program rttov_mie_params_aer

! Description:
!
! Generate a coefficient file containing aerosol optical parameters
!
! Usage: rttov_mie_params_aer.exe --rtcoef_file ... --opac_dir ...
!                                 [--nthreads n] [--temp_dir ...]
!
!        where rtcoef_file is a standard RTTOV op dep coefficient file
!              opac_dir is the path to the OPAC optdat/ directory
!              nthreads specifies the number of threads to use (if compiled
!                 with OpenMP), optional, defaults to 1
!              temp_dir is a directory to which, if supplied, temporary
!                 files will be written containing data for each particle
!                 type. This allows the processing to be restarted if it is
!                 interrupted, optional
!
! The new volcanic ash and Asian dust types are not from the OPAC database.
! The files containing the refractive indices are data/vapo00 and data/asdu00.
! They should be placed within the OPAC directory before running this program.
!
! For example, if the op dep coef file and optdat/ OPAC dir are in your current directory:
! $ rttov_mie_params_aer.exe --rtcoef_file rtcoef_msg_2_seviri.dat --opac_dir optdat
!
! The routine writes out an aerosol coef file with parameters for all channels
! in the input rtcoef file.
!
! The input file may be ASCII, binary or HDF5. The output aerosol coef file is the same
! format as the input file.
!
!
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 25 November 1998, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, KNMI and MeteoFrance.
!
!    Copyright 2012, EUMETSAT, All Rights Reserved.
!
! Method: Uses Mie theory to calculate aerosol optical parameters from
!         refractive indices taken from the OPAC database.
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.0    29/05/2012  Original code by Marco Matricardi, modified by JAH
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!

#include "throw.h"

use rttov_types, only : rttov_coefs

use rttov_const, only : &
  errorstatus_fatal,    &
  deg2rad,              &
  pi,                   &
  nphangle,             &
  phase_angles => phangle

use parkind1, only : jprb, jpim, jplm

use rttov_coef_io_mod, only : getlun, closelun

use rttov_getoptions, only : getoption

use rttov_unix_env, only: rttov_iargc, rttov_exit

use rttov_scattering_mod, only : &
  lognorm,    &
  gammadist,  &
  mie_sphere, &
  integrate,  &
  inter

#ifdef _RTTOV_HDF
use hdf5
use rttov_hdf_mod
#endif

implicit none

#include "rttov_read_ascii_coef.interface"
#include "rttov_read_binary_coef.interface"
#include "rttov_init_coef.interface"
#include "rttov_nullify_coef_scatt_ir.interface"
#include "rttov_write_ascii_scaercoef.interface"
#include "rttov_write_binary_scaercoef.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_errorreport.interface"
#include "rttov_cmpuc.interface"
#include "rttov_bpr_init.interface"
#include "rttov_bpr_calc.interface"
#include "rttov_bpr_dealloc.interface"

#ifdef _RTTOV_HDF
#include "rttov_hdf_load.interface"
#include "rttov_hdf_save.interface"
#endif

save
integer(kind=jpim), parameter :: iofileid  = 60    ! Logical unit used to read OPAC data
integer(kind=jpim), parameter :: nopacwavl = 61    ! Number of wavelengths in OPAC database
integer(kind=jpim), parameter :: nvapowavl = 52    ! Number of wavelengths in the new vol ash type file
integer(kind=jpim), parameter :: nasduwavl = 155   ! Number of wavelengths in the Asian dust type file
integer(kind=jpim), parameter :: maxnwavl  = max(nopacwavl, nvapowavl, nasduwavl)

integer(kind=jpim), parameter :: naercomp = 13     ! Number of aerosol types (components)
! integer(kind=jpim), parameter :: nphangle = 208    ! Number of phase angles
! integer(kind=jpim), parameter :: phangind90 = 118  ! Index of phase angle 90 degrees
integer(kind=jpim), parameter :: nrhum = 8         ! Number of relative humidities
integer(kind=jpim), parameter :: maxnradii = 10000 ! Maximum number of radius values
real(kind=jprb),    parameter :: rfac = 0.005      ! Multiplier for radius grid

character(len=4)   :: aernam(naercomp)      ! Aerosol component names in OPAC database
real(kind=jprb)    :: phangle(nphangle)     ! Phase angles
real(kind=jprb)    :: rhumval(nrhum)        ! Relative humidity values for water-affected types
character(len=2)   :: rhumstr(nrhum)        ! Strings of rel hum values
integer(kind=jpim) :: nrhumcomp(naercomp)   ! Number of relative humidities per component
integer(kind=jpim) :: compswitch(naercomp)  ! Switch to indicate which components to calculate

character(len=32)  :: format_in
character(len=256) :: rtcoef_file, f1
character(len=256) :: opac_dir
character(len=256) :: temp_dir
integer(kind=jpim) :: nthreads
character(len=256) :: fname
character(len=256) :: fnametemp
logical(kind=jplm) :: exists
integer(kind=jpim) :: err
integer(kind=jpim) :: file_id
character(len=32)  :: file_format
type(rttov_coefs)  :: coefs

type coef_phase_type
  real(kind=jprb), pointer :: phase(:,:,:)
endtype
type(coef_phase_type), allocatable :: coef_phase(:)

integer(kind=jpim) :: i, jr, ji, iaer, irhum, iwav, iang
integer(kind=jpim) :: nproc, iproc
integer(kind=jpim), allocatable :: iaerproc(:), irhumproc(:)
integer(kind=jpim), allocatable :: solar_chanlist(:) ! Temp array for identifying solar channels

real(kind=jprb)    :: sqrt2pi
integer(kind=jpim) :: nw                   ! Number of wavelengths in call to INTER
integer(kind=jpim) :: ntot                 ! Size of actual radius grid
real(kind=jprb)    :: rarr(maxnradii)      ! Radius grid for size dist.
real(kind=jprb)    :: n(maxnradii)         ! Size distribution
real(kind=jprb)    :: armin, armax         ! Size distribution parameters
real(kind=jprb)    :: rmin(naercomp,nrhum) ! Size distribution param.
real(kind=jprb)    :: rmax(naercomp,nrhum) ! Size distribution param.
real(kind=jprb)    :: rmod(naercomp,nrhum) ! Size distribution param.
real(kind=jprb)    :: sigma(naercomp)      ! Log-normal distribution param.
real(kind=jprb)    :: acoef(1)             ! Gamma distribution param.
real(kind=jprb)    :: alpha(1)             ! Gamma distribution param.
real(kind=jprb)    :: bcoef(1)             ! Gamma distribution param.
real(kind=jprb)    :: gamma(1)             ! Gamma distribution param.

integer(kind=jpim) :: iasdu                ! Variables for Asian dust type size dist. calculation:
real(kind=jprb)    :: nasdu(maxnradii)     ! Temporary number density array
integer(kind=jpim) :: asdu_opaci(3)        ! Indices of mineral components used
real(kind=jprb)    :: asdu_nd(3)           ! Number density scale factor for each mineral comp.

real(kind=jprb)    :: opacwavl(maxnwavl)  ! OPAC wavelengths
real(kind=jprb)    :: extc(maxnwavl)      ! OPAC parameters...
real(kind=jprb)    :: scac(maxnwavl)      ! These don't really need to be arrays
real(kind=jprb)    :: absc(maxnwavl)
real(kind=jprb)    :: sisca(maxnwavl)
real(kind=jprb)    :: asym(maxnwavl)
real(kind=jprb)    :: extnor(maxnwavl)
complex(kind=jprb) :: mopac(maxnwavl)     ! OPAC refractive index

real(kind=jprb)    :: refvapo(10)         ! Data read from new vol ash ref index file
integer(kind=jpim) :: nwr, nwi            ! Number of wavelengths for real/imag new vol ash ref indices
real(kind=jprb)    :: wavl                ! Wavelength read from new vol ash file
real(kind=jprb)    :: wavlr(nvapowavl)    ! Wavelengths for ref index real parts
real(kind=jprb)    :: wavli(nvapowavl)    ! Wavelengths for ref index imag parts
real(kind=jprb)    :: mopacr(nvapowavl)   ! Ref index real part for new vol ash
real(kind=jprb)    :: mopaci(nvapowavl)   ! Ref index imag part for new vol ash

real(kind=jprb)    :: refasdu(4)          ! Data read from Asian dust ref index file

real(kind=jprb)    :: minwavl             ! Smallest wavelength of ref index data for particle

real(kind=jprb)    :: mchanreal, mchanimag, hr  ! Arguments for INTER subroutine

real(kind=jprb), allocatable    :: chnwavl(:)   ! Instrument channel wavelengths
complex(kind=jprb), allocatable :: mchan(:)     ! Refrac. indices at chan wavelengths

real(kind=jprb) :: r, x, kl, r_angle, f, g      ! Used in calls to mie_sphere
real(kind=jprb) :: q_sct, q_ext, q_bsct, phasef
real(kind=jprb) :: xarr(maxnradii)
real(kind=jprb) :: yarr(maxnradii)
real(kind=jprb) :: zarr(maxnradii)
real(kind=jprb) :: parr(maxnradii)

integer(kind=jpim) :: ifail  ! Used in calls to NAG integration subroutine
real(kind=jprb)    :: error
real(kind=jprb)    :: ecoef
real(kind=jprb)    :: scoef
real(kind=jprb)    :: phfnc

! Variables used in bpr calculation
! integer(kind=jpim) :: j, k
! real(kind=jprb)    :: ang, ang1, arx(nphangle), ary(nphangle)
! real(kind=jprb)    :: xarr0(nphangle), xarr1(360), xarr2(phangind90), xarr3(nphangle-phangind90+1)
! real(kind=jprb)    :: yarr1(360), yarr2(phangind90), yarr3(nphangle-phangind90+1)
! real(kind=jprb)    :: intg(nphangle)
! real(kind=jprb)    :: cosT, phi, phasint, mu, mu1
! real(kind=jprb)    :: phas(nphangle,nphangle)

#ifndef RTTOV_NAG53
!$OMP THREADPRIVATE (mchan)
#endif


!-----Parameters used in the modified gamma distribution------------------------

data acoef /5461._jprb/
data alpha /1._jprb/
data bcoef /16._jprb/
data gamma /0.5_jprb/

!-------------------------------------------------------------------------------

!-----Parameters used in the log-normal distribution----------------------------

data rmin / 0.0050_jprb,  0.0050_jprb,  0.0050_jprb,  0.0050_jprb,  0.0050_jprb,  0.0050_jprb,  0.0050_jprb, &
            0.0050_jprb,  0.0200_jprb,  0.0050_jprb,  0.0050_jprb,  0.0050_jprb,  0.0100_jprb,  &
            0.0000_jprb,  0.0060_jprb,  0.0000_jprb,  0.0070_jprb,  0.0077_jprb,  0.0000_jprb,  0.0000_jprb, &
            0.0000_jprb,  0.0000_jprb,  0.0073_jprb,  0.0000_jprb,  0.0000_jprb,  0.0000_jprb,  &
            0.0000_jprb,  0.0064_jprb,  0.0000_jprb,  0.0085_jprb,  0.0085_jprb,  0.0000_jprb,  0.0000_jprb, &
            0.0000_jprb,  0.0000_jprb,  0.0079_jprb,  0.0000_jprb,  0.0000_jprb,  0.0000_jprb,  &
            0.0000_jprb,  0.0067_jprb,  0.0000_jprb,  0.0090_jprb,  0.0090_jprb,  0.0000_jprb,  0.0000_jprb, &
            0.0000_jprb,  0.0000_jprb,  0.0084_jprb,  0.0000_jprb,  0.0000_jprb,  0.0000_jprb,  &
            0.0000_jprb,  0.0071_jprb,  0.0000_jprb,  0.0099_jprb,  0.0099_jprb,  0.0000_jprb,  0.0000_jprb, &
            0.0000_jprb,  0.0000_jprb,  0.0090_jprb,  0.0000_jprb,  0.0000_jprb,  0.0000_jprb,  &
            0.0000_jprb,  0.0074_jprb,  0.0000_jprb,  0.0108_jprb,  0.0108_jprb,  0.0000_jprb,  0.0000_jprb, &
            0.0000_jprb,  0.0000_jprb,  0.0096_jprb,  0.0000_jprb,  0.0000_jprb,  0.0000_jprb,  &
            0.0000_jprb,  0.0078_jprb,  0.0000_jprb,  0.0115_jprb,  0.0115_jprb,  0.0000_jprb,  0.0000_jprb, &
            0.0000_jprb,  0.0000_jprb,  0.0101_jprb,  0.0000_jprb,  0.0000_jprb,  0.0000_jprb,  &
            0.0000_jprb,  0.0079_jprb,  0.0000_jprb,  0.0118_jprb,  0.0118_jprb,  0.0000_jprb,  0.0000_jprb, &
            0.0000_jprb,  0.0000_jprb,  0.0103_jprb,  0.0000_jprb,  0.0000_jprb,  0.0000_jprb   /

data rmax / 20.00_jprb,   20.00_jprb,   20.00_jprb,   20.00_jprb,   60.00_jprb,   20.00_jprb,   20.00_jprb, &
            60.00_jprb,    5.00_jprb,   20.00_jprb,   20.00_jprb,   20.00_jprb,   60.00_jprb,   &
             0.00_jprb,   25.00_jprb,    0.00_jprb,   32.20_jprb,   90.50_jprb,    0.00_jprb,    0.00_jprb, &
             0.00_jprb,    0.00_jprb,   30.20_jprb,    0.00_jprb,    0.00_jprb,    0.00_jprb,   &
             0.00_jprb,   27.40_jprb,    0.00_jprb,   36.30_jprb,  100.00_jprb,    0.00_jprb,    0.00_jprb, &
             0.00_jprb,    0.00_jprb,   33.60_jprb,    0.00_jprb,    0.00_jprb,    0.00_jprb,   &
             0.00_jprb,   29.80_jprb,    0.00_jprb,   39.90_jprb,  100.00_jprb,    0.00_jprb,    0.00_jprb, &
             0.00_jprb,    0.00_jprb,   36.40_jprb,    0.00_jprb,    0.00_jprb,    0.00_jprb,   &
             0.00_jprb,   35.00_jprb,    0.00_jprb,   47.80_jprb,  100.00_jprb,    0.00_jprb,    0.00_jprb, &
             0.00_jprb,    0.00_jprb,   42.30_jprb,    0.00_jprb,    0.00_jprb,    0.00_jprb,   &
             0.00_jprb,   42.80_jprb,    0.00_jprb,   58.50_jprb,  100.00_jprb,    0.00_jprb,    0.00_jprb, &
             0.00_jprb,    0.00_jprb,   50.20_jprb,    0.00_jprb,    0.00_jprb,    0.00_jprb,   &
             0.00_jprb,   58.50_jprb,    0.00_jprb,   78.30_jprb,  100.00_jprb,    0.00_jprb,    0.00_jprb, &
             0.00_jprb,    0.00_jprb,   65.50_jprb,    0.00_jprb,    0.00_jprb,    0.00_jprb,   &
             0.00_jprb,   74.90_jprb,    0.00_jprb,   98.50_jprb,  100.00_jprb,    0.00_jprb,    0.00_jprb, &
             0.00_jprb,    0.00_jprb,   81.80_jprb,    0.00_jprb,    0.00_jprb,    0.00_jprb    /

data rmod /0.4710_jprb,  0.0212_jprb,  0.0118_jprb,  0.2090_jprb,  1.7500_jprb,  0.0700_jprb,  0.3900_jprb, &
           1.9000_jprb,  0.5000_jprb,  0.0695_jprb,  0.0000_jprb,  0.610482_jprb,0.0000_jprb,  &
           0.0000_jprb,  0.0262_jprb,  0.0000_jprb,  0.3360_jprb,  2.8200_jprb,  0.0000_jprb,  0.0000_jprb, &
           0.0000_jprb,  0.0000_jprb,  0.0983_jprb,  0.0000_jprb,  0.0000_jprb,  0.0000_jprb,  &
           0.0000_jprb,  0.0285_jprb,  0.0000_jprb,  0.3780_jprb,  3.1700_jprb,  0.0000_jprb,  0.0000_jprb, &
           0.0000_jprb,  0.0000_jprb,  0.1090_jprb,  0.0000_jprb,  0.0000_jprb,  0.0000_jprb,  &
           0.0000_jprb,  0.0306_jprb,  0.0000_jprb,  0.4160_jprb,  3.4900_jprb,  0.0000_jprb,  0.0000_jprb, &
           0.0000_jprb,  0.0000_jprb,  0.1180_jprb,  0.0000_jprb,  0.0000_jprb,  0.0000_jprb,  &
           0.0000_jprb,  0.0348_jprb,  0.0000_jprb,  0.4970_jprb,  4.1800_jprb,  0.0000_jprb,  0.0000_jprb, &
           0.0000_jprb,  0.0000_jprb,  0.1350_jprb,  0.0000_jprb,  0.0000_jprb,  0.0000_jprb,  &
           0.0000_jprb,  0.0399_jprb,  0.0000_jprb,  0.6050_jprb,  5.1100_jprb,  0.0000_jprb,  0.0000_jprb, &
           0.0000_jprb,  0.0000_jprb,  0.1580_jprb,  0.0000_jprb,  0.0000_jprb,  0.0000_jprb,  &
           0.0000_jprb,  0.0476_jprb,  0.0000_jprb,  0.8010_jprb,  6.8400_jprb,  0.0000_jprb,  0.0000_jprb, &
           0.0000_jprb,  0.0000_jprb,  0.1950_jprb,  0.0000_jprb,  0.0000_jprb,  0.0000_jprb,  &
           0.0000_jprb,  0.0534_jprb,  0.0000_jprb,  0.9950_jprb,  8.5900_jprb,  0.0000_jprb,  0.0000_jprb, &
           0.0000_jprb,  0.0000_jprb,  0.2310_jprb,  0.0000_jprb,  0.0000_jprb,  0.0000_jprb   /

data sigma/2.510_jprb, 2.240_jprb, 2.000_jprb, 2.030_jprb, 2.030_jprb, 1.950_jprb, &
           2.000_jprb, 2.150_jprb, 2.200_jprb, 2.030_jprb, 0.000_jprb, 1.850_jprb, &
           0.000_jprb/

!-------------------------------------------------------------------------------

!-----Asian dust type size distribution data------------------------------------

data asdu_opaci /    6_jpim,     7_jpim,        8_jpim/
data asdu_nd    /0.862_jprb, 0.136_jprb, 0.217E-2_jprb/

!-------------------------------------------------------------------------------

!-----The name of the files stored in the OPAC database-------------------------

data aernam /'inso','waso','soot','ssam','sscm','minm','miam','micm','mitr',   &
             'suso','vola','vapo','asdu'/
!-------------------------------------------------------------------------------

!-----Relative humidity values--------------------------------------------------

data rhumval  /0.00_jprb, 50.00_jprb, 70.00_jprb, 80.00_jprb, 90.00_jprb, 95.00_jprb, 98.00_jprb, 99.00_jprb/

!-------------------------------------------------------------------------------

!-----Number of relative humidities for each aerosol type-----------------------

data nrhumcomp  /1, 8, 1, 8, 8, 1, 1, 1, 1, 8, 1, 1, 1/

!-------------------------------------------------------------------------------

!-----Switches to indicate which components to calculate------------------------

! Set to zero to disable computations: useful for testing
data compswitch /1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1/

!-------------------------------------------------------------------------------


TRY

!-------------------------------------------------------------------------------
! General info
!-------------------------------------------------------------------------------

! Each particle type requires a file containing refractive indices over a range
! of wavelengths and parameters defining a particle size distribution.

! All size distribution data is contained in this file (see above). All
! refractive index files (OPAC plus the two additional particle types) should
! be placed in a single directory whose path is passed as an argument to the
! executable.

! For particle types 1-11, the refractive indices come from the OPAC data base
! files. Types 1-10 assume a log-normal size distribution while type 11 assumes
! a modified gamma distribution.

! The size distribution for each particle/rel. humidity is calculated over the
! ranges of radii in rmin and rmax above. The mean and stdev of the log-normal
! distribution are contained in the rmod and sigma arrays.
! The parameters for the modified gamma distribution are contained in acoef,
! bcoef, alpha and gamma above.

! The vapo particle type (12) is based on a log-normal size distribution in the
! same way as types 1-11 and the refractive indices are contained in the file
! data/vapo00. The real and imaginary parts of the ref. index are interpolated
! over independent sets of wavelengths for vapo because they can contain missing
! real/imag values independently. Also, no values for wavelengths below 3um are
! available.

! The asdu particle type (13) is based on a combination of log-normal size
! distributions for particle types 6-8 (as in asdu_opaci) i.e. mineral nuc.,
! acc. and coarse modes. Each size distribution is calculated over the range
! 0.01-60um using the rmod and sigma values for types 6-8 and they are weighted
! according to the values in asdu_nd and summed. The refractive indices are
! contained in the file data/asdu00.

!-------------------------------------------------------------------------------
! Some initial set up
!-------------------------------------------------------------------------------

  if (rttov_iargc() == 0) then
    print *, 'Usage: --rtcoef_file   input rtcoef file name'
    print *, '       --opac_dir      path to directory containing OPAC data'
    print *, '       --nthreads      integer, optional'
    print *, '       --temp_dir      directory for temporary files, optional'
    stop
  endif

  call getoption("--rtcoef_file", rtcoef_file, mnd=.true._jplm)
  inquire(file=rtcoef_file, exist=exists)
  if (.not. exists) then
    print *, 'Cannot find rtcoef file: '//trim(rtcoef_file)
    stop
  endif

  call getoption("--opac_dir", opac_dir, mnd=.true._jplm)

  nthreads = -1
  call getoption("--nthreads", nthreads)
  if (nthreads <= 0) then
    nthreads = 1
  endif

  temp_dir = ''
  call getoption("--temp_dir", temp_dir)

  do i = 1, nrhum
    write(rhumstr(i), '(i2.2)') int(rhumval(i))
  enddo

  sqrt2pi = sqrt(2._jprb * pi)

  phangle(:) = phase_angles(:)
  call rttov_bpr_init(phangle, err)
  THROWM(err.ne.0, 'Error initialising bpr calculation tables')

!-------------------------------------------------------------------------------
! Read optical depth coefficient file
!-------------------------------------------------------------------------------

  call getlun(err, file_id, f1, file_format, .false._jplm, 'rtcoef', f=trim(rtcoef_file))
  THROW(err.ne.0)

  if (rttov_cmpuc(file_format, 'unformatted')) then
    call rttov_read_binary_coef(err, coefs%coef, file_id)
    THROWM(err.ne.0, 'Cannot open binary coefficient file '//trim(rtcoef_file))
    format_in = 'unformatted'
  else if (rttov_cmpuc(file_format, 'formatted')) then
    call rttov_read_ascii_coef(err, coefs%coef, file_id)
    THROWM(err.ne.0, 'Cannot open ASCII coefficient file '//trim(rtcoef_file))
    format_in = 'formatted'
  else if (rttov_cmpuc(file_format, 'hdf5')) then
#ifndef _RTTOV_HDF
    err =errorstatus_fatal
    THROWM(err.ne.0, 'This program is not compiled with HDF5 capability; use RTTOV_HDF=1 with Makefile.PL')
#else
    call open_hdf(.true., err)
    THROWM(err.ne.0, 'Error opening HDF5 interface')

    call rttov_hdf_load(err, f1, "/COEF", coef=coefs%coef)
    THROWM(err.ne.0, 'Cannot open HDF5 coefficient file '//trim(rtcoef_file))
    format_in = 'hdf5'

    call close_hdf(err)
    THROWM(err.ne.0, 'Error closing HDF5 interface')
#endif
  else
    err = errorstatus_fatal
    THROWM(err.ne.0, 'Unknown coefficient file format '//trim(file_format))
  endif

  call closelun(err, file_id)
  THROW(err.ne.0)

  call rttov_init_coef(err, coefs%coef)

!-------------------------------------------------------------------------------
! Populate the meta-data in the aerosol params structure
!-------------------------------------------------------------------------------

  call rttov_nullify_coef_scatt_ir(coefs%coef_scatt_ir)

  ! Number of channels
  coefs%coef_scatt_ir%fmv_aer_chn = coefs%coef%fmv_chn

  ! Create the index lookup for solar channels
  allocate(coefs%coef_scatt_ir%aer_pha_index(coefs%coef%fmv_chn))
  allocate(solar_chanlist(coefs%coef%fmv_chn))
  coefs%coef_scatt_ir%aer_pha_index(:) = 0
  coefs%coef_scatt_ir%fmv_aer_pha_chn = 0
  solar_chanlist(:) = 0

  do i = 1, coefs%coef%fmv_chn
    if (coefs%coef%ss_val_chn(i) > 0) then
      coefs%coef_scatt_ir%fmv_aer_pha_chn = coefs%coef_scatt_ir%fmv_aer_pha_chn + 1
      solar_chanlist(coefs%coef_scatt_ir%fmv_aer_pha_chn) = i
      coefs%coef_scatt_ir%aer_pha_index(i) = coefs%coef_scatt_ir%fmv_aer_pha_chn
    endif
  enddo

  coefs%coef_scatt_ir%fmv_aer_pha_ioff = 0
  if (coefs%coef_scatt_ir%fmv_aer_pha_chn > 0) then
    allocate(coefs%coef_scatt_ir%aer_pha_chanlist(coefs%coef_scatt_ir%fmv_aer_pha_chn))
    coefs%coef_scatt_ir%aer_pha_chanlist(:) = solar_chanlist(1:coefs%coef_scatt_ir%fmv_aer_pha_chn)
  endif
  deallocate(solar_chanlist)

  ! We now have:
  ! coefs%coef_scatt_ir%fmv_aer_pha_chn = number of solar channels (i.e. with phase functions)
  ! coefs%coef_scatt_ir%aer_pha_index(:) = index into phase array for each solar channel
  ! coefs%coef_scatt_ir%aer_pha_chanlist(:) = list of solar channel indexes

  ! Number of aerosol types
  coefs%coef_scatt_ir%fmv_aer_comp = naercomp
  allocate(coefs%coef_scatt_ir%fmv_aer_comp_name(naercomp))
  coefs%coef_scatt_ir%fmv_aer_comp_name(:) = aernam(1:naercomp)

  ! Allocate the optical params structure
  allocate(coefs%optp%optpaer(naercomp))

  do i = 1, naercomp
    call rttov_nullify_coef_scatt_ir(coefs%optp%optpaer(i))
  enddo

  ! Phase angles
  coefs%coef_scatt_ir%fmv_aer_ph = nphangle

  allocate(coefs%coef_scatt_ir%fmv_aer_ph_val(nphangle))
  coefs%coef_scatt_ir%fmv_aer_ph_val(:) = phangle(:)

  ! Relative humidities
  allocate(coefs%coef_scatt_ir%fmv_aer_rh(naercomp))
  coefs%coef_scatt_ir%fmv_aer_rh(:) = nrhumcomp(1:naercomp)

  do i = 1, naercomp
    allocate(coefs%optp%optpaer(i)%fmv_aer_rh_val(nrhumcomp(i)))
    coefs%optp%optpaer(i)%fmv_aer_rh_val(:) = rhumval(1:nrhumcomp(i))
  enddo

  ! Use a temporary array/struct for the phase functions as we need to calculate them
  ! for every channel, but we only keep them for the solar-affected channels.
  allocate(coef_phase(naercomp))

  ! Allocate arrays for optical parameters
  do i = 1, naercomp
    allocate(coefs%optp%optpaer(i)%abs(coefs%coef%fmv_chn,nrhumcomp(i)))
    allocate(coefs%optp%optpaer(i)%sca(coefs%coef%fmv_chn,nrhumcomp(i)))
    allocate(coefs%optp%optpaer(i)%bpr(coefs%coef%fmv_chn,nrhumcomp(i)))
    coefs%optp%optpaer(i)%abs(:,:)   = 0._jprb
    coefs%optp%optpaer(i)%sca(:,:)   = 0._jprb
    coefs%optp%optpaer(i)%bpr(:,:)   = 0._jprb
    if (coefs%coef_scatt_ir%fmv_aer_pha_chn > 0) then
      allocate(coefs%optp%optpaer(i)%pha(coefs%coef_scatt_ir%fmv_aer_pha_chn,nrhumcomp(i),nphangle))
      coefs%optp%optpaer(i)%pha(:,:,:) = 0._jprb
    endif
    allocate(coef_phase(i)%phase(coefs%coef%fmv_chn,nrhumcomp(i),nphangle))
    coef_phase(i)%phase(:,:,:) = 0._jprb
  enddo


  allocate(chnwavl(coefs%coef%fmv_chn))
  chnwavl(:) = 10000._jprb / coefs%coef%ff_cwn(:)

#ifndef RTTOV_NAG53
!$OMP PARALLEL NUM_THREADS(nthreads) DEFAULT(PRIVATE) SHARED(coefs)
#endif
  allocate(mchan(coefs%coef%fmv_chn))
#ifndef RTTOV_NAG53
!$OMP END PARALLEL
#endif

  ! Calculate number of particle type/rel hum combinations to process
  ! and set up arrays of indices for easy parallelisation.
  nproc = sum(nrhumcomp(1:naercomp) * compswitch(1:naercomp))
  allocate(iaerproc(nproc), irhumproc(nproc))
  i = 1
  do iaer = 1, naercomp
    if (compswitch(iaer) == 0) cycle
    do irhum = 1, nrhumcomp(iaer)
      iaerproc(i) = iaer
      irhumproc(i) = irhum
      i = i + 1
    enddo
  enddo

!-------------------------------------------------------------------------------
! For each aerosol type and rel hum, read and process the OPAC data
!-------------------------------------------------------------------------------
#ifndef RTTOV_NAG53
!$OMP PARALLEL DO NUM_THREADS(nthreads) DEFAULT(PRIVATE) SCHEDULE(DYNAMIC)                         &
!$OMP             SHARED(aernam, phangle, rhumval, rhumstr, nrhumcomp, coefs, opac_dir, chnwavl,   &
!$OMP                    rmod, rmin, rmax, sigma, acoef, bcoef, alpha, gamma, sqrt2pi, nproc,      &
!$OMP                    iaerproc, irhumproc, asdu_opaci, asdu_nd, coef_phase, temp_dir)
#endif
  do iproc = 1, nproc

    iaer = iaerproc(iproc)
    irhum = irhumproc(iproc)

    coefs%coef_scatt_ir%fmv_aer_comp_name(iaer) = aernam(iaer)

    ! Generate a unique file_id for this aerosol component and rel. hum.
    file_id = iofileid + iproc

    if (trim(temp_dir) /= '') then
      fnametemp = trim(temp_dir)//'/'//aernam(iaer)//rhumstr(irhum)//'.tmp'
      inquire(file=trim(fnametemp), exist=exists)
      if (exists) then
        open(file_id, file=trim(fnametemp), form='unformatted', status='old', action='read')
        read(file_id) coefs%optp%optpaer(iaer)%abs(:,irhum)
        read(file_id) coefs%optp%optpaer(iaer)%sca(:,irhum)
        read(file_id) coefs%optp%optpaer(iaer)%bpr(:,irhum)
        if (coefs%coef_scatt_ir%fmv_aer_pha_chn > 0) then
          read(file_id) coefs%optp%optpaer(iaer)%pha(:,irhum,:)
        endif
        close(file_id)
#ifndef RTTOV_NAG53
!$OMP CRITICAL
#endif
        print *,'Optical parameter data read from: '//trim(fnametemp)
#ifndef RTTOV_NAG53
!$OMP END CRITICAL
#endif
        cycle
      endif
    endif

    !------------------------------------------
    ! Open the OPAC file from current directory
    !------------------------------------------

    fname = trim(opac_dir)//'/'//aernam(iaer)//rhumstr(irhum)
    open(file_id, file=trim(fname), status='old', action='read', iostat=err)
#ifndef RTTOV_NAG53
!$OMP CRITICAL
#endif
    if (err == 0) then
      print *,'Reading ref. index file: '//trim(fname)
    else
      print *,'Cannot read ref. index file: '//trim(fname)//', skipping...'
    endif
#ifndef RTTOV_NAG53
!$OMP END CRITICAL
#endif
    if (err /= 0) cycle

    !------------------------------------------
    ! Read optical parameters from OPAC database
    !------------------------------------------

    if (iaer == 13) then  ! Asian dust (asdu)
      do i = 1, 2
        read(file_id, *)
      enddo

      nw = nasduwavl
      do i = 1, nasduwavl
        read (file_id, *) opacwavl(i), refasdu(:)
        mopac(i) = cmplx(refasdu(2),refasdu(3),jprb)
      enddo
      minwavl = minval(opacwavl)

    else if (iaer == 12) then  ! New vol ash (vapo)
      do i = 1, 6
        read(file_id, *)
      enddo

      ! The real and imaginary parts can be on different wavelengths as there can be missing
      ! data for one and not the other (vapo only)

      nwr = nvapowavl
      nwi = nvapowavl
      jr = 1
      ji = 1
      do i = 1, nvapowavl
        read (file_id, *) wavl, refvapo(:)
        if (refvapo(9) > -99) then
          wavlr(jr) = wavl
          mopacr(jr) = refvapo(9)
          jr = jr + 1
        else
          nwr = nwr - 1
        endif
        if (refvapo(10) > -99) then
          wavli(ji) = wavl
          mopaci(ji) = refvapo(10)
          ji = ji + 1
        else
          nwi = nwi - 1
        endif
      enddo
      minwavl = max(minval(wavlr(:)), minval(wavli(:)))

    else
      if (iaer == 11) then  ! OPAC vol ash (vola)
        do i = 1, 18
          read(file_id, *)
        enddo
      else                  ! All other types
        do i = 1, 17
          read(file_id, *)
        enddo
      endif

      nw = nopacwavl
      do i = 1, nopacwavl
        read (file_id, '(2x,7e10.3,2e11.3)') opacwavl(i), extc(i), scac(i), absc(i), &
                                                sisca(i), asym(i), extnor(i), mopac(i)
        mopac(i) = conjg(mopac(i))
      enddo
      minwavl = minval(opacwavl)
    endif

    close(file_id)

    !------------------------------------------
    ! Interpolate OPAC refractive indices to instrument wavelengths
    !------------------------------------------

    if (iaer == 12) then
      ! The real and imaginary parts can be on different wavelengths as there can be missing
      ! data for one and not the other (vapo only)

      do iwav = 1, coefs%coef%fmv_chn
        call INTER(nwr, nwr, 2, chnwavl(iwav), wavlr(1:nwr), mopacr(1:nwr), mchanreal, hr)
        call INTER(nwi, nwi, 2, chnwavl(iwav), wavli(1:nwi), mopaci(1:nwi), mchanimag, hr)
        mchan(iwav) = cmplx(mchanreal, mchanimag, kind=jprb)
      enddo

    else
      do iwav = 1, coefs%coef%fmv_chn
        call INTER(nw, nw, 2, chnwavl(iwav), opacwavl(1:nw), real(mopac(1:nw)),  mchanreal, hr)
        call INTER(nw, nw, 2, chnwavl(iwav), opacwavl(1:nw), aimag(mopac(1:nw)), mchanimag, hr)
        mchan(iwav) = cmplx(mchanreal, mchanimag, kind=jprb)
      enddo
    endif

    !------------------------------------------
    ! Generate grid of radius values over which to calculate the size dist.
    ! ntot is the number of values n the grid, rarr(:) contains the values
    !------------------------------------------

    armin  = rmin(iaer,irhum)
    armax  = rmax(iaer,irhum)

    rarr(1) = armin
    do i = 2, maxnradii
      rarr(i) = rarr(i-1) * 10**rfac
      if (rarr(i) > armax) then
        rarr(i) = armax
        exit
      endif
    enddo
    if (rarr(i) < armax) then
      print *, 'Increase maxnradii'
    endif
    ntot = i

    !------------------------------------------
    ! Compute size distribution
    !------------------------------------------

    if (iaer == 13) then
      ! The Asian dust type size distribution is a linear combination of the
      ! nuc., acc. and coarse mineral components. Calculate each size dist.
      ! and accumulate scaled values.

      n(1:ntot) = 0._jprb
      do iasdu = 1, size(asdu_nd)

        nasdu(:) = 0._jprb
        call lognorm(ntot, rarr(1:ntot), sigma(asdu_opaci(iasdu)), &
                     rmod(asdu_opaci(iasdu),1), sqrt2pi, nasdu(1:ntot))
        n(1:ntot) = n(1:ntot) + nasdu(1:ntot) * asdu_nd(iasdu)

      enddo

    else if (iaer <= 10 .or. iaer == 12) then
      ! For most aerosols the size distribution function is a log-normal distribution
      ! The size distribution function is computed for ntot points from rmin to rmax
      call lognorm(ntot, rarr(1:ntot), sigma(iaer), rmod(iaer,irhum), sqrt2pi, n(1:ntot))
    else
      ! For vola the size distribution function is a gamma distribution
      ! The size distribution function is computed for ntot points from rmin to rmax
      call gammadist(ntot, rarr(1:ntot), acoef(iaer-10), alpha(iaer-10), &
                     bcoef(iaer-10), gamma(iaer-10), n(1:ntot))
    endif

    !------------------------------------------
    ! Use the OPAC ref indices and the mie_sphere subroutine to calculate the optical params
    !------------------------------------------

    do iwav = 1, coefs%coef%fmv_chn

      if (chnwavl(iwav) < minwavl) then
#ifndef RTTOV_NAG53
!$OMP CRITICAL
#endif
        print *, 'Warning: wavelength too small for aer type '//aernam(iaer)
        print *, '         No params for channel ',iwav,' for this particle type'
#ifndef RTTOV_NAG53
!$OMP END CRITICAL
#endif
        coefs%optp%optpaer(iaer)%abs(iwav,irhum) = 0._jprb
        coefs%optp%optpaer(iaer)%sca(iwav,irhum) = 0._jprb
        coefs%optp%optpaer(iaer)%bpr(iwav,irhum) = 0._jprb
        coef_phase(iaer)%phase(iwav,irhum,:) = 0._jprb

      else

        do i = 1, ntot
          r  = rarr(i)
          x  = r * 2._jprb * pi / chnwavl(iwav)
          kl = 2._jprb * pi / chnwavl(iwav)
          r_angle = 0._jprb

          !---The optical parameters are computed using the Mie theory-------------------
          call mie_sphere(x, kl, r_angle, mchan(iwav), q_sct, q_ext, f, g, q_bsct, phasef)!, 0)
          !------------------------------------------------------------------------------

          ! Set up arrays to computed quantities for a sample of aerosols/cloud particles
          xarr(i) = r
          yarr(i) = q_ext          * n(i) * r**2 * pi
          zarr(i) = q_sct          * n(i) * r**2 * pi
        enddo

        !-Perform integration-----------------------------------------------------------
        call integrate(ntot, xarr(1:ntot), yarr(1:ntot), ecoef, error, ifail)
        call integrate(ntot, xarr(1:ntot), zarr(1:ntot), scoef, error, ifail)

        coefs%optp%optpaer(iaer)%abs(iwav,irhum) = (ecoef - scoef) * 0.001_jprb
        coefs%optp%optpaer(iaer)%sca(iwav,irhum) = scoef * 0.001_jprb
        !-------------------------------------------------------------------------------

        !------------------------------------------
        ! Calculate the phase function
        !------------------------------------------

        do iang = 1, nphangle

          r_angle = phangle(iang) * deg2rad

          !---The phase function is computed using the Mie theory ------------------------
          do i = 1, ntot
            r  = rarr(i)
            x  = r * 2._jprb * pi / chnwavl(iwav)
            kl = 2._jprb * pi / chnwavl(iwav)
            call mie_sphere(x, kl, r_angle, mchan(iwav), q_sct, q_ext, f, g, q_bsct, phasef)!, 1)
            parr(i) = phasef * n(i)
          enddo

          !-Perform integration-----------------------------------------------------------
          call integrate(ntot, xarr(1:ntot), parr(1:ntot), phfnc, error, ifail)

          coef_phase(iaer)%phase(iwav,irhum,iang) = 4._jprb * pi * phfnc / scoef
          !-------------------------------------------------------------------------------

        enddo

        !------------------------------------------
        ! Calculate the pseudo backscattering parameter bpr
        !------------------------------------------

        call rttov_bpr_calc(coef_phase(iaer)%phase(iwav,irhum,:), phangle, ecoef, err)
        coefs%optp%optpaer(iaer)%bpr(iwav,irhum) = ecoef

      endif !wavelength within range

    enddo !iwav

    ! coef_phase(:)%phase(:,:,:) contains phase functions for every channel. Copy just the
    ! phase functions for solar-affected channels into the coef pha(:,:,:) array
    if (coefs%coef_scatt_ir%fmv_aer_pha_chn > 0) then
      do i = 1, coefs%coef_scatt_ir%fmv_aer_pha_chn
        coefs%optp%optpaer(iaer)%pha(i,irhum,:) = coef_phase(iaer)%phase(coefs%coef_scatt_ir%aer_pha_chanlist(i),irhum,:)
      enddo
    endif

    if (trim(temp_dir) /= '') then
      open(file_id, file=trim(fnametemp), form='unformatted', status='replace', action='write', iostat=err)
      if (err == 0) then
        write(file_id) coefs%optp%optpaer(iaer)%abs(:,irhum)
        write(file_id) coefs%optp%optpaer(iaer)%sca(:,irhum)
        write(file_id) coefs%optp%optpaer(iaer)%bpr(:,irhum)
        if (coefs%coef_scatt_ir%fmv_aer_pha_chn > 0) then
          write(file_id) coefs%optp%optpaer(iaer)%pha(:,irhum,:)
        endif
        close(file_id)
      else
#ifndef RTTOV_NAG53
!$OMP CRITICAL
#endif
        print *, 'Failure to write temporary file: '//trim(fnametemp)
#ifndef RTTOV_NAG53
!$OMP END CRITICAL
#endif
      endif
    endif

  enddo !iproc
#ifndef RTTOV_NAG53
!$OMP END PARALLEL DO
#endif


!-------------------------------------------------------------------------------
! Write out aerosol coefficients
!-------------------------------------------------------------------------------

  call getlun(err, file_id, f1, file_format, .true._jplm, 'scaercoef', form=format_in, &
              instrument=(/coefs%coef%id_platform, coefs%coef%id_sat, coefs%coef%id_inst/))
  THROW(err.ne.0)

  if (trim(file_format) == 'formatted') then
    call rttov_write_ascii_scaercoef(err, coefs%coef, coefs%coef_scatt_ir, coefs%optp, &
                                     file_id, verbose=.true._jplm)
    THROWM(err.ne.0, 'Error writing ASCII aerosol coefficients')
  else if (trim(file_format) == 'unformatted') then
    call rttov_write_binary_scaercoef(err, coefs%coef, coefs%coef_scatt_ir, coefs%optp, &
                                      file_id, verbose=.true._jplm)
    THROWM(err.ne.0, 'Error writing binary aerosol coefficients')
  else if (trim(file_format) == 'hdf5') then
#ifndef _RTTOV_HDF
    err =errorstatus_fatal
    THROWM(err.ne.0, 'This program is not compiled with HDF5 capability; use RTTOV_HDF=1 with Makefile.PL')
#else
    call open_hdf(.true., err)
    THROWM(err.ne.0, 'Error opening HDF5 interface')

    call rttov_hdf_save(err, f1, "/SCAER", create=.true., &
                      & scaercoef=coefs%coef_scatt_ir, &
                      & optp=coefs%optp, compress=.true.)
    THROWM(err.ne.0, 'Error writing HDF5 aerosol coefficients')

    call close_hdf(err)
    THROWM(err.ne.0, 'Error closing HDF5 interface')
#endif
  else
    THROWM(err.ne.0, 'Unrecognised format for output')
  endif

  call closelun(err, file_id)
  THROW(err.ne.0)


!-------------------------------------------------------------------------------
! Clean up
!-------------------------------------------------------------------------------

  if (allocated(iaerproc))   deallocate(iaerproc)
  if (allocated(irhumproc))  deallocate(irhumproc)
#ifndef RTTOV_NAG53
!$OMP PARALLEL NUM_THREADS(nthreads) DEFAULT(PRIVATE)
#endif
  if (allocated(mchan))      deallocate(mchan)
#ifndef RTTOV_NAG53
!$OMP END PARALLEL
#endif
  if (allocated(chnwavl))    deallocate(chnwavl)

  if (allocated(coef_phase)) then
    do iaer = 1, naercomp
      if (associated(coef_phase(iaer)%phase)) deallocate(coef_phase(iaer)%phase)
    enddo
    deallocate(coef_phase)
  endif

  call rttov_bpr_dealloc(err)
  THROWM(err.ne.0, 'Error deallocating bpr calculation tables')

  call rttov_dealloc_coefs(err, coefs)
  THROWM(err.ne.0, 'Error deallocating coefficients')

PCATCH

end program rttov_mie_params_aer
