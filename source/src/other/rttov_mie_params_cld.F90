program rttov_mie_params_cld

! Description:
!
! Generate a coefficient file containing cloud optical parameters
!
! Usage: rttov_mie_params_cld.exe --rtcoef_file ... --opac_dir ...
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
! For example, if the op dep coef file and optdat/ OPAC dir are in your current directory:
! $ rttov_mie_params_cld.exe --rtcoef_file rtcoef_msg_2_seviri.dat --opac_dir optdat
!
! The routine writes out an cloud coef file with parameters for all channels
! in the input rtcoef file.
!
! The input file may be ASCII, binary or HDF5. The output cloud coef file is the same
! format as the input file.
!
! NB Currently does not calculate ice cloud parameters (water clouds only).
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
! Method: Uses Mie theory to calculate cloud optical parameters from
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
#include "rttov_write_ascii_sccldcoef.interface"
#include "rttov_write_binary_sccldcoef.interface"
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

integer(kind=jpim), parameter :: nopaccomp = 9     ! Number of cloud types in OPAC database
integer(kind=jpim), parameter :: nwclcomp = 5      ! Number of water cloud types
! integer(kind=jpim), parameter :: nphangle = 208    ! Number of phase angles
! integer(kind=jpim), parameter :: phangind90 = 118  ! Index of phase angle 90 degrees
integer(kind=jpim), parameter :: nrhum = 8         ! Number of relative humidities
integer(kind=jpim), parameter :: maxnradii = 10000 ! Maximum number of radius values
real(kind=jprb),    parameter :: rfac = 0.005      ! Multiplier for radius grid

character(len=4)   :: compnam(nopaccomp)    ! Cloud component names in OPAC database
real(kind=jprb)    :: phangle(nphangle)     ! Phase angles
real(kind=jprb)    :: rhumval(nrhum)        ! Relative humidity values for water-affected types
character(len=2)   :: rhumstr(nrhum)        ! Strings of rel hum values
integer(kind=jpim) :: nrhumcomp(nopaccomp)  ! Number of relative humidities per component
real(kind=jprb)    :: confac(nwclcomp)      ! LWC to number density conversion factors
integer(kind=jpim) :: compswitch(nwclcomp)  ! Switch to indicate which components to calculate

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
type(coef_phase_type), allocatable :: coef_phase_wcl(:)
! type(coef_phase_type), allocatable :: coef_phase_icl(:)

integer(kind=jpim) :: i, iwcl, irhum, iwav, iang
integer(kind=jpim) :: nproc, iproc
integer(kind=jpim), allocatable :: iwclproc(:), irhumproc(:)
integer(kind=jpim), allocatable :: solar_chanlist(:) ! Temp array for identifying solar channels

real(kind=jprb)    :: sqrt2pi
integer(kind=jpim) :: ntot                 ! Size of actual radius grid
real(kind=jprb)    :: rarr(maxnradii)      ! Radius grid for size dist.
real(kind=jprb)    :: n(maxnradii)         ! Size distribution
real(kind=jprb)    :: armin, armax         ! Size distribution parameters
real(kind=jprb)    :: rmin(nopaccomp)      ! Size distribution param.
real(kind=jprb)    :: rmax(nopaccomp)      ! Size distribution param.
real(kind=jprb)    :: acoef(5)             ! Gamma distribution param.
real(kind=jprb)    :: alpha(5)             ! Gamma distribution param.
real(kind=jprb)    :: bcoef(5)             ! Gamma distribution param.
real(kind=jprb)    :: gamma(5)             ! Gamma distribution param.

real(kind=jprb)    :: opacwavl(nopacwavl)  ! OPAC wavelengths
real(kind=jprb)    :: extc(nopacwavl)      ! OPAC parameters...
real(kind=jprb)    :: scac(nopacwavl)      ! These don't really need to be arrays
real(kind=jprb)    :: absc(nopacwavl)
real(kind=jprb)    :: sisca(nopacwavl)
real(kind=jprb)    :: asym(nopacwavl)
real(kind=jprb)    :: extnor(nopacwavl)
complex(kind=jprb) :: mopac(nopacwavl)     ! OPAC refractive index

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

data acoef/9.790E-03_jprb,3.820E-03_jprb,1.110E-03_jprb,8.120E-04_jprb,5.670E-05_jprb/
data alpha/5.000E+00_jprb,3.000E+00_jprb,5.000E+00_jprb,8.000E+00_jprb,4.000E+00_jprb/
data bcoef/9.380E-01_jprb,1.930E-01_jprb,7.820E-02_jprb,2.470E-01_jprb,7.130E-03_jprb/
data gamma/1.050E+00_jprb,1.300E+00_jprb,2.160E+00_jprb,2.150E+00_jprb,2.340E+00_jprb/

!-------------------------------------------------------------------------------

!-----Parameters used in the size distribution----------------------------------

data rmin /0.02_jprb, 0.02_jprb, 0.02_jprb, 0.02_jprb, 0.02_jprb, 0.02_jprb, 20.0_jprb, 20.0_jprb, 2.0_jprb/

data rmax /50.0_jprb, 50.0_jprb, 50.0_jprb, 50.0_jprb, 50.0_jprb, 50.0_jprb, 2000._jprb, 2000._jprb, 2000._jprb/

!-------------------------------------------------------------------------------

!-----The name of the files stored in the OPAC database-------------------------

data compnam /'stco','stma','cucc','cucp','cuma','fogr','cir1','cir2','cir3'/

!-------------------------------------------------------------------------------

!-----Relative humidity values--------------------------------------------------

data rhumval  /0.00_jprb, 50.00_jprb, 70.00_jprb, 80.00_jprb, 90.00_jprb, 95.00_jprb, 98.00_jprb, 99.00_jprb/

!-------------------------------------------------------------------------------

!-----Number of relative humidities for each cloud type-------------------------

data nrhumcomp  /1, 1, 1, 1, 1, 1, 1, 1, 1/

!-------------------------------------------------------------------------------

!-----LWC to number density conversion factor-----------------------------------

data confac /892.857_jprb, 266.667_jprb, 1538.461_jprb, 4347.826_jprb, 147.710_jprb/

!-------------------------------------------------------------------------------

!-----Switches to indicate which components to calculate------------------------

! Set to zero to disable computations: useful for testing
data compswitch /1, 1, 1, 1, 1/

!-------------------------------------------------------------------------------

TRY

!-------------------------------------------------------------------------------
! General info
!-------------------------------------------------------------------------------

! Each particle type requires a file containing refractive indices over a range
! of wavelengths and parameters defining a particle size distribution.

! All size distribution data is contained in this file (see above). All
! refractive index files (from the OPAC database) should be placed in a single
! directory whose path is passed as an argument to the executable.

! All water cloud types assume a modified gamma distribution.

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
! Populate the meta-data in the cloud params structure
!-------------------------------------------------------------------------------

  call rttov_nullify_coef_scatt_ir(coefs%coef_scatt_ir)

  ! Number of channels
  coefs%coef_scatt_ir%fmv_wcl_chn = coefs%coef%fmv_chn
  coefs%coef_scatt_ir%fmv_icl_chn = coefs%coef%fmv_chn

  ! Create the index lookup for solar channels
  allocate(coefs%coef_scatt_ir%wcl_pha_index(coefs%coef%fmv_chn))
  allocate(coefs%coef_scatt_ir%icl_pha_index(coefs%coef%fmv_chn))
  allocate(solar_chanlist(coefs%coef%fmv_chn))
  coefs%coef_scatt_ir%wcl_pha_index(:) = 0
  coefs%coef_scatt_ir%icl_pha_index(:) = 0
  coefs%coef_scatt_ir%fmv_wcl_pha_chn = 0
  coefs%coef_scatt_ir%fmv_icl_pha_chn = 0
  solar_chanlist(:) = 0

  do i = 1, coefs%coef%fmv_chn
    if (coefs%coef%ss_val_chn(i) > 0) then
      coefs%coef_scatt_ir%fmv_wcl_pha_chn = coefs%coef_scatt_ir%fmv_wcl_pha_chn + 1
      coefs%coef_scatt_ir%fmv_icl_pha_chn = coefs%coef_scatt_ir%fmv_icl_pha_chn + 1
      solar_chanlist(coefs%coef_scatt_ir%fmv_wcl_pha_chn) = i
      coefs%coef_scatt_ir%wcl_pha_index(i) = coefs%coef_scatt_ir%fmv_wcl_pha_chn
      coefs%coef_scatt_ir%icl_pha_index(i) = coefs%coef_scatt_ir%fmv_icl_pha_chn
    endif
  enddo

  coefs%coef_scatt_ir%fmv_wcl_pha_ioff = 0
  coefs%coef_scatt_ir%fmv_icl_pha_ioff = 0
  if (coefs%coef_scatt_ir%fmv_wcl_pha_chn > 0) then
    allocate(coefs%coef_scatt_ir%wcl_pha_chanlist(coefs%coef_scatt_ir%fmv_wcl_pha_chn))
    coefs%coef_scatt_ir%wcl_pha_chanlist(:) = solar_chanlist(1:coefs%coef_scatt_ir%fmv_wcl_pha_chn)
  endif
  if (coefs%coef_scatt_ir%fmv_icl_pha_chn > 0) then
    allocate(coefs%coef_scatt_ir%icl_pha_chanlist(coefs%coef_scatt_ir%fmv_icl_pha_chn))
    coefs%coef_scatt_ir%icl_pha_chanlist(:) = solar_chanlist(1:coefs%coef_scatt_ir%fmv_icl_pha_chn)
  endif
  deallocate(solar_chanlist)

  ! We now have:
  ! coefs%coef_scatt_ir%fmv_?cl_pha_chn = number of solar channels (i.e. with phase functions)
  ! coefs%coef_scatt_ir%?cl_pha_index(:) = index into phase array for each solar channel
  ! coefs%coef_scatt_ir%?cl_pha_chanlist(:) = list of solar channel indexes

  ! Number of cloud types
  coefs%coef_scatt_ir%fmv_wcl_comp = nwclcomp
  allocate(coefs%coef_scatt_ir%fmv_wcl_comp_name(nwclcomp))
  coefs%coef_scatt_ir%fmv_wcl_comp_name(:) = compnam(1:nwclcomp)

  ! Water cloud LWC-number density conversion factors
  allocate(coefs%coef_scatt_ir%confac(coefs%coef_scatt_ir%fmv_wcl_comp))
  coefs%coef_scatt_ir%confac(:) = 0._jprb

  ! Ice cloud data/arrays
  coefs%coef_scatt_ir%icl_nabs = 4
  coefs%coef_scatt_ir%icl_nsca = 4
  coefs%coef_scatt_ir%icl_nbpr = 4
  coefs%coef_scatt_ir%fmv_icl_comp = 30
  coefs%coef_scatt_ir%fmv_icl_ishp = 2
  allocate(coefs%coef_scatt_ir%fmv_icl_dg(coefs%coef_scatt_ir%fmv_icl_comp, coefs%coef_scatt_ir%fmv_icl_ishp))
  coefs%coef_scatt_ir%fmv_icl_dg(:,:) = 0._jprb
  allocate(coefs%coef_scatt_ir%fmv_icl_comp_name(coefs%coef_scatt_ir%fmv_icl_comp, coefs%coef_scatt_ir%fmv_icl_ishp))
  coefs%coef_scatt_ir%fmv_icl_comp_name(:,1) = 'Hexagonal'
  coefs%coef_scatt_ir%fmv_icl_comp_name(:,2) = 'Aggregate'

  ! Allocate the optical params structure
  allocate(coefs%optp%optpwcl(nwclcomp))
  allocate(coefs%optp%optpicl(coefs%coef_scatt_ir%fmv_icl_ishp))

  do i = 1, nwclcomp
    call rttov_nullify_coef_scatt_ir(coefs%optp%optpwcl(i))
  enddo
  do i = 1, coefs%coef_scatt_ir%fmv_icl_ishp
    call rttov_nullify_coef_scatt_ir(coefs%optp%optpicl(i))
  enddo

  ! Phase angles
  coefs%coef_scatt_ir%fmv_wcl_ph = nphangle
  coefs%coef_scatt_ir%fmv_icl_ph = nphangle

  allocate(coefs%coef_scatt_ir%fmv_wcl_ph_val(nphangle))
  allocate(coefs%coef_scatt_ir%fmv_icl_ph_val(nphangle))
  coefs%coef_scatt_ir%fmv_wcl_ph_val(:) = phangle(:)
  coefs%coef_scatt_ir%fmv_icl_ph_val(:) = phangle(:)

  ! Relative humidities
  allocate(coefs%coef_scatt_ir%fmv_wcl_rh(nwclcomp))
  coefs%coef_scatt_ir%fmv_wcl_rh(:) = nrhumcomp(1:nwclcomp)

  do i = 1, nwclcomp
    allocate(coefs%optp%optpwcl(i)%fmv_wcl_rh_val(nrhumcomp(i)))
    coefs%optp%optpwcl(i)%fmv_wcl_rh_val(:) = rhumval(1:nrhumcomp(i))
  enddo

  ! Use temporary arrays/structs for the phase functions as we need to calculate them
  ! for every channel, but we only keep them for the solar-affected channels.
  allocate(coef_phase_wcl(nwclcomp))
!   allocate(coef_phase_icl(coefs%coef_scatt_ir%fmv_icl_ishp))

  ! Allocate arrays for optical parameters
  do i = 1, nwclcomp
    allocate(coefs%optp%optpwcl(i)%abs(coefs%coef%fmv_chn,nrhumcomp(i)))
    allocate(coefs%optp%optpwcl(i)%sca(coefs%coef%fmv_chn,nrhumcomp(i)))
    allocate(coefs%optp%optpwcl(i)%bpr(coefs%coef%fmv_chn,nrhumcomp(i)))
    coefs%optp%optpwcl(i)%abs(:,:)   = 0._jprb
    coefs%optp%optpwcl(i)%sca(:,:)   = 0._jprb
    coefs%optp%optpwcl(i)%bpr(:,:)   = 0._jprb
    if (coefs%coef_scatt_ir%fmv_wcl_pha_chn > 0) then
      allocate(coefs%optp%optpwcl(i)%pha(coefs%coef_scatt_ir%fmv_wcl_pha_chn,nrhumcomp(i),nphangle))
      coefs%optp%optpwcl(i)%pha(:,:,:) = 0._jprb
    endif
    allocate(coef_phase_wcl(i)%phase(coefs%coef%fmv_chn,nrhumcomp(i),nphangle))
    coef_phase_wcl(i)%phase(:,:,:) = 0._jprb
  enddo

  do i = 1, coefs%coef_scatt_ir%fmv_icl_ishp
    allocate(coefs%optp%optpicl(i)%abs(coefs%coef%fmv_chn,coefs%coef_scatt_ir%icl_nabs))
    allocate(coefs%optp%optpicl(i)%sca(coefs%coef%fmv_chn,coefs%coef_scatt_ir%icl_nsca))
    allocate(coefs%optp%optpicl(i)%bpr(coefs%coef%fmv_chn,coefs%coef_scatt_ir%icl_nbpr))
    coefs%optp%optpicl(i)%abs(:,:)   = 0._jprb
    coefs%optp%optpicl(i)%sca(:,:)   = 0._jprb
    coefs%optp%optpicl(i)%bpr(:,:)   = 0._jprb
    if (coefs%coef_scatt_ir%fmv_icl_pha_chn > 0) then
      allocate(coefs%optp%optpicl(i)%pha(coefs%coef_scatt_ir%fmv_icl_pha_chn,coefs%coef_scatt_ir%fmv_icl_comp,nphangle))
      coefs%optp%optpicl(i)%pha(:,:,:) = 0._jprb
    endif
!     allocate(coef_phase_icl(i)%phase(coefs%coef%fmv_chn,coefs%coef_scatt_ir%fmv_icl_comp,nphangle))
!     coef_phase_icl(i)%phase(:,:,:) = 0._jprb
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
  nproc = sum(nrhumcomp(1:nwclcomp) * compswitch(1:nwclcomp))
  allocate(iwclproc(nproc), irhumproc(nproc))
  i = 1
  do iwcl = 1, nwclcomp
    if (compswitch(iwcl) == 0) cycle
    do irhum = 1, nrhumcomp(iwcl)
      iwclproc(i) = iwcl
      irhumproc(i) = irhum
      i = i + 1
    enddo
  enddo

!-------------------------------------------------------------------------------
! For each cloud type and rel hum, read and process the OPAC data
!-------------------------------------------------------------------------------
#ifndef RTTOV_NAG53
!$OMP PARALLEL DO NUM_THREADS(nthreads) DEFAULT(PRIVATE) SCHEDULE(DYNAMIC)                                &
!$OMP             SHARED(compnam, phangle, rhumval, rhumstr, nrhumcomp, confac, coefs, opac_dir, chnwavl, &
!$OMP                    rmin, rmax, acoef, bcoef, alpha, gamma, sqrt2pi, nproc, iwclproc, irhumproc,     &
!$OMP                    coef_phase_wcl, temp_dir) !, coef_phase_icl)
#endif
  do iproc = 1, nproc

    iwcl = iwclproc(iproc)
    irhum = irhumproc(iproc)

    coefs%coef_scatt_ir%fmv_wcl_comp_name(iwcl) = compnam(iwcl)
    coefs%coef_scatt_ir%confac(iwcl) = confac(iwcl)

    ! Generate a unique file_id for this cloud component and rel. hum.
    file_id = iofileid + iproc

    if (trim(temp_dir) /= '') then
      fnametemp = trim(temp_dir)//'/'//compnam(iwcl)//rhumstr(irhum)//'.tmp'
      inquire(file=trim(fnametemp), exist=exists)
      if (exists) then
        open(file_id, file=trim(fnametemp), form='unformatted', status='old', action='read')
        read(file_id) coefs%optp%optpwcl(iwcl)%abs(:,irhum)
        read(file_id) coefs%optp%optpwcl(iwcl)%sca(:,irhum)
        read(file_id) coefs%optp%optpwcl(iwcl)%bpr(:,irhum)
        if (coefs%coef_scatt_ir%fmv_wcl_pha_chn > 0) then
          read(file_id) coefs%optp%optpwcl(iwcl)%pha(:,irhum,:)
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

    fname = trim(opac_dir)//'/'//compnam(iwcl)//rhumstr(irhum)
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

    do i = 1, 18
      read(file_id, *)
    enddo

    do i = 1, nopacwavl
      read (file_id, '(2x,7e10.3,2e11.3)') opacwavl(i), extc(i), scac(i), absc(i), &
                                              sisca(i), asym(i), extnor(i), mopac(i)
      mopac(i) = conjg(mopac(i))
    enddo

    close(file_id)

    !------------------------------------------
    ! Interpolate OPAC refractive indices to instrument wavelengths
    !------------------------------------------

    do iwav = 1, coefs%coef%fmv_chn
      call INTER(nopacwavl, nopacwavl, 2, chnwavl(iwav), opacwavl, real(mopac),  mchanreal, hr)
      call INTER(nopacwavl, nopacwavl, 2, chnwavl(iwav), opacwavl, aimag(mopac), mchanimag, hr)
      mchan(iwav) = cmplx(mchanreal, mchanimag, kind=jprb)
    enddo

    !------------------------------------------
    ! Generate grid of radius values over which to calculate the size dist.
    ! ntot is the number of values n the grid, rarr(:) contains the values
    !------------------------------------------

    armin  = rmin(iwcl)
    armax  = rmax(iwcl)

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

    call gammadist(ntot, rarr(1:ntot), acoef(iwcl), alpha(iwcl), &
                   bcoef(iwcl), gamma(iwcl), n(1:ntot))

    !------------------------------------------
    ! Use the OPAC ref indices and the mie_sphere subroutine to calculate the optical params
    !------------------------------------------

    do iwav = 1, coefs%coef%fmv_chn

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

      coefs%optp%optpwcl(iwcl)%abs(iwav,irhum) = (ecoef - scoef) * 0.001_jprb
      coefs%optp%optpwcl(iwcl)%sca(iwav,irhum) = scoef * 0.001_jprb
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

        coef_phase_wcl(iwcl)%phase(iwav,irhum,iang) = 4._jprb * pi * phfnc / scoef
        !-------------------------------------------------------------------------------

      enddo

      !------------------------------------------
      ! Calculate the pseudo backscattering parameter bpr
      !------------------------------------------

      call rttov_bpr_calc(coef_phase_wcl(iwcl)%phase(iwav,irhum,:), phangle, ecoef, err)
      coefs%optp%optpwcl(iwcl)%bpr(iwav,irhum) = ecoef

    enddo !iwav

    ! coef_phase_wcl(:)%phase(:,:,:) contains phase functions for every channel. Copy just the
    ! phase functions for solar-affected channels into the coef pha(:,:,:) array
    if (coefs%coef_scatt_ir%fmv_wcl_pha_chn > 0) then
      do i = 1, coefs%coef_scatt_ir%fmv_wcl_pha_chn
        coefs%optp%optpwcl(iwcl)%pha(i,irhum,:) = coef_phase_wcl(iwcl)%phase(coefs%coef_scatt_ir%wcl_pha_chanlist(i),irhum,:)
      enddo
    endif

    if (trim(temp_dir) /= '') then
      open(file_id, file=trim(fnametemp), form='unformatted', status='replace', action='write', iostat=err)
      if (err == 0) then
        write(file_id) coefs%optp%optpwcl(iwcl)%abs(:,irhum)
        write(file_id) coefs%optp%optpwcl(iwcl)%sca(:,irhum)
        write(file_id) coefs%optp%optpwcl(iwcl)%bpr(:,irhum)
        if (coefs%coef_scatt_ir%fmv_wcl_pha_chn > 0) then
          write(file_id) coefs%optp%optpwcl(iwcl)%pha(:,irhum,:)
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
! Write out cloud coefficients
!-------------------------------------------------------------------------------

  call getlun(err, file_id, f1, file_format, .true._jplm, 'sccldcoef', form=format_in, &
              instrument=(/coefs%coef%id_platform, coefs%coef%id_sat, coefs%coef%id_inst/))
  THROW(err.ne.0)

  if (trim(file_format) == 'formatted') then
    call rttov_write_ascii_sccldcoef(err, coefs%coef, coefs%coef_scatt_ir, coefs%optp, &
                                     file_id, verbose=.true._jplm)
    THROWM(err.ne.0, 'Error writing ASCII cloud coefficients')
  else if (trim(file_format) == 'unformatted') then
    call rttov_write_binary_sccldcoef(err, coefs%coef, coefs%coef_scatt_ir, coefs%optp, &
                                      file_id, verbose=.true._jplm)
    THROWM(err.ne.0, 'Error writing binary cloud coefficients')
  else if (trim(file_format) == 'hdf5') then
#ifndef _RTTOV_HDF
    err =errorstatus_fatal
    THROWM(err.ne.0, 'This program is not compiled with HDF5 capability; use RTTOV_HDF=1 with Makefile.PL')
#else
    call open_hdf(.true., err)
    THROWM(err.ne.0, 'Error opening HDF5 interface')

    call rttov_hdf_save(err, f1, "/SCCLD", create=.true., &
                      & sccldcoef=coefs%coef_scatt_ir, &
                      & optp=coefs%optp, compress=.true.)
    THROWM(err.ne.0, 'Error writing HDF5 cloud coefficients')

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

  if (allocated(iwclproc))   deallocate(iwclproc)
  if (allocated(irhumproc))  deallocate(irhumproc)
#ifndef RTTOV_NAG53
!$OMP PARALLEL NUM_THREADS(nthreads) DEFAULT(PRIVATE)
#endif
  if (allocated(mchan))      deallocate(mchan)
#ifndef RTTOV_NAG53
!$OMP END PARALLEL
#endif
  if (allocated(chnwavl))    deallocate(chnwavl)

  if (allocated(coef_phase_wcl)) then
    do iwcl = 1, nwclcomp
      if (associated(coef_phase_wcl(iwcl)%phase)) deallocate(coef_phase_wcl(iwcl)%phase)
    enddo
    deallocate(coef_phase_wcl)
  endif
!   if (allocated(coef_phase_icl)) then
!     do iicl = 1, coefs%coef_scatt_ir%fmv_icl_ishp
!       if (associated(coef_phase_icl(iicl)%phase)) deallocate(coef_phase_icl(iicl)%phase)
!     enddo
!     deallocate(coef_phase_icl)
!   endif

  call rttov_bpr_dealloc(err)
  THROWM(err.ne.0, 'Error deallocating bpr calculation tables')

  call rttov_dealloc_coefs(err, coefs)
  THROWM(err.ne.0, 'Error deallocating coefficients')

PCATCH

end program rttov_mie_params_cld
