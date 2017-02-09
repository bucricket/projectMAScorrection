PROGRAM example_aer_param_fwd
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
  !    Copyright 2013, EUMETSAT, All Rights Reserved.
  !
  !     *************************************************************
  !
  !     TEST PROGRAM FOR RTTOV FORWARD MODEL WITH AEROSOLS
  !          RTTOV VERSION 11
  !
  ! Demonstrating how to supply aerosol optical parameters
  ! explicitly to RTTOV.
  !
  ! To run this program you must have the following files
  ! either resident in the same directory or set up as a
  ! symbolic link:
  !   the file containing input profiles (e.g. prof.dat)
  !   the file containing aerosol optical parameter profiles for each channel
  !   the RTTOV coefficient file
  !
  ! The script run_example_aer_param_fwd.sh may be used to run this program.
  ! The output is generated in a file called example_aer_param_fwd_output.dat.
  !
  ! NB This program is very similar to example_fwd.F90 with the addition
  !    of the aerosol scattering.
  !
  ! To adapt the code to their own applications users should
  ! edit the code highlighted like this:
  !     !================================
  !     !======Read =====start===========
  !          code to be modified
  !     !======Read ===== end ===========
  !     !================================
  !
  ! Current Code Owner: SAF NWP
  !
  ! Code Description:
  !   Language:          Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  !

  ! rttov_const contains useful RTTOV constants
  USE rttov_const, ONLY :     &
       & errorstatus_success, &
       & errorstatus_fatal,   &
       & platform_name,       &
       & inst_name

  ! rttov_types contains definitions of all RTTOV data types
  USE rttov_types, ONLY :     &
       & rttov_options,       &
       & rttov_coefs,         &
       & profile_type,        &
       & transmission_type,   &
       & radiance_type,       &
       & rttov_chanprof,      &
       & rttov_emissivity,    &
       & rttov_reflectance,   &
       & rttov_opt_param

  ! jpim, jprb and jplm are the RTTOV integer, real and logical KINDs
  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_unix_env, ONLY : rttov_exit

  IMPLICIT NONE

#include "rttov_direct.interface"
#include "rttov_parallel_direct.interface"
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_init_opt_param.interface"
#include "rttov_bpr_init.interface"
#include "rttov_bpr_calc.interface"
#include "rttov_bpr_dealloc.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_skipcommentline.interface"

  !--------------------------
  !
  INTEGER(KIND=jpim), PARAMETER :: iup   = 20   ! unit for input profile file
  INTEGER(KIND=jpim), PARAMETER :: ioout = 21   ! unit for output

  ! RTTOV variables/structures
  !====================
  TYPE(rttov_options)              :: opts                     ! Options structure
  TYPE(rttov_coefs)                :: coefs                    ! Coefficients structure
  TYPE(rttov_chanprof),    POINTER :: chanprof(:)    => NULL() ! Input channel/profile list
  LOGICAL(KIND=jplm),      POINTER :: calcemis(:)    => NULL() ! Flag to indicate calculation of emissivity within RTTOV
  TYPE(rttov_emissivity),  POINTER :: emissivity(:)  => NULL() ! Input/output surface emissivity
  LOGICAL(KIND=jplm),      POINTER :: calcrefl(:)    => NULL() ! Flag to indicate calculation of BRDF within RTTOV
  TYPE(rttov_reflectance), POINTER :: reflectance(:) => NULL() ! Input/output surface BRDF
  TYPE(profile_type),      POINTER :: profiles(:)    => NULL() ! Input profiles
  TYPE(rttov_opt_param)            :: aer_opt_param            ! Input aerosol optical parameters
  TYPE(transmission_type)          :: transmission             ! Output transmittances
  TYPE(radiance_type)              :: radiance                 ! Output radiances

  INTEGER(KIND=jpim)               :: errorstatus              ! Return error status of RTTOV subroutine calls

  INTEGER(KIND=jpim) :: alloc_status
  CHARACTER(LEN=21)  :: NameOfRoutine = 'example_aer_param_fwd'

  ! variables for input
  !====================
  CHARACTER(LEN=256) :: coef_filename
  CHARACTER(LEN=256) :: prof_filename
  CHARACTER(LEN=256) :: file_aer_opt_param
  INTEGER(KIND=jpim) :: nthreads
  INTEGER(KIND=jpim) :: dosolar
  INTEGER(KIND=jpim) :: nlevels
  INTEGER(KIND=jpim) :: nprof
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: nchanprof
  INTEGER(KIND=jpim), ALLOCATABLE :: channel_list(:)
  INTEGER(KIND=jpim) :: nphangle
  ! loop variables
  INTEGER(KIND=jpim) :: i, j, jch
  INTEGER(KIND=jpim) :: nch
  INTEGER(KIND=jpim) :: iprof, joff
  INTEGER            :: ios

  !- End of header --------------------------------------------------------

  ! The usual steps to take when running RTTOV are as follows:
  !   1. Specify required RTTOV options
  !   2. Read coefficients
  !   3. Allocate RTTOV input and output structures
  !   4. Set up the chanprof array with the channels/profiles to simulate
  !   5. Read input profile(s)
  !   6. Set up surface emissivity and/or reflectance
  !   7. Call rttov_direct and store results
  !   8. Deallocate all structures and arrays

  ! In this example we supply a profile of the optical parameters for each
  ! channel. Therefore an aerosol coefficient file is not required.

  ! The optical parameters required in the input file are:
  !   - the list of angles for which the phase function is defined
  ! and for each input profile:
  !   - absorption coefficients for all channels for each layer
  !   - scattering coefficients for all channels for each layer
  !   - the phase function over all angles for each channel for each layer

  ! The phase functions are used to calculate the "b parameters" which are
  ! required by the RTTOV scattering parameterisation. This calculation is
  ! relatively slow so you may want to consider calculating the b parameters
  ! off-line and storing them in a file if possible.

  ! The additional steps required when running with optical parameters are:
  !   1. Allocate the optical parameter structure
  !   2. Read the optical parameters from the input file
  !   3. Calculate the b parameters from the phase functions
  !   4. Call rttov_init_opt_params to pre-calculate some phase angle data
  !      (for solar calculations only)

  ! If nthreads is greater than 1 the parallel RTTOV interface is used.
  ! To take advantage of multi-threaded execution you must have compiled
  ! RTTOV with openmp enabled. See the user guide and the compiler flags.

  errorstatus = 0_jpim

  !=====================================================
  !========== Interactive inputs == start ==============

  WRITE(0,*) 'enter name of coefficient file (in current directory)'
  READ(*,*) coef_filename
  WRITE(0,*) 'enter name of file containing profile data (in current directory)'
  READ(*,*) prof_filename
  WRITE(0,*) 'enter name of aerosol optical parameter file (in current directory)'
  READ(*,*) file_aer_opt_param
  WRITE(0,*) 'enter number of angles for which phase function is defined'
  READ(*,*) nphangle
  WRITE(0,*) 'enter number of profiles'
  READ(*,*) nprof
  WRITE(0,*) 'enter number of profile levels'
  READ(*,*) nlevels
  WRITE(0,*) 'turn on solar simulations? (0=no, 1=yes)'
  READ(*,*) dosolar
  WRITE(0,*) 'enter number of channels to simulate per profile'
  READ(*,*) nchannels
  ALLOCATE(channel_list(nchannels))
  WRITE(0,*) 'enter space-separated channel list'
  READ(*,*,iostat=ios) channel_list(:)
  WRITE(0,*) 'enter number of threads to use'
  READ(*,*) nthreads


  ! --------------------------------------------------------------------------
  ! 1. Initialise RTTOV options structure
  ! --------------------------------------------------------------------------

  IF (dosolar == 1) THEN
    opts % rt_ir % addsolar = .TRUE.           ! Include solar radiation
  ELSE
    opts % rt_ir % addsolar = .FALSE.          ! Do not include solar radiation
  ENDIF
  opts % interpolation % addinterp   = .TRUE.  ! Allow interpolation of input profile
  opts % interpolation % interp_mode = 1       ! Set interpolation method
  opts % rt_all % addrefrac          = .TRUE.  ! Include refraction in path calc
  opts % rt_ir % addclouds           = .FALSE. ! Don't include cloud effects

  opts % rt_ir % addaerosl           = .TRUE.  ! Include aerosol effects
  opts % rt_ir % user_aer_opt_param  = .TRUE.  ! Supply optical parameters explictly

  opts % rt_ir % ozone_data          = .FALSE. ! Set the relevant flag to .TRUE.
  opts % rt_ir % co2_data            = .FALSE. !   when supplying a profile of the
  opts % rt_ir % n2o_data            = .FALSE. !   given trace gas (ensure the
  opts % rt_ir % ch4_data            = .FALSE. !   coef file supports the gas)
  opts % rt_ir % co_data             = .FALSE. !
  opts % rt_mw % clw_data            = .FALSE. !

  opts % config % verbose            = .TRUE.  ! Enable printing of warnings

  !========== Interactive inputs == end ==============
  !===================================================


  ! --------------------------------------------------------------------------
  ! 2. Read coefficients
  ! --------------------------------------------------------------------------
  CALL rttov_read_coefs(errorstatus, coefs, opts, file_coef=coef_filename)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'fatal error reading coefficients'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! Ensure input number of channels is not higher than number stored in coefficient file
  IF (nchannels > coefs % coef % fmv_chn) THEN
    nchannels = coefs % coef % fmv_chn
  ENDIF

  ! Ensure the options and coefficients are consistent
  CALL rttov_user_options_checkinput(errorstatus, opts, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'error in rttov options'
    CALL rttov_exit(errorstatus)
  ENDIF


  ! --------------------------------------------------------------------------
  ! 3. Allocate RTTOV input and output structures
  ! --------------------------------------------------------------------------

  ! Determine the total number of radiances to simulate (nchanprof).
  ! In this example we simulate all specified channels for each profile, but
  ! in general one can simulate a different number of channels for each profile.

  nchanprof = nchannels * nprof

  ! Allocate structures for rttov_direct
  CALL rttov_alloc_direct( &
      & errorstatus,                 &
      & 1_jpim,                      &  ! 1 => allocate
      & nprof,                       &
      & nchanprof,                   &
      & nlevels,                     &
      & chanprof,                    &
      & opts,                        &
      & profiles,                    &
      & coefs,                       &
      & transmission,                &
      & radiance,                    &
      & calcemis=calcemis,           &
      & emissivity=emissivity,       &
      & calcrefl=calcrefl,           &
      & reflectance=reflectance,     &
      & aer_nphangle=nphangle,       &
      & aer_opt_param=aer_opt_param, &
      & init=.TRUE._jplm)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'allocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF


  ! --------------------------------------------------------------------------
  ! 4. Build the list of profile/channel indices in chanprof
  ! --------------------------------------------------------------------------

  nch = 0_jpim
  DO j = 1, nprof
    DO jch = 1, nchannels
      nch = nch + 1_jpim
      chanprof(nch)%prof = j
      chanprof(nch)%chan = channel_list(jch)
    ENDDO
  ENDDO


  ! --------------------------------------------------------------------------
  ! 5. Read profile data
  ! --------------------------------------------------------------------------

  !===============================================
  !========== Read profiles == start =============

  OPEN(iup, file=TRIM(prof_filename), status='old', iostat=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'error opening profile file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF
  CALL rttov_skipcommentline(iup, errorstatus)

  ! Read gas units for profiles
  READ(iup,*) profiles(1) % gas_units
  profiles(:) % gas_units = profiles(1) % gas_units
  CALL rttov_skipcommentline(iup, errorstatus)

  ! Loop over all profiles and read data for each one
  DO iprof = 1, nprof

    ! Read pressure (hPa), temp (K), WV, O3 (gas units ppmv or kg/kg - as read above)
    READ(iup,*) profiles(iprof) % p(:)
    CALL rttov_skipcommentline(iup, errorstatus)
    READ(iup,*) profiles(iprof) % t(:)
    CALL rttov_skipcommentline(iup, errorstatus)
    READ(iup,*) profiles(iprof) % q(:)
    CALL rttov_skipcommentline(iup, errorstatus)
    ! Ozone profile is commented out in input profile data
!     READ(iup,*) profiles(iprof) % o3(:)
!     CALL rttov_skipcommentline(iup, errorstatus)

    ! 2 meter air variables
    READ(iup,*) profiles(iprof) % s2m % t, &
              & profiles(iprof) % s2m % q, &
              & profiles(iprof) % s2m % p, &
              & profiles(iprof) % s2m % u, &
              & profiles(iprof) % s2m % v, &
              & profiles(iprof) % s2m % wfetc
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Skin variables
    READ(iup,*) profiles(iprof) % skin % t, &
              & profiles(iprof) % skin % fastem   ! FASTEM only applies to MW instruments
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Surface type and water type
    READ(iup,*) profiles(iprof) % skin % surftype, &
              & profiles(iprof) % skin % watertype
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Elevation, latitude and longitude
    READ(iup,*) profiles(iprof) % elevation, &
              & profiles(iprof) % latitude,  &
              & profiles(iprof) % longitude
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Satellite and solar angles
    READ(iup,*) profiles(iprof) % zenangle,    &
              & profiles(iprof) % azangle,     &
              & profiles(iprof) % sunzenangle, &
              & profiles(iprof) % sunazangle
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Cloud variables for simple cloud scheme, set cfraction to 0. to turn this off (VIS/IR only)
    ! (Ignored for aerosol simulations)
    READ(iup,*) profiles(iprof) % ctp, &
              & profiles(iprof) % cfraction
    CALL rttov_skipcommentline(iup, errorstatus)

  ENDDO
  CLOSE(iup)

  ! --------------------------------------------------------------------------
  ! Read the aerosol optical parameter data for a single atmospheric profile
  ! --------------------------------------------------------------------------

  ! The parameters are supplied per channel so this file is specific to the instrument
  OPEN(iup, file=file_aer_opt_param, status='old', iostat=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'error opening profile file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF
  CALL rttov_skipcommentline(iup, errorstatus)

  ! Read in the angles over which the phase function is defined
  READ(iup,*) aer_opt_param % phangle(:)
  CALL rttov_skipcommentline(iup, errorstatus)

  DO iprof = 1, nprof
    joff = nchannels * (iprof-1) + 1

    ! Read absorption and scattering coefficients for each channel for each layer
    DO i = 1, nlevels-1_jpim
      READ(iup,*) aer_opt_param % ABS(joff:joff+nchannels-1,i)
    ENDDO
    CALL rttov_skipcommentline(iup, errorstatus)

    DO i = 1, nlevels-1_jpim
      READ(iup,*) aer_opt_param % sca(joff:joff+nchannels-1,i)
    ENDDO
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Read in phase function: each full phase function read for each channel for each layer
    DO i = 1, nlevels-1_jpim
      DO j = 1, nchannels
        READ(iup,*) aer_opt_param % pha(joff+j-1,i,:)
      ENDDO
    ENDDO
    CALL rttov_skipcommentline(iup, errorstatus)
  ENDDO
  CLOSE(iup)


  ! The "b parameters" must be calculated from the phase functions (this could be done
  ! off-line and the b parameters could be stored to avoid repeating this calculation
  ! for the same phase functions). NB if the phase function is the same for multiple
  ! layers it does not need re-calculating for every layer.

  WRITE(*,*) 'calculating b parameters...'

  ! Initialise the data for the bpr calculation
  CALL rttov_bpr_init(aer_opt_param % phangle(:), errorstatus)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'error calling rttov_bpr_init'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! Calculate the b parameter for each channel for each layer
  DO iprof = 1, nprof
    DO i = 1, nlevels-1_jpim
      DO j = 1, nchannels
        joff = nchannels * (iprof-1)

        ! Only do calculation for layers containing aerosols
        IF (ANY(aer_opt_param % pha(j+joff,i,:) > 0._jprb)) THEN

          ! If this phase function is same as previous layer copy the bpr
          IF (i > 1 .AND. ALL(aer_opt_param % pha(j+joff,i,:) == aer_opt_param % pha(j+joff,MAX(i-1,1),:))) THEN
            aer_opt_param % bpr(j+joff,i) = aer_opt_param % bpr(j+joff,i-1)
          ELSE
            CALL rttov_bpr_calc(aer_opt_param % pha(j+joff,i,:), &! input phase function
                              & aer_opt_param % phangle(:),      &! input phase angles
                              & aer_opt_param % bpr(j+joff,i),   &! output b parameter
                              & errorstatus)
            IF (errorstatus /= errorstatus_success) THEN
              WRITE(*,*) 'error calculating b parameters'
              CALL rttov_exit(errorstatus)
            ENDIF
          ENDIF

        ELSE
          aer_opt_param % bpr(j+joff,i) = 0._jprb
        ENDIF

      ENDDO
    ENDDO
  ENDDO

  ! Deallocate the data for the bpr calculation
  CALL rttov_bpr_dealloc(errorstatus)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'error calling rttov_bpr_dealloc'
  ENDIF

  WRITE(*,*) 'done calculating b parameters'

  ! If doing solar calculations pre-calculate some phase angle data for scattering calculations
  ! This needs to be re-run if the phase *angles* change (but not phase *functions*)
  IF (opts % rt_ir % addsolar) THEN
    CALL rttov_init_opt_param(errorstatus, opts, aer_opt_param)
    IF (errorstatus /= errorstatus_success) THEN
      WRITE(*,*) 'error initialising aer_opt_param'
      CALL rttov_exit(errorstatus)
    ENDIF
  ENDIF

  !========== READ profiles == end =============
  !=============================================


  ! --------------------------------------------------------------------------
  ! 6. Specify surface emissivity and reflectance
  ! --------------------------------------------------------------------------

  ! In this example we have no values for input emissivities
  emissivity(:) % emis_in = 0._jprb

  ! Calculate emissivity within RTTOV where the input emissivity value is
  ! zero or less (all channels in this case)
  calcemis(:) = (emissivity(:) % emis_in <= 0._jprb)

  ! In this example we have no values for input reflectances
  reflectance(:) % refl_in = 0._jprb

  ! Calculate BRDF within RTTOV where the input BRDF value is zero or less
  ! (all channels in this case)
  calcrefl(:) = (reflectance(:) % refl_in <= 0._jprb)

  ! Use default cloud top BRDF for simple cloud in VIS/NIR channels
  reflectance(:) % refl_cloud_top = 0._jprb


  ! --------------------------------------------------------------------------
  ! 7. Call RTTOV forward model
  ! --------------------------------------------------------------------------
  IF (nthreads <= 1) THEN
    CALL rttov_direct(                   &
          & errorstatus,                 &! out   error flag
          & chanprof,                    &! in    channel and profile index structure
          & opts,                        &! in    options structure
          & profiles,                    &! in    profile array
          & coefs,                       &! in    coefficients structure
          & transmission,                &! inout computed transmittances
          & radiance,                    &! inout computed radiances
          & calcemis      = calcemis,    &! in    flag for internal emissivity calcs
          & emissivity    = emissivity,  &! inout input/output emissivities per channel
          & calcrefl      = calcrefl,    &! in    flag for internal BRDF calcs
          & reflectance   = reflectance, &! inout input/output BRDFs per channel
          & aer_opt_param = aer_opt_param)! in    aerosol optical parameters
  ELSE
    CALL rttov_parallel_direct(           &
          & errorstatus,                  &! out   error flag
          & chanprof,                     &! in    channel and profile index structure
          & opts,                         &! in    options structure
          & profiles,                     &! in    profile array
          & coefs,                        &! in    coefficients structure
          & transmission,                 &! inout computed transmittances
          & radiance,                     &! inout computed radiances
          & calcemis      = calcemis,     &! in    flag for internal emissivity calcs
          & emissivity    = emissivity,   &! inout input/output emissivities per channel
          & calcrefl      = calcrefl,     &! in    flag for internal BRDF calcs
          & reflectance   = reflectance,  &! inout input/output BRDFs per channel
          & aer_opt_param = aer_opt_param,&! in    aerosol optical parameters
          & nthreads      = nthreads)      ! in    number of threads to use
  ENDIF

  IF (errorstatus /= errorstatus_success) THEN
    WRITE (*,*) 'rttov_direct error'
    CALL rttov_exit(errorstatus)
  ENDIF


  ! --- Output the results --------------------------------------------------

  ! Open output file where results are written
  OPEN(ioout, file='output_'//NameOfRoutine//'.dat', status='unknown', form='formatted', iostat=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'error opening the output file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF

  WRITE(ioout,*)' -----------------'
  WRITE(ioout,*)' Instrument ', inst_name(coefs % coef % id_inst)
  WRITE(ioout,*)' -----------------'
  WRITE(ioout,*)' '
  CALL rttov_print_opts(opts, lu=ioout)

  DO iprof = 1, nprof

    joff = (iprof-1_jpim) * nchannels

    WRITE(ioout,*)' '
    WRITE(ioout,*)' Profile ', iprof

    CALL rttov_print_profile(profiles(iprof), lu=ioout)

    WRITE(ioout,777)'CHANNELS PROCESSED FOR SAT ', platform_name(coefs % coef % id_platform), coefs % coef % id_sat
    WRITE(ioout,111) (chanprof(j) % chan, j = 1+joff, nchannels+joff)
    WRITE(ioout,*)' '
    WRITE(ioout,*)'CALCULATED BRIGHTNESS TEMPERATURES (K):'
    WRITE(ioout,222) (radiance % bt(j), j = 1+joff, nchannels+joff)
    IF (opts % rt_ir % addsolar) THEN
      WRITE(ioout,*)' '
      WRITE(ioout,*)'CALCULATED SATELLITE REFLECTANCES (BRF):'
      WRITE(ioout,444) (radiance % refl(j), j = 1+joff, nchannels+joff)
    ENDIF
    WRITE(ioout,*)' '
    WRITE(ioout,*)'CALCULATED RADIANCES (mW/m2/sr/cm-1):'
    WRITE(ioout,222) (radiance % total(j), j = 1+joff, nchannels+joff)
    WRITE(ioout,*)' '
    WRITE(ioout,*)'CALCULATED SURFACE EMISSIVITIES:'
    WRITE(ioout,444) (emissivity(j) % emis_out, j = 1+joff, nchannels+joff)
    IF (opts % rt_ir % addsolar) THEN
      WRITE(ioout,*)' '
      WRITE(ioout,*)'CALCULATED SURFACE BRDF:'
      WRITE(ioout,444) (reflectance(j) % refl_out, j = 1+joff, nchannels+joff)
    ENDIF
  ENDDO

  ! Close output file
  CLOSE(ioout, iostat=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'error closing the output file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF

  ! --- End of output section -----------------------------------------------


  ! --------------------------------------------------------------------------
  ! 8. Deallocate all RTTOV arrays and structures
  ! --------------------------------------------------------------------------
  DEALLOCATE (channel_list, stat=alloc_status)
  IF (alloc_status /= 0) THEN
    WRITE(*,*) 'mem dellocation error'
  ENDIF

  CALL rttov_alloc_direct( &
      & errorstatus,                 &
      & 0_jpim,                      &  ! 0 => deallocate
      & nprof,                       &
      & nchanprof,                   &
      & nlevels,                     &
      & chanprof,                    &
      & opts,                        &
      & profiles,                    &
      & coefs,                       &
      & transmission,                &
      & radiance,                    &
      & calcemis=calcemis,           &
      & emissivity=emissivity,       &
      & calcrefl=calcrefl,           &
      & reflectance=reflectance,     &
      & aer_nphangle=nphangle,       &
      & aer_opt_param=aer_opt_param)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'deallocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF

  CALL rttov_dealloc_coefs(errorstatus, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'coefs deallocation error'
  ENDIF


! Format definitions for output
111  FORMAT(1X,10I8)
222  FORMAT(1X,10F8.2)
444  FORMAT(1X,10F8.3)
777  FORMAT(/,A,A8,I3)

END PROGRAM example_aer_param_fwd
