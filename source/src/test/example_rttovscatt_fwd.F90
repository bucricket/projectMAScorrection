PROGRAM example_rttovscatt_fwd
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
  !    Copyright 2010, EUMETSAT, All Rights Reserved.
  !
  !     *************************************************************
  !
  !     TEST PROGRAM FOR RTTOV-SCATT FORWARD MODEL
  !          RTTOV VERSION 11
  !
  ! To run this program you must have the following files
  ! either resident in the same directory or set up as a
  ! symbolic link:
  !   the file containing input profiles (e.g. prof.dat)
  !   the file containing input cloud data
  !   the RTTOV coefficient file
  !   the RTTOV Mietable file
  !
  ! The script run_example_rttovscatt_fwd.sh may be used to run this program.
  ! The output is generated in a file called example_rttovscatt_fwd_output.dat.
  !
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
       & rttov_options_scatt, &
       & rttov_coefs,         &
       & rttov_scatt_coef,    &
       & profile_type,        &
       & profile_cloud_type,  &
       & transmission_type,   &
       & radiance_type,       &
       & rttov_chanprof,      &
       & rttov_emissivity

  ! jpim, jprb and jplm are the RTTOV integer, real and logical KINDs
  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_unix_env, ONLY : rttov_exit

  IMPLICIT NONE

#include "rttov_scatt.interface"
#include "rttov_alloc_scatt_prof.interface"
#include "rttov_read_scattcoeffs.interface"
#include "rttov_dealloc_scattcoeffs.interface"
#include "rttov_scatt_setupindex.interface"

#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_skipcommentline.interface"

  !--------------------------
  !
  INTEGER(KIND=jpim), PARAMETER :: iup   = 20   ! unit for input profile file
  INTEGER(KIND=jpim), PARAMETER :: ioout = 21   ! unit for output

  ! RTTOV variables/structures
  !====================
  TYPE(rttov_options)              :: opts                     ! Options structure - leave this set to defaults
  TYPE(rttov_options_scatt)        :: opts_scatt               ! RTTOV-SCATT options structure
  TYPE(rttov_coefs)                :: coefs                    ! Coefficients structure
  TYPE(rttov_scatt_coef)           :: coefs_scatt              ! RTTOV-SCATT coefficients structure
  TYPE(rttov_chanprof),    POINTER :: chanprof(:)    => NULL() ! Input channel/profile list
  INTEGER(KIND=jpim),      POINTER :: frequencies(:) => NULL() ! Channel indexes for Mietable lookup
  LOGICAL(KIND=jplm),      POINTER :: use_chan(:,:)  => NULL() ! Flags to specify channels to simulate
  LOGICAL(KIND=jplm),      POINTER :: calcemis(:)    => NULL() ! Flag to indicate calculation of emissivity within RTTOV
  TYPE(rttov_emissivity),  POINTER :: emissivity(:)  => NULL() ! Input/output surface emissivity
  TYPE(profile_type),      POINTER :: profiles(:)    => NULL() ! Input profiles
  TYPE(profile_cloud_type),POINTER :: cld_profiles(:)=> NULL() ! Input RTTOV-SCATT cloud/hydrometeor profiles
  TYPE(transmission_type)          :: transmission             ! Output transmittances
  TYPE(radiance_type)              :: radiance                 ! Output radiances

  INTEGER(KIND=jpim)               :: errorstatus              ! Return error status of RTTOV subroutine calls

  INTEGER(KIND=jpim) :: alloc_status
  CHARACTER(LEN=22)  :: NameOfRoutine = 'example_rttovscatt_fwd'

  ! variables for input
  !====================
  CHARACTER(LEN=256) :: coef_filename
  CHARACTER(LEN=256) :: prof_filename
  INTEGER(KIND=jpim) :: totalice, snowrain_units
  LOGICAL(KIND=jplm) :: use_totalice, mmr_snowrain
  INTEGER(KIND=jpim) :: nlevels
  INTEGER(KIND=jpim) :: nprof
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: nchanprof
  INTEGER(KIND=jpim), ALLOCATABLE :: channel_list(:)
  ! loop variables
  INTEGER(KIND=jpim) :: j
  INTEGER(KIND=jpim) :: ilev, nprint
  INTEGER(KIND=jpim) :: iprof, joff
  INTEGER            :: ios

  !- End of header --------------------------------------------------------

  ! The usual steps to take when running RTTOV-SCATT are as follows:
  !   1. Specify required RTTOV-SCATT options
  !   2. Read coefficients and Mietable file
  !   3. Allocate RTTOV input and output structures
  !   4. Set up the chanprof and frequencies arrays by calling rttov_scatt_setupindex
  !   5. Read input profile(s)
  !   6. Set up surface emissivity
  !   7. Call rttov_scatt and store results
  !   8. Deallocate all structures and arrays

  errorstatus = 0_jpim

  !=====================================================
  !========== Interactive inputs == start ==============

  WRITE(0,*) 'enter name of coefficient file (in current directory)'
  READ(*,*) coef_filename
  WRITE(0,*) 'enter name of file containing profile data (in current directory)'
  READ(*,*) prof_filename
  WRITE(0,*) 'enter number of profiles'
  READ(*,*) nprof
  WRITE(0,*) 'enter number of profile levels'
  READ(*,*) nlevels
  WRITE(0,*) 'use totalice? (0=no, 1=yes)'
  READ(*,*) totalice
  WRITE(0,*) 'snow/rain units? (0=kg/m2/s, 1=kg/kg)'
  READ(*,*) snowrain_units
  WRITE(0,*) 'enter number of channels to simulate per profile'
  READ(*,*) nchannels
  ALLOCATE(channel_list(nchannels))
  WRITE(0,*) 'enter space-separated channel list'
  READ(*,*,iostat=ios) channel_list(:)

  use_totalice = (totalice /= 0_jpim)
  mmr_snowrain = (snowrain_units /= 0_jpim)

  ! --------------------------------------------------------------------------
  ! 1. Initialise RTTOV-SCATT options structure
  ! --------------------------------------------------------------------------

  ! The rttov_options structure (opts) should be left with its default values.
  ! RTTOV-SCATT only allows access to a limited number of RTTOV options: these
  ! are set in the rttov_options_scatt structure (opts_scatt).

  ! For example:
  opts_scatt % interp_mode = 1                    ! Set interpolation method
  opts_scatt % config % verbose = .TRUE.          ! Enable printing of warnings

  ! See user guide for full list of RTTOV-SCATT options


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

  ! Read the RTTOV-SCATT Mietable file
  CALL rttov_read_scattcoeffs(errorstatus, coefs%coef, coefs_scatt)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'fatal error reading RTTOV-SCATT coefficients'
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
      & errorstatus,             &
      & 1_jpim,                  &  ! 1 => allocate
      & nprof,                   &
      & nchanprof,               &
      & nlevels,                 &
      & chanprof,                &
      & opts,                    &
      & profiles,                &
      & coefs,                   &
      & transmission,            &
      & radiance,                &
      & calcemis=calcemis,       &
      & emissivity=emissivity,   &
      & init=.TRUE._jplm)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'allocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! Allocate the RTTOV-SCATT cloud profiles structure
  ALLOCATE(cld_profiles(nprof), stat=alloc_status)
  IF (alloc_status /= 0) THEN
    WRITE(*,*) 'allocation error for cld_profiles'
    errorstatus = errorstatus_fatal
    CALL rttov_exit(errorstatus)
  ENDIF

  CALL rttov_alloc_scatt_prof(   &
      & nprof,                   &
      & cld_profiles,            &
      & nlevels,                 &
      & use_totalice,            &    ! false => separate ciw and snow; true => totalice
      & 1_jpim,                  &    ! 1 => allocate
      & init = .TRUE._jplm,      &
      & mmr_snowrain = mmr_snowrain)  ! snow/rain input units: false => kg/m2/s; true => kg/kg


  ! --------------------------------------------------------------------------
  ! 4. Populate chanprof and frequencies arrays
  ! --------------------------------------------------------------------------

  ! RTTOV-SCATT requires the frequencies array to be populated by a call to
  ! rttov_scatt_setupindex. This also populates the chanprof array. To specify
  ! only a subset of channels (i.e. those in channel_list) an array of flags is
  ! passed in (use_chan).

  ! use_chan array is dimensioned by the total number of instrument channels
  ALLOCATE(use_chan(nprof,coefs%coef%fmv_chn), &
           frequencies(nchanprof))

  ! Set use_chan to .TRUE. only for required channels
  use_chan(:,:) = .FALSE._jplm
  DO j = 1, nprof
    use_chan(j,channel_list(1:nchannels)) = .TRUE._jplm
  ENDDO

  ! Populate chanprof and frequencies arrays
  CALL rttov_scatt_setupindex ( &
      & nprof,              &
      & coefs%coef%fmv_chn, &
      & coefs%coef,         &
      & nchanprof,          &
      & chanprof,           &
      & frequencies,        &
      & use_chan)


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

    ! Read vertical profile data
    ! NB The bottom-most half pressure level is taken as the 2m pressure (see below)
    DO ilev = 1, nlevels
      IF (use_totalice) THEN
        READ(iup,*) &
          & profiles    (iprof) % p   (ilev),     &   ! full level pressure (hPa)
          & cld_profiles(iprof) % ph  (ilev),     &   ! half level pressure (hPa)
          & profiles    (iprof) % t   (ilev),     &   ! temperature (K)
          & profiles    (iprof) % q   (ilev),     &   ! specific humidity (ppmv or kg/kg - as read above)
          & cld_profiles(iprof) % cc  (ilev),     &   ! cloud cover (0-1)
          & cld_profiles(iprof) % clw (ilev),     &   ! liquid water (kg/kg)
          & cld_profiles(iprof) % totalice(ilev), &   ! combined ice water and snow (kg/kg)
          & cld_profiles(iprof) % rain(ilev)          ! rain (kg/kg)
      ELSE
        READ(iup,*) &
          & profiles    (iprof) % p   (ilev), &   ! full level pressure (hPa)
          & cld_profiles(iprof) % ph  (ilev), &   ! half level pressure (hPa)
          & profiles    (iprof) % t   (ilev), &   ! temperature (K)
          & profiles    (iprof) % q   (ilev), &   ! specific humidity (ppmv or kg/kg - as read above)
          & cld_profiles(iprof) % cc  (ilev), &   ! cloud cover (0-1)
          & cld_profiles(iprof) % clw (ilev), &   ! liquid water (kg/kg)
          & cld_profiles(iprof) % ciw (ilev), &   ! ice water (kg/kg)
          & cld_profiles(iprof) % rain(ilev), &   ! rain (kg/kg)
          & cld_profiles(iprof) % sp  (ilev)      ! frozen precip. (kg/kg)
      ENDIF
    ENDDO
    CALL rttov_skipcommentline(iup, errorstatus)

    ! 2 meter air variables
    READ(iup,*) profiles(iprof) % s2m % t, &
              & profiles(iprof) % s2m % q, &
              & profiles(iprof) % s2m % p, &
              & profiles(iprof) % s2m % u, &
              & profiles(iprof) % s2m % v
    CALL rttov_skipcommentline(iup, errorstatus)

    ! The bottom-most half pressure level is taken as the 2m pressure
    cld_profiles(iprof) % ph(nlevels+1) = profiles(iprof) % s2m % p

    ! Skin variables
    READ(iup,*) profiles(iprof) % skin % t, &
              & profiles(iprof) % skin % fastem
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

    ! Satellite angles
    READ(iup,*) profiles(iprof) % zenangle, &
              & profiles(iprof) % azangle
    CALL rttov_skipcommentline(iup, errorstatus)

  ENDDO
  CLOSE(iup)

  !========== READ profiles == end =============
  !=============================================


  ! --------------------------------------------------------------------------
  ! 6. Specify surface emissivity
  ! --------------------------------------------------------------------------

  ! In this example we have no values for input emissivities
  emissivity(:) % emis_in = 0._jprb

  ! Calculate emissivity within RTTOV where the input emissivity value is
  ! zero or less (all channels in this case)
  calcemis(:) = (emissivity(:) % emis_in <= 0._jprb)


  ! --------------------------------------------------------------------------
  ! 7. Call RTTOV-SCATT forward model
  ! --------------------------------------------------------------------------
  CALL rttov_scatt ( &
    & errorstatus,         &! out   error flag
    & opts_scatt,          &! in    RTTOV-SCATT options structure
    & nlevels,             &! in    number of profile levels
    & chanprof,            &! in    channel and profile index structure
    & frequencies,         &! in    channel indexes for Mietable lookup
    & profiles,            &! in    profile array
    & cld_profiles,        &! in    cloud/hydrometeor profile array
    & coefs,               &! in    coefficients structure
    & coefs_scatt,         &! in    Mietable structure
    & calcemis,            &! in    flag for internal emissivity calcs
    & emissivity,          &! inout input/output emissivities per channel
    & radiance)             ! inout computed radiances

  IF (errorstatus /= errorstatus_success) THEN
    WRITE (*,*) 'rttov_scatt error'
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

  DO iprof = 1, nprof

    joff = (iprof-1_jpim) * nchannels

    nprint = 1 + INT((nchannels-1)/10)
    WRITE(ioout,*)' '
    WRITE(ioout,*)' Profile ', iprof

    WRITE(ioout,777)'CHANNELS PROCESSED FOR SAT ', platform_name(coefs % coef % id_platform), coefs % coef % id_sat
    WRITE(ioout,111) (chanprof(j) % chan, j = 1+joff, nchannels+joff)
    WRITE(ioout,*)' '
    WRITE(ioout,*)'CALCULATED BRIGHTNESS TEMPERATURES (K):'
    WRITE(ioout,444) (radiance % bt(j), j = 1+joff, nchannels+joff)
    WRITE(ioout,*)' '
    WRITE(ioout,*)'CALCULATED SURFACE EMISSIVITIES:'
    WRITE(ioout,444) (emissivity(j) % emis_out, j = 1+joff, nchannels+joff)
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
  DEALLOCATE(channel_list, use_chan, frequencies, stat=alloc_status)
  IF (alloc_status /= 0) THEN
    WRITE(*,*) 'mem dellocation error'
  ENDIF

  CALL rttov_alloc_scatt_prof(   &
      & nprof,                   &
      & cld_profiles,            &
      & nlevels,                 &
      & use_totalice,            &  ! must match value used for allocation
      & 0_jpim)                     ! 0 => deallocate

  DEALLOCATE(cld_profiles, stat=alloc_status)
  IF (alloc_status /= 0) THEN
    WRITE(*,*) 'dellocation error for cld_profiles'
  ENDIF

  ! Deallocate structures for rttov_direct
  CALL rttov_alloc_direct( &
      & errorstatus,             &
      & 0_jpim,                  &  ! 0 => deallocate
      & nprof,                   &
      & nchanprof,               &
      & nlevels,                 &
      & chanprof,                &
      & opts,                    &
      & profiles,                &
      & coefs,                   &
      & transmission,            &
      & radiance,                &
      & calcemis=calcemis,       &
      & emissivity=emissivity)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'deallocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF

  CALL rttov_dealloc_scattcoeffs(coefs_scatt)

  CALL rttov_dealloc_coefs(errorstatus, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'coefs deallocation error'
  ENDIF


! Format definitions for output
111  FORMAT(1X,10I8)
444  FORMAT(1X,10F8.3)
777  FORMAT(/,A,A8,I3)

END PROGRAM example_rttovscatt_fwd
