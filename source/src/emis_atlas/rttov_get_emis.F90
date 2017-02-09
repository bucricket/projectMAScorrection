! Description:
!> @file
!!   Returns emissivity values for the given chanprof array.
!
!> @brief
!!   Returns emissivity values for the given chanprof array.
!!   Optionally, will return estimated errors or, for the
!!   TELSEM MW atlas, covariances.
!!
!! @details
!!    Optional arguments provide access to specific features
!!    of the two atlases. Note that not all arguments apply to
!!    all atlases.
!!
!!    The important elements of the profiles structure for the MW atlases are:
!!    skin\%surftype, latitude, longitude, zenangle.
!!
!!    The important elements of the profiles structure for the IR atlas are:
!!    skin\%surftype, latitude, longitude, snow_frac, zenangle, sunzenangle.
!!    The satellite and solar zenith angles are only required if the angle
!!    correction is being applied and sunzenangle is only used to distinguish
!!    between day (< 85 degrees) and night (> 85 degrees).
!!
!!
!! @param[out]  err          status on exit
!! @param[in]   opts         options to configure the simulations
!! @param[in]   chanprof     specifies channels and profiles to simulate
!! @param[in]   profiles     input atmospheric profiles and surface variables
!! @param[in]   coefs        coefficients structure for instrument to simulate
!! @param[in]   resolution   return emissivities at this resolution, TELSEM MW only,
!!                             units: degrees latitude/longitude, default: 0.25 degrees, optional
!! @param[out]  emissivity   emissivity values
!! @param[out]  emis_std     emissivity errors (standard deviations), IR and TELSEM MW only, optional
!! @param[out]  emis_cov     emissivity covariances, TELSEM MW only, dimensions are (nprof,nchan,nchan)
!!                             where nchan is the largest number of channels simulated per profile, optional
!! @param[out]  emis_flag    emissivity atlas flags, IR only, optional
!! @param[out]  pbats_veg    vegetation type, CNRM MW only, optional
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
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_get_emis(  &
              & err,        & ! out
              & opts,       & ! in
              & chanprof,   & ! in
              & profiles,   & ! in
              & coefs,      & ! in
              & resolution, & ! in, optional (MW atlas only)
              & emissivity, & ! out
              & emis_std,   & ! out, optional (not CNRM atlas)
              & emis_cov,   & ! out, optional (MW atlas only)
              & emis_flag,  & ! out, optional (IR atlas only)
              & pbats_veg)    ! out, optional (MW CNRM atlas only)
!INTF_OFF
#include "throw.h"
!INTF_ON

  USE parkind1, ONLY : jpim, jprb

  USE rttov_types, ONLY : &
          rttov_chanprof, &
          profile_type,   &
          rttov_coefs,    &
          rttov_options

!INTF_OFF
  USE rttov_const, ONLY :      &
          deg2rad,             &
          sensor_id_mw,        &
          sensor_id_po,        &
          surftype_land,       &
          surftype_seaice,     &
          pol_v,               &
          pol_h,               &
          errorstatus_success, &
          errorstatus_fatal

  USE mod_iratlas, ONLY :  &
          rttov_uwiremis,  &
          uw_ir_atlas_version => ir_atlas_version

  USE mod_mwatlas, ONLY :       &
          emis_interp_ind_mult, &
          emis_interp_int_mult, &
          telsem_mw_atlas_version =>  mw_atlas_version

  USE mod_cnrm_mw_atlas, ONLY :     &
          rttov_cnrmmwemis,         &
          cnrm_mw_atlas_version => mw_atlas_version

  USE mod_rttov_emis_atlas, ONLY : &
          ir_atlas_version,        &
          mw_atlas_version,        &
          ir_atlas_init,           &
          mw_atlas_init,           &
          ir_atlas_std_init

  USE yomhook, ONLY : &
          LHOOK,      &
          DR_HOOK
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim),   INTENT(OUT)           :: err
  TYPE(rttov_options),  INTENT(IN)            :: opts
  TYPE(rttov_chanprof), INTENT(IN)            :: chanprof(:)
  TYPE(profile_type),   INTENT(IN)            :: profiles(:)
  TYPE(rttov_coefs),    INTENT(IN)            :: coefs
  REAL(KIND=jprb),      INTENT(IN),  OPTIONAL :: resolution

  REAL(KIND=jprb),      INTENT(OUT)           :: emissivity(size(chanprof))
  REAL(KIND=jprb),      INTENT(OUT), OPTIONAL :: emis_std(size(chanprof))
  ! emis_cov dims are (nprofs, nchans, nchans) where nchans is the maximum # channels per profile
  REAL(KIND=jprb),      INTENT(OUT), OPTIONAL :: emis_cov(:,:,:)
  INTEGER(KIND=jpim),   INTENT(OUT), OPTIONAL :: emis_flag(size(chanprof))
  REAL(KIND=jprb),      INTENT(OUT), OPTIONAL :: pbats_veg(size(chanprof))

!INTF_END


  REAL(KIND=jprb)    :: ZHOOK_HANDLE
  CHARACTER(LEN=100) :: errormsg

  INTEGER(KIND=jpim) :: nchanprof                 ! number of channels/profiles
  INTEGER(KIND=jpim) :: prof                      ! index of current profile
  INTEGER(KIND=jpim) :: nchan                     ! number of channels in current profile
  INTEGER(KIND=jpim) :: chans(coefs%coef%fmv_chn) ! channel indexes for current profile

  INTEGER(KIND=jpim) :: i, j, k                   ! loop indexes
  INTEGER(KIND=jpim) :: lo, hi                    ! low/high chanprof indexes for current profile

  ! Local variables for IR atlas
  INTEGER(KIND=jpim)              :: instr_emis_flag
  REAL(KIND=jprb), ALLOCATABLE    :: instr_emis_cov(:)

  ! Local variables for MW atlas
  REAL(KIND=jprb), ALLOCATABLE    :: ev(:), eh(:)
  REAL(KIND=jprb), ALLOCATABLE    :: std(:,:)

  REAL(KIND=jprb)                 :: sinzen, sinview, sinview_sq, cosview_sq
  INTEGER(KIND=jpim), ALLOCATABLE :: pol_id(:)
  REAL(KIND=jprb)   , ALLOCATABLE :: emissfactor_h(:), emissfactor_v(:)

  REAL(KIND=jprb)  :: lpbats_veg
!----------------------------------------------------------------------------
  TRY

  IF (LHOOK) CALL DR_HOOK('RTTOV_GET_EMIS', 0_jpim, ZHOOK_HANDLE)

  !-----------------------------
  ! Initialise output arguments
  !-----------------------------

  err = errorstatus_success

  IF (PRESENT(emis_std)) emis_std(:) = 0.0_jprb
  IF (PRESENT(emis_cov)) emis_cov(:,:,:) = 0.0_jprb
  IF (PRESENT(emis_flag)) emis_flag(:) = 0_jpim
  IF (PRESENT(pbats_veg)) pbats_veg(:) = 0.0_jprb

  nchanprof = SIZE(chanprof(:))

  !-------------------------------------------
  ! Do array allocations before chanprof loop
  !-------------------------------------------

  ! Note that coefs%coef%fmv_chn is the largest possible number of channels per profile

  ! Also check optional arguments here and report errors when arguments are supplied
  ! which are not supported by the current atlas.

  IF ( coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po ) THEN
    ! MW sensor

    IF (.NOT. mw_atlas_init) THEN
      WRITE(errormsg,'(a)') 'MW atlas not initialised'
      err = errorstatus_fatal
      THROWM(err .NE. 0, errormsg)
    END IF

    IF (mw_atlas_version == telsem_mw_atlas_version) THEN

      ALLOCATE(ev(coefs%coef%fmv_chn),stat=err)
      THROWM(err .NE. 0, 'Allocation of ev')
      ALLOCATE(eh(coefs%coef%fmv_chn),stat=err)
      THROWM(err .NE. 0, 'Allocation of eh')

      IF ( PRESENT(emis_std) .OR. PRESENT(emis_cov) ) THEN
        ! Emissivity errors or covariances requested
        ALLOCATE(std(2*coefs%coef%fmv_chn,2*coefs%coef%fmv_chn),stat=err)
        THROWM(err .NE. 0, 'Allocation of std')
      END IF

      ALLOCATE(pol_id(coefs%coef%fmv_chn),stat=err)
      THROWM(err .NE. 0, 'Allocation of pol_id')
      ALLOCATE(emissfactor_v(coefs%coef%fmv_chn),stat=err)
      THROWM(err .NE. 0, 'Allocation of emissfactor_v')
      ALLOCATE(emissfactor_h(coefs%coef%fmv_chn),stat=err)
      THROWM(err .NE. 0, 'Allocation of emissfactor_h')

      IF ( PRESENT(emis_flag) .AND. opts%config%verbose ) THEN
        WARN('emis_flag argument not supported for TELSEM atlas')
      END IF
      IF ( PRESENT(pbats_veg) .AND. opts%config%verbose  ) THEN
        WARN('pbats_veg argument only supported for TELSEM atlas')
      END IF

    ELSEIF (mw_atlas_version == cnrm_mw_atlas_version) THEN

      IF ( PRESENT(resolution) .AND. opts%config%verbose ) THEN
        WARN('resolution argument not supported for CNRM MW atlas')
      END IF
      IF ( PRESENT(emis_std) .AND. opts%config%verbose ) THEN
        WARN('emis_std argument not supported for CNRM MW atlas')
      END IF
      IF ( PRESENT(emis_cov) .AND. opts%config%verbose ) THEN
        WARN('emis_cov argument not supported for CNRM MW atlas')
      END IF
      IF ( PRESENT(emis_flag) .AND. opts%config%verbose ) THEN
        WARN('emis_flag argument not supported for CNRM MW atlas')
      END IF

    ELSE
      WRITE(errormsg,'(a,i5)') 'Unknown MW atlas version: ', mw_atlas_version
      err = errorstatus_fatal
      THROWM(err .NE. 0, errormsg)
    END IF

  ELSE
    ! IR sensor

    IF (.NOT. ir_atlas_init) THEN
      WRITE(errormsg,'(a)') 'IR atlas not initialised'
      err = errorstatus_fatal
      THROWM(err .NE. 0, errormsg)
    END IF

    IF (ir_atlas_version == uw_ir_atlas_version) THEN

      ALLOCATE(instr_emis_cov(coefs%coef%fmv_chn),stat=err)
      THROWM(err .NE. 0, 'Allocation of instr_emis_cov')

      IF (.NOT. ir_atlas_std_init .AND. PRESENT(emis_std) .AND. opts%config%verbose) THEN
        WARN('IR atlas stddev not read, but emis_std argument present')
      END IF
      IF ( PRESENT(resolution) .AND. opts%config%verbose ) THEN
        WARN('resolution argument not supported for IR atlas')
      END IF
      IF ( PRESENT(emis_cov) .AND. opts%config%verbose ) THEN
        WARN('emis_cov argument not supported for IR atlas')
      END IF
      IF ( PRESENT(pbats_veg) .AND. opts%config%verbose ) THEN
        WARN('pbats_veg argument not supported for IR atlas')
      END IF

    ELSE
      WRITE(errormsg,'(a,i5)') 'Unknown IR atlas version: ', ir_atlas_version
      err = errorstatus_fatal
      THROWM(err .NE. 0, errormsg)
    END IF
  END IF

  !-----------
  ! Main loop
  !-----------

  ! Assume chanprof(:) lists all channels for profile 1, followed by all channels
  ! for profile 2, and so on. Loop through chanprof(:) until the profile number
  ! changes and in this way obtain the list of channels for the current profile.
  k = 1
  DO WHILE ( k <= nchanprof )
    prof = chanprof(k)%prof
    lo = k
    k  = lo+1
    DO
      IF ( k > nchanprof ) EXIT
      IF ( prof /= chanprof(k)%prof ) EXIT
      k = k+1
    END DO
    hi = k-1
    nchan = k-lo
    chans(1:nchan) = chanprof(lo:hi)%chan

    ! Profile index is 'prof'
    ! Number of channels for this profile is 'nchan'
    ! Channel indexes are 'chans(1:nchan)'

    IF ( coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po ) THEN
    !---------------------------
    ! MW atlas
    !---------------------------
      IF ( mw_atlas_version == telsem_mw_atlas_version ) THEN

        IF ( profiles(prof)%skin%surftype == surftype_land ) THEN

          IF ( PRESENT(emis_std) .OR. PRESENT(emis_cov) ) THEN
            ! Emissivity errors or covariances requested

            IF ( PRESENT(resolution) ) THEN
              ! User resolution
              CALL emis_interp_int_mult(                      &
                  & profiles(prof)%latitude,                  &! in
                  & profiles(prof)%longitude,                 &! in
                  & resolution,                               &! in
                  & profiles(prof)%zenangle,                  &! in
                  & coefs%coef%frequency_ghz(chans(1:nchan)), &! in
                  & nchan,                                    &! in
                  & ev(1:nchan),                              &! out
                  & eh(1:nchan),                              &! out
                  & std=std(1:(2*nchan),1:(2*nchan)),         &! out
                  & verb=0_jpim                               )! in
            ELSE
              ! Atlas resolution
              CALL emis_interp_ind_mult(                      &
                  & profiles(prof)%latitude,                  &! in
                  & profiles(prof)%longitude,                 &! in
                  & profiles(prof)%zenangle,                  &! in
                  & coefs%coef%frequency_ghz(chans(1:nchan)), &! in
                  & nchan,                                    &! in
                  & ev(1:nchan),                              &! out
                  & eh(1:nchan),                              &! out
                  & std=std(1:(2*nchan),1:(2*nchan)),         &! out
                  & verb=0_jpim                               )! in
            END IF

          ELSE
            ! No need to supply the std argument

            IF ( PRESENT(resolution) ) THEN
              ! User resolution
              CALL emis_interp_int_mult(                      &
                  & profiles(prof)%latitude,                  &! in
                  & profiles(prof)%longitude,                 &! in
                  & resolution,                               &! in
                  & profiles(prof)%zenangle,                  &! in
                  & coefs%coef%frequency_ghz(chans(1:nchan)), &! in
                  & nchan,                                    &! in
                  & ev(1:nchan),                              &! out
                  & eh(1:nchan),                              &! out
                  & verb=0_jpim                               )! in
            ELSE
              ! Atlas resolution
              CALL emis_interp_ind_mult(                      &
                  & profiles(prof)%latitude,                  &! in
                  & profiles(prof)%longitude,                 &! in
                  & profiles(prof)%zenangle,                  &! in
                  & coefs%coef%frequency_ghz(chans(1:nchan)), &! in
                  & nchan,                                    &! in
                  & ev(1:nchan),                              &! out
                  & eh(1:nchan),                              &! out
                  & verb=0_jpim                               )! in
            END IF
          END IF

          !------------------------------------------------------
          ! Combine H- and V-pol emissivities into output arrays
          !------------------------------------------------------

          sinzen     = SIN(profiles(prof)%zenangle * deg2rad)
          sinview    = sinzen / coefs%coef%ratoe
          sinview_sq = sinview * sinview
          cosview_sq = 1.0_jprb - sinview_sq

          pol_id(1:nchan) = coefs%coef%fastem_polar(chans(1:nchan)) + 1_jpim

          DO j = 1, nchan
            emissfactor_v(j) = pol_v(1,pol_id(j)) + pol_v(2,pol_id(j)) * sinview_sq + pol_v(3,pol_id(j)) * cosview_sq
            emissfactor_h(j) = pol_h(1,pol_id(j)) + pol_h(2,pol_id(j)) * sinview_sq + pol_h(3,pol_id(j)) * cosview_sq
          END DO
          emissivity(lo:hi) = ev(1:nchan)*emissfactor_v(1:nchan) + eh(1:nchan)*emissfactor_h(1:nchan)

          IF ( PRESENT(emis_std) ) THEN
            ! Generate the stdv

            ! Each dimension of std(:,:) has V-pol values for all channels followed by H-pol values
            ! for all channels.
            ! In this case std(:,:) contains covariances, so diagonal elements are variances (not stddev)
            ! but we want output emis_std(:) to be stddev

            DO j = 1, nchan
              emis_std(lo+j-1) = SQRT(emissfactor_v(j) * emissfactor_v(j) * std(j,j)            + &
                                      emissfactor_v(j) * emissfactor_h(j) * 2 * std(j,j+nchan)  + &
                                      emissfactor_h(j) * emissfactor_h(j) * std(j+nchan,j+nchan))
            END DO
          END IF

          IF ( PRESENT(emis_cov) ) THEN
            ! Generate the covariances

            ! std(:,:) contains covariances. Combine the H- and V-pol values for each pair of
            ! channels to find the covariance of the combined emissivities in the two channels.

            DO i = 1, nchan
              DO j = 1, nchan
                emis_cov(prof,i,j) = emissfactor_v(i) * emissfactor_v(j) * std(i,       j      ) + &
                                     emissfactor_v(i) * emissfactor_h(j) * std(i,       j+nchan) + &
                                     emissfactor_h(i) * emissfactor_v(j) * std(i+nchan, j      ) + &
                                     emissfactor_h(i) * emissfactor_h(j) * std(i+nchan, j+nchan)
              END DO
            END DO
          END IF

        ELSE
          ! Not land
          emissivity(lo:hi) = 0.0_jprb
        END IF

      ELSE IF (mw_atlas_version == cnrm_mw_atlas_version) THEN

        IF ( profiles(prof)%skin%surftype == surftype_land ) THEN

          CALL rttov_cnrmmwemis(                         &
                & err,                                   &! out
                & coefs%coef%id_inst,                    &! in
                & nchan,                                 &! in
                & profiles(prof)%latitude,               &! in
                & profiles(prof)%longitude,              &! in
                & profiles(prof)%zenangle,               &! in
                & coefs%coef%ff_ori_chn(chans(1:nchan)), &! in
                & emissivity(lo:hi),                     &! out
                & lpbats_veg)                             ! out
          THROWM(err .NE. 0, "error in CNRM atlas")
          IF (PRESENT(pbats_veg)) pbats_veg(lo:hi) = lpbats_veg

        ELSE
          ! Not land
          emissivity(lo:hi) = 0.0_jprb
        END IF

      ELSE
        WRITE(errormsg,'(a,i5)') 'Unknown MW atlas version: ', mw_atlas_version
        err = errorstatus_fatal
        THROWM(err .NE. 0, errormsg)
      END IF

    ELSE
    !---------------------------
    ! IR atlas
    !---------------------------
      IF (ir_atlas_version == uw_ir_atlas_version) THEN

        IF ( profiles(prof)%skin%surftype == surftype_land .OR. &
             profiles(prof)%skin%surftype == surftype_seaice    ) THEN

          CALL rttov_uwiremis(                        &
                & opts%config%verbose,                &! in
                & nchan,                              &! in
                & profiles(prof)%latitude,            &! in
                & profiles(prof)%longitude,           &! in
                & profiles(prof)%zenangle,            &! in
                & profiles(prof)%sunzenangle,         &! in
                & profiles(prof)%skin%surftype,       &! in
                & profiles(prof)%snow_frac,           &! in
                & coefs%coef%ff_cwn(chans(1:nchan)),  &! in
                & chans(1:nchan),                     &! in
                & emissivity(lo:hi),                  &! out
                & instr_emis_cov(1:nchan),            &! out
                & instr_emis_flag)                     ! out

          IF ( PRESENT(emis_std)  ) emis_std(lo:hi)  = instr_emis_cov(1:nchan)
          IF ( PRESENT(emis_flag) ) emis_flag(lo:hi) = instr_emis_flag

        ELSE
          ! Not land or seaice
          emissivity(lo:hi) = 0.0_jprb
        END IF

      END IF

    END IF ! coefs%coef%id_sensor

  END DO ! chanprof


  !------------------
  ! Deallocate arrays
  !------------------
  IF ( coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po ) THEN
    ! MW sensor
    IF (mw_atlas_version == telsem_mw_atlas_version) THEN
      IF ( ALLOCATED(ev) ) DEALLOCATE(ev)
      IF ( ALLOCATED(eh) ) DEALLOCATE(eh)

      IF ( ALLOCATED(std)   ) DEALLOCATE(std)

      IF ( ALLOCATED(pol_id)        ) DEALLOCATE(pol_id)
      IF ( ALLOCATED(emissfactor_v) ) DEALLOCATE(emissfactor_v)
      IF ( ALLOCATED(emissfactor_h) ) DEALLOCATE(emissfactor_h)

    ELSE  IF (mw_atlas_version == cnrm_mw_atlas_version) THEN
    END IF

  ELSE
    ! IR sensor

    IF ( ALLOCATED(instr_emis_cov) ) DEALLOCATE(instr_emis_cov)

  END IF

  IF (LHOOK) CALL DR_HOOK('RTTOV_GET_EMIS',1_jpim,ZHOOK_HANDLE)

  CATCH

  IF (LHOOK) CALL DR_HOOK('RTTOV_GET_EMIS',1_jpim,ZHOOK_HANDLE)

END SUBROUTINE rttov_get_emis
