! Description:
!> @file
!!   Returns BRDF values for the given chanprof array.
!
!> @brief
!!   Returns BRDF values for the given chanprof array.
!!   Optionally, will return bi-hemispherical (back-sky) albedo.
!!
!! @details
!!    The important elements of the profiles structure for the BRDF atlas are:
!!    skin\%surftype, skin\%watertype, latitude, longitude,
!!    zenangle, azangle, sunzenangle, sunazangle.
!!
!! @param[out]  err          status on exit
!! @param[in]   chanprof     specifies channels and profiles to simulate
!! @param[in]   profiles     input atmospheric profiles and surface variables
!! @param[in]   coefs        coefficients structure for instrument to simulate
!! @param[out]  brdf         BRDF values
!! @param[out]  brdf_flag    BRDF atlas flags, optional
!! @param[out]  bh_albedo    Bi-hemispherical (back-sky) albedo, optional
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
SUBROUTINE rttov_get_brdf(  &
              & err,        & ! out
              & chanprof,   & ! in
              & profiles,   & ! in
              & coefs,      & ! in
              & brdf,       & ! out
              & brdf_flag,  & ! out, optional
              & bh_albedo)    ! out, optional
!INTF_OFF
#include "throw.h"
!INTF_ON

  USE parkind1, ONLY : jpim, jprb

  USE rttov_types, ONLY : &
          rttov_chanprof, &
          profile_type,   &
          rttov_coefs

!INTF_OFF
  USE rttov_const, ONLY :      &
          surftype_seaice,     &
          errorstatus_success, &
          errorstatus_fatal

  USE mod_brdf_atlas, ONLY :      &
        & rttov_visnirbrdf, &
        & cms_vn_atlas_version => vn_atlas_version

  USE mod_rttov_brdf_atlas, ONLY : &
        & vn_atlas_version,        &
        & vn_atlas_init

  USE yomhook, ONLY : &
          LHOOK,      &
          DR_HOOK
!INTF_ON

  IMPLICIT NONE

  INTEGER(KIND=jpim),   INTENT(OUT)           :: err
  TYPE(rttov_chanprof), INTENT(IN)            :: chanprof(:)
  TYPE(profile_type),   INTENT(IN)            :: profiles(:)
  TYPE(rttov_coefs),    INTENT(IN)            :: coefs

  REAL(KIND=jprb),      INTENT(OUT)           :: brdf(size(chanprof))
  INTEGER(KIND=jpim),   INTENT(OUT), OPTIONAL :: brdf_flag(size(chanprof))
  REAL(KIND=jprb),      INTENT(OUT), OPTIONAL :: bh_albedo(size(chanprof))

!INTF_END


  REAL(KIND=jprb)    :: ZHOOK_HANDLE
  CHARACTER(LEN=100) :: errormsg

  INTEGER(KIND=jpim) :: nchanprof                 ! number of channels/profiles
  INTEGER(KIND=jpim) :: prof                      ! index of current profile
  INTEGER(KIND=jpim) :: nchan                     ! number of channels in current profile
  INTEGER(KIND=jpim) :: chans(coefs%coef%fmv_chn) ! channel indexes for current profile

  INTEGER(KIND=jpim) :: k                         ! loop indexes
  INTEGER(KIND=jpim) :: lo, hi                    ! low/high chanprof indexes for current profile

  INTEGER(KIND=jpim) :: instr_brdf_flag

!----------------------------------------------------------------------------
  TRY

  IF (LHOOK) CALL DR_HOOK('RTTOV_GET_BRDF', 0_jpim, ZHOOK_HANDLE)

  !-----------------------------
  ! Initialise output arguments
  !-----------------------------

  err = errorstatus_success

  IF (PRESENT(brdf_flag)) brdf_flag(:) = 0_jpim
  IF (PRESENT(bh_albedo)) bh_albedo(:) = 0._JPRB

  nchanprof = SIZE(chanprof(:))

  ! VIS/NIR sensor

  IF (.NOT. vn_atlas_init) THEN
    WRITE(errormsg,'(a)') 'vis/nir atlas not initialised'
    err = errorstatus_fatal
    THROWM(err .NE. 0, errormsg)
  END IF

  IF (vn_atlas_version == cms_vn_atlas_version) THEN
    !! need to be corrected

  ELSE
    WRITE(errormsg,'(a,i5)') 'Unknown vis/nir atlas version: ', vn_atlas_version
    err = errorstatus_fatal
    THROWM(err .NE. 0, errormsg)
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

    !---------------------------
    ! VIS/NIR atlas
    !---------------------------
    IF (vn_atlas_version == cms_vn_atlas_version) THEN

      IF ( profiles(prof)%skin%surftype == surftype_seaice ) THEN
!         WRITE(errormsg,'(a)') 'No BRDF data for sea ice surface'
!         err = errorstatus_fatal
!         THROWM(err .NE. 0, errormsg)
        brdf(lo:hi) = -999.0_jprb
      ELSE
        IF ( PRESENT(bh_albedo) ) THEN
          CALL rttov_visnirbrdf(                        &
                  & nchan,                              &! in
                  & profiles(prof)%latitude,            &! in
                  & profiles(prof)%longitude,           &! in
                  & profiles(prof)%skin%surftype,       &! in
                  & profiles(prof)%skin%watertype,      &! in
                  & profiles(prof)%zenangle,            &! in
                  & profiles(prof)%azangle,             &! in
                  & profiles(prof)%sunzenangle,         &! in
                  & profiles(prof)%sunazangle,          &! in
                  & coefs%coef%ff_cwn(chans(1:nchan)),  &! in
                  & chans(1:nchan),                     &! in
                  & brdf(lo:hi),                        &! out
                  & instr_brdf_flag,                    &! out
                  & bh_albedo(lo:hi))                    ! out
        ELSE
          CALL rttov_visnirbrdf(                        &
                  & nchan,                              &! in
                  & profiles(prof)%latitude,            &! in
                  & profiles(prof)%longitude,           &! in
                  & profiles(prof)%skin%surftype,       &! in
                  & profiles(prof)%skin%watertype,      &! in
                  & profiles(prof)%zenangle,            &! in
                  & profiles(prof)%azangle,             &! in
                  & profiles(prof)%sunzenangle,         &! in
                  & profiles(prof)%sunazangle,          &! in
                  & coefs%coef%ff_cwn(chans(1:nchan)),  &! in
                  & chans(1:nchan),                     &! in
                  & brdf(lo:hi),                        &! out
                  & instr_brdf_flag)                     ! out
        END IF

        IF ( PRESENT(brdf_flag) ) brdf_flag(lo:hi) = instr_brdf_flag

      END IF

    END IF

  END DO ! chanprof

  IF (LHOOK) CALL DR_HOOK('RTTOV_GET_BRDF',1_jpim,ZHOOK_HANDLE)

  CATCH

  IF (LHOOK) CALL DR_HOOK('RTTOV_GET_BRDF',1_jpim,ZHOOK_HANDLE)

END SUBROUTINE rttov_get_brdf
