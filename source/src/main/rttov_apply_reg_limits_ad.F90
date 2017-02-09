!
SUBROUTINE rttov_apply_reg_limits_ad( &
       & opts,        &
       & prof_user,   &
       & prof,        &
       & prof_ad,     &
       & coef,        &
       & coef_pccomp)
  ! Description:
  !   Clip profiles to coef file regression limits (AD).
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
  !    Copyright 2014, EUMETSAT, All Rights Reserved.
  !
  ! Current Code Owner: SAF NWP
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".

  USE rttov_types, ONLY :     &
         & rttov_coef,        &
         & rttov_options,     &
         & rttov_coef_pccomp, &
         & profile_type

!INTF_OFF
  USE parkind1, ONLY : jprb, jpim

  USE rttov_const, ONLY : &
      gas_id_watervapour, &
      gas_id_ozone,       &
      gas_id_co2,         &
      gas_id_co,          &
      gas_id_n2o,         &
      gas_id_ch4

  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON

  IMPLICIT NONE

  TYPE(rttov_options),     INTENT(IN)    :: opts
  TYPE(profile_type),      INTENT(IN)    :: prof_user(:) ! Profiles on user levels (only p(:) and 2m p used)
  TYPE(profile_type),      INTENT(IN)    :: prof(:)      ! Profiles on coef levels (ppmv dry)
  TYPE(profile_type),      INTENT(INOUT) :: prof_ad(:)   ! Profiles_ad on coef levels (ppmv dry)
  TYPE(rttov_coef),        INTENT(IN)    :: coef
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
!INTF_END

  REAL(KIND=jprb)    :: wind
  INTEGER(KIND=jpim) :: firstlevel, firstuserlevel, ilev
  INTEGER(KIND=jpim) :: ig
  INTEGER(KIND=jpim) :: nprofiles, iprof
  REAL(KIND=jprb)    :: ZHOOK_HANDLE

!- End of header --------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('RTTOV_APPLY_REG_LIMITS_AD',0_jpim,ZHOOK_HANDLE)

nprofiles = SIZE(prof)

DO iprof = 1, nprofiles

  ! If apply_reg_limits is true then check all profile values from TOA to surface.
  ! Otherwise, if reg_limits_extrap is true then check only the values above the
  !   top of the input profile.

  IF (opts%config%apply_reg_limits) THEN

    ! Find first user level at or below surface
    DO firstuserlevel = prof_user(iprof)%nlevels, 2, -1
      IF (prof_user(iprof)%p(firstuserlevel-1) < prof_user(iprof)%s2m%p) EXIT
    ENDDO

    ! Find first coef level at or below firstuserlevel
    DO firstlevel = prof(iprof)%nlevels, 2, -1
      IF (prof(iprof)%p(firstlevel-1) < prof_user(iprof)%p(firstuserlevel)) EXIT
    ENDDO

  ELSE IF (opts%interpolation%reg_limit_extrap) THEN

    IF (prof(iprof)%p(2) < prof_user(iprof)%p(1)) THEN

      ! Determine first coef level above user input profile top level
      DO firstlevel = 2, prof(iprof)%nlevels - 1
        IF (prof(iprof)%p(firstlevel+1) > prof_user(iprof)%p(1)) EXIT
      ENDDO

    ELSE
      CYCLE
    ENDIF

  ELSE

    ! Don't do anything if apply_reg_limits and reg_limit_extrap are both false
    EXIT

  ENDIF

  IF (.NOT. opts%rt_ir%pc%addpc) THEN

    DO ilev = 1, firstlevel
      IF (prof(iprof)%t(ilev) > coef%lim_prfl_tmax(ilev)) THEN
        prof_ad(iprof)%t(ilev) = 0._jprb
      ELSE IF (prof(iprof)%t(ilev) < coef%lim_prfl_tmin(ilev)) THEN
        prof_ad(iprof)%t(ilev) = 0._jprb
      ENDIF
    ENDDO

    ig = coef%fmv_gas_pos(gas_id_watervapour)
    DO ilev = 1, firstlevel
      IF (prof(iprof)%q(ilev) > coef%lim_prfl_gmax(ilev, ig)) THEN
        prof_ad(iprof)%q(ilev) = 0._jprb
      ELSE IF (prof(iprof)%q(ilev) < coef%lim_prfl_gmin(ilev, ig)) THEN
        prof_ad(iprof)%q(ilev) = 0._jprb
      ENDIF
    ENDDO

    IF (opts%rt_ir%ozone_data .AND. coef%nozone > 0) THEN
      ig = coef%fmv_gas_pos(gas_id_ozone)
      DO ilev = 1, firstlevel
        IF (prof(iprof)%o3(ilev) > coef%lim_prfl_gmax(ilev, ig)) THEN
          prof_ad(iprof)%o3(ilev) = 0._jprb
        ELSE IF (prof(iprof)%o3(ilev) < coef%lim_prfl_gmin(ilev, ig)) THEN
          prof_ad(iprof)%o3(ilev) = 0._jprb
        ENDIF
      ENDDO
    ENDIF

    IF (opts%rt_ir%co2_data .AND. coef%nco2 > 0) THEN
      ig = coef%fmv_gas_pos(gas_id_co2)
      DO ilev = 1, firstlevel
        IF (prof(iprof)%co2(ilev) > coef%lim_prfl_gmax(ilev, ig)) THEN
          prof_ad(iprof)%co2(ilev) = 0._jprb
        ELSE IF (prof(iprof)%co2(ilev) < coef%lim_prfl_gmin(ilev, ig)) THEN
          prof_ad(iprof)%co2(ilev) = 0._jprb
        ENDIF
      ENDDO
    ENDIF

    IF (opts%rt_ir%co_data .AND. coef%nco > 0) THEN
      ig = coef%fmv_gas_pos(gas_id_co)
      DO ilev = 1, firstlevel
        IF (prof(iprof)%co(ilev) > coef%lim_prfl_gmax(ilev, ig)) THEN
          prof_ad(iprof)%co(ilev) = 0._jprb
        ELSE IF (prof(iprof)%co(ilev) < coef%lim_prfl_gmin(ilev, ig)) THEN
          prof_ad(iprof)%co(ilev) = 0._jprb
        ENDIF
      ENDDO
    ENDIF

    IF (opts%rt_ir%n2o_data .AND. coef%nn2o > 0) THEN
      ig = coef%fmv_gas_pos(gas_id_n2o)
      DO ilev = 1, firstlevel
        IF (prof(iprof)%n2o(ilev) > coef%lim_prfl_gmax(ilev, ig)) THEN
          prof_ad(iprof)%n2o(ilev) = 0._jprb
        ELSE IF (prof(iprof)%n2o(ilev) < coef%lim_prfl_gmin(ilev, ig)) THEN
          prof_ad(iprof)%n2o(ilev) = 0._jprb
        ENDIF
      ENDDO
    ENDIF

    IF (opts%rt_ir%ch4_data .AND. coef%nch4 > 0) THEN
      ig = coef%fmv_gas_pos(gas_id_ch4)
      DO ilev = 1, firstlevel
        IF (prof(iprof)%ch4(ilev) > coef%lim_prfl_gmax(ilev, ig)) THEN
          prof_ad(iprof)%ch4(ilev) = 0._jprb
        ELSE IF (prof(iprof)%ch4(ilev) < coef%lim_prfl_gmin(ilev, ig)) THEN
          prof_ad(iprof)%ch4(ilev) = 0._jprb
        ENDIF
      ENDDO
    ENDIF

  ELSE ! addpc

    DO ilev = 1, firstlevel
      IF (prof(iprof)%t(ilev) > coef_pccomp%lim_pc_prfl_tmax(ilev)) THEN
        prof_ad(iprof)%t(ilev) = 0._jprb
      ELSE IF (prof(iprof)%t(ilev) < coef_pccomp%lim_pc_prfl_tmin(ilev)) THEN
        prof_ad(iprof)%t(ilev) = 0._jprb
      ENDIF
    ENDDO

    DO ilev = 1, firstlevel
      IF (prof(iprof)%q(ilev) > coef_pccomp%lim_pc_prfl_qmax(ilev)) THEN
        prof_ad(iprof)%q(ilev) = 0._jprb
      ELSE IF (prof(iprof)%q(ilev) < coef_pccomp%lim_pc_prfl_qmin(ilev)) THEN
        prof_ad(iprof)%q(ilev) = 0._jprb
      ENDIF
    ENDDO

    IF (opts%rt_ir%ozone_data .AND. coef%nozone > 0) THEN
      DO ilev = 1, firstlevel
        IF (prof(iprof)%o3(ilev) > coef_pccomp%lim_pc_prfl_ozmax(ilev)) THEN
          prof_ad(iprof)%o3(ilev) = 0._jprb
        ELSE IF (prof(iprof)%o3(ilev) < coef_pccomp%lim_pc_prfl_ozmin(ilev)) THEN
          prof_ad(iprof)%o3(ilev) = 0._jprb
        ENDIF
      ENDDO
    ENDIF

    IF (prof(iprof)%s2m%p < coef_pccomp%lim_pc_prfl_pmin) THEN
      prof_ad(iprof)%s2m%p = 0._jprb
    ELSE IF (prof(iprof)%s2m%p > coef_pccomp%lim_pc_prfl_pmax) THEN
      prof_ad(iprof)%s2m%p = 0._jprb
    ENDIF

    IF (prof(iprof)%s2m%t < coef_pccomp%lim_pc_prfl_tsmin) THEN
      prof_ad(iprof)%s2m%t = 0._jprb
    ELSE IF (prof(iprof)%s2m%t > coef_pccomp%lim_pc_prfl_tsmax) THEN
      prof_ad(iprof)%s2m%t = 0._jprb
    ENDIF

    IF (prof(iprof)%skin%t < coef_pccomp%lim_pc_prfl_skmin) THEN
      prof_ad(iprof)%skin%t = 0._jprb
    ELSE IF (prof(iprof)%skin%t > coef_pccomp%lim_pc_prfl_skmax) THEN
      prof_ad(iprof)%skin%t = 0._jprb
    ENDIF

    wind = SQRT(prof(iprof)%s2m%u * prof(iprof)%s2m%u + &
                prof(iprof)%s2m%v * prof(iprof)%s2m%v)

    IF (wind < coef_pccomp%lim_pc_prfl_wsmin) THEN
      prof_ad(iprof)%s2m%u = 0._jprb
      prof_ad(iprof)%s2m%v = 0._jprb
    ELSE IF (wind > coef_pccomp%lim_pc_prfl_wsmax) THEN
      prof_ad(iprof)%s2m%u = 0._jprb
      prof_ad(iprof)%s2m%v = 0._jprb
    ENDIF

  ENDIF ! addpc

ENDDO ! profiles

IF (LHOOK) CALL DR_HOOK('RTTOV_APPLY_REG_LIMITS_AD',1_jpim,ZHOOK_HANDLE)

END SUBROUTINE rttov_apply_reg_limits_ad
