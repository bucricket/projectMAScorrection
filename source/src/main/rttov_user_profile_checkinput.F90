! Description:
!> @file
!!   Check input profile variables are physically realistic.
!
!> @brief
!!   Check input profile variables are physically realistic and compare
!!   profile values against the regression limits used in coef training.
!!
!! @details
!!   Unphysical values return a fatal error status.
!!   If opts\%rt_config\%verbose is true then warnings are printed when
!!   regression limits are exceeded. If reg_limits_exceeded is passed in
!!   this is set to true if the limits are exceeded.
!!
!!   Comparisons for gas inputs are all carried out in units of ppmv
!!   over dry air as these are the units used in the coefficient
!!   generation. Where limits are exceeded values are reported in
!!   ppmv over dry air.
!!
!!   This subroutine tests a single profile on any pressure levels:
!!   profile levels are tested against regression limits for the nearest
!!   higher pressure coefficient level.
!!
!!   This can be used to check profiles before calling RTTOV instead of
!!   using the internal RTTOV profile checking which is controlled by
!!   opts\%config\%do_checkinput.
!!
!!   The subroutine also optionally checks aer/cld_opt_param structures:
!!   these are checked for unphysical values across *all* channels/profiles.
!!
!! @param[out]    err                  status on exit
!! @param[in]     opts                 options to configure the simulations
!! @param[in]     coefs                coefficients structure for instrument to simulate
!! @param[in]     prof                 input profile to check
!! @param[in]     aer_opt_param        input aerosol optical parameter profile to check, optional
!! @param[in]     cld_opt_param        input cloud optical parameter profile to check, optional
!! @param[in,out] reg_limits_exceeded  output flag indicating whether regression limits were exceeded, optional
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
SUBROUTINE rttov_user_profile_checkinput( &
       & err,                &
       & opts,               &
       & coefs,              &
       & prof,               &
       & aer_opt_param,      &
       & cld_opt_param,      &
       & reg_limits_exceeded )

!INTF_OFF
#include "throw.h"
!INTF_ON

  USE rttov_types, ONLY :  &
         & rttov_coefs,    &
         & rttov_options,  &
         & profile_type,   &
         & rttov_opt_param

  USE parkind1, ONLY : jpim, jplm

!INTF_OFF
  USE rttov_const, ONLY :     &
      ngases_unit,            &
      nsurftype,              &
      surftype_sea,           &
      nwatertype,             &
      wcl_id_stco,            &
      wcl_id_cucc,            &
      wcl_id_cucp,            &
      gas_id_watervapour,     &
      gas_id_ozone,           &
      gas_id_co2,             &
      gas_id_co,              &
      gas_id_ch4,             &
      gas_id_n2o,             &
      tmax, tmin,             &
      qmax, qmin,             &
      o3max, o3min,           &
      co2max, co2min,         &
      comax, comin,           &
      n2omax, n2omin,         &
      ch4max, ch4min,         &
      clwmax, clwmin,         &
      pmax, pmin,             &
      wmax,                   &
      zenmax, zenmaxv9,       &
      ctpmax, ctpmin,         &
      bemax, bemin,           &
      dgmin_hex, dgmax_hex,   &
      dgmin_agg, dgmax_agg,   &
      iwcmin_hex, iwcmax_hex, &
      iwcmin_agg, iwcmax_agg, &
      nish, nidg
  USE yomhook, ONLY : &
      LHOOK,&
      DR_HOOK
  USE parkind1, ONLY : jprb
  USE mod_rttov_baran2013_icldata, ONLY : &
      baran2013_iwc_min,  &
      baran2013_iwc_max,  &
      baran2013_temp_min, &
      baran2013_temp_max
!INTF_ON

  IMPLICIT NONE

  INTEGER(KIND=jpim),              INTENT(OUT) :: err                  ! status on exit
  TYPE(rttov_options),             INTENT(IN)  :: opts                 ! options to configure the simulations
  TYPE(rttov_coefs),               INTENT(IN)  :: coefs                ! coefficients structure for instrument to simulate
  TYPE(profile_type),              INTENT(IN)  :: prof                 ! input profile to check
  TYPE(rttov_opt_param), OPTIONAL, INTENT(IN)  :: aer_opt_param        ! input aerosol optical parameter profile to check
  TYPE(rttov_opt_param), OPTIONAL, INTENT(IN)  :: cld_opt_param        ! input cloud optical parameter profile to check
  LOGICAL(KIND=jplm),    OPTIONAL, INTENT(OUT) :: reg_limits_exceeded  ! output flag indicating whether reg. limits were exceeded

!INTF_END
#include "rttov_alloc_prof.interface"
#include "rttov_convert_gas_dry.interface"
#include "rttov_errorreport.interface"
  TYPE(rttov_options) :: opts_gas
  TYPE(profile_type)  :: prof_dry(1)
  REAL(KIND=jprb)     :: wind, zmax
  REAL(KIND=jprb)     :: iwcmin, iwcmax
  REAL(KIND=jprb)     :: dgmin, dgmax
  INTEGER(KIND=jpim)  :: firstlevel, firstcoeflevel, ilev, jlev
  INTEGER(KIND=jpim)  :: ilay
  INTEGER(KIND=jpim)  :: ig
  LOGICAL(KIND=jplm)  :: OK
  LOGICAL(KIND=jplm)  :: reg_limits_warn
  CHARACTER(256)      :: msg
  LOGICAL(KIND=jplm)  :: ltest(SIZE(prof%p))
  REAL(KIND=jprb)     :: ZHOOK_HANDLE

!- End of header --------------------------------------------------------

TRY

  !-------------
  ! Initialize
  !-------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_USER_PROFILE_CHECKINPUT',0_jpim,ZHOOK_HANDLE)

  reg_limits_warn = .FALSE._jplm

  ! We should check all levels which can contribute to output radiances.
  ! Exactly which user levels contribute depends on the user and coef
  ! levels, the surface pressure, and the interpolation mode used.
  ! In a future version we would like to ensure no input levels below
  ! the level immediately below the surface pressure contribute, but
  ! this is not necessarily the case currently.
  ! In general we should check down at least as far as the *coef*
  ! pressure level which lies on or below the surface pressure.

  ! Find first coef level at or below surface
  DO firstcoeflevel = coefs%coef%nlevels, 2, -1
    IF (coefs%coef%ref_prfl_p(firstcoeflevel-1) < prof%s2m%p) EXIT
  ENDDO

  ! Find first user level at or below firstcoeflevel
  DO firstlevel = prof%nlevels, 2, -1
    IF (prof%p(firstlevel-1) < coefs%coef%ref_prfl_p(firstcoeflevel)) EXIT
  ENDDO

  ! Create profile in units of ppmv dry
  opts_gas                 = opts
  opts_gas%rt_ir%addclouds = .FALSE.
  opts_gas%rt_ir%addaerosl = .FALSE.

  IF (prof%gas_units > ngases_unit) THEN
    err = errorstatus_fatal
    THROWM(err .NE. 0, "invalid gas units")
  ENDIF

  CALL rttov_alloc_prof(err, 1_jpim, prof_dry, prof%nlevels, opts_gas, 1_jpim)
  CALL rttov_convert_gas_dry(opts, (/ prof /), prof_dry)


  !------------------------------
  ! Check for unphysical values
  !------------------------------
  ! zenith angle
  IF (coefs%coef%fmv_model_ver == 9) THEN
    zmax = zenmaxv9
  ELSE
    zmax = zenmax
  ENDIF
  IF (prof%zenangle > zmax .OR. &
      prof%zenangle < 0._jprb) THEN
    err = errorstatus_fatal
    THROWM(err .NE. 0, "invalid zenith angle")
  ENDIF

  ! Solar zenith angle
  IF (opts%rt_ir%addsolar) THEN
    IF (prof%sunzenangle < 0._jprb) THEN
      err = errorstatus_fatal
      THROWM(err .NE. 0, "invalid solar zenith angle")
    ENDIF
  ENDIF

  ! Cloud Top Pressure
  IF (prof%ctp > ctpmax .OR. &
      prof%ctp < ctpmin) THEN
    err = errorstatus_fatal
    THROWM(err .NE. 0, "invalid cloud top pressure")
  ENDIF

  ! Cloud Fraction
  IF (prof%cfraction > 1._jprb .OR. &
      prof%cfraction < 0._jprb) THEN
    err = errorstatus_fatal
    THROWM(err .NE. 0, "invalid cloud fraction")
  ENDIF

  ! Zeeman variables
  IF (coefs%coef%inczeeman) THEN
    ! Magnetic field strength
    IF (prof%be > bemax .OR. &
        prof%be < bemin) THEN
      err = errorstatus_fatal
      THROWM(err .NE. 0, "invalid magnetic field strength")
    ENDIF

    ! Cosine of angle between path and mag. field
    IF (prof%cosbk > 1._jprb  .OR. &
        prof%cosbk < -1._jprb) THEN
      err = errorstatus_fatal
      THROWM(err .NE. 0, "invalid cosbk")
    ENDIF
  ENDIF

  !---------------------
  ! Surface variables
  !---------------------

  ! Pressure
  IF (prof%s2m%p > pmax .OR. &
      prof%s2m%p < pmin) THEN
    err = errorstatus_fatal
    THROWM(err .NE. 0, "invalid surface pressure")
  ENDIF

  ! 2m air temperature
  IF (prof%s2m%t > tmax .OR. &
      prof%s2m%t < tmin) THEN
    err = errorstatus_fatal
    THROWM(err .NE. 0, "invalid 2m air temperature")
  ENDIF

  ! 2m water vapour - only used if opts%rt_all%use_q2m is true
  IF (opts%rt_all%use_q2m) THEN
    IF (prof_dry(1)%s2m%q > qmax .OR. &
        prof_dry(1)%s2m%q < qmin) THEN
      err = errorstatus_fatal
      THROWM(err .NE. 0, "invalid 2m water vapour")
    ENDIF
  ENDIF

  ! surface wind speed
  wind = SQRT(prof%s2m%u * prof%s2m%u + &
              prof%s2m%v * prof%s2m%v)
  IF (wind > wmax .OR. &
      wind < 0._jprb) THEN
    err = errorstatus_fatal
    THROWM(err .NE. 0, "invalid 10m wind speed")
  ENDIF

  IF (opts%rt_ir%addsolar) THEN
    IF (prof%s2m%wfetc <= 0._jprb .AND. &
        prof%skin%surftype == surftype_sea) THEN
      err = errorstatus_fatal
      THROWM(err .NE. 0, "invalid wfetc")
    ENDIF
  ENDIF

  ! surface skin temperature
  IF (prof%skin%t > tmax .OR. &
      prof%skin%t < tmin) THEN
    err = errorstatus_fatal
    THROWM(err .NE. 0, "invalid skin surface temperature")
  ENDIF

  ! surface type
  IF (prof%skin%surftype < 0 .OR. &
      prof%skin%surftype > nsurftype) THEN
    err = errorstatus_fatal
    THROWM(err .NE. 0, "invalid surface type")
  ENDIF

  ! water type
  IF (prof%skin%watertype < 0 .OR. &
      prof%skin%watertype > nwatertype) THEN
    err = errorstatus_fatal
    THROWM(err .NE. 0, "invalid water type")
  ENDIF

  ! Foam fraction
  IF (opts%rt_mw%supply_foam_fraction) THEN
    IF (prof%skin%foam_fraction < 0._jprb .OR. &
        prof%skin%foam_fraction > 1._jprb) THEN
      err = errorstatus_fatal
      THROWM(err .NE. 0, "invalid foam fraction")
    ENDIF
  ENDIF

  ! snow fraction
  IF (prof%snow_frac < 0._jprb .OR. &
      prof%snow_frac > 1._jprb) THEN
    err = errorstatus_fatal
    THROWM(err .NE. 0, "invalid snow fraction")
  ENDIF


  !-------------------------
  ! Atmospheric variables
  !-------------------------

  ! Predictors are calculated on *all* levels so check
  ! the hard limits on every level to avoid errors.

  ! temperature
  ltest = (prof%t(:) > tmax)
  IF (ANY(ltest)) THEN
    err = errorstatus_fatal
    CALL print_info("Input temperature profile exceeds allowed maximum:", &
      (/tmax/), PACK(prof%p(:), mask=ltest), PACK(prof%t(:), mask=ltest))
  ENDIF

  ltest = (prof%t(:) < tmin)
  IF (ANY(ltest)) THEN
    err = errorstatus_fatal
    CALL print_info("Input temperature profile exceeds allowed minimum:", &
      (/tmin/), PACK(prof%p(:), mask=ltest), PACK(prof%t(:), mask=ltest))
  ENDIF
  THROWM(err .NE. 0, "some invalid atmospheric temperature")

  ! water vapour
  ltest = (prof_dry(1)%q(:) > qmax)
  IF (ANY(ltest)) THEN
    err = errorstatus_fatal
    CALL print_info("Input water vapour profile exceeds allowed maximum:", &
      (/qmax/), PACK(prof%p(:), mask=ltest), PACK(prof_dry(1)%q(:), mask=ltest))
  ENDIF

  ltest = (prof_dry(1)%q(:) < qmin)
  IF (ANY(ltest)) THEN
    err = errorstatus_fatal
    CALL print_info("Input water vapour profile exceeds allowed minimum:", &
      (/qmin/), PACK(prof%p(:), mask=ltest), PACK(prof_dry(1)%q(:), mask=ltest))
  ENDIF
  THROWM(err .NE. 0, "some invalid atmospheric water vapour")

  ! ozone
  IF (opts%rt_ir%ozone_Data .AND. coefs%coef%nozone > 0) THEN
    ltest = (prof_dry(1)%o3(:) > o3max)
    IF (ANY(ltest)) THEN
      err = errorstatus_fatal
      CALL print_info("Input ozone profile exceeds allowed maximum:", &
        (/o3max/), PACK(prof%p(:), mask=ltest), PACK(prof_dry(1)%o3(:), mask=ltest))
    ENDIF

    ltest = (prof_dry(1)%o3(:) < o3min)
    IF (ANY(ltest)) THEN
      err = errorstatus_fatal
      CALL print_info("Input ozone profile exceeds allowed minimum:", &
        (/o3min/), PACK(prof%p(:), mask=ltest), PACK(prof_dry(1)%o3(:), mask=ltest))
    ENDIF
    THROWM(err .NE. 0, "some invalid atmospheric ozone")
  ENDIF

  ! CO2
  IF (opts%rt_ir%co2_Data .AND. coefs%coef%nco2 > 0) THEN
    ltest = (prof_dry(1)%co2(:) > co2max)
    IF (ANY(ltest)) THEN
      err = errorstatus_fatal
      CALL print_info("Input CO2 profile exceeds allowed maximum:", &
        (/co2max/), PACK(prof%p(:), mask=ltest), PACK(prof_dry(1)%co2(:), mask=ltest))
    ENDIF

    ltest = (prof_dry(1)%co2(:) < co2min)
    IF (ANY(ltest)) THEN
      err = errorstatus_fatal
      CALL print_info("Input CO2 profile exceeds allowed minimum:", &
        (/co2min/), PACK(prof%p(:), mask=ltest), PACK(prof_dry(1)%co2(:), mask=ltest))
    ENDIF
    THROWM(err .NE. 0, "some invalid atmospheric CO2")
  ENDIF

  ! CO
  IF (opts%rt_ir%co_Data .AND. coefs%coef%nco > 0) THEN
    ltest = (prof_dry(1)%co(:) > comax)
    IF (ANY(ltest)) THEN
      err = errorstatus_fatal
      CALL print_info("Input CO profile exceeds allowed maximum:", &
        (/comax/), PACK(prof%p(:), mask=ltest), PACK(prof_dry(1)%co(:), mask=ltest))
    ENDIF

    ltest = (prof_dry(1)%co(:) < comin)
    IF (ANY(ltest)) THEN
      err = errorstatus_fatal
      CALL print_info("Input CO profile exceeds allowed minimum:", &
        (/comin/), PACK(prof%p(:), mask=ltest), PACK(prof_dry(1)%co(:), mask=ltest))
    ENDIF
    THROWM(err .NE. 0, "some invalid atmospheric CO")
  ENDIF

  ! N2O
  IF (opts%rt_ir%n2o_Data .AND. coefs%coef%nn2o > 0) THEN
    ltest = (prof_dry(1)%n2o(:) > n2omax)
    IF (ANY(ltest)) THEN
      err = errorstatus_fatal
      CALL print_info("Input N2O profile exceeds allowed maximum:", &
        (/n2omax/), PACK(prof%p(:), mask=ltest), PACK(prof_dry(1)%n2o(:), mask=ltest))
    ENDIF

    ltest = (prof_dry(1)%n2o(:) < n2omin)
    IF (ANY(ltest)) THEN
      err = errorstatus_fatal
      CALL print_info("Input N2O profile exceeds allowed minimum:", &
        (/n2omin/), PACK(prof%p(:), mask=ltest), PACK(prof_dry(1)%n2o(:), mask=ltest))
    ENDIF
    THROWM(err .NE. 0, "some invalid atmospheric N2O")
  ENDIF

  ! CH4
  IF (opts%rt_ir%ch4_Data .AND. coefs%coef%nch4 > 0) THEN
    ltest = (prof_dry(1)%ch4(:) > ch4max)
    IF (ANY(ltest)) THEN
      err = errorstatus_fatal
      CALL print_info("Input CH4 profile exceeds allowed maximum:", &
        (/ch4max/), PACK(prof%p(:), mask=ltest), PACK(prof_dry(1)%ch4(:), mask=ltest))
    ENDIF

    ltest = (prof_dry(1)%ch4(:) < ch4min)
    IF (ANY(ltest)) THEN
      err = errorstatus_fatal
      CALL print_info("Input CH4 profile exceeds allowed minimum:", &
        (/ch4min/), PACK(prof%p(:), mask=ltest), PACK(prof_dry(1)%ch4(:), mask=ltest))
    ENDIF
    THROWM(err .NE. 0, "some invalid atmospheric CH4")
  ENDIF

  ! cloud liquid water
  IF (opts%rt_mw%clw_Data) THEN
    ltest = (prof%clw(:) > clwmax)
    IF (ANY(ltest)) THEN
      err = errorstatus_fatal
      CALL print_info("Input cloud liquid water profile exceeds allowed maximum:", &
        (/clwmax/), PACK(prof%p(:), mask=ltest), PACK(prof%clw(:), mask=ltest))
    ENDIF

    ltest = (prof%clw(:) < clwmin)
    IF (ANY(ltest)) THEN
      err = errorstatus_fatal
      CALL print_info("Input cloud liquid water profile exceeds allowed minimum:", &
        (/clwmin/), PACK(prof%p(:), mask=ltest), PACK(prof%clw(:), mask=ltest))
    ENDIF
    THROWM(err .NE. 0, "some invalid atmospheric cloud liquid water")
  ENDIF

  ! Cloud input profile
  IF (opts%rt_ir%addclouds) THEN
    IF (.NOT. ASSOCIATED(prof%cfrac)) THEN
      err = errorstatus_fatal
      msg = "profiles structure not allocated for clouds; "// &
            "opts%rt_ir%addclouds must be true when calling rttov_alloc_prof"
      THROWM(err .NE. 0, msg)
    ENDIF

    IF (ANY( prof%cfrac(:) > 1._jprb ) .OR. &
        ANY( prof%cfrac(:) < 0._jprb )) THEN
      err = errorstatus_fatal
      THROWM(err .NE. 0, "some invalid cloud profile fraction (cfrac) (must be in [0,1])")
    ENDIF

    IF (.NOT. opts%rt_ir%user_cld_opt_param) THEN

      IF (prof%ish < 1_jpim .OR. prof%ish > nish) THEN
        err = errorstatus_fatal
        THROWM(err .NE. 0, "invalid ice particle shape (ish)")
      ENDIF

      IF (prof%ish < 3) THEN ! Not Baran scheme
        IF (prof%idg < 1_jpim .OR. prof%idg > nidg) THEN
          err = errorstatus_fatal
          THROWM(err .NE. 0, "invalid ice particle parameterisation (idg)")
        ENDIF
      ENDIF

      IF (ANY( prof%cloud(:,:) < 0._jprb )) THEN
        err = errorstatus_fatal
        THROWM(err .NE. 0, "some invalid cloud concentration (must be >=0)")
      ENDIF

      IF (ANY( prof%icede(:) < 0._jprb )) THEN
        err = errorstatus_fatal
        THROWM(err .NE. 0, "some invalid ice effective diameter (must be >=0)")
      ENDIF

      IF (opts%rt_ir%pc%addpc) THEN
        IF (ANY( prof%cloud(wcl_id_stco,:) > 0._jprb ) .OR. &
            ANY( prof%cloud(wcl_id_cucc,:) > 0._jprb ) .OR. &
            ANY( prof%cloud(wcl_id_cucp,:) > 0._jprb )) THEN
          err = errorstatus_fatal
          THROWM(err .NE. 0, "only maritime water cloud types are allowed for cloudy PC-RTTOV")
        ENDIF
      ENDIF

    ELSE IF (PRESENT(cld_opt_param)) THEN

      IF (ANY( cld_opt_param%abs(:,:) < 0._jprb )) THEN
        err = errorstatus_fatal
        THROWM(err .NE. 0, "some invalid cloud absorption coefficient (must be >=0)")
      ENDIF

      IF (ANY( cld_opt_param%sca(:,:) < 0._jprb )) THEN
        err = errorstatus_fatal
        THROWM(err .NE. 0, "some invalid cloud scattering coefficient (must be >=0)")
      ENDIF

      IF (ANY( cld_opt_param%bpr(:,:) < 0._jprb ) .OR. &
          ANY( cld_opt_param%bpr(:,:) > 1._jprb )) THEN
        err = errorstatus_fatal
        THROWM(err .NE. 0, "some invalid cloud bpr coefficient (must be in [0,1])")
      ENDIF

    ENDIF

  ENDIF

  ! Aerosol input profile
  IF (opts%rt_ir%addaerosl) THEN

    IF (.NOT. opts%rt_ir%user_aer_opt_param) THEN

      IF (.NOT. ASSOCIATED(prof%aerosols)) THEN
        err = errorstatus_fatal
        msg = "profiles structure not allocated for aerosols; " // &
              "opts%rt_ir%addaerosl must be true when calling rttov_alloc_prof"
        THROWM(err .NE. 0, msg)
      ENDIF

      IF (ANY( prof%aerosols(:,:) < 0._jprb )) THEN
        err = errorstatus_fatal
        THROWM(err .NE. 0, "some invalid aerosol concentration (must be >=0)")
      ENDIF

    ELSE IF (PRESENT(aer_opt_param)) THEN

      IF (ANY( aer_opt_param%abs(:,:) < 0._jprb )) THEN
        err = errorstatus_fatal
        THROWM(err .NE. 0, "some invalid aerosol absorption coefficient (must be >=0)")
      ENDIF

      IF (ANY( aer_opt_param%sca(:,:) < 0._jprb )) THEN
        err = errorstatus_fatal
        THROWM(err .NE. 0, "some invalid aerosol scattering coefficient (must be >=0)")
      ENDIF

      IF (ANY( aer_opt_param%bpr(:,:) < 0._jprb ) .OR. &
          ANY( aer_opt_param%bpr(:,:) > 1._jprb )) THEN
        err = errorstatus_fatal
        THROWM(err .NE. 0, "some invalid aerosol bpr coefficient (must be in [0,1])")
      ENDIF

    ENDIF

  ENDIF


  !---------------------------------
  ! Check against regression limits
  !---------------------------------

  ! Check cloud values for both PC and non-PC simulations

  IF (opts%config%verbose .OR. PRESENT(reg_limits_exceeded)) THEN

    IF (opts%rt_ir%addclouds .AND. .NOT. opts%rt_ir%user_cld_opt_param) THEN
      DO ilay = 1, firstlevel-1

        IF (prof%cloud(6_jpim,ilay) > 0._jprb) THEN

          IF (prof%ish == 1_jpim) THEN
            iwcmin = iwcmin_hex
            iwcmax = iwcmax_hex
          ELSE IF (prof%ish == 2_jpim) THEN
            iwcmin = iwcmin_agg
            iwcmax = iwcmax_agg
          ELSE IF (prof%ish >= 3_jpim) THEN
            iwcmin = baran2013_iwc_min
            iwcmax = baran2013_iwc_max
          ELSE
            CYCLE
          ENDIF

          IF (prof%cloud(6_jpim,ilay) < iwcmin) THEN
            reg_limits_warn = .TRUE.
            IF (opts%config%verbose) THEN
              CALL print_info_lev("Input ice water content profile exceeds lower coef limit", ilay, &
               &  iwcmin, prof%p(ilay), prof%cloud(6_jpim,ilay), .FALSE., .TRUE.)
            ENDIF
          ENDIF

          IF (prof%cloud(6_jpim,ilay) > iwcmax) THEN
            reg_limits_warn = .TRUE.
            IF (opts%config%verbose) THEN
              CALL print_info_lev("Input ice water content profile exceeds upper coef limit", ilay, &
               &  iwcmax, prof%p(ilay), prof%cloud(6_jpim,ilay), .FALSE., .TRUE.)
            ENDIF
          ENDIF

          IF (prof%ish >= 3_jpim) THEN
            IF (prof%t(ilay) < baran2013_temp_min) THEN
              reg_limits_warn = .TRUE.
              IF (opts%config%verbose) THEN
                CALL print_info_lev("Temperature profile exceeds lower coef limit for Baran parameterization", ilay, &
                  &  baran2013_temp_min, prof%p(ilay), prof%t(ilay), .FALSE., .FALSE.)
              ENDIF
            ENDIF
            IF (prof%t(ilay) > baran2013_temp_max) THEN
              reg_limits_warn = .TRUE.
              IF (opts%config%verbose) THEN
                CALL print_info_lev("Temperature profile exceeds upper coef limit for Baran parameterization", ilay, &
                  &  baran2013_temp_max, prof%p(ilay), prof%t(ilay), .FALSE., .FALSE.)
              ENDIF
            ENDIF
          ENDIF

          IF (prof%icede(ilay) > 0._jprb) THEN
            IF (prof%ish == 1_jpim) THEN
              dgmin = dgmin_hex
              dgmax = dgmax_hex
            ELSE IF (prof%ish == 2_jpim) THEN
              dgmin = dgmin_agg
              dgmax = dgmax_agg
            ELSE
              CYCLE
            ENDIF

            IF (prof%icede(ilay) < dgmin) THEN
              reg_limits_warn = .TRUE.
              IF (opts%config%verbose) THEN
                CALL print_info_lev("input effective diameter profile exceeds lower coef limit", ilay, &
                  &  dgmin, prof%p(ilay), prof%icede(ilay), .FALSE., .FALSE.)
              ENDIF
            ENDIF
            IF (prof%icede(ilay) > dgmax) THEN
              reg_limits_warn = .TRUE.
              IF (opts%config%verbose) THEN
                CALL print_info_lev("input effective diameter profile exceeds upper coef limit", ilay, &
                  &  dgmax, prof%p(ilay), prof%icede(ilay), .FALSE., .FALSE.)
              ENDIF
            ENDIF

          ENDIF
        ENDIF

      ENDDO
    ENDIF ! clouds

  ENDIF ! verbose or warn reg limits

  ! Check profiles against non-PC or PC limits

  IF (.NOT. opts%rt_ir%pc%addpc) THEN

    IF (opts%config%verbose .OR. PRESENT(reg_limits_exceeded)) THEN

      jlev = 1
      DO ilev = 1, firstlevel
        ! recherche des niveaux RTTOV qui encadrent le niveau user, les pressions sont
        ! arondies à 0.1Pa
        IF (jlev < coefs%coef%nlevels) THEN
          OK =  NINT(prof%p(ilev)*1000) >= NINT(coefs%coef%lim_prfl_p(jlev)*1000)   .AND. &
             &  NINT(prof%p(ilev)*1000) <  NINT(coefs%coef%lim_prfl_p(jlev+1)*1000)
        ENDIF

        DO WHILE (.NOT. OK .AND. jlev < coefs%coef%nlevels)
          jlev = jlev+1
          IF (jlev < coefs%coef%nlevels) THEN
            OK =  NINT(prof%p(ilev)*1000) >= NINT(coefs%coef%lim_prfl_p(jlev)*1000)   .AND. &
               &  NINT(prof%p(ilev)*1000) <  NINT(coefs%coef%lim_prfl_p(jlev+1)*1000)
          ELSE
            ! on est sur le dernier niveau des fichiers de coefs
            OK =  .TRUE.
          ENDIF
        ENDDO

        IF (prof%t(ilev) > coefs%coef%lim_prfl_tmax(jlev)) THEN
          reg_limits_warn = .TRUE.
          IF (opts%config%verbose) THEN
            CALL print_info_lev("Input temperature profile exceeds upper coef limit", ilev, &
              coefs%coef%lim_prfl_tmax(jlev), prof%p(ilev), prof%t(ilev), .TRUE., .FALSE.)
          ENDIF
        ENDIF

        IF (prof%t(ilev) < coefs%coef%lim_prfl_tmin(jlev)) THEN
          reg_limits_warn = .TRUE.
          IF (opts%config%verbose) THEN
            CALL print_info_lev("Input temperature profile exceeds lower coef limit", ilev, &
              coefs%coef%lim_prfl_tmin(jlev), prof%p(ilev), prof%t(ilev), .TRUE., .FALSE.)
          ENDIF
        ENDIF

        ig = coefs%coef%fmv_gas_pos( gas_id_watervapour )
        IF (prof_dry(1)%q(ilev) > coefs%coef%lim_prfl_gmax(jlev, ig)) THEN
          reg_limits_warn = .TRUE.
          IF (opts%config%verbose) THEN
            CALL print_info_lev("Input water vapour profile exceeds upper coef limit", ilev, &
              coefs%coef%lim_prfl_gmax(jlev, ig), prof%p(ilev), prof_dry(1)%q(ilev), .TRUE., .FALSE.)
          ENDIF
        ENDIF

        IF (prof_dry(1)%q(ilev) < coefs%coef%lim_prfl_gmin(jlev, ig)) THEN
          reg_limits_warn = .TRUE.
          IF (opts%config%verbose) THEN
            CALL print_info_lev("Input water vapour profile exceeds lower coef limit", ilev, &
              coefs%coef%lim_prfl_gmin(jlev, ig), prof%p(ilev), prof_dry(1)%q(ilev), .TRUE., .FALSE.)
          ENDIF
        ENDIF

        IF (opts%rt_ir%ozone_Data .AND. coefs%coef%nozone > 0) THEN
          ig = coefs%coef%fmv_gas_pos( gas_id_ozone )
          IF (prof_dry(1)%o3(ilev) > coefs%coef%lim_prfl_gmax(jlev, ig)) THEN
            reg_limits_warn = .TRUE.
            IF (opts%config%verbose) THEN
              CALL print_info_lev("Input ozone profile exceeds upper coef limit", ilev, &
                coefs%coef%lim_prfl_gmax(jlev, ig), prof%p(ilev), prof_dry(1)%o3(ilev), .TRUE., .FALSE.)
            ENDIF
          ENDIF

          IF (prof_dry(1)%o3(ilev) < coefs%coef%lim_prfl_gmin(jlev, ig)) THEN
            reg_limits_warn = .TRUE.
            IF (opts%config%verbose) THEN
              CALL print_info_lev("Input ozone profile exceeds lower coef limit", ilev, &
                coefs%coef%lim_prfl_gmin(jlev, ig), prof%p(ilev), prof_dry(1)%o3(ilev), .TRUE., .FALSE.)
            ENDIF
          ENDIF
        ENDIF

        IF (opts%rt_ir%co2_Data .AND. coefs%coef%nco2 > 0) THEN
          ig = coefs%coef%fmv_gas_pos( gas_id_co2 )
          IF (prof_dry(1)%co2(ilev) > coefs%coef%lim_prfl_gmax(jlev, ig)) THEN
            reg_limits_warn = .TRUE.
            IF (opts%config%verbose) THEN
              CALL print_info_lev("Input CO2 profile exceeds upper coef limit", ilev, &
                coefs%coef%lim_prfl_gmax(jlev, ig), prof%p(ilev), prof_dry(1)%co2(ilev), .TRUE., .FALSE.)
            ENDIF
          ENDIF

          IF (prof_dry(1)%co2(ilev) < coefs%coef%lim_prfl_gmin(jlev, ig)) THEN
            reg_limits_warn = .TRUE.
            IF (opts%config%verbose) THEN
              CALL print_info_lev("Input CO2 profile exceeds lower coef limit", ilev, &
                coefs%coef%lim_prfl_gmin(jlev, ig), prof%p(ilev), prof_dry(1)%co2(ilev), .TRUE., .FALSE.)
            ENDIF
          ENDIF
        ENDIF

        IF (opts%rt_ir%co_Data .AND. coefs%coef%nco > 0) THEN
          ig = coefs%coef%fmv_gas_pos( gas_id_co )
          IF (prof_dry(1)%co(ilev) > coefs%coef%lim_prfl_gmax(jlev, ig)) THEN
            reg_limits_warn = .TRUE.
            IF (opts%config%verbose) THEN
              CALL print_info_lev("Input CO profile exceeds upper coef limit", ilev, &
                coefs%coef%lim_prfl_gmax(jlev, ig), prof%p(ilev), prof_dry(1)%co(ilev), .TRUE., .FALSE.)
            ENDIF
          ENDIF

          IF (prof_dry(1)%co(ilev) < coefs%coef%lim_prfl_gmin(jlev, ig)) THEN
            reg_limits_warn = .TRUE.
            IF (opts%config%verbose) THEN
              CALL print_info_lev("Input CO profile exceeds lower coef limit", ilev, &
                coefs%coef%lim_prfl_gmin(jlev, ig), prof%p(ilev), prof_dry(1)%co(ilev), .TRUE., .FALSE.)
            ENDIF
          ENDIF
        ENDIF

        IF (opts%rt_ir%n2o_Data .AND. coefs%coef%nn2o > 0) THEN
          ig = coefs%coef%fmv_gas_pos( gas_id_n2o )
          IF (prof_dry(1)%n2o(ilev) > coefs%coef%lim_prfl_gmax(jlev, ig)) THEN
            reg_limits_warn = .TRUE.
            IF (opts%config%verbose) THEN
              CALL print_info_lev("Input N2O profile exceeds upper coef limit", ilev, &
                coefs%coef%lim_prfl_gmax(jlev, ig), prof%p(ilev), prof_dry(1)%n2o(ilev), .TRUE., .FALSE.)
            ENDIF
          ENDIF

          IF (prof_dry(1)%n2o(ilev) < coefs%coef%lim_prfl_gmin(jlev, ig)) THEN
            reg_limits_warn = .TRUE.
            IF (opts%config%verbose) THEN
              CALL print_info_lev("Input N2O profile exceeds lower coef limit", ilev, &
                coefs%coef%lim_prfl_gmin(jlev, ig), prof%p(ilev), prof_dry(1)%n2o(ilev), .TRUE., .FALSE.)
            ENDIF
          ENDIF
        ENDIF

        IF (opts%rt_ir%ch4_Data .AND. coefs%coef%nch4 > 0) THEN
          ig = coefs%coef%fmv_gas_pos( gas_id_ch4 )
          IF (prof_dry(1)%ch4(ilev) > coefs%coef%lim_prfl_gmax(jlev, ig)) THEN
            reg_limits_warn = .TRUE.
            IF (opts%config%verbose) THEN
              CALL print_info_lev("Input CH4 profile exceeds upper coef limit", ilev, &
                coefs%coef%lim_prfl_gmax(jlev, ig), prof%p(ilev), prof_dry(1)%ch4(ilev), .TRUE., .FALSE.)
            ENDIF
          ENDIF

          IF (prof_dry(1)%ch4(ilev) < coefs%coef%lim_prfl_gmin(jlev, ig)) THEN
            reg_limits_warn = .TRUE.
            IF (opts%config%verbose) THEN
              CALL print_info_lev("Input CH4 profile exceeds lower coef limit", ilev, &
                coefs%coef%lim_prfl_gmin(jlev, ig), prof%p(ilev), prof_dry(1)%ch4(ilev), .TRUE., .FALSE.)
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDIF ! verbose or warn reg limits

  ELSE IF (opts%rt_ir%pc%addpc) THEN

    IF (opts%config%verbose .OR. PRESENT(reg_limits_exceeded)) THEN

      IF ((prof%s2m%p < coefs%coef_pccomp%lim_pc_prfl_pmin) .OR. &
        & (prof%s2m%p > coefs%coef_pccomp%lim_pc_prfl_pmax)) THEN
        reg_limits_warn = .TRUE.
        IF (opts%config%verbose) THEN
          WARN("PC-RTTOV: surface pressure outside limits")
        ENDIF
      ENDIF

      IF ((prof%s2m%t < coefs%coef_pccomp%lim_pc_prfl_tsmin) .OR. &
        & (prof%s2m%t > coefs%coef_pccomp%lim_pc_prfl_tsmax)) THEN
        reg_limits_warn = .TRUE.
        IF (opts%config%verbose) THEN
          WARN("PC-RTTOV: surface temperature outside limits")
        ENDIF
      ENDIF

      IF ((prof%skin%t < coefs%coef_pccomp%lim_pc_prfl_skmin) .OR. &
        & (prof%skin%t > coefs%coef_pccomp%lim_pc_prfl_skmax)) THEN
        reg_limits_warn = .TRUE.
        IF (opts%config%verbose) THEN
          WARN("PC-RTTOV: skin temperature outside limits")
        ENDIF
      ENDIF

      wind = SQRT(&
         & prof%s2m%u * prof%s2m%u + &
         & prof%s2m%v * prof%s2m%v   )

      IF ((wind < coefs%coef_pccomp%lim_pc_prfl_wsmin) .OR. &
        & (wind > coefs%coef_pccomp%lim_pc_prfl_wsmax)) THEN
        reg_limits_warn = .TRUE.
        IF (opts%config%verbose) THEN
          WARN("PC-RTTOV: 10m wind speed outside limits")
        ENDIF
      ENDIF


      jlev = 1
      DO ilev = 1, firstlevel
        ! recherche des niveaux RTTOV qui encadrent le niveau user, les pressions sont
        ! arondies à 0.1Pa
        OK =  NINT(prof%p(ilev)*1000) >= NINT(coefs%coef%lim_prfl_p(jlev)*1000)   .AND. &
           &  NINT(prof%p(ilev)*1000) <  NINT(coefs%coef%lim_prfl_p(jlev+1)*1000)

        DO WHILE (.NOT. OK .AND. jlev < coefs%coef%nlevels)
          jlev = jlev+1
          IF (jlev < coefs%coef%nlevels) THEN
            OK =  NINT(prof%p(ilev)*1000) >= NINT(coefs%coef%lim_prfl_p(jlev)*1000)   .AND. &
               &  NINT(prof%p(ilev)*1000) <  NINT(coefs%coef%lim_prfl_p(jlev+1)*1000)
          ELSE
            ! on est sur le dernier niveau des fichiers de coefs
            OK = .TRUE.
          ENDIF
        ENDDO

        IF (prof%t(ilev) > coefs%coef_pccomp%lim_pc_prfl_tmax(jlev)) THEN
          reg_limits_warn = .TRUE.
          IF (opts%config%verbose) THEN
            CALL print_info_lev("PC-RTTOV: Input temperature profile exceeds upper coef limit", ilev, &
              coefs%coef_pccomp%lim_pc_prfl_tmax(jlev), prof%p(ilev), prof%t(ilev), .TRUE., .FALSE.)
          ENDIF
        ENDIF

        IF (prof%t(ilev) < coefs%coef_pccomp%lim_pc_prfl_tmin(jlev)) THEN
          reg_limits_warn = .TRUE.
          IF (opts%config%verbose) THEN
            CALL print_info_lev("PC-RTTOV: Input temperature profile exceeds lower coef limit", ilev, &
              coefs%coef_pccomp%lim_pc_prfl_tmin(jlev), prof%p(ilev), prof%t(ilev), .TRUE., .FALSE.)
          ENDIF
        ENDIF

        ig = coefs%coef%fmv_gas_pos( gas_id_watervapour )
        IF (prof_dry(1)%q(ilev) > coefs%coef_pccomp%lim_pc_prfl_qmax(jlev)) THEN
          reg_limits_warn = .TRUE.
          IF (opts%config%verbose) THEN
            CALL print_info_lev("PC-RTTOV: Input water vapour profile exceeds upper coef limit", ilev, &
              coefs%coef_pccomp%lim_pc_prfl_qmax(jlev), prof%p(ilev), prof_dry(1)%q(ilev), .TRUE., .FALSE.)
          ENDIF
        ENDIF

        IF (prof_dry(1)%q(ilev) < coefs%coef_pccomp%lim_pc_prfl_qmin(jlev)) THEN
          reg_limits_warn = .TRUE.
          IF (opts%config%verbose) THEN
            CALL print_info_lev("PC-RTTOV: Input water vapour profile exceeds lower coef limit", ilev, &
              coefs%coef_pccomp%lim_pc_prfl_qmin(jlev), prof%p(ilev), prof_dry(1)%q(ilev), .TRUE., .FALSE.)
          ENDIF
        ENDIF

        IF (opts%rt_ir%ozone_Data .AND. coefs%coef%nozone > 0) THEN
          ig = coefs%coef%fmv_gas_pos( gas_id_ozone )
          IF (prof_dry(1)%o3(ilev) > coefs%coef_pccomp%lim_pc_prfl_ozmax(jlev)) THEN
            reg_limits_warn = .TRUE.
            IF (opts%config%verbose) THEN
              CALL print_info_lev("PC-RTTOV: Input ozone profile exceeds upper coef limit", ilev, &
                coefs%coef_pccomp%lim_pc_prfl_ozmax(jlev), prof%p(ilev), prof_dry(1)%o3(ilev), .TRUE., .FALSE.)
            ENDIF
          ENDIF

          IF (prof_dry(1)%o3(ilev) < coefs%coef_pccomp%lim_pc_prfl_ozmin(jlev)) THEN
            reg_limits_warn = .TRUE.
            IF (opts%config%verbose) THEN
              CALL print_info_lev("PC-RTTOV: Input ozone profile exceeds lower coef limit", ilev, &
                coefs%coef_pccomp%lim_pc_prfl_ozmin(jlev), prof%p(ilev), prof_dry(1)%o3(ilev), .TRUE., .FALSE.)
            ENDIF
          ENDIF
        ENDIF

      ENDDO

    ENDIF ! verbose or warn reg limits

  ENDIF

  CALL rttov_alloc_prof(err, 1_jpim, prof_dry, prof%nlevels, opts_gas, 0_jpim)

  IF (PRESENT(reg_limits_exceeded)) reg_limits_exceeded = reg_limits_warn

IF (LHOOK) CALL DR_HOOK('RTTOV_USER_PROFILE_CHECKINPUT',1_jpim,ZHOOK_HANDLE)

CATCH

IF (LHOOK) CALL DR_HOOK('RTTOV_USER_PROFILE_CHECKINPUT',1_jpim,ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE print_info(msg1, limits, levels, values)
    CHARACTER(LEN=*), INTENT(IN) :: msg1
    REAL(KIND=jprb),  INTENT(IN) :: limits(:)
    REAL(KIND=jprb),  INTENT(IN) :: levels(:)
    REAL(KIND=jprb),  INTENT(IN) :: values(:)

    CHARACTER(LEN=256) :: msg2
    INTEGER(KIND=jpim) :: imax, lmax

    ! Print a maximum of 10 elements of the input arrays. Note limits(:) may be a
    ! single element array for hard limits so check size of this separately
    imax = MIN(10, SIZE(levels))
    lmax = MIN(imax, SIZE(limits))
    WARN(msg1)
    WRITE(msg2, '(a,10f10.4)') 'Limit   = ',limits(1:lmax)
    INFO(TRIM(msg2))
    WRITE(msg2, '(a,10f10.4)') 'p (hPa) = ',levels(1:imax)
    INFO(TRIM(msg2))
    WRITE(msg2, '(a,10f10.4)') 'Value   = ',values(1:imax)
    INFO(TRIM(msg2))
  END SUBROUTINE print_info

  SUBROUTINE print_info_lev(msg1, ilevel, limit, p, value, llevel, small)
    CHARACTER(LEN=*),   INTENT(IN) :: msg1    ! General warning message to print out
    INTEGER(KIND=jpim), INTENT(IN) :: ilevel  ! Level or layer index
    REAL(KIND=jprb),    INTENT(IN) :: limit   ! Limit value
    REAL(KIND=jprb),    INTENT(IN) :: p       ! Pressure of level
    REAL(KIND=jprb),    INTENT(IN) :: value   ! Input profile value
    LOGICAL(KIND=jplm), INTENT(IN) :: llevel  ! True/false for level/layer quantities resp.
    LOGICAL(KIND=jplm), INTENT(IN) :: small   ! True if profile values are small numbers

    CHARACTER(LEN=32)  :: s
    CHARACTER(LEN=256) :: msg2

    ! Append level/layer number to input message
    IF (llevel) THEN
      WRITE(s, '(" (level number = ",I4,")")') ilevel
    ELSE
      WRITE(s, '(" (layer number = ",I4,")")') ilevel
    ENDIF
    WARN(msg1//TRIM(s))

    ! Format the limit/pressure/value according to input logical small
    IF (small) THEN
      WRITE(msg2, '((a,e12.5,", "),(a,f10.4,", "),(a,e12.5,", "))') &
        'Limit   = ', limit, 'p (hPa) = ', p, 'Value   = ',value
    ELSE
      WRITE(msg2, '(3(a,f10.4,", "))') &
        'Limit   = ', limit, 'p (hPa) = ', p, 'Value   = ',value
    ENDIF
    INFO(TRIM(msg2))
  END SUBROUTINE print_info_lev

END SUBROUTINE rttov_user_profile_checkinput
