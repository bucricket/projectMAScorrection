! Description:
!> @file
!!   Convert input gas profiles to ppmv wrt dry air AD.
!
!> @brief
!!   Convert input gas profiles to ppmv wrt dry air AD.
!!
!! @details
!!   Convert input gas profiles to ppmv wrt dry air AD.
!!
!! @param[in]     opts             options to configure the simulations
!! @param[in]     profiles         profiles structure containing input gas profiles
!! @param[in,out] profiles_ad      profiles structure containing profile gas perturbations
!! @param[in]     profiles_dry     profiles structure containing converted profiles on exit
!! @param[in]     profiles_dry_ad  profiles structure containing converted profile perturbations
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
SUBROUTINE rttov_convert_gas_dry_ad(opts, profiles, profiles_ad, profiles_dry, profiles_dry_ad)

  USE rttov_types, ONLY : &
      rttov_options,      &
      profile_type
!INTF_OFF
  USE parkind1, ONLY : jpim, jprb

  USE rttov_const, ONLY : &
      gas_unit_specconc,  &
      gas_unit_ppmv,      &
      gas_id_watervapour, &
      gas_id_ozone,       &
      gas_id_co2,         &
      gas_id_n2o,         &
      gas_id_co,          &
      gas_id_ch4,         &
      mair,               &
      gas_mass

  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options), INTENT(IN)    :: opts
  TYPE(profile_type),  INTENT(IN)    :: profiles(:)
  TYPE(profile_type),  INTENT(INOUT) :: profiles_ad(SIZE(profiles))
  TYPE(profile_type),  INTENT(IN)    :: profiles_dry(SIZE(profiles))
  TYPE(profile_type),  INTENT(IN)    :: profiles_dry_ad(SIZE(profiles))
!INTF_END

  INTEGER(jpim) :: iprof, nprofiles
  REAL(jprb)    :: factor(profiles(1)%nlevels), factor2m
  REAL(jprb)    :: factor_ad(profiles(1)%nlevels), factor2m_ad
  REAL(jprb)    :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CONVERT_GAS_DRY_AD', 0_jpim, ZHOOK_HANDLE)

  nprofiles = SIZE(profiles)

  SELECT CASE (profiles(1)%gas_units)

  CASE (gas_unit_specconc)

!    vdry_gas = qwet_gas / (1._jprb - qwet_gas) * 1.E+06_jprb * Mair / gas_mass(gas_id)
!    vdry_gas_tl = qwet_gas_tl / (1._jprb - qwet_gas)**2 * 1.E+06_jprb * Mair / gas_mass(gas_id)
!    qwet_gas_ad = qwet_gas_ad + vdry_gas_ad / (1._jprb - qwet_gas)**2 * 1.E+06_jprb * Mair / gas_mass(gas_id)

    DO iprof = 1, nprofiles
      profiles_ad(iprof)%q = profiles_ad(iprof)%q + &
                             1.E+06_jprb * Mair / gas_mass(gas_id_watervapour) * &
                             profiles_dry_ad(iprof)%q / (1._jprb - profiles(iprof)%q)**2

      IF (opts%rt_all%use_q2m) profiles_ad(iprof)%s2m%q = profiles_ad(iprof)%s2m%q + &
                                 1.E+06_jprb * Mair / gas_mass(gas_id_watervapour) * &
                                 profiles_dry_ad(iprof)%s2m%q / (1._jprb - profiles(iprof)%s2m%q)**2

      IF (ASSOCIATED(profiles_dry_ad(iprof)%o3) .AND. ASSOCIATED(profiles_ad(iprof)%o3)) THEN
        profiles_ad(iprof)%o3 = profiles_ad(iprof)%o3 + &
                                1.E+06_jprb * Mair / gas_mass(gas_id_ozone) * &
                                profiles_dry_ad(iprof)%o3 / (1._jprb - profiles(iprof)%o3)**2
        IF (opts%rt_all%use_q2m) profiles_ad(iprof)%s2m%o = profiles_ad(iprof)%s2m%o + &
                                   1.E+06_jprb * Mair / gas_mass(gas_id_ozone) * &
                                   profiles_dry_ad(iprof)%s2m%o / (1._jprb - profiles(iprof)%s2m%o)**2
      ENDIF

      IF (ASSOCIATED(profiles_dry_ad(iprof)%co2) .AND. ASSOCIATED(profiles_ad(iprof)%co2)) THEN
        profiles_ad(iprof)%co2 = profiles_ad(iprof)%co2 + &
                                 1.E+06_jprb * Mair / gas_mass(gas_id_co2) * &
                                 profiles_dry_ad(iprof)%co2 / (1._jprb - profiles(iprof)%co2)**2
      ENDIF

      IF (ASSOCIATED(profiles_dry_ad(iprof)%co)  .AND. ASSOCIATED(profiles_ad(iprof)%co))  THEN
        profiles_ad(iprof)%co = profiles_ad(iprof)%co + &
                                1.E+06_jprb * Mair / gas_mass(gas_id_co) * &
                                profiles_dry_ad(iprof)%co / (1._jprb - profiles(iprof)%co)**2
      ENDIF

      IF (ASSOCIATED(profiles_dry_ad(iprof)%ch4) .AND. ASSOCIATED(profiles_ad(iprof)%ch4)) THEN
        profiles_ad(iprof)%ch4 = profiles_ad(iprof)%ch4 + &
                                 1.E+06_jprb * Mair / gas_mass(gas_id_ch4) * &
                                 profiles_dry_ad(iprof)%ch4 / (1._jprb - profiles(iprof)%ch4)**2
      ENDIF

      IF (ASSOCIATED(profiles_dry_ad(iprof)%n2o) .AND. ASSOCIATED(profiles_ad(iprof)%n2o)) THEN
        profiles_ad(iprof)%n2o = profiles_ad(iprof)%n2o + &
                                 1.E+06_jprb * Mair / gas_mass(gas_id_n2o) * &
                                 profiles_dry_ad(iprof)%n2o / (1._jprb - profiles(iprof)%n2o)**2
      ENDIF
    ENDDO


  CASE (gas_unit_ppmv)

!    vdry_gas = vwet_gas / (1._jprb - h2o * 1.E-06_jprb)
!    factor = 1. / (1._jprb - h2o * 1.E-06_jprb)
!    factor_tl = 1.E-06_jprb * h2o_tl * factor**2
!    vdry_gas_tl = vwet_gas_tl * factor + vwet_gas * factor_tl

!    vwet_gas_ad = vwet_gas_ad + factor * vdry_gas_ad
!    factor_ad = factor_ad + vwet_gas * vdry_gas_ad
!    h2o_ad = h2o_ad + 1.E-06_jprb * factor_ad * factor**2

    DO iprof = 1, nprofiles

      factor = 1._jprb / (1._jprb - profiles(iprof)%q * 1.E-06_jprb)
      IF (opts%rt_all%use_q2m) THEN
        factor2m = 1._jprb / (1._jprb - profiles(iprof)%s2m%q * 1.E-06_jprb)
      ENDIF

      profiles_ad(iprof)%q = profiles_ad(iprof)%q + factor * profiles_dry_ad(iprof)%q
      ! NB First use of factor_ad
      factor_ad = profiles(iprof)%q * profiles_dry_ad(iprof)%q

      IF (opts%rt_all%use_q2m) THEN
        profiles_ad(iprof)%s2m%q = profiles_ad(iprof)%s2m%q + factor2m * profiles_dry_ad(iprof)%s2m%q
        ! NB First use of factor2m_ad
        factor2m_ad = profiles(iprof)%s2m%q * profiles_dry_ad(iprof)%s2m%q
      ENDIF

      IF (ASSOCIATED(profiles_dry_ad(iprof)%o3) .AND. ASSOCIATED(profiles_ad(iprof)%o3)) THEN
        profiles_ad(iprof)%o3 = profiles_ad(iprof)%o3 + factor * profiles_dry_ad(iprof)%o3
        factor_ad = factor_ad + profiles(iprof)%o3 * profiles_dry_ad(iprof)%o3

        IF (opts%rt_all%use_q2m) THEN
          profiles_ad(iprof)%s2m%o = profiles_ad(iprof)%s2m%o + factor2m * profiles_dry_ad(iprof)%s2m%o
          factor2m_ad = factor2m_ad + profiles(iprof)%s2m%o * profiles_dry_ad(iprof)%s2m%o
        ENDIF
      ENDIF

      IF (ASSOCIATED(profiles_dry_ad(iprof)%co2) .AND. ASSOCIATED(profiles_ad(iprof)%co2)) THEN
        profiles_ad(iprof)%co2 = profiles_ad(iprof)%co2 + factor * profiles_dry_ad(iprof)%co2
        factor_ad = factor_ad + profiles(iprof)%co2 * profiles_dry_ad(iprof)%co2
      ENDIF

      IF (ASSOCIATED(profiles_dry_ad(iprof)%co)  .AND. ASSOCIATED(profiles_ad(iprof)%co))  THEN
        profiles_ad(iprof)%co = profiles_ad(iprof)%co + factor * profiles_dry_ad(iprof)%co
        factor_ad = factor_ad + profiles(iprof)%co * profiles_dry_ad(iprof)%co
      ENDIF

      IF (ASSOCIATED(profiles_dry_ad(iprof)%ch4) .AND. ASSOCIATED(profiles_ad(iprof)%ch4)) THEN
        profiles_ad(iprof)%ch4 = profiles_ad(iprof)%ch4 + factor * profiles_dry_ad(iprof)%ch4
        factor_ad = factor_ad + profiles(iprof)%ch4 * profiles_dry_ad(iprof)%ch4
      ENDIF

      IF (ASSOCIATED(profiles_dry_ad(iprof)%n2o) .AND. ASSOCIATED(profiles_ad(iprof)%n2o)) THEN
        profiles_ad(iprof)%n2o = profiles_ad(iprof)%n2o + factor * profiles_dry_ad(iprof)%n2o
        factor_ad = factor_ad + profiles(iprof)%n2o * profiles_dry_ad(iprof)%n2o
      ENDIF

      profiles_ad(iprof)%q = profiles_ad(iprof)%q + 1.E-06_jprb * factor_ad * factor**2
      IF (opts%rt_all%use_q2m) THEN
        profiles_ad(iprof)%s2m%q = profiles_ad(iprof)%s2m%q + 1.E-06_jprb * factor2m_ad * factor2m**2
      ENDIF
    ENDDO


  CASE DEFAULT ! Assume no conversion (i.e. compatibility mode or ppmv dry)

    DO iprof = 1, nprofiles
      profiles_ad(iprof)%q = profiles_ad(iprof)%q + profiles_dry_ad(iprof)%q
      profiles_ad(iprof)%s2m%q = profiles_ad(iprof)%s2m%q + profiles_dry_ad(iprof)%s2m%q
      IF (ASSOCIATED(profiles_dry_ad(iprof)%o3) .AND. ASSOCIATED(profiles_ad(iprof)%o3)) THEN
        profiles_ad(iprof)%o3 = profiles_ad(iprof)%o3 + profiles_dry_ad(iprof)%o3
        profiles_ad(iprof)%s2m%o = profiles_ad(iprof)%s2m%o + profiles_dry_ad(iprof)%s2m%o
      ENDIF
      IF (ASSOCIATED(profiles_dry_ad(iprof)%co2) .AND. ASSOCIATED(profiles_ad(iprof)%co2)) &
          profiles_ad(iprof)%co2 = profiles_ad(iprof)%co2 + profiles_dry_ad(iprof)%co2
      IF (ASSOCIATED(profiles_dry_ad(iprof)%co)  .AND. ASSOCIATED(profiles_ad(iprof)%co))  &
          profiles_ad(iprof)%co  = profiles_ad(iprof)%co  + profiles_dry_ad(iprof)%co
      IF (ASSOCIATED(profiles_dry_ad(iprof)%ch4) .AND. ASSOCIATED(profiles_ad(iprof)%ch4)) &
          profiles_ad(iprof)%ch4 = profiles_ad(iprof)%ch4 + profiles_dry_ad(iprof)%ch4
      IF (ASSOCIATED(profiles_dry_ad(iprof)%n2o) .AND. ASSOCIATED(profiles_ad(iprof)%n2o)) &
          profiles_ad(iprof)%n2o = profiles_ad(iprof)%n2o + profiles_dry_ad(iprof)%n2o
    ENDDO

  END SELECT

  IF (LHOOK) CALL DR_HOOK('RTTOV_CONVERT_GAS_DRY_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_convert_gas_dry_ad
