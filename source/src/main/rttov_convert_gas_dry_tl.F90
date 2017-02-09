! Description:
!> @file
!!   Convert input gas profiles to ppmv wrt dry air TL.
!
!> @brief
!!   Convert input gas profiles to ppmv wrt dry air TL.
!!
!! @details
!!   Convert input gas profiles to ppmv wrt dry air TL.
!!
!! @param[in]     opts             options to configure the simulations
!! @param[in]     profiles         profiles structure containing input gas profiles
!! @param[in]     profiles_tl      profiles structure containing input gas increments
!! @param[in]     profiles_dry     profiles structure containing converted profiles on exit
!! @param[in,out] profiles_dry_tl  profiles structure containing converted profile increments on exit
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
SUBROUTINE rttov_convert_gas_dry_tl(opts, profiles, profiles_tl, profiles_dry, profiles_dry_tl)

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
  TYPE(profile_type),  INTENT(IN)    :: profiles_tl(SIZE(profiles))
  TYPE(profile_type),  INTENT(IN)    :: profiles_dry(SIZE(profiles))
  TYPE(profile_type),  INTENT(INOUT) :: profiles_dry_tl(SIZE(profiles))
!INTF_END

  INTEGER(jpim) :: iprof, nprofiles
  REAL(jprb)    :: factor(profiles(1)%nlevels), factor2m
  REAL(jprb)    :: factor_tl(profiles(1)%nlevels), factor2m_tl
  REAL(jprb)    :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CONVERT_GAS_DRY_TL', 0_jpim, ZHOOK_HANDLE)

  nprofiles = SIZE(profiles)

  SELECT CASE (profiles(1)%gas_units)

  CASE (gas_unit_specconc)

!    vdry_gas = qwet_gas / (1._jprb - qwet_gas) * 1.E+06_jprb * Mair / gas_mass(gas_id)
!    vdry_gas_tl = qwet_gas_tl / (1._jprb - qwet_gas)**2 * 1.E+06_jprb * Mair / gas_mass(gas_id)

    DO iprof = 1, nprofiles
      profiles_dry_tl(iprof)%q = 1.E+06_jprb * Mair / gas_mass(gas_id_watervapour) * &
                                 profiles_tl(iprof)%q / (1._jprb - profiles(iprof)%q)**2

      IF (opts%rt_all%use_q2m) profiles_dry_tl(iprof)%s2m%q = &
                                 1.E+06_jprb * Mair / gas_mass(gas_id_watervapour) * &
                                 profiles_tl(iprof)%s2m%q / (1._jprb - profiles(iprof)%s2m%q)**2

      IF (ASSOCIATED(profiles_dry_tl(iprof)%o3) .AND. ASSOCIATED(profiles_tl(iprof)%o3)) THEN
        profiles_dry_tl(iprof)%o3 = 1.E+06_jprb * Mair / gas_mass(gas_id_ozone) * &
                                    profiles_tl(iprof)%o3 / (1._jprb - profiles(iprof)%o3)**2
        IF (opts%rt_all%use_q2m) profiles_dry_tl(iprof)%s2m%o = &
                                    1.E+06_jprb * Mair / gas_mass(gas_id_ozone) * &
                                    profiles_tl(iprof)%s2m%o / (1._jprb - profiles(iprof)%s2m%o)**2
      ENDIF

      IF (ASSOCIATED(profiles_dry_tl(iprof)%co2) .AND. ASSOCIATED(profiles_tl(iprof)%co2)) THEN
        profiles_dry_tl(iprof)%co2 = 1.E+06_jprb * Mair / gas_mass(gas_id_co2) * &
                                     profiles_tl(iprof)%co2 / (1._jprb - profiles(iprof)%co2)**2
      ENDIF

      IF (ASSOCIATED(profiles_dry_tl(iprof)%co)  .AND. ASSOCIATED(profiles_tl(iprof)%co))  THEN
        profiles_dry_tl(iprof)%co = 1.E+06_jprb * Mair / gas_mass(gas_id_co) * &
                                    profiles_tl(iprof)%co / (1._jprb - profiles(iprof)%co)**2
      ENDIF

      IF (ASSOCIATED(profiles_dry_tl(iprof)%ch4) .AND. ASSOCIATED(profiles_tl(iprof)%ch4)) THEN
        profiles_dry_tl(iprof)%ch4 = 1.E+06_jprb * Mair / gas_mass(gas_id_ch4) * &
                                     profiles_tl(iprof)%ch4 / (1._jprb - profiles(iprof)%ch4)**2
      ENDIF

      IF (ASSOCIATED(profiles_dry_tl(iprof)%n2o) .AND. ASSOCIATED(profiles_tl(iprof)%n2o)) THEN
        profiles_dry_tl(iprof)%n2o = 1.E+06_jprb * Mair / gas_mass(gas_id_n2o) * &
                                     profiles_tl(iprof)%n2o / (1._jprb - profiles(iprof)%n2o)**2
      ENDIF
    ENDDO


  CASE (gas_unit_ppmv)

!    vdry_gas = vwet_gas / (1._jprb - h2o * 1.E-06_jprb)
!    factor = 1. / (1._jprb - h2o * 1.E-06_jprb)
!    factor_tl = 1.E-06_jprb * h2o_tl * factor**2
!    vdry_gas_tl = vwet_gas_tl * factor + vwet_gas * factor_tl

    DO iprof = 1, nprofiles

      factor = 1._jprb / (1._jprb - profiles(iprof)%q * 1.E-06_jprb)
      factor_tl = 1.E-06_jprb * profiles_tl(iprof)%q * factor**2
      IF (opts%rt_all%use_q2m) THEN
        factor2m = 1._jprb / (1._jprb - profiles(iprof)%s2m%q * 1.E-06_jprb)
        factor2m_tl = 1.E-06_jprb * profiles_tl(iprof)%s2m%q * factor2m**2
      ENDIF

      profiles_dry_tl(iprof)%q = &
          profiles_tl(iprof)%q * factor + profiles(iprof)%q * factor_tl
      IF (opts%rt_all%use_q2m) profiles_dry_tl(iprof)%s2m%q = &
          profiles_tl(iprof)%s2m%q * factor2m + profiles(iprof)%s2m%q * factor2m_tl

      IF (ASSOCIATED(profiles_dry_tl(iprof)%o3) .AND. ASSOCIATED(profiles_tl(iprof)%o3)) THEN
        profiles_dry_tl(iprof)%o3 = &
            profiles_tl(iprof)%o3 * factor + profiles(iprof)%o3 * factor_tl
        IF (opts%rt_all%use_q2m) profiles_dry_tl(iprof)%s2m%o = &
            profiles_tl(iprof)%s2m%o * factor2m + profiles(iprof)%s2m%o * factor2m_tl
      ENDIF

      IF (ASSOCIATED(profiles_dry_tl(iprof)%co2) .AND. ASSOCIATED(profiles_tl(iprof)%co2)) THEN
        profiles_dry_tl(iprof)%co2 = &
            profiles_tl(iprof)%co2 * factor + profiles(iprof)%co2 * factor_tl
      ENDIF

      IF (ASSOCIATED(profiles_dry_tl(iprof)%co)  .AND. ASSOCIATED(profiles_tl(iprof)%co))  THEN
        profiles_dry_tl(iprof)%co = &
            profiles_tl(iprof)%co * factor + profiles(iprof)%co * factor_tl
      ENDIF

      IF (ASSOCIATED(profiles_dry_tl(iprof)%ch4) .AND. ASSOCIATED(profiles_tl(iprof)%ch4)) THEN
        profiles_dry_tl(iprof)%ch4 = &
            profiles_tl(iprof)%ch4 * factor + profiles(iprof)%ch4 * factor_tl
      ENDIF

      IF (ASSOCIATED(profiles_dry_tl(iprof)%n2o) .AND. ASSOCIATED(profiles_tl(iprof)%n2o)) THEN
        profiles_dry_tl(iprof)%n2o = &
            profiles_tl(iprof)%n2o * factor + profiles(iprof)%n2o * factor_tl
      ENDIF
    ENDDO


  CASE DEFAULT ! Assume no conversion (i.e. compatibility mode or ppmv dry)

    DO iprof = 1, nprofiles
      profiles_dry_tl(iprof)%q     = profiles_tl(iprof)%q
      profiles_dry_tl(iprof)%s2m%q = profiles_tl(iprof)%s2m%q
      IF (ASSOCIATED(profiles_dry_tl(iprof)%o3) .AND. ASSOCIATED(profiles_tl(iprof)%o3)) THEN
        profiles_dry_tl(iprof)%o3    = profiles_tl(iprof)%o3
        profiles_dry_tl(iprof)%s2m%o = profiles_tl(iprof)%s2m%o
      ENDIF
      IF (ASSOCIATED(profiles_dry_tl(iprof)%co2) .AND. ASSOCIATED(profiles_tl(iprof)%co2)) &
            profiles_dry_tl(iprof)%co2 = profiles_tl(iprof)%co2
      IF (ASSOCIATED(profiles_dry_tl(iprof)%co)  .AND. ASSOCIATED(profiles_tl(iprof)%co))  &
            profiles_dry_tl(iprof)%co  = profiles_tl(iprof)%co
      IF (ASSOCIATED(profiles_dry_tl(iprof)%ch4) .AND. ASSOCIATED(profiles_tl(iprof)%ch4)) &
            profiles_dry_tl(iprof)%ch4 = profiles_tl(iprof)%ch4
      IF (ASSOCIATED(profiles_dry_tl(iprof)%n2o) .AND. ASSOCIATED(profiles_tl(iprof)%n2o)) &
            profiles_dry_tl(iprof)%n2o = profiles_tl(iprof)%n2o
    ENDDO

  END SELECT

  IF (LHOOK) CALL DR_HOOK('RTTOV_CONVERT_GAS_DRY_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_convert_gas_dry_tl
