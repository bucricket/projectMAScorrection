! Description:
!> @file
!!   Convert input gas profiles to ppmv wrt dry air K.
!
!> @brief
!!   Convert input gas profiles to ppmv wrt dry air K.
!!
!! @details
!!   Convert input gas profiles to ppmv wrt dry air K.
!!   This does not correspond precisely to the AD for
!!   reasons of efficiency: multiplicative factors which
!!   only vary per-profile are recalculated only when the
!!   profile number changes in the loop over profiles_k.
!!
!! @param[in]     opts             options to configure the simulations
!! @param[in]     profiles         profiles structure containing input gas profiles
!! @param[in,out] profiles_k       profiles structure containing profile gas perturbations
!! @param[in]     profiles_dry     profiles structure containing converted profiles on exit
!! @param[in]     profiles_dry_k   profiles structure containing converted profile perturbations
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
SUBROUTINE rttov_convert_gas_dry_k(opts, chanprof, profiles, profiles_k, profiles_dry, profiles_dry_k)

  USE rttov_types, ONLY : &
      rttov_options,      &
      rttov_chanprof,     &
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

  TYPE(rttov_options),  INTENT(IN)    :: opts
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof(:)
  TYPE(profile_type),   INTENT(IN)    :: profiles(:)
  TYPE(profile_type),   INTENT(INOUT) :: profiles_k(SIZE(chanprof))
  TYPE(profile_type),   INTENT(IN)    :: profiles_dry(SIZE(profiles))
  TYPE(profile_type),   INTENT(IN)    :: profiles_dry_k(SIZE(chanprof))
!INTF_END

  INTEGER(jpim) :: ichan, nchannels, prof, lastprof
  REAL(jprb)    :: factor(profiles(1)%nlevels), factor2m
  REAL(jprb)    :: factor_k(profiles(1)%nlevels), factor2m_k
  REAL(jprb)    :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CONVERT_GAS_DRY_K', 0_jpim, ZHOOK_HANDLE)

! As noted above multiplicative factors which only vary per-profile are
! recalculated only when the profile number changes in the loop over
! profiles_k. This is much more efficient than re-calculating them in
! every loop, but it means this code differs superficially to the AD.

  nchannels = SIZE(profiles_k)

  SELECT CASE (profiles(1)%gas_units)

  CASE (gas_unit_specconc)

!    vdry_gas = qwet_gas / (1._jprb - qwet_gas) * 1.E+06_jprb * Mair / gas_mass(gas_id)
!    vdry_gas_tl = qwet_gas_tl / (1._jprb - qwet_gas)**2 * 1.E+06_jprb * Mair / gas_mass(gas_id)
!    qwet_gas_k = qwet_gas_k + vdry_gas_k / (1._jprb - qwet_gas)**2 * 1.E+06_jprb * Mair / gas_mass(gas_id)

    lastprof = -1
    DO ichan = 1, nchannels
      prof = chanprof(ichan)%prof

      IF (lastprof /= prof) THEN
        factor = 1.E+06_jprb * Mair / gas_mass(gas_id_watervapour) / &
                 (1._jprb - profiles(prof)%q)**2
        IF (opts%rt_all%use_q2m) factor2m = 1.E+06_jprb * Mair / gas_mass(gas_id_watervapour) / &
                                            (1._jprb - profiles(prof)%s2m%q)**2
        lastprof = prof
      ENDIF

      profiles_k(ichan)%q = profiles_k(ichan)%q + factor * profiles_dry_k(ichan)%q
      IF (opts%rt_all%use_q2m) profiles_k(ichan)%s2m%q = profiles_k(ichan)%s2m%q + &
                                                         factor2m * profiles_dry_k(ichan)%s2m%q
    ENDDO

    IF (ASSOCIATED(profiles_dry_k(1)%o3) .AND. ASSOCIATED(profiles_k(1)%o3)) THEN
      lastprof = -1
      DO ichan = 1, nchannels
        prof = chanprof(ichan)%prof

        IF (lastprof /= prof) THEN
          factor = 1.E+06_jprb * Mair / gas_mass(gas_id_ozone) / &
                   (1._jprb - profiles(prof)%o3)**2
          IF (opts%rt_all%use_q2m) factor2m = 1.E+06_jprb * Mair / gas_mass(gas_id_ozone) / &
                                              (1._jprb - profiles(prof)%s2m%o)**2
          lastprof = prof
        ENDIF

        profiles_k(ichan)%o3 = profiles_k(ichan)%o3 + factor * profiles_dry_k(ichan)%o3
        IF (opts%rt_all%use_q2m) profiles_k(ichan)%s2m%o = profiles_k(ichan)%s2m%o + &
                                                           factor2m * profiles_dry_k(ichan)%s2m%o
      ENDDO
    ENDIF

    IF (ASSOCIATED(profiles_dry_k(1)%co2) .AND. ASSOCIATED(profiles_k(1)%co2)) THEN
      lastprof = -1
      DO ichan = 1, nchannels
        prof = chanprof(ichan)%prof
        IF (lastprof /= prof) THEN
          factor = 1.E+06_jprb * Mair / gas_mass(gas_id_co2) / &
                   (1._jprb - profiles(prof)%co2)**2
          lastprof = prof
        ENDIF

        profiles_k(ichan)%co2 = profiles_k(ichan)%co2 + factor * profiles_dry_k(ichan)%co2
      ENDDO
    ENDIF

    IF (ASSOCIATED(profiles_dry_k(1)%co)  .AND. ASSOCIATED(profiles_k(1)%co))  THEN
      lastprof = -1
      DO ichan = 1, nchannels
        prof = chanprof(ichan)%prof
        IF (lastprof /= prof) THEN
          factor = 1.E+06_jprb * Mair / gas_mass(gas_id_co) / &
                   (1._jprb - profiles(prof)%co)**2
          lastprof = prof
        ENDIF

        profiles_k(ichan)%co = profiles_k(ichan)%co + factor * profiles_dry_k(ichan)%co
      ENDDO
    ENDIF

    IF (ASSOCIATED(profiles_dry_k(1)%ch4) .AND. ASSOCIATED(profiles_k(1)%ch4)) THEN
      lastprof = -1
      DO ichan = 1, nchannels
        prof = chanprof(ichan)%prof
        IF (lastprof /= prof) THEN
          factor = 1.E+06_jprb * Mair / gas_mass(gas_id_ch4) / &
                   (1._jprb - profiles(prof)%ch4)**2
          lastprof = prof
        ENDIF

        profiles_k(ichan)%ch4 = profiles_k(ichan)%ch4 + factor * profiles_dry_k(ichan)%ch4
      ENDDO
    ENDIF

    IF (ASSOCIATED(profiles_dry_k(1)%n2o) .AND. ASSOCIATED(profiles_k(1)%n2o)) THEN
      lastprof = -1
      DO ichan = 1, nchannels
        prof = chanprof(ichan)%prof
        IF (lastprof /= prof) THEN
          factor = 1.E+06_jprb * Mair / gas_mass(gas_id_n2o) / &
                   (1._jprb - profiles(prof)%n2o)**2
          lastprof = prof
        ENDIF

        profiles_k(ichan)%n2o = profiles_k(ichan)%n2o + factor * profiles_dry_k(ichan)%n2o
      ENDDO
    ENDIF


  CASE (gas_unit_ppmv)

!    vdry_gas = vwet_gas / (1._jprb - h2o * 1.E-06_jprb)
!    factor = 1. / (1._jprb - h2o * 1.E-06_jprb)
!    factor_tl = 1.E-06_jprb * h2o_tl * factor**2
!    vdry_gas_tl = vwet_gas_tl * factor + vwet_gas * factor_tl

!    vwet_gas_k = vwet_gas_k + factor * vdry_gas_k
!    factor_k = factor_k + vwet_gas * vdry_gas_k
!    h2o_k = h2o_k + 1.E-06_jprb * factor_k * factor**2

    lastprof = -1
    DO ichan = 1, nchannels
      prof = chanprof(ichan)%prof

      IF (lastprof /= prof) THEN
        factor = 1._jprb / (1._jprb - profiles(prof)%q * 1.E-06_jprb)
        IF (opts%rt_all%use_q2m) THEN
          factor2m = 1._jprb / (1._jprb - profiles(prof)%s2m%q * 1.E-06_jprb)
        ENDIF
        lastprof = prof
      ENDIF

      profiles_k(ichan)%q = profiles_k(ichan)%q + factor * profiles_dry_k(ichan)%q
      ! NB First use of factor_k
      factor_k = profiles(prof)%q * profiles_dry_k(ichan)%q

      IF (opts%rt_all%use_q2m) THEN
        profiles_k(ichan)%s2m%q = profiles_k(ichan)%s2m%q + factor2m * profiles_dry_k(ichan)%s2m%q
        ! NB First use of factor2m_k
        factor2m_k = profiles(prof)%s2m%q * profiles_dry_k(ichan)%s2m%q
      ENDIF

      IF (ASSOCIATED(profiles_dry_k(ichan)%o3) .AND. ASSOCIATED(profiles_k(ichan)%o3)) THEN
        profiles_k(ichan)%o3 = profiles_k(ichan)%o3 + factor * profiles_dry_k(ichan)%o3
        factor_k = factor_k + profiles(prof)%o3 * profiles_dry_k(ichan)%o3

        IF (opts%rt_all%use_q2m) THEN
          profiles_k(ichan)%s2m%o = profiles_k(ichan)%s2m%o + factor2m * profiles_dry_k(ichan)%s2m%o
          factor2m_k = factor2m_k + profiles(prof)%s2m%o * profiles_dry_k(ichan)%s2m%o
        ENDIF
      ENDIF

      IF (ASSOCIATED(profiles_dry_k(ichan)%co2) .AND. ASSOCIATED(profiles_k(ichan)%co2)) THEN
        profiles_k(ichan)%co2 = profiles_k(ichan)%co2 + factor * profiles_dry_k(ichan)%co2
        factor_k = factor_k + profiles(prof)%co2 * profiles_dry_k(ichan)%co2
      ENDIF

      IF (ASSOCIATED(profiles_dry_k(ichan)%co)  .AND. ASSOCIATED(profiles_k(ichan)%co))  THEN
        profiles_k(ichan)%co = profiles_k(ichan)%co + factor * profiles_dry_k(ichan)%co
        factor_k = factor_k + profiles(prof)%co * profiles_dry_k(ichan)%co
      ENDIF

      IF (ASSOCIATED(profiles_dry_k(ichan)%ch4) .AND. ASSOCIATED(profiles_k(ichan)%ch4)) THEN
        profiles_k(ichan)%ch4 = profiles_k(ichan)%ch4 + factor * profiles_dry_k(ichan)%ch4
        factor_k = factor_k + profiles(prof)%ch4 * profiles_dry_k(ichan)%ch4
      ENDIF

      IF (ASSOCIATED(profiles_dry_k(ichan)%n2o) .AND. ASSOCIATED(profiles_k(ichan)%n2o)) THEN
        profiles_k(ichan)%n2o = profiles_k(ichan)%n2o + factor * profiles_dry_k(ichan)%n2o
        factor_k = factor_k + profiles(prof)%n2o * profiles_dry_k(ichan)%n2o
      ENDIF

      profiles_k(ichan)%q = profiles_k(ichan)%q + 1.E-06_jprb * factor_k * factor**2
      IF (opts%rt_all%use_q2m) THEN
        profiles_k(ichan)%s2m%q = profiles_k(ichan)%s2m%q + 1.E-06_jprb * factor2m_k * factor2m**2
      ENDIF
    ENDDO


  CASE DEFAULT ! Assume no conversion (i.e. compatibility mode or ppmv dry)

    DO ichan = 1, nchannels
      prof = chanprof(ichan)%prof

      profiles_k(ichan)%q = profiles_k(ichan)%q + profiles_dry_k(ichan)%q
      profiles_k(ichan)%s2m%q = profiles_k(ichan)%s2m%q + profiles_dry_k(ichan)%s2m%q
      IF (ASSOCIATED(profiles_dry_k(ichan)%o3) .AND. ASSOCIATED(profiles_k(ichan)%o3)) THEN
        profiles_k(ichan)%o3 = profiles_k(ichan)%o3 + profiles_dry_k(ichan)%o3
        profiles_k(ichan)%s2m%o = profiles_k(ichan)%s2m%o + profiles_dry_k(ichan)%s2m%o
      ENDIF
      IF (ASSOCIATED(profiles_dry_k(ichan)%co2) .AND. ASSOCIATED(profiles_k(ichan)%co2)) &
          profiles_k(ichan)%co2 = profiles_k(ichan)%co2 + profiles_dry_k(ichan)%co2
      IF (ASSOCIATED(profiles_dry_k(ichan)%co)  .AND. ASSOCIATED(profiles_k(ichan)%co))  &
          profiles_k(ichan)%co  = profiles_k(ichan)%co  + profiles_dry_k(ichan)%co
      IF (ASSOCIATED(profiles_dry_k(ichan)%ch4) .AND. ASSOCIATED(profiles_k(ichan)%ch4)) &
          profiles_k(ichan)%ch4 = profiles_k(ichan)%ch4 + profiles_dry_k(ichan)%ch4
      IF (ASSOCIATED(profiles_dry_k(ichan)%n2o) .AND. ASSOCIATED(profiles_k(ichan)%n2o)) &
          profiles_k(ichan)%n2o = profiles_k(ichan)%n2o + profiles_dry_k(ichan)%n2o
    ENDDO

  END SELECT

  IF (LHOOK) CALL DR_HOOK('RTTOV_CONVERT_GAS_DRY_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_convert_gas_dry_k
