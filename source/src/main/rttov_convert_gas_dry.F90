! Description:
!> @file
!!   Convert input gas profiles to ppmv wrt dry air.
!
!> @brief
!!   Convert input gas profiles to ppmv wrt dry air.
!!
!! @details
!!   Convert input gas profiles to ppmv wrt dry air.
!!
!! @param[in]     opts           options to configure the simulations
!! @param[in]     profiles       profiles structure containing input gas profiles
!! @param[in,out] profiles_dry   profiles structure containing converted profiles on exit
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
SUBROUTINE rttov_convert_gas_dry(opts, profiles, profiles_dry)

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
  TYPE(profile_type),  INTENT(INOUT) :: profiles_dry(SIZE(profiles))
!INTF_END

  INTEGER(jpim) :: iprof, nprofiles
  REAL(jprb)    :: factor(profiles(1)%nlevels), factor2m
  REAL(jprb)    :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CONVERT_GAS_DRY', 0_jpim, ZHOOK_HANDLE)

! These conversion formulae could have been put into separate subroutines
! or functions, but the conversions (and especially the TL/AD/K) are more
! efficiently coded following the pattern below.

  nprofiles = SIZE(profiles)

  SELECT CASE (profiles(1)%gas_units)

  CASE (gas_unit_specconc)

!     vdry_gas = qwet_gas / (1._jprb - qwet_gas) * 1.E+06_jprb * Mair / gas_mass(gas_id)

    DO iprof = 1, nprofiles
      profiles_dry(iprof)%q = profiles(iprof)%q / (1._jprb - profiles(iprof)%q) * &
                              1.E+06_jprb * Mair / gas_mass(gas_id_watervapour)
      IF (opts%rt_all%use_q2m) profiles_dry(iprof)%s2m%q = &
                              profiles(iprof)%s2m%q / (1._jprb - profiles(iprof)%s2m%q) * &
                              1.E+06_jprb * Mair / gas_mass(gas_id_watervapour)

      IF (ASSOCIATED(profiles_dry(iprof)%o3) .AND. ASSOCIATED(profiles(iprof)%o3)) THEN
        profiles_dry(iprof)%o3 = profiles(iprof)%o3 / (1._jprb - profiles(iprof)%o3) * &
                                 1.E+06_jprb * Mair / gas_mass(gas_id_ozone)
        IF (opts%rt_all%use_q2m) profiles_dry(iprof)%s2m%o = &
                                 profiles(iprof)%s2m%o / (1._jprb - profiles(iprof)%s2m%o) * &
                                 1.E+06_jprb * Mair / gas_mass(gas_id_ozone)
      ENDIF

      IF (ASSOCIATED(profiles_dry(iprof)%co2) .AND. ASSOCIATED(profiles(iprof)%co2)) THEN
        profiles_dry(iprof)%co2 = profiles(iprof)%co2 / (1._jprb - profiles(iprof)%co2) * &
                                  1.E+06_jprb * Mair / gas_mass(gas_id_co2)
      ENDIF

      IF (ASSOCIATED(profiles_dry(iprof)%co)  .AND. ASSOCIATED(profiles(iprof)%co)) THEN
        profiles_dry(iprof)%co = profiles(iprof)%co / (1._jprb - profiles(iprof)%co) * &
                                 1.E+06_jprb * Mair / gas_mass(gas_id_co)
      ENDIF

      IF (ASSOCIATED(profiles_dry(iprof)%ch4) .AND. ASSOCIATED(profiles(iprof)%ch4)) THEN
        profiles_dry(iprof)%ch4 = profiles(iprof)%ch4 / (1._jprb - profiles(iprof)%ch4) * &
                                  1.E+06_jprb * Mair / gas_mass(gas_id_ch4)
      ENDIF

      IF (ASSOCIATED(profiles_dry(iprof)%n2o) .AND. ASSOCIATED(profiles(iprof)%n2o)) THEN
        profiles_dry(iprof)%n2o = profiles(iprof)%n2o / (1._jprb - profiles(iprof)%n2o) * &
                                  1.E+06_jprb * Mair / gas_mass(gas_id_n2o)
      ENDIF
    ENDDO


  CASE (gas_unit_ppmv)

!     vdry_gas = vwet_gas / (1._jprb - h2o * 1.E-06_jprb)

    DO iprof = 1, nprofiles

      factor = 1._jprb / (1._jprb - profiles(iprof)%q * 1.E-06_jprb)
      IF (opts%rt_all%use_q2m) factor2m = 1._jprb / (1._jprb - profiles(iprof)%s2m%q * 1.E-06_jprb)

      profiles_dry(iprof)%q = factor * profiles(iprof)%q
      IF (opts%rt_all%use_q2m) profiles_dry(iprof)%s2m%q = factor2m * profiles(iprof)%s2m%q

      IF (ASSOCIATED(profiles_dry(iprof)%o3) .AND. ASSOCIATED(profiles(iprof)%o3)) THEN
        profiles_dry(iprof)%o3 = factor * profiles(iprof)%o3
        IF (opts%rt_all%use_q2m) profiles_dry(iprof)%s2m%o = factor2m * profiles(iprof)%s2m%o
      ENDIF

      IF (ASSOCIATED(profiles_dry(iprof)%co2) .AND. ASSOCIATED(profiles(iprof)%co2)) THEN
        profiles_dry(iprof)%co2 = factor * profiles(iprof)%co2
      ENDIF

      IF (ASSOCIATED(profiles_dry(iprof)%co)  .AND. ASSOCIATED(profiles(iprof)%co)) THEN
        profiles_dry(iprof)%co = factor * profiles(iprof)%co
      ENDIF

      IF (ASSOCIATED(profiles_dry(iprof)%ch4) .AND. ASSOCIATED(profiles(iprof)%ch4)) THEN
        profiles_dry(iprof)%ch4 = factor * profiles(iprof)%ch4
      ENDIF

      IF (ASSOCIATED(profiles_dry(iprof)%n2o) .AND. ASSOCIATED(profiles(iprof)%n2o)) THEN
        profiles_dry(iprof)%n2o = factor * profiles(iprof)%n2o
      ENDIF
    ENDDO


  CASE DEFAULT ! Assume no conversion (i.e. compatibility mode or ppmv dry)

    DO iprof = 1, nprofiles
      profiles_dry(iprof)%q = profiles(iprof)%q
      IF (opts%rt_all%use_q2m) profiles_dry(iprof)%s2m%q = profiles(iprof)%s2m%q
      IF (ASSOCIATED(profiles_dry(iprof)%o3) .AND. ASSOCIATED(profiles(iprof)%o3)) THEN
        profiles_dry(iprof)%o3 = profiles(iprof)%o3
        IF (opts%rt_all%use_q2m) profiles_dry(iprof)%s2m%o = profiles(iprof)%s2m%o
      ENDIF
      IF (ASSOCIATED(profiles_dry(iprof)%co2) .AND. ASSOCIATED(profiles(iprof)%co2)) &
            profiles_dry(iprof)%co2 = profiles(iprof)%co2
      IF (ASSOCIATED(profiles_dry(iprof)%co)  .AND. ASSOCIATED(profiles(iprof)%co))  &
            profiles_dry(iprof)%co  = profiles(iprof)%co
      IF (ASSOCIATED(profiles_dry(iprof)%ch4) .AND. ASSOCIATED(profiles(iprof)%ch4)) &
            profiles_dry(iprof)%ch4 = profiles(iprof)%ch4
      IF (ASSOCIATED(profiles_dry(iprof)%n2o) .AND. ASSOCIATED(profiles(iprof)%n2o)) &
            profiles_dry(iprof)%n2o = profiles(iprof)%n2o
    ENDDO

  END SELECT

  IF (LHOOK) CALL DR_HOOK('RTTOV_CONVERT_GAS_DRY', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_convert_gas_dry
