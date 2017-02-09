SUBROUTINE rttov_add_prof( &
              profiles,     &
              profiles1,    &
              profiles2,    &
              lair,         &
              lground,      &
              profiles_gas)
! Description:
!   Adds two profile structures
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
  USE rttov_types, ONLY : profile_type
  USE parkind1, ONLY : jplm
!INTF_OFF
  USE parkind1, ONLY : jpim
!INTF_ON
  TYPE(profile_type), INTENT(INOUT)           :: profiles(:)
  TYPE(profile_type), INTENT(IN)              :: profiles1(SIZE(profiles))
  TYPE(profile_type), INTENT(IN)              :: profiles2(SIZE(profiles))
  LOGICAL(KIND=jplm), INTENT(IN), OPTIONAL    :: lair
  LOGICAL(KIND=jplm), INTENT(IN), OPTIONAL    :: lground
  TYPE(profile_type), INTENT(INOUT), OPTIONAL :: profiles_gas(SIZE(profiles))
!INTF_END
  LOGICAL(KIND=jplm) :: lair1
  LOGICAL(KIND=jplm) :: lground1
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: iprof
  lground1 = .TRUE.
  IF (PRESENT(lground)) lground1 = lground
  lair1 = .TRUE.
  IF (PRESENT(lair)) lair1 = lair
  nprofiles = SIZE(profiles1)
! Do not add:
!   profiles%id %date %time %gas_units
!
  DO iprof = 1, nprofiles
    IF (lground1) THEN
      profiles(iprof)%zenangle            = profiles1(iprof)%zenangle + profiles2(iprof)%zenangle
      profiles(iprof)%azangle             = profiles1(iprof)%azangle + profiles2(iprof)%azangle
      profiles(iprof)%sunzenangle         = profiles1(iprof)%sunzenangle + profiles2(iprof)%sunzenangle
      profiles(iprof)%sunazangle          = profiles1(iprof)%sunazangle + profiles2(iprof)%sunazangle
      profiles(iprof)%latitude            = profiles1(iprof)%latitude + profiles2(iprof)%latitude
      profiles(iprof)%longitude           = profiles1(iprof)%longitude + profiles2(iprof)%longitude
      profiles(iprof)%snow_frac           = profiles1(iprof)%snow_frac + profiles2(iprof)%snow_frac
      profiles(iprof)%soil_moisture       = profiles1(iprof)%soil_moisture + profiles2(iprof)%soil_moisture
      profiles(iprof)%ctp                 = profiles1(iprof)%ctp + profiles2(iprof)%ctp
      profiles(iprof)%cfraction           = profiles1(iprof)%cfraction + profiles2(iprof)%cfraction
      profiles(iprof)%elevation           = profiles1(iprof)%elevation + profiles2(iprof)%elevation
      profiles(iprof)%s2m%t               = profiles1(iprof)%s2m%t + profiles2(iprof)%s2m%t
      IF (PRESENT(profiles_gas)) THEN
        profiles_gas(iprof)%s2m%q         = profiles_gas(iprof)%s2m%q + profiles2(iprof)%s2m%q
        profiles_gas(iprof)%s2m%o         = profiles_gas(iprof)%s2m%o + profiles2(iprof)%s2m%o
      ELSE
        profiles(iprof)%s2m%q             = profiles1(iprof)%s2m%q + profiles2(iprof)%s2m%q
        profiles(iprof)%s2m%o             = profiles1(iprof)%s2m%o + profiles2(iprof)%s2m%o
      ENDIF
      profiles(iprof)%s2m%p               = profiles1(iprof)%s2m%p + profiles2(iprof)%s2m%p
      profiles(iprof)%s2m%u               = profiles1(iprof)%s2m%u + profiles2(iprof)%s2m%u
      profiles(iprof)%s2m%v               = profiles1(iprof)%s2m%v + profiles2(iprof)%s2m%v
      profiles(iprof)%s2m%wfetc           = profiles1(iprof)%s2m%wfetc + profiles2(iprof)%s2m%wfetc
      profiles(iprof)%skin%t              = profiles1(iprof)%skin%t + profiles2(iprof)%skin%t
      profiles(iprof)%skin%fastem(:)      = profiles1(iprof)%skin%fastem(:) + profiles2(iprof)%skin%fastem(:)
      profiles(iprof)%skin%salinity       = profiles1(iprof)%skin%salinity + profiles2(iprof)%skin%salinity
      profiles(iprof)%skin%foam_fraction  = profiles1(iprof)%skin%foam_fraction + &
                                            profiles2(iprof)%skin%foam_fraction
      profiles(iprof)%be                  = profiles1(iprof)%be + profiles2(iprof)%be
      profiles(iprof)%cosbk               = profiles1(iprof)%cosbk + profiles2(iprof)%cosbk
    ENDIF
    IF (lair1) THEN
      profiles(iprof)%p = profiles1(iprof)%p + profiles2(iprof)%p
      profiles(iprof)%t = profiles1(iprof)%t + profiles2(iprof)%t

      IF (PRESENT(profiles_gas)) THEN
        profiles_gas(iprof)%q            = profiles_gas(iprof)%q + profiles2(iprof)%q
        IF (ASSOCIATED(profiles2(iprof)%o3) .AND. ASSOCIATED(profiles_gas(iprof)%o3))   &
            profiles_gas(iprof)%o3       = profiles_gas(iprof)%o3 + profiles2(iprof)%o3
        IF (ASSOCIATED(profiles2(iprof)%co2) .AND. ASSOCIATED(profiles_gas(iprof)%co2)) &
            profiles_gas(iprof)%co2      = profiles_gas(iprof)%co2 + profiles2(iprof)%co2
        IF (ASSOCIATED(profiles2(iprof)%n2o) .AND. ASSOCIATED(profiles_gas(iprof)%n2o)) &
            profiles_gas(iprof)%n2o      = profiles_gas(iprof)%n2o + profiles2(iprof)%n2o
        IF (ASSOCIATED(profiles2(iprof)%co) .AND. ASSOCIATED(profiles_gas(iprof)%co))   &
            profiles_gas(iprof)%co       = profiles_gas(iprof)%co + profiles2(iprof)%co
        IF (ASSOCIATED(profiles2(iprof)%ch4) .AND. ASSOCIATED(profiles_gas(iprof)%ch4)) &
            profiles_gas(iprof)%ch4      = profiles_gas(iprof)%ch4 + profiles2(iprof)%ch4
      ELSE
        profiles(iprof)%q            = profiles1(iprof)%q + profiles2(iprof)%q
        IF (ASSOCIATED(profiles(iprof)%o3)  .AND. ASSOCIATED(profiles1(iprof)%o3) .AND. ASSOCIATED(profiles2(iprof)%o3))   &
            profiles(iprof)%o3       = profiles1(iprof)%o3 + profiles2(iprof)%o3
        IF (ASSOCIATED(profiles(iprof)%co2) .AND. ASSOCIATED(profiles1(iprof)%co2) .AND. ASSOCIATED(profiles2(iprof)%co2)) &
            profiles(iprof)%co2      = profiles1(iprof)%co2 + profiles2(iprof)%co2
        IF (ASSOCIATED(profiles(iprof)%n2o) .AND. ASSOCIATED(profiles1(iprof)%n2o) .AND. ASSOCIATED(profiles2(iprof)%n2o)) &
            profiles(iprof)%n2o      = profiles1(iprof)%n2o + profiles2(iprof)%n2o
        IF (ASSOCIATED(profiles(iprof)%co)  .AND. ASSOCIATED(profiles1(iprof)%co) .AND. ASSOCIATED(profiles2(iprof)%co))   &
            profiles(iprof)%co       = profiles1(iprof)%co + profiles2(iprof)%co
        IF (ASSOCIATED(profiles(iprof)%ch4) .AND. ASSOCIATED(profiles1(iprof)%ch4) .AND. ASSOCIATED(profiles2(iprof)%ch4)) &
            profiles(iprof)%ch4      = profiles1(iprof)%ch4 + profiles2(iprof)%ch4
      ENDIF

      IF (ASSOCIATED(profiles(iprof)%clw) .AND. ASSOCIATED(profiles1(iprof)%clw) .AND. ASSOCIATED(profiles2(iprof)%clw))   &
          profiles(iprof)%clw       = profiles1(iprof)%clw + profiles2(iprof)%clw
      IF (ASSOCIATED(profiles(iprof)%aerosols) .AND. ASSOCIATED(profiles1(iprof)%aerosols) .AND. &
          ASSOCIATED(profiles2(iprof)%aerosols))                                                 &
           profiles(iprof)%aerosols = profiles1(iprof)%aerosols + profiles2(iprof)%aerosols
      IF (ASSOCIATED(profiles(iprof)%cloud) .AND. ASSOCIATED(profiles1(iprof)%cloud) .AND.      &
          ASSOCIATED(profiles2(iprof)%cloud)         )                                          &
           profiles(iprof)%cloud    = profiles1(iprof)%cloud + profiles2(iprof)%cloud
      IF (ASSOCIATED(profiles(iprof)%cfrac) .AND. ASSOCIATED(profiles1(iprof)%cfrac) .AND.      &
          ASSOCIATED(profiles2(iprof)%cfrac)         )                                          &
           profiles(iprof)%cfrac    = profiles1(iprof)%cfrac + profiles2(iprof)%cfrac
      IF (ASSOCIATED(profiles(iprof)%icede) .AND. ASSOCIATED(profiles1(iprof)%icede) .AND.      &
          ASSOCIATED(profiles2(iprof)%icede)         )                                          &
           profiles(iprof)%icede    = profiles1(iprof)%icede + profiles2(iprof)%icede
    ENDIF
  ENDDO
END SUBROUTINE 
