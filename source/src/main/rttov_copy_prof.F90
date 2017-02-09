SUBROUTINE rttov_copy_prof( &
              profiles1,    &
              profiles2,    &
              larray,       &
              lscalar,      &
              profiles_gas)
! Description:
!   Copy profile structure
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
  TYPE(profile_type), INTENT(INOUT)        :: profiles1(:)
  TYPE(profile_type), INTENT(IN)           :: profiles2(SIZE(profiles1))
  LOGICAL(KIND=jplm), INTENT(IN), OPTIONAL :: larray
  LOGICAL(KIND=jplm), INTENT(IN), OPTIONAL :: lscalar
  TYPE(profile_type), INTENT(IN), OPTIONAL :: profiles_gas(SIZE(profiles1))
!INTF_END
  LOGICAL(KIND=jplm) :: larray1
  LOGICAL(KIND=jplm) :: lscalar1
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: iprof
  lscalar1 = .TRUE.
  IF (PRESENT(lscalar)) lscalar1 = lscalar
  larray1 = .TRUE.
  IF (PRESENT(larray)) larray1 = larray
  nprofiles = SIZE(profiles1)
  DO iprof = 1, nprofiles
    IF (lscalar1) THEN
      profiles1(iprof)%id            = profiles2(iprof)%id
      profiles1(iprof)%date          = profiles2(iprof)%date
      profiles1(iprof)%time          = profiles2(iprof)%time
      profiles1(iprof)%gas_units     = profiles2(iprof)%gas_units
      profiles1(iprof)%zenangle      = profiles2(iprof)%zenangle
      profiles1(iprof)%azangle       = profiles2(iprof)%azangle
      profiles1(iprof)%sunzenangle   = profiles2(iprof)%sunzenangle
      profiles1(iprof)%sunazangle    = profiles2(iprof)%sunazangle
      profiles1(iprof)%latitude      = profiles2(iprof)%latitude
      profiles1(iprof)%longitude     = profiles2(iprof)%longitude
      profiles1(iprof)%snow_frac     = profiles2(iprof)%snow_frac
      profiles1(iprof)%soil_moisture = profiles2(iprof)%soil_moisture
      profiles1(iprof)%ctp           = profiles2(iprof)%ctp
      profiles1(iprof)%cfraction     = profiles2(iprof)%cfraction
      profiles1(iprof)%ish           = profiles2(iprof)%ish
      profiles1(iprof)%idg           = profiles2(iprof)%idg
      profiles1(iprof)%elevation     = profiles2(iprof)%elevation
      profiles1(iprof)%s2m           = profiles2(iprof)%s2m
      IF (PRESENT(profiles_gas)) THEN
        profiles1(iprof)%s2m%q = profiles_gas(iprof)%s2m%q
        profiles1(iprof)%s2m%o = profiles_gas(iprof)%s2m%o
      ENDIF
      profiles1(iprof)%skin          = profiles2(iprof)%skin
      profiles1(iprof)%be            = profiles2(iprof)%be
      profiles1(iprof)%cosbk         = profiles2(iprof)%cosbk
    ENDIF
    IF (larray1) THEN
      profiles1(iprof)%nlevels = profiles2(iprof)%nlevels
      profiles1(iprof)%nlayers = profiles2(iprof)%nlayers
      profiles1(iprof)%p       = profiles2(iprof)%p
      profiles1(iprof)%t       = profiles2(iprof)%t
      IF (PRESENT(profiles_gas)) THEN
        profiles1(iprof)%q             = profiles_gas(iprof)%q
        IF (ASSOCIATED(profiles1(iprof)%o3) .AND. ASSOCIATED(profiles_gas(iprof)%o3))    &
            profiles1(iprof)%o3       = profiles_gas(iprof)%o3
        IF (ASSOCIATED(profiles1(iprof)%co2) .AND. ASSOCIATED(profiles_gas(iprof)%co2))  &
            profiles1(iprof)%co2      = profiles_gas(iprof)%co2
        IF (ASSOCIATED(profiles1(iprof)%n2o) .AND. ASSOCIATED(profiles_gas(iprof)%n2o))  &
            profiles1(iprof)%n2o      = profiles_gas(iprof)%n2o
        IF (ASSOCIATED(profiles1(iprof)%co) .AND. ASSOCIATED(profiles_gas(iprof)%co))    &
            profiles1(iprof)%co       = profiles_gas(iprof)%co
        IF (ASSOCIATED(profiles1(iprof)%ch4) .AND. ASSOCIATED(profiles_gas(iprof)%ch4))  &
            profiles1(iprof)%ch4      = profiles_gas(iprof)%ch4
      ELSE
        profiles1(iprof)%q             = profiles2(iprof)%q
        IF (ASSOCIATED(profiles1(iprof)%o3) .AND. ASSOCIATED(profiles2(iprof)%o3))       &
            profiles1(iprof)%o3       = profiles2(iprof)%o3
        IF (ASSOCIATED(profiles1(iprof)%co2) .AND. ASSOCIATED(profiles2(iprof)%co2))     &
            profiles1(iprof)%co2      = profiles2(iprof)%co2
        IF (ASSOCIATED(profiles1(iprof)%n2o) .AND. ASSOCIATED(profiles2(iprof)%n2o))     &
            profiles1(iprof)%n2o      = profiles2(iprof)%n2o
        IF (ASSOCIATED(profiles1(iprof)%co) .AND. ASSOCIATED(profiles2(iprof)%co))       &
            profiles1(iprof)%co       = profiles2(iprof)%co
        IF (ASSOCIATED(profiles1(iprof)%ch4) .AND. ASSOCIATED(profiles2(iprof)%ch4))     &
            profiles1(iprof)%ch4      = profiles2(iprof)%ch4
      ENDIF
      IF (ASSOCIATED(profiles1(iprof)%clw) .AND. ASSOCIATED(profiles2(iprof)%clw))           &
          profiles1(iprof)%clw      = profiles2(iprof)%clw
      IF (ASSOCIATED(profiles1(iprof)%aerosols) .AND. ASSOCIATED(profiles2(iprof)%aerosols)) &
          profiles1(iprof)%aerosols = profiles2(iprof)%aerosols
      IF (ASSOCIATED(profiles1(iprof)%cloud) .AND. ASSOCIATED(profiles2(iprof)%cloud))       &
          profiles1(iprof)%cloud    = profiles2(iprof)%cloud
      IF (ASSOCIATED(profiles1(iprof)%cfrac) .AND. ASSOCIATED(profiles2(iprof)%cfrac))       &
          profiles1(iprof)%cfrac    = profiles2(iprof)%cfrac
      IF (ASSOCIATED(profiles1(iprof)%icede) .AND. ASSOCIATED(profiles2(iprof)%icede))       &
          profiles1(iprof)%icede    = profiles2(iprof)%icede
    ENDIF
  ENDDO
END SUBROUTINE 
