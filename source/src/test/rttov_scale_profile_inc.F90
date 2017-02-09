
SUBROUTINE rttov_scale_profile_inc(profiles_inc, factor)
! Description:
!   Scales profile increment
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
! Method:
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
  USE parkind1, ONLY : jprb
  USE rttov_types, ONLY : profile_Type
!INTF_OFF
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE
  TYPE(profile_type), INTENT(INOUT)  :: profiles_inc(:)
  REAL(KIND=jprb), INTENT(IN)        :: factor
!INTF_END
  INTEGER(KIND=jpim) :: j, nprofiles

  nprofiles = size(profiles_inc)

  DO j = 1, nprofiles
! increments for atmospheric variables

    profiles_inc(j)%p = profiles_inc(j)%p * factor
    profiles_inc(j)%t = profiles_inc(j)%t * factor
    profiles_inc(j)%q = profiles_inc(j)%q * factor
    IF (associated(profiles_inc(j)%o3) ) profiles_inc(j)%o3  = profiles_inc(j)%o3 * factor
    IF (associated(profiles_inc(j)%co2)) profiles_inc(j)%co2 = profiles_inc(j)%co2 * factor
    IF (associated(profiles_inc(j)%co) ) profiles_inc(j)%co  = profiles_inc(j)%co * factor
    IF (associated(profiles_inc(j)%n2o)) profiles_inc(j)%n2o = profiles_inc(j)%n2o * factor
    IF (associated(profiles_inc(j)%ch4)) profiles_inc(j)%ch4 = profiles_inc(j)%ch4 * factor

    IF (associated(profiles_inc(j)%aerosols)) THEN
      profiles_inc(j)%aerosols =  profiles_inc(j)%aerosols * factor
    ENDIF

    IF (associated(profiles_inc(j)%cfrac)) THEN
      profiles_inc(j)%cfrac = profiles_inc(j)%cfrac * factor
    ENDIF

    IF (associated(profiles_inc(j)%cloud)) THEN
      profiles_inc(j)%cloud = profiles_inc(j)%cloud * factor
      IF (associated(profiles_inc(j)%icede)) THEN
        profiles_inc(j)%icede = profiles_inc(j)%icede * factor
      ENDIF
    ENDIF

    IF (associated(profiles_inc(j)%clw)) profiles_inc(j)%clw = profiles_inc(j)%clw * factor

    profiles_inc(j)%s2m%t              = profiles_inc(j)%s2m%t * factor
    profiles_inc(j)%s2m%q              = profiles_inc(j)%s2m%q * factor
    profiles_inc(j)%s2m%p              = profiles_inc(j)%s2m%p * factor
    profiles_inc(j)%s2m%u              = profiles_inc(j)%s2m%u * factor
    profiles_inc(j)%s2m%v              = profiles_inc(j)%s2m%v * factor
    profiles_inc(j)%s2m%o              = profiles_inc(j)%s2m%o * factor
    profiles_inc(j)%s2m%wfetc          = profiles_inc(j)%s2m%wfetc * factor
    profiles_inc(j)%skin%t             = profiles_inc(j)%skin%t * factor
    profiles_inc(j)%skin%fastem        = profiles_inc(j)%skin%fastem * factor
    profiles_inc(j)%skin%salinity      = profiles_inc(j)%skin%salinity * factor
    profiles_inc(j)%skin%foam_fraction = profiles_inc(j)%skin%foam_fraction * factor
    profiles_inc(j)%ctp                = profiles_inc(j)%ctp * factor
    profiles_inc(j)%cfraction          = profiles_inc(j)%cfraction * factor
    profiles_inc(j)%zenangle           = profiles_inc(j)%zenangle * factor
    profiles_inc(j)%azangle            = profiles_inc(j)%azangle * factor
    profiles_inc(j)%sunzenangle        = profiles_inc(j)%sunzenangle * factor
    profiles_inc(j)%sunazangle         = profiles_inc(j)%sunazangle * factor
    profiles_inc(j)%elevation          = profiles_inc(j)%elevation * factor
    profiles_inc(j)%latitude           = profiles_inc(j)%latitude * factor
    profiles_inc(j)%longitude          = profiles_inc(j)%longitude * factor
    profiles_inc(j)%snow_frac          = profiles_inc(j)%snow_frac * factor
    profiles_inc(j)%soil_moisture      = profiles_inc(j)%soil_moisture * factor
    profiles_inc(j)%be                 = profiles_inc(j)%be * factor
    profiles_inc(j)%cosbk              = profiles_inc(j)%cosbk * factor
  ENDDO

END SUBROUTINE

