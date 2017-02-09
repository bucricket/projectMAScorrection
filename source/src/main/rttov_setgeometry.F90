!
SUBROUTINE rttov_setgeometry( &
            & opts,           &
            & dosolar,        &
            & profiles,       &
            & aux,            &
            & coef,           &
            & angles,         &
            & raytracing,     &
            & profiles_dry,   &
            & do_pmc_calc)
! Description:
!   Compute all profile related viewing geometry
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
!    Copyright 2002, EUMETSAT, All Rights Reserved.
!
  USE rttov_types, ONLY :  &
       & rttov_coef,      &
       & profile_type,    &
       & profile_aux,     &
       & geometry_type,   &
       & raytracing_type, &
       & rttov_options
  USE parkind1, ONLY : jplm
!INTF_OFF
  USE rttov_const, ONLY : deg2rad
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options  ), INTENT(IN)              :: opts
  LOGICAL(jplm        ), INTENT(IN)              :: dosolar
  TYPE(profile_type   ), INTENT(IN)              :: profiles(:)
  TYPE(rttov_coef     ), INTENT(IN)              :: coef
  TYPE(profile_aux    ), INTENT(IN)   , OPTIONAL :: aux
  TYPE(geometry_type  ), INTENT(OUT)             :: angles(SIZE(profiles))
  TYPE(raytracing_type), INTENT(INOUT), OPTIONAL :: raytracing
  TYPE(profile_type   ), INTENT(IN)   , OPTIONAL :: profiles_dry(SIZE(profiles))
  LOGICAL(jplm        ), INTENT(IN)   , OPTIONAL :: do_pmc_calc
!INTF_END
#include "rttov_locpat.interface"

  INTEGER(KIND=jpim) :: i, nprofiles
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------

! JAH/PJR - 16/10/2013
! profiles%zenangle is geometrical (no refraction) zenith angle to satellite at surface
! angles%viewang is geometrical (no refraction) nadir angle to surface at satellite
! viewang calculated from zenangle by application of sine rule
! locpat calculates local path zenith angles due to refraction and curvature

!Notes on notation:
! zen  => zenith angle
!   (definition: angle at surface between view path to satellite and zenith)
! view => view angle
!   (definition: angle at the satellite between view path and nadir)
! _sq = square of given value
! _sqrt = square root of given value
! _minus1 = given value - 1
! trigonometric function abbreviations have their usual meanings
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETGEOMETRY', 0_jpim, ZHOOK_HANDLE)
  nprofiles = SIZE(profiles)
  DO i = 1, nprofiles
    angles(i)%sinzen           = SIN(profiles(i)%zenangle * deg2rad)
    angles(i)%sinzen_sq        = angles(i)%sinzen * angles(i)%sinzen
    angles(i)%coszen           = COS(profiles(i)%zenangle * deg2rad)
    angles(i)%coszen_sq        = angles(i)%coszen * angles(i)%coszen
    angles(i)%seczen           = 1.0_jprb / ABS(angles(i)%coszen)
    angles(i)%seczen_sq        = angles(i)%seczen * angles(i)%seczen
    angles(i)%seczen_sqrt      = SQRT(angles(i)%seczen)
    angles(i)%seczen_minus1    = angles(i)%seczen - 1.0_jprb
    angles(i)%seczen_minus1_sq = angles(i)%seczen_minus1 * angles(i)%seczen_minus1
    angles(i)%sinview          = angles(i)%sinzen / coef%ratoe
    angles(i)%sinview_sq       = angles(i)%sinview * angles(i)%sinview
    angles(i)%cosview_sq       = 1.0_jprb - angles(i)%sinview_sq
    angles(i)%normzen          = profiles(i)%zenangle / 60.0_jprb   ! normalized zenith angle for ISEM
    angles(i)%viewang          = ASIN(angles(i)%sinview) / deg2rad

    angles(i)%sinzen_sun       = SIN(profiles(i)%sunzenangle * deg2rad)

    angles(i)%sinlat           = SIN(profiles(i)%latitude * deg2rad)
    angles(i)%coslat           = COS(profiles(i)%latitude * deg2rad)
  ENDDO
  IF (PRESENT(raytracing) .AND. PRESENT(aux) .AND. PRESENT(profiles_dry)) THEN
    CALL rttov_locpat( &
          & opts,           &
          & dosolar,        &
          & profiles,       &
          & profiles_dry,   &
          & aux,            &
          & coef,           &
          & angles,         &
          & raytracing,     &
          & do_pmc_calc)
  ENDIF
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETGEOMETRY', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setgeometry
