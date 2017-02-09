!
SUBROUTINE rttov_setgeometry_ad( &
            & opts,           &
            & dosolar,        &
            & profiles,       &
            & profiles_ad,    &
            & aux,            &
            & coef,           &
            & angles,         &
            & raytracing,     &
            & raytracing_ad,  &
            & profiles_dry,   &
            & profiles_dry_ad,&
            & do_pmc_calc)
! Description:
!   Compute all profile related viewing geometry AD.
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
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options  ), INTENT(IN)    :: opts
  LOGICAL(jplm        ), INTENT(IN)    :: dosolar
  TYPE(profile_type   ), INTENT(IN)    :: profiles(:)
  TYPE(profile_type   ), INTENT(INOUT) :: profiles_ad(SIZE(profiles))
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  TYPE(profile_aux    ), INTENT(IN)    :: aux
  TYPE(geometry_type  ), INTENT(IN)    :: angles(SIZE(profiles))
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing_ad
  TYPE(profile_type   ), INTENT(IN)    :: profiles_dry(SIZE(profiles))
  TYPE(profile_type   ), INTENT(INOUT) :: profiles_dry_ad(SIZE(profiles))
  LOGICAL(jplm        ), INTENT(IN), OPTIONAL :: do_pmc_calc
!INTF_END
#include "rttov_locpat_ad.interface"

  REAL   (KIND=jprb) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('RTTOV_SETGEOMETRY_AD', 0_jpim, ZHOOK_HANDLE)
  CALL rttov_locpat_ad( &
        & opts,            &
        & dosolar,         &
        & profiles,        &
        & profiles_ad,     &
        & profiles_dry,    &
        & profiles_dry_ad, &
        & aux,             &
        & coef,            &
        & angles,          &
        & raytracing,      &
        & raytracing_ad,   &
        & do_pmc_calc)
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETGEOMETRY_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setgeometry_ad
