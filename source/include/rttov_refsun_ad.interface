Interface
SUBROUTINE rttov_refsun_ad( &
            & profiles,      &
            & profiles_ad,   &
            & coef,          &
            & aux,           &
            & sunglint,      &
            & sunglint_ad,   &
            & raytracing,    &
            & raytracing_ad)
  USE rttov_types, ONLY :  &
       & profile_aux,     &
       & profile_Type,    &
       & raytracing_type, &
       & sunglint_type,   &
       & rttov_coef
  IMPLICIT NONE
  TYPE(profile_type   ), INTENT(IN)    :: profiles   (:)
  TYPE(profile_type   ), INTENT(INOUT) :: profiles_ad(size(profiles))
  TYPE(profile_aux    ), INTENT(IN)    :: aux
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing_ad
  TYPE(sunglint_type  ), INTENT(IN)    :: sunglint
  TYPE(sunglint_type  ), INTENT(INOUT) :: sunglint_ad
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
End Subroutine
End Interface
