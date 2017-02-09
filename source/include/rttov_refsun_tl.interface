Interface
SUBROUTINE rttov_refsun_tl( &
            & profiles,      &
            & profiles_tl,   &
            & coef,          &
            & aux,           &
            & sunglint,      &
            & sunglint_tl,   &
            & raytracing,    &
            & raytracing_tl)
  USE rttov_types, ONLY :  &
       & profile_aux,     &
       & profile_Type,    &
       & raytracing_type, &
       & sunglint_type,   &
       & rttov_coef
  IMPLICIT NONE
  TYPE(profile_type   ), INTENT(IN)    :: profiles   (:)
  TYPE(profile_type   ), INTENT(IN)    :: profiles_tl(size(profiles))
  TYPE(profile_aux    ), INTENT(IN)    :: aux
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(IN)    :: raytracing_tl
  TYPE(sunglint_type  ), INTENT(IN)    :: sunglint
  TYPE(sunglint_type  ), INTENT(INOUT) :: sunglint_tl
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
End Subroutine
End Interface
