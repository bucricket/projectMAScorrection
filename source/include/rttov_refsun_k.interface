Interface
SUBROUTINE rttov_refsun_k( &
            & chanprof,     &
            & profiles,     &
            & profiles_k,   &
            & coef,         &
            & aux,          &
            & sunglint,     &
            & sunglint_k,   &
            & raytracing,   &
            & raytracing_k)
  USE rttov_types, ONLY :  &
       & rttov_chanprof,  &
       & profile_aux,     &
       & profile_Type,    &
       & raytracing_type, &
       & sunglint_type,   &
       & rttov_coef
  IMPLICIT NONE
  TYPE(rttov_chanprof ), INTENT(IN)    :: chanprof  (:)
  TYPE(profile_type   ), INTENT(IN)    :: profiles  (:)
  TYPE(profile_type   ), INTENT(INOUT) :: profiles_k(size(chanprof))
  TYPE(profile_aux    ), INTENT(IN)    :: aux
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing_k
  TYPE(sunglint_type  ), INTENT(IN)    :: sunglint
  TYPE(sunglint_type  ), INTENT(INOUT) :: sunglint_k
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
End Subroutine
End Interface
