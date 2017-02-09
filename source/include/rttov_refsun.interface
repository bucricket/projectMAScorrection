Interface
SUBROUTINE rttov_refsun( &
            & profiles,   &
            & coef,       &
            & aux,        &
            & sunglint,   &
            & raytracing)
  USE rttov_types, ONLY :  &
       & profile_aux,     &
       & profile_Type,    &
       & raytracing_type, &
       & sunglint_type,   &
       & rttov_coef
  IMPLICIT NONE
  TYPE(profile_type   ), INTENT(IN)    :: profiles(:)
  TYPE(profile_aux    ), INTENT(IN)    :: aux
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(sunglint_type  ), INTENT(INOUT) :: sunglint
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
End Subroutine
End Interface
