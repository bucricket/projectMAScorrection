Interface
SUBROUTINE rttov_profaux_ad( &
            & opts,    &
            & prof,    &
            & prof_ad, &
            & coef,    &
            & aux,     &
            & aux_ad)
  USE rttov_types, ONLY :  &
       & rttov_options, &
       & rttov_coef,    &
       & profile_Type,  &
       & profile_aux
  IMPLICIT NONE
  TYPE(rttov_options), INTENT(IN)    :: opts
  TYPE(profile_Type ), INTENT(IN)    :: prof   (:)
  TYPE(profile_Type ), INTENT(INOUT) :: prof_ad(:)
  TYPE(rttov_coef   ), INTENT(IN)    :: coef
  TYPE(profile_aux  ), INTENT(IN)    :: aux
  TYPE(profile_aux  ), INTENT(INOUT) :: aux_ad
End Subroutine
End Interface
