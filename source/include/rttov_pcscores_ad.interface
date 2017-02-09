Interface
SUBROUTINE rttov_pcscores_ad( &
            & opts,            &
            & chanprof,        &
            & chanprof_pc,     &
            & pccomp,          &
            & pccomp_ad,       &
            & coef_pccomp,     &
            & radiancedata_ad)
  USE rttov_types, ONLY :  &
       & rttov_options,     &
       & rttov_chanprof,    &
       & rttov_pccomp,      &
       & rttov_coef_pccomp, &
       & radiance_Type
  IMPLICIT NONE
  TYPE(rttov_options    ), INTENT(IN)    :: opts
  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof   (:) ! Channel indices
  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof_pc(:) ! Channel indices
  TYPE(rttov_pccomp     ), INTENT(INOUT) :: pccomp_ad
  TYPE(rttov_pccomp     ), INTENT(INOUT) :: pccomp
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  TYPE(radiance_Type    ), INTENT(INOUT) :: radiancedata_ad
End Subroutine
End Interface
