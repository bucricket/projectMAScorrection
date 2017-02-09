Interface
SUBROUTINE rttov_alloc_pc_dimensions( &
            & err,          &
            & opts,         &
            & npcscores,    &
            & nprofiles,    &
            & chanprof_in,  &
            & chanprof_pc,  &
            & asw,          &
            & channels_rec)
  USE parkind1, ONLY : jpim
  USE rttov_types, ONLY : rttov_chanprof, rttov_options
  IMPLICIT NONE
  INTEGER(KIND=jpim)  , INTENT(OUT)           :: err 
  TYPE(rttov_options ), INTENT(IN)            :: opts
  INTEGER(KIND=jpim)  , INTENT(IN)            :: npcscores
  INTEGER(KIND=jpim)  , INTENT(IN)            :: nprofiles
  TYPE(rttov_chanprof), POINTER               :: chanprof_pc(:)
  TYPE(rttov_chanprof), POINTER               :: chanprof_in(:)
  INTEGER(KIND=jpim)  , INTENT(IN)            :: asw
  INTEGER(KIND=jpim)  , INTENT(IN) , OPTIONAL :: channels_rec(:)
End Subroutine
End Interface
