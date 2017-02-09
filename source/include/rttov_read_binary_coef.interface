Interface
SUBROUTINE rttov_read_binary_coef( &
            & ERR,           &
            & coef,          &
            & file_lu,       &
            & channels)
#include "throw.h"
  USE rttov_types, ONLY :  &
       & rttov_coef
  USE parkind1, ONLY : jpim
  IMPLICIT NONE
  INTEGER(KIND=jpim)       , INTENT(IN)                :: file_lu       ! file logical unit number
  INTEGER(KIND=jpim)       , OPTIONAL     , INTENT(IN) :: channels   (:)! list of channels to extract
  TYPE(rttov_coef         ), INTENT(INOUT)             :: coef          ! coefficients
  INTEGER(KIND=jpim)       , INTENT(OUT)               :: ERR           ! return code
End Subroutine
End Interface
