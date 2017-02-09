Interface
SUBROUTINE rttov_read_binary_scaercoef( &
            & ERR,           &
            & coef,          &
            & coef_scatt_ir, &
            & optp,          &
            & file_lu,       &
            & channels)
#include "throw.h"
  USE rttov_types, ONLY :  &
       & rttov_coef,          &
       & rttov_optpar_ir,     &
       & rttov_coef_scatt_ir
  USE parkind1, ONLY : jpim
  IMPLICIT NONE
  INTEGER(KIND=jpim)       , INTENT(IN)                :: file_lu       ! file logical unit number
  INTEGER(KIND=jpim)       , OPTIONAL     , INTENT(IN) :: channels   (:)! list of channels to extract
  TYPE(rttov_coef         ), INTENT(INOUT)             :: coef          ! coefficients
  TYPE(rttov_optpar_ir    ), INTENT(INOUT)             :: optp
  TYPE(rttov_coef_scatt_ir), INTENT(INOUT)             :: coef_scatt_ir
  INTEGER(KIND=jpim)       , INTENT(OUT)               :: ERR           ! return code
End Subroutine
End Interface
