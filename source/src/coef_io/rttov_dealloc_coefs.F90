! Description:
!> @file
!!   Deallocate a coefficients structure.
!
!> @brief
!!   Deallocate a coefficients structure.
!!
!! @param[out]     err     status on exit
!! @param[in,out]  coefs   coefficients structure
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
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_dealloc_coefs(err, coefs)

#include "throw.h"

  USE rttov_types, ONLY : rttov_coefs
  USE parkind1, ONLY : jpim

  IMPLICIT NONE

  INTEGER(KIND=jpim), INTENT(OUT)   :: err
  TYPE(rttov_coefs),  INTENT(INOUT) :: coefs
!INTF_END

#include "rttov_errorreport.interface"
#include "rttov_dealloc_optpar_ir.interface"
#include "rttov_dealloc_coef_scatt_ir.interface"
#include "rttov_dealloc_coef_pccomp.interface"
#include "rttov_dealloc_coef.interface"
!- End of header --------------------------------------------------------

  TRY

  IF (ASSOCIATED(coefs%optp%optpaer) .OR. &
      ASSOCIATED(coefs%optp%optpwcl) .OR. &
      ASSOCIATED(coefs%optp%optpicl)) THEN
    CALL rttov_dealloc_optpar_ir(err, coefs%optp)
    THROW(err.NE.0)
  ENDIF

  IF(ASSOCIATED(coefs%coef_scatt_ir%fmv_aer_rh) .OR. &
     ASSOCIATED(coefs%coef_scatt_ir%fmv_wcl_rh) .OR. &
     ASSOCIATED(coefs%coef_scatt_ir%fmv_icl_dg)) THEN
    CALL rttov_dealloc_coef_scatt_ir (err, coefs%coef_scatt_ir)
    THROW(err.NE.0)
  ENDIF

  IF(ASSOCIATED(coefs%coef_pccomp%pcreg)) THEN
    CALL rttov_dealloc_coef_pccomp(err, coefs%coef_pccomp)
    THROW(err.NE.0)
  ENDIF

  CALL rttov_dealloc_coef(err, coefs%coef)
  THROW(err.NE.0)

  coefs%initialised = .FALSE.

  CATCH
END SUBROUTINE
