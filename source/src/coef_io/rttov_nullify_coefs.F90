SUBROUTINE rttov_nullify_coefs(coefs)
! Description:
!   Nullify the coefficients structure
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
!    Copyright 2010, EUMETSAT, All Rights Reserved.
!
! Method:
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
  USE rttov_types, ONLY : rttov_coefs
  IMPLICIT NONE
  TYPE(rttov_coefs), INTENT(INOUT) :: coefs
!INTF_END

#include "rttov_nullify_coef.interface"
#include "rttov_nullify_coef_pccomp.interface"
#include "rttov_nullify_coef_scatt_ir.interface"
#include "rttov_nullify_optpar_ir.interface"


CALL rttov_nullify_coef(coefs%coef)
CALL rttov_nullify_coef_pccomp(coefs%coef_pccomp)
CALL rttov_nullify_coef_scatt_ir(coefs%coef_scatt_ir)
CALL rttov_nullify_optpar_ir(coefs%optp)

END SUBROUTINE 
