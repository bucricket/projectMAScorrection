!
SUBROUTINE rttov_init_coefs(  ERR, opts, coefs )
! Description:
!   Initialise the coefficients structure.
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
! Imported Parameters:
#include "throw.h"
! Imported Type Definitions:
  USE rttov_types, ONLY :  &
       & rttov_coefs,   &
       & rttov_options
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE parkind1, ONLY : jprb
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(in):
! scalar arguments with intent(out):
  INTEGER(KIND=jpim)       , INTENT(OUT)               :: ERR          ! return code
  TYPE(rttov_options),       INTENT(IN)                :: opts
  TYPE(rttov_coefs         ), INTENT(INOUT)            :: coefs        ! coefficients
!INTF_END
  REAL(KIND=jprb) :: ZHOOK_HANDLE
#include "rttov_errorreport.interface"
#include "rttov_init_coef.interface"
#include "rttov_init_coef_optpar_ir.interface"
#include "rttov_init_coef_pccomp.interface"
! Local Scalars:
!- End of header --------------------------------------------------------
  TRY
!
  IF (LHOOK) CALL DR_HOOK('RTTOV_INIT_COEFS', 0_jpim, ZHOOK_HANDLE)

  IF (.NOT. coefs%initialised) THEN

    CALL rttov_init_coef(  ERR, coefs%coef )
    THROW(err.ne.0)

    IF( opts%rt_ir%addclouds ) THEN
      CALL rttov_init_coef_optpar_ir(   &
              & ERR,                 &
              & coefs%coef,          &
              & coefs%optp)
      THROW(err.ne.0)
    ENDIF

    IF( opts%rt_ir%pc%addpc ) THEN
      CALL rttov_init_coef_pccomp(   &
              & ERR,                 &
              & coefs%coef,          &
              & coefs%coef_pccomp)
      THROW(err.ne.0)
    ENDIF

    coefs%initialised = .TRUE.

  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_INIT_COEFS', 1_jpim, ZHOOK_HANDLE)
  CATCH
  IF (LHOOK) CALL DR_HOOK('RTTOV_INIT_COEFS', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE 
