!
SUBROUTINE rttov_calc_weighting_fn (err, p, opdep, weighting_fn)
  ! Description:
  ! Given a pressure profile and level to space optical depths
  !   this routine calculates the weighting function as
  !
  !   d(transmittance)/d(-log(p))
  !
  ! The input optical depths can be taken from
  !   traj % opdp_path % atm_level(:, ichan) for example.
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
  !    Copyright 2007, EUMETSAT, All Rights Reserved.
  !
  ! Method:
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  1.0       20/06/2011  Created (J Hocking)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:
  ! Imported Parameters:

  ! Imported Type Definitions:
  USE parkind1, ONLY : jprb, jpim

#include "throw.h"

!INTF_OFF
  USE rttov_const, ONLY :       &
        & errorstatus_success,  &
        & errorstatus_fatal
!INTF_ON

  IMPLICIT NONE

  ! subroutine arguments
  INTEGER(KIND=jpim), INTENT(OUT) :: err                      ! error status
  REAL(KIND=jprb),    INTENT(IN)  :: p(:)                     ! pressure levels
  REAL(KIND=jprb),    INTENT(IN)  :: opdep(size(p))           ! optical depth (levels to space)
  REAL(KIND=jprb),    INTENT(OUT) :: weighting_fn(size(p)-1)  ! calculated weighting function

!INTF_END

  ! local variables
  INTEGER(KIND=jpim) :: i
  INTEGER(KIND=jpim) :: nlevels

  TRY
  nlevels = size(p)
  
  DO i = 2, nlevels
    
    IF (p(i-1) >= p(i)) THEN
      err = errorstatus_fatal
      THROWM(err.ne.0, "Error in pressure levels")
    ELSE
      weighting_fn(i-1) = (exp(opdep(i)) - exp(opdep(i-1)))/(log(p(i-1)) - log(p(i)))
    END IF
  
  END DO
  
  CATCH
END SUBROUTINE rttov_calc_weighting_fn
