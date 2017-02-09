SUBROUTINE rttov_alloc_sunglint( &
            & err,       &
            & sunglint,  &
            & nprofiles, &
            & nomega,    &
            & asw,       &
            & init,      &
            & direct)
! Description:
!   Allocates/deallocates the sunglint structure
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

!INTF_OFF
#include "throw.h"
!INTF_ON
  USE parkind1, ONLY : jpim, jplm
  USE rttov_types, ONLY : sunglint_type
  INTEGER(KIND=jpim) , INTENT(INOUT)        :: err
  TYPE(sunglint_type), INTENT(INOUT)        :: sunglint
  INTEGER(KIND=jpim) , INTENT(IN)           :: nprofiles
  INTEGER(KIND=jpim) , INTENT(IN)           :: nomega
  INTEGER(KIND=jpim) , INTENT(IN)           :: asw
  LOGICAL(KIND=jplm) , INTENT(IN), OPTIONAL :: init
  LOGICAL(KIND=jplm) , INTENT(IN), OPTIONAL :: direct
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_init_sunglint.interface"

  LOGICAL(KIND=jplm) :: init1, direct1

  TRY

  init1 = .FALSE.
  IF (PRESENT(init)) init1 = init
  direct1 = .FALSE.
  IF (PRESENT(direct)) direct1 = direct
  IF (asw .EQ. 1) THEN
    NULLIFY (sunglint%s, sunglint%beta, sunglint%psi)
    ALLOCATE (sunglint%s(nprofiles), STAT = err)
    THROWM(err.ne.0,"Allocation of sunglint failed")
    IF (direct1) THEN
      ALLOCATE (sunglint%beta(nomega, nprofiles), sunglint%psi(nomega, nprofiles), STAT = err)
      THROWM(err.ne.0,"Allocation of sunglint failed")
    ENDIF
    IF (init1) THEN
      CALL rttov_init_sunglint(sunglint)
    ENDIF
  ENDIF
  IF (asw .EQ. 0) THEN
    IF (ASSOCIATED(sunglint%s)) DEALLOCATE (sunglint%s, STAT = err)
    THROWM(err.ne.0,"Deallocation of sunglint%s failed")
    IF (ASSOCIATED(sunglint%beta)) DEALLOCATE (sunglint%beta, STAT = err)
    THROWM(err.ne.0,"Deallocation of sunglint%beta failed")
    IF (ASSOCIATED(sunglint%psi)) DEALLOCATE (sunglint%psi, STAT = err)
    THROWM(err.ne.0,"Deallocation of sunglint%psi failed")
    NULLIFY (sunglint%s, sunglint%beta, sunglint%psi)
  ENDIF

  CATCH

END SUBROUTINE 
