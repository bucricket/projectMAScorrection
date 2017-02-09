SUBROUTINE rttov_nullify_optpar_ir(optp)
! Description:
!   Nullify the IR optical parameter coefficient structure
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
  USE rttov_types, ONLY : rttov_optpar_ir
  IMPLICIT NONE
  TYPE(rttov_optpar_ir), INTENT(INOUT) :: optp
!INTF_END
  NULLIFY (optp%optpaer)
  NULLIFY (optp%optpwcl)
  NULLIFY (optp%optpicl)
  NULLIFY (optp%optpiclb)
END SUBROUTINE 
