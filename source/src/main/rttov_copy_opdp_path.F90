!
SUBROUTINE rttov_copy_opdp_path(opts, opdp_pathA, opdp_pathB)
! Description:
! allocation/deallocation of a opdp_path structure
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
!  1.0       01/12/2002  New F90 code
!  1.1       11/10/2007  Add addclouds, addaerosl, init logicals
!                        nullify unused pointers P.Marguinaud
!  1.2       03/11/2009  Transmittances / optical depths on levels (A Geer)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
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
  USE rttov_types, ONLY : rttov_options, opdp_path_Type
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(out):
  TYPE(rttov_options),  INTENT(IN)    :: opts
  TYPE(opdp_path_Type), INTENT(INOUT) :: opdp_pathA
  TYPE(opdp_path_Type), INTENT(IN)    :: opdp_pathB
!INTF_END
!- End of header --------------------------------------------------------
  opdp_pathA%atm_level = opdp_pathB%atm_level
  IF (opts%rt_ir%addsolar) THEN
    opdp_pathA%sun_level_path2 = opdp_pathB%sun_level_path2
  ENDIF
END SUBROUTINE 
