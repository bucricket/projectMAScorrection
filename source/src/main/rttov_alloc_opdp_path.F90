!
SUBROUTINE rttov_alloc_opdp_path( &
            & ERR,       &
            & opts,      &
            & opdp_path, &
            & nlevels,   &
            & nchannels, &
            & asw,       &
            & init)
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
!INTF_OFF
#include "throw.h"
!INTF_ON
! Imported Type Definitions:
  USE rttov_types, ONLY : rttov_options, opdp_path_Type
  USE parkind1, ONLY : jpim, jplm
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(out):
  INTEGER(KIND=jpim)  , INTENT(OUT)             :: ERR      ! return code
  TYPE(rttov_options) , INTENT(IN)              :: opts
  INTEGER(KIND=jpim)  , INTENT(IN)              :: nlevels
  INTEGER(KIND=jpim)  , INTENT(IN)              :: nchannels
  TYPE(opdp_path_Type), INTENT(INOUT)           :: opdp_path
  INTEGER(KIND=jpim)  , INTENT(IN)              :: asw      ! 1=allocate, 0=deallocate
  LOGICAL(KIND=jplm)  , INTENT(IN)   , OPTIONAL :: init
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_init_opdp_path.interface"
! Local Arrays and Scalars:
  LOGICAL(KIND=jplm) :: init1
!- End of header --------------------------------------------------------
  TRY
  init1 = .FALSE.
  IF (PRESENT(init)) init1 = init

  IF (asw .EQ. 1) THEN
    NULLIFY (opdp_path%atm_level)
    NULLIFY (opdp_path%sun_level_path2)

    ALLOCATE (opdp_path%atm_level(nlevels, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of opdp_path%atm_level")
    IF (opts%rt_ir%addsolar) THEN
      ALLOCATE (opdp_path%sun_level_path2(nlevels, nchannels), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of opdp_path%sun_level_path2")
    ENDIF

    IF (init1) CALL rttov_init_opdp_path(opts, opdp_path)
  ENDIF

  IF (asw .EQ. 0) THEN
    DEALLOCATE (opdp_path%atm_level, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of opdp_path%atm_level")
    IF (opts%rt_ir%addsolar) THEN
      DEALLOCATE (opdp_path%sun_level_path2, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of opdp_path%sun_level_path2")
    ENDIF

    NULLIFY (opdp_path%atm_level)
    NULLIFY (opdp_path%sun_level_path2)
  ENDIF
  CATCH
END SUBROUTINE rttov_alloc_opdp_path
