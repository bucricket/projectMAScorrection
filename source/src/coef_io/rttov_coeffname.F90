!
SUBROUTINE rttov_coeffname (err, instrument, filetype, coeffname)
  ! Description:
  !
  ! Returns the file name of a coefficent file for the instrument given
  ! in argument excluding the file extension.
  ! Instrument refers to an array of 3 integers defining the satellite platform,
  ! satellite number and instrument number.
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
  !    Copyright 2002, EUMETSAT, All Rights Reserved.
  !
  ! Method:
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  1.0       01/12/2002  New F90 code with structures (P Brunel A Smith)
  !  2.0       02/12/2009  Introduced principal component capability. Marco Matricardi. ECMWF
  !
  ! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:
  ! Imported Parameters:

#include "throw.h"

  USE parkind1, ONLY : jpim

!INTF_OFF
  USE rttov_const, ONLY : &
      nplatforms,         &
      ninst,              &
      inst_name,          &
      platform_name
!INTF_ON

  IMPLICIT NONE

  INTEGER(KIND=jpim), INTENT(OUT) :: err            ! return code
  INTEGER(KIND=jpim), INTENT(IN)  :: instrument(3)  ! (platform, sat_id, inst) numbers
  CHARACTER(LEN=*),   INTENT(IN)  :: filetype       ! file type e.g. "rtcoef", "pccoef", etc
  CHARACTER(LEN=*),   INTENT(OUT) :: coeffname      ! filename of the coefficient file (excluding extension)
!INTF_END

#include "rttov_errorreport.interface"

  INTEGER(KIND=jpim) :: platform
  INTEGER(KIND=jpim) :: sat_id
  INTEGER(KIND=jpim) :: inst
  CHARACTER(LEN=2)   :: ch_sat_id
  !- End of header --------------------------------------------------------
TRY
  coeffname = 'no_name'
  err = errorstatus_success

  ! Expand instrument triplet
  platform = instrument(1)
  sat_id   = instrument(2)
  inst     = instrument(3)

  ! Test sat_id and convert to string
  IF (sat_id < 10 .AND. sat_id >= 0) THEN
    ! one digit
    WRITE(ch_sat_id,'(i1)') sat_id
  ELSE IF (sat_id >= 10 .AND. sat_id < 99) THEN
    ! two digits
    WRITE(ch_sat_id,'(i2)') sat_id
  ELSE
    ! ERROR and exit
    err = sat_id
    THROWM(err .NE. 0,"invalid sat_id")
  ENDIF

  ! Test platform number
  IF (platform <= 0 .OR. platform > nplatforms) err = errorstatus_fatal
  THROWM(err .NE. 0,"invalid platform number")

  ! Test instrument number  (0 is HIRS)
  IF (inst < 0 .OR. inst > ninst) err = errorstatus_fatal
  THROWM(err .NE. 0,"invalid instrument number")

  ! Create the file name e.g. "rtcoef_platform_satellite_inst"
  coeffname = TRIM(filetype) // '_'      // &
         & TRIM(platform_name(platform)) // &
         & '_'                           // &
         & TRIM(ch_sat_id)               // &
         & '_'                           // &
         & TRIM(inst_name(inst))

CATCH
END SUBROUTINE rttov_coeffname
