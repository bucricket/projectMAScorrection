!
Subroutine rttov_errorreport (error_status, error_message, name_of_routine)
  ! Description:
  !   Write out fatal error and informational messages to the unit
  !   specified by rttov_errorhandling.
  !
  ! Copyright:
  !
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
  !
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version  Date      Comment
  !
  !  1.0     01/04/11  Created, based on routine rttov_errorreport. (J Hocking)
  !
  ! Code Description:
  !   FORTRAN 90
  !
  ! Declarations
  !
  ! Modules used:

  USE parkind1, ONLY : jpim

!INTF_OFF
  USE rttov_const, ONLY :  &
      errorstatus_success
      
  USE rttov_global, ONLY : &
      err_init,            &
      error_unit
      
!INTF_ON

  IMPLICIT NONE

  ! Subroutine arguments
  
  INTEGER(KIND=jpim),           INTENT(IN) :: error_status    ! Success or fatal
  CHARACTER(LEN=*),             INTENT(IN) :: error_message   ! Message to output
  CHARACTER(LEN=*),   OPTIONAL, INTENT(IN) :: name_of_routine ! Routine calling this one
  
!INTF_END

#include "rttov_errorhandling.interface"

  ! local variables
  CHARACTER(LEN=8)   :: date
  CHARACTER(LEN=10)  :: time
  CHARACTER(LEN=21)  :: datetime
  
  !- End of header --------------------------------------------------------

  CALL DATE_AND_TIME(date, time)

  WRITE(datetime,"(1X,a4,2('/',a2),2x,2(a2,':'),a2)")   &
              & date(1:4), date(5:6), date(7:8),        &
              & time(1:2), time(3:4), time(5:6)
  
  ! If globlal variables not defined then use default values
  IF( .NOT. err_init ) CALL rttov_errorhandling(-1_jpim)

  IF ( error_status == errorstatus_success ) THEN

    ! Standard warning or informational message
    IF (PRESENT(name_of_routine)) THEN
      WRITE(error_unit,"(3a)") TRIM(datetime), "  ", TRIM(name_of_routine)
      WRITE(error_unit,"(5X,A)") TRIM(error_message)
    ELSE
      WRITE(error_unit,"(3a)") TRIM(datetime), "  ", TRIM(error_message)
    END IF
    
  ELSE
    
    ! fatal error message
    IF (PRESENT(name_of_routine)) THEN
      WRITE(error_unit,"(3a)") TRIM(datetime), "  fatal error in module ",TRIM(name_of_routine)
    ELSE
      WRITE(error_unit,"(2a)") TRIM(datetime), "  fatal error"
    END IF
    WRITE(error_unit,"(5X,A)") TRIM(error_message)
    
  END IF

End Subroutine rttov_errorreport
