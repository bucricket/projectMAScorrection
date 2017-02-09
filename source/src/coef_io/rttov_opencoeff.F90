!
SUBROUTINE rttov_opencoeff (&
       & err,        &
       & coeffname,  &
       & file_id,    &
       & for_output, &
       & lbinary     )
  ! Description:
  ! Opens a file given by the name "coeffname"  with  logical
  !   unit file_id for output (for_output= .true.) or input and returns
  !   the error status errorstatus.
  ! If file_id input is zero the routine uses the first free logical unit.
  !   The optional logical argument lbinary determines the expected data storage.
  ! If lbinary is false or not present the file is assumed as a sequential
  !   formatted, in other case it is sequential unformatted.
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

#include "throw.h"
  USE parkind1, ONLY : jpim, jplm

  IMPLICIT NONE

  INTEGER(KIND=jpim), INTENT(OUT)          :: err
  CHARACTER(LEN=*),   INTENT(IN)           :: coeffname  ! filename of the coefficient file
  INTEGER(KIND=jpim), INTENT(INOUT)        :: file_id
  LOGICAL(KIND=jplm), INTENT(IN), OPTIONAL :: for_output ! file access mode
  LOGICAL(KIND=jplm), INTENT(IN), OPTIONAL :: lbinary    ! if binary file wanted
!INTF_END

#include "rttov_errorreport.interface"

  CHARACTER(LEN=8)   :: file_status
  CHARACTER(LEN=8)   :: file_action
  CHARACTER(LEN=16)  :: file_form
  INTEGER(KIND=jpim) :: file_output  ! 1 for output; 0 for input
  LOGICAL            :: file_open
  LOGICAL            :: existence
  INTEGER(KIND=jpim) :: file_unit

  !- End of header --------------------------------------------------------
TRY

  file_unit = file_id

  ! Consider file_id argument to determine unit for open
  ! Be careful of the following loop for searching
  ! the first free logical unit. It has been observed that
  ! with some high level compiler options it can have some
  ! side effect, like returning file_id with 0 value.
  IF (file_id <= 0) THEN

    ! Find first free logical unit
    file_unit = 9
    file_open = .TRUE.
    DO
      file_unit = file_unit + 1
      INQUIRE(file_unit, OPENED=file_open)
      IF (.NOT. file_open) EXIT
    ENDDO
  ENDIF

  IF (file_id <= 0 .AND. file_unit >= 9) THEN
    file_id = file_unit
  ENDIF

  ! Consider lbinary option to create the option
  file_form = 'formatted'
  IF (PRESENT(lbinary)) THEN
    IF (lbinary) file_form = 'unformatted'
  ENDIF

  ! Access mode
  file_output = 0
  IF (PRESENT(for_output)) THEN
    IF (for_output) file_output = 1
  ENDIF


  ! Check data file existence

  INQUIRE(FILE=coeffname, EXIST=existence)
  IF (file_output == 0) THEN
    ! If data file does not exist, return an error
    IF (.NOT. existence) err = errorstatus_fatal
    THROWM(err .NE. 0, "Coefficient file"//TRIM(coeffname)//" not found")

    ! Set OPEN keywords for reading
    file_status = 'OLD   '
    file_action = 'READ '
  ELSE
    ! If data file does exist, output a warning message
    IF (existence) THEN
      INFO("Coefficient file"//TRIM(coeffname)//" will be overwritten")
    ENDIF

    ! Set OPEN keywords for writing
    file_status = 'REPLACE'
    file_action = 'WRITE'
  ENDIF


  ! Open the data file

  OPEN (file_id, FILE = coeffname,  &
        STATUS = TRIM(file_status), &
        ACTION = TRIM(file_action), &
        ACCESS = 'SEQUENTIAL',      &
        FORM   = TRIM(file_form),   &
        IOSTAT = err)
  THROWM(err .NE. 0, "Error opening "//TRIM(coeffname))

CATCH
END SUBROUTINE rttov_opencoeff
