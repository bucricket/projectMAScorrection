MODULE rttov_coef_io_mod
! Description:
!   Useful functions for I/O
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
USE parkind1, ONLY : jpim, jprb, jplm

#include "throw.h"

IMPLICIT NONE

#include "rttov_errorreport.interface"

CONTAINS

SUBROUTINE getlun(err, lun1, f1, form1, loutput, filetype, lun, f, form, instrument)
  ! Open a coefficient file for I/O: if lun was passed as input, use it, otherwise open the file
  ! whose name was passed in f, otherwise create the coefficient file name from the instrument
  ! definition

  USE rttov_const, ONLY :  &
       & rttov_magic_string

  INTEGER(KIND=jpim), INTENT(OUT)          :: err
  INTEGER(KIND=jpim), INTENT(OUT)          :: lun1     ! logical unit the file was opened with
  CHARACTER(LEN=*),   INTENT(OUT)          :: f1       ! filename of opened file
  CHARACTER(LEN=*),   INTENT(OUT)          :: form1    ! format (unformatted/formatted/hdf5) the file was opened with
  LOGICAL(KIND=jplm), INTENT(IN)           :: loutput  ! FALSE => open for input, TRUE => open for output
  CHARACTER(LEN=*),   INTENT(IN)           :: filetype ! file type = 'rtcoef', 'pccoef', ...
  INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: lun
  CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: f
  CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: form
  INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: instrument(3)

#include "rttov_cmpuc.interface"
#include "rttov_coeffname.interface"

  CHARACTER(LEN=8)   :: hdf5_format_signature
  ! DAR - removed kind=jplm from line below to comply with 2003 std which states
  !       that a logical must be of default type for an intrinsic (paraphrasing)
  LOGICAL            :: lopen, lexist 
  CHARACTER(LEN=4)   :: int4ch
  CHARACTER(LEN=8)   :: int8ch
  INTEGER(KIND=jpim) :: recl
  CHARACTER(LEN=16)  :: bin_check_string
  CHARACTER(LEN=8)   :: hdf5_check_signature
  CHARACTER(LEN=128) :: filestem

TRY

  hdf5_format_signature = CHAR(137)//'HDF'//CHAR(13)//CHAR(10)//CHAR(26)//CHAR(10)

  form1 = ''
  IF (PRESENT(form)) form1 = form

  ! If lun was passed as argument, then take it and return
  IF (PRESENT(lun)) THEN
    lun1 = lun

    ! If form was passed as argument then take it, otherwise make a guess
    IF (PRESENT(form)) THEN
      IF (form /= '') RETURN ! form1 = form (see above)
    ENDIF
    INQUIRE(lun1, form = form1)
    RETURN
  ENDIF

  ! Try to guess the filename
  IF (PRESENT(f)) THEN

    ! If f was passed as argument, then take it
    f1 = f
    IF (loutput .AND. (form1 == '')) THEN
      err = errorstatus_fatal
      THROWM(err.ne.0,'Format argument missing for output file')
    ENDIF

  ELSE IF (PRESENT(instrument)) THEN
    ! Guess the filename from the instrument triplet

    ! Get the filename from instrument triplet (excluding file extension)
    CALL rttov_coeffname(err, instrument, filetype, filestem)
    THROW(err.ne.0)

    IF (loutput .OR. form1 /= '') THEN

      ! Formatted output is the default
      IF (rttov_cmpuc('unformatted', form1)) THEN
        f1 = TRIM(filestem)//'.bin'
      ELSE IF (rttov_cmpuc('hdf5', form1)) THEN
        f1 = TRIM(filestem)//'.H5'
      ELSE
        form1 = 'formatted'
        f1 = TRIM(filestem)//'.dat'
      ENDIF

    ELSE

      ! Guess the input format (check for file existence)
      form1 = 'unformatted'
      f1 = TRIM(filestem)//'.bin'
      INQUIRE(file = f1, exist = lexist)

      IF (.NOT. lexist) THEN
        form1 = 'hdf5'
        f1 = TRIM(filestem)//'.h5'
        INQUIRE(file = f1, exist = lexist)
      ENDIF

      IF (.NOT. lexist) THEN
        form1 = 'hdf5'
        f1 = TRIM(filestem)//'.H5'
        INQUIRE(file = f1, exist = lexist)
      ENDIF

      IF (.NOT. lexist) THEN
        form1 = 'formatted'
        f1 = TRIM(filestem)//'.dat'
        INQUIRE(file = f1, exist = lexist)
      ENDIF

      IF (.NOT. lexist) THEN
        err = errorstatus_fatal
        THROWM(err.ne.0,'Cannot find any coefficient file named: '//TRIM(filestem)//'.dat/.bin/.h5/.H5')
      ENDIF

    ENDIF
  ELSE
    err = errorstatus_fatal
    THROWM(err.ne.0,'Opening the coefficient file requires the filename or the instrument ID')
  ENDIF

  ! Look for a free logical unit (not thread-safe)
  DO lun1 = 9, 99
    INQUIRE(lun1, opened = lopen)
    IF (.NOT. lopen) EXIT
  ENDDO
  IF (lopen) THEN
    err = errorstatus_fatal
    THROWM(err.ne.0,'Cannot find a free lun')
  ENDIF

  ! HDF5 files are not opened by this routine: just return file name,
  !   open will be done by load/save routines from rttov_hdf modules
  IF (loutput) THEN

    IF ( .NOT. rttov_cmpuc('hdf5', form1)) THEN
      OPEN(lun1, file = f1, form = form1, status = 'replace', action = 'write', iostat = err)
      THROWM(err.ne.0,'Cannot open '//TRIM(f1))
    ENDIF

  ELSE

    IF (form1 == '') THEN
      ! Guess the input file format with the embedded magic string

      ! Do not just open as unformatted and try reading the magic string because ifort takes the
      ! first 4/8 bytes as the record length and may allocate a huge amount of memory for the read
      ! if the file is actually formatted.
      form1 = 'formatted'

      INQUIRE(iolength = recl) int4ch, bin_check_string
      OPEN(lun1, file = f1, form = 'unformatted', access = 'direct',recl = recl, &
                  status = 'old', action = 'read', iostat = err)
      THROWM(err.ne.0,'Cannot open '//TRIM(f1))
      READ(lun1, rec = 1, iostat = err) int4ch, bin_check_string
      CLOSE(lun1, iostat = err)
      THROW(err.ne.0)
      IF (bin_check_string == rttov_magic_string) form1 = 'unformatted'

      INQUIRE(iolength = recl) int8ch, bin_check_string
      OPEN(lun1, file = f1, form = 'unformatted', access = 'direct',recl = recl, &
                  status = 'old', action = 'read', iostat = err)
      THROWM(err.ne.0,'Cannot open '//TRIM(f1))
      READ(lun1, rec = 1, iostat = err) int8ch, bin_check_string
      CLOSE(lun1, iostat = err)
      THROW(err.ne.0)
      IF (bin_check_string == rttov_magic_string) form1 = 'unformatted'

      INQUIRE(iolength = recl) hdf5_check_signature
      OPEN(lun1, file = f1, form = 'unformatted', access = 'direct',recl = recl, &
                  status = 'old', action = 'read', iostat = err)
      THROWM(err.ne.0,'Cannot open '//TRIM(f1))
      READ(lun1, rec = 1, iostat = err) hdf5_check_signature
      CLOSE(lun1, iostat = err)
      THROW(err.ne.0)
      IF (hdf5_check_signature == hdf5_format_signature) form1 = 'hdf5'

    ENDIF

    IF (.NOT. rttov_cmpuc('hdf5', form1)) THEN
      OPEN(lun1, file = f1, form = form1, status = 'old', action = 'read', iostat = err)
      THROWM(err.ne.0,'Cannot open '//TRIM(f1))
    ENDIF

  ENDIF

CATCH

END SUBROUTINE getlun

SUBROUTINE closelun(err, lun1, lun)
  ! Only close lun1 if lun is not present (so we can avoid
  ! avoid closing a logical unit that was opened by a user)
  INTEGER(KIND=jpim), INTENT(OUT)          :: err
  INTEGER(KIND=jpim), INTENT(IN)           :: lun1
  INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: lun

TRY
  IF (.NOT. PRESENT(lun)) THEN
    CLOSE(lun1, iostat = err)
    THROW(err.ne.0)
  ENDIF
CATCH
END SUBROUTINE closelun

END MODULE
