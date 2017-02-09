!
Subroutine rttov789_opencoeff (&
       & ERR,&
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
  ! Method:
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  1.0       01/12/2002  New F90 code with structures (P Brunel A Smith)
  !
  ! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Imported Parameters:
   
#include "throw.h"
  Use parkind1, Only : jpim, jplm


  Implicit None

  ! subroutine arguments
  ! scalar arguments with intent(in):
  Character (*), Intent (in) :: coeffname  ! filename of the coefficient file
  Logical(Kind=jplm), Optional, Intent (in) :: for_output     ! file access mode
  Logical(Kind=jplm), Optional, Intent (in) :: lbinary        ! if binary file wanted

  ! scalar arguments with intent(inout):
  Integer(Kind=jpim), Intent(inout) :: file_id

  ! scalar arguments with intent(out):
  Integer(Kind=jpim), Intent(out) :: ERR


!INTF_END

#include "rttov_errorreport.interface"



  ! Local Scalars:
  Character (len=8) :: file_status
  Character (len=8) :: file_action
  Character (len=16):: file_form
  Integer(Kind=jpim)           :: file_output      ! 1 for output; 0 for input
  Logical(Kind=jplm):: file_Open
  Logical(Kind=jplm):: existence
  Integer(Kind=jpim)           :: file_unit

  !- End of header --------------------------------------------------------
TRY

  file_unit = file_id

  ! Consider file_id argument to determine unit for open
  ! Be careful of the following loop for searching
  ! the first free logical unit. It has been observed that
  ! with some high level compiler options it can have some
  ! side effect, like returning file_id with 0 value.
  If( file_id <= 0 ) Then
     ! get the first free logical unit


     ! --------------------------------------------
     ! Initialise logical unit number and file_open
     ! --------------------------------------------

     file_unit = 9
     file_Open = .True.


     ! ------------------------------
     ! Start open loop for file_id search
     ! ------------------------------

     file_search: Do

        ! -- Increment logical unit number
        file_unit = file_unit + 1

        ! -- Check if file is open
        Inquire( file_unit, OPENED = file_Open )

        ! -- Is this file_id available?
        If ( .Not. file_Open ) Exit file_search

     End Do file_search

  Endif

  If( file_id <= 0 .and. file_unit >= 9) Then
     file_id = file_unit
  End If


  ! Consider lbinary option to create the option
  If(Present(lbinary)) Then
     If(lbinary) Then
        file_form = 'unformatted'
     Else
        file_form = 'formatted'
     Endif
  Else
     file_form = 'formatted'
  Endif

  ! mode access
  If(Present(for_output)) Then
     If(for_output) Then
        file_output = 1
     Else
        file_output = 0
     Endif
  Else
     file_output = 0
  Endif

  !#--------------------------------------------------------------------------#
  !#                      -- Check data file existence --                     #
  !#--------------------------------------------------------------------------#

  Inquire( FILE = coeffname, EXIST = existence )
  If ( file_output == 0 ) Then
     ! -- If data file does not exist, return an error

     If ( .Not. existence ) err = errorstatus_fatal
     THROWM( ERR .NE. 0, "Coefficient file"//TRIM(coeffname)//" not found" )

     ! -- Set OPEN keywords for reading
     file_status = 'OLD   '
     file_action = 'READ '

  Else

     ! -- If data file does exist, output a warning message
     If ( existence ) Then
        INFO( "Coefficient file"//TRIM(coeffname)//" will be overwritten" )
     End If

     ! -- Set OPEN keywords for writing
     file_status = 'REPLACE'
     file_action = 'WRITE'

  End If

  !#--------------------------------------------------------------------------#
  !#                        -- Open the data file --                          #
  !#--------------------------------------------------------------------------#

  Open( file_id, FILE   = coeffname, &
        & STATUS = Trim( file_status ), &
        & ACTION = Trim( file_action ), &
        & ACCESS = 'SEQUENTIAL', &
        & FORM   = Trim( file_form   ), &
        & IOSTAT = ERR ) 
  THROWM( ERR .NE. 0, "Error opening "//TRIM(coeffname) )


CATCH
End Subroutine rttov789_opencoeff
