!
Subroutine rttov789_readcoeffs  (&
       & errorstatus,   &! out
       & coef,          &! out
       & instrument,    &! in Optional
       & kmyproc,       &! in Optional
       & kioproc,       &! in Optional
       & file_id,       &! in Optional
       & channels      ) ! in Optional 
  ! Description:
  !
  ! Read an ASCII or binary coefficient file and allocate coeff structure
  !   arrays according to the optional list of channels.
  !!!!!!!
  ! This version can run in a distibuted mode :
  ! IO PE will read the data which are broadcasted to the other pes.
  ! Be careful any 'reading' modifications in rttov_readcoeffs_ascii or
  ! rttov_readcoeffs_binary may have to be reported in rttov_distribcoeffs
  !!!!!!!
  ! The optional arguments instrument and file_id determines whether the
  !   file is already opened or not
  ! if  "instrument" is present the routine will try to open
  !         the corresponding binary file  (extension .bin) in read only mode.
  !         If it fails then it tries to open the ASCII file (extension .dat)
  !         File is closed before return.
  ! if  "instrument" is not present but file_id is present the  routine will
  !  access to the coefficient file already opened with the logical unit file_id.
  !  The ASCII/binary test is performed by reading the first characters, binary
  !  files will always start by "%RTTOV_COEFF" characters. An ASCII file cannot
  !  contain such a string at the beginning of the file because it will be
  !  considered as a section name which will not be recognised. File is NOT
  !  closed on return.
  ! The user can provide an optional list of channels in "channels" argument
  !  array to reduce the output coefficient structure to this list. This
  ! can be important for reducing the memory allocation required when running
  ! with advanced IR sounders (e.g. AIRS or IASI). If the user
  !  wants all channels the "channels" argument shall not be present.
  !
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
  !  1.1       02/01/2003  A few comments added (R Saunders)
  !  1.2       03/05/2004  Add specific RTTOV8 CO2 variable (P. Brunel)
  !  1.3       02/06/2004  Change tests on id_comp_lvl == 7 by tests on fmv_model_ver (P. Brunel)
  !  1.4       08/09/2004  Change ascii/binary file test to use Inquire (J. Cameron)
  !  1.5       05/07/2007  Nullify coef_scatt_ir structure (P. Marguinaud)
  !  1.6       11/10/2007  Nullify other unused pointers (P. Marguinaud)
  !  1.7       23/94/2008  Add some initialisation of scalars (P. Brunel)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:
#include "throw.h"
  Use rttov_types, Only : &
        & rttov_coef
  Use parkind1, Only : jpim
!INTF_OFF
  Use rttov_const, Only :   &
        & version             ,&
        & release             ,&
        & minor_version       ,&
        & errorstatus_success ,&
        & errorstatus_fatal 
  USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
  Use parkind1, Only : jprb, jplm
!INTF_ON
  Implicit None


  ! subroutine arguments
  ! scalar arguments with intent(in):
  Integer(Kind=jpim), Optional, Intent(in) :: kmyproc  ! logical processor id
  Integer(Kind=jpim), Optional, Intent(in) :: kioproc  ! procs dedicated for io
  Integer(Kind=jpim), Optional, Intent (in) :: instrument(3)  ! (platform, satellite identification, instrument) number
  Integer(Kind=jpim), Optional, Intent (in) :: file_id      ! file logical unit number
  Integer(Kind=jpim), Optional, Intent (in) :: channels(:)      ! list of channels to extract


  ! scalar arguments with intent(out):
  Integer(Kind=jpim), Intent (out) :: errorstatus       ! return code
  Type( rttov_coef ), Intent (out) :: coef   ! coefficients

!INTF_END

#include "rttov789_coeffname.interface"
#include "rttov789_opencoeff.interface"
#include "rttov_errorreport.interface"
#include "rttov_nullify_coef.interface"
#include "rttov789_readcoeffs_ascii.interface"

  ! Local Scalars:
  Logical(Kind=jplm):: file_toclose
  Logical(Kind=jplm):: file_binary
  Logical(Kind=jplm):: existence
  Integer(Kind=jpim)   :: file_lu
  Integer(Kind=jpim)   :: imyproc,iioproc
  Integer(Kind=jpim)   :: err

  Character (len=256):: coeffname  ! file name for coefficient file
  Character (len=20) :: file_form
  Character (len=80) :: errMessage
REAL(KIND=JPRB) :: ZHOOK_HANDLE


  !- End of header --------------------------------------------------------
TRY
  ! 0 Initialise variables
  !---------------------------------------------
IF (LHOOK) CALL DR_HOOK('RTTOV789_READCOEFFS',0_jpim,ZHOOK_HANDLE)
  errorstatus     = errorstatus_success

  If ( .Not. Present (kmyproc) ) Then
     imyproc = 1
  Else
     imyproc = kmyproc
  Endif

  If ( .Not. Present (kioproc) ) Then
     iioproc = 1
  Else
     iioproc = kioproc
  Endif

  If (imyproc == iioproc ) then
     file_toclose = .False.
     file_binary  = .False.

     If ( .Not. Present (file_id) ) Then
        file_lu = 0
     Else
        file_lu = file_id
     Endif

     Write( errMessage, '( "RTTOV library version ",i2,1x,i1,".",i1 )' )&
               & version, release, minor_version 
     INFO(errMessage)

     ! 1 Beginning of coefficient opening sequence
     !---------------------------------------------

     ! test arguments instrument and file_id to decide whether to open
     ! the file or not.
     If ( Present (instrument ) ) Then

        ! Binary filename
        Call rttov789_coeffname ( err, instrument, coeffname, lbinary = .True._jplm  )
        THROW( err /= errorstatus_success )

        ! test existence of binary file
        Inquire( FILE = coeffname, EXIST = existence )
        If ( existence ) Then
           Write( errMessage, '( "open binary coefficient file ",a )' ) Trim(coeffname) 
           INFO(errMessage)
           ! Open binary file
           Call rttov789_opencoeff ( err, coeffname, file_lu, lbinary = .True._jplm )
           If ( err /= errorstatus_success ) Then
              ! Binary open fails, try ASCII access
              ! ASCII filename
              Call rttov789_coeffname ( err, instrument, coeffname )
              THROW( err /= errorstatus_success )
              ! Open ASCII file
              Call rttov789_opencoeff ( err, coeffname, file_lu)
              THROW( err /= errorstatus_success )
           Endif

        Else
           ! Try to open ASCII format
           ! ASCII filename
           Call rttov789_coeffname ( err, instrument, coeffname )
           THROW( err /= errorstatus_success )
           
           Write( errMessage, '( "open ASCII coefficient file ",a )' ) Trim(coeffname) 
           INFO(errMessage)
           ! Open ASCII file
           Call rttov789_opencoeff ( err, coeffname, file_lu)
           THROW( err /= errorstatus_success )

        End If
        file_toclose = .True.

     Else
        ! instrument argument missing
        If ( .Not. Present (file_id) ) Then
           ! file_id argument missing
           err = errorstatus_fatal
           Write( errMessage, '( "instrument and file_id missing arguments." )' )
           THROWM( err /= errorstatus_success, errMessage )
        Endif
     Endif

     ! Find out if the file is ascii or binary
     ! The inquire should work even if the file was opened externally
     INQUIRE(file_lu,FORM=file_form)
     IF ( file_form == 'FORMATTED' ) THEN
       file_binary = .FALSE.
     ELSEIF ( file_form == 'UNFORMATTED' ) THEN
       file_binary = .TRUE.
     ELSE
       err = errorstatus_fatal
       Write( errMessage, '(a)' ) 'Unknown file format: '//file_form
       THROWM( err /= errorstatus_success, errMessage )
       RETURN
     ENDIF

     ! End of coefficient opening sequence
     !-------------------------------------
  Endif


  ! 2 initialize coef structure
  !----------------------------
  Call rttov_nullify_coef(coef)

  If (imyproc == iioproc ) then
    ! 3 Read binary file
    !-------------------
    If( file_binary ) Then
       If( Present ( channels ) ) Then
       Else
       Endif

    ! 4 If no Binary file then read ASCII file
    !-----------------------------------------
    Else
       If( Present ( channels ) ) Then
          Call rttov789_readcoeffs_ascii  (&
                 & err,           &! out
                 & coef,          &! inout
                 & file_lu,       &! in
                 & channels = channels ) ! in Optional 
       Else
          Call rttov789_readcoeffs_ascii  (&
                 & err,           &! out
                 & coef,          &! inout
                 & file_lu       ) ! in 
       Endif

    Endif

    THROW( err /= errorstatus_success )

    If( file_toclose ) Then
       Close ( unit = file_lu )
    Endif

    Write( errMessage, '( "fast model version compatibility ",i2 )' )coef % fmv_model_ver
    INFO(errMessage)
  Endif

IF (LHOOK) CALL DR_HOOK('RTTOV789_READCOEFFS',1_jpim,ZHOOK_HANDLE)
CATCH
errorstatus = err
IF (LHOOK) CALL DR_HOOK('RTTOV789_READCOEFFS',1_jpim,ZHOOK_HANDLE)

End Subroutine rttov789_readcoeffs
