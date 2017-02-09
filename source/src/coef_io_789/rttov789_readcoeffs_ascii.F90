!
Subroutine rttov789_readcoeffs_ascii  (&
       & errorstatus,   &! out
       & coef,          &! inout
       & file_lu,       &! in
       & channels      ) ! in Optional 
  ! Description:
  !
  ! Read an ASCII coefficient file and fills coeff structure
  !   arrays according to the optional list of channels.
  !
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
  !  1.2       24/01/2003  Add return when section END encountered (P Brunel)
  !                        any I/O error is coded as fatal
  !                        Add GAZ_UNITS section
  !  1.3       02/06/2004  New format for FMV section with RTTOV8 (P. Brunel)
  !  1.4       15/06/2004  Corrected array dimension for coef % fmv_gas_pos (R Saunders)
  !  1.5        June 2005  Added Additional arrays for RTTOV8M Marco Matricardi
  !  1.6       05/12/2005  Corrected some array types for optical  const (R Saunders)
  !  1.7       19/03/2007  Reduced no of continuation lines to below 39 (R Saunders)
  !  1.8       12/12/2007  Add option of using top level (R Saunders)  
  !  1.9       01/11/2007  Removed hardcoded section length (A Geer)
  !  1.10      26/06/2008  Introduced the case where no channels are available for the 
  !                        phase function in the solar range (M. Matricardi)
  !  1.11      08/06/2009  Made interim fix to allocate cloud/aerosol arrays with right shape (R Saunders)
  !                        Ideally the channel order in all IR files needs to be in frequency not wavelen
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:
  ! Imported Parameters:
  Use rttov_types, Only : &
        & rttov_coef
  Use parkind1, Only : jpim
!INTF_OFF
  Use rttov_const, Only :   &
        & sensor_id_hi        , sensor_id_mw        , sensor_id_ir   ,sensor_id_po   ,&
        & errorstatus_success ,&
        & errorstatus_fatal   ,&
        & gas_id_mixed        ,&
        & gas_id_watervapour  ,&
        & gas_id_ozone        ,&
        & gas_id_wvcont       ,&
        & gas_id_co2          ,&
        & gas_id_n2o          ,&
        & gas_id_co           ,&
        & gas_id_ch4          ,&
        & sensor_name         , &
        & ngases_max          ,&
        & gas_name            , &
        & gas_unit_specconc  ,&
        & lensection
  Use parkind1, Only : jprb, jplm
  USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!INTF_ON
  Implicit None


  ! subroutine arguments
  ! scalar arguments with intent(in):
  Integer(Kind=jpim), Intent (in)           :: file_lu      ! file logical unit number
  Integer(Kind=jpim), Optional, Intent (in) :: channels(:)  ! list of channels to extract

  ! scalar arguments with intent(inout):
  Type( rttov_coef ), Intent (inout) :: coef   ! coefficients

  ! scalar arguments with intent(out):
  Integer(Kind=jpim), Intent (out) :: errorstatus       ! return code

!INTF_END

#include "rttov789_opencoeff.interface"
#include "rttov_errorreport.interface"
#include "rttov_skipcommentline.interface"
#include "rttov_deletecomment.interface"
#include "rttov_cmpuc.interface"
#include "rttov_findnextsection.interface"


  ! Local Scalars:
  Integer(Kind=jpim)           :: file_channels
  Logical(Kind=jplm)           :: for_output
  Integer(Kind=jpim)           :: file_lu_coef
  Logical(Kind=jplm)           :: all_channels
  Integer(Kind=jpim):: io_status
  Integer(Kind=jpim)           :: alloc_status(10)
  Real(Kind=jprb) :: pres
  Integer(Kind=jpim) :: i,j,k,l,n,inctop
  Integer(Kind=jpim) :: gas_id
  Integer(Kind=jpim) :: nbint
  Real(Kind=jprb) :: fc_speedl

  ! Placeholders for reading fastem coefficients
  Integer(Kind=jpim) :: fastem_coef_nb
  Real(Kind=jprb)    :: fastem_coef
  
  ! pointers for generic inputs
  Integer(Kind=jpim)          :: nvalues
  Complex (Kind=jprb),    Pointer :: values_c_0(:)
  Complex (Kind=jprb),    Pointer :: values_c_1(:)
  Real(Kind=jprb),    Pointer :: values0(:)
  Real(Kind=jprb),    Pointer :: values1(:)
  Real(Kind=jprb),    Pointer :: values2(:)
  Real(Kind=jprb),    Pointer :: values3(:)
  Real(Kind=jprb),    Pointer :: values4(:)
  Integer(Kind=jpim), Pointer :: ivalues0(:)
  Integer(Kind=jpim), Pointer :: ivalues1(:)
  Real(Kind=jprb),    Pointer :: coeffsarray(:,:,:)

  Character(len=36)  :: input_string
  Character(len=32)  :: gas_Type
  Character(len=lensection)  :: section
  Character (len=80) :: errMessage
  Character (len=25) :: NameOfRoutine = 'rttov789_readcoeffs_ascii'

  Logical(Kind=jplm) :: found
REAL(KIND=JPRB) :: ZHOOK_HANDLE

  !- End of header --------------------------------------------------------


  ! 0 Initialise variables
  !---------------------------------------------
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',0_jpim,ZHOOK_HANDLE)
  errorstatus     = errorstatus_success
  alloc_status(:) = 0

  ! test presence of channels argument
  If( Present ( channels ) ) Then
     all_channels = .False.
  Else
     all_channels = .True.
  Endif


  !read the file
  readfile: Do
     Call rttov_findnextsection( file_lu, io_status, section )
     If ( io_status < 0 ) Exit !end-of-file

     ! error message if any problem when reading
     errMessage = 'io status while reading section '//section

     Select Case( Trim(section) )


     Case( 'IDENTIFICATION' )
        ! Identification section
        ! 6 lines
        Call rttov_skipcommentline ( file_lu, io_status )
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

        Read( file_lu, *, iostat = io_status )&
               & coef % id_platform,&
               & coef % id_sat,&
               & coef % id_inst 
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

        Read( file_lu, '(a)', iostat = io_status ) coef % id_Common_name
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

        Read( file_lu, *, iostat = io_status ) input_string
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif
        Select Case (input_string)
        Case (sensor_name(sensor_id_ir))
           coef % id_sensor = sensor_id_ir   ! Infrared
        Case (sensor_name(sensor_id_mw))
           coef % id_sensor = sensor_id_mw   ! Micro Wave
        Case (sensor_name(sensor_id_hi))
           coef % id_sensor = sensor_id_hi   ! High resolution
        Case (sensor_name(sensor_id_po))
           coef % id_sensor = sensor_id_po   ! Polarimetric
        Case default
           coef % id_sensor = sensor_id_ir
        End Select

        Read( file_lu, *, iostat = io_status ) coef % id_comp_lvl
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif
        ! Error if the compatibility version of the coefficient file
        ! is not in the range defined by the constant module
        If(    coef % id_comp_lvl < 7  .Or.&
               & coef % id_comp_lvl > 9 ) Then 
           errorstatus = errorstatus_fatal
           Write( errMessage,&
                  & '( "Version of coefficient file is incompatible with RTTOV library")' ) 
           Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

        Read( file_lu,'(a)', iostat = io_status ) coef % id_creation
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

        Read( file_lu, *, iostat = io_status ) coef % id_creation_date
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

     Case ('LINE-BY-LINE')
       Call rttov_skipcommentline ( file_lu, io_status )
       LoopLBL : Do i = 1, Size( coef%line_by_line )
          Read( file_lu, '(a)', iostat = io_status ) coef%line_by_line(i)
          Do j = 1, Len( coef%line_by_line(i) )
            Select Case( coef%line_by_line(i)(j:j) )
              Case( '!' )
                Exit LoopLBL
              Case( ' ' )
                Cycle
              Case Default
                Exit
            End Select
          EndDo
        EndDo LoopLBL
        If( i .le. Size( coef%line_by_line ) ) &
          & coef%line_by_line(i) = 'xxxx'

     Case ('FAST_MODEL_VARIABLES')
        Call rttov_skipcommentline ( file_lu, io_status )
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
           IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

        ! fast model variables definition
        Read( file_lu, *, iostat = io_status ) coef % fmv_model_def
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
           IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

        If( coef % id_comp_lvl > 7 ) then
          ! fast model variables version
          Read( file_lu, *, iostat = io_status ) coef % fmv_model_ver
          If(io_status /= 0) Then
            Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
            errorstatus = errorstatus_fatal
            IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
            Return
          Endif
        Else
          coef % fmv_model_ver = 7
        Endif

        ! number of channels stored
        Read( file_lu, *, iostat = io_status ) coef % fmv_chn
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
           IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

        ! Take care of the user list of channels
        ! file_channels store the number of channels in the file
        ! coef % fmv_chn is the number of channels that the user requests
        file_channels = coef % fmv_chn
        If( .Not. all_channels ) Then
           coef % fmv_chn = Size( channels )
        Endif

        ! number of gases in file
        Read( file_lu, *, iostat = io_status ) coef % fmv_gas
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
           IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif
    
        ! Flag for treating top layer optimally (only for RTTOV-9 and later files)
        coef % IncZeeman =.False.
        If( coef % id_comp_lvl > 8 ) then
          ! Include top level properly flag 
          Read( file_lu, *, iostat = io_status ) inctop
          If(io_status /= 0) Then
             Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
             errorstatus = errorstatus_fatal
             IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
             Return
          Endif
        Endif

        ! allocate arrays of FAST_MODEL_VARIABLES section
        Allocate ( coef % fmv_gas_id ( coef % fmv_gas ), stat=alloc_status(1))
        Allocate ( coef % fmv_gas_pos( ngases_max ), stat=alloc_status(2))
        Allocate ( coef % fmv_var ( coef % fmv_gas ), stat=alloc_status(3))
        Allocate ( coef % fmv_lvl ( coef % fmv_gas ), stat=alloc_status(4))
        Allocate ( coef % fmv_coe ( coef % fmv_gas ), stat=alloc_status(5))
        If( Any(alloc_status /= 0) ) Then
           errorstatus = errorstatus_fatal
           Write( errMessage, '( "allocation of fmv coefs arrays")' )
           Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        End If

        If( Any(alloc_status /= 0) ) Then
           errorstatus = errorstatus_fatal
           Write( errMessage, '( "allocation of fmv coefs arrays")' )
           Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        End If

        coef % fmv_gas_id(:) = 0_jpim
        coef % fmv_gas_pos(:) = 0_jpim
        coef % fmv_var(:) = 0_jpim
        coef % fmv_lvl(:) = 0_jpim
        coef % fmv_coe(:) = 0_jpim

        Do n = 1, coef % fmv_gas
           ! gas id. number i gas_id list (fmv_gas)
           Read(file_lu,'(a)',iostat=io_status) gas_Type
           If(io_status /= 0) Then
              Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
              errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
              Return
           Endif
           Call Rttov_deletecomment(gas_Type)

           found = .False.
           Do i = 1, ngases_max
              If ( rttov_cmpuc( gas_Type , gas_name(i) ) ) Then
                 coef % fmv_gas_id(n) = i
                 found = .True.
                 Exit
              End If
           End Do
           If( .not. found ) Then
              errorstatus = errorstatus_fatal
              Write( errMessage, &
              & '( "Error: gas type ",a," not recognised")' )Trim(gas_Type)
              Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
              Return
           Endif

           ! store also the indice of this gas in the
           ! identification list
           ! so fmv_gas_pos(1) will give position of MxG in the file
           coef % fmv_gas_pos(coef % fmv_gas_id(n)) = n

           ! number of variables/predictors by gaz
           ! number of levels(pres/absorber) by gaz (fmv_gas

           if( coef % id_comp_lvl == 9 ) then
             Read(file_lu,* ,iostat=io_status)&
               & coef % fmv_var(n), coef % fmv_coe(n), coef % fmv_lvl(n)  
           else
             Read(file_lu,* ,iostat=io_status)&
               & coef % fmv_var(n), coef % fmv_lvl(n) 
           endif

           If(io_status /= 0) Then
              Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
              errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
              Return
           Endif
           ! Transfer information to some "classical" variables
           ! with more common names
           ! Note that the number of levels is taken from the Mixed Gases line
           if( coef % id_comp_lvl == 9 ) then
             Select Case( coef % fmv_gas_id(n) )
             Case( gas_id_mixed )
                coef % nmixed        = coef % fmv_var(n)
                coef % ncmixed       = coef % fmv_coe(n)
                coef % nlevels       = coef % fmv_lvl(n)
             Case( gas_id_watervapour )
                coef % nwater        = coef % fmv_var(n)
                coef % ncwater       = coef % fmv_coe(n)
             Case( gas_id_ozone )
                coef % nozone        = coef % fmv_var(n)
                coef % ncozone       = coef % fmv_coe(n)
             Case( gas_id_wvcont )
                coef % nwvcont       = coef % fmv_var(n)
                coef % ncwvcont      = coef % fmv_coe(n)
             Case( gas_id_co2 )
                coef % nco2          = coef % fmv_var(n)
                coef % ncco2         = coef % fmv_coe(n)
             Case( gas_id_n2o )
                coef % nn2o          = coef % fmv_var(n)
                coef % ncn2o         = coef % fmv_coe(n)
             Case( gas_id_co )
                coef % nco           = coef % fmv_var(n)
                coef % ncco          = coef % fmv_coe(n)
             Case( gas_id_ch4 )
                coef % nch4          = coef % fmv_var(n)
                coef % ncch4         = coef % fmv_coe(n)
             End Select
           else
             Select Case( coef % fmv_gas_id(n) )
             Case( gas_id_mixed )
                coef % nmixed        = coef % fmv_var(n)
                coef % ncmixed       = coef % nmixed
                coef % nlevels       = coef % fmv_lvl(n)
             Case( gas_id_watervapour )
                coef % nwater        = coef % fmv_var(n)
                coef % ncwater       = coef % nwater
             Case( gas_id_ozone )
                coef % nozone        = coef % fmv_var(n)
                coef % ncozone       = coef % nozone 
             Case( gas_id_wvcont )
                coef % nwvcont       = coef % fmv_var(n)
                coef % ncwvcont      = coef % nwvcont
             Case( gas_id_co2 )
                coef % nco2          = coef % fmv_var(n)
                coef % ncco2         = coef % nco2
             Case( gas_id_n2o )
                coef % nn2o          = coef % fmv_var(n)
                coef % ncn2o         = coef % nn2o
             Case( gas_id_co )
                coef % nco           = coef % fmv_var(n)
                coef % ncco          = coef % nco 
             Case( gas_id_ch4 )
                coef % nch4          = coef % fmv_var(n)
                coef % ncch4         = coef % nch4
             End Select
             coef % fmv_coe(n) = coef % fmv_var(n)
           endif
        End Do

        !
        ! compute number of layers
        ! for optical depths
        !
        coef % nlayers = coef % nlevels  

        ! Initialise the gaz units arary with defaults values
        ! (specific concetration kg/kg)
        Allocate ( coef % gaz_units  ( coef % fmv_gas ), stat=alloc_status(1) )
        If( Any(alloc_status /= 0) ) Then
           Write( errMessage, '( "allocation of gaz units coefs arrays")' )
           Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        End If
        coef % gaz_units( : ) = gas_unit_specconc


     Case ('GAS_SPECTRAL_INTERVAL')
        Call rttov_skipcommentline ( file_lu, io_status )
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif


        Do n = 1, coef % fmv_gas
           ! gas id. number i gas_id list (fmv_gas)
           Read(file_lu,'(a)',iostat=io_status) gas_Type
           If(io_status /= 0) Then
              Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
              errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
              Return
           Endif
           Call Rttov_deletecomment(gas_Type)

           gas_id = 0
           Do i = 1, ngases_max
              If ( rttov_cmpuc( gas_Type , gas_name(i) ) ) Then
                 gas_id = i
                 Exit
              End If
           End Do
           If ( gas_id == 0 ) Write(*,'(a)') &
                 & 'Error: gas type ' // Trim(gas_Type) // ' not recognised' 


           Read(file_lu,* ,iostat=io_status) nbint 
           If(io_status /= 0) Then
              Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
              errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
              Return
           Endif

           ! Transfer information to some "classical" variables
           ! with more common names
           Select Case( gas_id )

           Case( gas_id_mixed )
           coef%nintmixed = nbint
           allocate(coef%mixedgasint(2,coef % nintmixed),stat=alloc_status(1))
           If( Any(alloc_status /= 0) ) Then
             errorstatus = errorstatus_fatal
             Write( errMessage, '( "allocation of fmv coefs arrays")' )
             Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
             Return
           End If
             Do i=1,coef%nintmixed
               Read(file_lu,* ,iostat=io_status)coef%mixedgasint(1,i),coef%mixedgasint(2,i)
               If(io_status /= 0) Then
                 Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
                 errorstatus = errorstatus_fatal 
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                 Return
               Endif
             Enddo

           Case( gas_id_watervapour )
            coef%nintwater=nbint
            allocate(coef%watervapourint(2,coef % nintwater),stat=alloc_status(1))
            If( Any(alloc_status /= 0) ) Then
              errorstatus = errorstatus_fatal
              Write( errMessage, '( "allocation of fmv coefs arrays")' )
              Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
              Return
            End If
             Do i=1,coef%nintwater
               Read(file_lu,* ,iostat=io_status)coef%watervapourint(1,i),coef%watervapourint(2,i)
               If(io_status /= 0) Then
                 Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
                 errorstatus = errorstatus_fatal
               Endif
             Enddo

           Case( gas_id_ozone )
           coef%nintozone=nbint
           allocate(coef%ozoneint(2,coef % nintozone),stat=alloc_status(1))
           If( Any(alloc_status /= 0) ) Then
             errorstatus = errorstatus_fatal
             Write( errMessage, '( "allocation of fmv coefs arrays")' )
             Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
             Return
           End If
           Do i=1,coef%nintozone
               Read(file_lu,* ,iostat=io_status)coef%ozoneint(1,i),coef%ozoneint(2,i)
               If(io_status /= 0) Then
                 Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
                 errorstatus = errorstatus_fatal
               Endif
           Enddo

           Case( gas_id_wvcont )
           coef%nintwvcont=nbint
           allocate(coef%wvcontint(2,coef % nintwvcont),stat=alloc_status(1))
           If( Any(alloc_status /= 0) ) Then
             errorstatus = errorstatus_fatal
             Write( errMessage, '( "allocation of fmv coefs arrays")' )
             Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
             Return
           End If
           Do i=1,coef%nintwvcont
               Read(file_lu,* ,iostat=io_status)coef%wvcontint(1,i),coef%wvcontint(2,i)
               If(io_status /= 0) Then
                 Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
                 errorstatus = errorstatus_fatal
               Endif
           Enddo

           Case( gas_id_co2 )
           coef%nintco2=nbint
           allocate(coef%co2int(2,coef % nintco2),stat=alloc_status(1))
           If( Any(alloc_status /= 0) ) Then
             errorstatus = errorstatus_fatal
             Write( errMessage, '( "allocation of fmv coefs arrays")' )
             Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
             Return
           End If
           Do i=1,coef%nintco2
               Read(file_lu,* ,iostat=io_status)coef%co2int(1,i),coef%co2int(2,i)
               If(io_status /= 0) Then
                 Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
                 errorstatus = errorstatus_fatal
               Endif
           Enddo

           Case( gas_id_n2o )
           coef%nintn2o=nbint
           allocate(coef%n2oint(2,coef % nintn2o),stat=alloc_status(1))
           If( Any(alloc_status /= 0) ) Then
             errorstatus = errorstatus_fatal
             Write( errMessage, '( "allocation of fmv coefs arrays")' )
             Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
             Return
           End If
           Do i=1,coef%nintn2o
               Read(file_lu,* ,iostat=io_status)coef%n2oint(1,i),coef%n2oint(2,i)
               If(io_status /= 0) Then
                 Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
                 errorstatus = errorstatus_fatal
               Endif
           Enddo

           Case( gas_id_co )
           coef%nintco=nbint
           allocate(coef%coint(2,coef % nintco),stat=alloc_status(1))
           If( Any(alloc_status /= 0) ) Then
             errorstatus = errorstatus_fatal
             Write( errMessage, '( "allocation of fmv coefs arrays")' )
             Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
             Return
           End If
           Do i=1,coef%nintco
               Read(file_lu,* ,iostat=io_status)coef%coint(1,i),coef%coint(2,i)
               If(io_status /= 0) Then
                 Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
                 errorstatus = errorstatus_fatal
               Endif
           Enddo

           Case( gas_id_ch4 )
           coef%nintch4=nbint
           allocate(coef%ch4int(2,coef % nintch4),stat=alloc_status(1))
           If( Any(alloc_status /= 0) ) Then
             errorstatus = errorstatus_fatal
             Write( errMessage, '( "allocation of fmv coefs arrays")' )
             Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
             Return
           End If
           Do i=1,coef%nintch4
              Read(file_lu,* ,iostat=io_status)coef%ch4int(1,i),coef%ch4int(2,i)
              If(io_status /= 0) Then
                Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
                errorstatus = errorstatus_fatal
              Endif
           Enddo

           End Select
        End Do





     Case ('FILTER_FUNCTIONS')
        Call rttov_skipcommentline ( file_lu, io_status )
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

        !allocate FILTER_FUNCTIONS section  array size is fmv_chn
        Allocate ( coef % ff_ori_chn( coef % fmv_chn ), stat=alloc_status(1) )
        Allocate ( coef % ff_val_chn( coef % fmv_chn ), stat=alloc_status(2) )
        Allocate ( coef % ff_cwn( coef % fmv_chn ), stat=alloc_status(3) )
        Allocate ( coef % ff_bco( coef % fmv_chn ), stat=alloc_status(4) )
        Allocate ( coef % ff_bcs( coef % fmv_chn ), stat=alloc_status(5) )
        Allocate ( coef % ff_gam( coef % fmv_chn ), stat=alloc_status(6) )
        If( Any(alloc_status /= 0) ) Then
           errorstatus = errorstatus_fatal
           Write( errMessage, '( "allocation of ff coefs arrays")' )
           Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        End If


        If ( all_channels ) Then
           Do i = 1, coef % fmv_chn
              Read( file_lu, *, iostat = io_status )&
                     & coef % ff_ori_chn(i), &
                     & coef % ff_val_chn(i), &
                     & coef % ff_cwn(i),     &
                     & coef % ff_bco(i),     &
                     & coef % ff_bcs(i),     &
                     & coef % ff_gam(i) 
              If(io_status /= 0) Then
                 errorstatus = errorstatus_fatal
                 Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
                 errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                 Return
              Endif
           End Do
        Else
           Allocate ( ivalues0( file_channels ), stat=alloc_status(1) )
           Allocate ( ivalues1( file_channels ), stat=alloc_status(2) )
           Allocate (  values0( file_channels ), stat=alloc_status(3) )
           Allocate (  values1( file_channels ), stat=alloc_status(4) )
           Allocate (  values2( file_channels ), stat=alloc_status(5) )
           Allocate (  values3( file_channels ), stat=alloc_status(6) )
           If( Any(alloc_status /= 0) ) Then
              errorstatus = errorstatus_fatal
              Write( errMessage, '( "allocation of ff coefs arrays")' )
              Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
              Return
           End If
           Do i = 1, file_channels
              Read( file_lu, *, iostat = io_status )&
                     & ivalues0(i),&
                     & ivalues1(i),&
                     & values0(i) ,&
                     & values1(i) ,&
                     & values2(i) ,&
                     & values3(i) 
              If(io_status /= 0) Then
                 Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
                 errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                 Return
              Endif
           End Do
           coef % ff_ori_chn(:) = ivalues0 ( channels(:) )
           coef % ff_val_chn(:) = ivalues1 ( channels(:) )
           coef % ff_cwn(:)     =  values0 ( channels(:) )
           coef % ff_bco(:)     =  values1 ( channels(:) )
           coef % ff_bcs(:)     =  values2 ( channels(:) )
           coef % ff_gam(:)     =  values3 ( channels(:) )
           Deallocate ( ivalues0, stat=alloc_status(1) )
           Deallocate ( ivalues1, stat=alloc_status(2) )
           Deallocate (  values0, stat=alloc_status(3) )
           Deallocate (  values1, stat=alloc_status(4) )
           Deallocate (  values2, stat=alloc_status(5) )
           Deallocate (  values3, stat=alloc_status(6) )
           If( Any(alloc_status /= 0) ) Then
              errorstatus = errorstatus_fatal
              Write( errMessage, '( "deallocation of ff coefs arrays")' )
              Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
              Return
           End If
        Endif

!       coef % ff_gam(:)=1.

     Case('TRANSMITTANCE_TRESHOLD')
        Call rttov_skipcommentline ( file_lu, io_status )
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

        Allocate ( coef % tt_chn      ( coef % fmv_chn ), stat=alloc_status(1) )
        Allocate ( coef % tt_val_chn  ( coef % fmv_chn ), stat=alloc_status(2) )
        Allocate ( coef % tt_cwn      ( coef % fmv_chn ), stat=alloc_status(3) )
        Allocate ( coef % tt_a0       ( coef % fmv_chn ), stat=alloc_status(4) )
        Allocate ( coef % tt_a1       ( coef % fmv_chn ), stat=alloc_status(5) )
If( Any(alloc_status /= 0) ) Then
           errorstatus = errorstatus_fatal
           Write( errMessage, '( "allocation of transmittance treshold coefs arrays")' )
           Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        End If

        If( all_channels ) Then
           Do i = 1, coef % fmv_chn

              ! chan number
              ! validity of channel 
              ! central wave number
              ! transmittance treshold
              ! transmittance value            
              Read(file_lu,*,iostat=io_status)&
                      & coef % tt_chn(i)    ,  &
                      & coef % tt_val_chn(i),  &
                      & coef % tt_cwn(i)    ,  &
                      & coef % tt_a0(i)     ,  &
                      & coef % tt_a1(i)          
              If(io_status /= 0) Then
                 Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
                 errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                 Return
              Endif

           End Do
        Else
           Allocate ( ivalues0( file_channels ), stat=alloc_status(1) )
           Allocate ( ivalues1( file_channels ), stat=alloc_status(2) )
           Allocate (  values0( file_channels ), stat=alloc_status(3) )
           Allocate (  values1( file_channels ), stat=alloc_status(4) )
           Allocate (  values2( file_channels ), stat=alloc_status(5) )
           If( Any(alloc_status /= 0) ) Then
              errorstatus = errorstatus_fatal
              Write( errMessage, '( "allocation of transmittance treshold coefs arrays")' )
              Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
              Return
           End If
           Do i = 1, file_channels
              Read( file_lu, *, iostat = io_status )&
                      & ivalues0(i),&
                      & ivalues1(i) ,&
                      & values0(i) ,&
                      & values1(i) ,&
                      & values2(i)  
              If(io_status /= 0) Then
                 errorstatus = errorstatus_fatal
                 Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
                 errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                 Return
              Endif
           End Do
           coef % tt_chn(:)        = ivalues0 ( channels(:) )
           coef % tt_val_chn(:)    = ivalues1 ( channels(:) )
           coef % tt_cwn(:)        =  values0 ( channels(:) )
           coef % tt_a0(:)         =  values1 ( channels(:) )
           coef % tt_a1(:)         =  values2 ( channels(:) )
           Deallocate ( ivalues0, stat=alloc_status(1) )
           Deallocate ( ivalues1, stat=alloc_status(2) )
           Deallocate (  values0, stat=alloc_status(3) )
           Deallocate (  values1, stat=alloc_status(4) )
           Deallocate (  values2, stat=alloc_status(5) )
           If( Any(alloc_status /= 0) ) Then
              errorstatus = errorstatus_fatal
              Write( errMessage, '( "deallocation of transmittance treshold coefs arrays")' )
              Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
              Return
           End If
        Endif

     Case('SOLAR_SPECTRUM')
        Call rttov_skipcommentline ( file_lu, io_status )
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

        Allocate ( coef % ss_chn            ( coef % fmv_chn ), stat=alloc_status(1) )
        Allocate ( coef % ss_val_chn        ( coef % fmv_chn ), stat=alloc_status(2) )
        Allocate ( coef % ss_cwn            ( coef % fmv_chn ), stat=alloc_status(3) )
        Allocate ( coef % ss_solar_spectrum ( coef % fmv_chn ), stat=alloc_status(4) )
        If( Any(alloc_status /= 0) ) Then
           errorstatus = errorstatus_fatal
           Write( errMessage, '( "allocation of solar spectrum coefs arrays")' )
           Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        End If

        If( all_channels ) Then
           Do i = 1, coef % fmv_chn

              ! chan number
              ! validity of channel
              ! central wave number
              ! solar spectrum     
              Read(file_lu,*,iostat=io_status)&
                      & coef % ss_chn(i)    ,  &
                      & coef % ss_val_chn(i),  &
                      & coef % ss_cwn(i)    ,  &
                      & coef % ss_solar_spectrum(i) 
              If(io_status /= 0) Then
                 Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
                 errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                 Return
              Endif

           End Do
        Else
           Allocate ( ivalues0( file_channels ), stat=alloc_status(1) )
           Allocate ( ivalues1( file_channels ), stat=alloc_status(2) )
           Allocate (  values0( file_channels ), stat=alloc_status(3) )
           Allocate (  values1( file_channels ), stat=alloc_status(4) )
           If( Any(alloc_status /= 0) ) Then
              errorstatus = errorstatus_fatal
              Write( errMessage, '( "allocation of solar spectrum coefs arrays")' )
              Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
              Return
           End If
           Do i = 1, file_channels
              Read( file_lu, *, iostat = io_status )&
                      & ivalues0(i),&
                      & ivalues1(i) ,&
                      & values0(i) ,&
                      & values1(i) 
              If(io_status /= 0) Then
                 errorstatus = errorstatus_fatal
                 Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
                 errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                 Return
              Endif
           End Do
           coef % ss_chn(:)            = ivalues0 ( channels(:) )
           coef % ss_val_chn(:)        = ivalues1 ( channels(:) )
           coef % ss_cwn(:)            =  values0 ( channels(:) )
           coef % ss_solar_spectrum(:) =  values1 ( channels(:) )
           Deallocate ( ivalues0, stat=alloc_status(1) )
           Deallocate ( ivalues1, stat=alloc_status(2) )
           Deallocate (  values0, stat=alloc_status(3) )
           Deallocate (  values1, stat=alloc_status(4) )
           If( Any(alloc_status /= 0) ) Then
              errorstatus = errorstatus_fatal
              Write( errMessage, '( "deallocation of solar spectrum coefs arrays")' )
              Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
              Return
           End If
        Endif

     Case('WATER_OPTICAL_CONSTANT')
        Call rttov_skipcommentline ( file_lu, io_status )
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal 
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

        Allocate ( coef % woc_chn            ( coef % fmv_chn ), stat=alloc_status(1) )
        Allocate ( coef % woc_cwn            ( coef % fmv_chn ), stat=alloc_status(2) )
        Allocate ( coef % woc_waopc_ow      ( coef % fmv_chn ), stat=alloc_status(3) )
        Allocate ( coef % woc_waopc_fw      ( coef % fmv_chn ), stat=alloc_status(4) )
        If( Any(alloc_status /= 0) ) Then
           errorstatus = errorstatus_fatal
           Write( errMessage, '( "allocation of water optical constants coefs arrays")' )
           Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        End If

        If( all_channels ) Then
           Do i = 1, coef % fmv_chn

              ! chan number
              ! central wave number
              ! ocean water optical constants(real and imaginary part)
              ! fresh water optical constants(real and imaginary part)     
              Read(file_lu,*,iostat=io_status)&
                      & coef % woc_chn(i)    ,   &
                      & coef % woc_cwn(i)    ,   &
                      & coef % woc_waopc_ow(i)  ,&
                      & coef % woc_waopc_fw(i) 
              If(io_status /= 0) Then
                 Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
                 errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                 Return
              Endif

           End Do
        Else
           Allocate ( ivalues0( file_channels ), stat=alloc_status(1) )
           Allocate ( values1( file_channels ), stat=alloc_status(2) )
           Allocate ( values_c_0( file_channels ), stat=alloc_status(3) )
           Allocate ( values_c_1( file_channels ), stat=alloc_status(4) )
           If( Any(alloc_status /= 0) ) Then
              errorstatus = errorstatus_fatal
              Write( errMessage, '( "allocation of solar spectrum coefs arrays")' )
              Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
              Return
           End If
           Do i = 1, file_channels
              Read( file_lu, *, iostat = io_status )&
                      & ivalues0(i),&
                      & values1(i) ,&
                      & values_c_0(i) ,&
                      & values_c_1(i) 
              If(io_status /= 0) Then
                 errorstatus = errorstatus_fatal
                 Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
                 errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                 Return
              Endif
           End Do
           coef % woc_chn(:)              = ivalues0 ( channels(:) )
           coef % woc_cwn(:)              = values1 ( channels(:) )
           coef % woc_waopc_ow(:)         = values_c_0 ( channels(:) )
           coef % woc_waopc_fw(:)         = values_c_1 ( channels(:) )
           Deallocate ( ivalues0, stat=alloc_status(1) )
           Deallocate ( values1, stat=alloc_status(2) )
           Deallocate ( values_c_0, stat=alloc_status(3) )
           Deallocate ( values_c_1, stat=alloc_status(4) )
           If( Any(alloc_status /= 0) ) Then
              errorstatus = errorstatus_fatal
              Write( errMessage, '( "deallocation of solar spectrum coefs arrays")' )
              Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
              Return
           End If
        Endif

     Case('WAVE_SPECTRUM')
        Call rttov_skipcommentline ( file_lu, io_status )
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

     Read(file_lu,*,iostat=io_status)coef%ws_nomega

     Allocate ( coef % ws_k_omega ( coef % ws_nomega ), stat=alloc_status(1) )
     Allocate ( coef % ws_npoint  ( coef % ws_nomega ), stat=alloc_status(2) )

        If( Any(alloc_status /= 0) ) Then
           errorstatus = errorstatus_fatal
           Write( errMessage, '( "allocation of solar spectrum coefs arrays")' )
           Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        End If

     Do i = 1,coef%ws_nomega
       Read(file_lu,*,iostat=io_status)coef%ws_npoint(i),coef%ws_k_omega(i)
              If(io_status /= 0) Then
                 Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
                 errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                 Return
              Endif
     enddo

     Case ('FUNDAMENTAL_CONSTANTS')
        Call rttov_skipcommentline ( file_lu, io_status )
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

        ! speed of light (cm/s)
        Read(file_lu,*,iostat=io_status) fc_speedl
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

        ! first radiation constant (mW/(m2*sr*cm-4))
        ! second radiation constant (cm*K)
        Read(file_lu,*,iostat=io_status) coef % fc_planck_c1, coef % fc_planck_c2
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

        ! satellite nominal altitude (km)
        Read(file_lu,*,iostat=io_status) coef % fc_sat_height
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

     Case ('FASTEM')
        Call rttov_skipcommentline ( file_lu, io_status )
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

        ! fastem version number
        Read(file_lu,*,iostat=io_status) coef % fastem_ver
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

        ! number of coefficients
        Read(file_lu,*,iostat=io_status) fastem_coef_nb
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif


        Allocate ( coef % fastem_polar ( coef % fmv_chn ), stat=alloc_status(1) )
        If( Any(alloc_status /= 0) ) Then
           errorstatus = errorstatus_fatal
           Write( errMessage, '( "allocation of fastem polarisation array")' )
           Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        End If

        ! coefficients (fastem_coef_nb) - no need to store these
        Read(file_lu,*,iostat=io_status) (fastem_coef, i= 1, fastem_coef_nb)
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

        ! polarisation of each channel
        If( all_channels ) Then
           Read(file_lu,*,iostat=io_status) (coef % fastem_polar(i), i= 1, coef % fmv_chn)
        Else
           Allocate ( ivalues0( file_channels ), stat=alloc_status(1) )
           If( Any(alloc_status /= 0) ) Then
              errorstatus = errorstatus_fatal
              Write( errMessage, '( "allocation of fastem polarisation array")' )
              Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
              Return
           End If
           Read(file_lu,*,iostat=io_status) ( ivalues0(i), i= 1, file_channels)
           coef % fastem_polar(:) = ivalues0( channels (:) )
           Deallocate ( ivalues0, stat=alloc_status(1) )
           If( Any(alloc_status /= 0) ) Then
              errorstatus = errorstatus_fatal
              Write( errMessage, '( "deallocation of fastem polarisation array")' )
              Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
              Return
           End If
        Endif
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif


        !-------------------------------------------------------
     Case ('SSIREM')
        Call rttov_skipcommentline ( file_lu, io_status )
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

        ! version number
        Read(file_lu,*,iostat=io_status) coef % ssirem_ver
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

        Allocate ( coef % ssirem_chn ( coef % fmv_chn ), stat=alloc_status(1) )
        Allocate ( coef % ssirem_a0  ( coef % fmv_chn ), stat=alloc_status(2) )
        Allocate ( coef % ssirem_a1  ( coef % fmv_chn ), stat=alloc_status(3) )
        Allocate ( coef % ssirem_a2  ( coef % fmv_chn ), stat=alloc_status(4) )
        Allocate ( coef % ssirem_xzn1( coef % fmv_chn ), stat=alloc_status(5) )
        Allocate ( coef % ssirem_xzn2( coef % fmv_chn ), stat=alloc_status(6) )
        If( Any(alloc_status /= 0) ) Then
           errorstatus = errorstatus_fatal
           Write( errMessage, '( "allocation of ssirem coefs arrays")' )
           Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        End If

        If( all_channels ) Then
           Do i = 1, coef % fmv_chn

              ! original chan number
              ! constant coef
              ! first order coef
              ! second order coef
              ! 1st exponent on zenith angle
              ! 2nd exponent on zenith angle
              Read(file_lu,*,iostat=io_status)&
                     & coef % ssirem_chn(i), &
                     & coef % ssirem_a0(i),  &
                     & coef % ssirem_a1(i),  &
                     & coef % ssirem_a2(i),  &
                     & coef % ssirem_xzn1(i),&
                     & coef % ssirem_xzn2(i) 
              If(io_status /= 0) Then
                 Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
                 errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                 Return
              Endif

           End Do
        Else
           Allocate ( ivalues0( file_channels ), stat=alloc_status(1) )
           Allocate (  values0( file_channels ), stat=alloc_status(2) )
           Allocate (  values1( file_channels ), stat=alloc_status(3) )
           Allocate (  values2( file_channels ), stat=alloc_status(4) )
           Allocate (  values3( file_channels ), stat=alloc_status(5) )
           Allocate (  values4( file_channels ), stat=alloc_status(6) )
           If( Any(alloc_status /= 0) ) Then
              errorstatus = errorstatus_fatal
              Write( errMessage, '( "allocation of ssirem coefs arrays")' )
              Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
              Return
           End If
           Do i = 1, file_channels
              Read( file_lu, *, iostat = io_status )&
                     & ivalues0(i),&
                     & values0(i) ,&
                     & values1(i) ,&
                     & values2(i) ,&
                     & values3(i) ,&
                     & values4(i) 
              If(io_status /= 0) Then
                 errorstatus = errorstatus_fatal
                 Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
                 errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                 Return
              Endif
           End Do
           coef % ssirem_chn(:)    = ivalues0 ( channels(:) )
           coef % ssirem_a0(:)     =  values0 ( channels(:) )
           coef % ssirem_a1(:)     =  values1 ( channels(:) )
           coef % ssirem_a2(:)     =  values2 ( channels(:) )
           coef % ssirem_xzn1(:)   =  values3 ( channels(:) )
           coef % ssirem_xzn2(:)   =  values4 ( channels(:) )
           Deallocate ( ivalues0, stat=alloc_status(1) )
           Deallocate (  values0, stat=alloc_status(2) )
           Deallocate (  values1, stat=alloc_status(3) )
           Deallocate (  values2, stat=alloc_status(4) )
           Deallocate (  values3, stat=alloc_status(5) )
           Deallocate (  values4, stat=alloc_status(6) )
           If( Any(alloc_status /= 0) ) Then
              errorstatus = errorstatus_fatal
              Write( errMessage, '( "deallocation of ssirem coefs arrays")' )
              Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
              Return
           End If
        Endif

        !-------------------------------------------------------
     Case ('GAZ_UNITS')
        ! the array has already been allocated and initialised
        ! to specific concentration (kg/kg)
        !
        ! This section needs one input line per gaz
        ! in the same order as the gaz list defined inside
        !
        ! This is defining the units used for the sections
        ! REFERENCE_PROFILE and PROFILE_LIMITS
        !
        Call rttov_skipcommentline ( file_lu, io_status )
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

        Do n = 1, coef % fmv_gas
           Call rttov_skipcommentline (file_lu, io_status)
           If(io_status /= 0) Then
              errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
              Return
           Endif

           Read( file_lu, *, iostat=io_status ) coef % gaz_units( n )
           If(io_status /= 0) Then
              Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
              errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
              Return
           Endif

        End Do

        !-------------------------------------------------------
     Case ('REFERENCE_PROFILE')
        Call rttov_skipcommentline ( file_lu, io_status )
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

        Allocate ( coef % ref_prfl_p ( coef % fmv_lvl(gas_id_mixed) ), stat=alloc_status(1) )
        Allocate ( coef % ref_prfl_t ( coef % fmv_lvl(gas_id_mixed), coef % fmv_gas ), stat=alloc_status(2) )
        Allocate ( coef % ref_prfl_mr( coef % fmv_lvl(gas_id_mixed), coef % fmv_gas ), stat=alloc_status(3) )
        If( Any(alloc_status /= 0) ) Then
           Write( errMessage, '( "allocation of ref profile coefs arrays")' )
           Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        End If

        Do n = 1, coef % fmv_gas
           Call rttov_skipcommentline (file_lu, io_status)
           If(io_status /= 0) Then
              errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
              Return
           Endif

           ! units for reference gaz concentration is
           ! specified in GAZ_UNITS section (default is specific concentration (kg/kg))
           !
           Do i = 1, coef % nlevels
              Read( file_lu, *, iostat=io_status ) &
                     & pres                        ,&
                     & coef % ref_prfl_t  ( i, n )  ,&
                     & coef % ref_prfl_mr ( i, n ) 
              If(io_status /= 0) Then
                 Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
                 errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                 Return
              Endif
              If( coef % fmv_gas_id(n) == gas_id_mixed) coef % ref_prfl_p( i ) = pres
           End Do

        End Do

        !-------------------------------------------------------
     Case ('PROFILE_LIMITS')
        Call rttov_skipcommentline ( file_lu, io_status )
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

        Allocate ( coef % lim_prfl_p( coef % fmv_lvl(gas_id_mixed) ), stat=alloc_status(1) )
        Allocate ( coef % lim_prfl_tmax( coef % fmv_lvl(gas_id_mixed) ), stat=alloc_status(2) )
        Allocate ( coef % lim_prfl_tmin( coef % fmv_lvl(gas_id_mixed) ), stat=alloc_status(3) )
        Allocate ( coef % lim_prfl_gmin( coef % fmv_lvl(gas_id_mixed), coef % fmv_gas ), stat=alloc_status(4) )
        Allocate ( coef % lim_prfl_gmax( coef % fmv_lvl(gas_id_mixed), coef % fmv_gas ), stat=alloc_status(5) )
        If( Any(alloc_status /= 0) ) Then
           errorstatus = errorstatus_fatal
           Write( errMessage, '( "allocation of profile limits arrays")' )
           Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        End If

        Do l = 1, coef % nlevels
           ! pressure  (hPa)       (levels)
           ! max temperature (K)   (levels)
           ! min temperature (K)   (levels)
           Read(file_lu,*,iostat=io_status)&
                  & coef % lim_prfl_p(l), coef % lim_prfl_tmax(l), coef % lim_prfl_tmin(l) 
           If(io_status /= 0) Then
              Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
              errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
              Return
           Endif
        End Do

        Do n = 1, coef % fmv_gas
           Call rttov_skipcommentline ( file_lu, io_status )
           If(io_status /= 0) Then
              Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
              errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
              Return
           Endif

           Do i = 1, coef % nlevels
              ! max specific concentration (kg/kg) (levels, gases)
              ! min specific concentration (kg/kg) (levels, gases)
              ! or
              ! max volume mixing r (ppmv) (levels, gases)
              ! min volume mixing r (ppmv) (levels, gases)
              ! according to
              ! units specified in GAZ_UNITS section (default is specific concentration (kg/kg))
              Read(file_lu,*,iostat=io_status) &
                     & pres, coef % lim_prfl_gmax( i , n ), coef % lim_prfl_gmin( i , n) 
              If(io_status /= 0) Then
                 Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
                 errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                 Return
              Endif
           End Do

        End Do

        !-------------------------------------------------------
     Case ('FAST_COEFFICIENTS','COEF_SUB_FILES')
        Call rttov_skipcommentline ( file_lu, io_status )
        If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
           Return
        Endif

        ! If section is COEF_SUB_FILES then coefficients for each gaz is stored
        ! in separate files.
        ! This possibility could be used to store very large coefficient files
        ! (large number of channels or gases)
        ! Section contains 1 line per gaz in the same order as the
        ! FAST_MODEL_VARIABLES section
        ! line indicates the filename of the coefficient for that gas
        !
        ! No verification is done on the file.
        ! header lines starting with "!" sign are skipped

        ALLOCATE (coef%thermal, STAT = alloc_status(1))
        If( Any(alloc_status /= 0) ) Then
            errorstatus = errorstatus_fatal
            Write( errMessage, '( "allocation of thermal fast coefs")' )
            Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
            Return
        End If
        NULLIFY (coef%thermal%mixedgas)
        NULLIFY (coef%thermal%watervapour)
        NULLIFY (coef%thermal%ozone)
        NULLIFY (coef%thermal%wvcont)
        NULLIFY (coef%thermal%co2)
        NULLIFY (coef%thermal%n2o)
        NULLIFY (coef%thermal%co)
        NULLIFY (coef%thermal%ch4)

        ! loop on gases
        Do n = 1, coef % fmv_gas
           Call rttov_skipcommentline ( file_lu, io_status )
           If(io_status /= 0) Then
              Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
              errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
              Return
           Endif

           ! read dummy string of gas name or filename of the sub_coefficient file
           Read( file_lu, *, iostat=io_status ) input_string
           If(io_status /= 0) Then
              Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
              errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
              Return
           Endif

           ! Case of Sub coefficient files
           ! Open the file and skip the header
           If( Trim(section) == 'COEF_SUB_FILES') Then
              file_lu_coef = 0
              for_output = .False.
              Call rttov789_opencoeff (errorstatus, input_string, file_lu_coef,  for_output)
              If ( errorstatus /= 0 ) Then
                 errorstatus = errorstatus_fatal
                 Write( errMessage, '( "Error opening sub_coef file" )' )
                 Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                 Return
              Endif
              Call rttov_skipcommentline ( file_lu_coef, io_status )
              If(io_status /= 0) Then
                 Write( errMessage, &
                        & '( "I/O error while reading sub_coef file ",i5 )' ) io_status 
                 Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
                 errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                 Return
              Endif
           Else
              file_lu_coef = file_lu
           Endif

           Select Case( coef % fmv_gas_id(n) )

           Case(gas_id_mixed)
              nvalues      = coef % ncmixed
              Allocate ( coef % thermal % mixedgas ( coef % nlevels, coef % fmv_chn, coef % ncmixed), stat= alloc_status(1))
              If( Any(alloc_status /= 0) ) Then
                 errorstatus = errorstatus_fatal
                 Write( errMessage, '( "allocation of MxG coefs arrays")' )
                 Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                 Return
              End If
              If( all_channels ) Then
                 coeffsarray => coef % thermal % mixedgas
              Else
                 Allocate( coeffsarray( coef % nlevels, file_channels, coef % ncmixed ), stat= alloc_status(1))
                 If( Any(alloc_status /= 0) ) Then
                    errorstatus = errorstatus_fatal
                    Write( errMessage, '( "allocation of MxG coefs arrays")' )
                    Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                    Return
                 End If

              Endif

           Case(gas_id_watervapour)
              nvalues      = coef % ncwater
              Allocate ( coef % thermal % watervapour ( coef % nlevels, coef % fmv_chn, coef % ncwater), stat= alloc_status(1))
              If( Any(alloc_status /= 0) ) Then
                 errorstatus = errorstatus_fatal
                 Write( errMessage, '( "allocation of WV coefs arrays")' )
                 Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                 Return
              End If
              If( all_channels ) Then
                 coeffsarray => coef % thermal % watervapour
              Else
                 Allocate( coeffsarray( coef % nlevels, file_channels, coef % ncwater ), stat= alloc_status(1))
                 If( Any(alloc_status /= 0) ) Then
                    errorstatus = errorstatus_fatal
                    Write( errMessage, '( "allocation of WV coefs arrays")' )
                    Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                    Return
                 End If
              Endif

           Case(gas_id_ozone)
              nvalues      = coef % ncozone
              Allocate ( coef % thermal % ozone ( coef % nlevels, coef % fmv_chn, coef % ncozone), stat= alloc_status(1))
              If( Any(alloc_status /= 0) ) Then
                 errorstatus = errorstatus_fatal
                 Write( errMessage, '( "allocation of O3 coefs arrays")' )
                 Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                 Return
              End If
              If( all_channels ) Then
                 coeffsarray => coef % thermal % ozone
              Else
                 Allocate( coeffsarray( coef % nlevels, file_channels, coef % ncozone ), stat= alloc_status(1))
                 If( Any(alloc_status /= 0) ) Then
                    errorstatus = errorstatus_fatal
                    Write( errMessage, '( "allocation of O3 coefs arrays")' )
                    Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                    Return
                 End If
              Endif

           Case(gas_id_wvcont)
              nvalues      = coef % ncwvcont
              Allocate ( coef % thermal % wvcont ( coef % nlevels, coef % fmv_chn, coef % ncwvcont), stat= alloc_status(1))
              If( Any(alloc_status /= 0) ) Then
                 errorstatus = errorstatus_fatal
                 Write( errMessage, '( "allocation of WV continuum coefs arrays")' )
                 Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                 Return
              End If
              If( all_channels ) Then
                 coeffsarray => coef % thermal % wvcont
              Else
                 Allocate( coeffsarray( coef % nlevels, file_channels, coef % ncwvcont ), stat= alloc_status(1))
                 If( Any(alloc_status /= 0) ) Then
                    errorstatus = errorstatus_fatal
                    Write( errMessage, '( "allocation of WV continuum coefs arrays")' )
                    Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                    Return
                 End If
              Endif

           Case(gas_id_co2)
              nvalues      = coef % ncco2
              Allocate ( coef % thermal % co2 ( coef % nlevels, coef % fmv_chn, coef % ncco2), stat= alloc_status(1))
              If( Any(alloc_status /= 0) ) Then
                 errorstatus = errorstatus_fatal
                 Write( errMessage, '( "allocation of CO2 coefs arrays")' )
                 Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                 Return
              End If
              If( all_channels ) Then
                 coeffsarray => coef % thermal % co2
              Else
                 Allocate( coeffsarray( coef % nlevels, file_channels, coef % ncco2 ), stat= alloc_status(1))
                 If( Any(alloc_status /= 0) ) Then
                    errorstatus = errorstatus_fatal
                    Write( errMessage, '( "allocation of CO2 coefs arrays")' )
                    Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                    Return
                 End If
              Endif

           Case(gas_id_n2o)
              nvalues      = coef % ncn2o
              Allocate ( coef % thermal % n2o ( coef % nlevels, coef % fmv_chn, coef % ncn2o), stat= alloc_status(1))
              If( Any(alloc_status /= 0) ) Then
                 errorstatus = errorstatus_fatal
                 Write( errMessage, '( "allocation of N2O coefs arrays")' )
                 Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                 Return
              End If
              If( all_channels ) Then
                 coeffsarray => coef % thermal % n2o
              Else
                 Allocate( coeffsarray( coef % nlevels, file_channels, coef % ncn2o ), stat= alloc_status(1))
                 If( Any(alloc_status /= 0) ) Then
                    errorstatus = errorstatus_fatal
                    Write( errMessage, '( "allocation of N2O coefs arrays")' )
                    Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                    Return
                 End If
              Endif

           Case(gas_id_co)
              nvalues      = coef % ncco
              Allocate ( coef % thermal % co ( coef % nlevels, coef % fmv_chn, coef % ncco), stat= alloc_status(1))
              If( Any(alloc_status /= 0) ) Then
                 errorstatus = errorstatus_fatal
                 Write( errMessage, '( "allocation of CO coefs arrays")' )
                 Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                 Return
              End If
              If( all_channels ) Then
                 coeffsarray => coef % thermal % co
              Else
                 Allocate( coeffsarray( coef % nlevels, file_channels, coef % ncco ), stat= alloc_status(1))
                 If( Any(alloc_status /= 0) ) Then
                    errorstatus = errorstatus_fatal
                    Write( errMessage, '( "allocation of CO coefs arrays")' )
                    Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                    Return
                 End If
              Endif

           Case(gas_id_ch4)
              nvalues      = coef % ncch4
              Allocate ( coef % thermal % ch4 ( coef % nlevels, coef % fmv_chn, coef % ncch4), stat= alloc_status(1))
              If( Any(alloc_status /= 0) ) Then
                 errorstatus = errorstatus_fatal
                 Write( errMessage, '( "allocation of CH4 coefs arrays")' )
                 Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                 Return
              End If
              If( all_channels ) Then
                 coeffsarray => coef % thermal % ch4
              Else
                 Allocate( coeffsarray( coef % nlevels, file_channels, coef % ncch4 ), stat= alloc_status(1))
                 If( Any(alloc_status /= 0) ) Then
                    errorstatus = errorstatus_fatal
                    Write( errMessage, '( "allocation of CH4 coefs arrays")' )
                    Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                    Return
                 End If
              Endif


           End Select

           Read( file_lu_coef, *, iostat=io_status ) (((coeffsarray(i,j,k), &
                  & i = 1, coef % nlevels),&
                  & j = 1, file_channels),&
                  & k = 1, nvalues) 
           If(io_status /= 0) Then
              errmessage = 'erreur lecture'
              Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
              errorstatus = errorstatus_fatal
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
              Return
           Endif

           If( .Not. all_channels ) Then
              Select Case( coef % fmv_gas_id(n) )

              Case(gas_id_mixed)
                 coef % thermal % mixedgas( : , : , : ) = coeffsarray( : , channels(:) , : )
              Case(gas_id_watervapour)
                 coef % thermal % watervapour( : , : , : ) = coeffsarray( : , channels(:) , : )
              Case(gas_id_ozone)
                 coef % thermal % ozone( : , : , : ) = coeffsarray( : , channels(:) , : )
              Case(gas_id_wvcont)
                 coef % thermal % wvcont( : , : , : ) = coeffsarray( : , channels(:) , : )
              Case(gas_id_co2)
                 coef % thermal % co2( : , : , : ) = coeffsarray( : , channels(:) , : )
              Case(gas_id_n2o)
                 coef % thermal % n2o( : , : , : ) = coeffsarray( : , channels(:) , : )
              Case(gas_id_co)
                 coef % thermal % co( : , : , : ) = coeffsarray( : , channels(:) , : )
              Case(gas_id_ch4)
                 coef % thermal % ch4( : , : , : ) = coeffsarray( : , channels(:) , : )
                 !
              End Select
              Deallocate ( coeffsarray, stat=alloc_status(1) )
              If( Any(alloc_status /= 0) ) Then
                 errorstatus = errorstatus_fatal
                 Write( errMessage, '( "deallocation of intermediate coefs array")' )
                 Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
                 Return
              End If
           Endif

           ! For COEF_SUB_FILES close the intermediate coef file
           If( Trim(section) == 'COEF_SUB_FILES') Then
              Close( unit = file_lu_coef )
           Endif

        End Do

     Case ('END')
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)
        Return

     Case default

        Cycle readfile

     End Select

  End Do readfile
IF (LHOOK) CALL DR_HOOK('RTTOV_READCOEFFS_ASCII',1_jpim,ZHOOK_HANDLE)

End Subroutine rttov789_readcoeffs_ascii
