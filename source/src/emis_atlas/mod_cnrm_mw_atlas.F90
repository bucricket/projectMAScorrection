!
MODULE mod_cnrm_mw_atlas
  ! Description:
  !   Data and routines for MW emissivity atlas METEO-FRANCE CNRM.
  ! 
  !   Following hypothesis is done:
  !   Altas is ordered in latitude, longitude:
  !     all longitudes for a given latitude, next latitude etc...
  !
  ! Karbou, F., E. Gérard, and F. Rabier, 2006, 
  ! Microwave land emissivity and skin temperature for AMSU-A and -B assimilation over land,
  ! Q. J. R.Meteorol. Soc., vol. 132, No. 620, Part A, pp.2333-2355(23). 
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
  !  1.0      27/05/2010  Created, based on IFS/Arpege code (J. Hocking P. Brunel)
  !           22/11/2010  Change 30deg zenith angle threshold to 40deg (P. Brunel)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:
  !
#include "throw.h"

  USE parkind1, Only : jpim, jprb, jplm
  Use rttov_types, Only : rttov_coef
  Use rttov_const, Only :  &
        errorstatus_fatal, &
        inst_id_amsua,     &
        inst_id_amsub,     &
        inst_id_mhs,       &
        inst_id_ssmi,      &
        inst_id_ssmis,     &
        inst_name

  ! Disable implicit typing
  IMPLICIT NONE

#include "rttov_errorreport.interface"

  INTEGER(KIND=jpim), PARAMETER :: mw_atlas_version=200   ! Version of atlas
  
  ! Atlas constants (extracted from MODULE YOMEMIS of IFS/ARPEGE cy36)
   
  ! Emissivity from an atlas
  ! - low angles (an1) and high angles (an2)
  ! the parameters should be updated in case the emissivity atlases change
 
  INTEGER(KIND=JPIM) , PARAMETER :: TOT_EM_AN1 = 352207
  REAL(KIND=JPRB) ,DIMENSION(12,TOT_EM_AN1)   ::  DATA_EMIS_AN1
  integer(KIND=JPIM) :: lat_offset_an1(-90:90)
 
  INTEGER(KIND=JPIM) , PARAMETER :: TOT_EM_AN2 = 352207
  REAL(KIND=JPRB) ,DIMENSION(12,TOT_EM_AN2)   ::  DATA_EMIS_AN2
  integer(KIND=JPIM) :: lat_offset_an2(-90:90)

  INTEGER(KIND=JPIM) , PARAMETER :: TOT_EM_SSMIS = 186198
  REAL(KIND=JPRB) ,DIMENSION(20,TOT_EM_SSMIS)   ::  ATLAS_SSMIS
  integer(KIND=JPIM) :: lat_offset_ssmis(-90:90)

  INTEGER(KIND=JPIM) , PARAMETER :: TOT_EM_SSMI = 322156
  REAL(KIND=JPRB) ,DIMENSION(18,TOT_EM_SSMI)   ::  ATLAS_SSMI
  integer(KIND=JPIM) :: lat_offset_ssmi(-90:90)
 
CONTAINS

!------------------------------------------
! Routines for initialising database
!------------------------------------------
  
  SUBROUTINE rttov_cnrmmwemis_init( &
        &             path,       &! in
        &             imonth,     &! in
        &             coef,       &! in
        &             verbose,    &! in
        &             err         )! out
                   
    Implicit None
#include "rttov_opencoeff.interface"

    Character (len=*),  Intent(in)  :: path
    Integer(Kind=jpim), Intent(in)  :: imonth
    TYPE(rttov_coef),   INTENT(IN)  :: coef
    Logical(Kind=jplm), Intent(in)  :: verbose
    Integer(Kind=jpim), Intent(out) :: err

    
    Integer(Kind=jpim) ::  file_id
    Character (len=300) :: fn
    Character (len=2) :: cmonth

    Logical(Kind=jplm) :: file_exists
TRY
    write(cmonth,"(i2.2)") imonth
     
    !----------------------------------------------------------------------------

    IF ( coef%id_inst == inst_id_amsua .or. &
       & coef%id_inst == inst_id_amsub .or. &
       & coef%id_inst == inst_id_mhs   ) THEN

      fn=path//'AMSU_CNRM_'//cmonth//'_atlas_an1.dat'
      Inquire(FILE=fn, EXIST=file_exists)
      If ( .Not. file_exists ) err = errorstatus_fatal
      THROWM( ERR .NE. 0, "CNRM Atlas file"//TRIM(fn)//" not found" )

      file_id = 0
      Call rttov_opencoeff ( ERR, fn, file_id ) 
      THROWM( ERR .NE. 0, "Open CNRM Atlas file"//TRIM(fn) )

      If (verbose) INFO( "Using CNRM means "//TRIM(fn))

      READ(file_id,*, iostat=ERR) DATA_EMIS_AN1
      THROWM(err.ne.0,"io status while reading DATA_EMIS_AN1 ")

      CLOSE (unit=file_id)

      call rttov_cnrmmwemis_lat_limits ( data_emis_an1, lat_offset_an1 ) 


      fn=path//'AMSU_CNRM_'//cmonth//'_atlas_an2.dat'
      Inquire(FILE=fn, EXIST=file_exists)
      If ( .Not. file_exists ) err = errorstatus_fatal
      THROWM( ERR .NE. 0, "CNRM Atlas file"//TRIM(fn)//" not found" )

      file_id = 0
      Call rttov_opencoeff ( ERR, fn, file_id ) 
      THROWM( ERR .NE. 0, "Open CNRM Atlas file"//TRIM(fn) )

      INFO( "Using CNRM means "//TRIM(fn))

      READ(file_id,*, iostat=ERR) DATA_EMIS_AN2
      THROWM(err.ne.0,"io status while reading DATA_EMIS_AN2 ")

      CLOSE (unit=file_id)

      call rttov_cnrmmwemis_lat_limits ( data_emis_an2, lat_offset_an2 )
 

    ELSE IF ( coef%id_inst == inst_id_ssmi  ) THEN

      fn=path//'SSMI_CNRM_'//cmonth//'_atlas_ssmi.dat'
      Inquire(FILE=fn, EXIST=file_exists)
      If ( .Not. file_exists ) err = errorstatus_fatal
      THROWM( ERR .NE. 0, "CNRM Atlas file"//TRIM(fn)//" not found" )

      file_id = 0
      Call rttov_opencoeff ( ERR, fn, file_id ) 
      THROWM( ERR .NE. 0, "Open CNRM Atlas file"//TRIM(fn) )

      If (verbose) INFO( "Using CNRM means "//TRIM(fn))

      READ(file_id,*, iostat=ERR) ATLAS_SSMI
      THROWM(err.ne.0,"io status while reading ATLA_SSMI ")

      CLOSE (unit=file_id)

      call rttov_cnrmmwemis_lat_limits ( atlas_ssmi, lat_offset_ssmi )


!    ELSE IF ( coef%id_inst == inst_id_ssmis  ) THEN
!
!      fn=path//'SSMIS_CNRM_'//cmonth//'_atlas_ssmis.dat'
!      Inquire(FILE=fn, EXIST=file_exists)
!      If ( .Not. file_exists ) err = errorstatus_fatal
!      THROWM( ERR .NE. 0, "CNRM Atlas file"//TRIM(fn)//" not found" )
!
!      file_id = 0
!      Call rttov_opencoeff ( ERR, fn, file_id ) 
!      THROWM( ERR .NE. 0, "Open CNRM Atlas file"//TRIM(fn) )
!
!      If (verbose) INFO( "Using CNRM means "//TRIM(fn))
!
!      READ(file_id,*, iostat=ERR) ATLAS_SSMIS
!      THROWM(err.ne.0,"io status while reading ATLA_SSMIS ")
!
!      CLOSE (unit=file_id)
!
!      call rttov_cnrmmwemis_lat_limits ( atlas_ssmis, lat_offset_ssmis )


    ELSE

      ERR = 1
      THROWM(err.ne.0,"Instrument not supported "//Trim(inst_name(coef%id_inst)))

    END IF


CATCH
  END SUBROUTINE rttov_cnrmmwemis_init
  
  SUBROUTINE rttov_cnrmmwemis(    &
                & ERR,            &
                & id_inst,        &! in
                & nchan,          &! in 
                & latitude,       &! in 
                & longitude_in,   &! in
                & zenangle,       &! in
                & ori_chn,        &! in 
                & emis_out,       &! out 
                & pbats_veg)       ! out 
    Implicit None

    Integer(Kind=jpim), Intent(out) :: err         ! Error code
    Integer(Kind=jpim), Intent(in)  :: id_inst     ! ID instrument
    Integer(Kind=jpim), Intent(in)  :: nchan       ! number of channels
    Real(Kind=jprb),    Intent(in)  :: latitude    ! degrees
    Real(Kind=jprb),    Intent(in)  :: longitude_in   ! degrees
    Real(Kind=jprb),    Intent(in)  :: zenangle    ! zenith angle (degrees)
    Integer(Kind=jpim), Intent(in)  :: ori_chn(:)  ! instrument channel numbers
    Real(Kind=jprb),    Intent(out) :: emis_out(:) ! emissivity values 
    Real(Kind=jprb),    Intent(out) :: pbats_veg   ! vegetation
   
    Real(Kind=jprb)    :: longitude ! in ]-180,180]

    Real(Kind=jprb)    :: emis_amsua(15)
    Real(Kind=jprb)    :: emis_amsub(5)
    Real(Kind=jprb)    :: emis_ssmi(7)
    Real(Kind=jprb)    :: emis_ssmis(24)

    Real(Kind=jprb)    :: ze_ret3(4)
    Real(Kind=jprb)    :: ze_ret1
    Real(Kind=jprb)    :: zssmi_eret(7)
    Real(Kind=jprb)    :: zssmis_eret(8)
TRY

    IF (longitude_in .gt. 180.0_jprb) Then
      longitude = longitude_in - 360.0_jprb
    ELSEIF (longitude_in .le. -180.0_jprb) Then
      longitude = longitude_in + 360.0_jprb
    ELSE
      longitude = longitude_in
    ENDIF

    IF ( id_inst == inst_id_amsua ) THEN
      IF (zenangle > 40) THEN
        CALL LAND_AMSUA_AN2(latitude,longitude,ZE_RET3,PBATS_VEG)
      ELSE
        CALL LAND_AMSUA_AN1(latitude,longitude,ZE_RET3,PBATS_VEG)
      ENDIF

      emis_amsua(1)    = ZE_RET3(1)
      emis_amsua(2)    = ZE_RET3(2)
      emis_amsua(3:14) = ZE_RET3(3)
      emis_amsua(15)   = ZE_RET3(4)

      emis_out(1:nchan) = emis_amsua(ori_chn(1:nchan))

    ELSE IF ( id_inst == inst_id_amsub .or. &
            & id_inst == inst_id_mhs   ) THEN
      IF (zenangle > 40) THEN
        CALL LAND_AMSUB_AN2 (latitude,longitude,ZE_RET1,PBATS_VEG) 
      ELSE
        CALL LAND_AMSUB_AN1 (latitude,longitude,ZE_RET1,PBATS_VEG) 
      ENDIF

      emis_amsub(:) = ZE_RET1

      emis_out(1:nchan) = emis_amsub(ori_chn(1:nchan))

    ELSE IF ( id_inst == inst_id_ssmi  ) THEN
      CALL LAND_SSMI(latitude,longitude,ZSSMI_ERET,PBATS_VEG)

      emis_ssmi(:)  = ZSSMI_ERET(:)

      emis_out(1:nchan) = emis_ssmi(ori_chn(1:nchan))
      
    ELSE IF ( id_inst == inst_id_ssmis  ) THEN
      CALL LAND_SSMIS(latitude,longitude,ZSSMIS_ERET,PBATS_VEG)

      emis_ssmis(1:7)   = ZSSMIS_ERET(6)
      emis_ssmis(8:11)  = ZSSMIS_ERET(8)
      emis_ssmis(12)    = ZSSMIS_ERET(2)
      emis_ssmis(13)    = ZSSMIS_ERET(1)
      emis_ssmis(14)    = ZSSMIS_ERET(3)
      emis_ssmis(15)    = ZSSMIS_ERET(5)
      emis_ssmis(16)    = ZSSMIS_ERET(4)
      emis_ssmis(17)    = ZSSMIS_ERET(7)
      emis_ssmis(18)    = ZSSMIS_ERET(8)
      emis_ssmis(19:24) = ZSSMIS_ERET(6)

      emis_out(1:nchan) = emis_ssmis(ori_chn(1:nchan))

    ELSE

      ERR = 1
      THROWM(err.ne.0,"Instrument not supported "//Trim(inst_name(id_inst)))

    END IF

CATCH

  END SUBROUTINE rttov_cnrmmwemis

  SUBROUTINE rttov_cnrmmwemis_lat_limits ( data_emis, lat_offset )
  ! calculates the offset in the data_emis array
  ! for all interger latitude value
    Implicit None
    real(kind=jprb),    intent(in)  :: data_emis(:,:)
    integer(kind=jpim), intent(out) :: lat_offset(-90:90)
    integer :: last_lat, new_lat, jj

    lat_offset(:) = 1_jpim
    last_lat      = -90_jpim

    do jj = 1, SIZE(DATA_EMIS(1,:))
      new_lat = int(DATA_EMIS(2,jj) + 90.0_jprb) -90_jpim
      if( new_lat > last_lat ) then
        lat_offset(last_lat+1:new_lat) = jj
        last_lat = new_lat
      end if
    end do

    if (last_lat < 90_jpim ) then
    lat_offset(last_lat+1:) = lat_offset(last_lat) 
    end if

  END SUBROUTINE rttov_cnrmmwemis_lat_limits

!=======================================================================
!
! Following are IFS/ARPEGE routines cleaned for:
! Dr_HOOK 
! Unused variables and arguments
! Commented code
! efficient loops with latitude limits (use of lat_offset_xxx)
!
!=======================================================================
!
SUBROUTINE LAND_AMSUA_AN1 (PLAT,PLON,E_RET,PBATS_VEG)

! ROUTINE.
! --------
!    LAND_MF_AN1

! PURPOSE.
! --------
!     ROUTINE CALLED BY IFS OBSERVATION SCREENING MODULE TO 
!     INTERPOLATE SURFACE EMISSIVITY FROM AN ATLAS

! INPUT.
! ------
!     
!     PLAT - observation latitude
!     PLON - observation longitude


! OUTPUT.
! -------
!     E_RET- surface emissivity     
!     PBATS_VEG- BATS vegetation type

! F Karbou, N Bormann   09/05/07 : to caclucate land emissivities at mw frequencies   
! should be modified according to the emissivity atlas file format
    
USE PARKIND1  ,ONLY : JPIM     ,JPRB 
       
IMPLICIT NONE
      
REAL(KIND=JPRB), INTENT(IN)  :: PLAT
REAL(KIND=JPRB), INTENT(IN)  :: PLON
REAL(KIND=JPRB), INTENT(OUT) :: E_RET (4)
REAL(KIND=JPRB), INTENT(OUT) :: PBATS_VEG
 
INTEGER(KIND=JPIM) indice, jj
INTEGER(KIND=JPIM) lo, hi
REAL(KIND=JPRB) r1,r2,r3,r4,dx,dy
          
! find within AMSU mean emissivity map the closet emissivity sets to the observation
! 
          
r1 = 0.0
r2 = 0.0
r3 = 0.0
r4 = 0.0
indice=0

lo = 1
hi = SIZE(DATA_EMIS_AN1(1,:)) 
if ( nint(plat) > -90 ) lo = lat_offset_an1( nint(plat)-1 )
if ( nint(plat) <  90 ) hi = lat_offset_an1( nint(plat)+1 )

do jj = lo, hi

  dx = abs(DATA_EMIS_AN1(1,jj)-PLON)
  dy = abs(DATA_EMIS_AN1(2,jj)-PLAT)
        
  IF ((dx <= 0.5).and.(dy <= 0.5)) then
    IF ((DATA_EMIS_AN1(3,jj) > 0).AND.(DATA_EMIS_AN1(5,jj) > 0).AND. &
      & (DATA_EMIS_AN1(7,jj) > 0).AND.(DATA_EMIS_AN1(9,jj) > 0)) THEN
       indice = indice+1
       PBATS_VEG=DATA_EMIS_AN1(11,jj)
       IF (PBATS_VEG == 21) PBATS_VEG=14
       r1 = r1 + DATA_EMIS_AN1(3,jj)
       r2 = r2 + DATA_EMIS_AN1(5,jj)
       r3 = r3 + DATA_EMIS_AN1(7,jj)
       r4 = r4 + DATA_EMIS_AN1(9,jj)
    endif
 ENDIF

enddo
if (indice > 0) then
  E_RET(1)=(r1/indice)
  E_RET(2)=(r2/indice)
  E_RET(3)=(r3/indice)
  E_RET(4)=(r4/indice)   

else
  PBATS_VEG=-9.09
  if (abs(PLAT).gt.78) THEN
    E_RET(1)=0.75
    E_RET(2)=0.75
    E_RET(3)=0.75
    E_RET(4)=0.75
  else
    E_RET(1)=0.95005
    E_RET(2)=0.95005
    E_RET(3)=0.95005
    E_RET(4)=0.95005
  end if      
endif

END SUBROUTINE LAND_AMSUA_AN1

SUBROUTINE LAND_AMSUA_AN2 (PLAT,PLON,E_RET,PBATS_VEG)

! ROUTINE.
! --------
!    LAND_MF_AN1

! PURPOSE.
! --------
!     ROUTINE CALLED BY IFS OBSERVATION SCREENING MODULE TO 
!     INTERPOLATE SURFACE EMISSIVITY FROM AN ATLAS

! INPUT.
! ------
!     
!     PLAT - observation latitude
!     PLON - observation longitude


! OUTPUT.
! -------
!     E_RET- surface emissivity     
!     PBATS_VEG- BATS vegetation type

! F Karbou, N Bormann   09/05/07 : to caclucate land emissivities at mw frequencies   
! should be modified according to the emissivity atlas file format
         
USE PARKIND1  ,ONLY : JPIM     ,JPRB 
       
IMPLICIT NONE
      
REAL(KIND=JPRB), INTENT(IN)  :: PLAT
REAL(KIND=JPRB), INTENT(IN)  :: PLON
REAL(KIND=JPRB), INTENT(OUT) :: PBATS_VEG
REAL(KIND=JPRB), INTENT(OUT) :: E_RET (4)
     
INTEGER(KIND=JPIM) indice, jj
INTEGER(KIND=JPIM) lo, hi

REAL(KIND=JPRB) r1,r2,r3,r4,dx,dy


! find within AMSU mean emissivity map the closet emissivity sets to the observation
! 
r1 = 0.0
r2 = 0.0
r3 = 0.0
r4 = 0.0
indice=0

lo = 1
hi = SIZE(DATA_EMIS_AN2(1,:)) 
if ( nint(plat) > -90 ) lo = lat_offset_an1( nint(plat)-1 )
if ( nint(plat) <  90 ) hi = lat_offset_an1( nint(plat)+1 )

do jj = lo, hi
          
  dx = abs(DATA_EMIS_AN2(1,jj)-PLON)
  dy = abs(DATA_EMIS_AN2(2,jj)-PLAT)

  IF ((dx <= 0.5).and.(dy <= 0.5)) then
     IF ((DATA_EMIS_AN2(3,jj) > 0).AND.(DATA_EMIS_AN2(5,jj) > 0).AND. &
       & (DATA_EMIS_AN2(7,jj) > 0).AND.(DATA_EMIS_AN2(9,jj) > 0)) THEN
        indice = indice+1
        PBATS_VEG=DATA_EMIS_AN2(11,jj)
        IF (PBATS_VEG == 21) PBATS_VEG=14
        r1 = r1 + DATA_EMIS_AN2(3,jj)
        r2 = r2 + DATA_EMIS_AN2(5,jj)
        r3 = r3 + DATA_EMIS_AN2(7,jj)
        r4 = r4 + DATA_EMIS_AN2(9,jj)
      ENDIF
  ENDIF
enddo
if (indice > 0) then
  E_RET(1)=(r1/indice)
  E_RET(2)=(r2/indice)
  E_RET(3)=(r3/indice)
  E_RET(4)=(r4/indice)       
else
  PBATS_VEG=-9.09
  if (abs(PLAT).gt.78) THEN
    E_RET(1)=0.75
    E_RET(2)=0.75
    E_RET(3)=0.75
    E_RET(4)=0.75

  else
    E_RET(1)=0.95005
    E_RET(2)=0.95005
    E_RET(3)=0.95005
    E_RET(4)=0.95005

  end if                
endif         


END SUBROUTINE LAND_AMSUA_AN2

SUBROUTINE LAND_AMSUB_AN1 (PLAT,PLON,E_RET1,PBATS_VEG)

! ROUTINE.
! --------
!    LAND_AMSUB_AN1

! PURPOSE.
! --------
!     ROUTINE CALLED BY IFS OBSERVATION SCREENING MODULE TO 
!     INTERPOLATE SURFACE EMISSIVITY FROM AN ATLAS

! INPUT.
! ------
!     
!     PLAT - observation latitude
!     PLON - observation longitude


! OUTPUT.
! -------
!     E_RET- surface emissivity     
!     PBATS_VEG- BATS vegetation type

! F Karbou, N Bormann   09/05/07 : to caclucate land emissivities at mw frequencies   
! should be modified accoring to the emissivity atlas file format
          
USE PARKIND1  ,ONLY : JPIM     ,JPRB 

       
IMPLICIT NONE
      
REAL(KIND=JPRB), INTENT(IN)  :: PLAT
REAL(KIND=JPRB), INTENT(IN)  :: PLON

REAL(KIND=JPRB), INTENT(OUT) :: PBATS_VEG
REAL(KIND=JPRB), INTENT(OUT) :: E_RET1


INTEGER(KIND=JPIM) indice, jj
INTEGER(KIND=JPIM) lo, hi

REAL(KIND=JPRB) r1,dx,dy


r1 = 0.0
indice=0

lo = 1
hi = SIZE(DATA_EMIS_AN1(1,:)) 
if ( nint(plat) > -90 ) lo = lat_offset_an1( nint(plat)-1 )
if ( nint(plat) <  90 ) hi = lat_offset_an1( nint(plat)+1 )

do jj = lo, hi
          
  dx = abs(DATA_EMIS_AN1(1,jj)-PLON)
  dy = abs(DATA_EMIS_AN1(2,jj)-PLAT)
        
  if ((dx <= 0.5).and.(dy <= 0.5)) then
     IF (DATA_EMIS_AN1(9,jj) > 0) THEN
       indice = indice+1
       PBATS_VEG=DATA_EMIS_AN1(11,jj)
       IF (PBATS_VEG == 21) PBATS_VEG=14
       r1 = r1 + DATA_EMIS_AN1(9,jj)
     ENDIF
  endif

enddo
if (indice > 0) then
  E_RET1=(r1/indice)
else
  PBATS_VEG=-9.09
  if (abs(PLAT).gt.78) THEN
     E_RET1=0.75
  else
     E_RET1=0.95005
  end if
endif


END SUBROUTINE LAND_AMSUB_AN1

SUBROUTINE LAND_AMSUB_AN2 (PLAT,PLON,E_RET1,PBATS_VEG)

! ROUTINE.
! --------
!    LAND_AMSUA_AN2

! PURPOSE.
! --------
!     ROUTINE CALLED BY IFS OBSERVATION SCREENING MODULE TO 
!     INTERPOLATE SURFACE EMISSIVITY FROM AN ATLAS

! INPUT.
! ------
!     
!     PLAT - observation latitude
!     PLON - observation longitude


! OUTPUT.
! -------
!     E_RET- surface emissivity     
!     PBATS_VEG- BATS vegetation type

! F Karbou, N Bormann   09/05/07 : to caclucate land emissivities at mw frequencies   
! should be modified according to the emissivity atlas file format    
      
USE PARKIND1  ,ONLY : JPIM     ,JPRB 


IMPLICIT NONE
      
REAL(KIND=JPRB), INTENT(IN)  :: PLAT
REAL(KIND=JPRB), INTENT(IN)  :: PLON
REAL(KIND=JPRB), INTENT(OUT) :: PBATS_VEG
REAL(KIND=JPRB), INTENT(OUT) :: E_RET1

INTEGER(KIND=JPIM) indice, jj
INTEGER(KIND=JPIM) lo, hi

REAL(KIND=JPRB) r1,dx,dy
            

r1 = 0.0
indice=0

lo = 1
hi = SIZE(DATA_EMIS_AN2(1,:)) 
if ( nint(plat) > -90 ) lo = lat_offset_an1( nint(plat)-1 )
if ( nint(plat) <  90 ) hi = lat_offset_an1( nint(plat)+1 )

do jj = lo, hi
          
  dx = abs(DATA_EMIS_AN2(1,jj)-PLON)
  dy = abs(DATA_EMIS_AN2(2,jj)-PLAT)
        
  IF ((dx <= 0.5).and.(dy <= 0.5)) THEN
    IF (DATA_EMIS_AN2(9,jj) > 0) THEN
      indice = indice+1
      PBATS_VEG=DATA_EMIS_AN2(11,jj)
      IF (PBATS_VEG == 21) PBATS_VEG=14   
      r1 = r1 + DATA_EMIS_AN2(9,jj)                 
    ENDIF
  ENDIF
enddo
if (indice > 0) then
  E_RET1=(r1/indice)
else
  PBATS_VEG=-9.09
  if (abs(PLAT).gt.78) THEN
    E_RET1=0.75
  else
    E_RET1=0.95005
  end if
endif


END SUBROUTINE LAND_AMSUB_AN2

SUBROUTINE LAND_SSMI (PLAT,PLON,E_RET, PBATS_VEG)
 
! ROUTINE.
! --------
!    LAND_MF_AN1

! PURPOSE.
! --------
!     ROUTINE CALLED BY IFS OBSERVATION SCREENING MODULE TO 
!     INTERPOLATE SURFACE EMISSIVITY FROM AN ATLAS

! INPUT.
! ------
!     
!     PLAT - observation latitude
!     PLON - observation longitude


! OUTPUT.
! -------
!     E_RET- surface emissivity     
!     PBATS_VEG- BATS vegetation type

! F Karbou, N Bormann   09/05/07 : to caclucate land emissivities at mw frequencies
! should be modified according to the atlas file format
     
USE PARKIND1  ,ONLY : JPIM     ,JPRB 
       
IMPLICIT NONE
      
REAL(KIND=JPRB), INTENT(IN)  :: PLAT
REAL(KIND=JPRB), INTENT(IN)  :: PLON
REAL(KIND=JPRB), INTENT(OUT) :: PBATS_VEG
REAL(KIND=JPRB), INTENT(OUT) :: E_RET(7)

INTEGER(KIND=JPIM) indice, jj
INTEGER(KIND=JPIM) lo, hi

REAL(KIND=JPRB) r19h,r19v,dx,dy,r22v,r37v,r37h,r91v,r91h

r19v = 0.0
r19h = 0.0
r22v = 0.0
r37v = 0.0
r37h = 0.0
r91h = 0.0
r91v = 0.0
  
indice=0

lo = 1
hi = SIZE(ATLAS_SSMI(1,:)) 
if ( nint(plat) > -90 ) lo = lat_offset_an1( nint(plat)-1 )
if ( nint(plat) <  90 ) hi = lat_offset_an1( nint(plat)+1 )

do jj = lo, hi
          
  dx = abs(ATLAS_SSMI(1,jj)-PLON)
  dy = abs(ATLAS_SSMI(2,jj)-PLAT)
        
  if ((dx <= 0.5).and.(dy <= 0.5)) then
     indice = indice+1
     PBATS_VEG=ATLAS_SSMI(17,jj)
              
     r19v = r19v + ATLAS_SSMI(3,jj)
     r19h = r19h + ATLAS_SSMI(5,jj)
     r22v = r22v + ATLAS_SSMI(7,jj)
     r37v = r37v + ATLAS_SSMI(9,jj)
     r37h = r37h + ATLAS_SSMI(11,jj)
     r91v = r91v + ATLAS_SSMI(13,jj)
     r91h = r91h + ATLAS_SSMI(15,jj) 
                          
   endif
enddo
if (indice > 0) then
  E_RET(1)=r19v/indice
  E_RET(2)=r19h/indice
  E_RET(3)=r22v/indice
  E_RET(4)=r37v/indice
  E_RET(5)=r37h/indice
  E_RET(6)=r91v/indice
  E_RET(7)=r91h/indice     
else
  PBATS_VEG=-9.09
  E_RET(:)=-9.09       
endif


END SUBROUTINE LAND_SSMI

SUBROUTINE LAND_SSMIS (PLAT,PLON,E_RET, PBATS_VEG)
 
! ROUTINE.
! --------
!    LAND_MF_AN1

! PURPOSE.
! --------
!     ROUTINE CALLED BY IFS OBSERVATION SCREENING MODULE TO 
!     INTERPOLATE SURFACE EMISSIVITY FROM AN ATLAS

! INPUT.
! ------
!     
!     PLAT - observation latitude
!     PLON - observation longitude


! OUTPUT.
! -------
!     E_RET- surface emissivity     
!     PBATS_VEG- BATS vegetation type

! F Karbou, N Bormann   09/05/07 : to caclucate land emissivities at mw frequencies
! should be modified according to the atlas file format
   
USE PARKIND1  ,ONLY : JPIM     ,JPRB 
       
IMPLICIT NONE
      
REAL(KIND=JPRB), INTENT(IN)  :: PLAT
REAL(KIND=JPRB), INTENT(IN)  :: PLON
REAL(KIND=JPRB), INTENT(OUT) :: PBATS_VEG
REAL(KIND=JPRB), INTENT(OUT) :: E_RET(8)

INTEGER(KIND=JPIM) indice, jj
INTEGER(KIND=JPIM) lo, hi

REAL(KIND=JPRB) r19h,r19v,dx,dy,r22v,r37v,r37h,r50v,r91v,r91h


r19v = 0.0
r19h = 0.0
r22v = 0.0
r37v = 0.0
r37h = 0.0
r50v = 0.0
r91h = 0.0
r91v = 0.0
  
indice=0

lo = 1
hi = SIZE(ATLAS_SSMIS(1,:)) 
if ( nint(plat) > -90 ) lo = lat_offset_an1( nint(plat)-1 )
if ( nint(plat) <  90 ) hi = lat_offset_an1( nint(plat)+1 )

do jj = lo, hi
          
  dx = abs(ATLAS_SSMIS(1,jj)-PLON)
  dy = abs(ATLAS_SSMIS(2,jj)-PLAT)
        
  if ((dx <= 0.5).and.(dy <= 0.5)) then
     indice = indice+1
     PBATS_VEG=ATLAS_SSMIS(3,jj)
              
     r19v = r19v + ATLAS_SSMIS(4,jj)
     r19h = r19h + ATLAS_SSMIS(6,jj)
     r22v = r22v + ATLAS_SSMIS(8,jj)
     r37v = r37v + ATLAS_SSMIS(10,jj)
     r37h = r37h + ATLAS_SSMIS(12,jj)
     r50v = r50v + ATLAS_SSMIS(14,jj)
     r91v = r91v + ATLAS_SSMIS(16,jj)
     r91h = r91h + ATLAS_SSMIS(18,jj) 
                          
   endif
enddo
if (indice > 0) then
  E_RET(1)=r19v/indice
  E_RET(2)=r19h/indice
  E_RET(3)=r22v/indice
  E_RET(4)=r37v/indice
  E_RET(5)=r37h/indice
  E_RET(6)=r50v/indice
  E_RET(7)=r91v/indice
  E_RET(8)=r91h/indice     
else
  PBATS_VEG=-9.09
  E_RET(:)=-9.09       
endif


END SUBROUTINE LAND_SSMIS

  
END MODULE mod_cnrm_mw_atlas
