     MODULE RTTOV_HDF_COEFS
! Description
! Modules for the HDF5 handling
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
!    Copyright 2012, EUMETSAT, All Rights Reserved.
!
! Current Code Owner: SAF NWP
!
#include "throw.h"
      USE RTTOV_TYPES
      USE RTTOV_HDF_MOD
      USE HDF5, only : HID_T
      USE H5LT
!        
      IMPLICIT NONE
#include "rttov_errorreport.interface"

      PRIVATE
      PUBLIC :: RTTOV_HDF_COEF_WH,      RTTOV_HDF_COEF_RH
      PUBLIC :: RTTOV_HDF_SCCLDCOEF_WH, RTTOV_HDF_SCCLDCOEF_RH
      PUBLIC :: RTTOV_HDF_SCAERCOEF_WH, RTTOV_HDF_SCAERCOEF_RH
      PUBLIC :: RTTOV_HDF_PCCOEF_WH,    RTTOV_HDF_PCCOEF_RH

      CONTAINS

      SUBROUTINE RTTOV_HDF_COEF_WH(X,LUN,ERR,COMPRESS)
USE RTTOV_HDF_RTTOV_COEF_IO
USE RTTOV_HDF_RTTOV_FAST_COEF_IO
USE RTTOV_HDF_RTTOV_NLTE_COEF_IO
USE RTTOV_CONST, ONLY : SPEEDL

      TYPE(RTTOV_COEF),INTENT(IN)    ::X
      INTEGER(HID_T),INTENT(IN)      ::LUN
      INTEGER(KIND=JPIM),INTENT(OUT) ::ERR
      LOGICAL,INTENT(IN),OPTIONAL    ::COMPRESS

!
      INTEGER(HID_T) :: G_ID_SUB
!
TRY

        CALL RTTOV_HDF_RTTOV_COEF_WH( x, LUN, ERR, COMPRESS=COMPRESS )
        THROWM(ERR.NE.0,"CANNOT WRITE COEF")

        ! Write out speed of light to ensure backward compatibility
        CALL WRITE_ARRAY_HDF(LUN, 'FC_SPEEDL', &
          'speed of light - this value no longer used (see rttov_const.F90)', &
          ERR, R0=SPEEDL, UNITS='cm/s' )
        THROWM(ERR.NE.0,"CANNOT WRITE "//TRIM('FC_SPEEDL'))

        CALL H5LTSET_ATTRIBUTE_STRING_F(LUN, '.', "Description",   &
        "This is a RTTOV coefficient structure" // &
        CHAR(0), ERR )
        THROWM(ERR.NE.0,"CANNOT WRITE ATTRIBUTE")

        CALL MKPAR( LUN, "THERMAL", G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CREATE GROUP THERMAL")
      
        CALL RTTOV_HDF_RTTOV_FAST_COEF_WH( x%THERMAL, G_ID_SUB, ERR, COMPRESS=COMPRESS )
        THROWM(ERR.NE.0,"CANNOT WRITE COEF%THERMAL")

        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP THERMAL")
      
        IF( ASSOCIATED( x%SOLAR ) .AND. x%SOLARCOEF ) THEN
          CALL MKPAR( LUN, "SOLAR", G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT CREATE GROUP SOLAR")
      
          CALL RTTOV_HDF_RTTOV_FAST_COEF_WH( x%SOLAR, G_ID_SUB, ERR, COMPRESS=COMPRESS )
          THROWM(ERR.NE.0,"CANNOT WRITE COEF%SOLAR")

          CALL H5GCLOSE_F( G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT CLOSE GROUP SOLAR")
        ENDIF

        IF( x%NLTECOEF ) THEN
          CALL MKPAR( LUN, "NLTE_COEF", G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT CREATE GROUP NLTE_COEF")
      
          CALL RTTOV_HDF_RTTOV_NLTE_COEF_WH( x%NLTE_COEF, G_ID_SUB, ERR, COMPRESS=COMPRESS )
          THROWM(ERR.NE.0,"CANNOT WRITE COEF%NLTE_COEF")

          CALL H5GCLOSE_F( G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT CLOSE GROUP NLTE_COEF")
        ENDIF

CATCH
      END SUBROUTINE

      SUBROUTINE RTTOV_HDF_SCCLDCOEF_WH(X, Y, LUN,ERR,COMPRESS)
      TYPE(RTTOV_COEF_SCATT_IR), INTENT(IN)              :: X
      TYPE(RTTOV_OPTPAR_IR    ), INTENT(IN)              :: Y
      INTEGER(HID_T),INTENT(IN)      ::LUN
      INTEGER(KIND=JPIM),INTENT(OUT) ::ERR
      LOGICAL,INTENT(IN),OPTIONAL    ::COMPRESS
!
      INTEGER(HID_T) :: G_ID_SUB
!
TRY

        CALL H5LTSET_ATTRIBUTE_STRING_F(LUN, '.', "Description",   &
        "This is a RTTOV cloud coefficient structure" // &
        CHAR(0), ERR )
        THROWM(ERR.NE.0,"CANNOT WRITE ATTRIBUTE")
        
        CALL MKPAR( LUN, "WATERCLOUDS", G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CREATE GROUP WATERCLOUDS")
      
        CALL RTTOV_HDF_RTTOV_WATERCLOUDS_WH( x, y, G_ID_SUB, ERR, COMPRESS=COMPRESS )
        THROWM(ERR.NE.0,"CANNOT WRITE WATERCLOUDS COEFS")

        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP WATERCLOUDS")
      
        CALL MKPAR( LUN, "ICECLOUDS", G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CREATE GROUP ICECLOUDS")
      
        CALL RTTOV_HDF_RTTOV_ICECLOUDS_WH( x, y, G_ID_SUB, ERR, COMPRESS=COMPRESS )
        THROWM(ERR.NE.0,"CANNOT WRITE ICECLOUDS COEFS")

        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP ICECLOUDS")

CATCH
      END SUBROUTINE

      SUBROUTINE RTTOV_HDF_SCAERCOEF_WH(X, Y, LUN,ERR,COMPRESS)
USE RTTOV_CONST, ONLY: AER_NAME, NAER_MAX
USE RTTOV_UNIX_ENV, ONLY : RTTOV_UPPER_CASE

      TYPE(RTTOV_COEF_SCATT_IR), INTENT(IN)              :: X
      TYPE(RTTOV_OPTPAR_IR    ), INTENT(IN)              :: Y
      INTEGER(HID_T),INTENT(IN)      ::LUN
      INTEGER(KIND=JPIM),INTENT(OUT) ::ERR
      LOGICAL,INTENT(IN),OPTIONAL    ::COMPRESS

      CHARACTER(LEN=LENSH)  ::SNAME
      CHARACTER(LEN=LENSH)  ::GNAME 
      !CHARACTER(LEN=4)      ::VNAME(NAER_MAX)

      INTEGER(KIND=JPIM)    :: I
!
      INTEGER(HID_T) :: G_ID_SUB
!
TRY

      CALL H5LTSET_ATTRIBUTE_STRING_F(LUN, '.', "Description",   &
        & "This is a RTTOV aerosol coefficient structure" // &
        & CHAR(0), ERR )
      THROWM(ERR.NE.0,"CANNOT WRITE ATTRIBUTE")
        
      sname='FMV_AER_CHN'
      call write_array_hdf(lun,sname,&
        & 'number of channels for which optical parameters are stored',&
        & err,i0=x%FMV_AER_CHN )
      THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
      
      sname='FMV_AER_PHA_CHN'
      call write_array_hdf(lun,sname,&
        & 'Number of channels for which phase function values are stored',&
        & err,i0=x%FMV_AER_PHA_CHN)
      THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
      
      ! This value is always 0 in HDF files as we store the channels list
      sname='FMV_AER_PHA_IOFF'
      call write_array_hdf(lun,sname,&
        & 'index of channel for solar term',&
        & err,i0=0_jpim )
      THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
      
      sname='FMV_AER_COMP'
      call write_array_hdf(lun,sname,&
        & 'number of aerosols components',&
        & err,i0=x%FMV_AER_COMP )
      THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
      
      sname='FMV_AER_PH'
      call write_array_hdf(lun,sname,&
        & 'number of angles for phase function for aerosols',&
        & err,i0=x%FMV_AER_PH )
      THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
      
      if(associated(x%AER_PHA_CHANLIST))then
        sname='AER_PHA_CHANLIST'
        call write_array_hdf(lun,sname,&
          & 'The solar channel indexes',&
          & err,i1=x%aer_pha_chanlist )
        THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
      endif
      
      sname='FMV_AER_PH_VAL'
      if(associated(x%FMV_AER_PH_VAL))then
        call write_array_hdf(lun,sname,&
          & 'Angles for phase function for aerosols',&
          & err,r1=x%FMV_AER_PH_VAL , units = 'degrees')
        THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
      endif
      
      sname='FMV_AER_RH'
      if(associated(x%FMV_AER_RH))then
        call write_array_hdf(lun,sname,&
          & 'number of relative humidity for aerosols',&
          & err,i1=x%FMV_AER_RH )
        THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
      endif
      
      sname='FMV_AER_COMP_NAME'
      !CALL RTTOV_UPPER_CASE(VNAME, AER_NAME)
      call write_array_hdf(lun,sname,&
        & 'Aerosol names',&
        & err,c1=x%fmv_aer_comp_name)
      THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
      
      DO I = 1, x%FMV_AER_COMP
      
        CALL RTTOV_UPPER_CASE(GNAME, x%fmv_aer_comp_name(I))
      
        CALL MKPAR( LUN, TRIM(GNAME), G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CREATE GROUP "//TRIM(GNAME))
              
        sname='FMV_AER_RH_VAL'
        if(associated(y%optpaer(i)%FMV_AER_RH_VAL))then
          call write_array_hdf(g_id_sub,sname,&
            & 'Relative humidity for aerosols',&
            & err,r1=y%optpaer(i)%FMV_AER_RH_VAL , units = 'percent')
          THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
        endif
        
        sname='ABS'
        if(associated(y%optpaer(i)%ABS))then
          call write_array_hdf(g_id_sub,sname,&
            & 'Absorption (channels, humidity)',&
            & err,r2=y%optpaer(i)%ABS , units = 'm-1',compress=compress)
          THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
        endif
        
        sname='SCA'
        if(associated(y%optpaer(i)%SCA))then
          call write_array_hdf(g_id_sub,sname,&
            & 'Scattering (channels, humidity)',&
            & err,r2=y%optpaer(i)%SCA , units = 'm-1',compress=compress)
          THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
        endif
        
        sname='BPR'
        if(associated(y%optpaer(i)%BPR))then
          call write_array_hdf(g_id_sub,sname,&
            & 'Back scattering factor (channels, humidity or nbpr)',&
            & err,r2=y%optpaer(i)%BPR ,compress=compress)
          THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
        endif
        
        sname='PHA'
        if(associated(y%optpaer(i)%PHA))then
          call write_array_hdf(g_id_sub,sname,&
            & 'Phase function (channels, humidity, angles)',&
            & err,r3=y%optpaer(i)%PHA ,compress=compress)
          THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
        endif
        
        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP "//TRIM(GNAME))

ENDDO

CATCH
      END SUBROUTINE


      SUBROUTINE RTTOV_HDF_RTTOV_WATERCLOUDS_WH(X, Y, LUN,ERR,COMPRESS)
USE RTTOV_CONST, ONLY: WCL_NAME, NWCL_MAX
USE RTTOV_UNIX_ENV, ONLY : RTTOV_UPPER_CASE
      TYPE(RTTOV_COEF_SCATT_IR), INTENT(IN)              :: X
      TYPE(RTTOV_OPTPAR_IR    ), INTENT(IN)              :: Y
      INTEGER(HID_T),INTENT(IN)      ::LUN
      INTEGER(KIND=JPIM),INTENT(OUT) ::ERR
      LOGICAL,INTENT(IN),OPTIONAL    ::COMPRESS

      CHARACTER(LEN=LENSH)  ::SNAME
      CHARACTER(LEN=LENSH)  ::GNAME

      INTEGER(KIND=JPIM) :: I

!
      INTEGER(HID_T) :: G_ID_SUB
!
TRY

ERR=0_JPIM
sname='FMV_WCL_CHN'
call write_array_hdf(lun,sname,&
  & 'number of channels for which optical parameters are stored',&
  & err,i0=x%FMV_WCL_CHN )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='FMV_WCL_PHA_CHN'
call write_array_hdf(lun,sname,&
  & 'Number of channels for which phase function values are stored',&
  & err,i0=x%FMV_WCL_PHA_CHN )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

! This value is always 0 in HDF files as we store the channels list
sname='FMV_WCL_PHA_IOFF'
call write_array_hdf(lun,sname,&
  & 'index of channel for solar term',&
  & err,i0=0_jpim)
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='FMV_WCL_COMP'
call write_array_hdf(lun,sname,&
  & 'number of water cloud types',&
  & err,i0=x%FMV_WCL_COMP )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='FMV_WCL_PH'
call write_array_hdf(lun,sname,&
  & 'number of angles for phase function for water cloud types',&
  & err,i0=x%FMV_WCL_PH )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

if(associated(x%WCL_PHA_CHANLIST))then
  sname='WCL_PHA_CHANLIST'
  call write_array_hdf(lun,sname,&
    & 'The solar channel indexes',&
    & err,i1=x%wcl_pha_chanlist )
  THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
endif

sname='FMV_WCL_PH_VAL'
if(associated(x%FMV_WCL_PH_VAL))then
  call write_array_hdf(lun,sname,&
    & 'Angles for phase function for water cloud types',&
    & err,r1=x%FMV_WCL_PH_VAL , units = 'degrees')
  THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
endif

sname='FMV_WCL_RH'
if(associated(x%FMV_WCL_RH))then
  call write_array_hdf(lun,sname,&
    & 'number of relative humidity for water clouds',&
    & err,i1=x%FMV_WCL_RH )
  THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
endif

sname='CONFAC'
if(associated(x%CONFAC))then
  call write_array_hdf(lun,sname,&
    & 'Conversion from LWC to particule density',&
    & err,r1=x%CONFAC )
  THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
endif

sname='FMV_WCL_COMP_NAME'
!call RTTOV_UPPER_CASE(VNAME, x%fmv_wcl_comp_name)
call write_array_hdf(lun,sname,&
  & 'Cloud names',&
  & err,c1=x%fmv_wcl_comp_name)

DO I = 1, x%FMV_WCL_COMP

  CALL RTTOV_UPPER_CASE(GNAME, x%fmv_wcl_comp_name(I))

  CALL MKPAR( LUN, TRIM(GNAME), G_ID_SUB, ERR )
  THROWM(ERR.NE.0,"CANNOT CREATE GROUP "//TRIM(GNAME))

  sname='FMV_WCL_RH_VAL'
  if(associated(y%optpwcl(i)%FMV_WCL_RH_VAL))then
    call write_array_hdf(g_id_sub,sname,&
      & 'Relative humidity for water clouds',&
      & err,r1=y%optpwcl(i)%FMV_WCL_RH_VAL , units = 'percent')
    THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
  endif

  sname='ABS'
  if(associated(y%optpwcl(i)%ABS))then
    call write_array_hdf(g_id_sub,sname,&
      & 'Absorption (channels, humidity)',&
      & err,r2=y%optpwcl(i)%ABS , units = 'm-1',compress=compress)
    THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
  endif
  
  sname='SCA'
  if(associated(y%optpwcl(i)%SCA))then
    call write_array_hdf(g_id_sub,sname,&
      & 'Scattering (channels, humidity)',&
      & err,r2=y%optpwcl(i)%SCA , units = 'm-1',compress=compress)
    THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
  endif
  
  sname='BPR'
  if(associated(y%optpwcl(i)%BPR))then
    call write_array_hdf(g_id_sub,sname,&
      & 'Back scattering factor (channels, humidity)',&
      & err,r2=y%optpwcl(i)%BPR ,compress=compress)
    THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
  endif
  
  sname='PHA'
  if(associated(y%optpwcl(i)%PHA))then
    call write_array_hdf(g_id_sub,sname,&
      & 'Phase function (channels, humidity, angles)',&
      & err,r3=y%optpwcl(i)%PHA ,compress=compress)
    THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
  endif
  
  CALL H5GCLOSE_F( G_ID_SUB, ERR )
  THROWM(ERR.NE.0,"CANNOT CLOSE GROUP "//TRIM(GNAME))

ENDDO

err=0_jpim 

CATCH
      END SUBROUTINE


      SUBROUTINE RTTOV_HDF_RTTOV_ICECLOUDS_WH(X, Y, LUN,ERR,COMPRESS)
USE RTTOV_CONST, ONLY: WCL_NAME
USE RTTOV_UNIX_ENV, ONLY : RTTOV_UPPER_CASE
      TYPE(RTTOV_COEF_SCATT_IR), INTENT(IN)              :: X
      TYPE(RTTOV_OPTPAR_IR    ), INTENT(IN)              :: Y
      INTEGER(HID_T),INTENT(IN)      ::LUN
      INTEGER(KIND=JPIM),INTENT(OUT) ::ERR
      LOGICAL,INTENT(IN),OPTIONAL    ::COMPRESS

      CHARACTER(LEN=LENSH)  ::SNAME
      CHARACTER(LEN=LENSH)  ::GNAME
      CHARACTER(LEN=LENSH)  ::VNAME(x%FMV_ICL_ISHP)
      INTEGER(KIND=JPIM) :: I

!
      INTEGER(HID_T) :: G_ID_SUB
!
TRY

ERR=0_JPIM
sname='FMV_ICL_CHN'
call write_array_hdf(lun,sname,&
  & 'number of channels for which optical parameters are stored',&
  & err,i0=x%FMV_ICL_CHN )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='FMV_ICL_PHA_CHN'
call write_array_hdf(lun,sname,&
  & 'Number of channels for which phase function values are stored',&
  & err,i0=x%FMV_ICL_PHA_CHN )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

! This value is always 0 in HDF files as we store the channels list
sname='FMV_ICL_PHA_IOFF'
call write_array_hdf(lun,sname,&
  & 'index of channel for solar term ',&
  & err,i0=0_jpim )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='ICL_NABS'
call write_array_hdf(lun,sname,&
  & 'Number of regression coeffs for absorption optical depth',&
  & err,i0=x%icl_nabs )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='ICL_NSCA'
call write_array_hdf(lun,sname,&
  & 'Number of regression coeffs for scattering optical depth',&
  & err,i0=x%icl_nsca )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='ICL_NBPR'
call write_array_hdf(lun,sname,&
  & 'Number of regression coeffs for backscatter parameter',&
  & err,i0=x%icl_nbpr )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='FMV_ICL_ISHP'
call write_array_hdf(lun,sname,&
  & 'Number of shapes for which parameters are available',&
  & err,i0=x%FMV_ICL_ISHP )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='FMV_ICL_COMP'
call write_array_hdf(lun,sname,&
  & 'number effective diameter for each size distribution',&
  & err,i0=x%FMV_ICL_COMP )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='FMV_ICL_PH'
call write_array_hdf(lun,sname,&
  & 'number of angles for phase function for ice cloud types',&
  & err,i0=x%FMV_ICL_PH )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

if(associated(x%ICL_PHA_CHANLIST))then
  sname='ICL_PHA_CHANLIST'
  call write_array_hdf(lun,sname,&
    & 'The solar channel indexes',&
    & err,i1=x%icl_pha_chanlist )
  THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
endif

sname='FMV_ICL_PH_VAL'
if(associated(x%FMV_ICL_PH_VAL))then
  call write_array_hdf(lun,sname,&
    & 'Angles for phase function for ice cloud types',&
    & err,r1=x%FMV_ICL_PH_VAL , units = 'degrees')
  THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
endif

DO I = 1, x%FMV_ICL_ISHP
  IF( I .EQ. 1  ) THEN
    VNAME(I) = "HEXAGONAL"
  ELSEIF( I .EQ. 2  ) THEN
    VNAME(I) = "AGGREGATE"
  ELSE
    WRITE(VNAME(I),'(I2.2)') I
  ENDIF
ENDDO

sname='ICL_NAME'
call write_array_hdf(lun,sname,&
  & 'Ice  shape names',&
  & err,c1=VNAME)
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

DO I = 1, x%FMV_ICL_ISHP

  IF( I .EQ. 1  ) THEN
    GNAME = "HEXAGONAL"
  ELSEIF( I .EQ. 2  ) THEN
    GNAME = "AGGREGATE"
  ELSE
    WRITE(GNAME,'(I2.2)') I
  ENDIF
  
  CALL MKPAR( LUN, TRIM(GNAME), G_ID_SUB, ERR )
  THROWM(ERR.NE.0,"CANNOT CREATE GROUP "//TRIM(GNAME))
        
  sname='FMV_ICL_DG'
  if(associated(x%FMV_ICL_DG))then
    call write_array_hdf(g_id_sub,sname,&
      & 'effective diameter for each size distribution',&
      & err,r1=x%FMV_ICL_DG(:,I) , units = 'microns')
    THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
  endif

  sname='FMV_ICL_COMP_NAME'
  if(associated(x%FMV_ICL_COMP_NAME))then
    call write_array_hdf(g_id_sub,sname,&
      & 'ice particle type name',&
      & err,c1=x%FMV_ICL_COMP_NAME(:,I))
    THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
  endif

  sname='ABS'
  if(associated(y%optpicl(i)%ABS))then
    call write_array_hdf(g_id_sub,sname,&
      & 'Absorption (channels, nabs)',&
      & err,r2=y%optpicl(i)%ABS , units = 'm-1',compress=compress)
    THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
  endif
  
  sname='SCA'
  if(associated(y%optpicl(i)%SCA))then
    call write_array_hdf(g_id_sub,sname,&
      & 'Scattering (channels, nsca)',&
      & err,r2=y%optpicl(i)%SCA , units = 'm-1',compress=compress)
    THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
  endif
  
  sname='BPR'
  if(associated(y%optpicl(i)%BPR))then
    call write_array_hdf(g_id_sub,sname,&
      & 'Back scattering factor (channels, nbpr)',&
      & err,r2=y%optpicl(i)%BPR ,compress=compress)
    THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
  endif
  
  sname='PHA'
  if(associated(y%optpicl(i)%PHA))then
    call write_array_hdf(g_id_sub,sname,&
      & 'Phase function (channels, n_size_distrib_de, angles)',&
      & err,r3=y%optpicl(i)%PHA ,compress=compress)
    THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
  endif
  
  CALL H5GCLOSE_F( G_ID_SUB, ERR )
  THROWM(ERR.NE.0,"CANNOT CLOSE GROUP "//TRIM(GNAME))

ENDDO

err=0_jpim 

CATCH
      END SUBROUTINE
      
      SUBROUTINE RTTOV_HDF_PCCOEF_WH(X,LUN,ERR,COMPRESS)
USE RTTOV_HDF_RTTOV_COEF_PCC_IO
USE RTTOV_HDF_RTTOV_COEF_PCC1_IO
USE RTTOV_HDF_RTTOV_COEF_PCC2_IO

      TYPE(RTTOV_COEF_PCCOMP),INTENT(IN)    ::X
      INTEGER(HID_T),INTENT(IN)      ::LUN
      INTEGER(KIND=JPIM),INTENT(OUT) ::ERR
      LOGICAL,INTENT(IN),OPTIONAL    ::COMPRESS

      CHARACTER(LEN=LENSH)  :: GNAME, GNAME2
      INTEGER(KIND=JPIM)    :: I, J
!
      INTEGER(HID_T) :: G_ID_SUB, G_ID_SUB2, G_ID_SUB3
!
TRY

        CALL RTTOV_HDF_RTTOV_COEF_PCC_WH(X,LUN,ERR, COMPRESS=COMPRESS)
        THROWM(ERR.NE.0,"CANNOT WRITE PC COEF")
        
        CALL H5LTSET_ATTRIBUTE_STRING_F(LUN, '.', "Description",   &
        "This is a RTTOV coefficient structure PCCOMP" // &
        CHAR(0), ERR )
        THROWM(ERR.NE.0,"CANNOT WRITE ATTRIBUTE")

        ! PCREG structure
        CALL MKPAR( LUN, "PCREG", G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CREATE GROUP PCREG")
      
        DO J = 1, X%FMV_PC_BANDS
          WRITE(GNAME,'(I2.2)') J
          CALL MKPAR( G_ID_SUB, TRIM(GNAME), G_ID_SUB2, ERR )
          THROWM(ERR.NE.0,"CANNOT CREATE GROUP PCREG/"//TRIM(GNAME))

          DO I = 1, X%FMV_PC_SETS(J)

            WRITE(GNAME2,'(I2.2)') I
            CALL MKPAR( G_ID_SUB2, TRIM(GNAME2), G_ID_SUB3, ERR )
            THROWM(ERR.NE.0,"CANNOT CREATE GROUP PCREG/"//TRIM(GNAME)//'/'//TRIM(GNAME2))

            CALL RTTOV_HDF_RTTOV_COEF_PCC1_WH( X%PCREG(j,i), G_ID_SUB3, ERR , COMPRESS=COMPRESS)
            THROWM(ERR.NE.0,"CANNOT WRITE PCREG/"//TRIM(GNAME)//'/'//TRIM(GNAME2))

            CALL H5GCLOSE_F( G_ID_SUB3, ERR )
            THROWM(ERR.NE.0,"CANNOT CLOSE GROUP PCREG/"//TRIM(GNAME)//'/'//TRIM(GNAME2))

          END DO

          CALL H5GCLOSE_F( G_ID_SUB2, ERR )
          THROWM(ERR.NE.0,"CANNOT CLOSE GROUP PCREG/"//TRIM(GNAME))

        END DO

        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP PCREG")

        ! EIGEN structure
        CALL MKPAR( LUN, "EIGEN", G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CREATE GROUP EIGEN")

        DO I = 1, X%FMV_PC_BANDS

          WRITE(GNAME,'(I2.2)') I
          CALL MKPAR( G_ID_SUB, TRIM(GNAME), G_ID_SUB2, ERR )
          THROWM(ERR.NE.0,"CANNOT CREATE GROUP EIGEN/"//TRIM(GNAME))

          CALL RTTOV_HDF_RTTOV_COEF_PCC2_WH( X%EIGEN(i), G_ID_SUB2, ERR, COMPRESS=COMPRESS )
          THROWM(ERR.NE.0,"CANNOT WRITE EIGEN/"//TRIM(GNAME))

          CALL H5GCLOSE_F( G_ID_SUB2, ERR )
          THROWM(ERR.NE.0,"CANNOT CLOSE GROUP EIGEN/"//TRIM(GNAME))

        END DO

        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP EIGEN")

      
CATCH
      END SUBROUTINE
      
      SUBROUTINE RTTOV_HDF_COEF_RH(X,LUN,ERR)
      
USE RTTOV_CONST, ONLY :  &
       & gas_id_mixed,           &
       & gas_id_watervapour,     &
       & gas_id_ozone,           &
       & gas_id_wvcont,          &
       & gas_id_co2,             &
       & gas_id_n2o,             &
       & gas_id_co,              &
       & gas_id_ch4
USE RTTOV_HDF_RTTOV_COEF_IO
USE RTTOV_HDF_RTTOV_FAST_COEF_IO
USE RTTOV_HDF_RTTOV_NLTE_COEF_IO

      TYPE(RTTOV_COEF),INTENT(OUT)    ::X
      INTEGER(HID_T),INTENT(IN)       ::LUN
      INTEGER(KIND=JPIM),INTENT(OUT)  ::ERR
      
      LOGICAL :: LEXT
      INTEGER(HID_T) :: G_ID_SUB
      INTEGER(KIND=JPIM) :: n

TRY


        CALL RTTOV_HDF_RTTOV_COEF_RH( x, LUN, ERR )
        THROWM(ERR.NE.0,"CANNOT READ COEF")

        ! No channel selection with HDF5 so number of channels
        ! in file is same as number of channels extracted
        x%fmv_ori_nchn = x%fmv_chn

        DO n = 1, x%fmv_gas
          SELECT CASE (x%fmv_gas_id(n))
          CASE (gas_id_mixed)
            x%nmixed  = x%fmv_var(n)
            x%ncmixed = x%fmv_coe(n)
            x%nlevels = x%fmv_lvl(n)
          CASE (gas_id_watervapour)
            x%nwater  = x%fmv_var(n)
            x%ncwater = x%fmv_coe(n)
          CASE (gas_id_ozone)
            x%nozone  = x%fmv_var(n)
            x%ncozone = x%fmv_coe(n)
          CASE (gas_id_wvcont)
            x%nwvcont  = x%fmv_var(n)
            x%ncwvcont = x%fmv_coe(n)
          CASE (gas_id_co2)
            x%nco2  = x%fmv_var(n)
            x%ncco2 = x%fmv_coe(n)
          CASE (gas_id_n2o)
            x%nn2o  = x%fmv_var(n)
            x%ncn2o = x%fmv_coe(n)
          CASE (gas_id_co)
            x%nco  = x%fmv_var(n)
            x%ncco = x%fmv_coe(n)
          CASE (gas_id_ch4)
            x%nch4  = x%fmv_var(n)
            x%ncch4 = x%fmv_coe(n)
          END SELECT
        END DO
        x%NLAYERS = x%NLEVELS - 1
    

        NULLIFY( x%THERMAL)
        
        CALL H5GOPEN_F( LUN, "THERMAL", G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN GROUP THERMAL")
      
        ALLOCATE (x%THERMAL, STAT = ERR)
        THROWM( ERR .NE. 0, "ALLOCATION OF THERMAL FAST COEFS")
      
        CALL RTTOV_HDF_RTTOV_FAST_COEF_RH( x%THERMAL, G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT READ COEF%THERMAL")

        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP THERMAL")
        
        x%PMC_SHIFT = ASSOCIATED(x%PMC_COEF)

        x%SOLARCOEF = .FALSE.
        NULLIFY( x%SOLAR)

        CALL H5LEXISTS_F( LUN, "SOLAR", LEXT, ERR )
        THROWM(ERR.NE.0,"CANNOT TEST SOLAR ")
        IF( LEXT ) THEN

          CALL H5GOPEN_F( LUN, "SOLAR", G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT OPEN GROUP SOLAR")

          ALLOCATE (x%SOLAR, STAT = ERR)
          THROWM( ERR .NE. 0, "ALLOCATION OF SOLAR FAST COEFS")
          
          CALL RTTOV_HDF_RTTOV_FAST_COEF_RH( x%SOLAR, G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT READ COEF%SOLAR")

          CALL H5GCLOSE_F( G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT CLOSE GROUP SOLAR")
          
          x%SOLARCOEF = .TRUE.

        ENDIF


        x%NLTECOEF = .FALSE.
        NULLIFY( x%NLTE_COEF)

        CALL H5LEXISTS_F( LUN, "NLTE_COEF", LEXT, ERR )
        THROWM(ERR.NE.0,"CANNOT TEST NLTE_COEF ")
        IF( LEXT ) THEN

          CALL H5GOPEN_F( LUN, "NLTE_COEF", G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT OPEN GROUP NLTE_COEF")

          ALLOCATE (x%NLTE_COEF, STAT = ERR)
          THROWM( ERR .NE. 0, "ALLOCATION OF NLTE_COEF FAST COEFS")
          
          CALL RTTOV_HDF_RTTOV_NLTE_COEF_RH( x%NLTE_COEF, G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT READ COEF%NLTE_COEF")

          CALL H5GCLOSE_F( G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT CLOSE GROUP NLTE_COEF")
          
          x%NLTECOEF = .TRUE.

        ENDIF
          
CATCH
      END SUBROUTINE
      
      SUBROUTINE RTTOV_HDF_SCCLDCOEF_RH(X, Y, Z, LUN,ERR)
      TYPE(RTTOV_COEF_SCATT_IR), INTENT(INOUT)              :: X
      TYPE(RTTOV_OPTPAR_IR    ), INTENT(INOUT)              :: Y 
      TYPE(RTTOV_COEF    ), INTENT(IN)              :: Z

      INTEGER(HID_T),INTENT(IN)      ::LUN
      INTEGER(KIND=JPIM),INTENT(OUT) ::ERR
!
      INTEGER(HID_T) :: G_ID_SUB
!
TRY
        CALL H5GOPEN_F( LUN, "WATERCLOUDS", G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN GROUP WATERCLOUDS")
      
        CALL RTTOV_HDF_RTTOV_WATERCLOUDS_RH( x, y, z, G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT READ WATERCLOUDS COEFS")

        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP WATERCLOUDS")
      
        CALL H5GOPEN_F( LUN, "ICECLOUDS", G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN GROUP ICECLOUDS")
      
        CALL RTTOV_HDF_RTTOV_ICECLOUDS_RH( x, y, z, G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT READ ICECLOUDS COEFS")

        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP ICECLOUDS")

CATCH
      END SUBROUTINE
      
      SUBROUTINE RTTOV_HDF_RTTOV_WATERCLOUDS_RH(X, Y, Z, LUN,ERR)
USE RTTOV_CONST, ONLY: WCL_NAME, NWCL_MAX, DEG2RAD
USE RTTOV_UNIX_ENV, ONLY : RTTOV_UPPER_CASE
      TYPE(RTTOV_COEF_SCATT_IR), INTENT(INOUT)        :: X
      TYPE(RTTOV_OPTPAR_IR    ), INTENT(INOUT)        :: Y
      TYPE(RTTOV_COEF    ), INTENT(IN)                :: Z
      INTEGER(HID_T),INTENT(IN)      ::LUN
      INTEGER(KIND=JPIM),INTENT(OUT) ::ERR

      CHARACTER(LEN=LENSH)  ::SNAME
      CHARACTER(LEN=LENSH)  ::GNAME

      INTEGER(KIND=JPIM) :: I, ICOUNT, K, N
      LOGICAL  :: LEXT
!
      INTEGER(HID_T) :: G_ID_SUB
!
TRY

ERR=0_JPIM

sname='FMV_WCL_CHN'
call read_array_hdf(lun,sname,err,i0=x%FMV_WCL_CHN)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

IF (z%fmv_ori_nchn /= x%fmv_wcl_chn) THEN
    err = errorstatus_fatal
    THROWM(err.ne.0,"Incompatible channels between rtcoef and sccldcoef files")
ENDIF

! Read in variable n_phase_channels, FMV_WCL_PHA_CHN will be later affected
sname='FMV_WCL_PHA_CHN'
call read_array_hdf(lun,sname,err,i0=x%FMV_WCL_PHA_CHN)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

sname='FMV_WCL_PHA_IOFF'
call read_array_hdf(lun,sname,err,i0=x%FMV_WCL_PHA_IOFF)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

sname='FMV_WCL_COMP'
call read_array_hdf(lun,sname,err,i0=x%FMV_WCL_COMP)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

sname='FMV_WCL_PH'
call read_array_hdf(lun,sname,err,i0=x%FMV_WCL_PH)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))


! Sort out the solar channels/phase functions
IF (x%FMV_WCL_PHA_CHN > 0) THEN

    ! Always get solar channel indexes from HDF file
    sname='WCL_PHA_CHANLIST'
    call read_array_hdf(lun,sname,err,pi1=x%WCL_PHA_CHANLIST )
    THROWM(err.ne.0,"CANNOT READ "//trim(sname))


    IF (x%fmv_wcl_pha_chn > 0) THEN
        ! Copy solar channels to extract into correctly-sized array
        
        ! Create map from extracted channel list into pha array
        ALLOCATE (x%wcl_pha_index(z%fmv_chn))
        x%wcl_pha_index(:) = 0
        x%wcl_pha_index(x%wcl_pha_chanlist(1:x%fmv_wcl_pha_chn)) = &
            (/ (i, i = 1, x%fmv_wcl_pha_chn) /)
        
        ! Reset the solar channel offset if necessary
        IF (x%fmv_wcl_pha_ioff > 0) THEN
            x%fmv_wcl_pha_ioff = x%wcl_pha_chanlist(1)
        ENDIF
    ENDIF

ENDIF


sname='FMV_WCL_PH_VAL'
call read_array_hdf(lun,sname,err,pr1=x%FMV_WCL_PH_VAL)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

sname='FMV_WCL_RH'
call read_array_hdf(lun,sname,err,pi1=x%FMV_WCL_RH)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

sname='CONFAC'
call read_array_hdf(lun,sname,err,pr1=x%CONFAC)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

sname='FMV_WCL_COMP_NAME'
call read_array_hdf(lun,sname,err,pc1=x%fmv_wcl_comp_name)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

ALLOCATE (x%fmv_wcl_ph_val_cos(x%fmv_wcl_ph), STAT = ERR)
THROWM( ERR .NE. 0, "allocation of fmv_wcl_ph_val_cos")


    DO i = 1, x%fmv_wcl_ph
      x%fmv_wcl_ph_val_cos(i) = cos(x%fmv_wcl_ph_val(i) * deg2rad)
    ENDDO

    x%fmv_wcl_ph_val_min = 99999.0_JPRB

    DO i = 1, x%fmv_wcl_ph - 1

      IF (x%fmv_wcl_ph_val_min > x%fmv_wcl_ph_val(i + 1) - x%fmv_wcl_ph_val(i)) THEN
        x%fmv_wcl_ph_val_min = x%fmv_wcl_ph_val(i + 1) - x%fmv_wcl_ph_val(i)
      ENDIF

    ENDDO

    x%fmv_wcl_ph_val_min = x%fmv_wcl_ph_val_min / 2.0_JPRB
    icount = int (x%fmv_wcl_ph_val(x%fmv_wcl_ph) / x%fmv_wcl_ph_val_min, jpim)
    ALLOCATE (x%ifmv_wcl_ph_val(icount), STAT = ERR)
    THROWM( ERR .NE. 0, "allocation of ifmv_wcl_ph_val")

    k = 1

    DO i = 1, icount

      IF (x%fmv_wcl_ph_val(k) >= (i + 1) * x%fmv_wcl_ph_val_min) THEN
        x%ifmv_wcl_ph_val(i) = k
      ELSE
        k = k + 1
        x%ifmv_wcl_ph_val(i) = k
      ENDIF

    ENDDO


ALLOCATE (y%optpwcl(x%fmv_wcl_comp), STAT = ERR)
THROWM( ERR .NE. 0, "allocation of optpwcl")

DO n = 1, x%fmv_wcl_comp
  CALL rttov_hdf_nullify_coef_scatt_ir (y%optpwcl(n))
ENDDO


DO I = 1, x%FMV_WCL_COMP

  CALL RTTOV_UPPER_CASE(GNAME, x%fmv_wcl_comp_name(I))
!   write(0,*) "OPEN GROUP "//TRIM(GNAME)
  CALL H5GOPEN_F( LUN, TRIM(GNAME), G_ID_SUB, ERR )
  THROWM(ERR.NE.0,"CANNOT OPEN GROUP "//TRIM(GNAME))
  
  sname='FMV_WCL_RH_VAL'
  call read_array_hdf(G_ID_SUB,sname,err,pr1=y%optpwcl(i)%FMV_WCL_RH_VAL)
  THROWM(err.ne.0,"CANNOT READ "//trim(sname))
  
  sname='ABS'
  call read_array_hdf(G_ID_SUB,sname,err,pr2=y%optpwcl(i)%ABS)
  THROWM(err.ne.0,"CANNOT READ "//trim(sname))
  
  sname='SCA'
  call read_array_hdf(G_ID_SUB,sname,err,pr2=y%optpwcl(i)%SCA)
  THROWM(err.ne.0,"CANNOT READ "//trim(sname))
  
  sname='BPR'
  call read_array_hdf(G_ID_SUB,sname,err,pr2=y%optpwcl(i)%BPR)
  THROWM(err.ne.0,"CANNOT READ "//trim(sname))
  
  sname='PHA'
  call h5lexists_f( G_ID_SUB, sname, lext, err )
  if( lext ) then
  call read_array_hdf(G_ID_SUB,sname,err,pr3=y%optpwcl(i)%PHA)
  THROWM(err.ne.0,"CANNOT READ "//trim(sname))
  endif
  
  CALL H5GCLOSE_F( G_ID_SUB, ERR )
  THROWM(ERR.NE.0,"CANNOT CLOSE GROUP "//TRIM(GNAME))

ENDDO

err=0_jpim 

CATCH
      END SUBROUTINE 
      
      SUBROUTINE RTTOV_HDF_RTTOV_ICECLOUDS_RH(X, Y, Z,LUN,ERR)
USE RTTOV_CONST, ONLY: DEG2RAD
USE RTTOV_UNIX_ENV, ONLY : RTTOV_UPPER_CASE
      TYPE(RTTOV_COEF_SCATT_IR), INTENT(INOUT)        :: X
      TYPE(RTTOV_OPTPAR_IR    ), INTENT(INOUT)        :: Y
      TYPE(RTTOV_COEF    ), INTENT(IN)                :: Z
      INTEGER(HID_T),INTENT(IN)      ::LUN
      INTEGER(KIND=JPIM),INTENT(OUT) ::ERR

      CHARACTER(LEN=LENSH)  ::SNAME
      CHARACTER(LEN=LENSH)  ::GNAME
      CHARACTER(LEN=LENSH), POINTER  ::VNAME(:)
      INTEGER(KIND=JPIM) :: I, ICOUNT, K
      LOGICAL :: LEXT


!
      INTEGER(HID_T) :: G_ID_SUB
!
TRY

ERR=0_JPIM

sname='FMV_ICL_CHN'
call read_array_hdf(lun,sname,err,i0=x%FMV_ICL_CHN)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

IF (z%fmv_ori_nchn /= x%fmv_icl_chn) THEN
    err = errorstatus_fatal
    THROWM(err.ne.0,"Incompatible channels between rtcoef and sccldcoef files")
ENDIF

sname='FMV_ICL_PHA_CHN'
call read_array_hdf(lun,sname,err,i0=x%FMV_ICL_PHA_CHN)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

sname='FMV_ICL_PHA_IOFF'
call read_array_hdf(lun,sname,err,i0=x%FMV_ICL_PHA_IOFF)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

sname='ICL_NABS'
call read_array_hdf(lun,sname,err,i0=x%ICL_NABS)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

sname='ICL_NSCA'
call read_array_hdf(lun,sname,err,i0=x%ICL_NSCA)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

sname='ICL_NBPR'
call read_array_hdf(lun,sname,err,i0=x%ICL_NBPR)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

sname='FMV_ICL_ISHP'
call read_array_hdf(lun,sname,err,i0=x%FMV_ICL_ISHP)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

sname='FMV_ICL_COMP'
call read_array_hdf(lun,sname,err,i0=x%FMV_ICL_COMP)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

sname='FMV_ICL_PH'
call read_array_hdf(lun,sname,err,i0=x%FMV_ICL_PH)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

!sname='ICL_PHA_CHANLIST'
!call read_array_hdf(lun,sname,err,pi1=x%ICL_PHA_CHANLIST)
!THROWM(err.ne.0,"CANNOT READ "//trim(sname))

! Sort out the solar channels/phase functions
IF (x%FMV_ICL_PHA_CHN > 0) THEN

        ! Always get solar channel indexes from HDF file
    sname='ICL_PHA_CHANLIST'
    call read_array_hdf(lun,sname,err,pi1=x%ICL_PHA_CHANLIST )
    THROWM(err.ne.0,"CANNOT READ "//trim(sname))


    IF (x%fmv_icl_pha_chn > 0) THEN
        ! Copy solar channels to extract into correctly-sized array

        ! Create map from extracted channel list into pha array
        ALLOCATE (x%icl_pha_index(z%fmv_chn))
        x%icl_pha_index(:) = 0
        x%icl_pha_index(x%icl_pha_chanlist(1:x%fmv_icl_pha_chn)) = &
              (/ (i, i = 1, x%fmv_icl_pha_chn) /)

        ! Reset the solar channel offset if necessary
        IF (x%fmv_icl_pha_ioff > 0) THEN
          x%fmv_icl_pha_ioff = x%icl_pha_chanlist(1)
        ENDIF
    ENDIF

ENDIF




sname='FMV_ICL_PH_VAL'
call read_array_hdf(lun,sname,err,pr1=x%FMV_ICL_PH_VAL)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

sname='ICL_NAME'
call read_array_hdf(lun,sname,err,pc1=VNAME)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

ALLOCATE (x%fmv_icl_dg(x%fmv_icl_comp, x%fmv_icl_ishp), STAT = ERR)
THROWM( ERR .NE. 0, "allocation of fmv_icl_dg")

ALLOCATE (x%fmv_icl_comp_name(x%fmv_icl_comp, x%fmv_icl_ishp), STAT = ERR)
THROWM( ERR .NE. 0, "allocation of fmv_icl_comp_name")

ALLOCATE (x%fmv_icl_ph_val_cos(x%fmv_icl_ph), STAT = ERR)
THROWM( ERR .NE. 0, "allocation of fmv_icl_ph_val_cos")


    DO i = 1, x%fmv_icl_ph
      x%fmv_icl_ph_val_cos(i) = cos(x%fmv_icl_ph_val(i) * deg2rad)
    ENDDO

    x%fmv_icl_ph_val_min = 99999.0_JPRB

    DO i = 1, x%fmv_icl_ph - 1

      IF (x%fmv_icl_ph_val_min > x%fmv_icl_ph_val(i + 1) - x%fmv_icl_ph_val(i)) THEN
        x%fmv_icl_ph_val_min = x%fmv_icl_ph_val(i + 1) - x%fmv_icl_ph_val(i)
      ENDIF

    ENDDO

    x%fmv_icl_ph_val_min = x%fmv_icl_ph_val_min / 2.0_JPRB
    icount = int (x%fmv_icl_ph_val(x%fmv_icl_ph) / x%fmv_icl_ph_val_min, jpim)
    ALLOCATE (x%ifmv_icl_ph_val(icount), STAT = ERR)
    THROWM( ERR .NE. 0, "allocation of ifmv_icl_ph_val")

    k = 1

    DO i = 1, icount

      IF (x%fmv_icl_ph_val(k) >= (i + 1) * x%fmv_icl_ph_val_min) THEN
        x%ifmv_icl_ph_val(i) = k
      ELSE
        k = k + 1
        x%ifmv_icl_ph_val(i) = k
      ENDIF

    ENDDO
    
ALLOCATE (y%optpicl(x%fmv_icl_ishp), STAT = ERR)
THROWM( ERR .NE. 0, "allocation of optpicl")

DO i = 1, x%FMV_ICL_ISHP
  CALL rttov_hdf_nullify_coef_scatt_ir (y%optpicl(i))
ENDDO

DO I = 1, x%FMV_ICL_ISHP

  CALL RTTOV_UPPER_CASE(GNAME, VNAME(I))
        
  CALL H5GOPEN_F( LUN, TRIM(GNAME), G_ID_SUB, ERR )
  THROWM(ERR.NE.0,"CANNOT OPEN GROUP "//TRIM(GNAME))
  
  sname='FMV_ICL_DG'
  call read_array_hdf(G_ID_SUB,sname,err,r1=x%FMV_ICL_DG(:,I))
  THROWM(err.ne.0,"CANNOT READ "//trim(sname))
  
  sname='FMV_ICL_COMP_NAME'
  call read_array_hdf(G_ID_SUB,sname,err,c1=x%FMV_ICL_COMP_NAME(:,I))
  THROWM(err.ne.0,"CANNOT READ "//trim(sname))

  sname='ABS'
  call read_array_hdf(G_ID_SUB,sname,err,pr2=y%optpicl(i)%ABS)
  THROWM(err.ne.0,"CANNOT READ "//trim(sname))
  
  sname='SCA'
  call read_array_hdf(G_ID_SUB,sname,err,pr2=y%optpicl(i)%SCA)
  THROWM(err.ne.0,"CANNOT READ "//trim(sname))
  
  sname='BPR'
  call read_array_hdf(G_ID_SUB,sname,err,pr2=y%optpicl(i)%BPR)
  THROWM(err.ne.0,"CANNOT READ "//trim(sname))
  
  sname='PHA'
  call h5lexists_f( G_ID_SUB, sname, lext, err )
  if( lext ) then
  call read_array_hdf(G_ID_SUB,sname,err,pr3=y%optpicl(i)%PHA)
  THROWM(err.ne.0,"CANNOT READ "//trim(sname))
  endif
  
  CALL H5GCLOSE_F( G_ID_SUB, ERR )
  THROWM(ERR.NE.0,"CANNOT CLOSE GROUP "//TRIM(GNAME))

ENDDO

err=0_jpim 

CATCH
      END SUBROUTINE
      
      SUBROUTINE RTTOV_HDF_SCAERCOEF_RH(X, Y, Z, LUN,ERR)
USE RTTOV_CONST, ONLY: AER_NAME, NAER_MAX, DEG2RAD
USE RTTOV_UNIX_ENV, ONLY : RTTOV_UPPER_CASE

      TYPE(RTTOV_COEF_SCATT_IR), INTENT(INOUT)              :: X
      TYPE(RTTOV_OPTPAR_IR    ), INTENT(INOUT)              :: Y
      TYPE(RTTOV_COEF    ), INTENT(IN)                      :: Z
      INTEGER(HID_T),INTENT(IN)      ::LUN
      INTEGER(KIND=JPIM),INTENT(OUT) ::ERR

      CHARACTER(LEN=LENSH)  ::SNAME
      CHARACTER(LEN=LENSH)  ::GNAME 
      !CHARACTER(LEN=4), POINTER     ::VNAME(:)

      INTEGER(KIND=JPIM)    :: I, ICOUNT, K
      LOGICAL :: LEXT
!
      INTEGER(HID_T) :: G_ID_SUB
!
TRY
  
sname='FMV_AER_CHN'
call read_array_hdf(lun,sname,err,i0=x%FMV_AER_CHN)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

IF (z%fmv_ori_nchn /= x%fmv_aer_chn) THEN
    err = errorstatus_fatal
    THROWM(err.ne.0,"Incompatible channels between rtcoef and scaercoef files")
ENDIF

sname='FMV_AER_PHA_CHN'
call read_array_hdf(lun,sname,err,i0=x%FMV_AER_PHA_CHN)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

sname='FMV_AER_PHA_IOFF'
call read_array_hdf(lun,sname,err,i0=x%FMV_AER_PHA_IOFF)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

sname='FMV_AER_COMP'
call read_array_hdf(lun,sname,err,i0=x%FMV_AER_COMP)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

sname='FMV_AER_PH'
call read_array_hdf(lun,sname,err,i0=x%FMV_AER_PH)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

! Sort out the solar channels/phase functions
IF (x%FMV_AER_PHA_CHN > 0) THEN

        ! Always get solar channel indexes from HDF file
    sname='AER_PHA_CHANLIST'
    call read_array_hdf(lun,sname,err,pi1=x%aer_PHA_CHANLIST )
    THROWM(err.ne.0,"CANNOT READ "//trim(sname))


    IF (x%fmv_aer_pha_chn > 0) THEN
        ! Copy solar channels to extract into correctly-sized array

        ! Create map from extracted channel list into pha array
        ALLOCATE (x%aer_pha_index(z%fmv_chn))
        x%aer_pha_index(:) = 0
        x%aer_pha_index(x%aer_pha_chanlist(1:x%fmv_aer_pha_chn)) = &
              (/ (i, i = 1, x%fmv_aer_pha_chn) /)

        ! Reset the solar channel offset if necessary
        IF (x%fmv_aer_pha_ioff > 0) THEN
          x%fmv_aer_pha_ioff = x%aer_pha_chanlist(1)
        ENDIF
    ENDIF

ENDIF

sname='FMV_AER_PH_VAL'
call read_array_hdf(lun,sname,err,pr1=x%FMV_AER_PH_VAL)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

sname='FMV_AER_RH'
call read_array_hdf(lun,sname,err,pi1=x%FMV_AER_RH)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

ALLOCATE (x%fmv_aer_ph_val_cos(x%fmv_aer_ph), STAT = ERR)
THROWM( ERR .NE. 0, "allocation of fmv_aer_ph_val_cos")


    DO i = 1, x%fmv_aer_ph
      x%fmv_aer_ph_val_cos(i) = cos(x%fmv_aer_ph_val(i) * deg2rad)
    ENDDO

    x%fmv_aer_ph_val_min = 99999.0_JPRB

    DO i = 1, x%fmv_aer_ph - 1

      IF (x%fmv_aer_ph_val_min > x%fmv_aer_ph_val(i + 1) - x%fmv_aer_ph_val(i)) THEN
        x%fmv_aer_ph_val_min = x%fmv_aer_ph_val(i + 1) - x%fmv_aer_ph_val(i)
      ENDIF

    ENDDO

    x%fmv_aer_ph_val_min = x%fmv_aer_ph_val_min / 2.0_JPRB
    icount = int (x%fmv_aer_ph_val(x%fmv_aer_ph) / x%fmv_aer_ph_val_min, jpim)
    ALLOCATE (x%ifmv_aer_ph_val(icount), STAT = ERR)
    THROWM( ERR .NE. 0, "allocation of ifmv_aer_ph_val")

    k = 1

    DO i = 1, icount

      IF (x%fmv_aer_ph_val(k) >= (i + 1) * x%fmv_aer_ph_val_min) THEN
        x%ifmv_aer_ph_val(i) = k
      ELSE
        k = k + 1
        x%ifmv_aer_ph_val(i) = k
      ENDIF

    ENDDO
    
ALLOCATE (y%optpaer(x%fmv_aer_comp), STAT = ERR)
THROWM( ERR .NE. 0, "allocation of optaer")

DO i = 1, x%FMV_AER_COMP
  CALL rttov_hdf_nullify_coef_scatt_ir (y%optpaer(i))
ENDDO

sname='FMV_AER_COMP_NAME'
call read_array_hdf(lun,sname,err,pc1=x%FMV_AER_COMP_NAME)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

DO I = 1, x%FMV_AER_COMP

  CALL RTTOV_UPPER_CASE(GNAME, x%FMV_AER_COMP_NAME(I))
        
  CALL H5GOPEN_F( LUN, TRIM(GNAME), G_ID_SUB, ERR )
  THROWM(ERR.NE.0,"CANNOT OPEN GROUP "//TRIM(GNAME))
  
  sname='FMV_AER_RH_VAL'
  call read_array_hdf(G_ID_SUB,sname,err,pr1=y%optpaer(i)%FMV_AER_RH_VAL)
  THROWM(err.ne.0,"CANNOT READ "//trim(sname))
  
  sname='ABS'
  call read_array_hdf(G_ID_SUB,sname,err,pr2=y%optpaer(i)%ABS)
  THROWM(err.ne.0,"CANNOT READ "//trim(sname))
  
  sname='SCA'
  call read_array_hdf(G_ID_SUB,sname,err,pr2=y%optpaer(i)%SCA)
  THROWM(err.ne.0,"CANNOT READ "//trim(sname))
  
  sname='BPR'
  call read_array_hdf(G_ID_SUB,sname,err,pr2=y%optpaer(i)%BPR)
  THROWM(err.ne.0,"CANNOT READ "//trim(sname))
  
  sname='PHA'
  call h5lexists_f( G_ID_SUB, sname, lext, err )
  if( lext ) then
  call read_array_hdf(G_ID_SUB,sname,err,pr3=y%optpaer(i)%PHA)
  THROWM(err.ne.0,"CANNOT READ "//trim(sname))
  endif
  
  CALL H5GCLOSE_F( G_ID_SUB, ERR )
  THROWM(ERR.NE.0,"CANNOT CLOSE GROUP "//TRIM(GNAME))
ENDDO

CATCH
      END SUBROUTINE
      
      SUBROUTINE RTTOV_HDF_PCCOEF_RH(X, Y,LUN,ERR)
USE RTTOV_CONST, ONLY: ERRORSTATUS_FATAL, SENSOR_ID_HI
USE RTTOV_HDF_RTTOV_COEF_PCC_IO
USE RTTOV_HDF_RTTOV_COEF_PCC1_IO
USE RTTOV_HDF_RTTOV_COEF_PCC2_IO

      TYPE(RTTOV_COEF_PCCOMP),INTENT(INOUT)  ::X
      TYPE(RTTOV_COEF),INTENT(IN)            ::Y
      INTEGER(HID_T),INTENT(IN)      ::LUN
      INTEGER(KIND=JPIM),INTENT(OUT) ::ERR

      CHARACTER(LEN=LENSH)  :: GNAME, GNAME2
      INTEGER(KIND=JPIM)    :: I, J
!
      INTEGER(HID_T) :: G_ID_SUB, G_ID_SUB2, G_ID_SUB3
!
TRY

        CALL RTTOV_HDF_RTTOV_COEF_PCC_RH(X,LUN,ERR)
        THROWM(ERR.NE.0,"CANNOT READ COEF")

         IF ((y%id_sensor == sensor_id_hi)) THEN
          IF (y%id_comp_pc /= x%fmv_pc_comp_pc) THEN
            err = errorstatus_fatal
            THROWM( ERR .NE. 0, "Version of PC coef file is incompatible with RTTOV regression file")
          ENDIF        
        ENDIF        

        !PCREG structure
        CALL H5GOPEN_F( LUN, "PCREG", G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN GROUP "//TRIM(GNAME))
        
        ALLOCATE (X%PCREG(X%FMV_PC_BANDS,X%FMV_PC_MSETS), STAT = ERR)
        THROWM( ERR .NE. 0, "allocation of coef_pccomp %pcreg")
  
        DO J = 1, X%FMV_PC_BANDS

          WRITE(GNAME,'(I2.2)') J
          CALL H5GOPEN_F( G_ID_SUB, TRIM(GNAME), G_ID_SUB2, ERR )
          THROWM(ERR.NE.0,"CANNOT OPEN GROUP PCREG/"//TRIM(GNAME))

          DO I = 1, X%FMV_PC_SETS(J)

            WRITE(GNAME2,'(I2.2)') I
            CALL H5GOPEN_F( G_ID_SUB2, TRIM(GNAME2), G_ID_SUB3, ERR )
            THROWM(ERR.NE.0,"CANNOT OPEN GROUP PCREG/"//TRIM(GNAME)//'/'//TRIM(GNAME2))

            CALL RTTOV_HDF_RTTOV_COEF_PCC1_INIT(X%PCREG(j,i))

            CALL RTTOV_HDF_RTTOV_COEF_PCC1_RH( X%PCREG(j,i), G_ID_SUB3, ERR )
            THROWM(ERR.NE.0,"CANNOT READ PCREG/"//TRIM(GNAME)//'/'//TRIM(GNAME2))

            CALL H5GCLOSE_F( G_ID_SUB3, ERR )
            THROWM(ERR.NE.0,"CANNOT CLOSE GROUP PCREG/"//TRIM(GNAME)//'/'//TRIM(GNAME2))

          END DO

          CALL H5GCLOSE_F( G_ID_SUB2, ERR )
          THROWM(ERR.NE.0,"CANNOT CLOSE GROUP PCREG/"//TRIM(GNAME))

        END DO

        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP PCREG")

        !EIGEN structure
        CALL H5GOPEN_F( LUN, "EIGEN", G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN GROUP "//TRIM(GNAME))

        ALLOCATE (X%EIGEN(X%FMV_PC_BANDS), STAT = ERR)
        THROWM( ERR .NE. 0, "allocation of coef_pccomp %eigen")


        DO I = 1, X%FMV_PC_BANDS

          WRITE(GNAME,'(I2.2)') I
          CALL H5GOPEN_F( G_ID_SUB, TRIM(GNAME), G_ID_SUB2, ERR )
          THROWM(ERR.NE.0,"CANNOT OPEN GROUP EIGEN/"//TRIM(GNAME))

          CALL RTTOV_HDF_RTTOV_COEF_PCC2_INIT(X%EIGEN(i))

          CALL RTTOV_HDF_RTTOV_COEF_PCC2_RH( X%EIGEN(i), G_ID_SUB2, ERR )
          THROWM(ERR.NE.0,"CANNOT READ EIGEN/"//TRIM(GNAME))

          CALL H5GCLOSE_F( G_ID_SUB2, ERR )
          THROWM(ERR.NE.0,"CANNOT CLOSE GROUP EIGEN/"//TRIM(GNAME))

        END DO

        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP EIGEN")
      
CATCH
      END SUBROUTINE
      
      SUBROUTINE rttov_hdf_nullify_coef_scatt_ir(coef_scatt_ir)
! Description:
!   Nullify the IR scattering coefficient structure
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
  USE parkind1, ONLY : jpim, jprb
  USE rttov_types, ONLY : rttov_coef_scatt_ir
  IMPLICIT NONE
  TYPE(rttov_coef_scatt_ir), INTENT(INOUT) :: coef_scatt_ir
  coef_scatt_ir%fmv_aer_chn      = 0_JPIM
  coef_scatt_ir%fmv_wcl_chn      = 0_JPIM
  coef_scatt_ir%fmv_icl_chn      = 0_JPIM
  coef_scatt_ir%fmv_aer_pha_chn  = 0_JPIM
  coef_scatt_ir%fmv_wcl_pha_chn  = 0_JPIM
  coef_scatt_ir%fmv_icl_pha_chn  = 0_JPIM
  coef_scatt_ir%fmv_aer_comp     = 0_JPIM
  coef_scatt_ir%fmv_wcl_comp     = 0_JPIM
  coef_scatt_ir%fmv_icl_comp     = 0_JPIM
  coef_scatt_ir%fmv_icl_ishp     = 0_JPIM
  coef_scatt_ir%fmv_aer_pha_ioff = 0_JPIM
  coef_scatt_ir%fmv_wcl_pha_ioff = 0_JPIM
  coef_scatt_ir%fmv_icl_pha_ioff = 0_JPIM
  coef_scatt_ir%fmv_aer_ph       = 0_JPIM
  coef_scatt_ir%fmv_wcl_ph       = 0_JPIM
  coef_scatt_ir%fmv_icl_ph       = 0_JPIM
  coef_scatt_ir%icl_nabs         = 0_JPIM
  coef_scatt_ir%icl_nsca         = 0_JPIM
  coef_scatt_ir%icl_nbpr         = 0_JPIM
  NULLIFY (coef_scatt_ir%fmv_aer_comp_name)
  NULLIFY (coef_scatt_ir%fmv_wcl_comp_name)
  NULLIFY (coef_scatt_ir%fmv_icl_comp_name)
  NULLIFY (coef_scatt_ir%fmv_aer_rh)
  NULLIFY (coef_scatt_ir%fmv_wcl_rh)
  NULLIFY (coef_scatt_ir%fmv_aer_rh_val)
  NULLIFY (coef_scatt_ir%fmv_wcl_rh_val)
  NULLIFY (coef_scatt_ir%fmv_wcl_ph_val_cos)
  coef_scatt_ir%fmv_wcl_ph_val_min = 0._jprb
  NULLIFY (coef_scatt_ir%ifmv_wcl_ph_val)
  NULLIFY (coef_scatt_ir%fmv_aer_ph_val)
  NULLIFY (coef_scatt_ir%fmv_aer_ph_val_cos)
  NULLIFY (coef_scatt_ir%ifmv_aer_ph_val)
  coef_scatt_ir%fmv_aer_ph_val_min = 0._jprb
  NULLIFY (coef_scatt_ir%fmv_wcl_ph_val)
  NULLIFY (coef_scatt_ir%fmv_icl_ph_val)
  NULLIFY (coef_scatt_ir%fmv_icl_ph_val_cos)
  coef_scatt_ir%fmv_icl_ph_val_min = 0._jprb
  NULLIFY (coef_scatt_ir%ifmv_icl_ph_val)
  NULLIFY (coef_scatt_ir%fmv_icl_dg)
  NULLIFY (coef_scatt_ir%aer_pha_chanlist)
  NULLIFY (coef_scatt_ir%wcl_pha_chanlist)
  NULLIFY (coef_scatt_ir%icl_pha_chanlist)
  NULLIFY (coef_scatt_ir%aer_pha_index)
  NULLIFY (coef_scatt_ir%wcl_pha_index)
  NULLIFY (coef_scatt_ir%icl_pha_index)
  NULLIFY (coef_scatt_ir%abs)
  NULLIFY (coef_scatt_ir%sca)
  NULLIFY (coef_scatt_ir%bpr)
  NULLIFY (coef_scatt_ir%pha)
  NULLIFY (coef_scatt_ir%confac)
END SUBROUTINE 

END MODULE
