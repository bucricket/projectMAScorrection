      SUBROUTINE RTTOV_HDF_SAVE( ERR, F, PATH, CREATE, PROFILE, PROFILES, KMATRIX, OPTIONS, CHANPROF, &
         &  EMISSIVITY, REFLECTANCE, RADIANCE, RADIANCE2, PCCOMP, TRANSMISSION, COEF, SCCLDCOEF, SCAERCOEF, &
         &  OPTP, PCCOEF, &
         &  OPT_PARAM, OPTS, I0, R0, RM0, C0, L1, &
         &  IT1, IT2, IT3, IS1, IS2, IS3, I1, I2, I3, IL1,&
         &  R1, R2, R3, R4, RM1, RM2, RM3, RM4, &
         &  SNAME, COMMENT, UNITS, COMPRESS, FORCE_DOUBLE)
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
      USE PARKIND1
      USE RTTOV_TYPES
!INTF_OFF
      USE RTTOV_HDF_MOD
      USE HDF5
      USE H5LT
! 
      USE RTTOV_HDF_PROFILES
      USE RTTOV_HDF_COEFS

      USE RTTOV_HDF_PROFILE_IO
      USE RTTOV_HDF_S2M_IO
      USE RTTOV_HDF_SSKIN_IO
      USE RTTOV_HDF_OPTIONS_IO
      USE RTTOV_HDF_OPTIONS_CONFIG_IO
      USE RTTOV_HDF_OPTIONS_RT_ALL_IO
      USE RTTOV_HDF_OPTIONS_RT_IR_IO
      USE RTTOV_HDF_OPTIONS_RT_MW_IO
      USE RTTOV_HDF_OPTIONS_PC_IO
      USE RTTOV_HDF_OPTIONS_INTERP_IO
      USE RTTOV_HDF_PCCOMP_IO
      USE RTTOV_HDF_CHANPROF_IO
      USE RTTOV_HDF_EMISSIVITY_IO
      USE RTTOV_HDF_REFLECTANCE_IO
      USE RTTOV_HDF_RTTOV_COEF_IO
      USE RTTOV_HDF_RTTOV_FAST_COEF_IO   
      USE RTTOV_HDF_RADIANCE_IO
      USE RTTOV_HDF_RADIANCE2_IO
      USE RTTOV_HDF_TRANSMISSION_IO
      USE RTTOV_HDF_OPT_PARAM_IO
!INTF_ON
!
      IMPLICIT NONE
      CHARACTER(LEN=*),    INTENT(IN)  :: F
      CHARACTER(LEN=*),    INTENT(IN)  :: PATH
      INTEGER(KIND=JPIM),  INTENT(OUT) :: ERR
      LOGICAL,             INTENT(IN), OPTIONAL  :: CREATE
      TYPE(PROFILE_TYPE), INTENT(IN), OPTIONAL   :: PROFILE
      TYPE(PROFILE_TYPE), INTENT(IN), OPTIONAL   :: PROFILES(:)
      TYPE(PROFILE_TYPE), INTENT(IN), OPTIONAL   :: KMATRIX(:)
      TYPE(RTTOV_OPTIONS), INTENT(IN), OPTIONAL  :: OPTIONS
      TYPE(RTTOV_OPTIONS), INTENT(IN), OPTIONAL  :: OPTS             ! no output
      TYPE(RTTOV_CHANPROF),    INTENT(IN), OPTIONAL :: CHANPROF(:)
      TYPE(RTTOV_EMISSIVITY),  INTENT(IN), OPTIONAL :: EMISSIVITY(:)
      TYPE(RTTOV_REFLECTANCE), INTENT(IN), OPTIONAL :: REFLECTANCE(:)
      TYPE(RADIANCE_TYPE),     INTENT(IN), OPTIONAL :: RADIANCE
      TYPE(RADIANCE2_TYPE),    INTENT(IN), OPTIONAL :: RADIANCE2
      TYPE(TRANSMISSION_TYPE), INTENT(IN), OPTIONAL :: TRANSMISSION
      TYPE(RTTOV_PCCOMP),      INTENT(IN), OPTIONAL :: PCCOMP
      TYPE(RTTOV_COEF),        INTENT(IN), OPTIONAL :: COEF
      TYPE(RTTOV_COEF_SCATT_IR), INTENT(IN), OPTIONAL  :: SCCLDCOEF
      TYPE(RTTOV_COEF_SCATT_IR), INTENT(IN), OPTIONAL  :: SCAERCOEF
      TYPE(RTTOV_OPTPAR_IR),     INTENT(IN), OPTIONAL  :: OPTP
      TYPE(RTTOV_COEF_PCCOMP),   INTENT(IN), OPTIONAL  :: PCCOEF
      TYPE(RTTOV_OPT_PARAM),     INTENT(IN), OPTIONAL  :: OPT_PARAM

      INTEGER(KIND=JPIM),      INTENT(IN), OPTIONAL :: I0
      REAL(KIND=JPRB),         INTENT(IN), OPTIONAL :: R0
      REAL(KIND=JPRM),         INTENT(IN), OPTIONAL :: RM0
      CHARACTER(LEN=*),        INTENT(IN), OPTIONAL :: C0
      REAL(KIND=JPRB),         INTENT(IN), OPTIONAL :: R1(:)
      REAL(KIND=JPRB),         INTENT(IN), OPTIONAL :: R2(:,:)
      REAL(KIND=JPRB),         INTENT(IN), OPTIONAL :: R3(:,:,:)
      REAL(KIND=JPRB),         INTENT(IN), OPTIONAL :: R4(:,:,:,:)
      REAL(KIND=JPRM),         INTENT(IN), OPTIONAL :: RM1(:)
      REAL(KIND=JPRM),         INTENT(IN), OPTIONAL :: RM2(:,:)
      REAL(KIND=JPRM),         INTENT(IN), OPTIONAL :: RM3(:,:,:)
      REAL(KIND=JPRM),         INTENT(IN), OPTIONAL :: RM4(:,:,:,:)
      INTEGER(KIND=JPIT),      INTENT(IN), OPTIONAL :: IT1(:)
      INTEGER(KIND=JPIT),      INTENT(IN), OPTIONAL :: IT2(:,:)
      INTEGER(KIND=JPIT),      INTENT(IN), OPTIONAL :: IT3(:,:,:)
      INTEGER(KIND=JPIS),      INTENT(IN), OPTIONAL :: IS1(:)
      INTEGER(KIND=JPIS),      INTENT(IN), OPTIONAL :: IS2(:,:)
      INTEGER(KIND=JPIS),      INTENT(IN), OPTIONAL :: IS3(:,:,:)
      INTEGER(KIND=JPIM),      INTENT(IN), OPTIONAL :: I1(:)
      INTEGER(KIND=JPIM),      INTENT(IN), OPTIONAL :: I2(:,:)
      INTEGER(KIND=JPIM),      INTENT(IN), OPTIONAL :: I3(:,:,:)
      INTEGER(KIND=JPIB),      INTENT(IN), OPTIONAL :: IL1(:)
      LOGICAL,                 INTENT(IN), OPTIONAL :: L1(:)
      CHARACTER(LEN=*),        INTENT(IN), OPTIONAL :: SNAME
      CHARACTER(LEN=*),        INTENT(IN), OPTIONAL :: COMMENT
      CHARACTER(LEN=*),        INTENT(IN), OPTIONAL :: UNITS
      LOGICAL,                 INTENT(IN), OPTIONAL :: COMPRESS
      LOGICAL,                 INTENT(IN), OPTIONAL :: FORCE_DOUBLE


!INTF_END
#include "rttov_errorreport.interface"

      INTEGER(KIND=JPIM) :: I
      CHARACTER(LEN=4)   :: CHPROF
      CHARACTER(LEN=128) :: LCOMMENT
      LOGICAL            :: LCREATE
      LOGICAL            :: PATHROOT
      INTEGER(HID_T)     :: FILE_ID, G_ID, G_ID_SUB
!
TRY

      LCREATE = .FALSE.
      IF( PRESENT (CREATE) ) THEN
        LCREATE = CREATE
      ENDIF
      
      LCOMMENT=""
      IF( PRESENT (COMMENT) ) THEN
        LCOMMENT = COMMENT
      ELSE
        IF( PRESENT (SNAME) ) THEN
          LCOMMENT = SNAME
        ENDIF
      ENDIF
      
      PATHROOT = TRIM(PATH) .EQ. "/"
      
      IF( LCREATE ) THEN
        CALL RTTOV_HDF_OPEN_TRUNC( ERR, F, FILE_ID)
      ELSE
        CALL RTTOV_HDF_OPEN_RDWR( ERR, F, FILE_ID)      
      ENDIF
      THROWM(ERR.NE.0,"CANNOT OPEN FILE "//TRIM(F))
      
      IF( PATHROOT ) THEN
        G_ID=FILE_ID
      ELSE
        CALL MKPAR( FILE_ID, PATH, G_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CREATE GROUP "//TRIM(PATH)//" IN FILE "//TRIM(F))
      ENDIF
      
      IF( PRESENT(PROFILE) ) THEN
      
        CALL RTTOV_HDF_one_profile_WH( PROFILE, G_ID, ERR, COMPRESS=COMPRESS )
        THROWM(ERR.NE.0,"CANNOT WRITE PROFILE IN FILE "//TRIM(F))
           
      ENDIF
      
      IF( PRESENT(PROFILES) ) THEN
      
        DO I = 1, SIZE(PROFILES)
        
          WRITE(CHPROF,"(I4.4)") I
          CALL MKPAR( G_ID, CHPROF, G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT CREATE GROUP FOR PROFILE "//CHPROF)
      
          CALL RTTOV_HDF_one_profile_WH( PROFILES(I), G_ID_SUB, ERR, COMPRESS=COMPRESS )
          THROWM(ERR.NE.0,"CANNOT WRITE PROFILE"//CHPROF//" IN FILE "//TRIM(F))
          
          CALL H5GCLOSE_F( G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT CLOSE GROUP FOR PROFILE "//CHPROF)

        END DO
        
      ENDIF
    
      IF( PRESENT(KMATRIX) .AND. PRESENT(OPTS) ) THEN
      
        CALL RTTOV_HDF_KMATRIX_WH( KMATRIX, OPTS, G_ID, ERR, COMPRESS=COMPRESS )
        THROWM(ERR.NE.0,"CANNOT WRITE PROFILE IN FILE "//TRIM(F))
           
      ENDIF
      
      IF( PRESENT(OPTIONS) ) THEN
      
        CALL RTTOV_HDF_OPTIONS_WH( OPTIONS, G_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT WRITE OPTIONS IN FILE "//TRIM(F))

        CALL RTTOV_HDF_OPTIONS_CONFIG_WH( OPTIONS%CONFIG, G_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT WRITE OPTIONS%CONFIG IN FILE "//TRIM(F))

        CALL RTTOV_HDF_OPTIONS_RT_ALL_WH( OPTIONS%RT_ALL, G_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT WRITE OPTIONS%RT_ALL IN FILE "//TRIM(F))

        CALL RTTOV_HDF_OPTIONS_RT_IR_WH( OPTIONS%RT_IR, G_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT WRITE OPTIONS%RT_IR IN FILE "//TRIM(F))

        CALL RTTOV_HDF_OPTIONS_RT_MW_WH( OPTIONS%RT_MW, G_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT WRITE OPTIONS%RT_MW IN FILE "//TRIM(F))

        CALL RTTOV_HDF_OPTIONS_PC_WH( OPTIONS%RT_IR%PC, G_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT WRITE OPTIONS%RT_IR%PC IN FILE "//TRIM(F))

        CALL RTTOV_HDF_OPTIONS_INTERP_WH( OPTIONS%INTERPOLATION, G_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT WRITE OPTIONS%INTERPOLATION IN FILE "//TRIM(F))

      ENDIF
      
      IF( PRESENT(CHANPROF) ) THEN
      
        CALL RTTOV_HDF_CHANPROF_WH(CHANPROF,G_ID,ERR, COMPRESS=COMPRESS)
        THROWM(ERR.NE.0,"CANNOT WRITE CHANPROF IN FILE "//TRIM(F))

      ENDIF
      
      IF( PRESENT(EMISSIVITY) ) THEN
      
        CALL RTTOV_HDF_EMISSIVITY_WH(EMISSIVITY,G_ID,ERR, COMPRESS=COMPRESS)
        THROWM(ERR.NE.0,"CANNOT WRITE EMISSIVITY IN FILE "//TRIM(F))

      ENDIF
      
      IF( PRESENT(REFLECTANCE) ) THEN
      
        CALL RTTOV_HDF_REFLECTANCE_WH(REFLECTANCE,G_ID,ERR, COMPRESS=COMPRESS)
        THROWM(ERR.NE.0,"CANNOT WRITE REFLECTANCE IN FILE "//TRIM(F))

      ENDIF
      
      IF( PRESENT(RADIANCE) ) THEN
      
        CALL RTTOV_HDF_RADIANCE_WH(RADIANCE,G_ID,ERR, COMPRESS=COMPRESS)
        THROWM(ERR.NE.0,"CANNOT WRITE RADIANCE IN FILE "//TRIM(F))

      ENDIF
            
      IF( PRESENT(RADIANCE2) ) THEN

        CALL RTTOV_HDF_RADIANCE2_WH(RADIANCE2,G_ID,ERR, COMPRESS=COMPRESS)
        THROWM(ERR.NE.0,"CANNOT WRITE RADIANCE2 IN FILE "//TRIM(F))

      ENDIF
      
      IF( PRESENT(PCCOMP) ) THEN
      
        CALL RTTOV_HDF_PCCOMP_WH(PCCOMP,G_ID,ERR, COMPRESS=COMPRESS)
        THROWM(ERR.NE.0,"CANNOT WRITE PCCOMP IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(TRANSMISSION) ) THEN
      
        CALL RTTOV_HDF_TRANSMISSION_WH(TRANSMISSION,G_ID,ERR, COMPRESS=COMPRESS)
        THROWM(ERR.NE.0,"CANNOT WRITE TRANSMISSION IN FILE "//TRIM(F))

      ENDIF
      
      IF( PRESENT(COEF) ) THEN
      
        CALL RTTOV_HDF_COEF_WH(COEF,G_ID,ERR, COMPRESS=COMPRESS)
        THROWM(ERR.NE.0,"CANNOT WRITE COEF IN FILE "//TRIM(F))

      ENDIF
         
      IF( PRESENT(SCCLDCOEF) .AND. PRESENT(OPTP)) THEN
      
        CALL RTTOV_HDF_SCCLDCOEF_WH(SCCLDCOEF, OPTP, G_ID,ERR, COMPRESS=COMPRESS)
        THROWM(ERR.NE.0,"CANNOT WRITE CLOUDS COEF IN FILE "//TRIM(F))

      ENDIF
      
      IF( PRESENT(SCAERCOEF) .AND. PRESENT(OPTP)) THEN
      
        CALL RTTOV_HDF_SCAERCOEF_WH(SCAERCOEF, OPTP, G_ID,ERR, COMPRESS=COMPRESS)
        THROWM(ERR.NE.0,"CANNOT WRITE AEROSOL COEF IN FILE "//TRIM(F))

      ENDIF
      
      IF( PRESENT(PCCOEF) ) THEN
      
        CALL RTTOV_HDF_PCCOEF_WH(PCCOEF, G_ID,ERR, COMPRESS=COMPRESS)
        THROWM(ERR.NE.0,"CANNOT WRITE PC COEFS COEF IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(OPT_PARAM) ) THEN
      
        CALL RTTOV_HDF_OPT_PARAM_WH(OPT_PARAM, G_ID,ERR, COMPRESS=COMPRESS)
        THROWM(ERR.NE.0,"CANNOT WRITE OPT_PARAM IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(R1) .AND. PRESENT(SNAME) ) THEN
      
        call write_array_hdf(G_ID,SNAME,LCOMMENT,err,units=units,r1=r1, FORCE_DOUBLE=FORCE_DOUBLE)
        THROWM(ERR.NE.0,"CANNOT WRITE REAL ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(RM1) .AND. PRESENT(SNAME) ) THEN
      
        call write_array_hdf(G_ID,SNAME,LCOMMENT,err,units=units,rm1=rm1, FORCE_DOUBLE=FORCE_DOUBLE)
        THROWM(ERR.NE.0,"CANNOT WRITE REAL ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(R2) .AND. PRESENT(SNAME) ) THEN
      
        call write_array_hdf(G_ID,SNAME,LCOMMENT,err,units=units,r2=r2, COMPRESS=COMPRESS, FORCE_DOUBLE=FORCE_DOUBLE)
        THROWM(ERR.NE.0,"CANNOT WRITE REAL ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(RM2) .AND. PRESENT(SNAME) ) THEN
      
        call write_array_hdf(G_ID,SNAME,LCOMMENT,err,units=units,rm2=rm2, COMPRESS=COMPRESS, FORCE_DOUBLE=FORCE_DOUBLE)
        THROWM(ERR.NE.0,"CANNOT WRITE REAL ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(R3) .AND. PRESENT(SNAME) ) THEN
      
        call write_array_hdf(G_ID,SNAME,LCOMMENT,err,units=units,r3=r3, COMPRESS=COMPRESS, FORCE_DOUBLE=FORCE_DOUBLE)
        THROWM(ERR.NE.0,"CANNOT WRITE REAL ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(RM3) .AND. PRESENT(SNAME) ) THEN
      
        call write_array_hdf(G_ID,SNAME,LCOMMENT,err,units=units,rm3=rm3, COMPRESS=COMPRESS, FORCE_DOUBLE=FORCE_DOUBLE)
        THROWM(ERR.NE.0,"CANNOT WRITE REAL ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(R4) .AND. PRESENT(SNAME) ) THEN
      
        call write_array_hdf(G_ID,SNAME,LCOMMENT,err,units=units,r4=r4, COMPRESS=COMPRESS, FORCE_DOUBLE=FORCE_DOUBLE)
        THROWM(ERR.NE.0,"CANNOT WRITE REAL ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(RM4) .AND. PRESENT(SNAME) ) THEN
      
        call write_array_hdf(G_ID,SNAME,LCOMMENT,err,units=units,rm4=rm4, COMPRESS=COMPRESS, FORCE_DOUBLE=FORCE_DOUBLE)
        THROWM(ERR.NE.0,"CANNOT WRITE REAL ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(IT1) .AND. PRESENT(SNAME) ) THEN
      
        call write_array_hdf(G_ID,SNAME,LCOMMENT,err,units=units,it1=it1, COMPRESS=COMPRESS)
        THROWM(ERR.NE.0,"CANNOT WRITE INTEGER ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF
      
      IF( PRESENT(IT2) .AND. PRESENT(SNAME) ) THEN
      
        call write_array_hdf(G_ID,SNAME,LCOMMENT,err,units=units,it2=it2, COMPRESS=COMPRESS)
        THROWM(ERR.NE.0,"CANNOT WRITE INTEGER ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF
      
      IF( PRESENT(IT3) .AND. PRESENT(SNAME) ) THEN
      
        call write_array_hdf(G_ID,SNAME,LCOMMENT,err,units=units,it3=it3, COMPRESS=COMPRESS)
        THROWM(ERR.NE.0,"CANNOT WRITE INTEGER ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(IS1) .AND. PRESENT(SNAME) ) THEN
      
        call write_array_hdf(G_ID,SNAME,LCOMMENT,err,units=units,is1=is1, COMPRESS=COMPRESS)
        THROWM(ERR.NE.0,"CANNOT WRITE INTEGER ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF
      
      IF( PRESENT(IS2) .AND. PRESENT(SNAME) ) THEN
      
        call write_array_hdf(G_ID,SNAME,LCOMMENT,err,units=units,is2=is2, COMPRESS=COMPRESS)
        THROWM(ERR.NE.0,"CANNOT WRITE INTEGER ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF
      
      IF( PRESENT(IS3) .AND. PRESENT(SNAME) ) THEN
      
        call write_array_hdf(G_ID,SNAME,LCOMMENT,err,units=units,is3=is3, COMPRESS=COMPRESS)
        THROWM(ERR.NE.0,"CANNOT WRITE INTEGER ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(I1) .AND. PRESENT(SNAME) ) THEN
      
        call write_array_hdf(G_ID,SNAME,LCOMMENT,err,units=units,i1=i1, COMPRESS=COMPRESS)
        THROWM(ERR.NE.0,"CANNOT WRITE INTEGER ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF
      
      IF( PRESENT(I2) .AND. PRESENT(SNAME) ) THEN
      
        call write_array_hdf(G_ID,SNAME,LCOMMENT,err,units=units,i2=i2, COMPRESS=COMPRESS)
        THROWM(ERR.NE.0,"CANNOT WRITE INTEGER ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF
      
      IF( PRESENT(I3) .AND. PRESENT(SNAME) ) THEN
      
        call write_array_hdf(G_ID,SNAME,LCOMMENT,err,units=units,i3=i3, COMPRESS=COMPRESS)
        THROWM(ERR.NE.0,"CANNOT WRITE INTEGER ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF
      
      IF( PRESENT(IL1) .AND. PRESENT(SNAME) ) THEN
      
        call write_array_hdf(G_ID,SNAME,LCOMMENT,err,units=units,il1=il1, COMPRESS=COMPRESS)
        THROWM(ERR.NE.0,"CANNOT WRITE INTEGER ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(L1) .AND. PRESENT(SNAME) ) THEN
      
        call write_array_hdf(G_ID,SNAME,LCOMMENT,err,units=units,l1=l1)
        THROWM(ERR.NE.0,"CANNOT WRITE LOGICAL ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF
      
      IF( PRESENT(C0) .AND. PRESENT(SNAME) ) THEN
      
        call write_array_hdf(G_ID,SNAME,LCOMMENT,err,units=units,c0=c0)
        THROWM(ERR.NE.0,"CANNOT WRITE STRING "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF
      
      IF( PRESENT(I0) .AND. PRESENT(SNAME) ) THEN
      
        call write_array_hdf(G_ID,SNAME,LCOMMENT,err,units=units,i0=i0)
        THROWM(ERR.NE.0,"CANNOT WRITE INTEGER "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(R0) .AND. PRESENT(SNAME) ) THEN
      
        call write_array_hdf(G_ID,SNAME,LCOMMENT,err,units=units,r0=r0, FORCE_DOUBLE=FORCE_DOUBLE)
        THROWM(ERR.NE.0,"CANNOT WRITE REAL SCALAR "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(RM0) .AND. PRESENT(SNAME) ) THEN
      
        call write_array_hdf(G_ID,SNAME,LCOMMENT,err,units=units,rm0=rm0, FORCE_DOUBLE=FORCE_DOUBLE)
        THROWM(ERR.NE.0,"CANNOT WRITE REAL SCALAR "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF

      IF(.NOT. PATHROOT) THEN
        CALL H5GCLOSE_F( G_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP "//TRIM(PATH)//" IN FILE "//TRIM(F))
      ENDIF
        
      !INFO("avant close "//TRIM(F))
      !CALL RTTOV_HDF_LISTOBJ( ERR, FILE_ID)
      !THROWM(ERR.NE.0,"RTTOV_HDF_LISTOBJ on file :"//TRIM(F))


      CALL H5FCLOSE_F( FILE_ID, ERR )
      THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(F))


CATCH
      END SUBROUTINE
