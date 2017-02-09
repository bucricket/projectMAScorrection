!
SUBROUTINE rttov_read_ascii_pccoef( &
            & err,           &
            & opts,          &
            & coef,          &
            & coef_pccomp,   &
            & file_lu,       &
            & channels,      &
            & channels_rec)! in Optional
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
!  1.4       15/06/2004  Corrected array dimension for coef%fmv_gas_pos (R Saunders)
!  1.5        June 2005  Added Additional arrays for RTTOV8M Marco Matricardi
!  1.6       05/12/2005  Corrected some array types for optical  const (R Saunders)
!  1.7       19/03/2007  Reduced no of continuation lines to below 39 (R Saunders)
!  1.8       12/12/2007  Add option of using top level (R Saunders)
!  1.9       01/11/2007  Removed hardcoded section length (A Geer)
!  1.10      26/06/2008  Introduced the case where no channels are available for the
!                        phase function in the solar range (M. Matricardi)
!  1.11      27/02/2009  Profile levels to include ToA. Allocate coef arrays
!                        according to number of layers, not levels (P. Rayer)
!  1.12      06/03/2009  Separation of flags for IncZeeman and IncTop.
!                        Conditionals depending on coef%id_comp_lvl == 9
!                        extended to >= 9 (P. Rayer)
!  1.13      08/06/2009  Made interim fix to allocate cloud/aerosol arrays with right shape (R Saunders)
!                        Ideally the channel order in all IR files needs to be in frequency not wavelength
!  1.14      02/12/2009  Introduced principal component capability. Marco Matricardi. ECMWF
!            18/10/2010  Control on channels list against reference regression channels
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
! 2011/11 Introduced PC multiple band and cloudy computations, Marco Matricardi, ECMWF
!
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
! Imported Type Definitions:
  USE rttov_types, ONLY :  &
       & rttov_coef,       &
       & rttov_coef_pccomp,&
       & rttov_options
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE rttov_const, ONLY :  &
       & sensor_id_hi,           &
       & lensection
  USE parkind1, ONLY : jprb, jplm
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(in):
  INTEGER(KIND=jpim)       , INTENT(IN)                :: file_lu        ! file logical unit number
  INTEGER(KIND=jpim)       , OPTIONAL     , INTENT(IN) :: channels    (:)! list of channels to extract
  INTEGER(KIND=jpim)       , OPTIONAL     , INTENT(IN) :: channels_rec(:)
  TYPE(rttov_options      ), INTENT(IN)                :: opts
! scalar arguments with intent(inout):
  TYPE(rttov_coef         ), INTENT(INOUT)             :: coef          ! coefficients
  TYPE(rttov_coef_pccomp  ), INTENT(INOUT)             :: coef_pccomp
! scalar arguments with intent(out):
  INTEGER(KIND=jpim)       , INTENT(OUT)               :: err           ! return code
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_skipcommentline.interface"
#include "rttov_findnextsection.interface"
#include "rttov_nullify_coef_pccomp.interface"
! Local Scalars:
  INTEGER(KIND=jpim) :: file_channels
  INTEGER(KIND=jpim) :: file_channels_rec
  LOGICAL(KIND=jplm) :: all_channels
  LOGICAL(KIND=jplm) :: all_channels_rec
  INTEGER(KIND=jpim) :: io_status
  INTEGER(KIND=jpim) :: i, j, m, n
! pointers for generic inputs
  REAL   (KIND=jprb), POINTER :: values0      (:)
  REAL   (KIND=jprb), POINTER :: values1      (:)
  REAL   (KIND=jprb), POINTER :: values2      (:)
  REAL   (KIND=jprb), POINTER :: values3      (:)
  REAL   (KIND=jprb), POINTER :: values4      (:)
  REAL   (KIND=jprb), POINTER :: values5      (:)
  REAL   (KIND=jprb), POINTER :: values6      (:)
  REAL   (KIND=jprb), POINTER :: values7      (:)
  REAL   (KIND=jprb), POINTER :: values8      (:)
  INTEGER(KIND=jpim), POINTER :: ivalues0     (:)
  REAL   (KIND=jprb), POINTER :: eigenarray   (:,:)
  REAL   (KIND=jprb), POINTER :: noisearray   (:)
  REAL   (KIND=jprb), POINTER :: pcbcoarray   (:)
  REAL   (KIND=jprb), POINTER :: pcbcsarray   (:)
  REAL   (KIND=jprb), POINTER :: pccwnarray   (:)
  INTEGER(KIND=jpim), POINTER :: pcchnarray   (:)
  CHARACTER(LEN = lensection) :: section
!- End of header --------------------------------------------------------
  TRY
! 0 Initialise variables
!---------------------------------------------
! test presence of channels argument

  IF (coef%id_sensor /= sensor_id_hi) THEN
    err = errorstatus_fatal

    THROWM( ERR .NE. 0, "PC-RTTOV only compatible with hi-res sounders")
  ENDIF

  IF (Present(channels)) THEN
    all_channels = .FALSE.
  ELSE
    all_channels = .TRUE.
  ENDIF


  IF (Present(channels_rec)) THEN
    all_channels_rec = .FALSE.
  ELSE
    all_channels_rec = .TRUE.
  ENDIF

  CALL rttov_nullify_coef_pccomp(coef_pccomp)

!read the file

  readfile : DO
    CALL rttov_findnextsection(file_lu, io_status, section)
    IF (io_status < 0) EXIT!end-of-file

    SELECT CASE (Trim(section))
    CASE ('PRINCOMP_PREDICTORS')
      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM( ERR .NE. 0, 'io status while reading section '//section)


      READ (file_lu,  * , iostat=ERR)coef_pccomp%fmv_pc_comp_pc

      THROWM( ERR .NE. 0, 'io status while reading section '//section)


      IF (coef%id_comp_pc /= coef_pccomp%fmv_pc_comp_pc) THEN

        err = errorstatus_fatal

        THROWM( ERR .NE. 0, "Version of PC coef file is incompatible with RTTOV regression file")

      ENDIF

      READ (file_lu,  * , iostat=ERR)coef_pccomp%fmv_pc_cld

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

      IF (opts%rt_ir%pc%addpc .AND. opts%rt_ir%addclouds) THEN

        IF(coef_pccomp%fmv_pc_cld==0)THEN

          err = errorstatus_fatal

          THROWM( ERR .NE. 0, "Version of PC coef file is incompatible with cloudy computations")

        ENDIF

      ENDIF

      READ (file_lu,  * , iostat=ERR)coef_pccomp%fmv_pc_bands

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

      ALLOCATE (coef_pccomp%fmv_pc_sets(coef_pccomp%fmv_pc_bands), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of coef_pccomp%fmv_pc_sets array")

      ALLOCATE (coef_pccomp%eigen(coef_pccomp%fmv_pc_bands), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of coef_pccomp%eigen")
      

      READ (file_lu,  * , iostat=ERR)coef_pccomp%fmv_pc_msets

      THROWM( ERR .NE. 0, 'io status while reading section '//section)


      ALLOCATE (coef_pccomp%pcreg(coef_pccomp%fmv_pc_bands,coef_pccomp%fmv_pc_msets), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of coef_pccomp%pcreg arrays")


      READ (file_lu,  * , iostat=ERR)coef_pccomp%fmv_pc_mnum

      THROWM( ERR .NE. 0, 'io status while reading section '//section)


      READ (file_lu,  * , iostat=ERR)coef_pccomp%fmv_pc_mchn

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

! loop on bands 

      DO m = 1, coef_pccomp%fmv_pc_bands

        READ (file_lu,  * , iostat=ERR)coef_pccomp%fmv_pc_sets(m)

        THROWM( ERR .NE. 0, 'io status while reading section '//section)

! loop on predictor sets

        DO n = 1,  coef_pccomp%fmv_pc_sets(m)

          READ (file_lu,  * , iostat=ERR)coef_pccomp%pcreg(m,n)%fmv_pc_npred

          THROWM( ERR .NE. 0, 'io status while reading section '//section)


          CALL rttov_skipcommentline(file_lu, ERR)

          THROWM( ERR .NE. 0, 'io status while reading section '//section)

          ALLOCATE (coef_pccomp%pcreg(m,n)%predictindex(coef_pccomp%pcreg(m,n)%fmv_pc_npred), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of predictindex arrays")

          READ (file_lu,  * , iostat=ERR) &
            & (coef_pccomp%pcreg(m,n)%predictindex(i), i = 1, coef_pccomp%pcreg(m,n)%fmv_pc_npred)

          THROWM( ERR .NE. 0, "io status while reading predictindex arrays")

        ENDDO
      ENDDO



    CASE ('PRINCOMP_EIGENVECTORS')

      coef_pccomp%fmv_pc_nchn = coef_pccomp%fmv_pc_mchn

! Take care of the user list of channels
! file_channels store the number of channels in the file
! coef%fmv_chn is the number of channels that the user requests
      file_channels_rec = coef_pccomp%fmv_pc_nchn

      IF (.NOT. all_channels_rec) THEN
        coef_pccomp%fmv_pc_nchn = Size(channels_rec)
      ENDIF


      IF (.NOT. all_channels) THEN
        coef_pccomp%fmv_pc_nchn_noise = Size(channels)
      ELSE
        coef_pccomp%fmv_pc_nchn_noise = coef_pccomp%fmv_pc_mchn
      ENDIF

      DO m = 1,  coef_pccomp%fmv_pc_bands

        CALL rttov_skipcommentline(file_lu, ERR)

        THROWM( ERR .NE. 0, 'io status while reading section '//section)

        NULLIFY(coef_pccomp%eigen(m)%eigenvectors_t)

        ALLOCATE (coef_pccomp%eigen(m)%eigenvectors(coef_pccomp%fmv_pc_nchn, coef_pccomp%fmv_pc_mnum), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of eigenvectors")

        IF (all_channels_rec) THEN

          DO n = 1, coef_pccomp%fmv_pc_mnum
            READ (file_lu,  * , iostat=ERR)(coef_pccomp%eigen(m)%eigenvectors(i, n), i = 1, coef_pccomp%fmv_pc_nchn)

            THROWM( ERR .NE. 0, 'io status while reading section '//section)

          ENDDO

        ELSE

          ALLOCATE (eigenarray(file_channels_rec, coef_pccomp%fmv_pc_mnum), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of eigenarray")


          DO n = 1, coef_pccomp%fmv_pc_mnum
            READ (file_lu,  * , iostat=ERR)(eigenarray(i, n), i = 1, file_channels_rec)

            THROWM( ERR .NE. 0, 'io status while reading section '//section)

          ENDDO


          coef_pccomp%eigen(m)%eigenvectors(:,:) = eigenarray(channels_rec(:), :)
        ENDIF
      ENDDO

    CASE ('PRINCOMP_COEFFICIENTS')



      DO m=1,coef_pccomp%fmv_pc_bands

        DO n = 1, coef_pccomp%fmv_pc_sets(m)

          NULLIFY(coef_pccomp%pcreg(m,n)%coefficients_t)

          ALLOCATE (coef_pccomp%pcreg(m,n)%coefficients(coef_pccomp%pcreg(m,n)%fmv_pc_npred, coef_pccomp%fmv_pc_mnum), STAT =      &
            & ERR)

          THROWM( ERR .NE. 0, "allocation of pcreg(n)%coefficients")

          CALL rttov_skipcommentline(file_lu, ERR)

          THROWM( ERR .NE. 0, 'io status while reading section '//section)


          DO j = 1, coef_pccomp%fmv_pc_mnum

            READ (file_lu,  * , iostat=ERR)     &
              & (coef_pccomp%pcreg(m,n)%coefficients(i, j), i = 1, coef_pccomp%pcreg(m,n)%fmv_pc_npred)
            THROWM( ERR .NE. 0, 'io status while reading section '//section)

          ENDDO

        ENDDO
      ENDDO


    CASE ('EMISSIVITY_COEFFICIENTS')
      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

      READ (file_lu,  * , iostat=ERR)coef_pccomp%fmv_pc_nche

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

! Take care of the user list of channels
! file_channels store the number of channels in the file
! coef%fmv_chn is the number of channels that the user requests
      file_channels = coef_pccomp%fmv_pc_nche

      IF (.NOT. all_channels) THEN
        coef_pccomp%fmv_pc_nche = Size(channels)
      ENDIF

      ALLOCATE (coef_pccomp%emiss_chn(coef_pccomp%fmv_pc_nche), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of emiss_chn")

      ALLOCATE (coef_pccomp%emiss_c1(coef_pccomp%fmv_pc_nche), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of emiss_c1")

      ALLOCATE (coef_pccomp%emiss_c2(coef_pccomp%fmv_pc_nche), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of emiss_c2")

      ALLOCATE (coef_pccomp%emiss_c3(coef_pccomp%fmv_pc_nche), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of emiss_c3")

      ALLOCATE (coef_pccomp%emiss_c4(coef_pccomp%fmv_pc_nche), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of emiss_c4")

      ALLOCATE (coef_pccomp%emiss_c5(coef_pccomp%fmv_pc_nche), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of emiss_c5")

      ALLOCATE (coef_pccomp%emiss_c6(coef_pccomp%fmv_pc_nche), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of emiss_c6")

      ALLOCATE (coef_pccomp%emiss_c7(coef_pccomp%fmv_pc_nche), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of emiss_c7")

      ALLOCATE (coef_pccomp%emiss_c8(coef_pccomp%fmv_pc_nche), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of emiss_c8")

      ALLOCATE (coef_pccomp%emiss_c9(coef_pccomp%fmv_pc_nche), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of emiss_c9")


      IF (all_channels) THEN

        DO i = 1, coef_pccomp%fmv_pc_nche
          READ (file_lu,  * , iostat=ERR)coef_pccomp%emiss_chn(i), coef_pccomp%emiss_c1(i), coef_pccomp%emiss_c2(i), &
            & coef_pccomp%emiss_c3(i), coef_pccomp%emiss_c4(i), coef_pccomp%emiss_c5(i), coef_pccomp%emiss_c6(i),    &
            & coef_pccomp%emiss_c7(i), coef_pccomp%emiss_c8(i), coef_pccomp%emiss_c9(i)

          THROWM( ERR .NE. 0, 'io status while reading section '//section)

        ENDDO

      ELSE
        ALLOCATE (ivalues0(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%emiss iv0")

        ALLOCATE (values0(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%emiss v0")

        ALLOCATE (values1(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%emiss v1")

        ALLOCATE (values2(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%emiss v2")

        ALLOCATE (values3(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%emiss v3")

        ALLOCATE (values4(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%emiss v4")

        ALLOCATE (values5(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%emiss v5")

        ALLOCATE (values6(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%emiss v6")

        ALLOCATE (values7(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%emiss v7")

        ALLOCATE (values8(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%emiss v8")


        DO i = 1, file_channels
          READ (file_lu,  * , iostat=ERR)ivalues0(i), values0(i), values1(i), values2(i), values3(i), values4(i), &
            & values5(i), values6(i), values7(i), values8(i)

          THROWM(err.ne.0,"io status while reading section "//section)

        ENDDO

        coef_pccomp%emiss_chn(:) = ivalues0(channels(:))
        coef_pccomp%emiss_c1(:)  = values0(channels(:))
        coef_pccomp%emiss_c2(:)  = values1(channels(:))
        coef_pccomp%emiss_c3(:)  = values2(channels(:))
        coef_pccomp%emiss_c4(:)  = values3(channels(:))
        coef_pccomp%emiss_c5(:)  = values4(channels(:))
        coef_pccomp%emiss_c6(:)  = values5(channels(:))
        coef_pccomp%emiss_c7(:)  = values6(channels(:))
        coef_pccomp%emiss_c8(:)  = values7(channels(:))
        coef_pccomp%emiss_c9(:)  = values8(channels(:))
        DEALLOCATE (ivalues0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of iv0")

        DEALLOCATE (values0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of v0")

        DEALLOCATE (values1, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of v1")

        DEALLOCATE (values2, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of v2")

        DEALLOCATE (values3, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of v3")

        DEALLOCATE (values4, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of v4")

        DEALLOCATE (values5, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of v5")

        DEALLOCATE (values6, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of v6")

        DEALLOCATE (values7, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of v7")

        DEALLOCATE (values8, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of v8")

      ENDIF

    CASE ('PC_REFERENCE_PROFILE')
      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

      READ (file_lu,  * , iostat=ERR)coef_pccomp%fmv_pc_gas

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

      READ (file_lu,  * , iostat=ERR)coef_pccomp%fmv_pc_nlev

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

      ALLOCATE (coef_pccomp%ref_pc_prfl_p(coef_pccomp%fmv_pc_nlev), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of pccomp%ref_pc_prfl_p")

      ALLOCATE (coef_pccomp%ref_pc_prfl_mr(coef_pccomp%fmv_pc_nlev, coef_pccomp%fmv_pc_gas), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of pccomp%ref_pc_prfl_mr")


      DO n = 1, coef_pccomp%fmv_pc_gas
        CALL rttov_skipcommentline(file_lu, ERR)

        THROWM( ERR .NE. 0, 'io status while reading section '//section)


        DO i = 1, coef_pccomp%fmv_pc_nlev
          READ (file_lu,  * , iostat=ERR)coef_pccomp%ref_pc_prfl_p(i), coef_pccomp%ref_pc_prfl_mr(i, n)

          THROWM( ERR .NE. 0, 'io status while reading section '//section)

        ENDDO

      ENDDO

    CASE ('PC_PROFILE_LIMITS')
      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

      ALLOCATE (coef_pccomp%lim_pc_prfl_tmax(coef_pccomp%fmv_pc_nlev), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of pccomp%lim_pc_prfl_tmax")

      ALLOCATE (coef_pccomp%lim_pc_prfl_tmin(coef_pccomp%fmv_pc_nlev), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of pccomp%lim_pc_prfl_tmin")

      ALLOCATE (coef_pccomp%lim_pc_prfl_qmax(coef_pccomp%fmv_pc_nlev), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of pccomp%lim_pc_prfl_qmax")

      ALLOCATE (coef_pccomp%lim_pc_prfl_qmin(coef_pccomp%fmv_pc_nlev), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of pccomp%lim_pc_prfl_qmin")

      ALLOCATE (coef_pccomp%lim_pc_prfl_ozmax(coef_pccomp%fmv_pc_nlev), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of pccomp%lim_pc_prfl_ozmax")

      ALLOCATE (coef_pccomp%lim_pc_prfl_ozmin(coef_pccomp%fmv_pc_nlev), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of pccomp%lim_pc_prfl_ozmin")


      DO i = 1, coef_pccomp%fmv_pc_nlev
        READ (file_lu,  * , iostat=ERR)     &
          & coef_pccomp%ref_pc_prfl_p(i), coef_pccomp%lim_pc_prfl_tmin(i), coef_pccomp%lim_pc_prfl_tmax(i)

        THROWM( ERR .NE. 0, 'io status while reading section '//section)

      ENDDO

      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM( ERR .NE. 0, 'io status while reading section '//section)


      DO i = 1, coef_pccomp%fmv_pc_nlev
        READ (file_lu,  * , iostat=ERR)     &
          & coef_pccomp%ref_pc_prfl_p(i), coef_pccomp%lim_pc_prfl_qmin(i), coef_pccomp%lim_pc_prfl_qmax(i)

        THROWM( ERR .NE. 0, 'io status while reading section '//section)

      ENDDO

      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM( ERR .NE. 0, 'io status while reading section '//section)


      DO i = 1, coef_pccomp%fmv_pc_nlev
        READ (file_lu,  * , iostat=ERR)     &
          & coef_pccomp%ref_pc_prfl_p(i), coef_pccomp%lim_pc_prfl_ozmin(i), coef_pccomp%lim_pc_prfl_ozmax(i)

        THROWM( ERR .NE. 0, 'io status while reading section '//section)

      ENDDO

      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

      READ (file_lu,  * , iostat=ERR)coef_pccomp%lim_pc_prfl_pmin, coef_pccomp%lim_pc_prfl_pmax

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

      READ (file_lu,  * , iostat=ERR)coef_pccomp%lim_pc_prfl_tsmin, coef_pccomp%lim_pc_prfl_tsmax

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

      READ (file_lu,  * , iostat=ERR)coef_pccomp%lim_pc_prfl_skmin, coef_pccomp%lim_pc_prfl_skmax

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

      READ (file_lu,  * , iostat=ERR)coef_pccomp%lim_pc_prfl_wsmin, coef_pccomp%lim_pc_prfl_wsmax

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

    CASE ('INSTRUMENT_NOISE')
      ALLOCATE (coef_pccomp%noise_in(coef_pccomp%fmv_pc_nchn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of pccomp%noise_in")

      ALLOCATE (coef_pccomp%ff_ori_chn_in(coef_pccomp%fmv_pc_nchn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of pccomp%ff_ori_chn_in")

      ALLOCATE (coef_pccomp%ff_cwn_in(coef_pccomp%fmv_pc_nchn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of pccomp%ff_cwn_in")

      ALLOCATE (coef_pccomp%ff_bco_in(coef_pccomp%fmv_pc_nchn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of pccomp%ff_bco_in")

      ALLOCATE (coef_pccomp%ff_bcs_in(coef_pccomp%fmv_pc_nchn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of pccomp%ff_bcs_in")

      ALLOCATE (coef_pccomp%noise(coef_pccomp%fmv_pc_nchn_noise), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of pccomp%noise")

      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM( ERR .NE. 0, 'io status while reading section '//section)


      IF (all_channels_rec) THEN

        DO n = 1, coef_pccomp%fmv_pc_nchn

          READ (file_lu,  * , iostat=ERR)coef_pccomp%ff_ori_chn_in(n), coef_pccomp%ff_cwn_in(n),      &
            & coef_pccomp%ff_bco_in(n), coef_pccomp%ff_bcs_in(n), coef_pccomp%noise_in(n)

          THROWM( ERR .NE. 0, 'io status while reading section '//section)
  
        ENDDO

        IF (all_channels) THEN
          coef_pccomp%noise = coef_pccomp%noise_in
        ELSE
          coef_pccomp%noise = coef_pccomp%noise_in(coef%ff_ori_chn(:))
        ENDIF

      ELSE
        ALLOCATE (noisearray(file_channels_rec), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of noisearray")

        ALLOCATE (pcchnarray(file_channels_rec), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of pcchnarray")

        ALLOCATE (pccwnarray(file_channels_rec), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of pccwnarray")

        ALLOCATE (pcbcoarray(file_channels_rec), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of pcbcoarray")

        ALLOCATE (pcbcsarray(file_channels_rec), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of pcbcsarray")


        DO n = 1, file_channels_rec
          READ (file_lu,  * , iostat=ERR)pcchnarray(n), pccwnarray(n), pcbcoarray(n), pcbcsarray(n), noisearray(n)

          THROWM( ERR .NE. 0, 'io status while reading section '//section)

        ENDDO

        coef_pccomp%ff_cwn_in(:)     = pccwnarray(channels_rec(:))
        coef_pccomp%ff_ori_chn_in(:) = pcchnarray(channels_rec(:))
        coef_pccomp%ff_bco_in(:)     = pcbcoarray(channels_rec(:))
        coef_pccomp%ff_bcs_in(:)     = pcbcsarray(channels_rec(:))
        coef_pccomp%noise_in(:)      = noisearray(channels_rec(:))

        coef_pccomp%noise(:)         = noisearray(coef%ff_ori_chn(:))

        DEALLOCATE (noisearray, pccwnarray, pcchnarray, pcbcoarray, pcbcsarray, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of pcbcsarray")
      ENDIF

    CASE ('END')
      RETURN
    CASE DEFAULT
      CYCLE readfile
    END SELECT

  ENDDO readfile

  CATCH
END SUBROUTINE rttov_read_ascii_pccoef
