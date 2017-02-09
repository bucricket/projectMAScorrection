!
SUBROUTINE rttov_read_binary_pccoef( &
            & ERR,           &
            & opts,          &
            & coef,          &
            & coef_pccomp,   &
            & file_lu,       &
            & channels,      &
            & channels_rec)
! Description:
!
! Read an binary coefficient file and fills coeff structure
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
!  1.2       24/01/2003  add tests on all read statements (P Brunel)
!                        one record per channel for coefficients in binary format
!                        New header to allow checking R4<->R8
!  1.3       06/05/2003  Remove "optional" attribute of argument file_lu (P Brunel)
!  1.4       02/06/2004  New format for FMV section with RTTOV8  (P. Brunel)
!  1.5       15/06/2004  Corrected array dimension for coef%fmv_gas_pos (R Saunders)
!  1.6       14/05/2007  Updated for RTTOV-9 ( P Brunel)
!  1.7       12/12/2007  Option of additional layer at top flag read (R Saunders)
!  1.8       27/06/2008  Introduced the case when no channels are present for the
!                        phase function in the solar range (M. Matricardi)
!  1.9       27/02/2009  Profile levels to include ToA. Allocate coef arrays
!                        according to number of layers, not levels (P. Rayer)
!  1.10      06/03/2009  Separation of flags for IncZeeman and IncTop.
!                        Conditionals depending on coef%id_comp_lvl == 9
!                        extended to >= 9 (P. Rayer)
!  1.11      08/06/2009  Made interim fix to allocate cloud/aerosol arrays with right shape (R Saunders)
!                        Ideally the channel order in all IR files needs to be in frequency not wavelength
!  1.12      02/12/2009  Introduced principal component capability. Marco Matricardi. ECMWF
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
       & rttov_coef,          &
       & rttov_coef_pccomp,   &
       & rttov_options
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE rttov_const, ONLY :  &
       & rttov_magic_string,     &
       & rttov_magic_number,     &
       & sensor_id_hi
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
  INTEGER(KIND=jpim)       , INTENT(OUT)               :: ERR           ! return code
!INTF_END
#include "rttov_errorreport.interface"
! Local Scalars:
  INTEGER(KIND=jpim) :: file_channels
  INTEGER(KIND=jpim) :: file_channels_rec
  LOGICAL(KIND=jplm) :: all_channels
  LOGICAL(KIND=jplm) :: all_channels_rec
  INTEGER(KIND=jpim) :: i, j, m, n
! pointers for generic inputs
  REAL   (KIND=jprb), POINTER     :: values0         (:)
  REAL   (KIND=jprb), POINTER     :: values1         (:)
  REAL   (KIND=jprb), POINTER     :: values2         (:)
  REAL   (KIND=jprb), POINTER     :: values3         (:)
  REAL   (KIND=jprb), POINTER     :: values4         (:)
  REAL   (KIND=jprb), POINTER     :: values5         (:)
  REAL   (KIND=jprb), POINTER     :: values6         (:)
  REAL   (KIND=jprb), POINTER     :: values7         (:)
  REAL   (KIND=jprb), POINTER     :: values8         (:)
  INTEGER(KIND=jpim), POINTER     :: ivalues0        (:)
  REAL   (KIND=jprb), POINTER     :: eigenarray      (:,:)
  REAL   (KIND=jprb), POINTER     :: noisearray      (:)
  REAL   (KIND=jprb), POINTER     :: pcbcoarray      (:)
  REAL   (KIND=jprb), POINTER     :: pcbcsarray      (:)
  REAL   (KIND=jprb), POINTER     :: pccwnarray      (:)
  INTEGER(KIND=jpim), POINTER     :: pcchnarray      (:)
  CHARACTER(LEN = 16) :: bin_check_string
  REAL(KIND=jprb)     :: bin_check_number
  REAL(KIND=jprb)     :: bin_check_value
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

! 3 Read binary file
!-------------------
! Binary file
  READ (file_lu, iostat=ERR)bin_check_string, bin_check_number

  THROWM( ERR .NE. 0, 'io status while reading header')

! Verification of header string
  IF (bin_check_string /= rttov_magic_string) err = errorstatus_fatal

  THROWM( ERR .NE. 0,'Wrong header string in file')

! Verification of single/double precision using a 5 digit number
! with exponent 12, which is always Ok for single precision
  bin_check_value = 1._JPRB - abs(bin_check_number - rttov_magic_number)
  IF (bin_check_value > 1.01_JPRB .OR. bin_check_value < 0.99_JPRB) err = errorstatus_fatal

  THROWM( ERR .NE. 0,'File created with a different real precision (R4<->R8)')


! PRINCOMP_PREDICTORS

        READ (file_lu, iostat=ERR)coef_pccomp%fmv_pc_comp_pc

        THROWM( ERR .NE. 0, "reading coef_pccomp%fmv_pc_comp_pc")


        IF (coef%id_comp_pc /= coef_pccomp%fmv_pc_comp_pc) THEN
          err = errorstatus_fatal

          THROWM( ERR .NE. 0, "Version of PC coef file is incompatible with RTTOV regression file")

        ENDIF

        READ (file_lu, iostat=ERR)coef_pccomp%fmv_pc_cld

        THROWM( ERR .NE. 0, "reading coef_pccomp%fmv_pc_cld")

        IF (opts%rt_ir%pc%addpc .AND. opts%rt_ir%addclouds) THEN

          IF(coef_pccomp%fmv_pc_cld==0)THEN

            err = errorstatus_fatal

            THROWM( ERR .NE. 0, "Version of PC coef file is incompatible with cloudy computations")

          ENDIF

        ENDIF

        READ (file_lu, iostat=ERR)coef_pccomp%fmv_pc_bands

        THROWM( ERR .NE. 0, "reading coef_pccomp%fmv_pc_bands")

        ALLOCATE (coef_pccomp%fmv_pc_sets(coef_pccomp%fmv_pc_bands), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%fmv_pc_sets")

        ALLOCATE (coef_pccomp%eigen(coef_pccomp%fmv_pc_bands), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%eigen")

        READ (file_lu, iostat=ERR)coef_pccomp%fmv_pc_msets

        THROWM( ERR .NE. 0, "reading coef_pccomp%fmv_pc_msets")


        ALLOCATE (coef_pccomp%pcreg(coef_pccomp%fmv_pc_bands,coef_pccomp%fmv_pc_msets), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%pcreg")

        READ (file_lu, iostat=ERR)coef_pccomp%fmv_pc_mnum

        THROWM( ERR .NE. 0, "reading coef_pccomp%fmv_pc_mnum")

        READ (file_lu, iostat=ERR)coef_pccomp%fmv_pc_mchn

        THROWM( ERR .NE. 0, "reading coef_pccomp%fmv_pc_mchn")

! loop on bands
        DO m = 1, coef_pccomp%fmv_pc_bands

          READ (file_lu,iostat=ERR)coef_pccomp%fmv_pc_sets(m)

          THROWM( ERR .NE. 0, "reading coef_pccomp%fmv_pc_sets(m)")

!loop on predictor sets

          DO n = 1, coef_pccomp%fmv_pc_sets(m)

            READ (file_lu, iostat=ERR)coef_pccomp%pcreg(m,n)%fmv_pc_npred

            THROWM(err.ne.0,"reading coef_pccomp%pcreg(m,n)%fmv_pc_npred")

            ALLOCATE (coef_pccomp%pcreg(m,n)%predictindex(coef_pccomp%pcreg(m,n)%fmv_pc_npred), STAT = ERR)

            THROWM( ERR .NE. 0, "allocation of coef_pccomp%pcreg(m,n)%predictindex")

            READ (file_lu, iostat=ERR) &
              & (coef_pccomp%pcreg(m,n)%predictindex(i), i = 1, coef_pccomp%pcreg(m,n)%fmv_pc_npred)

            THROWM(err.ne.0,"reading coef_pccomp%pcreg(m,n)%predictindex")

          ENDDO
        ENDDO

! PRINCOMP_EIGENVECTORS
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
        DO m = 1, coef_pccomp%fmv_pc_bands

          NULLIFY(coef_pccomp%eigen(m)%eigenvectors_t)

          ALLOCATE (coef_pccomp%eigen(m)%eigenvectors(coef_pccomp%fmv_pc_nchn, coef_pccomp%fmv_pc_mnum), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of coef_pccomp%eigen(m)%eigenvectors")

          IF (all_channels_rec) THEN

            DO n = 1, coef_pccomp%fmv_pc_mnum
              READ (file_lu, iostat=ERR)(coef_pccomp%eigen(m)%eigenvectors(i, n), i = 1, coef_pccomp%fmv_pc_nchn)

              THROWM(err.ne.0,"reading coef_pccomp%eigen(m)%eigenvectors")

            ENDDO

          ELSE
            ALLOCATE (eigenarray(file_channels_rec, coef_pccomp%fmv_pc_mnum), STAT = ERR)

            THROWM( ERR .NE. 0, "allocation of eigenarray")


            DO n = 1, coef_pccomp%fmv_pc_mnum
              READ (file_lu, iostat=ERR)(eigenarray(i, n), i = 1, file_channels_rec)

              THROWM(err.ne.0,"reading eigenarray")

            ENDDO

            coef_pccomp%eigen(m)%eigenvectors(:,:) = eigenarray(channels_rec(:), :)
            DEALLOCATE (eigenarray, STAT = ERR)

            THROWM( ERR .NE. 0, "deallocation of eigenarray")

          ENDIF
        ENDDO

! PRINCOMP_COEFFICIENTS

        DO m = 1, coef_pccomp%fmv_pc_bands

          DO n = 1, coef_pccomp%fmv_pc_sets(m)
            NULLIFY(coef_pccomp%pcreg(m,n)%coefficients_t)

            ALLOCATE (coef_pccomp%pcreg(m,n)%coefficients(coef_pccomp%pcreg(m,n)%fmv_pc_npred, &
                                                          coef_pccomp%fmv_pc_mnum), STAT = ERR)

            THROWM( ERR .NE. 0, "allocation of coef_pccomp%pcreg(m,n)%coefficients")


            DO j = 1, coef_pccomp%fmv_pc_mnum
              READ (file_lu, iostat=ERR)     &
                & (coef_pccomp%pcreg(m,n)%coefficients(i, j), i = 1, coef_pccomp%pcreg(m,n)%fmv_pc_npred)

              THROWM(err.ne.0,"reading coef_pccomp%pcreg(m,n)%coefficients")

            ENDDO

          ENDDO
        ENDDO

!     EMISSIVITY_COEFFICIENTS
        READ (file_lu, iostat=ERR)coef_pccomp%fmv_pc_nche

        THROWM(err.ne.0,"reading coef_pccomp%fmv_pc_nche")

! Take care of the user list of channels
! file_channels store the number of channels in the file
! coef%fmv_chn is the number of channels that the user requests
        file_channels = coef_pccomp%fmv_pc_nche

        IF (.NOT. all_channels) THEN
          coef_pccomp%fmv_pc_nche = Size(channels)
        ENDIF

        ALLOCATE (coef_pccomp%emiss_chn(coef_pccomp%fmv_pc_nche), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%emiss_chn")

        ALLOCATE (coef_pccomp%emiss_c1(coef_pccomp%fmv_pc_nche), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%emiss_c1")

        ALLOCATE (coef_pccomp%emiss_c2(coef_pccomp%fmv_pc_nche), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%emiss_c2")

        ALLOCATE (coef_pccomp%emiss_c3(coef_pccomp%fmv_pc_nche), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%emiss_c3")

        ALLOCATE (coef_pccomp%emiss_c4(coef_pccomp%fmv_pc_nche), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%emiss_c4")

        ALLOCATE (coef_pccomp%emiss_c5(coef_pccomp%fmv_pc_nche), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%emiss_c5")

        ALLOCATE (coef_pccomp%emiss_c6(coef_pccomp%fmv_pc_nche), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%emiss_c6")

        ALLOCATE (coef_pccomp%emiss_c7(coef_pccomp%fmv_pc_nche), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%emiss_c7")

        ALLOCATE (coef_pccomp%emiss_c8(coef_pccomp%fmv_pc_nche), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%emiss_c8")

        ALLOCATE (coef_pccomp%emiss_c9(coef_pccomp%fmv_pc_nche), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%emiss_c9")


        IF (all_channels) THEN

          DO i = 1, coef_pccomp%fmv_pc_nche
            READ (file_lu, iostat=ERR)coef_pccomp%emiss_chn(i), coef_pccomp%emiss_c1(i), coef_pccomp%emiss_c2(i),   &
              & coef_pccomp%emiss_c3(i), coef_pccomp%emiss_c4(i), coef_pccomp%emiss_c5(i), coef_pccomp%emiss_c6(i), &
              & coef_pccomp%emiss_c7(i), coef_pccomp%emiss_c8(i), coef_pccomp%emiss_c9(i)

            THROWM(err.ne.0,"reading coef_pccomp%emiss_chn")

          ENDDO

        ELSE
          ALLOCATE (ivalues0(file_channels), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of ivalues0")

          ALLOCATE (values0(file_channels), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of values0")

          ALLOCATE (values1(file_channels), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of values1")

          ALLOCATE (values2(file_channels), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of values2")

          ALLOCATE (values3(file_channels), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of values3")

          ALLOCATE (values4(file_channels), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of values4")

          ALLOCATE (values5(file_channels), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of values5")

          ALLOCATE (values6(file_channels), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of values6")

          ALLOCATE (values7(file_channels), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of values7")

          ALLOCATE (values8(file_channels), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of values8")


          DO i = 1, file_channels
            READ (file_lu, iostat=ERR)ivalues0(i), values0(i), values1(i), values2(i), values3(i), values4(i), &
              & values5(i), values6(i), values7(i), values8(i)

            THROWM(err.ne.0,"reading ivalues0...")

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

          THROWM( ERR .NE. 0, "deallocation of ivalues0")

          DEALLOCATE (values0, STAT = ERR)

          THROWM( ERR .NE. 0, "deallocation of values0")

          DEALLOCATE (values1, STAT = ERR)

          THROWM( ERR .NE. 0, "deallocation of values1,")

          DEALLOCATE (values2, STAT = ERR)

          THROWM( ERR .NE. 0, "deallocation of values2")

          DEALLOCATE (values3, STAT = ERR)

          THROWM( ERR .NE. 0, "deallocation of values3")

          DEALLOCATE (values4, STAT = ERR)

          THROWM( ERR .NE. 0, "deallocation of values4")

          DEALLOCATE (values5, STAT = ERR)

          THROWM( ERR .NE. 0, "deallocation of values5")

          DEALLOCATE (values6, STAT = ERR)

          THROWM( ERR .NE. 0, "deallocation of values6")

          DEALLOCATE (values7, STAT = ERR)

          THROWM( ERR .NE. 0, "deallocation of values7")

          DEALLOCATE (values8, STAT = ERR)

          THROWM( ERR .NE. 0, "deallocation of values8")

        ENDIF

!     PC_REFERENCE_PROFILE
        READ (file_lu, iostat=ERR)coef_pccomp%fmv_pc_gas

        THROWM(err.ne.0,"reading coef_pccomp%fmv_pc_gas")

        READ (file_lu, iostat=ERR)coef_pccomp%fmv_pc_nlev

        THROWM(err.ne.0,"reading coef_pccomp%fmv_pc_nlev")

        ALLOCATE (coef_pccomp%ref_pc_prfl_p(coef_pccomp%fmv_pc_nlev), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%ref_pc_prfl_p")

        ALLOCATE (coef_pccomp%ref_pc_prfl_mr(coef_pccomp%fmv_pc_nlev, coef_pccomp%fmv_pc_gas), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%ref_pc_prfl_mr")


        DO n = 1, coef_pccomp%fmv_pc_gas

          DO i = 1, coef_pccomp%fmv_pc_nlev
            READ (file_lu, iostat=ERR)coef_pccomp%ref_pc_prfl_p(i), coef_pccomp%ref_pc_prfl_mr(i, n)

            THROWM(err.ne.0,"reading coef_pccomp%ref_pc_prfl_mr")

          ENDDO

        ENDDO

!     PC_PROFILE_LIMITS
        ALLOCATE (coef_pccomp%lim_pc_prfl_tmax(coef_pccomp%fmv_pc_nlev), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%lim_pc_prfl_tmax")

        ALLOCATE (coef_pccomp%lim_pc_prfl_tmin(coef_pccomp%fmv_pc_nlev), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%lim_pc_prfl_tmin")

        ALLOCATE (coef_pccomp%lim_pc_prfl_qmax(coef_pccomp%fmv_pc_nlev), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%lim_pc_prfl_qmax")

        ALLOCATE (coef_pccomp%lim_pc_prfl_qmin(coef_pccomp%fmv_pc_nlev), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%lim_pc_prfl_qmin")

        ALLOCATE (coef_pccomp%lim_pc_prfl_ozmax(coef_pccomp%fmv_pc_nlev), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%lim_pc_prfl_ozmax")

        ALLOCATE (coef_pccomp%lim_pc_prfl_ozmin(coef_pccomp%fmv_pc_nlev), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%lim_pc_prfl_ozmin")


        DO i = 1, coef_pccomp%fmv_pc_nlev
          READ (file_lu, iostat=ERR)     &
            & coef_pccomp%ref_pc_prfl_p(i), coef_pccomp%lim_pc_prfl_tmin(i), coef_pccomp%lim_pc_prfl_tmax(i)

          THROWM(err.ne.0,"reading coef_pccomp%ref_pc_prfl_tmin")

        ENDDO


        DO i = 1, coef_pccomp%fmv_pc_nlev
          READ (file_lu, iostat=ERR)     &
            & coef_pccomp%ref_pc_prfl_p(i), coef_pccomp%lim_pc_prfl_qmin(i), coef_pccomp%lim_pc_prfl_qmax(i)

          THROWM(err.ne.0,"reading coef_pccomp%lim_pc_prfl_qmin")

        ENDDO


        DO i = 1, coef_pccomp%fmv_pc_nlev
          READ (file_lu, iostat=ERR)     &
            & coef_pccomp%ref_pc_prfl_p(i), coef_pccomp%lim_pc_prfl_ozmin(i), coef_pccomp%lim_pc_prfl_ozmax(i)

          THROWM(err.ne.0,"reading coef_pccomp%lim_pc_prfl_ozmin")

        ENDDO

        READ (file_lu, iostat=ERR)coef_pccomp%lim_pc_prfl_pmin, coef_pccomp%lim_pc_prfl_pmax

        THROWM(err.ne.0,"reading coef_pccomp%lim_pc_prfl_pmin")

        READ (file_lu, iostat=ERR)coef_pccomp%lim_pc_prfl_tsmin, coef_pccomp%lim_pc_prfl_tsmax

        THROWM(err.ne.0,"reading coef_pccomp%lim_pc_prfl_tsmin")

        READ (file_lu, iostat=ERR)coef_pccomp%lim_pc_prfl_skmin, coef_pccomp%lim_pc_prfl_skmax

        THROWM(err.ne.0,"reading coef_pccomp%lim_pc_prfl_skmin")

        READ (file_lu, iostat=ERR)coef_pccomp%lim_pc_prfl_wsmin, coef_pccomp%lim_pc_prfl_wsmax

        THROWM(err.ne.0,"reading coef_pccomp%lim_pc_prfl_wsmin")

!     INSTRUMENT_NOISE
        ALLOCATE (coef_pccomp%noise_in(coef_pccomp%fmv_pc_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%noise_in")

        ALLOCATE (coef_pccomp%ff_ori_chn_in(coef_pccomp%fmv_pc_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%ff_ori_chn_in")

        ALLOCATE (coef_pccomp%ff_cwn_in(coef_pccomp%fmv_pc_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%ff_cwn_in")

        ALLOCATE (coef_pccomp%ff_bco_in(coef_pccomp%fmv_pc_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%ff_bco_in")

        ALLOCATE (coef_pccomp%ff_bcs_in(coef_pccomp%fmv_pc_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%ff_bcs_in")

        ALLOCATE (coef_pccomp%noise(coef_pccomp%fmv_pc_nchn_noise), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of coef_pccomp%noise")


        IF (all_channels_rec) THEN

          DO n = 1, coef_pccomp%fmv_pc_nchn
            READ (file_lu, iostat=ERR)coef_pccomp%ff_ori_chn_in(n), coef_pccomp%ff_cwn_in(n), coef_pccomp%ff_bco_in(n),      &
              & coef_pccomp%ff_bcs_in(n), coef_pccomp%noise_in(n)

            THROWM(err.ne.0,"reading instrument noise")

          ENDDO

          IF (all_channels) THEN
            coef_pccomp%noise = coef_pccomp%noise_in
          ELSE
            coef_pccomp%noise = coef_pccomp%noise_in(coef%ff_ori_chn(:))
          ENDIF

        ELSE
          ALLOCATE (noisearray(file_channels_rec), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of noisearray")

          ALLOCATE (pccwnarray(file_channels_rec), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of pccwnarray")

          ALLOCATE (pcchnarray(file_channels_rec), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of pcchnarray")

          ALLOCATE (pcbcoarray(file_channels_rec), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of pcbcoarray")

          ALLOCATE (pcbcsarray(file_channels_rec), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of pcbcsarray")


          DO n = 1, file_channels_rec
            READ (file_lu, iostat=ERR)pcchnarray(n), pccwnarray(n), pcbcoarray(n), pcbcsarray(n), noisearray(n)

            THROWM(err.ne.0,"reading instrument noise")

          ENDDO

          coef_pccomp%ff_ori_chn_in(:) = pcchnarray(channels_rec(:))
          coef_pccomp%ff_cwn_in(:)     = pccwnarray(channels_rec(:))
          coef_pccomp%ff_bco_in(:)     = pcbcoarray(channels_rec(:))
          coef_pccomp%ff_bcs_in(:)     = pcbcsarray(channels_rec(:))
          coef_pccomp%noise_in(:)      = noisearray(channels_rec(:))

          coef_pccomp%noise(:)         = noisearray(coef%ff_ori_chn(:))

          DEALLOCATE (noisearray, pccwnarray, pcchnarray, pcbcoarray, pcbcsarray, STAT = ERR)

          THROWM( ERR .NE. 0, "deallocation of pcbcsarray")
        ENDIF
!
! Here add reading of new sections for binary format in order to keep compatibility with
! previous versions
!

  CATCH
END SUBROUTINE rttov_read_binary_pccoef
