SUBROUTINE rttov_write_binary_pccoef( &
            & err,           &
            & coef_pccomp,   &
            & file_lu,       &
            & verbose)
! Description:
!   Write the PC coefficient structure
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
!
! Declarations:
! Modules used:
! Imported Parameters:
#include "throw.h"
! Imported Type Definitions:
  USE rttov_types, ONLY :  &
       & rttov_coef_pccomp
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE rttov_const, ONLY :  &
       & rttov_magic_string, &
       & rttov_magic_number
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(in):
  INTEGER(KIND=jpim)       , INTENT(IN)              :: file_lu      ! file logical unit number
  TYPE(rttov_coef_pccomp)  , INTENT(IN)              :: coef_pccomp
  LOGICAL(KIND=jplm)       , INTENT(IN), OPTIONAL    :: verbose      ! print out information messages or not
! scalar arguments with intent(in):
  INTEGER(KIND=jpim)       , INTENT(OUT)             :: err          ! return code
!INTF_END
#include "rttov_errorreport.interface"
  INTEGER(KIND=jpim)  :: n, m, i, j
  LOGICAL(KIND=jplm)  :: lverbose
  CHARACTER(LEN = 80) :: errMessage
!- End of header --------------------------------------------------------
  TRY
  IF (PRESENT(verbose)) THEN
    lverbose = verbose
  ELSE
    lverbose = .TRUE._jplm
  END IF

! Binary file
  IF (lverbose) THEN
    WRITE (errMessage, '( "write coefficient to file_id ", i2, " in binary format")')file_lu
    INFO(errMessage)
  END IF
! Write a string that could be displayed
! Write a real number to be able to check single/double precision
  WRITE (file_lu, iostat=err)rttov_magic_string, rttov_magic_number

  THROW(err.ne.0)

! PRINCOMP_PREDICTORS
    WRITE (file_lu, iostat=ERR)coef_pccomp%fmv_pc_comp_pc
    THROW(err.ne.0)

    WRITE (file_lu, iostat=ERR)coef_pccomp%fmv_pc_cld
    THROW(err.ne.0)

    WRITE (file_lu, iostat=ERR)coef_pccomp%fmv_pc_bands
    THROW(err.ne.0)

    WRITE (file_lu, iostat=ERR)coef_pccomp%fmv_pc_msets
    THROW(err.ne.0)

    WRITE (file_lu, iostat=ERR)coef_pccomp%fmv_pc_mnum
    THROW(err.ne.0)

    WRITE (file_lu, iostat=ERR)coef_pccomp%fmv_pc_mchn
    THROW(err.ne.0)

    DO m = 1, coef_pccomp%fmv_pc_bands

      WRITE (file_lu, iostat=ERR)coef_pccomp%fmv_pc_sets(m)
      THROW(err.ne.0)

      DO n = 1, coef_pccomp%fmv_pc_sets(m)

        WRITE (file_lu, iostat=ERR)coef_pccomp%pcreg(m,n)%fmv_pc_npred
        THROW(err.ne.0)

        WRITE (file_lu, iostat=ERR) &
         & (coef_pccomp%pcreg(m,n)%predictindex(i), i = 1, coef_pccomp%pcreg(m,n)%fmv_pc_npred)
       THROW(err.ne.0)

      ENDDO
    ENDDO

!PRINCOMP_EIGENVECTORS

    DO m = 1, coef_pccomp%fmv_pc_bands
      DO n = 1, coef_pccomp%fmv_pc_mnum

        WRITE (file_lu, iostat=ERR)(coef_pccomp%eigen(m)%eigenvectors(i, n), i = 1, coef_pccomp%fmv_pc_nchn)
        THROW(err.ne.0)

      ENDDO
    ENDDO

!PRINCOMP_COEFFICIENTS
    DO m = 1, coef_pccomp%fmv_pc_bands

      DO n = 1, coef_pccomp%fmv_pc_sets(m)

        DO j = 1, coef_pccomp%fmv_pc_mnum

          WRITE (file_lu, iostat=ERR)     &
            & (coef_pccomp%pcreg(m,n)%coefficients(i, j), i = 1, coef_pccomp%pcreg(m,n)%fmv_pc_npred)
          THROW(err.ne.0)

        ENDDO

      ENDDO
    ENDDO

!EMISSIVITY_COEFFICIENTS

    WRITE (file_lu, iostat=ERR)coef_pccomp%fmv_pc_nche
    THROW(err.ne.0)

    DO i = 1, coef_pccomp%fmv_pc_nche

      WRITE (file_lu, iostat=ERR)   &
        & coef_pccomp%emiss_chn(i), &
        & coef_pccomp%emiss_c1(i),  &
        & coef_pccomp%emiss_c2(i),  &
        & coef_pccomp%emiss_c3(i),  &
        & coef_pccomp%emiss_c4(i),  &
        & coef_pccomp%emiss_c5(i),  &
        & coef_pccomp%emiss_c6(i),  &
        & coef_pccomp%emiss_c7(i),  &
        & coef_pccomp%emiss_c8(i),  &
        & coef_pccomp%emiss_c9(i)
      THROW(err.ne.0)

    ENDDO

!PC_REFERENCE_PROFILE

    WRITE (file_lu, iostat=ERR)coef_pccomp%fmv_pc_gas
    THROW(err.ne.0)

    WRITE (file_lu, iostat=ERR)coef_pccomp%fmv_pc_nlev
    THROW(err.ne.0)

    DO n = 1, coef_pccomp%fmv_pc_gas

      DO i = 1, coef_pccomp%fmv_pc_nlev

        WRITE (file_lu, iostat=ERR)coef_pccomp%ref_pc_prfl_p(i), coef_pccomp%ref_pc_prfl_mr(i, n)
        THROW(err.ne.0)

      ENDDO

    ENDDO

!PC_PROFILE_LIMITS

    DO i = 1, coef_pccomp%fmv_pc_nlev

      WRITE (file_lu, iostat=ERR)          &
        & coef_pccomp%ref_pc_prfl_p(i),    &
        & coef_pccomp%lim_pc_prfl_tmin(i), &
        & coef_pccomp%lim_pc_prfl_tmax(i)
      THROW(err.ne.0)

    ENDDO

    DO i = 1, coef_pccomp%fmv_pc_nlev

      WRITE (file_lu, iostat=ERR)         &
        & coef_pccomp%ref_pc_prfl_p(i),   &
        & coef_pccomp%lim_pc_prfl_qmin(i),&
        & coef_pccomp%lim_pc_prfl_qmax(i)
  THROW(err.ne.0)
    ENDDO

    DO i = 1, coef_pccomp%fmv_pc_nlev

      WRITE (file_lu, iostat=ERR)           &
        & coef_pccomp%ref_pc_prfl_p(i),     &
        & coef_pccomp%lim_pc_prfl_ozmin(i), &
        & coef_pccomp%lim_pc_prfl_ozmax(i)
      THROW(err.ne.0)

    ENDDO

    WRITE (file_lu, iostat=ERR)coef_pccomp%lim_pc_prfl_pmin, coef_pccomp%lim_pc_prfl_pmax
    THROW(err.ne.0)

    WRITE (file_lu, iostat=ERR)coef_pccomp%lim_pc_prfl_tsmin, coef_pccomp%lim_pc_prfl_tsmax
    THROW(err.ne.0)

    WRITE (file_lu, iostat=ERR)coef_pccomp%lim_pc_prfl_skmin, coef_pccomp%lim_pc_prfl_skmax
    THROW(err.ne.0)

    WRITE (file_lu, iostat=ERR)coef_pccomp%lim_pc_prfl_wsmin, coef_pccomp%lim_pc_prfl_wsmax
    THROW(err.ne.0)

!INSTRUMENT_NOISE
    DO n = 1, coef_pccomp%fmv_pc_nchn

      WRITE (file_lu, iostat=ERR)       &
        & coef_pccomp%ff_ori_chn_in(n), &
        & coef_pccomp%ff_cwn_in(n),     &
        & coef_pccomp%ff_bco_in(n),     &
        & coef_pccomp%ff_bcs_in(n),     &
        & coef_pccomp%noise_in(n)
      THROW(err.ne.0)

    ENDDO

  IF (lverbose) INFO("end of write coefficient")
  CATCH
END SUBROUTINE 
