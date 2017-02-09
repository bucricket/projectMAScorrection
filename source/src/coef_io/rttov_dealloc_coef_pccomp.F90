!
SUBROUTINE rttov_dealloc_coef_pccomp (ERR, coef_pccomp)
! Description:
! de-allocation of a coefficient structure
! The allocation is done by the readcoef subroutine called by the user
! this subroutine should be called once per coef structure when
! all rttov calls are completed.
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
!  1.1       03/05/2004  Add specific RTTOV8 CO2 variable (P. Brunel)
!  1.2       02/06/2004  Update for RTTOV8 coefS (P. Brunel)
!  1.3       01/06/2005  Marco Matricardi (ECMWF):
!               --       N2O,CO and CH4 variables added
!  1.4       26/04/2007  Cloud/aerosol variables added (R Saunders)
!  1.5       11/10/2007  Nullify unused coef pointers P.Marguinaud
!  1.6       23/94/2008  Add some initialisation of scalars (P. Brunel)
!  1.7       06/03/2009  Deafault now coef % IncTop = .true. (P.Rayer)
!  1.8       02/12/2009  Add principal components (M. Matricardi, ECMWF)
!
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
! Imported Type Definitions:
#include "throw.h"
  USE rttov_types, ONLY :  &
       & rttov_coef_pccomp
  USE parkind1, ONLY : jpim
!INTF_OFF
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(out):
  INTEGER(KIND=jpim)       , INTENT(OUT)   :: ERR          ! return code
  TYPE(rttov_coef_pccomp  ), INTENT(INOUT) :: coef_pccomp  ! coefficients
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_nullify_coef_pccomp.interface"
! Local Arrays and Scalars:
  INTEGER(KIND=jpim) :: m, n
!- End of header --------------------------------------------------------
  TRY

    IF (Associated (coef_pccomp%emiss_chn   )) &
        DEALLOCATE (coef_pccomp%emiss_chn, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated (coef_pccomp%emiss_c1    )) &
        DEALLOCATE (coef_pccomp%emiss_c1, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated (coef_pccomp%emiss_c2    )) &
        DEALLOCATE (coef_pccomp%emiss_c2, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated (coef_pccomp%emiss_c3    )) &
        DEALLOCATE (coef_pccomp%emiss_c3, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated (coef_pccomp%emiss_c4    )) &
        DEALLOCATE (coef_pccomp%emiss_c4, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated (coef_pccomp%emiss_c5    )) &
        DEALLOCATE (coef_pccomp%emiss_c5, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated (coef_pccomp%emiss_c6    )) &
        DEALLOCATE (coef_pccomp%emiss_c6, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated (coef_pccomp%emiss_c7    )) &
        DEALLOCATE (coef_pccomp%emiss_c7, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated (coef_pccomp%emiss_c8    )) &
        DEALLOCATE (coef_pccomp%emiss_c8, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated (coef_pccomp%emiss_c9    )) &
        DEALLOCATE (coef_pccomp%emiss_c9, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated (coef_pccomp%ff_cwn_in   )) &
        DEALLOCATE (coef_pccomp%ff_cwn_in, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated (coef_pccomp%ff_bco_in   )) &
        DEALLOCATE (coef_pccomp%ff_bco_in, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated (coef_pccomp%ff_bcs_in   )) &
        DEALLOCATE (coef_pccomp%ff_bcs_in, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated (coef_pccomp%planck1_in  )) &
        DEALLOCATE (coef_pccomp%planck1_in, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated (coef_pccomp%planck2_in  )) &
        DEALLOCATE (coef_pccomp%planck2_in, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_pccomp%ff_ori_chn_in)) DEALLOCATE (coef_pccomp%ff_ori_chn_in, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_pccomp%noise_in)) DEALLOCATE (coef_pccomp%noise_in, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_pccomp%noise)) DEALLOCATE (coef_pccomp%noise, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_pccomp%noise_r)) DEALLOCATE (coef_pccomp%noise_r, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_pccomp%ref_pc_prfl_mr)) DEALLOCATE (coef_pccomp%ref_pc_prfl_mr, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_pccomp%ref_pc_prfl_p)) DEALLOCATE (coef_pccomp%ref_pc_prfl_p, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_pccomp%lim_pc_prfl_tmin)) DEALLOCATE (coef_pccomp%lim_pc_prfl_tmin, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_pccomp%lim_pc_prfl_tmax)) DEALLOCATE (coef_pccomp%lim_pc_prfl_tmax, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_pccomp%lim_pc_prfl_qmin)) DEALLOCATE (coef_pccomp%lim_pc_prfl_qmin, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_pccomp%lim_pc_prfl_qmax)) DEALLOCATE (coef_pccomp%lim_pc_prfl_qmax, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_pccomp%lim_pc_prfl_ozmin)) DEALLOCATE (coef_pccomp%lim_pc_prfl_ozmin, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_pccomp%lim_pc_prfl_ozmax)) DEALLOCATE (coef_pccomp%lim_pc_prfl_ozmax, STAT = ERR)

    THROW( ERR .NE. 0 )


    IF (Associated(coef_pccomp%co2_pc_ref)) &
      DEALLOCATE (coef_pccomp%co2_pc_ref, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_pccomp%co_pc_ref)) &
      DEALLOCATE (coef_pccomp%co_pc_ref, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_pccomp%ch4_pc_ref)) &
      DEALLOCATE (coef_pccomp%ch4_pc_ref, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_pccomp%n2o_pc_ref)) &
      DEALLOCATE (coef_pccomp%n2o_pc_ref, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_pccomp%eigen)) THEN
      DO m = 1, coef_pccomp%fmv_pc_bands
        IF (Associated(coef_pccomp%eigen(m)%eigenvectors)) &
          DEALLOCATE (coef_pccomp%eigen(m)%eigenvectors, STAT = ERR)
        IF (ASSOCIATED(coef_pccomp%eigen(m)%eigenvectors_t)) &
          DEALLOCATE (coef_pccomp%eigen(m)%eigenvectors_t, STAT = ERR)

        THROW( ERR .NE. 0 )
      ENDDO

      DEALLOCATE(coef_pccomp%eigen)

      THROW( ERR .NE. 0 )

    ENDIF

    IF (Associated(coef_pccomp%pcreg)) THEN

      DO m = 1, coef_pccomp%fmv_pc_bands
        DO n = 1, coef_pccomp%fmv_pc_sets(m)
          IF (Associated(coef_pccomp%pcreg(m,n)%coefficients)) &
            DEALLOCATE (coef_pccomp%pcreg(m,n)%coefficients, STAT = ERR)
          IF (ASSOCIATED(coef_pccomp%pcreg(m,n)%coefficients_t)) &
            DEALLOCATE (coef_pccomp%pcreg(m,n)%coefficients_t, STAT = ERR)

          THROW( ERR .NE. 0 )

          IF (Associated(coef_pccomp%pcreg(m,n)%predictindex)) &
            DEALLOCATE (coef_pccomp%pcreg(m,n)%predictindex, STAT = ERR)

          THROW( ERR .NE. 0 )

        ENDDO
      ENDDO

      DEALLOCATE (coef_pccomp%pcreg, STAT = ERR)

      THROW( ERR .NE. 0 )

    ENDIF

    IF (Associated(coef_pccomp%fmv_pc_sets)) &
      DEALLOCATE(coef_pccomp%fmv_pc_sets)


  CALL rttov_nullify_coef_pccomp (coef_pccomp)

  CATCH
END SUBROUTINE
