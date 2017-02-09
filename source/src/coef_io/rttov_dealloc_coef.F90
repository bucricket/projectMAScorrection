!
SUBROUTINE rttov_dealloc_coef(err, coef)
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
!  1.9       10/01/2013  Add PMC shifts (P Rayer)
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
       & rttov_coef
  USE parkind1, ONLY : jpim
!INTF_OFF
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(out):
  INTEGER(KIND=jpim), INTENT(OUT)   :: err          ! return code
  TYPE(rttov_coef),   INTENT(INOUT) :: coef         ! coefficients
!INTF_END
#include "rttov_nullify_coef.interface"
#include "rttov_errorreport.interface"
! Local Arrays and Scalars:
!- End of header --------------------------------------------------------

  TRY

  IF (ASSOCIATED(coef%fmv_gas_id)) DEALLOCATE (coef%fmv_gas_id, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%fmv_gas_pos)) DEALLOCATE (coef%fmv_gas_pos, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%fmv_var)) DEALLOCATE (coef%fmv_var, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%fmv_coe)) DEALLOCATE (coef%fmv_coe, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%fmv_lvl)) DEALLOCATE (coef%fmv_lvl, STAT = err)
  THROW( err .NE. 0 )


  IF (ASSOCIATED(coef%ff_ori_chn)) DEALLOCATE (coef%ff_ori_chn, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%ff_val_chn)) DEALLOCATE (coef%ff_val_chn, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%ff_cwn)) DEALLOCATE (coef%ff_cwn, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%ff_bco)) DEALLOCATE (coef%ff_bco, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%ff_bcs)) DEALLOCATE (coef%ff_bcs, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%ff_gam)) DEALLOCATE (coef%ff_gam, STAT = err)
  THROW( err .NE. 0 )


  IF (ASSOCIATED(coef%gaz_units)) DEALLOCATE (coef%gaz_units, STAT = err)
  THROW( err .NE. 0 )


  IF (ASSOCIATED(coef%fastem_polar)) DEALLOCATE (coef%fastem_polar, STAT = err)
  THROW( err .NE. 0 )


  IF (ASSOCIATED(coef%ssirem_chn)) DEALLOCATE (coef%ssirem_chn, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%ssirem_a0)) DEALLOCATE (coef%ssirem_a0, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%ssirem_a1)) DEALLOCATE (coef%ssirem_a1, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%ssirem_a2)) DEALLOCATE (coef%ssirem_a2, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%ssirem_xzn1)) DEALLOCATE (coef%ssirem_xzn1, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%ssirem_xzn2)) DEALLOCATE (coef%ssirem_xzn2, STAT = err)
  THROW( err .NE. 0 )


  IF (ASSOCIATED(coef%ref_prfl_p)) DEALLOCATE (coef%ref_prfl_p, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%ref_prfl_t)) DEALLOCATE (coef%ref_prfl_t, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%ref_prfl_mr)) DEALLOCATE (coef%ref_prfl_mr, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%lim_prfl_p)) DEALLOCATE (coef%lim_prfl_p, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%lim_prfl_tmax)) DEALLOCATE (coef%lim_prfl_tmax, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%lim_prfl_tmin)) DEALLOCATE (coef%lim_prfl_tmin, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%lim_prfl_gmin)) DEALLOCATE (coef%lim_prfl_gmin, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%lim_prfl_gmax)) DEALLOCATE (coef%lim_prfl_gmax, STAT = err)
  THROW( err .NE. 0 )


  IF (ASSOCIATED(coef%thermal)) THEN

    IF (ASSOCIATED(coef%thermal%mixedgas)) DEALLOCATE (coef%thermal%mixedgas, STAT = err)
    THROW( err .NE. 0 )

    IF (ASSOCIATED(coef%thermal%watervapour)) DEALLOCATE (coef%thermal%watervapour, STAT = err)
    THROW( err .NE. 0 )

    IF (ASSOCIATED(coef%thermal%ozone)) DEALLOCATE (coef%thermal%ozone, STAT = err)
    THROW( err .NE. 0 )

    IF (ASSOCIATED(coef%thermal%wvcont)) DEALLOCATE (coef%thermal%wvcont, STAT = err)
    THROW( err .NE. 0 )

    IF (ASSOCIATED(coef%thermal%co2)) DEALLOCATE (coef%thermal%co2, STAT = err)
    THROW( err .NE. 0 )

    IF (ASSOCIATED(coef%thermal%n2o)) DEALLOCATE (coef%thermal%n2o, STAT = err)
    THROW( err .NE. 0 )

    IF (ASSOCIATED(coef%thermal%co)) DEALLOCATE (coef%thermal%co, STAT = err)
    THROW( err .NE. 0 )

    IF (ASSOCIATED(coef%thermal%ch4)) DEALLOCATE (coef%thermal%ch4, STAT = err)
    THROW( err .NE. 0 )

    IF (ASSOCIATED(coef%thermal)) DEALLOCATE(coef%thermal, STAT = err)
    THROW( err .NE. 0 )

  ENDIF

  IF (coef%solarcoef) THEN

    IF (ASSOCIATED(coef%solar%mixedgas)) DEALLOCATE (coef%solar%mixedgas, STAT = err)
    THROW( err .NE. 0 )

    IF (ASSOCIATED(coef%solar%watervapour)) DEALLOCATE (coef%solar%watervapour, STAT = err)
    THROW( err .NE. 0 )

    IF (ASSOCIATED(coef%solar%ozone)) DEALLOCATE (coef%solar%ozone, STAT = err)
    THROW( err .NE. 0 )

    IF (ASSOCIATED(coef%solar%wvcont)) DEALLOCATE (coef%solar%wvcont, STAT = err)
    THROW( err .NE. 0 )

    IF (ASSOCIATED(coef%solar%co2)) DEALLOCATE (coef%solar%co2, STAT = err)
    THROW( err .NE. 0 )

    IF (ASSOCIATED(coef%solar%n2o)) DEALLOCATE (coef%solar%n2o, STAT = err)
    THROW( err .NE. 0 )

    IF (ASSOCIATED(coef%solar%co)) DEALLOCATE (coef%solar%co, STAT = err)
    THROW( err .NE. 0 )

    IF (ASSOCIATED(coef%solar%ch4)) DEALLOCATE (coef%solar%ch4, STAT = err)
    THROW( err .NE. 0 )

    IF (ASSOCIATED(coef%solar)) DEALLOCATE (coef%solar, STAT = err)
    THROW( err .NE. 0 )

  ENDIF

  IF (ASSOCIATED(coef%mixedgasint)) DEALLOCATE (coef%mixedgasint, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%watervapourint)) DEALLOCATE (coef%watervapourint, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%ozoneint)) DEALLOCATE (coef%ozoneint, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%wvcontint)) DEALLOCATE (coef%wvcontint, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%co2int)) DEALLOCATE (coef%co2int, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%n2oint)) DEALLOCATE (coef%n2oint, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%coint)) DEALLOCATE (coef%coint, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%ch4int)) DEALLOCATE (coef%ch4int, STAT = err)
  THROW( err .NE. 0 )


! NLTE coefficients
  IF (coef%nltecoef) THEN
    IF (ASSOCIATED(coef%nlte_coef%coef)) DEALLOCATE(coef%nlte_coef%coef, STAT = err)
    THROW( err .NE. 0 )

    IF (ASSOCIATED(coef%nlte_coef%sol_zen_angle)) DEALLOCATE(coef%nlte_coef%sol_zen_angle, STAT = err)
    THROW( err .NE. 0 )

    IF (ASSOCIATED(coef%nlte_coef%sat_zen_angle)) DEALLOCATE(coef%nlte_coef%sat_zen_angle, STAT = err)
    THROW( err .NE. 0 )

    IF (ASSOCIATED(coef%nlte_coef%cos_sol)) DEALLOCATE(coef%nlte_coef%cos_sol, STAT = err)
    THROW( err .NE. 0 )

    IF (ASSOCIATED(coef%nlte_coef%sec_sat)) DEALLOCATE(coef%nlte_coef%sec_sat, STAT = err)
    THROW( err .NE. 0 )

    IF (ASSOCIATED(coef%nlte_coef)) DEALLOCATE(coef%nlte_coef, STAT = err)
    THROW( err .NE. 0 )
  ENDIF

! PMC shift coefs
  IF (coef%pmc_shift) THEN
    IF (ASSOCIATED(coef%pmc_ppmc)) DEALLOCATE(coef%pmc_ppmc, STAT = err)
    THROW( err .NE. 0 )

    IF (ASSOCIATED(coef%pmc_coef)) DEALLOCATE(coef%pmc_coef, STAT = err)
    THROW( err .NE. 0 )

    IF (ASSOCIATED(coef%pmc_pnominal)) DEALLOCATE(coef%pmc_pnominal, STAT = err)
    THROW( err .NE. 0 )
  ENDIF

! planck variables
  IF (ASSOCIATED(coef%planck1)) DEALLOCATE (coef%planck1, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%planck2)) DEALLOCATE (coef%planck2, STAT = err)
  THROW( err .NE. 0 )

! frequency in GHz for MicroWaves

  IF (ASSOCIATED(coef%frequency_ghz)) DEALLOCATE (coef%frequency_ghz, STAT = err)
  THROW( err .NE. 0 )


  IF (ASSOCIATED(coef%dp)) DEALLOCATE (coef%dp, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%dpp)) DEALLOCATE (coef%dpp, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%tstar)) DEALLOCATE (coef%tstar, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%to3star)) DEALLOCATE (coef%to3star, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%wstar)) DEALLOCATE (coef%wstar, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%ostar)) DEALLOCATE (coef%ostar, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%co2star)) DEALLOCATE (coef%co2star, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%co2star)) DEALLOCATE (coef%co2star, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%n2ostar)) DEALLOCATE (coef%n2ostar, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%costar)) DEALLOCATE (coef%costar, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%ch4star)) DEALLOCATE (coef%ch4star, STAT = err)
  THROW( err .NE. 0 )


  IF (ASSOCIATED(coef%tt_chn)) DEALLOCATE (coef%tt_chn, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%tt_val_chn)) DEALLOCATE (coef%tt_val_chn, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%tt_cwn)) DEALLOCATE (coef%tt_cwn, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%tt_a0)) DEALLOCATE (coef%tt_a0, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%tt_a1)) DEALLOCATE (coef%tt_a1, STAT = err)
  THROW( err .NE. 0 )


  IF (ASSOCIATED(coef%pw_chn)) DEALLOCATE (coef%pw_chn, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%pw_val_chn)) DEALLOCATE (coef%pw_val_chn, STAT = err)
  THROW( err .NE. 0 )


  IF (ASSOCIATED(coef%ss_chn)) DEALLOCATE (coef%ss_chn, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%ss_val_chn)) DEALLOCATE (coef%ss_val_chn, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%ss_cwn)) DEALLOCATE (coef%ss_cwn, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%ss_solar_spectrum)) DEALLOCATE (coef%ss_solar_spectrum, STAT = err)
  THROW( err .NE. 0 )


  IF (ASSOCIATED(coef%woc_chn)) DEALLOCATE (coef%woc_chn, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%woc_cwn)) DEALLOCATE (coef%woc_cwn, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%woc_waopc_ow)) DEALLOCATE (coef%woc_waopc_ow, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%woc_waopc_fw)) DEALLOCATE (coef%woc_waopc_fw, STAT = err)
  THROW( err .NE. 0 )


  IF (ASSOCIATED(coef%ws_k_omega)) DEALLOCATE (coef%ws_k_omega, STAT = err)
  THROW( err .NE. 0 )

  IF (ASSOCIATED(coef%ws_npoint)) DEALLOCATE (coef%ws_npoint, STAT = err)
  THROW( err .NE. 0 )


  CALL rttov_nullify_coef(coef)

  CATCH
END SUBROUTINE rttov_dealloc_coef
