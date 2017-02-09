SUBROUTINE rttov_nullify_coef(coef)
! Description:
!   Nullify the coefficient structure
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
!INTF_OFF
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  USE rttov_types, ONLY : rttov_coef
  IMPLICIT NONE
  TYPE(rttov_coef), INTENT(INOUT) :: coef
!INTF_END
  coef%id_platform      = 0_JPIM
  coef%id_sat           = 0_JPIM
  coef%id_inst          = 0_JPIM
  coef%id_sensor        = 0_JPIM
  coef%id_comp_lvl      = 0_JPIM
  coef%id_creation_date = (/0_JPIM, 0_JPIM, 0_JPIM/)
  coef%line_by_line     = 'xxxx'
  coef%readme_srf       = 'xxxx'
  coef%id_creation      = 'xxxx'
  coef%id_Common_name   = 'xxxx'
  coef%fmv_model_def    = 'xxxx'
  coef%fmv_model_ver    = 0_JPIM
  coef%id_comp_pc       = 0_JPIM
  coef%inczeeman        = .FALSE.
  coef%fmv_ori_nchn     = 0_JPIM
  coef%fmv_chn          = 0_JPIM
  coef%fmv_gas          = 0_JPIM
  coef%nmixed           = 0_JPIM
  coef%nwater           = 0_JPIM
  coef%nozone           = 0_JPIM
  coef%nwvcont          = 0_JPIM
  coef%nco2             = 0_JPIM
  coef%nn2o             = 0_JPIM
  coef%nco              = 0_JPIM
  coef%nch4             = 0_JPIM
  coef%nlevels          = 0_JPIM
  coef%nlayers          = 0_JPIM
  coef%ncmixed          = 0_JPIM
  coef%ncwater          = 0_JPIM
  coef%ncozone          = 0_JPIM
  coef%ncwvcont         = 0_JPIM
  coef%ncco2            = 0_JPIM
  coef%ncn2o            = 0_JPIM
  coef%ncco             = 0_JPIM
  coef%ncch4            = 0_JPIM
  coef%nintmixed        = 0_JPIM
  coef%nintwater        = 0_JPIM
  coef%nintozone        = 0_JPIM
  coef%nintwvcont       = 0_JPIM
  coef%nintco2          = 0_JPIM
  coef%nintn2o          = 0_JPIM
  coef%nintco           = 0_JPIM
  coef%nintch4          = 0_JPIM
  coef%ws_nomega        = 1_JPIM                    ! used in automatic array declaration
  coef%fc_planck_c1     = 0._JPRB
  coef%fc_planck_c2     = 0._JPRB
  coef%fc_sat_height    = 0._JPRB
  coef%fastem_ver       = 0_JPIM
  coef%ssirem_ver       = 0_JPIM
  coef%ratoe            = 0._JPRB
  coef%mwcldtop         = 0_JPIM
  coef%solarcoef        = .FALSE.                   ! Solar coefs present
  coef%nltecoef         = .FALSE.
  coef%pmc_shift        = .FALSE.
  coef%pmc_lengthcell   = 0._JPRB
  coef%pmc_tempcell     = 0._JPRB
  coef%pmc_betaplus1    = 0._JPRB
  coef%pmc_nlay         = 0_JPIM
  coef%pmc_nvar         = 0_JPIM
  NULLIFY (coef%fmv_gas_id)
  NULLIFY (coef%fmv_gas_pos)
  NULLIFY (coef%fmv_var)
  NULLIFY (coef%fmv_coe)
  NULLIFY (coef%fmv_int)
  NULLIFY (coef%fmv_lvl)
  NULLIFY (coef%gaz_units)
  NULLIFY (coef%ff_ori_chn)
  NULLIFY (coef%ff_val_chn)
  NULLIFY (coef%ff_cwn)
  NULLIFY (coef%ff_bco)
  NULLIFY (coef%ff_bcs)
  NULLIFY (coef%ff_gam)
  NULLIFY (coef%tt_chn)
  NULLIFY (coef%tt_val_chn)
  NULLIFY (coef%tt_cwn)
  NULLIFY (coef%tt_a0)
  NULLIFY (coef%tt_a1)
  NULLIFY (coef%pw_chn)
  NULLIFY (coef%pw_val_chn)
  NULLIFY (coef%ss_chn)
  NULLIFY (coef%ss_val_chn)
  NULLIFY (coef%ss_cwn)
  NULLIFY (coef%ss_solar_spectrum)
  NULLIFY (coef%woc_chn)
  NULLIFY (coef%woc_cwn)
  NULLIFY (coef%woc_waopc_ow)
  NULLIFY (coef%woc_waopc_fw)
  NULLIFY (coef%ws_npoint)
  NULLIFY (coef%ws_k_omega)
  NULLIFY (coef%fastem_polar)
  NULLIFY (coef%ssirem_chn)
  NULLIFY (coef%ssirem_a0)
  NULLIFY (coef%ssirem_a1)
  NULLIFY (coef%ssirem_a2)
  NULLIFY (coef%ssirem_xzn1)
  NULLIFY (coef%ssirem_xzn2)
  NULLIFY (coef%ref_prfl_p)
  NULLIFY (coef%ref_prfl_t)
  NULLIFY (coef%ref_prfl_mr)
  NULLIFY (coef%lim_prfl_p)
  NULLIFY (coef%lim_prfl_tmax)
  NULLIFY (coef%lim_prfl_tmin)
  NULLIFY (coef%lim_prfl_gmax)
  NULLIFY (coef%lim_prfl_gmin)
  NULLIFY (coef%thermal)
  NULLIFY (coef%solar)
  NULLIFY (coef%mixedgasint)
  NULLIFY (coef%watervapourint)
  NULLIFY (coef%ozoneint)
  NULLIFY (coef%wvcontint)
  NULLIFY (coef%co2int)
  NULLIFY (coef%n2oint)
  NULLIFY (coef%coint)
  NULLIFY (coef%ch4int)
  NULLIFY (coef%planck1)
  NULLIFY (coef%planck2)
  NULLIFY (coef%frequency_ghz)
  NULLIFY (coef%dp)
  NULLIFY (coef%dpp)
  NULLIFY (coef%tstar)
  NULLIFY (coef%to3star)
  NULLIFY (coef%wstar)
  NULLIFY (coef%ostar)
  NULLIFY (coef%co2star)
  NULLIFY (coef%n2ostar)
  NULLIFY (coef%costar)
  NULLIFY (coef%ch4star)
  NULLIFY (coef%nlte_coef)
  NULLIFY (coef%pmc_coef)
  NULLIFY (coef%pmc_pnominal)
  NULLIFY (coef%pmc_ppmc)
END SUBROUTINE 
