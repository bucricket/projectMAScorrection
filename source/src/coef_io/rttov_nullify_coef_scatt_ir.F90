SUBROUTINE rttov_nullify_coef_scatt_ir(coef_scatt_ir)
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
!INTF_OFF
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  USE rttov_types, ONLY : rttov_coef_scatt_ir
  IMPLICIT NONE
  TYPE(rttov_coef_scatt_ir), INTENT(INOUT) :: coef_scatt_ir
!INTF_END
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
