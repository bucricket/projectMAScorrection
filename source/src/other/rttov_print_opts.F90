! Description:
!> @file
!!   Print out the contents of the rttov_options structure.
!
!> @brief
!!   Print out the contents of the rttov_options structure.
!!
!! @details
!!   If not supplied the output is written to the error_unit
!!   as set by rttov_errorhandling or the default if unset.
!!
!!   The optional text argument is printed at the top of the
!!   output.
!!
!! @param[in]   opts      options to configure the simulations
!! @param[in]   lu        logical unit for output, optional
!! @param[in]   text      additional text to print, optional
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
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_print_opts(opts, lu, text)

  USE rttov_types, ONLY : rttov_options
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE rttov_global, ONLY : error_unit
!INTF_ON

  IMPLICIT NONE

  TYPE(rttov_options), INTENT(IN)           :: opts  ! options
  INTEGER(KIND=jpim),  INTENT(IN), OPTIONAL :: lu    ! logical unit for print
  CHARACTER(LEN=*),    INTENT(IN), OPTIONAL :: text  ! additional text to print
!INTF_END

  INTEGER(KIND=jpim)  :: iu       ! logical unit for print
  CHARACTER(LEN=20)   :: tmp_text ! temporary string for formatting output

  iu = error_unit
  IF (PRESENT(lu)) iu = lu

  IF (PRESENT(text)) THEN
    WRITE(iu,'(/,a,a)') "RTTOV options structure: ", TRIM(text)
  ELSE
    WRITE(iu,'(/,a)') "RTTOV options structure"
  END IF

  WRITE(iu,'(a)')           "General config options"
  WRITE(iu,'(2x,a,l1)')     "apply_reg_limits     ", opts%config%apply_reg_limits
  WRITE(iu,'(2x,a,l1)')     "verbose              ", opts%config%verbose
  WRITE(iu,'(2x,a,l1)')     "do_checkinput        ", opts%config%do_checkinput
  WRITE(iu,'(a)',advance='yes')

  WRITE(iu,'(a)')           "General RT options"
  WRITE(iu,'(2x,a,l1)')     "addrefrac            ", opts%rt_all%addrefrac
  WRITE(iu,'(2x,a,l1)')     "switchrad            ", opts%rt_all%switchrad
  WRITE(iu,'(2x,a,l1)')     "use_q2m              ", opts%rt_all%use_q2m
  WRITE(iu,'(2x,a,l1)')     "do_lambertian        ", opts%rt_all%do_lambertian
  WRITE(iu,'(a)',advance='yes')

  WRITE(iu,'(a)')           "VIS/IR-only RT options"
  WRITE(iu,'(2x,a,l1)')     "addsolar             ", opts%rt_ir%addsolar
  WRITE(iu,'(2x,a,l1)')     "do_nlte_correction   ", opts%rt_ir%do_nlte_correction
  WRITE(iu,'(2x,a,l1)')     "addaerosl            ", opts%rt_ir%addaerosl
  WRITE(iu,'(2x,a,l1)')     "user_aer_opt_param   ", opts%rt_ir%user_aer_opt_param
  WRITE(iu,'(2x,a,l1)')     "addclouds            ", opts%rt_ir%addclouds
  WRITE(iu,'(2x,a,l1)')     "user_cld_opt_param   ", opts%rt_ir%user_cld_opt_param
  WRITE(tmp_text,'(f10.4)') opts%rt_ir%cldstr_threshold
  WRITE(iu,'(2x,a,a)')      "cldstr_threshold     ", TRIM(ADJUSTL(tmp_text))
  WRITE(iu,'(2x,a,l1)')     "cldstr_simple        ", opts%rt_ir%cldstr_simple
  WRITE(iu,'(2x,a,l1)')     "do_lambertian        ", opts%rt_ir%do_lambertian
  WRITE(iu,'(2x,a,l1)')     "ozone_data           ", opts%rt_ir%ozone_data
  WRITE(iu,'(2x,a,l1)')     "co2_data             ", opts%rt_ir%co2_data
  WRITE(iu,'(2x,a,l1)')     "n2o_data             ", opts%rt_ir%n2o_data
  WRITE(iu,'(2x,a,l1)')     "co_data              ", opts%rt_ir%co_data
  WRITE(iu,'(2x,a,l1)')     "ch4_data             ", opts%rt_ir%ch4_data
  WRITE(iu,'(a)',advance='yes')

  WRITE(iu,'(a)')           "PC-RTTOV options"
  WRITE(iu,'(2x,a,l1)')     "addpc                ", opts%rt_ir%pc%addpc
  WRITE(tmp_text,'(i3)') opts%rt_ir%pc%ipcbnd
  WRITE(iu,'(2x,a,a)')      "ipcbnd               ", TRIM(ADJUSTL(tmp_text))
  WRITE(tmp_text,'(i3)') opts%rt_ir%pc%ipcreg
  WRITE(iu,'(2x,a,a)')      "ipcreg               ", TRIM(ADJUSTL(tmp_text))
  WRITE(iu,'(2x,a,l1)')     "addradrec            ", opts%rt_ir%pc%addradrec
  WRITE(iu,'(a)',advance='yes')

  WRITE(iu,'(a)')           "MW-only RT options"
  WRITE(tmp_text,'(i3)') opts%rt_mw%fastem_version
  WRITE(iu,'(2x,a,a)')      "fastem_version       ", TRIM(ADJUSTL(tmp_text))
  WRITE(iu,'(2x,a,l1)')     "clw_data             ", opts%rt_mw%clw_data
  WRITE(iu,'(2x,a,l1)')     "do_lambertian        ", opts%rt_mw%do_lambertian
  WRITE(iu,'(2x,a,l1)')     "supply_foam_fraction ", opts%rt_mw%supply_foam_fraction
  WRITE(iu,'(a)',advance='yes')

  WRITE(iu,'(a)')           "Interpolation and vertical grid options"
  WRITE(iu,'(2x,a,l1)')     "addinterp            ", opts%interpolation%addinterp
  WRITE(tmp_text,'(i3)') opts%interpolation%interp_mode
  WRITE(iu,'(2x,a,a)')      "interp_mode          ", TRIM(ADJUSTL(tmp_text))
  WRITE(iu,'(2x,a,l1)')     "spacetop             ", opts%interpolation%spacetop
  WRITE(iu,'(2x,a,l1)')     "lgradp               ", opts%interpolation%lgradp
  WRITE(iu,'(2x,a,l1)')     "reg_limit_extrap     ", opts%interpolation%reg_limit_extrap
  WRITE(iu,'(a)',advance='yes')

END SUBROUTINE rttov_print_opts
