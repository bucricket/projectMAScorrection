FUNCTION rttov_opts_eq(opts1, opts2)
! Description:
!   Check equality of options structures
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

  USE parkind1, ONLY : jplm
  USE rttov_types, ONLY : rttov_options

  IMPLICIT NONE

  LOGICAL(KIND=jplm) :: rttov_opts_eq

  TYPE(rttov_options), INTENT(IN) :: opts1
  TYPE(rttov_options), INTENT(IN) :: opts2
!INTF_END

! Only options meaningful for calculations
  rttov_opts_eq = &
   ( opts1%config%apply_reg_limits        .EQV. opts2%config%apply_reg_limits        ) .AND. &
   ( opts1%config%do_checkinput           .EQV. opts2%config%do_checkinput           ) .AND. &
   ( opts1%rt_all%switchrad               .EQV. opts2%rt_all%switchrad               ) .AND. &
   ( opts1%rt_all%addrefrac               .EQV. opts2%rt_all%addrefrac               ) .AND. &
   ( opts1%rt_all%use_q2m                 .EQV. opts2%rt_all%use_q2m                 ) .AND. &
   ( opts1%rt_all%do_lambertian           .EQV. opts2%rt_all%do_lambertian           ) .AND. &
   ( opts1%rt_ir%addsolar                 .EQV. opts2%rt_ir%addsolar                 ) .AND. &
   ( opts1%rt_ir%do_nlte_correction       .EQV. opts2%rt_ir%do_nlte_correction       ) .AND. &
   ( opts1%rt_ir%addaerosl                .EQV. opts2%rt_ir%addaerosl                ) .AND. &
   ( opts1%rt_ir%addclouds                .EQV. opts2%rt_ir%addclouds                ) .AND. &
   ( opts1%rt_ir%user_aer_opt_param       .EQV. opts2%rt_ir%user_aer_opt_param       ) .AND. &
   ( opts1%rt_ir%user_cld_opt_param       .EQV. opts2%rt_ir%user_cld_opt_param       ) .AND. &
   ( opts1%rt_ir%cldstr_threshold          ==   opts2%rt_ir%cldstr_threshold         ) .AND. &
   ( opts1%rt_ir%cldstr_simple            .EQV. opts2%rt_ir%cldstr_simple            ) .AND. &
   ( opts1%rt_ir%do_lambertian            .EQV. opts2%rt_ir%do_lambertian            ) .AND. &
   ( opts1%rt_ir%ozone_data               .EQV. opts2%rt_ir%ozone_data               ) .AND. &
   ( opts1%rt_ir%co2_data                 .EQV. opts2%rt_ir%co2_data                 ) .AND. &
   ( opts1%rt_ir%n2o_data                 .EQV. opts2%rt_ir%n2o_data                 ) .AND. &
   ( opts1%rt_ir%co_data                  .EQV. opts2%rt_ir%co_data                  ) .AND. &
   ( opts1%rt_ir%ch4_data                 .EQV. opts2%rt_ir%ch4_data                 ) .AND. &
   ( opts1%rt_ir%pc%addpc                 .EQV. opts2%rt_ir%pc%addpc                 ) .AND. &
   ( opts1%rt_ir%pc%ipcbnd                 ==   opts2%rt_ir%pc%ipcbnd                ) .AND. &
   ( opts1%rt_ir%pc%ipcreg                 ==   opts2%rt_ir%pc%ipcreg                ) .AND. &
   ( opts1%rt_ir%pc%addradrec             .EQV. opts2%rt_ir%pc%addradrec             ) .AND. &
   ( opts1%rt_mw%fastem_version            ==   opts2%rt_mw%fastem_version           ) .AND. &
   ( opts1%rt_mw%clw_data                 .EQV. opts2%rt_mw%clw_data                 ) .AND. &
   ( opts1%rt_mw%do_lambertian            .EQV. opts2%rt_mw%do_lambertian            ) .AND. &
   ( opts1%rt_mw%supply_foam_fraction     .EQV. opts2%rt_mw%supply_foam_fraction     ) .AND. &
   ( opts1%interpolation%addinterp        .EQV. opts2%interpolation%addinterp        ) .AND. &
   ( opts1%interpolation%interp_mode       ==   opts2%interpolation%interp_mode      ) .AND. &
   ( opts1%interpolation%reg_limit_extrap .EQV. opts2%interpolation%reg_limit_extrap ) .AND. &
   ( opts1%interpolation%spacetop         .EQV. opts2%interpolation%spacetop         ) .AND. &
   ( opts1%interpolation%lgradp           .EQV. opts2%interpolation%lgradp           )

END FUNCTION
