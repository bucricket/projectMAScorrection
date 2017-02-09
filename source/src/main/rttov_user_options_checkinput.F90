! Description:
!> @file
!!   Check consistency of selected options and coefficient file.
!
!> @brief
!!   Check consistency of selected options and coefficient file,
!!   useful for "debugging" simulations.
!!
!! @param[out]    err                  status on exit
!! @param[in]     opts                 options to configure the simulations
!! @param[in]     coefs                coefficients structure for instrument to simulate
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
SUBROUTINE rttov_user_options_checkinput( &
       & err,                &
       & opts,               &
       & coefs               )

!INTF_OFF
#include "throw.h"
!INTF_ON

  USE rttov_types, ONLY : &
      rttov_coefs, &
      rttov_options

  USE parkind1, ONLY : jpim
!INTF_OFF
  USE yomhook, ONLY : &
      LHOOK,&
      DR_HOOK

  USE parkind1, ONLY : &
      jprb

  USE rttov_const, ONLY : &
      sensor_id_ir, &
      sensor_id_hi, &
      sensor_id_mw, &
      sensor_id_po, &
      ninterp_modes
!INTF_ON

  IMPLICIT NONE

  INTEGER(KIND=jpim),  INTENT(OUT) :: err   ! status on exit
  TYPE(rttov_options), INTENT(IN)  :: opts  ! options to configure the simulations
  TYPE(rttov_coefs),   INTENT(IN)  :: coefs ! coefficients structure for instrument to simulate
!INTF_END

#include "rttov_errorreport.interface"

  CHARACTER(LEN=256) :: msg
  REAL(KIND=jprb) :: ZHOOK_HANDLE

!- End of header --------------------------------------------------------

  TRY

  IF (LHOOK) CALL DR_HOOK('RTTOV_USER_OPTIONS_CHECKINPUT',0_jpim,ZHOOK_HANDLE)

  IF (opts%interpolation%addinterp) THEN
    IF (opts%interpolation%interp_mode < 1 .OR. &
        opts%interpolation%interp_mode > ninterp_modes) THEN
      err = errorstatus_fatal
      msg = 'Error in specified interpolation mode'
      THROWM(err .NE. 0, msg)
    ENDIF
  ELSE
    IF (opts%interpolation%lgradp) THEN
      err = errorstatus_fatal
      msg = 'Interpolation should be enabled if LGRADP is set to TRUE'
      THROWM(err .NE. 0, msg)
    ENDIF
  ENDIF

  IF (opts%rt_ir%ozone_data .AND. coefs%coef%nozone == 0_jpim) THEN
    err = errorstatus_fatal
    msg = 'Coefficient file does not allow variable O3'
    THROWM(err .NE. 0, msg)
  ENDIF

  IF (opts%rt_ir%co2_data .AND. coefs%coef%nco2 == 0_jpim) THEN
    err = errorstatus_fatal
    msg = 'Coefficient file does not allow variable CO2'
    THROWM(err .NE. 0, msg)
  ENDIF

  IF (opts%rt_ir%co_data .AND. coefs%coef%nco == 0_jpim) THEN
    err = errorstatus_fatal
    msg = 'Coefficient file does not allow variable CO'
    THROWM(err .NE. 0, msg)
  ENDIF

  IF (opts%rt_ir%n2o_data .AND. coefs%coef%nn2o == 0_jpim) THEN
    err = errorstatus_fatal
    msg = 'Coefficient file does not allow variable N2O'
    THROWM(err .NE. 0, msg)
  ENDIF

  IF (opts%rt_ir%ch4_data .AND. coefs%coef%nch4 == 0_jpim) THEN
    err = errorstatus_fatal
    msg = 'Coefficient file does not allow variable CH4'
    THROWM(err .NE. 0, msg)
  ENDIF


  IF (opts%rt_ir%addsolar .AND. coefs%coef%fmv_model_ver /= 9) THEN
    err = errorstatus_fatal
    msg = 'Coefficient file does not support solar calculations'
    THROWM(err .NE. 0, msg)
  ENDIF


  ! All IR instruments
  IF (coefs%coef%id_sensor == sensor_id_ir .OR. coefs%coef%id_sensor == sensor_id_hi) THEN
    IF (opts%rt_mw%clw_data) THEN
      err = errorstatus_fatal
      msg = 'Cloud liquid water only applicable to MW instruments'
      THROWM(err .NE. 0, msg)
    ENDIF
  ENDIF


  ! MW instruments
  IF (coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po) THEN
    IF (opts%rt_ir%addsolar) THEN
      err = errorstatus_fatal
      msg = 'Solar calculations not applicable to MW instruments'
      THROWM(err .NE. 0, msg)
    ENDIF
    IF (opts%rt_ir%addclouds) THEN
      err = errorstatus_fatal
      msg = 'addclouds not applicable to MW instruments: use rttov_scatt interface instead'
      THROWM(err .NE. 0, msg)
    ENDIF
    IF (opts%rt_ir%addaerosl) THEN
      err = errorstatus_fatal
      msg = 'Aerosol calculations not applicable to MW instruments'
      THROWM(err .NE. 0, msg)
    ENDIF
  ENDIF


  ! NLTE
  IF (opts%rt_ir%do_nlte_correction) THEN
    IF (.NOT. coefs%coef%nltecoef) THEN
      err = errorstatus_fatal
      msg = 'Coefficient file does not support NLTE correction or channel selection does not include any NLTE-affected channels'
      THROWM(err .NE. 0, msg)
    ENDIF
  ENDIF


  ! PC-RTTOV
  IF (opts%rt_ir%pc%addpc) THEN
    IF (coefs%coef%id_comp_pc == 0_jpim) THEN
      err = errorstatus_fatal
      msg = 'PC-RTTOV: optical depth coefficient file is not compatible with PCs'
      THROWM(err .NE. 0, msg)
    ENDIF

    IF (opts%rt_ir%addclouds .AND. coefs%coef_pccomp%fmv_pc_cld == 0_jpim) THEN
      err = errorstatus_fatal
      msg = 'PC-RTTOV: addclouds is TRUE, but PC coef file is not compatible with cloudy simulations'
      THROWM(err .NE. 0, msg)
    ENDIF

    IF (opts%rt_ir%pc%ipcbnd < 1_jpim .OR. &
      & opts%rt_ir%pc%ipcbnd > coefs%coef_pccomp%fmv_pc_bands) THEN
      err = errorstatus_fatal
      msg = 'PC-RTTOV: invalid spectral band (opts%rt_ir%pc%ipcbnd)'
      THROWM(err .NE. 0, msg)
    ENDIF

    IF (opts%rt_ir%pc%ipcreg < 1_jpim .OR. &
      & opts%rt_ir%pc%ipcreg > coefs%coef_pccomp%fmv_pc_sets(opts%rt_ir%pc%ipcbnd)) THEN
      err = errorstatus_fatal
      msg = 'PC-RTTOV: invalid predictor set (opts%rt_ir%pc%ipcreg) for selected spectral band (opts%rt_ir%pc%ipcbnd)'
      THROWM(err .NE. 0, msg)
    ENDIF

    IF (opts%rt_all%switchrad .AND. .NOT. opts%rt_ir%pc%addradrec) THEN
      err = errorstatus_fatal
      msg = 'PC-RTTOV: switchrad will have no effect unless addradrec is TRUE'
      THROWM(err .NE. 0, msg)
    ENDIF

    IF (opts%rt_ir%addaerosl) THEN
      err = errorstatus_fatal
      msg = 'PC-RTTOV: cannot be run with aerosol profiles'
      THROWM(err .NE. 0, msg)
    ENDIF

    IF (opts%rt_ir%addclouds .AND. opts%rt_ir%user_cld_opt_param) THEN
      err = errorstatus_fatal
      msg = 'PC-RTTOV: cannot be run with user input cloud optical parameters'
      THROWM(err .NE. 0, msg)
    ENDIF

    IF (opts%rt_ir%addsolar) THEN
      err = errorstatus_fatal
      msg = 'PC-RTTOV: cannot be run with addsolar TRUE'
      THROWM(err .NE. 0, msg)
    ENDIF

    IF (opts%rt_ir%do_nlte_correction) THEN
      err = errorstatus_fatal
      msg = 'PC-RTTOV: cannot be run with do_nlte_correction TRUE'
      THROWM(err .NE. 0, msg)
    ENDIF

    IF (opts%rt_all%do_lambertian .OR. opts%rt_ir%do_lambertian) THEN
      err = errorstatus_fatal
      msg = 'PC-RTTOV: cannot be run with do_lambertian TRUE'
      THROWM(err .NE. 0, msg)
    ENDIF
  ENDIF


  IF (LHOOK) CALL DR_HOOK('RTTOV_USER_OPTIONS_CHECKINPUT',1_jpim,ZHOOK_HANDLE)

  CATCH

  IF (LHOOK) CALL DR_HOOK('RTTOV_USER_OPTIONS_CHECKINPUT',1_jpim,ZHOOK_HANDLE)
END SUBROUTINE rttov_user_options_checkinput
