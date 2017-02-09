! Description:
!> @file
!!   Runs RTTOV tangent linear (TL) model
!
!> @brief
!!   Runs RTTOV tangent linear (TL) model
!!
!! @details
!!   Given an input profile and profile perturbation, computes
!!   the corresponding radiance perturbation from the tangent
!!   linear (TL) to the direct model evaluated at the input profile.
!!
!! @param[out]    errorstatus     status on exit
!! @param[in]     chanprof        specifies channels and profiles to simulate
!! @param[in]     opts            options to configure the simulations
!! @param[in]     profiles        input atmospheric profiles and surface variables
!! @param[in]     profiles_tl     input atmospheric profile and surface variable perturbations
!! @param[in]     coefs           coefficients structure for instrument to simulate
!! @param[in,out] transmission    output transmittances
!! @param[in,out] transmission_tl output transmittance perturbations
!! @param[in,out] radiancedata    output radiances and corresponding BTs and BRFs
!! @param[in,out] radiancedata_tl output radiance, BT and BRF perturbations
!! @param[in]     calcemis        flags for internal RTTOV surface emissivity calculation, optional
!! @param[in,out] emissivity      input/output surface emissivities, optional
!! @param[in,out] emissivity_tl   input/output surface emissivity perturbations, optional
!! @param[in]     calcrefl        flags for internal RTTOV surface BRDF calculation, optional
!! @param[in,out] reflectance     input/output surface BRDFs, input cloud top BRDF for simple cloud, optional
!! @param[in,out] reflectance_tl  input/output surface BRDF perturbations, optional
!! @param[in]     aer_opt_param   input aerosol optical parameters, optional
!! @param[in]     cld_opt_param   input cloud optical parameters, optional
!! @param[in,out] traj            RTTOV direct internal state, can be initialised outside RTTOV, optional
!! @param[in,out] traj_tl         RTTOV TL internal state, can be initialised outside RTTOV, optional
!! @param[in,out] pccomp          output PC scores and radiances from PC-RTTOV, optional
!! @param[in,out] pccomp_tl       output PC score and radiance perturbations from PC-RTTOV, optional
!! @param[in]     channels_rec    list of channels for which to calculate reconstructed radiances, optional
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
SUBROUTINE rttov_tl( &
            & errorstatus,       &
            & chanprof,          &
            & opts,              &
            & profiles,          &
            & profiles_tl,       &
            & coefs,             &
            & transmission,      &
            & transmission_tl,   &
            & radiancedata,      &
            & radiancedata_tl,   &
            & calcemis,          &
            & emissivity,        &
            & emissivity_tl,     &
            & calcrefl,          &
            & reflectance,       &
            & reflectance_tl,    &
            & aer_opt_param,     &
            & cld_opt_param,     &
            & traj,              &
            & traj_tl,           &
            & pccomp,            &
            & pccomp_tl,         &
            & channels_rec)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY :  &
       & rttov_options,     &
       & radiance_type,     &
       & transmission_type, &
       & profile_type,      &
       & rttov_coefs,       &
       & rttov_pccomp,      &
       & rttov_chanprof,    &
       & rttov_emissivity,  &
       & rttov_reflectance, &
       & rttov_opt_param,   &
       & rttov_traj
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE parkind1, ONLY : jprb
  USE rttov_const, ONLY :  &
         sensor_id_mw, &
         sensor_id_ir, &
         sensor_id_po, &
         sensor_id_hi, &
         ncldtyp
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE rttov_types, ONLY : &
       & rttov_traj_dyn,  &
       & rttov_traj_sta
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim)     , INTENT(OUT)                     :: errorstatus
  TYPE(profile_type)     , INTENT(IN)                      :: profiles(:)
  TYPE(rttov_chanprof)   , INTENT(IN)                      :: chanprof(:)
  TYPE(rttov_options )   , INTENT(IN)                      :: opts
  TYPE(profile_type  )   , INTENT(IN)                      :: profiles_tl(SIZE(profiles))
  TYPE(rttov_coefs   )   , INTENT(IN)   , TARGET           :: coefs
  TYPE(transmission_type), INTENT(INOUT)                   :: transmission
  TYPE(transmission_type), INTENT(INOUT)                   :: transmission_tl
  TYPE(radiance_type    ), INTENT(INOUT)                   :: radiancedata
  TYPE(radiance_type    ), INTENT(INOUT)                   :: radiancedata_tl
  LOGICAL(KIND=jplm)     , INTENT(IN)   , OPTIONAL         :: calcemis(SIZE(chanprof))
  TYPE(rttov_emissivity) , INTENT(INOUT), OPTIONAL         :: emissivity(SIZE(chanprof))
  TYPE(rttov_emissivity) , INTENT(INOUT), OPTIONAL         :: emissivity_tl(SIZE(chanprof))
  LOGICAL(KIND=jplm)     , INTENT(IN)   , OPTIONAL         :: calcrefl(SIZE(chanprof))
  TYPE(rttov_reflectance), INTENT(INOUT), OPTIONAL         :: reflectance(SIZE(chanprof))
  TYPE(rttov_reflectance), INTENT(INOUT), OPTIONAL         :: reflectance_tl(SIZE(chanprof))
  TYPE(rttov_opt_param)  , INTENT(IN)   , OPTIONAL         :: aer_opt_param
  TYPE(rttov_opt_param)  , INTENT(IN)   , OPTIONAL         :: cld_opt_param
  TYPE(rttov_traj       ), INTENT(INOUT), OPTIONAL, TARGET :: traj, traj_tl    ! target is needed here
  TYPE(rttov_pccomp     ), INTENT(INOUT), OPTIONAL         :: pccomp
  TYPE(rttov_pccomp     ), INTENT(INOUT), OPTIONAL         :: pccomp_tl
  INTEGER(KIND=jpim)     , INTENT(IN)   , OPTIONAL         :: channels_rec(:)
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_convert_gas_dry_tl.interface"
#include "rttov_apply_reg_limits_tl.interface"
#include "rttov_calcemis_ir_tl.interface"
#include "rttov_calcsurfrefl_tl.interface"
#include "rttov_profaux_tl.interface"
#include "rttov_setgeometry_tl.interface"
#include "rttov_setpredictors_7_tl.interface"
#include "rttov_setpredictors_8_tl.interface"
#include "rttov_setpredictors_9_tl.interface"
#include "rttov_opdep_tl.interface"
#include "rttov_opdep_9_tl.interface"
#include "rttov_transmit_tl.interface"
#include "rttov_transmit_9_solar_tl.interface"
#include "rttov_calcemis_mw_tl.interface"
#include "rttov_integrate_tl.interface"
#include "rttov_intavg_prof_tl.interface"
#include "rttov_intavg_chan_tl.interface"
#include "rttov_refsun_tl.interface"
#include "rttov_check_traj.interface"
#include "rttov_init_prof.interface"
#include "rttov_cldstr_tl.interface"
#include "rttov_opdpscattir_tl.interface"
#include "rttov_fresnel_tl.interface"
#include "rttov_init_raytracing.interface"
#include "rttov_pcscores_tl.interface"
#include "rttov_reconstruct_tl.interface"
#include "rttov_calcbt_pc_tl.interface"
#include "rttov_copy_raytracing.interface"
#include "rttov_copy_prof.interface"
#include "rttov_copy_aux_prof.interface"
#include "rttov_copy_opdp_path.interface"
#include "rttov_init_opdp_path.interface"
#include "rttov_alloc_traj_dyn.interface"
#include "rttov_alloc_traj_sta.interface"
#include "rttov_direct.interface"
#include "rttov_nlte_bias_correction_tl.interface"
#include "rttov_calcbt_tl.interface"
#include "rttov_calcsatrefl_tl.interface"

  INTEGER(KIND=jpim) :: nlevels
  LOGICAL(KIND=jplm) :: addcosmic! switch for adding temp of cosmic background

  TYPE(rttov_traj_dyn) :: traj0_dyn
  TYPE(rttov_traj_sta) :: traj0_sta
  TYPE(rttov_traj_dyn) :: traj0_tl_dyn

  LOGICAL(KIND=jplm)   :: ltraj_tl_dyn_dealloc

  TYPE(rttov_traj), TARGET  :: traj1
  TYPE(rttov_traj), POINTER :: traj0
  TYPE(rttov_traj), TARGET  :: traj1_tl
  TYPE(rttov_traj), POINTER :: traj0_tl
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim)  :: npcscores
  TYPE(rttov_options) :: opts_coef
  INTEGER(KIND=jpim)  :: err
  REAL   (KIND=JPRB)  :: ZHOOK_HANDLE

!- End of header --------------------------------------------------------
  TRY
!-------------
!0. initialize
!-------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_TL', 0_jpim, ZHOOK_HANDLE)
  nprofiles                 = SIZE(profiles)
  nchannels                 = SIZE(chanprof)
  nlevels                   = profiles(1)%nlevels
  opts_coef                 = opts
  opts_coef%rt_ir%addclouds = .FALSE.
  opts_coef%rt_ir%addaerosl = .FALSE.
  errorstatus               = errorstatus_success

  ltraj_tl_dyn_dealloc = .FALSE.

  IF (opts%rt_ir%pc%addpc) THEN
    npcscores = SIZE(pccomp%pcscores)
    IF (npcscores / nprofiles > coefs%coef_pccomp%fmv_pc_mnum) THEN
      err = errorstatus_fatal
      THROWM(err.ne.0, "npcscores is larger than maximum allowed by PC coefficient file")
    ENDIF
  ENDIF

  NULLIFY (traj0)
  NULLIFY (traj0_tl)
!-----------------------------------------------------------------------
!1. interpolator first call - input profiles from USER levs to COEF levs
!-----------------------------------------------------------------------
! done immediately in direct code - for quick escape if target values are out of bounds
  CALL rttov_check_traj( &
        & err,                 &
        & nprofiles,           &
        & nchannels,           &
        & opts,                &
        & nlevels,             &
        & coefs,               &
        & 1_jpim,              &
        & traj0 = traj0,       &
        & traj0_tl = traj0_tl, &
        & traj1 = traj1,       &
        & traj1_tl = traj1_tl, &
        & traj2 = traj,        &
        & traj2_tl = traj_tl)
  THROWM( err .NE. 0 , "rttov_check_traj fatal error")

  CALL rttov_direct( &
            & err,                           &
            & chanprof,                      &
            & opts,                          &
            & profiles,                      &
            & coefs,                         &
            & transmission,                  &
            & radiancedata,                  &
            & calcemis      = calcemis,      &
            & emissivity    = emissivity,    &
            & calcrefl      = calcrefl,      &
            & reflectance   = reflectance,   &
            & aer_opt_param = aer_opt_param, &
            & cld_opt_param = cld_opt_param, &
            & traj          = traj0,         &
            & traj_dyn      = traj0_dyn,     &
            & traj_sta      = traj0_sta,     &
            & pccomp        = pccomp,        &
            & channels_rec  = channels_rec)
  THROW(err.ne.0)

  IF (traj0_sta%dothermal .AND. .NOT. PRESENT(emissivity_tl)) THEN
    err = errorstatus_fatal
    THROWM(err.ne.0, "emissivity_tl parameter required")
  END IF
  IF (traj0_sta%dosolar .AND. .NOT. PRESENT(reflectance_tl)) THEN
    err = errorstatus_fatal
    THROWM(err.ne.0, "reflectance_tl parameter required")
  END IF

  CALL rttov_alloc_traj_dyn (err, traj0_tl_dyn, opts, nchannels, profiles(1)%nlayers, &
                             traj0_dyn%nstreams, ncldtyp, 1_jpim)
  THROW(err.ne.0)
  ltraj_tl_dyn_dealloc = .TRUE.
!
! Tangent Linear
!----------------

!-----------------------------------------------------------------------
! Convert input profiles to ppmv wrt dry air
!-----------------------------------------------------------------------

  CALL rttov_convert_gas_dry_tl( &
      opts,                      &
      profiles,                  &
      profiles_tl,               &
      traj0%profiles_dry,        &
      traj0_tl%profiles_dry)

!---------------------------------------------------------------------------------------
!1. interpolator first call - input profiles from USER levs to COEF levs
!---------------------------------------------------------------------------------------
  CALL rttov_init_prof(traj0_tl%profiles_coef)

  IF (opts%interpolation%addinterp) THEN
    CALL rttov_intavg_prof_tl( &
          & opts,                           &
          & nlevels,                        &
          & coefs%coef%nlevels,             &
          & profiles,                       &
          & profiles_tl,                    &
          & traj0%profiles_dry,             &
          & traj0_tl%profiles_dry,          &
          & traj0%profiles_coef,            &
          & traj0_tl%profiles_coef,         &
          & coefs%coef,                     &
          & coefs%coef_pccomp)
  ELSE
    CALL rttov_copy_prof( &
          & traj0_tl%profiles_coef, &
          & profiles_tl,            &
          & larray = .TRUE._jplm,   &
          & lscalar = .FALSE._jplm, &
          & profiles_gas = traj0_tl%profiles_dry)
  ENDIF

! complete profiles on COEF levs for surf, skin, cloud, aerosol, solar, angle
  CALL rttov_copy_prof( &
        & traj0_tl%profiles_coef, &
        & profiles_tl,            &
        & larray = .FALSE._jplm,  &
        & lscalar = .TRUE._jplm,  &
        & profiles_gas = traj0_tl%profiles_dry)
!------------------------------------------------------------------
!2. TL of apply_reg_limits - COEF levs
!------------------------------------------------------------------

  IF (opts%config%apply_reg_limits .OR. opts%interpolation%reg_limit_extrap) THEN
    CALL rttov_apply_reg_limits_tl( &
          & opts,                        &
          & profiles,                    &
          & traj0_sta%profiles_coef_ref, &
          & traj0_tl%profiles_coef,      &
          & coefs%coef,                  &
          & coefs%coef_pccomp)
  ENDIF

!------------------------------------------------------------------------
!3. determine cloud top, surface levels, ice cloud parameters
!-----------------------------------------------------------------------
!3.2 COEF levs
!-------------
!   CALL rttov_init_aux_prof(traj0_tl%aux_prof_coef)
  CALL rttov_profaux_tl( &
    & opts_coef,              &
    & traj0%profiles_coef,    &
    & traj0_tl%profiles_coef, &
    & coefs%coef,             &
    & traj0%aux_prof_coef,    &
    & traj0_tl%aux_prof_coef)

!3.1 USER levs
!-------------
  IF (opts%interpolation%addinterp .or. opts%rt_ir%addclouds) THEN
   CALL rttov_profaux_tl( &
     & opts,              &
     & profiles,          &
     & profiles_tl,       &
     & coefs%coef,        &
     & traj0%aux_prof,    &
     & traj0_tl%aux_prof)
 ELSE
   CALL rttov_copy_aux_prof(traj0_tl%aux_prof, traj0_tl%aux_prof_coef)
 ENDIF

! TL on geometry
!------------------------------------------------------------------------
!4. set up common geometric variables for the RT integration
!------------------------------------------------------------------------
!4.2 COEF levs
!-------------
 CALL rttov_setgeometry_tl( &
   & opts,                     &
   & traj0_sta%dosolar,        &
   & traj0%profiles_coef,      &
   & traj0_tl%profiles_coef,   &
   & traj0%aux_prof_coef,      &
   & coefs%coef,               &
   & traj0_sta%angles_coef,    &
   & traj0%raytracing_coef,    &
   & traj0_tl%raytracing_coef, &
   & traj0%profiles_coef,      &
   & traj0_tl%profiles_coef,   &
   & .TRUE._jplm)  ! do_pmc_calc on coef levels

!4.1 USER levs
!-------------
  CALL rttov_init_raytracing(traj0_sta%dosolar, traj0_tl%raytracing)
  IF (opts%interpolation%addinterp) THEN
    CALL rttov_setgeometry_tl(   &
        & opts,                  &
        & traj0_sta%dosolar,     &
        & profiles,              &
        & profiles_tl,           &
        & traj0%aux_prof,        &
        & coefs%coef,            &
        & traj0_sta%angles,      &
        & traj0%raytracing,      &
        & traj0_tl%raytracing,   &
        & traj0%profiles_dry,    &
        & traj0_tl%profiles_dry)
  ELSE
    traj0_sta%angles = traj0_sta%angles_coef
    ! copy raytracing_coef *in to* raytracing_user
    CALL rttov_copy_raytracing(traj0_sta%dosolar, traj0_tl%raytracing, traj0_tl%raytracing_coef)
  ENDIF

! TL of predictors
!---------------------------------------------------
!5. calculate atm/sol predictors - atm/sol COEF levs
!---------------------------------------------------

  IF (coefs%coef%fmv_model_ver == 7) THEN
    CALL rttov_setpredictors_7_tl( &
          & opts,                      &
          & traj0%profiles_coef,       &
          & traj0_tl%profiles_coef,    &
          & traj0_sta%angles_coef,     &
          & coefs%coef,                &
          & traj0%aux_prof_coef,       &
          & traj0_tl%aux_prof_coef,    &
          & traj0%predictors%path1,    &
          & traj0_tl%predictors%path1, &
          & traj0%raytracing_coef,     &
          & traj0_tl%raytracing_coef)
  ELSE IF (coefs%coef%fmv_model_ver == 8) THEN
    CALL rttov_setpredictors_8_tl( &
          & opts,                      &
          & traj0%profiles_coef,       &
          & traj0_tl%profiles_coef,    &
          & traj0_sta%angles_coef,     &
          & coefs%coef,                &
          & traj0%aux_prof_coef,       &
          & traj0%predictors%path1,    &
          & traj0_tl%predictors%path1, &
          & traj0%raytracing_coef,     &
          & traj0_tl%raytracing_coef)
  ELSE IF (coefs%coef%fmv_model_ver == 9) THEN
    IF (traj0_sta%dothermal) THEN
      CALL rttov_setpredictors_9_tl(            &
            & opts,                             &
            & traj0%profiles_coef,              &
            & traj0_tl%profiles_coef,           &
            & traj0_sta%angles_coef,            &
            & traj0%raytracing_coef%pathsat,    &
            & traj0_tl%raytracing_coef%pathsat, &
            & coefs%coef_pccomp,                &
            & coefs%coef,                       &
            & traj0%predictors%path1,           &
            & traj0_tl%predictors%path1)
    ENDIF
    IF (traj0_sta%dosolar) THEN
      CALL rttov_setpredictors_9_tl(            &
            & opts,                             &
            & traj0%profiles_coef,              &
            & traj0_tl%profiles_coef,           &
            & traj0_sta%angles_coef,            &
            & traj0%raytracing_coef%patheff,    &
            & traj0_tl%raytracing_coef%patheff, &
            & coefs%coef_pccomp,                &
            & coefs%coef,                       &
            & traj0%predictors%path2,           &
            & traj0_tl%predictors%path2)
    ENDIF
  ENDIF

!TL of optical depths
!--------------------------------------------------------------------
!6. predict atmospheric and solar optical depths - atm/sol COEF levs
!--------------------------------------------------------------------
  CALL rttov_init_opdp_path(opts, traj0_tl%opdp_path)
  CALL rttov_init_opdp_path(opts, traj0_tl%opdp_path_coef)

  IF (coefs%coef%fmv_model_ver == 9) THEN
    IF (traj0_sta%dothermal) THEN
      ! Calculate thermal path1 optical depths for all thermal channels
      CALL rttov_opdep_9_tl( &
            & coefs%coef%nlayers,                  &
            & chanprof,                            &
            & traj0_sta%thermal,                   &
            & traj0%predictors,                    &
            & traj0%predictors%path1,              &
            & traj0_tl%predictors%path1,           &
            & traj0%aux_prof_coef,                 &
            & traj0_tl%aux_prof_coef,              &
            & coefs%coef,                          &
            & coefs%coef%thermal,                  &
            & traj0%opdp_path_coef%atm_level,      &
            & traj0_tl%opdp_path_coef%atm_level,   &
            & traj0_sta%thermal_path1%opdp_ref_coef)
    ENDIF
    IF (traj0_sta%dosolar) THEN
      ! Calculate solar path2 optical depths for all solar channels
      CALL rttov_opdep_9_tl( &
            & coefs%coef%nlayers,                      &
            & chanprof,                                &
            & traj0_sta%solar,                         &
            & traj0%predictors,                        &
            & traj0%predictors%path2,                  &
            & traj0_tl%predictors%path2,               &
            & traj0%aux_prof_coef,                     &
            & traj0_tl%aux_prof_coef,                  &
            & coefs%coef,                              &
            & coefs%coef%solar,                        &
            & traj0%opdp_path_coef%sun_level_path2,    &
            & traj0_tl%opdp_path_coef%sun_level_path2, &
            & traj0_sta%solar_path2%opdp_ref_coef)
    ENDIF
  ELSE
    CALL rttov_opdep_tl( &
          & coefs%coef%nlayers,        &
          & chanprof,                  &
          & traj0%predictors,          &
          & traj0%predictors%path1,    &
          & traj0_tl%predictors%path1, &
          & traj0%aux_prof_coef,       &
          & traj0_tl%aux_prof_coef,    &
          & coefs%coef,                &
          & coefs%coef%thermal,        &
          & traj0%opdp_path_coef,      &
          & traj0_tl%opdp_path_coef,   &
          & traj0_sta%thermal_path1%opdp_ref_coef)
  ENDIF

!--------------------------------------------------------------------------
!7. interpolator second  call - optical depths from COEF levs to USER levs
!--------------------------------------------------------------------------

  IF (opts%interpolation%addinterp) THEN
    CALL rttov_intavg_chan_tl( &
          & opts,                           &
          & traj0_sta%thermal,              &
          & traj0_sta%solar,                &
          & coefs%coef%nlevels,             &
          & nlevels,                        &
          & chanprof,                       &
          & traj0%profiles_coef,            &
          & profiles,                       &
          & profiles_tl,                    &
          & traj0%opdp_path_coef,           &
          & traj0_tl%opdp_path_coef,        &
          & traj0%opdp_path,                &
          & traj0_tl%opdp_path)
  ELSE
    CALL rttov_copy_opdp_path(opts, traj0_tl%opdp_path, traj0_tl%opdp_path_coef)
  ENDIF

  IF (opts%interpolation%spacetop) THEN
    traj0_tl%opdp_path%atm_level(1,:) = 0._jprb
    IF (opts%rt_ir%addsolar) THEN
      traj0_tl%opdp_path%sun_level_path2(1,:) = 0._jprb
    ENDIF
  ENDIF

! TL of ircld
!--------------------------------------------------------------------
!8. If clouds are present,calculate the number of streams and
!   the cloud distribution  in each stream - remains on USER levs
!--------------------------------------------------------------------

  IF (opts%rt_ir%addclouds) THEN
    CALL rttov_cldstr_tl( &
          & opts%rt_ir,     &
          & profiles,       &
          & profiles_tl,    &
          & traj0%ircld,    &
          & traj0_tl%ircld)
  ELSE
    traj0_tl%ircld%xstrclr = 0._jprb
  ENDIF

! TL of opdpscattir
!----------------------------------------------------------------------------
!9. Calculate optical depths of aerosols and/or clouds - USER levs
!----------------------------------------------------------------------------

  IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
    CALL rttov_opdpscattir_tl( &
          & profiles(1)%nlayers,                           &
          & chanprof,                                      &
          & opts,                                          &
          & traj0%aux_prof,                                &
          & traj0_tl%aux_prof,                             &
          & profiles,                                      &
          & profiles_tl,                                   &
          & traj0_tl%profiles_dry,                         &
          & aer_opt_param,                                 &
          & cld_opt_param,                                 &
          & traj0_sta%dosolar,                             &
          & traj0_sta%solar,                               &
          & coefs%coef,                                    &
          & coefs%coef_scatt_ir,                           &
          & traj0%raytracing,                              &
          & traj0_tl%raytracing,                           &
          & traj0%transmission_scatt_ir,                   &
          & traj0_tl%transmission_scatt_ir,                &
          & traj0_dyn%transmission_scatt_ir_stream,        &
          & traj0_tl_dyn%transmission_scatt_ir_stream,     &
          & coefs%optp,                                    &
          & traj0%ircld,                                   &
          & traj0_tl%ircld)
  ENDIF

! TL of transmit
!----------------------------------------
!10. calculate transmittances - USER levs
!----------------------------------------

  IF (traj0_sta%dothermal) THEN
    CALL rttov_transmit_tl( &
          & opts%rt_ir%addaerosl,                           &
          & opts%rt_ir%addclouds,                           &
          & traj0_sta%do_lambertian,                        &
          & profiles(1)%nlayers,                            &
          & chanprof,                                       &
          & traj0_sta%thermal,                              &
          & traj0%aux_prof,                                 &
          & traj0_tl%aux_prof,                              &
          & coefs%coef,                                     &
          & traj0%ircld,                                    &
          & traj0_sta%angles,                               &
          & traj0%opdp_path%atm_level,                      &
          & traj0_tl%opdp_path%atm_level,                   &
          & traj0_sta%thermal_path1%od_level,               &
          & transmission%tau_levels,                        &
          & transmission%tau_total,                         &
          & transmission_tl%tau_levels,                     &
          & transmission_tl%tau_total,                      &
          & traj0_dyn%transmission_aux,                     &
          & traj0_dyn%transmission_aux%thermal_path1,       &
          & traj0_tl_dyn%transmission_aux%thermal_path1,    &
          & traj0%transmission_scatt_ir,                    &
          & traj0_tl%transmission_scatt_ir,                 &
          & traj0_dyn%transmission_scatt_ir_stream,         &
          & traj0_tl_dyn%transmission_scatt_ir_stream,      &
          & traj0_sta%thermal_path1%tau_ref,                &
          & traj0_sta%thermal_path1%tau_ref_surf,           &
          & traj0_sta%thermal_path1%tau_surf,               &
          & traj0_sta%thermal_path1%tau_level)
  ENDIF

  IF (traj0_sta%dosolar) THEN
    CALL rttov_transmit_9_solar_tl( &
          & opts%rt_ir%addaerosl,                          &
          & opts%rt_ir%addclouds,                          &
          & profiles(1)%nlayers,                           &
          & nprofiles,                                     &
          & chanprof,                                      &
          & traj0_sta%solar,                               &
          & traj0%aux_prof,                                &
          & traj0_tl%aux_prof,                             &
          & coefs%coef,                                    &
          & traj0%raytracing,                              &
          & traj0_tl%raytracing,                           &
          & traj0%ircld,                                   &
          & traj0%opdp_path,                               &
          & traj0_tl%opdp_path,                            &
          & traj0_sta%solar_path2,                         &
          & traj0_sta%solar_path1,                         &
          & transmission,                                  &
          & transmission_tl,                               &
          & traj0_dyn%transmission_aux,                    &
          & traj0_tl_dyn%transmission_aux,                 &
          & traj0%transmission_scatt_ir,                   &
          & traj0_tl%transmission_scatt_ir,                &
          & traj0_dyn%transmission_scatt_ir_stream,        &
          & traj0_tl_dyn%transmission_scatt_ir_stream)

  ENDIF

! TL of emissivity
!-------------------------------------------------
!11. calculate channel emissivity values - SURFACE
!-------------------------------------------------

  IF (traj0_sta%dothermal) THEN

    emissivity_tl%emis_out = emissivity_tl%emis_in

    traj0_tl%thermrefl =  - emissivity_tl%emis_out

    IF (ANY(calcemis)) THEN
      ! Calculate surface emissivity and traj0%thermrefl for selected channels

      IF (coefs%coef%id_sensor == sensor_id_ir) THEN
        ! Infrared
        traj0_tl%thermrefl(:) =  - emissivity_tl(:)%emis_out
      ELSE IF (coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po) THEN
        ! Microwave
        CALL rttov_calcemis_mw_tl( &
              & opts,                             &
              & profiles,                         &
              & profiles_tl,                      &
              & traj0_sta%angles,                 &
              & coefs%coef,                       &
              & chanprof,                         &
              & traj0_dyn%transmission_aux,       &
              & traj0_tl_dyn%transmission_aux,    &
              & calcemis,                         &
              & emissivity_tl%emis_out,           &
              & traj0_tl%thermrefl)
      ELSE
        ! Hires
        CALL rttov_calcemis_ir_tl(    &
              & profiles,             &
              & profiles_tl,          &
              & coefs%coef,           &
              & opts%rt_ir%pc%addpc,  &
              & coefs%coef_pccomp,    &
              & traj0_sta%thermal,    &
              & chanprof,             &
              & calcemis,             &
              & emissivity_tl%emis_out)
        traj0_tl%thermrefl(:) =  - emissivity_tl(:)%emis_out
      ENDIF

    ENDIF

  ELSE IF (PRESENT(emissivity_tl)) THEN

    emissivity_tl%emis_out = 0._jprb

  ENDIF

! TL of reflectance
!-------------------------------------------------------
!12. Compute surface reflectance
!-------------------------------------------------------

  IF (traj0_sta%dosolar) THEN

    reflectance_tl%refl_out = reflectance_tl%refl_in

    IF (ANY(calcrefl)) THEN

      CALL rttov_refsun_tl( &
            & profiles,            &
            & profiles_tl,         &
            & coefs%coef,          &
            & traj0%aux_prof,      &
            & traj0%sunglint,      &
            & traj0_tl%sunglint,   &
            & traj0%raytracing,    &
            & traj0_tl%raytracing)
      CALL rttov_fresnel_tl( &
            & chanprof,           &
            & profiles,           &
            & traj0_sta%solar,    &
            & coefs%coef,         &
            & traj0%sunglint,     &
            & traj0_tl%sunglint,  &
            & traj0%fresnrefl,    &
            & traj0_tl%fresnrefl)

      CALL rttov_calcsurfrefl_tl(     &
            & coefs%coef,             &
            & profiles,               &
            & traj0%sunglint,         &
            & traj0_tl%sunglint,      &
            & traj0%fresnrefl,        &
            & traj0_tl%fresnrefl,     &
            & traj0_sta%solar,        &
            & chanprof,               &
            & traj0_sta%refl_norm,    &
            & calcrefl,               &
            & emissivity,             &
            & emissivity_tl,          &
            & reflectance_tl%refl_out)

    ENDIF

  ELSE IF (PRESENT(reflectance_tl)) THEN

    reflectance_tl%refl_out = 0._jprb

  ENDIF

! TL of RTE
!---------------------------------------------------------
!13. integrate the radiative transfer equation - USER levs
!---------------------------------------------------------

  addcosmic = (coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po)
  CALL rttov_integrate_tl( &
        & addcosmic,                                    &
        & opts,                                         &
        & traj0_dyn%nstreams,                           &
        & chanprof,                                     &
        & emissivity,                                   &
        & emissivity_tl,                                &
        & reflectance,                                  &
        & reflectance_tl,                               &
        & traj0_sta%refl_norm,                          &
        & traj0%thermrefl,                              &
        & traj0_tl%thermrefl,                           &
        & traj0_sta%do_lambertian,                      &
        & traj0_sta%thermal,                            &
        & traj0_sta%dothermal,                          &
        & traj0_sta%solar,                              &
        & traj0_sta%dosolar,                            &
        & traj0_sta%solar_spec_esd,                     &
        & traj0_dyn%transmission_aux,                   &
        & traj0_tl_dyn%transmission_aux,                &
        & traj0%transmission_scatt_ir,                  &
        & traj0_tl%transmission_scatt_ir,               &
        & profiles,                                     &
        & profiles_tl,                                  &
        & traj0%profiles_dry,                           &
        & traj0_tl%profiles_dry,                        &
        & traj0%aux_prof,                               &
        & traj0_tl%aux_prof,                            &
        & coefs%coef,                                   &
        & traj0%raytracing,                             &
        & traj0_tl%raytracing,                          &
        & traj0%ircld,                                  &
        & traj0_tl%ircld,                               &
        & radiancedata,                                 &
        & traj0_sta%auxrad,                             &
        & traj0_dyn%auxrad_stream,                      &
        & traj0_tl_dyn%auxrad_stream,                   &
        & radiancedata_tl)


!---------------------------------------------------------
! 14. Calculate output BTs/reflectances/PCscores
!---------------------------------------------------------

  IF (opts%rt_ir%pc%addpc) THEN
    CALL rttov_pcscores_tl(          &
          & opts,                    &
          & chanprof,                &
          & traj0_sta%chanprof_pc,   &
          & pccomp,                  &
          & pccomp_tl,               &
          & coefs%coef_pccomp,       &
          & radiancedata_tl)

    IF (opts%rt_ir%pc%addradrec) THEN
      CALL rttov_reconstruct_tl(     &
            & opts,                  &
            & traj0_sta%chanprof_in, &
            & traj0_sta%chanprof_pc, &
            & pccomp,                &
            & pccomp_tl,             &
            & coefs%coef_pccomp)
    ENDIF
  ENDIF

! Do NLTE bias correction for hyperspectral instruments (not PC-RTTOV)
  IF (.NOT. opts%rt_ir%pc%addpc) THEN
    IF (coefs%coef%nltecoef) THEN
      IF (opts%rt_ir%do_nlte_correction .AND. &
          coefs%coef%id_sensor == sensor_id_hi) THEN
        CALL rttov_nlte_bias_correction_tl(coefs%coef, profiles, profiles_tl, &
                                           traj0_sta%angles, chanprof, radiancedata_tl)
      ENDIF
    ENDIF
  ENDIF

  IF (opts%rt_ir%pc%addpc) THEN

    IF (opts%rt_ir%pc%addradrec) &
      CALL rttov_calcbt_pc_tl(&
           traj0_sta%chanprof_in, &
           coefs%coef_pccomp,     &
           pccomp,                &
           pccomp_tl)

  ELSE

    IF (traj0_sta%dothermal) THEN
      CALL rttov_calcbt_tl(chanprof, coefs%coef, traj0_sta%thermal, radiancedata, radiancedata_tl)
    ELSE
      radiancedata_tl%bt_clear(:) = 0._jprb
      radiancedata_tl%bt(:) = 0._jprb
    ENDIF
    IF (traj0_sta%dosolar) THEN
      CALL rttov_calcsatrefl_tl(chanprof, profiles, traj0_sta%solar_spec_esd, traj0_sta%solar, radiancedata_tl)
    ELSE
      radiancedata_tl%refl_clear(:) = 0._jprb
      radiancedata_tl%refl(:) = 0._jprb
    ENDIF
  ENDIF

!--------------------
!15. deallocate memory
!--------------------

  CALL cleanup()

  IF (LHOOK) CALL DR_HOOK('RTTOV_TL', 1_jpim, ZHOOK_HANDLE)
  CATCH_C
  errorstatus = err
  IF (LHOOK) CALL DR_HOOK('RTTOV_TL', 1_jpim, ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE cleanup()
    INTEGER(KIND=jpim) :: error

    IF (traj0_dyn%nstreams >= 0) THEN
      CALL rttov_alloc_traj_dyn (error, traj0_dyn, opts, nchannels, profiles(1)%nlayers, &
                                 traj0_dyn%nstreams, ncldtyp, 0_jpim)
    ENDIF

    IF (ltraj_tl_dyn_dealloc) THEN
      CALL rttov_alloc_traj_dyn (error, traj0_tl_dyn, opts, nchannels, profiles(1)%nlayers, &
                                 traj0_dyn%nstreams, ncldtyp, 0_jpim)
    ENDIF

    IF (ASSOCIATED(traj0_sta%thermal)) THEN
      CALL rttov_alloc_traj_sta (error, traj0_sta, opts, coefs%coef, nlevels, nchannels, nprofiles, &
                                 0_jpim, npcscores, channels_rec)
    ENDIF

    IF (ASSOCIATED(traj0)) THEN
      CALL rttov_check_traj( &
            & error,               &
            & nprofiles,           &
            & nchannels,           &
            & opts,                &
            & nlevels,             &
            & coefs,               &
            & 0_jpim,              &
            & traj0 = traj0,       &
            & traj0_tl = traj0_tl, &
            & traj1 = traj1,       &
            & traj1_tl = traj1_tl, &
            & traj2 = traj,        &
            & traj2_tl = traj_tl)
    ENDIF
  END SUBROUTINE cleanup

END SUBROUTINE rttov_tl
