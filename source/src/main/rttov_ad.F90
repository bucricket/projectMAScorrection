! Description:
!> @file
!!   Runs RTTOV adjoint (AD) model
!
!> @brief
!!   Runs RTTOV adjoint (AD) model
!!
!! @details
!!   Given the gradient of a scalar function with respect to
!!   channel radiances, computes the gradient of the same scalar
!!   function with respect to the profile variables.
!!   This is the transpose of the tangent linear model.
!!
!!   All AD variables should be initialised to zero before calling.
!!   The input gradients should then be written to the radiancedata_ad
!!   structure. If opts\%rt_all\%switchrad is false the gradient should
!!   be wrt to radiance in radiancedata_ad%total, otherwise it is wrt
!!   to BT in radiancedata_ad\%bt. For visible/near-IR channels (at
!!   wavelengths less than 3 microns) the input gradients are always
!!   in radiancedata_ad\%total.
!!
!!   For PC-RTTOV the input gradients are specified in pccomp_ad instead.
!!   If opts\%rt_ir\%addpc\%addradrec is false the gradient should be wrt
!!   PC score in pccomp_ad\%pcscores. Otherwise it should be in
!!   pccomp_ad\%total_pccomp or pccomp_ad\%bt_pccomp according to the
!!   setting of opts\%rt_all\%switchrad (false/true respectively).
!!
!! @param[out]    errorstatus     status on exit
!! @param[in]     chanprof        specifies channels and profiles to simulate
!! @param[in]     opts            options to configure the simulations
!! @param[in]     profiles        input atmospheric profiles and surface variables
!! @param[in,out] profiles_ad     output gradient wrt atmospheric profile and surface variables
!! @param[in]     coefs           coefficients structure for instrument to simulate
!! @param[in,out] transmission    output transmittances
!! @param[in,out] transmission_ad input gradient wrt transmittances (usually zero)
!! @param[in,out] radiancedata    output radiances and corresponding BTs and BRFs
!! @param[in,out] radiancedata_ad input gradient wrt radiances or BTs
!! @param[in]     calcemis        flags for internal RTTOV surface emissivity calculation, optional
!! @param[in,out] emissivity      input/output surface emissivities, optional
!! @param[in,out] emissivity_ad   output gradient wrt surface emissivities, optional
!! @param[in]     calcrefl        flags for internal RTTOV surface BRDF calculation, optional
!! @param[in,out] reflectance     input/output surface BRDFs, input cloud top BRDF for simple cloud, optional
!! @param[in,out] reflectance_ad  output gradient wrt surface BRDFs, optional
!! @param[in]     aer_opt_param   input aerosol optical parameters, optional
!! @param[in]     cld_opt_param   input cloud optical parameters, optional
!! @param[in,out] traj            RTTOV direct internal state, can be initialised outside RTTOV, optional
!! @param[in,out] traj_ad         RTTOV AD internal state, can be initialised outside RTTOV, optional
!! @param[in,out] pccomp          output PC scores and radiances from PC-RTTOV, optional
!! @param[in,out] pccomp_ad       input gradient wrt PC scores, radiances or BTs for PC-RTTOV, optional
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
SUBROUTINE rttov_ad( &
            & errorstatus,       &
            & chanprof,          &
            & opts,              &
            & profiles,          &
            & profiles_ad,       &
            & coefs,             &
            & transmission,      &
            & transmission_ad,   &
            & radiancedata,      &
            & radiancedata_ad,   &
            & calcemis,          &
            & emissivity,        &
            & emissivity_ad,     &
            & calcrefl,          &
            & reflectance,       &
            & reflectance_ad,    &
            & aer_opt_param,     &
            & cld_opt_param,     &
            & traj,              &
            & traj_ad,           &
            & pccomp,            &
            & pccomp_ad,         &
            & channels_rec)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY :  &
       & rttov_coefs,       &
       & rttov_pccomp,      &
       & profile_type,      &
       & transmission_type, &
       & radiance_type,     &
       & rttov_options,     &
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
  USE rttov_types, ONLY :  &
       & rttov_traj_sta,   &
       & rttov_traj_dyn
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim)     , INTENT(OUT)                     :: errorstatus
  TYPE(profile_type  )   , INTENT(IN)                      :: profiles(:)
  TYPE(rttov_chanprof)   , INTENT(IN)                      :: chanprof(:)
  TYPE(rttov_options )   , INTENT(IN)                      :: opts
  TYPE(profile_type  )   , INTENT(INOUT)                   :: profiles_ad(SIZE(profiles))
  TYPE(rttov_coefs   )   , INTENT(IN)   , TARGET           :: coefs
  TYPE(transmission_type), INTENT(INOUT)                   :: transmission
  TYPE(transmission_type), INTENT(INOUT)                   :: transmission_ad
  TYPE(radiance_type    ), INTENT(INOUT)                   :: radiancedata
  TYPE(radiance_type    ), INTENT(INOUT)                   :: radiancedata_ad
  LOGICAL(KIND=jplm)     , INTENT(IN)   , OPTIONAL         :: calcemis(SIZE(chanprof))
  TYPE(rttov_emissivity) , INTENT(INOUT), OPTIONAL         :: emissivity(SIZE(chanprof))
  TYPE(rttov_emissivity) , INTENT(INOUT), OPTIONAL         :: emissivity_ad(SIZE(chanprof))
  LOGICAL(KIND=jplm)     , INTENT(IN)   , OPTIONAL         :: calcrefl(SIZE(chanprof))
  TYPE(rttov_reflectance), INTENT(INOUT), OPTIONAL         :: reflectance(SIZE(chanprof))
  TYPE(rttov_reflectance), INTENT(INOUT), OPTIONAL         :: reflectance_ad(SIZE(chanprof))
  TYPE(rttov_opt_param)  , INTENT(IN)   , OPTIONAL         :: aer_opt_param
  TYPE(rttov_opt_param)  , INTENT(IN)   , OPTIONAL         :: cld_opt_param
  TYPE(rttov_traj       ), INTENT(INOUT), OPTIONAL, TARGET :: traj, traj_ad       ! target is needed here
  TYPE(rttov_pccomp     ), INTENT(INOUT), OPTIONAL         :: pccomp
  TYPE(rttov_pccomp     ), INTENT(INOUT), OPTIONAL         :: pccomp_ad
  INTEGER(KIND=jpim)     , INTENT(IN)   , OPTIONAL         :: channels_rec(:)
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_convert_gas_dry_ad.interface"
#include "rttov_apply_reg_limits_ad.interface"
#include "rttov_setgeometry_ad.interface"
#include "rttov_calcemis_ir_ad.interface"
#include "rttov_calcsurfrefl_ad.interface"
#include "rttov_profaux_ad.interface"
#include "rttov_setpredictors_7_ad.interface"
#include "rttov_setpredictors_8_ad.interface"
#include "rttov_setpredictors_9_ad.interface"
#include "rttov_opdep_ad.interface"
#include "rttov_opdep_9_ad.interface"
#include "rttov_transmit_ad.interface"
#include "rttov_transmit_9_solar_ad.interface"
#include "rttov_calcemis_mw_ad.interface"
#include "rttov_integrate_ad.interface"
#include "rttov_intavg_prof_ad.interface"
#include "rttov_intavg_chan_ad.interface"
#include "rttov_refsun_ad.interface"
#include "rttov_init_prof.interface"
#include "rttov_check_traj.interface"
#include "rttov_init_predictor.interface"
#include "rttov_init_raytracing.interface"
#include "rttov_cldstr_ad.interface"
#include "rttov_opdpscattir_ad.interface"
#include "rttov_fresnel_ad.interface"
#include "rttov_pcscores_ad.interface"
#include "rttov_reconstruct_ad.interface"
#include "rttov_calcbt_pc_ad.interface"
#include "rttov_add_raytracing.interface"
#include "rttov_add_aux_prof.interface"
#include "rttov_add_prof.interface"
#include "rttov_init_sunglint.interface"
#include "rttov_init_ircld.interface"
#include "rttov_init_opdp_path.interface"
#include "rttov_init_transmission_aux.interface"
#include "rttov_add_opdp_path.interface"
#include "rttov_init_aux_prof.interface"
#include "rttov_init_auxrad_stream.interface"
#include "rttov_init_trans_scatt_ir.interface"
#include "rttov_alloc_traj_dyn.interface"
#include "rttov_alloc_traj_sta.interface"
#include "rttov_direct.interface"
#include "rttov_nlte_bias_correction_ad.interface"
#include "rttov_calcbt_ad.interface"
!#include "rttov_calcsatrefl_ad.interface"

  LOGICAL(KIND=jplm) :: addcosmic! switch for adding temp of cosmic background

  TYPE(rttov_traj_dyn) :: traj0_dyn
  TYPE(rttov_traj_sta) :: traj0_sta
  TYPE(rttov_traj_dyn) :: traj0_ad_dyn

  LOGICAL(KIND=jplm)   :: ltraj_ad_dyn_dealloc

  TYPE(rttov_traj), TARGET  :: traj1
  TYPE(rttov_traj), POINTER :: traj0
  TYPE(rttov_traj), TARGET  :: traj1_ad
  TYPE(rttov_traj), POINTER :: traj0_ad

  TYPE(rttov_options) :: opts_coef
  INTEGER(KIND=jpim)  :: nlevels
  INTEGER(KIND=jpim)  :: nprofiles
  INTEGER(KIND=jpim)  :: nchannels
  INTEGER(KIND=jpim)  :: npcscores
  INTEGER(KIND=jpim)  :: err
  REAL(KIND=JPRB)     :: ZHOOK_HANDLE

!- End of header --------------------------------------------------------
  TRY
  IF (LHOOK) CALL DR_HOOK('RTTOV_AD', 0_jpim, ZHOOK_HANDLE)
!-------------
!0. initialize
!-------------
  nprofiles                 = SIZE(profiles)
  nchannels                 = SIZE(chanprof)
  opts_coef                 = opts
  opts_coef%rt_ir%addaerosl = .FALSE.
  opts_coef%rt_ir%addclouds = .FALSE.
  nlevels                   = profiles(1)%nlevels
  errorstatus               = errorstatus_success

  ltraj_ad_dyn_dealloc = .FALSE.

  IF (opts%rt_ir%pc%addpc) THEN
    npcscores = SIZE(pccomp%pcscores)
    IF (npcscores / nprofiles > coefs%coef_pccomp%fmv_pc_mnum) THEN
      err = errorstatus_fatal
      THROWM(err.ne.0, "npcscores is larger than maximum allowed by PC coefficient file")
    ENDIF
  ENDIF

  NULLIFY (traj0)
  NULLIFY (traj0_ad)
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
        & traj0_ad = traj0_ad, &
        & traj1 = traj1,       &
        & traj1_ad = traj1_ad, &
        & traj2 = traj,        &
        & traj2_ad = traj_ad)
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

  IF (traj0_sta%dothermal .AND. .NOT. PRESENT(emissivity_ad)) THEN
    err = errorstatus_fatal
    THROWM(err.ne.0, "emissivity_ad parameter required")
  END IF
  IF (traj0_sta%dosolar .AND. .NOT. PRESENT(reflectance_ad)) THEN
    err = errorstatus_fatal
    THROWM(err.ne.0, "reflectance_ad parameter required")
  END IF

  CALL rttov_alloc_traj_dyn (err, traj0_ad_dyn, opts, nchannels, profiles(1)%nlayers, &
                             traj0_dyn%nstreams, ncldtyp, 1_jpim)
  THROW(err.ne.0)
  ltraj_ad_dyn_dealloc = .TRUE.
!
! Adjoint
!----------------
!
!---------------------------------------------
!0. allocate and initialize local AD variables
!---------------------------------------------
!0.1 profile AD-variables on COEF levels (USER levels done in tstrad_ad)
!------------------------------------------------------------------------
  CALL rttov_init_opdp_path(opts, traj0_ad%opdp_path)
  CALL rttov_init_opdp_path(opts, traj0_ad%opdp_path_coef)

  CALL rttov_init_ircld(traj0_ad%ircld)
  IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
    CALL rttov_init_trans_scatt_ir(traj0_ad%transmission_scatt_ir)
    CALL rttov_init_trans_scatt_ir(traj0_ad_dyn%transmission_scatt_ir_stream)
  ENDIF

  CALL rttov_init_auxrad_stream (traj0_ad_dyn%auxrad_stream)

  CALL rttov_init_transmission_aux (opts, traj0_ad_dyn%transmission_aux)

  CALL rttov_init_prof(traj0_ad%profiles_coef)
  CALL rttov_init_prof(traj0_ad%profiles_dry)
!0.2 auxillary profile AD-variables
!----------------------------------
! on USER levels
  CALL rttov_init_aux_prof(traj0_ad%aux_prof)
! on COEF levels
  CALL rttov_init_aux_prof(traj0_ad%aux_prof_coef)
!0.3 predictor AD-variables on COEF levels
!------------------------------------------
  CALL rttov_init_predictor(opts%rt_ir%addsolar, traj0_ad%predictors)
!0.4 raytracing AD-variables
!---------------------------
  CALL rttov_init_raytracing(traj0_sta%dosolar, traj0_ad%raytracing)
  CALL rttov_init_raytracing(traj0_sta%dosolar, traj0_ad%raytracing_coef)
!0.5 path opdep AD-variables on COEF levels
!--------------------------------------------------
! allocated with the direct arrays
! on USER levels
!0.6 transmission AD-variables on USER levels
!--------------------------------------------
! allocated with the direct arrays
!0.7 emissivity and traj0%thermrefl AD-variables - SURFACE
!-------------------------------------------------------
! emissivity_ad is init before calling
  traj0_ad%thermrefl(:) = 0._JPRB

  IF (traj0_sta%dosolar) THEN
    traj0_ad%fresnrefl(:) = 0._JPRB
    CALL rttov_init_sunglint(traj0_ad%sunglint)
  ENDIF

!--------------------------------------------------
! Set up radiance/BT inputs
!--------------------------------------------------

  IF (opts%rt_ir%pc%addpc) THEN

    IF (opts%rt_ir%pc%addradrec) THEN

      IF (opts%rt_all%switchrad) THEN
        CALL rttov_calcbt_pc_ad(         &
              & traj0_sta%chanprof_in,   &
              & coefs%coef_pccomp,       &
              & pccomp,                  &
              & pccomp_ad)
      ENDIF

      CALL rttov_reconstruct_ad( &
            & opts                 ,   &
            & traj0_sta%chanprof_in,   &
            & traj0_sta%chanprof_pc,   &
            & pccomp,                  &
            & pccomp_ad,               &
            & coefs%coef_pccomp)

    ENDIF

    CALL rttov_pcscores_ad( &
          & opts,                   &
          & chanprof,               &
          & traj0_sta%chanprof_pc,  &
          & pccomp,                 &
          & pccomp_ad,              &
          & coefs%coef_pccomp,      &
          & radiancedata_ad)

  ELSE

    ! Classical RTTOV
    IF (opts%rt_all%switchrad) THEN
      IF (traj0_sta%dothermal) &
        CALL rttov_calcbt_ad(chanprof, coefs%coef, traj0_sta%thermal, radiancedata, radiancedata_ad)
!       Input AD perturbation is always in radiance for pure-solar channels, never in reflectance
!       IF (traj0_sta%dosolar) &
!         CALL rttov_calcsatrefl_ad(chanprof, profiles, traj0_sta%solar_spec_esd, &
!                                   traj0_sta%thermal, traj0_sta%solar, radiancedata_ad)
    ELSE
      radiancedata_ad%clear(:) = 0._jprb
    ENDIF

  ENDIF


  ! Do NLTE bias correction for hyperspectral instruments (not PC-RTTOV)
  IF (.NOT. opts%rt_ir%pc%addpc) THEN
    IF (coefs%coef%nltecoef) THEN
      IF (opts%rt_ir%do_nlte_correction .AND. &
          coefs%coef%id_sensor == sensor_id_hi) THEN
        CALL rttov_nlte_bias_correction_ad(coefs%coef, profiles, profiles_ad, &
                                           traj0_sta%angles, chanprof, radiancedata_ad)
      ENDIF
    ENDIF
  ENDIF

!--------------------------------------------------
!1. AD of radiative transfer integration - USER levs
!--------------------------------------------------
  addcosmic = (coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po)

  CALL rttov_integrate_ad( &
        & addcosmic,                                 &
        & opts,                                      &
        & traj0_dyn%nstreams,                        &
        & chanprof,                                  &
        & emissivity,                                &
        & emissivity_ad,                             &
        & reflectance,                               &
        & reflectance_ad,                            &
        & traj0_sta%refl_norm,                       &
        & traj0%thermrefl,                           &
        & traj0_ad%thermrefl,                        &
        & traj0_sta%do_lambertian,                   &
        & traj0_sta%thermal,                         &
        & traj0_sta%dothermal,                       &
        & traj0_sta%solar,                           &
        & traj0_sta%dosolar,                         &
        & traj0_sta%solar_spec_esd,                  &
        & traj0_dyn%transmission_aux,                &
        & traj0_ad_dyn%transmission_aux,             &
        & traj0%transmission_scatt_ir,               &
        & traj0_ad%transmission_scatt_ir,            &
        & profiles,                                  &
        & profiles_ad,                               &
        & traj0%profiles_dry,                        &
        & traj0_ad%profiles_dry,                     &
        & traj0%aux_prof,                            &
        & traj0_ad%aux_prof,                         &
        & coefs%coef,                                &
        & traj0%raytracing,                          &
        & traj0_ad%raytracing,                       &
        & traj0%ircld,                               &
        & traj0_ad%ircld,                            &
        & radiancedata,                              &
        & traj0_sta%auxrad,                          &
        & traj0_dyn%auxrad_stream,                   &
        & traj0_ad_dyn%auxrad_stream,                &
        & radiancedata_ad)! inout  (output if converstion Bt -> rad)

!--------------------------------------------------------
!2. AD of channel reflectances and emissivities - SURFACE
!--------------------------------------------------------

  IF (traj0_sta%dosolar) THEN

    IF (ANY(calcrefl)) THEN

      CALL rttov_calcsurfrefl_ad(     &
            & coefs%coef,             &
            & profiles,               &
            & traj0%sunglint,         &
            & traj0_ad%sunglint,      &
            & traj0%fresnrefl,        &
            & traj0_ad%fresnrefl,     &
            & traj0_sta%solar,        &
            & chanprof,               &
            & traj0_sta%refl_norm,    &
            & calcrefl,               &
            & emissivity,             &
            & emissivity_ad,          &
            & reflectance_ad%refl_out)

      CALL rttov_fresnel_ad( &
            & chanprof,           &
            & profiles,           &
            & traj0_sta%solar,    &
            & coefs%coef,         &
            & traj0%sunglint,     &
            & traj0_ad%sunglint,  &
            & traj0%fresnrefl,    &
            & traj0_ad%fresnrefl)
      CALL rttov_refsun_ad( &
            & profiles,            &
            & profiles_ad,         &
            & coefs%coef,          &
            & traj0%aux_prof,      &
            & traj0%sunglint,      &
            & traj0_ad%sunglint,   &
            & traj0%raytracing,    &
            & traj0_ad%raytracing)

    ENDIF

    reflectance_ad%refl_in = reflectance_ad%refl_in + reflectance_ad%refl_out

  ENDIF


  IF (traj0_sta%dothermal) THEN

    IF (ANY(calcemis)) THEN
      ! Calculate surface emissivity and traj0%thermrefl for selected channels

      IF (coefs%coef%id_sensor == sensor_id_ir) THEN
        ! Infrared
        WHERE (calcemis)
          emissivity_ad(:)%emis_out = emissivity_ad(:)%emis_out - traj0_ad%thermrefl(:)
        ENDWHERE
      ELSE IF (coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po) THEN
        ! Microwave
        CALL rttov_calcemis_mw_ad( &
              & opts,                               &
              & profiles,                           &
              & profiles_ad,                        &
              & traj0_sta%angles,                   &
              & coefs%coef,                         &
              & chanprof,                           &
              & traj0_dyn%transmission_aux,         &
              & traj0_ad_dyn%transmission_aux,      &
              & calcemis,                           &
              & emissivity_ad%emis_out,             &
              & traj0_ad%thermrefl)
      ELSE
        ! Hires
        WHERE (calcemis)
          emissivity_ad(:)%emis_out = emissivity_ad(:)%emis_out - traj0_ad%thermrefl(:)
        ENDWHERE
        CALL rttov_calcemis_ir_ad(    &
              & profiles,             &
              & profiles_ad,          &
              & coefs%coef,           &
              & opts%rt_ir%pc%addpc,  &
              & coefs%coef_pccomp,    &
              & traj0_sta%thermal,    &
              & chanprof,             &
              & calcemis,             &
              & emissivity_ad%emis_out)
      ENDIF

    ENDIF

    IF (coefs%coef%id_sensor == sensor_id_ir .OR. &
        coefs%coef%id_sensor == sensor_id_hi .OR. &
        .NOT. ANY(calcemis)) THEN
      WHERE (.NOT. calcemis)
        emissivity_ad%emis_out = emissivity_ad%emis_out - traj0_ad%thermrefl
      ENDWHERE
    ENDIF

    emissivity_ad%emis_in = emissivity_ad%emis_in + emissivity_ad%emis_out

  ENDIF

!------------------------------------
!3. AD of transmittances - USER levs
!------------------------------------

  IF (traj0_sta%dosolar) THEN
    CALL rttov_transmit_9_solar_ad( &
          & opts%rt_ir%addaerosl,                      &
          & opts%rt_ir%addclouds,                      &
          & profiles(1)%nlayers,                       &
          & nprofiles,                                 &
          & chanprof,                                  &
          & traj0_sta%solar,                           &
          & traj0%aux_prof,                            &
          & traj0_ad%aux_prof,                         &
          & coefs%coef,                                &
          & traj0%raytracing,                          &
          & traj0_ad%raytracing,                       &
          & traj0%ircld,                               &
          & traj0%opdp_path,                           &
          & traj0_ad%opdp_path,                        &
          & traj0_sta%solar_path2,                     &
          & traj0_sta%solar_path1,                     &
          & transmission,                              &
          & transmission_ad,                           &
          & traj0_dyn%transmission_aux,                &
          & traj0_ad_dyn%transmission_aux,             &
          & traj0%transmission_scatt_ir,               &
          & traj0_ad%transmission_scatt_ir,            &
          & traj0_dyn%transmission_scatt_ir_stream,    &
          & traj0_ad_dyn%transmission_scatt_ir_stream)
  ENDIF

  IF (traj0_sta%dothermal) THEN
    CALL rttov_transmit_ad( &
          & opts%rt_ir%addaerosl,                        &
          & opts%rt_ir%addclouds,                        &
          & traj0_sta%do_lambertian,                     &
          & profiles(1)%nlayers,                         &
          & chanprof,                                    &
          & traj0_sta%thermal,                           &
          & traj0%aux_prof,                              &
          & traj0_ad%aux_prof,                           &
          & coefs%coef,                                  &
          & traj0%ircld,                                 &
          & traj0_sta%angles,                            &
          & traj0%opdp_path%atm_level,                   &
          & traj0_ad%opdp_path%atm_level,                &
          & traj0_sta%thermal_path1%od_level,            &
          & transmission%tau_levels,                     &
          & transmission%tau_total,                      &
          & transmission_ad%tau_levels,                  &
          & transmission_ad%tau_total,                   &
          & traj0_dyn%transmission_aux,                  &
          & traj0_dyn%transmission_aux%thermal_path1,    &
          & traj0_ad_dyn%transmission_aux%thermal_path1, &
          & traj0%transmission_scatt_ir,                 &
          & traj0_ad%transmission_scatt_ir,              &
          & traj0_dyn%transmission_scatt_ir_stream,      &
          & traj0_ad_dyn%transmission_scatt_ir_stream,   &
          & traj0_sta%thermal_path1%tau_ref,             &
          & traj0_sta%thermal_path1%tau_ref_surf,        &
          & traj0_sta%thermal_path1%tau_surf,            &
          & traj0_sta%thermal_path1%tau_level)
  ENDIF
!-------------------------------------------------------------
!4. AD of optical depths of aerosols and/or clouds - USER levs
!-------------------------------------------------------------

  IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
    CALL rttov_opdpscattir_ad( &
          & profiles(1)%nlayers,                       &
          & chanprof,                                  &
          & opts,                                      &
          & traj0%aux_prof,                            &
          & traj0_ad%aux_prof,                         &
          & profiles,                                  &
          & profiles_ad,                               &
          & traj0_ad%profiles_dry,                     &
          & aer_opt_param,                             &
          & cld_opt_param,                             &
          & traj0_sta%dosolar,                         &
          & traj0_sta%solar,                           &
          & coefs%coef,                                &
          & coefs%coef_scatt_ir,                       &
          & traj0%raytracing,                          &
          & traj0_ad%raytracing,                       &
          & traj0%transmission_scatt_ir,               &
          & traj0_ad%transmission_scatt_ir,            &
          & traj0_dyn%transmission_scatt_ir_stream,    &
          & traj0_ad_dyn%transmission_scatt_ir_stream, &
          & coefs%optp,                                &
          & traj0%ircld,                               &
          & traj0_ad%ircld)
  ENDIF

!----------------------------------------------------
!5. AD of cloud streams and distributions - USER levs
!----------------------------------------------------

  IF (opts%rt_ir%addclouds) THEN
    CALL rttov_cldstr_ad( &
          & opts%rt_ir,     &
          & profiles,       &
          & profiles_ad,    &
          & traj0%ircld,    &
          & traj0_ad%ircld)
  ENDIF

!------------------------------------
!6. AD of optical depth interpolation
!------------------------------------

  IF (opts%interpolation%spacetop) THEN
    traj0_ad%opdp_path%atm_level(1,:) = 0._jprb
    IF (opts%rt_ir%addsolar) THEN
      traj0_ad%opdp_path%sun_level_path2(1,:) = 0._jprb
    ENDIF
  ENDIF

  IF (opts%interpolation%addinterp) THEN
    CALL rttov_intavg_chan_ad( &
          & opts,                           &
          & traj0_sta%thermal,              &
          & traj0_sta%solar,                &
          & coefs%coef%nlevels,             &
          & nlevels,                        &
          & chanprof,                       &
          & traj0%profiles_coef,            &
          & profiles,                       &
          & profiles_ad,                    &
          & traj0%opdp_path_coef,           &
          & traj0_ad%opdp_path_coef,        &
          & traj0%opdp_path,                &
          & traj0_ad%opdp_path)
  ELSE
    CALL rttov_add_opdp_path(opts, traj0_ad%opdp_path_coef, traj0_ad%opdp_path_coef, traj0_ad%opdp_path)
  ENDIF

!---------------------------------------------------------
!7. AD of atmospheric and solar optical depths - COEF levs
!---------------------------------------------------------
! optical depth arrays allocated at start of AD-code

  IF (coefs%coef%fmv_model_ver == 9) THEN
    IF (traj0_sta%dosolar) THEN
      ! Calculate solar path2 optical depths for solar channels
      CALL rttov_opdep_9_ad( &
            & coefs%coef%nlayers,                      &
            & chanprof,                                &
            & traj0_sta%solar,                         &
            & traj0%predictors,                        &
            & traj0%predictors%path2,                  &
            & traj0_ad%predictors%path2,               &
            & traj0%aux_prof_coef,                     &
            & traj0_ad%aux_prof_coef,                  &
            & coefs%coef,                              &
            & coefs%coef%solar,                        &
            & traj0%opdp_path_coef%sun_level_path2,    &
            & traj0_ad%opdp_path_coef%sun_level_path2, &
            & traj0_sta%solar_path2%opdp_ref_coef)
    ENDIF
    IF (traj0_sta%dothermal) THEN
      ! Calculate thermal path1 optical depths for thermal channels
      CALL rttov_opdep_9_ad( &
            & coefs%coef%nlayers,                  &
            & chanprof,                            &
            & traj0_sta%thermal,                   &
            & traj0%predictors,                    &
            & traj0%predictors%path1,              &
            & traj0_ad%predictors%path1,           &
            & traj0%aux_prof_coef,                 &
            & traj0_ad%aux_prof_coef,              &
            & coefs%coef,                          &
            & coefs%coef%thermal,                  &
            & traj0%opdp_path_coef%atm_level,      &
            & traj0_ad%opdp_path_coef%atm_level,   &
            & traj0_sta%thermal_path1%opdp_ref_coef)
    ENDIF
  ELSE
    CALL rttov_opdep_ad( &
          & coefs%coef%nlayers,        &
          & chanprof,                  &
          & traj0%predictors,          &
          & traj0%predictors%path1,    &
          & traj0_ad%predictors%path1, &
          & traj0%aux_prof_coef,       &
          & traj0_ad%aux_prof_coef,    &
          & coefs%coef,                &
          & coefs%coef%thermal,        &
          & traj0%opdp_path_coef,      &
          & traj0_ad%opdp_path_coef,   &
          & traj0_sta%thermal_path1%opdp_ref_coef)
  ENDIF

!-------------------------------------------------------
!8. AD of RTTOV-7 RTTOV-8 RTTOV-9 predictors - COEF levs
!-------------------------------------------------------

  IF (coefs%coef%fmv_model_ver == 7) THEN
    CALL rttov_setpredictors_7_ad( &
          & opts,                      &
          & traj0%profiles_coef,       &
          & traj0_ad%profiles_coef,    &
          & traj0_sta%angles_coef,     &
          & coefs%coef,                &
          & traj0%aux_prof_coef,       &
          & traj0_ad%aux_prof_coef,       &
          & traj0%predictors%path1,    &
          & traj0_ad%predictors%path1, &
          & traj0%raytracing_coef,     &
          & traj0_ad%raytracing_coef)
  ELSE IF (coefs%coef%fmv_model_ver == 8) THEN
    CALL rttov_setpredictors_8_ad( &
          & opts,                      &
          & traj0%profiles_coef,       &
          & traj0_ad%profiles_coef,    &
          & traj0_sta%angles_coef,     &
          & coefs%coef,                &
          & traj0%aux_prof_coef,       &
          & traj0%predictors%path1,    &
          & traj0_ad%predictors%path1, &
          & traj0%raytracing_coef,     &
          & traj0_ad%raytracing_coef)
  ELSE IF (coefs%coef%fmv_model_ver == 9) THEN
    IF (traj0_sta%dosolar) THEN
      CALL rttov_setpredictors_9_ad(            &
            & opts,                             &
            & traj0%profiles_coef,              &
            & traj0_ad%profiles_coef,           &
            & traj0_sta%angles_coef,            &
            & traj0%raytracing_coef%patheff,    &
            & traj0_ad%raytracing_coef%patheff, &
            & coefs%coef_pccomp,                &
            & coefs%coef,                       &
            & traj0%predictors%path2,           &
            & traj0_ad%predictors%path2)
    ENDIF
    IF (traj0_sta%dothermal) THEN
      CALL rttov_setpredictors_9_ad(            &
            & opts,                             &
            & traj0%profiles_coef,              &
            & traj0_ad%profiles_coef,           &
            & traj0_sta%angles_coef,            &
            & traj0%raytracing_coef%pathsat,    &
            & traj0_ad%raytracing_coef%pathsat, &
            & coefs%coef_pccomp,                &
            & coefs%coef,                       &
            & traj0%predictors%path1,           &
            & traj0_ad%predictors%path1)
    ENDIF
  ENDIF

!-------------------------------------------------------
!9. AD of common geometric set-up for the RT integration
!-------------------------------------------------------

!9.1 USER levs
!-------------
  IF (opts%interpolation%addinterp) THEN
    CALL rttov_setgeometry_ad( &
      & opts,                &
      & traj0_sta%dosolar,   &
      & profiles,            &
      & profiles_ad,         &
      & traj0%aux_prof,      &
      & coefs%coef,          &
      & traj0_sta%angles,    &
      & traj0%raytracing,    &
      & traj0_ad%raytracing, &
      & traj0%profiles_dry,  &
      & traj0_ad%profiles_dry)
  ELSE
    CALL rttov_add_raytracing(traj0_sta%dosolar, traj0_ad%raytracing_coef, &
                              traj0_ad%raytracing_coef, traj0_ad%raytracing)
  ENDIF

!9.1 COEF levs
!-------------
  CALL rttov_setgeometry_ad( &
    & opts,                     &
    & traj0_sta%dosolar,        &
    & traj0%profiles_coef,      &
    & traj0_ad%profiles_coef,   &
    & traj0%aux_prof_coef,      &
    & coefs%coef,               &
    & traj0_sta%angles_coef,    &
    & traj0%raytracing_coef,    &
    & traj0_ad%raytracing_coef, &
    & traj0%profiles_coef,      &
    & traj0_ad%profiles_coef,   &
    & .TRUE._jplm)  ! do_pmc_calc on coef levels

!---------------------------------------------------------
!10. AD of cloud top, surface levels, ice cloud parameters
!---------------------------------------------------------
!10.2 USER levs
!-------------
  IF (opts%interpolation%addinterp .or. opts%rt_ir%addclouds) THEN
    CALL rttov_profaux_ad( &
      & opts,              &
      & profiles,          &
      & profiles_ad,       &
      & coefs%coef,        &
      & traj0%aux_prof,    &
      & traj0_ad%aux_prof)
  ELSE
    CALL rttov_add_aux_prof(traj0_ad%aux_prof_coef, traj0_ad%aux_prof_coef, traj0_ad%aux_prof)
  ENDIF

!10.1 COEF levs
!-------------
    CALL rttov_profaux_ad( &
          & opts_coef,              &
          & traj0%profiles_coef,    &
          & traj0_ad%profiles_coef, &
          & coefs%coef,             &
          & traj0%aux_prof_coef,    &
          & traj0_ad%aux_prof_coef)


!------------------------------------------------------------------
!11 AD of apply_reg_limits - COEF levs
!------------------------------------------------------------------

  IF (opts%config%apply_reg_limits .OR. opts%interpolation%reg_limit_extrap) THEN
    CALL rttov_apply_reg_limits_ad( &
          & opts,                        &
          & profiles,                    &
          & traj0_sta%profiles_coef_ref, &
          & traj0_ad%profiles_coef,      &
          & coefs%coef,                  &
          & coefs%coef_pccomp)
  ENDIF

!-----------------------------------------------------------------
!12. AD of profile variable interpolation - COEF levs to USER levs
!-----------------------------------------------------------------

  IF (opts%interpolation%addinterp) THEN
    CALL rttov_intavg_prof_ad( &
          & opts,                           &
          & nlevels,                        &
          & coefs%coef%nlevels,             &
          & profiles,                       &
          & profiles_ad,                    &
          & traj0%profiles_dry,             &
          & traj0_ad%profiles_dry,          &
          & traj0%profiles_coef,            &
          & traj0_ad%profiles_coef,         &
          & coefs%coef,                     &
          & coefs%coef_pccomp)
  ELSE
! COEF levs same as USER levs
    CALL rttov_add_prof( &
          & profiles_ad,            &
          & profiles_ad,            &
          & traj0_ad%profiles_coef, &
          & lair = .TRUE._jplm,     &
          & lground = .FALSE._jplm, &
          & profiles_gas = traj0_ad%profiles_dry)
  ENDIF

  CALL rttov_add_prof( &
        & profiles_ad,            &
        & profiles_ad,            &
        & traj0_ad%profiles_coef, &
        & lair = .FALSE._jplm,    &
        & lground = .TRUE._jplm, &
        & profiles_gas = traj0_ad%profiles_dry)

!-----------------------------------------------------------------------
! Convert input profiles to ppmv wrt dry air
!-----------------------------------------------------------------------

  CALL rttov_convert_gas_dry_ad( &
      opts,                      &
      profiles,                  &
      profiles_ad,               &
      traj0%profiles_dry,        &
      traj0_ad%profiles_dry)


!-----------------------------------------
!13. deallocate memory for local variables
!-----------------------------------------

  CALL cleanup()

  IF (LHOOK) CALL DR_HOOK('RTTOV_AD', 1_jpim, ZHOOK_HANDLE)
  CATCH_C
  errorstatus = err
  IF (LHOOK) CALL DR_HOOK('RTTOV_AD', 1_jpim, ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE cleanup()
    INTEGER(KIND=jpim) :: error

    IF (traj0_dyn%nstreams >= 0) THEN
      CALL rttov_alloc_traj_dyn (error, traj0_dyn, opts, nchannels, profiles(1)%nlayers, &
                                 traj0_dyn%nstreams, ncldtyp, 0_jpim)
    ENDIF

    IF (ltraj_ad_dyn_dealloc) THEN
      CALL rttov_alloc_traj_dyn (error, traj0_ad_dyn, opts, nchannels, profiles(1)%nlayers, &
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
            & traj0_ad = traj0_ad, &
            & traj1 = traj1,       &
            & traj1_ad = traj1_ad, &
            & traj2 = traj,        &
            & traj2_ad = traj_ad)
    ENDIF
  END SUBROUTINE cleanup

END SUBROUTINE rttov_ad
