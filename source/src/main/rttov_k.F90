! Description:
!> @file
!!   Runs RTTOV Jacobian (K) model
!
!> @brief
!!   Runs RTTOV Jacobian (K) model
!!
!! @details
!!   Given an input profile, calculates the Jacobian of the direct
!!   model with the outputs being the derivatives of the channel
!!   radiances or BTs with respect to each profile variable.
!!
!!   All K variables should be initialised to zero before calling.
!!   An input radiance or BT perturbation should be specified. This is
!!   usually 1.0: the value linearly scales the corresponding output
!!   Jacobians.
!!
!!   The input perturbations should then be written to the radiancedata_k
!!   structure. If opts\%rt_all\%switchrad is false the perturbations are
!!   in radiance in radiancedata_k%total, otherwise they are in BT in
!!   radiancedata_k\%bt. For visible/near-IR channels (at wavelengths less
!!   than 3 microns) the input perturbations are always in
!!   radiancedata_k\%total.
!!
!!   The output Jacobians are contained in profiles_k.
!!
!!   For PC-RTTOV the input perturbations are specified in pccomp_k instead.
!!   If opts\%rt_ir\%addpc\%addradrec is false the perturbations are in PC
!!   scores in pccomp_k\%pcscores. Otherwise they are in
!!   pccomp_k\%total_pccomp or pccomp_k\%bt_pccomp according to the
!!   setting of opts\%rt_all\%switchrad (false/true respectively).
!!
!!   For PC-RTTOV profiles_k contains the Jacobians of the PC predictor
!!   channels. The Jacobians of PC scores wrt profile variables are in
!!   profiles_k_pc while the Jacobians of reconstructed radiances or BTs
!!   are contained in profiles_k_rec.
!!
!! @param[out]    errorstatus     status on exit
!! @param[in]     chanprof        specifies channels and profiles to simulate
!! @param[in]     opts            options to configure the simulations
!! @param[in]     profiles        input atmospheric profiles and surface variables
!! @param[in,out] profiles_k      output Jacobians wrt atmospheric profile and surface variables
!! @param[in]     coefs           coefficients structure for instrument to simulate
!! @param[in,out] transmission    output transmittances
!! @param[in,out] transmission_k  input perturbations wrt transmittances (usually zero)
!! @param[in,out] radiancedata    output radiances and corresponding BTs and BRFs
!! @param[in,out] radiancedata_k  input perturbations wrt radiances or BTs
!! @param[in]     calcemis        flags for internal RTTOV surface emissivity calculation, optional
!! @param[in,out] emissivity      input/output surface emissivities, optional
!! @param[in,out] emissivity_k    output Jacobians wrt surface emissivities, optional
!! @param[in]     calcrefl        flags for internal RTTOV surface BRDF calculation, optional
!! @param[in,out] reflectance     input/output surface BRDFs, input cloud top BRDF for simple cloud, optional
!! @param[in,out] reflectance_k   output Jacobians wrt surface BRDFs, optional
!! @param[in]     aer_opt_param   input aerosol optical parameters, optional
!! @param[in]     cld_opt_param   input cloud optical parameters, optional
!! @param[in,out] traj            RTTOV direct internal state, can be initialised outside RTTOV, optional
!! @param[in,out] traj_k          RTTOV K internal state, can be initialised outside RTTOV, optional
!! @param[in,out] pccomp          output PC scores and radiances from PC-RTTOV, optional
!! @param[in,out] pccomp_k        input gradient in PC scores, radiances or BTs for PC-RTTOV, optional
!! @param[in,out] profiles_k_pc   output PC score Jacobians wrt atmospheric profile and surface variables, optional
!! @param[in,out] profiles_k_rec  output reconstructed radiance/BT Jacobians wrt atmospheric profile and surface
!!                                  variables, optional
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
SUBROUTINE rttov_k( &
            & errorstatus,      &
            & chanprof,         &
            & opts,             &
            & profiles,         &
            & profiles_k,       &
            & coefs,            &
            & transmission,     &
            & transmission_k,   &
            & radiancedata,     &
            & radiancedata_k,   &
            & calcemis,         &
            & emissivity,       &
            & emissivity_k,     &
            & calcrefl,         &
            & reflectance,      &
            & reflectance_k,    &
            & aer_opt_param,    &
            & cld_opt_param,    &
            & traj,             &
            & traj_k,           &
            & pccomp,           &
            & pccomp_k,         &
            & profiles_k_pc,    &
            & profiles_k_rec,   &
            & channels_rec)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY :  &
       & rttov_options,     &
       & rttov_coefs,       &
       & rttov_pccomp,      &
       & transmission_type, &
       & profile_type,      &
       & radiance_type,     &
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
       & rttov_traj_sta,  &
       & rttov_traj_dyn
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim)     , INTENT(OUT)                     :: errorstatus
  TYPE(profile_type  )   , INTENT(IN)                      :: profiles(:)
  TYPE(rttov_chanprof)   , INTENT(IN)                      :: chanprof(:)
  TYPE(rttov_options )   , INTENT(IN)                      :: opts
  TYPE(profile_type  )   , INTENT(INOUT)                   :: profiles_k(SIZE(chanprof))
  TYPE(rttov_coefs   )   , INTENT(IN)   , TARGET           :: coefs
  TYPE(transmission_type), INTENT(INOUT)                   :: transmission
  TYPE(transmission_type), INTENT(INOUT)                   :: transmission_k
  TYPE(radiance_type    ), INTENT(INOUT)                   :: radiancedata
  TYPE(radiance_type    ), INTENT(INOUT)                   :: radiancedata_k
  LOGICAL(KIND=jplm)     , INTENT(IN)   , OPTIONAL         :: calcemis(SIZE(chanprof))
  TYPE(rttov_emissivity) , INTENT(INOUT), OPTIONAL         :: emissivity(SIZE(chanprof))
  TYPE(rttov_emissivity) , INTENT(INOUT), OPTIONAL         :: emissivity_k(SIZE(chanprof))
  LOGICAL(KIND=jplm)     , INTENT(IN)   , OPTIONAL         :: calcrefl(SIZE(chanprof))
  TYPE(rttov_reflectance), INTENT(INOUT), OPTIONAL         :: reflectance(SIZE(chanprof))
  TYPE(rttov_reflectance), INTENT(INOUT), OPTIONAL         :: reflectance_k(SIZE(chanprof))
  TYPE(rttov_opt_param)  , INTENT(IN)   , OPTIONAL         :: aer_opt_param
  TYPE(rttov_opt_param)  , INTENT(IN)   , OPTIONAL         :: cld_opt_param
  TYPE(rttov_traj       ), INTENT(INOUT), OPTIONAL, TARGET :: traj, traj_k        ! Target is needed here
  TYPE(rttov_pccomp     ), INTENT(INOUT), OPTIONAL         :: pccomp
  TYPE(rttov_pccomp     ), INTENT(INOUT), OPTIONAL         :: pccomp_k
  TYPE(profile_type     ), INTENT(INOUT), OPTIONAL         :: profiles_k_pc(:)
  TYPE(profile_type     ), INTENT(INOUT), OPTIONAL         :: profiles_k_rec(:)
  INTEGER(KIND=jpim)     , INTENT(IN)   , OPTIONAL         :: channels_rec(:)
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_convert_gas_dry_k.interface"
#include "rttov_apply_reg_limits_k.interface"
#include "rttov_setgeometry_k.interface"
#include "rttov_calcemis_ir_k.interface"
#include "rttov_calcsurfrefl_k.interface"
#include "rttov_profaux_k.interface"
#include "rttov_setpredictors_7_k.interface"
#include "rttov_setpredictors_8_k.interface"
#include "rttov_setpredictors_9_k.interface"
#include "rttov_opdep_k.interface"
#include "rttov_opdep_9_k.interface"
#include "rttov_transmit_k.interface"
#include "rttov_transmit_9_solar_k.interface"
#include "rttov_calcemis_mw_k.interface"
#include "rttov_integrate_k.interface"
#include "rttov_intavg_prof_k.interface"
#include "rttov_intavg_chan_k.interface"
#include "rttov_refsun_k.interface"
#include "rttov_init_prof.interface"
#include "rttov_init_predictor.interface"
#include "rttov_init_raytracing.interface"
#include "rttov_check_traj.interface"
#include "rttov_cldstr_k.interface"
#include "rttov_opdpscattir_k.interface"
#include "rttov_fresnel_k.interface"
#include "rttov_mult_profiles_k.interface"
#include "rttov_pcscores_k.interface"
#include "rttov_pcscores_rec_k.interface"
#include "rttov_reconstruct_k.interface"
#include "rttov_calcbt_pc_ad.interface"
#include "rttov_add_raytracing.interface"
#include "rttov_add_aux_prof.interface"
#include "rttov_add_prof.interface"
#include "rttov_init_rad.interface"
#include "rttov_init_sunglint.interface"
#include "rttov_init_ircld.interface"
#include "rttov_init_opdp_path.interface"
#include "rttov_add_opdp_path.interface"
#include "rttov_init_aux_prof.interface"
#include "rttov_init_auxrad_stream.interface"
#include "rttov_init_trans_scatt_ir.interface"
#include "rttov_init_transmission_aux.interface"
#include "rttov_alloc_traj_dyn.interface"
#include "rttov_alloc_traj_sta.interface"
#include "rttov_direct.interface"
#include "rttov_nlte_bias_correction_k.interface"
#include "rttov_calcbt_ad.interface"
!#include "rttov_calcsatrefl_ad.interface"

  INTEGER(KIND=jpim) :: nlevels
  LOGICAL(KIND=jplm) :: addcosmic! switch for adding temp of cosmic background

  TYPE(rttov_traj_dyn) :: traj0_dyn
  TYPE(rttov_traj_sta) :: traj0_sta
  TYPE(rttov_traj_dyn) :: traj0_k_dyn

  LOGICAL(KIND=jplm)   :: ltraj_k_dyn_dealloc

  TYPE(rttov_traj), TARGET  :: traj1
  TYPE(rttov_traj), POINTER :: traj0
  TYPE(rttov_traj), TARGET  :: traj1_k
  TYPE(rttov_traj), POINTER :: traj0_k

  TYPE(rttov_options) :: opts_coef
  INTEGER(KIND=jpim) :: npcscores
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: nchannels_rec
  REAL(KIND=jprb)     , ALLOCATABLE :: total_k_pc (:, :, :), pcscores_k(:, :, :)
  INTEGER(KIND=jpim)  :: err
  REAL   (KIND=JPRB)  :: ZHOOK_HANDLE

!- End of header ------------------------------------------------------
  TRY
  IF (LHOOK) CALL DR_HOOK('RTTOV_K', 0_jpim, ZHOOK_HANDLE)
!-------------
!0. initialize
!-------------
  nprofiles                 = SIZE(profiles)
  nchannels                 = SIZE(chanprof)
  nlevels                   = profiles(1)%nlevels
  opts_coef                 = opts
  opts_coef%rt_ir%addclouds = .FALSE.
  opts_coef%rt_ir%addaerosl = .FALSE.
  errorstatus               = errorstatus_success

  ltraj_k_dyn_dealloc = .FALSE.

  IF (opts%rt_ir%pc%addpc) THEN
    npcscores = SIZE(pccomp%pcscores)
    IF (npcscores / nprofiles > coefs%coef_pccomp%fmv_pc_mnum) THEN
      err = errorstatus_fatal
      THROWM(err.ne.0, "npcscores is larger than maximum allowed by PC coefficient file")
    ENDIF
  ENDIF

  NULLIFY (traj0)
  NULLIFY (traj0_k)
!-----------------------------------------------------------------------
!1. interpolator first call - input profiles from USER levs to COEF levs
!-----------------------------------------------------------------------
! do immediately - for quick escape if target values are out of bounds
  CALL rttov_check_traj( &
        & err,               &
        & nprofiles,         &
        & nchannels,         &
        & opts,              &
        & nlevels,           &
        & coefs,             &
        & 1_jpim,            &
        & traj0 = traj0,     &
        & traj0_k = traj0_k, &
        & traj1 = traj1,     &
        & traj1_k = traj1_k, &
        & traj2 = traj,      &
        & traj2_k = traj_k)
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

  IF (traj0_sta%dothermal .AND. .NOT. PRESENT(emissivity_k)) THEN
    err = errorstatus_fatal
    THROWM(err.ne.0, "emissivity_k parameter required")
  END IF
  IF (traj0_sta%dosolar .AND. .NOT. PRESENT(reflectance_k)) THEN
    err = errorstatus_fatal
    THROWM(err.ne.0, "reflectance_k parameter required")
  END IF

  CALL rttov_alloc_traj_dyn (err, traj0_k_dyn, opts, nchannels, profiles(1)%nlayers, &
                             traj0_dyn%nstreams, ncldtyp, 1_jpim)
  THROW(err.ne.0)
  ltraj_k_dyn_dealloc = .TRUE.
! K matrix
!----------------
!---------------------------------------------
!0. allocate and initialize local AD variables
!---------------------------------------------
!0.1 auxillary profile K-variables
!----------------------------------
!0.6 path opdep K-variables on COEF levels
!-----------------------------------------
  CALL rttov_init_opdp_path(opts, traj0_k%opdp_path)
  CALL rttov_init_opdp_path(opts, traj0_k%opdp_path_coef)

  CALL rttov_init_ircld(traj0_k%ircld)
  IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
    CALL rttov_init_trans_scatt_ir(traj0_k%transmission_scatt_ir)
    CALL rttov_init_trans_scatt_ir(traj0_k_dyn%transmission_scatt_ir_stream)
  ENDIF

  CALL rttov_init_auxrad_stream (traj0_k_dyn%auxrad_stream)

  CALL rttov_init_transmission_aux (opts, traj0_k_dyn%transmission_aux)

  CALL rttov_init_prof(traj0_k%profiles_coef)
  CALL rttov_init_prof(traj0_k%profiles_dry)

!0.4 raytracing K-variables
!--------------------------
  CALL rttov_init_raytracing(traj0_sta%dosolar, traj0_k%raytracing)
  CALL rttov_init_raytracing(traj0_sta%dosolar, traj0_k%raytracing_coef)

! on USER levels
  CALL rttov_init_aux_prof(traj0_k%aux_prof)
! on COEF levels
  CALL rttov_init_aux_prof(traj0_k%aux_prof_coef)
!0.2 profile_all K-variables on USER levels
!------------------------------------------

  traj0_k%thermrefl(:) = 0._JPRB

  IF (traj0_sta%dosolar) THEN
    traj0_k%fresnrefl(:) = 0._JPRB
    CALL rttov_init_sunglint(traj0_k%sunglint)
  ENDIF

!--------------------------------------------------
! Set up radiance/BT inputs
!--------------------------------------------------
  IF (opts%rt_ir%pc%addpc) THEN

    CALL rttov_init_rad(radiancedata_k)
    radiancedata_k%total = 1._jprb

  ELSE

    ! Not PC-RTTOV
    IF (opts%rt_all%switchrad) THEN
      IF (traj0_sta%dothermal) &
        CALL rttov_calcbt_ad(chanprof, coefs%coef, traj0_sta%thermal, radiancedata, radiancedata_k)
!       Input K perturbation is always in radiance for pure-solar channels, never in reflectance
!       IF (traj0_sta%dosolar) &
!         CALL rttov_calcsatrefl_ad(chanprof, profiles, traj0_sta%solar_spec_esd, &
!                                   traj0_sta%thermal, traj0_sta%solar, radiancedata_k)
    ELSE
      radiancedata_k%clear(:) = 0._jprb
    ENDIF

  ENDIF


! Do NLTE bias correction for hyperspectral instruments (not PC-RTTOV)
  IF (.NOT. opts%rt_ir%pc%addpc) THEN
    IF (coefs%coef%nltecoef) THEN
      IF (opts%rt_ir%do_nlte_correction .AND. &
          coefs%coef%id_sensor == sensor_id_hi) THEN
        CALL rttov_nlte_bias_correction_k(coefs%coef, profiles, profiles_k, &
                                          traj0_sta%angles, chanprof, &
                                          radiancedata_k)
      ENDIF
    ENDIF
  ENDIF

!--------------------------------------------------
!1. K of radiative transfer integration - USER levs
!--------------------------------------------------
  addcosmic = (coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po)

  CALL rttov_integrate_k( &
        & addcosmic,                                &
        & opts,                                     &
        & traj0_dyn%nstreams,                       &
        & chanprof,                                 &
        & emissivity,                               &
        & emissivity_k,                             &
        & reflectance,                              &
        & reflectance_k,                            &
        & traj0_sta%refl_norm,                      &
        & traj0%thermrefl,                          &
        & traj0_k%thermrefl,                        &
        & traj0_sta%do_lambertian,                  &
        & traj0_sta%thermal,                        &
        & traj0_sta%dothermal,                      &
        & traj0_sta%solar,                          &
        & traj0_sta%dosolar,                        &
        & traj0_sta%solar_spec_esd,                 &
        & traj0_dyn%transmission_aux,               &
        & traj0_k_dyn%transmission_aux,             &
        & traj0%transmission_scatt_ir,              &
        & traj0_k%transmission_scatt_ir,            &
        & profiles,                                 &
        & profiles_k,                               &
        & traj0%profiles_dry,                       &
        & traj0_k%profiles_dry,                     &
        & traj0%aux_prof,                           &
        & traj0_k%aux_prof,                         &
        & coefs%coef,                               &
        & traj0%raytracing,                         &
        & traj0_k%raytracing,                       &
        & traj0%ircld,                              &
        & traj0_k%ircld,                            &
        & radiancedata,                             &
        & traj0_sta%auxrad,                         &
        & traj0_dyn%auxrad_stream,                  &
        & traj0_k_dyn%auxrad_stream,                &
        & radiancedata_k)

!--------------------------------------------------------
!2. K of channel reflectances and emissivities - SURFACE
!--------------------------------------------------------

  IF (traj0_sta%dosolar) THEN

    IF (ANY(calcrefl)) THEN

      CALL rttov_calcsurfrefl_k(     &
            & coefs%coef,            &
            & profiles,              &
            & traj0%sunglint,        &
            & traj0_k%sunglint,      &
            & traj0%fresnrefl,       &
            & traj0_k%fresnrefl,     &
            & traj0_sta%solar,       &
            & chanprof,              &
            & traj0_sta%refl_norm,   &
            & calcrefl,              &
            & emissivity,            &
            & emissivity_k,          &
            & reflectance_k%refl_out)

      CALL rttov_fresnel_k( &
            & chanprof,          &
            & profiles,          &
            & traj0_sta%solar,   &
            & coefs%coef,        &
            & traj0%sunglint,    &
            & traj0_k%sunglint,  &
            & traj0%fresnrefl,   &
            & traj0_k%fresnrefl)
      CALL rttov_refsun_k( &
            & chanprof,           &
            & profiles,           &
            & profiles_k,         &
            & coefs%coef,         &
            & traj0%aux_prof,     &
            & traj0%sunglint,     &
            & traj0_k%sunglint,   &
            & traj0%raytracing,   &
            & traj0_k%raytracing)

    ENDIF

    reflectance_k%refl_in = reflectance_k%refl_in + reflectance_k%refl_out

  ENDIF


  IF (traj0_sta%dothermal) THEN

    IF (ANY(calcemis)) THEN
      ! Calculate surface emissivity and traj0%thermrefl for selected channels

      IF (coefs%coef%id_sensor == sensor_id_ir) THEN
        ! Infrared
        WHERE (calcemis)
          emissivity_k(:)%emis_out = emissivity_k(:)%emis_out - traj0_k%thermrefl(:)
        ENDWHERE
      ELSE IF (coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po) THEN
        ! Microwave
        CALL rttov_calcemis_mw_k( &
              & opts,                           &
              & profiles,                       &
              & profiles_k,                     &
              & traj0_sta%angles,               &
              & coefs%coef,                     &
              & chanprof,                       &
              & traj0_dyn%transmission_aux,     &
              & traj0_k_dyn%transmission_aux,   &
              & calcemis,                       &
              & emissivity_k%emis_out,          &
              & traj0_k%thermrefl)
      ELSE
        ! Hires
        WHERE (calcemis)
          emissivity_k(:)%emis_out = emissivity_k(:)%emis_out - traj0_k%thermrefl(:)
        ENDWHERE
        CALL rttov_calcemis_ir_k( &
              & profiles,             &
              & profiles_k,           &
              & coefs%coef,           &
              & opts%rt_ir%pc%addpc,  &
              & coefs%coef_pccomp,    &
              & traj0_sta%thermal,    &
              & chanprof,             &
              & calcemis,             &
              & emissivity_k%emis_out)
      ENDIF

    ENDIF

    IF (coefs%coef%id_sensor == sensor_id_ir .OR. &
        coefs%coef%id_sensor == sensor_id_hi .OR. &
        .NOT. ANY(calcemis)) THEN
      WHERE (.NOT. calcemis)
        emissivity_k%emis_out = emissivity_k%emis_out - traj0_k%thermrefl
      ENDWHERE
    ENDIF

    emissivity_k%emis_in = emissivity_k%emis_in + emissivity_k%emis_out

  ENDIF
!----------------------------------
!3. K of transmittances - USER levs
!----------------------------------

  IF (traj0_sta%dosolar) THEN
    CALL rttov_transmit_9_solar_k( &
          & opts%rt_ir%addaerosl,                     &
          & opts%rt_ir%addclouds,                     &
          & profiles(1)%nlayers,                      &
          & chanprof,                                 &
          & traj0_sta%solar,                          &
          & traj0%aux_prof,                           &
          & traj0_k%aux_prof,                         &
          & coefs%coef,                               &
          & traj0%raytracing,                         &
          & traj0_k%raytracing,                       &
          & traj0%ircld,                              &
          & traj0%opdp_path,                          &
          & traj0_k%opdp_path,                        &
          & traj0_sta%solar_path2,                    &
          & traj0_sta%solar_path1,                    &
          & transmission,                             &
          & transmission_k,                           &
          & traj0_dyn%transmission_aux,               &
          & traj0_k_dyn%transmission_aux,             &
          & traj0%transmission_scatt_ir,              &
          & traj0_k%transmission_scatt_ir,            &
          & traj0_dyn%transmission_scatt_ir_stream,   &
          & traj0_k_dyn%transmission_scatt_ir_stream)
  ENDIF

  IF (traj0_sta%dothermal) THEN
    CALL rttov_transmit_k( &
          & opts%rt_ir%addaerosl,                       &
          & opts%rt_ir%addclouds,                       &
          & traj0_sta%do_lambertian,                    &
          & profiles(1)%nlayers,                        &
          & chanprof,                                   &
          & traj0_sta%thermal,                          &
          & traj0%aux_prof,                             &
          & traj0_k%aux_prof,                           &
          & coefs%coef,                                 &
          & traj0%ircld,                                &
          & traj0_sta%angles,                           &
          & traj0%opdp_path%atm_level,                  &
          & traj0_k%opdp_path%atm_level,                &
          & traj0_sta%thermal_path1%od_level,           &
          & transmission%tau_levels,                    &
          & transmission%tau_total,                     &
          & transmission_k%tau_levels,                  &
          & transmission_k%tau_total,                   &
          & traj0_dyn%transmission_aux,                 &
          & traj0_dyn%transmission_aux%thermal_path1,   &
          & traj0_k_dyn%transmission_aux%thermal_path1, &
          & traj0%transmission_scatt_ir,                &
          & traj0_k%transmission_scatt_ir,              &
          & traj0_dyn%transmission_scatt_ir_stream,     &
          & traj0_k_dyn%transmission_scatt_ir_stream,   &
          & traj0_sta%thermal_path1%tau_ref,            &
          & traj0_sta%thermal_path1%tau_ref_surf,       &
          & traj0_sta%thermal_path1%tau_surf,           &
          & traj0_sta%thermal_path1%tau_level)
  ENDIF
!------------------------------------------------------------
!4. K of optical depths of aerosols and/or clouds - USER levs
!------------------------------------------------------------

  IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
    CALL rttov_opdpscattir_k( &
          & profiles(1)%nlayers,                       &
          & chanprof,                                  &
          & opts,                                      &
          & traj0%aux_prof,                            &
          & traj0_k%aux_prof,                          &
          & profiles,                                  &
          & profiles_k,                                &
          & traj0_k%profiles_dry,                      &
          & aer_opt_param,                             &
          & cld_opt_param,                             &
          & traj0_sta%dosolar,                         &
          & traj0_sta%solar,                           &
          & coefs%coef,                                &
          & coefs%coef_scatt_ir,                       &
          & traj0%raytracing,                          &
          & traj0_k%raytracing,                        &
          & traj0%transmission_scatt_ir,               &
          & traj0_k%transmission_scatt_ir,             &
          & traj0_dyn%transmission_scatt_ir_stream,    &
          & traj0_k_dyn%transmission_scatt_ir_stream,  &
          & coefs%optp,                                &
          & traj0%ircld,                               &
          & traj0_k%ircld)
  ENDIF

!---------------------------------------------------
!5. K of cloud streams and distributions - USER levs
!---------------------------------------------------

  IF (opts%rt_ir%addclouds) THEN
    CALL rttov_cldstr_k( &
          & opts%rt_ir,     &
          & chanprof,       &
          & profiles,       &
          & profiles_k,     &
          & traj0%ircld,    &
          & traj0_k%ircld)
  ENDIF

!-----------------------------------
!6. K of optical depth interpolation
!-----------------------------------

  IF (opts%interpolation%spacetop) THEN
    traj0_k%opdp_path%atm_level(1,:) = 0._jprb
    IF (opts%rt_ir%addsolar) THEN
      traj0_k%opdp_path%sun_level_path2(1,:) = 0._jprb
    ENDIF
  ENDIF

  IF (opts%interpolation%addinterp) THEN
    CALL rttov_intavg_chan_k( &
          & opts,                           &
          & traj0_sta%thermal,              &
          & traj0_sta%solar,                &
          & coefs%coef%nlevels,             &
          & nlevels,                        &
          & chanprof,                       &
          & traj0%profiles_coef,            &
          & profiles,                       &
          & profiles_k,                     &
          & traj0%opdp_path_coef,           &
          & traj0_k%opdp_path_coef,         &
          & traj0%opdp_path,                &
          & traj0_k%opdp_path)
  ELSE
    CALL rttov_add_opdp_path(opts, traj0_k%opdp_path_coef, traj0_k%opdp_path_coef, traj0_k%opdp_path)
  ENDIF

!--------------------------------------------------------
!7. K of atmospheric and solar optical depths - COEF levs
!--------------------------------------------------------
! optical depth arrays allocated at start of K-code

  IF (coefs%coef%fmv_model_ver == 9) THEN
! DARFIX: TEST removing init stuff - only done for rttov7 predictors so far
    CALL rttov_init_predictor(opts%rt_ir%addsolar, traj0_k%predictors)
    IF (traj0_sta%dosolar) THEN
      ! Calculate solar path2 optical depths for solar channels
      CALL rttov_opdep_9_k( &
            & coefs%coef%nlayers,                     &
            & chanprof,                               &
            & traj0_sta%solar,                        &
            & traj0%predictors,                       &
            & traj0%predictors%path2,                 &
            & traj0_k%predictors%path2,               &
            & traj0%aux_prof_coef,                    &
            & traj0_k%aux_prof_coef,                  &
            & coefs%coef,                             &
            & coefs%coef%solar,                       &
            & traj0%opdp_path_coef%sun_level_path2,   &
            & traj0_k%opdp_path_coef%sun_level_path2, &
            & traj0_sta%solar_path2%opdp_ref_coef)
    ENDIF
    IF (traj0_sta%dothermal) THEN
      ! Calculate thermal path1 optical depths for thermal channels
      CALL rttov_opdep_9_k( &
            & coefs%coef%nlayers,                 &
            & chanprof,                           &
            & traj0_sta%thermal,                  &
            & traj0%predictors,                   &
            & traj0%predictors%path1,             &
            & traj0_k%predictors%path1,           &
            & traj0%aux_prof_coef,                &
            & traj0_k%aux_prof_coef,              &
            & coefs%coef,                         &
            & coefs%coef%thermal,                 &
            & traj0%opdp_path_coef%atm_level,     &
            & traj0_k%opdp_path_coef%atm_level,   &
            & traj0_sta%thermal_path1%opdp_ref_coef)
    ENDIF
  ELSE
    CALL rttov_opdep_k( &
          & coefs%coef%nlayers,       &
          & chanprof,                 &
          & traj0%predictors,         &
          & traj0%predictors%path1,   &
          & traj0_k%predictors%path1, &
          & traj0%aux_prof_coef,      &
          & traj0_k%aux_prof_coef,    &
          & coefs%coef,               &
          & coefs%coef%thermal,       &
          & traj0%opdp_path_coef,     &
          & traj0_k%opdp_path_coef,   &
          & traj0_sta%thermal_path1%opdp_ref_coef)
  ENDIF

!-------------------------------------------------------
!8. K of RTTOV-7 RTTOV-8 RTTOV-9 predictors - COEF levs
!-------------------------------------------------------

  IF (coefs%coef%fmv_model_ver == 7) THEN
    CALL rttov_setpredictors_7_k( &
          & opts,                     &
          & chanprof,                 &
          & traj0%profiles_coef,      &
          & traj0_k%profiles_coef,    &
          & traj0_sta%angles_coef,    &
          & coefs%coef,               &
          & traj0%aux_prof_coef,      &
          & traj0_k%aux_prof_coef,    &
          & traj0%predictors%path1,   &
          & traj0_k%predictors%path1, &
          & traj0%raytracing_coef,    &
          & traj0_k%raytracing_coef)
  ELSE IF (coefs%coef%fmv_model_ver == 8) THEN
    CALL rttov_setpredictors_8_k( &
          & opts,                     &
          & coefs%coef%nlayers,       &
          & traj0_sta%angles_coef,    &
          & chanprof,                 &
          & traj0%profiles_coef,      &
          & traj0_k%profiles_coef,    &
          & coefs%coef,               &
          & traj0%aux_prof_coef,      &
          & traj0%predictors%path1,   &
          & traj0_k%predictors%path1, &
          & traj0%raytracing_coef,    &
          & traj0_k%raytracing_coef)
  ELSE IF (coefs%coef%fmv_model_ver == 9) THEN
    IF (traj0_sta%dosolar) THEN
      CALL rttov_setpredictors_9_k(            &
            & opts,                            &
            & chanprof,                        &
            & coefs%coef%nlayers,              &
            & traj0%profiles_coef,             &
            & traj0_k%profiles_coef,           &
            & traj0_sta%angles_coef,           &
            & traj0%raytracing_coef%patheff,   &
            & traj0_k%raytracing_coef%patheff, &
            & coefs%coef_pccomp,               &
            & coefs%coef,                      &
            & traj0%predictors%path2,          &
            & traj0_k%predictors%path2)
    ENDIF
    IF (traj0_sta%dothermal) THEN
      CALL rttov_setpredictors_9_k(            &
            & opts,                            &
            & chanprof,                        &
            & coefs%coef%nlayers,              &
            & traj0%profiles_coef,             &
            & traj0_k%profiles_coef,           &
            & traj0_sta%angles_coef,           &
            & traj0%raytracing_coef%pathsat,   &
            & traj0_k%raytracing_coef%pathsat, &
            & coefs%coef_pccomp,               &
            & coefs%coef,                      &
            & traj0%predictors%path1,          &
            & traj0_k%predictors%path1)
    ENDIF
  ENDIF

!------------------------------------------------------
!9. K of common geometric set-up for the RT integration
!------------------------------------------------------
!9.1 USER levs
!-------------
  IF (opts%interpolation%addinterp) THEN
    CALL rttov_setgeometry_k( &
      & opts,               &
      & traj0_sta%dosolar,  &
      & chanprof,           &
      & profiles,           &
      & profiles_k,         &
      & traj0%aux_prof,     &
      & coefs%coef,         &
      & traj0_sta%angles,   &
      & traj0%raytracing,   &
      & traj0_k%raytracing, &
      & traj0%profiles_dry, &
      & traj0_k%profiles_dry)
  ELSE
    CALL rttov_add_raytracing(traj0_sta%dosolar, traj0_k%raytracing_coef, &
                              traj0_k%raytracing_coef, traj0_k%raytracing)
  ENDIF

!9.1 COEF levs
!-------------
    CALL rttov_setgeometry_k( &
          & opts,                    &
          & traj0_sta%dosolar,       &
          & chanprof,                &
          & traj0%profiles_coef,     &
          & traj0_k%profiles_coef,   &
          & traj0%aux_prof_coef,     &
          & coefs%coef,              &
          & traj0_sta%angles_coef,   &
          & traj0%raytracing_coef,   &
          & traj0_k%raytracing_coef, &
          & traj0%profiles_coef,     &
          & traj0_k%profiles_coef,   &
          & .TRUE._jplm)  ! do_pmc_calc on coef levels

!---------------------------------------------------------
!10. AD of cloud top, surface levels, ice cloud parameters
!---------------------------------------------------------
!10.1 COEF levs
!-------------

  IF (opts%interpolation%addinterp .or. opts%rt_ir%addclouds) THEN
    CALL rttov_profaux_k( &
      & opts,             &
      & chanprof,         &
      & profiles,         &
      & profiles_k,       &
      & coefs%coef,       &
      & traj0%aux_prof,   &
      & traj0_k%aux_prof)
  ELSE
    CALL rttov_add_aux_prof(traj0_k%aux_prof_coef, traj0_k%aux_prof_coef, traj0_k%aux_prof)
  ENDIF

!10.2 USER levs
!-------------
    CALL rttov_profaux_k( &
          & opts_coef,             &
          & chanprof,              &
          & traj0%profiles_coef,   &
          & traj0_k%profiles_coef, &
          & coefs%coef,            &
          & traj0%aux_prof_coef,   &
          & traj0_k%aux_prof_coef)

!------------------------------------------------------------------
!11 K of apply_reg_limits - COEF levs
!------------------------------------------------------------------

  IF (opts%config%apply_reg_limits .OR. opts%interpolation%reg_limit_extrap) THEN
    CALL rttov_apply_reg_limits_k( &
          & opts,                         &
          & chanprof,                     &
          & profiles,                     &
          & traj0_sta%profiles_coef_ref,  &
          & traj0_k%profiles_coef,        &
          & coefs%coef,                   &
          & coefs%coef_pccomp)
  ENDIF

!----------------------------------------------------------------
!12. K of profile variable interpolation - COEF levs to USER levs
!----------------------------------------------------------------

  IF (opts%interpolation%addinterp) THEN
    CALL rttov_intavg_prof_k( &
          & opts,                           &
          & chanprof,                       &
          & nlevels,                        &
          & coefs%coef%nlevels,             &
          & profiles,                       &
          & profiles_k,                     &
          & traj0%profiles_dry,             &
          & traj0_k%profiles_dry,           &
          & traj0%profiles_coef,            &
          & traj0_k%profiles_coef,          &
          & coefs%coef,                     &
          & coefs%coef_pccomp)
  ELSE
    CALL rttov_add_prof( &
          & profiles_k,             &
          & profiles_k,             &
          & traj0_k%profiles_coef,  &
          & lair = .TRUE._jplm,     &
          & lground = .FALSE._jplm, &
          & profiles_gas = traj0_k%profiles_dry)
  ENDIF

  CALL rttov_add_prof( &
        & profiles_k,            &
        & profiles_k,            &
        & traj0_k%profiles_coef, &
        & lair = .FALSE._jplm,   &
        & lground = .TRUE._jplm, &
        & profiles_gas = traj0_k%profiles_dry)

!-----------------------------------------------------------------------
! Convert input profiles to ppmv wrt dry air
!-----------------------------------------------------------------------

  CALL rttov_convert_gas_dry_k(  &
      opts,                      &
      chanprof,                  &
      profiles,                  &
      profiles_k,                &
      traj0%profiles_dry,        &
      traj0_k%profiles_dry)


  IF (opts%rt_ir%pc%addpc) THEN
    npcscores     = SIZE(traj0_sta%chanprof_pc)

    IF (opts%rt_ir%pc%addradrec) THEN
      nchannels_rec = SIZE(traj0_sta%chanprof_in)

      ALLOCATE ( &
        pcscores_k(npcscores / nprofiles, nchannels_rec / nprofiles, nprofiles),&
        total_k_pc(nchannels / nprofiles, nchannels_rec / nprofiles, nprofiles),&
        STAT = err)
      THROWM( err .NE. 0 , "allocation pc k arrays")

      IF (opts%rt_all%switchrad) THEN
        CALL rttov_calcbt_pc_ad( &
          & traj0_sta%chanprof_in,       &
          & coefs%coef_pccomp,           &
          & pccomp,                      &
          & pccomp_k)
      ENDIF

      CALL rttov_reconstruct_k( &
        & opts,                        &
        & traj0_sta%chanprof_in,       &
        & traj0_sta%chanprof_pc,       &
        & pccomp,                      &
        & pccomp_k,                    &
        & pcscores_k,                  &
        & coefs%coef_pccomp)

      CALL rttov_pcscores_rec_k(       &
        & opts,                        &
        & chanprof,                    &
        & pcscores_k,                  &
        & coefs%coef_pccomp,           &
        & total_k_pc)

! DAR: Don't need to initialise profiles if they're being overwritten immediately
!     CALL rttov_init_prof(profiles_k_rec)
      CALL rttov_mult_profiles_k(profiles_k_rec, profiles_k, total_k_pc, opts)

      DEALLOCATE (pcscores_k, total_k_pc, STAT = err)
      THROWM( err .NE. 0 , "deallocation pc k arrays")

    ELSE
      ALLOCATE ( &
        total_k_pc(nchannels / nprofiles, npcscores / nprofiles, nprofiles), &
        STAT = err)
      THROWM( err .NE. 0 , "allocation pc k arrays")

      CALL rttov_pcscores_k( &
            & opts,                        &
            & chanprof,                    &
            & traj0_sta%chanprof_pc,       &
            & pccomp,                      &
            & pccomp_k,                    &
            & coefs%coef_pccomp,           &
            & total_k_pc)

! DAR: Don't need to initialise profiles if they're being overwritten immediately
!     CALL rttov_init_prof(profiles_k_pc)
      CALL rttov_mult_profiles_k(profiles_k_pc, profiles_k, total_k_pc, opts)

      DEALLOCATE (total_k_pc, STAT = err)
      THROWM( err .NE. 0 , "deallocation pc k arrays")

    ENDIF
  ENDIF


!-----------------------------------------
!13. deallocate memory for local variables
!-----------------------------------------

  CALL cleanup()

  IF (LHOOK) CALL DR_HOOK('RTTOV_K', 1_jpim, ZHOOK_HANDLE)
  CATCH_C
  errorstatus = err
  IF (LHOOK) CALL DR_HOOK('RTTOV_K', 1_jpim, ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE cleanup()
    INTEGER(KIND=jpim) :: error

    IF (traj0_dyn%nstreams >= 0) THEN
      CALL rttov_alloc_traj_dyn (error, traj0_dyn, opts, nchannels, profiles(1)%nlayers, &
                                 traj0_dyn%nstreams, ncldtyp, 0_jpim)
    ENDIF

    IF (ltraj_k_dyn_dealloc) THEN
      CALL rttov_alloc_traj_dyn (error, traj0_k_dyn, opts, nchannels, profiles(1)%nlayers, &
                                 traj0_dyn%nstreams, ncldtyp, 0_jpim)
    ENDIF

    IF (ASSOCIATED(traj0_sta%thermal)) THEN
      CALL rttov_alloc_traj_sta (error, traj0_sta, opts, coefs%coef, nlevels, nchannels, nprofiles, &
                                 0_jpim, npcscores, channels_rec)
    ENDIF

    IF (ASSOCIATED(traj0)) THEN
      CALL rttov_check_traj( &
            & error,             &
            & nprofiles,         &
            & nchannels,         &
            & opts,              &
            & nlevels,           &
            & coefs,             &
            & 0_jpim,            &
            & traj0 = traj0,     &
            & traj0_k = traj0_k, &
            & traj1 = traj1,     &
            & traj1_k = traj1_k, &
            & traj2 = traj,      &
            & traj2_k = traj_k)
    ENDIF
  END SUBROUTINE cleanup

END SUBROUTINE rttov_k
