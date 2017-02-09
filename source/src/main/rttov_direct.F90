! Description:
!> @file
!!   Runs RTTOV direct model
!
!> @brief
!!   Runs RTTOV direct model
!!
!! @details
!!   Computes multi-channel level to space transmittances,
!!   top of atmosphere and level to space radiances, brightness
!!   temperatures and reflectances, and optionally surface
!!   emissivities and BRDFs for multiple profiles in a single call,
!!   for nadir-viewing satellite-based visible, infrared or microwave
!!   sensors. Requires a coefficient file for the sensor for which
!!   simulated radiances are requested.
!!
!!   The methodology is described in the following:
!!
!!   Eyre J.R. and H.M. Woolf  1988 Transmittance of atmospheric gases
!!   in the microwave region: a fast model. Applied Optics 27  3244-3249
!!
!!   Eyre J.R. 1991 A fast radiative transfer model for satellite sounding
!!   systems.  ECMWF Research Dept. Tech. Memo. 176 (available from the
!!   librarian at ECMWF).
!!
!!   Saunders R.W., M. Matricardi and P. Brunel 1999 An Improved Fast Radiative
!!   Transfer Model for Assimilation of Satellite Radiance Observations.
!!   QJRMS, 125, 1407-1425.
!!
!!   Matricardi, M., F. Chevallier and S. Tjemkes 2001 An improved general
!!   fast radiative transfer model for the assimilation of radiance
!!   observations. ECMWF Research Dept. Tech. Memo. 345
!!   (available from the librarian at ECMWF).
!!
!!   Matricardi, M. 2003 RTIASI-4, a new version of the ECMWF fast radiative
!!   transfer model for the infrared atmospheric sounding interferometer.
!!   ECMWF Research Dept. Tech. Memo. 425 (available from the librarian at ECMWF)
!!
!!   Rochon Y.J., L. Garand, D.S. Turner and S. Polavarapu Jacobian mapping
!!   between  vertical coordinate systems in data assimilation (submitted QJRMS
!!   June 2006)
!!
!!   Matricardi, M. 2009: An Observation operator for the assimilation of
!!   principal component scores into a NWP system. Available from EUMETSAT
!!
!! @param[out]    errorstatus    status on exit
!! @param[in]     chanprof       specifies channels and profiles to simulate
!! @param[in]     opts           options to configure the simulations
!! @param[in]     profiles       input atmospheric profiles and surface variables
!! @param[in]     coefs          coefficients structure for instrument to simulate
!! @param[in,out] transmission   output transmittances
!! @param[in,out] radiancedata   output radiances and corresponding BTs and BRFs
!! @param[in,out] radiancedata2  secondary output radiances, optional
!! @param[in]     calcemis       flags for internal RTTOV surface emissivity calculation, optional
!! @param[in,out] emissivity     input/output surface emissivities, optional
!! @param[in]     calcrefl       flags for internal RTTOV surface BRDF calculation, optional
!! @param[in,out] reflectance    input/output surface BRDFs, input cloud top BRDF for simple cloud, optional
!! @param[in]     aer_opt_param  input aerosol optical parameters, optional
!! @param[in]     cld_opt_param  input cloud optical parameters, optional
!! @param[in,out] traj           RTTOV internal state, can be initialised outside RTTOV, optional
!! @param[in,out] pccomp         output PC scores and radiances from PC-RTTOV, optional
!! @param[in]     channels_rec   list of channels for which to calculate reconstructed radiances, optional
!! @param[in,out] traj_sta       RTTOV internal state, optional, not intended for general use
!! @param[in,out] traj_dyn       RTTOV internal state, optional, not intended for general use
!! @param[in,out] lbl_check      used for coef verification, optional, not intended for general use
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
SUBROUTINE rttov_direct( &
            & errorstatus,    &
            & chanprof,       &
            & opts,           &
            & profiles,       &
            & coefs,          &
            & transmission,   &
            & radiancedata,   &
            & radiancedata2,  &
            & calcemis,       &
            & emissivity,     &
            & calcrefl,       &
            & reflectance,    &
            & aer_opt_param,  &
            & cld_opt_param,  &
            & traj,           &
            & traj_dyn,       &
            & traj_sta,       &
            & pccomp,         &
            & channels_rec,   &
            & lbl_check)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY :   &
       & rttov_coefs,       &
       & rttov_pccomp,      &
       & profile_type,      &
       & transmission_type, &
       & radiance_type,     &
       & radiance2_type,    &
       & rttov_options,     &
       & rttov_chanprof,    &
       & rttov_emissivity,  &
       & rttov_reflectance, &
       & rttov_opt_param,   &
       & rttov_traj,        &
       & rttov_traj_dyn,    &
       & rttov_traj_sta,    &
       & rttov_lbl_check
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE parkind1, ONLY : jprb
  USE rttov_const, ONLY :  &
       & sensor_id_mw, &
       & sensor_id_ir, &
       & sensor_id_hi, &
       & sensor_id_po, &
       & ncldtyp,      &
       & max_sol_zen,  &
       & surftype_sea, &
       & gas_unit_compatibility
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim)     , INTENT(OUT)                     :: errorstatus     ! return flag
  TYPE(profile_type  )   , INTENT(IN)                      :: profiles(:)     ! Atmospheric profiles supplied on
                                                                              ! user levels (nprofiles)
  TYPE(rttov_chanprof)   , INTENT(IN)                      :: chanprof(:)     ! Profile/channel indices (nchanprof)
  TYPE(rttov_options )   , INTENT(IN)                      :: opts
  TYPE(rttov_coefs   )   , INTENT(IN)   , TARGET           :: coefs           ! It is necessary to have "Target"
                                                                              ! attribute here
  TYPE(transmission_type), INTENT(INOUT)                   :: transmission    ! transmittances
  TYPE(radiance_type    ), INTENT(INOUT)                   :: radiancedata    ! radiances (mw/cm-1/ster/sq.m),
                                                                              ! BTs (degK) and reflectances (BRF)
  TYPE(radiance2_type   ), INTENT(INOUT), OPTIONAL         :: radiancedata2
  LOGICAL(KIND=jplm)     , INTENT(IN)   , OPTIONAL         :: calcemis(SIZE(chanprof))    ! switches for emis calcs
  TYPE(rttov_emissivity) , INTENT(INOUT), OPTIONAL         :: emissivity(SIZE(chanprof))  ! surface emis
  LOGICAL(KIND=jplm)     , INTENT(IN)   , OPTIONAL         :: calcrefl(SIZE(chanprof))    ! switches for refl calcs
  TYPE(rttov_reflectance), INTENT(INOUT), OPTIONAL         :: reflectance(SIZE(chanprof)) ! surface refl
  TYPE(rttov_opt_param)  , INTENT(IN)   , OPTIONAL         :: aer_opt_param
  TYPE(rttov_opt_param)  , INTENT(IN)   , OPTIONAL         :: cld_opt_param
  TYPE(rttov_traj  )     , INTENT(INOUT), OPTIONAL, TARGET :: traj            ! Target is *NEEDED* here (see
                                                                              ! rttov_check_traj)
  TYPE(rttov_traj_dyn)   , INTENT(INOUT), OPTIONAL, TARGET :: traj_dyn
  TYPE(rttov_traj_sta)   , INTENT(INOUT), OPTIONAL, TARGET :: traj_sta
  TYPE(rttov_pccomp)     , INTENT(INOUT), OPTIONAL         :: pccomp
  INTEGER(KIND=jpim)     , INTENT(IN)   , OPTIONAL         :: channels_rec(:)
  TYPE(rttov_lbl_check)  , INTENT(IN)   , OPTIONAL         :: lbl_check
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_convert_gas_dry.interface"
#include "rttov_cldstr.interface"
#include "rttov_apply_reg_limits.interface"
#include "rttov_checkinput.interface"
#include "rttov_profaux.interface"
#include "rttov_setgeometry.interface"
#include "rttov_setpredictors_7.interface"
#include "rttov_setpredictors_8.interface"
#include "rttov_setpredictors_9.interface"
#include "rttov_opdpscattir.interface"
#include "rttov_fresnel.interface"
#include "rttov_opdep.interface"
#include "rttov_opdep_9.interface"
#include "rttov_transmit.interface"
#include "rttov_transmit_9_solar.interface"
#include "rttov_calcemis_ir.interface"
#include "rttov_calcemis_mw.interface"
#include "rttov_calcsurfrefl.interface"
#include "rttov_integrate.interface"
#include "rttov_intavg_chan.interface"
#include "rttov_intavg_prof.interface"
#include "rttov_refsun.interface"
#include "rttov_check_traj.interface"
#include "rttov_calc_solar_spec_esd.interface"
#include "rttov_init_prof.interface"
#include "rttov_init_raytracing.interface"
#include "rttov_reconstruct.interface"
#include "rttov_pcscores.interface"
#include "rttov_calcbt_pc.interface"
#include "rttov_copy_raytracing.interface"
#include "rttov_copy_prof.interface"
#include "rttov_copy_aux_prof.interface"
#include "rttov_copy_opdp_path.interface"
#include "rttov_alloc_traj_dyn.interface"
#include "rttov_alloc_traj_sta.interface"
#include "rttov_checkpcchan.interface"
#include "rttov_nlte_bias_correction.interface"
#include "rttov_calcbt.interface"
#include "rttov_calcsatrefl.interface"

  INTEGER(KIND=jpim) :: i, j, prof, chan
  INTEGER(KIND=jpim) :: nlevels
  LOGICAL(KIND=jplm) :: addcosmic ! switch for adding temp of cosmic background

  TYPE(rttov_options) :: opts_coef

  TYPE(rttov_traj), TARGET  :: traj1
  TYPE(rttov_traj), POINTER :: traj0
  TYPE(rttov_traj_dyn), TARGET  :: traj1_dyn
  TYPE(rttov_traj_dyn), POINTER :: traj0_dyn
  TYPE(rttov_traj_sta), TARGET  :: traj1_sta
  TYPE(rttov_traj_sta), POINTER :: traj0_sta
  LOGICAL(KIND=jplm) :: ltraj_dyn_dealloc
  LOGICAL(KIND=jplm) :: ltraj_sta_dealloc
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: npcscores
  INTEGER(KIND=jpim) :: err

  REAL(KIND=jprb)     :: ZHOOK_HANDLE

!- End of header ------------------------------------------------------

  TRY
!-------------
!0. initialize
!-------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_DIRECT', 0_jpim, ZHOOK_HANDLE)

  nprofiles                 = SIZE(profiles)
  nchannels                 = SIZE(chanprof)
  nlevels                   = profiles(1)%nlevels
  opts_coef                 = opts
  opts_coef%rt_ir%addclouds = .FALSE._jplm
  opts_coef%rt_ir%addaerosl = .FALSE._jplm
  errorstatus               = errorstatus_success

  ltraj_dyn_dealloc   = .FALSE.
  ltraj_sta_dealloc   = .FALSE.
  NULLIFY (traj0, traj0_dyn, traj0_sta)

  IF (opts%rt_ir%user_aer_opt_param .AND. .NOT. PRESENT(aer_opt_param)) THEN
    err = errorstatus_fatal
    THROWM(err.ne.0, "if opts%rt_ir%user_aer_opt_param is TRUE then aer_opt_param is mandatory")
  ENDIF

  IF (opts%rt_ir%user_cld_opt_param .AND. .NOT. PRESENT(cld_opt_param)) THEN
    err = errorstatus_fatal
    THROWM(err.ne.0, "if opts%rt_ir%user_cld_opt_param is TRUE then cld_opt_param is mandatory")
  ENDIF

  IF (MAXVAL(chanprof(:)%chan) > coefs%coef%fmv_chn) THEN
    err = errorstatus_fatal
    THROWM(err.ne.0, "Channel index out of range for coefficients")
  ENDIF

  IF (ANY(profiles(:)%gas_units /= profiles(1)%gas_units)) THEN
    err = errorstatus_fatal
    THROWM(err.ne.0, "All profiles(:)%gas_units must be the same")
  ENDIF

  IF (opts%rt_ir%pc%addpc) THEN
    npcscores = SIZE(pccomp%pcscores)
    IF (npcscores / nprofiles > coefs%coef_pccomp%fmv_pc_mnum) THEN
      err = errorstatus_fatal
      THROWM(err.ne.0, "npcscores is larger than maximum allowed by PC coefficient file")
    ENDIF

    CALL rttov_checkpcchan( &
       & nprofiles,         &
       & nchannels,         &
       & opts,              &
       & chanprof,          &
       & coefs,             &
       & err                )
    THROWM(err .NE. 0, "rttov_checkpcchan fatal error")
  ENDIF

  CALL rttov_check_traj( &
        & err,           &
        & nprofiles,     &
        & nchannels,     &
        & opts,          &
        & nlevels,       &
        & coefs,         &
        & 1_jpim,        &
        & traj0 = traj0, &
        & traj1 = traj1, &
        & traj2 = traj)
  THROWM(err .NE. 0, "rttov_check_traj fatal error")

  IF (PRESENT(traj_sta)) THEN
    traj0_sta => traj_sta
  ELSE
    traj0_sta => traj1_sta
  ENDIF
  ltraj_sta_dealloc = .NOT. PRESENT(traj_sta)

  CALL rttov_alloc_traj_sta (err, traj0_sta, opts, coefs%coef, nlevels, nchannels, nprofiles, &
                            & 1_jpim, npcscores, channels_rec)
  THROW(err .NE. 0)

!----------------------------------------------------------------------------------
! Determine whether thermal and/or solar calculations are required for each channel
!----------------------------------------------------------------------------------

! Set flags to indicate if thermal and/or solar calculations should be performed for each channel.
! If RTTOV is called for solar-only channels and addsolar is FALSE, the resulting radiances will be zero.

  traj0_sta%thermal(:) = .TRUE._jplm
  traj0_sta%solar(:)   = .FALSE._jplm

  ! Turn off thermal calculations for any pure-solar channels
  IF (ASSOCIATED(coefs%coef%ss_val_chn)) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      traj0_sta%thermal(j) = .NOT. (coefs%coef%ss_val_chn(chan) == 2_jpim)
    ENDDO
  ENDIF

  traj0_sta%dothermal = ANY(traj0_sta%thermal(:))
  traj0_sta%dosolar   = .FALSE._jplm

  IF (opts%rt_ir%addsolar) THEN ! If addsolar is FALSE no solar calculations will be done at all

    IF ((coefs%coef%id_sensor == sensor_id_ir .OR. coefs%coef%id_sensor == sensor_id_hi) .AND. &
         coefs%coef%fmv_model_ver == 9) THEN

      IF (ASSOCIATED(coefs%coef%ss_val_chn)) THEN
        DO j = 1, nchannels
          prof = chanprof(j)%prof
          chan = chanprof(j)%chan

          IF (profiles(prof)%sunzenangle >= 0.0_jprb .AND. &
              profiles(prof)%sunzenangle < max_sol_zen) THEN

            ! ss_val_chn = 0 => thermal-only
            ! ss_val_chn = 1 => thermal + solar
            ! ss_val_chn = 2 => solar-only
            traj0_sta%solar(j) = (coefs%coef%ss_val_chn(chan) > 0_jpim)

          ENDIF
        ENDDO
      ENDIF

      traj0_sta%dosolar = ANY(traj0_sta%solar(:))

    ENDIF
  ENDIF

  IF (traj0_sta%dothermal .AND. .NOT. (PRESENT(calcemis) .AND. PRESENT(emissivity))) THEN
    err = errorstatus_fatal
    THROWM(err.ne.0, "calcemis and emissivity parameters required")
  END IF
  IF (traj0_sta%dosolar .AND. .NOT. (PRESENT(calcrefl) .AND. PRESENT(reflectance))) THEN
    err = errorstatus_fatal
    THROWM(err.ne.0, "calcrefl and reflectance parameters required")
  END IF

  IF (traj0_sta%dosolar) CALL rttov_calc_solar_spec_esd(coefs%coef, chanprof, profiles, traj0_sta%solar_spec_esd)

  ! Determine whether Lambertian option should be switched on for each channel
  ! For v11.3 we allow this for IR, but must retain the opts%rt_mw%do_lambertian option for compatibility.
  ! In v11.3 we add an rt_ir% option and an rt_all% option.
  ! In v12 opts%rt_all%do_lambertian will be the only option and the rt_ir% and rt_mw% options will be removed

  traj0_sta%do_lambertian(:) = .FALSE.
  IF (coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po) THEN
    IF (PRESENT(calcemis)) THEN
      IF (opts%rt_mw%do_lambertian .OR. opts%rt_all%do_lambertian) THEN
        ! FASTEM takes this into account already so do not use this option with FASTEM
        DO i = 1, nchannels
          traj0_sta%do_lambertian(i) = (profiles(chanprof(i)%prof)%skin%surftype /= surftype_sea) .OR. .NOT. calcemis(i)
        END DO
      END IF
    END IF
  ELSE
    IF (.NOT. opts%rt_ir%pc%addpc) THEN
      ! Lambertian option should not be used with PC-RTTOV
      traj0_sta%do_lambertian(:) = (opts%rt_ir%do_lambertian .OR. opts%rt_all%do_lambertian)
    ENDIF
  END IF

!-----------------------------------------------------------------------
! Convert input profiles to ppmv wrt dry air
!-----------------------------------------------------------------------

  CALL rttov_convert_gas_dry( &
      opts,                   &
      profiles,               &
      traj0%profiles_dry)

!-----------------------------------------------------------------------
!1. interpolator first call - input profiles from USER levs to COEF levs
!-----------------------------------------------------------------------

  CALL rttov_init_prof(traj0%profiles_coef, p = coefs%coef%ref_prfl_p)

  IF (opts%interpolation%addinterp) THEN
    CALL rttov_intavg_prof( &
          & opts,                           &
          & nlevels,                        &
          & coefs%coef%nlevels,             &
          & profiles,                       &
          & traj0%profiles_dry,             &
          & traj0%profiles_coef,            &
          & coefs%coef,                     &
          & coefs%coef_pccomp)
  ELSE
    CALL rttov_copy_prof( &
          & traj0%profiles_coef,    &
          & profiles,               &
          & larray = .TRUE._jplm,   &
          & lscalar = .FALSE._jplm, &
          & profiles_gas = traj0%profiles_dry)
  ENDIF

  CALL rttov_copy_prof( &
        & traj0%profiles_coef,   &
        & profiles,              &
        & larray = .FALSE._jplm, &
        & lscalar = .TRUE._jplm, &
        & profiles_gas = traj0%profiles_dry)

  WHERE (profiles(:)%gas_units /= gas_unit_compatibility)
    traj0%profiles_coef(:)%gas_units = -1  ! ppmv dry
  ENDWHERE

!------------------------------------------------------------------
!2. check input data is within suitable physical limits - COEF levs
!------------------------------------------------------------------

  IF (opts%config%do_checkinput) THEN
    CALL rttov_checkinput( &
           & opts,                   &
           & traj0%profiles_coef,    &
           & profiles,               &
           & coefs%coef,             &
           & coefs%coef_pccomp,      &
           & err)
    THROW(err .NE. 0)
  ENDIF

  ! Apply the regression limits: do this after checkinput so that unphysical
  ! inputs are reported as failures
  IF (opts%config%apply_reg_limits .OR. opts%interpolation%reg_limit_extrap) THEN
    CALL rttov_copy_prof(traj0_sta%profiles_coef_ref, traj0%profiles_coef)

    CALL rttov_apply_reg_limits( &
           & opts,                &
           & profiles,            &
           & traj0%profiles_coef, &
           & coefs%coef,          &
           & coefs%coef_pccomp)
  ENDIF


!------------------------------------------------------------------------
!3. determine cloud top, surface levels, ice cloud parameters
!------------------------------------------------------------------------

!3.2 COEF levs
!-------------

  CALL rttov_profaux( &
    & opts_coef,           &
    & traj0%profiles_coef, &
    & coefs%coef,          &
    & traj0%aux_prof_coef, &
    & on_coef_levels = LOGICAL((coefs%coef%fmv_model_ver == 7), KIND=jplm))

!3.1 USER levs
!-------------
  IF (opts%interpolation%addinterp .OR. opts%rt_ir%addclouds) THEN
    CALL rttov_profaux( &
      & opts,           &
      & profiles,       &
      & coefs%coef,     &
      & traj0%aux_prof, on_coef_levels = .FALSE._jplm) 
  ELSE
    CALL rttov_copy_aux_prof(traj0%aux_prof, traj0%aux_prof_coef)
  ENDIF

!------------------------------------------------------------------------
!4. set up common geometric variables for the RT integration
!------------------------------------------------------------------------

!4.2 COEF levs
!-------------
  CALL rttov_setgeometry( &
    & opts,                  &
    & traj0_sta%dosolar,     &
    & traj0%profiles_coef,   &
    & traj0%aux_prof_coef,   &
    & coefs%coef,            &
    & traj0_sta%angles_coef, &
    & traj0%raytracing_coef, &
    & traj0%profiles_coef,   &
    & .TRUE._jplm)  ! do_pmc_calc on coef levels

!4.1 USER levs
!-------------
  CALL rttov_init_raytracing(traj0_sta%dosolar, traj0%raytracing)
  IF (opts%interpolation%addinterp) THEN
    CALL rttov_setgeometry( &
      & opts,               &
      & traj0_sta%dosolar,  &
      & profiles,           &
      & traj0%aux_prof,     &
      & coefs%coef,         &
      & traj0_sta%angles,   &
      & traj0%raytracing,   &
      & traj0%profiles_dry)
  ELSE
    traj0_sta%angles = traj0_sta%angles_coef
    CALL rttov_copy_raytracing(traj0_sta%dosolar, traj0%raytracing, traj0%raytracing_coef)
  ENDIF


  IF (PRESENT(lbl_check)) THEN
!
! For lbl checks, we may want to test in a plane geometry
!
    IF( lbl_check%plane_geometry ) THEN
      DO i = 1, nprofiles
        traj0%raytracing%pathsat(:,i) = traj0_sta%angles(i)%seczen
        traj0%raytracing%zasat(:,i)   = traj0_sta%angles(i)%sinview
      ENDDO
      traj0%raytracing%pathsat_sqrt  = SQRT(traj0%raytracing%pathsat)
      traj0%raytracing%pathsat_rsqrt = 1.0_jprb/traj0%raytracing%pathsat_sqrt

      DO i = 1, nprofiles
        traj0%raytracing_coef%pathsat(:,i) = traj0_sta%angles_coef(i)%seczen
        traj0%raytracing_coef%zasat(:,i)   = traj0_sta%angles_coef(i)%sinview
      ENDDO
      traj0%raytracing_coef%pathsat_sqrt  = SQRT(traj0%raytracing_coef%pathsat)
      traj0%raytracing_coef%pathsat_rsqrt = 1.0_jprb/traj0%raytracing_coef%pathsat_sqrt

      IF (traj0_sta%dosolar) THEN
        DO i = 1, nprofiles
          traj0%raytracing%zasun(:,i)   = traj0_sta%angles(i)%sinzen_sun
          traj0%raytracing%pathsun(:,i) = &
            & 1.0_jprb/SQRT(1.0_jprb - traj0_sta%angles(i)%sinzen_sun**2_jpim)
          traj0%raytracing%patheff(:,i) = traj0%raytracing%pathsat(:,i) + traj0%raytracing%pathsun(:,i)
        ENDDO

        DO i = 1, nprofiles
          traj0%raytracing_coef%zasun(:,i)   = traj0_sta%angles_coef(i)%sinzen_sun
          traj0%raytracing_coef%pathsun(:,i) = &
            & 1.0_jprb/SQRT(1.0_jprb - traj0_sta%angles_coef(i)%sinzen_sun**2_jpim)
          traj0%raytracing_coef%patheff(:,i) = traj0%raytracing_coef%pathsat(:,i) + traj0%raytracing_coef%pathsun(:,i)
        ENDDO
      ENDIF
    ENDIF

  ENDIF

!---------------------------------------------------
!5. calculate atm/sol predictors - atm/sol COEF levs
!---------------------------------------------------

  IF (coefs%coef%fmv_model_ver == 7) THEN
    CALL rttov_setpredictors_7( &
          & opts,                   &
          & traj0%profiles_coef,    &
          & traj0_sta%angles_coef,  &
          & coefs%coef,             &
          & traj0%aux_prof_coef,    &
          & traj0%predictors,       &
          & traj0%predictors%path1, &
          & traj0%raytracing_coef)! inout  (in because of mem allocation)
  ELSE IF (coefs%coef%fmv_model_ver == 8) THEN
    CALL rttov_setpredictors_8( &
          & opts,                   &
          & traj0%profiles_coef,    &
          & traj0_sta%angles_coef,  &
          & coefs%coef,             &
          & traj0%aux_prof_coef,    &
          & traj0%predictors,       &
          & traj0%predictors%path1, &
          & traj0%raytracing_coef)! inout  (in because of mem allocation)
  ELSE IF (coefs%coef%fmv_model_ver == 9) THEN
    IF (traj0_sta%dothermal) THEN
      CALL rttov_setpredictors_9(            &
            & opts,                          &
            & traj0%profiles_coef,           &
            & traj0_sta%angles_coef,         &
            & traj0%raytracing_coef%pathsat, &
            & coefs%coef_pccomp,             &
            & coefs%coef,                    &
            & traj0%predictors,              &
            & traj0%predictors%path1)
    ENDIF
    IF (traj0_sta%dosolar) THEN
      CALL rttov_setpredictors_9(            &
            & opts,                          &
            & traj0%profiles_coef,           &
            & traj0_sta%angles_coef,         &
            & traj0%raytracing_coef%patheff, &
            & coefs%coef_pccomp,             &
            & coefs%coef,                    &
            & traj0%predictors,              &
            & traj0%predictors%path2)
    ENDIF
  ELSE
    err = errorstatus_fatal
    THROWM(err .NE. 0, "Unexpected RTTOV compatibility version number")
  ENDIF

!  End Do ! Profile loop
!--------------------------------------------------------------------
!6. predict atmospheric and solar optical depths - atm/sol COEF levs
!--------------------------------------------------------------------

  IF (coefs%coef%fmv_model_ver == 9) THEN
    IF (traj0_sta%dothermal) THEN
      ! Calculate thermal path1 optical depths for all thermal channels
      CALL rttov_opdep_9( &
            & coefs%coef%nlayers,                  &
            & chanprof,                            &
            & traj0_sta%thermal,                   &
            & traj0%predictors,                    &
            & traj0%predictors%path1,              &
            & traj0%aux_prof_coef,                 &
            & coefs%coef,                          &
            & coefs%coef%thermal,                  &
            & traj0%opdp_path_coef%atm_level,      &
            & traj0_sta%thermal_path1%opdp_ref_coef)
    ENDIF
    IF (traj0_sta%dosolar) THEN
      ! Calculate solar path2 optical depths for all solar channels
      CALL rttov_opdep_9( &
            & coefs%coef%nlayers,                   &
            & chanprof,                             &
            & traj0_sta%solar,                      &
            & traj0%predictors,                     &
            & traj0%predictors%path2,               &
            & traj0%aux_prof_coef,                  &
            & coefs%coef,                           &
            & coefs%coef%solar,                     &
            & traj0%opdp_path_coef%sun_level_path2, &
            & traj0_sta%solar_path2%opdp_ref_coef)
    ENDIF
  ELSE
    CALL rttov_opdep( &
          & coefs%coef%nlayers,     &
          & chanprof,               &
          & traj0%predictors,       &
          & traj0%predictors%path1, &
          & traj0%aux_prof_coef,    &
          & coefs%coef,             &
          & coefs%coef%thermal,     &
          & traj0%opdp_path_coef,   &
          & traj0_sta%thermal_path1%opdp_ref_coef)
  ENDIF

!--------------------------------------------------------------------------
!7. interpolator second  call - optical depths from COEF levs to USER levs
!--------------------------------------------------------------------------

  IF (opts%interpolation%addinterp) THEN
    CALL rttov_intavg_chan( &
          & opts,                           &
          & traj0_sta%thermal,              &
          & traj0_sta%solar,                &
          & coefs%coef%nlevels,             &
          & nlevels,                        &
          & chanprof,                       &
          & traj0%profiles_coef,            &
          & profiles,                       &
          & traj0%opdp_path_coef,           &
          & traj0%opdp_path)
  ELSE
    CALL rttov_copy_opdp_path(opts, traj0%opdp_path, traj0%opdp_path_coef)
  ENDIF

  IF (opts%interpolation%spacetop) THEN
    traj0%opdp_path%atm_level(1,:) = 0._jprb
    IF (opts%rt_ir%addsolar) THEN
      traj0%opdp_path%sun_level_path2(1,:) = 0._jprb
    ENDIF
  ENDIF

  ! This is for coefficient testing; we replace the optical depths
  ! with those from the line-by-line convoluated with ISRF
  IF (PRESENT(lbl_check)) THEN
    IF (ASSOCIATED(lbl_check%atm_layer)) THEN
      traj0%opdp_path%atm_level(1,:) = 0._jprb
      traj0%opdp_path%atm_level(2:,:) = lbl_check%atm_layer(:,:)
      IF (traj0_sta%dosolar) THEN
        traj0%opdp_path%sun_level_path2(1,:) = 0._jprb
        traj0%opdp_path%sun_level_path2(2:,:) = lbl_check%atm_layer_path2(:,:)
      ENDIF
    ENDIF
  ENDIF


!--------------------------------------------------------------------
!8. If clouds are present,calculate the number of streams and
!   the cloud distribution  in each stream - remains on USER levs
!--------------------------------------------------------------------
  IF (PRESENT(traj_dyn)) THEN
    traj0_dyn => traj_dyn
  ELSE
    traj0_dyn => traj1_dyn
  ENDIF
  ltraj_dyn_dealloc = .NOT. PRESENT(traj_dyn)

  traj0_dyn%nstreams = 0_jpim

  IF (opts%rt_ir%addclouds) THEN
    CALL rttov_cldstr( &
          & opts%rt_ir,                  &
          & profiles,                    &
          & traj0%ircld,                 &
          & traj0_dyn%nstreams)
  ELSE
    traj0%ircld%xstrclr = 1._jprb
    traj0%ircld%nstream = 0_jpim
    traj0%ircld%icldarr = 0_jpim
  ENDIF

  CALL rttov_alloc_traj_dyn (err, traj0_dyn, opts, nchannels, profiles(1)%nlayers, &
                             traj0_dyn%nstreams, ncldtyp, 1_jpim, .TRUE._jplm)
  THROW(err .NE. 0)

!----------------------------------------------------------------------------
!9. Calculate optical depths of aerosols and/or clouds - USER levs
!----------------------------------------------------------------------------

  IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
    CALL rttov_opdpscattir( &
          & profiles(1)%nlayers,                    &
          & chanprof,                               &
          & opts,                                   &
          & traj0%aux_prof,                         &
          & profiles,                               &
          & traj0%profiles_dry,                     &
          & aer_opt_param,                          &
          & cld_opt_param,                          &
          & traj0_sta%dosolar,                      &
          & traj0_sta%solar,                        &
          & coefs%coef,                             &
          & coefs%coef_scatt_ir,                    &
          & traj0%raytracing,                       &
          & traj0%transmission_scatt_ir,            &
          & traj0_dyn%transmission_scatt_ir_stream, &
          & coefs%optp,                             &
          & traj0%ircld)
  ENDIF

!----------------------------------------
!10. calculate transmittances - USER levs
!----------------------------------------

  IF (traj0_sta%dothermal) THEN
    ! Thermal path1 transmittances for all thermal channels
    CALL rttov_transmit( &
          & opts%rt_ir%addaerosl,                     &
          & opts%rt_ir%addclouds,                     &
          & traj0_sta%do_lambertian,                  &
          & profiles(1)%nlayers,                      &
          & chanprof,                                 &
          & traj0_sta%thermal,                        &
          & traj0%aux_prof,                           &
          & coefs%coef,                               &
          & traj0%ircld,                              &
          & traj0_sta%angles,                         &
          & traj0%opdp_path%atm_level,                &
          & traj0_sta%thermal_path1%od_level,         &
          & transmission%tau_levels,                  &
          & transmission%tau_total,                   &
          & traj0_dyn%transmission_aux,               &
          & traj0_dyn%transmission_aux%thermal_path1, &
          & traj0%transmission_scatt_ir,              &
          & traj0_dyn%transmission_scatt_ir_stream,   &
          & traj0_sta%thermal_path1%tau_ref,          &
          & traj0_sta%thermal_path1%tau_ref_surf,     &
          & traj0_sta%thermal_path1%tau_surf,         &
          & traj0_sta%thermal_path1%tau_level)
  ENDIF

  IF (traj0_sta%dosolar) THEN
    ! Solar path2 transmittances for all solar channels
    CALL rttov_transmit_9_solar( &
          & opts%rt_ir%addaerosl,                   &
          & opts%rt_ir%addclouds,                   &
          & profiles(1)%nlayers,                    &
          & chanprof,                               &
          & traj0_sta%solar,                        &
          & traj0%aux_prof,                         &
          & coefs%coef,                             &
          & traj0%raytracing,                       &
          & traj0%ircld,                            &
          & traj0%opdp_path,                        &
          & traj0_sta%solar_path2,                  &
          & traj0_sta%solar_path1,                  &
          & transmission,                           &
          & traj0_dyn%transmission_aux,             &
          & traj0%transmission_scatt_ir,            &
          & traj0_dyn%transmission_scatt_ir_stream)

  ENDIF

!-------------------------------------------------
!11. calculate channel emissivity values - SURFACE
!-------------------------------------------------

  IF (traj0_sta%dothermal) THEN

    WHERE (.NOT. calcemis)
      emissivity%emis_out = emissivity%emis_in
    ELSEWHERE
      emissivity%emis_out = 0._jprb
    ENDWHERE

    traj0%thermrefl = 1._jprb - emissivity%emis_out

    IF (ANY(calcemis)) THEN
      ! Calculate surface emissivity and traj0%thermrefl for selected channels

      IF (coefs%coef%id_sensor == sensor_id_ir) THEN
        ! Infrared
        CALL rttov_calcemis_ir( &
              & profiles,             &
              & traj0_sta%angles,     &
              & coefs%coef,           &
              & opts%rt_ir%pc%addpc,  &
              & coefs%coef_pccomp,    &
              & traj0_sta%thermal,    &
              & chanprof,             &
              & calcemis,             &
              & emissivity%emis_out,  &
              & err)
        THROWM(err .NE. 0, "calcemis_ir")
        traj0%thermrefl(:) = 1._jprb - emissivity(:)%emis_out

      ELSE IF (coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po) THEN
        ! Microwave
        CALL rttov_calcemis_mw( &
              & opts,                          &
              & profiles,                      &
              & traj0_sta%angles,              &
              & coefs%coef,                    &
              & chanprof,                      &
              & traj0_dyn%transmission_aux,    &
              & calcemis,                      &
              & emissivity%emis_out,           &
              & traj0%thermrefl,            &
              & err)
        THROWM(err .NE. 0, "calcemis_mw")

      ELSE
        ! Hi-res - separate from IR to match TL/AD/K
        CALL rttov_calcemis_ir( &
              & profiles,             &
              & traj0_sta%angles,     &
              & coefs%coef,           &
              & opts%rt_ir%pc%addpc,  &
              & coefs%coef_pccomp,    &
              & traj0_sta%thermal,    &
              & chanprof,             &
              & calcemis,             &
              & emissivity%emis_out,  &
              & err)
        THROWM(err .NE. 0, "hi-res calcemis_ir")
        traj0%thermrefl(:) = 1._jprb - emissivity(:)%emis_out

      ENDIF

    ENDIF ! calcemis

  ELSE IF (PRESENT(emissivity)) THEN

    emissivity%emis_out = 0._jprb

  ENDIF ! dothermal


!-------------------------------------------------------
!12. Compute surface reflectance 
!-------------------------------------------------------

  IF (traj0_sta%dosolar) THEN
  
    WHERE (.NOT. calcrefl)
      reflectance%refl_out = reflectance%refl_in
    ELSEWHERE
      reflectance%refl_out = 0._jprb
    ENDWHERE

    IF (ANY(calcrefl)) THEN  

      CALL rttov_refsun( &
            & profiles,         &
            & coefs%coef,       &
            & traj0%aux_prof,   &
            & traj0%sunglint,   &
            & traj0%raytracing)

      CALL rttov_fresnel( &
            & chanprof,        &
            & profiles,        &
            & traj0_sta%solar, &
            & coefs%coef,      &
            & traj0%sunglint,  &
            & traj0%fresnrefl)

    ENDIF

    ! rttov_calcsurfrefl populates refl_norm(:) so is called regardless of calcrefl
    CALL rttov_calcsurfrefl(     &
          & coefs%coef,          &
          & profiles,            &
          & traj0%sunglint,      &
          & traj0%fresnrefl,     &
          & traj0_sta%solar,     &
          & chanprof,            &
          & traj0_sta%refl_norm, &
          & calcrefl,            &
          & emissivity,          &
          & reflectance%refl_out)
    THROWM(err .NE. 0, "calcsurfrefl")

  ELSE IF (PRESENT(reflectance)) THEN

    reflectance%refl_out = 0._jprb

  ENDIF

!---------------------------------------------------------
!13. integrate the radiative transfer equation - USER levs
!---------------------------------------------------------
! for new modeltop - '% layer' renamed '% air' (use of 'layer' now misleading)

  addcosmic = (coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po)
  CALL rttov_integrate( &
        & addcosmic,                              &
        & opts,                                   &
        & traj0_dyn%nstreams,                     &
        & chanprof,                               &
        & emissivity,                             &
        & reflectance,                            &
        & traj0_sta%refl_norm,                    &
        & traj0%thermrefl,                        &
        & traj0_sta%do_lambertian,                &
        & traj0_sta%thermal,                      &
        & traj0_sta%dothermal,                    &
        & traj0_sta%solar,                        &
        & traj0_sta%dosolar,                      &
        & traj0_sta%solar_spec_esd,               &
        & traj0_dyn%transmission_aux,             &
        & traj0%transmission_scatt_ir,            &
        & profiles,                               &
        & traj0%profiles_dry,                     &
        & traj0%aux_prof,                         &
        & coefs%coef,                             &
        & traj0%raytracing,                       &
        & traj0%ircld,                            &
        & radiancedata,                           &
        & radiancedata2,                          &
        & traj0_sta%auxrad,                       &
        & traj0_dyn%auxrad_stream)


!---------------------------------------------------------
! 14. Calculate output BTs/reflectances/PCscores
!---------------------------------------------------------

  IF (opts%rt_ir%pc%addpc) THEN
    CALL rttov_pcscores( &
          & opts,                        &
          & chanprof,                    &
          & traj0_sta%chanprof_pc,       &
          & pccomp,                      &
          & coefs%coef_pccomp,           &
          & radiancedata)

    IF (opts%rt_ir%pc%addradrec) THEN
      CALL rttov_reconstruct( &
            & opts,                        &
            & traj0_sta%chanprof_in,       &
            & traj0_sta%chanprof_pc,       &
            & pccomp,                      &
            & coefs%coef_pccomp)
    ENDIF
  ENDIF

  ! Do NLTE bias correction for hyperspectral instruments (not PC-RTTOV)
  IF (.NOT. opts%rt_ir%pc%addpc) THEN
    IF (coefs%coef%nltecoef) THEN
      IF (opts%rt_ir%do_nlte_correction .AND. &
          coefs%coef%id_sensor == sensor_id_hi) THEN
        CALL rttov_nlte_bias_correction(coefs%coef, profiles, traj0_sta%angles, chanprof, radiancedata)
      ENDIF
    ENDIF
  ENDIF

  IF (opts%rt_ir%pc%addpc) THEN

    IF (opts%rt_ir%pc%addradrec) &
      CALL rttov_calcbt_pc(traj0_sta%chanprof_in, coefs%coef_pccomp, pccomp)

  ELSE

    IF (traj0_sta%dothermal) THEN
      CALL rttov_calcbt(chanprof, coefs%coef, traj0_sta%thermal, radiancedata)
    ELSE
      radiancedata%bt_clear(:) = 0._jprb
      radiancedata%bt(:) = 0._jprb
    ENDIF
    IF (traj0_sta%dosolar) THEN
      CALL rttov_calcsatrefl(chanprof, profiles, traj0_sta%solar_spec_esd, traj0_sta%solar, radiancedata)
    ELSE
      radiancedata%refl_clear(:) = 0._jprb
      radiancedata%refl(:) = 0._jprb
    ENDIF
  ENDIF

!---------------------
!15. deallocate memory
!---------------------

  CALL cleanup()

  IF (LHOOK) CALL DR_HOOK('RTTOV_DIRECT', 1_jpim, ZHOOK_HANDLE)
  CATCH_C
  errorstatus = err
  IF (LHOOK) CALL DR_HOOK('RTTOV_DIRECT', 1_jpim, ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE cleanup()
    INTEGER(KIND=jpim) :: error

    IF (ltraj_dyn_dealloc) THEN
      IF (ASSOCIATED(traj0_dyn)) THEN
        CALL rttov_alloc_traj_dyn (error, traj0_dyn, opts, nchannels, profiles(1)%nlayers, &
                                   traj0_dyn%nstreams, ncldtyp, 0_jpim, .TRUE._jplm)
      ENDIF
    ENDIF

    IF (ltraj_sta_dealloc) THEN
      IF (ASSOCIATED(traj0_sta)) THEN
        CALL rttov_alloc_traj_sta (error, traj0_sta, opts, coefs%coef, nlevels, nchannels, nprofiles, &
                                   0_jpim, npcscores, channels_rec)
      ENDIF
    ENDIF

    IF (ASSOCIATED(traj0)) THEN
      CALL rttov_check_traj( &
            & error,         &
            & nprofiles,     &
            & nchannels,     &
            & opts,          &
            & nlevels,       &
            & coefs,         &
            & 0_jpim,        &
            & traj0 = traj0, &
            & traj1 = traj1, &
            & traj2 = traj)
    ENDIF
  END SUBROUTINE cleanup

END SUBROUTINE rttov_direct
