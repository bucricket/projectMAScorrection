! Description:
!> @file
!!   Allocate/deallocate any/all arrays and structures for RTTOV tangent linear model.
!
!> @brief
!!   Allocate/deallocate any/all arrays and structures for RTTOV tangent linear model.
!!
!! @details
!!   Argument order follows that of rttov_tl where possible.
!!   This allows all necessary input/output arguments to rttov_tl to be
!!   allocated or deallocated in a single subroutine call: pass only those
!!   arguments you wish to allocate/deallocate.
!!
!! @param[out]    err             status on exit
!! @param[in]     asw             1_jpim => allocate; 0_jpim => deallocate
!! @param[in]     nprofiles       number of profiles being simulated
!! @param[in]     nchanprof       total number of channels being simulated over all profiles
!! @param[in]     nlevels         number of levels in input profiles
!! @param[in]     chanprof        specifies channels and profiles to simulate, optional
!! @param[in]     opts            options to configure the simulations
!! @param[in,out] profiles        input atmospheric profiles and surface variables, optional
!! @param[in,out] profiles_tl     input atmospheric profile and surface variable perturbations, optional
!! @param[in]     coefs           coefficients structure for instrument to simulate
!! @param[in,out] transmission    output transmittances, optional
!! @param[in,out] transmission_tl output transmittance perturbations, optional
!! @param[in,out] radiancedata    output radiances and corresponding BTs and BRFs, optional
!! @param[in,out] radiancedata_tl output radiance, BT and BRF perturbations, optional
!! @param[in]     calcemis        flags for internal RTTOV surface emissivity calculation, optional
!! @param[in,out] emissivity      input/output surface emissivities, optional
!! @param[in,out] emissivity_tl   input/output surface emissivity perturbations, optional
!! @param[in]     calcrefl        flags for internal RTTOV surface BRDF calculation, optional
!! @param[in,out] reflectance     input/output surface BRDFs, input cloud top BRDF for simple cloud, optional
!! @param[in,out] reflectance_tl  input/output surface BRDF perturbations, optional
!! @param[in]     aer_nphangle    number of angles for which phase fns are defined in aer_opt_param, optional
!! @param[in]     aer_opt_param   input aerosol optical parameters, optional
!! @param[in]     cld_nphangle    number of angles for which phase fns are defined in cld_opt_param, optional
!! @param[in]     cld_opt_param   input cloud optical parameters, optional
!! @param[in,out] traj            RTTOV direct internal state, can be initialised outside RTTOV, optional
!! @param[in,out] traj_tl         RTTOV TL internal state, can be initialised outside RTTOV, optional
!! @param[in]     npcscores       number of PC scores to calculate for PC-RTTOV, optional
!! @param[in]     nchannels_rec   number of channels for which to calculate reconstructed radiances, optional
!! @param[in,out] pccomp          output PC scores and radiances from PC-RTTOV, optional
!! @param[in,out] pccomp_tl       output PC score and radiance perturbations from PC-RTTOV, optional
!! @param[in,out] channels_rec    array for channels for which to calculate reconstructed radiances, optional
!! @param[in]     init            set .TRUE. to initialise newly allocated structures, optional
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
SUBROUTINE rttov_alloc_tl( &
              err,              &
              asw,              &
              nprofiles,        &
              nchanprof,        &
              nlevels,          &
              chanprof,         &
              opts,             &
              profiles,         &
              profiles_tl,      &
              coefs,            &
              transmission,     &
              transmission_tl,  &
              radiancedata,     &
              radiancedata_tl,  &
              calcemis,         &
              emissivity,       &
              emissivity_tl,    &
              calcrefl,         &
              reflectance,      &
              reflectance_tl,   &
              aer_nphangle,     &
              aer_opt_param,    &
              cld_nphangle,     &
              cld_opt_param,    &
              traj,             &
              traj_tl,          &
              npcscores,        &
              nchannels_rec,    &
              pccomp,           &
              pccomp_tl,        &
              channels_rec,     &
              init)

!INTF_OFF
#include "throw.h"
!INTF_ON

  USE parkind1, ONLY : jpim, jplm

  USE rttov_types, ONLY : &
      rttov_options,      &
      rttov_coefs,        &
      rttov_chanprof,     &
      rttov_emissivity,   &
      rttov_reflectance,  &
      profile_type,       &
      transmission_type,  &
      radiance_type,      &
      rttov_traj,         &
      rttov_pccomp,       &
      rttov_opt_param

  IMPLICIT NONE

  ! Input arguments
  INTEGER(jpim),           INTENT(OUT)             :: err               ! return code
  INTEGER(jpim),           INTENT(IN)              :: asw               ! 1=allocate, 0=deallocate
  INTEGER(jpim),           INTENT(IN)              :: nprofiles         ! number of profiles
  INTEGER(jpim),           INTENT(IN)              :: nchanprof         ! size of chanprof array
  INTEGER(jpim),           INTENT(IN)              :: nlevels           ! number of levels
  TYPE(rttov_options),     INTENT(IN)              :: opts              ! RTTOV options
  TYPE(rttov_coefs),       INTENT(IN),    TARGET   :: coefs             ! coefficients

  ! Output arguments
  TYPE(profile_type),      POINTER,       OPTIONAL :: profiles(:)       ! input profiles
  TYPE(profile_type),      POINTER,       OPTIONAL :: profiles_tl(:)    ! input profile perturbation
  TYPE(transmission_type), INTENT(INOUT), OPTIONAL :: transmission      ! output transmission
  TYPE(transmission_type), INTENT(INOUT), OPTIONAL :: transmission_tl   ! output transmission TL
  TYPE(radiance_type),     INTENT(INOUT), OPTIONAL :: radiancedata      ! output radiances
  TYPE(radiance_type),     INTENT(INOUT), OPTIONAL :: radiancedata_tl   ! output radiances TL

  TYPE(rttov_chanprof),    POINTER,       OPTIONAL :: chanprof(:)       ! channel and profile indices
  LOGICAL(jplm),           POINTER,       OPTIONAL :: calcemis(:)       ! switch for emissivity calculations
  TYPE(rttov_emissivity),  POINTER,       OPTIONAL :: emissivity(:)     ! input/output emissivity values
  TYPE(rttov_emissivity),  POINTER,       OPTIONAL :: emissivity_tl(:)  ! input/output emissivity perturbations
  LOGICAL(jplm),           POINTER,       OPTIONAL :: calcrefl(:)       ! switch for surface BRDF calculations
  TYPE(rttov_reflectance), POINTER,       OPTIONAL :: reflectance(:)    ! input/output surface BRDF values
  TYPE(rttov_reflectance), POINTER,       OPTIONAL :: reflectance_tl(:) ! input/output surface BRDF perturbations

  TYPE(rttov_traj),        INTENT(INOUT), OPTIONAL :: traj              ! RTTOV internal data
  TYPE(rttov_traj),        INTENT(INOUT), OPTIONAL :: traj_tl           ! RTTOV internal data TL

  ! Aerosol/cloud optical parameter arguments
  INTEGER(jpim),           INTENT(IN),    OPTIONAL :: aer_nphangle      ! number of phase angles for aerosols
  TYPE(rttov_opt_param),   INTENT(INOUT), OPTIONAL :: aer_opt_param     ! input aerosol optical parameters
  INTEGER(jpim),           INTENT(IN),    OPTIONAL :: cld_nphangle      ! number of phase angles for clouds
  TYPE(rttov_opt_param),   INTENT(INOUT), OPTIONAL :: cld_opt_param     ! input cloud optical parameters

  ! PC-RTTOV output arguments
  INTEGER(jpim),           INTENT(IN),    OPTIONAL :: npcscores         ! number of PC scores
  INTEGER(jpim),           INTENT(IN),    OPTIONAL :: nchannels_rec     ! number of reconstructed radiances
  TYPE(rttov_pccomp),      INTENT(INOUT), OPTIONAL :: pccomp            ! output PC scores
  TYPE(rttov_pccomp),      INTENT(INOUT), OPTIONAL :: pccomp_tl         ! output PC scores TL
  INTEGER(jpim),           POINTER,       OPTIONAL :: channels_rec(:)   ! reconstructed channel list

  ! Initialisation flag
  LOGICAL(jplm),           INTENT(IN),    OPTIONAL :: init              ! flag to select initialisation of structures
!INTF_END

#include "rttov_errorreport.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_alloc_traj.interface"

!- End of header --------------------------------------------------------

  TRY

  ! Check TL arguments before allocating anything
  IF (PRESENT(pccomp_tl) .AND. .NOT. PRESENT(npcscores)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'npscores required when (de)allocating pccomp_tl')
  ENDIF

  CALL rttov_alloc_direct( &
          err,                         &
          asw,                         &
          nprofiles,                   &
          nchanprof,                   &
          nlevels,                     &
          chanprof,                    &
          opts,                        &
          profiles,                    &
          coefs,                       &
          transmission,                &
          radiancedata,                &
          calcemis=calcemis,           &
          emissivity=emissivity,       &
          calcrefl=calcrefl,           &
          reflectance=reflectance,     &
          aer_nphangle=aer_nphangle,   &
          aer_opt_param=aer_opt_param, &
          cld_nphangle=cld_nphangle,   &
          cld_opt_param=cld_opt_param, &
          traj=traj,                   &
          npcscores=npcscores,         &
          nchannels_rec=nchannels_rec, &
          pccomp=pccomp,               &
          channels_rec=channels_rec,   &
          init=init)
  THROWM(err.NE.0, 'Error allocating direct structures')

  ! Use rttov_alloc_direct to allocate most of the TL variables
  CALL rttov_alloc_direct( &
          err,                          &
          asw,                          &
          nprofiles,                    &
          nchanprof,                    &
          nlevels,                      &
          opts=opts,                    &
          profiles=profiles_tl,         &
          coefs=coefs,                  &
          transmission=transmission_tl, &
          radiancedata=radiancedata_tl, &
          emissivity=emissivity_tl,     &
          reflectance=reflectance_tl,   &
          npcscores=npcscores,          &
          nchannels_rec=nchannels_rec,  &
          pccomp=pccomp_tl,             &
          init=init)
  THROWM(err.NE.0, 'Error allocating TL structures')

  ! traj_tl is (de)allocated here
  IF (PRESENT(traj_tl)) THEN
    CALL rttov_alloc_traj( &
            err,           &
            nprofiles,     &
            nchanprof,     &
            opts,          &
            nlevels,       &
            coefs,         &
            asw,           &
            traj_tl=traj_tl)
    THROWM(err.NE.0, 'Error (de)allocating traj_tl structure')
  ENDIF

  CATCH
END SUBROUTINE rttov_alloc_tl