!
SUBROUTINE rttov_k_ad( &
            & errorstatus,      &
            & chanprof,         &
            & opts,             &
            & profiles,         &
            & profiles_k,       &
            & coefs,            &
            & transmission,     &
            & radiancedata,     &
            & calcemis,         &
            & emissivity,       &
            & emissivity_k,     &
            & calcrefl,         &
            & reflectance,      &
            & reflectance_k,    &
            & aer_opt_param,    &
            & cld_opt_param,    &
            & pccomp,           &
            & profiles_k_pc,    &
            & profiles_k_rec,   &
            & channels_rec,     &
            & nthreads)
! Description:
!   Compute jacobian using adjoint code
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
#include "throw.h"
!INTF_ON

! Imported Type Definitions:
  USE rttov_types, ONLY :  &
       & rttov_coefs,       &
       & rttov_options,     &
       & profile_Type,      &
       & transmission_Type, &
       & rttov_chanprof,    &
       & rttov_emissivity,  &
       & rttov_reflectance, &
       & rttov_opt_param,   &
       & radiance_Type,     &
       & rttov_pccomp
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE parkind1, ONLY : jprb
  USE rttov_types, ONLY :  &
       & rttov_traj
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  INTEGER(KIND=jpim)     , INTENT(OUT)           :: errorstatus
  TYPE(profile_Type  )   , INTENT(IN)   , TARGET :: profiles(:)
  TYPE(rttov_chanprof)   , INTENT(IN)            :: chanprof(:)
  TYPE(profile_Type  )   , INTENT(INOUT), TARGET :: profiles_k(size(chanprof))
! for calls from RTTOV_CLD_K and RTTOV_SCATT_K.
  TYPE(rttov_options )   , INTENT(IN)            :: opts
  TYPE(rttov_coefs   )   , INTENT(IN)   , TARGET :: coefs
  TYPE(transmission_Type), INTENT(INOUT)         :: transmission                    ! in because of mem allocation
  TYPE(radiance_Type    ), INTENT(INOUT)         :: radiancedata                    ! in because of mem allocation
  LOGICAL(KIND=jplm)     , INTENT(IN)   , OPTIONAL         :: calcemis(size(chanprof))
  TYPE(rttov_emissivity) , INTENT(INOUT), OPTIONAL         :: emissivity(size(chanprof))
  TYPE(rttov_emissivity) , INTENT(INOUT), OPTIONAL         :: emissivity_k(size(chanprof))
  LOGICAL(KIND=jplm)     , INTENT(IN)   , OPTIONAL         :: calcrefl(size(chanprof))
  TYPE(rttov_reflectance), INTENT(INOUT), OPTIONAL         :: reflectance(size(chanprof))
  TYPE(rttov_reflectance), INTENT(INOUT), OPTIONAL         :: reflectance_k(size(chanprof))
  TYPE(rttov_opt_param  ), INTENT(IN)   , OPTIONAL         :: aer_opt_param
  TYPE(rttov_opt_param  ), INTENT(IN)   , OPTIONAL         :: cld_opt_param
  TYPE(rttov_pccomp     ), INTENT(INOUT), OPTIONAL         :: pccomp
  TYPE(profile_Type     ), INTENT(INOUT), OPTIONAL, TARGET :: profiles_k_pc(:)
  TYPE(profile_Type     ), INTENT(INOUT), OPTIONAL, TARGET :: profiles_k_rec(:)
  INTEGER(KIND=jpim)     , INTENT(IN)   , OPTIONAL         :: channels_rec(:)
  INTEGER(KIND=jpim)     , INTENT(IN)                      :: nthreads
!INTF_END

#include "rttov_ad.interface"
#include "rttov_parallel_ad.interface"
#include "rttov_alloc_prof.interface"
#include "rttov_alloc_rad.interface"
#include "rttov_alloc_pccomp.interface"
#include "rttov_alloc_traj.interface"
#include "rttov_init_prof.interface"
#include "rttov_copy_prof.interface"
#include "rttov_alloc_transmission.interface"
#include "rttov_init_pccomp.interface"
#include "rttov_init_rad.interface"
#include "rttov_init_transmission.interface"
#include "rttov_errorreport.interface"

  TYPE(rttov_traj)        :: traj, traj_ad
  TYPE(rttov_pccomp)      :: pccomp_ad
  TYPE(profile_Type), TARGET :: profiles_ad(size(profiles))
  TYPE(rttov_emissivity)     :: emissivity_ad(size(chanprof))
  TYPE(rttov_reflectance)    :: reflectance_ad(size(chanprof))
  TYPE(transmission_Type) :: transmission_ad
  TYPE(radiance_Type    ) :: radiancedata_ad
  INTEGER(KIND=jpim)      :: ichan          , err
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: npcscores
  INTEGER(KIND=jpim) :: nchannels_rec
  INTEGER(KIND=jpim) :: ipcscore
  INTEGER(KIND=jpim) :: ichannels_rec
  INTEGER(KIND=jpim) :: iprof

TRY


  nprofiles = size(profiles)
  nchannels = size(chanprof)

  npcscores = -1_jpim
  nchannels_rec = -1_jpim

  IF (opts%rt_ir%pc%addpc) THEN
    npcscores = size(pccomp%pcscores)
    IF (opts%rt_ir%pc%addradrec) THEN
      nchannels_rec = size(pccomp%bt_pccomp)
    ENDIF
  ENDIF

  IF (opts%rt_ir%pc%addpc) THEN
    CALL rttov_alloc_pccomp (err, pccomp_ad, npcscores, 1_jpim, nchannels_rec = nchannels_rec)
    THROW(err.ne.0)
  ENDIF

  CALL rttov_alloc_prof( &
        & err,                                                           &
        & nprofiles,                                                     &
        & profiles_ad,                                                   &
        & profiles(1)%nlevels,                                           &
        & opts,                                                          &
        & 1_jpim,                                                        &
        & coefs = coefs,                                                 &
        & init = .TRUE._jplm)
  THROW(err.ne.0)

  CALL rttov_alloc_rad( &
        & err,                          &
        & nchannels,                    &
        & radiancedata_ad,              &
        & profiles(1)%nlevels - 1_jpim, &
        & 1_jpim)
  THROW(err.ne.0)

  CALL rttov_alloc_transmission(err, transmission_ad, profiles(1)%nlevels-1_jpim, nchannels, 1_jpim)
  THROW(err.ne.0)

  IF (nthreads <= 0) THEN
    CALL rttov_alloc_traj(                          &
        & err,                  nprofiles,          &
        & nchannels,            opts,               &
        & profiles(1)%nlevels,  coefs,              &
        & 1_jpim,               traj = traj,        &
        & traj_ad = traj_ad )
    THROW(err.ne.0)
  ENDIF

!
! Regular channels
!
  DO ichan = 1, nchannels

    emissivity_ad%emis_in             = 0._jprb
    emissivity_ad%emis_out            = 0._jprb
    reflectance_ad%refl_in            = 0._jprb
    reflectance_ad%refl_out           = 0._jprb

    CALL rttov_init_transmission (transmission_ad)
    CALL rttov_init_rad (radiancedata_ad)
    IF (opts%rt_ir%pc%addpc) THEN
      radiancedata_ad%total(ichan) = 1._jprb
    ELSE IF (opts%rt_all%switchrad .and. &
             coefs%coef%ss_val_chn(chanprof(ichan)%chan) < 2) THEN
      radiancedata_ad%bt(ichan) = 1._jprb
    ELSE
      radiancedata_ad%total(ichan) = 1._jprb
    ENDIF

    CALL rttov_init_prof (profiles_ad)

    IF (opts%rt_ir%pc%addpc) THEN
      CALL rttov_init_pccomp (pccomp_ad)
    ENDIF

    IF (nthreads <= 0) THEN
      CALL rttov_ad(                     &
            & errorstatus,               &
            & chanprof,                  &
            & opts,                      &
            & profiles,                  &
            & profiles_ad,               &
            & coefs,                     &
            & transmission,              &
            & transmission_ad,           &
            & radiancedata,              &
            & radiancedata_ad,           &
            & calcemis,                  &
            & emissivity,                &
            & emissivity_ad,             &
            & calcrefl,                  &
            & reflectance,               &
            & reflectance_ad,            &
            & aer_opt_param,             &
            & cld_opt_param,             &
            & traj         = traj,       &
            & traj_ad      = traj_ad,    &
            & pccomp       = pccomp,     &
            & pccomp_ad    = pccomp_ad,  &
            & channels_rec = channels_rec)
    ELSE
      CALL rttov_parallel_ad(            &
            & errorstatus,               &
            & chanprof,                  &
            & opts,                      &
            & profiles,                  &
            & profiles_ad,               &
            & coefs,                     &
            & transmission,              &
            & transmission_ad,           &
            & radiancedata,              &
            & radiancedata_ad,           &
            & calcemis,                  &
            & emissivity,                &
            & emissivity_ad,             &
            & calcrefl,                  &
            & reflectance,               &
            & reflectance_ad,            &
            & aer_opt_param,             &
            & cld_opt_param,             &
            & pccomp       = pccomp,     &
            & pccomp_ad    = pccomp_ad,  &
            & channels_rec = channels_rec, &
            & nthreads     = nthreads)
    ENDIF

    IF (errorstatus /= 0) THEN
      err = errorstatus_fatal
      THROW(err.ne.0)
    ENDIF

    iprof = chanprof(ichan)%prof
    CALL rttov_copy_prof (profiles_k(ichan:ichan), profiles_ad(iprof:iprof))
    emissivity_k(ichan)%emis_in = emissivity_ad(ichan)%emis_in
    reflectance_k(ichan)%refl_in = reflectance_ad(ichan)%refl_in

  ENDDO

  IF (opts%rt_ir%pc%addpc) THEN

    IF (opts%rt_ir%pc%addradrec) THEN

!
! Reconstructed radiance
!

      DO ichannels_rec = 1, nchannels_rec
      
        emissivity_ad%emis_in             = 0._jprb
        emissivity_ad%emis_out            = 0._jprb
        reflectance_ad%refl_in            = 0._jprb
        reflectance_ad%refl_out           = 0._jprb

        CALL rttov_init_transmission (transmission_ad)
        CALL rttov_init_rad (radiancedata_ad)
      
        CALL rttov_init_pccomp (pccomp_ad)

        IF (opts%rt_all%switchrad) THEN
          pccomp_ad%bt_pccomp(ichannels_rec) = 1._jprb
        ELSE
          pccomp_ad%total_pccomp(ichannels_rec) = 1._jprb
        ENDIF
      
        CALL rttov_init_prof (profiles_ad)

        IF (nthreads <= 0) THEN
          CALL rttov_ad(                     &
                & errorstatus,               &
                & chanprof,                  &
                & opts,                      &
                & profiles,                  &
                & profiles_ad,               &
                & coefs,                     &
                & transmission,              &
                & transmission_ad,           &
                & radiancedata,              &
                & radiancedata_ad,           &
                & calcemis,                  &
                & emissivity,                &
                & emissivity_ad,             &
                & calcrefl,                  &
                & reflectance,               &
                & reflectance_ad,            &
                & aer_opt_param,             &
                & cld_opt_param,             &
                & traj         = traj,       &
                & traj_ad      = traj_ad,    &
                & pccomp       = pccomp,     &
                & pccomp_ad    = pccomp_ad,  &
                & channels_rec = channels_rec)
        ELSE
          CALL rttov_parallel_ad(            &
                & errorstatus,               &
                & chanprof,                  &
                & opts,                      &
                & profiles,                  &
                & profiles_ad,               &
                & coefs,                     &
                & transmission,              &
                & transmission_ad,           &
                & radiancedata,              &
                & radiancedata_ad,           &
                & calcemis,                  &
                & emissivity,                &
                & emissivity_ad,             &
                & calcrefl,                  &
                & reflectance,               &
                & reflectance_ad,            &
                & aer_opt_param,             &
                & cld_opt_param,             &
                & pccomp       = pccomp,     &
                & pccomp_ad    = pccomp_ad,  &
                & channels_rec = channels_rec, &
                & nthreads     = nthreads)
        ENDIF

        IF (errorstatus /= 0) THEN
          err = errorstatus_fatal
          THROW(err.ne.0)
        ENDIF
      
        iprof = 1 + ((nprofiles * (ichannels_rec-1))/nchannels_rec)
        CALL rttov_copy_prof (profiles_k_rec(ichannels_rec:ichannels_rec), profiles_ad(iprof:iprof))
      
      ENDDO

    ELSE
!
! PC-scores
!
      DO ipcscore = 1, npcscores
     
        emissivity_ad%emis_in             = 0._jprb
        emissivity_ad%emis_out            = 0._jprb
        reflectance_ad%refl_in            = 0._jprb
        reflectance_ad%refl_out           = 0._jprb

        CALL rttov_init_transmission (transmission_ad)
        CALL rttov_init_rad (radiancedata_ad)
      
        CALL rttov_init_pccomp (pccomp_ad)

        pccomp_ad%pcscores(ipcscore) = 1._jprb
      
        CALL rttov_init_prof (profiles_ad)

        IF (nthreads <= 0) THEN
          CALL rttov_ad(                     &
                & errorstatus,               &
                & chanprof,                  &
                & opts,                      &
                & profiles,                  &
                & profiles_ad,               &
                & coefs,                     &
                & transmission,              &
                & transmission_ad,           &
                & radiancedata,              &
                & radiancedata_ad,           &
                & calcemis,                  &
                & emissivity,                &
                & emissivity_ad,             &
                & calcrefl,                  &
                & reflectance,               &
                & reflectance_ad,            &
                & aer_opt_param,             &
                & cld_opt_param,             &
                & traj         = traj,       &
                & traj_ad      = traj_ad,    &
                & pccomp       = pccomp,     &
                & pccomp_ad    = pccomp_ad,  &
                & channels_rec = channels_rec)
        ELSE
          CALL rttov_parallel_ad(            &
                & errorstatus,               &
                & chanprof,                  &
                & opts,                      &
                & profiles,                  &
                & profiles_ad,               &
                & coefs,                     &
                & transmission,              &
                & transmission_ad,           &
                & radiancedata,              &
                & radiancedata_ad,           &
                & calcemis,                  &
                & emissivity,                &
                & emissivity_ad,             &
                & calcrefl,                  &
                & reflectance,               &
                & reflectance_ad,            &
                & aer_opt_param,             &
                & cld_opt_param,             &
                & pccomp       = pccomp,     &
                & pccomp_ad    = pccomp_ad,  &
                & channels_rec = channels_rec, &
                & nthreads = nthreads)
        ENDIF

        IF (errorstatus /= 0) THEN
          err = errorstatus_fatal
          THROW(err.ne.0)
        ENDIF
     
        iprof = 1 + ((nprofiles * (ipcscore-1))/npcscores)
        CALL rttov_copy_prof (profiles_k_pc(ipcscore:ipcscore), profiles_ad(iprof:iprof))
     
      ENDDO

    ENDIF

  ENDIF

  IF (opts%rt_ir%pc%addpc) THEN
    CALL rttov_alloc_pccomp (err, pccomp_ad, npcscores, 0_jpim, nchannels_rec = nchannels_rec)
    THROW(err.ne.0)
  ENDIF

  IF (nthreads <= 0) THEN
    CALL rttov_alloc_traj(                          &
        & err,                  nprofiles,          &
        & nchannels,            opts,               &
        & profiles(1)%nlevels,  coefs,              &
        & 0_jpim,               traj = traj,        &
        & traj_ad = traj_ad )
    THROW(err.ne.0)
  ENDIF

  CALL rttov_alloc_transmission(err, transmission_ad, profiles(1)%nlevels-1_jpim, nchannels, 0_jpim)
  THROW(err.ne.0)

  CALL rttov_alloc_rad( &
        & err,                          &
        & nchannels,                    &
        & radiancedata_ad,              &
        & profiles(1)%nlevels - 1_jpim, &
        & 0_jpim)
  THROW(err.ne.0)

  CALL rttov_alloc_prof( &
        & err,                 &
        & nprofiles,           &
        & profiles_ad,         &
        & profiles(1)%nlevels, &
        & opts,                &
        & 0_jpim)
  THROW(err.ne.0)

CATCH

  errorstatus = err

END SUBROUTINE rttov_k_ad

