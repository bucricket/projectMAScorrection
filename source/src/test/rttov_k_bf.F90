!
SUBROUTINE rttov_k_bf( &
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
!   Compute jacobian using brute force
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
  USE rttov_chain, ONLY :  &
       & chain,              &
       & pchain,             &
       & size_chain,         &
       & get_pointer_chain,  &
       & advance_chain
  USE rttov_test_k_mod, ONLY : &
    & make_chain_profile, &
    & free_chain_profile, &
    & assign_chain_profile
       
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  INTEGER(KIND=jpim)     , INTENT(OUT)           :: errorstatus
  TYPE(rttov_chanprof)   , INTENT(IN)            :: chanprof(:)
  TYPE(profile_Type  )   , INTENT(IN)   , TARGET :: profiles(:)
  TYPE(profile_Type  )   , INTENT(INOUT), TARGET :: profiles_k(size(chanprof))
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

#include "rttov_direct.interface"
#include "rttov_parallel_direct.interface"
#include "rttov_copy_prof.interface"
#include "rttov_alloc_prof.interface"
#include "rttov_alloc_pccomp.interface"
#include "rttov_alloc_traj.interface"
#include "rttov_alloc_rad.interface"
#include "rttov_alloc_transmission.interface"
#include "rttov_make_profile_inc.interface"
#include "rttov_errorreport.interface"
#include "rttov_alloc_pc_dimensions.interface"

  TYPE(rttov_traj)       :: traj
  TYPE(rttov_pccomp)     :: pccomp_bf
  TYPE(rttov_chanprof), POINTER :: chanprof_pc(:)
  TYPE(rttov_chanprof), POINTER :: chanprof_in(:)
  REAL(KIND=jprb),      POINTER :: w(:)
  REAL(KIND=jprb),      POINTER :: w_pc(:)
  REAL(KIND=jprb),      POINTER :: w_in(:)
  TYPE(profile_Type), TARGET :: profiles_bf (size(profiles))
  TYPE(profile_Type), TARGET :: profiles_inc(size(profiles))
  LOGICAL(KIND=jplm)      :: calcemis_bf      (size(chanprof))
  TYPE(rttov_emissivity)  :: emissivity_bf    (size(chanprof))
  LOGICAL(KIND=jplm)      :: calcrefl_bf      (size(chanprof))
  TYPE(rttov_reflectance) :: reflectance_bf   (size(chanprof))
  TYPE(transmission_Type) :: transmission_bf
  TYPE(radiance_Type    ) :: radiancedata_bf
  REAL(KIND=jprb)         :: out_bf(size(chanprof))
  TYPE(chain),  POINTER :: chain_profiles(:)
  TYPE(chain),  POINTER :: chain_profiles_bf(:)
  TYPE(chain),  POINTER :: chain_profiles_k(:)
  TYPE(chain),  POINTER :: chain_profiles_k_pc(:)
  TYPE(chain),  POINTER :: chain_profiles_k_rec(:)
  TYPE(chain),  POINTER :: chain_profiles_inc(:)
  TYPE(pchain), POINTER :: c(:)
  TYPE(pchain), POINTER :: c_bf(:)
  TYPE(pchain), POINTER :: c_k(:)
  TYPE(pchain), POINTER :: c_k_pc(:)
  TYPE(pchain), POINTER :: c_k_in(:)
  TYPE(pchain), POINTER :: c_inc(:)
  INTEGER(KIND=jpim) :: iprof, ichan, i
  REAL(KIND=jprb), POINTER :: x_profiles, x_profiles_inc, x_profiles_bf
  INTEGER(KIND=jpim) :: err, nsize
  REAL   (KIND=jprb) :: reflout(size(chanprof))
  REAL   (KIND=jprb) :: inc(size(profiles))
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: npcscores
  INTEGER(KIND=jpim) :: nchannels_rec
  
TRY




  nprofiles = size(profiles)
  nchannels = size(chanprof)

  CALL rttov_alloc_prof( &
        & err,                                                           &
        & nprofiles,                                                     &
        & profiles_bf,                                                   &
        & profiles(1)%nlevels,                                           &
        & opts,                                                          &
        & 1_jpim,                                                        &
        & coefs = coefs,                                                 &
        & init = .TRUE._jplm)
  THROW(err.ne.0)

  CALL rttov_copy_prof(profiles_bf, profiles)

  CALL rttov_alloc_prof( &
        & err,                                                           &
        & nprofiles,                                                     &
        & profiles_inc,                                                  &
        & profiles(1)%nlevels,                                           &
        & opts,                                                          &
        & 1_jpim,                                                        &
        & coefs = coefs,                                                 &
        & init = .TRUE._jplm)
  THROW(err.ne.0)

  CALL rttov_alloc_rad( &
        & err,                          &
        & nchannels,                    &
        & radiancedata_bf,              &
        & profiles(1)%nlevels - 1_jpim, &
        & 1_jpim)
  THROW(err.ne.0)

  CALL rttov_alloc_transmission(err, transmission_bf, profiles(1)%nlevels-1_jpim, nchannels, 1_jpim)
  THROW(err.ne.0)

  IF (nthreads <= 0) THEN
    CALL rttov_alloc_traj(                          &
        & err,                  nprofiles,          &
        & nchannels,            opts,               &
        & profiles(1)%nlevels,  coefs,              &
        & 1_jpim,               traj = traj)
    THROW(err.ne.0)
  ENDIF

  CALL rttov_make_profile_inc(profiles_inc, profiles, opts)

  IF (nthreads <= 0) THEN
    CALL rttov_direct(                     &
          & errorstatus,                   &
          & chanprof,                      &
          & opts,                          &
          & profiles,                      &
          & coefs,                         &
          & transmission,                  &
          & radiancedata,                  &
          & calcemis = calcemis,           &
          & emissivity = emissivity,       &
          & calcrefl = calcrefl,           &
          & reflectance = reflectance,     &
          & aer_opt_param = aer_opt_param, &
          & cld_opt_param = cld_opt_param, &
          & traj         = traj,           &
          & pccomp       = pccomp,         &
          & channels_rec = channels_rec)
  ELSE
    CALL rttov_parallel_direct(            &
          & errorstatus,                   &
          & chanprof,                      &
          & opts,                          &
          & profiles,                      &
          & coefs,                         &
          & transmission,                  &
          & radiancedata,                  &
          & calcemis = calcemis,           &
          & emissivity = emissivity,       &
          & calcrefl = calcrefl,           &
          & reflectance = reflectance,     &
          & aer_opt_param = aer_opt_param, &
          & cld_opt_param = cld_opt_param, &
          & pccomp       = pccomp,         &
          & channels_rec = channels_rec,   &
          & nthreads = nthreads)
  ENDIF

  IF (errorstatus /= 0) THEN
    err = errorstatus_fatal
    THROW(err.ne.0)
  ENDIF

  reflout(:) = reflectance(:)%refl_out
  
  CALL make_chain_profile (err, chain_profiles, c, "PROFILES", profiles, .false._jplm)
  THROW(err.ne.0)

  CALL make_chain_profile (err, chain_profiles_inc, c_inc, "PROFILES_INC", profiles_inc, .false._jplm)
  THROW(err.ne.0)

  CALL make_chain_profile (err, chain_profiles_bf, c_bf, "PROFILES_BF", profiles_bf, .false._jplm)
  THROW(err.ne.0)

  CALL make_chain_profile (err, chain_profiles_k, c_k, "PROFILES_K", profiles_k, .true._jplm)
  THROW(err.ne.0)

!
  NULLIFY (chain_profiles_k_rec, c_k_in)
  nchannels_rec = -1_jpim

  NULLIFY (chain_profiles_k_pc, c_k_pc)
  npcscores = -1_jpim

  NULLIFY (chanprof_in, chanprof_pc)
  NULLIFY (w, w_in, w_pc)

  IF (opts%rt_ir%pc%addpc) THEN

    npcscores = size(pccomp%pcscores)

    IF (opts%rt_ir%pc%addradrec) THEN
!
! Reconstructed radiances jacobians
!
      nchannels_rec = size(pccomp%bt_pccomp)
      CALL make_chain_profile (err, chain_profiles_k_rec, c_k_in, "PROFILES_K_REC", profiles_k_rec, .true._jplm)
      THROW(err.ne.0)

      ALLOCATE (w_in(nchannels_rec), STAT = err)
      THROWM(err.ne.0,"Cannot allocate w_in")

    ELSE
!
! PC-scores jacobians
!
      CALL make_chain_profile (err, chain_profiles_k_pc, c_k_pc, "PROFILES_K_PC", profiles_k_pc, .true._jplm)
      THROW(err.ne.0)

      ALLOCATE (w_pc(npcscores), STAT = err)
      THROWM(err.ne.0,"Cannot allocate w_pc")

    ENDIF

    CALL rttov_alloc_pc_dimensions(err, opts, npcscores, nprofiles, chanprof_in, chanprof_pc, &
                                   1_jpim, channels_rec = channels_rec)
    THROW(err.ne.0)


!
! Allocate pccomp structures
!

    CALL rttov_alloc_pccomp (err, pccomp_bf, npcscores, 1_jpim, nchannels_rec = nchannels_rec)
    THROW(err.ne.0)

  ENDIF

  ALLOCATE (w(nchannels), STAT = err)
  THROWM(err.ne.0,"Cannot allocate w")


  CALL size_chain(nsize, chain_profiles(1))
  DO i = 1, nsize

    DO iprof = 1, nprofiles
      CALL get_pointer_chain(c(iprof)%p, x_profiles)
      CALL get_pointer_chain(c_bf(iprof)%p, x_profiles_bf)
      CALL get_pointer_chain(c_inc(iprof)%p, x_profiles_inc)
      x_profiles_bf = x_profiles + x_profiles_inc
      inc(iprof)    = x_profiles_inc
    ENDDO

    IF (nthreads <= 0) THEN
      CALL rttov_direct(                     &
            & errorstatus,                   &
            & chanprof,                      &
            & opts,                          &
            & profiles_bf,                   &
            & coefs,                         &
            & transmission_bf,               &
            & radiancedata_bf,               &
            & calcemis = calcemis,           &
            & emissivity = emissivity,       &
            & calcrefl = calcrefl,           &
            & reflectance = reflectance,     &
            & aer_opt_param = aer_opt_param, &
            & cld_opt_param = cld_opt_param, &
            & traj         = traj,           &
            & pccomp       = pccomp_bf,      &
            & channels_rec = channels_rec)
    ELSE
      CALL rttov_parallel_direct(            &
            & errorstatus,                   &
            & chanprof,                      &
            & opts,                          &
            & profiles_bf,                   &
            & coefs,                         &
            & transmission_bf,               &
            & radiancedata_bf,               &
            & calcemis = calcemis,           &
            & emissivity = emissivity,       &
            & calcrefl = calcrefl,           &
            & reflectance = reflectance,     &
            & aer_opt_param = aer_opt_param, &
            & cld_opt_param = cld_opt_param, &
            & pccomp       = pccomp_bf,      &
            & channels_rec = channels_rec,   &
            & nthreads     = nthreads)
    ENDIF

    IF (errorstatus /= 0) THEN
      err = errorstatus_fatal
      THROW(err.ne.0)
    ENDIF

    DO iprof = 1, nprofiles
      CALL get_pointer_chain(c(iprof)%p, x_profiles)
      CALL get_pointer_chain(c_bf(iprof)%p, x_profiles_bf)
      x_profiles_bf = x_profiles
    ENDDO

    DO iprof = 1, nprofiles
      CALL advance_chain(c(iprof)%p)
      CALL advance_chain(c_bf(iprof)%p)
      CALL advance_chain(c_inc(iprof)%p)
    ENDDO


    IF (opts%rt_ir%pc%addpc) THEN
      IF (opts%rt_ir%pc%addradrec) THEN
        WHERE (inc(chanprof_in(:)%prof) .NE. 0._jprb)
          w_in = 1._jprb / inc(chanprof_in(:)%prof)
        ELSEWHERE
          w_in = 0._jprb
        ENDWHERE
      ELSE
        WHERE (inc(chanprof_pc(:)%prof) .NE. 0._jprb)
          w_pc = 1._jprb / inc(chanprof_pc(:)%prof)
        ELSEWHERE
          w_pc = 0._jprb
        ENDWHERE
      ENDIF
    END IF
    WHERE (inc(chanprof(:)%prof) .NE. 0._jprb)
      w = 1._jprb / inc(chanprof(:)%prof)
    ELSEWHERE
      w = 0._jprb
    ENDWHERE

    

    IF (opts%rt_all%switchrad) THEN
      IF (opts%rt_ir%pc%addpc) THEN
        CALL assign_chain_profile (c_k, (radiancedata_bf%total-radiancedata%total)*w)
        IF (opts%rt_ir%pc%addradrec) THEN
          CALL assign_chain_profile (c_k_in, (pccomp_bf%bt_pccomp-pccomp%bt_pccomp)*w_in)
        ELSE
          CALL assign_chain_profile (c_k_pc, (pccomp_bf%pcscores-pccomp%pcscores)*w_pc)
        ENDIF
      ELSE
        DO ichan = 1, nchannels
          IF (coefs%coef%ss_val_chn(chanprof(ichan)%chan) < 2) THEN
            out_bf(ichan) = (radiancedata_bf%bt(ichan)-radiancedata%bt(ichan))*w(ichan)
          ELSE
            out_bf(ichan) = (radiancedata_bf%total(ichan)-radiancedata%total(ichan))*w(ichan)
          ENDIF
        ENDDO
        CALL assign_chain_profile (c_k, out_bf)
      ENDIF
    ELSE
      CALL assign_chain_profile (c_k, (radiancedata_bf%total-radiancedata%total)*w)
      IF (opts%rt_ir%pc%addpc) THEN
        IF (opts%rt_ir%pc%addradrec) THEN
          CALL assign_chain_profile (c_k_in, (pccomp_bf%total_pccomp-pccomp%total_pccomp)*w_in)
        ELSE
          CALL assign_chain_profile (c_k_pc, (pccomp_bf%pcscores-pccomp%pcscores)*w_pc)
        ENDIF
      ENDIF
    ENDIF

  ENDDO

  CALL free_chain_profile (err, chain_profiles, c)
  THROW(err.ne.0)

  CALL free_chain_profile (err, chain_profiles_bf, c_bf)
  THROW(err.ne.0)

  CALL free_chain_profile (err, chain_profiles_inc, c_inc)
  THROW(err.ne.0)

  CALL free_chain_profile (err, chain_profiles_k, c_k)
  THROW(err.ne.0)

  IF (opts%rt_ir%pc%addpc) THEN
    IF (opts%rt_ir%pc%addradrec) THEN
      CALL free_chain_profile (err, chain_profiles_k_rec, c_k_in)
      THROW(err.ne.0)
    ELSE
      CALL free_chain_profile (err, chain_profiles_k_pc, c_k_pc)
      THROW(err.ne.0)
    ENDIF
  ENDIF

  calcemis_bf = .FALSE.
  DO ichan = 1, nchannels
    emissivity_bf%emis_in         = emissivity%emis_out
    emissivity_bf(ichan)%emis_in  = emissivity_bf(ichan)%emis_in - 0.01_jprb

    IF (nthreads <= 0) THEN
      CALL rttov_direct(                     &
            & errorstatus,                   &
            & chanprof,                      &
            & opts,                          &
            & profiles,                      &
            & coefs,                         &
            & transmission_bf,               &
            & radiancedata_bf,               &
            & calcemis = calcemis_bf,        &
            & emissivity = emissivity_bf,    &
            & calcrefl = calcrefl,           &
            & reflectance = reflectance,     &
            & aer_opt_param = aer_opt_param, &
            & cld_opt_param = cld_opt_param, &
            & traj         = traj,           &
            & pccomp       = pccomp_bf,      &
            & channels_rec = channels_rec)
    ELSE
      CALL rttov_parallel_direct(            &
            & errorstatus,                   &
            & chanprof,                      &
            & opts,                          &
            & profiles,                      &
            & coefs,                         &
            & transmission_bf,               &
            & radiancedata_bf,               &
            & calcemis = calcemis_bf,        &
            & emissivity = emissivity_bf,    &
            & calcrefl = calcrefl,           &
            & reflectance = reflectance,     &
            & aer_opt_param = aer_opt_param, &
            & cld_opt_param = cld_opt_param, &
            & pccomp       = pccomp_bf,      &
            & channels_rec = channels_rec,   &
            & nthreads     = nthreads)
    ENDIF

    IF (errorstatus /= 0) THEN
      err = errorstatus_fatal
      THROW(err.ne.0)
    ENDIF

    IF (opts%rt_all%switchrad .and. .not. opts%rt_ir%pc%addpc .and. &
        coefs%coef%ss_val_chn(chanprof(ichan)%chan) < 2) THEN
      emissivity_k(ichan)%emis_in = (radiancedata_bf%bt(ichan) - radiancedata%bt(ichan)) / ( - 0.01_jprb)
    ELSE
      emissivity_k(ichan)%emis_in = (radiancedata_bf%total(ichan) - radiancedata%total(ichan)) / ( - 0.01_jprb)
    ENDIF
  ENDDO

  IF (opts%rt_ir%addsolar) THEN
  
    ! Reflectances may depend on profile/emis values when calcrefl is true. Therefore we
    ! use the saved reflectances from the original call to rttov_direct.
    
    calcrefl_bf(:) = .FALSE.
    DO ichan = 1, nchannels
      reflectance_bf%refl_cloud_top = 0._jprb
      reflectance_bf%refl_in        = reflout
      reflectance_bf(ichan)%refl_in = reflectance_bf(ichan)%refl_in + 0.01_jprb

      IF (nthreads <= 0) THEN
        CALL rttov_direct(                     &
              & errorstatus,                   &
              & chanprof,                      &
              & opts,                          &
              & profiles,                      &
              & coefs,                         &
              & transmission_bf,               &
              & radiancedata_bf,               &
              & calcemis = calcemis,           &
              & emissivity = emissivity,       &
              & calcrefl = calcrefl_bf,        &
              & reflectance = reflectance_bf,  &
              & aer_opt_param = aer_opt_param, &
              & cld_opt_param = cld_opt_param, &
              & traj         = traj,           &
              & pccomp       = pccomp_bf,      &
              & channels_rec = channels_rec)
      ELSE
        CALL rttov_parallel_direct(            &
              & errorstatus,                   &
              & chanprof,                      &
              & opts,                          &
              & profiles,                      &
              & coefs,                         &
              & transmission_bf,               &
              & radiancedata_bf,               &
              & calcemis = calcemis,           &
              & emissivity = emissivity,       &
              & calcrefl = calcrefl_bf,        &
              & reflectance = reflectance_bf,  &
              & aer_opt_param = aer_opt_param, &
              & cld_opt_param = cld_opt_param, &
              & pccomp       = pccomp_bf,      &
              & channels_rec = channels_rec,   &
              & nthreads     = nthreads)
      ENDIF

      IF (errorstatus /= 0) THEN
        err = errorstatus_fatal
        THROW(err.ne.0)
      ENDIF
      
      IF (opts%rt_all%switchrad .and. &
        coefs%coef%ss_val_chn(chanprof(ichan)%chan) < 2) THEN
        reflectance_k(ichan)%refl_in = (radiancedata_bf%bt(ichan) - radiancedata%bt(ichan)) / ( +0.01_jprb)
      ELSE
        reflectance_k(ichan)%refl_in = (radiancedata_bf%total(ichan) - radiancedata%total(ichan)) / ( +0.01_jprb)
      ENDIF
    ENDDO

  ENDIF
  
    
  IF (opts%rt_ir%pc%addpc) THEN

    CALL rttov_alloc_pccomp (err, pccomp_bf, npcscores, 0_jpim, nchannels_rec = nchannels_rec)
    THROW(err.ne.0)

    CALL rttov_alloc_pc_dimensions(err, opts, npcscores, nprofiles, chanprof_in, chanprof_pc, &
                                   0_jpim)
    THROW(err.ne.0)

    IF (opts%rt_ir%pc%addradrec) THEN
      DEALLOCATE (w_in, STAT = err)
      THROWM(err.ne.0,"Cannot deallocate w_in")
    ELSE
      DEALLOCATE (w_pc, STAT = err)
      THROWM(err.ne.0,"Cannot deallocate w_pc")
    ENDIF

  ENDIF

  DEALLOCATE (w, STAT = err)
  THROWM(err.ne.0,"Cannot deallocate w")

  IF (nthreads <= 0) THEN
    CALL rttov_alloc_traj(                          &
        & err,                  nprofiles,          &
        & nchannels,            opts,               &
        & profiles(1)%nlevels,  coefs,              &
        & 0_jpim,               traj = traj)
    THROW(err.ne.0)
  ENDIF

  CALL rttov_alloc_transmission(err, transmission_bf, profiles(1)%nlevels-1_jpim, nchannels, 0_jpim)
  THROW(err.ne.0)

  CALL rttov_alloc_rad( &
        & err,                          &
        & nchannels,                    &
        & radiancedata_bf,              &
        & profiles(1)%nlevels - 1_jpim, &
        & 0_jpim)
  THROW(err.ne.0)

  CALL rttov_alloc_prof( &
        & err,                 &
        & nprofiles,           &
        & profiles_inc,        &
        & profiles(1)%nlevels, &
        & opts,                &
        & 0_jpim)
  THROW(err.ne.0)

  CALL rttov_alloc_prof( &
        & err,                 &
        & nprofiles,           &
        & profiles_bf,         &
        & profiles(1)%nlevels, &
        & opts,                &
        & 0_jpim)
  THROW(err.ne.0)

CATCH

  errorstatus = err


END SUBROUTINE rttov_k_bf
