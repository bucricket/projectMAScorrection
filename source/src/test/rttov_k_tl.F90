!
SUBROUTINE rttov_k_tl( &
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
!   Compute jacobian using tangent linear code
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

#include "rttov_tl.interface"
#include "rttov_parallel_tl.interface"
#include "rttov_alloc_prof.interface"
#include "rttov_alloc_pccomp.interface"
#include "rttov_alloc_traj.interface"
#include "rttov_alloc_rad.interface"
#include "rttov_alloc_transmission.interface"
#include "rttov_init_rad.interface"
#include "rttov_init_pccomp.interface"

  TYPE(rttov_traj)        :: traj, traj_tl
  TYPE(rttov_pccomp)      :: pccomp_tl
  TYPE(profile_Type), TARGET :: profiles_tl(size(profiles))
  TYPE(rttov_emissivity)     :: emissivity_tl(size(chanprof))
  TYPE(rttov_reflectance)    :: reflectance_tl(size(chanprof))
  TYPE(transmission_Type) :: transmission_tl
  TYPE(radiance_Type    ) :: radiancedata_tl
  REAL(KIND=jprb)         :: out_tl(size(chanprof))
  TYPE(chain),  POINTER :: chain_profiles_tl(:)
  TYPE(chain),  POINTER :: chain_profiles_k(:)
  TYPE(chain),  POINTER :: chain_profiles_k_pc(:)
  TYPE(chain),  POINTER :: chain_profiles_k_rec(:)
  TYPE(pchain), POINTER :: c_tl(:)
  TYPE(pchain), POINTER :: c_k(:)
  TYPE(pchain), POINTER :: c_k_pc(:)
  TYPE(pchain), POINTER :: c_k_in(:)
  INTEGER(KIND=jpim) :: err, i, nsize, ichan, iprof
  REAL(KIND=jprb), POINTER :: x
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: nchannels_rec
  INTEGER(KIND=jpim) :: npcscores
  LOGICAL(KIND=jplm) :: calcreflf(size(chanprof))
  TYPE(rttov_reflectance) :: reflectance_save(size(chanprof))
TRY

  nprofiles = size(profiles)
  nchannels = size(chanprof)

  CALL rttov_alloc_prof( &
        & err,                                                           &
        & nprofiles,                                                     &
        & profiles_tl,                                                   &
        & profiles(1)%nlevels,                                           &
        & opts,                                                          &
        & 1_jpim,                                                        &
        & coefs = coefs,                                                 &
        & init = .TRUE._jplm)
  THROW(err.ne.0)

  CALL rttov_alloc_rad( &
        & err,                          &
        & nchannels,                    &
        & radiancedata_tl,              &
        & profiles(1)%nlevels - 1_jpim, &
        & 1_jpim)
  THROW(err.ne.0)

  CALL rttov_alloc_transmission(err, transmission_tl, profiles(1)%nlevels-1_jpim, nchannels, 1_jpim)
  THROW(err.ne.0)

  IF (nthreads <= 0) THEN
    CALL rttov_alloc_traj(                          &
        & err,                  nprofiles,          &
        & nchannels,            opts,               &
        & profiles(1)%nlevels,  coefs,              &
        & 1_jpim,               traj = traj,        &
        & traj_tl = traj_tl )
    THROW(err.ne.0)
  ENDIF

!
! Chain input profiles
!

  CALL make_chain_profile (err, chain_profiles_tl, c_tl, "PROFILES_TL", profiles_tl, .false._jplm)
  THROW(err.ne.0)

!
! Regular radiances jacobians
!
  CALL make_chain_profile (err, chain_profiles_k, c_k, "PROFILES_K", profiles_k, .true._jplm)
  THROW(err.ne.0)

!
  NULLIFY (chain_profiles_k_rec, c_k_in)
  nchannels_rec = -1_jpim

  NULLIFY (chain_profiles_k_pc, c_k_pc)
  npcscores = -1_jpim

  IF (opts%rt_ir%pc%addpc) THEN

    npcscores = size(pccomp%pcscores)

    IF (opts%rt_ir%pc%addradrec) THEN
!
! Reconstructed radiances jacobians
!
      nchannels_rec = size(pccomp%bt_pccomp)
      CALL make_chain_profile (err, chain_profiles_k_rec, c_k_in, "PROFILES_K_REC", profiles_k_rec, .true._jplm)
      THROW(err.ne.0)
    ELSE
!
! PC-scores jacobians
!
      CALL make_chain_profile (err, chain_profiles_k_pc, c_k_pc, "PROFILES_K_PC", profiles_k_pc, .true._jplm)
      THROW(err.ne.0)

    ENDIF

!
! Allocate a pccomp_tl
!
    CALL rttov_alloc_pccomp (err, pccomp_tl, npcscores, 1_jpim, nchannels_rec = nchannels_rec)
    THROW(err.ne.0)
  ENDIF

  CALL size_chain(nsize, chain_profiles_tl(1))! assume all profiles have the same size

  DO i = 1, nsize
    emissivity_tl%emis_in      = 0._jprb
    emissivity_tl%emis_out     = 0._jprb
    reflectance_tl%refl_in     = 0._jprb
    reflectance_tl%refl_out    = 0._jprb

    CALL rttov_init_rad (radiancedata_tl)
    IF (opts%rt_ir%pc%addpc) THEN
      CALL rttov_init_pccomp (pccomp_tl)
    ENDIF

    DO iprof = 1, nprofiles
      CALL get_pointer_chain(c_tl(iprof)%p, x)
      x = 1._jprb
!      Print *, '-------', i
!      Call print_chain( 6, chain_profiles_tl(iprof) )
    ENDDO

    IF (nthreads <= 0) THEN
      CALL rttov_tl( &
            & errorstatus,               &
            & chanprof,                  &
            & opts,                      &
            & profiles,                  &
            & profiles_tl,               &
            & coefs,                     &
            & transmission,              &
            & transmission_tl,           &
            & radiancedata,              &
            & radiancedata_tl,           &
            & calcemis,                  &
            & emissivity,                &
            & emissivity_tl,             &
            & calcrefl,                  &
            & reflectance,               &
            & reflectance_tl,            &
            & aer_opt_param,             &
            & cld_opt_param,             &
            & traj         = traj,       &
            & traj_tl      = traj_tl,    &
            & pccomp       = pccomp,     &
            & pccomp_tl    = pccomp_tl,  &
            & channels_rec = channels_rec)
    ELSE
      CALL rttov_parallel_tl( &
            & errorstatus,               &
            & chanprof,                  &
            & opts,                      &
            & profiles,                  &
            & profiles_tl,               &
            & coefs,                     &
            & transmission,              &
            & transmission_tl,           &
            & radiancedata,              &
            & radiancedata_tl,           &
            & calcemis,                  &
            & emissivity,                &
            & emissivity_tl,             &
            & calcrefl,                  &
            & reflectance,               &
            & reflectance_tl,            &
            & aer_opt_param,             &
            & cld_opt_param,             &
            & pccomp       = pccomp,     &
            & pccomp_tl    = pccomp_tl,  &
            & channels_rec = channels_rec, &
            & nthreads     = nthreads)
    ENDIF

    IF (errorstatus /= 0) THEN
      err = errorstatus_fatal
      THROW(err.ne.0)
    ENDIF

    DO iprof = 1, nprofiles
      CALL get_pointer_chain(c_tl(iprof)%p, x)
      x = 0._jprb
      CALL advance_chain(c_tl(iprof)%p)
    ENDDO

    IF (opts%rt_all%switchrad) THEN
      IF (opts%rt_ir%pc%addpc) THEN
        CALL assign_chain_profile (c_k, radiancedata_tl%total)
        IF (opts%rt_ir%pc%addradrec) THEN
          CALL assign_chain_profile (c_k_in, pccomp_tl%bt_pccomp)
        ELSE
          CALL assign_chain_profile (c_k_pc, pccomp_tl%pcscores)
        ENDIF
      ELSE
        DO ichan = 1, nchannels
          IF (coefs%coef%ss_val_chn(chanprof(ichan)%chan) < 2) THEN
            out_tl(ichan) = radiancedata_tl%bt(ichan)
          ELSE
            out_tl(ichan) = radiancedata_tl%total(ichan)
          ENDIF
        ENDDO
        CALL assign_chain_profile (c_k, out_tl)
      ENDIF
    ELSE
      CALL assign_chain_profile (c_k, radiancedata_tl%total)
      IF (opts%rt_ir%pc%addpc) THEN
        IF (opts%rt_ir%pc%addradrec) THEN
          CALL assign_chain_profile (c_k_in, pccomp_tl%total_pccomp)
        ELSE
          CALL assign_chain_profile (c_k_pc, pccomp_tl%pcscores)
        ENDIF
      ENDIF
    ENDIF

   
  ENDDO

  CALL free_chain_profile (err, chain_profiles_tl, c_tl)
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

  ! Store any computed surface values
  reflectance_save%refl_in = reflectance%refl_out
  reflectance_save%refl_cloud_top = reflectance%refl_cloud_top

  DO ichan = 1, nchannels
    emissivity_tl%emis_in        = 0._jprb
    emissivity_tl(ichan)%emis_in = 1._jprb
    reflectance_tl%refl_in       = 0._jprb

    IF (nthreads <= 0) THEN
      CALL rttov_tl( &
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
            & calcemis,                  &
            & emissivity,                &
            & emissivity_tl,             &
            & calcrefl,                  &
            & reflectance,               &
            & reflectance_tl,            &
            & aer_opt_param,             &
            & cld_opt_param,             &
            & traj         = traj,       &
            & traj_tl      = traj_tl,    &
            & pccomp       = pccomp,     &
            & pccomp_tl    = pccomp_tl,  &
            & channels_rec = channels_rec)
    ELSE
      CALL rttov_parallel_tl( &
            & errorstatus,               &
            & chanprof,                  &
            & opts,                      &
            & profiles,                  &
            & profiles_tl,               &
            & coefs,                     &
            & transmission,              &
            & transmission_tl,           &
            & radiancedata,              &
            & radiancedata_tl,           &
            & calcemis,                  &
            & emissivity,                &
            & emissivity_tl,             &
            & calcrefl,                  &
            & reflectance,               &
            & reflectance_tl,            &
            & aer_opt_param,             &
            & cld_opt_param,             &
            & pccomp       = pccomp,     &
            & pccomp_tl    = pccomp_tl,  &
            & channels_rec = channels_rec, &
            & nthreads     = nthreads)
    ENDIF

    IF (errorstatus /= 0) THEN
      err = errorstatus_fatal
      THROW(err.ne.0)
    ENDIF

    IF (opts%rt_all%switchrad) THEN
      IF (coefs%coef%ss_val_chn(chanprof(ichan)%chan) < 2) THEN
        emissivity_k(ichan)%emis_in = radiancedata_tl%bt(ichan)
      ELSE
        emissivity_k(ichan)%emis_in = radiancedata_tl%total(ichan)
      ENDIF
    ELSE
      emissivity_k(ichan)%emis_in = radiancedata_tl%total(ichan)
    ENDIF
  ENDDO  

  
  IF (opts%rt_ir%addsolar) THEN

    calcreflf(:) = .FALSE.

    DO ichan = 1, nchannels
      emissivity_tl%emis_in         = 0._jprb
      reflectance_tl%refl_in        = 0._jprb
      reflectance_tl(ichan)%refl_in = 1._jprb

      IF (nthreads <= 0) THEN
        CALL rttov_tl( &
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
              & calcemis,                  &
              & emissivity,                &
              & emissivity_tl,             &
              & calcreflf,                 &
              & reflectance_save,          &
              & reflectance_tl,            &
              & aer_opt_param,             &
              & cld_opt_param,             &
              & traj         = traj,       &
              & traj_tl      = traj_tl,    &
              & pccomp       = pccomp,     &
              & pccomp_tl    = pccomp_tl,  &
              & channels_rec = channels_rec)
      ELSE
        CALL rttov_parallel_tl( &
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
              & calcemis,                  &
              & emissivity,                &
              & emissivity_tl,             &
              & calcreflf,                 &
              & reflectance_save,          &
              & reflectance_tl,            &
              & aer_opt_param,             &
              & cld_opt_param,             &
              & pccomp       = pccomp,     &
              & pccomp_tl    = pccomp_tl,  &
              & channels_rec = channels_rec, &
              & nthreads     = nthreads)
      ENDIF

      IF (errorstatus /= 0) THEN
        err = errorstatus_fatal
        THROW(err.ne.0)
      ENDIF
  
      IF (opts%rt_all%switchrad) THEN
        IF (coefs%coef%ss_val_chn(chanprof(ichan)%chan) < 2) THEN
          reflectance_k(ichan)%refl_in = radiancedata_tl%bt(ichan)
        ELSE
          reflectance_k(ichan)%refl_in = radiancedata_tl%total(ichan)
        ENDIF
      ELSE
        reflectance_k(ichan)%refl_in = radiancedata_tl%total(ichan)
      ENDIF
    ENDDO  
  
  ENDIF

  IF (opts%rt_ir%pc%addpc) THEN
    CALL rttov_alloc_pccomp (err, pccomp_tl, npcscores, 0_jpim, nchannels_rec = nchannels_rec)
    THROW(err.ne.0)
  ENDIF

  IF (nthreads <= 0) THEN
    CALL rttov_alloc_traj(                          &
        & err,                  nprofiles,          &
        & nchannels,            opts,               &
        & profiles(1)%nlevels,  coefs,              &
        & 0_jpim,               traj = traj,        &
        & traj_tl = traj_tl )
    THROW(err.ne.0)
  ENDIF

  CALL rttov_alloc_transmission(err, transmission_tl, profiles(1)%nlevels-1_jpim, nchannels, 0_jpim)
  THROWM(err.ne.0,"Deallocation of transmission_tl failed")

  CALL rttov_alloc_rad( &
        & err,                          &
        & nchannels,                    &
        & radiancedata_tl,              &
        & profiles(1)%nlevels - 1_jpim, &
        & 0_jpim)
  THROW(err.ne.0)

  CALL rttov_alloc_prof( &
        & err,                 &
        & nprofiles,           &
        & profiles_tl,         &
        & profiles(1)%nlevels, &
        & opts,                &
        & 0_jpim)
  THROW(err.ne.0)


CATCH

  errorstatus = err

END SUBROUTINE rttov_k_tl

