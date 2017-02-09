SUBROUTINE rttov_alloc_traj_sta (err, traj_sta, opts, coef, nlevels, nchannels, nprofiles, asw, npcscores, channels_rec)
! Description:
!   Allocates/deallocates static trajectory data
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
#include "throw.h"

USE rttov_types, ONLY : &
 & rttov_traj_sta, &
 & rttov_coef,     &
 & rttov_options

USE parkind1, ONLY : &
  & jpim

IMPLICIT NONE


INTEGER(KIND=jpim),   INTENT(OUT)   :: err
TYPE(rttov_traj_sta), INTENT(INOUT) :: traj_sta
TYPE(rttov_options),  INTENT(IN)    :: opts
TYPE(rttov_coef),     INTENT(IN)    :: coef
INTEGER(KIND=jpim),   INTENT(IN)    :: nlevels
INTEGER(KIND=jpim),   INTENT(IN)    :: nchannels
INTEGER(KIND=jpim),   INTENT(IN)    :: nprofiles
INTEGER(KIND=jpim),   INTENT(IN)    :: asw
INTEGER(KIND=jpim),   INTENT(IN), OPTIONAL :: npcscores
INTEGER(KIND=jpim),   INTENT(IN), OPTIONAL :: channels_rec(:)

!INTF_END

#include "rttov_alloc_prof.interface"
#include "rttov_alloc_auxrad.interface"
#include "rttov_alloc_pc_dimensions.interface"

INTEGER(KIND=jpim)  :: nlayers
TYPE(rttov_options) :: opts_coef


TRY

  nlayers = nlevels - 1

  opts_coef = opts
  opts_coef%rt_ir%addclouds = .FALSE.
  opts_coef%rt_ir%addaerosl = .FALSE.

  IF (asw == 1_jpim) THEN

    ALLOCATE (    &
            & traj_sta%do_lambertian(nchannels),                             &
            & traj_sta%thermal(nchannels),                                   &
            & traj_sta%solar(nchannels),                                     &
            & traj_sta%angles(nprofiles),                                    &
            & traj_sta%angles_coef(nprofiles),                               &
            ! These are passed into rttov_integrate so are required even when addsolar is false:
            & traj_sta%solar_spec_esd(nchannels),                            &
            & traj_sta%refl_norm(nchannels),                                 &
            & STAT = err)
    THROWM(err.NE.0,"Allocation of traj_sta failed")

    ! At this point we don't know if thermal calculations are required, but actually if not
    ! the thermal_path1 arrays are not necessary
    ALLOCATE(traj_sta%thermal_path1, STAT=err)
    THROWM(err.NE.0,"Allocation of traj_sta for thermal path1 failed")
    ALLOCATE (    &
            & traj_sta%thermal_path1%tau_ref(nlevels, nchannels),            &
            & traj_sta%thermal_path1%tau_ref_surf(nchannels),                &
            & traj_sta%thermal_path1%tau_surf(nchannels),                    &
            & traj_sta%thermal_path1%opdp_ref_coef(coef%nlayers, nchannels), &
            & traj_sta%thermal_path1%od_level(nlevels, nchannels),           &
            & traj_sta%thermal_path1%tau_level(nlevels, nchannels),          &
            & STAT = err)
    THROWM(err.NE.0,"Allocation of traj_sta for thermal path1 failed")

    IF (opts%rt_ir%addsolar) THEN
      ALLOCATE(traj_sta%solar_path2, STAT=err)
      THROWM(err.NE.0,"Allocation of traj_sta for solar path2 failed")
      ALLOCATE (  &
            & traj_sta%solar_path2%tau_ref(nlevels, nchannels),              &
            & traj_sta%solar_path2%tau_ref_surf(nchannels),                  &
            & traj_sta%solar_path2%tau_level(nlevels, nchannels),            &
            & traj_sta%solar_path2%tau_surf(nchannels),                      &
            & traj_sta%solar_path2%opdp_ref_coef(coef%nlayers , nchannels),  &
            & traj_sta%solar_path2%od_level(nlevels, nchannels),             &
            & traj_sta%solar_path2%od_singlelayer(nlayers, nchannels),       &
            & traj_sta%solar_path2%od_frac(nchannels),                       &
            & STAT = err)
      THROWM(err.NE.0,"Allocation of traj_sta for solar path2 failed")

      ALLOCATE(traj_sta%solar_path1, STAT=err)
      THROWM(err.NE.0,"Allocation of traj_sta for solar path1 failed")
      ALLOCATE (  &
            & traj_sta%solar_path1%tau_ref(nlevels, nchannels),              &
            & traj_sta%solar_path1%tau_ref_surf(nchannels),                  &
            & traj_sta%solar_path1%tau_level(nlevels, nchannels),            &
            & traj_sta%solar_path1%tau_surf(nchannels),                      &
            & STAT = err)
      THROWM(err.NE.0,"Allocation of traj_sta for solar path1 failed")
    ENDIF

    IF (opts%config%apply_reg_limits .OR. &
        opts%interpolation%reg_limit_extrap) THEN

      ALLOCATE(traj_sta%profiles_coef_ref(nprofiles), STAT = err)
      THROWM(err.NE.0,"Allocation of traj_sta%profiles_coef_ref failed")

      CALL rttov_alloc_prof(err, nprofiles, traj_sta%profiles_coef_ref, coef%nlevels, opts_coef, 1_jpim)
      THROW(err.NE.0)

    ENDIF

    CALL rttov_alloc_auxrad( &
          & err,                 &
          & traj_sta%auxrad,     &
          & nlevels,             &
          & nchannels,           &
          & 1_jpim)
    THROW(err.NE.0)

    IF (opts%rt_ir%pc%addpc) THEN

      CALL rttov_alloc_pc_dimensions( &
            & err,                   &
            & opts,                  &
            & npcscores,             &
            & nprofiles,             &
            & traj_sta%chanprof_in,  &
            & traj_sta%chanprof_pc,  &
            & 1_jpim,                &
            & channels_rec = channels_rec)

     THROW(err.NE.0)

    ENDIF

  ENDIF

  IF (asw == 0_jpim) THEN

    IF (opts%rt_ir%pc%addpc) THEN

      CALL rttov_alloc_pc_dimensions( &
            & err,                   &
            & opts,                  &
            & npcscores,             &
            & nprofiles,             &
            & traj_sta%chanprof_in,  &
            & traj_sta%chanprof_pc,  &
            & 0_jpim,                &
            & channels_rec = channels_rec)

     THROW(err.NE.0)

    ENDIF

    CALL rttov_alloc_auxrad( &
          & err,                 &
          & traj_sta%auxrad,     &
          & nlevels,             &
          & nchannels,           &
          & 0_jpim)
    THROW(err.NE.0)

    IF (opts%config%apply_reg_limits .OR. &
        opts%interpolation%reg_limit_extrap) THEN

      CALL rttov_alloc_prof (err, nprofiles, traj_sta%profiles_coef_ref, coef % nlevels, opts_coef, 0_jpim)
      THROW(err.NE.0)

      DEALLOCATE(traj_sta%profiles_coef_ref, STAT = err)
      THROWM(err.NE.0,"Deallocation of traj_sta%profiles_coef_ref failed")

    ENDIF

    DEALLOCATE (    &
            & traj_sta%do_lambertian,               &
            & traj_sta%thermal,                     &
            & traj_sta%solar,                       &
            & traj_sta%angles,                      &
            & traj_sta%angles_coef,                 &
            & traj_sta%solar_spec_esd,              &
            & traj_sta%refl_norm,                   &
            & STAT = err)
    THROWM(err.NE.0,"Deallocation of traj_sta failed")
    DEALLOCATE (    &
            & traj_sta%thermal_path1%tau_ref,       &
            & traj_sta%thermal_path1%tau_ref_surf,  &
            & traj_sta%thermal_path1%tau_surf,      &
            & traj_sta%thermal_path1%opdp_ref_coef, &
            & traj_sta%thermal_path1%od_level,      &
            & traj_sta%thermal_path1%tau_level,     &
            & STAT = err)
    THROWM(err.NE.0,"Deallocation of traj_sta for thermal path1 failed")
    DEALLOCATE(traj_sta%thermal_path1, STAT=err)
    THROWM(err.NE.0,"Deallocation of traj_sta for thermal path1 failed")

    IF (opts%rt_ir%addsolar) THEN
      DEALLOCATE (  &
            & traj_sta%solar_path2%tau_ref,        &
            & traj_sta%solar_path2%tau_ref_surf,   &
            & traj_sta%solar_path2%tau_level,      &
            & traj_sta%solar_path2%tau_surf,       &
            & traj_sta%solar_path2%opdp_ref_coef,  &
            & traj_sta%solar_path2%od_level,       &
            & traj_sta%solar_path2%od_singlelayer, &
            & traj_sta%solar_path2%od_frac,        &
            & STAT = err)
      THROWM(err.NE.0,"Deallocation of traj_sta for solar path2 failed")
      DEALLOCATE(traj_sta%solar_path2, STAT=err)
      THROWM(err.NE.0,"Deallocation of traj_sta for solar path2 failed")

      DEALLOCATE (  &
            & traj_sta%solar_path1%tau_ref,        &
            & traj_sta%solar_path1%tau_ref_surf,   &
            & traj_sta%solar_path1%tau_level,      &
            & traj_sta%solar_path1%tau_surf,       &
            & STAT = err)
      THROWM(err.NE.0,"Deallocation of traj_sta for solar path1 failed")
      DEALLOCATE(traj_sta%solar_path1, STAT=err)
      THROWM(err.NE.0,"Deallocation of traj_sta for solar path1 failed")
    ENDIF

  ENDIF

CATCH

END SUBROUTINE
