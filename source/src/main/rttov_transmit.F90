!
SUBROUTINE rttov_transmit( &
            & addaerosl,                    &
            & addclouds,                    &
            & do_lambertian,                &
            & nlayers,                      &
            & chanprof,                     &
            & chanflag,                     &
            & aux,                          &
            & coef,                         &
            & ircld,                        &
            & geometry,                     &
            & opdp_path,                    &
            & od_level,                     &
            & transmission_levels,          &
            & transmission_total,           &
            & transmission_aux,             &
            & transmission_aux_path,        &
            & transmission_scatt_ir,        &
            & transmission_scatt_ir_stream, &
            & tau_ref,                      &
            & tau_ref_surf,                 &
            & tau_surf,                     &
            & tau_level)
!
! Description:
!  To calculate optical depths for a number of channels
!  and profiles from every pressure level to space.
! To interpolate optical depths on to levels of radiative transfer model
! (which, at present, entails only surface transmittance, as
! other optical depths are on *rt* levels) and to convert
! optical depths to transmittances.
!
! Code based on OPDEP and RTTAU from previous versions of RTTOV
! Only one profile per call
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
!    Copyright 2002, EUMETSAT, All Rights Reserved.
!
! Method:
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0    01/06/2005  Marco Matricardi (ECMWF):
!            --       New routine based on rttov_transmit.F90.
!                     Variable trace gases, CO2, N2O,CO and CH4
!                     have been introduced for AIRS and IASI
!  1.1    29/01/2007  Removed polarisation R Saunders
!  1.2    22/08/2007  Optimised (D Salmond)
!  1.3    04/06/2008  Fix od_frac and od_frac_ac calculation near surface level (PB PM)
!  1.4    27/02/2009  Profile levels to include ToA. Distinguish between
!                     layer arrays and level arrays - SIZE, index
!                     labels, looping (P. Rayer)
!  1.5    03/11/2009  Transmittances on levels (A Geer)
!  1.6    08/03/2010  Missing levsurf assignment (N Bormann)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud

  USE rttov_types, ONLY :  &
       rttov_chanprof,             &
       rttov_coef,                 &
       transmission_type_aux,      &
       rttov_path_transmission,    &
       transmission_scatt_ir_type, &
       profile_aux,                &
       ircld_type,                 &
       geometry_type

  USE parkind1, ONLY : jpim, jprb, jplm
  USE rttov_math_mod
!INTF_OFF
  USE rttov_const, ONLY : sensor_id_hi, &
    max_optical_depth, min_od, min_tau, &
    sec_theta_eff

  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE

  LOGICAL(KIND=jplm)              , INTENT(IN)    :: addaerosl
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: addclouds
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: do_lambertian(:)
  INTEGER(KIND=jpim)              , INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof            ), INTENT(IN)    :: chanprof(:)
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: chanflag(:)
  TYPE(rttov_coef                ), INTENT(IN)    :: coef
  REAL(KIND=jprb)                 , INTENT(IN)    :: opdp_path(:,:)
  REAL(KIND=jprb)                 , INTENT(INOUT) :: transmission_levels(:,:)   ! Transmittances on levels
  REAL(KIND=jprb)                 , INTENT(INOUT) :: transmission_total(:)      ! Total transmittances
  TYPE(transmission_type_aux     ), INTENT(INOUT) :: transmission_aux           ! Top-level arrays
  TYPE(rttov_path_transmission   ), INTENT(INOUT) :: transmission_aux_path      ! Transmittances and single-layer od
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir_stream
  TYPE(profile_aux               ), INTENT(IN)    :: aux
  TYPE(ircld_type                ), INTENT(IN)    :: ircld
  TYPE(geometry_type),              INTENT(in)    :: geometry(:)
  REAL(KIND=jprb), INTENT(INOUT) :: tau_ref     (nlayers + 1, SIZE(chanprof))
  REAL(KIND=jprb), INTENT(INOUT) :: tau_ref_surf(SIZE(chanprof))
  REAL(KIND=jprb), INTENT(INOUT) :: od_level    (nlayers + 1, SIZE(chanprof))! sat to level optical depth
  REAL(KIND=jprb), INTENT(INOUT) :: tau_level   (nlayers + 1, SIZE(chanprof))! sat to level transmittance
  REAL(KIND=jprb), INTENT(INOUT) :: tau_surf    (SIZE(chanprof))! sat to surface transmittance
!INTF_END

  REAL   (KIND=jprb) :: od_surf(SIZE(chanprof))                 ! sat to surface optical depth
  REAL   (KIND=jprb) :: od_surf_ac(0:SIZE(transmission_aux_path%tau_level(1,:,1)), SIZE(chanprof))
  REAL   (KIND=jprb) :: od_frac(SIZE(chanprof)), od_frac_ac
  ! pgf90 doesn't like opdpac(nlayers+1) with do_lambertian + IR cld/aer simulations:
  REAL   (KIND=jprb) :: opdpac(SIZE(transmission_aux_path%tau_level(:,0,1)))
  REAL   (KIND=jprb) :: opdpsurfac
  REAL   (KIND=jprb) :: small_val
  REAL   (KIND=jprb) :: ref_power
  INTEGER(KIND=jpim) :: lev, lay, chan, i, j, prof, ist, isti, klevels
  INTEGER(KIND=jpim) :: nlevels
  INTEGER(KIND=jpim) :: levsurf
  INTEGER(KIND=jpim) :: nchannels

! allocate automatic arrays using SIZE(transmission_aux_path%tau_level(1,:,1) to keep pgf90 compiler happy
  REAL   (KIND=jprb) :: ztemp(nlayers+1, 0:SIZE(transmission_aux_path%tau_level(1,:,1)))
  LOGICAL(kind=jplm) :: all_tau_level_ge_zero

  REAL   (KIND=JPRB) :: ZHOOK_HANDLE

!- End of header --------------------------------------------------------

!--------------------------------------------------------------
!1. Assemble layer optical depths and convert to transmittances
!--------------------------------------------------------------
! - optical depths here are negative except for od_singlelayer
! - in rttov_opdep, already checked that value of opticaldepth is sensible
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT', 0_jpim, ZHOOK_HANDLE)
  small_val = (TINY(1._jprb)) ** (0.333333_jprb) ! XLF doesn't like 1/3
  nchannels = SIZE(chanprof)
  nlevels   = nlayers + 1
  DO j = 1, nchannels
    IF (chanflag(j)) THEN
      DO lev = 1, nlevels
        od_level(lev,j) = 0._jprb
      ENDDO
    ENDIF
  ENDDO
! single layer optical depths local variable
  DO j = 1, nchannels
    IF (chanflag(j)) THEN
      DO lay = 1, nlayers
        transmission_aux_path%od_singlelayer(:,lay,j) = - (opdp_path(lay+1,j) - opdp_path(lay,j))
      ENDDO
    ENDIF
  ENDDO
! level to space optical depths - local variable
! Introduce klevels/layers to stop NEC compiler merging next 2 loops
  klevels = nlevels
! gamma correction of local variables
  IF(ANY(coef%ff_gam(:) .NE. 1.0_jprb)) THEN
!cdir NODEP
!cdir COLLAPSE
    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        chan = chanprof(j)%chan
        DO lev = 1, klevels
          od_level(lev, j) = MAX(coef%ff_gam(chan) * opdp_path(lev, j), -max_optical_depth)
        ENDDO
        transmission_aux_path%od_singlelayer(:,:,j) = &
          coef%ff_gam(chan) * transmission_aux_path%od_singlelayer(:,:,j)
      ENDIF
    ENDDO
  ELSE
!cdir NODEP
!cdir COLLAPSE
    DO j=1,nchannels
      IF (chanflag(j)) THEN
        DO lev = 1, nlevels
          od_level(lev, j) = MAX(opdp_path(lev, j), -max_optical_depth)
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  IF (ALL(chanflag)) THEN
    CALL exponential(od_level(:, :), tau_ref(:, :))
  ELSE
    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        CALL exponential(od_level(:, j), tau_ref(:, j))
      ENDIF
    ENDDO
  ENDIF

  DO j = 1, nchannels
    IF (chanflag(j)) THEN
      DO lev = 1, nlevels
        tau_level(lev, j) = tau_ref(lev, j)
        transmission_levels(lev, j) = tau_ref(lev, j)
      ENDDO
    ENDIF
  ENDDO

  ! JAH - is this only required for hi-res sensors?
  IF (coef%id_sensor == sensor_id_hi .AND. coef%fmv_model_ver >= 9) THEN
    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        chan = chanprof(j)%chan
        IF (coef%tt_val_chn(chan) == 1) THEN
          DO lev = 1, nlevels
            IF (tau_level(lev, j) < coef%tt_a0(chan)) THEN
              tau_level(lev, j)               = coef%tt_a1(chan)
              transmission_levels(lev, j) = tau_level(lev, j)
            ENDIF
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF

!---------------------------------------------------------------------------------------
!2. Compute optical depth and transmittance at surface
!---------------------------------------------------------------------------------------
!   DO j = 1, nchannels
!     IF (chanflag(j)) THEN
!       prof       = chanprof(j)%prof
! ! as defined in rttov_profaux
!       levsurf    = aux%s(prof)%nearestlev_surf
! ! layer above this
! ! NB all od_level -ve
! ! if surface below nlevels, pfraction -ve
!       IF( (od_level(levsurf - 1, j) - od_level(levsurf, j) > 0._jprb)) THEN
!         od_surf(j) = od_level(levsurf, j) + aux%s(prof)%pfraction_surf * (od_level(levsurf - 1, j) - od_level(levsurf, j))
!         IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
!           od_frac(j) = MIN(od_surf(j) - od_level(levsurf - 1, j), 0._jprb)
!         ELSE
!           od_frac(j) = MIN(od_surf(j) - od_level(levsurf, j), 0._jprb)
!         ENDIF

!       ELSE
!         od_surf(j) = od_level(levsurf, j)
!         IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
!           od_frac(j) = MIN(od_surf(j) - od_level(levsurf - 1, j), 0._jprb)
!         ELSE
!           od_frac(j) = 0._jprb
!         ENDIF
!       ENDIF
!       tau_surf(j)     = EXP(od_surf(j))
!       tau_ref_surf(j) = tau_surf(j)
!     ENDIF
!   ENDDO

  DO j = 1, nchannels
    IF (chanflag(j)) THEN
      prof       = chanprof(j)%prof
! as defined in rttov_profaux
      levsurf    = aux%s(prof)%nearestlev_surf
! layer above this
! NB all od_level -ve
! if surface below nlevels, pfraction -ve
      od_surf(j) = od_level(levsurf, j) + aux%s(prof)%pfraction_surf * (od_level(levsurf - 1, j) - od_level(levsurf, j))
      IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
        od_frac(j) = od_surf(j) - od_level(levsurf - 1, j)
      ELSE
        od_frac(j) = od_surf(j) - od_level(levsurf, j)
      ENDIF
      tau_surf(j)     = EXP(od_surf(j))
      tau_ref_surf(j) = tau_surf(j)
    ENDIF
  ENDDO


  ! JAH - is this only required for hi-res sensors?
  IF (coef%id_sensor == sensor_id_hi .AND. coef%fmv_model_ver >= 9) THEN
    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        chan = chanprof(j)%chan
        IF (coef%tt_val_chn(chan) == 1) THEN
          IF (tau_surf(j) < coef%tt_a0(chan)) THEN
            tau_surf(j) = coef%tt_a1(chan)
          ENDIF
        ENDIF
      ENDIF
    ENDDO
  ENDIF
!---Loop over the streams-----------------------------------------------------------------
  IF (addaerosl .OR. addclouds) THEN
    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        prof = chanprof(j)%prof
        levsurf    = aux%s(prof)%nearestlev_surf
! all arrays here based on layers, not levels
        DO ist = 0, ircld%nstream(prof)
          od_surf_ac(ist,j) = transmission_scatt_ir_stream%opdpac(levsurf, ist, j) + aux%s(prof)%pfraction_surf * (     &
            & transmission_scatt_ir_stream%opdpac(levsurf - 1, ist, j) - &
            & transmission_scatt_ir_stream%opdpac(levsurf, ist, j))
          IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
            od_frac_ac = od_surf_ac(ist,j) - transmission_scatt_ir_stream%opdpac(levsurf - 1, ist, j)
          ELSE
            od_frac_ac = od_surf_ac(ist,j) - transmission_scatt_ir_stream%opdpac(levsurf, ist, j)
          ENDIF
          transmission_aux_path%tau_surf_ac(ist, j) = EXP(-od_surf_ac(ist,j))
          transmission_aux_path%od_sfrac(ist, j)    = - od_frac(j) + od_frac_ac
          IF (tau_surf(j) >= 0) THEN
            transmission_aux_path%tau_surf(ist, j) = tau_surf(j) * &
              & transmission_aux_path%tau_surf_ac(ist, j)
          ELSE
            transmission_aux_path%tau_surf(ist, j) = tau_surf(j)
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  IF (ANY(do_lambertian)) THEN
    DO j = 1, nchannels
      IF (chanflag(j) .AND. do_lambertian(j)) THEN
        prof = chanprof(j)%prof
        ref_power = sec_theta_eff * geometry(prof)%coszen
        ! Recall that od_level and od_surf are negative, but the aer/cld opdeps are positive
        IF (addclouds .OR. addaerosl) THEN
          DO ist = 0, ircld%nstream(prof)
            WHERE (transmission_scatt_ir_stream%opdpac(:,ist,j) < max_optical_depth)
              opdpac = transmission_scatt_ir_stream%opdpac(:,ist,j)
            ELSEWHERE
              opdpac = max_optical_depth
            ENDWHERE
            CALL exponential((od_level(:,j) - opdpac) * ref_power, &
                             transmission_aux_path%tau_level_p(:,ist,j))
            opdpsurfac = MIN(od_surf_ac(ist,j), max_optical_depth)
            transmission_aux_path%tau_surf_p(ist,j) = EXP((od_surf(j) - opdpsurfac) * ref_power)
          ENDDO
        ELSE
          CALL exponential(od_level(:,j) * ref_power, transmission_aux_path%tau_level_p(:,0,j))
          transmission_aux_path%tau_surf_p(0,j) = EXP(od_surf(j) * ref_power)
        ENDIF
        DO ist = 0, ircld%nstream(prof)
          CALL reciprocal(transmission_aux_path%tau_level_p(:,ist,j), transmission_aux_path%tau_level_p_r(:,ist,j))
          transmission_aux_path%tau_surf_p_r(ist,j) = 1._jprb / transmission_aux_path%tau_surf_p(ist,j)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
!-----------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!3. Store transmittances for other streams and single stream for o/p
!---------------------------------------------------------------------------------------
  IF (addaerosl .OR. addclouds) THEN
    DO j = 1, nchannels
      IF (chanflag(j)) transmission_total(j) = transmission_aux_path%tau_surf(0, j)
    ENDDO

    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        all_tau_level_ge_zero = .NOT. ANY((tau_level(:,j) < 0._jprb))
        prof = chanprof(j)%prof

        DO ist = 0, ircld%nstream(prof)
          CALL exponential(-transmission_scatt_ir_stream%opdpac(:, ist, j), &
            ztemp(:,ist))
        ENDDO

        WHERE (ztemp(:,0:ircld%nstream(prof)) < EXP(-max_optical_depth))
          ztemp(:,0:ircld%nstream(prof)) = EXP(-max_optical_depth)
        ENDWHERE

        IF (all_tau_level_ge_zero) THEN
          DO ist = 0, ircld%nstream(prof)
            transmission_aux_path%tau_level(:, ist, j) = &
              tau_level(:, j) * ztemp(:, ist)
          ENDDO
        ELSE
          DO lev = 1, nlevels
            IF (tau_level(lev, j) >= 0._jprb) THEN
              transmission_aux_path%tau_level(lev, 0:ircld%nstream(prof), j) = &
                tau_level(lev, j) * ztemp(lev, 0:ircld%nstream(prof))
            ELSE
              transmission_aux_path%tau_level(lev, 0:ircld%nstream(prof), j) = &
                tau_level(lev, j)
            ENDIF
          ENDDO
        ENDIF

        transmission_aux_path%od_singlelayer(:,:,j) = &
          transmission_aux_path%od_singlelayer(:,:,j) + transmission_scatt_ir%opdpacl(:,:,j)

        WHERE (transmission_aux_path%od_singlelayer(:,:,j) < small_val)
          transmission_aux_path%od_singlelayer(:,:,j) = small_val
        ENDWHERE

        transmission_aux_path%tau_surf(0:ircld%nstream(prof), j) = &
          MAX(small_val, transmission_aux_path%tau_surf(0:ircld%nstream(prof), j))

        CALL reciprocal(transmission_aux_path%tau_level(:,0:ircld%nstream(prof),j), &
                        transmission_aux_path%tau_level_r(:,0:ircld%nstream(prof),j))
        CALL reciprocal(transmission_aux_path%od_singlelayer(:,:,j), &
                        transmission_aux_path%od_singlelayer_r(:,:,j))
      ENDIF
    ENDDO

    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        DO lev = 1, nlevels
          IF (tau_level(lev, j) >= 0) THEN
            transmission_levels(lev, j) = transmission_aux_path%tau_level(lev, 0, j)
          ENDIF
        ENDDO
      ENDIF
    ENDDO

  ELSE ! no addclouds or add aerosol
    DO j = 1, nchannels
      IF (chanflag(j)) transmission_total(j) = tau_surf(j)
    ENDDO

    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        prof = chanprof(j)%prof
        ist = 0
        transmission_aux_path%od_sfrac(ist, j) = -od_frac(j)

! DAR stop changes in emissivity_ad/k for v. small tau and stop nag from complaining for under/overflows
        transmission_aux_path%tau_surf(ist, j) = MAX(small_val, tau_surf(j))

        transmission_aux_path%tau_level(:, ist, j) = tau_level(:, j)
        CALL reciprocal(transmission_aux_path%tau_level(:, ist, j), &
          transmission_aux_path%tau_level_r(:, ist, j))
        DO lay = 1, nlayers
          transmission_aux_path%od_singlelayer(0, lay, j) = &
              MAX(small_val, transmission_aux_path%od_singlelayer(0, lay, j))
        ENDDO

        CALL reciprocal(transmission_aux_path%od_singlelayer(0, :, j), &
          transmission_aux_path%od_singlelayer_r(0, :, j))
      ENDIF
    ENDDO
  ENDIF

! initi to zero means this is not allowed:
!transmission_aux_path%tau_surf_r = min(1e100_jprb, 1.0_jprb / transmission_aux_path%tau_surf)

  DO i = 1, nchannels
    IF (chanflag(i)) THEN
      prof = chanprof(i)%prof
      DO ist = 0, ircld%nstream(prof)
        transmission_aux_path%tau_surf_r(ist,i) = 1.0_jprb / transmission_aux_path%tau_surf(ist,i)
! stop floating overflow
        transmission_aux_path%od_sfrac_r(ist,i) = &
          1.0_jprb / (transmission_aux_path%od_sfrac(ist,i) + small_val)
      ENDDO
    ENDIF
  ENDDO
! Moved all of this code from rttov_integrate because it isn't to do with integrating!
! Because it's common in TL/AD/K code it's saved for later use.

  transmission_aux%anynegtau = -0.5_jprb

do i = 1, nchannels
  IF (chanflag(i)) THEN
    prof = chanprof(i)%prof
    do ist = 0, ircld%nstream(prof)
      do lay = 1, nlayers
        isti = ircld%icldarr(ist,lay,prof)
#ifdef RTTOV_XLF 
! IBM codepath uses fsel intrinsic that doesn't generated branched code

! DAR: if anynegtau > 0 then select a slower code path that does more tests (because something has gone wrong earlier?)
        transmission_aux%anynegtau = &
          transmission_aux%anynegtau + fsel(transmission_aux_path%tau_level(lay,ist,i), 0.0_jprb, 1.0_jprb)
! DAR: tau_level differences smaller than min optical depth. Is this condition necessary?
        transmission_aux%fac1(lay,ist,i) = fsel(transmission_aux_path%tau_level(lay, ist, i) - &
                transmission_aux_path%tau_level(lay+1, ist, i) - min_od, &
                1.0_jprb, 0.0_jprb)
! DAR: layer ods smaller than min optical depth. This should be done (is done?) in rttov_transmit as well.
        transmission_aux%fac1(lay,ist,i) = fsel(transmission_aux_path%od_singlelayer(isti,lay,i) - &
                min_od, transmission_aux%fac1(lay,ist,i), 0.0_jprb)

! DAR: is tau_level less than min_tau? THIS should probably be done elsewhere (rttov_transmit?)
!      is it a problem if tau_level is smaller than min_tau? Can we let this underflow gracefully? It's not being used as a divsor
!      consider removing this condition.
        transmission_aux_path%fac2(lay,ist,i) = &
          fsel(transmission_aux_path%tau_level(lay,ist,i) - min_tau, 1.0_jprb, 0.0_jprb)
#else
! Non-IBM codepath         
        lev = lay + 1

        if(transmission_aux_path%tau_level(lay,ist,i) < 0.0_jprb) then
          transmission_aux%anynegtau = 1.0_jprb
        endif

  ! DAR: separated fac1 and fac2 becuase fac2 can be used separately from fac1 in Phil Watts calcs

        if(transmission_aux_path%od_singlelayer(isti,lay,i) < min_od) then
          transmission_aux%fac1(lay,ist,i) = 0.0_jprb
        else
          if(transmission_aux_path%tau_level(lay, ist, i) - &
            transmission_aux_path%tau_level(lev, ist, i) < min_od) then
              transmission_aux%fac1(lay,ist,i) = 0.0_jprb
          else
              transmission_aux%fac1(lay,ist,i) = 1.0_jprb
          endif
        endif

        if(transmission_aux_path%tau_level(lay,ist,i) < min_tau) then
          transmission_aux_path%fac2(lay,ist,i) = 0.0_jprb
        else
          transmission_aux_path%fac2(lay,ist,i) = 1.0_jprb
        endif
#endif  
      enddo
      if(transmission_aux_path%tau_level(nlevels,ist,i) < min_tau) then
          transmission_aux_path%fac2(nlevels,ist,i) = 0.0_jprb
      else
          transmission_aux_path%fac2(nlevels,ist,i) = 1.0_jprb
      endif
        
    enddo
  ENDIF
enddo

do i = 1, nchannels
  IF (chanflag(i)) THEN
    prof = chanprof(i)%prof
    do ist = 0, ircld%nstream(prof)
#ifdef RTTOV_XLF
      transmission_aux%surf_fac(ist,i) = &
        & fsel(transmission_aux_path%tau_surf(ist, i) - min_tau, 1.0_jprb, 0.0_jprb)
#else
! Non-IBM codepath         
      if(transmission_aux_path%tau_surf(ist, i) > min_tau) then
        transmission_aux%surf_fac(ist,i) = 1.0_jprb
      else
        transmission_aux%surf_fac(ist,i) = 0.0_jprb
      endif
#endif
    enddo
  ENDIF
enddo

IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_transmit
