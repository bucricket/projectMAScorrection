SUBROUTINE rttov_transmit_ad( &
            & addaerosl,                       &
            & addclouds,                       &
            & do_lambertian,                   &
            & nlayers,                         &
            & chanprof,                        &
            & chanflag,                        &
            & aux,                             &
            & aux_ad,                          &
            & coef,                            &
            & ircld,                           &
            & geometry,                        &
            & opdp_path,                       &
            & opdp_path_ad,                    &
            & od_level,                        &
            & transmission_levels,             &
            & transmission_total,              &
            & transmission_levels_ad,          &
            & transmission_total_ad,           &
            & transmission_aux,                &
            & transmission_aux_path,           &
            & transmission_aux_path_ad,        &
            & transmission_scatt_ir,           &
            & transmission_scatt_ir_ad,        &
            & transmission_scatt_ir_stream,    &
            & transmission_scatt_ir_stream_ad, &
            & tau_ref,                         &
            & tau_ref_surf,                    &
            & tau_surf,                        &
            & tau_level)
! Description:
! Adjoint of rttov_transmit_tl
!  To calculate optical depths for a number of channels
!  and profiles from every pressure level to space.
! To interpolate optical depths on to levels of radiative transfer model
! (which, at present, entails only surface transmittance, as
! other optical depths are on *rt* levels) and to convert
! optical depths to transmittances.
!
! Code based on OPDEPAD and RTTAUAD from previous versions of RTTOV
! Only one profile per call
!
! Adjoint variables
! input transmission_aux_path_ad % tau_surf and transmission_aux_path_ad % tau_level
! set inside integrate_ad
!
! input/output aux_ad, aux_path_ad
!
! output predictors_ad initialised inside rttov_ad (need input
!    intent for memory allocation in calling routine)
!
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
!            --       New routine based on rttov_transmit_ad.F90.
!                     Variable trace gases, CO2, N2O,CO and CH4
!                     have been introduced for AIRS and IASI
!  1.1    06/02/2007  Removed polarisation R Saunders
!  1.2    11/12/2007  Input also adjoint of tau_total and tau_levels since
!                     these are the external outputs of RTTOV (A Geer)
!  1.3    04/06/2008  Fix od_frac and od_frac_ac calculation near surface level (PB PM)
!  1.4    15/08/2009  User defined ToA. Layers distinct from levels (P.Rayer)
!  1.5    03/11/2009  Transmittances / optical depths on levels (A Geer)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud

  USE rttov_types, ONLY :  &
       & rttov_chanprof,             &
       & rttov_coef,                 &
       & transmission_type_aux,      &
       & rttov_path_transmission,    &
       & transmission_scatt_ir_type, &
       & profile_aux,                &
       & ircld_type,                 &
       geometry_type
  USE parkind1, ONLY : jpim, jprb, jplm
  USE rttov_math_mod
!INTF_OFF
  USE rttov_const, ONLY : sensor_id_hi, sec_theta_eff
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE

  LOGICAL(KIND=jplm)              , INTENT(IN)    :: addaerosl
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: addclouds
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: do_lambertian(:)
  INTEGER(KIND=jpim)              , INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof            ), INTENT(IN)    :: chanprof(:)
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: chanflag(SIZE(chanprof))
  TYPE(rttov_coef                ), INTENT(IN)    :: coef
  TYPE(profile_aux               ), INTENT(IN)    :: aux
  TYPE(profile_aux               ), INTENT(INOUT) :: aux_ad
  TYPE(ircld_type                ), INTENT(IN)    :: ircld
  TYPE(geometry_type),              INTENT(in)    :: geometry(:)
  REAL(KIND=jprb)                 , INTENT(IN)    :: opdp_path(:,:)
  REAL(KIND=jprb)                 , INTENT(INOUT) :: opdp_path_ad(:,:)
  REAL(KIND=jprb)                 , INTENT(IN)    :: transmission_levels(:,:)
  REAL(KIND=jprb)                 , INTENT(IN)    :: transmission_total(:)
  REAL(KIND=jprb)                 , INTENT(INOUT) :: transmission_levels_ad(:,:)
  REAL(KIND=jprb)                 , INTENT(INOUT) :: transmission_total_ad(:)
  TYPE(transmission_type_aux     ), INTENT(IN)    :: transmission_aux
  TYPE(rttov_path_transmission   ), INTENT(IN)    :: transmission_aux_path
  TYPE(rttov_path_transmission   ), INTENT(INOUT) :: transmission_aux_path_ad
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_ad
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir_stream
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream_ad
  REAL(KIND=jprb)                 , INTENT(IN)    :: tau_ref(nlayers + 1, SIZE(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: tau_ref_surf(SIZE(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: od_level(nlayers + 1, SIZE(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: tau_surf(SIZE(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: tau_level(nlayers + 1, SIZE(chanprof))
!INTF_END

  REAL   (KIND=jprb) :: ref_power 
!   REAL   (KIND=jprb) :: od_surf(SIZE(chanprof))
  REAL   (KIND=jprb) :: od_level_ad(nlayers + 1, SIZE(chanprof))
  REAL   (KIND=jprb) :: od_surf_ad(SIZE(chanprof))
  REAL   (KIND=jprb) :: od_surf_ac_ad(0:SIZE(transmission_aux_path%tau_level(1,:,1)), SIZE(chanprof))
  REAL   (KIND=jprb) :: od_frac_ad(SIZE(chanprof)), od_frac_ac_ad
  REAL   (KIND=jprb) :: tau_surf_ad(SIZE(chanprof))
  REAL   (KIND=jprb) :: tau_level_ad(nlayers + 1, SIZE(chanprof))

  INTEGER(KIND=jpim) :: lev, lay, chan, j, levsurf, ist
  INTEGER(KIND=jpim) :: prof, nlevels
  INTEGER(KIND=jpim) :: nchannels

! allocate automatic arrays using SIZE(transmission_aux_path%tau_level(1,:,1) to keep pgf90 compiler happy
  REAL   (KIND=jprb) :: ztemp(nlayers+1, 0:SIZE(transmission_aux_path%tau_level(1,:,1)))
  REAL   (KIND=jprb) :: ztemp_ad(nlayers+1, 0:SIZE(transmission_aux_path%tau_level(1,:,1)))
  LOGICAL(kind=jplm) :: all_tau_level_ge_zero

  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_AD', 0_jpim, ZHOOK_HANDLE)
  nchannels              = SIZE(chanprof)
  nlevels                = nlayers + 1
!-----------------------------------------------------
!AD of store transmittances for other streams
!-----------------------------------------------------
  od_level_ad(:,:)       = 0._JPRB
  tau_level_ad(:,:)      = 0.0_JPRB
  tau_surf_ad(:)         = 0.0_JPRB
  od_frac_ad(:)          = 0.0_JPRB
  od_surf_ad(:)          = 0.0_JPRB
  od_surf_ac_ad          = 0.0_JPRB
  od_frac_ac_ad          = 0.0_JPRB
  
  IF (addaerosl .OR. addclouds) THEN
    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        all_tau_level_ge_zero = .NOT. ANY((tau_level(:,j) < 0._jprb))
        prof = chanprof(j)%prof

        DO ist = 0, ircld%nstream(prof)
          CALL exponential(-transmission_scatt_ir_stream%opdpac(:, ist, j), &
            ztemp(:,ist))
        ENDDO
!start AD

        transmission_scatt_ir_ad%opdpacl(:,:,j) = &
            transmission_scatt_ir_ad%opdpacl(:,:,j) + &
            transmission_aux_path_ad%od_singlelayer(:,:,j)

        IF (all_tau_level_ge_zero) THEN
          DO ist = 0, ircld%nstream(prof)
            tau_level_ad(:, j) = tau_level_ad(:, j) + &
              transmission_aux_path_ad%tau_level(:, ist, j) * ztemp(:, ist)

            ztemp_ad(:, ist) = &!ztemp_ad(:, ist) + &
              tau_level(:, j) * transmission_aux_path_ad%tau_level(:, ist, j)
          ENDDO
        ELSE
          DO lev = 1, nlevels
            IF (tau_level(lev, j) >= 0._jprb) THEN
              tau_level_ad(lev, j) = tau_level_ad(lev, j) + &
                SUM(transmission_aux_path_ad%tau_level(lev, 0:ircld%nstream(prof), j) * &
                    ztemp(lev, 0:ircld%nstream(prof)))

              ztemp_ad(lev, 0:ircld%nstream(prof)) = &!ztemp_ad(lev, 0:ircld%nstream(prof)) + &
                tau_level(lev, j) * &
                transmission_aux_path_ad%tau_level(lev, 0:ircld%nstream(prof), j)
            ELSE
              tau_level_ad(lev, j) = tau_level_ad(lev, j) + & 
                SUM(transmission_aux_path_ad%tau_level(lev, 0:ircld%nstream(prof), j))
            ENDIF
          ENDDO
        ENDIF

        DO ist = 0, ircld%nstream(prof)
          CALL exponential_ad(-ztemp(:,ist) , &
                             transmission_scatt_ir_stream_ad%opdpac(:, ist, j), &
                             ztemp_ad(:,ist), acc = .TRUE._jplm)
        ENDDO

      ENDIF
    ENDDO
  ELSE
    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        prof = chanprof(j)%prof
        ist = 0
        DO lev = nlevels, 1,  - 1
          tau_level_ad(lev, j) = tau_level_ad(lev, j) + &
            transmission_aux_path_ad%tau_level(lev, ist, j)
!            transmission_aux_path_ad%tau_level(lev, ist, j) = 0.0_JPRB
        ENDDO

        od_frac_ad(j) = od_frac_ad(j) - transmission_aux_path_ad%od_sfrac(ist, j)
!          transmission_aux_path_ad%od_sfrac(ist, j) = 0.0_JPRB
        tau_surf_ad(j) = tau_surf_ad(j) + transmission_aux_path_ad%tau_surf(ist, j)
!          transmission_aux_path_ad%tau_surf(ist, j) = 0.0_JPRB
        tau_level_ad(:, j) = tau_level_ad(:, j) + transmission_levels_ad(:, j)
        transmission_levels_ad(:, j) = 0
        tau_surf_ad(j) = tau_surf_ad(j) + transmission_total_ad(j)
        transmission_total_ad(j) = 0
      ENDIF
    ENDDO
  ENDIF

  IF (ANY(do_lambertian)) THEN
    DO j = 1, nchannels
      IF (chanflag(j) .AND. do_lambertian(j)) THEN
        prof = chanprof(j)%prof
        ref_power = sec_theta_eff * geometry(prof)%coszen
        DO ist = 0, ircld%nstream(prof)
          transmission_aux_path_ad%tau_surf_p(ist,j) = transmission_aux_path_ad%tau_surf_p(ist,j) - &
            transmission_aux_path_ad%tau_surf_p_r(ist,j) * transmission_aux_path%tau_surf_p_r(ist,j) ** 2_jpim

          CALL reciprocal_ad(transmission_aux_path%tau_level_p_r(:,ist,j), &
                             transmission_aux_path_ad%tau_level_p(:,ist,j), &
                             transmission_aux_path_ad%tau_level_p_r(:,ist,j), acc=.TRUE._jplm)
        ENDDO
        IF (addclouds .OR. addaerosl) THEN
          DO ist = 0, ircld%nstream(prof)
            od_surf_ad(j) = od_surf_ad(j) + &
              ref_power * transmission_aux_path%tau_surf_p(ist,j) * transmission_aux_path_ad%tau_surf_p(ist,j)

            od_surf_ac_ad(ist,j) = od_surf_ac_ad(ist,j) - &
              ref_power * transmission_aux_path%tau_surf_p(ist,j) * transmission_aux_path_ad%tau_surf_p(ist,j)

            ! must move ref_power to input because can't have it with intent inout
            CALL exponential_ad(transmission_aux_path%tau_level_p(:,ist,j) * ref_power, &
              od_level_ad(:,j), &
              transmission_aux_path_ad%tau_level_p(:,ist,j), acc=.TRUE._jplm)

            CALL exponential_ad(-transmission_aux_path%tau_level_p(:,ist,j) * ref_power, &
              transmission_scatt_ir_stream_ad%opdpac(:,ist,j), &
              transmission_aux_path_ad%tau_level_p(:,ist,j), acc=.TRUE._jplm)

          ENDDO
        ELSE
          od_surf_ad(j) = od_surf_ad(j) + &
            ref_power * transmission_aux_path%tau_surf_p(0,j) * transmission_aux_path_ad%tau_surf_p(0,j)

          ! must move ref_power to input because can't have it with intent inout
          CALL exponential_ad(transmission_aux_path%tau_level_p(:,0,j) * ref_power, &
            od_level_ad(:,j), &
            transmission_aux_path_ad%tau_level_p(:,0,j), acc=.TRUE._jplm)
        ENDIF
      ENDIF
    ENDDO
  ENDIF

!---------------------------------------------------------------------------------------
!AD of compute optical depth and transmittance at surface
!---------------------------------------------------------------------------------------
  IF (addaerosl .OR. addclouds) THEN
    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        prof    = chanprof(j)%prof
        levsurf = aux%s(prof)%nearestlev_surf
!---Loop over the streams-----------------------------------------------------------------
        DO ist = ircld%nstream(prof), 0,  - 1
          IF (tau_surf(j) >= 0) THEN
            transmission_aux_path_ad%tau_surf_ac(ist, j) =     &
              & transmission_aux_path_ad%tau_surf_ac(ist, j) + &
              & transmission_aux_path_ad%tau_surf(ist, j) * tau_surf(j)
            tau_surf_ad(j)                          =      &
              & tau_surf_ad(j) + transmission_aux_path_ad%tau_surf(ist, j) * &
              & transmission_aux_path%tau_surf_ac(ist, j)
          ELSE
            tau_surf_ad(j) = tau_surf_ad(j) + transmission_aux_path_ad%tau_surf(ist, j)
          ENDIF
          od_frac_ad(j) = od_frac_ad(j) - transmission_aux_path_ad%od_sfrac(ist, j)
          od_frac_ac_ad = od_frac_ac_ad + transmission_aux_path_ad%od_sfrac(ist, j)
          od_surf_ac_ad(ist,j) =      &
            & od_surf_ac_ad(ist,j) - transmission_aux_path_ad%tau_surf_ac(ist, j) * &
            & transmission_aux_path%tau_surf_ac(ist, j)
          IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
            od_surf_ac_ad(ist,j) = od_surf_ac_ad(ist,j) + od_frac_ac_ad
            transmission_scatt_ir_stream_ad%opdpac(levsurf - 1, ist, j) =     &
              & transmission_scatt_ir_stream_ad%opdpac(levsurf - 1, ist, j) - &
              & od_frac_ac_ad
          ELSE
            od_surf_ac_ad(ist,j) = od_surf_ac_ad(ist,j) + od_frac_ac_ad
            transmission_scatt_ir_stream_ad%opdpac(levsurf, ist, j) =     &
              & transmission_scatt_ir_stream_ad%opdpac(levsurf, ist, j) - &
              & od_frac_ac_ad
          ENDIF
          transmission_scatt_ir_stream_ad%opdpac(levsurf, ist, j)     =      &
            & transmission_scatt_ir_stream_ad%opdpac(levsurf, ist, j) +      &
            & od_surf_ac_ad(ist,j) * (1._jprb - aux%s(prof)%pfraction_surf)
          transmission_scatt_ir_stream_ad%opdpac(levsurf - 1, ist, j) =      &
            & transmission_scatt_ir_stream_ad%opdpac(levsurf - 1, ist, j) +  &
            & od_surf_ac_ad(ist,j) * aux%s(prof)%pfraction_surf
          aux_ad%s(prof)%pfraction_surf                               =      &
            & aux_ad%s(prof)%pfraction_surf + od_surf_ac_ad(ist,j) * (       &
            & transmission_scatt_ir_stream%opdpac(levsurf - 1, ist, j) -     &
            & transmission_scatt_ir_stream%opdpac(levsurf, ist, j))
!           od_surf_ac_ad(ist,j) = 0._jprb
          od_frac_ac_ad = 0._jprb
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  IF (coef%id_sensor == sensor_id_hi .AND. coef%fmv_model_ver >= 9) THEN
    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        chan = chanprof(j)%chan
        IF (coef%tt_val_chn(chan) == 1) THEN
          IF (tau_ref_surf(j) < coef%tt_a0(chan)) THEN
            tau_surf_ad(j) = 0._jprb
          ENDIF
        ENDIF
      ENDIF
    ENDDO
  ENDIF
  ! DO j = 1, nchannels
  !   prof          = chanprof(j)%prof
  !   levsurf       = aux%s(prof)%nearestlev_surf
  !   od_surf_ad(j) = od_surf_ad(j) + tau_surf_ad(j) * tau_ref_surf(j)

  !   IF( (od_level(levsurf - 1, j) - od_level(levsurf, j) > 0._jprb)) THEN
  !     od_surf(j) = od_level(levsurf, j) + aux%s(prof)%pfraction_surf * (od_level(levsurf - 1, j) - od_level(levsurf, j)) !FWD CODE
  !     IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
  !       IF(od_surf(j) - od_level(levsurf - 1, j) < 0._jprb) THEN
  !         od_surf_ad(j) = od_surf_ad(j) + od_frac_ad(j)
  !         od_level_ad(levsurf - 1,j) = od_level_ad(levsurf - 1,j) - od_frac_ad(j)
  !       ELSE
  !         od_frac_ad(j) = 0._jprb
  !       ENDIF
  !     ELSE
  !       IF(od_surf(j) - od_level(levsurf, j) < 0._jprb) THEN
  !         od_surf_ad(j) = od_surf_ad(j) + od_frac_ad(j)
  !         od_level_ad(levsurf,j) = od_level_ad(levsurf,j) - od_frac_ad(j)
  !       ELSE
  !         od_frac_ad(j) = 0._jprb
  !       ENDIF
  !     ENDIF
      
  !     od_level_ad(levsurf, j) = od_level_ad(levsurf, j) + od_surf_ad(j) * (1._jprb - aux%s(prof)%pfraction_surf)
  !     od_level_ad(levsurf-1, j) = od_level_ad(levsurf-1, j) + od_surf_ad(j) * aux%s(prof)%pfraction_surf
  !     aux_ad%s(prof)%pfraction_surf = aux_ad%s(prof)%pfraction_surf + &
  !       od_surf_ad(j) * (od_level(levsurf - 1, j) - od_level(levsurf, j))
  !   ELSE
  !     od_surf(j) = od_level(levsurf, j) !FWD CODE
  !     IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
  !       IF(od_surf(j) - od_level(levsurf - 1, j) < 0._jprb) THEN
  !         od_surf_ad(j) = od_surf_ad(j) + od_frac_ad(j)
  !         od_level_ad(levsurf - 1,j) = od_level_ad(levsurf - 1,j) - od_frac_ad(j)
  !       ELSE
  !         od_frac_ad(j) = 0._jprb
  !       ENDIF
  !     ELSE
  !       od_frac_ad(j) = 0._jprb
  !     ENDIF
  !     od_level_ad(levsurf, j) = od_level_ad(levsurf, j) + od_surf_ad(j)
  !   ENDIF
  !   od_surf_ad(j) = 0._jprb
  ! ENDDO

  DO j = 1, nchannels
    IF (chanflag(j)) THEN
      prof          = chanprof(j)%prof
      levsurf       = aux%s(prof)%nearestlev_surf
      od_surf_ad(j) = od_surf_ad(j) + tau_surf_ad(j) * tau_ref_surf(j)
      IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
        od_surf_ad(j)               = od_surf_ad(j) + od_frac_ad(j)
        od_level_ad(levsurf - 1, j) = od_level_ad(levsurf - 1, j) - od_frac_ad(j)
      ELSE
        od_surf_ad(j)           = od_surf_ad(j) + od_frac_ad(j)
        od_level_ad(levsurf, j) = od_level_ad(levsurf, j) - od_frac_ad(j)
      ENDIF
      od_level_ad(levsurf, j)       = od_level_ad(levsurf, j) + od_surf_ad(j) * (1._JPRB - aux%s(prof)%pfraction_surf)
      od_level_ad(levsurf - 1, j)   = od_level_ad(levsurf - 1, j) + od_surf_ad(j) * aux%s(prof)%pfraction_surf
      aux_ad%s(prof)%pfraction_surf =      &
        & aux_ad%s(prof)%pfraction_surf + od_surf_ad(j) * (od_level(levsurf - 1, j) - od_level(levsurf, j))
      od_surf_ad(j)                 = 0._JPRB
    ENDIF
  ENDDO
!-------------------------------------------
!AD of assemble layer optical depths
!-------------------------------------------
  IF (coef%id_sensor == sensor_id_hi .AND. coef%fmv_model_ver >= 9) THEN
    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        chan = chanprof(j)%chan
        IF (coef%tt_val_chn(chan) == 1) THEN
          DO lev = 1, nlevels
            IF (tau_ref(lev, j) < coef%tt_a0(chan)) THEN
              tau_level_ad(lev, j) = 0._JPRB
            ENDIF
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF
! ad of tau_level_tl
  DO j = 1, nchannels
    IF (chanflag(j)) THEN
      chan = chanprof(j)%chan
      DO lev = 1, nlevels
        od_level_ad(lev, j)  = od_level_ad(lev, j) + tau_level_ad(lev, j) * tau_ref(lev, j)
        od_level_ad(lev, j)  = coef%ff_gam(chan) * od_level_ad(lev, j)
        opdp_path_ad(lev, j) = opdp_path_ad(lev, j) + od_level_ad(lev, j)
      ENDDO
      transmission_aux_path_ad%od_singlelayer(:,:,j) = &
        coef%ff_gam(chan) * transmission_aux_path_ad%od_singlelayer(:,:,j)
    ENDIF
  ENDDO
! ad of level to space optical depths
! opdp_path_ad(:,:) = opdp_path_ad(:,:) + od_level_ad(:,:)
! ad of single layer optical depths
  DO j = 1, nchannels
    IF (chanflag(j)) THEN
      DO lay = nlayers, 1,  - 1
        opdp_path_ad(lay,j)   = &
          opdp_path_ad(lay,j) + SUM(transmission_aux_path_ad%od_singlelayer(:,lay,j))
        opdp_path_ad(lay+1,j) = &
          opdp_path_ad(lay+1,j) - SUM(transmission_aux_path_ad%od_singlelayer(:,lay,j))
      ENDDO
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_transmit_ad
