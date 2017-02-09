SUBROUTINE rttov_transmit_9_solar_ad( &
            & addaerosl,                       &
            & addclouds,                       &
            & nlayers,                         &
            & nprofiles,                       &
            & chanprof,                        &
            & solar,                           &
            & aux,                             &
            & aux_ad,                          &
            & coef,                            &
            & raytracing,                      &
            & raytracing_ad,                   &
            & ircld,                           &
            & opdp_path,                       &
            & opdp_path_ad,                    &
            & path2,                           &
            & path1,                           &
            & transmission,                    &
            & transmission_ad,                 &
            & transmission_aux,                &
            & transmission_aux_ad,             &
            & transmission_scatt_ir,           &
            & transmission_scatt_ir_ad,        &
            & transmission_scatt_ir_stream,    &
            & transmission_scatt_ir_stream_ad)
!
! Description:
! Adjoint of rttov_transmit_ad
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
! input transmission_aux_ad % solar_path2 % tau_surf and transmission_aux_ad % solar_path2 % tau_level
! set inside integrate_ad
!
! input/output aux_ad
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
!  1.1    06/02/2007  Removed polarisation index R Saunders
!  1.2    15/08/2009  User defined ToA. Layers distinct from levels (P.Rayer)
!  1.3    03/11/2009  Transmittances / optical depths on levels (A Geer)
!  1.4    02/12/2009  Pathsat, Pathsun and related quantities are now
!                     layer arrays (Marco Matricardi).
!  1.5    05/07/2010  Remove addsolar flag from profiles structure (J Hocking)
!  1.6    04/08/2010  Move addsolar check to calling routine (J Hocking)
!  1.7    14/12/2010  Use traj0_sta%solar array to flag channels for which solar calculations
!                     should be performed (J Hocking)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
! Declarations:
! Modules used:
! Imported Parameters:
  USE rttov_types, ONLY :  &
       & rttov_chanprof,             &
       & rttov_coef,                 &
       & opdp_path_type,             &
       & rttov_path_traj_trans,      &
       & transmission_type,          &
       & transmission_type_aux,      &
       & transmission_scatt_ir_type, &
       & profile_aux,                &
       & ircld_type,                 &
       & raytracing_type
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE parkind1, ONLY : jprb
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE rttov_const, ONLY : sensor_id_hi
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  LOGICAL(KIND=jplm), INTENT(IN)    :: addaerosl
  LOGICAL(KIND=jplm), INTENT(IN)    :: addclouds
  INTEGER(KIND=jpim)              , INTENT(IN)    :: nlayers
  INTEGER(KIND=jpim)              , INTENT(IN)    :: nprofiles
  TYPE(rttov_chanprof            ), INTENT(IN)    :: chanprof(:)
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: solar(SIZE(chanprof))
  TYPE(rttov_coef                ), INTENT(IN)    :: coef
  TYPE(profile_aux               ), INTENT(IN)    :: aux
  TYPE(profile_aux               ), INTENT(INOUT) :: aux_ad
  TYPE(ircld_type                ), INTENT(IN)    :: ircld
  TYPE(opdp_path_type            ), INTENT(IN)    :: opdp_path
  TYPE(opdp_path_type            ), INTENT(INOUT) :: opdp_path_ad
  TYPE(rttov_path_traj_trans     ), INTENT(IN)    :: path2
  TYPE(rttov_path_traj_trans     ), INTENT(IN)    :: path1
  TYPE(transmission_type         ), INTENT(IN)    :: transmission
  TYPE(transmission_type         ), INTENT(INOUT) :: transmission_ad
  TYPE(transmission_type_aux     ), INTENT(IN)    :: transmission_aux
  TYPE(transmission_type_aux     ), INTENT(INOUT) :: transmission_aux_ad
  TYPE(raytracing_type           ), INTENT(IN)    :: raytracing
  TYPE(raytracing_type           ), INTENT(INOUT) :: raytracing_ad
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_ad
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir_stream
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream_ad
!INTF_END
!local variables:
  REAL   (KIND=jprb) :: od_singlelayer_ad(nlayers,SIZE(chanprof))
  REAL   (KIND=jprb) :: od_surf_ac_ad(SIZE(chanprof))
  REAL   (KIND=jprb) :: od_level_ad(nlayers + 1, SIZE(chanprof))
  REAL   (KIND=jprb) :: od_surf_ad(SIZE(chanprof))
  REAL   (KIND=jprb) :: od_frac_ad(SIZE(chanprof))
  REAL   (KIND=jprb) :: tausun_level_ad(nlayers + 1, SIZE(chanprof))
  REAL   (KIND=jprb) :: tausun_surf_ad(SIZE(chanprof))
  REAL   (KIND=jprb) :: tau_level_ad(nlayers + 1, SIZE(chanprof))
  REAL   (KIND=jprb) :: tau_surf_ad(SIZE(chanprof))
  REAL   (KIND=jprb) :: pathsat_patheff_ad(nlayers, nprofiles)
  REAL   (KIND=jprb) :: pathsun_patheff_ad(nlayers, nprofiles)
  REAL   (KIND=jprb) :: od_surf, od_surf_ac, zpatheff, tauacsunpath1
  INTEGER(KIND=jpim) :: lev, lay, chan, j, prof, ist, isti, levsurf, nlevels
  INTEGER(KIND=jpim) :: nchannels
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_9_SOLAR_AD', 0_jpim, ZHOOK_HANDLE)
  nchannels              = SIZE(chanprof)
  nlevels                = nlayers + 1
!---------------------------------------------------------------------------------------
!AD of store transmittances for other polarisations
!---------------------------------------------------------------------------------------
  od_level_ad(:,:)       = 0._JPRB
  tausun_level_ad(:,:)   = 0.0_JPRB
  tau_level_ad(:,:)      = 0.0_JPRB
  od_singlelayer_ad(:,:) = 0.0_JPRB
  tausun_surf_ad(:)      = 0.0_JPRB
  tau_surf_ad(:)         = 0.0_JPRB
  od_frac_ad(:)          = 0.0_JPRB
  od_surf_ad(:)          = 0.0_JPRB
  od_surf_ac_ad(:)       = 0.0_JPRB
  pathsat_patheff_ad(:,:) = 0.0_JPRB
  pathsun_patheff_ad(:,:) = 0.0_JPRB
  DO j = nchannels, 1,  - 1
    prof = chanprof(j)%prof
    chan = chanprof(j)%chan
    IF (solar(j)) THEN

      IF (addaerosl .OR. addclouds) THEN
        DO lay = nlayers, 1, -1
          WHERE (transmission_scatt_ir%opdpext(:, lay, j) /= 0._jprb)
            transmission_scatt_ir_ad%opdpsca(:, lay, j) =      &
              & transmission_scatt_ir_ad%opdpsca(:, lay, j) +  &
              & transmission_scatt_ir_ad%ssa(:, lay, j) /      &
              & transmission_scatt_ir%opdpext(:, lay, j)
            transmission_scatt_ir_ad%opdpext(:, lay, j) =      &
              & transmission_scatt_ir_ad%opdpext(:, lay, j) -  &
              & transmission_scatt_ir_ad%ssa(:, lay, j) *      &
              & transmission_scatt_ir%opdpsca(:, lay, j) /     &
              & transmission_scatt_ir%opdpext(:, lay, j) ** 2
          ENDWHERE
          od_singlelayer_ad(lay, j)                               =  &
            & od_singlelayer_ad(lay, j) + SUM(transmission_scatt_ir_ad%opdpext(:, lay, j))
          transmission_scatt_ir_ad%opdpabs(:, lay, j)    =  &
            & transmission_scatt_ir_ad%opdpabs(:, lay, j) + &
            & transmission_scatt_ir_ad%opdpext(:, lay, j)
          transmission_scatt_ir_ad%opdpsca(:, lay, j)    =  &
            & transmission_scatt_ir_ad%opdpsca(:, lay, j) + &
            & transmission_scatt_ir_ad%opdpext(:, lay, j)
        ENDDO
      ENDIF

      IF (addaerosl .OR. addclouds) THEN
        DO isti = 0, 1
          DO lay = nlayers, 1,  - 1
            lev = lay + 1

            od_singlelayer_ad(lay, j) = od_singlelayer_ad(lay, j) + &
              & transmission_aux_ad%solar_path1%od_singlelayer(isti, lay, j) * &
              & raytracing%pathsat(lay, prof) / raytracing%patheff(lay, prof)
            transmission_scatt_ir_ad%opdpaclsun(isti, lay, j) = &
              & transmission_scatt_ir_ad%opdpaclsun(isti, lay, j) + &
              & transmission_aux_ad%solar_path1%od_singlelayer(isti, lay, j) * &
              & raytracing%pathsat(lay, prof) / raytracing%patheff(lay, prof)
            pathsat_patheff_ad(lay, prof) = pathsat_patheff_ad(lay, prof) + &
              & transmission_aux_ad%solar_path1%od_singlelayer(isti, lay, j) * &
              & (path2%od_singlelayer(lay, j) + transmission_scatt_ir%opdpaclsun(isti, lay, j))

            od_singlelayer_ad(lay, j) = od_singlelayer_ad(lay, j) + &
              & transmission_aux_ad%solar_path2%od_singlelayer(isti, lay, j) * &
              & raytracing%pathsun(lay, prof) / raytracing%patheff(lay, prof)
            transmission_scatt_ir_ad%opdpaclsun(isti, lay, j) = &
              & transmission_scatt_ir_ad%opdpaclsun(isti, lay, j) + &
              & transmission_aux_ad%solar_path2%od_singlelayer(isti, lay, j) * &
              & raytracing%pathsun(lay, prof) / raytracing%patheff(lay, prof)
            pathsun_patheff_ad(lay, prof) = pathsun_patheff_ad(lay, prof) + &
              & transmission_aux_ad%solar_path2%od_singlelayer(isti, lay, j) * &
              & (path2%od_singlelayer(lay, j) + transmission_scatt_ir%opdpaclsun(isti, lay, j))
          ENDDO
          IF (.NOT. addclouds) EXIT
        ENDDO
        DO ist = ircld%nstream(prof), 0,  - 1
          DO lev = nlevels, 1,  - 1
            IF (path1%tau_level(lev, j) >= 0) THEN
              IF (lev > 1) THEN
                lay = lev - 1

                tauacsunpath1 = EXP( - transmission_scatt_ir_stream%opdpacsun(lev, ist, j) * &
                  & raytracing%pathsat(lay, prof) / raytracing%patheff(lay, prof))

                tau_level_ad(lev, j) = tau_level_ad(lev, j) +      &
                  & tauacsunpath1 * transmission_aux_ad%solar_path1%tau_level(lev, ist, j)

                transmission_scatt_ir_stream_ad%opdpacsun(lev, ist, j) = &
                  & transmission_scatt_ir_stream_ad%opdpacsun(lev, ist, j) - &
                  & tauacsunpath1 * path1%tau_level(lev, j) * &
                  & transmission_aux_ad%solar_path1%tau_level(lev, ist, j) * &
                  & raytracing%pathsat(lay, prof) / raytracing%patheff(lay, prof)

                pathsat_patheff_ad(lay, prof) = pathsat_patheff_ad(lay, prof) - &
                  & tauacsunpath1 * path1%tau_level(lev, j) * &
                  & transmission_aux_ad%solar_path1%tau_level(lev, ist, j) * &
                  & transmission_scatt_ir_stream%opdpacsun(lev, ist, j)
              ELSE
                ! for lev == 1, opdpacsun(lev, ist, j) == 0.
                tau_level_ad(lev, j) = tau_level_ad(lev, j) + &
                  & transmission_aux_ad%solar_path1%tau_level(lev, ist, j)
              ENDIF
            ELSE
              tau_level_ad(lev, j) = tau_level_ad(lev, j) + &
                & transmission_aux_ad%solar_path1%tau_level(lev, ist, j)
            ENDIF

            IF (path2%tau_level(lev, j) >= 0) THEN
              tausun_level_ad(lev, j)                                = tausun_level_ad(lev, j) +      &
                & transmission_aux_ad%solar_path2%tau_level(lev, ist, j) *                            &
                & EXP( - transmission_scatt_ir_stream%opdpacsun(lev, ist, j))
              transmission_scatt_ir_stream_ad%opdpacsun(lev, ist, j) =                    &
                & transmission_scatt_ir_stream_ad%opdpacsun(lev, ist, j) -                &
                & transmission_aux_ad%solar_path2%tau_level(lev, ist, j) * path2%tau_level(lev, j) *  &
                & EXP( - transmission_scatt_ir_stream%opdpacsun(lev, ist, j))
            ELSE
              tausun_level_ad(lev, j) = tausun_level_ad(lev, j) + transmission_aux_ad%solar_path2%tau_level(lev, ist, j)
            ENDIF
          ENDDO
        ENDDO
      ELSE
        ist = 0
        transmission_ad%tausun_levels_path2(:, j) = 0
        transmission_ad%tausun_total_path2(j)     = 0
        transmission_ad%tausun_levels_path1(:, j) = 0
        transmission_ad%tausun_total_path1(j)     = 0

        tausun_level_ad(:, j) = tausun_level_ad(:, j) + transmission_aux_ad%solar_path2%tau_level(:, ist, j)
        tausun_surf_ad(j)     = tausun_surf_ad(j) + transmission_aux_ad%solar_path2%tau_surf(ist, j)
        tau_level_ad(:, j) = tau_level_ad(:, j) + transmission_aux_ad%solar_path1%tau_level(:, ist, j)
        tau_surf_ad(j)     = tau_surf_ad(j) + transmission_aux_ad%solar_path1%tau_surf(ist, j)
      ENDIF

      DO ist = 0, ircld%nstream(prof)
        transmission_aux_ad%solar_path2%tau_level(:, ist, j) = 0.0_JPRB
      ENDDO
      transmission_aux_ad%solar_path2%od_singlelayer(:, :, j) = 0.0_JPRB
    ENDIF
  ENDDO
!---------------------------------------------------------------------------------------
!AD of compute optical depth and transmittance at surface
!---------------------------------------------------------------------------------------
  DO j = 1, nchannels
    prof    = chanprof(j)%prof
    chan    = chanprof(j)%chan
! as defined in rttov_profaux
    levsurf = aux%s(prof)%nearestlev_surf
! layer above this
    IF (solar(j)) THEN
!---Loop over the streams-----------------------------------------------------------------
      IF (addaerosl .OR. addclouds) THEN
        DO ist = ircld%nstream(prof), 0,  - 1
          IF (path1%tau_surf(j) >= 0) THEN
            tau_surf_ad(j)                             =      &
              & tau_surf_ad(j) + transmission_aux_ad%solar_path1%tau_surf(ist, j) * &
              & transmission_aux%solar_path1%tau_surf_ac(ist, J)
            transmission_aux_ad%solar_path1%tau_surf_ac(ist, J) =     &
              & transmission_aux_ad%solar_path1%tau_surf_ac(ist, J) + &
              & transmission_aux_ad%solar_path1%tau_surf(ist, j) * path1%tau_surf(j)
          ELSE
            tau_surf_ad(j) = tau_surf_ad(j) + transmission_aux_ad%solar_path1%tau_surf(ist, j)
          ENDIF

          IF (path2%tau_surf(j) >= 0) THEN
            tausun_surf_ad(j)                          =      &
              & tausun_surf_ad(j) + transmission_aux_ad%solar_path2%tau_surf(ist, j) * &
              & transmission_aux%solar_path2%tau_surf_ac(ist, J)
            transmission_aux_ad%solar_path2%tau_surf_ac(ist, J) =     &
              & transmission_aux_ad%solar_path2%tau_surf_ac(ist, J) + &
              & transmission_aux_ad%solar_path2%tau_surf(ist, j) * path2%tau_surf(j)
          ELSE
            tausun_surf_ad(j) = tausun_surf_ad(j) + transmission_aux_ad%solar_path2%tau_surf(ist, j)
          ENDIF

          od_frac_ad(j)  = od_frac_ad(j) - &
            & transmission_aux_ad%solar_path1%od_sfrac(ist, J) * &
            & raytracing%pathsat(levsurf - 1, prof) / raytracing%patheff(levsurf - 1, prof)
          transmission_aux_ad%solar_path2%od_frac_ac(ist, J) = transmission_aux_ad%solar_path2%od_frac_ac(ist, J) +  &
            & transmission_aux_ad%solar_path1%od_sfrac(ist, J) * &
            & raytracing%pathsat(levsurf - 1, prof) / raytracing%patheff(levsurf - 1, prof)
          pathsat_patheff_ad(levsurf - 1, prof)     = pathsat_patheff_ad(levsurf - 1, prof) + &
            & ( - path2%od_frac(j) + transmission_aux%solar_path2%od_frac_ac(ist, J)) * &
            & transmission_aux_ad%solar_path1%od_sfrac(ist, J)

          od_frac_ad(j) = od_frac_ad(j) - &
            & transmission_aux_ad%solar_path2%od_sfrac(ist, J) * &
            & raytracing%pathsun(levsurf - 1, prof) / raytracing%patheff(levsurf - 1, prof)
          transmission_aux_ad%solar_path2%od_frac_ac(ist, J) = transmission_aux_ad%solar_path2%od_frac_ac(ist, J) +  &
            & transmission_aux_ad%solar_path2%od_sfrac(ist, J) * &
            & raytracing%pathsun(levsurf - 1, prof) / raytracing%patheff(levsurf - 1, prof)
          pathsun_patheff_ad(levsurf - 1, prof)     = pathsun_patheff_ad(levsurf - 1, prof) + &
            & ( - path2%od_frac(j) + transmission_aux%solar_path2%od_frac_ac(ist, J)) * &
            & transmission_aux_ad%solar_path2%od_sfrac(ist, J)

          od_surf_ac = transmission_scatt_ir_stream%opdpacsun(levsurf, ist, j) + aux%s(prof)%pfraction_surf * (     &
            & transmission_scatt_ir_stream%opdpacsun(levsurf - 1, ist, j) -                                         &
            & transmission_scatt_ir_stream%opdpacsun(levsurf, ist, j))

          pathsat_patheff_ad(levsurf - 1, prof) = pathsat_patheff_ad(levsurf - 1, prof) - &
            & od_surf_ac * transmission_aux%solar_path1%tau_surf_ac(ist, j) * &
            & transmission_aux_ad%solar_path1%tau_surf_ac(ist, j)

          od_surf_ac_ad(j) = od_surf_ac_ad(j) - &
            & transmission_aux%solar_path1%tau_surf_ac(ist, j) * &
            & transmission_aux_ad%solar_path1%tau_surf_ac(ist, j) * &
            & raytracing%pathsat(levsurf - 1, prof) / raytracing%patheff(levsurf - 1, prof)

          od_surf_ac_ad(j)                          =      &
            & od_surf_ac_ad(j) - transmission_aux_ad%solar_path2%tau_surf_ac(ist, j) * &
            & transmission_aux%solar_path2%tau_surf_ac(ist, j)

          IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
            od_surf_ac_ad(j)                                               =      &
              & od_surf_ac_ad(j) + transmission_aux_ad%solar_path2%od_frac_ac(ist, J)
            transmission_scatt_ir_stream_ad%opdpacsun(levsurf - 1, ist, j) =      &
              & transmission_scatt_ir_stream_ad%opdpacsun(levsurf - 1, ist, j) - &
              & transmission_aux_ad%solar_path2%od_frac_ac(ist, J)
          ELSE
            od_surf_ac_ad(j)                                           =      &
              & od_surf_ac_ad(j) + transmission_aux_ad%solar_path2%od_frac_ac(ist, J)
            transmission_scatt_ir_stream_ad%opdpacsun(levsurf, ist, j) =      &
              & transmission_scatt_ir_stream_ad%opdpacsun(levsurf, ist, j) - &
              & transmission_aux_ad%solar_path2%od_frac_ac(ist, J)
          ENDIF
          transmission_scatt_ir_stream_ad%opdpacsun(levsurf, ist, j)     =      &
            & transmission_scatt_ir_stream_ad%opdpacsun(levsurf, ist, j) +      &
            & od_surf_ac_ad(j) * (1._JPRB - aux%s(prof)%pfraction_surf)
          transmission_scatt_ir_stream_ad%opdpacsun(levsurf - 1, ist, j) =      &
            & transmission_scatt_ir_stream_ad%opdpacsun(levsurf - 1, ist, j) +  &
            & od_surf_ac_ad(j) * aux%s(prof)%pfraction_surf
          aux_ad%s(prof)%pfraction_surf                                  = aux_ad%s(prof)%pfraction_surf +      &
            & od_surf_ac_ad(j) * (transmission_scatt_ir_stream%opdpacsun(levsurf - 1, ist, j) -                 &
            & transmission_scatt_ir_stream%opdpacsun(levsurf, ist, j))
          od_surf_ac_ad(j)                                               = 0._JPRB
        ENDDO
      ENDIF
      ! JAH - I made this for hi-res only as in rttov_transmit
      IF (coef%id_sensor == sensor_id_hi) THEN
        IF (coef%tt_val_chn(chan) == 1) THEN
          IF (path2%tau_surf(j) < coef%tt_a0(chan)) THEN
            tausun_surf_ad(j) = 0._jprb
          ENDIF
          IF (path1%tau_ref_surf(j) < coef%tt_a0(chan)) THEN
            tau_surf_ad(j) = 0._jprb
          ENDIF
        ENDIF
      ENDIF
      od_surf = path2%od_level(levsurf, j) + &
        & aux%s(prof)%pfraction_surf * (path2%od_level(levsurf - 1, j) - path2%od_level(levsurf, j))
      pathsat_patheff_ad(levsurf - 1, prof) = pathsat_patheff_ad(levsurf - 1, prof) + &
        & path1%tau_ref_surf(j) * od_surf * tau_surf_ad(j)
      od_surf_ad(j) = od_surf_ad(j) + &
        & path1%tau_ref_surf(j) * tau_surf_ad(j) * &
        & raytracing%pathsat(levsurf - 1, prof) / raytracing%patheff(levsurf - 1, prof)
      od_surf_ad(j) = od_surf_ad(j) + tausun_surf_ad(j) * path2%tau_ref_surf(j)
      IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
        od_surf_ad(j)               = od_surf_ad(j) + od_frac_ad(J)
        od_level_ad(levsurf - 1, j) = od_level_ad(levsurf - 1, j) - od_frac_ad(J)
      ELSE
        od_surf_ad(j)           = od_surf_ad(j) + od_frac_ad(J)
        od_level_ad(levsurf, j) = od_level_ad(levsurf, j) - od_frac_ad(J)
      ENDIF
      od_level_ad(levsurf, j)       = od_level_ad(levsurf, j) + od_surf_ad(j) * (1._JPRB - aux%s(prof)%pfraction_surf)
      od_level_ad(levsurf - 1, j)   = od_level_ad(levsurf - 1, j) + od_surf_ad(j) * aux%s(prof)%pfraction_surf
      aux_ad%s(prof)%pfraction_surf =      &
        & aux_ad%s(prof)%pfraction_surf + od_surf_ad(j) * (path2%od_level(levsurf - 1, j) - path2%od_level(levsurf, j))
      od_surf_ad(j)                 = 0._JPRB
    ENDIF
  ENDDO
!-------------------------------------------
!AD of assemble layer optical depths
!-------------------------------------------
  ! JAH - I made this for hi-res only as in rttov_transmit
  IF (coef%id_sensor == sensor_id_hi) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      IF (solar(j)) THEN
        IF (coef%tt_val_chn(chan) == 1) THEN
          DO lev = 1, nlevels
            IF (path2%tau_ref(lev, j) < coef%tt_a0(chan)) THEN
              tausun_level_ad(lev, j) = 0._jprb
            ENDIF
            IF (path1%tau_ref(lev, j) < coef%tt_a0(chan)) THEN
              tau_level_ad(lev, j) = 0._jprb
            ENDIF
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF
  DO j = 1, nchannels
    IF (solar(j)) THEN
      prof = chanprof(j)%prof
      tau_level_ad(1, j) = 0._jprb
      DO lev = 2, nlevels
        lay = lev - 1

        pathsat_patheff_ad(lay, prof) = pathsat_patheff_ad(lay, prof) + &
          & path1%tau_ref(lev, j) * path2%od_level(lev, j) * tau_level_ad(lev, j)

        od_level_ad(lev, j) = od_level_ad(lev, j) + &
          & path1%tau_ref(lev, j) * tau_level_ad(lev, j) * &
          & raytracing%pathsat(lay, prof) / raytracing%patheff(lay, prof)
      ENDDO
    ENDIF
  ENDDO

  DO prof = 1, nprofiles
    DO lay = 1, nlayers
      zpatheff = 1._jprb / raytracing%patheff(lay, prof)

      raytracing_ad%pathsun(lay, prof) = raytracing_ad%pathsun(lay, prof) + &
          & pathsun_patheff_ad(lay, prof) * zpatheff

      raytracing_ad%patheff(lay, prof) = raytracing_ad%patheff(lay, prof) - &
          & raytracing%pathsun(lay, prof) * pathsun_patheff_ad(lay, prof) * zpatheff ** 2_jpim

      raytracing_ad%pathsat(lay, prof) = raytracing_ad%pathsat(lay, prof) + &
          & pathsat_patheff_ad(lay, prof) * zpatheff

      raytracing_ad%patheff(lay, prof) = raytracing_ad%patheff(lay, prof) - &
          & raytracing%pathsat(lay, prof) * pathsat_patheff_ad(lay, prof) * zpatheff ** 2_jpim
    ENDDO
  ENDDO

  DO j = 1, nchannels
    IF (solar(j)) THEN
      DO lev = 1, nlevels
        od_level_ad(lev, j) = od_level_ad(lev, j) + tausun_level_ad(lev, j) * path2%tau_ref(lev, j)
      ENDDO
    ENDIF
  ENDDO
  DO j = 1, nchannels
    chan = chanprof(j)%chan
    IF (solar(j)) THEN
      od_singlelayer_ad(:, j) = coef%ff_gam(chan) * od_singlelayer_ad(:, j)
      od_level_ad(:, j)       = coef%ff_gam(chan) * od_level_ad(:, j)
    ENDIF
  ENDDO
! ad of level to space optical depths
  opdp_path_ad%sun_level_path2(:,:) = opdp_path_ad%sun_level_path2(:,:) + od_level_ad(:,:)
! ad of single layer optical depths
  DO lay = nlayers, 1,  - 1
    opdp_path_ad%sun_level_path2(lay, :)     = opdp_path_ad%sun_level_path2(lay, :) + od_singlelayer_ad(lay, :)
    opdp_path_ad%sun_level_path2(lay + 1, :) = opdp_path_ad%sun_level_path2(lay + 1, :) - od_singlelayer_ad(lay, :)
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_9_SOLAR_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_transmit_9_solar_ad
