SUBROUTINE rttov_transmit_9_solar_k( &
            & addaerosl,                      &
            & addclouds,                      &
            & nlayers,                        &
            & chanprof,                       &
            & solar,                          &
            & aux,                            &
            & aux_k,                          &
            & coef,                           &
            & raytracing,                     &
            & raytracing_k,                   &
            & ircld,                          &
            & opdp_path,                      &
            & opdp_path_k,                    &
            & path2,                          &
            & path1,                          &
            & transmission,                   &
            & transmission_k,                 &
            & transmission_aux,               &
            & transmission_aux_k,             &
            & transmission_scatt_ir,          &
            & transmission_scatt_ir_k,        &
            & transmission_scatt_ir_stream,   &
            & transmission_scatt_ir_stream_k)
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
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0    01/06/2005  Marco Matricardi (ECMWF):
!            --       New routine based on rttov_transmit_k.F90.
!                     Variable trace gases, CO2, N2O,CO and CH4
!                     have been introduced for AIRS and IASI
!  2.0    01/07/2006  Marco Matricardi (ECMWF):
!            --       The contribution of aerosols and clouds
!                     has been added to the total transmission_aux.
!  3.0    06/02/2007  removed polarisation R Saunders
!  4.0    15/09/2009  User defined ToA. Layers distinct from levels (P.Rayer)
!  5.0    03/11/2009  Transmittances / optical depths on levels (A Geer)
!  6.0    02/12/2009  Pathsat, Pathsun and related quantities are now
!                     layer arrays (Marco Matricardi).
!  6.1    05/07/2010  Remove addsolar flag from profiles structure (J Hocking)
!  6.2    04/08/2010  Move addsolar check to calling routine (J Hocking)
!  6.3    14/12/2010  Use traj0_sta%solar array to flag channels for which solar calculations
!                     should be performed (J Hocking)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Adjoint variables
! input transmission_aux_k % solar_path2 % tau_surf and transmission_aux_k % solar_path2 % tau_level set inside integrate_k
!
! input/output aux_k
!
! output predictors_k initialised inside rttov_k (need input
!    intent for memory allocation in calling routine)
!
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
  TYPE(rttov_chanprof            ), INTENT(IN)    :: chanprof(:)
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: solar(SIZE(chanprof))
  TYPE(rttov_coef                ), INTENT(IN)    :: coef
  TYPE(profile_aux               ), INTENT(IN)    :: aux
  TYPE(profile_aux               ), INTENT(INOUT) :: aux_k
  TYPE(ircld_type                ), INTENT(IN)    :: ircld
  TYPE(opdp_path_type            ), INTENT(IN)    :: opdp_path
  TYPE(opdp_path_type            ), INTENT(INOUT) :: opdp_path_k
  TYPE(rttov_path_traj_trans     ), INTENT(IN)    :: path2
  TYPE(rttov_path_traj_trans     ), INTENT(IN)    :: path1
  TYPE(transmission_type         ), INTENT(IN)    :: transmission
  TYPE(transmission_type         ), INTENT(INOUT) :: transmission_k
  TYPE(transmission_type_aux     ), INTENT(IN)    :: transmission_aux
  TYPE(transmission_type_aux     ), INTENT(INOUT) :: transmission_aux_k
  TYPE(raytracing_type           ), INTENT(IN)    :: raytracing
  TYPE(raytracing_type           ), INTENT(INOUT) :: raytracing_k
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_k
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir_stream
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream_k
!INTF_END
!local variables:
  REAL   (KIND=jprb) :: od_singlelayer_k(nlayers, SIZE(chanprof))
  REAL   (KIND=jprb) :: od_surf_ac_k(SIZE(chanprof))
  REAL   (KIND=jprb) :: od_level_k(nlayers + 1, SIZE(chanprof))
  REAL   (KIND=jprb) :: od_surf_k(SIZE(chanprof))
  REAL   (KIND=jprb) :: od_frac_k(SIZE(chanprof))
  REAL   (KIND=jprb) :: tausun_level_k(nlayers + 1, SIZE(chanprof))
  REAL   (KIND=jprb) :: tausun_surf_k(SIZE(chanprof))
  REAL   (KIND=jprb) :: tau_level_k(nlayers + 1, SIZE(chanprof))
  REAL   (KIND=jprb) :: tau_surf_k(SIZE(chanprof))
  REAL   (KIND=jprb) :: pathsat_patheff_k(nlayers, SIZE(chanprof))
  REAL   (KIND=jprb) :: pathsun_patheff_k(nlayers, SIZE(chanprof))
  REAL   (KIND=jprb) :: od_surf, od_surf_ac, zpatheff, tauacsunpath1
  INTEGER(KIND=jpim) :: lev, lay, chan, j, prof, ist, isti, levsurf, nlevels
  INTEGER(KIND=jpim) :: nchannels
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!---------------------------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_9_SOLAR_K', 0_jpim, ZHOOK_HANDLE)
  nchannels             = SIZE(chanprof)
  nlevels               = nlayers + 1
!---------------------------------------------------------------------------------------
!K of store transmittances for other polarisations
!---------------------------------------------------------------------------------------
  od_level_k(:,:)       = 0._JPRB
  tausun_level_k(:,:)   = 0.0_JPRB
  tau_level_k(:,:)      = 0.0_JPRB
  od_singlelayer_k(:,:) = 0.0_JPRB
  tausun_surf_k(:)      = 0.0_JPRB
  tau_surf_k(:)         = 0.0_JPRB
  od_frac_k(:)          = 0.0_JPRB
  od_surf_k(:)          = 0.0_JPRB
  od_surf_ac_k(:)       = 0.0_JPRB
  pathsat_patheff_k(:,:) = 0.0_JPRB
  pathsun_patheff_k(:,:) = 0.0_JPRB
  DO j = nchannels, 1,  - 1
    prof = chanprof(j)%prof
    chan = chanprof(j)%chan
    IF (solar(j)) THEN

      IF (addaerosl .OR. addclouds) THEN
        DO lay = nlayers, 1, -1
          WHERE (transmission_scatt_ir%opdpext(:, lay, j) /= 0._jprb)
            transmission_scatt_ir_k%opdpsca(:, lay, j) =      &
              & transmission_scatt_ir_k%opdpsca(:, lay, j) +  &
              & transmission_scatt_ir_k%ssa(:, lay, j) /      &
              & transmission_scatt_ir%opdpext(:, lay, j)
            transmission_scatt_ir_k%opdpext(:, lay, j) =      &
              & transmission_scatt_ir_k%opdpext(:, lay, j) -  &
              & transmission_scatt_ir_k%ssa(:, lay, j) *      &
              & transmission_scatt_ir%opdpsca(:, lay, j) /     &
              & transmission_scatt_ir%opdpext(:, lay, j) ** 2
          ENDWHERE
          od_singlelayer_k(lay, j)                               =  &
            & od_singlelayer_k(lay, j) + SUM(transmission_scatt_ir_k%opdpext(:, lay, j))
          transmission_scatt_ir_k%opdpabs(:, lay, j)    =  &
            & transmission_scatt_ir_k%opdpabs(:, lay, j) + &
            & transmission_scatt_ir_k%opdpext(:, lay, j)
          transmission_scatt_ir_k%opdpsca(:, lay, j)    =  &
            & transmission_scatt_ir_k%opdpsca(:, lay, j) + &
            & transmission_scatt_ir_k%opdpext(:, lay, j)
        ENDDO
      ENDIF

      IF (addaerosl .OR. addclouds) THEN
        DO isti = 0, 1
          DO lay = nlayers, 1,  - 1
            lev = lay + 1

            od_singlelayer_k(lay, j) = od_singlelayer_k(lay, j) + &
              & transmission_aux_k%solar_path1%od_singlelayer(isti, lay, j) * &
              & raytracing%pathsat(lay, prof) / raytracing%patheff(lay, prof)
            transmission_scatt_ir_k%opdpaclsun(isti, lay, j) = &
              & transmission_scatt_ir_k%opdpaclsun(isti, lay, j) + &
              & transmission_aux_k%solar_path1%od_singlelayer(isti, lay, j) * &
              & raytracing%pathsat(lay, prof) / raytracing%patheff(lay, prof)
            pathsat_patheff_k(lay, j) = pathsat_patheff_k(lay, j) + &
              & transmission_aux_k%solar_path1%od_singlelayer(isti, lay, j) * &
              & (path2%od_singlelayer(lay, j) + transmission_scatt_ir%opdpaclsun(isti, lay, j))

            od_singlelayer_k(lay, j) = od_singlelayer_k(lay, j) + &
              & transmission_aux_k%solar_path2%od_singlelayer(isti, lay, j) * &
              & raytracing%pathsun(lay, prof) / raytracing%patheff(lay, prof)
            transmission_scatt_ir_k%opdpaclsun(isti, lay, j) = &
              & transmission_scatt_ir_k%opdpaclsun(isti, lay, j) + &
              & transmission_aux_k%solar_path2%od_singlelayer(isti, lay, j) * &
              & raytracing%pathsun(lay, prof) / raytracing%patheff(lay, prof)
            pathsun_patheff_k(lay, j) = pathsun_patheff_k(lay, j) + &
              & transmission_aux_k%solar_path2%od_singlelayer(isti, lay, j) * &
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

                tau_level_k(lev, j) = tau_level_k(lev, j) +      &
                  & tauacsunpath1 * transmission_aux_k%solar_path1%tau_level(lev, ist, j)

                transmission_scatt_ir_stream_k%opdpacsun(lev, ist, j) = &
                  & transmission_scatt_ir_stream_k%opdpacsun(lev, ist, j) - &
                  & tauacsunpath1 * path1%tau_level(lev, j) * &
                  & transmission_aux_k%solar_path1%tau_level(lev, ist, j) * &
                  & raytracing%pathsat(lay, prof) / raytracing%patheff(lay, prof)

                pathsat_patheff_k(lay, j) = pathsat_patheff_k(lay, j) - &
                  & tauacsunpath1 * path1%tau_level(lev, j) * &
                  & transmission_aux_k%solar_path1%tau_level(lev, ist, j) * &
                  & transmission_scatt_ir_stream%opdpacsun(lev, ist, j)
              ELSE
                ! for lev == 1, opdpacsun(lev, ist, j) == 0.
                tau_level_k(lev, j) = tau_level_k(lev, j) + &
                  & transmission_aux_k%solar_path1%tau_level(lev, ist, j)
              ENDIF
            ELSE
              tau_level_k(lev, j) = tau_level_k(lev, j) + &
                & transmission_aux_k%solar_path1%tau_level(lev, ist, j)
            ENDIF

            IF (path2%tau_level(lev, j) >= 0) THEN
              tausun_level_k(lev, j)                                = tausun_level_k(lev, j) +      &
                & transmission_aux_k%solar_path2%tau_level(lev, ist, j) *                                    &
                & EXP( - transmission_scatt_ir_stream%opdpacsun(lev, ist, j))
              transmission_scatt_ir_stream_k%opdpacsun(lev, ist, j) =                    &
                & transmission_scatt_ir_stream_k%opdpacsun(lev, ist, j) -                &
                & transmission_aux_k%solar_path2%tau_level(lev, ist, j) * path2%tau_level(lev, j) *  &
                & EXP( - transmission_scatt_ir_stream%opdpacsun(lev, ist, j))
            ELSE
              tausun_level_k(lev, j) = tausun_level_k(lev, j) + transmission_aux_k%solar_path2%tau_level(lev, ist, j)
            ENDIF
          ENDDO
        ENDDO
      ELSE
        ist = 0
        transmission_k%tausun_levels_path2(:, j) = 0
        transmission_k%tausun_total_path2(j)     = 0
        transmission_k%tausun_levels_path1(:, j) = 0
        transmission_k%tausun_total_path1(j)     = 0

        tausun_level_k(:, j) = tausun_level_k(:, j) + transmission_aux_k%solar_path2%tau_level(:, ist, j)
        tausun_surf_k(j)     = tausun_surf_k(j) + transmission_aux_k%solar_path2%tau_surf(ist, j)
        tau_level_k(:, j) = tau_level_k(:, j) + transmission_aux_k%solar_path1%tau_level(:, ist, j)
        tau_surf_k(j)     = tau_surf_k(j) + transmission_aux_k%solar_path1%tau_surf(ist, j)
      ENDIF

      DO ist = 0, ircld%nstream(prof)
        transmission_aux_k%solar_path2%tau_level(:, ist, j) = 0.0_JPRB
      ENDDO
      transmission_aux_k%solar_path2%od_singlelayer(:, :, j) = 0.0_JPRB
    ENDIF
  ENDDO
!---------------------------------------------------------------------------------------
!K of compute optical depth and transmittance at surface
!---------------------------------------------------------------------------------------
  DO j = 1, nchannels
    prof    = chanprof(j)%prof           ! Profile index
    chan    = chanprof(j)%chan
    levsurf = aux%s(prof)%nearestlev_surf
    IF (solar(j)) THEN
!---Loop over the streams-----------------------------------------------------------------
      IF (addaerosl .OR. addclouds) THEN
        DO ist = ircld%nstream(prof), 0,  - 1
          IF (path1%tau_surf(j) >= 0) THEN
            tau_surf_k(j)                             =      &
              & tau_surf_k(j) + transmission_aux_k%solar_path1%tau_surf(ist, j) * &
              & transmission_aux%solar_path1%tau_surf_ac(ist, J)
            transmission_aux_k%solar_path1%tau_surf_ac(ist, J) =     &
              & transmission_aux_k%solar_path1%tau_surf_ac(ist, J) + &
              & transmission_aux_k%solar_path1%tau_surf(ist, j) * path1%tau_surf(j)
          ELSE
            tau_surf_k(j) = tau_surf_k(j) + transmission_aux_k%solar_path1%tau_surf(ist, j)
          ENDIF

          IF (path2%tau_surf(j) >= 0) THEN
            tausun_surf_k(j)                          =      &
              & tausun_surf_k(j) + transmission_aux_k%solar_path2%tau_surf(ist, j) * &
              & transmission_aux%solar_path2%tau_surf_ac(ist, J)
            transmission_aux_k%solar_path2%tau_surf_ac(ist, J) =     &
              & transmission_aux_k%solar_path2%tau_surf_ac(ist, J) + &
              & transmission_aux_k%solar_path2%tau_surf(ist, j) * path2%tau_surf(j)
          ELSE
            tausun_surf_k(j) = tausun_surf_k(j) + transmission_aux_k%solar_path2%tau_surf(ist, j)
          ENDIF

          od_frac_k(j) = od_frac_k(j) - &
            & transmission_aux_k%solar_path1%od_sfrac(ist, J) * &
            & raytracing%pathsat(levsurf - 1, prof) / raytracing%patheff(levsurf - 1, prof)
          transmission_aux_k%solar_path2%od_frac_ac(ist, J) = transmission_aux_k%solar_path2%od_frac_ac(ist, J) +  &
            & transmission_aux_k%solar_path1%od_sfrac(ist, J) * &
            & raytracing%pathsat(levsurf - 1, prof) / raytracing%patheff(levsurf - 1, prof)
          pathsat_patheff_k(levsurf - 1, j) = pathsat_patheff_k(levsurf - 1, j) + &
            & ( - path2%od_frac(j) + transmission_aux%solar_path2%od_frac_ac(ist, J)) * &
            & transmission_aux_k%solar_path1%od_sfrac(ist, J)

          od_frac_k(j) = od_frac_k(j) - &
            & transmission_aux_k%solar_path2%od_sfrac(ist, J) * &
            & raytracing%pathsun(levsurf - 1, prof) / raytracing%patheff(levsurf - 1, prof)
          transmission_aux_k%solar_path2%od_frac_ac(ist, J) = transmission_aux_k%solar_path2%od_frac_ac(ist, J) + &
            & transmission_aux_k%solar_path2%od_sfrac(ist, J) * &
            & raytracing%pathsun(levsurf - 1, prof) / raytracing%patheff(levsurf - 1, prof)
          pathsun_patheff_k(levsurf - 1, j) = pathsun_patheff_k(levsurf - 1, j) + &
            & ( - path2%od_frac(j) + transmission_aux%solar_path2%od_frac_ac(ist, J)) * &
            & transmission_aux_k%solar_path2%od_sfrac(ist, J)

          od_surf_ac = transmission_scatt_ir_stream%opdpacsun(levsurf, ist, j) + aux%s(prof)%pfraction_surf * &
            & (transmission_scatt_ir_stream%opdpacsun(levsurf - 1, ist, j) - &
            & transmission_scatt_ir_stream%opdpacsun(levsurf, ist, j))

          pathsat_patheff_k(levsurf - 1, j) = pathsat_patheff_k(levsurf - 1, j) - &
            & od_surf_ac * transmission_aux%solar_path1%tau_surf_ac(ist, j) * &
            & transmission_aux_k%solar_path1%tau_surf_ac(ist, j)

          od_surf_ac_k(j) = od_surf_ac_k(j) - &
            & transmission_aux%solar_path1%tau_surf_ac(ist, j) * &
            & transmission_aux_k%solar_path1%tau_surf_ac(ist, j) * &
            & raytracing%pathsat(levsurf - 1, prof) / raytracing%patheff(levsurf - 1, prof)

          od_surf_ac_k(j)                          =      &
            & od_surf_ac_k(j) - transmission_aux_k%solar_path2%tau_surf_ac(ist, j) * &
            & transmission_aux%solar_path2%tau_surf_ac(ist, j)

          IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
            od_surf_ac_k(j) = od_surf_ac_k(j) + transmission_aux_k%solar_path2%od_frac_ac(ist, J)
            transmission_scatt_ir_stream_k%opdpacsun(levsurf - 1, ist, j) =     &
              & transmission_scatt_ir_stream_k%opdpacsun(levsurf - 1, ist, j) - &
              & transmission_aux_k%solar_path2%od_frac_ac(ist, J)
          ELSE
            od_surf_ac_k(j)                                           =      &
              & od_surf_ac_k(j) + transmission_aux_k%solar_path2%od_frac_ac(ist, J)
            transmission_scatt_ir_stream_k%opdpacsun(levsurf, ist, j) =     &
              & transmission_scatt_ir_stream_k%opdpacsun(levsurf, ist, j) - &
              & transmission_aux_k%solar_path2%od_frac_ac(ist, J)
          ENDIF
          transmission_scatt_ir_stream_k%opdpacsun(levsurf, ist, j)     =      &
            & transmission_scatt_ir_stream_k%opdpacsun(levsurf, ist, j) +      &
            & od_surf_ac_k(j) * (1._JPRB - aux%s(prof)%pfraction_surf)
          transmission_scatt_ir_stream_k%opdpacsun(levsurf - 1, ist, j) =      &
            & transmission_scatt_ir_stream_k%opdpacsun(levsurf - 1, ist, j) +  &
            & od_surf_ac_k(j) *  aux%s(prof)%pfraction_surf
          aux_k%s(j)%pfraction_surf = aux_k%s(j)%pfraction_surf + od_surf_ac_k(j) &
            &  * (transmission_scatt_ir_stream%opdpacsun(levsurf - 1, ist, j) - &
            & transmission_scatt_ir_stream%opdpacsun(levsurf, ist, j))
          od_surf_ac_k(j) = 0._JPRB
        ENDDO
      ENDIF
      ! JAH - I made this for hi-res only as in rttov_transmit
      IF (coef%id_sensor == sensor_id_hi) THEN
        IF (coef%tt_val_chn(chan) == 1) THEN
          IF (path2%tau_surf(j) < coef%tt_a0(chan)) THEN
            tausun_surf_k(j) = 0._jprb
          ENDIF
          IF (path1%tau_ref_surf(j) < coef%tt_a0(chan)) THEN
            tau_surf_k(j) = 0._jprb
          ENDIF
        ENDIF
      ENDIF
      od_surf = path2%od_level(levsurf, j) + &
        & aux%s(prof)%pfraction_surf * (path2%od_level(levsurf - 1, j) - path2%od_level(levsurf, j))
      pathsat_patheff_k(levsurf - 1, j) = pathsat_patheff_k(levsurf - 1, j) + &
        & path1%tau_ref_surf(j) * od_surf * tau_surf_k(j)
      od_surf_k(j) = od_surf_k(j) + &
        & path1%tau_ref_surf(j) * tau_surf_k(j) * &
        & raytracing%pathsat(levsurf - 1, prof) / raytracing%patheff(levsurf - 1, prof)
      od_surf_k(j) = od_surf_k(j) + tausun_surf_k(j) * path2%tau_ref_surf(j)
      IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
        od_surf_k(j)               = od_surf_k(j) + od_frac_k(J)
        od_level_k(levsurf - 1, j) = od_level_k(levsurf - 1, j) - od_frac_k(J)
      ELSE
        od_surf_k(j)           = od_surf_k(j) + od_frac_k(J)
        od_level_k(levsurf, j) = od_level_k(levsurf, j) - od_frac_k(J)
      ENDIF
      od_level_k(levsurf, j)     = od_level_k(levsurf, j) + od_surf_k(j) * (1._JPRB - aux%s(prof)%pfraction_surf)
      od_level_k(levsurf - 1, j) = od_level_k(levsurf - 1, j) + od_surf_k(j) * aux%s(prof)%pfraction_surf
      aux_k%s(j)%pfraction_surf  =      &
        & aux_k%s(j)%pfraction_surf + od_surf_k(j) * (path2%od_level(levsurf - 1, j) - path2%od_level(levsurf, j))
      od_surf_k(j)               = 0._JPRB
    ENDIF
  ENDDO
!-------------------------------------------
!K of assemble layer optical depths
!-------------------------------------------
  ! JAH - I made this for hi-res only as in rttov_transmit
  IF (coef%id_sensor == sensor_id_hi) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      IF (solar(j)) THEN
        IF (coef%tt_val_chn(chan) == 1) THEN
          DO lev = 1, nlevels
            IF (path2%tau_ref(lev, j) < coef%tt_a0(chan)) THEN
              tausun_level_k(lev, j) = 0._jprb
            ENDIF
            IF (path1%tau_ref(lev, j) < coef%tt_a0(chan)) THEN
              tau_level_k(lev, j) = 0._jprb
            ENDIF
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF
  DO j = 1, nchannels
    IF (solar(j)) THEN
      prof = chanprof(j)%prof
      tau_level_k(1, j) = 0._jprb
      DO lev = 2, nlevels
        lay = lev - 1

        pathsat_patheff_k(lay, j) = pathsat_patheff_k(lay, j) + &
          & path1%tau_ref(lev, j) * path2%od_level(lev, j) * tau_level_k(lev, j)

        od_level_k(lev, j) = od_level_k(lev, j) + &
          & path1%tau_ref(lev, j) * tau_level_k(lev, j) * &
          & raytracing%pathsat(lay, prof) / raytracing%patheff(lay, prof)
      ENDDO
    ENDIF
  ENDDO

  DO j = 1, nchannels
    IF (solar(j)) THEN
      DO lay = 1, nlayers
        prof = chanprof(j)%prof

        zpatheff = 1._jprb / raytracing%patheff(lay, prof)

        raytracing_k%pathsun(lay, j) = raytracing_k%pathsun(lay, j) + &
            & pathsun_patheff_k(lay, j) * zpatheff

        raytracing_k%patheff(lay, j) = raytracing_k%patheff(lay, j) - &
            & raytracing%pathsun(lay, prof) * pathsun_patheff_k(lay, j) * zpatheff ** 2_jpim

        raytracing_k%pathsat(lay, j) = raytracing_k%pathsat(lay, j) + &
            & pathsat_patheff_k(lay, j) * zpatheff

        raytracing_k%patheff(lay, j) = raytracing_k%patheff(lay, j) - &
            & raytracing%pathsat(lay, prof) * pathsat_patheff_k(lay, j) * zpatheff ** 2_jpim
      ENDDO
    ENDIF
  ENDDO

  DO j = 1, nchannels
    IF (solar(j)) THEN
      DO lev = 1, nlevels
        od_level_k(lev, j) = od_level_k(lev, j) + tausun_level_k(lev, j) * path2%tau_ref(lev, j)
      ENDDO
    ENDIF
  ENDDO
  DO j = 1, nchannels
    chan = chanprof(j)%chan
    IF (solar(j)) THEN
      od_singlelayer_k(:, j) = coef%ff_gam(chan) * od_singlelayer_k(:, j)
      od_level_k(:, j)       = coef%ff_gam(chan) * od_level_k(:, j)
    ENDIF
  ENDDO
! k of level to space optical depths
  opdp_path_k%sun_level_path2(:,:) = opdp_path_k%sun_level_path2(:,:) + od_level_k(:,:)
! k of single layer optical depths
  DO lay = nlayers, 1,  - 1
    opdp_path_k%sun_level_path2(lay, :)     = opdp_path_k%sun_level_path2(lay, :) + od_singlelayer_k(lay, :)
    opdp_path_k%sun_level_path2(lay + 1, :) = opdp_path_k%sun_level_path2(lay + 1, :) - od_singlelayer_k(lay, :)
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_9_SOLAR_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_transmit_9_solar_k
