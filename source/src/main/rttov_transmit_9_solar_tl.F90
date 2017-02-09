SUBROUTINE rttov_transmit_9_solar_tl( &
            & addaerosl,                       &
            & addclouds,                       &
            & nlayers,                         &
            & nprofiles,                       &
            & chanprof,                        &
            & solar,                           &
            & aux,                             &
            & aux_tl,                          &
            & coef,                            &
            & raytracing,                      &
            & raytracing_tl,                   &
            & ircld,                           &
            & opdp_path,                       &
            & opdp_path_tl,                    &
            & path2,                           &
            & path1,                           &
            & transmission,                    &
            & transmission_tl,                 &
            & transmission_aux,                &
            & transmission_aux_tl,             &
            & transmission_scatt_ir,           &
            & transmission_scatt_ir_tl,        &
            & transmission_scatt_ir_stream,    &
            & transmission_scatt_ir_stream_tl)
! Description:
! Tangent linear of rttov_transmit
!  To calculate optical depths for a number of channels
!  and profiles from every pressure level to space.
! To interpolate optical depths on to levels of radiative transfer model
! (which, at present, entails only surface transmittance, as
! other optical depths are on *rt* levels) and to convert
! optical depths to transmittances.
!
! Code based on OPDEPTL and RTTAUTL from previous versions of RTTOV
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
!            --       New routine based on rttov_transmit_tl.F90.
!                     Variable trace gases, CO2, N2O,CO and CH4
!                     have been introduced for AIRS and IASI
!  1.1    06/02/2007  Modified to remove polarisation (R Saunders)
!  1.2    15/07/2009  User defined ToA. Layers distinct from levels (P.Rayer)
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
  TYPE(rttov_chanprof)            , INTENT(IN)    :: chanprof(:)
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: solar(SIZE(chanprof))
  INTEGER(KIND=jpim)              , INTENT(IN)    :: nlayers
  INTEGER(KIND=jpim)              , INTENT(IN)    :: nprofiles
  TYPE(rttov_coef                ), INTENT(IN)    :: coef
  TYPE(profile_aux               ), INTENT(IN)    :: aux
  TYPE(profile_aux               ), INTENT(IN)    :: aux_tl
  TYPE(ircld_type                ), INTENT(IN)    :: ircld
  TYPE(opdp_path_type            ), INTENT(IN)    :: opdp_path
  TYPE(opdp_path_type            ), INTENT(INOUT) :: opdp_path_tl
  TYPE(rttov_path_traj_trans     ), INTENT(IN)    :: path2
  TYPE(rttov_path_traj_trans     ), INTENT(IN)    :: path1
  TYPE(transmission_type         ), INTENT(IN)    :: transmission
  TYPE(transmission_type         ), INTENT(INOUT) :: transmission_tl
  TYPE(transmission_type_aux     ), INTENT(IN)    :: transmission_aux
  TYPE(transmission_type_aux     ), INTENT(INOUT) :: transmission_aux_tl
  TYPE(raytracing_type           ), INTENT(IN)    :: raytracing
  TYPE(raytracing_type           ), INTENT(IN)    :: raytracing_tl
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_tl
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir_stream
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir_stream_tl
!INTF_END
!local variables:
  REAL   (KIND=jprb) :: od_singlelayer_tl(nlayers, SIZE(chanprof))
  REAL   (KIND=jprb) :: od_surf_ac_tl(SIZE(chanprof))
  REAL   (KIND=jprb) :: od_level_tl(nlayers + 1, SIZE(chanprof))
  REAL   (KIND=jprb) :: od_surf_tl(SIZE(chanprof))
  REAL   (KIND=jprb) :: od_frac_tl(SIZE(chanprof))
  REAL   (KIND=jprb) :: tausun_level_tl(nlayers + 1, SIZE(chanprof))
  REAL   (KIND=jprb) :: tausun_surf_tl(SIZE(chanprof))
  REAL   (KIND=jprb) :: tau_level_tl(nlayers + 1, SIZE(chanprof))
  REAL   (KIND=jprb) :: tau_surf_tl(SIZE(chanprof))
  REAL   (KIND=jprb) :: pathsat_patheff_tl(nlayers, nprofiles)
  REAL   (KIND=jprb) :: pathsun_patheff_tl(nlayers, nprofiles)
  REAL   (KIND=jprb) :: od_surf, od_surf_ac
  INTEGER(KIND=jpim) :: lev, lay, chan, j, prof, ist, levsurf, nlevels
  INTEGER(KIND=jpim) :: nchannels
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_9_SOLAR_TL', 0_jpim, ZHOOK_HANDLE)
  nchannels = SIZE(chanprof)
  nlevels   = nlayers + 1
!----------------------------------------
!2. Compute layer to space optical depths
!----------------------------------------
! notes: apply gamma correction; check value is sensible and constrain
! if necessary.
! note that optical depth in the calculations is negative
  DO j = 1, nchannels
    IF (solar(j)) THEN
      DO lay = 1, nlayers
        od_singlelayer_tl(lay, j) =  - (opdp_path_tl%sun_level_path2(lay + 1, j) - opdp_path_tl%sun_level_path2(lay, j))
      ENDDO
    ENDIF
  ENDDO
  DO j = 1, nchannels
    IF (solar(j)) THEN
      DO lev = 1, nlevels
        od_level_tl(lev, j) = opdp_path_tl%sun_level_path2(lev, j)
      ENDDO
    ENDIF
  ENDDO
  DO j = 1, nchannels
    chan = chanprof(j)%chan
    IF (solar(j)) THEN
      DO lev = 1, nlevels
        od_level_tl(lev, j) = coef%ff_gam(chan) * od_level_tl(lev, j)
      ENDDO
      DO lay = 1, nlayers
        od_singlelayer_tl(lay, j) = coef%ff_gam(chan) * od_singlelayer_tl(lay, j)
      ENDDO
    ENDIF
  ENDDO
! associated transmittances
  DO j = 1, nchannels
    IF (solar(j)) THEN
      DO lev = 1, nlevels
        tausun_level_tl(lev, j) = od_level_tl(lev, j) * path2%tau_ref(lev, j)
      ENDDO
    ENDIF
  ENDDO
  DO j = 1, nchannels
    IF (solar(j)) THEN
      prof = chanprof(j)%prof
      tau_level_tl(1, j) = 0._jprb
      DO lev = 2, nlevels
        lay = lev - 1

        pathsat_patheff_tl(lay, prof) = &
          & raytracing_tl%pathsat(lay, prof) / raytracing%patheff(lay, prof) - &
          & raytracing%pathsat(lay, prof) * raytracing_tl%patheff(lay, prof) / &
          & raytracing%patheff(lay, prof) ** 2_jpim

        pathsun_patheff_tl(lay, prof) = &
          & raytracing_tl%pathsun(lay, prof) / raytracing%patheff(lay, prof) - &
          & raytracing%pathsun(lay, prof) * raytracing_tl%patheff(lay, prof) / &
          & raytracing%patheff(lay, prof) ** 2_jpim

        tau_level_tl(lev, j) = path1%tau_ref(lev, j) * &
          & (path2%od_level(lev, j) * pathsat_patheff_tl(lay, prof) + od_level_tl(lev, j) * &
          & raytracing%pathsat(lay, prof) / raytracing%patheff(lay, prof))
      ENDDO
    ENDIF
  ENDDO
  ! JAH - I made this for hi-res only as in rttov_transmit
  IF (coef%id_sensor == sensor_id_hi) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      IF (solar(j)) THEN
        IF (coef%tt_val_chn(chan) == 1) THEN
          DO lev = 1, nlevels
            IF (path2%tau_ref(lev, j) < coef%tt_a0(chan)) THEN
              tausun_level_tl(lev, j) = 0._jprb
            ENDIF
            IF (path1%tau_ref(lev, j) < coef%tt_a0(chan)) THEN
              tau_level_tl(lev, j) = 0._jprb
            ENDIF
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF
!---------------------------------------------------------------------------------------
!2. Compute optical depth and transmittance at surface
!---------------------------------------------------------------------------------------
  DO j = 1, nchannels
    prof    = chanprof(j)%prof
    chan    = chanprof(j)%chan
! as defined in rttov_profaux
    levsurf = aux%s(prof)%nearestlev_surf
! layer above this
    IF (solar(j)) THEN
      od_surf_tl(j) = od_level_tl(levsurf, j) +                                                      &
        & aux_tl%s(prof)%pfraction_surf * (path2%od_level(levsurf - 1, j) - path2%od_level(levsurf, j)) +  &
        & aux%s(prof)%pfraction_surf * (od_level_tl(levsurf - 1, j) - od_level_tl(levsurf, j))
      IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
        od_frac_tl(J) = od_surf_tl(j) - od_level_tl(levsurf - 1, j)
      ELSE
        od_frac_tl(J) = od_surf_tl(j) - od_level_tl(levsurf, j)
      ENDIF
      tausun_surf_tl(j) = od_surf_tl(j) * path2%tau_ref_surf(j)
      od_surf = path2%od_level(levsurf, j) + &
        & aux%s(prof)%pfraction_surf * (path2%od_level(levsurf - 1, j) - path2%od_level(levsurf, j))
      tau_surf_tl(j) = path1%tau_ref_surf(j) * &
        & (od_surf * pathsat_patheff_tl(levsurf - 1, prof) + od_surf_tl(j) * &
        & raytracing%pathsat(levsurf - 1, prof) / raytracing%patheff(levsurf - 1, prof))
      ! JAH - I made this for hi-res only as in rttov_transmit
      IF (coef%id_sensor == sensor_id_hi) THEN
        IF (coef%tt_val_chn(chan) == 1) THEN
          IF (path2%tau_ref_surf(j) < coef%tt_a0(chan)) THEN
            tausun_surf_tl(j) = 0._jprb
          ENDIF
          IF (path1%tau_ref_surf(j) < coef%tt_a0(chan)) THEN
            tau_surf_tl(j) = 0._jprb
          ENDIF
        ENDIF
      ENDIF
!---Loop over the streams-----------------------------------------------------------------
      IF (addaerosl .OR. addclouds) THEN
        DO ist = 0, ircld%nstream(prof)
          od_surf_ac_tl(j) = transmission_scatt_ir_stream_tl%opdpacsun(levsurf, ist, j) + aux%s(prof)%pfraction_surf &
            &  * (transmission_scatt_ir_stream_tl%opdpacsun(levsurf - 1, ist, j) -                                   &
            & transmission_scatt_ir_stream_tl%opdpacsun(levsurf, ist, j)) + aux_tl%s(prof)%pfraction_surf * (        &
            & transmission_scatt_ir_stream%opdpacsun(levsurf - 1, ist, j) -                                          &
            & transmission_scatt_ir_stream%opdpacsun(levsurf, ist, j))
          IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
            transmission_aux_tl%solar_path2%od_frac_ac(ist, J) =      &
              & od_surf_ac_tl(j) - transmission_scatt_ir_stream_tl%opdpacsun(levsurf - 1, ist, j)
          ELSE
            transmission_aux_tl%solar_path2%od_frac_ac(ist, J) =      &
              & od_surf_ac_tl(j) - transmission_scatt_ir_stream_tl%opdpacsun(levsurf, ist, j)
          ENDIF

          transmission_aux_tl%solar_path2%tau_surf_ac(ist, j) =      &
            &  - od_surf_ac_tl(j) * transmission_aux%solar_path2%tau_surf_ac(ist, j)

          od_surf_ac = transmission_scatt_ir_stream%opdpacsun(levsurf, ist, j) + aux%s(prof)%pfraction_surf * (     &
            & transmission_scatt_ir_stream%opdpacsun(levsurf - 1, ist, j) -                                         &
            & transmission_scatt_ir_stream%opdpacsun(levsurf, ist, j))

          transmission_aux_tl%solar_path1%tau_surf_ac(ist, j) =      &
            & transmission_aux%solar_path1%tau_surf_ac(ist, j) * &
            & (- od_surf_ac * pathsat_patheff_tl(levsurf - 1, prof) - od_surf_ac_tl(j) * &
            & raytracing%pathsat(levsurf - 1, prof) / raytracing%patheff(levsurf - 1, prof))

          transmission_aux_tl%solar_path2%od_sfrac(ist, J)    =                                   &
            & ( - od_frac_tl(j) + transmission_aux_tl%solar_path2%od_frac_ac(ist, J)) * &
            &     raytracing%pathsun(levsurf - 1, prof) / raytracing%patheff(levsurf - 1, prof) +   &
            & ( - path2%od_frac(j) + transmission_aux%solar_path2%od_frac_ac(ist, J)) * &
            & pathsun_patheff_tl(levsurf - 1, prof)

          transmission_aux_tl%solar_path1%od_sfrac(ist, J)    =                                   &
            & ( - od_frac_tl(j) + transmission_aux_tl%solar_path2%od_frac_ac(ist, J)) * &
            &     raytracing%pathsat(levsurf - 1, prof) / raytracing%patheff(levsurf - 1, prof) +   &
            & ( - path2%od_frac(j) + transmission_aux%solar_path2%od_frac_ac(ist, J)) * &
            & pathsat_patheff_tl(levsurf - 1, prof)

          IF (path2%tau_surf(j) >= 0) THEN
            transmission_aux_tl%solar_path2%tau_surf(ist, j) = &
              & tausun_surf_tl(j) * transmission_aux%solar_path2%tau_surf_ac(ist, J) + &
              & transmission_aux_tl%solar_path2%tau_surf_ac(ist, J) * path2%tau_surf(j)
          ELSE
            transmission_aux_tl%solar_path2%tau_surf(ist, j) = tausun_surf_tl(j)
          ENDIF

          IF (path1%tau_surf(j) >= 0) THEN
            transmission_aux_tl%solar_path1%tau_surf(ist, j) = &
              & tau_surf_tl(j) * transmission_aux%solar_path1%tau_surf_ac(ist, J) + &
              & transmission_aux_tl%solar_path1%tau_surf_ac(ist, J) * path1%tau_surf(j)
          ELSE
            transmission_aux_tl%solar_path1%tau_surf(ist, j) = tau_surf_tl(j)
          ENDIF
        ENDDO
      ENDIF
!-----------------------------------------------------------------------------------------
    ENDIF
  ENDDO
!---------------------------------------------------------------------------------------
!3. Store transmittances for other streams
!---------------------------------------------------------------------------------------
  DO j = 1, nchannels
    prof = chanprof(j)%prof
    chan = chanprof(j)%chan
    IF (solar(j)) THEN
      IF (addaerosl .OR. addclouds) THEN
        DO ist = 0, ircld%nstream(prof)
          DO lev = 1, nlevels
            IF (path2%tau_level(lev, j) >= 0) THEN
              transmission_aux_tl%solar_path2%tau_level(lev, ist, j) = &
                & (tausun_level_tl(lev, j) - &
                & transmission_scatt_ir_stream_tl%opdpacsun(lev, ist, j) * path2%tau_level(lev, j)) * &
                & EXP( - transmission_scatt_ir_stream%opdpacsun(lev, ist, j))
            ELSE
              transmission_aux_tl%solar_path2%tau_level(lev, ist, j) = tausun_level_tl(lev, j)
            ENDIF

            IF (path1%tau_level(lev, j) >= 0) THEN
              IF (lev > 1) THEN
                lay = lev - 1
                transmission_aux_tl%solar_path1%tau_level(lev, ist, j) =      &
                  & EXP( - transmission_scatt_ir_stream%opdpacsun(lev, ist, j) * &
                  & raytracing%pathsat(lay, prof) / raytracing%patheff(lay, prof)) * &
                  & (tau_level_tl(lev, j) + path1%tau_level(lev, j) * &
                  & ( - transmission_scatt_ir_stream_tl%opdpacsun(lev, ist, j) * &
                  & raytracing%pathsat(lay, prof) / raytracing%patheff(lay, prof) + &
                  & ( - transmission_scatt_ir_stream%opdpacsun(lev, ist, j) * pathsat_patheff_tl(lay, prof))))
              ELSE
                ! for lev == 1, opdpacsun(lev, ist, j) == 0.
                transmission_aux_tl%solar_path1%tau_level(lev, ist, j) = tau_level_tl(lev, j)
              ENDIF
            ELSE
              transmission_aux_tl%solar_path1%tau_level(lev, ist, j) = tau_level_tl(lev, j)
            ENDIF
          ENDDO
        ENDDO
        DO lay = 1, nlayers
          lev = lay + 1
          transmission_aux_tl%solar_path2%od_singlelayer(:, lay, j)   =                                        &
            & (od_singlelayer_tl(lay, j) + transmission_scatt_ir_tl%opdpaclsun(:, lay, j)) *           &
            & raytracing%pathsun(lay, prof) / raytracing%patheff(lay, prof) +  &
            & (path2%od_singlelayer(lay, j) + transmission_scatt_ir%opdpaclsun(:, lay, j)) *           &
            & pathsun_patheff_tl(lay, prof)

          transmission_aux_tl%solar_path1%od_singlelayer(:, lay, j)   =                                        &
            & (od_singlelayer_tl(lay, j) + transmission_scatt_ir_tl%opdpaclsun(:, lay, j)) *           &
            & raytracing%pathsat(lay, prof) / raytracing%patheff(lay, prof) +  &
            & (path2%od_singlelayer(lay, j) + transmission_scatt_ir%opdpaclsun(:, lay, j)) *           &
            & pathsat_patheff_tl(lay, prof)
        ENDDO
      ELSE
        ist = 0
        transmission_tl%tausun_total_path2(j)     = tausun_surf_tl(j)
        transmission_tl%tausun_levels_path2(:, j) = tausun_level_tl(:, j)
        transmission_tl%tausun_total_path1(j)     = tau_surf_tl(j)
        transmission_tl%tausun_levels_path1(:, j) = tau_level_tl(:, j)
        transmission_aux_tl%solar_path2%tau_level(:, ist, j) = tausun_level_tl(:, j)
        transmission_aux_tl%solar_path2%tau_surf(ist, j)     = tausun_surf_tl(j)
        transmission_aux_tl%solar_path1%tau_level(:, ist, j) = tau_level_tl(:, j)
        transmission_aux_tl%solar_path1%tau_surf(ist, j)     = tau_surf_tl(j)
      ENDIF

      IF (addaerosl .OR. addclouds) THEN
        DO lay = 1, nlayers
          transmission_scatt_ir_tl%opdpext(:, lay, j) =                                  &
            & od_singlelayer_tl(lay, j) + transmission_scatt_ir_tl%opdpabs(:, lay, j) +  &
            & transmission_scatt_ir_tl%opdpsca(:, lay, j)
          WHERE (transmission_scatt_ir%opdpext(:, lay, j) /= 0._jprb)
            transmission_scatt_ir_tl%ssa(:, lay, j) =            &
              & transmission_scatt_ir_tl%opdpsca(:, lay, j) /    &
              & transmission_scatt_ir%opdpext(:, lay, j)         &
              & - transmission_scatt_ir_tl%opdpext(:, lay, j) *  &
              & transmission_scatt_ir%opdpsca(:, lay, j)         &
              &  / transmission_scatt_ir%opdpext(:, lay, j) ** 2
          ENDWHERE
        ENDDO
      ENDIF

    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_9_SOLAR_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_transmit_9_solar_tl
