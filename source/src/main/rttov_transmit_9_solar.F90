
SUBROUTINE rttov_transmit_9_solar( &
            & addaerosl,                    &
            & addclouds,                    &
            & nlayers,                      &
            & chanprof,                     &
            & solar,                        &
            & aux,                          &
            & coef,                         &
            & raytracing,                   &
            & ircld,                        &
            & opdp_path,                    &
            & path2,                        &
            & path1,                        &
            & transmission,                 &
            & transmission_aux,             &
            & transmission_scatt_ir,        &
            & transmission_scatt_ir_stream)
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
!  1.1    29/01/2007  Modified to remove polarisation (R Saunders)
!  1.2    28/08/2007  Optimised (D. Salmond)
!  1.3    27/02/2009  Profile levels to include ToA. Distinguish between
!                     layer arrays and level arrays - size, index
!                     labels, looping (P. Rayer)
!  1.4    03/11/2009  Transmittances / optical depths on levels (A Geer)
!  1.5    02/12/2009  Pathsat, Pathsun and related quantities are now
!                     layer arrays (Marco Matricardi).
!  1.6    05/07/2010  Remove addsolar flag from profiles structure (J Hocking)
!  1.7    04/08/2010  Move addsolar check to calling routine (J Hocking)
!  1.8    14/12/2010  Use traj0_sta%solar array to flag channels for which solar calculations
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
! Imported Type Definitions:
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
  USE rttov_const, ONLY : max_optical_depth, sensor_id_hi, min_tau
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  LOGICAL(KIND=jplm), INTENT(IN)    :: addaerosl
  LOGICAL(KIND=jplm), INTENT(IN)    :: addclouds
  INTEGER(KIND=jpim)              , INTENT(IN)    :: nlayers                  ! Number of pressure levels
  TYPE(rttov_chanprof            ), INTENT(IN)    :: chanprof(:)              ! Channel indices
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: solar(SIZE(chanprof))
  TYPE(rttov_coef                ), INTENT(IN)    :: coef                     ! Coefficients
  TYPE(opdp_path_type            ), INTENT(INOUT) :: opdp_path
  TYPE(rttov_path_traj_trans     ), INTENT(INOUT) :: path2
  TYPE(rttov_path_traj_trans     ), INTENT(INOUT) :: path1
  TYPE(transmission_type         ), INTENT(INOUT) :: transmission             ! Transmittances and single-layer od
  TYPE(transmission_type_aux     ), INTENT(INOUT) :: transmission_aux         ! Transmittances and single-layer od
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir_stream
  TYPE(profile_aux               ), INTENT(IN)    :: aux                      ! auxillary profiles informations
  TYPE(ircld_type                ), INTENT(IN)    :: ircld
  TYPE(raytracing_type           ), INTENT(IN)    :: raytracing
!INTF_END
!local variables:
  REAL   (KIND=jprb) :: od_surf_ac(SIZE(chanprof))
  REAL   (KIND=jprb) :: od_surf(SIZE(chanprof))                 ! sat to surface optical depth (path2)
  REAL   (KIND=jprb) :: od_level(nlayers + 1, SIZE(chanprof))   ! sat to level optical depth (path2)
  REAL   (KIND=jprb) :: od_singlelayer(nlayers,SIZE(chanprof))  ! single layer optical depth (path2)
  INTEGER(KIND=jpim) :: lev, lay, chan, j, prof, ist
  INTEGER(KIND=jpim) :: nlevels
! cloud liquid water local variables
  INTEGER(KIND=jpim) :: levsurf
  INTEGER(KIND=jpim) :: nchannels
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------

! For the most part, solar_path2 variables indicate optical depths and transmittances
! along the combined sun-surface-satellite path based on the effective path length
! given by (pathsat+pathsun). However, some solar_path2 variables are actually along
! the sun-surface path, in particular solar_path2%od_singlelayer and
! solar_path2%od_sfrac.

! All solar_path1 variables are along the surface-satellite path, but note that
! these are derived explicitly from the solar predictor calculations (on path2)
! so that they are consistent with the solar_path2 variables. Therefore these
! will always be slightly different to equivalent the thermal_path1 quantities.

!--------------------------------------------------------------
!1. Assemble layer optical depths and convert to transmittances
!--------------------------------------------------------------
! - optical depths here are negative except for od_singlelayer
! - in rttov_opdep, already checked that values are sensible
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_9_SOLAR', 0_jpim, ZHOOK_HANDLE)
  nchannels = SIZE(chanprof)
  nlevels   = nlayers + 1
  DO j = 1, nchannels
    IF (solar(j)) THEN
      DO lay = 1, nlayers
        od_singlelayer(lay, j) =  - (opdp_path%sun_level_path2(lay + 1, j) - opdp_path%sun_level_path2(lay, j))
      ENDDO
    ENDIF
  ENDDO
  DO j = 1, nchannels
    IF (solar(j)) THEN
      DO lev = 1, nlevels
        od_level(lev, j) = opdp_path%sun_level_path2(lev, j)
      ENDDO
    ENDIF
  ENDDO
  DO j = 1, nchannels
    IF (solar(j)) THEN
      chan = chanprof(j)%chan
      DO lev = 1, nlevels
        od_level(lev, j) = coef%ff_gam(chan) * od_level(lev, j)
      ENDDO
      DO lay = 1, nlayers
        od_singlelayer(lay, j) = coef%ff_gam(chan) * od_singlelayer(lay, j)
      ENDDO
    ENDIF
  ENDDO
! On some computers when optical depth is too thick
! there is an underflow during the conversion in
! transmittances. In that case uncomment following line
! and the declaration statement of max_optical_depth
  DO j = 1, nchannels
    IF (solar(j)) THEN
      DO lev = 1, nlevels
        od_level(lev, j) = MAX(od_level(lev, j),  - max_optical_depth)
      ENDDO
    ENDIF
  ENDDO
  DO j = 1, nchannels
    IF (solar(j)) THEN
      DO lev = 1, nlevels
        path2%tau_level(lev, j) = EXP(od_level(lev, j))
        path2%tau_ref(lev, j)   = path2%tau_level(lev, j)
        transmission%tausun_levels_path2(lev, j) = path2%tau_level(lev, j)
      ENDDO
    ENDIF
  ENDDO

  ! solar_path1 transmittances: the top-level-to-space transmittance is set to 1.0
  DO j = 1, nchannels
    IF (solar(j)) THEN
      prof = chanprof(j)%prof
      path1%tau_level(1, j) = 1.0_jprb
      transmission%tausun_levels_path1(1, j) = path1%tau_level(1, j)
      DO lev = 2, nlevels
        lay = lev - 1
        path1%tau_level(lev, j) = EXP(od_level(lev, j) * &
          & raytracing%pathsat(lay, prof) / raytracing%patheff(lay, prof))
        path1%tau_ref(lev, j)   = path1%tau_level(lev, j)
        transmission%tausun_levels_path1(lev, j) = path1%tau_level(lev, j)
      ENDDO
    ENDIF
  ENDDO

  ! JAH - I made this for hi-res only as in rttov_transmit
  IF (coef%id_sensor == sensor_id_hi) THEN
    DO j = 1, nchannels
      IF (solar(j)) THEN
        chan = chanprof(j)%chan
        IF (coef%tt_val_chn(chan) == 1) THEN
          DO lev = 1, nlevels
            IF (path2%tau_level(lev, j) < coef%tt_a0(chan)) THEN
              path2%tau_level(lev, j) = coef%tt_a1(chan)
              transmission%tausun_levels_path2(lev, j) = path2%tau_level(lev, j)
            ENDIF
            IF (path1%tau_level(lev, j) < coef%tt_a0(chan)) THEN
              path1%tau_level(lev, j) = coef%tt_a1(chan)
              transmission%tausun_levels_path1(lev, j) = path1%tau_level(lev, j)
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
    chan    = chanprof(j)%chan
    prof    = chanprof(j)%prof
! as defined in rttov_profaux
    levsurf = aux%s(prof)%nearestlev_surf
! layer above this
    IF (solar(j)) THEN
      od_surf(j) = od_level(levsurf, j) + &
        & aux%s(prof)%pfraction_surf * (od_level(levsurf - 1, j) - od_level(levsurf, j))
      IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
        path2%od_frac(J) = od_surf(j) - od_level(levsurf - 1, j)
      ELSE
        path2%od_frac(J) = od_surf(j) - od_level(levsurf, j)
      ENDIF
      path2%tau_surf(j)     = EXP(od_surf(j))
      path2%tau_ref_surf(j) = path2%tau_surf(j)
      path1%tau_surf(j)     = EXP(od_surf(j) * &
        & raytracing%pathsat(levsurf - 1, prof) / raytracing%patheff(levsurf - 1, prof))
      path1%tau_ref_surf(j) = path1%tau_surf(j)
      ! JAH - I made this for hi-res only as in rttov_transmit
      IF (coef%id_sensor == sensor_id_hi) THEN
        IF (coef%tt_val_chn(chan) == 1) THEN
          IF (path2%tau_surf(j) < coef%tt_a0(chan)) THEN
            path2%tau_surf(j) = coef%tt_a1(chan)
          ENDIF
          IF (path1%tau_surf(j) < coef%tt_a0(chan)) THEN
            path1%tau_surf(j) = coef%tt_a1(chan)
          ENDIF
        ENDIF
      ENDIF
!---Loop over the streams-----------------------------------------------------------------
      IF (addaerosl .OR. addclouds) THEN
        DO ist = 0, ircld%nstream(prof)
          od_surf_ac(j) = transmission_scatt_ir_stream%opdpacsun(levsurf, ist, j) + aux%s(prof)%pfraction_surf * (     &
            & transmission_scatt_ir_stream%opdpacsun(levsurf - 1, ist, j) -                                            &
            & transmission_scatt_ir_stream%opdpacsun(levsurf, ist, j))
          IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
            transmission_aux%solar_path2%od_frac_ac(ist, J) =      &
              & od_surf_ac(j) - transmission_scatt_ir_stream%opdpacsun(levsurf - 1, ist, j)
          ELSE
            transmission_aux%solar_path2%od_frac_ac(ist, J) =      &
              & od_surf_ac(j) - transmission_scatt_ir_stream%opdpacsun(levsurf, ist, j)
          ENDIF
          transmission_aux%solar_path2%tau_surf_ac(ist, j)     = EXP( - od_surf_ac(j))

          transmission_aux%solar_path1%tau_surf_ac(ist, j)     =     &
            & EXP( - od_surf_ac(j) * &
            & (raytracing%pathsat(levsurf - 1, prof) / raytracing%patheff(levsurf - 1, prof)))

          ! solar_path2%od_sfrac is actually along the sun-surface path (not "path2")
          transmission_aux%solar_path2%od_sfrac(ist, J)        =                  &
            & ( - path2%od_frac(j) + transmission_aux%solar_path2%od_frac_ac(ist, J)) * &
            & raytracing%pathsun(levsurf - 1, prof) / raytracing%patheff(levsurf - 1, prof)

          transmission_aux%solar_path1%od_sfrac(ist, J)        =  &
            & ( - path2%od_frac(j) + transmission_aux%solar_path2%od_frac_ac(ist, J)) * &
            & raytracing%pathsat(levsurf - 1, prof) / raytracing%patheff(levsurf - 1, prof)

          IF (path2%tau_surf(j) >= 0) THEN
            transmission_aux%solar_path2%tau_surf(ist, j) = path2%tau_surf(j) * &
              & transmission_aux%solar_path2%tau_surf_ac(ist, J)
          ELSE
            transmission_aux%solar_path2%tau_surf(ist, j) = path2%tau_surf(j)
          ENDIF

          IF (path1%tau_surf(j) >= 0) THEN
            transmission_aux%solar_path1%tau_surf(ist, j) = path1%tau_surf(j) * &
              & transmission_aux%solar_path1%tau_surf_ac(ist, J)
          ELSE
            transmission_aux%solar_path1%tau_surf(ist, j) = path1%tau_surf(j)
          ENDIF
          transmission_aux%solar_path1%tau_surf_r(ist, j) = 1.0_jprb / transmission_aux%solar_path1%tau_surf(ist, j)

        ENDDO
      ENDIF
!-----------------------------------------------------------------------------------------
    ENDIF
  ENDDO
!---------------------------------------------------------------------------------------
!3. Store transmittances for other streams
!---------------------------------------------------------------------------------------
  IF (addaerosl .OR. addclouds) THEN
    DO j = 1, nchannels
      IF (solar(j)) THEN
        transmission%tausun_total_path2(j) = transmission_aux%solar_path2%tau_surf(0, j)
        transmission%tausun_total_path1(j) = transmission_aux%solar_path1%tau_surf(0, j)
      ENDIF
    ENDDO
  ELSE
    DO j = 1, nchannels
      IF (solar(j)) THEN
        transmission%tausun_total_path2(j) = path2%tau_surf(j)
        transmission%tausun_total_path1(j) = path1%tau_surf(j)
      ENDIF
    ENDDO
  ENDIF
  DO j = 1, nchannels
    prof = chanprof(j)%prof! Profile index
    IF (solar(j)) THEN
      IF (addaerosl .OR. addclouds) THEN
        DO ist = 0, ircld%nstream(prof)
          DO lev = 1, nlevels
            IF (path2%tau_level(lev, j) >= 0) THEN
              transmission_aux%solar_path2%tau_level(lev, ist, j) =      &
                & path2%tau_level(lev, j) * EXP( - transmission_scatt_ir_stream%opdpacsun(lev, ist, j))
            ELSE
              transmission_aux%solar_path2%tau_level(lev, ist, j) = path2%tau_level(lev, j)
            ENDIF
            IF (path1%tau_level(lev, j) >= 0) THEN
              IF (lev > 1) THEN
                lay = lev - 1
                transmission_aux%solar_path1%tau_level(lev, ist, j) =      &
                  & path1%tau_level(lev, j) * EXP( - transmission_scatt_ir_stream%opdpacsun(lev, ist, j) * &
                  & raytracing%pathsat(lay, prof) / raytracing%patheff(lay, prof))
              ELSE
                ! for lev == 1, opdpacsun(lev, ist, j) == 0.
                transmission_aux%solar_path1%tau_level(lev, ist, j) = path1%tau_level(lev, j)
              ENDIF
            ELSE
              transmission_aux%solar_path1%tau_level(lev, ist, j) = path1%tau_level(lev, j)
            ENDIF
            transmission_aux%solar_path1%tau_level_r(lev, ist, j) = &
              & 1.0_jprb / transmission_aux%solar_path1%tau_level(lev, ist, j)
          ENDDO
        ENDDO
        DO lay = 1, nlayers
          ! solar_path2%od_singlelayer is actually along the sun-surface path (not "path2")
          transmission_aux%solar_path2%od_singlelayer(:, lay, j)   =                                   &
            & (od_singlelayer(lay, j) + transmission_scatt_ir%opdpaclsun(:, lay, j)) *  &
            & raytracing%pathsun(lay, prof) / raytracing%patheff(lay, prof)

          transmission_aux%solar_path1%od_singlelayer(:, lay, j)   =                                   &
            & (od_singlelayer(lay, j) + transmission_scatt_ir%opdpaclsun(:, lay, j)) *  &
            & raytracing%pathsat(lay, prof) / raytracing%patheff(lay, prof)
        ENDDO
        path2%od_level(:, j)       = od_level(:, j)
        path2%od_singlelayer(:, j) = od_singlelayer(:, j)
      ELSE
        ist = 0
        transmission_aux%solar_path2%tau_level(:, ist, j) = path2%tau_level(:, j)
        transmission_aux%solar_path2%tau_surf(ist, j)     = path2%tau_surf(j)
        path2%od_level(:, j)                        = od_level(:, j)
        path2%od_singlelayer(:, j)                  = od_singlelayer(:, j)
        transmission_aux%solar_path1%tau_level(:, ist, j) = path1%tau_level(:, j)
        transmission_aux%solar_path1%tau_surf(ist, j)     = path1%tau_surf(j)
      ENDIF

      ! Do non-cloud/cloud layer calculations

      IF (addaerosl .OR. addclouds) THEN
        DO lay = 1, nlayers
          transmission_scatt_ir%opdpext(:, lay, j) =                               &
            & od_singlelayer(lay, j) + transmission_scatt_ir%opdpabs(:, lay, j) +  &
            & transmission_scatt_ir%opdpsca(:, lay, j)
          WHERE (transmission_scatt_ir%opdpext(:, lay, j) /= 0._jprb)
            transmission_scatt_ir%ssa(:, lay, j) =      &
              & transmission_scatt_ir%opdpsca(:, lay, j) / transmission_scatt_ir%opdpext(:, lay, j)
          ENDWHERE
        ENDDO
      ENDIF

    ENDIF ! solar channel
  ENDDO

! See notes in rttov_transmit.F90 regarding fac2
  IF (addaerosl .OR. addclouds) THEN
    DO j = 1, nchannels
      IF (solar(j)) THEN
        prof = chanprof(j)%prof
        DO ist = 0, ircld%nstream(prof)
          DO lay = 1, nlayers
#ifdef RTTOV_XLF
! IBM codepath uses fsel intrinsic that doesn't generated branched code
            transmission_aux%solar_path1%fac2(lay,ist,j) = &
              fsel(transmission_aux%solar_path1%tau_level(lay,ist,j) - min_tau, 1.0_jprb, 0.0_jprb)
#else
! Non-IBM codepath
            lev = lay + 1
            IF (transmission_aux%solar_path1%tau_level(lay,ist,j) < min_tau) THEN
              transmission_aux%solar_path1%fac2(lay,ist,j) = 0.0_jprb
            ELSE
              transmission_aux%solar_path1%fac2(lay,ist,j) = 1.0_jprb
            ENDIF
#endif
          ENDDO
          IF (transmission_aux%solar_path1%tau_level(nlevels,ist,j) < min_tau) THEN
              transmission_aux%solar_path1%fac2(nlevels,ist,j) = 0.0_jprb
          ELSE
              transmission_aux%solar_path1%fac2(nlevels,ist,j) = 1.0_jprb
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_9_SOLAR', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_transmit_9_solar
