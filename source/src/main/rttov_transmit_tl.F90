SUBROUTINE rttov_transmit_tl( &
            & addaerosl,                       &
            & addclouds,                       &
            & do_lambertian,                   &
            & nlayers,                         &
            & chanprof,                        &
            & chanflag,                        &
            & aux,                             &
            & aux_tl,                          &
            & coef,                            &
            & ircld,                           &
            & geometry,                        &
            & opdp_path,                       &
            & opdp_path_tl,                    &
            & od_level,                        &
            & transmission_levels,             &
            & transmission_total,              &
            & transmission_levels_tl,          &
            & transmission_total_tl,           &
            & transmission_aux,                &
            & transmission_aux_path,           &
            & transmission_aux_path_tl,        &
            & transmission_scatt_ir,           &
            & transmission_scatt_ir_tl,        &
            & transmission_scatt_ir_stream,    &
            & transmission_scatt_ir_stream_tl, &
            & tau_ref,                         &
            & tau_ref_surf,                    &
            & tau_surf,                        &
            & tau_level)
!
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
!  1.1    29/01/2007  Removed polarisation R Saunders
!  1.2    11/12/2007  Return also TL of tau_total and tau_levels since
!                     these are the external outputs of RTTOV (A Geer)
!  1.3    04/06/2008  Fix od_frac and od_frac_ac calculation near surface level (PB PM)
!  1.4    15/07/2009  User defined ToA. Layers distinct from levels (P.Rayer)
!  1.5    03/11/2009  Transmittances / optical depths on levels (A Geer)
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
       & transmission_type_aux,      &
       & rttov_path_transmission,    &
       & transmission_scatt_ir_type, &
       & profile_aux,                &
       & ircld_type, &
       geometry_type
  USE parkind1, ONLY : jpim, jprb, jplm
  USE rttov_math_mod
!INTF_OFF
  USE rttov_const, ONLY : sensor_id_hi,     sec_theta_eff
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE

  LOGICAL(KIND=jplm)              , INTENT(IN)    :: addaerosl
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: addclouds
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: do_lambertian(:)
  TYPE(rttov_chanprof)            , INTENT(IN)    :: chanprof(:)
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: chanflag(SIZE(chanprof))
  INTEGER(KIND=jpim)              , INTENT(IN)    :: nlayers
  TYPE(rttov_coef                ), INTENT(IN)    :: coef
  TYPE(profile_aux               ), INTENT(IN)    :: aux
  TYPE(profile_aux               ), INTENT(IN)    :: aux_tl
  TYPE(ircld_type                ), INTENT(IN)    :: ircld
  TYPE(geometry_type),              INTENT(in)    :: geometry(:)
  REAL(KIND=jprb)                 , INTENT(IN)    :: opdp_path(:,:)
  REAL(KIND=jprb)                 , INTENT(IN)    :: opdp_path_tl(:,:)
  REAL(KIND=jprb)                 , INTENT(IN)    :: transmission_levels(:,:)
  REAL(KIND=jprb)                 , INTENT(IN)    :: transmission_total(:)
  REAL(KIND=jprb)                 , INTENT(INOUT) :: transmission_levels_tl(:,:)
  REAL(KIND=jprb)                 , INTENT(INOUT) :: transmission_total_tl(:)
  TYPE(transmission_type_aux     ), INTENT(IN)    :: transmission_aux
  TYPE(rttov_path_transmission   ), INTENT(IN)    :: transmission_aux_path
  TYPE(rttov_path_transmission   ), INTENT(INOUT) :: transmission_aux_path_tl
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir_tl
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir_stream
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir_stream_tl
  REAL(KIND=jprb)                 , INTENT(IN)    :: tau_ref(nlayers + 1, SIZE(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: tau_ref_surf(SIZE(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: od_level(nlayers + 1, SIZE(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: tau_level(nlayers + 1, SIZE(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: tau_surf(SIZE(chanprof))
!INTF_END

  REAL   (KIND=jprb) :: ref_power 
!   REAL   (KIND=jprb) :: od_surf(SIZE(chanprof))
  REAL   (KIND=jprb) :: od_level_tl(nlayers+1, SIZE(chanprof))
  REAL   (KIND=jprb) :: od_surf_tl(SIZE(chanprof))
  REAL   (KIND=jprb) :: od_surf_ac_tl(0:SIZE(transmission_aux_path%tau_level(1,:,1)), SIZE(chanprof))
  REAL   (KIND=jprb) :: od_frac_tl(SIZE(chanprof)), od_frac_ac_tl
  REAL   (KIND=jprb) :: tau_surf_tl(SIZE(chanprof))
  REAL   (KIND=jprb) :: tau_level_tl(nlayers+1, SIZE(chanprof))

  INTEGER(KIND=jpim) :: lev, lay, chan, j, levsurf
  INTEGER(KIND=jpim) :: prof, nlevels
  INTEGER(KIND=jpim) :: ist
  INTEGER(KIND=jpim) :: nchannels

! allocate automatic arrays using SIZE(transmission_aux_path%tau_level(1,:,1) to keep pgf90 compiler happy
  REAL   (KIND=jprb) :: ztemp(nlayers+1, 0:SIZE(transmission_aux_path%tau_level(1,:,1)))
  REAL   (KIND=jprb) :: ztemp_tl(nlayers+1, 0:SIZE(transmission_aux_path%tau_level(1,:,1)))
  LOGICAL(kind=jplm) :: all_tau_level_ge_zero

  REAL   (KIND=JPRB) :: ZHOOK_HANDLE

!- End of header --------------------------------------------------------
!--------------------------------------------------------------
!1. Assemble layer optical depths
!--------------------------------------------------------------
! - optical depths here are negative except for od_singlelayer
! - in rttov_opdep, already checked that values are sensible
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_TL', 0_jpim, ZHOOK_HANDLE)
  nchannels = SIZE(chanprof)
  nlevels   = nlayers + 1
! single layer optical depths - local variable
  DO j = 1, nchannels
    IF (chanflag(j)) THEN
      DO lay = 1, nlayers
        transmission_aux_path_tl%od_singlelayer(:,lay,j) =  &
          -(opdp_path_tl(lay + 1, j) - opdp_path_tl(lay, j))
      ENDDO
    ENDIF
  ENDDO
! level to space optical depths - local variable
! od_level_tl(:,:) = opdp_path_tl % atm_level(:,:)
! gamma correction of local variables
  DO j = 1, nchannels
    IF (chanflag(j)) THEN
      chan = chanprof(j)%chan
      DO lev = 1, nlevels
        od_level_tl(lev, j)  = coef%ff_gam(chan) * opdp_path_tl(lev, j)
        tau_level_tl(lev, j) = od_level_tl(lev, j) * tau_ref(lev, j)
      ENDDO
      transmission_aux_path_tl%od_singlelayer(:,:,j) = coef%ff_gam(chan) * transmission_aux_path_tl%od_singlelayer(:,:,j)
    ENDIF
  ENDDO
  IF (coef%id_sensor == sensor_id_hi .AND. coef%fmv_model_ver >= 9) THEN
    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        chan = chanprof(j)%chan
        IF (coef%tt_val_chn(chan) == 1) THEN
          DO lev = 1, nlevels
            IF (tau_ref(lev, j) < coef%tt_a0(chan)) THEN
              tau_level_tl(lev, j) = 0._jprb
            ENDIF
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF

!----------------------------------------------------------------------------------------
!2. Compute optical depth and transmittance at surface
!----------------------------------------------------------------------------------------
!   DO j = 1, nchannels
!     IF (chanflag(j)) THEN
!       prof          = chanprof(j)%prof
! ! as defined in rttov_profaux
!       levsurf       = aux%s(prof)%nearestlev_surf
! ! layer above this
! ! arrays here based on layers, not levels
!       IF( (od_level(levsurf - 1, j) - od_level(levsurf, j) > 0._jprb)) THEN
!         od_surf(j) = od_level(levsurf, j) + &
!           aux%s(prof)%pfraction_surf * (od_level(levsurf - 1, j) - od_level(levsurf, j)) !FWD CODE

!         od_surf_tl(j) = od_level_tl(levsurf, j) + &
!           aux%s(prof)%pfraction_surf * (od_level_tl(levsurf - 1, j) - od_level_tl(levsurf, j)) + &
!           aux_tl%s(prof)%pfraction_surf * (od_level(levsurf - 1, j) - od_level(levsurf, j))

!         IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
!           IF(od_surf(j) - od_level(levsurf - 1, j) < 0._jprb) THEN
!             od_frac_tl(j) = od_surf_tl(j) - od_level_tl(levsurf - 1, j)
!           ELSE
!             od_frac_tl(j) = 0._jprb
!           ENDIF
!         ELSE
!           IF(od_surf(j) - od_level(levsurf, j) < 0._jprb) THEN
!             od_frac_tl(j) = od_surf_tl(j) - od_level_tl(levsurf, j)
!           ELSE
!             od_frac_tl(j) = 0._jprb
!           ENDIF
!         ENDIF
!       ELSE
!         od_surf(j) = od_level(levsurf, j) ! FWD CODE
!         od_surf_tl(j) = od_level_tl(levsurf, j)
!         IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
!           IF(od_surf(j) - od_level(levsurf - 1, j) < 0._jprb) THEN
!             od_frac_tl(j) = od_surf_tl(j) - od_level_tl(levsurf - 1, j)
!           ELSE
!             od_frac_tl(j) = 0._jprb
!           ENDIF
!         ELSE
!           od_frac_tl(j) = 0._jprb
!         ENDIF
!       ENDIF
! ! associated transmittance
!       tau_surf_tl(j) = od_surf_tl(j) * tau_ref_surf(j)
!     ENDIF
!   ENDDO

  DO j = 1, nchannels
    IF (chanflag(j)) THEN
      prof          = chanprof(j)%prof
! as defined in rttov_profaux
      levsurf       = aux%s(prof)%nearestlev_surf
! layer above this
! arrays here based on layers, not levels
      od_surf_tl(j) =                                                                                                    &
        & od_level_tl(levsurf, j) + aux_tl%s(prof)%pfraction_surf * (od_level(levsurf - 1, j) - od_level(levsurf, j)) +  &
        & aux%s(prof)%pfraction_surf * (od_level_tl(levsurf - 1, j) - od_level_tl(levsurf, j))
      IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
        od_frac_tl(j) = od_surf_tl(j) - od_level_tl(levsurf - 1, j)
      ELSE
        od_frac_tl(j) = od_surf_tl(j) - od_level_tl(levsurf, j)
      ENDIF
! associated transmittance
      tau_surf_tl(j) = od_surf_tl(j) * tau_ref_surf(j)
    ENDIF
  ENDDO

  IF (coef%id_sensor == sensor_id_hi .AND. coef%fmv_model_ver >= 9) THEN
    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        chan = chanprof(j)%chan
        IF (coef%tt_val_chn(chan) == 1) THEN
          IF (tau_ref_surf(j) < coef%tt_a0(chan)) THEN
            tau_surf_tl(j) = 0._jprb
          ENDIF
        ENDIF
      ENDIF
    ENDDO
  ENDIF

!---Loop over the streams-----------------------------------------------------------------
  IF (addaerosl .OR. addclouds) THEN
    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        prof    = chanprof(j)%prof
        levsurf = aux%s(prof)%nearestlev_surf
! layer above this
        DO ist = 0, ircld%nstream(prof)
          od_surf_ac_tl(ist,j) = transmission_scatt_ir_stream_tl%opdpac(levsurf, ist, j) + aux%s(prof)%pfraction_surf * (     &
            & transmission_scatt_ir_stream_tl%opdpac(levsurf - 1, ist, j) -                                            &
            & transmission_scatt_ir_stream_tl%opdpac(levsurf, ist, j)) + aux_tl%s(prof)%pfraction_surf * (             &
            & transmission_scatt_ir_stream%opdpac(levsurf - 1, ist, j) -                                               &
            & transmission_scatt_ir_stream%opdpac(levsurf, ist, j))
          IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
            od_frac_ac_tl = od_surf_ac_tl(ist,j) - transmission_scatt_ir_stream_tl%opdpac(levsurf - 1, ist, j)
          ELSE
            od_frac_ac_tl = od_surf_ac_tl(ist,j) - transmission_scatt_ir_stream_tl%opdpac(levsurf, ist, j)
          ENDIF
          transmission_aux_path_tl%tau_surf_ac(ist, j) =  - od_surf_ac_tl(ist,j) * &
            & transmission_aux_path%tau_surf_ac(ist, j)
          transmission_aux_path_tl%od_sfrac(ist, j)    =  - od_frac_tl(j) + od_frac_ac_tl
          IF (tau_surf(j) >= 0) THEN
            transmission_aux_path_tl%tau_surf(ist, j) = tau_surf_tl(j) * &
              & transmission_aux_path%tau_surf_ac(ist, j) +      &
              & tau_surf(j) * transmission_aux_path_tl%tau_surf_ac(ist, j)
          ELSE
            transmission_aux_path_tl%tau_surf(ist, j) = tau_surf_tl(j)
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
        IF (addclouds .OR. addaerosl) THEN
          DO ist = 0, ircld%nstream(prof)
            CALL exponential_tl((od_level_tl(:,j) - transmission_scatt_ir_stream_tl%opdpac(:,ist,j)) * ref_power, &
              transmission_aux_path%tau_level_p(:,ist,j), &
              transmission_aux_path_tl%tau_level_p(:,ist,j))

            transmission_aux_path_tl%tau_surf_p(ist,j) = (od_surf_tl(j) - od_surf_ac_tl(ist,j)) * ref_power * &
                                                         transmission_aux_path%tau_surf_p(ist,j)
          ENDDO
        ELSE
          CALL exponential_tl(od_level_tl(:,j) * ref_power, &
            transmission_aux_path%tau_level_p(:,0,j), &
            transmission_aux_path_tl%tau_level_p(:,0,j))

          transmission_aux_path_tl%tau_surf_p(0,j) = od_surf_tl(j) * ref_power * transmission_aux_path%tau_surf_p(0,j)
        ENDIF
        DO ist = 0, ircld%nstream(prof)
          CALL reciprocal_tl(transmission_aux_path%tau_level_p_r(:,ist,j), &
                             transmission_aux_path_tl%tau_level_p(:,ist,j), &
                             transmission_aux_path_tl%tau_level_p_r(:,ist,j))

          transmission_aux_path_tl%tau_surf_p_r(ist,j) = -transmission_aux_path_tl%tau_surf_p(ist,j) * &
                                                          transmission_aux_path%tau_surf_p_r(ist,j) ** 2_jpim
        ENDDO
      ENDIF
    ENDDO
  ENDIF

!-----------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!3. Store transmittances for other streams
!---------------------------------------------------------------------------------------
  IF (addaerosl .OR. addclouds) THEN
    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        all_tau_level_ge_zero = .NOT. ANY((tau_level(:,j) < 0._jprb))        
        prof = chanprof(j)%prof! Profile index

        DO ist = 0, ircld%nstream(prof)
          CALL exponential(-transmission_scatt_ir_stream%opdpac(:, ist, j), &
            ztemp(:,ist))
          CALL exponential_tl(ztemp(:,ist) , &
                             -transmission_scatt_ir_stream_tl%opdpac(:, ist, j), &
                             ztemp_tl(:,ist))
        ENDDO

! DAR: replicate RTTOV 10.2 code
!        WHERE (ztemp(:,0:ircld%nstream(prof)) < EXP(-max_optical_depth))
!          ztemp = EXP(-max_optical_depth)
!          ztemp_tl = 0._jprb
!        ENDWHERE
       
        IF (all_tau_level_ge_zero) THEN
          DO ist = 0, ircld%nstream(prof)
            transmission_aux_path_tl%tau_level(:, ist, j) = &
              tau_level_tl(:, j) * ztemp(:, ist) + &
              tau_level(:, j) * ztemp_tl(:, ist)
          ENDDO
        ELSE
          DO lev = 1, nlevels
            IF (tau_level(lev, j) >= 0._jprb) THEN
              transmission_aux_path_tl%tau_level(lev, 0:ircld%nstream(prof), j) = &
                tau_level_tl(lev, j) * ztemp(lev, 0:ircld%nstream(prof)) + &
                tau_level(lev, j) * ztemp_tl(lev, 0:ircld%nstream(prof))
            ELSE
              transmission_aux_path_tl%tau_level(lev, 0:ircld%nstream(prof), j) = &
                tau_level_tl(lev, j)
            ENDIF
          ENDDO
        ENDIF

        transmission_aux_path_tl%od_singlelayer(:,:,j) = &
          transmission_aux_path_tl%od_singlelayer(:,:,j) + &
          transmission_scatt_ir_tl%opdpacl(:,:,j)

!        WHERE (transmission_aux_path%od_singlelayer(:,:,j) < small_val)
!          transmission_aux_path_tl%od_singlelayer(:,:,j) = 0._jprb
!        ENDWHERE

!        transmission_aux_path%tau_surf(:, j) = &
!          MAX(small_val, transmission_aux_path_tl%tau_surf(:, j))
      ENDIF
    ENDDO

  ELSE ! no addclouds or add aerosol

    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        prof = chanprof(j)%prof
        transmission_total_tl(j)     = tau_surf_tl(j)
        transmission_levels_tl(:, j) = tau_level_tl(:, j)

        ist = 0
        DO lev = 1, nlevels
          transmission_aux_path_tl%tau_level(lev, ist, j) = tau_level_tl(lev, j)
        ENDDO
        transmission_aux_path_tl%tau_surf(ist, j) = tau_surf_tl(j)
        transmission_aux_path_tl%od_sfrac(ist, j) =  - od_frac_tl(j)
      ENDIF
    ENDDO
  ENDIF
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_transmit_tl
