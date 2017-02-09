!
SUBROUTINE rttov_profaux( &
            & opts, &
            & prof, &
            & coef, &
            & aux,  &
            on_COEF_levels)
!
! Description:
! Calculates some variables related to the input profile.
! variables are nearest surface level, nearest cloud top level
! and Debye terms for MW
! The reason of having a separate structure for these
! variables is that the input profiles should be "read only"
! in RTTOV context.
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
!  1.0       01/12/2002  New F90 code with structures (P Brunel A Smith)
!  1.1       27/02/2009  Profile levels to include ToA. Distinguish between
!                        layer arrays and level arrays - size, index
!                        labels, looping (P. Rayer)
!  1.2       02/12/2009  Fixed a number of bugs due to the wrong assumption that cloud
!                        related quantities are defined on levels (thay are layer
!                        average quantities). Marco Matricardi
!  1.3       07/05/2013  Force ice crystal diameters to limits from Table 14 and Table 15 
!                        of RTTOV-9 SVR to avoid negative optical properties (J. Vidot)   
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
! debye coefficients
! Imported Type Definitions:
  USE rttov_types, ONLY :  &
       & rttov_options, &
       & rttov_coef,    &
       & profile_Type,  &
       & profile_aux
!INTF_OFF
  USE rttov_const, ONLY : sensor_id_mw, sensor_id_po, dcoeff, &
       & dgmin_hex, dgmax_hex, dgmin_agg, dgmax_agg 
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
  USE rttov_math_mod, ONLY: INVSQRT, reciprocal
!INTF_ON
  USE parkind1, ONLY : jplm
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options), INTENT(IN)    :: opts
  TYPE(profile_Type ), INTENT(IN)    :: prof(:)! profile
  TYPE(rttov_coef   ), INTENT(IN)    :: coef   ! coefficients
  TYPE(profile_aux  ), INTENT(INOUT) :: aux    ! auxilary profile info
  LOGICAL(jplm), INTENT(IN), OPTIONAL :: on_COEF_levels
!INTF_END
!local
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: iprof
  INTEGER(KIND=jpim) :: lev      , lay
  REAL   (KIND=jprb) :: dp
  REAL   (KIND=jprb) :: v                                           ! temperature ratio
  REAL   (KIND=jprb) :: amcfarq     , bmcfarq     , cmcfarq, zmcfarq
  REAL   (KIND=jprb) :: ztempc
  REAL   (KIND=jprb) :: zradipou_upp, zradipou_low
  REAL   (KIND=jprb) :: bwyser      , nft
  REAL(KIND=jprb), PARAMETER :: rtt      = 273.15_JPRB
  REAL(KIND=jprb), PARAMETER :: rtou_upp =  - 20._JPRB
  REAL(KIND=jprb), PARAMETER :: rtou_low =  - 60._JPRB
  REAL(KIND=JPRB) :: ZHOOK_HANDLE

  REAL(jprb) :: tstar_r(prof(1)%nlayers), ostar_r(prof(1)%nlayers)
  REAL(jprb) :: wstar_r(prof(1)%nlayers)
  REAL(jprb) :: sum1
  INTEGER(jpim) :: iv2lay, iv2lev, iv3lev
  INTEGER(jpim) :: nlayers, nlevels
  LOGICAL(jplm) :: on_COEF_levels1

  on_COEF_levels1 = .FALSE.
  IF(PRESENT(on_COEF_levels)) on_COEF_levels1 = on_COEF_levels
!- End of header --------------------------------------------------------
!-----------------------------------------
! determine cloud top and surface levels
!-----------------------------------------
! in line with coef % dp in rttov_initcoeffs
  IF (LHOOK) CALL DR_HOOK('RTTOV_PROFAUX', 0_jpim, ZHOOK_HANDLE)
  nprofiles = SIZE(prof)
  nlayers = prof(1)%nlayers
  nlevels = prof(1)%nlayers + 1

  DO iprof = 1, nprofiles
! nearest level above surface
    DO lev = prof(iprof)%nlevels - 1, 1,  - 1
      IF (prof(iprof)%s2m%p > prof(iprof)%p(lev)) EXIT
    ENDDO
! case-1: surf lies above lev=nlevels
!         at exit, lev is first level above surface
! case-2: surf lies below lev=nlevels
!         at exit, lev+1=nlevels, there is no level below surface
! case-1: first level below surface
! case-2: first level above surface
    aux%s(iprof)%nearestlev_surf = lev + 1
    dp = prof(iprof)%p(aux%s(iprof)%nearestlev_surf) - prof(iprof)%p(aux%s(iprof)%nearestlev_surf - 1)
    aux%s(iprof)%pfraction_surf  =      &
      & (prof(iprof)%p(aux%s(iprof)%nearestlev_surf) - prof(iprof)%s2m%p) / dp
! NB for case-2, aux % s(iprof) % pfraction_surf -ve
!nearest level above cloud top
    IF (coef%id_sensor /= sensor_id_mw .AND. coef%id_sensor /= sensor_id_po) THEN
      DO lev = prof(iprof)%nlevels - 1, 1,  - 1
        IF (prof(iprof)%ctp > prof(iprof)%p(lev)) EXIT
      ENDDO
      IF( lev > 1 ) THEN
        aux%s(iprof)%nearestlev_ctp = lev + 1
        dp = prof(iprof)%p(aux%s(iprof)%nearestlev_ctp) - prof(iprof)%p(aux%s(iprof)%nearestlev_ctp - 1)
        aux%s(iprof)%pfraction_ctp  =      &
          & (prof(iprof)%p(aux%s(iprof)%nearestlev_ctp) - prof(iprof)%ctp) / dp
        aux%s(iprof)%cfraction      = prof(iprof)%cfraction
      ELSE
        aux%s(iprof)%nearestlev_ctp = prof(iprof)%nlevels -1
        aux%s(iprof)%pfraction_ctp  = 0._JPRB
        aux%s(iprof)%cfraction      = 0._JPRB
      ENDIF
    ELSE
! for micro waves do not consider clouds in the RTTOV basis routines
      aux%s(iprof)%nearestlev_ctp = prof(iprof)%nlevels - 1
      aux%s(iprof)%pfraction_ctp  = 0._JPRB
      aux%s(iprof)%cfraction      = 0._JPRB
    ENDIF
!---------------------------------------------
! debye terms for pewrmittivity calculations
!---------------------------------------------
! Description:
!   To calculate individual debye terms for temperature
!   at each level. There are five debye terms. These
!   will be used in fastem and opdep to calculate
!   permittivity which is required for surface emissivity
!   and cloud modelling
!
! Method:
!   The model is a hybrid of LIEBE MPM 1993 and the PIOM laboratory
!   measurements reported by ELLISON et al. 1999.
    IF (coef%id_sensor == sensor_id_mw .OR. coef%id_sensor == sensor_id_po) THEN
      IF (opts%rt_mw%clw_Data) THEN
! cycle from from level corresponding to the parameter mwcldtp
        DO lev = coef%mwcldtop, prof(iprof)%nlevels
          v = 300.0_JPRB / prof(iprof)%t(lev) - 1.0_JPRB
          aux%debye_prof(1, lev, iprof) = dcoeff(1) - dcoeff(2) * v + dcoeff(3) * v * v
          aux%debye_prof(2, lev, iprof) = dcoeff(4) * aux%debye_prof(1, lev, iprof)
          aux%debye_prof(3, lev, iprof) = dcoeff(5) * v + dcoeff(6)
          aux%debye_prof(4, lev, iprof) = dcoeff(7) * aux%debye_prof(3, lev, iprof)
          aux%debye_prof(5, lev, iprof) = dcoeff(8)
        ENDDO
      ENDIF
    ENDIF
!-----------------------------------------------------------------------------------------
! If ice clouds, convert ice water content into effective generalized diameter
!-----------------------------------------------------------------------------------------
    IF (opts%rt_ir%addclouds .AND. .NOT. opts%rt_ir%user_cld_opt_param .AND. prof(iprof)%ish < 3) THEN
! Calculate upper and lower limits for Ou-Liou effective size
!
      zradipou_upp = 326.3_JPRB + rtou_upp * (12.42_JPRB + rtou_upp * (0.197_JPRB + rtou_upp * 0.0012_JPRB))
      zradipou_low = 326.3_JPRB + rtou_low * (12.42_JPRB + rtou_low * (0.197_JPRB + rtou_low * 0.0012_JPRB))
!
! and convert these to the "generalized" effective size used here (using McFarquhar et al 2003 equation),
! not forgetting the factor of 2 to convert from McFarquhar's radius to a diameter
!
      zradipou_upp =  - 1.56_JPRB + zradipou_upp * (0.388_JPRB + zradipou_upp * 0.00051_JPRB)
      zradipou_upp = 2.0_JPRB * zradipou_upp
      zradipou_low =  - 1.56_JPRB + zradipou_low * (0.388_JPRB + zradipou_low * 0.00051_JPRB)
      zradipou_low = 2.0_JPRB * zradipou_low
      DO lay = 1, prof(iprof)%nlayers
        lev = lay + 1
        IF (prof(iprof)%cloud(6, lay) > 0._jprb) THEN
          ! Use effective diameter from input profile if specified
          IF (prof(iprof)%icede(lay) > 0._jprb) THEN
            aux%dg(lay, iprof) = prof(iprof)%icede(lay)
          ELSE
            IF (prof(iprof)%idg == 1) THEN
  !Scheme by Ou and Liou, 1995, Atmos. Res., 35, 127-138.
              ztempc                  = prof(iprof)%t(lev) - rtt
  ! intermediate factors in calculating the generalized effective diameter
              aux%fac1_dg(lay, iprof) = 326.3_JPRB + ztempc * (12.42_JPRB + ztempc * (0.197_JPRB + ztempc * 0.0012_JPRB))
              aux%fac2_dg(lay, iprof) =      &
                &  - 1.56_JPRB + aux%fac1_dg(lay, iprof) * (0.388_JPRB + aux%fac1_dg(lay, iprof) * 0.00051_JPRB)
              aux%fac3_dg(lay, iprof) = 2.0_JPRB * aux%fac2_dg(lay, iprof)
  !
  ! Take Ou-Liou scheme as being valid only between -20C and -60C
  !
              aux%dg(lay, iprof)      = max(aux%fac3_dg(lay, iprof), zradipou_low)
              aux%dg(lay, iprof)      = min(aux%dg(lay, iprof), zradipou_upp)
            ELSE IF (prof(iprof)%idg == 2) THEN
  !Scheme by Wyser et al. (see McFarquhar et al. (2003))
              bwyser =  - 2.0_JPRB
              IF (prof(iprof)%t(lev) < 273._JPRB) THEN
                bwyser = bwyser + (0.001_JPRB * ((273._JPRB - prof(iprof)%t(lev)) ** 1.5_JPRB) * &
                  Log10(prof(iprof)%cloud(6, lay) / 50._JPRB))
              ENDIF
              aux%fac1_dg(lay, iprof) = 377.4_JPRB + bwyser * (203.3_JPRB + bwyser * (37.91_JPRB + bwyser * 2.3696_JPRB))
              nft = (sqrt(3._JPRB) + 4._JPRB) / (3._JPRB * sqrt(3._JPRB))
              aux%fac2_dg(lay, iprof) = aux%fac1_dg(lay, iprof) / nft
              aux%dg(lay, iprof)      = 2._JPRB * 4._JPRB * aux%fac2_dg(lay, iprof) * sqrt(3._JPRB) / 9._JPRB
            ELSE IF (prof(iprof)%idg == 3) THEN
  !Scheme by Boudala et al., 2002, Int. J. Climatol., 22, 1267-1284.
              ztempc             = prof(iprof)%t(lev) - rtt
              aux%dg(lay, iprof) = 53.005_JPRB * ((prof(iprof)%cloud(6, lay)) ** 0.06_JPRB) * exp(0.013_JPRB * ztempc)
            ELSE IF (prof(iprof)%idg == 4) THEN
  ! Scheme by McFarquhar et al. (2003)
              amcfarq                 = 1.78449_JPRB
              bmcfarq                 = 0.281301_JPRB
              cmcfarq                 = 0.0177166_JPRB
              zmcfarq                 = prof(iprof)%cloud(6, lay)
              aux%fac1_dg(lay, iprof) =      &
                & 10.0_JPRB ** (amcfarq + (bmcfarq * Log10(zmcfarq)) + (cmcfarq * Log10(zmcfarq) * Log10(zmcfarq)))
              aux%dg(lay, iprof)      = 2.0_JPRB * aux%fac1_dg(lay, iprof)
            ELSE
              aux%dg(lay, iprof) = 0._jprb
            ENDIF
          ENDIF

  ! test and force effective size betwwen values from Table 14 and Table 15 of RTTOV-9 SVR   
          IF ( ABS(aux%dg(lay, iprof)) > 0._JPRB ) THEN
            IF (prof(iprof)%ish == 1_JPIM) THEN
              IF ( aux%dg(lay, iprof) < dgmin_hex ) aux%dg(lay, iprof) = dgmin_hex
              IF ( aux%dg(lay, iprof) > dgmax_hex ) aux%dg(lay, iprof) = dgmax_hex
            ENDIF
            IF (prof(iprof)%ish == 2_JPIM) THEN
              IF ( aux%dg(lay, iprof) < dgmin_agg ) aux%dg(lay, iprof) = dgmin_agg
              IF ( aux%dg(lay, iprof) > dgmax_agg ) aux%dg(lay, iprof) = dgmax_agg
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDDO

  IF (on_COEF_levels1) THEN
!FWD ONLY
    aux%on_coef_levels = .TRUE.
! 1 profile layer quantities
!   the layer number agrees with the level number of its upper boundary
! layer N-1 lies between levels N-1 and N
    CALL reciprocal(coef%tstar, tstar_r)
    CALL reciprocal(coef%wstar, wstar_r)
    IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) CALL reciprocal(coef%ostar, ostar_r)

    !DAR add from set_predictors_7
    DO iprof = 1, nprofiles
      aux%t_layer(1:nlayers, iprof) = &
        (prof(iprof)%t(1:nlevels-1) + prof(iprof)%t(2:nlevels)) * 0.5_JPRB
      aux%w_layer(1:nlayers, iprof) = &
        (prof(iprof)%q(1:nlevels-1) + prof(iprof)%q(2:nlevels)) * 0.5_JPRB

      IF (opts%rt_all%use_q2m) THEN
        iv3lev = aux%s(iprof)%nearestlev_surf - 1! nearest level above surface
        iv2lev = aux%s(iprof)%nearestlev_surf    ! nearest level above surface

        IF (iv2lev <= coef%nlevels) THEN
          iv2lay       = iv2lev - 1
          aux%w_layer(iv2lay, iprof) = &
            (prof(iprof)%s2m%q + prof(iprof)%q(iv3lev)) * 0.5_JPRB
        ENDIF
      ENDIF

      IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN
!        aux%o3_layer(layer, iprof) = (prof(iprof)%o3(level - 1) + prof(iprof)%o3(level)) * 0.5_JPRB
        aux%o3_layer(1:nlayers, iprof) = &
          (prof(iprof)%o3(1:nlevels-1) + prof(iprof)%o3(2:nlevels)) * 0.5_JPRB
      ENDIF

! 2. calculate, for layers, deviations from reference profile
      aux%dt(:, iprof) = aux%t_layer(:, iprof)- coef%tstar(:)

      ! if no input O3 profile we still use the input temperature profile for dto
      IF (coef%nozone > 0) &
        aux%dto(:, iprof) = aux%t_layer(:, iprof) - coef%to3star(:)

! 3. calculate (profile / reference profile) ratios; tr wr or
      
      aux%tr(:, iprof) = aux%t_layer(:, iprof) * tstar_r(:)
      aux%wr(:, iprof) = aux%w_layer(:, iprof) * wstar_r(:)
      
! if no input O3 profile, set to reference value (or = 1)

      IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN
        aux%or(:, iprof) = aux%o3_layer(:, iprof) * ostar_r(:)
      ELSE
        aux%or(:, iprof) = 1._JPRB
      ENDIF
    ENDDO

! 4. calculate profile / reference profile sums: tw ww ow

    DO iprof = 1, nprofiles
      aux%tw(1, iprof) = 0._JPRB
      DO lay = 2, nlayers
        ! cumulate overlying layers: weighting tr relates to same layer as dpp
        ! do not need dpp(0) to start
        aux%tw(lay, iprof) = aux%tw(lay - 1, iprof) + &
          coef%dpp(lay - 1) * aux%tr(lay - 1, iprof)
      ENDDO
    ENDDO

!calc profile invariant normalisation coefficient (sum2) - FWD ONLY
    aux%sum(1,1) = coef%dpp(0) * coef%wstar(1)

    DO lay=2, nlayers
      aux%sum(lay,1) = aux%sum(lay-1,1) + coef%dpp(lay - 1) * coef%wstar(lay)
    ENDDO

    IF (coef%nozone > 0 .AND. opts%rt_ir%ozone_Data) THEN
      aux%sum(1,2) = coef%dpp(0) * coef%ostar(1)
      DO lay = 2, nlayers
        aux%sum(lay,2) = aux%sum(lay-1,2) + coef%dpp(lay - 1) * coef%ostar(lay)
      ENDDO
    ELSE
      aux%sum(:,2) = 1._jprb ! FILLER, not used but necessary for the reciprocal below
    ENDIF

    CALL reciprocal(aux%sum, aux%sum)

    DO iprof = 1, nprofiles
      ! cumulating column overlying layer and layer itself
      sum1 = 0._jprb
      DO lay = 1, nlayers
        ! cumulate overlying layers: weighting w or wstar relates to layer below dpp
        ! need dpp(0) to start
        sum1 = sum1 + coef%dpp(lay - 1) * aux%w_layer(lay, iprof)
        aux%ww(lay, iprof) = sum1 * aux%sum(lay,1)
      ENDDO
    ENDDO
    
    ! if no input O3 profile, set to reference value (ow =1)
    IF (coef%nozone > 0) THEN
      IF (opts%rt_ir%ozone_Data) THEN
        DO iprof = 1, nprofiles
          sum1 = 0._JPRB
          DO lay = 1, nlayers
            ! cumulate overlying layers: weighting o or ostar relates to layer below dpp
            ! need dpp(0) to start
            sum1 = sum1 + coef%dpp(lay - 1) * aux%o3_layer(lay, iprof) ! OK, dpp(0) defined 
            aux%ow(lay, iprof) = sum1 * aux%sum(lay,2)
          ENDDO
        ENDDO
      ELSE
        aux%ow(:, :) = 1._JPRB
      ENDIF
    ENDIF
!
    CALL reciprocal(aux%tr, aux%tr_r)

    aux%tw_sqrt = SQRT(aux%tw)
    aux%tw_sqrt(1,:) = 1e-100_jprb ! DARFIX used for division later...
    aux%tw_4rt = SQRT(aux%tw_sqrt)

    CALL INVSQRT(aux%wr, aux%wr_rsqrt)
    aux%wr_sqrt = aux%wr * aux%wr_rsqrt

    CALL reciprocal(aux%ww, aux%ww_r)

    IF (coef%nozone > 0) THEN
      CALL INVSQRT(aux%ow, aux%ow_rsqrt)
      aux%ow_r = aux%ow_rsqrt**2_jpim
      aux%or_sqrt = sqrt(aux%or)
      aux%ow_sqrt = aux%ow * aux%ow_rsqrt
    ENDIF
    
  ELSE
!FWD ONLY
    aux%on_COEF_levels = .FALSE.
  ENDIF
  
  IF (LHOOK) CALL DR_HOOK('RTTOV_PROFAUX', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_profaux
