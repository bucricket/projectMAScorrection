SUBROUTINE rttov_profaux_ad( &
            & opts,    &
            & prof,    &
            & prof_ad, &
            & coef,    &
            & aux,     &
            & aux_ad)
!
! Description:
! AD of rttov_profaux
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
!  1.0    07/10/2004 Added history
!  1.1    29/03/2005 Add end of header comment (J. Cameron)
!  1.2    04/02/2008 opts%lgradp option for TL/AD of pressure levels (N. Bormann)
!  1.3    15/08/2009 User defined ToA. Layers distinct from levels (P.Rayer)
!  1.4    02/12/2009 Fixed a number of bugs due to the wrong assumption that cloud
!                    related quantities are defined on levels (thay are layer
!                    average quantities). Marco Matricardi
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
       & rttov_options, &
       & rttov_coef,    &
       & profile_Type,  &
       & profile_aux
!INTF_OFF
  USE rttov_const, ONLY : sensor_id_mw, sensor_id_po, dcoeff, &
                          & inst_id_ssmis, inst_id_ssmisz, &
                          & dgmin_hex, dgmax_hex, dgmin_agg, dgmax_agg 
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb, jplm
  USE rttov_math_mod, ONLY: reciprocal, &
                            INVSQRT_AD, reciprocal_ad, sqrt_ad
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options), INTENT(IN)    :: opts
  TYPE(profile_Type ), INTENT(IN)    :: prof   (:)
  TYPE(profile_Type ), INTENT(INOUT) :: prof_ad(:)
  TYPE(rttov_coef   ), INTENT(IN)    :: coef
  TYPE(profile_aux  ), INTENT(IN)    :: aux
  TYPE(profile_aux  ), INTENT(INOUT) :: aux_ad

!INTF_END
! local
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: iprof
  INTEGER(KIND=jpim) :: lev      , lay
  REAL   (KIND=jprb) :: dp
  REAL   (KIND=jprb) :: dp_ad
  REAL   (KIND=jprb) :: v                                                       ! temperature ratio
  REAL   (KIND=jprb) :: v_ad                                                    ! temperature ratio
  REAL   (KIND=jprb) :: amcfarq     , bmcfarq     , cmcfarq, zmcfarq, zmcfarq_ad
  REAL   (KIND=jprb) :: ztempc      , ztempc_ad
  REAL   (KIND=jprb) :: zradipou_upp, zradipou_low
  REAL   (KIND=jprb) :: bwyser      , bwyser_ad   , nft
  REAL(KIND=jprb), PARAMETER :: rtt      = 273.15_JPRB
  REAL(KIND=jprb), PARAMETER :: rtou_upp =  - 20._JPRB
  REAL(KIND=jprb), PARAMETER :: rtou_low =  - 60._JPRB
  REAL(KIND=JPRB) :: ZHOOK_HANDLE

  REAL(jprb) :: tstar_r(prof(1)%nlayers), ostar_r(prof(1)%nlayers)
  REAL(jprb) :: wstar_r(prof(1)%nlayers)
  REAL(jprb) :: sum1
  INTEGER(jpim) :: iv2lay, iv2lev, iv3lev
  INTEGER(jpim) :: nlayers, nlevels

!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_PROFAUX_AD', 0_jpim, ZHOOK_HANDLE)
  nprofiles = size(prof)
  nlayers = prof(1)%nlayers
  nlevels = prof(1)%nlayers + 1

  DO iprof = 1, nprofiles
    zmcfarq_ad = 0._jprb
    ztempc_ad  = 0._jprb
    bwyser_ad  = 0._jprb
    dp_ad      = 0._jprb
!-----------------------------------------------------------------------------------------
! AD for ice water content into effective generalized diameter
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
      DO lay = prof(iprof)%nlayers, 1,  - 1
        lev = lay + 1
        IF (prof(iprof)%cloud(6, lay) > 0._jprb) THEN
          IF ( ABS(aux%dg(lay, iprof)) > 0._JPRB ) THEN
            IF (prof(iprof)%ish == 1_JPIM) THEN
              IF ( aux%dg(lay, iprof) < dgmin_hex ) THEN
                aux_ad%dg(lay, iprof) = 0._JPRB
              ENDIF
              IF ( aux%dg(lay, iprof) > dgmax_hex ) THEN
                aux_ad%dg(lay, iprof) = 0._JPRB
              ENDIF
            ENDIF
            IF (prof(iprof)%ish == 2_JPIM) THEN
              IF ( aux%dg(lay, iprof) < dgmin_agg ) THEN
                aux_ad%dg(lay, iprof) = 0._JPRB
              ENDIF
              IF ( aux%dg(lay, iprof) > dgmax_agg ) THEN
                aux_ad%dg(lay, iprof) = 0._JPRB
              ENDIF
            ENDIF
          ENDIF
          ! Use effective diameter from input profile if specified
          IF (prof(iprof)%icede(lay) > 0._jprb) THEN
            prof_ad(iprof)%icede(lay) = prof_ad(iprof)%icede(lay) + aux_ad%dg(lay, iprof)
          ELSE
            IF (prof(iprof)%idg == 1) THEN
  !Scheme by Ou and Liou, 1995, Atmos. Res., 35, 127-138.
              ztempc = prof(iprof)%t(lev) - rtt
  !
  ! Take Ou-Liou scheme as being valid only between 20C and 60C
  !
              IF (aux%fac3_dg(lay, iprof) < zradipou_low .OR. aux%fac3_dg(lay, iprof) > zradipou_upp) THEN
                aux_ad%dg(lay, iprof) = 0._JPRB
              ELSE
                aux_ad%fac3_dg(lay, iprof) = aux_ad%fac3_dg(lay, iprof) + aux_ad%dg(lay, iprof)
              ENDIF
              aux_ad%fac2_dg(lay, iprof) = aux_ad%fac2_dg(lay, iprof) + 2.0_JPRB * aux_ad%fac3_dg(lay, iprof)
              aux_ad%fac1_dg(lay, iprof) = aux_ad%fac1_dg(lay, iprof) + 0.388_JPRB * aux_ad%fac2_dg(lay, iprof)
              aux_ad%fac1_dg(lay, iprof) =      &
                & aux_ad%fac1_dg(lay, iprof) + aux_ad%fac2_dg(lay, iprof) * 0.00051_JPRB * 2 * aux%fac1_dg(lay, iprof)
              ztempc_ad                  = ztempc_ad + aux_ad%fac1_dg(lay, iprof) *      &
                & (12.42_JPRB + 2._JPRB * 0.197_JPRB * ztempc + 3._JPRB * 0.0012_JPRB * ztempc * ztempc)
              prof_ad(iprof)%t(lev)      = prof_ad(iprof)%t(lev) + ztempc_ad
              ztempc_ad                  = 0._jprb
            ELSE IF (prof(iprof)%idg == 2) THEN
  !Scheme by Wyser et al. (see McFarquhar et al. (2003))
              bwyser =  - 2.0_JPRB
              IF (prof(iprof)%t(lev) < 273._JPRB) THEN
                bwyser = bwyser + (0.001_JPRB * &
                  ((273._JPRB - prof(iprof)%t(lev)) ** 1.5_JPRB) * Log10(prof(iprof)%cloud(6, lay) / 50._JPRB))
              ENDIF
              nft = (sqrt(3._JPRB) + 4._JPRB) / (3._JPRB * sqrt(3._JPRB))
              aux_ad%fac2_dg(lay, iprof) =      &
                & aux_ad%fac2_dg(lay, iprof) + aux_ad%dg(lay, iprof) * 2._JPRB * 4._JPRB * sqrt(3._JPRB) / 9._JPRB
              aux_ad%fac1_dg(lay, iprof) = aux_ad%fac1_dg(lay, iprof) + aux_ad%fac2_dg(lay, iprof) / nft
              bwyser_ad                  = bwyser_ad +      &
                & aux_ad%fac1_dg(lay, iprof) * (203.3_JPRB + 2 * bwyser * 37.91_JPRB + 3 * bwyser * bwyser * 2.3696_JPRB)
              IF (prof(iprof)%t(lev) < 273._JPRB) THEN
                prof_ad(iprof)%t(lev)        = prof_ad(iprof)%t(lev) -                                    &
                  & bwyser_ad * 0.001_JPRB * 1.5_JPRB * ((273._JPRB - prof(iprof)%t(lev)) ** 0.5_JPRB) *  &
                  & Log10(prof(iprof)%cloud(6, lay) / 50._JPRB)
                prof_ad(iprof)%cloud(6, lay) = prof_ad(iprof)%cloud(6, lay) + &
                  bwyser_ad * 0.001_JPRB * ((273._JPRB - prof(iprof)%t(lev)) ** 1.5_JPRB) / &
                  prof(iprof)%cloud(6, lay) / log(10._JPRB)
              ENDIF
              bwyser_ad = 0._JPRB
            ELSE IF (prof(iprof)%idg == 3) THEN
  !Scheme by Boudala et al., 2002, Int. J. Climatol., 22, 1267-1284.
              ztempc = prof(iprof)%t(lev) - rtt
              prof_ad(iprof)%cloud(6, lay) = prof_ad(iprof)%cloud(6, lay) +                               &
                & aux_ad%dg(lay, iprof) * 53.005_JPRB * (prof(iprof)%cloud(6, lay) ** (0.06_JPRB - 1)) *  &
                & exp(0.013_JPRB * ztempc) * 0.06_JPRB
              ztempc_ad                    = ztempc_ad +                                                           &
                & aux_ad%dg(lay, iprof) * 53.005_JPRB * ((prof(iprof)%cloud(6, lay)) ** 0.06_JPRB) * 0.013_JPRB *  &
                & exp(0.013_JPRB * ztempc)
              prof_ad(iprof)%t(lev)        = prof_ad(iprof)%t(lev) + ztempc_ad
              ztempc_ad                    = 0._JPRB
            ELSE IF (prof(iprof)%idg == 4) THEN
  ! Scheme by McFarquhar et al. (2003)
              amcfarq                      = 1.78449_JPRB
              bmcfarq                      = 0.281301_JPRB
              cmcfarq                      = 0.0177166_JPRB
              zmcfarq                      = prof(iprof)%cloud(6, lay)
              aux_ad%fac1_dg(lay, iprof)   = aux_ad%fac1_dg(lay, iprof) + aux_ad%dg(lay, iprof) * 2.0_JPRB
              zmcfarq_ad                   = zmcfarq_ad + aux_ad%fac1_dg(lay, iprof) *                             &
                & 10.0_JPRB ** (amcfarq + bmcfarq * Log10(zmcfarq) + cmcfarq * Log10(zmcfarq) * Log10(zmcfarq)) *  &
                & (bmcfarq + 2._JPRB * cmcfarq * Log10(zmcfarq)) / zmcfarq
              prof_ad(iprof)%cloud(6, lay) = prof_ad(iprof)%cloud(6, lay) + zmcfarq_ad
              zmcfarq_ad                   = 0._jprb
            ELSE
              aux_ad%dg(lay, iprof) = 0._jprb
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDIF
!--------------------------------------------------
! AD for debye terms for pewrmittivity calculations
!--------------------------------------------------
    IF ((coef % id_sensor == sensor_id_mw .OR. coef % id_sensor == sensor_id_po) .AND. &
         opts%rt_mw%clw_data) THEN
      DO lev = prof(iprof)%nlevels, coef%mwcldtop,  - 1
        v = 300.0_JPRB / prof(iprof)%t(lev) - 1.0_JPRB
        aux_ad%debye_prof(3, lev, iprof) =      &
          & aux_ad%debye_prof(3, lev, iprof) + aux_ad%debye_prof(4, lev, iprof) * dcoeff(7)
        v_ad = aux_ad%debye_prof(3, lev, iprof) * dcoeff(5)
        aux_ad%debye_prof(1, lev, iprof) =      &
          & aux_ad%debye_prof(1, lev, iprof) + aux_ad%debye_prof(2, lev, iprof) * dcoeff(4)
        v_ad = v_ad + aux_ad%debye_prof(1, lev, iprof) * ( - dcoeff(2) + 2 * dcoeff(3) * v)
        prof_ad(iprof)%t(lev)            = prof_ad(iprof)%t(lev) + v_ad * ( - 300.0_JPRB / prof(iprof)%t(lev) ** 2)
      ENDDO
!aux_ad % debye_prof(:,:) = 0.
    ENDIF
!-----------------------------------------
! AD for cloud top and surface levels
!-----------------------------------------
    IF (coef%id_sensor /= sensor_id_mw .AND. coef%id_sensor /= sensor_id_po) THEN
!nearest level above cloud top
      dp = prof(iprof)%p(aux%s(iprof)%nearestlev_ctp) - prof(iprof)%p(aux%s(iprof)%nearestlev_ctp - 1)
      dp_ad = 0._jprb
      prof_ad(iprof)%cfraction = prof_ad(iprof)%cfraction + aux_ad%s(iprof)%cfraction
      prof_ad(iprof)%ctp       = prof_ad(iprof)%ctp - aux_ad%s(iprof)%pfraction_ctp / dp
      IF (opts%interpolation%lgradp) THEN
        prof_ad(iprof)%p(aux%s(iprof)%nearestlev_ctp) = prof_ad(iprof)%p(aux%s(iprof)%nearestlev_ctp) +      &
          & aux_ad%s(iprof)%pfraction_ctp / dp
        dp_ad = dp_ad - (prof(iprof)%p(aux%s(iprof)%nearestlev_ctp) - prof(iprof)%ctp) / dp ** 2 *  &
          & aux_ad%s(iprof)%pfraction_ctp
        prof_ad(iprof)%p(aux%s(iprof)%nearestlev_ctp) = prof_ad(iprof)%p(aux%s(iprof)%nearestlev_ctp) + dp_ad
        prof_ad(iprof)%p(aux%s(iprof)%nearestlev_ctp - 1) = prof_ad(iprof)%p(aux%s(iprof)%nearestlev_ctp - 1) - dp_ad
      ENDIF
!Else
! for micro waves do not consider clouds in the RTTOV basis routines
    ENDIF
!nearest level above surface
    dp = prof(iprof)%p(aux%s(iprof)%nearestlev_surf) - prof(iprof)%p(aux%s(iprof)%nearestlev_surf - 1)
    dp_ad = 0._jprb
    prof_ad(iprof)%s2m%p = prof_ad(iprof)%s2m%p - aux_ad%s(iprof)%pfraction_surf / dp
    IF (opts%interpolation%lgradp) THEN
      prof_ad(iprof)%p(aux%s(iprof)%nearestlev_surf) = prof_ad(iprof)%p(aux%s(iprof)%nearestlev_surf) +      &
        & aux_ad%s(iprof)%pfraction_surf / dp
      dp_ad = dp_ad - (prof(iprof)%p(aux%s(iprof)%nearestlev_surf) - prof(iprof)%s2m%p) / dp ** 2 *  &
        & aux_ad%s(iprof)%pfraction_surf
      prof_ad(iprof)%p(aux%s(iprof)%nearestlev_surf) = prof_ad(iprof)%p(aux%s(iprof)%nearestlev_surf) + dp_ad
      prof_ad(iprof)%p(aux%s(iprof)%nearestlev_surf - 1) = prof_ad(iprof)%p(aux%s(iprof)%nearestlev_surf - 1) - dp_ad
    ENDIF
  ENDDO

  IF(aux%on_COEF_levels) THEN
!FWD

    CALL reciprocal(coef%tstar, tstar_r)
    CALL reciprocal(coef%wstar, wstar_r)
    IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) CALL reciprocal(coef%ostar, ostar_r)

!AD
    CALL reciprocal_ad(aux%ww_r, aux_ad%ww, aux_ad%ww_r, acc = .TRUE._jplm)

    aux_ad%wr = aux_ad%wr + aux_ad%wr_sqrt * aux%wr_rsqrt

    aux_ad%wr_rsqrt = &!aux_ad%wr_rsqrt + 
      aux%wr * aux_ad%wr_sqrt

    CALL INVSQRT_AD(aux%wr_rsqrt, aux_ad%wr, aux_ad%wr_rsqrt, acc = .TRUE._jplm)

    CALL sqrt_ad(aux%tw_4rt, aux_ad%tw_sqrt, aux_ad%tw_4rt, acc = .FALSE._jplm)  
    CALL sqrt_ad(aux%tw_sqrt, aux_ad%tw, aux_ad%tw_sqrt, acc = .TRUE._jplm)  
    CALL reciprocal_ad(aux%tr_r, aux_ad%tr, aux_ad%tr_r, acc = .TRUE._jplm)
    
    IF (coef%nozone > 0) THEN
      CALL sqrt_ad(aux%or_sqrt, aux_ad%or, aux_ad%or_sqrt, acc = .TRUE._jplm)  
      aux_ad%ow_rsqrt = &!aux_ad%ow_rsqrt + 
        aux%ow * aux_ad%ow_sqrt

      aux_ad%ow = aux_ad%ow + aux_ad%ow_sqrt * aux%ow_rsqrt
      aux_ad%ow_rsqrt = aux_ad%ow_rsqrt + 2._jprb * aux%ow_rsqrt * aux_ad%ow_r
      CALL INVSQRT_AD(aux%ow_rsqrt, aux_ad%ow, aux_ad%ow_rsqrt, acc = .TRUE._jplm)
    ENDIF

! 4. calculate profile / reference profile sums: tw ww ow

    DO iprof = 1, nprofiles
      DO lay = nlayers, 2, -1
        ! cumulate overlying layers: weighting tr relates to same layer as dpp
        ! do not need dpp(0) to start
        aux_ad%tw(lay - 1, iprof) = aux_ad%tw(lay - 1, iprof) + aux_ad%tw(lay, iprof)
        aux_ad%tr(lay - 1, iprof) = aux_ad%tr(lay - 1, iprof) + coef%dpp(lay - 1) * aux_ad%tw(lay, iprof) 
      ENDDO
    ENDDO

    DO iprof = 1, nprofiles
      sum1 = 0._JPRB
      ! cumulating column overlying layer and layer itself
      DO lay = nlayers, 1, -1
        ! cumulate overlying layers: weighting w or wstar relates to layer below dpp
        ! need dpp(0) to start
        sum1 = sum1 + aux_ad%ww(lay, iprof) * aux%SUM(lay,1)
        aux_ad%w_layer(lay, iprof) = &!aux_ad%w_layer(lay, iprof) + 
          sum1 * coef%dpp(lay - 1)
      ENDDO
    ENDDO
 
    IF (coef%nozone > 0) THEN
      IF (opts%rt_ir%ozone_Data) THEN     ! if no input O3 profile, set to reference value (ow =1)
        DO iprof = 1, nprofiles
          sum1 = 0._JPRB
          DO lay = nlayers, 1, -1
            ! cumulate overlying layers: weighting o or ostar relates to layer below dpp
            ! need dpp(0) to start
            sum1 = sum1 + aux_ad%ow(lay, iprof) * aux%SUM(lay,2)
            aux_ad%o3_layer(lay, iprof) = &!aux_ad%o3_layer(lay, iprof) + 
              sum1 * coef%dpp(lay - 1) ! OK, dpp(0) defined 
          ENDDO
        ENDDO
      ELSE
        aux_ad%o3_layer(:, :) = 0._JPRB
      ENDIF
    ENDIF

! 3. calculate (profile / reference profile) ratios; tr wr or
    DO iprof = 1, nprofiles
      IF ((coef%id_inst == inst_id_ssmis .OR. coef%id_inst == inst_id_ssmisz) .AND. coef%IncZeeman) THEN
        aux_ad%t_layer(:, iprof) = aux_ad%t_layer(:, iprof) + aux_ad%tr(:, iprof) * tstar_r(:)
      ELSE
        aux_ad%t_layer(:, iprof) = aux_ad%tr(:, iprof) * tstar_r(:)
      ENDIF

      aux_ad%w_layer(:, iprof) = aux_ad%w_layer(:, iprof) + aux_ad%wr(:, iprof) * wstar_r(:)
      
! if no input O3 profile, set to reference value (or = 1)

      IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN
        aux_ad%o3_layer(:, iprof) = aux_ad%o3_layer(:, iprof) + aux_ad%or(:, iprof)  * ostar_r(:)
!      ELSE
!        aux_ad%or(:, iprof) = 0._JPRB
      ENDIF

! 2. calculate, for layers, deviations from reference profile
      aux_ad%t_layer(:, iprof) = aux_ad%t_layer(:, iprof) + aux_ad%dt(:, iprof) 

      ! if no input O3 profile we still use the input temperature profile for dto
      IF (coef%nozone > 0) &
        aux_ad%t_layer(:, iprof) = aux_ad%t_layer(:, iprof) + aux_ad%dto(:, iprof) 

! 1 profile layer quantities
!   the layer number agrees with the level number of its upper boundary
! layer N-1 lies between levels N-1 and N
    !DAR add from set_predictors_7
      IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN
        prof_ad(iprof)%o3(1) = prof_ad(iprof)%o3(1) + &
          0.5_JPRB * aux_ad%o3_layer(1, iprof)
        prof_ad(iprof)%o3(2:nlevels-1) = prof_ad(iprof)%o3(2:nlevels-1) + &
          0.5_JPRB * (aux_ad%o3_layer(1:nlevels-2, iprof) + aux_ad%o3_layer(2:nlevels-1, iprof))
        prof_ad(iprof)%o3(nlevels) = prof_ad(iprof)%o3(nlevels) + &
          0.5_JPRB * aux_ad%o3_layer(nlevels-1, iprof)
      ENDIF

      prof_ad(iprof)%t(1) = prof_ad(iprof)%t(1) + &
        0.5_JPRB * aux_ad%t_layer(1, iprof)
      prof_ad(iprof)%t(2:nlevels-1) = prof_ad(iprof)%t(2:nlevels-1) + &
        0.5_jprb * (aux_ad%t_layer(1:nlevels-2, iprof) + aux_ad%t_layer(2:nlevels-1, iprof))
      prof_ad(iprof)%t(nlevels) = prof_ad(iprof)%t(nlevels) + &
        0.5_JPRB * aux_ad%t_layer(nlevels-1, iprof)

      IF (opts%rt_all%use_q2m) THEN
        iv3lev = aux%s(iprof)%nearestlev_surf - 1! nearest level above surface
        iv2lev = aux%s(iprof)%nearestlev_surf    ! nearest level above surface

        IF (iv2lev <= coef%nlevels) THEN
          iv2lay       = iv2lev - 1
          prof_ad(iprof)%s2m%q =  prof_ad(iprof)%s2m%q + &
            aux_ad%w_layer(iv2lay, iprof) * 0.5_JPRB
          
          prof_ad(iprof)%q(iv2lev) = prof_ad(iprof)%q(iv2lev) - &
            aux_ad%w_layer(iv2lay, iprof) * 0.5_JPRB
        ENDIF
      ENDIF

      prof_ad(iprof)%q(1) = prof_ad(iprof)%q(1) + &
        0.5_JPRB * aux_ad%w_layer(1, iprof)
      prof_ad(iprof)%q(2:nlevels-1) = prof_ad(iprof)%q(2:nlevels-1) + &
        0.5_JPRB * (aux_ad%w_layer(1:nlevels-2, iprof) + aux_ad%w_layer(2:nlevels-1, iprof))
      prof_ad(iprof)%q(nlevels) = prof_ad(iprof)%q(nlevels) + &
        0.5_JPRB * aux_ad%w_layer(nlevels-1, iprof)

!       IF (opts%rt_all%use_q2m) THEN
!         iv3lev = aux%s(iprof)%nearestlev_surf - 1! nearest level above surface
!         iv2lev = aux%s(iprof)%nearestlev_surf    ! nearest level above surface

!         ! IF (iv2lev <= coef%nlevels) THEN
!         !   iv2lay       = iv2lev - 1
!         !   prof_ad(iprof)%s2m%q = prof_ad(iprof)%s2m%q + &
!         !     aux_ad%w_layer(iv2lay, iprof) * 0.5_JPRB

!         !   prof_ad(iprof)%q(iv3lev) = prof_ad(iprof)%q(iv3lev) + &
!         !      aux_ad%w_layer(iv2lay, iprof) * 0.5_JPRB
!         ! ENDIF

! !        IF (lev == iv2lev) THEN
!           iv2lay = iv2lev - 1
!           prof_ad(iprof)%s2m%q = prof_ad(iprof)%s2m%q + 0.5_JPRB * aux_ad%w_layer(iv2lay,iprof)
! !        ELSE
!           prof_ad(iprof)%q(iv2lev) = prof_ad(iprof)%q(iv2lev) - 0.5_JPRB * aux_ad%w_layer(iv2lay,iprof)
! !        ENDIF


!       ENDIF
    ENDDO
ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_PROFAUX_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_profaux_ad
