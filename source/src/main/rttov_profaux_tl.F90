SUBROUTINE rttov_profaux_tl( &
            & opts,    &
            & prof,    &
            & prof_tl, &
            & coef,    &
            & aux,     &
            & aux_tl)!
! Description:
! TL of rttov_profaux
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
!  1.0    07/10/2004  Added history
!  1.1    29/03/2005  Add end of header comment (J. Cameron)
!  1.2    04/02/2008  opts%lgradp option for TL/AD of pressure levels (N. Bormann)
!  1.3    15/07/2009  User defined ToA. Layers distinct from levels (P.Rayer)
!  1.4    02/12/2009  Fixed a number of bugs due to the wrong assumption that cloud
!                     related quantities are defined on levels (thay are layer
!                     average quantities). Marco Matricardi
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
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
       & dgmin_hex, dgmax_hex, dgmin_agg, dgmax_agg 
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
  USE rttov_math_mod, ONLY: reciprocal, &
                            INVSQRT_TL, reciprocal_tl, sqrt_tl
!INTF_ON

  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options), INTENT(IN)    :: opts
  TYPE(profile_Type ), INTENT(IN)    :: prof   (:)
  TYPE(profile_Type ), INTENT(IN)    :: prof_tl(:)
  TYPE(rttov_coef   ), INTENT(IN)    :: coef
  TYPE(profile_aux  ), INTENT(IN)    :: aux
  TYPE(profile_aux  ), INTENT(INOUT) :: aux_tl

!INTF_END
! local
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: iprof
  INTEGER(KIND=jpim) :: lev      , lay
  REAL   (KIND=jprb) :: dp
  REAL   (KIND=jprb) :: dp_tl
  REAL   (KIND=jprb) :: v                                                       ! temperature ratio
  REAL   (KIND=jprb) :: v_tl                                                    ! temperature ratio
  REAL   (KIND=jprb) :: amcfarq     , bmcfarq     , cmcfarq, zmcfarq, zmcfarq_tl
  REAL   (KIND=jprb) :: ztempc      , ztempc_tl
  REAL   (KIND=jprb) :: zradipou_upp, zradipou_low
  REAL   (KIND=jprb) :: bwyser      , bwyser_tl   , nft
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
!-----------------------------------------
! TL for cloud top and surface levels
!-----------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_PROFAUX_TL', 0_jpim, ZHOOK_HANDLE)
  nprofiles = size(prof)
  nlayers = prof(1)%nlayers
  nlevels = prof(1)%nlayers + 1

  DO iprof = 1, nprofiles
!nearest level above surface
    dp = prof(iprof)%p(aux%s(iprof)%nearestlev_surf) - prof(iprof)%p(aux%s(iprof)%nearestlev_surf - 1)
    IF (opts%interpolation%lgradp) THEN
      dp_tl = prof_tl(iprof)%p(aux%s(iprof)%nearestlev_surf) - prof_tl(iprof)%p(aux%s(iprof)%nearestlev_surf - 1)
      aux_tl%s(iprof)%pfraction_surf = &
        & (prof_tl(iprof)%p(aux%s(iprof)%nearestlev_surf) - prof_tl(iprof)%s2m%p) / dp -  &
        & (prof(iprof)%p(aux%s(iprof)%nearestlev_surf) - prof(iprof)%s2m%p) / dp ** 2 * dp_tl
    ELSE
      aux_tl%s(iprof)%pfraction_surf =  - prof_tl(iprof)%s2m%p / dp
    ENDIF
    IF (coef%id_sensor /= sensor_id_mw .AND. coef%id_sensor /= sensor_id_po) THEN
!nearest level above cloud top
      dp = prof(iprof)%p(aux%s(iprof)%nearestlev_ctp) - prof(iprof)%p(aux%s(iprof)%nearestlev_ctp - 1)
      IF (opts%interpolation%lgradp) THEN
        dp_tl = prof_tl(iprof)%p(aux%s(iprof)%nearestlev_ctp) - prof_tl(iprof)%p(aux%s(iprof)%nearestlev_ctp - 1)
        aux_tl%s(iprof)%pfraction_ctp = &
          & (prof_tl(iprof)%p(aux%s(iprof)%nearestlev_ctp) - prof_tl(iprof)%ctp) / dp -  &
          & (prof(iprof)%p(aux%s(iprof)%nearestlev_ctp) - prof(iprof)%ctp) / dp ** 2 * dp_tl
      ELSE
        aux_tl%s(iprof)%pfraction_ctp =  - prof_tl(iprof)%ctp / dp
      ENDIF
      aux_tl%s(iprof)%cfraction = prof_tl(iprof)%cfraction
    ELSE
! for micro waves do not consider clouds in the RTTOV basis routines
      aux_tl%s(iprof)%pfraction_ctp = 0._JPRB
      aux_tl%s(iprof)%cfraction     = 0._JPRB
    ENDIF
!--------------------------------------------------
! TL for debye terms for pewrmittivity calculations
!--------------------------------------------------
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
      DO lev = coef%mwcldtop, prof(iprof)%nlevels
        v = 300.0_JPRB / prof(iprof)%t(lev) - 1.0_JPRB
        v_tl =  - 300.0_JPRB * prof_tl(iprof)%t(lev) / prof(iprof)%t(lev) ** 2
        aux_tl%debye_prof(1, lev, iprof) =  - dcoeff(2) * v_tl + 2 * dcoeff(3) * v_tl * v
        aux_tl%debye_prof(2, lev, iprof) = dcoeff(4) * aux_tl%debye_prof(1, lev, iprof)
        aux_tl%debye_prof(3, lev, iprof) = dcoeff(5) * v_tl
        aux_tl%debye_prof(4, lev, iprof) = dcoeff(7) * aux_tl%debye_prof(3, lev, iprof)
        aux_tl%debye_prof(5, lev, iprof) = 0._JPRB
      ENDDO
    ENDIF
  ENDIF
!--------------------------------------------------------------
! TL for ice  water content into effective generalized diameter
!--------------------------------------------------------------
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
            aux_tl%dg(lay, iprof) = prof_tl(iprof)%icede(lay)
          ELSE
            IF (prof(iprof)%idg == 1) THEN
  !Scheme by Ou and Liou, 1995, Atmos. Res., 35, 127-138.
              ztempc = prof(iprof)%t(lev) - rtt
              ztempc_tl                  = prof_tl(iprof)%t(lev)
              aux_tl%fac1_dg(lay, iprof) =      &
                & ztempc_tl * (12.42_JPRB + 2._JPRB * 0.197_JPRB * ztempc + 3._JPRB * 0.0012_JPRB * ztempc * ztempc)
              aux_tl%fac2_dg(lay, iprof) = 0.388_JPRB * aux_tl%fac1_dg(lay, iprof) +      &
                & 0.00051_JPRB * 2 * aux%fac1_dg(lay, iprof) * aux_tl%fac1_dg(lay, iprof)
              aux_tl%fac3_dg(lay, iprof) = 2.0_JPRB * aux_tl%fac2_dg(lay, iprof)
  !
  ! Take Ou-Liou scheme as being valid only between 20C and 60C
  !
              IF (aux%fac3_dg(lay, iprof) < zradipou_low .OR. aux%fac3_dg(lay, iprof) > zradipou_upp) THEN
                aux_tl%dg(lay, iprof) = 0._JPRB
              ELSE
                aux_tl%dg(lay, iprof) = aux_tl%fac3_dg(lay, iprof)
              ENDIF
            ELSE IF (prof(iprof)%idg == 2) THEN
  !Scheme by Wyser et al. (see McFarquhar et al. (2003))
              bwyser    =  - 2.0_JPRB
              bwyser_tl = 0._JPRB
              IF (prof(iprof)%t(lev) < 273._JPRB) THEN
                bwyser    = bwyser + (0.001_JPRB * &
                  ((273._JPRB - prof(iprof)%t(lev)) ** 1.5_JPRB) * Log10(prof(iprof)%cloud(6, lay) / 50._JPRB))
                bwyser_tl =  -                                                                                        &
                  & 0.001_JPRB * 1.5_JPRB * ((273._JPRB - prof(iprof)%t(lev)) ** 0.5_JPRB) * prof_tl(iprof)%t(lev) *  &
                  & Log10(prof(iprof)%cloud(6, lay) / 50._JPRB) +                                                     &
                  & 0.001_JPRB * ((273._JPRB - prof(iprof)%t(lev)) ** 1.5_JPRB) * prof_tl(iprof)%cloud(6, lay) /      &
                  & prof(iprof)%cloud(6, lay) / log(10._JPRB)
              ENDIF
              aux_tl%fac1_dg(lay, iprof) =      &
                & (203.3_JPRB + 2 * bwyser * 37.91_JPRB + 3 * bwyser * bwyser * 2.3696_JPRB) * bwyser_tl
              nft = (sqrt(3._JPRB) + 4._JPRB) / (3._JPRB * sqrt(3._JPRB))
              aux_tl%fac2_dg(lay, iprof) = aux_tl%fac1_dg(lay, iprof) / nft
              aux_tl%dg(lay, iprof)      = 2._JPRB * 4._JPRB * aux_tl%fac2_dg(lay, iprof) * sqrt(3._JPRB) / 9._JPRB
            ELSE IF (prof(iprof)%idg == 3) THEN
  !Scheme by Boudala et al., 2002, Int. J. Climatol., 22, 1267-1284.
              ztempc                = prof(iprof)%t(lev) - rtt
              ztempc_tl             = prof_tl(iprof)%t(lev)
              aux_tl%dg(lay, iprof) = 53.005_JPRB * 0.06_JPRB *                                 &
                (prof(iprof)%cloud(6, lay) ** (0.06_JPRB - 1)) * prof_tl(iprof)%cloud(6, lay) * &
                exp(0.013_JPRB * ztempc) + 53.005_JPRB *                                        &
                ((prof(iprof)%cloud(6, lay)) ** 0.06_JPRB) * 0.013_JPRB * ztempc_tl * exp(0.013_JPRB * ztempc)
            ELSE IF (prof(iprof)%idg == 4) THEN
  ! Scheme by McFarquhar et al. (2003)
              amcfarq                    = 1.78449_JPRB
              bmcfarq                    = 0.281301_JPRB
              cmcfarq                    = 0.0177166_JPRB
              zmcfarq                    = prof(iprof)%cloud(6, lay)
              zmcfarq_tl                 = prof_tl(iprof)%cloud(6, lay)
              aux_tl%fac1_dg(lay, iprof) =                                                                         &
                & 10.0_JPRB ** (amcfarq + bmcfarq * Log10(zmcfarq) + cmcfarq * Log10(zmcfarq) * Log10(zmcfarq)) *  &
                & (bmcfarq + 2._JPRB * cmcfarq * Log10(zmcfarq)) / zmcfarq * zmcfarq_tl
              aux_tl%dg(lay, iprof)      = 2.0_JPRB * aux_tl%fac1_dg(lay, iprof)
            ELSE
              aux_tl%dg(lay, iprof) = 0._jprb
            ENDIF
          ENDIF
  ! test and force effective size betwwen values from Table 14 and Table 15 of RTTOV-9 SVR   
          IF ( ABS(aux%dg(lay, iprof)) > 0._JPRB ) THEN
            IF (prof(iprof)%ish == 1_JPIM) THEN
              IF ( aux%dg(lay, iprof) < dgmin_hex ) aux_tl%dg(lay, iprof) = 0._JPRB
              IF ( aux%dg(lay, iprof) > dgmax_hex ) aux_tl%dg(lay, iprof) = 0._JPRB
            ENDIF
            IF (prof(iprof)%ish == 2_JPIM) THEN
              IF ( aux%dg(lay, iprof) < dgmin_agg ) aux_tl%dg(lay, iprof) = 0._JPRB
              IF ( aux%dg(lay, iprof) > dgmax_agg ) aux_tl%dg(lay, iprof) = 0._JPRB
            ENDIF
          ENDIF          
        ENDIF
      ENDDO
    ENDIF
  ENDDO

  if(aux%on_COEF_levels) THEN
    ! 1 profile layer quantities
!   the layer number agrees with the level number of its upper boundary
! layer N-1 lies between levels N-1 and N
    CALL reciprocal(coef%tstar, tstar_r)
    CALL reciprocal(coef%wstar, wstar_r)
    IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) CALL reciprocal(coef%ostar, ostar_r)

    !DAR add from set_predictors_7
    DO iprof = 1, nprofiles
      aux_tl%t_layer(1:nlayers, iprof) = &
        (prof_tl(iprof)%t(1:nlevels-1) + prof_tl(iprof)%t(2:nlevels)) * 0.5_JPRB
      aux_tl%w_layer(1:nlayers, iprof) = &
        (prof_tl(iprof)%q(1:nlevels-1) + prof_tl(iprof)%q(2:nlevels)) * 0.5_JPRB

      IF (opts%rt_all%use_q2m) THEN
        iv3lev = aux%s(iprof)%nearestlev_surf - 1! nearest level above surface
        iv2lev = aux%s(iprof)%nearestlev_surf    ! nearest level above surface

        IF (iv2lev <= coef%nlevels) THEN
          iv2lay       = iv2lev - 1
          aux_tl%w_layer(iv2lay, iprof) = &
            (prof_tl(iprof)%s2m%q + prof_tl(iprof)%q(iv3lev)) * 0.5_JPRB
        ENDIF
      ENDIF

      IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN
!        aux_tl%o3_layer(layer, iprof) = (prof(iprof)%o3(level - 1) + prof(iprof)%o3(level)) * 0.5_JPRB
        aux_tl%o3_layer(1:nlayers, iprof) = &
          (prof_tl(iprof)%o3(1:nlevels-1) + prof_tl(iprof)%o3(2:nlevels)) * 0.5_JPRB
      ENDIF

! 2. calculate, for layers, deviations from reference profile
      aux_tl%dt(:, iprof) = aux_tl%t_layer(:, iprof)

      ! if no input O3 profile we still use the input temperature profile for dto
      IF (coef%nozone > 0) &
        aux_tl%dto(:, iprof) = aux_tl%t_layer(:, iprof)

! 3. calculate (profile / reference profile) ratios; tr wr or
      aux_tl%tr(:, iprof) = aux_tl%t_layer(:, iprof) * tstar_r(:)
      aux_tl%wr(:, iprof) = aux_tl%w_layer(:, iprof) * wstar_r(:)
      
! if no input O3 profile, set to reference value (or = 1)

      IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN
        aux_tl%or(:, iprof) = aux_tl%o3_layer(:, iprof) * ostar_r(:)
      ELSE
        aux_tl%or(:, iprof) = 0._JPRB
      ENDIF
    ENDDO

! 4. calculate profile / reference profile sums: tw ww ow

    DO iprof = 1, nprofiles
      aux_tl%tw(1, iprof) = 0._JPRB
      DO lay = 2, nlayers
        ! cumulate overlying layers: weighting tr relates to same layer as dpp
        ! do not need dpp(0) to start
        aux_tl%tw(lay, iprof) = aux_tl%tw(lay - 1, iprof) + &
          coef%dpp(lay - 1) * aux_tl%tr(lay - 1, iprof)
      ENDDO
    ENDDO

    DO iprof = 1, nprofiles
      sum1 = 0._JPRB
      ! cumulating column overlying layer and layer itself
      DO lay = 1, nlayers
        ! cumulate overlying layers: weighting w or wstar relates to layer below dpp
        ! need dpp(0) to start
        sum1 = sum1 + coef%dpp(lay - 1) * aux_tl%w_layer(lay, iprof)
        aux_tl%ww(lay, iprof) = sum1 * aux%sum(lay,1)
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
            sum1 = sum1 + coef%dpp(lay - 1) * aux_tl%o3_layer(lay, iprof) ! OK, dpp(0) defined 
            aux_tl%ow(lay, iprof) = sum1 * aux%sum(lay,2)
          ENDDO
        ENDDO
      ELSE
        aux_tl%ow(:, :) = 0._JPRB
      ENDIF
    ENDIF

    CALL reciprocal_tl(aux%tr_r, aux_tl%tr, aux_tl%tr_r)
    
    CALL sqrt_tl(aux%tw_sqrt, aux_tl%tw, aux_tl%tw_sqrt)  
    CALL sqrt_tl(aux%tw_4rt, aux_tl%tw_sqrt, aux_tl%tw_4rt)  

    CALL INVSQRT_TL(aux%wr_rsqrt, aux_tl%wr, aux_tl%wr_rsqrt)
    aux_tl%wr_sqrt = aux_tl%wr * aux%wr_rsqrt + &
                     aux%wr * aux_tl%wr_rsqrt

    CALL reciprocal_tl(aux%ww_r, aux_tl%ww, aux_tl%ww_r)

    IF (coef%nozone > 0) THEN
      CALL INVSQRT_TL(aux%ow_rsqrt, aux_tl%ow, aux_tl%ow_rsqrt)
      aux_tl%ow_sqrt = aux_tl%ow * aux%ow_rsqrt + &
                       aux%ow * aux_tl%ow_rsqrt
      aux_tl%ow_r = 2._jprb * aux%ow_rsqrt * aux_tl%ow_rsqrt

      CALL sqrt_tl(aux%or_sqrt, aux_tl%or, aux_tl%or_sqrt)  
    ENDIF   
  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_PROFAUX_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_profaux_tl
