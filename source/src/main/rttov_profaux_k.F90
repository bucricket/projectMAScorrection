SUBROUTINE rttov_profaux_k( &
            & opts,    &
            & chanprof,&
            & profiles,    &
            & profiles_k, &
            & coef,    &
            & aux,     &
            & aux_k)
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
       & rttov_chanprof, &
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
                            INVSQRT_K, reciprocal_k, sqrt_k
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options), INTENT(IN)    :: opts
  TYPE(rttov_chanprof), INTENT(IN)   :: chanprof(:)
  TYPE(profile_Type ), INTENT(IN)    :: profiles(:)
  TYPE(profile_Type ), INTENT(INOUT) :: profiles_k(:)
  TYPE(rttov_coef   ), INTENT(IN)    :: coef
  TYPE(profile_aux  ), INTENT(IN)    :: aux
  TYPE(profile_aux  ), INTENT(INOUT) :: aux_k

!INTF_END
! local
  INTEGER(KIND=jpim) :: nprofiles, nchannels
  INTEGER(KIND=jpim) :: i, prof
  INTEGER(KIND=jpim) :: lev      , lay
  REAL   (KIND=jprb) :: dp
  REAL   (KIND=jprb) :: dp_k
  REAL   (KIND=jprb) :: v                                                       ! temperature ratio
  REAL   (KIND=jprb) :: v_k                                                    ! temperature ratio
  REAL   (KIND=jprb) :: amcfarq     , bmcfarq     , cmcfarq, zmcfarq, zmcfarq_k
  REAL   (KIND=jprb) :: ztempc      , ztempc_k
  REAL   (KIND=jprb) :: zradipou_upp, zradipou_low
  REAL   (KIND=jprb) :: bwyser      , bwyser_k   , nft
  REAL(KIND=jprb), PARAMETER :: rtt      = 273.15_JPRB
  REAL(KIND=jprb), PARAMETER :: rtou_upp =  - 20._JPRB
  REAL(KIND=jprb), PARAMETER :: rtou_low =  - 60._JPRB
  REAL(KIND=JPRB) :: ZHOOK_HANDLE

  REAL(jprb) :: tstar_r(profiles(1)%nlayers), ostar_r(profiles(1)%nlayers)
  REAL(jprb) :: wstar_r(profiles(1)%nlayers)
  REAL(jprb) :: sum1
  INTEGER(jpim) :: iv2lay, iv2lev, iv3lev
  INTEGER(jpim) :: nlayers, nlevels
  INTEGER(jpim) :: map(SIZE(chanprof),2), prof_stat

!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_PROFAUX_K', 0_jpim, ZHOOK_HANDLE)
  nprofiles = size(profiles)
  nchannels = size(chanprof)
  nlayers = profiles(1)%nlayers
  nlevels = profiles(1)%nlayers + 1

  map(1,1) = chanprof(1)%prof
  map(1,2) = chanprof(1)%chan

  prof_stat = 1 ! assume profs are contiguous and monotonic
  DO i = 2, nchannels
    map(i,1) = chanprof(i)%prof
    map(i,2) = chanprof(i)%chan

    IF(map(i,1) < map(i-1,1)) THEN ! they are not.
      prof_stat = -1
    ENDIF
  ENDDO
  zmcfarq_k = 0._jprb
  ztempc_k  = 0._jprb
  bwyser_k  = 0._jprb
!------------------------------------------------------------
! K for ice water content into effective generalized diameter
!------------------------------------------------------------
  IF (opts%rt_ir%addclouds .AND. .NOT. opts%rt_ir%user_cld_opt_param) THEN
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
    DO i = 1, nchannels
      prof = chanprof(i)%prof
      IF (profiles(prof)%ish >= 3) CYCLE
      DO lay = profiles(prof)%nlayers, 1,  - 1
        lev = lay + 1
        IF (profiles(prof)%cloud(6, lay) > 0._jprb) THEN
          IF ( ABS(aux%dg(lay, prof)) > 0._JPRB ) THEN
            IF (profiles(prof)%ish == 1_JPIM) THEN
              IF ( aux%dg(lay, prof) < dgmin_hex ) THEN
                aux_k%dg(lay, i) = 0._JPRB
              ENDIF
              IF ( aux%dg(lay, prof) > dgmax_hex ) THEN
                aux_k%dg(lay, i) = 0._JPRB
              ENDIF
            ENDIF
            IF (profiles(prof)%ish == 2_JPIM) THEN
              IF ( aux%dg(lay, prof) < dgmin_agg ) THEN
                aux_k%dg(lay, i) = 0._JPRB
              ENDIF
              IF ( aux%dg(lay, prof) > dgmax_agg ) THEN
                aux_k%dg(lay, i) = 0._JPRB
              ENDIF
            ENDIF
          ENDIF
          ! Use effective diameter from input profile if specified
          IF (profiles(prof)%icede(lay) > 0._jprb) THEN
            profiles_k(i)%icede(lay) = profiles_k(i)%icede(lay) + aux_k%dg(lay, i)
          ELSE
            IF (profiles(prof)%idg == 1) THEN
  !Scheme by Ou and Liou, 1995, Atmos. Res., 35, 127-138.
              ztempc = profiles(prof)%t(lev) - rtt
  !
  ! Take Ou-Liou scheme as being valid only between 20C and 60C
  !
              IF (aux%fac3_dg(lay, prof) < zradipou_low .OR. aux%fac3_dg(lay, prof) > zradipou_upp) THEN
                aux_k%dg(lay, i) = 0._JPRB
              ELSE
                aux_k%fac3_dg(lay, i) = aux_k%fac3_dg(lay, i) + aux_k%dg(lay, i)
              ENDIF
              aux_k%fac2_dg(lay, i) = aux_k%fac2_dg(lay, i) + 2.0_JPRB * aux_k%fac3_dg(lay, i)
              aux_k%fac1_dg(lay, i) = aux_k%fac1_dg(lay, i) + 0.388_JPRB * aux_k%fac2_dg(lay, i)
              aux_k%fac1_dg(lay, i) =      &
                & aux_k%fac1_dg(lay, i) + aux_k%fac2_dg(lay, i) * 0.00051_JPRB * 2 * aux%fac1_dg(lay, prof)
              ztempc_k                   = ztempc_k + aux_k%fac1_dg(lay, i) *      &
                & (12.42_JPRB + 2._JPRB * 0.197_JPRB * ztempc + 3._JPRB * 0.0012_JPRB * ztempc * ztempc)
              profiles_k(i)%t(lev)       = profiles_k(i)%t(lev) + ztempc_k
              ztempc_k                   = 0._jprb
            ELSE IF (profiles(prof)%idg == 2) THEN
  !Scheme by Wyser et al. (see McFarquhar et al. (2003))
              bwyser =  - 2.0_JPRB
              IF (profiles(prof)%t(lev) < 273._JPRB) THEN
                bwyser = bwyser + (0.001_JPRB * &
                  ((273._JPRB - profiles(prof)%t(lev)) ** 1.5_JPRB) * Log10(profiles(prof)%cloud(6, lay) / 50._JPRB))
              ENDIF
              nft = (sqrt(3._JPRB) + 4._JPRB) / (3._JPRB * sqrt(3._JPRB))
              aux_k%fac2_dg(lay, i) =      &
                & aux_k%fac2_dg(lay, i) + aux_k%dg(lay, i) * 2._JPRB * 4._JPRB * sqrt(3._JPRB) / 9._JPRB
              aux_k%fac1_dg(lay, i) = aux_k%fac1_dg(lay, i) + aux_k%fac2_dg(lay, i) / nft
              bwyser_k                   = bwyser_k +      &
                & aux_k%fac1_dg(lay, i) * (203.3_JPRB + 2 * bwyser * 37.91_JPRB + 3 * bwyser * bwyser * 2.3696_JPRB)
              IF (profiles(prof)%t(lev) < 273._JPRB) THEN
                profiles_k(i)%t(lev)        = profiles_k(i)%t(lev) -                                     &
                  & bwyser_k * 0.001_JPRB * 1.5_JPRB * ((273._JPRB - profiles(prof)%t(lev)) ** 0.5_JPRB) *  &
                  & Log10(profiles(prof)%cloud(6, lay) / 50._JPRB)
                profiles_k(i)%cloud(6, lay) = profiles_k(i)%cloud(6, lay) +                                               &
                  & bwyser_k * 0.001_JPRB * ((273._JPRB - profiles(prof)%t(lev)) ** 1.5_JPRB) / profiles(prof)%cloud(6, lay) /  &
                  & log(10._JPRB)
              ENDIF
              bwyser_k = 0._JPRB
            ELSE IF (profiles(prof)%idg == 3) THEN
  !Scheme by Boudala et al., 2002, Int. J. Climatol., 22, 1267-1284.
              ztempc = profiles(prof)%t(lev) - rtt
              profiles_k(i)%cloud(6, lay) = profiles_k(i)%cloud(6, lay) +                                 &
                & aux_k%dg(lay, i) * 53.005_JPRB * (profiles(prof)%cloud(6, lay) ** (0.06_JPRB - 1)) *  &
                & exp(0.013_JPRB * ztempc) * 0.06_JPRB
              ztempc_k                    = ztempc_k +                                                             &
                & aux_k%dg(lay, i) * 53.005_JPRB * ((profiles(prof)%cloud(6, lay)) ** 0.06_JPRB) * 0.013_JPRB *  &
                & exp(0.013_JPRB * ztempc)
              profiles_k(i)%t(lev)        = profiles_k(i)%t(lev) + ztempc_k
              ztempc_k                    = 0._JPRB
            ELSE IF (profiles(prof)%idg == 4) THEN
  ! Scheme by McFarquhar et al. (2003)
              amcfarq                     = 1.78449_JPRB
              bmcfarq                     = 0.281301_JPRB
              cmcfarq                     = 0.0177166_JPRB
              zmcfarq                     = profiles(prof)%cloud(6, lay)
              aux_k%fac1_dg(lay, i)  = aux_k%fac1_dg(lay, i) + aux_k%dg(lay, i) * 2.0_JPRB
              zmcfarq_k                   = zmcfarq_k + aux_k%fac1_dg(lay, i) *                               &
                & 10.0_JPRB ** (amcfarq + bmcfarq * Log10(zmcfarq) + cmcfarq * Log10(zmcfarq) * Log10(zmcfarq)) *  &
                & (bmcfarq + 2._JPRB * cmcfarq * Log10(zmcfarq)) / zmcfarq
              profiles_k(i)%cloud(6, lay) = profiles_k(i)%cloud(6, lay) + zmcfarq_k
              zmcfarq_k                   = 0._jprb
            ELSE
              aux_k%dg(lay, i) = 0._jprb
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ENDIF
!-------------------------------------------------
! K for debye terms for pewrmittivity calculations
!-------------------------------------------------
  DO i = 1, nchannels
    prof = chanprof(i)%prof
    IF (opts%interpolation%lgradp) dp_k = 0._jprb
    IF ((coef % id_sensor == sensor_id_mw .or. coef % id_sensor == sensor_id_po) .and. &
         opts%rt_mw%clw_data) THEN
      DO lev = profiles(prof)%nlevels, coef%mwcldtop,  - 1
        v = 300.0_JPRB / profiles(prof)%t(lev) - 1.0_JPRB
        aux_k%debye_prof(3, lev, i) =      &
          & aux_k%debye_prof(3, lev, i) + aux_k%debye_prof(4, lev, i) * dcoeff(7)
        v_k = aux_k%debye_prof(3, lev, i) * dcoeff(5)
        aux_k%debye_prof(1, lev, i) =      &
          & aux_k%debye_prof(1, lev, i) + aux_k%debye_prof(2, lev, i) * dcoeff(4)
        v_k = v_k + aux_k%debye_prof(1, lev, i) * ( - dcoeff(2) + 2 * dcoeff(3) * v)
        profiles_k(i)%t(lev)             = profiles_k(i)%t(lev) + v_k * ( - 300.0_JPRB / profiles(prof)%t(lev) ** 2)
      ENDDO
!auxiles_k(i) % debye_prof(:,:) = 0.
    ENDIF
!-----------------------------------
! K for cloud top and surface levels
!-----------------------------------
    IF (coef%id_sensor /= sensor_id_mw .AND. coef%id_sensor /= sensor_id_po) THEN
!nearest level above cloud top
      dp = profiles(prof)%p(aux%s(prof)%nearestlev_ctp) - profiles(prof)%p(aux%s(prof)%nearestlev_ctp - 1)
      dp_k = 0._jprb
      profiles_k(i)%cfraction = profiles_k(i)%cfraction + aux_k%s(i)%cfraction
      profiles_k(i)%ctp       = profiles_k(i)%ctp - aux_k%s(i)%pfraction_ctp / dp
      IF (opts%interpolation%lgradp) THEN
        profiles_k(i)%p(aux%s(prof)%nearestlev_ctp) = profiles_k(i)%p(aux%s(prof)%nearestlev_ctp) +      &
          & aux_k%s(i)%pfraction_ctp / dp
        dp_k = dp_k - (profiles(prof)%p(aux%s(prof)%nearestlev_ctp) - profiles(prof)%ctp) / dp ** 2 *  &
          & aux_k%s(i)%pfraction_ctp
        profiles_k(i)%p(aux%s(prof)%nearestlev_ctp) = profiles_k(i)%p(aux%s(prof)%nearestlev_ctp) + dp_k
        profiles_k(i)%p(aux%s(prof)%nearestlev_ctp - 1) = profiles_k(i)%p(aux%s(prof)%nearestlev_ctp - 1) - dp_k
      ENDIF
!Else
! for micro waves do not consider clouds in the RTTOV basis routines
    ENDIF
!nearest level above surface
    dp = profiles(prof)%p(aux%s(prof)%nearestlev_surf) - profiles(prof)%p(aux%s(prof)%nearestlev_surf - 1)
    dp_k = 0._jprb
    profiles_k(i)%s2m%p = profiles_k(i)%s2m%p - aux_k%s(i)%pfraction_surf / dp
    IF (opts%interpolation%lgradp) THEN
      profiles_k(i)%p(aux%s(prof)%nearestlev_surf) = profiles_k(i)%p(aux%s(prof)%nearestlev_surf) +      &
        & aux_k%s(i)%pfraction_surf / dp
      dp_k = dp_k - (profiles(prof)%p(aux%s(prof)%nearestlev_surf) - profiles(prof)%s2m%p) / dp ** 2 *  &
        & aux_k%s(i)%pfraction_surf
      profiles_k(i)%p(aux%s(prof)%nearestlev_surf) = profiles_k(i)%p(aux%s(prof)%nearestlev_surf) + dp_k
      profiles_k(i)%p(aux%s(prof)%nearestlev_surf - 1) = profiles_k(i)%p(aux%s(prof)%nearestlev_surf - 1) - dp_k
    ENDIF
  ENDDO
! Channels loop

  IF(aux%on_COEF_levels) THEN

!FWD
    CALL reciprocal(coef%tstar, tstar_r)
    CALL reciprocal(coef%wstar, wstar_r)
    IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) CALL reciprocal(coef%ostar, ostar_r)
!AD
    CALL reciprocal_k(aux%ww_r, aux_k%ww, aux_k%ww_r, &
      acc = .TRUE._jplm, map = map, prof_stat = prof_stat)
    DO i = 1, nchannels
      prof = chanprof(i)%prof
      
      aux_k%wr(:,i) = aux_k%wr(:,i) + aux_k%wr_sqrt(:,i) * aux%wr_rsqrt(:,prof)

      aux_k%wr_rsqrt(:,i) = &!aux_k%wr_rsqrt + 
        aux%wr(:,prof) * aux_k%wr_sqrt(:,i)
    ENDDO

    CALL INVSQRT_K(aux%wr_rsqrt, aux_k%wr, aux_k%wr_rsqrt, &
      acc = .TRUE._jplm, map = map, prof_stat = prof_stat)

    CALL sqrt_k(aux%tw_4rt, aux_k%tw_sqrt, aux_k%tw_4rt, &
      acc = .FALSE._jplm, map = map, prof_stat = prof_stat)  
    CALL sqrt_k(aux%tw_sqrt, aux_k%tw, aux_k%tw_sqrt, &
      acc = .TRUE._jplm, map = map, prof_stat = prof_stat)  
    CALL reciprocal_k(aux%tr_r, aux_k%tr, aux_k%tr_r, &
      acc = .TRUE._jplm, map = map, prof_stat = prof_stat)
    
    IF (coef%nozone > 0) THEN
      CALL sqrt_k(aux%or_sqrt, aux_k%or, aux_k%or_sqrt, &
        acc = .TRUE._jplm, map = map, prof_stat = prof_stat)  
      DO i = 1, nchannels
        prof = chanprof(i)%prof
        aux_k%ow_rsqrt(:,i) = &!aux_k%ow_rsqrt + 
          aux%ow(:,prof) * aux_k%ow_sqrt(:,i)

        aux_k%ow(:,i) = aux_k%ow(:,i) + aux_k%ow_sqrt(:,i) * aux%ow_rsqrt(:,prof)
        aux_k%ow_rsqrt(:,i) = aux_k%ow_rsqrt(:,i) + &
          2._jprb * aux%ow_rsqrt(:,prof) * aux_k%ow_r(:,i)
      ENDDO
      CALL INVSQRT_K(aux%ow_rsqrt, aux_k%ow, aux_k%ow_rsqrt, &
        acc = .TRUE._jplm, map = map, prof_stat = prof_stat)
    ENDIF
    
! 4. calculate profile / reference profile sums: tw ww ow

    DO i = 1, nchannels
      DO lay = nlayers, 2, -1
        ! cumulate overlying layers: weighting tr relates to same layer as dpp
        ! do not need dpp(0) to start
        aux_k%tw(lay - 1, i) = aux_k%tw(lay - 1, i) + aux_k%tw(lay, i)
        aux_k%tr(lay - 1, i) = aux_k%tr(lay - 1, i) + coef%dpp(lay - 1) * aux_k%tw(lay, i) 
      ENDDO
    ENDDO

    DO i = 1, nchannels
      sum1 = 0._JPRB
      ! cumulating column overlying layer and layer itself
      DO lay = nlayers, 1, -1
        ! cumulate overlying layers: weighting w or wstar relates to layer below dpp
        ! need dpp(0) to start
        sum1 = sum1 + aux_k%ww(lay, i) * aux%SUM(lay,1)
        aux_k%w_layer(lay, i) = &!aux_k%w_layer(lay, i) + 
          sum1 * coef%dpp(lay - 1)
      ENDDO
    ENDDO
 
    IF (coef%nozone > 0) THEN
      IF (opts%rt_ir%ozone_Data) THEN     ! if no input O3 profile, set to reference value (ow =1)
        DO i = 1, nchannels
          sum1 = 0._JPRB
          DO lay = nlayers, 1, -1
            ! cumulate overlying layers: weighting o or ostar relates to layer below dpp
            ! need dpp(0) to start
            sum1 = sum1 + aux_k%ow(lay, i) * aux%SUM(lay,2)
            aux_k%o3_layer(lay, i) = &!aux_k%o3_layer(lay, i) + 
              sum1 * coef%dpp(lay - 1) ! OK, dpp(0) defined 
          ENDDO
        ENDDO
      ELSE
        aux_k%o3_layer(:, :) = 0._JPRB
      ENDIF
    ENDIF

! 3. calculate (profile / reference profile) ratios; tr wr or
    DO i = 1, nchannels
      prof = chanprof(i)%prof
      IF ((coef%id_inst == inst_id_ssmis .OR. coef%id_inst == inst_id_ssmisz) .AND. coef%IncZeeman) THEN
        aux_k%t_layer(:, i) = aux_k%t_layer(:, i) + aux_k%tr(:, i) * tstar_r(:)
      ELSE
        aux_k%t_layer(:, i) = aux_k%tr(:, i) * tstar_r(:)
      ENDIF

      aux_k%w_layer(:, i) = aux_k%w_layer(:, i) + aux_k%wr(:, i) * wstar_r(:)
      
! if no input O3 profile, set to reference value (or = 1)

      IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN
        aux_k%o3_layer(:, i) = aux_k%o3_layer(:, i) + aux_k%or(:, i)  * ostar_r(:)
      ENDIF

! 2. calculate, for layers, deviations from reference profile
      aux_k%t_layer(:, i) = aux_k%t_layer(:, i) + aux_k%dt(:, i) 

      ! if no input O3 profile we still use the input temperature profile for dto
      IF (coef%nozone > 0) &
        aux_k%t_layer(:, i) = aux_k%t_layer(:, i) + aux_k%dto(:, i) 

! 1 profile layer quantities
!   the layer number agrees with the level number of its upper boundary
! layer N-1 lies between levels N-1 and N
    !DAR add from set_predictors_7
      IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN
        profiles_k(i)%o3(1) = profiles_k(i)%o3(1) + &
          0.5_JPRB * aux_k%o3_layer(1, i)
        profiles_k(i)%o3(2:nlevels-1) = profiles_k(i)%o3(2:nlevels-1) + &
          0.5_JPRB * (aux_k%o3_layer(1:nlevels-2, i) + aux_k%o3_layer(2:nlevels-1, i))
        profiles_k(i)%o3(nlevels) = profiles_k(i)%o3(nlevels) + &
          0.5_JPRB * aux_k%o3_layer(nlevels-1, i)
      ENDIF

      profiles_k(i)%t(1) = profiles_k(i)%t(1) + &
        0.5_JPRB * aux_k%t_layer(1, i)
      profiles_k(i)%t(2:nlevels-1) = profiles_k(i)%t(2:nlevels-1) + &
        0.5_jprb * (aux_k%t_layer(1:nlevels-2, i) + aux_k%t_layer(2:nlevels-1, i))
      profiles_k(i)%t(nlevels) = profiles_k(i)%t(nlevels) + &
        0.5_JPRB * aux_k%t_layer(nlevels-1, i)

      profiles_k(i)%q(1) = profiles_k(i)%q(1) + &
        0.5_JPRB * aux_k%w_layer(1, i)
      profiles_k(i)%q(2:nlevels-1) = profiles_k(i)%q(2:nlevels-1) + &
        0.5_JPRB * (aux_k%w_layer(1:nlevels-2, i) + aux_k%w_layer(2:nlevels-1, i))
      profiles_k(i)%q(nlevels) = profiles_k(i)%q(nlevels) + &
        0.5_JPRB * aux_k%w_layer(nlevels-1, i)

      IF (opts%rt_all%use_q2m) THEN
        iv3lev = aux%s(prof)%nearestlev_surf - 1! nearest level above surface
        iv2lev = aux%s(prof)%nearestlev_surf    ! nearest level above surface

        IF (iv2lev <= coef%nlevels) THEN
          iv2lay       = iv2lev - 1
          profiles_k(i)%s2m%q =  profiles_k(i)%s2m%q + &
            aux_k%w_layer(iv2lay, i) * 0.5_JPRB
          
          profiles_k(i)%q(iv2lev) = profiles_k(i)%q(iv2lev) - &
            aux_k%w_layer(iv2lay, i) * 0.5_JPRB
        ENDIF
      ENDIF

    ENDDO
    
  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_PROFAUX_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_profaux_k
