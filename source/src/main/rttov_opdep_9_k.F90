SUBROUTINE rttov_opdep_9_k( &
            & nlayers,      &
            & chanprof,     &
            & chanflag,     &
            & pred_info,    &
            & predictors,   &
            & predictors_k, &
            & aux,          &
            & aux_k,        &
            & coef,         &
            & fast_coef,    &
            & opdp_path,    &
            & opdp_path_k,  &
            & opdp_ref)
! Description:
! K of rttov_opdep_9
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
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0    01/06/2005  Marco Matricardi (ECMWF):
!            --       New routine based on rttov_transmit_k.F90.
!                     Variable trace gases, CO2, N2O,CO and CH4
!                     have been introduced for AIRS and IASI
! 1.1      06/02/2007 Removed polaristion index R Saunders
! 1.2      07/04/2008 Marco Matricardi (ECMWF):
!                     modified the routine to be consistent
!                     with the LBLRTM coefficients
! 1.3      23/04/2008 Removed remaining use of fast_coef%watervapourint array P. Brunel
! 1.4      18/06/2009 Use rttov9 constant intervals from rttov_const
!                     and change test on upper bounds from <= to <
!                     Philippe Marguinaud
! 1.4.1    15/09/2009 User defined ToA. Layers distinct from levels (P.Rayer)
! 1.5      03/11/2009 Transmittances / optical depths on levels (A Geer)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
  USE rttov_types, ONLY : &
       & rttov_chanprof,  &
       & rttov_coef,      &
       & rttov_fast_coef, &
       & predictors_Type, &
       & rttov_path_pred, &
       & profile_aux
  USE parkind1, ONLY : jpim, jprb, jplm
!INTF_OFF
  USE rttov_const, ONLY :  &
       & sensor_id_mw,     &
       & rttov9_wv0690_50, &
       & rttov9_wv1050_00, &
       & rttov9_wv1095_25, &
       & rttov9_wv1100_25, &
       & rttov9_wv1350_25, &
       & rttov9_wv1750_25, &
       & rttov9_wv1900_25, &
       & rttov9_wv1995_00, &
       & rttov9_wv2000_00, &
       & rttov9_wv2250_00, &
       & rttov9_wv2295_25, &
       & rttov9_wv2380_25
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  INTEGER(KIND=jpim)   , INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof ), INTENT(IN)    :: chanprof(:)
  LOGICAL(KIND=jplm)   , INTENT(IN)    :: chanflag(SIZE(chanprof))
  TYPE(predictors_Type), INTENT(IN)    :: pred_info
  TYPE(rttov_path_pred), INTENT(IN)    :: predictors
  TYPE(rttov_path_pred), INTENT(INOUT) :: predictors_k
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  TYPE(rttov_fast_coef), INTENT(IN)    :: fast_coef
  TYPE(profile_aux    ), INTENT(IN)    :: aux
  TYPE(profile_aux    ), INTENT(INOUT) :: aux_k
  REAL(KIND=jprb)      , INTENT(IN)    :: opdp_path(:,:)
  REAL(KIND=jprb)      , INTENT(INOUT) :: opdp_path_k(:,:)
  REAL(KIND=jprb)      , INTENT(IN)    :: opdp_ref(nlayers, SIZE(chanprof))
!INTF_END
!local variables:
  REAL(KIND=jprb) :: opticaldepth_k(nlayers, SIZE(chanprof))
  REAL(KIND=jprb) :: chanx
  REAL(KIND=jprb), POINTER :: debye_prof  (:, :)
  REAL(KIND=jprb), POINTER :: debye_prof_k(:, :)
  INTEGER(KIND=jpim) :: lev, lay, chan, j, ii, nlevels
  INTEGER(KIND=jpim) :: prof
! cloud liquid water local variables
  REAL   (KIND=jprb) :: zf, zf_sq, z34_dif, z45_dif, z1_sq, z2_sq, z1_div, z2_div
  REAL   (KIND=jprb) :: z1_den, z2_den, zastar, z1_prod, z2_prod, z3_prod, z4_prod
  REAL   (KIND=jprb) :: zbstar, zbstar_sq, za2star, za2star_sq, zdiv, zgstar
  REAL   (KIND=jprb) :: z1f_sq_z1_sq, z2f_sq_z2_sq
  REAL   (KIND=jprb) :: z34_dif_k, z45_dif_k, z1_sq_k, z2_sq_k, z1_div_k , z2_div_k
  REAL   (KIND=jprb) :: z1_den_k, z2_den_k, zastar_k, z1_prod_k, z2_prod_k, z3_prod_k, z4_prod_k
  REAL   (KIND=jprb) :: zbstar_k, zbstar_sq_k, za2star_k, za2star_sq_k, zdiv_k, zgstar_k
  REAL   (KIND=jprb) :: z1f_sq_z1_sq_k, z2f_sq_z2_sq_k
  INTEGER(KIND=jpim) :: nchannels
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!---------------------------------------------------------------------------------------
!----------------------------------------
!2.Assemble layer optical depths
!----------------------------------------
! note that optical depth in the calculations is negative
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_9_K', 0_jpim, ZHOOK_HANDLE)
  nchannels           = SIZE(chanprof)
  nlevels             = nlayers + 1
  opticaldepth_k(:,:) = 0._jprb
  DO j = 1, nchannels
    IF (chanflag(j)) THEN
      DO lev = nlevels, 2,  - 1
        lay = lev - 1
        opdp_path_k(lev - 1, j) = opdp_path_k(lev - 1, j) + opdp_path_k(lev, j)
        opticaldepth_k(lay, j)  = opticaldepth_k(lay, j) + opdp_path_k(lev, j)
        opdp_path_k(lev, j)     = 0.0_jprb
      ENDDO
      opdp_path_k(1, j) = 0.0_jprb
    ENDIF
  ENDDO
  DO j = 1, nchannels
    IF (chanflag(j)) THEN
      WHERE (opdp_ref(:,j) > 0.0_jprb) opticaldepth_k(:,j) = 0.0_jprb
    ENDIF
  ENDDO
!--------------------
!1.9 add liquid water (MW only)
!--------------------
  IF (coef%id_sensor == sensor_id_mw) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      prof = chanprof(j)%prof
      debye_prof => aux%debye_prof(:, :, prof)
      debye_prof_k => aux_k%debye_prof(:, :, j)
      IF (pred_info%ncloud >= 1) THEN
        DO lay = nlayers, 1,  - 1
          lev = lay + 1
          IF (lev >= coef%mwcldtop) THEN
! Repeat direct code
            zf = coef%frequency_ghz(chan)
            zf_sq = zf * zf
            z1_sq = debye_prof(1, lay) * debye_prof(1, lay)
            z2_sq = debye_prof(2, lay) * debye_prof(2, lay)
            z34_dif                  = debye_prof(3, lay) - debye_prof(4, lay)
            z45_dif                  = debye_prof(4, lay) - debye_prof(5, lay)
            z1f_sq_z1_sq             = zf_sq + z1_sq
            z2f_sq_z2_sq             = zf_sq + z2_sq
            z1_div                   = 1.0_jprb / z1f_sq_z1_sq
            z2_div                   = 1.0_jprb / z2f_sq_z2_sq
            z1_den                   = z34_dif * z1_div
            z2_den                   = z45_dif * z2_div
            zastar                   = debye_prof(3, lay) - zf_sq * (z1_den + z2_den)
            z1_prod                  = z34_dif * debye_prof(1, lay)
            z2_prod                  = z1_prod * z1_div
            z3_prod                  = z45_dif * debye_prof(2, lay)
            z4_prod                  = z3_prod * z2_div
            zbstar                   =  - zf * (z2_prod + z4_prod)
            zbstar_sq                = zbstar * zbstar
            za2star                  = zastar + 2.0_jprb
            za2star_sq               = za2star * za2star
            zdiv = za2star_sq + zbstar_sq
            zgstar                   =  - 3.0_jprb * zbstar / zdiv
! Now compute Adjoint code
!opticaldepth_k(lay,j)= opticaldepth_k(lay,j)
            zgstar_k                 = opticaldepth_k(lay, j) * ( - 1.5_jprb * zf * predictors%clw(lay, prof))
            predictors_k%clw(lay, j) = predictors_k%clw(lay, j) + opticaldepth_k(lay, j) * ( - 1.5_jprb * zf * zgstar)
            zbstar_k                 =  - 3.0_jprb * zgstar_k / zdiv
            zdiv_k                   = 3.0_jprb * zgstar_k * zbstar / (zdiv * zdiv)
!zgstar_k =  0.
            za2star_sq_k             = zdiv_k
            zbstar_sq_k              = zdiv_k
!zdiv_k       = 0.
            za2star_k                = 2.0_jprb * za2star * za2star_sq_k
!za2star_sq_k = 0.
            zastar_k                 = za2star_k
!za2star_k    = 0.
            zbstar_k                 = zbstar_k + 2.0_jprb * zbstar * zbstar_sq_k
!zbstar_sq_k  = 0.
            z2_prod_k                =  - zf * zbstar_k
            z4_prod_k                =  - zf * zbstar_k
!zbstar_k     = 0.
            z3_prod_k                = z2_div * z4_prod_k
            z2_div_k                 = z3_prod * z4_prod_k
!z4_prod_k    = 0.
            z45_dif_k                = debye_prof(2, lay) * z3_prod_k
            debye_prof_k(2, lay)     = debye_prof_k(2, lay) + z45_dif * z3_prod_k
!z3_prod_k           = 0.
            z1_prod_k                = z1_div * z2_prod_k
            z1_div_k                 = z1_prod * z2_prod_k
!z2_prod_k    = 0.
            z34_dif_k                = debye_prof(1, lay) * z1_prod_k
            debye_prof_k(1, lay)     = debye_prof_k(1, lay) + z34_dif * z1_prod_k
!z1_prod_k           = 0.
            debye_prof_k(3, lay)     = debye_prof_k(3, lay) + zastar_k
            z1_den_k                 =  - zf_sq * zastar_k
            z2_den_k                 =  - zf_sq * zastar_k
!zastar_k            = 0.
            z2_div_k                 = z2_div_k + z45_dif * z2_den_k
            z45_dif_k                = z45_dif_k + z2_div * z2_den_k
!z2_den_k     = 0.
            z1_div_k                 = z1_div_k + z34_dif * z1_den_k
            z34_dif_k                = z34_dif_k + z1_div * z1_den_k
!z1_den_k     = 0.
            z2f_sq_z2_sq_k           =  - z2_div_k / (z2f_sq_z2_sq * z2f_sq_z2_sq)
!z2_div_k       = 0.
            z1f_sq_z1_sq_k           =  - z1_div_k / (z1f_sq_z1_sq * z1f_sq_z1_sq)
!z1_div_k       = 0.
            z2_sq_k                  = z2f_sq_z2_sq_k
!z2f_sq_z2_sq_k = 0.
            z1_sq_k                  = z1f_sq_z1_sq_k
!z1f_sq_z1_sq_k = 0.
            debye_prof_k(4, lay)     = debye_prof_k(4, lay) + z45_dif_k
            debye_prof_k(5, lay)     = debye_prof_k(5, lay) - z45_dif_k
!z45_dif_k           = 0.
            debye_prof_k(3, lay)     = debye_prof_k(3, lay) + z34_dif_k
            debye_prof_k(4, lay)     = debye_prof_k(4, lay) - z34_dif_k
!z34_dif_k           = 0.
            debye_prof_k(2, lay)     = debye_prof_k(2, lay) + z2_sq_k * 2.0_jprb * debye_prof(2, lay)
!z2_sq_k             = 0.
            debye_prof_k(1, lay)     = debye_prof_k(1, lay) + z1_sq_k * 2.0_jprb * debye_prof(1, lay)
!z1_sq_k             = 0.
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDIF
!-----------
!1.8 add CH4
!-----------
  IF (coef%nch4 > 0) THEN
    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        chan = chanprof(j)%chan
        DO lay = nlayers, 1,  - 1
          predictors_k%ch4(:, lay, j) = predictors_k%ch4(:, lay, j) + fast_coef%ch4(lay, chan, :) * opticaldepth_k(lay, j)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
!-----------
!1.7 add CO
!-----------
  IF (coef%nco > 0) THEN
    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        chan  = chanprof(j)%chan
        chanx = coef%ff_cwn(chan)
        prof  = chanprof(j)%prof
        IF ((chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2380_25) .OR. &
            (chanx < rttov9_wv2000_00)) THEN
          DO lay = nlayers, 1,  - 1
            predictors_k%co(1:11, lay, j) =      &
              & predictors_k%co(1:11, lay, j) + fast_coef%co(lay, chan, 1:11) * opticaldepth_k(lay, j)
          ENDDO
        ELSE
          DO lay = nlayers, 1,  - 1
            predictors_k%co(1:10, lay, j) =      &
              & predictors_k%co(1:10, lay, j) + fast_coef%co(lay, chan, 1:10) * opticaldepth_k(lay, j)
            predictors_k%co(12, lay, j) = predictors_k%co(12, lay, j) + fast_coef%co(lay, chan, 11) * opticaldepth_k(lay, j)
            predictors_k%co(13, lay, j) = predictors_k%co(13, lay, j) + fast_coef%co(lay, chan, 12) * opticaldepth_k(lay, j)
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF
!-----------
!1.6 add N2O
!-----------
  IF (coef%nn2o > 0) THEN
    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        chan  = chanprof(j)%chan
        chanx = coef%ff_cwn(chan)
        prof  = chanprof(j)%prof
        IF (chanx >= rttov9_wv1050_00 .AND. chanx < rttov9_wv1350_25) THEN
          DO lay = nlayers, 1,  - 1
            IF (coef%nch4 > 0) THEN
              predictors_k%ch4(10, lay, j) =      &
                & predictors_k%ch4(10, lay, j) + fast_coef%n2o(lay, chan, 12) * opticaldepth_k(lay, j)
              predictors_k%ch4(1, lay, j)  =      &
                & predictors_k%ch4(1, lay, j) + fast_coef%n2o(lay, chan, 11) * opticaldepth_k(lay, j)
            ENDIF
            predictors_k%n2o(1:10, lay, j) =      &
              & predictors_k%n2o(1:10, lay, j) + fast_coef%n2o(lay, chan, 1:10) * opticaldepth_k(lay, j)
          ENDDO
        ELSE IF (chanx >= rttov9_wv1995_00 .AND. chanx < rttov9_wv2000_00) THEN
          DO lay = nlayers, 1,  - 1
            IF (coef%nco > 0) THEN
              predictors_k%co(8, lay, j) = predictors_k%co(8, lay, j) -                                  &
                & fast_coef%n2o(lay, chan, 12) * opticaldepth_k(lay, j) * predictors%co(1, lay, prof) ** 3 /  &
                & predictors%co(8, lay, prof) ** 2
              predictors_k%co(1, lay, j) = predictors_k%co(1, lay, j) +                                      &
                & fast_coef%n2o(lay, chan, 12) * opticaldepth_k(lay, j) * 3 * predictors%co(1, lay, prof) ** 2 /  &
                & predictors%co(8, lay, prof)
              predictors_k%co(1, lay, j) = predictors_k%co(1, lay, j) + fast_coef%n2o(lay, chan, 11) * opticaldepth_k(lay, j)
            ENDIF
            predictors_k%n2o(1:10, lay, j) =      &
              & predictors_k%n2o(1:10, lay, j) + fast_coef%n2o(lay, chan, 1:10) * opticaldepth_k(lay, j)
          ENDDO
        ELSE IF (chanx < rttov9_wv2000_00) THEN
          DO lay = nlayers, 1,  - 1
            predictors_k%n2o(1:10, lay, j) =      &
              & predictors_k%n2o(1:10, lay, j) + fast_coef%n2o(lay, chan, 1:10) * opticaldepth_k(lay, j)
          ENDDO
        ELSE IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2250_00) THEN
          DO lay = nlayers, 1,  - 1
            predictors_k%n2o(13, lay, j) = predictors_k%n2o(13, lay, j) + fast_coef%n2o(lay, chan, 14) * opticaldepth_k(lay, j)
            predictors_k%n2o(12, lay, j) = predictors_k%n2o(12, lay, j) + fast_coef%n2o(lay, chan, 13) * opticaldepth_k(lay, j)
            IF (coef%nco > 0) THEN
              predictors_k%co(8, lay, j) = predictors_k%co(8, lay, j) -                                  &
                & fast_coef%n2o(lay, chan, 12) * opticaldepth_k(lay, j) * predictors%co(1, lay, prof) ** 3 /  &
                & predictors%co(8, lay, prof) ** 2
              predictors_k%co(1, lay, j) = predictors_k%co(1, lay, j) +                                      &
                & fast_coef%n2o(lay, chan, 12) * opticaldepth_k(lay, j) * 3 * predictors%co(1, lay, prof) ** 2 /  &
                & predictors%co(8, lay, prof)
              predictors_k%co(1, lay, j) = predictors_k%co(1, lay, j) + fast_coef%n2o(lay, chan, 11) * opticaldepth_k(lay, j)
            ENDIF
            predictors_k%n2o(10, lay, j)  =      &
              & predictors_k%n2o(10, lay, j) + fast_coef%n2o(lay, chan, 10) * opticaldepth_k(lay, j)
            predictors_k%n2o(11, lay, j)  = predictors_k%n2o(11, lay, j) + fast_coef%n2o(lay, chan, 9) * opticaldepth_k(lay, j)
            predictors_k%n2o(1:8, lay, j) =      &
              & predictors_k%n2o(1:8, lay, j) + fast_coef%n2o(lay, chan, 1:8) * opticaldepth_k(lay, j)
          ENDDO
        ELSE IF (chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2295_25) THEN
          DO lay = nlayers, 1,  - 1
            IF (coef%nco > 0) THEN
              predictors_k%co(8, lay, j) = predictors_k%co(8, lay, j) -                                  &
                & fast_coef%n2o(lay, chan, 12) * opticaldepth_k(lay, j) * predictors%co(1, lay, prof) ** 3 /  &
                & predictors%co(8, lay, prof) ** 2
              predictors_k%co(1, lay, j) = predictors_k%co(1, lay, j) +                                      &
                & fast_coef%n2o(lay, chan, 12) * opticaldepth_k(lay, j) * 3 * predictors%co(1, lay, prof) ** 2 /  &
                & predictors%co(8, lay, prof)
              predictors_k%co(1, lay, j) = predictors_k%co(1, lay, j) + fast_coef%n2o(lay, chan, 11) * opticaldepth_k(lay, j)
            ENDIF
            predictors_k%n2o(1:10, lay, j) =      &
              & predictors_k%n2o(1:10, lay, j) + fast_coef%n2o(lay, chan, 1:10) * opticaldepth_k(lay, j)
          ENDDO
        ELSE IF (chanx >= rttov9_wv2295_25 .AND. chanx < rttov9_wv2380_25) THEN
          DO lay = nlayers, 1,  - 1
            predictors_k%n2o(1:10, lay, j) =      &
              & predictors_k%n2o(1:10, lay, j) + fast_coef%n2o(lay, chan, 1:10) * opticaldepth_k(lay, j)
          ENDDO
        ELSE   ! chanx >= rttov9_wv2380_25
          DO lay = nlayers, 1,  - 1
            predictors_k%n2o(13, lay, j)  =      &
              & predictors_k%n2o(13, lay, j) + fast_coef%n2o(lay, chan, 12) * opticaldepth_k(lay, j)
            predictors_k%n2o(12, lay, j)  =      &
              & predictors_k%n2o(12, lay, j) + fast_coef%n2o(lay, chan, 11) * opticaldepth_k(lay, j)
            predictors_k%n2o(10, lay, j)  =      &
              & predictors_k%n2o(10, lay, j) + fast_coef%n2o(lay, chan, 10) * opticaldepth_k(lay, j)
            predictors_k%n2o(11, lay, j)  = predictors_k%n2o(11, lay, j) + fast_coef%n2o(lay, chan, 9) * opticaldepth_k(lay, j)
            predictors_k%n2o(1:8, lay, j) =      &
              & predictors_k%n2o(1:8, lay, j) + fast_coef%n2o(lay, chan, 1:8) * opticaldepth_k(lay, j)
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF
!-----------
!1.5 add CO2
!-----------
  IF (coef%nco2 > 0) THEN
    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        chan  = chanprof(j)%chan
        chanx = coef%ff_cwn(chan)
        prof  = chanprof(j)%prof
        IF (chanx >= rttov9_wv0690_50 .AND. chanx < rttov9_wv1100_25) THEN
          DO lay = nlayers, 1,  - 1
            predictors_k%co2(9, lay, j) = predictors_k%co2(9, lay, j) + fast_coef%co2(lay, chan, 9) * opticaldepth_k(lay, j)
            predictors_k%co2(8, lay, j) = predictors_k%co2(8, lay, j) + fast_coef%co2(lay, chan, 8) * opticaldepth_k(lay, j)
            predictors_k%co2(7, lay, j) = predictors_k%co2(7, lay, j) +      &
              & fast_coef%co2(lay, chan, 7) * 2 * predictors%co2(7, lay, prof) * opticaldepth_k(lay, j)
            predictors_k%co2(6, lay, j) = predictors_k%co2(6, lay, j) +      &
              & fast_coef%co2(lay, chan, 6) * 2 * predictors%co2(6, lay, prof) * opticaldepth_k(lay, j)
            predictors_k%co2(5, lay, j) = predictors_k%co2(5, lay, j) + fast_coef%co2(lay, chan, 5) * opticaldepth_k(lay, j)
            predictors_k%co2(3, lay, j) = predictors_k%co2(3, lay, j) +      &
              & fast_coef%co2(lay, chan, 4) * 2 * predictors%co2(3, lay, prof) * opticaldepth_k(lay, j)
            predictors_k%co2(5, lay, j) = predictors_k%co2(5, lay, j) +      &
              & fast_coef%co2(lay, chan, 3) * predictors%co2(6, lay, prof) ** 2 * opticaldepth_k(lay, j)
            predictors_k%co2(6, lay, j) = predictors_k%co2(6, lay, j) +                                     &
              & fast_coef%co2(lay, chan, 3) * 2 * predictors%co2(6, lay, prof) * predictors%co2(5, lay, prof) *  &
              & opticaldepth_k(lay, j)
            predictors_k%co2(2, lay, j) = predictors_k%co2(2, lay, j) + fast_coef%co2(lay, chan, 2) * opticaldepth_k(lay, j)
            predictors_k%co2(1, lay, j) = predictors_k%co2(1, lay, j) + fast_coef%co2(lay, chan, 1) * opticaldepth_k(lay, j)
          ENDDO
        ELSE IF (chanx >= rttov9_wv1995_00 .AND. chanx < rttov9_wv2000_00) THEN
          DO lay = nlayers, 1,  - 1
            predictors_k%co2(10, lay, j) = predictors_k%co2(10, lay, j) + fast_coef%co2(lay, chan, 11) * opticaldepth_k(lay, j)
            IF (coef%nco > 0) THEN
              predictors_k%co(1, lay, j) = predictors_k%co(1, lay, j) + fast_coef%co2(lay, chan, 10) * opticaldepth_k(lay, j)
            ENDIF
            predictors_k%co2(1:9, lay, j) =      &
              & predictors_k%co2(1:9, lay, j) + fast_coef%co2(lay, chan, 1:9) * opticaldepth_k(lay, j)
          ENDDO
        ELSE IF (chanx < rttov9_wv2000_00) THEN
          DO lay = nlayers, 1,  - 1
            predictors_k%co2(1:9, lay, j) =      &
              & predictors_k%co2(1:9, lay, j) + fast_coef%co2(lay, chan, 1:9) * opticaldepth_k(lay, j)
          ENDDO
        ELSE IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2250_00) THEN
          DO lay = nlayers, 1,  - 1
            predictors_k%co2(15, lay, j) = predictors_k%co2(15, lay, j) + fast_coef%co2(lay, chan, 14) * opticaldepth_k(lay, j)
            predictors_k%co2(14, lay, j) = predictors_k%co2(14, lay, j) + fast_coef%co2(lay, chan, 13) * opticaldepth_k(lay, j)
            predictors_k%co2(13, lay, j) = predictors_k%co2(13, lay, j) + fast_coef%co2(lay, chan, 12) * opticaldepth_k(lay, j)
            predictors_k%co2(12, lay, j) = predictors_k%co2(12, lay, j) + fast_coef%co2(lay, chan, 11) * opticaldepth_k(lay, j)
            predictors_k%co2(11, lay, j) = predictors_k%co2(11, lay, j) + fast_coef%co2(lay, chan, 10) * opticaldepth_k(lay, j)
            predictors_k%co2(10, lay, j) = predictors_k%co2(10, lay, j) + fast_coef%co2(lay, chan, 9) * opticaldepth_k(lay, j)
            IF (coef%nco > 0) THEN
              predictors_k%co(1, lay, j) = predictors_k%co(1, lay, j) + fast_coef%co2(lay, chan, 8) * opticaldepth_k(lay, j)
            ENDIF
            predictors_k%co2(8, lay, j)   = predictors_k%co2(8, lay, j) + fast_coef%co2(lay, chan, 7) * opticaldepth_k(lay, j)
            predictors_k%co2(7, lay, j)   = predictors_k%co2(7, lay, j) + fast_coef%co2(lay, chan, 6) * opticaldepth_k(lay, j)
            predictors_k%co2(1:5, lay, j) =      &
              & predictors_k%co2(1:5, lay, j) + fast_coef%co2(lay, chan, 1:5) * opticaldepth_k(lay, j)
          ENDDO
        ELSE IF (chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2295_25) THEN
          DO lay = nlayers, 1,  - 1
            predictors_k%co2(10, lay, j) = predictors_k%co2(10, lay, j) + fast_coef%co2(lay, chan, 11) * opticaldepth_k(lay, j)
            IF (coef%nco > 0) THEN
              predictors_k%co(1, lay, j) = predictors_k%co(1, lay, j) + fast_coef%co2(lay, chan, 10) * opticaldepth_k(lay, j)
            ENDIF
            predictors_k%co2(1:9, lay, j) =      &
              & predictors_k%co2(1:9, lay, j) + fast_coef%co2(lay, chan, 1:9) * opticaldepth_k(lay, j)
          ENDDO
        ELSE IF (chanx >= rttov9_wv2295_25 .AND. chanx < rttov9_wv2380_25) THEN
          DO lay = nlayers, 1,  - 1
            predictors_k%co2(1:10, lay, j) =      &
              & predictors_k%co2(1:10, lay, j) + fast_coef%co2(lay, chan, 1:10) * opticaldepth_k(lay, j)
          ENDDO
        ELSE   ! chanx >= rttov9_wv2380_25
          DO lay = nlayers, 1,  - 1
            predictors_k%co2(15, lay, j)  =      &
              & predictors_k%co2(15, lay, j) + fast_coef%co2(lay, chan, 13) * opticaldepth_k(lay, j)
            predictors_k%co2(14, lay, j)  =      &
              & predictors_k%co2(14, lay, j) + fast_coef%co2(lay, chan, 12) * opticaldepth_k(lay, j)
            predictors_k%co2(13, lay, j)  =      &
              & predictors_k%co2(13, lay, j) + fast_coef%co2(lay, chan, 11) * opticaldepth_k(lay, j)
            predictors_k%co2(12, lay, j)  =      &
              & predictors_k%co2(12, lay, j) + fast_coef%co2(lay, chan, 10) * opticaldepth_k(lay, j)
            predictors_k%co2(11, lay, j)  = predictors_k%co2(11, lay, j) + fast_coef%co2(lay, chan, 9) * opticaldepth_k(lay, j)
            predictors_k%co2(10, lay, j)  = predictors_k%co2(10, lay, j) + fast_coef%co2(lay, chan, 8) * opticaldepth_k(lay, j)
            predictors_k%co2(8, lay, j)   = predictors_k%co2(8, lay, j) + fast_coef%co2(lay, chan, 7) * opticaldepth_k(lay, j)
            predictors_k%co2(7, lay, j)   = predictors_k%co2(7, lay, j) + fast_coef%co2(lay, chan, 6) * opticaldepth_k(lay, j)
            predictors_k%co2(1:5, lay, j) =      &
              & predictors_k%co2(1:5, lay, j) + fast_coef%co2(lay, chan, 1:5) * opticaldepth_k(lay, j)
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF
!------------------------------
!1.4 add Water Vapour Continuum
!------------------------------
  IF (coef%nwvcont > 0) THEN
    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        chan = chanprof(j)%chan
        DO lay = nlayers, 1,  - 1
          predictors_k%wvcont(:, lay, j) =      &
            & predictors_k%wvcont(:, lay, j) + fast_coef%wvcont(lay, chan, :) * opticaldepth_k(lay, j)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
!-------------
!1.3 add ozone
!-------------
  IF (coef%nozone > 0) THEN
    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        chan  = chanprof(j)%chan
        chanx = coef%ff_cwn(chan)
        prof  = chanprof(j)%prof
        IF ((chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2380_25) .OR. &
            (chanx < rttov9_wv2000_00)) THEN
          DO lay = nlayers, 1,  - 1
            predictors_k%ozone(1:11, lay, j) =      &
              & predictors_k%ozone(1:11, lay, j) + fast_coef%ozone(lay, chan, 1:11) * opticaldepth_k(lay, j)
          ENDDO
        ELSE
          DO lay = nlayers, 1,  - 1
            predictors_k%ozone(15, lay, j) =      &
              & predictors_k%ozone(15, lay, j) + fast_coef%ozone(lay, chan, 13) * opticaldepth_k(lay, j)
            predictors_k%ozone(14, lay, j) =      &
              & predictors_k%ozone(14, lay, j) + fast_coef%ozone(lay, chan, 12) * opticaldepth_k(lay, j)
            predictors_k%ozone(11, lay, j) =      &
              & predictors_k%ozone(11, lay, j) + fast_coef%ozone(lay, chan, 11) * opticaldepth_k(lay, j)
            predictors_k%ozone(10, lay, j) =      &
              & predictors_k%ozone(10, lay, j) + fast_coef%ozone(lay, chan, 10) * opticaldepth_k(lay, j)
            predictors_k%ozone(9, lay, j)  =      &
              & predictors_k%ozone(9, lay, j) + fast_coef%ozone(lay, chan, 9) * opticaldepth_k(lay, j)
            predictors_k%ozone(13, lay, j) =      &
              & predictors_k%ozone(13, lay, j) + fast_coef%ozone(lay, chan, 8) * opticaldepth_k(lay, j)
            predictors_k%ozone(7, lay, j)  =      &
              & predictors_k%ozone(7, lay, j) + fast_coef%ozone(lay, chan, 7) * opticaldepth_k(lay, j)
            predictors_k%ozone(6, lay, j)  =      &
              & predictors_k%ozone(6, lay, j) + fast_coef%ozone(lay, chan, 6) * opticaldepth_k(lay, j)
            predictors_k%ozone(5, lay, j)  =      &
              & predictors_k%ozone(5, lay, j) + fast_coef%ozone(lay, chan, 5) * opticaldepth_k(lay, j)
            predictors_k%ozone(4, lay, j)  =      &
              & predictors_k%ozone(4, lay, j) + fast_coef%ozone(lay, chan, 4) * opticaldepth_k(lay, j)
            predictors_k%ozone(12, lay, j) =      &
              & predictors_k%ozone(12, lay, j) + fast_coef%ozone(lay, chan, 3) * opticaldepth_k(lay, j)
            predictors_k%ozone(2, lay, j)  =      &
              & predictors_k%ozone(2, lay, j) + fast_coef%ozone(lay, chan, 2) * opticaldepth_k(lay, j)
            predictors_k%ozone(1, lay, j)  =      &
              & predictors_k%ozone(1, lay, j) + fast_coef%ozone(lay, chan, 1) * opticaldepth_k(lay, j)
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF
!--------------------
!1.2 add water vapour
!--------------------
  DO j = 1, nchannels
    IF (chanflag(j)) THEN
      chan  = chanprof(j)%chan
      chanx = coef%ff_cwn(chan)
      prof  = chanprof(j)%prof
      IF (chanx >= rttov9_wv1095_25 .AND. chanx < rttov9_wv1750_25) THEN
        DO lay = nlayers, 1,  - 1
          IF (coef%nch4 > 0) THEN
            predictors_k%ch4(1, lay, j) = predictors_k%ch4(1, lay, j) +      &
              & opticaldepth_k(lay, j) * fast_coef%watervapour(lay, chan, 15) * predictors%ch4(3, lay, prof)
            predictors_k%ch4(3, lay, j) = predictors_k%ch4(3, lay, j) +      &
              & opticaldepth_k(lay, j) * fast_coef%watervapour(lay, chan, 15) * predictors%ch4(1, lay, prof)
            predictors_k%ch4(1, lay, j) =      &
              & predictors_k%ch4(1, lay, j) + opticaldepth_k(lay, j) * fast_coef%watervapour(lay, chan, 14)
          ENDIF
          predictors_k%watervapour(1:13, lay, j) =      &
            & predictors_k%watervapour(1:13, lay, j) + fast_coef%watervapour(lay, chan, 1:13) * opticaldepth_k(lay, j)
        ENDDO
      ELSE IF (chanx >= rttov9_wv1900_25 .AND. chanx < rttov9_wv1995_00) THEN
        DO lay = nlayers, 1,  - 1
          IF (coef%nco2 > 0) THEN
            predictors_k%co2(1, lay, j) =      &
              & predictors_k%co2(1, lay, j) + opticaldepth_k(lay, j) * fast_coef%watervapour(lay, chan, 14)
          ENDIF
          predictors_k%watervapour(1:13, lay, j) =      &
            & predictors_k%watervapour(1:13, lay, j) + fast_coef%watervapour(lay, chan, 1:13) * opticaldepth_k(lay, j)
        ENDDO
      ELSE IF (chanx >= rttov9_wv1995_00 .AND. chanx < rttov9_wv2000_00) THEN
        DO lay = nlayers, 1,  - 1
          IF (coef%nco > 0) THEN
            predictors_k%co(1, lay, j) =      &
              & predictors_k%co(1, lay, j) + fast_coef%watervapour(lay, chan, 15) * opticaldepth_k(lay, j)
          ENDIF
          IF (coef%nco2 > 0) THEN
            predictors_k%co2(1, lay, j) =      &
              & predictors_k%co2(1, lay, j) + fast_coef%watervapour(lay, chan, 14) * opticaldepth_k(lay, j)
          ENDIF
          DO ii = 1, 13
            predictors_k%watervapour(ii, lay, j) =      &
              & predictors_k%watervapour(ii, lay, j) + fast_coef%watervapour(lay, chan, ii) * opticaldepth_k(lay, j)
          ENDDO
        ENDDO
      ELSE IF (chanx < rttov9_wv2000_00) THEN
        DO lay = nlayers, 1,  - 1
          DO ii = 1, 13
            predictors_k%watervapour(ii, lay, j) =      &
              & predictors_k%watervapour(ii, lay, j) + fast_coef%watervapour(lay, chan, ii) * opticaldepth_k(lay, j)
          ENDDO
        ENDDO
      ELSE IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2250_00) THEN
        IF (coef%ss_val_chn(chan) > 0) THEN
          DO lay = nlayers, 1,  - 1
            IF (coef%nco > 0) THEN
              predictors_k%co(1, lay, j) =      &
                & predictors_k%co(1, lay, j) + fast_coef%watervapour(lay, chan, 15) * opticaldepth_k(lay, j)
            ENDIF
            IF (coef%nco2 > 0) THEN
              predictors_k%co2(1, lay, j) =      &
                & predictors_k%co2(1, lay, j) + fast_coef%watervapour(lay, chan, 14) * opticaldepth_k(lay, j)
            ENDIF
            predictors_k%watervapour(17, lay, j) =      &
              & predictors_k%watervapour(17, lay, j) + fast_coef%watervapour(lay, chan, 13) * opticaldepth_k(lay, j)
            predictors_k%watervapour(16, lay, j) =      &
              & predictors_k%watervapour(16, lay, j) + fast_coef%watervapour(lay, chan, 12) * opticaldepth_k(lay, j)
            predictors_k%watervapour(11, lay, j) =      &
              & predictors_k%watervapour(11, lay, j) + fast_coef%watervapour(lay, chan, 11) * opticaldepth_k(lay, j)
            predictors_k%watervapour(10, lay, j) =      &
              & predictors_k%watervapour(10, lay, j) + fast_coef%watervapour(lay, chan, 10) * opticaldepth_k(lay, j)
            predictors_k%watervapour(15, lay, j) =      &
              & predictors_k%watervapour(15, lay, j) + fast_coef%watervapour(lay, chan, 9) * opticaldepth_k(lay, j)
            predictors_k%watervapour(14, lay, j) =      &
              & predictors_k%watervapour(14, lay, j) + fast_coef%watervapour(lay, chan, 8) * opticaldepth_k(lay, j)
            DO ii = 1, 7
              predictors_k%watervapour(ii, lay, j) =      &
                & predictors_k%watervapour(ii, lay, j) + fast_coef%watervapour(lay, chan, ii) * opticaldepth_k(lay, j)
            ENDDO
          ENDDO
        ELSE
          DO lay = nlayers, 1,  - 1
            IF (coef%nco > 0) THEN
              predictors_k%co(1, lay, j) =      &
                & predictors_k%co(1, lay, j) + fast_coef%watervapour(lay, chan, 15) * opticaldepth_k(lay, j)
            ENDIF
            IF (coef%nco2 > 0) THEN
              predictors_k%co2(1, lay, j) =      &
                & predictors_k%co2(1, lay, j) + fast_coef%watervapour(lay, chan, 14) * opticaldepth_k(lay, j)
            ENDIF
            DO ii = 1, 13
              predictors_k%watervapour(ii, lay, j) =      &
                & predictors_k%watervapour(ii, lay, j) + fast_coef%watervapour(lay, chan, ii) * opticaldepth_k(lay, j)
            ENDDO
          ENDDO
        ENDIF
      ELSE IF (chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2295_25) THEN
        DO lay = nlayers, 1,  - 1
          IF (coef%nco > 0) THEN
            predictors_k%co(1, lay, j) =      &
              & predictors_k%co(1, lay, j) + fast_coef%watervapour(lay, chan, 15) * opticaldepth_k(lay, j)
          ENDIF
          IF (coef%nco2 > 0) THEN
            predictors_k%co2(1, lay, j) =      &
              & predictors_k%co2(1, lay, j) + fast_coef%watervapour(lay, chan, 14) * opticaldepth_k(lay, j)
          ENDIF
          DO ii = 1, 13
            predictors_k%watervapour(ii, lay, j) =      &
              & predictors_k%watervapour(ii, lay, j) + fast_coef%watervapour(lay, chan, ii) * opticaldepth_k(lay, j)
          ENDDO
        ENDDO
      ELSE IF (chanx >= rttov9_wv2295_25 .AND. chanx < rttov9_wv2380_25) THEN
        DO lay = nlayers, 1,  - 1
          DO ii = 1, 13
            predictors_k%watervapour(ii, lay, j) =      &
              & predictors_k%watervapour(ii, lay, j) + fast_coef%watervapour(lay, chan, ii) * opticaldepth_k(lay, j)
          ENDDO
        ENDDO
      ELSE   !chanx >= rttov9_wv2380_25
        DO lay = nlayers, 1,  - 1
          IF (coef%nch4 > 0) THEN
            predictors_k%ch4(3, lay, j) = predictors_k%ch4(3, lay, j) +      &
              & fast_coef%watervapour(lay, chan, 15) * predictors%ch4(1, lay, prof) ** 0.25 * opticaldepth_k(lay, j)
            predictors_k%ch4(1, lay, j) = predictors_k%ch4(1, lay, j) +                               &
              & fast_coef%watervapour(lay, chan, 15) * 0.25 * predictors%ch4(1, lay, prof) ** ( - 0.75) *  &
              & predictors%ch4(3, lay, prof) * opticaldepth_k(lay, j)
            predictors_k%ch4(1, lay, j) = predictors_k%ch4(1, lay, j) +      &
              & fast_coef%watervapour(lay, chan, 14) * 1.25 * predictors%ch4(1, lay, prof) ** 0.25 * opticaldepth_k(lay, j)
          ENDIF
          predictors_k%watervapour(19, lay, j) =      &
            & predictors_k%watervapour(19, lay, j) + fast_coef%watervapour(lay, chan, 13) * opticaldepth_k(lay, j)
          predictors_k%watervapour(16, lay, j) =      &
            & predictors_k%watervapour(16, lay, j) + fast_coef%watervapour(lay, chan, 12) * opticaldepth_k(lay, j)
          predictors_k%watervapour(11, lay, j) =      &
            & predictors_k%watervapour(11, lay, j) + fast_coef%watervapour(lay, chan, 11) * opticaldepth_k(lay, j)
          predictors_k%watervapour(18, lay, j) =      &
            & predictors_k%watervapour(18, lay, j) + fast_coef%watervapour(lay, chan, 10) * opticaldepth_k(lay, j)
          predictors_k%watervapour(15, lay, j) =      &
            & predictors_k%watervapour(15, lay, j) + fast_coef%watervapour(lay, chan, 9) * opticaldepth_k(lay, j)
          predictors_k%watervapour(14, lay, j) =      &
            & predictors_k%watervapour(14, lay, j) + fast_coef%watervapour(lay, chan, 8) * opticaldepth_k(lay, j)
          DO ii = 1, 7
            predictors_k%watervapour(ii, lay, j) =      &
              & predictors_k%watervapour(ii, lay, j) + fast_coef%watervapour(lay, chan, ii) * opticaldepth_k(lay, j)
          ENDDO
        ENDDO
      ENDIF
    ENDIF
  ENDDO
!--------------------------
!1.1 start with mixed gases
!--------------------------
  DO j = 1, nchannels
    IF (chanflag(j)) THEN
      chan  = chanprof(j)%chan
      chanx = coef%ff_cwn(chan)
      prof  = chanprof(j)%prof
      IF ((chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2380_25) .OR. &
          (chanx < rttov9_wv2000_00)) THEN
        DO lay = nlayers, 1,  - 1
          predictors_k%mixedgas(1:8, lay, j) =      &
            & predictors_k%mixedgas(1:8, lay, j) + fast_coef%mixedgas(lay, chan, 1:8) * opticaldepth_k(lay, j)
        ENDDO
      ELSE
        DO lay = nlayers, 1,  - 1
          predictors_k%mixedgas(10, lay, j) =      &
            & predictors_k%mixedgas(10, lay, j) + fast_coef%mixedgas(lay, chan, 10) * opticaldepth_k(lay, j)
          predictors_k%mixedgas(9, lay, j)  =      &
            & predictors_k%mixedgas(9, lay, j) + fast_coef%mixedgas(lay, chan, 9) * opticaldepth_k(lay, j)
          DO ii = 1, 8
            predictors_k%mixedgas(ii, lay, j) =      &
              & predictors_k%mixedgas(ii, lay, j) + fast_coef%mixedgas(lay, chan, ii) * opticaldepth_k(lay, j)
          ENDDO
        ENDDO
      ENDIF
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_9_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_opdep_9_k
