SUBROUTINE rttov_opdep_9_ad( &
            & nlayers,       &
            & chanprof,      &
            & chanflag,      &
            & pred_info,     &
            & predictors,    &
            & predictors_ad, &
            & aux,           &
            & aux_ad,        &
            & coef,          &
            & fast_coef,     &
            & opdp_path,     &
            & opdp_path_ad,  &
            & opdp_ref)
! Description:
! AD of rttov_opdep_9
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
!  1.1    06/02/2007  Removed polarisation R Saunders
!  1.2    07/04/2008  Marco Matricardi (ECMWF):
!                     modified the routine to be consistent
!                     with the LBLRTM coefficients
!  1.3    15/08/2009  User defined ToA. Layers distinct from levels (P.Rayer)
!  1.4    18/06/2009  Use rttov9 constant intervals from rttov_const
!                     and change test on upper bounds from <= to <
!                     Philippe Marguinaud
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
  USE rttov_types, ONLY : &
       & rttov_chanprof,  &
       & rttov_coef,      &
       & rttov_fast_coef, &
       & predictors_Type, &
       & rttov_path_pred, &
       & profile_aux
  USE parkind1, ONLY : jprb, jpim, jplm
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
  TYPE(rttov_path_pred), INTENT(INOUT) :: predictors_ad
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  TYPE(rttov_fast_coef), INTENT(IN)    :: fast_coef
  TYPE(profile_aux    ), INTENT(IN)    :: aux
  TYPE(profile_aux    ), INTENT(INOUT) :: aux_ad
  REAL(KIND=jprb)      , INTENT(IN)    :: opdp_path(:,:)
  REAL(KIND=jprb)      , INTENT(INOUT) :: opdp_path_ad(:,:)
  REAL(KIND=jprb)      , INTENT(IN)    :: opdp_ref(nlayers, SIZE(chanprof))
!INTF_END
!local variables:
  REAL(KIND=jprb) :: opticaldepth_ad(nlayers, SIZE(chanprof))
  REAL(KIND=jprb) :: chanx
  REAL(KIND=jprb), POINTER :: debye_prof   (:, :)
  REAL(KIND=jprb), POINTER :: debye_prof_ad(:, :)
  INTEGER(KIND=jpim) :: lev, lay, chan, j, nlevels
  INTEGER(KIND=jpim) :: prof
! cloud liquid water local variables
  REAL   (KIND=jprb) :: zf, zf_sq, z34_dif, z45_dif, z1_sq, z2_sq, z1_div, z2_div
  REAL   (KIND=jprb) :: z1_den, z2_den, zastar, z1_prod, z2_prod, z3_prod, z4_prod
  REAL   (KIND=jprb) :: zbstar, zbstar_sq, za2star, za2star_sq, zdiv, zgstar
  REAL   (KIND=jprb) :: z1f_sq_z1_sq, z2f_sq_z2_sq
  REAL   (KIND=jprb) :: z34_dif_ad, z45_dif_ad, z1_sq_ad, z2_sq_ad, z1_div_ad, z2_div_ad
  REAL   (KIND=jprb) :: z1_den_ad, z2_den_ad, zastar_ad, z1_prod_ad, z2_prod_ad, z3_prod_ad, z4_prod_ad
  REAL   (KIND=jprb) :: zbstar_ad, zbstar_sq_ad, za2star_ad, za2star_sq_ad, zdiv_ad, zgstar_ad
  REAL   (KIND=jprb) :: z1f_sq_z1_sq_ad, z2f_sq_z2_sq_ad
  INTEGER(KIND=jpim) :: nchannels
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
!----------------------------------------
!2.Assemble layer optical depths
!----------------------------------------
! note that optical depth in the calculations is negative
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_9_AD', 0_jpim, ZHOOK_HANDLE)
  nchannels            = SIZE(chanprof)
  nlevels              = nlayers + 1
  opticaldepth_ad(:,:) = 0._jprb
  DO j = 1, nchannels
    IF (chanflag(j)) THEN
      DO lev = nlevels, 2,  - 1
        lay = lev - 1
        opdp_path_ad(lev - 1, j) = opdp_path_ad(lev - 1, j) + opdp_path_ad(lev, j)
        opticaldepth_ad(lay, j)  = opticaldepth_ad(lay, j) + opdp_path_ad(lev, j)
        opdp_path_ad(lev, j)     = 0.0_jprb
      ENDDO
      opdp_path_ad(1, j) = 0.0_jprb
    ENDIF
  ENDDO
  DO j = 1, nchannels
    IF (chanflag(j)) THEN
      WHERE (opdp_ref(:,j) > 0.0_jprb) opticaldepth_ad(:,j) = 0.0_jprb
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
      debye_prof_ad => aux_ad%debye_prof(:, :, prof)
      IF (pred_info%ncloud >= 1) THEN
        DO lay = nlayers, 1,  - 1
          lev = lay + 1
          IF (lev >= coef%mwcldtop) THEN
! Repeat direct code
            zf = coef%frequency_ghz(chan)
            zf_sq = zf * zf
            z1_sq = debye_prof(1, lev) * debye_prof(1, lev)
            z2_sq = debye_prof(2, lev) * debye_prof(2, lev)
            z34_dif                      = debye_prof(3, lev) - debye_prof(4, lev)
            z45_dif                      = debye_prof(4, lev) - debye_prof(5, lev)
            z1f_sq_z1_sq                 = zf_sq + z1_sq
            z2f_sq_z2_sq                 = zf_sq + z2_sq
            z1_div = 1.0_jprb / z1f_sq_z1_sq
            z2_div = 1.0_jprb / z2f_sq_z2_sq
            z1_den = z34_dif * z1_div
            z2_den = z45_dif * z2_div
            zastar = debye_prof(3, lev) - zf_sq * (z1_den + z2_den)
            z1_prod                      = z34_dif * debye_prof(1, lev)
            z2_prod                      = z1_prod * z1_div
            z3_prod                      = z45_dif * debye_prof(2, lev)
            z4_prod                      = z3_prod * z2_div
            zbstar =  - zf * (z2_prod + z4_prod)
            zbstar_sq                    = zbstar * zbstar
            za2star                      = zastar + 2.0_jprb
            za2star_sq                   = za2star * za2star
            zdiv = za2star_sq + zbstar_sq
            zgstar =  - 3.0_jprb * zbstar / zdiv
! Now compute Adjoint code
!opticaldepth_ad(lay,j)= opticaldepth_ad(lay,j)
            zgstar_ad                    = opticaldepth_ad(lay, j) * ( - 1.5_jprb * zf * predictors%clw(lay, prof))
            predictors_ad%clw(lay, prof) =      &
              & predictors_ad%clw(lay, prof) + opticaldepth_ad(lay, j) * ( - 1.5_jprb * zf * zgstar)
            zbstar_ad                    =  - 3.0_jprb * zgstar_ad / zdiv
            zdiv_ad                      = 3.0_jprb * zgstar_ad * zbstar / (zdiv * zdiv)
!zgstar_ad =  0.
            za2star_sq_ad                = zdiv_ad
            zbstar_sq_ad                 = zdiv_ad
!zdiv_ad       = 0.
            za2star_ad                   = 2.0_jprb * za2star * za2star_sq_ad
!za2star_sq_ad = 0.
            zastar_ad                    = za2star_ad
!za2star_ad    = 0.
            zbstar_ad                    = zbstar_ad + 2.0_jprb * zbstar * zbstar_sq_ad
!zbstar_sq_ad  = 0.
            z2_prod_ad                   =  - zf * zbstar_ad
            z4_prod_ad                   =  - zf * zbstar_ad
!zbstar_ad     = 0.
            z3_prod_ad                   = z2_div * z4_prod_ad
            z2_div_ad                    = z3_prod * z4_prod_ad
!z4_prod_ad    = 0.
            z45_dif_ad                   = debye_prof(2, lev) * z3_prod_ad
            debye_prof_ad(2, lev)        = debye_prof_ad(2, lev) + z45_dif * z3_prod_ad
!z3_prod_ad           = 0.
            z1_prod_ad                   = z1_div * z2_prod_ad
            z1_div_ad                    = z1_prod * z2_prod_ad
!z2_prod_ad    = 0.
            z34_dif_ad                   = debye_prof(1, lev) * z1_prod_ad
            debye_prof_ad(1, lev)        = debye_prof_ad(1, lev) + z34_dif * z1_prod_ad
!z1_prod_ad           = 0.
            debye_prof_ad(3, lev)        = debye_prof_ad(3, lev) + zastar_ad
            z1_den_ad                    =  - zf_sq * zastar_ad
            z2_den_ad                    =  - zf_sq * zastar_ad
!zastar_ad            = 0.
            z2_div_ad                    = z2_div_ad + z45_dif * z2_den_ad
            z45_dif_ad                   = z45_dif_ad + z2_div * z2_den_ad
!z2_den_ad     = 0.
            z1_div_ad                    = z1_div_ad + z34_dif * z1_den_ad
            z34_dif_ad                   = z34_dif_ad + z1_div * z1_den_ad
!z1_den_ad     = 0.
            z2f_sq_z2_sq_ad              =  - z2_div_ad / (z2f_sq_z2_sq * z2f_sq_z2_sq)
!z2_div_ad       = 0.
            z1f_sq_z1_sq_ad              =  - z1_div_ad / (z1f_sq_z1_sq * z1f_sq_z1_sq)
!z1_div_ad       = 0.
            z2_sq_ad                     = z2f_sq_z2_sq_ad
!z2f_sq_z2_sq_ad = 0.
            z1_sq_ad                     = z1f_sq_z1_sq_ad
!z1f_sq_z1_sq_ad = 0.
            debye_prof_ad(4, lev)        = debye_prof_ad(4, lev) + z45_dif_ad
            debye_prof_ad(5, lev)        = debye_prof_ad(5, lev) - z45_dif_ad
!z45_dif_ad           = 0.
            debye_prof_ad(3, lev)        = debye_prof_ad(3, lev) + z34_dif_ad
            debye_prof_ad(4, lev)        = debye_prof_ad(4, lev) - z34_dif_ad
!z34_dif_ad           = 0.
            debye_prof_ad(2, lev)        = debye_prof_ad(2, lev) + z2_sq_ad * 2.0_jprb * debye_prof(2, lev)
!z2_sq_ad             = 0.
            debye_prof_ad(1, lev)        = debye_prof_ad(1, lev) + z1_sq_ad * 2.0_jprb * debye_prof(1, lev)
!z1_sq_ad             = 0.
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
        prof = chanprof(j)%prof
        DO lay = nlayers, 1,  - 1
          predictors_ad%ch4(:, lay, prof) =      &
            & predictors_ad%ch4(:, lay, prof) + fast_coef%ch4(lay, chan, :) * opticaldepth_ad(lay, j)
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
            predictors_ad%co(1:11, lay, prof) =      &
              & predictors_ad%co(1:11, lay, prof) + fast_coef%co(lay, chan, 1:11) * opticaldepth_ad(lay, j)
          ENDDO
        ELSE
          DO lay = nlayers, 1,  - 1
            predictors_ad%co(1:10, lay, prof) =      &
              & predictors_ad%co(1:10, lay, prof) + fast_coef%co(lay, chan, 1:10) * opticaldepth_ad(lay, j)
            predictors_ad%co(12, lay, prof)   =      &
              & predictors_ad%co(12, lay, prof) + fast_coef%co(lay, chan, 11) * opticaldepth_ad(lay, j)
            predictors_ad%co(13, lay, prof)   =      &
              & predictors_ad%co(13, lay, prof) + fast_coef%co(lay, chan, 12) * opticaldepth_ad(lay, j)
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
              predictors_ad%ch4(10, lay, prof) =      &
                & predictors_ad%ch4(10, lay, prof) + fast_coef%n2o(lay, chan, 12) * opticaldepth_ad(lay, j)
              predictors_ad%ch4(1, lay, prof)  =      &
                & predictors_ad%ch4(1, lay, prof) + fast_coef%n2o(lay, chan, 11) * opticaldepth_ad(lay, j)
            ENDIF
            predictors_ad%n2o(1:10, lay, prof) =      &
              & predictors_ad%n2o(1:10, lay, prof) + fast_coef%n2o(lay, chan, 1:10) * opticaldepth_ad(lay, j)
          ENDDO
        ELSE IF (chanx >= rttov9_wv1995_00 .AND. chanx < rttov9_wv2000_00) THEN
          DO lay = nlayers, 1,  - 1
            IF (coef%nco > 0) THEN
              predictors_ad%co(8, lay, prof) = predictors_ad%co(8, lay, prof) -                                &
                & fast_coef%n2o(lay, chan, 12) * opticaldepth_ad(lay, j) * predictors%co(1, lay, prof) ** 3 /  &
                & predictors%co(8, lay, prof) ** 2
              predictors_ad%co(1, lay, prof) = predictors_ad%co(1, lay, prof) +                                    &
                & fast_coef%n2o(lay, chan, 12) * opticaldepth_ad(lay, j) * 3 * predictors%co(1, lay, prof) ** 2 /  &
                & predictors%co(8, lay, prof)
              predictors_ad%co(1, lay, prof) =      &
                & predictors_ad%co(1, lay, prof) + fast_coef%n2o(lay, chan, 11) * opticaldepth_ad(lay, j)
            ENDIF
            predictors_ad%n2o(1:10, lay, prof) =      &
              & predictors_ad%n2o(1:10, lay, prof) + fast_coef%n2o(lay, chan, 1:10) * opticaldepth_ad(lay, j)
          ENDDO
        ELSE IF (chanx < rttov9_wv2000_00) THEN
          DO lay = nlayers, 1,  - 1
            predictors_ad%n2o(1:10, lay, prof) =      &
              & predictors_ad%n2o(1:10, lay, prof) + fast_coef%n2o(lay, chan, 1:10) * opticaldepth_ad(lay, j)
          ENDDO
        ELSE IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2250_00) THEN
          DO lay = nlayers, 1,  - 1
            predictors_ad%n2o(13, lay, prof) =      &
              & predictors_ad%n2o(13, lay, prof) + fast_coef%n2o(lay, chan, 14) * opticaldepth_ad(lay, j)
            predictors_ad%n2o(12, lay, prof) =      &
              & predictors_ad%n2o(12, lay, prof) + fast_coef%n2o(lay, chan, 13) * opticaldepth_ad(lay, j)
            IF (coef%nco > 0) THEN
              predictors_ad%co(8, lay, prof) = predictors_ad%co(8, lay, prof) -                                &
                & fast_coef%n2o(lay, chan, 12) * opticaldepth_ad(lay, j) * predictors%co(1, lay, prof) ** 3 /  &
                & predictors%co(8, lay, prof) ** 2
              predictors_ad%co(1, lay, prof) = predictors_ad%co(1, lay, prof) +                                    &
                & fast_coef%n2o(lay, chan, 12) * opticaldepth_ad(lay, j) * 3 * predictors%co(1, lay, prof) ** 2 /  &
                & predictors%co(8, lay, prof)
              predictors_ad%co(1, lay, prof) =      &
                & predictors_ad%co(1, lay, prof) + fast_coef%n2o(lay, chan, 11) * opticaldepth_ad(lay, j)
            ENDIF
            predictors_ad%n2o(10, lay, prof)  =      &
              & predictors_ad%n2o(10, lay, prof) + fast_coef%n2o(lay, chan, 10) * opticaldepth_ad(lay, j)
            predictors_ad%n2o(11, lay, prof)  =      &
              & predictors_ad%n2o(11, lay, prof) + fast_coef%n2o(lay, chan, 9) * opticaldepth_ad(lay, j)
            predictors_ad%n2o(1:8, lay, prof) =      &
              & predictors_ad%n2o(1:8, lay, prof) + fast_coef%n2o(lay, chan, 1:8) * opticaldepth_ad(lay, j)
          ENDDO
        ELSE IF (chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2295_25) THEN
          DO lay = nlayers, 1,  - 1
            IF (coef%nco > 0) THEN
              predictors_ad%co(8, lay, prof) = predictors_ad%co(8, lay, prof) -                                &
                & fast_coef%n2o(lay, chan, 12) * opticaldepth_ad(lay, j) * predictors%co(1, lay, prof) ** 3 /  &
                & predictors%co(8, lay, prof) ** 2
              predictors_ad%co(1, lay, prof) = predictors_ad%co(1, lay, prof) +                                    &
                & fast_coef%n2o(lay, chan, 12) * opticaldepth_ad(lay, j) * 3 * predictors%co(1, lay, prof) ** 2 /  &
                & predictors%co(8, lay, prof)
              predictors_ad%co(1, lay, prof) =      &
                & predictors_ad%co(1, lay, prof) + fast_coef%n2o(lay, chan, 11) * opticaldepth_ad(lay, j)
            ENDIF
            predictors_ad%n2o(1:10, lay, prof) =      &
              & predictors_ad%n2o(1:10, lay, prof) + fast_coef%n2o(lay, chan, 1:10) * opticaldepth_ad(lay, j)
          ENDDO
        ELSE IF (chanx >= rttov9_wv2295_25 .AND. chanx < rttov9_wv2380_25) THEN
          DO lay = nlayers, 1,  - 1
            predictors_ad%n2o(1:10, lay, prof) =      &
              & predictors_ad%n2o(1:10, lay, prof) + fast_coef%n2o(lay, chan, 1:10) * opticaldepth_ad(lay, j)
          ENDDO
        ELSE   ! chanx >= rttov9_wv2380_25
          DO lay = nlayers, 1,  - 1
            predictors_ad%n2o(13, lay, prof)  =      &
              & predictors_ad%n2o(13, lay, prof) + fast_coef%n2o(lay, chan, 12) * opticaldepth_ad(lay, j)
            predictors_ad%n2o(12, lay, prof)  =      &
              & predictors_ad%n2o(12, lay, prof) + fast_coef%n2o(lay, chan, 11) * opticaldepth_ad(lay, j)
            predictors_ad%n2o(10, lay, prof)  =      &
              & predictors_ad%n2o(10, lay, prof) + fast_coef%n2o(lay, chan, 10) * opticaldepth_ad(lay, j)
            predictors_ad%n2o(11, lay, prof)  =      &
              & predictors_ad%n2o(11, lay, prof) + fast_coef%n2o(lay, chan, 9) * opticaldepth_ad(lay, j)
            predictors_ad%n2o(1:8, lay, prof) =      &
              & predictors_ad%n2o(1:8, lay, prof) + fast_coef%n2o(lay, chan, 1:8) * opticaldepth_ad(lay, j)
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
        prof  = chanprof(j)%prof
        chanx = coef%ff_cwn(chan)
        IF (chanx >= rttov9_wv0690_50 .AND. chanx < rttov9_wv1100_25) THEN
          DO lay = nlayers, 1,  - 1
            predictors_ad%co2(9, lay, prof) =      &
              & predictors_ad%co2(9, lay, prof) + fast_coef%co2(lay, chan, 9) * opticaldepth_ad(lay, j)
            predictors_ad%co2(8, lay, prof) =      &
              & predictors_ad%co2(8, lay, prof) + fast_coef%co2(lay, chan, 8) * opticaldepth_ad(lay, j)
            predictors_ad%co2(7, lay, prof) = predictors_ad%co2(7, lay, prof) +      &
              & fast_coef%co2(lay, chan, 7) * 2 * predictors%co2(7, lay, prof) * opticaldepth_ad(lay, j)
            predictors_ad%co2(6, lay, prof) = predictors_ad%co2(6, lay, prof) +      &
              & fast_coef%co2(lay, chan, 6) * 2 * predictors%co2(6, lay, prof) * opticaldepth_ad(lay, j)
            predictors_ad%co2(5, lay, prof) =      &
              & predictors_ad%co2(5, lay, prof) + fast_coef%co2(lay, chan, 5) * opticaldepth_ad(lay, j)
            predictors_ad%co2(3, lay, prof) = predictors_ad%co2(3, lay, prof) +      &
              & fast_coef%co2(lay, chan, 4) * 2 * predictors%co2(3, lay, prof) * opticaldepth_ad(lay, j)
            predictors_ad%co2(5, lay, prof) = predictors_ad%co2(5, lay, prof) +      &
              & fast_coef%co2(lay, chan, 3) * predictors%co2(6, lay, prof) ** 2 * opticaldepth_ad(lay, j)
            predictors_ad%co2(6, lay, prof) = predictors_ad%co2(6, lay, prof) +                             &
              & fast_coef%co2(lay, chan, 3) * 2 * predictors%co2(6, lay, prof) * predictors%co2(5, lay, prof) *  &
              & opticaldepth_ad(lay, j)
            predictors_ad%co2(2, lay, prof) =      &
              & predictors_ad%co2(2, lay, prof) + fast_coef%co2(lay, chan, 2) * opticaldepth_ad(lay, j)
            predictors_ad%co2(1, lay, prof) =      &
              & predictors_ad%co2(1, lay, prof) + fast_coef%co2(lay, chan, 1) * opticaldepth_ad(lay, j)
          ENDDO
        ELSE IF (chanx >= rttov9_wv1995_00 .AND. chanx < rttov9_wv2000_00) THEN
          DO lay = nlayers, 1,  - 1
            predictors_ad%co2(10, lay, prof) =      &
              & predictors_ad%co2(10, lay, prof) + fast_coef%co2(lay, chan, 11) * opticaldepth_ad(lay, j)
            IF (coef%nco > 0) THEN
              predictors_ad%co(1, lay, prof) =      &
                & predictors_ad%co(1, lay, prof) + fast_coef%co2(lay, chan, 10) * opticaldepth_ad(lay, j)
            ENDIF
            predictors_ad%co2(1:9, lay, prof) =      &
              & predictors_ad%co2(1:9, lay, prof) + fast_coef%co2(lay, chan, 1:9) * opticaldepth_ad(lay, j)
          ENDDO
        ELSE IF (chanx < rttov9_wv2000_00) THEN
          DO lay = nlayers, 1,  - 1
            predictors_ad%co2(1:9, lay, prof) =      &
              & predictors_ad%co2(1:9, lay, prof) + fast_coef%co2(lay, chan, 1:9) * opticaldepth_ad(lay, j)
          ENDDO
        ELSE IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2250_00) THEN
          DO lay = nlayers, 1,  - 1
            predictors_ad%co2(15, lay, prof) =      &
              & predictors_ad%co2(15, lay, prof) + fast_coef%co2(lay, chan, 14) * opticaldepth_ad(lay, j)
            predictors_ad%co2(14, lay, prof) =      &
              & predictors_ad%co2(14, lay, prof) + fast_coef%co2(lay, chan, 13) * opticaldepth_ad(lay, j)
            predictors_ad%co2(13, lay, prof) =      &
              & predictors_ad%co2(13, lay, prof) + fast_coef%co2(lay, chan, 12) * opticaldepth_ad(lay, j)
            predictors_ad%co2(12, lay, prof) =      &
              & predictors_ad%co2(12, lay, prof) + fast_coef%co2(lay, chan, 11) * opticaldepth_ad(lay, j)
            predictors_ad%co2(11, lay, prof) =      &
              & predictors_ad%co2(11, lay, prof) + fast_coef%co2(lay, chan, 10) * opticaldepth_ad(lay, j)
            predictors_ad%co2(10, lay, prof) =      &
              & predictors_ad%co2(10, lay, prof) + fast_coef%co2(lay, chan, 9) * opticaldepth_ad(lay, j)
            IF (coef%nco > 0) THEN
              predictors_ad%co(1, lay, prof) =      &
                & predictors_ad%co(1, lay, prof) + fast_coef%co2(lay, chan, 8) * opticaldepth_ad(lay, j)
            ENDIF
            predictors_ad%co2(8, lay, prof)   =      &
              & predictors_ad%co2(8, lay, prof) + fast_coef%co2(lay, chan, 7) * opticaldepth_ad(lay, j)
            predictors_ad%co2(7, lay, prof)   =      &
              & predictors_ad%co2(7, lay, prof) + fast_coef%co2(lay, chan, 6) * opticaldepth_ad(lay, j)
            predictors_ad%co2(1:5, lay, prof) =      &
              & predictors_ad%co2(1:5, lay, prof) + fast_coef%co2(lay, chan, 1:5) * opticaldepth_ad(lay, j)
          ENDDO
        ELSE IF (chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2295_25) THEN
          DO lay = nlayers, 1,  - 1
            predictors_ad%co2(10, lay, prof) =      &
              & predictors_ad%co2(10, lay, prof) + fast_coef%co2(lay, chan, 11) * opticaldepth_ad(lay, j)
            IF (coef%nco > 0) THEN
              predictors_ad%co(1, lay, prof) =      &
                & predictors_ad%co(1, lay, prof) + fast_coef%co2(lay, chan, 10) * opticaldepth_ad(lay, j)
            ENDIF
            predictors_ad%co2(1:9, lay, prof) =      &
              & predictors_ad%co2(1:9, lay, prof) + fast_coef%co2(lay, chan, 1:9) * opticaldepth_ad(lay, j)
          ENDDO
        ELSE IF (chanx >= rttov9_wv2295_25 .AND. chanx < rttov9_wv2380_25) THEN
          DO lay = nlayers, 1,  - 1
            predictors_ad%co2(1:10, lay, prof) =      &
              & predictors_ad%co2(1:10, lay, prof) + fast_coef%co2(lay, chan, 1:10) * opticaldepth_ad(lay, j)
          ENDDO
        ELSE   ! chanx >= rttov9_wv2380_25
          DO lay = nlayers, 1,  - 1
            predictors_ad%co2(15, lay, prof)  =      &
              & predictors_ad%co2(15, lay, prof) + fast_coef%co2(lay, chan, 13) * opticaldepth_ad(lay, j)
            predictors_ad%co2(14, lay, prof)  =      &
              & predictors_ad%co2(14, lay, prof) + fast_coef%co2(lay, chan, 12) * opticaldepth_ad(lay, j)
            predictors_ad%co2(13, lay, prof)  =      &
              & predictors_ad%co2(13, lay, prof) + fast_coef%co2(lay, chan, 11) * opticaldepth_ad(lay, j)
            predictors_ad%co2(12, lay, prof)  =      &
              & predictors_ad%co2(12, lay, prof) + fast_coef%co2(lay, chan, 10) * opticaldepth_ad(lay, j)
            predictors_ad%co2(11, lay, prof)  =      &
              & predictors_ad%co2(11, lay, prof) + fast_coef%co2(lay, chan, 9) * opticaldepth_ad(lay, j)
            predictors_ad%co2(10, lay, prof)  =      &
              & predictors_ad%co2(10, lay, prof) + fast_coef%co2(lay, chan, 8) * opticaldepth_ad(lay, j)
            predictors_ad%co2(8, lay, prof)   =      &
              & predictors_ad%co2(8, lay, prof) + fast_coef%co2(lay, chan, 7) * opticaldepth_ad(lay, j)
            predictors_ad%co2(7, lay, prof)   =      &
              & predictors_ad%co2(7, lay, prof) + fast_coef%co2(lay, chan, 6) * opticaldepth_ad(lay, j)
            predictors_ad%co2(1:5, lay, prof) =      &
              & predictors_ad%co2(1:5, lay, prof) + fast_coef%co2(lay, chan, 1:5) * opticaldepth_ad(lay, j)
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
        prof = chanprof(j)%prof
        DO lay = nlayers, 1,  - 1
          predictors_ad%wvcont(:, lay, prof) =      &
            & predictors_ad%wvcont(:, lay, prof) + fast_coef%wvcont(lay, chan, :) * opticaldepth_ad(lay, j)
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
            predictors_ad%ozone(11, lay, prof) =      &
              & predictors_ad%ozone(11, lay, prof) + fast_coef%ozone(lay, chan, 11) * opticaldepth_ad(lay, j)
            predictors_ad%ozone(10, lay, prof) =      &
              & predictors_ad%ozone(10, lay, prof) + fast_coef%ozone(lay, chan, 10) * opticaldepth_ad(lay, j)
            predictors_ad%ozone(9, lay, prof)  =      &
              & predictors_ad%ozone(9, lay, prof) + fast_coef%ozone(lay, chan, 9) * opticaldepth_ad(lay, j)
            predictors_ad%ozone(8, lay, prof)  =      &
              & predictors_ad%ozone(8, lay, prof) + fast_coef%ozone(lay, chan, 8) * opticaldepth_ad(lay, j)
            predictors_ad%ozone(7, lay, prof)  =      &
              & predictors_ad%ozone(7, lay, prof) + fast_coef%ozone(lay, chan, 7) * opticaldepth_ad(lay, j)
            predictors_ad%ozone(6, lay, prof)  =      &
              & predictors_ad%ozone(6, lay, prof) + fast_coef%ozone(lay, chan, 6) * opticaldepth_ad(lay, j)
            predictors_ad%ozone(5, lay, prof)  =      &
              & predictors_ad%ozone(5, lay, prof) + fast_coef%ozone(lay, chan, 5) * opticaldepth_ad(lay, j)
            predictors_ad%ozone(4, lay, prof)  =      &
              & predictors_ad%ozone(4, lay, prof) + fast_coef%ozone(lay, chan, 4) * opticaldepth_ad(lay, j)
            predictors_ad%ozone(3, lay, prof)  =      &
              & predictors_ad%ozone(3, lay, prof) + fast_coef%ozone(lay, chan, 3) * opticaldepth_ad(lay, j)
            predictors_ad%ozone(2, lay, prof)  =      &
              & predictors_ad%ozone(2, lay, prof) + fast_coef%ozone(lay, chan, 2) * opticaldepth_ad(lay, j)
            predictors_ad%ozone(1, lay, prof)  =      &
              & predictors_ad%ozone(1, lay, prof) + fast_coef%ozone(lay, chan, 1) * opticaldepth_ad(lay, j)
          ENDDO
        ELSE
          DO lay = nlayers, 1,  - 1
            predictors_ad%ozone(15, lay, prof) =      &
              & predictors_ad%ozone(15, lay, prof) + fast_coef%ozone(lay, chan, 13) * opticaldepth_ad(lay, j)
            predictors_ad%ozone(14, lay, prof) =      &
              & predictors_ad%ozone(14, lay, prof) + fast_coef%ozone(lay, chan, 12) * opticaldepth_ad(lay, j)
            predictors_ad%ozone(11, lay, prof) =      &
              & predictors_ad%ozone(11, lay, prof) + fast_coef%ozone(lay, chan, 11) * opticaldepth_ad(lay, j)
            predictors_ad%ozone(10, lay, prof) =      &
              & predictors_ad%ozone(10, lay, prof) + fast_coef%ozone(lay, chan, 10) * opticaldepth_ad(lay, j)
            predictors_ad%ozone(9, lay, prof)  =      &
              & predictors_ad%ozone(9, lay, prof) + fast_coef%ozone(lay, chan, 9) * opticaldepth_ad(lay, j)
            predictors_ad%ozone(13, lay, prof) =      &
              & predictors_ad%ozone(13, lay, prof) + fast_coef%ozone(lay, chan, 8) * opticaldepth_ad(lay, j)
            predictors_ad%ozone(7, lay, prof)  =      &
              & predictors_ad%ozone(7, lay, prof) + fast_coef%ozone(lay, chan, 7) * opticaldepth_ad(lay, j)
            predictors_ad%ozone(6, lay, prof)  =      &
              & predictors_ad%ozone(6, lay, prof) + fast_coef%ozone(lay, chan, 6) * opticaldepth_ad(lay, j)
            predictors_ad%ozone(5, lay, prof)  =      &
              & predictors_ad%ozone(5, lay, prof) + fast_coef%ozone(lay, chan, 5) * opticaldepth_ad(lay, j)
            predictors_ad%ozone(4, lay, prof)  =      &
              & predictors_ad%ozone(4, lay, prof) + fast_coef%ozone(lay, chan, 4) * opticaldepth_ad(lay, j)
            predictors_ad%ozone(12, lay, prof) =      &
              & predictors_ad%ozone(12, lay, prof) + fast_coef%ozone(lay, chan, 3) * opticaldepth_ad(lay, j)
            predictors_ad%ozone(2, lay, prof)  =      &
              & predictors_ad%ozone(2, lay, prof) + fast_coef%ozone(lay, chan, 2) * opticaldepth_ad(lay, j)
            predictors_ad%ozone(1, lay, prof)  =      &
              & predictors_ad%ozone(1, lay, prof) + fast_coef%ozone(lay, chan, 1) * opticaldepth_ad(lay, j)
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
            predictors_ad%ch4(1, lay, prof) = predictors_ad%ch4(1, lay, prof) +      &
              & opticaldepth_ad(lay, j) * fast_coef%watervapour(lay, chan, 15) * predictors%ch4(3, lay, prof)
            predictors_ad%ch4(3, lay, prof) = predictors_ad%ch4(3, lay, prof) +      &
              & opticaldepth_ad(lay, j) * fast_coef%watervapour(lay, chan, 15) * predictors%ch4(1, lay, prof)
            predictors_ad%ch4(1, lay, prof) =      &
              & predictors_ad%ch4(1, lay, prof) + opticaldepth_ad(lay, j) * fast_coef%watervapour(lay, chan, 14)
          ENDIF
          predictors_ad%watervapour(13, lay, prof) =      &
            & predictors_ad%watervapour(13, lay, prof) + fast_coef%watervapour(lay, chan, 13) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(12, lay, prof) =      &
            & predictors_ad%watervapour(12, lay, prof) + fast_coef%watervapour(lay, chan, 12) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(11, lay, prof) =      &
            & predictors_ad%watervapour(11, lay, prof) + fast_coef%watervapour(lay, chan, 11) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(10, lay, prof) =      &
            & predictors_ad%watervapour(10, lay, prof) + fast_coef%watervapour(lay, chan, 10) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(9, lay, prof)  =      &
            & predictors_ad%watervapour(9, lay, prof) + fast_coef%watervapour(lay, chan, 9) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(8, lay, prof)  =      &
            & predictors_ad%watervapour(8, lay, prof) + fast_coef%watervapour(lay, chan, 8) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(7, lay, prof)  =      &
            & predictors_ad%watervapour(7, lay, prof) + fast_coef%watervapour(lay, chan, 7) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(6, lay, prof)  =      &
            & predictors_ad%watervapour(6, lay, prof) + fast_coef%watervapour(lay, chan, 6) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(5, lay, prof)  =      &
            & predictors_ad%watervapour(5, lay, prof) + fast_coef%watervapour(lay, chan, 5) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(4, lay, prof)  =      &
            & predictors_ad%watervapour(4, lay, prof) + fast_coef%watervapour(lay, chan, 4) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(3, lay, prof)  =      &
            & predictors_ad%watervapour(3, lay, prof) + fast_coef%watervapour(lay, chan, 3) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(2, lay, prof)  =      &
            & predictors_ad%watervapour(2, lay, prof) + fast_coef%watervapour(lay, chan, 2) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(1, lay, prof)  =      &
            & predictors_ad%watervapour(1, lay, prof) + fast_coef%watervapour(lay, chan, 1) * opticaldepth_ad(lay, j)
        ENDDO
      ELSE IF (chanx >= rttov9_wv1900_25 .AND. chanx < rttov9_wv1995_00) THEN
        DO lay = nlayers, 1,  - 1
          IF (coef%nco2 > 0) THEN
            predictors_ad%co2(1, lay, prof) =      &
              & predictors_ad%co2(1, lay, prof) + opticaldepth_ad(lay, j) * fast_coef%watervapour(lay, chan, 14)
          ENDIF
          predictors_ad%watervapour(13, lay, prof) =      &
            & predictors_ad%watervapour(13, lay, prof) + fast_coef%watervapour(lay, chan, 13) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(12, lay, prof) =      &
            & predictors_ad%watervapour(12, lay, prof) + fast_coef%watervapour(lay, chan, 12) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(11, lay, prof) =      &
            & predictors_ad%watervapour(11, lay, prof) + fast_coef%watervapour(lay, chan, 11) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(10, lay, prof) =      &
            & predictors_ad%watervapour(10, lay, prof) + fast_coef%watervapour(lay, chan, 10) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(9, lay, prof)  =      &
            & predictors_ad%watervapour(9, lay, prof) + fast_coef%watervapour(lay, chan, 9) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(8, lay, prof)  =      &
            & predictors_ad%watervapour(8, lay, prof) + fast_coef%watervapour(lay, chan, 8) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(7, lay, prof)  =      &
            & predictors_ad%watervapour(7, lay, prof) + fast_coef%watervapour(lay, chan, 7) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(6, lay, prof)  =      &
            & predictors_ad%watervapour(6, lay, prof) + fast_coef%watervapour(lay, chan, 6) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(5, lay, prof)  =      &
            & predictors_ad%watervapour(5, lay, prof) + fast_coef%watervapour(lay, chan, 5) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(4, lay, prof)  =      &
            & predictors_ad%watervapour(4, lay, prof) + fast_coef%watervapour(lay, chan, 4) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(3, lay, prof)  =      &
            & predictors_ad%watervapour(3, lay, prof) + fast_coef%watervapour(lay, chan, 3) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(2, lay, prof)  =      &
            & predictors_ad%watervapour(2, lay, prof) + fast_coef%watervapour(lay, chan, 2) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(1, lay, prof)  =      &
            & predictors_ad%watervapour(1, lay, prof) + fast_coef%watervapour(lay, chan, 1) * opticaldepth_ad(lay, j)
        ENDDO
      ELSE IF (chanx >= rttov9_wv1995_00 .AND. chanx < rttov9_wv2000_00) THEN
        DO lay = nlayers, 1,  - 1
          IF (coef%nco > 0) THEN
            predictors_ad%co(1, lay, prof) =      &
              & predictors_ad%co(1, lay, prof) + fast_coef%watervapour(lay, chan, 15) * opticaldepth_ad(lay, j)
          ENDIF
          IF (coef%nco2 > 0) THEN
            predictors_ad%co2(1, lay, prof) =      &
              & predictors_ad%co2(1, lay, prof) + fast_coef%watervapour(lay, chan, 14) * opticaldepth_ad(lay, j)
          ENDIF
          predictors_ad%watervapour(13, lay, prof) =      &
            & predictors_ad%watervapour(13, lay, prof) + fast_coef%watervapour(lay, chan, 13) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(12, lay, prof) =      &
            & predictors_ad%watervapour(12, lay, prof) + fast_coef%watervapour(lay, chan, 12) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(11, lay, prof) =      &
            & predictors_ad%watervapour(11, lay, prof) + fast_coef%watervapour(lay, chan, 11) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(10, lay, prof) =      &
            & predictors_ad%watervapour(10, lay, prof) + fast_coef%watervapour(lay, chan, 10) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(9, lay, prof)  =      &
            & predictors_ad%watervapour(9, lay, prof) + fast_coef%watervapour(lay, chan, 9) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(8, lay, prof)  =      &
            & predictors_ad%watervapour(8, lay, prof) + fast_coef%watervapour(lay, chan, 8) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(7, lay, prof)  =      &
            & predictors_ad%watervapour(7, lay, prof) + fast_coef%watervapour(lay, chan, 7) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(6, lay, prof)  =      &
            & predictors_ad%watervapour(6, lay, prof) + fast_coef%watervapour(lay, chan, 6) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(5, lay, prof)  =      &
            & predictors_ad%watervapour(5, lay, prof) + fast_coef%watervapour(lay, chan, 5) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(4, lay, prof)  =      &
            & predictors_ad%watervapour(4, lay, prof) + fast_coef%watervapour(lay, chan, 4) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(3, lay, prof)  =      &
            & predictors_ad%watervapour(3, lay, prof) + fast_coef%watervapour(lay, chan, 3) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(2, lay, prof)  =      &
            & predictors_ad%watervapour(2, lay, prof) + fast_coef%watervapour(lay, chan, 2) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(1, lay, prof)  =      &
            & predictors_ad%watervapour(1, lay, prof) + fast_coef%watervapour(lay, chan, 1) * opticaldepth_ad(lay, j)
        ENDDO
      ELSE IF (chanx < rttov9_wv2000_00) THEN
        DO lay = nlayers, 1,  - 1
          predictors_ad%watervapour(13, lay, prof) =      &
            & predictors_ad%watervapour(13, lay, prof) + fast_coef%watervapour(lay, chan, 13) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(12, lay, prof) =      &
            & predictors_ad%watervapour(12, lay, prof) + fast_coef%watervapour(lay, chan, 12) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(11, lay, prof) =      &
            & predictors_ad%watervapour(11, lay, prof) + fast_coef%watervapour(lay, chan, 11) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(10, lay, prof) =      &
            & predictors_ad%watervapour(10, lay, prof) + fast_coef%watervapour(lay, chan, 10) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(9, lay, prof)  =      &
            & predictors_ad%watervapour(9, lay, prof) + fast_coef%watervapour(lay, chan, 9) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(8, lay, prof)  =      &
            & predictors_ad%watervapour(8, lay, prof) + fast_coef%watervapour(lay, chan, 8) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(7, lay, prof)  =      &
            & predictors_ad%watervapour(7, lay, prof) + fast_coef%watervapour(lay, chan, 7) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(6, lay, prof)  =      &
            & predictors_ad%watervapour(6, lay, prof) + fast_coef%watervapour(lay, chan, 6) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(5, lay, prof)  =      &
            & predictors_ad%watervapour(5, lay, prof) + fast_coef%watervapour(lay, chan, 5) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(4, lay, prof)  =      &
            & predictors_ad%watervapour(4, lay, prof) + fast_coef%watervapour(lay, chan, 4) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(3, lay, prof)  =      &
            & predictors_ad%watervapour(3, lay, prof) + fast_coef%watervapour(lay, chan, 3) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(2, lay, prof)  =      &
            & predictors_ad%watervapour(2, lay, prof) + fast_coef%watervapour(lay, chan, 2) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(1, lay, prof)  =      &
            & predictors_ad%watervapour(1, lay, prof) + fast_coef%watervapour(lay, chan, 1) * opticaldepth_ad(lay, j)
        ENDDO
      ELSE IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2250_00) THEN
        IF (coef%ss_val_chn(chan) > 0) THEN
          DO lay = nlayers, 1,  - 1
            IF (coef%nco > 0) THEN
              predictors_ad%co(1, lay, prof) =      &
                & predictors_ad%co(1, lay, prof) + fast_coef%watervapour(lay, chan, 15) * opticaldepth_ad(lay, j)
            ENDIF
            IF (coef%nco2 > 0) THEN
              predictors_ad%co2(1, lay, prof) =      &
                & predictors_ad%co2(1, lay, prof) + fast_coef%watervapour(lay, chan, 14) * opticaldepth_ad(lay, j)
            ENDIF
            predictors_ad%watervapour(17, lay, prof) =      &
              & predictors_ad%watervapour(17, lay, prof) + fast_coef%watervapour(lay, chan, 13) * opticaldepth_ad(lay, j)
            predictors_ad%watervapour(16, lay, prof) =      &
              & predictors_ad%watervapour(16, lay, prof) + fast_coef%watervapour(lay, chan, 12) * opticaldepth_ad(lay, j)
            predictors_ad%watervapour(11, lay, prof) =      &
              & predictors_ad%watervapour(11, lay, prof) + fast_coef%watervapour(lay, chan, 11) * opticaldepth_ad(lay, j)
            predictors_ad%watervapour(10, lay, prof) =      &
              & predictors_ad%watervapour(10, lay, prof) + fast_coef%watervapour(lay, chan, 10) * opticaldepth_ad(lay, j)
            predictors_ad%watervapour(15, lay, prof) =      &
              & predictors_ad%watervapour(15, lay, prof) + fast_coef%watervapour(lay, chan, 9) * opticaldepth_ad(lay, j)
            predictors_ad%watervapour(14, lay, prof) =      &
              & predictors_ad%watervapour(14, lay, prof) + fast_coef%watervapour(lay, chan, 8) * opticaldepth_ad(lay, j)
            predictors_ad%watervapour(7, lay, prof)  =      &
              & predictors_ad%watervapour(7, lay, prof) + fast_coef%watervapour(lay, chan, 7) * opticaldepth_ad(lay, j)
            predictors_ad%watervapour(6, lay, prof)  =      &
              & predictors_ad%watervapour(6, lay, prof) + fast_coef%watervapour(lay, chan, 6) * opticaldepth_ad(lay, j)
            predictors_ad%watervapour(5, lay, prof)  =      &
              & predictors_ad%watervapour(5, lay, prof) + fast_coef%watervapour(lay, chan, 5) * opticaldepth_ad(lay, j)
            predictors_ad%watervapour(4, lay, prof)  =      &
              & predictors_ad%watervapour(4, lay, prof) + fast_coef%watervapour(lay, chan, 4) * opticaldepth_ad(lay, j)
            predictors_ad%watervapour(3, lay, prof)  =      &
              & predictors_ad%watervapour(3, lay, prof) + fast_coef%watervapour(lay, chan, 3) * opticaldepth_ad(lay, j)
            predictors_ad%watervapour(2, lay, prof)  =      &
              & predictors_ad%watervapour(2, lay, prof) + fast_coef%watervapour(lay, chan, 2) * opticaldepth_ad(lay, j)
            predictors_ad%watervapour(1, lay, prof)  =      &
              & predictors_ad%watervapour(1, lay, prof) + fast_coef%watervapour(lay, chan, 1) * opticaldepth_ad(lay, j)
          ENDDO
        ELSE
          DO lay = nlayers, 1,  - 1
            IF (coef%nco > 0) THEN
              predictors_ad%co(1, lay, prof) =      &
                & predictors_ad%co(1, lay, prof) + fast_coef%watervapour(lay, chan, 15) * opticaldepth_ad(lay, j)
            ENDIF
            IF (coef%nco2 > 0) THEN
              predictors_ad%co2(1, lay, prof) =      &
                & predictors_ad%co2(1, lay, prof) + fast_coef%watervapour(lay, chan, 14) * opticaldepth_ad(lay, j)
            ENDIF
            predictors_ad%watervapour(13, lay, prof) =      &
              & predictors_ad%watervapour(13, lay, prof) + fast_coef%watervapour(lay, chan, 13) * opticaldepth_ad(lay, j)
            predictors_ad%watervapour(12, lay, prof) =      &
              & predictors_ad%watervapour(12, lay, prof) + fast_coef%watervapour(lay, chan, 12) * opticaldepth_ad(lay, j)
            predictors_ad%watervapour(11, lay, prof) =      &
              & predictors_ad%watervapour(11, lay, prof) + fast_coef%watervapour(lay, chan, 11) * opticaldepth_ad(lay, j)
            predictors_ad%watervapour(10, lay, prof) =      &
              & predictors_ad%watervapour(10, lay, prof) + fast_coef%watervapour(lay, chan, 10) * opticaldepth_ad(lay, j)
            predictors_ad%watervapour(9, lay, prof)  =      &
              & predictors_ad%watervapour(9, lay, prof) + fast_coef%watervapour(lay, chan, 9) * opticaldepth_ad(lay, j)
            predictors_ad%watervapour(8, lay, prof)  =      &
              & predictors_ad%watervapour(8, lay, prof) + fast_coef%watervapour(lay, chan, 8) * opticaldepth_ad(lay, j)
            predictors_ad%watervapour(7, lay, prof)  =      &
              & predictors_ad%watervapour(7, lay, prof) + fast_coef%watervapour(lay, chan, 7) * opticaldepth_ad(lay, j)
            predictors_ad%watervapour(6, lay, prof)  =      &
              & predictors_ad%watervapour(6, lay, prof) + fast_coef%watervapour(lay, chan, 6) * opticaldepth_ad(lay, j)
            predictors_ad%watervapour(5, lay, prof)  =      &
              & predictors_ad%watervapour(5, lay, prof) + fast_coef%watervapour(lay, chan, 5) * opticaldepth_ad(lay, j)
            predictors_ad%watervapour(4, lay, prof)  =      &
              & predictors_ad%watervapour(4, lay, prof) + fast_coef%watervapour(lay, chan, 4) * opticaldepth_ad(lay, j)
            predictors_ad%watervapour(3, lay, prof)  =      &
              & predictors_ad%watervapour(3, lay, prof) + fast_coef%watervapour(lay, chan, 3) * opticaldepth_ad(lay, j)
            predictors_ad%watervapour(2, lay, prof)  =      &
              & predictors_ad%watervapour(2, lay, prof) + fast_coef%watervapour(lay, chan, 2) * opticaldepth_ad(lay, j)
            predictors_ad%watervapour(1, lay, prof)  =      &
              & predictors_ad%watervapour(1, lay, prof) + fast_coef%watervapour(lay, chan, 1) * opticaldepth_ad(lay, j)
          ENDDO
        ENDIF
      ELSE IF (chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2295_25) THEN
        DO lay = nlayers, 1,  - 1
          IF (coef%nco > 0) THEN
            predictors_ad%co(1, lay, prof) =      &
              & predictors_ad%co(1, lay, prof) + fast_coef%watervapour(lay, chan, 15) * opticaldepth_ad(lay, j)
          ENDIF
          IF (coef%nco2 > 0) THEN
            predictors_ad%co2(1, lay, prof) =      &
              & predictors_ad%co2(1, lay, prof) + fast_coef%watervapour(lay, chan, 14) * opticaldepth_ad(lay, j)
          ENDIF
          predictors_ad%watervapour(13, lay, prof) =      &
            & predictors_ad%watervapour(13, lay, prof) + fast_coef%watervapour(lay, chan, 13) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(12, lay, prof) =      &
            & predictors_ad%watervapour(12, lay, prof) + fast_coef%watervapour(lay, chan, 12) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(11, lay, prof) =      &
            & predictors_ad%watervapour(11, lay, prof) + fast_coef%watervapour(lay, chan, 11) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(10, lay, prof) =      &
            & predictors_ad%watervapour(10, lay, prof) + fast_coef%watervapour(lay, chan, 10) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(9, lay, prof)  =      &
            & predictors_ad%watervapour(9, lay, prof) + fast_coef%watervapour(lay, chan, 9) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(8, lay, prof)  =      &
            & predictors_ad%watervapour(8, lay, prof) + fast_coef%watervapour(lay, chan, 8) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(7, lay, prof)  =      &
            & predictors_ad%watervapour(7, lay, prof) + fast_coef%watervapour(lay, chan, 7) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(6, lay, prof)  =      &
            & predictors_ad%watervapour(6, lay, prof) + fast_coef%watervapour(lay, chan, 6) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(5, lay, prof)  =      &
            & predictors_ad%watervapour(5, lay, prof) + fast_coef%watervapour(lay, chan, 5) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(4, lay, prof)  =      &
            & predictors_ad%watervapour(4, lay, prof) + fast_coef%watervapour(lay, chan, 4) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(3, lay, prof)  =      &
            & predictors_ad%watervapour(3, lay, prof) + fast_coef%watervapour(lay, chan, 3) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(2, lay, prof)  =      &
            & predictors_ad%watervapour(2, lay, prof) + fast_coef%watervapour(lay, chan, 2) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(1, lay, prof)  =      &
            & predictors_ad%watervapour(1, lay, prof) + fast_coef%watervapour(lay, chan, 1) * opticaldepth_ad(lay, j)
        ENDDO
      ELSE IF (chanx >= rttov9_wv2295_25 .AND. chanx < rttov9_wv2380_25) THEN
        DO lay = nlayers, 1,  - 1
          predictors_ad%watervapour(13, lay, prof) =      &
            & predictors_ad%watervapour(13, lay, prof) + fast_coef%watervapour(lay, chan, 13) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(12, lay, prof) =      &
            & predictors_ad%watervapour(12, lay, prof) + fast_coef%watervapour(lay, chan, 12) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(11, lay, prof) =      &
            & predictors_ad%watervapour(11, lay, prof) + fast_coef%watervapour(lay, chan, 11) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(10, lay, prof) =      &
            & predictors_ad%watervapour(10, lay, prof) + fast_coef%watervapour(lay, chan, 10) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(9, lay, prof)  =      &
            & predictors_ad%watervapour(9, lay, prof) + fast_coef%watervapour(lay, chan, 9) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(8, lay, prof)  =      &
            & predictors_ad%watervapour(8, lay, prof) + fast_coef%watervapour(lay, chan, 8) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(7, lay, prof)  =      &
            & predictors_ad%watervapour(7, lay, prof) + fast_coef%watervapour(lay, chan, 7) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(6, lay, prof)  =      &
            & predictors_ad%watervapour(6, lay, prof) + fast_coef%watervapour(lay, chan, 6) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(5, lay, prof)  =      &
            & predictors_ad%watervapour(5, lay, prof) + fast_coef%watervapour(lay, chan, 5) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(4, lay, prof)  =      &
            & predictors_ad%watervapour(4, lay, prof) + fast_coef%watervapour(lay, chan, 4) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(3, lay, prof)  =      &
            & predictors_ad%watervapour(3, lay, prof) + fast_coef%watervapour(lay, chan, 3) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(2, lay, prof)  =      &
            & predictors_ad%watervapour(2, lay, prof) + fast_coef%watervapour(lay, chan, 2) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(1, lay, prof)  =      &
            & predictors_ad%watervapour(1, lay, prof) + fast_coef%watervapour(lay, chan, 1) * opticaldepth_ad(lay, j)
        ENDDO
      ELSE   !chanx >= rttov9_wv2380_25
        DO lay = nlayers, 1,  - 1
          IF (coef%nch4 > 0) THEN
            predictors_ad%ch4(3, lay, prof) = predictors_ad%ch4(3, lay, prof) +      &
              & fast_coef%watervapour(lay, chan, 15) * predictors%ch4(1, lay, prof) ** 0.25 * opticaldepth_ad(lay, j)
            predictors_ad%ch4(1, lay, prof) = predictors_ad%ch4(1, lay, prof) +                       &
              & fast_coef%watervapour(lay, chan, 15) * 0.25 * predictors%ch4(1, lay, prof) ** ( - 0.75) *  &
              & predictors%ch4(3, lay, prof) * opticaldepth_ad(lay, j)
            predictors_ad%ch4(1, lay, prof) = predictors_ad%ch4(1, lay, prof) +      &
              & fast_coef%watervapour(lay, chan, 14) * 1.25 * predictors%ch4(1, lay, prof) ** 0.25 * opticaldepth_ad(lay, j)
          ENDIF
          predictors_ad%watervapour(19, lay, prof) =      &
            & predictors_ad%watervapour(19, lay, prof) + fast_coef%watervapour(lay, chan, 13) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(16, lay, prof) =      &
            & predictors_ad%watervapour(16, lay, prof) + fast_coef%watervapour(lay, chan, 12) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(11, lay, prof) =      &
            & predictors_ad%watervapour(11, lay, prof) + fast_coef%watervapour(lay, chan, 11) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(18, lay, prof) =      &
            & predictors_ad%watervapour(18, lay, prof) + fast_coef%watervapour(lay, chan, 10) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(15, lay, prof) =      &
            & predictors_ad%watervapour(15, lay, prof) + fast_coef%watervapour(lay, chan, 9) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(14, lay, prof) =      &
            & predictors_ad%watervapour(14, lay, prof) + fast_coef%watervapour(lay, chan, 8) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(7, lay, prof)  =      &
            & predictors_ad%watervapour(7, lay, prof) + fast_coef%watervapour(lay, chan, 7) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(6, lay, prof)  =      &
            & predictors_ad%watervapour(6, lay, prof) + fast_coef%watervapour(lay, chan, 6) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(5, lay, prof)  =      &
            & predictors_ad%watervapour(5, lay, prof) + fast_coef%watervapour(lay, chan, 5) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(4, lay, prof)  =      &
            & predictors_ad%watervapour(4, lay, prof) + fast_coef%watervapour(lay, chan, 4) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(3, lay, prof)  =      &
            & predictors_ad%watervapour(3, lay, prof) + fast_coef%watervapour(lay, chan, 3) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(2, lay, prof)  =      &
            & predictors_ad%watervapour(2, lay, prof) + fast_coef%watervapour(lay, chan, 2) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour(1, lay, prof)  =      &
            & predictors_ad%watervapour(1, lay, prof) + fast_coef%watervapour(lay, chan, 1) * opticaldepth_ad(lay, j)
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
          predictors_ad%mixedgas(8, lay, prof) =      &
            & predictors_ad%mixedgas(8, lay, prof) + fast_coef%mixedgas(lay, chan, 8) * opticaldepth_ad(lay, j)
          predictors_ad%mixedgas(7, lay, prof) =      &
            & predictors_ad%mixedgas(7, lay, prof) + fast_coef%mixedgas(lay, chan, 7) * opticaldepth_ad(lay, j)
          predictors_ad%mixedgas(6, lay, prof) =      &
            & predictors_ad%mixedgas(6, lay, prof) + fast_coef%mixedgas(lay, chan, 6) * opticaldepth_ad(lay, j)
          predictors_ad%mixedgas(5, lay, prof) =      &
            & predictors_ad%mixedgas(5, lay, prof) + fast_coef%mixedgas(lay, chan, 5) * opticaldepth_ad(lay, j)
          predictors_ad%mixedgas(4, lay, prof) =      &
            & predictors_ad%mixedgas(4, lay, prof) + fast_coef%mixedgas(lay, chan, 4) * opticaldepth_ad(lay, j)
          predictors_ad%mixedgas(3, lay, prof) =      &
            & predictors_ad%mixedgas(3, lay, prof) + fast_coef%mixedgas(lay, chan, 3) * opticaldepth_ad(lay, j)
          predictors_ad%mixedgas(2, lay, prof) =      &
            & predictors_ad%mixedgas(2, lay, prof) + fast_coef%mixedgas(lay, chan, 2) * opticaldepth_ad(lay, j)
          predictors_ad%mixedgas(1, lay, prof) =      &
            & predictors_ad%mixedgas(1, lay, prof) + fast_coef%mixedgas(lay, chan, 1) * opticaldepth_ad(lay, j)
        ENDDO
      ELSE
        DO lay = nlayers, 1,  - 1
          predictors_ad%mixedgas(10, lay, prof) =      &
            & predictors_ad%mixedgas(10, lay, prof) + fast_coef%mixedgas(lay, chan, 10) * opticaldepth_ad(lay, j)
          predictors_ad%mixedgas(9, lay, prof)  =      &
            & predictors_ad%mixedgas(9, lay, prof) + fast_coef%mixedgas(lay, chan, 9) * opticaldepth_ad(lay, j)
          predictors_ad%mixedgas(8, lay, prof)  =      &
            & predictors_ad%mixedgas(8, lay, prof) + fast_coef%mixedgas(lay, chan, 8) * opticaldepth_ad(lay, j)
          predictors_ad%mixedgas(7, lay, prof)  =      &
            & predictors_ad%mixedgas(7, lay, prof) + fast_coef%mixedgas(lay, chan, 7) * opticaldepth_ad(lay, j)
          predictors_ad%mixedgas(6, lay, prof)  =      &
            & predictors_ad%mixedgas(6, lay, prof) + fast_coef%mixedgas(lay, chan, 6) * opticaldepth_ad(lay, j)
          predictors_ad%mixedgas(5, lay, prof)  =      &
            & predictors_ad%mixedgas(5, lay, prof) + fast_coef%mixedgas(lay, chan, 5) * opticaldepth_ad(lay, j)
          predictors_ad%mixedgas(4, lay, prof)  =      &
            & predictors_ad%mixedgas(4, lay, prof) + fast_coef%mixedgas(lay, chan, 4) * opticaldepth_ad(lay, j)
          predictors_ad%mixedgas(3, lay, prof)  =      &
            & predictors_ad%mixedgas(3, lay, prof) + fast_coef%mixedgas(lay, chan, 3) * opticaldepth_ad(lay, j)
          predictors_ad%mixedgas(2, lay, prof)  =      &
            & predictors_ad%mixedgas(2, lay, prof) + fast_coef%mixedgas(lay, chan, 2) * opticaldepth_ad(lay, j)
          predictors_ad%mixedgas(1, lay, prof)  =      &
            & predictors_ad%mixedgas(1, lay, prof) + fast_coef%mixedgas(lay, chan, 1) * opticaldepth_ad(lay, j)
        ENDDO
      ENDIF
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_9_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_opdep_9_ad
