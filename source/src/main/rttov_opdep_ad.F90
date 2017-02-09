SUBROUTINE rttov_opdep_ad( &
            & nlayers,       &
            & chanprof,      &
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
! Adjoint of rttov_transmit_tl
!  To calculate optical depths for a number of channels
!  and profiles from every pressure level to space.
! To interpolate optical depths on to levels of radiative transfer model
! (which, at present, entails only surface transmittance, as
! other optical depths are on *rt* levels) and to convert
! optical depths to transmittances.
!
! Code based on OPDEPAD and RTTAUAD from previous versions of RTTOV
! Only one profile per call
!
! input/output aux_ad
!
! output predictors_ad initialised inside rttov_ad (need input
!    intent for memory allocation in calling routine)
!
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
!  1.0    01/12/2002  New F90 code with structures (P Brunel A Smith)
!  1.1    29/01/2003  Add WV Continuum and CO2 capability (P Brunel)
!  1.2    04/12/2003  Optimisation (J Hague and D Salmond ECMWF)
!  1.3    26/09/2003  Modified to allow for multiple polarisations (S English)
!         06/09/2004  Mods. for Vectorisation (D Salmond ECMWF & B Carruthers, Cray)
!         28/02/2005  Improved vectorisation (D Dent)
!  1.4    01/06/2005  Marco Matricardi (ECMWF):
!            --       Computation ov variable OD_SFRAC added.
!  1.5    02/01/2007  Modified line in section 5 to reporduce rttov87
!  1.6    09/02/2007  Removed polarisation index (R Saunders)
!  1.7    15/08/2009  User defined ToA. Layers distinct from levels (P.Rayer)
!  1.8    03/11/2009  Transmittances / optical depths on levels (A Geer)
!  1.9    10/01/2013  Add PMC shifts (P Rayer)
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
       & rttov_chanprof,  &
       & rttov_coef,      &
       & rttov_fast_coef, &
       & predictors_Type, &
       & rttov_path_pred, &
       & opdp_path_Type,  &
       & profile_aux
  USE parkind1, ONLY : jpim, jprb
!INTF_OFF
  USE rttov_const, ONLY : sensor_id_mw
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  INTEGER(KIND=jpim)   , INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof ), INTENT(IN)    :: chanprof(:)
  TYPE(predictors_Type), INTENT(IN)    :: pred_info
  TYPE(rttov_path_pred), INTENT(IN)    :: predictors
  TYPE(rttov_path_pred), INTENT(INOUT) :: predictors_ad
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  TYPE(rttov_fast_coef), INTENT(IN)    :: fast_coef
  TYPE(profile_aux    ), INTENT(IN)    :: aux
  TYPE(profile_aux    ), INTENT(INOUT) :: aux_ad
  TYPE(opdp_path_Type ), INTENT(IN)    :: opdp_path                        ! optical depths
  TYPE(opdp_path_Type ), INTENT(INOUT) :: opdp_path_ad
  REAL(KIND=jprb)      , INTENT(IN)    :: opdp_ref(nlayers, size(chanprof))
!INTF_END
!local variables:
  REAL(KIND=jprb) :: opticaldepth_ad(nlayers, size(chanprof))
  REAL(KIND=jprb), POINTER :: debye_prof   (:, :)
  REAL(KIND=jprb), POINTER :: debye_prof_ad(:, :)
  INTEGER(KIND=jpim) :: lev, lay, chan      , j, nlevels
  INTEGER(KIND=jpim) :: prof
! cloud liquid water local variables
  REAL   (KIND=jprb) :: zf, zf_sq          , z34_dif   , z45_dif      , z1_sq     , z2_sq     , z1_div    , z2_div
  REAL   (KIND=jprb) :: z1_den         , z2_den         , zastar    , z1_prod      , z2_prod   , z3_prod   , z4_prod
  REAL   (KIND=jprb) :: zbstar         , zbstar_sq      , za2star   , za2star_sq   , zdiv      , zgstar
  REAL   (KIND=jprb) :: z1f_sq_z1_sq   , z2f_sq_z2_sq
  REAL   (KIND=jprb) :: z34_dif_ad     , z45_dif_ad     , z1_sq_ad  , z2_sq_ad     , z1_div_ad , z2_div_ad
  REAL   (KIND=jprb) :: z1_den_ad      , z2_den_ad      , zastar_ad , z1_prod_ad   , z2_prod_ad, z3_prod_ad, z4_prod_ad
  REAL   (KIND=jprb) :: zbstar_ad      , zbstar_sq_ad   , za2star_ad, za2star_sq_ad, zdiv_ad   , zgstar_ad
  REAL   (KIND=jprb) :: z1f_sq_z1_sq_ad, z2f_sq_z2_sq_ad, ztemp
  INTEGER(KIND=jpim) :: II
  INTEGER(KIND=jpim) :: nchannels
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_AD', 0_jpim, ZHOOK_HANDLE)
  nchannels = size(chanprof)
  nlevels   = nlayers + 1
!----------------------------------------
!2.Assemble layer optical depths
!----------------------------------------
! note that optical depth in the calculations is negative
  DO j = 1, nchannels
    ztemp = opdp_path_ad%atm_level(nlevels, j)
    DO lev = nlevels, 2,  - 1
      lay = lev - 1
      opticaldepth_ad(lay, j)        = ztemp
      ztemp = ztemp + opdp_path_ad%atm_level(lev - 1, j)
      opdp_path_ad%atm_level(lev, j) = 0.0_JPRB
    ENDDO
    opdp_path_ad%atm_level(1, j) = 0.0_JPRB
  ENDDO
  WHERE (opdp_ref(:,:) > 0.0_JPRB)
    opticaldepth_ad = 0.0_JPRB
  ENDWHERE

!------------------------------------------
!1.10 pressure-modulated cell (pmc) sensors
!------------------------------------------
  IF(coef%pmc_shift) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      prof = chanprof(j)%prof
      DO ii = 1, coef%pmc_nvar
        DO lay = 1, coef%pmc_nlay
          predictors_ad%pmc(ii, lay, prof, chan) = predictors_ad%pmc(ii, lay, prof, chan)                &
     &      + opticaldepth_ad(lay, j) * coef%pmc_coef(lay, chan, ii)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

!--------------------
!1.9 add liquid water (MW only)
!--------------------
  IF (coef%id_sensor == sensor_id_mw) THEN
    IF (pred_info%ncloud > 0) THEN
      DO j = 1, nchannels
        chan = chanprof(j)%chan
        prof = chanprof(j)%prof
        debye_prof => aux%debye_prof(:, :, prof)
        debye_prof_ad => aux_ad%debye_prof(:, :, prof)
!cdir nodep
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
            z1_div = 1.0_JPRB / z1f_sq_z1_sq
            z2_div = 1.0_JPRB / z2f_sq_z2_sq
            z1_den = z34_dif * z1_div
            z2_den = z45_dif * z2_div
            zastar = debye_prof(3, lev) - zf_sq * (z1_den + z2_den)
            z1_prod                      = z34_dif * debye_prof(1, lev)
            z2_prod                      = z1_prod * z1_div
            z3_prod                      = z45_dif * debye_prof(2, lev)
            z4_prod                      = z3_prod * z2_div
            zbstar =  - zf * (z2_prod + z4_prod)
            zbstar_sq                    = zbstar * zbstar
            za2star                      = zastar + 2.0_JPRB
            za2star_sq                   = za2star * za2star
            zdiv = za2star_sq + zbstar_sq
            zgstar =  - 3.0_JPRB * zbstar / zdiv
! Now compute Adjoint code
!opticaldepth_ad(lay,j)= opticaldepth_ad(lay,j)
            zgstar_ad                    = opticaldepth_ad(lay, j) * ( - 1.5_JPRB * zf * predictors%clw(lay, prof))
            predictors_ad%clw(lay, prof) =      &
              & predictors_ad%clw(lay, prof) + opticaldepth_ad(lay, j) * ( - 1.5_JPRB * zf * zgstar)
            zbstar_ad                    =  - 3.0_JPRB * zgstar_ad / zdiv
            zdiv_ad                      = 3.0_JPRB * zgstar_ad * zbstar / (zdiv * zdiv)
!zgstar_ad =  0.
            za2star_sq_ad                = zdiv_ad
            zbstar_sq_ad                 = zdiv_ad
!zdiv_ad       = 0.
            za2star_ad                   = 2.0_JPRB * za2star * za2star_sq_ad
!za2star_sq_ad = 0.
            zastar_ad                    = za2star_ad
!za2star_ad    = 0.
            zbstar_ad                    = zbstar_ad + 2.0_JPRB * zbstar * zbstar_sq_ad
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
            debye_prof_ad(2, lev)        = debye_prof_ad(2, lev) + z2_sq_ad * 2.0_JPRB * debye_prof(2, lev)
!z2_sq_ad             = 0.
            debye_prof_ad(1, lev)        = debye_prof_ad(1, lev) + z1_sq_ad * 2.0_JPRB * debye_prof(1, lev)
!z1_sq_ad             = 0.
          ENDIF
        ENDDO
      ENDDO
    ENDIF
  ENDIF
! N2O, CO and CH4 are not variable gases for v7 and v8 predictors.
! !-----------
! !1.8 add CH4
! !-----------
!   IF (coef%nch4 > 0) THEN
!     DO j = 1, nchannels
!       chan = chanprof(j)%chan
!       prof = chanprof(j)%prof
! !cdir nodep
!       DO lay = nlayers, 1,  - 1
!         predictors_ad%ch4(:, lay, prof) =      &
!           & predictors_ad%ch4(:, lay, prof) + fast_coef%ch4(lay, chan, :) * opticaldepth_ad(lay, j)
!       ENDDO
!     ENDDO
!   ENDIF
! !-----------
! !1.7 add CO
! !-----------
!   IF (coef%nco > 0) THEN
!     DO j = 1, nchannels
!       chan = chanprof(j)%chan
!       prof = chanprof(j)%prof
! !cdir nodep
!       DO lay = nlayers, 1,  - 1
!         predictors_ad%co(:, lay, prof) =      &
!           & predictors_ad%co(:, lay, prof) + fast_coef%co(lay, chan, :) * opticaldepth_ad(lay, j)
!       ENDDO
!     ENDDO
!   ENDIF
! !-----------
! !1.6 add N2O
! !-----------
!   IF (coef%nn2o > 0) THEN
!     DO j = 1, nchannels
!       chan = chanprof(j)%chan
!       prof = chanprof(j)%prof
! !cdir nodep
!       DO lay = nlayers, 1,  - 1
!         predictors_ad%n2o(:, lay, prof) =      &
!           & predictors_ad%n2o(:, lay, prof) + fast_coef%n2o(lay, chan, :) * opticaldepth_ad(lay, j)
!       ENDDO
!     ENDDO
!   ENDIF
!-----------
!1.5 add CO2
!-----------
  IF (coef%nco2 > 0) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      prof = chanprof(j)%prof
!cdir nodep
      DO lay = nlayers, 1,  - 1
        predictors_ad%co2(:, lay, prof) =      &
          & predictors_ad%co2(:, lay, prof) + fast_coef%co2(lay, chan, :) * opticaldepth_ad(lay, j)
      ENDDO
    ENDDO
  ENDIF
!------------------------------
!1.4 add Water Vapour Continuum
!------------------------------
  IF (coef%nwvcont > 0) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      prof = chanprof(j)%prof
      DO lay = nlayers, 1,  - 1
        predictors_ad%wvcont(:, lay, prof) =      &
          & predictors_ad%wvcont(:, lay, prof) + fast_coef%wvcont(lay, chan, :) * opticaldepth_ad(lay, j)
      ENDDO
    ENDDO
  ENDIF
!-------------
!1.3 add ozone
!-------------
  IF (coef%nozone > 0) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      prof = chanprof(j)%prof
      IF (coef%nozone == 11) THEN
!cdir nodep
        DO lay = nlayers, 1,  - 1
          predictors_ad%ozone(1, lay, prof)  =      &
            & predictors_ad%ozone(1, lay, prof) + fast_coef%ozone(lay, chan, 1) * opticaldepth_ad(lay, j)
          predictors_ad%ozone(2, lay, prof)  =      &
            & predictors_ad%ozone(2, lay, prof) + fast_coef%ozone(lay, chan, 2) * opticaldepth_ad(lay, j)
          predictors_ad%ozone(3, lay, prof)  =      &
            & predictors_ad%ozone(3, lay, prof) + fast_coef%ozone(lay, chan, 3) * opticaldepth_ad(lay, j)
          predictors_ad%ozone(4, lay, prof)  =      &
            & predictors_ad%ozone(4, lay, prof) + fast_coef%ozone(lay, chan, 4) * opticaldepth_ad(lay, j)
          predictors_ad%ozone(5, lay, prof)  =      &
            & predictors_ad%ozone(5, lay, prof) + fast_coef%ozone(lay, chan, 5) * opticaldepth_ad(lay, j)
          predictors_ad%ozone(6, lay, prof)  =      &
            & predictors_ad%ozone(6, lay, prof) + fast_coef%ozone(lay, chan, 6) * opticaldepth_ad(lay, j)
          predictors_ad%ozone(7, lay, prof)  =      &
            & predictors_ad%ozone(7, lay, prof) + fast_coef%ozone(lay, chan, 7) * opticaldepth_ad(lay, j)
          predictors_ad%ozone(8, lay, prof)  =      &
            & predictors_ad%ozone(8, lay, prof) + fast_coef%ozone(lay, chan, 8) * opticaldepth_ad(lay, j)
          predictors_ad%ozone(9, lay, prof)  =      &
            & predictors_ad%ozone(9, lay, prof) + fast_coef%ozone(lay, chan, 9) * opticaldepth_ad(lay, j)
          predictors_ad%ozone(10, lay, prof) =      &
            & predictors_ad%ozone(10, lay, prof) + fast_coef%ozone(lay, chan, 10) * opticaldepth_ad(lay, j)
          predictors_ad%ozone(11, lay, prof) =      &
            & predictors_ad%ozone(11, lay, prof) + fast_coef%ozone(lay, chan, 11) * opticaldepth_ad(lay, j)
        ENDDO
      ELSE
        DO ii = 1, coef%nozone
          DO lay = nlayers, 1,  - 1
            predictors_ad%ozone(ii, lay, prof) =      &
              & predictors_ad%ozone(ii, lay, prof) + fast_coef%ozone(lay, chan, ii) * opticaldepth_ad(lay, j)
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  ENDIF
!--------------------
!1.2 add water vapour
!--------------------
  DO j = 1, nchannels
    chan = chanprof(j)%chan
    prof = chanprof(j)%prof
    IF (coef%nwater == 15) THEN
!cdir nodep
      DO lay = nlayers, 1,  - 1
        predictors_ad%watervapour(1, lay, prof)  =      &
          & predictors_ad%watervapour(1, lay, prof) + fast_coef%watervapour(lay, chan, 1) * opticaldepth_ad(lay, j)
        predictors_ad%watervapour(2, lay, prof)  =      &
          & predictors_ad%watervapour(2, lay, prof) + fast_coef%watervapour(lay, chan, 2) * opticaldepth_ad(lay, j)
        predictors_ad%watervapour(3, lay, prof)  =      &
          & predictors_ad%watervapour(3, lay, prof) + fast_coef%watervapour(lay, chan, 3) * opticaldepth_ad(lay, j)
        predictors_ad%watervapour(4, lay, prof)  =      &
          & predictors_ad%watervapour(4, lay, prof) + fast_coef%watervapour(lay, chan, 4) * opticaldepth_ad(lay, j)
        predictors_ad%watervapour(5, lay, prof)  =      &
          & predictors_ad%watervapour(5, lay, prof) + fast_coef%watervapour(lay, chan, 5) * opticaldepth_ad(lay, j)
        predictors_ad%watervapour(6, lay, prof)  =      &
          & predictors_ad%watervapour(6, lay, prof) + fast_coef%watervapour(lay, chan, 6) * opticaldepth_ad(lay, j)
        predictors_ad%watervapour(7, lay, prof)  =      &
          & predictors_ad%watervapour(7, lay, prof) + fast_coef%watervapour(lay, chan, 7) * opticaldepth_ad(lay, j)
        predictors_ad%watervapour(8, lay, prof)  =      &
          & predictors_ad%watervapour(8, lay, prof) + fast_coef%watervapour(lay, chan, 8) * opticaldepth_ad(lay, j)
        predictors_ad%watervapour(9, lay, prof)  =      &
          & predictors_ad%watervapour(9, lay, prof) + fast_coef%watervapour(lay, chan, 9) * opticaldepth_ad(lay, j)
        predictors_ad%watervapour(10, lay, prof) =      &
          & predictors_ad%watervapour(10, lay, prof) + fast_coef%watervapour(lay, chan, 10) * opticaldepth_ad(lay, j)
        predictors_ad%watervapour(11, lay, prof) =      &
          & predictors_ad%watervapour(11, lay, prof) + fast_coef%watervapour(lay, chan, 11) * opticaldepth_ad(lay, j)
        predictors_ad%watervapour(12, lay, prof) =      &
          & predictors_ad%watervapour(12, lay, prof) + fast_coef%watervapour(lay, chan, 12) * opticaldepth_ad(lay, j)
        predictors_ad%watervapour(13, lay, prof) =      &
          & predictors_ad%watervapour(13, lay, prof) + fast_coef%watervapour(lay, chan, 13) * opticaldepth_ad(lay, j)
        predictors_ad%watervapour(14, lay, prof) =      &
          & predictors_ad%watervapour(14, lay, prof) + fast_coef%watervapour(lay, chan, 14) * opticaldepth_ad(lay, j)
        predictors_ad%watervapour(15, lay, prof) =      &
          & predictors_ad%watervapour(15, lay, prof) + fast_coef%watervapour(lay, chan, 15) * opticaldepth_ad(lay, j)
      ENDDO
    ELSE
      DO ii = 1, coef%nwater
        DO lay = nlayers, 1,  - 1
          predictors_ad%watervapour(ii, lay, prof) =      &
            & predictors_ad%watervapour(ii, lay, prof) + fast_coef%watervapour(lay, chan, ii) * opticaldepth_ad(lay, j)
        ENDDO
      ENDDO
    ENDIF
  ENDDO
!--------------------------
!1.1 start with mixed gases
!--------------------------
  DO j = 1, nchannels
    chan = chanprof(j)%chan
    prof = chanprof(j)%prof
    IF (coef%nmixed == 10) THEN
!cdir nodep
      DO lay = nlayers, 1,  - 1
        predictors_ad%mixedgas(1, lay, prof)  =      &
          & predictors_ad%mixedgas(1, lay, prof) + fast_coef%mixedgas(lay, chan, 1) * opticaldepth_ad(lay, j)
        predictors_ad%mixedgas(2, lay, prof)  =      &
          & predictors_ad%mixedgas(2, lay, prof) + fast_coef%mixedgas(lay, chan, 2) * opticaldepth_ad(lay, j)
        predictors_ad%mixedgas(3, lay, prof)  =      &
          & predictors_ad%mixedgas(3, lay, prof) + fast_coef%mixedgas(lay, chan, 3) * opticaldepth_ad(lay, j)
        predictors_ad%mixedgas(4, lay, prof)  =      &
          & predictors_ad%mixedgas(4, lay, prof) + fast_coef%mixedgas(lay, chan, 4) * opticaldepth_ad(lay, j)
        predictors_ad%mixedgas(5, lay, prof)  =      &
          & predictors_ad%mixedgas(5, lay, prof) + fast_coef%mixedgas(lay, chan, 5) * opticaldepth_ad(lay, j)
        predictors_ad%mixedgas(6, lay, prof)  =      &
          & predictors_ad%mixedgas(6, lay, prof) + fast_coef%mixedgas(lay, chan, 6) * opticaldepth_ad(lay, j)
        predictors_ad%mixedgas(7, lay, prof)  =      &
          & predictors_ad%mixedgas(7, lay, prof) + fast_coef%mixedgas(lay, chan, 7) * opticaldepth_ad(lay, j)
        predictors_ad%mixedgas(8, lay, prof)  =      &
          & predictors_ad%mixedgas(8, lay, prof) + fast_coef%mixedgas(lay, chan, 8) * opticaldepth_ad(lay, j)
        predictors_ad%mixedgas(9, lay, prof)  =      &
          & predictors_ad%mixedgas(9, lay, prof) + fast_coef%mixedgas(lay, chan, 9) * opticaldepth_ad(lay, j)
        predictors_ad%mixedgas(10, lay, prof) =      &
          & predictors_ad%mixedgas(10, lay, prof) + fast_coef%mixedgas(lay, chan, 10) * opticaldepth_ad(lay, j)
      ENDDO
    ELSE
      DO ii = 1, coef%nmixed
        DO lay = nlayers, 1,  - 1
          predictors_ad%mixedgas(ii, lay, prof) =      &
            & predictors_ad%mixedgas(ii, lay, prof) + fast_coef%mixedgas(lay, chan, ii) * opticaldepth_ad(lay, j)
        ENDDO
      ENDDO
    ENDIF
  ENDDO

  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_opdep_ad
