SUBROUTINE rttov_opdep_k( &
            & nlayers,      &
            & chanprof,     &
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
!  1.0
!  1.1    26/09/2003  Modified to allow for multiple polarisations (S English)
!  1.2    29/03/2005  Add end of header comment (J. Cameron)
!  1.3    01/06/2005  Marco Matricardi (ECMWF):
!            --       Computation ov variable OD_SFRAC added.
!  1.4    08/01/2007  Modified line in section 5 to reporduce rttov87 (RSaunders)
!  1.5    09/02/2007  Removed polarisation index.
!  1.6    15/09/2009  User defined ToA. Layers distinct from levels (P.Rayer)
!  1.7    03/11/2009  Transmittances / optical depths on levels (A Geer)
!  1.8    10/01/2013  Add PMC shifts (P Rayer)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Adjoint variables
! input transmission_k% tau_surf and transmission_k% tau_level set inside integrate_k
!
! input/output aux_k
!
! output predictors_k initialised inside rttov_k (need input
!    intent for memory allocation in calling routine)
!
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
  TYPE(rttov_path_pred), INTENT(INOUT) :: predictors_k
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  TYPE(rttov_fast_coef), INTENT(IN)    :: fast_coef
  TYPE(profile_aux    ), INTENT(IN)    :: aux
  TYPE(profile_aux    ), INTENT(INOUT) :: aux_k
  TYPE(opdp_path_Type ), INTENT(IN)    :: opdp_path                        ! optical depths
  TYPE(opdp_path_Type ), INTENT(INOUT) :: opdp_path_k
  REAL(KIND=jprb)      , INTENT(IN)    :: opdp_ref(nlayers, size(chanprof))
!INTF_END
!local variables:
  REAL(KIND=jprb) :: opticaldepth_k(nlayers, size(chanprof))
  REAL(KIND=jprb), POINTER :: debye_prof  (:, :)
  REAL(KIND=jprb), POINTER :: debye_prof_k(:, :)
  INTEGER(KIND=jpim) :: lev, lay, chan     , j, nlevels
  INTEGER(KIND=jpim) :: prof          , ii
! cloud liquid water local variables
  REAL   (KIND=jprb) :: zf, zf_sq         , z34_dif  , z45_dif     , z1_sq    , z2_sq    , z1_div   , z2_div
  REAL   (KIND=jprb) :: z1_den        , z2_den        , zastar   , z1_prod     , z2_prod  , z3_prod  , z4_prod
  REAL   (KIND=jprb) :: zbstar        , zbstar_sq     , za2star  , za2star_sq  , zdiv     , zgstar
  REAL   (KIND=jprb) :: z1f_sq_z1_sq  , z2f_sq_z2_sq
  REAL   (KIND=jprb) :: z34_dif_k     , z45_dif_k     , z1_sq_k  , z2_sq_k     , z1_div_k , z2_div_k
  REAL   (KIND=jprb) :: z1_den_k      , z2_den_k      , zastar_k , z1_prod_k   , z2_prod_k, z3_prod_k, z4_prod_k
  REAL   (KIND=jprb) :: zbstar_k      , zbstar_sq_k   , za2star_k, za2star_sq_k, zdiv_k   , zgstar_k
  REAL   (KIND=jprb) :: z1f_sq_z1_sq_k, z2f_sq_z2_sq_k, ztemp
  INTEGER(KIND=jpim) :: nchannels
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_K', 0_jpim, ZHOOK_HANDLE)
  nchannels = size(chanprof)
  nlevels   = nlayers + 1
!----------------------------------------
!2.Assemble layer optical depths
!----------------------------------------
! note that optical depth in the calculations is negative
  DO j = 1, nchannels
    ztemp = opdp_path_k%atm_level(nlevels, j)
    DO lev = nlevels, 2,  - 1
      lay = lev - 1
      opticaldepth_k(lay, j)        = ztemp
      ztemp = ztemp + opdp_path_k%atm_level(lev - 1, j)
      opdp_path_k%atm_level(lev, j) = 0.0_JPRB
    ENDDO
    opdp_path_k%atm_level(1, j) = 0.0_JPRB
  ENDDO
  WHERE (opdp_ref(:,:) > 0.0_JPRB)
    opticaldepth_k = 0.0_JPRB
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
          predictors_k%pmc(ii, lay, prof, chan) = &!predictors_k%pmc(ii, lay, prof, chan)                &
     &      + opticaldepth_k(lay, j) * coef%pmc_coef(lay, chan, ii)
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
        debye_prof_k => aux_k%debye_prof(:, :, j)
!cdir nodep
        DO lay = nlayers, 1,  - 1
          lev = lay + 1
          IF (lev >= coef%mwcldtop) THEN
! Repeat direct code
            zf = coef%frequency_ghz(chan)
            zf_sq = zf * zf
            z1_sq = debye_prof(1, lev) * debye_prof(1, lev)
            z2_sq = debye_prof(2, lev) * debye_prof(2, lev)
            z34_dif                  = debye_prof(3, lev) - debye_prof(4, lev)
            z45_dif                  = debye_prof(4, lev) - debye_prof(5, lev)
            z1f_sq_z1_sq             = zf_sq + z1_sq
            z2f_sq_z2_sq             = zf_sq + z2_sq
            z1_div                   = 1.0_JPRB / z1f_sq_z1_sq
            z2_div                   = 1.0_JPRB / z2f_sq_z2_sq
            z1_den                   = z34_dif * z1_div
            z2_den                   = z45_dif * z2_div
            zastar                   = debye_prof(3, lev) - zf_sq * (z1_den + z2_den)
            z1_prod                  = z34_dif * debye_prof(1, lev)
            z2_prod                  = z1_prod * z1_div
            z3_prod                  = z45_dif * debye_prof(2, lev)
            z4_prod                  = z3_prod * z2_div
            zbstar                   =  - zf * (z2_prod + z4_prod)
            zbstar_sq                = zbstar * zbstar
            za2star                  = zastar + 2.0_JPRB
            za2star_sq               = za2star * za2star
            zdiv = za2star_sq + zbstar_sq
            zgstar                   =  - 3.0_JPRB * zbstar / zdiv
! Now compute Adjoint code
!opticaldepth_k(lay,j)= opticaldepth_k(lay,j)
            zgstar_k                 = opticaldepth_k(lay, j) * ( - 1.5_JPRB * zf * predictors%clw(lay, prof))
            predictors_k%clw(lay, j) = &!predictors_k%clw(lay, j) + &
              opticaldepth_k(lay, j) * ( - 1.5_JPRB * zf * zgstar)
            zbstar_k                 =  - 3.0_JPRB * zgstar_k / zdiv
            zdiv_k                   = 3.0_JPRB * zgstar_k * zbstar / (zdiv * zdiv)
!zgstar_k =  0.
            za2star_sq_k             = zdiv_k
            zbstar_sq_k              = zdiv_k
!zdiv_k       = 0.
            za2star_k                = 2.0_JPRB * za2star * za2star_sq_k
!za2star_sq_k = 0.
            zastar_k                 = za2star_k
!za2star_k    = 0.
            zbstar_k                 = zbstar_k + 2.0_JPRB * zbstar * zbstar_sq_k
!zbstar_sq_k  = 0.
            z2_prod_k                =  - zf * zbstar_k
            z4_prod_k                =  - zf * zbstar_k
!zbstar_k     = 0.
            z3_prod_k                = z2_div * z4_prod_k
            z2_div_k                 = z3_prod * z4_prod_k
!z4_prod_k    = 0.
            z45_dif_k                = debye_prof(2, lev) * z3_prod_k
            debye_prof_k(2, lev)     = debye_prof_k(2, lev) + z45_dif * z3_prod_k
!z3_prod_k           = 0.
            z1_prod_k                = z1_div * z2_prod_k
            z1_div_k                 = z1_prod * z2_prod_k
!z2_prod_k    = 0.
            z34_dif_k                = debye_prof(1, lev) * z1_prod_k
            debye_prof_k(1, lev)     = debye_prof_k(1, lev) + z34_dif * z1_prod_k
!z1_prod_k           = 0.
            debye_prof_k(3, lev)     = debye_prof_k(3, lev) + zastar_k
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
            debye_prof_k(4, lev)     = debye_prof_k(4, lev) + z45_dif_k
            debye_prof_k(5, lev)     = debye_prof_k(5, lev) - z45_dif_k
!z45_dif_k           = 0.
            debye_prof_k(3, lev)     = debye_prof_k(3, lev) + z34_dif_k
            debye_prof_k(4, lev)     = debye_prof_k(4, lev) - z34_dif_k
!z34_dif_k           = 0.
            debye_prof_k(2, lev)     = debye_prof_k(2, lev) + z2_sq_k * 2.0_JPRB * debye_prof(2, lev)
!z2_sq_k             = 0.
            debye_prof_k(1, lev)     = debye_prof_k(1, lev) + z1_sq_k * 2.0_JPRB * debye_prof(1, lev)
!z1_sq_k             = 0.
          ELSE
            predictors_k%clw(lay, j) = 0._jprb ! initialisation moved to here from init_predictors
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
! !cdir nodep
!       DO lay = nlayers, 1,  - 1
!         predictors_k%ch4(:, lay, j) = predictors_k%ch4(:, lay, j) + fast_coef%ch4(lay, chan, :) * opticaldepth_k(lay, j)
!       ENDDO
!     ENDDO
!   ENDIF
! !-----------
! !1.7 add CO
! !-----------
!   IF (coef%nco > 0) THEN
!     DO j = 1, nchannels
!       chan = chanprof(j)%chan
! !cdir nodep
!       DO lay = nlayers, 1,  - 1
!         predictors_k%co(:, lay, j) = predictors_k%co(:, lay, j) + fast_coef%co(lay, chan, :) * opticaldepth_k(lay, j)
!       ENDDO
!     ENDDO
!   ENDIF
! !-----------
! !1.6 add N2O
! !-----------
!   IF (coef%nn2o > 0) THEN
!     DO j = 1, nchannels
!       chan = chanprof(j)%chan
! !cdir nodep
!       DO lay = nlayers, 1,  - 1
!         predictors_k%n2o(:, lay, j) = predictors_k%n2o(:, lay, j) + fast_coef%n2o(lay, chan, :) * opticaldepth_k(lay, j)
!       ENDDO
!     ENDDO
!   ENDIF
!-----------
!1.5 add CO2
!-----------
  IF (coef%nco2 > 0) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
!cdir nodep
      DO lay = nlayers, 1,  - 1
        predictors_k%co2(:, lay, j) = &!predictors_k%co2(:, lay, j) + 
          fast_coef%co2(lay, chan, :) * opticaldepth_k(lay, j)
      ENDDO
    ENDDO
  ENDIF
!------------------------------
!1.4 add Water Vapour Continuum
!------------------------------
  IF (coef%nwvcont > 0) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      DO lay = nlayers, 1,  - 1
        predictors_k%wvcont(:, lay, j) =      &
          ! predictors_k%wvcont(:, lay, j) + 
          fast_coef%wvcont(lay, chan, :) * opticaldepth_k(lay, j)
      ENDDO
    ENDDO
  ENDIF
!-------------
!1.3 add ozone
!-------------
  IF (coef%nozone > 0) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
        IF (coef%nozone == 11) THEN
!cdir nodep
        DO lay = 1, nlayers
          predictors_k%ozone(1, lay, j)  =      &
            ! predictors_k%ozone(1, lay, j) + 
            fast_coef%ozone(lay, chan, 1) * opticaldepth_k(lay, j)
          predictors_k%ozone(2, lay, j)  =      &
            ! predictors_k%ozone(2, lay, j) + 
            fast_coef%ozone(lay, chan, 2) * opticaldepth_k(lay, j)
          predictors_k%ozone(3, lay, j)  =      &
            ! predictors_k%ozone(3, lay, j) + 
            fast_coef%ozone(lay, chan, 3) * opticaldepth_k(lay, j)
          predictors_k%ozone(4, lay, j)  =      &
            ! predictors_k%ozone(4, lay, j) + 
            fast_coef%ozone(lay, chan, 4) * opticaldepth_k(lay, j)
          predictors_k%ozone(5, lay, j)  =      &
            ! predictors_k%ozone(5, lay, j) + 
            fast_coef%ozone(lay, chan, 5) * opticaldepth_k(lay, j)
          predictors_k%ozone(6, lay, j)  =      &
            ! predictors_k%ozone(6, lay, j) + 
            fast_coef%ozone(lay, chan, 6) * opticaldepth_k(lay, j)
          predictors_k%ozone(7, lay, j)  =      &
            ! predictors_k%ozone(7, lay, j) + 
            fast_coef%ozone(lay, chan, 7) * opticaldepth_k(lay, j)
          predictors_k%ozone(8, lay, j)  =      &
            ! predictors_k%ozone(8, lay, j) + 
            fast_coef%ozone(lay, chan, 8) * opticaldepth_k(lay, j)
          predictors_k%ozone(9, lay, j)  =      &
            ! predictors_k%ozone(9, lay, j) + 
            fast_coef%ozone(lay, chan, 9) * opticaldepth_k(lay, j)
          predictors_k%ozone(10, lay, j) =      &
            ! predictors_k%ozone(10, lay, j) + 
            fast_coef%ozone(lay, chan, 10) * opticaldepth_k(lay, j)
          predictors_k%ozone(11, lay, j) =      &
            ! predictors_k%ozone(11, lay, j) + 
            fast_coef%ozone(lay, chan, 11) * opticaldepth_k(lay, j)
        ENDDO
      ELSE
        DO ii = 1, coef%nozone
          DO lay = nlayers, 1,  - 1
            predictors_k%ozone(ii, lay, j) =      &
              ! predictors_k%ozone(ii, lay, j) + 
              fast_coef%ozone(lay, chan, ii) * opticaldepth_k(lay, j)
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
    IF (coef%nwater == 15) THEN
!cdir nodep
      DO lay = nlayers, 1,  - 1
        predictors_k%watervapour(1, lay, j)  =      &
!          & predictors_k%watervapour(1, lay, j) + &!
          fast_coef%watervapour(lay, chan, 1) * opticaldepth_k(lay, j)
        predictors_k%watervapour(2, lay, j)  =      &
!          & predictors_k%watervapour(2, lay, j) + &!
          fast_coef%watervapour(lay, chan, 2) * opticaldepth_k(lay, j)
        predictors_k%watervapour(3, lay, j)  =      &
!          & predictors_k%watervapour(3, lay, j) + &!
          fast_coef%watervapour(lay, chan, 3) * opticaldepth_k(lay, j)
        predictors_k%watervapour(4, lay, j)  =      &
!          & predictors_k%watervapour(4, lay, j) + &!
          fast_coef%watervapour(lay, chan, 4) * opticaldepth_k(lay, j)
        predictors_k%watervapour(5, lay, j)  =      &
!          & predictors_k%watervapour(5, lay, j) + &!
          fast_coef%watervapour(lay, chan, 5) * opticaldepth_k(lay, j)
        predictors_k%watervapour(6, lay, j)  =      &
!          & predictors_k%watervapour(6, lay, j) + &!
          fast_coef%watervapour(lay, chan, 6) * opticaldepth_k(lay, j)
        predictors_k%watervapour(7, lay, j)  =      &
!          & predictors_k%watervapour(7, lay, j) + &!
          fast_coef%watervapour(lay, chan, 7) * opticaldepth_k(lay, j)
        predictors_k%watervapour(8, lay, j)  =      &
!          & predictors_k%watervapour(8, lay, j) + &!
          fast_coef%watervapour(lay, chan, 8) * opticaldepth_k(lay, j)
        predictors_k%watervapour(9, lay, j)  =      &
!          & predictors_k%watervapour(9, lay, j) + &!
          fast_coef%watervapour(lay, chan, 9) * opticaldepth_k(lay, j)
        predictors_k%watervapour(10, lay, j) =      &
!          & predictors_k%watervapour(10, lay, j) + &!
          fast_coef%watervapour(lay, chan, 10) * opticaldepth_k(lay, j)
        predictors_k%watervapour(11, lay, j) =      &
!          & predictors_k%watervapour(11, lay, j) + &!
          fast_coef%watervapour(lay, chan, 11) * opticaldepth_k(lay, j)
        predictors_k%watervapour(12, lay, j) =      &
!          & predictors_k%watervapour(12, lay, j) + &!
          fast_coef%watervapour(lay, chan, 12) * opticaldepth_k(lay, j)
        predictors_k%watervapour(13, lay, j) =      &
!          & predictors_k%watervapour(13, lay, j) + &!
          fast_coef%watervapour(lay, chan, 13) * opticaldepth_k(lay, j)
        predictors_k%watervapour(14, lay, j) =      &
!          & predictors_k%watervapour(14, lay, j) + &!
          fast_coef%watervapour(lay, chan, 14) * opticaldepth_k(lay, j)
        predictors_k%watervapour(15, lay, j) =      &
!          & predictors_k%watervapour(15, lay, j) + &!
          fast_coef%watervapour(lay, chan, 15) * opticaldepth_k(lay, j)
      ENDDO
    ELSE
      DO ii = 1, coef%nwater
        DO lay = nlayers, 1,  - 1
          predictors_k%watervapour(ii, lay, j) =      &
          !  & predictors_k%watervapour(ii, lay, j) + 
          fast_coef%watervapour(lay, chan, ii) * opticaldepth_k(lay, j)
        ENDDO
      ENDDO
    ENDIF
  ENDDO
!--------------------------
!1.1 start with mixed gases
!--------------------------
  DO j = 1, nchannels
    chan = chanprof(j)%chan
    IF (coef%nmixed == 10) THEN
!cdir nodep
      DO lay = nlayers, 1,  - 1
        predictors_k%mixedgas(1, lay, j)  =      &
!          & predictors_k%mixedgas(1, lay, j) + &
          fast_coef%mixedgas(lay, chan, 1) * opticaldepth_k(lay, j)
        predictors_k%mixedgas(2, lay, j)  =      &
!          & predictors_k%mixedgas(2, lay, j) + &
          fast_coef%mixedgas(lay, chan, 2) * opticaldepth_k(lay, j)
        predictors_k%mixedgas(3, lay, j)  =      &
!          & predictors_k%mixedgas(3, lay, j) + &
          fast_coef%mixedgas(lay, chan, 3) * opticaldepth_k(lay, j)
        predictors_k%mixedgas(4, lay, j)  =      &
!          & predictors_k%mixedgas(4, lay, j) + &
          fast_coef%mixedgas(lay, chan, 4) * opticaldepth_k(lay, j)
        predictors_k%mixedgas(5, lay, j)  =      &
!          & predictors_k%mixedgas(5, lay, j) + &
          fast_coef%mixedgas(lay, chan, 5) * opticaldepth_k(lay, j)
        predictors_k%mixedgas(6, lay, j)  =      &
!          & predictors_k%mixedgas(6, lay, j) + &
          fast_coef%mixedgas(lay, chan, 6) * opticaldepth_k(lay, j)
        predictors_k%mixedgas(7, lay, j)  =      &
!          & predictors_k%mixedgas(7, lay, j) + &
          fast_coef%mixedgas(lay, chan, 7) * opticaldepth_k(lay, j)
        predictors_k%mixedgas(8, lay, j)  =      &
!          & predictors_k%mixedgas(8, lay, j) + &
          fast_coef%mixedgas(lay, chan, 8) * opticaldepth_k(lay, j)
        predictors_k%mixedgas(9, lay, j)  =      &
!          & predictors_k%mixedgas(9, lay, j) + &
          fast_coef%mixedgas(lay, chan, 9) * opticaldepth_k(lay, j)
        predictors_k%mixedgas(10, lay, j) =      &
!          & predictors_k%mixedgas(10, lay, j) + &
          fast_coef%mixedgas(lay, chan, 10) * opticaldepth_k(lay, j)
      ENDDO
    ELSE
      DO ii = 1, coef%nmixed
        DO lay = nlayers, 1,  - 1
          predictors_k%mixedgas(ii, lay, j) =      &
!            & predictors_k%mixedgas(ii, lay, j) + &
          fast_coef%mixedgas(lay, chan, ii) * opticaldepth_k(lay, j)
        ENDDO
      ENDDO
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_opdep_k
