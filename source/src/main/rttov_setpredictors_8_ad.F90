!
SUBROUTINE rttov_setpredictors_8_ad( &
            & opts,          &
            & prof,          &
            & prof_ad,       &
            & geom,          &
            & coef,          &
            & aux,           &
            & predictors,    &
            & predictors_ad, &
            & raytracing,    &
            & raytracing_ad)
!
! Description
! RTTOV-8 Model
! AD of rttov_setpredictors_8
! To calculate and store the profile variables (predictors) required
! in subsequent transmittance calculations.
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
! see RTTOV7 science and validation report pages 18/19
! variable names are close to the documentation
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0   29/01/2003  Original - copy of RTTOV7 model (P Brunel)
!                               as a template for RTTOV-8
!  1.1   30/09/2003  Added predictors for wv line and continuum and CO2 (R Saunders)
!  1.2   03/06/2004  Parkind parametrisation, correction of w_ad, t_ad calculation
!                    simplify AD relatted to predictor 9 for WVL (P. Brunel)
!  1.3   23/02/2005  Correction of Twr definition (P. Brunel)
!  1.4   29/03/2005  Add end of header comment (J. Cameron)
!  1.5   07/12/2005  Add surface humidity (R. Saunders)
!  1.6   03/03/2006  Marco Matricardi (ECMWF):
!           --       Introduced altitude dependent local zenith angle
!  1.7   14/03/2007  Corrected CO2 profile logic (R Saunders)
!  1.8   15/08/2009  User defined ToA. Layers distinct from levels (P.Rayer)
!  1.9   02/12/2009  Pathsat, Pathsun and related quantities are now layer arrays
!                    (Marco Matricardi).
!  1.10  10/11/2010  Remove rttov9_compat flag from code (J Hocking)
!  1.11  25/10/2011  Fix bug so that omitted trace gas profiles are treated 
!                    correctly (J Hocking)
!  1.12  12/10/2011  Added fix to stop ifort optimising array initialisation causing 
!                    problems (D Rundle)
!  1.13  23/11/2011  Fix bug in CO2 predictors (J Hocking)
!  1.14  10/01/2013  Add PMC shifts (P Rayer)
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
       & rttov_coef,      &
       & rttov_options,   &
       & profile_Type,    &
       & geometry_Type,   &
       & profile_aux,     &
       & rttov_path_pred, &
       & raytracing_type
!INTF_OFF
  USE rttov_const, ONLY : gravity, sensor_id_mw
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options  ), INTENT(IN)    :: opts
  TYPE(profile_Type   ), INTENT(IN)    :: prof(:)            ! profile (ppmv dry)
  TYPE(profile_Type   ), INTENT(INOUT) :: prof_ad(size(prof))
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  TYPE(geometry_Type  ), INTENT(IN)    :: geom(size(prof))
  TYPE(rttov_path_pred), INTENT(IN)    :: predictors
  TYPE(rttov_path_pred), INTENT(INOUT) :: predictors_ad
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing_ad
  TYPE(profile_aux    ), INTENT(IN)    :: aux                ! auxillary profiles info.
!INTF_END
!local variables:
  INTEGER(KIND=jpim) :: level , layer , iprof
  INTEGER(KIND=jpim) :: iv2lev, iv3lev, iv2lay
! user profile
  REAL   (KIND=Jprb) :: t(prof(1)%nlayers)
  REAL   (KIND=Jprb) :: w(prof(1)%nlayers)
  REAL   (KIND=Jprb) :: o(prof(1)%nlayers)
  REAL   (KIND=Jprb) :: co2  (prof(1)%nlayers)
! reference profile
  REAL   (KIND=Jprb) :: tr   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: wr   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: or   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: co2r (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: wwr  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: twr  (prof(1)%nlayers)
! user - reference
  REAL   (KIND=Jprb) :: dt   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: dtabs(prof(1)%nlayers)
! pressure weighted
  REAL   (KIND=Jprb) :: tw   (prof(1)%nlayers)
  REAL   (KIND=jprb) :: ww   (prof(1)%nlayers)
  REAL   (KIND=jprb) :: ow   (prof(1)%nlayers)
  REAL   (KIND=jprb) :: co2w (prof(1)%nlayers)
! intermediate variables
  REAL   (KIND=Jprb) :: sum1, sum2
  REAL   (KIND=Jprb) :: deltac     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sum2_ww    (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sum2_wwr   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sum2_ow    (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sum2_twr   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sum2_co2w  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tr_sq      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tr_sqrt    (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tr_4       (prof(1)%nlayers)
  REAL   (KIND=jprb) :: sec_wrwr   (prof(1)%nlayers)
  REAL   (KIND=jprb) :: sec_wr     (prof(1)%nlayers)
! AD variables
  REAL   (KIND=Jprb) :: t_ad       (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: w_ad       (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: o_ad       (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: co2_ad     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tr_ad      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: wr_ad      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: or_ad      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: wwr_ad     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: co2r_ad    (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: twr_ad     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: dt_ad      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: dto_ad     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tw_ad      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: ww_ad      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: ow_ad      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: co2w_ad    (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sec_or_ad  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sec_wr_ad  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sec_wrwr_ad(prof(1)%nlayers)
  INTEGER(KIND=jpim) :: nprofiles, nlayers

! pressure-moodulated cell (pmc) variables
  REAL   (KIND=Jprb) :: Lcel_cm
  REAL   (KIND=Jprb) :: Tcel
  REAL   (KIND=Jprb) :: betaplus1
  REAL   (KIND=Jprb) :: acm
  REAL   (KIND=Jprb) :: Pcel_Lev
  REAL   (KIND=Jprb) :: Pnom_LevM1
  REAL   (KIND=Jprb) :: Pcel_LevM1
  REAL   (KIND=Jprb) :: Pnom_Lev
  REAL   (KIND=Jprb) :: Pnom(coef%fmv_chn)
  REAL   (KIND=Jprb) :: Pcel(coef%fmv_chn)
  INTEGER(KIND=jpim) :: ichan
  INTEGER(KIND=jpim) :: nlay
  REAL   (KIND=Jprb) :: acm_AD
  REAL   (KIND=Jprb) :: Pcel_Lev_AD
  REAL   (KIND=Jprb) :: Pnom_LevM1_AD
  REAL   (KIND=Jprb) :: Pcel_LevM1_AD
  REAL   (KIND=Jprb) :: Pnom_Lev_AD

  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
!-------------------------------------------------------------------------------
! Recompute Direct variables
!-------------------------------------------------------------------------------
!1) Profile layer quantities
!-------------------------------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_8_AD', 0_jpim, ZHOOK_HANDLE)
  nprofiles = size(prof)
  nlayers = prof(1)%nlayers

  DO iprof = 1, nprofiles

    DO layer = 1, prof(1)%nlayers
      level    = layer + 1
      t(layer) = (prof(iprof)%t(level - 1) + prof(iprof)%t(level)) / 2._JPRB
      w(layer) = (prof(iprof)%q(level - 1) + prof(iprof)%q(level)) / 2._JPRB
    ENDDO

!

    IF (opts%rt_all%use_q2m) THEN
! include surface humidity
      iv3lev = aux%s(iprof)%nearestlev_surf - 1
      iv2lev = aux%s(iprof)%nearestlev_surf

      IF (iv2lev <= coef%nlevels) THEN
        iv2lay    = iv2lev - 1
        w(iv2lay) = (prof(iprof)%s2m%q + prof(iprof)%q(iv3lev)) / 2._JPRB
      ENDIF

    ENDIF

!

    IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN

      DO layer = 1, prof(1)%nlayers
        level    = layer + 1
        o(layer) = (prof(iprof)%o3(level - 1) + prof(iprof)%o3(level)) / 2._JPRB
      ENDDO

    ENDIF


    IF (opts%rt_ir%co2_Data .AND. coef%nco2 > 0) THEN

      DO layer = 1, prof(1)%nlayers
        level      = layer + 1
        co2(layer) = (prof(iprof)%co2(level - 1) + prof(iprof)%co2(level)) / 2._JPRB
      ENDDO

    ENDIF


!------------------------------------------------------------------------------
!2) calculate deviations from reference profile (layers)
!------------------------------------------------------------------------------
! All assignments 1:prof(1) % nlayers
    dt(:)      = t(:) - coef%tstar(:)
    dtabs(:)   = Abs(dt(:))
!------------------------------------------------------------------------------
!3) calculate (profile / reference profile) ratios; tr wr or
! if no input O3 profile, set to reference value (or =1)
!------------------------------------------------------------------------------
! All assignments 1:prof(1) % nlayers
    tr(:)      = t(:) / coef%tstar(:)
    tr_sq(:)   = tr(:) * tr(:)
    tr_4(:)    = tr_sq(:) * tr_sq(:)
    tr_sqrt(:) = Sqrt(tr(:))
    wr(:)      = w(:) / coef%wstar(:)

    IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN
      or(:) = o(:) / coef%ostar(:)
    ELSE
      or(:) = 1._JPRB
    ENDIF

    IF (coef%nco2 > 0 .AND. opts%rt_ir%co2_Data) THEN
      co2r(:) = co2(:) / coef%co2star(:)
    ELSE
      co2r(:) = 1._JPRB
    ENDIF
    
!-------------------------------------------------------------------
! 4. calculate profile / reference profile sums: tw wwr
!--------------------------------------------------------------------
    tw(1) = 0.

    DO layer = 2, prof(1)%nlayers
      tw(layer) = tw(layer - 1) + coef%dpp(layer - 1) * tr(layer - 1)
    ENDDO

    sum1 = 0._JPRB
    sum2 = 0._JPRB

    DO layer = 1, prof(1)%nlayers
      sum1            = sum1 + coef%dpp(layer - 1) * w(layer) * t(layer)
      sum2            = sum2 + coef%dpp(layer - 1) * coef%wstar(layer) * coef%tstar(layer)
      sum2_wwr(layer) = sum2
      wwr(layer)      = sum1 / sum2
    ENDDO

    sum1 = 0._JPRB
    sum2 = 0._JPRB

    DO layer = 1, prof(1)%nlayers
      sum1           = sum1 + coef%dpp(layer - 1) * w(layer)
      sum2           = sum2 + coef%dpp(layer - 1) * coef%wstar(layer)
      sum2_ww(layer) = sum2
      ww(layer)      = sum1 / sum2
    ENDDO


    IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN
      sum1 = 0._JPRB
      sum2 = 0._JPRB

      DO layer = 1, prof(1)%nlayers
        sum1           = sum1 + coef%dpp(layer - 1) * o(layer)
        sum2           = sum2 + coef%dpp(layer - 1) * coef%ostar(layer)
        sum2_ow(layer) = sum2
        ow(layer)      = sum1 / sum2
      ENDDO

    ELSE
      sum2_ow(:) = 0._JPRB
      ow(:)      = 1._JPRB
    ENDIF


    IF (opts%rt_ir%co2_Data .AND. coef%nco2 > 0) THEN
      sum1 = 0._JPRB
      sum2 = 0._JPRB

      DO layer = 1, prof(1)%nlayers
        sum1             = sum1 + coef%dpp(layer - 1) * co2(layer)
        sum2             = sum2 + coef%dpp(layer - 1) * coef%co2star(layer)
        sum2_co2w(layer) = sum2
        co2w(layer)      = sum1 / sum2
      ENDDO
    
    ELSE
      sum2_co2w(:) = 0._JPRB
      co2w(:)      = 1._JPRB
    ENDIF
    
    IF (coef%nco2 > 0) THEN
      sum1        = 0._JPRB
      sum2        = 0._JPRB
      sum2_twr(1) = 0._JPRB
      twr(1)      = 1._JPRB
      
      DO layer = 2, prof(1)%nlayers
        sum1            = sum1 + coef%dpp(layer - 1) * t(layer - 1)
        sum2            = sum2 + coef%dpp(layer - 1) * coef%tstar(layer - 1)
        sum2_twr(layer) = sum2
        twr(layer)      = sum1 / sum2
      ENDDO

    ENDIF
  

!-------------------------------------------------------------------------
! Adjoint code
!-------------------------------------------------------------------------
! DAR - 1:nlayers is seemingly necessary to stop ifort from optimising these lines out at -O3
    w_ad(1:nlayers)        = 0._JPRB
    wr_ad(1:nlayers)       = 0._JPRB
    ww_ad(1:nlayers)       = 0._JPRB
    wwr_ad(1:nlayers)      = 0._JPRB
    sec_wr_ad(1:nlayers)   = 0._JPRB
    sec_wrwr_ad(1:nlayers) = 0._JPRB
    dt_ad(1:nlayers)       = 0._JPRB
    t_ad(1:nlayers)        = 0._JPRB
    tr_ad(1:nlayers)       = 0._JPRB
    tw_ad(1:nlayers)       = 0._JPRB

!5.7  pressure-modulated cell (pmc) changes
!-------------------------------------------

      IF (COEF%PMC_SHIFT) THEN

! FWD
! Constants
        Lcel_cm=coef%pmc_lengthcell
        Tcel=coef%pmc_tempcell
        betaplus1=coef%pmc_betaplus1
        DO ichan=1,coef%fmv_chn
          Pnom(ichan)=coef%pmc_pnominal(ichan)
          Pcel(ichan)=coef%pmc_ppmc(ichan)
        ENDDO
        nlay=coef%pmc_nlay

! This is on coef levels so pressure is not an active variable in the AD

      DO layer = 1, prof(1)%nlayers
        level=layer+1
! nlay (for pmc) may be less than prof(1)%nlayers
        IF(LAYER .LE. NLAY) THEN
          acm = raytracing%co2_cm(iprof)/(2._JPRB*betaplus1*Lcel_cm)
          DO ichan=1, coef%fmv_chn
            Pcel_Lev= acm*raytracing%pathsat(layer,iprof)*prof(iprof)%p(level  )**2_JPIM + Pcel(ichan)**2_JPIM
            Pnom_LevM1=acm*raytracing%pathsat(layer,iprof)*prof(iprof)%p(level-1)**2_JPIM + Pnom(ichan)**2_JPIM
            Pcel_LevM1= acm*raytracing%pathsat(layer,iprof)*prof(iprof)%p(level-1)**2_JPIM + Pcel(ichan)**2_JPIM
            Pnom_Lev=acm*raytracing%pathsat(layer,iprof)*prof(iprof)%p(level  )**2_JPIM + Pnom(ichan)**2_JPIM
! AD
            Pcel_Lev_AD =  0._JPRB
            Pnom_LevM1_AD = 0._JPRB
            Pcel_LevM1_AD =  0._JPRB
            Pnom_Lev_AD = 0._JPRB
            acm_AD = 0._JPRB
! for PMC predictor 2
            predictors_AD%pmc(2,layer,iprof,ichan)=0._JPRB
! for PMC predictor 1
            Pcel_Lev_AD = Pcel_Lev_AD                                                    &
         &   + predictors_AD%pmc(1,layer,iprof,ichan) * Pnom_LevM1/(Pcel_Lev*Pnom_LevM1)
            Pnom_LevM1_AD = Pnom_LevM1_AD                                                  &
         &   + predictors_AD%pmc(1,layer,iprof,ichan) * Pcel_Lev/(Pcel_Lev*Pnom_LevM1)
            Pcel_LevM1_AD = Pcel_LevM1_AD                                                    &
         &   - predictors_AD%pmc(1,layer,iprof,ichan) * Pnom_Lev/(Pcel_LevM1*Pnom_Lev)
            Pnom_Lev_AD = Pnom_Lev_AD                                                  &
         &   - predictors_AD%pmc(1,layer,iprof,ichan) * Pcel_LevM1/(Pcel_LevM1*Pnom_Lev)
            predictors_AD%pmc(1,layer,iprof,ichan) = 0._JPRB

            acm_AD = acm_AD                                                              &
         &   + Pnom_Lev_AD * raytracing%pathsat(layer,iprof)*prof(iprof)%p(level  )**2_JPIM
            raytracing_AD%pathsat(layer,iprof) = raytracing_AD%pathsat(layer,iprof)      &
         &   + Pnom_Lev_AD * acm * prof(iprof)%p(level  )**2_JPIM
            prof_AD(iprof)%p(level  ) =  0._JPRB
            Pnom_Lev_AD = 0._JPRB

            acm_AD = acm_AD                                                              &
         &   + Pcel_LevM1_AD * raytracing%pathsat(layer,iprof)*prof(iprof)%p(level-1)**2_JPIM
            raytracing_AD%pathsat(layer,iprof) = raytracing_AD%pathsat(layer,iprof)      &
         &   + Pcel_LevM1_AD * acm * prof(iprof)%p(level-1)**2_JPIM
            prof_AD(iprof)%p(level-1) =  0._JPRB
            Pcel_LevM1_AD = 0._JPRB

            acm_AD = acm_AD                                                              &
         &   + Pnom_LevM1_AD * raytracing%pathsat(layer,iprof)*prof(iprof)%p(level-1)**2_JPIM
            raytracing_AD%pathsat(layer,iprof) = raytracing_AD%pathsat(layer,iprof)      &
         &   + Pnom_LevM1_AD * acm * prof(iprof)%p(level-1)**2_JPIM
            prof_AD(iprof)%p(level-1) =  0._JPRB
            Pnom_LevM1_AD = 0._JPRB

            acm_AD = acm_AD                                                              &
         &   + Pcel_Lev_AD * raytracing%pathsat(layer,iprof)*prof(iprof)%p(level  )**2_JPIM
            raytracing_AD%pathsat(layer,iprof) = raytracing_AD%pathsat(layer,iprof)      &
         &   + Pcel_Lev_AD * acm * prof(iprof)%p(level  )**2_JPIM
            prof_AD(iprof)%p(level-1) =  0._JPRB
            Pcel_Lev_AD = 0._JPRB

            raytracing_AD%co2_cm(iprof) = raytracing_AD%co2_cm(iprof)                    &
         &   + acm_AD/(2._JPRB*betaplus1*Lcel_cm)
            acm_AD = 0._JPRB

          ENDDO
        ENDIF
      ENDDO

      ENDIF
    

!5.6 CO2
!-------

    IF (coef%nco2 > 0) THEN
      co2r_ad(:) = 0._JPRB
      co2w_ad(:) = 0._JPRB
      twr_ad(:)  = 0._JPRB
      co2_ad(:)  = 0._JPRB
!

      DO layer = 1, prof(1)%nlayers
        level = layer + 1
        twr_ad(layer)                      =      &
          & twr_ad(layer) + predictors_ad%co2(10, layer, iprof) * raytracing%pathsat(layer, iprof) * tr_sqrt(layer)
        tr_ad(layer)                       =      &
          & tr_ad(layer) + predictors_ad%co2(10, layer, iprof) * 0.5_JPRB * predictors%co2(7, layer, iprof) / tr_sqrt(layer)
        raytracing_ad%pathsat(layer, iprof) =     &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%co2(10, layer, iprof) * twr(layer) * tr_sqrt(layer)
        twr_ad(layer)                      = twr_ad(layer) + &
          & predictors_ad%co2(9, layer, iprof) * 3._JPRB * raytracing%pathsat(layer, iprof) * &
          & predictors%co2(9, layer, iprof) / predictors%co2(7, layer, iprof)
        co2w_ad(layer)                     = co2w_ad(layer) + &
          & 2._JPRB * raytracing%pathsat(layer, iprof) * predictors_ad%co2(8, layer, iprof) * &
          & Sqrt(predictors%co2(8, layer, iprof))
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) + &
          & 2._JPRB * co2w(layer) * predictors_ad%co2(8, layer, iprof) * &
          & Sqrt(predictors%co2(8, layer, iprof))
        twr_ad(layer)                      =      &
          & twr_ad(layer) + predictors_ad%co2(7, layer, iprof) * raytracing%pathsat(layer, iprof)
        raytracing_ad%pathsat(layer, iprof) =     &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%co2(7, layer, iprof) * twr(layer)
        raytracing_ad%pathsat(layer, iprof) =     &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%co2(6, layer, iprof)
        tr_ad(layer)                       = tr_ad(layer) + predictors_ad%co2(5, layer, iprof)
        tr_ad(layer)                       =      &
          & tr_ad(layer) + 2._JPRB * predictors_ad%co2(4, layer, iprof) * predictors%co2(3, layer, iprof)
        raytracing_ad%pathsat(layer, iprof) =     &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%co2(4, layer, iprof) * predictors%co2(2, layer, iprof)
        tr_ad(layer)                       =      &
          & tr_ad(layer) + predictors_ad%co2(3, layer, iprof) * raytracing%pathsat(layer, iprof)
        raytracing_ad%pathsat(layer, iprof) =     &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%co2(3, layer, iprof) * tr(layer)
        tr_ad(layer)                       =      &
          & tr_ad(layer) + 2._JPRB * predictors_ad%co2(2, layer, iprof) * predictors%co2(5, layer, iprof)
        co2r_ad(layer)                     =      &
          & co2r_ad(layer) + predictors_ad%co2(1, layer, iprof) * raytracing%pathsat(layer, iprof)
        raytracing_ad%pathsat(layer, iprof) =     &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%co2(1, layer, iprof) * co2r(layer)
      ENDDO

    ENDIF

!5.5 cloud
!---------

    IF (coef%id_sensor == sensor_id_mw) THEN
      IF (opts%rt_mw%clw_Data) THEN

        DO layer = 1, prof(1)%nlayers
          deltac(layer) = 0.1820_JPRB * 100.0_JPRB * coef%dp(layer) / (4.3429_JPRB * gravity)
        ENDDO


        DO layer = 2, prof(1)%nlayers
          level = layer + 1
          prof_ad(iprof)%clw(level - 1)   =      &
            & prof_ad(iprof)%clw(level - 1) + 0.5_JPRB * predictors_ad%clw(layer, iprof) * deltac(layer) * geom(iprof)%seczen
          predictors_ad%clw(layer, iprof) = 0.5_JPRB * predictors_ad%clw(layer, iprof)
        ENDDO


        DO layer = 1, prof(1)%nlayers
          level = layer + 1
          prof_ad(iprof)%clw(level) =      &
            & prof_ad(iprof)%clw(level) + predictors_ad%clw(layer, iprof) * deltac(layer) * geom(iprof)%seczen
        ENDDO

      ENDIF
    ENDIF

!5.4 ozone
!---------

    IF (coef%nozone > 0) THEN
      o_ad(:)      = 0._JPRB
      or_ad(:)     = 0._JPRB
      wr_ad(:)     = 0._JPRB
      ow_ad(:)     = 0._JPRB
      dto_ad(:)    = 0._JPRB
      sec_or_ad(:) = 0._JPRB

      DO layer = 1, prof(1)%nlayers
        level = layer + 1
! One can pack all ow_ad lines in one longer statement
! same for sec_or_ad and dto_ad
        ow_ad(layer)                        = ow_ad(layer) +                                &
          & predictors_ad%ozone(11, layer, iprof) * 2._JPRB * raytracing%pathsat(layer, iprof) *  &
          & predictors%ozone(10, layer, iprof)
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%ozone(11, layer, iprof) * 2._JPRB * ow(layer) * predictors%ozone(10, layer, iprof)
        ow_ad(layer)                        =      &
          & ow_ad(layer) + predictors_ad%ozone(10, layer, iprof) * raytracing%pathsat(layer, iprof)
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%ozone(10, layer, iprof) * ow(layer)
        or_ad(layer)                        = or_ad(layer) +                                             &
          & predictors_ad%ozone(9, layer, iprof) * SQRT(raytracing%pathsat(layer, iprof) * ow(layer)) *  &
          & raytracing%pathsat(layer, iprof)
        ow_ad(layer)                        = ow_ad(layer) +                                  &
          & predictors_ad%ozone(9, layer, iprof) * predictors%ozone(1, layer, iprof) * 0.5_JPRB *  &
          & raytracing%pathsat(layer, iprof) / SQRT(raytracing%pathsat(layer, iprof) * ow(layer))
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%ozone(9, layer, iprof) * 1.5_JPRB * or(layer) * SQRT(predictors%ozone(10, layer, iprof))
        or_ad(layer)                        =      &
          & or_ad(layer) + predictors_ad%ozone(8, layer, iprof) * predictors%ozone(10, layer, iprof)
        ow_ad(layer)                        =      &
          & ow_ad(layer) + predictors_ad%ozone(8, layer, iprof) * raytracing%pathsat(layer, iprof) * or(layer)
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%ozone(8, layer, iprof) * or(layer) * ow(layer)
        or_ad(layer)                        =      &
          & or_ad(layer) + predictors_ad%ozone(7, layer, iprof) * 1.5_JPRB * predictors%ozone(2, layer, iprof) / ow(layer)
        ow_ad(layer)                        =      &
          & ow_ad(layer) - predictors_ad%ozone(7, layer, iprof) * predictors%ozone(7, layer, iprof) / ow(layer)
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%ozone(7, layer, iprof) * 0.5_JPRB * or(layer) ** 1.5_JPRB /              &
          & (raytracing%pathsat(layer, iprof) ** 0.5_JPRB * ow(layer))
        or_ad(layer)                        =      &
          & or_ad(layer) + predictors_ad%ozone(6, layer, iprof) * 2._JPRB * predictors%ozone(1, layer, iprof) * ow(layer)
        ow_ad(layer)                        = ow_ad(layer) +      &
          & predictors_ad%ozone(6, layer, iprof) * predictors%ozone(4, layer, iprof) / raytracing%pathsat(layer, iprof)
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%ozone(6, layer, iprof) * or(layer) ** 2 * ow(layer)
        sec_or_ad(layer)                    = sec_or_ad(layer) +                                   &
          & predictors_ad%ozone(5, layer, iprof) * 0.5_JPRB * predictors%ozone(3, layer, iprof) /  &
          & (predictors%ozone(1, layer, iprof) * predictors%ozone(2, layer, iprof))
        dto_ad(layer)                       =      &
          & dto_ad(layer) + predictors_ad%ozone(5, layer, iprof) * predictors%ozone(2, layer, iprof)
        sec_or_ad(layer)                    =      &
          & sec_or_ad(layer) + predictors_ad%ozone(4, layer, iprof) * 2._JPRB * predictors%ozone(1, layer, iprof)
        sec_or_ad(layer)                    = sec_or_ad(layer) +      &
          & predictors_ad%ozone(3, layer, iprof) * predictors%ozone(3, layer, iprof) / predictors%ozone(1, layer, iprof)
        dto_ad(layer)                       =      &
          & dto_ad(layer) + predictors_ad%ozone(3, layer, iprof) * predictors%ozone(1, layer, iprof)
        sec_or_ad(layer)                    =      &
          & sec_or_ad(layer) + predictors_ad%ozone(2, layer, iprof) * 0.5_JPRB / predictors%ozone(2, layer, iprof)
        sec_or_ad(layer)                    = sec_or_ad(layer) + predictors_ad%ozone(1, layer, iprof)
        or_ad(layer)                        = or_ad(layer) + sec_or_ad(layer) * raytracing%pathsat(layer, iprof)
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) + sec_or_ad(layer) * or(layer)
      ENDDO

    ENDIF

!5.3 Water Vapour Continuum based on RTIASI
!------------------------------------------

    IF (coef%nwvcont > 0) THEN

      DO layer = 1, prof(1)%nlayers
        level              = layer + 1
        sec_wr_ad(layer)   = sec_wr_ad(layer) + predictors_ad%wvcont(4, layer, iprof) / tr_sq(layer)
        tr_ad(layer)       = tr_ad(layer) -                                                              &
          & 2._JPRB * predictors_ad%wvcont(4, layer, iprof) * predictors%watervapour(7, layer, iprof) /  &
          & (tr_sq(layer) * tr(layer))
        sec_wr_ad(layer)   = sec_wr_ad(layer) + predictors_ad%wvcont(3, layer, iprof) / tr(layer)
        tr_ad(layer)       =      &
          & tr_ad(layer) - predictors_ad%wvcont(3, layer, iprof) * predictors%watervapour(7, layer, iprof) / tr_sq(layer)
        sec_wrwr_ad(layer) = sec_wrwr_ad(layer) + predictors_ad%wvcont(2, layer, iprof) / tr_4(layer)
        tr_ad(layer)       = tr_ad(layer) -      &
          & 4._JPRB * predictors_ad%wvcont(2, layer, iprof) * predictors%wvcont(1, layer, iprof) / tr_4(layer)
        sec_wrwr_ad(layer) = sec_wrwr_ad(layer) + predictors_ad%wvcont(1, layer, iprof) / tr(layer)
        tr_ad(layer)       =      &
          & tr_ad(layer) - predictors_ad%wvcont(1, layer, iprof) * predictors%wvcont(1, layer, iprof) / tr(layer)
      ENDDO

    ENDIF

!
!5.2 water vapour based on RTIASI
!--------------------------------

    DO layer = 1, prof(1)%nlayers
      level = layer + 1
      sec_wr = raytracing%pathsat(layer, iprof) * wr(layer)
      sec_wrwr = sec_wr(layer) * wr(layer)
      raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +                            &
        & predictors_ad%watervapour(12, layer, iprof) * 0.5_JPRB / Sqrt(raytracing%pathsat(layer, iprof)) *  &
        & (wr(layer) ** 1.5_JPRB / wwr(layer))
      wr_ad(layer)                        = wr_ad(layer) +                                                                   &
        & predictors_ad%watervapour(12, layer, iprof) * sqrt(raytracing%pathsat(layer, iprof)) * 1.5_JPRB * wr(layer) ** 0.5_JPRB &
        &  / wwr(layer)
      wwr_ad(layer)                       = wwr_ad(layer) -                                                               &
        & predictors_ad%watervapour(12, layer, iprof) * sqrt(raytracing%pathsat(layer, iprof)) * wr(layer) ** 1.5_JPRB /  &
        & wwr(layer) ** 2
      sec_wrwr_ad(layer)                  =      &
        & sec_wrwr_ad(layer) + predictors_ad%watervapour(11, layer, iprof) / wwr(layer)
      wwr_ad(layer)                       =      &
        & wwr_ad(layer) - predictors_ad%watervapour(11, layer, iprof) * sec_wrwr(layer) / wwr(layer) ** 2_JPIM
      dt_ad(layer)                        =      &
        & dt_ad(layer) + predictors_ad%watervapour(10, layer, iprof) * predictors%watervapour(5, layer, iprof)
      sec_wr_ad(layer)                    = sec_wr_ad(layer) +      &
        & 0.5_JPRB * predictors_ad%watervapour(10, layer, iprof) * dt(layer) / predictors%watervapour(5, layer, iprof)
      sec_wr_ad(layer)                    =      &
        & sec_wr_ad(layer) + predictors_ad%watervapour(9, layer, iprof) * dtabs(layer) * dt(layer)
      dt_ad(layer)                        = dt_ad(layer) +      &
        & 2.0_jprb * predictors_ad%watervapour(9, layer, iprof) * predictors%watervapour(7, layer, iprof) * dtabs(layer)
      sec_wr_ad(layer)                    =      &
        & sec_wr_ad(layer) + 3._JPRB * predictors_ad%watervapour(8, layer, iprof) * predictors%watervapour(1, layer, iprof)
      sec_wr_ad(layer)                    = sec_wr_ad(layer) + predictors_ad%watervapour(7, layer, iprof)
      sec_wr_ad(layer)                    = sec_wr_ad(layer) +      &
        & 0.25_JPRB * predictors_ad%watervapour(6, layer, iprof) / predictors%watervapour(6, layer, iprof) ** 3_JPIM
      sec_wr_ad(layer)                    =      &
        & sec_wr_ad(layer) + 0.5_JPRB * predictors_ad%watervapour(5, layer, iprof) / predictors%watervapour(5, layer, iprof)
      dt_ad(layer)                        =      &
        & dt_ad(layer) + predictors_ad%watervapour(4, layer, iprof) * predictors%watervapour(7, layer, iprof)
      sec_wr_ad(layer)                    = sec_wr_ad(layer) + predictors_ad%watervapour(4, layer, iprof) * dt(layer)
      ww_ad(layer)                        = ww_ad(layer) +                                            &
        & predictors_ad%watervapour(3, layer, iprof) * 2._JPRB * predictors%watervapour(2, layer, iprof) *  &
        & raytracing%pathsat(layer, iprof)
      raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
        & predictors_ad%watervapour(3, layer, iprof) * 2._JPRB * predictors%watervapour(2, layer, iprof) * ww(layer)
      ww_ad(layer)                        =      &
        & ww_ad(layer) + predictors_ad%watervapour(2, layer, iprof) * raytracing%pathsat(layer, iprof)
      raytracing_ad%pathsat(layer, iprof) =      &
        & raytracing_ad%pathsat(layer, iprof) + predictors_ad%watervapour(2, layer, iprof) * ww(layer)
      sec_wr_ad(layer)                    =      &
        & sec_wr_ad(layer) + 2._JPRB * predictors_ad%watervapour(1, layer, iprof) * predictors%watervapour(7, layer, iprof)
      sec_wr_ad(layer)                    = sec_wr_ad(layer) + sec_wrwr_ad(layer) * wr(layer)
      wr_ad(layer)                        = wr_ad(layer) + sec_wrwr_ad(layer) * predictors%watervapour(7, layer, iprof)
      wr_ad(layer)                        = wr_ad(layer) + sec_wr_ad(layer) * raytracing%pathsat(layer, iprof)
      raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) + sec_wr_ad(layer) * wr(layer)
    ENDDO

!5.1 mixed gases
!---------------

    DO layer = 2, prof(1)%nlayers
      level = layer + 1
! X10
      raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +         &
        & predictors_ad%mixedgas(10, layer, iprof) * 0.5_JPRB * tw(layer) ** 0.25_JPRB /  &
        & predictors%mixedgas(9, layer, iprof)
      tw_ad(layer)                        = tw_ad(layer) +      &
        & predictors_ad%mixedgas(10, layer, iprof) * 0.25_JPRB * predictors%mixedgas(9, layer, iprof) / tw(layer) ** 0.75_JPRB
    ENDDO


    DO layer = 1, prof_ad(1)%nlayers
      level = layer + 1
! X9
      raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
        & predictors_ad%mixedgas(9, layer, iprof) * 0.5_JPRB / predictors%mixedgas(9, layer, iprof)
! X8
      raytracing_ad%pathsat(layer, iprof) =      &
        & raytracing_ad%pathsat(layer, iprof) + predictors_ad%mixedgas(8, layer, iprof) * tw(layer) / tr(layer)
      tw_ad(layer)                        =      &
        & tw_ad(layer) + predictors_ad%mixedgas(8, layer, iprof) * raytracing%pathsat(layer, iprof) / tr(layer)
      tr_ad(layer)                        = tr_ad(layer) -      &
        & predictors_ad%mixedgas(8, layer, iprof) * tw(layer) * raytracing%pathsat(layer, iprof) / tr(layer) ** 2_JPIM
! X7
      raytracing_ad%pathsat(layer, iprof) =      &
        & raytracing_ad%pathsat(layer, iprof) + predictors_ad%mixedgas(7, layer, iprof) * tw(layer)
      tw_ad(layer)                        =      &
        & tw_ad(layer) + predictors_ad%mixedgas(7, layer, iprof) * raytracing%pathsat(layer, iprof)
! X6
      tr_ad(layer)                        =      &
        & tr_ad(layer) + predictors_ad%mixedgas(6, layer, iprof) * 2._JPRB * predictors%mixedgas(5, layer, iprof)
! X5
      tr_ad(layer)                        = tr_ad(layer) + predictors_ad%mixedgas(5, layer, iprof)
! X4
      tr_ad(layer)                        =      &
        & tr_ad(layer) + predictors_ad%mixedgas(4, layer, iprof) * 2._JPRB * predictors%mixedgas(3, layer, iprof)
      raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
        & predictors_ad%mixedgas(4, layer, iprof) * predictors%mixedgas(5, layer, iprof) ** 2_JPIM
! X3
      tr_ad(layer)                        =      &
        & tr_ad(layer) + predictors_ad%mixedgas(3, layer, iprof) * raytracing%pathsat(layer, iprof)
      raytracing_ad%pathsat(layer, iprof) =      &
        & raytracing_ad%pathsat(layer, iprof) + predictors_ad%mixedgas(3, layer, iprof) * tr(layer)
! X2
      raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
        & predictors_ad%mixedgas(2, layer, iprof) * 2._JPRB * raytracing%pathsat(layer, iprof)
! X1
      raytracing_ad%pathsat(layer, iprof) =      &
        & raytracing_ad%pathsat(layer, iprof) + predictors_ad%mixedgas(1, layer, iprof)
    ENDDO

!-------------------------------------------------------------------
!   calc adjoint of profile/reference sums
!-------------------------------------------------------------------

    IF (opts%rt_ir%co2_Data .AND. coef%nco2 > 0) THEN
      sum1 = 0._JPRB

      DO layer = prof(1)%nlayers, 1,  - 1
        sum1          = sum1 + co2w_ad(layer) / sum2_co2w(layer)
        co2_ad(layer) = co2_ad(layer) + sum1 * coef%dpp(layer - 1)
      ENDDO
    ENDIF

    IF (coef%nco2 > 0) THEN
      sum1 = 0._JPRB

      DO layer = prof(1)%nlayers, 2,  - 1
        sum1            = sum1 + twr_ad(layer) / sum2_twr(layer)
        t_ad(layer - 1) = t_ad(layer - 1) + sum1 * coef%dpp(layer - 1)
      ENDDO
    ENDIF

!
    sum1 = 0._JPRB

    IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN

      DO layer = prof(1)%nlayers, 1,  - 1
        sum1        = sum1 + ow_ad(layer) / sum2_ow(layer)
        o_ad(layer) = o_ad(layer) + sum1 * coef%dpp(layer - 1)
      ENDDO

    ELSE
      o_ad(:) = 0._JPRB
    ENDIF

!
    sum1 = 0._JPRB

    DO layer = prof(1)%nlayers, 1,  - 1
      sum1        = sum1 + wwr_ad(layer) / sum2_wwr(layer)
      w_ad(layer) = w_ad(layer) + sum1 * coef%dpp(layer - 1) * t(layer)
      t_ad(layer) = t_ad(layer) + sum1 * coef%dpp(layer - 1) * w(layer)
!    write(48,*) layer, w_ad(layer), sum1, coef % dpp( layer-1 ), t( layer )
    ENDDO

!
    sum1 = 0._JPRB

    DO layer = prof(1)%nlayers, 1,  - 1
      sum1        = sum1 + ww_ad(layer) / sum2_ww(layer)
      w_ad(layer) = w_ad(layer) + sum1 * coef%dpp(layer - 1)
    ENDDO


    DO layer = prof(1)%nlayers, 2,  - 1
      tw_ad(layer - 1) = tw_ad(layer - 1) + tw_ad(layer)
      tr_ad(layer - 1) = tr_ad(layer - 1) + tw_ad(layer) * coef%dpp(layer - 1)
    ENDDO

!-------------------------------------------------------------------
!   calc adjoint of profile deviations
!-------------------------------------------------------------------

    DO layer = 1, prof(1)%nlayers

      IF (opts%rt_ir%co2_Data .AND. coef%nco2 > 0) THEN
        co2_ad(layer) = co2_ad(layer) + co2r_ad(layer) / coef%co2star(layer)
      ENDIF


      IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN
        o_ad(layer) = o_ad(layer) + or_ad(layer) / coef%ostar(layer)
      ENDIF

      w_ad(layer) = w_ad(layer) + wr_ad(layer) / coef%wstar(layer)
      t_ad(layer) = t_ad(layer) + tr_ad(layer) / coef%tstar(layer)
!    write(48,*) layer, w_ad(layer), 'w_ad'

      IF (coef%nozone > 0) t_ad(layer) = t_ad(layer) + dto_ad(layer)

      t_ad(layer) = t_ad(layer) + dt_ad(layer)
    ENDDO

!-------------------------------------------------------------------
!   calc adjoint of profile layer means
!-------------------------------------------------------------------

    IF (opts%rt_ir%co2_Data .AND. coef%nco2 > 0) THEN

      DO level = 2, prof(1)%nlevels
        layer = level - 1
        prof_ad(iprof)%co2(level - 1) = prof_ad(iprof)%co2(level - 1) + 0.5_JPRB * co2_ad(layer)
        prof_ad(iprof)%co2(level)     = prof_ad(iprof)%co2(level) + 0.5_JPRB * co2_ad(layer)
      ENDDO

    ENDIF


    IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN

      DO level = 2, prof(1)%nlevels
        layer = level - 1
        prof_ad(iprof)%o3(level - 1) = prof_ad(iprof)%o3(level - 1) + 0.5_JPRB * o_ad(layer)
        prof_ad(iprof)%o3(level)     = prof_ad(iprof)%o3(level) + 0.5_JPRB * o_ad(layer)
      ENDDO

    ENDIF


! include adjoint surface humidity

    IF (opts%rt_all%use_q2m) THEN

      ! It is not immediately obvious that this is correct, but after careful comparison with
      ! the TL one can see that it is.
      DO level = 2, prof(1)%nlevels
        layer = level - 1
        prof_ad(iprof)%q(level - 1) = prof_ad(iprof)%q(level - 1) + 0.5_JPRB * w_ad(layer)

        IF (level == iv2lev) THEN
          iv2lay = iv2lev - 1
          prof_ad(iprof)%s2m%q = prof_ad(iprof)%s2m%q + 0.5_JPRB * w_ad(iv2lay)
        ELSE
          prof_ad(iprof)%q(level) = prof_ad(iprof)%q(level) + 0.5_JPRB * w_ad(layer)
        ENDIF
      ENDDO

    ELSE

      DO level = 2, prof(1)%nlevels
        layer = level - 1
        prof_ad(iprof)%q(level - 1) = prof_ad(iprof)%q(level - 1) + 0.5_JPRB * w_ad(layer)
        prof_ad(iprof)%q(level)     = prof_ad(iprof)%q(level) + 0.5_JPRB * w_ad(layer)
      ENDDO

    ENDIF


    DO level = 2, prof(1)%nlevels
      layer = level - 1
      prof_ad(iprof)%t(level - 1) = prof_ad(iprof)%t(level - 1) + 0.5_JPRB * t_ad(layer)
      prof_ad(iprof)%t(level)     = prof_ad(iprof)%t(level) + 0.5_JPRB * t_ad(layer)
    ENDDO

  ENDDO

  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_8_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setpredictors_8_ad
