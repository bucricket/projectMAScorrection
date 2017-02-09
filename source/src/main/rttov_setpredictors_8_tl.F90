!
SUBROUTINE rttov_setpredictors_8_tl( &
            & opts,          &
            & prof,          &
            & prof_tl,       &
            & geom,          &
            & coef,          &
            & aux,           &
            & predictors,    &
            & predictors_tl, &
            & raytracing,    &
            & raytracing_tl)
!
! Description
! RTTOV-8 Model
! TL of rttov_setpredictors_8
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
!  1.1   17/09/2003  Added predictors for wv line and continuum and CO2 (R Saunders)
!  1.2   03/06/2004  Parkind parametrisation, correction of wwr_tl calculation
!                    simplify TL of predictor 9 for WVL (P. Brunel)
!  1.3   23/02/2005  Correction of Twr definition (P. Brunel)
!  1.4   29/03/2005  Add end of header comment (J. Cameron)
!  1.5   07/12/2005  Add surface humidity (R. Saunders)
!  1.6   03/03/2006  Marco Matricardi (ECMWF):
!           --       Introduced altitude dependent local zenith angle
!  1.7   14/04/2007  Corrected CO2 profile logic (R Saunders)
!  1.8   15/07/2009  User defined ToA. Layers distinct from levels (P.Rayer)
!  1.9   02/12/2009  Pathsat, Pathsun and related quantities are now layer arrays
!                    (Marco Matricardi).
!  1.10  10/11/2010  Remove rttov9_compat flag from code (J Hocking)
!  1.11  25/10/2011  Fix bug so that omitted trace gas profiles are treated 
!                    correctly (J Hocking)
!  1.12  23/11/2011  Fix bug in CO2 predictors (J Hocking)
!  1.13  10/01/2013  Add PMC shifts (P Rayer)
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
  TYPE(profile_Type   ), INTENT(IN)    :: prof_tl(size(prof))
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  TYPE(geometry_Type  ), INTENT(IN)    :: geom(size(prof))
  TYPE(rttov_path_pred), INTENT(IN)    :: predictors
  TYPE(rttov_path_pred), INTENT(INOUT) :: predictors_tl      ! in because of mem allocation
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(IN)    :: raytracing_tl
  TYPE(profile_aux    ), INTENT(IN)    :: aux                ! auxillary profiles info.
!INTF_END
!local variables:
  INTEGER(KIND=jpim) :: level , layer , iprof
  INTEGER(KIND=jpim) :: iv2lev, iv3lev, iv2lay
! user profile
  REAL   (KIND=Jprb) :: t(prof(1)%nlayers)
  REAL   (KIND=Jprb) :: w(prof(1)%nlayers)
  REAL   (KIND=jprb) :: o(prof(1)%nlayers)
  REAL   (KIND=jprb) :: co2(prof(1)%nlayers)
! reference profile
  REAL   (KIND=Jprb) :: tr   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: wr   (prof(1)%nlayers)
  REAL   (KIND=jprb) :: or   (prof(1)%nlayers)
  REAL   (KIND=jprb) :: co2r (prof(1)%nlayers)
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
  REAL   (KIND=Jprb) :: tr_sq      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tr_sqrt    (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tr_4       (prof(1)%nlayers)
  REAL   (KIND=jprb) :: sec_wr     (prof(1)%nlayers)
  REAL   (KIND=jprb) :: sec_wrwr   (prof(1)%nlayers)
  REAL   (KIND=jprb) :: sec_or     (prof(1)%nlayers)
! TL variables
  REAL   (KIND=Jprb) :: t_tl       (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: w_tl       (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: o_tl       (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: co2_tl     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tr_tl      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: wr_tl      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: wwr_tl     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: or_tl      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: co2r_tl    (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: twr_tl     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: dt_tl      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: dto_tl     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tw_tl      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: ww_tl      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: ow_tl      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: co2w_tl    (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sec_or_tl  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sec_wr_tl  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sec_wrwr_tl(prof(1)%nlayers)
  INTEGER(KIND=jpim) :: nprofiles                     ! Number of profiles

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
  REAL   (KIND=Jprb) :: acm_TL
  REAL   (KIND=Jprb) :: Pcel_Lev_TL
  REAL   (KIND=Jprb) :: Pnom_LevM1_TL
  REAL   (KIND=Jprb) :: Pcel_LevM1_TL
  REAL   (KIND=Jprb) :: Pnom_Lev_TL

  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
!-------------------------------------------------------------------------------
! Recompute direct variables
!-------------------------------------------------------------------------------
! 1) profile layer mean quantities
!------------------------------------------------------------------------------
! layer N-1 lies between levels N-1 and N
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_8_TL', 0_jpim, ZHOOK_HANDLE)
  nprofiles = size(prof)

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
    dt(:)    = t(:) - coef%tstar(:)
    dtabs(:) = Abs(dt(:))
!------------------------------------------------------------------------------
!3) calculate (profile / reference profile) ratios; tr wr or
! if no input O3 profile, set to reference value (or =1)
!------------------------------------------------------------------------------
! All assignments 1:prof(1) % nlayers
    tr(:)    = t(:) / coef%tstar(:)
    wr(:)    = w(:) / coef%wstar(:)

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
! 4. calculate profile / reference profile sums: tw wwr twr
!--------------------------------------------------------------------
    tw(1) = 0.

    DO layer = 2, prof(1)%nlayers
      tw(layer) = tw(layer - 1) + coef%dpp(layer - 1) * tr(layer - 1)
    ENDDO

    sum1 = 0._JPRB
    sum2 = 0._JPRB

    DO layer = 1, prof(1)%nlayers
      sum1       = sum1 + coef%dpp(layer - 1) * w(layer) * t(layer)
      sum2       = sum2 + coef%dpp(layer - 1) * coef%wstar(layer) * coef%tstar(layer)
      wwr(layer) = sum1 / sum2
    ENDDO

    sum1 = 0._JPRB
    sum2 = 0._JPRB

    DO layer = 1, prof(1)%nlayers
      sum1      = sum1 + coef%dpp(layer - 1) * w(layer)
      sum2      = sum2 + coef%dpp(layer - 1) * coef%wstar(layer)
      ww(layer) = sum1 / sum2
    ENDDO


    IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN
      sum1 = 0._JPRB
      sum2 = 0._JPRB

      DO layer = 1, prof(1)%nlayers
        sum1      = sum1 + coef%dpp(layer - 1) * o(layer)
        sum2      = sum2 + coef%dpp(layer - 1) * coef%ostar(layer)
        ow(layer) = sum1 / sum2
      ENDDO

    ELSE
      ow(:) = 1._JPRB
    ENDIF

    IF (opts%rt_ir%co2_Data .AND. coef%nco2 > 0) THEN
      sum1 = 0._JPRB
      sum2 = 0._JPRB

      DO layer = 1, prof(1)%nlayers
        sum1        = sum1 + coef%dpp(layer - 1) * co2(layer)
        sum2        = sum2 + coef%dpp(layer - 1) * coef%co2star(layer)
        co2w(layer) = sum1 / sum2
      ENDDO

    ELSE
      co2w(:) = 1._JPRB
    ENDIF

    IF (coef%nco2 > 0) THEN
      sum1   = 0._JPRB
      sum2   = 0._JPRB
      twr(1) = 1._JPRB

      DO layer = 2, prof(1)%nlayers
        sum1       = sum1 + coef%dpp(layer - 1) * t(layer - 1)
        sum2       = sum2 + coef%dpp(layer - 1) * coef%tstar(layer - 1)
        twr(layer) = sum1 / sum2
      ENDDO
    ENDIF

!-------------------------------------------------------------------------------
! Now compute TL variables
!-------------------------------------------------------------------------------
! 1) profile layer mean quantities
!------------------------------------------------------------------------------

    DO layer = 1, prof(1)%nlayers
      level       = layer + 1
      t_tl(layer) = (prof_tl(iprof)%t(level - 1) + prof_tl(iprof)%t(level)) / 2._JPRB
      w_tl(layer) = (prof_tl(iprof)%q(level - 1) + prof_tl(iprof)%q(level)) / 2._JPRB
    ENDDO

! include tl surface humidity

    IF (opts%rt_all%use_q2m) THEN

      IF (iv2lev <= coef%nlevels) THEN
        w_tl(iv2lay) = (prof_tl(iprof)%s2m%q + prof_tl(iprof)%q(iv3lev)) / 2._JPRB
      ENDIF

    ENDIF


    IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN

      DO layer = 1, prof(1)%nlayers
        level       = layer + 1
        o_tl(layer) = (prof_tl(iprof)%o3(level - 1) + prof_tl(iprof)%o3(level)) / 2._JPRB
      ENDDO

    ENDIF


    IF (opts%rt_ir%co2_Data .AND. coef%nco2 > 0) THEN

      DO layer = 1, prof(1)%nlayers
        level         = layer + 1
        co2_tl(layer) = (prof_tl(iprof)%co2(level - 1) + prof_tl(iprof)%co2(level)) / 2._JPRB
      ENDDO

    ENDIF


!------------------------------------------------------------------------------
!2) calculate deviations from reference profile (layers)
!------------------------------------------------------------------------------
! All assignments 1:prof(1) % nlayers
    dt_tl(:) = t_tl(:)

    IF (coef%nozone > 0) dto_tl(:) = t_tl(:)

!------------------------------------------------------------------------------
!3) calculate (profile / reference profile) ratios; tr_tl wr_tl or_tl
!------------------------------------------------------------------------------
! All assignments 1:prof(1) % nlayers
    tr_tl(:) = t_tl(:) / coef%tstar(:)
    wr_tl(:) = w_tl(:) / coef%wstar(:)

    IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN
      or_tl(:) = o_tl(:) / coef%ostar(:)
    ELSE
      or_tl(:) = 0._JPRB
    ENDIF


    IF (opts%rt_ir%co2_Data .AND. coef%nco2 > 0) THEN
      co2r_tl(:) = co2_tl(:) / coef%co2star(:)
    ELSE
      co2r_tl(:) = 0._JPRB
    ENDIF

!-------------------------------------------------------------------------
! 4. calculate profile / reference profile sums: tw_tl ww_tl ow_tl co2w_tl
!-------------------------------------------------------------------------
    tw_tl(1) = 0._JPRB

    DO layer = 2, prof(1)%nlayers
      tw_tl(layer) = tw_tl(layer - 1) + coef%dpp(layer - 1) * tr_tl(layer - 1)
    ENDDO

    sum1 = 0._JPRB
    sum2 = 0._JPRB

    DO layer = 1, prof(1)%nlayers
      sum1         = sum1 + coef%dpp(layer - 1) * w_tl(layer)
      sum2         = sum2 + coef%dpp(layer - 1) * coef%wstar(layer)
      ww_tl(layer) = sum1 / sum2
    ENDDO

    sum1 = 0._JPRB
    sum2 = 0._JPRB

    DO layer = 1, prof(1)%nlayers
      sum1          = sum1 + coef%dpp(layer - 1) * (w_tl(layer) * t(layer) + w(layer) * t_tl(layer))
      sum2          = sum2 + coef%dpp(layer - 1) * coef%wstar(layer) * coef%tstar(layer)
      wwr_tl(layer) = sum1 / sum2
    ENDDO


    IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN
      sum1 = 0._JPRB
      sum2 = 0._JPRB

      DO layer = 1, prof(1)%nlayers
        sum1         = sum1 + coef%dpp(layer - 1) * o_tl(layer)
        sum2         = sum2 + coef%dpp(layer - 1) * coef%ostar(layer)
        ow_tl(layer) = sum1 / sum2
      ENDDO

    ELSE
      ow_tl(:) = 0._JPRB
    ENDIF


    IF (opts%rt_ir%co2_Data .AND. coef%nco2 > 0) THEN
      sum1 = 0._JPRB
      sum2 = 0._JPRB

      DO layer = 1, prof(1)%nlayers
        sum1           = sum1 + coef%dpp(layer - 1) * co2_tl(layer)
        sum2           = sum2 + coef%dpp(layer - 1) * coef%co2star(layer)
        co2w_tl(layer) = sum1 / sum2
      ENDDO
    
    ELSE
      co2w_tl(:) = 0._JPRB
    ENDIF

    IF (coef%nco2 > 0) THEN
      sum1      = 0._JPRB
      sum2      = 0._JPRB
      twr_tl(1) = 0._JPRB

      DO layer = 2, prof(1)%nlayers
        sum1          = sum1 + coef%dpp(layer - 1) * t_tl(layer - 1)
        sum2          = sum2 + coef%dpp(layer - 1) * coef%tstar(layer - 1)
        twr_tl(layer) = sum1 / sum2
      ENDDO
    ENDIF
    
! End of TL profile calcs
! ATTENTION
!  w_tl(:) = prof_tl(iprof) % q(:)
!5) set predictors for RTTOV-8 options
!--
!5.1 mixed gases
!---

    DO layer = 1, prof(1)%nlayers
      level = layer + 1! NB  raytracing % pathsat(layer,i) is angle at lower boundary of layer
      predictors_tl%mixedgas(1, layer, iprof) = raytracing_tl%pathsat(layer, iprof)
      predictors_tl%mixedgas(2, layer, iprof) =      &
        & raytracing_tl%pathsat(layer, iprof) * 2._JPRB * raytracing%pathsat(layer, iprof)
      predictors_tl%mixedgas(3, layer, iprof) =      &
        & tr_tl(layer) * raytracing%pathsat(layer, iprof) + raytracing_tl%pathsat(layer, iprof) * tr(layer)
      predictors_tl%mixedgas(4, layer, iprof) = 2._JPRB * tr_tl(layer) * predictors%mixedgas(3, layer, iprof) +      &
        & raytracing_tl%pathsat(layer, iprof) * predictors%mixedgas(5, layer, iprof) ** 2_JPIM
      predictors_tl%mixedgas(5, layer, iprof) = tr_tl(layer)
      predictors_tl%mixedgas(6, layer, iprof) = 2._JPRB * tr_tl(layer) * tr(layer)
      predictors_tl%mixedgas(7, layer, iprof) =      &
        & raytracing_tl%pathsat(layer, iprof) * tw(layer) + tw_tl(layer) * raytracing%pathsat(layer, iprof)
      predictors_tl%mixedgas(8, layer, iprof) = raytracing_tl%pathsat(layer, iprof) * tw(layer) / tr(layer) +      &
        & tw_tl(layer) * raytracing%pathsat(layer, iprof) / tr(layer) -                                            &
        & tr_tl(layer) * tw(layer) * raytracing%pathsat(layer, iprof) / tr(layer) ** 2_JPIM
      predictors_tl%mixedgas(9, layer, iprof) =      &
        & raytracing_tl%pathsat(layer, iprof) * 0.5_JPRB / predictors%mixedgas(9, layer, iprof)

      IF (layer == 1) THEN
        predictors_tl%mixedgas(10, layer, iprof) = 0._JPRB
      ELSE
        predictors_tl%mixedgas(10, layer, iprof) =                                                                         &
          & raytracing_tl%pathsat(layer, iprof) * 0.5_JPRB * tw(layer) ** 0.25_JPRB / predictors%mixedgas(9, layer, iprof) &
          &  + tw_tl(layer) * 0.25_JPRB * predictors%mixedgas(9, layer, iprof) / tw(layer) ** 0.75_JPRB
      ENDIF

!5.2 water vapour lines based on RTIASI
!--------------------------------------
      sec_wr(layer)                               = raytracing%pathsat(layer, iprof) * wr(layer)
      sec_wrwr(layer)                             = sec_wr(layer) * wr(layer)
      sec_wr_tl(layer)                            =      &
        & raytracing%pathsat(layer, iprof) * wr_tl(layer) + raytracing_tl%pathsat(layer, iprof) * wr(layer)
      sec_wrwr_tl(layer)                          =      &
        & sec_wr_tl(layer) * wr(layer) + predictors%watervapour(7, layer, iprof) * wr_tl(layer)
      predictors_tl%watervapour(:, layer, iprof)  = 0._JPRB
      predictors_tl%watervapour(1, layer, iprof)  = 2._JPRB * predictors%watervapour(7, layer, iprof) * sec_wr_tl(layer)
      predictors_tl%watervapour(2, layer, iprof)  =      &
        & raytracing%pathsat(layer, iprof) * ww_tl(layer) + raytracing_tl%pathsat(layer, iprof) * ww(layer)
      predictors_tl%watervapour(3, layer, iprof)  =                                                              &
        & 2._JPRB * predictors%watervapour(2, layer, iprof) * raytracing%pathsat(layer, iprof) * ww_tl(layer) +  &
        & 2._JPRB * predictors%watervapour(2, layer, iprof) * raytracing_tl%pathsat(layer, iprof) * ww(layer)
      predictors_tl%watervapour(4, layer, iprof)  =      &
        & predictors%watervapour(7, layer, iprof) * dt_tl(layer) + sec_wr_tl(layer) * dt(layer)
      predictors_tl%watervapour(5, layer, iprof)  =      &
        & 0.5_JPRB * sec_wr_tl(layer) / predictors%watervapour(5, layer, iprof)
      predictors_tl%watervapour(6, layer, iprof)  =      &
        & 0.25_JPRB * sec_wr_tl(layer) / predictors%watervapour(6, layer, iprof) ** 3_JPIM
      predictors_tl%watervapour(7, layer, iprof)  = sec_wr_tl(layer)
      predictors_tl%watervapour(8, layer, iprof)  = 3._JPRB * predictors%watervapour(1, layer, iprof) * sec_wr_tl(layer)
! NB can we sort this next one out?
      predictors_tl%watervapour(9, layer, iprof)  =      &
        & dtabs(layer) * (sec_wr_tl(layer) * dt(layer) + 2._JPRB * predictors%watervapour(7, layer, iprof) * dt_tl(layer))
      predictors_tl%watervapour(10, layer, iprof) = predictors%watervapour(5, layer, iprof) * dt_tl(layer) +      &
        & 0.5_JPRB * dt(layer) * sec_wr_tl(layer) / predictors%watervapour(5, layer, iprof)
      predictors_tl%watervapour(11, layer, iprof) =      &
        & sec_wrwr_tl(layer) / wwr(layer) - wwr_tl(layer) * sec_wrwr(layer) / wwr(layer) ** 2_JPIM
      predictors_tl%watervapour(12, layer, iprof) =                                                      &
        & 0.5_JPRB * raytracing_tl%pathsat(layer, iprof) / Sqrt(raytracing%pathsat(layer, iprof)) *           &
        & (wr(layer) ** 1.5_JPRB / wwr(layer)) +                                                         &
        & sqrt(raytracing%pathsat(layer, iprof)) * 1.5_JPRB * wr_tl(layer) * wr(layer) ** 0.5_JPRB / wwr(layer) -  &
        & sqrt(raytracing%pathsat(layer, iprof)) * wwr_tl(layer) * wr(layer) ** 1.5_JPRB / wwr(layer) ** 2_JPIM
!
! predictors_tl % watervapour(12,layer,iprof)  = 1.5 * sec_wr(layer)**0.5 * wr_tl(layer) / wwr(layer) - &
!                                    &  * wr * sec_wr(layer)**0.5 * wwr_tl(layer) / (wwr(layer) * wwr(layer))
!
!5.3 water vapour continuum transmittance based on RTIASI
!--------------------------------------------------------
!

      IF (coef%nwvcont > 0) THEN
        tr_sq(layer)                          = tr(layer) * tr(layer)
        tr_4(layer)                           = tr_sq(layer) * tr_sq(layer)
!predictors_tl % wvcont(:,layer,iprof)  = 0._JPRB
        predictors_tl%wvcont(1, layer, iprof) =      &
          & sec_wrwr_tl(layer) / tr(layer) - predictors%wvcont(1, layer, iprof) * tr_tl(layer) / tr(layer)
        predictors_tl%wvcont(2, layer, iprof) =      &
          & sec_wrwr_tl(layer) / tr_4(layer) - 4._JPRB * predictors%wvcont(1, layer, iprof) * tr_tl(layer) / tr_4(layer)
        predictors_tl%wvcont(3, layer, iprof) =      &
          & sec_wr_tl(layer) / tr(layer) - predictors%watervapour(7, layer, iprof) * tr_tl(layer) / tr_sq(layer)
        predictors_tl%wvcont(4, layer, iprof) = sec_wr_tl(layer) / tr_sq(layer) -      &
          & 2._JPRB * predictors%watervapour(7, layer, iprof) * tr_tl(layer) / (tr_sq(layer) * tr(layer))
      ENDIF

!
!5.4 ozone
!---------

      IF (coef%nozone > 0) THEN
        sec_or(layer)                         = raytracing%pathsat(layer, iprof) * or(layer)
        sec_or_tl(layer)                      =      &
          & raytracing_tl%pathsat(layer, iprof) * or(layer) + or_tl(layer) * raytracing%pathsat(layer, iprof)
        predictors_tl%ozone(1, layer, iprof)  = sec_or_tl(layer)
        predictors_tl%ozone(2, layer, iprof)  = 0.5_JPRB * sec_or_tl(layer) / predictors%ozone(2, layer, iprof)
        predictors_tl%ozone(3, layer, iprof)  =                                                         &
          & sec_or_tl(layer) * predictors%ozone(3, layer, iprof) / predictors%ozone(1, layer, iprof) +  &
          & predictors%ozone(1, layer, iprof) * dto_tl(layer)
        predictors_tl%ozone(4, layer, iprof)  = 2._JPRB * sec_or_tl(layer) * predictors%ozone(1, layer, iprof)
        predictors_tl%ozone(5, layer, iprof)  = 0.5_JPRB * sec_or_tl(layer) * predictors%ozone(3, layer, iprof) /      &
          & (predictors%ozone(1, layer, iprof) * predictors%ozone(2, layer, iprof)) +                                  &
          & predictors%ozone(2, layer, iprof) * dto_tl(layer)
        predictors_tl%ozone(6, layer, iprof)  =                                                    &
          & sec_or_tl(layer) * or(layer) * ow(layer) + or_tl(layer) * sec_or(layer) * ow(layer) +  &
          & ow_tl(layer) * sec_or(layer) * or(layer)
        predictors_tl%ozone(7, layer, iprof)  = predictors_tl%ozone(2, layer, iprof) * or(layer) / ow(layer) +      &
          & or_tl(layer) * predictors%ozone(2, layer, iprof) / ow(layer) -                                          &
          & ow_tl(layer) * predictors%ozone(2, layer, iprof) * or(layer) / ow(layer) ** 2_JPIM
        predictors_tl%ozone(8, layer, iprof)  = sec_or_tl(layer) * ow(layer) + sec_or(layer) * ow_tl(layer)
        predictors_tl%ozone(9, layer, iprof)  =                                                           &
          & raytracing%pathsat(layer, iprof) * or_tl(layer) * Sqrt(predictors%ozone(10, layer, iprof)) +  &
          & predictors%ozone(1, layer, iprof) * 0.5_JPRB * raytracing%pathsat(layer, iprof) * ow_tl(layer) /   &
          & Sqrt(predictors%ozone(10, layer, iprof)) +                                                    &
          & 1.5_JPRB * raytracing_tl%pathsat(layer, iprof) * or(layer) * Sqrt(predictors%ozone(10, layer, iprof))
        predictors_tl%ozone(10, layer, iprof) =      &
          & raytracing_tl%pathsat(layer, iprof) * ow(layer) + raytracing%pathsat(layer, iprof) * ow_tl(layer)
        predictors_tl%ozone(11, layer, iprof) =                                                         &
          & 2._JPRB * ow_tl(layer) * raytracing%pathsat(layer, iprof) * predictors%ozone(10, layer, iprof) +  &
          & 2._JPRB * raytracing_tl%pathsat(layer, iprof) * ow(layer) * predictors%ozone(10, layer, iprof)
      ENDIF

!
!5.5 cloud
!---------

      IF (coef%id_sensor == sensor_id_mw) THEN
        IF (opts%rt_mw%clw_Data) THEN
          deltac(layer)                               =      &
            & 0.1820_JPRB * 100.0_JPRB * coef%dp(layer) / (4.3429_JPRB * gravity)
          predictors_tl%clw(layer, iprof)             = deltac(layer) * prof_tl(iprof)%clw(level) * geom(iprof)%seczen
          predictors_tl%clw(2:prof(1)%nlayers, iprof) = 0.5_JPRB * (predictors_tl%clw(2:prof(1)%nlayers, iprof) +      &
            & deltac(2:prof(1)%nlayers) * prof_tl(iprof)%clw(1:prof(1)%nlevels - 1) * geom(iprof)%seczen)
        ELSE
          predictors_tl%clw    = 0._JPRB
        ENDIF
      ENDIF
!
!5.6 carbon dioxide transmittance based on RTIASI
!-------------------------------------------------
!

      IF (coef%nco2 > 0) THEN
        tr_sqrt(layer)                      = Sqrt(tr(layer))
        predictors_tl%co2(1, layer, iprof)  = raytracing%pathsat(layer, iprof) * co2r_tl(layer) + &
          & raytracing_tl%pathsat(layer, iprof) * co2r(layer)
        predictors_tl%co2(2, layer, iprof)  = 2._JPRB * tr_tl(layer) * predictors%co2(5, layer, iprof)
        predictors_tl%co2(3, layer, iprof)  = raytracing%pathsat(layer, iprof) * tr_tl(layer) + &
          & raytracing_tl%pathsat(layer, iprof) * tr(layer)
        predictors_tl%co2(4, layer, iprof)  = 2._JPRB * tr_tl(layer) * predictors%co2(3, layer, iprof) + &
          & raytracing_tl%pathsat(layer, iprof) * predictors%co2(2, layer, iprof)
        predictors_tl%co2(5, layer, iprof)  = tr_tl(layer)
        predictors_tl%co2(6, layer, iprof)  = raytracing_tl%pathsat(layer, iprof)
        predictors_tl%co2(7, layer, iprof)  = raytracing%pathsat(layer, iprof) * twr_tl(layer) + &
          & raytracing_tl%pathsat(layer, iprof) * twr(layer)
        predictors_tl%co2(8, layer, iprof)  =      &
          & 2._JPRB * Sqrt(predictors%co2(8, layer, iprof)) * &
          & (raytracing%pathsat(layer, iprof) * co2w_tl(layer) + raytracing_tl%pathsat(layer, iprof) * co2w(layer))
        predictors_tl%co2(9, layer, iprof)  = 3._JPRB * twr(layer) * twr(layer) * twr_tl(layer)
        predictors_tl%co2(10, layer, iprof) = raytracing%pathsat(layer, iprof) * tr_sqrt(layer) * twr_tl(layer) + &
          & 0.5_JPRB * predictors%co2(7, layer, iprof) * tr_tl(layer) / tr_sqrt(layer) + &
          & raytracing_tl%pathsat(layer, iprof) * twr(layer) * tr_sqrt(layer)
      ENDIF

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

! This is on coef levels so pressure is not an active variable in the TL

! nlay (for pmc) may be less than prof(1)%nlayers
        IF(LAYER .LE. NLAY) THEN
! FWD
          acm = raytracing%co2_cm(iprof)/(2._JPRB*betaplus1*Lcel_cm)
! TL
          acm_TL = raytracing_TL%co2_cm(iprof)/(2._JPRB*betaplus1*Lcel_cm)

          DO ichan=1, coef%fmv_chn
! FWD
            Pcel_Lev= acm*raytracing%pathsat(layer,iprof)*prof(iprof)%p(level  )**2_JPIM + Pcel(ichan)**2_JPIM
            Pnom_LevM1=acm*raytracing%pathsat(layer,iprof)*prof(iprof)%p(level-1)**2_JPIM + Pnom(ichan)**2_JPIM
            Pcel_LevM1= acm*raytracing%pathsat(layer,iprof)*prof(iprof)%p(level-1)**2_JPIM + Pcel(ichan)**2_JPIM
            Pnom_Lev=acm*raytracing%pathsat(layer,iprof)*prof(iprof)%p(level  )**2_JPIM + Pnom(ichan)**2_JPIM
! TL
! PMC Predictor-1
            Pcel_Lev_TL  = acm_TL*raytracing%pathsat(layer,iprof)*prof(iprof)%p(level  )**2_JPIM                    &
         &              +acm*raytracing_TL%pathsat(layer,iprof)*prof(iprof)%p(level  )**2_JPIM
            Pnom_LevM1_TL = acm_TL*raytracing%pathsat(layer,iprof)*prof(iprof)%p(level-1)**2_JPIM                   &
         &              +acm*raytracing_TL%pathsat(layer,iprof)*prof(iprof)%p(level-1)**2_JPIM
            Pcel_LevM1_TL  = acm_TL*raytracing%pathsat(layer,iprof)*prof(iprof)%p(level-1)**2_JPIM                  &
         &              +acm*raytracing_TL%pathsat(layer,iprof)*prof(iprof)%p(level-1)**2_JPIM
            Pnom_Lev_TL = acm_TL*raytracing%pathsat(layer,iprof)*prof(iprof)%p(level  )**2_JPIM                     &
         &              +acm*raytracing_TL%pathsat(layer,iprof)*prof(iprof)%p(level  )**2_JPIM
            predictors_TL%pmc(1,layer,iprof,ichan)=                                                                 &
         &     Pcel_Lev_TL  * Pnom_LevM1/(Pcel_Lev*Pnom_LevM1)                                                      &
         &   + Pnom_LevM1_TL * Pcel_Lev/(Pcel_Lev*Pnom_LevM1)                                                       &
         &   - Pcel_LevM1_TL  * Pnom_Lev/(Pcel_LevM1*Pnom_Lev)                                                      &
         &   - Pnom_Lev_TL * Pcel_LevM1/(Pcel_LevM1*Pnom_Lev)
! PMC Predictor-2 - use layer depth as lower bound for denominator
            predictors_TL%pmc(2,layer,iprof,ichan)=0._JPRB

          ENDDO
        ENDIF

      ENDIF

    ENDDO
! layers
  ENDDO
! profiles

  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_8_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setpredictors_8_tl
