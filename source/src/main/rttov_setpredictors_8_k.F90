!
SUBROUTINE rttov_setpredictors_8_k( &
            & opts,         &
            & nlayers,      &
            & angles,       &
            & chanprof,     &
            & prof,         &
            & prof_k,       &
            & coef,         &
            & aux,          &
            & predictors,   &
            & predictors_k, &
            & raytracing,   &
            & raytracing_k)
! Description
! RTTOV-9 Model
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
!    Copyright 2005, EUMETSAT, All Rights Reserved.
!
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0       22/06/2005  initial (P Brunel)
!                        based on version 1.4 (29/03/05) of AD code
!  1.1       07/12/2005  Add surface humidity (R. Saunders)
!  1.2       03/03/2006  Marco Matricardi (ECMWF):
!           --           Introduced altitude dependent local zenith angle
!  1.3       16/08/2006  Merged with v87 Roger Saunders
!  1.4       12/02/2007  Removed polarisation (R Saunders)
!  1.5       14/03/2007  Corrected CO2 profile logic (R Saunders)
!  1.6       09/03/2009  Changed index on line 328 from j to i (R.Saunders)
!  1.7       15/09/2009  User defined ToA. Layers distinct from levels (P.Rayer)
!  1.8       02/12/2009  Pathsat, Pathsun and related quantities are now layer arrays
!                        (Marco Matricardi).
!  1.9       14/09/2010  Bug fix (P Brunel)
!  1.10      10/11/2010  Remove rttov9_compat flag from code (J Hocking)
!  1.11      25/10/2011  Fix bug so that omitted trace gas profiles are treated 
!                        correctly (J Hocking)
!  1.12      23/11/2011  Fix bug in CO2 predictors (J Hocking)
!  1.13      10/01/2013  Add PMC shifts (P Rayer)
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
       & rttov_chanprof,  &
       & rttov_coef,      &
       & rttov_options,   &
       & profile_Type,    &
       & geometry_Type,   &
       & profile_aux,     &
       & rttov_path_pred, &
       & raytracing_type
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE rttov_const, ONLY : gravity, sensor_id_mw
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options)  , INTENT(IN)    :: opts
  INTEGER(KIND=jpim)   , INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof ), INTENT(IN)    :: chanprof(:)
  TYPE(profile_Type   ), INTENT(IN)    :: prof(:)                 ! profile (ppmv dry)
  TYPE(profile_Type   ), INTENT(INOUT) :: prof_k(size(chanprof))
  TYPE(geometry_Type  ), INTENT(IN)    :: angles(size(prof)    )
  TYPE(rttov_path_pred), INTENT(IN)    :: predictors
  TYPE(rttov_path_pred), INTENT(INOUT) :: predictors_k
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing_k
  TYPE(profile_aux    ), INTENT(IN)    :: aux                     ! auxillary profiles info.
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
!INTF_END
!local variables:
  INTEGER(KIND=jpim) :: level , layer
  INTEGER(KIND=jpim) :: i                                  ! channel indice
  INTEGER(KIND=jpim) :: j                                  ! profile indice
  INTEGER(KIND=jpim) :: k                                  ! chanprof%chan indice
  INTEGER(KIND=jpim) :: iv2lev , iv3lev
  INTEGER(KIND=jpim) :: iv2lay 
! user profile
  REAL   (KIND=Jprb) :: t(nlayers, size(prof))
  REAL   (KIND=Jprb) :: w(nlayers, size(prof))
  REAL   (KIND=Jprb) :: o(nlayers, size(prof))
  REAL   (KIND=Jprb) :: co2  (nlayers, size(prof))
! reference profile
  REAL   (KIND=Jprb) :: tr   (nlayers, size(prof))
  REAL   (KIND=Jprb) :: wr   (nlayers, size(prof))
  REAL   (KIND=Jprb) :: or   (nlayers, size(prof))
  REAL   (KIND=Jprb) :: co2r (nlayers, size(prof))
  REAL   (KIND=Jprb) :: wwr  (nlayers, size(prof))
  REAL   (KIND=Jprb) :: twr  (nlayers, size(prof))
! user - reference
  REAL   (KIND=Jprb) :: dt   (nlayers, size(prof))
  REAL   (KIND=Jprb) :: dtabs(nlayers, size(prof))
! pressure weighted
  REAL   (KIND=Jprb) :: tw   (nlayers, size(prof))
  REAL   (KIND=jprb) :: ww   (nlayers, size(prof))
  REAL   (KIND=jprb) :: ow   (nlayers, size(prof))
  REAL   (KIND=jprb) :: co2w (nlayers, size(prof))
! intermediate variables
  REAL   (KIND=Jprb) :: sum1, sum2
  REAL   (KIND=Jprb) :: deltac    (nlayers                )
  REAL   (KIND=Jprb) :: sum2_ww   (nlayers, size(prof)    )
  REAL   (KIND=Jprb) :: sum2_wwr  (nlayers, size(prof)    )
  REAL   (KIND=Jprb) :: sum2_ow   (nlayers, size(prof)    )
  REAL   (KIND=Jprb) :: sum2_twr  (nlayers, size(prof)    )
  REAL   (KIND=Jprb) :: sum2_co2w (nlayers, size(prof)    )
  REAL   (KIND=Jprb) :: tr_sq     (nlayers, size(prof)    )
  REAL   (KIND=Jprb) :: tr_sqrt   (nlayers, size(prof)    )
  REAL   (KIND=Jprb) :: tr_4      (nlayers, size(prof)    )
  REAL   (KIND=jprb) :: sec_wrwr  (nlayers, size(prof)    )
  REAL   (KIND=jprb) :: sec_wr    (nlayers, size(prof)    )
! K variables
  REAL   (KIND=Jprb) :: t_k       (nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: w_k       (nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: o_k       (nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: co2_k     (nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: tr_k      (nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: wr_k      (nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: or_k      (nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: wwr_k     (nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: co2r_k    (nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: twr_k     (nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: dt_k      (nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: dto_k     (nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: tw_k      (nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: ww_k      (nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: ow_k      (nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: co2w_k    (nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: sec_or_k  (nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: sec_wr_k  (nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: sec_wrwr_k(nlayers, size(chanprof))

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
  REAL   (KIND=Jprb) :: acm_K
  REAL   (KIND=Jprb) :: Pcel_Lev_K
  REAL   (KIND=Jprb) :: Pnom_LevM1_K
  REAL   (KIND=Jprb) :: Pcel_LevM1_K
  REAL   (KIND=Jprb) :: Pnom_Lev_K

  INTEGER(KIND=jpim) :: nprofiles                          ! Number of profiles
  INTEGER(KIND=jpim) :: nchannels                          ! Number of radiances computed (channels used * profiles)

  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_8_K', 0_jpim, ZHOOK_HANDLE)
  nprofiles = size(prof)
  nchannels = size(chanprof)
!-------------------------------------------------------------------------------
! Recompute Direct variables
!-------------------------------------------------------------------------------

  DO j = 1, nprofiles
!1) Profile layer quantities
!-------------------------------------------------------------------------------

    DO layer = 1, prof(1)%nlayers
      level       = layer + 1
      t(layer, j) = (prof(j)%t(level - 1) + prof(j)%t(level)) * 0.5_JPRB
      w(layer, j) = (prof(j)%q(level - 1) + prof(j)%q(level)) * 0.5_JPRB
    ENDDO

! include surface humidity

    IF (opts%rt_all%use_q2m) THEN
! include surface humidity
      iv3lev = aux%s(j)%nearestlev_surf - 1
      iv2lev = aux%s(j)%nearestlev_surf

      IF (iv2lev <= coef%nlevels) THEN
        iv2lay = iv2lev - 1
        w(iv2lay, j) = (prof(j)%s2m%q + prof(j)%q(iv3lev)) * 0.5_JPRB
      ENDIF

    ENDIF

!

    DO layer = 1, prof(1)%nlayers
      level = layer + 1

      IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN
        o(layer, j) = (prof(j)%o3(level - 1) + prof(j)%o3(level)) * 0.5_JPRB
      ELSE
        o(layer, j) = 0.0_JPRB
      ENDIF

    ENDDO


    DO layer = 1, prof(1)%nlayers
      level = layer + 1

      IF (opts%rt_ir%co2_Data .AND. coef%nco2 > 0) THEN
        co2(layer, j) = (prof(j)%co2(level - 1) + prof(j)%co2(level)) * 0.5_JPRB
      ELSE
        co2(layer, j) = 0.0_JPRB
      ENDIF

    ENDDO


!------------------------------------------------------------------------------
!2) calculate deviations from reference profile (layers)
!------------------------------------------------------------------------------
! All assignments 1:prof(1) % nlayers
    dt(:, j)      = t(:, j) - coef%tstar(:)
    dtabs(:, j)   = Abs(dt(:, j))
!------------------------------------------------------------------------------
!3) calculate (profile / reference profile) ratios; tr wr or
! if no input O3 profile, set to reference value (or =1)
!------------------------------------------------------------------------------
! All assignments 1:prof(1) % nlayers
    tr(:, j)      = t(:, j) / coef%tstar(:)
    tr_sq(:, j)   = tr(:, j) * tr(:, j)
    tr_4(:, j)    = tr_sq(:, j) * tr_sq(:, j)
    tr_sqrt(:, j) = Sqrt(tr(:, j))
    wr(:, j)      = w(:, j) / coef%wstar(:)

    IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN
      or(:, j) = o(:, j) / coef%ostar(:)
    ELSE
      or(:, j) = 1._JPRB
    ENDIF

    IF (coef%nco2 > 0 .AND. opts%rt_ir%co2_Data) THEN
      co2r(:, j) = co2(:, j) / coef%co2star(:)
    ELSE
      co2r(:, j) = 1._JPRB
    ENDIF

!-------------------------------------------------------------------
! 4. calculate profile / reference profile sums: tw wwr
!--------------------------------------------------------------------
    tw(1, j) = 0._JPRB

    DO layer = 2, prof(j)%nlayers
      tw(layer, j) = tw(layer - 1, j) + coef%dpp(layer - 1) * tr(layer - 1, j)
    ENDDO

    sum1 = 0._JPRB
    sum2 = 0._JPRB

    DO layer = 1, prof(j)%nlayers
      sum1 = sum1 + coef%dpp(layer - 1) * w(layer, j) * t(layer, j)
      sum2 = sum2 + coef%dpp(layer - 1) * coef%wstar(layer) * coef%tstar(layer)
      sum2_wwr(layer, j) = sum2
      wwr(layer, j)      = sum1 / sum2
    ENDDO

    sum1 = 0._JPRB
    sum2 = 0._JPRB

    DO layer = 1, prof(j)%nlayers
      sum1 = sum1 + coef%dpp(layer - 1) * w(layer, j)
      sum2 = sum2 + coef%dpp(layer - 1) * coef%wstar(layer)
      sum2_ww(layer, j) = sum2
      ww(layer, j)      = sum1 / sum2
    ENDDO


    IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN
      sum1 = 0._JPRB
      sum2 = 0._JPRB

      DO layer = 1, prof(j)%nlayers
        sum1 = sum1 + coef%dpp(layer - 1) * o(layer, j)
        sum2 = sum2 + coef%dpp(layer - 1) * coef%ostar(layer)
        sum2_ow(layer, j) = sum2
        ow(layer, j)      = sum1 / sum2
      ENDDO

    ELSE
      sum2_ow(:, j) = 0._JPRB
      ow(:, j)      = 1._JPRB
    ENDIF


    IF (opts%rt_ir%co2_Data .AND. coef%nco2 > 0) THEN
      sum1 = 0._JPRB
      sum2 = 0._JPRB

      DO layer = 1, prof(j)%nlayers
        sum1 = sum1 + coef%dpp(layer - 1) * co2(layer, j)
        sum2 = sum2 + coef%dpp(layer - 1) * coef%co2star(layer)
        sum2_co2w(layer, j) = sum2
        co2w(layer, j)      = sum1 / sum2
      ENDDO
    ELSE
      sum2_co2w(:, j) = 0._JPRB
      co2w(:, j)      = 1._JPRB
    ENDIF

    IF (coef%nco2 > 0) THEN
      sum1 = 0._JPRB
      sum2 = 0._JPRB
      sum2_twr(1, j) = 0._JPRB
      twr(1, j)      = 1._JPRB
      
      DO layer = 2, prof(j)%nlayers
        sum1 = sum1 + coef%dpp(layer - 1) * t(layer - 1, j)
        sum2 = sum2 + coef%dpp(layer - 1) * coef%tstar(layer - 1)
        sum2_twr(layer, j) = sum2
        twr(layer, j)      = sum1 / sum2
      ENDDO
    ENDIF

  ENDDO

!-------------------------------------------------------------------------
! K code
!-------------------------------------------------------------------------
  w_k(:,:)        = 0._JPRB
  wr_k(:,:)       = 0._JPRB
  ww_k(:,:)       = 0._JPRB
  wwr_k(:,:)      = 0._JPRB
  sec_wr_k(:,:)   = 0._JPRB
  sec_wrwr_k(:,:) = 0._JPRB
  dt_k(:,:)       = 0._JPRB
  t_k(:,:)        = 0._JPRB
  tr_k(:,:)       = 0._JPRB
  tw_k(:,:)       = 0._JPRB

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

! This is on coef levels so pressure is not an active variable in the K

    DO i = 1, nchannels
       j = chanprof(i)%prof
       k = chanprof(i)%chan

      DO layer = 1, prof(1)%nlayers
        level=layer+1
! nlay (for pmc) may be less than prof(1)%nlayers
        IF(LAYER .LE. NLAY) THEN
          acm = raytracing%co2_cm(j)/(2._JPRB*betaplus1*Lcel_cm)
            Pcel_Lev= acm*raytracing%pathsat(layer,j)*prof(j)%p(level  )**2_JPIM + Pcel(k)**2_JPIM
            Pnom_LevM1=acm*raytracing%pathsat(layer,j)*prof(j)%p(level-1)**2_JPIM + Pnom(k)**2_JPIM
            Pcel_LevM1= acm*raytracing%pathsat(layer,j)*prof(j)%p(level-1)**2_JPIM + Pcel(k)**2_JPIM
            Pnom_Lev=acm*raytracing%pathsat(layer,j)*prof(j)%p(level  )**2_JPIM + Pnom(k)**2_JPIM
! AD
            Pcel_Lev_K =  0._JPRB
            Pnom_LevM1_K = 0._JPRB
            Pcel_LevM1_K =  0._JPRB
            Pnom_Lev_K = 0._JPRB
            acm_K = 0._JPRB
! for PMC predictor 2
            predictors_K%pmc(2,layer,j,k)=0._JPRB
! for PMC predictor 1
            Pcel_Lev_K = Pcel_Lev_K                                                    &
         &   + predictors_K%pmc(1,layer,j,k) * Pnom_LevM1/(Pcel_Lev*Pnom_LevM1)
            Pnom_LevM1_K = Pnom_LevM1_K                                                  &
         &   + predictors_K%pmc(1,layer,j,k) * Pcel_Lev/(Pcel_Lev*Pnom_LevM1)
            Pcel_LevM1_K = Pcel_LevM1_K                                                    &
         &   - predictors_K%pmc(1,layer,j,k) * Pnom_Lev/(Pcel_LevM1*Pnom_Lev)
            Pnom_Lev_K = Pnom_Lev_K                                                  &
         &   - predictors_K%pmc(1,layer,j,k) * Pcel_LevM1/(Pcel_LevM1*Pnom_Lev)
            predictors_K%pmc(1,layer,j,k) = 0._JPRB

            acm_K = acm_K                                                              &
         &   + Pnom_Lev_K * raytracing%pathsat(layer,j)*prof(j)%p(level  )**2_JPIM
            raytracing_K%pathsat(layer,i) = raytracing_K%pathsat(layer,i)      &
         &   + Pnom_Lev_K * acm * prof(j)%p(level  )**2_JPIM
            prof_K(i)%p(level  ) =  0._JPRB
            Pnom_Lev_K = 0._JPRB

            acm_K = acm_K                                                              &
         &   + Pcel_LevM1_K * raytracing%pathsat(layer,j)*prof(j)%p(level-1)**2_JPIM
            raytracing_K%pathsat(layer,i) = raytracing_K%pathsat(layer,i)      &
         &   + Pcel_LevM1_K * acm * prof(j)%p(level-1)**2_JPIM
            prof_K(i)%p(level-1) =  0._JPRB
            Pcel_LevM1_K = 0._JPRB

            acm_K = acm_K                                                              &
         &   + Pnom_LevM1_K * raytracing%pathsat(layer,j)*prof(j)%p(level-1)**2_JPIM
            raytracing_K%pathsat(layer,i) = raytracing_K%pathsat(layer,i)      &
         &   + Pnom_LevM1_K * acm * prof(j)%p(level-1)**2_JPIM
            prof_K(i)%p(level-1) =  0._JPRB
            Pnom_LevM1_K = 0._JPRB

            acm_K = acm_K                                                              &
         &   + Pcel_Lev_K * raytracing%pathsat(layer,j)*prof(j)%p(level  )**2_JPIM
            raytracing_K%pathsat(layer,i) = raytracing_K%pathsat(layer,i)      &
         &   + Pcel_Lev_K * acm * prof(j)%p(level  )**2_JPIM
            prof_K(i)%p(level-1) =  0._JPRB
            Pcel_Lev_K = 0._JPRB

            raytracing_K%co2_cm(i) = raytracing_K%co2_cm(i)                    &
         &   + acm_K/(2._JPRB*betaplus1*Lcel_cm)
            acm_K = 0._JPRB

        ENDIF
      ENDDO
    ENDDO

  ENDIF


!5.6 CO2
!-------

  IF (coef%nco2 > 0) THEN
    co2r_k(:,:) = 0._JPRB
    co2w_k(:,:) = 0._JPRB
    twr_k(:,:)  = 0._JPRB
    co2_k(:,:)  = 0._JPRB

    DO i = 1, nchannels
      j = chanprof(i)%prof

      DO layer = 1, prof_k(i)%nlayers
        level = layer + 1
        twr_k(layer, i)               =      &
          & twr_k(layer, i) + predictors_k%co2(10, layer, i) * raytracing%pathsat(layer, j) * tr_sqrt(layer, j)
        tr_k(layer, i)                =      &
          & tr_k(layer, i) + predictors_k%co2(10, layer, i) * 0.5_JPRB * predictors%co2(7, layer, j) / tr_sqrt(layer, j)
        raytracing_k%pathsat(layer, i) =     &
          & raytracing_k%pathsat(layer, i) + predictors_k%co2(10, layer, i) * twr(layer, j) * tr_sqrt(layer, j)
        twr_k(layer, i)               = twr_k(layer, i) +                                                     &
          & predictors_k%co2(9, layer, i) * 3._JPRB * raytracing%pathsat(layer, j) * predictors%co2(9, layer, j) /  &
          & predictors%co2(7, layer, j)
        co2w_k(layer, i)              = co2w_k(layer, i) +      &
          & 2._JPRB * raytracing%pathsat(layer, j) * predictors_k%co2(8, layer, i) * Sqrt(predictors%co2(8, layer, j))
        raytracing_k%pathsat(layer, i) = raytracing_k%pathsat(layer, i) + &
          & 2._JPRB * co2w(layer, j) * predictors_k%co2(8, layer, i) * &
          & Sqrt(predictors%co2(8, layer, j))
        twr_k(layer, i)               = twr_k(layer, i) + predictors_k%co2(7, layer, i) * raytracing%pathsat(layer, j)
        raytracing_k%pathsat(layer, i) =     &
          & raytracing_k%pathsat(layer, i) + predictors_k%co2(7, layer, i) * twr(layer, j)
        raytracing_k%pathsat(layer, i) =     &
          & raytracing_k%pathsat(layer, i) + predictors_k%co2(6, layer, i)
        tr_k(layer, i)                = tr_k(layer, i) + predictors_k%co2(5, layer, i)
        tr_k(layer, i)                =      &
          & tr_k(layer, i) + 2._JPRB * predictors_k%co2(4, layer, i) * predictors%co2(3, layer, j)
        raytracing_k%pathsat(layer, i) =     &
          & raytracing_k%pathsat(layer, i) + predictors_k%co2(4, layer, i) * predictors%co2(2, layer, j)
        tr_k(layer, i)                = tr_k(layer, i) + predictors_k%co2(3, layer, i) * raytracing%pathsat(layer, j)
        raytracing_k%pathsat(layer, i) =     &
          & raytracing_k%pathsat(layer, i) + predictors_k%co2(3, layer, i) * tr(layer, j)
        tr_k(layer, i)                =      &
          & tr_k(layer, i) + 2._JPRB * predictors_k%co2(2, layer, i) * predictors%co2(5, layer, j)
        co2r_k(layer, i)              =      &
          & co2r_k(layer, i) + predictors_k%co2(1, layer, i) * raytracing%pathsat(layer, j)
        raytracing_k%pathsat(layer, i) =     &
          & raytracing_k%pathsat(layer, i) + predictors_k%co2(1, layer, i) * co2r(layer, j)
      ENDDO
! layers
    ENDDO
! channels
  ENDIF
! coefs CO2
!5.5 cloud
!---------

  IF (coef%id_sensor == sensor_id_mw) THEN

    DO layer = 1, prof(1)%nlayers
      deltac(layer) = 0.1820_JPRB * 100.0_JPRB * coef%dp(layer) / (4.3429_JPRB * gravity)
    ENDDO


    DO i = 1, nchannels
      j = chanprof(i)%prof

      IF (opts%rt_mw%clw_Data) THEN

        DO layer = 2, prof_k(i)%nlayers
          level = layer + 1
          prof_k(i)%clw(level - 1)   =      &
            & prof_k(i)%clw(level - 1) + 0.5_JPRB * predictors_k%clw(layer, i) * deltac(layer) * angles(j)%seczen
          predictors_k%clw(layer, i) = 0.5_JPRB * predictors_k%clw(layer, i)
        ENDDO


        DO layer = 1, prof_k(i)%nlayers
          level                = layer + 1
          prof_k(i)%clw(level) = prof_k(i)%clw(level) + predictors_k%clw(layer, i) * deltac(layer) * angles(j)%seczen
        ENDDO

      ENDIF

    ENDDO

  ENDIF

!5.4 ozone
!---------

  IF (coef%nozone > 0) THEN
    o_k(:,:)      = 0._JPRB
    or_k(:,:)     = 0._JPRB
    wr_k(:,:)     = 0._JPRB
    ow_k(:,:)     = 0._JPRB
    dto_k(:,:)    = 0._JPRB
    sec_or_k(:,:) = 0._JPRB

    DO i = 1, nchannels
      j = chanprof(i)%prof
! One can pack all ow_k lines in one longer statement
! same for sec_or_k and dto_k

      DO layer = 1, prof_k(i)%nlayers
        level = layer + 1
        ow_k(layer, i)                 = ow_k(layer, i) +      &
          & predictors_k%ozone(11, layer, i) * 2._JPRB * raytracing%pathsat(layer, j) * predictors%ozone(10, layer, j)
        raytracing_k%pathsat(layer, i) = raytracing_k%pathsat(layer, i) +      &
          & predictors_k%ozone(11, layer, i) * 2._JPRB * ow(layer, j) * predictors%ozone(10, layer, j)
        ow_k(layer, i)                 =      &
          & ow_k(layer, i) + predictors_k%ozone(10, layer, i) * raytracing%pathsat(layer, j)
        raytracing_k%pathsat(layer, i) =      &
          & raytracing_k%pathsat(layer, i) + predictors_k%ozone(10, layer, i) * ow(layer, j)
        or_k(layer, i)                 = or_k(layer, i) +                                          &
          & predictors_k%ozone(9, layer, i) * SQRT(raytracing%pathsat(layer, j) * ow(layer, j)) *  &
          & raytracing%pathsat(layer, j)
        ow_k(layer, i)                 = ow_k(layer, i) +                                                           &
          & predictors_k%ozone(9, layer, i) * predictors%ozone(1, layer, j) * 0.5_JPRB * raytracing%pathsat(layer, j) /  &
          & SQRT(raytracing%pathsat(layer, j) * ow(layer, j))
        raytracing_k%pathsat(layer, i) = raytracing_k%pathsat(layer, i) +      &
          & predictors_k%ozone(9, layer, i) * 1.5_JPRB * or(layer, j) * SQRT(predictors%ozone(10, layer, j))
        or_k(layer, i)                 =      &
          & or_k(layer, i) + predictors_k%ozone(8, layer, i) * predictors%ozone(10, layer, j)
        ow_k(layer, i)                 =      &
          & ow_k(layer, i) + predictors_k%ozone(8, layer, i) * raytracing%pathsat(layer, j) * or(layer, j)
        raytracing_k%pathsat(layer, i) =      &
          & raytracing_k%pathsat(layer, i) + predictors_k%ozone(8, layer, i) * or(layer, j) * ow(layer, j)
        or_k(layer, i)                 =      &
          & or_k(layer, i) + predictors_k%ozone(7, layer, i) * 1.5_JPRB * predictors%ozone(2, layer, j) / ow(layer, j)
        ow_k(layer, i)                 =      &
          & ow_k(layer, i) - predictors_k%ozone(7, layer, i) * predictors%ozone(7, layer, j) / ow(layer, j)
        raytracing_k%pathsat(layer, i) = raytracing_k%pathsat(layer, i) +      &
          & predictors_k%ozone(7, layer, i) * 0.5_JPRB * or(layer, j) ** 1.5_JPRB /      &
          & (raytracing%pathsat(layer, j) ** 0.5_JPRB * ow(layer, j))
        or_k(layer, i)                 =      &
          & or_k(layer, i) + predictors_k%ozone(6, layer, i) * 2._JPRB * predictors%ozone(1, layer, j) * ow(layer, j)
        ow_k(layer, i)                 =      &
          & ow_k(layer, i) + predictors_k%ozone(6, layer, i) * predictors%ozone(4, layer, j) / raytracing%pathsat(layer, j)
        raytracing_k%pathsat(layer, i) =      &
          & raytracing_k%pathsat(layer, i) + predictors_k%ozone(6, layer, i) * or(layer, j) ** 2 * ow(layer, j)
        sec_or_k(layer, i)             = sec_or_k(layer, i) +                             &
          & predictors_k%ozone(5, layer, i) * 0.5_JPRB * predictors%ozone(3, layer, j) /  &
          & (predictors%ozone(1, layer, j) * predictors%ozone(2, layer, j))
        dto_k(layer, i)                =      &
          & dto_k(layer, i) + predictors_k%ozone(5, layer, i) * predictors%ozone(2, layer, j)
        sec_or_k(layer, i)             =      &
          & sec_or_k(layer, i) + predictors_k%ozone(4, layer, i) * 2._JPRB * predictors%ozone(1, layer, j)
        sec_or_k(layer, i)             = sec_or_k(layer, i) +      &
          & predictors_k%ozone(3, layer, i) * predictors%ozone(3, layer, j) / predictors%ozone(1, layer, j)
        dto_k(layer, i)                =      &
          & dto_k(layer, i) + predictors_k%ozone(3, layer, i) * predictors%ozone(1, layer, j)
        sec_or_k(layer, i)             =      &
          & sec_or_k(layer, i) + predictors_k%ozone(2, layer, i) * 0.5_JPRB / predictors%ozone(2, layer, j)
        sec_or_k(layer, i)             = sec_or_k(layer, i) + predictors_k%ozone(1, layer, i)
        or_k(layer, i)                 = or_k(layer, i) + sec_or_k(layer, i) * raytracing%pathsat(layer, j)
        raytracing_k%pathsat(layer, i) = raytracing_k%pathsat(layer, i) + sec_or_k(layer, i) * or(layer, j)
      ENDDO
! layers
    ENDDO
! channels
  ENDIF
! Coef O3
!5.3 Water Vapour Continuum based on RTIASI
!------------------------------------------

  IF (coef%nwvcont > 0) THEN

    DO i = 1, nchannels
      j = chanprof(i)%prof

      DO layer = 1, prof_k(i)%nlayers
        level                = layer + 1
        sec_wr_k(layer, i)   = sec_wr_k(layer, i) + predictors_k%wvcont(4, layer, i) / tr_sq(layer, j)
        tr_k(layer, i)       = tr_k(layer, i) -      &
          & 2._JPRB * predictors_k%wvcont(4, layer, i) * predictors%watervapour(7, layer, j) / (tr_sq(layer, j) * tr(layer, j))
        sec_wr_k(layer, i)   = sec_wr_k(layer, i) + predictors_k%wvcont(3, layer, i) / tr(layer, j)
        tr_k(layer, i)       =      &
          & tr_k(layer, i) - predictors_k%wvcont(3, layer, i) * predictors%watervapour(7, layer, j) / tr_sq(layer, j)
        sec_wrwr_k(layer, i) = sec_wrwr_k(layer, i) + predictors_k%wvcont(2, layer, i) / tr_4(layer, j)
        tr_k(layer, i)       =      &
          & tr_k(layer, i) - 4._JPRB * predictors_k%wvcont(2, layer, i) * predictors%wvcont(1, layer, j) / tr_4(layer, j)
        sec_wrwr_k(layer, i) = sec_wrwr_k(layer, i) + predictors_k%wvcont(1, layer, i) / tr(layer, j)
        tr_k(layer, i)       =      &
          & tr_k(layer, i) - predictors_k%wvcont(1, layer, i) * predictors%wvcont(1, layer, j) / tr(layer, j)
      ENDDO

    ENDDO
! channels
  ENDIF
! Coefs WV Cont
!
!5.2 water vapour based on RTIASI
!--------------------------------

  DO i = 1, nchannels
    j = chanprof(i)%prof

    DO layer = 1, prof_k(i)%nlayers
      level = layer + 1
      sec_wr(layer, j)               = raytracing%pathsat(layer, j) * wr(layer, j)
      sec_wrwr(layer, j)             = sec_wr(layer, j) * wr(layer, j)
      raytracing_k%pathsat(layer, i) = raytracing_k%pathsat(layer, i) +                             &
        & predictors_k%watervapour(12, layer, i) * 0.5_JPRB / Sqrt(raytracing%pathsat(layer, j)) *  &
        & (wr(layer, j) ** 1.5_JPRB / wwr(layer, j))
      wr_k(layer, i)                 = wr_k(layer, i) +                                                                   &
        & predictors_k%watervapour(12, layer, i) * sqrt(raytracing%pathsat(layer, j)) * 1.5_JPRB * wr(layer, j) ** 0.5_JPRB /  &
        & wwr(layer, j)
      wwr_k(layer, i)                = wwr_k(layer, i) -                                                       &
        & predictors_k%watervapour(12, layer, i) * sqrt(raytracing%pathsat(layer, j)) * wr(layer, j) ** 1.5_JPRB /  &
        & wwr(layer, j) ** 2
      sec_wrwr_k(layer, i)           = sec_wrwr_k(layer, i) + predictors_k%watervapour(11, layer, i) / wwr(layer, j)
      wwr_k(layer, i)                =      &
        & wwr_k(layer, i) - predictors_k%watervapour(11, layer, i) * sec_wrwr(layer, j) / wwr(layer, j) ** 2_JPIM
      dt_k(layer, i)                 =      &
        & dt_k(layer, i) + predictors_k%watervapour(10, layer, i) * predictors%watervapour(5, layer, j)
      sec_wr_k(layer, i)             = sec_wr_k(layer, i) +      &
        & 0.5_JPRB * predictors_k%watervapour(10, layer, i) * dt(layer, j) / predictors%watervapour(5, layer, j)
      sec_wr_k(layer, i)             =      &
        & sec_wr_k(layer, i) + predictors_k%watervapour(9, layer, i) * dtabs(layer, j) * dt(layer, j)
      dt_k(layer, i)                 = dt_k(layer, i) +      &
        & 2._JPRB * predictors_k%watervapour(9, layer, i) * predictors%watervapour(7, layer, j) * dtabs(layer, j)
      sec_wr_k(layer, i)             =      &
        & sec_wr_k(layer, i) + 3._JPRB * predictors_k%watervapour(8, layer, i) * predictors%watervapour(1, layer, j)
      sec_wr_k(layer, i)             = sec_wr_k(layer, i) + predictors_k%watervapour(7, layer, i)
      sec_wr_k(layer, i)             =      &
        & sec_wr_k(layer, i) + 0.25_JPRB * predictors_k%watervapour(6, layer, i) / predictors%watervapour(6, layer, j) ** 3_JPIM
      sec_wr_k(layer, i)             =      &
        & sec_wr_k(layer, i) + 0.5_JPRB * predictors_k%watervapour(5, layer, i) / predictors%watervapour(5, layer, j)
      dt_k(layer, i)                 =      &
        & dt_k(layer, i) + predictors_k%watervapour(4, layer, i) * predictors%watervapour(7, layer, j)
      sec_wr_k(layer, i)             = sec_wr_k(layer, i) + predictors_k%watervapour(4, layer, i) * dt(layer, j)
      ww_k(layer, i)                 = ww_k(layer, i) +      &
        & predictors_k%watervapour(3, layer, i) * 2._JPRB * predictors%watervapour(2, layer, j) * raytracing%pathsat(layer, j)
      raytracing_k%pathsat(layer, i) = raytracing_k%pathsat(layer, i) +      &
        & predictors_k%watervapour(3, layer, i) * 2._JPRB * predictors%watervapour(2, layer, j) * ww(layer, j)
      ww_k(layer, i)                 =      &
        & ww_k(layer, i) + predictors_k%watervapour(2, layer, i) * raytracing%pathsat(layer, j)
      raytracing_k%pathsat(layer, i) =      &
        & raytracing_k%pathsat(layer, i) + predictors_k%watervapour(2, layer, i) * ww(layer, j)
      sec_wr_k(layer, i)             =      &
        & sec_wr_k(layer, i) + 2._JPRB * predictors_k%watervapour(1, layer, i) * predictors%watervapour(7, layer, j)
      sec_wr_k(layer, i)             = sec_wr_k(layer, i) + sec_wrwr_k(layer, i) * wr(layer, j)
      wr_k(layer, i)                 = wr_k(layer, i) + sec_wrwr_k(layer, i) * predictors%watervapour(7, layer, j)
      wr_k(layer, i)                 = wr_k(layer, i) + sec_wr_k(layer, i) * raytracing%pathsat(layer, j)
      raytracing_k%pathsat(layer, i) = raytracing_k%pathsat(layer, i) + sec_wr_k(layer, i) * wr(layer, j)
    ENDDO
! layers
  ENDDO
! channels
!5.1 mixed gases
!---------------

  DO i = 1, nchannels
    j = chanprof(i)%prof

    DO layer = 2, prof(1)%nlayers
      level = layer + 1
! X10
      raytracing_k%pathsat(layer, i) = raytracing_k%pathsat(layer, i) +      &
        & predictors_k%mixedgas(10, layer, i) * 0.5_JPRB * tw(layer, j) ** 0.25_JPRB / predictors%mixedgas(9, layer, j)
      tw_k(layer, i)                 = tw_k(layer, i) +      &
        & predictors_k%mixedgas(10, layer, i) * 0.25_JPRB * predictors%mixedgas(9, layer, j) / tw(layer, j) ** 0.75_JPRB
    ENDDO

  ENDDO


  DO i = 1, nchannels
    j = chanprof(i)%prof

    DO layer = 1, prof_k(i)%nlayers
      level = layer + 1
! X9
      raytracing_k%pathsat(layer, i) =      &
        & raytracing_k%pathsat(layer, i) + predictors_k%mixedgas(9, layer, i) * 0.5_JPRB / predictors%mixedgas(9, layer, j)
! X8
      raytracing_k%pathsat(layer, i) =      &
        & raytracing_k%pathsat(layer, i) + predictors_k%mixedgas(8, layer, i) * tw(layer, j) / tr(layer, j)
      tw_k(layer, i)                 =      &
        & tw_k(layer, i) + predictors_k%mixedgas(8, layer, i) * raytracing%pathsat(layer, j) / tr(layer, j)
      tr_k(layer, i)                 = tr_k(layer, i) -      &
        & predictors_k%mixedgas(8, layer, i) * tw(layer, j) * raytracing%pathsat(layer, j) / tr(layer, j) ** 2_JPIM
! X7
      raytracing_k%pathsat(layer, i) =      &
        & raytracing_k%pathsat(layer, i) + predictors_k%mixedgas(7, layer, i) * tw(layer, j)
      tw_k(layer, i)                 =      &
        & tw_k(layer, i) + predictors_k%mixedgas(7, layer, i) * raytracing%pathsat(layer, j)
! X6
      tr_k(layer, i)                 =      &
        & tr_k(layer, i) + predictors_k%mixedgas(6, layer, i) * 2._JPRB * predictors%mixedgas(5, layer, j)
! X5
      tr_k(layer, i)                 = tr_k(layer, i) + predictors_k%mixedgas(5, layer, i)
! X4
      tr_k(layer, i)                 =      &
        & tr_k(layer, i) + predictors_k%mixedgas(4, layer, i) * 2._JPRB * predictors%mixedgas(3, layer, j)
      raytracing_k%pathsat(layer, i) =      &
        & raytracing_k%pathsat(layer, i) + predictors_k%mixedgas(4, layer, i) * predictors%mixedgas(5, layer, j) ** 2_JPIM
! X3
      tr_k(layer, i)                 =      &
        & tr_k(layer, i) + predictors_k%mixedgas(3, layer, i) * raytracing%pathsat(layer, j)
      raytracing_k%pathsat(layer, i) =      &
        & raytracing_k%pathsat(layer, i) + predictors_k%mixedgas(3, layer, i) * tr(layer, j)
! X2
      raytracing_k%pathsat(layer, i) =      &
        & raytracing_k%pathsat(layer, i) + predictors_k%mixedgas(2, layer, i) * 2._JPRB * raytracing%pathsat(layer, j)
! X1
      raytracing_k%pathsat(layer, i) = raytracing_k%pathsat(layer, i) + predictors_k%mixedgas(1, layer, i)
    ENDDO
! layers
  ENDDO
! channels

  DO i = 1, nchannels
    j = chanprof(i)%prof
!-------------------------------------------------------------------
!   calc K of profile/reference sums
!-------------------------------------------------------------------

    IF (opts%rt_ir%co2_Data .AND. coef%nco2 > 0) THEN
      sum1 = 0._JPRB

      DO layer = prof_k(i)%nlayers, 1,  - 1
        sum1            = sum1 + co2w_k(layer, i) / sum2_co2w(layer, j)
        co2_k(layer, i) = co2_k(layer, i) + sum1 * coef%dpp(layer - 1)
      ENDDO
    ENDIF

    IF (coef%nco2 > 0) THEN
      sum1 = 0._JPRB

      DO layer = prof_k(i)%nlayers, 2,  - 1
        sum1 = sum1 + twr_k(layer, i) / sum2_twr(layer, j)
        t_k(layer - 1, i) = t_k(layer - 1, i) + sum1 * coef%dpp(layer - 1)
      ENDDO
    ENDIF

!
    sum1 = 0._JPRB

    IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN

      DO layer = prof_k(i)%nlayers, 1,  - 1
        sum1          = sum1 + ow_k(layer, i) / sum2_ow(layer, j)
        o_k(layer, i) = o_k(layer, i) + sum1 * coef%dpp(layer - 1)
      ENDDO

    ELSE
      o_k(:, i) = 0._JPRB
    ENDIF

!
    sum1 = 0._JPRB

    DO layer = prof_k(i)%nlayers, 1,  - 1
      sum1          = sum1 + wwr_k(layer, i) / sum2_wwr(layer, j)
      w_k(layer, i) = w_k(layer, i) + sum1 * coef%dpp(layer - 1) * t(layer, j)
      t_k(layer, i) = t_k(layer, i) + sum1 * coef%dpp(layer - 1) * w(layer, j)
    ENDDO

!
    sum1 = 0._JPRB

    DO layer = prof_k(i)%nlayers, 1,  - 1
      sum1          = sum1 + ww_k(layer, i) / sum2_ww(layer, j)
      w_k(layer, i) = w_k(layer, i) + sum1 * coef%dpp(layer - 1)
    ENDDO


    DO layer = prof_k(i)%nlayers, 2,  - 1
      tw_k(layer - 1, i) = tw_k(layer - 1, i) + tw_k(layer, i)
      tr_k(layer - 1, i) = tr_k(layer - 1, i) + tw_k(layer, i) * coef%dpp(layer - 1)
    ENDDO

!-------------------------------------------------------------------
!   calc adjoint of profile deviations
!-------------------------------------------------------------------

    DO layer = 1, prof(1)%nlayers

      IF (opts%rt_ir%co2_Data .AND. coef%nco2 > 0) THEN
        co2_k(layer, i) = co2_k(layer, i) + co2r_k(layer, i) / coef%co2star(layer)
      ENDIF


      IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN
        o_k(layer, i) = o_k(layer, i) + or_k(layer, i) / coef%ostar(layer)
      ENDIF

      w_k(layer, i) = w_k(layer, i) + wr_k(layer, i) / coef%wstar(layer)
      t_k(layer, i) = t_k(layer, i) + tr_k(layer, i) / coef%tstar(layer)

      IF (coef%nozone > 0) t_k(layer, i) = t_k(layer, i) + dto_k(layer, i)

      t_k(layer, i) = t_k(layer, i) + dt_k(layer, i)
    ENDDO

!-------------------------------------------------------------------
!   calc adjoint of profile layer means
!-------------------------------------------------------------------

    IF (opts%rt_ir%co2_Data .AND. coef%nco2 > 0) THEN

      DO level = 2, prof(1)%nlevels
        layer = level - 1
        prof_k(i)%co2(level - 1) = prof_k(i)%co2(level - 1) + 0.5_JPRB * co2_k(layer, i)
        prof_k(i)%co2(level)     = prof_k(i)%co2(level) + 0.5_JPRB * co2_k(layer, i)
      ENDDO

    ENDIF


    IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN

      DO level = 2, prof(1)%nlevels
        layer = level - 1
        prof_k(i)%o3(level - 1) = prof_k(i)%o3(level - 1) + 0.5_JPRB * o_k(layer, i)
        prof_k(i)%o3(level)     = prof_k(i)%o3(level) + 0.5_JPRB * o_k(layer, i)
      ENDDO

    ENDIF


! include K surface humidity

    IF (opts%rt_all%use_q2m) THEN
      iv2lev = aux%s(j)%nearestlev_surf

      ! It is not immediately obvious that this is correct, but after careful comparison with
      ! the TL one can see that it is.
      DO level = 2, prof(1)%nlevels
        layer = level - 1
        prof_k(i)%q(level - 1) = prof_k(i)%q(level - 1) + 0.5_JPRB * w_k(layer, i)

        IF (level == iv2lev) THEN
          iv2lay = iv2lev - 1
          prof_k(i)%s2m%q = prof_k(i)%s2m%q + 0.5_JPRB * w_k(iv2lay, i)
        ELSE
          prof_k(i)%q(level) = prof_k(i)%q(level) + 0.5_JPRB * w_k(layer, i)
        ENDIF
      ENDDO

    ELSE

      DO level = 2, prof(1)%nlevels
        layer = level - 1
        prof_k(i)%q(level - 1) = prof_k(i)%q(level - 1) + 0.5_JPRB * w_k(layer, i)
        prof_k(i)%q(level)     = prof_k(i)%q(level) + 0.5_JPRB * w_k(layer, i)
      ENDDO

    ENDIF


    DO level = 2, prof(1)%nlevels
      layer = level - 1
      prof_k(i)%t(level - 1) = prof_k(i)%t(level - 1) + 0.5_JPRB * t_k(layer, i)
      prof_k(i)%t(level)     = prof_k(i)%t(level) + 0.5_JPRB * t_k(layer, i)
    ENDDO

  ENDDO
! channels
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_8_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setpredictors_8_k
