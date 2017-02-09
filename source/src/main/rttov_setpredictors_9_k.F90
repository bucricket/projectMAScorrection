!
SUBROUTINE rttov_setpredictors_9_k( &
            & opts,         &
            & chanprof,     &
            & nlayers,      &
            & prof,         &
            & prof_k,       &
            & geom,         &
            & ray_path,     &
            & ray_path_k,   &
            & coef_pccomp,  &
            & coef,         &
            & predictors,   &
            & predictors_k)
! Description
! K of rttov_setpredictors_9
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
!  1.0   16/08/2006  Roger Saunders
!           --       New routine based on rttov_setpredictors_k.F90.
!           --       Variable trace gases CO2, N2O,CO and CH4
!           --       introduced for IASI and AIRS.
!           --       Altitude dependent local zenith angle also
!           --       introduced.
!  1.1   12/02/2007  Removed polarisation index (R Saunders)
!  1.2   15/09/2009  User defined ToA. Layers distinct from levels (P.Rayer)
!  1.3   02/12/2009  Introduced principal component capability. Pathsat, Pathsun and
!                    related quantities are now layer arrays (Marco Matricardi).
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
       & rttov_chanprof,    &
       & rttov_coef,        &
       & rttov_options,     &
       & rttov_coef_pccomp, &
       & profile_Type,      &
       & geometry_Type,     &
       & rttov_path_pred
  USE parkind1, ONLY : jpim, jprb
!INTF_OFF
  USE rttov_const, ONLY : gravity, sensor_id_mw
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options)    , INTENT(IN)    :: opts
  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof(:)
  INTEGER(KIND=jpim)     , INTENT(IN)    :: nlayers
  TYPE(profile_Type     ), INTENT(IN)    :: prof(:)         ! profile (ppmv dry)
  TYPE(profile_Type     ), INTENT(INOUT) :: prof_k(SIZE(chanprof))
  TYPE(geometry_Type    ), INTENT(IN)    :: geom(SIZE(prof))
  REAL(jprb             ), INTENT(IN)    :: ray_path(prof(1)%nlayers,SIZE(prof))
  REAL(jprb             ), INTENT(INOUT) :: ray_path_k(prof(1)%nlayers,SIZE(chanprof))
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  TYPE(rttov_coef       ), INTENT(IN)    :: coef
  TYPE(rttov_path_pred  ), INTENT(IN)    :: predictors
  TYPE(rttov_path_pred  ), INTENT(INOUT) :: predictors_k
!INTF_END
!local variables:
  INTEGER(KIND=jpim) :: level, layer
  INTEGER(KIND=jpim) :: i
  INTEGER(KIND=jpim) :: j
! user profile
  REAL   (KIND=jprb) :: t(nlayers, size(prof))
  REAL   (KIND=jprb) :: w(nlayers, size(prof))
  REAL   (KIND=jprb) :: o(nlayers, size(prof))
  REAL   (KIND=jprb) :: co2  (nlayers, size(prof))
  REAL   (KIND=jprb) :: co   (nlayers, size(prof))
  REAL   (KIND=jprb) :: n2o  (nlayers, size(prof))
  REAL   (KIND=jprb) :: ch4  (nlayers, size(prof))
! reference profile
  REAL   (KIND=jprb) :: tr   (nlayers, size(prof))
  REAL   (KIND=jprb) :: tro  (nlayers, size(prof))
  REAL   (KIND=jprb) :: wr   (nlayers, size(prof))
  REAL   (KIND=jprb) :: or   (nlayers, size(prof))
  REAL   (KIND=jprb) :: co2r (nlayers, size(prof))
  REAL   (KIND=jprb) :: n2or (nlayers, size(prof))
  REAL   (KIND=jprb) :: cor  (nlayers, size(prof))
  REAL   (KIND=jprb) :: ch4r (nlayers, size(prof))
  REAL   (KIND=jprb) :: wwr  (nlayers, size(prof))
! user - reference
  REAL   (KIND=jprb) :: dt   (nlayers, size(prof))
  REAL   (KIND=jprb) :: dto  (nlayers, size(prof))
  REAL   (KIND=jprb) :: dtabs(nlayers, size(prof))
! pressure weighted
  REAL   (KIND=jprb) :: tw   (nlayers, size(prof))
  REAL   (KIND=jprb) :: twr  (nlayers, size(prof))
  REAL   (KIND=jprb) :: tuw  (nlayers, size(prof))
  REAL   (KIND=jprb) :: tuwr (nlayers, size(prof))
  REAL   (KIND=jprb) :: ww   (nlayers, size(prof))
  REAL   (KIND=jprb) :: ow   (nlayers, size(prof))
  REAL   (KIND=jprb) :: co2w (nlayers, size(prof))
  REAL   (KIND=jprb) :: n2ow (nlayers, size(prof))
  REAL   (KIND=jprb) :: cow  (nlayers, size(prof))
  REAL   (KIND=jprb) :: ch4w (nlayers, size(prof))
  REAL   (KIND=jprb) :: n2owr(nlayers, size(prof))
  REAL   (KIND=jprb) :: cowr (nlayers, size(prof))
  REAL   (KIND=jprb) :: ch4wr(nlayers, size(prof))
! intermediate variables
  REAL   (KIND=jprb) :: sum1, sum2
  REAL   (KIND=jprb) :: deltac    (nlayers                )
  REAL   (KIND=jprb) :: sum2_ww   (nlayers, size(prof)    )
  REAL   (KIND=jprb) :: sum2_wwr  (nlayers, size(prof)    )
  REAL   (KIND=jprb) :: sum2_ow   (nlayers, size(prof)    )
  REAL   (KIND=jprb) :: sum2_twr  (nlayers, size(prof)    )
  REAL   (KIND=jprb) :: sum2_tuw  (nlayers, size(prof)    )
  REAL   (KIND=jprb) :: sum2_co2w (nlayers, size(prof)    )
  REAL   (KIND=jprb) :: sum2_ch4w (nlayers, size(prof)    )
  REAL   (KIND=jprb) :: sum2_ch4wr(nlayers, size(prof)    )
  REAL   (KIND=jprb) :: sum2_n2ow (nlayers, size(prof)    )
  REAL   (KIND=jprb) :: sum2_n2owr(nlayers, size(prof)    )
  REAL   (KIND=jprb) :: sum2_cow  (nlayers, size(prof)    )
  REAL   (KIND=jprb) :: sum2_cowr (nlayers, size(prof)    )
  REAL   (KIND=jprb) :: tr_sq     (nlayers, size(prof)    )
  REAL   (KIND=jprb) :: tr_4      (nlayers, size(prof)    )
  REAL   (KIND=jprb) :: sec_wrwr  (nlayers, size(prof)    )
  REAL   (KIND=jprb) :: sec_wr    (nlayers, size(prof)    )
! K variables
  REAL   (KIND=jprb) :: t_k       (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: w_k       (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: o_k       (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: co2_k     (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: n2o_k     (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: co_k      (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: ch4_k     (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: tr_k      (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: tro_k     (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: wr_k      (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: or_k      (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: wwr_k     (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: co2r_k    (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: n2or_k    (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: cor_k     (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: ch4r_k    (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: twr_k     (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: tuw_k     (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: tuwr_k    (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: dt_k      (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: dto_k     (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: tw_k      (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: ww_k      (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: ow_k      (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: co2w_k    (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: n2ow_k    (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: cow_k     (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: ch4w_k    (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: n2owr_k   (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: cowr_k    (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: ch4wr_k   (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: sec_or_k  (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: sec_wr_k  (nlayers, size(chanprof))
  REAL   (KIND=jprb) :: sec_wrwr_k(nlayers, size(chanprof))

! John.Hague@ecmwf.int optimisation start
  Real(Kind=jprb) :: tmp1(nlayers,size(prof))
  Real(Kind=jprb) :: tmp2(nlayers,size(prof))
  Real(Kind=jprb) :: tmp3(nlayers,size(prof))
  Real(Kind=jprb) :: tmp4(nlayers,size(prof))
  Real(Kind=jprb) :: tmp5(nlayers,size(prof))
  Real(Kind=jprb) :: tmp6(nlayers,size(prof))
  Real(Kind=jprb) :: tmp7(nlayers,size(prof))
  Real(Kind=jprb) :: tmp8(nlayers,size(prof))
! optimisation end

  INTEGER(KIND=jpim) :: nprofiles                          ! Number of profiles
  INTEGER(KIND=jpim) :: nchannels                          ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_9_K', 0_jpim, ZHOOK_HANDLE)
  nprofiles = size(prof)
  nchannels = size(chanprof)
!-------------------------------------------------------------------------------
! Recompute Direct variables
!-------------------------------------------------------------------------------
! layer N-1 lies between levels N-1 and N

  DO j = 1, nprofiles
!-------------------------------------------------------------------------------
!1) Profile layer quantities
!-------------------------------------------------------------------------------

    DO layer = 1, prof(j)%nlayers
      level       = layer + 1
!-Temperature--------------------------------------------------------------------
      t(layer, j) = (prof(j)%t(level - 1) + prof(j)%t(level)) / 2._jprb
!-H2O----------------------------------------------------------------------------
      w(layer, j) = (prof(j)%q(level - 1) + prof(j)%q(level)) / 2._jprb
!-O3-----------------------------------------------------------------------------
      IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) &
        & o(layer, j) = (prof(j)%o3(level - 1) + prof(j)%o3(level)) / 2._jprb
!-CO2----------------------------------------------------------------------------

      IF (opts%rt_ir%co2_Data .AND. coef%nco2 > 0) THEN

        IF (opts%rt_ir%pc%addpc) THEN
          co2(layer, j) = (coef_pccomp%co2_pc_ref(level - 1) + coef_pccomp%co2_pc_ref(level)) * 0.5_jprb
        ELSE
          co2(layer, j) = (prof(j)%co2(level - 1) + prof(j)%co2(level)) * 0.5_jprb
        ENDIF

      ENDIF

!-N2O----------------------------------------------------------------------------

      IF (opts%rt_ir%n2o_Data .AND. coef%nn2o > 0) THEN

        IF (opts%rt_ir%pc%addpc) THEN
          n2o(layer, j) = (coef_pccomp%n2o_pc_ref(level - 1) + coef_pccomp%n2o_pc_ref(level)) * 0.5_jprb
        ELSE
          n2o(layer, j) = (prof(j)%n2o(level - 1) + prof(j)%n2o(level)) * 0.5_jprb
        ENDIF

      ENDIF

!-CO-----------------------------------------------------------------------------

      IF (opts%rt_ir%co_Data .AND. coef%nco > 0) THEN

        IF (opts%rt_ir%pc%addpc) THEN
          co(layer, j) = (coef_pccomp%co_pc_ref(level - 1) + coef_pccomp%co_pc_ref(level)) * 0.5_jprb
        ELSE
          co(layer, j) = (prof(j)%co(level - 1) + prof(j)%co(level)) * 0.5_jprb
        ENDIF

      ENDIF

!-CH4----------------------------------------------------------------------------

      IF (opts%rt_ir%ch4_Data .AND. coef%nch4 > 0) THEN

        IF (opts%rt_ir%pc%addpc) THEN
          ch4(layer, j) = (coef_pccomp%ch4_pc_ref(level - 1) + coef_pccomp%ch4_pc_ref(level)) * 0.5_jprb
        ELSE
          ch4(layer, j) = (prof(j)%ch4(level - 1) + prof(j)%ch4(level)) * 0.5_jprb
        ENDIF

      ENDIF

    ENDDO


!------------------------------------------------------------------------------
!2) calculate deviations from reference profile (layers)
!------------------------------------------------------------------------------
! All assignments 1:prof(1) % nlayers
    dt(:, j)      = t(:, j) - coef%tstar(:)
    dtabs(:, j)   = ABS(dt(:, j))
!------------------------------------------------------------------------------
!3) calculate (profile / reference profile) ratios; tr wr or
! if no input O3 profile, set to reference value (or =1)
!------------------------------------------------------------------------------
! All assignments 1:prof(1) % nlayers
    tr(:, j)      = t(:, j) / coef%tstar(:)
    tr_sq(:, j)   = tr(:, j) * tr(:, j)
    tr_4(:, j)    = tr_sq(:, j) * tr_sq(:, j)
    wr(:, j)      = w(:, j) / coef%wstar(:)

    IF (coef%nozone > 0) THEN
      dto(:, j) = t(:, j) - coef%to3star(:)
      tro(:, j) = t(:, j) / coef%to3star(:)
      IF (opts%rt_ir%ozone_Data) THEN
        or(:, j)  = o(:, j) / coef%ostar(:)
      ELSE
        or(:, j)  = 1._jprb
      ENDIF
    ENDIF

!-CO2----------------------------------------------------------------------------

    IF (opts%rt_ir%co2_Data .AND. coef%nco2 > 0) THEN
      co2r(:, j) = co2(:, j) / coef%co2star(:)
    ELSE
      co2r(:, j) = 1._jprb
    ENDIF

!-N2O----------------------------------------------------------------------------

    IF (opts%rt_ir%n2o_Data .AND. coef%nn2o > 0) THEN
      n2or(:, j) = n2o(:, j) / coef%n2ostar(:)
    ELSE
      n2or(:, j) = 1._jprb
    ENDIF

!-CO----------------------------------------------------------------------------

    IF (opts%rt_ir%co_Data .AND. coef%nco > 0) THEN
      cor(:, j) = co(:, j) / coef%costar(:)
    ELSE
      cor(:, j) = 1._jprb
    ENDIF

!-CH4---------------------------------------------------------------------------

    IF (opts%rt_ir%ch4_Data .AND. coef%nch4 > 0) THEN
      ch4r(:, j) = ch4(:, j) / coef%ch4star(:)
    ELSE
      ch4r(:, j) = 1._jprb
    ENDIF

!-------------------------------------------------------------------
! 4. calculate profile / reference profile sums: tw wwr
!--------------------------------------------------------------------
    tw(1, j) = 0._jprb

    DO layer = 2, prof(j)%nlayers
      tw(layer, j) = tw(layer - 1, j) + coef%dpp(layer - 1) * tr(layer - 1, j)
    ENDDO

    sum1 = 0._jprb
    sum2 = 0._jprb

    DO layer = 1, prof(j)%nlayers
      sum1 = sum1 + t(layer, j)
      sum2 = sum2 + coef%tstar(layer)
      sum2_tuw(layer, j) = sum2
      tuw(layer, j)      = sum1 / sum2
    ENDDO

    tuwr(1, j)                 = coef%dpp(0) * t(1, j) / (coef%dpp(0) * coef%tstar(1))
    tuwr(2:prof(j)%nlayers, j) = tuw(2:prof(j)%nlayers, j)

    IF (coef%nco2 > 0) THEN
      sum1 = 0._jprb
      sum2 = 0._jprb

      DO layer = 1, prof(j)%nlayers
        sum1 = sum1 + coef%dpp(layer - 1) * t(layer, j)
        sum2 = sum2 + coef%dpp(layer - 1) * coef%tstar(layer)
        sum2_twr(layer, j) = sum2
        twr(layer, j)      = sum1 / sum2
      ENDDO
    ENDIF

    sum1 = 0._jprb
    sum2 = 0._jprb

    DO layer = 1, prof(j)%nlayers
      sum1 = sum1 + coef%dpp(layer - 1) * w(layer, j) * t(layer, j)
      sum2 = sum2 + coef%dpp(layer - 1) * coef%wstar(layer) * coef%tstar(layer)
      sum2_wwr(layer, j) = sum2
      wwr(layer, j)      = sum1 / sum2
    ENDDO

    sum1 = 0._jprb
    sum2 = 0._jprb

    DO layer = 1, prof(j)%nlayers
      sum1 = sum1 + coef%dpp(layer - 1) * w(layer, j)
      sum2 = sum2 + coef%dpp(layer - 1) * coef%wstar(layer)
      sum2_ww(layer, j) = sum2
      ww(layer, j)      = sum1 / sum2
    ENDDO


    IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN
      sum1 = 0._jprb
      sum2 = 0._jprb

      DO layer = 1, prof(j)%nlayers
        sum1 = sum1 + coef%dpp(layer - 1) * o(layer, j)
        sum2 = sum2 + coef%dpp(layer - 1) * coef%ostar(layer)
        sum2_ow(layer, j) = sum2
        ow(layer, j)      = sum1 / sum2
      ENDDO

    ELSE
      sum2_ow(:, j) = 0._jprb
      ow(:, j)      = 1._jprb
    ENDIF


    IF (opts%rt_ir%co2_Data .AND. coef%nco2 > 0) THEN
      sum1 = 0._jprb
      sum2 = 0._jprb

      DO layer = 1, prof(j)%nlayers
        sum1 = sum1 + coef%dpp(layer - 1) * co2(layer, j)
        sum2 = sum2 + coef%dpp(layer - 1) * coef%co2star(layer)
        sum2_co2w(layer, j) = sum2
        co2w(layer, j)      = sum1 / sum2
      ENDDO

    ELSE
      sum2_co2w(:, j) = 0._jprb
      co2w(:, j)      = 1._jprb
    ENDIF

!-N2O---------------------------------------------------------------------------

    IF (coef%nn2o > 0) THEN
      IF (opts%rt_ir%n2o_Data) THEN
        sum1 = 0._jprb
        sum2 = 0._jprb
  
        DO layer = 1, prof(j)%nlayers
          sum1 = sum1 + coef%dpp(layer - 1) * n2o(layer, j)
          sum2 = sum2 + coef%dpp(layer - 1) * coef%n2ostar(layer)
          sum2_n2ow(layer, j) = sum2
          n2ow(layer, j)      = sum1 / sum2
        ENDDO
  
        sum1 = 0._jprb
        sum2 = 0._jprb
  
        DO layer = 1, prof(j)%nlayers
          sum1 = sum1 + coef%dpp(layer - 1) * n2o(layer, j) * t(layer, j)
          sum2 = sum2 + coef%dpp(layer - 1) * coef%n2ostar(layer) * coef%tstar(layer)
          sum2_n2owr(layer, j) = sum2
          n2owr(layer, j)      = sum1 / sum2
        ENDDO
  
      ELSE
        sum2_n2ow(:, j)  = 0._jprb
        n2ow(:, j)       = 1._jprb

        sum1 = 0._jprb
        sum2 = 0._jprb
  
        DO layer = 1, prof(j)%nlayers
          sum1 = sum1 + coef%dpp(layer - 1) * coef%n2ostar(layer) * t(layer, j)
          sum2 = sum2 + coef%dpp(layer - 1) * coef%n2ostar(layer) * coef%tstar(layer)
          sum2_n2owr(layer, j) = sum2
          n2owr(layer, j)      = sum1 / sum2
        ENDDO

      ENDIF
    ENDIF

!-CO---------------------------------------------------------------------------

    IF (coef%nco > 0) THEN
      IF (opts%rt_ir%co_Data) THEN
        sum1 = 0._jprb
        sum2 = 0._jprb
  
        DO layer = 1, prof(j)%nlayers
          sum1 = sum1 + coef%dpp(layer - 1) * co(layer, j)
          sum2 = sum2 + coef%dpp(layer - 1) * coef%costar(layer)
          sum2_cow(layer, j) = sum2
          cow(layer, j)      = sum1 / sum2
        ENDDO
  
        sum1 = 0._jprb
        sum2 = 0._jprb
  
        DO layer = 1, prof(j)%nlayers
          sum1 = sum1 + coef%dpp(layer - 1) * co(layer, j) * t(layer, j)
          sum2 = sum2 + coef%dpp(layer - 1) * coef%costar(layer) * coef%tstar(layer)
          sum2_cowr(layer, j) = sum2
          cowr(layer, j)      = sum1 / sum2
        ENDDO
  
      ELSE
        sum2_cow(:, j)  = 0._jprb
        cow(:, j)       = 1._jprb
        
        sum1 = 0._jprb
        sum2 = 0._jprb
  
        DO layer = 1, prof(j)%nlayers
          sum1 = sum1 + coef%dpp(layer - 1) * coef%costar(layer) * t(layer, j)
          sum2 = sum2 + coef%dpp(layer - 1) * coef%costar(layer) * coef%tstar(layer)
          sum2_cowr(layer, j) = sum2
          cowr(layer, j)      = sum1 / sum2
        ENDDO        
      ENDIF
    ENDIF
    
!-CH4---------------------------------------------------------------------------

    IF (coef%nch4 > 0) THEN
      IF (opts%rt_ir%ch4_Data) THEN
        sum1 = 0._jprb
        sum2 = 0._jprb
  
        DO layer = 1, prof(j)%nlayers
          sum1 = sum1 + coef%dpp(layer - 1) * ch4(layer, j)
          sum2 = sum2 + coef%dpp(layer - 1) * coef%ch4star(layer)
          sum2_ch4w(layer, j) = sum2
          ch4w(layer, j)      = sum1 / sum2
        ENDDO
  
        sum1 = 0._jprb
        sum2 = 0._jprb
  
        DO layer = 1, prof(j)%nlayers
          sum1 = sum1 + coef%dpp(layer - 1) * ch4(layer, j) * t(layer, j)
          sum2 = sum2 + coef%dpp(layer - 1) * coef%ch4star(layer) * coef%tstar(layer)
          sum2_ch4wr(layer, j) = sum2
          ch4wr(layer, j)      = sum1 / sum2
        ENDDO
  
      ELSE
        sum2_ch4w(:, j)  = 0._jprb
        ch4w(:, j)       = 1._jprb
        
        sum1 = 0._jprb
        sum2 = 0._jprb
  
        DO layer = 1, prof(j)%nlayers
          sum1 = sum1 + coef%dpp(layer - 1) * coef%ch4star(layer) * t(layer, j)
          sum2 = sum2 + coef%dpp(layer - 1) * coef%ch4star(layer) * coef%tstar(layer)
          sum2_ch4wr(layer, j) = sum2
          ch4wr(layer, j)      = sum1 / sum2
        ENDDO
      ENDIF
    ENDIF

  ENDDO
! Loop on profiles
!-------------------------------------------------------------------------
! K code
!-------------------------------------------------------------------------
  w_k(:,:)        = 0._jprb
  wr_k(:,:)       = 0._jprb
  ww_k(:,:)       = 0._jprb
  wwr_k(:,:)      = 0._jprb
  sec_wr_k(:,:)   = 0._jprb
  sec_wrwr_k(:,:) = 0._jprb
  dt_k(:,:)       = 0._jprb
  dto_k(:,:)      = 0._jprb
  t_k(:,:)        = 0._jprb
  tr_k(:,:)       = 0._jprb
  tro_k(:,:)      = 0._jprb
  tw_k(:,:)       = 0._jprb
  twr_k(:,:)      = 0._jprb
  o_k(:,:)        = 0._jprb
  or_k(:,:)       = 0._jprb
  ow_k(:,:)       = 0._jprb
  sec_or_k(:,:)   = 0._jprb
  tuw_k(:,:)      = 0._jprb
  tuwr_k(:,:)     = 0._jprb
  ch4_k(:,:)      = 0._jprb
  ch4r_k(:,:)     = 0._jprb
  ch4w_k(:,:)     = 0._jprb
  ch4wr_k(:,:)    = 0._jprb
  co_k(:,:)       = 0._jprb
  cor_k(:,:)      = 0._jprb
  cow_k(:,:)      = 0._jprb
  cowr_k(:,:)     = 0._jprb
  n2o_k(:,:)      = 0._jprb
  n2or_k(:,:)     = 0._jprb
  n2ow_k(:,:)     = 0._jprb
  n2owr_k(:,:)    = 0._jprb
  co2_k(:,:)      = 0._jprb
  co2r_k(:,:)     = 0._jprb
  co2w_k(:,:)     = 0._jprb
!-------------------------------------------------
!5.9 ch4            transmittance based on RTIASI
!-------------------------------------------------
!

  IF (coef%nch4 > 0) THEN
! John.Hague@ecmwf.int optimisation start
    Do j = 1, nprofiles
      Do layer = 1,prof_k(1) % nlayers
        tmp1(layer,j) = 1.5_jprb*predictors%ch4(2,layer,j)/ch4w(layer,j)
        tmp2(layer,j) = predictors%ch4(2,layer,j)**3_jpim/(ray_path(layer, j)*ch4w(layer,j)**2_jpim)
        tmp3(layer,j) = 0.5_jprb*ch4r(layer,j)**1.5_jprb/(ray_path(layer, j)**0.5_jprb*ch4w(layer,j))
        tmp4(layer,j) = 0.25_jprb / predictors%ch4(6,layer,j)**3_jpim
        tmp5(layer,j) = 1._jprb / predictors%ch4(2,layer,j)
      Enddo
    Enddo
! optimisation end

    DO i = 1, nchannels
      j = chanprof(i)%prof

      DO layer = 1, prof_k(i)%nlayers
        level = layer + 1
        ch4r_k(layer, i)               = ch4r_k(layer, i) + predictors_k%ch4(11, layer, i) * tmp1(layer,j)
        ch4w_k(layer, i)               = ch4w_k(layer, i) - predictors_k%ch4(11, layer, i) * tmp2(layer,j)
        ray_path_k(layer, i) = ray_path_k(layer, i) + predictors_k%ch4(11, layer, i) * tmp3(layer,j)
        ch4w_k(layer, i)               =      &
          & ch4w_k(layer, i) + ray_path(layer, j) * predictors_k%ch4(10, layer, i)
        ray_path_k(layer, i) =      &
          & ray_path_k(layer, i) + predictors_k%ch4(10, layer, i) * ch4w(layer, j)
        ch4w_k(layer, i)               = ch4w_k(layer, i) +      &
          & predictors_k%ch4(9, layer, i) * predictors%ch4(10, layer, j) * 2._jprb * ray_path(layer, j)
        ray_path_k(layer, i) = ray_path_k(layer, i) +      &
          & predictors_k%ch4(9, layer, i) * 2._jprb * predictors%ch4(10, layer, j) * ch4w(layer, j)
        ch4wr_k(layer, i)              = ch4wr_k(layer, i) + predictors_k%ch4(8, layer, i)
        ch4wr_k(layer, i)              =      &
          & ch4wr_k(layer, i) + predictors_k%ch4(7, layer, i) * ray_path(layer, j)
        ray_path_k(layer, i) =      &
          & ray_path_k(layer, i) + predictors_k%ch4(7, layer, i) * ch4wr(layer, j)
        ch4r_k(layer, i)               = ch4r_k(layer, i) + &
          & predictors_k%ch4(6, layer, i) * ray_path(layer, j) * tmp4(layer, j)
        ray_path_k(layer, i) = ray_path_k(layer, i) +      &
          & predictors_k%ch4(6, layer, i) * ch4r(layer, j) * tmp4(layer,j)
        dt_k(layer, i)                 = dt_k(layer, i) + predictors_k%ch4(5, layer, i) * ch4r(layer, j)
        ch4r_k(layer, i)               = ch4r_k(layer, i) + predictors_k%ch4(5, layer, i) * dt(layer, j)
        ch4r_k(layer, i)               = ch4r_k(layer, i) +      &
          & predictors_k%ch4(4, layer, i) * predictors%ch4(1, layer, j) * ray_path(layer, j) * 2._jprb
        ray_path_k(layer, i) = ray_path_k(layer, i) +      &
          & predictors_k%ch4(4, layer, i) * 2 * predictors%ch4(1, layer, j) * ch4r(layer, j)
        dt_k(layer, i)                 = dt_k(layer, i) + predictors_k%ch4(3, layer, i) * predictors%ch4(1, layer, j)
        ch4r_k(layer, i)               =      &
          & ch4r_k(layer, i) + predictors_k%ch4(3, layer, i) * ray_path(layer, j) * dt(layer, j)
        ray_path_k(layer, i) =      &
          & ray_path_k(layer, i) + predictors_k%ch4(3, layer, i) * ch4r(layer, j) * dt(layer, j)
        ch4r_k(layer, i)               = ch4r_k(layer, i) +      &
          & predictors_k%ch4(2, layer, i) * 0.5_jprb * ray_path(layer, j) * tmp5(layer,j)
        ray_path_k(layer, i) = ray_path_k(layer, i) +      &
          & predictors_k%ch4(2, layer, i) * 0.5_jprb * ch4r(layer, j) * tmp5(layer,j)
        ch4r_k(layer, i)               = ch4r_k(layer, i) + predictors_k%ch4(1, layer, i) * ray_path(layer, j)
        ray_path_k(layer, i) = ray_path_k(layer, i) + predictors_k%ch4(1, layer, i) * ch4r(layer, j)
      ENDDO

    ENDDO

  ENDIF

!-------------------------------------------------
!5.8 CO             transmittance based on RTIASI
!-------------------------------------------------
!

  IF (coef%nco > 0) THEN
! John.Hague@ecmwf.int optimisation start
    Do j = 1, nprofiles
      Do layer = 1,prof_k(1) % nlayers
        tmp1(layer,j) = 0.25_jprb*ray_path(layer, j)*(ray_path(layer, j)*cowr(layer,j))**(-0.75_jprb)
        tmp2(layer,j) = 0.25_jprb*cowr(layer,j)*(ray_path(layer, j)*cowr(layer,j))**(-0.75_jprb)
        tmp3(layer,j) = ray_path(layer, j)*0.4_jprb*(ray_path(layer, j)*cowr(layer,j))**(-0.6_jprb)
        tmp4(layer,j) = cowr(layer,j)*0.4_jprb* (ray_path(layer, j)*cowr(layer,j))**(-0.6_jprb)
        tmp5(layer,j) = 2._jprb*predictors%co(1,layer,j)/cow(layer,j)**0.25_jprb
        tmp6(layer,j) = cor(layer,j)**2_jpim  /cow(layer,j)**0.25_jprb
        tmp7(layer,j) = 2._jprb*predictors%co(1,layer,j)/cow(layer,j)**0.5_jprb
        tmp8(layer,j) = 0.5_jprb*cor(layer,j)**1.5_jprb/(ray_path(layer, j)**0.5_jprb*cow(layer,j))
      Enddo
    Enddo
! optimisation end

    DO i = 1, nchannels
      j = chanprof(i)%prof

      DO layer = 1, prof_k(i)%nlayers
        level = layer + 1
        cowr_k(layer, i)               = cowr_k(layer, i) + predictors_k%co(13, layer, i) * tmp1(layer,j)
        ray_path_k(layer, i) = ray_path_k(layer, i) +  predictors_k%co(13, layer, i) * tmp2(layer,j)
        cowr_k(layer, i)               = cowr_k(layer, i) +                       &
          & predictors_k%co(12, layer, i) * tmp3(layer,j)
        ray_path_k(layer, i) = ray_path_k(layer, i) +      &
          & predictors_k%co(12, layer, i) * tmp4(layer,j)
        cor_k(layer, i)                =      &
          & cor_k(layer, i) + predictors_k%co(11, layer, i) * tmp5(layer,j)
        cow_k(layer, i)                =      &
          & cow_k(layer, i) - predictors_k%co(11, layer, i) * 0.25_jprb * predictors%co(11, layer, j) / cow(layer, j)
        ray_path_k(layer, i) =      &
          & ray_path_k(layer, i) + predictors_k%co(11, layer, i) * tmp6(layer,j)
        cor_k(layer, i)                =      &
          & cor_k(layer, i) + predictors_k%co(10, layer, i) * tmp7(layer,j)
        cow_k(layer, i)                =      &
          & cow_k(layer, i) - predictors_k%co(10, layer, i) * 0.5_jprb * predictors%co(10, layer, j) / cow(layer, j)
        ray_path_k(layer, i) =      &
          & ray_path_k(layer, i) + predictors_k%co(10, layer, i) * cor(layer, j) ** 2_jpim / SQRT(cow(layer, j))
        cor_k(layer, i)                =      &
          & cor_k(layer, i) + predictors_k%co(9, layer, i) * 1.5_jprb * predictors%co(2, layer, j) / cow(layer, j)
        cow_k(layer, i)                = cow_k(layer, i) -                    &
          & predictors_k%co(9, layer, i) * predictors%co(2, layer, j) ** 3_jpim /  &
          & (ray_path(layer, j) * cow(layer, j) ** 2_jpim)
        ray_path_k(layer, i) = ray_path_k(layer, i) +      &
          & predictors_k%co(9, layer, i) * tmp8(layer,j)
        cor_k(layer, i)                =      &
          & cor_k(layer, i) + predictors_k%co(8, layer, i) * 2._jprb * predictors%co(1, layer, j) / cow(layer, j)
        cow_k(layer, i)                =      &
          & cow_k(layer, i) - predictors_k%co(8, layer, i) * predictors%co(8, layer, j) / cow(layer, j)
        ray_path_k(layer, i) =      &
          & ray_path_k(layer, i) + predictors_k%co(8, layer, i) * cor(layer, j) ** 2_jpim / cow(layer, j)
        cor_k(layer, i)                =      &
          & cor_k(layer, i) + predictors_k%co(7, layer, i) * dtabs(layer, j) * ray_path(layer, j) * dt(layer, j)
        dt_k(layer, i)                 =      &
          & dt_k(layer, i) + predictors_k%co(7, layer, i) * dtabs(layer, j) * predictors%co(1, layer, j) * 2._jprb
        ray_path_k(layer, i) =      &
          & ray_path_k(layer, i) + predictors_k%co(7, layer, i) * dtabs(layer, j) * cor(layer, j) * dt(layer, j)
        cor_k(layer, i)                = cor_k(layer, i) +      &
          & predictors_k%co(6, layer, i) * ray_path(layer, j) * 0.25_jprb / predictors%co(6, layer, j) ** 3_jpim
        ray_path_k(layer, i) = ray_path_k(layer, i) +      &
          & predictors_k%co(6, layer, i) * cor(layer, j) * 0.25_jprb / predictors%co(6, layer, j) ** 3_jpim
        dt_k(layer, i)                 = dt_k(layer, i) + predictors_k%co(5, layer, i) * predictors%co(2, layer, j)
        cor_k(layer, i)                = cor_k(layer, i) +      &
          & predictors_k%co(5, layer, i) * 0.5_jprb * dt(layer, j) * ray_path(layer, j) / predictors%co(2, layer, j)
        ray_path_k(layer, i) = ray_path_k(layer, i) +      &
          & predictors_k%co(5, layer, i) * 0.5_jprb * dt(layer, j) * cor(layer, j) / predictors%co(2, layer, j)
        cor_k(layer, i)                =      &
          & cor_k(layer, i) + predictors_k%co(4, layer, i) * 2._jprb * predictors%co(1, layer, j) * ray_path(layer, j)
        ray_path_k(layer, i) =      &
          & ray_path_k(layer, i) + predictors_k%co(4, layer, i) * 2._jprb * predictors%co(1, layer, j) * cor(layer, j)
        dt_k(layer, i)                 = dt_k(layer, i) + predictors_k%co(3, layer, i) * predictors%co(1, layer, j)
        cor_k(layer, i)                =      &
          & cor_k(layer, i) + predictors_k%co(3, layer, i) * dt(layer, j) * ray_path(layer, j)
        ray_path_k(layer, i) =      &
          & ray_path_k(layer, i) + predictors_k%co(3, layer, i) * cor(layer, j) * dt(layer, j)
        cor_k(layer, i)                =      &
          & cor_k(layer, i) + predictors_k%co(2, layer, i) * 0.5_jprb * ray_path(layer, j) / predictors%co(2, layer, j)
        ray_path_k(layer, i) =      &
          & ray_path_k(layer, i) + predictors_k%co(2, layer, i) * 0.5_jprb * cor(layer, j) / predictors%co(2, layer, j)
        cor_k(layer, i)                = cor_k(layer, i) + predictors_k%co(1, layer, i) * ray_path(layer, j)
        ray_path_k(layer, i) = ray_path_k(layer, i) + predictors_k%co(1, layer, i) * cor(layer, j)
      ENDDO

    ENDDO

  ENDIF

!-------------------------------------------------
!5.7 N2O            transmittance based on RTIASI
!-------------------------------------------------
!

  IF (coef%nn2o > 0) THEN
! John.Hague@ecmwf.int optimisation start
    Do j = 1, nprofiles
      Do layer = 1,prof_k(1) % nlayers
        tmp1(layer,j)=0.5_jprb*n2or(layer,j)**1.5_jprb/(ray_path(layer, j)**0.5_jprb*n2ow(layer,j))
      Enddo
    Enddo
! optimisation end

    DO i = 1, nchannels
      j = chanprof(i)%prof

      DO layer = 1, prof_k(i)%nlayers
        level = layer + 1
        dt_k(layer, i)                 =      &
          & dt_k(layer, i) + predictors_k%n2o(13, layer, i) * ray_path(layer, j) ** 2_jpim * n2owr(layer, j)
        n2owr_k(layer, i)              =      &
          & n2owr_k(layer, i) + predictors_k%n2o(13, layer, i) * ray_path(layer, j) ** 2_jpim * dt(layer, j)
        ray_path_k(layer, i) = ray_path_k(layer, i) +      &
          & predictors_k%n2o(13, layer, i) * 2._jprb * ray_path(layer, j) * n2owr(layer, j) * dt(layer, j)
        n2owr_k(layer, i)              = n2owr_k(layer, i) +      &
          & predictors_k%n2o(12, layer, i) * 3._jprb * predictors%n2o(8, layer, j) ** 2_jpim * ray_path(layer, j)
        ray_path_k(layer, i) = ray_path_k(layer, i) +      &
          & predictors_k%n2o(12, layer, i) * 3._jprb * predictors%n2o(8, layer, j) ** 2_jpim * n2owr(layer, j)
        n2owr_k(layer, i)              = n2owr_k(layer, i) +      &
          & predictors_k%n2o(11, layer, i) * 2._jprb * predictors%n2o(8, layer, j) * ray_path(layer, j)
        ray_path_k(layer, i) = ray_path_k(layer, i) +      &
          & predictors_k%n2o(11, layer, i) * predictors%n2o(8, layer, j) * 2._jprb * n2owr(layer, j)
        n2or_k(layer, i)               =      &
          & n2or_k(layer, i) + predictors_k%n2o(10, layer, i) * 1.5_jprb * predictors%n2o(2, layer, j) / n2ow(layer, j)
        n2ow_k(layer, i)               = n2ow_k(layer, i) -                      &
          & predictors_k%n2o(10, layer, i) * predictors%n2o(2, layer, j) ** 3_jpim /  &
          & (ray_path(layer, j) * n2ow(layer, j) ** 2_jpim)
        ray_path_k(layer, i) = ray_path_k(layer, i) +      &
          & predictors_k%n2o(10, layer, i) * tmp1(layer,j)
        n2owr_k(layer, i)              = n2owr_k(layer, i) + predictors_k%n2o(9, layer, i)
        n2owr_k(layer, i)              =      &
          & n2owr_k(layer, i) + predictors_k%n2o(8, layer, i) * ray_path(layer, j)
        ray_path_k(layer, i) =      &
          & ray_path_k(layer, i) + predictors_k%n2o(8, layer, i) * n2owr(layer, j)
        n2ow_k(layer, i)               = n2ow_k(layer, i) + predictors_k%n2o(7, layer, i) * ray_path(layer, j)
        ray_path_k(layer, i) = ray_path_k(layer, i) + predictors_k%n2o(7, layer, i) * n2ow(layer, j)
        n2or_k(layer, i)               = n2or_k(layer, i) +      &
          & predictors_k%n2o(6, layer, i) * ray_path(layer, j) * 0.25_jprb / predictors%n2o(6, layer, j) ** 3_jpim
        ray_path_k(layer, i) = ray_path_k(layer, i) +      &
          & predictors_k%n2o(6, layer, i) * n2or(layer, j) * 0.25_jprb / predictors%n2o(6, layer, j) ** 3_jpim
        dt_k(layer, i)                 = dt_k(layer, i) + predictors_k%n2o(5, layer, i) * n2or(layer, j)
        n2or_k(layer, i)               = n2or_k(layer, i) + predictors_k%n2o(5, layer, i) * dt(layer, j)
        n2or_k(layer, i)               = n2or_k(layer, i) +      &
          & predictors_k%n2o(4, layer, i) * 2._jprb * predictors%n2o(1, layer, j) * ray_path(layer, j)
        ray_path_k(layer, i) = ray_path_k(layer, i) +      &
          & predictors_k%n2o(4, layer, i) * 2._jprb * predictors%n2o(1, layer, j) * n2or(layer, j)
        dt_k(layer, i)                 = dt_k(layer, i) + predictors_k%n2o(3, layer, i) * predictors%n2o(1, layer, j)
        n2or_k(layer, i)               =      &
          & n2or_k(layer, i) + predictors_k%n2o(3, layer, i) * ray_path(layer, j) * dt(layer, j)
        ray_path_k(layer, i) =      &
          & ray_path_k(layer, i) + predictors_k%n2o(3, layer, i) * n2or(layer, j) * dt(layer, j)
        n2or_k(layer, i)               = n2or_k(layer, i) +      &
          & predictors_k%n2o(2, layer, i) * 0.5_jprb * ray_path(layer, j) / predictors%n2o(2, layer, j)
        ray_path_k(layer, i) = ray_path_k(layer, i) +      &
          & predictors_k%n2o(2, layer, i) * 0.5_jprb * n2or(layer, j) / predictors%n2o(2, layer, j)
        n2or_k(layer, i)               = n2or_k(layer, i) + predictors_k%n2o(1, layer, i) * ray_path(layer, j)
        ray_path_k(layer, i) = ray_path_k(layer, i) + predictors_k%n2o(1, layer, i) * n2or(layer, j)
      ENDDO

    ENDDO

  ENDIF

!-----------------------------------------------------------------------------------------------
!5.6 CO2 transmittance
!-----------------------------------------------------------------------------------------------

  IF (coef%nco2 > 0) THEN
! John.Hague@ecmwf.int optimisation start
    Do j = 1, nprofiles
      Do layer = 1,prof_k(1) % nlayers
        tmp1(layer,j)=twr(layer,j)**2_jpim*tr(layer,j)**2_jpim*ray_path(layer, j)**0.5_jprb
        tmp2(layer,j)=twr(layer,j)*tr(layer,j)**0.5_jprb
      Enddo
    Enddo
! optimisation end

    DO i = 1, nchannels
      j = chanprof(i)%prof

      DO layer = 1, prof_k(i)%nlayers
        level = layer + 1
        tr_k(layer, i)                 =      &
          & tr_k(layer, i) + predictors_k%co2(15, layer, i) * 2._jprb * SQRT(predictors%co2(15, layer, j)) * twr(layer, j)
        twr_k(layer, i)                =      &
          & twr_k(layer, i) + predictors_k%co2(15, layer, i) * 2._jprb * SQRT(predictors%co2(15, layer, j)) * tr(layer, j)
        twr_k(layer, i)                = twr_k(layer, i) +                                 &
          & predictors_k%co2(14, layer, i) * 3._jprb * tmp1(layer,j)
        tr_k(layer, i)                 = tr_k(layer, i) +      &
          & predictors_k%co2(14, layer, i) * 2._jprb * tr(layer, j) * twr(layer, j) ** 3_jpim * SQRT(ray_path(layer, j))
        ray_path_k(layer, i) = ray_path_k(layer, i) +                    &
          & predictors_k%co2(14, layer, i) * tr(layer, j) ** 2_jpim * twr(layer, j) ** 3_jpim * 0.5_jprb /  &
          & SQRT(ray_path(layer, j))
        tr_k(layer, i)                 =      &
          & tr_k(layer, i) + predictors_k%co2(13, layer, i) * 3._jprb * tr(layer, j) ** 2_jpim * ray_path(layer, j)
        ray_path_k(layer, i) =      &
          & ray_path_k(layer, i) + predictors_k%co2(13, layer, i) * tr(layer, j) ** 3_jpim
        tr_k(layer, i)                 = tr_k(layer, i) + predictors_k%co2(12, layer, i) * 3 * tr(layer, j) ** 2_jpim
        co2r_k(layer, i)               = co2r_k(layer, i) +      &
          & predictors_k%co2(11, layer, i) * 0.5_jprb * ray_path(layer, j) / SQRT(predictors%co2(1, layer, j))
        ray_path_k(layer, i) = ray_path_k(layer, i) +      &
          & predictors_k%co2(11, layer, i) * 0.5_jprb * co2r(layer, j) / SQRT(predictors%co2(1, layer, j))
        twr_k(layer, i)                = twr_k(layer, i) +      &
          & predictors_k%co2(10, layer, i) * ray_path(layer, j) * SQRT(predictors%co2(5, layer, j))
        tr_k(layer, i)                 = tr_k(layer, i) +      &
          & predictors_k%co2(10, layer, i) * 0.5_jprb * predictors%co2(7, layer, j) / SQRT(predictors%co2(5, layer, j))
        ray_path_k(layer, i) =      &
          & ray_path_k(layer, i) + predictors_k%co2(10, layer, i) * tmp2(layer,j)
        twr_k(layer, i)                = twr_k(layer, i) + predictors_k%co2(9, layer, i) * 3._jprb * twr(layer, j) ** 2_jpim
        co2w_k(layer, i)               = co2w_k(layer, i) +      &
          & predictors_k%co2(8, layer, i) * 2._jprb * ray_path(layer, j) * SQRT(predictors%co2(8, layer, j))
        ray_path_k(layer, i) = ray_path_k(layer, i) +      &
          & predictors_k%co2(8, layer, i) * 2._jprb * co2w(layer, j) * SQRT(predictors%co2(8, layer, j))
        twr_k(layer, i)                = twr_k(layer, i) + predictors_k%co2(7, layer, i) * ray_path(layer, j)
        ray_path_k(layer, i) = ray_path_k(layer, i) + predictors_k%co2(7, layer, i) * twr(layer, j)
        ray_path_k(layer, i) = ray_path_k(layer, i) + predictors_k%co2(6, layer, i)
        tr_k(layer, i)                 = tr_k(layer, i) + predictors_k%co2(5, layer, i)
        tr_k(layer, i)                 =      &
          & tr_k(layer, i) + predictors_k%co2(4, layer, i) * 2._jprb * predictors%co2(3, layer, j)
        ray_path_k(layer, i) =      &
          & ray_path_k(layer, i) + predictors_k%co2(4, layer, i) * tr(layer, j) ** 2_jpim
        tr_k(layer, i)                 = tr_k(layer, i) + predictors_k%co2(3, layer, i) * ray_path(layer, j)
        ray_path_k(layer, i) = ray_path_k(layer, i) + predictors_k%co2(3, layer, i) * tr(layer, j)
        tr_k(layer, i)                 =      &
          & tr_k(layer, i) + predictors_k%co2(2, layer, i) * 2._jprb * predictors%co2(5, layer, j)
        co2r_k(layer, i)               = co2r_k(layer, i) + predictors_k%co2(1, layer, i) * ray_path(layer, j)
        ray_path_k(layer, i) = ray_path_k(layer, i) + predictors_k%co2(1, layer, i) * co2r(layer, j)
      ENDDO

    ENDDO

  ENDIF

!-------------------------------------------------------------------------------------
!5.5 cloud
!-------------------------------------------------------------------------------------

  IF (coef%id_sensor == sensor_id_mw) THEN

    IF (opts%rt_mw%clw_Data) THEN

      DO layer = 1, prof_k(1)%nlayers
        deltac(layer) = 0.1820_jprb * 100.0_jprb * coef%dp(layer) / (4.3429_jprb * gravity)
      ENDDO


      DO i = 1, nchannels
        j = chanprof(i)%prof

        DO layer = 2, prof_k(i)%nlayers
          level = layer + 1
          prof_k(i)%clw(level - 1)   =      &
            & prof_k(i)%clw(level - 1) + 0.5_jprb * predictors_k%clw(layer, i) * deltac(layer) * geom(j)%seczen
          predictors_k%clw(layer, i) = 0.5_jprb * predictors_k%clw(layer, i)
        ENDDO

        DO layer = 1, prof_k(i)%nlayers
          level                = layer + 1
          prof_k(i)%clw(level) = prof_k(i)%clw(level) + predictors_k%clw(layer, i) * deltac(layer) * geom(j)%seczen
        ENDDO

      ENDDO

    ENDIF

  ENDIF

!------------------------------------------------------------------------------------
!5.4 ozone
!------------------------------------------------------------------------------------

  IF (coef%nozone > 0) THEN
! John.Hague@ecmwf.int optimisation start
    Do j = 1, nprofiles
      Do layer = 1,prof_k(1) % nlayers
        tmp1(layer,j)=0.5_jprb*ray_path(layer, j)**(-0.5_jprb)* ow(layer,j)**2_jpim*dto(layer,j)
        tmp2(layer,j)=2._jprb*ow(layer,j)*ray_path(layer, j)**0.5_jprb *dto(layer,j)
        tmp3(layer,j)=ray_path(layer, j)**0.5_jprb *ow(layer,j)**2_jpim
        tmp4(layer,j)=1.75_jprb*ray_path(layer, j)*(ray_path(layer, j)*ow(layer,j))**0.75_jprb
        tmp5(layer,j)=1.75_jprb*ow(layer,j)*(ray_path(layer, j)*ow(layer,j))**0.75_jprb
        tmp6(layer,j)=0.5_jprb*or(layer,j)**1.5_jprb / (ray_path(layer, j)**0.5_jprb*ow(layer,j))
      Enddo
    Enddo
! optimisation end

    DO i = 1, nchannels
      j = chanprof(i)%prof

      DO layer = 1, prof_k(i)%nlayers
        level = layer + 1
! One can pack all ow_k lines in one longer statement
! same for sec_or_k and dto_k
        ray_path_k(layer, i) =      &
          & ray_path_k(layer, i) + predictors_k%ozone(15, layer, i) * tro(layer, j) ** 3_jpim
        tro_k(layer, i)                =      &
          & tro_k(layer, i) + predictors_k%ozone(15, layer, i) * 3._jprb * tro(layer, j) ** 2_jpim * ray_path(layer, j)
        ray_path_k(layer, i) = ray_path_k(layer, i) +                                            &
          & predictors_k%ozone(14, layer, i) * tmp1(layer,j)
        ow_k(layer, i)                 = ow_k(layer, i) +      &
          & predictors_k%ozone(14, layer, i) * tmp2(layer,j)
        dto_k(layer, i)                =      &
          & dto_k(layer, i) + predictors_k%ozone(14, layer, i) * tmp3(layer,j)
        ow_k(layer, i)                 = ow_k(layer, i) +                             &
          & predictors_k%ozone(13, layer, i) * tmp4(layer,j)
        ray_path_k(layer, i) = ray_path_k(layer, i) +      &
          & predictors_k%ozone(13, layer, i) * tmp5(layer,j)
        ray_path_k(layer, i) =      &
          & ray_path_k(layer, i) + predictors_k%ozone(12, layer, i) * or(layer, j) / ow(layer, j)
        or_k(layer, i)                 =      &
          & or_k(layer, i) + predictors_k%ozone(12, layer, i) * ray_path(layer, j) / ow(layer, j)
        ow_k(layer, i)                 = ow_k(layer, i) -      &
          & predictors_k%ozone(12, layer, i) * or(layer, j) * ray_path(layer, j) / ow(layer, j) ** 2_jpim
        ow_k(layer, i)                 = ow_k(layer, i) +      &
          & predictors_k%ozone(11, layer, i) * 2._jprb * ray_path(layer, j) * predictors%ozone(10, layer, j)
        ray_path_k(layer, i) = ray_path_k(layer, i) +      &
          & predictors_k%ozone(11, layer, i) * 2._jprb * ow(layer, j) * predictors%ozone(10, layer, j)
        ow_k(layer, i)                 =      &
          & ow_k(layer, i) + predictors_k%ozone(10, layer, i) * ray_path(layer, j)
        ray_path_k(layer, i) =      &
          & ray_path_k(layer, i) + predictors_k%ozone(10, layer, i) * ow(layer, j)
        or_k(layer, i)                 = or_k(layer, i) +                                          &
          & predictors_k%ozone(9, layer, i) * SQRT(ray_path(layer, j) * ow(layer, j)) *  &
          & ray_path(layer, j)
        ow_k(layer, i)                 = ow_k(layer, i) +                                                           &
          & predictors_k%ozone(9, layer, i) * predictors%ozone(1, layer, j) * 0.5_jprb * ray_path(layer, j) /  &
          & SQRT(ray_path(layer, j) * ow(layer, j))
        ray_path_k(layer, i) = ray_path_k(layer, i) +      &
          & predictors_k%ozone(9, layer, i) * 1.5_jprb * or(layer, j) * SQRT(predictors%ozone(10, layer, j))
        or_k(layer, i)                 =      &
          & or_k(layer, i) + predictors_k%ozone(8, layer, i) * predictors%ozone(10, layer, j)
        ow_k(layer, i)                 =      &
          & ow_k(layer, i) + predictors_k%ozone(8, layer, i) * ray_path(layer, j) * or(layer, j)
        ray_path_k(layer, i) =      &
          & ray_path_k(layer, i) + predictors_k%ozone(8, layer, i) * or(layer, j) * ow(layer, j)
        or_k(layer, i)                 =      &
          & or_k(layer, i) + predictors_k%ozone(7, layer, i) * 1.5_jprb * predictors%ozone(2, layer, j) / ow(layer, j)
        ow_k(layer, i)                 =      &
          & ow_k(layer, i) - predictors_k%ozone(7, layer, i) * predictors%ozone(7, layer, j) / ow(layer, j)
        ray_path_k(layer, i) = ray_path_k(layer, i) +      &
          & predictors_k%ozone(7, layer, i) * tmp6(layer,j)
        or_k(layer, i)                 =      &
          & or_k(layer, i) + predictors_k%ozone(6, layer, i) * 2._jprb * predictors%ozone(1, layer, j) * ow(layer, j)
        ow_k(layer, i)                 =      &
          & ow_k(layer, i) + predictors_k%ozone(6, layer, i) * predictors%ozone(4, layer, j) / ray_path(layer, j)
        ray_path_k(layer, i) =      &
          & ray_path_k(layer, i) + predictors_k%ozone(6, layer, i) * or(layer, j) ** 2_jpim * ow(layer, j)
        sec_or_k(layer, i)             = sec_or_k(layer, i) +                             &
          & predictors_k%ozone(5, layer, i) * 0.5_jprb * predictors%ozone(3, layer, j) /  &
          & (predictors%ozone(1, layer, j) * predictors%ozone(2, layer, j))
        dto_k(layer, i)                =      &
          & dto_k(layer, i) + predictors_k%ozone(5, layer, i) * predictors%ozone(2, layer, j)
        sec_or_k(layer, i)             =      &
          & sec_or_k(layer, i) + predictors_k%ozone(4, layer, i) * 2._jprb * predictors%ozone(1, layer, j)
        sec_or_k(layer, i)             = sec_or_k(layer, i) +      &
          & predictors_k%ozone(3, layer, i) * predictors%ozone(3, layer, j) / predictors%ozone(1, layer, j)
        dto_k(layer, i)                =      &
          & dto_k(layer, i) + predictors_k%ozone(3, layer, i) * predictors%ozone(1, layer, j)
        sec_or_k(layer, i)             =      &
          & sec_or_k(layer, i) + predictors_k%ozone(2, layer, i) * 0.5_jprb / predictors%ozone(2, layer, j)
        sec_or_k(layer, i)             = sec_or_k(layer, i) + predictors_k%ozone(1, layer, i)
        or_k(layer, i)                 = or_k(layer, i) + sec_or_k(layer, i) * ray_path(layer, j)
        ray_path_k(layer, i) = ray_path_k(layer, i) + sec_or_k(layer, i) * or(layer, j)
      ENDDO

    ENDDO

  ENDIF

!--------------------------------------------
!5.3 Water Vapour Continuum based on RTIASI
!--------------------------------------------

  IF (coef%nwvcont > 0) THEN

    DO i = 1, nchannels
      j = chanprof(i)%prof

      DO layer = 1, prof_k(i)%nlayers
        level                = layer + 1
        sec_wr_k(layer, i)   = sec_wr_k(layer, i) + predictors_k%wvcont(4, layer, i) / tr_sq(layer, j)
        tr_k(layer, i)       = tr_k(layer, i) -      &
          & 2 * predictors_k%wvcont(4, layer, i) * predictors%watervapour(7, layer, j) / (tr_sq(layer, j) * tr(layer, j))
        sec_wr_k(layer, i)   = sec_wr_k(layer, i) + predictors_k%wvcont(3, layer, i) / tr(layer, j)
        tr_k(layer, i)       =      &
          & tr_k(layer, i) - predictors_k%wvcont(3, layer, i) * predictors%watervapour(7, layer, j) / tr_sq(layer, j)
        sec_wrwr_k(layer, i) = sec_wrwr_k(layer, i) + predictors_k%wvcont(2, layer, i) / tr_4(layer, j)
        tr_k(layer, i)       =      &
          & tr_k(layer, i) - 4._jprb * predictors_k%wvcont(2, layer, i) * predictors%wvcont(1, layer, j) / tr_4(layer, j)
        sec_wrwr_k(layer, i) = sec_wrwr_k(layer, i) + predictors_k%wvcont(1, layer, i) / tr(layer, j)
        tr_k(layer, i)       =      &
          & tr_k(layer, i) - predictors_k%wvcont(1, layer, i) * predictors%wvcont(1, layer, j) / tr(layer, j)
      ENDDO

    ENDDO

  ENDIF

!---------------------------------
!5.2 water vapour based on RTIASI
!---------------------------------
! John.Hague@ecmwf.int optimisation start
  Do j = 1, nprofiles
    Do layer = 1,prof_k(1) % nlayers
      tmp1(layer,j)=predictors% watervapour(7,layer,j)**1.5_jprb
      tmp2(layer,j)=0.5_jprb*wr(layer,j)**1.5_jprb/ray_path(layer, j)**0.5_jprb
      tmp3(layer,j)=1.25_jprb*ww(layer,j)* (predictors % watervapour(2,layer,j))**0.25_jprb
      tmp4(layer,j)=0.5_jprb/SQRT(ray_path(layer, j))*(wr(layer,j)**1.5_jprb / wwr(layer,j))
      tmp5(layer,j)=SQRT(ray_path(layer, j))*1.5_jprb*wr(layer,j)**0.5_jprb/ wwr(layer,j)
      tmp6(layer,j)=SQRT(ray_path(layer, j))*wr(layer,j)**1.5_jprb / wwr(layer,j)**2_jpim
    Enddo
  Enddo
! optimisation end

  DO i = 1, nchannels
    j = chanprof(i)%prof

    DO layer = 1, prof_k(i)%nlayers
      level = layer + 1
      sec_wr(layer, j)               = ray_path(layer, j) * wr(layer, j)
      sec_wrwr(layer, j)             = sec_wr(layer, j) * wr(layer, j)
      wr_k(layer, i)                 =      &
        & wr_k(layer, i) + predictors_k%watervapour(19, layer, i) * 2._jprb * predictors%watervapour(7, layer, j) / ww(layer, j)
      ww_k(layer, i)                 = ww_k(layer, i) -                                   &
        & predictors_k%watervapour(19, layer, i) * predictors%watervapour(1, layer, j) /  &
        & (ray_path(layer, j) * ww(layer, j) ** 2_jpim)
      ray_path_k(layer, i) =      &
        & ray_path_k(layer, i) + predictors_k%watervapour(19, layer, i) * wr(layer, j) ** 2_jpim / ww(layer, j)
      dt_k(layer, i)                 =      &
        & dt_k(layer, i) + predictors_k%watervapour(18, layer, i) * tmp1(layer,j)
      wr_k(layer, i)                 = wr_k(layer, i) +                                                                     &
        & predictors_k%watervapour(18, layer, i) * 1.5_jprb * predictors%watervapour(5, layer, j) * ray_path(layer, j) &
        &  * dt(layer, j)
      ray_path_k(layer, i) = ray_path_k(layer, i) +      &
        & predictors_k%watervapour(18, layer, i) * 1.5_jprb * predictors%watervapour(5, layer, j) * wr(layer, j) * dt(layer, j)
      wr_k(layer, i)                 =      &
        & wr_k(layer, i) + predictors_k%watervapour(17, layer, i) * 1.5_jprb * predictors%watervapour(5, layer, j)
      ray_path_k(layer, i) = ray_path_k(layer, i) +      &
        & predictors_k%watervapour(17, layer, i) * tmp2(layer,j)
      ww_k(layer, i)                 = ww_k(layer, i) +                                   &
        & predictors_k%watervapour(16, layer, i) * 1.25_jprb * ray_path(layer, j) *  &
        & (predictors%watervapour(2, layer, j)) ** 0.25_jprb
      ray_path_k(layer, i) = ray_path_k(layer, i) +      &
        & predictors_k%watervapour(16, layer, i) * tmp3(layer,j)
      wr_k(layer, i)                 = wr_k(layer, i) +                                  &
        & predictors_k%watervapour(15, layer, i) * 1.5_jprb * ray_path(layer, j) *  &
        & SQRT(predictors%watervapour(7, layer, j))
      ray_path_k(layer, i) = ray_path_k(layer, i) +      &
        & predictors_k%watervapour(15, layer, i) * 1.5_jprb * wr(layer, j) * SQRT(predictors%watervapour(7, layer, j))
      ww_k(layer, i)                 = ww_k(layer, i) +                                  &
        & predictors_k%watervapour(14, layer, i) * 1.5_jprb * ray_path(layer, j) *  &
        & SQRT(predictors%watervapour(2, layer, j))
      ray_path_k(layer, i) = ray_path_k(layer, i) +      &
        & predictors_k%watervapour(14, layer, i) * 1.5_jprb * ww(layer, j) * SQRT(predictors%watervapour(2, layer, j))
      ray_path_k(layer, i) = ray_path_k(layer, i) +                        &
        & predictors_k%watervapour(13, layer, i) * tmp4(layer,j)
      wr_k(layer, i)                 = wr_k(layer, i) +                                                              &
        & predictors_k%watervapour(13, layer, i) * tmp5(layer,j)
      wwr_k(layer, i)                = wwr_k(layer, i) -                                                            &
        & predictors_k%watervapour(13, layer, i) * tmp6(layer,j)
      sec_wrwr_k(layer, i)           = sec_wrwr_k(layer, i) + predictors_k%watervapour(12, layer, i) / wwr(layer, j)
      wwr_k(layer, i)                =      &
        & wwr_k(layer, i) - predictors_k%watervapour(12, layer, i) * sec_wrwr(layer, j) / wwr(layer, j) ** 2_jpim
      dt_k(layer, i)                 =      &
        & dt_k(layer, i) + predictors_k%watervapour(11, layer, i) * predictors%watervapour(5, layer, j)
      sec_wr_k(layer, i)             = sec_wr_k(layer, i) +      &
        & 0.5_jprb * predictors_k%watervapour(11, layer, i) * dt(layer, j) / predictors%watervapour(5, layer, j)
      sec_wr_k(layer, i)             =      &
        & sec_wr_k(layer, i) + predictors_k%watervapour(10, layer, i) * dtabs(layer, j) * dt(layer, j)
      dt_k(layer, i)                 = dt_k(layer, i) +      &
        & 2 * predictors_k%watervapour(10, layer, i) * predictors%watervapour(7, layer, j) * dtabs(layer, j)
      sec_wr_k(layer, i)             =      &
        & sec_wr_k(layer, i) + predictors_k%watervapour(9, layer, i) * 4._jprb * predictors%watervapour(8, layer, j)
      sec_wr_k(layer, i)             =      &
        & sec_wr_k(layer, i) + 3._jprb * predictors_k%watervapour(8, layer, i) * predictors%watervapour(1, layer, j)
      sec_wr_k(layer, i)             = sec_wr_k(layer, i) + predictors_k%watervapour(7, layer, i)
      sec_wr_k(layer, i)             =      &
        & sec_wr_k(layer, i) + 0.25_jprb * predictors_k%watervapour(6, layer, i) / predictors%watervapour(6, layer, j) ** 3
      sec_wr_k(layer, i)             =      &
        & sec_wr_k(layer, i) + 0.5_jprb * predictors_k%watervapour(5, layer, i) / predictors%watervapour(5, layer, j)
      dt_k(layer, i)                 =      &
        & dt_k(layer, i) + predictors_k%watervapour(4, layer, i) * predictors%watervapour(7, layer, j)
      sec_wr_k(layer, i)             = sec_wr_k(layer, i) + predictors_k%watervapour(4, layer, i) * dt(layer, j)
      ww_k(layer, i)                 = ww_k(layer, i) +      &
        & predictors_k%watervapour(3, layer, i) * 2._jprb * predictors%watervapour(2, layer, j) * ray_path(layer, j)
      ray_path_k(layer, i) = ray_path_k(layer, i) +      &
        & predictors_k%watervapour(3, layer, i) * 2._jprb * predictors%watervapour(2, layer, j) * ww(layer, j)
      ww_k(layer, i)                 =      &
        & ww_k(layer, i) + predictors_k%watervapour(2, layer, i) * ray_path(layer, j)
      ray_path_k(layer, i) =      &
        & ray_path_k(layer, i) + predictors_k%watervapour(2, layer, i) * ww(layer, j)
      sec_wr_k(layer, i)             =      &
        & sec_wr_k(layer, i) + 2._jprb * predictors_k%watervapour(1, layer, i) * predictors%watervapour(7, layer, j)
      sec_wr_k(layer, i)             = sec_wr_k(layer, i) + sec_wrwr_k(layer, i) * wr(layer, j)
      wr_k(layer, i)                 = wr_k(layer, i) + sec_wrwr_k(layer, i) * predictors%watervapour(7, layer, j)
      wr_k(layer, i)                 = wr_k(layer, i) + sec_wr_k(layer, i) * ray_path(layer, j)
      ray_path_k(layer, i) = ray_path_k(layer, i) + sec_wr_k(layer, i) * wr(layer, j)
    ENDDO

  ENDDO

!-----------------------------------------------------------------------------------------------------
!5.1 mixed gases
!-----------------------------------------------------------------------------------------------------
! John.Hague@ecmwf.int optimisation start
  Do j = 1, nprofiles
    Do layer = 1,prof_k(1) % nlayers
      tmp1(layer,j)=ray_path(layer, j)**1.5_jprb/SQRT(tr(layer,j))
    Enddo
  Enddo
! optimisation end

  DO i = 1, nchannels
    j = chanprof(i)%prof

    DO layer = 1, prof_k(i)%nlayers
      level = layer + 1
! X10
      tr_k(layer, i)                 = tr_k(layer, i) +      &
        & predictors_k%mixedgas(10, layer, i) * 0.5_jprb * tmp1(layer,j)
      ray_path_k(layer, i) = ray_path_k(layer, i) +      &
        & predictors_k%mixedgas(10, layer, i) * 1.5_jprb * SQRT(predictors%mixedgas(3, layer, j))
! X9
      tr_k(layer, i)                 =      &
        & tr_k(layer, i) + predictors_k%mixedgas(9, layer, i) * 3._jprb * ray_path(layer, j) * tr(layer, j) ** 2_jpim
      ray_path_k(layer, i) =      &
        & ray_path_k(layer, i) + predictors_k%mixedgas(9, layer, i) * tr(layer, j) ** 3_jpim
! X8
      tuwr_k(layer, i)               =      &
        & tuwr_k(layer, i) + predictors_k%mixedgas(8, layer, i) * ray_path(layer, j)
      ray_path_k(layer, i) =      &
        & ray_path_k(layer, i) + predictors_k%mixedgas(8, layer, i) * tuwr(layer, j)
! X7
      tuw_k(layer, i)                =      &
        & tuw_k(layer, i) + predictors_k%mixedgas(7, layer, i) * ray_path(layer, j)
      ray_path_k(layer, i) =      &
        & ray_path_k(layer, i) + predictors_k%mixedgas(7, layer, i) * tuw(layer, j)
! X6
      tr_k(layer, i)                 =      &
        & tr_k(layer, i) + predictors_k%mixedgas(6, layer, i) * 2._jprb * predictors%mixedgas(5, layer, j)
! X5
      tr_k(layer, i)                 = tr_k(layer, i) + predictors_k%mixedgas(5, layer, i)
! X4
      tr_k(layer, i)                 =      &
        & tr_k(layer, i) + predictors_k%mixedgas(4, layer, i) * 2._jprb * predictors%mixedgas(3, layer, j)
      ray_path_k(layer, i) =      &
        & ray_path_k(layer, i) + predictors_k%mixedgas(4, layer, i) * predictors%mixedgas(5, layer, j) ** 2_jpim
! X3
      tr_k(layer, i)                 =      &
        & tr_k(layer, i) + predictors_k%mixedgas(3, layer, i) * ray_path(layer, j)
      ray_path_k(layer, i) =      &
        & ray_path_k(layer, i) + predictors_k%mixedgas(3, layer, i) * tr(layer, j)
! X2
      ray_path_k(layer, i) =      &
        & ray_path_k(layer, i) + predictors_k%mixedgas(2, layer, i) * 2._jprb * ray_path(layer, j)
! X1
      ray_path_k(layer, i) = ray_path_k(layer, i) + predictors_k%mixedgas(1, layer, i)
    ENDDO

  ENDDO

!-------------------------------------------------------------------
!   calc adjoint of profile/reference sums
!-------------------------------------------------------------------

  DO i = 1, nchannels
    j = chanprof(i)%prof

    IF (coef%nch4 > 0) THEN
      IF (opts%rt_ir%ch4_Data) THEN
        sum1 = 0._jprb

        DO layer = prof_k(i)%nlayers, 1,  - 1
          sum1            = sum1 + ch4w_k(layer, i) / sum2_ch4w(layer, j)
          ch4_k(layer, i) = ch4_k(layer, i) + sum1 * coef%dpp(layer - 1)
        ENDDO

        sum1 = 0._jprb

        DO layer = prof_k(i)%nlayers, 1,  - 1
          sum1            = sum1 + ch4wr_k(layer, i) / sum2_ch4wr(layer, j)
          ch4_k(layer, i) = ch4_k(layer, i) + sum1 * coef%dpp(layer - 1) * t(layer, j)
          t_k(layer, i)   = t_k(layer, i) + sum1 * coef%dpp(layer - 1) * ch4(layer, j)
        ENDDO

      ELSE
        ch4_k(:, i) = 0._jprb

        sum1 = 0._jprb

        DO layer = prof_k(i)%nlayers, 1,  - 1
          sum1            = sum1 + ch4wr_k(layer, i) / sum2_ch4wr(layer, j)
          t_k(layer, i)   = t_k(layer, i) + sum1 * coef%dpp(layer - 1) * coef%ch4star(layer)
        ENDDO

      ENDIF
    ENDIF


    IF (coef%nco > 0) THEN
      IF (opts%rt_ir%co_Data) THEN
        sum1 = 0._jprb

        DO layer = prof_k(i)%nlayers, 1,  - 1
          sum1           = sum1 + cow_k(layer, i) / sum2_cow(layer, j)
          co_k(layer, i) = co_k(layer, i) + sum1 * coef%dpp(layer - 1)
        ENDDO

        sum1 = 0._jprb

        DO layer = prof_k(i)%nlayers, 1,  - 1
          sum1           = sum1 + cowr_k(layer, i) / sum2_cowr(layer, j)
          co_k(layer, i) = co_k(layer, i) + sum1 * coef%dpp(layer - 1) * t(layer, j)
          t_k(layer, i)  = t_k(layer, i) + sum1 * coef%dpp(layer - 1) * co(layer, j)
        ENDDO

      ELSE
        co_k(:, i) = 0._jprb

        sum1 = 0._jprb

        DO layer = prof_k(i)%nlayers, 1,  - 1
          sum1           = sum1 + cowr_k(layer, i) / sum2_cowr(layer, j)
          t_k(layer, i)  = t_k(layer, i) + sum1 * coef%dpp(layer - 1) * coef%costar(layer)
        ENDDO

      ENDIF
    ENDIF

    IF (coef%nn2o > 0) THEN
      IF (opts%rt_ir%n2o_Data) THEN
        sum1 = 0._jprb

        DO layer = prof_k(i)%nlayers, 1,  - 1
          sum1            = sum1 + n2ow_k(layer, i) / sum2_n2ow(layer, j)
          n2o_k(layer, i) = n2o_k(layer, i) + sum1 * coef%dpp(layer - 1)
        ENDDO

        sum1 = 0._jprb

        DO layer = prof_k(i)%nlayers, 1,  - 1
          sum1            = sum1 + n2owr_k(layer, i) / sum2_n2owr(layer, j)
          n2o_k(layer, i) = n2o_k(layer, i) + sum1 * coef%dpp(layer - 1) * t(layer, j)
          t_k(layer, i)   = t_k(layer, i) + sum1 * coef%dpp(layer - 1) * n2o(layer, j)
        ENDDO

      ELSE
        n2o_k(:, i) = 0._jprb

        sum1 = 0._jprb

        DO layer = prof_k(i)%nlayers, 1,  - 1
          sum1            = sum1 + n2owr_k(layer, i) / sum2_n2owr(layer, j)
          t_k(layer, i)   = t_k(layer, i) + sum1 * coef%dpp(layer - 1) * coef%n2ostar(layer)
        ENDDO

      ENDIF
    ENDIF


    IF (opts%rt_ir%co2_Data .AND. coef%nco2 > 0) THEN
      sum1 = 0._jprb

      DO layer = prof_k(i)%nlayers, 1,  - 1
        sum1            = sum1 + co2w_k(layer, i) / sum2_co2w(layer, j)
        co2_k(layer, i) = co2_k(layer, i) + sum1 * coef%dpp(layer - 1)
      ENDDO

    ELSE
      co2_k(:, i) = 0._jprb
    ENDIF


    IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN
      sum1 = 0._jprb

      DO layer = prof_k(i)%nlayers, 1,  - 1
        sum1          = sum1 + ow_k(layer, i) / sum2_ow(layer, j)
        o_k(layer, i) = o_k(layer, i) + sum1 * coef%dpp(layer - 1)
      ENDDO

    ELSE
      o_k(:, i) = 0._jprb
    ENDIF

!
    sum1 = 0._jprb

    DO layer = prof_k(i)%nlayers, 1,  - 1
      sum1          = sum1 + wwr_k(layer, i) / sum2_wwr(layer, j)
      w_k(layer, i) = w_k(layer, i) + sum1 * coef%dpp(layer - 1) * t(layer, j)
      t_k(layer, i) = t_k(layer, i) + sum1 * coef%dpp(layer - 1) * w(layer, j)
    ENDDO

!
    sum1 = 0._jprb

    DO layer = prof_k(i)%nlayers, 1,  - 1
      sum1          = sum1 + ww_k(layer, i) / sum2_ww(layer, j)
      w_k(layer, i) = w_k(layer, i) + sum1 * coef%dpp(layer - 1)
    ENDDO

    sum1 = 0._jprb

    IF (coef%nco2 > 0) THEN
      DO layer = prof_k(i)%nlayers, 1,  - 1
        sum1          = sum1 + twr_k(layer, i) / sum2_twr(layer, j)
        t_k(layer, i) = t_k(layer, i) + sum1 * coef%dpp(layer - 1)
      ENDDO
    ENDIF

    tuw_k(2:prof(j)%nlayers, i) = tuw_k(2:prof(j)%nlayers, i) + tuwr_k(2:prof(j)%nlayers, i)
    t_k(1, i)                   = t_k(1, i) + tuwr_k(1, i) * coef%dpp(0) / (coef%dpp(0) * coef%tstar(1))
    sum1 = 0._jprb

    DO layer = prof_k(i)%nlayers, 1,  - 1
      sum1          = sum1 + tuw_k(layer, i) / sum2_tuw(layer, j)
      t_k(layer, i) = t_k(layer, i) + sum1
    ENDDO


    DO layer = prof_k(i)%nlayers, 2,  - 1
      tw_k(layer - 1, i) = tw_k(layer - 1, i) + tw_k(layer, i)
      tr_k(layer - 1, i) = tr_k(layer - 1, i) + tw_k(layer, i) * coef%dpp(layer - 1)
    ENDDO

!-------------------------------------------------------------------
!   calc adjoint of profile deviations
!-------------------------------------------------------------------
! All assignments 1:prof(1) % nlayers

    DO layer = 1, prof_k(i)%nlayers

      IF (opts%rt_ir%ch4_Data .AND. coef%nch4 > 0) THEN
        ch4_k(layer, i) = ch4_k(layer, i) + ch4r_k(layer, i) / coef%ch4star(layer)
      ELSE
        ch4_k(layer, i) = 0._jprb
      ENDIF


      IF (opts%rt_ir%co_Data .AND. coef%nco > 0) THEN
        co_k(layer, i) = co_k(layer, i) + cor_k(layer, i) / coef%costar(layer)
      ELSE
        co_k(layer, i) = 0._jprb
      ENDIF


      IF (opts%rt_ir%n2o_Data .AND. coef%nn2o > 0) THEN
        n2o_k(layer, i) = n2o_k(layer, i) + n2or_k(layer, i) / coef%n2ostar(layer)
      ELSE
        n2o_k(layer, i) = 0._jprb
      ENDIF


      IF (opts%rt_ir%co2_Data .AND. coef%nco2 > 0) THEN
        co2_k(layer, i) = co2_k(layer, i) + co2r_k(layer, i) / coef%co2star(layer)
      ELSE
        co2_k(layer, i) = 0._jprb
      ENDIF


      IF (coef%nozone > 0) THEN
        t_k(layer, i) = t_k(layer, i) + tro_k(layer, i) / coef%to3star(layer)
        IF (opts%rt_ir%ozone_Data) THEN
          o_k(layer, i) = o_k(layer, i) + or_k(layer, i) / coef%ostar(layer)
        ELSE
          o_k(layer, i) = 0._jprb
        ENDIF
      ENDIF

      w_k(layer, i) = w_k(layer, i) + wr_k(layer, i) / coef%wstar(layer)
      t_k(layer, i) = t_k(layer, i) + tr_k(layer, i) / coef%tstar(layer)

      IF (coef%nozone > 0) t_k(layer, i) = t_k(layer, i) + dto_k(layer, i)

      t_k(layer, i) = t_k(layer, i) + dt_k(layer, i)
    ENDDO

!-------------------------------------------------------------------
!   calc adjoint of profile layer means
!-------------------------------------------------------------------

    IF (opts%rt_ir%ch4_Data .AND. coef%nch4 > 0) THEN

      DO level = 2, prof_k(i)%nlevels

        IF (opts%rt_ir%pc%addpc) THEN
          layer           = level - 1
          ch4_k(layer, i) = 0._jprb
        ELSE
          layer = level - 1
          prof_k(i)%ch4(level - 1) = prof_k(i)%ch4(level - 1) + 0.5_jprb * ch4_k(layer, i)
          prof_k(i)%ch4(level)     = prof_k(i)%ch4(level) + 0.5_jprb * ch4_k(layer, i)
        ENDIF

      ENDDO

    ENDIF


    IF (opts%rt_ir%n2o_Data .AND. coef%nn2o > 0) THEN

      DO level = 2, prof_k(i)%nlevels

        IF (opts%rt_ir%pc%addpc) THEN
          layer           = level - 1
          n2o_k(layer, i) = 0._jprb
        ELSE
          layer = level - 1
          prof_k(i)%n2o(level - 1) = prof_k(i)%n2o(level - 1) + 0.5_jprb * n2o_k(layer, i)
          prof_k(i)%n2o(level)     = prof_k(i)%n2o(level) + 0.5_jprb * n2o_k(layer, i)
        ENDIF

      ENDDO

    ENDIF


    IF (opts%rt_ir%co_Data .AND. coef%nco > 0) THEN

      DO level = 2, prof_k(i)%nlevels

        IF (opts%rt_ir%pc%addpc) THEN
          layer          = level - 1
          co_k(layer, i) = 0._jprb
        ELSE
          layer = level - 1
          prof_k(i)%co(level - 1) = prof_k(i)%co(level - 1) + 0.5_jprb * co_k(layer, i)
          prof_k(i)%co(level)     = prof_k(i)%co(level) + 0.5_jprb * co_k(layer, i)
        ENDIF

      ENDDO

    ENDIF


    IF (opts%rt_ir%co2_Data .AND. coef%nco2 > 0) THEN

      DO level = 2, prof_k(i)%nlevels

        IF (opts%rt_ir%pc%addpc) THEN
          layer           = level - 1
          co2_k(layer, i) = 0._jprb
        ELSE
          layer = level - 1
          prof_k(i)%co2(level - 1) = prof_k(i)%co2(level - 1) + 0.5_jprb * co2_k(layer, i)
          prof_k(i)%co2(level)     = prof_k(i)%co2(level) + 0.5_jprb * co2_k(layer, i)
        ENDIF

      ENDDO

    ENDIF


    IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN

      DO level = 2, prof_k(i)%nlevels
        layer = level - 1
        prof_k(i)%o3(level - 1) = prof_k(i)%o3(level - 1) + 0.5_jprb * o_k(layer, i)
        prof_k(i)%o3(level)     = prof_k(i)%o3(level) + 0.5_jprb * o_k(layer, i)
      ENDDO

    ENDIF


    DO level = 2, prof_k(i)%nlevels
      layer = level - 1
      prof_k(i)%q(level - 1) = prof_k(i)%q(level - 1) + 0.5_jprb * w_k(layer, i)
      prof_k(i)%q(level)     = prof_k(i)%q(level) + 0.5_jprb * w_k(layer, i)
      prof_k(i)%t(level - 1) = prof_k(i)%t(level - 1) + 0.5_jprb * t_k(layer, i)
      prof_k(i)%t(level)     = prof_k(i)%t(level) + 0.5_jprb * t_k(layer, i)
    ENDDO

  ENDDO
! End of channel loop
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_9_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setpredictors_9_k
