!
SUBROUTINE rttov_setpredictors_9( &
            & opts,        &
            & prof,        &
            & geom,        &
            & ray_path,    &
            & coef_pccomp, &
            & coef,        &
            & pred_info,   &
            & predictors)
! Description
! To calculate and store the predictors required for the
! optical depth regression for emissive (thermal) or
! solar channels.
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
!  1.0   01/06/2005  Marco Matricardi (ECMWF):
!           --       New routine based on rttov_setpredictors.F90.
!           --       Variable trace gases CO2, N2O,CO and CH4
!           --       introduced for IASI and AIRS.
!           --       Altitude dependent local zenith angle also
!           --       introduced.
!  1.1   22/08/2007  Optimised (D Salmond)
!  1.2   27/02/2009  Profile levels to include ToA. Distinguish arrays
!                    in raytracing (on levels) from all others (on
!                    layers). Predictors prepared to maintain agreement
!                    with the RTTOV-9 prediction scheme (P. Rayer)
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
       & rttov_coef,        &
       & rttov_options,     &
       & rttov_coef_pccomp, &
       & profile_Type,      &
       & geometry_Type,     &
       & predictors_Type,   &
       & rttov_path_pred
  USE parkind1, ONLY : jprb
!INTF_OFF
  USE rttov_const, ONLY : gravity, sensor_id_mw
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options    ), INTENT(IN)    :: opts
  TYPE(profile_Type     ), INTENT(IN)    :: prof(:)         ! profile (ppmv dry)
  TYPE(geometry_Type    ), INTENT(IN)    :: geom(SIZE(prof))! geometry
  REAL(jprb             ), INTENT(IN)    :: ray_path(prof(1)%nlayers,SIZE(prof))
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  TYPE(rttov_coef       ), INTENT(IN)    :: coef            ! coefficients
  TYPE(predictors_Type  ), INTENT(INOUT) :: pred_info       ! predictor structure
  TYPE(rttov_path_pred  ), INTENT(INOUT) :: predictors      ! predictors
!INTF_END
!local variables:
  INTEGER(KIND=jpim) :: level, layer, iprof
! user profile
  REAL   (KIND=jprb) :: t(SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: w(SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: o(SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: co2  (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: n2o  (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: co   (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: ch4  (SIZE(prof(1)%p)-1)
! reference profile
  REAL   (KIND=jprb) :: tr   (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: tro  (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: wr   (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: or   (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: co2r (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: cor  (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: ch4r (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: n2or (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: twr  (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: tuw  (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: tuwr (SIZE(prof(1)%p)-1)
! user - reference
  REAL   (KIND=jprb) :: dt   (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: dto  (SIZE(prof(1)%p)-1)
! pressure weighted
  REAL   (KIND=jprb) :: tw   (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: ww   (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: ow   (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: co2w (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: cow  (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: n2ow (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: ch4w (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: wwr  (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: n2owr(SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: cowr (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: ch4wr(SIZE(prof(1)%p)-1)
! intermediate variables
  REAL   (KIND=jprb) :: sum1 , sum2 , sum3 , sum4
  REAL   (KIND=jprb) :: wwr_r
  REAL   (KIND=jprb) :: ztemp

  REAL   (KIND=jprb) :: deltac       (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: sec_or       (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: sec_wr       (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: sec_wrwr     (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: tr_sq        (SIZE(prof(1)%p)-1)
  REAL   (KIND=jprb) :: ray_path_sqrt(SIZE(prof(1)%p)-1)
  INTEGER(KIND=jpim) :: nprofiles
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!------------------------------------------------------------------------------
! 1 profile layer quantities
!------------------------------------------------------------------------------
! layer N-1 lies between levels N-1 and N
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_9', 0_jpim, ZHOOK_HANDLE)
  nprofiles = size(prof)

  IF (opts%rt_ir%pc%addpc) THEN
!-CO2--------------------------------------------------------------------------

    DO layer = 1, prof(1)%nlayers
      level    = layer + 1
      IF (opts%rt_ir%co2_Data .AND. coef%nco2 > 0) THEN
        co2(layer) = (coef_pccomp%co2_pc_ref(level - 1) + coef_pccomp%co2_pc_ref(level)) * 0.5_jprb
      ENDIF
    ENDDO

!-N2O--------------------------------------------------------------------------
    DO layer = 1, prof(1)%nlayers
      level    = layer + 1
      IF (opts%rt_ir%n2o_Data .AND. coef%nn2o > 0) THEN
        n2o(layer) = (coef_pccomp%n2o_pc_ref(level - 1) + coef_pccomp%n2o_pc_ref(level)) * 0.5_jprb
      ENDIF
    ENDDO

!-CO---------------------------------------------------------------------------
    DO layer = 1, prof(1)%nlayers
      level    = layer + 1
      IF (opts%rt_ir%co_Data .AND. coef%nco > 0) THEN
        co(layer) = (coef_pccomp%co_pc_ref(level - 1) + coef_pccomp%co_pc_ref(level)) * 0.5_jprb
      ENDIF
    ENDDO

!-CH4--------------------------------------------------------------------------
    DO layer = 1, prof(1)%nlayers
      level    = layer + 1
      IF (opts%rt_ir%ch4_Data .AND. coef%nch4 > 0) THEN
        ch4(layer) = (coef_pccomp%ch4_pc_ref(level - 1) + coef_pccomp%ch4_pc_ref(level)) * 0.5_jprb
      ENDIF
    ENDDO

  ENDIF ! addpc


  DO iprof = 1, nprofiles

    DO layer = 1, prof(1)%nlayers
      level    = layer + 1
!-Temperature  ----------------------------------------------------------------
      t(layer) = (prof(iprof)%t(level - 1) + prof(iprof)%t(level)) * 0.5_jprb
!-H2O--------------------------------------------------------------------------
      w(layer) = (prof(iprof)%q(level - 1) + prof(iprof)%q(level)) * 0.5_jprb
!-O3---------------------------------------------------------------------------
      IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0)     &
        &  o(layer) = (prof(iprof)%o3(level - 1) + prof(iprof)%o3(level)) * 0.5_jprb
    ENDDO

    IF (.not. opts%rt_ir%pc%addpc) THEN

      DO layer = 1, prof(1)%nlayers
        level    = layer + 1

!-CO2--------------------------------------------------------------------------

        IF (opts%rt_ir%co2_Data .AND. coef%nco2 > 0) THEN
          co2(layer) = (prof(iprof)%co2(level - 1) + prof(iprof)%co2(level)) * 0.5_jprb
        ENDIF

!-N2O--------------------------------------------------------------------------

        IF (opts%rt_ir%n2o_Data .AND. coef%nn2o > 0) THEN
          n2o(layer) = (prof(iprof)%n2o(level - 1) + prof(iprof)%n2o(level)) * 0.5_jprb
        ENDIF

!-CO---------------------------------------------------------------------------

        IF (opts%rt_ir%co_Data .AND. coef%nco > 0) THEN
          co(layer) = (prof(iprof)%co(level - 1) + prof(iprof)%co(level)) * 0.5_jprb
        ENDIF

!-CH4--------------------------------------------------------------------------

        IF (opts%rt_ir%ch4_Data .AND. coef%nch4 > 0) THEN
          ch4(layer) = (prof(iprof)%ch4(level - 1) + prof(iprof)%ch4(level)) * 0.5_jprb
        ENDIF

      ENDDO ! layers

    ENDIF ! not addpc


!------------------------------------------------------------------------------
! 2 calculate deviations from reference profile (layers)
! if no input O3 profile we still use the input temperature profile for dto
!-----------------------------------------------------------------------------
! All assignments 1:prof(1) % nlayers
    dt(:) = t(:) - coef%tstar(:)
    IF (coef%nozone > 0) dto(:) = t(:) - coef%to3star(:)

!------------------------------------------------------------------------------
! 3 calculate (profile / reference profile) ratios; tr,wr,or,co2r etc.
!------------------------------------------------------------------------------
! All assignments 1:prof(1) % nlayers
!-Temperature------------------------------------------------------------------
    tr(:) = t(:) / coef%tstar(:)
!-H2O--------------------------------------------------------------------------
    wr(:) = w(:) / coef%wstar(:)
!-Ozone------------------------------------------------------------------------
! if no input O3 profile, set or to reference value (or =1), but still use
! input temperature profile for tro

    IF (coef%nozone > 0) THEN
      tro(:) = t(:) / coef%to3star(:)
      IF (opts%rt_ir%ozone_Data) THEN
        or(:)  = o(:) / coef%ostar(:)
      ELSE
        or(:)  = 1._jprb
      ENDIF
    ENDIF

!-CO2--------------------------------------------------------------------------
! if no input CO2 profile, set to reference value (co2r=1)

    IF (opts%rt_ir%co2_Data .AND. coef%nco2 > 0) THEN
      co2r(:) = co2(:) / coef%co2star(:)
    ELSE
      co2r(:) = 1._jprb
    ENDIF

!-N2O--------------------------------------------------------------------------
! if no input N2O profile, set to reference value (n2or=1)

    IF (opts%rt_ir%n2o_Data .AND. coef%nn2o > 0) THEN
      n2or(:) = n2o(:) / coef%n2ostar(:)
    ELSE
      n2or(:) = 1._jprb
    ENDIF

!-CO---------------------------------------------------------------------------
! if no input CO profile, set to reference value (cor=1)

    IF (opts%rt_ir%co_Data .AND. coef%nco > 0) THEN
      cor(:) = co(:) / coef%costar(:)
    ELSE
      cor(:) = 1._jprb
    ENDIF

!-CH4--------------------------------------------------------------------------
! if no input CH4 profile, set to reference value (ch4r=1)

    IF (opts%rt_ir%ch4_Data .AND. coef%nch4 > 0) THEN
      ch4r(:) = ch4(:) / coef%ch4star(:)
    ELSE
      ch4r(:) = 1._jprb
    ENDIF

!------------------------------------------------------------------------------
! 4 calculate profile / reference profile sums: tw, ww, ow, twr, co2w etc.
!------------------------------------------------------------------------------
!-Temperature------------------------------------------------------------------
    tw(1) = 0._jprb

    DO layer = 2, prof(1)%nlayers
! cumulate overlying layers (tr relates to same layer as dpp)
! do not need dpp(0) to start
      tw(layer) = tw(layer - 1) + coef%dpp(layer - 1) * tr(layer - 1)
    ENDDO

    sum1 = 0._jprb
    sum2 = 0._jprb

    DO layer = 1, prof(1)%nlayers
! cumulating
      sum1       = sum1 + t(layer)
      sum2       = sum2 + coef%tstar(layer)
      tuw(layer) = sum1 / sum2
    ENDDO

! special case
    tuwr(1)                 = coef%dpp(0) * t(1) / (coef%dpp(0) * coef%tstar(1))
! otherwise
    tuwr(2:prof(1)%nlayers) = tuw(2:prof(1)%nlayers)

    
    ! twr used only in CO2 predictors
    IF (coef%nco2 > 0) THEN
      sum1 = 0._jprb
      sum2 = 0._jprb
  
      DO layer = 1, prof(1)%nlayers
! cumulate overlying layers (t,tstar relate to layer below dpp)
! need dpp(0) to start
        sum1       = sum1 + coef%dpp(layer - 1) * t(layer)
        sum2       = sum2 + coef%dpp(layer - 1) * coef%tstar(layer)
        twr(layer) = sum1 / sum2
      ENDDO
    ENDIF

!-H2O--------------------------------------------------------------------------
    sum1 = 0._jprb
    sum2 = 0._jprb
    sum3 = 0._jprb
    sum4 = 0._jprb

    DO layer = 1, prof(1)%nlayers
! cumulate overlying layers (w,wstar relate to layer below dpp)
! need dpp(0) to start
      sum1       = sum1 + coef%dpp(layer - 1) * w(layer)
      sum2       = sum2 + coef%dpp(layer - 1) * coef%wstar(layer)
      sum3       = sum3 + coef%dpp(layer - 1) * w(layer) * t(layer)
      sum4       = sum4 + coef%dpp(layer - 1) * coef%wstar(layer) * coef%tstar(layer)
      ww(layer)  = sum1 / sum2
      wwr(layer) = sum3 / sum4
    ENDDO

!-O3---------------------------------------------------------------------------
! if no input O3 profile, set to reference value (ow =1)

    IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN
      sum1 = 0._jprb
      sum2 = 0._jprb

      DO layer = 1, prof(1)%nlayers
! cumulate overlying layers (o,ostar relate to layer below dpp)
! need dpp(0) to start
        sum1      = sum1 + coef%dpp(layer - 1) * o(layer)
        sum2      = sum2 + coef%dpp(layer - 1) * coef%ostar(layer)
        ow(layer) = sum1 / sum2
      ENDDO

    ELSE
      ow(:) = 1._jprb
    ENDIF

!-CO2--------------------------------------------------------------------------
! if no input co2 profile, set to reference value (co2w=1)

    IF (opts%rt_ir%co2_Data .AND. coef%nco2 > 0) THEN
      sum1 = 0._jprb
      sum2 = 0._jprb

      DO layer = 1, prof(1)%nlayers
! cumulate overlying layers (co2,co2star relate to layer below dpp)
! need dpp(0) to start
        sum1        = sum1 + coef%dpp(layer - 1) * co2(layer)
        sum2        = sum2 + coef%dpp(layer - 1) * coef%co2star(layer)
        co2w(layer) = sum1 / sum2
      ENDDO

    ELSE
      co2w(:) = 1._jprb
    ENDIF

!-N2O--------------------------------------------------------------------------
! if no input n2o profile, set to reference value (n2ow=1)

    IF (coef%nn2o > 0) THEN
      IF (opts%rt_ir%n2o_Data) THEN
        sum1 = 0._jprb
        sum2 = 0._jprb
        sum3 = 0._jprb
        sum4 = 0._jprb
          
        DO layer = 1, prof(1)%nlayers
! cumulate overlying layers (n2o,n2ostar relate to layer below dpp)
! need dpp(0) to start
          sum1         = sum1 + coef%dpp(layer - 1) * n2o(layer)
          sum2         = sum2 + coef%dpp(layer - 1) * coef%n2ostar(layer)
          sum3         = sum3 + coef%dpp(layer - 1) * n2o(layer) * t(layer)
          sum4         = sum4 + coef%dpp(layer - 1) * coef%n2ostar(layer) * coef%tstar(layer)
          n2ow(layer)  = sum1 / sum2
          n2owr(layer) = sum3 / sum4
        ENDDO
  
      ELSE
        n2ow(:)  = 1._jprb
      
! no input n2o profile so use the n2o reference profile, 
! but still use the input temperature profile to calculate n2owr
        sum3 = 0._jprb
        sum4 = 0._jprb
    
        DO layer = 1, prof(1)%nlayers
! cumulate overlying layers (n2o,n2ostar relate to layer below dpp)
! need dpp(0) to start
          sum3         = sum3 + coef%dpp(layer - 1) * coef%n2ostar(layer) * t(layer)
          sum4         = sum4 + coef%dpp(layer - 1) * coef%n2ostar(layer) * coef%tstar(layer)
          n2owr(layer) = sum3 / sum4
        ENDDO
      ENDIF
    ENDIF
    
!-CO---------------------------------------------------------------------------
! if no input co profile, set to reference value (cow=1 )

    IF (coef%nco > 0) THEN
      IF (opts%rt_ir%co_Data) THEN
        sum1 = 0._jprb
        sum2 = 0._jprb
        sum3 = 0._jprb
        sum4 = 0._jprb
  
        DO layer = 1, prof(1)%nlayers
! cumulate overlying layers (co,costar relate to layer below dpp)
! need dpp(0) to start
          sum1        = sum1 + coef%dpp(layer - 1) * co(layer)
          sum2        = sum2 + coef%dpp(layer - 1) * coef%costar(layer)
          sum3        = sum3 + coef%dpp(layer - 1) * co(layer) * t(layer)
          sum4        = sum4 + coef%dpp(layer - 1) * coef%costar(layer) * coef%tstar(layer)
          cow(layer)  = sum1 / sum2
          cowr(layer) = sum3 / sum4
        ENDDO
  
      ELSE
        cow(:)  = 1._jprb
    
! no input co profile so use the co reference profile, 
! but still use the input temperature profile to calculate cowr
        sum3 = 0._jprb
        sum4 = 0._jprb
    
        DO layer = 1, prof(1)%nlayers
! cumulate overlying layers (co,costar relate to layer below dpp)
! need dpp(0) to start
          sum3        = sum3 + coef%dpp(layer - 1) * coef%costar(layer) * t(layer)
          sum4        = sum4 + coef%dpp(layer - 1) * coef%costar(layer) * coef%tstar(layer)
          cowr(layer) = sum3 / sum4
        ENDDO
      ENDIF
    ENDIF
    
!-CH4--------------------------------------------------------------------------
! if no input ch4 profile, set to reference value (ch4w=1 )

    IF (coef%nch4 > 0) THEN
      IF (opts%rt_ir%ch4_Data) THEN
        sum1 = 0._jprb
        sum2 = 0._jprb
        sum3 = 0._jprb
        sum4 = 0._jprb
  
        DO layer = 1, prof(1)%nlayers
! cumulate overlying layers (ch4,ch4star relate to layer below dpp)
! need dpp(0) to start
          sum1         = sum1 + coef%dpp(layer - 1) * ch4(layer)
          sum2         = sum2 + coef%dpp(layer - 1) * coef%ch4star(layer)
          sum3         = sum3 + coef%dpp(layer - 1) * ch4(layer) * t(layer)
          sum4         = sum4 + coef%dpp(layer - 1) * coef%ch4star(layer) * coef%tstar(layer)
          ch4w(layer)  = sum1 / sum2
          ch4wr(layer) = sum3 / sum4
        ENDDO
  
      ELSE
        ch4w(:)  = 1._jprb

! no input ch4 profile so use the ch4 reference profile, 
! but still use the input temperature profile to calculate ch4wr
        sum3 = 0._jprb
        sum4 = 0._jprb
    
        DO layer = 1, prof(1)%nlayers
! cumulate overlying layers (ch4,ch4star relate to layer below dpp)
! need dpp(0) to start
          sum3         = sum3 + coef%dpp(layer - 1) * coef%ch4star(layer) * t(layer)
          sum4         = sum4 + coef%dpp(layer - 1) * coef%ch4star(layer) * coef%tstar(layer)
          ch4wr(layer) = sum3 / sum4
        ENDDO
      ENDIF
    ENDIF
    
!5) set predictors for RTTOV-9 options (same for RTTOV-10)
!--

    DO layer = 1, prof(1)%nlayers
      level = layer + 1 ! raytracing % pathsat(layer,i) is angle at lower boundary of layer
      ray_path_sqrt(layer)                     = SQRT(ray_path(layer,iprof))
!5.1 mixed gases
!---
      tr_sq(layer)                             = tr(layer) * tr(layer)
      predictors%mixedgas(1, layer, iprof)     = ray_path(layer,iprof)
      predictors%mixedgas(2, layer, iprof)     = ray_path(layer,iprof) * ray_path(layer,iprof)
      predictors%mixedgas(3, layer, iprof)     = ray_path(layer,iprof) * tr(layer)
      predictors%mixedgas(4, layer, iprof)     = ray_path(layer,iprof) * tr_sq(layer)
      predictors%mixedgas(5, layer, iprof)     = tr(layer)
      predictors%mixedgas(6, layer, iprof)     = tr_sq(layer)
      predictors%mixedgas(7, layer, iprof)     = ray_path(layer,iprof) * tuw(layer)
      predictors%mixedgas(8, layer, iprof)     = ray_path(layer,iprof) * tuwr(layer)
      predictors%mixedgas(9, layer, iprof)     = ray_path(layer,iprof) * tr(layer) ** 3_jpim
      predictors%mixedgas(10, layer, iprof)    = ray_path(layer,iprof) * SQRT(ray_path(layer,iprof) * tr(layer))
!5.2 water vapour line transmittance based on RTIASI but with pred 9 removed
!----------------
      sec_wr(layer)                            = ray_path(layer,iprof) * wr(layer)
      sec_wrwr(layer)                          = sec_wr(layer) * wr(layer)
      wwr_r = 1.0_jprb / wwr(layer)
!predictors % watervapour(:,:,iprof) = 0._jprb
      predictors%watervapour(1, layer, iprof)  = sec_wr(layer) * sec_wr(layer)
      predictors%watervapour(2, layer, iprof)  = ray_path(layer,iprof) * ww(layer)
      predictors%watervapour(3, layer, iprof)  = (ray_path(layer,iprof) * ww(layer)) ** 2_jpim
      predictors%watervapour(4, layer, iprof)  = sec_wr(layer) * dt(layer)
      predictors%watervapour(5, layer, iprof)  = SQRT(sec_wr(layer))
      predictors%watervapour(6, layer, iprof)  = sec_wr(layer) ** 0.25_jprb
      predictors%watervapour(7, layer, iprof)  = sec_wr(layer)
      predictors%watervapour(8, layer, iprof)  = sec_wr(layer) ** 3_jpim
      predictors%watervapour(9, layer, iprof)  = sec_wr(layer) ** 4_jpim
      predictors%watervapour(10, layer, iprof) = sec_wr(layer) * dt(layer) * ABS(dt(layer))
      predictors%watervapour(11, layer, iprof) = SQRT(sec_wr(layer)) * dt(layer)
      predictors%watervapour(12, layer, iprof) = sec_wrwr(layer) * wwr_r
      predictors%watervapour(13, layer, iprof) = ray_path_sqrt(layer) * wr(layer) ** 1.5_jprb * wwr_r
      predictors%watervapour(14, layer, iprof) = (ray_path(layer,iprof) * ww(layer)) ** 1.5_jprb
      predictors%watervapour(15, layer, iprof) = sec_wr(layer) ** 1.5_jprb
      predictors%watervapour(16, layer, iprof) = (ray_path(layer,iprof) * ww(layer)) ** 1.25_jprb
      predictors%watervapour(17, layer, iprof) = SQRT(sec_wr(layer)) * wr(layer)
      predictors%watervapour(18, layer, iprof) = sec_wr(layer) ** 1.5_jprb * dt(layer)
      predictors%watervapour(19, layer, iprof) = ray_path(layer,iprof) * wr(layer) ** 2_jpim / ww(layer)
!5.3 water vapour continuum transmittance based on RTIASI
!----------------
!

      IF (coef%nwvcont > 0) THEN
!predictors % wvcont(:,:,iprof)  = 0._jprb
        ztemp = 1.0_jprb / tr(layer)
        predictors%wvcont(1, layer, iprof) = sec_wrwr(layer) * ztemp
        predictors%wvcont(3, layer, iprof) = sec_wr(layer) * ztemp
        ztemp = 1.0_jprb / tr_sq(layer)
        predictors%wvcont(2, layer, iprof) = sec_wrwr(layer) * ztemp * ztemp
        predictors%wvcont(4, layer, iprof) = sec_wr(layer) * ztemp
      ENDIF

!5.4 ozone
!---------
! if no input O3 profile, variables or, ow and dto have been set
! to the reference profile values (1, 1, 0)

      IF (coef%nozone > 0) THEN
        sec_or(layer)                      = or(layer) * ray_path(layer,iprof)
        predictors%ozone(1, layer, iprof)  = sec_or(layer)
        predictors%ozone(2, layer, iprof)  = SQRT(sec_or(layer))
        predictors%ozone(3, layer, iprof)  = sec_or(layer) * dto(layer)
        predictors%ozone(4, layer, iprof)  = sec_or(layer) * sec_or(layer)
        predictors%ozone(5, layer, iprof)  = SQRT(sec_or(layer)) * dto(layer)
        predictors%ozone(6, layer, iprof)  = sec_or(layer) * or(layer) * ow(layer)
        ztemp = 1.0_jprb / ow(layer)
        predictors%ozone(7, layer, iprof)  = SQRT(sec_or(layer)) * or(layer) * ztemp
        predictors%ozone(8, layer, iprof)  = sec_or(layer) * ow(layer)
        predictors%ozone(9, layer, iprof)  = sec_or(layer) * SQRT(ray_path(layer,iprof) * ow(layer))
        predictors%ozone(10, layer, iprof) = ray_path(layer,iprof) * ow(layer)
        predictors%ozone(11, layer, iprof) = ray_path(layer,iprof) * ow(layer) * ray_path(layer,iprof) * ow(layer)
        predictors%ozone(12, layer, iprof) = sec_or(layer) * ztemp
        predictors%ozone(13, layer, iprof) = (ray_path(layer,iprof) * ow(layer)) ** 1.75_jprb
        predictors%ozone(14, layer, iprof) = ray_path_sqrt(layer) * ow(layer) ** 2_jpim * dto(layer)
        predictors%ozone(15, layer, iprof) = ray_path(layer,iprof) * tro(layer) ** 3_jpim
      ENDIF

    ENDDO
! layers
!5.5 cloud
!---------

    IF (coef%id_sensor == sensor_id_mw) THEN
      IF (opts%rt_mw%clw_Data) THEN
        ztemp = 1.0_jprb / (4.3429_jprb * gravity)

        DO layer = 1, prof(1)%nlayers
          ! NB in RTTOV-10, layer 1 is the layer above level 2
          level = layer + 1
          deltac(layer)                = 0.1820_jprb * 100.0_jprb * coef%dp(layer) * ztemp
          predictors%clw(layer, iprof) = deltac(layer) * prof(iprof)%clw(level) * geom(iprof)%seczen
        ENDDO


        DO layer = 2, prof(1)%nlayers
          predictors%clw(layer, iprof) =      &
            & 0.5_jprb * (predictors%clw(layer, iprof) + deltac(layer) * prof(iprof)%clw(level - 1) * geom(iprof)%seczen)
          pred_info%ncloud            = 1_jpim
        ENDDO

      ELSE
        predictors%clw = 0._jprb
      ENDIF
    ENDIF
!5.6 carbon dioxide transmittance based on RTIASI
!-------------------------------------------------
!

    IF (coef%nco2 > 0) THEN

      DO layer = 1, prof(1)%nlayers
        level = layer + 1
        predictors%co2(1, layer, iprof)  = ray_path(layer,iprof) * co2r(layer)
        predictors%co2(2, layer, iprof)  = tr_sq(layer)
        predictors%co2(3, layer, iprof)  = ray_path(layer,iprof) * tr(layer)
        predictors%co2(4, layer, iprof)  = ray_path(layer,iprof) * tr_sq(layer)
        predictors%co2(5, layer, iprof)  = tr(layer)
        predictors%co2(6, layer, iprof)  = ray_path(layer,iprof)
        predictors%co2(7, layer, iprof)  = ray_path(layer,iprof) * twr(layer)
        predictors%co2(8, layer, iprof)  = (ray_path(layer,iprof) * co2w(layer)) ** 2_jpim
        predictors%co2(9, layer, iprof)  = twr(layer) * twr(layer) * twr(layer)
        predictors%co2(10, layer, iprof) = ray_path(layer,iprof) * twr(layer) * SQRT(tr(layer))
        predictors%co2(11, layer, iprof) = SQRT(ray_path(layer,iprof) * co2r(layer))
        predictors%co2(12, layer, iprof) = tr(layer) ** 3_jpim
        predictors%co2(13, layer, iprof) = ray_path(layer,iprof) * tr(layer) ** 3_jpim
        predictors%co2(14, layer, iprof) = ray_path_sqrt(layer) * tr(layer) ** 2_jpim * twr(layer) ** 3_jpim
        predictors%co2(15, layer, iprof) = tr(layer) ** 2_jpim * twr(layer) ** 2_jpim
      ENDDO

    ENDIF

!5.7 n2o transmittance based on RTIASI
!-------------------------------------
!

    IF (coef%nn2o > 0) THEN

      DO layer = 1, prof(1)%nlayers
        level = layer + 1
        predictors%n2o(1, layer, iprof)  = ray_path(layer,iprof) * n2or(layer)
        predictors%n2o(2, layer, iprof)  = SQRT(ray_path(layer,iprof) * n2or(layer))
        predictors%n2o(3, layer, iprof)  = ray_path(layer,iprof) * n2or(layer) * dt(layer)
        predictors%n2o(4, layer, iprof)  = (ray_path(layer,iprof) * n2or(layer)) ** 2_jpim
        predictors%n2o(5, layer, iprof)  = n2or(layer) * dt(layer)
        predictors%n2o(6, layer, iprof)  = (ray_path(layer,iprof) * n2or(layer)) ** 0.25_jprb
        predictors%n2o(7, layer, iprof)  = ray_path(layer,iprof) * n2ow(layer)
        predictors%n2o(8, layer, iprof)  = ray_path(layer,iprof) * n2owr(layer)
        predictors%n2o(9, layer, iprof)  = n2owr(layer)
        predictors%n2o(10, layer, iprof) = SQRT(ray_path(layer,iprof) * n2or(layer)) * n2or(layer) / n2ow(layer)
        predictors%n2o(11, layer, iprof) = (ray_path(layer,iprof) * n2owr(layer)) ** 2_jpim
        predictors%n2o(12, layer, iprof) = (ray_path(layer,iprof) * n2owr(layer)) ** 3_jpim
        predictors%n2o(13, layer, iprof) = ray_path(layer,iprof) ** 2_jpim * n2owr(layer) * dt(layer)
      ENDDO

    ENDIF

!5.8 co transmittance based on RTIASI
!------------------------------------

    IF (coef%nco > 0) THEN

      DO layer = 1, prof(1)%nlayers
        level = layer + 1
        predictors%co(1, layer, iprof)  = ray_path(layer,iprof) * cor(layer)
        predictors%co(2, layer, iprof)  = SQRT(ray_path(layer,iprof) * cor(layer))
        predictors%co(3, layer, iprof)  = ray_path(layer,iprof) * cor(layer) * dt(layer)
        predictors%co(4, layer, iprof)  = (ray_path(layer,iprof) * cor(layer)) ** 2_jpim
        predictors%co(5, layer, iprof)  = SQRT(ray_path(layer,iprof) * cor(layer)) * dt(layer)
        predictors%co(6, layer, iprof)  = (ray_path(layer,iprof) * cor(layer)) ** 0.25_jprb
        predictors%co(7, layer, iprof)  = ray_path(layer,iprof) * cor(layer) * dt(layer) * ABS(dt(layer))
        ztemp = 1.0_jprb / cow(layer)
        predictors%co(8, layer, iprof)  = ray_path(layer,iprof) * cor(layer) ** 2_jpim * ztemp
        predictors%co(9, layer, iprof)  = SQRT(ray_path(layer,iprof) * cor(layer)) * cor(layer) * ztemp
        predictors%co(10, layer, iprof) = ray_path(layer,iprof) * cor(layer) ** 2_jpim / SQRT(cow(layer))
        predictors%co(11, layer, iprof) = ray_path(layer,iprof) * cor(layer) ** 2_jpim / cow(layer) ** 0.25_jprb
        predictors%co(12, layer, iprof) = (ray_path(layer,iprof) * cowr(layer)) ** 0.4_jprb
        predictors%co(13, layer, iprof) = (ray_path(layer,iprof) * cowr(layer)) ** 0.25_jprb
      ENDDO

    ENDIF

!
!5.9 ch4             transmittance based on RTIASI
!-------------------------------------------------
!

    IF (coef%nch4 > 0) THEN

      DO layer = 1, prof(1)%nlayers
        level = layer + 1
        predictors%ch4(1, layer, iprof)  = ray_path(layer,iprof) * ch4r(layer)
        predictors%ch4(2, layer, iprof)  = SQRT(ray_path(layer,iprof) * ch4r(layer))
        predictors%ch4(3, layer, iprof)  = ray_path(layer,iprof) * ch4r(layer) * dt(layer)
        predictors%ch4(4, layer, iprof)  = (ray_path(layer,iprof) * ch4r(layer)) ** 2_jpim
        predictors%ch4(5, layer, iprof)  = ch4r(layer) * dt(layer)
        predictors%ch4(6, layer, iprof)  = (ray_path(layer,iprof) * ch4r(layer)) ** 0.25_jprb
        predictors%ch4(7, layer, iprof)  = ray_path(layer,iprof) * ch4wr(layer)
        predictors%ch4(8, layer, iprof)  = ch4wr(layer)
        predictors%ch4(9, layer, iprof)  = (ray_path(layer,iprof) * ch4w(layer)) ** 2_jpim
        predictors%ch4(10, layer, iprof) = ray_path(layer,iprof) * ch4w(layer)
        predictors%ch4(11, layer, iprof) = SQRT(ray_path(layer,iprof) * ch4r(layer)) * ch4r(layer) / ch4w(layer)
      ENDDO

    ENDIF

  ENDDO ! profiles

  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_9', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setpredictors_9
