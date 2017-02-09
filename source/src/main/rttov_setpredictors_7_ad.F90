!
SUBROUTINE rttov_setpredictors_7_ad( &
            & opts,          &
            & prof,          &
            & prof_ad,       &
            & geom,          &
            & coef,          &
            & aux,           &
            & aux_ad,        &
            & predictors,    &
            & predictors_ad, &
            & raytracing,    &
            & raytracing_ad)
!
! Description
! AD of rttov_setpredictors_7
! To calculate and store the profile variables (predictors) required
! in subsequent transmittance calculations.
! Code based on PRFTAU from previous versions of RTTOV
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
!  1.0       01/12/2002  New F90 code with structures (P Brunel A Smith)
!  1.1       04/12/2003  Optimisation (J Hague and D Salmond ECMWF)
!  1.2       29/03/2005  Add end of header comment (J. Cameron)
!  1.3       07/12/2005  Add surface humidity (R. Saunders)
!  1.4       03/03/2006  Marco Matricardi (ECMWF):
!               --       Altitude dependent local zenith angle introduced.
!  1.5       12/08/2008  Zeeman effect based on Yong Han (P. Rayer)
!  1.6       15/08/2009  User defined ToA. Layers distinct from levels (P.Rayer)
!  1.7       02/12/2009  Pathsat, Pathsun and related quantities are now layer arrays
!                        (Marco Matricardi).
!  1.8       17/06/2010  Combined non-Zeeman and Zeeman predictors for SSMIS for
!                        use with single coefficient file  (P Rayer)
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
! Imported Type Definitions:
  USE rttov_types, ONLY :  &
       & rttov_coef,      &
       & rttov_options,   &
       & profile_type,    &
       & geometry_type,   &
       & profile_aux,     &
       & rttov_path_pred, &
       & raytracing_type
!INTF_OFF
  USE rttov_const, ONLY : &
       & gravity,         &
       & sensor_id_mw,    &
       & inst_id_ssmis,   &
       & inst_id_ssmisz,  &
       & inst_id_amsua
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  USE rttov_math_mod
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options  ), INTENT(IN)    :: opts
  TYPE(profile_type   ), INTENT(IN)    :: prof(:)         ! profile (ppmv dry)
  TYPE(profile_type   ), INTENT(INOUT) :: prof_ad(:)
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  TYPE(geometry_type  ), INTENT(IN)    :: geom(:)
  TYPE(rttov_path_pred), INTENT(IN)    :: predictors
  TYPE(rttov_path_pred), INTENT(INOUT) :: predictors_ad
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing_ad
  TYPE(profile_aux    ), INTENT(IN)    :: aux
  TYPE(profile_aux    ), INTENT(INOUT) :: aux_ad
!INTF_END
!local variables:
  INTEGER(KIND=jpim) :: level, layer, lay
  INTEGER(KIND=jpim) :: i
  REAL   (KIND=jprb) :: deltac   (prof(1)%nlayers)

  REAL   (KIND=jprb) :: sec_wr   (prof(1)%nlayers), sec_wr_sqrt(prof(1)%nlayers) 
  REAL   (KIND=jprb) :: sec_or   (prof(1)%nlayers), sec_or_sqrt(prof(1)%nlayers) 
  REAL   (KIND=jprb) :: sec_or_ad(prof(1)%nlayers), sec_or_sqrt_ad(prof(1)%nlayers) 
  REAL   (KIND=jprb) :: sec_wr_ad(prof(1)%nlayers), sec_wr_sqrt_ad(prof(1)%nlayers) 
  REAL   (KIND=jprb) :: sec_wr_rsqrt(prof(1)%nlayers)

  REAL   (KIND=jprb) :: ztemp
  INTEGER(KIND=jpim) :: nprofiles, nlayers
  REAL   (KIND=jprb) :: ZHOOK_HANDLE

!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_7_AD', 0_jpim, ZHOOK_HANDLE)
  nprofiles = SIZE(prof)
  nlayers = prof(1)%nlayers

  IF(aux%on_coef_levels) raytracing_ad%pathsat = 0._jprb ! zero the COEF copy of pathsat

! cloud
!---------
  IF (coef%id_sensor == sensor_id_mw .AND. opts%rt_mw%clw_Data) THEN

    ztemp = 1.0_jprb / (4.3429_jprb * gravity)

    DO layer = 1, prof(1)%nlayers
      deltac(layer) = 0.1820_jprb * 100.0_jprb * coef%dp(layer) * ztemp
    ENDDO

!CDIR NOLOOPCHG

    DO layer = 2, prof(1)%nlayers
      level = layer + 1
      DO i = 1, nprofiles
        prof_ad(i)%clw(level - 1)   =      &
          & prof_ad(i)%clw(level - 1) + 0.5_jprb * predictors_ad%clw(layer, i) * deltac(layer) * geom(i)%seczen
        predictors_ad%clw(layer, i) = 0.5_jprb * predictors_ad%clw(layer, i)
      ENDDO
    ENDDO

!CDIR NOLOOPCHG

    DO layer = 1, prof(1)%nlayers
      level = layer + 1
      DO i = 1, nprofiles
        prof_ad(i)%clw(level) = prof_ad(i)%clw(level) + predictors_ad%clw(layer, i) * deltac(layer) * geom(i)%seczen
      ENDDO
    ENDDO
  ENDIF

!--
! mixed gases
!---

  IF (coef%id_inst == inst_id_ssmisz .AND.  coef%IncZeeman) THEN

      ! SSMISZ coefficient file - for testing only
      DO i = 1, nprofiles
        DO lay = 1, nlayers
          predictors_ad%mixedgas(10, lay, i) = predictors_ad%mixedgas(10, lay, i) + &
             predictors_ad%mixedgas(11, lay, i) * prof(i)%Be
          predictors_ad%mixedgas(3, lay, i) = predictors_ad%mixedgas(3, lay, i) + &
             predictors_ad%mixedgas(10, lay, i) * prof(i)%Be

          raytracing_ad%pathsat(lay, i) = raytracing_ad%pathsat(lay, i) + &
             prof(i)%Be ** 3_jpim * predictors_ad%mixedgas(9, lay, i)
          raytracing_ad%pathsat(lay, i) = raytracing_ad%pathsat(lay, i) + &
             prof(i)%Be * predictors_ad%mixedgas(8, lay, i)
          predictors_ad%mixedgas(6, lay, i) = predictors_ad%mixedgas(6, lay, i) + &
             predictors_ad%mixedgas(7, lay, i) / prof(i)%Be
          raytracing_ad%pathsat(lay, i) = raytracing_ad%pathsat(lay, i) + &
             predictors_ad%mixedgas(6, lay, i) / prof(i)%Be
          predictors_ad%mixedgas(2, lay, i) = predictors_ad%mixedgas(2, lay, i) + &
             predictors_ad%mixedgas(5, lay, i) * prof(i)%cosbk ** 2_jpim
          predictors_ad%mixedgas(2, lay, i) = predictors_ad%mixedgas(2, lay, i) + &
             predictors_ad%mixedgas(4, lay, i) / prof(i)%Be

          raytracing_ad%pathsat(lay, i) = raytracing_ad%pathsat(lay, i) + &
             prof(i)%cosbk ** 2_jpim * predictors_ad%mixedgas(3, lay, i)

          aux_ad%t_layer(lay, i) = &!aux_ad%t_layer(lay, i)
            -(predictors%mixedgas(2, lay, i) / &
             aux%t_layer(lay, i)) * predictors_ad%mixedgas(2, lay, i)

          raytracing_ad%pathsat(lay, i) = raytracing_ad%pathsat(lay, i) + &
             (300.0_jprb / aux%t_layer(lay, i)) * predictors_ad%mixedgas(2, lay, i)
          raytracing_ad%pathsat(lay, i) =  raytracing_ad%pathsat(lay, i) + predictors_ad%mixedgas(1, lay, i)

          ! Initialise some AD variables which are initialised in the ELSE clause for other instruments
          aux_ad%tr(lay, i) = 0._jprb
          aux_ad%tr_r(lay, i) = 0._jprb
          aux_ad%tw(lay, i) = 0._jprb
          aux_ad%tw_4rt(lay, i) = 0._jprb
          raytracing_ad%pathsat_sqrt(lay, i) = 0._jprb
        ENDDO
      ENDDO

  ELSE
    DO i = 1, nprofiles
      DO lay = 1, nlayers
        ! only effective for non-Zeeman chans 1-18,23-24 - coefficient file will have zeros for chan 19-22
        raytracing_ad%pathsat(lay, i) = raytracing_ad%pathsat(lay, i) + &
                                                 predictors_ad%mixedgas(1, lay, i) + &
          2._jprb * raytracing%pathsat(lay, i) * predictors_ad%mixedgas(2, lay, i) + &
          aux%tr(lay, i) *                       predictors_ad%mixedgas(3, lay, i) + &
          aux%tr(lay, i) * aux%tr(lay, i) *      predictors_ad%mixedgas(4, lay, i) + &
          aux%tw(lay, i) *                       predictors_ad%mixedgas(7, lay, i) + &
          aux%tw(lay, i) * aux%tr_r(lay, i) *    predictors_ad%mixedgas(8, lay, i)


        aux_ad%tr(lay, i) = &!aux_ad%tr(lay, i) + &
          raytracing%pathsat(lay, i) *                            predictors_ad%mixedgas(3, lay, i) + &
          2._jprb * aux%tr(lay, i) * raytracing%pathsat(lay, i) * predictors_ad%mixedgas(4, lay, i) + &
                                                                  predictors_ad%mixedgas(5, lay, i) + &
          2._jprb * aux%tr(lay, i) *                              predictors_ad%mixedgas(6, lay, i)


        aux_ad%tw(lay, i) = &!aux_ad%tw(lay, i) +
          raytracing%pathsat(lay, i) *                    predictors_ad%mixedgas(7, lay, i) + &
          raytracing%pathsat(lay, i) * aux%tr_r(lay, i) * predictors_ad%mixedgas(8, lay, i)

        aux_ad%tr_r(lay, i) = &!aux_ad%tr_r(lay, i) + &
          raytracing%pathsat(lay, i) * aux%tw(lay, i) * predictors_ad%mixedgas(8, lay, i)

        raytracing_ad%pathsat_sqrt(lay, i) = &!raytracing_ad%pathsat_sqrt(lay, i) +
                               predictors_ad%mixedgas(9, lay, i) + &
          aux%tw_4rt(lay, i) * predictors_ad%mixedgas(10, lay, i)

        aux_ad%tw_4rt(lay, i) = &!aux_ad%tw_4rt(lay, i) + &
          raytracing%pathsat_sqrt(lay, i) * predictors_ad%mixedgas(10, lay, i)
      ENDDO
    ENDDO

    IF (coef%id_inst == inst_id_ssmis .AND. coef%IncZeeman) THEN
      DO i = 1, nprofiles
        DO lay = 1, nlayers
          ! SSMIS with Zeeman coefficient file
          ! geomagnetic field variables (Be, cosbk) are part of the user input

          ! only effective for Zeeman chans 19-22 - coefficient file will have zeros for chan 1-18,23-24
          ! NB require prof(i) % Be >0. (divisor)

          ! X11 -> X21
          ! only effective for Zeeman chans 19-22 - coefficient file will have zeros for chan 1-18,23-24
          predictors_ad%mixedgas(20, lay, i) = predictors_ad%mixedgas(20, lay, i) + &
             predictors_ad%mixedgas(21, lay, i) * prof(i)%Be
          predictors_ad%mixedgas(13, lay, i) = predictors_ad%mixedgas(13, lay, i) + &
             predictors_ad%mixedgas(20, lay, i) * prof(i)%Be

          raytracing_ad%pathsat(lay, i) = raytracing_ad%pathsat(lay, i) + &
             prof(i)%Be ** 3_jpim * predictors_ad%mixedgas(19, lay, i)
          raytracing_ad%pathsat(lay, i) = raytracing_ad%pathsat(lay, i) + &
             prof(i)%Be * predictors_ad%mixedgas(18, lay, i)
          predictors_ad%mixedgas(16, lay, i) = predictors_ad%mixedgas(16, lay, i) + &
             predictors_ad%mixedgas(17, lay, i) / prof(i)%Be
          raytracing_ad%pathsat(lay, i) = raytracing_ad%pathsat(lay, i) + &
             predictors_ad%mixedgas(16, lay, i) / prof(i)%Be
          predictors_ad%mixedgas(12, lay, i) = predictors_ad%mixedgas(12, lay, i) + &
             predictors_ad%mixedgas(15, lay, i) * prof(i)%cosbk ** 2_jpim
          predictors_ad%mixedgas(12, lay, i) = predictors_ad%mixedgas(12, lay, i) + &
             predictors_ad%mixedgas(14, lay, i) / prof(i)%Be

          raytracing_ad%pathsat(lay, i) = raytracing_ad%pathsat(lay, i) + &
             prof(i)%cosbk ** 2_jpim * predictors_ad%mixedgas(13, lay, i)

          aux_ad%t_layer(lay, i) = &!aux_ad%t_layer(lay, i)
            -(predictors%mixedgas(12, lay, i) / &
             aux%t_layer(lay, i)) * predictors_ad%mixedgas(12, lay, i)

          raytracing_ad%pathsat(lay, i) = raytracing_ad%pathsat(lay, i) + &
             (300.0_jprb / aux%t_layer(lay, i)) * predictors_ad%mixedgas(12, lay, i)
          raytracing_ad%pathsat(lay, i) =  raytracing_ad%pathsat(lay, i) + predictors_ad%mixedgas(11, lay, i)
        ENDDO
      ENDDO
    ELSEIF (coef%id_inst == inst_id_amsua .AND. coef%IncZeeman) THEN
      ! AMSU-A with Zeeman coefficient file
      ! only effective for Zeeman chan 14 - coefficient file will have zeros for chan 1-13
      ! NB some of YH's original predictors omitted - effectively duplicated by predictors 1-4 above
      DO i = 1, nprofiles
        DO lay = 1, nlayers
          raytracing_ad%pathsat(lay, i) = raytracing_ad%pathsat(lay, i) + &
             prof(i)%cosbk ** 2_jpim * predictors_ad%mixedgas(11, lay, i)
          raytracing_ad%pathsat(lay, i) = raytracing_ad%pathsat(lay, i) +     &
            & 2.0_jprb * prof(i)%Be * raytracing%pathsat(lay, i) * predictors_ad%mixedgas(12, lay, i)
          raytracing_ad%pathsat(lay, i) = raytracing_ad%pathsat(lay, i) + &
             prof(i)%Be ** 3_jpim * predictors_ad%mixedgas(13, lay, i)
          raytracing_ad%pathsat(lay, i) = raytracing_ad%pathsat(lay, i) + &
            & 2.0_jprb * (prof(i)%cosbk * prof(i)%Be) ** 2_jpim * &
            raytracing%pathsat(lay, i) * predictors_ad%mixedgas(14, lay, i)
        ENDDO
      ENDDO
    ENDIF
  ENDIF

! water vapour - numbers in right hand are predictor numbers
! in the reference document for RTTOV7 (science and validation report)
!----------------

  DO i = 1, nprofiles
    sec_wr_sqrt(:) = raytracing%pathsat_sqrt(:, i) * aux%wr_sqrt(:, i)
    sec_wr_rsqrt(:) = raytracing%pathsat_rsqrt(:, i) * aux%wr_rsqrt(:, i)
    sec_wr(:) = sec_wr_sqrt(:)**2_jpim

    DO lay = 1, nlayers
      predictors_ad%watervapour(14, lay, i) = predictors_ad%watervapour(14, lay, i) + &
         aux%tr_r(lay, i)**3_jpim * predictors_ad%watervapour(15, lay, i)

      aux_ad%wr(lay, i) = &!aux_ad%wr(lay, i) + &
         sec_wr(lay) * aux%tr_r(lay, i) *      predictors_ad%watervapour(14, lay, i) + &
         sec_wr(lay) * aux%ww_r(lay, i) *      predictors_ad%watervapour(3, lay, i)  + &  !12
         sec_wr_sqrt(lay) * aux%ww_r(lay, i) * predictors_ad%watervapour(8, lay, i) !13

      aux_ad%tr_r(lay, i) = aux_ad%tr_r(lay, i) + &
         sec_wr(lay) * aux%wr(lay, i) * predictors_ad%watervapour(14, lay, i) + &
         3._jprb * aux%tr_r(lay, i)**2_jpim * predictors%watervapour(14, lay, i) * &
                                        predictors_ad%watervapour(15, lay, i)

      predictors_ad%watervapour(13, lay, i) = predictors_ad%watervapour(13, lay, i) + &
         2._jprb * predictors%watervapour(13, lay, i) * predictors_ad%watervapour(12, lay, i)  ! 3

      raytracing_ad%pathsat(lay, i) = raytracing_ad%pathsat(lay, i) + &
         2._jprb * raytracing%pathsat(lay, i) * aux%ww(lay, i) * aux%ww(lay, i) * & !12
                                        predictors_ad%watervapour(13, lay, i)

      aux_ad%ww(lay, i) = &!aux_ad%ww(lay, i) + &
         2._jprb * raytracing%pathsat(lay, i)**2_jpim * aux%ww(lay, i) * &
                                        predictors_ad%watervapour(13, lay, i)  !12

      aux_ad%dt(lay, i) = &!aux_ad%dt(lay, i) + &
         ABS(aux%dt(lay, i)) * sec_wr(lay) * 2._jprb * predictors_ad%watervapour(11, lay, i) + &      ! 10
         sec_wr(lay) *                                 predictors_ad%watervapour(4, lay, i)  + & !4
         sec_wr_sqrt(lay) *                            predictors_ad%watervapour(6, lay, i) !11
   ENDDO

    DO lay = 1, nlayers
      predictors_ad%watervapour(9, lay, i) = predictors_ad%watervapour(9, lay, i) + &
         (4._jprb / 3._jprb) * sec_wr(lay) *           predictors_ad%watervapour(10, lay, i)        ! 9
      predictors_ad%watervapour(5, lay, i) = predictors_ad%watervapour(5, lay, i) + &
         1.5_jprb * sec_wr(lay) *                      predictors_ad%watervapour(9, lay, i)

      sec_wr_ad(lay) = &!sec_wr_ad(lay) + &
                                               predictors_ad%watervapour(1, lay, i) + &  !  7
        aux%wr(lay, i) * aux%ww_r(lay, i) *    predictors_ad%watervapour(3, lay, i) + &  ! 12
        aux%dt(lay, i) *                       predictors_ad%watervapour(4, lay, i) + &  !  4
        2._jprb * sec_wr(lay) *                predictors_ad%watervapour(5, lay, i) + &  !  1
        ABS(aux%dt(lay, i)) * aux%dt(lay, i) * predictors_ad%watervapour(11, lay, i)+ &  ! 10
        aux%wr(lay, i) * aux%tr_r(lay, i) *    predictors_ad%watervapour(14, lay, i)     ! 14

      sec_wr_sqrt_ad(lay) = &!sec_wr_sqrt_ad(lay) + &
                                             predictors_ad%watervapour(2, lay, i) + &                  !  5
        aux%dt(lay, i) *                     predictors_ad%watervapour(6, lay, i) + &!11
        0.5_jprb * SQRT(sec_wr_rsqrt(lay)) * predictors_ad%watervapour(7, lay, i) + &! 6
        aux%wr(lay, i) * aux%ww_r(lay, i) *  predictors_ad%watervapour(8, lay, i) !13

      aux_ad%ww_r(lay, i) = &!aux_ad%ww_r(lay, i) + &
         sec_wr(lay) * aux%wr(lay, i) *      predictors_ad%watervapour(3, lay, i) + &!12
         sec_wr_sqrt(lay) * aux%wr(lay, i) * predictors_ad%watervapour(8, lay, i) ! 13
    ENDDO

    sec_wr_sqrt_ad(:) = sec_wr_sqrt_ad(:) + &
      2._jprb * sec_wr_sqrt(:) * sec_wr_ad(:)

    raytracing_ad%pathsat_sqrt(:, i) = raytracing_ad%pathsat_sqrt(:, i) + &
      sec_wr_sqrt_ad(:) * aux%wr_sqrt(:, i)

    aux_ad%wr_sqrt(:, i) = &!aux_ad%wr_sqrt(:, i) + &
      raytracing%pathsat_sqrt(:, i) * sec_wr_sqrt_ad(:)
  ENDDO

!
! ozone
!---------
! if no input O3 profile, variables or, ow and dto have been set
! to the reference profile values (1, 1, 0)

  IF (coef%nozone > 0) THEN
    DO i = 1, nprofiles
      sec_or_sqrt(:)             = raytracing%pathsat_sqrt(:, i) * aux%or_sqrt(:, i)
      sec_or(:)                  = sec_or_sqrt(:)**2_jpim
      DO lay = 1, prof(1)%nlayers
        predictors_ad%ozone(10, lay, i) = predictors_ad%ozone(10, lay, i) + &
           2._jprb * predictors%ozone(10, lay, i) * predictors_ad%ozone(11, lay, i)

        raytracing_ad%pathsat(lay, i) = raytracing_ad%pathsat(lay, i) + &
          aux%ow(lay, i) *                    predictors_ad%ozone(10, lay, i)

        raytracing_ad%pathsat_sqrt(lay, i) = raytracing_ad%pathsat_sqrt(lay, i) + &
          sec_or(lay) * aux%ow_sqrt(lay, i) * predictors_ad%ozone(9, lay, i)

        aux_ad%ow(lay, i) = &!aux_ad%ow(lay, i) + &
           sec_or(lay) * aux%or(lay, i) * predictors_ad%ozone(6, lay, i) + &
           sec_or(lay) *                  predictors_ad%ozone(8, lay, i) + &
           raytracing%pathsat(lay, i) *   predictors_ad%ozone(10, lay, i)

        sec_or_ad(lay) = &!sec_or_ad(lay) + &
                                                                  predictors_ad%ozone(1, lay, i) + &
          aux%dto(lay, i) *                                       predictors_ad%ozone(3, lay, i) + &
          2._jprb * sec_or(lay) *                                 predictors_ad%ozone(4, lay, i) + &
          aux%or(lay, i) * aux%ow(lay, i) *                       predictors_ad%ozone(6, lay, i) + &
          aux%ow(lay, i) *                                        predictors_ad%ozone(8, lay, i) + &
          raytracing%pathsat_sqrt(lay, i) * aux%ow_sqrt(lay, i) * predictors_ad%ozone(9, lay, i)

        sec_or_sqrt_ad(lay) = &!sec_or_sqrt_ad(lay) + &
                                              predictors_ad%ozone(2, lay, i) + &
          aux%dto(lay, i) *                   predictors_ad%ozone(5, lay, i) + &
          aux%or(lay, i) * aux%ow_r(lay, i) * predictors_ad%ozone(7, lay, i)

        aux_ad%dto(lay, i) = &!aux_ad%dto(lay, i) + &
          sec_or(lay) *      predictors_ad%ozone(3, lay, i) + &
          sec_or_sqrt(lay) * predictors_ad%ozone(5, lay, i)

        aux_ad%or(lay, i) = &!aux_ad%or(lay, i) + &
          sec_or(lay) * aux%ow(lay, i) *        predictors_ad%ozone(6, lay, i) + &
          sec_or_sqrt(lay) * aux%ow_r(lay, i) * predictors_ad%ozone(7, lay, i)

        aux_ad%ow_r(lay, i) = &!aux_ad%ow_r(lay, i) +
          sec_or_sqrt(lay) * aux%or(lay, i) * predictors_ad%ozone(7, lay, i)

        aux_ad%ow_sqrt(lay, i) = &!aux_ad%ow_sqrt(lay, i) + &
           sec_or(lay) * raytracing%pathsat_sqrt(lay, i) * predictors_ad%ozone(9, lay, i)
      ENDDO

      sec_or_sqrt_ad(:) = sec_or_sqrt_ad(:) + &
        2._jprb * sec_or_sqrt(:) * sec_or_ad(:)

      raytracing_ad%pathsat_sqrt(:, i) = raytracing_ad%pathsat_sqrt(:, i) + &
        aux%or_sqrt(:, i) * sec_or_sqrt_ad(:)
      aux_ad%or_sqrt(:, i) = &!aux_ad%or_sqrt(:, i) + &
        raytracing%pathsat_sqrt(:, i) * sec_or_sqrt_ad(:)
    ENDDO
  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_7_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setpredictors_7_ad
