
SUBROUTINE rttov_setpredictors_7_k( &
            & opts,          &
            & chanprof,      &
            & profiles,      &
            & profiles_k,    &
            & geom,          &
            & coef,          &
            & aux,           &
            & aux_k,         &
            & predictors,    &
            & predictors_k,  &
            & raytracing,    &
            & raytracing_k)
!
! Description
! K of rttov_setpredictors_7
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
       & rttov_chanprof,  &
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
  TYPE(rttov_chanprof ), INTENT(IN)    :: chanprof(:)
  TYPE(profile_type   ), INTENT(IN)    :: profiles(:)         ! profile (ppmv dry)
  TYPE(profile_type   ), INTENT(INOUT) :: profiles_k(:)
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  TYPE(geometry_type  ), INTENT(IN)    :: geom(:)
  TYPE(rttov_path_pred), INTENT(IN)    :: predictors
  TYPE(rttov_path_pred), INTENT(INOUT) :: predictors_k
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing_k
  TYPE(profile_aux    ), INTENT(IN)    :: aux
  TYPE(profile_aux    ), INTENT(INOUT) :: aux_k
!INTF_END
!local variables:
  INTEGER(KIND=jpim) :: lev, lay, prof, i
  REAL   (KIND=jprb) :: deltac   (profiles(1)%nlayers)

  REAL   (KIND=jprb) :: sec_or_k(profiles(1)%nlayers), sec_or_sqrt_k(profiles(1)%nlayers) 
  REAL   (KIND=jprb) :: sec_wr_sqrt(profiles(1)%nlayers), sec_or_sqrt(profiles(1)%nlayers)
  REAL   (KIND=jprb) :: sec_wr_k(profiles(1)%nlayers), sec_wr_sqrt_k(profiles(1)%nlayers) 

  REAL   (KIND=jprb) :: mgtemp(4,profiles(1)%nlayers), wtemp(16,profiles(1)%nlayers), o3temp(11,profiles(1)%nlayers)
  INTEGER(KIND=jpim) :: nlayers, nchannels
  INTEGER(KIND=jpim) :: last_prof
  REAL   (KIND=jprb) :: ZHOOK_HANDLE

 
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_7_K', 0_jpim, ZHOOK_HANDLE)

  nchannels = SIZE(chanprof)
  nlayers = profiles(1)%nlayers

  IF(aux%on_coef_levels) raytracing_k%pathsat = 0._jprb

! cloud
!---------
  IF (coef%id_sensor == sensor_id_mw .AND. opts%rt_mw%clw_Data) THEN
    DO lay = 1, profiles(1)%nlayers
      deltac(lay) = 0.1820_jprb * 100.0_jprb * coef%dp(lay) / (4.3429_jprb * gravity)
    ENDDO

!CDIR NOLOOPCHG
    DO lay = 2, profiles(1)%nlayers
      lev = lay + 1
      DO i = 1, nchannels
        prof = chanprof(i)%prof
        profiles_k(i)%clw(lev - 1) =      &
          & profiles_k(i)%clw(lev - 1) + 0.5_jprb * predictors_k%clw(lay, i) * deltac(lay) * geom(prof)%seczen
        predictors_k%clw(lay, i)   = 0.5_jprb * predictors_k%clw(lay, i)
      ENDDO
    ENDDO

!CDIR NOLOOPCHG
    DO lay = 1, profiles(1)%nlayers
      lev = lay + 1
      DO i = 1, nchannels
        prof = chanprof(i)%prof
        profiles_k(i)%clw(lev) =      &
          & profiles_k(i)%clw(lev) + predictors_k%clw(lay, i) * deltac(lay) * geom(prof)%seczen
      ENDDO
    ENDDO
  ENDIF


! mixed gases
!---

  IF (coef%id_inst == inst_id_ssmisz .AND.  coef%IncZeeman) THEN

      ! SSMISZ coefficient file - for testing only
      DO i = 1, nchannels
        prof = chanprof(i)%prof
        DO lay = 1, nlayers
          predictors_k%mixedgas(10, lay, i) = predictors_k%mixedgas(10, lay, i) + &
             predictors_k%mixedgas(11, lay, i) * profiles(prof)%Be
          predictors_k%mixedgas(3, lay, i) = predictors_k%mixedgas(3, lay, i) + &
             predictors_k%mixedgas(10, lay, i) * profiles(prof)%Be

          raytracing_k%pathsat(lay, i) = raytracing_k%pathsat(lay, i) + &
             profiles(prof)%Be ** 3_jpim * predictors_k%mixedgas(9, lay, i)
          raytracing_k%pathsat(lay, i) = raytracing_k%pathsat(lay, i) + &
             profiles(prof)%Be * predictors_k%mixedgas(8, lay, i)
          predictors_k%mixedgas(6, lay, i) = predictors_k%mixedgas(6, lay, i) + &
             predictors_k%mixedgas(7, lay, i) / profiles(prof)%Be
          raytracing_k%pathsat(lay, i) = raytracing_k%pathsat(lay, i) + &
             predictors_k%mixedgas(6, lay, i) / profiles(prof)%Be
          predictors_k%mixedgas(2, lay, i) = predictors_k%mixedgas(2, lay, i) + &
             predictors_k%mixedgas(5, lay, i) * profiles(prof)%cosbk ** 2_jpim
          predictors_k%mixedgas(2, lay, i) = predictors_k%mixedgas(2, lay, i) + &
             predictors_k%mixedgas(4, lay, i) / profiles(prof)%Be

          raytracing_k%pathsat(lay, i) = raytracing_k%pathsat(lay, i) + &
             profiles(prof)%cosbk ** 2_jpim * predictors_k%mixedgas(3, lay, i)

          aux_k%t_layer(lay, i) = &!aux_k%t_layer(lay, i)
            -(predictors%mixedgas(2, lay, prof) / &
             aux%t_layer(lay, prof)) * predictors_k%mixedgas(2, lay, i)

          raytracing_k%pathsat(lay, i) = raytracing_k%pathsat(lay, i) + &
             (300.0_jprb / aux%t_layer(lay, prof)) * predictors_k%mixedgas(2, lay, i)
          raytracing_k%pathsat(lay, i) =  raytracing_k%pathsat(lay, i) + predictors_k%mixedgas(1, lay, i)

          ! Initialise some K variables which are initialised in the ELSE clause for other instruments
          aux_k%tr(lay, i) = 0._jprb
          aux_k%tr_r(lay, i) = 0._jprb
          aux_k%tw(lay, i) = 0._jprb
          aux_k%tw_4rt(lay, i) = 0._jprb
          raytracing_k%pathsat_sqrt(lay, i) = 0._jprb
        ENDDO
      ENDDO

  ELSE
    last_prof = -1
    DO i = 1, nchannels
      prof = chanprof(i)%prof
      IF (prof .NE. last_prof) THEN
        DO lay=1, nlayers
          mgtemp(1,lay) = 2._jprb * raytracing%pathsat(lay, prof)
          mgtemp(2,lay) = mgtemp(1,lay) * aux%tr(lay, prof)
          mgtemp(3,lay) = aux%tw(lay, prof) * aux%tr_r(lay, prof)
          mgtemp(4,lay) = raytracing%pathsat(lay, prof) * aux%tw(lay, prof)
        ENDDO
        last_prof = prof
      ENDIF

      DO lay = 1, nlayers
        ! only effective for non-Zeeman chans 1-18,23-24 - coefficient file will have zeros for chan 19-22
        raytracing_k%pathsat(lay, i) = raytracing_k%pathsat(lay, i) + &
                                       predictors_k%mixedgas(1, lay, i) + &
          mgtemp(1,lay)   *            predictors_k%mixedgas(2, lay, i) + &
          aux%tr(lay, prof) *          predictors_k%mixedgas(3, lay, i) + &
          aux%tr(lay, prof) **2_jpim * predictors_k%mixedgas(4, lay, i) + &
          aux%tw(lay, prof) *          predictors_k%mixedgas(7, lay, i) + &
          mgtemp(3,lay)   *            predictors_k%mixedgas(8, lay, i)

        aux_k%tr(lay, i) = &!aux_k%tr(lay, i) + &
          raytracing%pathsat(lay, prof) * predictors_k%mixedgas(3, lay, i) + &
          mgtemp(2,lay) *                 predictors_k%mixedgas(4, lay, i) + &
                                          predictors_k%mixedgas(5, lay, i) + &
          2._jprb * aux%tr(lay, prof) *   predictors_k%mixedgas(6, lay, i)

        aux_k%tw(lay, i) = &!aux_k%tw(lay, i) +
          raytracing%pathsat(lay, prof) * (                   predictors_k%mixedgas(7, lay, i) + &
                                        aux%tr_r(lay, prof) * predictors_k%mixedgas(8, lay, i))

        aux_k%tr_r(lay, i) = &!aux_k%tr_r(lay, i) + &
          mgtemp(4,lay) * predictors_k%mixedgas(8, lay, i)

        raytracing_k%pathsat_sqrt(lay, i) = &!raytracing_k%pathsat_sqrt(lay, i) +
                                  predictors_k%mixedgas(9, lay, i) + &
          aux%tw_4rt(lay, prof) * predictors_k%mixedgas(10, lay, i)

        aux_k%tw_4rt(lay, i) = &!aux_k%tw_4rt(lay, i) + &
          raytracing%pathsat_sqrt(lay, prof) * predictors_k%mixedgas(10, lay, i)
      ENDDO
    ENDDO

    IF (coef%id_inst == inst_id_ssmis .AND. coef%IncZeeman) THEN
      DO i = 1, nchannels
        prof = chanprof(i)%prof
        DO lay = 1, nlayers
          ! SSMIS with Zeeman coefficient file
          ! geomagnetic field variables (Be, cosbk) are part of the user input

          ! only effective for Zeeman chans 19-22 - coefficient file will have zeros for chan 1-18,23-24
          ! NB require prof(i) % Be >0. (divisor)

          ! X11 -> X21
          ! only effective for Zeeman chans 19-22 - coefficient file will have zeros for chan 1-18,23-24
          predictors_k%mixedgas(20, lay, i) = predictors_k%mixedgas(20, lay, i) + &
             predictors_k%mixedgas(21, lay, i) * profiles(prof)%Be
          predictors_k%mixedgas(13, lay, i) = predictors_k%mixedgas(13, lay, i) + &
             predictors_k%mixedgas(20, lay, i) * profiles(prof)%Be

          raytracing_k%pathsat(lay, i) = raytracing_k%pathsat(lay, i) + &
             profiles(prof)%Be ** 3_jpim * predictors_k%mixedgas(19, lay, i)
          raytracing_k%pathsat(lay, i) = raytracing_k%pathsat(lay, i) + &
             profiles(prof)%Be * predictors_k%mixedgas(18, lay, i)
          predictors_k%mixedgas(16, lay, i) = predictors_k%mixedgas(16, lay, i) + &
             predictors_k%mixedgas(17, lay, i) / profiles(prof)%Be
          raytracing_k%pathsat(lay, i) = raytracing_k%pathsat(lay, i) + &
             predictors_k%mixedgas(16, lay, i) / profiles(prof)%Be
          predictors_k%mixedgas(12, lay, i) = predictors_k%mixedgas(12, lay, i) + &
             predictors_k%mixedgas(15, lay, i) * profiles(prof)%cosbk ** 2_jpim
          predictors_k%mixedgas(12, lay, i) = predictors_k%mixedgas(12, lay, i) + &
             predictors_k%mixedgas(14, lay, i) / profiles(prof)%Be

          raytracing_k%pathsat(lay, i) = raytracing_k%pathsat(lay, i) + &
             profiles(prof)%cosbk ** 2_jpim * predictors_k%mixedgas(13, lay, i)

          aux_k%t_layer(lay, i) = &!aux_k%t_layer(lay, i)
            -(predictors%mixedgas(12, lay, prof) / &
             aux%t_layer(lay, prof)) * predictors_k%mixedgas(12, lay, i)

          raytracing_k%pathsat(lay, i) = raytracing_k%pathsat(lay, i) + &
             (300.0_jprb / aux%t_layer(lay, prof)) * predictors_k%mixedgas(12, lay, i)
          raytracing_k%pathsat(lay, i) =  raytracing_k%pathsat(lay, i) + predictors_k%mixedgas(11, lay, i)
        ENDDO
      ENDDO
    ELSEIF (coef%id_inst == inst_id_amsua .AND. coef%IncZeeman) THEN
      ! AMSU-A with Zeeman coefficient file
      ! only effective for Zeeman chan 14 - coefficient file will have zeros for chan 1-13
      ! NB some of YH's original predictors omitted - effectively duplicated by predictors 1-4 above
      DO i = 1, nchannels
        prof = chanprof(i)%prof
        DO lay = 1, nlayers

          raytracing_k%pathsat(lay, i) = raytracing_k%pathsat(lay, i) + &
             profiles(prof)%cosbk ** 2_jpim * predictors_k%mixedgas(11, lay, i)
          raytracing_k%pathsat(lay, i) = raytracing_k%pathsat(lay, i) +     &
            & 2.0_jprb * profiles(prof)%Be * raytracing%pathsat(lay, prof) * predictors_k%mixedgas(12, lay, i)
          raytracing_k%pathsat(lay, i) = raytracing_k%pathsat(lay, i) + &
             profiles(prof)%Be ** 3_jpim * predictors_k%mixedgas(13, lay, i)
          raytracing_k%pathsat(lay, i) = raytracing_k%pathsat(lay, i) + &
            & 2.0_jprb * (profiles(prof)%cosbk * profiles(prof)%Be) ** 2_jpim * &
            raytracing%pathsat(lay, prof) * predictors_k%mixedgas(14, lay, i)
        ENDDO
      ENDDO
    ENDIF
  ENDIF

!--
! water vapour - numbers in right hand are predictor numbers
! in the reference document for RTTOV7 (science and validation report)
!----------------
  last_prof = -1
  DO i = 1, nchannels
    prof = chanprof(i)%prof
    IF (prof .NE. last_prof) THEN
      sec_wr_sqrt(:) = raytracing%pathsat_sqrt(:, prof) * aux%wr_sqrt(:, prof)
      wtemp(1,:) = sec_wr_sqrt(:)
      wtemp(2,:) = wtemp(1,:)**2_jpim

      DO lay=1, nlayers
        wtemp(3,lay) = aux%tr_r(lay, prof)**3_jpim
        wtemp(4,lay) = wtemp(2,lay) * aux%ww_r(lay, prof)
        wtemp(5,lay) = wtemp(1,lay) * aux%ww_r(lay, prof)
        wtemp(6,lay) = wtemp(2,lay) * aux%tr_r(lay, prof)
        wtemp(7,lay) = wtemp(2,lay) * aux%wr(lay, prof)
        wtemp(8,lay) = 2._jprb * raytracing%pathsat(lay, prof) * aux%ww(lay, prof)**2_jpim
        wtemp(9,lay) = 2._jprb * raytracing%pathsat(lay, prof)**2_jpim * aux%ww(lay, prof)
        wtemp(10,lay) = ABS(aux%dt(lay, prof)) * wtemp(2,lay) * 2._jprb
        wtemp(11,lay) = aux%wr(lay, prof) * aux%ww_r(lay, prof)
        wtemp(12,lay) = ABS(aux%dt(lay, prof)) * aux%dt(lay, prof)
        wtemp(13,lay) = aux%wr(lay, prof) * aux%tr_r(lay, prof)
        wtemp(15,lay) = wtemp(1,lay) * aux%wr(lay, prof)
        wtemp(16,lay) = 3._jprb * aux%tr_r(lay, prof)**2_jpim * predictors%watervapour(14, lay, prof)
      ENDDO

      wtemp(14,:) = 0.5_jprb * SQRT(raytracing%pathsat_rsqrt(:, prof) * aux%wr_rsqrt(:, prof))
      last_prof = prof
    ENDIF

    DO lay = 1, nlayers
      predictors_k%watervapour(14, lay, i) = predictors_k%watervapour(14, lay, i) + &
       wtemp(3,lay) * predictors_k%watervapour(15, lay, i)

      aux_k%wr(lay, i) = &!aux_k%wr(lay, i) + &
         wtemp(4,lay) * predictors_k%watervapour(3, lay, i)  + &  !12
         wtemp(5,lay) * predictors_k%watervapour(8, lay, i) + &!13
         wtemp(6,lay) * predictors_k%watervapour(14, lay, i)

      aux_k%tr_r(lay, i) = aux_k%tr_r(lay, i) + &!aux_k%tr_r(lay, i) + &
         wtemp(7,lay)  * predictors_k%watervapour(14, lay, i) + &
         wtemp(16,lay) * predictors_k%watervapour(15, lay, i)

      predictors_k%watervapour(13, lay, i) = predictors_k%watervapour(13, lay, i) + &
         2._jprb * predictors%watervapour(13, lay, prof) * predictors_k%watervapour(12, lay, i)  ! 3

      raytracing_k%pathsat(lay, i) = raytracing_k%pathsat(lay, i) + &
         wtemp(8,lay) * predictors_k%watervapour(13, lay, i) !12

      aux_k%ww(lay, i) = &!aux_k%ww(lay, i) + &
         wtemp(9,lay) * predictors_k%watervapour(13, lay, i)  !12

      aux_k%dt(lay, i) = &!aux_k%dt(lay, i) + &
         wtemp(2,lay) *  predictors_k%watervapour(4, lay, i) + & !4
         wtemp(1,lay) *  predictors_k%watervapour(6, lay, i) + & !11
         wtemp(10,lay) * predictors_k%watervapour(11, lay, i)    ! 10
   ENDDO

    DO lay = 1, nlayers
      predictors_k%watervapour(9, lay, i) = predictors_k%watervapour(9, lay, i) + &
         (4._jprb / 3._jprb) * wtemp(2,lay) *           predictors_k%watervapour(10, lay, i)        ! 9
      predictors_k%watervapour(5, lay, i) = predictors_k%watervapour(5, lay, i) + &
         1.5_jprb * wtemp(2,lay) *                      predictors_k%watervapour(9, lay, i) 

      sec_wr_k(lay) = &!sec_wr_k(lay) + &
                                 predictors_k%watervapour(1, lay, i) + &  !  7 
        wtemp(11,lay) *          predictors_k%watervapour(3, lay, i) + &  ! 12
        aux%dt(lay, prof) *      predictors_k%watervapour(4, lay, i) + &  !  4
        2._jprb * wtemp(2,lay) * predictors_k%watervapour(5, lay, i) + &  !  1
        wtemp(12,lay) *          predictors_k%watervapour(11, lay, i)+ &  ! 10
        wtemp(13,lay) *          predictors_k%watervapour(14, lay, i)     ! 14

      sec_wr_sqrt_k(lay) = &!sec_wr_sqrt_k(lay) + &
                            predictors_k%watervapour(2, lay, i) + &                  !  5 
        aux%dt(lay, prof) * predictors_k%watervapour(6, lay, i) + &!11
        wtemp(14,lay) *     predictors_k%watervapour(7, lay, i) + &! 6
        wtemp(11,lay) *     predictors_k%watervapour(8, lay, i) !13

      aux_k%ww_r(lay, i) = &!aux_k%ww_r(lay, i) + &
         wtemp(7,lay) *  predictors_k%watervapour(3, lay, i) + &!12
         wtemp(15,lay) * predictors_k%watervapour(8, lay, i) ! 13
    ENDDO

    sec_wr_sqrt_k(:) = sec_wr_sqrt_k(:) + &
      2._jprb * sec_wr_sqrt(:) * sec_wr_k(:)

    raytracing_k%pathsat_sqrt(:, i) = raytracing_k%pathsat_sqrt(:, i) + &
      sec_wr_sqrt_k(:) * aux%wr_sqrt(:, prof)

    aux_k%wr_sqrt(:, i) = &!aux_k%wr_sqrt(:, i) + &
      raytracing%pathsat_sqrt(:, prof) * sec_wr_sqrt_k(:) 
  ENDDO

! ozone
!---------
! if no input O3 profile, variables or, ow and dto have been set
! to the reference profile values (1, 1, 0)

  last_prof = -1
  IF (coef%nozone > 0) THEN
    DO i = 1, nchannels
      prof = chanprof(i)%prof
      IF(prof .NE. last_prof) THEN
        sec_or_sqrt(:) = raytracing%pathsat_sqrt(:, prof) * aux%or_sqrt(:, prof)
        o3temp(1,:) = sec_or_sqrt(:)
        o3temp(2,:) = sec_or_sqrt(:)**2_jpim

        DO lay=1, nlayers
          o3temp(3, lay)  = o3temp(2,lay) * aux%ow_sqrt(lay, prof)
          o3temp(4, lay)  = o3temp(2,lay) * aux%or(lay, prof)
          o3temp(5, lay)  = aux%ow(lay, prof) * aux%or(lay, prof)
          o3temp(6, lay)  = raytracing%pathsat_sqrt(lay, prof) * aux%ow_sqrt(lay, prof)
          o3temp(7, lay)  = aux%or(lay, prof) * aux%ow_r(lay, prof)
          o3temp(8, lay)  = o3temp(2,lay) * aux%ow(lay, prof)
          o3temp(9, lay)  = o3temp(1,lay) * aux%ow_r(lay, prof)
          o3temp(10, lay) = o3temp(1,lay) * aux%or(lay, prof)
          o3temp(11, lay) = o3temp(2,lay) * raytracing%pathsat_sqrt(lay, prof)
        ENDDO
        last_prof = prof
      ENDIF

      DO lay = 1, nlayers
        predictors_k%ozone(10, lay, i) = predictors_k%ozone(10, lay, i) + &
           2._jprb * predictors%ozone(10, lay, prof) * predictors_k%ozone(11, lay, i)

        raytracing_k%pathsat(lay, i) = raytracing_k%pathsat(lay, i) + &
          aux%ow(lay, prof) * predictors_k%ozone(10, lay, i)

        raytracing_k%pathsat_sqrt(lay, i) = raytracing_k%pathsat_sqrt(lay, i) + &
          o3temp(3,lay) * predictors_k%ozone(9, lay, i)

        aux_k%ow(lay, i) = &!aux_k%ow(lay, i) + &
           o3temp(4,lay) *                 predictors_k%ozone(6, lay, i) + &
           o3temp(2,lay) *                 predictors_k%ozone(8, lay, i) + & ! * sec_or(lay)
           raytracing%pathsat(lay, prof) * predictors_k%ozone(10, lay, i)

        sec_or_k(lay) = &!sec_or_k(lay) + &
                                    predictors_k%ozone(1, lay, i) + &
          aux%dto(lay, prof) *      predictors_k%ozone(3, lay, i) + &
          2._jprb * o3temp(2,lay) * predictors_k%ozone(4, lay, i) + &
          o3temp(5,lay) *           predictors_k%ozone(6, lay, i) + &
          aux%ow(lay, prof) *       predictors_k%ozone(8, lay, i)+ &
          o3temp(6,lay) *           predictors_k%ozone(9, lay, i)

        sec_or_sqrt_k(lay) = &!sec_or_sqrt_k(lay) + &
                               predictors_k%ozone(2, lay, i) + &
          aux%dto(lay, prof) * predictors_k%ozone(5, lay, i) + &
          o3temp(7,lay)    *   predictors_k%ozone(7, lay, i)

        aux_k%dto(lay, i) = &!aux_k%dto(lay, i) + &
          o3temp(2,lay) * predictors_k%ozone(3, lay, i) + &
          o3temp(1,lay) * predictors_k%ozone(5, lay, i)

        aux_k%or(lay, i) = &!aux_k%or(lay, i) + &
          o3temp(8,lay) * predictors_k%ozone(6, lay, i) + &
          o3temp(9,lay) * predictors_k%ozone(7, lay, i)

        aux_k%ow_r(lay, i) = &!aux_k%ow_r(lay, i) +
          o3temp(10,lay) * predictors_k%ozone(7, lay, i)

        aux_k%ow_sqrt(lay, i) = &!aux_k%ow_sqrt(lay, i) + &
          o3temp(11,lay) * predictors_k%ozone(9, lay, i)
      ENDDO

      sec_or_sqrt_k(:) = sec_or_sqrt_k(:) + &
        2._jprb * sec_or_sqrt(:) * sec_or_k(:)

      raytracing_k%pathsat_sqrt(:, i) = raytracing_k%pathsat_sqrt(:, i) + &
        aux%or_sqrt(:, prof) * sec_or_sqrt_k(:)

      aux_k%or_sqrt(:, i) = &!aux_k%or_sqrt(:, i) + &
        raytracing%pathsat_sqrt(:, prof) * sec_or_sqrt_k(:)
    ENDDO
  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_7_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setpredictors_7_k
