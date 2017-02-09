!
SUBROUTINE rttov_setpredictors_7_tl( &
            & opts,          &
            & prof,          &
            & prof_tl,       &
            & geom,          &
            & coef,          &
            & aux,           &
            & aux_tl,        &
            & predictors,    &
            & predictors_tl, &
            & raytracing,    &
            & raytracing_tl)
!
! Description
! TL of rttov_setpredictors_7
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
!  1.3       08/12/2005  Add surface humidity (R Saunders)
!  1.4       03/03/2006  Marco Matricardi (ECMWF):
!               --       Altitude dependent local zenith angle introduced.
!  1.5       12/08/2008  Zeeman effect based on Yong Han (P. Rayer)
!  1.6       15/07/2009  User defined ToA. Layers distinct from levels (P.Rayer)
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
! Imported Parameters:
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
  TYPE(profile_type   ), INTENT(IN)    :: prof(:)          ! profile (ppmv dry)
  TYPE(profile_type   ), INTENT(IN)    :: prof_tl(SIZE(prof))
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  TYPE(geometry_type  ), INTENT(IN)    :: geom(SIZE(prof))
  TYPE(rttov_path_pred), INTENT(IN)    :: predictors
  TYPE(rttov_path_pred), INTENT(INOUT) :: predictors_tl
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(IN)    :: raytracing_tl
  TYPE(profile_aux    ), INTENT(IN)    :: aux
  TYPE(profile_aux    ), INTENT(IN)    :: aux_tl
!INTF_END
!local variables:
  INTEGER(KIND=jpim) :: level, lay
  INTEGER(KIND=jpim) :: i
  REAL   (KIND=jprb) :: deltac   (prof(1)%nlayers)

  REAL   (KIND=jprb) :: sec_wr   (prof(1)%nlayers), sec_wr_sqrt(prof(1)%nlayers) 
  REAL   (KIND=jprb) :: sec_or   (prof(1)%nlayers), sec_or_sqrt(prof(1)%nlayers) 
  REAL   (KIND=jprb) :: sec_or_tl(prof(1)%nlayers), sec_or_sqrt_tl(prof(1)%nlayers) 
  REAL   (KIND=jprb) :: sec_wr_tl(prof(1)%nlayers), sec_wr_sqrt_tl(prof(1)%nlayers) 
  REAL   (KIND=jprb) :: sec_wr_rsqrt(prof(1)%nlayers)

  INTEGER(KIND=jpim) :: nprofiles, nlayers
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
! layer number agrees with the level number of its upper boundary
! layer N-1 lies between levels N-1 and N
! profile layer quantities
! Direct variables
!CDIR NOLOOPCHG
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_7_TL', 0_jpim, ZHOOK_HANDLE)
  nprofiles = SIZE(prof)
  nlayers = prof(1)%nlayers

! aux% variables are calculated in rttov_profaux

!--
! mixed gases
!---
  IF (coef%id_inst == inst_id_ssmisz .AND.  coef%IncZeeman) THEN

      ! SSMISZ coefficient file - for testing only
      DO i = 1, nprofiles
        DO lay = 1, nlayers
           predictors_tl%mixedgas(1, lay, i)  = raytracing_tl%pathsat(lay, i)
           predictors_tl%mixedgas(2, lay, i)  =  - (predictors%mixedgas(2, lay, i) / &
             aux%t_layer(lay, i)) * aux_tl%t_layer(lay, i) +      &
             & (300.0_jprb / aux%t_layer(lay, i)) * raytracing_tl%pathsat(lay, i)
           predictors_tl%mixedgas(3, lay, i)  = prof(i)%cosbk ** 2_jpim * raytracing_tl%pathsat(lay, i)
           predictors_tl%mixedgas(4, lay, i)  = predictors_tl%mixedgas(2, lay, i) / prof(i)%Be
           predictors_tl%mixedgas(5, lay, i)  = predictors_tl%mixedgas(2, lay, i) * prof(i)%cosbk ** 2_jpim
           predictors_tl%mixedgas(6, lay, i)  = raytracing_tl%pathsat(lay, i) / prof(i)%Be
           predictors_tl%mixedgas(7, lay, i)  = predictors_tl%mixedgas(6, lay, i) / prof(i)%Be
           predictors_tl%mixedgas(8, lay, i)  = prof(i)%Be * raytracing_tl%pathsat(lay, i)
           predictors_tl%mixedgas(9, lay, i)  = prof(i)%Be ** 3_jpim * raytracing_tl%pathsat(lay, i)
           predictors_tl%mixedgas(10, lay, i) = predictors_tl%mixedgas(3, lay, i) * prof(i)%Be
           predictors_tl%mixedgas(11, lay, i) = predictors_tl%mixedgas(10, lay, i) * prof(i)%Be
         ENDDO
      ENDDO

  ELSE
    DO i = 1, nprofiles
      DO lay = 1, nlayers
        ! only effective for non-Zeeman chans 1-18,23-24 - coefficient file will have zeros for chan 19-22
        predictors_tl%mixedgas(1, lay, i)  = raytracing_tl%pathsat(lay, i)
        predictors_tl%mixedgas(2, lay, i)  = 2._jprb * raytracing_tl%pathsat(lay, i) * raytracing%pathsat(lay, i)
        predictors_tl%mixedgas(3, lay, i)  = raytracing_tl%pathsat(lay, i) * aux%tr(lay, i) + &
                                             raytracing%pathsat(lay, i) * aux_tl%tr(lay, i)
        predictors_tl%mixedgas(4, lay, i)  = aux%tr(lay, i) * (&
                                                               raytracing_tl%pathsat(lay, i) * aux%tr(lay, i) + &
                                                               2._jprb * raytracing%pathsat(lay, i) * aux_tl%tr(lay, i))
        predictors_tl%mixedgas(5, lay, i)  = aux_tl%tr(lay, i)
        predictors_tl%mixedgas(6, lay, i)  = 2._jprb * aux%tr(lay, i) * aux_tl%tr(lay, i) !tr*tr
        predictors_tl%mixedgas(7, lay, i)  = raytracing_tl%pathsat(lay, i) * aux%tw(lay, i) + &
                                             raytracing%pathsat(lay, i) * aux_tl%tw(lay, i)
        predictors_tl%mixedgas(8, lay, i)  = raytracing_tl%pathsat(lay, i) * aux%tw(lay, i) * aux%tr_r(lay, i) + &
                                             raytracing%pathsat(lay, i) * (aux_tl%tw(lay, i) * aux%tr_r(lay, i) + &
                                                                           aux%tw(lay, i) * aux_tl%tr_r(lay, i))
        predictors_tl%mixedgas(9, lay, i)  = raytracing_tl%pathsat_sqrt(lay, i)
        predictors_tl%mixedgas(10, lay, i) = raytracing_tl%pathsat_sqrt(lay, i) * aux%tw_4rt(lay, i) + &
                                             raytracing%pathsat_sqrt(lay, i) * aux_tl%tw_4rt(lay, i)
      ENDDO
    ENDDO

    IF (coef%id_inst == inst_id_ssmis .AND. coef%IncZeeman) THEN
      DO i = 1, nprofiles
        DO lay = 1, nlayers
           ! SSMIS with Zeeman coefficient file
           ! geomagnetic field variables (Be, cosbk) are part of the user input

           ! only effective for Zeeman chans 19-22 - coefficient file will have zeros for chan 1-18,23-24
           ! NB require prof(i) % Be > 0. (divisor)

           ! X11 -> X21
           ! only effective for Zeeman chans 19-22 - coefficient file will have zeros for chan 1-18,23-24
           predictors_tl%mixedgas(11, lay, i)  = raytracing_tl%pathsat(lay, i)
           predictors_tl%mixedgas(12, lay, i)  =  - (predictors%mixedgas(12, lay, i) / &
             aux%t_layer(lay, i)) * aux_tl%t_layer(lay, i) +      &
             & (300.0_jprb / aux%t_layer(lay, i)) * raytracing_tl%pathsat(lay, i)
           predictors_tl%mixedgas(13, lay, i)  = prof(i)%cosbk ** 2_jpim * raytracing_tl%pathsat(lay, i)
           predictors_tl%mixedgas(14, lay, i)  = predictors_tl%mixedgas(12, lay, i) / prof(i)%Be
           predictors_tl%mixedgas(15, lay, i)  = predictors_tl%mixedgas(12, lay, i) * prof(i)%cosbk ** 2_jpim
           predictors_tl%mixedgas(16, lay, i)  = raytracing_tl%pathsat(lay, i) / prof(i)%Be
           predictors_tl%mixedgas(17, lay, i)  = predictors_tl%mixedgas(16, lay, i) / prof(i)%Be
           predictors_tl%mixedgas(18, lay, i)  = prof(i)%Be * raytracing_tl%pathsat(lay, i)
           predictors_tl%mixedgas(19, lay, i)  = prof(i)%Be ** 3_jpim * raytracing_tl%pathsat(lay, i)
           predictors_tl%mixedgas(20, lay, i) = predictors_tl%mixedgas(13, lay, i) * prof(i)%Be
           predictors_tl%mixedgas(21, lay, i) = predictors_tl%mixedgas(20, lay, i) * prof(i)%Be
         ENDDO
      ENDDO
    ELSEIF (coef%id_inst == inst_id_amsua .AND. coef%IncZeeman) THEN
      ! AMSU-A with Zeeman coefficient file
      ! only effective for Zeeman chan 14 - coefficient file will have zeros for chan 1-13
      ! NB some of YH's original predictors omitted - effectively duplicated by predictors 1-4 above
      DO i = 1, nprofiles
        DO lay = 1, nlayers
          predictors_tl%mixedgas(11, lay, i) = prof(i)%cosbk ** 2_jpim * raytracing_tl%pathsat(lay, i)
          predictors_tl%mixedgas(12, lay, i) =      &
            & 2.0_jprb * prof(i)%Be * raytracing%pathsat(lay, i) * raytracing_tl%pathsat(lay, i)
          predictors_tl%mixedgas(13, lay, i) = prof(i)%Be ** 3_jpim * raytracing_tl%pathsat(lay, i)
          predictors_tl%mixedgas(14, lay, i) =      &
            & 2.0_jprb * (prof(i)%cosbk * prof(i)%Be) ** 2_jpim * raytracing%pathsat(lay, i) * raytracing_tl%pathsat(lay, i)
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
    sec_wr_sqrt_tl(:) = raytracing_tl%pathsat_sqrt(:, i) * aux%wr_sqrt(:, i) + &
                        raytracing%pathsat_sqrt(:, i) * aux_tl%wr_sqrt(:, i)

    sec_wr_tl(:) = 2._jprb * sec_wr_sqrt(:) * sec_wr_sqrt_tl(:)

    DO lay = 1, nlayers
      predictors_tl%watervapour(1, lay, i)  = sec_wr_tl(lay)                                             !  7
      predictors_tl%watervapour(2, lay, i)  = sec_wr_sqrt_tl(lay)                                      !  5
      predictors_tl%watervapour(3, lay, i)  = sec_wr_tl(lay) * aux%wr(lay, i) * aux%ww_r(lay, i) + &
                                              sec_wr(lay) * (aux_tl%wr(lay, i) * aux%ww_r(lay, i) + &
                                                             aux%wr(lay, i) * aux_tl%ww_r(lay, i)) ! 12

      predictors_tl%watervapour(4, lay, i)  = sec_wr_tl(lay) * aux%dt(lay, i) + sec_wr(lay) * aux_tl%dt(lay, i)   !4
      predictors_tl%watervapour(5, lay, i)  = 2._jprb * sec_wr_tl(lay) * sec_wr(lay)                             !  1
      predictors_tl%watervapour(6, lay, i)  = sec_wr_sqrt_tl(lay) * aux%dt(lay, i) + &
                                              sec_wr_sqrt(lay) * aux_tl%dt(lay, i) !11

      predictors_tl%watervapour(7, lay, i)  = 0.5_jprb * SQRT(sec_wr_rsqrt(lay)) * sec_wr_sqrt_tl(lay)              ! 6

      predictors_tl%watervapour(8, lay, i)  = sec_wr_sqrt_tl(lay) * aux%wr(lay, i) * aux%ww_r(lay, i) + &
                                              sec_wr_sqrt(lay) * (aux_tl%wr(lay, i) * aux%ww_r(lay, i) + &
                                                                  aux%wr(lay, i) * aux_tl%ww_r(lay, i)) ! 13
      predictors_tl%watervapour(9, lay, i)  = 1.5_jprb * predictors_tl%watervapour(5, lay, i) * sec_wr(lay)
      predictors_tl%watervapour(10, lay, i) = (4._jprb / 3._jprb) * predictors_tl%watervapour(9, lay, i) * sec_wr(lay)       ! 9
    ENDDO

    DO lay = 1, nlayers
      predictors_tl%watervapour(11, lay, i) = ABS(aux%dt(lay, i)) * (&
                                                                     sec_wr_tl(lay) * aux%dt(lay, i) + &
                                                                     sec_wr(lay) * 2._jprb * aux_tl%dt(lay, i))     ! 10
      predictors_tl%watervapour(13, lay, i) = 2._jprb * (raytracing%pathsat(lay, i) * aux%ww(lay, i)) * &  ! 2
                                              (raytracing_tl%pathsat(lay, i) * aux%ww(lay, i) + &
                                               raytracing%pathsat(lay, i) * aux_tl%ww(lay, i))
      predictors_tl%watervapour(12, lay, i) = 2._jprb * predictors%watervapour(13, lay, i) * &
                                              predictors_tl%watervapour(13, lay, i) ! 3
      predictors_tl%watervapour(14, lay, i) = sec_wr_tl(lay) * aux%wr(lay, i) * aux%tr_r(lay, i) + &      ! 14
                                              sec_wr(lay) * (aux_tl%wr(lay, i) * aux%tr_r(lay, i) + &
                                                             aux%wr(lay, i) * aux_tl%tr_r(lay, i))

      predictors_tl%watervapour(15, lay, i) = aux%tr_r(lay, i)**2_jpim * ( &
                                                 3._jprb * predictors%watervapour(14, lay, i) * aux_tl%tr_r(lay, i) + &
                                                 predictors_tl%watervapour(14, lay, i) * aux%tr_r(lay, i))
    ENDDO
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
      sec_or_sqrt_tl(:)          = raytracing_tl%pathsat_sqrt(:, i) * aux%or_sqrt(:, i) + &
                                   raytracing%pathsat_sqrt(:, i) * aux_tl%or_sqrt(:, i)
      sec_or_tl(:)                  = 2._jprb * sec_or_sqrt(:) * sec_or_sqrt_tl(:)

      DO lay = 1, prof(1)%nlayers
        predictors_tl%ozone(1, lay, i)  = sec_or_tl(lay)
        predictors_tl%ozone(2, lay, i)  = sec_or_sqrt_tl(lay)
        predictors_tl%ozone(3, lay, i)  = sec_or_tl(lay) * aux%dto(lay, i) + sec_or(lay) * aux_tl%dto(lay, i)
        predictors_tl%ozone(4, lay, i)  = 2._jprb * sec_or(lay) * sec_or_tl(lay)
        predictors_tl%ozone(5, lay, i)  = sec_or_sqrt_tl(lay) * aux%dto(lay, i) + sec_or_sqrt(lay) * aux_tl%dto(lay, i)
        predictors_tl%ozone(6, lay, i)  = sec_or_tl(lay) * aux%or(lay, i) * aux%ow(lay, i) + &
                                          sec_or(lay) * (aux_tl%or(lay, i) * aux%ow(lay, i) + &
                                                         aux%or(lay, i) * aux_tl%ow(lay, i))
        predictors_tl%ozone(7, lay, i)  = sec_or_sqrt_tl(lay) * aux%or(lay, i) * aux%ow_r(lay, i) + &
                                          sec_or_sqrt(lay) * (aux_tl%or(lay, i) * aux%ow_r(lay, i) + &
                                                              aux%or(lay, i) * aux_tl%ow_r(lay, i))
        predictors_tl%ozone(8, lay, i)  = sec_or_tl(lay) * aux%ow(lay, i) + sec_or(lay) * aux_tl%ow(lay, i)
        predictors_tl%ozone(9, lay, i)  = sec_or_tl(lay) * raytracing%pathsat_sqrt(lay, i) * aux%ow_sqrt(lay, i) + &
                                          sec_or(lay) * (raytracing_tl%pathsat_sqrt(lay, i) * aux%ow_sqrt(lay, i) + &
                                                         raytracing%pathsat_sqrt(lay, i) * aux_tl%ow_sqrt(lay, i))
        predictors_tl%ozone(10, lay, i) = raytracing_tl%pathsat(lay, i) * aux%ow(lay, i) + &
                                          raytracing%pathsat(lay, i) * aux_tl%ow(lay, i)
        predictors_tl%ozone(11, lay, i) = 2._jprb * predictors%ozone(10, lay, i) * predictors_tl%ozone(10, lay, i)
      ENDDO
    ENDDO
  ENDIF

! cloud
!---------
!CDIR NOLOOPCHG
  IF (coef%id_sensor == sensor_id_mw) THEN
    IF (opts%rt_mw%clw_Data) THEN

      DO lay = 1, prof(1)%nlayers
        level         = lay + 1
        deltac(lay) = 0.1820_jprb * 100.0_jprb * coef%dp(lay) / (4.3429_jprb * gravity)
        DO i = 1, nprofiles
          predictors_tl%clw(lay, i) = deltac(lay) * prof_tl(i)%clw(level) * geom(i)%seczen
        ENDDO
      ENDDO

!CDIR NOLOOPCHG

      DO lay = 2, prof_tl(1)%nlayers
        level = lay + 1
        DO i = 1, nprofiles
          predictors_tl%clw(lay, i) =      &
            & 0.5_jprb * (predictors_tl%clw(lay, i) + deltac(lay) * prof_tl(i)%clw(level - 1) * geom(i)%seczen)
        ENDDO
      ENDDO

    ELSE
      predictors_tl%clw(:,:) = 0._jprb
    ENDIF
  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_7_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setpredictors_7_tl
