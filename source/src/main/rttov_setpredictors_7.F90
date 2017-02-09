!
SUBROUTINE rttov_setpredictors_7( &
            & opts,       &
            & prof,       &
            & geom,       &
            & coef,       &
            & aux,        &
            & pred_info,  &
            & predictors, &
            & raytracing)
!
! Description
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
!  1.4       16/01/2006  Marco Matricardi (ECMWF):
!               --       Altitude dependent local zenith angle introduced.
!  1.5       22/08/2007  Optimised (D Salmond)
!  1.6       12/08/2008  Zeeman effect based on Yong Han (P. Rayer)
!  1.7       27/02/2009  Profile levels to include ToA. Distinguish arrays
!                        in raytracing (on levels) from all others (on
!                        lays). Predictors prepared to maintain agreement
!                        with the RTTOV-7 scheme in RTTOV-9 (P. Rayer)
!  1.8       02/12/2009  Pathsat, Pathsun and related quantities are now lay arrays
!                        (Marco Matricardi).
!  1.9       17/06/2010  Combined non-Zeeman and Zeeman predictors for SSMIS for
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
       & predictors_type, &
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
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options  ), INTENT(IN)    :: opts
  TYPE(profile_type   ), INTENT(IN)    :: prof(:)         ! profile (ppmv dry)
  TYPE(rttov_coef     ), INTENT(IN)    :: coef            ! coefficients
  TYPE(geometry_type  ), INTENT(IN)    :: geom(SIZE(prof))! geometry
  TYPE(predictors_type), INTENT(INOUT) :: pred_info       ! predictors structure
  TYPE(rttov_path_pred), INTENT(INOUT) :: predictors      ! predictors
  TYPE(profile_aux    ), INTENT(IN)    :: aux             ! auxillary profiles info.
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing
!INTF_END
  INTEGER(KIND=jpim) :: level, lay
  INTEGER(KIND=jpim) :: i
  REAL   (KIND=jprb) :: deltac(prof(1)%nlayers)
  REAL   (KIND=jprb) :: ztemp
  REAL   (KIND=jprb) :: sec_or(prof(1)%nlayers), sec_or_sqrt(prof(1)%nlayers)
  REAL   (KIND=jprb) :: sec_wr(prof(1)%nlayers), sec_wr_sqrt(prof(1)%nlayers)
  INTEGER(KIND=jpim) :: nprofiles, nlayers
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_7', 0_jpim, ZHOOK_HANDLE)
  nprofiles = SIZE(prof)
  nlayers = prof(1)%nlayers

! aux% variables are calculated in rttov_profaux

! mixed gases
!---
  IF (coef%id_inst == inst_id_ssmisz .AND.  coef%IncZeeman) THEN
    ! The ssmisz instrument ID is for SSMIS coef files containing *only* Zeeman channels
    ! so we only calculate the Zeeman predictors (not the "regular" ones).
    ! This is intended for testing rather than for public coefs.
    DO i = 1, nprofiles
      DO lay = 1, nlayers
          predictors%mixedgas(1, lay, i) = raytracing%pathsat(lay, i)
          predictors%mixedgas(2, lay, i) = (300.0_jprb/aux%t_layer(lay, i)) * raytracing%pathsat(lay, i)
          predictors%mixedgas(3, lay, i) = prof(i)%cosbk**2_jpim * raytracing%pathsat(lay, i)
          predictors%mixedgas(4, lay, i) = predictors%mixedgas(2,lay, i) / prof(i)%Be
          predictors%mixedgas(5, lay, i) = predictors%mixedgas(2,lay, i) * prof(i)%cosbk**2_jpim
          predictors%mixedgas(6, lay, i) = raytracing%pathsat(lay, i) / prof(i)%Be
          predictors%mixedgas(7, lay, i) = predictors%mixedgas(6,lay, i) / prof(i)%Be
          predictors%mixedgas(8, lay, i) = prof(i)%Be * raytracing%pathsat(lay, i)
          predictors%mixedgas(9, lay, i) = prof(i)%Be**3_jpim * raytracing%pathsat(lay, i)
          predictors%mixedgas(10, lay, i) = predictors%mixedgas(3,lay, i) * prof(i)%Be
          predictors%mixedgas(11, lay, i) = predictors%mixedgas(10,lay, i) * prof(i)%Be
      ENDDO
    ENDDO

  ELSE
    DO i = 1, nprofiles
      DO lay = 1, nlayers
        ! only effective for non-Zeeman chans 1-18,23-24 - coefficient file will have zeros for chan 19-22
        predictors%mixedgas(1, lay, i)  = raytracing%pathsat(lay, i)
        predictors%mixedgas(2, lay, i)  = raytracing%pathsat(lay, i) ** 2_jpim
        predictors%mixedgas(3, lay, i)  = raytracing%pathsat(lay, i) * aux%tr(lay, i)
        predictors%mixedgas(4, lay, i)  = raytracing%pathsat(lay, i) * aux%tr(lay, i) * aux%tr(lay, i)
        predictors%mixedgas(5, lay, i)  = aux%tr(lay, i)
        predictors%mixedgas(6, lay, i)  = aux%tr(lay, i) * aux%tr(lay, i)!tr*tr
        predictors%mixedgas(7, lay, i)  = raytracing%pathsat(lay, i) * aux%tw(lay, i)
        predictors%mixedgas(8, lay, i)  = raytracing%pathsat(lay, i) * aux%tw(lay, i) * aux%tr_r(lay, i)
        predictors%mixedgas(9, lay, i)  = raytracing%pathsat_sqrt(lay, i)
        predictors%mixedgas(10, lay, i) = raytracing%pathsat_sqrt(lay, i) * aux%tw_4rt(lay, i)
      ENDDO
    ENDDO

    IF (coef%id_inst == inst_id_ssmis .AND. coef%IncZeeman) THEN
      DO i = 1, nprofiles
        DO lay = 1, nlayers
          ! SSMIS with Zeeman coefficient file
          ! geomagnetic field variables (Be, cosbk) are part of the user input

          ! only effective for Zeeman chans 19-22 - coefficient file will have zeros for chan 1-18,23-24
          ! NB require prof(i) % Be > 0. (divisor)
          predictors%mixedgas(11, lay, i) = raytracing%pathsat(lay, i)
          predictors%mixedgas(12, lay, i) = (300.0_jprb/aux%t_layer(lay, i)) * raytracing%pathsat(lay, i)
          predictors%mixedgas(13, lay, i) = prof(i)%cosbk**2_jpim * raytracing%pathsat(lay, i)
          predictors%mixedgas(14, lay, i) = predictors%mixedgas(12,lay, i) / prof(i)%Be
          predictors%mixedgas(15, lay, i) = predictors%mixedgas(12,lay, i) * prof(i)%cosbk**2_jpim
          predictors%mixedgas(16, lay, i) = raytracing%pathsat(lay, i) / prof(i)%Be
          predictors%mixedgas(17, lay, i) = predictors%mixedgas(16,lay, i) / prof(i)%Be
          predictors%mixedgas(18, lay, i) = prof(i)%Be * raytracing%pathsat(lay, i)
          predictors%mixedgas(19, lay, i) = prof(i)%Be**3_jpim * raytracing%pathsat(lay, i)
          predictors%mixedgas(20, lay, i) = predictors%mixedgas(13,lay, i) * prof(i)%Be
          predictors%mixedgas(21, lay, i) = predictors%mixedgas(20,lay, i) * prof(i)%Be
        ENDDO
      ENDDO
    ELSEIF (coef%id_inst == inst_id_amsua .AND. coef%IncZeeman) THEN
      ! AMSU-A with Zeeman coefficient file
      ! only effective for Zeeman chan 14 - coefficient file will have zeros for chan 1-13
      ! NB some of YH's original predictors omitted - effectively duplicated by predictors 1-4 above
      DO i = 1, nprofiles
        DO lay = 1, nlayers
          predictors%mixedgas(11, lay, i) = prof(i)%cosbk ** 2_jpim * raytracing%pathsat(lay, i)
          predictors%mixedgas(12, lay, i) = prof(i)%Be * raytracing%pathsat(lay, i) ** 2_jpim
          predictors%mixedgas(13, lay, i) = prof(i)%Be ** 3_jpim * raytracing%pathsat(lay, i)
          predictors%mixedgas(14, lay, i) = (prof(i)%cosbk * prof(i)%Be * raytracing%pathsat(lay, i)) ** 2_jpim
        ENDDO
      ENDDO
    ENDIF
  ENDIF
! water vapour - numbers in right hand are predictor numbers
! in the reference document for RTTOV7 (science and validation report)
!----------------

  DO i = 1, nprofiles
    sec_wr_sqrt(:) = raytracing%pathsat_sqrt(:, i) * aux%wr_sqrt(:, i)
    sec_wr(:) = sec_wr_sqrt(:)**2_jpim

    DO lay = 1, nlayers
      predictors%watervapour(1, lay, i)  = sec_wr(lay)                                             !  7
      predictors%watervapour(2, lay, i)  = sec_wr_sqrt(lay)                                      !  5
      predictors%watervapour(3, lay, i)  = sec_wr(lay) * aux%wr(lay, i) * aux%ww_r(lay, i)     ! 12
      predictors%watervapour(4, lay, i)  = sec_wr(lay) * aux%dt(lay, i)                              !  4
      predictors%watervapour(5, lay, i)  = sec_wr(lay) * sec_wr(lay)                             !  1
      predictors%watervapour(6, lay, i)  = sec_wr_sqrt(lay) * aux%dt(lay, i)        ! 11
      predictors%watervapour(7, lay, i)  = SQRT(predictors%watervapour(2, lay, i))                 ! 6
      predictors%watervapour(8, lay, i)  = predictors%watervapour(2, lay, i) * aux%wr(lay, i) * aux%ww_r(lay, i) ! 13
      predictors%watervapour(9, lay, i)  = predictors%watervapour(5, lay, i) * sec_wr(lay)       ! 8
      predictors%watervapour(10, lay, i) = predictors%watervapour(9, lay, i) * sec_wr(lay)       ! 9
    ENDDO

    DO lay = 1, nlayers
      predictors%watervapour(11, lay, i) = sec_wr(lay) * aux%dt(lay, i) * ABS(aux%dt(lay, i))          ! 10
      predictors%watervapour(13, lay, i) = (raytracing%pathsat(lay, i) * aux%ww(lay, i)) ** 2_jpim   ! 2
      predictors%watervapour(12, lay, i) = predictors%watervapour(13, lay, i) ** 2_jpim   ! 3    
      predictors%watervapour(14, lay, i) = sec_wr(lay) * aux%wr(lay, i) * aux%tr_r(lay, i)       ! 14
      predictors%watervapour(15, lay, i) = sec_wr(lay) * aux%wr(lay, i) * aux%tr_r(lay, i) ** 4_jpim
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

      DO lay = 1, prof(1)%nlayers
        predictors%ozone(1, lay, i)  = sec_or(lay)
        predictors%ozone(2, lay, i)  = sec_or_sqrt(lay)
        predictors%ozone(3, lay, i)  = sec_or(lay) * aux%dto(lay, i)
        predictors%ozone(4, lay, i)  = sec_or(lay) * sec_or(lay)
        predictors%ozone(5, lay, i)  = sec_or_sqrt(lay) * aux%dto(lay, i)
        predictors%ozone(6, lay, i)  = sec_or(lay) * aux%or(lay, i) * aux%ow(lay, i)
        predictors%ozone(7, lay, i)  = sec_or_sqrt(lay) * aux%or(lay, i) * aux%ow_r(lay, i)
        predictors%ozone(8, lay, i)  = sec_or(lay) * aux%ow(lay, i)
        predictors%ozone(9, lay, i)  = sec_or(lay) * raytracing%pathsat_sqrt(lay, i) * aux%ow_sqrt(lay, i)
        predictors%ozone(10, lay, i) = raytracing%pathsat(lay, i) * aux%ow(lay, i)
        predictors%ozone(11, lay, i) = (raytracing%pathsat(lay, i) * aux%ow(lay, i)) ** 2_jpim
      ENDDO
    ENDDO
  ENDIF

!
! cloud
!---------
  IF (coef%id_sensor == sensor_id_mw) THEN
    IF (opts%rt_mw%clw_Data) THEN

      pred_info%ncloud = 1

      ztemp = 1.0_jprb / (4.3429_jprb * gravity)

      DO i = 1, nprofiles
        DO lay = 1, prof(1)%nlayers
          ! NB in RTTOV-10, lay 1 is the lay above level 2
          level         = lay + 1
          deltac(lay) = 0.1820_jprb * 100.0_jprb * coef%dp(lay) * ztemp
          predictors%clw(lay, i) = deltac(lay) * prof(i)%clw(level) * geom(i)%seczen
        ENDDO
      ENDDO

      DO i = 1, nprofiles
        DO lay = 2, prof(1)%nlayers
          level = lay + 1
          predictors%clw(lay, i) =      &
            & 0.5_jprb * (predictors%clw(lay, i) + deltac(lay) * prof(i)%clw(level - 1) * geom(i)%seczen)
        ENDDO
      ENDDO
    ELSE
      predictors%clw(:,:) = 0._jprb
      pred_info%ncloud = 0_jpim
    ENDIF
  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_7', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setpredictors_7
