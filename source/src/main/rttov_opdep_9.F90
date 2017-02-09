!
SUBROUTINE rttov_opdep_9( &
            & nlayers,    &
            & chanprof,   &
            & chanflag,   &
            & pred_info,  &
            & predictors, &
            & aux,        &
            & coef,       &
            & fast_coef,  &
            & opdp_path,  &
            & opdp_ref)
!
! Description:
!  To calculate optical depths for a number of channels
!  and profiles from every pressure level to space.
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
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0    01/06/2005  Marco Matricardi (ECMWF):
!            --       New routine based on rttov_transmit.F90.
!                     Variable trace gases, CO2, N2O,CO and CH4
!                     have been introduced for AIRS and IASI
!  1.1    29/01/2007  Removed polarisation R Saunders
!  1.2    07/04/2008  Marco Matricardi (ECMWF):
!                     modified the routine to be consistent
!                     with the LBLRTM coefficients
!  1.3    27/02/2009  Profile levels to include ToA. Distinguish between
!                     layer arrays and level arrays - size, index labels,
!                     looping (P. Rayer)
!  1.4    18/06/2009  Use rttov9 constant intervals from rttov_const
!                     and change test on upper bounds from <= to <
!                     Philippe Marguinaud
!  1.5    03/11/2009  Transmittances / optical depths on levels (A Geer)
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
  USE rttov_types, ONLY : &
       & rttov_chanprof,  &
       & rttov_coef,      &
       & rttov_fast_coef, &
       & predictors_Type, &
       & rttov_path_pred, &
       & profile_aux
  USE parkind1, ONLY : jpim, jprb, jplm
!INTF_OFF
  USE rttov_const, ONLY :  &
       & sensor_id_mw,     &
       & rttov9_wv0690_50, &
       & rttov9_wv1050_00, &
       & rttov9_wv1095_25, &
       & rttov9_wv1100_25, &
       & rttov9_wv1350_25, &
       & rttov9_wv1750_25, &
       & rttov9_wv1900_25, &
       & rttov9_wv1995_00, &
       & rttov9_wv2000_00, &
       & rttov9_wv2250_00, &
       & rttov9_wv2295_25, &
       & rttov9_wv2380_25
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  INTEGER(KIND=jpim)   , INTENT(IN)    :: nlayers                          ! Number of pressure levels
  TYPE(rttov_chanprof ), INTENT(IN)    :: chanprof(:)                      ! Channel indices
  LOGICAL(KIND=jplm)   , INTENT(IN)    :: chanflag(SIZE(chanprof))         ! Channel calculation flag
  TYPE(predictors_Type), INTENT(IN)    :: pred_info                        ! Predictor structure
  TYPE(rttov_path_pred), INTENT(IN)    :: predictors                       ! Predictors
  TYPE(rttov_coef     ), INTENT(IN)    :: coef                             ! Coefficient structure
  TYPE(rttov_fast_coef), INTENT(IN)    :: fast_coef                        ! Fast coefficients
  TYPE(profile_aux    ), INTENT(IN)    :: aux                              ! auxillary profiles informations
  REAL(KIND=jprb)      , INTENT(INOUT) :: opdp_path(:,:)                   ! optical depths
  REAL(KIND=jprb)      , INTENT(INOUT) :: opdp_ref(nlayers, SIZE(chanprof))! layer optical depth
!INTF_END
!local variables:
  REAL(KIND=jprb) :: opticaldepth(nlayers, SIZE(chanprof)) ! raw layer optical depth
  REAL(KIND=jprb) :: chanx
  REAL(KIND=jprb), POINTER :: debye_prof(:, :) ! pointer on Debye profiles
  INTEGER(KIND=jpim) :: lev, lay, chan, j, prof, nlevels
! cloud liquid water local variables
  REAL   (KIND=jprb) :: zf, zf_sq, z34_dif, z45_dif, z1_sq, z2_sq, z1_div, z2_div
  REAL   (KIND=jprb) :: z1_den, z2_den, zastar, z1_prod, z2_prod, z3_prod, z4_prod
  REAL   (KIND=jprb) :: zbstar, zbstar_sq, za2star, za2star_sq, zdiv, zgstar
  REAL   (KIND=jprb) :: z1f_sq_z1_sq, z2f_sq_z2_sq
  INTEGER(KIND=jpim) :: ii
  REAL   (KIND=jprb) :: t1, t2
  INTEGER(KIND=jpim) :: nchannels
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_9', 0_jpim, ZHOOK_HANDLE)
  nchannels = SIZE(chanprof)
  nlevels   = nlayers + 1
  opticaldepth(:,:) = 0._jprb
!-----------------------------------------
!1. calculate layer gaseous optical depths
!-----------------------------------------

! Summary of which predictors are used (non-solar or solar):
!            WV < 2000.00 Non-solar predictors always used
! 2000.00 <= WV < 2250.00 Solar predictors always used EXCEPT for H2O where it depends on the value of ss_val_chn
! 2250.00 <= WV < 2380.25 Non-solar predictors always used
! 2380.25 <= WV           Solar predictors always used

!--------------------------
!1.1 start with mixed gases
!--------------------------
  DO j = 1, nchannels
    IF (chanflag(j)) THEN
      chan  = chanprof(j)%chan
      chanx = coef%ff_cwn(chan)
      prof  = chanprof(j)%prof
      IF ((chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2380_25) .OR. &
          (chanx < rttov9_wv2000_00)) THEN
        DO lay = 1, nlayers
          t1 = fast_coef%mixedgas(lay, chan, 1) * predictors%mixedgas(1, lay, prof) +      &
            & fast_coef%mixedgas(lay, chan, 2) * predictors%mixedgas(2, lay, prof) +       &
            & fast_coef%mixedgas(lay, chan, 3) * predictors%mixedgas(3, lay, prof) +       &
            & fast_coef%mixedgas(lay, chan, 4) * predictors%mixedgas(4, lay, prof)
          t2 = fast_coef%mixedgas(lay, chan, 5) * predictors%mixedgas(5, lay, prof) +      &
            & fast_coef%mixedgas(lay, chan, 6) * predictors%mixedgas(6, lay, prof) +       &
            & fast_coef%mixedgas(lay, chan, 7) * predictors%mixedgas(7, lay, prof) +       &
            & fast_coef%mixedgas(lay, chan, 8) * predictors%mixedgas(8, lay, prof)
          opticaldepth(lay, j) = t1 + t2
        ENDDO
      ELSE
        DO lay = 1, nlayers
          t1 = fast_coef%mixedgas(lay, chan, 1) * predictors%mixedgas(1, lay, prof) +      &
            & fast_coef%mixedgas(lay, chan, 2) * predictors%mixedgas(2, lay, prof) +       &
            & fast_coef%mixedgas(lay, chan, 3) * predictors%mixedgas(3, lay, prof) +       &
            & fast_coef%mixedgas(lay, chan, 4) * predictors%mixedgas(4, lay, prof) +       &
            & fast_coef%mixedgas(lay, chan, 5) * predictors%mixedgas(5, lay, prof)
          t2 = fast_coef%mixedgas(lay, chan, 6) * predictors%mixedgas(6, lay, prof) +      &
            & fast_coef%mixedgas(lay, chan, 7) * predictors%mixedgas(7, lay, prof) +       &
            & fast_coef%mixedgas(lay, chan, 8) * predictors%mixedgas(8, lay, prof) +       &
            & fast_coef%mixedgas(lay, chan, 9) * predictors%mixedgas(9, lay, prof) +       &
            & fast_coef%mixedgas(lay, chan, 10) * predictors%mixedgas(10, lay, prof)
          opticaldepth(lay, j) = t1 + t2
        ENDDO
      ENDIF
    ENDIF
  ENDDO
!--------------------
!1.2 add water vapour
!--------------------
  DO j = 1, nchannels
    IF (chanflag(j)) THEN
      chan  = chanprof(j)%chan
      chanx = coef%ff_cwn(chan)
      prof  = chanprof(j)%prof
      IF (chanx >= rttov9_wv1095_25 .AND. chanx < rttov9_wv1750_25) THEN
        DO lay = 1, nlayers
          t1 = fast_coef%watervapour(lay, chan, 1) * predictors%watervapour(1, lay, prof) +      &
            & fast_coef%watervapour(lay, chan, 2) * predictors%watervapour(2, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 3) * predictors%watervapour(3, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 4) * predictors%watervapour(4, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 5) * predictors%watervapour(5, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 6) * predictors%watervapour(6, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 7) * predictors%watervapour(7, lay, prof)
          t2 = fast_coef%watervapour(lay, chan, 8) * predictors%watervapour(8, lay, prof) +      &
            & fast_coef%watervapour(lay, chan, 9) * predictors%watervapour(9, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 10) * predictors%watervapour(10, lay, prof) +     &
            & fast_coef%watervapour(lay, chan, 11) * predictors%watervapour(11, lay, prof) +     &
            & fast_coef%watervapour(lay, chan, 12) * predictors%watervapour(12, lay, prof) +     &
            & fast_coef%watervapour(lay, chan, 13) * predictors%watervapour(13, lay, prof)
          IF (coef%nch4 > 0) THEN
            t2 = t2 + fast_coef%watervapour(lay, chan, 14) * predictors%ch4(1, lay, prof) +      &
              & fast_coef%watervapour(lay, chan, 15) * predictors%ch4(1, lay, prof) * predictors%ch4(3, lay, prof)
          ENDIF
          opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
        ENDDO
      ELSE IF (chanx >= rttov9_wv1900_25 .AND. chanx < rttov9_wv1995_00) THEN
        DO lay = 1, nlayers
          t1 = fast_coef%watervapour(lay, chan, 1) * predictors%watervapour(1, lay, prof) +      &
            & fast_coef%watervapour(lay, chan, 2) * predictors%watervapour(2, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 3) * predictors%watervapour(3, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 4) * predictors%watervapour(4, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 5) * predictors%watervapour(5, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 6) * predictors%watervapour(6, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 7) * predictors%watervapour(7, lay, prof)
          t2 = fast_coef%watervapour(lay, chan, 8) * predictors%watervapour(8, lay, prof) +      &
            & fast_coef%watervapour(lay, chan, 9) * predictors%watervapour(9, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 10) * predictors%watervapour(10, lay, prof) +     &
            & fast_coef%watervapour(lay, chan, 11) * predictors%watervapour(11, lay, prof) +     &
            & fast_coef%watervapour(lay, chan, 12) * predictors%watervapour(12, lay, prof) +     &
            & fast_coef%watervapour(lay, chan, 13) * predictors%watervapour(13, lay, prof)
          IF (coef%nco2 > 0) THEN
            t2 = t2 + fast_coef%watervapour(lay, chan, 14) * predictors%co2(1, lay, prof)
          ENDIF
          opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
        ENDDO
      ELSE IF (chanx >= rttov9_wv1995_00 .AND. chanx < rttov9_wv2000_00) THEN
        DO lay = 1, nlayers
          t1 = fast_coef%watervapour(lay, chan, 1) * predictors%watervapour(1, lay, prof) +      &
            & fast_coef%watervapour(lay, chan, 2) * predictors%watervapour(2, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 3) * predictors%watervapour(3, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 4) * predictors%watervapour(4, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 5) * predictors%watervapour(5, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 6) * predictors%watervapour(6, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 7) * predictors%watervapour(7, lay, prof)
          t2 = fast_coef%watervapour(lay, chan, 8) * predictors%watervapour(8, lay, prof) +      &
            & fast_coef%watervapour(lay, chan, 9) * predictors%watervapour(9, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 10) * predictors%watervapour(10, lay, prof) +     &
            & fast_coef%watervapour(lay, chan, 11) * predictors%watervapour(11, lay, prof) +     &
            & fast_coef%watervapour(lay, chan, 12) * predictors%watervapour(12, lay, prof) +     &
            & fast_coef%watervapour(lay, chan, 13) * predictors%watervapour(13, lay, prof)
          IF (coef%nco2 > 0) THEN
            t2 = t2 + fast_coef%watervapour(lay, chan, 14) * predictors%co2(1, lay, prof)
          ENDIF
          IF (coef%nco > 0) THEN
            t2 = t2 + fast_coef%watervapour(lay, chan, 15) * predictors%co(1, lay, prof)
          ENDIF
          opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
        ENDDO
      ELSE IF (chanx < rttov9_wv2000_00) THEN
        DO lay = 1, nlayers
          t1 = fast_coef%watervapour(lay, chan, 1) * predictors%watervapour(1, lay, prof) +      &
            & fast_coef%watervapour(lay, chan, 2) * predictors%watervapour(2, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 3) * predictors%watervapour(3, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 4) * predictors%watervapour(4, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 5) * predictors%watervapour(5, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 6) * predictors%watervapour(6, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 7) * predictors%watervapour(7, lay, prof)
          t2 = fast_coef%watervapour(lay, chan, 8) * predictors%watervapour(8, lay, prof) +      &
            & fast_coef%watervapour(lay, chan, 9) * predictors%watervapour(9, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 10) * predictors%watervapour(10, lay, prof) +     &
            & fast_coef%watervapour(lay, chan, 11) * predictors%watervapour(11, lay, prof) +     &
            & fast_coef%watervapour(lay, chan, 12) * predictors%watervapour(12, lay, prof) +     &
            & fast_coef%watervapour(lay, chan, 13) * predictors%watervapour(13, lay, prof)
          opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
        ENDDO
      ELSE IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2250_00) THEN
        IF (coef%ss_val_chn(chan) > 0) THEN
          DO lay = 1, nlayers
            t1 = fast_coef%watervapour(lay, chan, 1) * predictors%watervapour(1, lay, prof) +      &
              & fast_coef%watervapour(lay, chan, 2) * predictors%watervapour(2, lay, prof) +       &
              & fast_coef%watervapour(lay, chan, 3) * predictors%watervapour(3, lay, prof) +       &
              & fast_coef%watervapour(lay, chan, 4) * predictors%watervapour(4, lay, prof) +       &
              & fast_coef%watervapour(lay, chan, 5) * predictors%watervapour(5, lay, prof) +       &
              & fast_coef%watervapour(lay, chan, 6) * predictors%watervapour(6, lay, prof) +       &
              & fast_coef%watervapour(lay, chan, 7) * predictors%watervapour(7, lay, prof)
            t2 = fast_coef%watervapour(lay, chan, 8) * predictors%watervapour(14, lay, prof) +      &
              & fast_coef%watervapour(lay, chan, 9) * predictors%watervapour(15, lay, prof) +       &
              & fast_coef%watervapour(lay, chan, 10) * predictors%watervapour(10, lay, prof) +      &
              & fast_coef%watervapour(lay, chan, 11) * predictors%watervapour(11, lay, prof) +      &
              & fast_coef%watervapour(lay, chan, 12) * predictors%watervapour(16, lay, prof) +      &
              & fast_coef%watervapour(lay, chan, 13) * predictors%watervapour(17, lay, prof)
            IF (coef%nco2 > 0) THEN
              t2 = t2 + fast_coef%watervapour(lay, chan, 14) * predictors%co2(1, lay, prof)
            ENDIF
            IF (coef%nco > 0) THEN
              t2 = t2 + fast_coef%watervapour(lay, chan, 15) * predictors%co(1, lay, prof)
            ENDIF
            opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
          ENDDO
        ELSE
          DO lay = 1, nlayers
            t1 = fast_coef%watervapour(lay, chan, 1) * predictors%watervapour(1, lay, prof) +      &
              & fast_coef%watervapour(lay, chan, 2) * predictors%watervapour(2, lay, prof) +       &
              & fast_coef%watervapour(lay, chan, 3) * predictors%watervapour(3, lay, prof) +       &
              & fast_coef%watervapour(lay, chan, 4) * predictors%watervapour(4, lay, prof) +       &
              & fast_coef%watervapour(lay, chan, 5) * predictors%watervapour(5, lay, prof) +       &
              & fast_coef%watervapour(lay, chan, 6) * predictors%watervapour(6, lay, prof) +       &
              & fast_coef%watervapour(lay, chan, 7) * predictors%watervapour(7, lay, prof)
            t2 = fast_coef%watervapour(lay, chan, 8) * predictors%watervapour(8, lay, prof) +      &
              & fast_coef%watervapour(lay, chan, 9) * predictors%watervapour(9, lay, prof) +       &
              & fast_coef%watervapour(lay, chan, 10) * predictors%watervapour(10, lay, prof) +     &
              & fast_coef%watervapour(lay, chan, 11) * predictors%watervapour(11, lay, prof) +     &
              & fast_coef%watervapour(lay, chan, 12) * predictors%watervapour(12, lay, prof) +     &
              & fast_coef%watervapour(lay, chan, 13) * predictors%watervapour(13, lay, prof)
            IF (coef%nco2 > 0) THEN
              t2 = t2 + fast_coef%watervapour(lay, chan, 14) * predictors%co2(1, lay, prof)
            ENDIF
            IF (coef%nco > 0) THEN
              t2 = t2 + fast_coef%watervapour(lay, chan, 15) * predictors%co(1, lay, prof)
            ENDIF
            opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
          ENDDO
        ENDIF
      ELSE IF (chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2295_25) THEN
        DO lay = 1, nlayers
          t1 = fast_coef%watervapour(lay, chan, 1) * predictors%watervapour(1, lay, prof) +      &
            & fast_coef%watervapour(lay, chan, 2) * predictors%watervapour(2, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 3) * predictors%watervapour(3, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 4) * predictors%watervapour(4, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 5) * predictors%watervapour(5, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 6) * predictors%watervapour(6, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 7) * predictors%watervapour(7, lay, prof)
          t2 = fast_coef%watervapour(lay, chan, 8) * predictors%watervapour(8, lay, prof) +      &
            & fast_coef%watervapour(lay, chan, 9) * predictors%watervapour(9, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 10) * predictors%watervapour(10, lay, prof) +     &
            & fast_coef%watervapour(lay, chan, 11) * predictors%watervapour(11, lay, prof) +     &
            & fast_coef%watervapour(lay, chan, 12) * predictors%watervapour(12, lay, prof) +     &
            & fast_coef%watervapour(lay, chan, 13) * predictors%watervapour(13, lay, prof)
          IF (coef%nco2 > 0) THEN
            t2 = t2 + fast_coef%watervapour(lay, chan, 14) * predictors%co2(1, lay, prof)
          ENDIF
          IF (coef%nco > 0) THEN
            t2 = t2 + fast_coef%watervapour(lay, chan, 15) * predictors%co(1, lay, prof)
          ENDIF
          opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
        ENDDO
      ELSE IF (chanx >= rttov9_wv2295_25 .AND. chanx < rttov9_wv2380_25) THEN
        DO lay = 1, nlayers
          t1 = fast_coef%watervapour(lay, chan, 1) * predictors%watervapour(1, lay, prof) +      &
            & fast_coef%watervapour(lay, chan, 2) * predictors%watervapour(2, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 3) * predictors%watervapour(3, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 4) * predictors%watervapour(4, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 5) * predictors%watervapour(5, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 6) * predictors%watervapour(6, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 7) * predictors%watervapour(7, lay, prof)
          t2 = fast_coef%watervapour(lay, chan, 8) * predictors%watervapour(8, lay, prof) +      &
            & fast_coef%watervapour(lay, chan, 9) * predictors%watervapour(9, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 10) * predictors%watervapour(10, lay, prof) +     &
            & fast_coef%watervapour(lay, chan, 11) * predictors%watervapour(11, lay, prof) +     &
            & fast_coef%watervapour(lay, chan, 12) * predictors%watervapour(12, lay, prof) +     &
            & fast_coef%watervapour(lay, chan, 13) * predictors%watervapour(13, lay, prof)
          opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
        ENDDO
      ELSE   !chanx >= rttov9_wv2380_25
        DO lay = 1, nlayers
          t1 = fast_coef%watervapour(lay, chan, 1) * predictors%watervapour(1, lay, prof) +      &
            & fast_coef%watervapour(lay, chan, 2) * predictors%watervapour(2, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 3) * predictors%watervapour(3, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 4) * predictors%watervapour(4, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 5) * predictors%watervapour(5, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 6) * predictors%watervapour(6, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 7) * predictors%watervapour(7, lay, prof)
          t2 = fast_coef%watervapour(lay, chan, 8) * predictors%watervapour(14, lay, prof) +      &
            & fast_coef%watervapour(lay, chan, 9) * predictors%watervapour(15, lay, prof) +       &
            & fast_coef%watervapour(lay, chan, 10) * predictors%watervapour(18, lay, prof) +      &
            & fast_coef%watervapour(lay, chan, 11) * predictors%watervapour(11, lay, prof) +      &
            & fast_coef%watervapour(lay, chan, 12) * predictors%watervapour(16, lay, prof) +      &
            & fast_coef%watervapour(lay, chan, 13) * predictors%watervapour(19, lay, prof)
          IF (coef%nch4 > 0) THEN
            t2 = t2 + fast_coef%watervapour(lay, chan, 14) * predictors%ch4(1, lay, prof) ** 1.25 +      &
              & fast_coef%watervapour(lay, chan, 15) * predictors%ch4(1, lay, prof) ** 0.25 * predictors%ch4(3, lay, prof)
          ENDIF
          opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
        ENDDO
      ENDIF
    ENDIF
  ENDDO
!-------------
!1.3 add ozone
!-------------
  IF (coef%nozone > 0) THEN
    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        chan  = chanprof(j)%chan
        chanx = coef%ff_cwn(chan)
        prof  = chanprof(j)%prof
        IF ((chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2380_25) .OR. &
            (chanx < rttov9_wv2000_00)) THEN
          DO lay = 1, nlayers
            t1 = fast_coef%ozone(lay, chan, 1) * predictors%ozone(1, lay, prof) +      &
              & fast_coef%ozone(lay, chan, 2) * predictors%ozone(2, lay, prof) +       &
              & fast_coef%ozone(lay, chan, 3) * predictors%ozone(3, lay, prof) +       &
              & fast_coef%ozone(lay, chan, 4) * predictors%ozone(4, lay, prof) +       &
              & fast_coef%ozone(lay, chan, 5) * predictors%ozone(5, lay, prof)
            t2 = fast_coef%ozone(lay, chan, 6) * predictors%ozone(6, lay, prof) +      &
              & fast_coef%ozone(lay, chan, 7) * predictors%ozone(7, lay, prof) +       &
              & fast_coef%ozone(lay, chan, 8) * predictors%ozone(8, lay, prof) +       &
              & fast_coef%ozone(lay, chan, 9) * predictors%ozone(9, lay, prof) +       &
              & fast_coef%ozone(lay, chan, 10) * predictors%ozone(10, lay, prof) +     &
              & fast_coef%ozone(lay, chan, 11) * predictors%ozone(11, lay, prof)
            opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
          ENDDO
        ELSE
          DO lay = 1, nlayers
            t1 = fast_coef%ozone(lay, chan, 1) * predictors%ozone(1, lay, prof) +      &
              & fast_coef%ozone(lay, chan, 2) * predictors%ozone(2, lay, prof) +       &
              & fast_coef%ozone(lay, chan, 3) * predictors%ozone(12, lay, prof) +      &
              & fast_coef%ozone(lay, chan, 4) * predictors%ozone(4, lay, prof) +       &
              & fast_coef%ozone(lay, chan, 5) * predictors%ozone(5, lay, prof) +       &
              & fast_coef%ozone(lay, chan, 6) * predictors%ozone(6, lay, prof)
            t2 = fast_coef%ozone(lay, chan, 7) * predictors%ozone(7, lay, prof) +      &
              & fast_coef%ozone(lay, chan, 8) * predictors%ozone(13, lay, prof) +      &
              & fast_coef%ozone(lay, chan, 9) * predictors%ozone(9, lay, prof) +       &
              & fast_coef%ozone(lay, chan, 10) * predictors%ozone(10, lay, prof) +     &
              & fast_coef%ozone(lay, chan, 11) * predictors%ozone(11, lay, prof) +     &
              & fast_coef%ozone(lay, chan, 12) * predictors%ozone(14, lay, prof) +     &
              & fast_coef%ozone(lay, chan, 13) * predictors%ozone(15, lay, prof)
            opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF
!------------------------------
!1.4 add Water Vapour Continuum
!------------------------------
  IF (coef%nwvcont > 0) THEN
    IF (coef%nwvcont == 4) THEN
      DO j = 1, nchannels
        IF (chanflag(j)) THEN
          chan = chanprof(j)%chan
          prof = chanprof(j)%prof
          DO lay = 1, nlayers
            t1 = fast_coef%wvcont(lay, chan, 1) * predictors%wvcont(1, lay, prof) +      &
              & fast_coef%wvcont(lay, chan, 2) * predictors%wvcont(2, lay, prof)
            t2 = fast_coef%wvcont(lay, chan, 3) * predictors%wvcont(3, lay, prof) +      &
              & fast_coef%wvcont(lay, chan, 4) * predictors%wvcont(4, lay, prof)
            opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
          ENDDO
        ENDIF
      ENDDO
    ELSE
      DO j = 1, nchannels
        IF (chanflag(j)) THEN
          chan = chanprof(j)%chan
          prof = chanprof(j)%prof
          DO lay = 1, nlayers
            DO ii = 1, coef%nwvcont
              opticaldepth(lay, j) = opticaldepth(lay, j) + fast_coef%wvcont(lay, chan, ii) * predictors%wvcont(ii, lay, prof)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDIF
  ENDIF
!-----------
!1.5 add CO2
!-----------
  IF (coef%nco2 > 0) THEN
    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        chan  = chanprof(j)%chan
        chanx = coef%ff_cwn(chan)
        prof  = chanprof(j)%prof
        
        IF (chanx >= rttov9_wv0690_50 .AND. chanx < rttov9_wv1100_25) THEN
          DO lay = 1, nlayers
            t1 = fast_coef%co2(lay, chan, 1) * predictors%co2(1, lay, prof) + &
              & fast_coef%co2(lay, chan, 2) * predictors%co2(2, lay, prof) + &
              & fast_coef%co2(lay, chan, 3) * predictors%co2(5, lay, prof) * predictors%co2(6, lay, prof) ** 2 + &
              & fast_coef%co2(lay, chan, 4) * predictors%co2(3, lay, prof) ** 2
            t2 = fast_coef%co2(lay, chan, 5) * predictors%co2(5, lay, prof) +       &
              & fast_coef%co2(lay, chan, 6) * predictors%co2(6, lay, prof) ** 2 +  &
              & fast_coef%co2(lay, chan, 7) * predictors%co2(7, lay, prof) ** 2 +  &
              & fast_coef%co2(lay, chan, 8) * predictors%co2(8, lay, prof) +       &
              & fast_coef%co2(lay, chan, 9) * predictors%co2(9, lay, prof)
            opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
          ENDDO
        ELSE IF (chanx >= rttov9_wv1995_00 .AND. chanx < rttov9_wv2000_00) THEN
          DO lay = 1, nlayers
            t1 = fast_coef%co2(lay, chan, 1) * predictors%co2(1, lay, prof) + &
              & fast_coef%co2(lay, chan, 2) * predictors%co2(2, lay, prof) + &
              & fast_coef%co2(lay, chan, 3) * predictors%co2(3, lay, prof) + &
              & fast_coef%co2(lay, chan, 4) * predictors%co2(4, lay, prof) + &
              & fast_coef%co2(lay, chan, 5) * predictors%co2(5, lay, prof)
            t2 = fast_coef%co2(lay, chan, 6) * predictors%co2(6, lay, prof) + &
              & fast_coef%co2(lay, chan, 7) * predictors%co2(7, lay, prof) + &
              & fast_coef%co2(lay, chan, 8) * predictors%co2(8, lay, prof) + &
              & fast_coef%co2(lay, chan, 9) * predictors%co2(9, lay, prof)
            IF (coef%nco > 0) THEN
              t2 = t2 + fast_coef%co2(lay, chan, 10) * predictors%co(1, lay, prof)
            ENDIF
            t2 = t2 + fast_coef%co2(lay, chan, 11) * predictors%co2(10, lay, prof)
            opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
          ENDDO
        ELSE IF (chanx < rttov9_wv2000_00) THEN
          DO lay = 1, nlayers
            t1 = fast_coef%co2(lay, chan, 1) * predictors%co2(1, lay, prof) + &
              & fast_coef%co2(lay, chan, 2) * predictors%co2(2, lay, prof) + &
              & fast_coef%co2(lay, chan, 3) * predictors%co2(3, lay, prof) + &
              & fast_coef%co2(lay, chan, 4) * predictors%co2(4, lay, prof) + &
              & fast_coef%co2(lay, chan, 5) * predictors%co2(5, lay, prof)
            t2 = fast_coef%co2(lay, chan, 6) * predictors%co2(6, lay, prof) + &
              & fast_coef%co2(lay, chan, 7) * predictors%co2(7, lay, prof) + &
              & fast_coef%co2(lay, chan, 8) * predictors%co2(8, lay, prof) + &
              & fast_coef%co2(lay, chan, 9) * predictors%co2(9, lay, prof)
            opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
          ENDDO
        ELSE IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2250_00) THEN
          DO lay = 1, nlayers
            t1 = fast_coef%co2(lay, chan, 1) * predictors%co2(1, lay, prof) + &
                & fast_coef%co2(lay, chan, 2) * predictors%co2(2, lay, prof) + &
                & fast_coef%co2(lay, chan, 3) * predictors%co2(3, lay, prof) + &
                & fast_coef%co2(lay, chan, 4) * predictors%co2(4, lay, prof) + &
                & fast_coef%co2(lay, chan, 5) * predictors%co2(5, lay, prof) + &
                & fast_coef%co2(lay, chan, 6) * predictors%co2(7, lay, prof) + &
                & fast_coef%co2(lay, chan, 7) * predictors%co2(8, lay, prof)
            t2 = 0._jprb
            IF (coef%nco > 0) THEN
              t2 = t2 + fast_coef%co2(lay, chan, 8) * predictors%co(1, lay, prof)
            ENDIF
            t2 = t2 + fast_coef%co2(lay, chan, 9) * predictors%co2(10, lay, prof) +      &
              & fast_coef%co2(lay, chan, 10) * predictors%co2(11, lay, prof) +           &
              & fast_coef%co2(lay, chan, 11) * predictors%co2(12, lay, prof) +           &
              & fast_coef%co2(lay, chan, 12) * predictors%co2(13, lay, prof) +           &
              & fast_coef%co2(lay, chan, 13) * predictors%co2(14, lay, prof) +           &
              & fast_coef%co2(lay, chan, 14) * predictors%co2(15, lay, prof)
            opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
          ENDDO
        ELSE IF (chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2295_25) THEN
          DO lay = 1, nlayers
            t1 = fast_coef%co2(lay, chan, 1) * predictors%co2(1, lay, prof) + &
              & fast_coef%co2(lay, chan, 2) * predictors%co2(2, lay, prof) + &
              & fast_coef%co2(lay, chan, 3) * predictors%co2(3, lay, prof) + &
              & fast_coef%co2(lay, chan, 4) * predictors%co2(4, lay, prof) + &
              & fast_coef%co2(lay, chan, 5) * predictors%co2(5, lay, prof)
            t2 = fast_coef%co2(lay, chan, 6) * predictors%co2(6, lay, prof) + &
              & fast_coef%co2(lay, chan, 7) * predictors%co2(7, lay, prof) + &
              & fast_coef%co2(lay, chan, 8) * predictors%co2(8, lay, prof) + &
              & fast_coef%co2(lay, chan, 9) * predictors%co2(9, lay, prof)
            IF (coef%nco > 0) THEN
              t2 = t2 + fast_coef%co2(lay, chan, 10) * predictors%co(1, lay, prof)
            ENDIF
            t2 = t2 + fast_coef%co2(lay, chan, 11) * predictors%co2(10, lay, prof)
            opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
          ENDDO
        ELSE IF (chanx >= rttov9_wv2295_25 .AND. chanx < rttov9_wv2380_25) THEN
          DO lay = 1, nlayers
            t1 = fast_coef%co2(lay, chan, 1) * predictors%co2(1, lay, prof) + &
              & fast_coef%co2(lay, chan, 2) * predictors%co2(2, lay, prof) + &
              & fast_coef%co2(lay, chan, 3) * predictors%co2(3, lay, prof) + &
              & fast_coef%co2(lay, chan, 4) * predictors%co2(4, lay, prof) + &
              & fast_coef%co2(lay, chan, 5) * predictors%co2(5, lay, prof)
            t2 = fast_coef%co2(lay, chan, 6) * predictors%co2(6, lay, prof) + &
              & fast_coef%co2(lay, chan, 7) * predictors%co2(7, lay, prof) + &
              & fast_coef%co2(lay, chan, 8) * predictors%co2(8, lay, prof) + &
              & fast_coef%co2(lay, chan, 9) * predictors%co2(9, lay, prof) + &
              & fast_coef%co2(lay, chan, 10) * predictors%co2(10, lay, prof)
            opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
          ENDDO        
        ELSE   ! chanx >= rttov9_wv2380_25
          DO lay = 1, nlayers
            t1 = fast_coef%co2(lay, chan, 1) * predictors%co2(1, lay, prof) + &
                & fast_coef%co2(lay, chan, 2) * predictors%co2(2, lay, prof) + &
                & fast_coef%co2(lay, chan, 3) * predictors%co2(3, lay, prof) + &
                & fast_coef%co2(lay, chan, 4) * predictors%co2(4, lay, prof) + &
                & fast_coef%co2(lay, chan, 5) * predictors%co2(5, lay, prof) + &
                & fast_coef%co2(lay, chan, 6) * predictors%co2(7, lay, prof)
            t2 = fast_coef%co2(lay, chan, 7) * predictors%co2(8, lay, prof) + &
                & fast_coef%co2(lay, chan, 8) * predictors%co2(10, lay, prof) + &
                & fast_coef%co2(lay, chan, 9) * predictors%co2(11, lay, prof) + &
                & fast_coef%co2(lay, chan, 10) * predictors%co2(12, lay, prof) + &
                & fast_coef%co2(lay, chan, 11) * predictors%co2(13, lay, prof) + &
                & fast_coef%co2(lay, chan, 12) * predictors%co2(14, lay, prof) + &
                & fast_coef%co2(lay, chan, 13) * predictors%co2(15, lay, prof)
            opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF
!-----------
!1.6 add N2O
!-----------
  IF (coef%nn2o > 0) THEN
    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        chan  = chanprof(j)%chan
        chanx = coef%ff_cwn(chan)
        prof  = chanprof(j)%prof
        IF (chanx >= rttov9_wv1050_00 .AND. chanx < rttov9_wv1350_25) THEN
          DO lay = 1, nlayers
            t1 = fast_coef%n2o(lay, chan, 1) * predictors%n2o(1, lay, prof) + &
              & fast_coef%n2o(lay, chan, 2) * predictors%n2o(2, lay, prof) + &
              & fast_coef%n2o(lay, chan, 3) * predictors%n2o(3, lay, prof) + &
              & fast_coef%n2o(lay, chan, 4) * predictors%n2o(4, lay, prof) + &
              & fast_coef%n2o(lay, chan, 5) * predictors%n2o(5, lay, prof) + &
              & fast_coef%n2o(lay, chan, 6) * predictors%n2o(6, lay, prof)
            t2 = fast_coef%n2o(lay, chan, 7) * predictors%n2o(7, lay, prof) + &
              & fast_coef%n2o(lay, chan, 8) * predictors%n2o(8, lay, prof) + &
              & fast_coef%n2o(lay, chan, 9) * predictors%n2o(9, lay, prof) + &
              & fast_coef%n2o(lay, chan, 10) * predictors%n2o(10, lay, prof)
            IF (coef%nch4 > 0) THEN
              t2 = t2 + fast_coef%n2o(lay, chan, 11) * predictors%ch4(1, lay, prof) +      &
                & fast_coef%n2o(lay, chan, 12) * predictors%ch4(10, lay, prof) 
            ENDIF
            opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
          ENDDO
        ELSE IF (chanx >= rttov9_wv1995_00 .AND. chanx < rttov9_wv2000_00) THEN
          DO lay = 1, nlayers
            t1 = fast_coef%n2o(lay, chan, 1) * predictors%n2o(1, lay, prof) + &
              & fast_coef%n2o(lay, chan, 2) * predictors%n2o(2, lay, prof) + &
              & fast_coef%n2o(lay, chan, 3) * predictors%n2o(3, lay, prof) + &
              & fast_coef%n2o(lay, chan, 4) * predictors%n2o(4, lay, prof) + &
              & fast_coef%n2o(lay, chan, 5) * predictors%n2o(5, lay, prof) + &
              & fast_coef%n2o(lay, chan, 6) * predictors%n2o(6, lay, prof)
            t2 = fast_coef%n2o(lay, chan, 7) * predictors%n2o(7, lay, prof) + &
              & fast_coef%n2o(lay, chan, 8) * predictors%n2o(8, lay, prof) + &
              & fast_coef%n2o(lay, chan, 9) * predictors%n2o(9, lay, prof) + &
              & fast_coef%n2o(lay, chan, 10) * predictors%n2o(10, lay, prof)
            IF (coef%nco > 0) THEN
              t2 = t2 + fast_coef%n2o(lay, chan, 11) * predictors%co(1, lay, prof) +      &
                & fast_coef%n2o(lay, chan, 12) * predictors%co(1, lay, prof) ** 3 / predictors%co(8, lay, prof)
            ENDIF
            opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
          ENDDO
        ELSE IF (chanx < rttov9_wv2000_00) THEN
          DO lay = 1, nlayers
            t1 = fast_coef%n2o(lay, chan, 1) * predictors%n2o(1, lay, prof) + &
              & fast_coef%n2o(lay, chan, 2) * predictors%n2o(2, lay, prof) + &
              & fast_coef%n2o(lay, chan, 3) * predictors%n2o(3, lay, prof) + &
              & fast_coef%n2o(lay, chan, 4) * predictors%n2o(4, lay, prof) + &
              & fast_coef%n2o(lay, chan, 5) * predictors%n2o(5, lay, prof)
            t2 = fast_coef%n2o(lay, chan, 6) * predictors%n2o(6, lay, prof) + &
              & fast_coef%n2o(lay, chan, 7) * predictors%n2o(7, lay, prof) + &
              & fast_coef%n2o(lay, chan, 8) * predictors%n2o(8, lay, prof) + &
              & fast_coef%n2o(lay, chan, 9) * predictors%n2o(9, lay, prof) + &
              & fast_coef%n2o(lay, chan, 10) * predictors%n2o(10, lay, prof)
            opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
          ENDDO
        ELSE IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2250_00) THEN
          DO lay = 1, nlayers
            t1 = fast_coef%n2o(lay, chan, 1) * predictors%n2o(1, lay, prof) + &
              & fast_coef%n2o(lay, chan, 2) * predictors%n2o(2, lay, prof) + &
              & fast_coef%n2o(lay, chan, 3) * predictors%n2o(3, lay, prof) + &
              & fast_coef%n2o(lay, chan, 4) * predictors%n2o(4, lay, prof) + &
              & fast_coef%n2o(lay, chan, 5) * predictors%n2o(5, lay, prof) + &
              & fast_coef%n2o(lay, chan, 6) * predictors%n2o(6, lay, prof) + &
              & fast_coef%n2o(lay, chan, 7) * predictors%n2o(7, lay, prof)
            t2 = fast_coef%n2o(lay, chan, 8) * predictors%n2o(8, lay, prof) + &
              & fast_coef%n2o(lay, chan, 9) * predictors%n2o(11, lay, prof) + &
              & fast_coef%n2o(lay, chan, 10) * predictors%n2o(10, lay, prof)
            IF (coef%nco > 0) THEN
              t2 = t2 + fast_coef%n2o(lay, chan, 11) * predictors%co(1, lay, prof) +      &
                & fast_coef%n2o(lay, chan, 12) * predictors%co(1, lay, prof) ** 3 / predictors%co(8, lay, prof)
            ENDIF
            t2 = t2 + fast_coef%n2o(lay, chan, 13) * predictors%n2o(12, lay, prof) +      &
              & fast_coef%n2o(lay, chan, 14) * predictors%n2o(13, lay, prof)
            opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
          ENDDO
        ELSE IF (chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2295_25) THEN
          DO lay = 1, nlayers
            t1 = fast_coef%n2o(lay, chan, 1) * predictors%n2o(1, lay, prof) + &
              & fast_coef%n2o(lay, chan, 2) * predictors%n2o(2, lay, prof) + &
              & fast_coef%n2o(lay, chan, 3) * predictors%n2o(3, lay, prof) + &
              & fast_coef%n2o(lay, chan, 4) * predictors%n2o(4, lay, prof) + &
              & fast_coef%n2o(lay, chan, 5) * predictors%n2o(5, lay, prof) + &
              & fast_coef%n2o(lay, chan, 6) * predictors%n2o(6, lay, prof)
            t2 = fast_coef%n2o(lay, chan, 7) * predictors%n2o(7, lay, prof) + &
              & fast_coef%n2o(lay, chan, 8) * predictors%n2o(8, lay, prof) + &
              & fast_coef%n2o(lay, chan, 9) * predictors%n2o(9, lay, prof) + &
              & fast_coef%n2o(lay, chan, 10) * predictors%n2o(10, lay, prof)
            IF (coef%nco > 0) THEN
              t2 = t2 + fast_coef%n2o(lay, chan, 11) * predictors%co(1, lay, prof) +      &
                & fast_coef%n2o(lay, chan, 12) * predictors%co(1, lay, prof) ** 3 / predictors%co(8, lay, prof)
            ENDIF
            opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
          ENDDO
        ELSE IF (chanx >= rttov9_wv2295_25 .AND. chanx < rttov9_wv2380_25) THEN
          DO lay = 1, nlayers
            t1 = fast_coef%n2o(lay, chan, 1) * predictors%n2o(1, lay, prof) + &
              & fast_coef%n2o(lay, chan, 2) * predictors%n2o(2, lay, prof) + &
              & fast_coef%n2o(lay, chan, 3) * predictors%n2o(3, lay, prof) + &
              & fast_coef%n2o(lay, chan, 4) * predictors%n2o(4, lay, prof) + &
              & fast_coef%n2o(lay, chan, 5) * predictors%n2o(5, lay, prof) + &
              & fast_coef%n2o(lay, chan, 6) * predictors%n2o(6, lay, prof)
            t2 = fast_coef%n2o(lay, chan, 7) * predictors%n2o(7, lay, prof) + &
              & fast_coef%n2o(lay, chan, 8) * predictors%n2o(8, lay, prof) + &
              & fast_coef%n2o(lay, chan, 9) * predictors%n2o(9, lay, prof) + &
              & fast_coef%n2o(lay, chan, 10) * predictors%n2o(10, lay, prof)
            opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
          ENDDO
        ELSE   ! chanx >= rttov9_wv2380_25
          DO lay = 1, nlayers
            t1 = fast_coef%n2o(lay, chan, 1) * predictors%n2o(1, lay, prof) + &
              & fast_coef%n2o(lay, chan, 2) * predictors%n2o(2, lay, prof) + &
              & fast_coef%n2o(lay, chan, 3) * predictors%n2o(3, lay, prof) + &
              & fast_coef%n2o(lay, chan, 4) * predictors%n2o(4, lay, prof) + &
              & fast_coef%n2o(lay, chan, 5) * predictors%n2o(5, lay, prof) + &
              & fast_coef%n2o(lay, chan, 6) * predictors%n2o(6, lay, prof) + &
              & fast_coef%n2o(lay, chan, 7) * predictors%n2o(7, lay, prof)
            t2 = fast_coef%n2o(lay, chan, 8) * predictors%n2o(8, lay, prof) + &
              & fast_coef%n2o(lay, chan, 9) * predictors%n2o(11, lay, prof) + &
              & fast_coef%n2o(lay, chan, 10) * predictors%n2o(10, lay, prof) + &
              & fast_coef%n2o(lay, chan, 11) * predictors%n2o(12, lay, prof) + &
              & fast_coef%n2o(lay, chan, 12) * predictors%n2o(13, lay, prof)
            opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF
!-----------
!1.7 add CO
!-----------
  IF (coef%nco > 0) THEN
    DO j = 1, nchannels
      IF (chanflag(j)) THEN
        chan  = chanprof(j)%chan
        chanx = coef%ff_cwn(chan)
        prof  = chanprof(j)%prof
        IF ((chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2380_25) .OR. &
            (chanx < rttov9_wv2000_00)) THEN
          DO lay = 1, nlayers
            t1 = fast_coef%co(lay, chan, 1) * predictors%co(1, lay, prof) + &
              & fast_coef%co(lay, chan, 2) * predictors%co(2, lay, prof) + &
              & fast_coef%co(lay, chan, 3) * predictors%co(3, lay, prof) + &
              & fast_coef%co(lay, chan, 4) * predictors%co(4, lay, prof) + &
              & fast_coef%co(lay, chan, 5) * predictors%co(5, lay, prof) + &
              & fast_coef%co(lay, chan, 6) * predictors%co(6, lay, prof)
            t2 = fast_coef%co(lay, chan, 7) * predictors%co(7, lay, prof) + &
              & fast_coef%co(lay, chan, 8) * predictors%co(8, lay, prof) + &
              & fast_coef%co(lay, chan, 9) * predictors%co(9, lay, prof) + &
              & fast_coef%co(lay, chan, 10) * predictors%co(10, lay, prof) + &
              & fast_coef%co(lay, chan, 11) * predictors%co(11, lay, prof)
            opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
          ENDDO
        ELSE
          DO lay = 1, nlayers
            t1 = fast_coef%co(lay, chan, 1) * predictors%co(1, lay, prof) + &
              & fast_coef%co(lay, chan, 2) * predictors%co(2, lay, prof) + &
              & fast_coef%co(lay, chan, 3) * predictors%co(3, lay, prof) + &
              & fast_coef%co(lay, chan, 4) * predictors%co(4, lay, prof) + &
              & fast_coef%co(lay, chan, 5) * predictors%co(5, lay, prof) + &
              & fast_coef%co(lay, chan, 6) * predictors%co(6, lay, prof)
            t2 = fast_coef%co(lay, chan, 7) * predictors%co(7, lay, prof) + &
              & fast_coef%co(lay, chan, 8) * predictors%co(8, lay, prof) + &
              & fast_coef%co(lay, chan, 9) * predictors%co(9, lay, prof) + &
              & fast_coef%co(lay, chan, 10) * predictors%co(10, lay, prof) + &
              & fast_coef%co(lay, chan, 11) * predictors%co(12, lay, prof) + &
              & fast_coef%co(lay, chan, 12) * predictors%co(13, lay, prof)
            opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF
!-----------
!1.8 add CH4
!-----------
  IF (coef%nch4 > 0) THEN
    IF (coef%nch4 == 11) THEN
      DO j = 1, nchannels
        IF (chanflag(j)) THEN
          chan = chanprof(j)%chan
          prof = chanprof(j)%prof
          DO lay = 1, nlayers
            t1 = fast_coef%ch4(lay, chan, 1) * predictors%ch4(1, lay, prof) + &
              & fast_coef%ch4(lay, chan, 2) * predictors%ch4(2, lay, prof) + &
              & fast_coef%ch4(lay, chan, 3) * predictors%ch4(3, lay, prof) + &
              & fast_coef%ch4(lay, chan, 4) * predictors%ch4(4, lay, prof) + &
              & fast_coef%ch4(lay, chan, 5) * predictors%ch4(5, lay, prof)
            t2 = fast_coef%ch4(lay, chan, 6) * predictors%ch4(6, lay, prof) + &
              & fast_coef%ch4(lay, chan, 7) * predictors%ch4(7, lay, prof) + &
              & fast_coef%ch4(lay, chan, 8) * predictors%ch4(8, lay, prof) + &
              & fast_coef%ch4(lay, chan, 9) * predictors%ch4(9, lay, prof) + &
              & fast_coef%ch4(lay, chan, 10) * predictors%ch4(10, lay, prof) + &
              & fast_coef%ch4(lay, chan, 11) * predictors%ch4(11, lay, prof)
            opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
          ENDDO
        ENDIF
      ENDDO
    ELSE
      DO j = 1, nchannels
        IF (chanflag(j)) THEN
          chan = chanprof(j)%chan
          prof = chanprof(j)%prof
          DO lay = 1, nlayers
            DO ii = 1, coef%nch4
              opticaldepth(lay, j) = opticaldepth(lay, j) + fast_coef%ch4(lay, chan, ii) * predictors%ch4(ii, lay, prof)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDIF
  ENDIF
!--------------------
!1.9 add liquid water (MW only)
!--------------------
  IF (coef%id_sensor == sensor_id_mw) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      prof = chanprof(j)%prof
      debye_prof => aux%debye_prof(:, :, prof)
      IF (pred_info%ncloud >= 1) THEN
        DO lay = 1, nlayers
          lev = lay + 1
          IF (lev >= coef%mwcldtop) THEN
            zf = coef%frequency_ghz(chan)
            zf_sq                = zf * zf
            z1_sq                = debye_prof(1, lev) * debye_prof(1, lev)
            z2_sq                = debye_prof(2, lev) * debye_prof(2, lev)
            z34_dif              = debye_prof(3, lev) - debye_prof(4, lev)
            z45_dif              = debye_prof(4, lev) - debye_prof(5, lev)
            z1f_sq_z1_sq         = zf_sq + z1_sq
            z2f_sq_z2_sq         = zf_sq + z2_sq
            z1_div               = 1.0_jprb / z1f_sq_z1_sq
            z2_div               = 1.0_jprb / z2f_sq_z2_sq
            z1_den               = z34_dif * z1_div
            z2_den               = z45_dif * z2_div
            zastar               = debye_prof(3, lev) - zf_sq * (z1_den + z2_den)
            z1_prod              = z34_dif * debye_prof(1, lev)
            z2_prod              = z1_prod * z1_div
            z3_prod              = z45_dif * debye_prof(2, lev)
            z4_prod              = z3_prod * z2_div
            zbstar               =  - zf * (z2_prod + z4_prod)
            zbstar_sq            = zbstar * zbstar
            za2star              = zastar + 2.0_jprb
            za2star_sq           = za2star * za2star
            zdiv = za2star_sq + zbstar_sq
            zgstar               =  - 3.0_jprb * zbstar / zdiv
            opticaldepth(lay, j) = opticaldepth(lay, j) - 1.5_jprb * zf * zgstar * predictors%clw(lay, prof)
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDIF
!----------------------------------------
!2. Assemble layer optical depths
!----------------------------------------
! note that optical depth in the calculations is negative
! single layer
! store (neg) optical depth in reference array on COEF levs ...
! ... for TL, AD and K calculations
  DO j = 1, nchannels
    IF (chanflag(j)) THEN
      opdp_ref(:,j) = opticaldepth(:,j)
      WHERE (opticaldepth(:,j) > 0.0_jprb)
        opticaldepth(:,j) = 0.0_jprb
      ENDWHERE
    ENDIF
  ENDDO
! level to space optical depths
! these are path (i.e. cumulated layer) optical depths
  DO j = 1, nchannels
    IF (chanflag(j)) THEN
      opdp_path(1, j) = 0.0_jprb
      DO lev = 2, nlevels
        lay = lev - 1
        opdp_path(lev, j) = opdp_path(lev - 1, j) + opticaldepth(lay, j)
      ENDDO
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_9', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_opdep_9
