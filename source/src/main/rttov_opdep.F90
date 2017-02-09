SUBROUTINE rttov_opdep( &
            & nlayers,    &
            & chanprof,   &
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
! To interpolate optical depths on to levels of radiative transfer model
! (which, at present, entails only surface transmittance, as
! other optical depths are on *rt* levels) and to convert
! optical depths to transmittances.
!
! Code based on OPDEP and RTTAU from previous versions of RTTOV
! Only one profile per call
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
!  1.0    01/12/2002  New F90 code with structures (P Brunel A Smith)
!  1.1    29/01/2003  Add WV Continuum and CO2 capability (P Brunel)
!  1.2    04/12/2003  Optimisation (J Hague and D Salmond ECMWF)
!  1.3    26/09/2003  Modified to allow for multiple polarisations (S English)
!  1.4    17/08/2004  Bug fixed in setting transmission to 1 (S English)
!  1.5    28/02/2005  Improved vectorisation (D Dent)
!  1.6    01/06/2005  Marco Matricardi (ECMWF):
!            --       Computation of variable OD_SFRAC added.
!  1.7    24/01/2007  Removed polarisation index (R Saunders)
!  1.8    27/02/2009  Profile levels to include ToA. Distinguish between
!                     layer arrays and level arrays - size, index
!                     labels, looping (P. Rayer)
!  1.9    03/11/2009  Transmittances / optical depths on levels (A Geer)
!  1.10   10/01/2013  Add PMC shifts (P Rayer)
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
       & rttov_fast_coef, &
       & predictors_Type, &
       & rttov_path_pred, &
       & opdp_path_Type,  &
       & profile_aux
  USE parkind1, ONLY : jpim, jprb
!INTF_OFF
  USE rttov_const, ONLY : sensor_id_mw
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  INTEGER(KIND=jpim)   , INTENT(IN)    :: nlayers                          ! Number of pressure levels
  TYPE(rttov_chanprof ), INTENT(IN)    :: chanprof(:)                      ! Channel indices
  TYPE(predictors_Type), INTENT(IN)    :: pred_info                        ! Predictor structure
  TYPE(rttov_path_pred), INTENT(IN)    :: predictors                       ! Predictors
  TYPE(rttov_coef     ), INTENT(IN)    :: coef                             ! Coefficient structure
  TYPE(rttov_fast_coef), INTENT(IN)    :: fast_coef                        ! Fast coefficients
  TYPE(profile_aux    ), INTENT(IN)    :: aux                              ! auxillary profiles informations
  TYPE(opdp_path_Type ), INTENT(INOUT) :: opdp_path                        ! optical depths
  REAL(KIND=jprb)      , INTENT(OUT)   :: opdp_ref(nlayers, size(chanprof))! layer optical depth
!INTF_END
!local variables:
  REAL(KIND=jprb) :: opticaldepth(nlayers, size(chanprof))! raw layer optical depth
  REAL(KIND=jprb), POINTER :: debye_prof(:, :)! pointer on Debye profiles
  INTEGER(KIND=jpim) :: lev         , lay         , chan   , j, prof   , nlevels                  ! loop variables
! cloud liquid water local variables
  REAL   (KIND=jprb) :: zf, zf_sq       , z34_dif, z45_dif   , z1_sq  , z2_sq  , z1_div , z2_div
  REAL   (KIND=jprb) :: z1_den      , z2_den      , zastar , z1_prod   , z2_prod, z3_prod, z4_prod
  REAL   (KIND=jprb) :: zbstar      , zbstar_sq   , za2star, za2star_sq, zdiv   , zgstar
  REAL   (KIND=jprb) :: z1f_sq_z1_sq, z2f_sq_z2_sq, t1     , t2, ztemp
  INTEGER(KIND=jpim) :: ii
  INTEGER(KIND=jpim) :: nchannels
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP', 0_jpim, ZHOOK_HANDLE)
  nchannels = size(chanprof)
  nlevels   = nlayers + 1
!-----------------------------------------
!1. calculate layer gaseous optical depths
!-----------------------------------------
!--------------------------
!1.1 start with mixed gases
!--------------------------
  DO j = 1, nchannels
    chan = chanprof(j)%chan
    prof = chanprof(j)%prof
! unrolling reduces memory traffic and increases performance.
    IF (coef%nmixed == 10) THEN
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
    ELSE
      opticaldepth(1:nlayers, j) = 0._JPRB
      DO ii = 1, coef%nmixed
        DO lay = 1, nlayers
          opticaldepth(lay, j) =      &
            & opticaldepth(lay, j) + fast_coef%mixedgas(lay, chan, ii) * predictors%mixedgas(ii, lay, prof)
        ENDDO
      ENDDO
    ENDIF
  ENDDO
!--------------------
!1.2 add water vapour
!--------------------
  DO j = 1, nchannels
    chan = chanprof(j)%chan
    prof = chanprof(j)%prof
    IF (coef%nwater == 15) THEN
      DO lay = 1, nlayers
        t1 = fast_coef%watervapour(lay, chan, 1) * predictors%watervapour(1, lay, prof) +      &
          & fast_coef%watervapour(lay, chan, 2) * predictors%watervapour(2, lay, prof) +       &
          & fast_coef%watervapour(lay, chan, 3) * predictors%watervapour(3, lay, prof) +       &
          & fast_coef%watervapour(lay, chan, 4) * predictors%watervapour(4, lay, prof)
        t2 = fast_coef%watervapour(lay, chan, 5) * predictors%watervapour(5, lay, prof) +      &
          & fast_coef%watervapour(lay, chan, 6) * predictors%watervapour(6, lay, prof) +       &
          & fast_coef%watervapour(lay, chan, 7) * predictors%watervapour(7, lay, prof) +       &
          & fast_coef%watervapour(lay, chan, 8) * predictors%watervapour(8, lay, prof)
        opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
      ENDDO
! Split to stop register SPILLs (DJS)
      DO lay = 1, nlayers
        t1 = fast_coef%watervapour(lay, chan, 9) * predictors%watervapour(9, lay, prof) +      &
          & fast_coef%watervapour(lay, chan, 10) * predictors%watervapour(10, lay, prof) +     &
          & fast_coef%watervapour(lay, chan, 11) * predictors%watervapour(11, lay, prof) +     &
          & fast_coef%watervapour(lay, chan, 12) * predictors%watervapour(12, lay, prof)
        t2 = fast_coef%watervapour(lay, chan, 13) * predictors%watervapour(13, lay, prof) +      &
          & fast_coef%watervapour(lay, chan, 14) * predictors%watervapour(14, lay, prof) +       &
          & fast_coef%watervapour(lay, chan, 15) * predictors%watervapour(15, lay, prof)
        opticaldepth(lay, j) = opticaldepth(lay, j) + t1 + t2
      ENDDO
    ELSE
      DO ii = 1, coef%nwater
        DO lay = 1, nlayers
          opticaldepth(lay, j) =      &
            & opticaldepth(lay, j) + fast_coef%watervapour(lay, chan, ii) * predictors%watervapour(ii, lay, prof)
        ENDDO
      ENDDO
    ENDIF
  ENDDO
!-------------
!1.3 add ozone
!-------------
  IF (coef%nozone > 0) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      prof = chanprof(j)%prof
      IF (coef%nozone == 11) THEN
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
        DO ii = 1, coef%nozone
          DO lay = 1, nlayers
            opticaldepth(lay, j) = opticaldepth(lay, j) + fast_coef%ozone(lay, chan, ii) * predictors%ozone(ii, lay, prof)
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  ENDIF
!------------------------------
!1.4 add Water Vapour Continuum
!------------------------------
  IF (coef%nwvcont > 0) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      prof = chanprof(j)%prof
      DO ii = 1, coef%nwvcont
        DO lay = 1, nlayers
          opticaldepth(lay, j) = opticaldepth(lay, j) + fast_coef%wvcont(lay, chan, ii) * predictors%wvcont(ii, lay, prof)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
!-----------
!1.5 add CO2
!-----------
  IF (coef%nco2 > 0) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      prof = chanprof(j)%prof
      DO ii = 1, coef%nco2
        DO lay = 1, nlayers
          opticaldepth(lay, j) = opticaldepth(lay, j) + fast_coef%co2(lay, chan, ii) * predictors%co2(ii, lay, prof)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
! N2O, CO and CH4 are not variable gases for v7 and v8 predictors.
! !-----------
! !1.6 add N2O
! !-----------
!   DO j = 1, nchannels
!     chan = chanprof(j)%chan
!     prof = chanprof(j)%prof
!     DO ii = 1, coef%nn2o
!       DO lay = 1, nlayers
!         opticaldepth(lay, j) = opticaldepth(lay, j) + fast_coef%n2o(lay, chan, ii) * predictors%n2o(ii, lay, prof)
!       ENDDO
!     ENDDO
!   ENDDO
! !-----------
! !1.7 add CO
! !-----------
!   DO j = 1, nchannels
!     chan = chanprof(j)%chan
!     prof = chanprof(j)%prof
!     DO ii = 1, coef%nco
!       DO lay = 1, nlayers
!         opticaldepth(lay, j) = opticaldepth(lay, j) + fast_coef%co(lay, chan, ii) * predictors%co(ii, lay, prof)
!       ENDDO
!     ENDDO
!   ENDDO
! !-----------
! !1.8 add CH4
! !-----------
!   DO j = 1, nchannels
!     chan = chanprof(j)%chan
!     prof = chanprof(j)%prof
!     DO ii = 1, coef%nch4
!       DO lay = 1, nlayers
!         opticaldepth(lay, j) = opticaldepth(lay, j) + fast_coef%ch4(lay, chan, ii) * predictors%ch4(ii, lay, prof)
!       ENDDO
!     ENDDO
!   ENDDO
!--------------------
!1.9 add liquid water (MW only)
!--------------------
  IF (coef%id_sensor == sensor_id_mw) THEN
    IF (pred_info%ncloud > 0) THEN
      DO j = 1, nchannels
        chan = chanprof(j)%chan
        prof = chanprof(j)%prof
        debye_prof => aux%debye_prof(:, :, prof)
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
            z1_div               = 1.0_JPRB / z1f_sq_z1_sq
            z2_div               = 1.0_JPRB / z2f_sq_z2_sq
            z1_den               = z34_dif * z1_div
            z2_den               = z45_dif * z2_div
            zastar               = debye_prof(3, lev) - zf_sq * (z1_den + z2_den)
            z1_prod              = z34_dif * debye_prof(1, lev)
            z2_prod              = z1_prod * z1_div
            z3_prod              = z45_dif * debye_prof(2, lev)
            z4_prod              = z3_prod * z2_div
            zbstar               =  - zf * (z2_prod + z4_prod)
            zbstar_sq            = zbstar * zbstar
            za2star              = zastar + 2.0_JPRB
            za2star_sq           = za2star * za2star
            zdiv = za2star_sq + zbstar_sq
            zgstar               =  - 3.0_JPRB * zbstar / zdiv
            opticaldepth(lay, j) = opticaldepth(lay, j) - 1.5_JPRB * zf * zgstar * predictors%clw(lay, prof)
          ENDIF
        ENDDO
      ENDDO
    ENDIF
  ENDIF
!------------------------------------------
!1.10 pressure-modulated cell (pmc) sensors
!------------------------------------------
  IF (coef%pmc_shift) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      prof = chanprof(j)%prof
      DO ii = 1, coef%pmc_nvar
        DO lay = 1, coef%pmc_nlay
          opticaldepth(lay, j) = opticaldepth(lay, j) + coef%pmc_coef(lay, chan, ii) * predictors%pmc(ii, lay, prof, chan)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

!----------------------------------------
!2. Assemble layer optical depths
!----------------------------------------
! note that optical depth in the calculations is negative
! single layer
! store (neg) optical depth in reference array on COEF lays ...
! ... for TL, AD and K calculations
!  Replace array syntax with Do loops (DJS)
  DO j = 1, nchannels
    DO lay = 1, nlayers
      opdp_ref(lay, j) = opticaldepth(lay, j)
      IF (opticaldepth(lay, j) > 0.0_JPRB) THEN
        opticaldepth(lay, j) = 0.0_JPRB
      ENDIF
    ENDDO
  ENDDO
! level to space optical depths
  DO j = 1, nchannels
    opdp_path%atm_level(1, j) = 0.0_jprb
!  Introduce ztemp to stop store followed by load in this recursive loop (DJS)
    ztemp = 0.0_jprb
    DO lev = 2, nlevels
      lay = lev - 1
      ztemp = ztemp + opticaldepth(lay, j)
      opdp_path%atm_level(lev, j) = ztemp
    ENDDO
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_opdep
