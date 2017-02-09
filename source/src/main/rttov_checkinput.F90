!
SUBROUTINE rttov_checkinput( &
       & opts,        &
       & prof,        &
       & prof_user,   &
       & coef,        &
       & coef_pccomp, &
       & err          )
  ! Description:
  ! Check input profile/angles
  ! (i)  Are physically realistic
  ! (ii) Profile values are within the basis set used to
  !      generate the coefficients
  ! Unphysical values return a fatal error status
  ! Profile values outside the basis result in warnings being printed
  !
  ! Since the comparisons are made on coefficient levels
  ! gas abundances are compared in units of ppmv over dry air
  ! and where limits are exceeded values are reported in units
  ! of ppmv over dry air.
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
  ! Method: Check input profiles with fixed limits specified
  !         in constants and coeff file.
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  1.0       01/12/2002  New F90 code with structures (P Brunel A Smith)
  !  1.1       02/01/2003  More comments added (R Saunders)
  !  1.2       29/01/2003  More tests and add CO2 (P Brunel)
  !  1.3       27/06/2005  Uncommented water vapor and ozone profile checks  (R Saunders)
  !  1.4       23/01/2006  Changes from Marco (R Saunders)
  !  1.5       25/09/2007  Introduce a maximum number of warning messages (P Brunel)
  !  1.6       07/12/2007  Remove above; replace with a logical global variable... P.B.
  !  1.7       12/12/2007  Added limit checking for CH4, N2O, CO (R Saunders)
  !  1.8       16/01/2008  Added facility to apply regression limits (N Bormann)
  !  1.9       02/12/2009  Marco Matricardi: Added principal components
  ! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:

!INTF_OFF
#include "throw.h"
!INTF_ON

  USE rttov_types, ONLY : &
      rttov_coef,        &
      rttov_options,     &
      rttov_coef_pccomp, &
      profile_type

  USE parkind1, ONLY : jpim

!INTF_OFF
  USE rttov_const, ONLY : &
      nsurftype,          &
      nwatertype,         &
      surftype_sea,       &
      gas_id_watervapour, &
      gas_id_ozone,       &
      gas_id_co2,         &
      gas_id_co,          &
      gas_id_ch4,         &
      gas_id_n2o,         &
      tmax, tmin,         &
      qmax, qmin,         &
      o3max, o3min,       &
      co2max, co2min,     &
      comax, comin,       &
      n2omax, n2omin,     &
      ch4max, ch4min,     &
      clwmax, clwmin,     &
      pmax, pmin,         &
      wmax,               &
      zenmax, zenmaxv9,   &
      ctpmax, ctpmin,     &
      bemax, bemin,       &
      nish, nidg
  USE yomhook, ONLY : &
      LHOOK,&
      DR_HOOK
  USE parkind1, ONLY : &
      jprb, jplm
!INTF_ON

  IMPLICIT NONE

  TYPE(rttov_options),     INTENT(IN)    :: opts
  TYPE(profile_type),      INTENT(IN)    :: prof(:)
  TYPE(profile_type),      INTENT(IN)    :: prof_user(:)
  TYPE(rttov_coef),        INTENT(IN)    :: coef
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  INTEGER(KIND=jpim),      INTENT(OUT)   :: err

!INTF_END

#include "rttov_errorreport.interface"

  !local variables:
  REAL(KIND=jprb)    :: wind, zmax
  REAL(KIND=jprb)    :: dp(coef%nlevels)
  INTEGER(KIND=jpim) :: firstlevel, firstuserlevel, toplevel
  INTEGER(KIND=jpim) :: ig
  INTEGER(KIND=jpim) :: nprofiles, iprof
  CHARACTER(32)      :: sprof
  CHARACTER(128)     :: msg
  LOGICAL(KIND=jplm) :: ltest(SIZE(prof(1)%p))
  REAL(KIND=jprb)    :: ZHOOK_HANDLE

  !- End of header --------------------------------------------------------
TRY

  !-------------
  ! Initialize
  !-------------

IF (LHOOK) CALL DR_HOOK('RTTOV_CHECKINPUT',0_jpim,ZHOOK_HANDLE)

nprofiles = SIZE(prof)

DO iprof = 1, nprofiles
! For OpenMP only: Fortran I/O is not thread-safe
!$OMP CRITICAL
  WRITE(sprof,'(" (profile number = ",I8,")")') iprof
!$OMP END CRITICAL

  ! Compare profile levels and model levels
  IF (prof(iprof)%nlevels /= coef%nlevels) THEN
    err = errorstatus_fatal
    THROWM( err .NE. 0 , "invalid profile number of levels"//sprof)
  ENDIF

  dp(:) = ABS(coef%ref_prfl_p(:) - prof(iprof)%p(:)) / coef%ref_prfl_p(:)
  IF (ANY( dp > 0.01_jprb ) ) THEN
    err = errorstatus_fatal
    THROWM( err .NE. 0 , "invalid profile pressure levels"//sprof)
  ENDIF


  ! We should check all levels which can contribute to output radiances.
  ! Exactly which coef levels contribute depends on the user and coef
  ! levels, the surface pressure, and the interpolation mode used.
  ! We would like to avoid applying checks to levels which may contain
  ! extrapolated values below the surface.
  ! In general we should check down at least as far as the *input*
  ! pressure level which lies on or below the surface pressure.

  ! Find first user level at or below surface
  DO firstuserlevel = prof_user(iprof)%nlevels, 2, -1
    IF (prof_user(iprof)%p(firstuserlevel-1) < prof_user(iprof)%s2m%p) EXIT
  ENDDO

  ! Find first coef level at or below firstuserlevel
  DO firstlevel = prof(iprof)%nlevels, 2, -1
    IF (prof(iprof)%p(firstlevel-1) < prof_user(iprof)%p(firstuserlevel)) EXIT
  ENDDO


  ! When profile extrapolation is enabled we should ignore the coef levels
  ! above the top of the user profile because these will be clipped to the
  ! regression limits by rttov_apply_reg_limits

  IF (opts%interpolation%reg_limit_extrap .AND. &
      prof(iprof)%p(2) < prof_user(iprof)%p(1)) THEN

      ! Determine first coef level below user input profile top level
      DO toplevel = 3, prof(iprof)%nlevels
        IF (prof(iprof)%p(toplevel) > prof_user(iprof)%p(1)) EXIT
      ENDDO

  ELSE
    toplevel = 1
  ENDIF


  !------------------------------
  ! Check for unphysical values
  !------------------------------
  ! zenith angle
  IF (coef%fmv_model_ver == 9 ) THEN
    zmax = zenmaxv9
  Else
    zmax = zenmax
  ENDIF
  IF (prof(iprof)%zenangle > zmax .OR. &
      prof(iprof)%zenangle < 0._jprb) THEN
    err = errorstatus_fatal
    THROWM( err .NE. 0 , "invalid zenith angle"//sprof)
  ENDIF

  ! Solar zenith angle
  IF (opts%rt_ir%addsolar) THEN
    IF (prof(iprof)%sunzenangle < 0._jprb) THEN
      err = errorstatus_fatal
      THROWM( err .NE. 0 , "invalid solar zenith angle"//sprof)
    ENDIF
  ENDIF

  ! Cloud Fraction
  IF (prof(iprof)%cfraction > 1._jprb .OR. &
      prof(iprof)%cfraction < 0._jprb) THEN
    err = errorstatus_fatal
    THROWM( err .NE. 0 , "invalid cloud fraction"//sprof)
  ENDIF

  ! Cloud Top Pressure
  IF (prof(iprof)%cfraction .NE. 0 ) THEN
    IF (prof(iprof)%ctp > ctpmax .OR. &
        prof(iprof)%ctp < ctpmin) THEN
      err = errorstatus_fatal
      THROWM( err .NE. 0 , "invalid cloud top pressure"//sprof)
    ENDIF
  ENDIF

  ! Zeeman variables
  IF (coef%inczeeman) THEN
    ! Magnetic field strength
    IF (prof(iprof)%be > bemax .OR. &
        prof(iprof)%be < bemin) THEN
      err = errorstatus_fatal
      THROWM( err .NE. 0 , "invalid magnetic field strength"//sprof)
    ENDIF

    ! Cosine of angle between path and mag. field
    IF (prof(iprof)%cosbk > 1._jprb  .OR. &
        prof(iprof)%cosbk < -1._jprb) THEN
      err = errorstatus_fatal
      THROWM( err .NE. 0 , "invalid cosbk"//sprof)
    ENDIF
  ENDIF

  !1.1 surface variables
  !---------------------

  ! Pressure
  IF (prof(iprof)%s2m%p > pmax .OR. &
      prof(iprof)%s2m%p < pmin) THEN
    err = errorstatus_fatal
    THROWM( err .NE. 0 , "invalid surface pressure"//sprof)
  ENDIF

  ! 2m air temperature
  IF (prof(iprof)%s2m%t > tmax .OR. &
      prof(iprof)%s2m%t < tmin) THEN
    err = errorstatus_fatal
    THROWM( err .NE. 0 , "invalid 2m air temperature"//sprof)
  ENDIF

  ! 2m water vapour - only used IF opts%rt_all%use_q2m is TRUE
  IF (opts%rt_all%use_q2m ) THEN
    IF (prof(iprof)%s2m%q > qmax .OR. &
        prof(iprof)%s2m%q < qmin) THEN
      err = errorstatus_fatal
      THROWM( err .NE. 0 , "invalid 2m water vapour"//sprof)
    ENDIF
  ENDIF

  !  surface wind speed
  wind = SQRT(prof(iprof)%s2m%u * prof(iprof)%s2m%u + &
              prof(iprof)%s2m%v * prof(iprof)%s2m%v)
  IF (wind > wmax .OR. &
      wind < 0._jprb      ) THEN
    err = errorstatus_fatal
    THROWM( err .NE. 0 , "invalid 10m wind speed"//sprof)
  ENDIF

  IF (opts%rt_ir%addsolar) THEN
    IF (prof(iprof)%s2m%wfetc <= 0._jprb .AND. &
        prof(iprof)%skin%surftype == surftype_sea) THEN
      err = errorstatus_fatal
      THROWM( err .NE. 0 , "invalid wfetc"//sprof)
    ENDIF
  ENDIF

  ! surface skin temperature
  IF (prof(iprof)%skin%t > tmax .OR. &
      prof(iprof)%skin%t < tmin) THEN
    err = errorstatus_fatal
    THROWM( err .NE. 0 , "invalid skin surface temperature"//sprof)
  ENDIF

  ! surface type
  IF (prof(iprof)%skin%surftype < 0 .OR. &
      prof(iprof)%skin%surftype > nsurftype) THEN
    err = errorstatus_fatal
    THROWM( err .NE. 0 , "invalid surface type"//sprof)
  ENDIF

  ! water type
  IF (prof(iprof)%skin%watertype < 0 .OR. &
      prof(iprof)%skin%watertype > nwatertype) THEN
    err = errorstatus_fatal
    THROWM( err .NE. 0 , "invalid water type"//sprof)
  ENDIF

  ! Foam fraction
  IF (opts%rt_mw%supply_foam_fraction) THEN
    IF (prof(iprof)%skin%foam_fraction < 0._jprb .OR. &
        prof(iprof)%skin%foam_fraction > 1._jprb) THEN
      err = errorstatus_fatal
      THROWM(err .NE. 0, "invalid foam fraction")
    ENDIF
  ENDIF

  ! snow fraction
  IF (prof(iprof)%snow_frac < 0._jprb .OR. &
      prof(iprof)%snow_frac > 1._jprb) THEN
    err = errorstatus_fatal
    THROWM( err .NE. 0 , "invalid snow fraction"//sprof)
  ENDIF

  ! Atmospheric variables
  !-------------------------

  ! Predictors are calculated on *all* levels so check the hard
  ! limits on every level (to bottom of profile) to avoid errors.

  ltest = .FALSE.

  ! temperature
  ltest(toplevel:coef%nlevels) = (prof(iprof)%t(toplevel:coef%nlevels) > tmax)
  IF (ANY(ltest)) THEN
    err = errorstatus_fatal
    CALL print_info("Input temperature profile exceeds allowed maximum:", &
      (/tmax/), PACK(prof(iprof)%p(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)), &
                PACK(prof(iprof)%t(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)))
  ENDIF

  ltest(toplevel:coef%nlevels) = (prof(iprof)%t(toplevel:coef%nlevels) < tmin)
  IF (ANY(ltest)) THEN
    err = errorstatus_fatal
    CALL print_info("Input temperature profile exceeds allowed minimum:", &
      (/tmin/), PACK(prof(iprof)%p(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)), &
                PACK(prof(iprof)%t(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)))
  ENDIF
  THROWM( err .NE. 0 , "some invalid atmospheric temperature"//sprof)

  ! water vapour
  ltest(toplevel:coef%nlevels) = (prof(iprof)%q(toplevel:coef%nlevels) > qmax)
  IF (ANY(ltest)) THEN
    err = errorstatus_fatal
    CALL print_info("Input water vapour profile exceeds allowed maximum:", &
      (/qmax/), PACK(prof(iprof)%p(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)), &
                PACK(prof(iprof)%q(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)))
  ENDIF

  ltest(toplevel:coef%nlevels) = (prof(iprof)%q(toplevel:coef%nlevels) < qmin)
  IF (ANY(ltest)) THEN
    err = errorstatus_fatal
    CALL print_info("Input water vapour profile exceeds allowed minimum:", &
      (/qmin/), PACK(prof(iprof)%p(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)), &
                PACK(prof(iprof)%q(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)))
  ENDIF
  THROWM( err .NE. 0 , "some invalid atmospheric water vapour"//sprof)

  ! ozone
  IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN
    ltest(toplevel:coef%nlevels) = (prof(iprof)%o3(toplevel:coef%nlevels) > o3max)
    IF (ANY(ltest)) THEN
      err = errorstatus_fatal
      CALL print_info("Input ozone profile exceeds allowed maximum:", &
        (/o3max/), PACK(prof(iprof)%p(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)), &
                   PACK(prof(iprof)%o3(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)))
    ENDIF

    ltest(toplevel:coef%nlevels) = (prof(iprof)%o3(toplevel:coef%nlevels) < o3min)
    IF (ANY(ltest)) THEN
      err = errorstatus_fatal
      CALL print_info("Input ozone profile exceeds allowed minimum:", &
        (/o3min/), PACK(prof(iprof)%p(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)), &
                    PACK(prof(iprof)%o3(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)))
    ENDIF
    THROWM( err .NE. 0 , "some invalid atmospheric ozone"//sprof)
  ENDIF

  ! CO2
  IF (opts%rt_ir%co2_Data .AND. coef%nco2 > 0) THEN
    ltest(toplevel:coef%nlevels) = (prof(iprof)%co2(toplevel:coef%nlevels) > co2max)
    IF (ANY(ltest)) THEN
      err = errorstatus_fatal
      CALL print_info("Input CO2 profile exceeds allowed maximum:", &
        (/co2max/), PACK(prof(iprof)%p(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)), &
                    PACK(prof(iprof)%co2(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)))
    ENDIF

    ltest(toplevel:coef%nlevels) = (prof(iprof)%co2(toplevel:coef%nlevels) < co2min)
    IF (ANY(ltest)) THEN
      err = errorstatus_fatal
      CALL print_info("Input CO2 profile exceeds allowed minimum:", &
        (/co2min/), PACK(prof(iprof)%p(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)), &
                    PACK(prof(iprof)%co2(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)))
    ENDIF
    THROWM( err .NE. 0 , "some invalid atmospheric CO2"//sprof)
  ENDIF

  ! CO
  IF (opts%rt_ir%co_Data .AND. coef%nco > 0) THEN
    ltest(toplevel:coef%nlevels) = (prof(iprof)%co(toplevel:coef%nlevels) > comax)
    IF (ANY(ltest)) THEN
      err = errorstatus_fatal
      CALL print_info("Input CO profile exceeds allowed maximum:", &
        (/comax/), PACK(prof(iprof)%p(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)), &
                   PACK(prof(iprof)%co(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)))
    ENDIF

    ltest(toplevel:coef%nlevels) = (prof(iprof)%co(toplevel:coef%nlevels) < comin)
    IF (ANY(ltest)) THEN
      err = errorstatus_fatal
      CALL print_info("Input CO profile exceeds allowed minimum:", &
        (/comin/), PACK(prof(iprof)%p(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)), &
                   PACK(prof(iprof)%co(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)))
    ENDIF
    THROWM( err .NE. 0 , "some invalid atmospheric CO"//sprof)
  ENDIF

  ! N2O
  IF (opts%rt_ir%n2o_Data .AND. coef%nn2o > 0) THEN
    ltest(toplevel:coef%nlevels) = (prof(iprof)%n2o(toplevel:coef%nlevels) > n2omax)
    IF (ANY(ltest)) THEN
      err = errorstatus_fatal
      CALL print_info("Input N2O profile exceeds allowed maximum:", &
        (/n2omax/), PACK(prof(iprof)%p(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)), &
                    PACK(prof(iprof)%n2o(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)))
    ENDIF

    ltest(toplevel:coef%nlevels) = (prof(iprof)%n2o(toplevel:coef%nlevels) < n2omin)
    IF (ANY(ltest)) THEN
      err = errorstatus_fatal
      CALL print_info("Input N2O profile exceeds allowed minimum:", &
        (/n2omin/), PACK(prof(iprof)%p(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)), &
                    PACK(prof(iprof)%n2o(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)))
    ENDIF
    THROWM( err .NE. 0 , "some invalid atmospheric N2O"//sprof)
  ENDIF

  ! CH4
  IF (opts%rt_ir%ch4_Data .AND. coef%nch4 > 0) THEN
    ltest(toplevel:coef%nlevels) = (prof(iprof)%ch4(toplevel:coef%nlevels) > ch4max)
    IF (ANY(ltest)) THEN
      err = errorstatus_fatal
      CALL print_info("Input CH4 profile exceeds allowed maximum:", &
        (/ch4max/), PACK(prof(iprof)%p(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)), &
                    PACK(prof(iprof)%ch4(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)))
    ENDIF

    ltest(toplevel:coef%nlevels) = (prof(iprof)%ch4(toplevel:coef%nlevels) < ch4min)
    IF (ANY(ltest)) THEN
      err = errorstatus_fatal
      CALL print_info("Input CH4 profile exceeds allowed minimum:", &
        (/ch4min/), PACK(prof(iprof)%p(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)), &
                    PACK(prof(iprof)%ch4(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)))
    ENDIF
    THROWM( err .NE. 0 , "some invalid atmospheric CH4"//sprof)
  ENDIF

  ! cloud liquid water
  IF (opts%rt_mw%clw_Data ) THEN
    ltest(toplevel:coef%nlevels) = (prof(iprof)%clw(toplevel:coef%nlevels) > clwmax)
    IF (ANY(ltest)) THEN
      err = errorstatus_fatal
      CALL print_info("Input cloud liquid water profile exceeds allowed maximum:", &
        (/clwmax/), PACK(prof(iprof)%p(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)), &
                    PACK(prof(iprof)%clw(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)))
    ENDIF

    ltest(toplevel:coef%nlevels) = (prof(iprof)%clw(toplevel:coef%nlevels) < clwmin)
    IF (ANY(ltest)) THEN
      err = errorstatus_fatal
      CALL print_info("Input cloud liquid water profile exceeds allowed minimum:", &
        (/clwmin/), PACK(prof(iprof)%p(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)), &
                    PACK(prof(iprof)%clw(toplevel:coef%nlevels), mask=ltest(toplevel:coef%nlevels)))
    ENDIF
    THROWM( err .NE. 0 , "some invalid atmospheric cloud liquid water"//sprof)
  ENDIF


  ! Cloud parameterisations
  ! Cannot check cloud/cfrac/aerosol profiles as they are not interpolated onto coef levels
  IF (opts%rt_ir%addclouds .AND. .NOT. opts%rt_ir%user_cld_opt_param ) THEN
    IF (prof(iprof)%ish < 1_jpim .OR. prof(iprof)%ish > nish ) THEN
      err = errorstatus_fatal
      THROWM( err .NE. 0 , "invalid ice particle shape (ish)"//sprof)
    ENDIF

    IF (prof(iprof)%ish < 3 ) THEN ! Not Baran scheme
      IF (prof(iprof)%idg < 1_jpim .OR. prof(iprof)%idg > nidg ) THEN
        err = errorstatus_fatal
        THROWM( err .NE. 0 , "invalid ice particle parameterisation (idg)"//sprof)
      ENDIF
    ENDIF
  ENDIF

  !---------------------------------
  ! Check against regression limits
  !---------------------------------

  IF (opts%config%verbose .AND. .NOT. opts%config%apply_reg_limits) THEN

    ltest = .FALSE.

    IF(.NOT.opts%rt_ir%pc%addpc)then

      ltest(toplevel:firstlevel) = (prof(iprof)%t(toplevel:firstlevel) > coef%lim_prfl_tmax(toplevel:firstlevel))
      IF (ANY(ltest)) THEN
        CALL print_info("Input temperature profile exceeds upper coef limit"//sprof, &
          PACK(coef%lim_prfl_tmax(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
          PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
          PACK(prof(iprof)%t(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
      ENDIF

      ltest(toplevel:firstlevel) = (prof(iprof)%t(toplevel:firstlevel) < coef%lim_prfl_tmin(toplevel:firstlevel))
      IF (ANY(ltest)) THEN
        CALL print_info("Input temperature profile exceeds lower coef limit"//sprof, &
          PACK(coef%lim_prfl_tmin(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
          PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
          PACK(prof(iprof)%t(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
      ENDIF

      ig = coef%fmv_gas_pos( gas_id_watervapour )
      ltest(toplevel:firstlevel) = (prof(iprof)%q(toplevel:firstlevel) > coef%lim_prfl_gmax(toplevel:firstlevel, ig))
      IF (ANY(ltest)) THEN
        CALL print_info("Input water vapour profile exceeds upper coef limit"//sprof, &
          PACK(coef%lim_prfl_gmax(toplevel:firstlevel, ig), mask=ltest(toplevel:firstlevel)), &
          PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
          PACK(prof(iprof)%q(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
      ENDIF

      ltest(toplevel:firstlevel) = (prof(iprof)%q(toplevel:firstlevel) < coef%lim_prfl_gmin(toplevel:firstlevel, ig))
      IF (ANY(ltest)) THEN
        CALL print_info("Input water vapour profile exceeds lower coef limit"//sprof, &
          PACK(coef%lim_prfl_gmin(toplevel:firstlevel, ig), mask=ltest(toplevel:firstlevel)), &
          PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
          PACK(prof(iprof)%q(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
      ENDIF

      IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN
        ig = coef%fmv_gas_pos( gas_id_ozone )
        ltest(toplevel:firstlevel) = (prof(iprof)%o3(toplevel:firstlevel) > coef%lim_prfl_gmax(toplevel:firstlevel, ig))
        IF (ANY(ltest)) THEN
          CALL print_info("Input ozone profile exceeds upper coef limit"//sprof, &
            PACK(coef%lim_prfl_gmax(toplevel:firstlevel, ig), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%o3(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
        ENDIF

        ltest(toplevel:firstlevel) = (prof(iprof)%o3(toplevel:firstlevel) < coef%lim_prfl_gmin(toplevel:firstlevel, ig))
        IF (ANY(ltest)) THEN
          CALL print_info("Input ozone profile exceeds lower coef limit"//sprof, &
            PACK(coef%lim_prfl_gmin(toplevel:firstlevel, ig), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%o3(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
        ENDIF
      ENDIF

      IF (opts%rt_ir%co2_Data .AND. coef%nco2 > 0) THEN
        ig = coef%fmv_gas_pos( gas_id_co2 )
        ltest(toplevel:firstlevel) = (prof(iprof)%co2(toplevel:firstlevel) > coef%lim_prfl_gmax(toplevel:firstlevel, ig))
        IF (ANY(ltest)) THEN
          CALL print_info("Input CO2 profile exceeds upper coef limit"//sprof, &
            PACK(coef%lim_prfl_gmax(toplevel:firstlevel, ig), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%co2(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
        ENDIF

        ltest(toplevel:firstlevel) = (prof(iprof)%co2(toplevel:firstlevel) < coef%lim_prfl_gmin(toplevel:firstlevel, ig))
        IF (ANY(ltest)) THEN
          CALL print_info("Input CO2 profile exceeds lower coef limit"//sprof, &
            PACK(coef%lim_prfl_gmin(toplevel:firstlevel, ig), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%co2(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
        ENDIF
      ENDIF

      IF (opts%rt_ir%co_Data .AND. coef%nco > 0) THEN
        ig = coef%fmv_gas_pos( gas_id_co )
        ltest(toplevel:firstlevel) = (prof(iprof)%co(toplevel:firstlevel) > coef%lim_prfl_gmax(toplevel:firstlevel, ig))
        IF (ANY(ltest)) THEN
          CALL print_info("Input CO profile exceeds upper coef limit"//sprof, &
            PACK(coef%lim_prfl_gmax(toplevel:firstlevel, ig), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%co(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
        ENDIF

        ltest(toplevel:firstlevel) = (prof(iprof)%co(toplevel:firstlevel) < coef%lim_prfl_gmin(toplevel:firstlevel, ig))
        IF (ANY(ltest)) THEN
          CALL print_info("Input CO profile exceeds lower coef limit"//sprof, &
            PACK(coef%lim_prfl_gmin(toplevel:firstlevel, ig), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%co(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
        ENDIF
      ENDIF

      IF (opts%rt_ir%n2o_Data .AND. coef%nn2o > 0) THEN
        ig = coef%fmv_gas_pos( gas_id_n2o )
        ltest(toplevel:firstlevel) = (prof(iprof)%n2o(toplevel:firstlevel) > coef%lim_prfl_gmax(toplevel:firstlevel, ig))
        IF (ANY(ltest)) THEN
          CALL print_info("Input N2O profile exceeds upper coef limit"//sprof, &
            PACK(coef%lim_prfl_gmax(toplevel:firstlevel, ig), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%n2o(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
        ENDIF

        ltest(toplevel:firstlevel) = (prof(iprof)%n2o(toplevel:firstlevel) < coef%lim_prfl_gmin(toplevel:firstlevel, ig))
        IF (ANY(ltest)) THEN
          CALL print_info("Input N2O profile exceeds lower coef limit"//sprof, &
            PACK(coef%lim_prfl_gmin(toplevel:firstlevel, ig), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%n2o(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
        ENDIF
      ENDIF

      IF (opts%rt_ir%ch4_Data .AND. coef%nch4 > 0) THEN
        ig = coef%fmv_gas_pos( gas_id_ch4 )
        ltest(toplevel:firstlevel) = (prof(iprof)%ch4(toplevel:firstlevel) > coef%lim_prfl_gmax(toplevel:firstlevel, ig))
        IF (ANY(ltest)) THEN
          CALL print_info("Input CH4 profile exceeds upper coef limit"//sprof, &
            PACK(coef%lim_prfl_gmax(toplevel:firstlevel, ig), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%ch4(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
        ENDIF

        ltest(toplevel:firstlevel) = (prof(iprof)%ch4(toplevel:firstlevel) < coef%lim_prfl_gmin(toplevel:firstlevel, ig))
        IF (ANY(ltest)) THEN
          CALL print_info("Input CH4 profile exceeds lower coef limit"//sprof, &
            PACK(coef%lim_prfl_gmin(toplevel:firstlevel, ig), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%ch4(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
        ENDIF
      ENDIF

    ELSE ! addpc

      ltest(toplevel:firstlevel) = (prof(iprof)%t(toplevel:firstlevel) > coef_pccomp%lim_pc_prfl_tmax(toplevel:firstlevel))
      IF (ANY(ltest)) THEN
        CALL print_info("PC-RTTOV: Input temperature profile exceeds upper coef limit"//sprof, &
          PACK(coef_pccomp%lim_pc_prfl_tmax(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
          PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
          PACK(prof(iprof)%t(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
      ENDIF

      ltest(toplevel:firstlevel) = (prof(iprof)%t(toplevel:firstlevel) < coef_pccomp%lim_pc_prfl_tmin(toplevel:firstlevel))
      IF (ANY(ltest)) THEN
        CALL print_info("PC-RTTOV: Input temperature profile exceeds lower coef limit"//sprof, &
          PACK(coef_pccomp%lim_pc_prfl_tmin(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
          PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
          PACK(prof(iprof)%t(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
      ENDIF

      ig = coef%fmv_gas_pos( gas_id_watervapour )
      ltest(toplevel:firstlevel) = (prof(iprof)%q(toplevel:firstlevel) > coef_pccomp%lim_pc_prfl_qmax(toplevel:firstlevel))
      IF (ANY(ltest)) THEN
        CALL print_info("PC-RTTOV: Input water vapour profile exceeds upper coef limit"//sprof, &
          PACK(coef_pccomp%lim_pc_prfl_qmax(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
          PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
          PACK(prof(iprof)%q(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
      ENDIF

      ltest(toplevel:firstlevel) = (prof(iprof)%q(toplevel:firstlevel) < coef_pccomp%lim_pc_prfl_qmin(toplevel:firstlevel))
      IF (ANY(ltest)) THEN
        CALL print_info("PC-RTTOV: Input water vapour profile exceeds lower coef limit"//sprof, &
          PACK(coef_pccomp%lim_pc_prfl_qmin(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
          PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
          PACK(prof(iprof)%q(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
      ENDIF

      IF (opts%rt_ir%ozone_Data .AND. coef%nozone > 0) THEN
        ig = coef%fmv_gas_pos( gas_id_ozone )
        ltest(toplevel:firstlevel) = (prof(iprof)%o3(toplevel:firstlevel) > coef_pccomp%lim_pc_prfl_ozmax(toplevel:firstlevel))
        IF (ANY(ltest)) THEN
          CALL print_info("PC-RTTOV: Input ozone profile exceeds upper coef limit"//sprof, &
            PACK(coef_pccomp%lim_pc_prfl_ozmax(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%o3(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
        ENDIF

        ltest(toplevel:firstlevel) = (prof(iprof)%o3(toplevel:firstlevel) < coef_pccomp%lim_pc_prfl_ozmin(toplevel:firstlevel))
        IF (ANY(ltest)) THEN
          CALL print_info("PC-RTTOV: Input ozone profile exceeds lower coef limit"//sprof, &
            PACK(coef_pccomp%lim_pc_prfl_ozmin(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%o3(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
        ENDIF
      ENDIF

      IF ((prof(iprof)%s2m%p < coef_pccomp%lim_pc_prfl_pmin) .OR. &
          (prof(iprof)%s2m%p > coef_pccomp%lim_pc_prfl_pmax)) THEN
        msg = "PC-RTTOV: surface pressure outside limits"//sprof
        WARN(msg)
      ENDIF

      IF ((prof(iprof)%s2m%t < coef_pccomp%lim_pc_prfl_tsmin) .OR. &
          (prof(iprof)%s2m%t > coef_pccomp%lim_pc_prfl_tsmax)) THEN
        msg = "PC-RTTOV: surface temperature outside limits"//sprof
        WARN(msg)
      ENDIF

      IF ((prof(iprof)%skin%t < coef_pccomp%lim_pc_prfl_skmin) .OR.   &
          (prof(iprof)%skin%t > coef_pccomp%lim_pc_prfl_skmax)) THEN
        msg = "PC-RTTOV: skin temperature outside limits"//sprof
        WARN(msg)
      ENDIF

      ! wind was already calculated above
      IF ((wind < coef_pccomp%lim_pc_prfl_wsmin) .OR.  &
          (wind > coef_pccomp%lim_pc_prfl_wsmax)) THEN
        msg = "PC-RTTOV: 10m wind speed outside limits"//sprof
        WARN(msg)
      ENDIF

    ENDIF ! addpc
  ENDIF !verbose and not apply_reg_limits
ENDDO

IF (LHOOK) CALL DR_HOOK('RTTOV_CHECKINPUT',1_jpim,ZHOOK_HANDLE)

CATCH

IF (LHOOK) CALL DR_HOOK('RTTOV_CHECKINPUT',1_jpim,ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE print_info(msg1, limits, levels, values)
    CHARACTER(LEN=*), INTENT(IN) :: msg1
    REAL(KIND=jprb),  INTENT(IN) :: limits(:)
    REAL(KIND=jprb),  INTENT(IN) :: levels(:)
    REAL(KIND=jprb),  INTENT(IN) :: values(:)

    CHARACTER(LEN=256) :: msg2
    INTEGER(KIND=jpim) :: imax, lmax

    imax = MIN(10, SIZE(levels))
    lmax = MIN(imax, SIZE(limits))

! Replace warn/info macros from throw.h with in-line code because
! NAG v5.3 complains otherwise.

! For OpenMP only: Fortran I/O is not thread-safe
!$OMP CRITICAL
!     WARN(msg1)
    CALL rttov_errorreport(errorstatus_success, TRIM(msg1), 'rttov_checkinput.F90')
    WRITE(msg2, '(a,10f10.4)') 'Limit   = ',limits(1:lmax)
!     INFO(TRIM(msg2))
    CALL rttov_errorreport(errorstatus_success, TRIM(msg2))
    WRITE(msg2, '(a,10f10.4)') 'p (hPa) = ',levels(1:imax)
!     INFO(TRIM(msg2))
    CALL rttov_errorreport(errorstatus_success, TRIM(msg2))
    WRITE(msg2, '(a,10f10.4)') 'Value   = ',values(1:imax)
!     INFO(TRIM(msg2))
    CALL rttov_errorreport(errorstatus_success, TRIM(msg2))
!$OMP END CRITICAL

  END SUBROUTINE print_info

END SUBROUTINE rttov_checkinput
