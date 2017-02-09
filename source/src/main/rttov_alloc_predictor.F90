!
SUBROUTINE rttov_alloc_predictor( &
            & err,         &
            & npredictors, &
            & predictors,  &
            & coef,        &
            & asw,         &
            & addsolar,    &
            & init)
! Description:
!   Allocation/deallocation of a predictors structure
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
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
!INTF_OFF
#include "throw.h"
!INTF_ON

  USE rttov_types, ONLY : rttov_coef, predictors_type
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE rttov_const, ONLY : sensor_id_mw
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim)   , INTENT(OUT)             :: err         ! return code
  INTEGER(KIND=jpim)   , INTENT(IN)              :: npredictors
  TYPE(predictors_type), INTENT(INOUT)           :: predictors
  TYPE(rttov_coef     ), INTENT(IN)              :: coef
  INTEGER(KIND=jpim)   , INTENT(IN)              :: asw         ! 1=allocate, 0=deallocate
  LOGICAL(KIND=jplm)   , INTENT(IN)              :: addsolar
  LOGICAL(KIND=jplm)   , INTENT(IN)   , OPTIONAL :: init
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_init_predictor.interface"

  INTEGER(KIND=jpim) :: nlayers
  LOGICAL(KIND=jplm) :: init1
!- End of header --------------------------------------------------------
  TRY
  init1 = .FALSE.
  IF (PRESENT(init)) init1 = init
  nlayers = coef%nlayers

!  Allocate section
  IF (asw .EQ. 1) THEN
    CALL nullify_struct()
    predictors%nlevels = coef%nlevels
    predictors%nmixed  = coef%nmixed
    predictors%nwater  = coef%nwater
    predictors%nozone  = coef%nozone
    predictors%nwvcont = coef%nwvcont
    predictors%nco2    = coef%nco2
    predictors%nn2o    = coef%nn2o
    predictors%nco     = coef%nco
    predictors%nch4    = coef%nch4
    predictors%ncloud  = 0
    predictors%npmc    = coef%pmc_nvar

    ALLOCATE (predictors%path1%mixedgas(coef%nmixed, nlayers, npredictors), STAT = err)
    THROWM(err .NE. 0, "allocation of predictors%path1%mixedgas")
    ALLOCATE (predictors%path1%watervapour(coef%nwater, nlayers, npredictors), STAT = err)
    THROWM(err .NE. 0, "allocation of predictors%path1%watervapour")
    IF (coef%id_sensor == sensor_id_mw) THEN
      ALLOCATE (predictors%path1%clw(nlayers, npredictors), STAT = err)
      THROWM(err .NE. 0, "allocation of predictors%path1%clw")
    ENDIF
    IF (coef%nozone > 0) THEN
      ALLOCATE (predictors%path1%ozone(coef%nozone, nlayers, npredictors), STAT = err)
      THROWM(err .NE. 0, "allocation of predictors%path1%ozone")
    ENDIF
    IF (coef%nwvcont > 0) THEN
      ALLOCATE (predictors%path1%wvcont(coef%nwvcont, nlayers, npredictors), STAT = err)
      THROWM(err .NE. 0, "allocation of predictors%path1%wvcont")
    ENDIF
    IF (coef%nco2 > 0) THEN
      ALLOCATE (predictors%path1%co2(coef%nco2, nlayers, npredictors), STAT = err)
      THROWM(err .NE. 0, "allocation of predictors%path1%co2")
    ENDIF
    IF (coef%nn2o > 0) THEN
      ALLOCATE (predictors%path1%n2o(coef%nn2o, nlayers, npredictors), STAT = err)
      THROWM(err .NE. 0, "allocation of predictors%path1%n2o")
    ENDIF
    IF (coef%nco > 0) THEN
      ALLOCATE (predictors%path1%co(coef%nco, nlayers, npredictors), STAT = err)
      THROWM(err .NE. 0, "allocation of predictors%path1%co")
    ENDIF
    IF (coef%nch4 > 0) THEN
      ALLOCATE (predictors%path1%ch4(coef%nch4, nlayers, npredictors), STAT = err)
      THROWM(err .NE. 0, "allocation of predictors%path1%ch4")
    ENDIF
    IF (coef%pmc_shift) THEN
      ALLOCATE (predictors%path1%pmc(coef%pmc_nvar, nlayers+1, npredictors, coef%fmv_chn), STAT = err)
      THROWM(err .NE. 0, "allocation of predictors%path1%pmc")
    ENDIF

    IF (addsolar) THEN
      ALLOCATE (predictors%path2%mixedgas(coef%nmixed, nlayers, npredictors), STAT = err)
      THROWM(err .NE. 0, "allocation of predictors%path2%mixedgas")
      ALLOCATE (predictors%path2%watervapour(coef%nwater, nlayers, npredictors), STAT = err)
      THROWM(err .NE. 0, "allocation of predictors%path2%watervapour")
      IF (coef%nozone > 0) THEN
        ALLOCATE (predictors%path2%ozone(coef%nozone, nlayers, npredictors), STAT = err)
        THROWM(err .NE. 0, "allocation of predictors%path2%ozone")
      ENDIF
      IF (coef%nwvcont > 0) THEN
        ALLOCATE (predictors%path2%wvcont(coef%nwvcont, nlayers, npredictors), STAT = err)
        THROWM(err .NE. 0, "allocation of predictors%path2%wvcont")
      ENDIF
      IF (coef%nco2 > 0) THEN
        ALLOCATE (predictors%path2%co2(coef%nco2, nlayers, npredictors), STAT = err)
        THROWM(err .NE. 0, "allocation of predictors%path2%co2")
      ENDIF
      IF (coef%nn2o > 0) THEN
        ALLOCATE (predictors%path2%n2o(coef%nn2o, nlayers, npredictors), STAT = err)
        THROWM(err .NE. 0, "allocation of predictors%path2%n2o")
      ENDIF
      IF (coef%nco > 0) THEN
        ALLOCATE (predictors%path2%co(coef%nco, nlayers, npredictors), STAT = err)
        THROWM(err .NE. 0, "allocation of predictors%path2%co")
      ENDIF
      IF (coef%nch4 > 0) THEN
        ALLOCATE (predictors%path2%ch4(coef%nch4, nlayers, npredictors), STAT = err)
        THROWM(err .NE. 0, "allocation of predictors%path2%ch4")
      ENDIF
    ENDIF

    IF (init1) THEN
      CALL rttov_init_predictor(addsolar, predictors)
    ENDIF
  ENDIF

  IF (asw .EQ. 0) THEN
    DEALLOCATE (predictors%path1%mixedgas, STAT = err)
    THROWM(err .NE. 0, "deallocation of predictors%path1%mixedgas")
    DEALLOCATE (predictors%path1%watervapour, STAT = err)
    THROWM(err .NE. 0, "deallocation of predictors%path1%watervapour")
    IF (coef%id_sensor == sensor_id_mw) THEN
      DEALLOCATE (predictors%path1%clw, STAT = err)
      THROWM(err .NE. 0, "deallocation of predictors%path1%clw")
    ENDIF
    IF (coef%nozone > 0) THEN
      DEALLOCATE (predictors%path1%ozone, STAT = err)
      THROWM(err .NE. 0, "deallocation of predictors%path1%ozone")
    ENDIF
    IF (coef%nwvcont > 0) THEN
      DEALLOCATE (predictors%path1%wvcont, STAT = err)
      THROWM(err .NE. 0, "deallocation of predictors%path1%wvcont")
    ENDIF
    IF (coef%nco2 > 0) THEN
      DEALLOCATE (predictors%path1%co2, STAT = err)
      THROWM(err .NE. 0, "deallocation of predictors%path1%co2")
    ENDIF
    IF (coef%nn2o > 0) THEN
      DEALLOCATE (predictors%path1%n2o, STAT = err)
      THROWM(err .NE. 0, "deallocation of predictors%path1%n2o")
    ENDIF
    IF (coef%nco > 0) THEN
      DEALLOCATE (predictors%path1%co, STAT = err)
      THROWM(err .NE. 0, "deallocation of predictors%path1%co")
    ENDIF
    IF (coef%nch4 > 0) THEN
      DEALLOCATE (predictors%path1%ch4, STAT = err)
      THROWM(err .NE. 0, "deallocation of predictors%path1%ch4")
    ENDIF
    IF (coef%pmc_shift) THEN
      DEALLOCATE (predictors%path1%pmc, STAT = err)
      THROWM(err .NE. 0, "deallocation of predictors%path1%pmc")
    ENDIF

    IF (addsolar) THEN
      DEALLOCATE (predictors%path2%mixedgas, STAT = err)
      THROWM(err .NE. 0, "deallocation of predictors%path2%mixedgas")
      DEALLOCATE (predictors%path2%watervapour, STAT = err)
      THROWM(err .NE. 0, "deallocation of predictors%path2%watervapour")
      IF (coef%nozone > 0) THEN
        DEALLOCATE (predictors%path2%ozone, STAT = err)
        THROWM(err .NE. 0, "deallocation of predictors%path2%ozone")
      ENDIF
      IF (coef%nwvcont > 0) THEN
        DEALLOCATE (predictors%path2%wvcont, STAT = err)
        THROWM(err .NE. 0, "deallocation of predictors%path2%wvcont")
      ENDIF
      IF (coef%nco2 > 0) THEN
        DEALLOCATE (predictors%path2%co2, STAT = err)
        THROWM(err .NE. 0, "deallocation of predictors%path2%co2")
      ENDIF
      IF (coef%nn2o > 0) THEN
        DEALLOCATE (predictors%path2%n2o, STAT = err)
        THROWM(err .NE. 0, "deallocation of predictors%path2%n2o")
      ENDIF
      IF (coef%nco > 0) THEN
        DEALLOCATE (predictors%path2%co, STAT = err)
        THROWM(err .NE. 0, "deallocation of predictors%path2%co")
      ENDIF
      IF (coef%nch4 > 0) THEN
        DEALLOCATE (predictors%path2%ch4, STAT = err)
        THROWM(err .NE. 0, "deallocation of predictors%path2%ch4")
      ENDIF
    ENDIF

    CALL nullify_struct()

  ENDIF
  CATCH
CONTAINS
  SUBROUTINE nullify_struct()
    NULLIFY (predictors%path1%mixedgas)
    NULLIFY (predictors%path1%watervapour)
    NULLIFY (predictors%path1%ozone)
    NULLIFY (predictors%path1%wvcont)
    NULLIFY (predictors%path1%co2)
    NULLIFY (predictors%path1%n2o)
    NULLIFY (predictors%path1%co)
    NULLIFY (predictors%path1%ch4)
    NULLIFY (predictors%path1%clw)
    NULLIFY (predictors%path1%pmc)
    NULLIFY (predictors%path2%mixedgas)
    NULLIFY (predictors%path2%watervapour)
    NULLIFY (predictors%path2%ozone)
    NULLIFY (predictors%path2%wvcont)
    NULLIFY (predictors%path2%co2)
    NULLIFY (predictors%path2%n2o)
    NULLIFY (predictors%path2%co)
    NULLIFY (predictors%path2%ch4)
  END SUBROUTINE nullify_struct
END SUBROUTINE rttov_alloc_predictor
