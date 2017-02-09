SUBROUTINE rttov_init_predictor(addsolar, predictors)
! Description:
!   Initialise predictors structure
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
!    Copyright 2010, EUMETSAT, All Rights Reserved.
!
! Method:
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
  USE rttov_types, ONLY : predictors_type
  USE parkind1,    ONLY : jplm
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
  LOGICAL(KIND=jplm),    INTENT(IN)    :: addsolar
  TYPE(predictors_type), INTENT(INOUT) :: predictors
!INTF_END
  IF (Associated(predictors%path1%mixedgas)   ) predictors%path1%mixedgas    = 0._jprb
  IF (Associated(predictors%path1%watervapour)) predictors%path1%watervapour = 0._jprb
  IF (Associated(predictors%path1%ozone)      ) predictors%path1%ozone       = 0._jprb
  IF (Associated(predictors%path1%wvcont)     ) predictors%path1%wvcont      = 0._jprb
  IF (Associated(predictors%path1%co2)        ) predictors%path1%co2         = 0._jprb
  IF (Associated(predictors%path1%n2o)        ) predictors%path1%n2o         = 0._jprb
  IF (Associated(predictors%path1%co)         ) predictors%path1%co          = 0._jprb
  IF (Associated(predictors%path1%ch4)        ) predictors%path1%ch4         = 0._jprb
  IF (Associated(predictors%path1%clw)        ) predictors%path1%clw         = 0._jprb
  IF (Associated(predictors%path1%pmc)        ) predictors%path1%pmc         = 0._jprb
  
  IF (addsolar) THEN
    IF (Associated(predictors%path2%mixedgas)   ) predictors%path2%mixedgas    = 0._jprb
    IF (Associated(predictors%path2%watervapour)) predictors%path2%watervapour = 0._jprb
    IF (Associated(predictors%path2%ozone)      ) predictors%path2%ozone       = 0._jprb
    IF (Associated(predictors%path2%wvcont)     ) predictors%path2%wvcont      = 0._jprb
    IF (Associated(predictors%path2%co2)        ) predictors%path2%co2         = 0._jprb
    IF (Associated(predictors%path2%n2o)        ) predictors%path2%n2o         = 0._jprb
    IF (Associated(predictors%path2%co)         ) predictors%path2%co          = 0._jprb
    IF (Associated(predictors%path2%ch4)        ) predictors%path2%ch4         = 0._jprb
  ENDIF                                                                        
END SUBROUTINE 
