!
PROGRAM rttov_test_get_pc_predictindex
  ! Description:
  !   Test program for running rttov_get_pc_predictindex
  !   subroutine on commandline.
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
  ! Current Code Owner: SAF NWP
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  !
#include "throw.h"

  USE parkind1, ONLY : jpim
  USE rttov_types, ONLY : rttov_options
  USE rttov_getoptions
  USE rttov_unix_env, ONLY : rttov_iargc, rttov_exit

  IMPLICIT NONE

#include "rttov_errorreport.interface"
#include "rttov_get_pc_predictindex.interface"

  TYPE(rttov_options)  :: opts
  INTEGER(KIND=jpim)   :: ERR
  INTEGER(KIND=jpim)   :: ipcbnd, ipcreg
  CHARACTER(LEN=256)   :: f_pccoef_in = ""
  INTEGER(KIND=jpim), POINTER ::predictindex(:)

  !- End of header --------------------------------------------------------
TRY

  IF (rttov_iargc() == 0) THEN
    PRINT *, "Usage:  --pccoef-in  ...  --ipcbnd ... --ipcreg ... "
    STOP
  ENDIF

  CALL getoption( "--pccoef-in",  f_pccoef_in  )
  CALL getoption( "--ipcbnd",  ipcbnd )
  CALL getoption( "--ipcreg",  ipcreg )

  opts%rt_ir%pc%addpc  = f_pccoef_in .NE. ""
  opts%rt_ir%pc%ipcbnd = ipcbnd
  opts%rt_ir%pc%ipcreg = ipcreg

  CALL rttov_get_pc_predictindex( &
            & err,           &
            & opts,          &
            & predictindex,  &
            & file_pccoef   = f_pccoef_in)
  THROW( ERR .ne. 0)

  WRITE(*,*) "Number of channels ", SIZE(predictindex)
  WRITE(*,*) predictindex
  DEALLOCATE(predictindex)

PCATCH

END PROGRAM rttov_test_get_pc_predictindex
