!
SUBROUTINE rttov_init_trans_scatt_ir(transmission_scatt_ir)
! Description:
! allocation/deallocation of a  transmission_scatt_ir structure
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
!    Copyright 2007, EUMETSAT, All Rights Reserved.
!
! Method:
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0       22/20/2002  Creation
!  1.1       03/11/2009  Transmittances / optical depths on levels (A Geer)
!  1.2       02/12/2009  Introduced new variables for the mixed cloud scheme (Marco Matricardi)
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
  USE rttov_types, ONLY : transmission_scatt_ir_Type
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(out):
  TYPE(transmission_scatt_ir_Type), INTENT(INOUT) :: transmission_scatt_ir
!INTF_END
  IF (ASSOCIATED(transmission_scatt_ir%OPDPAAER)    ) transmission_scatt_ir%OPDPAAER     = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%OPDPSAER)    ) transmission_scatt_ir%OPDPSAER     = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%GPARAERA)    ) transmission_scatt_ir%GPARAERA     = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%PHASINTUPREF)) transmission_scatt_ir%PHASINTUPREF = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%PHASINTDOREF)) transmission_scatt_ir%PHASINTDOREF = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%AZPHAERUP)   ) transmission_scatt_ir%AZPHAERUP    = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%AZPHAERDO)   ) transmission_scatt_ir%AZPHAERDO    = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%GPARAER)     ) transmission_scatt_ir%GPARAER      = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%AZPHAERUPA)  ) transmission_scatt_ir%AZPHAERUPA   = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%AZPHAERDOA)  ) transmission_scatt_ir%AZPHAERDOA   = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%OPDPAERLA)   ) transmission_scatt_ir%OPDPAERLA    = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%AZPHUP)      ) transmission_scatt_ir%AZPHUP       = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%AZPHDO)      ) transmission_scatt_ir%AZPHDO       = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%AZPHUPCLS)   ) transmission_scatt_ir%AZPHUPCLS    = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%AZPHDOCLS)   ) transmission_scatt_ir%AZPHDOCLS    = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%AZPHUPTOT)   ) transmission_scatt_ir%AZPHUPTOT    = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%AZPHDOTOT)   ) transmission_scatt_ir%AZPHDOTOT    = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%OPDPA)       ) transmission_scatt_ir%OPDPA        = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%OPDPS)       ) transmission_scatt_ir%OPDPS        = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%GPAR)        ) transmission_scatt_ir%GPAR         = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%GPARTOT)     ) transmission_scatt_ir%GPARTOT      = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%OPDPACLS)    ) transmission_scatt_ir%OPDPACLS     = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%OPDPSCLS)    ) transmission_scatt_ir%OPDPSCLS     = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%GPARCLS)     ) transmission_scatt_ir%GPARCLS      = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%OPDPCLDLA)   ) transmission_scatt_ir%OPDPCLDLA    = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%AZPHACUP)    ) transmission_scatt_ir%AZPHACUP     = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%AZPHACDO)    ) transmission_scatt_ir%AZPHACDO     = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%OPDPABS)     ) transmission_scatt_ir%OPDPABS      = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%OPDPSCA)     ) transmission_scatt_ir%OPDPSCA      = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%OPDPEXT)     ) transmission_scatt_ir%OPDPEXT      = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%OPDPACLSUN)  ) transmission_scatt_ir%OPDPACLSUN   = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%OPDPACSUN)   ) transmission_scatt_ir%OPDPACSUN    = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%OPDPACL)     ) transmission_scatt_ir%OPDPACL      = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%OPDPAC)      ) transmission_scatt_ir%OPDPAC       = 0._jprb
  IF (ASSOCIATED(transmission_scatt_ir%SSA)         ) transmission_scatt_ir%SSA          = 0._jprb
END SUBROUTINE 
