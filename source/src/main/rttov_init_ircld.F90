!
SUBROUTINE rttov_init_ircld(irclds)
! Description:
!   Initialise ircld structure
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
!  1.1       02/12/2009  Increased the dimensions of cldtyp (Marco Matricardi)
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
  USE rttov_types, ONLY : ircld_Type
!INTF_OFF
  USE parkind1, ONLY : jprb
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(out):
  TYPE(ircld_Type), INTENT(INOUT) :: irclds
!INTF_END
  irclds%xstrclr = 0._jprb
  IF (Associated(irclds%icldarr)   ) irclds%icldarr    = 0_jpim
  IF (Associated(irclds%xstrref1)  ) irclds%xstrref1   = 0._jprb
  IF (Associated(irclds%xstrref2)  ) irclds%xstrref2   = 0._jprb
  IF (Associated(irclds%cldtyp)    ) irclds%cldtyp     = 0_jpim
  IF (Associated(irclds%indexstr)  ) irclds%indexstr   = 0_jpim
  IF (Associated(irclds%icount1ref)) irclds%icount1ref = 0_jpim
  IF (Associated(irclds%iloopin)   ) irclds%iloopin    = 0_jpim
  IF (Associated(irclds%iflag)     ) irclds%iflag      = 0_jpim
  IF (Associated(irclds%xstr)      ) irclds%xstr       = 0._jprb
  IF (Associated(irclds%xstrminref)) irclds%xstrminref = 0._jprb
  IF (Associated(irclds%xstrref)   ) irclds%xstrref    = 0._jprb
  IF (Associated(irclds%cldcfr)    ) irclds%cldcfr     = 0._jprb
  IF (Associated(irclds%maxcov)    ) irclds%maxcov     = 0._jprb
  IF (Associated(irclds%xstrmax)   ) irclds%xstrmax    = 0._jprb
  IF (Associated(irclds%xstrmin)   ) irclds%xstrmin    = 0._jprb
  IF (Associated(irclds%a)         ) irclds%a          = 0._jprb
  IF (Associated(irclds%ntotref)   ) irclds%ntotref    = 0._jprb
  IF (Associated(irclds%tave)      ) irclds%tave       = 0._jprb
  IF (Associated(irclds%wmixave)   ) irclds%wmixave    = 0._jprb
  IF (Associated(irclds%xpresave)  ) irclds%xpresave   = 0._jprb
  IF (Associated(irclds%ppv)       ) irclds%ppv        = 0._jprb
  IF (Associated(irclds%esw)       ) irclds%esw        = 0._jprb
  IF (Associated(irclds%esi)       ) irclds%esi        = 0._jprb
  IF (Associated(irclds%flag)      ) irclds%flag       = .FALSE.
END SUBROUTINE 
