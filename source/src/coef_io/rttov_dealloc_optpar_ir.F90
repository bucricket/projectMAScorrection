!
SUBROUTINE rttov_dealloc_optpar_ir(ERR, optp)
! Description:
! de-allocation of a coefficient structure
! The allocation is done by the readcoef subroutine called by the user
! this subroutine should be called once per coef structure when
! all rttov calls are completed.
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
!  1.0       01/12/2002  New F90 code with structures (P Brunel A Smith)
!  1.1       03/05/2004  Add specific RTTOV8 CO2 variable (P. Brunel)
!  1.2       02/06/2004  Update for RTTOV8 coefS (P. Brunel)
!  1.3       01/06/2005  Marco Matricardi (ECMWF):
!               --       N2O,CO and CH4 variables added
!  1.4       26/04/2007  Cloud/aerosol variables added (R Saunders)
!  1.5       11/10/2007  Nullify unused coef pointers P.Marguinaud
!  1.6       23/94/2008  Add some initialisation of scalars (P. Brunel)
!  1.7       06/03/2009  Deafault now coef % IncTop = .true. (P.Rayer)
!  1.8       02/12/2009  Add principal components (M. Matricardi, ECMWF)
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
#include "throw.h"
  USE rttov_types, ONLY :  &
       & rttov_optpar_ir
  USE parkind1, ONLY : jpim
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(out):
  INTEGER(KIND=jpim)       , INTENT(OUT)   :: ERR          ! return code
  TYPE(rttov_optpar_ir    ), INTENT(INOUT) :: optp
!INTF_END

#include "rttov_errorreport.interface"
#include "rttov_dealloc_coef_scatt_ir.interface"
#include "rttov_nullify_optpar_ir.interface"

! Local Arrays and Scalars:
  INTEGER(KIND=jpim) :: n
!- End of header --------------------------------------------------------
  TRY

    IF (Associated (optp%optpaer)) THEN

      DO n = 1, size (optp%optpaer)

        CALL rttov_dealloc_coef_scatt_ir (err, optp%optpaer(n))

        THROW( ERR .NE. 0 )

      ENDDO

      DEALLOCATE (optp%optpaer, STAT = ERR)

      THROW( ERR .NE. 0 )

    ENDIF


    IF (Associated(optp%optpwcl)) THEN

      DO n = 1, size(optp%optpwcl)

        CALL rttov_dealloc_coef_scatt_ir (err, optp%optpwcl(n))

        THROW( ERR .NE. 0 )

      ENDDO

      DEALLOCATE (optp%optpwcl, STAT = ERR)

      THROW( ERR .NE. 0 )

    ENDIF


    IF (Associated (optp%optpicl)) THEN

      DO n = 1, size(optp%optpicl)

        CALL rttov_dealloc_coef_scatt_ir (err, optp%optpicl(n))

        THROW( ERR .NE. 0 )

      ENDDO

      DEALLOCATE (optp%optpicl, STAT = ERR)

      THROW( ERR .NE. 0 )

    ENDIF


    IF (Associated (optp%optpiclb)) THEN

      DEALLOCATE (optp%optpiclb%iwn2013, STAT = ERR)
      THROW( ERR .NE. 0 )
      DEALLOCATE (optp%optpiclb%jwn2013, STAT = ERR)
      THROW( ERR .NE. 0 )
      DEALLOCATE (optp%optpiclb%dx_dwn2013, STAT = ERR)
      THROW( ERR .NE. 0 )
      NULLIFY (optp%optpiclb%iwn2013)
      NULLIFY (optp%optpiclb%jwn2013)
      NULLIFY (optp%optpiclb%dx_dwn2013)
      DEALLOCATE (optp%optpiclb%iwn2014, STAT = ERR)
      THROW( ERR .NE. 0 )
      DEALLOCATE (optp%optpiclb%jwn2014, STAT = ERR)
      THROW( ERR .NE. 0 )
      DEALLOCATE (optp%optpiclb%dx_dwn2014, STAT = ERR)
      THROW( ERR .NE. 0 )
      NULLIFY (optp%optpiclb%iwn2014)
      NULLIFY (optp%optpiclb%jwn2014)
      NULLIFY (optp%optpiclb%dx_dwn2014)
      DEALLOCATE (optp%optpiclb, STAT = ERR)
      THROW( ERR .NE. 0 )
      NULLIFY (optp%optpiclb)
    ENDIF


    CALL rttov_nullify_optpar_ir (optp)


  CATCH
END SUBROUTINE 
