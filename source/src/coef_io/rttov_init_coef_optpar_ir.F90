!
SUBROUTINE rttov_init_coef_optpar_ir(  ERR, coef, optp )
! Description:
!
!   IR cloud coef arrays initialisation
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
!    Copyright 2012, EUMETSAT, All Rights Reserved.
!
! Method:
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0       Oct 2012    Original (P. Brunel J. Vidot
!  1.1       Apr 2014    Add improved baran scheme (J. Vidot)
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
! Declarations:
! Modules used:
! Imported Parameters:
#include "throw.h"
! Imported Type Definitions:
  USE rttov_types, ONLY :  &
       & rttov_coef ,&
       & rttov_optpar_ir
  USE parkind1, ONLY : jpim
!INTF_OFF 
  USE mod_rttov_baran2013_icldata, ONLY : &
       & baran2013_wvn, &
       & n_baran2013_wn
  USE mod_rttov_baran2014_icldata, ONLY : &
       & baran2014_wvn, &
       & n_baran2014_wn
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim)       , INTENT(OUT)      :: ERR
  TYPE(rttov_coef         ), INTENT(IN)       :: coef
  TYPE(rttov_optpar_ir),     INTENT(INOUT)    :: optp

!INTF_END
#include "rttov_errorreport.interface"
! Local Scalars:
  INTEGER(KIND=jpim) :: ichn, iwn, jwn
  REAL(KIND=JPRB)    :: dx_dwn
  
!- End of header --------------------------------------------------------
  TRY
  
    ALLOCATE (optp%optpiclb, STAT = ERR)
    THROWM( ERR .NE. 0, "allocation of optp%optpiclb" )

    ALLOCATE (optp%optpiclb%iwn2013(coef%fmv_chn), STAT = ERR)
    THROWM( ERR .NE. 0, "allocation of optp%optpiclb%iwn2013" )
    
    ALLOCATE (optp%optpiclb%jwn2013(coef%fmv_chn), STAT = ERR)
    THROWM( ERR .NE. 0, "allocation of optp%optpiclb%jwn2013" )
    
    ALLOCATE (optp%optpiclb%dx_dwn2013(coef%fmv_chn), STAT = ERR)
    THROWM( ERR .NE. 0, "allocation of optp%optpiclb%dx_dwn2013" )

    ALLOCATE (optp%optpiclb%iwn2014(coef%fmv_chn), STAT = ERR)
    THROWM( ERR .NE. 0, "allocation of optp%optpiclb%iwn2014" )
    
    ALLOCATE (optp%optpiclb%jwn2014(coef%fmv_chn), STAT = ERR)
    THROWM( ERR .NE. 0, "allocation of optp%optpiclb%jwn2014" )
    
    ALLOCATE (optp%optpiclb%dx_dwn2014(coef%fmv_chn), STAT = ERR)
    THROWM( ERR .NE. 0, "allocation of optp%optpiclb%dx_dwn2014" )

    DO ichn = 1, coef%fmv_chn
      IF (baran2013_wvn(1_JPIM) .GE. coef%ff_cwn(ichn)) THEN
        iwn = 1_JPIM
        jwn = 1_JPIM
        dx_dwn  = 0.0_JPRB
      ELSEIF (baran2013_wvn(n_baran2013_wn) .LE. coef%ff_cwn(ichn)) THEN
        iwn = n_baran2013_wn
        jwn = n_baran2013_wn
        dx_dwn  = 0.0_JPRB
      ELSE
        iwn = 1_JPIM
        DO WHILE (baran2013_wvn(iwn) .LE. coef%ff_cwn(ichn))
         iwn = iwn + 1_JPIM
        ENDDO
        iwn = iwn - 1_JPIM
        jwn = iwn + 1_JPIM
        dx_dwn  = (coef%ff_cwn(ichn) - baran2013_wvn(iwn)) / (baran2013_wvn(jwn)  - baran2013_wvn(iwn))
      ENDIF
      
      optp%optpiclb%iwn2013(ichn) = iwn
      optp%optpiclb%jwn2013(ichn) = jwn
      optp%optpiclb%dx_dwn2013(ichn) = dx_dwn

      IF (baran2014_wvn(1_JPIM) .GE. coef%ff_cwn(ichn)) THEN
        iwn = 1_JPIM
        jwn = 1_JPIM
        dx_dwn  = 0.0_JPRB
      ELSEIF (baran2014_wvn(n_baran2014_wn) .LE. coef%ff_cwn(ichn)) THEN
        iwn = n_baran2014_wn
        jwn = n_baran2014_wn
        dx_dwn  = 0.0_JPRB
      ELSE
        iwn = 1_JPIM
        DO WHILE (baran2014_wvn(iwn) .LE. coef%ff_cwn(ichn))
         iwn = iwn + 1_JPIM
        ENDDO
        iwn = iwn - 1_JPIM
        jwn = iwn + 1_JPIM
        dx_dwn  = (coef%ff_cwn(ichn) - baran2014_wvn(iwn)) / (baran2014_wvn(jwn)  - baran2014_wvn(iwn))
      ENDIF
      
      optp%optpiclb%iwn2014(ichn) = iwn
      optp%optpiclb%jwn2014(ichn) = jwn
      optp%optpiclb%dx_dwn2014(ichn) = dx_dwn
      
      ! Interpolations will be done like this:
      ! value = value(iwn) + ( value(jwn) - value(iwn) ) * dx_dwn

    ENDDO
  CATCH
END SUBROUTINE 
