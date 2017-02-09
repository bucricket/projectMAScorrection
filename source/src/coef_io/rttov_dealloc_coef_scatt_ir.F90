!
SUBROUTINE rttov_dealloc_coef_scatt_ir (ERR, coef_scatt_ir)
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
       & rttov_coef_scatt_ir
  USE parkind1, ONLY : jpim
!INTF_OFF
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(out):
  INTEGER(KIND=jpim)       , INTENT(OUT)   :: ERR          ! return code
  TYPE(rttov_coef_scatt_ir), INTENT(INOUT) :: coef_scatt_ir
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_nullify_coef_scatt_ir.interface"

!- End of header --------------------------------------------------------
  TRY

    IF (Associated(coef_scatt_ir%fmv_aer_comp_name)) DEALLOCATE (coef_scatt_ir%fmv_aer_comp_name, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%fmv_wcl_comp_name)) DEALLOCATE (coef_scatt_ir%fmv_wcl_comp_name, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%fmv_icl_comp_name)) DEALLOCATE (coef_scatt_ir%fmv_icl_comp_name, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%fmv_aer_rh)) DEALLOCATE (coef_scatt_ir%fmv_aer_rh, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%fmv_wcl_rh)) DEALLOCATE (coef_scatt_ir%fmv_wcl_rh, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%fmv_aer_rh_val)) DEALLOCATE (coef_scatt_ir%fmv_aer_rh_val, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%fmv_wcl_rh_val)) DEALLOCATE (coef_scatt_ir%fmv_wcl_rh_val, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%fmv_wcl_ph_val_cos)) DEALLOCATE (coef_scatt_ir%fmv_wcl_ph_val_cos, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%ifmv_wcl_ph_val)) DEALLOCATE (coef_scatt_ir%ifmv_wcl_ph_val, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%fmv_aer_ph_val)) DEALLOCATE (coef_scatt_ir%fmv_aer_ph_val, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%fmv_aer_ph_val_cos)) DEALLOCATE (coef_scatt_ir%fmv_aer_ph_val_cos, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%ifmv_aer_ph_val)) DEALLOCATE (coef_scatt_ir%ifmv_aer_ph_val, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%fmv_wcl_ph_val)) DEALLOCATE (coef_scatt_ir%fmv_wcl_ph_val, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%fmv_icl_ph_val)) DEALLOCATE (coef_scatt_ir%fmv_icl_ph_val, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%fmv_icl_ph_val_cos)) DEALLOCATE (coef_scatt_ir%fmv_icl_ph_val_cos, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%ifmv_icl_ph_val)) DEALLOCATE (coef_scatt_ir%ifmv_icl_ph_val, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%fmv_icl_dg)) DEALLOCATE (coef_scatt_ir%fmv_icl_dg, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%aer_pha_chanlist)) DEALLOCATE (coef_scatt_ir%aer_pha_chanlist, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%wcl_pha_chanlist)) DEALLOCATE (coef_scatt_ir%wcl_pha_chanlist, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%icl_pha_chanlist)) DEALLOCATE (coef_scatt_ir%icl_pha_chanlist, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%aer_pha_index)) DEALLOCATE (coef_scatt_ir%aer_pha_index, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%wcl_pha_index)) DEALLOCATE (coef_scatt_ir%wcl_pha_index, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%icl_pha_index)) DEALLOCATE (coef_scatt_ir%icl_pha_index, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%abs)) DEALLOCATE (coef_scatt_ir%abs, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%sca)) DEALLOCATE (coef_scatt_ir%sca, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%bpr)) DEALLOCATE (coef_scatt_ir%bpr, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%pha)) DEALLOCATE (coef_scatt_ir%pha, STAT = ERR)

    THROW( ERR .NE. 0 )

    IF (Associated(coef_scatt_ir%confac)) DEALLOCATE (coef_scatt_ir%confac, STAT = ERR)

    THROW( ERR .NE. 0 )

    CALL rttov_nullify_coef_scatt_ir (coef_scatt_ir)

  CATCH
END SUBROUTINE
