!
SUBROUTINE rttov_init_coef_pccomp( &
            & err,           &
            & coef,          &
            & coef_pccomp)
! Description:
!
!   PC coef array initialisation
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
!  1.0       17/05/2004  Original (based on rttov_readcoeffs)
!  1.1       01/06/2005  Marco Matricardi (ECMWF):
!               --       Added variables for CO2,N2O, CO and CH4.
!  1.2       2/11/2007   Made mwcldtop variable (R Saunders)
!  1.3       2/11/2007   Included polarimetric type for GHz calc (RS)
!  1.4       27/02/2009  Profile levels to include ToA. Distinguish between
!                        layer arrays and level arrays - size, index labels,
!                        looping. Prepares predictor variables to agree
!                        with RTTOV-9 (P. Rayer)
!  1.5       16/7/2009   Remove trace gases for RTTOV-8 and 9 which are not supported (RS)
!  1.6       02/12/2009  Introduced principal component capability. Marco Matricardi. ECMWF
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
#include "throw.h"
! Imported Type Definitions:
  USE rttov_types, ONLY :  &
       & rttov_coef,          &
       & rttov_coef_pccomp
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(in):
! scalar arguments with intent(out):
  INTEGER(KIND=jpim),      INTENT(OUT)   :: err
  TYPE(rttov_coef),        INTENT(INOUT) :: coef
  TYPE(rttov_coef_pccomp), INTENT(INOUT) :: coef_pccomp
!INTF_END

  INTEGER(KIND=jpim)  :: i,j

#include "rttov_errorreport.interface"

!- End of header --------------------------------------------------------
  TRY

  ALLOCATE (coef_pccomp%planck1_in(coef_pccomp%fmv_pc_nchn), STAT = ERR)
  THROWM( ERR .NE. 0, "allocation of Planck1_in" )

  ALLOCATE (coef_pccomp%planck2_in(coef_pccomp%fmv_pc_nchn), STAT = ERR)
  THROWM( ERR .NE. 0, "allocation of Planck2_in" )

  coef_pccomp%planck1_in(:) = coef%fc_planck_c1 * coef_pccomp%ff_cwn_in(:) ** 3
  coef_pccomp%planck2_in(:) = coef%fc_planck_c2 * coef_pccomp%ff_cwn_in(:)

  ALLOCATE (coef_pccomp%co2_pc_ref(coef%nlevels), STAT = ERR)
  THROWM( ERR .NE. 0, "allocation of co2_pc_ref" )

  ALLOCATE (coef_pccomp%n2o_pc_ref(coef%nlevels), STAT = ERR)
  THROWM( ERR .NE. 0, "allocation of n2o_pc_ref" )

  ALLOCATE (coef_pccomp%co_pc_ref(coef%nlevels), STAT = ERR)
  THROWM( ERR .NE. 0, "allocation of co_pc_ref" )

  ALLOCATE (coef_pccomp%ch4_pc_ref(coef%nlevels), STAT = ERR)
  THROWM( ERR .NE. 0, "allocation of ch4_pc_ref" )

  IF (coef%nco2 > 0) THEN
    coef_pccomp%co2_pc_ref(:) = coef_pccomp%ref_pc_prfl_mr(:, 1)
  ENDIF

  IF (coef%nn2o > 0) THEN
    coef_pccomp%n2o_pc_ref(:) = coef_pccomp%ref_pc_prfl_mr(:, 2)
  ENDIF

  IF (coef%nco > 0) THEN
    coef_pccomp%co_pc_ref(:) = coef_pccomp%ref_pc_prfl_mr(:, 3)
  ENDIF

  IF (coef%nch4 > 0) THEN
    coef_pccomp%ch4_pc_ref(:) = coef_pccomp%ref_pc_prfl_mr(:, 4)
  ENDIF

! DAR: Added for RTTOV 11.2 - Allocate memory for transposed arrays for PC calcs
! We're doing it this way in order to maintain compatibility with the current
! coef files. The additional memory requirements for IASI are approx 40MB.
  DO j = 1, SIZE(coef_pccomp%pcreg(:,:), dim=2)
    DO i = 1, SIZE(coef_pccomp%pcreg(:,:), dim=1)
      ALLOCATE(coef_pccomp%pcreg(i,j)%coefficients_t( &
        SIZE(coef_pccomp%pcreg(i,j)%coefficients(1,:)), &
        SIZE(coef_pccomp%pcreg(i,j)%coefficients(:,1))))

      coef_pccomp%pcreg(i,j)%coefficients_t = &
        TRANSPOSE(coef_pccomp%pcreg(i,j)%coefficients)
    ENDDO
  ENDDO

  DO i = 1, SIZE(coef_pccomp%eigen(:))
    ALLOCATE(coef_pccomp%eigen(i)%eigenvectors_t( &
      SIZE(coef_pccomp%eigen(i)%eigenvectors(1,:)), &
      SIZE(coef_pccomp%eigen(i)%eigenvectors(:,1))))

    coef_pccomp%eigen(i)%eigenvectors_t = &
      TRANSPOSE(coef_pccomp%eigen(i)%eigenvectors)
  ENDDO

  ALLOCATE(coef_pccomp%noise_r(SIZE(coef_pccomp%noise)))
  coef_pccomp%noise_r = 1._jprb / coef_pccomp%noise

  CATCH
END SUBROUTINE 
