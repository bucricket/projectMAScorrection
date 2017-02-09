! Description:
!> @file
!!   Deallocate Mie coefficients.
!
!> @brief
!!   Deallocate Mie coefficients.
!!
!! @details
!!   Deallocates the Mie table structure passed as an argument.
!!
!! @param[in,out]  coef_scatt   Mie table structure to deallocate
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
Subroutine rttov_dealloc_scattcoeffs  (&
      & coef_scatt)     ! inout  
  ! History:
  ! version   date        comment
  ! -------   ----        -------
  !   1.0    04/2010      Initial version  (A Geer)

  ! Imported Type Definitions:
#include "throw.h"
  Use rttov_types, Only : rttov_scatt_coef
!INTF_OFF
  USE parkind1, ONLY : jpim, jprb
  USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!INTF_ON
  Implicit None

  ! subroutine arguments
  Type( rttov_scatt_coef ), Intent (inout) :: coef_scatt ! coefficients

!INTF_END

#include "rttov_errorreport.interface"

  ! local variables
  Integer(Kind=jpim)    :: ERR
  REAL(KIND=JPRB) :: ZHOOK_HANDLE

  !- End of header --------------------------------------------------------
  TRY

  IF (LHOOK) CALL DR_HOOK('RTTOV_DEALLOC_SCATTCOEFFS',0_jpim,ZHOOK_HANDLE)

  IF (Associated(coef_scatt % mie_freq)) DEALLOCATE (coef_scatt % mie_freq, STAT = ERR)
  THROW( ERR .NE. 0 )
  IF (Associated(coef_scatt % ext)) DEALLOCATE (coef_scatt % ext, STAT = ERR)
  THROW( ERR .NE. 0 )
  IF (Associated(coef_scatt % ssa)) DEALLOCATE (coef_scatt % ssa, STAT = ERR)
  THROW( ERR .NE. 0 )
  IF (Associated(coef_scatt % asp)) DEALLOCATE (coef_scatt % asp, STAT = ERR)
  THROW( ERR .NE. 0 )

  CATCH

  IF (LHOOK) CALL DR_HOOK('RTTOV_DEALLOC_SCATTCOEFFS',1_jpim,ZHOOK_HANDLE)
  End Subroutine rttov_dealloc_scattcoeffs
