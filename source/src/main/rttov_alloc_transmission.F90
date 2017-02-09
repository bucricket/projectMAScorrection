! Description:
!> @file
!!   Allocate/deallocate transmission structure.
!
!> @brief
!!   Allocate/deallocate transmission structure.
!!
!! @details
!!   The transmittance structure contains the calculated
!!   transmittances for the direct model and the output
!!   transmittance perturbations for the TL model.
!!   For the AD and K models the transmittance_ad/k
!!   should usually be initialised to zero.
!!
!! @param[out]    err            status on exit
!! @param[in,out] transmission   transmittances
!! @param[in]     nlayers        number of layers in input profiles (nlevels-1)
!! @param[in]     nchannels      total number of channels being simulated (SIZE(chanprof))
!! @param[in]     asw            1_jpim => allocate; 0_jpim => deallocate
!! @param[in]     init           set .TRUE. to initialise newly allocated structures, optional
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
SUBROUTINE rttov_alloc_transmission( &
            & err,          &
            & transmission, &
            & nlayers,      &
            & nchannels,    &
            & asw,          &
            & init)
!INTF_OFF
#include "throw.h"
!INTF_ON

  USE rttov_types, ONLY : transmission_type
  USE parkind1, ONLY : jpim, jplm

  IMPLICIT NONE

  INTEGER(KIND=jpim)     , INTENT(OUT)             :: err          ! return code
  TYPE(transmission_type), INTENT(INOUT)           :: transmission
  INTEGER(KIND=jpim)     , INTENT(IN)              :: nlayers
  INTEGER(KIND=jpim)     , INTENT(IN)              :: nchannels
  INTEGER(KIND=jpim)     , INTENT(IN)              :: asw          ! 1=allocate, 0=deallocate
  LOGICAL(KIND=jplm)     , INTENT(IN)   , OPTIONAL :: init
!INTF_END

#include "rttov_errorreport.interface"
#include "rttov_init_transmission.interface"

  INTEGER(KIND=jpim) :: nlevels
  LOGICAL(KIND=jplm) :: init1
!- End of header --------------------------------------------------------

  TRY
  nlevels = nlayers + 1
  init1   = .FALSE.
  IF (PRESENT(init)) init1 = init

  IF (asw .EQ. 1) THEN
    ALLOCATE (transmission%tau_levels(nlevels, nchannels), STAT = err)
    THROWM(err .NE. 0, "allocation of transmission%tau_levels")
    ALLOCATE (transmission%tau_total(nchannels), STAT = err)
    THROWM(err .NE. 0, "allocation of transmission%tau_total")
    ALLOCATE (transmission%tausun_levels_path2(nlevels, nchannels), STAT = err)
    THROWM(err .NE. 0, "allocation of transmission%tausun_levels_path2")
    ALLOCATE (transmission%tausun_total_path2(nchannels), STAT = err)
    THROWM(err .NE. 0, "allocation of transmission%tausun_total_path2")
    ALLOCATE (transmission%tausun_levels_path1(nlevels, nchannels), STAT = err)
    THROWM(err .NE. 0, "allocation of transmission%tausun_levels_path1")
    ALLOCATE (transmission%tausun_total_path1(nchannels), STAT = err)
    THROWM(err .NE. 0, "allocation of transmission%tausun_total_path1")

    IF (init1) THEN
      CALL rttov_init_transmission(transmission)
    ENDIF
  ENDIF

  IF (asw .EQ. 0) THEN
    DEALLOCATE (transmission%tau_levels, STAT = err)
    THROWM(err .NE. 0, "deallocation of transmission%tau_levels")
    DEALLOCATE (transmission%tau_total, STAT = err)
    THROWM(err .NE. 0, "deallocation of transmission%tau_total")
    DEALLOCATE (transmission%tausun_levels_path2, STAT = err)
    THROWM(err .NE. 0, "deallocation of transmission%tausun_levels_path2")
    DEALLOCATE (transmission%tausun_total_path2, STAT = err)
    THROWM(err .NE. 0, "deallocation of transmission%tausun_total_path2")
    DEALLOCATE (transmission%tausun_levels_path1, STAT = err)
    THROWM(err .NE. 0, "deallocation of transmission%tausun_levels_path1")
    DEALLOCATE (transmission%tausun_total_path1, STAT = err)
    THROWM(err .NE. 0, "deallocation of transmission%tausun_total_path1")
    NULLIFY (transmission%tau_levels)
    NULLIFY (transmission%tau_total)
    NULLIFY (transmission%tausun_levels_path2)
    NULLIFY (transmission%tausun_total_path2)
    NULLIFY (transmission%tausun_levels_path1)
    NULLIFY (transmission%tausun_total_path1)
  ENDIF
  CATCH
END SUBROUTINE rttov_alloc_transmission
