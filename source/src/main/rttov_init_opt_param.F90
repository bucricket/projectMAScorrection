! Description:
!> @file
!!   Initialise auxiliary phase function arrays in the user
!!   optical parameter structure.
!
!> @brief
!!   Initialise auxiliary phase function arrays in the user
!!   aerosol/cloud optical parameter structure.
!!
!! @details
!!   This subroutine should be called after allocating an
!!   rttov_opt_param structure for solar scattering simulations.
!!
!! @param[out]    err            status on exit
!! @param[in]     opts           RTTOV options
!! @param[in,out] opt_param      optical property structure
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
SUBROUTINE rttov_init_opt_param(err, opts, opt_param)

!INTF_OFF
#include "throw.h"
!INTF_ON

  USE rttov_types, ONLY : rttov_options, rttov_opt_param
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE parkind1, ONLY : jprb
  USE rttov_const, ONLY : deg2rad, errorstatus_fatal
!INTF_ON

  IMPLICIT NONE

  INTEGER(KIND=jpim),    INTENT(INOUT)          :: err
  TYPE(rttov_options),   INTENT(IN)             :: opts
  TYPE(rttov_opt_param), INTENT(INOUT)          :: opt_param

!INTF_END

  INTEGER(KIND=jpim)  :: nphangle, i, k, icount

#include "rttov_errorreport.interface"

  TRY

  IF (opts%rt_ir%addsolar) THEN

    IF (.NOT. ASSOCIATED(opt_param%phangle)) THEN
      err = errorstatus_fatal
      THROWM(err .NE. 0, "opt_param%phangle must be allocated: call rttov_alloc_opt_param first")
    ENDIF

    nphangle = size(opt_param%phangle)
    IF (ASSOCIATED(opt_param%cosphangle)) DEALLOCATE(opt_param%cosphangle, STAT = err)
    THROWM(err .NE. 0, "deallocation of aux phase array cosphangle")
    ALLOCATE(opt_param%cosphangle(nphangle), STAT = err)
    THROWM(err .NE. 0, "allocation of aux phase array cosphangle")

    DO i = 1, nphangle
      opt_param%cosphangle(i) = COS(opt_param%phangle(i) * deg2rad)
    ENDDO

    opt_param%minphadiff = opt_param%phangle(2) - opt_param%phangle(1)
    DO i = 2, nphangle - 1
      IF (opt_param%minphadiff > opt_param%phangle(i + 1) - opt_param%phangle(i)) THEN
        opt_param%minphadiff = opt_param%phangle(i + 1) - opt_param%phangle(i)
      ENDIF
    ENDDO

    opt_param%minphadiff = opt_param%minphadiff / 2.0_JPRB
    icount = opt_param%phangle(nphangle) / opt_param%minphadiff
    IF (ASSOCIATED(opt_param%iphangle)) DEALLOCATE(opt_param%iphangle, STAT = err)
    THROWM(err .NE. 0, "deallocation of aux phase array iphangle")
    ALLOCATE (opt_param%iphangle(icount), STAT = err)
    THROWM(err .NE. 0, "allocation of aux phase array iphangle")

    k = 1
    DO i = 1, icount
      IF (opt_param%phangle(k) >= (i + 1) * opt_param%minphadiff) THEN
        opt_param%iphangle(i) = k
      ELSE
        k = k + 1
        opt_param%iphangle(i) = k
      ENDIF
    ENDDO

  ENDIF

  CATCH
END SUBROUTINE rttov_init_opt_param
