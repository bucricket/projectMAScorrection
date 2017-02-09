! Description:
!> @file
!!   Allocate/deallocate the user optical parameter structure.
!
!> @brief
!!   Allocate/deallocate the user optical parameter structure for
!!   aerosols or clouds.
!!
!! @details
!!   This structure is used to pass in profiles of aerosol or cloud
!!   optical properties for each channel being simulated. For IR
!!   channels these consist of absorption and scattering coefficients
!!   and the back-scattering parameter. Phase functions are only
!!   required for solar single-scattering simulations.
!!
!!   Currently this input is not active in the TL/AD/K models so
!!   there are no TL/AD/K opt_param arguments.
!!
!! @param[out]    err            status on exit
!! @param[in,out] opt_param      optical property structure
!! @param[in]     nchanprof      total number of channels being simulated
!! @param[in]     nlayers        number of layers in input profiles
!! @param[in]     nphangle       number of angles over which phase functions are
!!                                 specified (for non-solar simulations this can be 1)
!! @param[in]     asw            1_jpim => allocate; 0_jpim => deallocate
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
SUBROUTINE rttov_alloc_opt_param( &
            & err,       &
            & opt_param, &
            & nchanprof, &
            & nlayers,   &
            & nphangle,  &
            & asw)

!INTF_OFF
#include "throw.h"
!INTF_ON

  USE parkind1, ONLY : jpim
  USE rttov_types, ONLY : rttov_opt_param

  IMPLICIT NONE

  INTEGER(KIND=jpim),    INTENT(OUT)   :: err
  TYPE(rttov_opt_param), INTENT(INOUT) :: opt_param
  INTEGER(KIND=jpim),    INTENT(IN)    :: nchanprof
  INTEGER(KIND=jpim),    INTENT(IN)    :: nlayers
  INTEGER(KIND=jpim),    INTENT(IN)    :: nphangle
  INTEGER(KIND=jpim),    INTENT(IN)    :: asw
!INTF_END

#include "rttov_errorreport.interface"

  TRY

  IF (asw == 1_jpim) THEN
    NULLIFY (opt_param%abs,      &
             opt_param%sca,      &
             opt_param%bpr,      &
             opt_param%pha,      &
             opt_param%phangle,  &
             opt_param%iphangle, &
             opt_param%cosphangle)
    ALLOCATE (opt_param%abs(nchanprof, nlayers),           &
              opt_param%sca(nchanprof, nlayers),           &
              opt_param%bpr(nchanprof, nlayers),           &
              opt_param%pha(nchanprof, nlayers, nphangle), &
              opt_param%phangle(nphangle),                 &
              STAT = err) !iphangle and cosphangle are allocated in init routine
    THROWM(err.ne.0,"Allocation of opt_param failed")
  ENDIF

  IF (asw == 0_jpim) THEN
    DEALLOCATE (opt_param%abs,        &
                opt_param%sca,        &
                opt_param%bpr,        &
                opt_param%pha,        &
                opt_param%phangle,    &
                STAT = err)
    THROWM(err.ne.0,"Deallocation of opt_param failed")
    IF (ASSOCIATED(opt_param%iphangle))   DEALLOCATE(opt_param%iphangle, STAT = err)
    THROWM(err.ne.0,"Deallocation of opt_param%iphangle failed")
    IF (ASSOCIATED(opt_param%cosphangle)) DEALLOCATE(opt_param%cosphangle, STAT = err)
    THROWM(err.ne.0,"Deallocation of opt_param%cosphangle failed")
    NULLIFY (opt_param%abs,      &
             opt_param%sca,      &
             opt_param%bpr,      &
             opt_param%pha,      &
             opt_param%phangle,  &
             opt_param%iphangle, &
             opt_param%cosphangle)
  ENDIF

  CATCH
END SUBROUTINE rttov_alloc_opt_param
