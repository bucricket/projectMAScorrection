! Description:
!> @file
!!   Allocate/deallocate radiance structure.
!
!> @brief
!!   Allocate/deallocate radiance structure.
!!
!! @details
!!   The radiance structure contains the output radiances, BTs and
!!   reflectances for the direct model, the output radiance
!!   perturbations for the TL model and the input gradients and
!!   perturbations for the AD and K models.
!!
!! @param[out]    err            status on exit
!! @param[in]     nchannels      total number of channels being simulated (SIZE(chanprof))
!! @param[in,out] radiance       radiances and corresponding BTs and BRFs
!! @param[in]     nlayers        number of layers in input profiles (nlevels-1)
!! @param[in]     asw            1_jpim => allocate; 0_jpim => deallocate
!! @param[in,out] radiance2      secondary radiances calculated by direct model only, optional
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
SUBROUTINE rttov_alloc_rad( &
            & err,       &
            & nchannels, &
            & radiance,  &
            & nlayers,   &
            & asw,       &
            & radiance2, &
            & init)
!INTF_OFF
#include "throw.h"
!INTF_ON

  USE rttov_types, ONLY : radiance_type, radiance2_type
  USE parkind1, ONLY : jpim, jplm

  IMPLICIT NONE

  INTEGER(KIND=jpim)  , INTENT(OUT)             :: err       ! return code
  INTEGER(KIND=jpim)  , INTENT(IN)              :: nchannels
  INTEGER(KIND=jpim)  , INTENT(IN)              :: nlayers
  TYPE(radiance_type) , INTENT(INOUT)           :: radiance  ! primary radiances
  INTEGER(KIND=jpim)  , INTENT(IN)              :: asw       ! 1=allocate, 0=deallocate
  TYPE(radiance2_type), OPTIONAL, INTENT(INOUT) :: radiance2 ! secondary radiances
  LOGICAL(KIND=jplm)  , OPTIONAL, INTENT(IN)    :: init
!INTF_END

#include "rttov_errorreport.interface"
#include "rttov_init_rad.interface"

  LOGICAL(KIND=jplm) :: init1
!- End of header --------------------------------------------------------

  TRY
  init1   = .FALSE.
  IF (PRESENT(init)) init1 = init

  IF (asw .EQ. 1) THEN
    ALLOCATE (radiance%clear(nchannels), STAT = err)
    THROWM(err .NE. 0, "allocation of radiance%clear")
    ALLOCATE (radiance%total(nchannels), STAT = err)
    THROWM(err .NE. 0, "allocation of radiance%total")
    ALLOCATE (radiance%bt_clear(nchannels), STAT = err)
    THROWM(err .NE. 0, "allocation of radiance%bt_clear")
    ALLOCATE (radiance%bt(nchannels), STAT = err)
    THROWM(err .NE. 0, "allocation of radiance%bt")
    ALLOCATE (radiance%refl_clear(nchannels), STAT = err)
    THROWM(err .NE. 0, "allocation of radiance%refl_clear")
    ALLOCATE (radiance%refl(nchannels), STAT = err)
    THROWM(err .NE. 0, "allocation of radiance%refl")
    ALLOCATE (radiance%cloudy(nchannels), STAT = err)
    THROWM(err .NE. 0, "allocation of radiance%cloudy")
    ALLOCATE (radiance%overcast(nlayers, nchannels), STAT = err)
    THROWM(err .NE. 0, "allocation of radiance%overcast")

    IF (PRESENT(radiance2)) THEN
      ALLOCATE (radiance2%upclear(nchannels), STAT = err)
      THROWM(err .NE. 0, "allocation of radiance2%upclear")
      ALLOCATE (radiance2%dnclear(nchannels), STAT = err)
      THROWM(err .NE. 0, "allocation of radiance2%dnclear")
      ALLOCATE (radiance2%refldnclear(nchannels), STAT = err)
      THROWM(err .NE. 0, "allocation of radiance2%refldnclear")
      ALLOCATE (radiance2%up(nlayers, nchannels), STAT = err)
      THROWM(err .NE. 0, "allocation of radiance2%up")
      ALLOCATE (radiance2%down(nlayers, nchannels), STAT = err)
      THROWM(err .NE. 0, "allocation of radiance2%down")
      ALLOCATE (radiance2%surf(nlayers, nchannels), STAT = err)
      THROWM(err .NE. 0, "allocation of radiance2%surf")
    ENDIF

    IF (init1) THEN
      CALL rttov_init_rad(radiance, radiance2)
    ENDIF
  ELSE

    DEALLOCATE (radiance%clear, STAT = err)
    THROWM(err .NE. 0, "deallocation of radiance%clear")
    DEALLOCATE (radiance%total, STAT = err)
    THROWM(err .NE. 0, "deallocation of radiance%total")
    DEALLOCATE (radiance%bt_clear, STAT = err)
    THROWM(err .NE. 0, "deallocation of radiance%bt_clear")
    DEALLOCATE (radiance%bt, STAT = err)
    THROWM(err .NE. 0, "deallocation of radiance%bt")
    DEALLOCATE (radiance%refl_clear, STAT = err)
    THROWM(err .NE. 0, "deallocation of radiance%refl_clear")
    DEALLOCATE (radiance%refl, STAT = err)
    THROWM(err .NE. 0, "deallocation of radiance%refl")
    DEALLOCATE (radiance%cloudy, STAT = err)
    THROWM(err .NE. 0, "deallocation of radiance%cloudy")
    DEALLOCATE (radiance%overcast, STAT = err)
    THROWM(err .NE. 0, "deallocation of radiance%overcast")

    IF (PRESENT(radiance2)) THEN
      DEALLOCATE (radiance2%upclear, STAT = err)
      THROWM(err .NE. 0, "deallocation of radiance2%upclear")
      DEALLOCATE (radiance2%dnclear, STAT = err)
      THROWM(err .NE. 0, "deallocation of radiance2%dnclear")
      DEALLOCATE (radiance2%refldnclear, STAT = err)
      THROWM(err .NE. 0, "deallocation of radiance2%refldnclear")
      DEALLOCATE (radiance2%up, STAT = err)
      THROWM(err .NE. 0, "deallocation of radiance2%up")
      DEALLOCATE (radiance2%down, STAT = err)
      THROWM(err .NE. 0, "deallocation of radiance2%down")
      DEALLOCATE (radiance2%surf, STAT = err)
      THROWM(err .NE. 0, "deallocation of radiance2%surf")
    ENDIF

    NULLIFY (radiance%clear)
    NULLIFY (radiance%total)
    NULLIFY (radiance%bt_clear)
    NULLIFY (radiance%bt)
    NULLIFY (radiance%refl_clear)
    NULLIFY (radiance%refl)
    NULLIFY (radiance%cloudy)
    NULLIFY (radiance%overcast)

    IF (PRESENT(radiance2)) THEN
      NULLIFY (radiance2%upclear)
      NULLIFY (radiance2%dnclear)
      NULLIFY (radiance2%refldnclear)
      NULLIFY (radiance2%up)
      NULLIFY (radiance2%down)
      NULLIFY (radiance2%surf)
    ENDIF
  ENDIF
  CATCH
END SUBROUTINE rttov_alloc_rad
