! Description:
!> @file
!!   Initialise an RTTOV-SCATT cloud profiles structure.
!
!> @brief
!!   Initialise an RTTOV-SCATT cloud profiles structure.
!!
!! @details
!!   The argument is an RTTOV-SCATT cloud profiles array: all
!!   array members will be initialised as well as cfrac.
!!
!! @param[in,out]  cld_profiles   Array of cloud profiles structures
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
SUBROUTINE rttov_init_scatt_prof(cld_profiles)

  USE rttov_types, ONLY : profile_cloud_type
!INTF_OFF
  USE parkind1, ONLY : jprb, jpim
!INTF_ON

  IMPLICIT NONE

  TYPE(profile_cloud_type), INTENT(INOUT) :: cld_profiles(:)
!INTF_END

  INTEGER(KIND=jpim) :: iprof, nprofiles

  nprofiles = SIZE(cld_profiles)

  DO iprof = 1, nprofiles
    cld_profiles(iprof)%cfrac = 0._jprb
    cld_profiles(iprof)%ph    = 0._jprb
    cld_profiles(iprof)%cc    = 0._jprb
    cld_profiles(iprof)%clw   = 0._jprb
    cld_profiles(iprof)%rain  = 0._jprb

    IF (ASSOCIATED(cld_profiles(iprof)%totalice)) cld_profiles(iprof)%totalice = 0._jprb
    IF (ASSOCIATED(cld_profiles(iprof)%sp))       cld_profiles(iprof)%sp       = 0._jprb
    IF (ASSOCIATED(cld_profiles(iprof)%ciw))      cld_profiles(iprof)%ciw      = 0._jprb
  ENDDO
END SUBROUTINE
