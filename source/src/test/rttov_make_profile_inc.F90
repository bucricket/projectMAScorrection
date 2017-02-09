
SUBROUTINE rttov_make_profile_inc(profiles_inc, profiles, opts)
! Description:
!   Computes a profile increment
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
!    Copyright 2010, EUMETSAT, All Rights Reserved.
!
! Method:
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
  USE rttov_types, ONLY : profile_Type, rttov_options
!INTF_OFF
  USE parkind1, ONLY : jprb,jpim
!INTF_ON
  IMPLICIT NONE
  TYPE(profile_type),  INTENT(INOUT) :: profiles_inc(:)
  TYPE(profile_type),  INTENT(IN)    :: profiles    (:)
  TYPE(rttov_options), INTENT(IN)    :: opts
!INTF_END
  INTEGER(KIND=jpim) :: j, nprofiles
  
  nprofiles = size(profiles_inc)
  DO j = 1, nprofiles
! increments for atmospheric variables
    If( opts%interpolation%addinterp .And. opts%interpolation%lgradp )Then
      profiles_inc(j)%p   =  0.00001_jprb * (Modulo(j,3_jpim)+1) * profiles(j)%p     ! tl on pressure levels
    Else
      profiles_inc(j)%p = 0._jprb
    End If
    profiles_inc(j)%t = - 0.005_jprb * profiles(j)%t ! - 0.5%
    profiles_inc(j)%q = - 0.01_jprb * profiles(j)%q  ! - 1%
    IF (associated(profiles_inc(j)%o3) ) profiles_inc(j)%o3  =  - 0.01_jprb * profiles(j)%o3(:) ! - 1%
    IF (associated(profiles_inc(j)%co2)) profiles_inc(j)%co2 =  - 0.01_jprb * profiles(j)%co2(:)! - 1%
    IF (associated(profiles_inc(j)%co) ) profiles_inc(j)%co  =  - 0.01_jprb * profiles(j)%co(:) ! - 1%
    IF (associated(profiles_inc(j)%n2o)) profiles_inc(j)%n2o =  - 0.01_jprb * profiles(j)%n2o(:)! - 1%
    IF (associated(profiles_inc(j)%ch4)) profiles_inc(j)%ch4 =  - 0.01_jprb * profiles(j)%ch4(:)! - 1%

    IF (associated(profiles_inc(j)%aerosols)) THEN
      profiles_inc(j)%aerosols =  - 0.01_jprb * profiles(j)%aerosols
    ENDIF

    IF (associated(profiles_inc(j)%cfrac)) THEN
      ! The cfrac increments must not alter the relative magnitudes of cfrac on consecutive
      ! layers (eg if cfrac(i) < cfrac(i+1) then this must hold for perturbed cfracs also).
      ! It also means that identical cfracs on adjacent levels must have identical increments.
      ! Increments must also not perturb cfrac beyond the range [0,1].
      profiles_inc(j)%cfrac(:) = - 0.01_jprb * profiles(j)%cfrac(:)
    ENDIF

    IF (associated(profiles_inc(j)%cloud)) THEN
      profiles_inc(j)%cloud =  - 0.0001_jprb * profiles(j)%cloud
      IF (associated(profiles_inc(j)%icede)) THEN
        profiles_inc(j)%icede = 0.000001_jprb * profiles(j)%icede
      ENDIF
    ENDIF

    IF (associated(profiles_inc(j)%clw)) profiles_inc(j)%clw = 0.001_jprb * profiles(j)%clw  ! 0.1%
!
    profiles_inc(j)%s2m%t       =  - 0.005_jprb * profiles(j)%s2m%t       ! -0.5% on T
    profiles_inc(j)%s2m%q       =  - 0.01_jprb  * profiles(j)%s2m%q       ! -1% on 2m wv
    profiles_inc(j)%s2m%p       =  - 0.00001_jprb * profiles(j)%s2m%p     ! -0.001% on pressure
    profiles_inc(j)%s2m%u       =    0.01_jprb  * profiles(j)%s2m%u       ! 1% of wind components
!     IF (profiles_inc(j)%s2m%u == 0._jprb) profiles_inc(j)%s2m%u = 0.01_jprb
    profiles_inc(j)%s2m%v       =    0.01_jprb  * profiles(j)%s2m%v       ! 1% of wind components
!     IF (profiles_inc(j)%s2m%v == 0._jprb) profiles_inc(j)%s2m%v = 0.01_jprb
    profiles_inc(j)%s2m%o       =  - 0.01_jprb  * profiles(j)%s2m%o
    profiles_inc(j)%s2m%wfetc   =  - 0.01_jprb  * profiles(j)%s2m%wfetc

    profiles_inc(j)%skin%t             =  - 0.005_jprb * profiles(j)%skin%t
    profiles_inc(j)%skin%fastem        =  - 0.001_jprb * profiles(j)%skin%fastem
    profiles_inc(j)%skin%salinity      =  0.005_jprb * profiles(j)%skin%salinity
    profiles_inc(j)%skin%foam_fraction =  - 0.01_jprb * profiles(j)%skin%foam_fraction

    profiles_inc(j)%ctp         =  - 0.00001_jprb * profiles(j)%ctp
    IF (profiles(j)%cfraction < 0.5_jprb) THEN   ! Ensure non-zero increment even when cfraction is zero
      profiles_inc(j)%cfraction =    0.005_jprb * (profiles(j)%cfraction + 0.005_jprb)
    ELSE
      profiles_inc(j)%cfraction =  - 0.005_jprb * profiles(j)%cfraction
    ENDIF
    profiles_inc(j)%zenangle    = 0._jprb
    profiles_inc(j)%azangle     = 0._jprb
    profiles_inc(j)%sunzenangle = 0._jprb
    profiles_inc(j)%sunazangle  = 0._jprb
    profiles_inc(j)%elevation   = 0._jprb
    profiles_inc(j)%latitude    = 0._jprb
    profiles_inc(j)%longitude   = 0._jprb
    profiles_inc(j)%snow_frac   = 0._jprb
    profiles_inc(j)%soil_moisture = 0._jprb
    profiles_inc(j)%be          = 0._jprb
    profiles_inc(j)%cosbk       = 0._jprb
  ENDDO

END SUBROUTINE

