MODULE rttov_wrapper_transfer
! Description
!   Helper subroutines for RTTOV wrapper: transfers data
!     between input/output arrays and RTTOV structures
!
!   rttov_copy_to_profiles - copy from arrays to profiles
!   rttov_copy_from_profiles_k - copy from profiles_k to arrays
!   rttov_copy_to_radiance_k - copy from arrays to radiance_k
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

#include "throw.h"

IMPLICIT NONE

CONTAINS

SUBROUTINE rttov_copy_to_profiles( &
    rth,                   &
    profiles,              &
    iprof1,                &
    datetimes,             &
    angles, surfgeom,      &
    surftype, skin, s2m,   &
    simplecloud, icecloud, &
    zeeman,                &
    p, t, gas_units, gas_id, gases)

! Description
!   Helper subroutine for RTTOV wrapper: copy profile data
!      from arrays into RTTOV profile structure

  USE parkind1, ONLY : jpim, jprb

  USE rttov_types, ONLY : &
      profile_type

  USE rttov_wrapper_handle

  USE rttov_const, ONLY : nwcl_max, naer_max

  IMPLICIT NONE

  TYPE(rttovwrapperhandle_type), INTENT(IN)    :: rth
  TYPE(profile_type),            INTENT(INOUT) :: profiles(:)
  INTEGER(jpim),                 INTENT(IN)    :: iprof1
  INTEGER(jpim),                 INTENT(IN)    :: datetimes(:,:)   !(6,nprofiles)
  REAL(jprb),                    INTENT(IN)    :: angles(:,:)      !(4,nprofiles)
  REAL(jprb),                    INTENT(IN)    :: surfgeom(:,:)    !(3,nprofiles)
  INTEGER(jpim),                 INTENT(IN)    :: surftype(:,:)    !(2,nprofiles)
  REAL(jprb),                    INTENT(IN)    :: skin(:,:)        !(9,nprofiles)
  REAL(jprb),                    INTENT(IN)    :: s2m(:,:)         !(6,nprofiles)
  REAL(jprb),                    INTENT(IN)    :: simplecloud(:,:) !(2,nprofiles)
  INTEGER(jpim),                 INTENT(IN)    :: icecloud(:,:)    !(2,nprofiles)
  REAL(jprb),                    INTENT(IN)    :: zeeman(:,:)      !(2,nprofiles)
  REAL(jprb),                    INTENT(IN)    :: p(:,:)           !(nlevels,nprofiles)
  REAL(jprb),                    INTENT(IN)    :: t(:,:)           !(nlevels,nprofiles)
  INTEGER(jpim),                 INTENT(IN)    :: gas_units
  INTEGER(jpim),                 INTENT(IN)    :: gas_id(:)        !(ngases)
  REAL(jprb),                    INTENT(IN)    :: gases(:,:,:)     !(nlevels,nprofiles,ngases)

  INTEGER(jpim) :: i, j, k, g, nprofiles, ngases
!------------------------------------------------------------------------------


  nprofiles = SIZE(profiles)
  ngases = SIZE(gas_id)

  DO i = 1, nprofiles
    j = iprof1 + i - 1

    profiles(i)%p(:) = p(:,j)
    profiles(i)%t(:) = t(:,j)

    profiles(i)%gas_units = gas_units

    DO g = 1, ngases
      IF (gas_id(g) == gas_id_q) THEN
        profiles(i)%q(:) = gases(:,j,g)
      ELSEIF (gas_id(g) == gas_id_o3 .AND. rth%opts%rt_ir%ozone_data) THEN
        profiles(i)%o3(:) = gases(:,j,g)
      ELSEIF (gas_id(g) == gas_id_co2 .AND. rth%opts%rt_ir%co2_data) THEN
        profiles(i)%co2(:) = gases(:,j,g)
      ELSEIF (gas_id(g) == gas_id_n2o .AND. rth%opts%rt_ir%n2o_data) THEN
        profiles(i)%n2o(:) = gases(:,j,g)
      ELSEIF (gas_id(g) == gas_id_co .AND. rth%opts%rt_ir%co_data) THEN
        profiles(i)%co(:) = gases(:,j,g)
      ELSEIF (gas_id(g) == gas_id_ch4 .AND. rth%opts%rt_ir%ch4_data) THEN
        profiles(i)%ch4(:) = gases(:,j,g)
      ELSEIF (gas_id(g) == gas_id_clw .AND. rth%opts%rt_mw%clw_data) THEN
        profiles(i)%clw(:) = gases(:,j,g)
      ENDIF

      IF (rth%opts%rt_ir%addclouds) THEN
        DO k = 1, nwcl_max
          IF (gas_id(g) == gas_id_lwc(k)) THEN
            profiles(i)%cloud(k,:) = gases(1:profiles(i)%nlayers,j,g)
          ENDIF
        ENDDO
        IF (gas_id(g) == gas_id_iwc) THEN
          profiles(i)%cloud(nwcl_max+1,:) = gases(1:profiles(i)%nlayers,j,g)
        ELSEIF (gas_id(g) == gas_id_icede) THEN
          profiles(i)%icede(:) = gases(1:profiles(i)%nlayers,j,g)
        ELSEIF (gas_id(g) == gas_id_cfrac) THEN
          profiles(i)%cfrac(:) = gases(1:profiles(i)%nlayers,j,g)
        ENDIF
      ENDIF

      IF (rth%opts%rt_ir%addaerosl) THEN
        DO k = 1, naer_max
          IF (gas_id(g) == gas_id_aer(k)) THEN
            profiles(i)%aerosols(k,:) = gases(1:profiles(i)%nlayers,j,g)
          ENDIF
        ENDDO
      ENDIF
    ENDDO

    profiles(i)%s2m%p              = s2m(1,j)
    profiles(i)%s2m%t              = s2m(2,j)
    profiles(i)%s2m%q              = s2m(3,j)
    profiles(i)%s2m%u              = s2m(4,j)
    profiles(i)%s2m%v              = s2m(5,j)
    profiles(i)%s2m%wfetc          = s2m(6,j)

    profiles(i)%skin%surftype      = surftype(1,j)
    profiles(i)%skin%watertype     = surftype(2,j)
    profiles(i)%skin%t             = skin(1,j)
    profiles(i)%skin%salinity      = skin(2,j)
    profiles(i)%snow_frac          = skin(3,j)
    profiles(i)%skin%foam_fraction = skin(4,j)
    profiles(i)%skin%fastem(:)     = skin(5:9,j)

    profiles(i)%zenangle           = angles(1,j)
    profiles(i)%azangle            = angles(2,j)
    profiles(i)%sunzenangle        = angles(3,j)
    profiles(i)%sunazangle         = angles(4,j)

    profiles(i)%latitude           = surfgeom(1,j)
    profiles(i)%longitude          = surfgeom(2,j)
    profiles(i)%elevation          = surfgeom(3,j)

    profiles(i)%date(:)            = datetimes(1:3,j)
    profiles(i)%time(:)            = datetimes(4:6,j)

    profiles(i)%ctp                = simplecloud(1,j)
    profiles(i)%cfraction          = simplecloud(2,j)

    profiles(i)%ish                = icecloud(1,j)
    profiles(i)%idg                = icecloud(2,j)

    profiles(i)%be                 = zeeman(1,j)
    profiles(i)%cosbk              = zeeman(2,j)

  ENDDO

END SUBROUTINE rttov_copy_to_profiles

SUBROUTINE rttov_copy_from_profiles_k( &
    rth,            &
    profiles_k,     &
    iprof1,         &
    skin_k, s2m_k,  &
    simplecloud_k,  &
    p_k, t_k,       &
    gas_id, gases_k)

! Description
!   Helper subroutine for RTTOV wrapper: copy profile_k structure
!      into wrapper output K arrays.

  USE parkind1, ONLY : jpim, jprb

  USE rttov_types, ONLY : &
      profile_type

  USE rttov_wrapper_handle

  USE rttov_const, ONLY : nwcl_max, naer_max

  IMPLICIT NONE

  TYPE(rttovwrapperhandle_type), INTENT(IN)    :: rth
  TYPE(profile_type),            INTENT(IN)    :: profiles_k(:)
  INTEGER(jpim),                 INTENT(IN)    :: iprof1
  REAL(jprb),                    INTENT(INOUT) :: skin_k(:,:,:)        !(9,nchannels,nprofiles)
  REAL(jprb),                    INTENT(INOUT) :: s2m_k(:,:,:)         !(6,nchannels,nprofiles)
  REAL(jprb),                    INTENT(INOUT) :: simplecloud_k(:,:,:) !(2,nchannels,nprofiles)
  REAL(jprb),                    INTENT(INOUT) :: p_k(:,:,:)           !(nlevels,nchannels,nprofiles)
  REAL(jprb),                    INTENT(INOUT) :: t_k(:,:,:)           !(nlevels,nchannels,nprofiles)
  INTEGER(jpim),                 INTENT(IN)    :: gas_id(:)            !(ngases)
  REAL(jprb),                    INTENT(INOUT) :: gases_k(:,:,:,:)     !(nlevels,nchannels,nprofiles,ngases)

  INTEGER(jpim) :: i, j, k, t, g, iprof, nchannels, nprofiles, ngases, nlayers
!------------------------------------------------------------------------------


  nchannels = SIZE(p_k, 2)
  nprofiles = SIZE(profiles_k) / nchannels
  ngases = SIZE(gas_id)
  nlayers = SIZE(p_k, 1) - 1

  DO i = 1, nprofiles
    iprof = iprof1 + i - 1

    DO j = 1, nchannels
      k = (i - 1) * nchannels + j

      p_k(:,j,iprof) = profiles_k(k)%p(:)
      t_k(:,j,iprof) = profiles_k(k)%t(:)

      DO g = 1, ngases
        IF (gas_id(g) == gas_id_q) THEN
          gases_k(:,j,iprof,g) = profiles_k(k)%q(:)
        ELSEIF (gas_id(g) == gas_id_o3 .AND. rth%opts%rt_ir%ozone_data) THEN
          gases_k(:,j,iprof,g) = profiles_k(k)%o3(:)
        ELSEIF (gas_id(g) == gas_id_co2 .AND. rth%opts%rt_ir%co2_data) THEN
          gases_k(:,j,iprof,g) = profiles_k(k)%co2(:)
        ELSEIF (gas_id(g) == gas_id_n2o .AND. rth%opts%rt_ir%n2o_data) THEN
          gases_k(:,j,iprof,g) = profiles_k(k)%n2o(:)
        ELSEIF (gas_id(g) == gas_id_co .AND. rth%opts%rt_ir%co_data) THEN
          gases_k(:,j,iprof,g) = profiles_k(k)%co(:)
        ELSEIF (gas_id(g) == gas_id_ch4 .AND. rth%opts%rt_ir%ch4_data) THEN
          gases_k(:,j,iprof,g) = profiles_k(k)%ch4(:)
        ELSEIF (gas_id(g) == gas_id_clw .AND. rth%opts%rt_mw%clw_data) THEN
          gases_k(:,j,iprof,g) = profiles_k(k)%clw(:)
        ENDIF

        IF (rth%opts%rt_ir%addclouds) THEN
          DO t = 1, nwcl_max
            IF (gas_id(g) == gas_id_lwc(t)) THEN
              gases_k(1:nlayers,j,iprof,g) = profiles_k(k)%cloud(t,:)
            ENDIF
          ENDDO
          IF (gas_id(g) == gas_id_iwc) THEN
            gases_k(1:nlayers,j,iprof,g) = profiles_k(k)%cloud(nwcl_max+1,:)
          ELSEIF (gas_id(g) == gas_id_icede) THEN
            gases_k(1:nlayers,j,iprof,g) = profiles_k(k)%icede(:)
          ELSEIF (gas_id(g) == gas_id_cfrac) THEN
            gases_k(1:nlayers,j,iprof,g) = profiles_k(k)%cfrac(:)
          ENDIF
        ENDIF

        IF (rth%opts%rt_ir%addaerosl) THEN
          DO t = 1, naer_max
            IF (gas_id(g) == gas_id_aer(t)) THEN
              gases_k(1:nlayers,j,iprof,g) = profiles_k(k)%aerosols(t,:)
            ENDIF
          ENDDO
        ENDIF
      ENDDO

      s2m_k(1,j,iprof) = profiles_k(k)%s2m%p
      s2m_k(2,j,iprof) = profiles_k(k)%s2m%t
      s2m_k(3,j,iprof) = profiles_k(k)%s2m%q
      s2m_k(4,j,iprof) = profiles_k(k)%s2m%u
      s2m_k(5,j,iprof) = profiles_k(k)%s2m%v
      s2m_k(6,j,iprof) = profiles_k(k)%s2m%wfetc

      skin_k(1,j,iprof) = profiles_k(k)%skin%t
      skin_k(2,j,iprof) = profiles_k(k)%skin%salinity
      skin_k(3,j,iprof) = 0._jprb  ! snow_frac is not active in TL/AD/K
      skin_k(4,j,iprof) = profiles_k(k)%skin%foam_fraction
      skin_k(5:9,j,iprof) = profiles_k(k)%skin%fastem(:)

      simplecloud_k(1,j,iprof) = profiles_k(k)%ctp
      simplecloud_k(2,j,iprof) = profiles_k(k)%cfraction

    ENDDO ! channels
  ENDDO ! profiles

END SUBROUTINE rttov_copy_from_profiles_k

SUBROUTINE rttov_copy_to_radiance_k( &
    radiance_k, &
    iprof1,     &
    nprofiles,  &
    bt_k, rad_k)

! Description
!   Helper subroutine for RTTOV wrapper: copy K perturbation data
!      from arrays into RTTOV radiance_k structure

  USE parkind1, ONLY : jpim, jprb

  USE rttov_types, ONLY : &
      radiance_type

  IMPLICIT NONE

  TYPE(radiance_type), INTENT(INOUT) :: radiance_k
  INTEGER(jpim),       INTENT(IN)    :: iprof1
  INTEGER(jpim),       INTENT(IN)    :: nprofiles
  REAL(jprb),          INTENT(IN)    :: bt_k(:,:)    !(nchannels,nprofiles)
  REAL(jprb),          INTENT(IN)    :: rad_k(:,:)   !(nchannels,nprofiles)

  INTEGER(jpim) :: i, j, lo, hi, nchannels
!------------------------------------------------------------------------------


  nchannels = SIZE(bt_k, 1)

  DO i = 1, nprofiles
    j = iprof1 + i - 1
    lo = (i - 1) * nchannels + 1
    hi = i * nchannels

    ! RTTOV K takes care of using the correct perturbation based on switchrad
    ! and channel type (IR/MW vs VIS/NIR) so copy both bt_k and rad_k into radiance_k

    radiance_k%bt(lo:hi) = bt_k(:,j)
    radiance_k%total(lo:hi) = rad_k(:,j)
  ENDDO

END SUBROUTINE rttov_copy_to_radiance_k

END MODULE rttov_wrapper_transfer
