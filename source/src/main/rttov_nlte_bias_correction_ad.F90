SUBROUTINE rttov_nlte_bias_correction_ad(coef,     &
                                         profiles,   profiles_ad, &
                                         angles,   &
                                         chanprof, &
                                                     rad_ad) 
                                           !IN         !INOUT
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 25 November 1998, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, KNMI and MeteoFrance.
!
!    Copyright 2012, EUMETSAT, All Rights Reserved.
!
! Current Code Owner: SAF NWP

  USE rttov_types, ONLY : rttov_coef, profile_type, geometry_type, &
                          rttov_chanprof, radiance_type
!INTF_OFF
  USE rttov_const, ONLY : deg2rad
  USE rttov_types, ONLY : rttov_nlte_coef
  USE parkind1, ONLY : jpim, jprb, jplm
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_coef),     INTENT(IN)    :: coef   ! nlte_coef structure 
  TYPE(profile_type),   INTENT(IN)    :: profiles(:) ! Atmospheric profiles
  TYPE(profile_type),   INTENT(INOUT) :: profiles_ad(:) 
  TYPE(geometry_type),  INTENT(IN)    :: angles(:)   ! geometry angles
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof(:)   ! geometry angles
  TYPE(radiance_type),  INTENT(INOUT) :: rad_ad
!INTF_END

  TYPE(rttov_nlte_coef), POINTER :: nlte_coef  ! nlte_coef structure 

  REAL(jprb)    :: p1, p2, p3, p
  REAL(jprb)    :: t_ad(coef%nlte_coef%ncoef,SIZE(profiles))
  REAL(jprb)    :: c(2), z(4,SIZE(profiles)) ! Bilinear interpolation coefs
  REAL(jprb)    :: d_rad(coef%nlte_coef%ncoef)
  INTEGER(jpim) :: bounds(2,2,SIZE(profiles))
  INTEGER(jpim) :: sol(2,SIZE(profiles)), sat(2,SIZE(profiles))
  INTEGER(jpim) :: isol, isat, i, j, k
  INTEGER(jpim) :: iprof, prof, chan
  LOGICAL(jplm) :: do_nlte(SIZE(profiles))

#define lb1 bounds(1,1,iprof)
#define ub1 bounds(2,1,iprof)
#define lb2 bounds(1,2,iprof)
#define ub2 bounds(2,2,iprof)

  nlte_coef => coef%nlte_coef
  do_nlte(:) = .TRUE.
  
  p1 = 0.005_jprb; p2 = 0.2_jprb; p3 = 52.0_jprb
  bounds = -99_jpim

!===============================================================================!
! Static code (same for FWD/TL/AD/K) - calc interpolation coefs per profile (z) !
!===============================================================================!
  DO iprof = 1, SIZE(profiles)
    ! Initialise to safe values - they will be changed if in range
    sol(1:2,iprof) = (/nlte_coef%nsol - 1_jpim, nlte_coef%nsol/)
    sat(1:2,iprof) = (/nlte_coef%nsat - 1_jpim, nlte_coef%nsat/)

    IF(profiles(iprof)%sunzenangle > nlte_coef%sol_zen_angle(nlte_coef%nsol) &
      .OR. profiles(iprof)%sunzenangle < 0._jprb) THEN ! sun below horizon do nothing
      do_nlte(iprof) = .FALSE.
      CYCLE
    ENDIF

    ! find bounds for solar zenith angle
    DO isol = 1, nlte_coef%nsol
      IF(nlte_coef%sol_zen_angle(isol) > profiles(iprof)%sunzenangle) THEN
        sol(2, iprof) = isol
        sol(1, iprof) = isol - 1_jpim ! safe because sunzenangle always gt 0 here
        EXIT
      ENDIF
    ENDDO

    ! find bounds for viewing angle
    DO isat = 1, nlte_coef%nsat
      IF(nlte_coef%sat_zen_angle(isat) > profiles(iprof)%zenangle) THEN
        sat(2, iprof) = isat
        sat(1, iprof) = isat - 1_jpim
        EXIT
      ENDIF
    ENDDO

! Simple linear interpolation just now  
    c(1) = (angles(iprof)%seczen - nlte_coef%sec_sat(sat(2,iprof))) / &
           (nlte_coef%sec_sat(sat(1,iprof)) - nlte_coef%sec_sat(sat(2,iprof)))

    c(2) = (COS(profiles(iprof)%sunzenangle * deg2rad) - &
           nlte_coef%cos_sol(sol(2,iprof))) / &
           (nlte_coef%cos_sol(sol(1,iprof)) - nlte_coef%cos_sol(sol(2,iprof)))

    z(1:4,iprof) = (/c(1) * c(2), (1 - c(1)) * c(2), &
                     c(1) * (1 - c(2)), (1 - c(1)) * (1 - c(2))/)

! find pressure bounds for finding average temperature
    p = 0._jprb; j = 1_jpim
    DO WHILE (p < p3)
      p = profiles(iprof)%p(j)
      
      IF( p < p1 ) THEN
        j = j + 1
        CYCLE
      ELSE IF(p >= p1 .AND. p < p2) THEN ! contribution to t1_avg
        IF(bounds(1,1,iprof) < 0) bounds(1,1,iprof) = j
        bounds(2,1,iprof) = j
      ELSE IF(p >= p2 .AND. p <= p3) THEN! contribution to t2_avg
        IF(bounds(1,2,iprof) < 0) bounds(1,2,iprof) = j
        bounds(2,2,iprof) = j
      ENDIF
      j = j + 1
    ENDDO
  ENDDO
! END static code

!=====================================================================!
!      Calculate radiance bias correction for required channels       !
!=====================================================================!
  t_ad = 0._jprb
  DO i = 1, SIZE(chanprof)
    prof = chanprof(i)%prof

    IF(do_nlte(prof)) THEN
      chan = chanprof(i)%chan
    
      IF(chan >= nlte_coef%start_chan .AND. chan <= nlte_coef%end_chan) THEN  
        ! corresponding array element for chan
        k = chan - nlte_coef%start_chan + 1 

        d_rad(:) = &
          z(1,prof) * nlte_coef%coef(:, sat(1,prof), sol(1,prof), k) + &
          z(2,prof) * nlte_coef%coef(:, sat(2,prof), sol(1,prof), k) + &
          z(3,prof) * nlte_coef%coef(:, sat(1,prof), sol(2,prof), k) + &
          z(4,prof) * nlte_coef%coef(:, sat(2,prof), sol(2,prof), k)

        ! only take into account contribution from rad_ad%total because 
        ! perturbation from rad_ad%clear = 0 by definition

        t_ad(1, prof) = t_ad(1, prof) + rad_ad%total(i) * d_rad(1)
        t_ad(2, prof) = t_ad(2, prof) + rad_ad%total(i) * d_rad(2)
        ! t_ad(3,:) later zero'd
!       t_ad(3, prof) = t_ad(3, prof) + rad_ad%total(i) * d_rad(3,prof)
      ENDIF
    ENDIF
  ENDDO

!================================!
! Calc t1 and t2 avg per profile !
!================================!
  DO iprof = 1, SIZE(profiles)
    profiles_ad(iprof)%t(lb1:ub1) = profiles_ad(iprof)%t(lb1:ub1) + t_ad(1,iprof) / (ub1 - lb1 + 1)
    profiles_ad(iprof)%t(lb2:ub2) = profiles_ad(iprof)%t(lb2:ub2) + t_ad(2,iprof) / (ub2 - lb2 + 1)
  ENDDO

END SUBROUTINE rttov_nlte_bias_correction_ad
