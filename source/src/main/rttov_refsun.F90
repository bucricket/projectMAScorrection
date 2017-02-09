!     Compute the fraction of solar radiance that is reflected by a wind
!     roughened water surface.
SUBROUTINE rttov_refsun( &
            & profiles,   &
            & coef,       &
            & aux,        &
            & sunglint,   &
            & raytracing)
!     Description:
!     The water surface is regarded as a collection of small mirror-like facets
!     each randomly tilted with respect to the local horizon. A time passes, the
!     tilt of a facet varies under the influence of the wind. When the open
!     ocean reflects the solar disk, these fluctuating facets produce a dancing
!     pattern known as sun glint.
!     The fraction of solar radiance reflected by a wind roughened water surface
!     is computed assuming that the slope of the facets obey to a
!     two-dimensional Gaussian random process  whose spectrum is specified by
!     the Joint North Sea Wave Project (JONSWAP) wave spectral model.
!     The total variance of the slope of the facet is obtained from the
!     frequency spectrum of the surface wave and the inverse function of the
!     dispersive relation of the full-gravity-capillary wave.
!     The shadowing of the surface of the facets on the backsides of the waves
!     and deep in the throughs between waves is considered.
!     A coordinate system is used where the average water surface lies in the
!     X-Y plane. The coordinate system is right-handed and is located at the
!     reflection point. The Z axis points toward the zenith and the X axis
!     points in the direction formed by the projection of the reflected ray
!     on the average water surface.
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
!     Method:
!     K. Yoshimori, K. Itoh and Y. Ichioka:'Optical charateristics of a
!     wind-roughened water surface:a two dimensional theory'.Applied Optics,
!     Vol.34,No.27,20 September 1995, pp.6236-6247.
!
!     Owner:
!     EUMETSAT
!
!     History:
!     Version      Date        Comment
!     1            15/07/2003  Oroginal code. RTIASI-4. Marco Matricardi. ECMWF.
!     2            27/02/2009  Determination of profile level index for the
!                              surface view angles simplified (P. Rayer)
!     3            02/12/2009  Pathsat, Pathsun and related quantities are now
!                              layer arrays (Marco Matricardi)
!     4            05/07/2010  Remove addsolar flag from profiles structure (J Hocking)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!
!     Module used:
  USE rttov_types, ONLY :  &
       & profile_aux,     &
       & profile_Type,    &
       & raytracing_type, &
       & sunglint_type,   &
       & rttov_coef
!INTF_OFF
  USE rttov_const, ONLY : deg2rad, pi,           &
       gravity, &
       surftype_sea, &
       max_sol_zen, min_exponent, max_exp_exponent
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  IMPLICIT NONE
!     Subroutine arguments:
!       Array arguments with intent in:
! the inclusion of solar
! radiation
  TYPE(profile_type   ), INTENT(IN)    :: profiles(:)
  TYPE(profile_aux    ), INTENT(IN)    :: aux
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(sunglint_type  ), INTENT(INOUT) :: sunglint
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
!INTF_END
!     End of subroutine arguments
#include "rttov_erfcx.interface"
!       Local arrays
  REAL   (KIND=jprb) :: ff      (coef%ws_nomega)  ! Working space
  REAL   (KIND=jprb) :: q_shad_a(size(profiles))
  REAL   (KIND=jprb) :: q_shad_b(size(profiles))
!       Local scalars:
  INTEGER(KIND=jpim) :: j, jj, i, np, ilevsur, ilevwrt
  REAL   (KIND=jprb) :: rqs                       ! Effective distribution function
  REAL   (KIND=jprb) :: theta_fi                  ! First order shadowing factor for the reflected ray
  REAL   (KIND=jprb) :: theta_csi                 ! First order shadowing factor for the incident ray surface
  REAL   (KIND=jprb) :: omega_1                   ! Frequency of the surface wave
  REAL   (KIND=jprb) :: k
  REAL   (KIND=jprb) :: sigma_a
  REAL   (KIND=jprb) :: sigma_b
  REAL   (KIND=jprb) :: sigma2
  REAL   (KIND=jprb) :: bb, aa, hinc, s2, s4, sa, sb, fac
  INTEGER(KIND=jpim) :: nprofiles

  REAL   (KIND=jprb) :: max_exp_exponent_4rt
  REAL   (KIND=jprb) :: omega_1_array(coef%ws_nomega)
  REAL   (KIND=jprb) :: omega_1_4_r(coef%ws_nomega)

  REAL   (KIND=JPRB) :: ZHOOK_HANDLE

  REAL   (KIND=JPRB) :: beta_exponent

!-----End of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_REFSUN', 0_jpim, ZHOOK_HANDLE)

  nprofiles = SIZE(profiles)
  fac = SQRT(2_jpim * pi)
  k = 3.3_jprb
  sigma_a = 0.07_jprb
  sigma_b = 0.09_jprb

  max_exp_exponent_4rt = SQRT(SQRT(max_exp_exponent))

  DO I = 1, coef%ws_nomega
    omega_1_array(i) = (i - 1) * 0.1_jprb
  ENDDO
  omega_1_4_r(2:coef%ws_nomega) = 1._jprb / (omega_1_array(2:coef%ws_nomega) ** 4_jpim)

  sunglint%s(:)%glint = 0._jprb
  sunglint%s(:)%omega = 0._jprb
  DO j = 1, nprofiles
    IF (profiles(j)%sunzenangle >= 0.0 .AND. &
        profiles(j)%sunzenangle < max_sol_zen) THEN
      IF (profiles(j)%skin%surftype == surftype_sea) THEN

!-----------Compute the wind speed----------------------------------------------
        sunglint%s(j)%windsp = sqrt(profiles(j)%s2m%u ** 2_jpim + profiles(j)%s2m%v ** 2_jpim)
!---------Compute the angle between the wind direction and the U axis-------
        IF (sunglint%s(j)%windsp > 0._jprb) THEN
          sunglint%s(j)%wangl = MOD(ATAN2(profiles(j)%s2m%v, profiles(j)%s2m%u), 2._jprb * pi)
        ELSE
          sunglint%s(j)%wangl = 0._jprb
        ENDIF
!-----------Compute the difference between the sun azimuth angle and the-------
!           azimuth angle of the direction formed by the projection on the
!           mean water surface of the surface-to-sun direction
! Note defn of azangle in appendix V of user guide.
!         sunglint%s(j)%dazng = profiles(j)%sunazangle - profiles(j)%azangle
        sunglint%s(j)%dazng = MOD(profiles(j)%azangle - profiles(j)%sunazangle, 360._jprb)
!-----------Find the index of the pressure level that is closest to the surface-
! two cases ...
! case-1: surf lies above profiles(j)%nlevels
! case-2: surf lies below profiles(j)%nlevels
        ilevwrt              = aux%s(j)%nearestlev_surf
! see rttov_profaux ...
! for case-1 this ilevwrt is first lev below surf
! for case-2 this ilevwrt = profiles(j)%nlevels (and is first lev above surf)
!-----------------------------------------------------
! first redundant section
!! case-2 is identified by -ve aux(j) %pfraction_surf
!             if (aux(j) %pfraction_surf < 0._jprb) then
!               ilevwrt=ilevwrt+1
!             end if
!-----------------------------------------------------
! second redundant section
!!! !           if(ilevwrt.gt.100)then
! should be
!             if(ilevwrt.gt.profiles(j)%nlevels)then
!               ilevsur=ilevwrt-1
!             else
!               ilevsur=ilevwrt
!             ENDIF
!-----------------------------------------------------
! level index for the raytracing angles to be used for surface
        ilevsur              = ilevwrt
!-----------Compute the secant of the satellite viewing angle at the surface----
        sunglint%s(j)%zensat = acos(1._jprb / raytracing%pathsat(ilevsur - 1, j))
!-----------Compute the secant of the sun zenith angle at the surface-----------
        sunglint%s(j)%zensun = acos(1._jprb / raytracing%pathsun(ilevsur - 1, j))
        IF (sunglint%s(j)%windsp == 0._jprb) THEN
!-----------In absence of wind the water surface is a mirror-like surface-------
          IF (sunglint%s(j)%zensat == sunglint%s(j)%zensun .AND. sunglint%s(j)%dazng == 180.0_jprb) THEN
            sunglint%s(j)%glint = 1._jprb
          ELSE
            sunglint%s(j)%glint = 0._jprb
          ENDIF
        ELSE
!-----------Set up some parameters used to define the JONSWAP frequency -------
!           spectrum of the surface wave                                       |
!------------------------------------------------------------------------------
          sunglint%s(j)%x_u     = profiles(j)%s2m%wfetc * gravity / sunglint%s(j)%windsp ** 2_jpim
          sunglint%s(j)%alfa1   = 0.076_jprb / sunglint%s(j)%x_u ** 0.22_jprb
          sunglint%s(j)%omega_m =      &
            & 2._jprb * pi * (3.5_jprb * (gravity / sunglint%s(j)%windsp) * (1._jprb / sunglint%s(j)%x_u ** 0.33_jprb))
!-----------Compute the total variance of the slope of the facet---------------
!                                                                              |
!           To compute the total variance of the slope of the facet (or the    |
!           total variance of the displacement of the water surface) the       |
!           knowledge of the frequency spectrum of the surface wave and the    |
!           inverse function of the dispersive relation of the full-gravity-   |
!           capillary wave is needed.                                          |
!                                                                              |
!           The dispersive relation of the full-gravity-capillary wave is      |
!           written as:                                                        |
!                                                                              |
!           omega_1(k)=k*sqrt((gravity/k+gamma*k/ro)*tanh(h*k))                   |
!                                                                              |
!           k    = wave-number                                                 |
!           gamma= surface-tension constant                                    |
!           ro   = density of water                                            |
!           h    = water depth                                                 |
!                                                                              |
!           The inverse function of the dispersive relation of the             |
!           full-gravity-capillary wave has been pre-computed numerically      |
!           using the following parameters:                                    |
!                                                                              |
!           gamma = 0.00735  [N/m]                                             |
!           ro    = 999.1    [kg/m^3]                                          |
!           h     = 50       [m]                                               |
!                                                                              |
!           The frequency spectrum (PSI) of the surface wave is computed       |
!           using the JONSWAP wave-spectral model.                             |
!                                                                              |
!           The total variance of the slope is then evaluated by integrating   |
!           the quantity (psi*k_omega**2):                                     |
!                                                                              |
!                     -infinity                                                |
!                    -                                                         |
!           gamma_sq=- (k_omega(omega_1)**2*psi(omega_1))  d omega_1           |
!                    -                                                         |
!                   -zero                                                      |
!                                                                              |
!           The integral is computed applying Simpson's rule using 302         |
!           quadrature points.                                                 |
!                                                                              |
!------------------------------------------------------------------------------
!-----------Compute the frequency spectrum of the surface wave------------------
          sunglint%psi(1, j) = 0._jprb
          FF(1) = 0._jprb
          DO I = 2, coef%ws_nomega
            omega_1 = omega_1_array(i)

            IF (sunglint%s(j)%omega_m < max_exp_exponent_4rt * omega_1) THEN ! stop exp overflow

              IF (omega_1 <= sunglint%s(j)%omega_m) THEN
                sigma2 = sigma_a ** 2_jpim
              ELSE
                sigma2 = sigma_b ** 2_jpim
              ENDIF

              beta_exponent = MIN(&
                   (omega_1 - sunglint%s(j)%omega_m) ** 2_jpim / &
                   (2._jprb * sigma2 * sunglint%s(j)%omega_m ** 2_jpim), &
                   max_exp_exponent)
              sunglint%beta(i, j) = EXP(-beta_exponent)

              sunglint%psi(i, j) = (sunglint%s(j)%alfa1 * gravity ** 2_jpim / omega_1 ** 5_jpim) * &
                   (EXP( - (5._jprb / 4._jprb) * sunglint%s(j)%omega_m ** 4_jpim * omega_1_4_r(i)))

              ! stop beta underflow
              IF (sunglint%beta(i, j) > min_exponent) sunglint%psi(i, j) = sunglint%psi(i, j) * (k ** sunglint%beta(i, j))

            ELSE
              sunglint%psi(i, j) = 0._jprb
            ENDIF
            FF(i) = sunglint%psi(i, j) * coef%ws_k_omega(i) ** 2_jpim
          ENDDO
!-----------Compute the Simpson's integral--------------------------------------
          bb   = 301._jprb
          aa   = 0._jprb
          np   = 1501_jpim
          hinc = (bb - aa) / REAL(2 * np, jprb)
          sa   = ff(1)
          sb   = ff(coef%ws_nomega)
          s4   = 0._jprb
          DO i = 1, nP
            jj = 2 * i - 1
            s4 = s4 + ff(jj)
          ENDDO
          s2 = 0._jprb
          DO i = 1, nP
            jj = 2 * i
            s2 = s2 + ff(jj)
          ENDDO
          sunglint%s(j)%gamma_sq = hinc / 3.0_jprb * (sa + sb + 4.0_jprb * s4 + 2.0_jprb * s2)
!-------------------------------------------------------------------------------
!-----------Compute the rms of the slope of the facet along the X (gamma_O)----
!           and Y (gamma_P) axis.                                              |
!------------------------------------------------------------------------------
! Note defn of azangle in appendix V of user guide.
!           sunglint%s(j)%wangl    = sunglint%s(j)%wangl - profiles(j)%azangle * deg2rad
          sunglint%s(j)%wangl    = &
              MOD(sunglint%s(j)%wangl - MOD(90.0_jprb-profiles(j)%azangle, 360.0_jprb) * deg2rad, 2._jprb * pi)
          sunglint%s(j)%gamma_o = sqrt((2 + cos(2._jprb * sunglint%s(j)%wangl)) * sunglint%s(j)%gamma_sq / 4._jprb)
          sunglint%s(j)%gamma_p = sqrt((2 - cos(2._jprb * sunglint%s(j)%wangl)) * sunglint%s(j)%gamma_sq / 4._jprb)
!-----------Obtain angles csi and alfa------------------------------------------
          sunglint%s(j)%csi     = atan(tan(sunglint%s(j)%zensun) * cos(PI - sunglint%s(j)%dazng * deg2rad))
          sunglint%s(j)%alfa    = asin(sin(sunglint%s(j)%zensun) * sin(PI - sunglint%s(j)%dazng * deg2rad))
          sunglint%s(j)%omega   = abs(sunglint%s(j)%csi + sunglint%s(j)%zensat) / 2._jprb
!-----------Compute the value of the function that represents the shadowing  --
!           of the incident and reflected ray                                  |
!------------------------------------------------------------------------------
!-----------First order shadowing (the slope of the facet is negative)---------
          sunglint%s(j)%gammax  = tan((sunglint%s(j)%csi - sunglint%s(j)%zensat) / 2._jprb)
          IF (abs(sunglint%s(j)%zensat) > 0._jprb) THEN
            IF ((1._jprb / tan(abs(sunglint%s(j)%zensat)) - &
                sunglint%s(j)%gammax * sunglint%s(j)%zensat / abs(sunglint%s(j)%zensat)) >= 0._jprb) THEN
              theta_fi = 1._jprb
            ELSE
              theta_fi = 0._jprb
            ENDIF
          ELSE
            theta_fi = 1._jprb
          ENDIF
          IF (abs(sunglint%s(j)%csi) > 0._jprb) THEN
            IF ((1._jprb / tan(abs(sunglint%s(j)%csi)) + &
                sunglint%s(j)%gammax * sunglint%s(j)%csi / abs(sunglint%s(j)%csi)) >= 0._jprb) THEN
              theta_csi = 1._jprb
            ELSE
              theta_csi = 0._jprb
            ENDIF
          ELSE
            theta_csi = 1._jprb
          ENDIF
          IF (abs(sunglint%s(j)%csi) > PI / 2._jprb) theta_csi = 0._jprb
!-----------Second order shadowing (the facet cannot be seen)-------------------
          IF ((sunglint%s(j)%zensat * sunglint%s(j)%csi) <= 0._jprb) THEN
            IF (sunglint%s(j)%zensat == 0._jprb .AND. sunglint%s(j)%csi == 0._jprb) THEN
              sunglint%s(j)%Q_shad = 1._jprb
            ELSE IF (sunglint%s(j)%zensat >= abs(sunglint%s(j)%csi)) THEN
              sunglint%s(j)%a_shad   = (1._jprb / tan(sunglint%s(j)%zensat)) / sunglint%s(j)%gamma_O
              sunglint%s(j)%lambda_a =                                                                           &
                & (1._jprb / (fac * sunglint%s(j)%a_shad)) * exp( - sunglint%s(j)%a_shad ** 2_jpim / 2._jprb) -  &
                & 0.5_jprb * RTTOV_ERFCX(sunglint%s(j)%a_shad / sqrt(2._jprb))
              Q_shad_A(j)            = 1._jprb / (1._jprb + sunglint%s(j)%lambda_a)
              sunglint%s(j)%Q_shad   = Q_shad_A(j) * theta_fi
            ELSE IF (abs(sunglint%s(j)%csi) > sunglint%s(j)%ZENSAT) THEN
              sunglint%s(j)%b_shad   = (1._jprb / tan(abs(sunglint%s(j)%csi))) / sunglint%s(j)%gamma_O
              sunglint%s(j)%lambda_b =                                                                           &
                & (1._jprb / (fac * sunglint%s(j)%b_shad)) * exp( - sunglint%s(j)%b_shad ** 2_jpim / 2._jprb) -  &
                & 0.5_jprb * RTTOV_ERFCX(sunglint%s(j)%b_shad / sqrt(2._jprb))
              Q_shad_B(j)            = 1._jprb / (1._jprb + sunglint%s(j)%lambda_b)
              sunglint%s(j)%Q_shad   = Q_shad_B(j) * theta_csi
            ENDIF
          ELSE
            sunglint%s(j)%a_shad   = (1._jprb / tan(sunglint%s(j)%ZENSAT)) / sunglint%s(j)%gamma_O
            sunglint%s(j)%b_shad   = (1._jprb / tan(abs(sunglint%s(j)%csi))) / sunglint%s(j)%gamma_O
            sunglint%s(j)%lambda_a =                                                                           &
              & (1._jprb / (fac * sunglint%s(j)%a_shad)) * exp( - sunglint%s(j)%a_shad ** 2_jpim / 2._jprb) -  &
              & 0.5_jprb * RTTOV_ERFCX(sunglint%s(j)%a_shad / sqrt(2._jprb))
            sunglint%s(j)%lambda_b =                                                                           &
              & (1._jprb / (fac * sunglint%s(j)%b_shad)) * exp( - sunglint%s(j)%b_shad ** 2_jpim / 2._jprb) -  &
              & 0.5_jprb * RTTOV_ERFCX(sunglint%s(j)%b_shad / sqrt(2._jprb))
            sunglint%s(j)%Q_shad   = 1._jprb / (1._jprb + sunglint%s(j)%lambda_a + sunglint%s(j)%lambda_b)
            sunglint%s(j)%Q_shad   = sunglint%s(j)%Q_shad * theta_fi * theta_csi
          ENDIF
!-------------------------------------------------------------------------------
!-----------Compute the probability density that the slope of the facet at ----
!           a certain point is gammax when the incident ray and the ray from   |
!           the reflection of the incident ray on the local surface do not     |
!           intersect with any other surface.                                  |
!------------------------------------------------------------------------------
          sunglint%s(j)%c_shad      = abs(cos(sunglint%s(j)%csi) + cos(sunglint%s(j)%ZENSAT)) * sunglint%s(j)%gamma_P
          sunglint%s(j)%p_prime     = 1._jprb / (fac * sunglint%s(j)%c_shad) *      &
            & exp( - 0.5_jprb * (tan(sunglint%s(j)%alfa) / sunglint%s(j)%c_shad) ** 2_jpim)
          sunglint%s(j)%g_shad      =      &
            & (1._jprb - tan((sunglint%s(j)%csi - sunglint%s(j)%ZENSAT) / 2._jprb) * tan(sunglint%s(j)%ZENSAT))
          sunglint%s(j)%fac1        = 2._jprb * (cos(sunglint%s(j)%alfa)) ** 2_jpim *      &
            & (cos((sunglint%s(j)%csi - sunglint%s(j)%ZENSAT) / 2)) ** 2_jpim
          sunglint%s(j)%Pxy_gammaxy = (1._jprb / (fac * sunglint%s(j)%gamma_O)) *      &
            & exp( - sunglint%s(j)%gammax ** 2_jpim / (2._jprb * sunglint%s(j)%gamma_O ** 2_jpim))
          rqs = sunglint%s(j)%Q_shad * sunglint%s(j)%g_shad * sunglint%s(j)%p_prime * sunglint%s(j)%Pxy_gammaxy /      &
            & sunglint%s(j)%fac1
!-------------------------------------------------------------------------------
          sunglint%s(j)%glint       = rqs! *pi*epsilon**2_jpim
        ENDIF
      ENDIF
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_REFSUN', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE RTTOV_REFSUN
