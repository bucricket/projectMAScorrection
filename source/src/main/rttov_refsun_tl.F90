SUBROUTINE rttov_refsun_tl( &
            & profiles,      &
            & profiles_tl,   &
            & coef,          &
            & aux,           &
            & sunglint,      &
            & sunglint_tl,   &
            & raytracing,    &
            & raytracing_tl)
!     Description:
!     Compute the fraction of solar radiance that is reflected by a wind
!     roughened water surface.
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
!     dispersive relation of the full-gravityity-capillary wave.
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
!     2            15/07/2009  Determination of profile level index for the
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
  USE rttov_const, ONLY :  &
       deg2rad, pi,           &
       gravity,      &
       surftype_sea, &
       max_sol_zen, max_exp_exponent, min_exponent

  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  IMPLICIT NONE
!     Subroutine arguments:
!       Array arguments with intent in:
! the inclusion of solar
! radiation
  TYPE(profile_type   ), INTENT(IN)    :: profiles   (:)
  TYPE(profile_type   ), INTENT(IN)    :: profiles_tl(size(profiles))
  TYPE(profile_aux    ), INTENT(IN)    :: aux
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(IN)    :: raytracing_tl
  TYPE(sunglint_type  ), INTENT(IN)    :: sunglint
  TYPE(sunglint_type  ), INTENT(INOUT) :: sunglint_tl
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
!INTF_END
!     end of subroutine arguments
!       Local arrays
  REAL   (KIND=jprb) :: psi_tl(coef%ws_nomega)              ! The frequency spectrum of the surface wave
  REAL   (KIND=jprb) :: ff_tl (coef%ws_nomega)              ! Working space
!       Local scalars:
  INTEGER(KIND=jpim) :: j, jj, i, np, ilevsur, ilevwrt
  REAL   (KIND=jprb) :: csi_tl                              ! Angle between the zenith and the projection of the
! incident ray on the x-z plane
  REAL   (KIND=jprb) :: alfa_tl                             ! Angle between the incident ray and the x-z plane
  REAL   (KIND=jprb) :: c_shad_tl                           ! The average magnitude of tan(alfa)
  REAL   (KIND=jprb) :: p_prime_tl                          ! The probability density of tan(alfa)
  REAL   (KIND=jprb) :: pxy_gammaxy_tl                      ! The joint  probability density of the along-view and
! cross view slope
  REAL   (KIND=jprb) :: rqs_tl                              ! Effective distribution function
  REAL   (KIND=jprb) :: gamma_sq_tl                         ! Total variance of the slope of the facet
  REAL   (KIND=jprb) :: gamma_o_tl                          ! The mean square of the along-view (X axis) slope
  REAL   (KIND=jprb) :: gamma_p_tl                          ! The mean square of the cross-view (Y axis) slope
  REAL   (KIND=jprb) :: g_shad_tl                           ! Normalization function
  REAL   (KIND=jprb) :: gammax_tl                           ! The x-slope of the facet
  REAL   (KIND=jprb) :: theta_fi                            ! First order shadowing factor for the reflected ray
  REAL   (KIND=jprb) :: theta_csi                           ! First ordesr shadowing factor for the incident ray
  REAL   (KIND=jprb) :: q_shad_tl                           ! Second order shadowing factor
  REAL   (KIND=jprb) :: q_shad_a_tl                         ! Second order shadowing factor
  REAL   (KIND=jprb) :: q_shad_b_tl                         ! Second order shadowing factor
  REAL   (KIND=jprb) :: zensat_tl                           ! Zenith angle of satellite viewing angle at surface
  REAL   (KIND=jprb) :: zensun_tl                           ! Zenith angle of sun at surface
  REAL   (KIND=jprb) :: omega_1                             ! Frequency of the surface wave
  REAL   (KIND=jprb) :: fac1_tl
  REAL   (KIND=jprb) :: a_shad_tl
  REAL   (KIND=jprb) :: b_shad_tl
  REAL   (KIND=jprb) :: lambda_a_tl
  REAL   (KIND=jprb) :: lambda_b_tl
  REAL   (KIND=jprb) :: x_u_tl
  REAL   (KIND=jprb) :: alfa1_tl
  REAL   (KIND=jprb) :: k
  REAL   (KIND=jprb) :: sigma_a
  REAL   (KIND=jprb) :: sigma_b
  REAL   (KIND=jprb) :: omega_m_tl
  REAL   (KIND=jprb) :: sigma2
  REAL   (KIND=jprb) :: beta_tl
  REAL   (KIND=jprb) :: windsp_tl
  REAL   (KIND=jprb) :: wangl_tl
  REAL   (KIND=jprb) :: s2_tl, s4_tl, sa_tl, sb_tl
  REAL   (KIND=jprb) :: bb, aa, hinc, fac
  INTEGER(KIND=jpim) :: nprofiles

  REAL   (KIND=jprb) :: max_exp_exponent_4rt, alfa1_r, logk
  REAL   (KIND=jprb) :: exp_max_exp_exponent_r, omega_m_3
  REAL   (KIND=jprb) :: omega_1_array(coef%ws_nomega)
  REAL   (KIND=jprb) :: omega_1_4_r(coef%ws_nomega)

  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!-----end of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_REFSUN_TL', 0_jpim, ZHOOK_HANDLE)
  nprofiles = SIZE(profiles)
  fac = SQRT(2_jpim * pi)
  k = 3.3_jprb
  sigma_a = 0.07_jprb
  sigma_b = 0.09_jprb

  max_exp_exponent_4rt = SQRT(SQRT(max_exp_exponent))
  logk = LOG(k) ! log(3.3)
  exp_max_exp_exponent_r = EXP(-max_exp_exponent)

  DO I = 1, coef%ws_nomega
    omega_1_array(i) = (i - 1) * 0.1_jprb
  ENDDO
  omega_1_4_r(2:coef%ws_nomega) = 1._jprb / (omega_1_array(2:coef%ws_nomega) ** 4_jpim)

  sunglint_tl%s(:)%glint = 0._jprb
  sunglint_tl%s(:)%omega = 0._jprb
  DO J = 1, nprofiles
    IF (profiles(j)%sunzenangle >= 0.0 .AND. &
        profiles(j)%sunzenangle < max_sol_zen) THEN
      IF (profiles(j)%skin%surftype == surftype_sea) THEN
!-----------Compute the wind speed----------------------------------------------
        IF (sunglint%s(j)%windsp > 0._jprb) THEN
          windsp_tl =      &
            & (profiles(j)%s2m%u * profiles_tl(j)%s2m%u + profiles(j)%s2m%v * profiles_tl(j)%s2m%v) / sunglint%s(j)%windsp
        ELSE
          windsp_tl = (profiles_tl(j)%s2m%u + profiles_tl(j)%s2m%v) / SQRT(2._jprb)
        ENDIF
!-----------Compute the angle between the wind direction and the U axis-------
        IF (sunglint%s(j)%windsp > 0._jprb) THEN
          wangl_tl = (1._jprb / (profiles(j)%s2m%u ** 2_jpim + profiles(j)%s2m%v ** 2_jpim)) *      &
            & (profiles(j)%s2m%u * profiles_tl(j)%s2m%v - profiles(j)%s2m%v * profiles_tl(j)%s2m%u)
        ELSE
          wangl_tl = 0._jprb
        ENDIF
!-----------Find the index of the pressure level that is closest to the surface-
        ilevwrt = aux%s(j)%nearestlev_surf
!------------------
!redundant sections
!            if (aux(j) %pfraction_surf < 0._jprb) then
!              ilevwrt=ilevwrt+1
!            end if
!
!            if(ilevwrt.gt.100)then
!              ilevsur=ilevwrt-1
!            else
!              ilevsur=ilevwrt
!            endif
!------------------
! level index for the raytracing angles to be used for surface
        ilevsur = ilevwrt
!-----------Compute the secant of the satellite viewing angle at the surface----
        IF (raytracing%pathsat(ilevsur - 1, j) .NE. 1) THEN
          zensat_tl = raytracing_tl%pathsat(ilevsur - 1, j) * (1._jprb / raytracing%pathsat(ilevsur - 1, j)) *      &
            & (1._jprb / SQRT(raytracing%pathsat(ilevsur - 1, j) ** 2_jpim - 1._jprb))
        ELSE
          zensat_tl = 0._jprb
        ENDIF
!-----------Compute the secant of the sun zenith angle at the surface-----------
        IF (raytracing%pathsun(ilevsur - 1, j) .NE. 1._jprb) THEN
          zensun_tl = raytracing_tl%pathsun(ilevsur - 1, j) * (1._jprb / raytracing%pathsun(ilevsur - 1, j)) *      &
            & (1._jprb / SQRT(raytracing%pathsun(ilevsur - 1, j) ** 2_jpim - 1._jprb))
        ELSE
          zensun_tl = 0._jprb
        ENDIF
        IF (sunglint%s(j)%windsp == 0._jprb) THEN
          IF (sunglint%s(j)%zensat == sunglint%s(j)%zensun .AND. sunglint%s(j)%dazng == 180.0_jprb) THEN
            sunglint_tl%s(j)%glint = 0._jprb
          ELSE
            sunglint_tl%s(j)%glint = 0._jprb
          ENDIF
        ELSE IF (sunglint%s(j)%windsp /= 0._jprb) THEN
!-----------Set up some parameters used to define the JONSWAP frequency -------
!           spectrum of the surface wave                                       |
!------------------------------------------------------------------------------
          x_u_tl     = profiles_tl(j)%s2m%wfetc * gravity / sunglint%s(j)%windsp ** 2_jpim -      &
            & 2._jprb * windsp_tl * profiles(j)%s2m%wfetc * gravity / sunglint%s(j)%windsp ** 3_jpim
          alfa1_tl   =  - 0.22_jprb * X_U_tl * sunglint%s(j)%alfa1 / sunglint%s(j)%x_u
          omega_m_tl =  - windsp_tl * sunglint%s(j)%omega_m / sunglint%s(j)%windsp - &
            & x_u_tl * 0.33_jprb * sunglint%s(j)%omega_m / sunglint%s(j)%x_u
!-----------Compute the total variance of the slope of the facet---------------
!                                                                              |
!           To compute the total variance of the slope of the facet (or the    |
!           total variance of the displacement of the water surface) the       |
!           knowledge of the frequency spectrum of the surface wave and the    |
!           inverse function of the dispersive relation of the full-gravityity-   |
!           capillary wave is needed.                                          |
!                                                                              |
!           The dispersive relation of the full-gravityity-capillary wave is      |
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
!           full-gravityity-capillary wave has been pre-computed numerically      |
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
          alfa1_r = 1._jprb / sunglint%s(j)%alfa1
          omega_m_3 = sunglint%s(j)%omega_m ** 3_jpim
          psi_tl(1) = 0._jprb
          FF_tl(1) = 0._jprb
          DO I = 2, coef%ws_nomega
            omega_1 = omega_1_array(i)

            IF (i == 1) THEN
              psi_tl(i) = 0._jprb
            ELSE
              IF (sunglint%s(j)%omega_m < max_exp_exponent_4rt * omega_1) THEN ! stop exp overflow

                IF (sunglint%beta(i, j) >= exp_max_exp_exponent_r) THEN
                  IF (omega_1 <= sunglint%s(j)%omega_m) THEN
                    sigma2 = sigma_a ** 2_jpim
                  ELSE
                    sigma2 = sigma_b ** 2_jpim
                  ENDIF

                  beta_tl = omega_m_tl * sunglint%beta(i, j) * omega_1 * &
                      (omega_1 - sunglint%s(j)%omega_m) / &
                      (sigma2 * omega_m_3)
                ELSE
                  beta_tl = 0._jprb
                ENDIF

                psi_tl(i) = alfa1_tl * sunglint%psi(i, j) * alfa1_r - &
                    omega_m_tl * sunglint%psi(i, j) * 5._jprb * omega_m_3 * omega_1_4_r(i)

                IF (sunglint%beta(i, j) > min_exponent) psi_tl(i) = psi_tl(i) + beta_tl * sunglint%psi(i, j) * logk

               ELSE
                 psi_tl(i) = 0._jprb
               ENDIF
            ENDIF
            FF_tl(i) = psi_tl(i) * coef%ws_k_omega(i) ** 2_jpim
          ENDDO
!-----------Compute the Simpson's integral--------------------------------------
          bb    = 301._jprb
          aa    = 0._jprb
          np    = 1501_jpim
          hinc  = (bb - aa) / REAL(2 * np, jprb)
          sa_tl = ff_tl(1)
          sb_tl = ff_tl(coef%ws_nomega)
          s4_tl = 0._jprb
          DO I = 1, NP
            jj    = 2 * I - 1
            s4_tl = s4_tl + ff_tl(jj)
          ENDDO
          S2_tl = 0._jprb
          DO I = 1, NP
            jj    = 2 * I
            s2_tl = s2_tl + ff_tl(jj)
          ENDDO
          gamma_sq_tl = hinc / 3.0_jprb * (sa_tl + sb_tl + 4.0_jprb * s4_tl + 2.0_jprb * s2_tl)
!-------------------------------------------------------------------------------
!-----------Compute the rms of the slope of the facet along the X (gamma_O)----
!           and Y (gamma_P) axis.                                              |
!------------------------------------------------------------------------------
          gamma_O_tl  = GAMMA_SQ_tl * (2._jprb + cos(2._jprb * sunglint%s(j)%wangl)) / (8._jprb * sunglint%s(j)%gamma_O)     &
            &  - wangl_tl * sin(2._jprb * sunglint%s(j)%wangl) * sunglint%s(j)%gamma_sq / (4._jprb * sunglint%s(j)%gamma_O)
          gamma_P_tl  = GAMMA_SQ_tl * (2._jprb - cos(2._jprb * sunglint%s(j)%wangl)) / (8._jprb * sunglint%s(j)%gamma_P) +      &
            & wangl_tl * sin(2._jprb * sunglint%s(j)%wangl) * sunglint%s(j)%gamma_sq / (4._jprb * sunglint%s(j)%gamma_P)
!-----------Obtain angles csi and alfa------------------------------------------
          csi_tl      =                                                                                                  &
            & zensun_tl * cos(PI - sunglint%s(j)%dazng * deg2rad) * (1._jprb / (cos(sunglint%s(j)%zensun)) ** 2_jpim) *  &
            & 1._jprb / (1._jprb + (tan(sunglint%s(j)%zensun) * cos(PI - sunglint%s(j)%dazng * deg2rad)) ** 2_jpim)
          IF ((sin(sunglint%s(j)%zensun) * sin(PI - sunglint%s(j)%dazng * deg2rad)) ** 2_jpim .NE. 1._jprb) THEN
            alfa_tl = zensun_tl * (cos(sunglint%s(j)%zensun) * sin(PI - sunglint%s(j)%dazng * deg2rad) *      &
              & (1._jprb / SQRT(1 - (sin(sunglint%s(j)%zensun) * sin(PI - sunglint%s(j)%dazng * deg2rad)) ** 2_jpim)))
          ELSE
            alfa_tl = 0._jprb
          ENDIF
          IF ((sunglint%s(j)%csi + sunglint%s(j)%zensat) >= 0._jprb) THEN
            sunglint_tl%s(j)%omega = (csi_tl + zensat_tl) / 2._jprb
          ELSE
            sunglint_tl%s(j)%omega =  - (csi_tl + zensat_tl) / 2._jprb
          ENDIF
!-----------Compute the value of the function that represents the shadowing  --
!           of the incident and reflected ray                                  |
!------------------------------------------------------------------------------
!-----------First order shadowing (the slope of the facet is negative)---------
          gammax_tl =                                                                                                     &
            & (1._jprb / (cos((sunglint%s(j)%csi - sunglint%s(j)%zensat) / 2._jprb)) ** 2_jpim) * (csi_tl - zensat_tl) /  &
            & 2._jprb
          IF (abs(sunglint%s(j)%zensat) > 0._jprb) THEN
            IF ((1._jprb / tan(abs(sunglint%s(j)%zensat)) - &
              & sunglint%s(j)%gammax * sunglint%s(j)%zensat / abs(sunglint%s(j)%zensat)) >= 0._jprb) THEN
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
              Q_shad_tl = 0._jprb
            ELSE IF (sunglint%s(j)%zensat >= abs(sunglint%s(j)%csi)) THEN
              a_shad_tl   =                                                                                               &
                &  - zensat_tl * (1._jprb / (sin(sunglint%s(j)%zensat)) ** 2_jpim) * (1._jprb / sunglint%s(j)%gamma_O) -  &
                & gamma_O_tl * (1._jprb / tan(sunglint%s(j)%zensat)) * (1._jprb / sunglint%s(j)%gamma_O ** 2_jpim)
              lambda_a_tl =  - a_shad_tl * (exp( - sunglint%s(j)%a_shad ** 2_jpim / 2._jprb)) /     &
                & (fac * sunglint%s(j)%a_shad ** 2_jpim)
              Q_shad_A_tl =  - 1._jprb / (1._jprb + sunglint%s(j)%lambda_a) ** 2_jpim * (lambda_a_tl)
              Q_shad_tl   = Q_shad_A_tl * theta_fi
            ELSE IF (abs(sunglint%s(j)%csi) > sunglint%s(j)%zensat) THEN
              IF (sunglint%s(j)%csi >= 0._jprb) THEN
                b_shad_tl = ( - 1._jprb / (sin(sunglint%s(j)%csi)) ** 2_jpim * 1._jprb / sunglint%s(j)%gamma_O) * csi_tl     &
                  &  - ((1._jprb / tan(sunglint%s(j)%csi)) / sunglint%s(j)%gamma_O ** 2_jpim) * gamma_O_tl
              ELSE
                b_shad_tl =                                                                                           &
                  & (1._jprb / (sin(sunglint%s(j)%csi)) ** 2_jpim * 1._jprb / sunglint%s(j)%gamma_O) * csi_tl +  &
                  & ((1._jprb / tan(sunglint%s(j)%csi)) / sunglint%s(j)%gamma_O ** 2_jpim) * gamma_O_tl
              ENDIF
              lambda_b_tl =  - b_shad_tl * (exp( - sunglint%s(j)%b_shad ** 2_jpim / 2._jprb)) /      &
                & (fac * sunglint%s(j)%b_shad ** 2_jpim)
              Q_shad_B_tl =  - 1._jprb / (1._jprb + sunglint%s(j)%lambda_b) ** 2_jpim * (lambda_b_tl)
              Q_shad_tl   = Q_shad_B_tl * theta_csi
            ENDIF
          ELSE
            a_shad_tl =                                                                                                 &
              &  - zensat_tl * (1._jprb / (sin(sunglint%s(j)%zensat)) ** 2_jpim) * (1._jprb / sunglint%s(j)%gamma_O) -  &
              & gamma_O_tl * (1._jprb / tan(sunglint%s(j)%zensat)) * (1._jprb / sunglint%s(j)%gamma_O ** 2_jpim)
            IF (sunglint%s(j)%csi >= 0._jprb) THEN
              b_shad_tl = ( - 1._jprb / (sin(sunglint%s(j)%csi)) ** 2_jpim * 1._jprb / sunglint%s(j)%gamma_O) * csi_tl     &
                &  - ((1._jprb / tan(sunglint%s(j)%csi)) / sunglint%s(j)%gamma_O ** 2_jpim) * gamma_O_tl
            ELSE
              b_shad_tl = (1._jprb / (sin(sunglint%s(j)%csi)) ** 2_jpim * 1._jprb / sunglint%s(j)%gamma_O) * csi_tl     &
                &  + ((1._jprb / tan(sunglint%s(j)%csi)) / sunglint%s(j)%gamma_O ** 2_jpim) * gamma_O_tl
            ENDIF
            lambda_a_tl =  - a_shad_tl * (exp( - sunglint%s(j)%a_shad ** 2_jpim / 2._jprb)) /      &
              & (fac * sunglint%s(j)%a_shad ** 2_jpim)
            lambda_b_tl =  - b_shad_tl * (exp( - sunglint%s(j)%b_shad ** 2_jpim / 2._jprb)) /      &
              & (fac * sunglint%s(j)%b_shad ** 2_jpim)
            Q_shad_tl   =  - 1._jprb / (1._jprb + sunglint%s(j)%lambda_a + sunglint%s(j)%lambda_b) ** 2_jpim *      &
              & (lambda_a_tl + lambda_b_tl)
            Q_shad_tl   = Q_shad_tl * theta_fi * theta_csi
          ENDIF
!-------------------------------------------------------------------------------
!-----------Compute the probability density that the slope of the facet at ----
!           a certain point is gammax when the incident ray and the ray from   |
!           the reflection of the incident ray on the local surface do not     |
!           intersect with any other surface.                                  |
!------------------------------------------------------------------------------
          IF ((cos(sunglint%s(j)%csi) + cos(sunglint%s(j)%zensat)) >= 0._jprb) THEN
            c_shad_tl = gamma_P_tl * (cos(sunglint%s(j)%csi) + cos(sunglint%s(j)%zensat)) -      &
              & csi_tl * sunglint%s(j)%gamma_P * sin(sunglint%s(j)%csi) -                           &
              & zensat_tl * sunglint%s(j)%gamma_P * sin(sunglint%s(j)%zensat)
          ELSE
            c_shad_tl = - gamma_P_tl * (cos(sunglint%s(j)%csi) + cos(sunglint%s(j)%zensat)) +      &
              & csi_tl * sunglint%s(j)%gamma_P * sin(sunglint%s(j)%csi) +                           &
              & zensat_tl * sunglint%s(j)%gamma_P * sin(sunglint%s(j)%zensat)
          ENDIF
          p_prime_tl             =  - alfa_tl * sunglint%s(j)%p_prime * tan(sunglint%s(j)%alfa) /               &
            & (cos(sunglint%s(j)%alfa) * sunglint%s(j)%c_shad) ** 2_jpim + c_shad_tl * sunglint%s(j)%p_prime *  &
            & ((tan(sunglint%s(j)%alfa)) ** 2_jpim / sunglint%s(j)%c_shad ** 3_jpim - 1._jprb / sunglint%s(j)%c_shad)
          g_shad_tl              = - 0.5_jprb * (csi_tl - zensat_tl) * tan(sunglint%s(j)%zensat) / &
            & cos(0.5_jprb * (sunglint%s(j)%csi - sunglint%s(j)%zensat)) ** 2_jpim - &
            & zensat_tl * tan(0.5_jprb * (sunglint%s(j)%csi - sunglint%s(j)%ZENSAT)) / &
            & cos(sunglint%s(j)%zensat) ** 2_jpim
          fac1_tl                =  - alfa_tl * 2._jprb * sunglint%s(j)%fac1 * tan(sunglint%s(j)%alfa) -           &
            & (csi_tl - zensat_tl) * sunglint%s(j)%fac1 * tan((sunglint%s(j)%csi - sunglint%s(j)%zensat) / 2._jprb)
          Pxy_gammaxy_tl         = gamma_O_tl * sunglint%s(j)%Pxy_gammaxy *                                           &
            & (sunglint%s(j)%gammax ** 2_jpim / sunglint%s(j)%gamma_O ** 3_jpim - 1._jprb / sunglint%s(j)%gamma_O) -  &
            & gammax_tl * sunglint%s(j)%Pxy_gammaxy * sunglint%s(j)%gammax / sunglint%s(j)%gamma_O ** 2_jpim
          RQS_tl                 =                                                                                         &
            & Q_shad_tl * sunglint%s(j)%g_shad * sunglint%s(j)%p_prime * sunglint%s(j)%Pxy_gammaxy / sunglint%s(j)%fac1 +  &
            & g_shad_tl * sunglint%s(j)%Q_shad * sunglint%s(j)%p_prime * sunglint%s(j)%Pxy_gammaxy / sunglint%s(j)%fac1 +  &
            & p_prime_tl * sunglint%s(j)%Q_shad * sunglint%s(j)%g_shad * sunglint%s(j)%Pxy_gammaxy / sunglint%s(j)%fac1 +  &
            & Pxy_gammaxy_tl * sunglint%s(j)%Q_shad * sunglint%s(j)%g_shad * sunglint%s(j)%p_prime / sunglint%s(j)%fac1 -  &
            & fac1_tl * sunglint%s(j)%Q_shad * sunglint%s(j)%g_shad * sunglint%s(j)%p_prime * sunglint%s(j)%Pxy_gammaxy /  &
            & sunglint%s(j)%fac1 ** 2_jpim
!-------------------------------------------------------------------------------
          sunglint_tl%s(j)%glint = rqs_tl! *PI*EPSILON**2
        ENDIF
      ENDIF
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_REFSUN_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_refsun_tl
