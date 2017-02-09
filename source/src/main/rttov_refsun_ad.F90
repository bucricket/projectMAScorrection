SUBROUTINE rttov_refsun_ad( &
            & profiles,      &
            & profiles_ad,   &
            & coef,          &
            & aux,           &
            & sunglint,      &
            & sunglint_ad,   &
            & raytracing,    &
            & raytracing_ad)
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
!     the Joint North Sea Wave ProJect (JONSWAP) wave spectral model.
!     The total variance of the slope of the facet is obtained from the
!     frequency spectrum of the surface wave and the inverse function of the
!     dispersive relation of the full-gravityity-capillary wave.
!     The shadowing of the surface of the facets on the backsides of the waves
!     and deep in the throughs between waves is considered.
!     A coordinate system is used where the average water surface lies in the
!     X-Y plane. The coordinate system is right-handed and is located at the
!     reflection point. The Z axis points toward the zenith and the X axis
!     points in the direction formed by the proJection of the reflected ray
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
!     2            15/08/2009  Determination of profile level index for the
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
  TYPE(profile_type   ), INTENT(IN)    :: profiles   (:)
  TYPE(profile_type   ), INTENT(INOUT) :: profiles_ad(size(profiles))
  TYPE(profile_aux    ), INTENT(IN)    :: aux
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing_ad
  TYPE(sunglint_type  ), INTENT(IN)    :: sunglint
  TYPE(sunglint_type  ), INTENT(INOUT) :: sunglint_ad
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
!INTF_END
!     End of subroutine arguments
!       Local arrays
  REAL   (KIND=jprb) :: psi_ad(coef%ws_nomega)              ! The frequency spectrum of the surface wave
  REAL   (KIND=jprb) :: ff_ad (coef%ws_nomega)              ! Working space
!       Local scalars:
  INTEGER(KIND=jpim) :: j, jj, i, np, ilevsur, ilevwrt
  REAL   (KIND=jprb) :: csi_ad                              ! Angle between the zenith and the proJection of the
! incident ray on the x-z plane
  REAL   (KIND=jprb) :: alfa_ad                             ! Angle between the incident ray and the x-z plane
  REAL   (KIND=jprb) :: c_shad_ad                           ! The average magnitude of tan(alfa)
  REAL   (KIND=jprb) :: p_prime_ad                          ! The probability density of tan(alfa)
  REAL   (KIND=jprb) :: pxy_gammaxy_ad                      ! The Joint  probability density of the along-view and
! cross view slope
  REAL   (KIND=jprb) :: rqs_ad                              ! Effective distribution function
  REAL   (KIND=jprb) :: gamma_sq_ad                         ! Total variance of the slope of the facet
  REAL   (KIND=jprb) :: gamma_o_ad                          ! The mean square of the along-view (X axis) slope
  REAL   (KIND=jprb) :: gamma_p_ad                          ! The mean square of the cross-view (Y axis) slope
  REAL   (KIND=jprb) :: g_shad_ad                           ! Normalization function
  REAL   (KIND=jprb) :: gammax_ad                           ! The x-slope of the facet
  REAL   (KIND=jprb) :: theta_fi                            ! First order shadowing factor for the reflected ray
  REAL   (KIND=jprb) :: theta_csi                           ! First ordesr shadowing factor for the incident ray
  REAL   (KIND=jprb) :: q_shad_ad                           ! Second order shadowing factor
  REAL   (KIND=jprb) :: q_shad_a_ad                         ! Second order shadowing factor
  REAL   (KIND=jprb) :: q_shad_b_ad                         ! Second order shadowing factor
  REAL   (KIND=jprb) :: zensat_ad                           ! Zenith angle of satellite viewing angle at surface
  REAL   (KIND=jprb) :: zensun_ad                           ! Zenith angle of sun at surface
  REAL   (KIND=jprb) :: omega_1                             ! Frequency of the surface wave
  REAL   (KIND=jprb) :: fac1_ad
  REAL   (KIND=jprb) :: a_shad_ad
  REAL   (KIND=jprb) :: b_shad_ad
  REAL   (KIND=jprb) :: lambda_a_ad
  REAL   (KIND=jprb) :: lambda_b_ad
  REAL   (KIND=jprb) :: x_u_ad
  REAL   (KIND=jprb) :: alfa1_ad
  REAL   (KIND=jprb) :: k
  REAL   (KIND=jprb) :: sigma_a
  REAL   (KIND=jprb) :: sigma_b
  REAL   (KIND=jprb) :: omega_m_ad
  REAL   (KIND=jprb) :: sigma2
  REAL   (KIND=jprb) :: beta_ad
  REAL   (KIND=jprb) :: windsp_ad
  REAL   (KIND=jprb) :: wangl_ad
  REAL   (KIND=jprb) :: s2_ad, s4_ad, sa_ad, sb_ad
  REAL   (KIND=jprb) :: bb, aa, hinc, fac
  INTEGER(KIND=jpim) :: nprofiles

  REAL   (KIND=jprb) :: max_exp_exponent_4rt, alfa1_r, logk
  REAL   (KIND=jprb) :: exp_max_exp_exponent_r, omega_m_3
  REAL   (KIND=jprb) :: omega_1_array(coef%ws_nomega)
  REAL   (KIND=jprb) :: omega_1_4_r(coef%ws_nomega)

  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_REFSUN_AD', 0_jpim, ZHOOK_HANDLE)
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

  DO J = 1, nprofiles
    IF (profiles(j)%sunzenangle >= 0.0 .AND. &
        profiles(j)%sunzenangle < max_sol_zen) THEN
      IF (profiles(j)%skin%surftype == surftype_sea) THEN
        windsp_ad      = 0._jprb
        wangl_ad       = 0._jprb
        OMEGA_M_ad     = 0._jprb
        RQS_ad         = 0._jprb
        Q_shad_ad      = 0._jprb
        Q_shad_A_ad    = 0._jprb
        Q_shad_B_ad    = 0._jprb
        g_shad_ad      = 0._jprb
        p_prime_ad     = 0._jprb
        Pxy_gammaxy_ad = 0._jprb
        fac1_ad        = 0._jprb
        gamma_O_ad     = 0._jprb
        gammax_ad      = 0._jprb
        alfa_ad        = 0._jprb
        csi_ad         = 0._jprb
        zensat_ad      = 0._jprb
        zensun_ad      = 0._jprb
        c_shad_ad      = 0._jprb
        gamma_p_ad     = 0._jprb
        gamma_sq_ad    = 0._jprb
        lambda_a_ad    = 0._jprb
        lambda_b_ad    = 0._jprb
        b_shad_ad      = 0._jprb
        a_shad_ad      = 0._jprb
        sa_ad          = 0._jprb
        sb_ad          = 0._jprb
        s4_ad          = 0._jprb
        s2_ad          = 0._jprb
        ff_ad(:)       = 0._jprb
        beta_ad        = 0._jprb
        alfa1_ad       = 0._jprb
        psi_ad(:)      = 0._jprb
        x_u_ad         = 0._jprb

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
        IF (sunglint%s(j)%windsp == 0._jprb) THEN
          IF (sunglint%s(j)%zensat == sunglint%s(j)%zensun .AND. sunglint%s(j)%dazng == 180.0_jprb) THEN
            sunglint_ad%s(j)%glint = 0._jprb
          ELSE
            sunglint_ad%s(j)%glint = 0._jprb
          ENDIF
        ELSE
          rqs_ad         = rqs_ad + sunglint_ad%s(j)%glint
!write(68,*)' rqs_ad, sunglint_ad%glint ',rqs_ad,sunglint_ad(j)%glint
!-----------Compute the probability density that the slope of the facet at ----
!           a certain point is gammax when the incident ray and the ray from   |
!           the reflection of the incident ray on the local surface do not     |
!           intersect with any other surface.                                  |
!------------------------------------------------------------------------------
          q_shad_ad      = q_shad_ad +      &
            & RQS_ad * sunglint%s(j)%g_shad * sunglint%s(j)%p_prime * sunglint%s(j)%Pxy_gammaxy / sunglint%s(j)%fac1
          g_shad_ad      = g_shad_ad +      &
            & RQS_ad * sunglint%s(j)%q_shad * sunglint%s(j)%p_prime * sunglint%s(j)%Pxy_gammaxy / sunglint%s(j)%fac1
          p_prime_ad     = p_prime_ad +      &
            & RQS_ad * sunglint%s(j)%q_shad * sunglint%s(j)%g_shad * sunglint%s(j)%Pxy_gammaxy / sunglint%s(j)%fac1
          Pxy_gammaxy_ad = Pxy_gammaxy_ad +      &
            & RQS_ad * sunglint%s(j)%q_shad * sunglint%s(j)%g_shad * sunglint%s(j)%p_prime / sunglint%s(j)%fac1
          fac1_ad        = fac1_ad -                                                                                      &
            & RQS_ad * sunglint%s(j)%q_shad * sunglint%s(j)%g_shad * sunglint%s(j)%p_prime * sunglint%s(j)%Pxy_gammaxy /  &
            & sunglint%s(j)%fac1 ** 2_jpim
          RQS_ad         = 0._JPRB

          gamma_O_ad     = gamma_O_ad + Pxy_gammaxy_ad * sunglint%s(j)%Pxy_gammaxy *      &
            & (sunglint%s(j)%gammax ** 2_jpim / sunglint%s(j)%gamma_O ** 3_jpim - 1._jprb / sunglint%s(j)%gamma_O)
          gammax_ad      = gammax_ad -      &
            & Pxy_gammaxy_ad * sunglint%s(j)%Pxy_gammaxy * sunglint%s(j)%gammax / sunglint%s(j)%gamma_O ** 2_jpim
          Pxy_gammaxy_ad = 0._JPRB

          alfa_ad        = alfa_ad - fac1_ad * 2._jprb * sunglint%s(j)%fac1 * tan(sunglint%s(j)%alfa)
          csi_ad         =      &
            & csi_ad - fac1_ad * sunglint%s(j)%fac1 * tan((sunglint%s(j)%csi - sunglint%s(j)%zensat) / 2._jprb)
          zensat_ad      =      &
            & zensat_ad + fac1_ad * sunglint%s(j)%fac1 * tan((sunglint%s(j)%csi - sunglint%s(j)%zensat) / 2._JPRB)
          fac1_ad        = 0._JPRB

          csi_ad         = csi_ad - g_shad_ad * 0.5_jprb * tan(sunglint%s(j)%zensat) /      &
            & cos(0.5_jprb * (sunglint%s(j)%csi - sunglint%s(j)%zensat)) ** 2_jpim
          zensat_ad      = zensat_ad + 0.5_jprb * g_shad_ad * tan(sunglint%s(j)%zensat) / &
            & cos(0.5_jprb * (sunglint%s(j)%csi - sunglint%s(j)%zensat)) ** 2_jpim
          zensat_ad      = zensat_ad - g_shad_ad * tan(0.5_jprb * (sunglint%s(j)%csi - sunglint%s(j)%ZENSAT)) / &
            & cos(sunglint%s(j)%zensat) ** 2_jpim
          g_shad_ad      = 0._JPRB
          alfa_ad        = alfa_ad - p_prime_ad * sunglint%s(j)%p_prime * tan(sunglint%s(j)%alfa) /      &
            & (cos(sunglint%s(j)%alfa) * sunglint%s(j)%c_shad) ** 2_jpim
          c_shad_ad      = c_shad_ad + p_prime_ad * sunglint%s(j)%p_prime *      &
            & ((tan(sunglint%s(j)%alfa)) ** 2_jpim / sunglint%s(j)%c_shad ** 3_jpim - 1._jprb / sunglint%s(j)%c_shad)
          p_prime_ad     = 0._jprb

          IF ((cos(sunglint%s(j)%csi) + cos(sunglint%s(j)%zensat)) >= 0._jprb) THEN
            zensat_ad  = zensat_ad - c_shad_ad * sunglint%s(j)%gamma_P * sin(sunglint%s(j)%zensat)
            csi_ad     = csi_ad - c_shad_ad * sunglint%s(j)%gamma_P * sin(sunglint%s(j)%csi)
            gamma_P_ad = gamma_P_ad + c_shad_ad * (cos(sunglint%s(j)%csi) + cos(sunglint%s(j)%zensat))
            c_shad_ad  = 0._jprb
          ELSE
            zensat_ad  = zensat_ad + c_shad_ad * sunglint%s(j)%gamma_P * sin(sunglint%s(j)%zensat)
            csi_ad     = csi_ad + c_shad_ad * sunglint%s(j)%gamma_P * sin(sunglint%s(j)%csi)
            gamma_P_ad = gamma_P_ad - c_shad_ad * (cos(sunglint%s(j)%csi) + cos(sunglint%s(j)%zensat))
            c_shad_ad  = 0._jprb
          ENDIF
!-----------Compute the value of the function that represents the shadowing  --
!           of the incident and reflected ray                                  |
!------------------------------------------------------------------------------
!-----------First order shadowing (the slope of the facet is negative)----------
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
              & sunglint%s(j)%gammax * sunglint%s(j)%csi / abs(sunglint%s(j)%csi)) >= 0._jprb) THEN
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
              Q_shad_ad = 0._jprb
            ELSE IF (sunglint%s(j)%zensat >= abs(sunglint%s(j)%csi)) THEN
              Q_shad_A_ad = Q_shad_A_ad + Q_shad_ad * theta_fi
              Q_shad_ad   = 0._jprb
              lambda_a_ad = lambda_a_ad - Q_shad_A_ad / (1._jprb + sunglint%s(j)%lambda_a) ** 2_jpim
              Q_shad_A_ad = 0._jprb
              a_shad_ad   = a_shad_ad - lambda_a_ad * (exp( - sunglint%s(j)%a_shad ** 2_jpim / 2._jprb)) /      &
                & (fac * sunglint%s(j)%a_shad ** 2_jpim)
              lambda_a_ad = 0._jprb
              zensat_ad   = zensat_ad -      &
                & a_shad_ad * (1._jprb / (sin(sunglint%s(j)%zensat)) ** 2_jpim) * (1._jprb / sunglint%s(j)%gamma_O)
              gamma_O_ad  = gamma_O_ad -      &
                & a_shad_ad * (1._jprb / tan(sunglint%s(j)%zensat)) * (1._jprb / sunglint%s(j)%gamma_O ** 2_jpim)
              a_shad_ad   = 0._jprb
            ELSE IF (abs(sunglint%s(j)%csi) > sunglint%s(j)%zensat) THEN
              Q_shad_B_ad = Q_shad_B_ad + Q_shad_ad * theta_csi
              Q_shad_ad   = 0._jprb
              lambda_b_ad = lambda_b_ad - Q_shad_B_ad / (1._jprb + sunglint%s(j)%lambda_b) ** 2_jpim
              Q_shad_B_ad = 0._jprb
              b_shad_ad   = b_shad_ad - lambda_b_ad * (exp( - sunglint%s(j)%b_shad ** 2_jpim / 2._jprb)) /      &
                & (fac * sunglint%s(j)%b_shad ** 2_jpim)
              lambda_b_ad = 0._jprb
              IF (sunglint%s(j)%csi >= 0._jprb) THEN
                csi_ad     =      &
                  & csi_ad + b_shad_ad * ( - 1._jprb / (sin(sunglint%s(j)%csi)) ** 2_jpim * 1._jprb / sunglint%s(j)%gamma_O)
                gamma_O_ad =      &
                  & gamma_O_ad - b_shad_ad * ((1._jprb / tan(sunglint%s(j)%csi)) / sunglint%s(j)%gamma_O ** 2_jpim)
                b_shad_ad  = 0._jprb
              ELSE IF (sunglint%s(j)%csi < 0._jprb) THEN
                csi_ad     =      &
                  & csi_ad + b_shad_ad * (1 / (sin(sunglint%s(j)%csi)) ** 2_jpim * 1._jprb / sunglint%s(j)%gamma_O)
                gamma_O_ad =      &
                  & gamma_O_ad + b_shad_ad * ((1._jprb / tan(sunglint%s(j)%csi)) / sunglint%s(j)%gamma_O ** 2_jpim)
                b_shad_ad  = 0._jprb
              ENDIF
            ENDIF
          ELSE
            q_shad_ad   = q_shad_ad * theta_fi * theta_csi
            lambda_a_ad =      &
              & lambda_a_ad - Q_shad_ad / (1._jprb + sunglint%s(j)%lambda_a + sunglint%s(j)%lambda_b) ** 2_jpim
            lambda_b_ad =      &
              & lambda_b_ad - Q_shad_ad / (1._jprb + sunglint%s(j)%lambda_a + sunglint%s(j)%lambda_b) ** 2_jpim
            Q_shad_ad   = 0._jprb
            b_shad_ad   = b_shad_ad - lambda_b_ad * (exp( - sunglint%s(j)%b_shad ** 2_jpim / 2._jprb)) /      &
              & (fac * sunglint%s(j)%b_shad ** 2_jpim)
            lambda_b_ad = 0._jprb
            a_shad_ad   = a_shad_ad - lambda_a_ad * (exp( - sunglint%s(j)%a_shad ** 2_jpim / 2._jprb)) /      &
              & (fac * sunglint%s(j)%a_shad ** 2_jpim)
            lambda_a_ad = 0._jprb
            IF (sunglint%s(j)%csi >= 0._jprb) THEN
              csi_ad     =      &
                & csi_ad + b_shad_ad * ( - 1._jprb / (sin(sunglint%s(j)%csi)) ** 2_jpim * 1._jprb / sunglint%s(j)%gamma_O)
              gamma_O_ad =      &
                & gamma_O_ad - b_shad_ad * ((1._jprb / tan(sunglint%s(j)%csi)) / sunglint%s(j)%gamma_O ** 2_jpim)
              b_shad_ad  = 0._jprb
            ELSE IF (sunglint%s(j)%csi < 0._jprb) THEN
              csi_ad     =      &
                & csi_ad + b_shad_ad * (1._jprb / (sin(sunglint%s(j)%csi)) ** 2_jpim * 1._jprb / sunglint%s(j)%gamma_O)
              gamma_O_ad =      &
                & gamma_O_ad + b_shad_ad * ((1._jprb / tan(sunglint%s(j)%csi)) / sunglint%s(j)%gamma_O ** 2_jpim)
              b_shad_ad  = 0._jprb
            ENDIF
            zensat_ad  = zensat_ad -      &
              & a_shad_ad * (1._jprb / (sin(sunglint%s(j)%zensat)) ** 2_jpim) * (1._jprb / sunglint%s(j)%gamma_O)
            gamma_O_ad =      &
              & gamma_O_ad - a_shad_ad * (1._jprb / tan(sunglint%s(j)%zensat)) * (1._jprb / sunglint%s(j)%gamma_O ** 2_jpim)
            a_shad_ad  = 0._jprb
          ENDIF
!-----------First order shadowing (the slope of the facet is negative)----------
          csi_ad    = csi_ad +      &
            & gammax_ad * (1._jprb / (cos((sunglint%s(j)%csi - sunglint%s(j)%zensat) / 2._jprb)) ** 2_jpim) / 2._jprb
          zensat_ad = zensat_ad -      &
            & gammax_ad * (1._jprb / (cos((sunglint%s(j)%csi - sunglint%s(j)%zensat) / 2._jprb)) ** 2_jpim) / 2._jprb
          gammax_ad = 0._jprb
!-----------Obtain angles csi and alfa------------------------------------------
          IF ((sunglint%s(j)%csi + sunglint%s(j)%zensat) >= 0._jprb) THEN
            csi_ad    = csi_ad + sunglint_ad%s(j)%omega / 2._jprb
            zensat_ad = zensat_ad + sunglint_ad%s(j)%omega / 2._jprb
          ELSE
            csi_ad    = csi_ad - sunglint_ad%s(j)%omega / 2._jprb
            zensat_ad = zensat_ad - sunglint_ad%s(j)%omega / 2._jprb
          ENDIF
          IF ((sin(sunglint%s(j)%zensun) * sin(PI - sunglint%s(j)%dazng * deg2rad)) ** 2_jpim .NE. 1._jprb) THEN
            zensun_ad = zensun_ad + alfa_ad * (cos(sunglint%s(j)%zensun) * sin(PI - sunglint%s(j)%dazng * deg2rad) *      &
              & (1._jprb / SQRT(1 - (sin(sunglint%s(j)%zensun) * sin(PI - sunglint%s(j)%dazng * deg2rad)) ** 2_jpim)))
          ELSE
            alfa_ad = 0._jprb
          ENDIF
          zensun_ad   = zensun_ad +                                                                                          &
            & csi_ad * cos(PI - sunglint%s(j)%dazng * deg2rad) * (1._jprb / (cos(sunglint%s(j)%zensun)) ** 2_jpim) * 1._jprb &
            &  / (1._jprb + (tan(sunglint%s(j)%zensun) * cos(PI - sunglint%s(j)%dazng * deg2rad)) ** 2_jpim)
          csi_ad      = 0._jprb
!-----------Compute the rms of the slope of the facet along the X (gamma_O)----
!           and Y (gamma_P) axis.                                              |
!-------------------------------------------------------------------------------
          gamma_sq_ad =      &
            & gamma_sq_ad + gamma_P_ad * (2._jprb - cos(2._jprb * sunglint%s(j)%wangl)) / (8._jprb * sunglint%s(j)%gamma_P)
          wangl_ad    = wangl_ad +      &
            & gamma_P_ad * sin(2._jprb * sunglint%s(j)%wangl) * sunglint%s(j)%gamma_sq / (4._jprb * sunglint%s(j)%gamma_P)
          gamma_P_ad  = 0._jprb
          gamma_sq_ad =      &
            & gamma_sq_ad + gamma_O_ad * (2._jprb + cos(2._jprb * sunglint%s(j)%wangl)) / (8._jprb * sunglint%s(j)%gamma_O)
          wangl_ad    = wangl_ad -      &
            & gamma_O_ad * sin(2._jprb * sunglint%s(j)%wangl) * sunglint%s(j)%gamma_sq / (4._jprb * sunglint%s(j)%gamma_O)
          gamma_O_ad  = 0._jprb
!-----------Compute the Simpson's integral--------------------------------------
          BB = 301._jprb
          AA = 0._jprb
          NP = 1501_jpim
          HINC        = (BB - AA) / REAL(2 * NP, jprb)
          SA_ad       = SA_ad + GAMMA_SQ_ad * HINC / 3.0_jprb
          SB_ad       = SB_ad + GAMMA_SQ_ad * HINC / 3.0_jprb
          S4_ad       = S4_ad + GAMMA_SQ_ad * HINC / 3.0_jprb * 4.0_jprb
          S2_ad       = S2_ad + GAMMA_SQ_ad * HINC / 3.0_jprb * 2.0_jprb
          DO I = 1, NP
            JJ = 2 * I
            ff_ad(JJ) = ff_ad(JJ) + S2_ad
          ENDDO
          S2_ad = 0._jprb
          DO I = 1, NP
            JJ = 2 * I - 1
            ff_ad(JJ) = ff_ad(JJ) + S4_ad
          ENDDO
          S4_ad = 0._jprb
          ff_ad(1)              = ff_ad(1) + SA_ad
          ff_ad(coef%ws_nomega) = ff_ad(coef%ws_nomega) + SB_ad
          SA_ad = 0._jprb
          SB_ad = 0._jprb
!-----------Compute the frequency spectrum of the surface wave------------------
          alfa1_r = 1._jprb / sunglint%s(j)%alfa1
          omega_m_3 = sunglint%s(j)%omega_m ** 3_jpim
          DO I = 2, coef%ws_nomega
            omega_1 = omega_1_array(i)

            psi_ad(i) = psi_ad(i) + FF_ad(i) * coef%ws_k_omega(i) ** 2_jpim

            IF (sunglint%s(j)%omega_m < max_exp_exponent_4rt * omega_1) THEN ! stop exp overflow

              IF (sunglint%beta(i, j) > min_exponent) beta_ad = beta_ad + psi_ad(i) * sunglint%psi(i, j) * logk

              omega_m_ad = omega_m_ad - &
                  psi_ad(i) * sunglint%psi(i, j) * 5._jprb * omega_m_3 * omega_1_4_r(i)
              alfa1_ad   = alfa1_ad + psi_ad(i) * sunglint%psi(i, j) * alfa1_r

              IF (sunglint%beta(i, j) >= exp_max_exp_exponent_r) THEN
                IF (omega_1 <= sunglint%s(j)%omega_m) THEN
                  sigma2 = sigma_a ** 2_jpim
                ELSE
                  sigma2 = sigma_b ** 2_jpim
                ENDIF

                omega_m_ad = omega_m_ad + beta_ad * sunglint%beta(i, j) * omega_1 * &
                    (omega_1 - sunglint%s(j)%omega_m) / &
                    (sigma2 * omega_m_3)
              ENDIF
            ENDIF
            beta_ad = 0._jprb
          ENDDO
!           FF_ad(:) = 0._jprb
!           psi_ad(:) = 0._jprb

!-----------Set up some parameters used to define the JONSWAP frequency -------
!           spectrum of the surface wave                                       |
!------------------------------------------------------------------------------
          WINDSP_ad                = WINDSP_ad - omega_m_ad * sunglint%s(j)%omega_m / sunglint%s(j)%windsp
          x_u_ad                   = x_u_ad - omega_m_ad * 0.33_jprb * sunglint%s(j)%omega_m / sunglint%s(j)%x_u
          omega_m_ad               = 0._jprb
          x_u_ad                   = x_u_ad - 0.22_jprb * alfa1_ad * sunglint%s(j)%alfa1 / sunglint%s(j)%x_u
          alfa1_ad                 = 0._jprb
          WINDSP_ad                =      &
            & WINDSP_ad - X_U_ad * 2._jprb * profiles(j)%s2m%wfetc * gravity / sunglint%s(j)%windsp ** 3_jpim
          profiles_ad(j)%s2m%wfetc = profiles_ad(j)%s2m%wfetc + X_U_ad * gravity / sunglint%s(j)%windsp ** 2_jpim
          X_U_ad                   = 0._jprb
!-----------Compute the secant of the sun zenith angle at the surface-----------
          IF (raytracing%pathsun(ilevsur - 1, j) .NE. 1._jprb) THEN
            raytracing_ad%pathsun(ilevsur - 1, j) = raytracing_ad%pathsun(ilevsur - 1, j) +      &
              & ZENSUN_ad * (1._jprb / raytracing%pathsun(ilevsur - 1, j)) *                     &
              & (1._jprb / SQRT(raytracing%pathsun(ilevsur - 1, j) ** 2_jpim - 1._jprb))
            ZENSUN_ad = 0._jprb
          ELSE
            ZENSUN_ad = 0._jprb
          ENDIF
!-----------Compute the secant of the satellite viewing angle at the surface----
          IF (raytracing%pathsat(ilevsur - 1, j) .NE. 1._jprb) THEN
            raytracing_ad%pathsat(ilevsur - 1, j) = raytracing_ad%pathsat(ilevsur - 1, j) +      &
              & ZENSAT_ad * (1._jprb / raytracing%pathsat(ilevsur - 1, j)) *                     &
              & (1._jprb / SQRT(raytracing%pathsat(ilevsur - 1, j) ** 2_jpim - 1._jprb))
            ZENSAT_ad = 0._jprb
          ELSE
            ZENSAT_ad = 0._jprb
          ENDIF
!-----------Compute the angle between the wind direction and the U axis-------
          IF (sunglint%s(j)%windsp > 0._jprb) THEN
            profiles_ad(j)%s2m%v = profiles_ad(j)%s2m%v +      &
              & (1._jprb / (profiles(j)%s2m%u ** 2_jpim + profiles(j)%s2m%v ** 2_jpim)) * profiles(j)%s2m%u * WANGL_ad
            profiles_ad(j)%s2m%u = profiles_ad(j)%s2m%u -      &
              & (1._jprb / (profiles(j)%s2m%u ** 2_jpim + profiles(j)%s2m%v ** 2_jpim)) * profiles(j)%s2m%v * WANGL_ad
          ELSE
            WANGL_ad = 0._jprb
          ENDIF
!-----------Compute the wind speed----------------------------------------------
          IF (sunglint%s(j)%windsp > 0._jprb) THEN
            profiles_ad(j)%s2m%u = profiles_ad(j)%s2m%u + WINDSP_ad * profiles(j)%s2m%u / sunglint%s(j)%windsp
            profiles_ad(j)%s2m%v = profiles_ad(j)%s2m%v + WINDSP_ad * profiles(j)%s2m%v / sunglint%s(j)%windsp
            WINDSP_ad            = 0._jprb
          ELSE
            profiles_ad(j)%s2m%u = profiles_ad(j)%s2m%u + WINDSP_ad / SQRT(2._jprb)
            profiles_ad(j)%s2m%v = profiles_ad(j)%s2m%v + WINDSP_ad / SQRT(2._jprb)
            WINDSP_ad            = 0._jprb
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDDO
  sunglint_ad%s(:)%GLINT = 0._jprb
  sunglint_ad%s(:)%OMEGA = 0._jprb
  IF (LHOOK) CALL DR_HOOK('RTTOV_REFSUN_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE RTTOV_REFSUN_AD
