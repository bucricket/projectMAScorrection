!     Compute the fraction of solar radiance that is reflected by a wind
!     roughened water surface.
SUBROUTINE rttov_refsun_k( &
            & chanprof,     &
            & profiles,     &
            & profiles_k,   &
            & coef,         &
            & aux,          &
            & sunglint,     &
            & sunglint_k,   &
            & raytracing,   &
            & raytracing_k)
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
!     1            15/07/2003  Oroginal code. RTIASI-4. MarcoMatricardi. ECMWF.
!     2            20/02/2007  Removed polarisation  R Saunders
!     3            26/3/2007   Added R*8 definition for all reals R Saunders
!     4            15/09/2009  Determination of profile level index for the
!                              surface view angles simplified (P. Rayer)
!     5            02/12/2009  Pathsat, Pathsun and related quantities are now
!                              layer arrays (Marco Matricardi)
!     6            05/07/2010  Remove addsolar flag from profiles structure (J Hocking)
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
       & rttov_chanprof,  &
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
!       Scalar arguments with intent in:
  TYPE(rttov_chanprof ), INTENT(IN)    :: chanprof  (:)
!       Array arguments with intent in:
  TYPE(profile_type   ), INTENT(IN)    :: profiles  (:)
  TYPE(profile_type   ), INTENT(INOUT) :: profiles_k(size(chanprof))
  TYPE(profile_aux    ), INTENT(IN)    :: aux
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing_k
  TYPE(sunglint_type  ), INTENT(IN)    :: sunglint
  TYPE(sunglint_type  ), INTENT(INOUT) :: sunglint_k
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
!INTF_END
!     End of subroutine arguments
!       Local arrays
  REAL   (KIND=jprb) :: psi_k(coef%ws_nomega)                            ! The frequency spectrum of the surface wave
  REAL   (KIND=jprb) :: ff_k (coef%ws_nomega)                            ! Working space
!       Local scalars:
  INTEGER(KIND=jpim) :: j, jj, i, np, ilevsur, ilevwrt, prof
  REAL   (KIND=jprb) :: csi_k                                            ! Angle between the zenith and the proJection of the
! incident ray on the x-z plane
  REAL   (KIND=jprb) :: alfa_k                                           ! Angle between the incident ray and the x-z plane
  REAL   (KIND=jprb) :: c_shad_k                                         ! The average magnitude of tan(alfa)
  REAL   (KIND=jprb) :: p_prime_k                                        ! The probability density of tan(alfa)
  REAL   (KIND=jprb) :: pxy_gammaxy_k                                    ! The Joint  probability density of the along-view and
! cross view slope
  REAL   (KIND=jprb) :: rqs_k                                            ! Effective distribution function
  REAL   (KIND=jprb) :: gamma_sq_k                                       ! Total variance of the slope of the facet
  REAL   (KIND=jprb) :: gamma_o_k                                        ! The mean square of the along-view (X axis) slope
  REAL   (KIND=jprb) :: gamma_p_k                                        ! The mean square of the cross-view (Y axis) slope
  REAL   (KIND=jprb) :: g_shad_k                                         ! Normalization function
  REAL   (KIND=jprb) :: gammax_k                                         ! The x-slope of the facet
  REAL   (KIND=jprb) :: theta_fi                                         ! First order shadowing factor for the reflected ray
  REAL   (KIND=jprb) :: theta_csi                                        ! First ordesr shadowing factor for the incident ray
  REAL   (KIND=jprb) :: q_shad_k                                         ! Second order shadowing factor
  REAL   (KIND=jprb) :: q_shad_a_k                                       ! Second order shadowing factor
  REAL   (KIND=jprb) :: q_shad_b_k                                       ! Second order shadowing factor
  REAL   (KIND=jprb) :: zensat_k                                         ! Zenith angle of satellite viewing angle at surface
  REAL   (KIND=jprb) :: zensun_k                                         ! Zenith angle of sun at surface
  REAL   (KIND=jprb) :: omega_1                                          ! Frequency of the surface wave
  REAL   (KIND=jprb) :: fac1_k
  REAL   (KIND=jprb) :: a_shad_k
  REAL   (KIND=jprb) :: b_shad_k
  REAL   (KIND=jprb) :: lambda_a_k
  REAL   (KIND=jprb) :: lambda_b_k
  REAL   (KIND=jprb) :: x_u_k
  REAL   (KIND=jprb) :: alfa1_k
  REAL   (KIND=jprb) :: k
  REAL   (KIND=jprb) :: sigma_a
  REAL   (KIND=jprb) :: sigma_b
  REAL   (KIND=jprb) :: omega_m_k
  REAL   (KIND=jprb) :: sigma2
  REAL   (KIND=jprb) :: beta_k
  REAL   (KIND=jprb) :: windsp_k
  REAL   (KIND=jprb) :: wangl_k
  REAL   (KIND=jprb) :: s2_k, s4_k, sa_k, sb_k
  REAL   (KIND=jprb) :: bb, aa, hinc, fac
  INTEGER(KIND=jpim) :: nchannels

  REAL   (KIND=jprb) :: max_exp_exponent_4rt, alfa1_r, logk
  REAL   (KIND=jprb) :: exp_max_exp_exponent_r, omega_m_3
  REAL   (KIND=jprb) :: omega_1_array(coef%ws_nomega)
  REAL   (KIND=jprb) :: omega_1_4_r(coef%ws_nomega)

  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_REFSUN_K', 0_jpim, ZHOOK_HANDLE)
  nchannels = SIZE(chanprof)
  fac = SQRT(2 * pi)
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

  DO J = 1, nchannels
    prof          = chanprof(j)%prof
    IF (profiles(prof)%sunzenangle >= 0.0 .AND. &
        profiles(prof)%sunzenangle < max_sol_zen) THEN
      IF (profiles(prof)%skin%surftype == surftype_sea) THEN
        windsp_k      = 0._jprb
        wangl_k       = 0._jprb
        OMEGA_M_k     = 0._jprb
        RQS_k         = 0._jprb
        Q_shad_k      = 0._jprb
        Q_shad_A_k    = 0._jprb
        Q_shad_B_k    = 0._jprb
        g_shad_k      = 0._jprb
        p_prime_k     = 0._jprb
        Pxy_gammaxy_k = 0._jprb
        fac1_k        = 0._jprb
        gamma_O_k     = 0._jprb
        gammax_k      = 0._jprb
        alfa_k        = 0._jprb
        csi_k         = 0._jprb
        zensat_k      = 0._jprb
        zensun_k      = 0._jprb
        c_shad_k      = 0._jprb
        gamma_p_k     = 0._jprb
        gamma_sq_k    = 0._jprb
        lambda_a_k    = 0._jprb
        lambda_b_k    = 0._jprb
        b_shad_k      = 0._jprb
        a_shad_k      = 0._jprb
        sa_k          = 0._jprb
        sb_k          = 0._jprb
        s4_k          = 0._jprb
        s2_k          = 0._jprb
        ff_k(:)       = 0._jprb
        beta_k        = 0._jprb
        alfa1_k       = 0._jprb
        psi_k(:)      = 0._jprb
        x_u_k         = 0._jprb

!-----------Find the index of the pressure level that is closest to the surface-
        ilevwrt = aux%s(prof)%nearestlev_surf
!------------------
!redundant sections
!            if (aux(prof) %pfraction_surf < 0._jprb) then
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
        IF (sunglint%s(prof)%windsp == 0._jprb) THEN
          IF (sunglint%s(prof)%zensat == sunglint%s(prof)%zensun .AND. sunglint%s(prof)%dazng == 180.0_jprb) THEN
            sunglint_k%s(j)%glint = 0._jprb
          ELSE
            sunglint_k%s(j)%glint = 0._jprb
          ENDIF
        ELSE
          RQS_k         = RQS_k + sunglint_k%s(j)%glint
!-----------Compute the probability density that the slope of the facet at ----
!           a certain point is gammax when the incident ray and the ray from   |
!           the reflection of the incident ray on the local surface do not     |
!           intersect with any other surface.                                  |
!------------------------------------------------------------------------------
          Q_shad_k      = Q_shad_k +                                                                       &
            & RQS_k * sunglint%s(prof)%g_shad * sunglint%s(prof)%p_prime * sunglint%s(prof)%Pxy_gammaxy /  &
            & sunglint%s(prof)%fac1
          g_shad_k      = g_shad_k +                                                                       &
            & RQS_k * sunglint%s(prof)%Q_shad * sunglint%s(prof)%p_prime * sunglint%s(prof)%Pxy_gammaxy /  &
            & sunglint%s(prof)%fac1
          p_prime_k     = p_prime_k +                                                                     &
            & RQS_k * sunglint%s(prof)%Q_shad * sunglint%s(prof)%g_shad * sunglint%s(prof)%Pxy_gammaxy /  &
            & sunglint%s(prof)%fac1
          Pxy_gammaxy_k = Pxy_gammaxy_k +      &
            & RQS_k * sunglint%s(prof)%Q_shad * sunglint%s(prof)%g_shad * sunglint%s(prof)%p_prime / sunglint%s(prof)%fac1
          fac1_k        = fac1_k - RQS_k * sunglint%s(prof)%Q_shad * sunglint%s(prof)%g_shad * sunglint%s(prof)%p_prime     &
            &  * sunglint%s(prof)%Pxy_gammaxy / sunglint%s(prof)%fac1 ** 2_jpim
          RQS_k         = 0._jprb

          gamma_O_k     = gamma_O_k + Pxy_gammaxy_k * sunglint%s(prof)%Pxy_gammaxy *      &
            & (sunglint%s(prof)%gammax ** 2_jpim / sunglint%s(prof)%gamma_O ** 3_jpim - 1._jprb / sunglint%s(prof)%gamma_O)
          gammax_k      = gammax_k -      &
            & Pxy_gammaxy_k * sunglint%s(prof)%Pxy_gammaxy * sunglint%s(prof)%gammax / sunglint%s(prof)%gamma_O ** 2_jpim
          Pxy_gammaxy_k = 0._jprb
          alfa_k        = alfa_k - fac1_k * 2._jprb * sunglint%s(prof)%fac1 * tan(sunglint%s(prof)%alfa)
          csi_k         = csi_k -      &
            & fac1_k * sunglint%s(prof)%fac1 * tan((sunglint%s(prof)%csi - sunglint%s(prof)%zensat) / 2._jprb)
          zensat_k      = zensat_k +      &
            & fac1_k * sunglint%s(prof)%fac1 * tan((sunglint%s(prof)%csi - sunglint%s(prof)%zensat) / 2._jprb)
          fac1_k        = 0._jprb

          csi_k         = csi_k - g_shad_k * 0.5_jprb * tan(sunglint%s(prof)%zensat) /      &
            & cos(0.5_jprb * (sunglint%s(prof)%csi - sunglint%s(prof)%zensat)) ** 2_jpim
          zensat_k      = zensat_k + 0.5_jprb * g_shad_k * tan(sunglint%s(prof)%zensat) / &
            & cos(0.5_jprb * (sunglint%s(prof)%csi - sunglint%s(prof)%zensat)) ** 2_jpim
          zensat_k      = zensat_k - g_shad_k * tan(0.5_jprb * (sunglint%s(prof)%csi - sunglint%s(prof)%ZENSAT)) / &
            & cos(sunglint%s(prof)%zensat) ** 2_jpim
          g_shad_k      = 0._jprb

          alfa_k        = alfa_k - p_prime_k * sunglint%s(prof)%p_prime * tan(sunglint%s(prof)%alfa) /      &
            & (cos(sunglint%s(prof)%alfa) * sunglint%s(prof)%c_shad) ** 2_jpim
          c_shad_k      = c_shad_k + p_prime_k * sunglint%s(prof)%p_prime *      &
            & ((tan(sunglint%s(prof)%alfa)) ** 2_jpim / sunglint%s(prof)%c_shad ** 3_jpim - &
            & 1._jprb / sunglint%s(prof)%c_shad)
          p_prime_k     = 0._jprb

          IF ((cos(sunglint%s(prof)%csi) + cos(sunglint%s(prof)%zensat)) >= 0._jprb) THEN
            zensat_k  = zensat_k - c_shad_k * sunglint%s(prof)%gamma_P * sin(sunglint%s(prof)%zensat)
            csi_k     = csi_k - c_shad_k * sunglint%s(prof)%gamma_P * sin(sunglint%s(prof)%csi)
            gamma_P_k = gamma_P_k + c_shad_k * (cos(sunglint%s(prof)%csi) + cos(sunglint%s(prof)%zensat))
            c_shad_k  = 0._jprb
          ELSE
            zensat_k  = zensat_k + c_shad_k * sunglint%s(prof)%gamma_P * sin(sunglint%s(prof)%zensat)
            csi_k     = csi_k + c_shad_k * sunglint%s(prof)%gamma_P * sin(sunglint%s(prof)%csi)
            gamma_P_k = gamma_P_k - c_shad_k * (cos(sunglint%s(prof)%csi) + cos(sunglint%s(prof)%zensat))
            c_shad_k  = 0._jprb
          ENDIF
!-----------Compute the value of the function that represents the shadowing  --
!           of the incident and reflected ray                                  |
!------------------------------------------------------------------------------
!-----------First order shadowing (the slope of the facet is negative)---------
          IF (abs(sunglint%s(prof)%zensat) > 0._jprb) THEN
            IF ((1._jprb / tan(abs(sunglint%s(prof)%zensat)) - &
                sunglint%s(prof)%gammax * sunglint%s(prof)%zensat / abs(sunglint%s(prof)%zensat)) >= 0._jprb) THEN
              theta_fi = 1._jprb
            ELSE
              theta_fi = 0._jprb
            ENDIF
          ELSE
            theta_fi = 1._jprb
          ENDIF
          IF (abs(sunglint%s(prof)%csi) > 0._jprb) THEN
            IF ((1._jprb / tan(abs(sunglint%s(prof)%csi)) + &
              & sunglint%s(prof)%gammax * sunglint%s(prof)%csi / abs(sunglint%s(prof)%csi)) >= 0._jprb) THEN
              theta_csi = 1._jprb
            ELSE
              theta_csi = 0._jprb
            ENDIF
          ELSE
            theta_csi = 1._jprb
          ENDIF
          IF (abs(sunglint%s(prof)%csi) > PI / 2._jprb) theta_csi = 0._jprb
!-----------Second order shadowing (the facet cannot be seen)-------------------
          IF ((sunglint%s(prof)%zensat * sunglint%s(prof)%csi) <= 0._jprb) THEN
            IF (sunglint%s(prof)%zensat == 0._jprb .AND. sunglint%s(prof)%csi == 0._jprb) THEN
              Q_shad_k = 0._jprb
            ELSE IF (sunglint%s(prof)%zensat >= abs(sunglint%s(prof)%csi)) THEN
              Q_shad_A_k = Q_shad_A_k + Q_shad_k * theta_fi
              Q_shad_k   = 0._jprb
              lambda_a_k = lambda_a_k - Q_shad_A_k / (1._jprb + sunglint%s(prof)%lambda_a) ** 2
              Q_shad_A_k = 0._jprb
              a_shad_k   = a_shad_k - lambda_a_k * (exp( - sunglint%s(prof)%a_shad ** 2 / 2._jprb)) /      &
                & (fac * sunglint%s(prof)%a_shad ** 2)
              lambda_a_k = 0._jprb
              zensat_k   = zensat_k -      &
                & a_shad_k * (1._jprb / (sin(sunglint%s(prof)%zensat)) ** 2) * (1._jprb / sunglint%s(prof)%gamma_O)
              gamma_O_k  = gamma_O_k -      &
                & a_shad_k * (1._jprb / tan(sunglint%s(prof)%zensat)) * (1._jprb / sunglint%s(prof)%gamma_O ** 2)
              a_shad_k   = 0._jprb
            ELSE IF (abs(sunglint%s(prof)%csi) > sunglint%s(prof)%zensat) THEN
              Q_shad_B_k = Q_shad_B_k + Q_shad_k * theta_csi
              Q_shad_k   = 0._jprb
              lambda_b_k = lambda_b_k - Q_shad_B_k / (1._jprb + sunglint%s(prof)%lambda_b) ** 2
              Q_shad_B_k = 0._jprb
              b_shad_k   = b_shad_k - lambda_b_k * (exp( - sunglint%s(prof)%b_shad ** 2 / 2._jprb)) /      &
                & (fac * sunglint%s(prof)%b_shad ** 2)
              lambda_b_k = 0._jprb
              IF (sunglint%s(prof)%csi >= 0._jprb) THEN
                csi_k     =      &
                  & csi_k + b_shad_k * ( - 1._jprb / (sin(sunglint%s(prof)%csi)) ** 2 * 1._jprb / sunglint%s(prof)%gamma_O)
                gamma_O_k =      &
                  & gamma_O_k - b_shad_k * ((1._jprb / tan(sunglint%s(prof)%csi)) / sunglint%s(prof)%gamma_O ** 2)
                b_shad_k  = 0._jprb
              ELSE
                csi_k     = csi_k +      &
                  & b_shad_k * (1._jprb / (sin(sunglint%s(prof)%csi)) ** 2 * 1._jprb / sunglint%s(prof)%gamma_O)
                gamma_O_k =      &
                  & gamma_O_k + b_shad_k * ((1._jprb / tan(sunglint%s(prof)%csi)) / sunglint%s(prof)%gamma_O ** 2)
                b_shad_k  = 0._jprb
              ENDIF
            ENDIF
          ELSE
            Q_shad_k   = Q_shad_k * theta_fi * theta_csi
            lambda_a_k = lambda_a_k - Q_shad_k / (1._jprb + sunglint%s(prof)%lambda_a + sunglint%s(prof)%lambda_b) ** 2
            lambda_b_k = lambda_b_k - Q_shad_k / (1._jprb + sunglint%s(prof)%lambda_a + sunglint%s(prof)%lambda_b) ** 2
            Q_shad_k   = 0._jprb
            b_shad_k   = b_shad_k - lambda_b_k * (exp( - sunglint%s(prof)%b_shad ** 2 / 2._jprb)) /      &
              & (fac * sunglint%s(prof)%b_shad ** 2)
            lambda_b_k = 0._jprb
            a_shad_k   = a_shad_k - lambda_a_k * (exp( - sunglint%s(prof)%a_shad ** 2 / 2._jprb)) /      &
              & (fac * sunglint%s(prof)%a_shad ** 2)
            lambda_a_k = 0._jprb
            IF (sunglint%s(prof)%csi >= 0._jprb) THEN
              csi_k     =      &
                & csi_k + b_shad_k * ( - 1._jprb / (sin(sunglint%s(prof)%csi)) ** 2 * 1._jprb / sunglint%s(prof)%gamma_O)
              gamma_O_k = gamma_O_k - b_shad_k * ((1._jprb / tan(sunglint%s(prof)%csi)) / sunglint%s(prof)%gamma_O ** 2)
              b_shad_k  = 0._jprb
            ELSE IF (sunglint%s(prof)%csi < 0._jprb) THEN
              csi_k     =      &
                & csi_k + b_shad_k * (1._jprb / (sin(sunglint%s(prof)%csi)) ** 2 * 1._jprb / sunglint%s(prof)%gamma_O)
              gamma_O_k =      &
                & gamma_O_k + b_shad_k * ((1._jprb / tan(sunglint%s(prof)%csi)) / sunglint%s(prof)%gamma_O ** 2)
              b_shad_k  = 0._jprb
            ENDIF
            ZENSAT_k  =      &
              & ZENSAT_k - a_shad_k * (1._jprb / (sin(sunglint%s(prof)%zensat)) ** 2) * (1._jprb / sunglint%s(prof)%gamma_O)
            gamma_O_k =      &
              & gamma_O_k - a_shad_k * (1._jprb / tan(sunglint%s(prof)%zensat)) * (1._jprb / sunglint%s(prof)%gamma_O ** 2)
            a_shad_k  = 0._jprb
          ENDIF
!-----------First order shadowing (the slope of the facet is negative)----------
          csi_k    = csi_k +      &
            & gammax_k * (1._jprb / (cos((sunglint%s(prof)%csi - sunglint%s(prof)%zensat) / 2._jprb)) ** 2) / 2._jprb
          zensat_k = zensat_k -      &
            & gammax_k * (1._jprb / (cos((sunglint%s(prof)%csi - sunglint%s(prof)%zensat) / 2._jprb)) ** 2) / 2._jprb
          gammax_k = 0._jprb
!-----------Obtain angles csi and alfa------------------------------------------
          IF ((sunglint%s(prof)%csi + sunglint%s(prof)%zensat) >= 0._jprb) THEN
            csi_k    = csi_k + sunglint_k%s(j)%omega / 2._jprb
            zensat_k = zensat_k + sunglint_k%s(j)%omega / 2._jprb
          ELSE
            csi_k    = csi_k - sunglint_k%s(j)%omega / 2._jprb
            zensat_k = zensat_k - sunglint_k%s(j)%omega / 2._jprb
          ENDIF
          IF ((sin(sunglint%s(prof)%zensun) * sin(PI - sunglint%s(prof)%dazng * deg2rad)) ** 2 .NE. 1._jprb) THEN
            ZENSUN_k = ZENSUN_k + alfa_k * (cos(sunglint%s(prof)%zensun) * sin(PI - sunglint%s(prof)%dazng * deg2rad) *      &
              & (1._jprb / SQRT(1._jprb - (sin(sunglint%s(prof)%zensun) * sin(PI - sunglint%s(prof)%dazng * deg2rad)) ** 2)) &
              & )
            alfa_k   = 0._jprb
          ELSE
            alfa_k = 0._jprb
          ENDIF
          ZENSUN_k   = ZENSUN_k +                                                                                            &
            & csi_k * cos(PI - sunglint%s(prof)%dazng * deg2rad) * (1._jprb / (cos(sunglint%s(prof)%zensun)) ** 2) * 1._jprb &
            &  / (1._jprb + (tan(sunglint%s(prof)%zensun) * cos(PI - sunglint%s(prof)%dazng * deg2rad)) ** 2)
          csi_k      = 0._jprb
!-----------Compute the rms of the slope of the facet along the X (gamma_O)----
!           and Y (gamma_P) axis.                                              |
!------------------------------------------------------------------------------
          GAMMA_SQ_k = GAMMA_SQ_k +      &
            & gamma_P_k * (2._jprb - cos(2._jprb * sunglint%s(prof)%WANGL)) / (8._jprb * sunglint%s(prof)%gamma_P)
          WANGL_k    = WANGL_k + gamma_P_k * sin(2._jprb * sunglint%s(prof)%WANGL) * sunglint%s(prof)%GAMMA_SQ /      &
            & (4._jprb * sunglint%s(prof)%gamma_P)
          gamma_P_k  = 0._jprb
          GAMMA_SQ_k = GAMMA_SQ_k +      &
            & gamma_O_k * (2._jprb + cos(2._jprb * sunglint%s(prof)%WANGL)) / (8._jprb * sunglint%s(prof)%gamma_O)
          WANGL_k    = WANGL_k - gamma_O_k * sin(2._jprb * sunglint%s(prof)%WANGL) * sunglint%s(prof)%GAMMA_SQ /      &
            & (4._jprb * sunglint%s(prof)%gamma_O)
          gamma_O_k  = 0._jprb
!-----------Compute the Simpson's integral--------------------------------------
          BB = 301._jprb
          AA = 0._jprb
          NP = 1501_jpim
          HINC       = (BB - AA) / REAL(2 * NP, jprb)
          SA_k       = SA_k + GAMMA_SQ_k * HINC / 3.0_jprb
          SB_k       = SB_k + GAMMA_SQ_k * HINC / 3.0_jprb
          S4_k       = S4_k + GAMMA_SQ_k * HINC / 3._jprb * 4._jprb
          S2_k       = S2_k + GAMMA_SQ_k * HINC / 3._jprb * 2._jprb
          GAMMA_SQ_k = 0._jprb
          DO I = 1, NP
            jj       = 2 * I
            ff_k(jj) = ff_k(jj) + S2_k
          ENDDO
          S2_k = 0._jprb
          DO I = 1, NP
            jj       = 2 * I - 1
            ff_k(jj) = ff_k(jj) + S4_k
          ENDDO
          S4_k = 0._jprb
          ff_k(1)              = ff_k(1) + SA_k
          ff_k(coef%ws_nomega) = ff_k(coef%ws_nomega) + SB_k
          SA_k = 0._jprb
          SB_k = 0._jprb
!-----------Compute the frequency spectrum of the surface wave------------------
          alfa1_r = 1._jprb / sunglint%s(prof)%alfa1
          omega_m_3 = sunglint%s(prof)%omega_m ** 3_jpim
          DO I = 2, coef%ws_nomega
            omega_1 = omega_1_array(i)

            psi_k(i) = psi_k(i) + FF_k(i) * coef%ws_k_omega(i) ** 2_jpim

            IF (sunglint%s(prof)%omega_m < max_exp_exponent_4rt * omega_1) THEN ! stop exp overflow

              IF (sunglint%beta(i, prof) > min_exponent) beta_k = beta_k + psi_k(i) * sunglint%psi(i, prof) * logk

              omega_m_k = omega_m_k - &
                  psi_k(i) * sunglint%psi(i, prof) * 5._jprb * omega_m_3 * omega_1_4_r(i)
              alfa1_k   = alfa1_k + psi_k(i) * sunglint%psi(i, prof) * alfa1_r

              IF (sunglint%beta(i, prof) >= exp_max_exp_exponent_r) THEN
                IF (omega_1 <= sunglint%s(prof)%omega_m) THEN
                  sigma2 = sigma_a ** 2_jpim
                ELSE
                  sigma2 = sigma_b ** 2_jpim
                ENDIF

                omega_m_k = omega_m_k + beta_k * sunglint%beta(i, prof) * omega_1 * &
                    (omega_1 - sunglint%s(prof)%omega_m) / &
                    (sigma2 * omega_m_3)
              ENDIF
            ENDIF
            beta_k = 0._jprb
          ENDDO
!           FF_k(:) = 0._jprb
!           psi_k(:) = 0._jprb

!-----------Set up some parameters used to define the JONSWAP frequency -------
!           spectrum of the surface wave                                       |
!------------------------------------------------------------------------------
          WINDSP_k                = WINDSP_k - omega_m_k * sunglint%s(prof)%omega_m / sunglint%s(prof)%windsp
          x_u_k                   = x_u_k - omega_m_k * 0.33_jprb * sunglint%s(prof)%omega_m / sunglint%s(prof)%x_u
          omega_m_k               = 0._jprb
          X_U_k = X_U_k - 0.22_jprb * alfa1_k * sunglint%s(prof)%alfa1 / sunglint%s(prof)%x_u
          alfa1_k                 = 0._jprb
          WINDSP_k                =      &
            & WINDSP_k - X_U_k * 2._jprb * profiles(prof)%s2m%wfetc * gravity / sunglint%s(prof)%windsp ** 3
          profiles_k(j)%s2m%wfetc = profiles_k(j)%s2m%wfetc + X_U_k * gravity / sunglint%s(prof)%windsp ** 2
          X_U_k = 0._jprb
!-----------Compute the secant of the sun zenith angle at the surface-----------
          IF (raytracing%pathsun(ilevsur - 1, prof) .NE. 1._jprb) THEN
            raytracing_k%pathsun(ilevsur - 1, j) = raytracing_k%pathsun(ilevsur - 1, j) +      &
              & ZENSUN_k * (1._jprb / raytracing%pathsun(ilevsur - 1, prof)) *                 &
              & (1._jprb / SQRT(raytracing%pathsun(ilevsur - 1, prof) ** 2 - 1._jprb))
            ZENSUN_k = 0._jprb
          ELSE
            ZENSUN_k = 0._jprb
          ENDIF
!-----------Compute the secant of the satellite viewing angle at the surface----
          IF (raytracing%pathsat(ilevsur - 1, prof) .NE. 1._jprb) THEN
            raytracing_k%pathsat(ilevsur - 1, j) = raytracing_k%pathsat(ilevsur - 1, j) +      &
              & ZENSAT_k * (1._jprb / raytracing%pathsat(ilevsur - 1, prof)) *                 &
              & (1._jprb / SQRT(raytracing%pathsat(ilevsur - 1, prof) ** 2 - 1._jprb))
            ZENSAT_k = 0._jprb
          ELSE
            ZENSAT_k = 0._jprb
          ENDIF
!-----------Compute the angle between the wind direction and the U axis-------
          IF (sunglint%s(prof)%windsp > 0._jprb) THEN
            profiles_k(j)%s2m%v = profiles_k(j)%s2m%v +      &
              & (1._jprb / (profiles(prof)%s2m%u ** 2 + profiles(prof)%s2m%v ** 2)) * profiles(prof)%s2m%u * WANGL_k
            profiles_k(j)%s2m%u = profiles_k(j)%s2m%u -      &
              & (1._jprb / (profiles(prof)%s2m%u ** 2 + profiles(prof)%s2m%v ** 2)) * profiles(prof)%s2m%v * WANGL_k
          ELSE
            WANGL_k = 0._jprb
          ENDIF
!-----------Compute the wind speed----------------------------------------------
          IF (sunglint%s(prof)%windsp > 0._jprb) THEN
            profiles_k(j)%s2m%u = profiles_k(j)%s2m%u + WINDSP_k * profiles(prof)%s2m%u / sunglint%s(prof)%windsp
            profiles_k(j)%s2m%v = profiles_k(j)%s2m%v + WINDSP_k * profiles(prof)%s2m%v / sunglint%s(prof)%windsp
            WINDSP_k            = 0._jprb
          ELSE
            profiles_k(j)%s2m%u = profiles_k(j)%s2m%u + WINDSP_k / SQRT(2._jprb)
            profiles_k(j)%s2m%v = profiles_k(j)%s2m%v + WINDSP_k / SQRT(2._jprb)
            WINDSP_k            = 0._jprb
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDDO
  sunglint_k%s(:)%GLINT = 0._jprb
  sunglint_k%s(:)%OMEGA = 0._jprb
  IF (LHOOK) CALL DR_HOOK('RTTOV_REFSUN_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_refsun_k
