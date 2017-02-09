!     Set up aerosols optical parameters for a climatological profile
SUBROUTINE rttov_cldstr_ad( &
            & opts_rt_ir,  &
            & profiles,    &
            & profiles_ad, &
            & ircld,       &
            & ircld_ad)
!     Description:
!     To set up profile-dependent variables for subsequent
!     rt calculations by other subroutines of RTIASI.
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
!
!     Owner:
!     EUMETSAT
!
!     History:
!     Version      Date        Comment
!     1           01/6/2004    Marco Matricardi. ECMWF.
!     2           24/03/2006   Marco Matricardi. ECMWF.
!                              A bug has been fixed tha caused
!                              the computation of a wrong number of
!                              streams in some situation.Also, a new feature
!                              has been introduced that allow to reduce the
!                              number of streams by considering only those
!                              streams whose weight is larger than
!                              cldstr_threshold. cldstr_threshold is pecified
!                              in the module rttov_const. By setting
!                              cldstr_threshold to a negative number,
!                              all the streams will be processed.
!                              This feature is used with caution.
!                              In fact, the sum of the weights of all streams
!                              (including the clear one) must be equal to 1.
!                              If some streams are not considered,
!                              this means that the weight of the clear stream
!                              has to be adjusted. As a consequence, if a too
!                              large value is used for cldstr_threshold,
!                              this can have serious implications for the
!                              accuracy of the results. In conclusion, only
!                              very small values should be used for
!                              cldstr_threshold just to remove the streams
!                              with a very small weight.
!     3           01/07/2006   Marco Matricardi. ECMWF.
!                              Rewritten for RTTOV
!     4           15/08/2009   User defined ToA. Layers distinct from levels (P.Rayer)
!     5           02/12/2009   All variables stored as layer arrays (Marco Matricardi)
!     6           03/12/2009   Fixed bug (Marco Matricardi)
!     7           07/06/2013   Single stream option to make memory/compute 
!                              requirements more reasonable (A Geer)
!     8           23/06/2014   ircld_tl%XSTR initialized to zero to avoid being undefined 
!                              in clear sky (S. Migliorini)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Module used:
  USE rttov_types, ONLY : profile_type, ircld_type, rttov_opts_rt_ir
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
  USE rttov_const, ONLY : cldstr_low_cloud_top
!INTF_ON
  IMPLICIT NONE
!       Arguments with intent in:
  TYPE(rttov_opts_rt_ir), INTENT(IN)    :: opts_rt_ir
  TYPE(profile_type),     INTENT(IN)    :: profiles(:)
!       Arguments with intent inout:
  TYPE(profile_type),     INTENT(INOUT) :: profiles_ad(SIZE(profiles))
  TYPE(ircld_type),       INTENT(INOUT) :: ircld
  TYPE(ircld_type),       INTENT(INOUT) :: ircld_ad
!INTF_END
!     End of subroutine arguments
!       Local scalars:
  INTEGER(KIND=jpim) :: i, j, istr, ijstr, ilay, ic
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: ibdy_layer(1), imax_cfrac(1)
  REAL   (KIND=jprb) :: cfrac_max, cloud_tot(profiles(1)%nlayers) 
!       Local tangent linear arrays:
  REAL   (KIND=jprb) :: cfrac_max_ad 
  REAL   (KIND=jprb) :: ntot(profiles(1)%nlevels)
!       Local direct arrays arrays:
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
!---------Compute number of streams and cloud distribution in each stream-------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CLDSTR_AD', 0_jpim, ZHOOK_HANDLE)
  nprofiles = SIZE(profiles)

  IF (opts_rt_ir%cldstr_simple) THEN

    ! A simpler, faster, but quite approximate "Cmax" single-stream approach. 
    ! Main benefit is that it is much more memory efficient. Intended mainly
    ! for mid- and upper-tropospheric channels.
    DO J = 1, NPROFILES

      ! Find maximum cloud fraction and cloud amount above the boundary layer
      ! NB pressure is on levels, cloud in layers, so p(1:nlevels-1) is the top of the layer
      ibdy_layer = MINLOC(ABS(profiles(j)%p(1:profiles(j)%nlevels-1) - cldstr_low_cloud_top))
      imax_cfrac = MAXLOC(profiles(j)%cfrac(1:ibdy_layer(1)))
      cfrac_max  = profiles(j)%cfrac(imax_cfrac(1))
      cloud_tot(:) = SUM(profiles(j)%cloud(:,:), DIM=1)

      cfrac_max_ad = 0.0_jprb

      ! Ignoring trivial amounts of cloud and precip
      IF (ANY(cloud_tot(1:ibdy_layer(1)) > 1E-6_jprb) .AND. cfrac_max > 1E-3) THEN

        ! Cloudy stream required
        cfrac_max_ad = cfrac_max_ad + ircld_ad%XSTR(2, J)
        cfrac_max_ad = cfrac_max_ad - ircld_ad%XSTRCLR(J)

      ENDIF

      profiles_ad(j)%cfrac(imax_cfrac(1)) = profiles_ad(j)%cfrac(imax_cfrac(1)) + cfrac_max_ad

    ENDDO

  ELSE

    DO J = 1, nprofiles
      loop2 : DO istr = 2, ircld%icount(J) - 1
        DO i = istr - 1, 1,  - 1
          ircld%xstrref1(istr, I, J) = ircld%xstrref(i, J)
          IF (ircld%xstrref(I, J) <= ircld%a(istr, J)) THEN
            ircld%xstrref(I + 1, J) = ircld%a(istr, J)
            CYCLE loop2
          ELSE
            ircld%xstrref(i + 1, J) = ircld%xstrref(i, J)
          ENDIF
        ENDDO
        i = 0
        ircld%xstrref(I + 1, J) = ircld%a(istr, J)
      ENDDO loop2
      ircld_ad%A(:, J)       = 0._jprb
      ircld_ad%MAXCOV(:, J)  = 0._jprb
      ircld_ad%CLDCFR(:, J)  = 0._jprb
      ircld_ad%XSTRMIN(:, J) = 0._jprb
      ircld_ad%XSTRMAX(:, J) = 0._jprb
      NTOT(:)                = 0._jprb
  !---------Consider only the streams whose weight is greater than cldstr_threshold-------
      IF (ircld%nstreamref(J) /= 0_jpim) THEN
        IF (ircld%icounstr(J) /= 0_jpim) THEN
          ircld_ad%xstr(ircld%icounstr(J) + 1, J) = ircld_ad%xstr(ircld%icounstr(J) + 1, J) - ircld_ad%xstrclr(j)
          ircld_ad%xstr(1, J)                     = ircld_ad%xstr(1, J) + ircld_ad%xstrclr(j)
          ircld_ad%xstrclr(J)                     = 0._jprb
          DO istr = ircld%icounstr(j) + 1, 1,  - 1
            IF (istr == 1) THEN
              IF (ircld%indexstr(istr, J) /= istr) THEN
                ircld_ad%xstr(ircld%indexstr(istr, J), J) =      &
                  & ircld_ad%xstr(ircld%indexstr(istr, J), J) + ircld_ad%xstr(istr, J)
                ircld_ad%xstr(istr, J)                    = 0._jprb
              ENDIF
            ELSE
              ircld_ad%xstr(istr - 1, J) = ircld_ad%xstr(istr - 1, j) + ircld_ad%xstr(istr, j)
              IF ((ircld%indexstr(istr - 1, J) + 1) /= istr) THEN
                ircld_ad%xstr(ircld%indexstr(istr - 1, J) + 1, J) =      &
                  & ircld_ad%xstr(ircld%indexstr(istr - 1, J) + 1, J) + ircld_ad%xstr(istr, j)
              ELSE
                ircld_ad%xstr(ircld%indexstr(istr - 1, J) + 1, J) = ircld_ad%xstr(istr, J)
              ENDIF
              IF (ircld%indexstr(istr - 1, J) /= istr) THEN
                ircld_ad%xstr(ircld%indexstr(istr - 1, j), j) =      &
                  & ircld_ad%xstr(ircld%indexstr(istr - 1, j), j) - ircld_ad%xstr(istr, j)
              ELSE
                ircld_ad%xstr(ircld%indexstr(istr - 1, j), j) =  - ircld_ad%xstr(istr, j)
              ENDIF
              IF (((ircld%indexstr(istr - 1, j) + 1) /= istr) .AND. (ircld%indexstr(istr - 1, J) /= istr)) THEN
                ircld_ad%xstr(istr, j) = 0._jprb
              ENDIF
            ENDIF
          ENDDO
        ELSE
          ircld_ad%XSTRCLR(J) = 0.0_JPRB
        ENDIF
      ENDIF
  !---------Compute the weight of the clear stream--------------------------------
      ircld_ad%xstr(ircld%nstreamref(j) + 1, j) = ircld_ad%xstr(ircld%nstreamref(j) + 1, j) - ircld_ad%xstrclr(j)
      ircld_ad%xstr(1, j)                       = ircld_ad%xstr(1, j) + ircld_ad%xstrclr(j)
  !---------Re-arrange the limits of each stream in ascending order---------------
      outer : DO i = ircld%iloop(j), 1,  - 1
        inner : DO istr = ircld%iloopin(i, j), 1,  - 1
          IF (ircld%xstrref2(i, istr, j) == ircld%xstrref2(i, istr + 1, j)) THEN
            IF (ircld%xstrref2(i, istr, j) /= 1.) THEN
              ircld%icount1ref(i, j) = ircld%icount1ref(i, j) - 1
              DO ijstr = ircld%icount1ref(i, j), istr,  - 1
                ircld_ad%xstr(ijstr + 1, j) = ircld_ad%xstr(ijstr + 1, j) + ircld_ad%xstr(ijstr, j)
                ircld_ad%xstr(ijstr, J)     = 0._jprb
              ENDDO
            ENDIF
          ENDIF
        ENDDO inner
      ENDDO outer
      loop1 : DO istr = ircld%icount(j) - 1, 2,  - 1
        IF (.NOT. ircld%flag(istr, J)) THEN
          i = 0
          ircld_ad%a(istr, J)     = ircld_ad%a(istr, j) + ircld_ad%xstr(i + 1, j)
          ircld_ad%xstr(i + 1, j) = 0._jprb                                      !
        ENDIF
        DO i = ircld%iflag(istr, j), istr - 1
          IF (ircld%xstrref1(istr, i, j) <= ircld%a(istr, j)) THEN
            ircld_ad%a(istr, j)     = ircld_ad%a(istr, j) + ircld_ad%xstr(i + 1, j)
            ircld_ad%xstr(i + 1, j) = 0._jprb
          ELSE
            ircld_ad%xstr(i, j)     = ircld_ad%xstr(i, j) + ircld_ad%xstr(i + 1, j)
            ircld_ad%xstr(i + 1, j) = 0._jprb
          ENDIF
        ENDDO
        ircld_ad%xstr(istr, j) = ircld_ad%xstr(istr, j) + ircld_ad%a(istr, j)
        ircld_ad%a(istr, j)    = 0._jprb
      ENDDO loop1
  !---------Determine the limits of each stream----------------------------------
      ic = ircld%icount(j)
      DO ilay = profiles(1)%nlayers, 1,  - 1
        IF (ircld%xstrmax(ilay, j) /= 0._jprb) THEN
          ic = ic - 1
          ircld_ad%xstrmax(ilay, j) = ircld_ad%xstrmax(ilay, j) + ircld_ad%xstr(ic, j)
          IC = IC - 1
          ircld_ad%xstrmin(ilay, j) = ircld_ad%xstrmin(ilay, j) + ircld_ad%xstr(ic, j)
        ENDIF
        IF (ircld%xstrminref(ilay, j) < 0._jprb) THEN
          ircld_ad%xstrmin(ilay, j) = 0._jprb
        ENDIF
        ntot(ilay)               = ntot(ilay) + ircld_ad%xstrmin(ilay, j)
        ircld_ad%cldcfr(ilay, j) = ircld_ad%cldcfr(ilay, j) - ircld_ad%xstrmin(ilay, j)
        ntot(ilay)               = ntot(ilay) + ircld_ad%xstrmax(ilay, j)
      ENDDO
  !---------Compute the cumulative cloud coverage usin the maximum-random---------
  !         overlap assumption
      DO ilay = profiles(1)%nlayers, 1,  - 1
        IF (ircld%cldcfr(ilay, j) /= 0._jprb) THEN
          ntot(ilay) =  - ntot(ilay)
        ELSE
          ntot(ilay) = 0._jprb
        ENDIF
      ENDDO
      DO ilay = profiles(1)%nlayers, 2,  - 1
        IF (ircld%ntotref(ilay - 1, j) == 0._jprb .AND. (1. - ircld%cldcfr(ilay - 1, J)) == 0._jprb) THEN
          ircld_ad%maxcov(ilay, J) = ircld_ad%maxcov(ilay, j) + ntot(ilay)
        ELSE
          ntot(ilay - 1)               =      &
            & ntot(ilay - 1) + ntot(ilay) * ircld%maxcov(ilay, j) / (1. - ircld%cldcfr(ilay - 1, j))
          ircld_ad%cldcfr(ilay - 1, J) = ircld_ad%cldcfr(ilay - 1, J) +      &
            & ntot(ilay) * ircld%ntotref(ilay - 1, j) * ircld%maxcov(ilay, J) / (1._jprb - ircld%cldcfr(ilay - 1, J)) ** 2
          ircld_ad%maxcov(ilay, J)     =      &
            & ircld_ad%maxcov(ilay, j) + ntot(ilay) * ircld%ntotref(ilay - 1, j) / (1._jprb - ircld%cldcfr(ilay - 1, J))
        ENDIF
        IF ((ircld%cldcfr(ilay - 1, J) > ircld%cldcfr(ilay, J))) THEN
          ircld_ad%cldcfr(ilay - 1, J) = ircld_ad%cldcfr(ilay - 1, J) - ircld_ad%maxcov(ilay, J)
          ircld_ad%MAXCOV(ilay, J)     = 0._jprb
        ELSE IF ((ircld%cldcfr(ilay - 1, J) < ircld%cldcfr(ilay, J))) THEN
          ircld_ad%cldcfr(ilay, J) = ircld_ad%cldcfr(ILAY, J) - ircld_ad%maxcov(ilay, J)
          ircld_ad%maxcov(ilay, J) = 0._jprb
        ELSE IF ((ircld%cldcfr(ilay - 1, J) == ircld%cldcfr(ilay, J))) THEN
          ircld_ad%cldcfr(ilay, j) = ircld_ad%cldcfr(ilay, J) - ircld_ad%maxcov(ilay, J)
          ircld_ad%maxcov(ilay, J) = 0._jprb
        ENDIF
      ENDDO
      ircld_ad%cldcfr(1, J) = ircld_ad%cldcfr(1, J) - ntot(1)
  !---------Compute number of streams and cloud distribution in each stream-------
      DO ilay = profiles(1)%nlayers, 1,  - 1
        IF (profiles(j)%cfrac(ilay) > 0.) THEN
          profiles_ad(j)%cfrac(ilay) = profiles_ad(j)%cfrac(ilay) + ircld_ad%cldcfr(ilay, j)
        ENDIF
      ENDDO
      ircld_ad%CLDCFR(:, j)  = 0._jprb
      ircld_ad%XSTRMIN(:, j) = 0._jprb
      ircld_ad%XSTRMAX(:, j) = 0._jprb
      ircld_ad%XSTR(:, j) = 0._jprb 
      ntot = 0._jprb
    ENDDO

  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_CLDSTR_AD', 1_jpim, ZHOOK_HANDLE)
!      ircld_ad(:)%xstrclr=0._jprb
END SUBROUTINE rttov_cldstr_ad
