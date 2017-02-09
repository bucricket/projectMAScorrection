!     Set up aerosols optical parameters for a climatological profile
SUBROUTINE rttov_cldstr_k( &
            & opts_rt_ir, &
            & chanprof,   &
            & profiles,   &
            & profiles_k, &
            & ircld,      &
            & ircld_k)
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
!     4           01/03/2007   Removed polarisation R Saunders
!     5           15/09/2009   User defined ToA. Layers distinct from levels (P.Rayer)
!     6           02/12/2009   All variables stored as layer arrays (Marco Matricardi)
!     7           03/12/2009   Fixed bug (Marco Matricardi)
!     8           07/06/2013   Single stream option to make memory/compute 
!                              requirements more reasonable (A Geer)
!     9           23/06/2014   ircld_k%XSTR initialized to zero to avoid being undefined 
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
  USE rttov_types, ONLY :  &
       & rttov_chanprof, &
       & profile_type,   &
       & rttov_opts_rt_ir, &
       & ircld_type
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
  USE rttov_const, ONLY : cldstr_low_cloud_top
!INTF_ON
  IMPLICIT NONE
!       Arguments with intent in:
  TYPE(rttov_opts_rt_ir), INTENT(IN)    :: opts_rt_ir
  TYPE(rttov_chanprof),   INTENT(IN)    :: chanprof(:)
  TYPE(profile_type),     INTENT(IN)    :: profiles(:)
!       Arguments with intent inout:
  TYPE(profile_type),     INTENT(INOUT) :: profiles_k(SIZE(chanprof))
  TYPE(ircld_type),       INTENT(INOUT) :: ircld
  TYPE(ircld_type),       INTENT(INOUT) :: ircld_k
!INTF_END
!     End of subroutine arguments
!       Local scalars:
  INTEGER(KIND=jpim) :: i, j, istr, ijstr, ilay, ic, jpk
  INTEGER(KIND=jpim) :: nchannels                             ! Number of radiances computed (channels used * profiles)
  INTEGER(KIND=jpim) :: ibdy_layer(1), imax_cfrac(1)
  REAL   (KIND=jprb) :: cfrac_max, cloud_tot(profiles(1)%nlayers) 
!       Local tangent linear arrays:
  REAL   (KIND=jprb) :: cfrac_max_k 
  REAL   (KIND=jprb) :: ntot(profiles(1)%nlevels)
!       Local direct arrays arrays:
  INTEGER(KIND=jpim) :: icount1ref(2*profiles(1)%nlayers)
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
!---------Compute number of streams and cloud distribution in each stream-------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CLDSTR_K', 0_jpim, ZHOOK_HANDLE)
  nchannels = SIZE(chanprof)

  IF (opts_rt_ir%cldstr_simple) THEN

    ! A simpler, faster, but quite approximate "Cmax" single-stream approach. 
    ! Main benefit is that it is much more memory efficient. Intended mainly
    ! for mid- and upper-tropospheric channels.
    DO J = 1, nchannels
      jpk = chanprof(j)%prof

      ! Find maximum cloud fraction and cloud amount above the boundary layer
      ! NB pressure is on levels, cloud in layers, so p(1:nlevels-1) is the top of the layer
      ibdy_layer = MINLOC(ABS(profiles(jpk)%p(1:profiles(jpk)%nlevels-1) - cldstr_low_cloud_top))
      imax_cfrac = MAXLOC(profiles(jpk)%cfrac(1:ibdy_layer(1)))
      cfrac_max  = profiles(jpk)%cfrac(imax_cfrac(1))
      cloud_tot(:) = SUM(profiles(jpk)%cloud(:,:), DIM=1)

      cfrac_max_k = 0.0_jprb

      ! Ignoring trivial amounts of cloud and precip
      IF (ANY(cloud_tot(1:ibdy_layer(1)) > 1E-6_jprb) .AND. cfrac_max > 1E-3) THEN

        ! Cloudy stream required
        cfrac_max_k = cfrac_max_k + ircld_k%XSTR(2, J)
        cfrac_max_k = cfrac_max_k - ircld_k%XSTRCLR(J)

      ENDIF

      profiles_k(j)%cfrac(imax_cfrac(1)) = profiles_k(j)%cfrac(imax_cfrac(1)) + cfrac_max_k

    ENDDO

  ELSE

    DO J = 1, nchannels
      jpk = chanprof(j)%prof
      loop2 : DO istr = 2, ircld%icount(jpk) - 1
        DO i = istr - 1, 1,  - 1
          ircld%xstrref1(istr, I, Jpk) = ircld%xstrref(i, Jpk)
          IF (ircld%xstrref(I, Jpk) <= ircld%a(istr, Jpk)) THEN
            ircld%xstrref(I + 1, Jpk) = ircld%a(istr, Jpk)
            CYCLE loop2
          ELSE
            ircld%xstrref(i + 1, Jpk) = ircld%xstrref(i, Jpk)
          ENDIF
        ENDDO
        i = 0
        ircld%xstrref(I + 1, Jpk) = ircld%a(istr, Jpk)
      ENDDO loop2
      ircld_k%A(:, J)       = 0._jprb
      ircld_k%MAXCOV(:, J)  = 0._jprb
      ircld_k%CLDCFR(:, J)  = 0._jprb
      ircld_k%XSTRMIN(:, J) = 0._jprb
      ircld_k%XSTRMAX(:, J) = 0._jprb
      NTOT(:)               = 0._jprb
  !---------Consider only the streams whose weight is greater than cldstr_threshold-------
      IF (ircld%nstreamref(Jpk) /= 0_jpim) THEN
        IF (ircld%icounstr(Jpk) /= 0_jpim) THEN
          ircld_k%xstr(ircld%icounstr(Jpk) + 1, J) = ircld_k%xstr(ircld%icounstr(Jpk) + 1, J) - ircld_k%xstrclr(j)
          ircld_k%xstr(1, J)                       = ircld_k%xstr(1, J) + ircld_k%xstrclr(j)
          ircld_k%xstrclr(J)                       = 0._jprb
          DO istr = ircld%icounstr(jpk) + 1, 1,  - 1
            IF (istr == 1) THEN
              IF (ircld%indexstr(istr, Jpk) /= istr) THEN
                ircld_k%xstr(ircld%indexstr(istr, Jpk), J) =      &
                  & ircld_k%xstr(ircld%indexstr(istr, Jpk), J) + ircld_k%xstr(istr, J)
                ircld_k%xstr(istr, J)                      = 0._jprb
              ENDIF
            ELSE
              ircld_k%xstr(istr - 1, J) = ircld_k%xstr(istr - 1, j) + ircld_k%xstr(istr, j)
              IF ((ircld%indexstr(istr - 1, Jpk) + 1) /= istr) THEN
                ircld_k%xstr(ircld%indexstr(istr - 1, Jpk) + 1, J) =      &
                  & ircld_k%xstr(ircld%indexstr(istr - 1, Jpk) + 1, J) + ircld_k%xstr(istr, j)
              ELSE
                ircld_k%xstr(ircld%indexstr(istr - 1, Jpk) + 1, J) = ircld_k%xstr(istr, J)
              ENDIF
              IF (ircld%indexstr(istr - 1, Jpk) /= istr) THEN
                ircld_k%xstr(ircld%indexstr(istr - 1, jpk), j) =      &
                  & ircld_k%xstr(ircld%indexstr(istr - 1, jpk), j) - ircld_k%xstr(istr, j)
              ELSE
                ircld_k%xstr(ircld%indexstr(istr - 1, jpk), j) =  - ircld_k%xstr(istr, j)
              ENDIF
              IF (((ircld%indexstr(istr - 1, jpk) + 1) /= istr) .AND. (ircld%indexstr(istr - 1, Jpk) /= istr)) THEN
                ircld_k%xstr(istr, j) = 0._jprb
              ENDIF
            ENDIF
          ENDDO
        ELSE
          ircld_k%XSTRCLR(J) = 0.0_JPRB
        ENDIF
      ENDIF
  !---------Compute the weight of the clear stream--------------------------------
      ircld_k%xstr(ircld%nstreamref(jpk) + 1, j) = ircld_k%xstr(ircld%nstreamref(jpk) + 1, j) - ircld_k%xstrclr(j)
      ircld_k%xstr(1, j)                         = ircld_k%xstr(1, j) + ircld_k%xstrclr(j)
  !---------Re-arrange the limits of each stream in ascending order---------------
      outer : DO i = ircld%iloop(jpk), 1,  - 1
        icount1ref(i) = ircld%icount1ref(i, jpk)
        inner : DO istr = ircld%iloopin(i, jpk), 1,  - 1
          IF (ircld%xstrref2(i, istr, jpk) == ircld%xstrref2(i, istr + 1, jpk)) THEN
            IF (ircld%xstrref2(i, istr, jpk) /= 1.) THEN
              icount1ref(i) = icount1ref(i) - 1
              DO ijstr = icount1ref(i), istr,  - 1
                ircld_k%xstr(ijstr + 1, j) = ircld_k%xstr(ijstr + 1, j) + ircld_k%xstr(ijstr, j)
                ircld_k%xstr(ijstr, J)     = 0._jprb
              ENDDO
            ENDIF
          ENDIF
        ENDDO inner
      ENDDO outer
      loop1 : DO istr = ircld%icount(jpk) - 1, 2,  - 1
        IF (.NOT. ircld%flag(istr, Jpk)) THEN
          i = 0
          ircld_k%a(istr, J)     = ircld_k%a(istr, j) + ircld_k%xstr(i + 1, j)
          ircld_k%xstr(i + 1, j) = 0._jprb                                    !
        ENDIF
        DO i = ircld%iflag(istr, jpk), istr - 1
          IF (ircld%xstrref1(istr, i, jpk) <= ircld%a(istr, jpk)) THEN
            ircld_k%a(istr, j)     = ircld_k%a(istr, j) + ircld_k%xstr(i + 1, j)
            ircld_k%xstr(i + 1, j) = 0._jprb
          ELSE
            ircld_k%xstr(i, j)     = ircld_k%xstr(i, j) + ircld_k%xstr(i + 1, j)
            ircld_k%xstr(i + 1, j) = 0._jprb
          ENDIF
        ENDDO
        ircld_k%xstr(istr, j) = ircld_k%xstr(istr, j) + ircld_k%a(istr, j)
        ircld_k%a(istr, j)    = 0._jprb
      ENDDO loop1
  !---------Determine the limits of each stream----------------------------------
      ic = ircld%icount(jpk)
      DO ilay = profiles(1)%nlayers, 1,  - 1
        IF (ircld%xstrmax(ilay, jpk) /= 0._jprb) THEN
          ic = ic - 1
          ircld_k%xstrmax(ilay, j) = ircld_k%xstrmax(ilay, j) + ircld_k%xstr(ic, j)
          IC = IC - 1
          ircld_k%xstrmin(ilay, j) = ircld_k%xstrmin(ilay, j) + ircld_k%xstr(ic, j)
        ENDIF
        IF (ircld%xstrminref(ilay, jpk) < 0._jprb) THEN
          ircld_k%xstrmin(ilay, j) = 0._jprb
        ENDIF
        ntot(ilay)              = ntot(ilay) + ircld_k%xstrmin(ilay, j)
        ircld_k%cldcfr(ilay, j) = ircld_k%cldcfr(ilay, j) - ircld_k%xstrmin(ilay, j)
        ntot(ilay)              = ntot(ilay) + ircld_k%xstrmax(ilay, j)
      ENDDO
  !---------Compute the cumulative cloud coverage usin the maximum-random---------
  !         overlap assumption
      DO ilay = profiles(1)%nlayers, 1,  - 1
        IF (ircld%cldcfr(ilay, jpk) /= 0._jprb) THEN
          ntot(ilay) =  - ntot(ilay)
        ELSE
          ntot(ilay) = 0._jprb
        ENDIF
      ENDDO
      DO ilay = profiles(1)%nlayers, 2,  - 1
        IF (ircld%ntotref(ilay - 1, jpk) == 0._jprb .AND. (1. - ircld%cldcfr(ilay - 1, Jpk)) == 0._jprb) THEN
          ircld_k%maxcov(ilay, J) = ircld_k%maxcov(ilay, j) + ntot(ilay)
        ELSE
          ntot(ilay - 1)              =      &
            & ntot(ilay - 1) + ntot(ilay) * ircld%maxcov(ilay, jpk) / (1. - ircld%cldcfr(ilay - 1, jpk))
          ircld_k%cldcfr(ilay - 1, J) = ircld_k%cldcfr(ilay - 1, J) +                &
            & ntot(ilay) * ircld%ntotref(ilay - 1, jpk) * ircld%maxcov(ilay, Jpk) /  &
            & (1._jprb - ircld%cldcfr(ilay - 1, Jpk)) ** 2
          ircld_k%maxcov(ilay, J)     =      &
            & ircld_k%maxcov(ilay, j) + ntot(ilay) * ircld%ntotref(ilay - 1, jpk) / (1._jprb - ircld%cldcfr(ilay - 1, Jpk))
        ENDIF
        IF ((ircld%cldcfr(ilay - 1, Jpk) > ircld%cldcfr(ilay, Jpk))) THEN
          ircld_k%cldcfr(ilay - 1, J) = ircld_k%cldcfr(ilay - 1, J) - ircld_k%maxcov(ilay, J)
          ircld_k%MAXCOV(ilay, J)     = 0._jprb
        ELSE IF ((ircld%cldcfr(ilay - 1, Jpk) < ircld%cldcfr(ilay, Jpk))) THEN
          ircld_k%cldcfr(ilay, J) = ircld_k%cldcfr(ILAY, J) - ircld_k%maxcov(ilay, J)
          ircld_k%maxcov(ilay, J) = 0._jprb
        ELSE IF ((ircld%cldcfr(ilay - 1, Jpk) == ircld%cldcfr(ilay, Jpk))) THEN
          ircld_k%cldcfr(ilay, j) = ircld_k%cldcfr(ilay, J) - ircld_k%maxcov(ilay, J)
          ircld_k%maxcov(ilay, J) = 0._jprb
        ENDIF
      ENDDO
      ircld_k%cldcfr(1, J) = ircld_k%cldcfr(1, J) - ntot(1)
  !---------Compute number of streams and cloud distribution in each stream-------
      DO ilay = profiles(1)%nlayers, 1,  - 1
        IF (profiles(jpk)%cfrac(ilay) > 0.) THEN
          profiles_k(j)%cfrac(ilay) = profiles_k(j)%cfrac(ilay) + ircld_k%cldcfr(ilay, j)
        ENDIF
      ENDDO
      ircld_k%CLDCFR(:, j)  = 0._jprb
      ircld_k%XSTRMIN(:, j) = 0._jprb
      ircld_k%XSTRMAX(:, j) = 0._jprb
      ircld_k%XSTR(:, j) = 0._jprb
      ntot = 0._jprb
    ENDDO
  ENDIF
  IF (LHOOK) CALL DR_HOOK('RTTOV_CLDSTR_K', 1_jpim, ZHOOK_HANDLE)
!      ircld_k(:)%xstrclr=0._jprb
END SUBROUTINE rttov_cldstr_k
