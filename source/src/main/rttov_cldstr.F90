!     Compute the number of streams,the area coverage of each stream and
!     the cloud distribution in each stream.
SUBROUTINE rttov_cldstr( &
            & opts_rt_ir,       &
            & profiles,         &
            & ircld,            &
            & nstreams)
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
!     1           30/7/2004    Marco Matricardi. ECMWF.
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
!     3           27/02/2009   Profile levels to include ToA. Cloud inputs on
!                              layers but stored as a level array in profiles.
!                              Distinguish between layer arrays and level
!                              arrays for index labels and looping (P. Rayer)
!     4           02/12/2009   All variables stored as layer arrays (Marco Matricardi)
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!     6           07/06/2013   Single stream option to make memory/compute 
!                              requirements more reasonable (A Geer)
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Module used:
  USE rttov_types, ONLY : profile_type, ircld_type, rttov_opts_rt_ir
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE parkind1, ONLY : jprb
  USE rttov_const, ONLY : ncldtyp, cldstr_low_cloud_top
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE
!       Arguments with intent out:
  INTEGER(KIND=jpim),     INTENT(OUT)   :: nstreams
!       Arguments with intent in:
  TYPE(rttov_opts_rt_ir), INTENT(IN)    :: opts_rt_ir
  TYPE(profile_type),     INTENT(IN)    :: profiles(:)
!       Arguments with intent inout:
  TYPE(ircld_type),       INTENT(INOUT) :: ircld
!INTF_END
!     End of subroutine arguments
!       Local scalars:
  INTEGER(KIND=jpim) :: i, j, istr, ijstr, ilay
  REAL   (KIND=jprb) :: delta_cfrac
  INTEGER(KIND=jpim) :: ibdy_layer(1), imax_cfrac(1)
  REAL   (KIND=jprb) :: cfrac_max, cloud_tot(profiles(1)%nlayers)
!       Local arrays:
  INTEGER(KIND=jpim) :: nprofiles
  REAL   (KIND=jprb) :: ntot(profiles(1)%nlevels)
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
!      ircld(:)%NSTREAM=0_jpim
!      ircld(:)%XSTRCLR=1._jprb
  IF (LHOOK) CALL DR_HOOK('RTTOV_CLDSTR', 0_jpim, ZHOOK_HANDLE)
  nprofiles = SIZE(profiles)

  ! Compute CLDTYP (not used here)
  IF (.NOT. opts_rt_ir%user_cld_opt_param) THEN
    DO J = 1, NPROFILES
      ircld%CLDTYP(:, :, J) = 0._jpim
      DO ILAY = 1, profiles(1)%nlayers
        DO I = 1, ncldtyp
          IF (profiles(j)%cloud(I, Ilay) /= 0.) THEN
            ircld%CLDTYP(I, ILAY, J) = I
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF (opts_rt_ir%cldstr_simple) THEN

    ! A simpler, faster, but quite approximate "Cmax" single-stream approach. 
    ! Main benefit is that it is much more memory efficient. Intended mainly
    ! for mid- and upper-tropospheric channels.
    DO J = 1, NPROFILES

      ircld%ICLDARR(:, :, j) = 0_jpim

      ! Find maximum cloud fraction and cloud amount above the boundary layer
      ! NB pressure is on levels, cloud in layers, so p(1:nlevels-1) is the top of the layer
      ibdy_layer = MINLOC(ABS(profiles(j)%p(1:profiles(j)%nlevels-1) - cldstr_low_cloud_top))
      imax_cfrac = MAXLOC(profiles(j)%cfrac(1:ibdy_layer(1)))
      cfrac_max  = profiles(j)%cfrac(imax_cfrac(1))
      cloud_tot(:) = SUM(profiles(j)%cloud(:,:), DIM=1)

      ! Ignoring trivial amounts of cloud and precip
      IF (ANY(cloud_tot(1:ibdy_layer(1)) > 1E-6_jprb) .AND. cfrac_max > 1E-3) THEN

        ! Cloudy stream required
        ircld%NSTREAM(j)     = 1
        ircld%XSTR(1, J)     = 0._jprb
        ircld%XSTR(2, J)     = cfrac_max
        ircld%XSTRCLR(J)     = 1._jprb - cfrac_max
        WHERE (cloud_tot > 1E-6_jprb) ircld%ICLDARR(1,:,J) = 1

      ELSE

        ! Pure clear-sky
        ircld%NSTREAM(j)     = 0
        ircld%XSTR(1:2, J)   = 0._jprb
        ircld%XSTRCLR(J)     = 1._jprb

      ENDIF

    ENDDO
    nstreams = MAXVAL(ircld%NSTREAM(1:nprofiles))

  ELSE

    nstreams  = 0_jpim
    delta_cfrac = 10.0_jprb * EPSILON(1.0_jprb)
    DO J = 1, NPROFILES
  !---------Compute number of streams and cloud distribution in each stream-------
      ircld%iflag(:, j)      = 1_jpim
      ircld%ICLDARR(:, :, j) = 0_jpim
      ircld%CLDCFR(:, j)     = 0._jprb
      ircld%XSTR(:, J)       = 0._jprb
      ircld%XSTRMIN(:, j)    = 0._jprb
      ircld%XSTRMAX(:, j)    = 0._jprb
      ircld%FLAG(:, j)       = .FALSE.
      NTOT = 0._jprb
      DO ILAY = 1, profiles(1)%nlayers
        IF (profiles(j)%cfrac(ILay) > 0.) THEN
          ircld%CLDCFR(ILAY, J) = profiles(j)%cfrac(ILay)
        ENDIF
          ENDDO

      ! Check for overcast layers and identical cfrac values on consecutive layers
      ! These need adjusting to ensure TL/AD/K models are correct
      DO ILAY = 2, profiles(1)%nlayers
        IF (ircld%CLDCFR(ILAY, J) > 0.) THEN

          ! Check for overcast layers
          IF (ircld%CLDCFR(ILAY, J) >= 1.0_jprb) THEN
            ircld%CLDCFR(ILAY, J) = 1.0_jprb - ilay*delta_cfrac
          ENDIF

          ! Check for identical adjacent cfrac (note that this won't always work if cldstr_threshold is +ve
          IF (ircld%CLDCFR(ILAY, J) == ircld%CLDCFR(ILAY-1, J)) THEN
            ircld%CLDCFR(ILAY, J) = ircld%CLDCFR(ILAY, J) - SIGN(delta_cfrac, ircld%CLDCFR(ILAY, J)-0.5_jprb)
          ENDIF

        ENDIF
      ENDDO

      ircld%ICOUNT(J)     = 0_jpim
      ircld%ICOUNSTR(J)   = 0_jpim
  !---------Compute the cumulative cloud coverage usin the maximum-random---------
  !         overlap assumption
      NTOT(1)             = 1. - ircld%CLDCFR(1, J)
      ircld%NTOTREF(1, J) = NTOT(1)
      DO ILAY = 2, profiles(1)%nlayers
        ircld%MAXCOV(ILAY, J) = (1. - MAX(ircld%CLDCFR(ILAY - 1, J), ircld%CLDCFR(ILAY, J)))
        IF (NTOT(ILAY - 1) == 0 .AND. (1. - ircld%CLDCFR(ILAY - 1, J)) == 0._JPRB) THEN
          NTOT(ILAY)             = ircld%MAXCOV(ILAY, J)
        ELSE
          NTOT(ILAY)             = NTOT(ILAY - 1) * ircld%MAXCOV(ILAY, J) / (1. - ircld%CLDCFR(ILAY - 1, J))
        ENDIF
          ircld%NTOTREF(ILAY, J) = NTOT(ILAY)
      ENDDO
      DO ILAY = 1, profiles(1)%nlayers
        IF (ircld%CLDCFR(ILAY, J) /= 0._JPRB) THEN
          NTOT(ILAY) = 1._JPRB - NTOT(ILAY)
        ELSE
          NTOT(ILAY) = 0._JPRB
        ENDIF
      ENDDO
  !---------Determine the limits of each stream----------------------------------
      ircld%ICOUNT(J) = 1_jpim
      DO ILAY = 1, profiles(1)%nlayers
        ircld%XSTRMAX(ILAY, J)    = NTOT(ILAY)
        ircld%XSTRMIN(ILAY, J)    = NTOT(ILAY) - ircld%CLDCFR(ILAY, J)
        ircld%XSTRMINREF(ILAY, J) = ircld%XSTRMIN(ILAY, J)
        IF (ircld%XSTRMIN(ILAY, J) < 0._JPRB) ircld%XSTRMIN(ILAY, J) = 0._JPRB
        IF (ircld%XSTRMAX(ILAY, J) /= 0._JPRB) THEN
          ircld%XSTR(ircld%ICOUNT(J), J) = ircld%XSTRMIN(ILAY, J)
          ircld%ICOUNT(J)                = ircld%ICOUNT(J) + 1
          ircld%XSTR(ircld%ICOUNT(J), J) = ircld%XSTRMAX(ILAY, J)
          ircld%ICOUNT(J)                = ircld%ICOUNT(J) + 1
        ENDIF
      ENDDO
      ircld%XSTRREF(:, J) = ircld%XSTR(:, J)
  !---------Re-arrange the limits of each stream in ascending order---------------
      LOOP1 : DO ISTR = 2, ircld%ICOUNT(J) - 1
        ircld%A(ISTR, J) = ircld%XSTR(ISTR, J)
        DO I = ISTR - 1, 1,  - 1
          IF (ircld%XSTR(I, J) <= ircld%A(ISTR, J)) THEN
            ircld%XSTR(I + 1, J) = ircld%A(ISTR, J)
            ircld%IFLAG(ISTR, J) = I
            ircld%FLAG(ISTR, J)  = .TRUE.
            CYCLE LOOP1
          ELSE
            ircld%XSTR(I + 1, J) = ircld%XSTR(I, J)
          ENDIF
        ENDDO
        I = 0
        ircld%XSTR(I + 1, J) = ircld%A(ISTR, J)
      ENDDO LOOP1
      ircld%ICOUNT1(J)    = ircld%ICOUNT(J) - 1
      ircld%ILOOP(J)      = 0_jpim
      ircld%ILOOPIN(:, J) = 0_jpim
      OUTER : DO
        ircld%ILOOP(J)                       = ircld%ILOOP(J) + 1
        ircld%ICOUNT1REF(ircld%ILOOP(J), J)  = ircld%ICOUNT1(J)
        ircld%XSTRREF2(ircld%ILOOP(J), :, J) = ircld%XSTR(:, J)
        INNER : DO ISTR = 1, ircld%ICOUNT1(J) - 1
          ircld%ILOOPIN(ircld%ILOOP(J), J) = ircld%ILOOPIN(ircld%ILOOP(J), J) + 1
          IF (ircld%XSTR(ISTR, J) == ircld%XSTR(ISTR + 1, J)) THEN
            IF (ircld%XSTR(ISTR, J) /= 1._jprb) THEN
              ircld%ICOUNT1(J) = ircld%ICOUNT1(J) - 1
              DO IJSTR = ISTR, ircld%ICOUNT1(J)
                ircld%XSTR(IJSTR, J) = ircld%XSTR(IJSTR + 1, J)
              ENDDO
              CYCLE OUTER
            ELSE
              EXIT OUTER
            ENDIF
          ENDIF
        ENDDO INNER
        EXIT OUTER
      ENDDO OUTER
      ircld%NSTREAM(J) = ircld%ICOUNT1(J) - 1
  !---------Compute the weight of the clear stream------------------------------
      DO ISTR = 1, ircld%NSTREAM(J)
        DO ILAY = 1, profiles(1)%nlayers
          IF (ircld%XSTRMIN(ILAY, J) <= ircld%XSTR(ISTR, J) .AND. ircld%XSTRMAX(ILAY, J) >= ircld%XSTR(ISTR + 1, J)) THEN
            ircld%ICLDARR(ISTR, ILAY, J) = 1_jpim
          ENDIF
        ENDDO
      ENDDO
      IF (ircld%NSTREAM(J) ==  - 1_jpim) THEN
        ircld%NSTREAM(J) = 0_jpim
      ENDIF
      ircld%XSTRCLR(J)    = 1._jprb - (ircld%XSTR(ircld%NSTREAM(J) + 1, J) - ircld%XSTR(1, J))
      ircld%NSTREAMREF(J) = ircld%NSTREAM(J)
  !---------Consider only the streams whose weight is greater than cldstr_threshold-------
      IF (ircld%NSTREAM(J) /= 0_jpim) THEN
        DO ISTR = 1, ircld%NSTREAM(J)
          IF ((ircld%XSTR(ISTR + 1, J) - ircld%XSTR(ISTR, J)) >= opts_rt_ir%cldstr_threshold) THEN
            ircld%ICOUNSTR(J)                    = ircld%ICOUNSTR(J) + 1
            ircld%INDEXSTR(ircld%ICOUNSTR(J), J) = ISTR
          ENDIF
        ENDDO
        IF (ircld%ICOUNSTR(J) /= 0_jpim) THEN
          DO ISTR = 1, ircld%ICOUNSTR(J)
            DO ILAY = 1, profiles(1)%nlayers
              ircld%ICLDARR(ISTR, ILAY, J) = ircld%ICLDARR(ircld%INDEXSTR(ISTR, J), ILAY, J)
            ENDDO
          ENDDO
          DO ISTR = 1, ircld%ICOUNSTR(J) + 1
            IF (ISTR == 1_jpim) THEN
              ircld%XSTR(ISTR, J) = ircld%XSTR(ircld%INDEXSTR(ISTR, J), J)
            ELSE
              ircld%XSTR(ISTR, J) = ircld%XSTR(ISTR - 1, J) +      &
                & (ircld%XSTR(ircld%INDEXSTR(ISTR - 1, J) + 1, J) - ircld%XSTR(ircld%INDEXSTR(ISTR - 1, J), J))
            ENDIF
          ENDDO
          ircld%XSTRCLR(J) = 1._jprb - (ircld%XSTR(ircld%ICOUNSTR(J) + 1, J) - ircld%XSTR(1, J))
          ircld%NSTREAM(J) = ircld%ICOUNSTR(J)
        ELSE IF (ircld%ICOUNSTR(J) == 0_jpim) THEN
          ircld%XSTRCLR(J) = 1._jprb
          ircld%NSTREAM(J) = 0_jpim
        ENDIF
      ELSE
        ircld%NSTREAM(J) = 0_jpim
      ENDIF
      nstreams = max(nstreams, ircld%NSTREAM(J))
    ENDDO

  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_CLDSTR', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_cldstr
