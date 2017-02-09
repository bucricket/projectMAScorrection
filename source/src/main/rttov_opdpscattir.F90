!!    Compute Optical parameters for aerosols and clouds
SUBROUTINE rttov_opdpscattir( &
            & nlayers,                      &
            & chanprof,                     &
            & opts,                         &
            & aux,                          &
            & profiles,                     &
            & profiles_dry,                 &
            & aer_opt_param,                &
            & cld_opt_param,                &
            & dosolar,                      &
            & solar,                        &
            & coef,                         &
            & coef_scatt_ir,                &
            & raytracing,                   &
            & transmission_scatt_ir,        &
            & transmission_scatt_ir_stream, &
            & optp,                         &
            & ircld)
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
! Current Code Owner: SAF NWP
!
!     Method:
!     See: User's manual and scientific report for RTIASI-5
!          (Available from EUMETSAT)
!
!     Owner:
!     EUMETSAT
!
!     History:
!     Version      Date        Comment
!     1           30/7/2004    RTIASI-5. Marco Matricardi. ECMWF.
!     2           27/02/2009   Profile levels to include ToA. Distinguish
!                              arrays in raytracing and profiles (on levels)
!                              from all others (on layers) - size, index
!                              labels, looping (P. Rayer)
!     3           09/06/2009   Corrected bug - to now assign ircld%xpresave from
!                              user pressure array, not refprf (R.Saunders)
!     4           03/11/2009   Transmittances / optical depths on levels (A Geer)
!     5           02/12/2009   Fixed a number of bugs due to the wrong assumption that aerosol/cloud
!                              related quantities are defined on levels (thay are layer
!                              average quantities). Marco Matricardi
!     6           02/12/2009   Introduced multiple cloud types in a single layer. pathsun, pathsat and
!                              related quantities are now layer arrays (Marco Matricardi)
!     7           05/07/2010   Remove addsolar flag from profiles structure (J Hocking)
!     8           14/12/2010   Use traj0_sta%solar array to flag channels for which solar calculations
!                              should be performed (J Hocking)
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
       & rttov_chanprof,             &
       & rttov_options,              &
       & rttov_coef,                 &
       & profile_Type,               &
       & rttov_opt_param,            &
       & raytracing_type,            &
       & transmission_scatt_ir_type, &
       & rttov_coef_scatt_ir,        &
       & rttov_optpar_ir,            &
       & profile_aux,                &
       & ircld_type
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jprb
  USE rttov_const, ONLY :  &
       & deg2rad, &
       & e00,     &
       & t00,     &
       & ti,      &
       & ncldtyp, &
       & phangle
!INTF_ON
  IMPLICIT NONE
  INTEGER(KIND=jpim)              , INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof            ), INTENT(IN)    :: chanprof(:)
  TYPE(rttov_options )            , INTENT(IN)    :: opts
  TYPE(profile_type              ), INTENT(IN)    :: profiles(:)
  TYPE(profile_type              ), INTENT(IN)    :: profiles_dry(SIZE(profiles))
  TYPE(rttov_opt_param), OPTIONAL , INTENT(IN)    :: aer_opt_param
  TYPE(rttov_opt_param), OPTIONAL , INTENT(IN)    :: cld_opt_param
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: dosolar
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: solar(SIZE(chanprof))
  TYPE(profile_aux               ), INTENT(INOUT) :: aux
  TYPE(rttov_coef                ), INTENT(IN)    :: coef
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream
  TYPE(rttov_coef_scatt_ir       ), INTENT(IN)    :: coef_scatt_ir
  TYPE(rttov_optpar_ir           ), INTENT(IN)    :: optp
  TYPE(ircld_type                ), INTENT(INOUT) :: ircld
  TYPE(raytracing_type           ), INTENT(INOUT) :: raytracing
!INTF_END
!     End of subroutine arguments
#include "rttov_baran2013_calc_optpar.interface"
#include "rttov_baran2014_calc_optpar.interface"
#include "rttov_baran_calc_phase.interface"
!       Local scalars:
  INTEGER(KIND=jpim) :: j, chan, chan1, i, prof, ish, ist, isti, k, ityp, iae
  INTEGER(KIND=jpim) :: lev, lay, lctyp
  INTEGER(KIND=jpim) :: icount, icl
  REAL   (KIND=jprb) :: opd, opdsun
  REAL   (KIND=jprb) :: absch
  REAL   (KIND=jprb) :: scach
  REAL   (KIND=jprb) :: bparh
  REAL   (KIND=jprb) :: afac
  REAL   (KIND=jprb) :: sfac
  REAL   (KIND=jprb) :: gfac
  REAL   (KIND=jprb) :: deltadg
  REAL   (KIND=jprb) :: delth
  REAL   (KIND=jprb) :: zdelth
  REAL   (KIND=jprb) :: frach
  REAL   (KIND=jprb) :: phasint
  REAL   (KIND=jprb) :: musat
  REAL   (KIND=jprb) :: musun
!       Local arrays:
  REAL   (KIND=jprb) :: OPDPAERL   (nlayers,SIZE(chanprof))
  REAL   (KIND=jprb) :: OPDPCLDL   (nlayers,SIZE(chanprof))
  REAL   (KIND=jprb) :: OPDPCLDLSUN(nlayers,SIZE(chanprof))
  REAL   (KIND=jprb) :: OPDPAERLSUN(nlayers,SIZE(chanprof))
  REAL   (KIND=jprb) :: pfac       (coef_scatt_ir%fmv_aer_ph)
  REAL   (KIND=jprb) :: phash      (coef_scatt_ir%fmv_aer_ph)
  REAL   (KIND=jprb) :: phasice    (coef_scatt_ir%fmv_icl_ph)
  REAL   (KIND=jprb) :: deltap     (coef_scatt_ir%fmv_icl_ph)
  REAL   (KIND=jprb) :: zminphadiff, relazi
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: nchannels
  
  REAL(KIND=jprb)  :: abso
  REAL(KIND=jprb)  :: sca
  REAL(KIND=jprb)  :: bpr
  REAL(KIND=jprb)  :: asym

  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDPSCATTIR', 0_jpim, ZHOOK_HANDLE)
  nprofiles   = SIZE(profiles)
  nchannels   = SIZE(chanprof)
  bparh       = 0._jprb
  OPDPAERL    = 0._jprb
  OPDPCLDL    = 0._jprb
  IF (dosolar) THEN
    OPDPAERLSUN = 0._jprb
    OPDPCLDLSUN = 0._jprb
  ENDIF
!-----------------------------------------------------------------------------------------
!         1.   CALCULATE OPTICAL DEPTHS OF AEROSOLS
!-----------------------------------------------------------------------------------------
  IF (opts%rt_ir%addaerosl .AND. .NOT. opts%rt_ir%user_aer_opt_param) THEN
    DO j = 1, nprofiles
      DO lay = 1, nlayers
        lev    = lay + 1
        icount = 0
        DO icl = 1, coef_scatt_ir%fmv_aer_comp
          IF (profiles(j)%aerosols(icl, lay) /= 0._jprb) THEN
            icount = icount + 1
            aux%iaertyp(icount,lay,j) = icl
          ENDIF
        ENDDO
        aux%iaernum(lay,j)    = icount
!-----Compute relative humidity-----------------------------------------------------------
        ircld%tave(lay,j)     = (profiles(j)%t(lev - 1) + profiles(j)%t(lev)) / 2._jprb
        ircld%wmixave(lay,j)  = (profiles_dry(j)%q(lev - 1) + profiles_dry(j)%q(lev)) / 2._jprb
        ircld%xpresave(lay,j) = (profiles(j)%p(lev - 1) + profiles(j)%p(lev)) / 2._jprb
! saturated vapour pressure
        ircld%esw(lay,j) = e00 * exp(17.502_jprb * (ircld%tave(lay,j) - T00) / (ircld%tave(lay,j) - 32.19_jprb))
        ircld%esi(lay,j) = E00 * exp(22.587_jprb * (ircld%tave(lay,j) - T00) / (ircld%tave(lay,j) + 0.7_jprb))
        IF (ircld%tave(lay,j) > t00) THEN
          ircld%ppv(lay,j) = ircld%esw(lay,j)! Water phase
        ELSE IF (ircld%tave(lay,j) > ti .AND. ircld%tave(lay,j) <= t00) THEN
          ircld%ppv(lay,j) = ircld%esi(lay,j) + &
            & (ircld%esw(lay,j) - ircld%esi(lay,j)) * ((ircld%tave(lay,j) - ti) / (t00 - ti)) ** 2! Mixed phase
        ELSE IF (ircld%tave(lay,j) <= ti) THEN
          ircld%ppv(lay,j) = ircld%esi(lay,j)! Ice phase
        ENDIF
        ircld%ppv(lay,j)     = ircld%ppv(lay,j) / 100._jprb
! layer average relative humidity
        aux%relhum(lay,j)    = 100._jprb * ircld%wmixave(lay,j) * &
          & 1.e-6_jprb * ircld%xpresave(lay,j) /      &
          & (ircld%ppv(lay,j) * (1._jprb + ircld%wmixave(lay,j) * 1.e-6_jprb))
        aux%relhumref(lay,j) = aux%relhum(lay,j)
        IF (aux%relhum(lay,j) > 99._jprb) THEN
          aux%relhum(lay,j) = 99._jprb
        ENDIF
!----------------------------------------------------------------------------------------
      ENDDO ! layers
    ENDDO ! profiles
  ENDIF ! opts%rt_ir%addaerosl

  IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
    transmission_scatt_ir%OPDPAAER = 0._jprb
    transmission_scatt_ir%OPDPSAER = 0._jprb
    IF (dosolar) THEN
      transmission_scatt_ir%AZPHAERUPA = 0._jprb
      transmission_scatt_ir%AZPHAERDOA = 0._jprb
    ENDIF
  ENDIF

  IF (opts%rt_ir%addaerosl) THEN
    transmission_scatt_ir%GPARAERA = 0._jprb
    transmission_scatt_ir%GPARAER  = 0._jprb
    IF (dosolar) THEN
      transmission_scatt_ir%AZPHAERUP = 0._jprb
      transmission_scatt_ir%AZPHAERDO = 0._jprb
    ENDIF

    DO j = 1, nchannels
      chan  = chanprof(j)%chan
      prof  = chanprof(j)%prof
      relazi  = profiles(prof)%azangle - profiles(prof)%sunazangle

      DO lay = 1, nlayers

        IF (opts%rt_ir%user_aer_opt_param) THEN
          transmission_scatt_ir%OPDPAAER(lay,j) = aer_opt_param%abs(j,lay) * &
                                                   raytracing%ltick(lay,prof)
          transmission_scatt_ir%OPDPSAER(lay,j) = aer_opt_param%sca(j,lay) * &
                                                   raytracing%ltick(lay,prof)
          transmission_scatt_ir%GPARAERA(lay,j) = aer_opt_param%sca(j,lay) * &
                                                   aer_opt_param%bpr(j,lay) * &
                                                   raytracing%ltick(lay,prof)

          IF (solar(j)) THEN
!-------------Average phase function for the upward scattered solar beam------------------
            musat   = 1._jprb / raytracing%pathsat(lay,prof)
            musun   =  - 1._jprb / raytracing%pathsun(lay,prof)
            zminphadiff = 1._jprb / (aer_opt_param%minphadiff * deg2rad)

            CALL int_phase_fn(musat, musun, 180.0_jprb - relazi, zminphadiff, &
                              aer_opt_param%pha(j,lay,:), aer_opt_param%cosphangle, &
                              aer_opt_param%iphangle, phasint)

            transmission_scatt_ir%PHASINTUPREF(1,lay,j) = phasint
            transmission_scatt_ir%AZPHAERUP(lay,j) = transmission_scatt_ir%AZPHAERUP(lay,j) + &
                                                      phasint * aer_opt_param%sca(j,lay) * &
                                                      raytracing%ltick(lay,prof)
!-------------Average phase function for the downward scattered solar beam----------------
            ! musun and zminphadiff are same as above; musat is -1 times the value above
            musat =  - 1._jprb * musat

            CALL int_phase_fn(musat, musun, relazi, zminphadiff,                    &
                              aer_opt_param%pha(j,lay,:), aer_opt_param%cosphangle, &
                              aer_opt_param%iphangle, phasint)

            transmission_scatt_ir%PHASINTDOREF(1,lay,j) = phasint
            transmission_scatt_ir%AZPHAERDO(lay,j) = transmission_scatt_ir%AZPHAERDO(lay,j) + &
                                                      phasint * aer_opt_param%sca(j,lay) * &
                                                      raytracing%ltick(lay,prof)
          ENDIF
        ELSE
          DO I = 1, aux%iaernum(lay,prof)
            iae = aux%iaertyp(i,lay,prof)
            IF (coef_scatt_ir%fmv_aer_rh(iae) /= 1) THEN
!---------------Interpolate scattering parameters to actual value of relative humidity-----
              DO k = 1, coef_scatt_ir%fmv_aer_rh(iae) - 1
                IF (aux%relhum(lay,prof) >= optp%optpaer(iae)%fmv_aer_rh_val(k) .AND.      &
                  & aux%relhum(lay,prof) <= optp%optpaer(iae)%fmv_aer_rh_val(k + 1)) THEN
                  delth  = (optp%optpaer(iae)%fmv_aer_rh_val(K + 1) - optp%optpaer(iae)%fmv_aer_rh_val(K))
                  zdelth = 1.0_jprb / delth
                  frach  = (aux%relhum(lay,prof) - optp%optpaer(iae)%fmv_aer_rh_val(K))
                  afac   = (optp%optpaer(iae)%abs(chan,k + 1) - optp%optpaer(iae)%abs(chan,k)) * zdelth
                  sfac   = (optp%optpaer(iae)%sca(chan,k + 1) - optp%optpaer(iae)%sca(chan,k)) * zdelth
                  gfac   = (optp%optpaer(iae)%bpr(chan,k + 1) - optp%optpaer(iae)%bpr(chan,k)) * zdelth
                  absch  = optp%optpaer(iae)%abs(chan,k) + afac * frach
                  scach  = optp%optpaer(iae)%sca(chan,k) + sfac * frach
                  bparh  = optp%optpaer(iae)%bpr(chan,k) + gfac * frach
                  IF (solar(j)) THEN
                    chan1 = coef_scatt_ir%aer_pha_index(chan)
                    pfac(1:coef_scatt_ir%fmv_aer_ph)  = &
                        & (optp%optpaer(iae)%pha(chan1,k + 1,1:coef_scatt_ir%fmv_aer_ph) - &
                        & optp%optpaer(iae)%pha(chan1,k,1:coef_scatt_ir%fmv_aer_ph)) * zdelth
                    phash(1:coef_scatt_ir%fmv_aer_ph) = &
                        & optp%optpaer(iae)%pha(chan1,k,1:coef_scatt_ir%fmv_aer_ph) + &
                        & pfac(1:coef_scatt_ir%fmv_aer_ph) * frach
                  ENDIF
                  EXIT
                ENDIF
              ENDDO
            ELSE
              absch = optp%optpaer(iae)%abs(chan,1)
              scach = optp%optpaer(iae)%sca(chan,1)
              bparh = optp%optpaer(iae)%bpr(chan,1)
              IF (solar(j)) THEN
                chan1 = coef_scatt_ir%aer_pha_index(chan)
                phash(1:coef_scatt_ir%fmv_aer_ph) = optp%optpaer(iae)%pha(chan1, 1, 1:coef_scatt_ir%fmv_aer_ph)
              ENDIF
            ENDIF
!-------------Compute optical parameters considering the contribution of------------------
!             all the aerosol components present in the layer
! For pre-defined particle types OD calculated as (aerosol amount in layer) * 
!   (aerosol ext coeff for given layer RH) * (layer depth as height diff of two upwardly successive levels)
! OD is accumulated for each aerosol type - arrays are zeroed above.
            transmission_scatt_ir%OPDPAAER(lay,j) = transmission_scatt_ir%OPDPAAER(lay,j) +      &
              & profiles(prof)%aerosols(iae,lay) * ABSCH * raytracing%ltick(lay,prof)
            transmission_scatt_ir%OPDPSAER(lay,j) = transmission_scatt_ir%OPDPSAER(lay,j) +      &
              & profiles(prof)%aerosols(iae,lay) * SCACH * raytracing%ltick(lay,prof)
            transmission_scatt_ir%GPARAERA(lay,j) = transmission_scatt_ir%GPARAERA(lay,j) +      &
              & profiles(prof)%aerosols(iae,lay) * SCACH * BPARH * raytracing%ltick(lay,prof)
!-------------If solar radiation is present,compute the azimuthally averaged--------------
!             value of the phase function for the given value of the viewing--------------
!             angle and solar zenith angle.
            IF (solar(j)) THEN
!-------------Average phase function for the upward scattered solar beam------------------
              musat   = 1._jprb / raytracing%pathsat(lay,prof)
              musun   =  - 1._jprb / raytracing%pathsun(lay,prof)
              zminphadiff = 1._jprb / (coef_scatt_ir%fmv_aer_ph_val_min * deg2rad)

              CALL int_phase_fn(musat, musun, 180.0_jprb - relazi, zminphadiff, &
                                phash, coef_scatt_ir%fmv_aer_ph_val_cos, &
                                coef_scatt_ir%ifmv_aer_ph_val, phasint)

              transmission_scatt_ir%PHASINTUPREF(i,lay,j) = phasint
              transmission_scatt_ir%AZPHAERUP(lay,j) = transmission_scatt_ir%AZPHAERUP(lay,j) + &
                & profiles(prof)%aerosols(iae,lay) * phasint * SCACH * raytracing%ltick(lay,prof)
!-------------Average phase function for the downward scattered solar beam----------------
              ! musun and zminphadiff are same as above; musat is -1 times the value above
              musat =  - 1._jprb * musat

              CALL int_phase_fn(musat, musun, relazi, zminphadiff,       &
                                phash, coef_scatt_ir%fmv_aer_ph_val_cos, &
                                coef_scatt_ir%ifmv_aer_ph_val, phasint)

              transmission_scatt_ir%PHASINTDOREF(i,lay,j) = phasint
              transmission_scatt_ir%AZPHAERDO(lay,j) = transmission_scatt_ir%AZPHAERDO(lay,j) +      &
                & profiles(prof)%aerosols(iae,lay) * phasint * SCACH * raytracing%ltick(lay,prof)
            ENDIF ! solar(j)
          ENDDO ! aer types
        ENDIF ! user_aer_opt_param
      ENDDO ! layers
!---------Compute final values for optical parameters-------------------------------------
      DO lay = 1, nlayers
        IF (transmission_scatt_ir%OPDPSAER(lay,j) /= 0._jprb) THEN
          transmission_scatt_ir%GPARAER(lay,j) =      &
            & transmission_scatt_ir%GPARAERA(lay,j) / transmission_scatt_ir%OPDPSAER(lay,j)
          IF (solar(j)) THEN
            transmission_scatt_ir%AZPHAERUPA(lay,j) =      &
              & transmission_scatt_ir%AZPHAERUP(lay,j) / transmission_scatt_ir%OPDPSAER(lay,j)
            transmission_scatt_ir%AZPHAERDOA(lay,j) =      &
              & transmission_scatt_ir%AZPHAERDO(lay,j) / transmission_scatt_ir%OPDPSAER(lay,j)
          ENDIF
        ENDIF
        transmission_scatt_ir%OPDPAERLA(lay,j) = transmission_scatt_ir%OPDPAAER(lay,j) +      &
          & transmission_scatt_ir%OPDPSAER(lay,j) * transmission_scatt_ir%GPARAER(lay,j)
        OPDPAERL(lay,j)                        =      &
          & transmission_scatt_ir%OPDPAERLA(lay,j) * raytracing%pathsat(lay,prof) * coef%ff_gam(chan)
        IF (solar(j)) THEN
          OPDPAERLSUN(lay,j) = &
            & transmission_scatt_ir%OPDPAERLA(lay,j) * raytracing%patheff(lay,prof) * coef%ff_gam(chan)
        ENDIF
      ENDDO ! layers
    ENDDO ! chanprof
  ENDIF ! opts%rt_ir%addaerosl
!-----------------------------------------------------------------------------------------
!         2.   CALCULATE OPTICAL DEPTHS OF CLOUDS
!-----------------------------------------------------------------------------------------
  IF (opts%rt_ir%addclouds) THEN
    transmission_scatt_ir%OPDPA   = 0._jprb
    transmission_scatt_ir%OPDPS   = 0._jprb
    transmission_scatt_ir%GPAR    = 0._jprb
    transmission_scatt_ir%GPARTOT = 0._jprb
    IF (dosolar) THEN
      transmission_scatt_ir%AZPHUPTOT = 0._jprb
      transmission_scatt_ir%AZPHDOTOT = 0._jprb
      transmission_scatt_ir%AZPHUP = 0._jprb
      transmission_scatt_ir%AZPHDO = 0._jprb
    ENDIF

    DO j = 1, nchannels
      chan  = chanprof(j)%chan
      prof  = chanprof(j)%prof
      ish  = profiles(prof)%ish
      relazi  = profiles(prof)%azangle - profiles(prof)%sunazangle
      DO lay = 1, nlayers

        IF (opts%rt_ir%user_cld_opt_param) THEN
          transmission_scatt_ir%OPDPA(lay,j) = cld_opt_param%abs(j,lay) * raytracing%ltick(lay,prof)
          transmission_scatt_ir%OPDPS(lay,j) = cld_opt_param%sca(j,lay) * raytracing%ltick(lay,prof)
          transmission_scatt_ir%GPARTOT(lay,j) = cld_opt_param%bpr(j,lay) * cld_opt_param%sca(j,lay) * &
                                               & raytracing%ltick(lay,prof)

          IF (solar(j)) THEN
!-------------Average phase function for the upward scattered solar beam------------------
            musat       = 1._jprb / raytracing%pathsat(lay,prof)
            musun       =  - 1._jprb / raytracing%pathsun(lay,prof)
            zminphadiff = 1._jprb / (cld_opt_param%minphadiff * deg2rad)

            CALL int_phase_fn(musat, musun, 180.0_jprb - relazi, zminphadiff, &
                              cld_opt_param%pha(j,lay,:), cld_opt_param%cosphangle, &
                              cld_opt_param%iphangle, phasint)

            transmission_scatt_ir%AZPHUPCLS(1,lay,j) = phasint
            transmission_scatt_ir%AZPHUPTOT(lay,j) = transmission_scatt_ir%AZPHUPTOT(lay,j) + &
                                                    & phasint * transmission_scatt_ir%OPDPS(lay,j)
!-------------Average phase function for the downward scattered solar beam----------------
            ! musun and zminphadiff are same as above; musat is -1 times the value above
            musat =  - 1._jprb * musat

            CALL int_phase_fn(musat, musun, relazi, zminphadiff,                    &
                              cld_opt_param%pha(j,lay,:), cld_opt_param%cosphangle, &
                              cld_opt_param%iphangle, phasint)

            transmission_scatt_ir%AZPHDOCLS(1,lay,j) = phasint
            transmission_scatt_ir%AZPHDOTOT(lay,j) = transmission_scatt_ir%AZPHDOTOT(lay,j) + &
                                                    & phasint * transmission_scatt_ir%OPDPS(lay,j)
          ENDIF
        ELSE
          DO lctyp = 1, ncldtyp
            IF (ircld%cldtyp(lctyp, lay,prof) /= 0) THEN
              ityp  = ircld%cldtyp(lctyp, lay,prof)
!---------------Compute cloud  optical parameters ----------------------------------------
              IF (ityp <= 5_jpim) THEN
!-----------------------------------------------------------------------------------------
!                 For water clouds use stored optical parameters
!-----------------------------------------------------------------------------------------
                transmission_scatt_ir%OPDPACLS(ityp,lay,j) =                                                      &
                  & profiles(prof)%cloud(ityp,lay) * coef_scatt_ir%confac(ityp) * optp%optpwcl(ityp)%abs(chan,1) * &
                  & raytracing%ltick(lay,prof)
                transmission_scatt_ir%OPDPA(lay,j)          =      &
                  & transmission_scatt_ir%OPDPA(lay,j) + transmission_scatt_ir%OPDPACLS(ityp,lay,j)
                transmission_scatt_ir%OPDPSCLS(ityp,lay,j) =  +                                                   &
                  & profiles(prof)%cloud(ityp,lay) * coef_scatt_ir%confac(ityp) * optp%optpwcl(ityp)%sca(chan,1) * &
                  & raytracing%ltick(lay,prof)
                transmission_scatt_ir%OPDPS(lay,j)          =      &
                  & transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSCLS(ityp,lay,j)
                transmission_scatt_ir%GPARCLS(ityp,lay,j)  = optp%optpwcl(ityp)%bpr(chan,1)
                transmission_scatt_ir%GPARTOT(lay,j)        = transmission_scatt_ir%GPARTOT(lay,j) +      &
                  & transmission_scatt_ir%GPARCLS(ityp,lay,j) * transmission_scatt_ir%OPDPSCLS(ityp,lay,j)
              ELSE
!-----------------------------------------------------------------------------------------
!                 For ice clouds optical parameters are computed using regression
!                 coefficients
!-----------------------------------------------------------------------------------------
                IF( ish .LE. 2_jpim ) THEN
                  transmission_scatt_ir%OPDPACLS(ityp,lay,j) = profiles(prof)%cloud(ityp,lay) * &
                    & (optp%optpicl(ish)%abs(chan,1) + optp%optpicl(ish)%abs(chan,2) * aux%dg(lay,prof) + &
                    & optp%optpicl(ish)%abs(chan,3) / aux%dg(lay,prof) + &
                    & optp%optpicl(ish)%abs(chan,4) / aux%dg(lay,prof) ** 2) * &
                    & raytracing%ltick(lay,prof)
                  transmission_scatt_ir%OPDPA(lay,j)          = &
                    & transmission_scatt_ir%OPDPA(lay,j) + transmission_scatt_ir%OPDPACLS(ityp,lay,j)
                  transmission_scatt_ir%OPDPSCLS(ityp,lay,j) = profiles(prof)%cloud(ityp,lay) * &
                    & (optp%optpicl(ish)%sca(chan,1) + optp%optpicl(ish)%sca(chan,2) * aux%dg(lay,prof) + &
                    & optp%optpicl(ish)%sca(chan,3) / aux%dg(lay,prof) + &
                    & optp%optpicl(ish)%sca(chan,4) / aux%dg(lay,prof) ** 2) * &
                    & raytracing%ltick(lay,prof)
                  transmission_scatt_ir%OPDPS(lay,j)          = &
                    & transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSCLS(ityp,lay,j)
                  transmission_scatt_ir%GPARCLS(ityp,lay,j)  =  + (                                     &
                    & optp%optpicl(ish)%bpr(chan,1) + optp%optpicl(ish)%bpr(chan,2) * aux%dg(lay,prof) +  &
                    & optp%optpicl(ish)%bpr(chan,3) * aux%dg(lay,prof) ** 2 +                             &
                    & optp%optpicl(ish)%bpr(chan,4) * aux%dg(lay,prof) ** 3)
                  transmission_scatt_ir%GPARTOT(lay,j)        = transmission_scatt_ir%GPARTOT(lay,j) +      &
                    & transmission_scatt_ir%GPARCLS(ityp,lay,j) * transmission_scatt_ir%OPDPSCLS(ityp,lay,j)
                ELSE
! ATTENTION MODIFIER LE NIVEAU DE TEMPERATURE
                  IF( ish == 3_jpim ) THEN
                    CALL rttov_baran2013_calc_optpar (optp, chan, &
                       & profiles(prof)%t(lay), profiles(prof)%cloud(ityp,lay), &
                       & abso, sca, bpr, asym)
                  ELSEIF( ish == 4_jpim ) THEN
                    CALL rttov_baran2014_calc_optpar (optp, chan, &
                       & profiles(prof)%t(lay), profiles(prof)%cloud(ityp,lay), &
                       & abso, sca, bpr, asym)
                  ENDIF
                  
                  transmission_scatt_ir%OPDPACLS(ityp,lay,j) = abso * raytracing%ltick(lay,prof)
                  transmission_scatt_ir%OPDPA(lay,j)          = &
                     & transmission_scatt_ir%OPDPA(lay,j) + transmission_scatt_ir%OPDPACLS(ityp,lay,j)
                  transmission_scatt_ir%OPDPSCLS(ityp,lay,j) = sca * raytracing%ltick(lay,prof)
                  transmission_scatt_ir%OPDPS(lay,j)          = &
                     & transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSCLS(ityp,lay,j)
                  transmission_scatt_ir%GPARCLS(ityp,lay,j)  = bpr
                  transmission_scatt_ir%GPARTOT(lay,j)        = transmission_scatt_ir%GPARTOT(lay,j) +      &
                    & transmission_scatt_ir%GPARCLS(ityp,lay,j) * transmission_scatt_ir%OPDPSCLS(ityp,lay,j)

                ENDIF
      
              ENDIF

!---------------If solar radiation is present,compute the azimuthally averaged------------
!               value of the phase function for the given value of the viewing------------
!               angle and solar zenith angle.
!-----------------------------------------------------------------------------------------
              IF (solar(j)) THEN
!-----------------If ice clouds are present the phase function for for the current value
!                 of the effective generalized diameter is obtained by linear
!                 interpolation and then the azimuthally averaged value is computed
!----------------------------------------------------------------------------------------
                IF (ityp == 6_jpim) THEN
                  IF( ish .LE. 2_jpim ) THEN
                    chan1 = coef_scatt_ir%icl_pha_index(chan)
                    DO k = 1, coef_scatt_ir%fmv_icl_comp - 1
                      IF (aux%dg(lay,prof) >= coef_scatt_ir%fmv_icl_dg(k, ish) .AND.      &
                        & aux%dg(lay,prof) <= coef_scatt_ir%fmv_icl_dg(k + 1, ish)) THEN
                        deltap     = (optp%optpicl(ish)%pha(chan1, k + 1, :) - optp%optpicl(ish)%pha(chan1, k, :))
                        deltadg    = (coef_scatt_ir%fmv_icl_dg(k + 1, ish) - coef_scatt_ir%fmv_icl_dg(k, ish))
                        phasice(:) = optp%optpicl(ish)%pha(chan1, k, :) +      &
                          & deltap(:) * (aux%dg(lay,prof) - coef_scatt_ir%fmv_icl_dg(k, ish)) / deltadg
                        EXIT
                      ENDIF
                    ENDDO
                  ELSE
                    CALL rttov_baran_calc_phase(asym, phangle, phasice)
                  ENDIF
!-----------------------Average phase function for the upward scattered solar beam----------
                  musat       = 1._jprb / raytracing%pathsat(lay,prof)
                  musun       =  - 1._jprb / raytracing%pathsun(lay,prof)
                  zminphadiff = 1._jprb / (coef_scatt_ir%fmv_icl_ph_val_min * deg2rad)

                  CALL int_phase_fn(musat, musun, 180.0_jprb - relazi, zminphadiff, &
                                    phasice, coef_scatt_ir%fmv_icl_ph_val_cos, &
                                    coef_scatt_ir%ifmv_icl_ph_val, phasint)

                  transmission_scatt_ir%AZPHUPCLS(ityp,lay,j) = phasint
                  transmission_scatt_ir%AZPHUPTOT(lay,j) = transmission_scatt_ir%AZPHUPTOT(lay,j) +      &
                    & transmission_scatt_ir%AZPHUPCLS(ityp,lay,j) * transmission_scatt_ir%OPDPSCLS(ityp,lay,j)
!---------------------Average phase function for the downward scattered solar beam--------
                  musat = -1._jprb * musat

                  CALL int_phase_fn(musat, musun, relazi, zminphadiff,         &
                                    phasice, coef_scatt_ir%fmv_icl_ph_val_cos, &
                                    coef_scatt_ir%ifmv_icl_ph_val, phasint)

                  transmission_scatt_ir%AZPHDOCLS(ityp,lay,j) = phasint
                  transmission_scatt_ir%AZPHDOTOT(lay,j) = transmission_scatt_ir%AZPHDOTOT(lay,j) +      &
                    & transmission_scatt_ir%AZPHDOCLS(ityp,lay,j) * transmission_scatt_ir%OPDPSCLS(ityp,lay,j)
                ELSE
!-----------------------------------------------------------------------------------------
!                     Water clouds
!-----------------------------------------------------------------------------------------
!---------------------Average phase function for the upward scattered solar beam----------
                  chan1 = coef_scatt_ir%wcl_pha_index(chan)
                  musat   = 1._jprb / raytracing%pathsat(lay,prof)
                  musun   =  - 1._jprb / raytracing%pathsun(lay,prof)
                  zminphadiff = 1._jprb / (coef_scatt_ir%fmv_wcl_ph_val_min * deg2rad)

                  CALL int_phase_fn(musat, musun, 180.0_jprb - relazi, zminphadiff, &
                                    optp%optpwcl(ityp)%pha(chan1,1,:),  &
                                    coef_scatt_ir%fmv_wcl_ph_val_cos,  &
                                    coef_scatt_ir%ifmv_wcl_ph_val, phasint)

                  transmission_scatt_ir%AZPHUPCLS(ityp,lay,j) = phasint
                  transmission_scatt_ir%AZPHUPTOT(lay,j) = transmission_scatt_ir%AZPHUPTOT(lay,j) +      &
                    & transmission_scatt_ir%AZPHUPCLS(ityp,lay,j) * transmission_scatt_ir%OPDPSCLS(ityp,lay,j)
!---------------------Average phase function for the downward scattered solar beam--------
                  musat = -1._jprb * musat

                  CALL int_phase_fn(musat, musun, relazi, zminphadiff, &
                                    optp%optpwcl(ityp)%pha(chan1,1,:),  &
                                    coef_scatt_ir%fmv_wcl_ph_val_cos,  &
                                    coef_scatt_ir%ifmv_wcl_ph_val, phasint)

                  transmission_scatt_ir%AZPHDOCLS(ityp,lay,j) = phasint
                  transmission_scatt_ir%AZPHDOTOT(lay,j) = transmission_scatt_ir%AZPHDOTOT(lay,j) +      &
                    & transmission_scatt_ir%AZPHDOCLS(ityp,lay,j) * transmission_scatt_ir%OPDPSCLS(ityp,lay,j)
                ENDIF ! ityp
              ENDIF ! solar(j)
            ENDIF
          ENDDO
        ENDIF ! user_cld_opt_param
        IF (transmission_scatt_ir%OPDPS(lay,j) /= 0._jprb) THEN
          transmission_scatt_ir%GPAR(lay,j) = &
            & transmission_scatt_ir%GPARTOT(lay,j) / transmission_scatt_ir%OPDPS(lay,j)
        ENDIF
        transmission_scatt_ir%OPDPCLDLA(lay,j) = transmission_scatt_ir%OPDPA(lay,j) + &
          & transmission_scatt_ir%GPAR(lay,j) * transmission_scatt_ir%OPDPS(lay,j)
        OPDPCLDL(lay,j)                        = &
          & transmission_scatt_ir%OPDPCLDLA(lay,j) * raytracing%pathsat(lay,prof) * coef%ff_gam(chan)
        IF (solar(j)) THEN
          IF (transmission_scatt_ir%OPDPS(lay,j) /= 0._jprb) THEN
            transmission_scatt_ir%AZPHUP(lay,j) =      &
              & transmission_scatt_ir%AZPHUPTOT(lay,j) / transmission_scatt_ir%OPDPS(lay,j)
            transmission_scatt_ir%AZPHDO(lay,j) =      &
              & transmission_scatt_ir%AZPHDOTOT(lay,j) / transmission_scatt_ir%OPDPS(lay,j)
          ENDIF
          OPDPCLDLSUN(lay,j) = &
            & transmission_scatt_ir%OPDPCLDLA(lay,j) * raytracing%patheff(lay,prof) * coef%ff_gam(chan)
        ENDIF
      ENDDO ! layers
    ENDDO ! channels
  ENDIF ! addclouds

!-----Compute optical parameters for each stream------------------------------------------
  IF (dosolar) transmission_scatt_ir%SSA = 0._jprb
  DO j = 1, nchannels
    chan = chanprof(j)%chan
    prof = chanprof(j)%prof

    ! For layer-specific quantities store just the non-cloudy and cloudy values
    ! When used later the code looks up the appropriate value for each cloud stream
    DO lay = 1, nlayers
      transmission_scatt_ir%OPDPACL(0,lay,j) = OPDPAERL(lay,j)
      IF (opts%rt_ir%addclouds) THEN
        transmission_scatt_ir%OPDPACL(1,lay,j) = OPDPAERL(lay,j) + &
                                                        OPDPCLDL(lay,j)
      ENDIF
    ENDDO

    IF (solar(j)) THEN
      DO lay = 1, nlayers

        transmission_scatt_ir%OPDPACLSUN(0,lay,j) = OPDPAERLSUN(lay,j)
        IF (opts%rt_ir%addclouds) THEN
          transmission_scatt_ir%OPDPACLSUN(1,lay,j) = OPDPAERLSUN(lay,j) + &
                                                      OPDPCLDLSUN(lay,j)
        ENDIF

        transmission_scatt_ir%AZPHACUP(0,lay,j) = transmission_scatt_ir%AZPHAERUPA(lay,j)
        transmission_scatt_ir%AZPHACDO(0,lay,j) = transmission_scatt_ir%AZPHAERDOA(lay,j)
        transmission_scatt_ir%OPDPABS(0,lay,j) = transmission_scatt_ir%OPDPAAER(lay,j) * &
                & raytracing%patheff(lay,prof) * coef%ff_gam(chan)
        transmission_scatt_ir%OPDPSCA(0,lay,j) = transmission_scatt_ir%OPDPSAER(lay,j) * &
                & raytracing%patheff(lay,prof) * coef%ff_gam(chan)

        IF (opts%rt_ir%addclouds) THEN
          IF (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j) /= 0._jprb) THEN
            transmission_scatt_ir%AZPHACUP(1,lay,j) = &
              & transmission_scatt_ir%AZPHAERUPA(lay,j) * transmission_scatt_ir%OPDPSAER(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) + &
              & transmission_scatt_ir%AZPHUP(lay,j) * transmission_scatt_ir%OPDPS(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j))
            transmission_scatt_ir%AZPHACDO(1,lay,j) = &
              & transmission_scatt_ir%AZPHAERDOA(lay,j) * transmission_scatt_ir%OPDPSAER(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) + &
              & transmission_scatt_ir%AZPHDO(lay,j) * transmission_scatt_ir%OPDPS(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j))
          ELSE
            transmission_scatt_ir%AZPHACUP(1,lay,j) = 0._jprb
            transmission_scatt_ir%AZPHACDO(1,lay,j) = 0._jprb
          ENDIF
          transmission_scatt_ir%OPDPABS(1,lay,j) = (transmission_scatt_ir%OPDPAAER(lay,j) + &
            & transmission_scatt_ir%OPDPA(lay,j)) * raytracing%patheff(lay,prof) * coef%ff_gam(chan)
          transmission_scatt_ir%OPDPSCA(1,lay,j) = (transmission_scatt_ir%OPDPSAER(lay,j) + &
            & transmission_scatt_ir%OPDPS(lay,j)) * raytracing%patheff(lay,prof) * coef%ff_gam(chan)
        ENDIF

      ENDDO ! layers
    ENDIF ! solar channel

    DO ist = 0, ircld%nstream(prof)
      IF (ist == 0) THEN
        IF (opts%rt_ir%addaerosl) THEN
          OPD = 0._jprb
          transmission_scatt_ir_stream%OPDPAC(1,ist,j) = 0._jprb
          IF (solar(j)) THEN
            OPDSUN = 0._jprb
            transmission_scatt_ir_stream%OPDPACSUN(1,ist,j) = 0._jprb
          ENDIF
          DO lay = 1, nlayers
            lev = lay + 1
            IF (solar(j)) THEN
              OPDSUN = OPDSUN + transmission_scatt_ir%OPDPACLSUN(0,lay,j)
              transmission_scatt_ir_stream%OPDPACSUN(lev,ist,j)  = OPDSUN
            ENDIF
            OPD = OPD + transmission_scatt_ir%OPDPACL(0,lay,j)
            transmission_scatt_ir_stream%OPDPAC(lev,ist,j)  = OPD
          ENDDO ! layers
        ELSE IF (opts%rt_ir%addclouds) THEN
          transmission_scatt_ir_stream%OPDPAC(:,ist,j)  = 0._jprb
          IF (solar(j)) THEN
            transmission_scatt_ir_stream%OPDPACSUN(:,ist,j)  = 0._jprb
          ENDIF
        ENDIF ! opts%rt_ir%addaerosl
      ELSE
        OPD = 0._jprb
        transmission_scatt_ir_stream%OPDPAC(1,ist,j) = 0._jprb
        IF (solar(j)) THEN
          OPDSUN = 0._jprb
          transmission_scatt_ir_stream%OPDPACSUN(1,ist,j) = 0._jprb
        ENDIF
        DO lay = 1, nlayers
          lev = lay + 1
          isti = ircld%icldarr(ist,lay,prof)
          IF (solar(j)) THEN
            OPDSUN = OPDSUN + transmission_scatt_ir%OPDPACLSUN(isti,lay,j)
            transmission_scatt_ir_stream%OPDPACSUN(lev,ist,j) = OPDSUN
          ENDIF
          OPD = OPD + transmission_scatt_ir%OPDPACL(isti,lay,j)
          transmission_scatt_ir_stream%OPDPAC(lev,ist,j) = OPD
        ENDDO ! layers
      ENDIF ! ist == 0
    ENDDO ! ist
  ENDDO ! channels

  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDPSCATTIR', 1_jpim, ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE int_phase_fn(musat, musun, relazi, zminphadiff, pha, cospha, ipha, phasint)
    ! Interpolate phase function to scattering angle
    REAL(KIND=jprb),    INTENT(IN)  :: musat, musun, relazi
    REAL(KIND=jprb),    INTENT(IN)  :: zminphadiff
    REAL(KIND=jprb),    INTENT(IN)  :: pha(:)
    REAL(KIND=jprb),    INTENT(IN)  :: cospha(:)
    INTEGER(KIND=jpim), INTENT(IN)  :: ipha(:)
    REAL(KIND=jprb),    INTENT(OUT) :: phasint

    INTEGER(KIND=jpim) :: ikk, kk
    REAL(KIND=jprb)    :: ztmpx, scattangle, deltap, delta

    phasint = 0._jprb
    ztmpx   = SQRT((1._jprb - musat ** 2) * (1._jprb - musun ** 2))
    scattangle = musat * musun + ztmpx * COS(relazi * deg2rad)
    ikk        = MAX(1_jpim, INT(ACOS(scattangle) * zminphadiff, jpim))
    kk         = ipha(ikk) - 1_jpim
    deltap     = pha(kk + 1) - pha(kk)
    delta      = cospha(kk) - cospha(kk + 1)
    phasint    = pha(kk) + deltap * (cospha(kk) - scattangle) / delta

  END SUBROUTINE int_phase_fn

END SUBROUTINE rttov_opdpscattir
