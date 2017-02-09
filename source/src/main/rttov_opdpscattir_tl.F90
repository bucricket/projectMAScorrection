!     Compute Optical parameters for aerosols and clouds
SUBROUTINE rttov_opdpscattir_tl( &
            & nlayers,                         &
            & chanprof,                        &
            & opts,                            &
            & aux,                             &
            & aux_tl,                          &
            & profiles,                        &
            & profiles_tl,                     &
            & profiles_dry_tl,                 &
            & aer_opt_param,                   &
            & cld_opt_param,                   &
            & dosolar,                         &
            & solar,                           &
            & coef,                            &
            & coef_scatt_ir,                   &
            & raytracing,                      &
            & raytracing_tl,                   &
            & transmission_scatt_ir,           &
            & transmission_scatt_ir_tl,        &
            & transmission_scatt_ir_stream,    &
            & transmission_scatt_ir_stream_tl, &
            & optp,                            &
            & ircld,                           &
            & ircld_tl)
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
!     1           30/07/2004   Marco Matricardi. ECMWF.
!     1.1         05/02/2007   Removed polarisation R Saunders
!     1.2         15/07/2009   User defined ToA. Layers distinct from levels (P.Rayer)
!     1.3         03/11/2009   Transmittances / optical depths on levels (A Geer)
!     1.4         02/12/2009   Fixed a number of bugs due to the wrong assumption that aerosol/cloud
!                              related quantities are defined on levels (thay are layer
!                              average quantities). Marco Matricardi
!     1.5         02/12/2009   Introduced multiple cloud types in a single layer. Pathsun, Pathsat and
!                              related quantities are now layer arrays (Marco Matricardi)
!     1.6         05/07/2010   Remove addsolar flag from profiles structure (J Hocking)
!     1.7         02/08/2010   Allow variable pressure levels (lgradp) (N Bormann)
!     1.8         14/12/2010   Use traj0_sta%solar array to flag channels for which solar calculations
!                              should be performed (J Hocking)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
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
  TYPE(profile_type              ), INTENT(IN)    :: profiles_tl(SIZE(profiles))
  TYPE(profile_type              ), INTENT(IN)    :: profiles_dry_tl(SIZE(profiles))
  TYPE(rttov_opt_param), OPTIONAL , INTENT(IN)    :: aer_opt_param
  TYPE(rttov_opt_param), OPTIONAL , INTENT(IN)    :: cld_opt_param
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: dosolar
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: solar(SIZE(chanprof))
  TYPE(profile_aux               ), INTENT(IN)    :: aux
  TYPE(profile_aux               ), INTENT(INOUT) :: aux_tl
  TYPE(rttov_coef                ), INTENT(IN)    :: coef
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_tl
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir_stream
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream_tl
  TYPE(rttov_coef_scatt_ir       ), INTENT(IN)    :: coef_scatt_ir
  TYPE(rttov_optpar_ir           ), INTENT(IN)    :: optp
  TYPE(ircld_type),                 INTENT(IN)    :: ircld
  TYPE(ircld_type),                 INTENT(INOUT) :: ircld_tl
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing_tl
!INTF_END
!     End of subroutine arguments
#include "rttov_baran2013_calc_optpar_tl.interface"
#include "rttov_baran2014_calc_optpar_tl.interface"
#include "rttov_baran_calc_phase_tl.interface"

!       Local scalars:
  INTEGER(KIND=jpim) :: j, chan, chan1, i, prof, ish, k, kk, ityp, iae
  INTEGER(KIND=jpim) :: lev, lay, lctyp
  INTEGER(KIND=jpim) :: ist, isti
  REAL   (KIND=jprb) :: OPD, OPDPSUN
  REAL   (KIND=jprb) :: absch_tl
  REAL   (KIND=jprb) :: scach_tl
  REAL   (KIND=jprb) :: bparh_tl
  REAL   (KIND=jprb) :: afac
  REAL   (KIND=jprb) :: sfac
  REAL   (KIND=jprb) :: gfac
  REAL   (KIND=jprb) :: deltaice
  REAL   (KIND=jprb) :: delth
  REAL   (KIND=jprb) :: frach_tl
  REAL   (KIND=jprb) :: phasint_tl
  REAL   (KIND=jprb) :: musat_tl
  REAL   (KIND=jprb) :: musun_tl
  REAL   (KIND=jprb) :: musat
  REAL   (KIND=jprb) :: musun
  REAL   (KIND=jprb) :: frach
  REAL   (KIND=jprb) :: absch
  REAL   (KIND=jprb) :: scach
  REAL   (KIND=jprb) :: bparh
!       Local arrays:
  REAL   (KIND=jprb) :: OPDPAERL   (nlayers,SIZE(chanprof))
  REAL   (KIND=jprb) :: OPDPCLDL   (nlayers,SIZE(chanprof))
  REAL   (KIND=jprb) :: OPDPCLDLSUN(nlayers,SIZE(chanprof))
  REAL   (KIND=jprb) :: OPDPAERLSUN(nlayers,SIZE(chanprof))
  REAL   (KIND=jprb) :: pfac       (coef_scatt_ir%fmv_aer_ph)
  REAL   (KIND=jprb) :: phash      (coef_scatt_ir%fmv_aer_ph)
  REAL   (KIND=jprb) :: phash_tl   (coef_scatt_ir%fmv_aer_ph)
  REAL   (KIND=jprb) :: phasice    (coef_scatt_ir%fmv_icl_ph)
  REAL   (KIND=jprb) :: phasice_tl (coef_scatt_ir%fmv_icl_ph)
  REAL   (KIND=jprb) :: deltap     (coef_scatt_ir%fmv_icl_ph)
  REAL   (KIND=jprb) :: zminphadiff, relazi
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: nchannels
  
  REAL(KIND=jprb)  :: abso, abso_tl
  REAL(KIND=jprb)  :: sca,  sca_tl
  REAL(KIND=jprb)  :: bpr,  bpr_tl
  REAL(KIND=jprb)  :: asym, asym_tl

  
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDPSCATTIR_TL', 0_jpim, ZHOOK_HANDLE)
  nprofiles   = SIZE(profiles)
  nchannels   = SIZE(chanprof)
  OPDPCLDL    = 0._jprb
  OPDPAERL    = 0._jprb
  IF (dosolar) THEN
    OPDPCLDLSUN = 0._jprb
    OPDPAERLSUN = 0._jprb
  ENDIF

!-----------------------------------------------------------------------------------------
!         1.   CALCULATE OPTICAL DEPTHS OF AEROSOLS
!-----------------------------------------------------------------------------------------
!-----Compute relative humidity-----------------------------------------------------------
  IF (opts%rt_ir%addaerosl .AND. .NOT. opts%rt_ir%user_aer_opt_param) THEN
    DO j = 1, nprofiles
      DO lay = 1, nlayers
        lev = lay + 1
        ircld_tl%tave(lay,j)     = (profiles_tl(j)%t(lev - 1) + profiles_tl(j)%t(lev)) / 2._jprb
        ircld_tl%wmixave(lay,j)  = (profiles_dry_tl(j)%q(lev - 1) + profiles_dry_tl(j)%q(lev)) / 2._jprb
        If(opts%interpolation%lgradp) ircld_tl%xpresave(lay,j) = &
                                  & (profiles_tl(j)%p(lev - 1) + profiles_tl(j)%p(lev)) / 2._jprb
!-----------Compute vater vapour partial pressure-----------------------------------------
        ircld_tl%ESW(lay,j)   = ircld_tl%TAVE(lay,j) * ircld%ESW(lay,j) * &
                               & 17.502_jprb * (T00 - 32.19_jprb) / (ircld%TAVE(lay,j) - 32.19_jprb) ** 2
        ircld_tl%ESI(lay,j)   = ircld_tl%TAVE(lay,j) * ircld%ESI(lay,j) * &
                               & 22.587_jprb * (T00 + 0.7_jprb) / (ircld%TAVE(lay,j) + 0.7_jprb) ** 2
        IF (ircld%TAVE(lay,j) > T00) THEN
          ircld_tl%PPV(lay,j) = ircld_tl%ESW(lay,j)! Water phase
        ELSE IF (ircld%TAVE(lay,j) > TI .AND. ircld%TAVE(lay,j) <= T00) THEN
          ircld_tl%PPV(lay,j) =                                                                             &
            & ircld_tl%ESI(lay,j) + ircld_tl%ESW(lay,j) * ((ircld%TAVE(lay,j) - TI) / (T00 - TI)) ** 2 -  &
            & ircld_tl%ESI(lay,j) * ((ircld%TAVE(lay,j) - TI) / (T00 - TI)) ** 2 +                         &
            & ircld_tl%TAVE(lay,j) * (ircld%ESW(lay,j) - ircld%ESI(lay,j)) * 2 *                          &
            & ((ircld%TAVE(lay,j) - TI) / (T00 - TI) ** 2)! Mixed phase
        ELSE IF (ircld%TAVE(lay,j) <= TI) THEN
          ircld_tl%PPV(lay,j) = ircld_tl%ESI(lay,j)! Ice phase
        ENDIF
        ircld_tl%PPV(lay,j)  = ircld_tl%PPV(lay,j) / 100._jprb
!-----------------------------------------------------------------------------------------
! layer average relative humidity
        aux_tl%RELHUM(lay,j) = ircld_tl%wmixave(lay,j) * 100._jprb *                       &
          & 1e-6_jprb * ircld%xpresave(lay,j)                                              &
          &  / (ircld%PPV(lay,j) * (1._jprb + ircld%wmixave(lay,j) * 1e-6_jprb)) -         &
          & ircld_tl%wmixave(lay,j) * 100._jprb * 1e-6_jprb ** 2 * ircld%xpresave(lay,j) * &
          & ircld%wmixave(lay,j) * ircld%PPV(lay,j) /                                      &
          & (ircld%PPV(lay,j) * (1._jprb + ircld%wmixave(lay,j) * 1.e-6_jprb)) ** 2 -      &
          & ircld_tl%PPV(lay,j) * 100._jprb * ircld%wmixave(lay,j) * 1.e-6_jprb *          &
          & ircld%xpresave(lay,j) * (1._jprb + ircld%wmixave(lay,j) * 1.e-6_jprb) /        &
          & (ircld%PPV(lay,j) * (1._jprb + ircld%wmixave(lay,j) * 1e-6_jprb)) ** 2
        IF(opts%interpolation%lgradp) aux_tl%RELHUM(lay,j) = aux_tl%RELHUM(lay,j) +        &
          &   100._jprb * ircld%wmixave(lay,j) * 1.e-6_jprb * ircld_tl%xpresave(lay,j) /   &
          &   (ircld%ppv(lay,j) * (1._jprb + ircld%wmixave(lay,j) * 1.e-6_jprb))
        IF (aux%RELHUMREF(lay,j) > 99._jprb) THEN
          aux_tl%RELHUM(lay,j) = 0._jprb
        ENDIF
      ENDDO ! layers
    ENDDO ! profiles
  ENDIF ! opts%rt_ir%addaerosl
!-----------------------------------------------------------------------------------------
  IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
    transmission_scatt_ir_tl%OPDPAAER = 0._jprb
    transmission_scatt_ir_tl%OPDPSAER = 0._jprb
    IF (dosolar) THEN
      transmission_scatt_ir_tl%AZPHAERUPA = 0._jprb
      transmission_scatt_ir_tl%AZPHAERDOA = 0._jprb
    ENDIF
  ENDIF

  IF (opts%rt_ir%addaerosl) THEN
    transmission_scatt_ir_tl%GPARAERA = 0._jprb
    transmission_scatt_ir_tl%GPARAER  = 0._jprb
    IF (dosolar) THEN
      transmission_scatt_ir_tl%AZPHAERUP = 0._jprb
      transmission_scatt_ir_tl%AZPHAERDO = 0._jprb
    ENDIF

    DO j = 1, nchannels
      chan = chanprof(j)%chan
      prof = chanprof(j)%prof
      relazi  = profiles(prof)%azangle - profiles(prof)%sunazangle
      DO lay = 1, nlayers
        
        IF (opts%rt_ir%user_aer_opt_param) THEN
          ! Input optical parameters are static

          transmission_scatt_ir_tl%OPDPAAER(lay,j) = aer_opt_param%abs(j,lay) * &
                                                      raytracing_tl%ltick(lay,prof)
          transmission_scatt_ir_tl%OPDPSAER(lay,j) = aer_opt_param%sca(j,lay) * &
                                                      raytracing_tl%ltick(lay,prof)
          transmission_scatt_ir_tl%GPARAERA(lay,j) = aer_opt_param%sca(j,lay) * &
                                                      aer_opt_param%bpr(j,lay) * &
                                                      raytracing_tl%ltick(lay,prof)

          IF (solar(j)) THEN
!-------------Average phase function for the upward scattered solar beam------------------
            musat       = 1._jprb / raytracing%pathsat(lay,prof)
            musun       =  - 1._jprb / raytracing%pathsun(lay,prof)
            musat_tl    =  - raytracing_tl%pathsat(lay,prof) / raytracing%pathsat(lay,prof) ** 2
            musun_tl    = raytracing_tl%pathsun(lay,prof) / raytracing%pathsun(lay,prof) ** 2
            zminphadiff = 1._jprb / (aer_opt_param%minphadiff * deg2rad)

            CALL int_phase_fn_tl(musat, musat_tl, musun, musun_tl, 180.0_jprb - relazi, zminphadiff, &
                                 aer_opt_param%pha(j,lay,:), aer_opt_param%cosphangle, &
                                 aer_opt_param%iphangle, phasint_tl)

            transmission_scatt_ir_tl%AZPHAERUP(lay,j) = transmission_scatt_ir_tl%AZPHAERUP(lay,j) + &
              & phasint_tl * aer_opt_param%sca(j,lay) * raytracing%ltick(lay,prof) +    &
              & raytracing_tl%ltick(lay,prof) * aer_opt_param%sca(j,lay) *              &
              & transmission_scatt_ir%PHASINTUPREF(1,lay,j)
!-------------Average phase function for the downward scattered solar beam------
            musat    = -1._jprb * musat
            musat_tl = -1._jprb * musat_tl

            CALL int_phase_fn_tl(musat, musat_tl, musun, musun_tl, relazi, zminphadiff, &
                                 aer_opt_param%pha(j,lay,:), aer_opt_param%cosphangle,  &
                                 aer_opt_param%iphangle, phasint_tl)

            transmission_scatt_ir_tl%AZPHAERDO(lay,j) = transmission_scatt_ir_tl%AZPHAERDO(lay,j) + &
              & phasint_tl * aer_opt_param%sca(j,lay) * raytracing%ltick(lay,prof) +    &
              & raytracing_tl%ltick(lay,prof) * aer_opt_param%sca(j,lay) *              &
              & transmission_scatt_ir%PHASINTDOREF(1,lay,j)
          ENDIF
        ELSE
          DO i = 1, aux%iaernum(lay,prof)
            iae = aux%iaertyp(i,lay,prof)
            IF (coef_scatt_ir%fmv_aer_rh(iae) /= 1) THEN
!---------------Interpolate scattering parameters to actual value of relative humidity------
              DO k = 1, coef_scatt_ir%fmv_aer_rh(iae) - 1
                IF (aux%relhum(lay,prof) >= optp%optpaer(iae)%fmv_aer_rh_val(k) .AND.      &
                  & aux%relhum(lay,prof) <= optp%optpaer(iae)%fmv_aer_rh_val(k + 1)) THEN
                  delth    = (optp%optpaer(iae)%fmv_aer_rh_val(K + 1) - optp%optpaer(iae)%fmv_aer_rh_val(K))
                  frach_tl = aux_tl%relhum(lay,prof)
                  afac     = (optp%optpaer(iae)%abs(chan,k + 1) - optp%optpaer(iae)%abs(chan,k)) / delth
                  sfac     = (optp%optpaer(iae)%sca(chan,k + 1) - optp%optpaer(iae)%sca(chan,k)) / delth
                  gfac     = (optp%optpaer(iae)%bpr(chan,k + 1) - optp%optpaer(iae)%bpr(chan,k)) / delth
                  frach    = (aux%relhum(lay,prof) - optp%optpaer(iae)%fmv_aer_rh_val(k))
                  absch    = optp%optpaer(iae)%abs(chan,k) + afac * frach
                  scach    = optp%optpaer(iae)%sca(chan,k) + sfac * frach
                  bparh    = optp%optpaer(iae)%bpr(chan,k) + gfac * frach
                  absch_tl = afac * frach_tl
                  scach_tl = sfac * frach_tl
                  bparh_tl = gfac * frach_tl
                  IF (solar(j)) THEN
                    chan1 = coef_scatt_ir%aer_pha_index(chan)
                    pfac(1:coef_scatt_ir%fmv_aer_ph) = (                                  &
                      & optp%optpaer(iae)%pha(chan1,k + 1,1:coef_scatt_ir%fmv_aer_ph) -  &
                      & optp%optpaer(iae)%pha(chan1,k,1:coef_scatt_ir%fmv_aer_ph)) / delth
                    phash(1:coef_scatt_ir%fmv_aer_ph) = &
                      & optp%optpaer(iae)%pha(chan1,k,1:coef_scatt_ir%fmv_aer_ph) + &
                      & pfac(1:coef_scatt_ir%fmv_aer_ph) * frach
                    phash_tl(1:coef_scatt_ir%fmv_aer_ph) = pfac(1:coef_scatt_ir%fmv_aer_ph) * frach_tl
                  ENDIF
                  EXIT
                ENDIF
              ENDDO
            ELSE
              absch = optp%optpaer(iae)%abs(chan,1)
              scach = optp%optpaer(iae)%sca(chan,1)
              bparh = optp%optpaer(iae)%bpr(chan,1)
              absch_tl = 0._jprb
              scach_tl = 0._jprb
              bparh_tl = 0._jprb
              IF (solar(j)) THEN
                chan1 = coef_scatt_ir%aer_pha_index(chan)
                phash(1:coef_scatt_ir%fmv_aer_ph) = optp%optpaer(iae)%pha(chan1, 1, 1:coef_scatt_ir%fmv_aer_ph)
                phash_tl(1:coef_scatt_ir%fmv_aer_ph) = 0._jprb
              ENDIF
            ENDIF
!-------------Compute optical parameters considering the contribution of------------------
!             all the aerosol components present in the layer
            transmission_scatt_ir_tl%OPDPAAER(lay,j) = transmission_scatt_ir_tl%OPDPAAER(lay,j) +              &
              & profiles_tl(prof)%aerosols(iae,lay) * absch * raytracing%ltick(lay,prof) +         &
              & absch_tl * profiles(prof)%aerosols(iae,lay) * raytracing%ltick(lay,prof) +         &
              & raytracing_tl%ltick(lay,prof) * absch * profiles(prof)%aerosols(iae,lay)
            transmission_scatt_ir_tl%OPDPSAER(lay,j) = transmission_scatt_ir_tl%OPDPSAER(lay,j) +              &
              & profiles_tl(prof)%aerosols(iae,lay) * scach * raytracing%ltick(lay,prof) +         &
              & scach_tl * profiles(prof)%aerosols(iae,lay) * raytracing%ltick(lay,prof) +         &
              & raytracing_tl%ltick(lay,prof) * scach * profiles(prof)%aerosols(iae,lay)
            transmission_scatt_ir_tl%GPARAERA(lay,j) = transmission_scatt_ir_tl%GPARAERA(lay,j) +              &
              & profiles_tl(prof)%aerosols(iae,lay) * scach * bparh *                                            &
              & raytracing%ltick(lay,prof) +                                                       &
              & scach_tl * profiles(prof)%aerosols(iae,lay) * raytracing%ltick(lay,prof) * bparh + &
              & raytracing_tl%ltick(lay,prof) * scach *                                            &
              & profiles(prof)%aerosols(iae,lay) * bparh +                                                       &
              & bparh_tl * profiles(prof)%aerosols(iae,lay) * raytracing%ltick(lay,prof) * scach
!-------------If solar radiation is present,compute the azimuthally averaged--------------
!             value of the phase function for the given value of the viewing--------------
!             angle and solar zenith angle.
            IF (solar(j)) THEN
              chan1 = coef_scatt_ir%aer_pha_index(chan)
!-------------Average phase function for the upward scattered solar beam------------------
              musat       = 1._jprb / raytracing%pathsat(lay,prof)
              musun       =  - 1._jprb / raytracing%pathsun(lay,prof)
              musat_tl    =  - raytracing_tl%pathsat(lay,prof) / raytracing%pathsat(lay,prof) ** 2
              musun_tl    = raytracing_tl%pathsun(lay,prof) / raytracing%pathsun(lay,prof) ** 2
              zminphadiff = 1._jprb / (coef_scatt_ir%fmv_aer_ph_val_min * deg2rad)

              CALL int_phase_fn_tl(musat, musat_tl, musun, musun_tl, 180.0_jprb - relazi, zminphadiff, &
                                   phash, coef_scatt_ir%fmv_aer_ph_val_cos, &
                                   coef_scatt_ir%ifmv_aer_ph_val, phasint_tl, phash_TL)

              transmission_scatt_ir_tl%AZPHAERUP(lay,j) = transmission_scatt_ir_tl%AZPHAERUP(lay,j) +            &
                & profiles_tl(prof)%aerosols(iae,lay) * transmission_scatt_ir%PHASINTUPREF(i,lay,j) * scach *    &
                & raytracing%ltick(lay,prof) +                                                       &
                & phasint_tl * profiles(prof)%aerosols(iae,lay) * scach *                                          &
                & raytracing%ltick(lay,prof) +                                                       &
                & scach_tl * profiles(prof)%aerosols(iae,lay) * transmission_scatt_ir%PHASINTUPREF(i,lay,j) *    &
                & raytracing%ltick(lay,prof) +                                                       &
                & raytracing_tl%ltick(lay,prof) * scach * profiles(prof)%aerosols(iae,lay) *         &
                & transmission_scatt_ir%PHASINTUPREF(i,lay,j)
!-------------Average phase function for the downward scattered solar beam------
              musat    = -1._jprb * musat
              musat_tl = -1._jprb * musat_tl

              CALL int_phase_fn_tl(musat, musat_tl, musun, musun_tl, relazi, zminphadiff, &
                                   phash, coef_scatt_ir%fmv_aer_ph_val_cos,               &
                                   coef_scatt_ir%ifmv_aer_ph_val, phasint_tl, phash_TL)

              transmission_scatt_ir_tl%AZPHAERDO(lay,j) = transmission_scatt_ir_tl%AZPHAERDO(lay,j) +            &
                & profiles_tl(prof)%aerosols(iae,lay) * transmission_scatt_ir%PHASINTDOREF(i,lay,j) * scach *    &
                & raytracing%ltick(lay,prof) +                                                       &
                & phasint_tl * profiles(prof)%aerosols(iae,lay) * scach *                                          &
                & raytracing%ltick(lay,prof) +                                                       &
                & scach_tl * profiles(prof)%aerosols(iae,lay) * transmission_scatt_ir%PHASINTDOREF(i,lay,j) *    &
                & raytracing%ltick(lay,prof) +                                                       &
                & raytracing_tl%ltick(lay,prof) * scach * profiles(prof)%aerosols(iae,lay) *         &
                & transmission_scatt_ir%PHASINTDOREF(i,lay,j)
            ENDIF
          ENDDO ! iaernum
        ENDIF ! opts%rt_ir%user_aer_opt_param
      ENDDO ! layers
!---------Compute final values for optical parameters-------------------------------------
      DO lay = 1, nlayers
        IF (transmission_scatt_ir%OPDPSAER(lay,j) /= 0._jprb) THEN
          transmission_scatt_ir_tl%GPARAER(lay,j) =                                                &
            & transmission_scatt_ir_tl%GPARAERA(lay,j) / transmission_scatt_ir%OPDPSAER(lay,j) -  &
            & transmission_scatt_ir_tl%OPDPSAER(lay,j) * transmission_scatt_ir%GPARAERA(lay,j) /  &
            & transmission_scatt_ir%OPDPSAER(lay,j) ** 2
          IF (solar(j)) THEN
            transmission_scatt_ir_tl%AZPHAERUPA(lay,j) =                                              &
              & transmission_scatt_ir_tl%AZPHAERUP(lay,j) / transmission_scatt_ir%OPDPSAER(lay,j) -  &
              & transmission_scatt_ir_tl%OPDPSAER(lay,j) * transmission_scatt_ir%AZPHAERUP(lay,j) /  &
              & transmission_scatt_ir%OPDPSAER(lay,j) ** 2
            transmission_scatt_ir_tl%AZPHAERDOA(lay,j) =                                              &
              & transmission_scatt_ir_tl%AZPHAERDO(lay,j) / transmission_scatt_ir%OPDPSAER(lay,j) -  &
              & transmission_scatt_ir_tl%OPDPSAER(lay,j) * transmission_scatt_ir%AZPHAERDO(lay,j) /  &
              & transmission_scatt_ir%OPDPSAER(lay,j) ** 2
          ENDIF
        ENDIF
        transmission_scatt_ir_tl%OPDPAERLA(lay,j) = transmission_scatt_ir_tl%OPDPAAER(lay,j) +      &
          & transmission_scatt_ir_tl%OPDPSAER(lay,j) * transmission_scatt_ir%GPARAER(lay,j) +       &
          & transmission_scatt_ir_tl%GPARAER(lay,j) * transmission_scatt_ir%OPDPSAER(lay,j)
        OPDPAERL(lay,j)                           =                                                        &
          & transmission_scatt_ir_tl%OPDPAERLA(lay,j) * raytracing%pathsat(lay,prof) * coef%ff_gam(chan) +  &
          & raytracing_tl%pathsat(lay,prof) * transmission_scatt_ir%OPDPAERLA(lay,j) * coef%ff_gam(chan)
        IF (solar(j)) THEN
          OPDPAERLSUN(lay,j) = &
            & transmission_scatt_ir_tl%OPDPAERLA(lay,j) * raytracing%patheff(lay,prof) * coef%ff_gam(chan) + &
            & raytracing_tl%patheff(lay,prof) * transmission_scatt_ir%OPDPAERLA(lay,j) * coef%ff_gam(chan)
        ENDIF
      ENDDO ! layers
    ENDDO ! channels
  ENDIF ! opts%rt_ir%addaerosl
!-------------------------------------------------------------------------------
!         2.   CALCULATE OPTICAL DEPTHS OF CLOUDS
!-------------------------------------------------------------------------------
  IF (opts%rt_ir%addclouds) THEN
    transmission_scatt_ir_tl%OPDPA   = 0._jprb
    transmission_scatt_ir_tl%OPDPS   = 0._jprb
    transmission_scatt_ir_tl%GPAR    = 0._jprb
    transmission_scatt_ir_tl%GPARTOT = 0._jprb
    IF (dosolar) THEN
      transmission_scatt_ir_tl%AZPHUPTOT = 0._jprb
      transmission_scatt_ir_tl%AZPHDOTOT = 0._jprb
      transmission_scatt_ir_tl%AZPHUP = 0._jprb
      transmission_scatt_ir_tl%AZPHDO = 0._jprb
    ENDIF

    DO j = 1, nchannels
      chan = chanprof(j)%chan
      prof = chanprof(j)%prof
      ish = profiles(prof)%ish
      relazi  = profiles(prof)%azangle - profiles(prof)%sunazangle
      DO lay = 1, nlayers
        IF (solar(j)) THEN
        ENDIF

        IF (opts%rt_ir%user_cld_opt_param) THEN
          ! Input optical parameters are static

          transmission_scatt_ir_tl%OPDPA(lay,j) = cld_opt_param%abs(j,lay) * &
                                                & raytracing_tl%ltick(lay,prof)
          transmission_scatt_ir_tl%OPDPS(lay,j) = cld_opt_param%sca(j,lay) * &
                                                & raytracing_tl%ltick(lay,prof)
          transmission_scatt_ir_tl%GPARTOT(lay,j) = cld_opt_param%bpr(j,lay) * cld_opt_param%sca(j,lay) * &
                                                  & raytracing_tl%ltick(lay,prof)
          IF (solar(j)) THEN
!-------------Average phase function for the upward scattered solar beam------------------
            musat       = 1._jprb / raytracing%pathsat(lay,prof)
            musun       =  - 1._jprb / raytracing%pathsun(lay,prof)
            musat_tl    =  - raytracing_tl%pathsat(lay,prof) / raytracing%pathsat(lay,prof) ** 2
            musun_tl    = raytracing_tl%pathsun(lay,prof) / raytracing%pathsun(lay,prof) ** 2
            zminphadiff = 1._jprb / (cld_opt_param%minphadiff * deg2rad)

            CALL int_phase_fn_tl(musat, musat_tl, musun, musun_tl, 180.0_jprb - relazi, zminphadiff, &
                                 cld_opt_param%pha(j,lay,:), cld_opt_param%cosphangle, &
                                 cld_opt_param%iphangle, phasint_tl)

            ! AZPHUPCLS(:,:,1) is the direct PHASINT
            transmission_scatt_ir_tl%AZPHUPTOT(lay,j) = transmission_scatt_ir_tl%AZPHUPTOT(lay,j) +      &
                                                       & phasint_tl * transmission_scatt_ir%OPDPS(lay,j) + &
                                                       & transmission_scatt_ir%AZPHUPCLS(1,lay,j) *      &
                                                       & transmission_scatt_ir_tl%OPDPS(lay,j)
!-------------Average phase function for the downward scattered solar beam------
            musat    = -1._jprb * musat
            musat_tl = -1._jprb * musat_tl

            CALL int_phase_fn_tl(musat, musat_tl, musun, musun_tl, relazi, zminphadiff, &
                                 cld_opt_param%pha(j,lay,:), cld_opt_param%cosphangle,  &
                                 cld_opt_param%iphangle, phasint_tl)

            ! AZPHDOCLS(:,:,1) is the direct PHASINT
            transmission_scatt_ir_tl%AZPHDOTOT(lay,j) = transmission_scatt_ir_tl%AZPHDOTOT(lay,j) +      &
                                                       & phasint_tl * transmission_scatt_ir%OPDPS(lay,j) + &
                                                       & transmission_scatt_ir%AZPHDOCLS(1,lay,j) *      &
                                                       & transmission_scatt_ir_tl%OPDPS(lay,j)
          ENDIF
        ELSE
          DO lctyp = 1, ncldtyp
            IF (ircld%cldtyp(lctyp, lay,prof) /= 0) THEN
              ityp  = ircld%cldtyp(lctyp, lay,prof)
!-----------------Compute cloud  optical parameters ----------------------------------------
              IF (ityp <= 5_jpim) THEN
!-----------------------------------------------------------------------------------------
!                   For water clouds use stored optical parameters
!-----------------------------------------------------------------------------------------
                transmission_scatt_ir_tl%OPDPACLS(ityp,lay,j) =                                             &
                  & profiles_tl(prof)%cloud(ityp,lay) * coef_scatt_ir%confac(ityp) *                         &
                  & optp%optpwcl(ityp)%abs(chan,1) * raytracing%ltick(lay,prof) +               &
                  & raytracing_tl%ltick(lay,prof) * profiles(prof)%cloud(ityp,lay) *            &
                  & coef_scatt_ir%confac(ityp) * optp%optpwcl(ityp)%abs(chan,1)
                transmission_scatt_ir_tl%OPDPA(lay,j)          =      &
                  & transmission_scatt_ir_tl%OPDPA(lay,j) + transmission_scatt_ir_tl%OPDPACLS(ityp,lay,j)
                transmission_scatt_ir_tl%OPDPSCLS(ityp,lay,j) =                                             &
                  & profiles_tl(prof)%cloud(ityp,lay) * coef_scatt_ir%confac(ityp) *                         &
                  & optp%optpwcl(ityp)%sca(chan,1) * raytracing%ltick(lay,prof) +               &
                  & raytracing_tl%ltick(lay,prof) * profiles(prof)%cloud(ityp,lay) *            &
                  & coef_scatt_ir%confac(ityp) * optp%optpwcl(ityp)%sca(chan,1)
                transmission_scatt_ir_tl%OPDPS(lay,j)          =      &
                  & transmission_scatt_ir_tl%OPDPS(lay,j) + transmission_scatt_ir_tl%OPDPSCLS(ityp,lay,j)
                transmission_scatt_ir_tl%GPARCLS(ityp,lay,j)  = 0._jprb
                transmission_scatt_ir_tl%GPARTOT(lay,j)        = transmission_scatt_ir_tl%GPARTOT(lay,j) +        &
                  & transmission_scatt_ir_tl%GPARCLS(ityp,lay,j) * transmission_scatt_ir%OPDPSCLS(ityp,lay,j) + &
                  & transmission_scatt_ir%GPARCLS(ityp,lay,j) * transmission_scatt_ir_tl%OPDPSCLS(ityp,lay,j)
              ELSE
!-----------------------------------------------------------------------------------------
!                   For ice clouds optical parameters are computed using regression
!                   coefficients
!-----------------------------------------------------------------------------------------
                IF( ish .LE. 2_jpim ) THEN
                  transmission_scatt_ir_tl%OPDPACLS(ityp,lay,j) = profiles_tl(prof)%cloud(ityp,lay) * (     &
                    & optp%optpicl(ish)%abs(chan,1) + optp%optpicl(ish)%abs(chan,2) * aux%dg(lay,prof) +      &
                    & optp%optpicl(ish)%abs(chan,3) / aux%dg(lay,prof) +                                      &
                    & optp%optpicl(ish)%abs(chan,4) / aux%dg(lay,prof) ** 2) *                                &
                    & raytracing%ltick(lay,prof) +                                              &
                    & aux_tl%dg(lay,prof) * optp%optpicl(ish)%abs(chan,2) * profiles(prof)%cloud(ityp,lay) *  &
                    & raytracing%ltick(lay,prof) -                                              &
                    & aux_tl%dg(lay,prof) * optp%optpicl(ish)%abs(chan,3) * profiles(prof)%cloud(ityp,lay) *  &
                    & raytracing%ltick(lay,prof) / aux%dg(lay,prof) ** 2 -                      &
                    & aux_tl%dg(lay,prof) * optp%optpicl(ish)%abs(chan,4) * profiles(prof)%cloud(ityp,lay) *  &
                    & raytracing%ltick(lay,prof) * 2._jprb / aux%dg(lay,prof) ** 3 +            &
                    & raytracing_tl%ltick(lay,prof) * profiles(prof)%cloud(ityp,lay) * (        &
                    & optp%optpicl(ish)%abs(chan,1) + optp%optpicl(ish)%abs(chan,2) * aux%dg(lay,prof) +      &
                    & optp%optpicl(ish)%abs(chan,3) / aux%dg(lay,prof) +                                      &
                    & optp%optpicl(ish)%abs(chan,4) / aux%dg(lay,prof) ** 2)
                  transmission_scatt_ir_tl%OPDPA(lay,j)          =      &
                    & transmission_scatt_ir_tl%OPDPA(lay,j) + transmission_scatt_ir_tl%OPDPACLS(ityp,lay,j)
                  transmission_scatt_ir_tl%OPDPSCLS(ityp,lay,j) = profiles_tl(prof)%cloud(ityp,lay) * (     &
                    & optp%optpicl(ish)%sca(chan,1) + optp%optpicl(ish)%sca(chan,2) * aux%dg(lay,prof) +      &
                    & optp%optpicl(ish)%sca(chan,3) / aux%dg(lay,prof) +                                      &
                    & optp%optpicl(ish)%sca(chan,4) / aux%dg(lay,prof) ** 2) *                                &
                    & raytracing%ltick(lay,prof) +                                              &
                    & aux_tl%dg(lay,prof) * optp%optpicl(ish)%sca(chan,2) * profiles(prof)%cloud(ityp,lay) *  &
                    & raytracing%ltick(lay,prof) -                                              &
                    & aux_tl%dg(lay,prof) * optp%optpicl(ish)%sca(chan,3) * profiles(prof)%cloud(ityp,lay) *  &
                    & raytracing%ltick(lay,prof) / aux%dg(lay,prof) ** 2 -                      &
                    & aux_tl%dg(lay,prof) * optp%optpicl(ish)%sca(chan,4) * profiles(prof)%cloud(ityp,lay) *  &
                    & raytracing%ltick(lay,prof) * 2._jprb / aux%dg(lay,prof) ** 3 +            &
                    & raytracing_tl%ltick(lay,prof) * profiles(prof)%cloud(ityp,lay) * (        &
                    & optp%optpicl(ish)%sca(chan,1) + optp%optpicl(ish)%sca(chan,2) * aux%dg(lay,prof) +      &
                    & optp%optpicl(ish)%sca(chan,3) / aux%dg(lay,prof) +                                      &
                    & optp%optpicl(ish)%sca(chan,4) / aux%dg(lay,prof) ** 2)
                  transmission_scatt_ir_tl%OPDPS(lay,j)          =      &
                    & transmission_scatt_ir_tl%OPDPS(lay,j) + transmission_scatt_ir_tl%OPDPSCLS(ityp,lay,j)
                  transmission_scatt_ir_tl%GPARCLS(ityp,lay,j)  =      &
                    & aux_tl%dg(lay,prof) * optp%optpicl(ish)%bpr(chan,2) +                               &
                    & aux_tl%dg(lay,prof) * 2._jprb * aux%dg(lay,prof) * optp%optpicl(ish)%bpr(chan,3) +  &
                    & aux_tl%dg(lay,prof) * 3._jprb * aux%dg(lay,prof) ** 2 * optp%optpicl(ish)%bpr(chan,4)
                  transmission_scatt_ir_tl%GPARTOT(lay,j)        = transmission_scatt_ir_tl%GPARTOT(lay,j) + &
                    & transmission_scatt_ir_tl%GPARCLS(ityp,lay,j) * &
                    & transmission_scatt_ir%OPDPSCLS(ityp,lay,j) +   &
                    & transmission_scatt_ir%GPARCLS(ityp,lay,j) *    &
                    & transmission_scatt_ir_tl%OPDPSCLS(ityp,lay,j)
                ELSE
                  ! ATTENTION MODIFIER LE NIVEAU DE TEMPERATURE
                  IF (ish == 3_JPIM) THEN 
                    CALL rttov_baran2013_calc_optpar_tl(optp, chan, &
                       & profiles(prof)%t(lay), profiles(prof)%cloud(ityp,lay), &
                       & profiles_tl(prof)%t(lay), profiles_tl(prof)%cloud(ityp,lay), &
                       & abso, sca, bpr, asym, &
                       & abso_tl, sca_tl, bpr_tl, asym_tl)
                  ELSEIF (ish == 4_JPIM) THEN 
                    CALL rttov_baran2014_calc_optpar_tl(optp, chan, &
                       & profiles(prof)%t(lay), profiles(prof)%cloud(ityp,lay), &
                       & profiles_tl(prof)%t(lay), profiles_tl(prof)%cloud(ityp,lay), &
                       & abso, sca, bpr, asym, &
                       & abso_tl, sca_tl, bpr_tl, asym_tl)
                  ENDIF
                  transmission_scatt_ir_tl%OPDPACLS(ityp,lay,j) = &
                           & abso_tl * raytracing%ltick(lay,prof) + &
                           & abso    * raytracing_tl%ltick(lay,prof)
                  transmission_scatt_ir_tl%OPDPA(lay,j)          =      &
                  & transmission_scatt_ir_tl%OPDPA(lay,j) + transmission_scatt_ir_tl%OPDPACLS(ityp,lay,j)
                  transmission_scatt_ir_tl%OPDPSCLS(ityp,lay,j) = &
                           & sca_tl * raytracing%ltick(lay,prof) + &
                           & sca    * raytracing_tl%ltick(lay,prof)
                  transmission_scatt_ir_tl%OPDPS(lay,j)          =      &
                  & transmission_scatt_ir_tl%OPDPS(lay,j) + transmission_scatt_ir_tl%OPDPSCLS(ityp,lay,j)
                  transmission_scatt_ir_tl%GPARCLS(ityp,lay,j)  = bpr_tl
                  transmission_scatt_ir_tl%GPARTOT(lay,j)        = transmission_scatt_ir_tl%GPARTOT(lay,j) +      &
                  & transmission_scatt_ir_tl%GPARCLS(ityp,lay,j) * transmission_scatt_ir%OPDPSCLS(ityp,lay,j) + &
                  & transmission_scatt_ir%GPARCLS(ityp,lay,j) * transmission_scatt_ir_tl%OPDPSCLS(ityp,lay,j)

                ENDIF
              ENDIF !ice/water
!-------------------If solar radiation is present,compute the azimuthally averaged----
!                   value of the phase function for the given value of the viewing----
!                   angle and solar zenith angle.
              IF (solar(j)) THEN
!-------------------If ice clouds are present the phase function for for the current value
!                   of the effective generalized diameter is obtained by linear
!                   interpolation and then the azimuthally averaged value is computed
!----------------------------------------------------------------------------------------
                IF (ityp == 6_jpim) THEN
                  IF( ish .LE. 2_jpim ) THEN
                    chan1 = coef_scatt_ir%icl_pha_index(chan)
                    DO k = 1, coef_scatt_ir%fmv_icl_comp - 1
                      IF (aux%dg(lay,prof) >= coef_scatt_ir%fmv_icl_dg(k, ish) .AND.      &
                        & aux%dg(lay,prof) <= coef_scatt_ir%fmv_icl_dg(k + 1, ish)) THEN
                        deltap       = (optp%optpicl(ish)%pha(chan1, k + 1, :) - optp%optpicl(ish)%pha(chan1, k, :))
                        deltaice     = (coef_scatt_ir%fmv_icl_dg(k + 1, ish) - coef_scatt_ir%fmv_icl_dg(k, ish))
                        phasice(:)   = optp%optpicl(ish)%pha(chan1, k, :) +      &
                          & deltap(:) * (aux%dg(lay,prof) - coef_scatt_ir%fmv_icl_dg(k, ish)) / deltaice
                        DO kk = 1, coef_scatt_ir%fmv_icl_ph
                          phasice_tl(kk) = deltap(kk) * aux_tl%dg(lay,prof) / deltaice
                        ENDDO
                        EXIT
                      ENDIF
                    ENDDO
                  ELSE
                    CALL rttov_baran_calc_phase_tl(asym, asym_tl, phangle, phasice, phasice_tl)
                  ENDIF
!-----------------------Average phase function for the upward scattered solar beam----------
                  musat       = 1._jprb / raytracing%pathsat(lay,prof)
                  musun       =  - 1._jprb / raytracing%pathsun(lay,prof)
                  musat_tl    =  - raytracing_tl%pathsat(lay,prof) / raytracing%pathsat(lay,prof) ** 2
                  musun_tl    = raytracing_tl%pathsun(lay,prof) / raytracing%pathsun(lay,prof) ** 2
                  zminphadiff = 1._jprb / (coef_scatt_ir%fmv_icl_ph_val_min * deg2rad)

                  CALL int_phase_fn_tl(musat, musat_tl, musun, musun_tl, 180.0_jprb - relazi, zminphadiff, &
                                       phasice, coef_scatt_ir%fmv_icl_ph_val_cos, &
                                       coef_scatt_ir%ifmv_icl_ph_val, phasint_tl, phasice_TL)

                  transmission_scatt_ir_tl%AZPHUPCLS(ityp,lay,j) = phasint_tl
                  transmission_scatt_ir_tl%AZPHUPTOT(lay,j) = transmission_scatt_ir_tl%AZPHUPTOT(lay,j) +       &
                                                             & transmission_scatt_ir%AZPHUPCLS(ityp,lay,j) *    &
                                                             & transmission_scatt_ir_tl%OPDPSCLS(ityp,lay,j) +  &
                                                             & transmission_scatt_ir_tl%AZPHUPCLS(ityp,lay,j) * &
                                                             & transmission_scatt_ir%OPDPSCLS(ityp,lay,j)

!---------------------Average phase function for the downward scattered solar beam--------
                  musat    = -1._jprb * musat
                  musat_tl = -1._jprb * musat_tl

                  CALL int_phase_fn_tl(musat, musat_tl, musun, musun_tl, relazi, zminphadiff, &
                                       phasice, coef_scatt_ir%fmv_icl_ph_val_cos,             &
                                       coef_scatt_ir%ifmv_icl_ph_val, phasint_tl, phasice_TL)

                  transmission_scatt_ir_tl%AZPHDOCLS(ityp,lay,j) = phasint_tl
                  transmission_scatt_ir_tl%AZPHDOTOT(lay,j) = transmission_scatt_ir_tl%AZPHDOTOT(lay,j) +       &
                                                             & transmission_scatt_ir%AZPHDOCLS(ityp,lay,j) *    &
                                                             & transmission_scatt_ir_tl%OPDPSCLS(ityp,lay,j) +  &
                                                             & transmission_scatt_ir_tl%AZPHDOCLS(ityp,lay,j) * &
                                                             & transmission_scatt_ir%OPDPSCLS(ityp,lay,j)
                ELSE
!-----------------------------------------------------------------------------------------
!                       Water clouds
!-----------------------------------------------------------------------------------------
                  chan1 = coef_scatt_ir%wcl_pha_index(chan)
!-----------------Average phase function for the upward scattered solar beam--------------
                  musat       = 1._jprb / raytracing%pathsat(lay,prof)
                  musun       =  - 1._jprb / raytracing%pathsun(lay,prof)
                  musat_tl    =  - raytracing_tl%pathsat(lay,prof) / raytracing%pathsat(lay,prof) ** 2
                  musun_tl    = raytracing_tl%pathsun(lay,prof) / raytracing%pathsun(lay,prof) ** 2
                  zminphadiff = 1._jprb / (coef_scatt_ir%fmv_wcl_ph_val_min * deg2rad)

                  CALL int_phase_fn_tl(musat, musat_tl, musun, musun_tl, 180.0_jprb - relazi, zminphadiff, &
                                       optp%optpwcl(ityp)%pha(chan1,1,:), coef_scatt_ir%fmv_wcl_ph_val_cos, &
                                       coef_scatt_ir%ifmv_wcl_ph_val, phasint_tl)

                  transmission_scatt_ir_tl%AZPHUPCLS(ityp,lay,j) = phasint_tl
                  transmission_scatt_ir_tl%AZPHUPTOT(lay,j) = transmission_scatt_ir_tl%AZPHUPTOT(lay,j) +       &
                                                             & transmission_scatt_ir%AZPHUPCLS(ityp,lay,j) *    &
                                                             & transmission_scatt_ir_tl%OPDPSCLS(ityp,lay,j) +  &
                                                             & transmission_scatt_ir_tl%AZPHUPCLS(ityp,lay,j) * &
                                                             & transmission_scatt_ir%OPDPSCLS(ityp,lay,j)

!-----------------------Average phase function for the downward scattered solar beam------------
                  musat    = -1._jprb * musat
                  musat_tl = -1._jprb * musat_tl

                  CALL int_phase_fn_tl(musat, musat_tl, musun, musun_tl, relazi, zminphadiff,              &
                                       optp%optpwcl(ityp)%pha(chan1,1,:), coef_scatt_ir%fmv_wcl_ph_val_cos, &
                                       coef_scatt_ir%ifmv_wcl_ph_val, phasint_tl)

                  transmission_scatt_ir_tl%AZPHDOCLS(ityp,lay,j) = phasint_tl
                  transmission_scatt_ir_tl%AZPHDOTOT(lay,j) = transmission_scatt_ir_tl%AZPHDOTOT(lay,j) +       &
                                                             & transmission_scatt_ir%AZPHDOCLS(ityp,lay,j) *    &
                                                             & transmission_scatt_ir_tl%OPDPSCLS(ityp,lay,j) +  &
                                                             & transmission_scatt_ir_tl%AZPHDOCLS(ityp,lay,j) * &
                                                             & transmission_scatt_ir%OPDPSCLS(ityp,lay,j)
                ENDIF !ityp
              ENDIF !solar
            ENDIF !cloud present
          ENDDO !ncldtyp
        ENDIF ! user_cld_opt_param
        IF (transmission_scatt_ir%OPDPS(lay,j) /= 0._jprb) THEN
          transmission_scatt_ir_tl%GPAR(lay,j) =                                               &
            & transmission_scatt_ir_tl%GPARTOT(lay,j) / transmission_scatt_ir%OPDPS(lay,j) -  &
            & transmission_scatt_ir_tl%OPDPS(lay,j) * transmission_scatt_ir%GPARTOT(lay,j) /  &
            & (transmission_scatt_ir%OPDPS(lay,j) ** 2)
        ENDIF
        transmission_scatt_ir_tl%OPDPCLDLA(lay,j) = transmission_scatt_ir_tl%OPDPA(lay,j) +      &
          & transmission_scatt_ir_tl%OPDPS(lay,j) * transmission_scatt_ir%GPAR(lay,j) +          &
          & transmission_scatt_ir_tl%GPAR(lay,j) * transmission_scatt_ir%OPDPS(lay,j)
        OPDPCLDL(lay,j)                           =                                                        &
          & transmission_scatt_ir_tl%OPDPCLDLA(lay,j) * raytracing%pathsat(lay,prof) * coef%ff_gam(chan) +  &
          & raytracing_tl%pathsat(lay,prof) * transmission_scatt_ir%OPDPCLDLA(lay,j) * coef%ff_gam(chan)
        IF (solar(j)) THEN
          IF (transmission_scatt_ir%OPDPS(lay,j) /= 0._jprb) THEN
            transmission_scatt_ir_tl%AZPHUP(lay,j) =                                               &
              & transmission_scatt_ir_tl%AZPHUPTOT(lay,j) / transmission_scatt_ir%OPDPS(lay,j) -  &
              & transmission_scatt_ir%AZPHUPTOT(lay,j) * transmission_scatt_ir_tl%OPDPS(lay,j) /  &
              & (transmission_scatt_ir%OPDPS(lay,j) ** 2)
            transmission_scatt_ir_tl%AZPHDO(lay,j) =                                               &
              & transmission_scatt_ir_tl%AZPHDOTOT(lay,j) / transmission_scatt_ir%OPDPS(lay,j) -  &
              & transmission_scatt_ir%AZPHDOTOT(lay,j) * transmission_scatt_ir_tl%OPDPS(lay,j) /  &
              & (transmission_scatt_ir%OPDPS(lay,j) ** 2)
          ENDIF
          OPDPCLDLSUN(lay,j) = &
            & transmission_scatt_ir_tl%OPDPCLDLA(lay,j) * raytracing%patheff(lay,prof) * coef%ff_gam(chan) +             &
            & raytracing_tl%patheff(lay,prof) * transmission_scatt_ir%OPDPCLDLA(lay,j) * coef%ff_gam(chan)
        ENDIF
      ENDDO !nlayers
    ENDDO !nchannels
  ENDIF !opts%rt_ir%addclouds
!-----Compute optical parameters for each stream--------------------------------
  IF (dosolar) transmission_scatt_ir_tl%SSA = 0._jprb
  DO j = 1, nchannels
    chan = chanprof(j)%chan
    prof = chanprof(j)%prof

    DO lay = 1, nlayers
      transmission_scatt_ir_tl%OPDPACL(0,lay,j) = OPDPAERL(lay,j)
      IF (opts%rt_ir%addclouds) THEN
        transmission_scatt_ir_tl%OPDPACL(1,lay,j) = OPDPAERL(lay,j) + &
                                                    OPDPCLDL(lay,j)
      ENDIF
    ENDDO

    IF (solar(j)) THEN
      DO lay = 1, nlayers

        transmission_scatt_ir_tl%OPDPACLSUN(0,lay,j) = OPDPAERLSUN(lay,j)
        IF (opts%rt_ir%addclouds) THEN
          transmission_scatt_ir_tl%OPDPACLSUN(1,lay,j) = OPDPAERLSUN(lay,j) + &
                                                         OPDPCLDLSUN(lay,j)
        ENDIF

        transmission_scatt_ir_tl%AZPHACUP(0,lay,j) = transmission_scatt_ir_tl%AZPHAERUPA(lay,j)
        transmission_scatt_ir_tl%AZPHACDO(0,lay,j) = transmission_scatt_ir_tl%AZPHAERDOA(lay,j)
        transmission_scatt_ir_tl%OPDPABS(0,lay,j) = &
          & transmission_scatt_ir_tl%OPDPAAER(lay,j) * raytracing%patheff(lay,prof) * coef%ff_gam(chan) + &
          & raytracing_tl%patheff(lay,prof) * transmission_scatt_ir%OPDPAAER(lay,j) * coef%ff_gam(chan)
        transmission_scatt_ir_tl%OPDPSCA(0,lay,j) = &
          & transmission_scatt_ir_tl%OPDPSAER(lay,j) * raytracing%patheff(lay,prof) * coef%ff_gam(chan) + &
          & raytracing_tl%patheff(lay,prof) * transmission_scatt_ir%OPDPSAER(lay,j) * coef%ff_gam(chan)

        IF (opts%rt_ir%addclouds) THEN

            IF (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j) /= 0._jprb) THEN
              transmission_scatt_ir_tl%AZPHACUP(1,lay,j) = &
                & transmission_scatt_ir_tl%AZPHAERUPA(lay,j) * transmission_scatt_ir%OPDPSAER(lay,j) / &
                & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) + &
                & transmission_scatt_ir_tl%OPDPSAER(lay,j) * transmission_scatt_ir%AZPHAERUPA(lay,j) / &
                & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) - &
                & transmission_scatt_ir_tl%OPDPS(lay,j) * &
                & transmission_scatt_ir%AZPHAERUPA(lay,j) * transmission_scatt_ir%OPDPSAER(lay,j) / &
                & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) ** 2 - &
                & transmission_scatt_ir_tl%OPDPSAER(lay,j) * transmission_scatt_ir%AZPHAERUPA(lay,j) * &
                & transmission_scatt_ir%OPDPSAER(lay,j) / &
                & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) ** 2 + &
                & transmission_scatt_ir_tl%AZPHUP(lay,j) * transmission_scatt_ir%OPDPS(lay,j) / &
                & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) + &
                & transmission_scatt_ir_tl%OPDPS(lay,j) * transmission_scatt_ir%AZPHUP(lay,j) / &
                & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) - &
                & transmission_scatt_ir_tl%OPDPS(lay,j) * &
                & transmission_scatt_ir%AZPHUP(lay,j) * transmission_scatt_ir%OPDPS(lay,j) / &
                & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) ** 2 - &
                & transmission_scatt_ir_tl%OPDPSAER(lay,j) * transmission_scatt_ir%AZPHUP(lay,j) * &
                & transmission_scatt_ir%OPDPS(lay,j) / &
                & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) ** 2
              transmission_scatt_ir_tl%AZPHACDO(1,lay,j) = &
                & transmission_scatt_ir_tl%AZPHAERDOA(lay,j) * transmission_scatt_ir%OPDPSAER(lay,j) / &
                & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) + &
                & transmission_scatt_ir_tl%OPDPSAER(lay,j) * transmission_scatt_ir%AZPHAERDOA(lay,j) / &
                & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) - &
                & transmission_scatt_ir_tl%OPDPS(lay,j) * &
                & transmission_scatt_ir%AZPHAERDOA(lay,j) * transmission_scatt_ir%OPDPSAER(lay,j) / &
                & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) ** 2 - &
                & transmission_scatt_ir_tl%OPDPSAER(lay,j) * transmission_scatt_ir%AZPHAERDOA(lay,j) * &
                & transmission_scatt_ir%OPDPSAER(lay,j) / &
                & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) ** 2 + &
                & transmission_scatt_ir_tl%AZPHDO(lay,j) * transmission_scatt_ir%OPDPS(lay,j) / &
                & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) + &
                & transmission_scatt_ir_tl%OPDPS(lay,j) * transmission_scatt_ir%AZPHDO(lay,j) / &
                & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) - &
                & transmission_scatt_ir_tl%OPDPS(lay,j) * &
                & transmission_scatt_ir%AZPHDO(lay,j) * transmission_scatt_ir%OPDPS(lay,j) / &
                & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) ** 2 - &
                & transmission_scatt_ir_tl%OPDPSAER(lay,j) * transmission_scatt_ir%AZPHDO(lay,j) * &
                & transmission_scatt_ir%OPDPS(lay,j) / &
                & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) ** 2
            ELSE
              transmission_scatt_ir_tl%AZPHACUP(1,lay,j) = 0._jprb
              transmission_scatt_ir_tl%AZPHACDO(1,lay,j) = 0._jprb
            ENDIF

          transmission_scatt_ir_tl%OPDPABS(1,lay,j) = &
            & raytracing_tl%patheff(lay,prof) * coef%ff_gam(chan) * &
            & (transmission_scatt_ir%OPDPAAER(lay,j) + transmission_scatt_ir%OPDPA(lay,j)) + &
            & (transmission_scatt_ir_tl%OPDPAAER(lay,j) + transmission_scatt_ir_tl%OPDPA(lay,j)) * &
            & raytracing%patheff(lay,prof) * coef%ff_gam(chan)
          transmission_scatt_ir_tl%OPDPSCA(1,lay,j) = &
            & raytracing_tl%patheff(lay,prof) * coef%ff_gam(chan) * &
            & (transmission_scatt_ir%OPDPSAER(lay,j) + transmission_scatt_ir%OPDPS(lay,j)) + &
            & (transmission_scatt_ir_tl%OPDPSAER(lay,j) + transmission_scatt_ir_tl%OPDPS(lay,j)) * &
            & raytracing%patheff(lay,prof) * coef%ff_gam(chan)
        ENDIF

      ENDDO ! layers
    ENDIF ! solar channel

    DO ist = 0, ircld%nstream(prof)
      IF (ist == 0) THEN
        IF (opts%rt_ir%addaerosl) THEN
          OPD = 0._jprb
          transmission_scatt_ir_stream_tl%OPDPAC(1,ist,j) = 0._jprb
          IF (solar(j)) THEN
            OPDPSUN = 0._jprb
            transmission_scatt_ir_stream_tl%OPDPACSUN(1,ist,j) = 0._jprb
          ENDIF
          DO lay = 1, nlayers
            lev = lay + 1
            IF (solar(j)) THEN
              OPDPSUN = OPDPSUN + transmission_scatt_ir_tl%OPDPACLSUN(0,lay,j)
              transmission_scatt_ir_stream_tl%OPDPACSUN(lev,ist,j) = OPDPSUN
            ENDIF
            OPD = OPD + transmission_scatt_ir_tl%OPDPACL(0,lay,j)
            transmission_scatt_ir_stream_tl%OPDPAC(lev,ist,j) = OPD
          ENDDO ! layers
        ELSE IF (opts%rt_ir%addclouds) THEN
          transmission_scatt_ir_stream_tl%OPDPAC(:,ist,j) = 0._jprb
          IF (solar(j)) THEN
            transmission_scatt_ir_stream_tl%OPDPACSUN(:,ist,j) = 0._jprb
          ENDIF
        ENDIF ! opts%rt_ir%addaerosl
      ELSE
        OPD = 0._jprb
        transmission_scatt_ir_stream_tl%OPDPAC(1,ist,j) = 0._jprb
        IF (solar(j)) THEN
          OPDPSUN = 0._jprb
          transmission_scatt_ir_stream_tl%OPDPACSUN(1,ist,j) = 0._jprb
        ENDIF
        DO lay = 1, nlayers
          lev = lay + 1
          isti = ircld%icldarr(ist,lay,prof)
          IF (solar(j)) THEN
            OPDPSUN = OPDPSUN + transmission_scatt_ir_tl%OPDPACLSUN(isti,lay,j)
            transmission_scatt_ir_stream_tl%OPDPACSUN(lev,ist,j) = OPDPSUN
          ENDIF
          OPD = OPD + transmission_scatt_ir_tl%OPDPACL(isti,lay,j)
          transmission_scatt_ir_stream_tl%OPDPAC(lev,ist,j) = OPD
        ENDDO ! layers
      ENDIF ! ist=0
    ENDDO ! istream
  ENDDO ! channels
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDPSCATTIR_TL', 1_jpim, ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE int_phase_fn_tl(musat, musat_tl, musun, musun_tl, relazi, zminphadiff, &
                             pha, cospha, ipha, phasint_tl, pha_tl)
    ! Interpolate phase function to scattering angle
    ! pha_tl may be omitted if the phase function is static
    REAL(KIND=jprb),           INTENT(IN)  :: musat, musat_tl, musun, musun_tl, relazi
    REAL(KIND=jprb),           INTENT(IN)  :: zminphadiff
    REAL(KIND=jprb),           INTENT(IN)  :: pha(:)
    REAL(KIND=jprb),           INTENT(IN)  :: cospha(:)
    INTEGER(KIND=jpim),        INTENT(IN)  :: ipha(:)
    REAL(KIND=jprb), OPTIONAL, INTENT(IN)  :: pha_tl(:)
    REAL(KIND=jprb),           INTENT(OUT) :: phasint_tl

    INTEGER(KIND=jpim) :: ikk, kk
    REAL(KIND=jprb)    :: ztmpx, ztmpx_tl
    REAL(KIND=jprb)    :: scattangle, scattangle_tl, deltap, deltap_tl, delta

    phasint_tl = 0._jprb

    ztmpx = SQRT((1._jprb - musat ** 2) * (1._jprb - musun ** 2))
    IF (ABS(musat) == 1._jprb .OR. ABS(musun) == 1._jprb) THEN
      ztmpx_tl = 0._jprb
    ELSE
      ztmpx_tl = - ((1._jprb - musun ** 2) * musat * musat_tl + &
                    (1._jprb - musat ** 2) * musun * musun_tl) / ztmpx
    ENDIF

    scattangle = musat * musun + ztmpx * COS(relazi * deg2rad)
    scattangle_tl = musat_tl * musun + musat * musun_tl + ztmpx_tl * COS(relazi * deg2rad)
    ikk        = MAX(1_jpim, INT(ACOS(scattangle) * zminphadiff, jpim))
    kk         = ipha(ikk) - 1_jpim
    deltap     = pha(kk + 1) - pha(kk)
    delta      = cospha(kk) - cospha(kk + 1)
    IF (PRESENT(pha_tl)) THEN
      deltap_tl   = pha_tl(kk + 1) - pha_tl(kk)
      phasint_tl  = pha_tl(kk) + deltap_tl * (cospha(kk) - scattangle) / delta - &
                                 deltap * scattangle_tl / delta
    ELSE
      phasint_tl  = - deltap * scattangle_tl / delta
    ENDIF

  END SUBROUTINE int_phase_fn_tl

END SUBROUTINE rttov_opdpscattir_tl
