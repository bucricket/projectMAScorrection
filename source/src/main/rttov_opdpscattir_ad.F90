!     Set up aerosols optical parameters for a climatological profile
SUBROUTINE rttov_opdpscattir_ad( &
            & nlayers,                         &
            & chanprof,                        &
            & opts,                            &
            & aux,                             &
            & aux_ad,                          &
            & profiles,                        &
            & profiles_ad,                     &
            & profiles_dry_ad,                 &
            & aer_opt_param,                   &
            & cld_opt_param,                   &
            & dosolar,                         &
            & solar,                           &
            & coef,                            &
            & coef_scatt_ir,                   &
            & raytracing,                      &
            & raytracing_ad,                   &
            & transmission_scatt_ir,           &
            & transmission_scatt_ir_ad,        &
            & transmission_scatt_ir_stream,    &
            & transmission_scatt_ir_stream_ad, &
            & optp,                            &
            & ircld,                           &
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
! Current Code Owner: SAF NWP
!
!     Method:
!
!     Owner:
!     EUMETSAT
!
!     History:
!     Version      Date        Comment
!     1           01/6/2004    Marco Matricardi. ECMWF.
!     1.1         05/02/2007   Removed polarisation index R Saunders
!     1.2         11/10/2007   iaernum & iaertyp in profile_aux
!     1.3         15/08/2009   User defined ToA. Layers distinct from levels (P.Rayer)
!     1.4         03/11/2009   Transmittances / optical depths on levels (A Geer)
!     1.5         02/12/2009   Fixed a number of bugs due to the wrong assumption that aerosol/cloud
!                              related quantities are defined on levels (thay are layer
!                              average quantities). Marco Matricardi
!     1.6         02/12/2009   Introduced multiple cloud types in a single layer. Pathsun, Pathsat and
!                              related quantities are now layer arrays (Marco Matricardi)
!     1.7         05/07/2010   Remove addsolar flag from profiles structure (J Hocking)
!     1.8         02/08/2010   Allow variable pressure levels (lgradp) (N Bormann)
!     1.9         14/12/2010   Use traj0_sta%solar array to flag channels for which solar calculations
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
  TYPE(rttov_chanprof            ), INTENT(IN)    :: chanprof   (:)
  TYPE(rttov_options )            , INTENT(IN)    :: opts
  TYPE(profile_type              ), INTENT(IN)    :: profiles   (:)
  TYPE(profile_type              ), INTENT(INOUT) :: profiles_ad(SIZE(profiles))
  TYPE(profile_type              ), INTENT(INOUT) :: profiles_dry_ad(SIZE(profiles))
  TYPE(rttov_opt_param), OPTIONAL , INTENT(IN)    :: aer_opt_param
  TYPE(rttov_opt_param), OPTIONAL , INTENT(IN)    :: cld_opt_param
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: dosolar
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: solar(SIZE(chanprof))
  TYPE(profile_aux               ), INTENT(IN)    :: aux
  TYPE(profile_aux               ), INTENT(INOUT) :: aux_ad
  TYPE(rttov_coef                ), INTENT(IN)    :: coef
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_ad
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir_stream
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream_ad
  TYPE(rttov_coef_scatt_ir       ), INTENT(IN)    :: coef_scatt_ir
  TYPE(rttov_optpar_ir           ), INTENT(IN)    :: optp
  TYPE(ircld_type                ), INTENT(IN)    :: ircld
  TYPE(ircld_type                ), INTENT(INOUT) :: ircld_ad
  TYPE(raytracing_type           ), INTENT(IN)    :: raytracing
  TYPE(raytracing_type           ), INTENT(INOUT) :: raytracing_ad
!INTF_END
!     End of subroutine arguments
#include "rttov_baran2013_calc_optpar.interface"
#include "rttov_baran2014_calc_optpar.interface"
#include "rttov_baran_calc_phase.interface"
#include "rttov_baran2013_calc_optpar_ad.interface"
#include "rttov_baran2014_calc_optpar_ad.interface"
#include "rttov_baran_calc_phase_ad.interface"

!       Local scalars:
  INTEGER(KIND=jpim) :: j, chan, chan1, i, prof, ish, k, kk, ityp, iae
  INTEGER(KIND=jpim) :: lev, lay, lctyp
  INTEGER(KIND=jpim) :: ist, isti
  REAL   (KIND=jprb) :: OPD, OPDSUN
  REAL   (KIND=jprb) :: absch_ad
  REAL   (KIND=jprb) :: scach_ad
  REAL   (KIND=jprb) :: bparh_ad
  REAL   (KIND=jprb) :: afac
  REAL   (KIND=jprb) :: sfac
  REAL   (KIND=jprb) :: gfac
  REAL   (KIND=jprb) :: deltaice
  REAL   (KIND=jprb) :: delth
  REAL   (KIND=jprb) :: frach_ad
  REAL   (KIND=jprb) :: phasint_ad
  REAL   (KIND=jprb) :: musat_ad
  REAL   (KIND=jprb) :: musun_ad
  REAL   (KIND=jprb) :: musat
  REAL   (KIND=jprb) :: musun
  REAL   (KIND=jprb) :: frach
  REAL   (KIND=jprb) :: absch
  REAL   (KIND=jprb) :: scach
  REAL   (KIND=jprb) :: bparh
  REAL   (KIND=jprb) :: zminphadiff, relazi
!       Local arrays:
  REAL   (KIND=jprb) :: OPDPAERL   (nlayers,SIZE(chanprof))
  REAL   (KIND=jprb) :: OPDPCLDL   (nlayers,SIZE(chanprof))
  REAL   (KIND=jprb) :: OPDPCLDLSUN(nlayers,SIZE(chanprof))
  REAL   (KIND=jprb) :: OPDPAERLSUN(nlayers,SIZE(chanprof))
  REAL   (KIND=jprb) :: pfac       (coef_scatt_ir%fmv_aer_ph)
  REAL   (KIND=jprb) :: phash      (coef_scatt_ir%fmv_aer_ph)
  REAL   (KIND=jprb) :: phash_ad   (coef_scatt_ir%fmv_aer_ph)
  REAL   (KIND=jprb) :: phasice    (coef_scatt_ir%fmv_icl_ph)
  REAL   (KIND=jprb) :: phasice_ad (coef_scatt_ir%fmv_icl_ph)
  REAL   (KIND=jprb) :: deltap     (coef_scatt_ir%fmv_icl_ph)
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: nchannels
  
  REAL(KIND=jprb)  :: abso, abso_ad
  REAL(KIND=jprb)  :: sca,  sca_ad
  REAL(KIND=jprb)  :: bpr,  bpr_ad
  REAL(KIND=jprb)  :: asym, asym_ad

  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!-----End of header-----------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDPSCATTIR_AD', 0_jpim, ZHOOK_HANDLE)
  nprofiles        = SIZE(profiles)
  nchannels        = SIZE(chanprof)
  frach_ad         = 0._jprb
  absch_ad         = 0._jprb
  scach_ad         = 0._jprb
  bparh_ad         = 0._jprb
  phasint_ad       = 0._jprb
  musat_ad         = 0._jprb
  musun_ad         = 0._jprb
  phash_ad         = 0._jprb
  phasice_ad       = 0._jprb
  abso_ad          = 0._jprb
  sca_ad           = 0._jprb
  bpr_ad           = 0._jprb
  asym_ad          = 0._jprb
  OPD              = 0._jprb
  OPDPAERL         = 0._jprb
  OPDPCLDL         = 0._jprb
  IF (dosolar) THEN
    OPDSUN      = 0._jprb
    OPDPAERLSUN = 0._jprb
    OPDPCLDLSUN = 0._jprb
  ENDIF

  ! transmission_scatt_ir_ad and transmission_scatt_ir_stream_ad are initialised in rttov_ad

!-----Compute optical parameters for each stream------------------------------------------
  DO j = nchannels, 1,  - 1
    chan = chanprof(j)%chan
    prof = chanprof(j)%prof
    DO ist = ircld%nstream(prof), 0,  - 1
      IF (ist == 0) THEN
        IF (opts%rt_ir%addaerosl) THEN
          DO lay = nlayers, 1,  - 1
            lev = lay + 1
            OPD = OPD + transmission_scatt_ir_stream_ad%OPDPAC(lev,ist,j)
            transmission_scatt_ir_ad%OPDPACL(0,lay,j) = &
              & transmission_scatt_ir_ad%OPDPACL(0,lay,j) + OPD
            IF (solar(j)) THEN
              OPDSUN = OPDSUN + transmission_scatt_ir_stream_ad%OPDPACSUN(lev,ist,j)
              transmission_scatt_ir_ad%OPDPACLSUN(0,lay,j) = &
                & transmission_scatt_ir_ad%OPDPACLSUN(0,lay,j) + OPDSUN
            ENDIF
          ENDDO
          OPD = 0._jprb
          transmission_scatt_ir_stream_ad%OPDPAC(1,ist,j) = 0._jprb
          IF (solar(j)) THEN
            OPDSUN = 0._jprb
            transmission_scatt_ir_stream_ad%OPDPACSUN(1,ist,j) = 0._jprb
          ENDIF
        ELSE IF (opts%rt_ir%addclouds) THEN
          IF (solar(j)) THEN
            transmission_scatt_ir_stream_ad%OPDPACSUN(:,ist,j) = 0._jprb
          ENDIF
          transmission_scatt_ir_stream_ad%OPDPAC(:,ist,j) = 0._jprb
        ENDIF
      ELSE ! ist/=0
        DO lay = nlayers, 1,  - 1
          lev = lay + 1
          isti = ircld%icldarr(ist,lay,prof)
          OPD = OPD + transmission_scatt_ir_stream_ad%OPDPAC(lev,ist,j)
          transmission_scatt_ir_ad%OPDPACL(isti,lay,j) = &
            & transmission_scatt_ir_ad%OPDPACL(isti,lay,j) + OPD
          IF (solar(j)) THEN
            OPDSUN = OPDSUN + transmission_scatt_ir_stream_ad%OPDPACSUN(lev,ist,j)
            transmission_scatt_ir_ad%OPDPACLSUN(isti,lay,j) = &
              & transmission_scatt_ir_ad%OPDPACLSUN(isti,lay,j) + OPDSUN
          ENDIF
        ENDDO
        OPD = 0._jprb
        transmission_scatt_ir_stream_ad%OPDPAC(1,ist,j) = 0._jprb
        IF (solar(j)) THEN
          OPDSUN = 0._jprb
          transmission_scatt_ir_stream_ad%OPDPACSUN(1,ist,j) = 0._jprb
        ENDIF
      ENDIF ! ist == 0
    ENDDO ! ist

    IF (solar(j)) THEN
      DO lay = nlayers, 1, -1

        IF (opts%rt_ir%addclouds) THEN
          raytracing_ad%patheff(lay,prof)             = raytracing_ad%patheff(lay,prof) + &
            & transmission_scatt_ir_ad%OPDPSCA(1,lay,j) * coef%ff_gam(chan) * &
            & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j))
          transmission_scatt_ir_ad%OPDPS(lay,j)      = transmission_scatt_ir_ad%OPDPS(lay,j) + &
            & transmission_scatt_ir_ad%OPDPSCA(1,lay,j) * &
            & raytracing%patheff(lay,prof) * coef%ff_gam(chan)
          transmission_scatt_ir_ad%OPDPSAER(lay,j)   = transmission_scatt_ir_ad%OPDPSAER(lay,j) + &
            & transmission_scatt_ir_ad%OPDPSCA(1,lay,j) * &
            & raytracing%patheff(lay,prof) * coef%ff_gam(chan)
          raytracing_ad%patheff(lay,prof)             = raytracing_ad%patheff(lay,prof) + &
            & transmission_scatt_ir_ad%OPDPABS(1,lay,j) * coef%ff_gam(chan) * &
            & (transmission_scatt_ir%OPDPA(lay,j) + transmission_scatt_ir%OPDPAAER(lay,j))
          transmission_scatt_ir_ad%OPDPA(lay,j)      = transmission_scatt_ir_ad%OPDPA(lay,j) + &
            & transmission_scatt_ir_ad%OPDPABS(1,lay,j) * &
            & raytracing%patheff(lay,prof) * coef%ff_gam(chan)
          transmission_scatt_ir_ad%OPDPAAER(lay,j)   = transmission_scatt_ir_ad%OPDPAAER(lay,j) + &
            & transmission_scatt_ir_ad%OPDPABS(1,lay,j) * &
            & raytracing%patheff(lay,prof) * coef%ff_gam(chan)

          IF ((transmission_scatt_ir%OPDPS(lay,j) +      &
            & transmission_scatt_ir%OPDPSAER(lay,j)) /= 0._jprb) THEN
            transmission_scatt_ir_ad%OPDPSAER(lay,j)   = transmission_scatt_ir_ad%OPDPSAER(lay,j) - &
              & transmission_scatt_ir_ad%AZPHACDO(1,lay,j) * transmission_scatt_ir%AZPHDO(lay,j) * &
              & transmission_scatt_ir%OPDPS(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) ** 2
            transmission_scatt_ir_ad%OPDPS(lay,j)      = transmission_scatt_ir_ad%OPDPS(lay,j) - &
              & transmission_scatt_ir_ad%AZPHACDO(1,lay,j) * &
              & transmission_scatt_ir%AZPHDO(lay,j) * transmission_scatt_ir%OPDPS(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) ** 2
            transmission_scatt_ir_ad%OPDPS(lay,j)      = transmission_scatt_ir_ad%OPDPS(lay,j) + &
              & transmission_scatt_ir_ad%AZPHACDO(1,lay,j) * &
              & transmission_scatt_ir%AZPHDO(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j))
            transmission_scatt_ir_ad%AZPHDO(lay,j)     = transmission_scatt_ir_ad%AZPHDO(lay,j) + &
              & transmission_scatt_ir_ad%AZPHACDO(1,lay,j) * transmission_scatt_ir%OPDPS(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j))
            transmission_scatt_ir_ad%OPDPSAER(lay,j)   = transmission_scatt_ir_ad%OPDPSAER(lay,j) - &
              & transmission_scatt_ir_ad%AZPHACDO(1,lay,j) * &
              & transmission_scatt_ir%AZPHAERDOA(lay,j) * transmission_scatt_ir%OPDPSAER(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) ** 2
            transmission_scatt_ir_ad%OPDPS(lay,j)      = transmission_scatt_ir_ad%OPDPS(lay,j) - &
              & transmission_scatt_ir_ad%AZPHACDO(1,lay,j) * &
              & transmission_scatt_ir%AZPHAERDOA(lay,j) * transmission_scatt_ir%OPDPSAER(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) ** 2
            transmission_scatt_ir_ad%OPDPSAER(lay,j)   = transmission_scatt_ir_ad%OPDPSAER(lay,j) + &
              & transmission_scatt_ir_ad%AZPHACDO(1,lay,j) * &
              & transmission_scatt_ir%AZPHAERDOA(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j))
            transmission_scatt_ir_ad%AZPHAERDOA(lay,j) = transmission_scatt_ir_ad%AZPHAERDOA(lay,j) + &
              & transmission_scatt_ir_ad%AZPHACDO(1,lay,j) * &
              & transmission_scatt_ir%OPDPSAER(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j))
            transmission_scatt_ir_ad%OPDPSAER(lay,j)   = transmission_scatt_ir_ad%OPDPSAER(lay,j) - &
              & transmission_scatt_ir_ad%AZPHACUP(1,lay,j) * transmission_scatt_ir%AZPHUP(lay,j) * &
              & transmission_scatt_ir%OPDPS(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) ** 2
            transmission_scatt_ir_ad%OPDPS(lay,j)      = transmission_scatt_ir_ad%OPDPS(lay,j) - &
              & transmission_scatt_ir_ad%AZPHACUP(1,lay,j) * &
              & transmission_scatt_ir%AZPHUP(lay,j) * transmission_scatt_ir%OPDPS(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) ** 2
            transmission_scatt_ir_ad%OPDPS(lay,j)      = transmission_scatt_ir_ad%OPDPS(lay,j) + &
              & transmission_scatt_ir_ad%AZPHACUP(1,lay,j) * &
              & transmission_scatt_ir%AZPHUP(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j))
            transmission_scatt_ir_ad%AZPHUP(lay,j)     = transmission_scatt_ir_ad%AZPHUP(lay,j) + &
              & transmission_scatt_ir_ad%AZPHACUP(1,lay,j) * transmission_scatt_ir%OPDPS(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j))
            transmission_scatt_ir_ad%OPDPSAER(lay,j)   = transmission_scatt_ir_ad%OPDPSAER(lay,j) - &
              & transmission_scatt_ir_ad%AZPHACUP(1,lay,j) * &
              & transmission_scatt_ir%AZPHAERUPA(lay,j) * transmission_scatt_ir%OPDPSAER(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) ** 2
            transmission_scatt_ir_ad%OPDPS(lay,j)      = transmission_scatt_ir_ad%OPDPS(lay,j) - &
              & transmission_scatt_ir_ad%AZPHACUP(1,lay,j) * &
              & transmission_scatt_ir%AZPHAERUPA(lay,j) * transmission_scatt_ir%OPDPSAER(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) ** 2
            transmission_scatt_ir_ad%OPDPSAER(lay,j)   = transmission_scatt_ir_ad%OPDPSAER(lay,j) + &
              & transmission_scatt_ir_ad%AZPHACUP(1,lay,j) * &
              & transmission_scatt_ir%AZPHAERUPA(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j))
            transmission_scatt_ir_ad%AZPHAERUPA(lay,j) = transmission_scatt_ir_ad%AZPHAERUPA(lay,j) + &
              & transmission_scatt_ir_ad%AZPHACUP(1,lay,j) * &
              & transmission_scatt_ir%OPDPSAER(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j))
          ELSE
            transmission_scatt_ir_ad%AZPHACUP(1,lay,j) = 0._jprb
            transmission_scatt_ir_ad%AZPHACDO(1,lay,j) = 0._jprb
          ENDIF
        ENDIF

        transmission_scatt_ir_ad%OPDPSAER(lay,j) = transmission_scatt_ir_ad%OPDPSAER(lay,j) + &
          & transmission_scatt_ir_ad%OPDPSCA(0,lay,j) * &
          & raytracing%patheff(lay,prof) * coef%ff_gam(chan)
        raytracing_ad%patheff(lay,prof) = raytracing_ad%patheff(lay,prof) + &
          & transmission_scatt_ir_ad%OPDPSCA(0,lay,j) * transmission_scatt_ir%OPDPSAER(lay,j) * &
          & coef%ff_gam(chan)
        transmission_scatt_ir_ad%OPDPAAER(lay,j) = transmission_scatt_ir_ad%OPDPAAER(lay,j) + &
          & transmission_scatt_ir_ad%OPDPABS(0,lay,j) * &
          & raytracing%patheff(lay,prof) * coef%ff_gam(chan)
        raytracing_ad%patheff(lay,prof) = raytracing_ad%patheff(lay,prof) + &
          & transmission_scatt_ir_ad%OPDPABS(0,lay,j) * transmission_scatt_ir%OPDPAAER(lay,j) * &
          & coef%ff_gam(chan)
        transmission_scatt_ir_ad%AZPHAERDOA(lay,j) = &
          & transmission_scatt_ir_ad%AZPHAERDOA(lay,j) + transmission_scatt_ir_ad%AZPHACDO(0,lay,j)
        transmission_scatt_ir_ad%AZPHAERUPA(lay,j) = &
          & transmission_scatt_ir_ad%AZPHAERUPA(lay,j) + transmission_scatt_ir_ad%AZPHACUP(0,lay,j)

        IF (opts%rt_ir%addclouds) THEN
          OPDPCLDLSUN(lay,j) = OPDPCLDLSUN(lay,j) + transmission_scatt_ir_ad%OPDPACLSUN(1,lay,j)
          OPDPAERLSUN(lay,j) = OPDPAERLSUN(lay,j) + transmission_scatt_ir_ad%OPDPACLSUN(1,lay,j)
        ENDIF
        OPDPAERLSUN(lay,j) = OPDPAERLSUN(lay,j) + transmission_scatt_ir_ad%OPDPACLSUN(0,lay,j)

      ENDDO ! layers
    ENDIF ! solar channel

    DO lay = nlayers, 1, -1
      IF (opts%rt_ir%addclouds) THEN
        OPDPCLDL(lay,j) = OPDPCLDL(lay,j) + transmission_scatt_ir_ad%OPDPACL(1,lay,j)
        OPDPAERL(lay,j) = OPDPAERL(lay,j) + transmission_scatt_ir_ad%OPDPACL(1,lay,j)
      ENDIF
      OPDPAERL(lay,j) = OPDPAERL(lay,j) + transmission_scatt_ir_ad%OPDPACL(0,lay,j)
    ENDDO

  ENDDO
!-------------------------------------------------------------------------------
!         2.   CALCULATE OPTICAL DEPTHS OF CLOUDS
!-------------------------------------------------------------------------------
  IF (opts%rt_ir%addclouds) THEN
    DO j = nchannels, 1,  - 1
      chan = chanprof(j)%chan
      prof = chanprof(j)%prof
      ish = profiles(prof)%ish
      relazi  = profiles(prof)%azangle - profiles(prof)%sunazangle
      DO lay = nlayers, 1,  - 1
        IF (solar(j)) THEN
!---------------Compute cloud  optical parameters ----------------------------------------
          transmission_scatt_ir_ad%OPDPCLDLA(lay,j) = transmission_scatt_ir_ad%OPDPCLDLA(lay,j) +    &
            & OPDPCLDLSUN(lay,j) * raytracing%patheff(lay,prof) * coef%ff_gam(chan)
          raytracing_ad%patheff(lay,prof)            = raytracing_ad%patheff(lay,prof) +      &
            & OPDPCLDLSUN(lay,j) * transmission_scatt_ir%OPDPCLDLA(lay,j) * coef%ff_gam(chan)
          IF (transmission_scatt_ir%OPDPS(lay,j) /= 0._jprb) THEN
            transmission_scatt_ir_ad%AZPHUPTOT(lay,j) = transmission_scatt_ir_ad%AZPHUPTOT(lay,j) +  &
              & transmission_scatt_ir_ad%AZPHUP(lay,j) / transmission_scatt_ir%OPDPS(lay,j)
            transmission_scatt_ir_ad%OPDPS(lay,j)     = transmission_scatt_ir_ad%OPDPS(lay,j) -      &
              & transmission_scatt_ir_ad%AZPHUP(lay,j) * transmission_scatt_ir%AZPHUPTOT(lay,j) /    &
              & (transmission_scatt_ir%OPDPS(lay,j) ** 2)
            transmission_scatt_ir_ad%AZPHDOTOT(lay,j) = transmission_scatt_ir_ad%AZPHDOTOT(lay,j) +  &
              & transmission_scatt_ir_ad%AZPHDO(lay,j) / transmission_scatt_ir%OPDPS(lay,j)
            transmission_scatt_ir_ad%OPDPS(lay,j)     = transmission_scatt_ir_ad%OPDPS(lay,j) -      &
              & transmission_scatt_ir_ad%AZPHDO(lay,j) * transmission_scatt_ir%AZPHDOTOT(lay,j) /    &
              & (transmission_scatt_ir%OPDPS(lay,j) ** 2)
          ENDIF
        ENDIF
        transmission_scatt_ir_ad%OPDPCLDLA(lay,j) = transmission_scatt_ir_ad%OPDPCLDLA(lay,j) + &
          & OPDPCLDL(lay,j) * raytracing%pathsat(lay,prof) * coef%ff_gam(chan)
        raytracing_ad%pathsat(lay,prof)            = raytracing_ad%pathsat(lay,prof) +            &
          & OPDPCLDL(lay,j) * transmission_scatt_ir%OPDPCLDLA(lay,j) * coef%ff_gam(chan)
        transmission_scatt_ir_ad%OPDPA(lay,j)     =      &
          & transmission_scatt_ir_ad%OPDPA(lay,j) + transmission_scatt_ir_ad%OPDPCLDLA(lay,j)
        transmission_scatt_ir_ad%OPDPS(lay,j)     = transmission_scatt_ir_ad%OPDPS(lay,j) +     &
          & transmission_scatt_ir_ad%OPDPCLDLA(lay,j) * transmission_scatt_ir%GPAR(lay,j)
        transmission_scatt_ir_ad%GPAR(lay,j)      = transmission_scatt_ir_ad%GPAR(lay,j) +      &
          & transmission_scatt_ir_ad%OPDPCLDLA(lay,j) * transmission_scatt_ir%OPDPS(lay,j)
        IF (transmission_scatt_ir%OPDPS(lay,j) /= 0._jprb) THEN
          transmission_scatt_ir_ad%GPARTOT(lay,j) = transmission_scatt_ir_ad%GPARTOT(lay,j) +   &
            & transmission_scatt_ir_ad%GPAR(lay,j) / transmission_scatt_ir%OPDPS(lay,j)
          transmission_scatt_ir_ad%OPDPS(lay,j)   = transmission_scatt_ir_ad%OPDPS(lay,j) -     &
            & transmission_scatt_ir_ad%GPAR(lay,j) * transmission_scatt_ir%GPARTOT(lay,j) /     &
            & (transmission_scatt_ir%OPDPS(lay,j) ** 2)
        ENDIF
!              transmission_scatt_ir_ad%GPAR (lay,j)    = 0._jprb

        IF (opts%rt_ir%user_cld_opt_param) THEN
          ! Input optical parameters are static

          IF (solar(j)) THEN
!-------------Average phase function for the downward scattered solar beam------------------
            musat       =  - 1._jprb / raytracing%pathsat(lay,prof)
            musun       =  - 1._jprb / raytracing%pathsun(lay,prof)
            zminphadiff = 1._jprb / (cld_opt_param%minphadiff * deg2rad)

            phasint_ad = phasint_ad + transmission_scatt_ir_ad%AZPHDOTOT(lay,j) * &
                       & cld_opt_param%sca(j,lay) * raytracing%ltick(lay,prof)
            ! AZPHDOCLS(:,:,1) is the direct PHASINT
            raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                       & transmission_scatt_ir_ad%AZPHDOTOT(lay,j) * cld_opt_param%sca(j,lay) *        &
                       & transmission_scatt_ir%AZPHDOCLS(1,lay,j)

            CALL int_phase_fn_ad(musat, musat_ad, musun, musun_ad, relazi, zminphadiff, &
                                 cld_opt_param%pha(j,lay,:), cld_opt_param%cosphangle,  &
                                 cld_opt_param%iphangle, phasint_ad)

            musat_ad = -1._jprb * musat_ad
            phasint_ad = 0._jprb
!-------------Average phase function for the upward scattered solar beam------
            musat    = -1._jprb * musat

            phasint_ad = phasint_ad + transmission_scatt_ir_ad%AZPHUPTOT(lay,j) * &
                       & cld_opt_param%sca(j,lay) * raytracing%ltick(lay,prof)
            ! AZPHUPCLS(:,:,1) is the direct PHASINT
            raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                       & transmission_scatt_ir_ad%AZPHUPTOT(lay,j) * cld_opt_param%sca(j,lay) *        &
                       & transmission_scatt_ir%AZPHUPCLS(1,lay,j)

            CALL int_phase_fn_ad(musat, musat_ad, musun, musun_ad, 180.0_jprb - relazi, zminphadiff, &
                                 cld_opt_param%pha(j,lay,:), cld_opt_param%cosphangle, &
                                 cld_opt_param%iphangle, phasint_ad)

            raytracing_ad%pathsat(lay,prof) = raytracing_ad%pathsat(lay,prof) - &
                                            & musat_ad / raytracing%pathsat(lay,prof) ** 2
            raytracing_ad%pathsun(lay,prof) = raytracing_ad%pathsun(lay,prof) + &
                                            & musun_ad / raytracing%pathsun(lay,prof) ** 2
            musat_ad = 0._jprb
            musun_ad = 0._jprb
            phasint_ad = 0._jprb
          ENDIF

          raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                                                      & transmission_scatt_ir_ad%GPARTOT(lay,j) *    &
                                                      & cld_opt_param%sca(j,lay) * cld_opt_param%bpr(j,lay)
          raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                                                      & transmission_scatt_ir_ad%OPDPS(lay,j) *      &
                                                      & cld_opt_param%sca(j,lay)
          raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                                                      & transmission_scatt_ir_ad%OPDPA(lay,j) *      &
                                                      & cld_opt_param%abs(j,lay)
        ELSE
!-----------------------------------------------------------------------------------------
!                 For water clouds use stored optical parameters
!-----------------------------------------------------------------------------------------
          DO lctyp = ncldtyp, 1, -1
            IF (ircld%cldtyp(lctyp, lay,prof) /= 0) THEN
              ityp  = ircld%cldtyp(lctyp, lay,prof)
              IF (ityp == 6_jpim .and. ish .GE. 3_jpim ) THEN
                IF (ish == 3_JPIM) THEN 
                   CALL rttov_baran2013_calc_optpar (optp, chan, &
                       & profiles(prof)%t(lay), profiles(prof)%cloud(ityp,lay), &
                       & abso, sca, bpr, asym)
                ELSEIF (ish == 4_JPIM) THEN 
                   CALL rttov_baran2014_calc_optpar (optp, chan, &
                       & profiles(prof)%t(lay), profiles(prof)%cloud(ityp,lay), &
                       & abso, sca, bpr, asym)
                ENDIF
              ENDIF
!-------------------Average phase function for the downward scattered solar beam------------
              IF (solar(j)) THEN
                IF (ityp == 6_jpim) THEN
!-----------------If ice clouds are present the phase function for for the current value
!                 of the effective generalized diameter is obtained by linear
!                 interpolation and then the azimuthally averaged value is computed
!-----------------------------------------------------------------------------------------
                  chan1 = coef_scatt_ir%icl_pha_index(chan)
                  IF( ish .LE. 2_jpim ) THEN
                    DO k = 1, coef_scatt_ir%fmv_icl_comp - 1
                      IF (aux%dg(lay,prof) >= coef_scatt_ir%fmv_icl_dg(k, ish) .AND.      &
                        & aux%dg(lay,prof) <= coef_scatt_ir%fmv_icl_dg(k + 1, ish)) THEN
                        deltap(:)    = (optp%optpicl(ish)%pha(chan1, k + 1, :) - optp%optpicl(ish)%pha(chan1, k, :))
                        deltaice     = (coef_scatt_ir%fmv_icl_dg(k + 1, ish) - coef_scatt_ir%fmv_icl_dg(k, ish))
                        phasice(:) = optp%optpicl(ish)%pha(chan1, k, :) +      &
                          & deltap(:) * (aux%dg(lay,prof) - coef_scatt_ir%fmv_icl_dg(k, ish)) / deltaice
                        EXIT
                      ENDIF
                    ENDDO
                  ELSE
                    CALL rttov_baran_calc_phase(asym, phangle, phasice)
                  ENDIF
!-----------------Average phase function for the downward scattered solar beam------------
                  musat       =  - 1._jprb / raytracing%pathsat(lay,prof)
                  musun       =  - 1._jprb / raytracing%pathsun(lay,prof)
                  zminphadiff = 1._jprb / (coef_scatt_ir%fmv_icl_ph_val_min * deg2rad)

                  transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) =               &
                              & transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) + &
                              & transmission_scatt_ir%AZPHDOCLS(ityp,lay,j) *   &
                              & transmission_scatt_ir_ad%AZPHDOTOT(lay,j)
                  transmission_scatt_ir_ad%AZPHDOCLS(ityp,lay,j) =               &
                              & transmission_scatt_ir_ad%AZPHDOCLS(ityp,lay,j) + &
                              & transmission_scatt_ir%OPDPSCLS(ityp,lay,j) *     &
                              & transmission_scatt_ir_ad%AZPHDOTOT(lay,j)
                  phasint_ad = phasint_ad + transmission_scatt_ir_ad%AZPHDOCLS(ityp,lay,j)

                  CALL int_phase_fn_ad(musat, musat_ad, musun, musun_ad, relazi, zminphadiff, &
                                       phasice, coef_scatt_ir%fmv_icl_ph_val_cos,             &
                                       coef_scatt_ir%ifmv_icl_ph_val, phasint_ad, phasice_ad)

                  musat_ad = -1._jprb * musat_ad
                  phasint_ad = 0._jprb
!-----------------Average phase function for the upward scattered solar beam--------------
                  musat    = -1._jprb * musat

                  transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) =               &
                              & transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) + &
                              & transmission_scatt_ir%AZPHUPCLS(ityp,lay,j) *   &
                              & transmission_scatt_ir_ad%AZPHUPTOT(lay,j)
                  transmission_scatt_ir_ad%AZPHUPCLS(ityp,lay,j) =               &
                              & transmission_scatt_ir_ad%AZPHUPCLS(ityp,lay,j) + &
                              & transmission_scatt_ir%OPDPSCLS(ityp,lay,j) *     &
                              & transmission_scatt_ir_ad%AZPHUPTOT(lay,j)
                  phasint_ad = phasint_ad + transmission_scatt_ir_ad%AZPHUPCLS(ityp,lay,j)

                  CALL int_phase_fn_ad(musat, musat_ad, musun, musun_ad, 180.0_jprb - relazi, zminphadiff, &
                                       phasice, coef_scatt_ir%fmv_icl_ph_val_cos, &
                                       coef_scatt_ir%ifmv_icl_ph_val, phasint_ad, phasice_ad)

                  raytracing_ad%pathsat(lay,prof) = raytracing_ad%pathsat(lay,prof) - &
                                                  & musat_ad / raytracing%pathsat(lay,prof) ** 2
                  raytracing_ad%pathsun(lay,prof) = raytracing_ad%pathsun(lay,prof) + &
                                                  & musun_ad / raytracing%pathsun(lay,prof) ** 2

                  musat_ad = 0._jprb
                  musun_ad = 0._jprb
                  phasint_ad = 0._jprb
!-----------------If ice clouds are present the phase function for for the current value
!                 of the effective generalized diameter is obtained by linear
!                 interpolation and then the azimuthally averaged value is computed
!----------------------------------------------------------------------------------------
                  IF( ISH .LE. 2_jpim ) THEN
                    DO k = 1, coef_scatt_ir%fmv_icl_comp - 1
                      IF (aux%dg(lay,prof) >= coef_scatt_ir%fmv_icl_dg(k, ish) .AND.      &
                        & aux%dg(lay,prof) <= coef_scatt_ir%fmv_icl_dg(k + 1, ish)) THEN
                        deltap(:) = (optp%optpicl(ish)%pha(chan1, k + 1, :) - optp%optpicl(ish)%pha(chan1, k, :))
                        deltaice  = (coef_scatt_ir%fmv_icl_dg(k + 1, ish) - coef_scatt_ir%fmv_icl_dg(k, ish))
                        DO kk = 1, coef_scatt_ir%fmv_icl_ph
                          aux_ad%dg(lay,prof) = aux_ad%dg(lay,prof) + phasice_ad(kk) * deltap(kk) / deltaice
                          phasice_ad(kk)      = 0._jprb
                        ENDDO
                        EXIT
                      ENDIF
                    ENDDO
                  ELSE
                    asym_ad = 0.0_jprb
                    CALL rttov_baran_calc_phase_ad(asym, asym_ad, phangle, phasice_ad)
                  ENDIF
                ELSE
!-----------------------------------------------------------------------------------------
!                         Water clouds
!-----------------------------------------------------------------------------------------
                  chan1 = coef_scatt_ir%wcl_pha_index(chan)
!-----------------------Average phase function for the downward scattered solar beam------------
                  musat       =  - 1._jprb / raytracing%pathsat(lay,prof)
                  musun       =  - 1._jprb / raytracing%pathsun(lay,prof)
                  zminphadiff = 1._jprb / (coef_scatt_ir%fmv_wcl_ph_val_min * deg2rad)

                  transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) =               &
                              & transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) + &
                              & transmission_scatt_ir%AZPHDOCLS(ityp,lay,j) *   &
                              & transmission_scatt_ir_ad%AZPHDOTOT(lay,j)
                  transmission_scatt_ir_ad%AZPHDOCLS(ityp,lay,j) =               &
                              & transmission_scatt_ir_ad%AZPHDOCLS(ityp,lay,j) + &
                              & transmission_scatt_ir%OPDPSCLS(ityp,lay,j) *     &
                              & transmission_scatt_ir_ad%AZPHDOTOT(lay,j)
                  phasint_ad = phasint_ad + transmission_scatt_ir_ad%AZPHDOCLS(ityp,lay,j)

                  CALL int_phase_fn_ad(musat, musat_ad, musun, musun_ad, relazi, zminphadiff,              &
                                       optp%optpwcl(ityp)%pha(chan1,1,:), coef_scatt_ir%fmv_wcl_ph_val_cos, &
                                       coef_scatt_ir%ifmv_wcl_ph_val, phasint_ad)

                  musat_ad = -1._jprb * musat_ad
                  phasint_ad = 0._jprb
!-------------------------Average phase function for the upward scattered solar beam--------------
                  musat    = -1._jprb * musat

                  transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) =               &
                              & transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) + &
                              & transmission_scatt_ir%AZPHUPCLS(ityp,lay,j) *   &
                              & transmission_scatt_ir_ad%AZPHUPTOT(lay,j)
                  transmission_scatt_ir_ad%AZPHUPCLS(ityp,lay,j) =               &
                              & transmission_scatt_ir_ad%AZPHUPCLS(ityp,lay,j) + &
                              & transmission_scatt_ir%OPDPSCLS(ityp,lay,j) *     &
                              & transmission_scatt_ir_ad%AZPHUPTOT(lay,j)
                  phasint_ad = phasint_ad + transmission_scatt_ir_ad%AZPHUPCLS(ityp,lay,j)

                  CALL int_phase_fn_ad(musat, musat_ad, musun, musun_ad, 180.0_jprb - relazi, zminphadiff, &
                                       optp%optpwcl(ityp)%pha(chan1,1,:), coef_scatt_ir%fmv_wcl_ph_val_cos, &
                                       coef_scatt_ir%ifmv_wcl_ph_val, phasint_ad)

                  raytracing_ad%pathsat(lay,prof) = raytracing_ad%pathsat(lay,prof) - &
                                                  & musat_ad / raytracing%pathsat(lay,prof) ** 2
                  raytracing_ad%pathsun(lay,prof) = raytracing_ad%pathsun(lay,prof) + &
                                                  & musun_ad / raytracing%pathsun(lay,prof) ** 2

                  musat_ad = 0._jprb
                  musun_ad = 0._jprb
                  phasint_ad = 0._jprb
                ENDIF !ityp
              ENDIF !solar
              IF (ityp <= 5_jpim) THEN
                transmission_scatt_ir_ad%GPAR(lay,j)           = 0._jprb
                transmission_scatt_ir_ad%GPARCLS(ityp,lay,j)  = transmission_scatt_ir_ad%GPARCLS(ityp,lay,j) +  &
                  & transmission_scatt_ir_ad%GPARTOT(lay,j) * transmission_scatt_ir%OPDPSCLS(ityp,lay,j)
                transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) = transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) + &
                  & transmission_scatt_ir_ad%GPARTOT(lay,j) * transmission_scatt_ir%GPARCLS(ityp,lay,j)
                transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) =      &
                  & transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) + transmission_scatt_ir_ad%OPDPS(lay,j)
                profiles_ad(prof)%cloud(ityp,lay)               = profiles_ad(prof)%cloud(ityp,lay) +      &
                  & transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) * coef_scatt_ir%confac(ityp) *        &
                  & optp%optpwcl(ityp)%sca(chan,1) * raytracing%ltick(lay,prof)
                raytracing_ad%ltick(lay,prof)     = raytracing_ad%ltick(lay,prof) +     &
                  & transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) * profiles(prof)%cloud(ityp,lay) *              &
                  & coef_scatt_ir%confac(ityp) * optp%optpwcl(ityp)%sca(chan,1)
                transmission_scatt_ir_ad%OPDPACLS(ityp,lay,j) =      &
                  & transmission_scatt_ir_ad%OPDPACLS(ityp,lay,j) + transmission_scatt_ir_ad%OPDPA(lay,j)
                profiles_ad(prof)%cloud(ityp,lay)               = profiles_ad(prof)%cloud(ityp,lay) +      &
                  & transmission_scatt_ir_ad%OPDPACLS(ityp,lay,j) * coef_scatt_ir%confac(ityp) *        &
                  & optp%optpwcl(ityp)%abs(chan,1) * raytracing%ltick(lay,prof)
                raytracing_ad%ltick(lay,prof)     = raytracing_ad%ltick(lay,prof) +     &
                  & transmission_scatt_ir_ad%OPDPACLS(ityp,lay,j) * profiles(prof)%cloud(ityp,lay) *              &
                  & coef_scatt_ir%confac(ityp) * optp%optpwcl(ityp)%abs(chan,1)
              ELSE
!-----------------------------------------------------------------------------------------
!                     For ice clouds optical parameters are computed using regression
!                     coefficients
!-----------------------------------------------------------------------------------------
                IF( ish .LE. 2_jpim ) THEN
                  transmission_scatt_ir_ad%GPARCLS(ityp,lay,j)  = &
                    & transmission_scatt_ir_ad%GPARCLS(ityp,lay,j) + &
                    & transmission_scatt_ir_ad%GPARTOT(lay,j) * transmission_scatt_ir%OPDPSCLS(ityp,lay,j)
                  transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) = &
                    & transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) + &
                    & transmission_scatt_ir_ad%GPARTOT(lay,j) * transmission_scatt_ir%GPARCLS(ityp,lay,j)
                  aux_ad%dg(lay,prof)                             = aux_ad%dg(lay,prof) + &
                    & transmission_scatt_ir_ad%GPARCLS(ityp,lay,j) * optp%optpicl(ish)%bpr(chan,2)
                  aux_ad%dg(lay,prof)                             = aux_ad%dg(lay,prof) + &
                    & transmission_scatt_ir_ad%GPARCLS(ityp,lay,j) *                    &
                    & 2 * aux%dg(lay,prof) * optp%optpicl(ish)%bpr(chan,3)
                  aux_ad%dg(lay,prof)                             = aux_ad%dg(lay,prof) + &
                    & transmission_scatt_ir_ad%GPARCLS(ityp,lay,j) *                    &
                    & 3 * aux%dg(lay,prof) ** 2 * optp%optpicl(ish)%bpr(chan,4)
                  transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) =      &
                    & transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) + transmission_scatt_ir_ad%OPDPS(lay,j)
                  raytracing_ad%ltick(lay,prof)     = raytracing_ad%ltick(lay,prof) + &
                    & transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) * profiles(prof)%cloud(ityp,lay) * (        &
                    & optp%optpicl(ish)%sca(chan,1) + optp%optpicl(ish)%sca(chan,2) * aux%dg(lay,prof) +          &
                    & optp%optpicl(ish)%sca(chan,3) / aux%dg(lay,prof) +                                          &
                    & optp%optpicl(ish)%sca(chan,4) / aux%dg(lay,prof) ** 2)
                  aux_ad%dg(lay,prof)                             = aux_ad%dg(lay,prof) -                    &
                    & transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) * optp%optpicl(ish)%sca(chan,4) *      &
                    & profiles(prof)%cloud(ityp,lay) * raytracing%ltick(lay,prof) * 2._jprb /  &
                    & aux%dg(lay,prof) ** 3
                  aux_ad%dg(lay,prof)                             = aux_ad%dg(lay,prof) -                &
                    & transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) * optp%optpicl(ish)%sca(chan,3) *  &
                    & profiles(prof)%cloud(ityp,lay) * raytracing%ltick(lay,prof) /        &
                    & aux%dg(lay,prof) ** 2
                  aux_ad%dg(lay,prof)                             = aux_ad%dg(lay,prof) +                &
                    & transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) * optp%optpicl(ish)%sca(chan,2) *  &
                    & profiles(prof)%cloud(ityp,lay) * raytracing%ltick(lay,prof)
                  profiles_ad(prof)%cloud(ityp,lay)               = profiles_ad(prof)%cloud(ityp,lay) +  &
                    & transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) * (                                &
                    & optp%optpicl(ish)%sca(chan,1) + optp%optpicl(ish)%sca(chan,2) * aux%dg(lay,prof) + &
                    & optp%optpicl(ish)%sca(chan,3) / aux%dg(lay,prof) +                                 &
                    & optp%optpicl(ish)%sca(chan,4) / aux%dg(lay,prof) ** 2) *                           &
                    & raytracing%ltick(lay,prof)
                  transmission_scatt_ir_ad%OPDPACLS(ityp,lay,j) =      &
                    & transmission_scatt_ir_ad%OPDPACLS(ityp,lay,j) + transmission_scatt_ir_ad%OPDPA(lay,j)
                  raytracing_ad%ltick(lay,prof)     = raytracing_ad%ltick(lay,prof) + &
                    & transmission_scatt_ir_ad%OPDPACLS(ityp,lay,j) * profiles(prof)%cloud(ityp,lay) * (        &
                    & optp%optpicl(ish)%abs(chan,1) + optp%optpicl(ish)%abs(chan,2) * aux%dg(lay,prof) +          &
                    & optp%optpicl(ish)%abs(chan,3) / aux%dg(lay,prof) +                                          &
                    & optp%optpicl(ish)%abs(chan,4) / aux%dg(lay,prof) ** 2)
                  aux_ad%dg(lay,prof)                             = aux_ad%dg(lay,prof) -                    &
                    & transmission_scatt_ir_ad%OPDPACLS(ityp,lay,j) * optp%optpicl(ish)%abs(chan,4) *      &
                    & profiles(prof)%cloud(ityp,lay) * raytracing%ltick(lay,prof) * 2._jprb /  &
                    & aux%dg(lay,prof) ** 3
                  aux_ad%dg(lay,prof)                             = aux_ad%dg(lay,prof) -                &
                    & transmission_scatt_ir_ad%OPDPACLS(ityp,lay,j) * optp%optpicl(ish)%abs(chan,3) *  &
                    & profiles(prof)%cloud(ityp,lay) * raytracing%ltick(lay,prof) /        &
                    & aux%dg(lay,prof) ** 2
                  aux_ad%dg(lay,prof)                             = aux_ad%dg(lay,prof) +                &
                    & transmission_scatt_ir_ad%OPDPACLS(ityp,lay,j) * optp%optpicl(ish)%abs(chan,2) *  &
                    & profiles(prof)%cloud(ityp,lay) * raytracing%ltick(lay,prof)
                  profiles_ad(prof)%cloud(ityp,lay)               = profiles_ad(prof)%cloud(ityp,lay) +  &
                    & transmission_scatt_ir_ad%OPDPACLS(ityp,lay,j) * (                                &
                    & optp%optpicl(ish)%abs(chan,1) + optp%optpicl(ish)%abs(chan,2) * aux%dg(lay,prof) + &
                    & optp%optpicl(ish)%abs(chan,3) / aux%dg(lay,prof) +                                 &
                    & optp%optpicl(ish)%abs(chan,4) / aux%dg(lay,prof) ** 2) *                           &
                    & raytracing%ltick(lay,prof)
                ELSE
                  !transmission_scatt_ir_tl%GPARTOT(lay,j)        = transmission_scatt_ir_tl%GPARTOT(lay,j) +        &
                  !  & transmission_scatt_ir_tl%GPARCLS(ityp,lay,j) * transmission_scatt_ir%OPDPSCLS(ityp,lay,j) + &
                  !  & transmission_scatt_ir%GPARCLS(ityp,lay,j) * transmission_scatt_ir_tl%OPDPSCLS(ityp,lay,j)

                  transmission_scatt_ir_ad%GPARCLS(ityp,lay,j)  =    &
                    & transmission_scatt_ir_ad%GPARCLS(ityp,lay,j) + &
                    & transmission_scatt_ir_ad%GPARTOT(lay,j) * transmission_scatt_ir%OPDPSCLS(ityp,lay,j)
                  transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) =     &
                    & transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) + &
                    & transmission_scatt_ir_ad%GPARTOT(lay,j) * transmission_scatt_ir%GPARCLS(ityp,lay,j)

                  !  transmission_scatt_ir_tl%GPARCLS(ityp,lay,j)  = bpr_tl
                  bpr_ad = 0.0_jprb
                  bpr_ad = bpr_ad + transmission_scatt_ir_ad%GPARCLS(ityp,lay,j)
                  !transmission_scatt_ir_ad%GPARCLS(ityp,lay,j) = 0.0_jprb


                  !  transmission_scatt_ir_tl%OPDPS(lay,j)          =      &
                  !  & transmission_scatt_ir_tl%OPDPS(lay,j) + transmission_scatt_ir_tl%OPDPSCLS(ityp,lay,j)
                  transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) = &
                      & transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) +transmission_scatt_ir_ad%OPDPS(lay,j)


                  !  transmission_scatt_ir_tl%OPDPSCLS(ityp,lay,j) = &
                  !           & sca_tl * raytracing%ltick(lay,prof) + &
                  !           & sca    * raytracing_tl%ltick(lay,prof)
                  sca_ad = 0.0_jprb
                  sca_ad = sca_ad + &
                     & transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j)* raytracing%ltick(lay,prof)
                  raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                     & sca * transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j)
                  !transmission_scatt_ir_ad%OPDPSCLS(ityp,lay,j) = 0.0_jprb

                  !  transmission_scatt_ir_tl%OPDPA(lay,j)          =      &
                  !  & transmission_scatt_ir_tl%OPDPA(lay,j) + transmission_scatt_ir_tl%OPDPACLS(ityp,lay,j)
                  transmission_scatt_ir_ad%OPDPACLS(ityp,lay,j) =     &
                    & transmission_scatt_ir_ad%OPDPACLS(ityp,lay,j) + &
                    & transmission_scatt_ir_ad%OPDPA(lay,j)

                  !  transmission_scatt_ir_tl%OPDPACLS(ityp,lay,j) = &
                  !           & abso_tl * raytracing%ltick(lay,prof) + &
                  !           & abso    * raytracing_tl%ltick(lay,prof)
                  abso_ad =0.0_jprb
                  abso_ad = abso_ad + &
                     & transmission_scatt_ir_ad%OPDPACLS(ityp,lay,j) * raytracing%ltick(lay,prof)
                  raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                    & transmission_scatt_ir_ad%OPDPACLS(ityp,lay,j) * abso
                  !transmission_scatt_ir_ad%OPDPACLS(ityp,lay,j) = 0.0_jprb

                  !  CALL rttov_baran_calc_optpar_tl (optp, chan,&
                  !       & profiles(prof)%t(lay), profiles(prof)%cloud(ityp,lay), &
                  !       & profiles_tl(prof)%t(lay), profiles_tl(prof)%cloud(ityp,lay), &
                  !       & abso, sca, bpr, asym ,&
                  !       & abso_tl, sca_tl, bpr_tl, asym_tl)
                  IF (ish == 3_JPIM) THEN
                    CALL rttov_baran2013_calc_optpar_ad (optp, chan, &
                         & profiles(prof)%t(lay), profiles(prof)%cloud(ityp,lay), &
                         & profiles_ad(prof)%t(lay), profiles_ad(prof)%cloud(ityp,lay), &
                         & abso_ad, sca_ad, bpr_ad, asym_ad)
                  ELSEIF (ish == 4_JPIM) THEN
                    CALL rttov_baran2014_calc_optpar_ad (optp, chan, &
                         & profiles(prof)%t(lay), profiles(prof)%cloud(ityp,lay), &
                         & profiles_ad(prof)%t(lay), profiles_ad(prof)%cloud(ityp,lay), &
                         & abso_ad, sca_ad, bpr_ad, asym_ad)
                  ENDIF
                  abso_ad =0.0_jprb
                  sca_ad =0.0_jprb
                  bpr_ad =0.0_jprb
                  asym_ad =0.0_jprb

                ENDIF
              ENDIF
            ENDIF ! cldtyp
          ENDDO ! ncldtyp
        ENDIF ! opts%rt_ir%user_cld_opt_param

        IF (solar(j)) THEN
          transmission_scatt_ir_ad%AZPHUP(lay,j) = 0._jprb
          transmission_scatt_ir_ad%AZPHDO(lay,j) = 0._jprb
        ENDIF
!         transmission_scatt_ir_ad%OPDPA(lay,j)  = 0._jprb
!         transmission_scatt_ir_ad%OPDPS(lay,j)  = 0._jprb
!         transmission_scatt_ir_ad%GPAR(lay,j)   = 0._jprb
      ENDDO !nlayers
    ENDDO !nchannels
    transmission_scatt_ir_ad%OPDPA       = 0._jprb
    transmission_scatt_ir_ad%OPDPS       = 0._jprb
!     transmission_scatt_ir_ad%OPDPACLS  = 0._jprb
!     transmission_scatt_ir_ad%OPDPSCLS  = 0._jprb
    transmission_scatt_ir_ad%GPAR        = 0._jprb
    transmission_scatt_ir_ad%GPARTOT     = 0._jprb
!     transmission_scatt_ir_ad%GPARCLS   = 0._jprb
    IF (dosolar) THEN
      transmission_scatt_ir_ad%AZPHUPCLS = 0._jprb
      transmission_scatt_ir_ad%AZPHDOCLS = 0._jprb
      transmission_scatt_ir_ad%AZPHUPTOT = 0._jprb
      transmission_scatt_ir_ad%AZPHDOTOT = 0._jprb
    ENDIF
  ENDIF !opts%rt_ir%addclouds
!-------------------------------------------------------------------------------
!         1.   CALCULATE OPTICAL DEPTHS OF AEROSOLS
!-------------------------------------------------------------------------------
  IF (opts%rt_ir%addaerosl) THEN
    DO J = nchannels, 1,  - 1
      chan = chanprof(j)%chan
      prof = chanprof(j)%prof
      relazi  = profiles(prof)%azangle - profiles(prof)%sunazangle
!---------Compute final values for optical parameters-------------------------------------
      DO lay = nlayers, 1,  - 1
        IF (solar(j)) THEN
          transmission_scatt_ir_ad%OPDPAERLA(lay,j) = transmission_scatt_ir_ad%OPDPAERLA(lay,j) +      &
            & OPDPAERLSUN(lay,j) * raytracing%patheff(lay,prof) * coef%ff_gam(chan)
          raytracing_ad%patheff(lay,prof)            = raytracing_ad%patheff(lay,prof) +      &
            & OPDPAERLSUN(lay,j) * transmission_scatt_ir%OPDPAERLA(lay,j) * coef%ff_gam(chan)
        ENDIF
        transmission_scatt_ir_ad%OPDPAERLA(lay,j) = transmission_scatt_ir_ad%OPDPAERLA(lay,j) +     &
          & OPDPAERL(lay,j) * raytracing%pathsat(lay,prof) * coef%ff_gam(chan)
        raytracing_ad%pathsat(lay,prof)            = raytracing_ad%pathsat(lay,prof) +     &
          & OPDPAERL(lay,j) * transmission_scatt_ir%OPDPAERLA(lay,j) * coef%ff_gam(chan)
        transmission_scatt_ir_ad%OPDPAAER(lay,j)  =      &
          & transmission_scatt_ir_ad%OPDPAAER(lay,j) + transmission_scatt_ir_ad%OPDPAERLA(lay,j)
        transmission_scatt_ir_ad%OPDPSAER(lay,j)  = transmission_scatt_ir_ad%OPDPSAER(lay,j) +      &
          & transmission_scatt_ir_ad%OPDPAERLA(lay,j) * transmission_scatt_ir%GPARAER(lay,j)
        transmission_scatt_ir_ad%GPARAER(lay,j)   = transmission_scatt_ir_ad%GPARAER(lay,j) +      &
          & transmission_scatt_ir_ad%OPDPAERLA(lay,j) * transmission_scatt_ir%OPDPSAER(lay,j)
        IF (transmission_scatt_ir%OPDPSAER(lay,j) /= 0._jprb) THEN
          IF (solar(j)) THEN
            transmission_scatt_ir_ad%AZPHAERDO(lay,j) = transmission_scatt_ir_ad%AZPHAERDO(lay,j) +      &
              & transmission_scatt_ir_ad%AZPHAERDOA(lay,j) / transmission_scatt_ir%OPDPSAER(lay,j)
            transmission_scatt_ir_ad%OPDPSAER(lay,j)  = transmission_scatt_ir_ad%OPDPSAER(lay,j) -      &
              & transmission_scatt_ir_ad%AZPHAERDOA(lay,j) * transmission_scatt_ir%AZPHAERDO(lay,j) /   &
              & transmission_scatt_ir%OPDPSAER(lay,j) ** 2
            transmission_scatt_ir_ad%AZPHAERUP(lay,j) = transmission_scatt_ir_ad%AZPHAERUP(lay,j) +      &
              & transmission_scatt_ir_ad%AZPHAERUPA(lay,j) / transmission_scatt_ir%OPDPSAER(lay,j)
            transmission_scatt_ir_ad%OPDPSAER(lay,j)  = transmission_scatt_ir_ad%OPDPSAER(lay,j) -      &
              & transmission_scatt_ir_ad%AZPHAERUPA(lay,j) * transmission_scatt_ir%AZPHAERUP(lay,j) /   &
              & transmission_scatt_ir%OPDPSAER(lay,j) ** 2
          ENDIF
          transmission_scatt_ir_ad%GPARAERA(lay,j) = transmission_scatt_ir_ad%GPARAERA(lay,j) +      &
            & transmission_scatt_ir_ad%GPARAER(lay,j) / transmission_scatt_ir%OPDPSAER(lay,j)
          transmission_scatt_ir_ad%OPDPSAER(lay,j) = transmission_scatt_ir_ad%OPDPSAER(lay,j) -      &
            & transmission_scatt_ir_ad%GPARAER(lay,j) * transmission_scatt_ir%GPARAERA(lay,j) /      &
            & transmission_scatt_ir%OPDPSAER(lay,j) ** 2
        ENDIF
      ENDDO
      DO lay = nlayers, 1,  - 1
        IF (opts%rt_ir%user_aer_opt_param) THEN
          IF (solar(j)) THEN
!-------------Average phase function for the downward scattered solar beam------------------
            musat       =  - 1._jprb / raytracing%pathsat(lay,prof)
            musun       =  - 1._jprb / raytracing%pathsun(lay,prof)
            zminphadiff = 1._jprb / (aer_opt_param%minphadiff * deg2rad)

            phasint_ad = phasint_ad + transmission_scatt_ir_ad%AZPHAERDO(lay,j) * &
                       & aer_opt_param%sca(j,lay) * raytracing%ltick(lay,prof)
            raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                       & transmission_scatt_ir_ad%AZPHAERDO(lay,j) * aer_opt_param%sca(j,lay) *        &
                       & transmission_scatt_ir%PHASINTDOREF(1,lay,j)

            CALL int_phase_fn_ad(musat, musat_ad, musun, musun_ad, relazi, zminphadiff, &
                                 aer_opt_param%pha(j,lay,:), aer_opt_param%cosphangle,  &
                                 aer_opt_param%iphangle, phasint_ad)

            musat_ad = -1._jprb * musat_ad
            phasint_ad = 0._jprb
!-------------Average phase function for the upward scattered solar beam------
            musat    = -1._jprb * musat

            phasint_ad = phasint_ad + transmission_scatt_ir_ad%AZPHAERUP(lay,j) * &
                       & aer_opt_param%sca(j,lay) * raytracing%ltick(lay,prof)
            raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                       & transmission_scatt_ir_ad%AZPHAERUP(lay,j) * aer_opt_param%sca(j,lay) *        &
                       & transmission_scatt_ir%PHASINTUPREF(1,lay,j)

            CALL int_phase_fn_ad(musat, musat_ad, musun, musun_ad, 180.0_jprb - relazi, zminphadiff, &
                                 aer_opt_param%pha(j,lay,:), aer_opt_param%cosphangle, &
                                 aer_opt_param%iphangle, phasint_ad)

            raytracing_ad%pathsat(lay,prof) = raytracing_ad%pathsat(lay,prof) - &
                                            & musat_ad / raytracing%pathsat(lay,prof) ** 2
            raytracing_ad%pathsun(lay,prof) = raytracing_ad%pathsun(lay,prof) + &
                                            & musun_ad / raytracing%pathsun(lay,prof) ** 2
            musat_ad = 0._jprb
            musun_ad = 0._jprb
            phasint_ad = 0._jprb
          ENDIF

          raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                                                      & transmission_scatt_ir_ad%GPARAERA(lay,j) *   &
                                                      & aer_opt_param%sca(j,lay) * aer_opt_param%bpr(j,lay)
          raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                                                      & transmission_scatt_ir_ad%OPDPSAER(lay,j) *   &
                                                      & aer_opt_param%sca(j,lay)
          raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                                                      & transmission_scatt_ir_ad%OPDPAAER(lay,j) *   &
                                                      & aer_opt_param%abs(j,lay)

        ELSE
          DO i = aux%iaernum(lay,prof), 1,  - 1
            iae = aux%iaertyp(i, lay,prof)
!-------------Repeat direct calculations--------------------------------------------------
            IF (coef_scatt_ir%fmv_aer_rh(iae) /= 1) THEN
              DO k = 1, coef_scatt_ir%fmv_aer_rh(iae) - 1
                IF (aux%relhum(lay,prof) >= optp%optpaer(iae)%fmv_aer_rh_val(k) .AND.      &
                  & aux%relhum(lay,prof) <= optp%optpaer(iae)%fmv_aer_rh_val(k + 1)) THEN
                  delth   = (optp%optpaer(iae)%fmv_aer_rh_val(K + 1) - optp%optpaer(iae)%fmv_aer_rh_val(K))
                  afac    = (optp%optpaer(iae)%abs(chan,k + 1) - optp%optpaer(iae)%abs(chan,k)) / delth
                  sfac    = (optp%optpaer(iae)%sca(chan,k + 1) - optp%optpaer(iae)%sca(chan,k)) / delth
                  gfac    = (optp%optpaer(iae)%bpr(chan,k + 1) - optp%optpaer(iae)%bpr(chan,k)) / delth
                  frach = (aux%relhum(lay,prof) - optp%optpaer(iae)%fmv_aer_rh_val(k))
                  absch = optp%optpaer(iae)%abs(chan,k) + afac * frach
                  scach = optp%optpaer(iae)%sca(chan,k) + sfac * frach
                  bparh = optp%optpaer(iae)%bpr(chan,k) + gfac * frach
                  IF (solar(j)) THEN
                    chan1 = coef_scatt_ir%aer_pha_index(chan)
                    pfac(1:coef_scatt_ir%fmv_aer_ph)    = (                               &
                      & optp%optpaer(iae)%pha(chan1,k + 1,1:coef_scatt_ir%fmv_aer_ph) -  &
                      & optp%optpaer(iae)%pha(chan1,k,1:coef_scatt_ir%fmv_aer_ph)) / delth
                    phash(1:coef_scatt_ir%fmv_aer_ph) = optp%optpaer(iae)%pha(chan1,k,1:coef_scatt_ir%fmv_aer_ph) + &
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
!-----------------------------------------------------------------------------------------
            IF (solar(j)) THEN
              chan1 = coef_scatt_ir%aer_pha_index(chan)
!-----------------Average phase function for the downward scattered solar beam------------
              musat       =  - 1._jprb / raytracing%pathsat(lay,prof)
              musun       =  - 1._jprb / raytracing%pathsun(lay,prof)
              zminphadiff = 1._jprb / (coef_scatt_ir%fmv_aer_ph_val_min * deg2rad)

              profiles_ad(prof)%aerosols(iae,lay)         = profiles_ad(prof)%aerosols(iae,lay) +              &
                & transmission_scatt_ir_ad%AZPHAERDO(lay,j) * transmission_scatt_ir%PHASINTDOREF(i,lay,j) * &
                & SCACH * raytracing%ltick(lay,prof)
              phasint_ad = phasint_ad +                                                                        &
                & transmission_scatt_ir_ad%AZPHAERDO(lay,j) * profiles(prof)%aerosols(iae,lay) * SCACH *      &
                & raytracing%ltick(lay,prof)
              scach_ad = scach_ad + transmission_scatt_ir_ad%AZPHAERDO(lay,j) * profiles(prof)%aerosols(iae,lay) * &
                & transmission_scatt_ir%PHASINTDOREF(i,lay,j) * raytracing%ltick(lay,prof)
              raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) +      &
                & transmission_scatt_ir_ad%AZPHAERDO(lay,j) * SCACH * profiles(prof)%aerosols(iae,lay) *      &
                & transmission_scatt_ir%PHASINTDOREF(i,lay,j)

              CALL int_phase_fn_ad(musat, musat_ad, musun, musun_ad, relazi, zminphadiff, &
                                   phash, coef_scatt_ir%fmv_aer_ph_val_cos,               &
                                   coef_scatt_ir%ifmv_aer_ph_val, phasint_ad, phash_ad)

              musat_ad = -1._jprb * musat_ad
              phasint_ad = 0._jprb
!-------------Average phase function for the upward scattered solar beam------------------
              musat    = -1._jprb * musat

              profiles_ad(prof)%aerosols(iae,lay)         = profiles_ad(prof)%aerosols(iae,lay) +              &
                & transmission_scatt_ir_ad%AZPHAERUP(lay,j) * transmission_scatt_ir%PHASINTUPREF(i,lay,j) * &
                & SCACH * raytracing%ltick(lay,prof)
              phasint_ad = phasint_ad +                                                                        &
                & transmission_scatt_ir_ad%AZPHAERUP(lay,j) * profiles(prof)%aerosols(iae,lay) * SCACH *      &
                & raytracing%ltick(lay,prof)
              scach_ad = scach_ad + transmission_scatt_ir_ad%AZPHAERUP(lay,j) * profiles(prof)%aerosols(iae,lay) * &
                & transmission_scatt_ir%PHASINTUPREF(i,lay,j) * raytracing%ltick(lay,prof)
              raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) +      &
                & transmission_scatt_ir_ad%AZPHAERUP(lay,j) * SCACH * profiles(prof)%aerosols(iae,lay) *      &
                & transmission_scatt_ir%PHASINTUPREF(i,lay,j)

              CALL int_phase_fn_ad(musat, musat_ad, musun, musun_ad, 180.0_jprb - relazi, zminphadiff, &
                                   phash, coef_scatt_ir%fmv_aer_ph_val_cos, &
                                   coef_scatt_ir%ifmv_aer_ph_val, phasint_ad, phash_ad)

              raytracing_ad%pathsat(lay,prof) = raytracing_ad%pathsat(lay,prof) - &
                                              & musat_ad / raytracing%pathsat(lay,prof) ** 2
              raytracing_ad%pathsun(lay,prof) = raytracing_ad%pathsun(lay,prof) + &
                                              & musun_ad / raytracing%pathsun(lay,prof) ** 2

              musat_ad = 0._jprb
              musun_ad = 0._jprb
              phasint_ad = 0._jprb
            ENDIF
!-----------------------------------------------------------------------------------------
            profiles_ad(prof)%aerosols(iae,lay)         = profiles_ad(prof)%aerosols(iae,lay) +      &
              & transmission_scatt_ir_ad%GPARAERA(lay,j) * SCACH * BPARH * raytracing%ltick(lay,prof)
            scach_ad = scach_ad + transmission_scatt_ir_ad%GPARAERA(lay,j) * profiles(prof)%aerosols(iae,lay) * &
              & raytracing%ltick(lay,prof) * BPARH
            raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) +      &
              & transmission_scatt_ir_ad%GPARAERA(lay,j) * SCACH * profiles(prof)%aerosols(iae,lay) * BPARH
            bparh_ad = bparh_ad + transmission_scatt_ir_ad%GPARAERA(lay,j) * profiles(prof)%aerosols(iae,lay) * &
              & raytracing%ltick(lay,prof) * SCACH
            profiles_ad(prof)%aerosols(iae,lay)         = profiles_ad(prof)%aerosols(iae,lay) +      &
              & transmission_scatt_ir_ad%OPDPSAER(lay,j) * SCACH * raytracing%ltick(lay,prof)
            scach_ad = scach_ad + transmission_scatt_ir_ad%OPDPSAER(lay,j) * profiles(prof)%aerosols(iae,lay) * &
              & raytracing%ltick(lay,prof)
            raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) +      &
              & transmission_scatt_ir_ad%OPDPSAER(lay,j) * SCACH * profiles(prof)%aerosols(iae,lay)
            profiles_ad(prof)%aerosols(iae,lay)         = profiles_ad(prof)%aerosols(iae,lay) +      &
              & transmission_scatt_ir_ad%OPDPAAER(lay,j) * ABSCH * raytracing%ltick(lay,prof)
            absch_ad = absch_ad + transmission_scatt_ir_ad%OPDPAAER(lay,j) * profiles(prof)%aerosols(iae,lay) * &
              & raytracing%ltick(lay,prof)
            raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) +      &
              & transmission_scatt_ir_ad%OPDPAAER(lay,j) * ABSCH * profiles(prof)%aerosols(iae,lay)
            IF (coef_scatt_ir%fmv_aer_rh(iae) /= 1_jpim) THEN
!---------------Interpolate scattering parameters to actual value of relative humidity--------
              DO k = coef_scatt_ir%fmv_aer_rh(iae) - 1, 1,  - 1
                IF (aux%relhum(lay,prof) >= optp%optpaer(iae)%fmv_aer_rh_val(k) .AND.      &
                  & aux%relhum(lay,prof) <= optp%optpaer(iae)%fmv_aer_rh_val(k + 1)) THEN
                  delth = (optp%optpaer(iae)%fmv_aer_rh_val(K + 1) - optp%optpaer(iae)%fmv_aer_rh_val(K))
                  afac  = (optp%optpaer(iae)%abs(chan,k + 1) - optp%optpaer(iae)%abs(chan,k)) / delth
                  sfac  = (optp%optpaer(iae)%sca(chan,k + 1) - optp%optpaer(iae)%sca(chan,k)) / delth
                  gfac  = (optp%optpaer(iae)%bpr(chan,k + 1) - optp%optpaer(iae)%bpr(chan,k)) / delth
                IF (solar(j)) THEN
                    chan1 = coef_scatt_ir%aer_pha_index(chan)
                    pfac(1:coef_scatt_ir%fmv_aer_ph) =                                         &
                      & (optp%optpaer(iae)%pha(chan1,k + 1,1:coef_scatt_ir%fmv_aer_ph) -      &
                      &  optp%optpaer(iae)%pha(chan1,k,1:coef_scatt_ir%fmv_aer_ph)) / delth
                    DO kk = coef_scatt_ir%fmv_aer_ph, 1,  - 1
                      frach_ad     = frach_ad + phash_ad(kk) * pfac(kk)
                      phash_ad(kk) = 0._jprb
                    ENDDO
                  ENDIF
                  frach_ad = frach_ad + bparh_ad * gfac
                  frach_ad = frach_ad + scach_ad * sfac
                  frach_ad = frach_ad + absch_ad * afac
                  absch_ad = 0._jprb
                  scach_ad = 0._jprb
                  bparh_ad = 0._jprb
                  aux_ad%relhum(lay,prof) = aux_ad%relhum(lay,prof) + frach_ad
                  frach_ad = 0._jprb
                  EXIT
                ENDIF
              ENDDO
            ELSE
              IF (solar(j)) THEN
                phash_ad(1:coef_scatt_ir%fmv_aer_ph) = 0._jprb
              ENDIF
              absch_ad = 0._jprb
              scach_ad = 0._jprb
              bparh_ad = 0._jprb
            ENDIF
          ENDDO ! aer types
        ENDIF ! opts%rt_ir%user_aer_opt_param
      ENDDO ! layers
    ENDDO ! channels
    transmission_scatt_ir_ad%OPDPAERLA = 0._jprb
    transmission_scatt_ir_ad%GPARAERA  = 0._jprb
    transmission_scatt_ir_ad%GPARAER   = 0._jprb
    IF (dosolar) THEN
      transmission_scatt_ir_ad%AZPHAERUP = 0._jprb
      transmission_scatt_ir_ad%AZPHAERDO = 0._jprb
    ENDIF
  ENDIF ! opts%rt_ir%addaerosl

  IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
    transmission_scatt_ir_ad%OPDPAAER = 0._jprb
    transmission_scatt_ir_ad%OPDPSAER = 0._jprb
    IF (dosolar) THEN
      transmission_scatt_ir_ad%AZPHAERUPA = 0._jprb
      transmission_scatt_ir_ad%AZPHAERDOA = 0._jprb
    ENDIF
  ENDIF
!-----Compute relative humidity-----------------------------------------------------------
  IF (opts%rt_ir%addaerosl .AND. .NOT. opts%rt_ir%user_aer_opt_param) THEN
    DO j = 1, nprofiles
      DO lay = nlayers, 1,  - 1
        IF (aux%RELHUMREF(lay,j) > 99._jprb) THEN
          aux_ad%RELHUM(lay,j) = 0._jprb
        ENDIF
        ircld_ad%wmixave(lay,j) = ircld_ad%wmixave(lay,j) +                                 &
          & aux_ad%RELHUM(lay,j) * 100._jprb * 1e-6_jprb * ircld%xpresave(lay,j) /          &
          & (ircld%PPV(lay,j) * (1._jprb + ircld%wmixave(lay,j) * 1e-6_jprb))
        ircld_ad%wmixave(lay,j) = ircld_ad%wmixave(lay,j) -                                 &
          & aux_ad%RELHUM(lay,j) * 100._jprb * 1e-6_jprb ** 2 * ircld%xpresave(lay,j) *     &
          & ircld%wmixave(lay,j) * ircld%PPV(lay,j) /                                       &
          & (ircld%PPV(lay,j) * (1._jprb + ircld%wmixave(lay,j) * 1e-6_jprb)) ** 2
        ircld_ad%PPV(lay,j)     = ircld_ad%PPV(lay,j) -                                     &
          & aux_ad%RELHUM(lay,j) * 100._jprb * ircld%wmixave(lay,j) *                       &
          & 1e-6_jprb * ircld%xpresave(lay,j) *                                             &
          & (1._jprb + ircld%wmixave(lay,j) * 1e-6_jprb) /                                  &
          & (ircld%PPV(lay,j) * (1._jprb + ircld%wmixave(lay,j) * 1e-6_jprb)) ** 2
        IF(opts%interpolation%lgradp) ircld_ad%xpresave(lay,j) = ircld_ad%xpresave(lay,j) + &
          &    100._jprb * ircld%wmixave(lay,j) * 1.e-6_jprb * aux_ad%relhum(lay,j) /       &
          &   (ircld%ppv(lay,j) * (1._jprb + ircld%wmixave(lay,j) * 1.e-6_jprb))
!---------Compute vater vapour partial pressure-------------------------------------------
        ircld_ad%PPV(lay,j)     = ircld_ad%PPV(lay,j) / 100._jprb
        IF (ircld%TAVE(lay,j) > T00) THEN
          ircld_ad%ESW(lay,j) = ircld_ad%ESW(lay,j) + ircld_ad%PPV(lay,j)
        ELSE IF (ircld%TAVE(lay,j) > TI .AND. ircld%TAVE(lay,j) <= T00) THEN
          ircld_ad%ESI(lay,j)  = ircld_ad%ESI(lay,j) + ircld_ad%PPV(lay,j)
          ircld_ad%ESW(lay,j)  =      &
            & ircld_ad%ESW(lay,j) + ircld_ad%PPV(lay,j) * ((ircld%TAVE(lay,j) - TI) / (T00 - TI)) ** 2
          ircld_ad%ESI(lay,j)  =      &
            & ircld_ad%ESI(lay,j) - ircld_ad%PPV(lay,j) * ((ircld%TAVE(lay,j) - TI) / (T00 - TI)) ** 2
          ircld_ad%TAVE(lay,j) = ircld_ad%TAVE(lay,j) +                           &
            & ircld_ad%PPV(lay,j) * (ircld%ESW(lay,j) - ircld%ESI(lay,j)) * 2 *  &
            & ((ircld%TAVE(lay,j) - TI) / (T00 - TI) ** 2)
        ELSE IF (ircld%TAVE(lay,j) <= TI) THEN
          ircld_ad%ESI(lay,j) = ircld_ad%ESI(lay,j) + ircld_ad%PPV(lay,j)
        ENDIF
        ircld_ad%TAVE(lay,j) = ircld_ad%TAVE(lay,j) +                                    &
          & ircld_ad%ESW(lay,j) * ircld%ESW(lay,j) * 17.502_jprb * (T00 - 32.19_jprb) /  &
          & (ircld%TAVE(lay,j) - 32.19_jprb) ** 2
        ircld_ad%TAVE(lay,j) = ircld_ad%TAVE(lay,j) +      &
          & ircld_ad%ESI(lay,j) * ircld%ESI(lay,j) * 22.587_jprb * (T00 + 0.7_jprb) /     &
          & (ircld%TAVE(lay,j) + 0.7_jprb) ** 2
!-----------------------------------------------------------------------------------------
      ENDDO
      DO lay = nlayers, 1,  - 1
        lev = lay + 1
        profiles_dry_ad(j)%q(lev - 1) = profiles_dry_ad(j)%q(lev - 1) + ircld_ad%wmixave(lay,j) / 2._jprb
        profiles_dry_ad(j)%q(lev)     = profiles_dry_ad(j)%q(lev) + ircld_ad%wmixave(lay,j) / 2._jprb
        profiles_ad(j)%t(lev - 1) = profiles_ad(j)%t(lev - 1) + ircld_ad%tave(lay,j) / 2._jprb
        profiles_ad(j)%t(lev)     = profiles_ad(j)%t(lev) + ircld_ad%tave(lay,j) / 2._jprb
        IF (opts%interpolation%lgradp) THEN
          profiles_ad(j)%p(lev - 1) = profiles_ad(j)%p(lev - 1) + ircld_ad%xpresave(lay,j) / 2._jprb
          profiles_ad(j)%p(lev)     = profiles_ad(j)%p(lev) + ircld_ad%xpresave(lay,j) / 2._jprb
        ENDIF
      ENDDO
    ENDDO
  ENDIF
!-----------------------------------------------------------------------------------------
  OPDPCLDL    = 0._jprb
  OPDPAERL    = 0._jprb
  OPDPCLDLSUN = 0._jprb
  OPDPAERLSUN = 0._jprb
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDPSCATTIR_AD', 1_jpim, ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE int_phase_fn_ad(musat, musat_ad, musun, musun_ad, relazi, zminphadiff, &
                             pha, cospha, ipha, phasint_ad, pha_ad)
    ! Interpolate phase function to scattering angle
    ! pha_ad may be omitted if the phase function is static
    REAL(KIND=jprb),           INTENT(IN)    :: musat, musun, relazi
    REAL(KIND=jprb),           INTENT(INOUT) :: musat_ad, musun_ad
    REAL(KIND=jprb),           INTENT(IN)    :: zminphadiff
    REAL(KIND=jprb),           INTENT(IN)    :: pha(:)
    REAL(KIND=jprb),           INTENT(IN)    :: cospha(:)
    INTEGER(KIND=jpim),        INTENT(IN)    :: ipha(:)
    REAL(KIND=jprb), OPTIONAL, INTENT(INOUT) :: pha_ad(:)
    REAL(KIND=jprb),           INTENT(INOUT) :: phasint_ad

    INTEGER(KIND=jpim) :: ikk, kk
    REAL(KIND=jprb)    :: ztmpx, ztmpx_ad
    REAL(KIND=jprb)    :: scattangle, scattangle_ad, deltap, deltap_ad, delta

    ztmpx = SQRT((1._jprb - musat ** 2) * (1._jprb - musun ** 2))

    ztmpx_ad      = 0._jprb
    scattangle_ad = 0._jprb
    deltap_ad     = 0._jprb

    scattangle    = musat * musun + ztmpx * COS(relazi * deg2rad)
    ikk           = MAX(1_jpim, INT(ACOS(scattangle) * zminphadiff, jpim))
    kk            = ipha(ikk) - 1_jpim
    deltap        = pha(kk + 1) - pha(kk)
    delta         = cospha(kk) - cospha(kk + 1)

    IF (PRESENT(pha_ad)) THEN
      pha_ad(kk) = pha_ad(kk) + phasint_ad
      deltap_ad = deltap_ad + phasint_ad * (cospha(kk) - scattangle) / delta
      scattangle_ad = scattangle_ad - deltap * phasint_ad / delta
      pha_ad(kk + 1) = pha_ad(kk + 1) + deltap_ad
      pha_ad(kk) = pha_ad(kk) - deltap_ad
    ELSE
      scattangle_ad = scattangle_ad - deltap * phasint_ad / delta
    ENDIF

    musat_ad = musat_ad + scattangle_ad * musun
    musun_ad = musun_ad + scattangle_ad * musat
    ztmpx_ad = ztmpx_ad + scattangle_ad * COS(relazi * deg2rad)

    IF (ABS(musat) == 1._jprb .OR. ABS(musun) == 1._jprb) THEN
      ztmpx_ad = 0._jprb
    ELSE
      musat_ad = musat_ad - (1._jprb - musun ** 2) * musat * ztmpx_ad / ztmpx
      musun_ad = musun_ad - (1._jprb - musat ** 2) * musun * ztmpx_ad / ztmpx
    ENDIF

  END SUBROUTINE int_phase_fn_ad

END SUBROUTINE rttov_opdpscattir_ad
