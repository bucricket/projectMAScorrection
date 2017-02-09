!    Set up aerosols optical parameters for a climatological profile
SUBROUTINE rttov_opdpscattir_k( &
            & nlayers,                        &
            & chanprof,                       &
            & opts,                           &
            & aux,                            &
            & aux_k,                          &
            & profiles,                       &
            & profiles_k,                     &
            & profiles_dry_k,                 &
            & aer_opt_param,                  &
            & cld_opt_param,                  &
            & dosolar,                        &
            & solar,                          &
            & coef,                           &
            & coef_scatt_ir,                  &
            & raytracing,                     &
            & raytracing_k,                   &
            & transmission_scatt_ir,          &
            & transmission_scatt_ir_k,        &
            & transmission_scatt_ir_stream,   &
            & transmission_scatt_ir_stream_k, &
            & optp,                           &
            & ircld,                          &
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
!     1.1         05/02/2007   Removed polarisation R Saunders
!     1.2         15/09/2009   User defined ToA. Layers distinct from levels (P.Rayer)
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
  TYPE(rttov_chanprof            ), INTENT(IN)    :: chanprof  (:)
  TYPE(rttov_options )            , INTENT(IN)    :: opts
  TYPE(profile_type              ), INTENT(IN)    :: profiles  (:)
  TYPE(profile_type              ), INTENT(INOUT) :: profiles_k(SIZE(chanprof))
  TYPE(profile_type              ), INTENT(INOUT) :: profiles_dry_k(SIZE(chanprof))
  TYPE(rttov_opt_param), OPTIONAL , INTENT(IN)    :: aer_opt_param
  TYPE(rttov_opt_param), OPTIONAL , INTENT(IN)    :: cld_opt_param
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: dosolar
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: solar(SIZE(chanprof))
  TYPE(profile_aux               ), INTENT(IN)    :: aux
  TYPE(profile_aux               ), INTENT(INOUT) :: aux_k
  TYPE(rttov_coef                ), INTENT(IN)    :: coef
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_k
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir_stream
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream_k
  TYPE(rttov_coef_scatt_ir       ), INTENT(IN)    :: coef_scatt_ir
  TYPE(rttov_optpar_ir           ), INTENT(IN)    :: optp
  TYPE(ircld_type),                 INTENT(IN)    :: ircld
  TYPE(ircld_type),                 INTENT(INOUT) :: ircld_k
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing_k
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
  REAL   (KIND=jprb) :: absch_k
  REAL   (KIND=jprb) :: scach_k
  REAL   (KIND=jprb) :: bparh_k
  REAL   (KIND=jprb) :: afac
  REAL   (KIND=jprb) :: sfac
  REAL   (KIND=jprb) :: gfac
  REAL   (KIND=jprb) :: deltaice
  REAL   (KIND=jprb) :: delth
  REAL   (KIND=jprb) :: frach_k
  REAL   (KIND=jprb) :: phasint_k
  REAL   (KIND=jprb) :: musat_k
  REAL   (KIND=jprb) :: musun_k
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
  REAL   (KIND=jprb) :: phash_k    (coef_scatt_ir%fmv_aer_ph)
  REAL   (KIND=jprb) :: phasice    (coef_scatt_ir%fmv_icl_ph)
  REAL   (KIND=jprb) :: phasice_k  (coef_scatt_ir%fmv_icl_ph)
  REAL   (KIND=jprb) :: deltap     (coef_scatt_ir%fmv_icl_ph)
  INTEGER(KIND=jpim) :: nchannels
  
  REAL(KIND=jprb)  :: abso, abso_k
  REAL(KIND=jprb)  :: sca,  sca_k
  REAL(KIND=jprb)  :: bpr,  bpr_k
  REAL(KIND=jprb)  :: asym, asym_k

  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!-----End of header-----------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDPSCATTIR_K', 0_jpim, ZHOOK_HANDLE)
  nchannels        = SIZE(chanprof)
  frach_k          = 0._jprb
  absch_k          = 0._jprb
  scach_k          = 0._jprb
  bparh_k          = 0._jprb
  phasint_k        = 0._jprb
  musat_k          = 0._jprb
  musun_k          = 0._jprb
  phash_k          = 0._jprb
  phasice_k        = 0._jprb
  abso_k           = 0._jprb
  sca_k            = 0._jprb
  bpr_k            = 0._jprb
  asym_k           = 0._jprb
  OPD              = 0._jprb
  OPDPAERL         = 0._jprb
  OPDPCLDL         = 0._jprb
  IF (dosolar) THEN
    OPDSUN      = 0._jprb
    OPDPAERLSUN = 0._jprb
    OPDPCLDLSUN = 0._jprb
  ENDIF

  ! transmission_scatt_ir_k and transmission_scatt_ir_k are initialised in rttov_k

!-----Compute optical parameters for each stream------------------------------------------
  DO j = nchannels, 1,  - 1
    chan = chanprof(j)%chan
    prof = chanprof(j)%prof
    DO ist = ircld%nstream(prof), 0,  - 1
      IF (ist == 0) THEN
        IF (opts%rt_ir%addaerosl) THEN
          DO lay = nlayers, 1,  - 1
            lev = lay + 1
            OPD = OPD + transmission_scatt_ir_stream_k%OPDPAC(lev,ist,j)
            transmission_scatt_ir_k%OPDPACL(0,lay,j) = &
              & transmission_scatt_ir_k%OPDPACL(0,lay,j) + OPD
            IF (solar(j)) THEN
              OPDSUN = OPDSUN + transmission_scatt_ir_stream_k%OPDPACSUN(lev,ist,j)
              transmission_scatt_ir_k%OPDPACLSUN(0,lay,j) = &
                & transmission_scatt_ir_k%OPDPACLSUN(0,lay,j) + OPDSUN
            ENDIF
          ENDDO
          OPD = 0._jprb
          transmission_scatt_ir_stream_k%OPDPAC(1,ist,j) = 0._jprb
          IF (solar(j)) THEN
            OPDSUN = 0._jprb
            transmission_scatt_ir_stream_k%OPDPACSUN(1,ist,j) = 0._jprb
          ENDIF
        ELSE IF (opts%rt_ir%addclouds) THEN
          IF (solar(j)) THEN
            transmission_scatt_ir_stream_k%OPDPACSUN(:,ist,j) = 0._jprb
          ENDIF
          transmission_scatt_ir_stream_k%OPDPAC(:,ist,j) = 0._jprb
        ENDIF
      ELSE ! ist/=0
        DO lay = nlayers, 1,  - 1
          lev = lay + 1
          isti = ircld%icldarr(ist,lay,prof)
          OPD = OPD + transmission_scatt_ir_stream_k%OPDPAC(lev,ist,j)
          transmission_scatt_ir_k%OPDPACL(isti,lay,j) = &
            & transmission_scatt_ir_k%OPDPACL(isti,lay,j) + OPD
          IF (solar(j)) THEN
            OPDSUN = OPDSUN + transmission_scatt_ir_stream_k%OPDPACSUN(lev,ist,j)
            transmission_scatt_ir_k%OPDPACLSUN(isti,lay,j) = &
              & transmission_scatt_ir_k%OPDPACLSUN(isti,lay,j) + OPDSUN
          ENDIF
        ENDDO
        OPD = 0._jprb
        transmission_scatt_ir_stream_k%OPDPAC(1,ist,j) = 0._jprb
        IF (solar(j)) THEN
          OPDSUN = 0._jprb
          transmission_scatt_ir_stream_k%OPDPACSUN(1,ist,j) = 0._jprb
        ENDIF
      ENDIF ! ist == 0
    ENDDO ! ist

    IF (solar(j)) THEN
      DO lay = nlayers, 1, -1

        IF (opts%rt_ir%addclouds) THEN
          raytracing_k%patheff(lay,j)             = raytracing_k%patheff(lay,j) + &
            & transmission_scatt_ir_k%OPDPSCA(1,lay,j) * coef%ff_gam(chan) * &
            & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j))
          transmission_scatt_ir_k%OPDPS(lay,j)      = transmission_scatt_ir_k%OPDPS(lay,j) + &
            & transmission_scatt_ir_k%OPDPSCA(1,lay,j) * &
            & raytracing%patheff(lay,prof) * coef%ff_gam(chan)
          transmission_scatt_ir_k%OPDPSAER(lay,j)   = transmission_scatt_ir_k%OPDPSAER(lay,j) + &
            & transmission_scatt_ir_k%OPDPSCA(1,lay,j) * &
            & raytracing%patheff(lay,prof) * coef%ff_gam(chan)
          raytracing_k%patheff(lay,j)             = raytracing_k%patheff(lay,j) + &
            & transmission_scatt_ir_k%OPDPABS(1,lay,j) * coef%ff_gam(chan) * &
            & (transmission_scatt_ir%OPDPA(lay,j) + transmission_scatt_ir%OPDPAAER(lay,j))
          transmission_scatt_ir_k%OPDPA(lay,j)      = transmission_scatt_ir_k%OPDPA(lay,j) + &
            & transmission_scatt_ir_k%OPDPABS(1,lay,j) * &
            & raytracing%patheff(lay,prof) * coef%ff_gam(chan)
          transmission_scatt_ir_k%OPDPAAER(lay,j)   = transmission_scatt_ir_k%OPDPAAER(lay,j) + &
            & transmission_scatt_ir_k%OPDPABS(1,lay,j) * &
            & raytracing%patheff(lay,prof) * coef%ff_gam(chan)

          IF ((transmission_scatt_ir%OPDPS(lay,j) +      &
            & transmission_scatt_ir%OPDPSAER(lay,j)) /= 0._jprb) THEN
            transmission_scatt_ir_k%OPDPSAER(lay,j)   = transmission_scatt_ir_k%OPDPSAER(lay,j) - &
              & transmission_scatt_ir_k%AZPHACDO(1,lay,j) * transmission_scatt_ir%AZPHDO(lay,j) * &
              & transmission_scatt_ir%OPDPS(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) ** 2
            transmission_scatt_ir_k%OPDPS(lay,j)      = transmission_scatt_ir_k%OPDPS(lay,j) - &
              & transmission_scatt_ir_k%AZPHACDO(1,lay,j) * &
              & transmission_scatt_ir%AZPHDO(lay,j) * transmission_scatt_ir%OPDPS(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) ** 2
            transmission_scatt_ir_k%OPDPS(lay,j)      = transmission_scatt_ir_k%OPDPS(lay,j) + &
              & transmission_scatt_ir_k%AZPHACDO(1,lay,j) * &
              & transmission_scatt_ir%AZPHDO(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j))
            transmission_scatt_ir_k%AZPHDO(lay,j)     = transmission_scatt_ir_k%AZPHDO(lay,j) + &
              & transmission_scatt_ir_k%AZPHACDO(1,lay,j) * transmission_scatt_ir%OPDPS(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j))
            transmission_scatt_ir_k%OPDPSAER(lay,j)   = transmission_scatt_ir_k%OPDPSAER(lay,j) - &
              & transmission_scatt_ir_k%AZPHACDO(1,lay,j) * &
              & transmission_scatt_ir%AZPHAERDOA(lay,j) * transmission_scatt_ir%OPDPSAER(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) ** 2
            transmission_scatt_ir_k%OPDPS(lay,j)      = transmission_scatt_ir_k%OPDPS(lay,j) - &
              & transmission_scatt_ir_k%AZPHACDO(1,lay,j) * &
              & transmission_scatt_ir%AZPHAERDOA(lay,j) * transmission_scatt_ir%OPDPSAER(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) ** 2
            transmission_scatt_ir_k%OPDPSAER(lay,j)   = transmission_scatt_ir_k%OPDPSAER(lay,j) + &
              & transmission_scatt_ir_k%AZPHACDO(1,lay,j) * &
              & transmission_scatt_ir%AZPHAERDOA(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j))
            transmission_scatt_ir_k%AZPHAERDOA(lay,j) = transmission_scatt_ir_k%AZPHAERDOA(lay,j) + &
              & transmission_scatt_ir_k%AZPHACDO(1,lay,j) * &
              & transmission_scatt_ir%OPDPSAER(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j))
            transmission_scatt_ir_k%OPDPSAER(lay,j)   = transmission_scatt_ir_k%OPDPSAER(lay,j) - &
              & transmission_scatt_ir_k%AZPHACUP(1,lay,j) * transmission_scatt_ir%AZPHUP(lay,j) * &
              & transmission_scatt_ir%OPDPS(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) ** 2
            transmission_scatt_ir_k%OPDPS(lay,j)      = transmission_scatt_ir_k%OPDPS(lay,j) - &
              & transmission_scatt_ir_k%AZPHACUP(1,lay,j) * &
              & transmission_scatt_ir%AZPHUP(lay,j) * transmission_scatt_ir%OPDPS(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) ** 2
            transmission_scatt_ir_k%OPDPS(lay,j)      = transmission_scatt_ir_k%OPDPS(lay,j) + &
              & transmission_scatt_ir_k%AZPHACUP(1,lay,j) * &
              & transmission_scatt_ir%AZPHUP(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j))
            transmission_scatt_ir_k%AZPHUP(lay,j)     = transmission_scatt_ir_k%AZPHUP(lay,j) + &
              & transmission_scatt_ir_k%AZPHACUP(1,lay,j) * transmission_scatt_ir%OPDPS(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j))
            transmission_scatt_ir_k%OPDPSAER(lay,j)   = transmission_scatt_ir_k%OPDPSAER(lay,j) - &
              & transmission_scatt_ir_k%AZPHACUP(1,lay,j) * &
              & transmission_scatt_ir%AZPHAERUPA(lay,j) * transmission_scatt_ir%OPDPSAER(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) ** 2
            transmission_scatt_ir_k%OPDPS(lay,j)      = transmission_scatt_ir_k%OPDPS(lay,j) - &
              & transmission_scatt_ir_k%AZPHACUP(1,lay,j) * &
              & transmission_scatt_ir%AZPHAERUPA(lay,j) * transmission_scatt_ir%OPDPSAER(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j)) ** 2
            transmission_scatt_ir_k%OPDPSAER(lay,j)   = transmission_scatt_ir_k%OPDPSAER(lay,j) + &
              & transmission_scatt_ir_k%AZPHACUP(1,lay,j) * &
              & transmission_scatt_ir%AZPHAERUPA(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j))
            transmission_scatt_ir_k%AZPHAERUPA(lay,j) = transmission_scatt_ir_k%AZPHAERUPA(lay,j) + &
              & transmission_scatt_ir_k%AZPHACUP(1,lay,j) * &
              & transmission_scatt_ir%OPDPSAER(lay,j) / &
              & (transmission_scatt_ir%OPDPS(lay,j) + transmission_scatt_ir%OPDPSAER(lay,j))
          ELSE
            transmission_scatt_ir_k%AZPHACUP(1,lay,j) = 0._jprb
            transmission_scatt_ir_k%AZPHACDO(1,lay,j) = 0._jprb
          ENDIF
        ENDIF

        transmission_scatt_ir_k%OPDPSAER(lay,j) = transmission_scatt_ir_k%OPDPSAER(lay,j) + &
          & transmission_scatt_ir_k%OPDPSCA(0,lay,j) * &
          & raytracing%patheff(lay,prof) * coef%ff_gam(chan)
        raytracing_k%patheff(lay,j) = raytracing_k%patheff(lay,j) + &
          & transmission_scatt_ir_k%OPDPSCA(0,lay,j) * transmission_scatt_ir%OPDPSAER(lay,j) * &
          & coef%ff_gam(chan)
        transmission_scatt_ir_k%OPDPAAER(lay,j) = transmission_scatt_ir_k%OPDPAAER(lay,j) + &
          & transmission_scatt_ir_k%OPDPABS(0,lay,j) * &
          & raytracing%patheff(lay,prof) * coef%ff_gam(chan)
        raytracing_k%patheff(lay,j) = raytracing_k%patheff(lay,j) + &
          & transmission_scatt_ir_k%OPDPABS(0,lay,j) * transmission_scatt_ir%OPDPAAER(lay,j) * &
          & coef%ff_gam(chan)
        transmission_scatt_ir_k%AZPHAERDOA(lay,j) = &
          & transmission_scatt_ir_k%AZPHAERDOA(lay,j) + transmission_scatt_ir_k%AZPHACDO(0,lay,j)
        transmission_scatt_ir_k%AZPHAERUPA(lay,j) = &
          & transmission_scatt_ir_k%AZPHAERUPA(lay,j) + transmission_scatt_ir_k%AZPHACUP(0,lay,j)

        IF (opts%rt_ir%addclouds) THEN
          OPDPCLDLSUN(lay,j) = OPDPCLDLSUN(lay,j) + transmission_scatt_ir_k%OPDPACLSUN(1,lay,j)
          OPDPAERLSUN(lay,j) = OPDPAERLSUN(lay,j) + transmission_scatt_ir_k%OPDPACLSUN(1,lay,j)
        ENDIF
        OPDPAERLSUN(lay,j) = OPDPAERLSUN(lay,j) + transmission_scatt_ir_k%OPDPACLSUN(0,lay,j)

      ENDDO ! layers
    ENDIF ! solar channel

    DO lay = nlayers, 1, -1
      IF (opts%rt_ir%addclouds) THEN
        OPDPCLDL(lay,j) = OPDPCLDL(lay,j) + transmission_scatt_ir_k%OPDPACL(1,lay,j)
        OPDPAERL(lay,j) = OPDPAERL(lay,j) + transmission_scatt_ir_k%OPDPACL(1,lay,j)
      ENDIF
      OPDPAERL(lay,j) = OPDPAERL(lay,j) + transmission_scatt_ir_k%OPDPACL(0,lay,j)
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
          transmission_scatt_ir_k%OPDPCLDLA(lay,j) = transmission_scatt_ir_k%OPDPCLDLA(lay,j) +    &
            & OPDPCLDLSUN(lay,j) * raytracing%patheff(lay,prof) * coef%ff_gam(chan)
          raytracing_k%patheff(lay,j)              = raytracing_k%patheff(lay,j) +      &
            & OPDPCLDLSUN(lay,j) * transmission_scatt_ir%OPDPCLDLA(lay,j) * coef%ff_gam(chan)
          IF (transmission_scatt_ir%OPDPS(lay,j) /= 0._jprb) THEN
            transmission_scatt_ir_k%AZPHUPTOT(lay,j) = transmission_scatt_ir_k%AZPHUPTOT(lay,j) +  &
              & transmission_scatt_ir_k%AZPHUP(lay,j) / transmission_scatt_ir%OPDPS(lay,j)
            transmission_scatt_ir_k%OPDPS(lay,j)     = transmission_scatt_ir_k%OPDPS(lay,j) -      &
              & transmission_scatt_ir_k%AZPHUP(lay,j) * transmission_scatt_ir%AZPHUPTOT(lay,j) /   &
              & (transmission_scatt_ir%OPDPS(lay,j) ** 2)
            transmission_scatt_ir_k%AZPHDOTOT(lay,j) = transmission_scatt_ir_k%AZPHDOTOT(lay,j) +  &
              & transmission_scatt_ir_k%AZPHDO(lay,j) / transmission_scatt_ir%OPDPS(lay,j)
            transmission_scatt_ir_k%OPDPS(lay,j)     = transmission_scatt_ir_k%OPDPS(lay,j) -      &
              & transmission_scatt_ir_k%AZPHDO(lay,j) * transmission_scatt_ir%AZPHDOTOT(lay,j) /   &
              & (transmission_scatt_ir%OPDPS(lay,j) ** 2)
          ENDIF
        ENDIF
        transmission_scatt_ir_k%OPDPCLDLA(lay,j) = transmission_scatt_ir_k%OPDPCLDLA(lay,j) +  &
          & OPDPCLDL(lay,j) * raytracing%pathsat(lay,prof) * coef%ff_gam(chan)
        raytracing_k%pathsat(lay,j)              = raytracing_k%pathsat(lay,j) +               &
          & OPDPCLDL(lay,j) * transmission_scatt_ir%OPDPCLDLA(lay,j) * coef%ff_gam(chan)
        transmission_scatt_ir_k%OPDPA(lay,j)     =      &
          & transmission_scatt_ir_k%OPDPA(lay,j) + transmission_scatt_ir_k%OPDPCLDLA(lay,j)
        transmission_scatt_ir_k%OPDPS(lay,j)     = transmission_scatt_ir_k%OPDPS(lay,j) +      &
          & transmission_scatt_ir_k%OPDPCLDLA(lay,j) * transmission_scatt_ir%GPAR(lay,j)
        transmission_scatt_ir_k%GPAR(lay,j)      = transmission_scatt_ir_k%GPAR(lay,j) +       &
          & transmission_scatt_ir_k%OPDPCLDLA(lay,j) * transmission_scatt_ir%OPDPS(lay,j)
        IF (transmission_scatt_ir%OPDPS(lay,j) /= 0._jprb) THEN
          transmission_scatt_ir_k%GPARTOT(lay,j) = transmission_scatt_ir_k%GPARTOT(lay,j) +    &
            & transmission_scatt_ir_k%GPAR(lay,j) / transmission_scatt_ir%OPDPS(lay,j)
          transmission_scatt_ir_k%OPDPS(lay,j)   = transmission_scatt_ir_k%OPDPS(lay,j) -      &
            & transmission_scatt_ir_k%GPAR(lay,j) * transmission_scatt_ir%GPARTOT(lay,j) /     &
            & (transmission_scatt_ir%OPDPS(lay,j) ** 2)
        ENDIF

        IF (opts%rt_ir%user_cld_opt_param) THEN
          ! Input optical parameters are static

          IF (solar(j)) THEN
!-------------Average phase function for the downward scattered solar beam------------------
            musat       =  - 1._jprb / raytracing%pathsat(lay,prof)
            musun       =  - 1._jprb / raytracing%pathsun(lay,prof)
            zminphadiff = 1._jprb / (cld_opt_param%minphadiff * deg2rad)

            phasint_k = phasint_k + transmission_scatt_ir_k%AZPHDOTOT(lay,j) * &
                      & cld_opt_param%sca(j,lay) * raytracing%ltick(lay,prof)
            ! AZPHDOCLS(:,:,1) is the direct PHASINT
            raytracing_k%ltick(lay,j) = raytracing_k%ltick(lay,j) + &
                       & transmission_scatt_ir_k%AZPHDOTOT(lay,j) * cld_opt_param%sca(j,lay) *   &
                       & transmission_scatt_ir%AZPHDOCLS(1,lay,j)

            CALL int_phase_fn_k(musat, musat_k, musun, musun_k, relazi, zminphadiff,  &
                                cld_opt_param%pha(j,lay,:), cld_opt_param%cosphangle, &
                                cld_opt_param%iphangle, phasint_k)

            musat_k = -1._jprb * musat_k
            phasint_k = 0._jprb
!-------------Average phase function for the upward scattered solar beam------
            musat     = -1._jprb * musat

            phasint_k = phasint_k + transmission_scatt_ir_k%AZPHUPTOT(lay,j) * &
                      & cld_opt_param%sca(j,lay) * raytracing%ltick(lay,prof)
            ! AZPHUPCLS(:,:,1) is the direct PHASINT
            raytracing_k%ltick(lay,j) = raytracing_k%ltick(lay,j) + &
                       & transmission_scatt_ir_k%AZPHUPTOT(lay,j) * cld_opt_param%sca(j,lay) *   &
                       & transmission_scatt_ir%AZPHUPCLS(1,lay,j)

            CALL int_phase_fn_k(musat, musat_k, musun, musun_k, 180.0_jprb - relazi, zminphadiff, &
                                cld_opt_param%pha(j,lay,:), cld_opt_param%cosphangle, &
                                cld_opt_param%iphangle, phasint_k)

            raytracing_k%pathsat(lay,j) = raytracing_k%pathsat(lay,j) - &
                                         & musat_k / raytracing%pathsat(lay,prof) ** 2
            raytracing_k%pathsun(lay,j) = raytracing_k%pathsun(lay,j) + &
                                         & musun_k / raytracing%pathsun(lay,prof) ** 2
            musat_k = 0._jprb
            musun_k = 0._jprb
            phasint_k = 0._jprb
          ENDIF

          raytracing_k%ltick(lay,j) = raytracing_k%ltick(lay,j) +   &
                                                   & transmission_scatt_ir_k%GPARTOT(lay,j) *    &
                                                   & cld_opt_param%sca(j,lay) * cld_opt_param%bpr(j,lay)
          raytracing_k%ltick(lay,j) = raytracing_k%ltick(lay,j) +   &
                                                   & transmission_scatt_ir_k%OPDPS(lay,j) *      &
                                                   & cld_opt_param%sca(j,lay)
          raytracing_k%ltick(lay,j) = raytracing_k%ltick(lay,j) +   &
                                                   & transmission_scatt_ir_k%OPDPA(lay,j) *      &
                                                   & cld_opt_param%abs(j,lay)
        ELSE
!-----------------------------------------------------------------------------------------
!                 For water clouds use stored optical parameters
!-----------------------------------------------------------------------------------------
          DO lctyp = 1, ncldtyp
            IF (ircld%cldtyp(lctyp, lay,prof) /= 0) THEN
              ityp  = ircld%cldtyp(lctyp, lay,prof)
              IF (ityp == 6_jpim .and. ish .GE. 3_jpim ) THEN
                IF (ish == 3_JPIM) THEN 
                  CALL rttov_baran2013_calc_optpar(optp, chan, &
                       & profiles(prof)%t(lay), profiles(prof)%cloud(ityp,lay), &
                       & abso, sca, bpr, asym)
                ELSEIF (ish == 4_JPIM) THEN 
                  CALL rttov_baran2014_calc_optpar(optp, chan, &
                       & profiles(prof)%t(lay), profiles(prof)%cloud(ityp,lay), &
                       & abso, sca, bpr, asym)
                ENDIF
              ENDIF

!-----------------Average phase function for the downward scattered solar beam------------
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
                        phasice(:)   = optp%optpicl(ish)%pha(chan1, k, :) +      &
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

                  transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) =               &
                              & transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) + &
                              & transmission_scatt_ir%AZPHDOCLS(ityp,lay,j) *  &
                              & transmission_scatt_ir_k%AZPHDOTOT(lay,j)
                  transmission_scatt_ir_k%AZPHDOCLS(ityp,lay,j) =               &
                              & transmission_scatt_ir_k%AZPHDOCLS(ityp,lay,j) + &
                              & transmission_scatt_ir%OPDPSCLS(ityp,lay,j) *    &
                              & transmission_scatt_ir_k%AZPHDOTOT(lay,j)
                  phasint_k = phasint_k + transmission_scatt_ir_k%AZPHDOCLS(ityp,lay,j)

                  CALL int_phase_fn_k(musat, musat_k, musun, musun_k, relazi, zminphadiff, &
                                      phasice, coef_scatt_ir%fmv_icl_ph_val_cos,           &
                                      coef_scatt_ir%ifmv_icl_ph_val, phasint_k, phasice_k)

                  musat_k = -1._jprb * musat_k
                  phasint_k = 0._jprb
!-----------------Average phase function for the upward scattered solar beam--------------
                  musat    = -1._jprb * musat

                  transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) =               &
                              & transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) + &
                              & transmission_scatt_ir%AZPHUPCLS(ityp,lay,j) *  &
                              & transmission_scatt_ir_k%AZPHUPTOT(lay,j)
                  transmission_scatt_ir_k%AZPHUPCLS(ityp,lay,j) =               &
                              & transmission_scatt_ir_k%AZPHUPCLS(ityp,lay,j) + &
                              & transmission_scatt_ir%OPDPSCLS(ityp,lay,j) *    &
                              & transmission_scatt_ir_k%AZPHUPTOT(lay,j)
                  phasint_k = phasint_k + transmission_scatt_ir_k%AZPHUPCLS(ityp,lay,j)

                  CALL int_phase_fn_k(musat, musat_k, musun, musun_k, 180.0_jprb - relazi, zminphadiff, &
                                      phasice, coef_scatt_ir%fmv_icl_ph_val_cos, &
                                      coef_scatt_ir%ifmv_icl_ph_val, phasint_k, phasice_k)

                  raytracing_k%pathsat(lay,j) = raytracing_k%pathsat(lay,j) - &
                                               & musat_k / raytracing%pathsat(lay,prof) ** 2
                  raytracing_k%pathsun(lay,j) = raytracing_k%pathsun(lay,j) + &
                                               & musun_k / raytracing%pathsun(lay,prof) ** 2

                  musat_k = 0._jprb
                  musun_k = 0._jprb
                  phasint_k = 0._jprb
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
                          aux_k%dg(lay,j) = aux_k%dg(lay,j) + phasice_k(kk) * deltap(kk) / deltaice
                          phasice_k(kk)    = 0._jprb
                        ENDDO
                        EXIT
                      ENDIF
                    ENDDO
                  ELSE
                    asym_k = 0.0_jprb
                    CALL rttov_baran_calc_phase_ad(asym, asym_k, phangle, phasice_k)
                  ENDIF

                ELSE
!-----------------------------------------------------------------------------------------
!                     Water clouds
!-----------------------------------------------------------------------------------------
                  chan1 = coef_scatt_ir%wcl_pha_index(chan)
!-----------------------Average phase function for the downward scattered solar beam------------
                  musat       =  - 1._jprb / raytracing%pathsat(lay,prof)
                  musun       =  - 1._jprb / raytracing%pathsun(lay,prof)
                  zminphadiff = 1._jprb / (coef_scatt_ir%fmv_wcl_ph_val_min * deg2rad)

                  transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) =               &
                              & transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) + &
                              & transmission_scatt_ir%AZPHDOCLS(ityp,lay,j) *  &
                              & transmission_scatt_ir_k%AZPHDOTOT(lay,j)
                  transmission_scatt_ir_k%AZPHDOCLS(ityp,lay,j) =               &
                              & transmission_scatt_ir_k%AZPHDOCLS(ityp,lay,j) + &
                              & transmission_scatt_ir%OPDPSCLS(ityp,lay,j) *    &
                              & transmission_scatt_ir_k%AZPHDOTOT(lay,j)
                  phasint_k = phasint_k + transmission_scatt_ir_k%AZPHDOCLS(ityp,lay,j)

                  CALL int_phase_fn_k(musat, musat_k, musun, musun_k, relazi, zminphadiff,                &
                                      optp%optpwcl(ityp)%pha(chan1,1,:), coef_scatt_ir%fmv_wcl_ph_val_cos, &
                                      coef_scatt_ir%ifmv_wcl_ph_val, phasint_k)

                  musat_k = -1._jprb * musat_k
                  phasint_k = 0._jprb
!-----------------Average phase function for the upward scattered solar beam--------------
                  musat    = -1._jprb * musat

                  transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) =               &
                              & transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) + &
                              & transmission_scatt_ir%AZPHUPCLS(ityp,lay,j) *  &
                              & transmission_scatt_ir_k%AZPHUPTOT(lay,j)
                  transmission_scatt_ir_k%AZPHUPCLS(ityp,lay,j) =               &
                              & transmission_scatt_ir_k%AZPHUPCLS(ityp,lay,j) + &
                              & transmission_scatt_ir%OPDPSCLS(ityp,lay,j) *    &
                              & transmission_scatt_ir_k%AZPHUPTOT(lay,j)
                  phasint_k = phasint_k + transmission_scatt_ir_k%AZPHUPCLS(ityp,lay,j)

                  CALL int_phase_fn_k(musat, musat_k, musun, musun_k, 180.0_jprb - relazi, zminphadiff, &
                                      optp%optpwcl(ityp)%pha(chan1,1,:), coef_scatt_ir%fmv_wcl_ph_val_cos, &
                                      coef_scatt_ir%ifmv_wcl_ph_val, phasint_k)

                  raytracing_k%pathsat(lay,j) = raytracing_k%pathsat(lay,j) - &
                                               & musat_k / raytracing%pathsat(lay,prof) ** 2
                  raytracing_k%pathsun(lay,j) = raytracing_k%pathsun(lay,j) + &
                                               & musun_k / raytracing%pathsun(lay,prof) ** 2

                  musat_k = 0._jprb
                  musun_k = 0._jprb
                  phasint_k = 0._jprb
                ENDIF
              ENDIF
              IF (ityp <= 5_jpim) THEN
                transmission_scatt_ir_k%GPAR(lay,j)           = 0._jprb
                transmission_scatt_ir_k%GPARCLS(ityp,lay,j)  = transmission_scatt_ir_k%GPARCLS(ityp,lay,j) +  &
                  & transmission_scatt_ir_k%GPARTOT(lay,j) * transmission_scatt_ir%OPDPSCLS(ityp,lay,j)
                transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) = transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) + &
                  & transmission_scatt_ir_k%GPARTOT(lay,j) * transmission_scatt_ir%GPARCLS(ityp,lay,j)
                transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) =      &
                  & transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) + transmission_scatt_ir_k%OPDPS(lay,j)
                profiles_k(j)%cloud(ityp,lay)                 = profiles_k(j)%cloud(ityp,lay) +      &
                  & transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) * coef_scatt_ir%confac(ityp) *     &
                  & optp%optpwcl(ityp)%sca(chan,1) * raytracing%ltick(lay,prof)
                raytracing_k%ltick(lay,j)       = raytracing_k%ltick(lay,j) +      &
                  & transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) * profiles(prof)%cloud(ityp,lay) *            &
                  & coef_scatt_ir%confac(ityp) * optp%optpwcl(ityp)%sca(chan,1)
                transmission_scatt_ir_k%OPDPACLS(ityp,lay,j) =      &
                  & transmission_scatt_ir_k%OPDPACLS(ityp,lay,j) + transmission_scatt_ir_k%OPDPA(lay,j)
                profiles_k(j)%cloud(ityp,lay)                 = profiles_k(j)%cloud(ityp,lay) +      &
                  & transmission_scatt_ir_k%OPDPACLS(ityp,lay,j) * coef_scatt_ir%confac(ityp) *     &
                  & optp%optpwcl(ityp)%abs(chan,1) * raytracing%ltick(lay,prof)
                raytracing_k%ltick(lay,j)       = raytracing_k%ltick(lay,j) +      &
                  & transmission_scatt_ir_k%OPDPACLS(ityp,lay,j) * profiles(prof)%cloud(ityp,lay) *            &
                  & coef_scatt_ir%confac(ityp) * optp%optpwcl(ityp)%abs(chan,1)
              ELSE
!-----------------------------------------------------------------------------------------
!                 For ice clouds optical parameters are computed using regression
!                 coefficients
!-----------------------------------------------------------------------------------------
                IF( ish .LE. 2_jpim ) THEN
                  transmission_scatt_ir_k%GPARCLS(ityp,lay,j)  =     &
                    & transmission_scatt_ir_k%GPARCLS(ityp,lay,j) +  &
                    & transmission_scatt_ir_k%GPARTOT(lay,j) * transmission_scatt_ir%OPDPSCLS(ityp,lay,j)
                  transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) =     &
                    & transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) + &
                    & transmission_scatt_ir_k%GPARTOT(lay,j) * transmission_scatt_ir%GPARCLS(ityp,lay,j)
                  aux_k%dg(lay,j)                               =     &
                    & aux_k%dg(lay,j) + transmission_scatt_ir_k%GPARCLS(ityp,lay,j) * &
                    & optp%optpicl(ish)%bpr(chan,2)
                  aux_k%dg(lay,j)                               = aux_k%dg(lay,j) +      &
                    & transmission_scatt_ir_k%GPARCLS(ityp,lay,j) * 2 * aux%dg(lay,prof) * &
                    & optp%optpicl(ish)%bpr(chan,3)
                  aux_k%dg(lay,j)                               = aux_k%dg(lay,j) +      &
                    & transmission_scatt_ir_k%GPARCLS(ityp,lay,j) * 3 * aux%dg(lay,prof) ** 2 * &
                    & optp%optpicl(ish)%bpr(chan,4)
                  transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) =     &
                    & transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) + transmission_scatt_ir_k%OPDPS(lay,j)
                  raytracing_k%ltick(lay,j)       = raytracing_k%ltick(lay,j) + &
                    & transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) * profiles(prof)%cloud(ityp,lay) * (     &
                    & optp%optpicl(ish)%sca(chan,1) + optp%optpicl(ish)%sca(chan,2) * aux%dg(lay,prof) +      &
                    & optp%optpicl(ish)%sca(chan,3) / aux%dg(lay,prof) +                                      &
                    & optp%optpicl(ish)%sca(chan,4) / aux%dg(lay,prof) ** 2)
                  aux_k%dg(lay,j)                               = aux_k%dg(lay,j) -                        &
                    & transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) * optp%optpicl(ish)%sca(chan,4) *       &
                    & profiles(prof)%cloud(ityp,lay) * raytracing%ltick(lay,prof) * 2._jprb /  &
                    & aux%dg(lay,prof) ** 3
                  aux_k%dg(lay,j)                               = aux_k%dg(lay,j) -                   &
                    & transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) * optp%optpicl(ish)%sca(chan,3) *  &
                    & profiles(prof)%cloud(ityp,lay) * raytracing%ltick(lay,prof) /       &
                    & aux%dg(lay,prof) ** 2
                  aux_k%dg(lay,j)                               = aux_k%dg(lay,j) +                   &
                    & transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) * optp%optpicl(ish)%sca(chan,2) *  &
                    & profiles(prof)%cloud(ityp,lay) * raytracing%ltick(lay,prof)
                  profiles_k(j)%cloud(ityp,lay)                 = profiles_k(j)%cloud(ityp,lay) +      &
                    & transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) * (                                 &
                    & optp%optpicl(ish)%sca(chan,1) + optp%optpicl(ish)%sca(chan,2) * aux%dg(lay,prof) + &
                    & optp%optpicl(ish)%sca(chan,3) / aux%dg(lay,prof) +                                 &
                    & optp%optpicl(ish)%sca(chan,4) / aux%dg(lay,prof) ** 2) *                           &
                    & raytracing%ltick(lay,prof)
                  transmission_scatt_ir_k%OPDPACLS(ityp,lay,j) =      &
                    & transmission_scatt_ir_k%OPDPACLS(ityp,lay,j) + transmission_scatt_ir_k%OPDPA(lay,j)
                  raytracing_k%ltick(lay,j)       = raytracing_k%ltick(lay,j) + &
                    & transmission_scatt_ir_k%OPDPACLS(ityp,lay,j) * profiles(prof)%cloud(ityp,lay) * (     &
                    & optp%optpicl(ish)%abs(chan,1) + optp%optpicl(ish)%abs(chan,2) * aux%dg(lay,prof) +      &
                    & optp%optpicl(ish)%abs(chan,3) / aux%dg(lay,prof) +                                      &
                    & optp%optpicl(ish)%abs(chan,4) / aux%dg(lay,prof) ** 2)
                  aux_k%dg(lay,j)                               = aux_k%dg(lay,j) -                        &
                    & transmission_scatt_ir_k%OPDPACLS(ityp,lay,j) * optp%optpicl(ish)%abs(chan,4) *       &
                    & profiles(prof)%cloud(ityp,lay) * raytracing%ltick(lay,prof) * 2._jprb /  &
                    & aux%dg(lay,prof) ** 3
                  aux_k%dg(lay,j)                               = aux_k%dg(lay,j) -                   &
                    & transmission_scatt_ir_k%OPDPACLS(ityp,lay,j) * optp%optpicl(ish)%abs(chan,3) *  &
                    & profiles(prof)%cloud(ityp,lay) * raytracing%ltick(lay,prof) /       &
                    & aux%dg(lay,prof) ** 2
                  aux_k%dg(lay,j)                               = aux_k%dg(lay,j) +                   &
                    & transmission_scatt_ir_k%OPDPACLS(ityp,lay,j) * optp%optpicl(ish)%abs(chan,2) *  &
                    & profiles(prof)%cloud(ityp,lay) * raytracing%ltick(lay,prof)
                  profiles_k(j)%cloud(ityp,lay)                 = profiles_k(j)%cloud(ityp,lay) +     &
                    & transmission_scatt_ir_k%OPDPACLS(ityp,lay,j) * (optp%optpicl(ish)%abs(chan,1) + &
                    & optp%optpicl(ish)%abs(chan,2) * aux%dg(lay,prof) +                                &
                    & optp%optpicl(ish)%abs(chan,3) / aux%dg(lay,prof) +                                &
                    & optp%optpicl(ish)%abs(chan,4) / aux%dg(lay,prof) ** 2) *                          &
                    & raytracing%ltick(lay,prof)
                ELSE
                  transmission_scatt_ir_k%GPARCLS(ityp,lay,j)  =    &
                    & transmission_scatt_ir_k%GPARCLS(ityp,lay,j) + &
                    & transmission_scatt_ir_k%GPARTOT(lay,j) * transmission_scatt_ir%OPDPSCLS(ityp,lay,j)
                  transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) =     &
                    & transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) + &
                    & transmission_scatt_ir_k%GPARTOT(lay,j) * transmission_scatt_ir%GPARCLS(ityp,lay,j)

                  bpr_k = 0.0_jprb
                  bpr_k = bpr_k + transmission_scatt_ir_k%GPARCLS(ityp,lay,j)

                  transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) = &
                      & transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j) +transmission_scatt_ir_k%OPDPS(lay,j)

                  sca_k = 0.0_jprb
                  sca_k = sca_k + &
                     & transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j)* raytracing%ltick(lay,prof)
                  raytracing_k%ltick(lay,j) = raytracing_k%ltick(lay,j) + &
                     & sca * transmission_scatt_ir_k%OPDPSCLS(ityp,lay,j)

                  transmission_scatt_ir_k%OPDPACLS(ityp,lay,j) =     &
                    & transmission_scatt_ir_k%OPDPACLS(ityp,lay,j) + &
                    & transmission_scatt_ir_k%OPDPA(lay,j)

                  abso_k =0.0_jprb
                  abso_k= abso_k + &
                     & transmission_scatt_ir_k%OPDPACLS(ityp,lay,j) * raytracing%ltick(lay,prof)
                  raytracing_k%ltick(lay,j) = raytracing_k%ltick(lay,j) + &
                    & transmission_scatt_ir_k%OPDPACLS(ityp,lay,j) * abso

                  IF (ish == 3_JPIM) THEN 
                    CALL rttov_baran2013_calc_optpar_ad(optp, chan, &
                         & profiles(prof)%t(lay), profiles(prof)%cloud(ityp,lay), &
                         & profiles_k(j)%t(lay), profiles_k(j)%cloud(ityp,lay), &
                         & abso_k, sca_k, bpr_k, asym_k)
                  ELSEIF (ish == 4_JPIM) THEN 
                    CALL rttov_baran2014_calc_optpar_ad(optp, chan, &
                         & profiles(prof)%t(lay), profiles(prof)%cloud(ityp,lay), &
                         & profiles_k(j)%t(lay), profiles_k(j)%cloud(ityp,lay), &
                         & abso_k, sca_k, bpr_k, asym_k)
                  ENDIF
                  
                  abso_k           = 0._jprb
                  sca_k            = 0._jprb
                  bpr_k            = 0._jprb
                  asym_k           = 0._jprb

                ENDIF
              ENDIF
            ENDIF ! cldtyp
          ENDDO ! ncldtyp
        ENDIF ! opts%rt_ir%user_cld_opt_param
        IF (solar(j)) THEN
          transmission_scatt_ir_k%AZPHUP(lay,j) = 0._jprb
          transmission_scatt_ir_k%AZPHDO(lay,j) = 0._jprb
        ENDIF
!         transmission_scatt_ir_k%OPDPA(lay,j)  = 0._jprb
!         transmission_scatt_ir_k%OPDPS(lay,j)  = 0._jprb
!         transmission_scatt_ir_k%GPAR(lay,j)   = 0._jprb
      ENDDO
    ENDDO
    transmission_scatt_ir_k%OPDPA      = 0._jprb
    transmission_scatt_ir_k%OPDPS      = 0._jprb
!     transmission_scatt_ir_k%OPDPACLS = 0._jprb
!     transmission_scatt_ir_k%OPDPSCLS = 0._jprb
    transmission_scatt_ir_k%GPAR       = 0._jprb
    transmission_scatt_ir_k%GPARTOT    = 0._jprb
!     transmission_scatt_ir_k%GPARCLS  = 0._jprb
    IF (dosolar) THEN
      transmission_scatt_ir_k%AZPHUPCLS = 0._jprb
      transmission_scatt_ir_k%AZPHDOCLS = 0._jprb
      transmission_scatt_ir_k%AZPHUPTOT = 0._jprb
      transmission_scatt_ir_k%AZPHDOTOT = 0._jprb
    ENDIF
  ENDIF
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
          transmission_scatt_ir_k%OPDPAERLA(lay,j) = transmission_scatt_ir_k%OPDPAERLA(lay,j) +      &
            & OPDPAERLSUN(lay,j) * raytracing%patheff(lay,prof) * coef%ff_gam(chan)
          raytracing_k%patheff(lay,j)              = raytracing_k%patheff(lay,j) +      &
            & OPDPAERLSUN(lay,j) * transmission_scatt_ir%OPDPAERLA(lay,j) * coef%ff_gam(chan)
        ENDIF
        transmission_scatt_ir_k%OPDPAERLA(lay,j) = transmission_scatt_ir_k%OPDPAERLA(lay,j) + &
          & OPDPAERL(lay,j) * raytracing%pathsat(lay,prof) * coef%ff_gam(chan)
        raytracing_k%pathsat(lay,j)              = raytracing_k%pathsat(lay,j) +              &
          & OPDPAERL(lay,j) * transmission_scatt_ir%OPDPAERLA(lay,j) * coef%ff_gam(chan)
        transmission_scatt_ir_k%OPDPAAER(lay,j)  =      &
          & transmission_scatt_ir_k%OPDPAAER(lay,j) + transmission_scatt_ir_k%OPDPAERLA(lay,j)
        transmission_scatt_ir_k%OPDPSAER(lay,j)  = transmission_scatt_ir_k%OPDPSAER(lay,j) +  &
          & transmission_scatt_ir_k%OPDPAERLA(lay,j) * transmission_scatt_ir%GPARAER(lay,j)
        transmission_scatt_ir_k%GPARAER(lay,j)   = transmission_scatt_ir_k%GPARAER(lay,j) +   &
          & transmission_scatt_ir_k%OPDPAERLA(lay,j) * transmission_scatt_ir%OPDPSAER(lay,j)
        IF (transmission_scatt_ir%OPDPSAER(lay,j) /= 0._jprb) THEN
          IF (solar(j)) THEN
            transmission_scatt_ir_k%AZPHAERDO(lay,j) = transmission_scatt_ir_k%AZPHAERDO(lay,j) +     &
              & transmission_scatt_ir_k%AZPHAERDOA(lay,j) / transmission_scatt_ir%OPDPSAER(lay,j)
            transmission_scatt_ir_k%OPDPSAER(lay,j)  = transmission_scatt_ir_k%OPDPSAER(lay,j) -      &
              & transmission_scatt_ir_k%AZPHAERDOA(lay,j) * transmission_scatt_ir%AZPHAERDO(lay,j) /  &
              & transmission_scatt_ir%OPDPSAER(lay,j) ** 2
            transmission_scatt_ir_k%AZPHAERUP(lay,j) = transmission_scatt_ir_k%AZPHAERUP(lay,j) +     &
              & transmission_scatt_ir_k%AZPHAERUPA(lay,j) / transmission_scatt_ir%OPDPSAER(lay,j)
            transmission_scatt_ir_k%OPDPSAER(lay,j)  = transmission_scatt_ir_k%OPDPSAER(lay,j) -      &
              & transmission_scatt_ir_k%AZPHAERUPA(lay,j) * transmission_scatt_ir%AZPHAERUP(lay,j) /  &
              & transmission_scatt_ir%OPDPSAER(lay,j) ** 2
          ENDIF
          transmission_scatt_ir_k%GPARAERA(lay,j) = transmission_scatt_ir_k%GPARAERA(lay,j) +      &
            & transmission_scatt_ir_k%GPARAER(lay,j) / transmission_scatt_ir%OPDPSAER(lay,j)
          transmission_scatt_ir_k%OPDPSAER(lay,j) = transmission_scatt_ir_k%OPDPSAER(lay,j) -      &
            & transmission_scatt_ir_k%GPARAER(lay,j) * transmission_scatt_ir%GPARAERA(lay,j) /     &
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

            phasint_k = phasint_k + transmission_scatt_ir_k%AZPHAERDO(lay,j) * &
                      & aer_opt_param%sca(j,lay) * raytracing%ltick(lay,prof)
            raytracing_k%ltick(lay,j) = raytracing_k%ltick(lay,j) + &
                       & transmission_scatt_ir_k%AZPHAERDO(lay,j) * aer_opt_param%sca(j,lay) *        &
                       & transmission_scatt_ir%PHASINTDOREF(1,lay,j)

            CALL int_phase_fn_k(musat, musat_k, musun, musun_k, relazi, zminphadiff,  &
                                aer_opt_param%pha(j,lay,:), aer_opt_param%cosphangle, &
                                aer_opt_param%iphangle, phasint_k)

            musat_k = -1._jprb * musat_k
            phasint_k = 0._jprb
!-------------Average phase function for the upward scattered solar beam------
            musat    = -1._jprb * musat

            phasint_k = phasint_k + transmission_scatt_ir_k%AZPHAERUP(lay,j) * &
                      & aer_opt_param%sca(j,lay) * raytracing%ltick(lay,prof)
            raytracing_k%ltick(lay,j) = raytracing_k%ltick(lay,j) + &
                       & transmission_scatt_ir_k%AZPHAERUP(lay,j) * aer_opt_param%sca(j,lay) *        &
                       & transmission_scatt_ir%PHASINTUPREF(1,lay,j)

            CALL int_phase_fn_k(musat, musat_k, musun, musun_k, 180.0_jprb - relazi, zminphadiff, &
                                aer_opt_param%pha(j,lay,:), aer_opt_param%cosphangle, &
                                aer_opt_param%iphangle, phasint_k)

            raytracing_k%pathsat(lay,j) = raytracing_k%pathsat(lay,j) - &
                                         & musat_k / raytracing%pathsat(lay,prof) ** 2
            raytracing_k%pathsun(lay,j) = raytracing_k%pathsun(lay,j) + &
                                         & musun_k / raytracing%pathsun(lay,prof) ** 2

            musat_k = 0._jprb
            musun_k = 0._jprb
            phasint_k = 0._jprb
          ENDIF

          raytracing_k%ltick(lay,j) = raytracing_k%ltick(lay,j) +    &
                                                   & transmission_scatt_ir_k%GPARAERA(lay,j) * &
                                                   & aer_opt_param%sca(j,lay) * aer_opt_param%bpr(j,lay)
          raytracing_k%ltick(lay,j) = raytracing_k%ltick(lay,j) +    &
                                                   & transmission_scatt_ir_k%OPDPSAER(lay,j) * &
                                                   & aer_opt_param%sca(j,lay)
          raytracing_k%ltick(lay,j) = raytracing_k%ltick(lay,j) +    &
                                                   & transmission_scatt_ir_k%OPDPAAER(lay,j) * &
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
                    phash(1:coef_scatt_ir%fmv_aer_ph) =                                   &
                      & optp%optpaer(iae)%pha(chan1,k,1:coef_scatt_ir%fmv_aer_ph) +      &
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

              profiles_k(j)%aerosols(iae,lay)         = profiles_k(j)%aerosols(iae,lay) +                   &
                & transmission_scatt_ir_k%AZPHAERDO(lay,j) * transmission_scatt_ir%PHASINTDOREF(i,lay,j) * &
                & SCACH * raytracing%ltick(lay,prof)
              phasint_k = phasint_k +                                                                         &
                & transmission_scatt_ir_k%AZPHAERDO(lay,j) * profiles(prof)%aerosols(iae,lay) * SCACH *      &
                & raytracing%ltick(lay,prof)
              scach_k = scach_k + transmission_scatt_ir_k%AZPHAERDO(lay,j) * profiles(prof)%aerosols(iae,lay) * &
                & transmission_scatt_ir%PHASINTDOREF(i,lay,j) * raytracing%ltick(lay,prof)
              raytracing_k%ltick(lay,j) = raytracing_k%ltick(lay,j) +           &
                & transmission_scatt_ir_k%AZPHAERDO(lay,j) * SCACH * profiles(prof)%aerosols(iae,lay) *      &
                & transmission_scatt_ir%PHASINTDOREF(i,lay,j)

              CALL int_phase_fn_k(musat, musat_k, musun, musun_k, relazi, zminphadiff, &
                                  phash, coef_scatt_ir%fmv_aer_ph_val_cos,             &
                                  coef_scatt_ir%ifmv_aer_ph_val, phasint_k, phash_k)

              musat_k = -1._jprb * musat_k
              phasint_k = 0._jprb
!-------------Average phase function for the upward scattered solar beam------------------
              musat    = -1._jprb * musat

              profiles_k(j)%aerosols(iae,lay)         = profiles_k(j)%aerosols(iae,lay) +                   &
                & transmission_scatt_ir_k%AZPHAERUP(lay,j) * transmission_scatt_ir%PHASINTUPREF(i,lay,j) * &
                & SCACH * raytracing%ltick(lay,prof)
              phasint_k = phasint_k +                                                                         &
                & transmission_scatt_ir_k%AZPHAERUP(lay,j) * profiles(prof)%aerosols(iae,lay) * SCACH *      &
                & raytracing%ltick(lay,prof)
              scach_k = scach_k + transmission_scatt_ir_k%AZPHAERUP(lay,j) * profiles(prof)%aerosols(iae,lay) * &
                & transmission_scatt_ir%PHASINTUPREF(i,lay,j) * raytracing%ltick(lay,prof)
              raytracing_k%ltick(lay,j) = raytracing_k%ltick(lay,j) +           &
                & transmission_scatt_ir_k%AZPHAERUP(lay,j) * SCACH * profiles(prof)%aerosols(iae,lay) *      &
                & transmission_scatt_ir%PHASINTUPREF(i,lay,j)

              CALL int_phase_fn_k(musat, musat_k, musun, musun_k, 180.0_jprb - relazi, zminphadiff, &
                                  phash, coef_scatt_ir%fmv_aer_ph_val_cos, &
                                  coef_scatt_ir%ifmv_aer_ph_val, phasint_k, phash_k)

              raytracing_k%pathsat(lay,j) = raytracing_k%pathsat(lay,j) - &
                                           & musat_k / raytracing%pathsat(lay,prof) ** 2
              raytracing_k%pathsun(lay,j) = raytracing_k%pathsun(lay,j) + &
                                           & musun_k / raytracing%pathsun(lay,prof) ** 2

              musat_k = 0._jprb
              musun_k = 0._jprb
              phasint_k = 0._jprb
            ENDIF
!-----------------------------------------------------------------------------------------
            profiles_k(j)%aerosols(iae,lay)         = profiles_k(j)%aerosols(iae,lay) +      &
              & transmission_scatt_ir_k%GPARAERA(lay,j) * SCACH * BPARH * raytracing%ltick(lay,prof)
            scach_k = scach_k + transmission_scatt_ir_k%GPARAERA(lay,j) * profiles(prof)%aerosols(iae,lay) *  &
              & raytracing%ltick(lay,prof) * BPARH
            raytracing_k%ltick(lay,j) = raytracing_k%ltick(lay,j) +      &
              & transmission_scatt_ir_k%GPARAERA(lay,j) * SCACH * profiles(prof)%aerosols(iae,lay) * BPARH
            bparh_k = bparh_k + transmission_scatt_ir_k%GPARAERA(lay,j) * profiles(prof)%aerosols(iae,lay) *  &
              & raytracing%ltick(lay,prof) * SCACH
            profiles_k(j)%aerosols(iae,lay)         = profiles_k(j)%aerosols(iae,lay) +      &
              & transmission_scatt_ir_k%OPDPSAER(lay,j) * SCACH * raytracing%ltick(lay,prof)
            scach_k = scach_k + transmission_scatt_ir_k%OPDPSAER(lay,j) * profiles(prof)%aerosols(iae,lay) *  &
              & raytracing%ltick(lay,prof)
            raytracing_k%ltick(lay,j) = raytracing_k%ltick(lay,j) +      &
              & transmission_scatt_ir_k%OPDPSAER(lay,j) * SCACH * profiles(prof)%aerosols(iae,lay)
            profiles_k(j)%aerosols(iae,lay)         = profiles_k(j)%aerosols(iae,lay) +      &
              & transmission_scatt_ir_k%OPDPAAER(lay,j) * ABSCH * raytracing%ltick(lay,prof)
            absch_k = absch_k + transmission_scatt_ir_k%OPDPAAER(lay,j) * profiles(prof)%aerosols(iae,lay) *  &
              & raytracing%ltick(lay,prof)
            raytracing_k%ltick(lay,j) = raytracing_k%ltick(lay,j) +      &
              & transmission_scatt_ir_k%OPDPAAER(lay,j) * ABSCH * profiles(prof)%aerosols(iae,lay)
            IF (coef_scatt_ir%fmv_aer_rh(iae) /= 1_jpim) THEN
!---------------Interpolate scattering parameters to actual value of relative humidity-------
              DO k = coef_scatt_ir%fmv_aer_rh(iae) - 1, 1,  - 1
                IF (aux%relhum(lay,prof) >= optp%optpaer(iae)%fmv_aer_rh_val(k) .AND.      &
                  & aux%relhum(lay,prof) <= optp%optpaer(iae)%fmv_aer_rh_val(k + 1)) THEN
                  delth = (optp%optpaer(iae)%fmv_aer_rh_val(K + 1) - optp%optpaer(iae)%fmv_aer_rh_val(K))
                  afac  = (optp%optpaer(iae)%abs(chan,k + 1) - optp%optpaer(iae)%abs(chan,k)) / delth
                  sfac  = (optp%optpaer(iae)%sca(chan,k + 1) - optp%optpaer(iae)%sca(chan,k)) / delth
                  gfac  = (optp%optpaer(iae)%bpr(chan,k + 1) - optp%optpaer(iae)%bpr(chan,k)) / delth
                  IF (solar(j)) THEN
                    chan1 = coef_scatt_ir%aer_pha_index(chan)
                    pfac(1:coef_scatt_ir%fmv_aer_ph) = &
                      & (optp%optpaer(iae)%pha(chan1,k + 1,1:coef_scatt_ir%fmv_aer_ph) &
                      &  - optp%optpaer(iae)%pha(chan1,k,1:coef_scatt_ir%fmv_aer_ph)) / delth
                    DO kk = coef_scatt_ir%fmv_aer_ph, 1,  - 1
                      frach_k     = frach_k + phash_k(kk) * pfac(kk)
                      phash_k(kk) = 0._jprb
                    ENDDO
                  ENDIF
                  frach_k                = frach_k + bparh_k * gfac
                  frach_k                = frach_k + scach_k * sfac
                  frach_k                = frach_k + absch_k * afac
                  absch_k                = 0._jprb
                  scach_k                = 0._jprb
                  bparh_k                = 0._jprb
                  aux_k%relhum(lay,j) = aux_k%relhum(lay,j) + frach_k
                  frach_k                = 0._jprb
                  EXIT
                ENDIF
              ENDDO
            ELSE
              IF (solar(j)) THEN
                phash_k(1:coef_scatt_ir%fmv_aer_ph) = 0._jprb
              ENDIF
              absch_k = 0._jprb
              scach_k = 0._jprb
              bparh_k = 0._jprb
            ENDIF
          ENDDO ! aer types
        ENDIF ! opts%rt_ir%user_aer_opt_param
      ENDDO ! layers
    ENDDO ! channels
    transmission_scatt_ir_k%OPDPAERLA = 0._jprb
    transmission_scatt_ir_k%GPARAERA  = 0._jprb
    transmission_scatt_ir_k%GPARAER   = 0._jprb
    IF (dosolar) THEN
      transmission_scatt_ir_k%AZPHAERUP = 0._jprb
      transmission_scatt_ir_k%AZPHAERDO = 0._jprb
    ENDIF
  ENDIF ! opts%rt_ir%addaerosl

  IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
    transmission_scatt_ir_k%OPDPAAER = 0._jprb
    transmission_scatt_ir_k%OPDPSAER = 0._jprb
    IF (dosolar) THEN
      transmission_scatt_ir_k%AZPHAERUPA = 0._jprb
      transmission_scatt_ir_k%AZPHAERDOA = 0._jprb
    ENDIF
  ENDIF
!-----Compute relative humidity-----------------------------------------------------------
  IF (opts%rt_ir%addaerosl .AND. .NOT. opts%rt_ir%user_aer_opt_param) THEN
    DO j = 1, nchannels
      prof = chanprof(j)%prof
      DO lay = nlayers, 1,  - 1
        lev = lay + 1
        IF (aux%RELHUMREF(lay,prof) > 99._jprb) THEN
          aux_k%RELHUM(lay,j) = 0._jprb
        ENDIF
        ircld_k%wmixave(lay,j) = ircld_k%wmixave(lay,j) +                                 &
          & aux_k%RELHUM(lay,j) * 100._jprb * 1e-6_jprb * ircld%xpresave(lay,prof) /      &
          & (ircld%PPV(lay,prof) * (1._jprb + ircld%wmixave(lay,prof) * 1e-6_jprb))
        ircld_k%wmixave(lay,j) = ircld_k%wmixave(lay,j) -                                 &
          & aux_k%RELHUM(lay,j) * 100._jprb * 1e-6_jprb ** 2 * ircld%xpresave(lay,prof) * &
          & ircld%wmixave(lay,prof) * ircld%PPV(lay,prof) /                               &
          & (ircld%PPV(lay,prof) * (1._jprb + ircld%wmixave(lay,prof) * 1e-6_jprb)) ** 2
        ircld_k%PPV(lay,j)     = ircld_k%PPV(lay,j) -                                     &
          & aux_k%RELHUM(lay,j) * 100._jprb * ircld%wmixave(lay,prof) *                   &
          & 1e-6_jprb * ircld%xpresave(lay,prof) *                                        &
          & (1._jprb + ircld%wmixave(lay,prof) * 1e-6_jprb) /                             &
          & (ircld%PPV(lay,prof) * (1._jprb + ircld%wmixave(lay,prof) * 1e-6_jprb)) ** 2
        IF(opts%interpolation%lgradp) ircld_k%xpresave(lay,j) = ircld_k%xpresave(lay,j) + &
          &    100._jprb * ircld%wmixave(lay,prof) * 1.e-6_jprb * aux_k%relhum(lay,j) /   &
          &   (ircld%ppv(lay,prof) * (1._jprb + ircld%wmixave(lay,prof) * 1.e-6_jprb))
!---------Compute vater vapour partial pressure-------------------------------------------
        ircld_k%PPV(lay,j)     = ircld_k%PPV(lay,j) / 100._jprb
        IF (ircld%TAVE(lay,prof) > T00) THEN
          ircld_k%ESW(lay,j) = ircld_k%ESW(lay,j) + ircld_k%PPV(lay,j)
        ELSE IF (ircld%TAVE(lay,prof) > TI .AND. ircld%TAVE(lay,prof) <= T00) THEN
          ircld_k%ESI(lay,j)  = ircld_k%ESI(lay,j) + ircld_k%PPV(lay,j)
          ircld_k%ESW(lay,j)  =      &
            & ircld_k%ESW(lay,j) + ircld_k%PPV(lay,j) * ((ircld%TAVE(lay,prof) - TI) / (T00 - TI)) ** 2
          ircld_k%ESI(lay,j)  =      &
            & ircld_k%ESI(lay,j) - ircld_k%PPV(lay,j) * ((ircld%TAVE(lay,prof) - TI) / (T00 - TI)) ** 2
          ircld_k%TAVE(lay,j) = ircld_k%TAVE(lay,j) +                                &
            & ircld_k%PPV(lay,j) * (ircld%ESW(lay,prof) - ircld%ESI(lay,prof)) * 2 *  &
            & ((ircld%TAVE(lay,prof) - TI) / (T00 - TI) ** 2)
        ELSE IF (ircld%TAVE(lay,prof) <= TI) THEN
          ircld_k%ESI(lay,j) = ircld_k%ESI(lay,j) + ircld_k%PPV(lay,j)
        ENDIF
        ircld_k%TAVE(lay,j) = ircld_k%TAVE(lay,j) +                                       &
          & ircld_k%ESW(lay,j) * ircld%ESW(lay,prof) * 17.502_jprb * (T00 - 32.19_jprb) /  &
          & (ircld%TAVE(lay,prof) - 32.19_jprb) ** 2
        ircld_k%TAVE(lay,j) = ircld_k%TAVE(lay,j) +                                     &
          & ircld_k%ESI(lay,j) * ircld%ESI(lay,prof) * 22.587_jprb * (T00 + 0.7_jprb) /  &
          & (ircld%TAVE(lay,prof) + 0.7_jprb) ** 2
!-----------------------------------------------------------------------------------------
        profiles_dry_k(j)%q(lev - 1) = profiles_dry_k(j)%q(lev - 1) + ircld_k%wmixave(lay,j) / 2._jprb
        profiles_dry_k(j)%q(lev)     = profiles_dry_k(j)%q(lev) + ircld_k%wmixave(lay,j) / 2._jprb
        profiles_k(j)%t(lev - 1) = profiles_k(j)%t(lev - 1) + ircld_k%tave(lay,j) / 2._jprb
        profiles_k(j)%t(lev)     = profiles_k(j)%t(lev) + ircld_k%tave(lay,j) / 2._jprb
        IF (opts%interpolation%lgradp) THEN
          profiles_k(j)%p(lev - 1) = profiles_k(j)%p(lev - 1) + ircld_k%xpresave(lay,j) / 2._jprb
          profiles_k(j)%p(lev)     = profiles_k(j)%p(lev) + ircld_k%xpresave(lay,j) / 2._jprb
        ENDIF
      ENDDO
    ENDDO
  ENDIF
!-----------------------------------------------------------------------------------------
  OPDPCLDL    = 0._jprb
  OPDPAERL    = 0._jprb
  OPDPCLDLSUN = 0._jprb
  OPDPAERLSUN = 0._jprb
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDPSCATTIR_K', 1_jpim, ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE int_phase_fn_k(musat, musat_k, musun, musun_k, relazi, zminphadiff, &
                            pha, cospha, ipha, phasint_k, pha_k)
    ! Interpolate phase function to scattering angle
    ! pha_k may be omitted if the phase function is static
    REAL(KIND=jprb),           INTENT(IN)    :: musat, musun, relazi
    REAL(KIND=jprb),           INTENT(INOUT) :: musat_k, musun_k
    REAL(KIND=jprb),           INTENT(IN)    :: zminphadiff
    REAL(KIND=jprb),           INTENT(IN)    :: pha(:)
    REAL(KIND=jprb),           INTENT(IN)    :: cospha(:)
    INTEGER(KIND=jpim),        INTENT(IN)    :: ipha(:)
    REAL(KIND=jprb), OPTIONAL, INTENT(INOUT) :: pha_k(:)
    REAL(KIND=jprb),           INTENT(INOUT) :: phasint_k

    INTEGER(KIND=jpim) :: ikk, kk
    REAL(KIND=jprb)    :: ztmpx, ztmpx_k
    REAL(KIND=jprb)    :: scattangle, scattangle_k, deltap, deltap_k, delta

    ztmpx = SQRT((1._jprb - musat ** 2) * (1._jprb - musun ** 2))

    ztmpx_k      = 0._jprb
    scattangle_k = 0._jprb
    deltap_k     = 0._jprb

    scattangle    = musat * musun + ztmpx * COS(relazi * deg2rad)
    ikk           = MAX(1_jpim, INT(ACOS(scattangle) * zminphadiff, jpim))
    kk            = ipha(ikk) - 1_jpim
    deltap        = pha(kk + 1) - pha(kk)
    delta         = cospha(kk) - cospha(kk + 1)

    IF (PRESENT(pha_k)) THEN
      pha_k(kk) = pha_k(kk) + phasint_k
      deltap_k = deltap_k + phasint_k * (cospha(kk) - scattangle) / delta
      scattangle_k = scattangle_k - deltap * phasint_k / delta
      pha_k(kk + 1) = pha_k(kk + 1) + deltap_k
      pha_k(kk) = pha_k(kk) - deltap_k
    ELSE
      scattangle_k = scattangle_k - deltap * phasint_k / delta
    ENDIF

    musat_k = musat_k + scattangle_k * musun
    musun_k = musun_k + scattangle_k * musat
    ztmpx_k = ztmpx_k + scattangle_k * COS(relazi * deg2rad)

    IF (ABS(musat) == 1._jprb .OR. ABS(musun) == 1._jprb) THEN
      ztmpx_k = 0._jprb
    ELSE
      musat_k = musat_k - (1._jprb - musun ** 2) * musat * ztmpx_k / ztmpx
      musun_k = musun_k - (1._jprb - musat ** 2) * musun * ztmpx_k / ztmpx
    ENDIF

  END SUBROUTINE int_phase_fn_k

END SUBROUTINE rttov_opdpscattir_k
