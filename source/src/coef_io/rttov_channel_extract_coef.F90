SUBROUTINE rttov_channel_extract_coef(err, coef1, coef2, channels)
! Description:
!
!   Given an rttov_coef structure containing coef1, extract
!   the data for the channels in the given list to a second
!   uninitialised coef structure.
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
!    Copyright 2014, EUMETSAT, All Rights Reserved.
!
#include "throw.h"
  USE parkind1, ONLY : jpim

  USE rttov_types, ONLY : &
    rttov_coef
!INTF_OFF
  USE rttov_const, ONLY : &
    gas_id_mixed
!INTF_ON
  IMPLICIT NONE

  INTEGER(jpim),    INTENT(OUT)   :: err
  TYPE(rttov_coef), INTENT(IN)    :: coef1
  TYPE(rttov_coef), INTENT(INOUT) :: coef2
  INTEGER(jpim),    INTENT(IN)    :: channels(:)
!INTF_END

  INTEGER(jpim)              :: i
  INTEGER(jpim)              :: nlte_start, nlte_count
  INTEGER(jpim), ALLOCATABLE :: nlte_chans(:)
! ----------------------------------------------------------------------------

TRY

  IF (MAXVAL(channels) > coef1%fmv_ori_nchn) THEN
    err = errorstatus_fatal
    THROWM(err .ne. 0, 'Channel index out of range for coefficient file')
  ENDIF
  ! Scalar variables

  ! The number of extracted channels is different...
  coef2%fmv_chn          = SIZE(channels)

  ! ... but everything else is the same
  coef2%id_platform      = coef1%id_platform
  coef2%id_sat           = coef1%id_sat
  coef2%id_inst          = coef1%id_inst
  coef2%id_sensor        = coef1%id_sensor
  coef2%id_comp_lvl      = coef1%id_comp_lvl
  coef2%id_creation_date = coef1%id_creation_date
  coef2%line_by_line     = coef1%line_by_line
  coef2%readme_srf       = coef1%readme_srf
  coef2%id_creation      = coef1%id_creation
  coef2%id_common_name   = coef1%id_common_name
  coef2%fmv_model_def    = coef1%fmv_model_def
  coef2%fmv_model_ver    = coef1%fmv_model_ver
  coef2%id_comp_pc       = coef1%id_comp_pc
  coef2%inczeeman        = coef1%inczeeman
  coef2%fmv_ori_nchn     = coef1%fmv_ori_nchn
  coef2%fmv_gas          = coef1%fmv_gas
  coef2%nmixed           = coef1%nmixed
  coef2%nwater           = coef1%nwater
  coef2%nozone           = coef1%nozone
  coef2%nwvcont          = coef1%nwvcont
  coef2%nco2             = coef1%nco2
  coef2%nn2o             = coef1%nn2o
  coef2%nco              = coef1%nco
  coef2%nch4             = coef1%nch4
  coef2%nlevels          = coef1%nlevels
  coef2%nlayers          = coef1%nlayers
  coef2%ncmixed          = coef1%ncmixed
  coef2%ncwater          = coef1%ncwater
  coef2%ncozone          = coef1%ncozone
  coef2%ncwvcont         = coef1%ncwvcont
  coef2%ncco2            = coef1%ncco2
  coef2%ncn2o            = coef1%ncn2o
  coef2%ncco             = coef1%ncco
  coef2%ncch4            = coef1%ncch4
!   coef2%nintmixed        = coef1%nintmixed
!   coef2%nintwater        = coef1%nintwater
!   coef2%nintozone        = coef1%nintozone
!   coef2%nintwvcont       = coef1%nintwvcont
!   coef2%nintco2          = coef1%nintco2
!   coef2%nintn2o          = coef1%nintn2o
!   coef2%nintco           = coef1%nintco
!   coef2%nintch4          = coef1%nintch4
  coef2%ws_nomega        = coef1%ws_nomega
  coef2%fc_planck_c1     = coef1%fc_planck_c1
  coef2%fc_planck_c2     = coef1%fc_planck_c2
  coef2%fc_sat_height    = coef1%fc_sat_height
  coef2%fastem_ver       = coef1%fastem_ver
  coef2%ssirem_ver       = coef1%ssirem_ver
  coef2%solarcoef        = coef1%solarcoef
  coef2%nltecoef         = coef1%nltecoef
  coef2%pmc_shift        = coef1%pmc_shift
  coef2%pmc_lengthcell   = coef1%pmc_lengthcell
  coef2%pmc_tempcell     = coef1%pmc_tempcell
  coef2%pmc_betaplus1    = coef1%pmc_betaplus1
  coef2%pmc_nlay         = coef1%pmc_nlay
  coef2%pmc_nvar         = coef1%pmc_nvar


  ! Fast model variables

  IF (ASSOCIATED(coef1%fmv_gas_id)) THEN

    ALLOCATE(coef2%fmv_gas_id(coef2%fmv_gas), stat=err)
    THROWM(err .ne. 0, 'allocation of fmv_chn')
    coef2%fmv_gas_id = coef1%fmv_gas_id

    ALLOCATE(coef2%fmv_gas_pos(SIZE(coef1%fmv_gas_pos)), stat=err)
    THROWM(err .ne. 0, 'allocation of fmv_gas_pos')
    coef2%fmv_gas_pos = coef1%fmv_gas_pos

    ALLOCATE(coef2%fmv_var(coef2%fmv_gas), stat=err)
    THROWM(err .ne. 0, 'allocation of fmv_var')
    coef2%fmv_var = coef1%fmv_var

    ALLOCATE(coef2%fmv_lvl(coef2%fmv_gas), stat=err)
    THROWM(err .ne. 0, 'allocation of fmv_lvl')
    coef2%fmv_lvl = coef1%fmv_lvl

    ALLOCATE(coef2%fmv_coe(coef2%fmv_gas), stat=err)
    THROWM(err .ne. 0, 'allocation of fmv_coe')
    coef2%fmv_coe = coef1%fmv_coe

  ENDIF


  ! Gas units

  IF (ASSOCIATED(coef1%gaz_units)) THEN

    ALLOCATE(coef2%gaz_units(coef2%fmv_gas), stat=err)
    THROWM(err .ne. 0, 'allocation of gaz_units')
    coef2%gaz_units = coef1%gaz_units

  ENDIF


!   ! Gas spectral interval
! 
!   IF (ASSOCIATED(coef1%mixedgasint)) THEN
! 
!     ALLOCATE(coef2%mixedgasint(2, coef2%nintmixed), stat=err)
!     THROWM(err .ne. 0, 'allocation of mixedgasint')
!     coef2%mixedgasint = coef1%mixedgasint
! 
!   ENDIF
! 
!   IF (ASSOCIATED(coef1%watervapourint)) THEN
! 
!     ALLOCATE(coef2%watervapourint(2, coef2%nintwater), stat=err)
!     THROWM(err .ne. 0, 'allocation of watervapourint')
!     coef2%watervapourint = coef1%watervapourint
! 
!   ENDIF
! 
!   IF (ASSOCIATED(coef1%ozoneint)) THEN
! 
!     ALLOCATE(coef2%ozoneint(2, coef2%nintozone), stat=err)
!     THROWM(err .ne. 0, 'allocation of ozoneint')
!     coef2%ozoneint = coef1%ozoneint
! 
!   ENDIF
! 
!   IF (ASSOCIATED(coef1%wvcontint)) THEN
! 
!     ALLOCATE(coef2%wvcontint(2, coef2%nintwvcont), stat=err)
!     THROWM(err .ne. 0, 'allocation of wvcontint')
!     coef2%wvcontint = coef1%wvcontint
! 
!   ENDIF
! 
!   IF (ASSOCIATED(coef1%co2int)) THEN
! 
!     ALLOCATE(coef2%co2int(2, coef2%nintco2), stat=err)
!     THROWM(err .ne. 0, 'allocation of co2int')
!     coef2%co2int = coef1%co2int
! 
!   ENDIF
! 
!   IF (ASSOCIATED(coef1%n2oint)) THEN
! 
!     ALLOCATE(coef2%n2oint(2, coef2%nintn2o), stat=err)
!     THROWM(err .ne. 0, 'allocation of n2oint')
!     coef2%n2oint = coef1%n2oint
! 
!   ENDIF
! 
!   IF (ASSOCIATED(coef1%coint)) THEN
! 
!     ALLOCATE(coef2%coint(2, coef2%nintco), stat=err)
!     THROWM(err .ne. 0, 'allocation of coint')
!     coef2%coint = coef1%coint
! 
!   ENDIF
! 
!   IF (ASSOCIATED(coef1%ch4int)) THEN
! 
!     ALLOCATE(coef2%ch4int(2, coef2%nintch4), stat=err)
!     THROWM(err .ne. 0, 'allocation of ch4int')
!     coef2%ch4int = coef1%ch4int
! 
!   ENDIF


  ! Filter

  IF (ASSOCIATED(coef1%ff_ori_chn)) THEN

    ALLOCATE(coef2%ff_ori_chn(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of ff_ori_chn')
    coef2%ff_ori_chn = coef1%ff_ori_chn(channels)

    ALLOCATE(coef2%ff_val_chn(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of ff_val_chn')
    coef2%ff_val_chn = coef1%ff_val_chn(channels)

    ALLOCATE(coef2%ff_cwn(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of ff_cwn')
    coef2%ff_cwn = coef1%ff_cwn(channels)

    ALLOCATE(coef2%ff_bco(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of ff_bco')
    coef2%ff_bco = coef1%ff_bco(channels)

    ALLOCATE(coef2%ff_bcs(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of ff_bcs')
    coef2%ff_bcs = coef1%ff_bcs(channels)

    ALLOCATE(coef2%ff_gam(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of ff_gam')
    coef2%ff_gam = coef1%ff_gam(channels)

  ENDIF


  ! Transmittance threshold

  IF (ASSOCIATED(coef1%tt_chn)) THEN

    ALLOCATE(coef2%tt_chn(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of tt_chn')
    coef2%tt_chn = coef1%tt_chn(channels)

    ALLOCATE(coef2%tt_val_chn(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of tt_val_chn')
    coef2%tt_val_chn = coef1%tt_val_chn(channels)

    ALLOCATE(coef2%tt_cwn(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of tt_cwn')
    coef2%tt_cwn = coef1%tt_cwn(channels)

    ALLOCATE(coef2%tt_a0(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of tt_a0')
    coef2%tt_a0 = coef1%tt_a0(channels)

    ALLOCATE(coef2%tt_a1(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of tt_a1')
    coef2%tt_a1 = coef1%tt_a1(channels)

  ENDIF


  ! Planck-weighted

  IF (ASSOCIATED(coef1%pw_chn)) THEN

    ALLOCATE(coef2%pw_chn(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of pw_chn')
    coef2%pw_chn = coef1%pw_chn(channels)

    ALLOCATE(coef2%pw_val_chn(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of pw_val_chn')
    coef2%pw_val_chn = coef1%pw_val_chn(channels)

  ENDIF


  ! Solar spectrum

  IF (ASSOCIATED(coef1%ss_chn)) THEN

    ALLOCATE(coef2%ss_chn(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of ss_chn')
    coef2%ss_chn = coef1%ss_chn(channels)

    ALLOCATE(coef2%ss_val_chn(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of ss_val_chn')
    coef2%ss_val_chn = coef1%ss_val_chn(channels)

    ALLOCATE(coef2%ss_cwn(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of ss_cwn')
    coef2%ss_cwn = coef1%ss_cwn(channels)

    ALLOCATE(coef2%ss_solar_spectrum(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of ss_solar_spectrum')
    coef2%ss_solar_spectrum = coef1%ss_solar_spectrum(channels)

  ENDIF


  ! Water optical constants

  IF (ASSOCIATED(coef1%woc_chn)) THEN

    ALLOCATE(coef2%woc_chn(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of woc_chn')
    coef2%woc_chn = coef1%woc_chn(channels)

    ALLOCATE(coef2%woc_cwn(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of woc_cwn')
    coef2%woc_cwn = coef1%woc_cwn(channels)

    ALLOCATE(coef2%woc_waopc_ow(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of woc_waopc_ow')
    coef2%woc_waopc_ow = coef1%woc_waopc_ow(channels)

    ALLOCATE(coef2%woc_waopc_fw(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of woc_waopc_fw')
    coef2%woc_waopc_fw = coef1%woc_waopc_fw(channels)

  ENDIF


  ! Wave spectrum

  IF (ASSOCIATED(coef1%ws_npoint)) THEN

    ALLOCATE(coef2%ws_npoint(coef2%ws_nomega), stat=err)
    THROWM(err .ne. 0, 'allocation of ws_npoint')
    coef2%ws_npoint = coef1%ws_npoint

    ALLOCATE(coef2%ws_k_omega(coef2%ws_nomega), stat=err)
    THROWM(err .ne. 0, 'allocation of ws_k_omega')
    coef2%ws_k_omega = coef1%ws_k_omega

  ENDIF


  ! FASTEM

  IF (ASSOCIATED(coef1%fastem_polar)) THEN

    ALLOCATE(coef2%fastem_polar(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of fastem_polar')
    coef2%fastem_polar = coef1%fastem_polar(channels)

  ENDIF


  ! SSIREM

  IF (ASSOCIATED(coef1%ssirem_chn)) THEN

    ALLOCATE(coef2%ssirem_chn(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of ssirem_chn')
    coef2%ssirem_chn = coef1%ssirem_chn(channels)

    ALLOCATE(coef2%ssirem_a0(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of ssirem_a0')
    coef2%ssirem_a0 = coef1%ssirem_a0(channels)

    ALLOCATE(coef2%ssirem_a1(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of ssirem_a1')
    coef2%ssirem_a1 = coef1%ssirem_a1(channels)

    ALLOCATE(coef2%ssirem_a2(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of ssirem_a2')
    coef2%ssirem_a2 = coef1%ssirem_a2(channels)

    ALLOCATE(coef2%ssirem_xzn1(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of ssirem_xzn1')
    coef2%ssirem_xzn1 = coef1%ssirem_xzn1(channels)

    ALLOCATE(coef2%ssirem_xzn2(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of ssirem_xzn2')
    coef2%ssirem_xzn2 = coef1%ssirem_xzn2(channels)

  ENDIF


  ! Reference profile

  IF (ASSOCIATED(coef1%ref_prfl_p)) THEN

    ALLOCATE(coef2%ref_prfl_p(coef2%fmv_lvl(gas_id_mixed)), stat=err)
    THROWM(err .ne. 0, 'allocation of ref_prfl_p')
    coef2%ref_prfl_p = coef1%ref_prfl_p

    ALLOCATE(coef2%ref_prfl_t(coef2%fmv_lvl(gas_id_mixed), coef2%fmv_gas), stat=err)
    THROWM(err .ne. 0, 'allocation of ref_prfl_t')
    coef2%ref_prfl_t = coef1%ref_prfl_t

    ALLOCATE(coef2%ref_prfl_mr(coef2%fmv_lvl(gas_id_mixed), coef2%fmv_gas), stat=err)
    THROWM(err .ne. 0, 'allocation of ref_prfl_mr')
    coef2%ref_prfl_mr = coef1%ref_prfl_mr

  ENDIF


  ! Profile limits

  IF (ASSOCIATED(coef1%lim_prfl_p)) THEN

    ALLOCATE(coef2%lim_prfl_p(coef2%fmv_lvl(gas_id_mixed)), stat=err)
    THROWM(err .ne. 0, 'allocation of lim_prfl_p')
    coef2%lim_prfl_p = coef1%lim_prfl_p

    ALLOCATE(coef2%lim_prfl_tmax(coef2%fmv_lvl(gas_id_mixed)), stat=err)
    THROWM(err .ne. 0, 'allocation of lim_prfl_tmax')
    coef2%lim_prfl_tmax = coef1%lim_prfl_tmax

    ALLOCATE(coef2%lim_prfl_tmin(coef2%fmv_lvl(gas_id_mixed)), stat=err)
    THROWM(err .ne. 0, 'allocation of lim_prfl_tmin')
    coef2%lim_prfl_tmin = coef1%lim_prfl_tmin

    ALLOCATE(coef2%lim_prfl_gmax(coef2%fmv_lvl(gas_id_mixed), coef2%fmv_gas), stat=err)
    THROWM(err .ne. 0, 'allocation of lim_prfl_gmax')
    coef2%lim_prfl_gmax = coef1%lim_prfl_gmax

    ALLOCATE(coef2%lim_prfl_gmin(coef2%fmv_lvl(gas_id_mixed), coef2%fmv_gas), stat=err)
    THROWM(err .ne. 0, 'allocation of lim_prfl_gmin')
    coef2%lim_prfl_gmin = coef1%lim_prfl_gmin

  ENDIF


  ! Thermal fast coefs

  IF (ASSOCIATED(coef1%thermal)) THEN

    ALLOCATE(coef2%thermal)
    NULLIFY (coef2%thermal%mixedgas)
    NULLIFY (coef2%thermal%watervapour)
    NULLIFY (coef2%thermal%ozone)
    NULLIFY (coef2%thermal%wvcont)
    NULLIFY (coef2%thermal%co2)
    NULLIFY (coef2%thermal%n2o)
    NULLIFY (coef2%thermal%co)
    NULLIFY (coef2%thermal%ch4)

    IF (ASSOCIATED(coef1%thermal%mixedgas)) THEN

      ALLOCATE(coef2%thermal%mixedgas(coef2%nlayers, coef2%fmv_chn, coef2%ncmixed), stat=err)
      THROWM(err .ne. 0, 'allocation of thermal%mixedgas')
      coef2%thermal%mixedgas(:,:,:) = coef1%thermal%mixedgas(:,channels,:)

    ENDIF

    IF (ASSOCIATED(coef1%thermal%watervapour)) THEN

      ALLOCATE(coef2%thermal%watervapour(coef2%nlayers, coef2%fmv_chn, coef2%ncwater), stat=err)
      THROWM(err .ne. 0, 'allocation of thermal%watervapour')
      coef2%thermal%watervapour(:,:,:) = coef1%thermal%watervapour(:,channels,:)

    ENDIF

    IF (ASSOCIATED(coef1%thermal%ozone)) THEN

      ALLOCATE(coef2%thermal%ozone(coef2%nlayers, coef2%fmv_chn, coef2%ncozone), stat=err)
      THROWM(err .ne. 0, 'allocation of thermal%ozone')
      coef2%thermal%ozone(:,:,:) = coef1%thermal%ozone(:,channels,:)

    ENDIF

    IF (ASSOCIATED(coef1%thermal%wvcont)) THEN

      ALLOCATE(coef2%thermal%wvcont(coef2%nlayers, coef2%fmv_chn, coef2%ncwvcont), stat=err)
      THROWM(err .ne. 0, 'allocation of thermal%wvcont')
      coef2%thermal%wvcont(:,:,:) = coef1%thermal%wvcont(:,channels,:)

    ENDIF

    IF (ASSOCIATED(coef1%thermal%co2)) THEN

      ALLOCATE(coef2%thermal%co2(coef2%nlayers, coef2%fmv_chn, coef2%ncco2), stat=err)
      THROWM(err .ne. 0, 'allocation of thermal%co2')
      coef2%thermal%co2(:,:,:) = coef1%thermal%co2(:,channels,:)

    ENDIF

    IF (ASSOCIATED(coef1%thermal%n2o)) THEN

      ALLOCATE(coef2%thermal%n2o(coef2%nlayers, coef2%fmv_chn, coef2%ncn2o), stat=err)
      THROWM(err .ne. 0, 'allocation of thermal%n2o')
      coef2%thermal%n2o(:,:,:) = coef1%thermal%n2o(:,channels,:)

    ENDIF

    IF (ASSOCIATED(coef1%thermal%co)) THEN

      ALLOCATE(coef2%thermal%co(coef2%nlayers, coef2%fmv_chn, coef2%ncco), stat=err)
      THROWM(err .ne. 0, 'allocation of thermal%co')
      coef2%thermal%co(:,:,:) = coef1%thermal%co(:,channels,:)

    ENDIF

    IF (ASSOCIATED(coef1%thermal%ch4)) THEN

      ALLOCATE(coef2%thermal%ch4(coef2%nlayers, coef2%fmv_chn, coef2%ncch4), stat=err)
      THROWM(err .ne. 0, 'allocation of thermal%ch4')
      coef2%thermal%ch4(:,:,:) = coef1%thermal%ch4(:,channels,:)

    ENDIF

  ENDIF


  ! Solar fast coefs

  IF (ASSOCIATED(coef1%solar)) THEN

    ALLOCATE(coef2%solar)
    NULLIFY (coef2%solar%mixedgas)
    NULLIFY (coef2%solar%watervapour)
    NULLIFY (coef2%solar%ozone)
    NULLIFY (coef2%solar%wvcont)
    NULLIFY (coef2%solar%co2)
    NULLIFY (coef2%solar%n2o)
    NULLIFY (coef2%solar%co)
    NULLIFY (coef2%solar%ch4)

    IF (ASSOCIATED(coef1%solar%mixedgas)) THEN

      ALLOCATE(coef2%solar%mixedgas(coef2%nlayers, coef2%fmv_chn, coef2%ncmixed), stat=err)
      THROWM(err .ne. 0, 'allocation of solar%mixedgas')
      coef2%solar%mixedgas(:,:,:) = coef1%solar%mixedgas(:,channels,:)

    ENDIF

    IF (ASSOCIATED(coef1%solar%watervapour)) THEN

      ALLOCATE(coef2%solar%watervapour(coef2%nlayers, coef2%fmv_chn, coef2%ncwater), stat=err)
      THROWM(err .ne. 0, 'allocation of solar%watervapour')
      coef2%solar%watervapour(:,:,:) = coef1%solar%watervapour(:,channels,:)

    ENDIF

    IF (ASSOCIATED(coef1%solar%ozone)) THEN

      ALLOCATE(coef2%solar%ozone(coef2%nlayers, coef2%fmv_chn, coef2%ncozone), stat=err)
      THROWM(err .ne. 0, 'allocation of solar%ozone')
      coef2%solar%ozone(:,:,:) = coef1%solar%ozone(:,channels,:)

    ENDIF

    IF (ASSOCIATED(coef1%solar%wvcont)) THEN

      ALLOCATE(coef2%solar%wvcont(coef2%nlayers, coef2%fmv_chn, coef2%ncwvcont), stat=err)
      THROWM(err .ne. 0, 'allocation of solar%wvcont')
      coef2%solar%wvcont(:,:,:) = coef1%solar%wvcont(:,channels,:)

    ENDIF

    IF (ASSOCIATED(coef1%solar%co2)) THEN

      ALLOCATE(coef2%solar%co2(coef2%nlayers, coef2%fmv_chn, coef2%ncco2), stat=err)
      THROWM(err .ne. 0, 'allocation of solar%co2')
      coef2%solar%co2(:,:,:) = coef1%solar%co2(:,channels,:)

    ENDIF

    IF (ASSOCIATED(coef1%solar%n2o)) THEN

      ALLOCATE(coef2%solar%n2o(coef2%nlayers, coef2%fmv_chn, coef2%ncn2o), stat=err)
      THROWM(err .ne. 0, 'allocation of solar%n2o')
      coef2%solar%n2o(:,:,:) = coef1%solar%n2o(:,channels,:)

    ENDIF

    IF (ASSOCIATED(coef1%solar%co)) THEN

      ALLOCATE(coef2%solar%co(coef2%nlayers, coef2%fmv_chn, coef2%ncco), stat=err)
      THROWM(err .ne. 0, 'allocation of solar%co')
      coef2%solar%co(:,:,:) = coef1%solar%co(:,channels,:)

    ENDIF

    IF (ASSOCIATED(coef1%solar%ch4)) THEN

      ALLOCATE(coef2%solar%ch4(coef2%nlayers, coef2%fmv_chn, coef2%ncch4), stat=err)
      THROWM(err .ne. 0, 'allocation of solar%ch4')
      coef2%solar%ch4(:,:,:) = coef1%solar%ch4(:,channels,:)

    ENDIF

  ENDIF


  ! NLTE coefs

  IF (ASSOCIATED(coef1%nlte_coef)) THEN

    ! For any monotonic channel selection we must find those selected channels
    ! which lie within the range of NLTE channels in the coef file. This
    ! constitutes another contiguous block of channels in the coef structure.

    ALLOCATE(nlte_chans(SIZE(channels))) ! Index of selected channels in nlte_coefs array in the file
    nlte_count = 0  ! Number of NLTE channels being read in
    nlte_start = 0  ! Index (in input channel list) of first NLTE channel being read in
    DO i = 1, SIZE(channels)
      IF (i > 1) THEN
        IF (channels(i) < channels(i-1)) THEN
          err = errorstatus_fatal
          THROWM( ERR .NE. 0, "non-monotonic channel selection incompatible with NLTE coefficients")
        ENDIF
      ENDIF
      IF (channels(i) >= coef1%nlte_coef%start_chan .AND. &
          channels(i) < coef1%nlte_coef%start_chan + coef1%nlte_coef%nchan) THEN
        nlte_count = nlte_count + 1
        nlte_chans(nlte_count) = channels(i) - coef1%nlte_coef%start_chan + 1
        IF (nlte_count == 1) nlte_start = i
      ENDIF
    ENDDO

    IF (nlte_count > 0) THEN
      ! We have some NLTE channels in the selection so go ahead and allocate space

      ALLOCATE(coef2%nlte_coef)

      coef2%nlte_coef%ncoef = coef1%nlte_coef%ncoef
      coef2%nlte_coef%nsol  = coef1%nlte_coef%nsol
      coef2%nlte_coef%nsat  = coef1%nlte_coef%nsat

      NULLIFY (coef2%nlte_coef%coef)
      NULLIFY (coef2%nlte_coef%sol_zen_angle, coef2%nlte_coef%cos_sol)
      NULLIFY (coef2%nlte_coef%sat_zen_angle, coef2%nlte_coef%sec_sat)

      coef2%nlte_coef%start_chan = nlte_start
      coef2%nlte_coef%nchan      = nlte_count
    ELSE

      ! No NLTE channels selected so no need for NLTE coef section
      coef2%nltecoef = .FALSE.

    ENDIF

    IF (coef2%nltecoef) THEN

      ALLOCATE(coef2%nlte_coef%sol_zen_angle(coef2%nlte_coef%nsol), stat=err)
      THROWM( ERR .NE. 0, "allocation of NLTE solar zenith angle array")
      coef2%nlte_coef%sol_zen_angle = coef1%nlte_coef%sol_zen_angle

      ALLOCATE(coef2%nlte_coef%sec_sat(coef2%nlte_coef%nsat), stat=err)
      THROWM( ERR .NE. 0, "allocation of NLTE satellite zenith angle array")
      coef2%nlte_coef%sec_sat = coef1%nlte_coef%sec_sat

      ALLOCATE(coef2%nlte_coef%coef(coef2%nlte_coef%ncoef, coef2%nlte_coef%nsat, &
               coef2%nlte_coef%nsol, coef2%nlte_coef%nchan), stat=err)
      THROWM( ERR .NE. 0, "allocation of NLTE coef array")
      coef2%nlte_coef%coef(:,:,:,:) = coef1%nlte_coef%coef(:,:,:,nlte_chans(1:nlte_count))

    ENDIF

    DEALLOCATE(nlte_chans)

  ENDIF


  ! PMC coefs

  IF (ASSOCIATED(coef1%pmc_coef)) THEN

    ALLOCATE(coef2%pmc_coef(coef2%pmc_nlay, coef2%fmv_chn, coef2%pmc_nvar), stat=err)
    THROWM(err .ne. 0, 'allocation of pmc_coef')
    coef2%pmc_coef(:,:,:) = coef1%pmc_coef(:,channels,:)

    ALLOCATE(coef2%pmc_pnominal(coef2%fmv_chn), stat=err)
    THROWM(err .ne. 0, 'allocation of pmc_pnominal')
    coef2%pmc_pnominal = coef1%pmc_pnominal(channels)

  ENDIF

CATCH
END SUBROUTINE rttov_channel_extract_coef
