SUBROUTINE rttov_channel_extract_pccoef(err, coef_pccomp1, coef_pccomp2, &
                                        channels, channels_rec)
! Description:
!
!   Given a PC coef structure, extract the data for the channels in
!   the given list to a second uninitialised PC coef structure.
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
    rttov_coef_pccomp
!INTF_OFF
  USE parkind1, ONLY : jplm
!INTF_ON
  IMPLICIT NONE

  INTEGER(jpim),           INTENT(OUT)   :: err
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp1
  TYPE(rttov_coef_pccomp), INTENT(INOUT) :: coef_pccomp2
  INTEGER(jpim),           INTENT(IN)    :: channels(:)
  INTEGER(jpim), OPTIONAL, INTENT(IN)    :: channels_rec(:)
!INTF_END

  LOGICAL(jplm) :: all_channels_rec
  INTEGER(jpim) :: m, n
! ----------------------------------------------------------------------------

TRY

  all_channels_rec = .NOT. PRESENT(channels_rec)

  ! Scalar variables

  ! A few are different...

  IF (all_channels_rec) THEN
    coef_pccomp2%fmv_pc_nchn = coef_pccomp1%fmv_pc_nchn
  ELSE
    coef_pccomp2%fmv_pc_nchn = SIZE(channels_rec)
  ENDIF

  coef_pccomp2%fmv_pc_nchn_noise = SIZE(channels)
  coef_pccomp2%fmv_pc_nche = SIZE(channels)

  ! ...but the rest are the same

  coef_pccomp2%fmv_pc_comp_pc    = coef_pccomp1%fmv_pc_comp_pc
  coef_pccomp2%fmv_pc_cld        = coef_pccomp1%fmv_pc_cld
  coef_pccomp2%fmv_pc_msets      = coef_pccomp1%fmv_pc_msets
  coef_pccomp2%fmv_pc_bands      = coef_pccomp1%fmv_pc_bands
  coef_pccomp2%fmv_pc_mnum       = coef_pccomp1%fmv_pc_mnum
  coef_pccomp2%fmv_pc_mchn       = coef_pccomp1%fmv_pc_mchn
  coef_pccomp2%fmv_pc_gas        = coef_pccomp1%fmv_pc_gas
  coef_pccomp2%fmv_pc_nlev       = coef_pccomp1%fmv_pc_nlev

  coef_pccomp2%lim_pc_prfl_pmin  = coef_pccomp1%lim_pc_prfl_pmin
  coef_pccomp2%lim_pc_prfl_pmax  = coef_pccomp1%lim_pc_prfl_pmax
  coef_pccomp2%lim_pc_prfl_tsmin = coef_pccomp1%lim_pc_prfl_tsmin
  coef_pccomp2%lim_pc_prfl_tsmax = coef_pccomp1%lim_pc_prfl_tsmax
  coef_pccomp2%lim_pc_prfl_skmin = coef_pccomp1%lim_pc_prfl_skmin
  coef_pccomp2%lim_pc_prfl_skmax = coef_pccomp1%lim_pc_prfl_skmax
  coef_pccomp2%lim_pc_prfl_wsmin = coef_pccomp1%lim_pc_prfl_wsmin
  coef_pccomp2%lim_pc_prfl_wsmax = coef_pccomp1%lim_pc_prfl_wsmax


  ! PC predictors

  ALLOCATE(coef_pccomp2%fmv_pc_sets(coef_pccomp2%fmv_pc_bands), stat=err)
  THROWM(err .ne. 0, 'allocation of fmv_pc_sets')
  coef_pccomp2%fmv_pc_sets = coef_pccomp1%fmv_pc_sets

  ALLOCATE(coef_pccomp2%pcreg(coef_pccomp2%fmv_pc_bands,coef_pccomp2%fmv_pc_msets), stat=err)
  THROWM(err .ne. 0, 'allocation of pcreg')

  DO m = 1, coef_pccomp2%fmv_pc_bands
    DO n = 1, coef_pccomp2%fmv_pc_sets(m)

      coef_pccomp2%pcreg(m,n)%fmv_pc_npred = coef_pccomp1%pcreg(m,n)%fmv_pc_npred

      ALLOCATE(coef_pccomp2%pcreg(m,n)%predictindex(coef_pccomp2%pcreg(m,n)%fmv_pc_npred), stat=err)
      THROWM(err .ne. 0, 'allocation of pcreg(m,n)%predictindex')
      coef_pccomp2%pcreg(m,n)%predictindex = coef_pccomp1%pcreg(m,n)%predictindex

    ENDDO
  ENDDO


  ! Eigenvectors

  ALLOCATE(coef_pccomp2%eigen(coef_pccomp2%fmv_pc_bands), stat=err)
  THROWM(err .ne. 0, 'allocation of eigen')

  DO m = 1, coef_pccomp2%fmv_pc_bands

    ALLOCATE(coef_pccomp2%eigen(m)%eigenvectors(coef_pccomp2%fmv_pc_nchn,coef_pccomp2%fmv_pc_mnum), stat=err)
    THROWM(err .ne. 0, 'allocation of eigen(m)%eigenvectors')

    IF (all_channels_rec) THEN
      coef_pccomp2%eigen(m)%eigenvectors = coef_pccomp1%eigen(m)%eigenvectors
    ELSE
      coef_pccomp2%eigen(m)%eigenvectors = coef_pccomp1%eigen(m)%eigenvectors(channels_rec,:)
    ENDIF

    NULLIFY(coef_pccomp2%eigen(m)%eigenvectors_t)

  ENDDO


  ! PC coefficients

  DO m = 1, coef_pccomp2%fmv_pc_bands
    DO n = 1, coef_pccomp2%fmv_pc_sets(m)

      ALLOCATE(coef_pccomp2%pcreg(m,n)%coefficients(coef_pccomp2%pcreg(m,n)%fmv_pc_npred, &
                                                    coef_pccomp2%fmv_pc_mnum), stat=err)
      THROWM(err .ne. 0, 'allocation of pcreg(m,n)%coefficients')
      coef_pccomp2%pcreg(m,n)%coefficients = coef_pccomp1%pcreg(m,n)%coefficients

      NULLIFY(coef_pccomp2%pcreg(m,n)%coefficients_t)

    ENDDO
  ENDDO


  ! Emissivity coefficients

  ALLOCATE(coef_pccomp2%emiss_chn(coef_pccomp2%fmv_pc_nche), stat=err)
  THROWM(err .ne. 0, 'allocation of emiss_chn')
  coef_pccomp2%emiss_chn = coef_pccomp1%emiss_chn(channels)

  ALLOCATE(coef_pccomp2%emiss_c1(coef_pccomp2%fmv_pc_nche), stat=err)
  THROWM(err .ne. 0, 'allocation of emiss_c1')
  coef_pccomp2%emiss_c1 = coef_pccomp1%emiss_c1(channels)

  ALLOCATE(coef_pccomp2%emiss_c2(coef_pccomp2%fmv_pc_nche), stat=err)
  THROWM(err .ne. 0, 'allocation of emiss_c2')
  coef_pccomp2%emiss_c2 = coef_pccomp1%emiss_c2(channels)

  ALLOCATE(coef_pccomp2%emiss_c3(coef_pccomp2%fmv_pc_nche), stat=err)
  THROWM(err .ne. 0, 'allocation of emiss_c3')
  coef_pccomp2%emiss_c3 = coef_pccomp1%emiss_c3(channels)

  ALLOCATE(coef_pccomp2%emiss_c4(coef_pccomp2%fmv_pc_nche), stat=err)
  THROWM(err .ne. 0, 'allocation of emiss_c4')
  coef_pccomp2%emiss_c4 = coef_pccomp1%emiss_c4(channels)

  ALLOCATE(coef_pccomp2%emiss_c5(coef_pccomp2%fmv_pc_nche), stat=err)
  THROWM(err .ne. 0, 'allocation of emiss_c5')
  coef_pccomp2%emiss_c5 = coef_pccomp1%emiss_c5(channels)

  ALLOCATE(coef_pccomp2%emiss_c6(coef_pccomp2%fmv_pc_nche), stat=err)
  THROWM(err .ne. 0, 'allocation of emiss_c6')
  coef_pccomp2%emiss_c6 = coef_pccomp1%emiss_c6(channels)

  ALLOCATE(coef_pccomp2%emiss_c7(coef_pccomp2%fmv_pc_nche), stat=err)
  THROWM(err .ne. 0, 'allocation of emiss_c7')
  coef_pccomp2%emiss_c7 = coef_pccomp1%emiss_c7(channels)

  ALLOCATE(coef_pccomp2%emiss_c8(coef_pccomp2%fmv_pc_nche), stat=err)
  THROWM(err .ne. 0, 'allocation of emiss_c8')
  coef_pccomp2%emiss_c8 = coef_pccomp1%emiss_c8(channels)

  ALLOCATE(coef_pccomp2%emiss_c9(coef_pccomp2%fmv_pc_nche), stat=err)
  THROWM(err .ne. 0, 'allocation of emiss_c9')
  coef_pccomp2%emiss_c9 = coef_pccomp1%emiss_c9(channels)


  ! PC reference profiles

  ALLOCATE(coef_pccomp2%ref_pc_prfl_p(coef_pccomp2%fmv_pc_nlev), stat=err)
  THROWM(err .ne. 0, 'allocation of ref_pc_prfl_p')
  coef_pccomp2%ref_pc_prfl_p = coef_pccomp1%ref_pc_prfl_p

  ALLOCATE(coef_pccomp2%ref_pc_prfl_mr(coef_pccomp2%fmv_pc_nlev,coef_pccomp2%fmv_pc_gas), stat=err)
  THROWM(err .ne. 0, 'allocation of ref_pc_prfl_mr')
  coef_pccomp2%ref_pc_prfl_mr = coef_pccomp1%ref_pc_prfl_mr

  ALLOCATE(coef_pccomp2%lim_pc_prfl_tmin(coef_pccomp2%fmv_pc_nlev), stat=err)
  THROWM(err .ne. 0, 'allocation of lim_pc_prfl_tmin')
  coef_pccomp2%lim_pc_prfl_tmin = coef_pccomp1%lim_pc_prfl_tmin

  ALLOCATE(coef_pccomp2%lim_pc_prfl_tmax(coef_pccomp2%fmv_pc_nlev), stat=err)
  THROWM(err .ne. 0, 'allocation of lim_pc_prfl_tmax')
  coef_pccomp2%lim_pc_prfl_tmax = coef_pccomp1%lim_pc_prfl_tmax

  ALLOCATE(coef_pccomp2%lim_pc_prfl_qmin(coef_pccomp2%fmv_pc_nlev), stat=err)
  THROWM(err .ne. 0, 'allocation of lim_pc_prfl_qmin')
  coef_pccomp2%lim_pc_prfl_qmin = coef_pccomp1%lim_pc_prfl_qmin

  ALLOCATE(coef_pccomp2%lim_pc_prfl_qmax(coef_pccomp2%fmv_pc_nlev), stat=err)
  THROWM(err .ne. 0, 'allocation of lim_pc_prfl_qmax')
  coef_pccomp2%lim_pc_prfl_qmax = coef_pccomp1%lim_pc_prfl_qmax

  ALLOCATE(coef_pccomp2%lim_pc_prfl_ozmin(coef_pccomp2%fmv_pc_nlev), stat=err)
  THROWM(err .ne. 0, 'allocation of lim_pc_prfl_ozmin')
  coef_pccomp2%lim_pc_prfl_ozmin = coef_pccomp1%lim_pc_prfl_ozmin

  ALLOCATE(coef_pccomp2%lim_pc_prfl_ozmax(coef_pccomp2%fmv_pc_nlev), stat=err)
  THROWM(err .ne. 0, 'allocation of lim_pc_prfl_ozmax')
  coef_pccomp2%lim_pc_prfl_ozmax = coef_pccomp1%lim_pc_prfl_ozmax


  ! Instrument noise

  ALLOCATE(coef_pccomp2%noise(coef_pccomp2%fmv_pc_nchn_noise), stat=err)
  THROWM(err .ne. 0, 'allocation of noise')
  coef_pccomp2%noise = coef_pccomp1%noise(channels)

  ALLOCATE(coef_pccomp2%noise_in(coef_pccomp2%fmv_pc_nchn), stat=err)
  THROWM(err .ne. 0, 'allocation of noise_in')

  ALLOCATE(coef_pccomp2%ff_ori_chn_in(coef_pccomp2%fmv_pc_nchn), stat=err)
  THROWM(err .ne. 0, 'allocation of ff_ori_chn_in')

  ALLOCATE(coef_pccomp2%ff_cwn_in(coef_pccomp2%fmv_pc_nchn), stat=err)
  THROWM(err .ne. 0, 'allocation of ff_cwn_in')

  ALLOCATE(coef_pccomp2%ff_bco_in(coef_pccomp2%fmv_pc_nchn), stat=err)
  THROWM(err .ne. 0, 'allocation of ff_bco_in')

  ALLOCATE(coef_pccomp2%ff_bcs_in(coef_pccomp2%fmv_pc_nchn), stat=err)
  THROWM(err .ne. 0, 'allocation of ff_bcs_in')

  IF (all_channels_rec) THEN
    coef_pccomp2%noise_in      = coef_pccomp1%noise_in
    coef_pccomp2%ff_ori_chn_in = coef_pccomp1%ff_ori_chn_in
    coef_pccomp2%ff_cwn_in     = coef_pccomp1%ff_cwn_in
    coef_pccomp2%ff_bco_in     = coef_pccomp1%ff_bco_in
    coef_pccomp2%ff_bcs_in     = coef_pccomp1%ff_bcs_in
  ELSE
    coef_pccomp2%noise_in      = coef_pccomp1%noise_in(channels_rec)
    coef_pccomp2%ff_ori_chn_in = coef_pccomp1%ff_ori_chn_in(channels_rec)
    coef_pccomp2%ff_cwn_in     = coef_pccomp1%ff_cwn_in(channels_rec)
    coef_pccomp2%ff_bco_in     = coef_pccomp1%ff_bco_in(channels_rec)
    coef_pccomp2%ff_bcs_in     = coef_pccomp1%ff_bcs_in(channels_rec)
  ENDIF

CATCH
END SUBROUTINE rttov_channel_extract_pccoef