SUBROUTINE rttov_channel_extract_scaercoef( &
              err,            &
              coef_scatt_ir1, &
              optp1,          &
              coef_scatt_ir2, &
              optp2,          &
              channels)
! Description:
!
!   Given an aerosol coef structure, extract the data
!   for the channels in the given list to a second
!   uninitialised aerosol coef structure.
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
    rttov_coef_scatt_ir,  &
    rttov_optpar_ir
!INTF_ON
  IMPLICIT NONE

  INTEGER(jpim),             INTENT(OUT)   :: err
  TYPE(rttov_coef_scatt_ir), INTENT(IN)    :: coef_scatt_ir1
  TYPE(rttov_optpar_ir),     INTENT(IN)    :: optp1
  TYPE(rttov_coef_scatt_ir), INTENT(INOUT) :: coef_scatt_ir2
  TYPE(rttov_optpar_ir),     INTENT(INOUT) :: optp2
  INTEGER(jpim),             INTENT(IN)    :: channels(:)
!INTF_END

#include "rttov_nullify_coef_scatt_ir.interface"

  INTEGER(jpim)              :: i, j, n
  INTEGER(jpim), ALLOCATABLE :: phase_channels(:)
  INTEGER(jpim), ALLOCATABLE :: phase_channels_ext(:)
  INTEGER(jpim), ALLOCATABLE :: phase_ext_index(:)
! ----------------------------------------------------------------------------

TRY

  ! Scalar variables

  ! A few of these are the same...
  coef_scatt_ir2%fmv_aer_comp       = coef_scatt_ir1%fmv_aer_comp
  coef_scatt_ir2%fmv_aer_ph         = coef_scatt_ir1%fmv_aer_ph
  coef_scatt_ir2%fmv_aer_ph_val_min = coef_scatt_ir1%fmv_aer_ph_val_min

  ! ... but we have to work out the solar channel numbers
  coef_scatt_ir2%fmv_aer_chn      = SIZE(channels)
  coef_scatt_ir2%fmv_aer_pha_ioff = 0   ! Don't use this offset anymore

  IF (coef_scatt_ir1%fmv_aer_pha_chn > 0) THEN
    ALLOCATE(phase_channels(coef_scatt_ir1%fmv_aer_pha_chn))

    ! Determine the solar channel numbers
    IF (coef_scatt_ir1%fmv_aer_pha_ioff > 0) THEN
      ! This maintains backward compatibility
      phase_channels(:) = (/ (i, i = coef_scatt_ir1%fmv_aer_pha_ioff, &
                                     coef_scatt_ir1%fmv_aer_pha_ioff + coef_scatt_ir1%fmv_aer_pha_chn - 1) /)
    ELSE
      phase_channels(:) = coef_scatt_ir1%aer_pha_chanlist(:)
    ENDIF

    ! Determine solar channels/phase functions to be extracted
    ALLOCATE(phase_channels_ext(coef_scatt_ir1%fmv_aer_pha_chn))
    ALLOCATE(phase_ext_index(coef_scatt_ir1%fmv_aer_pha_chn))
    coef_scatt_ir2%fmv_aer_pha_chn = 0 ! count up number of solar chans to extract as we go
    DO i = 1, coef_scatt_ir2%fmv_aer_chn
      j = 1  ! here j is the index into the phase_channels array (list of solar channels in file)
      DO WHILE (j <= coef_scatt_ir1%fmv_aer_pha_chn)
        IF (channels(i) == phase_channels(j)) THEN
          coef_scatt_ir2%fmv_aer_pha_chn = coef_scatt_ir2%fmv_aer_pha_chn + 1
          phase_channels_ext(coef_scatt_ir2%fmv_aer_pha_chn) = i
          phase_ext_index(coef_scatt_ir2%fmv_aer_pha_chn) = j
          EXIT
        ENDIF
        j = j + 1
      ENDDO
    ENDDO

    IF (coef_scatt_ir2%fmv_aer_pha_chn > 0) THEN
      ! Copy solar channels to extract into correctly-sized array
      ALLOCATE(coef_scatt_ir2%aer_pha_chanlist(coef_scatt_ir2%fmv_aer_pha_chn))
      coef_scatt_ir2%aer_pha_chanlist(:) = phase_channels_ext(1:coef_scatt_ir2%fmv_aer_pha_chn)

      ! Create map from extracted channel list into pha array
      ALLOCATE(coef_scatt_ir2%aer_pha_index(coef_scatt_ir2%fmv_aer_chn))
      coef_scatt_ir2%aer_pha_index(:) = 0
      coef_scatt_ir2%aer_pha_index(coef_scatt_ir2%aer_pha_chanlist(1:coef_scatt_ir2%fmv_aer_pha_chn)) = &
            (/ (i, i = 1, coef_scatt_ir2%fmv_aer_pha_chn) /)
    ENDIF
    DEALLOCATE(phase_channels_ext)

    ! See rttov_read_ascii_scaercoef.F90 for a description of what the various arrays contain.

    IF (ALLOCATED(phase_channels)) DEALLOCATE(phase_channels)
  ELSE
    ! What do we need to do if there are no phase fns at all?
    coef_scatt_ir2%fmv_aer_pha_chn = 0
  ENDIF

  ALLOCATE(coef_scatt_ir2%fmv_aer_ph_val(coef_scatt_ir2%fmv_aer_ph), stat=err)
  THROWM(err.ne.0, 'allocation of coef_scatt_ir2%fmv_aer_ph_val')
  coef_scatt_ir2%fmv_aer_ph_val = coef_scatt_ir1%fmv_aer_ph_val

  IF (coef_scatt_ir2%fmv_aer_pha_chn > 0) THEN

    ALLOCATE(coef_scatt_ir2%fmv_aer_ph_val_cos(coef_scatt_ir2%fmv_aer_ph), stat=err)
    THROWM(err.ne.0, 'allocation of coef_scatt_ir2%fmv_aer_ph_val_cos')
    coef_scatt_ir2%fmv_aer_ph_val_cos = coef_scatt_ir1%fmv_aer_ph_val_cos

    ALLOCATE(coef_scatt_ir2%ifmv_aer_ph_val(SIZE(coef_scatt_ir1%ifmv_aer_ph_val)), stat=err)
    THROWM(err.ne.0, 'allocation of coef_scatt_ir2%ifmv_aer_ph_val')
    coef_scatt_ir2%ifmv_aer_ph_val = coef_scatt_ir1%ifmv_aer_ph_val

  ENDIF

  ALLOCATE(coef_scatt_ir2%fmv_aer_comp_name(coef_scatt_ir2%fmv_aer_comp), stat=err)
  THROWM(err.ne.0, 'allocation of coef_scatt_ir2%fmv_aer_comp_name')
  coef_scatt_ir2%fmv_aer_comp_name = coef_scatt_ir1%fmv_aer_comp_name

  ALLOCATE(coef_scatt_ir2%fmv_aer_rh(coef_scatt_ir2%fmv_aer_comp), stat=err)
  THROWM(err.ne.0, 'allocation of coef_scatt_ir2%fmv_aer_rh')
  coef_scatt_ir2%fmv_aer_rh = coef_scatt_ir1%fmv_aer_rh

  ALLOCATE(optp2%optpaer(coef_scatt_ir2%fmv_aer_comp), stat=err)
  THROWM(err.ne.0, 'allocation of optpaer')

  DO n = 1, coef_scatt_ir2%fmv_aer_comp

    CALL rttov_nullify_coef_scatt_ir(optp2%optpaer(n))
    ALLOCATE(optp2%optpaer(n)%fmv_aer_rh_val(coef_scatt_ir2%fmv_aer_rh(n)), stat=err)
    THROWM(err.ne.0, 'allocation of optpaer(n)%fmv_aer_rh_val')
    optp2%optpaer(n)%fmv_aer_rh_val = optp1%optpaer(n)%fmv_aer_rh_val

    ALLOCATE(optp2%optpaer(n)%abs(coef_scatt_ir2%fmv_aer_chn,coef_scatt_ir2%fmv_aer_rh(n)), stat=err)
    THROWM(err.ne.0, 'allocation of optpaer(n)%abs')
    optp2%optpaer(n)%abs(:,:) = optp1%optpaer(n)%abs(channels,:)

    ALLOCATE(optp2%optpaer(n)%sca(coef_scatt_ir2%fmv_aer_chn,coef_scatt_ir2%fmv_aer_rh(n)), stat=err)
    THROWM(err.ne.0, 'allocation of optpaer(n)%sca')
    optp2%optpaer(n)%sca(:,:) = optp1%optpaer(n)%sca(channels,:)

    ALLOCATE(optp2%optpaer(n)%bpr(coef_scatt_ir2%fmv_aer_chn,coef_scatt_ir2%fmv_aer_rh(n)), stat=err)
    THROWM(err.ne.0, 'allocation of optpaer(n)%bpr')
    optp2%optpaer(n)%bpr(:,:) = optp1%optpaer(n)%bpr(channels,:)

    IF (coef_scatt_ir2%fmv_aer_pha_chn > 0) THEN

      ALLOCATE(optp2%optpaer(n)%pha(coef_scatt_ir2%fmv_aer_pha_chn,coef_scatt_ir2%fmv_aer_rh(n), &
                                    coef_scatt_ir2%fmv_aer_ph), stat=err)
      THROWM(err.ne.0, 'allocation of optpaer(n)%pha')
      optp2%optpaer(n)%pha(:,:,:) = optp1%optpaer(n)%pha(phase_ext_index(1:coef_scatt_ir2%fmv_aer_pha_chn),:,:)

    ENDIF

  ENDDO

  IF (ALLOCATED(phase_ext_index)) DEALLOCATE (phase_ext_index)

CATCH
END SUBROUTINE rttov_channel_extract_scaercoef