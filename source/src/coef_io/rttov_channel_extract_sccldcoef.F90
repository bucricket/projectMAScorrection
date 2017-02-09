SUBROUTINE rttov_channel_extract_sccldcoef( &
              err,            &
              coef_scatt_ir1, &
              optp1,          &
              coef_scatt_ir2, &
              optp2,          &
              channels)
! Description:
!
!   Given a cloud coef structure, extract the data
!   for the channels in the given list to a second
!   uninitialised cloud coef structure.
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

  ! --------------------------------------------------------------------------
  ! Water cloud
  ! --------------------------------------------------------------------------

  ! Scalar variables

  ! A few of these are the same...
  coef_scatt_ir2%fmv_wcl_comp       = coef_scatt_ir1%fmv_wcl_comp
  coef_scatt_ir2%fmv_wcl_ph         = coef_scatt_ir1%fmv_wcl_ph
  coef_scatt_ir2%fmv_wcl_ph_val_min = coef_scatt_ir1%fmv_wcl_ph_val_min

  ! ... but we have to work out the solar channel numbers
  coef_scatt_ir2%fmv_wcl_chn      = SIZE(channels)
  coef_scatt_ir2%fmv_wcl_pha_ioff = 0   ! Don't use this offset anymore

  IF (coef_scatt_ir1%fmv_wcl_pha_chn > 0) THEN
    ALLOCATE(phase_channels(coef_scatt_ir1%fmv_wcl_pha_chn))

    ! Determine the solar channel numbers
    IF (coef_scatt_ir1%fmv_wcl_pha_ioff > 0) THEN
      ! This maintains backward compatibility
      phase_channels(:) = (/ (i, i = coef_scatt_ir1%fmv_wcl_pha_ioff, &
                                     coef_scatt_ir1%fmv_wcl_pha_ioff + coef_scatt_ir1%fmv_wcl_pha_chn - 1) /)
    ELSE
      phase_channels(:) = coef_scatt_ir1%wcl_pha_chanlist(:)
    ENDIF

    ! Determine solar channels/phase functions to be extracted
    ALLOCATE(phase_channels_ext(coef_scatt_ir1%fmv_wcl_pha_chn))
    ALLOCATE(phase_ext_index(coef_scatt_ir1%fmv_wcl_pha_chn))
    coef_scatt_ir2%fmv_wcl_pha_chn = 0 ! count up number of solar chans to extract as we go
    DO i = 1, coef_scatt_ir2%fmv_wcl_chn
      j = 1  ! here j is the index into the phase_channels array (list of solar channels in file)
      DO WHILE (j <= coef_scatt_ir1%fmv_wcl_pha_chn)
        IF (channels(i) == phase_channels(j)) THEN
          coef_scatt_ir2%fmv_wcl_pha_chn = coef_scatt_ir2%fmv_wcl_pha_chn + 1
          phase_channels_ext(coef_scatt_ir2%fmv_wcl_pha_chn) = i
          phase_ext_index(coef_scatt_ir2%fmv_wcl_pha_chn) = j
          EXIT
        ENDIF
        j = j + 1
      ENDDO
    ENDDO

    IF (coef_scatt_ir2%fmv_wcl_pha_chn > 0) THEN
      ! Copy solar channels to extract into correctly-sized array
      ALLOCATE(coef_scatt_ir2%wcl_pha_chanlist(coef_scatt_ir2%fmv_wcl_pha_chn))
      coef_scatt_ir2%wcl_pha_chanlist(:) = phase_channels_ext(1:coef_scatt_ir2%fmv_wcl_pha_chn)

      ! Create map from extracted channel list into pha array
      ALLOCATE(coef_scatt_ir2%wcl_pha_index(coef_scatt_ir2%fmv_wcl_chn))
      coef_scatt_ir2%wcl_pha_index(:) = 0
      coef_scatt_ir2%wcl_pha_index(coef_scatt_ir2%wcl_pha_chanlist(1:coef_scatt_ir2%fmv_wcl_pha_chn)) = &
            (/ (i, i = 1, coef_scatt_ir2%fmv_wcl_pha_chn) /)
    ENDIF
    DEALLOCATE(phase_channels_ext)

    ! See rttov_read_ascii_sccldcoef.F90 for a description of what the various arrays contain.

    IF (ALLOCATED(phase_channels)) DEALLOCATE(phase_channels)
  ELSE
    ! What do we need to do if there are no phase fns at all?
    coef_scatt_ir2%fmv_wcl_pha_chn = 0
  ENDIF

  ALLOCATE(coef_scatt_ir2%fmv_wcl_ph_val(coef_scatt_ir2%fmv_wcl_ph), stat=err)
  THROWM(err.ne.0, 'allocation of coef_scatt_ir2%fmv_wcl_ph_val')
  coef_scatt_ir2%fmv_wcl_ph_val = coef_scatt_ir1%fmv_wcl_ph_val

  IF (coef_scatt_ir2%fmv_wcl_pha_chn > 0) THEN

    ALLOCATE(coef_scatt_ir2%fmv_wcl_ph_val_cos(coef_scatt_ir2%fmv_wcl_ph), stat=err)
    THROWM(err.ne.0, 'allocation of coef_scatt_ir2%fmv_wcl_ph_val_cos')
    coef_scatt_ir2%fmv_wcl_ph_val_cos = coef_scatt_ir1%fmv_wcl_ph_val_cos

    ALLOCATE(coef_scatt_ir2%ifmv_wcl_ph_val(SIZE(coef_scatt_ir1%ifmv_wcl_ph_val)), stat=err)
    THROWM(err.ne.0, 'allocation of coef_scatt_ir2%ifmv_wcl_ph_val')
    coef_scatt_ir2%ifmv_wcl_ph_val = coef_scatt_ir1%ifmv_wcl_ph_val

  ENDIF

  ALLOCATE(coef_scatt_ir2%fmv_wcl_comp_name(coef_scatt_ir2%fmv_wcl_comp), stat=err)
  THROWM(err.ne.0, 'allocation of coef_scatt_ir2%fmv_wcl_comp_name')
  coef_scatt_ir2%fmv_wcl_comp_name = coef_scatt_ir1%fmv_wcl_comp_name

  ALLOCATE(coef_scatt_ir2%fmv_wcl_rh(coef_scatt_ir2%fmv_wcl_comp), stat=err)
  THROWM(err.ne.0, 'allocation of coef_scatt_ir2%fmv_wcl_rh')
  coef_scatt_ir2%fmv_wcl_rh = coef_scatt_ir1%fmv_wcl_rh

  ALLOCATE (coef_scatt_ir2%confac(coef_scatt_ir2%fmv_wcl_comp), stat=err)
  THROWM(err.ne.0, 'allocation of coef_scatt_ir2%confac')
  coef_scatt_ir2%confac = coef_scatt_ir1%confac

  ALLOCATE(optp2%optpwcl(coef_scatt_ir2%fmv_wcl_comp), stat=err)
  THROWM(err.ne.0, 'allocation of optpwcl')

  DO n = 1, coef_scatt_ir2%fmv_wcl_comp

    CALL rttov_nullify_coef_scatt_ir(optp2%optpwcl(n))
    ALLOCATE(optp2%optpwcl(n)%fmv_wcl_rh_val(coef_scatt_ir2%fmv_wcl_rh(n)), stat=err)
    THROWM(err.ne.0, 'allocation of optpwcl(n)%fmv_wcl_rh_val')
    optp2%optpwcl(n)%fmv_wcl_rh_val = optp1%optpwcl(n)%fmv_wcl_rh_val

    ALLOCATE(optp2%optpwcl(n)%abs(coef_scatt_ir2%fmv_wcl_chn,coef_scatt_ir2%fmv_wcl_rh(n)), stat=err)
    THROWM(err.ne.0, 'allocation of optpwcl(n)%abs')
    optp2%optpwcl(n)%abs(:,:) = optp1%optpwcl(n)%abs(channels,:)

    ALLOCATE(optp2%optpwcl(n)%sca(coef_scatt_ir2%fmv_wcl_chn,coef_scatt_ir2%fmv_wcl_rh(n)), stat=err)
    THROWM(err.ne.0, 'allocation of optpwcl(n)%sca')
    optp2%optpwcl(n)%sca(:,:) = optp1%optpwcl(n)%sca(channels,:)

    ALLOCATE(optp2%optpwcl(n)%bpr(coef_scatt_ir2%fmv_wcl_chn,coef_scatt_ir2%fmv_wcl_rh(n)), stat=err)
    THROWM(err.ne.0, 'allocation of optpwcl(n)%bpr')
    optp2%optpwcl(n)%bpr(:,:) = optp1%optpwcl(n)%bpr(channels,:)

    IF (coef_scatt_ir2%fmv_wcl_pha_chn > 0) THEN

      ALLOCATE(optp2%optpwcl(n)%pha(coef_scatt_ir2%fmv_wcl_pha_chn,coef_scatt_ir2%fmv_wcl_rh(n), &
                                    coef_scatt_ir2%fmv_wcl_ph), stat=err)
      THROWM(err.ne.0, 'allocation of optpwcl(n)%pha')
      optp2%optpwcl(n)%pha(:,:,:) = optp1%optpwcl(n)%pha(phase_ext_index(1:coef_scatt_ir2%fmv_wcl_pha_chn),:,:)

    ENDIF

  ENDDO

  IF (ALLOCATED(phase_ext_index)) DEALLOCATE (phase_ext_index)


  ! --------------------------------------------------------------------------
  ! Ice cloud
  ! --------------------------------------------------------------------------

  ! Scalar variables

  ! A few of these are the same...
  coef_scatt_ir2%fmv_icl_comp       = coef_scatt_ir1%fmv_icl_comp
  coef_scatt_ir2%fmv_icl_ishp       = coef_scatt_ir1%fmv_icl_ishp
  coef_scatt_ir2%icl_nabs           = coef_scatt_ir1%icl_nabs
  coef_scatt_ir2%icl_nsca           = coef_scatt_ir1%icl_nsca
  coef_scatt_ir2%icl_nbpr           = coef_scatt_ir1%icl_nbpr
  coef_scatt_ir2%fmv_icl_ph         = coef_scatt_ir1%fmv_icl_ph
  coef_scatt_ir2%fmv_icl_ph_val_min = coef_scatt_ir1%fmv_icl_ph_val_min

  ! ... but we have to work out the solar channel numbers
  coef_scatt_ir2%fmv_icl_chn      = SIZE(channels)
  coef_scatt_ir2%fmv_icl_pha_ioff = 0   ! Don't use this offset anymore

  IF (coef_scatt_ir1%fmv_icl_pha_chn > 0) THEN
    ALLOCATE(phase_channels(coef_scatt_ir1%fmv_icl_pha_chn))

    ! Determine the solar channel numbers
    IF (coef_scatt_ir1%fmv_icl_pha_ioff > 0) THEN
      ! This maintains backward compatibility
      phase_channels(:) = (/ (i, i = coef_scatt_ir1%fmv_icl_pha_ioff, &
                                     coef_scatt_ir1%fmv_icl_pha_ioff + coef_scatt_ir1%fmv_icl_pha_chn - 1) /)
    ELSE
      phase_channels(:) = coef_scatt_ir1%icl_pha_chanlist(:)
    ENDIF

    ! Determine solar channels/phase functions to be extracted
    ALLOCATE(phase_channels_ext(coef_scatt_ir1%fmv_icl_pha_chn))
    ALLOCATE(phase_ext_index(coef_scatt_ir1%fmv_icl_pha_chn))
    coef_scatt_ir2%fmv_icl_pha_chn = 0 ! count up number of solar chans to extract as we go
    DO i = 1, coef_scatt_ir2%fmv_icl_chn
      j = 1  ! here j is the index into the phase_channels array (list of solar channels in file)
      DO WHILE (j <= coef_scatt_ir1%fmv_icl_pha_chn)
        IF (channels(i) == phase_channels(j)) THEN
          coef_scatt_ir2%fmv_icl_pha_chn = coef_scatt_ir2%fmv_icl_pha_chn + 1
          phase_channels_ext(coef_scatt_ir2%fmv_icl_pha_chn) = i
          phase_ext_index(coef_scatt_ir2%fmv_icl_pha_chn) = j
          EXIT
        ENDIF
        j = j + 1
      ENDDO
    ENDDO

    IF (coef_scatt_ir2%fmv_icl_pha_chn > 0) THEN
      ! Copy solar channels to extract into correctly-sized array
      ALLOCATE(coef_scatt_ir2%icl_pha_chanlist(coef_scatt_ir2%fmv_icl_pha_chn))
      coef_scatt_ir2%icl_pha_chanlist(:) = phase_channels_ext(1:coef_scatt_ir2%fmv_icl_pha_chn)

      ! Create map from extracted channel list into pha array
      ALLOCATE(coef_scatt_ir2%icl_pha_index(coef_scatt_ir2%fmv_icl_chn))
      coef_scatt_ir2%icl_pha_index(:) = 0
      coef_scatt_ir2%icl_pha_index(coef_scatt_ir2%icl_pha_chanlist(1:coef_scatt_ir2%fmv_icl_pha_chn)) = &
            (/ (i, i = 1, coef_scatt_ir2%fmv_icl_pha_chn) /)
    ENDIF
    DEALLOCATE(phase_channels_ext)

    ! See rttov_read_ascii_sccldcoef.F90 for a description of what the various arrays contain.

    IF (ALLOCATED(phase_channels)) DEALLOCATE(phase_channels)
  ELSE
    ! What do we need to do if there are no phase fns at all?
    coef_scatt_ir2%fmv_icl_pha_chn = 0
  ENDIF

  ALLOCATE(coef_scatt_ir2%fmv_icl_ph_val(coef_scatt_ir2%fmv_icl_ph), stat=err)
  THROWM(err.ne.0, 'allocation of coef_scatt_ir2%fmv_icl_ph_val')
  coef_scatt_ir2%fmv_icl_ph_val = coef_scatt_ir1%fmv_icl_ph_val

  IF (coef_scatt_ir2%fmv_icl_pha_chn > 0) THEN

    ALLOCATE(coef_scatt_ir2%fmv_icl_ph_val_cos(coef_scatt_ir2%fmv_icl_ph), stat=err)
    THROWM(err.ne.0, 'allocation of coef_scatt_ir2%fmv_icl_ph_val_cos')
    coef_scatt_ir2%fmv_icl_ph_val_cos = coef_scatt_ir1%fmv_icl_ph_val_cos

    ALLOCATE(coef_scatt_ir2%ifmv_icl_ph_val(SIZE(coef_scatt_ir1%ifmv_icl_ph_val)), stat=err)
    THROWM(err.ne.0, 'allocation of coef_scatt_ir2%ifmv_icl_ph_val')
    coef_scatt_ir2%ifmv_icl_ph_val = coef_scatt_ir1%ifmv_icl_ph_val

  ENDIF

  ALLOCATE(coef_scatt_ir2%fmv_icl_comp_name(coef_scatt_ir2%fmv_icl_comp,2), stat=err)
  THROWM(err.ne.0, 'allocation of coef_scatt_ir2%fmv_icl_comp_name')
  coef_scatt_ir2%fmv_icl_comp_name = coef_scatt_ir1%fmv_icl_comp_name

  ALLOCATE(coef_scatt_ir2%fmv_icl_dg(coef_scatt_ir2%fmv_icl_comp, coef_scatt_ir2%fmv_icl_ishp), stat=err)
  THROWM(err.ne.0, 'allocation of coef_scatt_ir2%fmv_icl_dg')
  coef_scatt_ir2%fmv_icl_dg = coef_scatt_ir1%fmv_icl_dg

  ALLOCATE (coef_scatt_ir2%fmv_icl_comp_name(coef_scatt_ir2%fmv_icl_comp,2), stat=err)
  THROWM(err.ne.0, 'allocation of coef_scatt_ir2%fmv_icl_comp_name')
  coef_scatt_ir2%fmv_icl_comp_name = coef_scatt_ir1%fmv_icl_comp_name

  ALLOCATE(optp2%optpicl(coef_scatt_ir2%fmv_icl_ishp), stat=err)
  THROWM(err.ne.0, 'allocation of optpicl')

  DO n = 1, coef_scatt_ir2%fmv_icl_ishp

    CALL rttov_nullify_coef_scatt_ir(optp2%optpicl(n))

    ALLOCATE(optp2%optpicl(n)%abs(coef_scatt_ir2%fmv_icl_chn,coef_scatt_ir2%icl_nabs), stat=err)
    THROWM(err.ne.0, 'allocation of optpicl(n)%abs')
    optp2%optpicl(n)%abs(:,:) = optp1%optpicl(n)%abs(channels,:)

    ALLOCATE(optp2%optpicl(n)%sca(coef_scatt_ir2%fmv_icl_chn,coef_scatt_ir2%icl_nsca), stat=err)
    THROWM(err.ne.0, 'allocation of optpicl(n)%sca')
    optp2%optpicl(n)%sca(:,:) = optp1%optpicl(n)%sca(channels,:)

    ALLOCATE(optp2%optpicl(n)%bpr(coef_scatt_ir2%fmv_icl_chn,coef_scatt_ir2%icl_nbpr), stat=err)
    THROWM(err.ne.0, 'allocation of optpicl(n)%bpr')
    optp2%optpicl(n)%bpr(:,:) = optp1%optpicl(n)%bpr(channels,:)

    IF (coef_scatt_ir2%fmv_icl_pha_chn > 0) THEN

      ALLOCATE(optp2%optpicl(n)%pha(coef_scatt_ir2%fmv_icl_pha_chn,coef_scatt_ir2%fmv_icl_comp, &
                                    coef_scatt_ir2%fmv_icl_ph), stat=err)
      THROWM(err.ne.0, 'allocation of optpicl(n)%pha')
      optp2%optpicl(n)%pha(:,:,:) = optp1%optpicl(n)%pha(phase_ext_index(1:coef_scatt_ir2%fmv_icl_pha_chn),:,:)

    ENDIF

  ENDDO

  IF (ALLOCATED(phase_ext_index)) DEALLOCATE (phase_ext_index)

CATCH
END SUBROUTINE rttov_channel_extract_sccldcoef