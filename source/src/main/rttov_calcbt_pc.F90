!
SUBROUTINE rttov_calcbt_pc(chanprof_in, coef_pccomp, pccomp)
! Description:
! To convert an array of radiances in many channels
! to equivalent black body brightness temperatures.
! Planck function is applied with a "central wave number"
! Temperature is corrected by "band corrections" coefficients
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
!    Copyright 2002, EUMETSAT, All Rights Reserved.
!
! Method: Uses band correction coefficients for each channel
!         read in from RT coefficient file.
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.0    02/12/2009  Marco Matricardi:Routine to be used in the radiance reconstruction
!                    when using principal components
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
! Declarations:
! Modules used:
! Imported Type Definitions:
  USE rttov_types, ONLY : rttov_chanprof, rttov_coef_pccomp, rttov_pccomp
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jprb
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof_in(:)
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp   ! Coefficients
  TYPE(rttov_pccomp     ), INTENT(INOUT) :: pccomp        ! input radiances and output BT
!INTF_END
! radiances are expressed in mw/cm-1/ster/sq.m
! and temperatures in Kelvin
!local variables:
  REAL   (KIND=jprb) :: tstore2
  INTEGER(KIND=jpim) :: chan, i
  INTEGER(KIND=jpim) :: nchannels_rec
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
  nchannels_rec = size(chanprof_in)
!- End of header ------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCBT_PC', 0_jpim, ZHOOK_HANDLE)
  DO i = 1, nchannels_rec
    chan = chanprof_in(i)%chan
    tstore2             = coef_pccomp%planck2_in(chan) / Log(1 + coef_pccomp%planck1_in(chan) / pccomp%total_pccomp(i))
    pccomp%bt_pccomp(i) = (tstore2 - coef_pccomp%ff_bco_in(chan)) / coef_pccomp%ff_bcs_in(chan)
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCBT_PC', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_calcbt_pc
