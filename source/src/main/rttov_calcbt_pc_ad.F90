SUBROUTINE rttov_calcbt_pc_ad( &
            & chanprof_in, &
            & coef_pccomp, &
            & pccomp,      &
            & pccomp_ad)
!
! Description:
! To convert an array of radiances in many channels
! to equivalent black body brightness temperatures.
! Planck function is applied with a "central wave number"
! Temperature is corrected by "band corrections" coefficients
! derivative of inverse Planck function with respect to radiance is
!
!                             C1 * C2 * Nu**4
! B-1'(R,Nu) = --------------------------------------------- dR
!                     (    C1 * Nu**3) (     C1 * Nu**3 )**2
!               R**2 *(1 + ----------) ( Ln( ---------- )
!                     (         R    ) (        R       )
!
! which can be reduced to the following, with
!  C1 = C1 * Nu**3
!  C2 = C2 * Nu
!
!                  C1 * B-1(R,Nu)**2
! B-1'(R,Nu) = ----------------------- dR
!                  C2 * R * (R + C1)
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
! 1.0    02/12/2009  Marco Matricardi: Routine to be used for radiance reconstruction
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
  USE parkind1, ONLY : jprb
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof_in(:)
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  TYPE(rttov_pccomp     ), INTENT(IN)    :: pccomp        ! rad%total rad%clear
  TYPE(rttov_pccomp     ), INTENT(INOUT) :: pccomp_ad
! output rad_ad%total AND clear
! input  rad_ad%bt AND bt_clear
!INTF_END
!local variables:
  REAL   (KIND=jprb) :: tstar
  REAL   (KIND=jprb) :: tstar_ad
  REAL   (KIND=jprb) :: radtotal
  INTEGER(KIND=jpim) :: chan, i
  INTEGER(KIND=jpim) :: nchannels_rec
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
  nchannels_rec = size(chanprof_in)
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCBT_PC_AD', 0_jpim, ZHOOK_HANDLE)
  DO i = 1, nchannels_rec
    chan = chanprof_in(i)%chan
! clear radiance
! T star for direct model
    tstar = coef_pccomp%ff_bco_in(chan) + coef_pccomp%ff_bcs_in(chan) * pccomp%bt_pccomp(i)
    radtotal                  = pccomp%total_pccomp(i)
! AD
    tstar_ad                  = pccomp_ad%bt_pccomp(i) / coef_pccomp%ff_bcs_in(chan)
    pccomp_ad%total_pccomp(i) = coef_pccomp%planck1_in(chan) * tstar ** 2 /      &
      & (coef_pccomp%planck2_in(chan) * radtotal * (radtotal + coef_pccomp%planck1_in(chan))) * tstar_ad
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCBT_PC_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_calcbt_pc_ad
