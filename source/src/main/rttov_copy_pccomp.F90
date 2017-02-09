SUBROUTINE rttov_copy_pccomp(pccomp1, pccomp2)
! Description:
!   Copy PC components structure
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
! Method:
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
  USE rttov_types, ONLY : rttov_pccomp
  TYPE(rttov_pccomp), INTENT(INOUT) :: pccomp1
  TYPE(rttov_pccomp), INTENT(IN)    :: pccomp2
!INTF_END
    IF (Associated(pccomp1%pcscores) &
    .AND. Associated(pccomp2%pcscores)) THEN
      pccomp1%pcscores = pccomp2%pcscores
    ENDIF
    IF (Associated(pccomp1%bt_pccomp) &
    .AND. Associated(pccomp2%bt_pccomp)) THEN
      pccomp1%bt_pccomp = pccomp2%bt_pccomp
    ENDIF
    IF (Associated(pccomp1%total_pccomp) &
    .AND. Associated(pccomp2%total_pccomp)) THEN
      pccomp1%total_pccomp = pccomp2%total_pccomp
    ENDIF
END SUBROUTINE 
