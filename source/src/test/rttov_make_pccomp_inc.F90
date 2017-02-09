SUBROUTINE rttov_make_pccomp_inc( &
            & pccomp_inc,       &
            & opts)
! Description:
!   Computes a sensible reconstructed radiance variation, either
!   in brightness temperature or in radiance, or as pcscore
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
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
!
  USE rttov_types, ONLY : rttov_pccomp, rttov_options
!INTF_OFF
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  IMPLICIT NONE
  TYPE(rttov_pccomp),  INTENT(INOUT) :: pccomp_inc
  TYPE(rttov_options), INTENT(IN)    :: opts
!INTF_END
  INTEGER(KIND=jpim) :: i

  If(opts%rt_ir%pc%addradrec)Then
    If (opts%rt_all%switchrad) Then
      Do i = 1, Size(pccomp_inc%bt_pccomp)
        pccomp_inc%bt_pccomp(i) = 0.01_jprb*(Modulo(i,113_jpim)+1)
      End Do
    Else
      Do i = 1, Size(pccomp_inc%total_pccomp)
        pccomp_inc%total_pccomp(i) = 0.01_jprb*(Modulo(i,113_jpim)+1)
      End Do
    End If
  Else
    Do i = 1, Size(pccomp_inc%pcscores)
      pccomp_inc%pcscores(i) = 0.01_jprb*(Modulo(i,113_jpim)+1)
    End Do
  End If
END SUBROUTINE 
