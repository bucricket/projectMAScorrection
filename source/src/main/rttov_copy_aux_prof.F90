SUBROUTINE rttov_copy_aux_prof( &
            & aux_prof1, &
            & aux_prof2)
! Description:
!   Copy auxiliary profile structure
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
  USE rttov_types, ONLY : profile_aux
  TYPE(profile_aux), INTENT(INOUT)             :: aux_prof1
  TYPE(profile_aux), INTENT(IN)                :: aux_prof2
!INTF_END
  aux_prof1%s = aux_prof2%s
  IF (Associated(aux_prof1%debye_prof) .AND. Associated(aux_prof2%debye_prof))     &
    &  aux_prof1%debye_prof = aux_prof2%debye_prof
  IF (Associated(aux_prof1%dg) .AND. Associated(aux_prof2%dg)                ) aux_prof1%dg         = aux_prof2%dg
  IF (Associated(aux_prof1%fac1_dg) .AND. Associated(aux_prof2%fac1_dg)      ) aux_prof1%fac1_dg    = aux_prof2%fac1_dg
  IF (Associated(aux_prof1%fac2_dg) .AND. Associated(aux_prof2%fac2_dg)      ) aux_prof1%fac2_dg    = aux_prof2%fac2_dg
  IF (Associated(aux_prof1%fac3_dg) .AND. Associated(aux_prof2%fac3_dg)      ) aux_prof1%fac3_dg    = aux_prof2%fac3_dg
  IF (Associated(aux_prof1%iaernum) .AND. Associated(aux_prof2%iaernum)      ) aux_prof1%iaernum    = aux_prof2%iaernum
  IF (Associated(aux_prof1%iaertyp) .AND. Associated(aux_prof2%iaertyp)      ) aux_prof1%iaertyp    = aux_prof2%iaertyp
  IF (Associated(aux_prof1%relhum) .AND. Associated(aux_prof2%relhum)        ) aux_prof1%relhum     = aux_prof2%relhum
  IF (Associated(aux_prof1%relhumref) .AND. Associated(aux_prof2%relhumref)  ) aux_prof1%relhumref  = aux_prof2%relhumref
END SUBROUTINE 
