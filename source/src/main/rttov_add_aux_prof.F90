SUBROUTINE rttov_add_aux_prof(aux_prof, aux_prof1, aux_prof2)
! Description:
!   Adds two auxiliary profile structures
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
  TYPE(profile_aux), INTENT(INOUT) :: aux_prof
  TYPE(profile_aux), INTENT(IN)    :: aux_prof1
  TYPE(profile_aux), INTENT(IN)    :: aux_prof2
!INTF_END
  aux_prof%s%pfraction_surf = aux_prof1%s%pfraction_surf + aux_prof2%s%pfraction_surf
  aux_prof%s%pfraction_ctp  = aux_prof1%s%pfraction_ctp + aux_prof2%s%pfraction_ctp
  aux_prof%s%cfraction      = aux_prof1%s%cfraction + aux_prof2%s%cfraction
  IF (associated(aux_prof%debye_prof) .AND. associated(aux_prof1%debye_prof) .AND. associated(aux_prof2%debye_prof))     &
    &  aux_prof%debye_prof = aux_prof1%debye_prof + aux_prof2%debye_prof
  IF (associated(aux_prof%dg) .AND. associated(aux_prof1%dg) .AND. associated(aux_prof2%dg)                        )     &
    &  aux_prof%dg         = aux_prof1%dg + aux_prof2%dg
  IF (associated(aux_prof%fac1_dg) .AND. associated(aux_prof1%fac1_dg) .AND. associated(aux_prof2%fac1_dg)         )     &
    &  aux_prof%fac1_dg    = aux_prof1%fac1_dg + aux_prof2%fac1_dg
  IF (associated(aux_prof%fac2_dg) .AND. associated(aux_prof1%fac2_dg) .AND. associated(aux_prof2%fac2_dg)         )     &
    &  aux_prof%fac2_dg    = aux_prof1%fac2_dg + aux_prof2%fac2_dg
  IF (associated(aux_prof%fac3_dg) .AND. associated(aux_prof1%fac3_dg) .AND. associated(aux_prof2%fac3_dg)         )     &
    &  aux_prof%fac3_dg    = aux_prof1%fac3_dg + aux_prof2%fac3_dg
  IF (associated(aux_prof%relhum) .AND. associated(aux_prof1%relhum) .AND. associated(aux_prof2%relhum)            )     &
    &  aux_prof%relhum     = aux_prof1%relhum + aux_prof2%relhum
  IF (associated(aux_prof%relhumref) .AND. associated(aux_prof1%relhumref) .AND. associated(aux_prof2%relhumref)   )     &
    &  aux_prof%relhumref  = aux_prof1%relhumref + aux_prof2%relhumref
END SUBROUTINE 
