SUBROUTINE rttov_init_aux_prof(aux_prof)
! Description:
!   Initialise auxiliary profile structure
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
!INTF_OFF
  USE parkind1, ONLY : jprb, jpim
!INTF_ON
  TYPE(profile_aux), INTENT(INOUT) :: aux_prof
!INTF_END
  aux_prof%s%nearestlev_surf = 0_jpim
  aux_prof%s%pfraction_surf  = 0._jprb
  aux_prof%s%nearestlev_ctp  = 0_jpim
  aux_prof%s%pfraction_ctp   = 0._jprb
  aux_prof%s%cfraction       = 0._jprb
  IF (associated(aux_prof%debye_prof)) aux_prof%debye_prof = 0._jprb
  IF (associated(aux_prof%dg)        ) aux_prof%dg         = 0._jprb
  IF (associated(aux_prof%fac1_dg)   ) aux_prof%fac1_dg    = 0._jprb
  IF (associated(aux_prof%fac2_dg)   ) aux_prof%fac2_dg    = 0._jprb
  IF (associated(aux_prof%fac3_dg)   ) aux_prof%fac3_dg    = 0._jprb
  IF (associated(aux_prof%iaernum)   ) aux_prof%iaernum    = 0_jpim
  IF (associated(aux_prof%iaertyp)   ) aux_prof%iaertyp    = 0_jpim
  IF (associated(aux_prof%relhum)    ) aux_prof%relhum     = 0._jprb
  IF (associated(aux_prof%relhumref) ) aux_prof%relhumref  = 0._jprb
END SUBROUTINE 
