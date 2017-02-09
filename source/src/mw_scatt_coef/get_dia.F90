subroutine get_dia ( dia_froz, i_type, ll_liu_dda, liu_habit)

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
!    Sets up the diameter range for calculations. Not appropriate in
!    the totalice case, where predict_psd_f07 provides it instead.
!
!    OUT: dia_froz   - particle diameter or maximum dimension [cm]
!    IN:  i_type     - hydrometeor type (see mod_mie.F90)   
!         ll_liu_dda - use Liu (2008) DDA results for snow
!         liu_habit  - habit for above, as defined in Liu (2008)


! Current Code Owner: SAF NWP

! History:
! Version   Date        Comment
! -------   ----        -------
!           02/03/2010  New function (Alan Geer)
!           15/03/2013  Fully-flexible Liu shapes (Alan Geer)

use parkind1, only: jprb, jpim, jplm
use mod_mie,  only: n_dia
!INTF_OFF
use mod_mie,  only: d_min, d_max, d_liu_max, d_liu_min
!INTF_ON
implicit none

! Interface
real (kind=jprb),    intent (  out) :: dia_froz(n_dia)
integer (kind=jpim), intent (in   ) :: i_type, liu_habit
logical (kind=jplm), intent (in   ) :: ll_liu_dda
!INTF_END

! Local variables
integer (kind=jpim) :: i_dia
real (kind=jprb)    :: d_dia_froz, d_max_use, d_min_use

!* size ranges (computed later in the totalice case)
if (ll_liu_dda) then
  ! Apply size limits appropriate to Liu shapes
  d_max_use = d_liu_max(liu_habit) 
  d_min_use = d_liu_min(liu_habit)
else
  d_max_use = d_max (i_type)
  d_min_use = d_min (i_type)
endif
d_dia_froz = (d_max_use - d_min_use) / n_dia
do i_dia = 1, n_dia
  dia_froz (i_dia) = d_min_use + (i_dia - 1) * d_dia_froz
enddo

return
end subroutine get_dia
