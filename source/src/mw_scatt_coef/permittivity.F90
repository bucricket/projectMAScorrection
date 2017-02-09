function permittivity (i_type, f_ghz, temp, density)

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

! Current Code Owner: SAF NWP

! History:
! Version   Date        Comment
! -------   ----        -------
!           10/03/2010  extracted from create_tables_spectra (Alan Geer)
!           03/03/2011  "ice-factor" density option
!           13/03/2013  Fully-flexible PSD, density and shape (Alan Geer)

use parkind1, only: jprb, jpim
! use mod_mie,  only: i_rain, i_snow, i_graupel, i_aggregate, i_clw, i_ciw, i_totalice, &
!                   & perm_air, density_air, density_ice, density_liq, pi
use mod_mie,  only: i_rain, i_clw, perm_air, density_air, density_ice
implicit none

! Interface
integer (kind=jpim), intent ( in) :: i_type
real    (kind=jprb), intent ( in) :: f_ghz
real    (kind=jprb), intent ( in) :: temp
real    (kind=jprb), intent ( in) :: density

! Function definition
complex (kind=jprb) :: permittivity

! Local variables
real    (kind=jprb) :: perm_re, perm_im
real (kind=jprb)    :: fv_air, fv_ice
complex (kind=jprb) :: perm_froz

! Interfaces of called functions
#include "perm_water.interface"
#include "perm_ice.interface"
#include "mg_ellips.interface"
! #include "density_all.interface"

!* permittivities
if (i_type == i_rain .or. i_type == i_clw) then

  !* Liquid
  call perm_water (f_ghz, temp, perm_re, perm_im)
  permittivity = cmplx (perm_re, perm_im, jprb)

else

  !* Frozen
  call perm_ice (f_ghz, temp, perm_re, perm_im)
  perm_froz = cmplx (perm_re, perm_im, jprb)

  fv_ice = (density - density_air) / (density_ice - density_air)
  fv_air = 1.0_jprb - fv_ice

  call mg_ellips (fv_air, perm_froz, perm_air, permittivity)

end if

return
end function permittivity
