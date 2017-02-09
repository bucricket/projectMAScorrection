subroutine mie_one_temp ( temp, wavelength, f_ghz, i_type, dens, &
  & psd, regime, ll_melt, ll_dsd, liu_habit, is_loaded, &
  & ext_tab, ssa_tab, asm_tab, zef_tab)

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
! Computes scattering parameters for a particular frequency, 
! temperature, and hydrometeor type, as a function of lwp
!
!    IN: temp       - temperature       [K]
!        wavelength -                   [cm]
!        f_ghz      - frequency         [GHz]
!        i_type     - hydrometeor type (see mod_mie.F90)   
!        ll_melt    - do melting layer of Bauer (2001) 
!        dens       - density parametrization 
!        psd        - PSD parametrization 
!        regime     - Field et al. (2007) PSD regime
!        ll_melt
!        ll_dsd     - Panegrossi et al. n0 vs T
!        liu_habit  - habit for Liu (2008) DDA shapes
!        is_loaded  - used by Liu DDA tables
!                     development work. Too slow for operational use.
!
!    OUT: ext_tab   - Bulk extinction               [km^-1]
!         ssa_tab   - Bulk single scattering albedo [ ]
!         asm_tab   - Bulk asymmetry paramter       [ ]
!         zef_tab   - Bulk radar reflectivity ?     [ ?? ] 
!

! Current Code Owner: SAF NWP

! History:
! Version   Date        Comment
! -------   ----        -------
!           10/03/2010  Basic routine for parallelisation (Alan Geer)
!           13/03/2013  Fully-flexible PSD, density and shape (Alan Geer)

use parkind1, only: jprb, jpim
use mod_mie, only: n_lwc
!INTF_OFF
use mod_mie, only: n_dia, get_lwc, cm2km, i_totalice
!INTF_ON

implicit none

real (kind=jprb),    intent(in   ) :: temp, wavelength, f_ghz
integer (kind=jpim), intent(in   ) :: i_type, liu_habit, dens, psd
logical,             intent(in   ) :: ll_melt, ll_dsd
character,           intent(in   ) :: regime
integer (kind=jpim), intent(inout) :: is_loaded
real (kind=jprb),    intent(  out) :: ext_tab(n_lwc), ssa_tab(n_lwc), asm_tab(n_lwc), zef_tab(n_lwc)
!INTF_END

logical :: ll_liu_dda, ll_variable_d
integer (kind=jpim) :: i_lwc!, i_habit
real (kind=jprb), dimension (n_dia) :: q_ext, q_sct, q_asm, q_bsct, nd, dia_froz
real (kind=jprb)                    :: ext, sct, asm, bsct, lwc

!* Interface blocks
#include "get_dia.interface"
#include "scattering.interface"
#include "set_spectra.interface"
#include "mie_one_wc.interface"

ll_liu_dda    = liu_habit >= 0
ll_variable_d = i_type == i_totalice .and. .not. ll_liu_dda

ext_tab = 0.0
ssa_tab = 0.0
asm_tab = 0.0
zef_tab = 0.0

if( .not. ll_variable_d ) then  

  ! Diameter range
  call get_dia( dia_froz, i_type, ll_liu_dda, liu_habit)

  ! Compute scattering parameters as a function of diameter. 
  call scattering( dia_froz, temp, wavelength, f_ghz, i_type, dens, &
    & ll_liu_dda, liu_habit, is_loaded, q_ext, q_sct, q_asm, q_bsct)

endif

do i_lwc = 1, n_lwc

  lwc = get_lwc(i_lwc)

  ! Get size distribution for this LWC 
  call set_spectra (i_type, lwc, temp, f_ghz, dia_froz, nd, psd, ll_dsd, dens, regime, &
    & ll_liu_dda, liu_habit, ll_variable_d)

  if( ll_variable_d ) then  

    ! Totalice diameters change according to the lwc, so must recalculate Mie each time - much slower...
    call scattering( dia_froz, temp, wavelength, f_ghz, i_type, dens, &
      & ll_liu_dda, liu_habit, is_loaded, q_ext, q_sct, q_asm, q_bsct)

  endif

  ! Integrate over sizes to get bulk layer properties
  call mie_one_wc( wavelength, f_ghz, i_type, (ll_melt .and. nint(temp) == 273), dia_froz, nd, &
    & q_ext, q_sct, q_asm, q_bsct, ext, sct, asm, bsct)

  ssa_tab(i_lwc) = ssa_tab(i_lwc) + sct / ext 
  ext_tab(i_lwc) = ext_tab(i_lwc) + ext * cm2km
  asm_tab(i_lwc) = asm_tab(i_lwc) + asm
  zef_tab(i_lwc) = zef_tab(i_lwc) + bsct

end do

end subroutine mie_one_temp
