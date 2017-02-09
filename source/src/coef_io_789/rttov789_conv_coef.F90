! Description:
!> @file
!!   Executable for converting RTTOV v7/v8/v9-compatible 
!!   coefficient files to v10/v11 format.
!
!> @brief
!!   Executable for converting RTTOV v7/v8/v9-compatible
!!   coefficient files to v10/v11 format.
!!
!! @details
!!   Usage:
!!   $ rttov789_conv_coef.exe \-\-coef-in ... \-\-coef-out ...
!!
!!   where \-\-coef-in specifies the input v7/v8/v9-compatible
!!   file and \-\-coef-out specifies the output v10/v11-
!!   compatible file. Input files must be ASCII format.
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
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
Program rttov789_conv_coef

#include "throw.h"

  Use rttov_types, Only :  &
        & rttov_coef

  Use rttov_getoptions


  Use rttov_unix_env, Only: rttov_iargc, rttov_exit

  Use parkind1, Only : jpim, jprb

  Use rttov_const, Only :     &
        & gas_id_mixed       , &
        & gas_id_watervapour , &
        & gas_id_ozone       , &
        & gas_id_wvcont      , &
        & gas_id_co2         , &
        & gas_id_n2o         , &
        & gas_id_co          , &
        & gas_id_ch4

  Implicit None

#include "rttov789_readcoeffs.interface"
#include "rttov_init_coef.interface"
#include "rttov_write_ascii_coef.interface"
#include "rttov_errorreport.interface"

  Type(rttov_coef)   :: coef

  Integer(Kind=jpim) :: file_id
  Integer(Kind=jpim) :: errorstatus
  Character(Len=256) :: f_coef_in = ""
  Character(Len=256) :: f_coef_out = ""
  Integer :: ioerr, err
  Real(Kind=Jprb), Pointer :: X1(:), X2(:,:)
  Integer(Kind=jpim) :: nlevels, n2, igaz
  Character(Len=2) :: dataset


TRY

  !- End of header --------------------------------------------------------



  If( rttov_iargc() .eq. 0 ) Then
    Print *, "Usage: "
    Print *, "       --coef-in ...  --coef-out ... "
    Stop
  EndIf


  Call getoption( "--coef-in",  f_coef_in  )

  Call getoption( "--coef-out",  f_coef_out  )

  If( f_coef_out .eq. "" .and. f_coef_in .ne. "" ) f_coef_out = Trim(f_coef_in)//'.bin'


  ! let the subroutine choose a logical unit for the file
  file_id      = 77

    Open( file_id,      File = f_coef_in, form = 'formatted', &
      Status = 'old', Iostat = ioerr )

  THROWM(ioerr.ne.0,"Cannot open "//Trim(f_coef_in))

  Call rttov789_readcoeffs ( errorstatus, coef, file_id = file_id )

  THROWM(errorstatus.ne.0,"Cannot read "//Trim(f_coef_in))


  coef%nlevels     = coef%nlevels + 1
  coef%fmv_lvl     = coef%nlevels
  coef%id_comp_pc  = 0
  coef%IncZeeman   = .false.
  coef%id_comp_lvl = 10

  nlevels = coef%nlevels

  Allocate(X1(nlevels)); X1(2:nlevels) = coef%ref_prfl_p; X1(1) = 0.
  deallocate(coef%ref_prfl_p); coef%ref_prfl_p => X1

  n2 = size(coef%ref_prfl_t,2)
  Allocate(X2(nlevels,n2)); X2(2:nlevels,:) = coef%ref_prfl_t; X2(1,:) = 0.
  deallocate(coef%ref_prfl_t); coef%ref_prfl_t => X2

  n2 = size(coef%ref_prfl_mr,2)
  Allocate(X2(nlevels,n2)); X2(2:nlevels,:) = coef%ref_prfl_mr; X2(1,:) = 0.
  deallocate(coef%ref_prfl_mr); coef%ref_prfl_mr => X2

  Allocate(X1(nlevels)); X1(2:nlevels) = coef%lim_prfl_p; X1(1) = 0.
  deallocate(coef%lim_prfl_p); coef%lim_prfl_p => X1

  Allocate(X1(nlevels)); X1(2:nlevels) = coef%lim_prfl_tmax; X1(1) = 0.
  deallocate(coef%lim_prfl_tmax); coef%lim_prfl_tmax => X1

  Allocate(X1(nlevels)); X1(2:nlevels) = coef%lim_prfl_tmin; X1(1) = 0.
  deallocate(coef%lim_prfl_tmin); coef%lim_prfl_tmin => X1

  n2 = size(coef%lim_prfl_gmax,2)
  Allocate(X2(nlevels,n2)); X2(2:nlevels,:) = coef%lim_prfl_gmax; X2(1,:) = 0.
  deallocate(coef%lim_prfl_gmax); coef%lim_prfl_gmax => X2

  n2 = size(coef%lim_prfl_gmin,2)
  Allocate(X2(nlevels,n2)); X2(2:nlevels,:) = coef%lim_prfl_gmin; X2(1,:) = 0.
  deallocate(coef%lim_prfl_gmin); coef%lim_prfl_gmin => X2


  Call rttov_init_coef( errorstatus, coef )

  dataset = grok_profile_dataset(coef%line_by_line)

  ! Top level of reference profiles
  Select Case (dataset)
    Case ("52")
      coef%ref_prfl_p(1) = 0.0050_jprb
      coef%ref_prfl_t(1,:) = 0.220285E+03_jprb
      Do igaz = 1, coef%fmv_gas
        Select case(coef%fmv_gas_id(igaz))
        Case (gas_id_mixed       )
          coef%ref_prfl_mr(1,igaz) = -.999900E+04_jprb
        Case (gas_id_watervapour )
          coef%ref_prfl_mr(1,igaz) = 0.160788E+01_jprb
        Case (gas_id_ozone       )
          coef%ref_prfl_mr(1,igaz) = 0.231238E+00_jprb
        Case (gas_id_wvcont      )
          coef%ref_prfl_mr(1,igaz) = 0.160788E+01_jprb
        Case (gas_id_co2         )
          coef%ref_prfl_mr(1,igaz) = 0.354357E+03_jprb
        Case (gas_id_n2o         )
          coef%ref_prfl_mr(1,igaz) = 0.731000E-02_jprb
        Case (gas_id_co          )
          coef%ref_prfl_mr(1,igaz) = 0.310768E+01_jprb
        Case (gas_id_ch4         )
          coef%ref_prfl_mr(1,igaz) = 0.145130E+00_jprb
        End Select
      EndDo

    Case("83")
      coef%ref_prfl_p(1) = 0.500000E-02_jprb
      coef%ref_prfl_t(1,:) = 0.191313E+03_jprb
      Do igaz = 1, coef%fmv_gas
        Select case(coef%fmv_gas_id(igaz))
        Case (gas_id_mixed       )
          coef%ref_prfl_mr(1,igaz) = -.999900E+04_jprb
        Case (gas_id_watervapour )
          coef%ref_prfl_mr(1,igaz) = 0.262772E+01_jprb
        Case (gas_id_ozone       )
          coef%ref_prfl_mr(1,igaz) = 0.296166E+00_jprb
        Case (gas_id_wvcont      )
          coef%ref_prfl_mr(1,igaz) = 0.262772E+01_jprb
        Case (gas_id_co2         )
          coef%ref_prfl_mr(1,igaz) = 0.374361E+03_jprb
        Case (gas_id_n2o         )
          coef%ref_prfl_mr(1,igaz) = 0.814357E-02_jprb
        Case (gas_id_co          )
          coef%ref_prfl_mr(1,igaz) = 0.844725E+00_jprb
        Case (gas_id_ch4         )
          coef%ref_prfl_mr(1,igaz) = 0.150528E+00_jprb
        End Select
      EndDo

    Case("43")
      coef%ref_prfl_p(1) = 0.500000E-02_jprb
      coef%ref_prfl_t(1,:) = coef%ref_prfl_t(2,:)
      Do igaz = 1, coef%fmv_gas
        Select case(coef%fmv_gas_id(igaz))
        Case (gas_id_mixed       )
          coef%ref_prfl_mr(1,igaz) = -.999900E+04_jprb
        Case (gas_id_watervapour )
          coef%ref_prfl_mr(1,igaz) = coef%ref_prfl_mr(2,igaz)
        Case (gas_id_ozone       )
          coef%ref_prfl_mr(1,igaz) = coef%ref_prfl_mr(2,igaz)
        Case (gas_id_wvcont      )
          coef%ref_prfl_mr(1,igaz) = coef%ref_prfl_mr(2,igaz)
        Case (gas_id_co2         )
          coef%ref_prfl_mr(1,igaz) = coef%ref_prfl_mr(2,igaz)
        Case (gas_id_n2o         )
          coef%ref_prfl_mr(1,igaz) = coef%ref_prfl_mr(2,igaz)
        Case (gas_id_co          )
          coef%ref_prfl_mr(1,igaz) = coef%ref_prfl_mr(2,igaz)
        Case (gas_id_ch4         )
          coef%ref_prfl_mr(1,igaz) = coef%ref_prfl_mr(2,igaz)
        End Select
      EndDo

    Case default
      coef%ref_prfl_p(1) = 0.500000E-02_jprb
      coef%ref_prfl_t(1,:) = coef%ref_prfl_t(2,:)
      Do igaz = 1, coef%fmv_gas
        Select case(coef%fmv_gas_id(igaz))
        Case (gas_id_mixed       )
          coef%ref_prfl_mr(1,igaz) = -.999900E+04_jprb
        Case (gas_id_watervapour )
          coef%ref_prfl_mr(1,igaz) = coef%ref_prfl_mr(2,igaz)
        Case (gas_id_ozone       )
          coef%ref_prfl_mr(1,igaz) = coef%ref_prfl_mr(2,igaz)
        Case (gas_id_wvcont      )
          coef%ref_prfl_mr(1,igaz) = coef%ref_prfl_mr(2,igaz)
        Case (gas_id_co2         )
          coef%ref_prfl_mr(1,igaz) = coef%ref_prfl_mr(2,igaz)
        Case (gas_id_n2o         )
          coef%ref_prfl_mr(1,igaz) = coef%ref_prfl_mr(2,igaz)
        Case (gas_id_co          )
          coef%ref_prfl_mr(1,igaz) = coef%ref_prfl_mr(2,igaz)
        Case (gas_id_ch4         )
          coef%ref_prfl_mr(1,igaz) = coef%ref_prfl_mr(2,igaz)
        End Select
      EndDo

  End Select


  ! Top level profile limits
  Select Case (dataset)

    Case ("52")
      coef%lim_prfl_p(1) = coef%ref_prfl_p(1)
      coef%lim_prfl_tmax(1) = 272.50_jprb
      coef%lim_prfl_tmin(1) = 130.81_jprb
      Do igaz = 1, coef%fmv_gas
        Select case(coef%fmv_gas_id(igaz))
        Case (gas_id_mixed       )
          coef%lim_prfl_gmax(1,igaz) = -.999900E+04_jprb
          coef%lim_prfl_gmin(1,igaz) = -.999900E+04_jprb
        Case (gas_id_watervapour )
          coef%lim_prfl_gmax(1,igaz) = 0.2868E+01_jprb
          coef%lim_prfl_gmin(1,igaz) = 0.1286E+01_jprb
        Case (gas_id_ozone       )
          coef%lim_prfl_gmax(1,igaz) = 0.1407E+01_jprb
          coef%lim_prfl_gmin(1,igaz) = 0.4828E-03_jprb
        Case (gas_id_wvcont      )
          coef%lim_prfl_gmax(1,igaz) = 0.2868E+01_jprb
          coef%lim_prfl_gmin(1,igaz) = 0.1286E+01_jprb
        Case (gas_id_co2         )
          coef%lim_prfl_gmax(1,igaz) = 0.4832E+03_jprb
          coef%lim_prfl_gmin(1,igaz) = 0.2491E+03_jprb
        Case default
          coef%lim_prfl_gmax(1,igaz) = coef%ref_prfl_mr(1,igaz) + (coef%lim_prfl_gmax(2,igaz) - coef%ref_prfl_mr(2,igaz))
          coef%lim_prfl_gmin(1,igaz) = coef%ref_prfl_mr(1,igaz) + (coef%lim_prfl_gmin(2,igaz) - coef%ref_prfl_mr(2,igaz))
        End Select
      EndDo

    Case ("83")
      coef%lim_prfl_p(1) = coef%ref_prfl_p(1)
      coef%lim_prfl_tmax(1) = 245.95_jprb
      coef%lim_prfl_tmin(1) = 143.65_jprb
      Do igaz = 1, coef%fmv_gas
        Select case(coef%fmv_gas_id(igaz))
        Case (gas_id_mixed       )
          coef%lim_prfl_gmax(1,igaz) = -.999900E+04_jprb
          coef%lim_prfl_gmin(1,igaz) = -.999900E+04_jprb
        Case (gas_id_watervapour )
          coef%lim_prfl_gmax(1,igaz) = 0.5241E+01_jprb
          coef%lim_prfl_gmin(1,igaz) = 0.9133E+00_jprb
        Case (gas_id_ozone       )
          coef%lim_prfl_gmax(1,igaz) = 0.1404E+01_jprb
          coef%lim_prfl_gmin(1,igaz) = 0.1402E-01_jprb
        Case (gas_id_wvcont      )
          coef%lim_prfl_gmax(1,igaz) = 0.5241E+01_jprb
          coef%lim_prfl_gmin(1,igaz) = 0.9133E+00_jprb
        Case (gas_id_co2         )
          coef%lim_prfl_gmax(1,igaz) = 0.4610E+03_jprb
          coef%lim_prfl_gmin(1,igaz) = 0.2933E+03_jprb
        Case (gas_id_n2o         )
          coef%lim_prfl_gmax(1,igaz) = 0.7563E-01_jprb
          coef%lim_prfl_gmin(1,igaz) = 0.2512E-04_jprb
        Case (gas_id_co          )
          coef%lim_prfl_gmax(1,igaz) = 0.8210E+01_jprb
          coef%lim_prfl_gmin(1,igaz) = 0.1216E-02_jprb
        Case (gas_id_ch4         )
          coef%lim_prfl_gmax(1,igaz) = 0.6238E+00_jprb
          coef%lim_prfl_gmin(1,igaz) = 0.6288E-03_jprb
        End Select
      EndDo

    Case default
      coef%lim_prfl_p(1) = coef%ref_prfl_p(1)
      coef%lim_prfl_tmax(1) = coef%ref_prfl_t(1,1) + (coef%lim_prfl_tmax(2) - coef%ref_prfl_t(2,1))
      coef%lim_prfl_tmin(1) = coef%ref_prfl_t(1,1) + (coef%lim_prfl_tmin(2) - coef%ref_prfl_t(2,1))
      coef%lim_prfl_gmax(1,:) = coef%ref_prfl_mr(1,:) + (coef%lim_prfl_gmax(2,:) - coef%ref_prfl_mr(2,:))
      coef%lim_prfl_gmin(1,:) = coef%ref_prfl_mr(1,:) + (coef%lim_prfl_gmin(2,:) - coef%ref_prfl_mr(2,:))
      Do igaz = 1, coef%fmv_gas
        Select case(coef%fmv_gas_id(igaz))
        Case (gas_id_mixed       )
          coef%lim_prfl_gmin(1,igaz) = -.999900E+04_jprb
          coef%lim_prfl_gmax(1,igaz) = -.999900E+04_jprb
        End Select
      EndDo

  End Select

  Close( file_id )

  Open (file_id, file = f_coef_out, iostat = ioerr, form = 'formatted')
  THROWM(ioerr.ne.0,"Cannot open "//Trim(f_coef_out))

  Call Rttov_write_ascii_coef( errorstatus, coef, file_id )
  THROWM(errorstatus.ne.0,"Cannot read "//Trim(f_coef_out))

  Close( file_id )

  Call rttov_dealloc_coef( errorstatus, coef )
  THROW(errorstatus.ne.0)

PCATCH

  Contains

  elemental subroutine rttov_lower_case(ous,ins)
  ! convert a word to lower case
  character (len=*) , intent(out) :: ous
  character (len=*) , intent(in) :: ins
  integer :: i,ic,nlen
  nlen = len(ins)
  ous = ''
  do i=1,nlen
     ic = ichar(ins(i:i))
     if (ic >= 65 .and. ic < 90) then
       ous(i:i) = char(ic+32)
     else
       ous(i:i) = ins(i:i)
     endif
  end do
  end subroutine rttov_lower_case

  Character(Len=2) Function grok_profile_dataset( lbl_id )
    Character(Len=132) :: lbl_id(20)
    Character(Len=132) :: lbl_id1(20)
    Integer(Kind=jpim) :: i

    grok_profile_dataset = "  "

    call rttov_lower_case(lbl_id1, lbl_id)

    Do i = 1, 20
      If( index(lbl_id1(i),"60l_sd  diverse_52profiles_60l") .ne. 0 ) Then
        grok_profile_dataset = "52"
        Return
      Endif
      If( index(lbl_id1(i),"tigr-43") .ne. 0 ) Then
        grok_profile_dataset = "43"
        Return
      Endif
      If( index(lbl_id1(i),"nwp-saf 52") .ne. 0 ) Then
        grok_profile_dataset = "52"
        Return
      Endif
      If( index(lbl_id1(i),"ecmwf_p83_91lev") .ne. 0 ) Then
        grok_profile_dataset = "83"
        Return
      Endif
      If( index(lbl_id1(i),"ecmwf_83p") .ne. 0 ) Then
        grok_profile_dataset = "83"
        Return
      Endif
    EndDo

    Do i = 1, 20
      If (.not. Trim(lbl_id(i)) == 'xxxx') Print *, Trim(lbl_id(i))
    EndDo
    Print *, "Cannot recognise profile dataset for top level values"

  End Function

End Program rttov789_conv_coef
