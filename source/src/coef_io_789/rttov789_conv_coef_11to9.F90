! Description:
!> @file
!!   Executable for converting RTTOV v10/v11-compatible
!!   coefficient files to v9 format.
!
!> @brief
!!   Executable for converting RTTOV v10/v11-compatible
!!   coefficient files to v9 format.
!!
!! @details
!!   Usage:
!!   $ rttov789_conv_coef_11to9.exe \-\-coef-in ... \-\-coef-out ...
!!
!!   where \-\-coef-in specifies the input v10/v11-compatible
!!   file and \-\-coef-out specifies the output v9-
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
Program rttov789_conv_coef_11to9

#include "throw.h"

  Use rttov_types, Only :  &
        & rttov_coef

  Use rttov_const, Only :  &
        & version,         &
        & release,         &
        & minor_version

  Use rttov_getoptions

  Use rttov_unix_env, Only: rttov_iargc, rttov_exit

  Use parkind1, Only : jpim, jprb

  Implicit None

#include "rttov_read_ascii_coef.interface"
#include "rttov789_write_ascii_coef.interface"
#include "rttov_errorreport.interface"

  Type(rttov_coef)   :: coef

  Integer(Kind=jpim) :: file_id
  Integer(Kind=jpim) :: errorstatus
  Character(Len=256) :: f_coef_in = ""
  Character(Len=256) :: f_coef_out = ""
  Character(Len=256) :: errmessage
  Integer :: ioerr, err
  Real(Kind=Jprb), Pointer :: X1(:), X2(:,:)
  Integer(Kind=jpim) :: nlevels, n2


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
  file_id = 77
  Open( file_id, File = f_coef_in, form = 'formatted', &
    Status = 'old', Iostat = ioerr )

  THROWM(ioerr.ne.0,"Cannot open "//Trim(f_coef_in))

  Call rttov_read_ascii_coef( errorstatus, coef, file_id )

  THROWM(errorstatus.ne.0,"Cannot read "//Trim(f_coef_in))

  Write( errMessage, '( "RTTOV library version ",i2,1x,i1,".",i1 )' )&
           & version, release, minor_version
  INFO(errMessage)
  Write( errMessage, '( "fast model version compatibility ",i2 )' ) coef%fmv_model_ver
  INFO(errMessage)

  ! Remove top level
  ! File format version is 9
  ! Gas units are unchanged (so will usually be ppmv)
  ! "IncTop" set to 0
  ! FASTEM version set to 3

  coef%nlevels     = coef%nlevels - 1
  coef%fmv_lvl     = coef%nlevels
  coef%id_comp_pc  = 0
  coef%id_comp_lvl = 9
  IF (coef%fastem_ver > 0) coef%fastem_ver  = 3

  nlevels = coef%nlevels

  Allocate(X1(nlevels)); X1 = coef%ref_prfl_p(2:)
  deallocate(coef%ref_prfl_p); coef%ref_prfl_p => X1

  n2 = size(coef%ref_prfl_t,2)
  Allocate(X2(nlevels,n2)); X2 = coef%ref_prfl_t(2:,:)
  deallocate(coef%ref_prfl_t); coef%ref_prfl_t => X2

  n2 = size(coef%ref_prfl_mr,2)
  Allocate(X2(nlevels,n2)); X2 = coef%ref_prfl_mr(2:,:)
  deallocate(coef%ref_prfl_mr); coef%ref_prfl_mr => X2

  Allocate(X1(nlevels)); X1 = coef%lim_prfl_p(2:)
  deallocate(coef%lim_prfl_p); coef%lim_prfl_p => X1

  Allocate(X1(nlevels)); X1 = coef%lim_prfl_tmax(2:)
  deallocate(coef%lim_prfl_tmax); coef%lim_prfl_tmax => X1

  Allocate(X1(nlevels)); X1 = coef%lim_prfl_tmin(2:)
  deallocate(coef%lim_prfl_tmin); coef%lim_prfl_tmin => X1

  n2 = size(coef%lim_prfl_gmax,2)
  Allocate(X2(nlevels,n2)); X2 = coef%lim_prfl_gmax(2:,:)
  deallocate(coef%lim_prfl_gmax); coef%lim_prfl_gmax => X2

  n2 = size(coef%lim_prfl_gmin,2)
  Allocate(X2(nlevels,n2)); X2 = coef%lim_prfl_gmin(2:,:)
  deallocate(coef%lim_prfl_gmin); coef%lim_prfl_gmin => X2

  Close( file_id )

  Open (file_id, file = f_coef_out, iostat = ioerr, form = 'formatted')
  THROWM(ioerr.ne.0,"Cannot open "//Trim(f_coef_out))

  Call rttov789_write_ascii_coef( errorstatus, coef, file_id )
  THROWM(errorstatus.ne.0,"Cannot read "//Trim(f_coef_out))

  Close( file_id )

  Call rttov_dealloc_coef( errorstatus, coef )
  THROW(errorstatus.ne.0)

PCATCH

End Program rttov789_conv_coef_11to9
