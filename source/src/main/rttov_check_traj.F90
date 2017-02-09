Subroutine rttov_check_traj( &
       & ERR,             &
       & nprofiles,       &
       & nchannels,       &
       & opts,            &
       & nlevels,         &
       & coefs,           &
       & asw,             &
       & traj0,    traj1,    traj2,    &
       & traj0_tl, traj1_tl, traj2_tl, &
       & traj0_ad, traj1_ad, traj2_ad, &
       & traj0_k,  traj1_k,  traj2_k   &
       & )
! Description:
!   Allocates/deallocates trajectory data for rttov_direct
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
!INTF_OFF
!INTF_OFF
#include "throw.h"
!INTF_ON
!INTF_ON
  Use parkind1, only : jpim

!INTF_OFF
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
  Use parkind1, only : jprb
!INTF_ON

  Use rttov_types, Only : &
          & rttov_coefs,    &
          & rttov_options,  &
          & rttov_traj

  Implicit None

  Integer(Kind=jpim),       Intent(in)    :: nprofiles
  Integer(Kind=jpim),       Intent(in)    :: nchannels

  Type(rttov_coefs),       Intent(in), Target :: coefs
  Type(rttov_options),     Intent(in)    :: opts
  Integer(Kind=jpim),       Intent(in)   :: nlevels

  Integer(Kind=jpim),       Intent(out)   :: ERR
  Integer(Kind=jpim),       Intent(in)    :: asw
  Type(rttov_traj),   Pointer, Optional                :: traj0, traj0_tl, traj0_ad, traj0_k
  Type(rttov_traj),   Target,  Optional, Intent(inout) :: traj1, traj1_tl, traj1_ad, traj1_k
  Type(rttov_traj),   Target,  Optional, Intent(inout) :: traj2, traj2_tl, traj2_ad, traj2_k

!INTF_END

#include "rttov_errorreport.interface"

REAL(KIND=JPRB) :: ZHOOK_HANDLE
TRY

IF (LHOOK) CALL DR_HOOK('RTTOV_CHECK_TEMP',0_jpim,ZHOOK_HANDLE)
  If( present( traj0 ) ) THEN
      Call check_traj(        &
       & ERR,                &
       & nprofiles,          &
       & nchannels,          &
       & opts,               &
       & nlevels,            &
       & coefs,              &
       & asw,                &
       & traj0,              &
       & traj1,              &
       & traj2 )
  THROWM( ERR .NE. 0 , "check traj0")
  ENDIF

  If( present( traj0_tl ) ) THEN
      Call check_traj(        &
       & ERR,                &
       & nprofiles,          &
       & nchannels,          &
       & opts,               &
       & nlevels,            &
       & coefs,              &
       & asw,                &
       & traj0_tl,           &
       & traj1_tl,           &
       & traj2_tl )
  THROWM( ERR .NE. 0 , "check traj0_tl")
  ENDIF


  If( present( traj0_ad ) ) THEN
      Call check_traj(       &
       & ERR,                &
       & nprofiles,          &
       & nchannels,          &
       & opts,               &
       & nlevels,            &
       & coefs,              &
       & asw,                &
       & traj0_ad,           &
       & traj1_ad,           &
       & traj2_ad )
  THROWM( ERR .NE. 0 , "check traj0_ad")
  ENDIF


  If( present( traj0_k ) )  THEN
      Call check_traj(       &
       & ERR,                &
       & nchannels,          &
       & nchannels,          &
       & opts,               &
       & nlevels,            &
       & coefs,              &
       & asw,                &
       & traj0_k,            &
       & traj1_k,            &
       & traj2_k )
   THROWM( ERR .NE. 0 , "check traj0_k")
   ENDIF

IF (LHOOK) CALL DR_HOOK('RTTOV_CHECK_TEMP',1_jpim,ZHOOK_HANDLE)

CATCH
IF (LHOOK) CALL DR_HOOK('RTTOV_CHECK_TEMP',1_jpim,ZHOOK_HANDLE)
  Contains

  Subroutine check_traj(  &
       & ERR,             &
       & nprofiles,       &
       & nchannels,       &
       & opts,            &
       & nlevels,         &
       & coefs,           &
       & asw,             &
       & traj0,           &
       & traj1,           &
       & traj2 )

  Integer(Kind=jpim),       Intent(in)    :: nprofiles
  Integer(Kind=jpim),       Intent(in)    :: nchannels
  Type(rttov_options),      Intent(in)    :: opts

  Integer(Kind=jpim),       Intent(in)    :: nlevels
  Type(rttov_coefs),        Intent(in), Target    :: coefs

  Integer(Kind=jpim),       Intent(out)   :: ERR
  Integer(Kind=jpim),       Intent(in)    :: asw
  Type(rttov_traj),   Pointer :: traj0
  Type(rttov_traj),   Target, Optional, Intent(inout) :: traj1
  Type(rttov_traj),   Target, Optional, Intent(in)    :: traj2

#include "rttov_alloc_traj.interface"
#include "rttov_opts_eq.interface"


REAL(KIND=JPRB) :: ZHOOK_HANDLE

TRY
IF (LHOOK) CALL DR_HOOK('CHECK_TEMP',0_jpim,ZHOOK_HANDLE)
  If( asw .eq. 1_jpim ) Then
    If( Present( traj2 ) ) Then

      If( ( .not. associated( traj2%coefs, coefs ) )              .or. &
         & ( traj2%nchannels .ne. nchannels )                     .or. &
         & ( .not. rttov_opts_eq (opts, traj2%opts) )             .or. &
         & ( traj2%nlevels   .ne. nlevels   ) ) Then
        err = errorstatus_fatal
        THROWM( ERR .NE. 0 , "rttov_check_traj fatal error dimensions mismatch")
      EndIf
      traj0 => traj2
    Else
      Call rttov_alloc_traj( &
         & ERR,               &
         & nprofiles,         &
         & nchannels,         &
         & opts,              &
         & nlevels,           &
         & coefs,             &
         & 1_jpim,            &
         & traj1 )
      traj0 => traj1
      THROWM( ERR .NE. 0 , "rttov_alloc_traj fatal error")

    EndIf
  Else
    If( Present( traj2 ) ) Then
      nullify(traj0)
    Else
      Call rttov_alloc_traj( &
         & ERR,               &
         & nprofiles,         &
         & nchannels,         &
         & opts,              &
         & nlevels,           &
         & coefs,             &
         & 0_jpim,            &
         & traj1 )
      THROWM( ERR .NE. 0 , "rttov_alloc_traj fatal error")
      traj0 => traj1
    EndIf
  EndIf
IF (LHOOK) CALL DR_HOOK('CHECK_TEMP',1_jpim,ZHOOK_HANDLE)

CATCH
IF (LHOOK) CALL DR_HOOK('CHECK_TEMP',1_jpim,ZHOOK_HANDLE)
  End Subroutine

End Subroutine
