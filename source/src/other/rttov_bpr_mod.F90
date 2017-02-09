MODULE RTTOV_BPR_MOD
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 25 November 1998, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, KNMI and MeteoFrance.
!
!    Copyright 2012, EUMETSAT, All Rights Reserved.
!
! Current Code Owner: SAF NWP

USE PARKIND1,    ONLY : JPRB, JPIM
USE RTTOV_CONST, ONLY : PI, DEG2RAD

IMPLICIT NONE

INTEGER(KIND=JPIM) :: NPHANGLE 
INTEGER(KIND=JPIM) :: PHANGIND90 ! INDEX OF PHASE ANGLE 90 DEGREES

LOGICAL :: PHASE_INIT=.FALSE.

INTEGER(KIND=JPIM),PARAMETER :: VALS_PER_DEG=100

REAL(KIND=JPRB), POINTER :: ARX(:)
REAL(KIND=JPRB), POINTER :: XARR0(:)
REAL(KIND=JPRB), POINTER :: CXARR0(:)
REAL(KIND=JPRB), POINTER :: MUX(:)
REAL(KIND=JPRB) :: CPHI(1:360_JPIM)
REAL(KIND=JPRB) :: XE(0:180*VALS_PER_DEG)

REAL(KIND=JPRB) :: TACOS1( 0:6400)
REAL(KIND=JPRB) :: TACOS2( 0:35250)
REAL(KIND=JPRB) :: TACOS3( 0:7420)
REAL(KIND=JPRB) :: TACOS4( 0:790)
REAL(KIND=JPRB) :: TACOS5( 0:100 )
REAL(KIND=JPRB) :: OFFACOS1, OFFACOS2,OFFACOS3,OFFACOS4,OFFACOS5

END MODULE
