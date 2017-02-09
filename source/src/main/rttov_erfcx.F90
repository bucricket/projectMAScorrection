!     Compute complementary error function
      FUNCTION RTTOV_ERFCX(X)

!     Description:
!     RTTOV_ERFCX - To compute complementary error function
!
!     Method:
!     Chebyshev fitting
!
!     Owner:
!     Argonne National Laboratory
!
!     History:
!     Version      Date        Comment
!     1            15/07/2003  Marco Matricardi. ECMWF.
!                              Based on the routine CALREF written by
!                              W.J.Cody, Mathematics and Computer Science
!                              Division Argonne National Laboratory
!                              Argonne, IL 60439
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code"


      Use parkind1   , Only :   &
!     Imported parameters:
     & jprb

!INTF_OFF
Use parkind1, Only : &
    jpim
!INTF_ON



      IMPLICIT NONE

!     Function arguments:

!       Scalar arguments with intent in:
      REAL (Kind=jprb)   ,INTENT(IN) :: X         ! Interpolation argument.
      REAL (Kind=jprb) :: RTTOV_ERFCX

!INTF_END

!     End of function arguments


!       Local arrays
      REAL    (Kind=jprb) :: A(5),B(4),C(9),D(8),P(6),Q(5)
!       Local scalars
      REAL    (Kind=jprb) :: Y,THRESH,YSQ,XSMALL,XNUM,XDEN,XBIG,SQRPI,DEL
      INTEGER (Kind=jpim) :: I


!-----End of header-------------------------------------------------------------

      DATA A/3.16112374387056560E00_JPRB,1.13864154151050156E02_JPRB,          &
           & 3.77485237685302021E02_JPRB,3.20937758913846947E03_JPRB,            &
           & 1.85777706184603153E-1_JPRB/
      DATA B/2.36012909523441209E01_JPRB,2.44024637934444173E02_JPRB,          &
           & 1.28261652607737228E03_JPRB,2.84423683343917062E03_JPRB/
      DATA C/5.64188496988670089E-1_JPRB,8.88314979438837594E0_JPRB,           &
           & 6.61191906371416295E01_JPRB,2.98635138197400131E02_JPRB,            &
           & 8.81952221241769090E02_JPRB,1.71204761263407058E03_JPRB,            &
           & 2.05107837782607147E03_JPRB,1.23033935479799725E03_JPRB,            &
           & 2.15311535474403846E-8_JPRB/
      DATA D/1.57449261107098347E01_JPRB,1.17693950891312499E02_JPRB,          &
           & 5.37181101862009858E02_JPRB,1.62138957456669019E03_JPRB,            &
           & 3.29079923573345963E03_JPRB,4.36261909014324716E03_JPRB,            &
           & 3.43936767414372164E03_JPRB,1.23033935480374942E03_JPRB/
      DATA P/3.05326634961232344E-1_JPRB,3.60344899949804439E-1_JPRB,          &
           & 1.25781726111229246E-1_JPRB,1.60837851487422766E-2_JPRB,            &
           & 6.58749161529837803E-4_JPRB,1.63153871373020978E-2_JPRB/
      DATA Q/2.56852019228982242E00_JPRB,1.87295284992346047E00_JPRB,          &
           & 5.27905102951428412E-1_JPRB,6.05183413124413191E-2_JPRB,            &
           & 2.33520497626869185E-3_JPRB/
      DATA THRESH /0.46875_JPRB/
      DATA XSMALL /5.96E-8_JPRB/
      DATA XBIG   /9.194_JPRB  /
      DATA SQRPI  /5.6418958354775628695E-1_JPRB/

      Y=ABS(X)

      IF (Y .LE. THRESH) THEN

        YSQ =0.

        IF (Y .GT. XSMALL)THEN
          YSQ = Y * Y
        ENDIF

        XNUM = A(5)*YSQ
        XDEN = YSQ

          DO I = 1, 3
            XNUM = (XNUM + A(I)) * YSQ
            XDEN = (XDEN + B(I)) * YSQ
          ENDDO

        RTTOV_ERFCX=1-X * (XNUM + A(4)) / (XDEN + B(4))

      ELSE IF (Y .LE. 4.) THEN

        XNUM = C(9)*Y
        XDEN = Y

        DO I = 1, 7
          XNUM = (XNUM + C(I)) * Y
          XDEN = (XDEN + D(I)) * Y
        ENDDO

        RTTOV_ERFCX= (XNUM + C(8)) / (XDEN + D(8))
        YSQ = AINT(Y*16.)/16.
        DEL = (Y-YSQ)*(Y+YSQ)
        RTTOV_ERFCX= EXP(-YSQ*YSQ) * EXP(-DEL) *RTTOV_ERFCX

      ELSE

        IF (Y .GE. XBIG) THEN
          RTTOV_ERFCX=0.
        ELSE
          YSQ = 1. / (Y * Y)
          XNUM = P(6)*YSQ
          XDEN = YSQ
            DO I = 1, 4
              XNUM = (XNUM + P(I)) * YSQ
              XDEN = (XDEN + Q(I)) * YSQ
            ENDDO
          RTTOV_ERFCX = YSQ *(XNUM + P(5)) / (XDEN + Q(5))
          RTTOV_ERFCX = (SQRPI -  RTTOV_ERFCX) / Y
          YSQ = AINT(Y*16.)/16.
          DEL = (Y-YSQ)*(Y+YSQ)
          RTTOV_ERFCX = EXP(-YSQ*YSQ) * EXP(-DEL) * RTTOV_ERFCX
        ENDIF
      ENDIF

      IF(X<0)THEN
        RTTOV_ERFCX=2.-RTTOV_ERFCX
      ENDIF



      END FUNCTION RTTOV_ERFCX
