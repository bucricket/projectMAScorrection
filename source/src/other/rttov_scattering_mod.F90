MODULE rttov_scattering_mod

  ! Description:
  !   Helper routines for scattering calculations.
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
  !    Copyright 2012, EUMETSAT, All Rights Reserved.
  !
  ! Method:
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date        Comment
  ! -------   ----        -------
  !  1.0      11/09/2012  Created (J Hocking)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:
  !
#include "throw.h"

IMPLICIT NONE

CONTAINS

subroutine integrate(n, x, f, intg, error, fail)

  ! This replaces the NAG library d01GAF subroutine to avoid dependency on an external
  ! library. The differences in the calculated quantities between this and d01GAF are
  ! in the least significant figure.

  ! Uses the method of Gill and Miller 1972 (I copied the code from the paper) and adds
  ! the estimated error to the output value for consistency with d01GAF (the NAG folks
  ! found this generally results in better agreement with the true value)

  use parkind1, only : jprb, jpim
  implicit none
  
  integer(kind=jpim), intent(in)  :: n
  real(kind=jprb),    intent(in)  :: x(1:n)
  real(kind=jprb),    intent(in)  :: f(1:n)
  real(kind=jprb),    intent(out) :: intg
  real(kind=jprb),    intent(out) :: error
  integer(kind=jpim), intent(out) :: fail
  
  integer(kind=jpim) :: i, j, k
  real(kind=jprb)    :: h1, h2, h3, h4, r1, r2, r3, r4, d1, d2, d3, c, s
  
  fail = 0
  
  if (n < 4) then
    fail = 1
    print *,'Requires at least 4 data points'
    stop
  endif
  
  intg  = 0._jprb
  error = 0._jprb
  s     = 0._jprb
  c     = 0._jprb
  r4    = 0._jprb
  
  j = 3
  k = n - 1
  do i = j, k
    if (i == j) then
      h2 = x(j-1) - x(j-2)
      d3 = (f(j-1) - f(j-2))/h2
      h3 = x(j) - x(j-1)
      d1 = (f(j) - f(j-1))/h3
      h1 = h2 + h3
      d2 = (d1 - d3)/h1
      h4 = x(j+1) - x(j)
      r1 = (f(j+1) - f(j))/h4
      r2 = (r1 - d1)/(h4 + h3)
      h1 = h1 + h4
      r3 = (r2 - d2)/h1
      intg = h2 * (f(1) + h2 * (d3/2._jprb - h2 * (d2/6._jprb - (h2 + 2._jprb * h3) * r3/12._jprb)))
      s    = -h2**3 * (h2 * (3._jprb * h2 + 5._jprb * h4) + 10._jprb * h3 * h1) / 60._jprb
    else
      h4 = x(i+1) - x(i)
      r1 = (f(i+1) - f(i))/h4
      r4 = h4 + h3
      r2 = (r1 - d1)/r4
      r4 = r4 + h2
      r3 = (r2 - d2)/r4
      r4 = r4 + h1
      r4 = (r3 - d3)/r4
    endif
    
    intg = intg + h3 * ((f(i) + f(i-1))/2._jprb - h3 * h3 * (d2 + r2 + (h2 - h4) * r3)/12._jprb)
    c    = h3**3 * (2._jprb * h3 * h3 +  5._jprb * (h3 * (h4 + h2) + 2._jprb * h4 * h2))/120._jprb
    error = error + (c + s) * r4
    if (i == j) then
      s = 2._jprb*c + s
    else
      s = c 
    endif
    
    if (i == k) then
      intg = intg + h4 * (f(n) - h4 * (r1/2._jprb + h4 * (r2/6._jprb + (2._jprb * h3 + h4) * r3/12._jprb)))
      error = error - h4**3 * r4 * (h4 * (3._jprb * h4 + 5._jprb * h2) + 10._jprb * h3 * (h2 + h3 + h4))/60._jprb
      error = error + s * r4
    else
      h1 = h2
      h2 = h3
      h3 = h4
      d1 = r1
      d2 = r2
      d3 = r3
    endif
  enddo
  intg = intg + error  ! The error is added to the answer in the NAG lib function
end subroutine integrate


subroutine mie_sphere(x, kl, r_angle, m, q_sct, q_ext,f, g, q_bsct, phasef)!, iphase)

!     Author: Marco Matricardi, ECMWF, 17-5-2012

!     Some lines not required by rttov_mie_params.F90 commented out

      use parkind1, only : jprb, jpim!, jplm
!       use rttov_const, only : pi
      implicit none

!-----Common variables
      real(kind=jprb)   , intent (in)  :: x, kl, r_angle             ! The size parameter
      complex(kind=jprb), intent (in)  :: m                          ! The refractive index
      real(kind=jprb)   , intent (out) :: q_sct, q_ext,f, g, q_bsct,phasef
!       integer(kind=jpim), intent (in)  :: iphase

!-----Local variables
      real(kind=jprb)    :: n1, n2, d0_re, d0_im, f_n, g_n
!       real(kind=jprb)    :: sum
      complex(kind=jprb) :: mc, mx, c_bsct, ss1,ss2
!       complex(kind=jprb) :: anbn, bnan, z
      integer(kind=jpim) :: n, n_end
!       integer(kind=jpim) :: nn
!       logical(kind=jplm) :: evenl

      complex(kind=jprb), allocatable :: dn (:), wn (:), an (:), bn (:),s1n(:),s2n(:)
      real(kind=jprb)   , allocatable :: pin(:),taun(:)
!       complex(kind=jprb), allocatable :: cm(:),c(:),dm(:),d(:)
!       real(kind=jprb)   , allocatable :: recip(:),bidel(:),pmom(:,:)
!       real(kind=jprb)   , allocatable :: am(:),bi(:)

!       integer(kind=jpim) :: ntrm,k,l,ld2,idel,i,mm,imax,mmax

!-----How many iterations
      n_end = anint (x + 4.05_jprb * x ** (1.0_jprb / 3.0_jprb) + 2.0_jprb)

      if (n_end <    5) n_end =    5
!      if (n_end > 1000) n_end = 1000

      allocate (dn (-1:n_end+1))
      allocate (wn (-1:n_end+1))
      allocate (an (-1:n_end+1))
      allocate (bn (-1:n_end+1))
      allocate (s1n(-1:n_end+1))
      allocate (s2n(-1:n_end+1))
      allocate (pin(-1:n_end+1))
      allocate (taun(-1:n_end+1))
!       allocate (cm (1:n_end+3))
!       allocate (c  (1:n_end+3))
!       allocate (dm (1:n_end+3))
!       allocate (d  (1:n_end+3))
!       allocate (recip(1:4*(n_end+3)+2))
!       allocate (am (0:n_end+3))
!       allocate (bi (0:n_end+3))
!       allocate (bidel (0:n_end+3))
!       allocate (pmom(0:3,4))

!       pmom(:,:)=0._jprb

      mc = conjg (m)
      mx = mc * x
      n1 = real(m) * x
      n2 = aimag(m) * x

      d0_re = sin  (n1) * cos  (n1) / (sin (n1) * sin (n1) + sinh              & 
              (n2) * sinh (n2))
      d0_im = sinh (n2) * cosh (n2) / (sin (n1) * sin (n1) + sinh              &
              (n2) * sinh (n2))

!-----Initialize

      pin(0)  =0._jprb
      pin(1)  =1._jprb
      pin(2)  =3._jprb*cos(r_angle)
      taun(0) =0._jprb
      taun(1) =cos(r_angle)
      taun(2) =3._jprb*cos(2._jprb*r_angle)


      dn ( 0) = cmplx (d0_re  ,             d0_im, kind=jprb)
      wn (-1) = cmplx (cos (x), -1.0_jprb * sin (x), kind=jprb)
      wn ( 0) = cmplx (sin (x),             cos (x), kind=jprb)

      q_ext  = 0.0_jprb
      q_sct  = 0.0_jprb
      q_bsct = 0.0_jprb
      f      = 0.0_jprb
      g      = 0.0_jprb
      c_bsct = cmplx (0.0_jprb,0.0_jprb,jprb)
      ss1    = (0._jprb,0._jprb)
      ss2    = (0._jprb,0._jprb)


      do n = 1, n_end + 1

        if(n>2)then
          pin(n)=((2._jprb*n-1._jprb)*cos(r_angle)*pin(n-1)/(n-1._jprb))-n*pin(n-2)/(n-1._jprb)
          taun(n)=n*cos(r_angle)*pin(n)-(n+1._jprb)*pin(n-1)
        endif


        f_n = 2.0_jprb * n - 1.0_jprb
        g_n = 2.0_jprb * n + 1.0_jprb

        wn(n) = f_n / x * wn(n - 1) - wn(n - 2)
        dn(n) = 1.0_jprb / (n / mx - dn(n - 1)) - n / mx

        an(n) = ((dn(n) / mc + n / x) * real(wn(n)) - real(wn(n - 1))) / &
                ((dn(n) / mc + n / x) * wn(n) - wn(n - 1))
        bn(n) = ((dn(n) * mc + n / x) * real(wn(n)) - real(wn(n - 1))) / &
                ((dn(n) * mc + n / x) * wn(n) - wn(n - 1))

        S1n(n)=(2._jprb*n+1._jprb)*(an(n)*pin(n)+bn(n)*taun(n))/(n*(n+1._jprb))
        S2n(n)=(2._jprb*n+1._jprb)*(bn(n)*pin(n)+an(n)*taun(n))/(n*(n+1._jprb))

        q_ext  = q_ext  + g_n * (real(an(n)) + real(bn(n)))
        q_sct  = q_sct  + g_n * (abs(an(n)) ** 2 + abs(bn(n)) ** 2)
        c_bsct = c_bsct + g_n * (-1) ** n * (an(n) - bn(n))

        ss1=ss1+S1n(n)
        ss2=ss2+S2n(n)

!         if (n > 1) then
!           nn = n - 1
! 
!           anbn = an (nn) * conjg (an (n)) + bn (nn) * conjg (bn (n))
!           bnan = an (nn) * conjg (bn (nn))
! 
!           g = g + nn * (nn + 2.0_jprb) / (nn + 1.0_jprb) * real(anbn) + (2.0_jprb * nn +     &
!               1.0_jprb) / (nn * (nn + 1.0_jprb)) * real(bnan)
!         end if
      end do

!       if(iphase==0)then
! 
! !-----Calculate Mueller C and D arrays-----------------------------------------
! 
!       do k=1,4*(n_end+3)+2
!         recip(k)=1./k
!       enddo
! 
!       ntrm=n_end+1
! 
!       cm( ntrm + 2 ) = ( 0., 0. )
!       dm( ntrm + 2 ) = ( 0., 0. )
!       cm( ntrm + 1 ) = ( 1. - recip( ntrm+1 ) ) * bn( ntrm )
!       dm( ntrm + 1 ) = ( 1. - recip( ntrm+1 ) ) * an( ntrm ) 
!       cm( ntrm ) = ( recip( ntrm ) + recip( ntrm+1 ) ) * an( ntrm ) +          &
!                    ( 1. - recip( ntrm ) )*bn( ntrm-1 ) 
!       dm( ntrm ) = ( recip( ntrm ) + recip( ntrm+1 ) ) * bn( ntrm ) +          &
!                    ( 1. - recip( ntrm ) )*an( ntrm-1 )
! 
!       do  k = ntrm-1, 2, -1
!         cm( k ) = cm( k+2 ) - ( 1. + recip(k+1) ) * bn( k+1 )                  &
!                             + ( recip(k) + recip(k+1) ) * an( k )              &
!                             + ( 1. - recip(k) ) * bn( k-1 )
!         dm( k ) = dm( k+2 ) - ( 1. + recip(k+1) ) * an( k+1 )                  &
!                             + ( recip(k) + recip(k+1) ) * bn( k )              &
!                             + ( 1. - recip(k) ) * an( k-1 )
!       enddO
!       cm( 1 ) = cm( 3 ) + 1.5 * ( an( 1 ) - bn( 2 ) )
!       dm( 1 ) = dm( 3 ) + 1.5 * ( bn( 1 ) - an( 2 ) )
! 
!       do k = 1, ntrm + 2
!         c( k ) = ( 2*k - 1 ) * cm( k )
!         d( k ) = ( 2*K - 1 ) * DM( K )
!       enddo
! 
! !------------------------------------------------------------------------------
! 
! 
! !-----Compute coefficients of Legendre series----------------------------------
! 
!       do l=0,2
! 
!         ld2=l/2
!         evenl=mod( l, 2 ).eq.0
! 
!           if( l.eq.0 ) then
! 
!             idel = 1
! 
!               do mm = 0, ntrm
!                 am( mm ) = 2.0 * recip( 2*mm + 1 )
!               enddo
! 
!             bi( 0 ) = 1.0
! 
!           else if( evenl ) then
! 
!             idel = 1
! 
!               do mm = ld2, ntrm
!                 am( mm ) = ( 1. + recip( 2*mm - l + 1 ) ) * am( mm )
!               enddo
! 
!               do i = 0, ld2 - 1
!                 bi( i ) = ( 1. - recip( l - 2*i ) ) * bi( i )
!               enddo
! 
!             bi( ld2 ) = ( 2. - recip( l ) ) * bi( ld2 - 1 )
! 
!           else
! 
!             idel = 2
! 
!               do  mm = ld2, ntrm
!                 am( mm ) = ( 1. - recip( 2*mm + l + 2 ) ) * am( mm )
!               enddo
! 
!               do i = 0, ld2
!                 bi( i ) = ( 1. - recip( l + 2*i + 1 ) ) * bi( I )
!               enddo
! 
!           end if
! 
! 
!       mmax = ntrm - idel
! 
!       mmax = mmax + 1                             !if( ipolzn.ge.0 )
!       imax = min( ld2, mmax - ld2 )
! 
!       if( imax.lt.0 ) go to  250
! 
!       do  i = 0, imax
!         bidel( i ) = bi( i )
!       enddo
! 
!       if( evenl ) bidel( 0 ) = 0.5*bidel( 0 )
! 
!       do i = 0, imax
! 
!         sum = 0.0
!           do mm = ld2, mmax - i
!              sum = sum + am( mm ) * real( c(mm-i+1) * conjg( c(mm+i+idel) ), kind=jprb )
!           enddo
! 
!         pmom( l, 1 ) = pmom( l, 1 ) + bidel( i ) * sum
!       enddo
! 
!       do i = 0, imax
! 
!         sum = 0.0
!         do mm = ld2, mmax - i
!           sum = sum + am(mm ) *real( d(mm-i+1) * conjg( d(mm+i+idel) ), kind=jprb )
!         enddo
! 
!         pmom( l, 2 ) = pmom( l, 2 ) + bidel( i ) * sum
!       enddo
! 
! 
!     enddo
! 
! 250 continue
! 
!     endif



      phasef = (abs (ss1) ** 2 + abs (ss2)** 2)

      q_ext  = q_ext * 2.0_jprb / (x * x)
      q_sct  = q_sct * 2.0_jprb / (x * x)

      if (q_sct > q_ext) q_sct = q_ext

!       g      = g  * 4.0_jprb / (x * x * q_sct)
!       f      = (pmom( 2, 2 )+pmom( 2, 1 ))*2_jprb/(x**2*q_sct)
!       q_bsct = abs (c_bsct) * abs (c_bsct) / (x * x)
!      phasef = phasef*2*pi/(kl**2)
       phasef = phasef/(2._jprb*kl**2)

      deallocate (an, bn, dn,s1n,s2n, wn)
      deallocate (pin,taun)
!       deallocate (cm,c,dm,d,recip,am,bi,bidel,pmom)

      return
end subroutine mie_sphere

subroutine lognorm(ntot, rarr, asigma, armod, sqrt2pi, n)

      use parkind1, only : jprb, jpim
      implicit none

      integer(kind=jpim), intent(in)  :: ntot
      real(kind=jprb),    intent(in)  :: rarr(1:ntot)
      real(kind=jprb),    intent(in)  :: asigma
      real(kind=jprb),    intent(in)  :: armod
      real(kind=jprb),    intent(in)  :: sqrt2pi
      real(kind=jprb),    intent(out) :: n(1:ntot)
      
      integer(kind=jpim) :: i
      real(kind=jprb)    :: r

      do i = 1, ntot
        r = rarr(i)
        n(i) = (1 / (sqrt2pi * r * log10(asigma) * log(10._jprb))) * &
               exp(-0.5 * ((log10(r) - log10(armod)) / log10(asigma))**2)
      enddo

      return
end subroutine lognorm

subroutine gammadist(ntot, rarr, aacoef, aalpha, abcoef, agamma, n)

      use parkind1, only : jprb, jpim
      implicit none
      
      integer(kind=jpim), intent(in)  :: ntot
      real(kind=jprb),    intent(in)  :: rarr(1:ntot)
      real(kind=jprb),    intent(in)  :: aacoef
      real(kind=jprb),    intent(in)  :: aalpha
      real(kind=jprb),    intent(in)  :: abcoef
      real(kind=jprb),    intent(in)  :: agamma
      real(kind=jprb),    intent(out) :: n(1:ntot)
      
      integer(kind=jpim) :: i
      real(kind=jprb)    :: r

      do i = 1, ntot
        r = rarr(i)
        n(i) = aacoef * (r**aalpha) * exp(-abcoef * (r**agamma))
      enddo

      return
end subroutine gammadist

!***********************************************************************
!  PROGRAM        INTER   SUBROUTINE
!-----------------------------------------------------------------------
!  PURPOSE        TO INTERPOLATE AT THE POINT ARG BETWEEN THE ARRAY ARX
!                 AND THE CORRESPONDING FUNCTION VALUES ARY.
!                 EXPONENTIAL INTERPOLATION IS USED FOR PRESSURE AND
!                 NUMBER DENSITY, LINEAR INTERPOLATION FOR TEMPERATURE.
!                 IF THE MODE NO IS <0 NEGATIVE NUMBERS ARE SET = 0.
!-----------------------------------------------------------------------
!  VERSION        4.0   D.P. EDWARDS   15/08/94
!                 Last changed 96/11/05
!-----------------------------------------------------------------------
!  ARGUMENTS      MXDIM   I*4  I/P  ARRAY DIMENSION
!                 NREC    I*4  I/P  NO OF ELEMENTS IN ARRAYS ARX AND ARY
!                 MODE    I*4  I/P  INTERPOLATION MODE
!                 ARG     R*4  I/P  INTERPOLATION ARGUMENT
!                 ARX     R*4  I/P  X VALUE ARRAY 
!                 ARY     R*4  I/P  Y VALUE FUNCTION ARRAY
!                 SS      R*4  O/P  INTERPOLATED FUNCTION VALUE AT ARG
!                 H       R*4  O/P  GRADIENT OR SCALE HEIGHT VALUE 
!***********************************************************************
       SUBROUTINE INTER(MXDIM,NREC,MODE,ARG,ARX,ARY,SS,H)
!-----------------------------------------------------------------------
       use parkind1, only : jprb, jpim
       
       INTEGER(KIND=JPIM)   MXDIM, NREC, MODE
       REAL(KIND=JPRB)      ARG
       REAL(KIND=JPRB)      ARX(NREC),ARY(NREC)
       REAL(KIND=JPRB)      SS, H
       
       INTEGER(KIND=JPIM)   IR, KL, NO, NOO, KS, KF, K
       REAL(KIND=JPRB)      AA, BB
!***********************************************************************
!
       IF (ARG .GE. ARX(1) .AND. ARG .LE. ARX(NREC)) THEN
         DO 10 IR=1,NREC-1
           IF (ARG .GE. ARX(IR) .AND. ARG .LT. ARX(IR+1)) KL = IR
   10    CONTINUE
         IF (ARG .EQ. ARX(NREC)) KL = NREC - 1
       ELSEIF (ARG .LT. ARX(1)) THEN
         KL = 1
       ELSEIF (ARG .GT. ARX(NREC)) THEN
         KL = NREC - 1
       ENDIF
!
!  INTERPOLATE FUNCTION VALUE AT ARG FROM DATA POINTS KL TO KL+1
!
!  EXPONENTIAL INTERPOLATION
!
       IF (ABS(MODE) .EQ. 1 .OR. ABS(MODE) .EQ. 4) THEN
         H = -(ARX(KL+1) - ARX(KL))/LOG(ARY(KL+1)/ARY(KL))
         SS = ARY(KL)*EXP(-(ARG-ARX(KL))/H)
         IF (ABS(MODE) .EQ. 4)  SS = 1.0/SS
         IF (MODE .LT. 0 .AND. SS .LT. 0.0) SS = 0.0
!
!  LINEAR INTERPOLATION
!
       ELSEIF (ABS(MODE) .EQ. 2 .OR. ABS(MODE) .EQ. 5) THEN
         H = (ARY(KL+1) - ARY(KL))/(ARX(KL+1) - ARX(KL))
         SS = ARY(KL) + H*(ARG - ARX(KL))
         IF (ABS(MODE) .EQ. 5) SS = 1.0/SS
         IF (MODE .LT. 0 .AND. SS .LT. 0.0) SS = 0.0
!
!  LINEAR-LOG INTERPOLATION
       ELSEIF (ABS(MODE) .EQ. 7) THEN
         H=(LOG(ARY(KL+1))-LOG(ARY(KL)))/(ARX(KL+1)-ARX(KL))
         SS = EXP(LOG(ARY(KL)) + H*(ARG - ARX(KL)))
         IF (MODE .LT. 0 .AND. SS .LT. 0.0) SS = 0.0
!
!  LOG-LOG INTERPOLATION
!
       ELSEIF (ABS(MODE) .EQ. 8) THEN
         H=(LOG(ARY(KL+1))-LOG(ARY(KL)))/(LOG(ARX(KL+1))-LOG(ARX(KL)))
         SS = EXP(LOG(ARY(KL)) + H*(LOG(ARG) - LOG(ARX(KL))))
         IF (MODE .LT. 0 .AND. SS .LT. 0.0) SS = 0.0
!
!  LOGARITHMIC INTERPOLATION
!
       ELSEIF (ABS(MODE) .EQ. 3) THEN
         AA = ARX(KL+1)/ARX(KL)
         BB = ARG/ARX(KL)
         IF (AA .EQ. BB) THEN
           SS = ARY(KL+1)
         ELSE
           H = (ARY(KL+1) - ARY(KL))/LOG(AA)
           SS = ARY(KL) + H*LOG(BB)
         ENDIF
         IF (MODE .LT. 0 .AND. SS .LT. 0.0) SS = 0.0
!
!  LAGRANGIAN INTERPOLATION
!
       ELSEIF (ABS(MODE) .EQ. 6) THEN
!
!  NUMBER OF DATA POINT POINTS TO INTERPOLATE BETWEEN
!
         NO = 4
!
!  FIND DATA POINTS BETWEEN WHICH TO INTERPOLATE
!
         NOO = NO
   20    IF (ARG .LT. ARX(1)) THEN
           NOO = 2
           KS = 1
           KF = 1 + NOO - 1
         ELSEIF (ARG .GT. ARX(NREC)) THEN
           NOO = 2
           KS = NREC - NOO + 1
           KF = NREC
         ELSE
           IF (MOD(NOO,2) .EQ. 0) THEN
             KS = KL - 0.5*NOO + 1
             KF = KL + 0.5*NOO
           ELSE
             KS = KL - 0.5*(NOO - 1) + 1
             KF = KL + 0.5*(NOO + 1)
           ENDIF
           IF (KS .LT. 1) KS = 1
           IF (KF .GT. NREC) KF = NREC
         ENDIF
!
!  INTERPOLATE FUNCTION VALUE AT ARG FROM DATA POINTS KS TO KF
!
         SS = 0.0
         DO 30 K=KS,KF
           SS = SS + XL(MXDIM,KS,KF,K,ARG,ARX)*ARY(K)     &
          /XL(MXDIM,KS,KF,K,ARX(K),ARX)
   30    CONTINUE
         H = NOO
!
!  IF INTERPOLATION HAS OVERSHOT, REDUCE ORDER AND TRY AGAIN
!
         IF (ARG .LT. ARX(1) .OR. ARG .GT. ARX(NREC)) THEN
           IF (SS .LT. 0.0) THEN
             IF (NOO .EQ. 2) THEN
               SS = 0.0
             ELSE
               NOO = NOO - 1
               GOTO 20
             ENDIF
           ENDIF
         ELSE
           IF (((ARY(KL) .LE. ARY(KL+1)) .AND.                       &
          (SS .LT. ARY(KL) .OR. SS .GT. ARY(KL+1))) .OR.            &
          ((ARY(KL) .GT. ARY(KL+1)) .AND.                           &
          (SS .GT. ARY(KL) .OR. SS .LT. ARY(KL+1)))) THEN
             NOO = NOO-1
             GOTO 20
           ENDIF
         ENDIF
!
       ENDIF
!
!-----------------------------------------------------------------------
       RETURN
       END SUBROUTINE
!***********************************************************************
!
!  PROGRAM        XL  FUNCTION
!
!  PURPOSE        TO COMPUTE LAGRANGE INTERPOLATION COEFFICIENTS
!
!  VERSION        3.0   D.P. EDWARDS   01/01/89
!
!  ARGUMENTS      MXDIM   I*4  I/P  ARRAY DIMENSION
!                 KS      I*4  I/P  LOWER LIMIT OF LAGRANGE SUM
!                 KF      I*4  I/P  UPPER LIMIT OF LAGRANGE SUM
!                 K       I*4  I/P  CURRENT INDEX OF LAGRANGE SUM
!                 ARG     R*4  I/P  INTERPOLATION ARGUMENT
!                 ARR     R*4  I/P  ARRAY TO INTERPOLATE BETWEEN
!
!***********************************************************************
!
       FUNCTION XL(MXDIM,KS,KF,K,ARG,ARR)
!-----------------------------------------------------------------------
       use parkind1, only : jprb, jpim
       
       INTEGER(KIND=JPIM) MXDIM
       INTEGER(KIND=JPIM) KS, KF, K
       REAL(KIND=JPRB)    ARG
       REAL(KIND=JPRB)    ARR(MXDIM)
       REAL(KIND=JPRB)    XL
       
       INTEGER(KIND=JPIM) J
       REAL(KIND=JPRB)    PROD
!-----------------------------------------------------------------------
!
       PROD = 1.0
       DO 10 J=KS,KF
         IF (J .NE. K) THEN
           PROD = PROD*(ARG - ARR(J))
         ENDIF
   10  CONTINUE
!
       XL = PROD
!
!-----------------------------------------------------------------------
       RETURN
       END FUNCTION


END MODULE rttov_scattering_mod