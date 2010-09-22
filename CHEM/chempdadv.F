
#include "ctrparam.h"

!	============================================================
!
!	PDADV.F:   Subroutines of Modified Bott advection scheme
!
!	------------------------------------------------------------
!
!	Author:		Chien Wang
!			MIT Joint Program on Science and Policy
!				of Global Change
!
!	----------------------------------------------------------
! 
!	Revision History:
!	
!	When	Who		What
!	----	----------	-------
!	080200	Chien Wang	repack based on CliChem3 & add cpp
!
!	==========================================================

C     **************************************
C     **************************************
      SUBROUTINE pdadv1(C,W4,W2,W1,N)
C     **************************************
C     **************************************
C
C     ******************************************************************
C     
C    This is a subroutine for the first part of Bott's advection scheme.
C
C     Andreas Bott  1989:  A Positive Definite Advection scheme obtained
C                   by Nonlinear Renormalization of the advective fluxes
C                                             Mon. Wea. Rev. 117 1006-15
C
C     Fourth Order: with coefficients from    Mon. Wea. Rev. 117 2633-36
C
C     Input: C=U*DT/DX[N+1]     Output: W4[3:N1,5],W2[2;3;n1;n,3] and 
C                                       W1[1;2;n;n+1,2]
C     On the Staggered Grid:      C(i')----Q(i)----C(i'+1)
C
C     ******************************************************************

      PARAMETER ( C0=1.0/1920.0,C1=1.0/384.00,C2=1.0/384.0
     &          , C3=1.0/768.00,C4=1.0/3840.0,EP=1.0E-15 )
c      parameter (cc0=1.,cc1=1./16.,cc2=1./48.)
      parameter (cc0=-1./24.,cc1=1./16.,cc2=1./48.)
c      parameter (cc0=-1./24.,cc1=1./16.,cc2=1./16.)
 
      DIMENSION   C(N+1),W4(N,5),W2(N,3),W1(4,2)

!	-----------------------------------------------------------

#if ( defined CPL_CHEM )

      n1=n-1
      n2=n-2
      n3=n-3

      do 1 i=1,n
      do 2 j=1,5
        w4(i,j)=0.0
2     continue
      do 3 j=1,3
        w2(i,j)=0.0
3     continue
1     continue

C
C     GET THE COEFFICIENTS DEPENDENT ON C ONLY
C 
      w1(1,1)=abs(c(1))
      w1(1,2)=0.0
      w1(2,1)=abs(c(2))
      w1(2,2)=2.0*w1(2,1)*(1.-w1(2,1))
      w1(3,1)=abs(c(n))
      w1(3,2)=2.0*w1(3,1)*(1.-w1(3,1))
      w1(4,1)=abs(c(n+1))
      w1(4,2)=0.0

      rr1=abs(c(2))
      rr2=1.-(rr1+rr1)
      r1=rr2**2
      r2=r1*rr2
      w2(2,1)=rr1*cc0
      w2(2,2)=(1.-r1)*cc1
      w2(2,3)=(1.-r2)*cc2

      rr1=abs(c(3))
      rr2=1.-(rr1+rr1)
      r1=rr2**2
      r2=r1*rr2
      w2(3,1)=rr1*cc0
      w2(3,2)=(1.-r1)*cc1
      w2(3,3)=(1.-r2)*cc2

      rr1=abs(c(4))
      rr2=1.-(rr1+rr1)
      r1=rr2**2
      r2=r1*rr2
      w2(4,1)=rr1*cc0
      w2(4,2)=(1.-r1)*cc1
      w2(4,3)=(1.-r2)*cc2

      rr1=abs(c(n2))
      rr2=1.-(rr1+rr1)
      r1=rr2**2
      r2=r1*rr2
      w2(n2,1)=rr1*cc0
      w2(n2,2)=(1.-r1)*cc1
      w2(n2,3)=(1.-r2)*cc2

      rr1=abs(c(n1))
      rr2=1.-(rr1+rr1)
      r1=rr2**2
      r2=r1*rr2
      w2(n1,1)=rr1*cc0
      w2(n1,2)=(1.-r1)*cc1
      w2(n1,3)=(1.-r2)*cc2

      rr1=abs(c(n))
      rr2=1.-(rr1+rr1)
      r1=rr2**2
      r2=r1*rr2
      w2(n,1)=rr1*cc0
      w2(n,2)=(1.-r1)*cc1
      w2(n,3)=(1.-r2)*cc2

      DO 100 I = 3 ,N1

        rr1  = ABS( C(I) )
        rr2  = 1.0 - (rr1+rr1)
        R1       = Rr2*Rr2
        R2       = R1*Rr2
        R3       = R2*Rr2
        R4       = R3*Rr2
 
        W4(I,1)  = rr1 *C0
        W4(I,2)  = (1.0-R1)*C1
        W4(I,3)  = (1.0-R2)*C2
        W4(I,4)  = (1.0-R3)*C3
        W4(I,5)  = (1.0-R4)*C4

100   CONTINUE
C

#endif

      return
      end


C     **************************************
C     **************************************
      SUBROUTINE pdadv2(C,Q,W4,W2,W1,ww,ww2,N,NOOS)
C     **************************************
C     **************************************
C
C     *************************************************************
C
C     This is a subroutine for the second part of Bott's advection
C     scheme.
C
C     Andreas Bott  1989:  A Positive Definite Advection scheme obtained
C                   by Nonlinear Renormalization of the advective fluxes
C                                             Mon. Wea. Rev. 117 1006-15
C
C     Fourth Order: with coefficients from    Mon. Wea. Rev. 117 2633-36
C
C     Input: C=U*DT/DX[N+1]  &  Q[N]      Output: Q[2 N-1]
C     On the Staggered Grid:      C(i')----Q(i)----C(i'+1)
C
C     NOSS     = 1: Perform non-oscillatory option
C
      PARAMETER ( C0=1.0/1920.0,C1=1.0/384.00,C2=1.0/384.0
     &          , C3=1.0/768.00,C4=1.0/3840.0,EP=1.0E-15 )
c      parameter ( cc0=1.,cc1=1./16.,cc2=1./24.)
      parameter ( cc0=-1./24.,cc1=1./16.,cc2=1./24.)
c      parameter ( cc0=-1./24.,cc1=1./16.,cc2=1./16.)

      DIMENSION   C(N+1),Q(N),W4(n,5),w2(n,3),w1(4,2),
     &            ww(n+1,5),ww2(n+1,5)
C

!	--------------------------------------------------------

#if ( defined CPL_CHEM )

      N1       = N-1
      N2       = N-2
      N3       = N-3

      do 1 i=1,(n+1)*5
        ww (i,1)=0.0
        ww2(i,1)=0.0
1     continue

C
C     FOR ANY POSITIVE-DEFINITE Q ADVECTION
C
C   1. First order scheme for i=2 and n:

      a0=q(1)
      a1=q(2)-q(1)
      ww(1,1)=a0
      ww(1,2)=a0*w1(1,1)
      ww(2,3)=a0*w1(2,1)+a1*w1(2,2)

      a0=q(n)
      a1=q(n)-q(n1)
      ww(n,1)=a0
      ww(n,2)=a0*w1(3,1)-a1*w1(3,2)
      ww(n+1,3)=a0*w1(4,1)

C   2. Second order scheme for i=2,3,n1,n:

      ww2(1,1)=ww(1,1)
      ww2(1,2)=ww(1,2)
      ww2(2,3)=ww(2,3)

      a0=q(3)-26.*q(2)+q(1)
      a1=q(3)-q(1)
      a2=q(3)-2.*q(2)+q(1)
      ww2(2,1)=cc0*a0+cc2*a2
      ww2(2,2)=a0*w2(2,1)-a1*w2(2,2)+a2*w2(2,3)
      ww2(3,3)=a0*w2(3,1)+a1*w2(3,2)+a2*w2(3,3) 

      a0=q(4)-26.*q(3)+q(2)
      a1=q(4)-q(2)
      a2=q(4)-2.*q(3)+q(2)
      ww2(3,1)=cc0*a0+cc2*a2
      ww2(3,2)=a0*w2(3,1)-a1*w2(3,2)+a2*w2(3,3)
      ww2(4,3)=a0*w2(4,1)+a1*w2(4,2)+a2*w2(4,3)
       
      a0=q(n1)-26.*q(n2)+q(n3)
      a1=q(n1)-q(n3)
      a2=q(n1)-2.0*q(n2)+q(n3)
      ww2(n2,1)=cc0*a0+cc2*a2
      ww2(n2,2)=a0*w2(n2,1)-a1*w2(n2,2)+a2*w2(n2,3)
      ww2(n1,3)=a0*w2(n1,1)+a1*w2(n1,2)+a2*w2(n1,3)

      a0=q(n)-26.*q(n1)+q(n2)
      a1=q(n)-q(n2)
      a2=q(n)-2.*q(n1)+q(n2)
      ww2(n1,1)=cc0*a0+cc2*a2
      ww2(n1,2)=a0*w2(n1,1)-a1*w2(n1,2)+a2*w2(n1,3)
      ww2(n,3) =a0*w2( n,1)+a1*w2( n,2)+a2*w2( n,3)

      ww2(n,1)  =ww(n,1)
      ww2(n,2)  =ww(n,2)
      ww2(n+1,3)=ww(n+1,3)

C   3. Fourth order scheme for i=3,n1:

      ww(2,1)=ww2(2,1)
      ww(2,2)=ww2(2,2)
      ww(3,3)=ww2(3,3)

      ww(n1,1)=ww2(n1,1)
      ww(n1,2)=ww2(n1,2)
      ww(n, 3)=ww2(n, 3)

      DO 200 I = 3 ,N2
      QL2      = Q(I-2)
      QL1      = Q(I-1)
      Q00      = Q(I)
      QR1      = Q(I+1)
      QR2      = Q(I+2)
      QP1      = QR1+QL1
      QP2      = QR2+QL2
      QM1      = QR1-QL1
      QM2      = QR2-QL2
C     COEFFICIENTS: AREA PRESERVING FLUX FORM
      A0       = 9.0*QP2 - 116.0*QP1 + 2134.0*Q00
      A1       =-5.0*QM2 +  34.0*QM1
      A2       =    -QP2 +  12.0*QP1 -   22.0*Q00
      A3       =     QM2 -   2.0*QM1
      A4       =     QP2 -   4.0*QP1 +    6.0*Q00
C     INTEGRALS: FOR THE USE OF IN/OUT FLUX OF THE GRID
      ww(I,1)   = C0*(A0+10.0*A2+A4)
c      ww(I,1)   = Q00
      ww(I,2) = A0*W4(I,1)-A1*W4(I,2)+A2*W4(I,3)
     &         - A3*W4(I,4)+A4*W4(I,5)
      ww(I+1,3)   = A0*W4(I+1,1)+A1*W4(I+1,2)+A2*W4(I+1,3)
     &         +A3*W4(I+1,4)+A4*W4(I+1,5)
200   CONTINUE
C
C     RESTRICT THE INTEGRALS TO PRESERVE THE SIGN
C
      I        = 1
      IF( C(I).GT.0.0 ) THEN
        ww(I,2) = 0.0
      ELSE IF( C(I).LT.0.0 ) THEN
        ww(I,2) = max( 0.0 , ww(I,2) )
      ENDIF
      DO 210 I = 2 ,N
      IF( C(I).GT.0.0 ) THEN
        ww(I,2) = 0.0
        ww(I,3) = max( 0.0 , ww(I,3) )
        ww2(i,2)= 0.0
        ww2(i,3)= max( 0.0, ww2(i,3))
      ELSE IF( C(I).LT.0.0 ) THEN
        ww(I,2) = max( 0.0 , ww(I,2) )
        ww(I,3) = 0.0
        ww2(i,2)= max( 0.0, ww2(i,2) )
        ww2(i,3)= 0.0
      ENDIF
210   CONTINUE
      I        = N+1
      IF( C(I).GT.0.0 ) THEN
        ww(I,3) = max( 0.0 , ww(I,3) )
      ELSE IF( C(I).LT.0.0 ) THEN
        ww(I,3) = 0.0
      ENDIF
      DO 220 I = 1 ,N
        ww(I,1)   = max( ww(I,2)+ww(I+1,3)+EP , ww(I,1) )
        ww2(i,1)  = max(ww2(i,2)+ww2(i+1,3)+ep,ww2(i,1))
220   CONTINUE
C
C     GET THE WEIGHTING FACTOR
C
      DO 230 I = 1 ,N
        ww(I,1)   = Q(I) / ww(I,1)
        ww2(i,1)  = q(i) /ww2(i,1)
230   CONTINUE
C                                         <= ww(I,2)
C     GET THE IN/OUT FLUX OF THE GRID  I --- I+1/2
C                                            ww(I,3) =>
      DO 250 I = 1 ,N+1
        if(i.ne.n+1) ww(I,2)   = ww(I,2)*ww(I,1)
        if(i.ne.1)   ww(I,3)   = ww(I,3)*ww(I-1,1)
        if(i.ne.n+1) ww2(i,2)  = ww2(i,2)*ww2(i,1)
        if(i.ne.1)   ww2(i,3)  = ww2(i,3)*ww2(i-1,1)
250   CONTINUE
C
      IF( NOOS.NE.1 ) THEN
C     COMPUTE THE TOTAL ADVECTION TENDENCY

c      DO 300 I = 2 ,N1
      q(2) =ww2(3,2)-ww2(3,3)-ww2(2,2) +ww2(2,3) 
      q(n1)=ww2(n,2)-ww2(n,3)-ww2(n1,2)+ww2(n1,3)
      DO 300 I = 3 ,N2
c        q(i) = ww(i+1,2)-ww(i+1,3)-ww(i,2)+ww(i,3)       !tendency
        q(i) = ww(i+1,2)-ww(i+1,3)-ww(i,2)+ww(i,3)+q(i) !value
300   CONTINUE
C
      ELSE
C
C     NON-OSCILLATORY OPTION: FCT LIMITER
C     P.K.Smolarkiewicz  &  W.W.Grabowski, 1990: The multidimensional
C     positive definite advection transport algorithm: Nonoscillatory
C     option,  J. Comput. Phys., 86, 355-375
C
C     GET THE DONOR-CELL FLUXES (Low-order)

      DO 400 I = 2 ,N
      IF( C(I).GT.0.0 ) THEN
        ww(I,1) = Q(I-1)
      ELSE
        ww(I,1) =-Q(I)
      ENDIF
400   CONTINUE

c      ww(1,1)=max(-q(1)*c(1),0.0)
      ww(1,1)=abs(q(1)*c(1))
      if(c(1).gt.0.0)then
        ww(1,4)=0.0
        ww(1,5)=ww(1,1)
      else
        ww(1,4)=ww(1,1)
        ww(1,5)=0.0
      endif

      DO 405 I = 2 ,N
        ww(I,1)   = ww(I,1) * C(I)
        ww(I,4)   = 0.0
        ww(I,5)  = 0.0
405   CONTINUE

c      ww(n+1,1)=max(q(n)*c(n+1),0.0)
      ww(n+1,1)=abs(q(n)*c(n+1))
      if(c(n+1).gt.0.0)then
        ww(n+1,4)=0.0
        ww(n+1,5)=ww(n+1,1)
      else
        ww(n+1,4)=ww(n+1,1)
        ww(n+1,5)=0.0
      endif

      DO 410 I = 2 ,N
      IF( C(I).GT.0.0 ) THEN
        ww(I,5)= ww(I,1)
      ELSE
        ww(I,4) = ww(I,1)
      ENDIF
410   CONTINUE

      DO 415 I = 1 ,N
        ww(I,1)   = ww(I+1,4) - ww(I+1,5) - ww(I,4) + ww(I,5)
415   CONTINUE

      DO 420 I = 1 ,N
        ww(I,1)   = ww(I,1) + Q(I)
c        ww(I,1)   = ww(I,1)
420   CONTINUE

C     GET THE A-FLUX = F(High-order)-F(Low-order)
      DO 430 I = 1 ,N
        ww(I,4)   = ww(I,2) - ww(I,4)                     
        ww(I,5)   = ww(I,3) - ww(I,5)                    
430   CONTINUE
      DO 435 I = 1 ,N
        ww(I,2)   = max( 0.0,ww(I,4) ) - min(0.0, ww(I,5) )
        ww(I,3)   = max( 0.0,ww(I,5) ) - min(0.0, ww(I,4) )
435   CONTINUE

      ww(1,4)=min(ww(1,1),ww(2,1),q(1),q(2))
      ww(1,5)=max(ww(1,1),ww(2,1),q(1),q(2))
      DO 440 I = 2 ,N1
      J        = I-1
      K        = I+1
      ww(I,4) = min(ww(J,1),ww(I,1),ww(K,1),Q(J),Q(I),Q(K)) 
      ww(I,5) = max(ww(J,1),ww(I,1),ww(K,1),Q(J),Q(I),Q(K))
440   CONTINUE
      ww(n,4)=min(ww(n1,1),ww(n,1),q(n1),q(n))
      ww(n,5)=max(ww(n1,1),ww(n,1),q(n1),q(n))

      DO 450 I = 1 ,N
      ww(I,4)  =(ww(I,1)-ww(I,4)) / (ww(I,2)+ww(I+1,3)+EP) 
      ww(I,5)  =(ww(I,5)-ww(I,1)) / (ww(I,3)+ww(I+1,2)+EP)  
      Q(I)     = ww(I,01)
450   CONTINUE

      DO 460 I = 2 ,N
        ww(I,1)   = min( 1.0,ww(I-1,5),ww(I,4) )
460   CONTINUE
      DO 465 I = 2 ,N
        ww(I,2)   = ww(I,2) * ww(I,1)                     
465   CONTINUE
      DO 470 I = 2 ,N
        ww(I,1)   = min( 1.0,ww(I-1,4),ww(I,5) )
470   CONTINUE
      DO 475 I = 2 ,N
        ww(I,3)   = ww(I,3) * ww(I,1)
475   CONTINUE
C     COMPUTE THE HIGH-ORDER ADVECTION TENDENCY
      DO 500 I = 2 ,N1
        ww(I,1)   = ww(I+1,2)-ww(I+1,3)-ww(I,2)+ww(I,3)
500   CONTINUE
C
C     COMPUTE THE TOTAL ADVECTION TENDENCY
C
      DO 600 I = 2 ,N1
c        q(i)     = ww(i,1)           !tendency
        q(i)     = ww(i,1)+q(i)     !value
600   CONTINUE

      ENDIF
 
#endif

      RETURN
      END

