
#include "ctrparam.h"

!       ==========================================================
!
!       MD2G04.F:  Lots of utility functions (some of which are unknown).
!
!       ----------------------------------------------------------
!
!       Revision History:
!
!       When    Who             What
!       ----    ----------      -------
!       073100  Chien Wang      repack based on CliChem3 & M24x11,
!                                 and add cpp.
!       082900  Andrei          some subroutines are taken out 
!
!       ==========================================================

      SUBROUTINE CHECKT (N)                                             8742.   
C****                                                                   8743.   
C**** THIS SUBROUTINE CHECKS WHETHER THE TEMPERATURES ARE REASONABLE    8744.   
C**** FOR DEBUGGING PURPOSES. IT IS TURNED ON BY SETTING IDACC(11)      8745.   
C**** TO BE POSITIVE.  REMEMBER TO SET IDACC(11) BACK TO ZERO AFTER     8746.   
C**** THE ERRORS ARE CORRECTED.                                         8747.   
C****                                                                   8748.   

#include "BD2G04.COM"                                                   8749.   

      COMMON U,V,T,P,Q                                                  8750.   
      IF(IDACC(11).LE.0) RETURN                                         8751.   
C****                                                                   8752.   
C**** CHECK WHETHER GDATA ARE REASONABLE AND CONSISTENT OVER EARTH      8753.   
C****                                                                   8754.   
      X=1.001                                                           8755.   
      DO 110 J=1,JM                                                     8756.   
      IMAX=IM                                                           8757.   
      IF((J.EQ.1).OR.(J.EQ.JM)) IMAX=1                                  8758.   
      DO 110 I=1,IMAX                                                   8759.   
      PEARTH=FDATA(I,J,2)*(1.-FDATA(I,J,3))                             8760.   
      IF(PEARTH.LE.0.) GO TO 110                                        8761.   
      IF(GDATA(I,J,2).GE.0..AND.GDATA(I,J,2)*GDATA(I,J,4).LE.0.)GO TO 508762.   
      WRITE (6,901) N,I,J,L,TAU,(GDATA(I,J,K),K=2,10)                   8763.   
      STOP 9021                                                         8764.   
   50 DO  70 L=1,2                                                      8765.   
      TGL=GDATA(I,J,4*L)                                                8766.   
      WTRL=GDATA(I,J,4*L+1)                                             8767.   
      ACEL=GDATA(I,J,4*L+2)                                             8768.   
      IF((TGL+60.)*(60.-TGL).GT.0.) GO TO 60                            8769.   
      WRITE (6,901) N,I,J,L,TAU,(GDATA(I,J,K),K=2,10)                   8770.   
   60 IF(WTRL.GE.0..AND.ACEL.GE.0..AND.TGL*WTRL.GE.0..AND.              8771.   
     *  TGL*ACEL.LE.0..AND.(WTRL+ACEL).LE.X*VDATA(I,J,8+L)) GO TO 70    8772.   
      WRITE (6,901) N,I,J,L,TAU,(GDATA(I,J,K),K=2,10)                   8773.   
      STOP 9021                                                         8774.   
   70 CONTINUE                                                          8775.   
  110 CONTINUE                                                          8776.   
      RETURN                                                            8777.   
C****                                                                   8778.   
  901 FORMAT ('0GDATA UNREASONABLE, N,I,J,L,TAU=',4I4,E14.5/            8779.   
     *  ' GDATA(2-10)=',9F12.4)                                         8780.   
      END                                                               8781.   
c     SUBROUTINE COSZR (IM,JM,SINL,COSL,SIND,COSD,ROT1,ROT2,COSZ,KR)    4795.5 
      SUBROUTINE COSZR (IM,JM,SIND,COSD,ROT1,ROT2,COSZ)
      DIMENSION COSZ(IM,JM),SINL(46),COSL(46)                           4795.51 
      COMMON/SCJL/SINL,COSL
      DATA TWOPI/6.283185/,ZERO/0./                                     4795.52 
      KR=0
      DO 200 J=1,JM                                                     4795.53 
      IF(SINL(J)*SIND-COSL(J)*COSD.GE.-ZERO)                            4795.54 
     *  GO TO 230                                                       4795.55 
      IF(SINL(J)*SIND+COSL(J)*COSD.LE.ZERO) GO TO 250                   4795.56 
      HALFD=ABS(ACOS(-SINL(J)*SIND/(COSL(J)*COSD)))                     4795.57 
      HALFS=HALFD*86400./TWOPI                                          4795.58 
      SINHD=SIN(HALFD)                                                  4795.59 
      TEMP=TWOPI*HALFS/86400.                                           4795.6  
      COSZ(1,J)=SINL(J)*SIND+COSL(J)*COSD*SINHD/TEMP                    4795.61 
      IF(KR.EQ.1) COSZ(1,J)=COSZ(1,J)*2.*HALFS/86400.                   4795.62 
      GO TO 200                                                         4795.63 
 230  COSZ(1,J)=SINL(J)*SIND                                            4795.64 
      GO TO 200                                                         4795.65 
 250  COSZ(1,J)=ZERO                                                    4795.66 
 200  CONTINUE                                                          4795.67 
      RETURN                                                            4795.68 
      END                                                               4795.69 
      SUBROUTINE FILTER                                                 8001.   
C*********                                                              8526.   
C     MFILTR=0  NO SHAPIRO FILTER                                       8527.   
C     MFILTR=1     SHAPIRO FILTER 1D  2TH ORDER ON SEA LEVEL PRESSURE   8528.   
C     MFILTR=2     SHAPIRO FILTER 1D  2TH ORDER ON SLP & THETA          8529.   
C     MFILTR=3     SHAPIRO FILTER 1D  2TH ORDER ON U-WIND               8530.   
C*********                                                              8531.   

#include "BD2G04.COM"                                                    8532.   

      COMMON U,V,T,P,Q                                                  8533.   
      COMMON/WORK2/X1JI(72,46),X2JI(72,46),X3JI(72,46),X1(72),X2(72),   8534.   
     *  X3(72),X4(72)                                                   8535.   
      LOGICAL HPRNT
      DATA IFIRST/1/                                                    8536.   
      IF ( IFIRST.NE.1 ) GO TO 100                                      8537.   
      IFIRST=0                                                          8538.   
C**** INITIALIZE CERTAIN QUANTITIES                                     8539.   
      CUTPLD=0.5                                                        8540.   
      BETAG=.0065/GRAV                                                  8541.   
      GBYRB=1./(RGAS*BETAG)                                             8542.   
      AKAP=KAPA-.205                                                    8543.   
      IMM1=IM-1                                                         8544.   
      IMBY2=IM/2                                                        8545.   
      NORDER=16                                                         8546.   
c     NORDER=8
      XXX=4.**NORDER                                                    8547.   
      WRITE(6,4000)    MFILTR,NORDER,XXX                                8548.   
 4000 FORMAT('  PERFORM SHAPIRO FILTER  MFILTR=',I2,'  NORDER=',I2,     8549.   
     *  '    XXX=',F15.0)                                               8550.   
      print *,' filter for ps, t, and q'
c     print *,' filter for q only'
  100 CONTINUE                                                          8551.   
      HPRNT=.FALSE.
C  
c     go to 250
C 
C**** PERFORM FILTER ON SEA LEVEL PRESSURE                              8552.   
c      IF (MFILTR.NE.1) GO TO 200                                       8553.   
      DO 110 J=1,JM                                                     8554.   
      DO 110 I=1,IM                                                     8555.   
      X2JI(I,J)=(1.+BETAG*FDATA(I,J,1)/BLDATA(I,J,2))**GBYRB            8556.   
      X1JI(I,J)=P(I,J)*X2JI(I,J)                                        8557.   
      X3JI(I,J)=X1JI(I,J)                                               8558.   
  110 CONTINUE                                                          8559. 
      CALL SHAP2D (MFILTR,NORDER,XXX,IM,JM,1,1)                         8560.   
      DO 44 I=1,IM                                                      8561.   
      DO 44 J=1,JM                                                      8562.   
      P(I,J)=X1JI(I,J)/X2JI(I,J)                                        8563.   
   44 CONTINUE                                                          8564.   
      DOPK=1.                                                           8565.   
  200 continue
c 200 IF(MFILTR.NE.2) GO TO 300                                         8566.   
C**** PERFORM FILTER ON POTENTIAL TEMPERATURE                           8567.   
      DO 240 L=1,LM                                                     8568.   
      DO 210 J=1,JM                                                     8569.   
      DO 210 I=1,IM                                                     8570.   
      X2JI(I,J)=1.                                                      8571.   
      X1JI(I,J)=T(I,J,L)*X2JI(I,J)                                      8572.   
      X3JI(I,J)=X1JI(I,J)                                               8573.   
  210 CONTINUE                                                          8574.   
      CALL SHAP2D (MFILTR,NORDER,XXX,IM,JM,1,1)                         8575.   
      DO 240 J=1,JM                                                     8576.   
      DO 240 I=1,IM                                                     8577.   
  240 T(I,J,L)=X1JI(I,J)/X2JI(I,J)                                      8578.   
  250 continue
C**** PERFORM FILTER ON SPECIFIC HUMIDITY                               8578.01 
      DO 280 L=1,LM                                                     8578.02 
      DO 260 J=1,JM                                                     8578.03 
      DO 260 I=1,IM                                                     8578.04 
      X2JI(I,J)=1.                                                      8578.05 
      X1JI(I,J)=Q(I,J,L)*X2JI(I,J)                                      8578.06 
      X3JI(I,J)=X1JI(I,J)                                               8578.07 
  260 CONTINUE                                                          8578.08 
      CALL SHAP2D (MFILTR,NORDER,XXX,IM,JM,1,1)                         8578.09 
      DO 280 J=1,JM                                                     8578.1  
      DO 280 I=1,IM                                                     8578.11 
      QTEMP=X1JI(I,J)/X2JI(I,J)                                         8578.12 
      IF(QTEMP.LE.0.) QTEMP=0.                                          8578.13 
  280 Q(I,J,L)=QTEMP                                                    8578.15 
  300 IF (MFILTR.NE.3) RETURN                                           8579.   
      print *,' WRONG FILTER'
C**** PERFORM FILTER ON U-WIND                                          8580.   
      DO 340 L=1,LM                                                     8581.   
      DO 310 J=2,JM                                                     8582.   
      X1JI(1,J)=U(1,J,L)                                                8583.   
  310 X3JI(1,J)=U(1,J,L)                                                8584.   
      CALL SHAP2D (MFILTR,NORDER,XXX,IM,JM,2,3)                         8585.   
      DO 340 J=2,JM                                                     8586.   
  340 U(1,J,L)=X1JI(1,J)                                                8587.   
C**** PERFORM FILTER ON V-WIND        
      if(HPRNT)then
       print *,' from filter 1'
       print *,' V(J,3)'
       print *,(V(1,J,3),J=1,JM)
      endif
      DO 440 L=1,LM                  
      DO 410 J=2,JM                 
      X1JI(1,J)=V(1,J,L)           
  410 X3JI(1,J)=V(1,J,L)          
      CALL SHAP2D (MFILTR,NORDER,XXX,IM,JM,2,3)     
      DO 440 J=2,JM                               
  440 V(1,J,L)=X1JI(1,J)     
      if(HPRNT)then
       print *,' from filter 2'
       print *,' V(J,3)'
       print *,(V(1,J,3),J=1,JM)
      endif
      RETURN                                                            8588.   
      END                                                               8589.   
      SUBROUTINE DTDX1D(XNOW,XREF,XDT0,SDT0,KFOR,NFOR)                  9112.4  
      DIMENSION XNOW(5),XREF(5),XDT0(5),KFOR(5)                         9112.5  
C     FCO2(X)=DLOG(1.0+1.3365*X+9.29E-03*X*X+2.1766E-06*X**3)           9112.6  
C     FFNC(X,Y)=1.54*DLOG(1.0+1.11*X**0.77) - 0.068*DLOG(1.0+X*Y)       9112.7  
C    +          +0.231*Y**0.65/(1.0+0.113*Y**0.71)                      9112.8  
C                                                                       9112.9  
      FCO2(X)=DLOG(1.0+0.942*X/(1.0+6.2E-04*X)+8.8E-03*X*X              9113.   
     +        +3.26E-06*X**3+0.156*X**1.3*EXP(-X/760.0))                9113.1  
C                                                                       9113.2  
      FFNC(X,Y)=1.556*DLOG(1.0+1.098*X**0.77*(1.0+0.032*X)              9113.3  
     +                                       /(1.0+0.0014*X*X))         9113.4  
     +          +(0.394*Y**0.66+0.16*Y*EXP(-1.6*Y))                     9113.5  
     +                         /(1.0+0.169*Y**0.62)                     9113.6  
     +          -0.14*DLOG(1.0+0.636*(X*Y)**0.75+0.007*Y*(X*Y)**1.52)   9113.7  
C                                                                       9113.8  
      SUM=0.0                                                           9113.9  
      DO 200 K=1,NFOR                                                   9114.   
      XDT0(K)=0.0                                                       9114.1  
      XX=XNOW(K)                                                        9114.2  
      XR=XREF(K)                                                        9114.3  
      IF(KFOR(K).EQ.0) GO TO 200                                        9114.4  
      GO TO (110,120,130,140,150),K                                     9114.5  
      GO TO 200                                                         9114.6  
 110  XDTX=FCO2(XX)                                                     9114.7  
      XDTR=FCO2(XR)                                                     9114.8  
      GO TO 190                                                         9114.9  
 120  YY=XNOW(3)                                                        9115.   
      YR=XREF(3)                                                        9115.1  
      XDTX=FFNC(XX,YR)                                                  9115.2  
      XDTR=FFNC(XR,YR)                                                  9115.3  
      GO TO 190                                                         9115.4  
 130  XR=XREF(2)                                                        9115.5  
      YY=XNOW(3)                                                        9115.6  
      YR=XREF(3)                                                        9115.7  
      XDTX=FFNC(XR,YY)                                                  9115.8  
      XDTR=FFNC(XR,YR)                                                  9115.9  
      GO TO 190                                                         9116.   
 140  XX=XNOW(4)                                                        9116.1  
      XR=XREF(4)                                                        9116.2  
      XDTX=0.066*XX                                                     9116.3  
      XDTR=0.066*XR                                                     9116.4  
      GO TO 190                                                         9116.5  
 150  XX=XNOW(5)                                                        9116.6  
      XR=XREF(5)                                                        9116.7  
      XDTX=0.084*XX                                                     9116.8  
      XDTR=0.084*XR                                                     9116.9  
      GO TO 190                                                         9117.   
 190  XDT0(K)=XDTX-XDTR                                                 9117.1  
      SUM=SUM+XDT0(K)                                                   9117.2  
 200  CONTINUE                                                          9117.3  
      SDT0=SUM                                                          9117.4  
      RETURN                                                            9117.5  
      END                                                               9117.6  
      SUBROUTINE DTDX3D(XNOW,XREF,XDT0,SDT0,KFOR,NFOR)                  9126.5  
      DIMENSION XNOW(5),XREF(5),XDT0(5),KFOR(5)                         9126.6  
C     FCO2(X)=DLOG(1.0 + 0.71*X**1.1 + 0.0084*X**1.9 + 2.8E-07*X**3.3)  9126.7  
      FCO2(X)=DLOG(1.0 + 1.2*X + 0.005*X**2 + 1.4E-06*X**3)             9126.8  
      FFNC(X,Y)=3.1515*X/(1.0+0.21448*X**0.87) - 0.056*X*Y              9126.9  
     +          +0.3234*Y/(1.0+0.052*Y**0.84)                           9127.   
C                                                                       9127.1  
      SUM=0.0                                                           9127.2  
      DO 200 K=1,NFOR                                                   9127.3  
      XDT0(K)=0.0                                                       9127.4  
      IF(KFOR(K).EQ.0) GO TO 200                                        9127.5  
      GO TO (110,120,130,140,150),K                                     9127.6  
      GO TO 200                                                         9127.7  
 110  XX=XNOW(1)                                                        9127.8  
      XR=XREF(1)                                                        9127.9  
      XDTX=FCO2(XX)                                                     9128.   
      XDTR=FCO2(XR)                                                     9128.1  
      GO TO 190                                                         9128.2  
 120  XX=XNOW(2)                                                        9128.3  
      XR=XREF(2)                                                        9128.4  
      YR=XREF(3)                                                        9128.5  
      XDTX=FFNC(XX,YR)                                                  9128.6  
      XDTR=FFNC(XR,YR)                                                  9128.7  
      GO TO 190                                                         9128.8  
 130  XR=XREF(2)                                                        9128.9  
      YY=XNOW(3)                                                        9129.   
      YR=XREF(3)                                                        9129.1  
      XDTX=FFNC(XR,YY)                                                  9129.2  
      XDTR=FFNC(XR,YR)                                                  9129.3  
      GO TO 190                                                         9129.4  
 140  XX=XNOW(4)                                                        9129.5  
      XR=XREF(4)                                                        9129.6  
      XDTX=0.066*XX                                                     9129.7  
      XDTR=0.066*XR                                                     9129.8  
      GO TO 190                                                         9129.9  
 150  XX=XNOW(5)                                                        9130.   
      XR=XREF(5)                                                        9130.1  
      XDTX=0.084*XX                                                     9130.2  
      XDTR=0.084*XR                                                     9130.3  
      GO TO 190                                                         9130.4  
 190  XDT0(K)=XDTX-XDTR                                                 9130.5  
      SUM=SUM+XDT0(K)                                                   9130.6  
 200  CONTINUE                                                          9130.7  
      SDT0=SUM                                                          9130.8  
      RETURN                                                            9130.9  
      END                                                               9131.   
      SUBROUTINE DXDT1D(XNOW,XREF,XDT0,SDT0,KFOR,NFOR)                  9117.7  
      DIMENSION XNOW(5),XREF(5),XDT0(5),KFOR(5)                         9117.8  
C     FCO2(X)=DLOG(1.0+1.3365*X+9.29E-03*X*X+2.1766E-06*X**3)           9117.9  
C     FFNC(X,Y)=1.54*DLOG(1.0+1.11*X**0.77) - 0.068*DLOG(1.0+X*Y)       9118.   
C    +          +0.231*Y**0.65/(1.0+0.113*Y**0.71)                      9118.1  
C                                                                       9118.2  
      FCO2(X)=DLOG(1.0+0.942*X/(1.0+6.2E-04*X)+8.8E-03*X*X              9118.3  
     +        +3.26E-06*X**3+0.156*X**1.3*EXP(-X/760.0))                9118.4  
C                                                                       9118.5  
      FFNC(X,Y)=1.556*DLOG(1.0+1.098*X**0.77*(1.0+0.032*X)              9118.6  
     +                                       /(1.0+0.0014*X*X))         9118.7  
     +          +(0.394*Y**0.66+0.16*Y*EXP(-1.6*Y))                     9118.8  
     +                         /(1.0+0.169*Y**0.62)                     9118.9  
     +          -0.14*DLOG(1.0+0.636*(X*Y)**0.75+0.007*Y*(X*Y)**1.52)   9119.   
C                                                                       9119.1  
      DO 200 K=1,NFOR                                                   9119.2  
      XN=XNOW(K)                                                        9119.3  
      XR=XREF(K)                                                        9119.4  
      ITER=0                                                            9119.5  
      DTK=XDT0(K)                                                       9119.6  
      IF(XN.EQ.XR) GO TO 200                                            9119.7  
      IF(KFOR(K).EQ.0) GO TO 200                                        9119.8  
      GO TO (110,120,130,140,150),K                                     9119.9  
      GO TO 200                                                         9120.   
 110  TR=FCO2(XR)                                                       9120.1  
 111  TN=FCO2(XN)                                                       9120.2  
      DTNR=TN-TR                                                        9120.3  
      ITER=ITER+1                                                       9120.4  
      SLOPE=DTNR/(XN-XR)                                                9120.5  
      DT=DTK-DTNR                                                       9120.6  
      DX=DT/SLOPE                                                       9120.7  
      XN=XN+DX                                                          9120.8  
      IF(ITER.LT.15.AND.ABS(DT).GT.1.E-05) GO TO 111                    9120.9  
      GO TO 190                                                         9121.   
  120 YR=XREF(3)                                                        9121.1  
      TR=FFNC(XR,YR)                                                    9121.2  
  121 TN=FFNC(XN,YR)                                                    9121.3  
      DTNR=TN-TR                                                        9121.4  
      ITER=ITER+1                                                       9121.5  
      SLOPE=DTNR/(XN-XR)                                                9121.6  
      DT=DTK-DTNR                                                       9121.7  
      DX=DT/SLOPE                                                       9121.8  
      XN=XN+DX                                                          9121.9  
      IF(ITER.LT.15.AND.ABS(DT).GT.1.E-05) GO TO 121                    9122.   
      GO TO 190                                                         9122.1  
 130  XR=XREF(2)                                                        9122.2  
      YN=XNOW(3)                                                        9122.3  
      YR=XREF(3)                                                        9122.4  
      TR=FFNC(XR,YR)                                                    9122.5  
  131 TN=FFNC(XR,YN)                                                    9122.6  
      DTNR=TN-TR                                                        9122.7  
      ITER=ITER+1                                                       9122.8  
      SLOPE=DTNR/(YN-YR)                                                9122.9  
      DT=DTK-DTNR                                                       9123.   
      DY=DT/SLOPE                                                       9123.1  
      YN=YN+DY                                                          9123.2  
      IF(ITER.LT.15.AND.ABS(DT).GT.1.E-05) GO TO 131                    9123.3  
      XN=YN                                                             9123.4  
      GO TO 190                                                         9123.5  
 140  XN=XNOW(4)                                                        9123.6  
      XR=XREF(4)                                                        9123.7  
  141 TN=0.066*XN                                                       9123.8  
      TR=0.066*XR                                                       9123.9  
      DTNR=TN-TR                                                        9124.   
      ITER=ITER+1                                                       9124.1  
      SLOPE=DTNR/(XN-XR)                                                9124.2  
      DT=DTK-DTNR                                                       9124.3  
      DX=DT/SLOPE                                                       9124.4  
      XN=XN+DX                                                          9124.5  
      IF(ITER.LT.15.AND.ABS(DT).GT.1.E-05) GO TO 141                    9124.6  
      TR=0.066*XR                                                       9124.7  
      GO TO 190                                                         9124.8  
 150  XN=XNOW(5)                                                        9124.9  
      XR=XREF(5)                                                        9125.   
      TR=0.084*XR                                                       9125.1  
  151 TN=0.084*XN                                                       9125.2  
      DTNR=TN-TR                                                        9125.3  
      ITER=ITER+1                                                       9125.4  
      SLOPE=DTNR/(XN-XR)                                                9125.5  
      DT=DTK-DTNR                                                       9125.6  
      DX=DT/SLOPE                                                       9125.7  
      XN=XN+DX                                                          9125.8  
      IF(ITER.LT.15.AND.ABS(DT).GT.1.E-05) GO TO 151                    9125.9  
      GO TO 190                                                         9126.   
 190  XNOW(K)=XN                                                        9126.1  
 200  CONTINUE                                                          9126.2  
      RETURN                                                            9126.3  
      END                                                               9126.4  
      SUBROUTINE DXDT3D(XNOW,XREF,XDT0,SDT0,KFOR,NFOR)                  9131.1  
      DIMENSION XNOW(5),XREF(5),XDT0(5),KFOR(5)                         9131.2  
C     FCO2(X)=DLOG(1.0 + 0.71*X**1.1 + 0.0084*X**1.9 + 2.8E-07*X**3.3)  9131.3  
      FCO2(X)=DLOG(1.0 + 1.2*X + 0.005*X**2 + 1.4E-06*X**3)             9131.4  
      FFNC(X,Y)=3.1515*X/(1.0+0.21448*X**0.87) - 0.056*X*Y              9131.5  
     +          +0.3234*Y/(1.0+0.052*Y**0.84)                           9131.6  
C                                                                       9131.7  
      DO 200 K=1,NFOR                                                   9131.8  
      XN=XNOW(K)                                                        9131.9  
      XR=XREF(K)                                                        9132.   
      ITER=0                                                            9132.1  
      DTK=XDT0(K)                                                       9132.2  
      IF(XN.EQ.XR) GO TO 200                                            9132.3  
      IF(KFOR(K).EQ.0) GO TO 200                                        9132.4  
      GO TO (110,120,130,140,150),K                                     9132.5  
      GO TO 200                                                         9132.6  
 110  TR=FCO2(XR)                                                       9132.7  
 111  TN=FCO2(XN)                                                       9132.8  
      DTNR=TN-TR                                                        9132.9  
      ITER=ITER+1                                                       9133.   
      SLOPE=DTNR/(XN-XR)                                                9133.1  
      DT=DTK-DTNR                                                       9133.2  
      DX=DT/SLOPE                                                       9133.3  
      XN=XN+DX                                                          9133.4  
      IF(ITER.LT.15.AND.ABS(DT).GT.1.E-05) GO TO 111                    9133.5  
      GO TO 190                                                         9133.6  
  120 YR=XREF(3)                                                        9133.7  
      TR=FFNC(XR,YR)                                                    9133.8  
  121 TN=FFNC(XN,YR)                                                    9133.9  
      DTNR=TN-TR                                                        9134.   
      ITER=ITER+1                                                       9134.1  
      SLOPE=DTNR/(XN-XR)                                                9134.2  
      DT=DTK-DTNR                                                       9134.3  
      DX=DT/SLOPE                                                       9134.4  
      XN=XN+DX                                                          9134.5  
      IF(ITER.LT.15.AND.ABS(DT).GT.1.E-05) GO TO 121                    9134.6  
      GO TO 190                                                         9134.7  
 130  XR=XREF(2)                                                        9134.8  
      YN=XNOW(3)                                                        9134.9  
      YR=XREF(3)                                                        9135.   
      TR=FFNC(XR,YR)                                                    9135.1  
  131 TN=FFNC(XR,YN)                                                    9135.2  
      DTNR=TN-TR                                                        9135.3  
      ITER=ITER+1                                                       9135.4  
      SLOPE=DTNR/(YN-YR)                                                9135.5  
      DT=DTK-DTNR                                                       9135.6  
      DY=DT/SLOPE                                                       9135.7  
      YN=YN+DY                                                          9135.8  
      IF(ITER.LT.15.AND.ABS(DT).GT.1.E-05) GO TO 131                    9135.9  
      XN=YN                                                             9136.   
      GO TO 190                                                         9136.1  
 140  XN=XNOW(4)                                                        9136.2  
      XR=XREF(4)                                                        9136.3  
  141 TN=0.066*XN                                                       9136.4  
      TR=0.066*XR                                                       9136.5  
      DTNR=TN-TR                                                        9136.6  
      ITER=ITER+1                                                       9136.7  
      SLOPE=DTNR/(XN-XR)                                                9136.8  
      DT=DTK-DTNR                                                       9136.9  
      DX=DT/SLOPE                                                       9137.   
      XN=XN+DX                                                          9137.1  
      IF(ITER.LT.15.AND.ABS(DT).GT.1.E-05) GO TO 141                    9137.2  
      TR=0.066*XR                                                       9137.3  
      GO TO 190                                                         9137.4  
 150  XN=XNOW(5)                                                        9137.5  
      XR=XREF(5)                                                        9137.6  
      TR=0.084*XR                                                       9137.7  
  151 TN=0.084*XN                                                       9137.8  
      DTNR=TN-TR                                                        9137.9  
      ITER=ITER+1                                                       9138.   
      SLOPE=DTNR/(XN-XR)                                                9138.1  
      DT=DTK-DTNR                                                       9138.2  
      DX=DT/SLOPE                                                       9138.3  
      XN=XN+DX                                                          9138.4  
      IF(ITER.LT.15.AND.ABS(DT).GT.1.E-05) GO TO 151                    9138.5  
      GO TO 190                                                         9138.6  
 190  XNOW(K)=XN                                                        9138.7  
 200  CONTINUE                                                          9138.8  
      RETURN                                                            9138.9  
      END                                                               9139.   

      FUNCTION FUNC1(KX,KY,KZ)                                          2666.5  
      CONST=0.5                                                         2667.   
      KA=KX-1                                                           2667.5  
      KB=KY-1                                                           2668.   
      KC=KZ-1                                                           2668.5  
      K1=KA+KB                                                          2669.   
      K2=KA+KC                                                          2669.5  
      K3=KB+KC                                                          2670.   
      FUNC1=0.                                                          2670.5  
      IF((KA.EQ.K3).OR.(KB.EQ.K2).OR.(KC.EQ.K1)) FUNC1=CONST            2671.   
      RETURN                                                            2671.5  
C                                                                       2672.   
      ENTRY FUNC2(KX,KY,KZ)                                             2672.5  
      CONST=0.5                                                         2673.   
      KA=KX-1                                                           2673.5  
      KB=KY-1                                                           2674.   
      KC=KZ-1                                                           2674.5  
      K1=KA+KB                                                          2675.   
      K2=KA+KC                                                          2675.5  
      K3=KB+KC                                                          2676.   
      FUNC2=0.                                                          2676.5  
      IF((KA.GT.KB).AND.(KA.GT.KC).AND.(KA.EQ.K3)) FUNC2 =-CONST        2677.   
      IF((KA.LT.KB).AND.(KB.EQ.K2)) FUNC2=CONST                         2677.5  
      IF((KA.LT.KC).AND.(KC.EQ.K1)) FUNC2=CONST                         2678.   
      RETURN                                                            2678.5  
      END                                                               2679.   

c     FUNCTION RANDU (X)                                                3801.   
      SUBROUTINE RANDUU(RX,X)
      integer IX,IY
C**** THIS FUNCTION GENERATES RANDOM NUMBERS ON AN IBM 360 OR 370       3802.   
   10 IY=IX*65539                                                       3803.   
      IF(IY) 20,40,30                                                   3804.   
   20 IY=(IY+2147483647)+1                                              3805.   
   30 IX=IY                                                             3806.   
c     RANDU=FLOAT(IY)*.4656613E-9                                       3807.   
      RX=FLOAT(IY)*.4656613E-9
      RETURN                                                            3808.   
   40 IX=1                                                              3809.   
      GO TO 10                                                          3810.   
      ENTRY RINIT (INIT)                                                3811.   
      IX=INIT                                                           3812.   
      RETURN                                                            3813.   
      ENTRY RFINAL (IFINAL)                                             3814.   
      IFINAL=IX                                                         3815.   
      RETURN                                                            3816.   
      END                                                               3817.   
      SUBROUTINE SHAP2D (MFILTR,NORDER,XXX,IM,JM,J1,ITYPE)              8590.   
      COMMON/WORK2/X1JI(72,46),X2JI(72,46),X3JI(72,46),X1(72),X2(72),   8591.   
     *  X3(72),X4(72),XM1(72),XJMP1(72)                                 8592.   
C      VARIABLE ITYPE DETERMINES TYPE OF BOUNDARY CONDITIONS
C      ITYPE=1 FOR PS,T AND Q ( XM1=X2)
C      ITYPE=2 FOR U (XM1=X1)
C      ITYPE=3 FOR V (XM1=-X1)
      JMM1=JM-1                                                         8593.   
      J2=J1+1                                                           8594.   
      IMBY2=1                                                           8595.   
      DO 145 N=1,NORDER                                                 8596.   
      DO 146 K=1,IM                                                     8597.   
      X1(K)=X1JI(K,J1)                                                  8598.   
      X2(K)=X1JI(K,J2)                                                  8599.   
      X3(K)=X1JI(K,JMM1)                                                8600.   
      X4(K)=X1JI(K,JM)                                                  8601.   
      IF(ITYPE.EQ.1)THEN
        XM1(K)=X1JI(K,J2)
        XJMP1(K)=X1JI(K,JMM1)
      ELSEIF(ITYPE.EQ.2)THEN
        XM1(K)=X1JI(K,J1)
        XJMP1(K)=X1JI(K,JM)
      ELSE
        XM1(K)=-X1JI(K,J1)
        XJMP1(K)=-X1JI(K,JM)
      ENDIF
  146 CONTINUE                                                          8602.   
      DO 142 I=1,IM                                                     8603.   
      X1IM1=X1JI(I,J1)                                                  8604.   
      DO 142 J=J2,JMM1                                                  8605.   
      X1I=X1JI(I,J)                                                     8606.   
      X1JI(I,J)=X1IM1-X1I-X1I+X1JI(I,J+1)                               8607.   
      X1IM1=X1I                                                         8608.   
 142  CONTINUE                                                          8609.   
      SUM1=0.                                                           8610.   
      SUMJM=0.                                                          8611.   
      DO 144 K=1,IMBY2                                                  8612.   
ccc   SUM1 =SUM1 +X2(K)-X1(K)-X1(K)+X2(K)                               8613.   
      SUM1 =SUM1 +XM1(K)-X1(K)-X1(K)+X2(K)
ccc   SUMJM=SUMJM+X3(K)-X4(K)-X4(K)+X3(K)                               8614.   
      SUMJM=SUMJM+X3(K)-X4(K)-X4(K)+XJMP1(K)
  144 CONTINUE                                                          8615.   
      X1SUM     =SUM1 /IMBY2                                            8616.   
      XJMSUM    =SUMJM/IMBY2                                            8617.   
      DO 147 K=1,IM                                                     8618.   
      X1JI(K,JM)=XJMSUM                                                 8619.   
  147 X1JI(K,J1)= X1SUM                                                 8620.   
 145  CONTINUE                                                          8621.   
      DO 160 I=1,IM                                                     8622.   
      DO 160 J=J1,JM                                                    8623.   
      X1JI(I,J)=(X3JI(I,J)-X1JI(I,J)/XXX)                               8624.   
  160 CONTINUE                                                          8625.   
      RETURN                                                            8626.   
      END                                                               8627.   
