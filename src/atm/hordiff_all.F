
#include "ctrparam.h"

!       ==========================================================
!
!       HORDIFF.F:  subroutine for calculating horizontal
!                       diffusion of Q.
!
!       ----------------------------------------------------------
!
!       Revision History:
!
!       When    Who             What
!       ----    ----------      -------
!       080100  Chien Wang      repack based on CliChem3 & M24x11,
!                                 and add cpp.
!
!       ==========================================================

      SUBROUTINE HORDIFFALL(DTDIF)
#include "BD2G04.COM"
      COMMON U,V,T,P,Q  
      DIMENSION VT(IM0,JM0,LM0),TT(IM0,JM0,LM0),PT(IM0,JM0), 
     &  QT(IM0,JM0,LM0),PU(IM0,JM0),FD(IM0,JM0),DQDY(JM0,LM0)
     &  ,DTDY(JM0,LM0),DUDY(JM0,LM0),DVDY(JM0,LM0)
     & ,UT(IM0,JM0,LM0)
      COMMON/HDFLUX/VQHD(JM0,LM0),VTHD(JM0,LM0),VUHD(JM0,LM0),
     &             VVHD(JM0,LM0)
      logical first
      data first /.true./
      I=1
      JMM1=JM0-1
      FDIFU=5.0E5
      FDIFU=2.50E5

      FDIFF=5.0E5

      FDIFFQ=5.0E5
      FDIFFQ=1.00E6

      if(first)then
        print *,' HOR DIFF for Q and T'
        print *,' FDIFF=',FDIFF,' FDIFFQ=',FDIFFQ
        print *,' HOR DIFF for U and V'
        print *,' FDIFU=',FDIFU
        print *,'IM0=',IM0,' JM0=',JM0,' LM0=',LM0
        print *,'IM=',IM,' JM=',JM,' LM=',LM
        first=.false.
      endif
      DO 50 J=1,JM0                                            
 50   FD(I,J)=P(I,J)*DXYP(J)    
      DO 57 L=1,LM0         
      DO 57 J=1,JM0            
       AJL(J,L,56)=AJL(J,L,56)-Q (I,J,L)*P(I,J)
      TT (I,J,L)=T (I,J,L)*FD(I,J) 
 57   QT (I,J,L)=Q (I,J,L)*FD(I,J) 
      DO 100 L=1,LM0
      DO 100 J=2,JM0
        DQDY(J,L)=(Q (1,J,L)-Q (1,J-1,L))/DYV(J)
        DTDY(J,L)=(T (1,J,L)-T (1,J-1,L))/DYV(J)
  100 CONTINUE
      DO 200 L=1,LM0
       PSAV=0.5*(P(1,1)+P(1,2))
       !CAS FLUXQL=FDIFF*DQDY(2,L)*DXV(2)*PSAV*DTDIF
       FLUXQL=FDIFFQ*DQDY(2,L)*DXV(2)*PSAV*DTDIF
       FLUXQL=DMAX1( -0.5*QT(1,2,L), DMIN1(0.5*QT(1,1,L),FLUXQL))
       FLUXTL=FDIFF*DTDY(2,L)*DXV(2)*PSAV*DTDIF
       QT (1,1,L)=QT (1,1,L)+FLUXQL
       TT (1,1,L)=TT (1,1,L)+FLUXTL
       DO 210 J=2,JMM1
         PSAV=0.5*(P(1,J)+P(1,J+1))
         !CAS FLUXQR=FDIFF*DQDY(J+1,L)*DXV(J+1)*PSAV*DTDIF
         FLUXQR=FDIFFQ*DQDY(J+1,L)*DXV(J+1)*PSAV*DTDIF
         FLUXQR=DMAX1( -0.5*QT(1,J+1,L), DMIN1(0.5*QT(1,J,L),FLUXQR))
         FLUXTR=FDIFF*DTDY(J+1,L)*DXV(J+1)*PSAV*DTDIF
         QT (1,J,L)=QT (1,J,L)+(FLUXQR-FLUXQL)
         TT (1,J,L)=TT (1,J,L)+(FLUXTR-FLUXTL)
	 VQHD(J,L)=-FLUXQL/(DXV(J)*0.5*(P(1,J)+P(1,J-1))*DTDIF) 
	 VTHD(J,L)=-FLUXTL/(DXV(J)*0.5*(P(1,J)+P(1,J-1))*DTDIF) 
         FLUXQL=FLUXQR
         FLUXTL=FLUXTR
  210  CONTINUE
         QT (1,JM,L)=QT (1,JM,L)-FLUXQR
         TT (1,JM,L)=TT (1,JM,L)-FLUXTR
	 J=JM
	 VQHD(J,L)=-FLUXQL/(DXV(J)*0.5*(P(1,J)+P(1,J-1))*DTDIF) 
	 VTHD(J,L)=-FLUXTL/(DXV(J)*0.5*(P(1,J)+P(1,J-1))*DTDIF) 
  200 CONTINUE
      DO 300 L=1,LM0
      DO 300 J=1,JM0
        Q (I,J,L)=QT (I,J,L)/FD(I,J)
        T (I,J,L)=TT (I,J,L)/FD(I,J)
       AJL(J,L,56)=AJL(J,L,56)+Q (I,J,L)*P(I,J)
  300 CONTINUE
      DOPK=1.  
      FD(I,1)=2.*FD(I,1)      
      FD(I,JM)=2.*FD(I,JM)     
      DO 65 J=2,JM           
      FDU=.5*(FD(I,J)+FD(I,J-1))   
      DO 65 L=1,LM
      VT(I,J,L)=V(I,J,L)*FDU  
 65   UT(I,J,L)=U(I,J,L)*FDU  
      DO 110 L=1,LM0
      DO 110 J=3,JM0
        DUDY(J,L)=(U (1,J,L)-U (1,J-1,L))/DYP(J)
        DVDY(J,L)=(V (1,J,L)-V (1,J-1,L))/DYP(J)
  110 CONTINUE
      DO 400 L=1,LM0
       PSAV=P(1,2)
       FLUXVL=FDIFU*DVDY(3,L)*DXP(2)*PSAV*DTDIF
       FLUXUL=FDIFU*DUDY(3,L)*DXP(2)*PSAV*DTDIF
       VT (1,2,L)=VT (1,2,L)+FLUXVL
       UT (1,2,L)=UT (1,2,L)+FLUXUL
       DO 410 J=3,JMM1
         PSAV=P(1,J)
         FLUXVR=FDIFU*DVDY(J+1,L)*DXP(J)*PSAV*DTDIF
         FLUXUR=FDIFU*DUDY(J+1,L)*DXP(J)*PSAV*DTDIF
         VT (1,J,L)=VT (1,J,L)+(FLUXVR-FLUXVL)
         UT (1,J,L)=UT (1,J,L)+(FLUXUR-FLUXUL)
	 VVHD(J-1,L)=-FLUXVL/(DXP(J-1)*P(1,J-1)*DTDIF) 
	 VUHD(J-1,L)=-FLUXUL/(DXP(J-1)*P(1,J-1)*DTDIF) 
         FLUXVL=FLUXVR
         FLUXUL=FLUXUR
  410  CONTINUE
         VT (1,JM,L)=VT (1,JM,L)-FLUXVR
         UT (1,JM,L)=UT (1,JM,L)-FLUXUR
	 J=JM
	 VVHD(J-1,L)=-FLUXVL/(DXP(J-1)*P(1,J-1)*DTDIF) 
	 VUHD(J-1,L)=-FLUXUL/(DXP(J-1)*P(1,J-1)*DTDIF) 
  400 CONTINUE
      DO 75 J=2,JM           
      RFDU=2./(FD(I,J)+FD(I,J-1))   
      DO 75 L=1,LM0
      V(I,J,L)=VT(I,J,L)*RFDU  
      U(I,J,L)=UT(I,J,L)*RFDU  
#if ( defined HR_DATA )
       if(L.le.4)then
        uyzhr(L,J)=U(I,J,L)
        vyzhr(L,J)=V(I,J,L)
       endif
#endif
 75   CONTINUE
      RETURN 
      END
