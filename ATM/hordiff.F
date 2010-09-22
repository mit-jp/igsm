C $Header: /home/igsm-repository/ATM/hordiff.F,v 1.1 2010/08/30 15:36:35 jscott Exp $
C $Name:  $

#include "ctrparam.h"

!	==========================================================
!
!	HORDIFF.F:  subroutine for calculating horizontal 
!			diffusion of Q.
!
!	----------------------------------------------------------
!
!	Revision History:
!	
!	When	Who		What
!	----	----------	-------	
!	080100	Chien Wang	repack based on CliChem3 & M24x11,
!				  and add cpp.
!
!	==========================================================

      SUBROUTINE HORDIFF(DTDIF)

#include "BD2G04.COM"
 
      COMMON U,V,T,P,Q  
      DIMENSION VT(IM0,JM0,LM0),TT(IM0,JM0,LM0),PT(IM0,JM0), 
     *  QT(IM0,JM0,LM0),PU(IM0,JM0),FD(IM0,JM0),DQDY(JM0,LM0)
      COMMON/HDFLUX/VQHD(JM0,LM0)
      logical first
      data first /.true./
      I=1
      JMM1=JM0-1
      FDIFF=2.5E5
      FDIFF=5.0E5
      if(first)then
        print *,' HOR DIFF for Q'
        print *,' FDIFF=',FDIFF
        print *,'IM0=',IM0,' JM0=',JM0,' LM0=',LM0
        first=.false.
      endif
      DO 50 J=1,JM0                                            
 50   FD(I,J)=P(I,J)*DXYP(J)    
      DO 57 L=1,LM0         
      DO 57 J=1,JM0            
       AJL(J,L,56)=AJL(J,L,56)-Q (I,J,L)*P(I,J)
 57   QT (I,J,L)=Q (I,J,L)*FD(I,J) 
      DO 100 L=1,LM0
      DO 100 J=2,JM0
        DQDY(J,L)=(Q (1,J,L)-Q (1,J-1,L))/DYV(J)
  100 CONTINUE
      DO 200 L=1,LM0
       PSAV=0.5*(P(1,1)+P(1,2))
       FLUXL=FDIFF*DQDY(2,L)*DXV(2)*PSAV*DTDIF
       FLUXL=DMAX1( -0.5*QT(1,2,L), DMIN1(0.5*QT(1,1,L),FLUXL))
c      FLUXR=0.
       QT (1,1,L)=QT (1,1,L)+FLUXL
       DO 210 J=2,JMM1
         PSAV=0.5*(P(1,J)+P(1,J+1))
         FLUXR=FDIFF*DQDY(J+1,L)*DXV(J+1)*PSAV*DTDIF
         FLUXR=DMAX1( -0.5*QT(1,J+1,L), DMIN1(0.5*QT(1,J,L),FLUXR))
         QT (1,J,L)=QT (1,J,L)+(FLUXR-FLUXL)
	 VQHD(J,L)=-FLUXL/(DXV(J)*0.5*(P(1,J)+P(1,J-1))*DTDIF) 
         FLUXL=FLUXR
  210  CONTINUE
         QT (1,JM,L)=QT (1,JM,L)-FLUXR
	 J=JM
	 VQHD(J,L)=-FLUXL/(DXV(J)*0.5*(P(1,J)+P(1,J-1))*DTDIF) 
  200 CONTINUE
      DO 300 L=1,LM0
      DO 300 J=1,JM0
        Q (I,J,L)=QT (I,J,L)/FD(I,J)
       AJL(J,L,56)=AJL(J,L,56)+Q (I,J,L)*P(I,J)
  300 CONTINUE
      DOPK=1.  
      RETURN 
      END

