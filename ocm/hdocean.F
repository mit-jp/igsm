#include "ctrparam.h"

      SUBROUTINE HDOCEAN(R,AR,DX,DY,DXY,DML,FDIFF,DTDIF)
#include "ODIFF.COM"
      PARAMETER(JM0=N_LAT)
      DIMENSION R(JM0,LMO),DRDY(JM0+1,LMO),FDIFF(LMO)
     &   ,AR(JM0),RT(JM0,LMO),FD(JM0),DX(JM0),DY(JM0)
     &   ,DXY(JM0),DXO(JM0),DML(JM0)
      logical first
      data first /.true./
      JM=JM0
      JMM1=JM-1
      if(first)then
        print *,' HOR DIFF for R'
        print *,' JM0=',JM0,' LMO=',LMO
        print *,' FDIFF=',FDIFF
c       print *,'DX'
c       print *,DX
c       print *,'DY'
c       print *,DY
c       print *,'AR'
c       print *,AR
c       print *,'DXY'
c       print *,DXY
      DO 50 J=1,JM
 50   FD(J)=DXY(J)*AR(J)
      DO J=2,JM
      DXO(J)=DX(J)*0.5*(AR(J-1)+AR(J))
      END DO
c       print *,'FD'
c       print *,FD
        first=.false.
      endif	! first
      DO L=1,LMO
      DO 57 J=1,JM0
 57   RT (J,L)=R (J,L)*FD(J)
c     DRDY(3,L)=0.
      DO 100 J=1,JM
        if(AR(J-1).le.0.0.or.AR(J).le.0.0)then
          DRDY(j,L)=0.
        else
          DRDY(J,L)=(R (J,L)-R (J-1,L))/DY(J)
        endif
  100 CONTINUE
      DRDY(JM+1,L)=0.
      END DO
      DO 58 J=1,JM0
 58   RT (J,1)=RT(J,1)*DML(J)
c     print *,'DRDY'
c     print *,(DRDY(j,1),j=1,jm)
      RTT=0.0
      do j=1,jm
       if(AR(J).gt.0.0)then
       RTT=RTT+RT (J,1)
       endif
      enddo
      FLTOT=0.0
      DO 200 L=1,LMO
	cfr=1.
	cfm=1.
       FLUXL=0.0
       DO 210 J=1,JM
        if(AR(J).gt.0.0)then
c        print *,'210 j=',j
  	 if(L.eq.1)then
	  cfr=0.5*(DML(J+1)+DML(J))
c  cfm=DML(J)
	 end if
         if(J.lt.jm)then
          FLUXR=FDIFF(L)*DRDY(J+1,L)*DXO(J+1)*DTDIF
          FLUXR=FLUXR*cfr
         else
          FLUXR=0.0
         endif
c        FLUXR=DMAX1( -0.5*RT(1,J+1,L), DMIN1(0.5*RT(1,J,L),FLUXR))
c 	 if(L.eq.1)then
         FLTOT=FLTOT+(FLUXR-FLUXL)/cfm
c        endif
         RT (J,L)=RT (J,L)+(FLUXR-FLUXL)/cfm
         FLUXL=FLUXR
        endif
  210  CONTINUE
  200 CONTINUE
      RTT=0.0
      do j=1,jm
       if(AR(J).gt.0.0)then
       RTT=RTT+RT (J,1)
       endif
      enddo
      DO 300 J=1,JM
        if(AR(J).gt.0.0)then
         DO 301 L=1,LMO
           R (J,L)=RT (J,L)/FD(J)
  301    CONTINUE
           R (J,1)=R (J,1)/DML(J)
        endif
  300 CONTINUE
      RETURN
      END
