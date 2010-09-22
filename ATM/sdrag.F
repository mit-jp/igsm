C $Header: /home/igsm-repository/ATM/sdrag.F,v 1.1 2010/08/30 15:36:35 jscott Exp $
C $Name:  $

#include "ctrparam.h"

      SUBROUTINE SDRAG(WLMMAX,JWLMMAX)                                  7801.
C****                                                                   7802.
C**** THIS SUBROUTINE PUTS A DRAG ON THE WINDS ON THE TOP LAYER OF      7803.
C**** THE ATMOSPHERE                                                    7804.
C****                                                                   7805.
#include "BD2G04.COM"                                                   7806.

      COMMON U,V,T,P,Q                                                  7807.
      COMMON/WORK1/CONV(IM0,JM0,LM0),PK(IM0,JM0,LM0)                    7808.
      real TAUJ(JM0+1)
      DO 100 J=2,JM                                                     7809.
      I=IM                                                              7810.
      DO 100 IP1=1,IM                                                   7811.
      SP=.25*(P(I,J-1)+P(IP1,J-1)+P(I,J)+P(IP1,J))                      7812.
      WLM=SQRT(U(I,J,LM)*U(I,J,LM)+V(I,J,LM)*V(I,J,LM))                 7813.
      if(WLM.gt.WLMMAX)then
         WLMMAX=WLM
         JWLMMAX=J
      endif
      RHO=PTOP/(RGAS*T(I,J,LM)*PK(I,J,LM))                              7814.
      CDN=DUMMY1(1)+DUMMY1(2)*WLM                                       7815.
      TAUJ(J)=CDN*100.*RHO*WLM*U(I,J,LM)
      X=NDYN*DT*RHO*CDN*WLM*GRAV/(SP*DSIG(LM))                          7816.
      U(I,J,LM)=U(I,J,LM)*(1.-X)                                        7817.
      V(I,J,LM)=V(I,J,LM)*(1.-X)                                        7818.
  100 I=IP1                                                             7819.
C****                                                                   5973.
      TAUJ(1)=TAUJ(2)
      TAUJ(JM+1)=TAUJ(JM)
      do J=1,JM
      PLAND=FDATA(I,J,2)                                                5974.
      PWATER=1.-PLAND                                                   5975.
      PLICE=FDATA(I,J,3)*PLAND                                          5976.
      PEARTH=PLAND-PLICE                                                5977.
      POICE=ODATA(I,J,2)*PWATER                                         5978.
      POCEAN=PWATER-POICE                                               5979.
      if(POCEAN.LE.1.E-5)then
         POCEAN=0.
         POICE=PWATER
      endif
      TTOFR=PEARTH+PLICE+POICE+POCEAN
      if(abs(TTOFR-1.).gt.1.e-3)then
        print *,'From sdrag TTOFR=',TTOFR
        print *,' J=',J,' PLAND=',PLAND,' POCEAN=',POCEAN
        print *,' I=',I,' PWATER=',PWATER,' POICE=',POICE
        print *, 'ODATA(I,J,2)=',ODATA(I,J,2)
        stop
      end if
C  Fric
         AJ(J,38)=AJ(J,38)+(TAUJ(j)+TAUJ(j+1))*POCEAN
         BJ(J,38)=BJ(J,38)+(TAUJ(j)+TAUJ(j+1))*PLAND
         CJ(J,38)=CJ(J,38)+(TAUJ(j)+TAUJ(j+1))*POICE
       enddo
C  Fric
      RETURN                                                            7820.
      END                                                               7821.
