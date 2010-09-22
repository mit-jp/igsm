C $Header: /home/igsm-repository/OCN_ML/qflux_rel.F,v 1.1 2010/09/03 19:37:33 jscott Exp $
C $Name:  $

#include "ctrparam.h"

C Routine used for mixed layer ocean model

      subroutine qflux_rel

#include "BD2G04.COM"
#include "AGRID.h"
#include "ODIFF.COM"
#include "BDOCEAN.COM"

      COMMON /HF2ML/ HFML(IO0,JM0)
      dimension qfav(jm0),socean(jm0)
      DATA IFIRST/1/
      logical PRNT
      if(IFIRST.eq.1)then
      XLHM=334000.
      CW=1000.*4185.
      print *,'DFLDSST=',DFLDSST
c     TAUREL=45*24*3600.
c     print *,'TAUREL=',TAUREL
c     print *,'DFLDSST=',CW*100./TAUREL
      endif
c     QFGAV=0.0
c     SO=0.0
      DO 650 J=1,JM
       qfav(j)=0.0
       socean(j)=0.0
       DO 650 I=1,IO
      PRNT=I.eq.5.and.J.eq.37
      PRNT=.false.
      PRNT=J.eq.36.and.((I.ge.21.and.I.le.27).or.(I.ge.64.and.I.le.-66))
     &   .and.JDATE.eq.16.and.JYEAR.eq.1
      PRNT=.false.
       if(CLAND_OCEAN(i,j).lt.1.)then
C      Used in the run 1043.03
c       if(J.eq.37.or.J.eq.36)then
c           ODATAT(i,j,2)=0.5*ODATAT(i,j,2)
c       endif       
c       if(J.eq.37.and.I.ge.64.and.I.le.66)then
c          if(ODATAT(I,J,1).lt.-1.0)ODATAT(I,J,1)=-1.0
c       endif       
C      Used in the run 1040.03
        OICEM(i,j)=ODATA_OCEAN(i,j,3)*ODATA_OCEAN(i,j,2)
        OICEMT(i,j)=ODATAT(i,j,3)*ODATAT(i,j,2)
        if(OICEMT(i,j).gt.0.0.and.ODATAT(I,J,1).gt.0.0)then
         OICEMT(i,j)=0.0
        endif
        QICE=XLHM*(OICEM(i,j)-OICEMT(i,j))/(CW*Z1O(I,J))
        QTEM=ODATAT(I,J,1)-ODATA_OCEAN(I,J,1)
        QFLUXT(I,J)=DFLDSST*(QTEM+QICE)
c for taurel
c       QICE=XLHM*(OICEM(i,j)-OICEMT(i,j))/TAUREL
c       QTEM=CW*Z1O(I,J)*(ODATAT(I,J,1)-ODATA_OCEAN(I,J,1))/TAUREL
c       QFLUXT(I,J)=(QTEM+QICE)
c for taurel
       else
          QICE=0.0
          QTEM=0.0
          QFLUXT(I,J)=0.0
       endif
          if(prnt)then
           write(77,*),'JDATE=',JDATE,' MONTH=',MONTH
           write(77,*),'QFLUX_REL ',j,i,' Focean=',1.-CLAND_OCEAN(i,j)
           write(77,*),'ICEOBS=',OICEMT(i,j),' ICE=',OICEM(i,j)
           write(77,*),'SSTOBS=',ODATAT(I,J,1),
     &       ' SST=',ODATA_OCEAN(I,J,1)
c          write(77,*),'DICE=',QICE,' DSST=',QTEM
c          write(77,*),'Z1O(I,J)=',Z1O(I,J)
           write(77,*),'QFICE=',DFLDSST*QICE,' QFSST=',DFLDSST*QTEM
           write(77,*),'QFLUXT(I,J)=',QFLUXT(I,J)
          endif
          POCEAN=1.-CLAND_OCEAN(i,j)
          HFMLAV=HFMLAV+POCEAN*DXYP(J)*HFML(I,J)
          QFGAV=QFGAV+POCEAN*DXYP(J)*QFLUXT(I,J)
          QFTAV=QFTAV+DFLDSST*QTEM*POCEAN*DXYP(J)
          QFIAV=QFIAV+DFLDSST*QICE*POCEAN*DXYP(J)
  650   CONTINUE
        IFIRST=0
        return
        end
