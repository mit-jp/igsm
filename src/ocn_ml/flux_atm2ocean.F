
#include "ctrparam.h"

C ML Ocean Model

      subroutine flux_atm2ocean

#include "BD2G04.COM"
#include "AGRID.h"
#include "ODIFF.COM"
#include "BDOCEAN.COM"

      logical prnt
      STBO=.5672573E-7
      SPERDAY=24.*3600.
       do j=1,JM
         do i=1,IO
            PLAND=CLAND_OCEAN(I,J)
            POCEAN= (1.-PLAND)
            POICE= ODATA_OCEAN(I,J,2)*(1.-PLAND)
            POCEAN=1.-PLAND
            PWATER=1.-PLAND-POICE
          prnt=i.eq.21.and.j.eq.37
          prnt=.false.
          if(POCEAN.gt.0.0)then
          if(PWATER.gt.0.0)then
           E0(I,J,1)=hfluxo(j)-dhfodtg(j)*(osst(j)-ODATA_OCEAN(i,j,1))
           E0(I,J,1)=DTOCEAN*E0(I,J,1)
           EVAPOR(I,J,1)=evao(j)+
     &     devodtg(j)*(osst(j)-ODATA_OCEAN(i,j,1))
           EVAPOR(I,J,1)=DTOCEAN*EVAPOR(I,J,1)/SPERDAY
          endif
          if(POICE.gt.0.0)then
           E0(I,J,2)=hfluxi(j)-dhfidtg(j)*(tice1(j)-GDATA_ICE(i,j,2))
           E0(I,J,2)=DTOCEAN*E0(I,J,2)
           EVAPOR(I,J,2)=evai(j)+
     &     devidtg(j)*(tice1(j)-GDATA_ICE(i,j,2))
           EVAPOR(I,J,2)=DTOCEAN*EVAPOR(I,J,2)/SPERDAY
          endif
          PREC(I,J)=DTOCEAN*precip(j)/SPERDAY
          TPREC(I,J)=tempr(j)
            HFAV=HFAV+(E0(I,J,1)*PWATER+E0(I,J,2)*POICE)*DXYP(J)
            HFOAV=HFOAV+E0(I,J,1)*PWATER*DXYP(J)
            HFIAV=HFIAV+E0(I,J,2)*POICE*DXYP(J)
            EVIAV=EVIAV+EVAPOR(I,J,2)*POICE*DXYP(J)
            PRIAV=PRIAV+PREC(I,J)*POICE*DXYP(J)
            TPRIAV=TPRIAV+TPREC(I,J)*POICE*DXYP(J)
            PWTRAV1=PWTRAV1+PWATER*DXYP(J)
          if(prnt)then
           write(77,*) ,' '
           write(77,*) ,'FROM FLUX  i=',i,' j=',j
           write(77,*),'DAY=',JDAY,' TAUO=',TAUO
            write(77,*),'PLAND=',PLAND
           write(77,*) ,' PWATER=',PWATER
           write(77,*) ,' POICE=',POICE
           if(PWATER.gt.0.0)then
           write(77,*) ,'hfluxo(j)=',hfluxo(j)
     &     ,' dhfodtg(j)=',dhfodtg(j)
           write(77,*) ,'E0(I,J,1)=',E0(I,J,1)/DTOCEAN
           write(77,*) ,'osst(j)=',osst(j)
           write(77,*) ,'sst(i,j)=',ODATA_OCEAN(i,j,1)
           write(77,*) ,'evao(j) mm/day =',evao(j)
     &     ,' devodtg(j)=',devodtg(j)
           write(77,*) ,'evao(j)=',DTOCEAN*evao(j)/SPERDAY
           write(77,*) ,'EVAPOR(I,J,1)=',EVAPOR(I,J,1)
           write(77,*) ,'PREC(I,J)=',PREC(I,J)
           write(77,*) ,'TPREC(I,J)=',TPREC(I,J)
           endif
           if(POICE.gt.0.0)then
           write(77,*) ,'OCEAN ICE',' POICE=',POICE
           write(77,*) ,'hfluxi(j)=',hfluxi(j)
     &     ,' dhfidtg(j)=',dhfidtg(j)
           write(77,*) ,'E0(I,J,2)=',E0(I,J,2)/DTOCEAN
           write(77,*) ,'tice1(j)=',tice1(j)
           write(77,*) ,'tice1(i,j)=',GDATA_ICE(i,j,2)
           write(77,*) ,'evai(j) mm/day =',evai(j)
     &     ,' devidtg(j)=',devidtg(j)
           write(77,*) ,'evai(j)=',DTOCEAN*evai(j)/SPERDAY
           write(77,*) ,'EVAPOR(I,J,2)=',EVAPOR(I,J,2)
           endif
          endif
          endif 	! PWATER > 0
         enddo
         enddo
         return
         end
