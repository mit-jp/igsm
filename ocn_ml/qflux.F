
#include "ctrparam.h"

C Routine used for mixed layer ocean model

      subroutine qflux

#include "BD2G04.COM"
#include "AGRID.h"
#include "ODIFF.COM"
#include "BDOCEAN.COM"

      LOGICAL PRNT
      dimension qfav(jm0),socean(jm0)
      DATA IFIRST/1/
      if(IFIRST.eq.1)then
        read(649),ARUNQ,QFLUXMN
        print *,'Q-flux from run ',arunq
        QFGAVI=0.0
        SWI=0.0
        do j=1,jm
c          print *,'DXYP(J)=',DXYP(J)
        do i=1,io
          POCEAN=1.-CLAND_OCEAN(i,j)
          if(POCEAN.gt.0.0)then
          SWI=SWI+POCEAN*DXYP(J)
          do nn=1,12
          QFGAVI=QFGAVI+POCEAN*DXYP(J)*QFLUXMN(I,J,NN)
          enddo
          endif
        enddo
        enddo
        print *,' QFGAVI=',QFGAVI/(12.*SWI)
        IFIRST=0
      endif

      if(JDATE.le.16)then
       do 723 j=1,JM
       do 723 i=1,IO
       if(CLAND_OCEAN(i,j).lt.1.)then
         QFLUXT(i,j)=((16-JDATE)*QFLUXMN(i,j,MONTH-1)+
     *            (JDATE+15)*QFLUXMN(i,j,MONTH))/31.
       else
          QFLUXT(I,J)=0.0
       endif
          POCEAN=1.-CLAND_OCEAN(i,j)
          QFGAV=QFGAV+POCEAN*DXYP(J)*QFLUXT(I,J)
          NQFGAV=NQFGAV+1
  723  continue
      else
       do 724 j=1,JM
       do 724 i=1,IO
       if(CLAND_OCEAN(i,j).lt.1.)then
        QFLUXT(i,j)=((JDATE-16)*QFLUXMN(i,j,MONTH+1)+
     *            (31-JDATE+16)*QFLUXMN(i,j,MONTH))/31.
       else
          QFLUXT(I,J)=0.0
       endif
          POCEAN=1.-CLAND_OCEAN(i,j)
          QFGAV=QFGAV+POCEAN*DXYP(J)*QFLUXT(I,J)
          NQFGAV=NQFGAV+1
  724  continue
      endif
      PRNT=TAU.gt.17520.0
      PRNT=.true.
      PRNT=.false.
      if(PRNT)then
       write(77, *),' QFLUX',' IDAY=',IDAY,' JYEAR=',JYEAR
       write(77, *),' JDATE=',JDATE,' JMONTH=',JMONTH
       write(77, *),' JDAY=',JDAY,' TOFDAY=',TOFDAY
       write(77, *),' MONTH=',MONTH,' TAU=',TAU
       write(77,*) ,'QFLUXT(5,37)=',QFLUXT(5,37)
      endif
c      write(78, *),' QFLUX  JDATE=',JDATE,' JMONTH=',JMONTH
c      write(78, *),' JDAY=',JDAY,' TOFDAY=',TOFDAY
c      write(78, *),' MONTH=',MONTH,' NQFGAV=',NQFGAV

      return
      end
