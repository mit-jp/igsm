
#include "ctrparam.h"

C routine used in mixed layer ocean carbon model

        subroutine ocean_ml(dtocn,dtatm)

#      include "BD2G04.COM"
#      include "AGRID.h"

       dimension ndmonth(12),amonth(12)
       character *4 amonth
       character *4 anmonth
       logical dayli
       DATA ndmonth/31,28,31,30,31,30,31,31,30,31,30,31/
       DATA AMONTH/'JAN','FEB','MAR','APR','MAY','JUNE','JULY'
     &  ,'AUG','SEP','OCT','NOV','DEC'/


         KN=KOCEAN
!        print *,' dtocn from ocean_ml ',dtocn
!        print *,' dtatm from ocean_ml ',dtatm
         call ML_OCEAN(DTOCN,DTATM,KN)
         do j=1,jm0
         do i=1,IO0
           ODATA(i,j,1)=osst(j)
           ODATA(i,j,2)=foice(j)
           ODATA(i,j,3)=aoice(j)
           ODATA(i,j,4)=osst2(j)
           ODATA(i,j,5)=osst3(j)
           GDATA(i,j,1)=snowice(j)
           GDATA(i,j,3)=tice1(j)
           GDATA(i,j,7)=tice2(j)
         enddo
         enddo
c        JPR=36
c        write(78, *)'OCEAN_ML JPR=',JPR
c        write(78, *),osst(JPR),aoice(JPR),foice(JPR)
c        write(78, *),snowice(JPR),tice1(JPR),tice2(JPR)
c        JPR=37
c        write(78, *)'OCEAN_ML JPR=',JPR
c        write(78, *),osst(JPR),aoice(JPR),foice(JPR)
c        write(78, *),snowice(JPR),tice1(JPR),tice2(JPR)
        go to 555
       print *,' Everything from N.P. to S.P.'
        print *,' '
c       print *,'Month=',anmonth,'DAY=',jdatoc
        print *,'DAY=',jdayoc
        print *,' SST'
        print '(12f7.1)',(osst(j),j=jm0,1,-1)
        print *,' ICE'
        print '(12f7.1)',(aoice(j)/916.,j=jm0,1,-1)
        print *,' SNOW_ICE'
        print '(12f7.1)',(snowice(j),j=jm0,1,-1)
        print *,' FICE'
        print '(12f7.1)',(foice(j),j=jm0,1,-1)
        print *,' TEMICE_1'
        print '(12f7.1)',(tice1(j),j=jm0,1,-1)
c       print *,' TEMICE_2'
c       print '(12f7.1)',(tice2(j),j=jm0,1,-1)
 555    continue
       return
       end
