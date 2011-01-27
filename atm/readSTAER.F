
#include "ctrparam.h"

	subroutine read_staer (nydata,amn)
#include "BD2G04.COM"
        parameter( naery=2000)
	character * 40 title
	character * 80 title1
	character * 120 fl_volaer
        character *5 xxx
        character *41 fordir
        dimension amn(jm0,12,naery),amnglav(12,naery)
        common/vaerosol/fl_volaer
        SS=0.0
        do j=1,jm
          SS=SS+DXYP(j)
        enddo
        do ii =1,25
        print *,fl_volaer
        enddo
        open (560,file=fl_volaer,form='formatted',status='old')
c       read(560,'(a40')title
c       print *,title
c       read(560,'(a80')title1
c       print *,title1
        do ny=1,naery
        do j=1,jm
         read (560,101,end=500),(amn(j,m,ny),m=1,12)
        enddo
        enddo
  500   continue
        nydata=ny-1
        print *,'Volcanic aerosol data for year =',1850+ny-2
        do n=1,nydata
        do nm=1,12
          amnglav(nm,n)=0.0
          do j=1,jm
          amnglav(nm,n)=amnglav(nm,n)+amn(j,nm,n)*DXYP(j)
          enddo
        enddo
        print 102,1849+n,(amnglav(nm,n)/SS,nm=1,12)
        enddo
        close (560)
  101	format(12f7.4)
  102	format(i5,12f7.4)
        return
        end

