#include "ctrparam.h"

	subroutine read_orbital 
#include "BD2G04.COM"
#include "ORB.h"
	character * 40 title
	character * 80 title1
	character * 256 fl_orbital
        character *5 xxx
        common/orbitalfl/fl_orbital
        SS=0.0
        do j=1,jm
          SS=SS+DXYP(j)
        enddo
!       do ii =1,25
        print *,fl_orbital
!       enddo
        open (560,file=fl_orbital,form='formatted',status='old')
        do ii =1,ncomments
        read(560,'(a80)'),title1
        enddo
        print *,title1
        do ny=1,nyorbit
         read (560,*,end=500), yorbit(ny),eccenh(ny),obliqh(ny)
     &     ,perih(ny)
!       print *,yorbit(ny),eccenh(ny),obliqh(ny)
        enddo
  500   continue
        close (560)
  101	format(12f7.4)
  102	format(i5,12f7.4)
        return
        end

