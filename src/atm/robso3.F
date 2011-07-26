
#include "ctrparam.h"

      subroutine  robso3(iyearr1)
#include "BD2G04.COM" 
#include "chem_para"
#include "TEM.h"
      character * 120 flo3,o3datadir
      character * 4 year
      character * 6 name
      data ifirst /1/
      iyearr=iyearr1
      if(ifirst.eq.1) then
       print *,'FROM ROBSO3'
       print *,'O3DATADIR=',o3datadir
       if(CLIMO3) print *,'With fixed O3 for 1860'
c      stop
       id2=index(o3datadir," ")
       print *,'O3DATADIR=',o3datadir(1:id2-1)
       ifirst=0
      endif
        if(CLIMO3) iyearr=1860
        write (year,'(i4)'),iyearr
        flo3=o3datadir(1:id2-1)//'o3bau'//year//'.dat46'
        print *,flo3
        open(561,file=flo3,status='old')
        do mn=1,12
        read (561,*),iiyear,nm1,name
        print *,'O3 for ',iiyear,nm1
        do j=1,jm0
         read (561,*),(obso3(n,j,mn),n=1,8)
        enddo
        enddo
        close (561)
      return
      end
