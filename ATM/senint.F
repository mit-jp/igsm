
#include "ctrparam.h"

       subroutine senint(s)
       parameter(nsen=200)
       common/cldfdb/coefcl(3),CLDFEED,SEN_dat
       logical CLDFEED
       character * 120 SEN_dat
       dimension sens(nsen),cfcl(nsen)
       data ifisrt /0/
       if(ifirst.eq.0) then
         open (543,file=SEN_dat,form='formatted',status='old')
         do i=1,nsen
           read (543,*,end=500),cfcl(i),sens(i)
c          read (543,100),cfcl(i),sens(i)
         enddo
 500     continue
         nsenr=i-1
         ifirst=1
         print *,' Number of points in SEN_dat file ',nsenr
       endif
       if(s.lt.sens(1))then
        ip=2
        go to 200
       elseif(s.ge.sens(nsenr))then
         ip=nsenr
         go to 200
       endif
       do i=1,nsenr
         if(s.lt.sens(i))then
            ip=i
            go to 200
         endif
       enddo
  200  continue
       im=ip-1
       cf=(cfcl(im)*(sens(ip)-s)+cfcl(ip)*(s-sens(im)))
     &      /(sens(ip)-sens(im))
       do i=1,3
         coefcl(i)=cf
       enddo
  100  format(f9.6,f5.2)
       return
       end
