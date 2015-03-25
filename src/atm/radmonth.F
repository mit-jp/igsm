#include "ctrparam.h"
      subroutine radmonth 
#include "BD2G04.COM"
      common /RADMN/ nradmn,radmonthly(JM0,LM0+3,9)
      do k=1,9
      do l=1,lm+3
      do j=1,jm
       radmonthly(j,l,k)=radmonthly(j,l,k)/nradmn
      enddo
      enddo
      enddo
      write (469)radmonthly
       radmonthly = 0.0
       nradmn=0
      return
      end
