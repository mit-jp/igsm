C $Header: /home/igsm-repository/OCN_ML/wrrstrt_ocean.F,v 1.1 2010/09/03 19:37:33 jscott Exp $
C $Name:  $

#include "ctrparam.h"

C write mixed layer ocean model restart file

      subroutine wrrstrt_ocean

#include "BD2G04.COM"   
#include "AGRID.h"
#include "ODIFF.COM"
#include "BDOCEAN.COM"   

      rewind 658
      print *,'Write data for ocean restart'
      print *,'AEXP=',AEXP
      print *,tofday,jday,jdate,jmonth,jyear
      write (658),AEXP
      write (658),tofday,jday,jdate,jmonth,jyear
      write (658),ODATA_OCEAN,GDATA_ICE,Z1OOLD,OLDICEM
      print *,' ODATA_OCEAN(i,7,1)'
      print *,(ODATA_OCEAN(i,7,1),i=1,io0,6)
       if(TRRUN)then
      write (658),RTGO2D
      endif
      rewind 658
      return
      end
