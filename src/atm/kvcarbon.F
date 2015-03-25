#include "ctrparam.h"
      subroutine kvcarbon(iyr)
#include "ATMSIZE.h"
#include "DRIVER.h"
      common/Garydiff/depthml(jm0),edzcar(jm0)
      common /Garykvct/cfkvct,edzcart(jm0)
      common/cfocarb/KVCT
      logical KVCT
      real depthml,edzcar,cfkvct,edzcart
      REAL pi, yy1
      pi=4.*atan(1.)
      cfkvct=1.0
      if(KVCT) then
      if (iyr.ge.1991) then
         yy1=(iyr-1990)/40.
         cfkvct=sqrt((0.5*pi-atan(yy1))/(0.5*pi))
      endif
      endif
      print *,'iyr=',iyr,' cfkvct=',cfkvct
       do j=1,jm0
         edzcart(j)=cfkvct*edzcar(j)
      enddo
      return
      end
