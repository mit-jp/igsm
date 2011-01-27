
#include "ctrparam.h"
      subroutine kvcarbon(iyr)
#include "BD2G04.COM"
#include "DRIVER.h"
      common/Garydiff/depthml(jm0),edzcar(jm0)
      common /Garykvct/cfkvct,edzcart(jm0)
      real depthml,edzcar,cfkvct,edzcart
      REAL pi, yy
      pi=4.*atan(1.)
      cfkvct=1.0
      if (iyr.ge.1991) then
         yy=(iyr-1990)/40.
         cfkvct=sqrt((0.5*pi-atan(yy))/(0.5*pi))
      endif
      print *,'iyr=',iyr,' cfkvct=',cfkvct
       do j=1,jm0
         edzcart(j)=cfkvct*edzcar(j)
      enddo
      return
      end
