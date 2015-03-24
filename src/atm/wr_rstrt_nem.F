#include "ctrparam.h"
	subroutine wr_rstrt_nem
#include "BD2G04.COM"
#include "TEM.h"

C  unit 877 is last_clm file (really nem pickup)
! This file is not used anymore

       if (1.eq.0) then
	write(877),co24tem,
     &            o34tem,
     &            temp4tem,
     &            dtem4tem,
     &            sws4tem,
     &            pre4tem,
     &            strmdur,
     &            qstrm,
     &            aet,
     &            sh2o1m,
     &            sh2o2m,
     &            swe,sfr,drn,
     &            daytsoil,
     &            daysho,
     &            rsh2o
c    &            istorms,
c    &            istrmdry
cprint *,'AET for J=23'
c       print *,(aet(nc,23),nc=1,ncoh)
        rewind 877
        endif

C  unit 876 is last_nep file (older code, unit 367)
        write(876),adupt,temco2
     &          ,temch4,temn2o
        rewind 876

	return
	end
