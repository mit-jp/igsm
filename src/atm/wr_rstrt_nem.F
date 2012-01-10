#include "ctrparam.h"
	subroutine wr_rstrt_nem
#include "BD2G04.COM"
#include "TEM.h"

C  unit 877 is last_clm file (really nem pickup)
! This file is not used anymore

C  unit 876 is last_nep file (older code, unit 367)
        write(876),adupt,temco2
     &          ,temch4,temn2o
        rewind 876

	return
	end
