
#include "ctrparam.h"

!	==========================================================
!
!	VFUNCT.F:  "unknown"
!
!	----------------------------------------------------------
!
!	Revision History:
!	
!	When	Who		What
!	-----	----------	-------	
!	080100	Chien Wang	repack based on M24x11, add cpp
!
!	==========================================================

	subroutine vfunct(v,z)

	zm=11000.
	d=2.71364
	a=1.27
	c=50900.
	zz=abs(z-zm)
	v=d*exp(-(zz**a/c))

	return 
	end
