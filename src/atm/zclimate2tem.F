
#include "ctrparam.h"
! 	This subroutine is called one a month by driver
!	to zero out  climate data for TEM
	subroutine zclimate2tem
#include "BD2G04.COM"
#include "TEM.h"        
	do j=1,jm0
	  co24tem(j)=0.0
	  temp4tem(j)=0.0
	  sws4tem(j)=0.0
	  pre4tem(j)=0.0
	  do n=1,no3
	   o34tem(n,j)=0.0
	  enddo
	enddo
	return
	end
