#include "ctrparam.h"
! 	This subroutine is called one a month by climate2tem
!	to adjust vegfraction given to TEM to conserve land area
	subroutine checkvegfraction
#include "BD2G04.COM"
#include "TEM.h"        
#include "CLM.h"        
!
#if ( defined CLM35 )
!
       do j=1,jm0
       do nc=1,ncoh
          if (Incohortarea(nc,j).ne.cohortarea(nc,j)) then
             print *,'From checkvegfraction'
             print '(2i3,2i20)',j,nc,Incohortarea(nc,j),cohortarea(nc,j)
             stop
          endif
          enddo
          enddo
        return
#endif
	end
