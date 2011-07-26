
#include "ctrparam.h"

!	============================================================
!
!	CHEM_FILTER.F:   Subroutine for applying Shapiro filter to
!			   tracer gases in MIT Global Chemistry Model
!
!	------------------------------------------------------------
!
!	Author:		Chien Wang
!			MIT Joint Program on Science and Policy
!				of Global Change
!
!	----------------------------------------------------------
! 
!	Revision History:
!	
!	When	Who		What
!	----	----------	-------
!	050994	Chien Wang	rev.	
!	080100	Chien Wang	repack based on CliChem3 & add cpp
!	051804	Chien Wang	rev.
!
!	==========================================================

      Subroutine chemfilter(x,kstart,kend,norder)

#include "chem_para"
#include "chem_com"
#include "BD2G04.COM"

      common /WORK2/x1ji(72,46),x2ji(72,46),x3ji(72,46),
     &              x1  (72),   x2  (72),   x3  (72),   x4 (72)

      dimension x(nlon,nlat,nlev)

#if ( defined CPL_CHEM )

c      norder =16
c      norder =6
      xxx    =4.**norder
c-------------------------------------------------------
c  Applying the Shapiro filter to CFC11:
c

      do k=kstart,kend

      do j=1,n2dh
        x2ji(1,j)=1.0
        x1ji(1,j)=x(1,j,k)*x2ji(1,j)
        x3ji(1,j)=x1ji(1,j)
      end do

      call chemshap2d (mfilter, norder, xxx, nlon, nlat,1,1)

      do j=1,n2dh
        xtemp   = x1ji(1,j)/x2ji(1,j)
        xtemp   = max(0.0,xtemp)
        x(1,j,k)= xtemp
      end do

      end do

#endif

      return
       end
