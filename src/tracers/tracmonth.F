
#include "ctrparam.h"

!	============================================================
!
!	CHEMMONTH1.F:   Subroutine for accumulating tracer mixing
!				ratios prepared for monthly averaging
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
!	052000	Chien Wang	rev.	
!	080200	Chien Wang	repack based on CliChem3 & add cpp
!	091901	Chien Wang	fix "write(180,*),"
!	092001	Chine Wang	add bc and oc
!	051804	Chien Wang	rev.
!
!	==========================================================

!	=====================
 	subroutine tracmonth1
!	=====================

#include "chem_para"
#include "chem_com"
#include "BD2G04.COM"
#include "tracer_com.h"


!	------------------------------------------------------


      do n=1,ntracers
       i = 1
       do k=1,nlev
       do j=1,nlat
        tracmonth(i,j,k,n)= tracmonth(i,j,k,n)+
     &           tracers(i,j,k,n)
       end do
       end do
      end do



      navtrac=navtrac+1


      return
       end

!	=====================
	subroutine tracmonth2
!	=====================

!	=============================================================
!
!	CHEMMONTH2.F:   Subroutine for calculating monthly averaged
!				mixing ratios of tracers
!	-------------------------------------------------------------
!	Author:		Chien Wang
!			MIT Joint Program on Science and Policy
!			    of Global Change
!	Last Revised:	June 29, 1999
!
!	=============================================================

#include "chem_para"
#include "chem_com"
#include "BD2G04.COM"
#include "tracer_com.h"


!	---------------------------------------------------------


	scale=1./float(navtrac)
      do n=1,ntracers
       i = 1
       do k=1,nlev
       do j=1,nlat
        tracmonth(i,j,k,n)= tracmonth(i,j,k,n)*scale
       end do
       end do
      end do



 	write(269)JYEAR,JDAY,tracmonth

	
      globalmass = 0.0
      do j=1,nlat
      do k=1,nlev
        globalmass = globalmass + airmass(1,j,k)
      enddo
      enddo

      do n=1,ntracers
       tracglob(n)= 0.0
       i = 1
       do k=1,nlev
       do j=1,nlat
        tracglob(n)= tracglob(n)+
     &          + airmass(i,j,k)
     &          * tracmonth(i,j,k,n)
       end do
       end do
       tracglob(n)= tracglob(n)/globalmass
      end do





      write(276,101)JMONTH,JDAY,tracglob
101   format(a6,i6,5E12.4)

        tracmonth = 0.0

        navtrac=0


      return
       end

     
