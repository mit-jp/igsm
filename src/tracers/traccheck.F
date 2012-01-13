
#include "ctrparam.h"

!	============================================================
!
!	CHEMCHECK.F:	Subroutine for checking negative values and
!			    forcing them to be passitive ones
!				   of MIT Global Chemistry Model
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
!	081795	Chien Wang	rev.	
!	080100	Chien Wang	repack based on CliChem3 & add cpp
!	051804	Chien Wang	rev.
!
!	==========================================================

      Subroutine chemcheck(x00,x11,dta)

#include "chem_para"
#include "BD2G04.COM"

      dimension x00  (nlon,nlat,nlev)

!	----------------------------------------------------------

c ===
c === 081295: check for negative values:
c ===

	oneover3 = 1./3.

	do i=1,nlon
	do j=1,nlat

	  if(x00(i,j,1).lt.0.0) then
	    xmean =  max(0.0,
     &                  (x00(i,j,1) + x00(i,j,2))*0.5 )
            x00(i,j,1) = xmean
            x00(i,j,2) = xmean
	  endif

	do k =2,nlev1
	  if(x00(i,j,k).lt.0.0)then
	     xmeanv = max(0.0,
     &                (x00(i,j,k) 
     &               + x00(i,j,k-1) + x00(i,j,k+1))*oneover3 )
	     x00(i,j,k-1) = xmeanv
	     x00(i,j,k)   = xmeanv
	     x00(i,j,k+1) = xmeanv
	  endif
	end do

	  if(x00(i,j,nlev).lt.0.0) then
	    xmean = max(0.0,
     &                       (x00(i,j,nlev)
     &                      + x00(i,j,nlev1))*0.5 )
            x00(i,j,nlev)  = xmean
            x00(i,j,nlev1) = xmean
	  endif

	end do
	end do


      return
       end

