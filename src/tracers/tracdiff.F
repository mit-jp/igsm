#include "ctrparam.h"

!	============================================================
!
!      CHEMDIFF.F:    Subroutine for calculating horizontal
!       		diffusion of MIT Global Chemistry Model
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
!	013096	Chien Wang	rev.	
!	080100	Chien Wang	repack based on CliChem3 & add cpp
!	051804	Chien Wang	rev.
!
!	==========================================================

      Subroutine tracdiff(ifdiff,x00,x11,dta)

#include "chem_para"
#include "chem_com"
#include "BD2G04.COM"

      dimension x00  (nlon,nlat,nlev)
      dimension x11  (nlon,nlat,nlev)

      dimension dcdy(nlat,nlev)


c-------------------------------------------------------
c   Definitions of parameters:
c
      istart=1
      iend  =nlon
	
c
c 013096
c   fktdif span from 2.e6 in the first three years to
c	5.e5 or 1.e6 in twenty years and maintain this value
c	thereafter:
c
c	xxx    = float(myyear - 3)/20.0
c	xxx    = amin1(1.0,amax1(0.0,xxx))
c	fktdif = (20.0 - xxx * 10.0)*1.e5	! m2/s

c      fktdif = 4.e6 	 !m2/s

c	if(ifdiff.eq.1)then
c	  fktdif = 2.e6
c	else if(ifdiff.eq.2)then
c	  fktdif = 3.e6
c	endif

c 111596:
c	fktdif = float(ifdiff)*1.e6
	fktdif = float(ifdiff)*1.e5

c=====
c   Calculate dcdy:
c
      do i=istart,iend
      do j=2,nlat
      do k=1,nlev
        dcdy(j,k)=(x11(i,j,k)-x11(i,j-1,k))
     &           /dyv(j)
      end do
      end do
      end do

c=====
c   Calculate meridional eddy diffusion:
c
      do k=1,nlev
        paver  = 0.5*(p00(1,1)+p00(1,2))
        fluxl  =-fktdif
     &         /dyv(2)*dcdy(2,k)*dta 
     &         * paver
        fluxl=max(-0.5*x00(1,2,k), min(0.5*x00(1,1,k),fluxl))
      do j=2,nlat1
        paver  = 0.5*(p00(1,j)+p00(1,j+1))
        fluxr     =-fktdif
     &            /dyv(j+1)*dcdy(j+1,k)*dta
     &         * paver
        fluxr=max(-0.5*x00(1,j+1,k),min(0.5*x00(1,j,k),fluxr))
        x00(1,j,k)=x00(1,j,k)-(fluxr-fluxl)
        fluxl=fluxr
       end do
       end do

c	call traccheck(x00)


	return
	 end

