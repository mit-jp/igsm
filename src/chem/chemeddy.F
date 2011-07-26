
#include "ctrparam.h"

!	============================================================
!
!	CHEMEDDY.F:    Subroutine for calculating zonal-average eddy
!			 diffusion of MIT Global Chemistry Model
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
!	051804	Chien Wang	rev. for 46x11
!
!	==========================================================

      Subroutine chemeddy(ifdiff,x00,x11,dta)

#include "chem_para"
#include "chem_com"
#include "BD2G04.COM"

      dimension x00  (nlon,nlat,nlev)
      dimension x11  (nlon,nlat,nlev)

      dimension vc   (nlat,nlev)
      dimension beta5(nlat)

      dimension dcdy(nlat,nlev)
      dimension dcdz(nlat,nlev)
      dimension dcdc(nlat,nlev)

!	----------------------------------------------------------

#if ( defined CPL_CHEM )

c-------------------------------------------------------
c   Definitions of parameters:
c

      istart=1
      iend  =nlon

      beta5(1)=0.0
      do j=2,nlat1
        beta5(j)=0.573*sqrt(beta2(j))/(1.-0.427*beta2(j)**0.302)
      end do
      beta5(nlat)=0.0

c=====
c   Calculate dcdy and dcdz:
c
      do 5  i=istart,iend
      do 5  j=2,nlat
      do 5  k=1,nlev
        dcdy(j,k)=(x11(i,j,k)-x11(i,j-1,k))/dyv(j)
5     continue

      do 6 i=istart,iend
      do 6 j=1,nlat
      do 6 k=1,nlev1
        dcdz(j,k)=-(x11(i,j,k+1)-x11(i,j,k))*deltap(j,k)
6     continue

      do 61 i=istart,iend
      do 61 j=1,nlat
      do 62 k=2,nlev
        dcdz(j,k)=dcdz(j,k-1)*dp2dz(j,k-1)
62    continue
      dcdz(j,1)=dcdz(j,2)
61    continue

      do 7 j=2,nlat1
      do 7 k=2,nlev1
        dcdc(j,k)=dcdz(j,k)*4.0
     &           /(dcdy(j,k)  +dcdy(j+1,k)
     &            +dcdy(j,k+1)+dcdy(j+1,k+1)+1.e-20)
7     continue

      do 8 j=2,nlat1
        alamor   =beta5(j)/beta2(j)
        alamor2  =alamor  /beta2(j)
        oneoalam1=1./(1.+beta5(j))
      do 8 k=2,nlev1
        dcdc(j,k)=oneoalam1*beta1(j)*beta3(j,k)
     &           *(1.0+alamor
     &            +beta3(j,k)*0.25*beta1(j)*dcdc(j,k)
     &            *(1.+alamor2))
8     continue      

c=====
c   Calculate meridional eddy diffusion:
c
      do 10 k=1,nlev
        paver  = 0.5*(p00(1,1)+p00(1,2))
        fluxl  =-fkt(2,k)
     &         /dyv(2)*dcdy(2,k)*dta 
     &         * paver
        fluxl=max(-0.5*x00(1,2,k),min(0.5*x00(1,1,  k),fluxl))
        vc(2,k)=fluxl/(paver+1.e-20)
      do 11 j=2,nlat1
        paver  = 0.5*(p00(1,j)+p00(1,j+1))
        fluxr     =-fkt(j+1,k)
     &            /dyv(j+1)*dcdy(j+1,k)*dta
     &         * paver
        fluxr=max(-0.5*x00(1,j+1,k),min(0.5*x00(1,j,k),fluxr))
        vc (j+1,k)=fluxr/(paver+1.e-20)
        x00(1,j,k)=x00(1,j,k)-(fluxr-fluxl)
        fluxl=fluxr
11    continue
10    continue

c=====
c   Calculate vertical eddy diffusion:
c
c 112696 changed also in eddypa.f for beta4
c
      do 12 j=2,nlat1
        fluxb=0.0
      do 14 k=1,n_tropopause	! ktrop = 7 for both 9 and 11 layer model
        fluxt=0.25*(vc(j,k)+vc(j,k+1)+vc(j+1,k)+vc(j+1,k+1))
     &       *dcdc(j,k+1)
     &       *beta4(j,k+1)
     &       *p00(1,j)

c        fluxt=max(-0.5*x00(1,j,k)  *dsig(k),
c     &        min( 0.5*x00(1,j,k+1)*dsig(k+1),fluxt))
c        if(fluxt*dcdz(j,k+1).lt.0.0) fluxt=0.0
c        x00(1,j,k)=x00(1,j,k)+(fluxt-fluxb)/dsig(k)

        fluxt=max(-0.5*x00(1,j,k+1)*dsig(k+1),
     &        min( 0.5*x00(1,j,k)  *dsig(k),fluxt))
        if(fluxt*dcdz(j,k+1).gt.0.0) fluxt=0.0
        x00(1,j,k)=x00(1,j,k)-(fluxt-fluxb)/dsig(k)
        fluxb=fluxt
14    continue
12    continue

c	write(6,*)"FKT = "
c     	write(6,*)fkt
c	write(6,*)"VC = "
c	write(6,*)vc
c	write(6,*)"DCDY = "
c     	write(6,*)dcdy
c
c  040895 test:
c

1996	continue

c ======
c 013096 
c Apply horizontal diffussion to some tracers
c	to reduce initialization errors in the
c	global distribution:
c
      if(ifdiff.ne.0) call chemdiff(ifdiff,x00,x11,dta)

	call chemcheck(x00)

#endif

      return
       end

