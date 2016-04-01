
#include "ctrparam.h"

!	============================================================
!
!	CHEMADV.F:   Calculating advection and eddy
!			transport of chemical species
!
!			chemadv0
!			chemadv
!			chemadv2
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
!	062298	Chien Wang	rev.	
!	080200	Chien Wang	repack based on CliChem3 & add cpp
!	101800	Chien Wang	replaced if_3gases with cpp
!	092001	Chien Wang	add bc and oc
!	051804	Chien Wang	rev. for 46x11
!
!	==========================================================

!	=======================
	subroutine chemadv0 (dt)
!	=======================


#include "chem_para"
#include "chem_com"

	dimension tracert(nlon,nlat,nlev)


#if ( defined CPL_CHEM )
c
c-----------------------------------
c Calculating advection and eddy diffusion:
c

c
c CFC11:
c

      do i=1,n3d
        tracert(i,1,1)=cfc11(i,1,1)
      enddo

      call chemadv(8,
     &             tracert,cfc11,cfc110,dt,ddepref)
c 
c CFC12:
c
      do i=1,n3d
        tracert(i,1,1) = cfc12(i,1,1)
      enddo

      call chemadv(8,
     &             tracert,cfc12,cfc12,dt,ddepref)
c 
c N2O:
c
      do i=1,n3d
        tracert(i,1,1) = xn2o(i,1,1)
      enddo

      call chemadv(8,
     &             tracert,xn2o,xn2o,dt,ddepref)

! === if hfc, pfc, and sf6 are included:
#if ( defined INC_3GASES )

! === 032698
! === HFC134a:
! ===

      do i=1,n3d
        tracert(i,1,1) = hfc134a(i,1,1)
      enddo

      call chemadv(8,
     &             tracert,hfc134a,hfc134a,dt,ddepref)

! ===
! === PFC:
! ===

      do i=1,n3d
        tracert(i,1,1) = pfc(i,1,1)
      enddo

      call chemadv(8,
     &             tracert,pfc,pfc,dt,ddepref)

! ===
! === SF6:
! ===

      do i=1,n3d
        tracert(i,1,1) = sf6(i,1,1)
      enddo

      call chemadv(8,
     &             tracert,sf6,sf6,dt,ddepref)

#endif

!
! === Black Carbon:
!

      do i=1,n3d
        tracert(i,1,1) = bcarbon(i,1,1)
      enddo

      call chemadv(8,
     &             tracert,bcarbon,bcarbon,dt,ddepbc)

!
! === Organic Carbon:
!

      do i=1,n3d
        tracert(i,1,1) = ocarbon(i,1,1)
      enddo

      call chemadv(8,
     &             tracert,ocarbon,ocarbon,dt,ddepoc)


! ===

c 
c O3:
c
      do i=1,n3d
        tracert(i,1,1) = o3(i,1,1)
      enddo

c 051295 use chemadv2 with different top vbc:

      call chemadv2(8,
     &              tracert,o3,o3,dt,ddepo3)
c      call chemadv(tracert,o3,o3,dt)

c 051698 use prescribed top o3:

	do k=n_tropopause+1,nlev
	do j=1,nlat
	  o3(1,j,k) = o3top(j,k-n_tropopause,mymonth)
	end do
	end do
c 
c CO:
c
      do i=1,n3d
        tracert(i,1,1) = co(i,1,1)
      enddo

      call chemadv(8,
     &             tracert,co,co,dt,ddepref)

c 
c CO2:
c
      do i=1,n3d
        tracert(i,1,1) = zco2(i,1,1)
      enddo

      call chemadv(8,
     &             tracert,zco2,zco2,dt,ddepref)

c 
c NO:
c
      do i=1,n3d
        tracert(i,1,1) = xno(i,1,1)
      enddo

      call chemadv(8,
     &             tracert,xno,xno,dt,ddepno)

c 
c NO2:
c
      do i=1,n3d
        tracert(i,1,1) = xno2(i,1,1)
      enddo

      call chemadv(8,
     &             tracert,xno2,xno2,dt,ddepno2)

c 
c N2O5:
c
      do i=1,n3d
        tracert(i,1,1) = xn2o5(i,1,1)
      enddo

      call chemadv(8,
     &             tracert,xn2o5,xn2o5,dt,ddepn2o5)

c 
c HNO3:
c
      do i=1,n3d
        tracert(i,1,1) = hno3(i,1,1)
      enddo

      call chemadv(8,
     &             tracert,hno3,hno3,dt,ddephno3)

c 
c CH4:
c
      do i=1,n3d
        tracert(i,1,1) = ch4(i,1,1)
      enddo

      call chemadv(8,
     &             tracert,ch4,ch4,dt,ddepref)

c 
c CH2O:
c
      do i=1,n3d
        tracert(i,1,1) = ch2o(i,1,1)
      enddo

      call chemadv(8,
     &             tracert,ch2o,ch2o,dt,ddepref)

c 
c SO2:
c
      do i=1,n3d
        tracert(i,1,1) = so2(i,1,1)
      enddo

      call chemadv(8,
     &             tracert,so2,so2,dt,ddepref)

c 
c H2SO4:
c
      do i=1,n3d
        tracert(i,1,1) = h2so4(i,1,1)
      enddo

      call chemadv(8,
     &             tracert,h2so4,h2so4,dt,ddepref)

#endif

	return
	 end


c     =====================================================
      Subroutine chemadv(ifdiff,
     &                   x00,x11,xinit,dt1,ddepspd)
c     =====================================================

c==================================================================c
c                                                                  c
c   CHEMADV.F:    Subroutine for calculating advection and eddy    c
c                   transport of MIT Global Chemistry Model        c
c        -------------------------------------------------         c
c   Author:		Chien Wang                                 c
c			MIT Joint Program on Science and Policy    c
c			    of Global Change                       c
c   Last Revised on:    January 30, 1996                           c
c                                                                  c
c==================================================================c

      parameter(xxx1 = 1./6., xxx2 = 4.0*xxx1)
      parameter(yyy3 = 1./36.,yyy2 =10.0*yyy3, yyy1 =25.0*yyy3)

#include "chem_para"
#include "chem_com"
#include "BD2G04.COM"

      common /WORK1/pit(nlon,nlat),sd(nlon,nlat,nlev1)

      dimension x00  (nlon,nlat,nlev)
      dimension x11  (nlon,nlat,nlev)
      dimension xinit(nlon,nlat,nlev)

c 062095 dry deposition speed, in sigma/second and positive
c        speed is updraft
c
      dimension ddepspd(nlon,nlat)

      dimension c(nlat+1),x(nlat),w1(4,2),w2(nlat,3),
     &         w4(nlat,5),ww(nlat+1,5),ww2(nlat+1,5)

!	---------------------------------------------------------

#if ( defined CPL_CHEM )

c-------------------------------------------------------
c   Definitions of parameters:
c
c   Basic time step for advection
c     dta =dt1*3.0      ! dta=1 hr.
      dta =dt1           ! dt1=20 min in GISS therefore here 
c
c 111596:
c
	dt1hr = 3600.

      dt2 =dta*0.5

      istart=1
      iend  =nlon

c-------------------------------
c  Start do loop for tracers:

      do 1000 ntime=1,ndyn

c-------------------------------------------------------
c   Scaling mixing ratio with PAI:
c
      do 1 k=1,nlev
      do 1 i=1,n2dh
        x00(i,1,k)=x00(i,1,k)*p00(i,1)
1     continue

      do 11 k=1,nlev
        pvv(1,1,k)=0.0
11    continue

      do 12 j=1,nlat
        pww(1,j,1)=0.0
12    continue

c-------------------------------------------------------
c   Calculating meridional advection:
c
	i = 1
       do 61 k=1,nlev
c       do 61 i=istart,iend

c       c(1)=0.0
       c(1)= c(2)
       do 62 j=2,nlat
         c(j)   =pvv(i,j,k)/dyv(j)*dta
62     continue
c       c(nlat+1)=0.0
       c(nlat+1)= c(nlat)
      
       call pdadv1(c,w4,w2,w1,nlat)
        
       do 63 j=1,nlat
         x(j)=x00(i,j,k)
63     continue

       call pdadv2(c,x,w4,w2,w1,ww,ww2,nlat,1)

c---------------------------
c  Lateral BC:
c
c  South pole:
c
        fluxl=pvv(i,2,k)*(x00(i,2,k)+x00(i,1,k))/dyv(2)
     &       *dta*0.5

        fluxl=max(-x00(i,2,k),
     &        min( x00(i,1,k),fluxl))

        fluxr=pvv(i,3,k)*(x00(i,3,k)+x00(i,2,k))/dyv(3)
     &       *dta*0.5

        fluxr=max(-x00(i,3,k),
     &        min( x00(i,2,k),fluxr))

        x00(i,2,k)=x00(i,2,k)
     &            -(fluxr-fluxl)

        fluxlbc = 
     &            -min(0.0,pvv(i,2,k))
     &            *(x11(i,2,k)-x11(i,1,k))/dyv(2)
     &            *(p00(i,2)+p00(i,1))*0.5
     &            *dta
        fluxlbc =
     &            max(-x00(i,1,k),
     &            min( x00(i,2,k),fluxlbc))

        x00(i,1,k)=x00(i,1,k)
     &            +fluxlbc

c
c  North pole:
c
        fluxl=pvv(i,nlat1,k)*(x00(i,nlat1,k)
     &                       +x00(i,nlat2,k))/dyv(nlat1)
     &       *dta*0.5

        fluxl=max(-x00(i,nlat1,k),
     &        min( x00(i,nlat2,k),fluxl))

        fluxr=pvv(i,nlat, k)*(x00(i,nlat, k)
     &                       +x00(i,nlat1,k))/dyv(nlat)
     &       *dta*0.5

        fluxr=
     &        max(-x00(i,nlat, k),
     &        min( x00(i,nlat1,k),fluxr))

        x00(i,nlat1,k)=x00(i,nlat1,k)
     &                -(fluxr-fluxl)

        fluxlbc =
     &                -max(0.0,pvv(i,nlat,k))
     &                *(x11(i,nlat,k)-x11(i,nlat1,k))/dyv(nlat)
     &                *(p00(i,nlat)+p00(i,nlat1))*0.5
     &                *dta
        fluxlbc =
     &                max(-x00(i,nlat,k),
     &                min( x00(i,nlat1,k),fluxlbc))

        x00(i,nlat,k) =x00(i,nlat,k)
     &                +fluxlbc

c---
c   Adjustment of momentum equation
c
c       do 64 j=2,nlat1
       do 64 j=3,nlat2
         if(k.ne.1
     & .and.k.ne.nlev
     &                )then
	 deltac = max(-1.0,
     &            min(+1.0,c(j+1)-c(j)))
         x00(i,j,k)=x(j)
     &             +x00(i,j,k)*deltac
         endif
64     continue

c 051595:

       if(  k.ne.1
     & .and.k.ne.nlev
     &                               )then
        x00(i,2,k)=x00(i,2,k)
     &            *(1.0+(pvv(i,3,k)/dyv(3)
     &                  -pvv(i,2,k)/dyv(2))*dta)
        x00(i,1,k)=x00(i,1,k)
     &            *(1.0+pvv(i,2,k)/dyv(2)
     &                                   *dta)
        x00(i,nlat1,k)=x00(i,nlat1,k)
     &                *(1.0+(pvv(i,nlat, k)/dyv(nlat)
     &                      -pvv(i,nlat1,k)/dyv(nlat1))*dta)
        x00(i,nlat,k) =x00(i,nlat,k)
     &                *(1.0-pvv(i,nlat,k)/dyv(nlat)
     &                                            *dta)
       endif
c=====

61     continue
       
c-------------------------------------------------------
c   Calculating vertical advection:
c
c       do 66 i=istart,iend
	i = 1
       do 66 j=1,nlat
c       do 66 j=2,nlat1

       c(1)     =0.0
       do 67 k=2,nlev1
         c(k)   =-pww(i,j,k)/dsig(k)*dta
67     continue
       c(nlev)  =-pww(i,j,nlev1)/dsig(nlev1)*dta
       c(nlev+1)=0.0
      
       call pdadv1(c,w4,w2,w1,nlev)
        
       do 68 k=1,nlev
         x(k)=x00(i,j,k)
68     continue

       call pdadv2(c,x,w4,w2,w1,ww,ww2,nlev,1)

c---
c  VBC:
c
         fluxt=pww(i,j,3)*(x00(i,j,3)+x00(i,j,2))
     &        /dsig(3)
     &        *dta*0.5

c 112596:
c         fluxt=-max(-x00(i,j,3),
c     &          min( x00(i,j,2),-fluxt))
         fluxt=-max(-x00(i,j,3)*0.5,
     &          min( x00(i,j,2)*0.5,-fluxt))

         fluxb=pww(i,j,2)*(x00(i,j,2)+x00(i,j,1))
     &        /dsig(2)
     &        *dta*0.5

c 112596:
c         fluxb=-max(-x00(i,j,2),
c     &          min( x00(i,j,1),-fluxb))
         fluxb=-max(-x00(i,j,2)*0.5,
     &          min( x00(i,j,1)*0.5,-fluxb))

         x00(i,j,2)=x00(i,j,2)
     &             +(fluxt-fluxb)

         x00(i,j,1)=x00(i,j,1)
cc     &             +fluxb
c     &             +pww(i,j,2)
c 062095 add dry deposition:
!     &             +max(0.0,pww(i,j,2)-ddepspd(i,j))
     &             +(pww(i,j,2)-ddepspd(i,j))
     &             *(x11(i,j,2)-x11(i,j,1))/dsig(1)
     &             *p00(i,j)
     &             *dta

         fluxt=pww(i,j,nlev)*(x00(i,j,nlev)
     &                      +x00(i,j,nlev1))
     &        /dsig(nlev)
     &        *dta*0.5

c 112596:
c         fluxt=-max(-x00(i,j,nlev),
c     &          min( x00(i,j,nlev1),-fluxt))
         fluxt=-max(-x00(i,j,nlev)*0.5,
     &          min( x00(i,j,nlev1)*0.5,-fluxt))

         fluxb=pww(i,j,nlev1)*(x00(i,j,nlev1)
     &                       +x00(i,j,nlev2))
     &        /dsig(nlev1)
     &        *dta*0.5

c 112596
c         fluxb=-max(-x00(i,j,nlev1),
c     &          min( x00(i,j,nlev2),-fluxb))
         fluxb=-max(-x00(i,j,nlev1)*0.5,
     &          min( x00(i,j,nlev2)*0.5,-fluxb))

         x00(i,j,nlev1)=x00(i,j,nlev1)
     &                 +(fluxt-fluxb)

         x00(i,j,nlev)=x00(i,j,nlev)
c     &                -fluxb
     &                +min(0.0,pww(i,j,nlev))
     &                *(x11(i,j,nlev)-x11(i,j,nlev1))/dsig(nlev)
     &                *p00(i,j)
     &                *dta

c---
c
c       do 69 k=2,nlev1
       do 69 k=3,nlev2
         if(j.ne.1.and.j.ne.nlat)then
	 deltac = max(-1.0,
     &            min(+1.0,c(k+1)-c(k)))
         x00(i,j,k)=x(k)
     &             +x00(i,j,k)*deltac
         endif
69     continue

c 051595:

       if(  j.ne.1.and.j.ne.nlat
c     & .and.j.ne.2.and.j.ne.nlat1
     &                            )then
c ===
c === 081295: set limitation of xyz
c ===
	xyz = (pww(i,j,3)/dsig(3)
     &        -pww(i,j,2)/dsig(2))*dta
	xyz = min( 1.0,
     &        max(-1.0 ,xyz))

         x00(i,j,2)=x00(i,j,2)
     &             *(1.0-xyz)

        xyz = (pww(i,j,nlev) /dsig(nlev)
     &        -pww(i,j,nlev1)/dsig(nlev1))*dta
	xyz = min( 1.0,
     &        max(-1.0 ,xyz))

         x00(i,j,nlev1)=x00(i,j,nlev1)
     &                 *(1.0-xyz)
       endif
c=====

66     continue

c       goto 2001
c9000   continue
       
c--------------------------------------------------------
c   Rescaling mixing ratio with PAI:
c

      do 200 k=1,nlev
      do 200 i=1,n2dh
        x00(i,1,k)=x00(i,1,k)/p11(i,1)

c 012797: limit error 
        deltazu   =1.2*x11(i,1,k)
        deltazl   =0.8*x11(i,1,k)
	x00(i,1,k) = max(deltazl, min(deltazu, x00(i,1,k)))

        deltax    =abs(x00(i,1,k)-x11(i,1,k))
        deltay    =1.e-10*x11(i,1,k)
        if(deltax.gt.deltay)then 
          x11(i,1,k)=x00(i,1,k)
        else
          x00(i,1,k)=x11(i,1,k)
        endif
200   continue

c----------------------------------------------------
c   Corner points:
c

      x00(1,1,1)       = xxx1*(x00(1,1,2) + x00(1,2,1))
     &                 + xxx2* x00(1,1,1)
      x11(1,1,1)       = x00(1,1,1)

      x00(1,nlat,1)    = xxx1*(x00(1,nlat,2) + x00(1,nlat1,1))
     &                 + xxx2* x00(1,nlat,1)
      x11(1,nlat,1)    = x00(1,nlat,1)

      x00(1,1,nlev)    = xxx1*(x00(1,1,nlev1) + x00(1,2,nlev))
     &                 + xxx2* x00(1,1,nlev)
      x11(1,1,nlev)    = x00(1,1,nlev)

      x00(1,nlat,nlev) = xxx1*(x00(1,nlat,nlev1) + x00(1,nlat1,nlev))
     &                 + xxx2* x00(1,nlat,nlev)
      x11(1,nlat,nlev) = x00(1,nlat,nlev)

c-----------------------------------------------------
c  LBC smoothing:
c

      do k=1,nlev
c        x00(1,2,k)     = xxx1*(x00(1,1,k) + x00(1,3,k))
c     &                 + xxx2* x00(1,2,k)

c        x00(1,nlat1,k) = xxx1*(x00(1,nlat2,k) + x00(1,nlat,k))
c     &                 + xxx2* x00(1,nlat1,k)

        x00(1,1,k)     = yyy1*x00(1,1,k) + yyy2*x00(1,2,k)
     &                 + yyy3*x00(1,3,k)

        x00(1,nlat,k)  = yyy1*x00(1,nlat,k) + yyy2*x00(1,nlat1,k)
     &                 + yyy3*x00(1,nlat2,k)

        x11(1,1,k)     = x00(1,1,k)
c        x11(1,2,k)     = x00(1,2,k)
        x11(1,nlat,k)  = x00(1,nlat,k)
c        x11(1,nlat1,k) = x00(1,nlat1,k)

      enddo

c------------------------
c  Artificial diffusion
c    for top two levels:
c
      atfk=1.e-3

      i=1
      k=nlev1
      do 606 j=1,nlat
        x00(i,j,k)=x00(i,j,k)
     &            +atfk
     &           *(x00(i,j,k+1)+x00(i,j,k-1)-2.0*x00(i,j,k))
        x11(i,j,k)=x00(i,j,k)
606   continue

      k=nlev
      do 607 j=1,nlat
        x00(i,j,k)=x00(i,j,k)
     &            +atfk
     &           *min(0.0,(x00(i,j,k-1)-x00(i,j,k)))
        x11(i,j,k)=x00(i,j,k)
607   continue

c-------------------------
c  Artificial diffusion
c    for bottom two levels:
c
      atfk = 1.e-6	!1.e-5	!1.e-4

c      k=2
c      do 616 j=1,nlat
c        x00(i,j,k)=x00(i,j,k)
c     &            +atfk
c     &           *(x00(i,j,k+1)+x00(i,j,k-1)-2.0*x00(i,j,k))
c        x11(i,j,k)=x00(i,j,k)
c616   continue
c
c      k=1
c      do 617 j=1,nlat
c        x00(i,j,k)=x00(i,j,k)
c     &            +atfk
c     &           *(x00(i,j,k+1)-x00(i,j,k))
c        x11(i,j,k)=x00(i,j,k)
c617   continue

1000  continue

c-----------------------
c   Calculate eddy diffusion
c     and re-scale mass:
c
      do 501 k=1,nlev
      do 501 i=1,n2dh
        x00(i,1,k)=x00(i,1,k)*p00(i,1)
501   continue

c 111596:
c
       dteddy = dt1hr            !GISS model calculate eddy diffusion
                    	         !at every 1 hr in new version

      if(meddy1.eq.1) call chemeddy(ifdiff,x00,x11,dteddy)

      do 502 k=1,nlev
      do 502 i=1,n2dh
          x00(i,1,k)=x00(i,1,k)
     &              +x11(i,1,k)*(p11(i,1)-p00(i,1))
        x00(i,1,k)=x00(i,1,k)/p11(i,1)

c 012797: limit error
        deltazu   =1.2*x11(i,1,k)
        deltazl   =0.8*x11(i,1,k)
        x00(i,1,k) = max(deltazl, min(deltazu, x00(i,1,k)))

        x11(i,1,k)=x00(i,1,k)
502   continue

c
c-----------------------

#endif

      return
       end


c     ======================================================
      Subroutine chemadv2(ifdiff,
     &                    x00,x11,xinit,dt1,ddepspd)
c     ======================================================

c==================================================================c
c                                                                  c
c   CHEMADV.F:    Subroutine for calculating advection and eddy    c
c                   transport of MIT Global Chemistry Model        c
c        -------------------------------------------------         c
c   Author:		Chien Wang                                 c
c			MIT Joint Program on Science and Policy    c
c			    of Global Change                       c
c   Last Revised on:    January 30, 1996                           c
c                                                                  c
c==================================================================c

      parameter(xxx1 = 1./6., xxx2 = 4.0*xxx1)
      parameter(yyy3 = 1./36.,yyy2 =10.0*yyy3, yyy1 =25.0*yyy3)

#include "chem_para"
#include "chem_com"
#include "BD2G04.COM"

      common /WORK1/pit(nlon,nlat),sd(nlon,nlat,nlev1)

      dimension x00  (nlon,nlat,nlev)
      dimension x11  (nlon,nlat,nlev)
      dimension xinit(nlon,nlat,nlev)

c 062095 dry deposition speed, in sigma/second and positive
c        speed is updraft
c
      dimension ddepspd(nlon,nlat)

      dimension c(nlat+1),x(nlat),w1(4,2),w2(nlat,3),
     &         w4(nlat,5),ww(nlat+1,5),ww2(nlat+1,5)

!	---------------------------------------------------------

#if ( defined CPL_CHEM )

c-------------------------------------------------------
c   Definitions of parameters:
c
c   Basic time step for advection
c     dta =dt1*3.0      ! dta=1 hr.
      dta =dt1           ! dt1=20 min in GISS therefore here 
c
c 111596:
c
        dt1hr = 3600.

      dt2 =dta*0.5

      istart=1
      iend  =nlon

c-------------------------------
c  Start do loop for tracers:

      do 1000 ntime=1,ndyn

c-------------------------------------------------------
c   Scaling mixing ratio with PAI:
c
      do 1 k=1,nlev
      do 1 i=1,n2dh
        x00(i,1,k)=x00(i,1,k)*p00(i,1)
1     continue

      do 11 k=1,nlev
        pvv(1,1,k)=0.0
11    continue

      do 12 j=1,nlat
        pww(1,j,1)=0.0
12    continue

c-------------------------------------------------------
c   Calculating meridional advection:
c
	i = 1
       do 61 k=1,nlev
c       do 61 i=istart,iend

c       c(1)=0.0
       c(1)= c(2)
       do 62 j=2,nlat
         c(j)   =pvv(i,j,k)/dyv(j)*dta
62     continue
c       c(nlat+1)=0.0
       c(nlat+1)= c(nlat)
      
       call pdadv1(c,w4,w2,w1,nlat)
        
       do 63 j=1,nlat
         x(j)=x00(i,j,k)
63     continue

       call pdadv2(c,x,w4,w2,w1,ww,ww2,nlat,1)

c---------------------------
c  Lateral BC:
c
c  South pole:
c
        fluxl=pvv(i,2,k)*(x00(i,2,k)+x00(i,1,k))/dyv(2)
     &       *dta*0.5

        fluxl=max(-x00(i,2,k),
     &        min( x00(i,1,k),fluxl))

        fluxr=pvv(i,3,k)*(x00(i,3,k)+x00(i,2,k))/dyv(3)
     &       *dta*0.5

        fluxr=max(-x00(i,3,k),
     &        min( x00(i,2,k),fluxr))

        x00(i,2,k)=x00(i,2,k)
     &            -(fluxr-fluxl)

        fluxlbc = 
     &            -min(0.0,pvv(i,2,k))
     &            *(x11(i,2,k)-x11(i,1,k))/dyv(2)
     &            *(p00(i,2)+p00(i,1))*0.5
     &            *dta
        fluxlbc =
     &            max(-x00(i,1,k),
     &            min( x00(i,2,k),fluxlbc))

        x00(i,1,k)=x00(i,1,k)
     &            +fluxlbc

c
c  North pole:
c
        fluxl=pvv(i,nlat1,k)*(x00(i,nlat1,k)
     &                       +x00(i,nlat2,k))/dyv(nlat1)
     &       *dta*0.5

        fluxl=max(-x00(i,nlat1,k),
     &        min( x00(i,nlat2,k),fluxl))

        fluxr=pvv(i,nlat, k)*(x00(i,nlat, k)
     &                       +x00(i,nlat1,k))/dyv(nlat)
     &       *dta*0.5

        fluxr=
     &        max(-x00(i,nlat, k),
     &        min( x00(i,nlat1,k),fluxr))

        x00(i,nlat1,k)=x00(i,nlat1,k)
     &                -(fluxr-fluxl)

        fluxlbc =
     &                -max(0.0,pvv(i,nlat,k))
     &                *(x11(i,nlat,k)-x11(i,nlat1,k))/dyv(nlat)
     &                *(p00(i,nlat)+p00(i,nlat1))*0.5
     &                *dta
        fluxlbc =
     &                max(-x00(i,nlat,k),
     &                min( x00(i,nlat1,k),fluxlbc))

        x00(i,nlat,k) =x00(i,nlat,k)
     &                +fluxlbc

c---
c   Adjustment of momentum equation
c
c       do 64 j=2,nlat1
       do 64 j=3,nlat2
         if(k.ne.1
     & .and.k.ne.nlev
     &                )then
	 deltac = max(-1.0,
     &            min(+1.0,c(j+1)-c(j)))
         x00(i,j,k)=x(j)
     &             +x00(i,j,k)*deltac
         endif
64     continue

c 051595:

       if(  k.ne.1
     & .and.k.ne.nlev
     &                               )then
        x00(i,2,k)=x00(i,2,k)
     &            *(1.0+(pvv(i,3,k)/dyv(3)
     &                  -pvv(i,2,k)/dyv(2))*dta)
        x00(i,1,k)=x00(i,1,k)
     &            *(1.0+pvv(i,2,k)/dyv(2)
     &                                   *dta)
        x00(i,nlat1,k)=x00(i,nlat1,k)
     &                *(1.0+(pvv(i,nlat, k)/dyv(nlat)
     &                      -pvv(i,nlat1,k)/dyv(nlat1))*dta)
        x00(i,nlat,k) =x00(i,nlat,k)
     &                *(1.0-pvv(i,nlat,k)/dyv(nlat)
     &                                            *dta)
       endif
c=====

61     continue
       
c-------------------------------------------------------
c   Calculating vertical advection:
c
c       do 66 i=istart,iend
	i = 1
       do 66 j=1,nlat
c       do 66 j=2,nlat1

       c(1)     =0.0
       do 67 k=2,nlev1
         c(k)   =-pww(i,j,k)/dsig(k)*dta
#ifdef O3_TOP_SPARC
         if (k.ge.n_tropopause) then
           c(k)   = c(k)*0.25
!          c(k)   = c(k)*0.5
!          c(k)   = 0.0
         endif 
#endif
67     continue
       c(nlev)  =-pww(i,j,nlev1)/dsig(nlev1)*dta
       c(nlev+1)=0.0
      
       call pdadv1(c,w4,w2,w1,nlev)
        
       do 68 k=1,nlev
         x(k)=x00(i,j,k)
68     continue

       call pdadv2(c,x,w4,w2,w1,ww,ww2,nlev,1)

c---
c  VBC:
c
         fluxt=pww(i,j,3)*(x00(i,j,3)+x00(i,j,2))
     &        /dsig(3)
     &        *dta*0.5

c 112796
c         fluxt=-max(-x00(i,j,3),
c     &          min( x00(i,j,2),-fluxt))
         fluxt=-max(-x00(i,j,3)*0.5,
     &          min( x00(i,j,2)*0.5,-fluxt))

         fluxb=pww(i,j,2)*(x00(i,j,2)+x00(i,j,1))
     &        /dsig(2)
     &        *dta*0.5

c 112796
c         fluxb=-max(-x00(i,j,2),
c     &          min( x00(i,j,1),-fluxb))
         fluxb=-max(-x00(i,j,2)*0.5,
     &          min( x00(i,j,1)*0.5,-fluxb))

         x00(i,j,2)=x00(i,j,2)
     &             +(fluxt-fluxb)

         x00(i,j,1)=x00(i,j,1)
cc     &             +fluxb
c     &             +pww(i,j,2)
c 062095 add dry deposition:
!     &             +max(0.0,pww(i,j,2)-ddepspd(i,j))
     &             +(pww(i,j,2)-ddepspd(i,j))
     &             *(x11(i,j,2)-x11(i,j,1))/dsig(1)
     &             *p00(i,j)
     &             *dta

         fluxt=pww(i,j,nlev)*(x00(i,j,nlev)
     &                      +x00(i,j,nlev1))
     &        /dsig(nlev)
     &        *dta*0.5

c 112796
c         fluxt=-max(-x00(i,j,nlev),
c     &          min( x00(i,j,nlev1),-fluxt))
         fluxt=-max(-x00(i,j,nlev)*0.5,
     &          min( x00(i,j,nlev1)*0.5,-fluxt))

         fluxb=pww(i,j,nlev1)*(x00(i,j,nlev1)
     &                       +x00(i,j,nlev2))
     &        /dsig(nlev1)
     &        *dta*0.5

c 112796
c         fluxb=-max(-x00(i,j,nlev1),
c     &          min( x00(i,j,nlev2),-fluxb))
         fluxb=-max(-x00(i,j,nlev1)*0.5,
     &          min( x00(i,j,nlev2)*0.5,-fluxb))

         x00(i,j,nlev1)=x00(i,j,nlev1)
     &                 +(fluxt-fluxb)

         x00(i,j,nlev)=x00(i,j,nlev)
c     &                -fluxb
     &                +min(0.0,pww(i,j,nlev))
     &                *(x11(i,j,nlev)-x11(i,j,nlev1))/dsig(nlev)
     &                *p00(i,j)
     &                *dta


c---
c
c       do 69 k=2,nlev1
       do 69 k=3,nlev2
         if(j.ne.1.and.j.ne.nlat)then
	 deltac = max(-1.0,
     &            min(+1.0,c(k+1)-c(k)))
         x00(i,j,k)=x(k)
     &             +x00(i,j,k)*deltac
         endif
69     continue

c 051595:

       if(  j.ne.1.and.j.ne.nlat
c     & .and.j.ne.2.and.j.ne.nlat1
     &                            )then
c ===
c === 081295: set limitation of xyz
c ===
	xyz = (pww(i,j,3)/dsig(3)
     &        -pww(i,j,2)/dsig(2))*dta
	xyz = min( 1.0,
     &        max(-1.0 ,xyz))

         x00(i,j,2)=x00(i,j,2)
     &             *(1.0-xyz)

        xyz = (pww(i,j,nlev) /dsig(nlev)
     &        -pww(i,j,nlev1)/dsig(nlev1))*dta
	xyz = min( 1.0,
     &        max(-1.0 ,xyz))

         x00(i,j,nlev1)=x00(i,j,nlev1)
     &                 *(1.0-xyz)
       endif
c=====

66     continue

c       goto 2001
c9000   continue
       
c--------------------------------------------------------
c   Rescaling mixing ratio with PAI:
c

      do 200 k=1,nlev
      do 200 i=1,n2dh
        x00(i,1,k)=x00(i,1,k)/p11(i,1)

c 012797: limit error
        deltazu   =1.2*x11(i,1,k)
        deltazl   =0.8*x11(i,1,k)
        x00(i,1,k) = max(deltazl, min(deltazu, x00(i,1,k)))

        deltax    =abs(x00(i,1,k)-x11(i,1,k))
        deltay    =1.e-10*x11(i,1,k)
        if(deltax.gt.deltay)then 
          x11(i,1,k)=x00(i,1,k)
        else
          x00(i,1,k)=x11(i,1,k)
        endif
200   continue

c----------------------------------------------------
c   Corner points:
c

      x00(1,1,1)       = xxx1*(x00(1,1,2) + x00(1,2,1))
     &                 + xxx2* x00(1,1,1)
      x11(1,1,1)       = x00(1,1,1)

      x00(1,nlat,1)    = xxx1*(x00(1,nlat,2) + x00(1,nlat1,1))
     &                 + xxx2* x00(1,nlat,1)
      x11(1,nlat,1)    = x00(1,nlat,1)

      x00(1,1,nlev)    = xxx1*(x00(1,1,nlev1) + x00(1,2,nlev))
     &                 + xxx2* x00(1,1,nlev)
      x11(1,1,nlev)    = x00(1,1,nlev)

      x00(1,nlat,nlev) = xxx1*(x00(1,nlat,nlev1) + x00(1,nlat1,nlev))
     &                 + xxx2* x00(1,nlat,nlev)
      x11(1,nlat,nlev) = x00(1,nlat,nlev)

c-----------------------------------------------------
c  LBC smoothing:
c

      do k=1,nlev
c        x00(1,2,k)     = xxx1*(x00(1,1,k) + x00(1,3,k))
c     &                 + xxx2* x00(1,2,k)

c        x00(1,nlat1,k) = xxx1*(x00(1,nlat2,k) + x00(1,nlat,k))
c     &                 + xxx2* x00(1,nlat1,k)

        x00(1,1,k)     = yyy1*x00(1,1,k) + yyy2*x00(1,2,k)
     &                 + yyy3*x00(1,3,k)

        x00(1,nlat,k)  = yyy1*x00(1,nlat,k) + yyy2*x00(1,nlat1,k)
     &                 + yyy3*x00(1,nlat2,k)

        x11(1,1,k)     = x00(1,1,k)
c        x11(1,2,k)     = x00(1,2,k)
        x11(1,nlat,k)  = x00(1,nlat,k)
c        x11(1,nlat1,k) = x00(1,nlat1,k)

      enddo

c------------------------
c  Artificial diffusion
c    for top two levels:
c
c ===== 091495:
c      atfk=1.e-4       !138 d
c     atfk=2.e-5        !690 d
      atfk= 1.e-10	!1.e-6 better - 031696    !20*690 d

      i=1
      k=nlev1
      do 606 j=1,nlat
        x00(i,j,k)=x00(i,j,k)
     &            +atfk
     &           *(x00(i,j,k+1)+x00(i,j,k-1)-2.0*x00(i,j,k))
        x11(i,j,k)=x00(i,j,k)
606   continue

c      k=nlev
c      do 607 j=1,nlat
c        x00(i,j,k)=x00(i,j,k)
c     &            +atfk
c     &           *min(0.0,(x00(i,j,k-1)-x00(i,j,k)))
c        x11(i,j,k)=x00(i,j,k)
c607   continue

c-------------------------
c  Artificial diffusion
c    for bottom two levels:
c
      atfk = 1.e-6	!1.e-5	!1.e-4

c      k=2
c      do 616 j=1,nlat
c        x00(i,j,k)=x00(i,j,k)
c     &            +atfk
c     &           *(x00(i,j,k+1)+x00(i,j,k-1)-2.0*x00(i,j,k))
c        x11(i,j,k)=x00(i,j,k)
c616   continue
c
c      k=1
c      do 617 j=1,nlat
c        x00(i,j,k)=x00(i,j,k)
c     &            +atfk
c     &           *(x00(i,j,k+1)-x00(i,j,k))
c        x11(i,j,k)=x00(i,j,k)
c617   continue

1000  continue

c-----------------------
c   Calculate eddy diffusion
c     and re-scale mass:
c
      do 501 k=1,nlev
      do 501 i=1,n2dh
        x00(i,1,k)=x00(i,1,k)*p00(i,1)
501   continue

c 111596:
c
       dteddy = dt1hr            !GISS model calculate eddy diffusion
                                 !at every 1 hr in new version

      if(meddy1.eq.1) call chemeddy(ifdiff,x00,x11,dteddy)

      do 502 k=1,nlev
      do 502 i=1,n2dh
          x00(i,1,k)=x00(i,1,k)
     &              +x11(i,1,k)*(p11(i,1)-p00(i,1))
        x00(i,1,k)=x00(i,1,k)/p11(i,1)

c 012797: limit error
        deltazu   =1.2*x11(i,1,k)
        deltazl   =0.8*x11(i,1,k)
        x00(i,1,k) = max(deltazl, min(deltazu, x00(i,1,k)))

        x11(i,1,k)=x00(i,1,k)
502   continue

c
c-----------------------
#endif

      return
       end

