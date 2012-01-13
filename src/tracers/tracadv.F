
#include "ctrparam.h"
!       =======================
        subroutine tracadv0 (dt)
!       =======================


#include "chem_para"
#include "chem_com"
#include "tracer_com.h"

        dimension tracer0(nlon,nlat,nlev)
        dimension tracer1(nlon,nlat,nlev)
        dimension ddepspd(nlon,nlat)
          
!        print *,'DT from tracadv0=',dt
!        print *,'PVV from tracadv0='
!        print *,pvv
!        print *,'PWW from tracadv0='
!        print *,pww


         ifdiff = 8
         do ntr=1,ntracers
            tracer0 = tracers (:,:,:,ntr)
            tracer1 = tracers (:,:,:,ntr)
!           print *,' TRACER1 ', ntr
!           print *,tracer1
            if (trtype(ntr).eq.3.or.trtype(ntr).eq.4) then
              ddepspd = -0.1150E-06
            else
              ddepspd = 0.0
           endif

           call tracadv(ifdiff,tracer0,tracer1,dt,ddepspd)

            tracers (:,:,:,ntr) = tracer1
         enddo
       end

c     =====================================================
      Subroutine tracadv(ifdiff,
     &                   x00,x11,dt1,ddepspd)
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

c 062095 dry deposition speed, in sigma/second and positive
c        speed is updraft
c
      dimension ddepspd(nlon,nlat)

      dimension c(nlat+1),x(nlat),w1(4,2),w2(nlat,3),
     &         w4(nlat,5),ww(nlat+1,5),ww2(nlat+1,5)

!	---------------------------------------------------------


c-------------------------------------------------------
c   Definitions of parameters:
c
c   Basic time step for advection
      dta =dt1           ! dt1 dynamical time step
c                        ! 20 min for JM=24 and 10 for JM=46
c 111596:
c
	dt1hr = 3600.

!        print *,'DT from tracadv=',dt1
!        print *,'ndyn from tracadv=',ndyn
!     print *,'X11 before  anything'
!     print *,x11

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

       do 62 j=2,nlat
         c(j)   =pvv(i,j,k)/dyv(j)*dta
62     continue
c       c(1)=0.0
       c(1)= c(2)
c       c(nlat+1)=0.0
       c(nlat+1)= c(nlat)
      
!      print *,'FROM meridional advection'
!      print *,c
       call pdadvt1(c,w4,w2,w1,nlat)
        
       do 63 j=1,nlat
         x(j)=x00(i,j,k)
63     continue

       call pdadvt2(c,x,w4,w2,w1,ww,ww2,nlat,1)

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
!      print *,'FROM vertical advection'
!      print *,c
      
       call pdadvt1(c,w4,w2,w1,nlev)
        
       do 68 k=1,nlev
         x(k)=x00(i,j,k)
68     continue

       call pdadvt2(c,x,w4,w2,w1,ww,ww2,nlev,1)

c---
c  VBC:
c
         fluxt=pww(i,j,3)*(x00(i,j,3)+x00(i,j,2))
     &        /dsig(3)
     &        *dta*0.5

         fluxt=-max(-x00(i,j,3)*0.5,
     &          min( x00(i,j,2)*0.5,-fluxt))

         fluxb=pww(i,j,2)*(x00(i,j,2)+x00(i,j,1))
     &        /dsig(2)
     &        *dta*0.5

         fluxb=-max(-x00(i,j,2)*0.5,
     &          min( x00(i,j,1)*0.5,-fluxb))

         x00(i,j,2)=x00(i,j,2)
     &             +(fluxt-fluxb)

         x00(i,j,1)=x00(i,j,1)
cc     &             +fluxb
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
        if (x11(i,1,k).gt.0.0) then
        deltazu   =1.2*x11(i,1,k)
        deltazl   =0.8*x11(i,1,k)
	x00(i,1,k) = max(deltazl, min(deltazu, x00(i,1,k)))
        endif

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

!     print *,'X11 before  LBC smoothing'
!     print *,x11
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
!     print *,'X11 before eddy'
!     print *,x11

!     if(meddy1.eq.1) call traceddy(ifdiff,x00,x11,dteddy)
      call traceddy(ifdiff,x00,x11,dteddy)

      do 502 k=1,nlev
      do 502 i=1,n2dh
          x00(i,1,k)=x00(i,1,k)
     &              +x11(i,1,k)*(p11(i,1)-p00(i,1))
        x00(i,1,k)=x00(i,1,k)/p11(i,1)

c 012797: limit error
        if (x11(i,1,k).gt.0.0) then
        deltazu   =1.2*x11(i,1,k)
        deltazl   =0.8*x11(i,1,k)
        x00(i,1,k) = max(deltazl, min(deltazu, x00(i,1,k)))
        endif

        x11(i,1,k)=x00(i,1,k)
502   continue

c
c-----------------------


      return
       end


c     ======================================================
