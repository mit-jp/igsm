
#include "ctrparam.h"
!       =======================
        subroutine tracdrydep0 
!       =======================


#include "chem_para"
#include "chem_com"
#include "tracer_com.h"
#include "BD2G04.COM"
        COMMON U,V,T,P,Q

        dimension tracer0(nlon,nlat,nlev)
        dimension tracer1(nlon,nlat,nlev)
        dimension ddepspd(nlon,nlat)
          
!        print *,'From tracdrydep0='
!        print *,'PVV from tracadv0='
!        print *,pvv
!        print *,'PWW from tracadv0='
!        print *,pww


         do ntr=1,ntracers
            if (trtype(ntr).eq.3.or.trtype(ntr).eq.4) then
            tracer0 = tracers (:,:,:,ntr)
            tracer1 = tracers (:,:,:,ntr)
!           print *,' TRACER1 ', ntr
!           print *,tracer1

           call tracdrydep(tracer0,tracer1,dt,ddepspd)

            tracers (:,:,:,ntr) = tracer1
           endif
         enddo
!        stop
       end

c     =====================================================
      Subroutine tracdrydep(x00,x11,dt1,ddepspd)
c     =====================================================

c==================================================================c
c                                                                  c
c   CHEMADV.F:    Subroutine for calculating dry deposition        c
c                    of tracers MIT Global  Model                  c
c        -------------------------------------------------         c
c   Last Revised on:    November 14, 2011                           c
c                                                                  c
c==================================================================c


#include "chem_para"
#include "chem_com"
#include "BD2G04.COM"
        COMMON U,V,T,P,Q

      common /WORK1/pit(nlon,nlat),sd(nlon,nlat,nlev1)

      dimension x00  (nlon,nlat,nlev)
      dimension x11  (nlon,nlat,nlev)
      dimension pl(nlev),pk(nlev),tl(nlev)
      dimension hge(0:nlev),hgm(nlev),ro(nlev)
      dimension tausec(nlev)

c 062095 dry deposition speed, in sigma/second and positive
c        speed is updraft
c

!	---------------------------------------------------------


c-------------------------------------------------------
c   Definitions of parameters:
c
c   Basic time step for advection
      dta =dt1           ! dt1 dynamical time step
c 111596:
c
	dt1hr = 3600.
        wdep=1.e-3	! m/sec

!        print *,'DT from tracadv=',dt1
!        print *,'ndyn from tracadv=',ndyn
!     print *,'X11 before  anything'
!     print *,x11


c-------------------------------
c  Start do loop for tracers:

      i=1
      do j=1,nlat
!     print *,'j from drydep ',j
      DO 340 L=1,nlev       
      PL(L)=SIG(L)*P(I,J)+PTOP  
C**** TEMPERATURES                   
      PK(L)=EXPBYK(SIG(L)*P(I,J)+PTOP)
      TL(L)=T(I,J,L)*PK(L)
  340 CONTINUE  
c
      hge(0)=0.0
      do  l=1,nlev
          hge(l)=hge(l-1)+TL(L)*RGAS/sig(l)
     &    *dsig(l)/9.81
      enddo
      do  l=1,nlev
          hgm(l)=0.5*(hge(l)+hge(l-1))
          ro(l)=pl(l)/(TL(L)*RGAS)*100.
          tausec(l)=hgm(l)/wdep
!     print '(i3,3f10.1)',l,tl(l),hgm(l),tausec(l)/dt1hr/24
      enddo
      do  l=1,nlev
        x11(i,j,l)=x00(i,j,l)*(1.-1./(tausec(l)/dt1hr))
      enddo
      enddo	! j
c-----------------------


      return
       end


c     ======================================================
