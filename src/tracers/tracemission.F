
#include "ctrparam.h"

!	============================================================
!
!	TRASEMISSION.F:  Subroutine for calculating surface emission
!				of tracers of MIT Global Chemistry Model
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
!
!	==========================================================

!
	subroutine tracemission 
!	===================================

#include "chem_para"
#include "chem_com"
#include "BD2G04.COM"
#include "tracer_com.h"


        COMMON U,V,T,P,Q
        dimension  emi_accri(nlev),x11(nlev),emi_trac(nlev_accri)
	
	real emi_new

!	----------------------------------------------------------
!
!        print *,'First year of emissions from chememi=',myyear


c
       i=1
       do  j=1,nlat


        tmass    = 0.0                !Total air mass in kg
        tmassx   = 0.0
        do k=1,nlev
          tmass    = tmass + airmass(i,j,k)
        enddo



      do n=1,ntracers

c
      do k=1,nlev
        x11(k) = tracers (i,j,k,n)
      enddo
      do k=1,nlev_accri
        emi_trac(k)=tracemis (i,j,k,n)
      enddo
       if(j.eq. -39.and.n.eq.3)then
          print *,' From TRACEMISS'
          print *,' PS=',ps4accri(i,j)
          print *,' N=',N
         endif
       call int_trac_emi(j,ps4accri(i,j),emi_trac,emi_accri)

       if(j.eq.-39.and.n.eq.3)then
      do k=1,nlev_accri
        print*,k,tracemis (i,j,k,n)
      enddo
      print *,' '
      print *,' '
      do k=1,nlev
          print *,k,emi_accri(k)
      enddo

         endif
      call tracemissionA(1,j,nlev,
     &                   emi_accri,x11)

      do k=1,nlev
         tracers (i,j,k,n) = x11(k)
      enddo

      enddo ! ntracers
      enddo ! J


      return
       end

c=====================================================

      subroutine tracemissionA(i,j,ktop,
     &                         emi_accri,x11)
c     ===========================================

c======================================================c
c  A subroutine for adding ACCRI emission               c
c     into layers from k=1 to ktop (nlev)    c
c======================================================c

#include "chem_para"
#include "chem_com"

      dimension x11(ktop)
      dimension emi_accri(ktop)

!     print *,'  '
!       xmass    	!Total tracer mass in 10^-9 kg
      do k=1,ktop
        xmass    = 
     &             x11(k)
        xemi     = ( xmass 
     &           + emi_accri(k) / airmass(i,j,k) )	! ppbm

!     print *,(x11(k)-xemi),emi_accri(k)

        x11(k)  = xemi
      enddo
!     print *,'  '

      return
       end

c=====================================================
!  Subroutine for interpolating ACCRI emissions
      subroutine int_trac_emi(j,psrf,xxaccri,emi_accri)
#include "chem_para"
#include "chem_com"
#include "BD2G04.COM"

      dimension  emi_accri(nlev),xxaccri(nlev_accri)
      dimension  phb(0:nlev),dp_accri(nlev_accri)
      dimension ndmonth(12)
      logical first
      real totaccri,totigsm
      integer ncall
      save ncall,dp_accri
      data first /.true./
      if(first) then
      print *,'DP_ACCRI'
      do l=1,nlev_accri
        dp_accri(l)=p_accri(l-1)-p_accri(l)
        print *,l,dp_accri(l)
      enddo
      ncall=1
      first = .false.
      endif
!     if(J.eq.23) then
!      print *,'From int_accri_emi'
!      print *,month,iiyear,psrf
!     endif
      totaccri=0.0
      totigsm=0.0
      do l=1,nlev_accri
       totaccri=totaccri+xxaccri(l)
      enddo
      do l=0,nlev
       phb(l)=psrf*SIGE(l+1)*1013./984.
      enddo
      lb=1
      p0=1013
      emi0=0.0
      do l=1,nlev
       emi_accri(l)=emi0
       emi0=0.0
  333  continue
       if(p_accri(lb).ge.phb(l)) then
         emi_accri(l)=emi_accri(l)+xxaccri(lb)
         lb=lb+1
         if(lb.lt.nlev_accri)go to 333
       else
         dp=p_accri(lb-1)-phb(l)
         emi=xxaccri(lb)*dp/dp_accri(lb)
         emi_accri(l)=emi_accri(l)+emi
         emi0=xxaccri(lb)-emi
         lb=lb+1
       endif
       if(lb.gt.nlev_accri)go to 334
      enddo
  334  continue
!      if(J.eq.23) then
!      print *,'scale_accri=',scale_accri(iiyear)
!      endif
      do l=1,nlev
      totigsm=totigsm+emi_accri(l)
       if(J.eq.-23) then
       print *,emi_accri(l)
       endif
      enddo
       if(J.eq.-39) then
       print *,totigsm,totaccri
       endif
      
      if(totaccri.gt.1e-3) then
!     if( abs(totaccri-totigsm).ge.1e-3) then
      if( abs(totaccri-totigsm)/totaccri.ge.1e-3) then
         print *,'error in accri emissions interpolation'
         print *,j,month,psrf,ncall
         print *,totaccri,totigsm,totaccri-totigsm
         print *,p_accri
         print *,phb
      do l=1,nlev
       print *,l,emi_accri(l)
      enddo
      do l=1,nlev_accri
       print *,l,xxaccri(l)
      enddo
         stop
       endif
       endif
      ncall=ncall+1
      return
      end


