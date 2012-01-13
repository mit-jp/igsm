#include "ctrparam.h"

!  Subroutine for interpolating ACCRI emissions
      subroutine int_accri_emi(j,psrf,xxaccri,emi_accri,month,iiyear)
#include "chem_para"
#include "chem_com"
#include "BD2G04.COM"
#if ( defined ACCRI )

      dimension  emi_accri(nlev),xxaccri(nlon,nlat,nlev_accri,12)
      dimension  phb(0:nlev),dp_accri(nlev_accri)
      dimension ndmonth(12)
      logical first
      real totaccri,totigsm
      integer ncall
      save ncall,dp_accri
      data first /.true./
      DATA ndmonth/31,28,31,30,31,30,31,31,30,31,30,31/
      if(first) then
      do l=1,nlev_accri
        dp_accri(l)=p_accri(l-1)-p_accri(l)
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
       totaccri=totaccri+xxaccri(1,j,l,month)
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
         emi_accri(l)=emi_accri(l)+xxaccri(1,j,lb,month)
         lb=lb+1
         if(lb.lt.nlev_accri)go to 333
       else
         dp=p_accri(lb-1)-phb(l)
         emi=xxaccri(1,j,lb,month)*dp/dp_accri(lb)
         emi_accri(l)=emi_accri(l)+emi
         emi0=xxaccri(1,j,lb,month)-emi
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
       emi_accri(l)=scale_accri(iiyear)*emi_accri(l)
     &   /ndmonth(month)
       if(J.eq.-23) then
       print *,emi_accri(l)
       endif
      enddo
      
      if(totaccri.gt.1e-3) then
!     if( abs(totaccri-totigsm).ge.1e-3) then
      if( abs(totaccri-totigsm)/totaccri.ge.1e-3) then
         print *,'error in accri emissions interpolation'
         print *,j,month,psrf,ncall
         print *,totaccri,totigsm,totaccri-totigsm
         print *,p_accri
         print *,phb
      do l=1,nlev
       print *,l,emi_accri(l)/scale_accri(iiyear)
     &   *ndmonth(month)
      enddo
      do l=1,nlev_accri
       print *,l,xxaccri(1,j,l,month)
      enddo
         stop
       endif
       endif
      ncall=ncall+1
#endif
      return
      end
