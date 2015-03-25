#include "ctrparam.h"
! 	This subroutine is called one a month by climate2tem
!	to adjust vegfraction given to TEM to conserve land area
	subroutine vegfracadjustment
#include "BD2G04.COM"
#include "TEM.h"        
#include "CLM.h"        
!
#if ( defined CLM35 )
        logical first
!
        integer vegarea_igsm(jm0),vegarea_tem(jm0)
        logical vegsurf
        vegsurf(n)=((n.gt.1.and.n.le.17).or.n.eq.34)
!       vegsurf(n)=((n.gt.0.and.n.le.17).or.n.eq.34)
!
        data first /.true./
        if (first) then
!
!          print *, ' LAND AREAS'
          do j=1,jm0
           vegarea_tem(j)=0
           do nc=1,ncoh
             if(vegsurf(nc)) then
               vegarea_tem(j)=vegarea_tem(j)+ Incohortarea(nc,j)
             endif
!            if(.not.vegsurf(nc)) then
!              nvegfrac(j)=nvegfrac(j)+vegfrac(nc,j)
!            endif
           enddo
          enddo
!          print *, ' '
          first = .false.
        endif
           print *,'From vegadju'
!
        do j=1,jm0
         do nc=1,ncoh
          if(.not.vegsurf(nc)) then
          if (Incohortarea(nc,j).ne.cohortarea(nc,j)) then
             print *,'From checkvegfraction'
             print '(2i3,2i20)',j,nc,Incohortarea(nc,j),cohortarea(nc,j)
             stop
          endif
          endif
         enddo
        enddo
 
        ivegtot_tem=0
        ivegtot_igsm=0
        do j=1,jm0
         vegarea_igsm(j)= 0
         do nc=1,ncoh
           if(vegsurf(nc)) then
             vegarea_igsm(j)=vegarea_igsm(j)+ cohortarea(nc,j)
           endif
         enddo
!
!        Adjustment
         if(vegarea_igsm(j).ne.vegarea_tem(j)) then
          idveg=vegarea_igsm(j)-vegarea_tem(j)
!           print *,'Adjustent'
!           print *,j,idveg
          ivegmax=0
          do nc=1,ncoh
           if(vegsurf(nc)) then
             if(cohortarea(nc,j).gt.ivegmax) then
               ivegmax=cohortarea(nc,j)
               ncmax=nc
             endif
           endif
          enddo
          if(idveg.gt.ivegmax) then
           print *,'DVEG gt VEGMAX'
           print *,idveg,ivegmax
           stop
          else
!           print *,ncmax,coharea(ncmax,j),vegfrac(ncmax,j)
            cohortarea(ncmax,j)=cohortarea(ncmax,j)-idveg
          endif
         endif
!        Adjustment
!
         vegarea_igsm(j)= 0
         do nc=1,ncoh
           if(vegsurf(nc)) then
             vegarea_igsm(j)=vegarea_igsm(j)+ cohortarea(nc,j)
           endif
         enddo
         ivegtot_tem=ivegtot_tem+vegarea_tem(j)
         ivegtot_igsm=ivegtot_tem+vegarea_igsm(j)
         print '(F6.0,2i10)',lat(j)*180./3.14,
     &         vegarea_tem(j),vegarea_igsm(j)
          enddo
         print '(A15,2i10)','GLOBAL VEG AREA',
     &         ivegtot_tem,ivegtot_igsm



!
#endif
	return
	end
