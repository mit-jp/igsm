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
        integer coharea(ncoh,jm0),vegarea(jm0),vegarea_tem(jm0)
     &    ,landarea_igsm(jm0),vegarea_igsm(jm0)
        real nvegfrac(jm0),vegfrac0(ncoh,jm0)
        logical vegsurf
        vegsurf(n)=((n.gt.1.and.n.le.17).or.n.eq.34)
!       vegsurf(n)=((n.gt.0.and.n.le.17).or.n.eq.34)
!
        data first /.true./
        if (first) then
!
            vegfrac0 = vegfrac
!
           print *, ' LAND AREAS'
          do j=1,jm0
!          landarea(j)=0
           vegarea_tem(j)=0
           nvegfrac(j)=0.0
           vegarea_igsm(j)=0
           do nc=1,ncoh
             landarea_igsm(j)=nint(cellarea(j)*landfrac(j))
!            landarea(j)=landarea(j)+ cohortarea(nc,j)
             if(vegsurf(nc)) then
               vegarea_igsm(j)=vegarea_igsm(j)+
     &          nint(vegfrac(nc,j)*landarea_igsm(j))
               vegarea_tem(j)=vegarea_tem(j)+ cohortarea(nc,j)
             endif
!            if(nc.gt.17.and.nc.ne.34) then
             if(.not.vegsurf(nc)) then
               nvegfrac(j)=nvegfrac(j)+vegfrac(nc,j)
             endif
           enddo
           print '(F5.0,f15.3,i15,f10.5)'
     &     ,lat(j)*180./3.14,
     &     cellarea(j)*landfrac(j),landarea_igsm(j)
     &      ,nvegfrac(j)
          enddo
           print *, ' '
          first = .false.
        endif
!
!       do j=1,jm0
!        do nc=1,ncoh
!         if(.not.vegsurf(nc)) then
!            vegfrac(nc,j)=vegfrac0(nc,j)
!         endif
!        enddo
!       enddo
!
        svegtotg=0.0
        vegtotg=0.0
        ivegtot1g=0
        ivegtot2g=0
        ivegtot3g=0
        do j=1,jm0
         vegarea_igsm(j)= 0
         vegarea(j)=0
         do nc=1,ncoh
           coharea(nc,j)= 0
         enddo
         vegtot=0.0
         svegtot=0.0
         do nc=1,ncoh
           coharea(nc,j)=coharea(nc,j) + 
     &       nint(vegfrac(nc,j)*landarea_igsm(j))
           if(vegsurf(nc)) then
             vegarea(j)=vegarea(j)+ coharea(nc,j)
           endif
!          svegtot=svegtot+vegfrac(nc,j)*
!    &       cellarea(j)*landfrac(j)
           if(vegsurf(nc)) then
            vegtot=vegtot+vegfrac(nc,j)*
     &        cellarea(j)*landfrac(j)
           endif
         enddo
!
!        Adjustment
         if(vegarea(j).ne.vegarea_tem(j)) then
          idveg=vegarea(j)-vegarea_tem(j)
!           print *,'Adjustent'
!           print *,j,idveg
          ivegmax=0
          do nc=1,ncoh
           if(vegsurf(nc)) then
             if(coharea(nc,j).gt.ivegmax) then
               ivegmax=coharea(nc,j)
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
            coharea(ncmax,j)=coharea(ncmax,j)-idveg
            vegfrac(ncmax,j)=float(coharea(ncmax,j))/landarea_igsm(j)
!           print *,vegfrac(ncmax,j),coharea(ncmax,j), 
!    &           nint(vegfrac(ncmax,j)*landarea_igsm(j))
            chnew=nint(vegfrac(ncmax,j)*landarea_igsm(j))
!           print *,'Adjustent'
           if(coharea(ncmax,j).ne.chnew) then
             print *,'Adjustent did not work'
           endif
          endif
         endif
!        Adjustment
!
         do nc=1,ncoh
           if(vegsurf(nc)) then
             vegarea_igsm(j)=vegarea_igsm(j)+ coharea(nc,j)
           endif
         enddo
         vegtotg=vegtotg+vegtot
         ivegtot1g=ivegtot1g+vegarea_tem(j)
         ivegtot2g=ivegtot2g+vegarea(j)
         ivegtot3g=ivegtot3g+vegarea_igsm(j)
         print '(F6.0,f15.3,3i10)',lat(j)*180./3.14,
     &         vegtot,vegarea_tem(j),vegarea(j),vegarea_igsm(j)
          enddo
         print '(A6,f15.3,3i10)','GLOBAL',
     &         vegtotg,ivegtot1g,ivegtot2g,ivegtot3g



!
          do j=1,jm0
!           print *,'vegfrac',j,vegfrac(1,j),vegfrac(31,j)
           do nc=1,ncoh
            if(.not.vegsurf(nc)) then
              if(vegfrac(nc,j).ne.vegfrac0(nc,j)) then
               print *,j,nc,vegfrac(nc,j),vegfrac0(nc,j)
              endif
            endif
!           if(cohortarea(nc,j)-coharea(nc,j).ne.0) then
!             print ('2i5,2i15'),j,nc,cohortarea(nc,j),coharea(nc,j)
!           endif
           enddo
          enddo
!       stop
#endif
	return
	end
