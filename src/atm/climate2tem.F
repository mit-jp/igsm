#include "ctrparam.h"
! 	This subroutine is called one a month by driver
!	to normalize monthly climate data for TEM
	subroutine climate2tem(month,ndays)
#include "BD2G04.COM"
#include "TEM.h"        
#include "CLM.h"        
        logical first,fixedLC
!
        real wetlands(35,jm0),wetlat(jm0)
!
        data first /.true./
        if (first) then
        fixedLC=.not.dyn_pft
!       fixedLC=.false.
!       rampyr_dynpft=0
!
#if ( defined CLM35 )
        print *,'climate4tem for CLM35'
        if(fixedLC) then
          print *,'Run with fixed land cover'
        endif
          open (661,
     &      file=wetfrac_data,status='old',form='unformatted')
          read (661),wetlat,wetlands
          close(661)
!
          cellarea = nint(DXYP*1.e-6)
          landfrac =  FDATA(1,:,2)
!         print *,'CELLAREA'
!         print *,CELLAREA
!         print *,'LANDFRAC'
!         print *,LANDFRAC
          print *,'LANDAREA FOR TEM'
          do j=1,jm0
           print '(F5.0,f13.3,f7.3,f15.3)'
     &     ,lat(j)*180./3.14,cellarea(j),landfrac(j)
     &     ,cellarea(j)*landfrac(j)
          enddo
#if ( defined DATA4TEM )
          open (931,
     &      file='lndfrac4tem',status='new',form='unformatted')
          write(931),cellarea,landfrac
          close (931)
#endif
#else
        print *,'climate4tem for CLM21'
#endif
          first = .false.
          return
        endif
	do j=1,jm0
	  temp4tem(j)=temp4tem(j)/ndays
	  sws4tem(j)=sws4tem(j)/ndays
	  pre4tem(j)=pre4tem(j)/ndays
#if ( !defined CPL_CHEM )
	  do n=1,no3
	   o34tem(n,j)=obso3(n,j,month)
	  enddo
#endif
	enddo
        print *,'Month=',month,' CO2(23) for TEM=',co24tem(jm0/2)
#if ( defined CLM35 )
c       print *,'O3 for TEM J=23'
c       print *,(o34tem(n,23),n=1,no3)
!  For crop and pasture
!   in CLM3.5 crop - 16, pasture - 34
!   in TEM44c crop1 - 16, crop2 -17, pasture - 34
!  Defining grop1 and crop2
       if (month.eq.1) then
        croptotal =vegfrac(16,:)
        if(fixedLC)then
!
         if(INYEAR.le.1950) then
           crop1=1.0
         elseif(INYEAR.lt.1990) then
           crop1=1.0-1./40.*(INYEAR-1950)
         else
           crop1=0.0
         endif
         iycrop=INYEAR
!
        else
!
         if(rampyr_dynpft.eq.0) then
!
          if(JYEAR.le.1950) then
            crop1=1.0
          elseif(JYEAR.lt.1990) then
            crop1=1.0-1./40.*(JYEAR-1950)
          else
            crop1=0.0
          endif
          iycrop=JYEAR
!
         else
!
          if(rampyr_dynpft.le.1950) then
            crop1=1.0
          elseif(rampyr_dynpft.lt.1990) then
            crop1=1.0-1./40.*(rampyr_dynpft-1950)
          else
            crop1=0.0
          endif
          iycrop=rampyr_dynpft
!
         endif
!
        endif
!
          vegfrac(16,:) = crop1* croptotal
          vegfrac(17,:) = (1.-crop1)*croptotal
        print *,' CROP1 for year ',iycrop,crop1
!
        
!   in CLM3.5 total wetland area - 18
!   in TEM44c wetlands - 18 - 30
!      
        do j=1,jm0
         totwet=vegfrac(18,j)
        do nc=18,30
         vegfrac(nc,j)=totwet*wetlands(nc,j)
        enddo
        enddo
!
         call vegfracadjustment

        endif
#endif
#if ( defined DATA4TEM )
        write(937),month,ndays,co24tem,temp4tem,pre4tem,sws4tem,
     &		o34tem,dtem4tem
!       write(938),month,ndays,strmdur,qstrm,pet,aet,sh2o1m,sh2o2m,
!    &  swe,sfr,drn,vegfrac
!       write(939),month,ndays,daytsoil,daysho,rsh2o
#endif
	return
	end
