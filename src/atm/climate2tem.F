#include "ctrparam.h"
! 	This subroutine is called one a month by driver
!	to normalize monthly climate data for TEM
	subroutine climate2tem(month,ndays)
#include "BD2G04.COM"
#include "TEM.h"        
#include "CLM.h"        
!       logical first,fixedLC,wr88,fixedLC4TEM
!       fixedLC4TEM move to TEM.h
        logical first,fixedLC,wr88
        logical tem44c
!
        real wetlands(35,jm0),wetlat(jm0)
!
        COMMON/CO2TRND/ALFFOR,CO2TR,YEARGT,CO2IN,INYRAD
        data first /.true./
        wr88=.true.
        wr88=.false.
        tem44c=.false.
!       fixedLC4TEM=.true.
!       moved to input.F
        if (first) then
         if(.not.dyn_pft.or.rampyr_dynpft.ne.0) then
          fixedLC=.true.
         else
          fixedLC=.false.
         endif
!TEST
       if(tem44c) then
       print *,' use data from TEMDATA for TEM44d'
       flsotem='/home/sokolov/IGSM2/TEMDATA/igsmsoilchar4x360b1.mit'
       fllatem='/home/sokolov/IGSM2/TEMDATA/igsmslayerchar4x360b1.mit'
       open(394,file=flsotem,status='old')
       open(395,file=fllatem,status='old')
       do j=1,jm0
        read(394,*),lontem(j),lattem(j),xx,areatem(j),sandtem(j),
     &  silttem(j),claytem(j),wsoiltem(j),phtem(j),yy,zz
        do k=1,mnlay1
        read(395,*),lontem(j),lattem(j),xx,nk,zthicktem(k,j),
     &  portem(k,j),ksattem(k,j),yy,zz
      enddo
      enddo
      close(394)
      close(395)
      endif
!TEST
#if ( defined CLM35 )
        print *,'climate4tem for CLM35'
        if(fixedLC) then
          print *,'Run with fixed land cover'
        else 
          print *,'Run with changing land cover'
        endif
        if(fixedLC4TEM) then
           print *,'Run with fixed land cover for TEM'
        else 
          print *,'Run with changing land cover for TEM'
        endif
          open (661,
     &      file=wetfrac_data,status='old',form='unformatted')
          read (661),wetlat,wetlands
          close(661)
!
          cellarea = DXYP*1.e-6
          landfrac =  FDATA(1,:,2)
          do j=1,jm0
             landarea(j)=nint(nint(cellarea(j))*landfrac(j))
          enddo
!         print *,'CELLAREA'
!         print *,CELLAREA
!         print *,'LANDFRAC'
!         print *,LANDFRAC
!         print *,'LANDAREA FOR TEM'
          if(wr88)then
           write(88, *),'LANDAREA FOR TEM'
           do j=1,jm0
!           print '(F5.0,f13.3,f20.16,i13)'
            write (88, '(F5.0,f13.3,f20.16,i13)')
     &      ,lat(j)*180./3.14,cellarea(j),landfrac(j)
     &      ,nint(cellarea(j)*landfrac(j))
           enddo
          endif
!TEST
       if(tem44c) then
! use data from TEMDATA
          do j=1,jm0
            pctsand4tem(j)=sandtem(j)
            pctsilt4te(j)=silttem(j)
            pctclay4tem(j)=claytem(j)
            soilph4tem(j)=phtem(j)
           do k=1,mnlay1
            zthick4tem(k,j)=zthicktem(k,j)
            por4tem(k,j)=portem(k,j)
            ksat4tem(k,j)=ksattem(k,j)
           enddo
          enddo
       endif
!TEST
!         print *,'TEXTURE'
!         do j=1,jm0
!          print '(F5.0,f13.3,3f6.2)'
!    &     ,lat(j)*180./3.14,cellarea(j)
!    &     , pctsand4tem(j),pctsilt4te(j),pctclay4tem(j)
!         enddo
!         print *,'SOILAYER'
!         do j=1,jm0
!          do k=1,mnlay1
!          do k=1,1
!          print '(F5.0,i3,3f15.5)'
!    &     ,lat(j)*180./3.14,k,zthick4tem(k,j)
!    &     , por4tem(k,j),ksat4tem(k,j)
!          enddo
!         enddo
#if ( defined DATA4TEM )
          open (931,
     &      file='lndfrac4tem',status='new',form='unformatted')
          write(931),cellarea,landarea
     & , pctclay4tem, pctsand4tem, pctsilt4te
     & ,soilph4tem, ksat4tem, por4tem, zthick4tem
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
c       print *,'O3 for TEM J=23'
c       print *,(o34tem(n,23),n=1,no3)
#if ( defined CLM35 )

        if(fixedLC4TEM) then
        do j=1,jm0
        do nc=1,35
         cohortarea(nc,j)=Incohortarea(nc,j)
        enddo
        enddo

        else
!  For crop and pasture
!   in CLM3.5 crop - 16, pasture - 34
!   in TEM44c crop1 - 16, crop2 -17, pasture - 34
!  Defining grop1 and crop2
       if (month.eq.1) then
        croptotal =vegfrac(16,:)
        if(fixedLC)then
         if(.not.dyn_pft) then
          iycrop=YEARGT
         else
          iycrop=rampyr_dynpft
         endif
        else
          iycrop=JYEAR
        endif
        if(iycrop.le.1950) then
          crop1=1.0
        elseif(iycrop.lt.1990) then
          crop1=1.0-1./40.*(iycrop-1950)
        else
          crop1=0.0
        endif
!      endif
!
          vegfrac(16,:) = crop1* croptotal
          vegfrac(17,:) = (1.-crop1)*croptotal
        print *,' CROP1 for year ',iycrop,crop1
        
!   in CLM3.5 total wetland area - 18
!   in TEM44c wetlands - 18 - 30
!      
        do j=1,jm0
         totwet=vegfrac(18,j)
        do nc=18,30
         vegfrac(nc,j)=totwet*wetlands(nc,j)
        enddo
        enddo

        do j=1,jm0
        do nc=1,35
         cohortarea(nc,j)=
     &   nint(landarea(j)*vegfrac(nc,j))
        enddo
        enddo


      if(wr88)then
       write(88, *),' From climate2tem '
       do j=1,jm0
!      print '(i3, f12.2,f20.16,i20)'
       write( 88, '(i3,2f12.2,i20,i20)')
     &    ,j,lat(j)*180./3.14,cellarea(j),landarea(j)
       do nc=1,ncoh
!      print '(i3, f20.16,2i20)'
       write( 88, '(i3, f30.26,3i20)')
     & ,nc,vegfrac(nc,j),Incohortarea(nc,j),
     &   cohortarea(nc,j)
     &  ,Incohortarea(nc,j)-cohortarea(nc,j)
       enddo
       enddo
      endif

      if(.not.fixedLC)then
       call vegfracadjustment

       if(wr88)then
        write(88, *),' After vegfracadjustment '
        do j=1,jm0
        write( 88, '(i3,2f12.2,f20.16,i20)')
     &    ,j,lat(j)*180./3.14,cellarea(j),landfrac(j)
        do nc=1,ncoh
        write( 88, '(i3, f30.26,3i20)')
     &  ,nc,vegfrac(nc,j),Incohortarea(nc,j),
     &    cohortarea(nc,j)
     &   ,Incohortarea(nc,j)- cohortarea(nc,j)
        enddo
        enddo
       endif
      else
!       call checkvegfraction
        call vegfracadjustment 
      endif

!
        endif  ! fixedLC4TEM
        endif  ! nmoth.eq.1
#endif
#if ( defined DATA4TEM )
        write(937),month,ndays,co24tem,temp4tem,pre4tem,sws4tem,
     &		o34tem,dtem4tem
#if ( defined CLM35 )
        write(938),month,ndays,strmdur,qstrm,pet,aet,sh2o1m,sh2o2m,
     &  swe,sfr,drn,vegfrac
#if ( defined DATA4NEM )
        write(939),month,ndays,daytsoil,daysho,rsh2o
#endif
        jpr=-10
        ncpr=12
!       do j=1,jm0
        j=jpr
        if(vegfrac(ncpr,jpr).gt.0.0)then
        print *,'For lat=',lat(jpr)*180./3.14,' cohort=',ncpr
           coharea=vegfrac(ncpr,j)*
     &     cellarea(j)*landfrac(j)
!       print ('(1p5e12.3)'),coharea,co24tem(jpr),temp4tem(jpr),
        print *,coharea,co24tem(jpr),temp4tem(jpr),
     &    sws4tem(jpr),pre4tem(jpr)
        print *,pet(ncpr,jpr),aet(ncpr,jpr),
!       print ('(1p4e12.3)'),pet(ncpr,jpr),aet(ncpr,jpr),
     &     sh2o1m(ncpr,jpr),
     &     sh2o2m(ncpr,jpr)
        print *,swe(ncpr,jpr),sfr(ncpr,jpr),drn(ncpr,jpr)
!       print ('(1p4e12.3)'),swe(ncpr,jpr),sfr(ncpr,jpr),drn(ncpr,jpr)
        endif
!       enddo
#endif
#endif
!       stop
	return
	end
