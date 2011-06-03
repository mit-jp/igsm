
#include "ctrparam.h"

      subroutine stbtrnd(YEAR,mndriver)
      parameter (nyd=2000,nyd1=nyd+1)
C
C
C     For IPCC stabilisation simulations
C     CO2 concentrations  from file co2_data
C     Concentrations of other GHGs are
c     those for 1765.
C     DATA are annual means, but they are
C       used as a data for the middle of the coresponding year
C
C
!     parameter(jm0=46)
#include "BD2G04.COM"
#include "TEM.h"
#include "RADCOM.COM"
      common /Garyflux/pC_atm(jm0),wind_amp,fluxco2(jm0)
      real XGASS(NGAS)
      dimension xemi(nyd),ighgyr(nyd),xxemi(nyd1)
      real emisdata(5)
      dimension NDAYMN(12)
      data NDAYMN /31,28,31,30,31,30,31,31,30,31,30,31/

      common/boxmod/YEAR0
      common /bmtrdata/co2_data
      common/ BACKGRGHG/GHGBGR(5)
      character * 64 co2_data
      character * 10 name,scen(5)
      logical first
      data first /.true./
      data ncall /1/
      if(first)then
         nsc=KTREN-KTREN/10*10
         open (unit=561,file=co2_data,
     *     status='OLD',form='formatted')
         read (561,*)name
         read (561,*)name,scen
         print *,' Simulation with emission scenario ',scen(nsc)
       do 50 i=1,nyd
         read (561,*,end=500),ighgyr(i),xemi(i)
   50 continue
  500 continue
      close (561)
      ny=i-2
      ny1=ny+1
      ny2=ny+2
      YEAR0=ighgyr(1)-1
        print *,'From FORSET EMISSIPCC'
        print *,' year=',year
        print *,'YEAR0=',YEAR0
  501    format(i4,2f8.2)
         xxemi(1)=2.*xemi(1)-xemi(2)
         do 10 i=1,ny1
           xxemi(i+1)=xemi(i)
   10    continue
         first=.false.
      endif
       do j=1,jm
         OCUPT=OCUPT+fluxco2(j)
         TEMUPTANN=TEMUPTANN+1.e-3*temco2(j)/NDAYMN(mndriver)
       enddo

        call ghgint(xx,xxemi,ny2,YEAR)
        yy0=FULGAS(2)*GHGBGR(1)
        if(ncall.eq.2)then
           print *,'IPPCCEMIS '
           print *,' year=',year
           print *,'YEAR0=',YEAR0
           print *,'XGASS'
           print *,XGASS
        endif
        if(ncall.lt.3)ncall=ncall+1
      return
      end

#include "ctrparam.h"
C
      SUBROUTINE emissipcc_mn(nepmn,ocumn,xco2ann)
C     |==========================================================|
C     | This is moved out of the main driver code - uses BD2G04  |
C     | Called only at end-of-month period                       |
C     \==========================================================/
!     IMPLICIT NONE

#include "BD2G04.COM"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
      REAL    nepmn
      REAL    ocumn
      REAL    xco2ann

C     LOCAL VARIABLES
!     REAL*4 tnow

      tnow= jyear + (jday-.5)/365.
      print *,'TNOW from emissipcc_mn=',TNOW     
      print *,'CO2=',co2
      print *,nepmn,ocumn,xco2ann

      CALL emissipcc(tnow,nepmn,ocumn,co2,xco2ann)

      RETURN
      END
C

#include "ctrparam.h"
      SUBROUTINE emissipcc_yr
C     |==========================================================|
C     | This is moved out of the main driver code - uses BD2G04  |
C     | Called only at end-of-year period.                       |
C     \==========================================================/
!     IMPLICIT NONE

#include "BD2G04.COM"

      COMMON/ BACKGRGHG/ghgbgr
       REAL   ghgbgr(5)

      REWIND 861
      WRITE (861,*)co2*ghgbgr(1)

      RETURN
      END

         
      subroutine emissipcc(YEAR,tcu,ocu,xco2,xco2ann)
      parameter (nyd=2000,nyd1=nyd+1)
C        
C      
C     For IPCC stabilisation simulations with CO2 emissions
C     CO2 emissions  from file co2_data
C     Concentrations of other GHGs are
c     those for 1765.
C     DATA are annual means, but they are
C       used as a data for the middle of the coresponding year
C     
C     
#include "BD2G04.COM"
!include "TEM.h"
!include "RADCOM.COM"
c     common /Garyflux/pC_atm(jm0),wind_amp,fluxco2(jm0)
!     parameter(jm0=46)
      common /ATCO2/atm_co2(jm0)
      dimension xemi(nyd),eghgyr(nyd),xxemi(nyd1),iyeare(nyd1)
      real emidata(4),co2stocker(9),xxco2(nyd)
      dimension NDAYMN(12)
      data NDAYMN /31,28,31,30,31,30,31,31,30,31,30,31/
  
      common/boxmod/YEAR0
      common /bmtrdata/co2_data
      common /ipccdata/init_co2
      common/ BACKGRGHG/GHGBGR(5)
      character * 120 co2_data
      character * 30 name
      logical first
      data first /.true./
      data ncall /1/
      if(first)then
         id2=index(co2_data," ")
         open (unit=561,file=co2_data,
     *     status='OLD',form='formatted')
         read (561,*)name
         read (561,*)name
!        print *,' Simulation with emission scenario ',co2_data(1:id2-1)
         print *,' Simulation with emission scenario ',name
       do 50 i=1,nyd
         read (561,*,end=500),eghgyr(i),xemi(i)
         print *,eghgyr(i),xemi(i)
   50 continue
  500 continue
      close (561)
      ny=i-2
      ny1=ny+1
      ny2=ny+2
      YEAR0=eghgyr(1)-1
        print *,'From FORSET EMISSIPCC'
        print *,' year=',year
        print *,'YEAR0=',YEAR0
        print *,'ny=',ny
  501    format(i4,2f8.2)
         xxemi(1)=xemi(1)
         iyeare(1)=YEAR0
         do 10 i=1,ny1
           iyeare(i+1)=eghgyr(i)
           xxemi(i+1)=xemi(i)
   10    continue
         print *,(iyeare(i),i=1,ny2)
         print *,(xxemi(i),i=1,ny2)
         dtemi=1./12.
         first=.false.
      endif

c        call ghgint(xx,xxemi,ny2,YEAR)
         call ghgint1(xx,xxemi,iyeare,ny2,YEAR)
        yy0=xco2*GHGBGR(1)
        xco2ann=xco2ann+yy0
c       yy0 and yy CO2 in PPM, xx in GtC/year, ocu and tcu in GtC/month
c       dtemi=1./12.
c       print *,'ocu=',ocu,' tcu=',tcu
c       ocu=1.25*ocu
c       tcu=1.8*tcu
        print *,'ocu=',ocu,' tcu=',tcu
        print *,'emi=',xx*dtemi
        yy=yy0+(xx*dtemi-ocu-tcu)/2.1
        xco2=yy/GHGBGR(1)
c       if(ncall.eq.2)then
           print *,'EMISSIPCC '
           print *,' year=',year
           print *,'YEAR0=',YEAR0
           print *,yy0,xx
           print *,ocu,tcu
           print *,yy,xco2
         do 30 j=1,jm0
          atm_co2(j)=yy
   30    continue
c          stop
c       endif
        if(ncall.lt.3)ncall=ncall+1
      return
      end


      subroutine ghgint1(xx,xxf,iyeare,ny,year)
      common/boxmod/YEAR0
      dimension xxf(ny),iyeare(ny)
!     print *,ny,iyeare(1),xxf(1)
      do i=1,ny
       year1=iyeare(i)+0.5
       year2=iyeare(i+1)+0.5
c      if(year.gt.iyeare(i).and.year.le.iyeare(i+1))go to 100
       if(year.ge.year1.and.year.le.year2)go to 100
      enddo
       print *,' Wrong year year=',year
       print *,'YEAR=',YEAR,year1,year2
       stop
  100  continue
c        print *,'YEAR=',YEAR,year1,year2
c       x=(iyeare(i+1)-year)/(iyeare(i+1)-iyeare(i))
        x=(year2-year)/(year2-year1)
         xx=xxf(i+1)*(1.-x)+xxf(i)*x
c       print *,i,x,xxf(i),xxf(i+1),xx
      return
      end
