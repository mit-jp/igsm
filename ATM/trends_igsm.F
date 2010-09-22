C $Header: /home/igsm-repository/ATM/trends_igsm.F,v 1.1 2010/08/30 15:36:35 jscott Exp $
C $Name:  $

#include "ctrparam.h"

      SUBROUTINE GTREND(XGAS,YEART,NGAS)                                9155.9
C                                                                       9156.
C-----------------------------------------------------------------------9156.1
C     T-GAS SCENARIO  G:                                                9156.2
C                                                                       9156.3
C                                                                       9156.4
C                   CO2:                                                9156.5
C                          0.01/YR INCREASE  form 1977                  9156.6
C                                                                       9156.7
C                                                                       9156.8
C                           N2O, CH4 and CFC for 1977
C                                                                       9156.9
C                   N2O:   WEISS FORMULA THRU 1978, 0.035/YR INCREASE   9157.
C                          IN 1978 INCREMENT THRU 1990, 0.025/YR        9157.1
C                          THRU 2000, 0.015/YR THRU 2010, 0.005/YR      9157.2
C                          INCREASE IN INCREMENT BEYOND 2010            9157.3
C                                                                       9157.4
C                   CH4:   OBSERVED RATES THRU 1980,                    9157.5
C                          0.015/YR INCREASE THRU 1990, 0.010/YR        9157.6
C                          INCREASE THRU 2000, 0.005/YR BEYOND 2000     9157.7
C                                                                       9157.8
C                   CFC:   OBSERVED CONCENTRATION THRU 1984,            9157.9
C                          0.03 EMISSION RATE INCREASE TO 1990,         9158.
C                          0.02 EMISSION RATE INCREASE TO 2000,         9158.1
C                          0.01 EMISSION RATE INCREASE TO 2010,         9158.2
C                          CONSTANT EMISSION RATE BEYOND 2010           9158.3
C                                                                       9158.4
C-----------------------------------------------------------------------9158.5
C                                                                       9158.6
      DIMENSION XGAS(5)                                                 9158.7
      DIMENSION XXF11(45),XXF12(45),XXCO2(27)                           9158.8
      COMMON/CO2TRND/ALFFOR,CO2TR,YEARGT,CO2IN,INYRAD
      common/ BACKGRGHG/GHGBGR(5)
      save
      data ifirst /0/
      DATA XXF11/                                                       9158.9
     10.00001,0.00001,0.00001,0.00001,0.00001,0.00004,0.00010,0.00023,  9159.
     20.00041,0.00066,0.00100,0.00148,0.00215,0.00297,0.00397,0.00523,  9159.1
     30.00662,0.00791,0.00921,0.01093,0.01316,0.01596,0.01939,0.02345,  9159.2
     40.02807,0.03321,0.03904,0.04569,0.05343,0.06224,0.07186,0.08267,  9159.3
     50.09505,0.10867,0.12174,0.13422,0.14665,0.15820,0.16871,0.17841,  9159.4
     60.18919,0.19983,0.21032,0.22068,0.23090/                          9159.5
      DATA XXF12/                                                       9159.6
     10.00001,0.00001,0.00001,0.00040,0.00112,0.00198,0.00298,0.00411,  9159.7
     20.00539,0.00679,0.00836,0.00998,0.01181,0.01387,0.01619,0.01888,  9159.8
     30.02195,0.02513,0.02868,0.03292,0.03766,0.04312,0.04954,0.05702,  9159.9
     40.06546,0.07485,0.08545,0.09756,0.11130,0.12634,0.14250,0.16011,  9160.
     50.17962,0.20082,0.22144,0.24105,0.25951,0.27630,0.29247,0.30831,  9160.1
     60.32585,0.34328,0.36059,0.37778,0.39486/                          9160.2
      DATA XXCO2/                                                       9160.3
     1 315.00, 315.69, 316.66, 317.27, 318.24, 318.93, 319.33, 319.95,  9160.4
     2 320.80, 321.44, 322.25, 323.44, 324.82, 326.01, 326.97, 328.86,  9160.5
     3 330.35, 330.71, 331.69, 332.83, 334.58, 335.95, 337.65, 339.04,  9160.6
     4 340.01, 341.56, 343.81/                                          9160.7
C                                                                       9160.8
      if (ifirst.eq.0)then
       if(ALFFOR.gt.0.0)then
         CO2IN=GHGBGR(1)
       endif
       print *,' gtrend.f',100*ALFFOR,'% per year from',CO2IN
       print *,'N2O, CH4 and CFCs '
       print '(4e10.5)',(GHGBGR(ii),ii=2,5)
c      ifirst=1
      endif
      YEAR=YEART
      DT=YEAR-INYRAD
      XX0=CO2IN
      XX=XX0*(1.+ALFFOR)**DT
      XGAS(1)=XX                                                        9182.3
      CO2TR=XGAS(1)
       do ii=2,5
         XGAS(ii)=GHGBGR(ii)
        enddo
      if (ifirst.eq.0)then
       print *,'From GTREND'
       print '(5e12.4)',(XGAS(ii),ii=1,5)
       ifirst=1
      endif
      RETURN                                                            9182.5
      END                                                               9182.6
!	----------------------------------------------------------

      subroutine bmtrnd(XGASS,YEAR,NGAS)
      parameter (nyd=2000,nyd1=nyd+1)
C
C
C     For simulations with observed forcing using equivalent
C     CO2 concentrations  from file co2_data
C     Concentrations of other GHGs are
c     those for 1958.
C     DATA are annual means, but they are
C       used as a data for the middle of the coresponding year
C
C
      real XGASS(NGAS)
      dimension XG1958(5), xf11(nyd),xf12(nyd),xco2(nyd),xn2o(nyd),
     *xch4(nyd),ighgyr(nyd)
      real xxf11(nyd1),xxf12(nyd1),xxco2(nyd1),xxn2o(nyd1),xxch4(nyd1)
      common/boxmod/YEAR0
      common /bmtrdata/co2_data
      common/ BACKGRGHG/GHGBGR(5)
      character * 64 co2_data
      logical first
      data first /.true./
      data ncall /1/
      if(first)then
         print *,' BMTREND for observed forcing'
         open (unit=561,file=co2_data,
     *     status='OLD',form='formatted')
       do 50 i=1,nyd
c        read (561,501,end=500),ighgyr(i),tmp,xco2(i)
         read (561,*,end=500),ighgyr(i),xco2(i)
   50 continue
  500 continue
      close (561)
      do ii=1,5
       XG1958(ii)=GHGBGR(ii)
      enddo
      ny=i-2
      ny1=ny+1
      ny2=ny+2
      YEAR0=ighgyr(1)-1
        print *,'From FORSET'
        print *,' year=',year
        print *,'YEAR0=',YEAR0
  501    format(i4,2f8.2)
c        xxco2(1)=XG1958(1)
         xxco2(1)=2.*xco2(1)-xco2(2)
         xxn2o(1)=XG1958(2)
         xxch4(1)=XG1958(3)
         xxf11(1)=XG1958(4)*1000.
         xxf12(1)=XG1958(5)*1000.
         do 10 i=1,ny1
           xxco2(i+1)=xco2(i)
   10    continue
         first=.false.
      endif
        if(ncall.eq.2)then
           print *,'BMTRND from FORGET'
           print *,' year=',year
           print *,'YEAR0=',YEAR0
        endif
        call ghgint(xx,xxco2,ny2,YEAR)
        XGASS(1)=xx
        XGASS(2)=XG1958(2)                                                     
        XGASS(3)=XG1958(3)                                                     
        XGASS(4)=XG1958(4)*1000.                                                
        XGASS(5)=XG1958(5)*1000.                                            
        if(ncall.lt.3)ncall=ncall+1
      return
      end

      subroutine ghgint(xx,xxf,ny,year)
      common/boxmod/YEAR0
      real xxf(ny)
      i=year+0.5
      dl=year+.5-i
      i=i-YEAR0
      if(i.lt.1.or.i.ge.ny)then
         print *,' error in ghgint i=',i
         stop 25
      endif
c     print *,'YEAR=',YEAR
      xx=xxf(i)*(1.-dl)+xxf(i+1)*dl
c     print *,i,dl,xxf(i),xxf(i+1),xx
      return
      end
      subroutine bmtrndmg(XGASS,YEAR,NGAS)
      parameter (nyd=2000,nyd1=nyd+1)
c
c     Multi-gas version of box-model trend routine
c     - added 22 Aug 2002 - ceforest
C     - greenhouse gas concentrations taken from GISS dataset on web
c     - CFC-11 includes other long-lived gases
C
C     For simulations with observed forcing using equivalent
C     CO2 concentrations  from file co2_data
C     Concentrations of other GHGs are
c     those for 1958.
C     DATA are annual means, but they are
C       used as a data for the middle of the coresponding year
C
C
      real XGASS(NGAS)
      dimension xf11(nyd),xf12(nyd),xco2(nyd),xn2o(nyd),
     *xch4(nyd),ighgyr(nyd)
      real xxf11(nyd1),xxf12(nyd1),xxco2(nyd1),xxn2o(nyd1),xxch4(nyd1)
      common/boxmod/YEAR0
      common /bmtrdata/co2_data
      common/ BACKGRGHG/GHGBGR(5)
      character * 64 co2_data
      logical first
      data first /.true./
      data ncall /1/
      if(first)then
         print *,' BMTRENDMG for observed forcing'
         open (unit=561,file=co2_data,
     *     status='OLD',form='formatted')
       do 50 i=1,nyd
c        read (561,501,end=500),ighgyr(i),tmp,xco2(i)
c     f11 includes other long-lived gases
         read (561,*,end=500),ighgyr(i),xco2(i),xn2o(i),
     *          xch4(i), xf11(i),xf12(i)
   50 continue
  500 continue
      close (561)
      do ii=1,3
       XGASS(ii)=GHGBGR(ii)
      enddo
      do ii=4,5
       XGASS(ii)=GHGBGR(ii)*1000.
      enddo
!     CFC11 and CFC12 are converted in ppb 
      ny=i-2
      ny1=ny+1
      ny2=ny+2
      YEAR0=ighgyr(1)-1
        print *,'From FORSET'
        print *,' TREF=',YEAR
        print *,' XREF'
        print *,XGASS
        print *,' year=',year
        print *,' ghgdata for ',ny,' years from ', ighgyr(1), ' to '
     &  , ighgyr(ny)
  501    format(i4,2f8.2)
         xxco2(1)=2.*xco2(1)-xco2(2)
         xxn2o(1)=2.*xn2o(1)-xn2o(2)
         xxch4(1)=2.*xch4(1)-xch4(2)
         xxf11(1)=(2.*xf11(1)-xf11(2))*1000.
         xxf12(1)=(2.*xf12(1)-xf12(2))*1000.
         do 10 i=1,ny1
           xxco2(i+1)=xco2(i)
           xxn2o(i+1)=xn2o(i)
           xxch4(i+1)=xch4(i)
           xxf11(i+1)=xf11(i)*1000.
           xxf12(i+1)=xf12(i)*1000.
   10    continue
         first=.false.
         return
      endif
        call ghgint(xx,xxco2,ny2,YEAR)
        XGASS(1)=xx
        call ghgint(xx,xxn2o,ny2,YEAR)
        XGASS(2)=xx
        call ghgint(xx,xxch4,ny2,YEAR)
        XGASS(3)=xx
        call ghgint(xx,xxf11,ny2,YEAR)
        XGASS(4)=xx
        call ghgint(xx,xxf12,ny2,YEAR)
        XGASS(5)=xx
        if(ncall.eq.2)then
           print *,'BMTRND from FORGET'
           print *,' year=',year
           print *,'YEAR0=',YEAR0
        print *,' XGASS'
        print *,XGASS
        endif
        if(ncall.lt.3)ncall=ncall+1
      return
      end
