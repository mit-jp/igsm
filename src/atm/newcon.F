
#include "ctrparam.h"

!	==========================================================
!
!	NEWCON.F:  A subroutine for convert predicted concentrations
!			of chemical species to GISS radiation
!			scales.
!
!	----------------------------------------------------------
!
!	Revision History:
!	
!	When	Who		What
!	----	----------	-------	
!	080100	Chien Wang	repack based on CliChem3 & M24x11,
!				  and add cpp.
!
!	==========================================================

      subroutine newcon(co2n,co2o,xn,yn,xo,yo,co20,x0,y0)
!       x - N2O; y - CH4

      FFN3(X,Y)=3.1515*X/(1.0+0.21448*X**0.87) - 0.056*X*Y              9131.5
     +          +0.3234*Y/(1.0+0.052*Y**0.84)
      FFN1(X,Y)=1.556*LOG(1.0+1.098*X**0.77*(1.0+0.032*X)               9113.3
     +                                       /(1.0+0.0014*X*X))         9113.4
     +          +(0.394*Y**0.66+0.16*Y*EXP(-1.6*Y))                     9113.5
     +                         /(1.0+0.169*Y**0.62)                     9113.6
     +          -0.14*LOG(1.0+0.636*(X*Y)**0.75+0.007*Y*(X*Y)**1.52)    9113.7
#if ( !defined IPCC_FORCING )
      FCH4(X,Y)=FFN1(X,Y)
      FN2O(X,Y)=FFN1(X,Y)
#else
!     FCO2(X)=CF*5.35*LOG(X)
!     Instanteneous forcing
!     FCO2(X)=1.136*CF*5.35*LOG(X)
!      Adjustment for 3D CO2 fit
      FCO2(X)=1.1024*CF*5.35*LOG(X)
      GI(X,Y)=0.47*DLOG(1.+2.01e-5*(X*Y)**0.75+5.31e-15*Y*(X*Y)**1.52)
      FCH4(X,Y)=CF*(0.036*SQRT(Y)-GI(X,Y))
      FN2O(X,Y)=CF*(0.12*SQRT(X)-GI(X,Y))
      FCO23(X)=DLOG(1.0 + 1.2*X + 0.005*X**2 + 1.4E-06*X**3)
#endif
        data ifirst /1/
        data inprnt /26/
        if(ifirst.eq.1)then
         CF=1.22/4.3
         sfac=1.
#if ( defined IPCC_FORCING )
         sfac=1000.
#endif
c        print *,x0,y0,xo,yo
         frCH41=FCH4(x0*sfac,y0*sfac)
         frN2O1=FN2O(x0*sfac,y0*sfac)
         frCH43=ffn3(x0,y0)
         frN2O3=ffn3(x0,y0)
#if ( defined IPCC_FORCING )
         frCO21=FCO2(co20)
         frCO23=FCO23(co20)
#endif
         ifirst=0
        endif
        fx1=FN2O(xo*sfac,y0*sfac)-frN2O1
        fy1=FCH4(x0*sfac,yo*sfac)-frCH41
        x=xo
        iter=0
  100   continue
        fx3=ffn3(x,y0)-frN2O3
!APS 11.23.13
        fx3=1.5*fx3
!   N2O forcing is 1.5 larger that given by fit
        if(inprnt.lt.25) then
        print *,' N2O iter=',iter
        print *,' xo=',xo,' fx1=',fx1
        print *,' x=',x,' fx3=',fx3
        endif
        dfdx=fx3/(x-x0)
        df=fx1-fx3
        x=x+df/dfdx
        iter=iter+1
        if(iter.lt.15.and.abs(df).gt.1.e-5) go to 100
        xn=x
        iter=0
        y=yo
  200   continue
        fy3=ffn3(x0,y)-frCH43
!APS 11.23.13
        fy3=1.6*fy3
!   CH4 forcing is 1.6 larger that given by fit
        if(inprnt.lt.25)then
        print *,' iter=',iter
        print *,' yo=',yo,' fy1=',fy1
        print *,' y=',y,' fy3=',fy3
        endif
        dfdy=fy3/(y-y0)
        df=fy1-fy3
        y=y+df/dfdy
        iter=iter+1
        if(iter.lt.15.and.abs(df).gt.1.e-5) go to 200
        yn=y
#if ( defined IPCC_FORCING )
        c=co2o
        iter=0
        fc1=FCO2(c)-frCO21
  300   continue
        fc3=FCO23(c)-frCO23
        dfdx=fc3/(c-co20)
        df=fc1-fc3
        c=c+df/dfdx
        iter=iter+1
        if(iter.lt.15.and.abs(df).gt.1.e-5) go to 300
        co2n=c
#else
         co2n=co2o
#endif
!      print *,'From NEWCON'
!      print *,'N2O ',x0,xo,xn
!      print *,'CH4 ',y0,yo,yn
        inprnt=inprnt+1
      return
      end

!	==========================================================
!
!	NEWCON_ETMINAN.F:  A subroutine for convert predicted concentrations
!			of chemical species to GISS radiation
!			scales.
!
!	----------------------------------------------------------
!
!	Revision History:
!	
!	When	Who		What
!	----	----------	-------	
!      06.02.2020  APS
!      Radiative forcing of carbon dioxide, methane, and nitrous
!      oxide: A significant revision of the methane radiative
!      forcing.   Etminan et al GRL 2016
!
!	==========================================================

      subroutine newcon_ETMINAN(co2n,co2o,xn,yn,xo,yo,co20,x0,y0)
!       x - N2O; y - CH4
!       XN - N2O; XM - CH4

#if ( defined ETMINAN_FORCING )
!   Fits for 3D forcing
      FFN3(X,Y)=3.1515*X/(1.0+0.21448*X**0.87) - 0.056*X*Y              9131.5
     +          +0.3234*Y/(1.0+0.052*Y**0.84)
      FCO23(X)=DLOG(1.0 + 1.2*X + 0.005*X**2 + 1.4E-06*X**3)



      FCO2I(C,C0,XN,XN0)=(a1*(C-C0)**2+b1*(C-C0)+c1*0.5*(XN+XN0)
     &   +5.36)*LOG(C/C0)
!    Adjustment for 3D CO2 fit
     &   *1.1024*CF
#ifdef ADJUSTED_CO2
     &   *4.0/4.9
#endif

      FCH4I(XM,XM0,XN,XN0)=(a3*0.5*(XM+XM0)+b3*0.5*(XN+XN0)
     &   +0.043)*(SQRT(XM)-SQRT(XM0))
     &   *CF
      FN2OI(XM,XM0,XN,XN0,C,C0)=(a2*0.5*(C+C0)+b2*0.5*(XN+XN0)
     &   +c2*0.5*(XM+XM0)+0.117)*(SQRT(XN)-SQRT(XN0))
     &   *CF
      FF11I(X)=CF*0.25*X
      FF12I(X)=CF*0.32*X
        data ifirst /1/
        data inprnt /26/
        if(ifirst.eq.1)then
        print *,'newcon for ETMINAN_FORCING'
      a1=-2.4e-7
      b1=7.2e-4
      c1=-2.1e-4
      a2=-8.e-6
      b2=4.2e-6
      c2=-4.9e-6
      a3=-1.3e-6
      b3=-8.2e-6
         CF=1.22/4.3
         sfac=1000.
         frCH43=ffn3(x0,y0)
         frN2O3=ffn3(x0,y0)
         frCO23=FCO23(co20)
         ifirst=0
        endif
        fx1=FN2OI(yo*sfac,y0*sfac,xo*sfac,x0*sfac,co2o,co20)
        fy1=FCH4I(yo*sfac,y0*sfac,xo*sfac,x0*sfac)
        x=xo
        iter=0
  100   continue
        fx3=ffn3(x,y0)-frN2O3
!APS 11.23.13
        fx3=1.5*fx3
!   N2O forcing is 1.5 larger that given by fit
        if(inprnt.lt.25) then
        print *,' N2O iter=',iter
        print *,' xo=',xo,' fx1=',fx1
        print *,' x=',x,' fx3=',fx3
        endif
        dfdx=fx3/(x-x0)
        df=fx1-fx3
        x=x+df/dfdx
        iter=iter+1
        if(iter.lt.15.and.abs(df).gt.1.e-5) go to 100
        xn=x
        iter=0
        y=yo
  200   continue
        fy3=ffn3(x0,y)-frCH43
!APS 11.23.13
        fy3=1.6*fy3
!   CH4 forcing is 1.6 larger that given by fit
        if(inprnt.lt.25)then
        print *,' iter=',iter
        print *,' yo=',yo,' fy1=',fy1
        print *,' y=',y,' fy3=',fy3
        endif
        dfdy=fy3/(y-y0)
        df=fy1-fy3
        y=y+df/dfdy
        iter=iter+1
        if(iter.lt.15.and.abs(df).gt.1.e-5) go to 200
        yn=y
!      CO2 forcing
        c=co2o
        iter=0
        fc1=FCO2I(c,co20,xo,x0)
  300   continue
        fc3=FCO23(c)-frCO23
        dfdx=fc3/(c-co20)
        df=fc1-fc3
        c=c+df/dfdx
        iter=iter+1
        if(iter.lt.15.and.abs(df).gt.1.e-5) go to 300
        co2n=c
!      print *,'From NEWCON'
!      print *,'N2O ',x0,xo,xn
!      print *,'CH4 ',y0,yo,yn
        inprnt=inprnt+1
#endif
      return
      end
