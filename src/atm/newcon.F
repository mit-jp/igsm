
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

      subroutine newcon(xn,yn,xo,yo,x0,y0)

      FFN3(X,Y)=3.1515*X/(1.0+0.21448*X**0.87) - 0.056*X*Y              9131.5
     +          +0.3234*Y/(1.0+0.052*Y**0.84)
      FFN1(X,Y)=1.556*LOG(1.0+1.098*X**0.77*(1.0+0.032*X)               9113.3
     +                                       /(1.0+0.0014*X*X))         9113.4
     +          +(0.394*Y**0.66+0.16*Y*EXP(-1.6*Y))                     9113.5
     +                         /(1.0+0.169*Y**0.62)                     9113.6
     +          -0.14*LOG(1.0+0.636*(X*Y)**0.75+0.007*Y*(X*Y)**1.52)    9113.7
        data ifirst /1/
        if(ifirst.eq.1)then
c       print *,x0,y0,xo,yo
        fr1=ffn1(x0,y0)
        fr3=ffn3(x0,y0)
        ifirst=0
        inprnt=0
        endif
        fx1=ffn1(xo,y0)-fr1
        fy1=ffn1(x0,yo)-fr1
        x=xo
        iter=0
  100   continue
        fx3=ffn3(x,y0)-fr3
        dfdx=fx3/(x-x0)
        df=fx1-fx3
        x=x+df/dfdx
        iter=iter+1
        if(iter.lt.15.and.abs(df).gt.1.e-5) go to 100
        xn=x
c       if(inprnt.lt.25) then
c       print *,' iter=',iter
c       print *,' xo=',xo,' fx1=',fx1,' x=',x,' fx3=',fx3
c       endif
        iter=0
        y=yo
  200   continue
        fy3=ffn3(x0,y)-fr3
        dfdy=fy3/(y-y0)
        df=fy1-fy3
        y=y+df/dfdy
        iter=iter+1
        if(iter.lt.15.and.abs(df).gt.1.e-5) go to 200
        yn=y
c       if(inprnt.lt.25)then
c       print *,' iter=',iter
c       print *,' yo=',yo,' fy1=',fy1,' y=',y,' fy3=',fy3
c       endif
c       inprnt=inprnt+1
      return
      end
