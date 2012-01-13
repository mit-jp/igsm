
#include "ctrparam.h"
#if ( defined OLD_CHEM )

!	============================================================
!
!	CHEMTROP.F:	Subroutine for calculating chemical
!				reactions in troposphere
!				of MIT Global Chemistry Model
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
!	081898	Chien Wang	rev.	
!	080200	Chien Wang	repack based on CliChem3 & add cpp
!
!	==========================================================

!
	subroutine chemtrop (dtr, ifss, ktrop, Temp, qv, den)
!	=====================================================

!	==========================================================
!	Brief Description of Dummy Variables:
!
!	dtr:	time step of integration in second
!	ifss:	1 for using steady-state assumption
!		0 for not using
!	ktrop:  layer label of tropopause
!	Temp:	temperature in Kelvin
!	p:	air pressure in hPa
!	qv:	water vapor density ratio in kg/m^3
!	den:	air density in kg/m^3
!
!	==========================================================	

#include "chem_para"
#include "chem_const1"
#include "chem_com"

	dimension Temp (nlon,nlat,nlev)
	dimension qv   (nlon,nlat,nlev)
	dimension den  (nlon,nlat,nlev)

	dimension ychem    (nvaria)
	dimension awchem   (nvaria)
	dimension revawchem(nvaria)
	dimension ychem0   (2)

	real ychem_hfc134a

!	------------------------------------------------------------

#if ( defined CPL_CHEM )
	
	awchem(1) = awO3
	awchem(2) = awCO
	awchem(3) = awCO2
	awchem(4) = awHO
	awchem(5) = awHO2
	awchem(6) = awH2O2
	awchem(7) = awNO
	awchem(8) = awNO2
	awchem(9) = awCH4
	awchem(10)= awCH2O
	awchem(11)= awCH3O2H
	awchem(12)= awSO2
	awchem(13)= awN2O5
	awchem(14)= awHNO3
	awchem(15)= awH2SO4

!	write(6,*)"Here are solarflux & coszangle"
!	write(6,*)solarflux
!	write(6,*)coszangle
!	write(6,*)

       do ind = 1,15
         revawchem(ind) = 1./awchem(ind)
       enddo
       
       yyy       = Avogadro*1.e-3

       yyy1      = yyy/28.964 		!dry air
       yyy2      = yyy/18.0152		!water vapor

       do 1 k=1,ktrop
         patm = airpress(k)/1013.27
       do 2 j=1,n2dh

       Rflux     = solarflux(1,j,k)
     &           * 0.5678               !surface ratio of
                                        ! Act(<= 800nm)/Act(total)
                                        ! WMO 1981

	if(Rflux.le.50.0) goto 2           !no photochemistry
 					   !for nighttime

! === Convert air density from kg/m^3 to mols/cm^3:

       ychem0(1) = den(1,j,k)*yyy1	   	!dry air
       ychem0(2) = qv (1,j,k)*yyy2		!water vapor

! === Convert mixing ratio from ppb(m) to mols/cm^3:

       ychem( 1) = o3    (1,j,k)
       ychem( 2) = co    (1,j,k)
       ychem( 3) = zco2  (1,j,k)
       ychem( 4) = ho    (1,j,k)
       ychem( 5) = ho2   (1,j,k)
       ychem( 6) = h2o2  (1,j,k)
       ychem( 7) = xno   (1,j,k)
       ychem( 8) = xno2  (1,j,k)
       ychem( 9) = ch4   (1,j,k)
       ychem(10) = ch2o  (1,j,k)
       ychem(11) = ch3o2h(1,j,k)
       ychem(12) = so2   (1,j,k)

       ychem(13) = xn2o5 (1,j,k)
       ychem(14) = hno3  (1,j,k)
       ychem(15) = h2so4 (1,j,k)

       do iii = 1,nvaria
         ychem(iii) = max(0.0, ychem(iii))
       enddo

       call ppbm2mcm (nvaria, ychem, den(1,j,k), revawchem)

#ifdef INC_3GASES       
        ! === if hfc, pfc, and sf6 are included:
	! === 
	! === 033098 add hfc134a:
	  ychem_hfc134a = hfc134a (1,j,k)
	       
	  revawchemhfc134a = 1./awHFC134a
	  call ppbm2mcm(1,ychem_hfc134a,den(1,j,k),
     &		          revawchemhfc134a)
#endif

! === Calculate reaction rates:

        ntemp = nint(max(1.0,
     &               min(300.0,(Temp(1,j,k) - 200.)*2.0)))

! ===
! === 042596 add zenith angle function:

       call rateph (Temp(1,j,k), ntemp, Rflux, coszangle(1,j))

       call ratebm (Temp(1,j,k), ntemp, patm)

       call ratetm (Temp(1,j,k), ntemp, ychem0)

!       if(ifss.eq.1) call chemsteady(nvaria, ychem0,ychem)

! =====
! ===== New scheme for calculating NOy and S(VI)
! =====     Chien Wang, 092195
! =====

!
! === weighting d[N2O5] and d[HNO3]:
!
	xxxn2o5 = rk(14)*ychem(1)	!  k14[O3]
	xxxhno3 = rk(13)*ychem(4)	!  k13[HO]
	xxxtotal= 2*xxxn2o5 + xxxhno3	! [NO2] = 1./
					! (2xxxn2o5 + xxxhno3)
!
! === save initial NOx and SO2:
!
	ynox0 = ychem(7) + ychem(8)
	yso20 = ychem(12)

!	-------------------------------------------
!       call tropreact1(dtr, nvaria, ychem0, ychem)
       call tropreact1(dtr, 12, ychem0, ychem)
!	-------------------------------------------

!
! === Derive [N2O5], [HNO3], and [H2SO4]; 092195:
!
	dnox = min(0.0, ychem(7) + ychem(8) - ynox0)
	if(xxxtotal.gt.0.0)then
	ychem(13) = ychem(13) - xxxn2o5/xxxtotal * dnox 
	ychem(14) = ychem(14) - xxxhno3/xxxtotal * dnox
	endif

	ychem(15) = ychem(15) - min(0.0, ychem(12) - yso20)

#ifdef INC_3GASES
        ! === if hfc, pfc, and sf6 are included:
	! ===
	! === 033098 add hfc134a oxidation: CF3CH2F + OH -> CF3CHF + H2O
	! ===	here derive exp{-k[OH]t} first:

	  xxxhfc134a    = exp(- 1.7e-12 
     &                  * exp(-1750.0/Temp(1,j,k))
     &		        * dtr*ychem(4))
     	  ychem_hfc134a = max(0.0, ychem_hfc134a*xxxhfc134a)
	
	  call mcm2ppbm (1,ychem_hfc134a,den(1,j,k),awHFC134a)
	
	  hfc134a(1,j,k)= ychem_hfc134a
#endif
	
! === Convert concentration to ppb(m):

       call mcm2ppbm (nvaria, ychem, den(1,j,k), awchem)

       do 10 iii = 1,nvaria
         ychem(iii) = max(0.0, ychem(iii))
10     continue

! === Give value back to the common block:

       o3    (1,j,k) = ychem(01)
       co    (1,j,k) = ychem(02)
       zco2  (1,j,k) = ychem(03)
       ho    (1,j,k) = ychem(04)
       ho2   (1,j,k) = ychem(05)
       h2o2  (1,j,k) = ychem(06)
       xno   (1,j,k) = ychem(07)
       xno2  (1,j,k) = ychem(08)
       ch4   (1,j,k) = ychem(09)
       ch2o  (1,j,k) = ychem(10)
       ch3o2h(1,j,k) = ychem(11)
       so2   (1,j,k) = ychem(12)
       xn2o5 (1,j,k) = ychem(13)
       hno3  (1,j,k) = ychem(14)
       h2so4 (1,j,k) = ychem(15)

2      continue
1      continue

#endif

       return
        end

!
! --- won't need this portion for run excludes chemistry
!
#if ( defined CPL_CHEM )
 
c
      subroutine tropreact1 (tout, neq, c00, y)
c     ==============================================

c==================================================================c
c                                                                  c
c   TROPREACT1.F:	Subroutine for calculating chemical        c
c			  reactions in troposphere                 c
c                             of MIT Global Chemistry Model        c
c                         partly based on LSODES                   c
c        -------------------------------------------------         c
c   Author:		Chien Wang                                 c
c			MIT Joint Program on Science and Policy    c
c			    of Global Change                       c             
c   Last Revised on:    January 31, 1995                           c
c                                                                  c
c==================================================================c

      external fex, jex

      dimension y(neq), rwork(872), iwork(45)
      dimension c00(2)

      common /c00tmp/c00tmp(2)

c      data lrw/500/, liw/30/

      c00tmp(1) = c00(1)
      c00tmp(2) = c00(2)

      t = 0.
      itol = 1
      rtol = 1.e-3
      atol = 1.e+3	!molecules/cm3

c     atol(1) = 1.e-6
c     atol(2) = 1.e-10
c     atol(3) = 1.e-6

      itask  = 1
      istate = 1
      iopt   = 0
      lrw    = 872
      liw    = 45
      mf     = 21

        call lsodenew(fex,neq,y,t,tout,itol,rtol,atol,itask,istate,
     &     iopt,rwork,lrw,iwork,liw,jex,mf)
c        call lsode(fex,neq,y,t,tout,itol,rtol,atol,itask,istate,
c     &     iopt,rwork,lrw,iwork,liw,jex,mf)

      return
       end

c
      subroutine fex (neq, t, y, ydot)
c     ================================
c
c  Subroutine for defining of ydot
c
     
#include "chem_para"
#include "chem_com"
      
      common /c00tmp/c00(2)

      dimension y(neq), ydot(neq)

        bbb = 0.7809*c00(1)
        ccc = 0.2095*c00(1)

c === redesigned scheme, 092195
c

	xbeta1 = rk(3)*bbb + rk(4)*ccc	!k3[N2] + k4[O2]
	xbeta2 = rk(2)*c00(2)		!k2[H2O]
	xbeta  = xbeta2/(xbeta1 + xbeta2)

	rr1 = rk(1) *y(1)		!j1[O3]
	rr5 = rk(5) *y(2) *y(4)		!k5[CO][HO]
	rr7 = rk(7) *y(5) *y(7)		!k7[HO2][NO]
	rr8 = rk(8) *y(8)		!k8[NO2]
	rr10= rk(10)*y(1) *y(5)		!k10[O3][HO2]
	rr11= rk(11)*y(1) *y(4)		!k11[O3][HO]
	rr12= rk(12)*y(1) *y(7)		!k12[O3][NO]
	rr13= rk(13)*y(4) *y(8)		!k13[HO][NO2]
	rr14= rk(14)*y(1) *y(8)		!k14[O3][NO2]
	rr16= rk(16)*y(5) *y(5)		!k16[HO2][HO2]
	rr17= rk(17)*y(6)		!j17[H2O2]
	rr18= rk(18)*y(4) *y(6)		!k18[HO][H2O2]
	rr19= rk(19)*y(4) *y(9)		!k19[HO][CH4]
	rr21= rk(21)*y(7)		!k21[NO] - for ss only
	rr23= rk(23)*y(5)		!k23[HO2]- for ss only
	rr24= rk(24)*y(11)		!j24[CH3O2H]
	rr25= rk(25)*y(4) *y(11)	!k25[HO][CH3O2H]
	rr26= rk(26)*y(10)		!j26[CH2O]
	rr27= rk(27)*y(4) *y(10)	!k27[HO][CH2O]
	rr29= rk(29)*y(4) *y(12)	!k29[HO][SO2]

c New reactions - 062995:

	rr32= rk(32)*y(4)*y(5)		!k32[HO][HO2]
	rr33= rk(33)*y(4)*y(4)		!k33[HO][HO]
	rr34= rk(34)*y(4)*y(4)		!k34[HO][HO]

	xxx   = (rr21 + rr23)

	gamax = 0.0
	if(rr21.gt.1.e-12)gamax = rr21/xxx

	gamay = 0.0
	if(xxx.gt.0.0)gamay = rr23/xxx

      ydot(1)  = rr8 + rr33				!O3
     &	       -(xbeta*rr1 + rr10 + rr11 + rr12 + rr14)		

      ydot(2)  = rr26 + rr27				!CO
     &         - rr5

      ydot(3)  = rr5 					!CO2

      ydot(4)  = 2.0*xbeta*rr1 + rr7 + rr10 		!HO
     &	       + 2.0*rr17+ rr24
     &         -(rr5 + rr11 + rr13 + rr18 + rr19 
     &          +rr25+ rr27 + rr29
     &	        +rr32+ rr33 + rr34 )

      ydot(5)  = rr5 + rr11 + rr18 + rr27 + rr29	!HO2
     &         + rr24+ 2.0*rr26
     &         +(2.0*gamax - 1.0)*(rr19 + rr25)
     &         -(rr7 + rr10 + rr16 + rr32)

      ydot(6)  = rr16 + rr34				!H2O2
     &         -(rr17+ rr18)

      ydot(7)  = rr8 					!NO
     &	       -(rr7 + rr12 + gamax*(rr19 + rr25))

      ydot(8)  = rr7 + rr12 + gamax*(rr19 + rr25)	!NO2
     &	       -(rr8 + rr13+ rr14+ rr14)

      ydot(9)  =-rr19					!CH4

      ydot(10) = rr24 + gamax*(rr19 + rr25)		!CH2O
     &	       -(rr26 + rr27)

      ydot(11) = gamay*(rr19 + rr25)			!CH3O2H
     &	       -(rr24 + rr25) 

      ydot(12) =-rr29 					!SO2

      return
       end

c
      subroutine jex (neq, t, y, ml, mu, pd, nrpd)
c     ==========================================
c
c  Subroutine for calculating the Jacobian for f(i)
c
c     P(m,n) here represents df(m)/dn
c                         or d ydot(m)/dn
c

#include "chem_para"
#include "chem_com"

      common /c00tmp/c00(2)

      dimension y(neq), pd(nrpd,neq)

      bbb = 0.7809*c00(1)
      ccc = 0.2095*c00(1)

        xbeta1 = rk(3)*bbb + rk(4)*ccc   !k3[N2] + k4[O2]
        xbeta2 = rk(2)*c00(2)            !k2[H2O]
        xbeta  = xbeta2/(xbeta1 + xbeta2)

        rr21= rk(21)*y(7)               !k21[NO] - for ss only
        rr23= rk(23)*y(5)               !k23[HO2]- for ss only
        xxx   = (rr21 + rr23)

        gamax = 0.0
        if(rr21.gt.1.e-12)gamax = rr21/xxx

        gamay = 0.0
        if(xxx.gt.0.0)gamay = rr23/xxx

      do 1 j=1,neq
      do 1 i=1,neq
        pd(i,j) = 0.0
1     continue

c =====
c ===== pd(i,j) = d ydot(i)/dy(j):
c =====
	
	rk5y2    = rk(5)*y(2)
	rk5y4    = rk(5)*y(4)
	rk7y5    = rk(7)*y(5)
	rk7y7    = rk(7)*y(7)
	rk10y1   = rk(10)*y(1)
	rk10y5   = rk(10)*y(5)
	rk11y1   = rk(11)*y(1)
	rk11y4   = rk(11)*y(4)
	rk12y1   = rk(12)*y(1)
	rk12y7   = rk(12)*y(7)
	rk13y4   = rk(13)*y(4)
	rk13y8   = rk(13)*y(8)
	rk14y1   = rk(14)*y(1)
	rk14y8   = rk(14)*y(8)
	rk18y4   = rk(18)*y(4)
	rk18y6   = rk(18)*y(6)
	rk19y4   = rk(19)*y(4)
	rk19y9   = rk(19)*y(9)
	rk25y4   = rk(25)*y(4)
	rk25y11  = rk(25)*y(11)
	rk27y4   = rk(27)*y(4)
	rk27y10  = rk(27)*y(10)
	rk29y4   = rk(29)*y(4)
	rk29y12  = rk(29)*y(12)
	rk32y4   = rk(32)*y(4)
	rk32y5   = rk(32)*y(5)

	pd(1,1)  = -(xbeta*rk(1)   + rk10y5 + rk11y4
     &	                  +rk12y7  + rk14y8)
	pd(4,1)  = 2.0*xbeta*rk(1) + rk10y5 - rk11y4
	pd(5,1)  = rk11y4 - rk10y5
	pd(8,1)  = rk12y7 - 2.*rk14y8

	pd(2,2)  =-rk5y4
	pd(3,2)  = rk5y4
	pd(4,2)  =-rk5y4
	pd(5,2)  = rk5y4

	pd(1,4)  = 2.0*rk(33)*y(4) - rk11y1
	pd(2,4)  = rk27y10    - rk5y2
	pd(3,4)  = rk5y2
	pd(4,4)  =-( 
     &	 	   rk5y2  + rk11y1 + rk(13)*y(8)
     &		 + rk18y6 + rk19y9 + rk25y11
     &		 + rk27y10+ rk29y12+ rk32y5
     &		 +(rk(33) + rk(34))*y(4)
     &              )
	pd(5,4)  =(rk5y2 + rk11y1 + rk18y6
     &           + rk27y10+ rk29y12)
     &           +(2.0*gamax - 1.0)*(rk19y9 + rk25y11)
     &		 - rk32y5
	pd(6,4)  =-rk18y6
	pd(7,4)  =-gamax*(rk19y9 + rk25y11)
	pd(8,4)  = gamax*(rk19y9 + rk25y11)
     &		 - rk13y8
	pd(9,4)  =-rk19y9
	pd(10,4) = gamax*(rk19y9 + rk25y11)
     &		 - rk27y10
	pd(11,4) = gamay*(rk19y9 + rk25y11)
     &		 - rk25y11
	pd(12,4) =-rk29y12

	pd(1,5)  = rk11y1
	pd(4,5)  = rk7y7 + rk10y1
     &           - rk32y4
	pd(5,5)  =
     &		 -(rk7y7 + rk10y1 + rk(16)*y(5) + rk32y4)
	pd(6,5)  = 2.0*   (rk(16) + rk(34))*y(5)
	pd(7,5)  =-rk7y7
	pd(8,5)  = rk7y7
c 	pd(11,5)

	pd(4,6)  = 2.0*rk(17)
     &           - rk18y4
	pd(5,6)  = rk18y4
	pd(6,6)  =-(rk(17) + rk18y4)

	pd(1,7)  =-rk12y1
	pd(4,7)  = rk7y5
	pd(5,7)  =-rk7y5
	pd(7,7)  =-(rk7y5 + rk12y1)
	pd(8,7)  = (rk7y5 + rk12y1)
	
	pd(1,8)  = rk(8)
     &		 - rk14y1
	pd(4,8)  =-rk13y4
	pd(7,8)  = rk(8)
	pd(8,8)  =-(rk(8) + rk13y4 + 2.0*rk14y1)

	pd(4,9)  =-rk19y4
	pd(5,9)  =(2.0*gamax - 1.0)*rk19y4
	pd(7,9)  =-gamax*rk19y4
	pd(8,9)  = gamax*rk19y4
	pd(9,9)  =-rk19y4
	pd(10,9) = gamax*rk19y4
	pd(11,9) = gamay*rk19y4

	pd(2,10) = rk(26) + rk27y4
	pd(4,10) =-rk27y4
	pd(5,10) = rk27y4 + 2.0*rk(26)
	pd(10,10)=-(rk(26)+ rk27y4)

	pd(4,11) = rk(24) - rk25y4
	pd(5,11) = rk(24) + (2.0*gamax - 1.0)*rk25y4
	pd(7,11) =-gamax*rk25y4
	pd(8,11) = gamax*rk25y4
	pd(10,11)= rk(24) + gamax*rk25y4
	pd(11,11)=-(rk(24) + rk25y4*max(0.0,1.0 - gamay))

	pd(4,12) =-rk(29)*y(4)
	pd(5,12) = rk29y4
	pd(12,12)=-rk29y4

      return
      end
 

c==================================================================c
c                                                                  c
c   CHEMFUNCPACK.F:  Subroutine and function pack                  c
c                             of MIT Global Chemistry Model        c
c        -------------------------------------------------         c
c   Author:		Chien Wang                                 c
c			MIT Joint Program on Science and Policy    c
c			    of Global Change                       c             
c   Last Revised on:    February 22, 1995                          c
c                                                                  c
c==================================================================c

c
      subroutine ppbm2mcm (np, xx, den, revaw) 
c     ==========================================

c------------------------------------------------
c   A program for convert mixing ratio in
c     ppb(m) to concentration in molecules/cm^3
c------------------------------------------------

      dimension xx    (np)
      dimension revaw (np)

      ddd = 6.02217e+11*den
      do 1 i=1,np
        xx (i) = xx (i)
     &          *ddd*revaw(i)
1     continue

      return
       end

c
      subroutine mcm2ppbm (np, xx, den, aw)
c     ==========================================

c------------------------------------------------
c   A program for convert concentration in
c     molecules/cm^3 to mixing ratio in ppb(m)
c------------------------------------------------

      dimension xx (np)
      dimension aw (np)

      ddd = 1./(6.02217e+11*den)
      do 1 i=1,np
        xx (i) = xx (i)
     &          *ddd*aw(i)
1     continue

      return
       end

c
      subroutine rateph (T, ntemp, Rflux, coszagcm)
c     ===========================

c------------------------------------------------
c   A program for calculating rate constants
c     of photochemical reactions
c         ----------------------------
c   The time unit for all reaction rates is second
c
c       New version revised at June 3, 1995
c
c------------------------------------------------

#include "chem_para"
#include "chem_com"
     
	coszagcm = max(cosza4rk(10),
     &             min(cosza4rk(1),coszagcm))

	do 1 nband=1,9
	  if(coszagcm.le.cosza4rk(nband).and.
     &       coszagcm.gt.cosza4rk(nband+1))then

	    wgtxx = 1.0/(cosza4rk(nband) - cosza4rk(nband+1))
	    wgt1  = (cosza4rk(nband) - coszagcm)  *wgtxx
	    wgt2  = (coszagcm - cosza4rk(nband+1))*wgtxx
c
c	+--- wgt2 ---+--- wgt1 ---+
c	nband --- coszagcm --- nband+1

c-- R8: NO2 + hv -> NO + O
c   Atkinson et al., 1992:
c
c        rk(8) = (18096.7314929
c     &          -3.9704*(T-273.))
c     &        *5.0e-10
c     &        *Rflux

c        rk(8) = rktable1(08,ntemp)*Rflux

        rk(8) = ((rk08gama(nband)*wgt2 + rk08gama(nband+1)*wgt1)
     &        +  (rk08aaa (nband)*wgt2 + rk08aaa (nband+1)*wgt1)
     &          *(T-273.))
     &        *Rflux

c-- R1: O3 + hv -> O(1D) + O2
c   Atkinson et al., 1992:
c
c        rk(1)  = 0.0028*rk(8)
cc       rk(1) = 3.467625419939578E-08*Rflux

c	rk(1) = (rk01table(nband)*wgt2 
c     &         + rk01table(nband+1)*wgt1)*Rflux
        rk(1)  = 0.0028*rk(8)

c-- R17: H2O2 + hv -> 2OH
c   Atkinson et al., 1992:
c
c       rk(17)= 6.263464249748237E-09*Rflux

	rk(17) = (rk17table(nband)*wgt2 
     &         +  rk17table(nband+1)*wgt1)*Rflux

c-- R24: CH3O2H + hv -> CH3O + HO
c   Atkinson et al., 1992:
c
c       rk(24)= 4.474208962739173E-09*Rflux

	rk(24) = (rk24table(nband)*wgt2 
     &         +  rk24table(nband+1)*wgt1)*Rflux

c-- R26: CH2O + hv -> CHO + H
c   Atkinson et al., 1992:
c
c       rk(26)= 9.410107812688823E-08*Rflux

	rk(26) = (rk26table(nband)*wgt2 
     &         +  rk26table(nband+1)*wgt1)*Rflux

	goto 2
	endif

1	continue

c = For znith angle.ge.86:

        rk(8) = (rk08gama(10)
     &        +  rk08aaa (10)
     &          *(T-273.))
     &        *Rflux

c	rk(1) = rk01table(10)
c     &        * Rflux
        rk(1)  = 0.0028*rk(8)

	rk(17) = rk17table(10)
     &         * Rflux

	rk(24) = rk24table(10)
     &         * Rflux

	rk(26) = rk26table(10)
     &         *Rflux

2	continue

      return
       end

c
      subroutine ratebm (T, ntemp, patm)
c     ===========================

c------------------------------------------------
c   A program for calculating rate constants
c     of bimolecular reactions
c         ----------------------------
c   The time unit for all reaction rates is second
c
c------------------------------------------------

#include "chem_para"
#include "chem_com"

	Trev = 1./T

c-- R2: O(1D) + H2O -> 2OH

	rk(2)  = 2.2e-10

c-- R3: O(1D) + N2 -> O + N2

c	rk(3)  = 1.8e-11 * exp( 107.0*Trev)
	rk(3)  = rktable1(03,ntemp)

c-- R4: O(1D) + O2 -> O + O2

c	rk(4)  = 3.2e-11 * exp(  67.0*Trev)
	rk(4)  = rktable1(04,ntemp)

c-- R5: CO + HO -> H + CO2

	rk(5)  = 1.5e-13 * (1.0 + 0.6*patm)

c-- R7: HO2 + NO -> HO + NO2

c	rk(7)  = 3.7e-12 * exp( 240.0*Trev)
	rk(7)  = rktable1(07,ntemp)

c-- R10: HO2 + O3 -> HO + 2O2

c	rk(10) = 1.1e-14 * exp(-500.0*Trev)
	rk(10) = rktable1(10,ntemp)

c-- R11: HO + O3 -> HO2 + O2

c	rk(11) = 1.6e-12 * exp(-940.0*Trev)
	rk(11) = rktable1(11,ntemp)
c	rk(11) = 0.0

c-- R12: NO + O3 -> NO2 + O2

c	rk(12) = 2.0e-12 * exp(-1400.0*Trev)
	rk(12) = rktable1(12,ntemp)

c-- R14: NO2 + O3 -> NO3 + O2

c	rk(14) = 1.2e-13 * exp(-2450.0*Trev)
	rk(14) = rktable1(14,ntemp)

c-- R16: HO2 + HO2 -> H2O2 + O2

c	rk(16) = 2.3e-13 * exp( 600.0*Trev)

c-- R18: H2O2 + HO -> HO2 + H2O

c	rk(18) = 2.9e-12 * exp(-160.0*Trev)
	rk(18) = rktable1(18,ntemp)

c-- R19: CH4 + HO -> CH3 + H2O

c	rk(19) = 2.65e-12 * exp(-1800.0*Trev)
	rk(19) = rktable1(19,ntemp)

c-- R21: CH3O2 + NO -> CH3O + NO2

c	rk(21) = 4.2e-12 * exp( 180.0*Trev)
	rk(21) = rktable1(21,ntemp)

c-- R22: CH3O + O2 -> CH2O + HO2

c	rk(22) = 3.9e-14 * exp(-900.0*Trev)
	rk(22) = rktable1(22,ntemp)

c-- R23: CH3O2 + HO2 -> CH3O2H + O2

c	rk(23) = 3.8e-13 * exp( 780.0*Trev)
	rk(23) = rktable1(23,ntemp)

c-- R25: CH3O2H + HO -> CH3O2 + H2O

c	rk(25) = 1.9e-12 * exp( 190.0*Trev)
	rk(25) = rktable1(25,ntemp)

c-- R27: CH2O + HO -> CHO + H2O

	rk(27) = 1.0e-11

c-- R28: CHO + O2 -> CO + HO2

c	rk(28) = 3.5e-12 * exp( 140.0*Trev)
	rk(28) = rktable1(28,ntemp)

c-- R30: HOSO2 + O2 -> HO2 + SO3

c	rk(30) = 1.3e-12 * exp(-330.0*Trev)
	rk(30) = rktable1(30,ntemp)

c-- R31: SO3 + H2O -> H2SO4

	rk(31) = 2.4e-15

c New reactions - 062995:

c-- R32: HO + HO2 -> H2O + O2

c	rk(32) = 4.8e-11 * exp( 250.0*Trev)
	rk(32) = rktable1(32,ntemp)

c-- R33: HO + HO -> H2O + O

c	rk(33) = 4.2e-12 * exp(-240.0*Trev)
	rk(33) = rktable1(33,ntemp)

      return
       end

c
      subroutine ratetm (T, ntemp, ym)
c     =========================

c------------------------------------------------
c   A program for calculating rate constants
c     of termolecular reactions
c         ----------------------------
c   The time unit for all reaction rates is second
c
c------------------------------------------------

#include "chem_para"
#include "chem_com"

        dimension ym(2)

        T300 = T/300.
        xm   = ym(1)
        xn2  = 0.7809*xm
        h2o  = ym(2)

c-- R6:  H + O2 + M -> HO2 + M

c	rkt0   = 6.2e-32 * xm * T300**(-1.6)
	rkt0   = rktable1(06,ntemp) * xm
	rkt00  = 7.5e-11
	rk(6)  = calkmt(rkt0, rkt00)

c-- R9:  O + O2 + M -> O3 + M  (no high pres. limit)

c	rkt0   = 5.6e-34 * xm * T300**(-2.8)
	rkt0   = rktable1(09,ntemp) * xm
	rk(9)  = rkt0

c-- R13: NO2 + HO + M -> HNO3 + M

c	rkt0   = 2.6e-30 * xm * T300**(-3.2)
c	rkt00  = 2.4e-11 * T300**(-1.3)
	rkt0   = rktable1(13,ntemp) * xm
	rkt00  = rktable2(1, ntemp)
	rk(13) = calkmt(rkt0, rkt00)

c-- R15: NO3 + NO2 + M -> N2O5 + M

c	rkt0   = 2.2e-30 * xm * T300**(-3.9)
c	rkt00  = 1.5e-12 * T300**(-0.7)
	rkt0   = rktable1(15,ntemp) * xm
	rkt00  = rktable2(2, ntemp)
	rk(15) = calkmt(rkt0, rkt00)

c-- R16: HO2 + HO2 + M -> H2O2 + O2
c  062895:

c	rk(16) = 1.7e-33 * xn2 * exp(1000.0/T)
	rk(16) = rktable1(16,ntemp) * xn2
     &         *(1.+ 1.4e-21*h2o*exp(2200.0/T))

c-- R20: CH3 + O2 + M -> CH3O2 + M

c	rkt0   = 4.5e-31 * xm * T300**(-3.0)
c	rkt00  = 1.8e-12 * T300**(-1.7)
	rkt0   = rktable1(20,ntemp) * xm
	rkt00  = rktable2(3, ntemp)
	rk(20) = calkmt(rkt0, rkt00)

c-- R29: SO2 + HO + M -> HOSO2 + M

c	rkt0   = 3.0e-31 * xm * T300**(-3.3)
	rkt0   = rktable1(29,ntemp) * xm
	rkt00  = 1.5e-12
	rk(29) = calkmt(rkt0, rkt00)

c New reactions, 063095:

c-- R34: HO + HO + M -> H2O2 + M

c        rkt0   = 6.9e-31 * xm * T300**(-0.8)
        rkt0   = rktable1(34,ntemp) * xm
        rkt00  = 1.5e-11
        rk(34) = calkmt(rkt0, rkt00)


      return
       end

c
      function calkmt (rkt0, rkt00)
c     =============================

c
c  A function for calculate k(M,T) of termolecular
c    reaction rate from low and high limit
c

      aaa = rkt0/rkt00
      bbb = log10(aaa)**2
      ccc = 1./(1.+bbb)

      calkmt = rkt0/(1.+aaa) * 0.6**ccc

      return
       end


c
	Block data rktable3dat
c	======================

c----------------------------------------c
c  Block data for holding rktable3 data  c
c----------------------------------------c

#include "chem_para"
#include "chem_com"

	data cosza4rk/
     &      0.100000E+01,  0.984808E+00,  0.939693E+00,  0.866025E+00,
     &      0.766044E+00,  0.642788E+00,  0.500000E+00,  0.342020E+00,
     &      0.207912E+00,  0.697565E-01/

	data rk08gama/
     &      0.958260E-05,  0.956540E-05,  0.947786E-05,  0.933214E-05,
     &      0.904837E-05,  0.853986E-05,  0.790570E-05,  0.644646E-05,
     &      0.540693E-05,  0.640000E-05/

	data rk08aaa/
     &     -0.213357E-08, -0.212833E-08, -0.210322E-08, -0.206239E-08,
     &     -0.198520E-08, -0.185147E-08, -0.169861E-08, -0.133796E-08,
     &     -0.109912E-08, -0.131294E-08/

	data rk01table/
     &      0.553973E-07,  0.541612E-07,  0.501300E-07,  0.436050E-07,
     &      0.346763E-07,  0.241176E-07,  0.141663E-07,  0.593145E-08,
     &      0.238321E-08,  0.148548E-08/

	data rk17table/
     &      0.761430E-08,  0.755268E-08,  0.730677E-08,  0.691454E-08,
     &      0.626346E-08,  0.536943E-08,  0.430112E-08,  0.288139E-08,
     &      0.199044E-08,  0.209279E-08/

	data rk24table/
     &      0.535338E-08,  0.531281E-08,  0.515418E-08,  0.489785E-08,
     &      0.447421E-08,  0.387878E-08,  0.318768E-08,  0.217960E-08,
     &      0.155503E-08,  0.170283E-08/

	data rk26table/
     &      0.109106E-06,  0.108473E-06,  0.105771E-06,  0.101587E-06,
     &      0.941011E-07,  0.831775E-07,  0.701839E-07,  0.494384E-07,
     &      0.364853E-07,  0.410051E-07/

	end

#endif

#endif

#if ( !defined OLD_CHEM )

!	============================================================
!
!	CHEMTROP.F:	Subroutine for calculating chemical
!				reactions in troposphere
!				of MIT Global Chemistry Model
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
!	081898	Chien Wang	rev.	
!	080200	Chien Wang	repack based on CliChem3 & add cpp
!
!	==========================================================

!
	subroutine chemtrop (dtr, ifss, ktrop, Temp, qv, den)
!	=====================================================

!	==========================================================
!	Brief Description of Dummy Variables:
!
!	dtr:	time step of integration in second
!	ifss:	1 for using steady-state assumption
!		0 for not using
!	ktrop:  layer label of tropopause
!	Temp:	temperature in Kelvin
!	p:	air pressure in hPa
!	qv:	water vapor density ratio in kg/m^3
!	den:	air density in kg/m^3
!
!	==========================================================	

#include "chem_para"
#include "chem_const1"
#include "chem_com"

	dimension Temp (nlon,nlat,nlev)
	dimension qv   (nlon,nlat,nlev)
	dimension den  (nlon,nlat,nlev)

	dimension ychem    (nvaria)
	dimension awchem   (nvaria)
	dimension revawchem(nvaria)
	dimension ychem0   (2)

	real ychem_hfc134a

!	------------------------------------------------------------

#if ( defined CPL_CHEM )
	
	awchem(1) = awO3
	awchem(2) = awCO
	awchem(3) = awCO2
	awchem(4) = awHO
	awchem(5) = awHO2
	awchem(6) = awH2O2
	awchem(7) = awNO
	awchem(8) = awNO2
	awchem(9) = awCH4
	awchem(10)= awCH2O
	awchem(11)= awCH3O2H
	awchem(12)= awSO2
	awchem(13)= awN2O5
	awchem(14)= awHNO3
	awchem(15)= awH2SO4

!	write(6,*)"Here are solarflux & coszangle"
!	write(6,*)solarflux
!	write(6,*)coszangle
!	write(6,*)

       do ind = 1,15
         revawchem(ind) = 1./awchem(ind)
       enddo
       
       yyy       = Avogadro*1.e-3

       yyy1      = yyy/28.964 		!dry air
       yyy2      = yyy/18.0152		!water vapor

       do 1 k=1,ktrop
         patm = airpress(k)/1013.27
       do 2 j=1,n2dh

       Rflux     = solarflux(1,j,k)
     &           * 0.5678               !surface ratio of
                                        ! Act(<= 800nm)/Act(total)
                                        ! WMO 1981

	if(Rflux.le.50.0) goto 2           !no photochemistry
 					   !for nighttime

! === Convert air density from kg/m^3 to mols/cm^3:

       ychem0(1) = den(1,j,k)*yyy1	   	!dry air
       ychem0(2) = qv (1,j,k)*yyy2		!water vapor

! === Convert mixing ratio from ppb(m) to mols/cm^3:

       ychem( 1) = o3    (1,j,k)
       ychem( 2) = co    (1,j,k)
       ychem( 3) = zco2  (1,j,k)
       ychem( 4) = ho    (1,j,k)
       ychem( 5) = ho2   (1,j,k)
       ychem( 6) = h2o2  (1,j,k)
       ychem( 7) = xno   (1,j,k)
       ychem( 8) = xno2  (1,j,k)
       ychem( 9) = ch4   (1,j,k)
       ychem(10) = ch2o  (1,j,k)
       ychem(11) = ch3o2h(1,j,k)
       ychem(12) = so2   (1,j,k)

       ychem(13) = xn2o5 (1,j,k)
       ychem(14) = hno3  (1,j,k)
       ychem(15) = h2so4 (1,j,k)

       do iii = 1,nvaria
         ychem(iii) = max(0.0, ychem(iii))
       enddo

       call ppbm2mcm (nvaria, ychem, den(1,j,k), revawchem)

#ifdef INC_3GASES       
        ! === if hfc, pfc, and sf6 are included:
	! === 
	! === 033098 add hfc134a:
	  ychem_hfc134a = hfc134a (1,j,k)
	       
	  revawchemhfc134a = 1./awHFC134a
	  call ppbm2mcm(1,ychem_hfc134a,den(1,j,k),
     &		          revawchemhfc134a)
#endif

! === Calculate reaction rates:

        ntemp = nint(max(1.0,
     &               min(300.0,(Temp(1,j,k) - 200.)*2.0)))

! ===
! === 042596 add zenith angle function:

       call rateph (Temp(1,j,k), ntemp, Rflux, coszangle(1,j))

       call ratebm (Temp(1,j,k), ntemp, patm)

       call ratetm (Temp(1,j,k), ntemp, ychem0)

!       if(ifss.eq.1) call chemsteady(nvaria, ychem0,ychem)

! =====
! ===== New scheme for calculating NOy and S(VI)
! =====     Chien Wang, 092195
! =====

!
! === weighting d[N2O5] and d[HNO3]:
!
	xxxn2o5 = rk(14)*ychem(1)	!  k14[O3]
	xxxhno3 = rk(13)*ychem(4)	!  k13[HO]
	xxxtotal= 2*xxxn2o5 + xxxhno3	! [NO2] = 1./
					! (2xxxn2o5 + xxxhno3)
!
! === save initial NOx and SO2:
!
	ynox0 = ychem(7) + ychem(8)
	yso20 = ychem(12)

!	-------------------------------------------
!       call tropreact1(dtr, nvaria, ychem0, ychem)
       call tropreact1(dtr, 12, ychem0, ychem)
!	-------------------------------------------

!
! === Derive [N2O5], [HNO3], and [H2SO4]; 092195:
!
	dnox = min(0.0, ychem(7) + ychem(8) - ynox0)
	if(xxxtotal.gt.0.0)then
	ychem(13) = ychem(13) - xxxn2o5/xxxtotal * dnox 
	ychem(14) = ychem(14) - xxxhno3/xxxtotal * dnox
	endif

	ychem(15) = ychem(15) - min(0.0, ychem(12) - yso20)

#ifdef INC_3GASES
        ! === if hfc, pfc, and sf6 are included:
	! ===
	! === 033098 add hfc134a oxidation: CF3CH2F + OH -> CF3CHF + H2O
	! ===	here derive exp{-k[OH]t} first:

	  xxxhfc134a    = exp(- 1.7e-12 
     &                  * exp(-1750.0/Temp(1,j,k))
     &		        * dtr*ychem(4))
     	  ychem_hfc134a = max(0.0, ychem_hfc134a*xxxhfc134a)
	
	  call mcm2ppbm (1,ychem_hfc134a,den(1,j,k),awHFC134a)
	
	  hfc134a(1,j,k)= ychem_hfc134a
#endif
	
! === Convert concentration to ppb(m):

       call mcm2ppbm (nvaria, ychem, den(1,j,k), awchem)

       do 10 iii = 1,nvaria
         ychem(iii) = max(0.0, ychem(iii))
10     continue

! === Give value back to the common block:

       o3    (1,j,k) = ychem(01)
       co    (1,j,k) = ychem(02)
       zco2  (1,j,k) = ychem(03)
       ho    (1,j,k) = ychem(04)
       ho2   (1,j,k) = ychem(05)
       h2o2  (1,j,k) = ychem(06)
       xno   (1,j,k) = ychem(07)
       xno2  (1,j,k) = ychem(08)
       ch4   (1,j,k) = ychem(09)
       ch2o  (1,j,k) = ychem(10)
       ch3o2h(1,j,k) = ychem(11)
       so2   (1,j,k) = ychem(12)
       xn2o5 (1,j,k) = ychem(13)
       hno3  (1,j,k) = ychem(14)
       h2so4 (1,j,k) = ychem(15)

2      continue
1      continue

#endif

       return
        end

!
! --- won't need this portion for run excludes chemistry
!
#if ( defined CPL_CHEM )
 
c
      subroutine tropreact1 (tout, neq, c00, y)
c     ==============================================

c==================================================================c
c                                                                  c
c   TROPREACT1.F:	Subroutine for calculating chemical        c
c			  reactions in troposphere                 c
c                             of MIT Global Chemistry Model        c
c                         partly based on LSODES                   c
c        -------------------------------------------------         c
c   Author:		Chien Wang                                 c
c			MIT Joint Program on Science and Policy    c
c			    of Global Change                       c             
c   Last Revised on:    January 31, 1995                           c
c                                                                  c
c==================================================================c

      external fex, jex

      dimension y(neq), rwork(872), iwork(45)
      dimension c00(2)

      common /c00tmp/c00tmp(2)

c      data lrw/500/, liw/30/

      c00tmp(1) = c00(1)
      c00tmp(2) = c00(2)

      t = 0.
      itol = 1
      rtol = 1.e-3
      atol = 1.e+3	!molecules/cm3

c     atol(1) = 1.e-6
c     atol(2) = 1.e-10
c     atol(3) = 1.e-6

      itask  = 1
      istate = 1
      iopt   = 0
      lrw    = 872
      liw    = 45
      mf     = 21

        call lsodenew(fex,neq,y,t,tout,itol,rtol,atol,itask,istate,
     &     iopt,rwork,lrw,iwork,liw,jex,mf)
c        call lsode(fex,neq,y,t,tout,itol,rtol,atol,itask,istate,
c     &     iopt,rwork,lrw,iwork,liw,jex,mf)

      return
       end

c
      subroutine fex (neq, t, y, ydot)
c     ================================
c
c  Subroutine for defining of ydot
c
     
#include "chem_para"
#include "chem_com"
      
      common /c00tmp/c00(2)

      dimension y(neq), ydot(neq)

        bbb = 0.7809*c00(1)
        ccc = 0.2095*c00(1)

c === redesigned scheme, 092195
c

	xbeta1 = rk(3)*bbb + rk(4)*ccc	!k3[N2] + k4[O2]
	xbeta2 = rk(2)*c00(2)		!k2[H2O]
	xbeta  = xbeta2/(xbeta1 + xbeta2)

	rr1 = rk(1) *y(1)		!j1[O3]
	rr5 = rk(5) *y(2) *y(4)		!k5[CO][HO]
	rr7 = rk(7) *y(5) *y(7)		!k7[HO2][NO]
	rr8 = rk(8) *y(8)		!k8[NO2]
	rr10= rk(10)*y(1) *y(5)		!k10[O3][HO2]
	rr11= rk(11)*y(1) *y(4)		!k11[O3][HO]
	rr12= rk(12)*y(1) *y(7)		!k12[O3][NO]
	rr13= rk(13)*y(4) *y(8)		!k13[HO][NO2]
	rr14= rk(14)*y(1) *y(8)		!k14[O3][NO2]
	rr16= rk(16)*y(5) *y(5)		!k16[HO2][HO2]
	rr17= rk(17)*y(6)		!j17[H2O2]
	rr18= rk(18)*y(4) *y(6)		!k18[HO][H2O2]
	rr19= rk(19)*y(4) *y(9)		!k19[HO][CH4]
	rr21= rk(21)*y(7)		!k21[NO] - for ss only
	rr23= rk(23)*y(5)		!k23[HO2]- for ss only
	rr24= rk(24)*y(11)		!j24[CH3O2H]
	rr25= rk(25)*y(4) *y(11)	!k25[HO][CH3O2H]
	rr26= rk(26)*y(10)		!j26[CH2O]
	rr27= rk(27)*y(4) *y(10)	!k27[HO][CH2O]
	rr29= rk(29)*y(4) *y(12)	!k29[HO][SO2]

c New reactions - 062995:

	rr32= rk(32)*y(4)*y(5)		!k32[HO][HO2]
	rr33= rk(33)*y(4)*y(4)		!k33[HO][HO]
	rr34= rk(34)*y(4)*y(4)		!k34[HO][HO]

	xxx   = (rr21 + rr23)

	gamax = 0.0
	if(rr21.gt.1.e-12)gamax = rr21/xxx

	gamay = 0.0
	if(xxx.gt.0.0)gamay = rr23/xxx

      ydot(1)  = rr8 + rr33				!O3
     &	       -(xbeta*rr1 + rr10 + rr11 + rr12 + rr14)		

      ydot(2)  = rr26 + rr27				!CO
     &         - rr5

      ydot(3)  = rr5 					!CO2

      ydot(4)  = 2.0*xbeta*rr1 + rr7 + rr10 		!HO
     &	       + 2.0*rr17+ rr24
     &         -(rr5 + rr11 + rr13 + rr18 + rr19 
     &          +rr25+ rr27 + rr29
     &	        +rr32+ rr33 + rr34 )

      ydot(5)  = rr5 + rr11 + rr18 + rr27 + rr29	!HO2
     &         + rr24+ 2.0*rr26
     &         +(2.0*gamax - 1.0)*(rr19 + rr25)
     &         -(rr7 + rr10 + rr16 + rr32)

      ydot(6)  = rr16 + rr34				!H2O2
     &         -(rr17+ rr18)

      ydot(7)  = rr8 					!NO
     &	       -(rr7 + rr12 + gamax*(rr19 + rr25))

      ydot(8)  = rr7 + rr12 + gamax*(rr19 + rr25)	!NO2
     &	       -(rr8 + rr13+ rr14+ rr14)

      ydot(9)  =-rr19					!CH4

      ydot(10) = rr24 + gamax*(rr19 + rr25)		!CH2O
     &	       -(rr26 + rr27)

      ydot(11) = gamay*(rr19 + rr25)			!CH3O2H
     &	       -(rr24 + rr25) 

      ydot(12) =-rr29 					!SO2

      return
       end

c
      subroutine jex (neq, t, y, ml, mu, pd, nrpd)
c     ==========================================
c
c  Subroutine for calculating the Jacobian for f(i)
c
c     P(m,n) here represents df(m)/dn
c                         or d ydot(m)/dn
c

#include "chem_para"
#include "chem_com"

      common /c00tmp/c00(2)

      dimension y(neq), pd(nrpd,neq)

      bbb = 0.7809*c00(1)
      ccc = 0.2095*c00(1)

        xbeta1 = rk(3)*bbb + rk(4)*ccc   !k3[N2] + k4[O2]
        xbeta2 = rk(2)*c00(2)            !k2[H2O]
        xbeta  = xbeta2/(xbeta1 + xbeta2)

        rr21= rk(21)*y(7)               !k21[NO] - for ss only
        rr23= rk(23)*y(5)               !k23[HO2]- for ss only
        xxx   = (rr21 + rr23)

        gamax = 0.0
        if(rr21.gt.1.e-12)gamax = rr21/xxx

        gamay = 0.0
        if(xxx.gt.0.0)gamay = rr23/xxx

      do 1 j=1,neq
      do 1 i=1,neq
        pd(i,j) = 0.0
1     continue

c =====
c ===== pd(i,j) = d ydot(i)/dy(j):
c =====
	
	rk5y2    = rk(5)*y(2)
	rk5y4    = rk(5)*y(4)
	rk7y5    = rk(7)*y(5)
	rk7y7    = rk(7)*y(7)
	rk10y1   = rk(10)*y(1)
	rk10y5   = rk(10)*y(5)
	rk11y1   = rk(11)*y(1)
	rk11y4   = rk(11)*y(4)
	rk12y1   = rk(12)*y(1)
	rk12y7   = rk(12)*y(7)
	rk13y4   = rk(13)*y(4)
	rk13y8   = rk(13)*y(8)
	rk14y1   = rk(14)*y(1)
	rk14y8   = rk(14)*y(8)
	rk18y4   = rk(18)*y(4)
	rk18y6   = rk(18)*y(6)
	rk19y4   = rk(19)*y(4)
	rk19y9   = rk(19)*y(9)
	rk25y4   = rk(25)*y(4)
	rk25y11  = rk(25)*y(11)
	rk27y4   = rk(27)*y(4)
	rk27y10  = rk(27)*y(10)
	rk29y4   = rk(29)*y(4)
	rk29y12  = rk(29)*y(12)
	rk32y4   = rk(32)*y(4)
	rk32y5   = rk(32)*y(5)

	pd(1,1)  = -(xbeta*rk(1)   + rk10y5 + rk11y4
     &	                  +rk12y7  + rk14y8)
	pd(4,1)  = 2.0*xbeta*rk(1) + rk10y5 - rk11y4
	pd(5,1)  = rk11y4 - rk10y5
	pd(7,1)  = rk(12)*y(7)				! 012507 added
	pd(8,1)  = rk12y7 - 2.*rk14y8

	pd(2,2)  =-rk5y4
	pd(3,2)  = rk5y4
	pd(4,2)  =-rk5y4
	pd(5,2)  = rk5y4

	pd(1,4)  = 2.0*rk(33)*y(4) - rk11y1
	pd(2,4)  = rk27y10    - rk5y2
	pd(3,4)  = rk5y2
	pd(4,4)  =-( 
     &	 	   rk5y2  + rk11y1 + rk(13)*y(8)
     &		 + rk18y6 + rk19y9 + rk25y11
     &		 + rk27y10+ rk29y12+ rk32y5
     &		 +2.0*(rk(33) + rk(34))*y(4)		! 012507 added 2.0
     &              )
	pd(5,4)  =(rk5y2 + rk11y1 + rk18y6
     &           + rk27y10+ rk29y12)
     &           +(2.0*gamax - 1.0)*(rk19y9 + rk25y11)
     &		 - rk32y5
	pd(6,4)  = 2.0*rk(34)*y(4) -rk18y6		! 012507 fixed
	pd(7,4)  =-gamax*(rk19y9 + rk25y11)
	pd(8,4)  = gamax*(rk19y9 + rk25y11)
     &		 - rk13y8
	pd(9,4)  =-rk19y9
	pd(10,4) = gamax*(rk19y9 + rk25y11)
     &		 - rk27y10
	pd(11,4) = gamay*(rk19y9 + rk25y11)
     &		 - rk25y11
	pd(12,4) =-rk29y12

	pd(1,5)  = -rk10y1	! fixed 012507 originally rk11y1
	pd(4,5)  = rk7y7 + rk10y1
     &           - rk32y4
	pd(5,5)  =
     &		 -(rk7y7 + rk10y1 + rk(16)*y(5) + rk32y4)
	pd(6,5)  = 2.0*rk(16)*y(5)			! 012507 fixed
	pd(7,5)  =-rk7y7
	pd(8,5)  = rk7y7
c 	pd(11,5)

	pd(4,6)  = 2.0*rk(17)
     &           - rk18y4
	pd(5,6)  = rk18y4
	pd(6,6)  =-(rk(17) + rk18y4)

	pd(1,7)  =-rk12y1
	pd(4,7)  = rk7y5
	pd(5,7)  =-rk7y5
	pd(7,7)  =-(rk7y5 + rk12y1)
	pd(8,7)  = (rk7y5 + rk12y1)
	
	pd(1,8)  = rk(8)
     &		 - rk14y1
	pd(4,8)  =-rk13y4
	pd(7,8)  = rk(8)
	pd(8,8)  =-(rk(8) + rk13y4 + 2.0*rk14y1)

	pd(4,9)  =-rk19y4
	pd(5,9)  =(2.0*gamax - 1.0)*rk19y4
	pd(7,9)  =-gamax*rk19y4
	pd(8,9)  = gamax*rk19y4
	pd(9,9)  =-rk19y4
	pd(10,9) = gamax*rk19y4
	pd(11,9) = gamay*rk19y4

	pd(2,10) = rk(26) + rk27y4
	pd(4,10) =-rk27y4
	pd(5,10) = rk27y4 + 2.0*rk(26)
	pd(10,10)=-(rk(26)+ rk27y4)

	pd(4,11) = rk(24) - rk25y4
	pd(5,11) = rk(24) + (2.0*gamax - 1.0)*rk25y4
	pd(7,11) =-gamax*rk25y4
	pd(8,11) = gamax*rk25y4
	pd(10,11)= rk(24) + gamax*rk25y4
	pd(11,11)=-(rk(24) + rk25y4*max(0.0,1.0 - gamay))

	pd(4,12) =-rk(29)*y(4)
	pd(5,12) = rk29y4
	pd(12,12)=-rk29y4

      return
      end
 

c==================================================================c
c                                                                  c
c   CHEMFUNCPACK.F:  Subroutine and function pack                  c
c                             of MIT Global Chemistry Model        c
c        -------------------------------------------------         c
c   Author:		Chien Wang                                 c
c			MIT Joint Program on Science and Policy    c
c			    of Global Change                       c             
c   Last Revised on:    February 22, 1995                          c
c                                                                  c
c==================================================================c

c
      subroutine ppbm2mcm (np, xx, den, revaw) 
c     ==========================================

c------------------------------------------------
c   A program for convert mixing ratio in
c     ppb(m) to concentration in molecules/cm^3
c------------------------------------------------

      dimension xx    (np)
      dimension revaw (np)

      ddd = 6.02217e+11*den
      do 1 i=1,np
        xx (i) = xx (i)
     &          *ddd*revaw(i)
1     continue

      return
       end

c
      subroutine mcm2ppbm (np, xx, den, aw)
c     ==========================================

c------------------------------------------------
c   A program for convert concentration in
c     molecules/cm^3 to mixing ratio in ppb(m)
c------------------------------------------------

      dimension xx (np)
      dimension aw (np)

      ddd = 1./(6.02217e+11*den)
      do 1 i=1,np
        xx (i) = xx (i)
     &          *ddd*aw(i)
1     continue

      return
       end

c
      subroutine rateph (T, ntemp, Rflux, coszagcm)
c     ===========================

c------------------------------------------------
c   A program for calculating rate constants
c     of photochemical reactions
c         ----------------------------
c   The time unit for all reaction rates is second
c
c       New version revised at June 3, 1995
c
c------------------------------------------------

#include "chem_para"
#include "chem_com"
     
	coszagcm = max(cosza4rk(10),
     &             min(cosza4rk(1),coszagcm))

	do 1 nband=1,9
	  if(coszagcm.le.cosza4rk(nband).and.
     &       coszagcm.gt.cosza4rk(nband+1))then

	    wgtxx = 1.0/(cosza4rk(nband) - cosza4rk(nband+1))
	    wgt1  = (cosza4rk(nband) - coszagcm)  *wgtxx
	    wgt2  = (coszagcm - cosza4rk(nband+1))*wgtxx
c
c	+--- wgt2 ---+--- wgt1 ---+
c	nband --- coszagcm --- nband+1

c-- R8: NO2 + hv -> NO + O
c   Atkinson et al., 1992:
c
c        rk(8) = (18096.7314929
c     &          -3.9704*(T-273.))
c     &        *5.0e-10
c     &        *Rflux

c        rk(8) = rktable1(08,ntemp)*Rflux

        rk(8) = ((rk08gama(nband)*wgt2 + rk08gama(nband+1)*wgt1)
     &        +  (rk08aaa (nband)*wgt2 + rk08aaa (nband+1)*wgt1)
     &          *(T-273.))
     &        *Rflux

c-- R1: O3 + hv -> O(1D) + O2
c   Atkinson et al., 1992:
c
c        rk(1)  = 0.0028*rk(8)
cc       rk(1) = 3.467625419939578E-08*Rflux

c	rk(1) = (rk01table(nband)*wgt2 
c     &         + rk01table(nband+1)*wgt1)*Rflux
        rk(1)  = 0.0028*rk(8)

c-- R17: H2O2 + hv -> 2OH
c   Atkinson et al., 1992:
c
c       rk(17)= 6.263464249748237E-09*Rflux

	rk(17) = (rk17table(nband)*wgt2 
     &         +  rk17table(nband+1)*wgt1)*Rflux

c-- R24: CH3O2H + hv -> CH3O + HO
c   Atkinson et al., 1992:
c
c       rk(24)= 4.474208962739173E-09*Rflux

	rk(24) = (rk24table(nband)*wgt2 
     &         +  rk24table(nband+1)*wgt1)*Rflux

c-- R26: CH2O + hv -> CHO + H
c   Atkinson et al., 1992:
c
c       rk(26)= 9.410107812688823E-08*Rflux

	rk(26) = (rk26table(nband)*wgt2 
     &         +  rk26table(nband+1)*wgt1)*Rflux

	goto 2
	endif

1	continue

c = For znith angle.ge.86:

        rk(8) = (rk08gama(10)
     &        +  rk08aaa (10)
     &          *(T-273.))
     &        *Rflux

c	rk(1) = rk01table(10)
c     &        * Rflux
        rk(1)  = 0.0028*rk(8)

	rk(17) = rk17table(10)
     &         * Rflux

	rk(24) = rk24table(10)
     &         * Rflux

	rk(26) = rk26table(10)
     &         *Rflux

2	continue

      return
       end

c
      subroutine ratebm (T, ntemp, patm)
c     ===========================

c------------------------------------------------
c   A program for calculating rate constants
c     of bimolecular reactions
c         ----------------------------
c   The time unit for all reaction rates is second
c
c------------------------------------------------

#include "chem_para"
#include "chem_com"

	Trev = 1./T

c-- R2: O(1D) + H2O -> 2OH

	rk(2)  = 2.2e-10

c-- R3: O(1D) + N2 -> O + N2

c	rk(3)  = 1.8e-11 * exp( 107.0*Trev)
	rk(3)  = rktable1(03,ntemp)

c-- R4: O(1D) + O2 -> O + O2

c	rk(4)  = 3.2e-11 * exp(  67.0*Trev)
	rk(4)  = rktable1(04,ntemp)

c-- R5: CO + HO -> H + CO2

	rk(5)  = 1.5e-13 * (1.0 + 0.6*patm)

c-- R7: HO2 + NO -> HO + NO2

c	rk(7)  = 3.7e-12 * exp( 240.0*Trev)
	rk(7)  = rktable1(07,ntemp)

c-- R10: HO2 + O3 -> HO + 2O2

c	rk(10) = 1.1e-14 * exp(-500.0*Trev)
	rk(10) = rktable1(10,ntemp)

c-- R11: HO + O3 -> HO2 + O2

c	rk(11) = 1.6e-12 * exp(-940.0*Trev)
	rk(11) = rktable1(11,ntemp)
c	rk(11) = 0.0

c-- R12: NO + O3 -> NO2 + O2

c	rk(12) = 2.0e-12 * exp(-1400.0*Trev)
	rk(12) = rktable1(12,ntemp)

c-- R14: NO2 + O3 -> NO3 + O2

c	rk(14) = 1.2e-13 * exp(-2450.0*Trev)
	rk(14) = rktable1(14,ntemp)

c-- R16: HO2 + HO2 -> H2O2 + O2

c	rk(16) = 2.3e-13 * exp( 600.0*Trev)

c-- R18: H2O2 + HO -> HO2 + H2O

c	rk(18) = 2.9e-12 * exp(-160.0*Trev)
	rk(18) = rktable1(18,ntemp)

c-- R19: CH4 + HO -> CH3 + H2O

c	rk(19) = 2.65e-12 * exp(-1800.0*Trev)
	rk(19) = rktable1(19,ntemp)

c-- R21: CH3O2 + NO -> CH3O + NO2

c	rk(21) = 4.2e-12 * exp( 180.0*Trev)
	rk(21) = rktable1(21,ntemp)

c-- R22: CH3O + O2 -> CH2O + HO2

c	rk(22) = 3.9e-14 * exp(-900.0*Trev)
	rk(22) = rktable1(22,ntemp)

c-- R23: CH3O2 + HO2 -> CH3O2H + O2

c	rk(23) = 3.8e-13 * exp( 780.0*Trev)
	rk(23) = rktable1(23,ntemp)

c-- R25: CH3O2H + HO -> CH3O2 + H2O

c	rk(25) = 1.9e-12 * exp( 190.0*Trev)
	rk(25) = rktable1(25,ntemp)

c-- R27: CH2O + HO -> CHO + H2O

	rk(27) = 1.0e-11

c-- R28: CHO + O2 -> CO + HO2

c	rk(28) = 3.5e-12 * exp( 140.0*Trev)
	rk(28) = rktable1(28,ntemp)

c-- R30: HOSO2 + O2 -> HO2 + SO3

c	rk(30) = 1.3e-12 * exp(-330.0*Trev)
	rk(30) = rktable1(30,ntemp)

c-- R31: SO3 + H2O -> H2SO4

	rk(31) = 2.4e-15

c New reactions - 062995:

c-- R32: HO + HO2 -> H2O + O2

c	rk(32) = 4.8e-11 * exp( 250.0*Trev)
	rk(32) = rktable1(32,ntemp)

c-- R33: HO + HO -> H2O + O

c	rk(33) = 4.2e-12 * exp(-240.0*Trev)
	rk(33) = rktable1(33,ntemp)

      return
       end

c
      subroutine ratetm (T, ntemp, ym)
c     =========================

c------------------------------------------------
c   A program for calculating rate constants
c     of termolecular reactions
c         ----------------------------
c   The time unit for all reaction rates is second
c
c------------------------------------------------

#include "chem_para"
#include "chem_com"

        dimension ym(2)

        T300 = T/300.
        xm   = ym(1)
        xn2  = 0.7809*xm
        h2o  = ym(2)

c-- R6:  H + O2 + M -> HO2 + M

c	rkt0   = 6.2e-32 * xm * T300**(-1.6)
	rkt0   = rktable1(06,ntemp) * xm
	rkt00  = rktable2(01,ntemp)
	rk(6)  = calkmt(rkt0, rkt00)

c-- R9:  O + O2 + M -> O3 + M  (no high pres. limit)

c	rkt0   = 5.6e-34 * xm * T300**(-2.8)
	rkt0   = rktable1(09,ntemp) * xm
	rk(9)  = rkt0

c-- R13: NO2 + HO + M -> HNO3 + M

c	rkt0   = 2.6e-30 * xm * T300**(-3.2)
c	rkt00  = 2.4e-11 * T300**(-1.3)
	rkt0   = rktable1(13,ntemp) * xm
	rkt00  = rktable2(02,ntemp)
	rk(13) = calkmt(rkt0, rkt00)

c-- R15: NO3 + NO2 + M -> N2O5 + M

c	rkt0   = 2.2e-30 * xm * T300**(-3.9)
c	rkt00  = 1.5e-12 * T300**(-0.7)
	rkt0   = rktable1(15,ntemp) * xm
	rkt00  = rktable2(03,ntemp)
	rk(15) = calkmt(rkt0, rkt00)

c-- R16: HO2 + HO2 + M -> H2O2 + O2
c  062895:

c	rk(16) = 1.7e-33 * xn2 * exp(1000.0/T)
	rk(16) = rktable1(16,ntemp) * xn2
     &         *(1.+ 1.4e-21*h2o*exp(2200.0/T))

c-- R20: CH3 + O2 + M -> CH3O2 + M

c	rkt0   = 4.5e-31 * xm * T300**(-3.0)
c	rkt00  = 1.8e-12 * T300**(-1.7)
	rkt0   = rktable1(20,ntemp) * xm
	rkt00  = rktable2(04,ntemp)
	rk(20) = calkmt(rkt0, rkt00)

c-- R29: SO2 + HO + M -> HOSO2 + M

c	rkt0   = 3.0e-31 * xm * T300**(-3.3)
	rkt0   = rktable1(29,ntemp) * xm
	rkt00  = rktable2(05,ntemp)
	rk(29) = calkmt(rkt0, rkt00)

c New reactions, 063095:

c-- R34: HO + HO + M -> H2O2 + M

c        rkt0   = 6.9e-31 * xm * T300**(-0.8)
        rkt0   = rktable1(34,ntemp) * xm
        rkt00  = rktable2(06,ntemp)
        rk(34) = calkmt(rkt0, rkt00)


      return
       end

c
      function calkmt (rkt0, rkt00)
c     =============================

c
c  A function for calculate k(M,T) of termolecular
c    reaction rate from low and high limit
c

      aaa = rkt0/rkt00
      bbb = dlog10(aaa)**2
      ccc = 1./(1.+bbb)

      calkmt = rkt0/(1.+aaa) * 0.6**ccc

      return
       end


c
	Block data rktable3dat
c	======================

c----------------------------------------c
c  Block data for holding rktable3 data  c
c----------------------------------------c

#include "chem_para"
#include "chem_com"

	data cosza4rk/
     &      0.100000E+01,  0.984808E+00,  0.939693E+00,  0.866025E+00,
     &      0.766044E+00,  0.642788E+00,  0.500000E+00,  0.342020E+00,
     &      0.207912E+00,  0.697565E-01/

	data rk08gama/
     &      0.958260E-05,  0.956540E-05,  0.947786E-05,  0.933214E-05,
     &      0.904837E-05,  0.853986E-05,  0.790570E-05,  0.644646E-05,
     &      0.540693E-05,  0.640000E-05/

	data rk08aaa/
     &     -0.213357E-08, -0.212833E-08, -0.210322E-08, -0.206239E-08,
     &     -0.198520E-08, -0.185147E-08, -0.169861E-08, -0.133796E-08,
     &     -0.109912E-08, -0.131294E-08/

	data rk01table/
     &      0.553973E-07,  0.541612E-07,  0.501300E-07,  0.436050E-07,
     &      0.346763E-07,  0.241176E-07,  0.141663E-07,  0.593145E-08,
     &      0.238321E-08,  0.148548E-08/

	data rk17table/
     &      0.761430E-08,  0.755268E-08,  0.730677E-08,  0.691454E-08,
     &      0.626346E-08,  0.536943E-08,  0.430112E-08,  0.288139E-08,
     &      0.199044E-08,  0.209279E-08/

	data rk24table/
     &      0.535338E-08,  0.531281E-08,  0.515418E-08,  0.489785E-08,
     &      0.447421E-08,  0.387878E-08,  0.318768E-08,  0.217960E-08,
     &      0.155503E-08,  0.170283E-08/

	data rk26table/
     &      0.109106E-06,  0.108473E-06,  0.105771E-06,  0.101587E-06,
     &      0.941011E-07,  0.831775E-07,  0.701839E-07,  0.494384E-07,
     &      0.364853E-07,  0.410051E-07/

	end

#endif
#endif
