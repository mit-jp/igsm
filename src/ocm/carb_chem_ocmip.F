
#include "ctrparam.h"

!	==========================================================
!
!	CARB_CHEM_OCMIP.F:  Subroutine for computing pCO2 in
!                           mixed layer ocean for the IGSM2.
!                           carbon coefficientsnd CO2* found using
!                           OCMIP2 protocols (which follow:
!                         DOE Methods Handbook 1994 Ver.2, 
C ORNL/CDIAC-74, dickson and Goyet, eds. (Ch 2 p 10, Eq A.49) and
!                           Millero).
!                           pH found using Follows et al (2006)
!
! NOTE ON UNITS: IGSM2.2 wants units of mol/l, calc_pco2_follows.F
!                wants mol/m3 (as input - has its own conversion
!                              to mol/kg internally)
!	     -------------------------------------------------
! using structure of carb_mxdlyr_chem.F (by G. Holian)
!
!	Author: 	    S. Dutkiewicz
!                           Oct 2006
!	----------------------------------------------------------
! 
!	Revision History:
!	
!	When	Who		What
!	-----	----------	-------
!	==========================================================

c
      subroutine carb_chem_ocmip(focean)
c     ===========================

!     implicit none
#include "ODIFF.COM"
#include "com_carb"

      REAL surfphos(nlat), surfsi(nlat), AlkT(nlat), S(nlat)

      real focean(nlat)
      common /SSS/Salt(nlat)

      data ifirst /1/
!     save

      real   ff(nlat)    !fugactity
      real   ak0(nlat)   !solubility
      real   ak1(nlat), ak2(nlat)   !carbonate eq. coeff
      real   akb(nlat), aksi(nlat)  !borate and silicate eq. coeff
      real   ak1p(nlat), ak2p(nlat), ak3p(nlat)  ! phos eq. coeff
      real   akw(nlat)   !dissociation in water
      real   bt(nlat)    !diss org. borate
      
      real  dic(nlat), co2s(nlat), co3(nlat) !in units mol/m3

      real  sss
      real  ta, pt, sit
      real  tk, tk100, tk1002, dlogtk, invtk
      real  sqrtis, sqrts, s15, scl, s2, is, is2

      real  convert    !mol/m3 to mol/l

      integer it, nit
!	-----------------------------------------------------------

      convert=1000.     !liters to m3

      if (ifirst.eq.1) then

        read (668),AlkT         !mol eq/m3
        print *,'AlkT',iyear
        print '12F7.2',(AlkT(j),j=1,nlat)

        read (669),S

!       read (670),DICEQ        !unit??

        DO j=1,nlat
          surfphos(j)= 5.1225e-4  !mol/m3
          surfsi(j)= 7.6838e-3    !mol/m3
          Salt(j)=S(j)
          volume_ml(j) = area_ml(j)*depth_ml(j)
!         DICEQ(j)=(DICEQ(j)*volume_ml(j))    !mol/m3 -> mol   
!         TCO2_ml(j)=(DICEQ(j)/volume_ml(j))/convert  !mol/m3 -> mol/l
          TCO2_ml(j)=RCO2(j,1)*1.e-3
          co2s(j)=0.e0
          co3(j)=0.e0
        ENDDO
        print *,' Initial Hguess'
        print *,Hguess
        print *,' Initial TCO2_ml'
        print *,TCO2_ml
      endif	! ifirst


c**** Compute Temperature Dependent Equilibrium Constants 
c*** note that units will be in mol/kg

      DO j=1,nlat
       if(focean(j).gt.0.0) then
        tk=TEMP(j)  ! in kelvin
        tk100 = tk/100.0
        tk1002=tk100*tk100
        invtk=1.0/tk
        dlogtk=log(tk)
        sss=S(j)
        is=19.924*sss/(1000.-1.005*sss)
        is2=is*is
        sqrtis=sqrt(is)
        s2=sss*sss
        sqrts=sqrt(sss)
        s15=sss**1.5
        scl=sss/1.80655

C------------------------------------------------------------------------
C f = k0(1-pH2O)*correction term for non-ideality (fugacity)
C Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)
        ff(j)=exp(-162.8301 + 218.2968/tk100  +
     &          90.9241*log(tk100) - 1.47696*tk1002 +
     &          sss * (.025695 - .025225*tk100 +
     &          0.0049867*tk1002))
c ------------------------------------------------------------------------
C K0 from Weiss 1974
           ak0(j) = exp(93.4517/tk100 - 60.2409 +
     &        23.3585 * log(tk100) +
     &        sss * (0.023517 - 0.023656*tk100 +
     &        0.0047036*tk1002))
C------------------------------------------------------------------------
C k1 = [H][HCO3]/[H2CO3]
C k2 = [H][CO3]/[HCO3]
C Millero p.664 (1995) using Mehrbach et al. data on seawater scale 
           ak1(j)=10**(-1*(3670.7*invtk -
     &          62.008 + 9.7944*dlogtk -
     &          0.0118 * sss + 0.000116*s2))
           ak2(j)=10**(-1*(1394.7*invtk + 4.777 -
     &          0.0184*sss + 0.000118*s2))
C------------------------------------------------------------------------
C kb = [H][BO2]/[HBO2]
C Millero p.669 (1995) using data from dickson (1990)
           akb(j)=exp((-8966.90 - 2890.53*sqrts - 77.942*sss +
     &          1.728*s15 - 0.0996*s2)*invtk +
     &          (148.0248 + 137.1942*sqrts + 1.62142*sss) +
     &          (-24.4344 - 25.085*sqrts - 0.2474*sss) *
     &          dlogtk + 0.053105*sqrts*tk)
C------------------------------------------------------------------------
C k1p = [H][H2PO4]/[H3PO4]
C DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
           ak1p(j) = exp(-4576.752*invtk + 115.525 -
     &          18.453*dlogtk + 
     &          (-106.736*invtk + 0.69171)*sqrts +
     &          (-0.65643*invtk - 0.01844)*sss)
C------------------------------------------------------------------------
C k2p = [H][HPO4]/[H2PO4]
C DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
           ak2p(j) = exp(-8814.715*invtk + 172.0883 -
     &          27.927*dlogtk + 
     &          (-160.340*invtk + 1.3566) * sqrts +
     &          (0.37335*invtk - 0.05778) * sss)
C------------------------------------------------------------------------
C k3p = [H][PO4]/[HPO4]
!DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
           ak3p(j) = exp(-3070.75*invtk - 18.141 +
     &          (17.27039*invtk + 2.81197) *
     &          sqrts + (-44.99486*invtk - 0.09984) * sss)
C------------------------------------------------------------------------
C ksi = [H][SiO(OH)3]/[Si(OH)4]
C Millero p.671 (1995) using data from Yao and Millero (1995)
           aksi(j) = exp(-8904.2*invtk + 117.385 -
     &          19.334*dlogtk +
     &          (-458.79*invtk + 3.5913) * sqrtis +
     &          (188.74*invtk - 1.5998) * is +
     &          (-12.1652*invtk + 0.07871) * is2 +
     &          log(1.0-0.001005*sss))
C------------------------------------------------------------------------
C kw = [H][OH]
C Millero p.670 (1995) using composite data
           akw(j) = exp(-13847.26*invtk + 148.9652 -
     &          23.6521*dlogtk +
     &          (118.67*invtk - 5.977 + 1.0495 * dlogtk) *
     &          sqrts - 0.01615 * sss)
C------------------------------------------------------------------------
C Calculate concentrations for borate, sulfate, and fluoride
C Uppstrom (1974)
           bt(j) = 0.000232 * scl/10.811
C------------------------------------------------------------------------
         else
            ff(j)=0.0
            ak0(j)= 0.0
            ak1(j)= 0.0
            ak2(j)= 0.0
            akb(j)= 0.0
            ak1p(j) = 0.0
            ak2p(j) = 0.0
            ak3p(j) = 0.0
            aksi(j) = 0.0
            akw(j) = 0.0
            bt(j) = 0.0
         endif
c----------------------------------------------------------
      ENDDO

        if (ifirst.eq.1) then
             nit=10    !if start of run, do multiple iterations
             ifirst=0
        else
            nit=1
        endif

            nit=10
!           nit=1
!     print *,'nit=',nit
      DO j=1,nlat

       if (focean(j).gt.0.0) then
        do it=1,nit
!          print *,j,it
!           TCO2_ml(j)=RCO2(j,1)*1.e-3
            pH(j)=-log10(Hguess(j))
            dic(j)=TCO2_ml(j)*convert               !mol/l -> mol/m3
            CALL CALC_PCO2_FOLLOWS( 
     I        dic(j), surfphos(j),
     I        surfsi(j),AlkT(j),
     I        ak1(j),ak2(j),
     I        ak1p(j),ak2p(j),ak3p(j),
     I        akb(j),akw(j),
     I        aksi(j),ff(j),bt(j), 
     U        pH(j),co2s(j),co3(j) )
            Hguess(j)=10.0**(-pH(j))
        enddo 
!        print *,j,
!    I        dic(j), surfphos(j),
!    I        surfsi(j),AlkT(j),
!    I        ak1(j),ak2(j),
!    I        ak1p(j),ak2p(j),ak3p(j),
!    I        akb(j),akw(j),
!    I        aksi(j),ff(j),bt(j), 
!    U        pH(j),co2s(j),co3(j)
c not completely sure about the next line (but it isn't needed anyway)
        CO2_ml(j) =co2s(j)/convert             ! mol/m3 -> mol/l
        CO3_ml(j) =co3(j)/convert              ! mol/m3 -> mol/l
        HCO3_ml(j) = TCO2_ml(j)-CO3_ml(j)-CO2_ml(j)
        pC_mla(j) = (CO2_ml(j)/ff(j))       ! QQ units atm?? 
        H(j)=Hguess(j)
       else
        pH(j)=0.0
        dic(j)=0.0
        CO2_ml(j)=0.0
        CO3_ml(j) = 0.0
        HCO3_ml(j) = 0.0
        pC_mla(j) = 0.0
        Hguess(j)=0.0
        H(j)=0.0

       endif !focean

!      print *,j,ff(j),dic(j),CO2_ml(j)
!      print *,pC_mla(j)
       if(j.eq.-27)then
        print *,'MLCHEM'
        print *,CO2_ml(j),ff(j)
        print *,pH(j)
        print *,'temp=',TEMP(j)
       endif

      ENDDO
!     print *,'Before return'
      return
      end

