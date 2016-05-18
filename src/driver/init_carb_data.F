
#include "ctrparam.h"

c
      subroutine init_carb_data(co2eq)
c     =========================

c==================================================================c
c                                                                  c
c   INIT_CARB_DATA.F:   Subroutine for computing the Carbon        c
c                      Dioxide Data for the Ocean in the MIT       c
c                          Global Climate/Chemistry Model          c
c        -------------------------------------------------         c
c   Author:             Gary L. Holian                             c
c                       MIT Joint Program on Science and Policy    c
c                           of Global Change                       c
c   Last Revised on:    July 15, 1996                              c
c                                                                  c
c==================================================================c
 

!#include "chem_para"
#include "ODIFF.COM"
#include "com_carb"


c**** Compute Data Vectors

      

      DO 10 j=1,nlat

c     pC_atm(j) = 278.000000E-6
      pC_atm(j) = co2eq*1.e-6
      volume_ml(j) = area_ml(j)*depth_ml(j)

 10   continue


      DO 15 j=1,nlat

      pC_mla(j) = 275.000000E-6
      pC_mla(j) = (co2eq-3.)*1.e-6
      TCO2_ml(j)=pC_mla(j)*3.2*2.

 15   continue

!     read (667)DIC_ml



      DO 100 j=1,nlat
      
      if(volume_ml(j).gt.0.0)then
!       DIC_ml(j)=DIC_ml(j)*volume_ml(j)
        DIC_ml(j)=TCO2_ml(j)*volume_ml(j)*1000.
      else
        DIC_ml(j)=0.0
      end if

      DO 90 k=1,LMO
 
      RCO2(j,k) = 0.0
      if(volume_ml(j).gt.0.0)then
        RCO2(j,k) = DIC_ml(j)/volume_ml(j)
      end if

 90   continue

 100  continue
	print *,' from INIT'
cprint *,' area_srf'
cprint*,area_srf
	print *,' depth_ml'
	print *,depth_ml
cprint *,' RCO2'
cprint *,RCO2
cprint *,' pC_mla'
cprint *,pC_mla
	print *,' DIC'
	print *,DIC_ml
cprint *,' from INIT'
cprint *,' ************'
cprint *,'wind_amp=',wind_amp
cprint *,' ************'
c       print *,'pC_atm*10**4'
c       print 5001,(pC_atm(j)*1.e6-300.,j=1,nlat)
c       print *,'pC_mla*10**4'
c       print 5001,(pC_mla(j)*1.e6-300.,j=1,nlat)
c       print *,'DIC_ml*10**(-15)'
c       print 5001,(DIC_ml(j)*1.e-15,j=1,nlat)
c       print *,'TCO2_ml'
c       print 5001,(TCO2_ml(j)*1.e-6,j=1,nlat)
c       print *,' RCO2'
c       print 5001,((RCO2(j,k),j=1,nlat),k=1,LMO)
 5001 format(24f5.1)


      return

      end
