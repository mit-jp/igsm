
#include "ctrparam.h"

!	==========================================================
!
!	CARB_AIRSEA_FLX.F:  Subroutine for computing the Air to
!			   Sea Flux of Carbon Dioxide in the MIT
!			       Global Climate/Chemistry Model
!	     -------------------------------------------------
!
!	Author: 	    Gary L. Holian
!			    MIT Joint Program on Science and Policy
!				of Global Change
!
!	----------------------------------------------------------
! 
!	Revision History:
!	
!	When	Who		What
!	-----	----------	-------
!	102996	G. Holian	rev.
!	080200	Chien Wang	repack based on CliChem3 and add cpp.
!
!	==========================================================

c
      subroutine carb_airsea_flx(dtco2)
c     ==========================

!#include "chem_para"
#include "ODIFF.COM"
#include "com_carb"

      common /SSS/Salt(nlat)
      real pC_atml(nlat)

!	----------------------------------------------------------
      sca1 = 2073.1
      sca2 = -125.62
      sca3 =    3.6276
      sca4 =   -0.043219


C**** Change in DIC_ml Following Diffusion

      DO 800 j=1,nlat
c
      volume_ml(j) = area_ml(j)*depth_ml(j)

      DIC_ml(j) = RCO2(j,1)*volume_ml(j)
!     DIC_ml(j) = DICEQ(j)+RCO2(j,1)*volume_ml(j)

 800  continue


c**** Compute the Piston velocity for the Air to Sea Exchange of CO2
c**** as a function of the current surface wind speeds. Units of Vpiston
c**** are moles/sec*microatm*m^2. 
   
      
      DO 10 j=1,nlat

         tempc=TEMP(j)-273.15
         pisvel = 0.337*wind(j)**2/3.6d5
c****     m/s
     &    *dtco2
c****   convertion  to m/xxx

!calculate SCHMIDT NO. for CO2
              SchmidtNoDIC =
     &            sca1
     &          + sca2 * tempc
     &          + sca3 * tempc*tempc
     &          + sca4 * tempc*tempc
     &                *tempc

         Vpiston(j) = pisvel/
     &    sqrt(SchmidtNoDIC/660.0)
c**** convertion moles/sec*microatm*m^2. 


 10   continue


c**** Compute the CO2 concentration differences between the Atmosphere
c**** and the Surface Ocean Mixed Layer. 
!         Units are microatm (10^-6 atm) ?


      DO 20 j=1,nlat

      tk100=TEMP(j)/100.
      tk1002=tk100*tk100

      ak0 = exp(93.4517/tk100 - 60.2409 +
     &        23.3585 * log(tk100) +
     &        salt(j) * (0.023517 - 0.023656*tk100 +
     &        0.0047036*tk1002))

           ff = exp(-162.8301 + 218.2968/tk100  +
     &          90.9241*log(tk100) - 1.47696*tk1002 +
     &          salt(j) * (.025695 - .025225*tk100 +
     &          0.0049867*tk1002))
      
       CO2_atm=pC_atm(j)*ak0
!      pC_mla(j)=CO2_ml(j)/ff

!      if(j.eq.23) then
!        print *,'OCM'
!        print *,'ATM'
!        print *,pC_atm(j),CO2_atm
!        print *,'OCEAN'
!        print *,CO2_ml(j),pC_mla(j),CO2_ml(j)/ak0
!      endif

       dpCO2_atm_ml(j) = (CO2_atm - CO2_ml(j))

!      dpCO2_atm_ml(j) = (CO2_atm - CO2_ml(j)/ff*ak0)

!      ddp=(1E6)*(pC_atm(j) - pC_mla(j))

!      if(j.eq.23) then
!        print *,dpCO2_atm_ml(j),Vpiston(j)
!        print *,ddp,Vpiston(j)*ak0
!      endif

 20   continue


c**** Compute the CO2 flux to the Surface Ocean Mixed Layer arising from
c**** a concentration difference between air and sea. Units are moles/time.


      DO 30 j=1,nlat

      if(area_srf(j).gt.0.)then
      Flux_atm_ml(j) = Vpiston(j)*area_srf(j)*dpCO2_atm_ml(j)
     &  *1024.0
!      density
!      if(j.eq.23)then
!        print *,pC_atm(j)*1.e6,pC_mla(j)*1.e6
!        print *,Vpiston(j),area_srf(j),wind(j)
!        print *,Flux_atm_ml(j),Vpiston(j)*ak0**area_srf(j)*ddp
!     endif
      end if

 30   continue


c**** Add the CO2 flux to the existing pool of Dissolved Inorganic Carbon 
c**** in the Surface Ocean Mixed Layer. Units of DICml are moles.


      DO 40 j=1,nlat

      if(volume_ml(j).gt.0.0)then

      DIC_ml(j) = DIC_ml(j) + Flux_atm_ml(j)

	fluxco2(j) = Flux_atm_ml(j)

      end if


 40   continue



c**** Compute the new concentration of Dissolved Inorganic Carbon in the
c**** Surface Ocean Mixed Layer. Units of TCO2 are moles/L.


      DO 50 j=1,nlat

      if(volume_ml(j).gt.0.0)then
        TCO2_ml(j) = (DIC_ml(j)/volume_ml(j))/1000.
        RCO2(j,1)= DIC_ml(j)/volume_ml(j)
!       RCO2(j,1)= (DIC_ml(j)-DICEQ(j))/volume_ml(j)
      end if


 50   continue

      return
      end

