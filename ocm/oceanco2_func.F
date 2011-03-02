
#include "ctrparam.h"

!	==========================================================
!
!	OCEANCO2_FUNC.F:  Utility subroutines of ocean CO2 model
!
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
!	080200	Chien Wang	repack based on CliChem3 and add cpp.
!
!	==========================================================

!	==================
	subroutine avegary
!	==================

!#include "chem_para"
#include "ODIFF.COM"
#include "com_carb"


        ncallgary=ncallgary+1
        do j=1,nlat
	pC_atmav(j)=pC_atmav(j)+pC_atm(j)
	pC_mlaav(j)=pC_mlaav(j)+pC_mla(j)
	Flux_atm_mlav(j)=Flux_atm_mlav(j)+Flux_atm_ml(j)
	Hav(j)=Hav(j)+H(j)
	DIC_mlav(j)=DIC_mlav(j)+DIC_ml(j)
	TCO2_mlav(j)=TCO2_mlav(j)+TCO2_ml(j)
	Alkcav(j)=Alkcav(j)+AlkC(j)
	AlkBav(j)=AlkBav(j)+AlkB(j)
	AlkSav(j)=AlkSav(j)+AlkS(j)
	AlkWav(j)=AlkWav(j)+AlkW(j)
	CO2av(j)=CO2av(j)+CO2_ml(j)
	CO3av(j)=CO3av(j)+CO3_ml(j)
         do l=1,lmo
           Rco2av(j,l)=Rco2av(j,l)+rco2(j,l)
         end do !       l
        end do  !       j

cprint *,'from avegary ncallgary=',ncallgary

	return
	end

!	==================================
	subroutine wrgary
!	==================================

!#include "chem_para"
#include "ODIFF.COM"
#include "com_carb"
!print *,'from wrgary ncallgary=',ncallgary
        do j=1,nlat
	pC_atmav(j)=pC_atmav(j)/ncallgary
	pC_mlaav(j)=pC_mlaav(j)/ncallgary
	Flux_atm_mlav(j)=Flux_atm_mlav(j)/ncallgary
	Hav(j)=Hav(j)/ncallgary
	DIC_mlav(j)=DIC_mlav(j)/ncallgary
	TCO2_mlav(j)=TCO2_mlav(j)/ncallgary
	CO2av(j)=CO2av(j)/ncallgary
	CO3av(j)=CO3av(j)/ncallgary
	Alkcav(j)=Alkcav(j)/ncallgary
	AlkBav(j)=AlkBav(j)/ncallgary
	AlkSav(j)=AlkSav(j)/ncallgary
	AlkWav(j)=AlkWav(j)/ncallgary
         do l=1,lmo
           RCO2av(j,l)=RCO2av(j,l)/ncallgary
         end do !       l
        end do  !       j
	write (602)pC_atmav,pC_mlaav,Flux_atm_mlav
	write (602)DIC_mlav,TCO2_mlav,AlkBav,AlkSav,AlkWav
	write(602)Hav,Alkcav,CO2av,Co3av
	write(602)RCO2av
 5001 format(24f5.1)
 5002 format(24e10.2)


	return
	end

!	==============================
	subroutine zerogary
!	==============================

!#include "chem_para"
#include "ODIFF.COM"
#include "com_carb"


cprint *,' from zerogary ncallgary=',ncallgary,' nlat=',nlat
        do j=1,nlat
	pC_atmav(j)=0.
	pC_mlaav(j)=0.
	Flux_atm_mlav(j)=0.
	Hav(j)=0.
	DIC_mlav(j)=0.
	TCO2_mlav(j)=0.
        AlkCav(j)=0.
        AlkBav(j)=0.
        AlkSav(j)=0.
        AlkWav(j)=0.
        CO2av(j)=0.
        CO3av(j)=0.
         do l=1,lmo
           RCO2av(j,l)=0.
         end do !       l
        end do  !       j
	ncallgary=0.


	return
	end
