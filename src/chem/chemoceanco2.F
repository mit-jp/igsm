
#include "ctrparam.h"

!	==========================================================
!
!	CHEMOCEANCO2.F:  Subroutine for deduce CO2 amount caused by
!			   ocean uptake of MIT Global Chemistry Model
!
!	----------------------------------------------------------
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
!	102596	Chien Wang	rev.
!	080100	Chien Wang	repack based on CliChem3 & add cpp
!
!	==========================================================

!	==============================
	subroutine chemoceanco2 (fco2)
!	==============================

#include "chem_para"
#include "chem_com"
#include "BD2G04.COM"

	dimension fco2(nlat)	!flux toward ocean of co2 in mole(C)

!	------------------------------------------------------------

#if ( defined CPL_CHEM) && ( defined CPL_OCEANCO2 )

	xconv = 44.e6         !from mole(C) to 10^-9 kg(CO2)
c----
	i     = 1

c 052295:
	ktop  = 2
     
#if ( N_LAT == 24 )
	do  j=3,nlat	! throwout 2 south pole points
#endif
#if ( N_LAT == 46 )
	do  j=4,nlat	! throwout 3 south pole points
#endif
c
c  Calculate total mass of air and CO2:
c
	  tmass    = 0.0                !Total air mass in kg
	  xmass    = 0.0                !Total co2 mass in 10^-9 kg
	do k=1,ktop
	  tmass    = tmass + airmass(i,j,k)
	  xmass    = xmass
     &             + airmass(i,j,k)*zco2(i,j,k)
	enddo
	  tmass1   = 1./tmass
c
c Convert flux from mole(C) to ppb(CO2) in mass:
c
	  yamount = fco2(j)*xconv
c
c Map CO2 sink into two layer:
c
	  xemi     = max(0.0, (xmass - yamount))
     &                            * tmass1     !ppbm

	do k=1,ktop
	  zco2(i,j,k) = xemi
	enddo

	enddo

#endif

	return
	 end

