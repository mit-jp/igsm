C $Header: /home/igsm-repository/ATM/chemtemco2.F,v 1.1 2010/08/30 15:36:35 jscott Exp $
C $Name:  $

#include "ctrparam.h"

!	==========================================================
!
!	CHEMTEMCO2.F:  Subroutine for deduce CO2 amount caused by
!                  land (TEM)  uptake of MIT Global Chemistry Model
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
!       101200                  developed from chemoceanco2
!
!	==========================================================

!	==============================
	subroutine chemtemco2 (fnep)
!	==============================

#include "chem_para"
#include "chem_com"
#include "BD2G04.COM"

	dimension fnep(nlat)	!flux toward land of carbon in Gt

!	------------------------------------------------------------

#if ( defined CPL_CHEM) && ( defined CPL_TEM )

	xconv = 44.e6           !from mole(C) to 10^-9 kg(CO2)
     &   *1.e15/12.            !from Gt (C) to mole(C)
c----
	i     = 1

c 052295:
	ktop  = 2
     
	do 1 j=3,nlat	! throwout 2 south pole points
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
	  yamount = fnep(j)*xconv
c
c Map CO2 sink into two layer:
c
	  xemi     = max(0.0, (xmass - yamount))
     &                            * tmass1     !ppbm

	do k=1,ktop
	  zco2(i,j,k) = xemi
	enddo

1	continue

#endif

	return
	 end

