
#include "ctrparam.h"

!	============================================================
!
!	CHEMBUDGET.F:  Subroutine for calculating total amount of
!			 tracers in MIT Global Chemistry Model
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
!	022999	Chien Wang	rev.	
!	080100	Chien Wang	repack based on CliChem3 & add cpp
!	092001	Chien Wang	add bc and oc
!	051804	Chien Wang	rev. for 46x11
!
!	==========================================================

      subroutine chembudget (p)

#include "chem_para"
#include "chem_com"
#include "BD2G04.COM"

      dimension tlayer(nlev)
	
	real :: totalair, tropair, totalf11, totalch4, totalco
	real :: totalo3,  totaloh, totalso2, totalnox
	real :: totalbc,  totaloc

!	---------------------------------------------------------
#if ( defined CPL_CHEM )

      totalair = 0.0
      tropair  = 0.0
      totalf11 = 0.0
	
c add ch4 and co budget 022996:

	totalch4 = 0.0
	totalco  = 0.0

c and others 030996:

	totalo3  = 0.0
	totaloh  = 0.0
	totalso2 = 0.0
	totalnox = 0.0

! === add bc and oc
	totalbc  = 0.0
	totaloc  = 0.0

      i = 1 

c Troposphere first:

      do 111 k=1,n_tropopause
      do 112 j=1,nlat
        amass    = airmass0(i,j,k)*p(i,j)
        totalair = totalair
     &           + amass 

        totalf11 = totalf11 
     &           + amass * cfc11(i,j,k)
        totalch4 = totalch4 
     &           + amass * ch4(i,j,k)
        totalco  = totalco 
     &           + amass * co(i,j,k)
	totalo3  = totalo3
     &           + amass * o3(i,j,k)
	totaloh  = totaloh
     &           + amass * ho(i,j,k)
	totalso2 = totalso2
     &           + amass * so2(i,j,k)
	totalnox = totalnox
     &           + amass * 14.0*(xno(i,j,k)/30.0
     &		               + xno2(i,j,k)/46.0)  !Convert to N
        totalbc  = totalbc 
     &           + amass * bcarbon(i,j,k)
	totaloc  = totaloc
     &           + amass * ocarbon(i,j,k)
112   continue
111   continue

	tropair = totalair

c Whole Atmosphere:

      do 115 k=n_tropopause+1,nlev
      do 116 j=1,nlat
        amass    = airmass0(i,j,k)*p(i,j)
        totalair = totalair
     &           + amass 

        totalf11 = totalf11 
     &           + amass * cfc11(i,j,k)
        totalch4 = totalch4 
     &           + amass * ch4(i,j,k)
        totalco  = totalco 
     &           + amass * co(i,j,k)
        totalbc  = totalbc 
     &           + amass * bcarbon(i,j,k)
	totaloc  = totaloc
     &           + amass * ocarbon(i,j,k)
116   continue
115   continue

	totalf11 = totalf11 * 1.e-15		!10^6 kg
	totalch4 = totalch4 * 1.e-18		!TgCH4
	totalco  = totalco  * 1.e-18		!TgCO
	totalbc  = totalbc  * 1.e-18		!TgC
	totaloc  = totaloc  * 1.e-18		!TgC

c troposphere only:
	totalo3  = totalo3  * 1.e-18		!TgO3
	totaloh  = totaloh  			!10^-9 kg
	totalso2 = totalso2 * 1.e-18		!TGSO2
	totalnox = totalnox * 1.e-18		!TGN

      write(167,100) iyearchem, monthchem,
     &               totalair,tropair,
     &               totalf11,
     &               totalch4,
     &               totalco,
     &               totalo3,
     &               totaloh,
     &               totalso2,
     &               totalnox,
     &		     totalbc,
     &		     totaloc

c100   format(" Yr ",i4," Mon ",i3,
c     &       " Total air mass = ",e12.3," kg",/
c     &       " Total F11 = ",f12.3," 10^6 kg",/
c     &       " Total CH4 = ",f12.3," TgCH4",/
c     &       " Total CO  = ",f12.3," TgCO")

100    format(i5,i5,2e12.4,3f12.3,/ 
     &        10x,f12.3,e12.4,2f12.3,2e12.4)

#endif

      return
       end
