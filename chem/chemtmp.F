
#include "ctrparam.h"

!	============================================================
!
!	CHEMTMP1.F:   Subroutine for calculating some temperary
!			parameters for monthly averaging
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
!	103194	Chien Wang	rev.	
!	080200	Chien Wang	repack based on CliChem3 & add cpp
!
!	==========================================================

	subroutine chemtmp1 (dth,airmass0,p,pvv,xx1)

#include "chem_para"
#include "chem_tmp"
#include "BD2G04.COM"

	dimension airmass0(nlon,nlat,nlev)
	dimension pvv     (nlon,nlat,nlev)
	dimension xx1     (nlon,nlat,nlev)

!	----------------------------------------------------------

#if ( defined CPL_CHEM )

c
c  CFC11 N-S transport amount in 10^-9 kg
c
      i = 1
      do 1 k=1,nlev
        cfcnsf(k) = cfcnsf(k)                
     &           +( max(0.0,pvv(i,13,k))
     &             *airmass0(i,12,k)*p(i,12)*xx1(i,12,k)
     &            + min(0.0,pvv(i,13,k))
     &             *airmass0(i,13,k)*p(i,13)*xx1(i,13,k) )
     &             /dyv(13)*dth   *1.e-9      !kg
c     &            +( (pvv(i,14,k)+pvv(i,13,k))/dyv(13)
c     &               *airmass0(i,13,k)*p(i,13)*xx1(i,13,k)
c     &              -(pvv(i,13,k)+pvv(i,12,k))/dyv(12)
c     &               *airmass0(i,12,k)*p(i,12)*xx1(i,12,k) )
c     &            *dth*0.5e-9                 !kg
1     continue

#endif

      return
       end

c====================================================

      subroutine chemtmp2

c==================================================================c
c                                                                  c
c   CHEMTMP1.F:   Subroutine for calculating some temperary        c
c                   parameters for monthly averaging               c
c        -------------------------------------------------         c
c   Author:		Chien Wang                                 c
c			MIT Joint Program on Science and Policy    c
c			    of Global Change                       c             
c   Last Revised on:    October 31, 1994                           c
c                                                                  c
c==================================================================c

#include "chem_para"
#include "chem_tmp"
#include "BD2G04.COM"

#if ( defined CPL_CHEM )

      write(168)cfcnsf

      do k=1,nlev
        cfcnsf(k) = 0.0
      end do

#endif

      return
       end

     
