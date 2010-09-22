
#include "ctrparam.h"

!	============================================================
!
!	CHEMAIRMASS.F:  Subroutine for calculating air mass
!			in MIT Global Chemistry Model
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
!	122095	Chien Wang	rev.	
!	080200	Chien Wang	repack based on CliChem3 & add cpp
!	051804	Chien Wang	rev.
!
!	==========================================================

!	=========================
	Subroutine chemairmass(p)
!	=========================

#include "chem_para"
#include "chem_com"
#include "BD2G04.COM"

!	----------------------------------------------------------

#if ( defined CPL_CHEM )

c-------------------------------
c  Calculate air mass:
c
      do k=1,nlev
      do i=1,n2dh
        airmass(i,1,k) = airmass0(i,1,k)*p(i,1)
      end do
      end do

#endif

      return
       end

!	================================
	Subroutine chemmass1(x11,xtotal)
!	================================

c==================================================================c
c                                                                  c
c   CHEMMASS1.F:  Subroutine for calculating total amount (mass)   c
c                   of tracers in MIT Global Chemistry Model       c
c        -------------------------------------------------         c
c   Author:		Chien Wang                                 c
c			MIT Joint Program on Science and Policy    c
c			    of Global Change                       c
c   Last Revised on:    August 8, 1995                             c
c                                                                  c
c==================================================================c

#include "chem_para"
#include "chem_com"

      dimension x11     (nlon,nlat,nlev)

!	--------------------------------------------------------

#if ( defined CPL_CHEM )

c-------------------------------
c  Calculate total amount of tracer:
c
      xtotal = 0.0
      do i=1,n3d
        xtotal  = xtotal
     &          + airmass(i,1,1)
     &          * x11(i,1,1)
      end do

#endif

      return
       end

!	=======================================
	Subroutine chemmass2(adjcoe,x11,xtotal)
!	=======================================

c==================================================================c
c                                                                  c
c   CHEMMASS2.F:  Subroutine for conpensating mass loss during     c
c		    simulation based on unified mapping 	   c
c                   in MIT Global Chemistry Model 	           c
c        -------------------------------------------------         c
c   Author:		Chien Wang                                 c
c			MIT Joint Program on Science and Policy    c
c			    of Global Change                       c
c   Last Revised on:    September 15, 1995                         c
c                                                                  c
c==================================================================c

#include "chem_para"
#include "chem_com"

      dimension x11     (nlon,nlat,nlev)

!	-------------------------------------------------

#if ( defined CPL_CHEM )

c---------------------------
c  Readjust tracer's mass:
c
      xtotal2 = 0.0
      do i=1,n3d
        xtotal2      = xtotal2
     &               + airmass(i,1,1)
     &               * x11(i,1,1)
      end do

      xgain  = (xtotal-xtotal2)
        if(xgain.gt.0) xgain = xgain*adjcoe

      xratio = xgain
     &       /float(nlat*nlev)

      do i=1,n3d
        x11(i,1,1) = max(0.0,x11(i,1,1)
     &             + xratio/airmass(i,1,1))
      end do

#endif

      return
       end
	
!	=======================================
	Subroutine chemmass3(adjcoe,x11,xtotal)
!	=======================================

c==================================================================c
c                                                                  c
c   CHEMMASS3.F:  Subroutine for conpensating mass loss during     c
c		    simulation based on unified mapping 	   c
c                   in MIT Global Chemistry Model 	           c
c		  Old chemmass2
c        -------------------------------------------------         c
c   Author:		Chien Wang                                 c
c			MIT Joint Program on Science and Policy    c
c			    of Global Change                       c
c   Last Revised on:    September 15, 1995                         c
c                                                                  c
c==================================================================c

#include "chem_para"
#include "chem_com"

      dimension x11     (nlon,nlat,nlev)

!	-------------------------------------------------------

#if ( defined CPL_CHEM )

c---------------------------
c  Readjust tracer's mass:
c
      xtotal2 = 0.0
      do i=1,n3d
        xtotal2      = xtotal2
     &               + airmass(i,1,1)
     &               * x11(i,1,1)
      end do

      xgain  = (xtotal-xtotal2)*adjcoe
c        if(xgain.gt.0) xgain = xgain*adjcoe

      xratio = xgain
     &       /float(nlat*nlev)

      do i=1,n3d
        x11(i,1,1) = max(0.0,x11(i,1,1)
     &             + xratio/airmass(i,1,1))
      end do

#endif

      return
       end

!	=======================================	
	Subroutine chemmass4(adjcoe,x11,xtotal)
!	=======================================

c==================================================================c
c                                                                  c
c   CHEMMASS4.F:  Subroutine for conpensating mass loss during     c
c                   simulation based on unified mass mapping       c
c		    plus N.H. extra share due to lbc loss          c
c                   in MIT Global Chemistry Model                  c
c        -------------------------------------------------         c
c   Author:             Chien Wang                                 c
c                       MIT Joint Program on Science and Policy    c
c                           of Global Change                       c
c   Last Revised on:    August 8, 1995                             c
c                                                                  c
c==================================================================c

#include "chem_para"
#include "chem_com"

      dimension x11     (nlon,nlat,nlev)

!	----------------------------------------------------

#if ( defined CPL_CHEM )

c---------------------------
c  Readjust tracer's mass:
c
      xtotal2 = 0.0
      do i=1,n3d
        xtotal2      = xtotal2
     &               + airmass(i,1,1)
     &               * x11(i,1,1)
      end do

      xxx    = (xtotal-xtotal2)
     &       /float(nlat*nlev)

      xratio = xxx
     &       *adjcoe

      xratio2= xxx
     &       *(1.0-adjcoe)/144.

      do i=1,n3d
        x11(i,1,1) = x11(i,1,1)
     &             + xratio/airmass(i,1,1)
      end do

      i=1
      ntropics = nlat/2
      do k=1,n_tropopause
      do j=ntropics+1,nlat1
        x11(i,j,k) 
     &             = x11(i,j,k)
     &             +(float(ntropics-j))**2
     &             *xratio2/airmass(i,j,k)
      end do
      end do

#endif

      return
       end

!	=============================================
	Subroutine chemmass6(adjyr,adjcoe,x11,xtotal)
!	=============================================

c==================================================================c
c                                                                  c
c   CHEMMASS2.F:  Subroutine for conpensating mass loss during     c
c		    simulation based on unified mapping 	   c
c		     and take away lossed mass directly according  c
c			to the tropospheric life time              c
c                   in MIT Global Chemistry Model 	           c
c        -------------------------------------------------         c
c   Author:		Chien Wang                                 c
c			MIT Joint Program on Science and Policy    c
c			    of Global Change                       c
c   Last Revised on:    September 15, 1995                         c
c                                                                  c
c==================================================================c

#include "chem_para"
#include "chem_com"

      dimension x11     (nlon,nlat,nlev)

!	----------------------------------------------------

#if ( defined CPL_CHEM )

c---------------------------
c  Readjust tracer's mass:
c
      xtotal2 = 0.0
      do i=1,n3d
        xtotal2      = xtotal2
     &               + airmass(i,1,1)
     &               * x11(i,1,1)
      end do

      xloss  = xtotal/(adjyr*8760.)		!8760 = 365d x 24h
						! and adjyr is in yr
      xgain  = (xtotal-xtotal2)
	if(xgain.gt.0) xgain = xgain*adjcoe

      xratio =(  xgain
     &         - xloss )
     &       /float(nlat*nlev)

      do i=1,n3d
        x11(i,1,1) = max(0.0, x11(i,1,1)
     &             + xratio/airmass(i,1,1))
      end do

#endif

      return
       end

!	==============================================
	Subroutine chemmass66(adjyr,adjcoe,x11,xtotal)
!	==============================================

c==================================================================c
c                                                                  c
c   CHEMMASS2.F:  Subroutine for conpensating mass loss during     c
c		    simulation based on unified mapping 	   c
c		     and take away lossed mass directly according  c
c			to the tropospheric life time              c
c                   in MIT Global Chemistry Model 	           c
c        -------------------------------------------------         c
c   Author:		Chien Wang                                 c
c			MIT Joint Program on Science and Policy    c
c			    of Global Change                       c
c   Last Revised on:    July 23, 1997                              c
c                                                                  c
c==================================================================c

#include "chem_para"
#include "chem_com"

      dimension x11     (nlon,nlat,nlev)

!	----------------------------------------------------

#if ( defined CPL_CHEM )

c---------------------------
c  Readjust tracer's mass:
c
c ======
c 102596
c  close this ocean sink to use ocean model calculate co2 uptake
!
! --- 080200 reopen under the cpp parameter
!
#if ( !defined CPL_OCEANCO2 )
      tropmass = 28.97296245*1.e-3/44.009
      co2ref   = 275./tropmass
     &		/(adjyr*8760.)		! 8760 = 365d x 24h
					!       and adjyr is in yr
#endif

      xtotal2 = 0.0
      do i=1,n3d
        xtotal2      = xtotal2
     &               + airmass(i,1,1)
     &               * x11(i,1,1)
      end do

c 072397:
c 101300:
#if ( !defined CPL_TEM )
      xloss  = bio_uptake *4.185692e17	! GTC/yr to 10-9kg/hr
#else
      xloss  = 0.0
#endif

#if ( !defined CPL_OCEANCO2 )
     &       + xtotal/(adjyr*8760.)
#endif

      xgain  = (xtotal-xtotal2)
	if(xgain.gt.0) xgain = xgain*adjcoe

      xratio =(  xgain
     &         - xloss )
     &       /float(nlat*nlev)

      do i=1,n3d
        x11(i,1,1) = x11(i,1,1)
     &             + xratio/airmass(i,1,1)
#if ( !defined CPL_OCEANCO2 )
     &             + co2ref
#endif
	if ( x11(i,1,1) .le. 0.0 ) x11(i,1,1) = 0.0

      end do

#endif

      return
       end

