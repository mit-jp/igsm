C $Header: /home/igsm-repository/CHEM/chemstrat.F,v 1.1 2010/09/02 17:56:19 jscott Exp $
C $Name:  $

#include "ctrparam.h"

!	============================================================
!
!	CHEMEMISSION.F:  Subroutine for calculating stratospheric
!			  decuction of MIT Global Chemistry Model
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
!	102694	Chien Wang	rev.	
!	080200	Chien Wang	repack based on CliChem3 & add cpp
!
!	==========================================================

      subroutine chemstrat (dtb)

#include "chem_para"
#include "chem_com"
#include "BD2G04.COM"

!	----------------------------------------------------------

#if ( defined CPL_CHEM )

c------------------------
c  Deduction of CFC11
c    in stratosphere:
c
      dta = dtb *3.0
      i=1
      k=nlev

c 4/19/95:

      tf11=1./3.1536e8*dta	!10 yr
c      tf12=1./1.173139e9*dta
c      tn2o=1./1.110067e9*dta

      tf12=tf11
      tn2o=tf11

      do j=1,nlat
c        cfc11sd     = cfc11sd 
c     &              + airmass0(i,j,k)*p(i,j)
c     &              * tf11 * cfc11(i,j,k)
        cfc11(i,j,k)= cfc11(i,j,k)
     &              * (1.-tf11)

c        cfc12sd     = cfc12sd 
c     &              + airmass0(i,j,k)*p(i,j)
c     &              * tf12 * cfc12(i,j,k)
        cfc12(i,j,k)= cfc12(i,j,k)
     &              * (1.-tf12)

c        xn2osd      = xn2osd 
c     &              + airmass0(i,j,k)*p(i,j)
c     &              * tn2o * xn2o(i,j,k)
        xn2o (i,j,k)= xn2o(i,j,k)
     &              * (1.-tn2o)
      enddo

#endif

      return
       end
