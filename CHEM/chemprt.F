C $Header: /home/igsm-repository/CHEM/chemprt.F,v 1.1 2010/09/02 17:56:19 jscott Exp $
C $Name:  $

#include "ctrparam.h"

!	============================================================
!
!	CHEMPRT.F:  	Subroutine for writing data into files
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
!	030995	Chien Wang	rev.	
!	080200	Chien Wang	repack based on CliChem3 & add cpp
!
!	==========================================================

c
      subroutine chemprt
c     ====================

#include "chem_para"
#include "chem_com"

#if ( defined CPL_CHEM )

      write(166)cfc11
      write(166)cfc12
      write(166)xn2o

      write(171)o3
      write(171)o1d
      write(171)ho
      write(171)ho2
      write(171)h2o2
       
      write(172)xno
      write(172)xno2
      write(172)xno3
      write(172)xn2o5
      write(172)hno3

      write(173)co
      write(173)zco2
      write(173)ch4
      write(173)ch2o

      write(174)so2
      write(174)hoso2
      write(174)so3
      write(174)h2so4

c        write(6,*)atomo(1,1,1),
c     & atomh (1,1,1),
c     & ch3   (1,1,1),
c     & cho   (1,1,1),
c     & ch3o  (1,1,1),
c     & ch3o2 (1,1,1),
c     & ch3o2h(1,1,1),

c	do j=1,nlat
c	write(6,*)"sviod = ",sviod(1,j,1)
c	enddo

#endif

       return
        end

