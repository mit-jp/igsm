C $Header: /home/igsm-repository/ATM/chemglobal.F,v 1.1 2010/08/30 15:36:35 jscott Exp $
C $Name:  $

#include "ctrparam.h"

!	==========================================================
!
!	CHEMGLOBAL.F:	Subroutine for calculating global averaged
!		  	 vertical profiles of trace gases
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
!	080100	Chien Wang	repack based on CliChem3 & add cpp
!
!	==========================================================

	subroutine chemglobal(p)
!	========================

#include "chem_para"
#include "chem_com"

	dimension p(nlon,nlat)
!	----------------------------------------------

	ptotal = 0.0
	do j=1,nlat
	  ptotal = ptotal + p(1,j)
	enddo
	xxx = 1./ptotal

	do 1 k=1,nlev

      	  glbgas(k,1) = 0.0		!CO2
      	  glbgas(k,2) = 0.0		!N2O
      	  glbgas(k,3) = 0.0		!CH4
      	  glbgas(k,4) = 0.0		!F11
      	  glbgas(k,5) = 0.0		!F12

      do j=1,nlat
        glbgas(k,1) = glbgas(k,1)
     &              + zco2  (1,j,k)
     &                    *p(1,j)
        glbgas(k,2) = glbgas(k,2)
     &              + xn2o  (1,j,k)
     &                    *p(1,j)
        glbgas(k,3) = glbgas(k,3)
     &              + ch4   (1,j,k)
     &                    *p(1,j)
        glbgas(k,4) = glbgas(k,4)
     &              + cfc11 (1,j,k)
     &                    *p(1,j)
        glbgas(k,5) = glbgas(k,5)
     &              + cfc12 (1,j,k)
     &                    *p(1,j)
      enddo
	
	glbgas(k,1) = glbgas(k,1)*xxx
	glbgas(k,2) = glbgas(k,2)*xxx
	glbgas(k,3) = glbgas(k,3)*xxx
	glbgas(k,4) = glbgas(k,4)*xxx
	glbgas(k,5) = glbgas(k,5)*xxx

1	continue

      return
       end

     
