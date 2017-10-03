
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

#include "BD2G04.COM"
#include "chem_para"
#include "chem_com"


!dimension p(nlon,nlat),dxyp(nlat)
!	----------------------------------------------
        common/ BACKGRGHG/GHGBGR(5)

!         print *,'From chemglobal '
!         print *,'P'
!         print *,p
!         print *,'dxyp'
!         print *,dxyp
	ptotal = 0.0
	do j=1,nlat
	  ptotal = ptotal + p(1,j)*dxyp(j)
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
     &                    *p(1,j)*dxyp(j)
        glbgas(k,2) = glbgas(k,2)
     &              + xn2o  (1,j,k)
     &                    *p(1,j)*dxyp(j)
        glbgas(k,3) = glbgas(k,3)
     &              + ch4   (1,j,k)
     &                    *p(1,j)*dxyp(j)
        glbgas(k,4) = glbgas(k,4)
     &              + cfc11 (1,j,k)
     &                    *p(1,j)*dxyp(j)
        glbgas(k,5) = glbgas(k,5)
     &              + cfc12 (1,j,k)
     &                    *p(1,j)*dxyp(j)
      enddo
	
	glbgas(k,1) = glbgas(k,1)*xxx
	glbgas(k,2) = glbgas(k,2)*xxx
	glbgas(k,3) = glbgas(k,3)*xxx
	glbgas(k,4) = glbgas(k,4)*xxx
	glbgas(k,5) = glbgas(k,5)*xxx

1	continue
!       print *,'GLBGAS'
!       do l=1,nlev
!        print *,(glbgas(l,i),i=1,5)
!       enddo

        do i=1,5
         totgas(i)=0.0
         adjfactor(i)=1.0
           do k=1,nlev
            totgas(i)=totgas(i)+glbgas(k,i)*dsig(k)
           enddo
        enddo
        tropmass = 28.97296245*1.e-3
!       print *,'TROPMASS=',tropmass
!       print *,'TOTGAS'
!        print *,(totgas(i),i=1,5)
          co2o=totgas(1)
     &      *tropmass/44.0098
          xxxo=totgas(2)
     &      *tropmass/44.0000
          yyyo=totgas(3)
     &      *tropmass/16.0426
          call newcon(co2n,co2o,xxxn,yyyn,xxxo,yyyo,
     &         GHGBGR(1),GHGBGR(2),GHGBGR(3))
          adjfactor(1)=co2n/co2o
          adjfactor(2)=xxxn/xxxo
          adjfactor(3)=yyyn/yyyo
!         print *,'N2O ',xxxn,xxxo,GHGBGR(2)
!         print *,'CH4 ',yyyn,yyyo,GHGBGR(3)
!         print *,'CO2 ',co2n,co2o,GHGBGR(1)

      return
       end

     
