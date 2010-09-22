C $Header: /home/igsm-repository/ATM/rtgases.F,v 1.1 2010/08/30 15:36:35 jscott Exp $
C $Name:  $

#include "ctrparam.h"

!	==========================================================
!
!	RTGASES.F:   Subroutine for reading monthly-mean tracer mixing
!			ratios from climate chemistry run for 
!                             a given year
!           This subroutine is called if READGHG.eq.1
!	----------------------------------------------------------
!
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

       subroutine rtgases(CO2,JMONTH)

#include "chem_para"
#include "chem_com"

      DIMENSION AMONTH(12),cff11(3)
      CHARACTER*4 AMONTH,JMONTH,RMONTH
      DATA AMONTH/'JAN','FEB','MAR','APR','MAY','JUNE','JULY','AUG',
     *  'SEP','OCT','NOV','DEC'/  
      data nread /0/
      data ifirst/1/
      data cff11 /1.04,1.36,2.75/

!
#ifdef PREDICTED_GASES
!
c============================
c  Unit here for gases is
c	10^-9 kg(c)/kg(air), or ppb(m):
c
!     if(CO2.le.0.0.or.CO2.gt.2.5)then
!       print *,' rtgases CO2=',CO2
!       stop
!     endif
! 100 continue
      if(nread.eq.0)then
        read(569),aexp1,iy1
        read(669),aexp2,iy2
        if(abs(aexp1-aexp2).gt.0.1)then
         print *,'Wrong montly files'
         print *,'aexp1=',aexp1,' aexp2=',aexp2
         stop
        else
         if(ifirst.eq.1)then
           print *,' Run with GHGs concentrations for '
           print *,' year ',iy1,' from run ', aexp1
           ifirst=0
         endif
        endif
      endif
      nread=nread+1
      print *,' nread=',nread
      read(569)cfc11
      read(569)cfc12
      read(569)xn2o
      read(569)o3
      read(569)co
      read(569)zco2
      read(569)ho		!useless
      read(569)xno		!useless
      read(569)xno		!useless
      read(569)ch4
      read(569)so2		!useless
      read(569)h2so4
      read(569)sviod		!optical depth of S(VI)
      read(569)bcm
      read(569)ocm
      read(569)bcodm
      read(569)ocodm

      read(669)hfc134a
      read(669)pfc
      read(669)sf6
      do l=1,nlev
      do j=1,nlat
      do i=1,nlon
       cfc11(i,j,l)=cfc11(i,j,l)+cff11(1)*hfc134a(i,j,l)
     &             +cff11(2)*pfc(i,j,l)+cff11(3)*sf6(i,j,l)
      enddo
      enddo
      enddo


      RMONTH=AMONTH(nread)
      if(RMONTH.ne.JMONTH)then
        print *,' RTGASES NREAD=',nread,' ',RMONTH,JMONTH
        stop
      endif
      if(nread.eq.12)then
          rewind (569)
          rewind (669)
          nread=0
      endif
!     do 10 k=1,nlev
!       do 10 j=1,nlat
!        do 10 i=1,nlon
!          zco2(i,j,k)=CO2*zco2(i,j,k)
!  10 continue


!
#endif
!
        return
        end
