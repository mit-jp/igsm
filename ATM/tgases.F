
#include "ctrparam.h"

!	==========================================================
!
!	TGASES.F:   Subroutine for reading monthly-mean tracer mixing
!			ratios to carryout offline calculation
!			
!			Version allows reading continually.
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

       subroutine tgases (CO2,JMONTH)

#include "chem_para"
#include "chem_com"

      DIMENSION AMONTH(12),cff11(3)
      CHARACTER*4 AMONTH,JMONTH,RMONTH
      DATA AMONTH/'JAN','FEB','MAR','APR','MAY','JUNE','JULY','AUG',
     *  'SEP','OCT','NOV','DEC'/
      data nread /0/
      data ifirst/1/
      data cff11 /1.04,1.36,2.75/


c============================
c  Unit here for gases is
c	10^-9 kg(c)/kg(air), or ppb(m):
c
      if(CO2.le.0.0.or.CO2.gt.1.5)then
        print *,' tgases CO2=',CO2
        stop
      endif
  100 continue
      nread=nread+1
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

      read(679)hfc134a
      read(679)pfc
      read(679)sf6
#ifdef PREDICTED_GASES
      do l=1,nlev
      do j=1,nlat
      do i=1,nlon
       cfc11(i,j,l)=cfc11(i,j,l)+cff11(1)*hfc134a(i,j,l)
     &             +cff11(2)*pfc(i,j,l)+cff11(3)*sf6(i,j,l)
      enddo
      enddo
      enddo
#endif


      nr1=nread/12
      nr=nread-nr1*12
      if(nr.eq.0)nr=12
      RMONTH=AMONTH(nr)
      if(RMONTH.ne.JMONTH)then
         print *,' DISAGREMENT in TGASES'
         print *,nread,nr,RMONTH,JMONTH
         stop
      endif
      if(ifirst.eq.1)then
        print *,' TGASES NREAD=',nread,nr,RMONTH,JMONTH
        print *,'CO2'
c       print *,(zco2(1,12,l),l=1,nlev)
c       ifirst=0
      endif
ctropmass = 28.97296245

c        cfc11 = cfc11*tropmass/137.3675*1.e3 !pptv
c        cfc12 = cfc12*tropmass/120.9054*1.e3 !pptv
c        xn2o  = xn2o *tropmass/44.0000       !ppbv
c        o3    = o3   *tropmass/47.9982       !ppbv
c        zco2  = zco2 *tropmass/44.0098*1.e-3 !ppmv
c        co    = co   *tropmass/28.0104       !ppbv
c        ch4   = ch4  *tropmass/16.0426*1.e-3 !ppmv

        return
        end
