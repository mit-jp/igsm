
#include "ctrparam.h"

!	============================================================
!
!	CHEMEMISSION.F:  Subroutine for calculating surface emission
!				of tracers of MIT Global Chemistry Model
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
!	031097	Chien Wang	brought in NEM
!	052200	Chien Wang	rev.	
!	080100	Chien Wang	repack based on CliChem3 & add cpp
!	101800	Chien Wang	delete if_3gases & use cpp
!	121800	Chien Wang	move xnoxltnt init to cheminit
!	092001	Chien Wang	add bc and oc
!	051804	Chien Wang	rev. for 46x11
!
!	==========================================================

!
	subroutine chememission (tinterval)
!	===================================

#include "chem_para"
#include "chem_com"
#include "BD2G04.COM"
#include "chem_meta"

#if ( defined CPL_NEM )
C For Emission
c === 031097
      common/EMFORCHIEN/ECH4CHIEN(nlat),EN2OCHIEN(nlat)
c For Emission
#endif

	dimension x11(nlev)
#if ( defined ACCRI )
        COMMON U,V,T,P,Q
        dimension  emi_accri(nlev)
#endif
	
	real emi_new
        data ifirst /0/

!	----------------------------------------------------------
!
#if ( defined CPL_CHEM )

! --- call meta model
#if ( defined CPL_META )
	call chemmeta
#endif
        if (ifirst.eq.0) then
#if ( defined ACCRI )
        print *,'ACCRI is true in chem with H2O emissions'
#endif
         print *,'First year of emissions from chememi=',myyear
	 print *,'CO2 emissions for ',myyear
         print *,(edailyco2(1,j,myyear),j=1,nlat)
! print *,'CO2 emissions for ',myyear-1
!        print *,(edailyco2(1,j,myyear-1),j=1,nlat)
!        print *,'n_urban for ',myyear
!        do k=1,3
!         print *,k
!         print *,(n_urban(k,j,myyear),j=1,nlat)
!        enddo
         ifirst=1
        endif
	do k=1,nlev
	  x11(k) = 0.0
	enddo

	xconv   = tinterval*1.e15	!time & unit conversion
c
c NO by lightning 062795:
c
	yearemi = xnoxltnt(myyear)
     &          * xnoxltnm(mymonth)
	ynoltn  = yearemi*xconv
     &          * 12.0			!forgot to set it as an annual

c----
      i    = 1

c 052295:
      ktop = 2
c
      kup  = ktop + 1
     
      do 1 j=2,nlat1                     !Only apply to emission point
c      do 1 j=1,nlat

        eratio = ehpbl(i,j)

        tmass    = 0.0                !Total air mass in kg
        tmassx   = 0.0
      do k=1,ktop
        tmass    = tmass + airmass(i,j,k)
      enddo
        tmassx   = tmass + airmass(i,j,kup)*eratio

        tmass1   = 1./tmass
        tmass1x  = 1./tmassx

c
c 1. N2O:
c
      do k=1,kup
        x11(k) = xn2o (i,j,k)
      enddo

#if ( defined CPL_NEM )
C For Emission
c === 031097
        xemin2o = edailyn2o(i,j,myyear)
c       conversion of NEM emission from kg of N to ppbm of N2O
     &          + EN2OCHIEN(j) *1.e9 *44.0/28.0  !convert to ppb(m)
cold &          + EN2OCHIEN(j) *1.e9 *44.0/28.0*9.333/13.  !convert to ppb(m)

      call chememission2(1,j,ktop,0.0,
     &                   tmass1,xemin2o,x11)
#else
!APS 062210
        xemin2o = edailyn2o(i,j,myyear)
c       fixed  NEM of 10.8 *1.e9  kg of N 
     &          + 10.8/365.*1.e9 *44.0/28.0  !convert to ppb(m)
      call chememission2(1,j,ktop,0.0,
!!   &                   tmass1,edailyn2o(i,j,myyear),x11)
     &                   tmass1,xemin2o,x11)
#endif

      do k=1,kup
        xn2o(i,j,k) = x11(k)
      enddo

c
c 2. CH4:
c
      do k=1,kup
        x11(k) = ch4 (i,j,k)
      enddo

#if ( defined CPL_NEM )
C For Emission
c === 031097
        xemich4 = edailych4(i,j,myyear)
c       conversion of NEM emission from kg of CH4 to ppbm of CH4
     &          + ECH4CHIEN(j) *1.e9          !convert to ppb(m)
cold &          + ECH4CHIEN(j) *1.e9*126./136.          !convert to ppb(m)

      call chememission2(1,j,ktop,0.0,
     &                   tmass1,xemich4,x11)
#else
!APS 062210
        xemich4 = edailych4(i,j,myyear)
c       fixed NEM emission  144 *1.e9 kg of CH4
     &          + 144./365.*1.e9          !convert to ppb(m)
      call chememission2(1,j,ktop,0.0,
     &                   tmass1,xemich4,x11)
!!!  &                   tmass1,edailych4(i,j,myyear),x11)
#endif

      do k=1,kup
        ch4(i,j,k) = x11(k)
      enddo

c
c 3. CO:
c
      do k=1,kup
        x11(k) = co (i,j,k)
      enddo
#if ( defined CPL_META )
      ! --- Convert flux in kg/day into
      ! --- 	daily emission in 10^-9 kg
      ! --- Note the flux has been multiplied by n
      ! ---	already in chemmeta
      emi_new = 0.0
      do nn = 1,3
        emi_new = emi_new + results_meta(11,nn,j)
      end do
!      emi_new = edailyco(i,j,myyear)*(1.0 - alpha_co(j,myyear))
      emi_new = edailyco(i,j,myyear) - edailyuco(i,j,myyear)
     &        + emi_new*1.e9
     
      call chememission2(1,j,ktop,0.0,
     &                   tmass1,emi_new,x11)
#else
      call chememission2(1,j,ktop,0.0,
     &                   tmass1,edailyco(i,j,myyear),x11)
#endif

      do k=1,kup
        co(i,j,k) = x11(k)
      enddo

c
c 4. NOx:
c
      do k=1,kup
        x11(k) = xno  (i,j,k)
c        x11(k) = xno2 (i,j,k)
      enddo

#if ( defined CPL_META )
      ! --- Convert flux in kg/km^2/day into
      ! --- 	daily emission in 10^-9 kg
      ! --- Note the flux has been multiplied by n
      ! ---	already in chemmeta
      emi_new = 0.0
      do nn = 1,3
        emi_new = emi_new + results_meta(2,nn,j)
      end do
!      emi_new = edailynox(i,j,myyear)*(1.0 - alpha_nox(j,myyear))
      emi_new = edailynox(i,j,myyear) - edailyunox(i,j,myyear)
     &        + emi_new*1.e9
     
      call chememission2(1,j,ktop,0.0,
     &                   tmass1,emi_new,x11)
#else
      call chememission2(1,j,ktop,0.0,
     &                   tmass1,edailynox(i,j,myyear),x11)
#endif

      do k=1,kup
        xno(i,j,k) = x11(k)
c        xno2(i,j,k) = x11(k)
      enddo

c
c  NO - Lightning contribution:
c
      do k=1,n_tropopause  ! note: ktrop = 7 in 9 & 11 layer model
        x11(k) = xno  (i,j,k)
      enddo
      xfunc    = xnoxltnd(i,j,mymonth)

      call chememiltn(1,j,7,
     &                   xfunc,ynoltn,x11)

      do k=1,n_tropopause
        xno(i,j,k) = x11(k)
      enddo

c-- Use chememission3 below:

c
c 5. CFC11:
c
      do k=1,kup
        x11(k) = cfc11 (i,j,k)
      enddo

      call chememission3(1,j,ktop,eratio,
     &                   tmass1x,edailyf11(i,j,myyear),x11)

      do k=1,kup
        cfc11(i,j,k) = x11(k)
      enddo
c
c 6. CFC12:
c
      do k=1,kup
        x11(k) = cfc12 (i,j,k)
      enddo

      call chememission3(1,j,ktop,eratio,
     &                   tmass1x,edailyf12(i,j,myyear),x11)

      do k=1,kup
        cfc12(i,j,k) = x11(k)
      enddo
c
c 7. SO2:
c
      do k=1,kup
        x11(k) = so2 (i,j,k)
      enddo

      xxx=edailyso2(i,j,myyear)
      call chememission3(1,j,ktop,eratio,
     &                   tmass1x,xxx,x11)

      do k=1,kup
        so2(i,j,k) = x11(k)
      enddo

c
c 8. CO2:
c
      do k=1,kup
        x11(k) = zco2 (i,j,k)
      enddo

      call chememission2(1,j,ktop,0.0,
     &                   tmass1,edailyco2(i,j,myyear),x11)

      do k=1,kup
        zco2(i,j,k) = x11(k)
      enddo

! === if hfc, pfc, and sf6 are included:
#if ( defined INC_3GASES )

! === 032698
! === 9. HFC134a:
! ===
      do k=1,kup
        x11(k) = hfc134a(i,j,k)
      enddo

      call chememission2(1,j,ktop,0.0,
     &                   tmass1,edailyhfc134a(i,j,myyear),x11)

      do k=1,kup
        hfc134a(i,j,k) = x11(k)
      enddo

! === 
! === 10. PFC:
! ===
      do k=1,kup
        x11(k) = pfc(i,j,k)
      enddo

      call chememission2(1,j,ktop,0.0,
     &                   tmass1,edailypfc(i,j,myyear),x11)

      do k=1,kup
        pfc(i,j,k) = x11(k)
      enddo

! === 
! === 11. SF6:
! ===
      do k=1,kup
        x11(k) = sf6(i,j,k)
      enddo

      call chememission2(1,j,ktop,0.0,
     &                   tmass1,edailysf6(i,j,myyear),x11)

      do k=1,kup
        sf6(i,j,k) = x11(k)
      enddo

#endif
! ===

! 
! === 12. Black Carbon:
!
      do k=1,kup
        x11(k) = bcarbon(i,j,k)
      enddo

      call chememission2(1,j,ktop,0.0,
     &                   tmass1,edailybc(i,j,myyear),x11)

      do k=1,kup
        bcarbon(i,j,k) = x11(k)
      enddo

! 
! === 13. Organic Carbon:
!
      do k=1,kup
        x11(k) = ocarbon(i,j,k)
      enddo

      call chememission2(1,j,ktop,0.0,
     &                   tmass1,edailyoc(i,j,myyear),x11)

      do k=1,kup
        ocarbon(i,j,k) = x11(k)
      enddo

#if ( defined ACCRI )
!       print *,"ACCRI emissions"
c
c 1. CO2:
c
      do k=1,nlev
        x11(k) = zco2 (i,j,k)
      enddo
       call int_accri_emi(j,ps4accri(i,j),co2accri,emi_accri,
     &                     mymonth,myyear)

      call chememissionA(1,j,nlev,
     &                   emi_accri,x11)

      if(j.eq.-23) then
      print *,'CO2 before and after chememissionA'
      do k=1,nlev
        print *,k,(zco2 (i,j,k)-x11(k)),emi_accri(k)
      enddo
      endif
      do k=1,nlev
        zco2(i,j,k) = x11(k)
      enddo
c
c 2. CO:
c
      do k=1,nlev
        x11(k) =  co  (i,j,k)
      enddo
       call int_accri_emi(j,ps4accri(i,j),coaccri,emi_accri,
     &                     mymonth,myyear)

      call chememissionA(1,j,nlev,
     &                   emi_accri,x11)
      do k=1,nlev
         co (i,j,k) = x11(k)
      enddo
c
c
c 3. NO:
c
      do k=1,nlev
        x11(k) =  xno  (i,j,k)
      enddo
       call int_accri_emi(j,ps4accri(i,j),NOaccri,emi_accri,
     &                     mymonth,myyear)

      call chememissionA(1,j,nlev,
     &                   emi_accri,x11)
      do k=1,nlev
         xno (i,j,k) = x11(k)
      enddo
c
c
c 4. NO2:
c
      do k=1,nlev
        x11(k) =  xno2  (i,j,k)
      enddo
       call int_accri_emi(j,ps4accri(i,j),NO2accri,emi_accri,
     &                     mymonth,myyear)

      call chememissionA(1,j,nlev,
     &                   emi_accri,x11)
      do k=1,nlev
         xno2 (i,j,k) = x11(k)
      enddo
c
c
c 5. BC:
c
      do k=1,nlev
        x11(k) =  bcarbon  (i,j,k)
      enddo
       call int_accri_emi(j,ps4accri(i,j),BCaccri,emi_accri,
     &                     mymonth,myyear)

      call chememissionA(1,j,nlev,
     &                   emi_accri,x11)
      do k=1,nlev
         bcarbon (i,j,k) = x11(k)
      enddo
c
c
c 6. OC:
c
      do k=1,nlev
        x11(k) =  ocarbon  (i,j,k)
      enddo
       call int_accri_emi(j,ps4accri(i,j),OCaccri,emi_accri,
     &                     mymonth,myyear)

      call chememissionA(1,j,nlev,
     &                   emi_accri,x11)
      do k=1,nlev
         ocarbon (i,j,k) = x11(k)
      enddo
c
c
c 7. SO2:
c
      do k=1,nlev
        x11(k) =  so2  (i,j,k)
      enddo
       call int_accri_emi(j,ps4accri(i,j),SO2accri,emi_accri,
     &                     mymonth,myyear)

      call chememissionA(1,j,nlev,
     &                   emi_accri,x11)
      do k=1,nlev
         so2 (i,j,k) = x11(k)
      enddo
c
c
c 8. CH2O:
c
      do k=1,nlev
        x11(k) =  ch2o  (i,j,k)
      enddo
       call int_accri_emi(j,ps4accri(i,j),CH2Oaccri,emi_accri,
     &                     mymonth,myyear)

      call chememissionA(1,j,nlev,
     &                   emi_accri,x11)
      do k=1,nlev
         ch2o (i,j,k) = x11(k)
      enddo
c
c 9. H2O:
c
!     print *,' ACCRI H2o'
      do k=1,nlev
        x11(k) =  Q  (i,j,k)
      enddo
       call int_accri_emi(j,ps4accri(i,j),H2Oaccri,emi_accri,
     &                     mymonth,myyear)

      call chememissionA(1,j,nlev,
     &                   emi_accri,x11)
      do k=1,nlev
         Q (i,j,k) = x11(k)
      enddo
c

#endif

1     continue

#endif

      return
       end

#if ( defined ACCRI )
c=====================================================

      subroutine chememissionA(i,j,ktop,
     &                         emi_accri,x11)
c     ===========================================

c======================================================c
c  A subroutine for adding ACCRI emission               c
c     into layers from k=1 to ktop (nlev)    c
c======================================================c

#include "chem_para"
#include "chem_com"

      dimension x11(ktop)
      dimension emi_accri(ktop)

#if ( defined CPL_CHEM )
!     print *,'  '
!       xmass    	!Total tracer mass in 10^-9 kg
      do k=1,ktop
        xmass    = 
!    &             airmass(i,j,k)*x11(k)
     &             x11(k)
        xemi     = ( xmass 
!    &           + emi_accri(k) ) / airmass(i,j,k)	! ppbm
     &           + emi_accri(k) / airmass(i,j,k) )	! ppbm

!     print *,(x11(k)-xemi),emi_accri(k)

        x11(k)  = xemi
      enddo
!     print *,'  '
#endif

      return
       end

#endif
c=====================================================

      subroutine chememission2(i,j,ktop,eratio,tmass1,
     &                         yamount,x11)
c     ===========================================

c======================================================c
c  A subroutine for calculating emission               c
c    by mixing tracers into from k=1 to ktop layers    c
c======================================================c

#include "chem_para"
#include "chem_com"

      dimension x11(ktop+1)

#if ( defined CPL_CHEM )
        xmass    = 0.0	!Total tracer mass in 10^-9 kg
      do k=1,ktop
        xmass    = xmass
     &           + airmass(i,j,k)*x11(k)
      enddo
        xemi     = ( xmass 
     &           + yamount ) * tmass1	! ppbm

      do k=1,ktop
        x11(k)  = xemi
      enddo
#endif

      return
       end

c=====================================================

      subroutine chememission3(i,j,ktop,eratio,tmass1,
     &                         yamount,x11)
c     ===========================================

c======================================================c
c  A subroutine for calculating emission               c
c    by mixing tracers into one extra layer than       c
c    chemiemission2                                    c
c======================================================c

#include "chem_para"
#include "chem_com"

      dimension x11(ktop+1)

#if ( defined CPL_CHEM )

      kup = ktop + 1  

        xmass    = 0.0                !Total tracer mass in kg
      do k=1,ktop
        xmass    = xmass
     &           + airmass(i,j,k)*x11(k)
      enddo
        xmass    = xmass
     &           + airmass(i,j,kup)*x11(kup)
     &            *eratio
     &           + yamount 
        xemi     = xmass * tmass1     !ppbm

      do k=1,ktop
        x11(k)  = xemi
      enddo
        x11(kup)= max(0.0,x11(kup) 
     &          + eratio*(xemi - x11(kup)))

#endif

      return
       end

c======================================================
      subroutine chememiltn(i,j,ktop,
     &                         xfunc,yamount,x11)
c     ===========================================

c======================================================c
c  A subroutine for calculating lightning-produced     c
c    NO from k=1 to ktop                               c
c======================================================c

#include "chem_para"
#include "chem_com"

      dimension x11(ktop)

#if ( defined CPL_CHEM )

c first 4 layers share 70% total production:

      tshare = yamount*xfunc
      xshare = 0.175*tshare		!0.7/4

      do k=1,4
        x11(k) = x11(k) + xshare/airmass(i,j,k)
      enddo

c then the rest 3 layers share 30%

      xshare = 0.1*tshare

      do k=5,n_tropopause
        x11(k) = x11(k) + xshare/airmass(i,j,k)
      enddo

#endif

      return
       end


c===============================================================

      Block data xemission

c==============================================================c
c     Global emission amount of tracers in 10^6 kg and         c
c            start from 1931 current end in 1991.              c
c==============================================================c

#include "chem_para"
#include "chem_com"

#if ( defined CPL_CHEM )

c
c percentage of air mass of the third layer assumed
c   to be involved in emission calculation:
c   April 10, 1995:
c

#if ( N_LAT == 24 )
      data ehpbl/
     &     0.000, 0.000, 0.000, 0.000, 0.000,
     &     0.000, 0.000, 0.000, 0.096, 0.221,    
     &     0.059, 0.203, 0.122, 0.199, 0.606,    
     &     1.000, 1.000, 1.000, 1.000, 1.000,    
     &     1.000, 0.531, 0.000, 0.000
     &					/
#endif
#if ( N_LAT == 46 )
      data    ehpbl/
     &      0.000, 0.000, 0.000, 0.000, 0.000,
     &      0.000, 0.000, 0.000, 0.000, 0.000,
     &      0.000, 0.000, 0.000, 0.000, 0.015,
     &      0.064, 0.118, 0.182, 0.189, 0.106,
     &      0.091, 0.165, 0.183, 0.142, 0.143,
     &      0.182, 0.317, 0.525, 0.729, 0.930,
     &      1.000, 1.000, 1.000, 1.000, 1.000,
     &      1.000, 1.000, 1.000, 1.000, 1.000,
     &      0.792, 0.552, 0.283, 0.012, 0.000,
     &      0.000
     &     /
#endif

c---
c Lightning produced NO emission data, 062795:
c   based on Kumar et al., 1995:
c
! --- 121800: moved initialization to cheminit.F
!       data xnoxltnt/
c     &      61*1.072e+4
!     &      124*1.072e+4
!     &                 /      
       data xnoxltnm/
     &      0.0776,  0.0776,  0.0866,  0.0866,  0.0866,  0.0697,
     &      0.0697,  0.0697,  0.0995,  0.0995,  0.0995,  0.0776/

#if ( N_LAT == 24 )
       data xnoxltnd/
     &      0.0000,  0.0000,  0.0000,  0.0025,  0.0124,  0.0248,
     &      0.0413,  0.0661,  0.1321,  0.1131,  0.1024,  0.0941,
     &      0.0983,  0.0661,  0.0520,  0.0537,  0.0528,  0.0413,
     &      0.0256,  0.0165,  0.0050,  0.0000,  0.0000,  0.0000,

     &      0.0000,  0.0000,  0.0000,  0.0025,  0.0124,  0.0248,
     &      0.0413,  0.0661,  0.1321,  0.1131,  0.1024,  0.0941,
     &      0.0983,  0.0661,  0.0520,  0.0537,  0.0528,  0.0413,
     &      0.0256,  0.0165,  0.0050,  0.0000,  0.0000,  0.0000,

     &      0.0000,  0.0000,  0.0000,  0.0089,  0.0192,  0.0296,
     &      0.0400,  0.0570,  0.0725,  0.0814,  0.0918,  0.1066,
     &      0.0807,  0.0777,  0.0777,  0.0718,  0.0563,  0.0666,
     &      0.0370,  0.0148,  0.0104,  0.0000,  0.0000,  0.0000,

     &      0.0000,  0.0000,  0.0000,  0.0089,  0.0192,  0.0296,
     &      0.0400,  0.0570,  0.0725,  0.0814,  0.0918,  0.1066,
     &      0.0807,  0.0777,  0.0777,  0.0718,  0.0563,  0.0666,
     &      0.0370,  0.0148,  0.0104,  0.0000,  0.0000,  0.0000,

     &      0.0000,  0.0000,  0.0000,  0.0089,  0.0192,  0.0296,
     &      0.0400,  0.0570,  0.0725,  0.0814,  0.0918,  0.1066,
     &      0.0807,  0.0777,  0.0777,  0.0718,  0.0563,  0.0666,
     &      0.0370,  0.0148,  0.0104,  0.0000,  0.0000,  0.0000,

     &      0.0000,  0.0000,  0.0000,  0.0028,  0.0110,  0.0221,
     &      0.0331,  0.0414,  0.0423,  0.0643,  0.0855,  0.0827,
     &      0.1048,  0.1232,  0.1232,  0.1048,  0.0763,  0.0404,
     &      0.0239,  0.0129,  0.0055,  0.0000,  0.0000,  0.0000,

     &      0.0000,  0.0000,  0.0000,  0.0028,  0.0110,  0.0221,
     &      0.0331,  0.0414,  0.0423,  0.0643,  0.0855,  0.0827,
     &      0.1048,  0.1232,  0.1232,  0.1048,  0.0763,  0.0404,
     &      0.0239,  0.0129,  0.0055,  0.0000,  0.0000,  0.0000,

     &      0.0000,  0.0000,  0.0000,  0.0028,  0.0110,  0.0221,
     &      0.0331,  0.0414,  0.0423,  0.0643,  0.0855,  0.0827,
     &      0.1048,  0.1232,  0.1232,  0.1048,  0.0763,  0.0404,
     &      0.0239,  0.0129,  0.0055,  0.0000,  0.0000,  0.0000,

     &      0.0000,  0.0000,  0.0000,  0.0077,  0.0071,  0.0180,
     &      0.0258,  0.0444,  0.0657,  0.0811,  0.0927,  0.1037,
     &      0.1050,  0.1017,  0.0985,  0.0882,  0.0670,  0.0476,
     &      0.0296,  0.0122,  0.0039,  0.0000,  0.0000,  0.0000,

     &      0.0000,  0.0000,  0.0000,  0.0077,  0.0071,  0.0180,
     &      0.0258,  0.0444,  0.0657,  0.0811,  0.0927,  0.1037,
     &      0.1050,  0.1017,  0.0985,  0.0882,  0.0670,  0.0476,
     &      0.0296,  0.0122,  0.0039,  0.0000,  0.0000,  0.0000,

     &      0.0000,  0.0000,  0.0000,  0.0077,  0.0071,  0.0180,
     &      0.0258,  0.0444,  0.0657,  0.0811,  0.0927,  0.1037,
     &      0.1050,  0.1017,  0.0985,  0.0882,  0.0670,  0.0476,
     &      0.0296,  0.0122,  0.0039,  0.0000,  0.0000,  0.0000,

     &      0.0000,  0.0000,  0.0000,  0.0025,  0.0124,  0.0248,
     &      0.0413,  0.0661,  0.1321,  0.1131,  0.1024,  0.0941,
     &      0.0983,  0.0661,  0.0520,  0.0537,  0.0528,  0.0413,
     &      0.0256,  0.0165,  0.0050,  0.0000,  0.0000,  0.0000/

#endif
#if ( N_LAT == 46 )
      data xnoxltnd/
     &      0.0000, 0.0000, 0.0000, 0.0000, 0.0001, 0.0014,
     &      0.0032, 0.0082, 0.0135, 0.0198, 0.0266, 0.0351,
     &      0.0446, 0.0573, 0.0764, 0.1101, 0.1287, 0.1190,
     &      0.1110, 0.1055, 0.1006, 0.0963, 0.0951, 0.0973,
     &      0.0897, 0.0733, 0.0620, 0.0548, 0.0525, 0.0534,
     &      0.0534, 0.0529, 0.0487, 0.0428, 0.0354, 0.0273,
     &      0.0220, 0.0173, 0.0116, 0.0058, 0.0028, 0.0002,
     &      0.0000, 0.0000, 0.0000, 0.0000,

     &      0.0000, 0.0000, 0.0000, 0.0000, 0.0001, 0.0014,
     &      0.0032, 0.0082, 0.0135, 0.0198, 0.0266, 0.0351,
     &      0.0446, 0.0573, 0.0764, 0.1101, 0.1287, 0.1190,
     &      0.1110, 0.1055, 0.1006, 0.0963, 0.0951, 0.0973,
     &      0.0897, 0.0733, 0.0620, 0.0548, 0.0525, 0.0534,
     &      0.0534, 0.0529, 0.0487, 0.0428, 0.0354, 0.0273,
     &      0.0220, 0.0173, 0.0116, 0.0058, 0.0028, 0.0002,
     &      0.0000, 0.0000, 0.0000, 0.0000,

     &      0.0000, 0.0000, 0.0000, 0.0000, 0.0004, 0.0049,
     &      0.0096, 0.0149, 0.0201, 0.0254, 0.0308, 0.0361,
     &      0.0423, 0.0510, 0.0594, 0.0673, 0.0741, 0.0786,
     &      0.0835, 0.0888, 0.0951, 0.1027, 0.1003, 0.0870,
     &      0.0799, 0.0784, 0.0777, 0.0777, 0.0759, 0.0728,
     &      0.0666, 0.0587, 0.0600, 0.0652, 0.0554, 0.0403,
     &      0.0281, 0.0168, 0.0129, 0.0107, 0.0058, 0.0005,
     &      0.0000, 0.0000, 0.0000, 0.0000,

     &      0.0000, 0.0000, 0.0000, 0.0000, 0.0004, 0.0049,
     &      0.0096, 0.0149, 0.0201, 0.0254, 0.0308, 0.0361,
     &      0.0423, 0.0510, 0.0594, 0.0673, 0.0741, 0.0786,
     &      0.0835, 0.0888, 0.0951, 0.1027, 0.1003, 0.0870,
     &      0.0799, 0.0784, 0.0777, 0.0777, 0.0759, 0.0728,
     &      0.0666, 0.0587, 0.0600, 0.0652, 0.0554, 0.0403,
     &      0.0281, 0.0168, 0.0129, 0.0107, 0.0058, 0.0005,
     &      0.0000, 0.0000, 0.0000, 0.0000,

     &      0.0000, 0.0000, 0.0000, 0.0000, 0.0004, 0.0049,
     &      0.0096, 0.0149, 0.0201, 0.0254, 0.0308, 0.0361,
     &      0.0423, 0.0510, 0.0594, 0.0673, 0.0741, 0.0786,
     &      0.0835, 0.0888, 0.0951, 0.1027, 0.1003, 0.0870,
     &      0.0799, 0.0784, 0.0777, 0.0777, 0.0759, 0.0728,
     &      0.0666, 0.0587, 0.0600, 0.0652, 0.0554, 0.0403,
     &      0.0281, 0.0168, 0.0129, 0.0107, 0.0058, 0.0005,
     &      0.0000, 0.0000, 0.0000, 0.0000,

     &      0.0000, 0.0000, 0.0000, 0.0000, 0.0001, 0.0016,
     &      0.0033, 0.0075, 0.0120, 0.0177, 0.0233, 0.0289,
     &      0.0342, 0.0384, 0.0415, 0.0420, 0.0462, 0.0575,
     &      0.0685, 0.0794, 0.0849, 0.0834, 0.0881, 0.0994,
     &      0.1097, 0.1191, 0.1232, 0.1232, 0.1175, 0.1081,
     &      0.0953, 0.0807, 0.0635, 0.0452, 0.0342, 0.0257,
     &      0.0195, 0.0139, 0.0098, 0.0060, 0.0031, 0.0002,
     &      0.0000, 0.0000, 0.0000, 0.0000,

     &      0.0000, 0.0000, 0.0000, 0.0000, 0.0001, 0.0016,
     &      0.0033, 0.0075, 0.0120, 0.0177, 0.0233, 0.0289,
     &      0.0342, 0.0384, 0.0415, 0.0420, 0.0462, 0.0575,
     &      0.0685, 0.0794, 0.0849, 0.0834, 0.0881, 0.0994,
     &      0.1097, 0.1191, 0.1232, 0.1232, 0.1175, 0.1081,
     &      0.0953, 0.0807, 0.0635, 0.0452, 0.0342, 0.0257,
     &      0.0195, 0.0139, 0.0098, 0.0060, 0.0031, 0.0002,
     &      0.0000, 0.0000, 0.0000, 0.0000,

     &      0.0000, 0.0000, 0.0000, 0.0000, 0.0001, 0.0016,
     &      0.0033, 0.0075, 0.0120, 0.0177, 0.0233, 0.0289,
     &      0.0342, 0.0384, 0.0415, 0.0420, 0.0462, 0.0575,
     &      0.0685, 0.0794, 0.0849, 0.0834, 0.0881, 0.0994,
     &      0.1097, 0.1191, 0.1232, 0.1232, 0.1175, 0.1081,
     &      0.0953, 0.0807, 0.0635, 0.0452, 0.0342, 0.0257,
     &      0.0195, 0.0139, 0.0098, 0.0060, 0.0031, 0.0002,
     &      0.0000, 0.0000, 0.0000, 0.0000,

     &      0.0000, 0.0000, 0.0000, 0.0000, 0.0003, 0.0043,
     &      0.0077, 0.0074, 0.0081, 0.0136, 0.0189, 0.0229,
     &      0.0283, 0.0378, 0.0477, 0.0586, 0.0684, 0.0763,
     &      0.0834, 0.0893, 0.0951, 0.1008, 0.1040, 0.1047,
     &      0.1041, 0.1024, 0.1008, 0.0991, 0.0953, 0.0900,
     &      0.0811, 0.0703, 0.0601, 0.0502, 0.0408, 0.0316,
     &      0.0226, 0.0137, 0.0087, 0.0045, 0.0022, 0.0002,
     &      0.0000, 0.0000, 0.0000, 0.0000,

     &      0.0000, 0.0000, 0.0000, 0.0000, 0.0003, 0.0043,
     &      0.0077, 0.0074, 0.0081, 0.0136, 0.0189, 0.0229,
     &      0.0283, 0.0378, 0.0477, 0.0586, 0.0684, 0.0763,
     &      0.0834, 0.0893, 0.0951, 0.1008, 0.1040, 0.1047,
     &      0.1041, 0.1024, 0.1008, 0.0991, 0.0953, 0.0900,
     &      0.0811, 0.0703, 0.0601, 0.0502, 0.0408, 0.0316,
     &      0.0226, 0.0137, 0.0087, 0.0045, 0.0022, 0.0002,
     &      0.0000, 0.0000, 0.0000, 0.0000,

     &      0.0000, 0.0000, 0.0000, 0.0000, 0.0003, 0.0043,
     &      0.0077, 0.0074, 0.0081, 0.0136, 0.0189, 0.0229,
     &      0.0283, 0.0378, 0.0477, 0.0586, 0.0684, 0.0763,
     &      0.0834, 0.0893, 0.0951, 0.1008, 0.1040, 0.1047,
     &      0.1041, 0.1024, 0.1008, 0.0991, 0.0953, 0.0900,
     &      0.0811, 0.0703, 0.0601, 0.0502, 0.0408, 0.0316,
     &      0.0226, 0.0137, 0.0087, 0.0045, 0.0022, 0.0002,
     &      0.0000, 0.0000, 0.0000, 0.0000,

     &      0.0000, 0.0000, 0.0000, 0.0000, 0.0001, 0.0014,
     &      0.0032, 0.0082, 0.0135, 0.0198, 0.0266, 0.0351,
     &      0.0446, 0.0573, 0.0764, 0.1101, 0.1287, 0.1190,
     &      0.1110, 0.1055, 0.1006, 0.0963, 0.0951, 0.0973,
     &      0.0897, 0.0733, 0.0620, 0.0548, 0.0525, 0.0534,
     &      0.0534, 0.0529, 0.0487, 0.0428, 0.0354, 0.0273,
     &      0.0220, 0.0173, 0.0116, 0.0058, 0.0028, 0.0002,
     &      0.0000, 0.0000, 0.0000, 0.0000
     &     /
#endif
#endif

      end
