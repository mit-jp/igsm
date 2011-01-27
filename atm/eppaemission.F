
#include "ctrparam.h"

        subroutine eppaemission (iyr)
!       ===================================
#if ( defined CPL_CHEM )
#include "chem_para"
#include "chem_com"
#include "BD2G04.COM"
#include "chem_meta"
      parameter (neppa=124)
      character * 120 emiss_data,SO2ERATIO
      namelist /EPPA/ LYEAREM,emiss_data,SO2ERATIO
      common /SO2EMIN/SO2EM
      data ifirst / 0 /
      common /EPPAEMISS/  edailyf11eppa     (nlon,nlat,neppa),
     &           edailyf12eppa     (nlon,nlat,neppa),
     &           edailyn2oeppa     (nlon,nlat,neppa),
     &           edailycoeppa      (nlon,nlat,neppa),
     &           edailynoxeppa     (nlon,nlat,neppa),
     &           edailych4eppa     (nlon,nlat,neppa),
     &           edailyso2eppa     (nlon,nlat,neppa),
     &           edailyco2eppa     (nlon,nlat,neppa),
     &           edailyhfc134aeppa (nlon,nlat,neppa),
     &           edailypfceppa     (nlon,nlat,neppa),
     &           edailysf6eppa     (nlon,nlat,neppa),
     &           edailyucoeppa     (nlon,nlat,neppa),
     &           edailyunmveppa    (nlon,nlat,neppa),
     &           edailyunoxeppa    (nlon,nlat,neppa),
     &           edailyusoxeppa    (nlon,nlat,neppa),
     &           edailybceppa      (nlon,nlat,neppa),
     &           edailynh3eppa     (nlon,nlat,neppa),
     &           edailyoceppa      (nlon,nlat,neppa),
     &           edailyubceppa     (nlon,nlat,neppa),
     &           edailyunh3eppa    (nlon,nlat,neppa),
     &           edailyuoceppa     (nlon,nlat,neppa),
     &           n_total_urbaneppa (nlat,neppa),
     &           n_urbaneppa       (3,nlat,neppa)

      if (ifirst.eq.0) then
      open(535,file='eppaemis.dat')
      read (535,NML=EPPA)
      close (535)
      open(164,file=emiss_data,
     &        form='unformatted',
     &        status='old')
      open(664,file=SO2ERATIO,
     &        form='formatted',
     &        status='old')

      read(664,'(f10.6)')SO2EM
       print *,' SO2EM=',SO2EM

      read(164)edailyf11eppa,
     &         edailyf12eppa,
     &         edailyn2oeppa,
     &         edailycoeppa,
     &         edailynoxeppa,
     &         edailych4eppa,
     &         edailyso2eppa,
     &         edailyco2eppa,
     &         edailyhfc134aeppa,
     &         edailypfceppa,
     &         edailysf6eppa,
     &         edailyucoeppa,
     &         edailyunmveppa,
     &         edailyunoxeppa,
     &         edailyusoxeppa,
     &         edailybceppa,
     &         edailynh3eppa,
     &         edailyoceppa,
     &         edailyubceppa,
     &         edailyunh3eppa,
     &         edailyuoceppa,
     &         n_total_urbaneppa,
     &         n_urbaneppa

        close (164)
        close (664)
        print *,'Emissioms for ',neppa,' years'
        ifirst = 1
	return
	endif 	! first
        iyremi  = iyr-1976       !year from starting point
        myyearlast  = min(LYEAREM-1976,neppa)  !last year of emission
        print *,' EPPAEMISION myyearlast=',myyearlast
        iyremi = min(iyremi,myyearlast)
        print *,' EPPAEMISION iyr=',iyr,' iyremi=',iyremi
        print *,' EPPAEMISION nlat=',nlat
         print *,'CO2 emissions for ',iyremi
         print *,(edailyco2eppa(1,j,iyremi),j=1,nlat)
        do j=1,nlat
         edailyf11(1,j,1)=edailyf11eppa(1,j,iyremi)
         edailyf12(1,j,1)=edailyf12eppa(1,j,iyremi)
         edailyn2o(1,j,1)=edailyn2oeppa(1,j,iyremi)
         edailyco(1,j,1)=edailycoeppa(1,j,iyremi)
         edailynox(1,j,1)=edailynoxeppa(1,j,iyremi)
         edailych4(1,j,1)=edailych4eppa(1,j,iyremi)
         edailyso2(1,j,1)=edailyso2eppa(1,j,iyremi)
         edailyco2(1,j,1)=edailyco2eppa(1,j,iyremi)
         edailyhfc134a(1,j,1)=edailyhfc134aeppa(1,j,iyremi)
         edailypfc(1,j,1)=edailypfceppa(1,j,iyremi)
         edailysf6(1,j,1)=edailysf6eppa(1,j,iyremi)
         edailyuco(1,j,1)=edailyucoeppa(1,j,iyremi)
         edailyunmv(1,j,1)=edailyunmveppa(1,j,iyremi)
         edailyunox(1,j,1)=edailyunoxeppa(1,j,iyremi)
         edailyusox(1,j,1)=edailyusoxeppa(1,j,iyremi)
         edailybc(1,j,1)=edailybceppa(1,j,iyremi)
         edailynh3(1,j,1)=edailynh3eppa(1,j,iyremi)
         edailyoc(1,j,1)=edailyoceppa(1,j,iyremi)
         edailyubc(1,j,1)=edailyubceppa(1,j,iyremi)
         edailyunh3(1,j,1)=edailyunh3eppa(1,j,iyremi)
         edailyuoc(1,j,1)=edailyuoceppa(1,j,iyremi)
         n_total_urban (j,1)=n_total_urbaneppa(j,iyremi)
         do k=1,3
            n_urban (k,j,1)= n_urbaneppa(k,j,iyremi)
         enddo
       enddo
   
#endif
	return
	end
