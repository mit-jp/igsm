#include "ctrparam.h"

        subroutine eppaemission (iyr)
!       ===================================
#if ( defined CPL_CHEM )
#include "chem_para"
#include "chem_com"
#include "BD2G04.COM"
#include "chem_meta"
#if ( defined OLD_EMISS )
      parameter (neppa=124,naccri=110)
!     parameter (neppa=324,naccri=110)
!     parameter (neppa=524,naccri=110)
#else
      parameter (neppa=1,naccri=110)
#endif
#endif

      logical FORMSO2
      character * 120 emiss_data,SO2ERATIO
#if ( defined ACCRI )
     &  ,accri_data,accri_scaling
#endif
      namelist /EPPA/ LYEAREM,emiss_data,SO2ERATIO,FORMSO2
     & ,IYRSO2EMI
#if ( defined ACCRI )
     &  ,accri_data,accri_scaling
#endif
      common /SO2EMIN/IYRSO2EMI,SO2EM
      data ifirst / 0 /
#if ( defined CPL_CHEM )
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
#if ( defined ACCRI )
      real scale_accri_eppa(neppa)
#endif
#endif

      if (ifirst.eq.0) then
!     FORMSO2=.TRUE.
      open(535,file='eppaemis.dat')
       print *,' CALL EPPAEMISSION iyr=',iyr
      IYRSO2EMI=1997
      read (535,NML=EPPA)
      print *,' FORMSO2=',FORMSO2
      close (535)
#if ( defined CPL_CHEM )
      open(164,file=emiss_data,
     &        form='unformatted',
     &        status='old')
#endif
      if(FORMSO2) then
      open(664,file=SO2ERATIO,
     &        form='formatted',
     &        status='old')
       read(664,'(f10.6)')SO2EM
      else
      open(664,file=SO2ERATIO,
     &        form='unformatted',
     &        status='old')
       read(664)SO2EM
      endif
       print *,' SO2EM=from eppaemissions ',SO2EM
        close (664)
#if ( defined CPL_CHEM )
#if ( defined OLD_EMISS )
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
        print *,'Emissioms for ',neppa,' years'
#else
      read(164)iyemif,iyemil
      print *,'EPPA emissions for years ',iyemif,iyemil
      if (iyemil.lt.LYEAREM) then
       print *,'Wrong iyemil and LYEAREM'
       print *,'iyemil=',iyemil,' LYEAREM=',LYEAREM
       stop
      endif
      if (iyemif.gt.iyr) then
       print *,'Wrong iyemif and iyr'
       print *,'iyemif=',iyemif,' iyr=',iyr
       stop
      else
       iyr1=min(iyr-1,LYEAREM)
       do iy=iyemif,iyr1
        read(164)iyemi
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
        print *,'iyemi=',iyemi
       enddo
      endif

#endif
#if ( defined ACCRI )
      print *,' Run with ACCRI emissions with H2O'
      print *,accri_data
      open(264,file=accri_data,
     &        form='unformatted',
     &        status='old')

      open(265,file=accri_scaling,
     &        form='formatted',
     &        status='old')

      read(264)p_accri
      print *,'p_accri'
      print *,p_accri
#endif
#endif

        ifirst = 1
	return
	endif 	! first
#if ( defined CPL_CHEM )
#if ( defined OLD_EMISS )
        iyremi  = iyr-1976       !year from starting point
        myyearlast  = min(LYEAREM-1976,neppa)  !last year of emission
        print *,' EPPAEMISION myyearlast=',myyearlast
        iyremi = min(iyremi,myyearlast)
        print *,' EPPAEMISION iyr=',iyr,' iyremi=',iyremi
        print *,' EPPAEMISION nlat=',nlat
         print *,'CO2 emissions for ',iyremi
         print *,(edailyco2eppa(1,j,iyremi),j=1,nlat)
#else
        iyremi  = 1
      if(iyemi.lt.iyemil.and.iyemi.lt.LYEAREM) then
        read(164)iyemi
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
         endif

         print *,iyemi,' emissions used for ',iyr
         print *,'CO2 emissions for ',iyemi
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
#if ( defined ACCRI )
      print *,'Reading ACCRI emissions for year ',iyr
      read(264)COaccri,
     &        NOaccri,
     &        NO2accri,
     &        BCaccri,
     &        OCaccri,
     &        CO2accri,
     &        SO2accri,
!    &        CH2Oaccri
     &        CH2Oaccri,
     &        H2Oaccri
      print*,'COaccri=',COaccri(1,23,1,1),
     &        'NOaccri=',NOaccri(1,23,1,1),
     &        'NO2accri=',NO2accri(1,23,1,1),
     &        'BCaccri=',BCaccri(1,23,1,1),
     &        'OCaccri=',OCaccri(1,23,1,1),
     &        'CO2accri=',CO2accri(1,23,1,1),
     &        'SO2accri=',SO2accri(1,23,1,1),
     &        'CH2Oaccri=',CH2Oaccri(1,23,1,1),
     &        'H2Oaccri=',H2Oaccri(1,23,1,1)
!       print *,'H2Oaccri'
!       print *,(H2Oaccri(1,23,1,nm),nm=1,12)
         scale_accri(1)=1.0
#endif
#endif
   
	return
	end
