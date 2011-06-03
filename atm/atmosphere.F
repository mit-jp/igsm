
#include "ctrparam.h"

!	==========================================================
!
!	Atmosphere.F:	 Former main program of the MIT Global Climate and
!			Biogeochemistry Model.
!
!	----------------------------------------------------------
!
!	Repacking Note: This version is combined with main.f of  
!		several originally separated model versions such as 
!		CliChem 3.0, CliChemNem, as well as MODEL24x11.
!
!		The chemistry module, ocean CO2 module,
!		and TEM module are all controlled by cpp now.
!
!		Chien Wang
!		MIT Joint Program for the Science and Policy
!			of Global Change
!
!	----------------------------------------------------------
! 
!	Revision History:
!	
!	When	Who		What
!	----	----------	-------	
!	073100	Chien Wang	Created from CliChem3.0 and MODEL24x11
!       100500  Andrei          Converted to subroutine. When called first
!                               time does only initialization. 
!                               Creates monthly data for TEM
!
!	==========================================================


      SUBROUTINE ATMOSPHERE(DTATM,mndriver)
C**** MD2G04 BD2G04 MD2G04                01/02/93                         0.1  
C**** OPT(3)                                                               0.2  
C****                                                                      0.3  
C**** Dynamics and physics programs for 2-D model.                         0.4  
C**** Like D2G04 but run on work station.                                  0.5 
 
#include "BD2G04.COM"					     1.   
#include "ODIFF.COM"
#include "run.COM"					     1.   

#include "chem_para"
#include "chem_com"

#if ( defined CPL_CHEM )

#include "chem_tmp"
	integer hrcnt, cnt3hr(8)	! for ozone impact
	real  sfc3hro3(nlat,8)
creal*4  sfc3hro3(nlat,8)

#endif

      COMMON/INTA/COE1(01,01,01),COE2(01,01,01)                            1.1  
c     DIMENSION ACO(36,01),BSI(36,01),CCO(36,01),DSI(36,01)                1.2  
      COMMON/SPEC1/                                                        1.3  
     *  XA(IM0,JM0),XB(IM0,JM0),YA(IM0,JM0),YB(IM0,JM0),ZA(IM0,JM0)
     &,ZB(IM0,JM0)                                                         1.4  
      COMMON/SPEC2/KM,KINC,COEK,C3LAND(IO0,JM0),C3OICE(IO0,JM0)            1.5  
     *  ,C3LICE(IO0,JM0),WMGE(IO0,JM0),TSSFC(1,JM0,4)                      1.6  
      COMMON/EPARA/VTH(JM0,LM0),WTH(JM0,LM0),VU(JM0,LM0),VV(JM0,LM0)
     &,DQSDT(JM0,LM0)                                                      1.7  
     *  ,DWV(JM0),PHIT(JM0,LM0),TPRIM2(JM0,LM0),WU(JM0,LM0),CKS,CKN        1.8  
     *  ,WQ(JM0,LM0),VQ(JM0,LM0),MRCHT                                     1.9  
      COMMON U,V,T,P,Q                                                     2.   
      COMMON/WORK1/WORKX(IM0,JM0,LM0),UT(IM0,JM0,LM0),VT(IM0,JM0,LM0),     3.   
     *  TT(IM0,JM0,LM0),PT(IM0,JM0),QT(IM0,JM0,LM0)                        4.   
      COMMON/WORK2/UX(IM0,JM0,LM0),VX(IM0,JM0,LM0),TX(IM0,JM0,LM0)
     &,PX(IM0,JM0)                                                         5.   
      COMMON/OLDZO/ZMLOLD(IO0,JM0)
C     COMMON/KEYS/KEYNR(42,50)                                             8.   
      CHARACTER C*4,CYEAR*4,CMMND*80                                       8.1  
      DIMENSION C(39),JC(100),RC(161)                                      8.2  
      EQUIVALENCE (JC(1),IM),(C(1),XLABEL(1),LABEL1),(RC(1),TAU)           8.3  
      CHARACTER*8 LABSSW,LABEL1,OFFSSW/'XXXXXXXX'/                         8.4  
      LOGICAL EVENT,wr25,HPRNT,TRANSR,LHORDIF                              9.   
     & ,CONTRR,OBSFOR,FIRST,NOCLM
      common/conprn/HPRNT,JPR,LPR
      common/wrcom/wr25,TRANSR,CONTRR,OBSFOR
      COMMON/CWMG/WMGEA(JM0),NWMGEA(JM0),RIGA(JM0),DTAV(JM0),DQAV(JM0),
     *Z0AV(JM0),WSAV(JM0),WS0AV(JM0),TAUAV(JM0)
      COMMON/OCN/TG3M(1,JM0,12),RTGO(1,JM0,lmo),STG3(1,JM0),DTG3(1,JM0)
      common/SURRAD/TRSURF(JM0,4),SRSURF(JM0,4)
      dimension RTGOAV(JM0,lmo)
      common/tprmtg/tprmg(JM0),ntprmg(JM0)
      common/aexpc/AEXP,ISTRT1,ISTRTCHEM
      common/mixlr/Z1OAV(JM0),NZ1OAV(JM0)
      common/flxio/FLIO(JM0),NFLIO(JM0)
      common/surps/srps(JM0+3),nsrps
c     character *19 buf
c     character *8 buf1
      character *120 file1,file2,plotfl,nwrfl,qffile,clfile
      CHARACTER*4 AMONTH(12),JMONTHPR
      common/files/file1,file2,flotfl,nwrfl,qffile,clfile
      common/PRNT1/NCOMP
      common/Dscale/DWAV0(JM0)
      COMMON/CO2TRND/ALFFOR,CO2TR
      COMMON/FRMIC/ FRMDICE(JM0)
      common/ BACKGRGHG/GHGBGR(5)
      DATA AMONTH/'JAN','FEB','MAR','APR','MAY','JUNE','JULY'
     &  ,'AUG','SEP','OCT','NOV','DEC'/
      dimension fluxnep(jm0)
      logical odifcarbon

#if ( defined CLM )
#  include "CLM.h"
#  if ( defined CPL_TEM )      
C  For TEM
#    include "TEM.h"
#  endif
#endif

#if ( defined  CPL_OCEANCO2 )
#include "OCM.h"
      common /Garyflux/pC_atm(jm0),wind_amp,fluxco2(jm0)
#  if ( defined ML_2D )
      common/Garyclim/tggary(jm0),wsgary(jm0),areaml(jm0),arsrf(jm0)
      common/Garydiff/depthml(jm0),edzcar(jm0),dzg(lmo),dzog(lmo-1),
     &Rco2(jm0,lmo),edohd(lmo),zg(lmo),focean(jm0)
      common /Garychem/Hg(jm0)
      common/qfl/QFLUX(JM0,0:13),ZOAV(JM0),QFLUXT(JM0)
      common /Garyvdif/iyearocm,vdfocm,acvdfc
      common /Garyvlog/odifcarbon,ocarcont
      common /Garykvct/cfkvct,edzcart(jm0)
# endif
#endif

      INTEGER dtatm, mndriver   !routine arguments jrs

#if ( defined OCEAN_3D || defined ML_2D )
#include "AGRID.h"
        dimension oimeltt(jm0),dhdtav(jm0),devdtav(jm0)
#endif

C **** CLEAR SKY
      common/clrsk/CLEAR(JM0),NCLR(JM0),AJCLR(JM0,12),BJCLR(JM0,12),
     * CJCLR(JM0,12)
!     common/TSUR/TSURFC(JM0,0:13),TSURFT(JM0),TSURFD(JM0),DTSURF(JM0)
#include "TSRF.COM"
      real TSURFW(JM0),TLANDW(JM0)
      integer CLEAR

      common /ATCO2/atm_co2(jm0)

#if ( defined CPL_NEM )      
C     For Emission
c === 031097
      real  ECH4COR(JM0),ECH4PY(JM0,12), ECH4OUT(JM0),EPJT(JM0)
     &,ECH4CTR(JM0)
      real  EN2OCOR(JM0),EN2OPY(JM0,12), EN2OOUT(JM0),EPJTN2O(JM0)
     &,EN2OCTR(JM0)
c     ECH4CHIEN and EN2OCHIEN are used in chemistry model
      common/EMFORCHIEN/ECH4CHIEN(JM0),EN2OCHIEN(JM0)
C     For Emission
#endif
      dimension NDAYMN(12)
      data NDAYMN /31,28,31,30,31,30,31,31,30,31,30,31/
      INTEGER JDOFM(13)
      DATA JDOFM/0,31,59,90,120,151,181,212,243,273,304,334,365/
      DATA FIRST/.TRUE./
      SAVE

C     AJCLR
C            1  SW INC AT P0            RD (AJ(1))
C            2  SW ABS BELOW P0         RD (AJ(2))
C            3  SW ABS BELOW P1         RD (AJ(3))
C            4  SW ABS AT Z0            RD (AJ(6))
C            5  SW INC AT Z0            RD (AJ(5))
C            6  LW INC AT Z0            RD (AJ(67))
C            7  NET LW AT Z0            SF (AJ(9))
C            8  NET LW AT P0            RD (AJ(7))
C            9  NET LW AT P1            RD (AJ(8))
C           10  NET RAD AT P0           DG (AJ(10))
C           11  NET RAD AT P1           DG (AJ(11))
C           12  NET RAD AT Z0           DG (AJ(12))
C **** CLEAR SKY
      INTFX(XTAU)=INT(XTAU*XINT+.5)                                       10.   
      EVENT(XTAU)=MOD(ITAU,INTFX(XTAU)).LT.IDTHR                          11. 

!
! --- assign input and output files
!	Note: Due to historical reasons, no all files are
!		assigned here - in case you want to search
!		something use
!		grep -i "needed characters" *.F
!
!		You have my sympathy.  
!
!	Chien 080400
!

c
      if(FIRST) then
c
      CALL CLOCKS (MNOW)                                                  12.   
C     CALL ERRSET (206,1,0,0,0,301)                                       13.   
C     CALL ERRSET (208,256,-1)                                            14.   
      MBEGIN=MNOW                                                         14.1  
c     CALL HARMO(36,1,24,ACO,BSI,CCO,DSI,1)                               14.5  
      IPFLAG=0                                                            15.   
C     CALL ENQJOB                                                         16.   
      CALL INPUT                                                          17.   
      
#if ( defined CPL_CHEM )
!
! --- Chemistry model
! --- Set year and month index:
!
      myyear  = 1       !year from starting point
      iyearchem  = 1       !year from starting point
      mymonth = 1       !month

      ihaha = 1
      ievenodd = 0      ! even hour 0
                        ! odd  hour 1

      call chembudget (p)
!
! --- Set cfcnsf = 0.0
!
      do k=1,nlev
        cfcnsf(k) = 0.0
      enddo
!
#endif
 
C
      KDISK0=500+KDISK
	ndaa=3
c       LHORDIF=.false.
	LHORDIF=.true.
	if(.not.LHORDIF)print *,' NO HOR. DIFFUSION for Q'
	if(LHORDIF)print *,' HOR. DIFFUSION for Q after COND'
	print *,' RADIATION EVERY ',NRAD/NDYN,' HOURS'
	
       odifcarbon=.false.
#if ( defined CPL_OCEANCO2 && defined ML_2D )
      odifcarbon=.true.
      wind_amp=1.
      dtco2=3600.*24.
c     ncallgary=0
      do j=1,jm
        areaml(j)=dxyp(j)*(1-FDATA(1,J,2))
        focean(j)=(1-FDATA(1,J,2))
        DEPTHML(j)=ZOAV(j)
      end do    !       j
        print *,' RCO2'
!       print 5001,((Rco2(j,k),j=1,jm),k=1,LMO)
        print 5001,((Rco2(j,k)*1.e2,j=1,jm),k=1,LMO)
        dzog(1)=10./SQRT(1.7010587)
        dzg(2)=10.
        do l=2,lmo-1
          dzog(L)=dzog(L-1)*1.7010587
          dzg(L+1)=dzg(L)*1.7010587
        end do
        zg(1)=50.
        dzg(1)=100.
        do l=2,lmo
        zg(l)=zg(l-1)+0.5*(dzg(l-1)+dzg(l))
        end do
        do l=1,lmo
c         edohd(l)=2.5e4/(zg(l)/zg(1))**1.0
C New coefficients for horizontal diffusion 11/16/00
          edohd(l)=1.55e4-9.231e3*(atan((zg(l)-300)/50))
        end do
#endif


#if ( defined CPL_CHEM )
!
! --- Initialization of chemistry model:
!
c     call cheminit(ISTRT1,T,q)
c     11/30/2000 ISTRTCHEM for restart of chemistry model
      call cheminit(ISTRTCHEM,T,q)
c     print *,'H2SO4 after cheminit ',h2so4(1,33,1)


	!
	! --- tmp output
	!
copen(124,file='OUTPUT/hro3',form='unformatted',
c    &		 status='unknown')
     	hrcnt         = 1
     	cnt3hr(:)     = 0
	sfc3hro3(:,:) = 0.0

#endif

       print *,' IRAND=',IRAND
       print *,' NCNDS=',NCNDS
        print *,' after INPUT MRCHT=',MRCHT
      JDAY00=JDAY-1
      AEXP4=AEXP
      TAU4=TAU
      print *,' ISTRT1=',ISTRT1
      if(ISTART.eq.2.or.ISTRT1.eq.0)then
         nwr=1
         WRITE (546) AEXP4,TAU4,XLABEL
         write(547)AEXP,nwr
      elseif(ISTART.eq.10)then
         read(547)AEXPX,nwr
         print *,' NWR=',nwr
         if(abs(AEXPX-AEXP).gt.0.05)then
          print *,' DISAGREEMENT BETWEEN AEXPX AND AEXP FILE 47'
          print *,' AEXPX=',AEXPX,' AEXP=',AEXP
          stop
         endif
         READ (546) AEXP4
         if(abs(AEXP4-AEXP).gt.0.05)then
          print *,' DISAGREEMENT BETWEEN AEXP4 AND AEXP FILE 46'
          print *,' AEXP4=',AEXP4,' AEXP=',AEXP
          stop
         endif
C***
         do 245 nr=1,nwr-1
          READ (546)  
  245    continue
      endif
c     CALL FRTR0(IO)                                                      18.   
      KBGN=KINC+1                                                         18.5  
      KM2=KM*2-1                                                          18.51 
      KM3=KM2                                                             18.52 
      KB1=3                                                               18.53 
      KB2=5                                                               18.54 
      KB3=5                                                               18.55 
      IS=IM                                                               18.56 
      IF(KM.EQ.1) IS=1                                                    18.57 
      FIO=IO                                                              18.58 
      JMM2=JM-2                                                           18.59 
      HR24=24.                                                            18.6  
      HR12=12.                                                            18.61 
      MSTART=MNOW+MDYN+MCNDS+MRAD+MSURF+MDIAG+MELSE                       19.   
      LMT3P1=LM*3+1                                                       20.   
C**** INITIALIZE TIME PARAMETERS                                          21.   
      DTHR=DT/3600.                                                       22.   
      IDTHR=INTFX(DTHR)                                                   23.   
      DTFS=DT*2./3.                                                       24.   
      DTLF=2.*DT                                                          25.   
      NDYNO=MOD(NDYN,2)                                                   26.   
      I24=INTFX(24.)                                                      27.   
      NSTEP0=.5+TAUI/DTHR                                                 28.   
      NSTEP=INT(.5+TAU/DTHR)-NSTEP0                                       29.   
      NSTEP1=NSTEP                                                        29.5  
      NSTEP2=NSTEP                                                        29.6  
      MRCHT=0.                                                            29.7  
      ITAU=(NSTEP+NSTEP0)*IDTHR                                           30.   
cjrs changed to dfloat 8/2/07
      TAU=DFLOAT(ITAU)/XINT                                               31.   
      IDAY=1+ITAU/I24                                                     32.   
      TOFDAY=(ITAU-(IDAY-1)*I24)/XINT                                     33.   
!     if(ISTART.eq.2.or.ISTRT1.eq.0.and..not.CONTRR)then
!     do 458 j=1,JM
!      TSURFD(j)=0.
!      TSURFT(j)=0.
! 458 continue
!     endif
      if(JDATE.eq.100)then
       print *,JDATE,JMONTH,JYEAR
       print *,' main before daily0'
       print *,' T1 ocean'
       print 5001,(ODATA(1,j,1),j=1,JM0)
       print *,' T2 ocean'
       print 5001,(ODATA(1,j,4),j=1,JM0)
       print *,' T3 ocean'
       print 5001,(ODATA(1,j,5),j=1,JM0)
       print *,' sea ice'
       print 5002,(ODATA(1,j,2),j=1,JM0)
      endif
      CALL DAILY_NEW0                                                         34.   
      print *,' Main after DAILYNEW0 JYEAR=',JYEAR
       print *,"DT2MGL"
       print *,DT2MGL
       print *,"DT2MLD"
       print *,DT2MLD
#if( !defined OCEAN_3D&& !defined ML_2D )
      CALL DAILY_OCEAN
      print *,' AFTER DAILY_OCEAN IDAY=',IDAY,' IYEAR=',IYEAR
      print *,' JYEAR=',JYEAR,' JDAY=',JDAY
      print *,' JDATE=',JDATE,' JMONTH=',JMONTH
#endif
      if(JDATE.eq.100)then
       print *,JDATE,JMONTH,JYEAR
       print *,' main after daily0'
       print *,' T1 ocean'
       print 5001,(ODATA(1,j,1),j=1,JM0)
       print *,' T2 ocean'
       print 5001,(ODATA(1,j,4),j=1,JM0)
       print *,' T3 ocean'
       print 5001,(ODATA(1,j,5),j=1,JM0)
       print *,' sea ice'
       print 5002,(ODATA(1,j,2),j=1,JM0)
      endif
c     print *,' Z12O'
c     print *,(Z12O(1,j),j=1,JM0)
c99   CONTINUE                                                            34.993
      CALL CLOCKS (MLAST)                                                 35.   
      MINC=MNOW-MLAST                                                     36.   
      MELSE=MELSE+MINC                                                    37.   
      PERCNT=100.*MELSE/(MSTART-MLAST+1.E-5)                              38.   
c     WRITE (6,901) IDAY,TOFDAY,JDATE,JMONTH,MINC,MELSE,PERCNT,TAU        39.   
      DOPK=1.                                                             40.   
      MODD5K=1000                                                         41.   
      IF(TAU.GE.TAUE) GO TO 820                                           42.   
      HPRNT=.TRUE.
      HPRNT=.FALSE.
      JPR=7
      JPR=1
      JPR=14
      LPR=1
      print *,' MAIN MRCHT=',MRCHT
c
c     print *,' TAUE=',TAUE
      TAUEM=TAUE
      FIRST=.FALSE.
C
c     TAU for coupler
c
      TAUATM=TAU
      MONTHATM=JMONTH
      JDATEATM=JDATE
      JYEARATM=JYEAR
C
#if ( defined  CPL_OCEANCO2 )
       do j=1,jm
        fluxco2(j)=0.0
       enddo
#endif
#if ( defined CPL_CHEM) && ( defined CPL_TEM )
C  For TEM
      if(ISTRT1.eq.0) then
c      New run
c      Reading from flin_nep
        read(537)adupt,temco2
     &     ,temch4,temn2o
      else
c      Restart of the run
c      Reading from last_nep
cjrs   file previously opened in input.F
       read(876)adupt,temco2
     &    ,temch4,temn2o
C       CLOSE(876) 
       rewind 876
      endif

      adupt= 0.0

       aduptd=adupt/(365.*JM)
       temnepgl=0.0
        do j=1,jm
         temnepgl=temnepgl+temco2(j)
        enddo
        print *,'ADNEP=',adupt
        print *,'Initial NEP=',adupt+temnepgl*1.e-3
        temup0=0.0
#endif
#if (!defined OCEAN_3D && !defined ML_2D)
      if(TRANSR)then
       if(LMO.eq.11) then
         call ODIFS
       elseif(LMO.eq.12) then
         call ODIFS12
       else
         Print *,' Wrong LMO',LMO
         stop
       endif
      endif
#endif
!#if (defined PREDICTED_GASES)
#if (defined CPL_TEM || defined CPL_OCEANCO2 )
      if(OBSFOR) then
       call obsco2(iyear,imontha)
       mnobco2=imonth
      endif
#endif
!#endif

#ifdef ML_2D
         do j=1,jm
           do i=1,io
            CLAND4OCEAN(i,j)= C3LAND(I,J)
           enddo
          enddo
         IDAYM=IDAY
         JDAYM=JDAY
         JDATEM=JDATE
         JMONTHM=JMONTH
         JYEARM=JYEAR
         TAUML=TAU
         TOFDAYML=TOFDAY
         KOCEANM=KOCEAN
#endif
#if (defined OCEAN_3D || defined ML_2D)
         do l=1,lm
          sigfl(l)=sig(l)
         enddo
         print *,'SIGFL'
         print *,sigfl
#endif
#if ( defined CPL_CHEM )
        do j=1,jm
          atm_co2(j)=zco2(1,j,1)
     &           *28.97296245/44.0*1.e-9
                        !ppb(m) to kg per volume base
     &           *1.e6
        enddo
#else
      if(.not.OBSFOR) then
        CFF=1.
        if (CO2.gt.0.0)CFF=CO2
        do j=1,jm
          atm_co2(j)=CFF*GHGBGR(1)
        enddo
      endif
#endif
#if (defined CPL_TEM || defined CPL_OCEANCO2 )
        print *,'ATM_CO2'
        print *,atm_co2
#endif
      JDAYLAST=-1
      ncallclm=0
      NOCLM=.true.
#if ( defined CLM )
      NOCLM=.false.
#endif
      print *,' atmosphere DTATM=',DTATM
      print *,' It is running'
      print *,'End of atmospheric model initialization'
      print *,'  '
      print *,'  '
      print *,' JMONTHM= ',JMONTHM
      print *,' TOFDAYML= ',TOFDAYML
      print *,'  '
      return
      endif	! first
C****                                                                     43.   
C**** MAIN LOOP                                                           44.   
C****                                                                     45.   
C
cprint *,' atmosphere TAU=',tau
      TAUE=TAU+DTATM
c     HPRNT=TAU.ge.17520.0.and.TAU.lt.17545.0
c     print *,' TAUE=',TAUE
#if ( defined OCEAN_3D || defined ML_2D)
         do j=1,jm0
           tauu(j)=0.
           tauv(j)=0.
           precip(j)=0.
           evao(j)=0.
           evai(j)=0.
           hfluxo(j)=0.
           hfluxi(j)=0.
           dhfodtg(j)=0.
           devodtg(j)=0.
           dhfidtg(j)=0.
           devidtg(j)=0.
           dhfodtgeq(j)=0.
           devodtgeq(j)=0.
           dhfidtgeq(j)=0.
           devidtgeq(j)=0.
           tempr(j)=0.
cjrs change var name to arunoff
           arunoff(j)=0.
           solarinc_ice(j)=0.
           solarnet_ice(j)=0.
           solarinc_ocean(j)=0.
           solarnet_ocean(j)=0.
           naveo(j)=0.
           navei(j)=0.
           navrad(j)=0.
           navrado(j)=0.
c
           ps4ocean(j)=0.
           do l=1,lm
            qyz4ocean(j,l)=0.
            tyz4ocean(j,l)=0.
           enddo
c
         enddo
#endif
#ifdef OCEAN_3D
C get data from atm-ocean common block
       do j=1,jm0
         ODATA(1,j,1)=mmsst(j)
         ODATA(1,j,2)=mmfice(j)
         GDATA(1,j,3)=mmtice(j)
         GDATA(1,j,1)=mmsnowm(j)
         ODATA(1,j,3)=mmicem(j)
         GDATA(1,j,7)=0.5*(mmtice2(j)+mmtice1(j))
#  ifdef CPL_OCEANCO2 
         fluxco2(j)=fluxco2(j) + dtatm*3600.*mmco2flux(j)
#  endif
       enddo
#endif
      WLMMAX=0.0
C
  100 IF(.NOT.EVENT(TAUT)) GO TO 200                                      46.   
c     HPRNT=TAU.ge.17520.00
      NSTEP1=NSTEP                                                        46.5  
C**** WRITE RESTART INFORMATION ONTO DISK                                 47.   
  120 CALL RFINAL (IRAND)                                                 48.   
      IF(NSTEP.EQ.NSTEP2) GO TO 116                                       48.3  
      DO 115 K=1,22                                                       48.5  
      DO 115 J=1,JM                                                       48.6  
      DO 115 I=2,IO                                                       48.7  
 115  GDATA(I,J,K)=GDATA(1,J,K)                                           48.9  
 116  CONTINUE                                                            48.91 
      if(wr25.and.ISTART.eq.2)then
       print *,' main write'
       print *,' T1 ocean'
       print 5001,(ODATA(1,j,1),j=1,JM0)
       print *,' T2 ocean'
       print 5001,(ODATA(1,j,4),j=1,JM0)
       print *,' T3 ocean'
       print 5001,(ODATA(1,j,5),j=1,JM0)
      REWIND KDISK0                                                       49.   
      if(TRANSR)then
      WRITE(KDISK0) AEXP,TAU,JC,C,RC,KEYNR,U,V,T,P,Q,ODATA,GDATA,BLDATA,  50.
     *  RQT,SRHR,TRHR,(AJ(K,1),K=1,KACC),TAU,TSSFC,CKS,CKN,WMGE,TPRIM2    51.
     *  ,MRCHT,TRSURF,SRSURF,TLANDW,TSURFW,DWAV0
     *  ,TG3M,RTGO,STG3,DTG3
      print *,' STG3'
      print 5001,(STG3(1,j),j=1,JM0)
      print *,' DTG3/356'
      print 5001,(DTG3(1,j)/365.,j=1,JM0)
      print *,' RTGO'
      print 5001,((RTGO(1,j,k),j=1,JM0),k=1,lmo)
 5001 format(23f6.1)
 5002 format(23f6.3)
      else
      WRITE(KDISK0) AEXP,TAU,JC,C,RC,KEYNR,U,V,T,P,Q,ODATA,GDATA,BLDATA,  50.   
     *  RQT,SRHR,TRHR,(AJ(K,1),K=1,KACC),TAU,TSSFC,CKS,CKN,WMGE,TPRIM2    51.   
     *  ,MRCHT,TRSURF,SRSURF,TLANDW,TSURFW,DWAV0
      endif
      REWIND KDISK0                                                       52.   
      end if	! ISTART.eq.2
      CALL CLOCKS (MNOW)                                                  53.   
      MINC=MLAST-MNOW                                                     54.   
      MELSE=MELSE+MINC                                                    55.   
      PERCNT=100.*MELSE/(MSTART-MNOW+1.E-5)                               56.   
      MLAST=MNOW                                                          59.   
C**** TEST FOR TERMINATION OF RUN                                         60.   
  200 continue
c     HPRNT=TAU.gt.45.0.and.TAU.lt.60.0
c     HPRNT=TAU.gt.470.0.and.TAU.lt.550.0
      NCOMP=0
      IF(TAU+.06125.GE.TAUE) GO TO 820                                    63.   
      JDAY00=JDAY
C**** IF TIME TO ZERO OUT DIAGNOSTIC ACCUMULATING ARRAYS, DO SO           64.   
C****   (NORMALLY DONE AT THE BEGINNING OF A MONTH)                       65.   
      IF (TAU.EQ.TAUI) GO TO 260                                          66.   
      IF(.NOT.EVENT(24.)) GO TO 290                                       67.   
      DO 250 K=1,13                                                       68.   
      IF(JDAY.EQ.NDZERO(K)) GO TO 260                                     69.   
  250 CONTINUE                                                            70.   
      GO TO 290                                                           71.   
  260 CONTINUE
      TAU0=TAU                                                            72.   
      IDAY0=IDAY                                                          73.   
      TOFDY0=TOFDAY                                                       74.   
      JDATE0=JDATE                                                        75.   
      JMNTH0=JMONTH                                                       76.   
      JYEAR0=JYEAR                                                        77.   
      DO 270 I=1,12                                                       78.   
  270 IDACC(I)=0                                                          79.   
      NODIFS=0
      nsrps=0.
      DO 280 K=1,KACC                                                     80.   
  280 AJ(K,1)=0.                                                          81.   
      do 5280 j=1,JM
         tprmg(j)=0.
         ntprmg(j)=0.
         Z1OAV(j)=0.
         NZ1OAV(j)=0
         NFLIO(J)=0
         FLIO(J)=0.
         NCLR(J)=0
         do 5282 k=1,lmo
            RTGOAV(j,k)=0.
 5282    continue
         do 5281 n=1,12
              AJCLR(J,n)=0.
              BJCLR(J,n)=0.
              CJCLR(J,n)=0.
 5281    continue 
 5280 continue
      do 5286 j=1,jm+3
         srps(j)=0.
 5286 continue

c#if ( defined CPL_OCEANCO2 && defined ML_2D )
c       call zerogary(ncallgary)
c       call zerogary
c#endif

         CALL DIAG9A (1)                                                  82.   

#if ( defined CPL_NEM )
C For Emission
c === 031097
      DO 5655 MONTH=1,12
      IF(JDAY.LE.JDOFM(MONTH+1)) GO TO 5656
 5655 CONTINUE
 5656 MNHTEM=MONTH-1
      if(MNHTEM.eq.0)MNHTEM=12
      do J=1,JM
c      ECH4CHIEN(J)=temch4(J)/NDAYMN(MNHTEM)/1000.
c      EN2OCHIEN(J)=temn2o(J)/NDAYMN(MNHTEM)/1000.

       ECH4CHIEN(J)=temch4(J)/NDAYMN(MNHTEM)
       EN2OCHIEN(J)=temn2o(J)/NDAYMN(MNHTEM)

      enddo

C For Emission
#endif

  290 CONTINUE                                                            83.   
C****                                                                     84.   
C**** INTEGRATE DYNAMIC TERMS                                             85.   
      if(HPRNT)then
      print *,' main before comp1 1',' TAU=',TAU,JDATE,JMONTH
#include "PRNT.COM"
      print *,' GDATA(1,7,5)=',GDATA(1,7,5),' GDATA(1,7,6)='
     *  ,GDATA(1,7,6)
      endif
C****                                                                     86.   
         MODD5D=MOD(NSTEP,NDA5D)                                          87.   
         IF(MODD5D.EQ.0) CALL DIAG5A (2,0)                                88.   
c     IF(NDYNO.EQ.1) then
c     print *,NDYNO
c     print *,jm,im,LMT3P1
c     endif
      DO 310 J=1,JM                                                       89.   
      DO 300 L=1,LMT3P1                                                   90.   
      DO 300 I=1,IM                                                       91.   
      UX(I,J,L)=U(I,J,L)                                                  92.   
  300 UT(I,J,L)=U(I,J,L)                                                  93.   
      DO 310 L=1,LM                                                       94.   
      DO 310 I=1,IM                                                       95.   
  310 QT(I,J,L)=Q(I,J,L)                                                  96.   
      if(HPRNT)then
      print *,' main before comp1 2',' TAU=',TAU
#include "PRNT.COM"
      print *,' GDATA(1,7,5)=',GDATA(1,7,5),' GDATA(1,7,6)='
     *  ,GDATA(1,7,6)
      endif
C**** INITIAL FORWARD STEP, QX = Q + .667*DT*F(Q)                         97.   
      NS=0                                                                98.   
      MRCH=0                                                              99.   
      CALL COMP1 (UX,VX,TX,PX,Q,U,V,T,P,Q,DTFS,NS)                       100.   
      if(HPRNT)then
      print *,' main after comp1',' TAU=',TAU,' MRCH=',MRCH      
#include "PRNT.COM"
      endif
C
      IF(NDYNO.EQ.1) GO TO 320                                           101.   
C**** INITIAL BACKWARD STEP IS ODD, QT = QT + DT*F(QX)                   102.   
      MRCH=-1                                                            103.   
      CALL COMP1 (UT,VT,TT,PT,QT,UX,VX,TX,PX,Q,DT,NS)                    104. 
      if(HPRNT)then
      print *,' main after comp1',' TAU=',TAU,' MRCH=',MRCH
#include "PRNT.COM"
      endif
C
      GO TO 360                                                          105.   
C**** INITIAL BACKWARD STEP IS EVEN, Q = Q + DT*F(QX)                    106.   
  320 NS=1                                                               107.   
         MODD5K=MOD(NSTEP+NS-NDYN+NDA5K,NDA5K)                           108.   
      MRCH=1                                                             109.   
      CALL COMP1 (U,V,T,P,Q,UX,VX,TX,PX,QT,DT,NS)                        110.  
      if(HPRNT)then
      print *,' main after comp1',' TAU=',TAU,' MRCH=',MRCH
#include "PRNT.COM"
      endif
C
CD       DIAGA SHOULD BE CALLED HERE BUT THEN ARRAYS MUST BE CHANGED     111.   
c
c
C**** ODD LEAP FROG STEP, QT = QT + 2*DT*F(Q)                            112.   
  340 MRCH=-2                                                            113.   
      CALL COMP1 (UT,VT,TT,PT,QT,U,V,T,P,Q,DTLF,NS)                      114.   
      if(HPRNT)then
      print *,' main after comp1',' TAU=',TAU,' MRCH=',MRCH
#include "PRNT.COM"
      endif
C
C**** EVEN LEAP FROG STEP, Q = Q + 2*DT*F(QT)                            115.   
  360 NS=NS+2                                                            116.   
         MODD5K=MOD(NSTEP+NS-NDYN+NDA5K,NDA5K)                           117.   
      MRCH=2                                                             118.   
      CALL COMP1 (U,V,T,P,Q,UT,VT,TT,PT,QT,DTLF,NS)                      119.   
      if(HPRNT)then
      print *,' main after comp1',' TAU=',TAU,' MRCH=',MRCH
#include "PRNT.COM"
      endif
C
      IF(NS.LT.NDYN) GO TO 340                                           122.   
c     IF(MOD(NSTEP+NS-NDYN+NDAA,NDAA).LT.MRCH) THEN                      120.   
        CALL DIAGA (UT,VT,TT,PT,QT,NOCLM)                                121.   
        if(HPRNT)then
         print *,' main after DIAGA',' TAU=',TAU,' MRCH=',MRCH
#include "PRNT.COM"
        endif
c     ENDIF
c     IF(NS.LT.NDYN) GO TO 340                                           122.   
      DOPK=1.                                                            123.   
      CALL CLOCKS (MNOW)                                                 124.   
      MINC=MLAST-MNOW                                                    125.   
      MDYN=MDYN+MINC                                                     126.   
      MLAST=MNOW                                                         127.   
      PERCNT=100.*MDYN/(MSTART-MNOW+1.E-5)                               128.   
C                                                                        130.   
         CALL DIAG9A (2)                                                 131.   
C****                                                                    133.   
C**** INTEGRATE SOURCE TERMS                                             134.   
C****                                                                    135.   
#if (!defined PREDICTED_GASES)
#if (defined CPL_TEM || defined CPL_OCEANCO2 )
      if(OBSFOR) then
      if(JMONTH.ne.AMONTH(mnobco2)) then
         mnobco2=mnobco2+1
         if(mnobco2.eq.13)mnobco2=1
         call obsco2(jyear,mnobco2)
       endif
       endif
#endif
#endif
C
C
      MODRD=MOD(NSTEP,NRAD)                                              136.   
         MODD5S=MOD(NSTEP,NDA5S)                                         137.   
         IF(MODD5S.EQ.0) IDACC(8)=IDACC(8)+1                             138.   
C                                                                        139.   
C**** CONDENSTATION, SUPER SATURATION AND MOIST CONVECTION               140.   
      if(HPRNT)then
      print *,' main before conds',' TAU=',TAU
#include "PRNT.COM"
      endif


#if ( defined CPL_CHEM )
!
! ===== Calculate airmass at grid point
! =====   Chien Wang, 092395
!

	call chemairmass(p4chem0)

!	
! === Calculating total mass of tracers
! === 	with long residence times:
!
	call chemmass1(cfc11, cfc11mass)
	call chemmass1(cfc12, cfc12mass)
	call chemmass1(xn2o,  xn2omass)
	call chemmass1(zco2,  zco2mass)
	call chemmass1(ch4,   ch4mass)

        ! === if hfc, pfc, and sf6 are included:
#ifdef INC_3GASES
	! === 032698
	call chemmass1(hfc134a, hfc134amass)
	call chemmass1(pfc,     pfcmass)
	call chemmass1(sf6,     sf6mass)
	! ===
#endif

!
! === Calculating advection and eddy diffusion:
!
	call chemadv0 (dt)

!
! === Calculate total n-s transport amount
! === 	of cfc11 - temperary:
!
	dth = 3600.0	! for relexible setting of dt dt * 3.0
!       call chemtmp1 (dth,airmass0,p,pvv,cfc11)

! === Readjust mass of tracers, 1:
!
! 	-----------------------------
!  Use tropospheric life time (yr) to calculate mass
!	loss of tracers and at the same time compensate
!	all the numerical loss back to tracer's mass
!	which equavalent to use adjcoe = 1.0 for 
!	chemmass2.f
!
!	Chien Wang,  September 12,1995
!	------------------------------

! === 092595 update p

	call chemairmass(p4chem1)

	call chemmass6(46.0, 1.0, cfc11,cfc11mass)
	call chemmass6(120.0,1.0, cfc12,cfc12mass)
 
! === 102596
! === close tau type of ocean uptake of co2:

	call chemmass66(1.0, 1.0,zco2,zco2mass)

!call chemmass6(150.0,1.0,xn2o,xn2omass)
        call chemmass6(120.0,1.0,xn2o,xn2omass)
	call chemmass2(1.0,ch4, ch4mass )

        ! === if hfc, pfc, and sf6 are included:
#ifdef INC_3GASES
	! === 032698
!	call chemmass6(14.6,   1.0,hfc134a, hfc134amass)
	call chemmass2(1.0, hfc134a,   hfc134amass)
	call chemmass6(10000.0,1.0,pfc,pfcmass)
	call chemmass6(3200.0, 1.0,sf6,sf6mass)
	! ===
#endif

!
! === Calculate tropospheric gaseous reactions
! === 	every nhr_for_chem hours:
!
       dt_chem_h = dth*float(nhr_for_chem)
       if(ievenodd.eq.0) then
       call chemtrop0(0, T, q, dt_chem_h, 1)

c     print *,'H2SO4 after chemtrop0 ',h2so4(1,33,1)
c     print *,'SVIOD after chemtrop0 ',SVIOD(1,33,1)

	!
	! --- tmp output
	!
	cnt3hr(hrcnt)	       = cnt3hr(hrcnt) + 1
	sfc3hro3(1:nlat,hrcnt) = sfc3hro3(1:nlat,hrcnt) 
     &		        + (o3(1,1:nlat,1))*29.0/48.0
     	
	hrcnt = hrcnt + 1
	if (hrcnt .gt. 8 ) hrcnt = 1
	end if

       ievenodd = ievenodd + 1
       if(ievenodd.ge.nhr_for_chem) ievenodd = 0

!
! === Calculating stratospheric processes:
! === 092595
! === adjust startospheric loss to whole global:
!
      call chemmass1(cfc11, cfc11mass)
      call chemmass1(cfc12, cfc12mass)
      call chemmass1(xn2o,  xn2omass)

      call chemstrat (dt)

      call chemmass2(1.0,cfc11,cfc11mass)
      call chemmass2(1.0,cfc12,cfc12mass)
      call chemmass2(1.0,xn2o, xn2omass)

!
! === Get total mass of chemically active species:
!
!	call chemairmass(p)
!	write(6,*)"P before 2nd", p
c     print *,'H2SO4 after chemairmass ',h2so4(1,33,1)

      call chemmass1(cfc11, cfc11mass)
      call chemmass1(cfc12, cfc12mass)
      call chemmass1(xn2o,  xn2omass)
      call chemmass1(zco2,  zco2mass)
      call chemmass1(ch4,   ch4mass)

        ! === if hfc, pfc, and sf6 are included:
#ifdef INC_3GASES
	! === 032698
        call chemmass1(hfc134a, hfc134amass)
        call chemmass1(pfc,     pfcmass)
        call chemmass1(sf6,     sf6mass)
	! ===
#endif

#endif

c     print *,'H2SO4 before CONDSE ',h2so4(1,33,1)
      CALL CONDSE(mndriver)                                              141.   
c     print *,'H2SO4 after CONDSE ',h2so4(1,33,1)

#if ( defined CPL_CHEM )
!
! === Calculate emission once per hour:
!
!      timeinhr=1./(365.*24.)      !hourly  emission
!      call chememission(timeinhr)
!
!
! === Print hourly
!
!      call chemprt
!
! ==============================================
#endif

      if(HPRNT)then
      print *,' main after conds',' TAU=',TAU
#include "PRNT.COM"
      endif
#if ( !defined CLM )
      CALL PRECIP_LAND(mndriver)                                         142.   
      if(HPRNT)then
      print *,' main after preci',' TAU=',TAU
#include "PRNT.COM"
      endif
#endif
      CALL CLOCKS (MNOW)                                                 143.   
      MINC=MLAST-MNOW                                                    144.   
      MCNDS=MCNDS+MLAST-MNOW                                             145.   
      MLAST=MNOW                                                         146.   
C                                                                        147.   
         CALL DIAG9A (3)                                                 148.   
C**** RADIATION, SOLAR AND THERMAL                                       149.   
      if(HPRNT)then
      print *,' main before radia',' TAU=',TAU
#include "PRNT.COM"
      print *,' IRAND=',IRAND
      endif

#if ( defined PREDICTED_GASES || defined PREDICTED_AEROSOL )

#    if ( defined FIXED_FOR )
         CALL RADIA
#    else
      call radia_chem
#    endif

#else

#    if ( defined SVI_ALBEDO )
         CALL RADIAGSO
#    else
       if(OBSFOR) then
#        if ( defined FIXED_FOR )
          CALL RADIA
#        else
          CALL RADIAGSO
#        endif
       else
         CALL RADIA                   
       endif
#    endif

#endif

      if(HPRNT)then
      print *,' main after radia',' TAU=',TAU
#include "PRNT.COM"
      endif
      CALL CLOCKS (MNOW)                                                 151.   
      MINC=MINC+MLAST-MNOW                                               152.   
      MRAD=MRAD+MLAST-MNOW                                               153.   
      MLAST=MNOW                                                         154.   
      if(HPRNT)then
      print *,' main before diag9a',' TAU=',TAU
#include "PRNT.COM"
      endif
C                                                                        155.   
         CALL DIAG9A (4)                                                 156.   
      if(HPRNT)then
      print *,' main after diag9a',' TAU=',TAU
#include "PRNT.COM"
      endif
C**** SURFACE INTERACTION AND GROUND CALCULATION                         157. 

#if ( defined CLM )
      if(HPRNT)then
      print *,' main before surf4clm',' TAU=',TAU
#include "PRNT.COM"
      endif
      CALL SUR4CLM
      if(HPRNT)then
      print *,' main after surf4clm',' TAU=',TAU
#include "PRNT.COM"
      endif
      i=1
      do j=1,jm
         pcpl4clm(i,j)=pcpl4clm(i,j)*prlnd2total(j,mndriver)
     &     *3600./(NDYN*DT)
         pcpc4clm(i,j)=pcpc4clm(i,j)*prlnd2total(j,mndriver)
     &     *3600./(NDYN*DT)
      enddo
!     print *,' main after surf4clm',' TAU=',TAU
!       print  ('2(12f7.2,/,11f7.2,/)'),ps4clm,pcpl4clm,
!    & pcpc4clm,tpr4clm,
!    &  tsl4clm,
!    &  qs4clm,ws4clm
!    & ,us4clm,vs4clm,
!    &  dsw4clm,
!    &  dlw4clm,pco24clm
!    &  ,swinr4clm,swvis4clm
#if ( defined DATA4TEM )
c     if(JYEAR.gt.20)then
        write (935),ps4clm,pcpl4clm,
     & pcpc4clm,tpr4clm,
     &  tsl4clm,
     &  qs4clm,ws4clm
     & ,us4clm,vs4clm,
     &  dsw4clm,
     &  dlw4clm,pco24clm
     &  ,swinr4clm,swvis4clm
c      endif
#endif

       ncallclm=ncallclm+1
c       print *,'before clm4mit2d ncallclm=',ncallclm
      if(HPRNT)then
      print *,' main before clm4mit2d',' TAU=',TAU
#include "PRNT.COM"
      endif
       call clm4mit2d
      if(HPRNT)then
      print *,' main after clm4mit2d',' TAU=',TAU
#include "PRNT.COM"
      endif

c     if(JYEAR.gt.20)then
!      write (934),tau,tsoiclm,snwdclm,snwcclm,
!    & lwuclm,tref2mclm,tflxclm,tgndclm,
!    & lhfclm,shfclm,tauxclm,tauyclm,
!    & asdirclm,aldirclm,asdifclm,aldifclm,
!    & sroclm,ssrclm,glrclm
!    &,h2olclm,h2oiclm
c     endif
!     print *,' main after clm4mit2d',' TAU=',TAU
!       print  ('2(12f7.2,/,11f7.2,/)'),tsoiclm,snwdclm,snwcclm,
!    & lwuclm,tref2mclm,tflxclm,tgndclm,
!    & lhfclm,shfclm,tauxclm,tauyclm,
!    & asdirclm,aldirclm,asdifclm,aldifclm

      CALL SURF_CLM
      CALL SURF_OCEAN
      CALL GR_CLM
      if(HPRNT)then
         print *,' main after surfc',' TAU=',TAU
#include "PRNT.COM"
      endif
#else
!     CALL SURF_LAND
!     CALL SURF_OCEAN
      CALL SURFCE       ! 07/14/2006
      if(HPRNT)then
         print *,' main after surfc',' TAU=',TAU
#include "PRNT.COM"
      endif
      CALL GRLAND                                                        159.   
#endif

#if ( defined OCEAN_3D || defined ML_2D)
      CALL GRFOROCEAN
#else
      CALL GROCEAN(mndriver)                                             159.   
#endif
c
      if(HPRNT)then
      print *,' main after groun',' TAU=',TAU
#include "PRNT.COM"
      endif
c     print *,'H2SO4 before DRYCNV ',h2so4(1,33,1)
      CALL DRYCNV                                                        160.   
c     print *,'H2SO4 after DRYCNV ',h2so4(1,33,1)
      if(HPRNT)then
      print *,' main after drycn',' TAU=',TAU
#include "PRNT.COM"
      endif
      CALL CLOCKS (MNOW)                                                 161.   
      MINC=MINC+MLAST-MNOW                                               162.   
      MSURF=MSURF+MLAST-MNOW                                             163.   
      MLAST=MNOW                                                         164.   
         CALL DIAG9A (5)                                                 165.   
C**** STRATOSPHERIC MOMENTUM DRAG                                        166. 
      CALL SDRAG(WLMMAX,JWLMMAX)                                         167.   
      if(HPRNT)then
      print *,' main after sdrag',' TAU=',TAU
#include "PRNT.COM"
      endif
      CALL CLOCKS (MNOW)                                                 168.   
      MINC=MINC+MLAST-MNOW                                               169.   
      MSURF=MSURF+MLAST-MNOW                                             170.   
      MLAST=MNOW                                                         171.   
C                                                                        172.   
         CALL DIAG9A (6)                                                 173.   
      MSRCE=MCNDS+MRAD+MSURF                                             174.   
      PERCNT=100.*MSRCE/(MSTART-MNOW+1.E-5)                              175.   
C**** SEA LEVEL PRESSURE FILTER                                          177.   
      NFILTR=  NDYN                                                      177.5  
      IF(MFILTR.LE.0.OR.MOD(NSTEP,NFILTR).NE.0) GO TO 500                178.   
         IDACC(10)=IDACC(10)+1                                           179.   
      if(HPRNT)then
      print *,' main before filtr',' TAU=',TAU
#include "PRNT.COM"
      endif
C
C                                                                        180.   
C ****************
      if(LHORDIF)then
        DTDIF=3600.
       if(JM.eq.24)then
        CALL HORDIFF(DTDIF)
       else if(JM.eq.46)then
        CALL HORDIFFALL(DTDIF)
       else
        print *,' Wromg JM=',JM
        stop
       endif
      end if
C ****************
c     CALL FILTER                                                        181.   
C
#if ( defined CPL_CHEM )
!
! === Readjust total mass of tracers:
!
      call chemairmass(p)

      call chemmass2(1.00,cfc11,cfc11mass)
      call chemmass2(1.00,cfc12,cfc12mass)
      call chemmass2(1.00,xn2o ,xn2omass )
      call chemmass2(1.00,zco2,zco2mass)
      call chemmass2(1.00,ch4, ch4mass )

        ! === if hfc, pfc, and sf6 are included:
#ifdef INC_3GASES
	! === 032698
        call chemmass2(1.0,hfc134a, hfc134amass)
        call chemmass2(1.0,pfc,     pfcmass)
        call chemmass2(1.0,sf6,     sf6mass)
	! ===
#endif

!
! === Accumulative calculation prepared for
! ===  carrying out monthly average:
!
      call chemmonth1
!
#endif

      CALL CLOCKS (MNOW)                                                 182.   
      MDYN=MDYN+MLAST-MNOW                                               183.   
      MLAST=MNOW                                                         184.   
C                                                                        185.   
         CALL DIAG9A (7)                                                 186.   
C****                                                                    187.   
C**** UPDATE MODEL TIME AND CALL DAILY IF REQUIRED                       188.   
C****                                                                    189.   
  500 NSTEP=NSTEP+NDYN                                                   190.   
      ITAU=(NSTEP+NSTEP0)*IDTHR                                          191.   
cJRS fix to DFLOAT 8/2/07
      TAU=DFLOAT(ITAU)/XINT                                               192.   
      IDAY=1+ITAU/I24                                                    193.   
      TOFDAYPR=TOFDAY+1.00
      TOFDAY=(ITAU-(IDAY-1)*I24)/XINT                                    194.   
      IF(.NOT.EVENT(24.)) GO TO 550                                      195.   
C                                                                        196.   
      do J=1,JM0
        TSURFW(J)=TSURFD(J)
        TLANDW(J)=TLANDD(J)
      enddo

      JDATECLM=JDATE
      JDATEPR=JDATE
      JMONTHPR=JMONTH
      JYEARPR=JYEAR
      CALL DAILY_NEW                                                     197.   
!     print *,' AFTER DAILY_NEW IDAY=',IDAY,' IYEAR=',IYEAR
!     print *,' JYEAR=',JYEAR,' JDAY=',JDAY
#if( !defined OCEAN_3D && !defined ML_2D )
      CALL DAILY_OCEAN
c     print *,' AFTER DAILY_OCEAN IDAY=',IDAY,' IYEAR=',IYEAR
c     print *,' JYEAR=',JYEAR,' JDAY=',JDAY
c     print *,' JDATE=',JDATE,' JMONTH=',JMONTH
#endif
      if(JDATE.eq.100)then
       print *,JDATE,JMONTH,JYEAR
       print *,' main after daily'
       print *,' T1 ocean'
       print 5001,(ODATA(1,j,1),j=1,JM0)
       print *,' T2 ocean'
       print 5001,(ODATA(1,j,4),j=1,JM0)
       print *,' T3 ocean'
       print 5001,(ODATA(1,j,5),j=1,JM0)
       print *,' sea ice'
       print 5002,(ODATA(1,j,2),j=1,JM0)
      endif
      CALL CLOCKS (MNOW)                                                 198.   
      MELSE=MELSE+(MLAST-MNOW)                                           199.   
      MLAST=MNOW                                                         200.   
         NDAILY=SDAY/DT                                                  201.   

#if ( defined CPL_CHEM )
!
! === Calculate air mass, second step:
!
      i=1
      call chemairmass(p)

!
! === Calculate emission once per day:
!
c      zxy=0.
c      do j=1,jm
c       zxy=zxy+zco2(1,j,1)
c    &           *28.97296245/44.0*1.e-3
c      enddo
c      print *,' CO2 before emission ',zxy/jm

      timeinday=1./(365.)          !daily emission
#if ( defined ACCRI )
       do j=1,jm
      ps4accri(1,j) = p(1,j)
       enddo
#endif
      call chememission(timeinday)

c      zxy=0.
c      do j=1,jm
c       zxy=zxy+zco2(1,j,1)
c    &           *28.97296245/44.0*1.e-3
c      enddo
c      print *,' CO2 after emission ',zxy/jm
!
#endif

C                                                                        202.   
         CALL DIAG9A (8)                                                 203.   
#if ( !defined OCEAN_3D && !defined ML_2D )
      IF(KOCEAN.EQ.1) THEN
         DO 540 J=1,JM                                                   203.11 
         DO 540 I=1,IM                                                   203.12 
         AIJ(I,J,59)=AIJ(I,J,59)+ODATA(I,J,4)                            203.13 
  540    AIJ(I,J,60)=AIJ(I,J,60)+ODATA(I,J,5)                            203.14 
      if(TRANSR)then
       if(LMO.eq.11) then
         CALL ODIFS
       elseif(LMO.eq.12) then
         CALL ODIFS12
       else
         Print *,' Wromng LMO',LMO
         stop
       endif
       NODIFS=NODIFS+1
        do 5283 j=1,JM0
          do 5283 k=1,lmo
            RTGOAV(j,k)=RTGOAV(j,k)+RTGO(1,j,k)
 5283  continue
      endif
C**** RESTRUCTURE THE OCEAN LAYERS AND ELIMINATE SMALL ICE BERGS         203.15 
      CALL OSTRUC                                                        203.16 
      if(JDATE.eq.41)then
       print *,JDATE,JMONTH,JYEAR
       print *,' main after ostruc'
       print *,'Z1O'
       print 5001,(Z1O(1,j),j=1,JM0)
       print *,' T1 ocean'
       print 5001,(ODATA(1,j,1),j=1,JM0)
       print *,' T2 top ocean'
       print 5001,((2.*ODATA(1,j,4)-ODATA(1,j,5)),j=1,JM0)
       print *,' T2 ocean'
       print 5001,(ODATA(1,j,4),j=1,JM0)
       print *,' T3 ocean'
       print 5001,(ODATA(1,j,5),j=1,JM0)
      endif
      ENDIF	! KOCEAN
#endif
#if ( defined ML_2D )
      if(TRANSR)then
       NODIFS=NODIFS+1
        do 5283 j=1,JM0
          do 5283 k=1,lmo
            RTGOAV(j,k)=RTGOAV(j,k)+RTGO(1,j,k)
 5283  continue
      endif
#endif

#if ( defined CPL_OCEANCO2 )

C    For OCM or 3D ocean model with carbon cycle
#ifdef PREDICTED_GASES
c -------
c 102596
c  air co2 mixing ratio goes to ocean carbon model:
c
      do j=1,jm
        pC_atm(j)=zco2(1,j,1)
     &           *28.97296245/44.0*1.e-9	
			!ppb(m) to kg per volume base

          atm_co2(j)=pC_atm(j)*1.e6

      enddo    !       j
c
c -------
#else

      do j=1,jm
         pC_atm(j)=atm_co2(j)*1.e-6
      enddo    !       j

#endif

#if ( defined ML_2D )
CB Gary CO2 uptake by the ocean

      do j=1,jm
        if(NWMGEA(J).gt.0)then
          WSAV(J)=WSAV(J)/NWMGEA(J)
          NWMGEA(J)=0.
        else
          WSAV(J)=0.
        end if
      end do    ! j

      do j=1,jm
        tggary(j)=ODATA(1,j,1)+273.16
        wsgary(j)=WSAV(J)
        WSAV(j)=0.
        arsrf(j)=areaml(j)*(1.-ODATA(1,j,2))
        DEPTHML(j)=ZOAV(J)
        co24ocnan(j)=co24ocnan(j)+pC_atm(j)*1.e6
      enddo    !       j
c     print *,'CO2 for 2D ocean'
c     print *,JYEAR,JMONTH
c     print *,'co2=',pC_atm(27)*1.e6,' ws=',wsgary(27)
c     print *,'tem=',tggary(27)
c     print '12f7.1,/,2(11f7.1,/,),12f7.1',(pC_atm(j)*1.e6,j=1,jm)
c     print '12f7.1,/,2(11f7.1,/,),12f7.1',(rco2(j,1),j=1,jm)
c     ncallgary=ncallgary+1
!     10/28/06
!     call carb_mxdlyr_chem(focean)
!     call carb_airsea_flx
!
!     3D ocean chemistry
      call carb_chem_ocmip(focean)
      call carb_airsea_flx(dtco2)
!     3D ocean chemistry
!     10/28/06
c     print *,'FCO2 ncallgary=',ncallgary
c     print '12f7.1,/,2(11f7.1,/,),12f7.1',
c    &       (fluxco2(j)*12.e-15*365.,j=1,jm)
#endif

C    For ocean carbon model
c  Annual oceanic CO2 uptake
       do j=1,jm
         OCUPT=OCUPT+fluxco2(j)
          ocumn=ocumn+fluxco2(j)
          fluxco2mn(j)=fluxco2mn(j)+fluxco2(j)
       enddo
       print *,' OCUPT=',OCUPT*12.e-15

#if ( defined CPL_CHEM )
!
! === Calculate ocean uptake of CO2
! === 	once per day:
!
      i=1
      call chemairmass(p)	!update airmass

      call chemoceanco2(fluxco2)
!
#endif

       do j=1,jm
        fluxco2(j)=0.0
       enddo
#endif

#if ( defined CPL_TEM )
!#if ( defined CLM )
c     print *,'JDATE for TEM=',JDATECLM
      do j=1,jm
        if(npred4tem(j).gt.0)then
c        pred4tem(j)=pred4tem(j)/npred4tem(j)
         ewvd4tem(j)=ewvd4tem(j)/npred4tem(j)
         pre4tem(J)=pre4tem(J)+pred4tem(j)
        endif
c
        if(nt2md4tem(j).gt.0)then
         t2md4tem(j)=t2md4tem(j)/nt2md4tem(j)
        endif
        temp4tem(j)=temp4tem(j)+t2md4tem(j)
        dtem4tem(JDATECLM,j)=t2md4tem(j)
c
        if(nradd4tem(j).gt.0)then
         cldd4tem(j)=cldd4tem(j)/ncldd4tem(j)
         swtd4tem(j)=swtd4tem(j)/nradd4tem(j)
         swsd4tem(j)=swsd4tem(j)/nradd4tem(j)
        endif
        sws4tem(j)=sws4tem(j)+swsd4tem(j)
c
      enddo
c      
       do j=1,jm
         pred4tem(j)=0.0
         ewvd4tem(j)=0.0
         t2md4tem(j)=0.0
         cldd4tem(j)=0.0
         swtd4tem(j)=0.0
         swsd4tem(j)=0.0
         npred4tem(j)=0
         ncldd4tem(j)=0
         nradd4tem(j)=0
         nt2md4tem(j)=0.
      enddo
#endif

#if ( defined CPL_OCEANCO2 && defined ML_2D )
C    For OCM 

!       dtco2=3600.*24.
!       cfkvct=1.0
!       if (JYEAR.ge.1991) then begin
!        if (JYEAR.le.2100) then begin
!         cfkvct=(1.0*(2100-JYEAR)+0.25*(JYEAR-1990))/110.
!        esle
!         cfkvct=0.25
!        endif
!       endif
!       do j=1,jm
!        edzcatr(j)=cfkvct*edzcar(j)
!       enddo
        call diffusco2(lmo,jm,dtco2,0.5,edzcart,depthml,focean,
     &                dzg,dzog,rco2)
        call hdocean(rco2,focean,dxv,dyv,DXYP,depthml,edohd,dtco2)
        call avegary
CB Gary CO2 uptake by the ocean
#endif

#if ( defined CPL_CHEM) && ( defined CPL_TEM )
C take into accout land uptake form TEM for previous  month
      do j=1,jm
         fluxnep(j)=aduptd+1.e-3*temco2(j)/NDAYMN(mndriver)
c  Annual TEM CO2 uptake
         TEMUPTANN=TEMUPTANN+fluxnep(j)
       enddo
        if(jdate.eq.1)then
          print *,'Monthly TEM uptake'
c         print *,mndriver,adupt,temuptann-temup0
          temup0=temuptann
        endif
C
c
           i=1
      call chemairmass(p)       !update airmass

      call chemtemco2(fluxnep)
C
c
!
#endif

c     End of month
      if(JDATE.eq.1)then

#if ( defined CPL_CHEM )
!
! === Calculating monthly averaged mixing ratios:
!
      call chemmonth2
c     print *,'Atmosphere after chemmonth2'

! === Calculate and print monthly n-s transport
! ===   of cfc11:
!
!      call chemtmp2
!
        do nhr = 1,8
          sfc3hro3(1:nlat,nhr) = sfc3hro3(1:nlat,nhr)
     &                             /float(cnt3hr(nhr))
#if ( defined CPL_TEM )
        o34tem(nhr,1:nlat)=sfc3hro3(1:nlat,nhr)
#endif
        end do
cwrite(124)sfc3hro3
        cnt3hr(1:8)          = 0
        sfc3hro3(1:nlat,1:8) = 0.0

! === Writing rawdata every month:
!
!      call chemprt     !closed 032697

      call chembudget (p)
      print *,' Atmosphre after chembudget mymonth=',mymonth
! === 09/26/94
! ===   Reset year and month index:
!
      mymonth = mymonth + 1
      if(mymonth.gt.12)then
        iyearchem  = iyearchem +1
        mymonth = 1
!     endif ! 27/8/2005

! === 092295
! === write rawdata for possible renew run
! ===   at end of each month:
! ===   at end of each year: 27/8/2005
!
        rewind 178
        print *,'For chem restart ',iyearchem,mymonth
        write(178)iyearchem,mymonth,airmass,
     &              cfc11,cfc110,
     &              cfc11m,
     &              cfc11sd,
     &              cfc12,cfc12m,
     &              cfc12sd,
     &              xn2o ,xn2om ,
     &              xn2osd,
     &              hfc134a,hfc134am,
     &              pfc   ,pfcm,
     &              sf6   ,sf6m,
     &              bcarbon,bcm,
     &              ocarbon,ocm,
     &              atomo ,
     &              o1d   ,
     &              o3    ,o3m  ,
     &              co    ,com  ,
     &              zco2  ,zco2m,
     &              atomh ,
     &              ho    ,
     &              ho2   ,hoxm ,
     &              h2o2  ,
     &              xno   ,
     &              xno2  ,xnoxm,
     &              xno3  ,
     &              xn2o5 ,xnoym,
     &              hno3  ,
     &              ch4   ,ch4m ,
     &              ch3   ,
     &              cho   ,
     &              ch2o  ,
     &              ch3o  ,
     &              ch3o2 ,
     &              ch3o2h,
     &              so2   ,so2m ,
     &              hoso2 ,
     &              so3   ,
     &              h2so4 ,h2so4m,
     &              sviod ,sviodm
        rewind 178

      endif
#endif

#if ( defined CPL_TEM )
!#if ( defined CLM )
      do J=1,JM
c
#ifdef PREDICTED_GASES
        co24tem(j)=zco2(1,j,1)
     &           *28.97296245/44.0*1.e-3
                        !ppm(m) to kg per volume base
#else
         co24tem(j)=atm_co2(j)

#endif
c
      enddo
#endif

      endif	! end of month
 550  continue
C END of DAilY


C     CALL CHECKT (11)                                                   203.17 
      CALL CLOCKS (MNOW)                                                 203.18 
      MSURF=MSURF+(MLAST-MNOW)                                           203.19 
      MLAST=MNOW                                                         203.2  
C****                                                                    204.   
C**** WRITE INFORMATION ONTO A TAPE EVERY USET HOURS                     205.   
C****                                                                    206.   
      IF(USET.LE.0.) GO TO 600                                           207.   
      IF(.NOT.EVENT(USET)) GO TO 600                                     208.   
C     COMPUTATIONS FOR XXXXXX                                            209.   
      WRITE (520) TAU,XXXXXX                                             210.   
      CALL CLOCKS (MNOW)                                                 211.   
      MINC=MLAST-MNOW                                                    212.   
      MELSE=MELSE+MINC                                                   213.   
      PERCNT=100.*MELSE/(MSTART-MNOW+1.E-5)                              214.   
c     WRITE (6,910) MINC,MELSE,PERCNT,TAU                                215.   
C****                                                                    216.   
C**** CALL DIAGNOSTIC ROUTINES                                           217.   
C****                                                                    218.   
  600 IF(MOD(NSTEP-NDYN,NDA4).EQ.0) CALL DIAG4A                          219.   
C                                                                        220.   
      IF(NDPRNT(1).GE.0) GO TO 610                                       221.   
C**** PRINT CURRENT DIAGNOSTICS (INCLUDING THE INITIAL CONDITIONS)       222.   
      IF(KDIAG(1).LT.9) CALL DIAG1(NOCLM)                                223.   
      IF(KDIAG(2).LT.9) CALL DIAG2                                       224.   
      IF(KDIAG(7).LT.9) CALL DIAG7P                                      225.   
      IF(KDIAG(3).LT.9) CALL DIAG3                                       226.   
C                                                                        227.   
C                                                                        228.   
      IF(KDIAG(4).LT.9) CALL DIAG4                                       229.   
      IF(TAU.LE.TAUI+DTHR*(NDYN+.5)) CALL DIAGKN                         230.   
      NDPRNT(1)=NDPRNT(1)+1                                              231.   
C  690 changed to 691 02/21/2003
  610 IF(.NOT.EVENT(24.)) GO TO 691                                      232.   
C**** PRINT DIAGNOSTIC TIME AVERAGED QUANTITIES ON NDPRNT-TH DAY OF RUN  233.   
      DO 620 K=1,13                                                      234.   
      IF (JDAY.EQ.NDPRNT(K)) GO TO 630                                   235.   
  620 CONTINUE                                                           236.   
      GO TO 640                                                          237.   
c 630 WRITE (6,920)                                                      238.   
  630 continue
      IF(KDIAG(1).LT.9) CALL DIAG1(NOCLM)                                239.   
      IF(KDIAG(2).LT.9) CALL DIAG2                                       240.   
      IF(KDIAG(7).LT.9) CALL DIAG7P                                      241.   
      IF(KDIAG(3).LT.9) CALL DIAG3                                       242.   
C                                                                        243.   
C                                                                        244.   
C     IF(KDIAG(6).LT.9) CALL DIAG6                                       245.   
      IF(KDIAG(4).LT.9) CALL DIAG4                                       246.   
      IF(KDIAG(8).LT.9) CALL DIAG8 (IPFLAG)                              247.   
C**** THINGS TO DO BEFORE ZEROING OUT THE ACCUMULATING ARRAYS            248.   
C****   (NORMALLY DONE AT THE END OF A MONTH)                            249.   
  640 DO 650 K=1,13                                                      250.   
      IF(JDAY.EQ.NDZERO(K)) GO TO 660                                    251.   
  650 CONTINUE                                                           252.   
      GO TO 690                                                          253.   
C**** PRINT THE KEY DIAGNOSTICS                                          254.   
  660 CONTINUE                                                           255.   
      CALL DIAGKN                                                        255.   
C**** PRINT AND ZERO OUT THE TIMING NUMBERS                              256.   
      CALL CLOCKS (MNOW)                                                 257.   
      MDIAG=MDIAG+(MLAST-MNOW)                                           258.   
      MLAST=MNOW                                                         259.   
      TOTALT=.01*(MSTART-MNOW)                                           260.   
      PDYN=MDYN/TOTALT                                                   261.   
      PCDNS=MCNDS/TOTALT                                                 262.   
      PRAD=MRAD/TOTALT                                                   263.   
      PSURF=MSURF/TOTALT                                                 264.   
      PDIAG=MDIAG/TOTALT                                                 265.   
      PELSE=MELSE/TOTALT                                                 266.   
      DTIME=24.*TOTALT/(60.*(TAU-TAU0))                                  267.   
c     WRITE (6,909) DTIME,PDYN,PCDNS,PRAD,PSURF,PDIAG,PELSE              268.   
      MDYN=0                                                             269.   
      MCNDS=0                                                            270.   
      MRAD=0                                                             271.   
      MSURF=0                                                            272.   
      MDIAG=0                                                            273.   
      MELSE=0                                                            274.   
      MSTART=MNOW                                                        275.   
      if(TRANSR)then
      do 5284 j=1,JM0
        do 5284 k=1,lmo
          RTGOAV(j,k)=RTGOAV(j,k)/NODIFS
 5284 continue
c     print *,'ATM RTGOAV monthly'
c     print *,(RTGOAV(J,1),j=1,jm)
       endif	! TRANSR
       SPGAV=srps(jm+3)/nsrps+PTOP
       do 5285 j=1,JM+3
        GBUDG(j,38,1)=(srps(j)/nsrps+PTOP)*1013./SPGAV
 5285 continue
c      do 5287 j=1,JM+3
c       GBUDG(j,38,1)=GBUDG(j,37,1)*1013./GBUDG(jm+3,37,1)
c5287 continue
!      print *,'FRMDICE'
!      print '6(1PE12.4)',FRMDICE
       ENKE=0.0
       ENPT=0.0
      do ii=1,4
       ENKE=ENKE+SPECA(1,19,ii)
       ENPT=ENPT+SPECA(1,20,ii)
      enddo
c     print *,'ENKE=',ENKE,' ENTP=',ENPT,' ENTT=',ENKE+ENPT
      print *,'ENKE=',QTABLE(JMP3,19,10),' ENTP=',QTABLE(JMP3,20,10),
     &     ' ENTT=',QTABLE(JMP3,19,10)+QTABLE(JMP3,20,10)
      print *,'WLMMAX=',WLMMAX,' JWLMMAX=',JWLMMAX
c     print *,'AJ(*,37)'
c     print *,(AJ(J,37),j=1,jm)
c     print *,'AJ(*,28)'
c     print *,(AJ(J,28),j=1,jm)
c     IF(USEP.LE.0.) GO TO 680                                           276.   
C**** WRITE SELECTED DIAGNOSTICS ONTO A DISK DATA SET FOR PLOTTING       277.   
c     IF(TAU.LE.TAUI+1080.) GO TO 675                                    278.   
c 670 READ (16) TAUX                                                     279.   
c     IF(TAU.GT.TAUX+1080.) GO TO 670                                    280.   
  675 WRITE (546) AEXP4,JDATE,JMONTH,JYEAR,JDATE0,JMNTH0,JYEAR0,         281.   
     *  GBUDG,QMAPS,QTABLE,INQTAB,J1QT,INQMAP,RTGOAV
      print *,'From atm write(546) ',JMNTH0,JYEAR0,' ',JMONTH,JYEAR
c     print *,0.1*GBUDG(1,26,2),0.1*GBUDG(1,35,2)
      nwr=nwr+1
  690 continue
      if(JDAY.eq.1)then
      rewind 547
      write(547)AEXP,nwr
      rewind 547
      print *,'From atm write(547) ',AEXP,nwr
      endif
      if(wr25.and.JDAY.eq.1)then
c      Write a restart file once a year
       print *,'Write a restart file once a year.'

#if ( defined CPL_OCEANCO2 && defined ML_2D )
C     Data for possible restart for OCM
        write(369)jyear-1,vdfocm
        write(369)Hg
        write(369)Rco2
        rewind 369
#endif

       print *,' KDISK0=',KDISK0
       CALL RFINAL (IRAND)
       REWIND KDISK0                                                    
       if(TRANSR)then
        WRITE(KDISK0) AEXP,TAU,JC,C,RC,KEYNR,U,V,T,P,Q,ODATA,GDATA,
     *   BLDATA,
     *  RQT,SRHR,TRHR,(AJ(K,1),K=1,KACC),TAU,TSSFC,CKS,CKN,WMGE,TPRIM2 
     *  ,MRCHT,TRSURF,SRSURF,TLANDW,TSURFW,DWAV0
     *  ,TG3M,RTGO,STG3,DTG3
       else
        WRITE(KDISK0) AEXP,TAU,JC,C,RC,KEYNR,U,V,T,P,Q,ODATA,GDATA,
     *   BLDATA,
     *  RQT,SRHR,TRHR,(AJ(K,1),K=1,KACC),TAU,TSSFC,CKS,CKN,WMGE,TPRIM2 
     *  ,MRCHT,TRSURF,SRSURF,TLANDW,TSURFW,DWAV0
c       print *,' TSURFT'
c       print 5001,TSURFT
c       print *,' TSURFW'
c       print 5001,TSURFW
       endif
       REWIND KDISK0                                                   
       KDISK=3-KDISK
       KDISK0=500+KDISK
      end if
C  690 changed to 691 02/21/2003
  680    IF(KCOPY.LE.0) GO TO 691                                        284.   
C**** WRITE A COPY OF THE FINAL RESTART DATA SET ONTO DISK               285.   
         CALL RFINAL (IRAND)                                             285.5  
         print *,' after 680'
         print *,' TAU=',TAU,' IRAND=',IRAND
         IF(KCOPY.GT.99) GO TO 687                                       286.   
  683    READ (KCOPY) TAUX                                               286.5  
         IF(TAU.GT.TAUX+3240.) GO TO 683                                 287.   
  685    WRITE (KCOPY) TAU,JC,C,RC,KEYNR,U,V,T,P,Q,ODATA,GDATA,BLDATA,   287.5  
     *   RQT,SRHR,TRHR,(AJ(K,1),K=1,KACC),TAU,TSSFC,CKS,CKN              288.   
     *   ,WMGE,TPRIM2,MRCHT,TRSURF,SRSURF,TLANDW,TSURFW,DWAV0
     *   ,TG3M,RTGO,STG3,DTG3
         REWIND KCOPY                                                    288.5  
C  690 changed to 691 02/21/2003
         GO TO 691                                                       289.   
  687    KCOPY=KCOPY-100                                                 289.5  
         GO TO 685                                                       289.6  
C**** TIME FOR CALLING DIAGNOSTICS                                       290.   
C  690 changed to 691 02/21/2003
  691 CALL CLOCKS (MNOW)                                                 291.   
      MDIAG=MDIAG+(MLAST-MNOW)                                           292.   
      MLAST=MNOW                                                         293.   
  780 IF(TAU.LE.TAUI+DTHR*(NDYN+.5)) GO TO 120                           294.   
      GO TO 100                                                          295.   
C****                                                                    296.   
C**** END OF MAIN LOOP                                                   297.   
C****                                                                    298.   
C**** RUN TERMINATED BECAUSE SENSE SWITCH 6 WAS TURNED ON                299.   
  800 WRITE (6,904)                                                      300.   
      IF(EVENT(TAUT)) GO TO 820                                          301.   
      CALL RFINAL (IRAND)                                                302.   
      print *,' after 800'
      print *,' TAU=',TAU,' IRAND=',IRAND
      if(wr25) then
      REWIND KDISK0                                                      303.   
      WRITE(KDISK0) AEXP,TAU,JC,C,RC,KEYNR,U,V,T,P,Q,ODATA,GDATA,BLDATA, 304.   
     *  RQT,SRHR,TRHR,(AJ(K,1),K=1,KACC),TAU,TSSFC,CKS,CKN,WMGE,TPRIM2   305.   
     *  ,MRCHT,TRSURF,SRSURF,TLANDW,TSURFW,DWAV0
      end if
C     WRITE (6,908)                                                      306.   
C**** RUN TERMINATED BECAUSE IT REACHED TAUE (OR SS6 WAS TURNED ON)      307.   
c 820 WRITE (6,905) TAU,IDAY,TOFDAY                                      308.   
  820 continue

#if ( defined OCEAN_3D || defined ML_2D )
C     DTATM time step of atm model in hours
C     precip and evap in mm/day or kg/m**2/day
         do j=1,jm0
Cjrs #if ( defined OCEAN_3D && defined CPL_OCEANCO2 ) 
#ifdef OCEAN_3D
!jrs           ncallatm=ncallatm+1
! 020107
!          co24ocean(j)=pC_atm(j)*1.e6
! jrs give CO2 even if ocn carbon off
           co24ocean(j)=atm_co2(j)
#  ifdef CPL_OCEANCO2
           co24ocnan(j)=co24ocnan(j)+co24ocean(j)
#  endif
#endif
#ifdef ML_2D
cjrs block only MD_2D
           rseaice(j)=ODATA(1,J,2)
#endif
           tauu(j)=tauu(j)/(NSURF*DTATM)
           tauv(j)=tauv(j)/(NSURF*DTATM)
           tempr(j)=tempr(j)/DTATM
           precip(j)=precip(j)/(DTATM/24.)
           fland=FDATA(1,J,2)
           if (fland.lt.1.0)then
           precip(j)=precip(j)*(1.-fland*prlnd2total(j,mndriver))
     &      /(1.-fland)
           endif
Cjrs            surfpr(j)=surfpr(j)/(DTATM/24.)
c
           if(naveo(j).gt.0)then
             hfluxo(j)=NSURF*hfluxo(j)/(NDYN*DT*naveo(j))
             dhfodtg(j)=dhfodtg(j)/naveo(j)
             dhfodtgeq(j)=dhfodtgeq(j)/naveo(j)
             evao(j)=NSURF*evao(j)/(NDYN*DT*naveo(j))
             devodtg(j)=devodtg(j)/naveo(j)
             devodtgeq(j)=devodtgeq(j)/naveo(j)
C           From mm/sec to mm/day
             evao(j)=24.*3600.*evao(j)
             devodtg(j)=24.*3600.*devodtg(j)
             devodtgeq(j)=24.*3600.*devodtgeq(j)
           endif
C
           if(navei(j).gt.0)then
             hfluxi(j)=NSURF*hfluxi(j)/(NDYN*DT*navei(j))
             dhfidtg(j)=dhfidtg(j)/navei(j)
             dhfidtgeq(j)=dhfidtgeq(j)/navei(j)
             evai(j)=NSURF*evai(j)/(NDYN*DT*navei(j))
             devidtg(j)=devidtg(j)/navei(j)
             devidtgeq(j)=devidtgeq(j)/navei(j)
C           From mm/sec to mm/day
             evai(j)=24.*3600.*evai(j)
             devidtg(j)=24.*3600.*devidtg(j)
             devidtgeq(j)=24.*3600.*devidtgeq(j)
           endif
C
           if(navrad(j).gt.0)then
             solarinc_ice(j)=solarinc_ice(j)/navrad(j)
             solarnet_ice(j)=solarnet_ice(j)/navrad(j)
           endif
           if(navrado(j).gt.0)then
             solarinc_ocean(j)=solarinc_ocean(j)/navrado(j)
             solarnet_ocean(j)=solarnet_ocean(j)/navrado(j)
           endif
c          Runoff is a flux of water from land in mm/day
c          not for m**2
cjrs change runoff to new name arunoff
           arunoff(j)=arunoff(j)/(DTATM/24.)*FDATA(1,j,2)
     &          *DXYP(J)
           if(NWMGEA(J).gt.0)then
             wsocean(J)=WSAV(J)/NWMGEA(J)
             NWMGEA(J)=0.
           else
             wsocean(J)=0.
           end if
           WSAV(J)=0.
c          dhdtav(j)=dhdtav(j)+dhfdtg(j)
c          devdtav(j)=devdtav(j)+devdtg(j)
c
c
           ps4ocean(j)=ps4ocean(j)/DTATM
           do l=1,lm
            qyz4ocean(j,l)=qyz4ocean(j,l)/DTATM
            tyz4ocean(j,l)=tyz4ocean(j,l)/DTATM
           enddo
c
c
        end do    ! j
           rungl=0.0
           runn=0.0
           runt=0.0
           runs=0.0
           SLAND=0.0
           CLAT=20.*TWOPI/360.
           do j=1,jm
             SLAND=SLAND+FDATA(1,j,2)*DXYP(J)
c jrs runoff->arunoff
             rungl=rungl+arunoff(j)
             if(LAT(J).lt.-CLAT)then
               runs=runs+arunoff(j)
             else if(LAT(J).lt.CLAT)then
               runt=runt+arunoff(j)
             else
               runn=runn+arunoff(j)
             endif
          enddo
c          print *,'RUNOFF TOFDAY=',TOFDAY
c          print *,rungl/SLAND,rungl,runs,runt,runn
c       nmonth=JMNTH0
#ifdef ML_2D
c jrs only ML_2D
        nmonth=AMONTH(mndriver)
#endif
        jdatefl=jdate-1
c       if(JMONTH.ne.JMNTH0)jdatefl=NDAYMN(mndriver)
c       if(JMONTH.ne.JMNTH0)then
c       print *,'OCEAN FLUXS'
c       print *,JMONTH,JMNTH0,JDATE
c       print *,'Month=',AMONTH(mndriver),' day=',jdatefl
c       print *,' TAIR_OCEAN'
c       print *,(tairo(j),j=1,jm0)
c       print *,' TAIR_ICE'
c       print *,(tairi(j),j=1,jm0)
c       print *,' TAUU'
c       print *,(tauu(j),j=1,jm0)
c       print *,' TAUV'
c       print *,(tauv(j),j=1,jm0)
c       print *,' EVA-E OCEAN'
c       print *,(evao(j),j=1,jm0)
c       print *,' EVA-I OCEAN'
c       print *,(evai(j),j=1,jm0)
c       print *,' P-E ICE'
c       print *,(pmei(j),j=1,jm0)
c       print *,' HEAT FLUX OCEAN'
c       print *,(hfluxo(j),j=1,jm0)
c       print *,' DHEATFLUX_OCEAN/DTG'
c       print *,(dhfodtg(j),j=1,jm0)
c       print *,' EVA_OCEAN/DTG'
c       print *,(devodtg(j),j=1,jm0)
c       print *,' EVA_ICE/DTG'
c       print *,(devidtg(j),j=1,jm0)
c       print *,' HEAT FLUX ICE'
c       print *,(hfluxi(j),j=1,jm0)
c       print *,' DHEATFLUX_ICE/DTG'
c       print *,(dhfidtg(j),j=1,jm0)
c       print *,' DLH_OCEAN/DTG'
c       print *,(devidtg(j),j=1,jm0)
c       print *,'PS4OCEAN'
c       print *,(ps4ocean(j),j=1,jm0)
c       print *,'QYZ4OCEAN'
c       do l=1,lm
c        print *,(qyz4ocean(j,l),j=1,jm0)
c       enddo
c       endif
c       go to 587
c       write(893),nmonth,jdatefl,tempr,tauu,tauv,precip,evao,
c    &  evai,hfluxo,dhfodtg,devodtg,hfluxi,dhfidtg,devidtg,
c    &  solarinc_ice,solarnet_ice,rseaice
#ifdef ML_2D
c jrs only ML_2D
         do j=1,jm
           osst(j)=ODATA(1,j,1)
           aoice(j)=ODATA(1,j,3)
           foice(j)=ODATA(1,j,2)
           snowice(j)=GDATA(1,j,1)
           tice1(j)=GDATA(1,j,3)
           tice2(j)=GDATA(1,j,7)
          enddo
#endif
c        write (894),nmonth,jdatefl,osst,aoice,foice,snowice,tice1,tice2
 587     continue
#endif

C
c     TAU for coupler
c
      TAUATM=TAU
      MONTHATM=JMONTH
      JDATEATM=JDATE
      JYEARATM=JYEAR
C
#ifdef ML_2D
Cjrs change this block to only ML_2D
         IDAYM=IDAY
         JDAYM=JDAY
         JDATEM=JDATE
         JMONTHM=JMONTH
         JYEARM=JYEAR
         TAUML=TAU
         TOFDAYML=TOFDAY
#endif
C
      if(JDAY.ne.JDAYLAST)then
c      print *,'co24ocean=',co24ocean(jm/2)
       if(JDATEPR.ne.0)then
c      WRITE (6,905) TOFDAY,JDATE,JMONTH,JYEAR
       WRITE (6,905) TOFDAYPR,JDATEPR,JMONTHPR,JYEARPR
        endif
cjrs        print *,'ncallclm=',ncallclm
       JDAYLAST=JDAY
c      if(ncallclm.gt.6) stop
c      stop
      endif
      return
C     CALL ENQJOB                                                        309.   
C     CALL ENQJOB                                                        310.   
      IF(IPFLAG.EQ.0) STOP 13                                            311.   
      STOP 1                                                             312.   
C****                                                                    313.   
  901 FORMAT ('0CLIMATE MODEL STARTED UP',14X,'DAY',I5,', HR',F6.2,I6,   314.   
     *  A5,I27,I7,F7.1,'   TAU',F9.2)                                    315.   
  902 FORMAT ('   DYNAMIC TERMS INTEGRATED, MRCH=',I1,4X,'DAY',I5,       316.   
     *  ', HR',F6.2,I6,A5,2I7,F7.1,23X,'TAU',F9.2)                       317.   
  903 FORMAT (' SOURCE TERMS INTEGRATED',64X,2I7,F7.1)                   318.   
  904 FORMAT ('0SENSE SWITCH 6 HAS BEEN TURNED ON.')                     319.   
  905 FORMAT (/1(1X,33('****')/)/                                        320.   
c    *  ' PROGRAM HAS TERMINATED NORMALLY.  TAU,IDAY,TOFDAY=',F9.2,      321.   
     *  ' ATM HAS TERMINATED NORMALLY. AT ',F9.2,I3,A5,i5, 
     *   /1(1X,33('****')/))                                             322.
  906 FORMAT (' OUTPUT RECORD WRITTEN ON UNIT',I3,55X,2I7,F7.1,          323.   
     *  '   TAU',F9.2,' ON ',A4)                                         324.   
  908 FORMAT (' OUTPUT RECORD WRITTEN ON UNIT',I3,79X,'TAU',F9.2,        325.   
     *  ' ON ',A4)                                                       326.   
  909 FORMAT (/'0TIME',F7.2,'(MINUTES)    DYNAMICS',F5.1,                327.   
     *  '    CONDENSATION',F5.1,'    RADIATION',F5.1,'    SURFACE',      328.   
     *  F5.1,'    DIAGNOSTICS',F5.1,'    OTHER',F5.1//)                  329.   
  910 FORMAT (' INFORMATION WRITTEN ON UNIT 20',57X,2I7,F7.1,            330.   
     *  '   TAU',F9.2,' ON TAPE')                                        331.   
  920 FORMAT ('1'/64(1X/))                                               332.   
      END                                                                333.   
