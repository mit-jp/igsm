C $Header: /home/igsm-repository/ATM/input.F,v 1.3 2010/09/02 22:02:10 jscott Exp $
C $Name:  $

#include "ctrparam.h"

!	==========================================================
!
!	INPUT.F: THIS SUBROUTINE SETS THE PARAMETERS IN THE 
!		   C ARRAY, READS IN THE INITIAL CONDITIONS,
!		   AND CALCULATES THE DISTANCE PROJECTION ARRAYS
!
!	----------------------------------------------------------
!
!	Author of Chemistry Modules:	Chien Wang
!
!	----------------------------------------------------------
! 
!	Revision History:
!	
!	When	Who		What
!	----	----------	-------	
!	073100	Chien Wang	repack based on CliChem3 and add cpp
!
!	==========================================================

      SUBROUTINE INPUT                                                  1501.   
C****                                                                   1502.   
C**** THIS SUBROUTINE SETS THE PARAMETERS IN THE C ARRAY, READS IN THE  1503.   
C**** INITIAL CONDITIONS, AND CALCULATES THE DISTANCE PROJECTION ARRAYS 1504.   
C****                                                                   1505.   

#if ( defined CPL_CHEM )
!
#include "chem_para"
#include "chem_com"
!
#endif

#include "ODIFF.COM"
#include "BD2G04.COM"
#include "RADCOM.COM"
#include "run.COM" 
#include "DRIVER.h"

#if ( defined OCEAN_3D )
#include "AGRID.h"
#endif

cjrs done in driver.h #if ( defined CPL_TEM )
cjrs#include "TEM.h"
cjrs#endif
!
! === Chien Wang 062904
!
	character(100) :: cfname
!
      COMMON/OCN/TG3M(1,JM0,12),RTGO(1,JM0,lmo),STG3(1,JM0),DTG3(1,JM0)
      COMMON/SPEC2/KM,KINC,COEK,C3LAND(IO0,JM0),C3OICE(IO0,JM0)         1506.1  
     *  ,C3LICE(IO0,JM0),WMGE(IO0,JM0),TSSFC(IM0,JM0,4)                 1506.2  
      COMMON/EPARA/VTH(JM0,LM0),WTH(JM0,LM0),VU(JM0,LM0),VV(JM0,LM0)
     &   ,DQSDT(JM0,LM0)                                                1506.3  
     *  ,DWV(JM0),PHIT(JM0,LM0),TPRIM2(JM0,LM0),WU(JM0,LM0),CKS,CKN     1506.4  
     *  ,WQ(JM0,LM0),VQ(JM0,LM0),MRCHT                                  1506.5  
      CHARACTER*4 XLABL1                                                1506.6  
      COMMON U,V,T,P,Q                                                  1507.   
C     COMMON/KEYS/KEYNR(42,50)                                          1508.   
c      COMMON/RADCOM/VADATA(11,4,3)
c     CHARACTER*8 RECORD,ANDEND,NLREC*80                                1510.   
      CHARACTER*12 RECORD,ANDEND,NLREC*80                               1510.   
      CHARACTER*4 C,C1,NAMD60,DISK,RUNID                                1510.1  
      CHARACTER*5 TSCNTR
      DIMENSION RECORD(10)                                              1510.2  
      DIMENSION JC(100),C(39),RC(161),JC1(100),C1(39),RC1(161)          1511.   
      EQUIVALENCE (JC(1),IM),(C(1),XLABEL(1)),(RC(1),TAU)               1511.1  
      DIMENSION IDAYS0(13),NAMD60(4),SIG0(36),SIGE0(37)                 1512.   
      DIMENSION XA(1,JM0),XB(1,JM0),XLABL1(33)                          1512.5  
      DIMENSION JDOFM(13),VMASK(JM0)
     & ,DSIGF(LM0),DSIGH(LM0)
      character *120 file1,file2,plotfl,nwrfl
      character * 120  t3file,tsfile,zmfile,qffile,clfile,wrcldf
     &	  ,ochemfile,deepco2in
cjrs DRIVER.h caruptfile  & ,oco2file,co2rfile,caruptfile,flrco2av
     &	  ,oco2file,co2rfile,flrco2av
     &    ,ghg_monthly,ghg_monthly2,co2_data,o3_data
     &    ,bgrghg_data
      character * 120 sulf1986,sulf2050,sulfamp,SO2_EM,
     &  S0C_data,
     &  dirdat1,dirdat2
     &  ,bc_data
      character * 120 chemdata,chemout,init_4nem,pov_deepo
cjrs DRIVER.h last_nep character * 120 chemdata,chemout,last_nep,init_4nem,pov_deepo
     & ,flin_nep,last_clm,SEN_dat
cjrs DRIVER.h fnememiss &  ,fl_init_alkt,fl_init_salt,fl_dic_eq,fnememiss,
     &  ,fl_init_alkt,fl_init_salt,fl_dic_eq,
     & chem_init,chem_init2,chem_init4pfc,chemrstfl
      common/files/file1,file2,plotfl,nwrfl,qffile,clfile,wrcldf
     *,t3file,tsfile,zmfile,ochemfile,deepco2in
      character * 120 ghostfile
c==== 012201
      common/nemdata/nemdatdir
      character *120 nemdatdir

      common /bmtrdata/co2_data
      common /bghgdata/bgrghg_data
      common /sulfdata/sulf1986,sulf2050,sulfamp,SO2_EM
      common /o3data/o3_data
      common /solardata/S0C_data
      common/aexpc/AEXP,ISTRT1,ISTRTCHEM
      common/qfl/QFLUX(JM0,0:13),ZOAV(JM0),QFLUXT(JM0)
      common/fixcld/cldssm(JM0,LM0,0:13),cldmcm(JM0,LM0,0:13)
     &  ,CLDSST(JM0,LM0),CLDMCT(JM0,LM0)
      common/SURRAD/TRSURF(JM0,4),SRSURF(JM0,4)
      NAMELIST/INPUTZ/IM,JM,LM,LS1,LBLM,LMCM,LSSM,KOCEAN,ISTART,KDISK,  1513.   
     *  TAUP,TAUNI,TAUE,TAUT,TAUO,NDYN,NCNDS,NRAD,NSURF,NGRND,NFILTR,   1514.   
     *  NDAA,NDA5D,NDA5K,NDA5S,NDA4,NDASF,DT,TAU,XINT,INYEAR,           1515.   
     *  LHE,LHM,LHS,RADIUS,GRAV,RGAS,KAPA,OMEGA,RHMAX,ETA,S0X,CO2,SRCOR,1516.   
     *  PTOP,PSF,PSL,PTRUNC,DLAT,DLON,AREAG,IRAND,IJRA,MFILTR,NDIFS,    1517.   
     *  KACC,KEYCT,SKIPSE,USESLP,USEP,USET,KCOPY,DUMMY1,IDACC,KDIAG,    1518.   
     *  NDZERO,NDPRNT,IJD6,NAMD6,SIG,SIGE,                              1519.   
     *  KM,KINC,COEK,INDAY,IMONTH,LDAY,LMONTH,LYEAR,AEXP,               1519.5  
     *  READGHG,wr25,LFR,ISTRT1,PCLOUD,QFCOR,TRANSR,WRCLD,NWRCLD,CONTRR,
     *  ISTWRC,CLDFEED,OBSFOR,ALFFOR,YEARGT,CO2IN,ISTRTCHEM,
     *  AERFOR,AERF4BC,
     *  S0RATE,CFS0X,
     *  CFAEROSOL,CFBC,
     &  cfocdif,rkv,diffcar0,ocarcont,ocarindata,
!       Kvc=diffcar0+cfocdif*Kvh
     *  file1,file2,plotfl,nwrfl,qffile,clfile,wrcldf,clmsen,cfdif0,
     *  t3file,tsfile,zmfile,ochemfile,deepco2in,
     &  fl_init_alkt,fl_init_salt,fl_dic_eq,
     *  ghg_monthly,ghg_monthly2,co2_data,o3_data,
     &    bgrghg_data,
     *  sulf1986,sulf2050,sulfamp,SO2_EM,
     &  S0C_data,cfvolaer,
     &  dirdat1,dirdat2
     &  ,chemdata,chemout,last_nep,init_4nem,pov_deepo,fnememiss
     &  ,chem_init,chem_init2,chem_init4pfc,chemrstfl
     &	  ,oco2file,co2rfile,caruptfile,flrco2av
     &  ,flin_nep,last_clm,SEN_dat,nemdatdir
     &  ,GHSFALB,GHSF,ALBCF,FVOLADD,ghostfile,fl_volaer
     &  ,STRARFOR,GSOEQ,CO2FOR,CO2F,FORSULF,FORBC,S0FOR,FORVOL
     &  ,VEGCH,vegfile,TRVEG
     &  ,fclmlice,fbaresoil,fwmax,fprratio,o3datadir,CLIMO3
     &  ,OCNGEOM,ocngmfile,ocndata4atm
     &  ,bc_data
     &  ,PRTREND
Cjrs     &  ,dtatm,dtocn
#if ( defined IPCC_EMI )
     &  ,init_co2
      character * 120 init_co2
#endif
      character * 120 ocngmfile,ocndata4atm,fl_volaer
      character * 120 vegfile
     &  ,fclmlice,fbaresoil,fwmax,fprratio,o3datadir
      common/wrcom/wr25,TRANSR,CONTRR,OBSFOR
c jrs      common/TIMESTEPS/dtatm,dtocn
      LOGICAL LFR,NLFR,wr25,TRANSR,WRCLD,CONTRR,CLDFEED,OBSFOR
     &,GHSF,VEGCH,TRVEG,GSOEQ,OCNGEOM,GHSFALB,STRARFOR,CO2FOR
     & ,FORSULF,FORBC,S0FOR,FORVOL
      common/FORAERSOL/FORSULF,FORBC,FORVOL
      common/vaerosol/fl_volaer
      common/eqgso/GSOEQ
      common/ghstfor/GHSFALB,GHSF,ALBCF,FVOLADD,STRARFOR,S0FOR,CO2FOR,
     &  CO2F,ghostfv(LM0+1),ghostf(LM0+1,JM0)
      common/veg/TRVEG,IYVEG
      common/COMCLD/READGHG,PCLOUD,WRCLD,NWRCLD,NWRCL,INYEAR,JNDAY
     &,CFAEROSOL,ALFA,CFBC,cfvolaer
      common/ BACKGRGHG/GHGBGR(5)
      COMMON/CO2TRND/ALFFOR,CO2TR,YEARGT,CO2IN,INYRAD
      common/ S0XR/S0RATE,CFS0X
      common/cldfdb/coefcl(3),CLDFEED,SEN_dat
      common/diff/cfdiff,rkv
      common/Dscale/DWAV0(JM0)
      dimension fland_temp(jm0)
      common/atmos_lo/fland_atm(jm0)  !jrs not sure this does anything

#if ( defined CLM )
cjrs alreadyin DRIVER.h#include "CLM.COM"
      dimension clmlice(jm0),baresoil(jm0),
     &  w1maxclm(jm0),w2maxclm(jm0),vmaskclm(jm0)
      character * 120 lineclm
#endif
#if ( defined CPL_OCEANCO2 && defined ML_2D )
      common/Garyclim/tggary(jm0),wsgary(jm0),areaml(jm0),arsrf(jm0)
!     common/Garydiff/depthml(jm0),edzon(jm0,lmo),dzg(lmo),dzog(lmo-1),
!    &Rco2(jm0,lmo),edohd(lmo),zg(lmo),focean(jm0)
      common/Garydiff/depthml(jm0),edzon(jm0),dzg(lmo),dzog(lmo-1),
     &Rco2(jm0,lmo),edohd(lmo),zg(lmo),focean(jm0)
      common /Garychem/Hg(jm0)
      common /Garyvdif/iyearocm,vdfocm,acvdfc,cfocdif,diffcar0
      common /Garyvlog/odifcarbon,ocarcont
      real  Rco2in(jm0,lmo),Hgin(jm0)
      logical odifcarbon,ocarcont,ocarindata
#endif

      common /SO2EMIN/SO2EM
#if (!defined PREDICTED_GASES)                                               
#if (defined CPL_TEM || defined CPL_OCEANCO2 )
      common /ATCO2/atm_co2(jm0),oco2file,co2rfile
#endif
#endif

      integer PCLOUD
!     common/TSUR/TSURFC(JM0,0:13),TLANDD(JM0),TSURFD(JM0),DTSURF(JM0)
!     common/TSLD/TLANDC(JM0,0:13),TLANDT(JM0),TLANDD(JM0),DTLAND(JM0)
#include "TSRF.COM"
      DATA DISK/'DISK'/,ANDEND/' &END   '/                              1520.   
      DATA IDAYS0/0,1,32,60,91,121,152,182,213,244,274,305,335/         1521.   
      DATA NAMD60/'AUSD','MWST','SAHL','EPAC'/                          1522.   
      DATA EDPERD/1./,EDPERY/365./                                      1527.   
      DATA JDOFM/0,31,59,90,120,151,181,212,243,273,304,334,365/
      DATA JDPERY/365/
!     dimension GHGBGR1860(5),GHGBGR1958(5),GHGBGR1977(5),
!    &  GHGBGR1980(5), GHGBGR1765(5),GHGBGR1991(5),GHGBGR2000(5)
!    &  ,GHGBGR1990(5)
!     DATA GHGBGR1765/280.0,0.275,0.791,0.00E-6,00.0E-6/
!     DATA GHGBGR1860/286.4,0.276,0.805,0.00E-6,00.0E-6/
!     DATA GHGBGR1958/314.9,0.291,1.224,7.60E-6,29.6E-6/
!     DATA GHGBGR1977/331.8,0.292,1.613,13.2e-6,25.2e-6/
!     DATA GHGBGR1980/337.9,0.301,1.547,166.6e-6,300.0e-6/
!     DATA GHGBGR1990/351.0,0.308,1.67,500.6e-6,470.0e-6/
c     DATA GHGBGR1991/355.7,0.310,1.704,268.6e-6,492.0e-6/
C     average for years 1990 and 1991 from  ghgdata.GISS.modified.dat
!     DATA GHGBGR1991/352.7,0.3091,1.681,5.21E-04,4.855E-04/
!     average for years 2000 and 2001 from  ghgdata.GISS.dat
!     DATA GHGBGR2000/368.75,0.316,1.735,5.825E-04,5.35E-04/


C                                                                       1527.5  
C     DEFINE THE COMMON BLOCK /SPEC2/                                   1527.51 
C                                                                       1527.52 
      KM=1                                                              1527.53 
      KINC=1                                                            1527.54 
      COEK=2.                                                           1527.55 
C                                                                       1527.56 
C**** SET PARAMETER DEFAULTS                                            1528.   
      DO 10 K=1,100                                                     1529.   
   10 JC(K)=0                                                           1530.   
      DO 15 K=1,161                                                     1530.1  
   15 RC(K)=0.                                                          1530.2  
      JM=JM0                                                            1531.   
      IM=IM0
      IO=IO0                                                            1532.5  
      LM=LM0                                                            1533.   
C     LS1 is a lowest stratospheric layer, thus LS1=8 means that
C     there are 2 layers in the strosphere k=8 and 9 (for LM=9)
C     and four for LM=11.
C     and four for LM=11.
      LS1=8                                                             1534.   
C
      LBLM=2                                                            1535.   
      ISTART=10                                                         1536.   
      CONTRR=.false.
      CLDFEED=.false.
      CFAEROSOL=1.0
      CFBC=1.0
      cfvolaer=1.0
      ALFFOR=0.
      S0RATE=0.
      YEARGT=1958.
      CO2IN=280.
      INYRAD=1.
      coefcl(1)=0.
      coefcl(2)=0.
      coefcl(3)=0.
      cfdiff=1.
      rkv=1.
      GHSF=.false.
      GHSFALB=.false.
      STRARFOR=.false.
      FORSULF=.true.
      FORBC=.false.
      FORVOL=.false.
      CO2FOR=.false.
      S0FOR=.false.
      CO2F=2.0
      GSOEQ=.false.
      VEGCH=.false.
      TRVEG=.false.
      PRTREND=.false.
      OCNGEOM=.false.
      do L=1,LM+1
        ghostfv(L)=0.
        do j=1,jm
         ghostf(l,j)=0.0
        enddo
      enddo
      cfocdif=0.42
      cfocdif=1.375
      cfocdif=0.6
      cfocdif=3.0
      diffcar0=2.85
      diffcar0=1.00
cc    ALFA = 8.0*1.e3
      AERF4BC=-0.35
      TRANSR=.false.
      WRCLD=.FALSE.
      clfile = 'undefined'
      zmfile = 'undefined'
      t3file = 'undefined'
      tsfile = 'undefined'
      qffile = 'undefined'
      wrcldf = 'undefined'
      oco2file = 'undefined'
      co2rfile = 'undefined'
      fl_volaer = 'undefined'
      SO2_EM = 'undefined'
      sulf1986 = 'undefined'
      sulf2050 = 'undefined'
      sulfamp = 'undefined'
      caruptfile = 'undefined'
      flrco2av = 'undefined'
      ISTRT1=0
      ISTRTCHEM=0
      ISTWRC=0
      QFCOR=.FALSE.
      READGHG=0.
cjrs      dtatm=1
cjrs      dtocn=1
      ghg_monthly =  'undefined'
      ghg_monthly2 =  'undefined'
      co2_data =  'undefined'
      bgrghg_data =  'undefined'
      bc_data =  'undefined'
      S0C_data =  'undefined'
      o3_data =  'undefined'
      ochemfile =  'undefined'
      deepco2in =  'undefined'
      fl_init_alkt = 'undefined'
      fl_init_salt = 'undefined'
      fl_dic_eq = 'undefined'
      ocarcont=.true.
      ocarindata=.false.
      vegfile =  'undefined'
      fclmlice =  'undefined'
      fprratio =  'undefined'
      fbaresoil =  'undefined'
      o3datadir =  'undefined'
      CLIMO3=.false.
      fwmax =  'undefined'
      ocngmfile =  'undefined'
      ocndata4atm =  'undefined'
      ghostfile =  'undefined'
      NWRCLD=0
      TAUNI=0.
      LFR=.TRUE.
      wr25=.true.
      chemout = 'DUMP'
      chemdata = 'DATA'
      chem_init = 'init-data_46x11_1991'
      chem_init2 = 'init-data2_46x11_1991'
      chem_init4pfc = 'undefined'
      chemrstfl = 'undefined'
      nemdatdir = 'TEMDATA'
      last_nep = 'undefined'
      last_clm = 'undefined'
      flin_nep = 'undefined'
      init_4nem = 'undefined'
      fnememiss = 'undefined'
      pov_deepo = 'undefined'
      SKIPSE=1.
      TAUT=6.                                                           1537.   
      TAUT=24.
      KOCEAN=1                                                          1537.1  
      KDISK=1                                                           1538.   
      DT=900.                                                           1539.   
      DT=1200.
      XINT=120.                                                         1540.   
      XINT=24.
      NDYN=4                                                            1541.   
      NDYN=3
#if ( defined CLM )
      NSURF=1
#else
      NSURF=2                                                           1542.   
#endif
!     NSURF=1	! 07/17/2006
      NGRND=1                                                           1543.   
      TAUP=-1.                                                          1544.   
      TAUI=-1.                                                          1545.   
      TAUE=1.E30                                                        1546.   
      TAUO=1.E30                                                        1547.   
      IYEAR=1976                                                        1548.   
      TWOPI=8.*atan(1.)
      TWOPI=6.283185                                                    1549.   
      SDAY=86400.                                                       1550.   
      LHE=2500000.                                                      1551.   
      LHM=334000.                                                       1552.   
      LHS=2834000.                                                      1553.   
      RADIUS=6375000.                                                   1554.   
      GRAV=9.81                                                         1555.   
      RGAS=287.                                                         1556.   
      KAPA=.286                                                         1557.   
      PTOP=10.                                                          1558.   
      PSF=984.                                                          1559.   
      PSL=1000.                                                         1560.   
      PTRUNC=1./8192.                                                   1561.   
      S0X=1.                                                            1561.1  
      CFS0X=1.   
      CO2=1.                                                            1561.2  
      SRCOR=1.                                                          1561.3  
      ED=.1                                                             1562.   
      EDM=.1                                                            1563.   
      ETA=0.                                                            1564.   
      ETA=1.
      RHMAX=100000.                                                     1565.   
      RHMAX=150000.
      CDX=1.                                                            1566.   
      IRAND=123456789                                                   1567.   
      IJRA=1                                                            1568.   
      MFILTR=1                                                          1569.   
      MFILTR=2
      KEYCT=1                                                           1570.   
      CKN=1.00                                                          1570.5  
      CKS=1.00                                                          1570.6  
      DUMMY1(1)=.0005                                                   1571.   
      DUMMY1(2)=.00005                                                  1572.   
      IJD6(1,1)=32                                                      1573.   
      IJD6(2,1)=9                                                       1574.   
      IJD6(1,2)=9                                                       1575.   
      IJD6(2,2)=18                                                      1576.   
      IJD6(1,3)=19                                                      1577.   
      IJD6(2,3)=14                                                      1578.   
      IJD6(1,4)=7                                                       1579.   
      IJD6(2,4)=12                                                      1580.   
      DO 20 KR=1,4                                                      1581.   
   20 NAMD6(KR)=NAMD60(KR)                                              1582.   
      NDPRNT(1)=-1                                                      1583.   
      DO 30 K=2,13                                                      1584.   
      NDPRNT(K)=IDAYS0(K)                                               1585.   
   30 NDZERO(K)=IDAYS0(K)                                               1586.   
      DO 40 K=1,50                                                      1587.   
   40 KEYNR(2,K)=0                                                      1588.   
      DO 45 K=1,12                                                      1589.   
   45 KDIAG(K)=10                                                       1590.   
      KDIAG(1)=0.
      KDIAG(2)=0.
c     KDIAG(3)=0.
        DST=.02053388
        DSB=1.-.948665
        if(LM.eq.9)then
          CALL MESH09(LM,DST,DSB,SIG0,SIGE0,DSIGH,DSIGF)
	elseif(LM.eq.11)then
          CALL MESH11(LM,DST,DSB,SIG0,SIGE0,DSIGH,DSIGF)
	else
	 print *,' wrong LM LM=',LM
	 stop
	endif
      DO  L=1,LM 
	SIGE(l)=SIGE0(LM+2-l)
	SIG(l)=SIG0(LM+1-l)
      ENDDO
c     DO 50 L=1,LM                                                      1591.   
c     SIG(L)=SIG0(L)                                                    1592.   
c  50 SIGE(L)=SIGE0(L)                                                  1593.   
      SIGE(LM+1)=0.                                                     1594.   
      WRITE (6,901)                                                     1595.   
      open(535,file='name.dat')
      READ (535,902) XLABEL                                             1596.   
      open(514,file='name.tmp')
      DO 51 I=1,33                                                      1596.1  
   51 XLABL1(I)=XLABEL(I)                                               1596.2  
      XLABEL(33)=DISK                                                   1597.   
      WRITE (6,903) XLABEL                                              1598.   
C**** COPY INPUTZ NAMELIST ONTO CORE TAPE AND TITLE PAGE                1599.   
   60 CONTINUE
      READ (535,904) RECORD                                             1600.   
      WRITE (514,904) RECORD                                            1601.   
      WRITE (6,905) RECORD                                              1602.   
      IF(RECORD(1).NE.ANDEND) GO TO 60                                  1603.   
      rewind 514
      read (UNIT=514,NML=INPUTZ)
      REWIND 514                                                        1606.   
C JRS ignore name.dat values, start Jan 1. with couple.nml years  
      inyear = startYear
      lyear = endYear +1
      INDAY    = 1
      IMONTH   = 1
      LDAY     = 1
      LMONTH   = 1
Cjrs      dtatmo=dtatm
cjrs      dtocno=dtocn
cb    open statments
c
c     File which depend on resolution 
      id2=index(dirdat2," ")
c
      open( unit=519,file=dirdat2(1:id2-1)//'FILE19',
     *    status='OLD',form='unformatted')
      open( unit=515,file=dirdat2(1:id2-1)//'FILE15',
     *    status='OLD',form='unformatted')
      open( unit=523,file=dirdat2(1:id2-1)//'FILE23',
     *    status='OLD',form='unformatted')
      open( unit=526,file=dirdat2(1:id2-1)//'FILE26',
     *    status='OLD',form='unformatted')
      open( unit=562,file=dirdat2(1:id2-1)//'FILE62',
     *    status='OLD',form='unformatted')
c
      if(VEGCH.or.TRVEG)then
        close(523)
        open( unit=523,file=vegfile,
     &    status='OLD',form='unformatted')
      endif
c     File which do not depend on resolution 
      id1=index(dirdat1," ")
c
      open( unit=509,file=dirdat1(1:id1-1)//'FILE09',
     * status='OLD',form='unformatted')
      open( unit=507,file=dirdat1(1:id1-1)//'FILE07',
     * status='OLD',form='unformatted')
      open( unit=516,file=dirdat1(1:id1-1)//'FILE16',
     * status='OLD',form='unformatted')
      open( unit=517,file=dirdat1(1:id1-1)//'FILE17',
     * status='OLD',form='unformatted')
      open( unit=522,file=dirdat1(1:id1-1)//'FILE22',
     * status='OLD',form='unformatted')
      open( unit=521,file=dirdat1(1:id1-1)//'FILE21',
     * status='OLD',form='unformatted')
c
      if(LMO.eq.12) then
        open( unit=593,file=dirdat1(1:id1-1)//'FOCEAN_12',
     *   status='OLD',form='unformatted')
      endif
c
      if(GHSF)then
        print *,ghostfile
        open(unit=599,file=ghostfile,
     &   status='OLD',form='unformatted')
         read(599) ghostf
         close(599)
      endif
#ifndef CPL_CHEM 
#ifdef PREDICTED_BC
      if(READGHG.eq.0)then
!      data for BC only
       print *,'Data for black carbon'
       open(769,file=bc_data,
     &  status='old',form='unformatted')
      endif
#endif
#endif

#ifdef CPL_CHEM
!
!Chemistry-model
!
! 1. Open initial atmospheric condition to start simulation:

! 2. Open a given data set to initialize simulation:
      print *,' From assign',chemdata,chemout
      id2=index(chemdata," ")
      open(120,file=chemdata(1:id2-1)//'chem.nml',status='old')
C     open(120 moved from chemint.F
C     file chem.nml moved to chemdata

#if ( N_LAT == 24 ) && ( N_LEV == 9 )
      open(163,file=chemdata(1:id2-1)//'init-data_24x9',
     &        form='unformatted',
     &        status='old')

      open(160,file=chemdata(1:id2-1)//'init-data2_24x9',
     &        form='unformatted',
     &        status='old')

      open(162,file=chemdata(1:id2-1)//'topo3_24x9',
     &        form='unformatted',
     &        status='old')
#endif
#if ( N_LAT == 46 ) && ( N_LEV == 11 )
      id1=index(chem_init," ")
      print *,'INITIAL DATA FOR CHEMISTRY'
      print *,chem_init
c     open(163,file=chemdata(1:id2-1)//'init-data_46x11',
      open(163,file=chemdata(1:id2-1)//chem_init(1:id1-1),
     &        form='unformatted',
     &        status='old')

      id1=index(chem_init2," ")
      print *,chem_init2
c     open(160,file=chemdata(1:id2-1)//'init-data2_46x11',
      open(160,file=chemdata(1:id2-1)//chem_init2(1:id1-1),
     &        form='unformatted',
     &        status='old')

!     New init file hfc, pfc, and sf6 for Dec 1990
      id1=index(chem_init4pfc," ")
      print *,chem_init2
      open(161,file=chemdata(1:id2-1)//chem_init4pfc(1:id1-1),
     &        form='unformatted',
     &        status='old')

      open(162,file=chemdata(1:id2-1)//'topo3_46x11',
     &        form='unformatted',
     &        status='old')
#endif
c     open(164,file=chemdata(1:id2-1)//'edaily.dat',
c    &        form='unformatted',
c    &        status='old')

!     open(164,file=emiss_data,
!    &        form='unformatted',
!    &        status='old')

!     open(664,file=SO2ERATIO,
!    &        form='formatted',
!    &        status='old')
!      moved to input.F

      open(165,file=chemdata(1:id2-1)//'rktable.dat',
     &        form='unformatted',
     &        status='old')

! 3. Open a data file to save rawdata:

c closed 032697:

c      open(166,file='OUTPUT/cfcn2o',
c     &        form='unformatted',
c     &        status='new')
c      open(171,file='OUTPUT/hox',
c     &        form='unformatted',
c     &        status='new')
c      open(172,file='OUTPUT/nox',
c     &        form='unformatted',
c     &        status='new')
c      open(173,file='OUTPUT/hcs',
c     &        form='unformatted',
c     &        status='new')
c      open(174,file='OUTPUT/sox',
c     &        form='unformatted',
c     &        status='new')

      cfname = chemout(1:len_trim(chemout))//'/'//'globalmean'
      open(176,file=cfname,
     &        status='new')

      cfname = chemout(1:len_trim(chemout))//'/'//'globalmean2'
      open(180,file=cfname,
     &         status='new')
     
      cfname = chemout(1:len_trim(chemout))//'/'//'meta'
      open(181,file=cfname,
     &        form='unformatted',
     &        status='new')

      cfname = chemout(1:len_trim(chemout))//'/'//'budget'
      open(167, file=cfname,status='unknown')

      cfname = chemout(1:len_trim(chemout))//'/'//'monthly'
      open(169,file=cfname,
     &        form='unformatted',
     &        status='new')

      cfname = chemout(1:len_trim(chemout))
     &          //'/'//'monthly2'
      open(179,file=cfname,
     &         form='unformatted',
     &         status='new')
! ===

      cfname = chemout(1:len_trim(chemout))//'/'//'photodiag'
      open(177,file=cfname,	!032697
     &        form='unformatted',
     &        status='new')

c 5. Open a data file to save diagnostic output:

!      open(168,file='OUTPUT/cfc11nsf_monthly',
!     &        form='unformatted',
!     &        status='new')

        cfname = chemout(1:len_trim(chemout))
     &           //'/'//'chemrenew.dat'
	open(178,file=cfname,
     &        form='unformatted',
     &        status='unknown')

        if(ISTRTCHEM.eq.1)then

C  open statment moved from cheminit.F file number change to 279
        cfname = chemrstfl(1:len_trim(chemrstfl))
          open(279,file=cfname,
     &             form='unformatted',
     &             status='unknown')

        endif

#endif

c     open file for carbon uptake
#if ( defined CPL_TEM || defined CPL_OCEANCO2 )
c     open(333,file=caruptfile,status='new',form='formatted')
      open(333,file=caruptfile,status='replace',form='formatted')
      close(333)
#endif

#if ( defined CPL_OCEANCO2 && defined ML_2D)
      open(668,file=fl_init_alkt,
     &   form='unformatted',status='old')
      open(669,file=fl_init_salt,
     &   form='unformatted',status='old')
!     open(670,file=fl_dic_eq,
!    &   form='unformatted',status='old')
      open(602,file=flrco2av,status='new',form='unformatted')
#endif

ce    open statments
      call bgrghg(YEARGT)
#if ( defined IPCC_EMI )
!     if(YEARGT.eq.1765)then
!      GHGBGR(1)=277.6
       open (unit=861,file=init_co2,
     &     status='OLD',form='formatted')
       read (861,*)xco2init
       CO2=xco2init/GHGBGR(1)
       print *,'IPCC EMI CO2=',CO2
!     else
!        print *,' Wrong YEARGT ', YEARGT
!        stop
!     endif
#endif
      print *,'Background GHGs for year ',YEARGT
      print '(5E12.4)',GHGBGR
      if(CLDFEED)then
C     Calculate coefcl for given clmsen
        print *,'Climate sensitivity=',abs(clmsen)
        call senint(abs(clmsen))
C       NEW
        if(clmsen.gt.0.0)then
         print *,'With coefficients of different signs for clouds '
         print *, '        of diffrent types'
         coefcl(2)=-coefcl(2)
         coefcl(3)=-coefcl(3)
        else
C       OLD
         print *,'With the same coefficient for clouds '
         print *, '        of all types'
        endif
        print *,'coefcl=',coefcl
      else
        print *,'No additional cloud feedback'
        print *,'coefcl=',coefcl
      endif
      if(TRANSR)then
        cfdiff=cfdif0/2.5
        print *,'cfdiff=',cfdiff
       print *,' Weight for old diffusion coefficeints=',rkv
       print *,' Weight for new diffusion coefficeints=',1.-rkv
      else
        print *,'No diffusion into deep ocean'
      endif
      if (PRTREND) then
         print *,' With changes in  stochastic precip'
      else
         print *,' Without changes in  stochastic precip'
      endif
#if ( defined PREDICTED_AEROSOL )
!#if ( defined CPL_CHEM )
Cold  AFBYCF=0.6725
Cold  SO2EREF=123.57
Cigsm1AFBYCF=0.6054
Cigsm1SO2EREF=135.272
!     AFBYCF=1.101897
!     SO2EREF=152.3631
!      read(664,'(f10.6)')SO2EM
!      SO2ER=SO2EM/SO2EREF
!      CFAEROSOL=-AERFOR/(AFBYCF*SO2ER)
Cold  CFAEROSOL=(-AERFOR/AFBYCF)**1.035/(SO2ER**0.7248)
!      CFAEROSOL=(-AERFOR/AFBYCF)**1.035/(SO2ER**1.0391)


      SO2EREF=147.375
!     open(664,file=SO2ERATIO,
!    &        form='formatted',
!    &        status='old')
!     read(664,'(f10.6)')SO2EM
       print *,' SO2EM from input ',SO2EM
      SO2ER=SO2EM/SO2EREF
!     F90BYF80=0.948
!     AERFOR90=AERFOR*F90BYF80
!     AFBYCF90=1.04 
!     print *,'AFBYCF90=',AFBYCF90
!     12/21/2006 (from runs 265x.06)
!     Adjustment for BC forcing
       print *,' AERF4BC=',AERF4BC
      AERFOR=AERFOR+AERF4BC
!     Adjustment for BC forcing
      AFBYCF=1.05
      print *,'SO2ER=',SO2ER
      print *,'AFBYCF=',AFBYCF
      CFAEROSOL=(-AERFOR/AFBYCF)**1.21
      CFAEROSOL=CFAEROSOL/(SO2ER**1.01)
      print *,'AERFOR=',AERFOR,'CFAEROSOL=',CFAEROSOL
#if ( defined OCEAN_3D )
       CFAEROSOL=CFAEROSOL/1.35
       print *,'CFAEROSOL_3D=',CFAEROSOL
#endif
#endif
#if ( defined SVI_ALBEDO )
      ALFA=-16.7*AERFOR*1.e3
      print *,'AERFOR=',AERFOR,' ALFA=',ALFA
#endif
C**** SET DEPENDENT QUANTITIES                                          1608.   
   80 DLON=TWOPI/IM                                                     1609.   
      DLAT=.5*TWOPI/(JM-1)                                              1610.   
      JMM1=JM-1                                                         1611.   
      FIM=IM                                                            1612.   
      FIO=IO                                                            1612.5  
      LMM1=LM-1                                                         1613.   
      LMP1=LM+1                                                         1614.   
      LTM=LS1-1                                                         1615.   
      LSSM=LM                                                           1616.   
      LMCM=LTM                                                          1617.   
c     LMCM=LTM+2
      NCNDS=NDYN                                                        1618.   
      NRAD=5*NDYN                                                       1619.   
      NFILTR=2*NDYN                                                     1620.   
      NFILTR=0
      NDAA=3*NDYN+2                                                     1621.   
      NDA5D=NDYN                                                        1622.   
      NDA5K=NDAA                                                        1623.   
	ndaa=3
      NDA5S=3*NDYN                                                      1624.   
      NDA4=24*NDYN                                                      1625.   
      NDASF=2*NSURF-1                                                   1626.   
      KACC=KACC0
      IF(SKIPSE.GE.1.) KACC=KACC-IM*JM*LM*3+6                           1630.   
	print *,' KACC0=',KACC0,' KACC=',KACC
      IF(ISTART.GE.4) GO TO 90 
      KACC=JM*80*3 + JM*80 + JM*3 + JM*LM*59 + JM*3*4 + IM*JM*75        1627.   
     *  + IM*LM*16 + IM*JM*LM*3 + 20*100 + JM*36 + (IM/2+1)*20*8 + 8*2  1628.   
     *  +24*50*4 + 2                                                    1629.   
      IF(SKIPSE.GE.1.) KACC=KACC-IM*JM*LM*3+6                           1630.   
	print *,' KACC=',KACC
   90 continue
#if( !defined OCEAN_3D)
      open( unit=525,file=zmfile,
     * status='OLD',form='unformatted')
#endif
      if(ISTART.eq.2)then
          open( unit=501,file=file1,
     *    status='new',form='unformatted')
          open( unit=502,file=file2,
     *    status='new',form='unformatted')
          open( unit=546,file=plotfl,
     *    status='new',form='unformatted')
          open( unit=547,file=nwrfl,
     *    status='new',form='unformatted')
      elseif(ISTART.eq.10)then
          open( unit=501,file=file1,
     *     status='OLD',form='unformatted')
          open( unit=502,file=file2,
     *     status='OLD',form='unformatted')
        if(ISTRT1.eq.0)then
          open( unit=546,file=plotfl,
     *    status='new',form='unformatted')
          open( unit=547,file=nwrfl,
     *    status='new',form='unformatted')
#if ( defined CPL_CHEM) && ( defined CPL_TEM )
          open( unit=537,file=flin_nep,
     *    status='OLD',form='unformatted')
#endif
#if ( defined CPL_TEM )
          open (876,file=last_nep,form='unformatted',status='new')
!         open (877,file=last_clm,form='unformatted',status='new')
c         file last_clm contains data for posible restart of NEM
c         this file is writen at the end of the run
#if ( defined CPL_NEM )
!         open (368,file=init_4nem,form='unformatted',status='old')
c         file init_4nem contains data for the restart of NEM
c         from the results of a previous run
          open (277,file=fnememiss,form='unformatted',status='replace')
          close(277)
#endif
#endif
#if ( defined CPL_OCEANCO2 && defined ML_2D )
          open (369,file=pov_deepo,form='unformatted',status='new')
#endif
        else
C      For restart of the run
          open( unit=546,file=plotfl,
     *    status='OLD',form='unformatted')
          open( unit=547,file=nwrfl,
     *    status='OLD',form='unformatted')
#if ( defined CPL_TEM )
          open (876,file=last_nep,form='unformatted',status='old')
!         open (877,file=last_clm,form='unformatted',status='old')
#if ( defined CPL_NEM )
!         open (368,file=init_4nem,form='unformatted',status='old')
c         file init_4nem contains data for the restart of NEM
c         from the results of a previous run
          open (277,file=fnememiss,form='unformatted',status='replace')
          close(277)
#endif
#endif
#if ( defined CPL_OCEANCO2 && defined ML_2D )
          open (369,file=pov_deepo,form='unformatted',status='old')
#endif
        endif
       endif
C****                                                                   1719.   
C**** RESTART ON DATA SETS 1 OR 2, ISTART=10-13                         1720.   
C****                                                                   1721.   
C**** CHOOSE DATA SET TO RESTART ON                                     1722.   
  400 TAU1=-1.                                                          1723.   
      READ (501,ERR=410) AEXPX1,TAU1                                     1724.   
  410 REWIND 501                                                        1725.   
      TAU2=-1.                                                          1726.   
      READ (502,ERR=420) AEXPX2,TAU2                                     1727.   
  420 REWIND 502                                                        1728.   
	print *,' TAU1=',TAU1,' TAU2=',TAU2
      KDISK=1                                                           1729.   
      IF(TAU1+TAU2.LE.-2.) GO TO 850                                    1730.   
      IF(TAU2.GT.TAU1) KDISK=2                                          1731.   
      if(KDISK.eq.1)AEXPX=AEXPX1
      if(KDISK.eq.2)AEXPX=AEXPX2
      IF(ISTART.GE.13) KDISK=3-KDISK                                    1732.   
      GO TO 450                                                         1733.   
  440 KDISK=ISTART-10                                                   1734.   
C**** RESTART ON UNIT KDISK                                             1735.   
  450 ISTAR0=ISTART                                                     1736.   
      KDISK0=KDISK+500                                                  1737.   
      if(ISTRT1.eq.0) then
C *****
C      FOR ISTRT1 = 0
C *******
c      print *,' Form input '
       if(.not.CONTRR)then
        READ (KDISK0,ERR=840)AEXPX,TAUX,JC1,C1,RC1,KEYNR,U,V,T,P,Q,
     &  ODATA,
     *  GDATA,BLDATA,RQT,SRHR,TRHR,(AJ(K,1),K=1,KACC),TAUY,TSSFC,CKS,   1739.   
     *  CKN,WMGE,TPRIM2,MRCHT,TRSURF,SRSURF,TLANDD,TSURFD,DWAV0
       else
!       READ (KDISK0,ERR=840)AEXPX,TAUX,JC1,C1,RC1,KEYNR,U,V,T,P,Q,
        READ (KDISK0,ERR=840)AEXPX,TAUX,JC,C,RC,KEYNR,U,V,T,P,Q,
     &  ODATA,
     *   GDATA,BLDATA,RQT,SRHR,TRHR,(AJ(K,1),K=1,KACC),TAUY,TSSFC,CKS,  1739.
     *   CKN,WMGE,TPRIM2,MRCHT,TRSURF,SRSURF,TLANDD,TSURFD,DWAV0,
     *   TG3M,RTGO,STG3,DTG3
        endif
       if(.not.CONTRR)then
	JC(16)=JC1(16)
	C(38)=C1(38)
	C(39)=C1(39)
	do i=41,50
c   print *,i,JC1(i)
	  JC(i)=JC1(i)
	end do
	end if
c       print *,' NCNDS=',NCNDS
c       print *,' WITH DEEP MIXED LAYER'
 5001 format(24f5.1)
        print *,' START OF RUN ', AEXP
        print *,' INDAY=',INDAY,' IMONTH=',IMONTH
        print *,' INYEAR=',INYEAR
        print *,'INITIAL CONDITIONs FROM EXP.=',AEXPX
        print *,'JDAY=',JDAY,' JDATE=',JDATE,' JMONTH=',JMONTH
        print *,' JYEAR=',JYEAR
        AEXPX=AEXP
        TAU=0.
        JNDAY=INDAY+JDOFM(IMONTH)
        IYEAR=INYEAR
        IDAY=INDAY+JDOFM(IMONTH)
         TAUI=(IDAY-1)*24.
         TAU=TAUI
         TAUX=TAU
         TAUY=TAU
         JYEAR=INYEAR
         JYEAR0=INYEAR
#if( !defined OCEAN_3D && !defined ML_2D )
         if(TRANSR.and..not.CONTRR)then
         open( unit=575,file=t3file,
     *   status='OLD',form='unformatted')
         read(575)AEXTG3
         read(575)TG3M
         print *,' TG3 from ',AEXTG3
            do 5368 j=1,JM
              STG3(1,j)=0.
              DTG3(1,j)=0.
              do 5368 k=1,lmo
                RTGO(1,j,k)=0.
 5368    continue
        end if
#endif

#if ( defined CPL_OCEANCO2 && defined ML_2D )
       if(ocarcont) then
        if(ocarindata)then
          print *,'Wrong setting of ocarcont and ocarindata'
          print *,ocarcont,ocarindata
          stop
        endif
        open(116,file=deepco2in,
     *   status='old',form='unformatted')
        print *,' AFTER OPEN INIT. data for ocean chem.'
        print *,deepco2in
        read(116)iyearocm,vdfocm
        print *,' iyearocm=',iyearocm
        print *,'Vertical diffusion coefficeint for carbon=',vdfocm
        if(iyearocm.ne.JYEAR-1) then
          print *,'Data for ocean carbon for wrong year'
          print *,' iyearocm=',iyearocm
          print *,' JYEAR=',JYEAR
          stop
        endif
        read(116)Hgin
        read(116)Rco2in
        do k=1,lmo
        do j=1,jm0
         if(k.eq.1)Hg(j)=Hgin(j)
         Rco2(j,k)=Rco2in(j,k)
        end do
        end do
       else
        if(ocarindata)then
         print *,' Reading initial data for ocean carbon'
         open(116,file=deepco2in,
     *    status='old',form='unformatted')
         read(116)
         read(116)Hgin
         read(116)Rco2in
         print *,(Rco2in(j,1),j=1,jm)
         do k=1,lmo
         do j=1,jm0
!         if(k.eq.1)Hg(j)=Hgin(j)
          if(k.eq.1)Hg(j)=1.0e-8
          Rco2(j,k)=Rco2in(j,k)
         end do
         end do
         close(116)
        else
         do k=1,lmo
         do j=1,jm0
          if(k.eq.1)Hg(j)=1.0e-8
          Rco2(j,k)=0.0
         end do
         end do
        endif
       endif
#endif

      else 
C *****
C      FOR ISTRT1 = 1
C *******
        if(TRANSR)then
         READ (KDISK0,ERR=840)AEXPX,TAUX,JC,C,RC,KEYNR,U,V,T,P,Q,ODATA, 1738.
     *   GDATA,BLDATA,RQT,SRHR,TRHR,(AJ(K,1),K=1,KACC),TAUY,TSSFC,CKS,  1739.
     *   CKN,WMGE,TPRIM2,MRCHT,TRSURF,SRSURF,TLANDD,TSURFD,DWAV0,
     *   TG3M,RTGO,STG3,DTG3
        else
         READ (KDISK0,ERR=840)AEXPX,TAUX,JC,C,RC,KEYNR,U,V,T,P,Q,ODATA, 1738.
     *   GDATA,BLDATA,RQT,SRHR,TRHR,(AJ(K,1),K=1,KACC),TAUY,TSSFC,CKS,  1739.
     *   CKN,WMGE,TPRIM2,MRCHT,TRSURF,SRSURF,TLANDD,TSURFD,DWAV0
        endif
C
         if(abs(AEXPX-AEXP).gt.0.05)then
          print *,' DISAGREEMENT BETWEEN AEXPX AND AEXP '
          print *,' FILE ',KDISK0
          print *,' AEXPX=',AEXPX,' AEXP=',AEXP
          stop
         else
          print *,' RESTART OF EXP. ',AEXP
          print *,'JDAY=',JDAY,' JDATE=',JDATE,' JMONTH=',JMONTH
          print *,' JYEAR=',JYEAR
         endif

#if ( defined CPL_OCEANCO2 && defined ML_2D )

        print *,' AFTER OPEN INIT. data for ocean chem.'
        read(369)iyearocm,vdfocm
        print *,' iyearocm=',iyearocm
        print *,'Vertical diffusion coefficeint for carbon=',vdfocm
        if(iyearocm.ne.JYEAR-1) then
c       if(iyearocm.ne.JYEAR) then
          print *,'Data for ocean carbon for wrong year'
          print *,' iyearocm=',iyearocm
          print *,' JYEAR=',JYEAR
          stop
        endif
        read(369)Hgin
        read(369)Rco2in
        print *,(Rco2in(j,1),j=1,jm)
        do k=1,lmo
        do j=1,jm0
         if(k.eq.1)Hg(j)=Hgin(j)
         Rco2(j,k)=Rco2in(j,k)
        end do
        end do
#endif

      endif        ! endif for ISTRT1
#if ( defined CPL_OCEANCO2 && defined ML_2D )
       print *,'Rco2 from input'
!        print *,(Rco2(j,1),j=1,jm)
         print *,Rco2
#endif

      if(KOCEAN.eq.1) then
       print *,' T1 ocean'
       print 5001,(ODATA(1,j,1),j=1,JM)
       print *,' T2 ocean'
       print 5001,(ODATA(1,j,4),j=1,JM)
       print *,' T3 ocean'
       print 5001,(ODATA(1,j,5),j=1,JM)
       endif

c      if(TRANSR)then
c     print *,' STG3'
c     print 5001,(STG3(1,j),j=1,JM)
c     print *,' DTG3/356'
c     print 5001,(DTG3(1,j)/365.,j=1,JM)
c     print *,' RTGO'
c     print 5001,((RTGO(1,j,k),j=1,JM),k=1,lmo)
c     endif
      REWIND KDISK0                                                     1740.   
      ISTART=ISTAR0                                                     1741.   
      KDISK=KDISK0-500                                                  1742.   
      IF(TAUX.NE.TAUY) GO TO 860                                        1743.   
      DO 451 I=1,33                                                     1744.1  
  451 XLABEL(I)=XLABL1(I)                                               1744.2  
      TAU=TAUX                                                          1745.   
      TAUP=TAUX                                                         1746.   
C**** UPDATE C ARRAY FROM INPUTZ                                        1747.   
  500 READ (514,INPUTZ)                                                 1748.   
!      nrad=NDYN
#if ( defined IPCC_EMI )
       CO2=xco2init/GHGBGR(1)
#endif
      INYRAD=INYEAR
      JNDAY=INDAY+JDOFM(IMONTH)
      WRITE (6,907) KDISK,TAUX,AEXP
      ITAUX=TAUX
      NHY=24*365
c     print *,' NHY=',NHY,' ITAUX=',ITAUX,' IYEAR=',IYEAR
      IYEARX=ITAUX/NHY+IYEAR
c     print *,' IYEARX=',IYEARX
      ITAUX=TAUX-(IYEARX-IYEAR)*NHY
c     print *,' ITAUX=',ITAUX
      DO 871 IMNTHX=1,12
       ITM=ITAUX-JDOFM(IMNTHX)*24
       IT=ITAUX-JDOFM(IMNTHX+1)*24
       if(IT.lt.0)go to 872
  871 CONTINUE
  872 CONTINUE
c     print *,' ITAUX=',ITAUX
c     print *,' ITM=',ITM,' IT=',IT
      NWR10Y=20*365*24/5
c     print *,' NWR10Y=',NWR10Y
      IDAYX=ITM/24+1
c     print *,' IYEARX=',IYEARX,' IMNTHX=',IMNTHX,' IDAYX=',IDAYX
      WRCL=((IYEARX-INYEAR)*365.+(JDOFM(IMNTHX)-JDOFM(IMONTH))
     * +IDAYX-INDAY)*24./5.
      NWRCL=(WRCL+0.99)
      INWR=NWRCL/NWR10Y
      NWRCL=NWRCL-NWR10Y*INWR
c     print *,' NWRCL=',NWRCL
      NWRGHG=(IYEARX-INYEAR)*12.+(IMNTHX-IMONTH)
c     print *,' NWRGHG=',NWRGHG
      if(CLDFEED)then
       open( unit=576,file=tsfile,
     * status='OLD',form='unformatted')
       read(576)EXPTSF,TSCNTR
       print *,TSCNTR,' form EXP=',EXPTSF
!      read(576)EXPTSF
!      print *,' TSURF form EXP=',EXPTSF
#if ( defined CLM )
        read(576)TSURFC,TLANDC
#else
        read(576)TSURFC
#endif
!      read(576)TLANDC
      endif
        CORSR=1.0
      if(KOCEAN.eq.0)then
        print *,' fixed SST'
      elseif(KOCEAN.eq.1)then
        print *,' interactive ocean'
#if( !defined OCEAN_3D && !defined ML_2D )
      open( unit=527,file=qffile,
     * status='OLD',form='unformatted')
        READ (527) ANEXPQ,QFLUX,ZOAV
        if(QFCOR)READ (527) CORSR,QFLUX
        do 475 j=1,JM
         QFLUXT(j)=0.
         do 476 n=1,12
           QFLUXT(j)=QFLUXT(j)+QFLUX(j,n)/12.
  476    continue
  475   continue
        print *,' Q-flux from EXP=',ANEXPQ
        print *,(QFLUXT(J),J=1,JM)
#endif
#if(  defined CPL_OCEANCO2 && defined ML_2D )
      open( unit=527,file=qffile,
     * status='OLD',form='unformatted')
        READ (527) ZOAV
      print *,'ZOAV for OCM'
#endif
      else
        print *,' value of KOCEAN is wrong'
        stop
      endif
      NLFR=.NOT.LFR
c     PRINT *,' NLFR=',NLFR
      if(LFR) then
      print *,'  with LAND fractions'
      print *,'  with LAND fractions'
      else
      print *,'  without LAND fractions'
      print *,'  without LAND fractions'
      endif
      if(READGHG.eq.2)then
        open( unit=569,file=ghg_monthly,
     *  status='OLD',form='unformatted')
        open( unit=679,file=ghg_monthly2,
     *  status='OLD',form='unformatted')
	print *,' GHGs from ',ghg_monthly
        if(ISTRT1.eq.1)then
          do 369 i=1,NWRGHG
           do 369 ii=1,13
           read(569)
           if(ii.le.3)read(679)
  369     continue
      endif
      endif
      if(READGHG.eq.1)then
c       open( unit=569,file=dirdat1(1:id1-1)//'ghgsm77',
c    *  status='OLD',form='unformatted')
        open( unit=569,file=ghg_monthly,
     *  status='OLD',form='unformatted')
        open( unit=679,file=ghg_monthly2,
     *  status='OLD',form='unformatted')
        print *,' GHGs from ',ghg_monthly
      endif
      if(PCLOUD.eq.1)then
        print *,' prescribed clouds from GISS GCM'
      elseif(abs(PCLOUD-3.).lt.1.5)then
        open( unit=585,file=clfile,
     *  status='OLD',form='unformatted')
        print *,' fixed  clouds from ',clfile,' for each 5 h.'
        if(ISTRT1.eq.1)then
          do 367 i=1,NWRCL
           read(585)
  367     continue
        endif
        if(PCLOUD.eq.2)then
          print *,' fixed  MC and SS clouds '
        elseif(PCLOUD.eq.3)then
          print *,' fixed SS clouds '
          print *,' interactive MC clouds new scheme'
        elseif(PCLOUD.eq.4)then
          print *,' fixed MC clouds '
          print *,' interactive SS clouds new scheme'
        endif
      elseif(PCLOUD.eq.0)then
        print *,' interactive clouds new scheme'
      elseif(PCLOUD.eq.5)then
        print *,' interactive clouds old scheme'
      elseif(PCLOUD.eq.6)then
      open( unit=528,file=clfile,
     * status='OLD',form='unformatted')
        read (528) EXPCL,cldssm,cldmcm
        print *,' fixed  clouds form EXP=',EXPCL
      else
        print *,' CLOUDS ARE NOT ASSIGNED'
        stop
      endif
      if(WRCLD)then
        if(ISTWRC.eq.0)then
          open( unit=581,file=wrcldf,
     *    status='new',form='unformatted')
        else 
          open( unit=581,file=wrcldf,
     *    status='old',form='unformatted')
        do 368 i=1,NWRCL
          read(581)
  368   continue
        endif
      endif
      TAUE=((LYEAR-IYEAR)*365.+(JDOFM(LMONTH)-JDOFM(IMONTH))+
     *      LDAY-INDAY)*24.+TAUI
      if(ISTART.eq.2)then      
         TAUE=8017.
         if(KOCEAN.eq.1)then
          print *,' SST is not assined'
          stop
         endif
       endif
      REWIND 514                                                        1749.   
      close (514)
      IF (TAU.LT.TAUP-.06125) GO TO 900                                 1750.   
      IF(USET.LE.0.) GO TO 600                                          1751.   
C**** REPOSITION THE OUTPUT TAPE ON UNIT 20 FOR RESTARTING              1752.   
      IF(TAU.LE.TAUO+.06125) GO TO 600                                  1753.   
  520 READ (520,ERR=870,END=880) TAUZ                                   1754.   
      IF(TAU.GE.TAUZ+USET-.06125) GO TO 520                             1755.   
      WRITE (6,908) TAUZ                                                1756.   
C****                                                                   1757.   
C**** CONSTANT ARRAYS TO BE CALCULATED OR READ IN EACH RUN              1758.   
C****                                                                   1759.   
C**** CALCULATE SPHERICAL GEOMETRY                                      1760.   
  600 continue
      TWOPI=8.*atan(1.)
      TWOPI=6.283185                                                    1549.   
      DLON=TWOPI/float(IM)
      DLAT=.5*TWOPI/float(JM-1)
      LAT(1)=-.25*TWOPI                                                 1761.   
      LAT(JM)=-LAT(1)                                                   1762.   
      SINP(1)=-1.                                                       1763.   
      SINP(JM)=1.                                                       1764.   
      COSP(1)=0.                                                        1765.   
      COSP(JM)=0.                                                       1766.   
      DXP(1)=0.                                                         1767.   
      DXP(JM)=0.                                                        1768.   
      DO 620 J=2,JMM1                                                   1769.   
      LAT(J)=LAT(J-1)+DLAT                                              1770.   
      SINP(J)=SIN(LAT(J))                                               1771.   
      COSP(J)=COS(LAT(J))                                               1772.   
  620 DXP(J)=RADIUS*DLON*COSP(J)                                        1773.   
c     print *,(360./TWOPI*acos(COSP(J)),J=1,JM)
c     print *,' COSP'
c     print *,(COSP(J),J=1,JM)
      DO 640 J=2,JM                                                     1774.   
      COSV(J)=.5*(COSP(J-1)+COSP(J))                                    1775.   
      DXV(J)=.5*(DXP(J-1)+DXP(J))                                       1776.   
  640 DYV(J)=RADIUS*(LAT(J)-LAT(J-1))                                   1777.   
c     print *,' DXV(JM/2+1)=',DXV(JM/2+1)
      COSV(JM/2+1)=1.
      DXV(JM/2+1)=RADIUS*DLON
c     print *,' DXV(JM/2+1)=',DXV(JM/2+1)
      print *,' YV'
      print *,(360./TWOPI*acos(COSV(J)),J=2,JM)
      print *,' YP'
      print *,(360./TWOPI*acos(COSP(J)),J=1,JM)
c     print *,' COSV'
c     print *,(COSV(J),J=2,JM)
      DYP(1)=.5*DYV(2)                                                  1778.   
      DYP(JM)=.5*DYV(JM)                                                1779.   
      DXYP(1)=.5*DXV(2)*DYP(1)                                          1780.   
      DXYP(JM)=.5*DXV(JM)*DYP(JM)                                       1781.   
      DXYS(1)=0.                                                        1782.   
      DXYS(JM)=DXYP(JM)                                                 1783.   
      DXYN(1)=DXYP(1)                                                   1784.   
      DXYN(JM)=0.                                                       1785.   
      AREAG=DXYP(1)+DXYP(JM)                                            1786.   
      DO 660 J=2,JMM1                                                   1787.   
      DYP(J)=.5*(DYV(J)+DYV(J+1))                                       1788.   
      DXYP(J)=.5*(DXV(J)+DXV(J+1))*DYP(J)                               1789.   
      DXYS(J)=.5*DXYP(J)                                                1790.   
      DXYN(J)=.5*DXYP(J)                                                1791.   
  660 AREAG=AREAG+DXYP(J)                                               1792.   
      print *,' DXYP'
      print *,(DXYP(J),J=1,JM)
       SS=0.
       SN=0.
       do 578 j=1,12
         SS=SS+DXYP(j)/(RADIUS**2*DLON*DLAT)
         SN=SN+DXYP(j+12)/(RADIUS**2*DLON*DLAT)
  578 continue
      ATMMASS=(ss+sn)*984.*100./9.81
      AREAG=AREAG*FIM                                                   1793.   
      RAVPS(1)=0.                                                       1794.   
      RAVPN(JM)=0.                                                      1795.   
      DO 680 J=2,JM                                                     1796.   
      DXYV(J)=DXYN(J-1)+DXYS(J)                                         1797.   
      RAPVS(J)=.5*DXYS(J)/DXYV(J)                                       1798.   
      RAPVN(J-1)=.5*DXYN(J-1)/DXYV(J)                                   1799.   
      RAVPS(J)=.5*DXYS(J)/DXYP(J)                                       1800.   
  680 RAVPN(J-1)=.5*DXYN(J-1)/DXYP(J-1)                                 1801.   
cprint *,DXP(1),DXP(2),DXV(2)
C**** CALCULATE CORIOLIS PARAMETER                                      1802.   
      OMEGA=TWOPI*(EDPERD+EDPERY)/(EDPERD*EDPERY*SDAY)                  1803.   
      F(1)=-RADIUS*OMEGA*.5*COSP(2)*DXV(2)                              1804.   
      F(JM)=-F(1)                                                       1805.   
      DO 720 J=2,JMM1                                                   1806.   
  720 F(J)=OMEGA*(DXV(J)*DXV(J)-DXV(J+1)*DXV(J+1))/DLON                 1807.   
C**** CALCULATE DSIG AND DSIGO                                          1808.   
      DO 740 L=1,LM                                                     1809.   
  740 DSIG(L)=SIGE(L)-SIGE(L+1)                                         1810.   
      DO 760 L=1,LMM1                                                   1811.   
  760 DSIGO(L)=SIG(L)-SIG(L+1)                                          1812.   

#if ( defined CPL_CHEM )
!
! --- Calculate air mass, First step
! ---   (need to time surface pressure p(i,j) :
!
      i=1
      do 112 k=1,nlev
      do 112 j=1,nlat
        airmass0(i,j,k)=dsig(k)*dxyp(j)*100.
     &                /grav
112   continue

!	open(122,file='airmass0.dat',form='unformatted',
!     &		 status='unknown')
!	write(122)airmass0
!	stop
!
#endif

C**** READ IN FDATA: PHIS, PLAND AND RLICE                              1813.   
      READ (526) FDATA                                                  1814.   
      REWIND 526                                                        1815.   
      print *,' NLFR=',NLFR,' IO=',IO
      DO 283 J=1,JM                                                     1815.5  
      DO 283 I=1,IO                                                     1815.51 
        FDATA(I,J,1)=0.
#if ( defined ML_2D)
      if(FDATA(I,J,2).ge.0.94)then
         FDATA(I,J,2)=1.00
      endif
#endif
      C3LICE(I,J)=FDATA(I,J,2)*FDATA(I,J,3)                             1815.52 
      C3LAND(I,J)=FDATA(I,J,2)                                          1815.53 
      if(NLFR)FDATA(I,J,2)=0.                                           1815.54
 283  continue
      do 284 J=1,JM
        ILAND=0.
        IICE=0.
        CONT1=0.
        CONT2=0.
        do 285 I=1,IO
           PLAND=FDATA(I,J,2)
           PICE=FDATA(I,J,3)
           CONT1=CONT1+PLAND 
           CONT2=CONT2+PICE
           ILAND=ILAND+1
           IF(PLAND.GT.0.)IICE=IICE+1
  285   continue
        do 286 I=1,IO
           IF(ILAND.GT.0)FDATA(I,J,2)=CONT1/ILAND
c          IF(FDATA(I,J,2).LT.0.01)FDATA(I,J,2)=0.
           IF(IICE.GT.0)FDATA(I,J,3)=CONT2/IICE
  286   continue
        fland_temp(j)=FDATA(1,J,2)
  284 continue
#if( defined OCEAN_3D)
Cjrs        if(jmocean.ne.jm0-2)then
C          print *,"Wrong jm or jmocean"
C          stop
C       endif
       OCNGEOM=.false.
        print *,'With land/ocean fractions directly from 3D ocean model'
        print *,cflan
        do i=1,IO
CJRS         FDATA(I,1,2)=cflan(1)
C         do j=2,jm0-1
C            FDATA(I,J,2)=cflan(j-1)
C         enddo
         do j=1,jm0
            FDATA(I,J,2)=cflan(j)
         enddo
CJRS         FDATA(I,JM0,2)=cflan(jmocean)
       enddo
      do j=1,jm0
       fland_atm(j)=FDATA(1,J,2)
       enddo
       if(ISTRT1.eq.0)then
        open (505,file=ocndata4atm,form='unformatted',
     &    status='new') 
       else
        open (505,file=ocndata4atm,form='unformatted',
     &    status='old') 
       endif
#else
       if(OCNGEOM)then
        print *,'With land/ocean fractions as in 3D ocean model'
        open (626,file=ocngmfile,
     & status='old')
        do j=1,jm0
           read (626,*),iii,fo3d
           print *,360./TWOPI*acos(COSP(J)),fo3d
           do i=1,IO
            FDATA(I,J,2)=1.0-fo3d
           enddo
        enddo
      endif
#endif
#if ( defined CLM )
       open (767,file=fclmlice,status='old')
       do j=1,jm
        read(767,*),clmlice(j)
        FDATA(1,J,3)=0.01*clmlice(j)
       enddo
       close (767)
       open (767,file=fbaresoil,status='old')
       do j=1,jm
        read(767,*),baresoil(j),baresoil(j)
       enddo
       close (767)
       open (767,file=fwmax,status='old')
       read (767,*),lineclm
       read(767,*),(w1maxclm(j),j=1,jm)
       read (767,*),lineclm
       read(767,*),(w2maxclm(j),j=1,jm)
       read (767,*),lineclm
       read(767,*),(vmaskclm(j),j=jm,1,-1)
       close (767)
!      open (767,file=fprratio,status='old')
!      do j=1,jm
!       read(767,*),(prlnd2total(j,n),n=1,12)
!      enddo
!      close (767)
#else
!      do j=1,jm
!      do n=1,12
!       prlnd2total(j,n)=1.0
!      enddo
!      enddo
#endif
       open (767,file=fprratio,status='old')
       do j=1,jm
        read(767,*),(prlnd2total(j,n),n=1,12)
       enddo
       close (767)
       print *,'Ratio of land precipitation to total'
       do j=1,jm
        print('12f7.4'),(prlnd2total(j,n),n=1,12)
       enddo
      print *,' FDATA(1,J,2)='
      print *,(FDATA(1,J,2),J=1,JM)
      print *,' FDATA(1,J,3)='
      print *,(FDATA(1,J,3),J=1,JM)

#if( !defined OCEAN_3D)
C**** READ IN MAXIMUM MIXED LAYER DEPTHS FOR PREDICTED OCEAN RUNS       1815.6  
      IF(KOCEAN.NE.1) GO TO 764                                         1815.61 
      READ (525) Z12O                                                   1815.62 
      REWIND 525                                                        1815.63 
      DO 628 J=1,JM                                                     1815.641
      SUM2=0.                                                           1815.643
      CONT1=0.                                                          1815.644
      DO 626 I=1,IO                                                     1815.645
      PWATER=1.-C3LAND(I,J)                                             1815.647
      IF(PWATER.LE.0.) GO TO 626                                        1815.648
      CONT1=CONT1+PWATER                                                1815.649
      SUM2=SUM2+Z12O(I,J)*PWATER                                        1815.651
  626 CONTINUE                                                          1815.652
      IF(CONT1.LE.0.) GO TO 628                                         1815.653
      IF(J.EQ.1.OR.J.EQ.JM) GO TO 628                                   1815.654
      SUM2=SUM2/CONT1                                                   1815.656
      DO 627 I=1,IO                                                     1815.657
  627 Z12O(I,J)=SUM2                                                    1815.659
  628 CONTINUE                                                          1815.66 
      DO 629 J=1,2                                                      1815.661
      DO 629 I=1,IO                                                     1815.662
  629 Z12O(I,J)=Z12O(I,3)                                               1815.664
  764 CONTINUE
#endif

C**** READ IN EARTH RATIOS FOR THE 8 VEGETATION TYPES AND THE VADATA :  1816.   
C     VADATA(TYPE,SEASON,1)=GROUND ALBEDO FOR A GIVEN TYPE AND SEASON   1817.   
C         1       2       3       4       5       6       7       8     1818.   
C       DESRT   TNDRA   GRASS   SHRUB   TREES   DECID   EVERG   RAINF   1819.   
C SPRN   0.35,   0.12,   0.16,   0.16,   0.14,   0.18,   0.12,   0.11,  1820.   
C SUMR   0.35,   0.12,   0.20,   0.18,   0.14,   0.12,   0.12,   0.11,  1821.   
C FALL   0.35,   0.17,   0.20,   0.25,   0.17,   0.15,   0.15,   0.11,  1822.   
C WNTR   0.35,   0.15,   0.18,   0.20,   0.12,   0.12,   0.11,   0.11/  1823.   
C                                                                       1824.   
C     VADATA(TYPE,SEASON,2)=RATIO OF NEAR IR ALBEDO TO VIS ALBEDO FOR...1825.   
C         1       2       3       4       5       6       7       8     1826.   
C       DESRT   TNDRA   GRASS   SHRUB   TREES   DECID   EVERG   RAINF   1827.   
C SPRN    1.0,    3.0,    3.0,    3.0,    3.0,    3.0,    3.0,    3.0,  1828.   
C SUMR    1.0,    3.3,    3.5,    3.0,    3.3,    4.0,    3.0,    3.0,  1829.   
C FALL    1.0,    3.5,    4.0,    3.0,    3.5,    5.0,    3.0,    3.0,  1830.   
C WNTR    1.0,    3.2,    3.5,    3.0,    3.2,    4.0,    3.0,    3.0/  1831.   
C                                                                       1832.   
C     VADATA(TYPE,1,3)=MASKING DEPTH FOR A GIVEN TYPE                   1833.   
C                                                                       1834.   
C         1       2       3       4       5       6       7       8     1835.   
C       DESRT   TNDRA   GRASS   SHRUB   TREES   DECID   EVERG   RAINF   1836.   
C        10.,    20.,    20.,    50.,   200.,   500.,   1000.,  2500.,  1837.   
C                                                                       1838.   
C     VADATA(TYPE,1+K,3)=WATER FIELD CAPACITY FOR K-TH GROUND LAYER     1839.   
C                                                                       1840.   
C    1   10.,    30.,    30.,    30.,    30.,    30.,     30.,   200.,  1841.   
C    2   10.,   200.,   200.,   300.,   300.,   450.,    450.,   450.,  1842.   
C   (3)   0.,     0.,     0.,     0.,     0.,     0.,      0.,    0./   1843.   
C                                                                       1844.   
#if ( !defined CLM )
      if(VEGCH.or.TRVEG)then
         READ(523) IYVEG
c        if(.not.TRVEG)then
          print *,'VDATA for year=',iyveg
c        endif
      IYVEGIN=IYVEG
      endif
      READ (523) (((VDATA(I,J,K),I=1,IO),J=1,JM),K=1,8),                1845.   
     *  (((VADATA(I,J,K),I=1,8),J=1,4),K=1,3)                           1845.1  
        print *,'VADATA'
        do k=1,3
        print *,' K=',k
        print '(8f7.2)',((VADATA(I,J,K),I=1,8),J=1,4)
        enddo
       if(ISTRT1.eq.1.and.TRVEG)then
        print *,'Restart with TRVEG'
        if(JYEAR.le.1992)then
          JYEARV=JYEAR
        else
          JYEARV=1992
          print *,' End of vegfile has been reached'
          print *,' VDATA for year 1992 are used for the rest of run'
        endif
        do ii=IYVEGIN,JYEARV-1
          READ(523) IYVEG
          READ (523) (((VDATA(I,J,K),I=1,IO),J=1,JM),K=1,8)
        enddo
        if(IYVEG.ne.JYEARV)then
         print *,' Wrong IYVEG'
         print *,' IYVEG=',IYVEG,' JYEARV=',JYEARV
         stop
        endif
        print *,'VDATA for year=',JYEARV
      endif
c     REWIND 523                                                        1846.   
C**** MODIFY THE VADATA IF DESIRED                                      1847.   
C     NO MODIFICATIONS                                                  1848.   
C**** COMPUTE WATER FIELD CAPACITIES FOR GROUND LAYERS 1 AND 2          1849.   
      IOFF=0                                                            1849.1  
      IF(VADATA(4,2,3).LT.100.) IOFF=1                                  1849.2  
      ERROR=.001                                                        1849.3  
      DEFLT=24.                                                         1850.   
      DO 785 L=1,2                                                      1851.   
      DO 780 J=1,JM                                                     1852.   
      DO 780 I=1,IO                                                     1853.   
      WFCIJL=0.                                                         1854.   
      DO 770 K=1,8                                                      1855.   
  770 WFCIJL=WFCIJL+VDATA(I,J,K)*VADATA(K,L+IOFF,3)                     1856.   
      IF (WFCIJL.LT.1.) WFCIJL=DEFLT                                    1857.   
      IF(ISTART.NE.2) GO TO 780
      IF(GDATA(I,J,4*L+1)+GDATA(I,J,4*L+2).LE.WFCIJL) GO TO 780         1858.   
      X=WFCIJL/(GDATA(I,J,4*L+1)+GDATA(I,J,4*L+2)+1.E-3)                1859.   
      GDATA(I,J,4*L+1)=GDATA(I,J,4*L+1)*X                               1860.   
      GDATA(I,J,4*L+2)=GDATA(I,J,4*L+2)*X                               1861.   
  780 VDATA(I,J,L+8)=WFCIJL                                             1862.   
      DEFLT=60.                                                         1863.   
  785 CONTINUE                                                          1864.   
      DO 765 K=1,10                                                     1864.5  
      DO 765 J=2,JMM1                                                   1864.51 
      CONT1=0.                                                          1864.52 
      SUM1=0.                                                           1864.53 
      DO 766 I=1,IO                                                     1864.54 
      PEARTH=C3LAND(I,J)-C3LICE(I,J)                                    1864.55 
      CONT1=CONT1+PEARTH                                                1864.56 
 766  SUM1=SUM1+PEARTH*VDATA(I,J,K)                                     1864.57 
      IF (CONT1.LE.0.) GO TO 765                                        1864.58 
      SUM1=SUM1/CONT1                                                   1864.59 
      DO 767 I=1,IO                                                     1864.6  
 767  VDATA(I,J,K)=SUM1                                                 1864.61 
 765  CONTINUE                                                          1864.62 
c     print *,' BEAR LAND'
c     print '(12f7.2,/,11f7.2)',(VDATA(1,j,1),j=1,JM)
      print *,' INPUT'
      print *,' WMAX1'
      print '(12f7.2,/,11f7.2)',(VDATA(1,j,9),j=1,JM)
      print *,' WMAX2'
      print '(12f7.2,/,11f7.2)',(VDATA(1,j,10),j=1,JM)
C *************
      print *,'Vadata'
        print '(8f7.2)',(VADATA(K,4,3),K=1,8)
      DO  K=1,8
c     VADATA(K,4,3)=0.1*VADATA(K,4,3)
      VADATA(K,4,3)=VADATA(K,3,3)
      ENDDO
        print '(8f7.2)',(VADATA(K,4,3),K=1,8)
      DO J=1,JM 
      CONT1=0.   
      SUM1=0.    
      DO  I=1,IO
      WFCIJL=0.   
      PEARTH=C3LAND(I,J)-C3LICE(I,J)
      CONT1=CONT1+PEARTH
c     SUM1=SUM1+PEARTH*WFCIJL
      DO  K=1,8
      WFCIJL=WFCIJL+VDATA(I,J,K)*VADATA(K,4,3)
      ENDDO     ! K
      SUM1=SUM1+PEARTH*WFCIJL
      ENDDO     ! I
      IF (CONT1.LE.0.) GO TO 865            
      SUM1=SUM1/CONT1
      VMASK(J)=SUM1
 865  CONTINUE
      ENDDO     ! J
      print *,' VMASK form NP to SP in meters of water'
      print '(12f7.2,/11f7.2)',(VMASK(jm-j+1),j=1,JM)
C ************
#else
      READ (523) (((VDATA(I,J,K),I=1,IO),J=1,JM),K=1,8),                1845.   
     *  (((VADATA(I,J,K),I=1,8),J=1,4),K=1,3)                           1845.1  
        print *,'VADATA'
        do k=1,3
        print *,' K=',k
        print '(8f7.2)',((VADATA(I,J,K),I=1,8),J=1,4)
        enddo
      do j=1,jm
        VDATA(1,j,1)=0.01*baresoil(j)
        VDATA(1,j,2)=1.-0.01*baresoil(j)
        do k=3,8
         VDATA(1,j,k)=0.0
        enddo
        VDATA(1,j,9)=w1maxclm(j)
        VDATA(1,j,10)=w2maxclm(j)
      enddo
      print *,'Vadata'
        print '(8f7.2)',(VADATA(K,4,3),K=1,8)
      DO  K=1,8
c     VADATA(K,4,3)=0.1*VADATA(K,4,3)
      VADATA(K,4,3)=VADATA(K,3,3)
      ENDDO
        print '(8f7.2)',(VADATA(K,4,3),K=1,8)
#endif

      CALL RINIT (IRAND)                                                1865.   
C     CALL IJSET (IM,JM,FDATA(1,1,2))                                   1866.   
      WRITE (6,INPUTZ)                                                  1867.   
C     information for main program
      AEXPA=AEXP
      INDAYA=INDAY
      IMONTHA=IMONTH
      INYEARA=INYEAR
      LDAYA=LDAY
      LMONTHA=LMONTH
      LYEARA=LYEAR
      IYEARA=JYEAR
      IRESTART=ISTRT1
      RETURN                                                            1868.   
C****                                                                   1869.   
C**** TERMINATE BECAUSE OF IMPROPER PICK-UP                             1870.   
C****                                                                   1871.   
  800 WRITE (6,910) ISTART                                              1872.   
      STOP 3                                                            1873.   
  810 WRITE (6,911) TAUP,TAUX                                           1874.   
      STOP 3                                                            1875.   
  820 WRITE (6,912) TAUP,TAUX                                           1876.   
      STOP 3                                                            1877.   
  830 WRITE (6,913)                                                     1878.   
      STOP 3                                                            1879.   
  840 IF(3-KDISK.EQ.KLAST) GO TO 850                                    1880.   
      REWIND KDISK                                                      1881.   
      KLAST=KDISK                                                       1882.   
      KDISK=3-KDISK                                                     1883.   
      WRITE (6,914) KLAST,KDISK                                         1884.   
      GO TO 450                                                         1885.   
  850 WRITE (6,915)                                                     1886.   
      STOP 3                                                            1887.   
  860 WRITE (6,916) TAUX,TAUY                                           1888.   
      STOP 3                                                            1889.   
  870 WRITE (6,917) TAUZ,TAU                                            1890.   
      STOP 3                                                            1891.   
  880 WRITE (6,918) TAUZ,TAU                                            1892.   
      STOP 3                                                            1893.   
  890 WRITE (6,919) ISTART                                              1894.   
      STOP 3                                                            1895.   
  900 WRITE (6,920) TAUP,TAU                                            1896.   
      STOP 3                                                            1897.   
C****                                                                   1898.   
  901 FORMAT ('0',40X,'GISS N LAYER WEATHER MODEL'/)                    1899.   
  902 FORMAT (20A4/11A4,A2,30X,A4)                                      1900.   
  903 FORMAT ('0',31A4,A3,A4/)                                          1901.   
  904 FORMAT (10A12)                                                     1902.   
  905 FORMAT (35X,10A12)                                                 1903.   
  906 FORMAT ('0ATMOSPHERIC I.C.  ISTART,TAUX=',I4,F10.2,3X,20A4)       1904.   
  907 FORMAT ('0RESTART DISK READ ON UNIT',I2,', TAUX=',F9.2,'AEXP=',
     *           F9.2,/,3X,20A4,A6)                                     1905.   
  908 FORMAT ('0OUTPUT TAPE REPOSITIONED.  LAST TAU READ WAS',F9.2)     1906.   
  910 FORMAT ('0ERROR ENCOUNTERED READING I.C. ON UNIT 9.  ISTART=',I4) 1907.   
  911 FORMAT ('0EOF ON UNIT 9.  LATER I.C. NEEDED.  TAUP,TAUX=',2F10.2) 1908.   
  912 FORMAT ('0EARLIER I.C. NEEDED ON UNIT 9.  TAUP,TAUX=',2F10.2)     1909.   
  913 FORMAT ('0ERROR ENCOUNTERED READING GROUND CONDITIONS ON UNIT 7.')1910.   
  914 FORMAT ('0ERROR ENCOUNTERED READING RESTART TAPE ON UNIT',I3/,    1911.   
     *  '  TRY TO RESTART THE JOB WITH ISTART=3,KDISK=',I1)             1912.   
  915 FORMAT ('0ERRORS ON BOTH RESTART DATA SETS.')                     1913.   
  916 FORMAT ('0TAUX,TAUY=',2F10.2/'0DISK RESTART FILE DESTROYED, TRY T 1914.   
     * RESTART THE JOB WITH ISTART=99, OR TERMINATE THE JOB.')          1915.   
  917 FORMAT ('0ERROR ENCOUNTERED REPOSITIONING TAPE ON UNIT 27.  TAUZ,T1916.   
     *AU=',2F10.2)                                                      1917.   
  918 FORMAT ('0EOF ON UNIT 20 WHILE REPOSITIONING TAPE.  TAUZ,TAU=',   1918.   
     *  2F10.2)                                                         1919.   
  919 FORMAT ('0INCORRECT VALUE OF ISTART',I5)                          1920.   
  920 FORMAT ('0PREVIOUS TAUE=',F10.2,' WAS NOT YET REACHED. TAU=',     1921.   
     *  F10.2,' RESUBMIT THE JOB WITH AN EARLIER TAUE CARD')            1922.   
      END                                                               1923.   
