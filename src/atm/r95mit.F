
#include "ctrparam.h"

!	==========================================================
!
!	R95MIT.F:  Model II radiation: 1958 Atmosphere and mean 
!			strat aeros (.012) Zenith angle dependence 
!			for aerosols not used.
!
!	----------------------------------------------------------
!
!	Author of Chemistry Modules:	Chien Wang
!
!	----------------------------------------------------------
!
!	Important Note: Because the original components of chemistry 
!		module in this file are used by some runs not using
!		interactive chemistry-climate model, therefore, the 
!		cpp header CPL_CHEM is barely applied. Instead, 
!		PREDICTED_GASES is appearing at related places.
!
!		Chien Wang
!		080100
!		 
!	Revision History:
!	
!	When	Who		What
!	----	----------	-------	
!	073100	Chien Wang	repack based on CliChem3 & M24x11,
!				  add cpp, and float -> dble
!	093001  Chien Wang      add bc & oc and rewrote aerosol/radiation
!				  interface including S(VI)
!	062604	Chien Wang	merge with current igsm module
!
!	==========================================================


C**** R83XX B83XX R83ZA                      02/4/93                       0.1  
C**** OPT(3)                                                               0.2  
C****                                                                      0.3  
C**** Model II radiation: 1958 Atmosphere and mean strat aeros (.012)      0.4  
C**** Zenith angle dependence for aerosols not used.                       0.5  
***** R83ZA B83XX R83ZA                      12/23/91                      0.1  
***** OPT(3)                                                               0.2  
*****                                                                      0.3  
***** Model II radiation with 1958 Atmosphere, mean strat aeros (.012).    0.4  
***** Aerosols: Zenith angle dependence and other changes implemented      0.5  
C     SUBROUTINE RCOMP1(NFTTTR,NFTTSR,NFTFOR)                              1.   
c     6/20/2005
      SUBROUTINE RCOMP1(NFTTTR,NFTTSR,NFTFOR,KTREND) 

#include "chem_para"
#include "chem_com"
#include "B83XX.COM"						          1012.

c     DOUBLE PRECISION PFOFTK,TKOFPF,WAVNA,WAVNB,TK,PFWI                  64.   
      DATA WAVNA/850.0/,WAVNB/900.0/                                      64.5  
C                                                                         65.   
      NKSR=6                                                              66.   
C               -----------------------------------------------------     67.   
C               READ IN GAS TAU TABLE AND DISTRIBUTED PLANCK FUNCTION     68.   
C               ALSO THERMAL RAD AEROSOL, CLOUD & SURFACE ALBEDO DATA     69.   
C                                                                         70.   
C               IF(NFTTTR.GE.1) TAU TABLE DATA ARE READ (UNIT=NFTTTR)     71.   
C                               WINDOW FLUX B TEMP CONVERSION FACTORS     72.   
C                               ARE ALSO COMPUTED AT THIS TIME            73.   
C                                                                         74.   
C               IF(NFTTTR.LT.1) TAU TABLE DATA ARE NOT READ FROM DISK     75.   
C                               WINDOW FLUX B TEMP CONVERSION FACTORS     76.   
C                               ARE NOT COMPUTED                          77.   
C                               COMMON/RADCOM/PARAMETERS CAN BE RESET     78.   
C                               MORE CONVENIENTLY                         79.   
C               -----------------------------------------------------     80.   
C                                                                         81.   
      IF(NFTTTR.LT.1) GO TO 110                                           82.   
C                                                                         83.   
C$    REWIND NFTTTR                                                       84.   
C$    READ  (NFTTTR) ITRHDR,TAUTBL,PLANCK,TRAQEX,TRAQSC,TRACOS            85.   
C$   +              ,TRCQEX,TRCQSC,TRCCOS,AOCEAN,AGSIDV,CLDALB            86.   
C$   +              ,TRACEG                                               87.   
      REWIND NFTTTR                                                       88.   
      READ  (NFTTTR) TAUTBL                                               89.   
      REWIND NFTTTR                                                       89.5  
C                                                                         90.   
      NFTTTP=NFTTTR+1                                                     91.   
      REWIND NFTTTP                                                       92.   
      READ  (NFTTTP) PLANCK                                               93.   
      REWIND NFTTTP                                                       93.5  
C                                                                         94.   
C                                                                         95.   
C$    IF(NFTTSR.GT.1) REWIND NFTTSR                                       96.   
C$    IF(NFTTSR.GT.1) READ (NFTTSR) ISRHDR,SRTBL                          97.   
C                                                                         98.   
      ID5(1)=8304                                                         99.   
      ID5(2)=8106                                                        100.   
      ID5(3)=8106                                                        101.   
C                                                                        102.   
      NKTR  =25                                                          103.   
      IT0   =123                                                         104.   
      ITNEXT=250                                                         105.   
C                                                                        106.   
C     ---------------------------------------------------------------    107.   
C     DEFINE WINDOW FLUX TO BRIGHTNESS TEMPERATURE CONVERSION FACTORS    108.   
C     ---------------------------------------------------------------    109.   
C                                                                        110.   
      DO 100 I=1,630                                                     111.   
      PFWI=0.001*I                                                       112.   
      IF(I.GT.100) PFWI=(0.1+0.01*(I-100))                               113.   
      IF(I.GT.190) PFWI=(1.0+0.10*(I-190))                               114.   
  100 TKPFW(I)=TKOFPF(WAVNA,WAVNB,PFWI)                                  115.   
  110 CONTINUE                                                           116.   
C                 ---------------------------------------------------    117.   
C                 SET ALBEDO,GAS,AEROSOL DISTRIBUTIONS & COEFFICIENTS    118.   
C                 ALSO CALLED ARE ALBDAY,O3DDAY,O3DLAT FOR  JDAY,JLAT    119.   
C--------------------------------------------------------------------    120.   
C                                                                        121.   
      IF(KFORCE.GT.0) CALL SETFOR(NFTFOR)                                122.   
      IF(LASTVC.GE.0) CALL SETATM                                        123.   
C **************** Print *******
c       print *,' from RCOMP1'
c       print *,'JMLAT=',JMLAT
c       print *,'DLAT'
c       print *,DLAT
C **************** Print *******
      CALL  SETALB                                                       124.   
c     CALL  SETGAS                                                       125.   
c  6/202005
!     print *,' Before CALL  SETGAS'
      CALL  SETGAS(KTREND)
!     print *,' After CALL  SETGAS'
      CALL  SETAER                                                       126.   

C-----------------                                                       127.   
      RETURN                                                             128.   
C                                                                        129.   
C----------------------------------------------------------------------- 130.   
C                     RESET SEASON (JDAY) DEPENDENT QUANTITIES AS NEEDED 131.   
      ENTRY RCOMPT                                                       132.   
c     print *,' from RCOMPT JDAY=',JDAY
c     print *,'PPMV58 from RCOMPT'
c      print *,PPMV58
C----------------------------------------------------------------------- 133.   
C                                                                        134.   
      IF(KFORCE.GT.0) CALL GETFOR                                        135.   
      IF(LAPGAS.EQ.2) CALL SETLAP                                        136.   
      CALL  ALBDAY                                                       137.   
      CALL  O3DDAY                                                       138.   

#if ( !defined CPL_CHEM ) && ( defined O3_FORCED )
C--Addition by CEForest (18Mar98)
!call getforcedozone(jyear,jday)
	call getforcedozone(jyearo3,jday)
C--End of Addition
#endif

      RETURN                                                             139.   
C                                                                        140.   
C----------------------------------------------------------------------- 141.   
C                   RESET LATITUDE (JLAT) DEPENDENT QUANTITIES AS NEEDED 142.   
      ENTRY RCOMPJ                                                       143.   
C----------------------------------------------------------------------- 144.   
      CALL  O3DLAT                                                       145.   
      RETURN                                                             146.   
C                                                                        147.   
C----------------------------------------------------------------------- 148.   
C                  GET ALBEDO,GAS AEROSOL DATA THEN COMPUTE THERML/SOLAR 149.   
      ENTRY RCOMPX                                                       150.   
C----------------------------------------------------------------------- 151.   
      CALL  GETALB                                                       152.   
      CALL  GETGAS                                                       153.   
      CALL  GETAER                                                       154.   
C                        --------------------------------------------    155.   
C                        SPECIFY SURFACE LAYER GAS ABSORPTION AMOUNTS    156.   
C                        --------------------------------------------    157.   
      DO 350 K=1,11                                                      158.   
      TAUSL(K)=RATQSL*FRACSL*TAUN(1+K*NL-NL)                             159.   
  350 TAUN(1+K*NL-NL)=TAUN(1+K*NL-NL)-TAUSL(K)                           160.   
      DO 360 K=12,NKTR                                                   161.   
      TAUSL(K)=       FRACSL*TAUN(1+K*NL-NL)                             162.   
  360 TAUN(1+K*NL-NL)=TAUN(1+K*NL-NL)-TAUSL(K)                           163.   
C-----------------                                                       164.   
      CALL  THERML                                                       165.   
C-----------------                                                       166.   
      IF(KGASSR.GT.0) CALL SOLGAS                                        167.   
C                                                                        168.   
C$ ****************************COMMENTED OUT CARDS INTERPOLATE SOLAR TAU 169.   
C$    DO 300 I=1,600                                                     170.   
C$300 SRTAU(I)=0.                                                        171.   
C$    CALL SRGAS (NL,PL,DPL,TLM,ULGAS,SRTAU,SRTBL,1,3)                   172.   
C----------------                                                        173.   
      CALL  SOLAR                                                        174.   
C----------------                                                        175.   
C                                                                        176.   

      RETURN                                                             177.   
      END                                                                178.   
      SUBROUTINE SETALB                                                  179.   
      COMMON/ADDALB/BVSURFA,XVSURFA,BNSURFA,XNSURFA

#include "B83XX.COM"                                                     180.   
#include "chem_para"
#include "chem_com"
#if ( defined OCEAN_3D )
      PARAMETER (IM0=01,JM0=N_LAT,LM0=N_LEV,IO0=N_LON0,KAIJ0=75)
#include "AGRID.h"
#endif
#if ( defined CLM )
#include "CLM.h"
#endif

      DIMENSION ALVISK(11,4),ALNIRK(11,4),FIELDC(11,3),VTMASK(11)        241.   
C$$   DIMENSION ALMEAN(11,4),RATIRV(11,4)                                242.   
C                                                                        243.   
      EQUIVALENCE                                                        244.   
     +          (VADATA(1,1,1),ALVISK(1,1)),(VADATA(1,1,2),ALNIRK(1,1))  245.   
     +,         (VADATA(1,1,3),FIELDC(1,1)),(VADATA(1,4,3),VTMASK(1))    246.   
C$$  +          (VADATA(1,1,1),ALMEAN(1,1)),(VADATA(1,1,2),RATIRV(1,1))  247.   
C$$  +,         (VADATA(1,1,3),VTMASK(1)),(VADATA(1,2,3),FIELDC(1,1))    248.   
C                                                                        249.   
      EQUIVALENCE                                                        250.   
     +          (FEMTRA(1),ECLTRA),  (FZASRA(1),ZCLSRA)                  251.   
     +,         (FEMTRA(2),EOCTRA),  (FZASRA(2),ZOCSRA)                  252.   
     +,         (FEMTRA(3),ESNTRA),  (FZASRA(3),ZSNSRA)                  253.   
     +,         (FEMTRA(4),EICTRA),  (FZASRA(4),ZICSRA)                  254.   
     +,         (FEMTRA(5),EDSTRA),  (FZASRA(5),ZDSSRA)                  255.   
     +,         (FEMTRA(6),EVGTRA),  (FZASRA(6),ZVGSRA)                  256.   
C                                                                        257.   
      EQUIVALENCE                                                        258.   
     + (BXA(1),BOCVIS),(BXA(5),BEAVIS),(BXA( 9),BOIVIS),(BXA(13),BLIVIS) 259.   
     +,(BXA(2),BOCNIR),(BXA(6),BEANIR),(BXA(10),BOINIR),(BXA(14),BLINIR) 260.   
     +,(BXA(3),XOCVIS),(BXA(7),XEAVIS),(BXA(11),XOIVIS),(BXA(15),XLIVIS) 261.   
     +,(BXA(4),XOCNIR),(BXA(8),XEANIR),(BXA(12),XOINIR),(BXA(16),XLINIR) 262.   
C                                                                        263.   
     +,               (BXA(17),EXPSNE),(BXA(18),EXPSNO),(BXA(19),EXPSNL) 264.   
C                                                                        265.   
     +,               (BXA(20),BSNVIS),(BXA(21),BSNNIR)                  266.   
     +,               (BXA(22),XSNVIS),(BXA(23),XSNNIR)                  267.   
C                                                                        268.   
      EQUIVALENCE (ASNALB(1),ASNVIS),(ASNALB(2),ASNNIR)                  269.   
      EQUIVALENCE (AOIALB(1),AOIVIS),(AOIALB(2),AOINIR)                  270.   
      EQUIVALENCE (ALIALB(1),ALIVIS),(ALIALB(2),ALINIR)                  271.   
      DIMENSION SRBALB(6),SRXALB(6)                                      272.   
C                                                                        273.   
      EQUIVALENCE (SRBXAL(1,1),SRBALB(1)),(SRBXAL(1,2),SRXALB(1))        274.   
C                                                                        275.   
C                          1       2       3       4                     276.   
C                        WINTER  SPRING  SUMMER  AUTUMN                  277.   
      REAL SEASON(4)/  15.00,  105.0,  196.0,  288.0/                    278.   
C                                                                        279.   
C----------------------------------------------------------------------- 280.   
C     SOLAR: OCEAN ALBEDO DEPENDENCE ON ZENITH ANGLE & WIND SPEED        281.   
C                                                                        282.   
      BVH2O(WMAG)=.0488+.0974/(5.679+WMAG)+.0004/(.3333+WMAG)            283.   
      XVH2O(WMAG,X)=.021+X*X*(.0421+X*(.1283+X*(-.04+X*(3.117/           284.   
     +              (5.679+WMAG)+X*.025/(.3333+WMAG)))))                 285.   
C----------------------------------------------------------------------- 286.   
C                                                                        287.   
      data IFFF0 /1/
      data ICLM /1/
C **************** Print *******
cprint *,' from SETALB'
cprint *,'DLAT'
cprint *,DLAT
cprint *,'VTMASK'
cprint *,VTMASK
cprint *,'FIELDC'
cprint *,FIELDC
C **************** Print *******
      JNORTH=JMLAT/2+1                                                   288.   
      VISNIR=MEANAL                                                      289.   
C                                                                        290.   
C**** FOR OLD ALBEDO FILES COMPUTE VISUAL AND NEAR-IR ALBEDOS            290.8  
c     print *,' SETALB NV=',NV
c     print *,' ALVISK=',ALVISK
c     print *,' ALNIRK=',ALNIRK
      IF (VADATA(4,2,3).GT.100.) GO TO 101                               290.9  
      DO 100 L=1,4                                                       291.   
      DO 100 K=1,8                                                       292.   
      ALMEAN=ALVISK(K,L)                                                 292.1  
      RATIRV=ALNIRK(K,L)                                                 292.2  
      ALVISK(K,L)=ALMEAN/(0.6+0.4*RATIRV)                                293.   
  100 ALNIRK(K,L)=ALMEAN/(0.4+0.6/RATIRV)                                294.   
  101 CONTINUE                                                           294.1  
C                                                                        295.   
C----------------------------------------------------------------------- 296.   
C                   DEFINE SEASONAL ALBEDO (ALVISD,ALNIRD) FOR VEG TYPES 297.   
      ENTRY ALBDAY                                                       298.   
C----------------------------------------------------------------------- 299.   
C                                                                        300.   
      XJDAY=JDAY                                                         301.   
c     print *,'from ALBDAY XJDAY=',XJDAY
      SEASN1=-77.0                                                       302.   
      DO 110 K=1,4                                                       303.   
      SEASN2=SEASON(K)                                                   304.   
      IF(XJDAY.LE.SEASN2) GO TO 120                                      305.   
  110 SEASN1=SEASN2                                                      306.   
      K=1                                                                307.   
      SEASN2=380.0                                                       308.   
  120 WT2=(XJDAY-SEASN1)/(SEASN2-SEASN1)                                 309.   
      WT1=1.-WT2                                                         310.   
      KS1=1+MOD(K,4)                                                     311.   
      KS2=1+MOD(K+1,4)                                                   312.   
      KN1=1+MOD(K+2,4)                                                   313.   
      KN2=K                                                              314.   
      DO 130 K=1,NV                                                      315.   
C------------------------                                                316.   
C     SOUTHERN HEMISPHERE                                                317.   
C------------------------                                                318.   
      ALVISD(K   )=WT1*ALVISK(K,KS1)+WT2*ALVISK(K,KS2)                   319.   
      ALNIRD(K   )=WT1*ALNIRK(K,KS1)+WT2*ALNIRK(K,KS2)                   320.   
C------------------------                                                321.   
C     NORTHERN HEMISPHERE                                                322.   
C------------------------                                                323.   
      ALVISD(K+NV)=WT1*ALVISK(K,KN1)+WT2*ALVISK(K,KN2)                   324.   
  130 ALNIRD(K+NV)=WT1*ALNIRK(K,KN1)+WT2*ALNIRK(K,KN2)                   325.   
c     print *,' ALVISK=',ALVISK
c     print *,' ALNIRK=',ALNIRK
c     print *,' ALVISD=',ALVISD
c     print *,' ALNIRD=',ALNIRD
      RETURN                                                             326.   
C                                                                        327.   
C----------------------------------------------------------------------- 328.   
C                 ALBEDO,THERMAL FLUX,FLUX DERIVATIVE FOR EACH SURF TYPE 329.   
      ENTRY GETALB                                                       330.   
C----------------------------------------------------------------------- 331.   
C                                                                        332.   
      LATHEM=NV                                                          333.   
      IF(JLAT.LT.JNORTH) LATHEM=0                                        334.   
c      print *,'From GETALB JLAT=',JLAT
c      print *,POCEAN,PEARTH,POICE,PLICE
C                                                                        335.   
C                                            -------------------------   336.   
C                                            SNOW ALBEDO SPECIFICATION   337.   
C                                            -------------------------   338.   
Ccc   ASNAGE=0.35*EXP(-0.2*AGESN)                                        339.   
      if(IFFF0.eq.1)then
       print *,' FRSNALB=',FRSNALB
	print *,' ASNALB(1)=',ASNALB(1),' ASNALB(2)=',ASNALB(2)
	print *,' AOIALB(1)=',AOIALB(1),' AOIALB(2)=',AOIALB(2)
	print *,' ALIALB(1)=',ALIALB(1),' ALIALB(2)=',ALIALB(2)
       IFFF0=0
      endif
      ASNAGE=FRSNALB*EXP(-0.2*AGESN)
      BSNVIS=ASNVIS+ASNAGE                                               340.   
      BSNNIR=ASNNIR+ASNAGE                                               341.   
      XSNVIS=BSNVIS                                                      342.   
      XSNNIR=BSNNIR                                                      343.   
C                                                                        344.   
      EXPSNE=1.                                                          345.   
      EXPSNO=1.                                                          346.   
      EXPSNL=1.                                                          347.   
C                                                                        348.   
      DO 200 I=1,16                                                      349.   
  200 BXA(I)=0.                                                          350.   
C                                                                        351.   
      DO 210 K=1,NKTR                                                    352.   
      TRGALB(K)=0.                                                       353.   
      BGFEMD(K)=0.                                                       354.   
  210 BGFEMT(K)=0.                                                       355.   
C                                                                        356.   
      BOCSUM=0.                                                          357.   
      BEASUM=0.                                                          358.   
      BOISUM=0.                                                          359.   
      BLISUM=0.                                                          360.   
C                                                                        361.   
      DO 220 K=1,4                                                       362.   
  220 DTRUFG(K)=0.                                                       363.   
C                                                                        364.   
C                                           --------------------------   365.   
C                                           OCEAN ALBEDO SPECIFICATION   366.   
C                                           --------------------------   367.   
      IF(POCEAN.LT.1.E-04) GO TO 400                                     368.   
      X=0.5+(0.5-COSZ)*ZOCSRA                                            369.   
      BOCVIS=BVH2O(WMAG)                                                 370.   
      XOCVIS=XVH2O(WMAG,X)                                               371.   
      BOCNIR=BOCVIS                                                      372.   
      XOCNIR=XOCVIS                                                      373.   
C                                                                        374.   
      X=1./(1.+WMAG)                                                     375.   
      AV=(-.0147087*X*X+.0292266*X-.0081079)*EOCTRA                      376.   
      BV=(1.01673-0.0083652*WMAG)*EOCTRA                                 377.   
C                                                                        378.   
      ITOC=TGO                                                           379.   
      WTOC=TGO-ITOC                                                      380.   
      ITOC=ITOC-IT0                                                      381.   
      BOCSUM=0.                                                          382.   
      BOCM=0.                                                            383.   
      BOCP=0.                                                            384.   
C                                                                        385.   
      DO 310 K=1,NKTR                                                    386.   
      TRAPOC=AV+BV*AOCEAN(K)                                             387.   
      BOCM1 =(PLANCK(ITOC-1)-(PLANCK(ITOC-1)-PLANCK(ITOC  ))*WTOC)       388.   
     +      *(1.-TRAPOC)                                                 389.   
      BOCM  =BOCM+BOCM1                                                  390.   
      BOCP1 =(PLANCK(ITOC+1)-(PLANCK(ITOC+1)-PLANCK(ITOC+2))*WTOC)       391.   
     +      *(1.-TRAPOC)                                                 392.   
      BOCP  =BOCP+BOCP1                                                  393.   
      BOC   =(PLANCK(ITOC  )-(PLANCK(ITOC  )-PLANCK(ITOC+1))*WTOC)       394.   
     +      *(1.-TRAPOC)                                                 395.   
      BOCSUM=BOCSUM+BOC                                                  396.   
      ITOC=ITOC+ITNEXT                                                   397.   
C                                                                        398.   
      TRGALB(K)=TRGALB(K)+POCEAN*TRAPOC                                  399.   
      BGFEMD(K)=BGFEMD(K)+POCEAN*(BOCP1-BOCM1)                           400.   
  310 BGFEMT(K)=BGFEMT(K)+POCEAN*BOC                                     401.   
      DTRUFG(1)=0.5*(BOCP-BOCM)                                          402.   
C                                        -----------------------------   403.   
C                                        SOIL/VEG ALBEDO SPECIFICATION   404.   
C                                        -----------------------------   405.   
  400 DSFRAC=PVT(1)                                                      406.   
      VGFRAC=1.-DSFRAC                                                   407.   
      IF(PEARTH.LT.1.E-04) GO TO 500                                     408.   
      IF(SNOWE .GT.1.E-04) GO TO 420                                     409.   
C                                                                        410.   
      BEAVIS=PVT(1)*ALVISD(1+LATHEM)*(1.0-0.5*WEARTH*WETSRA)             411.   
      BEANIR=PVT(1)*ALNIRD(1+LATHEM)*(1.0-0.5*WEARTH*WETSRA)             412.   
      DO 410 K=2,NV                                                      413.   
      BEAVIS=BEAVIS+PVT(K)*ALVISD(K+LATHEM)                              414.   
  410 BEANIR=BEANIR+PVT(K)*ALNIRD(K+LATHEM)                              415.   
      GO TO 440                                                          416.   
  420 VTFRAC=PVT(1)*EXP(-SNOWE/VTMASK(1))                                417.   
      EXPSNE=VTFRAC                                                      418.   
      DSFRAC=VTFRAC                                                      419.   
C$    BEAVIS=VTFRAC*ALVISD(1+LATHEM)*(1.-0.5*WEARTH*WETSRA)              420.   
      BEAVIS=PVT(1)*ALVISD(1+LATHEM)*(1.-0.5*WEARTH*WETSRA)              421.   
     +      *(1.-VTFRAC)                                                 422.   
C$    BEANIR=VTFRAC*ALNIRD(1+LATHEM)*(1.-0.5*WEARTH*WETSRA)              423.   
      BEANIR=PVT(1)*ALNIRD(1+LATHEM)*(1.-0.5*WEARTH*WETSRA)              424.   
     +      *(1.-VTFRAC)                                                 425.   
      DO 430 K=2,NV                                                      426.   
      VTFRAC=PVT(K)*EXP(-SNOWE/VTMASK(K))                                427.   
      BEAVIS=BEAVIS+PVT(K)*ALVISD(K+LATHEM)*(1.-VTFRAC)                  428.   
C$    BEAVIS=BEAVIS+VTFRAC*ALVISD(K+LATHEM)   *******************CORRECT 429.   
      BEANIR=BEANIR+PVT(K)*ALNIRD(K+LATHEM)*(1.-VTFRAC)                  430.   
C$    BEANIR=BEANIR+VTFRAC*ALNIRD(K+LATHEM)   *******************CORRECT 431.   
  430 EXPSNE=EXPSNE+VTFRAC                                               432.   
C                                                                        433.   
  440 XEAVIS=BEAVIS                                                      434.   
      XEANIR=BEANIR                                                      435.   
C$    BEAVIS=BEAVIS+BSNVIS*(1.-EXPSNE)                                   436.   
C$    BEANIR=BEANIR+BSNNIR*(1.-EXPSNE)                                   437.   
C$    XEAVIS=XEAVIS+XSNVIS*(1.-EXPSNE)                                   438.   
C$    XEANIR=XEANIR+XSNNIR*(1.-EXPSNE)                                   439.   
      BEAVIS=BEAVIS*EXPSNE+BSNVIS*(1.-EXPSNE)                            440.   
      BEANIR=BEANIR*EXPSNE+BSNNIR*(1.-EXPSNE)                            441.   
      XEAVIS=XEAVIS*EXPSNE+XSNVIS*(1.-EXPSNE)                            442.   
      XEANIR=XEANIR*EXPSNE+XSNNIR*(1.-EXPSNE)                            443.   
      VGFRAC=EXPSNE-DSFRAC                                               444.   

#if ( defined CLM )
c     if(ncallclm.ge.1)then
      i=1
      BEAVIS=0.7*asdirclm(i,JLAT)+0.3*asdifclm(i,JLAT)
      BEANIR=0.7*aldirclm(i,JLAT)+0.3*aldifclm(i,JLAT)
      XEAVIS=BEAVIS
      XEANIR=BEANIR
c     endif
c     if(ncallclm.eq.0)then
c       print *,JLAT,BEAVIS,BEANIR
c     endif
#endif

C                                                                        445.   
      ITEA=TGE                                                           446.   
      WTEA=TGE-ITEA                                                      447.   
      ITEA=ITEA-IT0                                                      448.   
      BEASUM=0.                                                          449.   
      BEAM=0.                                                            450.   
      BEAP=0.                                                            451.   
C                                                                        452.   
C                                                                        467.   
      DO 450 K=1,NKTR                                                    453.   
      TRAPEA=AGSIDV(K,1)*(1.-EXPSNE)                                     454.   
     +      +AGSIDV(K,3)*DSFRAC*(1.-WETTRA*WEARTH)                       455.   
     +      +AGSIDV(K,4)*VGFRAC                                          456.   
      BEAM1 =(PLANCK(ITEA-1)-(PLANCK(ITEA-1)-PLANCK(ITEA  ))*WTEA)       457.   
     +      *(1.-TRAPEA)                                                 458.   
      BEAM  =BEAM+BEAM1                                                  459.   
      BEAP1 =(PLANCK(ITEA+1)-(PLANCK(ITEA+1)-PLANCK(ITEA+2))*WTEA)       460.   
     +      *(1.-TRAPEA)                                                 461.   
      BEAP  =BEAP+BEAP1                                                  462.   
      BEA   =(PLANCK(ITEA  )-(PLANCK(ITEA  )-PLANCK(ITEA+1))*WTEA)       463.   
     +      *(1.-TRAPEA)                                                 464.   
      BEASUM=BEASUM+BEA                                                  465.   
      ITEA=ITEA+ITNEXT                                                   466.   
      TRGALB(K)=TRGALB(K)+PEARTH*TRAPEA                                  468.   
      BGFEMD(K)=BGFEMD(K)+PEARTH*(BEAP1-BEAM1)                           469.   
  450 BGFEMT(K)=BGFEMT(K)+PEARTH*BEA                                     470.   
      DTRUFG(2)=0.5*(BEAP-BEAM)                                          471.   
      if(ncallclm.eq.-1)then
        print *,'471 JLAT=',JLAT
      print *,(ITEA-1),(ITEA),(ITEA+1)
      print *,PLANCK(ITEA-1),PLANCK(ITEA),PLANCK(ITEA+1)
      print *,' VGFRAC=',VGFRAC,' DSFRAC=',DSFRAC
      print *,' WTEA=',WTEA,' WEARTH=',WEARTH
      print *,' SNOWE=',SNOWE,' EXPSNE=',EXPSNE
c     print *,JLAT,' BEAVIS=',BEAVIS,' BEANIR=',BEANIR
      endif
C                                                                        472.   
C                                       ------------------------------   473.   
C                                       OCEAN ICE ALBEDO SPECIFICATION   474.   
C                                       ------------------------------   475.   
  500 CONTINUE                                                           476.   
      IF(POICE.LT.1.E-04) GO TO 600                                      477.   
      EXPSNO=EXP(-SNOWOI/DMOICE)                                         478.   
      BOIVIS=AOIVIS*EXPSNO+BSNVIS*(1.-EXPSNO)                            479.   
      BOINIR=AOINIR*EXPSNO+BSNNIR*(1.-EXPSNO)                            480.   
Cnext lines - use albedos from ice model
#if ( defined OCEAN_3D )
      BOIVIS=mmsAlb(JLAT)
      BOINIR=mmsAlbNIR(JLAT)
#endif
      XOIVIS=BOIVIS                                                      481.   
      XOINIR=BOINIR                                                      482.   
C                                                                        483.   
      ITOI=TGOI                                                          484.   
      WTOI=TGOI-ITOI                                                     485.   
      ITOI=ITOI-IT0                                                      486.   
      BOISUM=0.                                                          487.   
      BOIM=0.                                                            488.   
      BOIP=0.                                                            489.   
C                                                                        490.   
      DO 510 K=1,NKTR                                                    491.   
      TRAPOI=AGSIDV(K,1)*ESNTRA*(1.-EXPSNO)                              492.   
     +      +AGSIDV(K,2)*EICTRA*EXPSNO                                   493.   
      BOIM1 =(PLANCK(ITOI-1)-(PLANCK(ITOI-1)-PLANCK(ITOI  ))*WTOI)       494.   
     +      *(1.-TRAPOI)                                                 495.   
      BOIM  =BOIM+BOIM1                                                  496.   
      BOIP1 =(PLANCK(ITOI+1)-(PLANCK(ITOI+1)-PLANCK(ITOI+2))*WTOI)       497.   
     +      *(1.-TRAPOI)                                                 498.   
      BOIP  =BOIP+BOIP1                                                  499.   
      BOI   =(PLANCK(ITOI  )-(PLANCK(ITOI  )-PLANCK(ITOI+1))*WTOI)       500.   
     +      *(1.-TRAPOI)                                                 501.   
      BOISUM=BOISUM+BOI                                                  502.   
      ITOI=ITOI+ITNEXT                                                   503.   
C                                                                        504.   
      TRGALB(K)=TRGALB(K)+POICE*TRAPOI                                   505.   
      BGFEMD(K)=BGFEMD(K)+POICE*(BOIP1-BOIM1)                            506.   
  510 BGFEMT(K)=BGFEMT(K)+POICE*BOI                                      507.   
      DTRUFG(3)=0.5*(BOIP-BOIM)                                          508.   
C                                                                        509.   
C                                        -----------------------------   510.   
C                                        LAND ICE ALBEDO SPECIFICATION   511.   
C                                        -----------------------------   512.   
  600 CONTINUE                                                           513.   
      IF(PLICE.LT.1.E-04) GO TO 700                                      514.   
      EXPSNL=EXP(-SNOWLI/DMLICE)                                         515.   
      BLIVIS=ALIVIS*EXPSNL+BSNVIS*(1.-EXPSNL)                            516.   
      BLINIR=ALINIR*EXPSNL+BSNNIR*(1.-EXPSNL)                            517.   

#if ( defined CLM )
c     if(ncallclm.ge.1)then
      i=1
      BLIVIS=0.7*asdirclm(i,JLAT)+0.3*asdifclm(i,JLAT)
      BLINIR=0.7*aldirclm(i,JLAT)+0.3*aldifclm(i,JLAT)
c     endif
#endif

      XLIVIS=BLIVIS                                                      518.   
      XLINIR=BLINIR                                                      519.   
C                                                                        520.   
      ITLI=TGLI                                                          521.   
      WTLI=TGLI-ITLI                                                     522.   
      ITLI=ITLI-IT0                                                      523.   
C                                                                        524.   
      BLISUM=0.                                                          525.   
      BLIM=0.                                                            526.   
      BLIP=0.                                                            527.   
      BGF=0.                                                             528.   
C                                                                        529.   
      DO 610 K=1,NKTR                                                    530.   
      TRAPLI=AGSIDV(K,1)*ESNTRA*(1.-EXPSNL)                              531.   
     +      +AGSIDV(K,2)*EICTRA*EXPSNL                                   532.   
      BLIM1 =(PLANCK(ITLI-1)-(PLANCK(ITLI-1)-PLANCK(ITLI  ))*WTLI)       533.   
     +      *(1.-TRAPLI)                                                 534.   
      BLIM  =BLIM+BLIM1                                                  535.   
      BLIP1 =(PLANCK(ITLI+1)-(PLANCK(ITLI+1)-PLANCK(ITLI+2))*WTLI)       536.   
     +      *(1.-TRAPLI)                                                 537.   
      BLIP  =BLIP+BLIP1                                                  538.   
      BLI   =(PLANCK(ITLI  )-(PLANCK(ITLI  )-PLANCK(ITLI+1))*WTLI)       539.   
     +      *(1.-TRAPLI)                                                 540.   
      BLISUM=BLISUM+BLI                                                  541.   
      ITLI=ITLI+ITNEXT                                                   542.   
C                                                                        543.   
      TRGALB(K)=TRGALB(K)+PLICE*TRAPLI                                   544.   
      BGFEMD(K)=BGFEMD(K)+PLICE*(BLIP1-BLIM1)                            545.   
  610 BGFEMT(K)=BGFEMT(K)+PLICE*BLI                                      546.   
      DTRUFG(4)=0.5*(BLIP-BLIM)                                          547.
C                                                                        548.   
  700 CONTINUE                                                           549.   
      BVSURF=POCEAN*BOCVIS +PEARTH*BEAVIS +POICE*BOIVIS +PLICE*BLIVIS    550.   
      XVSURF=POCEAN*XOCVIS +PEARTH*XEAVIS +POICE*XOIVIS +PLICE*XLIVIS    551.   
      BNSURF=POCEAN*BOCNIR +PEARTH*BEANIR +POICE*BOINIR +PLICE*BLINIR    552.   
      XNSURF=POCEAN*XOCNIR +PEARTH*XEANIR +POICE*XOINIR +PLICE*XLINIR    553.   

#if ( !defined CPL_CHEM ) && ( (defined SVI_ALBEDO || defined GHS_ALB) )
         IF(COSZ.GE.0.01) then
            XALBEDO=0.6*XVSURF+0.4*XNSURF
            SECZ=1./COSZ
            if(JLAT.le.-2)then
               print *,' JLAT=',JLAT
               print *,' COSZ=',COSZ
               print*,POCEAN,PEARTH,POICE,PLICE
               print *,' XALBEDO=',XALBEDO
               print *,BVSURF,XVSURF,BNSURF,XNSURF
            endif
            BVSURF=BVSURF+BVSURFA*(1.-XALBEDO)**2*SECZ
            XVSURF=XVSURF+XVSURFA*(1.-XALBEDO)**2*SECZ
            BNSURF=BNSURF+BNSURFA*(1.-XALBEDO)**2*SECZ
            XNSURF=XNSURF+XNSURFA*(1.-XALBEDO)**2*SECZ
            if(JLAT.eq.-10)then
               print *,' After add'
               print *,'BVSURFA=',BVSURFA
               print *,'DAsrf=',BVSURFA*(1.-XALBEDO)**2*SECZ
               print *,BVSURF,XVSURF,BNSURF,XNSURF
            endif
         endif
#endif

C     ----------------------------------------------------------------   554.   
C     SPECTRAL DISTRIBUTION ASSUMES THAT:  AMEAN = 0.6*AVIS + 0.4*ANIR   555.   
C     ----------------------------------------------------------------   556.   
C                                                                        557.   
      IF(KEEPAL.EQ.1) GO TO 800                                          558.   
      SRBALB(6)=BVSURF+0.4*VISNIR*(BNSURF-BVSURF)                        559.   
      SRXALB(6)=XVSURF+0.4*VISNIR*(XNSURF-XVSURF)                        560.   
      DO 710 I=1,5                                                       561.   
      SRBALB(I)=BNSURF-0.6*VISNIR*(BNSURF-BVSURF)                        562.   
  710 SRXALB(I)=XNSURF-0.6*VISNIR*(XNSURF-XVSURF)                        563.   
      IF(KALVIS.EQ.0) GO TO 800                                          564.   
      SRBALB(4)=SRBALB(6)                                                565.   
      SRXALB(4)=SRXALB(6)                                                566.   
C                                                                        567.   
C--------------------------------------------------------------------    568.   
C     DEFINE SURFACE FLUX FACTORS, FLUX DERIVATIVES FOR EACH SURFTYPE    569.   
C--------------------------------------------------------------------    570.   
  800 BGF=0.                                                             571.   
      DO 810 K=1,NKTR                                                    572.   
      BGFEMD(K)=BGFEMD(K)*0.5                                            573.   
  810 BGF=BGF+BGFEMT(K)                                                  574.   
C                                                                        575.   
      BGM=BOCM*POCEAN+BEAM*PEARTH+BOIM*POICE+BLIM*PLICE                  576.   
      BGP=BOCP*POCEAN+BEAP*PEARTH+BOIP*POICE+BLIP*PLICE                  577.   
      TTRUFG=0.5*(BGP-BGM)                                               578.   
C                                                                        579.   
      FTRUFG(1)=BOCSUM/BGF                                               580.   
      FTRUFG(2)=BEASUM/BGF                                               581.   
      FTRUFG(3)=BOISUM/BGF                                               582.   
      FTRUFG(4)=BLISUM/BGF                                               583.   
C                                                                        584.   
      RETURN                                                             585.   
      END                                                                586.   
c     SUBROUTINE SETGAS                                                  587.   
c    20/06/2005
      SUBROUTINE SETGAS(KTREND)

#include "B83XX.COM"                                                     588.   
#include "chem_para"
#include "chem_com"

      EQUIVALENCE (IMG(1),IMGAS1),(IMG(2),IMGAS2)                        649.   
      EQUIVALENCE (ILG(1),ILGAS1),(ILG(2),ILGAS2)                        650.   
C                                                                        651.   
C                                                                        652.   
C----------------------------------------------------------------------  653.   
C     GLOBAL   U.S. (1976) STANDARD ATMOSPHERE   P, T, GEO H  PARAMETERS 654.   
C----------------------------------------------------------------------  655.   
C                                                                        656.   
      COMMON/O3GLOB/ PLB0(40),TLM0(40),U0GAS3(40)                        656.11 
      DIMENSION SPLB(8),STLB(8),SHLB(8),SDLB(8)                          657.   
      DATA SPLB/1013.25,226.32,54.748,8.6801,1.109,.66938,.039564        658.   
     +         ,3.7338E-03/                                              659.   
      DATA STLB/288.15,216.65,216.65,228.65,270.65,270.65,214.65,186.87/ 660.   
      DATA SHLB/0.0,11.0,20.0,32.0,47.0,51.0,71.0,84.852/                661.   
      DATA SDLB/-6.5,0.0,1.0,2.8,0.0,-2.8,-2.0,0.0/                      662.   
      DATA HPCON/34.16319/                                               663.   
      DATA PI/3.1415926/                                                 664.   
      DATA P0/1013.25/                                                   665.   
C                                                                        666.   
      DIMENSION KGAS(9,3)                                                667.   
      DATA KGAS/ 1, 2, 3, 0, 0, 9, 11, 12, 13                            668.   
     +         , 4, 6, 8, 0, 0,10,  0,  0,  0                            669.   
     +         , 5, 7, 0, 0, 0, 0,  0,  0,  0/                           670.   
C                                                                        671.   
C                -----------------------------------------------------   672.   
C                USE PLB TO FIX STANDARD HEIGHTS FOR GAS DISTRIBUTIONS   673.   
C                -----------------------------------------------------   674.   
C                                                                        675.   
c     print *,'FROM SETGAS PREDICTED_GASES=',PREDICTED_GASES
c     6/20/2005
      if(KTREND.le.0)then
C     assign background GHGs
      PPMV58(2)=GHGBGR(1)	! CO2
      PPMV58(6)=GHGBGR(2)	! N2O
      PPMV58(7)=GHGBGR(3)	! CH4
      PPMV58(8)=GHGBGR(4)	! F11
      PPMV58(9)=GHGBGR(5)	! F12
      endif
      print *,'PPMV58 from SETGAS'
       print *,PPMV58
      NLP=NL+1                                                           676.   
      NLMOD=NLP-LAYRAD                                                   677.   
      PS0=PLB(1)                                                         678.   
      PTOP=PLB(NLP-LAYRAD)                                               679.   
C                                                                        680.   
      DO 100 L=1,NL                                                      681.   
      DPL(L)=PLB(L)-PLB(L+1)                                             682.   
  100 PL(L)=(PLB(L)+PLB(L+1))*0.5                                        683.   
      NLNKTR=NL*NKTR                                                     684.   
C                                                                        685.   
      IF(LASTVC.GE.0) GO TO 107                                          686.   
C                                                                        687.   
      DO 105 L=1,NL                                                      688.   
      P=PLB(L)                                                           689.   
      DO 101 N=2,8                                                       690.   
      IF(P.GT.SPLB(N)) GO TO 102                                         691.   
  101 CONTINUE                                                           691.5  
      N=9                                                                692.   
  102 N=N-1                                                              693.   
      IF(ABS(SDLB(N)).LT.1.E-04) GO TO 103                               694.   
      H=SHLB(N)+STLB(N)/SDLB(N)*((SPLB(N)/P)**(SDLB(N)/HPCON)-1.)        695.   
      GO TO 104                                                          696.   
C ALOG
  103 H=SHLB(N)+STLB(N)/HPCON*LOG(SPLB(N)/P)                             697.   
C ALOG
  104 T=STLB(N)+SDLB(N)*(H-SHLB(N))                                      698.   
      TLB(L)=T                                                           699.   
  105 HLB(L)=H                                                           700.   
!     print *,' After 105'
      HLB(1)=1.E-10                                                      701.   
      HLB(NL+1)=99.99                                                    702.   
      TLB(NL+1)=STLB(8)                                                  703.   
      DO 106 L=1,NL                                                      704.   
      TLT(L)=TLB(L+1)                                                    705.   
  106 TLM(L)=0.5*(TLB(L)+TLT(L))                                         706.   
      TLB(NL+1)=TLT(NL)                                                  707.   
C                                                                        708.   
  107 NLAY=LASTVC/100000                                                 709.   
      NATM=(LASTVC-NLAY*100000)/10000                                    710.   
      IF(NATM.GT.0) GO TO 112                                            711.   
C                                                                        712.   
C---------------------------------------------------------------------   713.   
C     DEFINE GLOBAL MEAN GAS AMOUNTS FOR TRACEGAS & OVERLAP ABSORPTION   714.   
C---------------------------------------------------------------------   715.   
C                                                                        716.   
C                                         ----------------------------   717.   
C                                         GLOBAL MEAN H2O DISTRIBUTION   718.   
C                                         ----------------------------   719.   
      RHP=0.77                                                           720.   
      EST=10.0**(9.4051-2353.0/TLB(1))                                   721.   
      FWB=0.662*RHP*EST/(PLB(1)-RHP*EST)                                 722.   
      DO 111 L=1,NL                                                      723.   
      PLT=PLB(L+1)                                                       724.   
      DP=PLB(L)-PLT                                                      725.   
      RHP=0.77*(PLT/P0-0.02)/.98                                         726.   
      EST=10.0**(9.4051-2353.0/TLT(L))                                   727.   
      FWT=0.662*RHP*EST/(PLT-RHP*EST)                                    728.   
      IF(FWT.GT.3.E-06) GO TO 110                                        729.   
      FWT=3.E-06                                                         730.   
      RHP=FWT*PLT/(EST*(FWT+0.662))                                      731.   
  110 ULGASL=0.5*(FWB+FWT)*DP*1270.                                      732.   
C$110 ULGASL=0.5*(FWB+FWT)*DP*1268.75                                    733.   
      U0GAS(L,1)=ULGASL                                                  734.   
      SHL(L)=ULGASL/(ULGASL+1268.75*DP)                                  735.   
      EQ=0.5*(PLB(L)+PLT)*SHL(L)/(0.662+0.378*SHL(L))                    736.   
      ES=10.**(9.4051-2353./TLM(L))                                      737.   
      RHL(L)=EQ/ES                                                       738.   
  111 FWB=FWT                                                            739.   
  112 CONTINUE                                                           740.   
C                                         ----------------------------   741.   
C                                         GLOBAL MEAN  O3 DISTRIBUTION   742.   
C----------------                         ----------------------------   743.   
!     print *,' Before SETO3D'
      CALL SETO3D                                                        744.   
!     print *,' After SETO3D'
C----------------                                                        745.   
      JJLAT=JLAT                                                         746.   
C                             IF(JDAY.LT.1)  KEEP SETATM O3 DISTRIBUTION 747.   
C                             ------------------------------------------ 748.   
      IF(JDAY.LT.1) GO TO 125                                            749.   
C----------------                                                        750.   
!     print *,' Before O3DDAY'
      CALL O3DDAY                                                        751.   
!     print *,' After O3DDAY'
C----------------                                                        752.   
C                                                                        753.   
      DO 120 J=1,JMLAT                                                   754.   
      RADLAT=PI*DLAT(J)/180.                                             755.   
  120 COSLAT(J)=0.5+0.5*SIN(RADLAT)                                      756.   
C                                                                        757.   
      DO 121 N=1,NL                                                      758.   
  121 UO3L(N)=0.                                                         759.   
      DO 123 JLAT=1,JMLAT                                                760.   
C----------------                                                        761.   
!     print *,' Before O3DLAT'
      CALL O3DLAT                                                        762.   
!     print *,' After O3DLAT'
C----------------                                                        763.   
      JB=JLAT+1                                                          764.   
      JA=JLAT-1                                                          765.   
      IF(JB.GT.JMLAT) JB=JMLAT                                           766.   
      IF(JA.LT.1    ) JA=1                                               767.   
      WTLAT=0.5*(COSLAT(JB)-COSLAT(JA))                                  768.   
      DO 122 N=1,NL                                                      769.   
  122 UO3L(N)=UO3L(N)+U0GAS(N,3)*WTLAT                                   770.   
  123 CONTINUE                                                           771.   
      DO 124 N=1,NL                                                      772.   
  124 U0GAS(N,3)=UO3L(N)                                                 773.   
  125 JLAT=JJLAT                                                         774.   
!     print *,' After 774'
      XXXX=SETAO3(OCM)                                                   775.   
!     print *,' After 775'
C                                                                        775.11 
C     SAVE GLOBAL MEAN P,T,O3 FOR UPDATING LAPGAS TAU TABLE IN SETLAP    775.12 
C     ---------------------------------------------------------------    775.13 
C                                                                        775.14 
      DO 126 N=1,NL                                                      775.15 
      PLB0(N)=PLB(N)                                                     775.16 
      TLM0(N)=TLM(N)                                                     775.17 
  126 U0GAS3(N)=U0GAS(N,3)                                               775.18 
      PLB0(NLP)=PLB(NLP)                                                 775.19 
C                                         ----------------------------   776.   
C                                         GLOBAL MEAN NO2 DISTRIBUTION   777.   
C                                         ----------------------------   778.   
!     print *,' After 778'
      ACM=0.0                                                            779.   
      HI=0.0                                                             780.   
      FI=CMANO2(1)                                                       781.   
      HL=HLB(2)                                                          782.   
      L=1                                                                783.   
      J=1                                                                784.   
  130 J=J+1                                                              785.   
      IF(J.GT.42) GO TO 133                                              786.   
      HJ=HI+2.0                                                          787.   
      FJ=CMANO2(J)                                                       788.   
  131 DH=HJ-HI                                                           789.   
      IF(HJ.GT.HL) GO TO 132                                             790.   
      ACM=ACM+(FI+FJ)*DH*0.5                                             791.   
      HI=HJ                                                              792.   
      FI=FJ                                                              793.   
      GO TO 130                                                          794.   
  132 FF=FI+(FJ-FI)*(HL-HI)/DH                                           795.   
      DH=HL-HI                                                           796.   
      ACM=ACM+(FI+FJ)*DH*0.5                                             797.   
      U0GAS(L,5)=ACM                                                     798.   
      ACM=0.0                                                            799.   
      HI=HL                                                              800.   
      FI=FF                                                              801.   
      IF(L.EQ.NL) GO TO 133                                              802.   
      L=L+1                                                              803.   
      HL=HLB(L+1)                                                        804.   
      GO TO 131                                                          805.   
  133 U0GAS(L,5)=ACM                                                     806.   
      ACM=0.0                                                            807.   
      L=L+1                                                              808.   
      IF(L.LT.NLP) GO TO 133                                             809.   
!     print *,' After 809'
C                            -----------------------------------------   810.   
C                            (CO2,O2) UNIFORMLY MIXED GAS DISTRIBUTION   811.   
C                            -----------------------------------------   812.   
      DO 141 K=2,4,2                                                     813.   
      DO 140 N=1,NL                                                      814.   
  140 U0GAS(N,K)=PPMV58(K)*0.8*DPL(N)/P0                                 815.   
  141 CONTINUE                                                           816.   
C  PRINT
       print *,' CO2',PPMV58(2)
c     print *,'NLMOD=',NLMOD
c     print *,'PSIG'
c     print *,(PSIG(L),L=1,NLMOD+1)
c     print *,'PLB'
c     print *,(PLB(L),L=1,NLMOD+1)
c      print *,(U0GAS(n,2),n=1,nl)
C PRINT
C                -----------------------------------------------------   817.   
C                (N20,CH4,F11,F12) SPECIFIED VERTICAL GAS DISTRIBUTION   818.   
C                -----------------------------------------------------   819.   
      DO 151 K=6,9                                                       820.   
      DO 150 N=1,NL                                                      821.   
      U0GAS(N,K)=PPMV58(K)*0.8*DPL(N)/P0                                 822.   
      ZT=(HLB(N+1)-Z0(K))/ZH(K)                                          823.   
      IF(ZT.LE.0.) GO TO 150                                             824.   
      ZB=(HLB(N)-Z0(K))/ZH(K)                                            825.   
      EXPZT=EXP(-ZT)                                                     826.   
      EXPZB=EXP(-ZB)                                                     827.   
      IF(ZB.LT.0.) EXPZB=1.-ZB                                           828.   
      U0GAS(N,K)=U0GAS(N,K)*(EXPZB-EXPZT)/(ZT-ZB)                        829.   
  150 CONTINUE                                                           830.   
  151 CONTINUE                                                           831.   
C                     ------------------------------------------------   832.   
C                     SPECIFIED GAS AMOUNTS (INCLUDING SCALING FACTOR)   833.   
C                     ------------------------------------------------   834.   
C                                                                        835.   
      DO 161 K=1,9                                                       836.   
      DO 160 N=1,NL                                                      837.   
  160 ULGAS(N,K)=U0GAS(N,K)*FULGAS(K)                                    838.   
  161 CONTINUE                                                           839.   
C PRINT 
!     print *,' after 161'
!     print *,(ULGAS(n,2),n=1,nl)
C PRINT
C                                                                        840.   
C-------------------------------                                         841.   
      CALL SETAO2(ULGAS(1,4),NL)                                         842.   
C-------------------------------                                         843.   
C                                                                        844.   
C       --------------------------------------------------------------   845.   
C       OVERLAP ABSORPTION (ILGAS1,ILGAS2) FOR GLOBAL MEAN GAS AMOUNTS   846.   
C       --------------------------------------------------------------   847.   
      DO 170 K=1,30                                                      848.   
  170 MLGAS(K)=0                                                         849.   
      IF(LAPGAS.LT.1) GO TO 174                                          850.   
      DO 172 L=1,3                                                       851.   
      DO 171 K=ILGAS1,ILGAS2                                             852.   
      M=KGAS(K,L)                                                        853.   
      IF(M.GT.3) MLGAS(M)=1                                              854.   
  171 CONTINUE                                                           855.   
  172 CONTINUE                                                           856.   
      DO 173 K=1,15                                                      857.   
  173 MLGAS(15+K)=MLGAS(K)                                               858.   
  174 CONTINUE                                                           859.   
C                                                                        860.   
C     ----------------------------------------------------------------   861.   
C     TAULAP=OVERLAP ABSORPTION KEPT AS INITIALIZED (NO CHANGES LATER)   862.   
C     ----------------------------------------------------------------   863.   
C                                                                        864.   
      DO 180 I=1,1000                                                    865.   
      TAULAP(I)=0.                                                       866.   
  180 TAUN(I)=0.                                                         867.   
C                                                                        868.   
C--------------------------------                                        869.   
      IF(LAPGAS.GT.0) CALL TAUGAS                                        870.   
C--------------------------------                                        871.   
C                                                                        872.   
      DO 181 I=1,NLNKTR                                                  873.   
  181 TAULAP(I)=TAUN(I)                                                  874.   
C                                                                        875.   
C           ----------------------------------------------------------   876.   
C           MAIN GAS (IMGAS1,IMGAS2) ABSORPTION INTERPOLATED AS NEEDED   877.   
C           ----------------------------------------------------------   878.   
C                                                                        879.   
      DO 191 L=1,3                                                       880.   
      DO 190 K=IMGAS1,IMGAS2                                             881.   
      M=KGAS(K,L)                                                        882.   
      IF(M.GT.0) MLGAS(M)=1                                              883.   
  190 CONTINUE                                                           884.   
  191 CONTINUE                                                           885.   
      DO 192 K=1,13                                                      886.   
  192 MLGAS(K)=MLGAS(K)*(MLGAS(K)-MLGAS(K+15))                           887.   
      IF(IMGAS1.EQ.1) MLGAS(14)=1                                        888.   
      IF(KWVCON.EQ.1) MLGAS(15)=1                                        889.   
      DO 193 K=1,30                                                      890.   
  193 MLLAP(K)=MLGAS(K)                                                  891.   
C                                                                        892.   
      RETURN                                                             893.   
C                                                                        894.   
C----------------------------------------------------------------------- 895.   
C                REDEFINE TAULAP TABLE: GET ABSORPTION FROM TAUGAS TABLE 896.   
      ENTRY SETLAP                                                       897.   
C----------------------------------------------------------------------- 898.   
C                                                                        899.   
      IF(LAPGAS.EQ.1) RETURN                                             900.   
C                                                                        901.   
      DO 200 I=1,1000                                                    902.   
  200 TAULAP(I)=0.                                                       903.   
      IF(LAPGAS.EQ.0) RETURN                                             904.   
C                                                                        905.   
      DO 210 K=1,15                                                      906.   
  210 MLGAS(K)=MLLAP(K+15)                                               907.   
C                                                                        908.   
      DO 220 I=1,NLNKTR                                                  909.   
  220 TAUN(I)=TAULAP(I)                                                  910.   
C                                                                        911.   
      DO 230 L=1,NL                                                      912.   
      DPL(L)=PLB0(L)-PLB0(L+1)                                           912.11 
      PL(L)=(PLB0(L)+PLB0(L+1))*0.5                                      912.12 
      TLM(L)=TLM0(L)                                                     912.13 
      U0GAS(L,3)=U0GAS3(L)                                               912.14 
C                                                                        912.15 
      ULGAS(L,1)=U0GAS(L,1)*FULGAS(1)                                    913.   
      ULGAS(L,3)=U0GAS(L,3)*FULGAS(3)                                    914.   
  230 ULGAS(L,5)=U0GAS(L,5)*FULGAS(5)                                    915.   
C                                                                        916.   
c
      tropmass = 28.97296245*1.e-3*0.8/P0
      trpm=tropmass*1.e3
      DO 240 L=1,nlev
      ULGAS(L,4)=U0GAS(L,4)*FULGAS(4)

#if ( defined PREDICTED_GASES  && !defined FIXED_FOR )
        pxxx = dpl(l)

        ULGAS(L,2)=glbgas(l,1)*tropmass/44.0098
     &            *pxxx
        ULGAS(L,6)=glbgas(l,2)*tropmass/44.0000
     &            *pxxx
        ULGAS(L,7)=glbgas(l,3)*tropmass/16.0426
     &            *pxxx
        ULGAS(L,8)=glbgas(l,4)*tropmass/137.3675
     &            *pxxx
        ULGAS(L,9)=glbgas(l,5)*tropmass/120.9054
     &            *pxxx
#else
	!
	!prescribed greenhouse
	!  gas profiles
	!
      ULGAS(L,2)=U0GAS(L,2)*FULGAS(2)                                    918.   
      ULGAS(L,6)=U0GAS(L,6)*FULGAS(6)                                    920.   
      ULGAS(L,7)=U0GAS(L,7)*FULGAS(7)                                    921.   
      ULGAS(L,8)=U0GAS(L,8)*FULGAS(8)                                    922.   
      ULGAS(L,9)=U0GAS(L,9)*FULGAS(9)
#endif
 240  continue
      ll=nlev
      do 2240 l=nlev+1,NL
      ULGAS(L,4)=U0GAS(L,4)*FULGAS(4)
#if ( defined PREDICTED_GASES  && !defined FIXED_FOR )
        pxxx = dpl(l)

        ULGAS(L,2)=glbgas(ll,1)*tropmass/44.0098
     &            *pxxx
        ULGAS(L,6)=glbgas(ll,2)*tropmass/44.0000
     &            *pxxx
        ULGAS(L,7)=glbgas(ll,3)*tropmass/16.0426
     &            *pxxx
        ULGAS(L,8)=glbgas(ll,4)*tropmass/137.3675
     &            *pxxx
        ULGAS(L,9)=glbgas(ll,5)*tropmass/120.9054
     &            *pxxx
#else 
      ULGAS(L,2)=U0GAS(L,2)*FULGAS(2)    
      ULGAS(L,6)=U0GAS(L,6)*FULGAS(6)   
      ULGAS(L,7)=U0GAS(L,7)*FULGAS(7)           
      ULGAS(L,8)=U0GAS(L,8)*FULGAS(8)         
      ULGAS(L,9)=U0GAS(L,9)*FULGAS(9)
#endif
 2240  continue
C PRINT
c     print *,' after 240'
c     print *,(ULGAS(n,2),n=1,nl)
C PRINT
C                                                                        924.   
C-----------------                                                       925.   
      CALL  TAUGAS                                                       926.   
C-----------------                                                       927.   
C                                                                        928.   
      DO 250 I=1,NLNKTR                                                  929.   
  250 TAULAP(I)=TAUN(I)                                                  930.   
C                                                                        931.   
      DO 260 K=1,15                                                      932.   
  260 MLGAS(K)=MLLAP(K)                                                  933.   
C                                                                        934.   
      RETURN                                                             935.   
C                                                                        936.   
C----------------------------------------------------------------------- 937.   
C                SPECIFY ULGAS: GET MAINGAS ABSORPTION FROM TAUGAS TABLE 938.   
      ENTRY GETGAS                                                       939.   
C----------------------------------------------------------------------- 940.   
C                                                                        941.   
C-----------------                                                       942.   
      CALL  O3DLON                                                       943.   
C-----------------                                                       944.   
C                                                                        945.   
      DO 300 L=1,NL                                                      946.   
      DPL(L)=PLB(L)-PLB(L+1)                                             947.   
  300 PL(L)=(PLB(L)+PLB(L+1))*0.5                                        948.   
C                                                                        949.   
      IF(KEEPRH.EQ.1) GO TO 311                                          950.   
      DO 310 L=1,NL                                                      951.   
  310 U0GAS(L,1)=1268.75*DPL(L)*SHL(L)                                   952.   
C$310 U0GAS(L,1)=1268.75*DPL(L)*SHL(L)/(1.-SHL(L))   ********CORRECT     953.   
      GO TO 313                                                          954.   
  311 CONTINUE                                                           955.   
      DO 312 L=1,NL                                                      956.   
      ES=10.0**(9.4051-2353.0/TLM(L))                                    957.   
      SHL(L)=0.622*(RHL(L)*ES)/(PL(L)-0.378*(RHL(L)*ES))                 958.   
  312 U0GAS(L,1)=1268.75*DPL(L)*SHL(L)                                   959.   
C$312 U0GAS(L,1)=1268.75*DPL(L)*SHL(L)/(1.-SHL(L))   ********CORRECT     960.   
  313 CONTINUE                                                           961.   
C                                                                        962.   
      DO 320 I=1,NLNKTR                                                  963.   
  320 TAUN(I)=TAULAP(I)                                                  964.   
C                                                                        965.   
      DO 330 L=1,NL                                                      966.   
      ULGAS(L,1)=U0GAS(L,1)*FULGAS(1)                                    967.   
      ULGAS(L,3)=U0GAS(L,3)*FULGAS(3)                                    968.   
  330 ULGAS(L,5)=U0GAS(L,5)*FULGAS(5)                                    969.   
C                                                                        970.   
      PART=(PLB(1)-PTOP)/(PS0-PTOP)                                      971.   

!
! --- Chemistry model patch  080895
!
! --- Note: most of the modifications in following 
!	sections were made originally as a part of chemistry
!	module ( PREDICTED_GASES == CPL_CHEM ). However,
!	they can be used by non-interactive
!	chemistry-climate runs now, as far as the prescribed
!	profiles of chemical species and aerosols are
!	available.
!
!	Chien Wang
!	080100
!

c ===
c Prescribed gaseous profiles:
c
c      DO 340 L=1,NL                                                      972.   
c      IF(L.EQ.NLMOD) PART=1.                                             973.   
c      ULGAS(L,2)=U0GAS(L,2)*FULGAS(2)*PART                               974.   
c      ULGAS(L,4)=U0GAS(L,4)*FULGAS(4)*PART                               975.   
c      ULGAS(L,6)=U0GAS(L,6)*FULGAS(6)*PART                               976.   
c      ULGAS(L,7)=U0GAS(L,7)*FULGAS(7)*PART                               977.   
c      ULGAS(L,8)=U0GAS(L,8)*FULGAS(8)*PART                               978.   
c      ULGAS(L,9)=U0GAS(L,9)*FULGAS(9)*PART                               979.   
c340   continue
c	goto 9341
c
c ===

!
! --- Use predicted gaseous profiles:
!
      tropmass = 28.97296245*1.e-3*0.8/P0
      trpm=tropmass*1.e3

!
! --- Use internal point to avoid possible unstable
! ---     problem related to LBC:
!
	jyyy = max(3, min(nlat2,JLAT))
!

      do 2340 l=1,nlev
	IF(L.EQ.NLMOD) PART=1.

      ULGAS(L,4)=U0GAS(L,4)*FULGAS(4)*PART	!O2

#if ( defined PREDICTED_GASES  && !defined FIXED_FOR )
!
! --- predicted greenhouse gas profiles
!  
	pxxx = dpl(l)*part

c       if (JLAT.eq.12) then 
c         print *,'zco2=',zco2(1,jlat,l)
c       endif
	ULGAS(L,2)=dmax1(0.0,zco2(ILON,jyyy,l))/44.0098
     &            *pxxx*tropmass
! --- Chem adjustmen of CO2 for IPCC forcing
     &            *adjfactor(1)
!       if (JLAT.eq.14) then 
!         print *,'l=',L,' ULGAS(L,2)=',ULGAS(L,2),
!    &              zco2(ILON,jyyy,l)
!       endif
     
#ifdef O3_RAD
!
!	=== Chien Wang 121797 then 062498 ===
!	===  add to use predicted ozone   ===
!	===  in troposphere only          ===
	if(l.le.n_tropopause)
     &	ULGAS(L,3)=dmax1(0.0,o3(ILON,jyyy,l))/48.0
     &            *pxxx*tropmass
#endif

!       if (JLAT.eq.14) then 
!         print *,'l=',L,' ULGAS(L,3)=',ULGAS(L,3)
!       endif
   
! --- Chem adjustmen of N2O and CH4 concentrations
!
        xxxo=dmax1(0.0,xn2o(ILON,jyyy,l))
     &      *tropmass/44.0000*1.25*P0
        yyyo=dmax1(0.0,ch4(ILON,jyyy,l))
     &      *tropmass/16.0426*1.25*P0
!       call newcon(xxxn,yyyn,xxxo,yyyo,ppmv58(6),ppmv58(7))
         xxxn=adjfactor(2)*xxxo
         yyyn=adjfactor(3)*yyyo

        ULGAS(L,6)=xxxn*0.8*pxxx/P0
        ULGAS(L,7)=yyyn*0.8*pxxx/P0

#if ( defined INC_3GASES  && !defined FIXED_FOR )
!
! === if hfc, pfc, and sf6 are included:
!
! === 032698
! === add hfc134a, pfc and sf6 to equivilent f11:
! ===
	  equi_cfc11 = cfc11(ILON,jyyy,l)
     &	  	     + hfc134a(ilon,jyyy,l)*dhfc134a_df11
     &		     + pfc    (ilon,jyyy,l)*dpfmethane_df11
     &		     + sf6    (ilon,jyyy,l)*dsf6_df11
#else 
	  equi_cfc11 = cfc11(ILON,jyyy,l)
#endif
	! ===
        ULGAS(L,8)=dmax1(0.0,equi_cfc11)
     &            *tropmass/137.3675
     &            *pxxx
        ULGAS(L,9)=dmax1(0.0,cfc12(ILON,jyyy,l))
     &            *tropmass/120.9054
     &            *pxxx
     
#else		
!
! --- prescribed greenhouse gas profiles
!  
      ULGAS(L,2)=U0GAS(L,2)*FULGAS(2)*PART
      ULGAS(L,6)=U0GAS(L,6)*FULGAS(6)*PART
      ULGAS(L,7)=U0GAS(L,7)*FULGAS(7)*PART
      ULGAS(L,8)=U0GAS(L,8)*FULGAS(8)*PART
      ULGAS(L,9)=U0GAS(L,9)*FULGAS(9)*PART
#endif

#if ( !defined CPL_CHEM ) && ( defined O3_FORCED )
C-- Added ozone forcing from external source. 
C--   changed 18Mar98 CEForest
C     NB. ozone is updated daily
C      o3 = ppb(m)
C      48 = mol weight of o3
C      ULGAS = cm^3 (STP)/cm^2
C
C     15JAN03  CEForest
C       changed to use total ozone, rather than anomalies, from GISS data
C
      pxxx = dpl(l)*part
      ULGAS(L,3)= o3dev(ILON,jlat,l)/48.0
     &       *pxxx*tropmass
C              if (l.eq.9) print *,'Ozone Dev: ',jlat,o3dev(ILON,jlat,l)
C--   end of change 18Mar98
#endif

 2340 continue

      ll=nlev
      do 2342 l=nlev+1,NL
	IF(L.EQ.NLMOD) PART=1.

      ULGAS(L,4)=U0GAS(L,4)*FULGAS(4)*PART	!O2

#if ( defined PREDICTED_GASES  && !defined FIXED_FOR )
!
! --- predicted greenhouse gas profiles
!  
	pxxx = dpl(l)*part

        ULGAS(L,2)=dmax1(0.0,zco2(ILON,jyyy,ll))
     &            *tropmass/44.0098
     &            *pxxx
!
! --- Chem adjustmen of N2O and CH4 concentrations
!
        xxxo=dmax1(0.0,xn2o(ILON,jyyy,ll))
     &      *tropmass/44.0000*1.25*P0
        yyyo=dmax1(0.0,ch4(ILON,jyyy,ll))
     &      *tropmass/16.0426*1.25*P0
!       call newcon(xxxn,yyyn,xxxo,yyyo,ppmv58(6),ppmv58(7))
         xxxn=adjfactor(2)*xxxo
         yyyn=adjfactor(3)*yyyo

        ULGAS(L,6)=xxxn*0.8*pxxx/P0
        ULGAS(L,7)=yyyn*0.8*pxxx/P0

#if ( defined INC_3GASES  && !defined FIXED_FOR )
!
! === if hfc, pfc, and sf6 are included:
!
! === 032698
! === add hfc134a, pfc and sf6 to equivilent f11:
! ===
	  equi_cfc11 = cfc11(ILON,jyyy,ll)
     &		     + hfc134a(ilon,jyyy,ll)*dhfc134a_df11
     &		     + pfc    (ilon,jyyy,ll)*dpfmethane_df11
     &		     + sf6    (ilon,jyyy,ll)*dsf6_df11
#else
	  equi_cfc11 = cfc11(ILON,jyyy,ll)
#endif
        ULGAS(L,8)=dmax1(0.0,equi_cfc11)
     &            *tropmass/137.3675
     &            *pxxx
        ULGAS(L,9)=dmax1(0.0,cfc12(ILON,jyyy,ll))
     &            *tropmass/120.9054
     &            *pxxx
#else
!
! --- prescribed greenhouse gas profiles
!
      ULGAS(L,2)=U0GAS(L,2)*FULGAS(2)*PART
      ULGAS(L,6)=U0GAS(L,6)*FULGAS(6)*PART
      ULGAS(L,7)=U0GAS(L,7)*FULGAS(7)*PART
      ULGAS(L,8)=U0GAS(L,8)*FULGAS(8)*PART
      ULGAS(L,9)=U0GAS(L,9)*FULGAS(9)*PART
#endif

#if ( !defined CPL_CHEM ) && ( defined O3_FORCED )
C-- Added ozone forcing from external source. 
C--   changed 18Mar98 CEForest
C     NB. ozone is updated daily
C      o3 = ppb(m)
C      48 = mol weight of o3
C      ULGAS = cm^3 (STP)/cm^2
C
C     15JAN03  CEForest
C       changed to use total ozone, rather than anomalies, from GISS data
C
C       added adjustment to layers (nlev+1:nlev+3) above dynamics layers
      pxxx = dpl(l)*part
      ULGAS(L,3)= o3dev(ILON,jlat,l)/48.0
     &        *pxxx*tropmass
C              if (l.eq.9) print *,'Ozone Dev: ',jlat,o3dev(ILON,jlat,l)
C--   end of change 18Mar98
#endif


 2342 continue

c
c-------------------------------------------------------

C-----------------                                                       981.   
      CALL  TAUGAS                                                       982.   
C-----------------                                                       983.   
C                                                                        984.   
      RETURN                                                             985.   
C                                                                        986.   
C----------------------------------------------------------------------- 987.   
C                IF(KGASSR.GT.0)  REDEFINE ULGAS FOR SOLAR FULGAS VALUES 988.   
      ENTRY SOLGAS                                                       989.   
C----------------------------------------------------------------------- 990.   
C                                                                        991.   
C                                                                        992.   
      DO 400 L=1,NL                                                      993.   
      ULGAS(L,1)=U0GAS(L,1)*FULGAS(1+9)                                  994.   
      ULGAS(L,3)=U0GAS(L,3)*FULGAS(3+9)                                  995.   
  400 ULGAS(L,5)=U0GAS(L,5)*FULGAS(5+9)                                  996.   
C                                                                        997.   
      PART=(PLB(1)-PTOP)/(PS0-PTOP)                                      998.   
      DO 410 L=1,NL                                                      999.   
      IF(L.EQ.NLMOD) PART=1.                                            1000.   
      ULGAS(L,2)=U0GAS(L,2)*FULGAS(2+9)*PART                            1001.   
      ULGAS(L,4)=U0GAS(L,4)*FULGAS(4+9)*PART                            1002.   
      ULGAS(L,6)=U0GAS(L,6)*FULGAS(6+9)*PART                            1003.   
      ULGAS(L,7)=U0GAS(L,7)*FULGAS(7+9)*PART                            1004.   
      ULGAS(L,8)=U0GAS(L,8)*FULGAS(8+9)*PART                            1005.   
  410 ULGAS(L,9)=U0GAS(L,9)*FULGAS(9+9)*PART                            1006.   
C                                                                       1007.   
C                                                                       1008.   
      RETURN                                                            1009.   
      END                                                               1010.   
      SUBROUTINE SETAER                                                 1011.   

#include "chem_para"
#include "chem_com"
#include "B83XX.COM"                                                    1012.   

C                                                                       1073.   
      EQUIVALENCE (FEMTRA(1),ECLTRA)                                    1074.   
      EQUIVALENCE (ISPARE(2),NEWAQA)                                    1074.1  
      EQUIVALENCE (ISPARE(3),NEWCQA)                                    1074.2  
C                                                                       1075.   
      DIMENSION SRAX(40,6,5),SRAS(40,6,5),SRAC(40,6,5)                  1076.   
C                                                                       1077.   
C-----------------------------------------------------------------------1078.   
C     THERMAL: SET (5) AEROSOL TYPE COMPOSITIONS & VERTICAL DISTRIBUTION1079.   
C-----------------------------------------------------------------------1080.   
C                                                                       1081.   
      DO 100 J=1,NGOLDH                                                 1082.   
      DO 100 K=1,NKTR                                                   1083.   
      DO 100 L=1,NL                                                     1084.   
  100 TRAX(L,K,J)=0.                                                    1085.   
C                                                                       1086.   
      DO 103 I=1,NAERO                                                  1087.   
      DO 103 J=1,NGOLDH                                                 1088.   
      IF(AGOLDH(I,J).LT.1.E-06) GO TO 103                               1089.   
      C=CGOLDH(I,J)                                                     1090.   
      BC=EXP(-BGOLDH(I,J)/C)                                            1091.   
      ABC=AGOLDH(I,J)*(1.0+BC)                                          1092.   
C                                                                       1093.   
      DO 102 L=1,NL                                                     1094.   
C AMIN
      ABCD=ABC/(1.0+BC*EXP(DMIN1(HLB(L  )/C,80.)))                      1095.
     +    -ABC/(1.0+BC*EXP(DMIN1(HLB(L+1)/C,80.)))                      1096.
C AMIN
      DO 101 K=1,NKTR                                                   1097.   
      TRANEW=TRACOS(K,I)                                                1097.5  
      IF(NEWAQA.GT.0) TRANEW=1.0                                        1097.6  
  101 TRAX(L,K,J)=TRAX(L,K,J)+ABCD*(TRAQEX(K,I)-TRANEW*TRAQSC(K,I))     1098.   
  102 CONTINUE                                                          1099.   
  103 CONTINUE                                                          1100.   
C                                                                       1101.   
      DO 104 J=1,2                                                      1102.   
      DO 104 K=1,NKTR                                                   1103.   
      TRCNEW=TRCCOS(K,J)                                                1103.5  
      IF(NEWCQA.GT.0) TRCNEW=1.0                                        1103.6  
  104 TRCX(K,J)=TRCQEX(K,J)-TRCNEW*TRCQSC(K,J)                          1104.   
C                                                                       1105.   
C-----------------------------------------------------------------------1106.   
C     SOLAR:   SET (5) AEROSOL TYPE COMPOSITIONS & VERTICAL DISTRIBUTION1107.   
C-----------------------------------------------------------------------1108.   
C                                                                       1109.   
      DO 110 J=1,NGOLDH                                                 1110.   
      DO 110 K=1,NKSR                                                   1111.   
      DO 110 L=1,NL                                                     1112.   
      SRAX(L,K,J)=1.E-30                                                1113.   
      SRAS(L,K,J)=1.E-31                                                1114.   
  110 SRAC(L,K,J)=0.                                                    1115.   
C                                                                       1116.   
      DO 113 I=1,NAERO                                                  1117.   
      DO 113 J=1,NGOLDH                                                 1118.   
      IF(AGOLDH(I,J).LT.1.E-06) GO TO 113                               1119.   
      C=CGOLDH(I,J)                                                     1120.   
      BC=EXP(-BGOLDH(I,J)/C)                                            1121.   
      ABC=AGOLDH(I,J)*(1.0+BC)                                          1122.   
C                                                                       1123.   
      DO 112 L=1,NL                                                     1124.   
C AMIN
      ABCD=ABC/(1.0+BC*EXP(DMIN1(HLB(L  )/C,80.)))                      1095.
     +    -ABC/(1.0+BC*EXP(DMIN1(HLB(L+1)/C,80.)))                      1096.
C AMIN
      DO 111 K=1,NKSR                                                   1127.   
      SRAX(L,K,J)=SRAX(L,K,J)+ABCD*SRAQEX(K,I)                          1128.   
      SRAS(L,K,J)=SRAS(L,K,J)+ABCD*SRAQSC(K,I)                          1129.   
  111 SRAC(L,K,J)=SRAC(L,K,J)+ABCD*SRACOS(K,I)*SRAQSC(K,I)              1130.   
  112 CONTINUE                                                          1131.   
  113 CONTINUE                                                          1132.   
C                                                                       1133.   
      DO 114 J=1,NGOLDH                                                 1134.   
      DO 114 K=1,NKSR                                                   1135.   
      DO 114 L=1,NL                                                     1136.   
  114 SRAC(L,K,J)=SRAC(L,K,J)/SRAS(L,K,J)                               1137.   
C                                                                       1138.   
C-----------------                                                      1139.   
      ENTRY GETAER                                                      1140.   
C-----------------                                                      1141.   
C                                                                       1142.   
C-----------------------------------------------------------------------1143.   
C                            GET CLOUD & AEROSOL AMOUNTS & DISTRIBUTIONS1144.   
C-----------------------------------------------------------------------1145.   
      LBOTCL=0                                                          1146.   
      LTOPCL=0                                                          1147.   
      DO 203 L=1,NL                                                     1148.   
      KCLD=1                                                            1149.   
      IF(TLM(L).LT.TKCICE) KCLD=2                                       1150.   
      IF(CLDTAU(NLP-L).GT.0.1) LTOPCL=NLP-L                             1151.   
C$    IF(CLDTAU(NLP-L).GT.0.1) LBOTCL=NLP-L   *******************CORRECT1152.   
      IF(CLDTAU(    L).GT.0.1) LBOTCL=L                                 1153.   
C$    IF(CLDTAU(    L).GT.0.1) LTOPCL=L   ***********************CORRECT1154.   
C                                                (THERMAL)              1155.   
C                                                ---------              1156.   
      DO 202 K=1,NKTR                                                   1157.   
      SUMEXT=1.E-30                                                     1158.   
      DO 201 J=1,NGOLDH                                                 1159.   
  201 SUMEXT=SUMEXT+FGOLDH(J)*TRAX(L,K,J)                               1160.   
      TRAEXT(L,K)=SUMEXT+CLDTAU(L)*TRCX(K,KCLD)*FCLDTR                  1161.   
  202 TAUN(L+(K-1)*NL)=TAUN(L+(K-1)*NL)+TRAEXT(L,K)                     1162.   
  203 CONTINUE                                                          1163.   
C                                                                       1164.   
C-----------------------------------------------------------------------1165.   
C                         CLOUD ALBEDO & SURFACE LAYER FOG SPECIFICATION1166.   
C-----------------------------------------------------------------------1167.   
C                                                                       1168.   
      DO 204 K=1,NKTR                                                   1169.   
  204 FTAUSL(K)=FOGTSL*TRCX(K,1)*FCLDTR                                 1170.   
      IF(LTOPCL.GT.0) GO TO 206                                         1171.   
      DO 205 K=1,NKTR                                                   1172.   
  205 TRCALB(K)=0.                                                      1173.   
      GO TO 210                                                         1174.   
  206 KCLD=1                                                            1175.   
      IF(TLM(LTOPCL).LT.TKCICE) KCLD=2                                  1176.   
      DO 207 K=1,NKTR                                                   1177.   
  207 TRCALB(K)=(1.0-EXP(-CLDTAU(LTOPCL)*TRCX(K,KCLD)))*CLDALB(K,KCLD)  1178.   
     +         *ECLTRA*FCLDTR                                           1179.   
  210 CONTINUE                                                          1180.   
C                                                (SOLAR)                1181.   
C                                                -------                1182.   
      KSR=9*KAERSR                                                      1183.   
      DO 9212 K=1,NKSR                                                  1184.   
      DO 212 L=1,NL                                                     1185.   
      EXTSUM=1.E-30                                                     1186.   
      SCTSUM=1.E-31                                                     1187.   
      COSSUM=0.                                                         1188.   
      DO 211 J=1,NGOLDH                                                 1189.   
      EXTSUM=EXTSUM+FGOLDH(J+KSR)*SRAX(L,K,J)                           1190.   
      SCTSUM=SCTSUM+FGOLDH(J+KSR)*SRAS(L,K,J)                           1191.   
  211 COSSUM=COSSUM+FGOLDH(J+KSR)*SRAS(L,K,J)*SRAC(L,K,J)               1192.   

#if ( (defined PREDICTED_BC || defined PREDICTED_AEROSOL) && !defined FIXED_FOR )
!
! --- Chemistry model patch, 092901
!
! ===  Chien Wang
! ===   (1) add to type 3 aerosol with
! ===    chemistry model predicted S(VI);
! ===   (2) add type 11 aerosol with
! ===    chemistry model predicted bcarbon
! ===
       if ( L .le. nlev1 ) then
         !
         ! === add as global aerosol
         !     Note: if needed the AGOLDH for prescribed
         !       tropospheric S(VI), SLFT1 & SLFT2, can be
         !       set to zero in later part of the code
         !
         !           FAERSOL/svi_intensity is added for using
         !       FAERSOL to switch between diagnostic/prognostic loops
         !       while normalize it to 1 in prognostic loop
         !       FBC added for black carbon  7/22/04
         !
                dsviod = 0.0
                dbcod  = 0.0

#if ( defined PREDICTED_AEROSOL  && !defined FIXED_FOR )
          dsviod = max(0.0,
     &            (sviod(1,jlat,L) - sviod(1,jlat,L+1))
     &             *FAERSOL )
#endif

#if ( defined PREDICTED_BC  && !defined FIXED_FOR )
          dbcod  = max(0.0,
     &            (bcod(1,jlat,L) - bcod(1,jlat,L+1))
     &             *FBC )
#endif

          EXTSUM = EXTSUM
     &           + dsviod*SRAQEX(K,3)
     &           + dbcod*SRAQEX(K,11)
          SCTSUM = SCTSUM
     &           + dsviod*SRAQSC(K,3)
     &           + dbcod*SRAQSC(K,11)
          COSSUM = COSSUM
     &           + dsviod*SRAQSC(K,3)*SRACOS(K,3)
     &           + dbcod*SRAQSC(K,11)*SRACOS(K,11)

            if(jlat.eq.-22.or.jlat.eq.-33)then
             if(L.eq.1.and.k.eq.1)then
             print *,'From r95 jlat=',jlat,' L=',L
c            print *,' LATHEM=',LATHEM, ' JNORTH=',JNORTH
c            print *,'FAERSOL=',FAERSOL,' FBC=',FBC
             print *,sviod(1,jlat,L),sviod(1,jlat,L+1)
c            print *,dsviod,SRAQEX(K,3)
             print *,bcod(1,jlat,L),bcod(1,jlat,L+1)
c            print *,dbcod,SRAQEX(K,11)
c            print *,SRAQSC(K,11),SRACOS(K,11)
             endif
            endif
       end if
#endif

      EXTAER(L,K)=EXTSUM                                                1193.   
      SCTAER(L,K)=SCTSUM                                                1194.   
      COSAER(L,K)=COSSUM/SCTSUM                                         1195.   

212   continue
9212  continue
c
c ======================================================

      IF(NTRACE.GT.0) GO TO 300                                         1196.   
C                                                                       1197.   
C-----------                                                            1198.   
      RETURN                                                            1199.   
C-----------                                                            1200.   
C                                                                       1201.   
  300 CONTINUE                                                          1202.   
C-----------------------------------------------------------------------1203.   
C                 ADD TRACER AEROSOL THERMAL & SOLAR CONTRIBUTIONS      1204.   
C-----------------------------------------------------------------------1205.   
      DO 303 JJ=1,NTRACE                                                1206.   
      J=NGOLDH+JJ                                                       1207.   
      I=ITR(JJ)                                                         1208.   
C                                                (THERMAL)              1209.   
C                                                ---------              1210.   
      DO 302 K=1,NKTR                                                   1211.   
C$    SUMEXT=FGOLDH(J+KSR)*(TRAQEX(K,I)-TRACOS(K,I)*TRAQSC(K,I))        1212.   
      SUMEXT=FGOLDH(J+KSR)*(TRAQEX(K,I)-TRAQSC(K,I))                    1212.11 
      DO 301 L=1,NL                                                     1213.   
  301 TAUN(L+(K-1)*NL)=TAUN(L+(K-1)*NL)+SUMEXT*TRACER(L,JJ)             1214.   
  302 CONTINUE                                                          1215.   
  303 CONTINUE                                                          1216.   
C                                                                       1217.   
C                                                (SOLAR)                1218.   
C                                                -------                1219.   
      DO 305 K=1,NKSR                                                   1220.   
      DO 305 L=1,NL                                                     1221.   
      EXTSUM=EXTAER(L,K)                                                1222.   
      SCTSUM=SCTAER(L,K)                                                1223.   
      COSSUM=COSAER(L,K)*SCTAER(L,K)                                    1224.   
      DO 304 JJ=1,NTRACE                                                1225.   
      J=NGOLDH+JJ                                                       1226.   
      I=ITR(JJ)                                                         1227.   
      EXTSUM=EXTSUM+FGOLDH(J+KSR)*TRACER(L,JJ)*SRAQEX(K,I)              1228.   
      SCTSUM=SCTSUM+FGOLDH(J+KSR)*TRACER(L,JJ)*SRAQSC(K,I)              1229.   
  304 COSSUM=COSSUM+FGOLDH(J+KSR)*TRACER(L,JJ)*SRAQSC(K,I)*SRACOS(K,I)  1230.   
      EXTAER(L,K)=EXTSUM                                                1231.   
      SCTAER(L,K)=SCTSUM                                                1232.   
  305 COSAER(L,K)=COSSUM/SCTSUM                                         1233.   
      RETURN                                                            1234.   
      END                                                               1235.   
      SUBROUTINE TAUGAS                                                 1236.   

#include "B83XX.COM"

C     TAUGAS INPUT REQUIRES:  NL,TLM,ULGAS,TRACEG,PL,DPL,TAUTBL,MLGAS   1295.11 
C     TAUGAS OUTPUT DATA IS:  TAUN                                      1295.12 
C                                                                       1296.   
      DIMENSION IGASX(11),KGX(11),NUX(11),IGUX(11),NGX(3),IG1X(3)       1297.   
      DIMENSION ULOX(165),DUX(165),PX(15),H2OCON(25)                    1298.   
C                                                                       1299.   
      DATA NTX/8/, TLOX/181./,DTX/23./                                  1300.   
      DATA NPX/15/, PX/1000., 975., 910., 800., 645.,                   1301.   
     *                  480., 330., 205., 110.,  40.,                   1302.   
     *                   7.5,  3.5,  1.0,  0.1, .001/                   1303.   
C                                                                       1304.   
      DATA NGUX/652/, NPUX/15/                                          1305.   
      DATA NGX/10,10,04/, IG1X/2,12,22/                                 1306.   
      DATA                                                              1307.   
     *       IGASX/  1,  2,  3,  1,  1,  2,  2,  3,  6,  6,  7/,        1308.   
     *         KGX/  1,  2,  3,  2,  3,  1,  3,  2,  1,  2,  1/,        1309.   
     *         NUX/ 25,  9,  9,  9,  9,  5,  5,  5,  1,  1,  1/,        1310.   
     *        IGUX/  0,250,340,376,466,502,552,572,622,632,642/         1311.   
C                                                                       1312.   
C                                                                       1313.   
      DATA ULOX/    .25E+2,.25E+2,.50E+2,.50E+2,.25E+2,.50E+1,.10E+1,   1314.   
     *.25E+0,.10E+0,.50E-1,.10E-1,.10E-1,.10E-3,.10E-5,.10E-5,          1315.   
     *.50E+1,.50E+1,.80E+1,.10E+2,.20E+2,.20E+2,.10E+1,.10E+2,.80E+1,   1316.   
     *.10E+1,.10E+1,.25E+0,.50E-1,.50E-2,.50E-3,       .10E-3,.10E-3,   1317.   
     *.40E-3,.60E-3,.10E-2,.24E-2,.48E-2,.48E-2,.64E-2,.64E-2,.64E-2,   1318.   
     *.40E-2,.10E-4,.80E-7,.40E-7,       .25E+2,.25E+2,.50E+2,.50E+2,   1319.   
     *.25E+2,.50E+1,.10E+1,.25E+0,.10E+0,.50E-1,.10E-1,.10E-1,.10E-3,   1320.   
     *.10E-5,.10E-5,       .25E+2,.25E+2,.50E+2,.50E+2,.25E+2,.50E+1,   1321.   
     *.10E+1,.25E+0,.10E+0,.50E-1,.10E-1,.10E-1,.10E-3,.10E-5,.10E-5,   1322.   
     *       .50E+1,.50E+1,.80E+1,.10E+2,.20E+2,.20E+2,.10E+2,.10E+2,   1323.   
     *.80E+1,.10E+1,.10E+1,.25E+0,.50E-1,.50E-2,.50E-3,       .50E+1,   1324.   
     *.50E+1,.80E+1,.10E+2,.20E+2,.20E+2,.10E+2,.10E+2,.80E+1,.10E+1,   1325.   
     *.10E+1,.25E+0,.50E-1,.50E-2,.50E-3,       .10E-3,.10E-3,.40E-3,   1326.   
     *.60E-3,.10E-2,.24E-2,.48E-2,.48E-2,.64E-2,.64E-2,.64E-2,.40E-2,   1327.   
     *.10E-4,.80E-7,.40E-7,       .11E-1,.11E-1,.18E-1,.31E-1,.37E-1,   1328.   
     *.35E-1,.31E-1,.24E-1,.18E-1,.13E-1,.11E-2,.66E-3,.44E-3,.44E-4,   1329.   
     *.44E-6,       .11E-1,.11E-1,.18E-1,.31E-1,.37E-1,.35E-1,.31E-1,   1330.   
     *.24E-1,.18E-1,.13E-1,.11E-2,.66E-3,.44E-3,.44E-4,.44E-6,          1331.   
     *.64E-1,.64E-1,.10E+0,.18E+0,.22E+0,.20E+0,.18E+0,.14E+0,.10E+0,   1332.   
     *.77E-1,.64E-2,.38E-2,.26E-2,.26E-3,.26E-5/                        1333.   
C                                                                       1334.   
      DATA  DUX/    .75E+2,.75E+2,.10E+3,.10E+3,.75E+2,.50E+2,.10E+2,   1335.   
     *.20E+1,.20E+0,.10E+0,.50E-1,.10E-1,.40E-2,.40E-3,.40E-4,          1336.   
     *.50E+1,.50E+1,.80E+1,.10E+2,.10E+2,.10E+2,.10E+2,.10E+2,.80E+1,   1337.   
     *.50E+1,.35E+1,.25E+0,.25E+0,.10E+0,.10E-1,       .30E-3,.30E-3,   1338.   
     *.50E-3,.80E-3,.10E-2,.16E-2,.64E-2,.16E-2,.25E-1,.25E-1,.25E-1,   1339.   
     *.45E-2,.25E-2,.10E-2,.25E-4,       .24E+3,.24E+3,.30E+3,.30E+3,   1340.   
     *.24E+3,.15E+3,.30E+2,.60E+1,.60E+0,.30E+0,.15E+0,.30E-1,.12E-1,   1341.   
     *.12E-2,.12E-3,       .24E+3,.24E+3,.30E+3,.30E+3,.24E+3,.15E+3,   1342.   
     *.30E+2,.60E+1,.60E+0,.30E+0,.15E+0,.30E-1,.12E-1,.12E-2,.12E-3,   1343.   
     *       .10E+2,.10E+2,.16E+2,.20E+2,.20E+2,.20E+2,.20E+2,.20E+2,   1344.   
     *.16E+2,.10E+2,.70E+1,.50E+0,.50E+0,.20E+0,.20E-1,       .10E+2,   1345.   
     *.10E+2,.16E+2,.20E+2,.20E+2,.20E+2,.20E+2,.20E+2,.16E+2,.10E+2,   1346.   
     *.70E+1,.50E+0,.50E+0,.20E+0,.20E-1,       .60E-3,.60E-3,.10E-2,   1347.   
     *.16E-2,.20E-2,.32E-2,.13E-1,.32E-1,.50E-1,.50E-1,.50E-1,.90E-2,   1348.   
     *.50E-2,.20E-2,.50E-4,       45*0./                                1349.   
C                                                                       1350.   
      DATA H2OCON/       .767116, .322401,  .572299,.58537, .48869,     1351.   
     *   .43539, .44322, .64072,  .89293,  1.12733,1.65550, .865210,    1352.   
     *  1.38403,1.80159,1.99196, 2.03403,  2.20561,2.42859,2.56883,     1353.   
     *  2.67157,2.71888, .45534,  .44735,   .44534, .44365/             1354.   
C                                                                       1355.   
C--------------------------------------------------------------------   1356.   
C        ABSORPTION (TAU) INTERPOLATION FOR GAS AMOUNTS IN ULGAS(N,K)   1357.   
C--------------------------------------------------------------------   1358.   
C                                                                       1359.   
      IPX=2                                                             1360.   
      DO 100 IP=1,NL                                                    1361.   
C                                                                       1362.   
   20 WPB = (PL(IP)-PX(IPX))/(PX(IPX-1)-PX(IPX))                        1363.   
      IF(WPB.GE.0. .OR. IPX.GE.NPX) GO TO 30                            1364.   
      IPX = IPX+1                                                       1365.   
      GO TO 20                                                          1366.   
C                                                                       1367.   
   30 WTB = (TLM(IP)-TLOX)/DTX                                          1368.   
      ITX = MIN0(MAX0(INT(WTB),0),NTX-2)                                1369.   
      WTB = WTB-FLOAT(ITX)                                              1370.   
C                                                                       1371.   
      WBB = WPB*WTB                                                     1372.   
      WBA = WPB-WBB                                                     1373.   
      WAB = WTB-WBB                                                     1374.   
      WAA = 1.-(WBB+WBA+WAB)                                            1375.   
C                                                                       1376.   
      IAA = NGUX*(ITX+NTX*(IPX-1))                                      1377.   
      IBA = IAA-NGUX*NTX                                                1378.   
C                                                                       1379.   
      DO 90 IGAS=1,11                                                   1380.   
      IF(MLGAS(IGAS).LT.1) GO TO 90                                     1381.   
C                                                                       1382.   
      UGAS = ULGAS(IP,IGASX(IGAS))                                      1383.   
      IF(UGAS.LT.1.E-10) GO TO 90                                       1384.   
C                                                                       1385.   
      IU = IPX + NPUX*(IGAS-1)                                          1386.   
      NU = NUX(IGAS)                                                    1387.   
      IF(NU.GT.1) GO TO 40                                              1388.   
      XUA = 0.                                                          1389.   
      XUB = 0.                                                          1390.   
      GO TO 50                                                          1391.   
   40 XUA = (UGAS-ULOX(IU))/DUX(IU)                                     1392.   
      XUB = (UGAS-ULOX(IU-1))/DUX(IU-1)                                 1393.   
   50 IUA = INT(XUA)                                                    1394.   
      IUB = INT(XUB)                                                    1395.   
C                                                                       1396.   
      QAA = 1.                                                          1397.   
      QAB = 1.                                                          1398.   
      IF(XUA.GT.0. .AND. IUA.LT.NU-1) GO TO 60                          1399.   
c     XUA = DMIN1(DMAX1(XUA,0.),FLOAT(NU-1))                            1400.   
      XUA = DMIN1(DMAX1(XUA,0.),dble(NU-1))                                1400.   
      IUA = MIN0(INT(XUA),NU-2)                                         1401.   
      QAA = UGAS/(ULOX(IU)+DUX(IU)*FLOAT(IUA))                          1402.   
      QAB = UGAS/(ULOX(IU)+DUX(IU)*FLOAT(IUA+1))                        1403.   
C                                                                       1404.   
   60 QBA = 1.                                                          1405.   
      QBB = 1.                                                          1406.   
      IF(XUB.GT.0. .AND. IUB.LT.NU-1) GO TO 70                          1407.   
c     XUB = DMIN1(DMAX1(XUB,0.),FLOAT(NU-1))                            1408.   
      XUB = DMIN1(DMAX1(XUB,0.),dble(NU-1))                              1408.   
      IUB = MIN0(INT(XUB),NU-2)                                         1409.   
      QBA = UGAS/(ULOX(IU-1)+DUX(IU-1)*FLOAT(IUB))                      1410.   
      QBB = UGAS/(ULOX(IU-1)+DUX(IU-1)*FLOAT(IUB+1))                    1411.   
C                                                                       1412.   
   70 UAB = XUA-FLOAT(IUA)                                              1413.   
      UBB = XUB-FLOAT(IUB)                                              1414.   
      UAA = 1.-UAB                                                      1415.   
      UBA = 1.-UBB                                                      1416.   
C                                                                       1417.   
C                                                                       1418.   
      WAAA = WAA*UAA*QAA                                                1419.   
      WAAB = WAA*UAB*QAB                                                1420.   
      WABA = WAB*UAA*QAA                                                1421.   
      WABB = WAB*UAB*QAB                                                1422.   
      WBAA = WBA*UBA*QBA                                                1423.   
      WBAB = WBA*UBB*QBB                                                1424.   
      WBBA = WBB*UBA*QBA                                                1425.   
      WBBB = WBB*UBB*QBB                                                1426.   
C                                                                       1427.   
      NG = NGX(KGX(IGAS))                                               1428.   
      IAAA = IAA+IGUX(IGAS) + NG*IUA                                    1429.   
      IAAB = IAAA+NG                                                    1430.   
      IABA = IAAA+NGUX                                                  1431.   
      IABB = IABA+NG                                                    1432.   
      IBAA = IBA+IGUX(IGAS) + NG*IUB                                    1433.   
      IBAB = IBAA+NG                                                    1434.   
      IBBA = IBAA+NGUX                                                  1435.   
      IBBB = IBBA+NG                                                    1436.   
C                                                                       1437.   
C                                                                       1438.   
      IPG = IP+NL*(IG1X(KGX(IGAS))-1)                                   1439.   
      DO 80 IG=1,NG                                                     1440.   
      TAUN(IPG) = TAUN(IPG)                                             1441.   
     *     + WAAA*TAUTBL(IAAA+IG)                                       1442.   
     *     + WAAB*TAUTBL(IAAB+IG)                                       1443.   
     *     + WABA*TAUTBL(IABA+IG)                                       1444.   
     *     + WABB*TAUTBL(IABB+IG)                                       1445.   
     *     + WBAA*TAUTBL(IBAA+IG)                                       1446.   
     *     + WBAB*TAUTBL(IBAB+IG)                                       1447.   
     *     + WBBA*TAUTBL(IBBA+IG)                                       1448.   
     *     + WBBB*TAUTBL(IBBB+IG)                                       1449.   
   80 IPG = IPG+NL                                                      1450.   
   90 CONTINUE                                                          1451.   
  100 CONTINUE                                                          1452.   
C                                                                       1453.   
      IF(MLGAS(12).LT.1) GO TO 110                                      1454.   
C-------------------------------------------------------------------    1455.   
C                                    PICK UP CCL3F1 (F11) ABSORPTION    1456.   
C-------------------------------------------------------------------    1457.   
C                                                                       1458.   
      DO 102 K=1,25                                                     1459.   
      XKPCMA=TRACEG(K,1)                                                1460.   
      IF(XKPCMA.LT.1.E-10) GO TO 102                                    1461.   
      DO 101 N=1,NL                                                     1462.   
      NK=N+(K-1)*NL                                                     1463.   
  101 TAUN(NK)=TAUN(NK)+ULGAS(N,8)*XKPCMA                               1464.   
  102 CONTINUE                                                          1465.   
C                                                                       1466.   
  110 IF(MLGAS(13).LT.1) GO TO 120                                      1467.   
C-------------------------------------------------------------------    1468.   
C                                    PICK UP CCL2F2 (F12) ABSORPTION    1469.   
C-------------------------------------------------------------------    1470.   
C                                                                       1471.   
      DO 112 K=1,25                                                     1472.   
      XKPCMA=TRACEG(K,2)                                                1473.   
      IF(XKPCMA.LT.1.E-10) GO TO 112                                    1474.   
      DO 111 N=1,NL                                                     1475.   
      NK=N+(K-1)*NL                                                     1476.   
  111 TAUN(NK)=TAUN(NK)+ULGAS(N,9)*XKPCMA                               1477.   
  112 CONTINUE                                                          1478.   
C                                                                       1479.   
  120 IF(MLGAS(14).LT.1) GO TO 130                                      1480.   
C-------------------------------------------------------------------    1481.   
C                              PICK UP WINDOW H2O GASEOUS ABSORPTION    1482.   
C-------------------------------------------------------------------    1483.   
C                                                                       1484.   
      DO 121 N=1,NL                                                     1485.   
      TAUN(N) = TAUN(N)                                                 1486.   
  121 CONTINUE                                                          1487.   
  130 CONTINUE                                                          1488.   
C-------------------------------------------------------------------    1489.   
C                                   PICK UP H2O CONTINUUM ABSORPTION    1490.   
C-------------------------------------------------------------------    1491.   
C                                                                       1492.   
      IF(MLGAS(15).LT.1) GO TO 140                                      1493.   
      DO 131 N=1,NL                                                     1494.   
      TAUN(N) = TAUN(N) + 2.21866E-11*                                  1495.   
     *          PL(N)*ULGAS(N,1)*EXP(1800./TLM(N))*                     1496.   
     *          (ULGAS(N,1)/DPL(N)+.808563)                             1497.   
  131 CONTINUE                                                          1498.   
C                                                                       1499.   
C$ ********************************REMOVE FOLLOWING STATEMENT TO CORRECT1500.   
      IF(NL.GT.0) RETURN                                                1501.   
      DO 133 N=1,NL                                                     1502.   
      PH2O=12.38E-4*ULGAS(N,1)*PL(N)/DPL(N)                             1503.   
      TH2O=EXP(1800./TLM(N)-6.081081)                                   1504.   
      COEC=PH2O*TH2O+.0015*(PL(N)-PH2O)                                 1505.   
      DO 132 K=2,25                                                     1506.   
      COEF=H2OCON(K)*1.E-5                                              1507.   
      NK=N+(K-1)*NL                                                     1508.   
  132 TAUN(NK)=TAUN(NK)+ULGAS(N,1)*COEC*COEF                            1509.   
  133 CONTINUE                                                          1510.   
  140 CONTINUE                                                          1511.   
C                                                                       1512.   
      RETURN                                                            1513.   
      END                                                               1514.   
      SUBROUTINE THERML                                                 1515.   

#include "B83XX.COM"
#if ( defined CLM )
#include "CLM.h"
#endif
  
      DATA R6,R24/.1666667,4.166667E-02/                                1577.   
      DATA A,B,C/0.3825,0.5742,0.0433/                                  1578.   
C                                                                       1579.   
C-----------------------------------------------------------------------1580.   
C                                   LAYER EDGE TEMPERATURE INTERPOLATION1581.   
C-----------------------------------------------------------------------1582.   
      IF(TLGRAD.LT.0.) GO TO 103                                        1583.   
      TA=TLM(1)                                                         1584.   
      TB=TLM(2)                                                         1585.   
      P1=PLB(1)                                                         1586.   
      P2=PLB(2)                                                         1587.   
      P3=PLB(3)                                                         1588.   
      DT1CPT=0.5*TA*(EXPBYK(PLB(1))-EXPBYK(PLB(2)))/EXPBYK(PL(1))       1589.   
      DTHALF=(TA-TB)*(P1-P2)/(P1-P3)                                    1590.   
      IF(DTHALF.GT.DT1CPT) DTHALF=DT1CPT                                1591.   
      TLB(1)=TA+DTHALF*TLGRAD                                           1592.   
      TLT(1)=TA-DTHALF*TLGRAD                                           1593.   
      DO 101 L=3,NL                                                     1594.   
      TC=TLM(L)                                                         1595.   
      P4=PLB(L+1)                                                       1596.   
      DTHALF=0.5*((TA-TB)/(P1-P3)+(TB-TC)/(P2-P4))*(P2-P3)*TLGRAD       1597.   
      TLB(L-1)=TB+DTHALF                                                1598.   
      TLT(L-1)=TB-DTHALF                                                1599.   
      TA=TB                                                             1600.   
      TB=TC                                                             1601.   
      P1=P2                                                             1602.   
      P2=P3                                                             1603.   
  101 P3=P4                                                             1604.   
      DTHALF=(TA-TB)*(P2-P3)/(P1-P3)*TLGRAD                             1605.   
      TLB(NL)=TC+DTHALF                                                 1606.   
      TLT(NL)=TC-DTHALF                                                 1607.   
      L=NLP                                                             1608.   
      DO 102 N=1,NL                                                     1609.   
      L=L-1                                                             1610.   
      IF(PLB(L).GT.PTLISO) GO TO 103                                    1611.   
      TLT(L)=TLM(L)                                                     1612.   
  102 TLB(L)=TLM(L)                                                     1613.   
  103 CONTINUE                                                          1614.   
C-----------------------------------------------------------------------1615.   
C                   WEIGHT ASSIGNMENTS FOR PLANCK FUNCTION INTERPOLATION1616.   
C-----------------------------------------------------------------------1617.   
      DO 104 L=1,NL                                                     1618.   
      ITL=TLB(L)                                                        1619.   
      WTLB(L)=TLB(L)-ITL                                                1620.   
      ITLB(L)=ITL-IT0                                                   1621.   
      ITL=TLT(L)                                                        1622.   
      WTLT(L)=TLT(L)-ITL                                                1623.   
  104 ITLT(L)=ITL-IT0                                                   1624.   
      ITS=TSL                                                           1625.   
      WTS=TSL-ITS                                                       1626.   
      ITS=ITS-IT0                                                       1627.   
C                                                                       1628.   
C     ------------------------------------------------------------------1629.   
C                                         WINDOW REGION FLUX COMPUTATION1630.   
C     ------------------------------------------------------------------1631.   
C                                                          DOWNWARD FLUX1632.   
C     ------------------------------------------------------------------1633.   
      K=1                                                               1634.   
      BG=BGFEMT(K)                                                      1635.   
c     print *,'1635 K=',k,' PEARTH=',PEARTH
c     print *,'BG=',BG
      WTS1=1.-WTS                                                       1636.   
      TRSLTS=0.                                                         1637.   
      TRSLTG=0.                                                         1638.   
      TRSLWV=0.                                                         1639.   
      TRSLBS=0.                                                         1640.   
      DNA=0.                                                            1641.   
      DNB=0.                                                            1642.   
      DNC=0.                                                            1643.   
      NLK0=0                                                            1644.   
      NLK=NL                                                            1645.   
      TRDFLB(NLP)=0.                                                    1646.   
  100 TAUA=TAUN(NLK)                                                    1647.   
      IF(TAUA.GT.1.E-05) GO TO 120                                      1648.   
      TRDFLB(NLK)=0.                                                    1649.   
      NLK=NLK-1                                                         1650.   
      IF(NLK.GT.NLK0) GO TO 100                                         1651.   
  110 NLK=NLK+1                                                         1652.   
      TRUFLB(NLK)=BG                                                    1653.   
      IF(NLK.LT.NLP) GO TO 110                                          1654.   
      TRUFG=BG                                                          1655.   
      TRDFG=0.                                                          1656.   
      TRUFGW=BG                                                         1657.   
      TRUFGW=0.                                                         1658.   
      TRUFTW=TRUFLB(NLP)                                                1659.   
      GO TO 200                                                         1660.   
  120 N=NLK                                                             1661.   
  130 ITL=ITLT(N)                                                       1662.   
      BTOP=PLANCK(ITL)-(PLANCK(ITL)-PLANCK(ITL+1))*WTLT(N)              1663.   
      ITL=ITLB(N)                                                       1664.   
      BBOT=PLANCK(ITL)-(PLANCK(ITL)-PLANCK(ITL+1))*WTLB(N)              1665.   
      TAUA=TAUN(N)                                                      1666.   
      TAUB=TAUA+TAUA                                                    1667.   
      TAUC=10.*TAUA                                                     1668.   
      IF(TAUA.GT.1.E-01) GO TO 140                                      1669.   
      IF(TAUA.LT.1.E-03) GO TO 135                                      1670.   
      TAU2=TAUA*TAUA                                                    1671.   
      BDIF=BBOT-BTOP                                                    1672.   
      BBTA=BDIF/TAUA                                                    1673.   
      BBTB=BDIF/TAUB                                                    1674.   
      BBTC=BDIF/TAUC                                                    1675.   
      TRAN=1.-TAUA+(0.5-R6*TAUA+R24*TAU2)*TAU2                          1676.   
      GO TO 145                                                         1677.   
  135 BDIF=.5*(BTOP+BBOT)                                               1678.   
      TRA(N)=1.-TAUA                                                    1679.   
      ENA(N)=BDIF*TAUA                                                  1680.   
      DNA=DNA*TRA(N)+ENA(N)                                             1681.   
      TRB(N)=1.-TAUB                                                    1682.   
      ENB(N)=BDIF*TAUB                                                  1683.   
      DNB=DNB*TRB(N)+ENB(N)                                             1684.   
      TRC(N)=1.-TAUC                                                    1685.   
      ENC(N)=BDIF*TAUC                                                  1686.   
      DNC=DNC*TRC(N)+ENC(N)                                             1687.   
      GO TO 160                                                         1688.   
  140 BDIF=BBOT-BTOP                                                    1689.   
      BBTA=BDIF/TAUA                                                    1690.   
      BBTB=BDIF/TAUB                                                    1691.   
      BBTC=BDIF/TAUC                                                    1692.   
      IF(TAUA.GT.7.) GO TO 150                                          1693.   
      TRAN=EXP(-TAUA)                                                   1694.   
  145 TRA(N)=TRAN                                                       1695.   
      ENA(N)=BTOP+BBTA-(BBOT+BBTA)*TRAN                                 1696.   
      DNA   =BBOT-BBTA-(BTOP-BBTA-DNA)*TRAN                             1697.   
      TRBN=TRAN*TRAN                                                    1698.   
      TRB(N)=TRBN                                                       1699.   
      ENB(N)=BTOP+BBTB-(BBOT+BBTB)*TRBN                                 1700.   
      DNB   =BBOT-BBTB-(BTOP-BBTB-DNB)*TRBN                             1701.   
      TRCN=(TRBN*TRBN*TRAN)**2                                          1702.   
      TRC(N)=TRCN                                                       1703.   
      ENC(N)=BTOP+BBTC-(BBOT+BBTC)*TRCN                                 1704.   
      DNC   =BBOT-BBTC-(BTOP-BBTC-DNC)*TRCN                             1705.   
      GO TO 160                                                         1706.   
 150  TRA(N)=0.                                                         1707.   
      TRB(N)=0.                                                         1708.   
      TRC(N)=0.                                                         1709.   
      ENA(N)=BTOP+BBTA                                                  1710.   
      ENB(N)=BTOP+BBTB                                                  1711.   
      ENC(N)=BTOP+BBTC                                                  1712.   
      DNA=BBOT-BBTA                                                     1713.   
      DNB=BBOT-BBTB                                                     1714.   
      DNC=BBOT-BBTC                                                     1715.   
  160 TRDFLB(N)=A*DNA+B*DNB+C*DNC                                       1716.   
      N=N-1                                                             1717.   
      IF(N.GT.0) GO TO 130                                              1718.   
      IF(LTOPCL.LT.1) GO TO 165                                         1719.   
      ENA(LTOPCL)=ENA(LTOPCL)*(1.0-TRCALB(K))+TRCALB(K)*TRDFLB(LTOPCL+1)1720.   
      ENB(LTOPCL)=ENB(LTOPCL)*(1.0-TRCALB(K))+TRCALB(K)*TRDFLB(LTOPCL+1)1721.   
      ENC(LTOPCL)=ENC(LTOPCL)*(1.0-TRCALB(K))+TRCALB(K)*TRDFLB(LTOPCL+1)1722.   
  165 CONTINUE                                                          1723.   
C     ------------------------------------------------------------------1724.   
C                                         SURFACE LAYER FLUX COMPUTATION1725.   
C     ------------------------------------------------------------------1726.   
      N=1                                                               1727.   
      TRDFG=TRDFLB(1)                                                   1728.   
      TAUA=TAUSL(1)+FTAUSL(1)                                           1729.   
      IF(TAUA.GT.1.E-05) GO TO 170                                      1730.   
      BG=BG+TRDFG*TRGALB(K)                                             1731.   
      UNB=BG                                                            1733.   
      UNC=BG                                                            1734.   
      FUNABC=BG                                                         1735.   
      GO TO 180                                                         1736.   
  170 BS=PLANCK(ITS)*WTS1+PLANCK(ITS+1)*WTS                             1737.   
      TA=EXP(-TAUA)                                                     1738.   
      TB=TA*TA                                                          1739.   
      TC=(TB*TB*TA)**2                                                  1740.   
      DNA=(DNA-BS)*TA+BS                                                1741.   
      DNB=(DNB-BS)*TB+BS                                                1742.   
      DNC=(DNC-BS)*TC+BS                                                1743.   
      TRDFG=A*DNA+B*DNB+C*DNC                                           1744.   
      BG=BG+TRDFG*TRGALB(K)                                             1745.   
      UNA=(BG-BS)*TA+BS                                                 1746.   
      UNB=(BG-BS)*TB+BS                                                 1747.   
      UNC=(BG-BS)*TC+BS                                                 1748.   
      FUNABC=A*UNA+B*UNB+C*UNC                                          1749.   
      BSP=PLANCK(ITS+1)*WTS1+PLANCK(ITS+2)*WTS                          1750.   
      BSM=PLANCK(ITS-1)*WTS1+PLANCK(ITS  )*WTS                          1751.   
      SLABS=1.-A*TA-B*TB-C*TC                                           1752.   
      TRSLTS=TRSLTS+(BSP-BSM)*SLABS                                     1753.   
      TRSLTG=TRSLTG+BGFEMD(K)*SLABS                                     1754.   
      TRSLBS=TRSLBS+BS*SLABS                                            1755.   
C     ------------------------------------------------------------------1756.   
C                                                UPWARD FLUX COMPUTATION1757.   
C     ------------------------------------------------------------------1758.   
  180 TRUFLB(N)=FUNABC                                                  1759.   
      IF(N.GT.NLK) GO TO 190                                            1760.   
      UNA=UNA*TRA(N)+ENA(N)                                             1761.   
      UNB=UNB*TRB(N)+ENB(N)                                             1762.   
      UNC=UNC*TRC(N)+ENC(N)                                             1763.   
      FUNABC=A*UNA+B*UNB+C*UNC                                          1764.   
  190 N=N+1                                                             1765.   
      IF(N.LT.NLP) GO TO 180                                            1766.   
      TRUFLB(N)=FUNABC                                                  1767.   
      TRUFTW=FUNABC                                                     1768.   
      TRDFGW=TRDFG                                                      1769.   
      TRUFGW=BG                                                         1770.   
      TRUFG=BG                                                          1771.   
      DO 195 L=1,NLP                                                    1772.   
      DFLB(L,1)=TRDFLB(L)                                               1773.   
  195 UFLB(L,1)=TRUFLB(L)                                               1774.   
      DFSL(1)=TRDFLB(1)                                                 1775.   
      UFSL(1)=TRUFLB(1)                                                 1776.   
      DFLB(1,1)=TRDFGW                                                  1777.   
      UFLB(1,1)=TRUFGW                                                  1778.   
c     print *,' 1778 TRUFLB(1)=',TRUFLB(1)
C     ------------------------------------------------------------------1779.   
C         END WINDOW REGION FLUX COMPUTATION;       CONTINUE INTEGRATION1780.   
C     ------------------------------------------------------------------1781.   
C     ------------------------------------------------------------------1782.   
C          DOWNWARD FLUX COMPUTATION                                    1783.   
C     ------------------------------------------------------------------1784.   
 200  ITK0=K*ITNEXT                                                     1785.   
      K=K+1                                                             1786.   
      IF(K.GT.NKTR) GO TO 300                                           1787.   
      DFLB(NLP,K)=0.                                                    1788.   
      BG=BGFEMT(K)                                                      1789.   
      ITS=ITS+ITNEXT                                                    1790.   
      NLK0=NLK0+NL                                                      1791.   
      NLK=NLK0+NL                                                       1792.   
      NLL=NL                                                            1793.   
  210 TAUA=TAUN(NLK)                                                    1794.   
      IF(TAUA.GT.1.E-05) GO TO 220                                      1795.   
      DFLB(NLL,K)=0.                                                    1796.   
      NLK=NLK-1                                                         1797.   
      NLL=NLL-1                                                         1798.   
      IF(NLL.GT.0) GO TO 210                                            1799.   
      TRUFG=TRUFG+BG                                                    1800.   
      DO 215 N=1,NLP                                                    1801.   
      UFLB(N,K)=BG                                                      1802.   
  215 TRUFLB(N)=TRUFLB(N)+BG                                            1803.   
      GO TO 200                                                         1804.   
  220 N=NLL                                                             1805.   
      DNA=0.                                                            1806.   
      DNB=0.                                                            1807.   
      DNC=0.                                                            1808.   
  230 ITL=ITLT(N)+ITK0                                                  1809.   
      BTOP=PLANCK(ITL)-(PLANCK(ITL)-PLANCK(ITL+1))*WTLT(N)              1810.   
      ITL=ITLB(N)+ITK0                                                  1811.   
      BBOT=PLANCK(ITL)-(PLANCK(ITL)-PLANCK(ITL+1))*WTLB(N)              1812.   
      TAUA=TAUN(NLK)                                                    1813.   
      TAUB=TAUA+TAUA                                                    1814.   
      TAUC=10.*TAUA                                                     1815.   
      IF(TAUA.GT.1.E-01) GO TO 240                                      1816.   
      IF(TAUA.LT.1.E-03) GO TO 235                                      1817.   
      TAU2=TAUA*TAUA                                                    1818.   
      BDIF=BBOT-BTOP                                                    1819.   
      BBTA=BDIF/TAUA                                                    1820.   
      BBTB=BDIF/TAUB                                                    1821.   
      BBTC=BDIF/TAUC                                                    1822.   
      TRAN=1.-TAUA+(0.5-R6*TAUA+R24*TAU2)*TAU2                          1823.   
      GO TO 245                                                         1824.   
  235 BDIF=.5*(BTOP+BBOT)                                               1825.   
      TRA(N)=1.-TAUA                                                    1826.   
      ENA(N)=BDIF*TAUA                                                  1827.   
      DNA=DNA*TRA(N)+ENA(N)                                             1828.   
      TRB(N)=1.-TAUB                                                    1829.   
      ENB(N)=BDIF*TAUB                                                  1830.   
      DNB=DNB*TRB(N)+ENB(N)                                             1831.   
      TRC(N)=1.-TAUC                                                    1832.   
      ENC(N)=BDIF*TAUC                                                  1833.   
      DNC=DNC*TRC(N)+ENC(N)                                             1834.   
      GO TO 260                                                         1835.   
  240 BDIF=BBOT-BTOP                                                    1836.   
      BBTA=BDIF/TAUA                                                    1837.   
      BBTB=BDIF/TAUB                                                    1838.   
      BBTC=BDIF/TAUC                                                    1839.   
      IF(TAUA.GT.7.) GO TO 250                                          1840.   
      TRAN=EXP(-TAUA)                                                   1841.   
  245 TRA(N)=TRAN                                                       1842.   
      ENA(N)=BTOP+BBTA-(BBOT+BBTA)*TRAN                                 1843.   
      DNA   =BBOT-BBTA-(BTOP-BBTA-DNA)*TRAN                             1844.   
      TRBN=TRAN*TRAN                                                    1845.   
      TRB(N)=TRBN                                                       1846.   
      ENB(N)=BTOP+BBTB-(BBOT+BBTB)*TRBN                                 1847.   
      DNB   =BBOT-BBTB-(BTOP-BBTB-DNB)*TRBN                             1848.   
      TRCN=(TRBN*TRBN*TRAN)**2                                          1849.   
      TRC(N)=TRCN                                                       1850.   
      ENC(N)=BTOP+BBTC-(BBOT+BBTC)*TRCN                                 1851.   
      DNC   =BBOT-BBTC-(BTOP-BBTC-DNC)*TRCN                             1852.   
      GO TO 260                                                         1853.   
 250  TRA(N)=0.                                                         1854.   
      TRB(N)=0.                                                         1855.   
      TRC(N)=0.                                                         1856.   
      ENA(N)=BTOP+BBTA                                                  1857.   
      ENB(N)=BTOP+BBTB                                                  1858.   
      ENC(N)=BTOP+BBTC                                                  1859.   
      DNA=BBOT-BBTA                                                     1860.   
      DNB=BBOT-BBTB                                                     1861.   
      DNC=BBOT-BBTC                                                     1862.   
  260 FDNABC=A*DNA+B*DNB+C*DNC                                          1863.   
      TRDFLB(N)=TRDFLB(N)+FDNABC                                        1864.   
      DFLB(N,K)=FDNABC                                                  1865.   
      N=N-1                                                             1866.   
      NLK=NLK-1                                                         1867.   
      IF(N.GT.0) GO TO 230                                              1868.   
      DFSL(K)=FDNABC                                                    1869.   
      IF(LTOPCL.LT.1) GO TO 265                                         1870.   
      ENA(LTOPCL)=ENA(LTOPCL)*(1.0-TRCALB(K))+TRCALB(K)*DFLB(LTOPCL+1,K)1871.   
      ENB(LTOPCL)=ENB(LTOPCL)*(1.0-TRCALB(K))+TRCALB(K)*DFLB(LTOPCL+1,K)1872.   
      ENC(LTOPCL)=ENC(LTOPCL)*(1.0-TRCALB(K))+TRCALB(K)*DFLB(LTOPCL+1,K)1873.   
  265 CONTINUE                                                          1874.   
C     ------------------------------------------------------------------1875.   
C                                         SURFACE LAYER FLUX COMPUTATION1876.   
C     ------------------------------------------------------------------1877.   
      N=1                                                               1878.   
      TAUA=TAUSL(K)+FTAUSL(K)                                           1879.   
      IF(TAUA.GT.1.E-05) GO TO 270                                      1880.   
      BG=BG+FDNABC*TRGALB(K)                                            1881.   
      UNA=BG                                                            1882.   
      UNB=BG                                                            1883.   
      UNC=BG                                                            1884.   
      FUNABC=BG                                                         1885.   
      GO TO 280                                                         1886.   
  270 BS=PLANCK(ITS)*WTS1+PLANCK(ITS+1)*WTS                             1887.   
      TA=EXP(-TAUA)                                                     1888.   
      TB=TA*TA                                                          1889.   
      TC=(TB*TB*TA)**2                                                  1890.   
      DNA=(DNA-BS)*TA+BS                                                1891.   
      DNB=(DNB-BS)*TB+BS                                                1892.   
      DNC=(DNC-BS)*TC+BS                                                1893.   
      FDNABC=A*DNA+B*DNB+C*DNC                                          1894.   
      BG=BGFEMT(K)+FDNABC*TRGALB(K)                                     1895.   
      UNA=(BG-BS)*TA+BS                                                 1896.   
      UNB=(BG-BS)*TB+BS                                                 1897.   
      UNC=(BG-BS)*TC+BS                                                 1898.   
      FUNABC=A*UNA+B*UNB+C*UNC                                          1899.   
      BSP=PLANCK(ITS+1)*WTS1+PLANCK(ITS+2)*WTS                          1900.   
      BSM=PLANCK(ITS-1)*WTS1+PLANCK(ITS  )*WTS                          1901.   
      SLABS=1.-A*TA-B*TB-C*TC                                           1902.   
      TRSLTS=TRSLTS+(BSP-BSM)*SLABS                                     1903.   
      TRSLTG=TRSLTG+BGFEMD(K)*SLABS                                     1904.   
      TRSLBS=TRSLBS+BS*SLABS                                            1905.   
C     ------------------------------------------------------------------1906.   
C                                                UPWARD FLUX COMPUTATION1907.   
C     ------------------------------------------------------------------1908.   
  280 TRUFLB(N)=TRUFLB(N)+FUNABC                                        1909.   
      UFLB(N,K)=FUNABC                                                  1910.   
      IF(N.GT.NLL) GO TO 290                                            1911.   
      UNA=UNA*TRA(N)+ENA(N)                                             1912.   
      UNB=UNB*TRB(N)+ENB(N)                                             1913.   
      UNC=UNC*TRC(N)+ENC(N)                                             1914.   
      FUNABC=A*UNA+B*UNB+C*UNC                                          1915.   
  290 N=N+1                                                             1916.   
      IF(N.LT.NLP) GO TO 280                                            1917.   
      TRUFLB(NLP)=TRUFLB(NLP)+FUNABC                                    1918.   
      UFLB(NLP,K)=FUNABC                                                1919.   
      UFSL(K)=UFLB(1,K)                                                 1920.   
      TRDFG=TRDFG+FDNABC                                                1921.   
      DFLB(1,K)=FDNABC                                                  1922.   
      TRUFG=TRUFG+BG                                                    1923.   
      UFLB(1,K)=BG                                                      1924.   
      IF(K.EQ.11) TRSLWV=TRUFLB(1)-TRDFLB(1)+TRDFG-TRUFG                1925.   
      GO TO 200                                                         1926.   
 300  CONTINUE                                                          1927.   
c     if(PEARTH.gt.0.0.or.PLICE.gt.0.0)then
c     print * ,'1927 JLAT=',JLAT,PEARTH,PLICE
c     print *,' TRUFLB(1)=',TRUFLB(1),' TRUFG=',TRUFG
c     endif

#if ( defined CLM)
c     if(ncallclm.ge.1)then
c     if(PEARTH.gt.0.0.or.PLICE.gt.0.0)then
c     TRUFG=-lwuclm(ILON,JLAT)
c     print *,' CLM TRUFG=',TRUFG
c     endif
c     endif
#endif
C     ------------------------------------------------------------------1928.   
C                                                   END FLUX COMPUTATION1929.   
C     ------------------------------------------------------------------1930.   
      TRSLCR=TRUFLB(1)-TRDFLB(1)+TRDFG-TRUFG                            1931.   
      TRDFSL=TRDFLB(1)                                                  1932.   
      TRDFLB(1)=TRDFG                                                   1933.   
      TRUFSL=TRUFLB(1)                                                  1934.   
      TRUFLB(1)=TRUFG                                                   1935.   
      DO 310 L=1,NLP                                                    1936.   
  310 TRNFLB(L)=TRUFLB(L)-TRDFLB(L)                                     1937.   
      DO 320 L=1,NL                                                     1938.   
  320 TRFCRL(L)=TRNFLB(L+1)-TRNFLB(L)                                   1939.   
      PFW=10.*TRUFTW                                                    1940.   
      IPF=PFW                                                           1941.   
      IF(IPF.LT.10) GO TO 330                                           1942.   
      DPF=PFW-IPF                                                       1943.   
      IPF=IPF+180                                                       1944.   
      GO TO 350                                                         1945.   
  330 PFW=10.*PFW                                                       1946.   
      IPF=PFW                                                           1947.   
      IF(IPF.LT.10) GO TO 340                                           1948.   
      DPF=PFW-IPF                                                       1949.   
      IPF=IPF+90                                                        1950.   
      GO TO 350                                                         1951.   
  340 PFW=10.*PFW                                                       1952.   
      IPF=PFW                                                           1953.   
      IF(IPF.LT.1) IPF=1                                                1954.   
  350 BTEMPW=TKPFW(IPF)+DPF*(TKPFW(IPF+1)-TKPFW(IPF))                   1955.   
      RETURN                                                            1956.   
      END                                                               1957.   
      SUBROUTINE SOLAR                                                  1958.   
C-----------------------------------------------------------------------1959.   
C          SOLAR RETURNS                                                1960.   
C-----------------------------------------------------------------------1961.   
C                       SRDFLB   SOLAR DOWNWARD FLUX AT LAYER BOTTOM    1962.   
C                       SRUFLB   SOLAR UPWARD FLUX AT LAYER BOTTOM EDGE 1963.   
C                       SRNFLB   SOLAR NET (DOWNWARD) FLUX (WATTS/M**2) 1964.   
C                       SRFHRL   SOLAR HEATING RATE : FLUX (WATTS/M**2) 1965.   
C                       SRRVIS   VISALB OF ATMOSPHERE (AS IF RSURFX=0.) 1966.   
C                       SRTATM   ATMOS. TRANSMISSIVITY (TOTAL SPECTRUM) 1967.   
C                       PLAVIS   PLANETARY ALBEDO 0.2-0.7 MICRON REGION 1968.   
C                       ALBVIS   ALBEDO AT GROUND 0.2-0.7 MICRON REGION 1969.   
C                       PLANIR   PLANETARY ALBEDO WAV>0.7 MICRON REGION 1970.   
C                       ALBNIR   ALBEDO AT GROUND WAV>0.7 MICRON REGION 1971.   
C-----------------------------------------------------------------------1972.   
C                COMMENT                                                1973.   
C-----------------------------------------------------------------------1974.   
C                       SOLAR DATA IS RETURNED IN RADCOM LINES:  N,O,P,Q1975.   
C                       NORMS0=1 FLUXES ARE NORMALIZED BY SOLAR CONSTANT1976.   
C                       VERTICAL FLUX DISTRIBUTIONS CONTAIN SOLAR ZENITH1977.   
C                       ANGLE (COSZ) DEPENDENCE                         1978.   
C                       RETURNED SOLAR FLUX VALUES SHOULD BE MULTIPLIED 1979.   
C                       BY COSZ WHEN COMPUTING ATMOSPHERIC HEATING RATE 1980.   
C-----------------------------------------------------------------------1981.   

#include "B83XX.COM"
 
      DIMENSION PFR(52),PFRI(52), PI0C(14),DKS0(14)                     2036.   
      DATA PFR/                                                         2037.   
     1.4144,.4917,.5265,.5530,.5757,.5966,.6159,.6345,.6522,.6689,.6849,2038.   
     2.7003,.7152,.7293,.7428,.7557,.7680,.7796,.7905,.8008,.8105,.8198,2039.   
     3.8286,.8368,.8444,.8515,.8581,.8642,.8699,.8750,.8798,.8843,.8886,2040.   
     4.8928,.8968,.9005,.9040,.9072,.9101,.9129,.9153,.9174,.9193,.9212,2041.   
     5.9227,.9242,.9254,.9266,.9275,.9284,.864245    ,.864245    /      2042.   
      DATA PFRI/                                                        2043.   
     1.4950,.5300,.5620,.5882,.6088,.6302,.6537,.6763,.6969,.7157,.7332,2044.   
     2.7499,.7658,.7806,.7945,.8074,.8194,.8306,.8409,.8504,.8592,.8674,2045.   
     3.8751,.8822,.8886,.8946,.9000,.9050,.9097,.9139,.9177,.9210,.9246,2046.   
     4.9280,.9313,.9343,.9371,.9394,.9415,.9438,.9458,.9475,.9488,.9500,2047.   
     5.9507,.9515,.9529,.9532,.9538,.9541,.876178    ,.876178    /      2048.   
      DATA PI0C/.66,.91,.975,.99,.995,.999,.999,.999,.999,.999,.999,    2049.   
     +          .999,.9999,.99999/                                      2050.   
      DATA DKS0/.01,.03,.04,.04,.04,.002,.004,.013,.002,.003,.003,      2051.   
     +          .072,.20,.53/                                           2052.   
      DIMENSION DBLN(20), KSLAM(14), CPFFL(40)                          2053.   
      DATA DBLN/2.,4.,8.,16.,32.,64.,128.,256.,512.,1024.,2048.,4096.,  2054.   
     +    8192.,16384.,32768.,65536.,131072.,262144.,524288.,1048576./  2055.   
      DATA NKSLAM/14/, KSLAM/1,1,2,2,5,5,5,5,1,1,1,3,4,6/               2056.   
      DATA XCMNO2/5.465/                                                2057.   
      DATA XCMO3/.0399623/                                              2058.   
      DATA TOTRAY/0.000155/                                             2059.   
C                                                                       2060.   
      DIMENSION SRBALB(6),SRXALB(6)                                     2061.   
      EQUIVALENCE (SRBXAL(1,1),SRBALB(1)),(SRBXAL(1,2),SRXALB(1))       2062.   
C                                                                       2063.   
      EQUIVALENCE                                                       2064.   
     + (BXA(1),BOCVIS),(BXA(5),BEAVIS),(BXA( 9),BOIVIS),(BXA(13),BLIVIS)2065.   
     +,(BXA(2),BOCNIR),(BXA(6),BEANIR),(BXA(10),BOINIR),(BXA(14),BLINIR)2066.   
     +,(BXA(3),XOCVIS),(BXA(7),XEAVIS),(BXA(11),XOIVIS),(BXA(15),XLIVIS)2067.   
     +,(BXA(4),XOCNIR),(BXA(8),XEANIR),(BXA(12),XOINIR),(BXA(16),XLINIR)2068.   
     +,               (BXA(17),EXPSNE),(BXA(18),EXPSNO),(BXA(19),EXPSNL)2069.   
     +,               (BXA(20),BSNVIS),(BXA(21),BSNNIR)                 2070.   
     +,               (BXA(22),XSNVIS),(BXA(23),XSNNIR)                 2071.   
C                                                                       2072.   
      EQUIVALENCE (SRXATM(1),SRXVIS),(SRXATM(2),SRXNIR)                 2073.   
      EQUIVALENCE (SRXATM(3),XXAVIS),(SRXATM(4),XXANIR)                 2074.   
C                                                                       2075.   
      EQUIVALENCE (ISPARE(1),NEWASZ)                                    2075.5  
C                                                                       2076.   
C-----------------------------------------------------------------------2077.   
C     SOLAR: NET FLUX AT GROUND FOR FRACTIONAL GRID SURFACE ALBEDOS     2078.   
C                                                                       2079.   
C     PFNFG(DT,XA,RSA,RX,RB)=(DT*(1.-RB)-XA*(RX-RB)*(1.-RSA))           2080.   
C    +                      /(1.-RSA*RB)                                2081.   
C-----------------------------------------------------------------------2082.   
C                                                                       2083.   
C                                                                       2084.   
C     O3ABS(X)= 1.08173*X/(1.00+                                        2085.   
C    $ 138.57*X)**0.805 + 0.0658*X/(1.00+(103.63*X)**3)                 2086.   
C                                                                       2087.   
      S0COSZ=S0                                                         2088.   
      IF(NORMS0.EQ.0) S0COSZ=S0*COSZ                                    2089.   
C                                                                       2090.   
      DO 10 N=1,NLP                                                     2091.   
      SRNFLB(N)=0.                                                      2092.   
      SRDFLB(N)=0.                                                      2093.   
      SRUFLB(N)=0.                                                      2094.   
      SRFHRL(N)=0.                                                      2095.   
   10 CONTINUE                                                          2096.   
      SRIVIS=0.                                                         2097.   
      SROVIS=0.                                                         2098.   
      SRINIR=0.                                                         2099.   
      SRONIR=0.                                                         2100.   
      SRDVIS=0.                                                         2101.   
      SRUVIS=0.                                                         2102.   
      SRDNIR=0.                                                         2103.   
      SRUNIR=0.                                                         2104.   
      SRTVIS=0.                                                         2105.   
      SRAVIS=0.                                                         2106.   
      SRTNIR=0.                                                         2107.   
      SRANIR=0.                                                         2108.   
      SRSLHR=0.                                                         2109.   
      PLAVIS=1.                                                         2110.   
      PLANIR=1.                                                         2111.   
      ALBVIS=1.                                                         2112.   
      ALBNIR=1.                                                         2113.   
      SRRVIS=1.                                                         2114.   
      SRRNIR=0.                                                         2115.   
      SRTNIR=0.                                                         2116.   
      SRXVIS=0.                                                         2117.   
      SRXNIR=0.                                                         2118.   
C                                                                       2119.   
      XXVIS=.53/(1.-SRBALB(6))                                          2120.   
      XXNIR=.47/(1.-SRBALB(5))                                          2121.   
      DO 20 N=1,4                                                       2122.   
   20 FSRNFG(N)=XXVIS*(1.-BXA(4*N-3))+XXNIR*(1.-BXA(4*N-2))             2123.   
C                                                                       2124.   
      IF(COSZ.LT.0.01) RETURN                                           2125.   
      COSMAG=35.0/SQRT(1224.*COSZ*COSZ+1.0)                             2126.   
      TAURAY=TOTRAY*FRAYLE                                              2127.   
      CPF=49.999/COSMAG                                                 2128.   
      IPF=CPF                                                           2129.   
      DPF=CPF-IPF                                                       2130.   
      IF(ISOSCT.EQ.1) IPF=51                                            2131.   
      CPFF=(1.0-DPF)*PFR(IPF)+DPF*PFR(IPF+1)                            2132.   
      CPFFI=(1.0-DPF)*PFRI(IPF)+DPF*PFRI(IPF+1)                         2133.   
      SECZ=1./COSZ                                                      2134.   
      DO 100 N=1,NL                                                     2135.   
      CPFFL(N)=CPFF                                                     2136.   
      IF(TLM(N).LT.TKCICE) CPFFL(N)=CPFFI                               2137.   
 100  CONTINUE                                                          2138.   
C                                                                       2139.   
      K = 0                                                             2140.   
  300 K = K+1                                                           2141.   
C                                                                       2142.   
      KLAM=KSLAM(K)                                                     2143.   
      DKS0K=DKS0(K)                                                     2144.   
      DKS0X=DKS0K*S0COSZ                                                2145.   
      RBNB=SRBALB(KLAM)                                                 2146.   
      RBNX=SRXALB(KLAM)                                                 2147.   
      RCNB=0.0                                                          2148.   
      RCNX=0.0                                                          2149.   
C                                                                       2150.   
      N = 0                                                             2151.   
  200 N = N+1                                                           2152.   
C                                                                       2153.   
      CPFF=CPFFL(N)                                                     2154.   
      SRB(N)=RBNB                                                       2155.   
      SRX(N)=RBNX                                                       2156.   
      TLN=TLM(N)                                                        2157.   
      PLN=PL(N)                                                         2158.   
      ULN=ULGAS(N,1)                                                    2159.   
      RTAU=1.E-06                                                       2160.   
      GO TO (101,102,103,104,105,106,107,108,109,110,111,112,113,114),K 2161.   
 101  CONTINUE                                                          2162.   
C--------K=6-------H2O       DS0=.01                                    2163.   
      TERMA=(35.66+TLN*(.0416-.0004622*TLN+.001057*PLN))*(1.+.04286*PLN)2164.   
      TERMB=(1.+.00171*ULN)*(1.+PLN*(189.088+.1316*PLN))                2165.   
      TAU1  =TERMA/TERMB                                                2166.   
      IF(TAU1.GT.0.02343) TAU1=0.02343                                  2167.   
      TAU=TAU1*ULN                                                      2168.   
      GO TO 120                                                         2169.   
 102  CONTINUE                                                          2170.   
C--------K=5-------H2O       DS0=.03                                    2171.   
      TERMA=(2.792+TLN*(.0914-.0002848*TLN+.0003395*PLN))               2172.   
     +     *(1.+.02964*PLN)                                             2173.   
      TERMB=(1.0+.000657*ULN)*(1.+PLN*(240.70+.13847*PLN))              2174.   
      TAU1  =TERMA/TERMB                                                2175.   
      IF(TAU1.GT.0.00520) TAU1=0.00520                                  2176.   
      TAU=TAU1*ULN                                                      2177.   
      GO TO 120                                                         2178.   
 103  CONTINUE                                                          2179.   
C--------K=4-------H2O       DS0=.04                                    2180.   
      TERMA=(.4768+.467E-04*PLN*TLN)*(1.+TLN*(.00191-.719E-05*TLN))     2181.   
      TERMB=(1.+.717E-04*ULN)*(1.+PLN*(130.56+.0876*PLN))/(1.+.0266*PLN)2182.   
      TAU1  =TERMA/TERMB                                                2183.   
      IF(TAU1.GT.0.00150) TAU1=0.0015                                   2184.   
      TAU=TAU1*ULN                                                      2185.   
      GO TO 120                                                         2186.   
 104  CONTINUE                                                          2187.   
C--------K=3-------H2O       DS0=.04                                    2188.   
      TERMA=(.000247*TLN-.091+PLN*(.00035+.78E-06*TLN))*(1.+.2847*PLN)  2189.   
      TERMB=(1.+.2066E-04*ULN)*(1.+PLN*(137.17+.16132*PLN))             2190.   
      TAU   =(TERMA/TERMB)*ULN                                          2191.   
      GO TO 120                                                         2192.   
 105  CONTINUE                                                          2193.   
C--------K=2-------H2O       DS0=.04                                    2194.   
      TERMA=(PLN*(1.974/TLN+.0001117*TLN)-10.713)*(1.+.005788*TLN)      2195.   
     +     *(1.+.001517*PLN)                                            2196.   
      TERMB=(1.+.3218E-04*ULN)*(1.+PLN*(863.44+.2048*PLN))              2197.   
      TAU   =(TERMA/TERMB)*ULN                                          2198.   
      GO TO 120                                                         2199.   
 106  CONTINUE                                                          2200.   
C--------K=4-------O2        DS0=.002                                   2201.   
      ULN=ULGAS(N,4)                                                    2202.   
      TERMA=(.2236E-05-.1181E-09*TLN)*(1.+PLN*(.6364E-05*PLN+.001168))  2203.   
      TERMB=1.+.1521E-05*ULN                                            2204.   
      TAU   =(TERMA/TERMB)*ULN                                          2205.   
      GO TO 120                                                         2206.   
 107  CONTINUE                                                          2207.   
C--------K=3-------O2        DS0=.004                                   2208.   
      ULN=ULGAS(N,4)                                                    2209.   
      TERMA=(.3179E-06-.9263E-11*TLN)*(1.+PLN*(.8832E-05*PLN+.0005292)) 2210.   
      TERMB=1.+.1968E-06*ULN                                            2211.   
      TAU   =(TERMA/TERMB)*ULN                                          2212.   
      GO TO 120                                                         2213.   
 108  CONTINUE                                                          2214.   
C--------K=2-------O2        DS0=.013                                   2215.   
      ULN=ULGAS(N,4)                                                    2216.   
      TERMA=(.2801E-07-.1638E-12*TLN)*(1.+PLN*(.1683E-04*PLN-.001721))  2217.   
      TERMB=1.+.8097E-07*ULN                                            2218.   
      TAU   =(TERMA/TERMB)*ULN                                          2219.   
      GO TO 120                                                         2220.   
 109  CONTINUE                                                          2221.   
C--------K=4-------CO2       DS0=.002                                   2222.   
      ULN=ULGAS(N,2)                                                    2223.   
      TERMA=(50.73-.03155*TLN-PLN*(.5543+.00091*TLN))*(1.-.1004*PLN)    2224.   
      TERMB=(1.+.006468*ULN)*(1.+PLN*(49.51+.8285*PLN))                 2225.   
      TAU   =(TERMA/TERMB)*ULN                                          2226.   
      IF(PLN.LT.175.0) TAU=(.00018*PLN+0.00001)*ULN                     2227.   
      GO TO 120                                                         2228.   
 110  CONTINUE                                                          2229.   
C--------K=3-------CO2       DS0=.003                                   2230.   
      ULN=ULGAS(N,2)                                                    2231.   
      TERMA=(1.+.01319*TLN)*(PLN*(.008001*ULN+.4589E-03)-.8396*ULN)     2232.   
      TERMB=ULN*(PLN+295.7+1.967*ULN)+.15126*PLN                        2233.   
      TAU   =(TERMA/TERMB)*ULN                                          2234.   
      GO TO 120                                                         2235.   
 111  CONTINUE                                                          2236.   
C--------K=2-------CO2       DS0=.003                                   2237.   
      ULN=ULGAS(N,2)                                                    2238.   
      TERMA=(1.+.02257*TLN)*(PLN*(.002295*ULN-.5489E-04)-.7571*ULN)     2239.   
      TERMB=ULN*(PLN+803.9+2.477*ULN)-.09899*PLN                        2240.   
      TAU   =(TERMA/TERMB)*ULN                                          2241.   
      GO TO 120                                                         2242.   
 112  CONTINUE                                                          2243.   
      TAU=0.0                                                           2244.   
      GO TO 120                                                         2245.   
 113  CONTINUE                                                          2246.   
      TAU=0.0                                                           2247.   
      GO TO 120                                                         2248.   
 114  CONTINUE                                                          2249.   
      TAU=XCMNO2*ULGAS(N,5)+XCMO3*ULGAS(N,3)                            2250.   
      RTAU=TAURAY*(PLB(N)-PLB(N+1))                                     2251.   
 120  CONTINUE                                                          2252.   
      IF(TAU.LT.0.0) TAU=0.0                                            2253.   
      CTAU=CLDTAU(N)*FCLDSR                                             2254.   
      CPI0=PI0C(K)                                                      2255.   
      ATAU=EXTAER(N,KLAM)                                               2256.   
      TAU=TAU+CTAU+ATAU+RTAU                                            2257.   
      IF(TAU.LT.TAUMIN) GO TO 180                                       2258.   
      CTAUSC=CPI0*CTAU                                                  2259.   
      ATAUSC=SCTAER(N,KLAM)                                             2260.   
      TAUSCT=CTAUSC+ATAUSC+RTAU                                         2261.   
      PIZERO=TAUSCT/TAU                                                 2262.   
      IF(PIZERO.GT.0.001) GO TO 130                                     2263.   
      GO TO 180                                                         2264.   
 130  CONTINUE                                                          2265.   
      APFF=COSAER(N,KLAM)                                               2266.   
      APFF0=APFF                                                        2266.1  
      IF(NEWASZ.GT.0) CALL HGAER1(COSZ,ATAUSC,APFF0,APFF)               2266.2  
      PFF=(CPFF*CTAUSC+APFF*ATAUSC)/TAUSCT                              2267.   
      IF(ISOSCT.GT.1) GO TO 131                                         2268.   
      GO TO 132                                                         2269.   
  131 TAU=TAU-TAUSCT*PFF                                                2270.   
      PIZERO=PIZERO*(1.-PFF)/(1.-PIZERO*PFF)                            2271.   
      PFF=0.                                                            2272.   
  132 CONTINUE                                                          2273.   
      PR=1.0-PFF                                                        2274.   
      PT=1.0+PFF                                                        2275.   
      IF(TAU.LT.0.015625) GO TO 140                                     2276.   
C ALOG
      DBLS=7.001+1.44269*LOG(TAU)                                       2277.   
C ALOG
      NDBLS=DBLS                                                        2278.   
      TAU=TAU/DBLN(NDBLS)                                               2279.   
      GO TO 150                                                         2280.   
 140  XANB=EXP(-TAU-TAU)                                                2281.   
      XANX=EXP(-TAU*SECZ)                                               2282.   
      TANB=PT*XANB                                                      2283.   
      XXT=(SECZ-2.0)*TAU                                                2284.   
      TANX=PT*SECZ*(.5+XXT*(.25+XXT*(.0833333+XXT*(.0208333+XXT))))*XANX2285.   
      RASB=PR*(1.0-TAU*(2.0-2.66667*TAU*(1.0-TAU)))                     2286.   
      XXT=(SECZ+2.0)*TAU                                                2287.   
      RASX=PR*SECZ*(.5-XXT*(.25-XXT*(.0833333-XXT*(.0208333-XXT))))     2288.   
      BNORM=(1.0-XANB)/(RASB+TANB)*PIZERO                               2289.   
      XNORM=(1.0-XANX)/(RASX+TANX)*PIZERO                               2290.   
      RASB=RASB*BNORM                                                   2291.   
      RASX=RASX*XNORM                                                   2292.   
      TANB=TANB*BNORM                                                   2293.   
      TANX=TANX*XNORM                                                   2294.   
      GO TO 170                                                         2295.   
 150  XANB=EXP(-TAU-TAU)                                                2296.   
      XANX=EXP(-TAU*SECZ)                                               2297.   
      TANB=PT*XANB                                                      2298.   
      XXT=(SECZ-2.0)*TAU                                                2299.   
      TANX=PT*SECZ*(.5+XXT*(.25+XXT*(.0833333+XXT*(.0208333+XXT))))*XANX2300.   
      RASB=PR*(1.0-TAU*(2.0-2.66667*TAU*(1.0-TAU)))                     2301.   
      XXT=(SECZ+2.0)*TAU                                                2302.   
      RASX=PR*SECZ*(.5-XXT*(.25-XXT*(.0833333-XXT*(.0208333-XXT))))     2303.   
      BNORM=(1.0-XANB)/(RASB+TANB)*PIZERO                               2304.   
      XNORM=(1.0-XANX)/(RASX+TANX)*PIZERO                               2305.   
      RASB=RASB*BNORM                                                   2306.   
      RASX=RASX*XNORM                                                   2307.   
      TANB=TANB*BNORM                                                   2308.   
      TANX=TANX*XNORM                                                   2309.   
      DO 160 NN=1,NDBLS                                                 2310.   
      RARB=RASB*RASB                                                    2311.   
      RARX=XANX*RASX                                                    2312.   
      XATB=XANB+TANB                                                    2313.   
      DENOM=1.0-RARB                                                    2314.   
      DB=(TANB+XANB*RARB)/DENOM                                         2315.   
      DX=(TANX+RARX*RASB)/DENOM                                         2316.   
      UB=RASB*(XANB+DB)                                                 2317.   
      UX=RARX+RASB*DX                                                   2318.   
      RASB=RASB+XATB*UB                                                 2319.   
      RASX=RASX+XATB*UX                                                 2320.   
      TANB=XANB*TANB+XATB*DB                                            2321.   
      TANX=XANX*TANX+XATB*DX                                            2322.   
      XANB=XANB*XANB                                                    2323.   
      XANX=XANX*XANX                                                    2324.   
 160  CONTINUE                                                          2325.   
 170  RARB=RASB*RBNB                                                    2326.   
      RARX=RASB*RBNX                                                    2327.   
      XATB=XANB+TANB                                                    2328.   
      DENOM=1.0-RARB                                                    2329.   
      DB=(TANB+XANB*RARB)/DENOM                                         2330.   
      DX=(TANX+XANX*RARX)/DENOM                                         2331.   
      UB=RBNB*(XANB+DB)                                                 2332.   
      UX=RBNX*XANX+RBNB*DX                                              2333.   
      RBNB=RASB+XATB*UB                                                 2334.   
      RBNX=RASX+XATB*UX                                                 2335.   
      XATC=XATB/(1.0-RASB*RCNB)                                         2336.   
      RCNX=RASX+(XANX*RCNX+TANX*RCNB)*XATC                              2337.   
      RCNB=RASB+RCNB*XATB*XATC                                          2338.   
      GO TO 190                                                         2339.   
 180  RASB=0.0                                                          2340.   
      RASX=0.0                                                          2341.   
      TANB=0.0                                                          2342.   
      TANX=0.0                                                          2343.   
      XANB=EXP(-TAU-TAU)                                                2344.   
      XANX=EXP(-TAU*SECZ)                                               2345.   
      DX=0.0                                                            2346.   
      UX=RBNX*XANX                                                      2347.   
      RBNB=RBNB*XANB*XANB                                               2348.   
      RBNX=UX*XANB                                                      2349.   
      RCNB=RCNB*XANB*XANB                                               2350.   
      RCNX=RCNX*XANX*XANB                                               2351.   
 190  RNB(N)=RASB                                                       2352.   
      RNX(N)=RASX                                                       2353.   
      TNB(N)=TANB                                                       2354.   
      TNX(N)=TANX                                                       2355.   
      XNB(N)=XANB                                                       2356.   
      XNX(N)=XANX                                                       2357.   
      IF(N.LT.NL) GO TO 200                                             2358.   
C                                                                       2359.   
      IF(K.EQ.NKSLAM) GO TO 301                                         2360.   
      SRDFLB(NLP)=SRDFLB(NLP)+DKS0X                                     2361.   
      SRUFLB(NLP)=SRUFLB(NLP)+DKS0X*RBNX                                2362.   
      SRDFLB(NL)=SRDFLB(NL)+DKS0X*(XANX+DX)                             2363.   
      SRUFLB(NL)=SRUFLB(NL)+DKS0X*UX                                    2364.   
      RMEAN=RBNX                                                        2365.   
      DO 230 M=2,NL                                                     2366.   
      N=NLP-M                                                           2367.   
      XBNB=XNB(N)                                                       2368.   
      XBNX=XNX(N)                                                       2369.   
      RBNX=RNX(N)                                                       2370.   
      IF(RBNX.GT.1.E-05) GO TO 210                                      2371.   
      RASB=RASB*XBNB*XBNB                                               2372.   
      TANX=TANX*XBNB                                                    2373.   
      GO TO 220                                                         2374.   
 210  RBNB=RNB(N)                                                       2375.   
      TBNB=TNB(N)                                                       2376.   
      TBNX=TNX(N)                                                       2377.   
      RARB=RASB*RBNB                                                    2378.   
      XBTB=XBNB+TBNB                                                    2379.   
      DENOM=1.0-RARB                                                    2380.   
      TANX=TBNX*XANX+XBTB*(TANX+XANX*RBNX*RASB)/DENOM                   2381.   
      RASB=RBNB+XBTB*XBTB*RASB/DENOM                                    2382.   
 220  XANX=XANX*XBNX                                                    2383.   
      RBNB=SRB(N)                                                       2384.   
      RBNX=SRX(N)                                                       2385.   
      DX=(TANX+XANX*RBNX*RASB)/(1.0-RASB*RBNB)                          2386.   
      UX=RBNX*XANX+RBNB*DX                                              2387.   
      SRUFLB(N)=SRUFLB(N)+DKS0X*UX                                      2388.   
  230 SRDFLB(N)=SRDFLB(N)+DKS0X*(XANX+DX)                               2389.   
      SRRNIR=SRRNIR+DKS0K*RCNX                                          2390.   
      SRTNIR=SRTNIR+DKS0K*(TANX+XANX)                                   2391.   
      SRXNIR=SRXNIR+DKS0K*XANX                                          2392.   
      GO TO 300                                                         2393.   
C                                                                       2394.   
  301 CONTINUE                                                          2395.   
      SRTNIR=SRTNIR/0.459                                               2396.   
      SRRNIR=SRRNIR/0.459                                               2397.   
      SRXNIR=SRXNIR/0.459                                               2398.   
      SRANIR=1.0-SRTNIR-SRRNIR                                          2399.   
C                                                                       2400.   
      VRD(NLP)=DKS0X                                                    2401.   
      VRU(NLP)=DKS0X*RBNX                                               2402.   
      O3PATH=(1.9+XANX*(COSMAG-1.9))*ULGAS(NL,3)                        2403.   
      ATOP=0.                                                           2404.   
      ABOT=O3ABS(O3PATH)                                                2405.   
      ASUM=(ABOT-ATOP)*XANX                                             2406.   
      O3A(NL)=ASUM*S0COSZ                                               2407.   
      ATOP=ABOT                                                         2408.   
      VRD(NL)=DKS0X*(XANX+DX)-ASUM*S0COSZ                               2409.   
      VRU(NL)=DKS0X*UX                                                  2410.   
      FAC(NL)=UX                                                        2411.   
      RMEAN=RBNX                                                        2412.   
      N=NL                                                              2413.   
  305 N=N-1                                                             2414.   
      XBNB=XNB(N)                                                       2415.   
      XBNX=XNX(N)                                                       2416.   
      RBNX=RNX(N)                                                       2417.   
      IF(RBNX.GT.1.E-05) GO TO 310                                      2418.   
      RASB=RASB*XBNB*XBNB                                               2419.   
      TANX=TANX*XBNB                                                    2420.   
      GO TO 320                                                         2421.   
 310  RBNB=RNB(N)                                                       2422.   
      TBNB=TNB(N)                                                       2423.   
      TBNX=TNX(N)                                                       2424.   
      RARB=RASB*RBNB                                                    2425.   
      XBTB=XBNB+TBNB                                                    2426.   
      DENOM=1.0-RARB                                                    2427.   
      TANX=TBNX*XANX+XBTB*(TANX+XANX*RBNX*RASB)/DENOM                   2428.   
      RASB=RBNB+XBTB*XBTB*RASB/DENOM                                    2429.   
 320  XANX=XANX*XBNX                                                    2430.   
      RBNB=SRB(N)                                                       2431.   
      RBNX=SRX(N)                                                       2432.   
      DX=(TANX+XANX*RBNX*RASB)/(1.0-RASB*RBNB)                          2433.   
      UX=RBNX*XANX+RBNB*DX                                              2434.   
      FAC(N)=UX                                                         2435.   
      VRU(N)=DKS0X*UX                                                   2436.   
      O3PATH=O3PATH+(1.9+XANX*(COSMAG-1.9))*ULGAS(N,3)                  2437.   
      ABOT=O3ABS(O3PATH)                                                2438.   
      ASUM=ASUM+(ABOT-ATOP)*XANX                                        2439.   
      ATOP=ABOT                                                         2440.   
      VRD(N)=DKS0X*(XANX+DX)-ASUM*S0COSZ                                2441.   
      O3A(N)=ASUM*S0COSZ                                                2442.   
      IF(N.GT.1) GO TO 305                                              2443.   
C                                                                       2444.   
      O3SUM=0.                                                          2445.   
      DO 324 I=1,NL                                                     2446.   
  324 O3SUM=O3SUM+ULGAS(I,3)                                            2447.   
      SRXVIS=XANX*(1.-O3ABS(COSMAG*O3SUM)/0.53)                         2448.   
      SRTVIS=TANX+XANX-ASUM/DKS0K                                       2449.   
      RGRND=UX/(XANX+DX+1.E-05)                                         2450.   
      IF(RGRND.GT.1.0) RGRND=1.0                                        2451.   
      ASUM=ASUM*RGRND                                                   2452.   
      VRU(N)=VRU(N)-ASUM*S0COSZ                                         2453.   
  325 CONTINUE                                                          2454.   
      O3PATH=O3PATH+1.9*ULGAS(N,3)                                      2455.   
      ATOP=O3ABS(O3PATH)                                                2456.   
      ASUM=ASUM+(ATOP-ABOT)*FAC(N)                                      2457.   
      ABOT=ATOP                                                         2458.   
      N=N+1                                                             2459.   
      VRU(N)=VRU(N)-ASUM*S0COSZ                                         2460.   
      IF(N.LT.NLP) GO TO 325                                            2461.   
      SRRVIS=RCNX-ASUM/DKS0K                                            2462.   
      SRAVIS=1.0-SRRVIS-SRTVIS                                          2463.   
      TFU=VRU(NLP)                                                      2464.   
      BFU=VRU(1)                                                        2465.   
      IF(BFU.GE.0.) GO TO 327                                           2466.   
      DO 326 N=1,NLP                                                    2467.   
  326 VRU(N)=(VRU(N)-BFU)*(TFU/(TFU-BFU))                               2468.   
      BFU=VRU(1)                                                        2469.   
  327 BFD=VRD(1)                                                        2470.   
      IF(BFD.GT.BFU) GO TO 329                                          2471.   
      TFD=VRD(NLP)                                                      2472.   
      BFUD=BFU/TFD                                                      2473.   
      TFDD=TFD/(TFD-BFD)                                                2474.   
      DO 328 N=1,NLP                                                    2475.   
  328 VRD(N)=(VRD(N)*(1.-BFUD)-BFD+BFUD*TFD)*TFDD                       2476.   
  329 SRDVIS=VRD(1)                                                     2477.   
      SRUVIS=VRU(1)                                                     2478.   
      ALBVIS=SRUVIS/(SRDVIS+1.E-10)                                     2479.   
      TAU1=0.                                                           2480.   
      SRIVIS=VRD(NLP)                                                   2481.   
      SROVIS=VRU(NLP)                                                   2482.   
      PLAVIS=SROVIS/SRIVIS                                              2483.   
C                                                                       2484.   
      TAU2=0.                                                           2485.   
      TAU3=0.                                                           2486.   
      TRN1=0.                                                           2487.   
      TRN2=0.                                                           2488.   
      TRN3=0.                                                           2489.   
      N=NLP                                                             2490.   
C                                                                       2491.   
C             THE FOLLOWING IS CONSIDERED PART OF THE NEAR-IR SPECTRUM  2492.   
C             --------------------------------------------------------  2493.   
      DO 330 M=1,NL                                                     2494.   
      N=N-1                                                             2495.   
      PLN=PL(N)                                                         2496.   
      ULN=ULGAS(N,2)*SECZ                                               2497.   
      ULX=ULN                                                           2498.   
      IF(ULN.GT.7.0) ULN=7.0                                            2499.   
C--------K=5-------CO2       DS0=.002                                   2500.   
      TERMA=.003488*PLN*(1.+39.59*EXP(-8.769*ULN/(1.+4.419*ULN)))       2501.   
     +     *(1.+ULN*(.001938*PLN-.00503*ULN))                           2502.   
      TERMB=(1.+.04712*PLN*(1.+.4877*ULN))                              2503.   
      TAU=TERMA/TERMB                                                   2504.   
      IF(TAU.LT.1.E-06) TAU=1.E-06                                      2505.   
      TAU1=TAU1+TAU*ULX                                                 2506.   
      ULN=ULGAS(N,1)*SECZ                                               2507.   
C--------K=7-------H2O       DS0=.01(DS0=.008 + DS0=.002 CO2 OVERLAP)   2508.   
      TERMA=.001582*PLN*(1.+6.769*EXP(-9.59*ULN/(1.+5.026*ULN)))        2509.   
     +     *(1.+ULN*(.2757E-03*PLN+.001429*ULN))                        2510.   
      TERMB=(1.+.003683*PLN*(1.+1.187*ULN))                             2511.   
      TAU2=TAU2+(TERMA/TERMB)*ULN                                       2512.   
      ULN=ULGAS(N,4)*SECZ                                               2513.   
C--------K=5-------O2        DS0=.001                                   2514.   
      TERMA=(.1366E-03-.2203E-07*TLN)*(1.+PLN*(.1497E-06*ULN+.001261))  2515.   
      TERMB=(1.+.3867E-03*ULN)/(1.+.2075E-04*ULN)                       2516.   
      TAU3=TAU3+(TERMA/TERMB)*ULN                                       2517.   
      IF(TAU1.LT.10.0) TRN1=EXP(-TAU1)                                  2518.   
      IF(TAU2.LT.10.0) TRN2=EXP(-TAU2)                                  2519.   
      IF(TAU3.LT.10.0) TRN3=EXP(-TAU3)                                  2520.   
      FAC(N)=.004358*TRN1+.01743*TRN2+.00218*TRN3                       2521.   
  330 SRDFLB(N)=SRDFLB(N)+SRDFLB(N)*FAC(N)                              2522.   
      FAC(NLP)=.023968                                                  2523.   
      SRDFLB(NLP)=SRDFLB(NLP)+SRDFLB(NLP)*FAC(NLP)                      2524.   
      DO 340 N=1,NLP                                                    2525.   
  340 SRUFLB(N)=SRUFLB(N)+SRUFLB(N)*FAC(1)                              2526.   
      SRINIR=SRDFLB(NLP)                                                2527.   
      SRONIR=SRUFLB(NLP)                                                2528.   
      PLANIR=SRONIR/SRINIR                                              2529.   
      SRDNIR=SRDFLB(1)                                                  2530.   
      SRUNIR=SRUFLB(1)                                                  2531.   
      ALBNIR=SRUNIR/(SRDNIR+1.E-10)                                     2532.   
      DO 350 N=1,NLP                                                    2533.   
      SRDFLB(N)=SRDFLB(N)+VRD(N)                                        2534.   
      SRUFLB(N)=SRUFLB(N)+VRU(N)                                        2535.   
  350 SRNFLB(N)=SRDFLB(N)-SRUFLB(N)                                     2536.   
      DO 360 N=1,NL                                                     2537.   
  360 SRFHRL(N)=SRNFLB(N+1)-SRNFLB(N)                                   2538.   
      SRSLHR=FRACSL*SRFHRL(1)                                           2539.   
C                                                                       2540.   
C---------------------------------                                      2541.   
      CALL O2HEAT(FAC,COSZ,S0COSZ)                                      2542.   
C---------------------------------                                      2543.   
C                                                                       2544.   
      DO 500 L=1,NL                                                     2545.   
  500 SRFHRL(L)=SRFHRL(L)+FAC(L)                                        2546.   
      L=NLP                                                             2547.   
      DO 510 N=1,NL                                                     2548.   
      L=L-1                                                             2549.   
      IF(PLB(L).GT.0.09) GO TO 520                                      2550.   
  510 SRFHRL(L)=FAC(L)+O3A(L)                                           2551.   
  520 CONTINUE                                                          2552.   
C     I=NLP+1-II                                                        2553.   
C                                                                       2554.   
C-----------------------------------------------------------------------2555.   
C     SOLAR NET FLUX (SRNFLB(1)) DISTRIBUTION ACCORDING TO SURFACE TYPE 2556.   
CR    NOT USED AND NOT SAFE (CAUSES DIVIDE CHECKS)                      2556.1  
C-----------------------------------------------------------------------2557.   
CR    FSRVIS=0.53                                                       2558.   
CR    FSRNIR=0.47                                                       2559.   
C                                                                       2560.   
CR    RASVIS=0.                                                         2561.   
CR    IF(SRUVIS.GT.1.E-03) RASVIS=(SRDVIS-SRTVIS*SRIVIS)/SRUVIS         2562.   
CR    XXAVIS=0.                                                         2563.   
CR    DENOM=SRIVIS*(SRXALB(6)-SRBALB(6))                                2564.   
CR    IF(ABS(DENOM).GT.1.E-03) XXAVIS=(SRUVIS-SRDVIS*SRBALB(6))/DENOM   2565.   
C$    PNFVIS=(SRDVIS-SRUVIS)/SRIVIS                                     2566.   
CR    IF(SRIVIS.GT.1.E-03) PNFVIS=(SRDVIS-SRUVIS)/SRIVIS                2566.11 
CR    RASNIR=0.                                                         2567.   
CR    IF(PNFVIS.LT.1.E-03) RETURN                                       2568.   
CR    IF(SRUNIR.GT.1.E-03) RASNIR=(SRDNIR-SRTNIR*SRINIR)/SRUNIR         2569.   
CR    XXANIR=0.                                                         2570.   
CR    DENOM=SRINIR*(SRXALB(5)-SRBALB(5))                                2571.   
CR    IF(ABS(DENOM).GT.1.E-03) XXANIR=(SRUNIR-SRDNIR*SRBALB(5))/DENOM   2572.   
C$    PNFNIR=(SRDNIR-SRUNIR)/SRINIR                                     2573.   
CR    IF(SRINIR.GT.1.E-03) PNFNIR=(SRDNIR-SRUNIR)/SRINIR                2573.11 
CR    IF(PNFNIR.LT.1.E-03) RETURN                                       2574.   
C                                                                       2575.   
CR    FNSROC=0.                                                         2576.   
CR    IF(POCEAN.LT.1.E-04) GO TO 601                                    2577.   
CR    POCVIS=PFNFG(SRTVIS,SRXVIS,RASVIS,XOCVIS,BOCVIS)                  2578.   
CR    POCNIR=PFNFG(SRTNIR,SRXNIR,RASNIR,XOCVIS,BOCVIS)                  2579.   
CR    FNSROC=(FSRVIS*POCVIS/PNFVIS+FSRNIR*POCNIR/PNFNIR)                2580.   
C                                                                       2581.   
CR601 FNSREA=0.                                                         2582.   
CR    IF(PEARTH.LT.1.E-04) GO TO 602                                    2583.   
CR    PEAVIS=PFNFG(SRTVIS,SRXVIS,RASVIS,XEAVIS,BEAVIS)                  2584.   
CR    PEANIR=PFNFG(SRTNIR,SRXNIR,RASNIR,XEANIR,BEANIR)                  2585.   
CR    FNSREA=(FSRVIS*PEAVIS/PNFVIS+FSRNIR*PEANIR/PNFNIR)                2586.   
C                                                                       2587.   
CR602 FNSROI=0.                                                         2588.   
CR    IF(POICE .LT.1.E-04) GO TO 603                                    2589.   
CR    POIVIS=PFNFG(SRTVIS,SRXVIS,RASVIS,XOIVIS,BOIVIS)                  2590.   
CR    POINIR=PFNFG(SRTNIR,SRXNIR,RASNIR,XOINIR,BOINIR)                  2591.   
CR    FNSROI=(FSRVIS*POIVIS/PNFVIS+FSRNIR*POINIR/PNFNIR)                2592.   
C                                                                       2593.   
CR603 FNSRLI=0.                                                         2594.   
CR    IF(PLICE .LT.1.E-04) GO TO 604                                    2595.   
CR    PLIVIS=PFNFG(SRTVIS,SRXVIS,RASVIS,XLIVIS,BLIVIS)                  2596.   
CR    PLINIR=PFNFG(SRTNIR,SRXNIR,RASNIR,XLINIR,BLINIR)                  2597.   
CR    FNSRLI=(FSRVIS*PLIVIS/PNFVIS+FSRNIR*PLINIR/PNFNIR)                2598.   
C                                                                       2599.   
CR604 FNORM=FNSROC*POCEAN+FNSREA*PEARTH+FNSROI*POICE+FNSRLI*PLICE       2600.   
C                                                                       2601.   
CR    FSRNFG(1)=FNSROC/FNORM                                            2602.   
CR    FSRNFG(2)=FNSREA/FNORM                                            2603.   
CR    FSRNFG(3)=FNSROI/FNORM                                            2604.   
CR    FSRNFG(4)=FNSRLI/FNORM                                            2605.   
C                                                                       2606.   
      RETURN                                                            2607.   
      END                                                               2608.   
      SUBROUTINE SETAO2(O2CMA,NL)                                       2609.   
      DIMENSION O2CMA(40),O2FHRL(40)                                    2610.   
      DIMENSION SFWM2(18),SIGMA(18,6)                                   2611.   
      DATA SFWM2/                                                       2612.   
     A 2.196E-03, 0.817E-03, 1.163E-03, 1.331E-03, 1.735E-03, 1.310E-03,2613.   
     B 1.311E-03, 2.584E-03, 2.864E-03, 4.162E-03, 5.044E-03, 6.922E-03,2614.   
     C 6.906E-03,10.454E-03, 5.710E-03, 6.910E-03,14.130E-03,18.080E-03/2615.   
      DATA SIGMA/                                                       2616.   
     A     2.74E-19, 2.74E-19, 2.74E-19, 2.74E-19, 2.74E-19, 2.74E-19,  2617.   
     B     4.33E-21, 4.89E-21, 6.63E-21, 1.60E-20, 7.20E-20, 1.59E-18,  2618.   
     C     2.10E-21, 2.32E-21, 3.02E-21, 6.30E-21, 3.46E-20, 7.52E-19,  2619.   
     D     5.95E-22, 9.72E-22, 2.53E-21, 7.57E-21, 7.38E-20, 7.44E-19,  2620.   
     E     3.33E-22, 1.02E-22, 4.09E-21, 1.63E-20, 8.79E-20, 3.81E-19,  2621.   
     F     1.09E-21, 1.16E-21, 1.45E-21, 3.32E-21, 2.00E-20, 4.04E-19,  2622.   
     G     1.15E-21, 1.30E-21, 1.90E-21, 4.89E-21, 2.62E-20, 4.08E-19,  2623.   
     H     3.90E-22, 4.90E-22, 9.49E-22, 3.33E-21, 2.14E-20, 2.39E-19,  2624.   
     I     1.29E-22, 2.18E-22, 8.28E-22, 3.46E-21, 1.94E-20, 1.06E-19,  2625.   
     J     6.26E-23, 7.80E-23, 2.62E-22, 1.83E-21, 1.25E-20, 3.95E-20,  2626.   
     K     2.74E-23, 3.58E-23, 8.64E-23, 4.03E-22, 2.13E-21, 1.95E-20,  2627.   
     L     1.95E-23, 2.44E-23, 4.89E-23, 2.87E-22, 1.95E-21, 1.36E-20,  2628.   
     M     1.84E-23, 1.96E-23, 2.71E-23, 8.52E-23, 6.48E-22, 3.89E-21,  2629.   
     N     1.80E-23, 1.81E-23, 1.87E-23, 2.69E-23, 1.34E-22, 1.52E-21,  2630.   
     O     1.80E-23, 1.80E-23, 1.82E-23, 2.40E-23, 5.71E-23, 5.70E-22,  2631.   
     P     1.76E-23, 1.76E-23, 1.76E-23, 1.76E-23, 1.76E-23, 3.50E-23,  2632.   
     Q     1.71E-23, 1.71E-23, 1.71E-23, 1.71E-23, 1.71E-23, 2.68E-23,  2633.   
     R     1.00E-23, 1.00E-23, 1.00E-23, 1.00E-23, 1.00E-23, 1.00E-23/  2634.   
      REAL WTKO2(6)/0.05,0.20,0.25,0.25,0.20,0.05/                      2635.   
C                                                                       2636.   
      DATA STPMOL/2.68714E+19/,S00/1367.0/                              2637.   
      DATA NW/18/,NZ/11/,NKO2/6/                                        2638.   
      DIMENSION ZTABLE(40,11)                                           2639.   
      DIMENSION ZCOSJ(11)                                               2640.   
      NLP=NL+1                                                          2641.   
      FSUM=0.0                                                          2642.   
      DO 100 I=1,NW                                                     2643.   
  100 FSUM=FSUM+SFWM2(I)                                                2644.   
      DO 110 J=1,NZ                                                     2645.   
  110 ZTABLE(NLP,J)=FSUM                                                2646.   
      SUMMOL=0.0                                                        2647.   
      DO 150 N=1,NL                                                     2648.   
      L=NLP-N                                                           2649.   
      SUMMOL=SUMMOL+O2CMA(L)*STPMOL                                     2650.   
      DO 140 J=1,NZ                                                     2651.   
      ZCOS=0.01*(1/J)+0.1*(J-1)                                         2652.   
      ZCOSJ(J)=ZCOS                                                     2653.   
      FSUM=0.0                                                          2654.   
      DO 130 I=1,NW                                                     2655.   
      WSUM=0.0                                                          2656.   
      DO 120 K=1,NKO2                                                   2657.   
      TAU=SIGMA(I,K)*SUMMOL/ZCOS                                        2658.   
      IF(TAU.GT.30.0) TAU=30.0                                          2659.   
  120 WSUM=WSUM+WTKO2(K)*EXP(-TAU)                                      2660.   
  130 FSUM=FSUM+WSUM*SFWM2(I)                                           2661.   
  140 ZTABLE(L,J)=FSUM                                                  2662.   
  150 CONTINUE                                                          2663.   
      DO 170 J=1,NZ                                                     2664.   
      DO 160 L=1,NL                                                     2665.   
  160 ZTABLE(L,J)=ZTABLE(L+1,J)-ZTABLE(L,J)                             2666.   
  170 CONTINUE                                                          2667.   
      RETURN                                                            2668.   
C                                                                       2669.   
C---------------------------------                                      2670.   
      ENTRY O2HEAT(O2FHRL,COSZ,S0)                                      2671.   
C---------------------------------                                      2672.   
C                                                                       2673.   
      ZCOS=1.0+10.0*COSZ                                                2674.   
      JI=ZCOS                                                           2675.   
      IF(JI.GT.10) JI=10                                                2676.   
      JJ=JI+1                                                           2677.   
      WTJ=ZCOS-JI                                                       2678.   
      WTI=1.0-WTJ                                                       2679.   
      DO 200 L=1,NLP-1                                                  2680.   
  200 O2FHRL(L)=(WTI*ZTABLE(L,JI)+WTJ*ZTABLE(L,JJ))*S0/S00              2681.   
      RETURN                                                            2682.   
      END                                                               2683.   
      FUNCTION O3ABS(OCM)                                               2684.   
c     DOUBLE PRECISION O3UVAB                                           2684.1  
      DIMENSION  AO3(460)                                               2685.   
C                                                                       2686.   
      IP=0                                                              2687.   
      XX=OCM*1.E+04                                                     2688.   
      IX=XX                                                             2689.   
      IF(IX.GT.99) GO TO 110                                            2690.   
      IF(IX.LT.1 ) GO TO 130                                            2691.   
      GO TO 120                                                         2692.   
  110 IP=IP+90                                                          2693.   
      XX=XX*0.1                                                         2694.   
      IX=XX                                                             2695.   
      IF(IX.GT.99) GO TO 110                                            2696.   
  120 DX=XX-IX                                                          2697.   
      IX=IX+IP                                                          2698.   
      O3ABS=AO3(IX)+DX*(AO3(IX+1)-AO3(IX))                              2699.   
      RETURN                                                            2700.   
  130 O3ABS=XX*AO3(1)                                                   2701.   
      RETURN                                                            2702.   
C                                                                       2703.   
C----------------------                                                 2704.   
      ENTRY SETAO3(OCM)                                                 2705.   
C----------------------                                                 2706.   
C                                                                       2707.   
!     print *,'After 2707'
      DO 140 I=1,460                                                    2708.   
      II=(I-10)/90-4                                                    2709.   
      XX=I-((I-10)/90)*90                                               2710.   
!     print *,i,ii,xx
!     OCM=XX*10.**II                                                    2711.   
!     05/14/2006
      OCM=XX*10.**float(II)
!     print *,ocm
!     05/14/2006
  140 AO3(I)=O3UVAB(OCM)                                                2712.   
!     print *,'After 2712'
      O3ABS=1.                                                          2713.   
      RETURN                                                            2714.   
      END                                                               2715.   
      FUNCTION O3UVAB(OCM)                                              2716.   
c     IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               2717.   
c     REAL OCM                                                          2718.   
C-----------------------------------------------------------------------2719.   
C**** OZONE ABSORPTION COEFFICIENT DATA FROM HANDBOOK OF GEOPHYSICS 19612720.   
C****                   T = -44 DEG CENTR.                              2721.   
C-----------------------------------------------------------------------2722.   
      DIMENSION X(226),F(226)                                           2723.   
      DIMENSION OWMUV2(115),OWMUV3(111),OKEUV2(115),OKEUV3(111)         2724.   
      EQUIVALENCE (X(1),OWMUV2(1)),(X(116),OWMUV3(1)),                  2725.   
     *(F(1),OKEUV2(1)),(F(116),OKEUV3(1))                               2726.   
      DATA OWMUV2/.2002,.2012,.2022,.2032,.2042,.2052,.2062,.2072,.2082,2727.   
     $.2092,.2102,.2112,.2122,.2132,.2142,.2152,.2162,.2172,.2182,.2192,2728.   
     $.2202,.2212,.2222,.2232,.2242,.2252,.2262,.2272,.2282,.2292,.2302,2729.   
     $.2312,.2322,.2332,.2342,.2352,.2362,.2372,.2382,.2392,.2400,.2402,2730.   
     $.2412,.2422,.2432,.2438,.2444,.2452,.2458,.2463,.2472,.2478,.2482,2731.   
     $.2490,.2492,.2500,.2508,.2519,.2527,.2539,.2543,.2553,.2562,.2566,2732.   
     $.2571,.2575,.2579,.2587,.2597,.2604,.2617,.2624,.2635,.2643,.2650,2733.   
     $.2654,.2662,.2669,.2675,.2682,.2692,.2695,.2702,.2712,.2718,.2722,2734.   
     $.2732,.2742,.2746,.2752,.2762,.2772,.2782,.2792,.2802,.2812,.2822,2735.   
     $.2830,.2842,.2852,.2862,.2872,.2882,.2892,.2902,.2912,.2922,.2932,2736.   
     $.2942,.2952,.2962,.2972,.2982,.2992,.2998/                        2737.   
      DATA OWMUV3/.3004,.3016,.3021,.3029,.3036,.3037,.3051,.3053,.3059,2738.   
     $.3061,.3066,.3075,.3077,.3083,.3085,.3092,.3098,.3100,.3104,.3106,2739.   
     $.3109,.3112,.3130,.3135,.3146,.3148,.3151,.3154,.3167,.3170,.3173,2740.   
     $.3176,.3190,.3194,.3199,.3200,.3209,.3210,.3216,.3220,.3223,.3226,2741.   
     $.3239,.3242,.3245,.3248,.3253,.3255,.3269,.3272,.3275,.3279,.3292,2742.   
     $.3295,.3299,.3303,.3309,.3312,.3328,.3332,.3334,.3338,.3357,.3365,2743.   
     $.3369,.3372,.3391,.3395,.3398,.3401,.3417,.3421,.3426,.3430,.3437,2744.   
     $.3439,.3451,.3455,.3460,.3463,.3466,.3472,.3481,.3485,.3489,.3493,2745.   
     $.3499,.3501,.3506,.3514,.3521,.3523,.3546,.3550,.3554,.3556,.3561,2746.   
     $.3567,.3572,.3573,.3588,.3594,.3599,.3600,.3604,.3606,.3639,.3647,2747.   
     $.3650,.3654,.3660/                                                2748.   
      DATA OKEUV2/  8.3,  8.3,  8.1,  8.3,  8.6,  9.0,  9.7, 10.8, 11.7,2749.   
     $ 13.0, 14.3, 16.0, 18.0, 20.6, 23.0, 26.1, 29.3, 32.6, 36.9, 40.8,2750.   
     $ 46.9, 51.4, 56.7, 63.4, 69.1, 76.6, 84.0, 91.4, 99.9,110.0,118.0,2751.   
     $126.0,136.0,145.0,154.0,164.0,175.0,186.0,192.0,201.0,210.0,212.0,2752.   
     $221.0,230.0,239.0,248.0,250.0,259.0,264.0,264.0,273.0,277.0,275.0,2753.   
     $283.0,283.0,290.0,283.0,297.0,290.0,300.0,290.0,302.0,295.0,283.0,2754.   
     $293.0,290.0,286.0,297.0,281.0,280.0,271.0,275.0,254.0,264.0,250.0,2755.   
     $248.0,242.0,228.0,230.0,216.0,213.0,211.0,199.0,188.0,188.0,178.0,2756.   
     $169.0,153.0,155.0,148.0,136.0,127.0,117.0,108.0, 97.0, 88.7, 81.3,2757.   
     $ 78.7, 67.9, 61.4, 54.3, 49.6, 43.1, 38.9, 34.6, 30.2, 27.5, 23.9,2758.   
     $ 21.0, 18.6, 16.2, 14.2, 12.3, 10.7,  9.5/                        2759.   
      DATA OKEUV3/8.880,7.520,6.960,6.160,5.810,5.910,4.310,4.430,4.130,2760.   
     $4.310,4.020,3.330,3.390,3.060,3.100,2.830,2.400,2.490,2.330,2.320,2761.   
     $2.120,2.200,1.436,1.595,1.074,1.138,1.068,1.262,0.818,0.948,0.860,2762.   
     $1.001,0.543,0.763,0.665,0.781,0.382,0.406,0.373,0.608,0.484,0.601,2763.   
     $0.209,0.276,0.259,0.470,0.319,0.354,0.131,0.223,0.185,0.339,0.080,2764.   
     $0.093,0.079,0.184,0.139,0.214,0.053,0.074,0.068,0.152,0.038,0.070,2765.   
     $.0540000,.1030000,.0240000,.0382500,.0292500,.0550000,.0135000,   2766.   
     $.0155250,.0127500,.0188250,.0167250,.0262500,.0115500,.0140250,   2767.   
     $.0099750,.0115500,.0081000,.0104250,.0050100,.0057000,.0046650,   2768.   
     $.0073425,.0051825,.0055275,.0040575,.0077700,.0048900,.0054600,   2769.   
     $.0015375,.0017775,.0013275,.0014100,.0011550,.0023325,.0018825,   2770.   
     $.0019650,.0009600,.0013650,.0011925,.0013200,.0008925,.0009825,   2771.   
     $.0001350,.0006300,.0004500,.0006225,0.0/                          2772.   
C                                                                       2773.   
C        THEKAERAKA SOLAR FLUX                                          2774.   
C                                                                       2775.   
      DIMENSION Y(190),H(190)                                           2776.   
      DATA     H/.007,.900,.007,.007,.030,.070,.230,.630,1.25,2.71,10.7,2777.   
     1 22.9,57.5,64.9,66.7,59.3,63.0,72.3,70.4,104.,130.,185.,232.,204.,2778.   
     2    222.,315.,482.,584.,514.,603.,689.,764.,830.,975.,1059.,1081.,2779.   
     31074.,1069.,1093.,1083.,1068.,1132.,1181.,1157.,1120.,1098.,1098.,2780.   
     41189.,1429.,1644.,1751.,1774.,1747.,1693.,1639.,1663.,1810.,1922.,2781.   
     52006.,2057.,2066.,2048.,2033.,2044.,2074.,1976.,1950.,1960.,1942.,2782.   
     61920.,1882.,1833.,1833.,1852.,1842.,1818.,1783.,1754.,1725.,1720.,2783.   
     71695.,1705.,1712.,1719.,1715.,1712.,1700.,1682.,1666.,1647.,1635.,2784.   
     81602.,1570.,1544.,1511.,1486.,1456.,1427.,1402.,1389.,1344.,1314.,2785.   
     91290.,1260.,1235.,1211.,1185.,1159.,1134.,1109.,1085.,1060.,1036.,2786.   
     A1013.,990.,968.,947.,926.,908.,891.,880.,869.,858.,847.,837.,820.,2787.   
     B 803.,785.,767.,748.,668.,593.,535.,485.,438.,397.,358.,337.,312.,2788.   
     C 288.,267.,245.,223.,202.,180.,159.,142.,126.,114.,103., 90., 79.,2789.   
     D 69.0,62.0,55.0,48.0,43.0,39.0,35.0,31.0,26.0,22.6,19.2,16.6,14.6,2790.   
     E 13.5,12.3,11.1,10.3, 9.5,8.70,7.80,7.10,6.50,5.92,5.35,4.86,4.47,2791.   
     F   4.11,3.79,1.82,0.99,.585,.367,.241,.165,.117,.0851,.0634,.0481/2792.   
      DATA     Y/.115,.120,.125,.130,.140,.150,.160,.170,.180,.190,.200,2793.   
     1 .210,.220,.225,.230,.235,.240,.245,.250,.255,.260,.265,.270,.275,2794.   
     2      .280,.285,.290,.295,.300,.305,.310,.315,.320,.325,.330,.335,2795.   
     3           .340,.345,.350,.355,.360,.365,.370,.375,.380,.385,.390,2796.   
     4           .395,.400,.405,.410,.415,.420,.425,.430,.435,.440,.445,2797.   
     5           .450,.455,.460,.465,.470,.475,.480,.485,.490,.495,.500,2798.   
     6           .505,.510,.515,.520,.525,.530,.535,.540,.545,.550,.555,2799.   
     7           .560,.565,.570,.575,.580,.585,.590,.595,.600,.605,.610,2800.   
     8           .620,.630,.640,.650,.660,.670,.680,.690,.700,.710,.720,2801.   
     9           .730,.740,.750,.760,.770,.780,.790,.800,.810,.820,.830,2802.   
     A .840,.850,.860,.870,.880,.890,.900,.910,.920,.930,.940,.950,.960,2803.   
     B 0.97,0.98,0.99,1.00,1.05,1.10,1.15,1.20,1.25,1.30,1.35,1.40,1.45,2804.   
     C 1.50,1.55,1.60,1.65,1.70,1.75,1.80,1.85,1.90,1.95,2.00,2.10,2.20,2805.   
     D 2.30,2.40,2.50,2.60,2.70,2.80,2.90,3.00,3.10,3.20,3.30,3.40,3.50,2806.   
     E 3.60,3.70,3.80,3.90,4.00,4.10,4.20,4.30,4.40,4.50,4.60,4.70,4.80,2807.   
     F      4.9, 5.0, 6.0, 7.0, 8.0, 9.0,10.0,11.0,12.0,13.0,14.0,15.00/2808.   
      NH=190                                                            2809.   
      NG=226                                                            2810.   
      XA=X(1)                                                           2811.   
      XB=X(NG)                                                          2812.   
      SOLCON=0.1353D0                                                   2813.   
      ABINT=0.D0                                                        2814.   
      X2=DMIN1(X(NG),Y(NH))                                             2815.   
      IF(XA.GE.X2) GO TO 160                                            2816.   
      X1=DMAX1(X(1),Y(1))                                               2817.   
      IF(XB.LE.X1) GO TO 160                                            2818.   
      YA=XA                                                             2819.   
      IF(XA.LT.X1) YA=X1                                                2820.   
      YB=XB                                                             2821.   
      IF(YB.GT.X2) YB=X2                                                2822.   
      DO 100 JG=2,NG                                                    2823.   
      XJ=X(JG)                                                          2824.   
      IF(XJ.GT.YA) GO TO 110                                            2825.   
100   CONTINUE                                                          2825.1  
      JG=NG+1                                                           2825.2  
110   IG=JG-1                                                           2826.   
      XI=X(IG)                                                          2827.   
      TAU=F(IG)*OCM                                                     2828.   
      IF(TAU.GT.35.D0) TAU=35.D0                                        2829.   
      GI=1.D0-DEXP(-TAU)                                                2830.   
      TAU=F(JG)*OCM                                                     2831.   
      IF(TAU.GT.35.D0) TAU=35.D0                                        2832.   
      GJ=1.D0-DEXP(-TAU)                                                2833.   
      B=(GJ-GI)/(XJ-XI)                                                 2834.   
      A=GJ-B*XJ                                                         2835.   
      DO 120 JH=2,NH                                                    2836.   
      YJ=Y(JH)                                                          2837.   
      IF(YJ.GT.YA) GO TO 130                                            2838.   
120   CONTINUE                                                          2838.1  
      JH=NH+1                                                           2838.2  
130   IH=JH-1                                                           2839.   
      YI=Y(IH)                                                          2840.   
      HI=H(IH)/10000.D0                                                 2841.   
      HJ=H(JH)/10000.D0                                                 2842.   
      D=(HJ-HI)/(YJ-YI)                                                 2843.   
      C=HJ-D*YJ                                                         2844.   
      X2=YA                                                             2845.   
140   X1=X2                                                             2846.   
      X2=DMIN1(XJ,YJ)                                                   2847.   
      DELTA=(XJ-YJ)/(XJ+YJ)                                             2848.   
      IF(X2.GT.YB) X2=YB                                                2849.   
      DINT=(X2-X1)*(A*C+0.5D0*(B*C+A*D)*(X2+X1)+B*D*(X2*(X2+X1)+X1*X1)/ 2850.   
     $3.D0)                                                             2851.   
      ABINT=ABINT+DINT                                                  2852.   
      IF(X2.GE.YB) GO TO 160                                            2853.   
      IF(DELTA.GT.1.D-14) GO TO 150                                     2854.   
      XI=XJ                                                             2855.   
      GI=GJ                                                             2856.   
      JG=JG+1                                                           2857.   
      XJ=X(JG)                                                          2858.   
      TAU=F(JG)*OCM                                                     2859.   
      IF(TAU.GT.35.D0) TAU=35.D0                                        2860.   
      GJ=1.D0-DEXP(-TAU)                                                2861.   
      B=(GJ-GI)/(XJ-XI)                                                 2862.   
      A=GJ-B*XJ                                                         2863.   
      IF(DABS(DELTA).LE.1.D-14) GO TO 150                               2864.   
      GO TO 140                                                         2865.   
150   YI=YJ                                                             2866.   
      HI=HJ                                                             2867.   
      JH=JH+1                                                           2868.   
      YJ=Y(JH)                                                          2869.   
      HJ=H(JH)/10000.D0                                                 2870.   
      D=(HJ-HI)/(YJ-YI)                                                 2871.   
      C=HJ-D*YJ                                                         2872.   
      GO TO 140                                                         2873.   
160   O3UVAB=ABINT/SOLCON                                               2874.   
      RETURN                                                            2875.   
      END                                                               2876.   
      SUBROUTINE SETO3D                                                 2877.   

#include "B83XX.COM"  

C-----------------------------------------------------------------------2915.   
C                                                                       2916.   
C        LONDON ET AL (1976) JUL,1957-DEC,1970 NCAR ATLAS OF TOTAL OZONE2917.   
C                                                                       2918.   
C        AVERAGE GLOBAL COLUMN AMOUNT -- O3AVE(MONTH,LATITUDE,LONGITUDE)2919.   
C                                                                       2920.   
C                                  MONTH=1-12  JAN,FEB,...,DEC          2921.   
C                                  LAT  =1-18  -85,-75,..., 85          2922.   
C                                                                       2923.   
C-----------------------------------------------------------------------2924.   
      REAL O3AVEA(216),O3AVEB(216),O3AVEC(216),O3AVED(216),O3AVEE(216)  2925.   
      REAL O3AVEF(216),O3AVEG(216),O3AVEH(216),O3AVEI(216),O3AVEJ(216)  2926.   
      REAL O3AVEK(216),O3AVEL(216),O3AVEM(216),O3AVEN(216),O3AVEO(216)  2927.   
      REAL O3AVEP(216),O3AVEQ(216),O3AVER(216),O3AVE(12,18,18)          2928.   
      EQUIVALENCE (O3AVE(1,1,10),O3AVEA(1)),(O3AVE(1,1,11),O3AVEB(1))   2929.   
     1           ,(O3AVE(1,1,12),O3AVEC(1)),(O3AVE(1,1,13),O3AVED(1))   2930.   
     2           ,(O3AVE(1,1,14),O3AVEE(1)),(O3AVE(1,1,15),O3AVEF(1))   2931.   
     3           ,(O3AVE(1,1,16),O3AVEG(1)),(O3AVE(1,1,17),O3AVEH(1))   2932.   
     4           ,(O3AVE(1,1,18),O3AVEI(1)),(O3AVE(1,1,01),O3AVEJ(1))   2933.   
     5           ,(O3AVE(1,1,02),O3AVEK(1)),(O3AVE(1,1,03),O3AVEL(1))   2934.   
     6           ,(O3AVE(1,1,04),O3AVEM(1)),(O3AVE(1,1,05),O3AVEN(1))   2935.   
     7           ,(O3AVE(1,1,06),O3AVEO(1)),(O3AVE(1,1,07),O3AVEP(1))   2936.   
     8           ,(O3AVE(1,1,08),O3AVEQ(1)),(O3AVE(1,1,09),O3AVER(1))   2937.   
      DATA O3AVEA/                                                      2938.   
     A     .317,.295,.291,.292,.293,.298,.300,.305,.313,.324,.369,.355, 2939.   
     B     .319,.300,.296,.292,.291,.300,.301,.304,.314,.322,.358,.350, 2940.   
     C     .312,.301,.295,.287,.286,.298,.302,.305,.316,.322,.343,.335, 2941.   
     D     .299,.291,.285,.280,.279,.290,.295,.300,.307,.319,.327,.316, 2942.   
     E     .281,.275,.279,.268,.266,.278,.282,.290,.295,.306,.306,.296, 2943.   
     F     .266,.261,.259,.256,.252,.261,.267,.277,.280,.289,.285,.277, 2944.   
     G     .252,.249,.248,.246,.240,.249,.252,.262,.264,.273,.265,.258, 2945.   
     H     .240,.238,.240,.242,.237,.242,.240,.249,.252,.258,.251,.245, 2946.   
     I     .232,.230,.238,.241,.240,.238,.234,.241,.241,.245,.239,.236, 2947.   
     J     .235,.235,.244,.252,.253,.244,.236,.237,.232,.230,.230,.232, 2948.   
     K     .249,.256,.264,.269,.267,.261,.245,.245,.238,.234,.233,.237, 2949.   
     L     .278,.289,.294,.300,.294,.284,.265,.265,.256,.249,.248,.261, 2950.   
     M     .318,.338,.343,.351,.342,.324,.300,.296,.287,.275,.279,.299, 2951.   
     N     .347,.368,.383,.383,.370,.351,.335,.319,.304,.288,.296,.321, 2952.   
     O     .364,.394,.418,.410,.402,.371,.358,.340,.312,.298,.302,.325, 2953.   
     P     .356,.388,.421,.414,.394,.360,.337,.319,.299,.285,.292,.313, 2954.   
     Q     .364,.403,.431,.426,.398,.358,.328,.303,.292,.287,.297,.324, 2955.   
     R     .373,.421,.447,.440,.408,.355,.323,.295,.289,.291,.305,.329/ 2956.   
      DATA O3AVEB/                                                      2957.   
     A     .318,.295,.291,.293,.293,.299,.301,.305,.314,.326,.372,.358, 2958.   
     B     .321,.300,.295,.293,.291,.301,.301,.306,.314,.326,.361,.353, 2959.   
     C     .315,.302,.296,.291,.288,.300,.303,.306,.318,.328,.348,.340, 2960.   
     D     .307,.296,.291,.284,.278,.298,.299,.305,.314,.326,.335,.324, 2961.   
     E     .294,.285,.286,.272,.270,.286,.288,.296,.302,.315,.315,.304, 2962.   
     F     .278,.271,.265,.260,.258,.270,.273,.283,.287,.298,.293,.284, 2963.   
     G     .262,.259,.254,.250,.247,.255,.259,.268,.270,.282,.274,.266, 2964.   
     H     .247,.246,.244,.245,.239,.245,.247,.255,.255,.266,.257,.250, 2965.   
     I     .235,.235,.239,.244,.240,.238,.236,.244,.244,.249,.244,.239, 2966.   
     J     .233,.234,.243,.251,.249,.240,.234,.235,.232,.231,.231,.231, 2967.   
     K     .247,.254,.263,.267,.262,.253,.242,.240,.237,.232,.232,.237, 2968.   
     L     .279,.287,.296,.282,.286,.275,.260,.257,.253,.246,.246,.258, 2969.   
     M     .320,.336,.345,.348,.325,.309,.293,.282,.279,.267,.272,.294, 2970.   
     N     .346,.369,.379,.377,.348,.330,.317,.299,.286,.280,.288,.312, 2971.   
     O     .368,.406,.412,.401,.373,.345,.332,.312,.293,.284,.293,.316, 2972.   
     P     .366,.409,.423,.418,.386,.349,.326,.307,.290,.278,.295,.312, 2973.   
     Q     .366,.407,.428,.429,.396,.352,.323,.296,.287,.282,.298,.318, 2974.   
     R     .372,.420,.446,.441,.407,.352,.320,.292,.286,.290,.305,.327/ 2975.   
      DATA O3AVEC/                                                      2976.   
     A     .319,.296,.292,.294,.294,.299,.302,.306,.316,.328,.372,.359, 2977.   
     B     .321,.300,.295,.297,.293,.303,.305,.309,.319,.332,.367,.359, 2978.   
     C     .322,.309,.302,.297,.293,.309,.309,.314,.326,.338,.362,.353, 2979.   
     D     .324,.313,.303,.294,.295,.314,.311,.318,.330,.342,.353,.343, 2980.   
     E     .315,.308,.296,.286,.287,.305,.306,.314,.326,.335,.338,.326, 2981.   
     F     .294,.290,.281,.271,.273,.287,.290,.299,.307,.319,.312,.303, 2982.   
     G     .274,.272,.264,.258,.258,.268,.272,.281,.286,.297,.290,.281, 2983.   
     H     .254,.254,.251,.248,.248,.254,.257,.263,.267,.276,.271,.262, 2984.   
     I     .240,.239,.241,.245,.241,.243,.244,.250,.251,.256,.250,.246, 2985.   
     J     .230,.231,.238,.249,.246,.237,.234,.233,.234,.233,.230,.228, 2986.   
     K     .238,.244,.251,.258,.253,.244,.236,.235,.233,.228,.228,.230, 2987.   
     L     .259,.269,.276,.279,.268,.254,.246,.241,.238,.235,.237,.246, 2988.   
     M     .289,.305,.312,.306,.289,.270,.261,.255,.249,.246,.252,.268, 2989.   
     N     .321,.347,.354,.343,.315,.291,.281,.273,.262,.259,.268,.285, 2990.   
     O     .351,.394,.396,.384,.353,.315,.300,.288,.275,.271,.282,.296, 2991.   
     P     .363,.414,.422,.415,.382,.333,.313,.292,.281,.276,.292,.306, 2992.   
     Q     .366,.415,.430,.433,.398,.346,.313,.288,.282,.280,.299,.317, 2993.   
     R     .372,.421,.445,.441,.406,.348,.316,.289,.285,.289,.306,.327/ 2994.   
      DATA O3AVED/                                                      2995.   
     A     .320,.296,.293,.294,.295,.300,.303,.308,.317,.330,.374,.361, 2996.   
     B     .322,.300,.297,.299,.296,.307,.310,.314,.323,.339,.373,.366, 2997.   
     C     .329,.313,.310,.304,.302,.320,.318,.326,.338,.352,.373,.367, 2998.   
     D     .343,.330,.318,.306,.315,.333,.329,.337,.354,.366,.370,.366, 2999.   
     E     .334,.324,.311,.299,.312,.326,.329,.333,.352,.357,.354,.342, 3000.   
     F     .304,.300,.291,.279,.285,.302,.308,.315,.324,.328,.325,.312, 3001.   
     G     .277,.276,.268,.262,.266,.279,.283,.289,.296,.303,.299,.283, 3002.   
     H     .256,.257,.253,.249,.252,.259,.266,.269,.274,.278,.273,.263, 3003.   
     I     .242,.243,.243,.248,.247,.251,.255,.256,.258,.260,.253,.249, 3004.   
     J     .231,.234,.238,.250,.255,.251,.250,.246,.248,.244,.237,.229, 3005.   
     K     .235,.241,.248,.257,.259,.257,.248,.246,.245,.244,.233,.230, 3006.   
     L     .256,.261,.267,.270,.269,.262,.251,.247,.247,.248,.239,.248, 3007.   
     M     .293,.304,.306,.302,.288,.272,.259,.256,.256,.256,.254,.269, 3008.   
     N     .327,.344,.356,.346,.319,.291,.272,.270,.264,.267,.270,.285, 3009.   
     O     .356,.392,.402,.388,.359,.312,.289,.281,.276,.281,.285,.297, 3010.   
     P     .368,.416,.424,.415,.388,.328,.304,.285,.279,.284,.295,.309, 3011.   
     Q     .370,.418,.436,.436,.402,.338,.306,.283,.278,.284,.301,.320, 3012.   
     R     .373,.422,.446,.441,.407,.345,.312,.286,.275,.291,.307,.328/ 3013.   
      DATA O3AVEE/                                                      3014.   
     A     .319,.295,.293,.295,.296,.300,.304,.309,.318,.332,.375,.362, 3015.   
     B     .325,.301,.300,.302,.300,.309,.313,.319,.328,.345,.378,.370, 3016.   
     C     .332,.314,.312,.310,.310,.327,.329,.335,.347,.362,.381,.375, 3017.   
     D     .348,.334,.324,.312,.328,.346,.366,.352,.372,.381,.377,.373, 3018.   
     E     .337,.327,.318,.303,.322,.335,.342,.347,.363,.366,.358,.344, 3019.   
     F     .301,.297,.292,.282,.291,.307,.314,.321,.331,.332,.324,.309, 3020.   
     G     .275,.271,.269,.264,.270,.279,.286,.292,.299,.301,.293,.281, 3021.   
     H     .255,.253,.252,.251,.253,.258,.265,.269,.275,.277,.268,.262, 3022.   
     I     .245,.244,.246,.250,.249,.253,.254,.257,.259,.260,.252,.249, 3023.   
     J     .240,.239,.245,.255,.256,.260,.256,.253,.253,.251,.243,.237, 3024.   
     K     .247,.248,.252,.263,.270,.268,.258,.256,.256,.252,.244,.238, 3025.   
     L     .263,.263,.268,.277,.282,.276,.261,.259,.259,.258,.251,.251, 3026.   
     M     .299,.304,.309,.309,.302,.291,.269,.266,.268,.269,.269,.275, 3027.   
     N     .346,.358,.365,.353,.335,.307,.276,.272,.276,.283,.289,.300, 3028.   
     O     .379,.400,.414,.401,.373,.319,.286,.280,.283,.293,.303,.314, 3029.   
     P     .382,.421,.437,.427,.398,.323,.293,.280,.280,.293,.308,.321, 3030.   
     Q     .375,.424,.444,.440,.405,.334,.298,.278,.276,.290,.306,.326, 3031.   
     R     .374,.424,.448,.443,.406,.345,.310,.284,.281,.292,.309,.328/ 3032.   
      DATA O3AVEF/                                                      3033.   
     A     .318,.294,.294,.295,.298,.301,.304,.311,.320,.333,.377,.361, 3034.   
     B     .324,.298,.300,.304,.305,.310,.315,.323,.331,.348,.383,.371, 3035.   
     C     .337,.311,.314,.313,.317,.330,.333,.344,.354,.369,.386,.377, 3036.   
     D     .350,.330,.324,.317,.332,.349,.351,.362,.378,.390,.380,.372, 3037.   
     E     .333,.322,.314,.307,.323,.339,.345,.358,.369,.372,.357,.340, 3038.   
     F     .300,.292,.286,.284,.294,.307,.316,.327,.335,.334,.323,.307, 3039.   
     G     .275,.269,.264,.263,.269,.277,.285,.292,.300,.303,.290,.279, 3040.   
     H     .254,.251,.250,.251,.254,.256,.261,.267,.271,.276,.266,.261, 3041.   
     I     .243,.242,.242,.247,.248,.250,.247,.251,.252,.258,.253,.247, 3042.   
     J     .237,.239,.243,.253,.255,.255,.246,.243,.244,.245,.239,.236, 3043.   
     K     .246,.247,.253,.263,.265,.265,.253,.245,.247,.247,.239,.238, 3044.   
     L     .265,.265,.276,.283,.284,.280,.261,.254,.253,.258,.250,.250, 3045.   
     M     .306,.309,.321,.316,.318,.292,.273,.259,.265,.271,.273,.277, 3046.   
     N     .365,.369,.381,.363,.347,.313,.278,.264,.275,.290,.302,.307, 3047.   
     O     .396,.416,.431,.415,.405,.322,.282,.271,.288,.303,.321,.328, 3048.   
     P     .397,.433,.455,.436,.404,.322,.287,.273,.276,.302,.320,.333, 3049.   
     Q     .382,.429,.451,.442,.408,.331,.297,.274,.273,.295,.311,.333, 3050.   
     R     .375,.427,.450,.445,.407,.343,.309,.283,.280,.295,.311,.330/ 3051.   
      DATA O3AVEG/                                                      3052.   
     A     .317,.293,.293,.295,.299,.299,.305,.311,.320,.335,.378,.360, 3053.   
     B     .323,.296,.300,.304,.306,.310,.317,.325,.334,.353,.385,.367, 3054.   
     C     .335,.307,.310,.312,.318,.328,.335,.347,.357,.376,.390,.372, 3055.   
     D     .346,.324,.320,.317,.332,.349,.354,.367,.384,.393,.384,.368, 3056.   
     E     .331,.318,.311,.305,.324,.339,.349,.365,.378,.377,.360,.339, 3057.   
     F     .301,.293,.286,.285,.296,.309,.321,.334,.344,.339,.325,.309, 3058.   
     G     .276,.270,.266,.267,.271,.280,.287,.295,.303,.308,.294,.282, 3059.   
     H     .257,.253,.250,.252,.254,.257,.261,.266,.271,.279,.268,.261, 3060.   
     I     .240,.241,.241,.246,.246,.250,.246,.249,.253,.259,.254,.248, 3061.   
     J     .234,.238,.245,.256,.258,.259,.244,.243,.241,.243,.237,.235, 3062.   
     K     .244,.249,.259,.271,.274,.274,.257,.251,.248,.248,.238,.237, 3063.   
     L     .270,.272,.289,.297,.298,.294,.277,.267,.260,.262,.251,.254, 3064.   
     M     .329,.338,.353,.338,.333,.313,.296,.275,.273,.282,.281,.296, 3065.   
     N     .401,.414,.424,.392,.369,.329,.298,.272,.282,.303,.321,.341, 3066.   
     O     .420,.451,.461,.432,.389,.331,.291,.272,.279,.313,.343,.358, 3067.   
     P     .411,.451,.468,.447,.403,.320,.289,.271,.277,.308,.334,.349, 3068.   
     Q     .386,.434,.456,.443,.404,.332,.297,.273,.273,.300,.317,.339, 3069.   
     R     .378,.430,.453,.446,.407,.342,.310,.282,.279,.296,.314,.332/ 3070.   
      DATA O3AVEH/                                                      3071.   
     A     .315,.292,.293,.295,.299,.297,.303,.311,.320,.334,.378,.358, 3072.   
     B     .320,.294,.298,.303,.306,.308,.316,.325,.337,.355,.387,.362, 3073.   
     C     .330,.304,.307,.311,.315,.323,.334,.345,.360,.381,.389,.366, 3074.   
     D     .339,.318,.312,.314,.328,.344,.355,.368,.388,.401,.384,.360, 3075.   
     E     .325,.313,.302,.300,.318,.339,.354,.369,.381,.380,.360,.337, 3076.   
     F     .299,.291,.285,.284,.296,.313,.326,.340,.350,.343,.328,.312, 3077.   
     G     .277,.271,.269,.269,.272,.281,.288,.296,.308,.311,.298,.289, 3078.   
     H     .257,.253,.252,.254,.253,.257,.262,.267,.272,.281,.272,.265, 3079.   
     I     .241,.241,.241,.246,.245,.248,.246,.248,.253,.260,.255,.250, 3080.   
     J     .234,.236,.242,.256,.260,.260,.246,.244,.240,.241,.237,.237, 3081.   
     K     .243,.246,.257,.273,.279,.276,.261,.258,.251,.246,.238,.238, 3082.   
     L     .270,.269,.288,.299,.308,.299,.283,.276,.269,.263,.252,.257, 3083.   
     M     .327,.339,.358,.349,.351,.337,.313,.292,.288,.280,.284,.302, 3084.   
     N     .407,.419,.432,.407,.390,.356,.324,.298,.300,.304,.327,.368, 3085.   
     O     .421,.455,.459,.439,.393,.333,.306,.287,.289,.311,.345,.377, 3086.   
     P     .408,.452,.465,.443,.399,.323,.296,.276,.279,.309,.338,.362, 3087.   
     Q     .387,.437,.459,.444,.404,.334,.301,.276,.277,.302,.320,.345, 3088.   
     R     .379,.433,.455,.447,.408,.343,.313,.282,.279,.298,.315,.336/ 3089.   
      DATA O3AVEI/                                                      3090.   
     A     .313,.291,.291,.293,.299,.296,.302,.310,.319,.333,.379,.354, 3091.   
     B     .316,.292,.295,.300,.307,.306,.315,.322,.333,.354,.384,.354, 3092.   
     C     .322,.302,.301,.307,.309,.319,.331,.340,.357,.379,.385,.356, 3093.   
     D     .328,.310,.301,.306,.316,.332,.347,.359,.380,.397,.379,.348, 3094.   
     E     .315,.304,.293,.296,.308,.328,.345,.360,.374,.376,.356,.329, 3095.   
     F     .292,.285,.277,.278,.288,.304,.318,.330,.340,.340,.324,.306, 3096.   
     G     .271,.266,.262,.263,.266,.277,.283,.291,.301,.307,.293,.284, 3097.   
     H     .253,.249,.249,.252,.250,.256,.261,.267,.271,.278,.267,.263, 3098.   
     I     .240,.238,.240,.247,.244,.248,.247,.250,.254,.258,.251,.249, 3099.   
     J     .233,.236,.243,.254,.259,.258,.248,.246,.241,.243,.238,.238, 3100.   
     K     .242,.246,.256,.268,.273,.271,.260,.255,.250,.244,.240,.239, 3101.   
     L     .258,.266,.278,.290,.295,.288,.277,.269,.265,.257,.253,.256, 3102.   
     M     .294,.308,.325,.326,.322,.308,.297,.284,.278,.271,.277,.287, 3103.   
     N     .338,.368,.383,.371,.357,.329,.316,.294,.287,.288,.303,.324, 3104.   
     O     .375,.420,.429,.411,.382,.328,.312,.293,.287,.299,.322,.354, 3105.   
     P     .388,.440,.454,.437,.396,.328,.307,.285,.282,.305,.330,.359, 3106.   
     Q     .386,.439,.457,.444,.404,.338,.309,.283,.280,.304,.321,.349, 3107.   
     R     .379,.435,.456,.448,.408,.345,.316,.286,.281,.300,.317,.337/ 3108.   
      DATA O3AVEJ/                                                      3109.   
     A     .313,.290,.290,.291,.298,.294,.301,.309,.318,.331,.378,.353, 3110.   
     B     .313,.291,.291,.296,.304,.302,.311,.318,.330,.348,.382,.350, 3111.   
     C     .315,.297,.294,.300,.306,.310,.325,.334,.348,.364,.378,.346, 3112.   
     D     .316,.301,.292,.297,.305,.317,.334,.346,.360,.371,.366,.335, 3113.   
     E     .304,.293,.283,.286,.295,.313,.330,.344,.356,.359,.346,.316, 3114.   
     F     .284,.276,.268,.271,.279,.297,.309,.320,.325,.330,.317,.296, 3115.   
     G     .265,.258,.254,.257,.261,.273,.280,.288,.289,.296,.287,.274, 3116.   
     H     .250,.245,.244,.249,.247,.255,.260,.265,.268,.273,.263,.257, 3117.   
     I     .237,.235,.238,.246,.246,.249,.247,.249,.251,.257,.249,.247, 3118.   
     J     .234,.236,.245,.256,.259,.255,.248,.249,.244,.245,.242,.238, 3119.   
     K     .244,.249,.259,.271,.273,.270,.258,.256,.253,.247,.243,.242, 3120.   
     L     .261,.273,.283,.291,.292,.284,.271,.269,.263,.257,.254,.257, 3121.   
     M     .289,.305,.319,.321,.315,.301,.287,.281,.273,.268,.272,.282, 3122.   
     N     .321,.347,.364,.358,.344,.319,.305,.293,.282,.281,.291,.313, 3123.   
     O     .357,.400,.409,.397,.373,.332,.314,.295,.286,.293,.309,.333, 3124.   
     P     .377,.429,.442,.429,.396,.338,.317,.294,.287,.302,.321,.351, 3125.   
     Q     .385,.439,.458,.443,.407,.345,.318,.292,.284,.304,.322,.349, 3126.   
     R     .380,.437,.458,.449,.408,.348,.319,.289,.283,.301,.319,.340/ 3127.   
      DATA O3AVEK/                                                      3128.   
     A     .311,.289,.289,.290,.298,.293,.300,.308,.317,.329,.377,.352, 3129.   
     B     .308,.290,.288,.291,.301,.296,.307,.315,.326,.340,.377,.344, 3130.   
     C     .305,.291,.287,.293,.297,.302,.315,.325,.335,.346,.369,.333, 3131.   
     D     .299,.289,.281,.287,.293,.302,.317,.327,.335,.344,.353,.318, 3132.   
     E     .287,.279,.272,.277,.281,.295,.309,.320,.325,.332,.331,.301, 3133.   
     F     .272,.264,.259,.262,.268,.281,.292,.300,.300,.309,.305,.282, 3134.   
     G     .257,.249,.246,.250,.254,.264,.271,.278,.279,.285,.278,.263, 3135.   
     H     .246,.239,.239,.245,.245,.252,.255,.261,.262,.267,.259,.250, 3136.   
     I     .234,.231,.239,.245,.245,.248,.245,.249,.248,.254,.246,.243, 3137.   
     J     .235,.237,.247,.258,.260,.257,.250,.250,.245,.246,.241,.240, 3138.   
     K     .248,.254,.264,.276,.276,.272,.262,.258,.255,.250,.248,.246, 3139.   
     L     .267,.278,.289,.300,.296,.286,.272,.270,.263,.258,.258,.262, 3140.   
     M     .292,.310,.325,.329,.319,.302,.288,.280,.273,.268,.274,.281, 3141.   
     N     .323,.346,.365,.365,.347,.320,.305,.291,.282,.281,.292,.305, 3142.   
     O     .352,.390,.405,.398,.378,.338,.316,.300,.290,.294,.309,.330, 3143.   
     P     .376,.424,.440,.431,.404,.350,.323,.303,.293,.303,.321,.349, 3144.   
     Q     .386,.442,.462,.448,.411,.354,.324,.298,.289,.306,.325,.349, 3145.   
     R     .381,.441,.459,.452,.410,.352,.322,.293,.286,.301,.320,.342/ 3146.   
      DATA O3AVEL/                                                      3147.   
     A     .309,.290,.288,.288,.295,.292,.299,.307,.315,.327,.375,.350, 3148.   
     B     .306,.289,.287,.288,.298,.293,.304,.311,.320,.333,.372,.340, 3149.   
     C     .298,.286,.282,.288,.290,.294,.308,.316,.322,.332,.362,.325, 3150.   
     D     .289,.280,.274,.281,.282,.290,.304,.312,.317,.325,.342,.309, 3151.   
     E     .276,.269,.264,.268,.271,.281,.293,.300,.304,.313,.318,.290, 3152.   
     F     .262,.256,.253,.255,.258,.267,.278,.283,.283,.293,.294,.272, 3153.   
     G     .250,.245,.241,.245,.246,.255,.261,.267,.265,.282,.272,.256, 3154.   
     H     .240,.235,.236,.243,.240,.245,.249,.254,.253,.260,.254,.247, 3155.   
     I     .232,.229,.239,.245,.244,.247,.241,.245,.241,.246,.243,.241, 3156.   
     J     .235,.236,.247,.258,.258,.254,.246,.246,.239,.240,.238,.240, 3157.   
     K     .248,.253,.263,.273,.271,.267,.256,.253,.245,.243,.243,.244, 3158.   
     L     .265,.274,.287,.293,.290,.281,.267,.262,.256,.251,.253,.258, 3159.   
     M     .293,.307,.324,.323,.315,.298,.284,.275,.268,.263,.271,.278, 3160.   
     N     .326,.348,.370,.363,.347,.320,.304,.290,.281,.278,.291,.306, 3161.   
     O     .357,.391,.412,.404,.380,.347,.322,.303,.296,.296,.313,.334, 3162.   
     P     .381,.431,.447,.439,.412,.363,.331,.311,.301,.308,.331,.353, 3163.   
     Q     .389,.449,.470,.456,.417,.363,.329,.306,.296,.308,.331,.354, 3164.   
     R     .382,.441,.462,.454,.413,.354,.325,.296,.289,.301,.319,.343/ 3165.   
      DATA O3AVEM/                                                      3166.   
     A     .309,.290,.288,.289,.293,.292,.299,.306,.313,.325,.374,.350, 3167.   
     B     .306,.289,.286,.285,.296,.291,.300,.308,.316,.326,.369,.339, 3168.   
     C     .297,.284,.281,.285,.288,.290,.302,.308,.315,.324,.355,.323, 3169.   
     D     .287,.278,.272,.275,.277,.284,.295,.300,.306,.316,.333,.304, 3170.   
     E     .273,.266,.261,.263,.267,.274,.284,.288,.292,.302,.311,.286, 3171.   
     F     .260,.253,.250,.252,.253,.261,.268,.273,.275,.284,.288,.269, 3172.   
     G     .247,.244,.241,.245,.243,.250,.254,.260,.260,.270,.268,.254, 3173.   
     H     .238,.234,.235,.242,.239,.243,.244,.250,.249,.255,.253,.245, 3174.   
     I     .231,.231,.238,.244,.242,.246,.238,.242,.239,.243,.242,.239, 3175.   
     J     .236,.238,.247,.257,.254,.253,.245,.244,.237,.235,.235,.236, 3176.   
     K     .250,.254,.263,.270,.266,.264,.254,.250,.244,.239,.237,.243, 3177.   
     L     .270,.279,.289,.290,.285,.279,.267,.261,.256,.250,.251,.258, 3178.   
     M     .301,.317,.329,.322,.314,.298,.285,.277,.270,.263,.270,.282, 3179.   
     N     .342,.367,.380,.369,.351,.326,.309,.294,.286,.284,.295,.314, 3180.   
     O     .380,.412,.424,.411,.388,.357,.331,.311,.303,.302,.325,.347, 3181.   
     P     .398,.448,.457,.449,.419,.373,.343,.318,.309,.314,.341,.366, 3182.   
     Q     .396,.456,.480,.466,.424,.370,.338,.311,.303,.311,.336,.363, 3183.   
     R     .384,.442,.464,.456,.414,.358,.327,.297,.290,.302,.322,.344/ 3184.   
      DATA O3AVEN/                                                      3185.   
     A     .311,.291,.287,.288,.293,.292,.297,.305,.312,.325,.373,.350, 3186.   
     B     .307,.290,.286,.285,.293,.292,.300,.305,.315,.326,.366,.341, 3187.   
     C     .300,.287,.283,.282,.288,.292,.300,.306,.313,.324,.351,.323, 3188.   
     D     .290,.281,.274,.276,.279,.285,.293,.298,.303,.315,.330,.308, 3189.   
     E     .276,.272,.265,.264,.267,.274,.281,.287,.288,.302,.309,.289, 3190.   
     F     .263,.259,.254,.253,.257,.262,.267,.272,.274,.285,.287,.273, 3191.   
     G     .252,.247,.244,.248,.247,.252,.254,.260,.262,.270,.268,.259, 3192.   
     H     .243,.238,.239,.244,.241,.245,.245,.251,.251,.257,.253,.249, 3193.   
     I     .236,.233,.238,.244,.244,.246,.238,.243,.242,.245,.243,.242, 3194.   
     J     .237,.241,.247,.256,.255,.254,.245,.245,.242,.234,.234,.236, 3195.   
     K     .252,.259,.266,.271,.269,.269,.257,.256,.251,.242,.240,.245, 3196.   
     L     .277,.286,.296,.298,.292,.290,.276,.275,.267,.259,.259,.267, 3197.   
     M     .323,.342,.352,.339,.333,.319,.303,.298,.288,.280,.285,.296, 3198.   
     N     .374,.403,.413,.392,.376,.351,.332,.319,.306,.303,.317,.340, 3199.   
     O     .408,.448,.448,.433,.410,.375,.351,.330,.317,.318,.343,.368, 3200.   
     P     .418,.467,.473,.464,.426,.383,.347,.328,.316,.319,.347,.376, 3201.   
     Q     .402,.459,.482,.474,.426,.374,.343,.313,.306,.313,.338,.368, 3202.   
     R     .384,.440,.463,.458,.415,.360,.328,.299,.291,.301,.319,.344/ 3203.   
      DATA O3AVEO/                                                      3204.   
     A     .313,.291,.288,.288,.292,.292,.298,.305,.312,.324,.364,.351, 3205.   
     B     .311,.294,.289,.286,.294,.293,.302,.306,.316,.326,.358,.345, 3206.   
     C     .308,.296,.291,.286,.294,.297,.303,.310,.316,.330,.354,.331, 3207.   
     D     .301,.292,.284,.282,.286,.295,.301,.307,.310,.326,.334,.318, 3208.   
     E     .290,.283,.274,.273,.276,.286,.291,.297,.299,.314,.314,.302, 3209.   
     F     .280,.272,.266,.263,.264,.272,.277,.283,.286,.297,.295,.286, 3210.   
     G     .267,.261,.256,.254,.255,.260,.263,.268,.272,.280,.276,.271, 3211.   
     H     .254,.250,.249,.249,.247,.251,.251,.256,.259,.264,.261,.258, 3212.   
     I     .242,.242,.243,.245,.244,.248,.242,.247,.248,.252,.248,.248, 3213.   
     J     .237,.242,.249,.256,.255,.255,.245,.244,.243,.237,.236,.236, 3214.   
     K     .253,.256,.267,.271,.270,.270,.259,.258,.252,.245,.242,.248, 3215.   
     L     .279,.283,.296,.296,.294,.292,.280,.279,.269,.260,.260,.268, 3216.   
     M     .327,.339,.357,.345,.338,.328,.319,.309,.293,.284,.285,.302, 3217.   
     N     .386,.409,.421,.405,.388,.363,.346,.332,.314,.311,.319,.348, 3218.   
     O     .419,.450,.459,.445,.418,.384,.361,.338,.322,.320,.340,.373, 3219.   
     P     .419,.461,.473,.468,.423,.358,.358,.331,.316,.319,.343,.376, 3220.   
     Q     .401,.453,.477,.469,.423,.375,.345,.314,.307,.312,.333,.361, 3221.   
     R     .382,.437,.461,.455,.415,.361,.329,.299,.291,.301,.316,.341/ 3222.   
      DATA O3AVEP/                                                      3223.   
     A     .314,.293,.289,.290,.292,.294,.299,.305,.312,.323,.363,.352, 3224.   
     B     .315,.298,.293,.290,.294,.299,.303,.307,.316,.324,.365,.350, 3225.   
     C     .315,.303,.296,.291,.300,.306,.311,.316,.323,.336,.360,.341, 3226.   
     D     .308,.301,.293,.291,.297,.308,.312,.318,.324,.337,.345,.329, 3227.   
     E     .299,.292,.284,.283,.285,.299,.306,.311,.317,.326,.327,.314, 3228.   
     F     .285,.280,.272,.272,.274,.284,.293,.296,.301,.308,.306,.297, 3229.   
     G     .272,.266,.262,.261,.262,.269,.275,.280,.283,.289,.284,.280, 3230.   
     H     .256,.253,.251,.251,.251,.255,.256,.264,.266,.271,.267,.263, 3231.   
     I     .241,.242,.244,.245,.245,.248,.245,.251,.251,.255,.252,.251, 3232.   
     J     .236,.239,.247,.253,.253,.251,.242,.244,.239,.237,.235,.236, 3233.   
     K     .248,.250,.262,.267,.264,.262,.254,.250,.244,.240,.235,.239, 3234.   
     L     .268,.270,.286,.287,.284,.278,.267,.264,.256,.250,.245,.256, 3235.   
     M     .301,.308,.329,.322,.317,.300,.297,.281,.272,.264,.263,.279, 3236.   
     N     .351,.362,.380,.372,.360,.337,.320,.305,.295,.285,.287,.316, 3237.   
     O     .383,.406,.427,.415,.391,.365,.345,.324,.310,.304,.310,.342, 3238.   
     P     .393,.428,.450,.441,.404,.373,.353,.324,.310,.310,.321,.356, 3239.   
     Q     .387,.435,.461,.456,.412,.370,.341,.313,.303,.306,.321,.353, 3240.   
     R     .381,.432,.457,.452,.413,.361,.328,.299,.291,.298,.314,.338/ 3241.   
      DATA O3AVEQ/                                                      3242.   
     A     .315,.293,.289,.291,.293,.295,.298,.305,.312,.323,.362,.354, 3243.   
     B     .316,.301,.295,.291,.294,.300,.303,.307,.316,.322,.361,.350, 3244.   
     C     .318,.305,.297,.292,.298,.306,.311,.314,.324,.334,.354,.340, 3245.   
     D     .309,.301,.292,.289,.295,.305,.312,.317,.326,.335,.343,.326, 3246.   
     E     .295,.288,.279,.279,.284,.297,.305,.305,.316,.321,.324,.310, 3247.   
     F     .279,.272,.266,.269,.272,.281,.289,.291,.299,.303,.305,.293, 3248.   
     G     .263,.259,.254,.257,.259,.266,.273,.276,.281,.285,.284,.277, 3249.   
     H     .247,.246,.244,.248,.247,.252,.253,.261,.265,.269,.267,.259, 3250.   
     I     .235,.236,.239,.244,.243,.246,.243,.247,.251,.253,.249,.246, 3251.   
     J     .231,.234,.243,.250,.251,.247,.240,.238,.233,.234,.232,.233, 3252.   
     K     .242,.244,.257,.262,.260,.255,.247,.243,.235,.235,.228,.233, 3253.   
     L     .257,.263,.278,.280,.275,.269,.258,.252,.242,.239,.235,.243, 3254.   
     M     .280,.288,.308,.307,.299,.287,.274,.267,.255,.250,.246,.259, 3255.   
     N     .309,.319,.348,.340,.332,.309,.293,.286,.273,.264,.261,.282, 3256.   
     O     .339,.357,.388,.376,.360,.334,.320,.305,.289,.282,.279,.306, 3257.   
     P     .365,.393,.424,.411,.386,.355,.340,.316,.300,.303,.297,.329, 3258.   
     Q     .375,.415,.445,.439,.404,.365,.336,.310,.298,.299,.306,.338, 3259.   
     R     .379,.428,.453,.447,.412,.360,.326,.298,.291,.296,.310,.335/ 3260.   
      DATA O3AVER/                                                      3261.   
     A     .316,.295,.291,.292,.292,.296,.299,.305,.313,.323,.361,.355, 3262.   
     B     .317,.301,.296,.292,.292,.300,.302,.305,.314,.319,.358,.348, 3263.   
     C     .316,.303,.295,.289,.291,.301,.306,.307,.317,.324,.348,.336, 3264.   
     D     .303,.294,.286,.283,.285,.296,.304,.304,.313,.322,.333,.318, 3265.   
     E     .283,.277,.272,.272,.273,.284,.290,.296,.302,.309,.314,.299, 3266.   
     F     .265,.262,.259,.259,.259,.268,.274,.282,.286,.293,.293,.279, 3267.   
     G     .252,.249,.248,.249,.247,.253,.258,.265,.272,.277,.273,.265, 3268.   
     H     .241,.238,.240,.242,.241,.244,.246,.252,.257,.260,.256,.249, 3269.   
     I     .231,.229,.238,.241,.241,.242,.237,.242,.244,.247,.242,.239, 3270.   
     J     .231,.233,.242,.249,.251,.246,.237,.235,.230,.230,.229,.230, 3271.   
     K     .241,.250,.257,.265,.262,.257,.245,.243,.234,.230,.229,.231, 3272.   
     L     .260,.273,.281,.285,.280,.272,.257,.256,.245,.238,.237,.245, 3273.   
     M     .285,.302,.312,.314,.305,.294,.278,.277,.262,.252,.251,.262, 3274.   
     N     .310,.331,.347,.346,.336,.320,.303,.298,.281,.267,.267,.283, 3275.   
     O     .331,.354,.383,.378,.364,.342,.324,.315,.293,.278,.279,.297, 3276.   
     P     .350,.379,.414,.398,.381,.343,.335,.317,.299,.287,.285,.311, 3277.   
     Q     .367,.404,.436,.428,.399,.361,.332,.307,.295,.293,.298,.327, 3278.   
     R     .376,.424,.450,.442,.409,.358,.326,.296,.290,.294,.306,.332/ 3279.   
C                                                                       3280.   
      DIMENSION AO3AVE(18,12),SO3JF(11,19),SO3SO(11,19)                 3281.   
      DATA AO3AVE/            .3148,.3160,.3171,.3159,.3027,.2824,.2645,3282.   
     A.2493,.2376,.2344,.2455,.2667,.3038,.3467,.3753,.3842,.3817,.3780,3283.   
     B.2926,.2959,.3008,.3035,.2943,.2763,.2600,.2463,.2366,.2366,.2500,3284.   
     C.2735,.3166,.3661,.4076,.4270,.4310,.4309,.2904,.2937,.2974,.2959,3285.   
     D.2869,.2704,.2561,.2454,.2403,.2443,.2590,.2844,.3293,.3803,.4210,3286.   
     E.4439,.4534,.4539,.2918,.2943,.2965,.2940,.2834,.2687,.2561,.2476,3287.   
     F.2450,.2538,.2676,.2888,.3259,.3692,.4077,.4325,.4454,.4476,.2951,3288.   
     G.2979,.2994,.3001,.2904,.2731,.2575,.2467,.2441,.2548,.2675,.2873,3289.   
     H.3181,.3517,.3828,.4002,.4080,.4096,.2960,.3012,.3084,.3132,.3044,3290.   
     I.2852,.2660,.2515,.2465,.2521,.2641,.2802,.3023,.3257,.3417,.3457,3291.   
     J.3521,.3517,.3008,.3070,.3153,.3211,.3127,.2934,.2714,.2545,.2437,3292.   
     K.2440,.2528,.2665,.2875,.3064,.3191,.3222,.3210,.3201,.3074,.3126,3293.   
     L.3221,.3276,.3211,.3015,.2783,.2603,.2478,.2431,.2499,.2624,.2784,3294.   
     M.2928,.3024,.3017,.2954,.2914,.3156,.3224,.3326,.3391,.3300,.3071,3295.   
     N.2827,.2632,.2489,.2399,.2455,.2566,.2720,.2854,.2939,.2931,.2889,3296.   
     O.2854,.3282,.3354,.3456,.3504,.3368,.3124,.2899,.2692,.2532,.2389,3297.   
     P.2415,.2521,.2672,.2844,.2967,.3003,.2986,.2966,.3723,.3713,.3661,3298.   
     Q.3538,.3332,.3072,.2826,.2626,.2481,.2359,.2373,.2489,.2700,.2936,3299.   
     R.3113,.3172,.3154,.3130,.3554,.3533,.3467,.3353,.3146,.2925,.2723,3300.   
     S.2562,.2450,.2350,.2387,.2554,.2828,.3140,.3331,.3406,.3408,.3351/3301.   
C                                                                       3302.   
      DATA SO3JF/                                                       3303.   
     A     13.0,12.3,11.7,10.5,8.90,6.20,4.50,3.30,2.20,1.80,1.00,      3304.   
     B     13.6,12.9,11.9,10.3,8.30,6.10,4.45,3.40,2.50,1.85,1.00,      3305.   
     C     14.8,13.9,12.8,10.3,8.00,6.00,4.55,3.60,2.70,1.90,1.00,      3306.   
     D     16.6,15.1,14.0,11.0,7.95,6.00,4.65,3.70,2.95,1.95,1.00,      3307.   
     E     18.1,16.0,14.6,12.0,8.00,6.00,4.80,3.75,3.00,1.98,1.00,      3308.   
     F     18.3,16.3,14.8,12.6,8.20,6.15,4.80,3.80,3.05,2.00,1.00,      3309.   
     G     17.3,16.1,14.7,12.7,9.10,6.10,4.70,3.75,3.00,2.00,1.00,      3310.   
     H     16.3,15.5,14.5,12.6,9.00,6.00,4.55,3.65,2.95,1.98,1.00,      3311.   
     I     15.7,14.9,14.1,12.4,8.70,5.90,4.40,3.45,2.80,1.96,1.00,      3312.   
     J     15.3,14.1,13.5,12.2,8.30,5.85,4.25,3.40,2.75,1.95,1.00,      3313.   
     K     15.6,14.9,14.0,12.4,9.00,6.10,4.55,3.50,2.85,1.96,1.00,      3314.   
     L     17.4,16.6,16.0,14.0,10.0,7.30,5.10,3.90,3.00,1.97,1.00,      3315.   
     M     17.6,18.3,17.8,15.8,12.3,9.00,6.05,4.40,3.20,1.97,1.00,      3316.   
     N     16.0,16.9,17.8,16.8,15.2,12.0,7.90,5.10,3.65,1.97,1.00,      3317.   
     O     12.3,13.8,15.7,16.2,16.2,14.8,10.0,6.00,4.00,1.96,1.00,      3318.   
     P     12.0,11.9,12.0,13.8,14.3,14.3,12.0,6.80,4.30,1.95,1.00,      3319.   
     Q     11.9,11.8,11.7,11.6,11.8,12.0,10.3,7.20,4.50,1.90,1.00,      3320.   
     R     11.6,11.5,11.4,11.2,11.0,10.4,9.00,7.20,4.15,1.85,1.00,      3321.   
     S     11.2,10.9,10.7,10.5,10.0,9.75,8.60,7.00,3.80,1.80,1.00/      3322.   
      DATA SO3SO/                                                       3323.   
     A     10.5,10.5,10.5,10.6,10.5,10.3,8.20,4.80,3.10,1.90,1.00,      3324.   
     B     11.5,11.5,11.6,12.1,12.1,10.8,8.05,4.95,3.40,1.92,1.00,      3325.   
     C     12.7,13.8,14.0,14.1,12.9,10.9,7.95,5.10,3.70,1.96,1.00,      3326.   
     D     15.4,15.9,16.0,15.4,13.2,10.7,7.40,5.15,3.85,1.98,1.00,      3327.   
     E     17.9,18.0,17.4,16.1,13.0,10.0,6.70,4.90,3.80,1.99,1.00,      3328.   
     F     18.3,18.6,17.8,16.1,12.1,9.10,5.95,4.80,3.70,2.00,1.00,      3329.   
     G     18.6,18.5,17.8,15.9,11.1,8.00,5.55,4.40,3.45,2.00,1.00,      3330.   
     H     18.2,18.1,17.2,15.1,10.3,7.40,5.10,4.00,3.10,1.99,1.00,      3331.   
     I     17.5,16.8,16.2,14.0,9.90,7.00,4.90,3.85,2.95,1.98,1.00,      3332.   
     J     16.5,15.8,15.0,12.9,9.40,6.65,4.80,3.70,2.90,1.96,1.00,      3333.   
     K     16.3,15.8,15.0,12.9,9.20,6.80,5.00,3.85,2.95,1.96,1.00,      3334.   
     L     16.4,16.2,15.8,14.0,9.80,7.10,5.10,3.95,3.00,1.96,1.00,      3335.   
     M     16.6,16.5,16.2,14.8,10.8,7.75,5.50,4.05,3.05,1.97,1.00,      3336.   
     N     16.5,16.6,16.5,16.0,12.1,9.00,6.00,4.40,3.10,1.97,1.00,      3337.   
     O     15.8,16.2,16.4,16.1,14.2,10.9,6.60,4.50,3.20,1.97,1.00,      3338.   
     P     12.2,14.2,15.5,15.3,14.7,12.4,7.40,4.70,3.10,1.96,1.00,      3339.   
     Q     11.6,11.9,12.1,14.0,13.9,12.3,8.00,4.40,2.95,1.90,1.00,      3340.   
     R     11.2,11.2,11.4,11.6,11.8,10.9,8.00,3.95,2.60,1.87,1.00,      3341.   
     S     11.0,10.8,10.5,10.3,10.1,9.70,7.00,3.65,2.20,1.80,1.00/      3342.   
C                                                                       3343.   
      DIMENSION XJDMO(14),HKMSPR(14),HKMAUT(14)                         3344.   
      DIMENSION CNCAUT(14),CNCSPR(14),DEGLAT(14)                        3345.   
      DATA DEGLAT/-85.0,-71.0,-59.0,-47.0,-35.0,-22.0,-9.0,             3346.   
     +            9.0,22.0,35.0,47.0,59.0,71.0,85.0/                    3347.   
      DATA XJDMO/-15.0,16.0,45.0,75.0,105.0,136.0,166.0,197.0,228.0     3348.   
     +           ,258.0,289.0,319.0,350.0,381.0/                        3349.   
      DATA HKMSPR/18.5,18.5,19.0,23.5,24.0,24.5,26.5,                   3350.   
     +            26.5,25.0,22.5,21.0,20.0,18.5,16.5/                   3351.   
      DATA HKMAUT/16.5,18.5,20.0,21.0,22.5,25.0,26.5,                   3352.   
     +            26.5,24.5,24.0,23.5,19.0,18.5,18.5/                   3353.   
      DATA CNCSPR/0.0181,0.0212,0.0187,0.0167,0.0162,0.0183,0.0175,     3354.   
     +            0.0187,0.0200,0.0196,0.0225,0.0291,0.0287,0.0300/     3355.   
      DATA CNCAUT/0.0300,0.0287,0.0291,0.0225,0.0196,0.0200,0.0187,     3356.   
     +            0.0175,0.0183,0.0162,0.0167,0.0187,0.0212,0.0181/     3357.   
C                                                                       3358.   
      DIMENSION PLBSO3(11),SOJDAY(6),PMLAT(6)                           3359.   
      DATA PLBSO3/10.0,7.0,5.0,3.0,2.0,1.5,1.0,0.7,0.5,0.3,0.1/         3360.   
      DATA SOJDAY/-91.,31.,92.,213.,274.,396./                          3361.   
      DATA PMLAT/1.,1.,-1.,-1.,1.,1./                                   3362.   
      DIMENSION AO3JIM(144),O3LB(40),PLB0(40)                           3363.   
      DIMENSION CONCS(144),CONCA(144),BHKMS(144),BHKMA(144)             3364.   
      DIMENSION WTJLAT(144),WTJLON(144),ILATIJ(144),ILONIJ(144)         3365.   
      DIMENSION WTLSEP(144),WTLJAN(144),LSEPJ(144),LJANJ(144)           3366.   
      DATA ACMMGG/2.37251E-4/,ACMPKM/7.1509E-4/,H10MB/31.05467/         3367.   
      DATA A,B,C,D/0.331,23.0,4.553,5.23/                               3368.   
      LOGICAL SKIPI                                                     3369.   
C                                                                       3370.   
C-----------------------------------------------------------------------3371.   
C----SET O3 VERTICAL PROFILE PARAMETERS FOR LATITUDE GCM GRID POINTS    3372.   
C-----------------------------------------------------------------------3373.   
      SKIPI =.FALSE.                                                    3374.   
      IF(ABS(FLONO3).LT.1.E-04) SKIPI =.TRUE.                           3375.   
      DO 100 L=1,NL                                                     3376.   
  100 PLB0(L)=PLB(L)                                                    3377.   
      DO 103 J=1,JMLAT                                                  3378.   
      DLATJ=DLAT(J)                                                     3379.   
      ILATI=(DLATJ+95.001)/10.                                          3380.   
      IF(ILATI.LT. 1) ILATI= 1                                          3381.   
      IF(ILATI.GT.17) ILATI=17                                          3382.   
      ILATIJ(J)=ILATI                                                   3383.   
      LATD=ILATI*10-95                                                  3384.   
      WTJL=(DLATJ-LATD)*0.1                                             3385.   
      WTJLAT(J)=WTJL                                                    3386.   
      DO 101 JJ=2,14                                                    3387.   
      II=JJ-1                                                           3388.   
      IF(DLATJ.LE.DEGLAT(JJ)) GO TO 102                                 3389.   
  101 CONTINUE                                                          3389.1  
      JJ=14                                                             3390.   
  102 WTJJ=(DLATJ-DEGLAT(II))/(DEGLAT(JJ)-DEGLAT(II))                   3391.   
      WTII=1.-WTJJ                                                      3392.   
      CONCS(J)=WTII*CNCSPR(II)+WTJJ*CNCSPR(JJ)                          3393.   
      CONCA(J)=WTII*CNCAUT(II)+WTJJ*CNCAUT(JJ)                          3394.   
      BHKMS(J)=WTII*HKMSPR(II)+WTJJ*HKMSPR(JJ)                          3395.   
  103 BHKMA(J)=WTII*HKMAUT(II)+WTJJ*HKMAUT(JJ)                          3396.   
C                                                                       3397.   
      DO 104 I=1,IMLON                                                  3398.   
      DLONI=DLON(I)                                                     3399.   
      ILONG=DLONI/20.0                                                  3400.   
      WTJLG=(DLONI-ILONG*20)/20.0                                       3401.   
      WTJLON(I)=WTJLG                                                   3402.   
      WTILG=1.-WTJLG                                                    3403.   
      ILONG=ILONG+1                                                     3404.   
      JLONG=ILONG+1                                                     3405.   
      IF(ILONG.GT.18) ILONG=18                                          3406.   
      IF(ILONG.GT.17) JLONG=1                                           3407.   
  104 ILONIJ(I)=ILONG                                                   3408.   
      NLAY=LASTVC/100000                                                3409.   
      NATM=(LASTVC-NLAY*100000)/10000                                   3410.   
      IF(NATM.GT.0) GO TO 106                                           3411.   
C                                                                       3412.   
      O3B=0.343                                                         3413.   
      DO 105 L=1,NL                                                     3414.   
      HLT=HLB(L+1)                                                      3415.   
      O3T=A*(1.0+EXP(-B/C))/(1.0+EXP((HLT-B)/C))+(0.343-A)*EXP(-HLT/D)  3416.   
      U0GAS(L,3)=(O3B-O3T)                                              3417.   
  105 O3B=O3T                                                           3418.   
C                                                                       3419.   
  106 AO3J=0.0                                                          3420.   
      RETURN                                                            3421.   
C-----------------------------------------------------------------------3422.   
      ENTRY O3DDAY                                                      3423.   
C-----------------------------------------------------------------------3424.   
      XJDAY=JDAY                                                        3425.   
      WTAUT=(XJDAY-91.)/213.                                            3426.   
      IF(XJDAY.LT. 91.) WTAUT=( 91.-XJDAY)/152.                         3427.   
      IF(XJDAY.GT.304.) WTAUT=(456.-XJDAY)/152.                         3428.   
      WTSPR=1.-WTAUT                                                    3429.   
      DO 200 JMO=1,14                                                   3430.   
      XJDMJ=XJDMO(JMO)                                                  3431.   
      IF(XJDAY.LT.XJDMJ) GO TO 201                                      3432.   
  200 XJDMI=XJDMJ                                                       3433.   
      XJDMI=XJDMO(13)                                                   3434.   
  201 DAYMO=XJDMJ-XJDMI                                                 3435.   
      WTJM=(XJDAY-XJDMI)/DAYMO                                          3436.   
      WTIM=1.-WTJM                                                      3437.   
      JMO=JMO-1                                                         3438.   
      IMO=JMO-1                                                         3439.   
      IF(IMO.LT.1) IMO=12                                               3440.   
      IF(JMO.GT.12) JMO=1                                               3441.   
      JJDAY=1                                                           3442.   
      SJDAY=SOJDAY(JJDAY)                                               3443.   
  202 JJDAY=JJDAY+1                                                     3444.   
      SIDAY=SJDAY                                                       3445.   
      SJDAY=SOJDAY(JJDAY)                                               3446.   
      IF(XJDAY.GT.SJDAY) GO TO 202                                      3447.   
      WTJAN=(XJDAY-SIDAY)/(SJDAY-SIDAY)                                 3448.   
      IF(JJDAY.EQ.3.OR.JJDAY.EQ.5) WTJAN=1.-WTJAN                       3449.   
      WTSEP=1.0-WTJAN                                                   3450.   
      DO 203 J=1,JMLAT                                                  3451.   
      DLATJ=DLAT(J)                                                     3452.   
      DLSEP=10.0+0.099999*DLATJ*PMLAT(JJDAY)                            3453.   
      DLJAN=10.0+0.099999*DLATJ*PMLAT(JJDAY-1)                          3454.   
      LSEP=DLSEP                                                        3455.   
      LJAN=DLJAN                                                        3456.   
      LJANJ(J)=LJAN                                                     3457.   
      LSEPJ(J)=LSEP                                                     3458.   
      WTLSEP(J)=DLSEP-LSEP                                              3459.   
  203 WTLJAN(J)=DLJAN-LJAN                                              3460.   
      IF(AO3J.GT.1.E-10) GO TO 400                                      3461.   
C                                                                       3462.   
C-----------------------------------------------------------------------3463.   
      ENTRY O3DLAT                                                      3464.   
C-----------------------------------------------------------------------3465.   
      ILATI=ILATIJ(JLAT)                                                3466.   
      WTJL=WTJLAT(JLAT)                                                 3467.   
      WTIL=1.-WTJL                                                      3468.   
      JLATI=ILATI+1                                                     3469.   
      LSEP=LSEPJ(JLAT)                                                  3470.   
      LJAN=LJANJ(JLAT)                                                  3471.   
      WTLS=WTLSEP(JLAT)                                                 3472.   
      WTLJ=WTLJAN(JLAT)                                                 3473.   
      AO3J=WTIM*(WTIL*AO3AVE(ILATI,IMO)+WTJL*AO3AVE(JLATI,IMO))         3474.   
     +    +WTJM*(WTIL*AO3AVE(ILATI,JMO)+WTJL*AO3AVE(JLATI,JMO))         3475.   
      BHKMJ=WTSPR*BHKMS(JLAT)+WTAUT*BHKMA(JLAT)                         3476.   
      CONCJ=WTSPR*CONCS(JLAT)+WTAUT*CONCA(JLAT)                         3477.   
      AO3JJ=AO3J                                                        3478.   
      IF(SKIPI) GO TO 400                                               3479.   
      DO 300 I=1,IMLON                                                  3480.   
      ILONG=ILONIJ(I)                                                   3481.   
      JLONG=ILONG+1                                                     3482.   
      IF(JLONG.GT.18) JLONG=1                                           3483.   
      WTJLG=WTJLON(I)                                                   3484.   
      WTILG=1.0-WTJLG                                                   3485.   
      AO3J=WTIM*(WTIL*(WTILG*O3AVE(IMO,ILATI,ILONG)                     3486.   
     +                +WTJLG*O3AVE(IMO,ILATI,JLONG))                    3487.   
     +          +WTJL*(WTILG*O3AVE(IMO,JLATI,ILONG)                     3488.   
     +                +WTJLG*O3AVE(IMO,JLATI,JLONG)))                   3489.   
     +    +WTJM*(WTIL*(WTILG*O3AVE(JMO,ILATI,ILONG)                     3490.   
     +                +WTJLG*O3AVE(JMO,ILATI,JLONG))                    3491.   
     +          +WTJL*(WTILG*O3AVE(JMO,JLATI,ILONG)                     3492.   
     +                +WTJLG*O3AVE(JMO,JLATI,JLONG)))                   3493.   
  300 AO3JIM(I)=AO3J                                                    3494.   
      AO3J=AO3JJ                                                        3495.   
C                                                                       3496.   
C-----------------------------------------------------------------------3497.   
      ENTRY O3DLON                                                      3498.   
C-----------------------------------------------------------------------3499.   
C                                                                       3500.   
      IF(SKIPI) RETURN                                                  3501.   
      AO3J=AO3JJ+ABS((AO3JIM(ILON)-AO3JJ))*FLONO3                       3502.   
C                                                                       3503.   
  400 CKMJ=0.25*AO3J/CONCJ                                              3504.   
      GTOP=0.0                                                          3505.   
      POI=0.0                                                           3506.   
      FI=0.0                                                            3507.   
      L=NL                                                              3508.   
      PLL=PLB0(L)                                                       3509.   
      J=12                                                              3510.   
  401 J=J-1                                                             3511.   
      IF(J.LT.1) GO TO 404                                              3512.   
      POJ=PLBSO3(J)                                                     3513.   
      FJ=WTSEP*(WTLS*SO3SO(J,LSEP+1)+(1.-WTLS)*SO3SO(J,LSEP))           3514.   
     +  +WTJAN*(WTLJ*SO3JF(J,LJAN+1)+(1.-WTLJ)*SO3JF(J,LJAN))           3515.   
  402 DP=POJ-POI                                                        3516.   
      IF(POJ.GT.PLL) GO TO 403                                          3517.   
      GTOP=GTOP+(FI+FJ)*DP*ACMMGG                                       3518.   
      POI=POJ                                                           3519.   
      FI=FJ                                                             3520.   
      GO TO 401                                                         3521.   
  403 FF=(FJ-FI)/DP                                                     3522.   
      DP=PLL-POI                                                        3523.   
      FF=FI+FF*DP                                                       3524.   
      GTOP=GTOP+(FI+FF)*DP*ACMMGG                                       3525.   
      POI=PLL                                                           3526.   
      FI=FF                                                             3527.   
      O3LB(L)=GTOP                                                      3528.   
      L=L-1                                                             3529.   
      PLL=PLB0(L)                                                       3530.   
      GO TO 402                                                         3531.   
  404 FI=FJ*ACMPKM                                                      3532.   
      HI=H10MB                                                          3533.   
      HJ=BHKMJ+CKMJ                                                     3534.   
      XPBC=EXP(-BHKMJ/CKMJ)                                             3535.   
      XPHC=EXP(HJ/CKMJ)                                                 3536.   
      DTERM=1.0+XPHC*XPBC                                               3537.   
      ATERM=(1.0+XPBC)/DTERM                                            3538.   
      FTERM=ATERM/DTERM*XPHC*XPBC/CKMJ                                  3539.   
      TTERM=AO3J-GTOP-FI*(HI-HJ)*0.5                                    3540.   
      AA=TTERM/(FTERM*(HI-HJ)*0.5+1.0-ATERM)                            3541.   
      FJ=AA*FTERM                                                       3542.   
      GTOPBC=GTOP+(FI+FJ)*(HI-HJ)*0.5-AA*ATERM                          3543.   
      TOP=AA*(1.0+XPBC)                                                 3544.   
      GO TO 406                                                         3545.   
  405 DH=HI-HJ                                                          3546.   
      FF=(FJ-FI)/DH                                                     3547.   
      DH=HI-H                                                           3548.   
      FF=FI+FF*DH                                                       3549.   
      GTOP=GTOP+(FI+FF)*DH*0.5                                          3550.   
      HI=H                                                              3551.   
      FI=FF                                                             3552.   
      O3LB(L)=GTOP                                                      3553.   
      L=L-1                                                             3554.   
  406 CONTINUE                                                          3555.   
      H=HLB(L)                                                          3556.   
      IF(H.GT.HJ) GO TO 405                                             3557.   
      O3LB(L)=TOP/(1.+XPBC*EXP(H/CKMJ))+GTOPBC                          3558.   
      L=L-1                                                             3559.   
      IF(L.GT.0) GO TO 406                                              3560.   
      O3LB(NLP)=0.                                                      3561.   
      DO 407 L=1,NL                                                     3562.   
  407 U0GAS(L,3)=(O3LB(L)-O3LB(L+1))                                    3563.   
      RETURN                                                            3564.   
      END                                                               3565.   
      BLOCK DATA                                                        3566.   

#include "B83XX.COM" 

C-----------------------------------------------------------------------3597.   
C              SEASONAL ALBEDOS FOR  11  VEGETATION TYPES               3598.   
C-----------------------------------------------------------------------3599.   
C                                                                       3600.   
      EQUIVALENCE                                                       3601.   
     +          (VADATA(1,1,1),ALVISK(1,1)),(VADATA(1,1,2),ALNIRK(1,1)) 3602.   
     +,         (VADATA(1,1,3),FIELDC(1,1)),(VADATA(1,4,3),VTMASK(1))   3603.   
C$$  +          (VADATA(1,1,1),ALMEAN(1,1)),(VADATA(1,1,2),RATIRV(1,1)) 3604.   
C$$  +,         (VADATA(1,1,3),VTMASK(1)),(VADATA(1,2,3),FIELDC(1,1))   3605.   
C                                                                       3606.   
      EQUIVALENCE                                                       3607.   
     +          (FEMTRA(1),ECLTRA),  (FZASRA(1),ZCLSRA)                 3608.   
     +,         (FEMTRA(2),EOCTRA),  (FZASRA(2),ZOCSRA)                 3609.   
     +,         (FEMTRA(3),ESNTRA),  (FZASRA(3),ZSNSRA)                 3610.   
     +,         (FEMTRA(4),EICTRA),  (FZASRA(4),ZICSRA)                 3611.   
     +,         (FEMTRA(5),EDSTRA),  (FZASRA(5),ZDSSRA)                 3612.   
     +,         (FEMTRA(6),EVGTRA),  (FZASRA(6),ZVGSRA)                 3613.   
C                                                                       3614.   
      EQUIVALENCE (IMG(1),IMGAS1),(IMG(2),IMGAS2)                       3615.   
      EQUIVALENCE (ILG(1),ILGAS1),(ILG(2),ILGAS2)                       3616.   
C                                                                       3617.   
      EQUIVALENCE (ID5(1),IDPROG),(ID5(2),ID2TRD),(ID5(3),ID3SRD)       3618.   
      EQUIVALENCE (ID5(4),ID4VEG),(ID5(5),ID5FOR)                       3619.   
C                                                                       3620.   
      EQUIVALENCE (FRC(1),FRACCC),(FRC(2),  FCHI),(FRC(3),  FCMI)       3621.   
     +                           ,(FRC(4),  FCLO),(FRC(5),  FCOV)       3622.   
C                                                                       3623.   
      DIMENSION ALVISK(11,4),ALNIRK(11,4)                               3624.   
C$$   DIMENSION ALMEAN(11,4),RATIRV(11,4)                               3625.   
      DIMENSION FIELDC(11,3),VTMASK(11)                                 3626.   
C                                                                       3627.   
C                          1       2       3       4                    3628.   
C                        WINTER  SPRING  SUMMER  AUTUMN                 3629.   
C                                                                       3630.   
      DATA ALVISK/                                                      3631.   
C        1     2     3     4     5     6     7     8     9    10    11  3632.   
C      DESRT TNDRA GRASS SHRUB TREES DECID EVERG RAINF ROCKS CROPS ALGAE3633.   
     1 .350, .067, .089, .089, .078, .100, .067, .061, .100, .070, .001,3634.   
     2 .350, .063, .100, .100, .073, .055, .067, .061, .100, .070, .001,3635.   
     3 .350, .085, .091, .139, .085, .058, .083, .061, .100, .070, .001,3636.   
     4 .350, .080, .090, .111, .064, .055, .061, .061, .100, .070, .001/3637.   
C                                                                       3638.   
      DATA ALNIRK/                                                      3639.   
C        1     2     3     4     5     6     7     8     9    10    11  3640.   
C      DESRT TNDRA GRASS SHRUB TREES DECID EVERG RAINF ROCKS CROPS ALGAE3641.   
     1 .350, .200, .267, .267, .233, .300, .200, .183, .100, .070, .001,3642.   
     2 .350, .206, .350, .300, .241, .218, .200, .183, .100, .070, .001,3643.   
     3 .350, .298, .364, .417, .298, .288, .250, .183, .100, .070, .001,3644.   
     4 .350, .255, .315, .333, .204, .218, .183, .183, .100, .070, .001/3645.   
C                                                                       3646.   
C$$   DATA ALMEAN/                                                      3647.   
C        1     2     3     4     5     6     7     8     9    10    11  3648.   
C      DESRT TNDRA GRASS SHRUB TREES DECID EVERG RAINF ROCKS CROPS ALGAE3649.   
C$$  1 .350, .120, .160, .160, .140, .180, .120, .110, .100, .070, .001,3650.   
C$$  2 .350, .120, .200, .180, .140, .120, .120, .110, .100, .070, .001,3651.   
C$$  3 .350, .170, .200, .250, .170, .150, .150, .110, .100, .070, .001,3652.   
C$$  4 .350, .150, .180, .200, .120, .120, .110, .110, .100, .070, .001/3653.   
C                                                                       3654.   
C$$   DATA RATIRV/                                                      3655.   
C        1     2     3     4     5     6     7     8     9    10    11  3656.   
C      DESRT TNDRA GRASS SHRUB TREES DECID EVERG RAINF ROCKS CROPS ALGAE3657.   
C$$  1 1.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 1.00, 3.50, 1.50,3658.   
C$$  2 1.00, 3.30, 3.50, 3.00, 3.30, 4.00, 3.00, 3.00, 1.00, 3.50, 1.50,3659.   
C$$  3 1.00, 3.50, 4.00, 3.00, 3.50, 5.00, 3.00, 3.00, 1.00, 3.50, 1.50,3660.   
C$$  4 1.00, 3.20, 3.50, 3.00, 3.20, 4.00, 3.00, 3.00, 1.00, 3.50, 1.50/3661.   
C                                                                       3662.   
      DATA FIELDC/                                                      3663.   
C          (KG/M**2)                                                    3664.   
C        1     2     3     4     5     6     7     8     9    10    11  3665.   
C      DESRT TNDRA GRASS SHRUB TREES DECID EVERG RAINF ROCKS CROPS ALGAE3666.   
     1 10.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 200., 10.0, 30.0, 999.,3667.   
     2 10.0, 200., 200., 300., 300., 450., 450., 450., 10.0, 200., 999.,3668.   
     3 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0/3669.   
C                                                                       3670.   
      DATA VTMASK/                                                      3671.   
C          (KG/M**2)                                                    3672.   
C        1     2     3     4     5     6     7     8     9    10    11  3673.   
C      DESRT TNDRA GRASS SHRUB TREES DECID EVERG RAINF ROCKS CROPS ALGAE3674.   
     4 10.0, 20.0, 20.0, 50.0, 200., 500.,1000.,2500., 10.0, 30.0, .001/3676.   
C                                                                       3677.   
C                                                                       3678.   
      DATA DLAT/                                                        3679.   
     +-90.000000,-82.173913,-74.347826,-66.521739,-58.695652,-50.869565,3680.   
     +-43.043478,-35.217391,-27.391304,-19.565217,-11.739130,- 3.913043,3681.   
     +  3.913043, 11.739130, 19.565217, 27.391304, 35.217391, 43.043478,3682.   
     + 50.869565, 58.695652, 66.521739, 74.347826, 82.173913, 90.000000,3683.   
     + 22*0.0000/                                                       3684.   
C                                                                       3685.   
      DATA DLON/                                                        3686.   
     +   0.0,  10.0,  20.0,  30.0,  40.0,  50.0,  60.0,  70.0,  80.0,   3687.   
     +  90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0,   3688.   
     + 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 260.0,   3689.   
     + 270.0, 280.0, 290.0, 300.0, 310.0, 320.0, 330.0, 340.0, 350.0,   3690.   
     +36*0.0/                                                           3691.   
C                                                                       3692.   
C-----------------------------------------------------------------------3693.   
C     TRACE GAS REFERENCE AMOUNTS & DISTRIBUTIONS ARE DEFINED IN  SETGAS3694.   
C-----------------------------------------------------------------------3695.   
C                                                                       3696.   
C                                                                       3697.   
C                  H2O  CO2   O3   O2  NO2  N2O  CH4  F11  F12          3698.   
C                   1    2    3    4    5    6    7    8    9           3699.   
      DATA FULGAS/ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0          3700.   
     +           , 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/         3701.   
C                                                                       3702.   
C               GLOBAL OCEAN LAND DESERT  HAZE   TR1  TR2  TR3  TR4     3703.   
C                   1    2    3     4      5      6    7    8    9      3704.   
C                                                                       3705.   
      DATA FGOLDH/ 1.0, .68, .32, 1.E-20,1.E-20, 0.0, 0.0, 0.0, 0.0     3706.   
     +           , 1.0, .68, .32, 1.E-20,1.E-20, 0.0, 0.0, 0.0, 0.0/    3707.   
C                                                                       3708.   
      DATA    LASTVC/-123456/,  KFORCE/-123456789/                      3709.   
C                                                                       3710.   
C                                                                       3711.   
      DATA   TAUMIN/1.0E-04/, TLGRAD/ 1.0/, EOCTRA/1.0/, ZOCSRA/1.0/    3712.   
      DATA   FRACSL/1.0E-02/, TKCICE/258./, ESNTRA/1.0/, ZSNSRA/1.0/    3713.   
      DATA   RATQSL/1.0    /, FLONO3/ 0.0/, EICTRA/1.0/, ZICSRA/1.0/    3714.   
      DATA   FOGTSL/0.0    /, ECLTRA/1.00/, EDSTRA/1.0/, ZDSSRA/1.0/    3715.   
      DATA   PTLISO/2.5E+00/, ZCLSRA/1.00/, EVGTRA/1.0/, ZVGSRA/1.0/    3716.   
C                                                                       3717.   
      DATA   FMARCL/0.50/,    FCLDTR/1.0/,  NTRACE/0/,   IDPROG/0/      3718.   
      DATA   WETTRA/1.00/,    FCLDSR/1.0/,  ITR(1)/0/,   ID2TRD/0/      3719.   
      DATA   WETSRA/1.00/,    FALGAE/1.0/,  ITR(2)/0/,   ID3SRD/0/      3720.   
      DATA   DMOICE/10.0/,    FRAYLE/1.0/,  ITR(3)/0/,   ID4VEG/0/      3721.   
      DATA   DMLICE/10.0/,    LICETK/  0/,  ITR(4)/0/,   ID5FOR/0/      3722.   
C                                                                       3723.   
      DATA                                                   NV/ 8/     3724.   
      DATA   IMGAS1/1/,       KEEPRH/0/,    KGASSR/0/,   LAYRAD/ 3/     3725.   
      DATA   IMGAS2/3/,       KEEPAL/0/,    KAERSR/0/,       NL/12/     3726.   
      DATA   ILGAS1/2/,       ISOSCT/0/,    KFRACC/0/,      NLP/13/     3727.   
      DATA   ILGAS2/9/,       IHGSCT/0/,    MARCLD/0/,    JMLAT/24/     3728.   
      DATA   KWVCON/1/,       LAPGAS/1/,    NORMS0/1/,    IMLON/36/     3729.   
C                                                                       3730.   
      DATA     JYEAR/1958/,   JLAT/18/,     S0/1367.0/                  3731.   
      DATA      JDAY/   0/,   ILON/18/,   COSZ/0.5000/                  3732.   
C                                                                       3733.   
      DATA   POCEAN/0.700/,   TGO/288.15/,   AGESN/1.00/,    WMAG/2.00/ 3734.   
      DATA   PEARTH/0.100/,   TGE/288.15/,   SNOWE/0.30/,  WEARTH/0.00/ 3735.   
      DATA    POICE/0.100/,  TGOI/288.15/,  SNOWOI/0.10/,   ZOICE/10.0/ 3736.   
      DATA    PLICE/0.100/,  TGLI/288.15/,  SNOWLI/0.20/,  FRACCC/0.00/ 3737.   
      DATA                    TSL/288.15/                               3738.   
C                                                                       3739.   
      DATA PLB/                                                         3740.   
     +     1013.2500, 961.7485, 879.3460, 741.3219, 566.2166, 401.4117, 3741.   
     +      262.3575, 154.2043,  71.8018,  10.0000,   5.0000,   2.0000, 3742.   
     +        1.E-05,  27*0.00/                                         3743.   
C                                                                       3744.   
      DATA HLB/                                                         3745.   
     +         1.E-10,1.0,2.0,3.0,6.0,11.0,20.0,32.0,47.0,51.0,71.0     3746.   
     +        ,84.852,99.99,27*99.999/                                  3747.   
C                                                                       3748.   
      DATA TLB/40*250./                                                 3749.   
      DATA TLT/40*250./                                                 3750.   
      DATA TLM/40*250./                                                 3751.   
C                                                                       3752.   
      DATA U0GAS/360*0./                                                3753.   
      DATA ULGAS/360*0./                                                3754.   
C                                                                       3755.   
      DATA TRACER/160*0./                                               3756.   
      DATA CLDTAU/ 40*0./                                               3757.   
C                                                                       3758.   
      DATA SHL/40*0./                                                   3759.   
      DATA RHL/40*0./                                                   3760.   
C                                                                       3761.   
      DATA PVT/8*0.125,3*0.0/                                           3762.   
C                                                                       3763.   
      DATA SRBXAL/30*0./                                                3764.   
      DATA BXA/153*0./                                                  3765.   
C                                                                       3766.   
      DATA LUXGAS/1/                                                    3767.   
      DATA KALVIS/0/                                                    3768.   
      DATA MEANAL/0/                                                    3769.   
C                                                                       3770.   
C-----------------------------------------------------------------------3771.   
C  AEROSOL RADIATIVE PROPERTIES,COMPOSITION,TYPE & VERTICAL DISTRIBUTION3772.   
C-----------------------------------------------------------------------3773.   
C                                                                       3774.   
C                                BLOCKD INITIALIZED DEFAULT DATA        3775.   
C                                                                       3776.   
C                                                                       3785.   
      DIMENSION     QACID1(25),QACID2(25),QSLFT1(25),QSLFT2(25)         3786.   
     T             ,QBSLT1(25),QBSLT2(25),QSSALT(25),QDUST1(25)         3787.   
     T             ,QDUST2(25),QCARB1(25),QCARB2(25)                    3788.   
     T             ,SACID1(25),SACID2(25),SSLFT1(25),SSLFT2(25)         3789.   
     T             ,SBSLT1(25),SBSLT2(25),SSSALT(25),SDUST1(25)         3790.   
     T             ,SDUST2(25),SCARB1(25),SCARB2(25)                    3791.   
     T             ,CACID1(25),CACID2(25),CSLFT1(25),CSLFT2(25)         3792.   
     T             ,CBSLT1(25),CBSLT2(25),CSSALT(25),CDUST1(25)         3793.   
     T             ,CDUST2(25),CCARB1(25),CCARB2(25)                    3794.   
     T             ,QWATER(25),QICE25(25),SWATER(25),SICE25(25)         3795.   
     T             ,CWATER(25),CICE25(25)                               3796.   
C                                                                       3797.   
     S      ,XACID1(6),XACID2(6),XSLFT1(6),XSLFT2(6),XBSLT1(6),XBSLT2(6)3798.   
     S      ,XSSALT(6),XDUST1(6),XDUST2(6),XCARB1(6),XCARB2(6)          3799.   
     S      ,YACID1(6),YACID2(6),YSLFT1(6),YSLFT2(6),YBSLT1(6),YBSLT2(6)3800.   
     S      ,YSSALT(6),YDUST1(6),YDUST2(6),YCARB1(6),YCARB2(6)          3801.   
     S      ,ZACID1(6),ZACID2(6),ZSLFT1(6),ZSLFT2(6),ZBSLT1(6),ZBSLT2(6)3802.   
     S      ,ZSSALT(6),ZDUST1(6),ZDUST2(6),ZCARB1(6),ZCARB2(6)          3803.   
     S      ,XWATER(6),XICE25(6),YWATER(6),YICE25(6),ZWATER(6),ZICE25(6)3804.   
C                                                                       3805.   
      EQUIVALENCE (TRAQEX(1, 1),QACID1(1)),(TRAQEX(1, 2),QACID2(1))     3806.   
     1           ,(TRAQEX(1, 3),QSLFT1(1)),(TRAQEX(1, 4),QSLFT2(1))     3807.   
     2           ,(TRAQEX(1, 5),QBSLT1(1)),(TRAQEX(1, 6),QBSLT2(1))     3808.   
     3           ,(TRAQEX(1, 7),QSSALT(1)),(TRAQEX(1, 8),QDUST1(1))     3809.   
     4           ,(TRAQEX(1, 9),QDUST2(1)),(TRAQEX(1,10),QCARB1(1))     3810.   
     5           ,(TRAQEX(1,11),QCARB2(1))                              3811.   
C                                                                       3812.   
      EQUIVALENCE (TRAQSC(1, 1),SACID1(1)),(TRAQSC(1, 2),SACID2(1))     3813.   
     1           ,(TRAQSC(1, 3),SSLFT1(1)),(TRAQSC(1, 4),SSLFT2(1))     3814.   
     2           ,(TRAQSC(1, 5),SBSLT1(1)),(TRAQSC(1, 6),SBSLT2(1))     3815.   
     3           ,(TRAQSC(1, 7),SSSALT(1)),(TRAQSC(1, 8),SDUST1(1))     3816.   
     4           ,(TRAQSC(1, 9),SDUST2(1)),(TRAQSC(1,10),SCARB1(1))     3817.   
     5           ,(TRAQSC(1,11),SCARB2(1))                              3818.   
C                                                                       3819.   
      EQUIVALENCE (TRACOS(1, 1),CACID1(1)),(TRACOS(1, 2),CACID2(1))     3820.   
     1           ,(TRACOS(1, 3),CSLFT1(1)),(TRACOS(1, 4),CSLFT2(1))     3821.   
     2           ,(TRACOS(1, 5),CBSLT1(1)),(TRACOS(1, 6),CBSLT2(1))     3822.   
     3           ,(TRACOS(1, 7),CSSALT(1)),(TRACOS(1, 8),CDUST1(1))     3823.   
     4           ,(TRACOS(1, 9),CDUST2(1)),(TRACOS(1,10),CCARB1(1))     3824.   
     5           ,(TRACOS(1,11),CCARB2(1))                              3825.   
C                                                                       3826.   
      EQUIVALENCE (TRCQEX(1, 1),QWATER(1)),(TRCQEX(1, 2),QICE25(1))     3827.   
      EQUIVALENCE (TRCQSC(1, 1),SWATER(1)),(TRCQSC(1, 2),SICE25(1))     3828.   
      EQUIVALENCE (TRCCOS(1, 1),CWATER(1)),(TRCCOS(1, 2),CICE25(1))     3829.   
                                                                        3830.   
C                                                                       3831.   
      EQUIVALENCE (SRAQEX(1, 1),XACID1(1)),(SRAQEX(1, 2),XACID2(1))     3832.   
     1           ,(SRAQEX(1, 3),XSLFT1(1)),(SRAQEX(1, 4),XSLFT2(1))     3833.   
     2           ,(SRAQEX(1, 5),XBSLT1(1)),(SRAQEX(1, 6),XBSLT2(1))     3834.   
     3           ,(SRAQEX(1, 7),XSSALT(1)),(SRAQEX(1, 8),XDUST1(1))     3835.   
     4           ,(SRAQEX(1, 9),XDUST2(1)),(SRAQEX(1,10),XCARB1(1))     3836.   
     5           ,(SRAQEX(1,11),XCARB2(1))                              3837.   
C                                                                       3838.   
      EQUIVALENCE (SRAQSC(1, 1),YACID1(1)),(SRAQSC(1, 2),YACID2(1))     3839.   
     1           ,(SRAQSC(1, 3),YSLFT1(1)),(SRAQSC(1, 4),YSLFT2(1))     3840.   
     2           ,(SRAQSC(1, 5),YBSLT1(1)),(SRAQSC(1, 6),YBSLT2(1))     3841.   
     3           ,(SRAQSC(1, 7),YSSALT(1)),(SRAQSC(1, 8),YDUST1(1))     3842.   
     4           ,(SRAQSC(1, 9),YDUST2(1)),(SRAQSC(1,10),YCARB1(1))     3843.   
     5           ,(SRAQSC(1,11),YCARB2(1))                              3844.   
C                                                                       3845.   
      EQUIVALENCE (SRACOS(1, 1),ZACID1(1)),(SRACOS(1, 2),ZACID2(1))     3846.   
     1           ,(SRACOS(1, 3),ZSLFT1(1)),(SRACOS(1, 4),ZSLFT2(1))     3847.   
     2           ,(SRACOS(1, 5),ZBSLT1(1)),(SRACOS(1, 6),ZBSLT2(1))     3848.   
     3           ,(SRACOS(1, 7),ZSSALT(1)),(SRACOS(1, 8),ZDUST1(1))     3849.   
     4           ,(SRACOS(1, 9),ZDUST2(1)),(SRACOS(1,10),ZCARB1(1))     3850.   
     5           ,(SRACOS(1,11),ZCARB2(1))                              3851.   
C                                                                       3852.   
      EQUIVALENCE (SRCQEX(1, 1),XWATER(1)),(SRCQEX(1, 2),XICE25(1))     3853.   
      EQUIVALENCE (SRCQSC(1, 1),YWATER(1)),(SRCQSC(1, 2),YICE25(1))     3854.   
      EQUIVALENCE (SRCCOS(1, 1),ZWATER(1)),(SRCCOS(1, 2),ZICE25(1))     3855.   
                                                                        3856.   
C                                                                       3857.   
      DATA NGOLDH/5/,NAERO/11/                                          3858.   
C                                                                       3859.   
C-----------------------------------------------------------------------3860.   
C      COMPOSITION & VERTICAL DISTRIBUTION FOR 5 SPECIFIED AEROSOL TYPES3861.   
C-----------------------------------------------------------------------3862.   
C TYPE                                                                  3863.   
C    1   STRATOSPHERIC GLOBAL AEROSOL  A,B,C ARE GLOBAL AVERAGE VALUES  3864.   
C    2    TROPOSPHERIC  OCEAN AEROSOL  A,B,C ARE GLOBAL AVERAGE VALUES  3865.   
C    3    TROPOSPHERIC   LAND AEROSOL  A,B,C ARE GLOBAL AVERAGE VALUES  3866.   
C    4    TROPOSPHERIC DESERT AEROSOL  A,B,C ARE  LOCAL AVERAGE VALUES  3867.   
C    5    TROPOSPHERIC   HAZE AEROSOL  A,B,C ARE  LOCAL AVERAGE VALUES  3868.   
C                                                                       3869.   
C        1     2     3     4     5     6     7     8     9    10    11  3870.   
C      ACID1 OCT82 SLFT1 SLFT2 BSLT1 BSLT2 SSALT DUST1 DUST2 MAY82 CARB23871.   
      DATA AGOLDH/                                                      3872.   
     1 .012,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,3873.   
     2   .0,   .0, .018, .033, .012, .023, .011,   .0,   .0,   .0,   .0,3874.   
     3   .0,   .0, .031, .057, .021, .042,   .0,   .0,   .0,   .0, .018,3875.   
     4   .0,   .0,   .0,   .0,   .0,   .0,   .0, .300, .300,   .0,   .0,3876.   
     5   .0, .250,   .0,   .0,   .0,   .0,   .0, .300,   .0,   .0,   .0/3877.   
      DATA BGOLDH/                                                      3878.   
     1 20.0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,3879.   
     2   .0,   .0, 4.00, 0.00, 4.00, 1.00, 0.00,   .0,   .0,   .0,   .0,3880.   
     3   .0,   .0, 4.00, 0.00, 4.00, 0.00,   .0,   .0,   .0,   .0, 0.00,3881.   
     4   .0,   .0,   .0,   .0,   .0,   .0,   .0, 3.50, 0.00,   .0,   .0,3882.   
     5   .0, 0.00,   .0,   .0,   .0,   .0,   .0, 3.50,   .0,   .0,   .0/3883.   
      DATA CGOLDH/                                                      3884.   
     1 3.00,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,3885.   
     2   .0,   .0, 3.00, 1.00, 3.00,  0.5, 1.00,   .0,   .0,   .0,   .0,3886.   
     3   .0,   .0, 3.00, 1.00, 3.00, 1.00,   .0,   .0,   .0,   .0, 1.00,3887.   
     4   .0,   .0,   .0,   .0,   .0,   .0,   .0, 1.00, 1.00,   .0,   .0,3888.   
     5   .0, 1.00,   .0,   .0,   .0,   .0,   .0, 1.00,   .0,   .0,   .0/3889.   
C                                                                       3890.   
C-----------------------------------------------------------------------3891.   
C THERMAL RADIATION 25 K-INTERVAL MERGED AEROSOL DATA FOR QEXT,QSCA,COSB3892.   
C-----------------------------------------------------------------------3893.   
      DATA QACID1/                                                      3894.   
     +        0.04052,0.05895,0.08506,0.06673,0.05160,0.04437,0.03864,  3895.   
     +        0.02719,0.01668,0.01146,0.00705,0.03286,0.02449,0.03017,  3896.   
     +        0.03198,0.02891,0.02634,0.02366,0.02300,0.02271,0.02159,  3897.   
     +        0.08516,0.08825,0.08982,0.09284/                          3898.   
      DATA SACID1/                                                      3899.   
     +        0.00095,0.00361,0.00273,0.00226,0.00150,0.00141,0.00131,  3900.   
     +        0.00090,0.00049,0.00029,0.00014,0.00072,0.00049,0.00031,  3901.   
     +        0.00023,0.00023,0.00022,0.00020,0.00019,0.00018,0.00018,  3902.   
     +        0.00183,0.00201,0.00205,0.00207/                          3903.   
      DATA CACID1/                                                      3904.   
     +        0.11030,0.17256,0.17138,0.19696,0.19510,0.18945,0.18874,  3905.   
     +        0.18795,0.18313,0.17814,0.17075,0.10583,0.09756,0.08388,  3906.   
     +        0.07246,0.07266,0.07099,0.06873,0.06754,0.06661,0.06674,  3907.   
     +        0.11197,0.11068,0.10998,0.10852/                          3908.   
C                                                                       3909.   
      DATA QACID2/                                                      3910.   
     +        0.05764,0.15189,0.06264,0.04527,0.03973,0.03646,0.03375,  3911.   
     +        0.02163,0.01337,0.00979,0.00724,0.04076,0.03631,0.04273,  3912.   
     +        0.04072,0.03752,0.03290,0.03012,0.02968,0.02914,0.02763,  3913.   
     +        0.10731,0.12510,0.12901,0.13232/                          3914.   
      DATA SACID2/                                                      3915.   
     +        0.00367,0.00752,0.00264,0.00172,0.00188,0.00221,0.00225,  3916.   
     +        0.00134,0.00066,0.00034,0.00012,0.00237,0.00121,0.00084,  3917.   
     +        0.00080,0.00081,0.00074,0.00069,0.00067,0.00065,0.00064,  3918.   
     +        0.00674,0.00807,0.00825,0.00837/                          3919.   
      DATA CACID2/                                                      3920.   
     +        0.05720,0.11171,0.11850,0.11443,0.12325,0.13171,0.13500,  3921.   
     +        0.13575,0.13419,0.12666,0.10961,0.05186,0.04026,0.03219,  3922.   
     +        0.03060,0.03105,0.03041,0.02959,0.02911,0.02884,0.02901,  3923.   
     +        0.07145,0.07168,0.07134,0.07096/                          3924.   
C                                                                       3925.   
      DATA QSLFT1/                                                      3926.   
     +        0.15555,0.16333,0.16406,0.16396,0.16070,0.14074,0.11920,  3927.   
     +        0.09140,0.07341,0.06645,0.05871,0.15301,0.13456,0.15809,  3928.   
     +        0.16264,0.14805,0.12798,0.10588,0.09960,0.09604,0.08844,  3929.   
     +        0.35895,0.27430,0.26964,0.27183/                          3930.   
      DATA SSLFT1/                                                      3931.   
     +        0.13162,0.13152,0.11642,0.12932,0.10550,0.08323,0.07081,  3932.   
     +        0.05079,0.03287,0.02458,0.01871,0.12787,0.11183,0.09490,  3933.   
     +        0.08739,0.08716,0.08022,0.07182,0.06899,0.06700,0.06496,  3934.   
     +        0.13067,0.12933,0.12878,0.12808/                          3935.   
      DATA CSLFT1/                                                      3936.   
     +        0.52508,0.48102,0.59654,0.66259,0.66566,0.70224,0.71546,  3937.   
     +        0.69308,0.62819,0.55963,0.45811,0.52840,0.54500,0.51620,  3938.   
     +        0.50685,0.52475,0.54985,0.58351,0.59484,0.60203,0.61652,  3939.   
     +        0.45926,0.47060,0.47243,0.47178/                          3940.   
C                                                                       3941.   
      DATA QSLFT2/                                                      3942.   
     +        0.44109,0.37065,0.38095,0.40554,0.37738,0.32564,0.27970,  3943.   
     +        0.21687,0.17752,0.16154,0.14952,0.43239,0.38517,0.39512,  3944.   
     +        0.39098,0.36978,0.32960,0.28406,0.27042,0.26204,0.24771,  3945.   
     +        0.63665,0.59084,0.58844,0.59078/                          3946.   
      DATA SSLFT2/                                                      3947.   
     +        0.37818,0.31549,0.29505,0.33810,0.28074,0.22692,0.19562,  3948.   
     +        0.14289,0.09653,0.07449,0.06008,0.36685,0.33089,0.28296,  3949.   
     +        0.26185,0.26286,0.24369,0.22019,0.21220,0.20647,0.20093,  3950.   
     +        0.31870,0.30963,0.30762,0.30507/                          3951.   
      DATA CSLFT2/                                                      3952.   
     +        0.54586,0.50074,0.62826,0.69007,0.69596,0.73443,0.74600,  3953.   
     +        0.71846,0.64430,0.57291,0.47311,0.54977,0.56612,0.53939,  3954.   
     +        0.53105,0.54799,0.57221,0.60426,0.61497,0.62179,0.63518,  3955.   
     +        0.51454,0.52095,0.52268,0.52316/                          3956.   
C                                                                       3957.   
      DATA QBSLT1/                                                      3958.   
     +        0.19787,0.15206,0.14808,0.15505,0.14132,0.12508,0.10931,  3959.   
     +        0.07946,0.05659,0.04675,0.03801,0.20081,0.15823,0.15732,  3960.   
     +        0.15377,0.14273,0.13163,0.12005,0.11684,0.11523,0.11121,  3961.   
     +        0.36601,0.39099,0.39240,0.39274/                          3962.   
      DATA SBSLT1/                                                      3963.   
     +        0.09892,0.12369,0.09780,0.11017,0.08914,0.08577,0.07794,  3964.   
     +        0.05688,0.03912,0.03069,0.02440,0.09492,0.08277,0.05817,  3965.   
     +        0.04773,0.04970,0.04568,0.04058,0.03865,0.03717,0.03641,  3966.   
     +        0.07710,0.08232,0.08235,0.08163/                          3967.   
      DATA CBSLT1/                                                      3968.   
     +        0.54090,0.49369,0.59375,0.67539,0.69444,0.71623,0.71674,  3969.   
     +        0.69425,0.63125,0.57379,0.48766,0.54072,0.57272,0.57215,  3970.   
     +        0.57655,0.59243,0.60616,0.62323,0.62911,0.63253,0.63934,  3971.   
     +        0.51632,0.50380,0.50414,0.50666/                          3972.   
C                                                                       3973.   
      DATA QBSLT2/                                                      3974.   
     +        0.49004,0.35700,0.34009,0.38146,0.35476,0.32874,0.29258,  3975.   
     +        0.21726,0.16067,0.13571,0.11451,0.48169,0.40550,0.37263,  3976.   
     +        0.35312,0.33842,0.31466,0.28850,0.28051,0.27574,0.26813,  3977.   
     +        0.59495,0.63654,0.63850,0.63742/                          3978.   
      DATA SBSLT2/                                                      3979.   
     +        0.26833,0.30862,0.25309,0.29334,0.24644,0.24238,0.22164,  3980.   
     +        0.16459,0.11742,0.09480,0.07809,0.26006,0.23936,0.17265,  3981.   
     +        0.14418,0.15103,0.13960,0.12488,0.11925,0.11488,0.11275,  3982.   
     +        0.19766,0.20963,0.20969,0.20807/                          3983.   
      DATA CBSLT2/                                                      3984.   
     +        0.57850,0.51330,0.62334,0.70306,0.72063,0.74166,0.74111,  3985.   
     +        0.71466,0.64442,0.58410,0.49911,0.58174,0.60690,0.60535,  3986.   
     +        0.60954,0.62353,0.63716,0.65423,0.66019,0.66381,0.67030,  3987.   
     +        0.58670,0.57707,0.57759,0.58014/                          3988.   
C                                                                       3989.   
      DATA QSSALT/                                                      3990.   
     +        0.27651,0.36950,0.40122,0.39669,0.34286,0.33458,0.29978,  3991.   
     +        0.26075,0.26470,0.26660,0.28507,0.27114,0.23752,0.18761,  3992.   
     +        0.16890,0.17532,0.17705,0.17827,0.17801,0.17743,0.17914,  3993.   
     +        0.34241,0.33620,0.33607,0.33681/                          3994.   
      DATA SSSALT/                                                      3995.   
     +        0.27651,0.36950,0.40121,0.39659,0.34226,0.33245,0.29555,  3996.   
     +        0.22360,0.16290,0.13425,0.11177,0.27114,0.23751,0.18755,  3997.   
     +        0.16883,0.17526,0.17700,0.17823,0.17797,0.17739,0.17911,  3998.   
     +        0.34241,0.33620,0.33607,0.33681/                          3999.   
      DATA CSSALT/                                                      4000.   
     +        0.66858,0.50298,0.60372,0.65282,0.66694,0.67041,0.66666,  4001.   
     +        0.62258,0.52248,0.44732,0.32878,0.66866,0.66680,0.66404,  4002.   
     +        0.66252,0.66281,0.66265,0.66244,0.66232,0.66223,0.66226,  4003.   
     +        0.67338,0.67406,0.67410,0.67408/                          4004.   
C                                                                       4005.   
      DATA QDUST1/                                                      4006.   
     +        0.60958,0.65996,0.59890,0.73030,0.64827,0.55835,0.48157,  4007.   
     +        0.34847,0.23144,0.18097,0.13460,0.59012,0.47533,0.39938,  4008.   
     +        0.36575,0.35808,0.33834,0.31587,0.30849,0.30369,0.29821,  4009.   
     +        0.91360,1.14613,1.16193,1.16619/                          4010.   
      DATA SDUST1/                                                      4011.   
     +        0.32015,0.60541,0.49800,0.59591,0.46651,0.39745,0.34242,  4012.   
     +        0.23468,0.13039,0.08473,0.04350,0.29084,0.23940,0.16410,  4013.   
     +        0.13070,0.13267,0.12095,0.10691,0.10167,0.09788,0.09578,  4014.   
     +        0.39128,0.54469,0.55555,0.55942/                          4015.   
      DATA CDUST1/                                                      4016.   
     +        0.50425,0.49645,0.57736,0.63615,0.63373,0.66224,0.67205,  4017.   
     +        0.67034,0.65137,0.61767,0.53600,0.49640,0.47921,0.43825,  4018.   
     +        0.40760,0.41364,0.41120,0.40706,0.40418,0.40149,0.40315,  4019.   
     +        0.47280,0.39308,0.38801,0.38670/                          4020.   
C                                                                       4021.   
      DATA QDUST2/                                                      4022.   
     +        0.95483,0.71515,0.77676,0.91847,0.93699,0.89565,0.82979,  4023.   
     +        0.74871,0.70959,0.69272,0.68748,0.94632,0.90846,0.85600,  4024.   
     +        0.83350,0.83544,0.82317,0.80807,0.80270,0.79879,0.79577,  4025.   
     +        1.02427,1.12417,1.13054,1.13169/                          4026.   
      DATA SDUST2/                                                      4027.   
     +        0.49885,0.58157,0.55165,0.64038,0.59140,0.55222,0.50136,  4028.   
     +        0.42019,0.36087,0.33502,0.31667,0.49026,0.47989,0.42207,  4029.   
     +        0.39751,0.40487,0.39774,0.38819,0.38426,0.38107,0.38027,  4030.   
     +        0.49780,0.59147,0.59817,0.60013/                          4031.   
      DATA CDUST2/                                                      4032.   
     +        0.74352,0.54594,0.68229,0.72513,0.73598,0.75710,0.75041,  4033.   
     +        0.70723,0.65024,0.61702,0.58021,0.74556,0.74741,0.75647,  4034.   
     +        0.76384,0.76647,0.77599,0.78746,0.79136,0.79400,0.79700,  4035.   
     +        0.71874,0.62817,0.62224,0.62062/                          4036.   
C                                                                       4037.   
      DATA QCARB1/                                                      4038.   
     +        0.44718,0.51882,0.26055,0.20526,0.19295,0.18655,0.17520,  4039.   
     +        0.11120,0.06749,0.04893,0.03537,0.32912,0.25261,0.24973,  4040.   
     +        0.23947,0.22883,0.20424,0.18781,0.18400,0.18032,0.17370,  4041.   
     +        0.57200,0.64430,0.65267,0.65790/                          4042.   
      DATA SCARB1/                                                      4043.   
     +        0.17857,0.12659,0.06506,0.05088,0.05317,0.05712,0.05562,  4044.   
     +        0.03310,0.01705,0.01009,0.00493,0.13908,0.08683,0.06332,  4045.   
     +        0.06114,0.06260,0.05755,0.05319,0.05155,0.05032,0.04981,  4046.   
     +        0.19594,0.21003,0.20967,0.20853/                          4047.   
      DATA CCARB1/                                                      4048.   
     +        0.40490,0.48729,0.43960,0.40824,0.46236,0.51422,0.53366,  4049.   
     +        0.53211,0.51283,0.46211,0.32882,0.40923,0.35984,0.30817,  4050.   
     +        0.30468,0.31306,0.31215,0.30857,0.30555,0.30388,0.30644,  4051.   
     +        0.43102,0.40748,0.40436,0.40208/                          4052.   
C                                                                       4053.   
      DATA QCARB2/                                                      4054.   
     +        0.09591,0.22971,0.21603,0.21745,0.17928,0.17061,0.15202,  4055.   
     +        0.10846,0.06721,0.04817,0.03076,0.09456,0.08428,0.07093,  4056.   
     +        0.06589,0.06737,0.06766,0.06782,0.06771,0.06754,0.06792,  4057.   
     +        0.12455,0.12130,0.12121,0.12155/                          4058.   
      DATA SCARB2/                                                      4059.   
     +        0.00748,0.06133,0.05031,0.04978,0.03714,0.03448,0.03065,  4060.   
     +        0.02099,0.01137,0.00688,0.00291,0.00728,0.00544,0.00350,  4061.   
     +        0.00276,0.00291,0.00290,0.00288,0.00285,0.00282,0.00286,  4062.   
     +        0.01420,0.01327,0.01324,0.01332/                          4063.   
      DATA CCARB2/                                                      4064.   
     +        0.14117,0.25269,0.27090,0.30506,0.29845,0.28974,0.28880,  4065.   
     +        0.28843,0.28603,0.28395,0.29112,0.14128,0.12741,0.11121,  4066.   
     +        0.09892,0.09935,0.09786,0.09604,0.09517,0.09448,0.09466,  4067.   
     +        0.18297,0.17686,0.17658,0.17696/                          4068.   
C                                                                       4069.   
      DATA QWATER/                                                      4070.   
     +        0.82334,0.89509,1.13254,1.20762,1.24075,1.18580,1.07585,  4071.   
     +        0.95283,0.89542,0.86914,0.85864,0.87834,0.94021,1.03878,  4072.   
     +        1.07876,1.06927,1.06987,1.07153,1.07327,1.07505,1.07280,  4073.   
     +        1.20709,1.20194,1.20383,1.20978/                          4074.   
      DATA SWATER/                                                      4075.   
     +        0.34695,0.68566,0.86748,0.89010,0.83121,0.75556,0.65338,  4076.   
     +        0.51441,0.40925,0.36469,0.31873,0.39396,0.39368,0.43707,  4077.   
     +        0.45625,0.44997,0.45039,0.45146,0.45251,0.45357,0.45227,  4078.   
     +        0.85537,0.85478,0.85718,0.86370/                          4079.   
      DATA CWATER/                                                      4080.   
     +        0.91848,0.65450,0.79206,0.82335,0.83709,0.84869,0.84338,  4081.   
     +        0.77907,0.68419,0.62521,0.54076,0.91355,0.89224,0.85667,  4082.   
     +        0.84557,0.85029,0.85229,0.85399,0.85411,0.85389,0.85524,  4083.   
     +        0.91095,0.91472,0.91488,0.91467/                          4084.   
C                                                                       4085.   
      DATA QICE25/                                                      4086.   
     +        1.15210,0.81551,0.98885,1.10325,1.17652,1.14217,1.07777,  4087.   
     +        1.08252,1.14496,1.16939,1.22006,1.16194,1.16781,1.19342,  4088.   
     +        1.20279,1.19736,1.19435,1.19146,1.19097,1.19095,1.18924,  4089.   
     +        1.19321,1.21794,1.21959,1.21942/                          4090.   
      DATA SICE25/                                                      4091.   
     +        0.57392,0.45452,0.57278,0.68806,0.74580,0.69171,0.64662,  4092.   
     +        0.62884,0.64120,0.64892,0.66105,0.59403,0.60241,0.67853,  4093.   
     +        0.70399,0.68299,0.66547,0.64731,0.64301,0.64122,0.63321,  4094.   
     +        0.71867,0.77122,0.77524,0.77622/                          4095.   
      DATA CICE25/                                                      4096.   
     +        0.93634,0.72920,0.86084,0.88431,0.87489,0.88472,0.86613,  4097.   
     +        0.82078,0.79850,0.79041,0.78539,0.93377,0.91036,0.85751,  4098.   
     +        0.84228,0.85220,0.86089,0.87036,0.87263,0.87355,0.87810,  4099.   
     +        0.94697,0.94840,0.94812,0.94714/                          4100.   
C                                                                       4101.   
C-----------------------------------------------------------------------4102.   
C    SOLAR RADIATION 6 K-INTERVAL MERGED AEROSOL DATA FOR QEXT,QSCA,COSB4103.   
C-----------------------------------------------------------------------4104.   
C                                                                       4105.   
      DATA XACID1/    0.05776,0.10033,0.19099,0.36614,0.55931,1.04703/  4106.   
      DATA YACID1/    0.01880,0.09956,0.19090,0.36613,0.55931,1.04703/  4107.   
      DATA ZACID1/    0.36054,0.51871,0.57276,0.62068,0.65273,0.68988/  4108.   
C                                                                       4109.   
      DATA XACID2/0.13360,0.33875,0.51498,0.68359,0.79939,0.94494/      4110.   
      DATA YACID2/0.07420,0.33691,0.51483,0.68358,0.79939,0.94494/      4111.   
C$    DATA ZACID2/0.40248,0.62259,0.68524,0.71328,0.71195,0.72894/      4112.   
      DATA ZACID2/0.39821,0.54835,0.60846,0.63637,0.63503,0.65221/      4112.1  
C                                                                       4113.   
      DATA XSLFT1/    0.31035,0.44757,0.54238,0.66756,0.78260,1.04454/  4114.   
      DATA YSLFT1/    0.24589,0.44490,0.54224,0.66755,0.78260,1.04454/  4115.   
      DATA ZSLFT1/    0.70591,0.67557,0.66832,0.66438,0.66199,0.66008/  4116.   
C                                                                       4117.   
      DATA XSLFT2/    0.60959,0.74888,0.81124,0.87560,0.92632,1.00936/  4118.   
      DATA YSLFT2/    0.50477,0.74262,0.81090,0.87556,0.92631,1.00935/  4119.   
      DATA ZSLFT2/    0.74067,0.70281,0.69748,0.69922,0.70070,0.70754/  4120.   
C                                                                       4121.   
      DATA XBSLT1/    0.30419,0.46195,0.54908,0.66403,0.77732,1.02644/  4122.   
      DATA YBSLT1/    0.28732,0.44765,0.53358,0.64786,0.76063,1.00769/  4123.   
      DATA ZBSLT1/    0.67768,0.66588,0.66785,0.66932,0.66671,0.66818/  4124.   
C                                                                       4125.   
      DATA XBSLT2/    0.62145,0.76377,0.81783,0.87743,0.92782,1.00765/  4126.   
      DATA YBSLT2/    0.58466,0.73120,0.78367,0.84258,0.89259,0.96944/  4127.   
      DATA ZBSLT2/    0.70368,0.69767,0.70313,0.70847,0.70983,0.71935/  4128.   
C                                                                       4129.   
      DATA XSSALT/    0.64091,0.78294,0.83066,0.87490,0.92554,1.00414/  4130.   
      DATA YSSALT/    0.64091,0.78294,0.83066,0.87490,0.92554,1.00413/  4131.   
      DATA ZSSALT/    0.67233,0.68272,0.68718,0.69084,0.69334,0.69627/  4132.   
C                                                                       4133.   
      DATA XDUST1/    1.17571,1.20282,1.13894,1.08190,1.04572,0.99864/  4134.   
      DATA YDUST1/    1.04642,1.12320,1.04442,0.97057,0.93288,0.78720/  4135.   
      DATA ZDUST1/    0.72235,0.68164,0.69516,0.72361,0.74315,0.80409/  4136.   
C                                                                       4137.   
      DATA XDUST2/    1.09335,1.12888,1.09512,1.05217,1.02411,1.00081/  4138.   
      DATA YDUST2/    0.83740,0.93590,0.88162,0.81721,0.78602,0.68767/  4139.   
      DATA ZDUST2/    0.78776,0.76447,0.77511,0.79364,0.80840,0.85594/  4140.   
C                                                                       4141.   
      DATA XCARB1/0.74444,1.11851,1.14599,1.09902,1.05179,1.00292/      4142.   
      DATA YCARB1/0.53412,1.11290,1.14544,1.09899,1.05179,1.00292/      4143.   
C$    DATA ZCARB1/0.75767,0.74553,0.72950,0.71977,0.71968,0.74073/      4144.   
      DATA ZCARB1/0.71248,0.66984,0.65284,0.64292,0.64282,0.66426/      4144.1  
C                                                                       4145.   
      DATA XCARB2/    0.54418,0.82500,0.91922,0.97919,1.00345,0.99476/  4146.   
      DATA YCARB2/    0.19636,0.34820,0.40558,0.44719,0.46860,0.48132/  4147.   
      DATA ZCARB2/    0.45878,0.59691,0.65112,0.70444,0.74341,0.79820/  4148.   
C                                                                       4149.   
      DATA XWATER/    1.10372,1.05381,1.03792,1.02265,1.01285,0.99989/  4150.   
      DATA YWATER/    0.84758,1.03190,1.02896,1.02226,1.01282,0.99988/  4151.   
      DATA ZWATER/    0.87621,0.84587,0.84884,0.85323,0.85888,0.86321/  4152.   
C                                                                       4153.   
      DATA XICE25/    1.05394,1.02884,1.02030,1.01257,1.00706,0.99981/  4154.   
      DATA YICE25/    0.75677,0.96035,1.00797,1.01184,1.00702,0.99981/  4155.   
      DATA ZICE25/    0.92708,0.88645,0.87975,0.87906,0.87391,0.87623/  4156.   
C                                                                       4157.   
C-----------------------------------------------------------------------4158.   
C THERMAL RADIATION 25 K-INTERVAL MERGED CLOUD & SURFACE ALBEDO DATA    4159.   
C-----------------------------------------------------------------------4160.   
      DATA AGSIDV/                                                      4161.   
     S        0.01407,0.01653,0.03230,0.08764,0.10055,0.09095,0.08892,  4162.   
     S        0.07985,0.06411,0.05926,0.05398,0.01576,0.02449,0.05091,  4163.   
     S        0.05680,0.05325,0.04652,0.03809,0.03574,0.03451,0.03082,  4164.   
     S        0.01757,0.02022,0.02059,0.02082,                          4165.   
     I        0.01407,0.01653,0.03230,0.08764,0.10055,0.09095,0.08892,  4166.   
     I        0.07985,0.06411,0.05926,0.05398,0.01576,0.02449,0.05091,  4167.   
     I        0.05680,0.05325,0.04652,0.03809,0.03574,0.03451,0.03082,  4168.   
     I        0.01757,0.02022,0.02059,0.02082,                          4169.   
     D     0.04500,0.10414,0.06739,0.08448,0.08516,0.06283,0.05230,     4170.   
     D     0.03382,0.01901,0.01542,0.01178,0.05142,0.04835,0.05505,     4171.   
     D     0.05600,0.05310,0.04603,0.03731,0.03472,0.03328,0.03000,     4172.   
     D     0.16159,0.17592,0.17812,0.17927,                             4173.   
     V     25*0.0/                                                      4174.   
      DATA AOCEAN/                                                      4175.   
     +        0.04000,0.05965,0.06124,0.08339,0.09235,0.09510,0.09908,  4176.   
     +        0.11117,0.12263,0.12577,0.12931,0.04700,0.06894,0.08970,  4177.   
     +        0.09574,0.09565,0.09619,0.09672,0.09703,0.09723,0.09700,  4178.   
     +        0.04645,0.04487,0.04482,0.04493/                          4179.   
C                                                                       4180.   
      DATA CLDALB/                                                      4181.   
     +        0.01332,0.08190,0.07036,0.05082,0.04486,0.04673,0.04770,  4182.   
     +        0.05130,0.05240,0.05251,0.05259,0.01558,0.01763,0.02410,  4183.   
     +        0.02571,0.02514,0.02448,0.02366,0.02347,0.02340,0.02294,  4184.   
     +        0.04566,0.04499,0.04518,0.04544,                          4185.   
     +        0.01407,0.01653,0.03230,0.08764,0.10055,0.09095,0.08892,  4186.   
     +        0.07985,0.06411,0.05926,0.05398,0.01576,0.02449,0.05091,  4187.   
     +        0.05680,0.05325,0.04652,0.03809,0.03574,0.03451,0.03082,  4188.   
     +        0.01757,0.02022,0.02059,0.02082/                          4189.   
C                                                                       4190.   
      DATA ASNALB/0.600,0.350,13*0.0/                                   4191.   
C&&   DATA ASNALB/0.550,0.300,13*0.0/
C                                                                       4192.   
C&&   DATA AOIALB/0.550,0.300,13*0.0/                                   4193.   
      DATA AOIALB/0.600,0.350,13*0.0/
C                                                                       4194.   
      DATA ALIALB/0.600,0.350,13*0.0/                                   4195.   
C                                                                       4196.   
C-----------------------------------------------------------------------4197.   
C            TRACE GAS VERTICAL DISTRIBUTION & 1958 MEAN CONCENTRATION  4198.   
C-----------------------------------------------------------------------4199.   
C                                                                       4200.   
      DATA CMANO2/                                                      4201.   
     1  8.66E-06,5.15E-06,2.85E-06,1.50E-06,9.89E-07,6.91E-07,7.17E-07, 4202.   
     2  8.96E-07,3.67E-06,4.85E-06,5.82E-06,6.72E-06,7.77E-06,8.63E-06, 4203.   
     3  8.77E-06,8.14E-06,6.91E-06,5.45E-06,4.00E-06,2.67E-06,1.60E-06, 4204.   
     4  8.36E-07,3.81E-07,1.58E-07,6.35E-08,2.57E-08,1.03E-08,4.18E-09, 4205.   
     5  1.66E-09,6.57E-10,2.58E-10,1.02E-10,4.11E-11,1.71E-11,7.73E-12, 4206.   
     6  9.07E-12,4.63E-12,2.66E-12,1.73E-12,1.28E-12,1.02E-12,1.00E-30/ 4207.   
C                                                                       4208.   
C                                                                       4209.   
C     GAS  NUMBER  1     2    3       4  5     6     7       8      9   4210.   
C                 H2O   CO2  O3      O2 NO2   N2O   CH4   CCL3F1 CCL2F2 4211.   
C     DATA FULGAS/1.0,  1.0,1.0,    1.0,1.0,  1.0,  1.0,    1.0,    1.0/4212.   
c     DATA PPMV58/0.0,315.0,0.0,210000.,0.0,0.295,1.400,8.00E-6,25.0E-6/4213.   
      DATA PPMV58/0.0,  0.0,0.0,210000.,0.0,4*0.0/
C$    DATA Z0/    0.0,  0.0,0.0,    0.0,0.0, 15.0, 10.0,   12.0,   12.0/4214.   
      DATA Z0/    0.0,  0.0,0.0,    0.0,0.0,915.0,910.0,   12.0,   12.0/4215.   
      DATA ZH/    8.0,  8.0,8.0,    8.0,8.0, 10.0, 30.0,    3.0,    3.0/4216.   
C                                                                       4217.   
C-----------------------------------------------------------------------4218.   
C                     TRACE GAS ABSORPTION COEFFICIENTS FOR  F11 & F12  4219.   
C-----------------------------------------------------------------------4220.   
C                                                                       4221.   
      DIMENSION F11PCM(25),F12PCM(25)                                   4222.   
      EQUIVALENCE (TRACEG(1,1),F11PCM(1)),(TRACEG(1,2),F12PCM(1))       4223.   
C                                                                       4224.   
C                                                                       4225.   
      DATA F11PCM/                                                      4226.   
     +        13.6000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,  4227.   
     +         0.0000, 0.0000, 0.0000, 0.0000,11.9504, 2.5138, 0.5054,  4228.   
     +         0.1086, 0.0308, 0.0178, 0.0054, 0.0000, 0.0000, 0.0000,  4229.   
     +         2.5220, 1.1731, 0.8627, 0.7445/                          4230.   
C                                                                       4231.   
      DATA F12PCM/                                                      4232.   
     +         5.4900, 1.3339, 0.7739, 0.1304, 0.0286, 0.0051, 0.0000,  4233.   
     +         0.0000, 0.0000, 0.0000, 0.0000, 9.0745, 2.3577, 0.4135,  4234.   
     +         0.0575, 0.0000, 0.2507, 0.6215, 0.7262, 0.7972, 0.9150,  4235.   
     +        13.1663, 1.1564, 0.0388, 0.0082/                          4236.   
C                                                                       4236.11 
C     ------------------------------------------------------------------4236.12 
C     DECEMBER 4, 1991 UPDATE   PROVIDES FOR THE FOLLOWING IMPROVEMENTS:4236.13 
C     ------------------------------------------------------------------4236.14 
C     IF(NEWASZ.GT.0) ALL AEROSOL SOLAR ZENITH ANGLE DEPENDENCE IMPROVED4236.15 
C     IF(NEWAQA.GT.0) ALL AERSOL THERMAL CROSS-SECTIONS ARE Q-ABSORPTION4236.16 
C             (TRACER AEROSOLS ALREADY USE Q-ABSORPTION IN XRAD83XX)    4236.17 
C     IF(NEWCQA.GT.0) ALL CLOUDS THERMAL CROSS-SECTIONS ARE Q-ABSORPTION4236.18 
C     ------------------------------------------------------------------4236.21 
C                                                                       4236.22 
      EQUIVALENCE (ISPARE(1),NEWASZ)                                    4236.23 
      EQUIVALENCE (ISPARE(2),NEWAQA)                                    4236.24 
      EQUIVALENCE (ISPARE(3),NEWCQA)                                    4236.25 
C                                                                       4236.26 
      DATA NEWASZ/0/,  NEWAQA/0/,  NEWCQA/0/                            4236.27 
C                                                                       4236.28 
      END                                                               4237.   
      SUBROUTINE PHDATM(P,H,D,T,O,Q,S,OCM,WCM,NPHD,NATM)                4238.   
C                                                                       4239.   
C     ------------------------------------------------------------------4240.   
C     -------------     MCCLATCHY (1972) ATMOSPHERE DATA     -----------4241.   
C     ------------------------------------------------------------------4242.   
C                                                                       4243.   
C        INPUT DATA                                                     4244.   
C------------------                                                     4245.   
C                  NATM=0  GIVES ABREVIATED DATA FOR  STANDARD ATMOSPHER4246.   
C                                (INPUT: P OR H) (RETURNS: H OR P & D,T)4247.   
C                                                                       4248.   
C                  NATM=1  GIVES ATMOSPHERE DATA FOR  TROPICAL LATITUDES4249.   
C                  NATM=2  GIVES ATMOSPHERE DATA FOR  MIDLATITUDE SUMMER4250.   
C                  NATM=3  GIVES ATMOSPHERE DATA FOR  MIDLATITUDE WINTER4251.   
C                  NATM=4  GIVES ATMOSPHERE DATA FOR  SUBARCTIC SUMMER  4252.   
C                  NATM=5  GIVES ATMOSPHERE DATA FOR  SUBARCTIC WINTER  4253.   
C                  NATM=6  GIVES ATMOSPHERE DATA FOR  STANDARD ATMOSPHER4254.   
C                                                                       4255.   
C                  NPHD=1  RETURNS H,D,T,O,Q,S DATA FOR GIVEN PRESSURE P4256.   
C                  NPHD=2  RETURNS P,D,T,O,Q,S DATA FOR GIVEN   HEIGHT H4257.   
C                  NPHD=3  RETURNS P,H,T,O,Q,S DATA FOR GIVEN  DENSITY D4258.   
C                                                                       4259.   
C       OUTPUT DATA                                                     4260.   
C------------------                                                     4261.   
C                  P = PRESSURE IN MILLIBARS                            4262.   
C                  H = HEIGHT IN KILOMETERS                             4263.   
C                  D = DENSITY IN GRAMS/METER**3                        4264.   
C                  T = TEMPERATURE (ABSOLUTE)                           4265.   
C                  O = OZONE MIXING RATIO (GRAMS OZONE)/(GRAMS AIR)     4266.   
C                  Q = SPECIFIC HUMIDITY (GRAMS WATER VAPOR)/(GRAMS AIR)4267.   
C                  S = SATURATION RATIO (GRAMS WATER VAPOR)/(GRAMS AIR) 4268.   
C                  OCM = OZONE (CM-STP) ABOVE GIVEN HEIGHT              4269.   
C                  WCM = WATER VAPOR (CM-STP) ABOVE GIVEN HEIGHT        4270.   
C                                                                       4271.   
C           REMARKS                                                     4272.   
C------------------                                                     4273.   
C                  INPUT P,H,D PARAMETERS ARE NOT ALTERED               4274.   
C                  P,D INTERPOLATION IS EXPONENTIAL WITH HEIGHT         4275.   
C                  NO EXTRAPOLATION IS MADE OUTSIDE 0-100 KM INTERVAL   4276.   
C                  S  IS NOT COMPUTED ABOVE 40 KM (FORMULA NOT ACCURATE)4277.   
C                                                                       4278.   
C                  R = Q/S          GIVES RELATIVE HUMIDITY             4279.   
C                  W = Q/(1-Q)      GIVES WATER VAPOR MIXING RATIO      4280.   
C                  N = D*2.079E 16  GIVES NUMBER DENSITY PER CM**3      4281.   
C                                                                       4282.   
C                                                                       4283.   
C                                                                       4284.   
C                                                                       4285.   
C                                                                       4286.   
      DIMENSION    PRS1(33),PRS2(33),PRS3(33),PRS4(33),PRS5(33),PRS6(33)4287.   
     1            ,DNS1(33),DNS2(33),DNS3(33),DNS4(33),DNS5(33),DNS6(33)4288.   
     2            ,TMP1(33),TMP2(33),TMP3(33),TMP4(33),TMP5(33),TMP6(33)4289.   
     3            ,WVP1(33),WVP2(33),WVP3(33),WVP4(33),WVP5(33),WVP6(33)4290.   
     4            ,OZO1(33),OZO2(33),OZO3(33),OZO4(33),OZO5(33),OZO6(33)4291.   
      DIMENSION   PRES(33,6),DENS(33,6),TEMP(33,6),WVAP(33,6),OZON(33,6)4292.   
C                                                                       4293.   
      EQUIVALENCE                                                       4294.   
     +     (PRES(1,1),PRS1(1)),(DENS(1,1),DNS1(1)),(TEMP(1,1),TMP1(1))  4295.   
     +    ,(PRES(1,2),PRS2(1)),(DENS(1,2),DNS2(1)),(TEMP(1,2),TMP2(1))  4296.   
     +    ,(PRES(1,3),PRS3(1)),(DENS(1,3),DNS3(1)),(TEMP(1,3),TMP3(1))  4297.   
     +    ,(PRES(1,4),PRS4(1)),(DENS(1,4),DNS4(1)),(TEMP(1,4),TMP4(1))  4298.   
     +    ,(PRES(1,5),PRS5(1)),(DENS(1,5),DNS5(1)),(TEMP(1,5),TMP5(1))  4299.   
     +    ,(PRES(1,6),PRS6(1)),(DENS(1,6),DNS6(1)),(TEMP(1,6),TMP6(1))  4300.   
      EQUIVALENCE (WVAP(1,1),WVP1(1)),(OZON(1,1),OZO1(1))               4301.   
      EQUIVALENCE (WVAP(1,2),WVP2(1)),(OZON(1,2),OZO2(1))               4302.   
      EQUIVALENCE (WVAP(1,3),WVP3(1)),(OZON(1,3),OZO3(1))               4303.   
      EQUIVALENCE (WVAP(1,4),WVP4(1)),(OZON(1,4),OZO4(1))               4304.   
      EQUIVALENCE (WVAP(1,5),WVP5(1)),(OZON(1,5),OZO5(1))               4305.   
      EQUIVALENCE (WVAP(1,6),WVP6(1)),(OZON(1,6),OZO6(1))               4306.   
C                                                                       4307.   
C                                                                       4308.   
      DIMENSION HTKM(33)                                                4309.   
      DATA HTKM/1.0E-09, 1., 2., 3., 4., 5., 6., 7., 8., 9.,10.,11.     4310.   
     1         ,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.     4311.   
     2         ,25.,30.,35.,40.,45.,50.,70.,99.9/                       4312.   
C                                                                       4313.   
C                                                                       4314.   
C---------------------------------------------------------------------- 4315.   
C0000 GLOBAL   U.S. (1976) STANDARD ATMOSPHERE   P, T, GEO H  PARAMETERS4316.   
C---------------------------------------------------------------------- 4317.   
C                                                                       4318.   
      DIMENSION SPLB(8),STLB(8),SHLB(8),SDLB(8)                         4319.   
      DATA SPLB/1013.25,226.32,54.748,8.6801,1.109,.66938,.039564       4320.   
     +         ,3.7338E-03/                                             4321.   
      DATA STLB/288.15,216.65,216.65,228.65,270.65,270.65,214.65,186.87/4322.   
      DATA SHLB/0.0,11.0,20.0,32.0,47.0,51.0,71.0,84.852/               4323.   
      DATA SDLB/-6.5,0.0,1.0,2.8,0.0,-2.8,-2.0,0.0/                     4324.   
      DATA HPCON/34.16319/                                              4325.   
C                                                                       4326.   
C                                                                       4327.   
C-----------------------------------------------------------------------4328.   
C1111 TROPICAL LATITUDES      MCCLATCHY (1972) ATMOSPHERE DATA VS HEIGHT4329.   
C-----------------------------------------------------------------------4330.   
C                                                                       4331.   
      DATA PRS1/      1.013E 03,9.040E 02,8.050E 02,7.150E 02,6.330E 02,4332.   
     1      5.590E 02,4.920E 02,4.320E 02,3.780E 02,3.290E 02,2.860E 02,4333.   
     2      2.470E 02,2.130E 02,1.820E 02,1.560E 02,1.320E 02,1.110E 02,4334.   
     3      9.370E 01,7.890E 01,6.660E 01,5.650E 01,4.800E 01,4.090E 01,4335.   
     4      3.500E 01,3.000E 01,2.570E 01,1.220E 01,6.000E 00,3.050E 00,4336.   
     5      1.590E 00,8.540E-01,5.790E-02,3.000E-04/                    4337.   
      DATA DNS1/      1.167E 03,1.064E 03,9.689E 02,8.756E 02,7.951E 02,4338.   
     1      7.199E 02,6.501E 02,5.855E 02,5.258E 02,4.708E 02,4.202E 02,4339.   
     2      3.740E 02,3.316E 02,2.929E 02,2.578E 02,2.260E 02,1.972E 02,4340.   
     3      1.676E 02,1.382E 02,1.145E 02,9.515E 01,7.938E 01,6.645E 01,4341.   
     4      5.618E 01,4.763E 01,4.045E 01,1.831E 01,8.600E 00,4.181E 00,4342.   
     5      2.097E 00,1.101E 00,9.210E-02,5.000E-04/                    4343.   
      DATA TMP1/  300.0,294.0,288.0,284.0,277.0,270.0,264.0,257.0,250.0,4344.   
     1244.0,237.0,230.0,224.0,217.0,210.0,204.0,197.0,195.0,199.0,203.0,4345.   
     2207.0,211.0,215.0,217.0,219.0,221.0,232.0,243.0,254.0,265.0,270.0,4346.   
     3  219.0,210.0/                                                    4347.   
      DATA WVP1/1.9E 01,1.3E 01,9.3E 00,4.7E 00,2.2E 00,1.5E 00,8.5E-01,4348.   
     1  4.7E-01,2.5E-01,1.2E-01,5.0E-02,1.7E-02,6.0E-03,1.8E-03,1.0E-03,4349.   
     2  7.6E-04,6.4E-04,5.6E-04,5.0E-04,4.9E-04,4.5E-04,5.1E-04,5.1E-04,4350.   
     3  5.4E-04,6.0E-04,6.7E-04,3.6E-04,1.1E-04,4.3E-05,1.9E-05,6.3E-06,4351.   
     4  1.4E-07,1.0E-09/                                                4352.   
      DATA OZO1/5.6E-05,5.6E-05,5.4E-05,5.1E-05,4.7E-05,4.5E-05,4.3E-05,4353.   
     1  4.1E-05,3.9E-05,3.9E-05,3.9E-05,4.1E-05,4.3E-05,4.5E-05,4.5E-05,4354.   
     2  4.7E-05,4.7E-05,6.9E-05,9.0E-05,1.4E-04,1.9E-04,2.4E-04,2.8E-04,4355.   
     3  3.2E-04,3.4E-04,3.4E-04,2.4E-04,9.2E-05,4.1E-05,1.3E-05,4.3E-06,4356.   
     4  8.6E-08,4.3E-11/                                                4357.   
C                                                                       4358.   
C-----------------------------------------------------------------------4359.   
C2222 MIDLATITUDE SUMMER      MCCLATCHY (1972) ATMOSPHERE DATA VS HEIGHT4360.   
C-----------------------------------------------------------------------4361.   
C                                                                       4362.   
      DATA PRS2/      1.013E 03,9.020E 02,8.020E 02,7.100E 02,6.280E 02,4363.   
     1      5.540E 02,4.870E 02,4.260E 02,3.720E 02,3.240E 02,2.810E 02,4364.   
     2      2.430E 02,2.090E 02,1.790E 02,1.530E 02,1.300E 02,1.110E 02,4365.   
     3      9.500E 01,8.120E 01,6.950E 01,5.950E 01,5.100E 01,4.370E 01,4366.   
     4      3.760E 01,3.220E 01,2.770E 01,1.320E 01,6.520E 00,3.330E 00,4367.   
     5      1.760E 00,9.510E-01,6.710E-02,3.000E-04/                    4368.   
      DATA DNS2/      1.191E 03,1.080E 03,9.757E 02,8.846E 02,7.998E 02,4369.   
     1      7.211E 02,6.487E 02,5.830E 02,5.225E 02,4.669E 02,4.159E 02,4370.   
     2      3.693E 02,3.269E 02,2.882E 02,2.464E 02,2.104E 02,1.797E 02,4371.   
     3      1.535E 02,1.305E 02,1.110E 02,9.453E 01,8.056E 01,6.872E 01,4372.   
     4      5.867E 01,5.014E 01,4.288E 01,1.322E 01,6.519E 00,3.330E 00,4373.   
     5      1.757E 00,9.512E-01,6.706E-02,5.000E-04/                    4374.   
      DATA TMP2/  294.0,290.0,285.0,279.0,273.0,267.0,261.0,255.0,248.0,4375.   
     1242.0,235.0,229.0,222.0,216.0,216.0,216.0,216.0,216.0,216.0,217.0,4376.   
     2218.0,219.0,220.0,222.0,223.0,224.0,234.0,245.0,258.0,270.0,276.0,4377.   
     3  218.0,210.0/                                                    4378.   
      DATA WVP2/1.4E 01,9.3E 00,5.9E 00,3.3E 00,1.9E 00,1.0E 00,6.1E-01,4379.   
     1  3.7E-01,2.1E-01,1.2E-01,6.4E-02,2.2E-02,6.0E-03,1.8E-03,1.0E-03,4380.   
     2  7.6E-04,6.4E-04,5.6E-04,5.0E-04,4.9E-04,4.5E-04,5.1E-04,5.1E-04,4381.   
     3  5.4E-04,6.0E-04,6.7E-04,3.6E-04,1.1E-04,4.3E-05,1.9E-05,6.3E-06,4382.   
     4  1.4E-07,1.0E-09/                                                4383.   
      DATA OZO2/6.0E-05,6.0E-05,6.0E-05,6.2E-05,6.4E-05,6.6E-05,6.9E-05,4384.   
     1  7.5E-05,7.9E-05,8.6E-05,9.0E-05,1.1E-04,1.2E-04,1.5E-04,1.8E-04,4385.   
     2  1.9E-04,2.1E-04,2.4E-04,2.8E-04,3.2E-04,3.4E-04,3.6E-04,3.6E-04,4386.   
     3  3.4E-04,3.2E-04,3.0E-04,2.0E-04,9.2E-05,4.1E-05,1.3E-05,4.3E-06,4387.   
     4  8.6E-08,4.3E-11/                                                4388.   
C                                                                       4389.   
C-----------------------------------------------------------------------4390.   
C3333 MIDLATITUDE WINTER      MCCLATCHY (1972) ATMOSPHERE DATA VS HEIGHT4391.   
C-----------------------------------------------------------------------4392.   
C                                                                       4393.   
      DATA PRS3/      1.018E 03,8.973E 02,7.897E 02,6.938E 02,6.081E 02,4394.   
     1      5.313E 02,4.627E 02,4.016E 02,3.473E 02,2.992E 02,2.568E 02,4395.   
     2      2.199E 02,1.882E 02,1.610E 02,1.378E 02,1.178E 02,1.007E 02,4396.   
     3      8.610E 01,7.350E 01,6.280E 01,5.370E 01,4.580E 01,3.910E 01,4397.   
     4      3.340E 01,2.860E 01,2.430E 01,1.110E 01,5.180E 00,2.530E 00,4398.   
     5      1.290E 00,6.820E-01,4.670E-02,3.000E-04/                    4399.   
      DATA DNS3/      1.301E 03,1.162E 03,1.037E 03,9.230E 02,8.282E 02,4400.   
     1      7.411E 02,6.614E 02,5.886E 02,5.222E 02,4.619E 02,4.072E 02,4401.   
     2      3.496E 02,2.999E 02,2.572E 02,2.206E 02,1.890E 02,1.620E 02,4402.   
     3      1.388E 02,1.188E 02,1.017E 02,8.690E 01,7.421E 01,6.338E 01,4403.   
     4      5.415E 01,4.624E 01,3.950E 01,1.783E 01,7.924E 00,3.625E 00,4404.   
     5      1.741E 00,8.954E-01,7.051E-02,5.000E-04/                    4405.   
      DATA TMP3/  272.2,268.7,265.2,261.7,255.7,249.7,243.7,237.7,231.7,4406.   
     1225.7,219.7,219.2,218.7,218.2,217.7,217.2,216.7,216.2,215.7,215.2,4407.   
     2215.2,215.2,215.2,215.2,215.2,215.2,217.4,227.8,243.2,258.5,265.7,4408.   
     3  230.7,210.2/                                                    4409.   
      DATA WVP3/3.5E 00,2.5E 00,1.8E 00,1.2E 00,6.6E-01,3.8E-01,2.1E-01,4410.   
     1  8.5E-02,3.5E-02,1.6E-02,7.5E-03,6.9E-03,6.0E-03,1.8E-03,1.0E-03,4411.   
     2  7.6E-04,6.4E-04,5.6E-04,5.0E-04,4.9E-04,4.5E-04,5.1E-04,5.1E-04,4412.   
     3  5.4E-04,6.0E-04,6.7E-04,3.6E-04,1.1E-04,4.3E-05,1.9E-05,6.3E-06,4413.   
     4  1.4E-07,1.0E-09/                                                4414.   
      DATA OZO3/6.0E-05,5.4E-05,4.9E-05,4.9E-05,4.9E-05,5.8E-05,6.4E-05,4415.   
     1  7.7E-05,9.0E-05,1.2E-04,1.6E-04,2.1E-04,2.6E-04,3.0E-04,3.2E-04,4416.   
     2  3.4E-04,3.6E-04,3.9E-04,4.1E-04,4.3E-04,4.5E-04,4.3E-04,4.3E-04,4417.   
     3  3.9E-04,3.6E-04,3.4E-04,1.9E-04,9.2E-05,4.1E-05,1.3E-05,4.3E-06,4418.   
     4  8.6E-08,4.3E-11/                                                4419.   
C                                                                       4420.   
C-----------------------------------------------------------------------4421.   
C4444 SUBARCTIC SUMMER        MCCLATCHY (1972) ATMOSPHERE DATA VS HEIGHT4422.   
C-----------------------------------------------------------------------4423.   
C                                                                       4424.   
      DATA PRS4/      1.010E 03,8.960E 02,7.929E 02,7.000E 02,6.160E 02,4425.   
     1      5.410E 02,4.730E 02,4.130E 02,3.590E 02,3.107E 02,2.677E 02,4426.   
     2      2.300E 02,1.977E 02,1.700E 02,1.460E 02,1.250E 02,1.080E 02,4427.   
     3      9.280E 01,7.980E 01,6.860E 01,5.890E 01,5.070E 01,4.360E 01,4428.   
     4      3.750E 01,3.227E 01,2.780E 01,1.340E 01,6.610E 00,3.400E 00,4429.   
     5      1.810E 00,9.870E-01,7.070E-02,3.000E-04/                    4430.   
      DATA DNS4/      1.220E 03,1.110E 03,9.971E 02,8.985E 02,8.077E 02,4431.   
     1      7.244E 02,6.519E 02,5.849E 02,5.231E 02,4.663E 02,4.142E 02,4432.   
     2      3.559E 02,3.059E 02,2.630E 02,2.260E 02,1.943E 02,1.671E 02,4433.   
     3      1.436E 02,1.235E 02,1.062E 02,9.128E 01,7.849E 01,6.750E 01,4434.   
     4      5.805E 01,4.963E 01,4.247E 01,1.338E 01,6.614E 00,3.404E 00,4435.   
     5      1.817E 00,9.868E-01,7.071E-02,5.000E-04/                    4436.   
      DATA TMP4/  287.0,282.0,276.0,271.0,266.0,260.0,253.0,246.0,239.0,4437.   
     1232.0,225.0,225.0,225.0,225.0,225.0,225.0,225.0,225.0,225.0,225.0,4438.   
     2225.0,225.0,225.0,225.0,226.0,228.0,235.0,247.0,262.0,274.0,277.0,4439.   
     3  216.0,210.0/                                                    4440.   
      DATA WVP4/9.1E 00,6.0E 00,4.2E 00,2.7E 00,1.7E 00,1.0E 00,5.4E-01,4441.   
     1  2.9E-01,1.3E-02,4.2E-02,1.5E-02,9.4E-03,6.0E-03,1.8E-03,1.0E-03,4442.   
     2  7.6E-04,6.4E-04,5.6E-04,5.0E-04,4.9E-04,4.5E-04,5.1E-04,5.1E-04,4443.   
     3  5.4E-04,6.0E-04,6.7E-04,3.6E-04,1.1E-04,4.3E-05,1.9E-05,6.3E-06,4444.   
     4  1.4E-07,1.0E-09/                                                4445.   
      DATA OZO4/4.9E-05,5.4E-05,5.6E-05,5.8E-05,6.0E-05,6.4E-05,7.1E-05,4446.   
     1  7.5E-05,7.9E-05,1.1E-04,1.3E-04,1.8E-04,2.1E-04,2.6E-04,2.8E-04,4447.   
     2  3.2E-04,3.4E-04,3.9E-04,4.1E-04,4.1E-04,3.9E-04,3.6E-04,3.2E-04,4448.   
     3  3.0E-04,2.8E-04,2.6E-04,1.4E-04,9.2E-05,4.1E-05,1.3E-05,4.3E-06,4449.   
     4  8.6E-08,4.3E-11/                                                4450.   
C                                                                       4451.   
C-----------------------------------------------------------------------4452.   
C5555 SUBARCTIC WINTER        MCCLATCHY (1972) ATMOSPHERE DATA VS HEIGHT4453.   
C-----------------------------------------------------------------------4454.   
C                                                                       4455.   
      DATA PRS5/      1.013E 03,8.878E 02,7.775E 02,6.798E 02,5.932E 02,4456.   
     1      5.158E 02,4.467E 02,3.853E 02,3.308E 02,2.829E 02,2.418E 02,4457.   
     2      2.067E 02,1.766E 02,1.510E 02,1.291E 02,1.103E 02,9.431E 01,4458.   
     3      8.058E 01,6.882E 01,5.875E 01,5.014E 01,4.277E 01,3.647E 01,4459.   
     4      3.109E 01,2.649E 01,2.256E 01,1.020E 01,4.701E 00,2.243E 00,4460.   
     5      1.113E 00,5.719E-01,4.016E-02,3.000E-04/                    4461.   
      DATA DNS5/      1.372E 03,1.193E 03,1.058E 03,9.366E 02,8.339E 02,4462.   
     1      7.457E 02,6.646E 02,5.904E 02,5.226E 02,4.538E 02,3.879E 02,4463.   
     2      3.315E 02,2.834E 02,2.422E 02,2.071E 02,1.770E 02,1.517E 02,4464.   
     3      1.300E 02,1.113E 02,9.529E 01,8.155E 01,6.976E 01,5.966E 01,4465.   
     4      5.100E 01,4.358E 01,3.722E 01,1.645E 01,7.368E 00,3.330E 00,4466.   
     5      1.569E 00,7.682E-01,5.695E-02,5.000E-04/                    4467.   
      DATA TMP5/  257.1,259.1,255.9,252.7,247.7,240.9,234.1,227.3,220.6,4468.   
     1217.2,217.2,217.2,217.2,217.2,217.2,217.2,216.6,216.0,215.4,214.8,4469.   
     2214.1,213.6,213.0,212.4,211.8,211.2,216.0,222.2,234.7,247.0,259.3,4470.   
     3  245.7,210.0/                                                    4471.   
      DATA WVP5/1.2E 00,1.2E 00,9.4E-01,6.8E-01,4.1E-01,2.0E-01,9.8E-02,4472.   
     1  5.4E-02,1.1E-02,8.4E-03,5.5E-03,3.8E-03,2.6E-03,1.8E-03,1.0E-03,4473.   
     2  7.6E-04,6.4E-04,5.6E-04,5.0E-04,4.9E-04,4.5E-04,5.1E-04,5.1E-04,4474.   
     3  5.4E-04,6.0E-04,6.7E-04,3.6E-04,1.1E-04,4.3E-05,1.9E-05,6.3E-06,4475.   
     4  1.4E-07,1.0E-09/                                                4476.   
      DATA OZO5/4.1E-05,4.1E-05,4.1E-05,4.3E-05,4.5E-05,4.7E-05,4.9E-05,4477.   
     1  7.1E-05,9.0E-05,1.6E-04,2.4E-04,3.2E-04,4.3E-04,4.7E-04,4.9E-04,4478.   
     2  5.6E-04,6.2E-04,6.2E-04,6.2E-04,6.0E-04,5.6E-04,5.1E-04,4.7E-04,4479.   
     3  4.3E-04,3.6E-04,3.2E-04,1.5E-04,9.2E-05,4.1E-05,1.3E-05,4.3E-06,4480.   
     4  8.6E-08,4.3E-11/                                                4481.   
C                                                                       4482.   
C---------------------------------------------------------------------- 4483.   
C6666 GLOBAL   U.S. (1976) STANDARD ATMOSPHERE   P, T, GEO H  PARAMETERS4484.   
C---------------------------------------------------------------------- 4485.   
C                                                                       4486.   
      DATA PRS6/    1.01325E+03,8.987E+02,7.950E+02,7.011E+02,6.164E+02,4487.   
     1      5.402E+02,4.718E+02,4.106E+02,3.560E+02,3.074E+02,2.644E+02,4488.   
     2      2.263E+02,1.933E+02,1.651E+02,1.410E+02,1.204E+02,1.029E+02,4489.   
     3      8.787E+01,7.505E+01,6.410E+01,5.475E+01,4.678E+01,4.000E+01,4490.   
     4      3.422E+01,2.931E+01,2.511E+01,1.172E+01,5.589E+00,2.775E+00,4491.   
     5      1.431E+00,7.594E-01,4.634E-02,2.384E-04/                    4492.   
      DATA DNS6/      1.225E+03,1.112E+03,1.006E+03,9.091E+02,8.191E+02,4493.   
     1      7.361E+02,6.597E+02,5.895E+02,5.252E+02,4.663E+02,4.127E+02,4494.   
     2      3.639E+02,3.108E+02,2.655E+02,2.268E+02,1.937E+02,1.654E+02,4495.   
     3      1.413E+02,1.207E+02,1.031E+02,8.803E+01,7.487E+01,6.373E+01,4496.   
     4      5.428E+01,4.627E+01,3.947E+01,1.801E+01,8.214E+00,3.851E+00,4497.   
     5      1.881E+00,9.775E-01,7.424E-02,4.445E-04/                    4498.   
      DATA TMP6/                                                        4499.   
     1         288.150,281.650,275.150,268.650,262.150,255.650,249.150, 4500.   
     2         242.650,236.150,229.650,223.150,216.650,216.650,216.650, 4501.   
     3         216.650,216.650,216.650,216.650,216.650,216.650,216.650, 4502.   
     4         217.650,218.650,219.650,220.650,221.650,226.650,237.050, 4503.   
     5         251.050,265.050,270.650,217.450,186.870/                 4504.   
      DATA WVP6/      1.083E+01,6.323E+00,3.612E+00,2.015E+00,1.095E+00,4505.   
     1      5.786E-01,2.965E-01,1.469E-01,7.021E-02,3.226E-02,1.419E-02,4506.   
     2      5.956E-03,5.002E-03,4.186E-03,3.490E-03,2.896E-03,2.388E-03,4507.   
     3      1.954E-03,1.583E-03,1.267E-03,9.967E-04,8.557E-04,7.104E-04,4508.   
     4      5.600E-04,4.037E-04,2.406E-04,5.404E-05,2.464E-05,1.155E-05,4509.   
     5      5.644E-06,2.932E-06,2.227E-07,1.334E-09/                    4510.   
      DATA OZO6/      7.526E-05,3.781E-05,6.203E-05,3.417E-05,5.694E-05,4511.   
     1      3.759E-05,5.970E-05,4.841E-05,7.102E-05,6.784E-05,9.237E-05,4512.   
     2      9.768E-05,1.251E-04,1.399E-04,1.715E-04,1.946E-04,2.300E-04,4513.   
     3      2.585E-04,2.943E-04,3.224E-04,3.519E-04,3.714E-04,3.868E-04,4514.   
     4      3.904E-04,3.872E-04,3.728E-04,2.344E-04,9.932E-05,3.677E-05,4515.   
     5      1.227E-05,4.324E-06,5.294E-08,1.262E-10/                    4516.   
C                                                                       4517.   
C                                                                       4518.   
      IF(NATM.GT.0) GO TO 200                                           4519.   
      O=1.E-10                                                          4520.   
      Q=1.E-10                                                          4521.   
      S=1.E-10                                                          4522.   
      OCM=1.E-10                                                        4523.   
      WCM=1.E-10                                                        4524.   
      IF(NPHD.LT.2) GO TO 150                                           4525.   
      DO 110 N=2,8                                                      4526.   
      IF(H.LT.SHLB(N)) GO TO 120                                        4527.   
 110  CONTINUE                                                          4528.   
      N=9                                                               4529.   
 120  N=N-1                                                             4530.   
      IF(ABS(SDLB(N)).LT.1.E-04) GO TO 130                              4531.   
      P=SPLB(N)*(1.+SDLB(N)/STLB(N)*(H-SHLB(N)))**(-HPCON/SDLB(N))      4532.   
      GO TO 140                                                         4533.   
 130  P=SPLB(N)*EXP(-HPCON/STLB(N)*(H-SHLB(N)))                         4534.   
 140  T=STLB(N)+SDLB(N)*(H-SHLB(N))                                     4535.   
      D=P/T*28.9644E 05/8.31432E 03                                     4536.   
      RETURN                                                            4537.   
C                                                                       4538.   
 150  CONTINUE                                                          4539.   
      DO 160 N=2,8                                                      4540.   
 160  IF(P.GT.SPLB(N)) GO TO 170                                        4541.   
      N=9                                                               4542.   
 170  N=N-1                                                             4543.   
      IF(ABS(SDLB(N)).LT.1.E-04) GO TO 180                              4544.   
      H=SHLB(N)+STLB(N)/SDLB(N)*((SPLB(N)/P)**(SDLB(N)/HPCON)-1.)       4545.   
      GO TO 190                                                         4546.   
C ALOG
 180  H=SHLB(N)+STLB(N)/HPCON*LOG(SPLB(N)/P)                            4547.   
C ALOG
 190  T=STLB(N)+SDLB(N)*(H-SHLB(N))                                     4548.   
      D=P/T*28.9644E 05/8.31432E 03                                     4549.   
      RETURN                                                            4550.   
C                                                                       4551.   
 200  CONTINUE                                                          4552.   
      IF(NPHD.EQ.1) GO TO 240                                           4553.   
      IF(NPHD.EQ.2) GO TO 220                                           4554.   
      XX=D                                                              4555.   
      XI=DENS(1,NATM)                                                   4556.   
      IF(D.GT.XI) XX=XI                                                 4557.   
      IF(D.LT.5.0E-04) GO TO 280                                        4558.   
      DO 210 J=2,33                                                     4559.   
      XJ=DENS(J,NATM)                                                   4560.   
      IF(XX.GT.XJ) GO TO 260                                            4561.   
 210  XI=XJ                                                             4562.   
 220  XX=H                                                              4563.   
      XI=HTKM(1)                                                        4564.   
      IF(H.LT.XI) XX=XI                                                 4565.   
      IF(H.GT.99.9) GO TO 280                                           4566.   
      DO 230 J=2,33                                                     4567.   
      XJ=HTKM(J)                                                        4568.   
      IF(XX.LT.XJ) GO TO 260                                            4569.   
 230  XI=XJ                                                             4570.   
 240  XX=P                                                              4571.   
      XI=PRES(1,NATM)                                                   4572.   
      IF(P.GT.XI) XX=XI                                                 4573.   
      IF(P.LT.3.0E-04) GO TO 280                                        4574.   
      DO 250 J=2,33                                                     4575.   
      XJ=PRES(J,NATM)                                                   4576.   
      IF(XX.GT.XJ) GO TO 260                                            4577.   
 250  XI=XJ                                                             4578.   
 260  DELTA=(XX-XI)/(XJ-XI)                                             4579.   
      I=J-1                                                             4580.   
C ALOG
      IF(NPHD.NE.2) H=HTKM(I)+(HTKM(J)-HTKM(I))*LOG(XX/XI)/LOG(XJ/XI)   4581.   
C ALOG
      PI=PRES(I,NATM)                                                   4582.   
      PJ=PRES(J,NATM)                                                   4583.   
      DI=DENS(I,NATM)                                                   4584.   
      DJ=DENS(J,NATM)                                                   4585.   
      IF(NPHD.NE.1) P=PI+DELTA*(PJ-PI)                                  4586.   
      IF(NPHD.NE.3) D=DI+DELTA*(DJ-DI)                                  4587.   
      T=TEMP(I,NATM)+DELTA*(TEMP(J,NATM)-TEMP(I,NATM))                  4588.   
      O=OZON(I,NATM)/DI+DELTA*(OZON(J,NATM)/DJ-OZON(I,NATM)/DI)         4589.   
      Q=WVAP(I,NATM)/DI+DELTA*(WVAP(J,NATM)/DJ-WVAP(I,NATM)/DI)         4590.   
      ES=10.**(9.4051-2353./T)                                          4591.   
      IF(P.LT.PI) PI=P                                                  4592.   
      S=1.E+06                                                          4593.   
      RS=(PI-ES+0.622*ES)/(0.622*ES)                                    4594.   
      IF(RS.GT.1.E-06) S=1./RS                                          4595.   
      OI=O                                                              4596.   
      QI=Q                                                              4597.   
      OCM=0.                                                            4598.   
      WCM=0.                                                            4599.   
      DO 270 K=J,33                                                     4600.   
      PJ=PRES(K,NATM)                                                   4601.   
      DJ=DENS(K,NATM)                                                   4602.   
      OJ=OZON(K,NATM)/DJ                                                4603.   
      QJ=WVAP(K,NATM)/DJ                                                4604.   
      DP=PI-PJ                                                          4605.   
      OCM=OCM+0.5*(OI+OJ)*DP                                            4606.   
      WCM=WCM+0.5*(QI+QJ)*DP                                            4607.   
      OI=OJ                                                             4608.   
      QI=QJ                                                             4609.   
 270  PI=PJ                                                             4610.   
      WCM=WCM/0.980*22420.7/18.0                                        4611.   
      OCM=OCM/0.980*22420.7/48.0                                        4612.   
      RETURN                                                            4613.   
 280  T=210.0                                                           4614.   
      IF(NATM.EQ.6) T=186.87                                            4615.   
      O=1.E-10                                                          4616.   
      Q=1.E-10                                                          4617.   
      S=1.E-10                                                          4618.   
      OCM=1.E-10                                                        4619.   
      WCM=1.E-10                                                        4620.   
      IF(NPHD.NE.1) P=1.E-05                                            4621.   
      IF(NPHD.NE.2) H=99.99                                             4622.   
      IF(NPHD.NE.3) D=2.E-05                                            4623.   
      RETURN                                                            4624.   
      END                                                               4625.   
      FUNCTION PFOFTK(WAVNA,WAVNB,TK)                                   4626.   
C     ------------------------------------------------------------------4627.   
C                                                                       4628.   
C        INPUT DATA                                                     4629.   
C                  WAVNA,WAVNB  SPECLTRAL INTERVAL IN WAVENUMBERS       4630.   
C                               (ORDER OF WAVNA,WAVNB NOT IMPORTANT)    4631.   
C                                                                       4632.   
C                  TK           ABSOLUTE TEMPERATURE IN DEGREES KELVIN  4633.   
C                                                                       4634.   
C       OUTPUT DATA                                                     4635.   
C                  PFOFTK       PLANCK FLUX (W/M**2)                    4636.   
C                                                                       4637.   
C                                                                       4638.   
C           REMARKS                                                     4639.   
C                   PLANCK INTENSITY (W/M**2/STER) IS GIVEN BY PFOFTK/PI4640.   
C                                                                       4641.   
C     ------------------------------------------------------------------4642.   
c     IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               4643.   
      DIMENSION BN(21),BD(21)                                           4644.   
      DATA BN/1.D0,-1.D0,1.D0,-1.D0,1.D0,-1.D0,5.D0,-691.D0,7.D0        4645.   
     1,-3617.D0,43867.D0,-174611.D0,854513.D0,-236364091.D0             4646.   
     2,8553103.D0,-23749461029.D0,8615841276005.D0,-7709321041217.D0    4647.   
     3,2577687858367.D0,-2631527155305348D 04,2929993913841559D0/       4648.   
      DATA BD/1.D0,2.D0,6.D0,30.D0,42.D0,30.D0,66.D0,2730.D0,6.D0       4649.   
     1,510.D0,798.D0,330.D0,138.D0,2730.D0,6.D0,870.D0,14322.D0         4650.   
     2,510.D0,6.D0,1919190.D0,6.D0/                                     4651.   
      DATA PI4/97.40909103400244D0/                                     4652.   
C     DATA  PI/3.141592653589793D0/                                     4653.   
      DATA HCK/1.43879D0/                                               4654.   
      DATA DGXLIM/1.D-06/                                               4655.   
      PFOFTK=0.D0                                                       4656.   
      IF(TK.LT.1.D-06) RETURN                                           4657.   
      DO 160 II=1,2                                                     4658.   
      IF(II.EQ.1) X=HCK*WAVNA/TK                                        4659.   
      IF(II.EQ.2) X=HCK*WAVNB/TK                                        4660.   
      IF(X.GT.2.3D0) GO TO 120                                          4661.   
      XX=X*X                                                            4662.   
      GSUM=1.D0/3.D0-X/8.D0+XX/60.D0                                    4663.   
      NB=3                                                              4664.   
      XNF=XX/2.D0                                                       4665.   
      DO 100 N=4,38,2                                                   4666.   
      NB=NB+1                                                           4667.   
      NNB=NB                                                            4668.   
      B=BN(NB)/BD(NB)                                                   4669.   
      XN3=N+3                                                           4670.   
      XNM=N*(N-1)                                                       4671.   
      XNF=XNF*(XX/XNM)                                                  4672.   
      DG=B/XN3*XNF                                                      4673.   
      GSUM=GSUM+DG                                                      4674.   
      DGB=DG                                                            4675.   
      IF(DABS(DG).LT.DGXLIM) GO TO 110                                  4676.   
 100  CONTINUE                                                          4677.   
 110  GX=GSUM*XX*X                                                      4678.   
      GO TO 150                                                         4679.   
 120  GSUM=PI4/15.D0                                                    4680.   
      DO 130 N=1,20                                                     4681.   
      NNB=N                                                             4682.   
      XN=N                                                              4683.   
      XNN=XN*XN                                                         4684.   
      XNX=XN*X                                                          4685.   
      IF(XNX.GT.100.D0) GO TO 140                                       4686.   
      GTERM=(X*X*(3.D0+XNX)+6.D0*(1.D0+XNX)/XNN)/XNN                    4687.   
      DG=GTERM*DEXP(-XNX)                                               4688.   
      GSUM=GSUM-DG                                                      4689.   
      DGB=DG                                                            4690.   
      IF(DG.LT.DGXLIM) GO TO 140                                        4691.   
 130  CONTINUE                                                          4692.   
 140  GX=GSUM                                                           4693.   
 150  CONTINUE                                                          4694.   
      IF(II.EQ.1) GXA=GX                                                4695.   
      IF(II.EQ.2) GXB=GX                                                4696.   
 160  CONTINUE                                                          4697.   
      PNORM=15.D0/PI4                                                   4698.   
      PFOFTK=DABS(GXB-GXA)*PNORM                                        4699.   
      PFOFTK=PFOFTK*5.6692D-08*TK**4                                    4700.   
      RETURN                                                            4701.   
      END                                                               4702.   
      FUNCTION TKOFPF(WAVNA,WAVNB,FLUXAB)                               4703.   
C     ------------------------------------------------------------------4704.   
C                                                                       4705.   
C        INPUT DATA                                                     4706.   
C------------------                                                     4707.   
C                  WAVNA,WAVNB  SPECLTRAL INTERVAL IN WAVENUMBERS       4708.   
C                               (ORDER OF WAVNA,WAVNB NOT IMPORTANT)    4709.   
C                  FLUXAB       PLANCK FLUX (W/M**2) IN INTERVAL        4710.   
C                                                       (WAVNA,WAVNB)   4711.   
C                                                                       4712.   
C       OUTPUT DATA                                                     4713.   
C------------------                                                     4714.   
C                  TK           BRIGHTNESS TEMPERATURE IN DEGREES KELVIN4715.   
C                                                                       4716.   
C                                                                       4717.   
C           REMARKS                                                     4718.   
C------------------                                                     4719.   
C                   TKOFPF IS INVERSE FUNCTION OF PFOFTK(WAVNA,WAVNB,TK)4720.   
C                   THE OUTPUT OF TKOFPF SATISFIES THE IDENTITY         4721.   
C                                 FLUXAB=PFOFTK(WAVNA,WAVNB,TK)         4722.   
C                   (UNITS FOR FLUXAB AND PFOFTK MUST BE IDENTICAL)     4723.   
C                                                                       4724.   
C     ------------------------------------------------------------------4725.   
c     IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               4726.   
      LOGICAL LOGFIT                                                    4727.   
      DATA DELFIT/1.D-06/                                               4728.   
      DATA NMAX/20/                                                     4729.   
C     IF(FLUXAB.LE.0.D0) RETURN                                         4730.   
      LOGFIT=.FALSE.                                                    4731.   
      NFIT=0                                                            4732.   
      PF=FLUXAB                                                         4733.   
      XA=0.D0                                                           4734.   
      YA=0.D0                                                           4735.   
      XB=250.D0                                                         4736.   
      YB=PFOFTK(WAVNA,WAVNB,XB)                                         4737.   
      XX=PF*XB/YB                                                       4738.   
      YY=PFOFTK(WAVNA,WAVNB,XX)                                         4739.   
      IF(DABS(YY-PF).LT.DELFIT) GO TO 200                               4740.   
      IF((YY/PF).LT.0.5D0) GO TO 150                                    4741.   
      IF((YY/PF).GT.2.0D0) GO TO 170                                    4742.   
      IF(XX.GT.XB) GO TO 110                                            4743.   
      XC=XB                                                             4744.   
      YC=YB                                                             4745.   
      XB=XX                                                             4746.   
      YB=YY                                                             4747.   
      GO TO 120                                                         4748.   
  110 XC=XX                                                             4749.   
      YC=YY                                                             4750.   
  120 XBA=XB-XA                                                         4751.   
      XCA=XC-XA                                                         4752.   
      XBC=XB-XC                                                         4753.   
      YBA=YB-YA                                                         4754.   
      YCA=YC-YA                                                         4755.   
      YBC=YB-YC                                                         4756.   
      NFIT=NFIT+1                                                       4757.   
      IF(NFIT.GT.NMAX) GO TO 200                                        4758.   
      YXBA=YBA/XBA                                                      4759.   
      YXCA=YCA/XCA                                                      4760.   
      C=(YXBA-YXCA)/XBC                                                 4761.   
      B=YXBA-(XB+XA)*C                                                  4762.   
      A=YA-XA*(B+XA*C)                                                  4763.   
      ROOT=DSQRT(B*B+4.D0*C*(PF-A))                                     4764.   
      XX=0.5D0*(ROOT-B)/C                                               4765.   
      IF(XX.LT.XA.OR.XX.GT.XC) XX=-0.5D0*(ROOT+B)/C                     4766.   
      YY=PFOFTK(WAVNA,WAVNB,XX)                                         4767.   
      IF(LOGFIT) YY=DLOG(YY)                                            4768.   
      IF(DABS(YY-PF).LT.DELFIT) GO TO 200                               4769.   
      IF(XX.GT.XB) GO TO 130                                            4770.   
      XC=XB                                                             4771.   
      YC=YB                                                             4772.   
      GO TO 140                                                         4773.   
  130 XA=XB                                                             4774.   
      YA=YB                                                             4775.   
  140 XB=XX                                                             4776.   
      YB=YY                                                             4777.   
      GO TO 120                                                         4778.   
  150 XA=XX                                                             4779.   
      YA=YY                                                             4780.   
  160 XC=XB                                                             4781.   
      YC=YB                                                             4782.   
      XB=XB/2.D0                                                        4783.   
      YB=PFOFTK(WAVNA,WAVNB,XB)                                         4784.   
      IF(YB.LT.YA) GO TO 190                                            4785.   
      IF(YB.GT.PF) GO TO 160                                            4786.   
      XA=XB                                                             4787.   
      YA=YB                                                             4788.   
      GO TO 190                                                         4789.   
  170 XC=XX                                                             4790.   
      YC=YY                                                             4791.   
  180 XA=XB                                                             4792.   
      YA=YB                                                             4793.   
      XB=XB*2.D0                                                        4794.   
      YB=PFOFTK(WAVNA,WAVNB,XB)                                         4795.   
      IF(YB.GT.YC) GO TO 190                                            4796.   
      IF(YB.LT.PF) GO TO 180                                            4797.   
      XC=XB                                                             4798.   
      YC=YB                                                             4799.   
  190 XB=XA+(PF-YA)*(XC-XA)/(YC-YA)                                     4800.   
      YB=PFOFTK(WAVNA,WAVNB,XB)                                         4801.   
      XX=XB                                                             4802.   
      IF(DABS(YB-PF).LT.DELFIT) GO TO 200                               4803.   
      PF=DLOG(PF)                                                       4804.   
      YA=DLOG(YA)                                                       4805.   
      YB=DLOG(YB)                                                       4806.   
      YC=DLOG(YC)                                                       4807.   
      LOGFIT=.TRUE.                                                     4808.   
      GO TO 120                                                         4809.   
  200 TKOFPF=XX                                                         4810.   
      RETURN                                                            4811.   
      END                                                               4812.   
      SUBROUTINE WRITER(INDEX,KPAGE)                                    4813.   

#include "B83XX.COM"
 
      DIMENSION SRAOC(15),SRAEA(15),SRAOI(15),SRALI(15),SRASN(15)       4875.   
C                                                                       4876.   
      DIMENSION SRBALB(6),SRXALB(6)                                     4877.   
      EQUIVALENCE (SRBXAL(1,1),SRBALB(1)),(SRBXAL(1,2),SRXALB(1))       4878.   
C                                                                       4879.   
     +,               (BXA(20),BSNVIS),(BXA(21),BSNNIR)                 4880.   
     +,               (BXA(22),XSNVIS),(BXA(23),XSNNIR)                 4881.   
C                                                                       4882.   
      EQUIVALENCE (ASNALB(1),ASNVIS),(ASNALB(2),ASNNIR)                 4883.   
      EQUIVALENCE (AOIALB(1),AOIVIS),(AOIALB(2),AOINIR)                 4884.   
      EQUIVALENCE (ALIALB(1),ALIVIS),(ALIALB(2),ALINIR)                 4885.   
C                                                                       4886.   
      EQUIVALENCE                                                       4887.   
     +          (FEMTRA(1),ECLTRA),  (FZASRA(1),ZCLSRA)                 4888.   
     +,         (FEMTRA(2),EOCTRA),  (FZASRA(2),ZOCSRA)                 4889.   
     +,         (FEMTRA(3),ESNTRA),  (FZASRA(3),ZSNSRA)                 4890.   
     +,         (FEMTRA(4),EICTRA),  (FZASRA(4),ZICSRA)                 4891.   
     +,         (FEMTRA(5),EDSTRA),  (FZASRA(5),ZDSSRA)                 4892.   
     +,         (FEMTRA(6),EVGTRA),  (FZASRA(6),ZVGSRA)                 4893.   
C                                                                       4894.   
      EQUIVALENCE (IMG(1),IMGAS1),(IMG(2),IMGAS2)                       4895.   
      EQUIVALENCE (ILG(1),ILGAS1),(ILG(2),ILGAS2)                       4896.   
C                                                                       4897.   
      EQUIVALENCE (ID5(1),IDPROG),(ID5(2),ID2TRD),(ID5(3),ID3SRD)       4898.   
      EQUIVALENCE (ID5(4),ID4VEG),(ID5(5),ID5FOR)                       4899.   
C                                                                       4900.   
      EQUIVALENCE (PVT( 1),DESRT),(PVT( 2),TNDRA),(PVT( 3),GRASS)       4901.   
     +           ,(PVT( 4),SHRUB),(PVT( 5),TREES),(PVT( 6),DECID)       4902.   
     +           ,(PVT( 7),EVERG),(PVT( 8),RAINF),(PVT( 9),ROCKS)       4903.   
     +           ,(PVT(10),CROPS),(PVT(11),ALGAE)                       4904.   
C                                                                       4905.   
      EQUIVALENCE (FRC(1),FRACCC),(FRC(2),  FCHI),(FRC(3),  FCMI)       4906.   
     +                           ,(FRC(4),  FCLO),(FRC(5),  FCOV)       4907.   
C                                                                       4908.   
C                                                                       4909.   
      CHARACTER*8 FTYPE                                                 4910.   
      DIMENSION BGFLUX(25),BGFRAC(25),TAUSUM(25)                        4911.   
      DIMENSION SUM0(15),SUM1(40),SUM2(40),SUM3(40),FTYPE(5),AUXGAS(4)  4912.   
      DATA FTYPE/'DOWNWARD','  UPWARD','UPWD NET','COOLRATE','FRACTION'/4913.   
      DATA AUXGAS/1H0,1HL,1HX,1HX/                                      4914.   
      DATA P0/1013.25/                                                  4915.   
C                                                                       4916.   
      INDJ=MOD(INDEX,10)                                                4917.   
      IF(INDJ.LT.1) INDJ=10                                             4918.   
      INDI=1                                                            4919.   
      IF(INDEX.LT.11) INDI=INDJ                                         4920.   
      DO 9999 INDX=INDI,INDJ                                            4921.   
C                                                                       4922.   
      IF(INDEX.EQ.0) GO TO 10                                           4923.   
      GO TO (100,200,300,400,500,600,700,800,900,1000),INDX             4924.   
C                                                                       4925.   
C-------------                                                          4926.   
   10 CONTINUE                                                          4927.   
C-------------                                                          4928.   
C                                                                       4929.   
      NPAGE=1                                                           4930.   
      WRITE(6,6001) NPAGE                                               4931.   
 6001 FORMAT(1I1,'(1) RADCOM M/R: (CONTROL/INPUT PARAMETERS)'           4932.   
     +      ,' DEFAULT VALUES & MODIFICATIONS'/)                        4933.   
      WRITE(6,6002)                                                     4934.   
 6002 FORMAT(20X,'PARAMETER/VALUE',5X,'COMMENTS RE PARAMETER DEFAULT'   4935.   
     +      ,' VALUE AND PARAMETER RANGE AND EFFECT'/10X,'AEROSOLS')    4936.   
      WRITE(6,6003)                                                     4937.   
 6003 FORMAT(20X,'FGOLDH(1) = 1.0',5X,'STRATOSPHERIC AEROSOL, GLOBAL'   4938.   
     +      ,' BACKGROUND - TAU(.55) = 0.005'                           4939.   
     +      /20X,'FGOLDH(2) = 1.0',5X,' TROPOSPHERIC AEROSOL OVER LAND' 4940.   
     +      ,' BACKGROUND:  TAU(.55) = 0.125'                           4941.   
     +      /20X,'FGOLDH(3) = 0.0',5X,' TROPOSPHERIC AEROSOL OVER LAND' 4942.   
     +      ,' BACKGROUND:  TAU(.55) = 0.125 (FOR FGOLDH(3)=1.0'        4943.   
     +      /)                                                          4944.   
      GO TO 9999                                                        4945.   
C                                                                       4946.   
C-------------                                                          4947.   
  100 CONTINUE                                                          4948.   
C-------------                                                          4949.   
C                                                                       4950.   
C                                                                       4951.   
      NPAGE=1                                                           4952.   
      IF(INDEX.LT.11) NPAGE=KPAGE                                       4953.   
      WRITE(6,6101) NPAGE,LASTVC,KFORCE                                 4954.   
      WRITE(6,6102)                                                     4955.   
      IDPROG=ID5(1)                                                     4956.   
      ID2TRD=ID5(2)                                                     4957.   
      ID3SRD=ID5(3)                                                     4958.   
      ID4VEG=ID5(4)                                                     4959.   
      ID5FOR=ID5(5)                                                     4960.   
      FACTOR=P0/(PLB(1)-PLB(2))*1.25                                    4961.   
      PPMCO2=ULGAS(1,2)*FACTOR                                          4962.   
      PPMO2 =ULGAS(1,4)*FACTOR                                          4963.   
      PPMN2O=ULGAS(1,6)*FACTOR                                          4964.   
      PPMCH4=ULGAS(1,7)*FACTOR                                          4965.   
      PPMF11=ULGAS(1,8)*FACTOR                                          4966.   
      PPMF12=ULGAS(1,9)*FACTOR                                          4967.   
      WRITE(6,6103) (FULGAS(I),I=1,9),(FGOLDH(I),I=1,5)                 4968.   
      IF(KGASSR.GT.0.OR.KAERSR.GT.0)                                    4969.   
     +WRITE(6,6104) (FULGAS(I+9),I=1,9),(FGOLDH(I+9),I=1,5)             4970.   
!
! === Chien Wang 121797
!
#if ( defined CPL_CHEM  && !defined FIXED_FOR )
	WRITE(6,6105) PPMCO2,PPMO3,PPMO2,PPMN2O,PPMCH4,PPMF11,PPMF12
#else
      WRITE(6,6105) PPMCO2,PPMO2,PPMN2O,PPMCH4,PPMF11,PPMF12            4971.   
#endif
     +             ,(FGOLDH(I),I=6,9),NV                                4972.   
      WRITE(6,6106) TAUMIN,TLGRAD,EOCTRA,ZOCSRA,FMARCL,FCLDTR,NTRACE    4973.   
     +             ,IDPROG,IMGAS1,KEEPRH,KGASSR,LAYRAD                  4974.   
      WRITE(6,6107) FRACSL,TKCICE,ESNTRA,ZSNSRA,WETTRA,FCLDSR,ITR(1)    4975.   
     +             ,ID2TRD,IMGAS2,KEEPAL,KAERSR,NL                      4976.   
      WRITE(6,6108) RATQSL,FLONO3,EICTRA,ZICSRA,WETSRA,FALGAE,ITR(2)    4977.   
     +             ,ID3SRD,ILGAS1,ISOSCT,KFRACC,NLP                     4978.   
      WRITE(6,6109) FOGTSL,ECLTRA,EDSTRA,ZDSSRA,DMOICE,FRAYLE,ITR(3)    4979.   
     +             ,ID4VEG,ILGAS2,IHGSCT,MARCLD,JMLAT                   4980.   
      WRITE(6,6110) PTLISO,ZCLSRA,EVGTRA,ZVGSRA,DMLICE,LICETK,ITR(4)    4981.   
     +             ,ID5FOR,KWVCON,LAPGAS,NORMS0,IMLON                   4982.   
C                                                                       4983.   
 6101 FORMAT(1I1,'(1) RADCOM 1/F: (CONTROL/INPUT PARAMETERS)'           4984.   
     +      ,' (GAS/AEROSOL REFERENCE AMOUNT SCALE FACTORS,'            4985.   
     +      ,' DEFAULTS & OPTIONS IN FORCE)  LASTVC=',I7                4986.   
     +      /1X,113('-'),'  KFORCE=',I10)                               4987.   
 6102 FORMAT(4X,'GAS: ','H2O',5X,'CO2',7X,'O3',6X,'O2',5X,'NO2'         4988.   
     +      ,5X,'N2O',5X,'CH4',6X,'CCL3F1',3X,'CCL2F2'                  4989.   
     +      ,3X,'AERSOL: GLOBAL    OCEAN     LAND  DESERT    HAZE')     4990.   
 6103 FORMAT(1X,'FULGAS=',1P,1E7.1,1P,2E8.1,1P,1E9.1,1P,3E8.1,1P,2E9.1  4991.   
     +      ,3X,'FGOLDH=',1P,1E7.1,1P,2E9.2,1P,2E8.1)                   4992.   
 6104 FORMAT(1H+,T84,'T'                                                4993.   
     +      /1X,'FULGAS=',1P,1E7.1,1P,2E8.1,1P,1E9.1,1P,3E8.1,1P,2E9.1  4994.   
     +      ,' S',1X,'FGOLDH=',1P,1E7.1,1P,2E9.2,1P,2E8.1)              4995.   
!
! === Chien Wang 121797
!
#if ( defined CPL_CHEM  && !defined FIXED_FOR )
 6105 FORMAT(1X,'PPM(1)=$',6X,0P,1F8.3,f8.4,F8.0,8X,F8.4,F8.4,1X,F8.7
#else
 6105 FORMAT(1X,'PPM(1)=$',6X,0P,1F8.3,9X,F8.0,8X,F8.4,F8.4,1X,F8.7     4996.   
#endif
     +      ,1X,F8.7,3X,'TRACER=',1P,E7.1,1P,2E9.2,1P,E8.1,'   NV=',I2) 4997.   
 6106 FORMAT(1X,'TAUMIN=',1P,E7.1,1X,'TLGRAD=',0P,F4.1,' EOCTRA=',F3.1  4998.   
     +      ,1X,'ZOCSRA=',  F3.1,1X,'FMARCL=',  F4.2,1X,'FCLDTR=',F3.1  4999.   
     +      ,1X,'NTRACE=',    I2,3X,'IDPROG=',    I4,1X,'IMGAS1=',  I1  5000.   
     +      ,1X,'KEEPRH=',    I1,1X,'KGASSR=',    I1,1X,'LAYRAD=',  I2) 5001.   
 6107 FORMAT(1X,'FRACSL=',1P,E7.1,1X,'TKCICE=',0P,F4.0,' ESNTRA=',F3.1  5002.   
     +      ,1X,'ZSNSRA=',  F3.1,1X,'WETTRA=',  F4.2,1X,'FCLDSR=',F3.1  5003.   
     +      ,1X,'ITR(1)=',    I2,3X,'ID2TRD=',    I4,1X,'IMGAS2=',  I1  5004.   
     +      ,1X,'KEEPAL=',    I1,1X,'KAERSR=',    I1,1X,'    NL=',  I2) 5005.   
 6108 FORMAT(1X,'RATQSL=',  F4.2,4X,'FLONO3=',  F4.1,1X,'EICTRA=',F3.1  5006.   
     +      ,1X,'ZICSRA=',  F3.1,1X,'WETSRA=',  F4.2,1X,'FALGAE=',F3.1  5007.   
     +      ,1X,'ITR(2)=',    I2,3X,'ID3SRD=',    I4,1X,'ILGAS1=',  I1  5008.   
     +      ,1X,'ISOSCT=',    I1,1X,'KFRACC=',    I1,1X,'   NLP=',  I2) 5009.   
 6109 FORMAT(1X,'FOGTSL=',  F4.2,4X,'ECLTRA=',  F4.2,1X,'EDSTRA=',F3.1  5010.   
     +      ,1X,'ZDSSRA=',  F3.1,1X,'DMOICE=',  F4.1,1X,'FRAYLE=',F3.1  5011.   
     +      ,1X,'ITR(3)=',    I2,3X,'ID4VEG=',    I4,1X,'ILGAS2=',  I1  5012.   
     +      ,1X,'IHGSCT=',    I1,1X,'MARCLD=',    I1,1X,' JMLAT=',  I2) 5013.   
 6110 FORMAT(1X,'PTLISO=',1PE7.1,1X,'ZCLSRA=',0PF4.2,1X,'EVGTRA=',F3.1  5014.   
     +      ,1X,'ZVGSRA=',  F3.1,1X,'DMLICE=',  F4.1,1X,'LICETK=',  I3  5015.   
     +      ,1X,'ITR(4)=',    I2,3X,'ID5FOR=',    I4,1X,'KWVCON=',  I1  5016.   
     +      ,1X,'LAPGAS=',    I1,1X,'NORMS0=',    I1,1X,'IMLON=',  I3)  5017.   
      GO TO 9999                                                        5018.   
C                                                                       5019.   
C-------------                                                          5020.   
  200 CONTINUE                                                          5021.   
C-------------                                                          5022.   
C                                                                       5023.   
      NPAGE=0                                                           5024.   
      IF(INDEX.LT.11) NPAGE=KPAGE                                       5025.   
      WRITE(6,6201) NPAGE,AUXGAS(LUXGAS+1),S0,COSZ                      5026.   
      DO 202 K=1,9                                                      5027.   
      DO 201 L=1,NL                                                     5028.   
      IF(LUXGAS.EQ.0) UXGAS(L,K)=U0GAS(L,K)                             5029.   
  201 IF(LUXGAS.EQ.1) UXGAS(L,K)=ULGAS(L,K)                             5030.   
  202 CONTINUE                                                          5031.   
      IF(LUXGAS.LT.2) GO TO 205                                         5032.   
      LGS=(LUXGAS-2)*9                                                  5033.   
      DO 203 L=1,NL                                                     5034.   
      UXGAS(L,1)=U0GAS(L,1)*FULGAS(1+LGS)                               5035.   
      UXGAS(L,3)=U0GAS(L,3)*FULGAS(3+LGS)                               5036.   
  203 UXGAS(L,5)=U0GAS(L,5)*FULGAS(5+LGS)                               5037.   
C                                                                       5038.   
      DO 204 L=1,NL                                                     5039.   
      UXGAS(L,2)=U0GAS(L,2)*FULGAS(2+LGS)                               5040.   
      UXGAS(L,4)=U0GAS(L,4)*FULGAS(4+LGS)                               5041.   
      UXGAS(L,6)=U0GAS(L,6)*FULGAS(6+LGS)                               5042.   
      UXGAS(L,7)=U0GAS(L,7)*FULGAS(7+LGS)                               5043.   
      UXGAS(L,8)=U0GAS(L,8)*FULGAS(8+LGS)                               5044.   
  204 UXGAS(L,9)=U0GAS(L,9)*FULGAS(9+LGS)                               5045.   
  205 CONTINUE                                                          5046.   
      DO 206 N=1,NL                                                     5047.   
      L=NLP-N                                                           5048.   
      WRITE(6,6202) L,PLB(L),HLB(L),TLB(L),TLT(L),TLM(L)                5049.   
     +     ,(UXGAS(L,K),K=1,9),CLDTAU(L),SHL(L),RHL(L)                  5050.   
  206 CONTINUE                                                          5051.   
      DO 207 I=1,15                                                     5052.   
  207 SUM0(I)=0.                                                        5053.   
      DO 210 L=1,NL                                                     5054.   
      DO 208 I=1,9                                                      5055.   
  208 SUM0(I)=SUM0(I)+ULGAS(L,I)                                        5056.   
      DO 209 I=1,4                                                      5057.   
  209 SUM0(11+I)=SUM0(11+I)+TRACER(L,I)                                 5058.   
  210 SUM0(10)=SUM0(10)+CLDTAU(L)                                       5059.   
      DO 212 J=1,NGOLDH                                                 5060.   
      TAU55=0.                                                          5061.   
      DO 211 I=1,NAERO                                                  5062.   
  211 TAU55=TAU55+AGOLDH(I,J)*FGOLDH(J)                                 5063.   
  212 SUM0(11)=SUM0(11)+TAU55                                           5064.   
      TGMEAN=POCEAN*TGO**4+PEARTH*TGE**4+PLICE*TGLI**4+POICE*TGOI**4    5065.   
      TGMEAN=SQRT(TGMEAN)                                               5066.   
      TGMEAN=SQRT(TGMEAN)                                               5067.   
      WRITE(6,6203) SUM0(11),(SUM0(I),I=1,10)                           5068.   
C     WRITE(6,6204) POCEAN,   TGO, AGESN, ZOICE,LASTVC, DESRT, DECID    5069.   
C    +             ,SRAOC(1),SRAEA(1),SRAOI(1),SRALI(1),SRASN(1)        5070.   
C    +             ,SRDALB(1),SRXALB(1)                                 5071.   
C     WRITE(6,6205) PEARTH,   TGE, SNOWE,WEARTH, PSIG0, TNDRA, EVERG    5072.   
C     WRITE(6,6206)  POICE,  TGOI,SNOWOI,FRACCC, ALGAE, GRASS, RAINF    5073.   
C     WRITE(6,6207)  PLICE,  TGLI,SNOWLI, JYEAR,TRACR1, SHRUB, ROCKS    5074.   
C     WRITE(6,6208) MEANAL,TGMEAN,EXSNEA,  JDAY,TRACR2, TREES, CROPS    5075.   
C     WRITE(6,6209) KALVIS,   TSL,EXSNOI,  JLAT,TRACR3,  FCHI,  FCLO    5076.   
C     WRITE(6,6210) LUXGAS,  WMAG,EXSNLI,  ILON,TRACR4,  FCMI,  FCOV    5077.   
C                                                                       5078.   
      WRITE(6,6204) POCEAN,TGO,AGESN,WMAG,SUM0(12),JYEAR,BSNVIS,BSNNIR  5079.   
     +             ,LASTVC                                              5080.   
      WRITE(6,6205) PEARTH,TGE,SNOWE,WEARTH,SUM0(13),JDAY,XSNVIS,XSNNIR 5081.   
      WRITE(6,6206) POICE,TGOI,SNOWOI,ZOICE,SUM0(14),JLAT               5082.   
     +             ,(SRBALB(I),I=1,6)                                   5083.   
      WRITE(6,6207) PLICE,TGLI,SNOWLI,FRC(5),SUM0(15),ILON              5084.   
     +             ,(SRXALB(I),I=1,6)                                   5085.   
      WRITE(6,6208) TGMEAN,LUXGAS,PSUM,TSL,MEANAL,KALVIS,(PVT(I),I=1,11)5086.   
      WRITE(6,6209) (BXA(I),I=1,19)                                     5087.   
 6201 FORMAT(1I1,'(2) RADCOM G/L: (INPUT DATA)'                         5088.   
     +      ,T41,'    ABSORBER AMOUNT PER LAYER:'                       5089.   
     +      ,'  U',1A1,'GAS(L,K) IN CM**3(STP)/CM**2'                   5090.   
     +      ,T109,'S0=',F8.3,3X,'COSZ=',F6.4/1X,132('-')                5091.   
     +      /' LN     PLB   HLB    TLB    TLT    TLM      '             5092.   
     +      ,'H2O     CO2    O3       O2    NO2      N2O    CH4'        5093.   
     +      ,'   CCL3F1   CCL2F2  CLDTAU   SHL    RHL ')                5094.   
 6202 FORMAT(1X,I2,F9.3,F6.2,3F7.2,F9.3,F8.3,1X,F6.5,F8.0,1P,1E9.2      5095.   
     +      ,1X,0P,1F6.5,F7.4,1P,2E9.2,0P,1F7.2,1X,F7.6,1X,F5.4)        5096.   
 6203 FORMAT( 1X,'$SUM AERSOL=',F5.3,7X,'$COLUMN AMOUNT',F9.3           5097.   
     +  ,F8.3,1X,F6.5,F8.0,1P,1E9.2,1X,0P,1F6.5,F7.4,1P,2E9.2,0P,F7.2)  5098.   
 6204 FORMAT(/1X,'POCEAN=',F6.4,'    TGO=' ,F6.2,1X,' AGESN=',F6.3      5099.   
     +      , 1X,'  WMAG=',F6.3,' TRACER 1=',F5.3,' JYEAR=',I4          5100.   
     +      , 3X,'BSNVIS=',F6.4,' BSNNIR=' ,F6.4,7X,'LASTVC=',I7)       5101.   
 6205 FORMAT(    ' PEARTH=',F6.4,'    TGE=',F6.2,'  SNOWE=',F6.3        5102.   
     +      ,    ' WEARTH=',F6.3,' $SUMS: 2=',F5.3                      5103.   
     +      ,     '  JDAY=',I4  ,2X,' XSNVIS=',F6.4,' XSNNIR=',F6.4     5104.   
     +      , 8X,'NIRALB VISALB')                                       5105.   
 6206 FORMAT(    '  POICE=',F6.4,'   TGOI=',F6.2,' SNOWOI=',F6.3        5106.   
     +      ,    '  ZOICE=',F6.3,'        3=',F5.3                      5107.   
     +      ,     '  JLAT=',I4,  2X,' SRBALB=',F6.4                     5108.   
     +      ,4F7.4,F7.4)                                                5109.   
 6207 FORMAT(    '  PLICE=',F6.4,'   TGLI=',F6.2,' SNOWLI=',F6.3        5110.   
     +      ,    ' FRC(5)=',F6.3,'        4=',F5.3                      5111.   
     +      ,     '  ILON=',I4,  2X,' SRXALB=',F6.4                     5112.   
     +      ,4F7.4,F7.4)                                                5113.   
 6208 FORMAT( 1X,13('-'),'$TGMEAN=',F6.2,14X,' LUXGAS=',I1,5X           5114.   
     +      ,1X,'DESERT TUNDRA GRASSL SHRUBS TREES  DECIDF'             5115.   
     +       ,' EVERGF',' RAINF','  ROCKS','  CROPS','  ALGAE'          5116.   
     +      /    '  $PSUM=',F6.4,'    TSL=',F6.2,' MEANAL=',I1          5117.   
     +       ,5X,' KALVIS=',I1,T54,'PVT=',F6.4,10F7.4)                  5118.   
 6209 FORMAT(' BOCVIS BOCNIR XOCVIS XOCNIR|BEAVIS BEANIR XEAVIS XEANIR' 5119.   
     +      ,'|BOIVIS BOINIR XOIVIS XOINIR|BLIVIS BLINIR XLIVIS XLINIR' 5120.   
     +      ,'|EXPSNE|EXPSNO|EXPSNL'/1X,F6.4,18F7.4)                    5121.   
      GO TO 9999                                                        5122.   
C                                                                       5123.   
C-------------                                                          5124.   
  300 CONTINUE                                                          5125.   
C-------------                                                          5126.   
C                                                                       5127.   
      NPAGE=0                                                           5128.   
      IF(INDEX.LT.11) NPAGE=KPAGE                                       5129.   
      IF(NL.GT.13) NPAGE=1                                              5130.   
      L=NLP                                                             5131.   
      STNFLB=SRNFLB(L)-TRNFLB(L)                                        5132.   
      WRITE(6,6301) NPAGE,NORMS0                                        5133.   
      WRITE(6,6302) L,PLB(L),HLB(L),TLB(L)                              5134.   
     +             ,TRDFLB(L),TRUFLB(L),TRNFLB(L)                       5135.   
     +             ,SRDFLB(L),SRUFLB(L),SRNFLB(L),STNFLB                5136.   
      DO 301 N=1,NL                                                     5137.   
      L=NLP-N                                                           5138.   
      CRHRF=8.4167/(PLB(L)-PLB(L+1))                                    5139.   
      STNFLB=SRNFLB(L)-TRNFLB(L)                                        5140.   
      STFHR =SRFHRL(L)-TRFCRL(L)                                        5141.   
      TRDCR =TRFCRL(L)*CRHRF                                            5142.   
      SRDHR =SRFHRL(L)*CRHRF                                            5143.   
      STDHR=STFHR*CRHRF                                                 5144.   
      SRALB =SRUFLB(L)/(SRDFLB(L)+1.E-10)                               5145.   
      SRXVIS=SRXATM(1)                                                  5146.   
      SRXNIR=SRXATM(2)                                                  5147.   
      WRITE(6,6303) L,PLB(L),HLB(L),TLB(L),TLT(L)                       5148.   
     +             ,TRDFLB(L),TRUFLB(L),TRNFLB(L),TRFCRL(L)             5149.   
     +             ,SRDFLB(L),SRUFLB(L),SRNFLB(L),SRFHRL(L)             5150.   
     +             ,STNFLB,STFHR,STDHR,TRDCR,SRDHR,SRALB                5151.   
  301 CONTINUE                                                          5152.   
C                                                                       5153.   
      WRITE(6,6304) BTEMPW,TRUFTW,SRIVIS,SROVIS,PLAVIS,SRINIR,SRONIR    5154.   
     +             ,PLANIR                                              5155.   
      WRITE(6,6305) TRDFGW,TRUFGW,SRDVIS,SRUVIS,ALBVIS,SRDNIR,SRUNIR    5156.   
     +             ,ALBNIR                                              5157.   
      WRITE(6,6306) SRXVIS,SRXNIR,SRTVIS,SRRVIS,SRAVIS,SRTNIR,SRRNIR    5158.   
     +             ,SRANIR                                              5159.   
      WRITE(6,6307) TRDFSL,TRUFSL,TRSLCR,TRSLTS,TRSLTG,TRSLWV,TRSLBS    5160.   
     +             ,SRSLHR                                              5161.   
C                                                                       5162.   
      WRITE(6,6308) (FSRNFG(I),I=1,4),LTOPCL,JLAT,JYEAR                 5163.   
      WRITE(6,6309) (FTRUFG(I),I=1,4),LBOTCL,ILON,JDAY                  5164.   
      WRITE(6,6310) (DTRUFG(I),I=1,4),TTRUFG,COSZ                       5165.   
C                                                                       5166.   
 6301 FORMAT(1I1,'(3) RADCOM M/S: (OUTPUT DATA)'                        5167.   
     +      ,T37,'THERMAL FLUXES (W/M**2)',4X,'SOLAR FLUXES (W/M**2)'   5168.   
     +      ,1X,'NORMS0=',I1,'  ENERGY  INPUT  HEAT/COOL DEG/DAY ALB'   5169.   
     +      ,'DO'/1X,31('-'),2X,9('---'),2X,10('---'),1X,'$',7('-')     5170.   
     +      ,'$',5('-'),1X,'$',5('-'),'$',5('-'),'$',5('-'),1X,'$----'  5171.   
     +      /' LN     PLB   HLB    TLB    TLT '                         5172.   
     +      ,'  TRDFLB TRUFLB TRNFLB TRFCRL   SRDFLB  SRUFLB  SRNFLB'   5173.   
     +      ,' SRFHRL  STNFLB  STFHR  STDHR TRDCR SRDHR SRALB')         5174.   
 6302 FORMAT(1X,I2,F9.3,F6.2,1X,F6.2,8X,3F7.2,8X,3F8.2,7X,F8.2)         5175.   
 6303 FORMAT(1X,I2,F9.3,F6.2,2F7.2,1X,3F7.2,F7.2,1X,3F8.2,F7.2,1X,F7.2  5176.   
     +      ,1X,F6.2,1X,3F6.2,1X,F5.4)                                  5177.   
 6304 FORMAT(/1X,'AT ATM TOP:    ',' BTEMPW=',F6.2,1X,' TRUFTW=',F6.3   5178.   
     +      , 2X,' SRIVIS=',F6.2,' SROVIS=',F6.2,   ' PLAVIS=',F6.4     5179.   
     +      , 2X,' SRINIR=',F6.2,' SRONIR=',F6.2,   ' PLANIR=',F6.4)    5180.   
 6305 FORMAT( 1X,'AT GROUND :    ',' TRDFGW=',F6.3,1X,' TRUFGW=',F6.3   5181.   
     +      , 2X,' SRDVIS=',F6.2,' SRUVIS=',F6.2,   ' ALBVIS=',F6.4     5182.   
     +      , 2X,' SRDNIR=',F6.2,' SRUNIR=',F6.2,   ' ALBNIR=',F6.4)    5183.   
 6306 FORMAT( 1X,'ATMOSPHERE:    ',' SRXVIS=',F6.4,1X,' SRXNIR=',F6.4   5184.   
     +      , 2X,' SRTVIS=',F6.4,' SRRVIS=',F6.4,   ' SRAVIS=',F6.4     5185.   
     +      , 2X,' SRTNIR=',F6.4,' SRRNIR=',F6.4,   ' SRANIR=',F6.4)    5186.   
 6307 FORMAT( 1X,'SURF LAYER:    ',' TRDRSL=',F6.2,1X,' TRUFSL=',F6.2   5187.   
     +      , 2X,' TRSLCR=',F6.4,'+TRSLTS=',F6.4,   '-TRSLTG=',F6.4     5188.   
     +      , 2X,' TRSLWV=',F6.4,' TRSLBS=',F6.3,   ' SRSLHR=',F6.4)    5189.   
 6308 FORMAT(/1X,'FSRNFG(I)=>     FRAC SRNFLB(1) EACH SURFTYPE'         5190.   
     +      ,'   OCEAN=',F7.4,' EARTH=',F7.4,'  OICE=',F7.4,'   LICE='  5191.   
     +      ,F7.4,1X,' LTOPCL=',I2,' JLAT=',I2,' JYEAR',I4)             5192.   
 6309 FORMAT( 1X,'FTRUFG(I)=>     FRAC TRUFLB(1) EACH SURFTYPE'         5193.   
     +      ,'   OCEAN=',F7.4,' EARTH=',F7.4,'  OICE=',F7.4,'   LICE='  5194.   
     +      ,F7.4,1X,' LBOTCL=',I2,' ILON=',I2,' JDAY=',I4)             5195.   
 6310 FORMAT( 1X,'DTRUFG(I)=>    DERIV TRUFLB(1) EACH SURFTYPE'         5196.   
     +      ,'   OCEAN=',F7.4,' EARTH=',F7.4,'  OICE=',F7.4,'   LICE='  5197.   
     +      ,F7.4,  '=>TTRUFG=',F6.4,'   COSZ=',F6.4)                   5198.   
      GO TO 9999                                                        5199.   
C                                                                       5200.   
C-------------                                                          5201.   
  400 CONTINUE                                                          5202.   
C-------------                                                          5203.   
      GO TO 9999                                                        5204.   
C                                                                       5205.   
C-------------                                                          5206.   
  500 CONTINUE                                                          5207.   
C-------------                                                          5208.   
C                                                                       5209.   
      NPAGE=1                                                           5210.   
      IF(INDEX.LT.11) NPAGE=KPAGE                                       5211.   
      SIGMA=5.6697D-08                                                  5212.   
      TGMEAN=POCEAN*TGO**4+PEARTH*TGE**4+PLICE*TGLI**4+POICE*TGOI**4    5213.   
      TGMEAN=SQRT(TGMEAN)                                               5214.   
      TGMEAN=SQRT(TGMEAN)                                               5215.   
      SIGT4=SIGMA*TGMEAN**4                                             5216.   
      ITG=TGMEAN                                                        5217.   
      WTG=TGMEAN-ITG                                                    5218.   
      ITG=ITG-IT0                                                       5219.   
      SUMK=0.0                                                          5220.   
      DO 501 K=1,NKTR                                                   5221.   
      BGFLUX(K)=PLANCK(ITG)-(PLANCK(ITG)-PLANCK(ITG+1))*WTG             5222.   
      BGFRAC(K)=BGFLUX(K)/SIGT4                                         5223.   
      SUMK=SUMK+BGFLUX(K)                                               5224.   
      ITG=ITG+ITNEXT                                                    5225.   
  501 CONTINUE                                                          5226.   
      WRITE(6,6501) NPAGE                                               5227.   
      WRITE(6,6502) (K,K=1,11)                                          5228.   
      DO 502 N=1,NL                                                     5229.   
      L=NLP-N                                                           5230.   
      LI=L                                                              5231.   
      LL=NL*10+L                                                        5232.   
      WRITE(6,6503) L,PL(L),DPL(L),TLM(L),(TAULAP(I),I=LI,LL,NL)        5233.   
  502 CONTINUE                                                          5234.   
      LK=0                                                              5235.   
      DO 504 K=1,NKTR                                                   5236.   
      TAUSUM(K)=0.                                                      5237.   
      DO 503 L=1,NL                                                     5238.   
      LK=LK+1                                                           5239.   
  503 TAUSUM(K)=TAUSUM(K)+TAULAP(LK)                                    5240.   
  504 CONTINUE                                                          5241.   
      WRITE(6,6504) (TAUSUM(K),K=1,11)                                  5242.   
      WRITE(6,6505)                                                     5243.   
      WRITE(6,6506)  SUMK,(BGFLUX(K),K=1,11)                            5244.   
      WRITE(6,6507) TGMEAN,SIGT4,(BGFRAC(K),K=1,11)                     5245.   
      NPAGE=0                                                           5246.   
      IF(NL.GT.13) NPAGE=1                                              5247.   
      WRITE(6,6508) NPAGE                                               5248.   
      WRITE(6,6509) (K,K=12,25)                                         5249.   
      DO 505 N=1,NL                                                     5250.   
      L=NLP-N                                                           5251.   
      LI=NL*11+L                                                        5252.   
      LL=NL*24+L                                                        5253.   
      WRITE(6,6510) L,(TAULAP(I),I=LI,LL,NL)                            5254.   
  505 CONTINUE                                                          5255.   
      WRITE(6,6511) (TAUSUM(K),K=12,NKTR)                               5256.   
      WRITE(6,6512) (BGFLUX(K),K=12,NKTR)                               5257.   
      WRITE(6,6513) (BGFRAC(K),K=12,NKTR)                               5258.   
C                                                                       5259.   
 6501 FORMAT(1I1,'(5) TAULAP TABLE FOR THERMAL RADIATION: INCLUDES'     5260.   
     +      ,' WEAK OVERLAPPING GAS ABSORPTION BY'                      5261.   
     +      ,' H2O, CO2, O3, N2O, CH4',T117,'LIST:  TAULAP(LK)'/        5262.   
     +      ,/1X,'K-DISTRIBUTION BREAKDOWN:'                            5263.   
     +       ,T31,'WINDOW',T65,'WATER VAPOR:  PRINCIPAL ABSORBER REGION'5264.   
     +      ,/T30,8('-'),3X,93('-'))                                    5265.   
 6502 FORMAT(' LN     PL     DPL    TLM    K='                          5266.   
     +      ,I4,5X,'K=',I4,I10,5I9,3I10)                                5267.   
 6503 FORMAT(1X,I2,2F8.3,F7.2,F11.6,F11.6,F10.6,5F9.5,3F10.5)           5268.   
 6504 FORMAT(/13X,'COLUMN AMOUNT=',F10.6,F11.6,F10.6,5F9.5,3F10.5)      5269.   
 6505 FORMAT(' K-INTERVAL CONTRIBUTIONS:'/' COMPARE WITH GROUND FLUX:') 5270.   
 6506 FORMAT( 1X,'PLANCK FLUX W/M**2=',F6.2,2F11.3,F10.3,5F9.3,3F10.3)  5271.   
 6507 FORMAT( 1X,'FRAC(TG=',F6.2,')**4=',F6.2,2X                        5272.   
     +      ,F9.5,F11.5,F10.5,5F9.5,3F10.5)                             5273.   
 6508 FORMAT(1I1/T25,'CARBON DIOXIDE:   PRINCIPAL ABSORBER REGION'      5274.   
     +      ,T100,'OZONE:   PRINCIPAL ABSORBER REGION'                  5275.   
     +      /4X,92('-'),3X,34('-'))                                     5276.   
 6509 FORMAT(1X,'LN K=',I4,I10,6I9,2I10,5X,'K=',I3,3I9)                 5277.   
 6510 FORMAT( 1X,  I2,F9.6,F10.6,6F9.5,2F10.5,F10.5,3F9.5)              5278.   
 6511 FORMAT(/1X,'CA',F9.6,F10.6,6F9.5,2F10.5,F10.5,3F9.5)              5279.   
 6512 FORMAT(/1X,'PF',F9.3,F10.3,6F9.3,2F10.3,F10.3,3F9.3)              5280.   
 6513 FORMAT( 1X,'FR',F9.5,F10.5,6F9.5,2F10.5,F10.5,3F9.5)              5281.   
      GO TO 9999                                                        5282.   
C                                                                       5283.   
C-------------                                                          5284.   
  600 CONTINUE                                                          5285.   
C-------------                                                          5286.   
C                                                                       5287.   
      NPAGE=1                                                           5288.   
      IF(INDEX.LT.11) NPAGE=KPAGE                                       5289.   
      SIGMA=5.6697D-08                                                  5290.   
      TGMEAN=POCEAN*TGO**4+PEARTH*TGE**4+PLICE*TGLI**4+POICE*TGOI**4    5291.   
      TGMEAN=SQRT(TGMEAN)                                               5292.   
      TGMEAN=SQRT(TGMEAN)                                               5293.   
      SIGT4=SIGMA*TGMEAN**4                                             5294.   
      ITG=TGMEAN                                                        5295.   
      WTG=TGMEAN-ITG                                                    5296.   
      ITG=ITG-IT0                                                       5297.   
      SUMK=0.0                                                          5298.   
      DO 601 K=1,NKTR                                                   5299.   
      BGFLUX(K)=PLANCK(ITG)-(PLANCK(ITG)-PLANCK(ITG+1))*WTG             5300.   
      BGFRAC(K)=BGFLUX(K)/SIGT4                                         5301.   
      SUMK=SUMK+BGFLUX(K)                                               5302.   
      ITG=ITG+ITNEXT                                                    5303.   
  601 CONTINUE                                                          5304.   
      WRITE(6,6601) NPAGE                                               5305.   
      WRITE(6,6602) (K,K=1,11)                                          5306.   
      DO 602 N=1,NL                                                     5307.   
      L=NLP-N                                                           5308.   
      LI=L                                                              5309.   
      LL=NL*10+L                                                        5310.   
      WRITE(6,6603) L,PL(L),DPL(L),TLM(L),(TAUN(I),I=LI,LL,NL)          5311.   
  602 CONTINUE                                                          5312.   
      LK=0                                                              5313.   
      DO 604 K=1,NKTR                                                   5314.   
      TAUSUM(K)=TAUSL(K)                                                5315.   
      DO 603 L=1,NL                                                     5316.   
      LK=LK+1                                                           5317.   
  603 TAUSUM(K)=TAUSUM(K)+TAUN(LK)                                      5318.   
  604 CONTINUE                                                          5319.   
      WRITE(6,6604) (TAUSL(K),K=1,11)                                   5320.   
      WRITE(6,6605) (TAUSUM(K),K=1,11)                                  5321.   
      WRITE(6,6606)  SUMK,(BGFLUX(K),K=1,11)                            5322.   
      WRITE(6,6607) TGMEAN,SIGT4,(BGFRAC(K),K=1,11)                     5323.   
      NPAGE=0                                                           5324.   
      IF(NL.GT.13) NPAGE=1                                              5325.   
      WRITE(6,6608) NPAGE                                               5326.   
      WRITE(6,6609) (K,K=12,25)                                         5327.   
      DO 605 N=1,NL                                                     5328.   
      L=NLP-N                                                           5329.   
      LI=NL*11+L                                                        5330.   
      LL=NL*24+L                                                        5331.   
      WRITE(6,6610) L,(TAUN(I),I=LI,LL,NL)                              5332.   
  605 CONTINUE                                                          5333.   
      WRITE(6,6611) ( TAUSL(K),K=12,NKTR)                               5334.   
      WRITE(6,6612) (TAUSUM(K),K=12,NKTR)                               5335.   
      WRITE(6,6613) (BGFLUX(K),K=12,NKTR)                               5336.   
      WRITE(6,6614) (BGFRAC(K),K=12,NKTR)                               5337.   
C                                                                       5338.   
 6601 FORMAT(1I1,'(6) TAU TABLE FOR THERMAL RADIATION: INCLUDES ANY'    5339.   
     +      ,' SPECIFIED OVERLAP, CLOUD & AEROSOL ABSORPTION'           5340.   
     +      ,T117,'TAUN(LK),TAUSL(L)'/                                  5341.   
     +      ,/1X,'K-DISTRIBUTION BREAKDOWN:'                            5342.   
     +       ,T31,'WINDOW',T65,'WATER VAPOR:  PRINCIPAL ABSORBER REGION'5343.   
     +      ,/T30,8('-'),3X,93('-'))                                    5344.   
 6602 FORMAT(' LN     PL     DPL    TLM    K='                          5345.   
     +      ,I4,5X,'K=',I4,I10,5I9,3I10)                                5346.   
 6603 FORMAT(1X,I2,2F8.3,F7.2,F11.6,F11.6,F10.6,5F9.5,3F10.5)           5347.   
 6604 FORMAT(/13X,'SURFACE LAYER=',F10.6,F11.6,F10.6,5F9.5,3F10.5)      5348.   
 6605 FORMAT(/13X,'COLUMN AMOUNT=',F10.3,F11.3,F10.3,5F9.3,3F10.3)      5349.   
 6606 FORMAT(/1X,'PLANCK FLUX W/M**2=',F6.2,2F11.3,F10.3,5F9.3,3F10.3)  5350.   
 6607 FORMAT( 1X,'FRAC(TG=',F6.2,')**4=',F6.2,2X                        5351.   
     +      ,F9.5,F11.5,F10.5,5F9.5,3F10.5)                             5352.   
 6608 FORMAT(1I1/T25,'CARBON DIOXIDE:   PRINCIPAL ABSORBER REGION'      5353.   
     +      ,T100,'OZONE:   PRINCIPAL ABSORBER REGION'                  5354.   
     +      /4X,92('-'),3X,34('-'))                                     5355.   
 6609 FORMAT(1X,'LN K=',I4,I10,6I9,2I10,5X,'K=',I3,3I9)                 5356.   
 6610 FORMAT( 1X,  I2,F9.6,F10.6,6F9.5,2F10.5,F10.5,3F9.5)              5357.   
 6611 FORMAT(/1X,'SL',F9.6,F10.6,6F9.5,2F10.5,F10.5,3F9.5)              5358.   
 6612 FORMAT(/1X,'CA',F9.3,F10.3,6F9.3,2F10.3,F10.3,3F9.3)              5359.   
 6613 FORMAT(/1X,'PF',F9.3,F10.3,6F9.3,2F10.3,F10.3,3F9.3)              5360.   
 6614 FORMAT( 1X,'FR',F9.5,F10.5,6F9.5,2F10.5,F10.5,3F9.5)              5361.   
      GO TO 9999                                                        5362.   
C                                                                       5363.   
C-------------                                                          5364.   
  700 CONTINUE                                                          5365.   
C-------------                                                          5366.   
C                                                                       5367.   
      NPAGE=1                                                           5368.   
      IF(INDEX.LT.11) NPAGE=KPAGE                                       5369.   
      SIGMA=5.6697D-08                                                  5370.   
      TGMEAN=POCEAN*TGO**4+PEARTH*TGE**4+PLICE*TGLI**4+POICE*TGOI**4    5371.   
      TGMEAN=SQRT(TGMEAN)                                               5372.   
      TGMEAN=SQRT(TGMEAN)                                               5373.   
      SIGT4=SIGMA*TGMEAN**4                                             5374.   
      ITG=TGMEAN                                                        5375.   
      WTG=TGMEAN-ITG                                                    5376.   
      ITG=ITG-IT0                                                       5377.   
      SUMK=0.0                                                          5378.   
      DO 701 K=1,NKTR                                                   5379.   
      BGFLUX(K)=PLANCK(ITG)-(PLANCK(ITG)-PLANCK(ITG+1))*WTG             5380.   
      BGFRAC(K)=BGFLUX(K)/SIGT4                                         5381.   
      SUMK=SUMK+BGFLUX(K)                                               5382.   
      ITG=ITG+ITNEXT                                                    5383.   
  701 CONTINUE                                                          5384.   
      WRITE(6,6701) NPAGE                                               5385.   
      WRITE(6,6702) (K,K=1,11)                                          5386.   
      DO 702 N=1,NL                                                     5387.   
      L=NLP-N                                                           5388.   
      WRITE(6,6703) L,PL(L),DPL(L),TLM(L),(TRAEXT(L,K),K=1,11)          5389.   
  702 CONTINUE                                                          5390.   
      DO 704 K=1,NKTR                                                   5391.   
      TAUSUM(K)=0.                                                      5392.   
      DO 703 L=1,NL                                                     5393.   
  703 TAUSUM(K)=TAUSUM(K)+TRAEXT(L,K)                                   5394.   
  704 CONTINUE                                                          5395.   
      WRITE(6,6704) (TAUSUM(K),K=1,11)                                  5396.   
      WRITE(6,6705)                                                     5397.   
      WRITE(6,6706)         SUMK,(BGFLUX(K),K=1,11)                     5398.   
      WRITE(6,6707) TGMEAN,SIGT4,(BGFRAC(K),K=1,11)                     5399.   
      NPAGE=0                                                           5400.   
      IF(NL.GT.13) NPAGE=1                                              5401.   
      WRITE(6,6708) NPAGE                                               5402.   
      WRITE(6,6709) (K,K=12,25)                                         5403.   
      DO 705 N=1,NL                                                     5404.   
      L=NLP-N                                                           5405.   
      WRITE(6,6710) L,(TRAEXT(L,K),K=12,NKTR)                           5406.   
  705 CONTINUE                                                          5407.   
      WRITE(6,6711) (TAUSUM(K),K=12,NKTR)                               5408.   
      WRITE(6,6712) (BGFLUX(K),K=12,NKTR)                               5409.   
      WRITE(6,6713) (BGFRAC(K),K=12,NKTR)                               5410.   
C                                                                       5411.   
 6701 FORMAT(1I1,'(7) AEROSOL TAU TABLE FOR THERMAL RADIATION:'         5412.   
     +      ,'  CLOUD & AEROSOL ABSORPTION'                             5413.   
     +      ,T116,'LIST:  TRAEXT(L,K)'/                                 5414.   
     +      ,/1X,'K-DISTRIBUTION BREAKDOWN:'                            5415.   
     +       ,T31,'WINDOW',T65,'WATER VAPOR:  PRINCIPAL ABSORBER REGION'5416.   
     +      ,/T30,8('-'),3X,93('-'))                                    5417.   
 6702 FORMAT(' LN     PL     DPL    TLM    K='                          5418.   
     +      ,I4,5X,'K=',I4,I10,5I9,3I10)                                5419.   
 6703 FORMAT(1X,I2,2F8.3,F7.2,F11.6,F11.6,F10.6,5F9.5,3F10.5)           5420.   
 6704 FORMAT(/13X,'COLUMN AMOUNT=',F10.6,F11.6,F10.6,5F9.5,3F10.5)      5421.   
 6705 FORMAT(' K-INTERVAL CONTRIBUTIONS:'/' COMPARE WITH GROUND FLUX:') 5422.   
 6706 FORMAT( 1X,'PLANCK FLUX W/M**2=',F6.2,2F11.3,F10.3,5F9.3,3F10.3)  5423.   
 6707 FORMAT( 1X,'FRAC(TG=',F6.2,')**4=',F6.2,2X                        5424.   
     +      ,F9.5,F11.5,F10.5,5F9.5,3F10.5)                             5425.   
 6708 FORMAT(1I1/T25,'CARBON DIOXIDE:   PRINCIPAL ABSORBER REGION'      5426.   
     +      ,T100,'OZONE:   PRINCIPAL ABSORBER REGION'                  5427.   
     +      /4X,92('-'),3X,34('-'))                                     5428.   
 6709 FORMAT(1X,'LN K=',I4,I10,6I9,2I10,5X,'K=',I3,3I9)                 5429.   
 6710 FORMAT( 1X,  I2,F9.6,F10.6,6F9.5,2F10.5,F10.5,3F9.5)              5430.   
 6711 FORMAT(/1X,'CA',F9.6,F10.6,6F9.5,2F10.5,F10.5,3F9.5)              5431.   
 6712 FORMAT(/1X,'PF',F9.3,F10.3,6F9.3,2F10.3,F10.3,3F9.3)              5432.   
 6713 FORMAT( 1X,'FR',F9.5,F10.5,6F9.5,2F10.5,F10.5,3F9.5)              5433.   
      GO TO 9999                                                        5434.   
C                                                                       5435.   
C-------------                                                          5436.   
  800 CONTINUE                                                          5437.   
C-------------                                                          5438.   
C                                                                       5439.   
      NPAGE=1                                                           5440.   
      IF(INDEX.LT.11) NPAGE=KPAGE                                       5441.   
      WRITE(6,6801) NPAGE                                               5442.   
      DO 802 K=1,NKSR                                                   5443.   
      SUM1(K)=0.                                                        5444.   
      SUM2(K)=0.                                                        5445.   
      SUM3(K)=0.                                                        5446.   
      DO 801 L=1,NL                                                     5447.   
      SUM1(K)=SUM1(K)+EXTAER(L,K)                                       5448.   
      SUM2(K)=SUM2(K)+SCTAER(L,K)                                       5449.   
      SUM3(K)=SUM3(K)+SCTAER(L,K)*COSAER(L,K)                           5450.   
  801 PI0AER(L,K)=SCTAER(L,K)/(EXTAER(L,K)+1.E-10)                      5451.   
      SUM3(K)=SUM3(K)/(SUM2(K)+1.E-10)                                  5452.   
      SUM0(K)=SUM2(K)/(SUM1(K)+1.E-10)                                  5453.   
  802 CONTINUE                                                          5454.   
      WRITE(6,6802) (K,K=1,6),(K,K=1,6)                                 5455.   
      DO 803 N=1,NL                                                     5456.   
      L=NLP-N                                                           5457.   
      WRITE(6,6803) L,PLB(L),HLB(L)                                     5458.   
     +              ,(EXTAER(L,J),J=1,6),(SCTAER(L,J),J=1,6)            5459.   
  803 CONTINUE                                                          5460.   
      WRITE(6,6804) (SUM1(K),K=1,NKSR),(SUM2(K),K=1,NKSR)               5461.   
      NPAGE=0                                                           5462.   
      IF(NL.GT.13) NPAGE=1                                              5463.   
      WRITE(6,6805) NPAGE                                               5464.   
      WRITE(6,6806) (K,K=1,6),(K,K=1,6)                                 5465.   
      DO 804 N=1,NL                                                     5466.   
      L=NLP-N                                                           5467.   
      WRITE(6,6807) L,PL(L),DPL(L)                                      5468.   
     +              ,(COSAER(L,J),J=1,6),(PI0AER(L,J),J=1,6)            5469.   
  804 CONTINUE                                                          5470.   
      WRITE(6,6808) (SUM3(K),K=1,NKSR),(SUM0(K),K=1,NKSR)               5471.   
      WRITE(6,6809) (SRBALB(K),K=1,NKSR)                                5472.   
      WRITE(6,6810) (SRXALB(K),K=1,NKSR)                                5473.   
      WRITE(6,6811)                                                     5474.   
      SUM=0.                                                            5475.   
      DO 806 J=1,5                                                      5476.   
      TAU55=0.                                                          5477.   
      DO 805 I=1,NAERO                                                  5478.   
  805 TAU55=TAU55+AGOLDH(I,J)*FGOLDH(J)                                 5479.   
      WRITE(6,6812) J,FGOLDH(J),TAU55                                   5480.   
  806 SUM=SUM+TAU55                                                     5481.   
      WRITE(6,6813) SUM                                                 5482.   
C                                                                       5483.   
 6801 FORMAT(1I1,'(8) AEROSOL INPUT FOR SOLAR RADIATION:'               5484.   
     +      ,'  AEROSOL RADIATIVE PROPERTIES'                           5485.   
     +      ,T81,'LIST: EXTAER(L,K),SCTAER(L,K),COSAER(L,K),PIZERO(L,K)'5486.   
     +      //T42,'TAU -- EXTINCTION',T99,'TAU -- SCATTERING'           5487.   
     +      ,/T24,53('-'),4X,53('-'))                                   5488.   
 6802 FORMAT(' LN    PLB     HLB     K=',I3,5I9,7X,'K=',I3,5I9)         5489.   
 6803 FORMAT(1X,I2,2F8.3,3X,6F9.6,3X,6F9.6)                             5490.   
 6804 FORMAT(/1X,T7,'COLUMN AMOUNT=',2X,6F9.6,3X,6F9.6)                 5491.   
 6805 FORMAT(1I1/T48,'COSBAR',T105,'PIZERO'                             5492.   
     +      ,/T24,53('-'),4X,53('-'))                                   5493.   
 6806 FORMAT(' LN     PL     DPL     K=',I3,5I9,7X,'K=',I3,5I9)         5494.   
 6807 FORMAT(1X,I2,2F8.3,3X,6F9.6,3X,6F9.6)                             5495.   
 6808 FORMAT(/1X,T7,'COLUMN   MEAN=',2X,6F9.6,3X,6F9.6)                 5496.   
 6809 FORMAT(/1X,T7,'ALBEDO RSURFB=',2X,6F9.6,3X,6F9.6)                 5497.   
 6810 FORMAT( 1X,T7,'ALBEDO RSURFX=',2X,6F9.6,3X,6F9.6)                 5498.   
      GO TO 9999                                                        5499.   
 6811 FORMAT(///T44,'AEROSOL COMPOSITION AND TYPE MIX:'                 5500.   
     +      ,T81,'FACTOR',6X,'VALUE',T107,'TAU(0.55)'/)                 5501.   
 6812 FORMAT(T81,'FGOLDH(',I1,') =',1P,E9.2,5X,0P,F7.4)                 5502.   
 6813 FORMAT(/T81,'SUM COLUMN TAU(0.55) =',F10.4)                       5503.   
C                                                                       5504.   
C-------------                                                          5505.   
  900 CONTINUE                                                          5506.   
C-------------                                                          5507.   
C                                                                       5508.   
      SIGMA=5.6697D-08                                                  5509.   
      TGMEAN=POCEAN*TGO**4+PEARTH*TGE**4+PLICE*TGLI**4+POICE*TGOI**4    5510.   
      TGMEAN=SQRT(TGMEAN)                                               5511.   
      TGMEAN=SQRT(TGMEAN)                                               5512.   
      SIGT4=SIGMA*TGMEAN**4                                             5513.   
      ITG=TGMEAN                                                        5514.   
      WTG=TGMEAN-ITG                                                    5515.   
      ITG=ITG-IT0                                                       5516.   
      DO 901 K=1,NKTR                                                   5517.   
      BGFLUX(K)=PLANCK(ITG)-(PLANCK(ITG)-PLANCK(ITG+1))*WTG             5518.   
      BGFRAC(K)=BGFLUX(K)/SIGT4                                         5519.   
      ITG=ITG+ITNEXT                                                    5520.   
  901 CONTINUE                                                          5521.   
      DO 910 NW=1,5                                                     5522.   
      DO 903 K=1,NKTR                                                   5523.   
      DO 902 L=1,NLP                                                    5524.   
      IF(NW.EQ.1) WFLB(L,K)=DFLB(L,K)                                   5525.   
      IF(NW.EQ.2) WFLB(L,K)=UFLB(L,K)                                   5526.   
      IF(NW.EQ.3) WFLB(L,K)=UFLB(L,K)-DFLB(L,K)                         5527.   
      IF(NW.GT.3.AND.L.GT.NL) GO TO 902                                 5528.   
      IF(NW.EQ.4) WFLB(L,K)=WFLB(L+1,K)-WFLB(L,K)                       5529.   
      IF(NW.EQ.5.AND.ABS(TRFCRL(L)).LT.1.E-10) WFLB(L,K)=1.E-30         5530.   
      IF(NW.EQ.5) WFLB(L,K)=WFLB(L,K)/(ABS(TRFCRL(L))+1.E-10)           5531.   
  902 CONTINUE                                                          5532.   
      IF(NW.EQ.1) WFSL(K)=DFSL(K)                                       5533.   
      IF(NW.EQ.2) WFSL(K)=UFSL(K)                                       5534.   
      IF(NW.EQ.3) WFSL(K)=UFSL(K)-DFSL(K)                               5535.   
      IF(NW.EQ.4) WFSL(K)=WFSL(K)-UFLB(1,K)+DFLB(1,K)                   5536.   
      IF(NW.EQ.5.AND.ABS(TRSLCR).LT.1.E-10) WFSL(K)=1.E-30              5537.   
      IF(NW.EQ.5) WFSL(K)=WFSL(K)/(ABS(TRSLCR)+1.E-10)                  5538.   
  903 CONTINUE                                                          5539.   
      DO 907 L=1,NLP                                                    5540.   
      IF(L.GT.NL.AND.NW.GT.3) GO TO 907                                 5541.   
      ASUM1=0.                                                          5542.   
      BSUM1=0.                                                          5543.   
      CSUM1=0.                                                          5544.   
      DSUM1=0.                                                          5545.   
      ESUM1=0.                                                          5546.   
      FSUM1=0.                                                          5547.   
      SUM=0.                                                            5548.   
      DO 904 K=2,11                                                     5549.   
      ASUM1=ASUM1+  WFSL(K)                                             5550.   
      BSUM1=BSUM1+ BGFEMT(K)                                            5551.   
      CSUM1=CSUM1+BGFLUX(K)                                             5552.   
      DSUM1=DSUM1+BGFRAC(K)                                             5553.   
      ESUM1=ESUM1+TRCALB(K)                                             5554.   
      FSUM1=FSUM1+ TRGALB(K)                                            5555.   
  904 SUM=SUM+WFLB(L,K)                                                 5556.   
      SUM1(L)=SUM                                                       5557.   
      ASUM2=0.                                                          5558.   
      BSUM2=0.                                                          5559.   
      CSUM2=0.                                                          5560.   
      DSUM2=0.                                                          5561.   
      ESUM2=0.                                                          5562.   
      FSUM2=0.                                                          5563.   
      SUM=0.                                                            5564.   
      DO 905 K=12,21                                                    5565.   
      ASUM2=ASUM2+  WFSL(K)                                             5566.   
      BSUM2=BSUM2+ BGFEMT(K)                                            5567.   
      CSUM2=CSUM2+BGFLUX(K)                                             5568.   
      DSUM2=DSUM2+BGFRAC(K)                                             5569.   
      ESUM2=ESUM2+TRCALB(K)                                             5570.   
      FSUM2=FSUM2+ TRGALB(K)                                            5571.   
  905 SUM=SUM+WFLB(L,K)                                                 5572.   
      SUM2(L)=SUM                                                       5573.   
      ASUM3=0.                                                          5574.   
      BSUM3=0.                                                          5575.   
      CSUM3=0.                                                          5576.   
      DSUM3=0.                                                          5577.   
      ESUM3=0.                                                          5578.   
      FSUM3=0.                                                          5579.   
      SUM=0.                                                            5580.   
      DO 906 K=22,NKTR                                                  5581.   
      ASUM3=ASUM3+  WFSL(K)                                             5582.   
      BSUM3=BSUM3+ BGFEMT(K)                                            5583.   
      CSUM3=CSUM3+BGFLUX(K)                                             5584.   
      DSUM3=DSUM3+BGFRAC(K)                                             5585.   
      ESUM3=ESUM3+TRCALB(K)                                             5586.   
      FSUM3=FSUM3+ TRGALB(K)                                            5587.   
  906 SUM=SUM+WFLB(L,K)                                                 5588.   
      SUM3(L)=SUM                                                       5589.   
  907 CONTINUE                                                          5590.   
C                                                                       5591.   
      NPAGE=1                                                           5592.   
      WRITE(6,6901) NPAGE,NW,FTYPE(NW)                                  5593.   
      WRITE(6,6902) (K,K=1,11)                                          5594.   
      DO 908 N=1,NLP                                                    5595.   
      L=NLP+1-N                                                         5596.   
      IF(L.GT.NL.AND.NW.GT.3) GO TO 908                                 5597.   
      SUML=SUM1(L)+SUM2(L)+SUM3(L)+WFLB(L,1)                            5598.   
      WRITE(6,6903) L,SUML,SUM1(L),SUM2(L),SUM3(L),(WFLB(L,K),K=1,11)   5599.   
  908 CONTINUE                                                          5600.   
      SUMA=ASUM1+ASUM2+ASUM3+  WFSL(1)                                  5601.   
      SUMB=BSUM1+BSUM2+BSUM3+ BGFEMT(1)                                 5602.   
      SUMC=CSUM1+CSUM2+CSUM3+BGFLUX(1)                                  5603.   
      SUMD=DSUM1+DSUM2+DSUM3+BGFRAC(1)                                  5604.   
      SUME=ESUM1+ESUM2+ESUM3+TRCALB(1)                                  5605.   
      SUMF=FSUM1+FSUM2+FSUM3+TRGALB(1)                                  5606.   
      WRITE(6,6904) SUMA,ASUM1,ASUM2,ASUM3,(  WFSL(K),K=1,11)           5607.   
      WRITE(6,6905) SUMB,BSUM1,BSUM2,BSUM3,( BGFEMT(K),K=1,11)          5608.   
      WRITE(6,6906) SUMC,CSUM1,CSUM2,CSUM3,(BGFLUX(K),K=1,11)           5609.   
      WRITE(6,6907) SUMD,DSUM1,DSUM2,DSUM3,(BGFRAC(K),K=1,11)           5610.   
      WRITE(6,6908) SUME,ESUM1,ESUM2,ESUM3,(TRCALB(K),K=1,11)           5611.   
      WRITE(6,6909) SUMF,FSUM1,FSUM2,FSUM3,(TRGALB(K),K=1,11)           5612.   
      NPAGE=0                                                           5613.   
      IF(NL.GT.13)  NPAGE=1                                             5614.   
      WRITE(6,6910) NPAGE                                               5615.   
      WRITE(6,6911) (K,K=12,25)                                         5616.   
      DO 909 N=1,NLP                                                    5617.   
      L=NLP+1-N                                                         5618.   
      IF(L.GT.NL.AND.NW.GT.3) GO TO 909                                 5619.   
      WRITE(6,6912) L,(WFLB(L,K),K=12,NKTR)                             5620.   
  909 CONTINUE                                                          5621.   
      WRITE(6,6913) (  WFSL(K),K=12,NKTR)                               5622.   
      WRITE(6,6914) ( BGFEMT(K),K=12,NKTR)                              5623.   
      WRITE(6,6915) (BGFLUX(K),K=12,NKTR)                               5624.   
      WRITE(6,6916) (BGFRAC(K),K=12,NKTR)                               5625.   
      WRITE(6,6917) (TRCALB(K),K=12,NKTR)                               5626.   
      WRITE(6,6918) ( TRGALB(K),K=12,NKTR)                              5627.   
  910 CONTINUE                                                          5628.   
C                                                                       5629.   
 6901 FORMAT(1I1,'(9.',I1,') THERMAL RADIATION: K-DISTRIBUTION'         5630.   
     +      ,' BREAKDOWN FOR  ',1A8,' FLUX'/                            5631.   
     +       /T8,'SUM   PRINCIPAL REGION SUM',4X                        5632.   
     +       ,'WINDOW',T66,'WATER VAPOR:  PRINCIPAL ABSORBER REGION'    5633.   
     +      ,/T7,'-----',2X,20('-'),4X,6('-'),3X,87('-'))               5634.   
 6902 FORMAT(1X,'LN   TOTAL    H2O    CO2     O3     K='                5635.   
     +      ,I2,5X,'K=',I2,9I9)                                         5636.   
 6903 FORMAT( 1X,I2,F8.2,1X,3F7.2,F10.3,10F9.3)                         5637.   
 6904 FORMAT(/' SL',F8.2,1X,3F7.2,F10.3,10F9.3)                         5638.   
 6905 FORMAT(/' BG',F8.2,1X,3F7.2,F10.3,10F9.3)                         5639.   
 6906 FORMAT( ' PF',F8.2,1X,3F7.2,F10.3,10F9.3)                         5640.   
 6907 FORMAT( ' FR',F8.4,1X,3F7.4,F10.5,10F9.5)                         5641.   
 6908 FORMAT(/' AC',F8.2,1X,3F7.2,F10.3,10F9.3)                         5642.   
 6909 FORMAT( ' AG',F8.2,1X,3F7.2,F10.3,10F9.3)                         5643.   
 6910 FORMAT(1I1/T26,'CARBON DIOXIDE:  PRINCIPAL ABSORBER REGION'       5644.   
     +      ,T100,'OZONE:   PRINCIPAL ABSORBER REGION'                  5645.   
     +      /5X,89('-'),5X,34('-'))                                     5646.   
 6911 FORMAT(1X,'LN  K=',I4,9I9,7X,'K=',I3,3I9)                         5647.   
 6912 FORMAT( 1X,I2,1X,10F9.3,3X,4F9.3)                                 5648.   
 6913 FORMAT(/' SL',1X,10F9.3,3X,4F9.3)                                 5649.   
 6914 FORMAT(/' BG',1X,10F9.3,3X,4F9.3)                                 5650.   
 6915 FORMAT( ' PF',1X,10F9.3,3X,4F9.3)                                 5651.   
 6916 FORMAT( ' FR',1X,10F9.5,3X,4F9.5)                                 5652.   
 6917 FORMAT(/' AC',1X,10F9.3,3X,4F9.3)                                 5653.   
 6918 FORMAT( ' AG',1X,10F9.3,3X,4F9.3)                                 5654.   
      RETURN                                                            5655.   
C                                                                       5656.   
C-------------                                                          5657.   
 1000 CONTINUE                                                          5658.   
C-------------                                                          5659.   
C                                                                       5660.   
      NPAGE=1                                                           5661.   
      IF(INDEX.LT.11) NPAGE=KPAGE                                       5662.   
      WRITE(6,7001) NPAGE                                               5663.   
 7001 FORMAT(1I1,'(10) BLOCK DATA AEROSOL PROPERTY SPECIFICATION:')     5664.   
 9999 CONTINUE                                                          5665.   
      RETURN                                                            5666.   
      END                                                               5667.   
      SUBROUTINE SOLARZ(NG,KWRITE)                                      5668.   
#include "B83XX.COM"                                                    5669.   
      DIMENSION SRDATA(187),ZRDATA(187)                                 5730.   
      EQUIVALENCE (SRDFLB(1),SRDATA(1))                                 5731.   
c     DOUBLE PRECISION XMU(50),WT(50)                                   5732.   
      dimension XMU(50),WT(50)
      DATA NSRD/187/                                                    5733.   
      DIMENSION NOFLUX(7)                                               5734.   
      DATA NOFLUX/164,167,168,169,170,171,174/                          5735.   
C                                                                       5736.   
C-------------------------------------                                  5737.   
      CALL GAUSST(NG,0.D0,1.D0,XMU,WT)                                  5738.   
C-------------------------------------                                  5739.   
      DO 100 J=1,NG                                                     5740.   
  100 WT(J)=WT(J)*2.D0*XMU(J)                                           5741.   
C                                                                       5742.   
      DO 110 I=1,NSRD                                                   5743.   
  110 ZRDATA(I)=0.                                                      5744.   
C                                                                       5745.   
      NORM=NORMS0                                                       5746.   
      ZCOS=COSZ                                                         5747.   
C                                                                       5748.   
      DO 130 J=1,NG                                                     5749.   
      COSZ=XMU(J)                                                       5750.   
      NORMS0=1                                                          5751.   
C---------------                                                        5752.   
      CALL SOLAR                                                        5753.   
C---------------                                                        5754.   
      DO 120 I=1,NSRD                                                   5755.   
  120 ZRDATA(I)=ZRDATA(I)+SRDATA(I)*WT(J)                               5756.   
      KPAGE=J-(J/2)*2                                                   5757.   
      IF(KWRITE.GT.1) CALL WRITER(3,KPAGE)                              5758.   
  130 CONTINUE                                                          5759.   
C                                                                       5760.   
      DO 150 I=1,NSRD                                                   5761.   
      FACTOR=0.25                                                       5762.   
      DO 140 K=1,7                                                      5763.   
      IF(I.EQ.NOFLUX(K)) FACTOR=1.                                      5764.   
  140 CONTINUE                                                          5765.   
      IF(I.GT.176) FACTOR=1.                                            5766.   
  150 SRDATA(I)=ZRDATA(I)*FACTOR                                        5767.   
      COSZ=NG                                                           5768.   
      IF(NG.GT.9) COSZ=.1*NG                                            5769.   
      COSZ=COSZ+NG/1000.                                                5770.   
      KPAGE=1                                                           5771.   
C                                                                       5772.   
      NORMS0=100                                                        5773.   
C                                                                       5774.   
      IF(KWRITE.GT.0) CALL WRITER(13,KPAGE)                             5775.   
C                                                                       5776.   
      COSZ=ZCOS                                                         5777.   
      NORMS0=NORM                                                       5778.   
C                                                                       5779.   
      RETURN                                                            5780.   
      END                                                               5781.   
      SUBROUTINE GAUSST(NG,X1,X2,XP,WT)                                 5782.   
c     IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               5783.   
      DIMENSION XP(1),WT(1)                                             5784.   
	real*8 pi, ps, dxl
      DATA PI,PS,DXL/3.141592653589793D0,1.013211836423378D-01,1.D-16/  5785.   
      XMID=(X2+X1)/2.D0                                                 5786.   
      XDIF=X2-X1                                                        5787.   
      XHAF=XDIF/2.D0                                                    5788.   
      DNG=NG                                                            5789.   
      NN=NG/2                                                           5790.   
      N2=NN*2                                                           5791.   
      IF(N2.EQ.NG) GO TO 110                                            5792.   
      XP(NN+1)=XMID                                                     5793.   
      WT(NN+1)=XDIF                                                     5794.   
      IF(NG.LT.2) RETURN                                                5795.   
      PN=1.D0                                                           5796.   
      N=0                                                               5797.   
 100  N=N+2                                                             5798.   
      DN=N                                                              5799.   
      DM=DN-1.D0                                                        5800.   
      PN=PN*(DM/DN)                                                     5801.   
      IF(N.LT.N2) GO TO 100                                             5802.   
      WT(NN+1)=XDIF/(DNG*PN)**2                                         5803.   
 110  I=0                                                               5804.   
      C=PI/DSQRT(DNG*(DNG+1.D0)+0.5D0-PS)/105.D0                        5805.   
 120  I=I+1                                                             5806.   
      DI=I                                                              5807.   
      Z=PS/(4.D0*DI-1.D0)**2                                            5808.   
      ZZ=(105.D0+Z*(210.D0-Z*(2170.D0-Z*(105812.D0-12554474.D0*Z))))    5809.   
      X=DCOS(ZZ*C*(DI-0.25D0))                                          5810.   
 130  N=1                                                               5811.   
      DM=1.D0                                                           5812.   
      PNI=1.D0                                                          5813.   
      PNJ=X                                                             5814.   
 140  N=N+1                                                             5815.   
      DN=N                                                              5816.   
      PNK=((DM+DN)*X*PNJ-DM*PNI)/DN                                     5817.   
      PNI=PNJ                                                           5818.   
      PNJ=PNK                                                           5819.   
      DM=DN                                                             5820.   
      IF(N.LT.NG) GO TO 140                                             5821.   
      DX=PNJ*(1.D0-X*X)/DNG/(PNI-X*PNJ)                                 5822.   
      X=X-DX                                                            5823.   
      IF(DABS(DX).GT.DXL) GO TO 130                                     5824.   
      J=NG+1-I                                                          5825.   
      XP(I)=XMID-XHAF*X                                                 5826.   
      XP(J)=XMID+XHAF*X                                                 5827.   
      WT(I)=XDIF*(1.D0-X*X)/(DNG*PNI)**2                                5828.   
      WT(J)=WT(I)                                                       5829.   
      IF(I.LT.NN) GO TO 120                                             5830.   
      RETURN                                                            5831.   
      END                                                               5832.   
      SUBROUTINE SETATM                                                 5833.   
#include "B83XX.COM"                                                    5834.   
      DIMENSION NL4(4),PLB4(40,4)                                       5877.   
      DATA NL4/12,12,24,35/                                             5878.   
      DATA PLB4/                                                        5879.   
     1     1013.2500, 961.7485, 879.3460, 741.3219, 566.2166, 401.4117, 5880.   
     1      262.3575, 154.2043,  71.8018,  10.0000,   5.0000,   2.0000, 5881.   
     1        1.E-05,  27*0.,                                           5882.   
C                                                                       5883.   
     2      984.0000, 934.0000, 854.0000, 720.0000, 550.0000, 390.0000, 5884.   
     2      255.0000, 150.0000,  70.0000,  10.0000,   5.0000,   2.0000, 5885.   
     2        1.E-05,  27*0.,                                           5886.   
C                                                                       5887.   
     3     1013.2500, 988.8846, 956.9068, 910.2775, 820.4963, 683.6775, 5888.   
     3      521.6665, 356.3138, 209.4467, 102.9552,  47.7944,  22.1797, 5889.   
     3      10.29439,  4.77932,  2.21785,  1.01932,  0.46761,  0.21156, 5890.   
     3      0.092671, 0.047500, 0.021885, 0.010000, 0.005000, 0.002000, 5891.   
     3      1.00E-05, 15*0.0,                                           5892.   
C                                                                       5893.   
     4     1013.2500,1000.0000, 950.0000, 900.0000, 850.0000, 800.0000, 5894.   
     4      750.0000, 700.0000, 650.0000, 600.0000, 550.0000, 500.0000, 5895.   
     4      450.0000, 400.0000, 350.0000, 300.0000, 250.0000, 200.0000, 5896.   
     4      150.0000, 100.0000,  50.0000,  20.0000,  10.0000,   5.0000, 5897.   
     4        2.0000,   1.0000,   0.5000,   0.2000,   0.1000,   0.0500, 5898.   
     4        0.0200,   0.0100,   0.0050,   0.0020,   0.0010,   1.E-05, 5899.   
     4        4*0./                                                     5900.   
C                                                                       5901.   
      LAST=LASTVC                                                       5902.   
      LMAG=100000                                                       5903.   
C                          ------------------------------------------   5904.   
C                          NLAY:   ATMOSPHERIC LAYERING SPECIFICATION   5905.   
C                          ------------------------------------------   5906.   
      NLAY=LAST/LMAG                                                    5907.   
      LAST=LAST-LMAG*NLAY                                               5908.   
      LMAG=LMAG/10                                                      5909.   
C                                                                       5910.   
      KSCALE=0                                                          5911.   
      IF(NLAY.GT.9) KSCALE=1                                            5912.   
      IF(NLAY.GT.9) NLAY=NLAY-10                                        5913.   
C                                                                       5914.   
      IF(NLAY.LT.1.OR.NLAY.GT.8)  GO TO 20                              5915.   
      GO TO (10,10,10,10,12,14,16,18),NLAY                              5916.   
   10 NL=NL4(NLAY)                                                      5917.   
      NLP=NL+1                                                          5918.   
C                            (1-4)=(12,12,24,35 PRESSURE SPECIFICATIONS)5919.   
C                            -------------------------------------------5920.   
      DO 11 N=1,NLP                                                     5921.   
   11 PLB(N)=PLB4(N,NLAY)                                               5922.   
      GO TO 20                                                          5923.   
C                                    (5)=(1-D MODEL LAYER SPECIFICATION)5924.   
C                                    -----------------------------------5925.   
   12 NL=18                                                             5926.   
      DO 13 N=1,NL                                                      5927.   
      HLB(N)=N-1+2*(N/7)                                                5928.   
      IF(N.GT. 8) HLB(N)=4*N-24-N/11-N/12                               5929.   
   13 IF(N.GT.13) HLB(N)=30+(N-14)*5                                    5930.   
      HLB( 1)=1.0E-10                                                   5931.   
      HLB(19)=99.99                                                     5932.   
      GO TO 20                                                          5933.   
C                                 (6)=(LINE-BY-LINE LAYER SPECIFICATION)5934.   
C                                 --------------------------------------5935.   
   14 NL=30                                                             5936.   
      DO 15 N=1,NL                                                      5937.   
      HLB(N)=N-1+(N-17)*(N/17)                                          5938.   
   15 IF(N.GT.20) HLB(N)=20+(N-20)*5                                    5939.   
      HLB( 1)=1.0E-10                                                   5940.   
      HLB(31)=99.99                                                     5941.   
      GO TO 20                                                          5942.   
C                                   (7)=(MCCLATCHEY LAYER SPECIFICATION)5943.   
C                                   ------------------------------------5944.   
   16 NL=32                                                             5945.   
      DO 17 N=1,NL                                                      5946.   
      HLB(N)=N-1                                                        5947.   
   17 IF(N.GT.25) HLB(N)=25+5*(N-26)                                    5948.   
      HLB( 1)=1.0E-10                                                   5949.   
      HLB(32)=70.00                                                     5950.   
      HLB(33)=99.99                                                     5951.   
      GO TO 20                                                          5952.   
C                                       (8)=(HI-RES LAYER SPECIFICATION)5953.   
C                                       --------------------------------5954.   
   18 NL=39                                                             5955.   
      DO 19 N=1,NL                                                      5956.   
      HLB(N)=N-1                                                        5957.   
      IF(N.GT.21) HLB(N)=20+(N-21)*2                                    5958.   
      IF(N.GT.31) HLB(N)=40+(N-31)*5                                    5959.   
   19 IF(N.GT.37) HLB(N)=70+(N-37)*10                                   5960.   
      HLB( 1)=1.0E-10                                                   5961.   
      HLB(40)=99.99                                                     5962.   
C                                                                       5963.   
C                         -------------------------------------------   5964.   
C                         NATM:   ATMOSPHERIC STRUCTURE SPECIFICATION   5965.   
C                         -------------------------------------------   5966.   
   20 NATM=LAST/LMAG                                                    5967.   
      LAST=LAST-LMAG*NATM                                               5968.   
      LMAG=LMAG/10                                                      5969.   
C                                                                       5970.   
      IF(KSCALE.NE.1) GO TO 24                                          5971.   
C                                                                       5972.   
C                   SIGMA LEVEL RESCALING OF PRESSURES RELATIVE TO PSIG05973.   
C                   ----------------------------------------------------5974.   
C                                                                       5975.   
      NLMOD=NL-LAYRAD                                                   5976.   
      IF(NLAY.GT.4) GO TO 22                                            5977.   
      PTOP=PLB(NLMOD+1)                                                 5978.   
      PBOT=PLB(1)                                                       5979.   
      DO 21 L=1,NLMOD                                                   5980.   
      PSIG(L)=(PLB(L)-PTOP)/(PBOT-PTOP)                                 5981.   
   21 PLB(L) =PSIG(L)*(PSIG0-PTOP)+PTOP                                 5982.   
      PSIG(NLMOD+1)=0.                                                  5983.   
      GO TO 24                                                          5984.   
C                                                                       5985.   
C                     SIGMA LEVEL RESCALING OF HEIGHTS RELATIVE TO PSIG05986.   
C                     --------------------------------------------------5987.   
   22 HTOP=HLB(NLMOD+1)                                                 5988.   
      HBOT=HLB(1)                                                       5989.   
      DO 23 L=1,NLMOD                                                   5990.   
      PSIG(L)=(HLB(L)-HTOP)/(HBOT-HTOP)                                 5991.   
   23 HLB(L) =PSIG(L)*(PSIG0-HTOP)+HTOP                                 5992.   
      PSIG(NLMOD+1)=0.                                                  5993.   
   24 CONTINUE                                                          5994.   
C                                                                       5995.   
      NLP=NL+1                                                          5996.   
      NPHD=1+NLAY/5                                                     5997.   
      N=1                                                               5998.   
      IF(NPHD.EQ.1) P=PLB(N)                                            5999.   
      IF(NPHD.EQ.2) H=HLB(N)                                            6000.   
      CALL PHDATM(P,H,D,T,O,Q,S,OCM,WCM,NPHD,NATM)                      6001.   
      IF(NPHD.EQ.1) HLB(N)=H                                            6002.   
      IF(NPHD.EQ.2) PLB(N)=P                                            6003.   
      PB=P                                                              6004.   
      TB=T                                                              6005.   
      OB=OCM                                                            6006.   
      WB=WCM                                                            6007.   
      DO 25 N=1,NL                                                      6008.   
      IF(NPHD.EQ.1) P=PLB(N+1)                                          6009.   
      IF(NPHD.EQ.2) H=HLB(N+1)                                          6010.   
      CALL PHDATM(P,H,D,T,O,Q,S,OCM,WCM,NPHD,NATM)                      6011.   
      IF(NPHD.EQ.1) HLB(N+1)=H                                          6012.   
      IF(NPHD.EQ.2) PLB(N+1)=P                                          6013.   
      TLB(N)=TB                                                         6014.   
      TLT(N)=T                                                          6015.   
      TLM(N)=0.5*(T+TB)                                                 6016.   
      U0GAS(N,1)=WB-WCM                                                 6017.   
      U0GAS(N,3)=OB-OCM                                                 6018.   
      SHL(N)=U0GAS(N,1)/(U0GAS(N,1)+1268.75*(PB-P))                     6019.   
      EQ=0.5*(PB+P)*SHL(N)/(0.662+0.338*SHL(N))                         6020.   
C$    EQ=0.5*(PB+P)*SHL(N)/(0.622+0.338*SHL(N))                         6021.   
      ES=10.0**(9.4051-2353.0/TLM(N))                                   6022.   
      RHL(N)=EQ/ES                                                      6023.   
      PB=P                                                              6024.   
      TB=T                                                              6025.   
      OB=OCM                                                            6026.   
   25 WB=WCM                                                            6027.   
      TLB(NLP)=TLT(NL)                                                  6028.   
      TSL=TLB(1)                                                        6029.   
      TGO=TLB(1)                                                        6030.   
      TGE=TLB(1)                                                        6031.   
      TGOI=TGO-5.                                                       6032.   
      TGLI=TGE-5.                                                       6033.   
C                                 ----------------------------------    6034.   
C                                 NSUR:   SURFACE TYPE SPECIFICATION    6035.   
C                                 ----------------------------------    6036.   
   30 NSUR=LAST/LMAG                                                    6037.   
      LAST=LAST-LMAG*NSUR                                               6038.   
      LMAG=LMAG/10                                                      6039.   
C                                                                       6040.   
      IF(NSUR.EQ.0) GO TO 40                                            6041.   
      POCEAN=0.                                                         6042.   
      PEARTH=0.                                                         6043.   
      POICE =0.                                                         6044.   
      PLICE =0.                                                         6045.   
      AGESN =0.                                                         6046.   
      SNOWE =0.                                                         6047.   
      SNOWOI=0.                                                         6048.   
      SNOWLI=0.                                                         6049.   
C                                                                       6050.   
      IF(NSUR.EQ.1) POCEAN=1.                                           6051.   
      IF(NSUR.EQ.2) PEARTH=1.                                           6052.   
      IF(NSUR.EQ.3) POICE =1.                                           6053.   
      IF(NSUR.EQ.4) PLICE =1.                                           6054.   
      IF(NSUR.EQ.5) PEARTH=1.                                           6055.   
      IF(NSUR.EQ.5) SNOWE =1.                                           6056.   
      IF(NSUR.GT.5) PLICE =1.                                           6057.   
      IF(NSUR.EQ.6) SNOWLI=1.                                           6058.   
      IF(NSUR.LT.7) GO TO 40                                            6059.   
      BXAVIS=0.                                                         6060.   
      BXANIR=0.                                                         6061.   
      IF(NSUR.EQ.7) BXAVIS=1.                                           6062.   
      IF(NSUR.GT.7) BXANIR=1.                                           6063.   
      IF(NSUR.EQ.9) BXAVIS=1.                                           6064.   
      DO 31 I=1,5                                                       6065.   
      SRBXAL(I,1)=BXANIR                                                6066.   
   31 SRBXAL(I,2)=BXANIR                                                6067.   
      SRBXAL(6,1)=BXAVIS                                                6068.   
      SRBXAL(6,2)=BXAVIS                                                6069.   
      IF(KALVIS.GT.0) SRBXAL(4,1)=SRBXAL(6,1)                           6070.   
      IF(KALVIS.GT.0) SRBXAL(4,2)=SRBXAL(6,2)                           6071.   
C                                                                       6072.   
C                            ----------------------------------------   6073.   
C                            NTRA:   TRACER COMPOSITION SPECIFICATION   6074.   
C                            ----------------------------------------   6075.   
   40 NTRA=LAST/LMAG                                                    6076.   
      LAST=LAST-LMAG*NTRA                                               6077.   
      LMAG=LMAG/10                                                      6078.   
C                                                                       6079.   
      TAUT55=1.0                                                        6080.   
      NTRACE=1                                                          6081.   
      IF(NTRA.LT.1) TAUT55=0.                                           6082.   
      IF(NTRA.LT.1) NTRACE=0                                            6083.   
      ITR(1)=NTRA                                                       6084.   
      DO 41 L=1,NL                                                      6085.   
   41 TRACER(L,1)=TAUT55*(PLB(L)-PLB(L+1))/PLB(1)                       6086.   
C                                                                       6087.   
C                              -------------------------------------    6088.   
C                              NVEG:   VEGETATION TYPE SPECIFICATION    6089.   
C                              -------------------------------------    6090.   
   50 NVEG=LAST/LMAG                                                    6091.   
      LAST=LAST-LMAG*NVEG                                               6092.   
      LMAG=LMAG/10                                                      6093.   
C                                                                       6094.   
      DO 51 K=1,11                                                      6095.   
   51 PVT(K)=0.                                                         6096.   
      IF(NVEG.LT.1) GO TO 60                                            6097.   
      PVT(NVEG)=1.                                                      6098.   
C                               -------------------------------------   6099.   
C                               NCLD:   CLOUD LAYER,TAU SPECIFICATION   6100.   
C                               -------------------------------------   6101.   
   60 NCLD=LAST                                                         6102.   
      DO 61 L=1,NL                                                      6103.   
   61 CLDTAU(L)=0.                                                      6104.   
      IF(NCLD.GT.0) CLDTAU(NCLD)=64./2**NCLD                            6105.   
      RETURN                                                            6106.   
      END                                                               6107.   
      SUBROUTINE SETFOR(NFTFOR)                                         6108.   
#include "B83XX.COM"                                                    6109.   
C     COMMON/TMINOR/FCO2,FN2O,FCH4,FF11,FF12,FVOL,FSUN                  6150.   
C                                                                       6151.   
C-----------------------------------------------------------------------6152.   
C     EXTERNAL FORCING FOR  CO2,N2O,CH4,F11,F12,VOLCANIC AER,SOLAR CONST6153.   
C     STARTING FROM  JAN 1,1880  PROJECTED THROUGH  DEC 31,2100         6154.   
C     INPUT FORCING DATA READ IN FROM DISK DATA  DSN=CLIM.RUN.FORCING   6155.   
C                                                                       6156.   
C     CALL SETFOR  TO READ IN AND/OR INITIALIZE DATA AND/OR RESET PARAMS6157.   
C                                                                       6158.   
C     IF(NFTFOR.GT.0)  FORCING DATA WILL BE READ IN FROM DISKUNIT=NFTFOR6159.   
C     IF(NFTFOR.EQ.0)  NO DATA READ, SELECT CONSTITUENTS FOR EXT FORCING6160.   
C     IF(NFTFOR.LT.0)  NO DATA READ, RESET ONLY SOL CONST REFERENCE VALU6161.   
C-----------------------------------------------------------------------6162.   
C                                                                       6163.   
      DIMENSION YEAR(221),SCO2(221),SCH4(221),SN2O(221)                 6164.   
      DIMENSION SF11(221),SF12(221),UPPM(221)                           6165.   
      DIMENSION TAUS(12,221),TAUM(2652)                                 6166.   
      EQUIVALENCE (TAUS(1,1),TAUM(1))                                   6167.   
C                                                                       6168.   
      DIMENSION INDEX(9),INFOR(9)                                       6169.   
      EQUIVALENCE (INFOR(1),KVOL),(INFOR(2),KCO2),(INFOR(3),KXXX)       6170.   
      EQUIVALENCE (INFOR(4),KSUN),(INFOR(5),KYYY),(INFOR(6),KN2O)       6171.   
      EQUIVALENCE (INFOR(7),KCH4),(INFOR(8),KF11),(INFOR(9),KF12)       6172.   
C                                                                       6173.   
      DIMENSION DMO(12),JDY(12)                                         6174.   
      DATA DMO/31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31./         6175.   
      DATA JDY/31,59,90,120,151,181,212,243,278,304,334,365/            6176.   
C                                                                       6177.   
      IF(NFTFOR.LT.0) GO TO 150                                         6178.   
      IF(NFTFOR.LT.1) GO TO 110                                         6179.   
C                                                                       6180.   
      REWIND NFTFOR                                                     6181.   
      READ  (NFTFOR) NOUT,NEND,KFS,KCS,(YEAR(L),SCO2(L),SCH4(L),SN2O(L) 6182.   
     +              ,SF11(L),SF12(L),UPPM(L),(TAUS(K,L),K=1,12),L=1,221)6183.   
     +              ,IDATE                                              6184.   
      REWIND NFTFOR                                                     6184.5  
C                                                                       6185.   
      ID5(5)=IDATE+10*KFS+KCS                                           6186.   
C                                                                       6187.   
C-----------------------------------------------------------------------6188.   
C     REFERENCE YEAR IS (1958) WHERE FULGAS(K)=1 FOR CO2,N2O,CH4,F11,F126189.   
C     MEAN 1958 BACKGROUND CO2=315 N2O=.295 CH4=1.4 F11=8.E-6 F12=25.E-66190.   
C     GAS PPM IS LINEARLY INTERPOLATED (MEAN ANNUAL PPM OCCURS JDAY=183)6191.   
C                                                                       6192.   
C     BACKGROUND TAU STRATAER=0.012  (VOLCANIC CONTRIBUTION IS ADDITIVE)6193.   
C                                                                       6194.   
C           KFS=IDENTIFIER FOR F11,F12 ABUNDANCE SCENARIOS              6195.   
C           KCS=IDENTIFIER FOR CO2 ABUNDANCE SCENARIOS                  6196.   
C           ID5(5)=IDATE+10*KFS+KCS IS THE FORCING DATA SET IDENTIFIER  6197.   
C-----------------------------------------------------------------------6198.   
C                                                                       6199.   
      RRCO2=PPMV58(2)                                                   6200.   
      RCH4=PPMV58(7)                                                    6201.   
      RN2O=PPMV58(6)                                                    6202.   
C                           (F11,F12 EXTERNAL FORCING DATA ARE IN PPM)  6203.   
      RF11=PPMV58(8)*1000.                                              6204.   
      RF12=PPMV58(9)*1000.                                              6205.   
C                                                                       6206.   
      RVOL=AGOLDH(1,1)                                                  6207.   
C-----------------------------------------------------------------------6208.   
C                                                                       6209.   
C     SELECT CONSTITUENTS FOR WHICH EXTERNAL FORCING WILL BE IMPLEMENTED6210.   
C                                                                       6211.   
C     KFORCE  IS AN INTEGER UP TO NINE DIGITS LONG, SUCH THAT EACH DIGIT6212.   
C             IS AN ON/OFF SWITCH FOR IMPLEMENTING EXTERNAL FORCING FOR:6213.   
C                                                                       6214.   
C             (1)      (2)  (4)      (6)  (7)  (8)  (9)   CODED DIGITS  6215.   
C             VOL-AER, CO2, SOL-CON, N2O, CH4, F11, F12,  RESPECTIVELY. 6216.   
C                                 (THE DIGITS  (3) & (5)...ARE NOT USED)6217.   
C                                                                       6218.   
C             EXAMPLE:  1206789 SELECTS FORCING FOR ALL EXCEPT SOL CONST6219.   
C                       (ORDER OR REPETITION OF DIGITS IS NOT IMPORTANT)6220.   
C-----------------------------------------------------------------------6221.   
  110 KFOR=KFORCE                                                       6222.   
      KMAG=100000000                                                    6223.   
      DO 120 K=1,9                                                      6224.   
      KF=KFOR/KMAG                                                      6225.   
      INDEX(K)=KF                                                       6226.   
      KFOR=KFOR-KF*KMAG                                                 6227.   
  120 KMAG=KMAG/10                                                      6228.   
      DO 130 K=1,9                                                      6229.   
  130 INFOR(K)=0                                                        6230.   
      DO 140 K=1,9                                                      6231.   
      IF(INDEX(K).EQ.0) GO TO 140                                       6232.   
      INFOR(INDEX(K))=1                                                 6233.   
  140 CONTINUE                                                          6234.   
C                                                                       6235.   
C-----------------------------------------------------------------------6236.   
C     SELECT REFERENCE SOLAR CONSTANT (S0)  AS PASSED IN  COMMON/RADCOM/6237.   
C-----------------------------------------------------------------------6238.   
C                                                                       6239.   
  150 S00=S0                                                            6240.   
      RETURN                                                            6241.   
C                                                                       6242.   
C-----------------                                                      6243.   
      ENTRY GETFOR                                                      6244.   
C-----------------                                                      6245.   
C                                                                       6246.   
C-----------------------------------------------------------------------6247.   
C     EXTERNAL FORCING RETURNED FOR CONSTITUENTS PRESELECTED IN   SETFOR6248.   
C                                                                       6249.   
C            RADCOM  INPUT DATA:   JYEAR, JDAY                          6250.   
C                                                                       6251.   
C            RADCOM OUTPUT DATA:   FULGAS(K),K=2,6,7,8,9;  FGOLDH(1), S06252.   
C                                                                       6253.   
C-----------------------------------------------------------------------6254.   
C                                                                       6255.   
      JDM=JDAY                                                          6256.   
      DO 210 JMONTH=1,12                                                6257.   
      IF(JDAY.GT.JDY(JMONTH)) GO TO 210                                 6258.   
      GO TO 220                                                         6259.   
  210 JDM=JDAY-JDY(JMONTH)                                              6260.   
      JMONTH=12                                                         6261.   
  220 MO=JMONTH+(JYEAR-1880)*12                                         6262.   
      IF(MO.LT.   1) MO=1                                               6263.   
      IF(MO.GT.2651) MO=2651                                            6264.   
C                                                                       6265.   
      FRACYR=(JDAY-183)/365.                                            6266.   
      FRACMO=JDM/DMO(JMONTH)                                            6267.   
C                                                                       6268.   
      NY=JYEAR-1880+1                                                   6269.   
      IF(JDAY.LT.183) NY=NY-1                                           6270.   
      IF(JDAY.LT.183) FRACYR=FRACYR+0.5                                 6271.   
      IF(NY.LT.  1) NY=1                                                6272.   
      IF(NY.GT.220) NY=220                                              6273.   
      FCO2=SCO2(NY)+(SCO2(NY+1)-SCO2(NY))*FRACYR                        6274.   
      FCH4=SCH4(NY)+(SCH4(NY+1)-SCH4(NY))*FRACYR                        6275.   
      FN2O=SN2O(NY)+(SN2O(NY+1)-SN2O(NY))*FRACYR                        6276.   
      FF11=SF11(NY)+(SF11(NY+1)-SF11(NY))*FRACYR                        6277.   
      FF12=SF12(NY)+(SF12(NY+1)-SF12(NY))*FRACYR                        6278.   
      FSUN=UPPM(NY)+(UPPM(NY+1)-UPPM(NY))*FRACYR                        6279.   
      FVOL=TAUM(MO)+(TAUM(MO+1)-TAUM(MO))*FRACMO                        6280.   
C                                                                       6281.   
C-----------------------------------------------------------------------6282.   
C                                                    OUTPUT FORCING DATA6283.   
C-----------------------------------------------------------------------6284.   
C                                                                       6285.   
      IF(KCO2.GT.0) FULGAS(2)=FCO2/RRCO2                                6286.   
      IF(KN2O.GT.0) FULGAS(6)=FN2O/RN2O                                 6287.   
      IF(KCH4.GT.0) FULGAS(7)=FCH4/RCH4                                 6288.   
      IF(KF11.GT.0) FULGAS(8)=FF11/RF11                                 6289.   
      IF(KF12.GT.0) FULGAS(9)=FF12/RF12                                 6290.   
      IF(KVOL.GT.0) FGOLDH(1)=(RVOL+FVOL)/RVOL                          6291.   
      IF(KSUN.GT.0) S0=S00+S00*0.03*(FSUN-0.2)                          6292.   
C                                                                       6293.   
      RETURN                                                            6294.   
      END                                                               6295.   
      SUBROUTINE HGAER1(XMU,TAU,G,GG)                                   6301.   
C                                                                       6302.   
      DIMENSION C05T00(51),C06T00(51),C07T00(51),C08T00(51),C09T00(51)  6303.   
      DIMENSION C05T01(51),C06T01(51),C07T01(51),C08T01(51),C09T01(51)  6304.   
      DIMENSION C05T02(51),C06T02(51),C07T02(51),C08T02(51),C09T02(51)  6305.   
      DIMENSION C05T03(51),C06T03(51),C07T03(51),C08T03(51),C09T03(51)  6306.   
      DIMENSION C05T04(51),C06T04(51),C07T04(51),C08T04(51),C09T04(51)  6307.   
      DIMENSION C05T05(51),C06T05(51),C07T05(51),C08T05(51),C09T05(51)  6308.   
      DIMENSION C05T06(51),C06T06(51),C07T06(51),C08T06(51),C09T06(51)  6309.   
      DIMENSION C05T07(51),C06T07(51),C07T07(51),C08T07(51),C09T07(51)  6310.   
      DIMENSION C05T08(51),C06T08(51),C07T08(51),C08T08(51),C09T08(51)  6311.   
      DIMENSION C05T09(51),C06T09(51),C07T09(51),C08T09(51),C09T09(51)  6312.   
      DIMENSION C05T10(51),C06T10(51),C07T10(51),C08T10(51),C09T10(51)  6313.   
C                                                                       6314.   
      DIMENSION C05TAU(51,11),C06TAU(51,11),C07TAU(51,11),C08TAU(51,11) 6315.   
      DIMENSION C09TAU(51,11)                                           6316.   
C                                                                       6317.   
      DIMENSION GTAU(51,11,5)                                           6318.   
C                                                                       6319.   
      EQUIVALENCE (C05TAU(1, 1),C05T00(1)),(C05TAU(1, 2),C05T01(1))     6320.   
      EQUIVALENCE (C05TAU(1, 3),C05T02(1)),(C05TAU(1, 4),C05T03(1))     6321.   
      EQUIVALENCE (C05TAU(1, 5),C05T04(1)),(C05TAU(1, 6),C05T05(1))     6322.   
      EQUIVALENCE (C05TAU(1, 7),C05T06(1)),(C05TAU(1, 8),C05T07(1))     6323.   
      EQUIVALENCE (C05TAU(1, 9),C05T08(1)),(C05TAU(1,10),C05T09(1))     6324.   
      EQUIVALENCE (C05TAU(1,11),C05T10(1))                              6325.   
C                                                                       6326.   
      EQUIVALENCE (C06TAU(1, 1),C06T00(1)),(C06TAU(1, 2),C06T01(1))     6327.   
      EQUIVALENCE (C06TAU(1, 3),C06T02(1)),(C06TAU(1, 4),C06T03(1))     6328.   
      EQUIVALENCE (C06TAU(1, 5),C06T04(1)),(C06TAU(1, 6),C06T05(1))     6329.   
      EQUIVALENCE (C06TAU(1, 7),C06T06(1)),(C06TAU(1, 8),C06T07(1))     6330.   
      EQUIVALENCE (C06TAU(1, 9),C06T08(1)),(C06TAU(1,10),C06T09(1))     6331.   
      EQUIVALENCE (C06TAU(1,11),C06T10(1))                              6332.   
C                                                                       6333.   
      EQUIVALENCE (C07TAU(1, 1),C07T00(1)),(C07TAU(1, 2),C07T01(1))     6334.   
      EQUIVALENCE (C07TAU(1, 3),C07T02(1)),(C07TAU(1, 4),C07T03(1))     6335.   
      EQUIVALENCE (C07TAU(1, 5),C07T04(1)),(C07TAU(1, 6),C07T05(1))     6336.   
      EQUIVALENCE (C07TAU(1, 7),C07T06(1)),(C07TAU(1, 8),C07T07(1))     6337.   
      EQUIVALENCE (C07TAU(1, 9),C07T08(1)),(C07TAU(1,10),C07T09(1))     6338.   
      EQUIVALENCE (C07TAU(1,11),C07T10(1))                              6339.   
C                                                                       6340.   
      EQUIVALENCE (C08TAU(1, 1),C08T00(1)),(C08TAU(1, 2),C08T01(1))     6341.   
      EQUIVALENCE (C08TAU(1, 3),C08T02(1)),(C08TAU(1, 4),C08T03(1))     6342.   
      EQUIVALENCE (C08TAU(1, 5),C08T04(1)),(C08TAU(1, 6),C08T05(1))     6343.   
      EQUIVALENCE (C08TAU(1, 7),C08T06(1)),(C08TAU(1, 8),C08T07(1))     6344.   
      EQUIVALENCE (C08TAU(1, 9),C08T08(1)),(C08TAU(1,10),C08T09(1))     6345.   
      EQUIVALENCE (C08TAU(1,11),C08T10(1))                              6346.   
C                                                                       6347.   
      EQUIVALENCE (C09TAU(1, 1),C09T00(1)),(C09TAU(1, 2),C09T01(1))     6348.   
      EQUIVALENCE (C09TAU(1, 3),C09T02(1)),(C09TAU(1, 4),C09T03(1))     6349.   
      EQUIVALENCE (C09TAU(1, 5),C09T04(1)),(C09TAU(1, 6),C09T05(1))     6350.   
      EQUIVALENCE (C09TAU(1, 7),C09T06(1)),(C09TAU(1, 8),C09T07(1))     6351.   
      EQUIVALENCE (C09TAU(1, 9),C09T08(1)),(C09TAU(1,10),C09T09(1))     6352.   
      EQUIVALENCE (C09TAU(1,11),C09T10(1))                              6353.   
C                                                                       6354.   
      EQUIVALENCE (C05TAU(1,1),GTAU(1,1,1)),(C06TAU(1,1),GTAU(1,1,2))   6355.   
      EQUIVALENCE (C07TAU(1,1),GTAU(1,1,3)),(C08TAU(1,1),GTAU(1,1,4))   6356.   
      EQUIVALENCE (C09TAU(1,1),GTAU(1,1,5))                             6357.   
C                                                                       6358.   
C                                                                       6359.   
      DATA C05T00/0.0,                                                  6360.   
     1    .0179,.0379,.0574,.0767,.0958,.1147,.1334,.1520,.1703,.1884,  6361.   
     2    .2062,.2238,.2410,.2580,.2747,.2910,.3070,.3226,.3380,.3530,  6362.   
     3    .3675,.3819,.3958,.4094,.4227,.4355,.4481,.4603,.4722,.4838,  6363.   
     4    .4950,.5059,.5166,.5269,.5370,.5468,.5563,.5655,.5745,.5832,  6364.   
     5    .5917,.5999,.6079,.6157,.6233,.6306,.6378,.6445,.6513,.6578/  6365.   
C                                                                       6366.   
      DATA C05T01/0.0,                                                  6367.   
     1    .0000,.0226,.0463,.0679,.0885,.1084,.1278,.1469,.1655,.1838,  6368.   
     2    .2018,.2194,.2367,.2537,.2704,.2866,.3026,.3182,.3335,.3484,  6369.   
     3    .3630,.3773,.3911,.4047,.4180,.4308,.4433,.4556,.4675,.4791,  6370.   
     4    .4904,.5014,.5121,.5224,.5326,.5424,.5520,.5613,.5703,.5792,  6371.   
     5    .5877,.5961,.6041,.6120,.6197,.6271,.6344,.6414,.6483,.6550/  6372.   
C                                                                       6373.   
      DATA C05T02/0.0,                                                  6374.   
     1    .0000,.0207,.0434,.0649,.0856,.1057,.1252,.1444,.1632,.1816,  6375.   
     2    .1996,.2173,.2346,.2516,.2683,.2845,.3005,.3161,.3313,.3463,  6376.   
     3    .3608,.3750,.3889,.4024,.4156,.4284,.4410,.4532,.4651,.4767,  6377.   
     4    .4880,.4990,.5097,.5201,.5303,.5401,.5497,.5591,.5682,.5771,  6378.   
     5    .5857,.5941,.6022,.6102,.6179,.6254,.6327,.6398,.6467,.6535/  6379.   
C                                                                       6380.   
      DATA C05T03/0.0,                                                  6381.   
     1    .0095,.0317,.0517,.0712,.0904,.1095,.1283,.1469,.1651,.1832,  6382.   
     2    .2009,.2184,.2355,.2523,.2688,.2849,.3008,.3162,.3313,.3461,  6383.   
     3    .3605,.3747,.3885,.4019,.4151,.4278,.4403,.4525,.4643,.4759,  6384.   
     4    .4872,.4981,.5089,.5192,.5294,.5392,.5488,.5582,.5673,.5762,  6385.   
     5    .5848,.5932,.6013,.6093,.6170,.6246,.6319,.6391,.6460,.6528/  6386.   
C                                                                       6387.   
      DATA C05T04/0.0,                                                  6388.   
     1    .0260,.0472,.0656,.0833,.1008,.1183,.1359,.1534,.1709,.1882,  6389.   
     2    .2053,.2223,.2389,.2554,.2715,.2873,.3029,.3181,.3330,.3476,  6390.   
     3    .3619,.3759,.3895,.4028,.4158,.4284,.4408,.4529,.4647,.4762,  6391.   
     4    .4873,.4982,.5089,.5192,.5293,.5391,.5487,.5580,.5671,.5759,  6392.   
     5    .5845,.5929,.6010,.6090,.6167,.6243,.6316,.6388,.6457,.6525/  6393.   
C                                                                       6394.   
      DATA C05T05/0.0,                                                  6395.   
     1    .0428,.0635,.0812,.0978,.1140,.1302,.1465,.1629,.1793,.1958,  6396.   
     2    .2121,.2284,.2444,.2603,.2760,.2914,.3066,.3214,.3360,.3504,  6397.   
     3    .3643,.3781,.3915,.4046,.4175,.4299,.4422,.4541,.4657,.4771,  6398.   
     4    .4882,.4990,.5095,.5197,.5298,.5395,.5490,.5583,.5673,.5761,  6399.   
     5    .5846,.5930,.6011,.6090,.6167,.6243,.6316,.6387,.6457,.6524/  6400.   
C                                                                       6401.   
      DATA C05T06/0.0,                                                  6402.   
     1    .0590,.0796,.0969,.1129,.1283,.1435,.1588,.1741,.1896,.2051,  6403.   
     2    .2206,.2360,.2514,.2667,.2818,.2967,.3114,.3258,.3401,.3541,  6404.   
     3    .3677,.3812,.3943,.4072,.4198,.4321,.4441,.4559,.4673,.4786,  6405.   
     4    .4895,.5002,.5106,.5207,.5306,.5403,.5497,.5589,.5678,.5766,  6406.   
     5    .5850,.5934,.6014,.6093,.6170,.6244,.6317,.6388,.6458,.6525/  6407.   
C                                                                       6408.   
      DATA C05T07/0.0,                                                  6409.   
     1    .0742,.0948,.1120,.1277,.1427,.1572,.1716,.1861,.2007,.2153,  6410.   
     2    .2300,.2447,.2594,.2740,.2885,.3028,.3171,.3310,.3448,.3584,  6411.   
     3    .3717,.3849,.3977,.4103,.4227,.4347,.4465,.4581,.4693,.4804,  6412.   
     4    .4912,.5017,.5120,.5220,.5318,.5413,.5506,.5597,.5686,.5772,  6413.   
     5    .5856,.5939,.6019,.6097,.6173,.6247,.6320,.6390,.6459,.6526/  6414.   
C                                                                       6415.   
      DATA C05T08/0.0,                                                  6416.   
     1    .0885,.1090,.1263,.1418,.1565,.1705,.1844,.1982,.2121,.2260,  6417.   
     2    .2400,.2540,.2680,.2819,.2958,.3096,.3233,.3368,.3502,.3633,  6418.   
     3    .3763,.3890,.4015,.4138,.4259,.4377,.4493,.4606,.4717,.4825,  6419.   
     4    .4931,.5035,.5136,.5235,.5331,.5425,.5517,.5607,.5695,.5780,  6420.   
     5    .5864,.5945,.6024,.6102,.6177,.6251,.6323,.6393,.6461,.6528/  6421.   
C                                                                       6422.   
      DATA C05T09/0.0,                                                  6423.   
     1    .1017,.1223,.1395,.1550,.1695,.1833,.1968,.2101,.2234,.2367,  6424.   
     2    .2501,.2634,.2768,.2902,.3035,.3167,.3299,.3429,.3558,.3686,  6425.   
     3    .3811,.3935,.4057,.4176,.4295,.4409,.4523,.4634,.4742,.4849,  6426.   
     4    .4952,.5054,.5154,.5251,.5346,.5439,.5530,.5618,.5705,.5789,  6427.   
     5    .5871,.5952,.6031,.6107,.6182,.6255,.6326,.6396,.6464,.6530/  6428.   
C                                                                       6429.   
      DATA C05T10/0.0,                                                  6430.   
     1    .1139,.1346,.1518,.1673,.1817,.1953,.2086,.2216,.2344,.2472,  6431.   
     2    .2600,.2728,.2857,.2985,.3113,.3240,.3367,.3492,.3617,.3740,  6432.   
     3    .3862,.3982,.4100,.4217,.4332,.4444,.4554,.4663,.4769,.4873,  6433.   
     4    .4975,.5075,.5173,.5268,.5362,.5453,.5543,.5630,.5715,.5799,  6434.   
     5    .5880,.5960,.6037,.6113,.6187,.6259,.6330,.6399,.6466,.6532/  6435.   
C                                                                       6436.   
      DATA C06T00/0.0,                                                  6437.   
     1    .0250,.0525,.0792,.1056,.1316,.1572,.1823,.2070,.2311,.2547,  6438.   
     2    .2776,.3000,.3217,.3427,.3631,.3827,.4019,.4201,.4378,.4550,  6439.   
     3    .4713,.4872,.5024,.5170,.5312,.5446,.5576,.5701,.5820,.5936,  6440.   
     4    .6047,.6153,.6257,.6354,.6450,.6541,.6628,.6713,.6794,.6873,  6441.   
     5    .6948,.7021,.7091,.7159,.7224,.7287,.7348,.7407,.7462,.7516/  6442.   
C                                                                       6443.   
      DATA C06T01/0.0,                                                  6444.   
     1    .0000,.0339,.0652,.0941,.1216,.1480,.1737,.1987,.2229,.2466,  6445.   
     2    .2694,.2918,.3134,.3344,.3548,.3744,.3935,.4118,.4295,.4467,  6446.   
     3    .4632,.4792,.4945,.5092,.5236,.5372,.5504,.5631,.5753,.5871,  6447.   
     4    .5984,.6093,.6198,.6299,.6396,.6490,.6580,.6667,.6751,.6832,  6448.   
     5    .6909,.6984,.7056,.7126,.7194,.7259,.7322,.7382,.7441,.7498/  6449.   
C                                                                       6450.   
      DATA C06T02/0.0,                                                  6451.   
     1    .0000,.0307,.0608,.0893,.1168,.1433,.1690,.1941,.2183,.2420,  6452.   
     2    .2648,.2871,.3087,.3296,.3500,.3696,.3887,.4070,.4247,.4420,  6453.   
     3    .4584,.4745,.4898,.5047,.5191,.5328,.5461,.5590,.5713,.5832,  6454.   
     4    .5947,.6057,.6164,.6266,.6365,.6460,.6552,.6641,.6726,.6808,  6455.   
     5    .6887,.6964,.7038,.7110,.7178,.7245,.7309,.7371,.7431,.7489/  6456.   
C                                                                       6457.   
      DATA C06T03/0.0,                                                  6458.   
     1    .0130,.0424,.0692,.0953,.1210,.1462,.1709,.1952,.2188,.2420,  6459.   
     2    .2645,.2865,.3078,.3285,.3486,.3680,.3870,.4051,.4228,.4399,  6460.   
     3    .4563,.4723,.4877,.5025,.5169,.5306,.5440,.5569,.5692,.5812,  6461.   
     4    .5927,.6038,.6146,.6248,.6348,.6444,.6537,.6626,.6712,.6796,  6462.   
     5    .6876,.6954,.7028,.7101,.7170,.7238,.7303,.7366,.7427,.7486/  6463.   
C                                                                       6464.   
      DATA C06T04/0.0,                                                  6465.   
     1    .0314,.0594,.0842,.1080,.1315,.1549,.1781,.2012,.2238,.2461,  6466.   
     2    .2678,.2892,.3099,.3302,.3499,.3690,.3876,.4055,.4230,.4399,  6467.   
     3    .4561,.4720,.4872,.5019,.5163,.5299,.5432,.5561,.5684,.5804,  6468.   
     4    .5918,.6029,.6137,.6240,.6340,.6436,.6529,.6619,.6705,.6790,  6469.   
     5    .6870,.6948,.7023,.7096,.7167,.7235,.7300,.7364,.7425,.7485/  6470.   
C                                                                       6471.   
      DATA C06T05/0.0,                                                  6472.   
     1    .0503,.0777,.1014,.1237,.1456,.1673,.1889,.2105,.2319,.2531,  6473.   
     2    .2739,.2944,.3145,.3341,.3533,.3718,.3901,.4076,.4247,.4413,  6474.   
     3    .4573,.4730,.4880,.5025,.5167,.5302,.5434,.5562,.5684,.5803,  6475.   
     4    .5917,.6028,.6135,.6238,.6338,.6434,.6527,.6617,.6703,.6787,  6476.   
     5    .6868,.6946,.7021,.7095,.7165,.7233,.7299,.7363,.7425,.7485/  6477.   
C                                                                       6478.   
      DATA C06T06/0.0,                                                  6479.   
     1    .0686,.0956,.1188,.1403,.1611,.1814,.2017,.2220,.2421,.2622,  6480.   
     2    .2820,.3016,.3208,.3397,.3582,.3762,.3939,.4110,.4276,.4439,  6481.   
     3    .4596,.4749,.4897,.5040,.5180,.5313,.5443,.5569,.5690,.5808,  6482.   
     4    .5921,.6031,.6138,.6240,.6339,.6435,.6527,.6617,.6703,.6787,  6483.   
     5    .6868,.6946,.7021,.7094,.7165,.7233,.7300,.7364,.7425,.7485/  6484.   
C                                                                       6485.   
      DATA C06T07/0.0,                                                  6486.   
     1    .0859,.1128,.1357,.1567,.1767,.1961,.2154,.2345,.2535,.2725,  6487.   
     2    .2913,.3099,.3283,.3464,.3642,.3816,.3987,.4153,.4315,.4473,  6488.   
     3    .4626,.4776,.4920,.5061,.5198,.5329,.5457,.5582,.5701,.5818,  6489.   
     4    .5930,.6038,.6144,.6245,.6344,.6439,.6530,.6620,.6705,.6789,  6490.   
     5    .6869,.6947,.7022,.7095,.7166,.7234,.7300,.7364,.7426,.7486/  6491.   
C                                                                       6492.   
      DATA C06T08/0.0,                                                  6493.   
     1    .1022,.1290,.1517,.1723,.1919,.2107,.2291,.2473,.2654,.2834,  6494.   
     2    .3013,.3191,.3366,.3539,.3710,.3877,.4042,.4202,.4360,.4513,  6495.   
     3    .4662,.4808,.4950,.5087,.5221,.5350,.5476,.5598,.5715,.5830,  6496.   
     4    .5941,.6048,.6152,.6252,.6350,.6444,.6535,.6624,.6709,.6792,  6497.   
     5    .6872,.6949,.7024,.7097,.7167,.7235,.7301,.7365,.7427,.7486/  6498.   
C                                                                       6499.   
      DATA C06T09/0.0,                                                  6500.   
     1    .1173,.1440,.1666,.1871,.2063,.2246,.2425,.2600,.2773,.2945,  6501.   
     2    .3116,.3285,.3453,.3619,.3783,.3943,.4102,.4257,.4409,.4558,  6502.   
     3    .4703,.4845,.4982,.5116,.5248,.5374,.5497,.5617,.5732,.5845,  6503.   
     4    .5954,.6060,.6163,.6262,.6358,.6451,.6541,.6629,.6713,.6796,  6504.   
     5    .6875,.6952,.7026,.7099,.7168,.7236,.7302,.7365,.7427,.7487/  6505.   
C                                                                       6506.   
      DATA C06T10/0.0,                                                  6507.   
     1    .1314,.1581,.1806,.2009,.2199,.2379,.2553,.2722,.2889,.3055,  6508.   
     2    .3218,.3381,.3541,.3700,.3858,.4012,.4165,.4315,.4462,.4606,  6509.   
     3    .4746,.4884,.5018,.5148,.5277,.5400,.5520,.5638,.5751,.5862,  6510.   
     4    .5969,.6073,.6175,.6272,.6367,.6459,.6548,.6635,.6719,.6800,  6511.   
     5    .6879,.6955,.7029,.7101,.7170,.7237,.7303,.7366,.7427,.7487/  6512.   
C                                                                       6513.   
      DATA C07T00/0.0,                                                  6514.   
     1    .0360,.0751,.1129,.1498,.1858,.2209,.2546,.2873,.3183,.3484,  6515.   
     2    .3767,.4040,.4296,.4540,.4773,.4990,.5199,.5392,.5577,.5753,  6516.   
     3    .5916,.6073,.6220,.6358,.6492,.6615,.6733,.6845,.6950,.7051,  6517.   
     4    .7147,.7237,.7324,.7406,.7484,.7559,.7630,.7698,.7762,.7824,  6518.   
     5    .7883,.7940,.7994,.8046,.8096,.8144,.8190,.8234,.8276,.8317/  6519.   
C                                                                       6520.   
      DATA C07T01/0.0,                                                  6521.   
     1    .0000,.0500,.0929,.1323,.1696,.2052,.2391,.2719,.3029,.3329,  6522.   
     2    .3612,.3886,.4144,.4390,.4625,.4845,.5058,.5256,.5445,.5626,  6523.   
     3    .5795,.5957,.6109,.6253,.6392,.6521,.6644,.6762,.6872,.6979,  6524.   
     4    .7079,.7174,.7266,.7351,.7434,.7513,.7587,.7659,.7727,.7793,  6525.   
     5    .7855,.7915,.7971,.8026,.8079,.8129,.8177,.8223,.8268,.8310/  6526.   
C                                                                       6527.   
      DATA C07T02/0.0,                                                  6528.   
     1    .0000,.0433,.0845,.1233,.1604,.1958,.2296,.2623,.2932,.3232,  6529.   
     2    .3515,.3788,.4047,.4294,.4530,.4753,.4967,.5168,.5360,.5544,  6530.   
     3    .5715,.5881,.6037,.6184,.6327,.6459,.6586,.6707,.6821,.6931,  6531.   
     4    .7034,.7133,.7228,.7316,.7402,.7484,.7561,.7636,.7706,.7774,  6532.   
     5    .7839,.7901,.7960,.8017,.8071,.8123,.8173,.8221,.8267,.8311/  6533.   
C                                                                       6534.   
      DATA C07T03/0.0,                                                  6535.   
     1    .0139,.0544,.0915,.1272,.1620,.1958,.2284,.2601,.2903,.3197,  6536.   
     2    .3475,.3745,.4001,.4246,.4481,.4703,.4918,.5119,.5311,.5496,  6537.   
     3    .5669,.5836,.5993,.6142,.6287,.6420,.6550,.6673,.6789,.6901,  6538.   
     4    .7006,.7107,.7204,.7294,.7382,.7465,.7545,.7621,.7693,.7763,  6539.   
     5    .7829,.7893,.7953,.8012,.8067,.8121,.8172,.8221,.8269,.8314/  6540.   
C                                                                       6541.   
      DATA C07T04/0.0,                                                  6542.   
     1    .0339,.0723,.1065,.1393,.1714,.2028,.2336,.2637,.2927,.3210,  6543.   
     2    .3480,.3743,.3993,.4234,.4465,.4684,.4897,.5096,.5288,.5471,  6544.   
     3    .5644,.5811,.5968,.6118,.6263,.6398,.6528,.6652,.6769,.6882,  6545.   
     4    .6988,.7090,.7188,.7280,.7369,.7454,.7534,.7612,.7685,.7756,  6546.   
     5    .7823,.7888,.7950,.8009,.8066,.8120,.8173,.8223,.8271,.8317/  6547.   
C                                                                       6548.   
      DATA C07T05/0.0,                                                  6549.   
     1    .0546,.0920,.1246,.1553,.1852,.2144,.2432,.2715,.2990,.3260,  6550.   
     2    .3519,.3772,.4015,.4249,.4474,.4689,.4897,.5093,.5283,.5464,  6551.   
     3    .5635,.5801,.5957,.6106,.6251,.6386,.6516,.6640,.6757,.6871,  6552.   
     4    .6978,.7080,.7179,.7272,.7361,.7447,.7528,.7606,.7680,.7752,  6553.   
     5    .7820,.7886,.7948,.8008,.8065,.8121,.8174,.8224,.8273,.8320/  6554.   
C                                                                       6555.   
      DATA C07T06/0.0,                                                  6556.   
     1    .0749,.1117,.1434,.1728,.2010,.2284,.2554,.2820,.3079,.3335,  6557.   
     2    .3582,.3825,.4058,.4284,.4502,.4711,.4914,.5106,.5292,.5470,  6558.   
     3    .5639,.5802,.5957,.6105,.6248,.6382,.6511,.6635,.6752,.6865,  6559.   
     4    .6972,.7075,.7174,.7267,.7357,.7442,.7524,.7603,.7677,.7750,  6560.   
     5    .7818,.7884,.7947,.8008,.8065,.8121,.8174,.8226,.8275,.8322/  6561.   
C                                                                       6562.   
      DATA C07T07/0.0,                                                  6563.   
     1    .0943,.1306,.1617,.1902,.2173,.2434,.2689,.2940,.3185,.3427,  6564.   
     2    .3662,.3893,.4117,.4334,.4545,.4747,.4944,.5131,.5312,.5486,  6565.   
     3    .5651,.5812,.5964,.6110,.6252,.6384,.6512,.6635,.6752,.6864,  6566.   
     4    .6971,.7073,.7172,.7265,.7355,.7440,.7522,.7601,.7676,.7748,  6567.   
     5    .7817,.7883,.7946,.8007,.8065,.8121,.8175,.8227,.8276,.8324/  6568.   
C                                                                       6569.   
      DATA C07T08/0.0,                                                  6570.   
     1    .1125,.1486,.1793,.2071,.2334,.2585,.2828,.3066,.3299,.3529,  6571.   
     2    .3753,.3973,.4186,.4395,.4597,.4792,.4982,.5164,.5340,.5510,  6572.   
     3    .5672,.5829,.5978,.6122,.6261,.6392,.6518,.6640,.6755,.6867,  6573.   
     4    .6973,.7074,.7172,.7265,.7354,.7440,.7522,.7600,.7675,.7748,  6574.   
     5    .7816,.7883,.7946,.8007,.8066,.8122,.8176,.8228,.8277,.8325/  6575.   
C                                                                       6576.   
      DATA C07T09/0.0,                                                  6577.   
     1    .1296,.1655,.1958,.2232,.2489,.2732,.2966,.3194,.3416,.3635,  6578.   
     2    .3848,.4058,.4262,.4462,.4656,.4844,.5028,.5203,.5374,.5539,  6579.   
     3    .5697,.5850,.5997,.6137,.6274,.6403,.6527,.6647,.6761,.6872,  6580.   
     4    .6977,.7077,.7175,.7267,.7356,.7441,.7522,.7601,.7675,.7748,  6581.   
     5    .7816,.7883,.7946,.8007,.8066,.8122,.8176,.8228,.8278,.8326/  6582.   
C                                                                       6583.   
      DATA C07T10/0.0,                                                  6584.   
     1    .1456,.1813,.2114,.2384,.2635,.2872,.3099,.3319,.3532,.3742,  6585.   
     2    .3946,.4147,.4342,.4533,.4720,.4901,.5078,.5248,.5413,.5573,  6586.   
     3    .5727,.5876,.6019,.6156,.6290,.6417,.6539,.6657,.6770,.6879,  6587.   
     4    .6982,.7082,.7178,.7270,.7358,.7442,.7523,.7602,.7676,.7748,  6588.   
     5    .7816,.7883,.7946,.8007,.8066,.8122,.8176,.8228,.8278,.8326/  6589.   
C                                                                       6590.   
      DATA C08T00/0.0,                                                  6591.   
     1    .0568,.1172,.1747,.2295,.2813,.3300,.3748,.4169,.4547,.4903,  6592.   
     2    .5220,.5517,.5784,.6030,.6257,.6460,.6652,.6825,.6985,.7134,  6593.   
     3    .7269,.7396,.7513,.7621,.7723,.7816,.7904,.7987,.8064,.8137,  6594.   
     4    .8204,.8268,.8329,.8385,.8439,.8490,.8538,.8584,.8627,.8668,  6595.   
     5    .8707,.8744,.8780,.8814,.8846,.8877,.8906,.8934,.8961,.8987/  6596.   
C                                                                       6597.   
      DATA C08T01/0.0,                                                  6598.   
     1    .0045,.0786,.1413,.1980,.2505,.2994,.3445,.3870,.4255,.4620,  6599.   
     2    .4948,.5257,.5538,.5798,.6039,.6258,.6464,.6650,.6823,.6985,  6600.   
     3    .7132,.7270,.7398,.7516,.7629,.7730,.7826,.7917,.8000,.8080,  6601.   
     4    .8153,.8223,.8289,.8350,.8408,.8463,.8514,.8564,.8610,.8654,  6602.   
     5    .8696,.8736,.8773,.8809,.8843,.8876,.8907,.8937,.8965,.8992/  6603.   
C                                                                       6604.   
      DATA C08T02/0.0,                                                  6605.   
     1    .0000,.0639,.1239,.1794,.2314,.2799,.3249,.3675,.4063,.4431,  6606.   
     2    .4766,.5081,.5370,.5637,.5888,.6115,.6330,.6525,.6707,.6878,  6607.   
     3    .7032,.7179,.7314,.7440,.7559,.7667,.7769,.7865,.7954,.8038,  6608.   
     4    .8117,.8190,.8260,.8325,.8387,.8445,.8499,.8551,.8600,.8647,  6609.   
     5    .8690,.8733,.8772,.8810,.8845,.8880,.8912,.8943,.8973,.9001/  6610.   
C                                                                       6611.   
      DATA C08T03/0.0,                                                  6612.   
     1    .0129,.0725,.1266,.1778,.2266,.2730,.3165,.3580,.3962,.4326,  6613.   
     2    .4659,.4975,.5265,.5536,.5790,.6021,.6241,.6441,.6628,.6804,  6614.   
     3    .6964,.7116,.7256,.7386,.7510,.7622,.7728,.7828,.7921,.8009,  6615.   
     4    .8090,.8167,.8240,.8307,.8372,.8432,.8489,.8543,.8594,.8642,  6616.   
     5    .8688,.8731,.8772,.8811,.8848,.8884,.8917,.8949,.8980,.9009/  6617.   
C                                                                       6618.   
      DATA C08T04/0.0,                                                  6619.   
     1    .0338,.0901,.1399,.1870,.2320,.2754,.3165,.3561,.3930,.4283,  6620.   
     2    .4609,.4920,.5207,.5477,.5730,.5962,.6184,.6385,.6575,.6753,  6621.   
     3    .6916,.7071,.7214,.7347,.7474,.7589,.7698,.7801,.7896,.7987,  6622.   
     4    .8071,.8150,.8225,.8294,.8361,.8423,.8481,.8537,.8589,.8639,  6623.   
     5    .8686,.8731,.8773,.8813,.8851,.8887,.8922,.8955,.8986,.9016/  6624.   
C                                                                       6625.   
      DATA C08T05/0.0,                                                  6626.   
     1    .0561,.1105,.1578,.2017,.2435,.2838,.3224,.3597,.3948,.4287,  6627.   
     2    .4602,.4904,.5185,.5450,.5699,.5930,.6150,.6351,.6541,.6720,  6628.   
     3    .6884,.7040,.7185,.7319,.7448,.7565,.7676,.7781,.7877,.7970,  6629.   
     4    .8056,.8136,.8213,.8284,.8352,.8416,.8476,.8533,.8586,.8637,  6630.   
     5    .8685,.8731,.8774,.8815,.8854,.8891,.8926,.8960,.8991,.9022/  6631.   
C                                                                       6632.   
      DATA C08T06/0.0,                                                  6633.   
     1    .0782,.1314,.1770,.2187,.2581,.2958,.3319,.3670,.4002,.4324,  6634.   
     2    .4626,.4917,.5189,.5447,.5691,.5918,.6134,.6334,.6522,.6700,  6635.   
     3    .6864,.7020,.7165,.7300,.7430,.7548,.7660,.7766,.7864,.7957,  6636.   
     4    .8044,.8126,.8204,.8276,.8345,.8410,.8471,.8529,.8583,.8635,  6637.   
     5    .8684,.8731,.8774,.8816,.8856,.8893,.8929,.8963,.8996,.9027/  6638.   
C                                                                       6639.   
      DATA C08T07/0.0,                                                  6640.   
     1    .0994,.1518,.1962,.2363,.2739,.3095,.3436,.3765,.4080,.4385,  6641.   
     2    .4673,.4951,.5213,.5463,.5700,.5921,.6134,.6329,.6515,.6691,  6642.   
     3    .6854,.7009,.7154,.7289,.7418,.7536,.7649,.7755,.7854,.7948,  6643.   
     4    .8036,.8118,.8197,.8270,.8340,.8405,.8467,.8526,.8581,.8634,  6644.   
     5    .8683,.8731,.8775,.8817,.8857,.8896,.8932,.8967,.8999,.9031/  6645.   
C                                                                       6646.   
      DATA C08T08/0.0,                                                  6647.   
     1    .1197,.1714,.2148,.2538,.2899,.3238,.3562,.3874,.4172,.4461,  6648.   
     2    .4735,.5001,.5253,.5493,.5722,.5937,.6144,.6335,.6518,.6691,  6649.   
     3    .6852,.7005,.7148,.7283,.7412,.7529,.7642,.7748,.7847,.7942,  6650.   
     4    .8030,.8113,.8192,.8265,.8336,.8402,.8464,.8524,.8579,.8632,  6651.   
     5    .8682,.8730,.8775,.8818,.8858,.8897,.8934,.8969,.9002,.9034/  6652.   
C                                                                       6653.   
      DATA C08T09/0.0,                                                  6654.   
     1    .1387,.1899,.2326,.2705,.3055,.3382,.3691,.3988,.4271,.4546,  6655.   
     2    .4808,.5061,.5302,.5533,.5754,.5962,.6163,.6350,.6528,.6698,  6656.   
     3    .6855,.7007,.7148,.7281,.7409,.7526,.7638,.7744,.7843,.7937,  6657.   
     4    .8025,.8109,.8188,.8262,.8333,.8399,.8462,.8521,.8577,.8631,  6658.   
     5    .8681,.8730,.8775,.8818,.8859,.8898,.8935,.8971,.9004,.9036/  6659.   
C                                                                       6660.   
      DATA C08T10/0.0,                                                  6661.   
     1    .1567,.2073,.2493,.2865,.3206,.3522,.3819,.4104,.4374,.4636,  6662.   
     2    .4886,.5128,.5359,.5580,.5793,.5994,.6188,.6370,.6544,.6710,  6663.   
     3    .6864,.7013,.7152,.7283,.7410,.7526,.7637,.7742,.7840,.7935,  6664.   
     4    .8023,.8106,.8185,.8259,.8330,.8396,.8459,.8519,.8576,.8630,  6665.   
     5    .8680,.8729,.8774,.8818,.8859,.8899,.8936,.8972,.9006,.9038/  6666.   
C                                                                       6667.   
      DATA C09T00/0.0,                                                  6668.   
     1    .1151,.2302,.3312,.4172,.4903,.5514,.6016,.6447,.6796,.7102,  6669.   
     2    .7355,.7578,.7769,.7935,.8085,.8212,.8330,.8432,.8524,.8609,  6670.   
     3    .8683,.8752,.8815,.8872,.8926,.8974,.9019,.9061,.9100,.9136,  6671.   
     4    .9170,.9201,.9231,.9258,.9284,.9309,.9332,.9354,.9374,.9394,  6672.   
     5    .9412,.9430,.9446,.9462,.9477,.9492,.9506,.9519,.9531,.9543/  6673.   
C                                                                       6674.   
      DATA C09T01/0.0,                                                  6675.   
     1    .0245,.1526,.2576,.3468,.4239,.4902,.5461,.5952,.6357,.6717,  6676.   
     2    .7017,.7283,.7513,.7712,.7891,.8043,.8183,.8304,.8413,.8512,  6677.   
     3    .8599,.8680,.8753,.8818,.8880,.8934,.8985,.9032,.9075,.9116,  6678.   
     4    .9153,.9187,.9220,.9250,.9278,.9305,.9329,.9353,.9375,.9396,  6679.   
     5    .9415,.9434,.9451,.9468,.9484,.9499,.9513,.9527,.9540,.9552/  6680.   
C                                                                       6681.   
      DATA C09T02/0.0,                                                  6682.   
     1    .0057,.1184,.2173,.3044,.3816,.4494,.5078,.5598,.6035,.6428,  6683.   
     2    .6758,.7053,.7309,.7532,.7733,.7904,.8062,.8197,.8320,.8432,  6684.   
     3    .8529,.8619,.8700,.8772,.8841,.8901,.8956,.9008,.9055,.9099,  6685.   
     4    .9139,.9177,.9212,.9244,.9274,.9302,.9329,.9354,.9377,.9399,  6686.   
     5    .9419,.9439,.9457,.9475,.9491,.9507,.9521,.9535,.9549,.9561/  6687.   
C                                                                       6688.   
      DATA C09T03/0.0,                                                  6689.   
     1    .0177,.1190,.2077,.2880,.3610,.4269,.4847,.5372,.5820,.6227,  6690.   
     2    .6574,.6886,.7157,.7396,.7612,.7796,.7967,.8113,.8246,.8367,  6691.   
     3    .8472,.8570,.8657,.8735,.8809,.8873,.8933,.8989,.9039,.9086,  6692.   
     4    .9129,.9168,.9205,.9239,.9271,.9301,.9329,.9355,.9379,.9402,  6693.   
     5    .9423,.9444,.9462,.9481,.9497,.9514,.9529,.9543,.9557,.9570/  6694.   
C                                                                       6695.   
      DATA C09T04/0.0,                                                  6696.   
     1    .0383,.1335,.2145,.2879,.3553,.4173,.4729,.5241,.5685,.6094,  6697.   
     2    .6446,.6766,.7046,.7294,.7519,.7713,.7891,.8046,.8186,.8314,  6698.   
     3    .8425,.8529,.8621,.8704,.8782,.8850,.8913,.8972,.9025,.9074,  6699.   
     4    .9119,.9161,.9200,.9235,.9269,.9300,.9328,.9356,.9381,.9405,  6700.   
     5    .9427,.9448,.9467,.9486,.9503,.9520,.9535,.9550,.9564,.9577/  6701.   
C                                                                       6702.   
      DATA C09T05/0.0,                                                  6703.   
     1    .0614,.1528,.2288,.2967,.3590,.4167,.4692,.5181,.5613,.6013,  6704.   
     2    .6363,.6684,.6966,.7219,.7449,.7648,.7832,.7993,.8138,.8271,  6705.   
     3    .8387,.8495,.8591,.8678,.8759,.8830,.8896,.8958,.9013,.9064,  6706.   
     4    .9111,.9154,.9195,.9232,.9266,.9298,.9328,.9356,.9382,.9407,  6707.   
     5    .9429,.9451,.9471,.9490,.9508,.9525,.9541,.9556,.9570,.9583/  6708.   
C                                                                       6709.   
      DATA C09T06/0.0,                                                  6710.   
     1    .0849,.1736,.2461,.3098,.3680,.4217,.4710,.5172,.5586,.5974,  6711.   
     2    .6316,.6632,.6913,.7166,.7398,.7599,.7787,.7951,.8100,.8236,  6712.   
     3    .8355,.8467,.8566,.8656,.8740,.8813,.8882,.8945,.9002,.9055,  6713.   
     4    .9104,.9148,.9190,.9228,.9264,.9297,.9328,.9356,.9383,.9408,  6714.   
     5    .9431,.9454,.9474,.9494,.9512,.9529,.9545,.9561,.9575,.9589/  6715.   
C                                                                       6716.   
      DATA C09T07/0.0,                                                  6717.   
     1    .1078,.1944,.2643,.3249,.3797,.4300,.4764,.5199,.5594,.5965,  6718.   
     2    .6296,.6605,.6881,.7132,.7362,.7565,.7753,.7918,.8069,.8208,  6719.   
     3    .8330,.8443,.8545,.8637,.8723,.8799,.8869,.8934,.8992,.9047,  6720.   
     4    .9097,.9143,.9186,.9225,.9262,.9295,.9327,.9356,.9384,.9409,  6721.   
     5    .9433,.9456,.9477,.9497,.9515,.9533,.9549,.9565,.9579,.9593/  6722.   
C                                                                       6723.   
      DATA C09T08/0.0,                                                  6724.   
     1    .1297,.2146,.2824,.3405,.3927,.4402,.4839,.5250,.5625,.5979,  6725.   
     2    .6298,.6597,.6866,.7113,.7340,.7541,.7729,.7895,.8046,.8186,  6726.   
     3    .8309,.8424,.8528,.8621,.8709,.8786,.8858,.8924,.8984,.9040,  6727.   
     4    .9091,.9138,.9182,.9222,.9259,.9294,.9326,.9356,.9384,.9410,  6728.   
     5    .9434,.9457,.9479,.9499,.9518,.9536,.9552,.9568,.9583,.9597/  6729.   
C                                                                       6730.   
      DATA C09T09/0.0,                                                  6731.   
     1    .1505,.2340,.2999,.3561,.4060,.4512,.4927,.5315,.5672,.6009,  6732.   
     2    .6315,.6603,.6865,.7105,.7328,.7526,.7713,.7878,.8029,.8169,  6733.   
     3    .8293,.8409,.8513,.8608,.8697,.8775,.8848,.8916,.8976,.9033,  6734.   
     4    .9085,.9133,.9178,.9219,.9257,.9292,.9325,.9356,.9384,.9411,  6735.   
     5    .9435,.9459,.9480,.9501,.9520,.9538,.9555,.9571,.9586,.9600/  6736.   
C                                                                       6737.   
      DATA C09T10/0.0,                                                  6738.   
     1    .1702,.2523,.3168,.3712,.4193,.4625,.5020,.5390,.5729,.6050,  6739.   
     2    .6344,.6620,.6873,.7108,.7325,.7520,.7703,.7867,.8017,.8157,  6740.   
     3    .8281,.8397,.8502,.8597,.8687,.8766,.8840,.8908,.8970,.9028,  6741.   
     4    .9080,.9129,.9174,.9216,.9254,.9290,.9324,.9355,.9384,.9411,  6742.   
     5    .9436,.9460,.9482,.9502,.9522,.9540,.9557,.9574,.9589,.9603/  6743.   
C                                                                       6744.   
C                                                                       6745.   
      IF(TAU.GT.1.0) THEN                                               6746.   
      CALL HGCLD1(XMU,TAU,G,GG)                                         6747.   
      GO TO 130                                                         6748.   
      ENDIF                                                             6749.   
C                                                                       6750.   
C     ----------------------------------------------------------------  6751.   
C     COSBAR ADJUSTMENT TO REPRODUCE THE SOLAR ZENITH ANGLE DEPENDENCE  6752.   
C     FOR AEROSOL ALBEDOS FOR OPTICAL THICKNESSES OF (0.0 < TAU < 1.0)  6753.   
C     ----------------------------------------------------------------  6754.   
C                                                                       6755.   
C                                                                       6756.   
C                          -------------------------------------------  6757.   
C                          XMU (COSZ) SOLAR ZENITH ANGLE INTERPOLATION  6758.   
C                          DATA INTERVAL:  0.02  ON  (0.0 < XMU < 1.0)  6759.   
C                          -------------------------------------------  6760.   
C                                                                       6761.   
      XI=XMU*50.0+0.9999                                                6762.   
      IX=XI                                                             6763.   
      IF(IX.LT.1) IX=1                                                  6764.   
      JX=IX+1                                                           6765.   
      WXJ=XI-IX                                                         6766.   
      WXI=1.0-WXJ                                                       6767.   
C                                                                       6768.   
C                                            -------------------------  6769.   
C                                            AEROSOL TAU INTERPOLATION  6770.   
C                                            0.10 ON (0.0 < XMU < 1.0)  6771.   
C                                            -------------------------  6772.   
C                                                                       6773.   
      TI=TAU*10.0+0.9999                                                6774.   
      IT=TI                                                             6775.   
      IF(IT.LT.1) IT=1                                                  6776.   
      IF(IT.GT.11) IT=11                                                6777.   
      JT=IT+1                                                           6778.   
      IF(JT.GT.11) JT=11                                                6779.   
      WTJ=TI-IT                                                         6780.   
      WTI=1.0-WTJ                                                       6781.   
C                                                                       6782.   
C                                      -------------------------------  6783.   
C                                      COSBAR DEPENDENCE INTERPOLATION  6784.   
C                                         0.10 ON (0.5 < COSBAR < 0.9)  6785.   
C                                      LINEAR FOR (0.0 < COSBAR < 0.5)  6786.   
C                                      -------------------------------  6787.   
C                                                                       6788.   
      GI=G*10.0                                                         6789.   
      IF(GI.GT.5.0) GO TO 110                                           6790.   
      JG=1                                                              6791.   
      GG=G*(WTI*(WXI*GTAU(IX,IT,JG)+WXJ*GTAU(JX,IT,JG))                 6792.   
     +  +   WTJ*(WXI*GTAU(IX,JT,JG)+WXJ*GTAU(JX,JT,JG)))                6793.   
      GG=GG+GG                                                          6794.   
      GO TO 130                                                         6795.   
C                                                                       6796.   
  110 IG=GI                                                             6797.   
      WGJ=GI-IG                                                         6798.   
      WGI=1.0-WGJ                                                       6799.   
      IG=IG-4                                                           6800.   
      JG=IG+1                                                           6801.   
      IF(IG.GT.4) GO TO 120                                             6802.   
C                                                                       6803.   
      GG=WGI*(WTI*(WXI*GTAU(IX,IT,IG)+WXJ*GTAU(JX,IT,IG))               6804.   
     +      + WTJ*(WXI*GTAU(IX,JT,IG)+WXJ*GTAU(JX,JT,IG)))              6805.   
     +  +WGJ*(WTI*(WXI*GTAU(IX,IT,JG)+WXJ*GTAU(JX,IT,JG))               6806.   
     +      + WTJ*(WXI*GTAU(IX,JT,JG)+WXJ*GTAU(JX,JT,JG)))              6807.   
      GO TO 130                                                         6808.   
C                                                                       6809.   
  120 IG=5                                                              6810.   
C                                                                       6811.   
      GG=WGI*(WTI*(WXI*GTAU(IX,IT,IG)+WXJ*GTAU(JX,IT,IG))               6812.   
     +      + WTJ*(WXI*GTAU(IX,JT,IG)+WXJ*GTAU(JX,JT,IG)))              6813.   
     +  +WGJ                                                            6814.   
C                                                                       6815.   
  130 CONTINUE                                                          6816.   
C                                                                       6817.   
      RETURN                                                            6818.   
      END                                                               6819.   
      SUBROUTINE HGCLD1(XMU,TAU,G,GG)                                   6820.   
C                                                                       6821.   
      DIMENSION C05T01(51),C06T01(51),C07T01(51),C08T01(51),C09T01(51)  6822.   
      DIMENSION C05T02(51),C06T02(51),C07T02(51),C08T02(51),C09T02(51)  6823.   
      DIMENSION C05T03(51),C06T03(51),C07T03(51),C08T03(51),C09T03(51)  6824.   
      DIMENSION C05T04(51),C06T04(51),C07T04(51),C08T04(51),C09T04(51)  6825.   
      DIMENSION C05T05(51),C06T05(51),C07T05(51),C08T05(51),C09T05(51)  6826.   
      DIMENSION C05T06(51),C06T06(51),C07T06(51),C08T06(51),C09T06(51)  6827.   
      DIMENSION C05T07(51),C06T07(51),C07T07(51),C08T07(51),C09T07(51)  6828.   
      DIMENSION C05T08(51),C06T08(51),C07T08(51),C08T08(51),C09T08(51)  6829.   
      DIMENSION C05T09(51),C06T09(51),C07T09(51),C08T09(51),C09T09(51)  6830.   
      DIMENSION C05T10(51),C06T10(51),C07T10(51),C08T10(51),C09T10(51)  6831.   
      DIMENSION C05T99(51),C06T99(51),C07T99(51),C08T99(51),C09T99(51)  6832.   
C                                                                       6833.   
      DIMENSION C05TAU(51,11),C06TAU(51,11),C07TAU(51,11),C08TAU(51,11) 6834.   
      DIMENSION C09TAU(51,11)                                           6835.   
C                                                                       6836.   
      DIMENSION GTAU(51,11,5)                                           6837.   
C                                                                       6838.   
      EQUIVALENCE (C05TAU(1, 1),C05T01(1)),(C05TAU(1, 2),C05T02(1))     6839.   
      EQUIVALENCE (C05TAU(1, 3),C05T03(1)),(C05TAU(1, 4),C05T04(1))     6840.   
      EQUIVALENCE (C05TAU(1, 5),C05T05(1)),(C05TAU(1, 6),C05T06(1))     6841.   
      EQUIVALENCE (C05TAU(1, 7),C05T07(1)),(C05TAU(1, 8),C05T08(1))     6842.   
      EQUIVALENCE (C05TAU(1, 9),C05T09(1)),(C05TAU(1,10),C05T10(1))     6843.   
      EQUIVALENCE (C05TAU(1,11),C05T99(1))                              6844.   
C                                                                       6845.   
      EQUIVALENCE (C06TAU(1, 1),C06T01(1)),(C06TAU(1, 2),C06T02(1))     6846.   
      EQUIVALENCE (C06TAU(1, 3),C06T03(1)),(C06TAU(1, 4),C06T04(1))     6847.   
      EQUIVALENCE (C06TAU(1, 5),C06T05(1)),(C06TAU(1, 6),C06T06(1))     6848.   
      EQUIVALENCE (C06TAU(1, 7),C06T07(1)),(C06TAU(1, 8),C06T08(1))     6849.   
      EQUIVALENCE (C06TAU(1, 9),C06T09(1)),(C06TAU(1,10),C06T10(1))     6850.   
      EQUIVALENCE (C06TAU(1,11),C06T99(1))                              6851.   
C                                                                       6852.   
      EQUIVALENCE (C07TAU(1, 1),C07T01(1)),(C07TAU(1, 2),C07T02(1))     6853.   
      EQUIVALENCE (C07TAU(1, 3),C07T03(1)),(C07TAU(1, 4),C07T04(1))     6854.   
      EQUIVALENCE (C07TAU(1, 5),C07T05(1)),(C07TAU(1, 6),C07T06(1))     6855.   
      EQUIVALENCE (C07TAU(1, 7),C07T07(1)),(C07TAU(1, 8),C07T08(1))     6856.   
      EQUIVALENCE (C07TAU(1, 9),C07T09(1)),(C07TAU(1,10),C07T10(1))     6857.   
      EQUIVALENCE (C07TAU(1,11),C07T99(1))                              6858.   
C                                                                       6859.   
      EQUIVALENCE (C08TAU(1, 1),C08T01(1)),(C08TAU(1, 2),C08T02(1))     6860.   
      EQUIVALENCE (C08TAU(1, 3),C08T03(1)),(C08TAU(1, 4),C08T04(1))     6861.   
      EQUIVALENCE (C08TAU(1, 5),C08T05(1)),(C08TAU(1, 6),C08T06(1))     6862.   
      EQUIVALENCE (C08TAU(1, 7),C08T07(1)),(C08TAU(1, 8),C08T08(1))     6863.   
      EQUIVALENCE (C08TAU(1, 9),C08T09(1)),(C08TAU(1,10),C08T10(1))     6864.   
      EQUIVALENCE (C08TAU(1,11),C08T99(1))                              6865.   
C                                                                       6866.   
      EQUIVALENCE (C09TAU(1, 1),C09T01(1)),(C09TAU(1, 2),C09T02(1))     6867.   
      EQUIVALENCE (C09TAU(1, 3),C09T03(1)),(C09TAU(1, 4),C09T04(1))     6868.   
      EQUIVALENCE (C09TAU(1, 5),C09T05(1)),(C09TAU(1, 6),C09T06(1))     6869.   
      EQUIVALENCE (C09TAU(1, 7),C09T07(1)),(C09TAU(1, 8),C09T08(1))     6870.   
      EQUIVALENCE (C09TAU(1, 9),C09T09(1)),(C09TAU(1,10),C09T10(1))     6871.   
      EQUIVALENCE (C09TAU(1,11),C09T99(1))                              6872.   
C                                                                       6873.   
      EQUIVALENCE (C05TAU(1,1),GTAU(1,1,1)),(C06TAU(1,1),GTAU(1,1,2))   6874.   
      EQUIVALENCE (C07TAU(1,1),GTAU(1,1,3)),(C08TAU(1,1),GTAU(1,1,4))   6875.   
      EQUIVALENCE (C09TAU(1,1),GTAU(1,1,5))                             6876.   
C                                                                       6877.   
C                                                                       6878.   
      DATA C05T01/0.0,                                                  6879.   
     1    .1139,.1346,.1518,.1673,.1817,.1953,.2086,.2216,.2344,.2472,  6880.   
     2    .2600,.2728,.2857,.2985,.3113,.3240,.3367,.3492,.3617,.3740,  6881.   
     3    .3862,.3982,.4100,.4217,.4332,.4444,.4554,.4663,.4769,.4873,  6882.   
     4    .4975,.5075,.5173,.5268,.5362,.5453,.5543,.5630,.5715,.5799,  6883.   
     5    .5880,.5960,.6037,.6113,.6187,.6259,.6330,.6399,.6466,.6532/  6884.   
C                                                                       6885.   
      DATA C05T02/0.0,                                                  6886.   
     1    .1981,.2188,.2361,.2514,.2656,.2788,.2912,.3031,.3145,.3255,  6887.   
     2    .3362,.3466,.3569,.3669,.3768,.3865,.3962,.4057,.4151,.4244,  6888.   
     3    .4337,.4428,.4519,.4609,.4698,.4785,.4872,.4958,.5043,.5127,  6889.   
     4    .5209,.5290,.5371,.5450,.5528,.5604,.5679,.5753,.5826,.5898,  6890.   
     5    .5968,.6037,.6105,.6171,.6237,.6301,.6364,.6425,.6486,.6545/  6891.   
C                                                                       6892.   
      DATA C05T03/0.0,                                                  6893.   
     1    .2435,.2639,.2809,.2960,.3099,.3227,.3348,.3463,.3571,.3676,  6894.   
     2    .3777,.3874,.3969,.4060,.4150,.4237,.4323,.4407,.4489,.4570,  6895.   
     3    .4650,.4728,.4806,.4882,.4957,.5031,.5104,.5177,.5248,.5319,  6896.   
     4    .5388,.5457,.5525,.5592,.5659,.5724,.5788,.5852,.5915,.5977,  6897.   
     5    .6038,.6098,.6157,.6215,.6273,.6330,.6385,.6440,.6494,.6547/  6898.   
C                                                                       6899.   
      DATA C05T04/0.0,                                                  6900.   
     1    .2714,.2914,.3081,.3229,.3365,.3491,.3608,.3719,.3824,.3925,  6901.   
     2    .4022,.4115,.4205,.4292,.4377,.4459,.4540,.4618,.4694,.4769,  6902.   
     3    .4842,.4914,.4985,.5054,.5122,.5189,.5255,.5320,.5384,.5447,  6903.   
     4    .5509,.5570,.5631,.5690,.5749,.5807,.5865,.5921,.5977,.6033,  6904.   
     5    .6087,.6141,.6194,.6246,.6298,.6349,.6399,.6448,.6497,.6545/  6905.   
C                                                                       6906.   
      DATA C05T05/0.0,                                                  6907.   
     1    .2900,.3097,.3262,.3408,.3541,.3664,.3778,.3887,.3989,.4088,  6908.   
     2    .4181,.4272,.4358,.4442,.4524,.4602,.4680,.4754,.4827,.4898,  6909.   
     3    .4967,.5035,.5101,.5166,.5230,.5293,.5354,.5415,.5474,.5533,  6910.   
     4    .5590,.5647,.5703,.5757,.5812,.5865,.5918,.5970,.6021,.6071,  6911.   
     5    .6121,.6171,.6219,.6267,.6315,.6361,.6407,.6453,.6498,.6542/  6912.   
C                                                                       6913.   
      DATA C05T06/0.0,                                                  6914.   
     1    .3033,.3228,.3390,.3534,.3665,.3786,.3898,.4005,.4105,.4201,  6915.   
     2    .4292,.4380,.4465,.4546,.4625,.4701,.4776,.4848,.4918,.4986,  6916.   
     3    .5053,.5118,.5182,.5244,.5305,.5364,.5423,.5480,.5537,.5592,  6917.   
     4    .5646,.5700,.5753,.5804,.5855,.5905,.5955,.6004,.6052,.6099,  6918.   
     5    .6146,.6192,.6237,.6282,.6326,.6370,.6413,.6456,.6498,.6539/  6919.   
C                                                                       6920.   
      DATA C05T07/0.0,                                                  6921.   
     1    .3133,.3325,.3485,.3627,.3757,.3876,.3987,.4092,.4190,.4284,  6922.   
     2    .4374,.4460,.4543,.4622,.4700,.4774,.4846,.4916,.4984,.5051,  6923.   
     3    .5115,.5178,.5240,.5300,.5359,.5416,.5472,.5528,.5582,.5635,  6924.   
     4    .5687,.5738,.5789,.5838,.5887,.5935,.5982,.6029,.6074,.6119,  6925.   
     5    .6164,.6208,.6251,.6293,.6335,.6377,.6418,.6458,.6498,.6537/  6926.   
C                                                                       6927.   
      DATA C05T08/0.0,                                                  6928.   
     1    .3210,.3400,.3559,.3699,.3827,.3945,.4054,.4158,.4255,.4348,  6929.   
     2    .4436,.4521,.4602,.4680,.4756,.4829,.4900,.4968,.5034,.5099,  6930.   
     3    .5162,.5224,.5284,.5342,.5400,.5455,.5510,.5564,.5616,.5667,  6931.   
     4    .5718,.5767,.5816,.5864,.5911,.5957,.6003,.6047,.6091,.6135,  6932.   
     5    .6177,.6219,.6261,.6302,.6342,.6381,.6421,.6459,.6497,.6535/  6933.   
C                                                                       6934.   
      DATA C05T09/0.0,                                                  6935.   
     1    .3271,.3460,.3618,.3757,.3883,.4000,.4108,.4211,.4306,.4398,  6936.   
     2    .4485,.4569,.4649,.4726,.4800,.4872,.4941,.5008,.5074,.5137,  6937.   
     3    .5199,.5259,.5318,.5375,.5431,.5486,.5539,.5591,.5642,.5693,  6938.   
     4    .5742,.5790,.5837,.5884,.5930,.5974,.6018,.6062,.6104,.6146,  6939.   
     5    .6188,.6228,.6268,.6308,.6347,.6385,.6423,.6460,.6497,.6533/  6940.   
C                                                                       6941.   
      DATA C05T10/0.0,                                                  6942.   
     1    .3321,.3509,.3665,.3803,.3929,.4045,.4152,.4253,.4348,.4439,  6943.   
     2    .4525,.4607,.4686,.4762,.4836,.4906,.4975,.5041,.5105,.5168,  6944.   
     3    .5229,.5288,.5345,.5401,.5457,.5510,.5562,.5614,.5664,.5713,  6945.   
     4    .5761,.5808,.5854,.5900,.5944,.5988,.6031,.6073,.6115,.6156,  6946.   
     5    .6196,.6236,.6275,.6313,.6351,.6388,.6425,.6461,.6497,.6532/  6947.   
C                                                                       6948.   
      DATA C05T99/0.0,                                                  6949.   
     1    .3759,.3933,.4078,.4204,.4320,.4425,.4522,.4614,.4699,.4781,  6950.   
     2    .4857,.4930,.5000,.5067,.5131,.5192,.5252,.5309,.5364,.5417,  6951.   
     3    .5469,.5519,.5568,.5615,.5661,.5705,.5749,.5791,.5832,.5873,  6952.   
     4    .5912,.5950,.5988,.6024,.6060,.6095,.6130,.6164,.6196,.6229,  6953.   
     5    .6260,.6292,.6322,.6352,.6381,.6410,.6439,.6467,.6494,.6521/  6954.   
C                                                                       6955.   
      DATA C06T01/0.0,                                                  6956.   
     1    .1314,.1581,.1806,.2009,.2199,.2379,.2553,.2722,.2889,.3055,  6957.   
     2    .3218,.3381,.3541,.3700,.3858,.4012,.4165,.4315,.4462,.4606,  6958.   
     3    .4746,.4884,.5018,.5148,.5277,.5400,.5520,.5638,.5751,.5862,  6959.   
     4    .5969,.6073,.6175,.6272,.6367,.6459,.6548,.6635,.6719,.6800,  6960.   
     5    .6879,.6955,.7029,.7101,.7170,.7237,.7303,.7366,.7427,.7487/  6961.   
C                                                                       6962.   
      DATA C06T02/0.0,                                                  6963.   
     1    .2301,.2561,.2779,.2973,.3151,.3317,.3472,.3620,.3761,.3897,  6964.   
     2    .4028,.4155,.4279,.4399,.4518,.4633,.4747,.4858,.4968,.5076,  6965.   
     3    .5182,.5287,.5389,.5490,.5589,.5686,.5781,.5875,.5967,.6057,  6966.   
     4    .6144,.6230,.6315,.6397,.6478,.6556,.6633,.6708,.6781,.6853,  6967.   
     5    .6922,.6991,.7057,.7121,.7184,.7246,.7306,.7364,.7421,.7476/  6968.   
C                                                                       6969.   
      DATA C06T03/0.0,                                                  6970.   
     1    .2848,.3100,.3311,.3497,.3668,.3825,.3971,.4110,.4240,.4365,  6971.   
     2    .4484,.4599,.4710,.4816,.4921,.5021,.5119,.5214,.5308,.5399,  6972.   
     3    .5488,.5575,.5661,.5745,.5828,.5908,.5988,.6066,.6142,.6217,  6973.   
     4    .6291,.6364,.6435,.6505,.6574,.6641,.6707,.6772,.6835,.6898,  6974.   
     5    .6959,.7019,.7077,.7135,.7191,.7246,.7300,.7353,.7404,.7455/  6975.   
C                                                                       6976.   
      DATA C06T04/0.0,                                                  6977.   
     1    .3189,.3434,.3639,.3819,.3983,.4134,.4273,.4406,.4529,.4647,  6978.   
     2    .4759,.4867,.4970,.5069,.5165,.5258,.5348,.5435,.5519,.5602,  6979.   
     3    .5682,.5761,.5837,.5912,.5985,.6057,.6127,.6196,.6263,.6330,  6980.   
     4    .6395,.6459,.6521,.6583,.6644,.6703,.6761,.6819,.6875,.6931,  6981.   
     5    .6985,.7039,.7091,.7143,.7194,.7243,.7292,.7340,.7387,.7433/  6982.   
C                                                                       6983.   
      DATA C06T05/0.0,                                                  6984.   
     1    .3420,.3660,.3859,.4034,.4193,.4339,.4474,.4601,.4720,.4833,  6985.   
     2    .4940,.5043,.5141,.5235,.5326,.5413,.5498,.5579,.5658,.5736,  6986.   
     3    .5810,.5883,.5954,.6023,.6091,.6157,.6221,.6285,.6346,.6407,  6987.   
     4    .6466,.6525,.6582,.6638,.6693,.6747,.6800,.6853,.6904,.6955,  6988.   
     5    .7004,.7053,.7101,.7148,.7194,.7240,.7285,.7329,.7372,.7415/  6989.   
C                                                                       6990.   
      DATA C06T06/0.0,                                                  6991.   
     1    .3586,.3821,.4016,.4187,.4342,.4484,.4615,.4739,.4854,.4964,  6992.   
     2    .5067,.5166,.5260,.5350,.5438,.5521,.5602,.5680,.5755,.5829,  6993.   
     3    .5899,.5968,.6036,.6101,.6165,.6227,.6287,.6347,.6405,.6462,  6994.   
     4    .6517,.6571,.6625,.6677,.6729,.6779,.6828,.6877,.6925,.6972,  6995.   
     5    .7018,.7063,.7108,.7152,.7195,.7237,.7279,.7320,.7360,.7400/  6996.   
C                                                                       6997.   
      DATA C06T07/0.0,                                                  6998.   
     1    .3711,.3942,.4133,.4301,.4453,.4592,.4720,.4841,.4953,.5060,  6999.   
     2    .5160,.5256,.5348,.5435,.5520,.5600,.5678,.5753,.5826,.5896,  7000.   
     3    .5964,.6031,.6095,.6157,.6219,.6278,.6336,.6392,.6447,.6501,  7001.   
     4    .6554,.6606,.6657,.6706,.6755,.6802,.6849,.6895,.6940,.6985,  7002.   
     5    .7028,.7071,.7113,.7154,.7195,.7235,.7274,.7313,.7351,.7388/  7003.   
C                                                                       7004.   
      DATA C06T08/0.0,                                                  7005.   
     1    .3808,.4036,.4224,.4390,.4539,.4676,.4801,.4920,.5029,.5134,  7006.   
     2    .5232,.5326,.5415,.5500,.5582,.5660,.5736,.5809,.5880,.5948,  7007.   
     3    .6014,.6078,.6140,.6200,.6259,.6316,.6372,.6427,.6480,.6532,  7008.   
     4    .6582,.6632,.6681,.6728,.6775,.6820,.6865,.6909,.6952,.6994,  7009.   
     5    .7036,.7077,.7117,.7156,.7195,.7233,.7270,.7307,.7343,.7379/  7010.   
C                                                                       7011.   
      DATA C06T09/0.0,                                                  7012.   
     1    .3886,.4111,.4297,.4460,.4607,.4742,.4865,.4982,.5089,.5192,  7013.   
     2    .5288,.5380,.5467,.5551,.5631,.5708,.5782,.5853,.5922,.5988,  7014.   
     3    .6052,.6115,.6175,.6234,.6291,.6347,.6401,.6454,.6505,.6555,  7015.   
     4    .6604,.6652,.6699,.6745,.6790,.6834,.6877,.6920,.6961,.7002,  7016.   
     5    .7042,.7081,.7119,.7157,.7195,.7231,.7267,.7303,.7337,.7372/  7017.   
C                                                                       7018.   
      DATA C06T10/0.0,                                                  7019.   
     1    .3949,.4172,.4356,.4517,.4663,.4796,.4917,.5032,.5138,.5239,  7020.   
     2    .5334,.5424,.5510,.5592,.5671,.5746,.5819,.5888,.5955,.6021,  7021.   
     3    .6083,.6144,.6203,.6261,.6317,.6371,.6424,.6475,.6525,.6574,  7022.   
     4    .6622,.6668,.6714,.6759,.6802,.6845,.6887,.6928,.6968,.7008,  7023.   
     5    .7046,.7085,.7122,.7159,.7195,.7230,.7265,.7299,.7333,.7366/  7024.   
C                                                                       7025.   
      DATA C06T99/0.0,                                                  7026.   
     1    .4509,.4707,.4871,.5013,.5141,.5256,.5362,.5461,.5551,.5638,  7027.   
     2    .5718,.5794,.5866,.5934,.6000,.6062,.6122,.6178,.6233,.6286,  7028.   
     3    .6336,.6386,.6433,.6478,.6523,.6565,.6607,.6647,.6686,.6724,  7029.   
     4    .6761,.6797,.6832,.6866,.6900,.6932,.6964,.6995,.7025,.7055,  7030.   
     5    .7084,.7112,.7140,.7167,.7194,.7220,.7245,.7270,.7295,.7319/  7031.   
C                                                                       7032.   
      DATA C07T01/0.0,                                                  7033.   
     1    .1456,.1813,.2114,.2384,.2635,.2872,.3099,.3319,.3532,.3742,  7034.   
     2    .3946,.4147,.4342,.4533,.4720,.4901,.5078,.5248,.5413,.5573,  7035.   
     3    .5727,.5876,.6019,.6156,.6290,.6417,.6539,.6657,.6770,.6879,  7036.   
     4    .6982,.7082,.7178,.7270,.7358,.7442,.7523,.7602,.7676,.7748,  7037.   
     5    .7816,.7883,.7946,.8007,.8066,.8122,.8176,.8228,.8278,.8326/  7038.   
C                                                                       7039.   
      DATA C07T02/0.0,                                                  7040.   
     1    .2601,.2939,.3219,.3466,.3691,.3898,.4090,.4272,.4442,.4606,  7041.   
     2    .4762,.4912,.5057,.5198,.5334,.5466,.5596,.5721,.5843,.5963,  7042.   
     3    .6078,.6192,.6302,.6410,.6515,.6616,.6715,.6811,.6904,.6995,  7043.   
     4    .7083,.7168,.7251,.7331,.7409,.7483,.7556,.7626,.7694,.7760,  7044.   
     5    .7824,.7885,.7945,.8002,.8058,.8111,.8163,.8214,.8262,.8309/  7045.   
C                                                                       7046.   
      DATA C07T03/0.0,                                                  7047.   
     1    .3256,.3578,.3842,.4074,.4283,.4473,.4648,.4813,.4966,.5111,  7048.   
     2    .5248,.5379,.5504,.5624,.5740,.5851,.5959,.6063,.6163,.6262,  7049.   
     3    .6357,.6450,.6540,.6628,.6715,.6798,.6880,.6960,.7037,.7113,  7050.   
     4    .7187,.7259,.7330,.7398,.7465,.7530,.7594,.7656,.7716,.7774,  7051.   
     5    .7831,.7887,.7940,.7993,.8044,.8093,.8141,.8188,.8233,.8278/  7052.   
C                                                                       7053.   
      DATA C07T04/0.0,                                                  7054.   
     1    .3675,.3983,.4235,.4455,.4652,.4831,.4995,.5149,.5290,.5424,  7055.   
     2    .5550,.5670,.5783,.5892,.5996,.6096,.6192,.6284,.6374,.6461,  7056.   
     3    .6544,.6626,.6705,.6781,.6857,.6929,.7000,.7070,.7137,.7204,  7057.   
     4    .7268,.7331,.7393,.7453,.7512,.7569,.7625,.7680,.7734,.7786,  7058.   
     5    .7837,.7887,.7936,.7983,.8030,.8075,.8119,.8163,.8205,.8246/  7059.   
C                                                                       7060.   
      DATA C07T05/0.0,                                                  7061.   
     1    .3963,.4260,.4503,.4714,.4902,.5073,.5228,.5374,.5507,.5634,  7062.   
     2    .5752,.5864,.5970,.6071,.6168,.6260,.6349,.6434,.6516,.6596,  7063.   
     3    .6672,.6746,.6818,.6888,.6956,.7022,.7086,.7149,.7210,.7270,  7064.   
     4    .7328,.7384,.7440,.7494,.7547,.7599,.7650,.7699,.7748,.7796,  7065.   
     5    .7842,.7887,.7932,.7976,.8018,.8060,.8101,.8141,.8180,.8218/  7066.   
C                                                                       7067.   
      DATA C07T06/0.0,                                                  7068.   
     1    .4172,.4461,.4696,.4900,.5082,.5246,.5395,.5535,.5662,.5783,  7069.   
     2    .5895,.6001,.6102,.6198,.6289,.6376,.6460,.6540,.6617,.6691,  7070.   
     3    .6763,.6832,.6899,.6964,.7028,.7089,.7148,.7206,.7263,.7318,  7071.   
     4    .7371,.7424,.7475,.7525,.7574,.7622,.7668,.7714,.7759,.7803,  7072.   
     5    .7846,.7888,.7929,.7969,.8009,.8048,.8086,.8123,.8159,.8195/  7073.   
C                                                                       7074.   
      DATA C07T07/0.0,                                                  7075.   
     1    .4331,.4613,.4842,.5040,.5216,.5375,.5520,.5654,.5777,.5893,  7076.   
     2    .6001,.6104,.6200,.6291,.6379,.6462,.6542,.6618,.6691,.6762,  7077.   
     3    .6830,.6896,.6959,.7021,.7081,.7138,.7194,.7249,.7302,.7354,  7078.   
     4    .7404,.7453,.7502,.7548,.7594,.7639,.7683,.7726,.7768,.7809,  7079.   
     5    .7849,.7888,.7927,.7965,.8002,.8038,.8074,.8109,.8143,.8177/  7080.   
C                                                                       7081.   
      DATA C07T08/0.0,                                                  7082.   
     1    .4455,.4731,.4955,.5148,.5320,.5475,.5616,.5747,.5866,.5979,  7083.   
     2    .6083,.6182,.6275,.6363,.6448,.6528,.6605,.6678,.6748,.6816,  7084.   
     3    .6881,.6944,.7005,.7064,.7121,.7176,.7230,.7282,.7332,.7382,  7085.   
     4    .7430,.7476,.7522,.7566,.7610,.7652,.7694,.7735,.7774,.7813,  7086.   
     5    .7851,.7889,.7925,.7961,.7996,.8030,.8064,.8097,.8130,.8162/  7087.   
C                                                                       7088.   
      DATA C07T09/0.0,                                                  7089.   
     1    .4555,.4826,.5046,.5235,.5404,.5555,.5692,.5820,.5936,.6046,  7090.   
     2    .6147,.6244,.6334,.6420,.6502,.6579,.6654,.6725,.6793,.6859,  7091.   
     3    .6921,.6982,.7041,.7098,.7153,.7206,.7257,.7308,.7356,.7404,  7092.   
     4    .7449,.7494,.7538,.7581,.7622,.7663,.7703,.7742,.7780,.7817,  7093.   
     5    .7853,.7889,.7924,.7958,.7992,.8024,.8057,.8088,.8119,.8150/  7094.   
C                                                                       7095.   
      DATA C07T10/0.0,                                                  7096.   
     1    .4637,.4903,.5120,.5306,.5471,.5620,.5754,.5879,.5993,.6101,  7097.   
     2    .6200,.6294,.6382,.6466,.6546,.6621,.6694,.6763,.6829,.6893,  7098.   
     3    .6954,.7013,.7070,.7125,.7179,.7230,.7280,.7328,.7375,.7421,  7099.   
     4    .7465,.7509,.7551,.7592,.7632,.7672,.7710,.7747,.7784,.7820,  7100.   
     5    .7855,.7889,.7923,.7956,.7988,.8020,.8051,.8081,.8111,.8140/  7101.   
C                                                                       7102.   
      DATA C07T99/0.0,                                                  7103.   
     1    .5366,.5590,.5770,.5924,.6060,.6180,.6289,.6389,.6480,.6565,  7104.   
     2    .6643,.6717,.6785,.6850,.6912,.6969,.7025,.7077,.7127,.7175,  7105.   
     3    .7220,.7264,.7306,.7347,.7386,.7423,.7460,.7495,.7529,.7562,  7106.   
     4    .7594,.7625,.7655,.7684,.7712,.7740,.7767,.7793,.7818,.7843,  7107.   
     5    .7867,.7891,.7914,.7937,.7959,.7981,.8002,.8022,.8043,.8062/  7108.   
C                                                                       7109.   
      DATA C08T01/0.0,                                                  7110.   
     1    .1567,.2073,.2493,.2865,.3206,.3522,.3819,.4104,.4374,.4636,  7111.   
     2    .4886,.5128,.5359,.5580,.5793,.5994,.6188,.6370,.6544,.6710,  7112.   
     3    .6864,.7013,.7152,.7283,.7410,.7526,.7637,.7742,.7840,.7935,  7113.   
     4    .8023,.8106,.8185,.8259,.8330,.8396,.8459,.8519,.8576,.8630,  7114.   
     5    .8680,.8729,.8774,.8818,.8859,.8899,.8936,.8972,.9006,.9038/  7115.   
C                                                                       7116.   
      DATA C08T02/0.0,                                                  7117.   
     1    .2878,.3342,.3718,.4041,.4329,.4588,.4824,.5045,.5249,.5442,  7118.   
     2    .5623,.5797,.5962,.6120,.6272,.6417,.6559,.6693,.6823,.6949,  7119.   
     3    .7069,.7186,.7298,.7405,.7509,.7606,.7701,.7792,.7879,.7963,  7120.   
     4    .8042,.8118,.8191,.8260,.8327,.8390,.8451,.8509,.8564,.8617,  7121.   
     5    .8667,.8716,.8762,.8806,.8848,.8888,.8926,.8963,.8998,.9032/  7122.   
C                                                                       7123.   
      DATA C08T03/0.0,                                                  7124.   
     1    .3656,.4087,.4432,.4725,.4984,.5215,.5422,.5614,.5789,.5954,  7125.   
     2    .6106,.6251,.6387,.6517,.6641,.6758,.6872,.6981,.7085,.7187,  7126.   
     3    .7283,.7378,.7468,.7555,.7641,.7722,.7801,.7878,.7951,.8022,  7127.   
     4    .8091,.8157,.8221,.8282,.8342,.8399,.8454,.8507,.8558,.8608,  7128.   
     5    .8655,.8700,.8744,.8786,.8826,.8865,.8903,.8939,.8973,.9006/  7129.   
C                                                                       7130.   
      DATA C08T04/0.0,                                                  7131.   
     1    .4167,.4573,.4895,.5167,.5405,.5616,.5805,.5979,.6136,.6283,  7132.   
     2    .6419,.6547,.6668,.6781,.6890,.6992,.7091,.7184,.7274,.7361,  7133.   
     3    .7444,.7525,.7602,.7677,.7750,.7820,.7888,.7954,.8018,.8080,  7134.   
     4    .8139,.8197,.8254,.8308,.8361,.8412,.8462,.8510,.8556,.8601,  7135.   
     5    .8645,.8687,.8728,.8767,.8805,.8842,.8877,.8912,.8945,.8977/  7136.   
C                                                                       7137.   
      DATA C08T05/0.0,                                                  7138.   
     1    .4528,.4913,.5218,.5473,.5696,.5893,.6069,.6230,.6375,.6511,  7139.   
     2    .6635,.6752,.6862,.6965,.7063,.7156,.7245,.7329,.7409,.7487,  7140.   
     3    .7561,.7633,.7703,.7769,.7834,.7896,.7957,.8015,.8072,.8127,  7141.   
     4    .8180,.8232,.8283,.8332,.8379,.8426,.8470,.8514,.8556,.8598,  7142.   
     5    .8638,.8677,.8714,.8751,.8787,.8821,.8855,.8887,.8919,.8950/  7143.   
C                                                                       7144.   
      DATA C08T06/0.0,                                                  7145.   
     1    .4795,.5164,.5454,.5697,.5909,.6095,.6261,.6412,.6548,.6675,  7146.   
     2    .6791,.6901,.7003,.7098,.7190,.7275,.7357,.7435,.7509,.7581,  7147.   
     3    .7648,.7714,.7778,.7838,.7898,.7954,.8009,.8063,.8115,.8165,  7148.   
     4    .8214,.8261,.8307,.8352,.8395,.8437,.8479,.8519,.8558,.8596,  7149.   
     5    .8633,.8669,.8704,.8738,.8772,.8804,.8836,.8866,.8896,.8925/  7150.   
C                                                                       7151.   
      DATA C08T07/0.0,                                                  7152.   
     1    .5000,.5356,.5635,.5868,.6070,.6248,.6406,.6550,.6679,.6800,  7153.   
     2    .6909,.7013,.7109,.7199,.7285,.7365,.7442,.7515,.7584,.7651,  7154.   
     3    .7715,.7776,.7835,.7892,.7947,.7999,.8051,.8100,.8148,.8195,  7155.   
     4    .8240,.8284,.8327,.8368,.8408,.8448,.8486,.8523,.8560,.8595,  7156.   
     5    .8630,.8663,.8696,.8728,.8759,.8790,.8820,.8849,.8877,.8905/  7157.   
C                                                                       7158.   
      DATA C08T08/0.0,                                                  7159.   
     1    .5162,.5507,.5777,.6002,.6197,.6368,.6519,.6657,.6781,.6896,  7160.   
     2    .7001,.7100,.7191,.7277,.7359,.7435,.7508,.7577,.7643,.7706,  7161.   
     3    .7766,.7824,.7880,.7933,.7986,.8035,.8083,.8130,.8175,.8219,  7162.   
     4    .8261,.8302,.8343,.8381,.8419,.8456,.8492,.8527,.8561,.8595,  7163.   
     5    .8627,.8659,.8690,.8720,.8750,.8778,.8806,.8834,.8861,.8887/  7164.   
C                                                                       7165.   
      DATA C08T09/0.0,                                                  7166.   
     1    .5293,.5629,.5891,.6109,.6298,.6464,.6610,.6743,.6862,.6974,  7167.   
     2    .7074,.7169,.7257,.7340,.7418,.7491,.7561,.7627,.7690,.7750,  7168.   
     3    .7807,.7863,.7916,.7967,.8016,.8063,.8109,.8154,.8196,.8238,  7169.   
     4    .8278,.8317,.8356,.8392,.8428,.8463,.8497,.8531,.8563,.8595,  7170.   
     5    .8625,.8656,.8685,.8714,.8742,.8769,.8796,.8822,.8847,.8872/  7171.   
C                                                                       7172.   
      DATA C08T10/0.0,                                                  7173.   
     1    .5401,.5729,.5985,.6197,.6381,.6542,.6684,.6813,.6929,.7036,  7174.   
     2    .7134,.7226,.7311,.7390,.7466,.7536,.7604,.7667,.7728,.7786,  7175.   
     3    .7841,.7894,.7945,.7994,.8042,.8087,.8131,.8173,.8214,.8254,  7176.   
     4    .8292,.8330,.8366,.8401,.8436,.8469,.8502,.8534,.8564,.8595,  7177.   
     5    .8624,.8653,.8681,.8708,.8735,.8761,.8787,.8812,.8836,.8860/  7178.   
C                                                                       7179.   
      DATA C08T99/0.0,                                                  7180.   
     1    .6384,.6631,.6821,.6978,.7111,.7227,.7328,.7420,.7501,.7576,  7181.   
     2    .7644,.7707,.7765,.7819,.7870,.7918,.7963,.8005,.8045,.8084,  7182.   
     3    .8120,.8154,.8187,.8219,.8250,.8278,.8307,.8334,.8360,.8385,  7183.   
     4    .8409,.8432,.8455,.8477,.8498,.8519,.8539,.8559,.8578,.8596,  7184.   
     5    .8614,.8632,.8648,.8665,.8681,.8697,.8712,.8728,.8742,.8757/  7185.   
C                                                                       7186.   
      DATA C09T01/0.0,                                                  7187.   
     1    .1702,.2523,.3168,.3712,.4193,.4625,.5020,.5390,.5729,.6050,  7188.   
     2    .6344,.6620,.6873,.7108,.7325,.7520,.7703,.7867,.8017,.8157,  7189.   
     3    .8281,.8397,.8502,.8597,.8687,.8766,.8840,.8908,.8970,.9028,  7190.   
     4    .9080,.9129,.9174,.9216,.9254,.9290,.9324,.9355,.9384,.9411,  7191.   
     5    .9436,.9460,.9482,.9502,.9522,.9540,.9557,.9574,.9589,.9603/  7192.   
C                                                                       7193.   
      DATA C09T02/0.0,                                                  7194.   
     1    .3174,.3895,.4438,.4879,.5256,.5583,.5872,.6136,.6374,.6597,  7195.   
     2    .6802,.6995,.7175,.7345,.7506,.7655,.7798,.7930,.8055,.8173,  7196.   
     3    .8281,.8385,.8481,.8570,.8655,.8731,.8804,.8872,.8935,.8994,  7197.   
     4    .9049,.9099,.9148,.9191,.9233,.9271,.9307,.9341,.9373,.9402,  7198.   
     5    .9430,.9456,.9480,.9503,.9524,.9544,.9563,.9581,.9598,.9613/  7199.   
C                                                                       7200.   
      DATA C09T03/0.0,                                                  7201.   
     1    .4078,.4729,.5209,.5592,.5915,.6191,.6431,.6649,.6842,.7022,  7202.   
     2    .7185,.7339,.7481,.7614,.7741,.7859,.7972,.8078,.8178,.8274,  7203.   
     3    .8364,.8451,.8532,.8608,.8682,.8750,.8815,.8877,.8934,.8989,  7204.   
     4    .9040,.9089,.9135,.9177,.9218,.9256,.9292,.9326,.9358,.9388,  7205.   
     5    .9416,.9443,.9468,.9491,.9514,.9535,.9554,.9573,.9591,.9607/  7206.   
C                                                                       7207.   
      DATA C09T04/0.0,                                                  7208.   
     1    .4692,.5288,.5723,.6066,.6353,.6597,.6807,.6997,.7163,.7318,  7209.   
     2    .7457,.7588,.7708,.7821,.7927,.8026,.8121,.8210,.8295,.8376,  7210.   
     3    .8452,.8525,.8595,.8661,.8724,.8784,.8841,.8896,.8948,.8998,  7211.   
     4    .9044,.9089,.9132,.9172,.9210,.9247,.9281,.9314,.9345,.9374,  7212.   
     5    .9402,.9429,.9453,.9477,.9500,.9521,.9541,.9560,.9579,.9596/  7213.   
C                                                                       7214.   
      DATA C09T05/0.0,                                                  7215.   
     1    .5136,.5690,.6090,.6404,.6666,.6886,.7076,.7246,.7394,.7532,  7216.   
     2    .7655,.7771,.7877,.7976,.8069,.8156,.8239,.8316,.8390,.8461,  7217.   
     3    .8528,.8592,.8653,.8711,.8767,.8820,.8871,.8920,.8967,.9012,  7218.   
     4    .9054,.9095,.9134,.9171,.9207,.9241,.9274,.9305,.9335,.9363,  7219.   
     5    .9390,.9416,.9440,.9464,.9486,.9507,.9527,.9546,.9565,.9582/  7220.   
C                                                                       7221.   
      DATA C09T06/0.0,                                                  7222.   
     1    .5473,.5993,.6366,.6658,.6900,.7102,.7277,.7432,.7568,.7693,  7223.   
     2    .7805,.7910,.8006,.8095,.8179,.8257,.8332,.8401,.8468,.8531,  7224.   
     3    .8591,.8648,.8703,.8755,.8806,.8853,.8899,.8944,.8986,.9027,  7225.   
     4    .9066,.9103,.9140,.9174,.9207,.9239,.9270,.9299,.9327,.9354,  7226.   
     5    .9380,.9405,.9429,.9451,.9473,.9494,.9514,.9533,.9551,.9568/  7227.   
C                                                                       7228.   
      DATA C09T07/0.0,                                                  7229.   
     1    .5737,.6230,.6581,.6855,.7081,.7271,.7433,.7577,.7703,.7819,  7230.   
     2    .7922,.8019,.8107,.8189,.8266,.8338,.8406,.8470,.8530,.8588,  7231.   
     3    .8643,.8695,.8745,.8793,.8839,.8883,.8925,.8966,.9004,.9042,  7232.   
     4    .9078,.9113,.9146,.9178,.9209,.9239,.9268,.9295,.9322,.9348,  7233.   
     5    .9372,.9396,.9419,.9441,.9462,.9482,.9502,.9520,.9538,.9555/  7234.   
C                                                                       7235.   
      DATA C09T08/0.0,                                                  7236.   
     1    .5950,.6420,.6754,.7013,.7226,.7405,.7557,.7693,.7811,.7919,  7237.   
     2    .8016,.8106,.8188,.8265,.8337,.8403,.8466,.8525,.8582,.8635,  7238.   
     3    .8686,.8734,.8781,.8825,.8868,.8908,.8947,.8985,.9021,.9056,  7239.   
     4    .9089,.9121,.9153,.9183,.9212,.9240,.9267,.9293,.9318,.9343,  7240.   
     5    .9366,.9389,.9411,.9432,.9452,.9472,.9490,.9509,.9526,.9543/  7241.   
C                                                                       7242.   
      DATA C09T09/0.0,                                                  7243.   
     1    .6125,.6576,.6894,.7142,.7345,.7514,.7659,.7787,.7899,.8001,  7244.   
     2    .8093,.8177,.8255,.8327,.8394,.8457,.8516,.8572,.8624,.8675,  7245.   
     3    .8722,.8767,.8811,.8852,.8892,.8930,.8966,.9002,.9035,.9068,  7246.   
     4    .9100,.9130,.9159,.9187,.9215,.9241,.9267,.9292,.9316,.9339,  7247.   
     5    .9361,.9383,.9404,.9424,.9443,.9462,.9481,.9498,.9515,.9532/  7248.   
C                                                                       7249.   
      DATA C09T10/0.0,                                                  7250.   
     1    .6272,.6706,.7012,.7249,.7443,.7605,.7743,.7866,.7972,.8069,  7251.   
     2    .8156,.8236,.8310,.8378,.8442,.8501,.8558,.8610,.8660,.8708,  7252.   
     3    .8752,.8795,.8836,.8875,.8913,.8949,.8983,.9016,.9048,.9079,  7253.   
     4    .9109,.9137,.9165,.9192,.9218,.9243,.9267,.9291,.9314,.9336,  7254.   
     5    .9357,.9378,.9398,.9417,.9436,.9454,.9472,.9489,.9506,.9522/  7255.   
C                                                                       7256.   
      DATA C09T99/0.0,                                                  7257.   
     1    .7681,.7934,.8109,.8243,.8350,.8439,.8514,.8579,.8636,.8687,  7258.   
     2    .8732,.8774,.8812,.8847,.8880,.8910,.8938,.8964,.8989,.9013,  7259.   
     3    .9035,.9056,.9076,.9095,.9113,.9130,.9147,.9163,.9178,.9193,  7260.   
     4    .9207,.9221,.9234,.9247,.9260,.9271,.9283,.9294,.9305,.9316,  7261.   
     5    .9326,.9336,.9346,.9355,.9364,.9373,.9382,.9390,.9398,.9406/  7262.   
C                                                                       7263.   
C                                                                       7264.   
C     ----------------------------------------------------------------  7265.   
C     COSBAR ADJUSTMENT TO REPRODUCE THE SOLAR ZENITH ANGLE DEPENDENCE  7266.   
C     FOR CLOUD ALBEDOS FOR OPTICAL THICKNESS FROM (1.0 < TAU < 99.0)   7267.   
C     ----------------------------------------------------------------  7268.   
C                                                                       7269.   
C                                                                       7270.   
C                          -------------------------------------------  7271.   
C                          XMU (COSZ) SOLAR ZENITH ANGLE INTERPOLATION  7272.   
C                          DATA INTERVAL:  0.02  ON  (0.0 < XMU < 1.0)  7273.   
C                          -------------------------------------------  7274.   
C                                                                       7275.   
      XI=XMU*50.0+0.9999                                                7276.   
      IX=XI                                                             7277.   
      IF(IX.LT.1) IX=1                                                  7278.   
      JX=IX+1                                                           7279.   
      WXJ=XI-IX                                                         7280.   
      WXI=1.0-WXJ                                                       7281.   
C                                                                       7282.   
C                                              -----------------------  7283.   
C                                              CLOUD TAU INTERPOLATION  7284.   
C                                              1.0 OVER (1 < TAU < 10)  7285.   
C                                              LINEAR (10 < TAU < 100)  7286.   
C                                              -----------------------  7287.   
C                                                                       7288.   
      TI=TAU                                                            7289.   
      IT=TI                                                             7290.   
      IF(IT.LT.1) IT=1                                                  7291.   
      WTJ=TI-IT                                                         7292.   
      IF(IT.GT.9) THEN                                                  7293.   
      WTJ=(TAU-10.0)/90.0                                               7294.   
      IT=10                                                             7295.   
      ENDIF                                                             7296.   
      WTI=1.0-WTJ                                                       7297.   
      JT=IT+1                                                           7298.   
C                                                                       7299.   
C                                      -------------------------------  7300.   
C                                      COSBAR DEPENDENCE INTERPOLATION  7301.   
C                                         0.10 ON (0.5 < COSBAR < 0.9)  7302.   
C                                      LINEAR FOR (0.0 < COSBAR < 0.5)  7303.   
C                                      -------------------------------  7304.   
C                                                                       7305.   
      GI=G*10.0                                                         7306.   
      IF(GI.GT.5.0) GO TO 110                                           7307.   
      JG=1                                                              7308.   
      GG=G*(WTI*(WXI*GTAU(IX,IT,JG)+WXJ*GTAU(JX,IT,JG))                 7309.   
     +  +   WTJ*(WXI*GTAU(IX,JT,JG)+WXJ*GTAU(JX,JT,JG)))                7310.   
      GG=GG+GG                                                          7311.   
      GO TO 130                                                         7312.   
C                                                                       7313.   
  110 IG=GI                                                             7314.   
      WGJ=GI-IG                                                         7315.   
      WGI=1.0-WGJ                                                       7316.   
      IG=IG-4                                                           7317.   
      JG=IG+1                                                           7318.   
      IF(IG.GT.4) GO TO 120                                             7319.   
C                                                                       7320.   
      GG=WGI*(WTI*(WXI*GTAU(IX,IT,IG)+WXJ*GTAU(JX,IT,IG))               7321.   
     +      + WTJ*(WXI*GTAU(IX,JT,IG)+WXJ*GTAU(JX,JT,IG)))              7322.   
     +  +WGJ*(WTI*(WXI*GTAU(IX,IT,JG)+WXJ*GTAU(JX,IT,JG))               7323.   
     +      + WTJ*(WXI*GTAU(IX,JT,JG)+WXJ*GTAU(JX,JT,JG)))              7324.   
      GO TO 130                                                         7325.   
C                                                                       7326.   
  120 IG=5                                                              7327.   
C                                                                       7328.   
      GG=WGI*(WTI*(WXI*GTAU(IX,IT,IG)+WXJ*GTAU(JX,IT,IG))               7329.   
     +      + WTJ*(WXI*GTAU(IX,JT,IG)+WXJ*GTAU(JX,JT,IG)))              7330.   
     +  +WGJ                                                            7331.   
C                                                                       7332.   
  130 CONTINUE                                                          7333.   
C                                                                       7334.   
      RETURN                                                            7335.   
      END                                                               7336.   
