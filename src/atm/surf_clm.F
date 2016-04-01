
#include "ctrparam.h"

!	==========================================================
!
!	SURFACE.F:  THIS SUBROUTINE CALCULATES THE SURFACE FLUXES 
!		    WHICH INCLUDE SENSIBLE HEAT, EVAPORATION, 
!		    THERMAL RADIATION, AND MOMENTUM DRAG.  IT ALSO 
!		    CALCULATES INSTANTANEOUSLY SURFACE TEMPERATURE,
!		    SURFACE SPECIFIC HUMIDITY, AND SURFACE WIND 
!		    COMPONENTS.
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
!	092301	Chien Wang	add bc and oc
!
!	==========================================================

      SUBROUTINE SURF_CLM

C****                                                                   5802.   
C**** THIS SUBROUTINE CALCULATES THE SURFACE FLUXES WHICH INCLUDE       5803.   
C**** SENSIBLE HEAT, EVAPORATION, THERMAL RADIATION, AND MOMENTUM       5804.   
C**** DRAG.  IT ALSO CALCULATES INSTANTANEOUSLY SURFACE TEMPERATURE,    5805.   
C**** SURFACE SPECIFIC HUMIDITY, AND SURFACE WIND COMPONENTS.           5806.   
C****                                                                   5807.   

#if ( defined CLM )
#if ( defined CPL_CHEM )
!
#include "chem_para"
#include "chem_com"
!
#endif

#include "BD2G04.COM"

#include "CLM.h"

      COMMON/SPEC2/KM,KINC,COEK,C3LAND(IO0,JM0),C3OICE(IO0,JM0)         5808.1  
     *  ,C3LICE(IO0,JM0),WMGE(IO0,JM0),TSSFC(1,JM0,4)                   5808.2  
      COMMON U,V,T,P,Q                                                  5809.   
      COMMON/WORK1/CONV(IM0,JM0,LM0),PK(IM0,JM0,LM0),PREC(IM0,JM0),
     &  TPREC(IM0,JM0),                                                 5810.   
     *  COSZ1(IO0,JM0)                                                  5811.   
      COMMON/WORK2/UT(IM0,JM0,LM0),VT(IM0,JM0,LM0),DU1(IO0,JM0),
     &  DV1(IO0,JM0),                                                   5812.   
     *  RA(8),ID(8),UMS(8)                                              5813.   
      COMMON/WORK3/E0(IO0,JM0,4),E1(IO0,JM0,4),EVAPOR(IO0,JM0,4),       5814.   
     *  TGRND(IO0,JM0,4)                                                5814.1  
      COMMON/RDATA/ROUGHL(IO0,JM0)                                      5815.   
      DIMENSION SINI(72),COSI(72)                                       5816.   
      LOGICAL POLE,PRNT,HPRNT
      common/conprn/HPRNT,JPR
!     common/TSUR/TSURFC(JM0,0:13),TSURFT(JM0),TSURFD(JM0),DTEMSR(JM0)
#include "TSRF.COM"
      common/SURRAD/TRSURF(JM0,4),SRSURF(JM0,4)
c     REAL*8 B,TGV,TKV,TSV0,TSV1,TSV                                    5818.   
      COMMON/CWMG/WMGEA(JM0),NWMGEA(JM0),CHAVER(JM0),DTAV(JM0),DQAV(JM0)
     & ,Z0AV(JM0),WSAV(JM0),WS0AV(JM0),TAUAV(JM0)
C
      COMMON/SURFLAND/ DUL1(JM0),DVL1(JM0),DT1L(JM0),DQ1L(JM0),
     & WSSL(JM0),T2ML(JM0),
     & TSSL(JM0),QSSL(JM0),USSL(JM0),VSSL(JM0),TAUSL(JM0),BLJ(JM0,50)
     & ,ELHTG(JM0),SHTG(JM0),TAUXG(JM0),TAUYG(JM0)
c
      DATA RVAP/461.5/                                                  5819.   
      DATA SHV/0./,SHW/4185./,SHI/2060./,RHOW/1000./,RHOI/916.6/,       5820.   
     *  ALAMI/2.1762/,STBO/.5672573E-7/,TF/273.16/,TFO/-1.56/           5821.   
      DATA Z1I/.1/,Z2LI/2.9/,Z1E/.1/,Z2E/4./,RHOS/91.66/,ALAMS/.35/     5822.   
      QSAT(TM,PR,QLH)=3.797915*EXP(QLH*(7.93252E-6-2.166847E-3/TM))/PR  5836.   
      DLQSDT(TM,QLH)=QLH*2.166847E-3/(TM*TM)
      DATA IFIRST/1/                                                    5838.   
      ROSNOW(X)=0.54*X/LOG(1.+0.54*X/275.)
      ALSNOW(X)=2.8E-6*X**2
C****                                                                   5839.   
C**** FDATA  2  LAND COVERAGE (1)                                       5840.   
C****        3  RATIO OF LAND ICE COVERAGE TO LAND COVERAGE (1)         5841.   
C****                                                                   5842.   
C**** ODATA  1  OCEAN TEMPERATURE (C)                                   5843.   
C****        2  RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)       5844.   
C****        3  OCEAN ICE AMOUNT OF SECOND LAYER (KG/M**2)              5845.   
C****                                                                   5846.   
C**** GDATA  1  OCEAN ICE SNOW AMOUNT (KG/M**2)                         5847.   
C****        2  EARTH SNOW AMOUNT (KG/M**2)                             5848.   
C****        3  OCEAN ICE TEMPERATURE OF FIRST LAYER (C)                5849.   
C****        4  EARTH TEMPERATURE OF FIRST LAYER (C)                    5850.   
C****        5  EARTH WATER OF FIRST LAYER (KG/M**2)                    5851.   
C****        6  EARTH ICE OF FIRST LAYER (KG/M**2)                      5852.   
C****        7  OCEAN ICE TEMPERATURE OF SECOND LAYER (C)               5853.   
C****        8  EARTH TEMPERATURE OF SECOND LAYER (C)                   5854.   
C****        9  EARTH WATER OF SECOND LAYER (KG/M**2)                   5855.   
C****       10  EARTH ICE OF SECOND LAYER (KG/M**2)                     5856.   
C****       12  LAND ICE SNOW AMOUNT (KG/M**2)                          5857.   
C****       13  LAND ICE TEMPERATURE OF FIRST LAYER (C)                 5858.   
C****       14  LAND ICE TEMPERATURE OF SECOND LAYER (C)                5859.   
C****                                                                   5860.   
C**** BLDATA 1  COMPOSITE SURFACE WIND MAGNITUDE (M/S)                  5861.   
C****        2  COMPOSITE SURFACE AIR TEMPERATURE (K)                   5862.   
C****        3  COMPOSITE SURFACE AIR SPECIFIC HUMIDITY (1)             5863.   
C****        4  LAYER TO WHICH DRY CONVECTION MIXES (1)                 5864.   
C****        5  FREE                                                    5865.   
C****        6  COMPOSITE SURFACE U WIND                                5866.   
C****        7  COMPOSITE SURFACE V WIND                                5867.   
C****        8  COMPOSITE SURFACE MOMENTUM TRANSFER (TAU)               5868.   
C****                                                                   5869.   
C**** VDATA  9  WATER FIELD CAPACITY OF FIRST LAYER (KG/M**2)           5870.   
C****       10  WATER FIELD CAPACITY OF SECOND LAYER (KG/M**2)          5871.   
C****                                                                   5872.   
C**** ROUGHL    LOG(ZGS/ROUGHNESS LENGTH) (LOGARITHM TO BASE 10)        5873.   
C**** ROUGHL will be ROUGHNESS LENGTH
C****                                                                   5874.   
c     print *,'surface TAU=',TAU
      NSTEPS=NSURF*NSTEP/NDYN                                           5875.   
      IF(IFIRST.NE.1) GO TO 50                                          5876.   
      print *,' SURFACE FOR LAND AFTER CLM'
      IFIRST=0                                                          5877.   
      COEFSN=100./ROSNOW(10.)
      COEFSN=1.
      DTSURF=NDYN*DT/NSURF                                              5884.   
      NGRNDZ=NGRND
      DTGRND=DTSURF/NGRNDZ                                              6143.   
      NCLMPERDAY=(24.*3600.)/DTSURF
      print *,' DTSURF=',DTSURF
      print *,' DTGRND=',DTGRND
      print *,' NCLMPERDAY=',NCLMPERDAY
      SHA=RGAS/KAPA                                                     5886.   
      RVX=0.                                                            5887.   
      ZS1CO=.5*DSIG(1)*RGAS/GRAV                                        5896.   
      P1000K=EXPBYK(1000.)                                              5897.   
      COEFS=GRAV/(100.*DSIG(1))                                         5898.   
      COEF1=(1.-SIG(2))/DSIGO(1)                                        5899.   
      COEF2=(SIG(1)-1.)/DSIGO(1)                                        5900.   
   50 CONTINUE
C**** ZERO OUT ENERGY AND EVAPORATION FOR GROUND AND INITIALIZE TGRND   5906.   
      DO 70 J=1,JM                                                      5907.   
      DO 70 I=1,IM                                                      5908.   
      TGRND(I,J,3)=GDATA(I,J,13)                                        5910.   
      TGRND(I,J,4)=GDATA(I,J,4)                                         5911.   
      DO 70 K=3,4                                                       5912.
      EVAPOR(I,J,K)=0.  
      E1(I,J,K)=0.                                                      5913.
      E0(I,J,K)=0.                                                      5913.
   70 CONTINUE
         IHOUR=1.5+TOFDAY                                               5914.   
C****                                                                   5915.   
C**** OUTSIDE LOOP OVER TIME STEPS, EXECUTED NSURF TIMES EVERY HOUR     5916.   
C****                                                                   5917.   
C**** ZERO OUT LAYER 1 WIND INCREMENTS                                  5922.   
      DO 60 J=1,JM                                                      5923.   
      DUL1(J)=0. 
   60 DVL1(J)=0.
C****                                                                   5927.   
C**** OUTSIDE LOOP OVER J AND I, EXECUTED ONCE FOR EACH GRID POINT      5928.   
C****                                                                   5929.   
!     print *,'From surf_clm TAU=',TAU,' HPRNT=',HPRNT,JPR
      DO 7000 J=1,JM                                                    5930.   
      if(PRNT)then
       if(ns.eq.1)then
       write(78,*) ,'    '
       write(78,*) ,'TAU=',TAU
       endif
        write(78,*),'NS=',ns
      endif
  100    CONTINUE
         BTRHDT=0.                                                      5958.   
         BSHDT=0.                                                       5961.   
         BEVHDT=0.                                                      5964.   
         BT2=0.                                                         5967.   
         BDT2=0.                                                         5967.   
         BTAUL=0.
         BTAUF=0.
      IMAX=IM                                                            5969.   
      DO 6000 I=1,IMAX                                                  5970.   
C****                                                                   5971.   
C**** DETERMINE SURFACE CONDITIONS                                      5972.   
C****                                                                   5973.   
      PLAND=FDATA(I,J,2)                                                5974.   
      if(PLAND.gt.0.0)then
      SP=P(I,J)                                                         5980.   
      PS=SP+PTOP                                                        5981.   
      PSK=EXPBYK(PS)                                                    5982.   
      P1=SIG(1)*SP+PTOP                                                 5983.   
      P1K=EXPBYK(P1)                                                    5984.   
C     surface fluxes from radiation
      TRHT0=TRSURF(J,2)
      SRHEAT=SRSURF(J,2)*COSZ1(I,J)
C     surface fluxes from radiation

      RMBYA=100.*SP*DSIG(1)/GRAV                                        6001.   
C**** ZERO OUT QUANTITIES TO BE SUMMED OVER SURFACE TYPES               6002.   
      TAUS=0.                                                           6008.   
      T2MS=0.
C****                                                                   6032.   
      SNOW=snwdclm(i,j)
      ACE1=h2oiclm(i,j,1)
      WTR1=h2olclm(i,j,1)
      IF(SNOW.GT.0.) THEN
        ELHX=LHS  
      ELSE
        PFROZN=ACE1/(WTR1+ACE1+1.E-20)
        ELHX=LHE+LHM*PFROZN          
      ENDIF
 
      SRHDT=SRHEAT*DTSURF
      TKV=THV1*PSK                                                      6137.   
      ZS1=ZS1CO*TKV*SP/PS                                               6138.   
      P1=SIG(1)*SP+PTOP                                                 6139.   
      SHDT=0.                                                           6144.   
      EVHDT=0.                                                          6145.   
      TRHDT=0.                                                          6146.   
      F1DT=0.                                                           6147.   
C**** LOOP OVER GROUND TIME STEPS                                       6148.   
C**** CALCULATE FLUXES OF SENSIBLE HEAT, LATENT HEAT, THERMAL           6478.   
C****   RADIATION, AND CONDUCTION HEAT (WATTS/M**2)                     6479.   
      SHEAT=shfclm(i,j)
      EVHEAT=lhfclm(i,j)
#if ( defined CLM35 )
      TG=(abs(lwuclm(i,j))/STBOLAND(J))**(1./4.)
#else
      TG=(abs(lwuclm(i,j))/STBO)**(1./4.)
#endif
      TG=(abs(lwuclm(i,j))/STBO)**(1./4.)

      TG1=TG-TF
c     print *,'From surf_clm TAU=',TAU,' J=',j
c     print *,LHS,LHE,LHM,ELHX
c     print *,shfclm(i,j),lhfclm(i,j)
c     TRHEAT=TRHT0-STBO*(TG*TG)*(TG*TG)
      TRHEAT=TRHT0+lwuclm(i,j)
      SHDT=DTSURF*SHEAT                                                 6505.
      EVHDT=DTSURF*EVHEAT                                               6506.
      TRHDT=DTSURF*TRHEAT
!      if(TAU.gt.3623.0.and.TAU.lt.3744.0) then
       if(HPRNT)then
         if(J.eq.jpr)then
          print *,'Form surf_clm J=',j,' TAU=',TAU
          print *,'TG1=',TG1,' lwuclm(i,j)=',lwuclm(i,j)
          print '(A5,3F9.2)',' TGL ',TAU,TOFDAY,TG1
          print *,'TRHT0=',TRHT0,' TRHEAT=',TRHEAT
!         print *,'Form surf_clm J=',j,' TAU=',TAU
           print *,'shfclm(i,j)=',shfclm(i,j),
     &     'lhfclm(i,j)=',lhfclm(i,j)
           print *,'tref2mclm(i,j)=',tref2mclm(i,j),
     &     'tgndclm(i,j)=',tgndclm(i,j)
           print *,'tauxclm(i,j)=',tauxclm(i,j),
     &     'tauyclm(i,j)=',tauyclm(i,j)
          print *,'TG1=',TG1,' lwuclm(i,j)=',lwuclm(i,j)
      print *,SHEAT,EVHEAT
      print *,TRHEAT,SRHEAT
!         write (329) ,TAU,TOFDAY,TG1,TRHT0,lwuclm(i,j),
!    &      shfclm(i,j),lhfclm(i,j)
         endif
         if(J.eq.-30)then
          write (330) ,TAU,TOFDAY,TG1,TRHT0,lwuclm(i,j),
     &      shfclm(i,j),lhfclm(i,j)
         endif
        endif

c     DQ1=EVHDT/(ELHX*RMBYA)
c     EVAP=-EVHDT/ELHX
      EVAP=vetclm(i,j)+sevclm(i,j)+cevclm(i,j)
!!    EVAP=vetclm(i,j)
      EVAPV=vetclm(i,j)
      EVAPS=sevclm(i,j)
      EVAPC=cevclm(i,j)
      EVAP=EVAP*DTSURF
      DQ1=-EVAP/RMBYA
c     print *,EVHDT,SHDT,DQ1,EVAP

      F0=SRHEAT+TRHEAT+SHEAT+EVHEAT                                     6487.   

      TAUL=tauxclm(i,j)
      TAUF=tauyclm(i,j)
      WR=SQRT(VSSL(J)**2+USSL(J)**2)/WSSL(J)
      TAUL=WR*TAUL
      TAUF=WR*TAUF
      DUL1(J)=DUL1(J)+PLAND*DTGRND*TAUL*COEFS/SP  
      DVL1(J)=DVL1(J)+PLAND*DTGRND*TAUF*COEFS/SP 
      TAUYG(J)=TAUL
      TAUXG(J)=TAUF
c     print *,tauxclm(i,j),tauyclm(i,j)
c     print *,PLAND,DTGRND,COEFS,SP
c     print *,DUL1(J),DVL1(J)


c     TH2M=tref2mclm(i,j)
c     t2md4tem(j)=t2md4tem(j)+TH2M
      T2M=tref2mclm(i,j)
      t2md4tem(j)=t2md4tem(j)+T2M
      nt2md4tem(j)=nt2md4tem(j)+1
      DT2M=DT2MLD(j)
      F0DT=CORSR*SRHDT+TRHDT+SHDT+EVHDT                                 6510. 
      if (J.eq.-23) then
      print *,' From surf_clm TAU=',TAU
      print *,j,i,' T2M=',t2m
      print *,tsoiclm(i,j,1),TG,tgndclm(i,j)
!     TG=(abs(lwuclm(i,j))/STBO)**(1./4.)
      print *,SHEAT,EVHEAT
      print *,TRHEAT,SRHEAT
      endif
c     print *,'From surface ',TAU,CORSR,SRHDT,TRHDT,SHDT,EVHDT
C**** ACCUMULATE SURFACE FLUXES AND PROGNOSTIC AND DIAGNOSTIC QUANTITIES6517.   
      do ITYPE=3,4
      E0(I,J,ITYPE)=E0(I,J,ITYPE)+F0DT                                  6518.   
      E1(I,J,ITYPE)=E1(I,J,ITYPE)+F1DT                                  6519.   
      EVAPOR(I,J,ITYPE)=EVAPOR(I,J,ITYPE)+EVAP                          6520.   
      TGRND(I,J,ITYPE)=TG1                                              6521.   
      enddo

      DTH1=-SHDT*PLAND/(SHA*RMBYA*P1K)
      DQQ1=-DQ1*PLAND


      TAUS=TAUS+SQRT(TAUL**2+TAUF**2)*PLAND
      T2MS=T2MS+T2M*PLAND
      BSHDT=BSHDT+SHDT*PLAND
      BEVHDT=BEVHDT+EVHDT*PLAND
      BTRHDT=BTRHDT+TRHDT*PLAND
c     BT2=BT2+(TH2M-TF)*PLAND
      BT2=BT2+(T2M-TF)*PLAND
      BDT2=BDT2+DT2M*PLAND
      BTAUL=BTAUL+TAUL*PLAND
      BTAUF=BTAUF+TAUF*PLAND

 5000 CONTINUE
      DT1L(J)=DTH1       
      DQ1L(J)=DQQ1       
      TAUSL(J)=TAUS
      T2ML(J)=T2MS
      TLANDD(J)=TLANDD(J)+(T2M-273.16)/NCLMPERDAY
C****                                                                   6596.   
C**** ACCUMULATE DIAGNOSTICS                                            6597.   
C****                                                                   6598.   
      endif
 6000 CONTINUE
C**** QUANTITIES ACCUMULATED FOR SURFACE TYPE TABLES IN DIAG1           6663.   
         BLJ(J,9)=BTRHDT        
         BLJ(J,13)=BSHDT      
         BLJ(J,14)=BEVHDT    
         BLJ(J,32)=BTAUL
         BLJ(J,33)=BTAUF
         BLJ(J,38)=BTAUL
         BLJ(J,26)=BT2 
         BLJ(J,27)=BDT2 
         BJ(J,41)=BJ(J,41)+EVAPV*PLAND
         BJ(J,42)=BJ(J,42)+EVAPS*PLAND
         BJ(J,16)=BJ(J,16)+EVAPC*PLAND
         if(J.eq.-23)then
          print *,'TAU=',TAU,' EVHEAT=',EVHEAT
          print *,'EVAP=',EVAP,' EVAP1=',-EVHDT/ELHX
         endif
 7000 CONTINUE                                                          6677.   
!     do J=1,JM
!      TLANDD(J)=TLANDD(J)+(T2ML(J)-273.16)/NCLMPERDAY
!     enddo
!     print *,' From surf_clm T2ML TAU=',TAU
!     print *,T2ML
!     print *,'TLANDD'
!     print *,TLANDD
!     print *,'NCLMPERDAY=',NCLMPERDAY
c      write(935),TAU,ELHTG,SHTG,TAUXG,TAUYG
C****                                                                   6678.   
#endif
      RETURN                                                            6795.   
  990 FORMAT ('0PPBL',3I4,14F8.2)                                       6818.   
  991 FORMAT ('0SURFACE ',4I4,5F10.4,3F11.7)                            6819.   
  992 FORMAT ('0',I2,10F10.4/23X,4F10.4,10X,2F10.4/                     6820.   
     *  33X,3F10.4,10X,2F10.4)                                          6821.   
  993 FORMAT ('0',I2,10F10.4/23X,7F10.4/33X,7F10.4)                     6822.   
  994 FORMAT ('0',I2,11F10.4)                                           6823.   
      END                                                               6824.   
