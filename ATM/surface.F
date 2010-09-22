C $Header: /home/igsm-repository/ATM/surface.F,v 1.1 2010/08/30 15:36:35 jscott Exp $
C $Name:  $

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

      SUBROUTINE SURFCE                                                 5801.   

C****                                                                   5802.   
C**** THIS SUBROUTINE CALCULATES THE SURFACE FLUXES WHICH INCLUDE       5803.   
C**** SENSIBLE HEAT, EVAPORATION, THERMAL RADIATION, AND MOMENTUM       5804.   
C**** DRAG.  IT ALSO CALCULATES INSTANTANEOUSLY SURFACE TEMPERATURE,    5805.   
C**** SURFACE SPECIFIC HUMIDITY, AND SURFACE WIND COMPONENTS.           5806.   
C****                                                                   5807.   

#if ( defined CPL_CHEM )
!
#include "chem_para"
#include "chem_com"
!
#endif

#include "BD2G04.COM"

#if ( defined CLM )
#include "CLM.h"
#endif

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
      DIMENSION WMGMINO(JM0)
      LOGICAL POLE,PRNT,HPRNT
      common/conprn/HPRNT          
!     common/TSUR/TSURFC(JM0,0:13),TSURFT(JM0),TSURFD(JM0),DTEMSR(JM0)
#include "TSRF.COM"
      common/SURRAD/TRSURF(JM0,4),SRSURF(JM0,4)
c     REAL*8 B,TGV,TKV,TSV0,TSV1,TSV                                    5818.   
      COMMON/CWMG/WMGEA(JM0),NWMGEA(JM0),CHAVER(JM0),DTAV(JM0),DQAV(JM0)
     & ,Z0AV(JM0),WSAV(JM0),WS0AV(JM0),TAUAV(JM0)
C
#if ( defined OCEAN_3D || defined ML_2D)
#include "AGRID.h"
#endif
c
      DATA RVAP/461.5/                                                  5819.   
      DATA SHV/0./,SHW/4185./,SHI/2060./,RHOW/1000./,RHOI/916.6/,       5820.   
     *  ALAMI/2.1762/,STBO/.5672573E-7/,TF/273.16/,TFO/-1.56/           5821.   
      DATA Z1I/.1/,Z2LI/2.9/,Z1E/.1/,Z2E/4./,RHOS/91.66/,ALAMS/.35/     5822.   
      DIMENSION AROUGH(20),BROUGH(20),CROUGH(20),DROUGH(20),EROUGH(20)  5823.   
      DATA  AROUGH/16.59,13.99,10.4,7.35,5.241,3.926,3.126,2.632,2.319, 5824.   
     *2.116,1.982,1.893,1.832,1.788,1.757,1.733,1.714,1.699,1.687,1.677/5825.   
      DATA BROUGH/3.245,1.733,0.8481,0.3899,0.1832,0.9026E-1,0.4622E-1, 5826.   
     * .241E-1,.1254E-1,.6414E-2,.3199E-2,.1549E-2,.7275E-3,.3319E-3,   5827.   
     * .1474E-3,.6392E-4,.2713E-4,.1130E-4,.4630E-5,.1868E-5/           5828.   
      DATA  CROUGH/5.111,3.088,1.682,.9239,.5626,.3994,.3282,.3017,.299 5829.   
     *,.3114,.3324,.3587,.3881,.4186,.4492,.4792,.5082,.5361,.5627,     5830.   
     * .5882/                                                           5831.   
      DATA DROUGH/1.24,1.02,0.806,0.682,0.661,0.771,0.797,0.895,0.994,  5832.   
     * 1.09,1.18,1.27,1.35,1.43,1.50,1.58,1.65,1.71,1.78,1.84/          5833.   
      DATA EROUGH/0.128,0.130,0.141,0.174,0.238,0.330,0.438,0.550,0.660,5834.   
     * 0.766,0.866,0.962,1.05,1.14,1.22,1.30,1.37,1.45,1.52,1.58/       5835.   
      QSAT(TM,PR,QLH)=3.797915*EXP(QLH*(7.93252E-6-2.166847E-3/TM))/PR  5836.   
      DLQSDT(TM,QLH)=QLH*2.166847E-3/(TM*TM)
c     TLOG(Z0)=ALOG(.36*SQRTT/(FMAG*Z0))+2.302585*ROUGH-.08             5837.   
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
C****        5   SURFACE MOMENTUM TRANSFER (TAU) OCEAN                  5865.   
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
      print *,' SURFACE CORSR=',CORSR
      print *,' ZGS=30 m for LAND '
      IFIRST=0                                                          5877.   
      WMGMINL = 5.0 
      print *,'WMGMIN 4 LAND=',WMGMINL
      print *,'over land  WMG=max(WMG0,WMGMIN)'
      WMGM0=8.0
      WMGM45=25.
      print *,'  WMGM0=', WMGM0,' WMGM45=',WMGM45
      WMGMAV=0.5*( WMGM0+WMGM45)
      DWGM=0.5*( WMGM0-WMGM45)
      do j = 1,jm0
        rhrad = 3.14159*(-90.+4.*(j-1))/180.
        WMGMINO(j) =  WMGMAV+DWGM*cos(4.*rhrad)
      enddo
      print *,'WMGMIN 4 OCEAN is a function of latitude'
      print 258,(WMGMINO(J),J=1,JM)
      print *,' WMGE'
        print 258,(WMGE(1,J),J=1,JM)
  258 format(12f5.1)
!     print *,'ODATA(1,7,2)=',ODATA(1,7,2)
      COEFSN=100./ROSNOW(10.)
      COEFSN=1.
      print *,' COEFSN=',COEFSN
      do 2567 J=1,JM
        NWMGEA(J)=0
        WMGEA(J)=0.
        CHAVER(J)=0.
        DTAV(J)=0.
        DQAV(J)=0.
        Z0AV(J)=0.
        WSAV(J)=0.
        WS0AV(J)=0. 
	TAUAV(J)=0.
 2567 CONTINUE
      READ (519) ((ROUGHL(I,J),I=1,IO),J=1,JM)                          5878.   
c     DO 10 J=2,JMM1                                                    5878.01 
C *************
      DO 10 J=1,JM
C *************
      ILAND=0.
      SUM1=0.                                                           5878.02 
      CONT1=0.                                                          5878.03 
      CONT2=0.
      DO 11 I=1,IO                                                      5878.04 
      PLAND=C3LAND(I,J)                                                 5878.05 
      CONT1=CONT1+PLAND                                                 5878.06 
       ROUGHL(I,J)=10**(log10(30.)-ROUGHL(I,J))
C**** ROUGHL IS NOW ROUGHNESS LENGTH 
 11   SUM1=SUM1+PLAND*ROUGHL(I,J)                                       5878.07 
      IF(CONT1.LE.0.) GO TO 10                                          5878.08 
      SUM1=SUM1/CONT1                                                   5878.09 
      DO 12 I=1,IO                                                      5878.1  
 12   ROUGHL(I,J)=SUM1                                                  5878.11 
 10   CONTINUE                                                          5878.12 
C     SRCORX=1.                                                         5878.13 
      CIAX=0.3
      print *,' surfacen  CIAX=',CIAX
      print *,' QS=Q1, TS=T1'
      print *,' WS=sqrt(0.75*W1+WGEM) '
	print *,' ROUGHL'
	print *,(ROUGHL(1,J),J=1,jm)
      REWIND 519                                                        5879.   
      LBLMM1=LBLM-1                                                     5880.   
      IQ1=IM/4+1                                                        5881.   
      IQ2=IM/2+1                                                        5882.   
      IQ3=3*IM/4+1                                                      5883.   
      DTSURF=NDYN*DT/NSURF                                              5884.   
      print *,' DTSURF=',DTSURF
         DTSRCE=DT*NDYN                                                 5885.   
      SHA=RGAS/KAPA                                                     5886.   
      RVX=0.                                                            5887.   
      ACE1I=Z1I*RHOI                                                    5888.   
      HC1I=ACE1I*SHI                                                    5889.   
      HC2LI=Z2LI*RHOI*SHI                                               5890.   
      HC1DE=Z1E*1129950.                                                5891.   
      HC2DE=Z2E*1129950.+3.5*.125*RHOW*3100.                            5892.   
      Z1IBYL=Z1I/ALAMI                                                  5893.   
      Z2LI3L=Z2LI/(3.*ALAMI)                                            5894.   
      BYRSL=1./(RHOS*ALAMS)                                             5895.   
      ZS1CO=.5*DSIG(1)*RGAS/GRAV                                        5896.   
      P1000K=EXPBYK(1000.)                                              5897.   
      COEFS=GRAV/(100.*DSIG(1))                                         5898.   
      COEF1=(1.-SIG(2))/DSIGO(1)                                        5899.   
      COEF2=(SIG(1)-1.)/DSIGO(1)                                        5900.   
      DO 20 I=1,IM                                                      5901.   
      SINI(I)=SIN((I-1)*TWOPI/FIM)                                      5902.   
   20 COSI(I)=COS((I-1)*TWOPI/FIM)                                      5903.   
   50 S0=S0X*1367./RSDIST                                               5904.   
      BYS0=RSDIST/1367.                                                 5905.   
C**** ZERO OUT ENERGY AND EVAPORATION FOR GROUND AND INITIALIZE TGRND   5906.   
      DO 70 J=1,JM                                                      5907.   
      DO 70 I=1,IM                                                      5908.   
      TGRND(I,J,2)=GDATA(I,J,3)                                         5909.   
      TGRND(I,J,3)=GDATA(I,J,13)                                        5910.   
      TGRND(I,J,4)=GDATA(I,J,4)                                         5911.   
      DO 70 K=1,12                                                      5912.   
   70 E0(I,J,K)=0.                                                      5913.   
         IHOUR=1.5+TOFDAY                                               5914.   
C****                                                                   5915.   
C**** OUTSIDE LOOP OVER TIME STEPS, EXECUTED NSURF TIMES EVERY HOUR     5916.   
C****                                                                   5917.   
      DO 9000 NS=1,NSURF                                                5918.   
         MODDSF=MOD(NSTEPS+NS-1,NDASF)                                  5919.   
         IF(MODDSF.EQ.0) IDACC(3)=IDACC(3)+1                            5920.   
         MODD6=MOD(IDAY+NS,NSURF)                                       5921.   
C**** ZERO OUT LAYER 1 WIND INCREMENTS                                  5922.   
      DO 60 J=1,JM                                                      5923.   
      DO 60 I=1,IM                                                      5924.   
      DU1(I,J)=0.                                                       5925.   
   60 DV1(I,J)=0.                                                       5926.   
C****                                                                   5927.   
C**** OUTSIDE LOOP OVER J AND I, EXECUTED ONCE FOR EACH GRID POINT      5928.   
C****                                                                   5929.   
      JPR=-7
      DO 7000 J=1,JM                                                    5930.   
      PRNT=j.eq.8
      PRNT=.FALSE.
      if(PRNT)then
       if(ns.eq.1)then
       write(78,*) ,'    '
       write(78,*) ,'TAU=',TAU
       endif
        write(78,*),'NS=',ns
      endif
      HEMI=1.                                                           5931.   
      IF(J.LE.JM/2) HEMI=-1.                                            5932.   
      FCOR=2.*OMEGA*SINP(J)                                             5933.   
      FMAG=FCOR*HEMI                                                    5934.   
      ROOT2F=SQRT(2.*FMAG)                                              5935.   
      IF(J.EQ.1) GO TO 80                                               5936.   
      IF(J.EQ.JM) GO TO 90                                              5937.   
      WMG0=.5*(WMGE(1,J)+WMGE(1,J+1))+.001                              5937.5  
      POLE=.FALSE.                                                      5938.   
      IMAX=IM                                                           5939.   
      GO TO 100                                                         5940.   
C**** CONDITIONS AT THE SOUTH POLE                                      5941.   
   80 POLE=.TRUE.                                                       5942.   
      IMAX=1                                                            5943.   
      JVPO=2                                                            5944.   
      RAPO=2.*RAPVN(1)                                                  5945.   
      U1=.25*(U(1,2,1)+V(IQ1,2,1)-U(IQ2,2,1)-V(IQ3,2,1))                5946.   
      V1=.25*(V(1,2,1)-U(IQ1,2,1)-V(IQ2,2,1)+U(IQ3,2,1))                5947.   
      WMG0=WMGE(1,2)+.001                                               5947.5  
      GO TO 100                                                         5948.   
C**** CONDITIONS AT THE NORTH POLE                                      5949.   
   90 POLE=.TRUE.                                                       5950.   
      IMAX=1                                                            5951.   
      JVPO=JM                                                           5952.   
      RAPO=2.*RAPVS(JM)                                                 5953.   
      U1=.25*(U(1,JM,1)-V(IQ1,JM,1)-U(IQ2,JM,1)+V(IQ3,JM,1))            5954.   
      V1=.25*(V(1,JM,1)+U(IQ1,JM,1)-V(IQ2,JM,1)-U(IQ3,JM,1))            5955.   
      WMG0=WMGE(1,JM)+.001                                              5955.5  
C**** ZERO OUT SURFACE DIAGNOSTICS WHICH WILL BE SUMMED OVER LONGITUDE  5956.   
  100    ATRHDT=0.                                                      5957.   
         BTRHDT=0.                                                      5958.   
         CTRHDT=0.                                                      5959.   
         ASHDT=0.                                                       5960.   
         BSHDT=0.                                                       5961.   
         CSHDT=0.                                                       5962.   
         AEVHDT=0.                                                      5963.   
         BEVHDT=0.                                                      5964.   
         CEVHDT=0.                                                      5965.   
         ATS=0.                                                         5966.   
         BTS=0.                                                         5967.   
         CTS=0.                                                         5968.   
         AT2=0.                                                         5966.   
         BT2=0.                                                         5967.   
         CT2=0.                                                         5968.   
         ATAUL=0.
         ATAUF=0.
         BTAUL=0.
         BTAUF=0.
         CTAUL=0.
         CTAUF=0.
	AWS=0.
	BWS=0.
	CWS=0.
        AWMG=0.
        BWMG=0.
        CWMG=0.
	ACH=0.
	BCH=0.
	CCH=0.
      IM1=IM                                                            5969.   
#if ( defined CLM )
         if(NS.eq.1)then
           i=1
           tsl4clm(i,j)=0.0
           qs4clm(i,j)=0.0
           ps4clm(i,j)=0.0
           ws4clm(i,j)=0.0
           us4clm(i,j)=0.0
           vs4clm(i,j)=0.0
         endif
#endif
      DO 6000 I=1,IMAX                                                  5970.   
C****                                                                   5971.   
C**** DETERMINE SURFACE CONDITIONS                                      5972.   
C****                                                                   5973.   
      PLAND=FDATA(I,J,2)                                                5974.   
      PWATER=1.-PLAND                                                   5975.   
      PLICE=FDATA(I,J,3)*PLAND                                          5976.   
      PEARTH=PLAND-PLICE                                                5977.   
      POICE=ODATA(I,J,2)*PWATER                                         5978.   
      POCEAN=PWATER-POICE                                               5979.   
      if(POCEAN.LE.1.E-5)then
         POCEAN=0.
         POICE=PWATER
      endif
      TTOFR=PEARTH+PLICE+POICE+POCEAN 
      if(abs(TTOFR-1).gt.1.e-3)then
        print *,' From surface TTOFR=',TTOFR
        print *,' J=',J,' PLAND=',PLAND,' POCEAN=',POCEAN
        print *,'POICE=',POICE,' ODATA(I,J,2)=',ODATA(I,J,2)
        stop
      end if
      SP=P(I,J)                                                         5980.   
      PS=SP+PTOP                                                        5981.   
      PSK=EXPBYK(PS)                                                    5982.   
      P1=SIG(1)*SP+PTOP                                                 5983.   
      P1K=EXPBYK(P1)                                                    5984.   
      WSOLD=BLDATA(I,J,1)                                               5985.   
      USOLD=BLDATA(I,J,6)                                               5986.   
      VSOLD=BLDATA(I,J,7)                                               5987.   
      TOLD=BLDATA(I,J,8)                                                5988.   
      SQRTT=SQRT(TOLD)                                                  5989.   
      GKBYFW=.1296*GRAV/(FCOR*FMAG*WSOLD+1.E-20)                        5990.   
      COSWS=GKBYFW*USOLD                                                5991.   
      SINWS=GKBYFW*VSOLD                                                5992.   
      IF(POLE) GO TO 1200                                               5993.   
      U1=.25*(U(IM1,J,1)+U(I,J,1)+U(IM1,J+1,1)+U(I,J+1,1))              5994.   
      V1=.25*(V(IM1,J,1)+V(I,J,1)+V(IM1,J+1,1)+V(I,J+1,1))              5995.   
      if(J.eq.JPR.or.J.eq.-12)then
        print *,' J=',J
        print *,'TAU=',TAU,' NS=',NS,' ITYPE=',ITYPE
       print *,'U(I,J,1)=',U(I,J,1),' V(I,J,1)=',V(I,J,1)
       print *,'U(I,J+1,1)=',U(I,J+1,1),' V(I,J+1,1)=',V(I,J+1,1)
       print *,'U(IM1,J,1)=',U(IM1,J,1),' V(IM1,J,1)=',V(IM1,J,1)
       print *,'U(IM1,J+1,1)=',U(IM1,J+1,1),
     &    ' V(IM1,J+1,1)=',V(IM1,J+1,1)
      endif
 1200 TH1=T(I,J,1)                                                      5996.   
      Q1=Q(I,J,1)                                                       5997.   
      DTH1=0.0
      DQQ1=0.0
      THV1=TH1*(1.+Q1*RVX)                                              5998.   
c     SRHEAT=SRHR(I,J,1)*COSZ1(I,J)*SRCOR                               5999.   
c     SRHDT=SRHEAT*DTSURF                                               6000.   
      RMBYA=100.*SP*DSIG(1)/GRAV                                        6001.   
C**** ZERO OUT QUANTITIES TO BE SUMMED OVER SURFACE TYPES               6002.   
      USS=0.                                                            6003.   
      VSS=0.                                                            6004.   
      WSS=0.                                                            6005.   
      TSS=0.                                                            6006.   
      QSS=0.                                                            6007.   
      TAUS=0.                                                           6008.   
         SINAPS=0.                                                      6009.   
         COSAPS=0.                                                      6010.   
         JR=J
         DXYPJ=DXYP(J)                                                  6012.   
         TG1S=0.                                                        6013.   
         QGS=0.                                                         6014.   
         BETAS=0.                                                       6015.   
         TRHDTS=0.                                                      6016.   
         SHDTS=0.                                                       6017.   
         EVHDTS=0.                                                      6018.   
         UGS=0.                                                         6019.   
         VGS=0.                                                         6020.   
         WGS=0.                                                         6021.   
         USRS=0.                                                        6022.   
         VSRS=0.                                                        6023.   
         RIS1S=0.                                                       6024.   
         RIGSS=0.                                                       6025.   
         CDMS=0.                                                        6026.   
         CDHS=0.                                                        6027.   
         DGSS=0.                                                        6028.   
         EDS1S=0.                                                       6029.   
         PPBLS=0.                                                       6030.   
         EVAPS=0.                                                       6031.   
C****                                                                   6032.   
      IF(POCEAN.LE.0.) GO TO 2200                                       6033.   
C****                                                                   6034.   
C**** OCEAN                                                             6035.   
C****                                                                   6036.   
      ITYPE=1                                                           6037.   
      PTYPE=POCEAN                                                      6038.   
c     formula charnoka
      TOCEAN=BLDATA(I,J,5)
      ROUGH=MAX(0.018*TOCEAN/GRAV,1.5E-4)
c     ROUGH=MAX(0.035*TOCEAN/GRAV,2.5E-4) used in 008.03
      ZGS=10.                                                           6041.   
!     WMGMIN=8.
      WMGMIN=WMGMINO(J)
      NGRNDZ=1                                                          6043.   
      TG1=ODATA(I,J,1)                                                  6044.   
      BETA=1.                                                           6045.   
      ELHX=LHE                                                          6046.   
      TRHT0=TRSURF(J,1)
      SRHEAT=SRSURF(J,1)*COSZ1(I,J)*SRCOR
      GO TO 3000                                                        6047.   
C****                                                                   6048.   
 2200 IF(POICE.LE.0.) GO TO 2400                                        6049.   
C****                                                                   6050.   
C**** OCEAN ICE                                                         6051.   
C****                                                                   6052.   
      ITYPE=2                                                           6053.   
      PTYPE=POICE                                                       6054.   
      NGRNDZ=NGRND                                                      6055.   
      SNOW=GDATA(I,J,1)                                                 6056.   
      TG1=TGRND(I,J,2)                                                  6057.   
      TG2=GDATA(I,J,7)                                                  6058.   
      ACE2=ODATA(I,J,3)                                                 6059.   
      Z2=ACE2/RHOI                                                      6060.   
      Z2BY4L=Z2/(4.*ALAMI)                                              6061.   
      if (SNOW.gt.10.)then
         RHOS0=ROSNOW(SNOW)
      else
         RHOS0=275.
      endif
      RHOS=COEFSN*RHOS0
      ALAMS=ALSNOW(RHOS0)
      BYRSL=1./(RHOS*ALAMS)
c     Z1BY6L=(Z1IBYL+SNOW*BYRSL)*.1666667                               6062.   
c     CDTERM=1.5*TG2-.5*TFO                                             6063.   
      CDTERM=TG2
c     CDENOM=1./(2.*Z1BY6L+Z2BY4L)                                      6064.   
      Z1BY2L=(Z1IBYL+SNOW*BYRSL)*0.5
      CDENOM=1./(Z1BY2L+2.*Z2BY4L)
      ROUGH=10**(log10(10.)-4.37)
      ZGS=10.                                                           6067.   
!     WMGMIN=8.
      WMGMIN=WMGMINO(J)
      HC1=HC1I+SNOW*SHI                                                 6069.   
      BETA=1.                                                           6070.   
      ELHX=LHS                                                          6071.   
      TRHT0=TRSURF(J,3)
      SRHEAT=SRSURF(J,3)*COSZ1(I,J)*SRCOR
      GO TO 3000                                                        6072.   
C****                                                                   6073.   
 2400 IF(PLAND.LE.0.) GO TO 5000                                        6074.   
      NGRNDZ=NGRND                                                      6075.   
      ROUGH=ROUGHL(I,J)                                                 6076.   
      ZGS=30.                                                           6078.   
      WMGMIN=WMGMINL
      TRHT0=TRSURF(J,2)
      SRHEAT=SRSURF(J,2)*COSZ1(I,J)*SRCOR
      IF(PLICE.LE.0.) GO TO 2600                                        6080.   
C****                                                                   6081.   
C**** LAND ICE                                                          6082.   
C****                                                                   6083.   
      ITYPE=3                                                           6084.   
      PTYPE=PLICE                                                       6085.   
      SNOW=GDATA(I,J,12)                                                6086.   
      TG1=TGRND(I,J,3)                                                  6087.   
      TG2=GDATA(I,J,14)                                                 6088.   
      if (SNOW.gt.10.)then
         RHOS0=ROSNOW(SNOW)
      else
         RHOS0=275.
      endif
      RHOS=COEFSN*RHOS0
      ALAMS=ALSNOW(RHOS0)
      BYRSL=1./(RHOS*ALAMS)
c     Z1BY6L=(Z1IBYL+SNOW*BYRSL)*.1666667                               6089.   
      CDTERM=TG2                                                        6090.   
c     CDENOM=1./(2.*Z1BY6L+Z2LI3L)                                      6091.   
      Z1BY2L=(Z1IBYL+SNOW*BYRSL)*0.5
      CDENOM=1./(Z1BY2L+3.*Z2LI3L/2.)
      HC1=HC1I+SNOW*SHI                                                 6092.   
      BETA=1.                                                           6093.   
      ELHX=LHS                                                          6094.   
      GO TO 3000                                                        6095.   
C****                                                                   6096.   
 2600 IF(PEARTH.LE.0.) GO TO 5000                                       6097.   
C****                                                                   6098.   
C**** EARTH                                                             6099.   
C****                                                                   6100.   
      ITYPE=4                                                           6101.   
      PTYPE=PEARTH                                                      6102.   
      SNOW=GDATA(I,J,2)                                                 6103.   
      TG1=TGRND(I,J,4)                                                  6104.   
      WTR1=GDATA(I,J,5)                                                 6105.   
      ACE1=GDATA(I,J,6)                                                 6106.   
      TG2=GDATA(I,J,8)                                                  6107.   
      WTR2=GDATA(I,J,9)                                                 6108.   
      ACE2=GDATA(I,J,10)                                                6109.   
      WFC1=VDATA(I,J,9)                                                 6110.   
      WFC2=VDATA(I,J,10)                                                6111.   
      WTR1DRY=0.025*WFC1
      HC1=HC1DE+WTR1*SHW+(ACE1+SNOW)*SHI                                6112.   
      ALAM1D=2.+.5*(1.+2.*WTR1/WFC1)                                    6113.   
      ALAM2D=4.                                                         6114.   
      RMULCH=1.                                                         6115.   
      IF((SINP(J).GT..5).AND.(JDAY-91)*(273-JDAY).LT.0) RMULCH=.25      6116.   
      IF((SINP(J).LT.-.5).AND.(JDAY-91)*(273-JDAY).GE.0) RMULCH=.25     6117.   
      ALAM1V=RMULCH*(.4185+1.2555*WTR1/WFC1+ALAMI*ACE1/(Z1E*RHOI))      6118.   
      ALAM3V=.8370                                                      6119.   
      IF(TG2.LT.0.) ALAM3V=.4185+ALAMI*.15                              6120.   
      ALAM2V=.125*(.4185+1.2555*WTR2/WFC2+ALAMI*ACE2/(5.*Z1E*RHOI))     6121.   
     *  +.875*ALAM3V                                                    6122.   
      ALAM1E=VDATA(I,J,1)*ALAM1D+(1.-VDATA(I,J,1))*ALAM1V               6123.   
      ALAM2E=VDATA(I,J,1)*ALAM2D+(1.-VDATA(I,J,1))*ALAM2V               6124.   
      if (SNOW.gt.10.)then
         RHOS0=ROSNOW(SNOW)
      else
         RHOS0=275.
      endif
      RHOS=COEFSN*RHOS0
      ALAMS=ALSNOW(RHOS0)
      BYRSL=1./(RHOS*ALAMS)
c     Z1BY6L=(Z1E/ALAM1E+SNOW*BYRSL)*.1666667                           6125.   
      Z1BY2L=(Z1E/ALAM1E+SNOW*BYRSL)*0.5
      CDTERM=TG2                                                        6126.   
c     CDENOM=1./(2.*Z1BY6L+Z2E/(3.*ALAM2E))                             6127.   
      CDENOM=1./(Z1BY2L+Z2E/(2.*ALAM2E))  
      BETA=1.                                                           6128.   
      ELHX=LHS                                                          6129.   
      IF(SNOW.GT.0.) GO TO 3000                                         6130.   
      BETA=(WTR1+ACE1)/WFC1                                             6131.   
      BETA=max(((WTR1+ACE1-WTR1DRY)/WFC1),0.0)
      PFROZN=ACE1/(WTR1+ACE1+1.E-20)                                    6132.   
      ELHX=LHE+LHM*PFROZN                                               6133.   
       HC2E=HC2DE+WTR2*SHW+ACE2*SHI
C****                                                                   6134.   
C**** BOUNDARY LAYER INTERACTION                                        6135.   
C****                                                                   6136.   
 3000 continue
      SRHDT=SRHEAT*DTSURF
      TKV=THV1*PSK                                                      6137.   
      ZS1=ZS1CO*TKV*SP/PS                                               6138.   
      P1=SIG(1)*SP+PTOP                                                 6139.   
      DTGRND=DTSURF/NGRNDZ                                              6143.   
      SHDT=0.                                                           6144.   
      EVHDT=0.                                                          6145.   
      TRHDT=0.                                                          6146.   
      F1DT=0.                                                           6147.   
C**** LOOP OVER GROUND TIME STEPS                                       6148.   
      DO 3600 NG=1,NGRNDZ                                               6149.   
      TG=TG1+TF                                                         6150.   
      QG=QSAT(TG,PS,ELHX)                                               6151.   
      TGV=TG*(1.+QG*RVX)                                                6152.   
      W1=SQRT(U1*U1+V1*V1)
      WS0=W1
c     WS=SQRT(W1*W1+0.8*WMG)
!     WMG=WMG0+WMGMIN
!     07/17/2006
      if(ITYPE.le.2)then
        WMG=WMG0+WMGMIN
      else
        WMG=max(WMG0,WMGMIN)
      endif
!     07/17/2006
      WS=SQRT((0.75*W1)**2+WMG)
      if(J.eq.JPR)then
        print *,' '
        print *,'TAU=',TAU,' NS=',NS,' ITYPE=',ITYPE
        print *,'TG=',TG,' QG=',QG
        print *,'RVX=',RVX,' TG1=',TG1
      endif
!      if(TAU.gt.3623.0.and.TAU.lt.3744.0) then
!        if(ITYPE.eq.4)then
!         if(J.eq.29)then
!           write (329) ,TAU,TOFDAY,TG1
!         endif
!         endif
!       endif

#if ( defined CPL_OCEANCO2  || defined OCEAN_3D)
      if(ITYPE.eq.1)then
      NWMGEA(J)=NWMGEA(J)+1
      WSAV(J)=WSAV(J)+WS
      end if
#endif
 
      WG=WS
      THS=TH1
      QS=Q1
      TSV=THS*PSK
      Z0=ROUGH
      ROUGH=log10(ZGS/ROUGH)
      CDN=.0231/(ROUGH*ROUGH)
c     if(ITYPE.eq.1)then 
c         CDN=.00075+.000067*WSOLD
c         ROUGH=7.126-1.068*LOG(WSOLD+1.E-12)
c     endif
      LR=ROUGH*2.-.5  
      IF(LR.GT.20) LR=20
      IF(LR.LT.1) LR=1
      RIGS=ZGS*GRAV*(TSV-TGV)/(TGV*WS*WS)
      SINAP=0.
      COSAP=1.
      IF(RIGS.LE.0) THEN
C     surface layer has unstable stratification
         CIA=TWOPI*0.0625/(1.+WS*CIAX)
         DM=SQRT((1.-AROUGH(LR)*RIGS)*(1.-BROUGH(LR)*RIGS)/ 
     *     (1.-CROUGH(LR)*RIGS)) 
         DH=1.35*SQRT((1.-DROUGH(LR)*RIGS)/(1.-EROUGH(LR)*RIGS)) 
      ELSE
C     surface layer has stable stratification
         CIA=TWOPI*(0.09375-0.03125/(1.+4*RIGS**2))/(1.+WS*CIAX)
         DM=1./(1.+(11.238+89.9*RIGS)*RIGS)
         DH=1.35/(1.+1.93*RIGS)
      END IF
      CDH=CDN*DM*DH 
      if(J.eq.JPR)then
        print *,'TAU=',TAU,' NS=',NS,' ITYPE=',ITYPE
        print *,'WS=',WS,' ZGS=',ZGS
        print *,'DM=',DM,' DH=',DH
        print *,'RIGS=',RIGS,' TGV=',TGV
      endif
      USR=COS(CIA)
      VSR=SIN(CIA)*HEMI
      UG=U1
      VG=V1
      US=(USR*UG-VSR*VG)
      VS=(VSR*UG+USR*VG)
      RCDHWS=CDH*WS*100.*PS/(RGAS*TSV)
      if(J.eq.JPR)then
c       print *,' '
        print *,'TAU=',TAU,' NS=',NS,' ITYPE=',ITYPE
        print *,'CDH=',CDH,' RGAS=',RGAS
        print *,'PS=',PS,' TSV=',TSV
        print *,'WS=',WS,' RCDHWS=',RCDHWS
      endif
      TS=TSV/(1.+QS*RVX)                                                6467.   
      QSATS=QSAT(TS,PS,ELHX)                                            6468.   
c dLQS/dTs
      DLQSDTS=DLQSDT(TS,ELHX)
c dLQS/dTs
      IF(QS.LE.QSATS) GO TO 3500                                        6469.   
      DQSSDT=QSATS*ELHX/(RVAP*TS*TS)                                    6470.   
      X=(QS-QSATS)/(DQSSDT+(SHA/ELHX))                                  6471.   
      TS=TS+X                                                           6472.   
      QS=QS+X*(SHA/ELHX)                                                6473.   
C**** CALCULATE RHOS*CDM*WS AND RHOS*CDH*WS                             6474.   
 3500 CDM=CDN*DM                                                        6475.   
      RCDMWS=CDM*WS*100.*PS/(RGAS*TS)                                   6476.   
C**** CALCULATE FLUXES OF SENSIBLE HEAT, LATENT HEAT, THERMAL           6478.   
C****   RADIATION, AND CONDUCTION HEAT (WATTS/M**2)                     6479.   
      SHEAT=SHA*RCDHWS*(TS-TG)                                          6480.   
      BETAUP=BETA                                                       6481.   
      IF(QS.GT.QG) BETAUP=1.                                            6482.   
      EVHEAT=(LHE+TG1*SHV)*BETAUP*RCDHWS*(QS-QG)                        6483.   
c     TRHEAT=TRHR(I,J,1)-STBO*(TG*TG)*(TG*TG)                           6484.   
      TRHEAT=TRHT0-STBO*(TG*TG)*(TG*TG)
#if ( defined CLM )
         if(NS.eq.1)then
          if(ITYPE.EQ.4.or.ITYPE.EQ.3)then
           tsl4clm(i,j)=tsl4clm(i,j)+TS*PTYPE/PLAND
           qs4clm(i,j)=qs4clm(i,j)+QS*PTYPE/PLAND
           ps4clm(i,j)=ps4clm(i,j)+PS*PTYPE/PLAND
           ws4clm(i,j)=ws4clm(i,j)+WS*PTYPE/PLAND
           us4clm(i,j)=us4clm(i,j)+US*PTYPE/PLAND
           vs4clm(i,j)=vs4clm(i,j)+VS*PTYPE/PLAND
         endif
        endif
#endif
      if(J.eq.JPR)then
c       print *,' '
        print *,'TAU=',TAU,' NS=',NS,' ITYPE=',ITYPE
        print *,'TRHT0=',TRHT0,' STBO=',STBO
        print *,'TG=',TG,' TS=',TS
        print *,'TRHEAT=',TRHEAT
        print *,'SHA=',SHA,' RCDHWS=',RCDHWS
        print *,'SHEAT=',SHEAT
      endif
      TG1OLD=TG1
      SHEATOLD=SHEAT
#if ( defined OCEAN_3D )
      IF(ITYPE.EQ.1 .or. ITYPE.EQ.2) GO TO 3620  
#else
      IF(ITYPE.EQ.1) GO TO 3620                                         6485.   
#endif
C**** CALCULATE FLUXES USING IMPLICIT TIME STEP FOR NON-OCEAN POINTS    6486.   
      F0=SRHEAT+TRHEAT+SHEAT+EVHEAT                                     6487.   
c     F1=(TG1-CDTERM-F0*Z1BY6L)*CDENOM                                  6488.   
      F1=(TG1-CDTERM)*CDENOM
      DSHDTG=-RCDHWS*SHA 
      DQGDTG=QG*ELHX/(RVAP*TG*TG)                                       6490.   
      DEVDTG=-RCDHWS*LHE*BETAUP*DQGDTG   
      DTRDTG=-4.*STBO*TG*TG*TG                                          6492.   
      DF0DTG=DSHDTG+DEVDTG+DTRDTG                                       6493.   
c     DFDTG=DF0DTG-(1.-DF0DTG*Z1BY6L)*CDENOM                            6493.5  
      DFDTG=DF0DTG-CDENOM
c     DF1DTG=(1.-DF0DTG*Z1BY6L)*CDENOM
      DF1DTG=CDENOM
      DTG=(F0-F1)*DTGRND/(HC1-DTGRND*DFDTG)                             6494.   
      SHDT=SHDT+DTGRND*(SHEAT+DTG*DSHDTG)                               6495.   
      EVHDT=EVHDT+DTGRND*(EVHEAT+DTG*DEVDTG)                            6496.   
      TRHDT=TRHDT+DTGRND*(TRHEAT+DTG*DTRDTG)                            6497.   
      TG1=TG1+DTG
c     F1DT=F1DT+DTGRND*(TG1-CDTERM-(F0+DTG*DF0DTG)*Z1BY6L)*CDENOM       6498.   
      F1DT=F1DT+DTGRND*(TG1-CDTERM)*CDENOM
      DU1(I,J)=DU1(I,J)+PTYPE*DTGRND*RCDMWS*US*COEFS/SP                 6499.   
      DV1(I,J)=DV1(I,J)+PTYPE*DTGRND*RCDMWS*VS*COEFS/SP                 6500.   
c     TG1=TG1+DTG                                                       6501.   
 3600 CONTINUE                                                          6502.   
      GO TO 3700                                                        6503.   
C**** CALCULATE FLUXES USING EXPLICIT TIME STEP FOR OCEAN POINTS        6504.   
 3620 SHDT=DTSURF*SHEAT                                                 6505.   
      EVHDT=DTSURF*EVHEAT                                               6506.   
      TRHDT=DTSURF*TRHEAT                                               6507.   
      DU1(I,J)=DU1(I,J)+PTYPE*DTSURF*RCDMWS*US*COEFS/SP                 6508.   
      DV1(I,J)=DV1(I,J)+PTYPE*DTSURF*RCDMWS*VS*COEFS/SP                 6509.   
 3700 CONTINUE
      EPS=1.D-8
c     print *,'FROM SURFACE NS=',NS
c     print *,'J=',J,' ITYPE=',ITYPE
c     print *,RCDMWS,WS
      WWS=max(W1,1.D-4)
c     RO=SP*100/(RGAS*TG)
c     print *,'RO=',RO
c     USTAR=SQRT(RCDMWS*WS/RO)
c     TSTAR=SHEATOLD/(0.35*1007.*RO*USTAR)
c     ALPHAH=DH
c     TT2M=TG+TSTAR/ALPHAH*LOG(2.0/Z0)
c     TT2M=TG+TSTAR/ALPHAH*LOG(ZGS/Z0)
c       print *,'RIGS=',RIGS,' Z0=',Z0
c       print *,'CDN=',CDN
c       print *,'H=',SHDT/DTSURF,' TGM=',RCDMWS*WS
c       print *,' SHEATOLD=',SHEATOLD
c     print *,' USTAR=',USTAR,' TSTAR=',TSTAR
c     print *,' ALPHAH=',ALPHAH,' TT2M=',TT2M
c     print *,' TT2M=',TT2M
      ZTEM=ZGS
      ZTEM=2.0
c     print *,'ZTEM=',ZTEM
      CALL TZM(T2M,TH2M,ZTEM,Z0,ZGS,SP,TG,TS,RIGS,WS,
     &     -SHEATOLD,RCDMWS*WS,LR,EPS)
c     print *,'FROM SURFACE'
c     print *,'TS=',TS,' TG=',TG
c     print *,' T2M=',T2M,' TH2M=',TH2M
      F0DT=CORSR*SRHDT+TRHDT+SHDT+EVHDT                                 6510. 
      if(J.eq.JPR)then
        print *,'TAU=',TAU,' NS=',NS,' ITYPE=',ITYPE
        print *,'DTSURF=',DTSURF,' CORSR=',CORSR
        print *,'SRHDT=',SRHDT,' TRHDT=',TRHDT
        print *,'SHDT=',SHDT,' EVHDT=',EVHDT
        print *,'F0DT=',F0DT
        print *,'US=',US,' VS=',VS
        print *,'COEFS=',COEFS,' SP=',SP
      endif
c     print *,'From surface ',TAU,CORSR,SRHDT,TRHDT,SHDT,EVHDT
C**** CALCULATE EVAPORATION                                             6511.   
CCC   DQ1=EVHDT/((LHE+TG1*SHV)*RMBYA)                                   6512.   
      DQ1=EVHDT/(ELHX*RMBYA)
      IF(DQ1*PTYPE.LE.Q1) GO TO 3720                                    6513.   
      DQ1=Q1/PTYPE                                                      6514.   
CCC   EVHDT=DQ1*(LHE+TG1*SHV)*RMBYA                                     6515.   
      EVHDT=DQ1*ELHX*RMBYA
 3720 EVAP=-DQ1*RMBYA                                                   6516.   
C**** ACCUMULATE SURFACE FLUXES AND PROGNOSTIC AND DIAGNOSTIC QUANTITIES6517.   
      E0(I,J,ITYPE)=E0(I,J,ITYPE)+F0DT                                  6518.   
      E1(I,J,ITYPE)=E1(I,J,ITYPE)+F1DT                                  6519.   
      EVAPOR(I,J,ITYPE)=EVAPOR(I,J,ITYPE)+EVAP                          6520.   
      if(PRNT)then
c      write(78,*) ,'    '
c      write(78,*) ,'TAU=',TAU
       write(78,*) ,'J=',j,' ITYPE=',ITYPE,' PTYPE=',PTYPE
       write(78,*) ,'TS=',TS,' TG=',TG,' QS=',QS
       write(78,*) ,'TG1=',TG1,' TG1OLD=',TG1OLD
       write(78,*) ,'TG2=',TG2
       write(78,*) ,'SHEAT=',SHEAT,' EVHEAT=',EVHEAT
       write(78,*) ,'TRHEAT=',TRHEAT,' SRHEAT=',SRHEAT
       write(78,*) ,'EVAP mm/day=',24.*3600.*EVAP/DTSURF
       write(78,*) ,'EVAP=',EVAP,
     &  ' F0DT=',F0DT/DTSURF,' F1DT=',F1DT/DTSURF
      endif
#if ( defined OCEAN_3D || defined ML_2D )
C       For ocean model
c       if(NS.eq.2)then
#if ( defined ML_2D )
        if(ITYPE.eq.1)then
#endif
C        DNetHeat by DTG
          DSHDTG=-RCDHWS*SHA
          DQGDTG=QG*ELHX/(RVAP*TG*TG)
          DEVDTG=-RCDHWS*LHE*BETAUP*DQGDTG
          DTRDTG=-4.*STBO*TG*TG*TG
          DF0DTG=DSHDTG+DEVDTG+DTRDTG
          if(EVHEAT.lt.0.0)then
           DEVDTGEQ=EVHEAT*DLQSDTS
          else
           DEVDTGEQ=0.0
          endif
C        DNetHeat by DTG
#if ( defined OCEAN_3D )
        if(ITYPE.eq.1)then
#endif
          dhfodtg(j)=dhfodtg(j)+DF0DTG
          devodtg(j)=devodtg(j)-DEVDTG/LHE
          dhfodtgeq(j)=dhfodtgeq(j)+DEVDTGEQ
          devodtgeq(j)=devodtgeq(j)-DEVDTGEQ/LHE
           evao(j)=evao(j)+EVAP
           hfluxo(j)=hfluxo(j)+F0DT
           naveo(j)=naveo(j)+1
        endif
        if(ITYPE.eq.2)then
           evai(j)=evai(j)+EVAP
           hfluxi(j)=hfluxi(j)+F0DT
           dhfidtg(j)=dhfidtg(j)+DF0DTG
           devidtg(j)=devidtg(j)-DEVDTG/LHE
           dhfidtgeq(j)=dhfidtgeq(j)+DEVDTGEQ
           devidtgeq(j)=devidtgeq(j)-DEVDTGEQ/LHE
c          tairi(j)=tairi(j)+TS
           navei(j)=navei(j)+1
        endif
c       endif 	! NS
         tauu(j)=tauu(j)+RCDMWS*US*PTYPE
         tauv(j)=tauv(j)+RCDMWS*VS*PTYPE
C       For ocean model
#endif
      TGRND(I,J,ITYPE)=TG1                                              6521.   
       TSSFC(I,J,ITYPE)=TS                                              6521.5  

c     TH1=TH1-SHDT*PTYPE/(SHA*RMBYA*P1K)                                6522.   
c     Q1=Q1-DQ1*PTYPE                                                   6523.   
      
      DTH1=DTH1-SHDT*PTYPE/(SHA*RMBYA*P1K)    
      DQQ1=DQQ1-DQ1*PTYPE                    

      USS=USS+US*PTYPE                                                  6524.   
      VSS=VSS+VS*PTYPE                                                  6525.   
      WSS=WSS+WS*PTYPE                                                  6526.   
      TSS=TSS+TS*PTYPE                                                  6527.   
      QSS=QSS+QS*PTYPE                                                  6528.   
      TAUS=TAUS+CDM*WS*W1*PTYPE                                         6529.   
         SINAPS=SINAPS+SINAP*PTYPE                                      6530.   
         COSAPS=COSAPS+COSAP*PTYPE                                      6531.   
         TG1S=TG1S+TG1*PTYPE                                            6532.   
         QGS=QGS+QG*PTYPE                                               6533.   
         BETAS=BETAS+BETA*PTYPE                                         6534.   
         TRHDTS=TRHDTS+TRHDT*PTYPE                                      6535.   
         SHDTS=SHDTS+SHDT*PTYPE                                         6536.   
         EVHDTS=EVHDTS+EVHDT*PTYPE                                      6537.   
         UGS=UGS+UG*PTYPE                                               6538.   
         VGS=VGS+VG*PTYPE                                               6539.   
         WGS=WGS+WG*PTYPE                                               6540.   
         USRS=USRS+USR*PTYPE                                            6541.   
         VSRS=VSRS+VSR*PTYPE                                            6542.   
c        RIS1S=RIS1S+RIS1*PTYPE                                         6543.   
         RIGSS=RIGSS+RIGS*PTYPE                                         6544.   
         CDMS=CDMS+CDM*PTYPE                                            6545.   
         CDHS=CDHS+CDH*PTYPE                                            6546.   
c        DGSS=DGSS+DGS*PTYPE                                            6547.   
c        EDS1S=EDS1S+EDS1*PTYPE                                         6548.   
c        PPBLS=PPBLS+PPBL*PTYPE                                         6549.   
         EVAPS=EVAPS+EVAP*PTYPE                                         6550.   
      GO TO (4000,4100,4400,4600),ITYPE                                 6551.   
C****                                                                   6552.   
C**** OCEAN                                                             6553.   
C****                                                                   6554.   
 4000    ASHDT=ASHDT+SHDT*POCEAN                                        6555.   
         AEVHDT=AEVHDT+EVHDT*POCEAN                                     6556.   
         ATRHDT=ATRHDT+TRHDT*POCEAN                                     6557.   
         ATS=ATS+(TS-TF)*POCEAN                                         6558.   
         AT2=AT2+(TH2M-TF)*POCEAN 
         BLDATA(I,J,5)=CDM*WS*W1
         ATAUL=ATAUL+RCDMWS*US*POCEAN
         ATAUF=ATAUF+RCDMWS*VS*POCEAN
         AWS=AWS+WS*POCEAN
         AWMG=AWMG+SQRT(WMG)*POCEAN
         ACH=ACH+RCDHWS*POCEAN
      GO TO 2200                                                        6559.   
C****                                                                   6560.   
C**** OCEAN ICE                                                         6561.   
C****                                                                   6562.   
 4100    CSHDT=CSHDT+SHDT*POICE                                         6563.   
         CEVHDT=CEVHDT+EVHDT*POICE                                      6564.   
         CTRHDT=CTRHDT+TRHDT*POICE                                      6565.   
         CTS=CTS+(TS-TF)*POICE                                          6566.   
         CT2=CT2+(TH2M-TF)*POICE                                          6566.   
         CTAUL=CTAUL+RCDMWS*US*POICE
         CTAUF=CTAUF+RCDMWS*VS*POICE
         CWS=CWS+WS*POICE 
         CWMG=CWMG+SQRT(WMG)*POICE
         CCH=CCH+RCDHWS*POICE 
      GO TO 2400                                                        6567.   
C****                                                                   6568.   
C**** LAND ICE                                                          6569.   
C****                                                                   6570.   
 4400    BSHDT=BSHDT+SHDT*PLICE                                         6571.   
         BEVHDT=BEVHDT+EVHDT*PLICE                                      6572.   
         BTRHDT=BTRHDT+TRHDT*PLICE                                      6573.   
         BTS=BTS+(TS-TF)*PLICE                                          6574.   
         BT2=BT2+(TH2M-TF)*PLICE  
         BTAUL=BTAUL+RCDMWS*US*PLICE
         BTAUF=BTAUF+RCDMWS*VS*PLICE
         BWS=BWS+WS*PLICE
         BWMG=BWMG+SQRT(WMG)*PLICE
         BCH=BCH+RCDHWS*PLICE
      GO TO 2600                                                        6575.   
C****                                                                   6576.   
C**** EARTH                                                             6577.   
C****                                                                   6578.   
 4600    BSHDT=BSHDT+SHDT*PEARTH                                        6579.   
         BEVHDT=BEVHDT+EVHDT*PEARTH                                     6580.   
         BTRHDT=BTRHDT+TRHDT*PEARTH                                     6581.   
         BTS=BTS+(TS-TF)*PEARTH                                         6582.   
         BT2=BT2+(TH2M-TF)*PEARTH
         BTAUL=BTAUL+RCDMWS*US*PEARTH
         BTAUF=BTAUF+RCDMWS*VS*PEARTH
         BWS=BWS+WS*PEARTH
         BWMG=BWMG+SQRT(WMG)*PEARTH
         BCH=BCH+RCDHWS*PEARTH
C**** NON-OCEAN POINTS WHICH ARE NOT MELTING OR FREEZING WATER USE      6583.   
C****   IMPLICIT TIME STEPS                                             6584.   
C****                                                                   6585.   
C**** UPDATE SURFACE AND FIRST LAYER QUANTITIES                         6586.   
C****                                                                   6587.   
 5000 CONTINUE
      T(I,J,1)=TH1                                                      6588.   
     & +DTH1
      Q(I,J,1)=Q1                                                       6589.   
     & +DQQ1
      BLDATA(I,J,1)=WSS                                                 6590.   
      BLDATA(I,J,2)=TSS                                                 6591.   
      BLDATA(I,J,3)=QSS                                                 6592.   
      BLDATA(I,J,6)=USS                                                 6593.   
      BLDATA(I,J,7)=VSS                                                 6594.   
      BLDATA(I,J,8)=TAUS                                                6595.   
c     print *,j,T(I,J,1),Q(I,J,1)
c     print *,(TGRND(I,J,k),k=1,4)
c     print *,(EVAPOR(I,J,k),k=1,4)
c     print *,(E0(I,J,k),k=1,4)
c     print *,(E1(I,J,k),k=1,4)
c     print *,j,DU1(1,j),DV1(1,j)
C****                                                                   6596.   
C**** ACCUMULATE DIAGNOSTICS                                            6597.   
C****                                                                   6598.   
C**** QUANTITIES ACCUMULATED FOR REGIONS IN DIAG1                       6599.   
         IF(JR.EQ.JM) GO TO 5700                                        6600.   
         DJ(JR,9)=DJ(JR,9)+TRHDTS*DXYPJ                                 6601.   
         DJ(JR,13)=DJ(JR,13)+SHDTS*DXYPJ                                6602.   
         DJ(JR,14)=DJ(JR,14)+EVHDTS*DXYPJ                               6603.   
         DJ(JR,19)=DJ(JR,19)+EVAPS*DXYPJ                                6604.   
         IF(MODDSF.NE.0) GO TO 5700                                     6605.   
         DJ(JR,23)=DJ(JR,23)+(TSS-TF)*DXYPJ                             6606.   
 5700    CONTINUE
 6000 IM1=I                                                             6662.   
C**** QUANTITIES ACCUMULATED FOR SURFACE TYPE TABLES IN DIAG1           6663.   
         AJ(J,9)=AJ(J,9)+ATRHDT                                         6664.   
         BJ(J,9)=BJ(J,9)+BTRHDT                                         6665.   
         CJ(J,9)=CJ(J,9)+CTRHDT                                         6666.   
         AJ(J,13)=AJ(J,13)+ASHDT                                        6667.   
         BJ(J,13)=BJ(J,13)+BSHDT                                        6668.   
         CJ(J,13)=CJ(J,13)+CSHDT                                        6669.   
         AJ(J,14)=AJ(J,14)+AEVHDT                                       6670.   
         BJ(J,14)=BJ(J,14)+BEVHDT                                       6671.   
         CJ(J,14)=CJ(J,14)+CEVHDT                                       6672.   
         AJ(J,32)=AJ(J,32)+ATAUL
         BJ(J,32)=BJ(J,32)+BTAUL
         CJ(J,32)=CJ(J,32)+CTAUL
         AJ(J,33)=AJ(J,33)+ATAUF
         BJ(J,33)=BJ(J,33)+BTAUF
         CJ(J,33)=CJ(J,33)+CTAUF
         AJ(J,37)=AJ(J,37)+AWS
         BJ(J,37)=BJ(J,37)+BWS
         CJ(J,37)=CJ(J,37)+CWS
         AJ(J,28)=AJ(J,28)+AWMG
         BJ(J,28)=BJ(J,28)+BWMG
         CJ(J,28)=CJ(J,28)+CWMG
         AJ(J,38)=AJ(J,38)+ATAUL/NSURF
         BJ(J,38)=BJ(J,38)+BTAUL/NSURF
         CJ(J,38)=CJ(J,38)+CTAUL/NSURF
         IF(MODDSF.NE.0) GO TO 7000                                     6673.   
         AJ(J,23)=AJ(J,23)+ATS                                          6674.   
         BJ(J,23)=BJ(J,23)+BTS                                          6675.   
         CJ(J,23)=CJ(J,23)+CTS                                          6676.   
         AJ(J,26)=AJ(J,26)+AT2                                          6674.   
         BJ(J,26)=BJ(J,26)+BT2                                          6675.   
         CJ(J,26)=CJ(J,26)+CT2                                          6676.   
c     print *,j,'ATS=',ATS,' AT2=',AT2
c     print *,'BLDATA'
c     print *,(BLDATA(1,j,k),k=1,3)
c     print *,(BLDATA(1,j,k),k=6,8)
 
 7000 CONTINUE                                                          6677.   
C****                                                                   6678.   
C**** ADD IN SURFACE FRICTION TO FIRST LAYER WIND                       6679.   
C****                                                                   6680.   
      DO 7600 I=1,IM                                                    6681.   
      U(I,2,1)=U(I,2,1)-2.*(DU1(1,1)*COSI(I)-DV1(1,1)*SINI(I))*RAPVN(1) 6682.   
      V(I,2,1)=V(I,2,1)-2.*(DV1(1,1)*COSI(I)+DU1(1,1)*SINI(I))*RAPVN(1) 6683.   
      U(I,JM,1)=U(I,JM,1)                                               6684.   
     *  -2.*(DU1(1,JM)*COSI(I)+DV1(1,JM)*SINI(I))*RAPVS(JM)             6685.   
 7600 V(I,JM,1)=V(I,JM,1)                                               6686.   
     *  -2.*(DV1(1,JM)*COSI(I)-DU1(1,JM)*SINI(I))*RAPVS(JM)             6687.   
      DO 7700 J=2,JMM1                                                  6688.   
      I=IM                                                              6689.   
      DO 7700 IP1=1,IM                                                  6690.   
      if(J.eq.JPR.or.J.eq.-12)then
        print *,' J=',J,' before'
       print *,'U(I,J,1)=',U(I,J,1),' V(I,J,1)=',V(I,J,1)
       print *,'U(I,J+1,1)=',U(I,J+1,1),' V(I,J+1,1)=',V(I,J+1,1)
       print *,'DU1(I,J)=',DU1(I,J),' DU1(IP1,J)=',DU1(IP1,J)
      endif
      U(I,J,1)=U(I,J,1)-(DU1(I,J)+DU1(IP1,J))*RAPVS(J)                  6691.   
      V(I,J,1)=V(I,J,1)-(DV1(I,J)+DV1(IP1,J))*RAPVS(J)                  6692.   
      U(I,J+1,1)=U(I,J+1,1)-(DU1(I,J)+DU1(IP1,J))*RAPVN(J)              6693.   
      V(I,J+1,1)=V(I,J+1,1)-(DV1(I,J)+DV1(IP1,J))*RAPVN(J)              6694.   
      if(J.eq.JPR.or.J.eq.-12)then
        print *,' J=',J,' after'
       print *,'U(I,J,1)=',U(I,J,1),' V(I,J,1)=',V(I,J,1)
       print *,'U(I,J+1,1)=',U(I,J+1,1),' V(I,J+1,1)=',V(I,J+1,1)
       print *,'DU1(I,J)=',DU1(I,J),' DU1(IP1,J)=',DU1(IP1,J)
      endif
 7700 I=IP1                                                             6695.   
c     print *,'U V'
c     do j=1,jm
c     print *,j,U(I,J,1),v(I,J,1)
c     enddo
C****                                                                   6696.   
C**** DRY CONVECTION ORIGINATING FROM THE FIRST LAYER                   6697.   
C****                                                                   6698.   
C**** LOAD U,V INTO UT,VT.  UT,VT WILL BE FIXED DURING DRY CONVECTION   6699.   
C****   WHILE U,V WILL BE UPDATED.                                      6700.   
      DO 8050 L=1,LM                                                    6701.   
      DO 8050 J=2,JM                                                    6702.   
      DO 8050 I=1,IM                                                    6703.   
      UT(I,J,L)=U(I,J,L)                                                6704.   
 8050 VT(I,J,L)=V(I,J,L)                                                6705.   
C**** OUTSIDE LOOPS OVER J AND I                                        6706.   
      DO 8500 J=1,JM                                                    6707.   
      POLE=.FALSE.                                                      6708.   
      IF(J.EQ.1.OR.J.EQ.JM) POLE=.TRUE.                                 6709.   
      IMAX=IM                                                           6710.   
      IF(POLE) IMAX=IM                                                  6711.   
      DO 8120 K=1,2                                                     6712.   
      RA(K)=RAPVS(J)                                                    6713.   
 8120 RA(K+2)=RAPVN(J)                                                  6714.   
      IM1=IM                                                            6715.   
      DO 8500 I=1,IMAX                                                  6716.   
      BLDATA(I,J,4)=1.                                                  6717.   
      IF(T(I,J,1)*(1.+Q(I,J,1)*RVX).LE.                                 6718.   
     *   T(I,J,2)*(1.+Q(I,J,2)*RVX)) GO TO 8500                         6719.   
C**** MIX HEAT AND MOISTURE THROUGHOUT THE BOUNDARY LAYER               6720.   
      PKMS=PK(I,J,1)*DSIG(1)+PK(I,J,2)*DSIG(2)                          6721.   
      THPKMS=T(I,J,1)*(PK(I,J,1)*DSIG(1))+T(I,J,2)*(PK(I,J,2)*DSIG(2))  6722.   
      QMS=Q(I,J,1)*DSIG(1)+Q(I,J,2)*DSIG(2)                             6723.   
      TVMS=T(I,J,1)*(1.+Q(I,J,1)*RVX)*(PK(I,J,1)*DSIG(1))               6724.   
     *    +T(I,J,2)*(1.+Q(I,J,2)*RVX)*(PK(I,J,2)*DSIG(2))               6725.   
      THETA=TVMS/PKMS                                                   6726.   

#if ( defined CPL_CHEM )
!
! --- 03/23/95
!
      cfc11ms = cfc11(i,j,1)*dsig(1) + cfc11(i,j,2)*dsig(2)

      cfc12ms = cfc12(i,j,1)*dsig(1) + cfc12(i,j,2)*dsig(2)

      xn2oms  = xn2o (i,j,1)*dsig(1) + xn2o (i,j,2)*dsig(2)

      o3ms    = o3   (i,j,1)*dsig(1) + o3   (i,j,2)*dsig(2)

      coms    = co   (i,j,1)*dsig(1) + co   (i,j,2)*dsig(2)

      zco2ms  = zco2 (i,j,1)*dsig(1) + zco2 (i,j,2)*dsig(2)

      xnoms   = xno  (i,j,1)*dsig(1) + xno  (i,j,2)*dsig(2)

      xno2ms  = xno2 (i,j,1)*dsig(1) + xno2 (i,j,2)*dsig(2)

      xn2o5ms = xn2o5(i,j,1)*dsig(1) + xn2o5(i,j,2)*dsig(2)

      hno3ms  = hno3 (i,j,1)*dsig(1) + hno3 (i,j,2)*dsig(2)

      ch4ms   = ch4  (i,j,1)*dsig(1) + ch4  (i,j,2)*dsig(2)

      ch2oms  = ch2o (i,j,1)*dsig(1) + ch2o (i,j,2)*dsig(2)

      so2ms   = so2  (i,j,1)*dsig(1) + so2  (i,j,2)*dsig(2)

      h2so4ms = h2so4(i,j,1)*dsig(1) + h2so4(i,j,2)*dsig(2)

        ! === if hfc, pfc, and sf6 are included:
#ifdef INC_3GASES
	! === 032698
      hfc134ams = hfc134a(i,j,1)*dsig(1) 
     &	        + hfc134a(i,j,2)*dsig(2)

      pfcms     = pfc(i,j,1)*dsig(1) 
     &	        + pfc(i,j,2)*dsig(2)

      sf6ms     = sf6(i,j,1)*dsig(1) 
     &	        + sf6(i,j,2)*dsig(2)
	! ===
#endif

      bcms  = bcarbon (i,j,1)*dsig(1) + bcarbon (i,j,2)*dsig(2)
      ocms  = ocarbon (i,j,1)*dsig(1) + ocarbon (i,j,2)*dsig(2)

c 062295
c      h2o2ms  = h2o2 (i,j,1)*dsig(1) + h2o2 (i,j,2)*dsig(2)

!
#endif

      DO 8140 L=3,LM                                                    6727.   
      IF(THETA.LT.T(I,J,L)*(1.+Q(I,J,L)*RVX)) GO TO 8160                6728.   
      PKMS=PKMS+(PK(I,J,L)*DSIG(L))                                     6729.   
      THPKMS=THPKMS+T(I,J,L)*(PK(I,J,L)*DSIG(L))                        6730.   
      QMS=QMS+Q(I,J,L)*DSIG(L)                                          6731.   
      TVMS=TVMS+T(I,J,L)*(1.+Q(I,J,L)*RVX)*(PK(I,J,L)*DSIG(L))          6732.   

#if ( defined CPL_CHEM )
!
! --- 03/23/95
!
        cfc11ms = cfc11ms + cfc11(i,j,l)*dsig(l)

        cfc12ms = cfc12ms + cfc12(i,j,l)*dsig(l)

        xn2oms  = xn2oms  + xn2o (i,j,l)*dsig(l)

        o3ms    = o3ms    + o3   (i,j,l)*dsig(l)

        coms    = coms    + co   (i,j,l)*dsig(l)

        zco2ms  = zco2ms  + zco2 (i,j,l)*dsig(l)

        xnoms   = xnoms   + xno  (i,j,l)*dsig(l)

        xno2ms  = xno2ms  + xno2 (i,j,l)*dsig(l)

        xn2o5ms = xn2o5ms + xn2o5(i,j,l)*dsig(l)

        hno3ms  = hno3ms  + hno3 (i,j,l)*dsig(l)

        ch4ms   = ch4ms   + ch4  (i,j,l)*dsig(l)

        ch2oms  = ch2oms  + ch2o (i,j,l)*dsig(l)

        so2ms   = so2ms   + so2  (i,j,l)*dsig(l)

        h2so4ms = h2so4ms + h2so4(i,j,l)*dsig(l)

        ! === if hfc, pfc, and sf6 are included:
#ifdef INC_3GASES
	! === 032698
        hfc134ams = hfc134ams 
     &		  + hfc134a(i,j,l)*dsig(l)

        pfcms     = pfcms 
     &		  + pfc(i,j,l)*dsig(l)

        sf6ms     = sf6ms 
     &		  + sf6(i,j,l)*dsig(l)
	! ===
#endif

	bcms  = bcms  + bcarbon (i,j,l)*dsig(l)
	ocms  = ocms  + ocarbon (i,j,l)*dsig(l)

c 062295
c        h2o2ms  = h2o2ms  + h2o2 (i,j,l)*dsig(l)
!
#endif

 8140 THETA=TVMS/PKMS                                                   6733.   
      L=LM+1                                                            6734.   
 8160 LMAX=L-1                                                          6735.   
      RDSIGS=1./(SIGE(1)-SIGE(LMAX+1))                                  6736.   
      THM=THPKMS/PKMS                                                   6737.   
      QMS=QMS*RDSIGS                                                    6738.   

#if ( defined CPL_CHEM )
!
! --- 03/23/95
!
      cfc11ms = cfc11ms*rdsigs

      cfc12ms = cfc12ms*rdsigs

      xn2oms  = xn2oms *rdsigs

      o3ms    = o3ms   *rdsigs

      coms    = coms   *rdsigs

      zco2ms  = zco2ms *rdsigs

      xnoms   = xnoms  *rdsigs

      xno2ms  = xno2ms *rdsigs

      xn2o5ms = xn2o5ms*rdsigs

      hno3ms  = hno3ms *rdsigs

      ch4ms   = ch4ms  *rdsigs

      ch2oms  = ch2oms *rdsigs

      so2ms   = so2ms  *rdsigs

      h2so4ms = h2so4ms*rdsigs

        ! === if hfc, pfc, and sf6 are included:
#ifdef INC_3GASES
	! === 032698
      hfc134ams = hfc134ams*rdsigs

      pfcms     = pfcms*rdsigs

      sf6ms     = sf6ms*rdsigs
	! ===
#endif

      bcms  = bcms*rdsigs
      ocms  = ocms*rdsigs

c 062295
c      h2o2ms  = h2o2ms*rdsigs
c
!
#endif

      BLDATA(I,J,4)=LMAX                                                6739.   
      DO 8180 L=1,LMAX                                                  6740.   
         AJL(J,L,12)=AJL(J,L,12)+(THM-T(I,J,L))*PK(I,J,L)*P(I,J)        6741.   
      T(I,J,L)=THM                                                      6742.   

#if ( defined CPL_CHEM )
!
! --- 03/23/95
!
        cfc11(i,j,l) = cfc11ms

        cfc12(i,j,l) = cfc12ms

        xn2o(i,j,l)  = xn2oms

        o3(i,j,l)    = o3ms

        co(i,j,l)    = coms

        zco2(i,j,l)  = zco2ms

        xno(i,j,l)   = xnoms

        xno2(i,j,l)  = xno2ms

        xn2o5(i,j,l) = xn2o5ms

        hno3(i,j,l)  = hno3ms

        ch4(i,j,l)   = ch4ms

        ch2o(i,j,l)  = ch2oms

        so2(i,j,l)   = so2ms

        h2so4(i,j,l) = h2so4ms

        ! === if hfc, pfc, and sf6 are included:
#ifdef INC_3GASES
! === 032698
        hfc134a(i,j,l) = hfc134ams

        pfc(i,j,l)     = pfcms

        sf6(i,j,l)     = sf6ms
	! ===
#endif

	bcarbon(i,j,l)  = bcms
	ocarbon(i,j,l)  = ocms

c 062295
c        h2o2(i,j,l)  = h2o2ms
c
!
#endif

 8180 Q(I,J,L)=QMS                                                      6743.   
      IF(POLE) GO TO 8300                                               6744.   
C**** MIX MOMENTUM THROUGHOUT THE BOUNDARY LAYER AT NON-POLAR GRID BOXES6745.   
      ID(1)=I+(J-1)*IM                                                  6748.   
      ID(2)=ID(1)+IM*JM*LM                                              6749.   
      ID(3)=I+J*IM                                                      6752.   
      ID(4)=ID(3)+IM*JM*LM                                              6753.   
      if(J.eq.JPR)then
       print *,'ID for J=',j
       print *,(ID(k),k=1,4)
       print *,'RA for J=',j
       print *,(RA(k),k=1,4)
      endif
      DO 8240 K=1,4                                                     6754.   
      UMS(K)=0.                                                         6755.   
      DO 8220 L=1,LMAX                                                  6756.   
 8220 UMS(K)=UMS(K)+UT(ID(K),1,L)*DSIG(L)                               6757.   
 8240 UMS(K)=UMS(K)*RDSIGS                                              6758.   
      DO 8260 L=1,LMAX                                                  6759.   
         AJL(J,L,38)=AJL(J,L,38)+(UMS(1)-UT(I,J,L))*.5*                 6760.   
     *     P(I,J)*RA(1)                                                 6761.   
         AJL(J+1,L,38)=AJL(J+1,L,38)+(UMS(3)-                           6762.   
     *      UT(I,J+1,L))*P(I,J)*RA(3)*.5                                6763.   
      DO 8260 K=1,4                                                     6764.   
      if(J.eq.JPR)then
       print *,'L=',L,' K=',K
       print *,'ID(K)=',ID(K),' RA(K)=',RA(K)
       print *,'UMS(K)=',UMS(K),' UT(ID(K),1,L)=',UT(ID(K),1,L)
      endif
 8260 U(ID(K),1,L)=U(ID(K),1,L)+(UMS(K)-UT(ID(K),1,L))*RA(K)            6765.   
      GO TO 8400                                                        6766.   
C**** MIX MOMENTUM THROUGHOUT THE BOUNDARY LAYER AT POLAR GRID BOXES    6767.   
 8300 JVPO=2                                                            6768.   
      IF(J.EQ.JM) JVPO=JM                                               6769.   
      RAPO=2.*RAPVN(1)                                                  6770.   
      DO 8360 IPO=1,IM                                                  6771.   
      UMSPO=0.                                                          6772.   
      VMSPO=0.                                                          6773.   
      DO 8320 L=1,LMAX                                                  6774.   
      UMSPO=UMSPO+UT(IPO,JVPO,L)*DSIG(L)                                6775.   
 8320 VMSPO=VMSPO+VT(IPO,JVPO,L)*DSIG(L)                                6776.   
      UMSPO=UMSPO*RDSIGS                                                6777.   
      VMSPO=VMSPO*RDSIGS                                                6778.   
      DO 8340 L=1,LMAX                                                  6779.   
      U(IPO,JVPO,L)=U(IPO,JVPO,L)+(UMSPO-UT(IPO,JVPO,L))*RAPO           6780.   
      V(IPO,JVPO,L)=V(IPO,JVPO,L)+(VMSPO-VT(IPO,JVPO,L))*RAPO           6781.   
 8340    AJL(JVPO,L,38)=AJL(JVPO,L,38)                                  6782.   
     *  +(UMSPO-UT(IPO,JVPO,L))*P(1,J)*RAPO                             6783.   
 8360 CONTINUE                                                          6784.   
C**** ACCUMULATE BOUNDARY LAYER DIAGNOSTICS                             6785.   
 8400    IF(MODD6.NE.0) GO TO 8500                                      6786.   
         DO 8420 KR=1,4                                                 6787.   
         IF(I.EQ.IJD6(1,KR).AND.J.EQ.IJD6(2,KR)) GO TO 8440             6788.   
 8420    CONTINUE                                                       6789.   
         GO TO 8500                                                     6790.   
 8440    ADAILY(IHOUR,47,KR)=ADAILY(IHOUR,47,KR)+1.                     6791.   
         ADAILY(IHOUR,48,KR)=ADAILY(IHOUR,48,KR)+LMAX                   6792.   
 8500 IM1=I                                                             6793.   
      do j=1,jm
      I=1
      if(J.eq.JPR.or.J.eq.-12)then
        print *,' J=',J,' after dry convection'
       print *,'U(I,J,1)=',U(I,J,1),' V(I,J,1)=',V(I,J,1)
       print *,'U(I,J+1,1)=',U(I,J+1,1),' V(I,J+1,1)=',V(I,J+1,1)
      endif
      enddo
 9000 CONTINUE                                                          6794.   
      do 9001 J=1,JM
       TSURFD(J)=TSURFD(J)+(BLDATA(1,J,2)-273.16)/24.
 9001 continue
c     write (935) ,ps4clm,
c    &  tsl4clm,
c    &  qs4clm,ws4clm
c    & ,us4clm,vs4clm
      RETURN                                                            6795.   
  990 FORMAT ('0PPBL',3I4,14F8.2)                                       6818.   
  991 FORMAT ('0SURFACE ',4I4,5F10.4,3F11.7)                            6819.   
  992 FORMAT ('0',I2,10F10.4/23X,4F10.4,10X,2F10.4/                     6820.   
     *  33X,3F10.4,10X,2F10.4)                                          6821.   
  993 FORMAT ('0',I2,10F10.4/23X,7F10.4/33X,7F10.4)                     6822.   
  994 FORMAT ('0',I2,11F10.4)                                           6823.   
      END                                                               6824.   
