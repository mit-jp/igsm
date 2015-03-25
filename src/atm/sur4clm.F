
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

      SUBROUTINE SUR4CLM

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
      LOGICAL POLE,PRNT,HPRNT
      common/conprn/HPRNT,JPR,LPR
!     common/TSUR/TSURFC(JM0,0:13),TSURFT(JM0),TSURFD(JM0),DTEMSR(JM0)
#include "TSRF.COM"
      common/SURRAD/TRSURF(JM0,4),SRSURF(JM0,4)
c     REAL*8 B,TGV,TKV,TSV0,TSV1,TSV                                    5818.   
      integer IQ1,IQ2,IQ3
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
C****                                                                   5874.   
      save
c     print *,'sur4clm TAU=',TAU
      NSTEPS=NSURF*NSTEP/NDYN                                           5875.   
      IF(IFIRST.NE.1) GO TO 50                                          5876.   
      ZGSL = 100.
      ZGSL = 30.
      print *,' SURFACE FOR CLM'
!     print *,' ZGS=30 m for LAND '
      print *,' ZGS=',ZGSL,' m for LAND '
      WMGMIN=0.
#if ( defined CLM35 )
      WMGMIN=3.
#else
      WMGMIN=8.
#endif
      print *,'WMGMIN 4 LAND=',WMGMIN
      IFIRST=0                                                          5877.   
      print *,' WMGE'
        print 258,(WMGE(1,J)+WMGMIN,J=1,JM)
  258 format(12f5.1)
C     SRCORX=1.                                                         5878.13 
      CIAX=0.3
      print *,' surfacen  CIAX=',CIAX
      print *,' QS=Q1, TS=T1'
      print *,' WS=sqrt(0.75*W1+WGEM) '
      IQ1=IM/4+1                                                        5881.   
      IQ2=IM/2+1                                                        5882.   
      IQ3=3*IM/4+1                                                      5883.   
!     DTSURF=NDYN*DT/NSURF                                              5884.   
!     print *,' From SRF4CLM DTSURF=',DTSURF
!        DTSRCE=DT*NDYN                                                 5885.   
      SHA=RGAS/KAPA                                                     5886.   
      RVX=0.                                                            5887.   
   50 CONTINUE
C**** ZERO OUT ENERGY AND EVAPORATION FOR GROUND AND INITIALIZE TGRND   5906.   
      DO 70 J=1,JM                                                      5907.   
      DO 70 I=1,IM                                                      5908.   
      TGRND(I,J,3)=GDATA(I,J,13)                                        5910.   
      TGRND(I,J,4)=GDATA(I,J,4)                                         5911.   
   70 CONTINUE
c     print *,'After 70'
         IHOUR=1.5+TOFDAY                                               5914.   
C****                                                                   5915.   
C**** OUTSIDE LOOP OVER TIME STEPS, EXECUTED NSURF TIMES EVERY HOUR     5916.   
C****                                                                   5917.   
C****                                                                   5927.   
C**** OUTSIDE LOOP OVER J AND I, EXECUTED ONCE FOR EACH GRID POINT      5928.   
C****                                                                   5929.   
      JPR=-7
      DO 7000 J=1,JM                                                    5930.   
c     print *,'After 7000 J=',J
c     print *,IQ3
      PRNT=.FALSE.
      HEMI=1.                                                           5931.   
      IF(J.LE.JM/2) HEMI=-1.                                            5932.   
      IF(J.EQ.1) GO TO 80                                               5936.   
      IF(J.EQ.JM) GO TO 90                                              5937.   
      WMG0=.5*(WMGE(1,J)+WMGE(1,J+1))+.001                              5937.5  
      POLE=.FALSE.                                                      5938.   
      IMAX=IM                                                           5939.   
      GO TO 100                                                         5940.   
C**** CONDITIONS AT THE SOUTH POLE                                      5941.   
   80 POLE=.TRUE.                                                       5942.   
c     print *,'After 80'
c     print *,IQ1,IQ2,IQ3
      IMAX=1                                                            5943.   
      JVPO=2                                                            5944.   
      RAPO=2.*RAPVN(1)                                                  5945.   
c     print *,' RAPO=', RAPO
c     II1=IQ1
c     II2=IQ2
c     II3=IQ3
c     print *,II1,II2,II3
c     print *,' III=',III
c     print *,' U(IQ3,2,1)=',U(IQ3,2,1)
      U1=.25*(U(1,2,1)+V(IQ1,2,1)-U(IQ2,2,1)-V(IQ3,2,1))                5946.   
      V1=.25*(V(1,2,1)-U(IQ1,2,1)-V(IQ2,2,1)+U(IQ3,2,1))                5947.   
      WMG0=WMGE(1,2)                                                    5947.5  
      GO TO 100                                                         5948.   
C**** CONDITIONS AT THE NORTH POLE                                      5949.   
   90 POLE=.TRUE.                                                       5950.   
      IMAX=1                                                            5951.   
      JVPO=JM                                                           5952.   
      RAPO=2.*RAPVS(JM)                                                 5953.   
      U1=.25*(U(1,JM,1)-V(IQ1,JM,1)-U(IQ2,JM,1)+V(IQ3,JM,1))            5954.   
      V1=.25*(V(1,JM,1)+U(IQ1,JM,1)-V(IQ2,JM,1)-U(IQ3,JM,1))            5955.   
      WMG0=WMGE(1,JM)                                                   5955.5  
C**** ZERO OUT SURFACE DIAGNOSTICS WHICH WILL BE SUMMED OVER LONGITUDE  5956.   
  100 CONTINUE
c     print *,'After 100'
      BTS=0.     
      BWS=0.
      BWMG=0.
      IM1=IM                                                            5969.   
       i=1
       tsl4clm(i,j)=0.0
       qs4clm(i,j)=0.0
       ps4clm(i,j)=0.0
       ws4clm(i,j)=0.0
       us4clm(i,j)=0.0
       vs4clm(i,j)=0.0
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
      IF(POLE) GO TO 1200                                               5993.   
      U1=.25*(U(IM1,J,1)+U(I,J,1)+U(IM1,J+1,1)+U(I,J+1,1))              5994.   
      V1=.25*(V(IM1,J,1)+V(I,J,1)+V(IM1,J+1,1)+V(I,J+1,1))              5995.   
 1200 TH1=T(I,J,1)                                                      5996.   
      Q1=Q(I,J,1)                                                       5997.   
      THV1=TH1*(1.+Q1*RVX)                                              5998.   
C**** ZERO OUT QUANTITIES TO BE SUMMED OVER SURFACE TYPES               6002.   
      USS=0.                                                            6003.   
      VSS=0.                                                            6004.   
      WSS=0.                                                            6005.   
      TSS=0.                                                            6006.   
      QSS=0.                                                            6007.   
C****                                                                   6032.   
 2400 IF(PLAND.LE.0.) GO TO 5000                                        6074.   
!     ZGS=30.                                                           6078.   
      ZGS=ZGSL
      IF(PLICE.LE.0.) GO TO 2600                                        6080.   
C****                                                                   6081.   
C**** LAND ICE                                                          6082.   
C****                                                                   6083.   
      ITYPE=3                                                           6084.   
      PTYPE=PLICE                                                       6085.   
      TG1=TGRND(I,J,3)                                                  6087.   
c     ELHX=LHS                                                          6094.   
      if (TG1.le.0.0)ELHX=LHS
      if (TG1.gt.0.0)ELHX=LHE
      GO TO 3000                                                        6095.   
C****                                                                   6096.   
 2600 IF(PEARTH.LE.0.) GO TO 5000                                       6097.   
C****                                                                   6098.   
C**** EARTH                                                             6099.   
C****                                                                   6100.   
      ITYPE=4                                                           6101.   
      PTYPE=PEARTH                                                      6102.   
      TG1=TGRND(I,J,4)                                                  6104.
      if (TG1.le.0.0)ELHX=LHS
      if (TG1.gt.0.0)ELHX=LHE
C****                                                                   6134.   
C**** BOUNDARY LAYER INTERACTION                                        6135.   
C****                                                                   6136.   
 3000 continue
      if(J.eq.JPR)then
        print *,' after 3000'
        print *,'TAU=',TAU,' NS=',NS,' ITYPE=',ITYPE
        print *,'CDH=',CDH,' RGAS=',RGAS
        print *,'U1=',U1,' V1=',V1
        print *,'WMGO=',WMGO
      endif
      TKV=THV1*PSK                                                      6137.   
C**** LOOP OVER GROUND TIME STEPS                                       6148.   
      TG=TG1+TF                                                         6150.   
      QG=QSAT(TG,PS,ELHX)                                               6151.   
      TGV=TG*(1.+QG*RVX)                                                6152.   
      UG=0.75*U1
      VG=0.75*V1
      W1=SQRT(UG*UG+VG*VG)
      WMG=WMG0+WMGMIN
      WS=SQRT(W1**2+WMG)
      RW=1.0
      if(W1.ne.0.0)RW=WS/W1
      THS=TH1
      QS=Q1
      TSV=THS*PSK
      RIGS=ZGS*GRAV*(TSV-TGV)/(TGV*WS*WS)
      IF(RIGS.LE.0) THEN
C     surface layer has unstable stratification
         CIA=TWOPI*0.0625/(1.+WS*CIAX)
      ELSE
C     surface layer has stable stratification
         CIA=TWOPI*(0.09375-0.03125/(1.+4*RIGS**2))/(1.+WS*CIAX)
      END IF
      if(J.eq.JPR)then
        print *,'TAU=',TAU,' NS=',NS,' ITYPE=',ITYPE
        print *,'WS=',WS,' ZGS=',ZGS
        print *,'RIGS=',RIGS,' TGV=',TGV
      endif
      USR=COS(CIA)
      VSR=SIN(CIA)*HEMI
      US=(USR*UG-VSR*VG)
      VS=(VSR*UG+USR*VG)
      if(J.eq.JPR)then
c       print *,' '
        print *,'TAU=',TAU,' NS=',NS,' ITYPE=',ITYPE
        print *,'CDH=',CDH,' RGAS=',RGAS
        print *,'PS=',PS,' TSV=',TSV
        print *,'WS=',WS
      endif
      TS=TSV/(1.+QS*RVX)                                                6467.   
#if ( defined CLM35 )
!!!   QSATS=QSAT(TS,PS,ELHX)                                            6468.   
!!!   IF(QS.LE.QSATS) GO TO 3500                                        6469.   
!!!   DQSSDT=QSATS*ELHX/(RVAP*TS*TS)                                    6470.   
!!!   X=(QS-QSATS)/(DQSSDT+(SHA/ELHX))                                  6471.   
!     print *,j,TS,QS
!!!   TS=TS+X                                                           6472.   
!!!   QS=QS+X*(SHA/ELHX)                                                6473.   
!     print *,TS,QS
 3500 CONTINUE
#endif

          if(ITYPE.EQ.4.or.ITYPE.EQ.3)then
           tsl4clm(i,j)=tsl4clm(i,j)+TS*PTYPE/PLAND
           qs4clm(i,j)=qs4clm(i,j)+QS*PTYPE/PLAND
           ps4clm(i,j)=ps4clm(i,j)+PS*PTYPE/PLAND
           ws4clm(i,j)=ws4clm(i,j)+WS*PTYPE/PLAND
!          us4clm(i,j)=us4clm(i,j)+US*PTYPE/PLAND
!          vs4clm(i,j)=vs4clm(i,j)+VS*PTYPE/PLAND
           us4clm(i,j)=us4clm(i,j)+RW*US*PTYPE/PLAND
           vs4clm(i,j)=vs4clm(i,j)+RW*VS*PTYPE/PLAND
         endif

      TSSFC(I,J,ITYPE)=TS                                               6521.5  
      USS=USS+US*PTYPE                                                  6524.   
      VSS=VSS+VS*PTYPE                                                  6525.   
      WSS=WSS+WS*PTYPE                                                  6526.   
      TSS=TSS+TS*PTYPE                                                  6527.   
      QSS=QSS+QS*PTYPE                                                  6528.   
      GO TO (5000,5000,4400,4600),ITYPE                                 6551.   
C****                                                                   6552.   
C**** LAND ICE                                                          6569.   
C****                                                                   6570.   
 4400    CONTINUE
         BTS=BTS+(TS-TF)*PLICE                                          6574.   
         BWS=BWS+WS*PLICE
         BWMG=BWMG+SQRT(WMG)*PLICE
      GO TO 2600                                                        6575.   
C****                                                                   6576.   
C**** EARTH                                                             6577.   
C****                                                                   6578.   
 4600    CONTINUE
         BTS=BTS+(TS-TF)*PEARTH                                         6582.   
         BWS=BWS+WS*PEARTH
         BWMG=BWMG+SQRT(WMG)*PEARTH

C**** NON-OCEAN POINTS WHICH ARE NOT MELTING OR FREEZING WATER USE      6583.   
C****   IMPLICIT TIME STEPS                                             6584.   
C****                                                                   6585.   
C**** UPDATE SURFACE AND FIRST LAYER QUANTITIES                         6586.   
C****                                                                   6587.   
 5000 CONTINUE
C****                                                                   6596.   
C**** ACCUMULATE DIAGNOSTICS                                            6597.   
C****                                                                   6598.   
 6000 IM1=I                                                             6662.   
      WSSL(J)=WSS     
      TSSL(J)=TSS    
      QSSL(J)=QSS   
      USSL(J)=USS  
      VSSL(J)=VSS   
C**** QUANTITIES ACCUMULATED FOR SURFACE TYPE TABLES IN DIAG1           6663.   
         BLJ(J,37)=BWS
         BLJ(J,28)=BWMG
         BLJ(J,23)=BTS  
 7000 CONTINUE                                                          6677.   
      if(HPRNT) then
      print *,' From SRF4CLM TAU=',TAU,' JPR=',JPR
!          print *,'TS(1),TS(2)'
           print *,tsl4clm(1,jpr),qs4clm(1,jpr),ps4clm(1,jpr)
           print *,vs4clm(i,jpr),vs4clm(i,jpr),ws4clm(i,jpr)
           print *,dsw4clm(i,jpr),dlw4clm(i,jpr)
           print *,swinr4clm(i,jpr), swvis4clm(i,jpr)
      endif
C****                                                                   6678.   
#endif
      RETURN                                                            6795.   
      END                                                               6824.   
