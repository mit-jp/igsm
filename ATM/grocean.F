
#include "ctrparam.h"

!	==========================================================
!
!	GROUND.F: THIS SUBROUTINE USES THE SURFACE FLUXES TO 
!		    PREDICT IN TIME THE GROUND TEMPERATURE, GROUND 
!		    WATER AND ICE, AND SNOW MELTING.
!
!	----------------------------------------------------------
! 
!	Revision History:
!	
!	When	Who		What
!	----	----------	-------	
!	073100	Chien Wang	add cpp & dmax dmin -> max, min
!
!	==========================================================

      SUBROUTINE GROCEAN(mndriver)                                      7001.   
C****                                                                   7001.5  
C**** THIS SUBROUTINE USES THE SURFACE FLUXES TO PREDICT IN TIME THE    7002.   
C**** GROUND TEMPERATURE, GROUND WATER AND ICE, AND SNOW MELTING.       7002.5  
C****                                                                   7003.   

#include "BD2G04.COM"                                      		7003.5  

      COMMON U,V,T,P,Q                                                  7004.   
      COMMON/WORK1/CONV(IM0,JM0,LM0),PK(IM0,JM0,LM0),PREC(IM0,JM0),
     &   TPREC(IM0,JM0)                                                 4005.5
      COMMON/WORK3/E0(IO0,JM0,4),E1(IO0,JM0,4),EVAPOR(IO0,JM0,4)        7004.5  
      COMMON/OT/OTA(IO0,JM0),OTB(IO0,JM0),OTC(IO0,JM0)                  7005.   
      COMMON/SPEC2/KM,KINC,COEK,C3LAND(IO0,JM0),C3OICE(IO0,JM0)         7005.1  
     *  ,C3LICE(IO0,JM0),WMGE(IO0,JM0),TSSFC(1,JM0,4)                   7005.2  
      common/qfl/QFLUX(JM0,0:13),ZOAV(JM0),QFLUXT(JM0)
      COMMON/OLDZO/Z1OOLD(IO0,JM0)
      DIMENSION FWATER(JM0),TOBS(JM0),TF68(JM0),DELR(JM0)               7005.3  
      DATA SHV/0./,SHW/4185./,SHI/2060./,RHOW/1000./,RHOI/916.6/,       7005.9  
     *  ALAMI/2.1762/,TFO/-1.56/,Z1I/.1/,Z2LI/2.9/,Z1E/.1/,Z2E/4./      7006.   
      INTEGER JDOFM(13)                                                 7006.1  
      DATA JDOFM/0,31,59,90,120,151,181,212,243,273,304,334,365/        7006.2  
      DATA Z2OIM/0.9/,Z2OIX/4.9/                                        7006.5  
      DATA TTRUNC/0./                                                   7007.   
      DATA IFIRST/1/                                                    7007.5  
C****                                                                   7008.   
C**** FDATA  2  LAND COVERAGE (1)                                       7008.5  
C****        3  RATIO OF LAND ICE COVERAGE TO LAND COVERAGE (1)         7009.   
C****                                                                   7009.5  
C**** ODATA  1  OCEAN TEMPERATURE (C)                                   7010.   
C****        2  RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)       7010.5  
C****        3  OCEAN ICE AMOUNT OF SECOND LAYER (KG/M**2)              7011.   
C****                                                                   7011.5  
C**** GDATA  1  OCEAN ICE SNOW AMOUNT (KG/M**2)                         7012.   
C****        2  EARTH SNOW AMOUNT (KG/M**2)                             7012.5  
C****        3  OCEAN ICE TEMPERATURE OF FIRST LAYER (C)                7013.   
C****        4  EARTH TEMPERATURE OF FIRST LAYER (C)                    7013.5  
C****        5  EARTH WATER OF FIRST LAYER (KG/M**2)                    7014.   
C****        6  EARTH ICE OF FIRST LAYER (KG/M**2)                      7014.5  
C****        7  OCEAN ICE TEMPERATURE OF SECOND LAYER (C)               7015.   
C****        8  EARTH TEMPERATURE OF SECOND LAYER (C)                   7015.5  
C****        9  EARTH WATER OF SECOND LAYER (KG/M**2)                   7016.   
C****       10  EARTH ICE OF SECOND LAYER (KG/M**2)                     7016.5  
C****       12  LAND ICE SNOW AMOUNT (KG/M**2)                          7017.   
C****       13  LAND ICE TEMPERATURE OF FIRST LAYER (C)                 7017.5  
C****       14  LAND ICE TEMPERATURE OF SECOND LAYER (C)                7018.   
C****                                                                   7018.5  
C**** VDATA  9  WATER FIELD CAPACITY OF FIRST LAYER (KG/M**2)           7019.   
C****       10  WATER FIELD CAPACITY OF SECOND LAYER (KG/M**2)          7019.5  
C****                                                                   7020.   
      IF(IFIRST.NE.1) GO TO 50                                          7020.5  
      IFIRST=0                                                          7021.   
      FIO=IO                                                            7021.1  
      JDAYPR=0
      IF(KOCEAN.NE.1) GO TO 10                                          7021.5  
      DO 30 ii=1,12
      DO 30 J=1,JM                                                      7022.88
      SUMJ=SUMJ+QFLUX(J,ii)*DXYP(J)/12.
   30 continue
      SUMJ=SUMJ/AREAG                                                   7022.9
      WRITE(6,31) SUMJ                                                  7022.91
   31 FORMAT(1X,'GLOBAL ANNUAL MEAN QFLUX ',E15.4)   
   10 DTSRCE=NDYN*DT                                                    7023.   
      ACE1I=Z1I*RHOI                                                    7023.5  
      AC2OIM=Z2OIM*RHOI                                                 7024.   
      ATRUNC=2.**(-13)                                                  7024.5  
      BYZICX=1./(Z1I+Z2OIX)                                             7024.6  
      HC1I=ACE1I*SHI                                                    7025.   
      HC2LI=Z2LI*RHOI*SHI                                               7025.5  
      HC1DE=Z1E*1129950.                                                7026.   
      HC2DE=Z2E*1129950.+3.5*.125*RHOW*3100.                            7026.5  
      DIFFUS=DTSRCE/SDAY                                                7027.   
C     OTCOR=-.927E18                                                    7027.1  
      print *,'From GROCEAN PREC over ocean'
      print *,'E0(1,j,1)'
      print *,(E0(1,j,1),j=1,JM0)
   50 ANGLE=TWOPI*JDAY/365.                                             7027.5  
      DO 55 MONTH=1,12                                                  7027.6  
      IF(JDAY.LE.JDOFM(MONTH+1)) GO TO 56                               7027.7  
   55 CONTINUE                                                          7027.8  
   56 CONTINUE                                                          7027.9  
      SINANG=SIN(ANGLE)                                                 7028.   
      COSANG=COS(ANGLE)                                                 7028.5  
C****                                                                   7029.   
C**** OUTSIDE LOOP OVER J AND I, EXECUTED ONCE FOR EACH GRID POINT      7029.5  
C****                                                                   7030.   
      JRPR=0
      DO 980 J=1,JM                                                     7030.5  
      if(J.eq.-3.and.TAU.ge. 17520.00)then
         print *,'ODATA(1,J,1)=',ODATA(1,J,1)
         print *,'ODATA(1,J,2)=',ODATA(1,J,2)
         print *,'ODATA(1,J,3)=',ODATA(1,J,3)
      endif
      IMAX=IM                                                           7031.   
      IF((J.EQ.1).OR.(J.EQ.JM)) IMAX=1                                  7031.5  
         BF1DT=0.                                                       7032.   
         CF1DT=0.                                                       7032.5  
         AOTDT=0.                                                       7033.   
         COTDT=0.                                                       7033.5  
         AEFO=0.                                                        7034.   
         CEFI=0.                                                        7034.5  
         BEDIFS=0.                                                      7035.   
         CEDIFS=0.                                                      7035.5  
         BERUN0=0.                                                      7036.   
         CF2DT=0.                                                       7036.5  
         BERUN2=0.                                                      7037.   
         CERUN2=0.                                                      7037.5  
         AERUN4=0.                                                      7038.   
         CERUN4=0.                                                      7038.5  
         ATG1=0.                                                        7039.   
         BTG1=0.                                                        7039.5  
         CTG1=0.                                                        7040.   
         ATG2=0.                                                        7040.5  
         BTG2=0.                                                        7041.   
         CTG2=0.                                                        7041.5  
         ATG3=0.                                                        7042.   
         AEVAP=0.                                                       7042.5  
         BEVAP=0.                                                       7043.   
         CEVAP=0.                                                       7043.5  
         BDIFS=0.                                                       7044.   
         CDIFS=0.                                                       7044.5  
         AIFO=0.                                                        7045.   
         CIFI=0.                                                        7045.5  
         BRUN0=0.                                                       7046.   
         CRUN0=0.                                                       7046.5  
         BRUN2=0.                                                       7047.   
         CRUN2=0.                                                       7047.5  
         ARUN4=0.                                                       7048.   
         CRUN4=0.                                                       7048.5  
         BWTR1=0.                                                       7049.   
         BACE1=0.                                                       7049.5  
         BWTR2=0.                                                       7050.   
         BACE2=0.                                                       7050.5  
         CACE2=0.                                                       7051.   
         BSNOW=0.                                                       7051.5  
         CSNOW=0.                                                       7052.   
         CICOV=0.                                                       7052.5  
      DO 960 I=1,IMAX                                                   7053.   
C****                                                                   7053.5  
C**** DETERMINE SURFACE CONDITIONS                                      7054.   
C****                                                                   7054.5  
      PLAND=FDATA(I,J,2)                                                7055.   
      PWATER=1.-PLAND                                                   7055.5  
      PLICE=FDATA(I,J,3)*PLAND                                          7056.   
      PEARTH=PLAND-PLICE                                                7056.5  
      ROICE=ODATA(I,J,2)                                                7057.   
      POICE=ROICE*PWATER                                                7057.5  
      POCEAN=PWATER-POICE                                               7058.   
      if(POCEAN.LE.1.E-5)then
         POCEAN=0.
         POICE=PWATER
      endif
         JR=J
         DXYPJ=DXYP(J)                                                  7059.   
         SNOWS=0.                                                       7059.5  
         WTR1S=0.                                                       7060.   
         ACE1S=0.                                                       7060.5  
         WTR2S=0.                                                       7061.   
         ACE2S=0.                                                       7061.5  
         TG1S=0.                                                        7062.   
         TG2S=0.                                                        7062.5  
         EVAPS=0.                                                       7063.   
         RUN0S=0.                                                       7063.5  
         DIFSS=0.                                                       7064.   
C**** CALCULATE PRECIPITATION HEAT FLUX (FALLS AT 0 DEGREES CENTIGRADE) 4044.5
      PRCP=PREC(I,J)                                                    4045.
          if (fland.lt.1.0)then
           PRCP=PRCP*(1.-fland*prlnd2total(j,mndriver))
     &      /(1.-fland)
           endif
      TPRCP=TPREC(I,J)                                                  4045.5
      IF(TPRCP.LT.0.) GO TO 530                                         4046.
C     EPRCP=PRCP*TPRCP*SHW                                              4046.5
      EPRCP=0.                                                          4047.
      ENRGP=EPRCP                                                       4047.5
      GO TO 550                                                         4048.
C     EPRCP=PRCP*TPRCP*SHI                                              4048.5
  530 EPRCP=0.                                                          4049.
      ENRGP=EPRCP-PRCP*LHM                                              4049.5
         AIJ(I,J,70)=AIJ(I,J,70)+PRCP                                   4050.
  550 CONTINUE
c     if(TAU.ge.17520.0.and.TAU.lt.17524.0)then
c      print *,'GROCEAN TAU=',TAU
c      print *,'PRCP=',PRCP,' TPRCP=',TPRCP
c      print *,'ENRGP=',ENRGP
c     endif
C****
C****                                                                   7064.5  
      IF(PWATER.LE.0.0) GO TO 400                                       7065.   
C****                                                                   7065.5  
C**** OCEAN                                                             7066.   
C****                                                                   7066.5  
      EVAP=EVAPOR(I,J,1)                                                7067.   
         ATG1=ATG1+ODATA(I,J,1)*POCEAN                                  7067.5  
         TG1S=TG1S+ODATA(I,J,1)*POCEAN                                  7068.   
         AEVAP=AEVAP+EVAP*POCEAN                                        7068.5  
         EVAPS=EVAPS+EVAP*POCEAN                                        7069.   
         AIJ(I,J,61)=AIJ(I,J,61)+EVAP*POCEAN                            7069.5  
      IF(KOCEAN.EQ.1) GO TO 60                                          7070.   
         ATG2=ATG2+ODATA(I,J,1)*POCEAN                                  7070.5  
         TG2S=TG2S+ODATA(I,J,1)*POCEAN                                  7071.   
      IF(POICE.GT.0.) GO TO 110                                         7071.5  
      GO TO 400                                                         7072.   
   60 TGW=ODATA(I,J,1)                                                  7072.5  
      OTDT=DTSRCE*QFLUXT(J)
         AIJ(I,J,57)=AIJ(I,J,57)+TGW                                    7073.   
      WTRO=Z1O(I,J)*RHOW                                                7073.5  
      ENRGO0=WTRO*TGW*SHW                                               7074.   
      EOFRZ=WTRO*TFO*SHW                                                7074.5  
      ENRGFO=0.
      ENRGO=0.
      ACEFO=0.
      if(POCEAN.LE.0.0)go to 100
      F0DT=E0(I,J,1)                                                    7075.   
         AIJ(I,J,65)=AIJ(I,J,65)+F0DT*POCEAN                            7075.5  
         ATG2=ATG2+ODATA(I,J,4)*POCEAN                                  7076.5  
         TG2S=TG2S+ODATA(I,J,4)*POCEAN                                  7077.   
         ATG3=ATG3+ODATA(I,J,5)*POCEAN                                  7077.5  
         AOTDT=AOTDT+OTDT*POCEAN                                        7078.   
CPRE  RUN4=-EVAP                                                        7078.5  
      RUN4=PRCP-EVAP                                                    7078.5  
      if(J.eq.-3.and.TAU.ge. 17520.00)then
         print *,'EVAP=',EVAP,' PRCP=',PRCP
      endif
      ERUN4=RUN4*TGW*SHW                                                7079.   
         AERUN4=AERUN4+ERUN4*POCEAN                                     7079.5  
         ARUN4=ARUN4+RUN4*POCEAN                                        7080.   
CPRE  ENRGO=F0DT+OTDT-ERUN4                                             7080.5  
      ENRGO=F0DT+OTDT+ENRGP-ERUN4                                       7080.5  
      IF(ENRGO0+ENRGO.LT.EOFRZ) GO TO 80                                7081.   
C**** FLUXES RECOMPUTE TGO WHICH IS ABOVE FREEZING POINT FOR OCEAN      7081.5  
      ENRGFO=0.                                                         7082.   
      ACEFO=0.                                                          7082.5  
      IF(ROICE.GT.0.) GO TO 100                                         7083.   
      ODATA(I,J,1)=TGW+(ENRGO/(WTRO*SHW)+TTRUNC)                        7083.5  
      if(J.eq.-3.and.TAU.ge. 17520.00)then
         print *,'1 ODATA(i,j,1)=',ODATA(i,j,1)
         print *,'TGW=',TGW,' WTRO=',WTRO
         print *,'RHOW=',RHOW,' Z1O(i,j)=',Z1O(i,j)
         print *,'SHW=',SHW
      endif
      GO TO 400                                                         7084.   
C**** FLUXES COOL TGO TO FREEZING POINT FOR OCEAN AND FORM SOME ICE     7084.5  
   80 ACEFO=(ENRGO0+ENRGO-EOFRZ)/(TFO*(SHI-SHW)-LHM)                    7085.   
      ENRGFO=ACEFO*(TFO*SHI-LHM)                                        7085.5  
         AEFO=AEFO-ENRGFO*POCEAN                                        7086.   
         AIFO=AIFO-ACEFO*POCEAN                                         7086.5  
      IF(ROICE.GT.0.) GO TO 100                                         7087.   
      ROICE=ACEFO/(ACE1I+AC2OIM)                                        7087.5  
      ODATA(I,J,1)=TFO                                                  7088.   
c     if(J.eq.8.and.TAU.ge. 17520.00)then
c        print *,'2 ODATA(i,j,1)=',ODATA(i,j,1)
c     endif
      if(J.eq.-3.and.TAU.ge. 17520.00)then
         print *,'ROICE==',ROICE
      endif
      ODATA(I,J,2)=ROICE                                                7088.5  
      GDATA(I,J,1)=0.                                                   7089.   
      GDATA(I,J,3)=TFO                                                  7089.5  
      GDATA(I,J,7)=TFO                                                  7090.   
      ODATA(I,J,3)=AC2OIM                                               7090.5  
      GO TO 400                                                         7091.   
C****                                                                   7091.5  
  100 ACE2F=0.                                                          7092.   
      ACE2M=0.                                                          7092.5  
C****                                                                   7093.   
C**** OCEAN ICE                                                         7093.5  
C****                                                                   7094.   
  110 SNOW=GDATA(I,J,1)                                                 7094.5  
      TG1=GDATA(I,J,3)                                                  7095.   
      TG2=GDATA(I,J,7)                                                  7095.5  
      ACE2=ODATA(I,J,3)                                                 7096.   
C        AIJ(I,J,1)=AIJ(I,J,1)+POICE                                    7096.5  
         AIJ(I,J,58)=AIJ(I,J,58)+ACE2*POICE                             7097.   
      F0DT=E0(I,J,2)                                                    7097.5  
         AIJ(I,J,66)=AIJ(I,J,66)+F0DT*POICE                             7098.   
      F1DT=E1(I,J,2)                                                    7098.5  
      EVAP=EVAPOR(I,J,2)                                                7099.   
         AIJ(I,J,62)=AIJ(I,J,62)+EVAP*POICE                             7099.5  
      Z2=ACE2/RHOI                                                      7100.   
         IF(KOCEAN.NE.1) GO TO 120                                      7100.5  
      WTRI0=WTRO-(SNOW+ACE1I+ACE2)                                      7101.   
      EIW0=WTRI0*TGW*SHW                                                7101.5  
      WTRW0=WTRO-ROICE*(SNOW+ACE1I+ACE2)                                7102.   
      ENRGW0=WTRW0*TGW*SHW                                              7102.5  
      RUN0=0.                                                           7103.   
      DIFSI=0.                                                          7103.5  
      EDIFSI=0.                                                         7104.   
      RUN4=-EVAP                                                        7104.5  
      ERUN4=TGW*RUN4*SHW                                                7105.   
         CERUN4=CERUN4+ERUN4*POICE                                      7105.5  
         CRUN4=CRUN4+RUN4*POICE                                         7106.   
C****                                                                   7106.5  
C**** OCEAN ICE, CALCULATE TG1                                          7107.   
C****                                                                   7107.5  
  120 SNANDI=SNOW+ACE1I-EVAP                                            7108.   
      if(TPRCP.lt.0.0)then
        SNANDI=SNANDI+PRCP
      endif
      HC1=SNANDI*SHI                                                    7108.5  
      ENRG1=F0DT+EVAP*(TG1*SHI-LHM)-F1DT                                7109.   
C PRE
      ENRG1=ENRG1+ENRGP
      if(J.eq.-3.and.TAU.ge. 17520.00)then
         print *,'ENRG1=',ENRG1,' ENRGP=',ENRGP
         print *,'ENRGEV=',EVAP*(TG1*SHI-LHM)
         print *,'F0DT=',F0DT,' F1DT=',F1DT
      endif
C PRE
      IF(ENRG1.LE.-TG1*HC1) GO TO 130                                   7109.5  
C**** FLUXES HEAT UP TG1 TO FREEZING POINT AND MELT SOME SNOW AND ICE   7110.   
      RUN0=(ENRG1+TG1*HC1)/LHM                                          7110.5  
      TG1=0.                                                            7111.   
      SNANDI=SNANDI-RUN0                                                7111.5  
         CRUN0=CRUN0+RUN0*POICE                                         7112.   
         RUN0S=RUN0S+RUN0*POICE                                         7112.5  
      GO TO 140                                                         7113.   
C**** FLUXES RECOMPUTE TG1 WHICH IS BELOW FREEZING POINT                7113.5  
  130 TG1=TG1+ENRG1/HC1                                                 7114.   
  140 IF(SNANDI.GE.ACE1I) GO TO 160                                     7114.5  
C**** SOME ICE HAS MELTED OR EVAPORATED, TAKE IT FROM G2                7115.   
      SNOW=0.                                                           7115.5  
      DIFS=SNANDI-ACE1I                                                 7116.   
      TG1=(TG1*SNANDI-TG2*DIFS)/ACE1I                                   7116.5  
      EDIFS=DIFS*(TG2*SHI-LHM)                                          7117.   
      IF(KOCEAN.EQ.1) GO TO 150                                         7117.5  
         CEDIFS=CEDIFS+EDIFS*POICE                                      7118.   
         CDIFS=CDIFS+DIFS*POICE                                         7118.5  
         CERUN2=CERUN2+EDIFS*POICE                                      7119.   
         CRUN2=CRUN2+DIFS*POICE                                         7119.5  
         DIFSS=DIFSS+DIFS*POICE                                         7120.   
      GO TO 200                                                         7120.5  
  150 ACE2=ACE2+(DIFS+ATRUNC)                                           7121.   
      DIFSI=ROICE*DIFS                                                  7121.5  
      EDIFSI=ROICE*EDIFS                                                7122.   
      GO TO 210                                                         7122.5  
  160 SNOW=SNANDI-ACE1I                                                 7123.   
C FROM PREC
      if(SNOW.GT.ACE1I)then
C****  SNOW IS COMPACTED INTO ICE, ICE FROM LAYER 1 GOES DOWN TO LAYER 2
       DIFS=SNOW-.9*ACE1I                                               4095.5
       SNOW=.9*ACE1I                                                    4096.
       EDIFS=DIFS*(TG1*SHI-LHM)                                         4096.5
       IF(KOCEAN.EQ.0) THEN
         ERUN2=DIFS*(TG2*SHI-LHM)                                       4097.5
          CEDIFS=CEDIFS+EDIFS*POICE                                     4098.5
          CDIFS=CDIFS+DIFS*POICE                                        4099.
          DIFSS=DIFSS+DIFS*POICE                                        4099.5
          CERUN2=CERUN2+ERUN2*POICE                                     4100.
          CRUN2=CRUN2+DIFS*POICE                                        4100.5
       ENDIF
C****  DIFFUSION CHANGES ICE AMOUNT AND TEMPERATURE OF SECOND LAYER     4101.5
       TG2=TG2+(TG1-TG2)*DIFS/(ACE2+DIFS)
      endif
C FROM PREC
      IF(KOCEAN.EQ.1) GO TO 210                                         7123.5  
C****                                                                   7124.   
C**** OCEAN ICE, CALCULATE TG2                                          7124.5  
C****                                                                   7125.   
  200 F2DT=DTSRCE*(TG2-TFO)*ALAMI*2./Z2
      TG2=TG2+(F1DT-F2DT)/(ACE2*SHI)                                    7126.   
      GO TO 370                                                         7126.5  
  210 F2DT=DTSRCE*(TG2-TGW)*ALAMI*2./Z2
      ENRG2=F1DT-F2DT                                                   7127.5  
      ENRGIW=F2DT+OTDT-ERUN4                                            7128.   
      ENRGFI=0.                                                         7128.5  
      HC2=ACE2*SHI                                                      7129.   
      IF(ENRG2.LE.-TG2*HC2) GO TO 220                                   7129.5  
C**** FLUXES HEAT UP TG2 TO FREEZING POINT AND MELT SOME ICE            7130.   
      ACE2M=(ENRG2+TG2*HC2)/LHM                                         7130.5  
      TG2=0.                                                            7131.   
      ACE2=ACE2+(-ACE2M+ATRUNC)                                         7131.5  
         AIFI=AIFI+ACE2M*POICE                                          7132.   
      GO TO 300                                                         7132.5  
C**** CALCULATE THE ENERGY OF THE WATER BELOW THE ICE AT THE FREEZING   7133.   
C****   POINT AND TEST WHETHER NEW ICE MUST BE FORMED                   7133.5  
  220 WTRI1=WTRO-(SNOW+ACE1I+ACE2)                                      7134.   
      EFIW=WTRI1*TFO*SHW                                                7134.5  
      IF(EIW0+ENRGIW.LT.EFIW) GO TO 240                                 7135.   
C**** FLUXES RECOMPUTE TG2, THE WATER BELOW THE ICE IS ABOVE FREEZING   7135.5  
      TG2=TG2+ENRG2/HC2                                                 7136.   
      GO TO 300                                                         7136.5  
C**** FLUXES WOULD COOL TGIW TO BELOW FREEZING, FREEZE SOME MORE ICE    7137.   
  240 ACE2F=(EIW0+ENRGIW-EFIW)/(TFO*(SHI-SHW)-LHM)                      7137.5  
      ENRGFI=ACE2F*(TFO*SHI-LHM)                                        7138.   
         CEFI=CEFI-ENRGFI*POICE                                         7138.5  
         CIFI=CIFI-ACE2F*POICE                                          7139.   
      ACE2=ACE2+(ACE2F+ATRUNC)                                          7139.5  
      TG2=TG2+(ENRG2+ACE2F*(TFO-TG2)*SHI)/(ACE2*SHI)                    7140.   
C****                                                                   7140.5  
C**** CALCULATE COMPOSITE TEMPERATURES AND REDISTRIBUTION OF ICE        7141.   
C****                                                                   7141.5  
  300 WTRW=WTRW0-(1.-ROICE)*ACEFO+ROICE*(RUN0-RUN4+ACE2M-ACE2F)         7142.   
      ENRGW=ENRGW0+(1.-ROICE)*(ENRGO-ENRGFO)+ROICE*(ENRGIW-ENRGFI)      7142.5  
      TGW=ENRGW/(WTRW*SHW)+TTRUNC                                       7143.   
      if(J.eq.-3.and.TAU.ge. 17520.00)then
         print *,' ENRGW=',ENRGW
         print *,' ENRGW0=',ENRGW0
         print *,' ENRGO=',ENRGO
         print *,' ENRGFO=',ENRGFO
         print *,' ENRGIW=',ENRGIW
         print *,' ENRGFI=',ENRGFI
         print *,' ROICE=',ROICE
      endif
      IF(ACEFO.LE.0.) GO TO 310                                         7143.5  
C**** NEW ICE FORMED ON THE OCEAN SURFACE                               7144.   
      DRO=(1.-ROICE)*ACEFO/(ACE1I+AC2OIM)                               7144.5  
      TG1=TG1+(TFO-TG1)*DRO*ACE1I/(ROICE*(SNOW+ACE1I)+DRO*ACE1I)        7145.   
      TG2=TG2+(TFO-TG2)*DRO*AC2OIM/(ROICE*ACE2+DRO*AC2OIM)              7145.5  
      SNOW=SNOW*ROICE/(ROICE+DRO)                                       7146.   
      ROICE=ROICE+DRO                                                   7146.5  
      ACE2=ACE2+(DRO*(AC2OIM-ACE2)/ROICE+ATRUNC)                        7147.   
  310 IF(ACE2.GE.AC2OIM) GO TO 320                                      7147.5  
C**** ICE IS TO THIN, COMPRESS IT HORIZONTALLY                          7148.   
      ROICEN=ROICE*(ACE1I+ACE2)/(ACE1I+AC2OIM)                          7148.5  
      GO TO 340                                                         7149.   
  320 OPNOCN=.06*(RHOI/(ACE1I+ACE2)-BYZICX)                             7149.5  
      IF(1.-ROICE.GT.OPNOCN) GO TO 360                                  7149.6  
C**** TOO LITTLE OPEN OCEAN, COMPRESS THE ICE HORIZONTALLY              7150.   
      ROICEN=1.-OPNOCN                                                  7150.5  
  340 DRI=ROICE-ROICEN                                                  7151.   
      DIFS=DRI*ACE1I/ROICE                                              7151.5  
      SNOW=SNOW*(ROICE/ROICEN)                                          7152.   
      TG2=TG2+(TG1-TG2)*DIFS/(ACE2+DIFS)                                7152.5  
      ACE2=ACE2+(DRI*(ACE1I+ACE2)/ROICEN+ATRUNC)                        7153.   
      EDIFSI=EDIFSI+ROICE*DIFS*(TG1*SHI-LHM)                            7153.5  
      DIFSI=DIFSI+ROICE*DIFS                                            7154.   
      ROICE=ROICEN                                                      7154.5  
C**** RESAVE PROGNOSTIC QUANTITIES                                      7155.   
  360 ODATA(I,J,1)=TGW                                                  7155.5  
      if(J.eq.-3.and.TAU.ge. 17520.00)then
         print *,'3 ODATA(i,j,1)=',ODATA(i,j,1)
      endif
      if(ROICE.lt.0.0.or.ACE2.lt.AC2OIM)then
        print *,'From ground'
        print *,'J=',' ROICE=',ROICE,' ACE2=',ACE2
        stop
      endif
      if(J.eq.-3.and.TAU.ge. 17520.00)then
         print *,'after 360 ROICE==',ROICE
         print *,' TG1=',TG1,' TG2=',TG2
      endif
      ODATA(I,J,2)=ROICE                                                7156.   
      ODATA(I,J,3)=ACE2                                                 7156.5  
         COTDT=COTDT+OTDT*POICE                                         7157.   
         CEDIFS=CEDIFS+EDIFSI*PWATER                                    7157.5  
         CDIFS=CDIFS+DIFSI*PWATER                                       7158.   
         DIFSS=DIFSS+DIFSI*PWATER                                       7158.5  
  370 GDATA(I,J,1)=SNOW                                                 7159.   
      GDATA(I,J,3)=TG1                                                  7159.5  
      GDATA(I,J,7)=TG2                                                  7160.   
         CSNOW=CSNOW+SNOW*POICE                                         7160.5  
         CTG1=CTG1+TG1*POICE                                            7161.   
         CTG2=CTG2+TG2*POICE                                            7161.5  
         CACE2=CACE2+ACE2*POICE                                         7162.   
         CF1DT=CF1DT+F1DT*POICE                                         7162.5  
         CF2DT=CF2DT+F2DT*POICE                                         7163.   
         CEVAP=CEVAP+EVAP*POICE                                         7163.5  
         CICOV=CICOV+POICE                                              7164.   
         SNOWS=SNOWS+SNOW*POICE                                         7164.5  
         TG1S=TG1S+TG1*POICE                                            7165.   
         ACE1S=ACE1S+ACE1I*POICE                                        7165.5  
         ACE2S=ACE2S+ACE2*POICE                                         7166.   
         TG2S=TG2S+TG2*POICE                                            7166.5  
         EVAPS=EVAPS+EVAP*POICE                                         7167.   
C****                                                                   7167.5  
  400 continue
         AIJ(I,J,7)=AIJ(I,J,7)+(WTR1+ACE1)/WFC1                         7305.5  
         AIJ(I,J,50)=AIJ(I,J,50)+(WTR1+ACE1+WTR2+ACE2)                  7306.   
C****                                                                   7306.5  
C**** ACCUMULATE DIAGNOSTICS                                            7307.   
C****                                                                   7307.5  
C**** QUANTITIES ACCUMULATED FOR LATITUDE-LONGITUDE MAPS IN DIAGIJ      7314.   
  950    AIJ(I,J,6)=AIJ(I,J,6)+EVAPS                                    7314.5  
         AIJ(I,J,28)=AIJ(I,J,28)+TG1S                                   7315.   
  960 CONTINUE                                                          7315.5  
C**** LONGITUDINALLY INTEGRATED QUANTITIES FOR DIAGJ                    7316.   
         CJ(J,15)=CJ(J,15)+CF2DT                                        7316.5  
         AJ(J,17)=AJ(J,17)+ATG2                                         7317.   
         CJ(J,17)=CJ(J,17)+CTG2                                         7318.   
         AJ(J,18)=AJ(J,18)+ATG1                                         7318.5  
         CJ(J,18)=CJ(J,18)+CTG1                                         7319.5  
         AJ(J,19)=AJ(J,19)+AEVAP                                        7320.   
         CJ(J,19)=CJ(J,19)+CEVAP                                        7321.   
         CJ(J,30)=CJ(J,30)+CICOV                                        7321.5  
CCC      AJ(J,33)=AJ(J,33)+AOTDT                                        7322.   
CCC      CJ(J,33)=CJ(J,33)+COTDT                                        7322.5  
         AJ(J,34)=AJ(J,34)+ATG3                                         7323.   
         CJ(J,41)=CJ(J,41)+CEDIFS                                       7324.5  
         CJ(J,42)=CJ(J,42)+CF1DT                                        7325.5  
         AJ(J,43)=AJ(J,43)+AEFO                                         7326.   
         CJ(J,43)=CJ(J,43)+(CERUN2+CEFI)                                7327.   
         CJ(J,45)=CJ(J,45)+CDIFS                                        7328.   
         AJ(J,46)=AJ(J,46)+AIFO                                         7328.5  
         CJ(J,46)=CJ(J,46)+(CRUN2+CIFI)                                 7329.5  
         AJ(J,47)=AJ(J,47)+ARUN4                                        7330.   
         CJ(J,47)=CJ(J,47)+CRUN4                                        7330.5  
         AJ(J,48)=AJ(J,48)+AERUN4                                       7331.   
         CJ(J,48)=CJ(J,48)+CERUN4                                       7331.5  
         CJ(J,52)=CJ(J,52)+CACE2                                        7334.   
         CJ(J,53)=CJ(J,53)+CSNOW                                        7335.   
         CJ(J,54)=CJ(J,54)+CRUN0                                        7336.   
c     if(J.eq.8.and.TAU.ge. 17520.00)then
c        print *,'4 ODATA(1,j,1)=',ODATA(1,j,1)
c     endif
  980 CONTINUE                                                          7336.5  
      RETURN                                                            7337.   
      END                                                               7337.5  
