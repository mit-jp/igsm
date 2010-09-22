
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

      SUBROUTINE GRLAND                                                 7001.   
C****                                                                   7001.5  
C**** THIS SUBROUTINE USES THE SURFACE FLUXES TO PREDICT IN TIME THE    7002.   
C**** GROUND TEMPERATURE, GROUND WATER AND ICE, AND SNOW MELTING.       7002.5  
C****                                                                   7003.   

#include "BD2G04.COM"                                      		7003.5  

#if ( defined OCEAN_3D || defined ML_2D )
#include "AGRID.h"
#endif

      COMMON U,V,T,P,Q                                                  7004.   
      COMMON/WORK1/CONV(IM0,JM0,LM0),PK(IM0,JM0,LM0),PREC(IM0,JM0),
     &   TPREC(IM0,JM0)
      COMMON/WORK3/E0(IO0,JM0,4),E1(IO0,JM0,4),EVAPOR(IO0,JM0,4)        7004.5  
      COMMON/OT/OTA(IO0,JM0),OTB(IO0,JM0),OTC(IO0,JM0)                  7005.   
      COMMON/SPEC2/KM,KINC,COEK,C3LAND(IO0,JM0),C3OICE(IO0,JM0)         7005.1  
     *  ,C3LICE(IO0,JM0),WMGE(IO0,JM0),TSSFC(1,JM0,4)                   7005.2  
      common/qfl/QFLUX(JM0,0:13),ZOAV(JM0),QFLUXT(JM0)
      COMMON/OLDZO/Z1OOLD(IO0,JM0)
      DIMENSION FWATER(JM0),TOBS(JM0),TF68(JM0),DELR(JM0)               7005.3  
      COMMON/FRMIC/ FRMDICE(JM0)
      DATA SHV/0./,SHW/4185./,SHI/2060./,RHOW/1000./,RHOI/916.6/,       7005.9  
     *  ALAMI/2.1762/,TFO/-1.56/,Z1I/.1/,Z2LI/2.9/,Z1E/.1/,Z2E/4./      7006.   
      INTEGER JDOFM(13)                                                 7006.1  
      DATA JDOFM/0,31,59,90,120,151,181,212,243,273,304,334,365/        7006.2  
      DATA Z2OIM/0.9/,Z2OIX/4.9/                                        7006.5  
      DATA TTRUNC/0./                                                   7007.   
      DATA IFIRST/1/                                                    7007.5  
      ROSNOW(X)=0.54*X/LOG(1.+0.54*X/275.)
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
      print *,'E0(1,j,1)'
      print *,(E0(1,j,1),j=1,JM0)
      COEFSN=1.
      print *,' COEFSN=',COEFSN
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
      IMAX=IM                                                           7031.   
      IF((J.EQ.1).OR.(J.EQ.JM)) IMAX=1                                  7031.5  
C     DELR(J)=(TF68(J)-TOBS(J))*Z1O(1,J)*RHOW*SHW/(365.*SDAY)           7031.6  
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
	 BRUNS0=0.
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
         RUNS0=0.                                                       7063.5  
         DIFSS=0.                                                       7064.   
C****                                                                   7064.5  
C****                                                                   7167.5  
  400 IF(PLICE.LE.1.E-5) GO TO 600                                      7168.   
C****                                                                   7168.5  
C**** LAND ICE                                                          7169.   
C****                                                                   7169.5  
      SNOW=GDATA(I,J,12)                                                7170.   
      TG1=GDATA(I,J,13)                                                 7170.5  
      TG2=GDATA(I,J,14)                                                 7171.   
      F0DT=E0(I,J,3)                                                    7171.5  
         AIJ(I,J,67)=AIJ(I,J,67)+F0DT                                   7172.   
      F1DT=E1(I,J,3)                                                    7172.5  
      EVAP=EVAPOR(I,J,3)                                                7173.   
         AIJ(I,J,63)=AIJ(I,J,63)+EVAP                                   7173.5  
C**** CALCULATE TG1                                                     7174.   
      SNANDI=SNOW+ACE1I-EVAP                                            7174.5  
c     if (SNANDI.le.0.0)then
c       print *,' After 7174.5 TAU=',TAU,' J=',J
c     endif
      HC1=SNANDI*SHI                                                    7175.   
      ENRG1=F0DT+EVAP*(TG1*SHI-LHM)-F1DT                                7175.5  
      IF(ENRG1.LE.-TG1*HC1) GO TO 420                                   7176.   
C**** FLUXES HEAT UP TG1 TO FREEZING POINT AND MELT SOME SNOW AND ICE   7176.5  
      RUN0=(ENRG1+TG1*HC1)/LHM                                          7177.   
c	RUNS0 does not include runoff due melting of land ice
c	RUNS0 includes runoff due melting of land ice formed from snow
c       during the run.
	RUNS0=DMIN1(RUN0,SNOW+FRMDICE(J))
        if(RUN0.gt.SNOW)then
          FRMDICE(j)=FRMDICE(j)-(RUN0-SNOW)
          if( FRMDICE(j).lt.0.0) FRMDICE(j)=0.0
        endif
c
      TG1=0.                                                            7177.5  
      SNANDI=SNANDI-RUN0                                                7178.   
c     if (SNANDI.le.0.0)then
c       print *,' After 7178 TAU=',TAU,' J=',J
c     endif
         BRUN0=BRUN0+RUN0*PLICE                                         7178.5  
	 BRUNS0=BRUNS0+RUNS0*PLICE
         RUN0S=RUN0S+RUN0*PLICE                                         7179.   
         AIJ(I,J,33)=AIJ(I,J,33)+RUN0                                   7179.5  
      GO TO 440                                                         7180.   
C**** FLUXES RECOMPUTE TG1 WHICH IS BELOW FREEZING POINT                7180.5  
  420 TG1=TG1+ENRG1/HC1                                                 7181.   
  440 IF(SNANDI.GE.ACE1I) GO TO 460                                     7181.5  
C**** SOME ICE HAS MELTED OR EVAPORATED, TAKE IT FROM G2                7182.   
      SNOW=0.                                                           7182.5  
      DIFS=SNANDI-ACE1I                                                 7183.   
      TG1=(TG1*SNANDI-TG2*DIFS)/ACE1I                                   7183.5  
c     if (DIFS.le.0.0)then
c       print *,' After 7183.5 TAU=',TAU,' J=',J
c       print *,'SNANDI=',SNANDI,' DIFS=',DIFS
c       print *,'TG1=',TG1,' TG2=',TG2
c     endif
      EDIFS=DIFS*(TG2*SHI-LHM)                                          7184.   
         BEDIFS=BEDIFS+EDIFS*PLICE                                      7184.5  
         AIJ(I,J,69)=AIJ(I,J,69)+EDIFS                                  7185.   
         BDIFS=BDIFS+DIFS*PLICE                                         7185.5  
         DIFSS=DIFSS+DIFS*PLICE                                         7186.   
         BERUN2=BERUN2+EDIFS*PLICE                                      7186.5  
C        AIJ(I,J,72)=AIJ(I,J,72)+ERUN2                                  7187.   
         AIJ(I,J,72)=AIJ(I,J,72)+EDIFS                                  7187.1  
         BRUN2=BRUN2+DIFS*PLICE                                         7187.5  
      GO TO 500                                                         7188.   
  460 SNOW=SNANDI-ACE1I                                                 7188.5  
C**** CALCULATE TG2                                                     7189.   
  500 TG2=TG2+F1DT/HC2LI                                                7189.5  
C**** RESAVE PROGNOSTIC QUANTITIES                                      7190.   
      GDATA(I,J,12)=SNOW                                                7190.5  
      GDATA(I,J,13)=TG1                                                 7191.   
      GDATA(I,J,14)=TG2                                                 7191.5  
         if (SNOW.gt.10.)then
            RHOS0=ROSNOW(SNOW)
         else
            RHOS0=275.
         endif
         RHOS=COEFSN*RHOS0
         BSNOW=BSNOW+100*SNOW/RHOS*PLICE  
c        BSNOW=BSNOW+SNOW*PLICE                                         7192.   
         BTG1=BTG1+TG1*PLICE                                            7192.5  
         BTG2=BTG2+TG2*PLICE                                            7193.   
         BF1DT=BF1DT+F1DT*PLICE                                         7193.5  
         AIJ(I,J,69)=AIJ(I,J,69)+F1DT                                   7194.   
         BEVAP=BEVAP+EVAP*PLICE                                         7194.5  
         SNOWS=SNOWS+SNOW*PLICE                                         7195.   
         TG1S=TG1S+TG1*PLICE                                            7195.5  
         ACE1S=ACE1S+ACE1I*PLICE                                        7196.   
         ACE2S=ACE2S+Z2LI*RHOI*PLICE                                    7196.5  
         TG2S=TG2S+TG2*PLICE                                            7197.   
         EVAPS=EVAPS+EVAP*PLICE                                         7197.5  
C****                                                                   7198.   
  600 IF(PEARTH.LE.1.E-5) GO TO 940                                     7198.5  
C****                                                                   7199.   
C**** EARTH                                                             7199.5  
C****                                                                   7200.   
      SNOW=GDATA(I,J,2)                                                 7200.5  
      TG1=GDATA(I,J,4)                                                  7201.   
      WTR1=GDATA(I,J,5)                                                 7201.5  
      ACE1=GDATA(I,J,6)                                                 7202.   
      TG2=GDATA(I,J,8)                                                  7202.5  
      WTR2=GDATA(I,J,9)                                                 7203.   
      ACE2=GDATA(I,J,10)                                                7203.5  
      WFC1=VDATA(I,J,9)                                                 7204.   
      WFC2=VDATA(I,J,10)                                                7204.5  
      CHI1=(WTR1+ACE1)/WFC1                                             7205.   
      GFAC=WFC2/WFC1                                                    7205.5  
      HC1=HC1DE+WTR1*SHW+(ACE1+SNOW)*SHI                                7206.   
      F0DT=E0(I,J,4)                                                    7206.5  
         AIJ(I,J,68)=AIJ(I,J,68)+F0DT                                   7207.   
      F1DT=E1(I,J,4)                                                    7207.5  
      EVAP=EVAPOR(I,J,4)                                                7208.   
         EVAPS=EVAPS+EVAP*PEARTH                                        7208.5  
         AIJ(I,J,64)=AIJ(I,J,64)+EVAP                                   7209.   
C****                                                                   7209.5  
C**** EARTH, DETERMINE EVAPORATION AND DIFFUSION OF WATER               7210.   
C****                                                                   7210.5  
      ENRG1=0.                                                          7211.   
      IF(SNOW.LE.0.) GO TO 625                                          7211.5  
      IF(EVAP.GT.SNOW) GO TO 620                                        7212.   
C**** SOME SNOW EVAPORATES                                              7212.5  
  610 ENRG1=EVAP*(TG1*SHI-LHM)                                          7213.   
      SNOW=SNOW-EVAP                                                    7213.5  
      GO TO 660                                                         7214.   
C**** ALL SNOW EVAPORATES                                               7214.5  
  620 EVAP=EVAP-SNOW                                                    7215.   
      ENRG1=SNOW*(TG1*SHI-LHM)                                          7215.5  
      SNOW=0.                                                           7216.   
      GO TO 630                                                         7216.5  
  625 IF(WTR1+ACE1-EVAP.LE.WFC1) GO TO 630                              7217.   
C**** DEW+WTR+ACE EXCEEDS WFC, DO NOT BOTHER TO CORRECT ALL DIAGNOSTICS 7217.5  
C     IF(TG1.LE.0.) GO TO 610                                           7218.   
C     RUN0=WTR1+ACE1-EVAP-WFC1                                          7218.5  
C     ERUN0=RUN0*TG1*SHW                                                7219.   
C     WTR1=WFC1-ACE1                                                    7219.5  
C     ENRG1=-ERUN0                                                      7220.   
      GO TO 660                                                         7220.5  
  630 IF(EVAP.LE.WTR1+ACE1) GO TO 640                                   7221.   
C**** ALL WATER AND ICE EVAPORATES,RECOMPUTE EVHDT AND EVAP             7221.5  
      DEVAP=EVAP-(WTR1+ACE1)                                            7222.   
      DEVHDT=-DEVAP*(LHE+TG1*SHV)                                       7222.5  
      ENRG1=ENRG1+WTR1*TG1*SHW+ACE1*(TG1*SHI-LHM)                       7223.   
      EVAP=WTR1+ACE1                                                    7223.5  
      WTR1=0.                                                           7224.   
      ACE1=0.                                                           7224.5  
      F0DT=F0DT-DEVHDT                                                  7225.   
         BJ(J,14)=BJ(J,14)-DEVHDT*PEARTH                                7225.5  
         AIJ(I,J,23)=AIJ(I,J,23)-DEVHDT*PEARTH                          7226.   
         AIJ(I,J,68)=AIJ(I,J,68)-DEVHDT                                 7226.5  
         AIJ(I,J,64)=AIJ(I,J,64)-DEVAP                                  7227.   
         EVAPS=EVAPS-DEVAP*PEARTH                                       7227.5  
      GO TO 660                                                         7228.   
C**** EVAPORATION FROM WATER AND ICE                                    7228.5  
  640 IF(EVAP.GT.WTR1+ACE1) EVAP=WTR1+ACE1                              7228.6  
      DWET=EVAP/(WTR1+ACE1+1.E-20)                                      7229.   
  650 ENRG1=ENRG1+DWET*(WTR1*TG1*SHW+ACE1*(TG1*SHI-LHM))                7229.5  
      WTR1=WTR1*(1.-DWET)                                               7230.   
      ACE1=ACE1*(1.-DWET)                                               7230.5  
C**** DETERMINE DIFFUSION OF WATER                                      7231.   
  660 X=(WTR1+ACE1)/(WTR2+ACE2+1.E-20)                                  7231.5  
      IF(1..LT.X*GFAC) GO TO 670                                        7232.   
      GROW=1.                                                           7232.5  
      IF((SINP(J).GT..5).AND.(JDAY-121)*(243-JDAY).LT.0) GROW=0.        7233.   
      IF((SINP(J).LT.-.5).AND.(JDAY-60)*(304-JDAY).GE.0) GROW=0.        7233.5  
      DIFS=GROW*(1.-VDATA(I,J,1))*WTR2*(X*GFAC-1.)/(1.+GFAC)            7234.   
      EDIFS=TG2*DIFS*SHW                                                7234.5  
c     if(JDAY.GE.186.and.(J.eq.43.or.J.eq.44))then
c      if(JRPR.ne.J) then
c       print *,JDAY,J,SINP(J)
c       print *,WTR1,WTR2,ACE1,ACE2
c       print *,TG1,TG2,GROW,VDATA(I,J,1)
c       print *,DIFS
c      endif
c     endif
      GO TO 690                                                         7235.   
  670 DIFS=DIFFUS*WTR1*(GFAC-1./X)/(1.+GFAC)                            7235.5  
      EDIFS=TG1*DIFS*SHW                                                7236.   
  690    BEDIFS=BEDIFS+EDIFS*PEARTH                                     7236.5  
         BDIFS=BDIFS+DIFS*PEARTH                                        7237.   
         DIFSS=DIFSS+DIFS*PEARTH                                        7237.5  
C****                                                                   7238.   
C**** EARTH, CALCULATE TG1                                              7238.5  
C****                                                                   7239.   
      ENRG1=ENRG1+F0DT-F1DT-EDIFS                                       7239.5  
      IF(TG1) 710,740,750                                               7240.   
C**** FREEZE THE WATER THAT DIFFUSES INTO G1                            7240.5  
  710 ENRG1=ENRG1+DIFS*(TG1*SHI-LHM)                                    7241.   
      ACE1=ACE1-DIFS                                                    7241.5  
      HC1=HC1DE+(ACE1+SNOW)*SHI                                         7242.   
      IF(ENRG1.LE.-TG1*HC1) GO TO 780                                   7242.5  
C**** FLUXES HEAT UP TG1 TO FREEZING POINT                              7243.   
      ENRG1=ENRG1+TG1*HC1                                               7243.5  
      TG1=0.                                                            7244.   
  720 IF(ENRG1.LE.(ACE1+SNOW)*LHM) GO TO 730                            7244.5  
C**** SNOW AND GROUND ICE MELTS, RECOMPUTE TG1                          7245.   
      RUN0=MAX(SNOW*.5*CHI1,SNOW+WTR1+ACE1-WFC1)                        7245.5  
      WTR1=WTR1+ACE1+SNOW-RUN0                                          7246.   
      TG1=(ENRG1-(ACE1+SNOW)*LHM)/(HC1DE+WTR1*SHW)                      7246.5  
      ACE1=0.                                                           7247.   
      SNOW=0.                                                           7247.5  
      GO TO 790                                                         7248.   
C**** SOME SNOW AND GROUND ICE MELTS, TG1 IS AT FREEZING POINT          7248.5  
  730 DWATER=ENRG1/LHM                                                  7249.   
      DSNOW=MIN(SNOW,DWATER)                                            7249.5  
      RUN0=MAX(DSNOW*.5*CHI1,DSNOW+WTR1+ACE1-WFC1)                      7250.   
      SNOW=SNOW-DSNOW                                                   7250.5  
      ACE1=ACE1-(DWATER-DSNOW)                                          7251.   
      WTR1=WTR1+(DWATER-RUN0)                                           7251.5  
      GO TO 790                                                         7252.   
C**** TG1 IS AT FREEZING POINT, SUBTRACT WATER THAT DIFFUSES OUT OF G1  7252.5  
  740 WTR1=WTR1-DIFS                                                    7253.   
      HC1=HC1DE+WTR1*SHW+(ACE1+SNOW)*SHI                                7253.5  
      IF(ENRG1.GT.0.) GO TO 720                                         7254.   
      GO TO 760                                                         7254.5  
C**** THE WATER THAT DIFFUSES OUT OF G1 IS ABOVE THE FREEZING POINT     7255.   
  750 ENRG1=ENRG1+TG1*DIFS*SHW                                          7255.5  
      WTR1=WTR1-DIFS                                                    7256.   
      HC1=HC1DE+WTR1*SHW                                                7256.5  
      IF(-ENRG1.LE.TG1*HC1) GO TO 780                                   7257.   
C**** FLUXES COOL TG1 TO FREEZING POINT                                 7257.5  
      ENRG1=ENRG1+TG1*HC1                                               7258.   
      TG1=0.                                                            7258.5  
  760 IF(-ENRG1.LE.WTR1*LHM) GO TO 770                                  7259.   
C**** GROUND WATER FREEZES, RECOMPUTE TG1                               7259.5  
      ACE1=WTR1+ACE1                                                    7260.   
      TG1=(ENRG1+WTR1*LHM)/(HC1DE+(ACE1+SNOW)*SHI)                      7260.5  
      WTR1=0.                                                           7261.   
      GO TO 800                                                         7261.5  
C**** SOME GROUND WATER FREEZES, TG1 IS AT FREEZING POINT               7262.   
  770 DICE=-ENRG1/LHM                                                   7262.5  
      WTR1=WTR1-DICE                                                    7263.   
      ACE1=ACE1+DICE                                                    7263.5  
      GO TO 800                                                         7264.   
C**** FLUXES DO NOT CAUSE TG1 TO CROSS THE FREEZING POINT               7264.5  
  780 TG1=TG1+ENRG1/HC1                                                 7265.   
      GO TO 800                                                         7265.5  
  790    BRUN0=BRUN0+RUN0*PEARTH                                        7266.   
	 RUNS0=RUN0
	 BRUNS0=BRUNS0+RUNS0*PEARTH
         RUN0S=RUN0S+RUN0*PEARTH                                        7266.5  
         AIJ(I,J,32)=AIJ(I,J,32)+RUN0                                   7267.   
C****                                                                   7267.5  
C**** EARTH, CALCULATE TG2                                              7268.   
C****                                                                   7268.5  
  800 ENRG2=F1DT+EDIFS                                                  7269.   
      HC2=HC2DE+WTR2*SHW+ACE2*SHI                                       7269.5  
      IF(TG2) 810,840,850                                               7270.   
C**** FREEZE THE WATER THAT DIFFUSES AND PERCOLATES INTO G2             7270.5  
  810 ENRG2=ENRG2-DIFS*(TG2*SHI-LHM)                                    7271.   
      ACE2=ACE2+DIFS                                                    7271.5  
      HC2=HC2DE+ACE2*SHI                                                7272.   
      IF(ENRG2.LE.-TG2*HC2) GO TO 880                                   7272.5  
C**** FLUXES HEAT UP TG2 TO FREEZING POINT                              7273.   
      ENRG2=ENRG2+TG2*HC2                                               7273.5  
      TG2=0.                                                            7274.   
  820 IF(ENRG2.LE.ACE2*LHM) GO TO 830                                   7274.5  
C**** GROUND ICE MELTS, RECOMPUTE TG2                                   7275.   
      WTR2=WTR2+ACE2                                                    7275.5  
      TG2=(ENRG2-ACE2*LHM)/(HC2DE+WTR2*SHW)                             7276.   
      ACE2=0.                                                           7276.5  
      GO TO 890                                                         7277.   
C**** SOME GROUND ICE MELTS, TG2 IS AT FREEZING POINT                   7277.5  
  830 DWATER=ENRG2/LHM                                                  7278.   
      WTR2=WTR2+DWATER                                                  7278.5  
      ACE2=ACE2-DWATER                                                  7279.   
      GO TO 890                                                         7279.5  
C**** TG2 IS AT FREEZING POINT, ADD IN WATER THAT DIFFUSES OR PERCOLATE 7280.   
  840 WTR2=WTR2+DIFS                                                    7280.5  
      HC2=HC2DE+WTR2*SHW+ACE2*SHI                                       7281.   
      IF(ENRG2.GT.0.) GO TO 820                                         7281.5  
      GO TO 860                                                         7282.   
C**** WATER THAT DIFFUSES OR PERCOLATES IS ABOVE FREEZING POINT         7282.5  
  850 ENRG2=ENRG2-TG2*DIFS*SHW                                          7283.   
      WTR2=WTR2+DIFS                                                    7283.5  
      HC2=HC2DE+WTR2*SHW                                                7284.   
      IF(-ENRG2.LE.TG2*HC2) GO TO 880                                   7284.5  
C**** FLUXES COOL TG2 TO FREEZING POINT                                 7285.   
      ENRG2=ENRG2+TG2*HC2                                               7285.5  
      TG2=0.                                                            7286.   
  860 IF(-ENRG2.LE.WTR2*LHM) GO TO 870                                  7286.5  
C**** GROUND WATER FREEZES, RECOMPUTE TG2                               7287.   
      ACE2=WTR2+ACE2                                                    7287.5  
      TG2=(ENRG2+WTR2*LHM)/(HC2DE+ACE2*SHI)                             7288.   
      WTR2=0.                                                           7288.5  
      GO TO 890                                                         7289.   
C**** SOME GROUND WATER FREEZES, TG2 IS AT FREEZING POINT               7289.5  
  870 DICE=-ENRG2/LHM                                                   7290.   
      WTR2=WTR2-DICE                                                    7290.5  
      ACE2=ACE2+DICE                                                    7291.   
      GO TO 890                                                         7291.5  
C**** FLUXES DO NOT CAUSE TG2 TO CROSS THE FREEZING POINT               7292.   
  880 TG2=TG2+ENRG2/HC2                                                 7292.5  
  890 CONTINUE                                                          7293.   
C**** RESAVE PROGNOSTIC QUANTITIES                                      7293.5  
      GDATA(I,J,2)=SNOW                                                 7294.   
      GDATA(I,J,4)=TG1                                                  7294.5  
c     print *,'From GRLAND TG1=',TG1
C     13/11/03
c     IF(WTR1+ACE1.GT.WFC1) WTR1=.99999*WTR1                            7294.6  
c     IF(WTR1+ACE1.GT.WFC1) ACE1=.99999*ACE1                            7294.7  
      IF(WTR1+ACE1.GT.WFC1) THEN
       if(TG1.gt.0.0)then
         print *,'From ground WFC1'
         print *,JDAY,J
         print *,WFC2,WTR2,TG2,ACE2
         print *,WFC1,WTR1,TG1,ACE1
        WTRT=WFC1-ACE1
        RUNS0=WTR1-WTRT
        BRUNS0=BRUNS0+RUNS0*PEARTH
        WTR1=max(WTRT,0.0)
       endif
      ENDIF
      GDATA(I,J,5)=WTR1                                                 7295.   
      GDATA(I,J,6)=ACE1                                                 7295.5  
      GDATA(I,J,8)=TG2                                                  7296.   
         JDAYPR=JDAY
c     if((WTR2.lt.1.e-7.and.ACE2.lt.1.e-7)
c    &  .or.(WTR1.lt.1.e-7.and.ACE1.lt.1.e-7) )then
c     if(WTR2.lt.1.e-3*WFC2.and.ACE2.lt.1.e-5)then
c       if(JDAY.ne.JDAYPR.and.JRPR.ne.J)then
c        JDAYPR=JDAY
c        JRPR=J
c        print *,'From ground'
c        print *,JDAY,J
c        print *,WFC2,WTR2,TG2,ACE2
c        print *,WFC1,WTR1,TG1,ACE1
c        print *,'WTR2 was not  changed'
c       endif
c       WTR2=1.e-7
c       WTR2=1.e-3*WFC2
c     endif
      IF(WTR2+ACE2.GT.WFC2) THEN
       if(TG2.gt.0.0)then
         print *,'From ground WFC2'
         print *,JDAY,J
         print *,WFC2,WTR2,TG2,ACE2
         print *,WFC1,WTR1,TG1,ACE1
        WTRT=WFC2-ACE2
        RUNS0=WTR2-WTRT
        BRUNS0=BRUNS0+RUNS0*PEARTH
        WTR2=max(WTRT,0.0)
       endif
      ENDIF
      GDATA(I,J,9)=WTR2                                                 7296.5  
      GDATA(I,J,10)=ACE2                                                7297.   
         if (SNOW.gt.10.)then
            RHOS0=ROSNOW(SNOW)
         else
            RHOS0=275.
         endif
         RHOS=COEFSN*RHOS0
         BSNOW=BSNOW+100*SNOW/RHOS*PEARTH  
c        BSNOW=BSNOW+SNOW*PEARTH                                        7297.5  
         BTG1=BTG1+TG1*PEARTH                                           7298.   
         BTG2=BTG2+TG2*PEARTH                                           7298.5  
         BWTR1=BWTR1+WTR1*PEARTH                                        7299.   
         BACE1=BACE1+ACE1*PEARTH                                        7299.5  
         BWTR2=BWTR2+WTR2*PEARTH                                        7300.   
         BACE2=BACE2+ACE2*PEARTH                                        7300.5  
         BF1DT=BF1DT+F1DT*PEARTH                                        7301.   
         BEVAP=BEVAP+EVAP*PEARTH                                        7301.5  
         SNOWS=SNOWS+SNOW*PEARTH                                        7302.   
         TG1S=TG1S+TG1*PEARTH                                           7302.5  
         WTR1S=WTR1S+WTR1*PEARTH                                        7303.   
         ACE1S=ACE1S+ACE1*PEARTH                                        7303.5  
         WTR2S=WTR2S+WTR2*PEARTH                                        7304.   
         ACE2S=ACE2S+ACE2*PEARTH                                        7304.5  
         TG2S=TG2S+TG2*PEARTH                                           7305.   
         AIJ(I,J,7)=AIJ(I,J,7)+(WTR1+ACE1)/WFC1                         7305.5  
         AIJ(I,J,50)=AIJ(I,J,50)+(WTR1+ACE1+WTR2+ACE2)                  7306.   
C****                                                                   7306.5  
C**** ACCUMULATE DIAGNOSTICS                                            7307.   
C****                                                                   7307.5  
C**** QUANTITIES ACCUMULATED FOR REGIONS IN DIAGJ                       7308.   
  940    IF(JR.EQ.JM0) GO TO 950                                        7308.5  
         DJ(JR,17)=DJ(JR,17)+TG2S*DXYPJ                                 7309.   
         DJ(JR,18)=DJ(JR,18)+TG1S*DXYPJ                                 7309.5  
         DJ(JR,45)=DJ(JR,45)+DIFSS*DXYPJ                                7310.5  
         DJ(JR,49)=DJ(JR,49)+WTR1S*DXYPJ                                7311.   
         DJ(JR,50)=DJ(JR,50)+ACE1S*DXYPJ                                7311.5  
         DJ(JR,51)=DJ(JR,51)+WTR2S*DXYPJ                                7312.   
         DJ(JR,52)=DJ(JR,52)+ACE2S*DXYPJ                                7312.5  
         DJ(JR,53)=DJ(JR,53)+SNOWS*DXYPJ                                7313.   
         DJ(JR,54)=DJ(JR,54)+RUN0S*DXYPJ                                7313.5  
C**** QUANTITIES ACCUMULATED FOR LATITUDE-LONGITUDE MAPS IN DIAGIJ      7314.   
  950    AIJ(I,J,6)=AIJ(I,J,6)+EVAPS                                    7314.5  
         AIJ(I,J,28)=AIJ(I,J,28)+TG1S                                   7315.   
  960 CONTINUE                                                          7315.5  
C**** LONGITUDINALLY INTEGRATED QUANTITIES FOR DIAGJ                    7316.   
         BJ(J,17)=BJ(J,17)+BTG2                                         7317.5  
         BJ(J,18)=BJ(J,18)+BTG1                                         7319.   
         BJ(J,19)=BJ(J,19)+BEVAP                                        7320.5  
         BJ(J,40)=BJ(J,40)+BERUN0                                       7323.5  
         BJ(J,41)=BJ(J,41)+BEDIFS                                       7324.   
         BJ(J,42)=BJ(J,42)+BF1DT                                        7325.   
!        BJ(J,43)=BJ(J,43)+BERUN2                                       7326.5  
!        BJ(J,45)=BJ(J,45)+BDIFS                                        7327.5  
         BJ(J,46)=BJ(J,46)+BRUN2                                        7329.   
         BJ(J,49)=BJ(J,49)+BWTR1                                        7332.   
         BJ(J,50)=BJ(J,50)+BACE1                                        7332.5  
         BJ(J,51)=BJ(J,51)+BWTR2                                        7333.   
         BJ(J,52)=BJ(J,52)+BACE2                                        7333.5  
         BJ(J,53)=BJ(J,53)+BSNOW                                        7334.5  
         BJ(J,47)=BJ(J,47)+BRUN0                                        7335.5  
	 BJ(J,54)=BJ(J,54)+BRUNS0
#if ( defined OCEAN_3D || defined ML_2D )
C	 Runoff from first layer of soil does not include ice melting
	 if(PLICE+PEARTH.gt.0.0)then
           arunoff(j)=arunoff(j)+BRUNS0/(PLICE+PEARTH)
	 end if
#endif
  980 CONTINUE                                                          7336.5  
      RETURN                                                            7337.   
      END                                                               7337.5  
