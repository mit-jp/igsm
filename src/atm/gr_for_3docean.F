
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

      SUBROUTINE GRFOROCEAN                                                7001.   
C****                                                                   7001.5  
C**** THIS SUBROUTINE USES THE SURFACE FLUXES TO PREDICT IN TIME THE    7002.   
C**** GROUND TEMPERATURE, GROUND WATER AND ICE, AND SNOW MELTING.       7002.5  
C****                                                                   7003.   

#include "BD2G04.COM"                                      		7003.5  

      COMMON U,V,T,P,Q                                                  7004.   
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
c     print *,'E0(1,j,1)'
c     print *,(E0(1,j,1),j=1,JM0)
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
c     DO 960 I=1,IMAX                                                   7053.   
       I=1
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
         ATG2=ATG2+ODATA(I,J,4)*POCEAN                                  7070.5  
         TG2S=TG2S+ODATA(I,J,4)*POCEAN                                  7071.   
         ATG3=ATG3+ODATA(I,J,5)*POCEAN                                  7070.5  
      IF(POICE.GT.0.) GO TO 110                                         7071.5  
      GO TO 400                                                         7072.   
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
C****                                                                   7124.   
C****                                                                   7125.   
  200 F2DT=DTSRCE*(TG2-TFO)*ALAMI*2./Z2
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
C****                                                                   7306.5  
C**** ACCUMULATE DIAGNOSTICS                                            7307.   
C****                                                                   7307.5  
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
  980 CONTINUE                                                          7336.5  
      RETURN                                                            7337.   
      END                                                               7337.5  
