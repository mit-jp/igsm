C $Header: /home/igsm-repository/ATM/gr_clm.F,v 1.1 2010/08/30 15:36:35 jscott Exp $
C $Name:  $

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

      SUBROUTINE GR_CLM
C****                                                                   7001.5  
C**** THIS SUBROUTINE USES THE SURFACE FLUXES TO PREDICT IN TIME THE    7002.   
C**** GROUND TEMPERATURE, GROUND WATER AND ICE, AND SNOW MELTING.       7002.5  
C****                                                                   7003.   

#if ( defined CLM )
#include "BD2G04.COM"                                      		7003.5  
#if ( defined OCEAN_3D || defined ML_2D )
#include "AGRID.h"
#endif
#include "CLM.h"
      dimension DZSOI(nlevsoi)


      COMMON U,V,T,P,Q                                                  7004.   
      COMMON/WORK1/CONV(IM0,JM0,LM0),PK(IM0,JM0,LM0),PREC(IM0,JM0),
     &   TPREC(IM0,JM0)
      COMMON/WORK3/E0(IO0,JM0,4),E1(IO0,JM0,4),EVAPOR(IO0,JM0,4)        7004.5  
     *  ,TGRND(IO0,JM0,4)
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
      DATA DZSOI / 0.0175128179162552, 
     &  0.0275789692596763,  
     &  0.0454700332424132, 
     &  0.0749674109862086,
     &  0.123600365102281,
     &  0.203782551010432,  
     &  0.335980626448433, 
     &  0.553938405368685,
     &  0.913290031589061,
     &  1.1369718045125 /
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
      ZS1=0.0
      ZS2=0.0
      ZS3=0.0
      do k=1,nlevsoi
       if(k.le.3)then
       ZS1=ZS1+DZSOI(k)
       else if(k.le.8)then
       ZS2=ZS2+DZSOI(k)
       endif
       ZS3=ZS3+DZSOI(k)
      enddo
      print *,'Soil depth ',ZS1,ZS2,ZS3
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
         BF1DT=0.                                                       7032.   
         BEDIFS=0.                                                      7035.   
         BERUN0=0.                                                      7036.   
         BERUN2=0.                                                      7037.   
         BTG1=0.                                                        7039.5  
         BTG2=0.                                                        7041.   
         BEVAP=0.                                                       7043.   
         BDIFS=0.                                                       7044.   
         BRUN0=0.                                                       7046.   
	 BRUNS0=0.
         BRUN2=0.                                                       7047.   
         BWTR1=0.                                                       7049.   
         BACE1=0.                                                       7049.5  
         BWTR2=0.                                                       7050.   
         BACE2=0.                                                       7050.5  
         BSNOW=0.                                                       7051.5  
         BSCOV=0.
         BLWUP=0.
         BSWUP=0.
      DO 960 I=1,IMAX                                                   7053.   
C****                                                                   7053.5  
C****                                                                   7054.5  
C     assign CLM results to GDATA
      TG11=0.0
      WTR1=0.0
      ACE1=0.0
      do k=1,3
       TG11=TG11+tsoiclm(i,j,k)*DZSOI(k)
c      WTR1=WTR1+h2olclm(i,j,k)*DZSOI(k)
c      ACE1=ACE1+h2oiclm(i,j,k)*DZSOI(k)
       WTR1=WTR1+h2olclm(i,j,k)
       ACE1=ACE1+h2oiclm(i,j,k)
      enddo
      TG11=TG11/ZS1
c     WTR1=0.1*WTR1/ZS1
c     ACE1=0.1*ACE1/ZS1
      TG2=0.0
      WTR2=0.0
      ACE2=0.0
      do k=4,8
       TG2=TG2+tsoiclm(i,j,k)*DZSOI(k)
c      WTR2=WTR2+h2olclm(i,j,k)*DZSOI(k)
c      ACE2=ACE2+h2oiclm(i,j,k)*DZSOI(k)
       WTR2=WTR2+h2olclm(i,j,k)
       ACE2=ACE2+h2oiclm(i,j,k)
      enddo
      TG2=TG2/ZS2-273.16
      WTR2=WTR2/ZS2
      ACE2=ACE2/ZS2
c     TG1=tgndclm(i,j)-273.16
      TG1=TGRND(I,J,3)
      SNOW=100.*snwdclm(i,j)
C        snwdclm in meters
      SCOVL=snwcclm(i,j)
      RUNS0=DTSRCE*(sroclm(i,j)+ssrclm(i,j)+glrclm(i,j))
!     RUN0=DTSRCE*(sroclm(i,j)+glrclm(i,j))
      RUN0=DTSRCE*+glrclm(i,j)
      RUN2=DTSRCE*ssrclm(i,j)
c     print *,'From gr_clm J=',j,' TAU=',TAU
c     print *,tsoiclm(i,j,1),tsoiclm(i,j,nlevsoi)
c     print *,TG1,TG2,TG11
c     print *,WTR1,WTR2,SNOW
      albinr=0.7*aldirclm(i,j)+0.3*aldifclm(i,j)
      albvis=0.7*asdirclm(i,j)+0.3*asdifclm(i,j)
      swuclm=albinr*swinr4clm(i,j)+albvis*swvis4clm(i,j)
!        if(J.eq.29)then
!          if(TAU.ge.3644.0.and.TAU.le.3743.0) then
!         print *,'Form gr_clm J=',j,' TAU=',TAU
!         print *,' albinr=',albinr,' albvis=',albvis
!         print *,' swuclm=',swuclm
!         endif
!        endif

      GDATA(I,J,2)=SNOW
      GDATA(I,J,4)=TG1                                                  7294.5  
      GDATA(I,J,5)=WTR1                                                 7295.   
      GDATA(I,J,6)=ACE1                                                 7295.5  
      GDATA(I,J,8)=TG2                                                  7296.   
      GDATA(I,J,9)=WTR2                                                 7296.5  
      GDATA(I,J,10)=ACE2                                                7297.   
      GDATA(I,J,12)=SNOW
      GDATA(I,J,13)=TG1                                                 7191.   
      GDATA(I,J,14)=TG2                                                 7191.5  

      PLAND=FDATA(I,J,2)                                                7055.   
C****                                                                   7167.5  
      F0DT=E0(I,J,3)                                                    7171.5  
      F1DT=E1(I,J,3)                                                    7172.5  
      EVAP=EVAPOR(I,J,3)                                                7173.   


C**** RESAVE PROGNOSTIC QUANTITIES                                      7293.5  
         BSNOW=BSNOW+SNOW*PLAND
         BSCOV=BSCOV+SCOVL*PLAND
         BTG1=BTG1+TG1*PLAND
         BTG2=BTG2+TG2*PLAND
         BWTR1=BWTR1+WTR1*PLAND
         BACE1=BACE1+ACE1*PLAND
         BWTR2=BWTR2+WTR2*PLAND
         BACE2=BACE2+ACE2*PLAND
         BF1DT=BF1DT+F1DT*PLAND
         BEVAP=BEVAP+EVAP*PLAND
         BRUNS0=BRUNS0+RUNS0*PLAND
         BRUN0=BRUN0+RUN0*PLAND
         BRUN2=BRUN2+RUN2*PLAND
         BLWUP=BLWUP+lwuclm(i,j)*PLAND
         BSWUP=BSWUP+swuclm*PLAND
C****                                                                   7306.5  
C**** ACCUMULATE DIAGNOSTICS                                            7307.   
C****                                                                   7307.5  
  960 CONTINUE                                                          7315.5  
C**** LONGITUDINALLY INTEGRATED QUANTITIES FOR DIAGJ                    7316.   
         BJ(J,17)=BJ(J,17)+BTG2                                         7317.5  
         BJ(J,18)=BJ(J,18)+BTG1                                         7319.   
         BJ(J,19)=BJ(J,19)+BEVAP                                        7320.5  
c        BJ(J,40)=BJ(J,40)+BERUN0                                       7323.5  
c        BJ(J,41)=BJ(J,41)+BEDIFS                                       7324.   
!        BJ(J,42)=BJ(J,42)+BF1DT                                        7325.   
c        BJ(J,43)=BJ(J,43)+BERUN2                                       7326.5  
         BJ(J,43)=BJ(J,43)+BSWUP
c        BJ(J,45)=BJ(J,45)+BDIFS                                        7327.5  
         BJ(J,45)=BJ(J,45)+BLWUP
         BJ(J,49)=BJ(J,49)+BWTR1                                        7332.   
         BJ(J,50)=BJ(J,50)+BACE1                                        7332.5  
         BJ(J,51)=BJ(J,51)+BWTR2                                        7333.   
         BJ(J,52)=BJ(J,52)+BACE2                                        7333.5  
         BJ(J,53)=BJ(J,53)+BSNOW                                        7334.5  
c        moved from dd2g04.F
         BJ(J,31)=BJ(J,31)+BSCOV
C     From prland
c        BJ(J,39)=BJ(J,39)+BENRGP                                       4200.
c        BJ(J,40)=BJ(J,40)+BERUN0                                       4201.
c        BJ(J,41)=BJ(J,41)+BEDIFS                                       4201.5
c        BJ(J,43)=BJ(J,43)+BERUN2                                       4203.
c        BJ(J,45)=BJ(J,45)+BDIFS                                        4204.
C        Runoff from first layer of soil including ice melting
         BJ(J,47)=BJ(J,47)+BRUN0
C        Runoff from first layer of soil does not include ice melting
         BJ(J,54)=BJ(J,54)+BRUNS0
         BJ(J,46)=BJ(J,46)+BRUN2                                        
C     From prland
#if ( defined OCEAN_3D || defined ML_2D )
COO	 Runoff from first layer of soil does not include ice melting
C	 Runoff from first layer of soil including glacial ice melting
           arunoff(j)=arunoff(j)+RUNS0
#endif
  980 CONTINUE                                                          7336.5  
#endif
      RETURN                                                            7337.   
      END                                                               7337.5  
