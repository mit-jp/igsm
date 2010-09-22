C $Header: /home/igsm-repository/ATM/prland.F,v 1.1 2010/08/30 15:36:35 jscott Exp $
C $Name:  $

#include "ctrparam.h"

      SUBROUTINE PRECIP_LAND(mndriver)                                  4001.   
C****                                                                   4001.5  
C**** THIS SUBROUTINE USES THE PRECIPITATION TO CALCULATE THE GROUND    4002.   
C**** WATER, GROUND ICE, SNOW COVER, AND RUNOFF                         4002.5  
C****                                                                   4003.   
C**** RUN1 IS NOT ACUMULATED IN ADAILY FOR DIAG6                        4003.5  
C****                                                                   4004.   

#include "BD2G04.COM"                                                   4004.5  
#if ( defined OCEAN_3D || defined ML_2D )
#include "AGRID.h"
#endif

      COMMON U,V,T,P,Q                                                  4005.   
      COMMON/WORK1/CONV(IM0,JM0,LM0),PK(IM0,JM0,LM0),PREC(IM0,JM0),
     &   TPREC(IM0,JM0)                                                 4005.5  
      COMMON/FRMIC/ FRMDICE(JM0)
      DATA SHW/4185./,SHI/2060./,RHOI/916.6/                            4006.   
      DATA Z1I/.1/,Z1E/.1/,Z2LI/2.9/                                    4006.5  
      DATA RHOW/1000./,Z2OIM/0.9/,TFO/-1.56/                            4007.   
      DATA TTRUNC/0./                                                   4007.5  
      DATA IFIRST/1/                                                    4008.   
C****                                                                   4008.5  
C**** FDATA  2  LAND COVERAGE (1)                                       4009.   
C****        3  RATIO OF LAND ICE COVERAGE TO LAND COVERAGE (1)         4009.5  
C****                                                                   4010.   
C**** ODATA  1  OCEAN TEMPERATURE (C)                                   4010.5  
C****        2  RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)       4011.   
C****        3  OCEAN ICE AMOUNT OF SECOND LAYER (KG/M**2)              4011.5  
C****                                                                   4012.   
C**** GDATA  1  OCEAN ICE SNOW AMOUNT (KG/M**2)                         4012.5  
C****        2  EARTH SNOW AMOUNT (KG/M**2)                             4013.   
C****        3  OCEAN ICE TEMPERATURE OF FIRST LAYER (C)                4013.5  
C****        4  EARTH TEMPERATURE OF FIRST LAYER (C)                    4014.   
C****        5  EARTH WATER OF FIRST LAYER (KG/M**2)                    4014.5  
C****        6  EARTH ICE OF FIRST LAYER (KG/M**2)                      4015.   
C****        7  OCEAN ICE TEMPERATURE OF SECOND LAYER (C)               4015.5  
C****       12  LAND ICE SNOW AMOUNT (KG/M**2)                          4016.   
C****       13  LAND ICE TEMPERATURE OF FIRST LAYER (C)                 4016.5  
C****       14  LAND ICE TEMPERATURE OF SECOND LAYER (C)                4017.   
C****                                                                   4017.5  
C**** VDATA  9  WATER FIELD CAPACITY OF FIRST LAYER (KG/M**2)           4018.   
C****                                                                   4018.5  
      IF(IFIRST.NE.1) GO TO 10                                          4019.   
      IFIRST=0                                                          4019.5  
      SHA=RGAS/KAPA                                                     4020.   
      ACE1I=Z1I*RHOI                                                    4020.5  
      AC2OIM=Z2OIM*RHOI                                                 4021.   
      ATRUNC=2.**(-13)                                                  4021.5  
      ACE2LI=Z2LI*RHOI                                                  4022.   
      HC1I=ACE1I*SHI                                                    4022.5  
      HC1DE=Z1E*1129950.                                                4023.   
      DO J=1,JM  
        FRMDICE(j)=0.0
      ENDDO
      print *,'SNOWMAX=4*ACE1I (1m)'
C****                                                                   4023.5  
C**** OUTSIDE LOOP OVER J AND I, EXECUTED ONCE FOR EACH GRID POINT      4024.   
C****                                                                   4024.5  
   10 DO 980 J=1,JM                                                     4025.   
      IMAX=IM                                                           4025.5  
      IF(J.EQ.1.OR.J.EQ.JM) IMAX=1                                      4026.   
         BENRGP=0.                                                      4027.   
         BEDIFS=0.                                                      4028.   
         BERUN0=0.                                                      4029.5  
         BERUN2=0.                                                      4030.   
         BDIFS=0.                                                       4032.   
         BRUN0=0.                                                       4033.5  
	 BRUNS0=0.
         BRUN2=0.                                                       4034.5  
      DO 960 I=1,IMAX                                                   4036.5  
      IF(PREC(I,J).LE.0.) GO TO 960                                     4037.   
C****                                                                   4037.5  
C**** DETERMINE SURFACE CONDITIONS                                      4038.   
C****                                                                   4038.5  
      PLAND=FDATA(I,J,2)                                                4039.   
      PWATER=1.-PLAND                                                   4039.5  
      PLICE=FDATA(I,J,3)*PLAND                                          4040.   
      PEARTH=PLAND-PLICE                                                4040.5  
      ROICE=ODATA(I,J,2)                                                4041.   
      POICE=ROICE*PWATER                                                4041.5  
      POCEAN=PWATER-POICE                                               4042.   
         JR=J
         DXYPJ=DXYP(J)                                                  4043.   
         RUN0S=0.                                                       4043.5  
         DIFSS=0.                                                       4044.   
C**** CALCULATE PRECIPITATION HEAT FLUX (FALLS AT 0 DEGREES CENTIGRADE) 4044.5  
!     PRCP=PREC(I,J)                                                     4045.   
!    07.18.2006
      PRCP=PREC(I,J)*prlnd2total(j,mndriver)  
      TPRCP=TPREC(I,J)                                                  4045.5  
      IF(TPRCP.LT.0.) GO TO 30                                          4046.   
C     EPRCP=PRCP*TPRCP*SHW                                              4046.5  
      EPRCP=0.                                                          4047.   
      ENRGP=EPRCP                                                       4047.5  
      GO TO 50                                                          4048.   
C     EPRCP=PRCP*TPRCP*SHI                                              4048.5  
   30 EPRCP=0.                                                          4049.   
      ENRGP=EPRCP-PRCP*LHM                                              4049.5  
         AIJ(I,J,70)=AIJ(I,J,70)+PRCP                                   4050.   
C****                                                                   4050.5  
   50 CONTINUE
C****                                                                   4114.   
  400 IF(PLICE.LE.0.) GO TO 600                                         4114.5  
C****                                                                   4115.   
C**** LAND ICE                                                          4115.5  
C****                                                                   4116.   
      SNOW=GDATA(I,J,12)                                                4116.5  
      TG1=GDATA(I,J,13)                                                 4117.   
      TG2=GDATA(I,J,14)                                                 4117.5  
         BENRGP=BENRGP+ENRGP*PLICE                                      4118.   
         AIJ(I,J,67)=AIJ(I,J,67)+ENRGP                                  4118.5  
      HC1=HC1I+SNOW*SHI                                                 4119.   
      RUN0=0.                                                           4119.5  
        if(j.eq.-42)then
        print *,' J=',J
        print *,TPRCP,EPRCP,-TG1*HC1
        endif
      IF(TPRCP.LT.0.) GO TO 480                                         4120.   
      IF(EPRCP.LT.-TG1*HC1) GO TO 460                                   4120.5  
C**** RAIN HEATS UP TG1 TO FREEZING POINT AND MELTS SOME SNOW OR ICE    4121.   
      DWATER=(TG1*HC1+EPRCP)/LHM                                        4121.5  
      TG1=0.                                                            4122.   
      RUN0=DWATER+PRCP                                                  4122.5  
c
c	RUNS0 does not include runoff due melting of land ice
	RUNS0=DMIN1(DWATER,SNOW+FRMDICE(J))+PRCP
        if(j.eq.-42)then
        print *,'FRMDICE J=',J
        print *,FRMDICE(J),SNOW,DWATER
        endif
        IF(DWATER.GT.SNOW) THEN
          FRMDICE(j)=FRMDICE(j)-(DWATER-SNOW)
          if( FRMDICE(j).lt.0.0) FRMDICE(j)=0.0
        if(j.eq.-42)then
        print *,'After melting',FRMDICE(J)
        endif
        ENDIF
c
      IF(DWATER.GT.SNOW) GO TO 440                                      4123.   
C**** RAIN MELTS SOME SNOW                                              4123.5  
      SNOW=SNOW-DWATER                                                  4124.   
      GO TO 580                                                         4124.5  
C**** RAIN MELTS ALL SNOW AND SOME ICE, ICE MOVES UP THROUGH THE LAYERS 4125.   
  440 DIFS=SNOW-DWATER                                                  4125.5  
      SNOW=0.                                                           4126.   
      TG1=-TG2*DIFS/ACE1I                                               4126.5  
      EDIFS=DIFS*(TG2*SHI-LHM)                                          4127.   
      ERUN2=EDIFS                                                       4127.5  
      GO TO 560                                                         4128.   
C**** RAIN COOLS TO FREEZING POINT AND HEATS UP TG1                     4128.5  
  460 TG1=TG1+EPRCP/HC1                                                 4129.   
      RUN0=PRCP                                                         4129.5  
        if(j.eq.-42)then
        print *,'After 460 TG1=',TG1,'PRCP=',PRCP
        endif
c     For runoff added on 7/30/03
         RUNS0=RUN0
      GO TO 590                                                         4130.   
C**** SNOW INCREASES SNOW AMOUNT AND SNOW TEMPERATURE RECOMPUTES TG1    4130.5  
  480 TG1=(TG1*HC1+EPRCP)/(HC1+PRCP*SHI)                                4131.   
      SNOW=SNOW+PRCP                                                    4131.5  
        if(j.eq.-42)then
        print *,'After 480 TG1=',TG1,'SNOW=',SNOW
        endif
c     GO TO 580
c     IF(SNOW.LE.ACE1I) GO TO 580                                       4132.   
      IF(SNOW.LE.4.0*ACE1I) GO TO 580 
c     4.*ACE1I=360 kg/m^2 of snow, for used function for show density
c     this corresponds to 1 m deep show.
C**** SNOW IS COMPACTED INTO ICE, ICE MOVES DOWN THROUGH THE LAYERS     4132.5  
c     DIFS=SNOW-.9*ACE1I                                                4133.   
c     SNOW=.9*ACE1I                                                     4133.5  
      DIFS=SNOW-4.*ACE1I  
      SNOW=4.*ACE1I      
      FRMDICE(j)=FRMDICE(j)+DIFS
      FRMDICE(j)=0.0
        if(j.eq.-42)then
        print *,'Before 560 DIFS=',DIFS,'SNOW=',SNOW
        print *,'FRMDICE(j)=',FRMDICE(j)
        endif
      EDIFS=DIFS*(TG1*SHI-LHM)                                          4134.   
      ERUN2=DIFS*(TG2*SHI-LHM)                                          4134.5  
      GDATA(I,J,14)=TG2+(TG1-TG2)*DIFS/ACE2LI                           4135.   
  560    BEDIFS=BEDIFS+EDIFS*PLICE                                      4135.5  
         AIJ(I,J,69)=AIJ(I,J,69)+EDIFS                                  4136.   
         BDIFS=BDIFS+DIFS*PLICE                                         4136.5  
         DIFSS=DIFSS+DIFS*PLICE                                         4137.   
         BERUN2=BERUN2+ERUN2*PLICE                                      4137.5  
         AIJ(I,J,72)=AIJ(I,J,72)+ERUN2                                  4138.   
         BRUN2=BRUN2+DIFS*PLICE                                         4138.5  
  580 GDATA(I,J,12)=SNOW                                                4139.   
  590 GDATA(I,J,13)=TG1                                                 4139.5  
         BRUN0=BRUN0+RUN0*PLICE                                         4140.   
	 BRUNS0=BRUNS0+RUNS0*PLICE
        if(j.eq.-42)then
        print *,'After 590 RUNS0=',RUNS0,'RUN0=',RUN0
        endif
c
         RUN0S=RUN0S+RUN0*PLICE                                         4140.5  
         AIJ(I,J,33)=AIJ(I,J,33)+RUN0                                   4141.   
C****                                                                   4141.5  
  600 IF(PEARTH.LE.0.) GO TO 940                                        4142.   
C****                                                                   4142.5  
C**** EARTH                                                             4143.   
C****                                                                   4143.5  
      SNOW=GDATA(I,J,2)                                                 4144.   
      TG1=GDATA(I,J,4)                                                  4144.5  
      WTR1=GDATA(I,J,5)                                                 4145.   
      ACE1=GDATA(I,J,6)                                                 4145.5  
         BENRGP=BENRGP+ENRGP*PEARTH                                     4146.   
         AIJ(I,J,68)=AIJ(I,J,68)+ENRGP                                  4146.5  
      WFC1=VDATA(I,J,9)                                                 4147.   
      WFC2=VDATA(I,J,10) 
      CHI1=(WTR1+ACE1)/WFC1                                             4147.5  
      HC1=HC1DE+WTR1*SHW+(ACE1+SNOW)*SHI                                4148.   
      RUN0=0.                                                           4148.5  
      ERUN0=0.                                                          4149.   
      IF(TPRCP.LT.0.) GO TO 660                                         4149.5  
      IF(TG1.LE.0.) GO TO 620                                           4150.   
C**** RAIN ON GROUND ABOVE FREEZING POINT, RECOMPUTE TG1                4150.5  
      TG1=(TG1*HC1+EPRCP)/(HC1+PRCP*SHW)                                4151.   
      RUN0=DMAX1(PRCP*.5*CHI1,PRCP+WTR1-WFC1)                           4151.5  
      WTR1=WTR1+(PRCP-RUN0)                                             4152.   
      ERUN0=TG1*RUN0*SHW                                                4152.5  
      GO TO 890                                                         4153.   
  620 IF(EPRCP.LT.-TG1*HC1) GO TO 640                                   4153.5  
C**** RAIN HEATS UP TG1 TO FREEZING POINT                               4154.   
      EPRCP=EPRCP+TG1*HC1                                               4154.5  
      TG1=0.                                                            4155.   
      IF(EPRCP.LT.(ACE1+SNOW)*LHM) GO TO 630                            4155.5  
C**** RAIN MELTS SNOW AND ICE AND HEATS UP TG1 ABOVE FREEZING POINT     4156.   
      RUN0=DMAX1((PRCP+SNOW)*.5*CHI1,PRCP+SNOW+WTR1+ACE1-WFC1)          4156.5  
      WTR1=WTR1+ACE1+SNOW+(PRCP-RUN0)                                   4157.   
      TG1=(EPRCP-(ACE1+SNOW)*LHM)/(HC1DE+(WTR1+RUN0)*SHW)               4157.5  
      ACE1=0.                                                           4158.   
      SNOW=0.                                                           4158.5  
      ERUN0=TG1*RUN0*SHW                                                4159.   
      GO TO 880                                                         4159.5  
C**** RAIN MELTS SOME SNOW AND ICE, TG1 IS AT FREEZING POINT            4160.   
  630 DWATER=EPRCP/LHM                                                  4160.5  
      DSNOW=DMIN1(SNOW,DWATER)                                          4161.   
      RUN0=DMAX1((PRCP+DSNOW)*.5*CHI1,PRCP+DSNOW+WTR1+ACE1-WFC1)        4161.5  
      WTR1=WTR1+DWATER+PRCP-RUN0                                        4162.   
      IF(WTR1.LT.0.) WTR1=0.                                            4162.1  
      SNOW=SNOW-DSNOW                                                   4162.5  
      ACE1=ACE1-DWATER+DSNOW                                            4163.   
      GO TO 880                                                         4163.5  
C**** RAIN COOLS TO FREEZING POINT AND HEATS UP TG1                     4164.   
  640 TG1=TG1+EPRCP/HC1                                                 4164.5  
      RUN0=DMAX1(PRCP*.5*CHI1,PRCP+ACE1-WFC1)                           4165.   
      PRCP=PRCP-RUN0                                                    4165.5  
      IF(PRCP*LHM.LT.-TG1*HC1) GO TO 650                                4166.   
C**** SOME RAIN FREEZES AND TG1 HEATS UP TO FREEZING POINT              4166.5  
      DICE=-TG1*HC1/LHM                                                 4167.   
      TG1=0.                                                            4167.5  
      ACE1=ACE1+DICE                                                    4168.   
      WTR1=PRCP-DICE                                                    4168.5  
      GO TO 890                                                         4169.   
C**** RAIN FREEZES AND HEATS UP TG1, BUT STILL BELOW FREEZING POINT     4169.5  
  650 TG1=(TG1*HC1+PRCP*LHM)/(HC1+PRCP*SHI)                             4170.   
      ACE1=ACE1+PRCP                                                    4170.5  
      GO TO 890                                                         4171.   
  660 IF(TG1.LE.0.) GO TO 690                                           4171.5  
      IF(-EPRCP.LT.TG1*HC1) GO TO 670                                   4172.   
C**** NEW SNOW HEATS UP AND COOLS TG1 TO FREEZING POINT                 4172.5  
      EPRCP=EPRCP+TG1*HC1                                               4173.   
      TG1=0.                                                            4173.5  
      SNOW=PRCP                                                         4174.   
      GO TO 700                                                         4174.5  
C**** NEW SNOW HEATS UP TO FREEZING POINT AND COOLS TG1                 4175.   
  670 TG1=TG1+EPRCP/HC1                                                 4175.5  
      IF(PRCP*LHM.LT.TG1*HC1) GO TO 680                                 4176.   
C**** SOME NEW SNOW MELTS UNTIL TG1 COOLS TO FREEZING POINT             4176.5  
      DWATER=TG1*HC1/LHM                                                4177.   
      TG1=0.                                                            4177.5  
      SNOW=PRCP-DWATER                                                  4178.   
      RUN0=DMAX1(DWATER*.5*CHI1,DWATER+WTR1-WFC1)                       4178.5  
      WTR1=WTR1+(DWATER-RUN0)                                           4179.   
      GO TO 880                                                         4179.5  
C**** ALL NEW SNOW MELTS, RECOMPUTE TG1                                 4180.   
  680 TG1=(TG1*HC1-PRCP*LHM)/(HC1+PRCP*SHW)                             4180.5  
      RUN0=DMAX1(PRCP*.5*CHI1,PRCP+WTR1-WFC1)                           4181.   
      WTR1=WTR1+(PRCP-RUN0)                                             4181.5  
      ERUN0=TG1*RUN0*SHW                                                4182.   
      GO TO 890                                                         4182.5  
  690 SNOW=SNOW+PRCP                                                    4183.   
C     Restriction of SNOW cover
c     if(SNOW.gt.ACE1I)then
c       SNOW=ACE1I
c     endif
C
      IF(WTR1.GT.0.) GO TO 700                                          4183.5  
C**** NEW SNOW INCREASES SNOW AMOUNT AND SNOW TEMP RECOMPUTES TG1       4184.   
      TG1=(TG1*HC1+EPRCP)/(HC1+PRCP*SHI)                                4184.5  
      GO TO 880                                                         4185.   
  700 IF(-EPRCP.LT.WTR1*LHM) GO TO 710                                  4185.5  
C**** GROUND WATER FREEZES, RECOMPUTE TG1                               4186.   
      ACE1=ACE1+WTR1                                                    4186.5  
      HC1=HC1DE+(ACE1+SNOW)*SHI                                         4187.   
      TG1=(EPRCP+WTR1*LHM)/HC1                                          4187.5  
      WTR1=0.                                                           4188.   
      GO TO 880                                                         4188.5  
C**** SOME GROUND WATER FREEZES UNTIL SNOW TEMP HEATS TO FREEZING POINT 4189.   
  710 DICE=-EPRCP/LHM                                                   4189.5  
      WTR1=WTR1-DICE                                                    4190.   
      ACE1=ACE1+DICE                                                    4190.5  
      IF(WTR1+ACE1.GT.WFC1) WTR1=.99999*WTR1                            4190.6  
      IF(WTR1+ACE1.GT.WFC1) ACE1=.99999*ACE1                            4190.7  
  880 GDATA(I,J,2)=SNOW                                                 4191.   
  890 GDATA(I,J,4)=TG1                                                  4191.5  
      GDATA(I,J,5)=WTR1                                                 4192.   
      GDATA(I,J,6)=ACE1                                                 4192.5  
         BERUN0=BERUN0+ERUN0*PEARTH                                     4193.   
         BRUN0=BRUN0+RUN0*PEARTH                                        4193.5  
	 RUNS0=RUN0
	 BRUNS0=BRUNS0+RUNS0*PEARTH
         RUN0S=RUN0S+RUN0*PEARTH                                        4194.   
         AIJ(I,J,32)=AIJ(I,J,32)+RUN0                                   4194.5  
C****                                                                   4195.   
C**** ACCUMULATE DIAGNOSTICS                                            4195.5  
C****                                                                   4196.   
  940    DJ(JR,39)=DJ(JR,39)+ENRGP*DXYPJ                                4196.5  
         DJ(JR,45)=DJ(JR,45)+DIFSS*DXYPJ                                4197.   
         DJ(JR,54)=DJ(JR,54)+RUN0S*DXYPJ                                4197.5  
         AIJ(I,J,5)=AIJ(I,J,5)+PREC(I,J)*prlnd2total(j,mndriver)        4198.   
         AIJ(I,J,23)=AIJ(I,J,23)+ENRGP                                  4198.5  
  960 CONTINUE                                                          4199.   
         BJ(J,39)=BJ(J,39)+BENRGP                                       4200.   
         BJ(J,40)=BJ(J,40)+BERUN0                                       4201.   
         BJ(J,41)=BJ(J,41)+BEDIFS                                       4201.5  
!        BJ(J,43)=BJ(J,43)+BERUN2                                       4203.   
!        BJ(J,45)=BJ(J,45)+BDIFS                                        4204.   
C	 Runoff from first layer of soil including ice melting
 	 BJ(J,47)=BJ(J,47)+BRUN0
C	 Runoff from first layer of soil does not include ice melting
         BJ(J,54)=BJ(J,54)+BRUNS0
         BJ(J,46)=BJ(J,46)+BRUN2                                        4206.5  
#if ( defined OCEAN_3D || defined ML_2D )
C	 Runoff from first layer of soil does not include ice melting
	 if(PLICE+PEARTH.gt.0.0)then
           arunoff(j)=arunoff(j)+BRUNS0/(PLICE+PEARTH)
	 endif
#endif`
  980 CONTINUE                                                          4209.5  
      RETURN                                                            4210.   
      END                                                               4210.5  
