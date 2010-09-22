
#include "ctrparam.h"
!	==========================================================
!
!	MD2G04.F:  Lots of utility functions.
!
!	----------------------------------------------------------
! 
!	Revision History:
!	
!	When	Who		What
!	----	----------	-------	
!	073100	Chien Wang	repack based on CliChem3 & M24x11,
!				  and add cpp.
!
!	==========================================================


      SUBROUTINE DAILY_FOR_OCEAN
C****                                                                   1002.
C**** THIS SUBROUTINE PERFORMS THOSE FUNCTIONS OF THE PROGRAM WHICH     1003.
C**** TAKE PLACE AT THE BEGINNING OF A NEW DAY.                         1004.
C****                                                                   1005.

#include "BD2G04.COM"
#include "AGRID.h"
#include "ODIFF.COM"
#include "BDOCEAN.COM"

      DIMENSION AMONTH(12),JDOFM(13)
      CHARACTER*4 AMONTH                                                1009.1
      LOGICAL PRNT
      DIMENSION XA(1,JM0),XB(1,JM0),OI(IO0,JM0),XOI(IO0,JM0)            1009.5
      DATA AMONTH/'JAN','FEB','MAR','APR','MAY','JUNE','JULY','AUG',    1010.
     *  'SEP','OCT','NOV','DEC'/                                        1011.
      DATA JDOFM/0,31,59,90,120,151,181,212,243,273,304,334,365/        1012.
c     DATA JDPERY/365/,JMPERY/12/,EDPERY/365./,Z1I/.1/,RHOI/916.6/      1013.
c     DATA IFIRST/1/,Z2OIM/0.4/    
      DATA IFIRST/1/
#include "OCNDATA.COM"
C****                                                                   1016.
C**** CALCULATE THE DAILY CALENDAR                                      1035.
C****                                                                   1036.
      if(IFIRST.eq.1)then
      XLHM=334000.
      CW=1000.*4185.
      MICE2MIN=RHOI*Z2OIM
c     print *,'MICE2MIN=',MICE2MIN
      print *,'DAILY_OCN IDAY=',IDAY
      print *,' JYEAR=',JYEAR,' JDAY=',JDAY
      print *,'before 220'
      print *,' JDATE=',JDATE,' JMONTH=',JMONTH
      IFIRST=0
      endif
      DO 210 MONTH=1,JMPERY                                             1039.
      IF(JDAY.LE.JDOFM(MONTH+1)) GO TO 220                              1040.
  210 CONTINUE                                                          1041.
  220 JDATE=JDAY-JDOFM(MONTH)                                           1042.
      JMONTH=AMONTH(MONTH)                                              1043.
      PRNT=TAU.gt.17520.0
      PRNT=.false.
      if(PRNT)then
      write(77, *),' DAILY_OCEAN KOCEAN=',KOCEAN
      write(77, *),' JYEAR=',JYEAR,' JDAY=',JDAY
      write(77, *),' JDATE=',JDATE,' JMONTH=',JMONTH
      write(77, *),' TAU=',TAU,' TOFDAY=',TOFDAY
      endif
c
C****                                                                   1049.
C**** CALCULATE DAILY OCEAN DATA FROM CLIMATOLOGY                       1050.
C****                                                                   1051.
C**** ODATA  1  OCEAN TEMPERATURE (C)                                   1052.
C****        2  RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)       1053.
C****        3  OCEAN ICE AMOUNT OF SECOND LAYER (KG/M**2)              1054.
C****                                                                   1055.
C**** READ IN TWO MONTHS OF OCEAN DATA                                  1056.
      do 385 j=1,JM
      do 385 i=1,IO
       OLDICEM(i,j)=ODATAT(i,j,3)*ODATAT(i,j,2)
  385 continue
      REWIND 515 
      IF(JDAY.GE.16) GO TO 310                                          1057.
      MD=JDATE+15                                                       1058.
      GO TO 320                                                         1059.
  310 IF(JDAY.LE.350) GO TO 340                                         1060.
      MD=JDATE-16                                                       1061.
  320 READ (515) M,XO                                                   1062.
      MDMAX=31                                                          1063.
      DO 330 MX=1,10                                                    1064.
  330 READ (515) M                                                      1065.
      READ (515) M,(((ODATAT(I,J,K),I=1,IO),J=1,JM),K=1,3)              1066.
      GO TO 400                                                         1067.
  340 DO 350 MX=1,12                                                    1068.
      READ (515) M,(((ODATAT(I,J,K),I=1,IO),J=1,JM),K=1,3)              1069.
      IF(M.EQ.MONTH) GO TO 360                                          1070.
      IF(M+1.EQ.MONTH.AND.JDATE.LT.16) GO TO 370                        1071.
  350 CONTINUE                                                          1072.
      STOP 2                                                            1073.
  360 IF(JDATE.EQ.16) GO TO 480                                         1074.
      MDMAX=JDOFM(MONTH+1)-JDOFM(MONTH)                                 1075.
      MD=JDATE-16                                                       1076.
      GO TO 380                                                         1077.
  370 MDMAX=JDOFM(MONTH)-JDOFM(MONTH-1)                                 1078.
      MD=MDMAX+JDATE-16                                                 1079.
  380 READ (515) M,XO                                                   1080.
C**** INTERPOLATE OCEAN DATA TO CURRENT DAY                             1081.
  400 X1=FLOAT(MDMAX-MD)/MDMAX                                          1082.
      X2=1.-X1                                                          1083.
      DO 420 K=1,3                                                      1084.
      DO 420 J=1,JM                                                     1085.
      DO 420 I=1,IO                                                     1086.
  420 ODATAT(I,J,K)=X1*ODATAT(I,J,K)+X2*XO(I,J,K)                       1087.
  480 REWIND 515                                                        1088.
C**** WHEN TGO IS NOT DEFINED, MAKE IT A REASONALBE VALUE               1089.
C**** ZERO OUT SNOWOI, TG1OI, TG2OI AND ACE2OI IF THERE IS NO OCEAN ICE 1100.
      DO 426 J=1,JM                                                     1090.
      DO 426 I=1,IO                                                     1091.
      PLAND=CLAND_OCEAN(I,J)
      if(PLAND.ge.1.0)THEN
        do k=1,3
          ODATAT(I,J,K)=0.0
        enddo  
      GDATA_ICE(I,J,1)=0.  
      GDATA_ICE(I,J,3)=0.   
      GDATA_ICE(I,J,2)=0.    
      ELSE
        IF(ODATAT(I,J,1).LT.-1.56) ODATAT(I,J,1)=-1.56 
        IF(ODATAT(I,J,2).LE.0.) ODATAT(I,J,3)=0.
      ENDIF
  426 CONTINUE                                                          1093.
C**** REDUCE THE RATIO OF OCEAN ICE TO WATER BY .1*RHOI/ACEOI           1094.
c     DO 490 J=1,JM                                                     1095.
c     DO 490 I=1,IO                                                     1096.
c     IF(ODATAT(I,J,2).LE.0.) GO TO 490                                 1097.
c     BYZICE=RHOI/(Z1I*RHOI+ODATAT(I,J,3))                              1097.1
c     ODATAT(I,J,2)=ODATAT(I,J,2)*(1.-.06*(BYZICE-1./5.))               1098.
c 490 CONTINUE                                                          1099.
      IF(KOCEAN.EQ.0.) THEN
        print *,'DAILY KOCEAN=',KOCEAN
        DO 640 J=1,JM 
        DO 640 I=1,IO  
        do nn=1,5
          ODATA_OCEAN(I,J,nn)=ODATAT(I,J,nn)
        enddo
!!!       ODATA_OCEAN(I,J,1)=ODATAT(I,J,1)+4.
c       GDATA_ICE(I,J,3)=-1.56
c       GDATA_ICE(I,J,2)=-1.56
  640   CONTINUE     
c     print *,'ICE'
c     print *,(ODATA_OCEAN(1,j,2),j=1,jm)
c     write(78,*),(ODATA_OCEAN(i,35,2),i=1,im)
c     write(78,*),(GDATA_ICE(i,35,2),i=1,im)
       RETURN                                                           1108.
      ENDIF
C****                                                                   1108.01
C**** CALCULATE DAILY OCEAN MIXED LAYER DEPTHS FROM CLIMATOLOGY         1108.02
C****                                                                   1108.03
C**** SAVE PREVIOUS DAY'S MIXED LAYER DEPTH IN WORK2                    1108.04
  500 DO 510 J=1,JM                                                     1108.05
      DO 510 I=1,IO                                                     1108.06
      ZMLOLD(I,J)=Z1O(I,J)
      Z1OOLD(I,J)=Z1O(I,J)                                              1108.07
  510 CONTINUE
c      go to 626
C**** READ IN TWO MONTHS OF OCEAN DATA                                  1108.08
      IF(JDAY.GE.16) GO TO 520                                          1108.09
      MD=JDATE+15                                                       1108.1
      GO TO 530                                                         1108.11
  520 IF(JDAY.LE.350) GO TO 550                                         1108.12
      MD=JDATE-16                                                       1108.13
  530 READ (515) M,XZO,XOI,XZO,XZO                                      1108.14
      MDMAX=31                                                          1108.15
      DO 540 MX=1,10                                                    1108.16
  540 READ (515) M                                                      1108.17
      READ (515) M,Z1O,OI,Z1O,Z1O                                       1108.18
      GO TO 600                                                         1108.19
  550 DO 560 MX=1,12                                                    1108.2
      READ (515) M,Z1O,OI,Z1O,Z1O                                       1108.21
      IF(M.EQ.MONTH) GO TO 570                                          1108.22
      IF(M+1.EQ.MONTH.AND.JDATE.LT.16) GO TO 580                        1108.23
  560 CONTINUE                                                          1108.24
      STOP 2                                                            1108.25
  570 IF(JDATE.EQ.16) GO TO 625                                         1108.26
      MDMAX=JDOFM(MONTH+1)-JDOFM(MONTH)                                 1108.27
      MD=JDATE-16                                                       1108.28
      GO TO 590                                                         1108.29
  580 MDMAX=JDOFM(MONTH)-JDOFM(MONTH-1)                                 1108.3
      MD=MDMAX+JDATE-16                                                 1108.31
  590 READ (515) M,XZO,XOI,XZO,XZO                                      1108.32
C**** INTERPOLATE OCEAN DATA TO CURRENT DAY                             1108.33
  600 X1=FLOAT(MDMAX-MD)/MDMAX                                          1108.34
      X2=1.-X1                                                          1108.35
      DO 610 J=1,JM                                                     1108.36
      DO 610 I=1,IO                                                     1108.37
      PWATER=1.-CLAND_OCEAN(I,J)
      if(PWATER.gt.0.0)then
      OI(I,J)=X1*OI(I,J)+X2*XOI(I,J)                                    1108.371
      IF(OI(I,J).GT.0.) OI(I,J)=OI(I,J)*                                1108.373
     *  (1.-.1*RHOI/(Z1I*RHOI+ODATA_OCEAN(I,J,3)))                      1108.374
      Z1O(I,J)=X1*Z1O(I,J)+X2*XZO(I,J)                                  1108.38
      Z1OMIN=.09166+.001*(GDATA_ICE(I,J,1)+ODATA_OCEAN(I,J,3)) 
      IF(Z1O(I,J).LT.Z1OMIN) Z1O(I,J)=Z1OMIN                            1108.391
      IF(Z1OMIN.GT.Z12O(I,J)-.1) WRITE(6,605)I,J,MONTH,Z1OMIN,XZO(I,J)  1108.4
  605 FORMAT (' OCEAN ICE CLOSE TO MLD AT I,J,MONTH',3I3,2F10.3)        1108.41
      IF(Z1OMIN.GT.Z12O(I,J)-.1) STOP 8148                              1108.42
C**** PREVENT Z1O, THE MIXED LAYER DEPTH, FROM EXCEEDING Z12O           1108.491
      IF(Z1O(I,J).GT.Z12O(I,J)) Z1O(I,J)=Z12O(I,J)                      1108.494
      endif
  610 CONTINUE                                                          1108.43
  625 REWIND 515                                                        1108.44
       go to 627
  626 CONTINUE
c     print *,' Z1O and Z12O = Z1OAV'
      DO 630 J=1,JM                                                     1108.492
      DO 630 I=1,IO                                                     1108.493
c     Z1O(I,J)=Z1OAV(I,J)
c     Z12O(I,J)=Z1OAV(I,J)
c     IF(Z1O(I,J).GT.Z12O(I,J)-.01) Z1O(I,J)=Z12O(I,J)                  1108.494
      Z1O(I,J)=50.00
      Z12O(I,J)=50.00
  630 CONTINUE                                                          1108.495
      print *,' Z1O and Z12O = 50.0 meters'
  627 CONTINUE
      if(PRNT)then
       write(77,*) ,'Z1O(5,7)=',Z1O(5,7)
      endif
      RETURN                                                            1108.5
C****                                                                   1109.
  901 FORMAT ('0PRESSURE ADDED IN GMP IS',F10.6/)                       1114.
  902 FORMAT ('0MEAN SURFACE PRESSURE OF THE ATMOSPHERE IS',F10.4)      1115.
  910 FORMAT('1',33A4/)                                                 1116.
  915 FORMAT (47X,'DAY',I5,', HR',I3,' (',I2,A5,I5,')',F8.1)            1117.
  920 FORMAT('1')                                                       1118.
      END                                                               1119.
