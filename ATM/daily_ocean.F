
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


      SUBROUTINE DAILY_OCEAN                                            1001.
C****                                                                   1002.
C**** THIS SUBROUTINE PERFORMS THOSE FUNCTIONS OF THE PROGRAM WHICH     1003.
C**** TAKE PLACE AT THE BEGINNING OF A NEW DAY.                         1004.
C****                                                                   1005.

#include "BD2G04.COM"                                                   1006.

      COMMON/SPEC2/KM,KINC,COEK,C3LAND(IO0,JM0),C3OICE(IO0,JM0)         1006.1
     *  ,C3LICE(IO0,JM0),WMGE(IO0,JM0),TSSFC(IM0,JM0,4)                 1006.2
      COMMON U,V,T,P,Q                                                  1007.
      COMMON/WORK2/Z1OOLD(IO0,JM0),XO(IO0,JM0,3),XZO(IO0,JM0)           1008.
      COMMON/OLDZO/ZMLOLD(IO0,JM0)
      DIMENSION AMONTH(12),JDOFM(13)                                    1009.
      CHARACTER*4 AMONTH                                                1009.1
      DIMENSION XA(1,JM0),XB(1,JM0),OI(IO0,JM0),XOI(IO0,JM0)            1009.5
      dimension sst1(JM0,3),sst2(JM0,3),dsst(JM0,3),intem(3),
     &  sstmin(12,2)
     &       ,miceo(JM0)
      common/qfl/QFLUX(JM0,0:13),ZOAV(JM0),QFLUXT(JM0)
!     common/TSUR/TSURFC(JM0,0:13),TSURFT(JM0),TSURFD(JM0),DTSURF(JM0)
#include "TSRF.COM"
      common/fixcld/cldssm(JM0,LM0,0:13),cldmcm(JM0,LM0,0:13),
     &  CLDSST(JM0,LM0),
     &  CLDMCT(JM0,LM0)
      common/surps/srps(JM0+3),nsrps
      LOGICAL HPRNT
      common/conprn/HPRNT,JPR,LPR
      data ifirst /1/
      data intem /1,4,5/
      data sstmin /-1.56,-1.56,-0.75,6*0.0,2*-0.75,-1.56,
     * 3*0.0,2*-0.75,3*-1.56,-0.75,3*0.0/
      DATA AMONTH/'JAN','FEB','MAR','APR','MAY','JUNE','JULY','AUG',    1010.
     *  'SEP','OCT','NOV','DEC'/                                        1011.
      DATA JDOFM/0,31,59,90,120,151,181,212,243,273,304,334,365/        1012.
      DATA JDPERY/365/,JMPERY/12/,EDPERY/365./,Z1I/.1/,RHOI/916.6/      1013.
C**** ORBITAL PARAMETERS FOR EARTH FOR YEAR 2000 A.D.                   1014.
      DATA SOLS/173./,APHEL/186./,OBLIQ/23.44/,ECCN/.0167/              1015.
c     DATA SOLS/173./,APHEL/186./,OBLIQ/25.00/,ECCN/.0167/              1015.
C****                                                                   1016.
C**** CALCULATE THE DAILY CALENDAR                                      1035.
C****                                                                   1036.
  200 JYEAR=IYEAR+(IDAY-1)/JDPERY                                       1037.
      JDAY=IDAY-(JYEAR-IYEAR)*JDPERY                                    1038.
      DO 210 MONTH=1,JMPERY                                             1039.
      IF(JDAY.LE.JDOFM(MONTH+1)) GO TO 220                              1040.
  210 CONTINUE                                                          1041.
  220 JDATE=JDAY-JDOFM(MONTH)                                           1042.
      JMONTH=AMONTH(MONTH)                                              1043.
      if(ifirst.eq.1.or.HPRNT)then
      print *,' DAILY_OCEAN IDAY=',IDAY,' IYEAR=',IYEAR
      print *,' JYEAR=',JYEAR,' JDAY=',JDAY
      print *,' JDATE=',JDATE,' JMONTH=',JMONTH
      if(KOCEAN.eq.1)ifirst=0
      endif
c
      IF(KOCEAN.EQ.1) GO TO 500                                         1048.1
C****                                                                   1049.
C**** CALCULATE DAILY OCEAN DATA FROM CLIMATOLOGY                       1050.
C****                                                                   1051.
C**** ODATA  1  OCEAN TEMPERATURE (C)                                   1052.
C****        2  RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)       1053.
C****        3  OCEAN ICE AMOUNT OF SECOND LAYER (KG/M**2)              1054.
C****                                                                   1055.
C**** READ IN TWO MONTHS OF OCEAN DATA                                  1056.
      do 385 j=1,JM
       miceo(j)=ODATA(1,j,3)*ODATA(1,j,2)
  385 continue
      IF(JDAY.GE.16) GO TO 310                                          1057.
      MD=JDATE+15                                                       1058.
      GO TO 320                                                         1059.
  310 IF(JDAY.LE.350) GO TO 340                                         1060.
      MD=JDATE-16                                                       1061.
  320 READ (515) M,XO                                                   1062.
      MDMAX=31                                                          1063.
      DO 330 MX=1,10                                                    1064.
  330 READ (515) M                                                      1065.
      READ (515) M,(((ODATA(I,J,K),I=1,IO),J=1,JM),K=1,3)               1066.
      GO TO 400                                                         1067.
  340 DO 350 MX=1,12                                                    1068.
      READ (515) M,(((ODATA(I,J,K),I=1,IO),J=1,JM),K=1,3)               1069.
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
  420 ODATA(I,J,K)=X1*ODATA(I,J,K)+X2*XO(I,J,K)                         1087.
  480 REWIND 515                                                        1088.
      DO 255 J=1,JM                                                     1088.5
      SUM1=0.                                                           1088.51
      SUM2=0.                                                           1088.511
      SUM3=0.                                                           1088.512
      CONT1=0.                                                          1088.52
      DO 256 I=1,IO                                                     1088.53
      PLAND=C3LAND(I,J)                                                 1088.54
      POICE= ODATA(I,J,2)*(1.-PLAND)                                    1088.55
      C3OICE(I,J)=POICE                                                 1088.56
      PWATER=1.-PLAND                                                   1088.57
      IF(PWATER.LE.0.) GO TO 256                                        1088.58
      CONT1=CONT1+PWATER                                                1088.59
      SUM1=SUM1+PWATER*ODATA(I,J,1)                                     1088.6
      SUM2=SUM2+PWATER*ODATA(I,J,4)                                     1088.601
      SUM3=SUM3+PWATER*ODATA(I,J,5)                                     1088.602
 256  CONTINUE                                                          1088.61
      IF(CONT1.EQ.0.) GO TO 255                                         1088.62
      IF (J.EQ.1.OR.J.EQ.JM) GO TO 255                                  1088.63
      SUM1=SUM1/CONT1                                                   1088.64
      SUM2=SUM2/CONT1                                                   1088.642
      SUM3=SUM3/CONT1                                                   1088.643
      DO 258 I=1,IO                                                     1088.65
      ODATA(I,J,4)=SUM2                                                 1088.651
      ODATA(I,J,5)=SUM3                                                 1088.652
 258  ODATA(I,J,1)=SUM1                                                 1088.66
 255  CONTINUE                                                          1088.67
      DO 257 J=2,JMM1                                                   1088.68
      SUM1=0.                                                           1088.69
      SUM2=0.                                                           1088.7
      CONT1=0.                                                          1088.71
      DO 254 I=1,IO                                                     1088.72
      POICE=ODATA(I,J,2)*(1.-C3LAND(I,J))                               1088.73
      SUM1=SUM1+POICE                                                   1088.74
      SUM2=SUM2+POICE*ODATA(I,J,3)                                      1088.75
 254  CONT1=CONT1+(1.-C3LAND(I,J))                                      1088.76
      IF(SUM1.LE.0.) GO TO 425                                          1088.77
      SUM2=SUM2/SUM1                                                    1088.78
      DO 423 I=1,IO                                                     1088.79
 423  ODATA(I,J,3)=SUM2                                                 1088.8
 425  CONTINUE                                                          1088.81
      IF(CONT1.LE.0.) GO TO 257                                         1088.82
      RATIO=SUM1/CONT1                                                  1088.83
      DO 253 I=1,IO                                                     1088.84
 253  ODATA(I,J,2)=RATIO                                                1088.85
 257  CONTINUE                                                          1088.86
      go to 678
      DO 251 J=1,2                                                      1088.87
      DO 251 I=1,IO                                                     1088.88
      ODATA(I,J,1)=ODATA(I,3,1)                                         1088.881
  251 ODATA(I,J,2)=1.                                                   1088.89
      DO 428 J=1,2                                                      1088.9
      DO 428 I=1,IO                                                     1088.91
      ODATA(I,J,4)=ODATA(I,3,4)                                         1088.911
      ODATA(I,J,5)=ODATA(I,3,5)                                         1088.912
  428 ODATA(I,J,3)=ODATA(I,3,3)                                         1088.92
  678 continue
c     print *,'ICE FractionS'
c     print *,(ODATA(1,J,2),j=1,jm)
C     Skip adjustment
c     go to 950
      if(JDATE.eq.46)then
      print *,' before'
      do 567 M=1,5
        print *,'  '
        print *,' ODATA ',M
        print *,(ODATA(1,J,M),J=1,JM)
  567 continue
      endif
      do 558 J=1,JM
        do 559 ntem=1,3
           ITEM=intem(ntem)
           sst1(J,ntem)=ODATA(1,J,ITEM)
           sst2(J,ntem)=ODATA(1,J,ITEM)
           dsst(J,ntem)=0.
  559   continue
        if(ODATA(1,J,2).ge.0.2)then
        dmice=ODATA(1,J,3)*ODATA(1,j,2)-miceo(J)
        if(dmice.ge.0.0)then
          do 561 ntem=1,3
           if(sst1(J,ntem).gt.-1.56)then
             sst2(J,ntem)=-1.56
             dsst(J,ntem)=-1.56-sst1(J,ntem)
           endif
  561     continue
        else
          do 569 ntem=1,3
           if(sst1(J,ntem).gt.0.0)then
             sst2(J,ntem)=0.0
             dsst(J,ntem)=-sst1(J,ntem)
           endif
  569     continue
        endif
        else
          ODATA(1,J,2)=0.
          ODATA(1,J,3)=0.
        endif
  558 continue
      do 562 j=2,JM/2
        jnr=JM-j+1
        do 563 ntem=1,3
         if(dsst(j+1,ntem).eq.0.0.and.dsst(j,ntem).eq.0.0
     *          .and.dsst(j-1,ntem).ne.0.0) then
           sst2(j,ntem)=sst1(j,ntem)+0.5*dsst(j-1,ntem)
           sst2(j+1,ntem)=sst1(j+1,ntem)+0.25*dsst(j-1,ntem)
         endif
         if(dsst(jnr-1,ntem).eq.0.0.and.dsst(jnr,ntem).eq.0.0
     *           .and.dsst(jnr+1,ntem).ne.0.0) then
           sst2(jnr,ntem)=sst1(jnr,ntem)+0.5*dsst(jnr+1,ntem)
           sst2(jnr-1,ntem)=sst1(jnr-1,ntem)+0.25*dsst(jnr+1,ntem)
         endif
  563   continue
  562 continue
      do 663 J=1,JM
        do 664 ntem=1,3
           ITEM=intem(ntem)
           ODATA(1,J,ITEM)=sst2(J,ntem)
  664   continue
  663 continue
      if(JDATE.eq.46)then
      print *,' after'
      do 557 M=1,5
        print *,'  '
        print *,' ODATA ',M
        print *,(ODATA(1,J,M),J=1,JM)
  557 continue
      endif
C
      go to 955
  950 continue
      if(ifirst.eq.1)then
      print *,' Adjustment of SST and sea ice is skiped'
      print *,' Adjustment of SST and sea ice is skiped'
      print *,' Adjustment of SST and sea ice is skiped'
      ifirst=0
      endif
  955 continue
      if(ifirst.eq.1)then
      print *,' With adjustment of SST and sea ice '
      print *,' With adjustment of SST and sea ice '
      print *,' With adjustment of SST and sea ice '
        ifirst=0
      endif
c     JDAY=JDSAVE                                                       1088.93
c     JDATE=JDATES                                                      1088.94
c     MONTH=MONSAV                                                      1088.95
C**** WHEN TGO IS NOT DEFINED, MAKE IT A REASONALBE VALUE               1089.
      DO 426 J=1,JM                                                     1090.
      DO 426 I=1,IO                                                     1091.
      IF(ODATA(I,J,1).LT.-10.) ODATA(I,J,1)=-10.                        1092.
  426 CONTINUE                                                          1093.
C**** REDUCE THE RATIO OF OCEAN ICE TO WATER BY .1*RHOI/ACEOI           1094.
      DO 490 J=1,JM                                                     1095.
      DO 490 I=1,IO                                                     1096.
      IF(ODATA(I,J,2).LE.0.) GO TO 490                                  1097.
      BYZICE=RHOI/(Z1I*RHOI+ODATA(I,J,3))                               1097.1
      ODATA(I,J,2)=ODATA(I,J,2)*(1.-.06*(BYZICE-1./5.))                 1098.
  490 CONTINUE                                                          1099.
C**** ZERO OUT SNOWOI, TG1OI, TG2OI AND ACE2OI IF THERE IS NO OCEAN ICE 1100.
      DO 620 J=1,JM                                                     1101.
      DO 620 I=1,IO                                                     1102.
      IF(ODATA(I,J,2).GT.0.) GO TO 620                                  1103.
      GDATA(I,J,1)=0.                                                   1104.
      GDATA(I,J,3)=0.                                                   1105.
      GDATA(I,J,7)=0.                                                   1106.
  620 CONTINUE                                                          1107.
      RETURN                                                            1108.
C****                                                                   1108.01
C**** CALCULATE DAILY OCEAN MIXED LAYER DEPTHS FROM CLIMATOLOGY         1108.02
C****                                                                   1108.03
C**** SAVE PREVIOUS DAY'S MIXED LAYER DEPTH IN WORK2                    1108.04
  500 DO 510 J=1,JM                                                     1108.05
      DO 510 I=1,IO                                                     1108.06
      ZMLOLD(I,J)=Z1O(I,J)
  510 Z1OOLD(I,J)=Z1O(I,J)                                              1108.07
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
      OI(I,J)=X1*OI(I,J)+X2*XOI(I,J)                                    1108.371
      IF(OI(I,J).GT.0.) OI(I,J)=OI(I,J)*                                1108.373
     *  (1.-.1*RHOI/(Z1I*RHOI+ODATA(I,J,3)))                            1108.374
      Z1O(I,J)=X1*Z1O(I,J)+X2*XZO(I,J)                                  1108.38
      Z1OMIN=.09166+.001*(GDATA(I,J,1)+ODATA(I,J,3))                    1108.39
      IF(Z1O(I,J).LT.Z1OMIN) Z1O(I,J)=Z1OMIN                            1108.391
      IF(Z1OMIN.GT.Z12O(I,J)-.1) WRITE(6,605)I,J,MONTH,Z1OMIN,XZO(I,J)  1108.4
  605 FORMAT (' OCEAN ICE CLOSE TO MLD AT I,J,MONTH',3I3,2F10.3)        1108.41
      IF(Z1OMIN.GT.Z12O(I,J)-.1) STOP 8148                              1108.42
  610 CONTINUE                                                          1108.43
  625 REWIND 515                                                        1108.44
      DO 628 J=1,JM                                                     1108.441
      SUM1=0.                                                           1108.442
      CONT1=0.                                                          1108.444
      DO 626 I=1,IO                                                     1108.445
      C3OICE(I,J)=OI(I,J)*(1.-C3LAND(I,J))                              1108.446
      PWATER=1.-C3LAND(I,J)                                             1108.447
      IF(PWATER.LE.0.) GO TO 626                                        1108.448
      CONT1=CONT1+PWATER                                                1108.449
      SUM1=SUM1+Z1O(I,J)*PWATER                                         1108.45
  626 CONTINUE                                                          1108.452
      IF(CONT1.LE.0.) GO TO 628                                         1108.453
      IF(J.EQ.1.OR.J.EQ.JM) GO TO 628                                   1108.454
      SUM1=SUM1/CONT1                                                   1108.455
      DO 627 I=1,IO                                                     1108.457
      Z1O(I,J)=SUM1                                                     1108.458
  627 CONTINUE                                                          1108.459
  628 CONTINUE                                                          1108.46
      DO 629 J=1,2                                                      1108.461
      DO 629 I=1,IO                                                     1108.462
      Z1O(I,J)=Z1O(I,3)                                                 1108.463
  629 CONTINUE                                                          1108.464
C**** PREVENT Z1O, THE MIXED LAYER DEPTH, FROM EXCEEDING Z12O           1108.491
      DO 630 J=1,JM                                                     1108.492
      DO 630 I=1,IO                                                     1108.493
CCC   Z1O(I,J)=ZOAV(J)
      IF(Z1O(I,J).GT.Z12O(I,J)-.01) Z1O(I,J)=Z12O(I,J)                  1108.494
  630 CONTINUE                                                          1108.495
c     print *,' DAILY JDATE=',JDATE,' MONTH=',MONTH
c     print *,'TSURFD'
c     print *,TSURFD
c     print *,'TSURFT'
c     print *,TSURFT
      do 725 j=1,JM
!       DTSURF(j)=TSURFD(j)-TSURFT(j)
        DT2MGL(j)=TSURFD(j)-TSURFT(j)
        TSURFD(j)=0.
  725 continue
      if(JDATE.le.16)then
       do 723 j=1,JM
        QFLUXT(j)=((16-JDATE)*QFLUX(j,MONTH-1)+
     *            (JDATE+15)*QFLUX(j,MONTH))/31.
        TSURFT(j)=((16-JDATE)*TSURFC(j,MONTH-1)+
     *            (JDATE+15)*TSURFC(j,MONTH))/31.
  723  continue
      else
       do 724 j=1,JM
        QFLUXT(j)=((JDATE-16)*QFLUX(j,MONTH+1)+
     *            (31-JDATE+16)*QFLUX(j,MONTH))/31.
        TSURFT(j)=((JDATE-16)*TSURFC(j,MONTH+1)+
     *            (31-JDATE+16)*TSURFC(j,MONTH))/31.
  724  continue
      endif
c     print *,' NEW TSURFT'
c     print *,TSURFT
      RETURN                                                            1108.5
C****                                                                   1109.
  901 FORMAT ('0PRESSURE ADDED IN GMP IS',F10.6/)                       1114.
  902 FORMAT ('0MEAN SURFACE PRESSURE OF THE ATMOSPHERE IS',F10.4)      1115.
  910 FORMAT('1',33A4/)                                                 1116.
  915 FORMAT (47X,'DAY',I5,', HR',I3,' (',I2,A5,I5,')',F8.1)            1117.
  920 FORMAT('1')                                                       1118.
      END                                                               1119.
