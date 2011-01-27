
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


      SUBROUTINE DAILY_NEW                                              1001.
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
#if ( defined OCEAN_3D || defined ML_2D)
#include "AGRID.h"
#endif
#if ( defined CLM )
#include "CLM.h"
#endif
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
      DATA OMEGT/282.9/
C****                                                                   1016.
C**** THE GLOBAL MEAN PRESSURE IS KEPT CONSTANT AT PSF MILLIBARS        1017.
C****                                                                   1018.
C**** CALCULATE THE CURRENT GLOBAL MEAN PRESSURE                        1019.
  100 SMASS=0.                                                          1020.
c       print *,' from Daily KOCEAN=',KOCEAN
      nsrps=nsrps+1
      DO 120 J=1,JM                                                     1021.
      SPRESS=0.                                                         1022.
      DO 110 I=1,IM                                                     1023.
  110 SPRESS=SPRESS+P(I,J)                                              1024.
       srps(J)=srps(J)+P(1,J)
      SMASS=SMASS+SPRESS*DXYP(J)                                        1025.
       if(J.EQ.JM/2)PBARSH=SMASS
  120 continue
      PBAR=SMASS/AREAG+PTOP                                             1026.
      PBARNH=2.*(SMASS-PBARSH)/AREAG
      PBARSH=2.*PBARSH/AREAG
      srps(JM+1)=srps(JM+1)+PBARSH
      srps(JM+2)=srps(JM+2)+PBARNH
      srps(JM+3)=srps(JM+3)+PBAR-PTOP
#if ( defined OCEAN_3D)
Cjrs      do j=1,jm
Cjrs       surfpr(j)=surfpr(j)+P(1,J)
Cjrs      enddo
#endif
C**** CORRECT PRESSURE FIELD FOR ANY LOSS OF MASS BY TRUNCATION ERROR   1027.
  130 DELTAP=PSF-PBAR                                                   1028.
      if(DELTAP.gt.1.)then
        print *,' from Daily DELTAP=',DELTAP
        print *,' PBAR=',PBAR,' PBARNH=',PBARNH,' PBARSH=',PBARSH
      endif
c     GO TO 1140
      DO 140 J=1,JM                                                     1029.
      DO 140 I=1,IM                                                     1030.
  140 P(I,J)=P(I,J)+DELTAP                                              1031.
      DOPK=1.                                                           1032.
 1140 continue
C     WRITE (6,901) DELTAP                                              1033.
C****                                                                   1034.
C**** CALCULATE THE DAILY CALENDAR                                      1035.
C****                                                                   1036.
  200 JYEAR=IYEAR+(IDAY-1)/JDPERY                                       1037.
      JDAY=IDAY-(JYEAR-IYEAR)*JDPERY                                    1038.
      DO 210 MONTH=1,JMPERY                                             1039.
      IF(JDAY.LE.JDOFM(MONTH+1)) GO TO 220                              1040.
  210 CONTINUE                                                          1041.
  220 JDATE=JDAY-JDOFM(MONTH)                                           1042.
      JMONTH=AMONTH(MONTH)                                              1043.
C**** CALCULATE SOLAR ANGLES AND ORBIT POSITION                         1044.
       if(ifirst.eq.1.or.HPRNT)then
      print *,' DAILY_ATM IDAY=',IDAY,' IYEAR=',IYEAR
      print *,' JYEAR=',JYEAR,' JDAY=',JDAY
      print *,' JDATE=',JDATE,' JMONTH=',JMONTH
      print *,'OBLIQ=',OBLIQ
         ifirst=0
       endif
      JDSAVE=JDAY                                                       1044.5
      JDATES=JDATE                                                      1044.51
      MONSAV=MONTH                                                      1044.52
c     JDAY=197                                                          1044.53
c     JDATE=16                                                          1044.54
c     MONTH=7                                                           1044.55
!     RSDIST=(1.+ECCN*COS(TWOPI*(JDAY-APHEL)/EDPERY))**2                1045.
!     DEC=COS(TWOPI*(JDAY-SOLS)/EDPERY)*OBLIQ*TWOPI/360.                1046.
!     SIND=SIN(DEC)                                                     1047.
!     COSD=COS(DEC)                                                     1048.
!     03/03/06
!     Fixed calculation of incoming solar radiation
      CALL ORBIT (OBLIQ,ECCN,OMEGT,JDAY-0.5,RSDIST,SIND,COSD,LAMBDA)
      if(JDATE.le.16)then
       do 7231 j=1,JM
       do 7231 L=1,LM
         CLDSST(j,L)=((16-JDATE)*cldssm(j,L,MONTH-1)+
     *               (JDATE+15)*cldssm(j,L,MONTH))/31.
         CLDMCT(j,L)=((16-JDATE)*cldmcm(j,L,MONTH-1)+
     *               (JDATE+15)*cldmcm(j,L,MONTH))/31.
 7231  continue
      else
       do 7241 j=1,JM
       do 7241 L=1,LM
         CLDSST(j,L)=((JDATE-16)*cldssm(j,L,MONTH+1)+
     *               (31-JDATE+16)*cldssm(j,L,MONTH))/31.
         CLDMCT(j,L)=((JDATE-16)*cldmcm(j,L,MONTH+1)+
     *               (31-JDATE+16)*cldmcm(j,L,MONTH))/31.
 7241  continue
      endif
#if (defined OCEAN_3D || defined ML_2D)
       if(JDATE.le.16)then
       do 723 j=1,JM
         TSURFT(j)=((16-JDATE)*TSURFC(j,MONTH-1)+
     *   (JDATE+15)*TSURFC(j,MONTH))/31.
         TLANDT(j)=((16-JDATE)*TLANDC(j,MONTH-1)+
     *   (JDATE+15)*TLANDC(j,MONTH))/31.
 723  continue
      else
       do 724 j=1,JM
        TSURFT(j)=((JDATE-16)*TSURFC(j,MONTH+1)+
     *            (31-JDATE+16)*TSURFC(j,MONTH))/31.
        TLANDT(j)=((JDATE-16)*TLANDC(j,MONTH+1)+
     *            (31-JDATE+16)*TLANDC(j,MONTH))/31.
  724  continue
      endif
!       print *,'From daily_new TSURFD'
!       print *,TSURFD
!       print *,'TSURFT'
!       print *,TSURFT
!       print *,'From daily_new TLANDD'
!       print *,TLANDD
!       print *,'TLANDT'
!       print *,TLANDT
      do 725 j=1,JM
       DT2MGL(j)=TSURFD(j)-TSURFT(j)
       DT2MLD(j)=TLANDD(j)-TLANDT(j)
       TSURFD(j)=0.
       TLANDD(j)=0.
  725  continue
#if ( defined CLM )
        DT2MLAND=0.
        if(PRTREND) then
        AREAL=0.
        do j=1,jm
         DT2MLAND=DT2MLAND+DT2MLD(J)*DXYP(j)*FDATA(1,j,2)
         AREAL=AREAL+DXYP(j)*FDATA(1,j,2)
        end do  !j
        DT2MLAND=DT2MLAND/AREAL
!       print *,'DT2MLD'
!       print *,DT2MLD
       if(JDATE.eq.1)then
        print *,'JDATE=',JDATE,' DT2MLAND=',DT2MLAND
       endif
       endif
!       print *,'AREAL=',AREAL
#endif

#endif

      RETURN                                                            1108.5
C****                                                                   1109.
      ENTRY DAILY_NEW0                                                  1110.
c     IF(TAU.GT.TAUI+DT/7200.) GO TO 200                                1111.
c     GO TO 100                                                         1112.
      go to 200
C*****                                                                  1113.
  901 FORMAT ('0PRESSURE ADDED IN GMP IS',F10.6/)                       1114.
  902 FORMAT ('0MEAN SURFACE PRESSURE OF THE ATMOSPHERE IS',F10.4)      1115.
  910 FORMAT('1',33A4/)                                                 1116.
  915 FORMAT (47X,'DAY',I5,', HR',I3,' (',I2,A5,I5,')',F8.1)            1117.
  920 FORMAT('1')                                                       1118.
      END                                                               1119.
      SUBROUTINE ORBIT (OBLIQ,ECCN,OMEGT,DAY,SDIST,SIND,COSD,LAMBDA)    8201.
C****                                                                   8202.
C**** ORBIT receives the orbital parameters and time of year, and       8203.
C**** returns the distance from the sun and its declination angle.      8204.
C**** The reference for the following caculations is: V.M.Blanco        8205.
C**** and S.W.McCuskey, 1961, "Basic Physics of the Solar System",      8206.
C**** pages 135 - 151.                                                  8207.
C****                                                                   8208.
C**** Program authors: Gary L. Russell and Robert J. Suozzo, 12/13/85   8209.
C****                                                                   8210.
C****        All computations are in double-precision;                  8211.
C****        but the arguments are single-precision.                    8212.
C**** Input: OBLIQ = latitude of tropics in degrees                     8213.
C****        ECCEN = eccentricity of the orbital ellipse                8214.
C****        OMEGT = angle from vernal equinox to perihelion in degrees 8215.
C****        DAY   = day of the year in days; 0 = Jan 1, hour 0         8216.
C****                                                                   8217.
C**** Constants: EDAYPY = Earth days per year = 365                     8218.
C****            VERQNX = occurence of vernal equinox = day 79 = Mar 21 8219.
C****                                                                   8220.
C**** Intermediate quantities:                                          8221.
C****    PERIHE = perihelion during the year in temporal radians        8222.
C****    MA     = mean anomaly in temporal radians = 2J DAY/365 - PERIHE8223.
C****    EA     = eccentric anomaly in radians                          8224.
C****    TA     = true anomaly in radians                               8225.
C****    BSEMI  = semi minor axis in units of the semi major axis       8226.
C****    GREENW = longitude of Greenwich in the Earth's reference frame 8227.
C****                                                                   8228.
C**** Output: DIST = distance to the sun in units of the semi major axis8229.
C****        SDIST = square of DIST                                     8229.5
C****         SIND = sine of the declination angle                      8230.
C****         COSD = cosine of the declination angle                    8231.
C****       LAMBDA = sun longitude in Earth's rotating reference frame  8232.
C****                                                                   8233.
      IMPLICIT REAL*8 (A-H,O-Z)                                         8234.
      REAL*8 MA                                                         8235.
C     REAL*4 SIND,COSD,SDIST,LAMBDA,OBLIQ,ECCN,OMEGT,DAY                8236.
C****                                                                   8237.
      PI = 3.14159265358979D0                                           8238.
      EDAYPY = 365.                                                     8239.
      VERQNX = 79.                                                      8240.
      OMEGA=OMEGT*(PI/180.D0)                                           8241.
      DOBLIQ=OBLIQ*(PI/180.D0)                                          8242.
      ECCEN=ECCN                                                        8243.
C****                                                                   8244.
C**** Determine time of perihelion using Kepler's equation:             8245.
C**** PERIHE-VERQNX = OMEGA - ECCEN sin(OMEGA)                          8246.
C****                                                                   8247.
      PERIHE = OMEGA-ECCEN*SIN(OMEGA)+VERQNX*2.*PI/365.                 8248.
C     PERIHE = DMOD(PERIHE,2.*PI)                                       8249.
      MA = 2.*PI*DAY/365.-PERIHE                                        8250.
      MA = DMOD(MA,2.*PI)                                               8251.
C****                                                                   8252.
C**** Numerically solve Kepler's equation: MA = EA - ECCEN sin(EA)      8253.
C****                                                                   8254.
      EA = MA+ECCEN*(SIN(MA)+ECCEN*SIN(2.*MA)/2.)                       8255.
  110 DEA = (MA-EA+ECCEN*SIN(MA))/(1.-ECCEN*COS(EA))                    8256.
      EA = EA+DEA                                                       8257.
      IF (DABS(DEA).GT.1.D-8)  GO TO 110                                8258.
C****                                                                   8259.
C**** Calculate the distance to the sun and the true anomaly            8260.
C****                                                                   8261.
      BSEMI = DSQRT(1.-ECCEN*ECCEN)                                     8262.
      COSEA = COS(EA)                                                   8263.
      SINEA = SIN(EA)                                                   8264.
      SDIST  = (1.-ECCEN*COSEA)*(1.-ECCEN*COSEA)                        8265.
      TA = DATAN2(SINEA*BSEMI,COSEA-ECCEN)                              8266.
C****                                                                   8267.
C**** Change the reference frame to be the Earth's equatorial plane     8268.
C**** with the Earth at the center and the positive x axis parallel to  8269.
C**** the ray from the sun to the Earth were it at vernal equinox.      8270.
C**** The distance from the current Earth to that ray (or x axis) is:   8271.
C**** DIST sin(TA+OMEGA).  The sun is located at:                       8272.
C****                                                                   8273.
C**** SUN    = (-DIST cos(TA+OMEGA),                                    8274.
C****           -DIST sin(TA+OMEGA) cos(OBLIQ),                         8275.
C****            DIST sin(TA+OMEGA) sin(OBLIQ))                         8276.
C**** SIND   = sin(TA+OMEGA) sin(OBLIQ)                                 8277.
C**** COSD   = sqrt(1-SIND**2)                                          8278.
C**** LAMBDA = atan[tan(TA+OMEGA) cos(OBLIQ)] - GREENW                  8279.
C**** GREENW = 2*3.14159 DAY (EDAYPY-1)/EDAYPY                          8280.
C****                                                                   8281.
      SINDD = SIN(TA+OMEGA)*SIN(DOBLIQ)                                 8282.
      COSD = DSQRT(1.-SINDD*SINDD)                                      8283.
      SIND = SINDD                                                      8284.
C     GREENW = 2.*PI*(DAY-VERQNX)*(EDAYPY+1.)/EDAYPY                    8285.
C     SUNX = -COS(TA+OMEGA)                                             8286.
C     SUNY = -SIN(TA+OMEGA)*COS(DOBLIQ)                                 8287.
C     LAMBDA = DATAN2(SUNY,SUNX)-GREENW                                 8288.
C     LAMBDA = DMOD(LAMBDA,2.*PI)                                       8289.
C****                                                                   8290.
      RETURN                                                            8291.
      END                                                               8292.
