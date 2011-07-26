
#include "ctrparam.h"

!	==========================================================
!
!	TRVDATA.F: THIS SUBROUTINE CALCULATS TRANSIENT VDATA
!	
!	==========================================================

      SUBROUTINE TRVDATA
C****                                                                   1502.   
C****                                                                   1505.   


#include "ODIFF.COM"
#include "BD2G04.COM"
#include "RADCOM.COM"
#include "run.COM" 


      COMMON/SPEC2/KM,KINC,COEK,C3LAND(IO0,JM0),C3OICE(IO0,JM0)         1506.1  
     *  ,C3LICE(IO0,JM0),WMGE(IO0,JM0),TSSFC(IM0,JM0,4)                 1506.2  
      common/veg/TRVEG,IYVEG
      logical TRVEG
      dimension VMASK(JM0)
      JM=JM0                                                            1531.   
      IM=IM0
      IO=IO0                                                            1532.5  
      LM=LM0                                                            1533.   
C**** READ IN EARTH RATIOS FOR THE 8 VEGETATION TYPES AND THE VADATA :  1816.   
C     VADATA(TYPE,SEASON,1)=GROUND ALBEDO FOR A GIVEN TYPE AND SEASON   1817.   
C         1       2       3       4       5       6       7       8     1818.   
C       DESRT   TNDRA   GRASS   SHRUB   TREES   DECID   EVERG   RAINF   1819.   
C SPRN   0.35,   0.12,   0.16,   0.16,   0.14,   0.18,   0.12,   0.11,  1820.   
C SUMR   0.35,   0.12,   0.20,   0.18,   0.14,   0.12,   0.12,   0.11,  1821.   
C FALL   0.35,   0.17,   0.20,   0.25,   0.17,   0.15,   0.15,   0.11,  1822.   
C WNTR   0.35,   0.15,   0.18,   0.20,   0.12,   0.12,   0.11,   0.11/  1823.   
C                                                                       1824.   
C     VADATA(TYPE,SEASON,2)=RATIO OF NEAR IR ALBEDO TO VIS ALBEDO FOR...1825.   
C         1       2       3       4       5       6       7       8     1826.   
C       DESRT   TNDRA   GRASS   SHRUB   TREES   DECID   EVERG   RAINF   1827.   
C SPRN    1.0,    3.0,    3.0,    3.0,    3.0,    3.0,    3.0,    3.0,  1828.   
C SUMR    1.0,    3.3,    3.5,    3.0,    3.3,    4.0,    3.0,    3.0,  1829.   
C FALL    1.0,    3.5,    4.0,    3.0,    3.5,    5.0,    3.0,    3.0,  1830.   
C WNTR    1.0,    3.2,    3.5,    3.0,    3.2,    4.0,    3.0,    3.0/  1831.   
C                                                                       1832.   
C     VADATA(TYPE,1,3)=MASKING DEPTH FOR A GIVEN TYPE                   1833.   
C                                                                       1834.   
C         1       2       3       4       5       6       7       8     1835.   
C       DESRT   TNDRA   GRASS   SHRUB   TREES   DECID   EVERG   RAINF   1836.   
C        10.,    20.,    20.,    50.,   200.,   500.,   1000.,  2500.,  1837.   
C                                                                       1838.   
C     VADATA(TYPE,1+K,3)=WATER FIELD CAPACITY FOR K-TH GROUND LAYER     1839.   
C                                                                       1840.   
C    1   10.,    30.,    30.,    30.,    30.,    30.,     30.,   200.,  1841.   
C    2   10.,   200.,   200.,   300.,   300.,   450.,    450.,   450.,  1842.   
C   (3)   0.,     0.,     0.,     0.,     0.,     0.,      0.,    0./   1843.   
C                                                                       1844.   
      READ (523) IYVEG
  764 READ (523) (((VDATA(I,J,K),I=1,IO),J=1,JM),K=1,8)
        print *,' From trvdata iyveg=',iyveg
c       print *,'VADATA'
c       do k=1,3
c       print *,' K=',k
c       print '(8f7.2)',((VADATA(I,J,K),I=1,8),J=1,4)
c       enddo
C**** MODIFY THE VADATA IF DESIRED                                      1847.   
C     NO MODIFICATIONS                                                  1848.   
C**** COMPUTE WATER FIELD CAPACITIES FOR GROUND LAYERS 1 AND 2          1849.   
      IOFF=0                                                            1849.1  
      IF(VADATA(4,2,3).LT.100.) IOFF=1                                  1849.2  
      ERROR=.001                                                        1849.3  
      DEFLT=24.                                                         1850.   
      DO 785 L=1,2                                                      1851.   
      DO 780 J=1,JM                                                     1852.   
      DO 780 I=1,IO                                                     1853.   
      WFCIJL=0.                                                         1854.   
      DO 770 K=1,8                                                      1855.   
  770 WFCIJL=WFCIJL+VDATA(I,J,K)*VADATA(K,L+IOFF,3)                     1856.   
      IF (WFCIJL.LT.1.) WFCIJL=DEFLT                                    1857.   
      IF(ISTART.NE.2) GO TO 780
      IF(GDATA(I,J,4*L+1)+GDATA(I,J,4*L+2).LE.WFCIJL) GO TO 780         1858.   
      X=WFCIJL/(GDATA(I,J,4*L+1)+GDATA(I,J,4*L+2)+1.E-3)                1859.   
      GDATA(I,J,4*L+1)=GDATA(I,J,4*L+1)*X                               1860.   
      GDATA(I,J,4*L+2)=GDATA(I,J,4*L+2)*X                               1861.   
  780 VDATA(I,J,L+8)=WFCIJL                                             1862.   
      DEFLT=60.                                                         1863.   
  785 CONTINUE                                                          1864.   
      DO 765 K=1,10                                                     1864.5  
      DO 765 J=2,JMM1                                                   1864.51 
      CONT1=0.                                                          1864.52 
      SUM1=0.                                                           1864.53 
      DO 766 I=1,IO                                                     1864.54 
      PEARTH=C3LAND(I,J)-C3LICE(I,J)                                    1864.55 
      CONT1=CONT1+PEARTH                                                1864.56 
 766  SUM1=SUM1+PEARTH*VDATA(I,J,K)                                     1864.57 
      IF (CONT1.LE.0.) GO TO 765                                        1864.58 
      SUM1=SUM1/CONT1                                                   1864.59 
      DO 767 I=1,IO                                                     1864.6  
 767  VDATA(I,J,K)=SUM1                                                 1864.61 
 765  CONTINUE                                                          1864.62 
c     print *,' BEAR LAND'
c     print '(12f7.2,/,11f7.2)',(VDATA(1,j,1),j=1,JM)
c     print *,' TRVDATA'
c     print *,' WMAX1'
c     print '(12f7.2,/,11f7.2)',(VDATA(1,j,9),j=1,JM)
c     print *,' WMAX2'
c     print '(12f7.2,/,11f7.2)',(VDATA(1,j,10),j=1,JM)
C *************
      DO  K=1,8
c     VADATA(K,4,3)=0.1*VADATA(K,4,3)
      VADATA(K,4,3)=VADATA(K,3,3)
      ENDDO
      DO J=1,JM 
      CONT1=0.   
      SUM1=0.    
      DO  I=1,IO
      WFCIJL=0.   
      PEARTH=C3LAND(I,J)-C3LICE(I,J)
      CONT1=CONT1+PEARTH
      SUM1=SUM1+PEARTH*WFCIJL
      DO  K=1,8
      WFCIJL=WFCIJL+VDATA(I,J,K)*VADATA(K,4,3)
      ENDDO     ! K
      SUM1=SUM1+PEARTH*WFCIJL
      ENDDO     ! I
      IF (CONT1.LE.0.) GO TO 865            
      SUM1=SUM1/CONT1
      VMASK(J)=SUM1
 865  CONTINUE
      ENDDO     ! J
c     print *,' VMASK form NP to SP in meters of water'
c     print '(12f7.2,/11f7.2)',(VMASK(jm-j+1),j=1,JM)
C ************
      return
      END                                                               1923.   
