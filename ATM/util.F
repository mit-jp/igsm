
#include "ctrparam.h"

!	==========================================================
!
!	UTIL.F:  Some utility functions for the climate model.
!
!	----------------------------------------------------------
!
!	Revision History:
!	
!	When	Who		What
!	-----	----------	-------	
!	080200	Chien Wang	repack based on CliChem3 & M24x11,
!				  and add cpp.
!
!	==========================================================

C CMS SYSTEM ROUTINES EMULATION FOR IBM RS/6000
C
      SUBROUTINE CLOCKS(IHSC)
C THIS VERSION OF CLOCKS RETURNS PROCESS TIME OF USER AND
C SYSTEM TIME OF CHILD PROCESSES
C NOTE: MCLOCK IS REALLY IN HUNDREDTHS OF A SECOND, NOT SIXTIETHS.
CCC      IHSC=-MCLOCK()
      logical first
      real *4 zero,a
      data first /.true./
      if(first) then
      zero=0.0
      a=secnds(zero)
      first=.false.
      end if
      IHSC=100*secnds(a)
c     IHSC=0.
      RETURN
      END
      FUNCTION THBAR (X,Y)
      REAL   A,B,C,D,E,F,G,Q,AL
CC    DOUBLE PRECISION A,B,C,D,E,F,G,Q,AL
      DATA A,B,C,D,E,F,G/113.4977618974100,438.5012518098521,
     *  88.49964112645850,-11.50111432385882,
     *  30.00033943846368,299.9975118132485,299.9994728900967/
      Q=X/Y
      AL=(A+Q*(B+Q*(C+Q*(D+Q))))/(E+Q*(F+G*Q))
      THBAR=X*AL
      RETURN
      END
!
      FUNCTION EXPBYK (X)
      EXPBYK=X**.286
      RETURN
      END
!
      FUNCTION EXPBYKOLD (X)
C     EXPBYK=X**.286
c     DOUBLE PRECISION  A(7),B(7),C(7),D(7),E(7),F(7),G(7),H(7)
      REAL A(7),B(7),C(7),D(7),E(7),F(7),G(7),H(7)
c     DOUBLE PRECISION  TOP, BOT
C
      DATA A(1)     /.3910084705257427D12/
      DATA B(1)     /.1323236271112985D11/
      DATA C(1)     /.4866245535199495D8/
      DATA D(1)     /.2825751070482957D5/
      DATA E(1)     /.2021679763023094D12/
      DATA F(1)     /.3219813576002469D10/
      DATA G(1)     /.7026939414893149D7/
      DATA H(1)     /.2245905505347945D4/
CUT1     DC    E'.272E3'
      DATA A(2)     /.1527376839478999D10/
      DATA B(2)     /.2067556673365934D9/
      DATA C(2)     /.3041403459557123D7/
      DATA D(2)     /.7064377675658254D4/
      DATA E(2)     /.1173982318753656D10/
      DATA F(2)     /.7478937617614235D8/
      DATA G(2)     /.6528830353471617D6/
      DATA H(2)     /.8346812273642854D3/
CUT2     DC    E'.68E2'
      DATA A(3)     /.5966315773653637D7/
      DATA B(3)     /.3230557302541690D7/
      DATA C(3)     /.1900877162187466D6/
      DATA D(3)     /.1766094418795829D4/
      DATA E(3)     /.6817273980203066D7/
      DATA F(3)     /.1737197094305671D7/
      DATA G(3)     /.6066030070508742D5/
      DATA H(3)     /.3102057273820339D3/
CUT3     DC    E'.17E2'
      DATA A(4)     /.2330592097816283D5/
      DATA B(4)     /.5047745784578189D5/
      DATA C(4)     /.1188048226473864D5/
      DATA D(4)     /.4415236046666724D3/
      DATA E(4)     /.3958766989520573D5/
      DATA F(4)     /.4035136939373286D5/
      DATA G(4)     /.5636035678286606D4/
      DATA H(4)     /.1152866390241182D3/
CUT4     DC    E'.425E1'
      DATA A(5)     /.9103875388060500D2/
      DATA B(5)     /.7887102788893237D3/
      DATA C(5)     /.7425301414937417D3/
      DATA D(5)     /.1103809011715425D3/
      DATA E(5)     /.2298842049606823D3/
      DATA F(5)     /.9372759240517918D3/
      DATA G(5)     /.5236521711877195D3/
      DATA H(5)     /.4284578897126435D2/
CUT5     DC    E'.10625E1'
      DATA A(6)     /.3556201323540056D0/
      DATA B(6)     /.1232359810733637D2/
      DATA C(6)     /.4640813384330962D2/
      DATA D(6)     /.2759522529901992D2/
      DATA E(6)     /.1334929482333115D1/
      DATA F(6)     /.2177091313247190D2/
      DATA G(6)     /.4865327546221521D2/
      DATA H(6)     /.1592345521834502D2/
CUT6     DC    E'.265625E0'
      DATA A(7)     /.1389141141676475D-2/
      DATA B(7)     /.1925562204367069D0/
      DATA C(7)     /.2900508365211162D1/
      DATA D(7)     /.6898806323349460D1/
      DATA E(7)     /.7751888493733211D-2/
      DATA F(7)     /.5056917033590845D0/
      DATA G(7)     /.4520445715351951D1/
      DATA H(7)     /.5917884392240980D1/
C
C
      IF(X.LT.272.)  GO TO 10
      K=1
      GO TO 100
   10 IF(X.LT.68.)  GO TO 20
      K=2
      GO TO 100
   20 IF(X.LT.17.)  GO TO 30
      K=3
      GO TO 100
   30 IF(X.LT.4.25) GO TO 40
      K=4
      GO TO 100
   40 IF(X.LT.1.0625) GO TO 50
      K=5
      GO TO 100
   50 IF(X.LT.0.265625) GO TO 60
      K=6
      GO TO 100
   60 K=7
C
  100 CONTINUE
      IF(X.LT.272.   .AND. X.GE.68.)       K=2
      IF(X.LT.68.    .AND. X.GE.17.)       K=3
      IF(X.LT.17.    .AND. X.GE.4.25)      K=4
      IF(X.LT.4.25   .AND. X.GE.1.0625)    K=5
      IF(X.LT.1.0625 .AND. X.GE.0.265625)  K=6
      IF(X.LT.0.265625)                    K=7
C
      TOP = X
      BOT = H(K)
      TOP = TOP + D(K)
      BOT = BOT * X
      TOP = TOP * X
      BOT = BOT + G(K)
      TOP = TOP + C(K)
      BOT = BOT * X
      TOP = TOP * X
      BOT = BOT + F(K)
      TOP = TOP + B(K)
      BOT = BOT * X
      TOP = TOP * X
      BOT = BOT + E(K)
      TOP = TOP + A(K)
C
      TOP = TOP / BOT
      EXPBYK = TOP
      RETURN
      END
C     SUBROUTINE READT (IUNIT,NSKIP,AIN,LENGTH,AOUT,IPOS)
C****
C**** READ IN TITLE & REAL*4 ARRAY
C****
C     REAL*4 AIN(LENGTH),X
C     REAL*4 AOUT(LENGTH)
C     CHARACTER*80 TITLE
C     DO 10 N=1,IPOS-1
C  10 READ (IUNIT,END=920)
C     READ (IUNIT,ERR=910,END=920,NUM=LEN) TITLE,(X,N=1,NSKIP),AIN
C     IF(LEN.LT.4*(20+NSKIP+LENGTH)) GO TO 930
C     DO 100 N=1,LENGTH
C 100 AOUT(N)=AIN(N)
C     WRITE(6,'('' Read from Unit '',I2,'':'',A80)') IUNIT,TITLE
C     RETURN
C 910 WRITE(6,*) 'READ ERROR ON UNIT',IUNIT
C     STOP 'READ ERROR'
C 920 WRITE(6,*) 'END OF FILE ENCOUNTERED ON UNIT',IUNIT
C     STOP 'NO DATA TO READ'
C 930 WRITE(6,*) LEN/4,' RATHER THAN',20+NSKIP+LENGTH,' WORDS ON UNIT',
C    *  IUNIT
C     STOP 'NOT ENOUGH DATA FOUND'
C     END
