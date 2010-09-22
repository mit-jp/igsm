
#include "ctrparam.h"

!	==========================================================
!
!	VWEIGAV.F: 
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

      SUBROUTINE VWEIGHAV (GRAV,VIN,VOUT,PHI,PHIS,ZSTARAV,LB,LE,JB,JE
     *  ,PTOP,SCLH,P,SIGE,DSIGO,CKS,CKN,IO,IM,JM,LM)                    2595.5
C                                                                       2596.
      DIMENSION P(IM,JM),SCLH(JM,LM),SIGE(37),DSIGO(35),
     &     PHIS(IO,JM)                                                  2596.5
      DIMENSION VIN(JM,LM),VOUT(JM),PHI(IM,JM,LM),ZSTARAV(JM)  
c     JHALF=12                                                          2597.5
      JHALF=JM/2
      DO 10 J=JB,JE                                                     2598.
      COEF=.5                                                           2598.5
      IF(J.EQ.JB.OR.J.EQ.JE) COEF=1.                                    2599.
      COEKD=CKS                                                         2599.5
      IF(J.GT.JHALF) COEKD=CKN                                          2600.
      ZTEM=0.                                                           2600.5
      CONT=0.                                                           2601.
      TEM1=0.                                                           2602.
      TEM2=0.                                                           2602.5
      IF(J.NE.JB) TEM1=ZSTARAV(J) 
      IF(J.NE.JE) TEM2=ZSTARAV(J+1)
      DO 20 L=LB,LE                                                     2601.5
      TEM=P(1,J)*SCLH(J,L)/(P(1,J)*SIGE(L+1)+PTOP)                      2604.
      TEM=TEM*DSIGO(L)*EXP(-.5*(PHI(1,J,L)+PHI(1,J,L+1)-2.*PHIS(1,J))/  2604.5
     *  GRAV/(COEF*(TEM1+TEM2)*COEKD+1.E-20))                           2605.
      CONT=CONT+TEM                                                     2605.5
 20   ZTEM=ZTEM+VIN(J,L)*TEM                                            2606.
      VOUT(J)=ZTEM/(CONT+1.E-20)                                        2606.5
 10   CONTINUE                                                          2607.
C                                                                       2607.5
      RETURN                                                            2608.
      END                                                               2608.5
      SUBROUTINE VWEI1AV  (GRAV,VIN,VOUT,PHI,PHIS,ZSTARAV,LB,LE,JB,JE 
     *  ,PTOP,SCLH,P,SIGE,DSIGO,CKS,CKN,IO,IM,JM,LM)                    2609.5
C                                                                       2610.
      DIMENSION P(IM,JM),SCLH(JM,LM),SIGE(37),DSIGO(35),PHIS(IO,JM)     2610.5
      DIMENSION VIN(JM,LM),VOUT(JM),PHI(IM,JM,LM),ZSTARAV(JM) 
c     JHALF=12                                                          2611.5
      JHALF=JM/2
      DO 10 J=JB,JE                                                     2612.
      COEKD=CKS                                                         2612.5
      IF(J.GT.JHALF) COEKD=CKN                                          2613.
      ZTEM=0.                                                           2613.5
      CONT=0.                                                           2614.
      DO 20 L=LB,LE                                                     2614.5
      TEM=(P(1,J)+P(1,J-1))*.5*(SCLH(J,L)+SCLH(J-1,L))/                 2615.
     *  ((P(1,J)+P(1,J-1))*SIGE(L+1)+2.*PTOP)                           2615.5
      TEM=TEM*DSIGO(L)*EXP(-.25*(PHI(1,J,L)+PHI(1,J-1,L)+PHI(1,J,L+1)+  2616.
     *  PHI(1,J-1,L+1)-2.*PHIS(1,J)-2.*PHIS(1,J-1))/GRAV/               2616.5
     *  (ZSTARAV(J)*COEKD+1.E-20))        
      CONT=CONT+TEM                                                     2617.5
 20   ZTEM=ZTEM+VIN(J,L)*TEM                                            2618.
      VOUT(J)=ZTEM/(CONT+1.E-20)                                        2618.5
 10   CONTINUE                                                          2619.
C                                                                       2619.5
      RETURN                                                            2620.
      END                                                               2620.5
