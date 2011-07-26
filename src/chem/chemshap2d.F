
#include "ctrparam.h"

!	============================================================
!
!	CHEMSHAP2D.F:   A revised version of SHAP2D.F which is a
!			subroutine for applying Shapiro (2d)
!			smoothing of MIT Global Chemistry Model
!
!	------------------------------------------------------------
!
!	Author:		Chien Wang
!			MIT Joint Program on Science and Policy
!				of Global Change
!
!	----------------------------------------------------------
! 
!	Revision History:
!	
!	When	Who		What
!	----	----------	-------
!	080494	Chien Wang	rev.	
!	080200	Chien Wang	repack based on CliChem3 & add cpp
!
!	==========================================================

      subroutine chemshap2d (MFILTR,NORDER,XXX,IM,JM,J1,ITYPE)              8590.   

      COMMON/WORK2/X1JI(72,46),X2JI(72,46),X3JI(72,46),X1(72),X2(72),
     *  X3(72),X4(72),XM1(72),XJMP1(72)                             

!	----------------------------------------------------------

#if ( defined CPL_CHEM )

C      VARIABLE ITYPE DETERMINES TYPE OF BOUNDARY CONDITIONS
C      ITYPE=1 FOR PS,T AND Q ( XM1=X2)
C      ITYPE=2 FOR U (XM1=X1)
C      ITYPE=3 FOR V (XM1=-X1)

      JMM1=JM-1                                       
      J2=J1+1                                        
      IMBY2=1              
      DO 145 N=1,NORDER   

      DO 146 K=1,IM      
      X1(K)=X1JI(K,J1)  
      X2(K)=X1JI(K,J2) 
      X3(K)=X1JI(K,JMM1)   
      X4(K)=X1JI(K,JM)    
      IF(ITYPE.EQ.1)THEN
        XM1(K)=X1JI(K,J2)
        XJMP1(K)=X1JI(K,JMM1)
      ELSEIF(ITYPE.EQ.2)THEN
        XM1(K)=X1JI(K,J1)
        XJMP1(K)=X1JI(K,JM)
      ELSE
        XM1(K)=-X1JI(K,J1)
        XJMP1(K)=-X1JI(K,JM)
      ENDIF
  146 CONTINUE

      DO 142 I=1,IM 
      X1IM1=X1JI(I,J1) 
      DO 142 J=J2,JMM1
      X1I=X1JI(I,J)  
      X1JI(I,J)=X1IM1-X1I-X1I+X1JI(I,J+1)
      X1IM1=X1I                         
 142  CONTINUE                         

      SUM1=0.                         
      SUMJM=0.                       
      DO 144 K=1,IMBY2              
ccc   SUM1 =SUM1 +X2(K)-X1(K)-X1(K)+X2(K)  
      SUM1 =SUM1 +XM1(K)-X1(K)-X1(K)+X2(K)
ccc   SUMJM=SUMJM+X3(K)-X4(K)-X4(K)+X3(K) 
      SUMJM=SUMJM+X3(K)-X4(K)-X4(K)+XJMP1(K)
  144 CONTINUE                             

      X1SUM     =SUM1 /IMBY2              
      XJMSUM    =SUMJM/IMBY2             
c      DO 147 K=1,IM                     
c      X1JI(K,JM)=XJMSUM                
c  147 X1JI(K,J1)= X1SUM               

 145  CONTINUE      

      DO 160 I=1,IM
c      DO 160 J=J1,JM
      do 160 j=j2,jmm1
      X1JI(I,J)=(X3JI(I,J)-X1JI(I,J)/XXX)
  160 CONTINUE                          

#endif

      RETURN                           
      END                             
