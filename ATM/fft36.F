
#include "ctrparam.h"

!	==========================================================
!
!	FFT36.F:  FFT functions.
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

      SUBROUTINE FRTR0 (KM)                                                1.   
C**** INITIALIZATION ENTRY TO CALCULATE SIN VALUES AND CHECK THAT KM=36    2.   
      COMMON/FCOM/BYKM,BYKMH,BYKM2,SIN10,SIN20,SIN30,SIN40,SIN50,SIN60,    3.   
     *  SIN70,SIN80                                                        4.   
CV    REAL*8 TWOPI/6.283185307179586477/                                   5.   
c     DOUBLE PRECISION TWOPI                                               6.   
      DATA   TWOPI/6.283185307179586477/                                   7.   
      IF(KM.NE.36) GO TO 220                                               8.   
      BYKM=1./KM                                                           9.   
      BYKMH=2./KM                                                         10.   
      BYKM2=1./(2.*KM)                                                    11.   
      SIN10=DSIN(TWOPI/36.)                                               12.   
      SIN20=DSIN(TWOPI/18.)                                               13.   
      SIN30=1.D0/2.                                                       14.   
      SIN40=DSIN(TWOPI/9.)                                                15.   
      SIN50=DCOS(TWOPI/9.)                                                16.   
      SIN60=DSQRT(3.D0)/2.                                                17.   
      SIN70=DCOS(TWOPI/18.)                                               18.   
      SIN80=DCOS(TWOPI/36.)                                               19.   
CYDBL SIN10=DSIN(TWOPI/36.)                                               20.   
CYDBL SIN20=DSIN(TWOPI/18.)                                               21.   
CYDBL SIN30=1.D0/2.                                                       22.   
CYDBL SIN40=DSIN(TWOPI/9.)                                                23.   
CYDBL SIN50=DCOS(TWOPI/9.)                                                24.   
CYDBL SIN60=DSQRT(3.D0)/2.                                                25.   
CYDBL SIN70=DCOS(TWOPI/18.)                                               26.   
CYDBL SIN80=DCOS(TWOPI/36.)                                               27.   
C     SIN10=SIN(TWOPI/36.)                                                28.   
C     SIN20=SIN(TWOPI/18.)                                                29.   
C     SIN30=.5                                                            30.   
C     SIN40=SIN(TWOPI/9.)                                                 31.   
C     SIN50=COS(TWOPI/9.)                                                 32.   
C     SIN60=SQRT(3.)/2.                                                   33.   
C     SIN70=COS(TWOPI/18.)                                                34.   
C     SIN80=COS(TWOPI/36.)                                                35.   
      RETURN                                                              36.   
  220 WRITE (6,901) KM                                                    37.   
      STOP                                                                38.   
  901 FORMAT ('0THIS FOURT SUBROUTINE NOT SUITED FOR KM = ',I8)           39.   
      END                                                                 40.   
      SUBROUTINE FRTR (F)                                                101.   
C**** THIS SUBROUTINE PERFORMS A FOURIER ANALYSIS ON THE ONE DIMENSIONAL 102.   
C**** ARRAY F WHICH MUST BE DIMENSIONED 36.  IT RETURNS IN F THE ENERGY  103.   
C**** ASSOCIATED WITH EACH WAVE NUMBER.  UPON ENTERING THIS ROUTINE,     104.   
C**** THE TOTAL ENERGY IS                                                105.   
C****   .5*SUM(F(K)*F(K))                                                106.   
C**** WITH THE SUM BEING TAKEN OVER ALL K FROM 1 TO 36.  UPON LEAVING    107.   
C**** THIS ROUTINE, THE TOTAL ENERGY IS                                  108.   
C****   SUM(F(N+1))                                                      109.   
C**** WITH THE SUM BEING TAKEN OVER ALL WAVE NUMBERS FROM 0 TO 18.       110.   
      COMMON/FCOM/BYKM,BYKMH,BYKM2,SIN10,SIN20,SIN30,SIN40,SIN50,SIN60,  111.   
     *  SIN70,SIN80                                                      112.   
      DIMENSION F(36)                                                    113.   
   10 CC00=F(12)+F(24)+F(36)                                             114.   
      CC01=F(36)-(F(12)+F(24))*SIN30                                     115.   
      CC10=F(1)+F(13)+F(25)                                              116.   
      CC11=F(1)*SIN80-F(13)*SIN40-F(25)*SIN20                            117.   
      CC20=F(2)+F(14)+F(26)                                              118.   
      CC21=F(2)*SIN70-F(14)*SIN50-F(26)*SIN10                            119.   
      CC30=F(3)+F(15)+F(27)                                              120.   
      CC31=(F(3)-F(15))*SIN60                                            121.   
      CC40=F(4)+F(16)+F(28)                                              122.   
      CC41=F(4)*SIN50-F(16)*SIN70+F(28)*SIN10                            123.   
      CC50=F(5)+F(17)+F(29)                                              124.   
      CC51=F(5)*SIN40-F(17)*SIN80+F(29)*SIN20                            125.   
      CC60=F(6)+F(18)+F(30)                                              126.   
      CC61=(F(6)+F(30))*SIN30-F(18)                                      127.   
      CC70=F(7)+F(19)+F(31)                                              128.   
      CC71=F(7)*SIN20-F(19)*SIN80+F(31)*SIN40                            129.   
      CC80=F(8)+F(20)+F(32)                                              130.   
      CC81=F(8)*SIN10-F(20)*SIN70+F(32)*SIN50                            131.   
      CC90=F(9)+F(21)+F(33)                                              132.   
      CC91=(F(33)-F(21))*SIN60                                           133.   
      CCA0=F(10)+F(22)+F(34)                                             134.   
      CCA1=F(34)*SIN70-F(10)*SIN10-F(22)*SIN50                           135.   
      CCB0=F(11)+F(23)+F(35)                                             136.   
      CCB1=F(35)*SIN80-F(11)*SIN20-F(23)*SIN40                           137.   
      SC01=(F(12)-F(24))*SIN60                                           138.   
      SC11=F(1)*SIN10+F(13)*SIN50-F(25)*SIN70                            139.   
      SC21=F(2)*SIN20+F(14)*SIN40-F(26)*SIN80                            140.   
      SC31=(F(3)+F(15))*SIN30-F(27)                                      141.   
      SC41=F(4)*SIN40+F(16)*SIN20-F(28)*SIN80                            142.   
      SC51=F(5)*SIN50+F(17)*SIN10-F(29)*SIN70                            143.   
      SC61=(F(6)-F(30))*SIN60                                            144.   
      SC71=F(7)*SIN70-F(19)*SIN10-F(31)*SIN50                            145.   
      SC81=F(8)*SIN80-F(20)*SIN20-F(32)*SIN40                            146.   
      SC91=F(9)-(F(21)+F(33))*SIN30                                      147.   
      SCA1=F(10)*SIN80-F(22)*SIN40-F(34)*SIN20                           148.   
      SCB1=F(11)*SIN70-F(23)*SIN50-F(35)*SIN10                           149.   
C**** CALCULATE EXPRESSIONS SUMMED BY INCREMENTS OF 4                    150.   
      C400=CC00+CC40+CC80                                                151.   
      C401=CC01+CC41+CC81                                                152.   
      C403=CC00-(CC40+CC80)*SIN30                                        153.   
      C402=(CC01-(CC41+CC81)*SIN30)+((SC41-SC81)*SIN60)                  154.   
      C404=(CC01-(CC41+CC81)*SIN30)-((SC41-SC81)*SIN60)                  155.   
      C410=CC10+CC50+CC90                                                156.   
      C411=CC11+CC51+CC91                                                157.   
      C413=(CC10-CC50)*SIN60                                             158.   
      C412=((CC11-CC51)*SIN60)+((SC11+SC51)*SIN30-SC91)                  159.   
      C414=((CC11-CC51)*SIN60)-((SC11+SC51)*SIN30-SC91)                  160.   
      C420=CC20+CC60+CCA0                                                161.   
      C421=CC21+CC61+CCA1                                                162.   
      C423=(CC20+CCA0)*SIN30-CC60                                        163.   
      C422=((CC21+CCA1)*SIN30-CC61)+((SC21-SCA1)*SIN60)                  164.   
      C424=((CC21+CCA1)*SIN30-CC61)-((SC21-SCA1)*SIN60)                  165.   
      C430=CC30+CC70+CCB0                                                166.   
      C431=CC31+CC71+CCB1                                                167.   
      C433=(CCB0-CC70)*SIN60                                             168.   
      C432=((CCB1-CC71)*SIN60)+(SC31-(SC71+SCB1)*SIN30)                  169.   
      C434=((CCB1-CC71)*SIN60)-(SC31-(SC71+SCB1)*SIN30)                  170.   
      S401=SC01+SC41+SC81                                                171.   
      S403=(CC40-CC80)*SIN60                                             172.   
      S402=((CC41-CC81)*SIN60)+((SC41+SC81)*SIN30-SC01)                  173.   
      S404=((CC41-CC81)*SIN60)-((SC41+SC81)*SIN30-SC01)                  174.   
      S411=SC11+SC51+SC91                                                175.   
      S413=(CC10+CC50)*SIN30-CC90                                        176.   
      S412=((CC11+CC51)*SIN30-CC91)+((SC51-SC11)*SIN60)                  177.   
      S414=((CC11+CC51)*SIN30-CC91)-((SC51-SC11)*SIN60)                  178.   
      S421=SC21+SC61+SCA1                                                179.   
      S423=(CC20-CCA0)*SIN60                                             180.   
      S422=((CC21-CCA1)*SIN60)+(SC61-(SC21+SCA1)*SIN30)                  181.   
      S424=((CC21-CCA1)*SIN60)-(SC61-(SC21+SCA1)*SIN30)                  182.   
      S431=SC31+SC71+SCB1                                                183.   
      S433=CC30-(CC70+CCB0)*SIN30                                        184.   
      S432=(CC31-(CC71+CCB1)*SIN30)+((SC71-SCB1)*SIN60)                  185.   
      S434=(CC31-(CC71+CCB1)*SIN30)-((SC71-SCB1)*SIN60)                  186.   
C**** CALCULATE EXPRESSIONS SUMMED BY INCREMENTS OF 2                    187.   
      C200=C400+C420                                                     188.   
      C201=C401+C421                                                     189.   
      C202=C402+C422                                                     190.   
      C203=C403+C423                                                     191.   
      C204=C404+C424                                                     192.   
      C205=C404-C424                                                     193.   
      C206=C403-C423                                                     194.   
      C207=C402-C422                                                     195.   
      C208=C401-C421                                                     196.   
C     C209=C400-C420                                                     197.   
      C210=C410+C430                                                     198.   
      C211=C411+C431                                                     199.   
      C212=C412+C432                                                     200.   
      C213=C413+C433                                                     201.   
      C214=C414+C434                                                     202.   
      C215=S414-S434                                                     203.   
      C216=S413-S433                                                     204.   
      C217=S412-S432                                                     205.   
      C218=S411-S431                                                     206.   
C     C219=0                                                             207.   
C     S200=0                                                             208.   
      S201=S401+S421                                                     209.   
      S202=S402+S422                                                     210.   
      S203=S403+S423                                                     211.   
      S204=S404+S424                                                     212.   
      S205=S424-S404                                                     213.   
      S206=S423-S403                                                     214.   
      S207=S422-S402                                                     215.   
      S208=S421-S401                                                     216.   
C     S209=0                                                             217.   
C     S210=0                                                             218.   
      S211=S411+S431                                                     219.   
      S212=S412+S432                                                     220.   
      S213=S413+S433                                                     221.   
      S214=S414+S434                                                     222.   
      S215=C414-C434                                                     223.   
      S216=C413-C433                                                     224.   
      S217=C412-C432                                                     225.   
      S218=C411-C431                                                     226.   
C     S219=C410-C430                                                     227.   
C**** CALCULATE THE SQUARE OF THE MAGNITUDE OF G(1,N)+I*G(2,N)           228.   
   20 F(1)=(C200+C210)*(C200+C210)*BYKM2                                 229.   
      F(2)=((C201+C211)*(C201+C211)+(S201+S211)*(S201+S211))*BYKM        230.   
      F(3)=((C202+C212)*(C202+C212)+(S202+S212)*(S202+S212))*BYKM        231.   
      F(4)=((C203+C213)*(C203+C213)+(S203+S213)*(S203+S213))*BYKM        232.   
      F(5)=((C204+C214)*(C204+C214)+(S204+S214)*(S204+S214))*BYKM        233.   
      F(6)=((C205+C215)*(C205+C215)+(S205+S215)*(S205+S215))*BYKM        234.   
      F(7)=((C206+C216)*(C206+C216)+(S206+S216)*(S206+S216))*BYKM        235.   
      F(8)=((C207+C217)*(C207+C217)+(S207+S217)*(S207+S217))*BYKM        236.   
      F(9)=((C208+C218)*(C208+C218)+(S208+S218)*(S208+S218))*BYKM        237.   
      F(10)=((C400-C420)*(C400-C420)+(C410-C430)*(C410-C430))*BYKM       238.   
      F(11)=((C208-C218)*(C208-C218)+(S218-S208)*(S218-S208))*BYKM       239.   
      F(12)=((C207-C217)*(C207-C217)+(S217-S207)*(S217-S207))*BYKM       240.   
      F(13)=((C206-C216)*(C206-C216)+(S216-S206)*(S216-S206))*BYKM       241.   
      F(14)=((C205-C215)*(C205-C215)+(S215-S205)*(S215-S205))*BYKM       242.   
      F(15)=((C204-C214)*(C204-C214)+(S214-S204)*(S214-S204))*BYKM       243.   
      F(16)=((C203-C213)*(C203-C213)+(S213-S203)*(S213-S203))*BYKM       244.   
      F(17)=((C202-C212)*(C202-C212)+(S212-S202)*(S212-S202))*BYKM       245.   
      F(18)=((C201-C211)*(C201-C211)+(S211-S201)*(S211-S201))*BYKM       246.   
      F(19)=(C200-C210)*(C200-C210)*BYKM2                                247.   
      RETURN                                                             248.   
      END                                                                249.   
      SUBROUTINE GETAN (F,G)                                             301.   
C**** GETAN RETRIEVES THE FOURIER COEFFICIENTS CONTAINED IN AN           302.   
C**** ARRAY G DIMENSIONED 2 BY 19 AND DEFINED BY                         303.   
C****   G(1,N+1)+I*G(2,N+1)=SUM(F(K)*EXP(-2*PI*I*N*K/KM))/KMH            304.   
C**** WITH THE SUM TAKEN OVER ALL K FROM 1 TO KM.  KMH = KM FOR N = 0    305.   
C**** OR 18, OTHERWISE KMH = KM/2.  THE INTERNAL NOTATION CPQN MEANS     306.   
C****   CPQN = SUM(F(K)*COS(2*PI*N*K/KM))                                307.   
C**** WITH THE SUM BEING TAKEN OVER ALL K FROM 1 TO KM WHICH ARE EQUAL   308.   
C**** TO Q MODULO(P).  SPQN IS THE SAME BUT WITH COS REPLACED BY SIN.    309.   
C**** THE NOTATION A=10, B=11, ETC. IS USED FOR  P, Q AND N.             310.   
      COMMON/FCOM/BYKM,BYKMH,BYKM2,SIN10,SIN20,SIN30,SIN40,SIN50,SIN60,  311.   
     *  SIN70,SIN80                                                      312.   
      DIMENSION F(36),G(2,19)                                            313.   
C**** CALCULATE EXPRESSIONS SUMMED BY INCREMENTS OF 12                   314.   
   10 CC00=F(12)+F(24)+F(36)                                             315.   
      CC01=F(36)-(F(12)+F(24))*SIN30                                     316.   
      CC10=F(1)+F(13)+F(25)                                              317.   
      CC11=F(1)*SIN80-F(13)*SIN40-F(25)*SIN20                            318.   
      CC20=F(2)+F(14)+F(26)                                              319.   
      CC21=F(2)*SIN70-F(14)*SIN50-F(26)*SIN10                            320.   
      CC30=F(3)+F(15)+F(27)                                              321.   
      CC31=(F(3)-F(15))*SIN60                                            322.   
      CC40=F(4)+F(16)+F(28)                                              323.   
      CC41=F(4)*SIN50-F(16)*SIN70+F(28)*SIN10                            324.   
      CC50=F(5)+F(17)+F(29)                                              325.   
      CC51=F(5)*SIN40-F(17)*SIN80+F(29)*SIN20                            326.   
      CC60=F(6)+F(18)+F(30)                                              327.   
      CC61=(F(6)+F(30))*SIN30-F(18)                                      328.   
      CC70=F(7)+F(19)+F(31)                                              329.   
      CC71=F(7)*SIN20-F(19)*SIN80+F(31)*SIN40                            330.   
      CC80=F(8)+F(20)+F(32)                                              331.   
      CC81=F(8)*SIN10-F(20)*SIN70+F(32)*SIN50                            332.   
      CC90=F(9)+F(21)+F(33)                                              333.   
      CC91=(F(33)-F(21))*SIN60                                           334.   
      CCA0=F(10)+F(22)+F(34)                                             335.   
      CCA1=F(34)*SIN70-F(10)*SIN10-F(22)*SIN50                           336.   
      CCB0=F(11)+F(23)+F(35)                                             337.   
      CCB1=F(35)*SIN80-F(11)*SIN20-F(23)*SIN40                           338.   
      SC01=(F(12)-F(24))*SIN60                                           339.   
      SC11=F(1)*SIN10+F(13)*SIN50-F(25)*SIN70                            340.   
      SC21=F(2)*SIN20+F(14)*SIN40-F(26)*SIN80                            341.   
      SC31=(F(3)+F(15))*SIN30-F(27)                                      342.   
      SC41=F(4)*SIN40+F(16)*SIN20-F(28)*SIN80                            343.   
      SC51=F(5)*SIN50+F(17)*SIN10-F(29)*SIN70                            344.   
      SC61=(F(6)-F(30))*SIN60                                            345.   
      SC71=F(7)*SIN70-F(19)*SIN10-F(31)*SIN50                            346.   
      SC81=F(8)*SIN80-F(20)*SIN20-F(32)*SIN40                            347.   
      SC91=F(9)-(F(21)+F(33))*SIN30                                      348.   
      SCA1=F(10)*SIN80-F(22)*SIN40-F(34)*SIN20                           349.   
      SCB1=F(11)*SIN70-F(23)*SIN50-F(35)*SIN10                           350.   
C**** CALCULATE EXPRESSIONS SUMMED BY INCREMENTS OF 4                    351.   
      C400=CC00+CC40+CC80                                                352.   
      C401=CC01+CC41+CC81                                                353.   
      C403=CC00-(CC40+CC80)*SIN30                                        354.   
      C402=(CC01-(CC41+CC81)*SIN30)+((SC41-SC81)*SIN60)                  355.   
      C404=(CC01-(CC41+CC81)*SIN30)-((SC41-SC81)*SIN60)                  356.   
      C410=CC10+CC50+CC90                                                357.   
      C411=CC11+CC51+CC91                                                358.   
      C413=(CC10-CC50)*SIN60                                             359.   
      C412=((CC11-CC51)*SIN60)+((SC11+SC51)*SIN30-SC91)                  360.   
      C414=((CC11-CC51)*SIN60)-((SC11+SC51)*SIN30-SC91)                  361.   
      C420=CC20+CC60+CCA0                                                362.   
      C421=CC21+CC61+CCA1                                                363.   
      C423=(CC20+CCA0)*SIN30-CC60                                        364.   
      C422=((CC21+CCA1)*SIN30-CC61)+((SC21-SCA1)*SIN60)                  365.   
      C424=((CC21+CCA1)*SIN30-CC61)-((SC21-SCA1)*SIN60)                  366.   
      C430=CC30+CC70+CCB0                                                367.   
      C431=CC31+CC71+CCB1                                                368.   
      C433=(CCB0-CC70)*SIN60                                             369.   
      C432=((CCB1-CC71)*SIN60)+(SC31-(SC71+SCB1)*SIN30)                  370.   
      C434=((CCB1-CC71)*SIN60)-(SC31-(SC71+SCB1)*SIN30)                  371.   
      S401=SC01+SC41+SC81                                                372.   
      S403=(CC40-CC80)*SIN60                                             373.   
      S402=((CC41-CC81)*SIN60)+((SC41+SC81)*SIN30-SC01)                  374.   
      S404=((CC41-CC81)*SIN60)-((SC41+SC81)*SIN30-SC01)                  375.   
      S411=SC11+SC51+SC91                                                376.   
      S413=(CC10+CC50)*SIN30-CC90                                        377.   
      S412=((CC11+CC51)*SIN30-CC91)+((SC51-SC11)*SIN60)                  378.   
      S414=((CC11+CC51)*SIN30-CC91)-((SC51-SC11)*SIN60)                  379.   
      S421=SC21+SC61+SCA1                                                380.   
      S423=(CC20-CCA0)*SIN60                                             381.   
      S422=((CC21-CCA1)*SIN60)+(SC61-(SC21+SCA1)*SIN30)                  382.   
      S424=((CC21-CCA1)*SIN60)-(SC61-(SC21+SCA1)*SIN30)                  383.   
      S431=SC31+SC71+SCB1                                                384.   
      S433=CC30-(CC70+CCB0)*SIN30                                        385.   
      S432=(CC31-(CC71+CCB1)*SIN30)+((SC71-SCB1)*SIN60)                  386.   
      S434=(CC31-(CC71+CCB1)*SIN30)-((SC71-SCB1)*SIN60)                  387.   
C**** CALCULATE EXPRESSIONS SUMMED BY INCREMENTS OF 2                    388.   
      C200=C400+C420                                                     389.   
      C201=C401+C421                                                     390.   
      C202=C402+C422                                                     391.   
      C203=C403+C423                                                     392.   
      C204=C404+C424                                                     393.   
      C205=C404-C424                                                     394.   
      C206=C403-C423                                                     395.   
      C207=C402-C422                                                     396.   
      C208=C401-C421                                                     397.   
C     C209=C400-C420                                                     398.   
      C210=C410+C430                                                     399.   
      C211=C411+C431                                                     400.   
      C212=C412+C432                                                     401.   
      C213=C413+C433                                                     402.   
      C214=C414+C434                                                     403.   
      C215=S414-S434                                                     404.   
      C216=S413-S433                                                     405.   
      C217=S412-S432                                                     406.   
      C218=S411-S431                                                     407.   
C     C219=0                                                             408.   
C     S200=0                                                             409.   
      S201=S401+S421                                                     410.   
      S202=S402+S422                                                     411.   
      S203=S403+S423                                                     412.   
      S204=S404+S424                                                     413.   
      S205=S424-S404                                                     414.   
      S206=S423-S403                                                     415.   
      S207=S422-S402                                                     416.   
      S208=S421-S401                                                     417.   
C     S209=0                                                             418.   
C     S210=0                                                             419.   
      S211=S411+S431                                                     420.   
      S212=S412+S432                                                     421.   
      S213=S413+S433                                                     422.   
      S214=S414+S434                                                     423.   
      S215=C414-C434                                                     424.   
      S216=C413-C433                                                     425.   
      S217=C412-C432                                                     426.   
      S218=C411-C431                                                     427.   
C     S219=C410-C430                                                     428.   
C**** CALCULATE FINAL COEFFICIENTS OF FOURIER EXPANSION                  429.   
      G(1,1)=(C200+C210)*BYKM                                            430.   
      G(1,2)=(C201+C211)*BYKMH                                           431.   
      G(1,3)=(C202+C212)*BYKMH                                           432.   
      G(1,4)=(C203+C213)*BYKMH                                           433.   
      G(1,5)=(C204+C214)*BYKMH                                           434.   
      G(1,6)=(C205+C215)*BYKMH                                           435.   
      G(1,7)=(C206+C216)*BYKMH                                           436.   
      G(1,8)=(C207+C217)*BYKMH                                           437.   
      G(1,9)=(C208+C218)*BYKMH                                           438.   
      G(1,10)=(C400-C420)*BYKMH                                          439.   
      G(1,11)=(C208-C218)*BYKMH                                          440.   
      G(1,12)=(C207-C217)*BYKMH                                          441.   
      G(1,13)=(C206-C216)*BYKMH                                          442.   
      G(1,14)=(C205-C215)*BYKMH                                          443.   
      G(1,15)=(C204-C214)*BYKMH                                          444.   
      G(1,16)=(C203-C213)*BYKMH                                          445.   
      G(1,17)=(C202-C212)*BYKMH                                          446.   
      G(1,18)=(C201-C211)*BYKMH                                          447.   
      G(1,19)=(C200-C210)*BYKM                                           448.   
      G(2,1)=0.                                                          449.   
      G(2,2)=(S201+S211)*BYKMH                                           450.   
      G(2,3)=(S202+S212)*BYKMH                                           451.   
      G(2,4)=(S203+S213)*BYKMH                                           452.   
      G(2,5)=(S204+S214)*BYKMH                                           453.   
      G(2,6)=(S205+S215)*BYKMH                                           454.   
      G(2,7)=(S206+S216)*BYKMH                                           455.   
      G(2,8)=(S207+S217)*BYKMH                                           456.   
      G(2,9)=(S208+S218)*BYKMH                                           457.   
      G(2,10)=(C410-C430)*BYKMH                                          458.   
      G(2,11)=(S218-S208)*BYKMH                                          459.   
      G(2,12)=(S217-S207)*BYKMH                                          460.   
      G(2,13)=(S216-S206)*BYKMH                                          461.   
      G(2,14)=(S215-S205)*BYKMH                                          462.   
      G(2,15)=(S214-S204)*BYKMH                                          463.   
      G(2,16)=(S213-S203)*BYKMH                                          464.   
      G(2,17)=(S212-S202)*BYKMH                                          465.   
      G(2,18)=(S211-S201)*BYKMH                                          466.   
      G(2,19)=0.                                                         467.   
      RETURN                                                             468.   
      END                                                                469.   
