
#include "ctrparam.h"

!	==========================================================
!
!	DD2G04.F:  Diagnostic program for 2-D model.
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

C**** DD2G04 BD2G04 DD2G04                 01/02/93                        0.1  
C**** OPT(3)                                                               0.2  
C****                                                                      0.3  
C**** Diagnostic program for 2-D model.                                    0.4  
C**** Like D2G04, but run on work station.                                 0.5  
c     SUBROUTINE DIAGA (U,V,T,P,Q,TC)                                      1.   
      SUBROUTINE DIAGA (U,V,T,P,Q,NOCLM) 
C****                                                             IDACC    2.   
C**** CONTENTS OF AJ(J,N)  (SUM OVER LONGITUDE AND TIME OF)                3.   
C****   1  SRINCP0 (W/M**2)                                        2 RD    4.   
C****   2  SRNFP0 (W/M**2)                                         2 RD    5.   
C****   3  SRNFP1 (W/M**2)                                         2 RD    6.   
C****   4  SRABSATM=AJ(2)-AJ(6) (W/M**2)                           2 D1    7.   
C****   5  SRINCG (W/M**2)                                         2 RD    8.   
C****   6  SRNFG (W/M**2)                                          2 RD    9.   
C****   7  TRNFP0=AJ(74)+A2BYA1*AJ(9)/DTSRCE (W/M**2)              2 D1   10.   
C****   8  TRNFP1=AJ(75)+A2BYA1*AJ(9)/DTSRCE (W/M**2)              2 D1   11.   
C****   9  TRHDT (J/M**2)                                          1 SF   12.   
C****  10  RNFP0=AJ(2)+AJ(7) (W/M**2)                              2 D1   13.   
C****  11  RNFP1=AJ(3)+AJ(8) (W/M**2)                              2 D1   14.   
C****  12  RHDT=A1BYA2*AJ(6)*DTSRCE+AJ(9) (J/M**2)                 1 D1   15.   
C****  13  QDT (J/M**2)                                            1 SF   16.   
C****  14  EVHDT (J/M**2)                                          1 SF   17.   
C****  15  F2DT (J/M**2)                                           1 GD   18.   
C****  16  HEATZ1=AJ(41)+AJ(42)                                    1 D1   19.   
C****  17  TG2 (K-273.16)                                          1 GD   20.   
C****  18  TG1 (K-273.16)                                          1 GD   21.   
C****  19  EVAP (KG/M**2)                                          1 GD   22.   
C****  20  PRCP=AJ(61)+AJ(62) (100 PA)                             1 D1   23.   
C****  21  TX (K-273.16)  (INTEGRAL OVER ATMOSPHERE OF)            4 DA   24.   
C****  22  TX1 (K-273.16)                                          4 DA   25.   
C****  23  TS (K-273.16)                                           3 SF   26.   
C****  24  DTH/DPHI  (STRATOSPHERE)                                4 DA   27.   
C****  25  DTH/DPHI  (TROPOSPHERE)                                 4 DA   28.   
C****  26  .0625*DTH*DLNP/(DU*DU+DV*DV)  (STRATOSPHERE)            4 DA   29.   
C****  27  .0625*DTH*DLNP/(DU*DU+DV*DV)  (TROPOSPHERE)             4 DA   30.   
C****  28  4*UMAX/(DX*SINJ)  (STRATOSPHERE) replaced               4 DA   31.   
C****  29  4*UMAX/(DX*SINJ)  (TROPOSPHERE)  replaced               4 DA   32.   
C****  28  Wsfr variance
C****  29  Sea level pressure
C****  30  POICE (1)                                               1 GD   33.   
C****  31  PSNOW (1)                                               4 DA   34.   
c      32 TAUL
c      33 TAUF
C****  34  TG3 ocean only
C****  35  T(J+1)-T(J-1)  (SUM OVER STRATOSPHERE OF)               4 DA   36.   
C****  36  T(J+1)-T(J-1)  (SUM OVER TROPOSPHERE OF)                4 DA   37.   
C****  37  SQRT(DTH/DLNP)/SINJ  (STRATOSPHERE)  replaced           4 DA   38.   
C****  38  SQRT(DTH/DLNP)/SINJ  (TROPOSPHERE)   replaced           4 DA   39.   
C****  37  Surf Wind
C****  38  Total stress
C****  39  ENERGP (J/M**2)                                         1 PR   40.   
C****  40  ERUN1 (J/M**2)                                          1 GP   41.   
C****  41  EDIFS (J/M**2)                                          1 GP   42.   
C****  42  F1DT (J/M**2)                                           1 GD   43.   
C****  43  ERUN2 (J/M**2)                                          1 GP   44.   
C****  44  HEATZ0=AJ(12)+AJ(13)+AJ(14)+AJ(39)-AJ(40) (J/M**2)      1 D1   45.   
C****  45  DIFS (KG/M**2)                                          1 GP   46.   
C****  46  DWTR1=AJ(20)*SCALE(20)-(AJ(19)+AJ(45)+AJ(54))*SCALE(19) 1 D1   47.   
C****  47  RUN2 (KG/M**2)                                          1 GP   48.   
C****  48  DWTR2=AJ(45)-AJ(47) (KG/M**2)                           1 D1   49.   
C****  49  WTR1 (KG/M**2)                                          1 GD   50.   
C****  50  ACE1 (KG/M**2)                                          1 GD   51.   
C****  51  WTR2 (KG/M**2)                                          1 GD   52.   
C****  52  ACE2 (KG/M**2)                                          1 GD   53.   
C****  53  SNOW (KG/M**2)                                          1 GD   54.   
C****  54  RUN1 (KG/M**2)                                          1 GP   55.   
C****  55  BTEMPW-TF                                               2 RD   56.   
C****  56  HEATZ2=AJ(15)+AJ(43) (J/M**2)                           1 D1   57.   
C****  57  PCLDSS (1)  (COMPOSITE OVER ATMOSPHERE)                 2 RD   58.   
C****  58  PCLDMC (1)  (COMPOSITE OVER ATMOSPHERE)                 2 RD   59.   
C****  59  PCLD (1)  (COMPOSITE OVER ATMOSPHERE)                   2 RD   60.   
C****  60  CLDTOPMC=AJ(80)/AJ(58) (100 PA)                         0 D1   61.   
C****  61  PRCPSS (100 PA)                                         1 CN   62.   
C****  62  PRCPMC (100 PA)                                         1 CN   63.   
C****  63  Q*P (100 PA)  (INTEGRAL OVER ATMOSPHERE OF)             4 DA   64.   
C****  64  GAM  (K/M)  (*SIG(TROPOSPHERE)/GRAV)                    4 DA   65.   
C****  65  GAMM  (K-S**2/M**2)  (SIG(TROPOSPHERE)/GAMD)            4 DA   66.   
C****  66  GAMC  (K/M)                                             4 DA   67.   
C****  67  TRINCG (W/M**2)                                         2 RD   68.   
C****  68-69  FREE                                                        69.   
C****  70  TRNFP0-TRNFG (W/M**2)                                   2 RD   70.   
C****  71  TRNFP1-TRNFG (W/M**2)                                   2 RD   71.   
C****  72  PLAVIS*S0*COSZ (W/M**2)                                 2 RD   72.   
C****  73  PLANIR*S0*COSZ (W/M**2)                                 2 RD   73.   
C****  74  ALBVIS*S0*COSZ (W/M**2)                                 2 RD   74.   
C****  75  ALBNIR*S0*COSZ (W/M**2)                                 2 RD   75.   
C****  76  SRRVIS*S0*COSZ (W/M**2)                                 2 RD   76.   
C****  77  SRRNIR*S0*COSZ (W/M**2)                                 2 RD   77.   
C****  78  SRAVIS*S0*COSZ (W/M**2)                                 2 RD   78.   
C****  79  SRANIR*S0*COSZ (W/M**2)                                 2 RD   79.   
C****  80  PBOTMC-PTOPMC (100 PA)                                  2 RD   80.   
C****                                                                     81.   
C**** CONTENTS OF APJ(J,N)  (SUM OVER LONGITUDE AND TIME OF)              82.   
C****   1  P (100 PA)                                              4 DA   83.   
C****   2  4*P4I (100 PA)  (UV GRID)                               4 DA   84.   
C****                                                                     85.   
C**** CONTENTS OF AJL(J,L,N)  (SUM OVER LONGITUDE AND TIME OF)            86.   
C****   1  (TX-273.16)*P (100 K*PA)                                4 DA   87.   
C****   2  PHI*P (100 N/S**2)                                      4 DA   88.   
C****   3  Q*P (100 PA)                                            4 DA   89.   
C****   4  4*PU4I (100 PA*M/S)  (UV GRID)                          4 DA   90.   
C****   5  4*PV4I (100 PA*M/S)  (UV GRID)                          4 DA   91.   
C****   6  SD (100 N/S)                                            4 DA   92.   
C****   7  -(PWAI*DA-PWI*SPAI/SPI*DA)  (P GRID)                    4 DA   93.   
C****   8  FMX(MC)*P (100 PA)                                      1 CN   94.   
C****   9  SRHR (W/M**2)                                           2 RD   95.   
C****  10  TRHR (W/M**2)                                           2 RD   96.   
C****  11  DTX(SS)*P (100 K*PA)                                    1 CN   97.   
C****  12  DT(DC)*P                                                  C3   98.   
C****  13  DT(MC)*P (100 PA*K)  DRY HEATING                          C3   99.   
C****  14  4*(PU4I*PU4I+PV4I*PV4I)/P4I (100 N/S**2)  (UV GRID)     4 DA  100.   
C****  15  4*PWW4I (100 N/S**2)  (UV GRID)                         4 DA  101.   
C****  16  (TH*SQRT(P)-THGM)**2/GMEAN(PR**(1-KAPA)*DTH/DPR)        4 DA  102.   
C****  17  TH*P                                                    4 DA  103.   
C****  18  RH*P                                                    4 DA  104.   
C****  19  PCLD*P (TOTAL)                                            C3  105.   
C****  20  16*PT16I*PV4I/P4I (100 PA*K*M/S)  (UV GRID)             4 DA  106.   
C****  21  16*PTV16I (100 PA*K*M/S)  (UV GRID)                     4 DA  107.   
C****  22  16*PZ16I*PV4I/P4I (100 W/S**2)  (UV GRID)               4 DA  108.   
C****  23  16*PZV16I (100 W/S**2)  (UV GRID)                       4 DA  109.   
C****  24  16*PQ16I*PV4I/P4I (100 PA*M/S)  (UV GRID)               4 DA  110.   
C****  25  16*PQV16I (100 PA*M/S)  (UV GRID)                       4 DA  111.   
C****  26  4*PWW4I*PV4I/P4I (100 W/S**2)  (UV GRID)                4 DA  112.   
C****  27  4*PWWV4I (100 W/S**2)  (UV GRID)                        4 DA  113.   
C****  28  PCLD*P (SS)                                               C3  114.   
C****  29  PCLD*P (MC)                                               C3  115.   
C****  30  2*(SHA*T+PHI) * SDMEAN                                  4 DA  116.   
C****  31  2*(SHA*T+PHI) * SD                                      4 DA  117.   
C****  32  2*Q * SDMEAN                                            4 DA  118.   
C****  33  2*Q * SD                                                4 DA  119.   
C****  34  2*PHI * (SD-SDMEAN)                                     4 DA  120.   
C****  35  16*(U*U+V*V) * SD  (UV GRID)                            4 DA  121.   
C****  36  8*(U*SD - U*SDMEAN)  (UV GRID)                          4 DA  122.   
C****  37  8*(U+R*OMEGA*COSJ) * SD  (UV GRID)                      4 DA  123.   
C****  38  DU(DC)*P  (UV GRID)                                       GD  124.   
C****  39  DU(MC)*P (100 N/M/S)  (UV GRID)                         1 CN  125.   
C****  40  DU(ED)*P*(DTSURF*DSIG*ED/DZ**2)  (UV GRID)                SF  126.   
C****  41  U  (SUM OVER I FROM 5 TO 9)  (PV GRID) (COMMENTED OUT)  4 DA  127.   
C****  41  P*V*((TH-THMEAN) * (DU/DP) / (DTH/DP) - U+UMEAN )       4 DA  128.   
C****  42  SRNFLB (W/M**2)                                         2 RD 
C****  43  TRNFLB (W/M**2)                                         2 RD
C****  44  U  (SUM OVER I FROM 35 TO 3)  (PV GRID) (COMMENTED OUT) 4 DA  131.   
C****  44  (2F-2D(UDX))*16PV(TH-THMEAN)/(DTH/DSIG)+(SD-SDMEAN)*8U  4 DA  132.   
C****  45  SRNFLB CLEAR SKY (W/M**2)                               2 RD 
C****  46  TRNFLB CLEAR SKY (W/M**2)                               2 RD
C****  47  V-V*  =D((V-VI)*(T-TI)/DTHDP)/DP                        4 DA  135.   
C****  48  4*PU4I*PV4I/P4I (100 N/S**2)  (UV GRID)                 4 DA  136.   
C****  49  4*PUV4I (100 N/S**2)  (UV GRID)                         4 DA  137.   
C****  50  DT(MC)*P (100 PA*K)  CHANGE OF PHASE                    1 CN  138.   
C****  51  FREE                                                          139.   
C****  52  PV =DTH/DP*(DUCOS(LAT)/COS(LAT)DLAT-F)-DTH/DLAT*DU/DP   4 DA  140.   
C****  53  VDVT - VPDA4*DVT/PDA4                                   4 DA  141.   
C****  54  SIGMA  (VARIANCE FOR MOIST CONVECTION)                  4 DA  142.   
C      55  NORT. TRANSPORT of Q by hor. diff.
C      56  DQ due to hori. diff.
C      57  DQ due to MC
C****                                                                    143.   
C**** CONTENTS OF ASJL(J,L,N)  (SUM OVER LONGITUDE AND TIME OF)          144.   
C****   1  TX (K)                                                  4 DA  145.   
C****   2  PHI (M**2/S**2)                                         4 DA  146.   
C****   3  SRHR (W/M**2)                                           2 RD  147.   
C****   4  TRHR (W/M**2)                                           2 RD  148.   
C****                                                                    149.   
C**** CONTENTS OF AIJ(I,J,N)  (SUM OVER TIME OF)                         150.   
C****   1  POICE (1)                                               4 DA  151.   
C****   2  PSNOW (1)                                               4 DA  152.   
C****   3  SNOW (KG/M**2)                                          4 DA  153.   
C****   4  QDT (J/M**2)                                            1 SF  154.   
C****   5  PREC (KG/M**2)                                          1 PR  155.   
C****   6  EVAP (KG/M**2)                                          1 SF  156.   
C****   7  BETA (1)                                                3 SF  157.   
C****   8  4*P4 (100 PA)  (UV GRID)                                4 DA  158.   
C****   9  PHI1000 (M**2/S**2)                                     4 DA  159.   
C****  10  PHI850 (M**2/S**2-1500*GRAV)                            4 DA  160.   
C****  11  PHI700-3000*GRAV                                        4 DA  161.   
C****  12  PHI500-5600*GRAV                                        4 DA  162.   
C****  13  PHI300-9500*GRAV                                        4 DA  163.   
C****  14  PHI100-16400*GRAV                                       4 DA  164.   
C****  15  PHI30-24000*GRAV                                        4 DA  165.   
C****  16  T700-273.16 (K-273.16)*GRAV)                            4 DA  166.   
C****  17  PCLDMC (1)  (COMPOSITE OVER ATMOSPHERE)                 2 RD  167.   
C****  18  PBOTMC-PTOPMC (100 PA)                                  2 RD  168.   
C****  19  PCLD (1)  (COMPOSITE OVER ATMOSPHERE)                   2 RD  169.   
C****  20  16*P4*(SHA*T4+Z4)*V1*DSIG*DXV (100 W*M/S**2)  (UV GRID) 4 DA  170.   
C****  21  TRNFP0 (W/M**2)                                         2 RS  171.   
C****  22  SRHDT+TRHDT (J/M**2)                                    1 SF  172.   
C****  23  SRHDT+TRHDT+QDT+EVHDT+ENRGP (J/M**2)                    1 SP  173.   
C****  24  SRNFP0 (W/M**2)                                         2 RD  174.   
C****  25  SRINCP0 (W/M**2)                                        2 RD  175.   
C****  26  SRNFG (W/M**2)                                          2 RD  176.   
C****  27  SRINCG (W/M**2)                                         2 RD  177.   
C****  28  TG1 (K-273.16)                                          1 GD  178.   
C****  29  PLICE+PEARTH*PFROZEN (1)                                4 DA  179.   
C****  30  TG2 (K-273.16)                                          1 GD  180.   
C****  31  DTH/DPHI  (TROPOSPHERE)                                 4 DA  181.   
C****  32  RUN1 (KG/M**2)                                          1 SF  182.   
C****  33  TS (K-273.16)  (USING LAPSE RATE FROM TX1)              4 DA  183.   
C****  34  CDM (1)                                                 3 SF  184.   
C****  35  TS (K-273.16)                                           3 SF  185.   
C****  36  US (M/S)                                                3 SF  186.   
C****  37  VS (M/S)                                                3 SF  187.   
C****  38  PSL (100 PA-1000)  (USING TS)                           4 DA  188.   
C****  39  UJET (M/S)                                              4 DA  189.   
C****  40  VJET (M/S)                                              4 DA  190.   
C****  41  PCLD(LOW) (1)                                           2 RD  191.   
C****  42  PCLD(MID) (1)                                           2 RD  192.   
C****  43  PCLD(HIGH) (1)                                          2 RD  193.   
C****  44  BTEMPW-TF (K-273.16)                                    2 RD  194.   
C****  45  PLAVIS*S0*COSZ (W/M**2)                                 2 RD  195.   
C****  46  ALPHA0 (1)                                              1 SF  196.   
C****                                                                    197.   
C**** CONTENTS OF AIL(I,L,N)  (SUM OVER TIME OF)                         198.   
C**** WE ARE NOT TAKING INTO ACCOUNT THE VARIATION OF MASS               199.   
C****   1  U (M/S)  (SUM OVER J FROM 11 TO 13)  (PU GRID)          4 DA  200.   
C****   2  V (M/S)  (SUM OVER J FROM 11 TO 13)  (PU GRID)          4 DA  201.   
C****   3  SD (100 N/S)  (SUM OVER J FROM 11 TO 13)                4 DA  202.   
C****   4  TX (K-273.16)  (SUM OVER J FROM 11 TO 13)               4 DA  203.   
C****   5  RH (1)  (SUM OVER J FROM 11 TO 13)                      4 DA  204.   
C****   6  DTX(MC)*P*DA (100 K*N)  (SUM OVER J FROM 11 TO 13)      1 CN  205.   
C****   7  (SRHR+TRHR)*DA (W)  (SUM OVER J FROM 11 TO 13)          2 RD  206.   
C****   9  SD (100 N/S)  (AT J=19)                                 4 DA  207.   
C****  10  TX-273.16  (AT J = 19)                                  4 DA  208.   
C****  11  SR+TR  (AT J = 19)                                      2 RD  209.   
C****  12  2*U  (AT J=19)                                          4 DA  210.   
C****  13  SD  (AT J = 21)                                         4 DA  211.   
C****  14  TX-273.16  (AT J = 21)                                  4 DA  212.   
C****  15  SR+TR  (AT J = 21)                                      2 RD  213.   
C****  16  2*U  (AT J=21)                                          4 DA  214.   
C****                                                                    215.   
C**** CONTENTS OF AIJL(I,J,L,N)  (SUM OVER TIME OF)                      216.   
C****   1  4*P4*U1 (100 PA*M/S)  (UV GRID)                         4 DA  217.   
C****   2  4*P4*V1 (100 PA*M/S)  (UV GRID)                         4 DA  218.   
C****   3  16*P4*(SHA*T4+Z4) (100 N/S**2)  (UV GRID)               4 DA  219.   
C****                                                                    220.   
C**** CONTENTS OF IDACC(N), NUMBER OF ACCUMULATION TIMES OF              221.   
C****   1  SOURCE TERMS  (DETERMINED BY NDYN)                            222.   
C****   2  RADIATION SOURCE TERMS  (DETERMINED BY NRAD)                  223.   
C****   3  SURFACE INTERACTION SOURCE TERMS  (DETERMINED BY NDASF)       224.   
C****   4  QUANTITIES IN DIAGA  (DETERMINED BY NDAA)                     225.   
C****   5  ENERGY NUMBERS IN DIAG4  (DEYERMINED BY NDA4)                 226.   
C****   6  KINETIC ENERGY IN DIAG5 FROM DYNAMICS  (DETERMINED BY NDA5K)  227.   
C****   7  ENERGY IN DIAG5 FROM DYNAMICS  (DETERMINED BY NDA5D)          228.   
C****   8  ENERGY IN DIAG5 FROM SOURCES  (DETERMINED BY NDA5S)           229.   
C****   9  WAVE ENERGY IN DIAG7  (EVERY 12 HOURS)                        230.   
C****  10  ENERGY IN DIAG5 FROM FILTER  (DETERMINED BY NFILTR)           231.   
C****                                                                    232.   

#include "BD2G04.COM"                                                     233.   

      COMMON/SPEC2/KMT,KINC,COEK                                         233.1  
      DIMENSION TSZM(JM0),THSZM(JM0),                                    233.2  
     *  TEMZ(JM0,LM0),COEKY(JM0),DTZDY(JM0,LM0),DTZDZ(JM0,LM0)
     & ,TA(JM0,LM0),PA(JM0)                                              233.3  
      COMMON/EPARA/VTH(JM0,LM0),WTH(JM0,LM0),VU(JM0,LM0),VV(JM0,LM0),
     &  DQSDT(JM0,LM0)                                                   233.4  
     *  ,DWV(JM0),PHIT(JM0,LM0),TPRIM2(JM0,LM0),WU(JM0,LM0),CKS,CKN      233.5  
     *  ,WQ(JM0,LM0),VQ(JM0,LM0)                                         233.6  
c     COMMON/HDFLUX/VQHD(JM0,LM0)
      COMMON/HDFLUX/VQHD(JM0,LM0),VTHD(JM0,LM0),VUHD(JM0,LM0),
     &             VVHD(JM0,LM0)
      LOGICAL POLE                                                       234.   
      LOGICAL NOCLM
      COMMON/WORK1/PIT(IM0,JM0),SD(IM0,JM0,1)                            235.   
      COMMON/WORK2/PK(IM0,JM0,LM0),W(IM0,JM0,LM0),PHIE(IM0,JM0,LM0-1),
c    &  FM(36,9),                                                        236.   
     *  GMEAN(36),THJL(46,36),THSQJL(46,36),SDMEAN(46,35),               237.   
     *  UA(36),SQRTP(72),PDA(72),TRI(3),                                 238.   
     *   PDAN(72)                                                        239.   
      COMMON/WORK3/PHI(IM0,JM0,LM0),TX(IM0,JM0,LM0),THSEC(72),PSEC(72),  240.   
     *  D2SIG(36),SHETH(36)                                              241.   
      DIMENSION LUPA(36),LDNA(36)                                        242.   
      COMMON/WORK5/DUT(IO0,JM0,LM0),DVT(IO0,JM0,LM0)                     243.   
      CHARACTER*16 TITLE                                                 244.   
      DIMENSION PMB(7),GHT(7)                                            245.   
      DATA PMB/1000.,850.,700.,500.,300.,100.,30./                       246.   
      DATA GHT/0.,1500.,3000.,5600.,9500.,16400.,24000./                 247.   
      DATA IFIRST/1/                                                     248.   
C**** QSAT=(RAIR/RVAPOR)*6.1071*EXP((L/RVAPOR)*(1/TF-1/T))/P             249.   
      QSAT(TM,PR,QL)=3.797915*EXP(QL*(7.93252E-6-2.166847E-3/TM))/PR     250.   
      CALL CLOCKS (MBEGIN)                                               251.   
      IDACC(4)=IDACC(4)+1                                                252.   
      IF(IFIRST.NE.1) GO TO 50                                           253.   
      IFIRST=0                                                           254.   
C**** INITIALIZE CERTAIN QUANTITIES                                      255.   
      P1000=1000.**KAPA                                                  255.5  
      HLAT=LHE                                                           255.6  
      DQDTX=.622*HLAT/RGAS                                               255.7  
      JET=LTM                                                            256.   
      BYIM=1./FIM                                                        257.   
      SHA=RGAS/KAPA                                                      258.   
      BETA=.0065                                                         259.   
      BBYG=BETA/GRAV                                                     260.   
      RBBYG=RGAS*BETA/GRAV                                               261.   
      GBYRB=GRAV/(RGAS*BETA)                                             262.   
      EPSLON=1.                                                          263.   
      PTOPK=EXPBYK(PTOP)                                                 264.   
      KM=0                                                               265.   
      DO 5 K=1,7                                                         266.   
      IF(PTOP.GT.PMB(K)) GO TO 6                                         267.   
    5 KM=KM+1                                                            268.   
    6 JEQ=2.+.5*JMM1                                                     269.   
      PRQ1=.75*PTOP                                                      270.   
      DLNP12=DLOG(.75/.35)                                               271.   
      DLNP23=DLOG(.35/.1)                                                272.   
      DO 10 L=1,LM                                                       273.   
      LUPA(L)=L+1                                                        274.   
   10 LDNA(L)=L-1                                                        275.   
      LDNA(1)=1                                                          276.   
      LUPA(LM)=LM                                                        277.   
   50 CONTINUE                                                           278.   
c     print *,' PRINT FROM DIAGA'
C****                                                                    279.   
C**** FILL IN HUMIDITY AND SIGMA DOT ARRAYS AT THE POLES                 280.   
C****                                                                    281.   
C****                                                                    286.   
C**** CALCULATE PK AND TX, THE REAL TEMPERATURE                          287.   
C****                                                                    288.   
      DO 80 L=1,LM                                                       289.   
      PK(1,1,L)=EXPBYK(SIG(L)*P(1,1)+PTOP)                               290.   
      TX(1,1,L)=T(1,1,L)*PK(1,1,L)                                       291.   
      PK(1,JM,L)=EXPBYK(SIG(L)*P(1,JM)+PTOP)                             292.   
      TX(1,JM,L)=T(1,JM,L)*PK(1,JM,L)                                    293.   
      DO 80 J=2,JMM1                                                     301.   
      DO 80 I=1,IM                                                       302.   
      PK(I,J,L)=EXPBYK(SIG(L)*P(I,J)+PTOP)                               303.   
   80 TX(I,J,L)=T(I,J,L)*PK(I,J,L)                                       304.   
C C C C C C C C C C C C C C C                                            304.5  
      JMM2=JM-2                                                          304.51 
      IF(COEK.EQ.0.) GO TO 89                                            304.52 
      DO 74 J=1,JM                                                       304.53 
      PA(J)=0.                                                           304.54 
      DO 71 I=1,IM                                                       304.55 
 71   PA(J)=PA(J)+P(I,J)                                                 304.56 
      DO 73 L=1,LM                                                       304.57 
      TA(J,L)=0.                                                         304.58 
      DO 72 I=1,IM                                                       304.59 
 72   TA(J,L)=TA(J,L)+T(I,J,L)                                           304.6  
 73   TA(J,L)=TA(J,L)/FIM                                                304.61 
 74   PA(J)=PA(J)/FIM                                                    304.62 
 89   CONTINUE                                                           304.63 
C****                                                                    305.   
C**** OUTSIDE J LOOP FOR ALL PRIMARY GRID ROWS                           306.   
C****                                                                    307.   
      DO 490 J=1,JM                                                      308.   
      JM1=J-1                                                            309.   
      IF(J.EQ.1) JM1=1                                                   310.   
      JP1=J+1                                                            311.   
      IF(J.EQ.JM) JP1=JM                                                 312.   
      POLE=.FALSE.                                                       313.   
      IF(J.EQ.1.OR.J.EQ.JM) POLE=.TRUE.                                  314.   
      IMAX=IM                                                            315.   
      IF(POLE) IMAX=1                                                    316.   
      CMAX=IMAX                                                          316.5  
      DXYPJ=DXYP(J)                                                      317.   
C**** NUMBERS ACCUMULATED FOR A SINGLE LEVEL                             318.   
      AT1=0.                                                             319.   
      BT1=0.                                                             320.   
      CT1=0.                                                             321.   
      BSCOV=0.                                                           322.   
      CSCOV=0.                                                           323.   
      SPI=0.                                                             324.   
      DO 120 I=1,IMAX                                                    325.   
      JR=J
      PLAND=FDATA(I,J,2)                                                 327.   
      POICE=ODATA(I,J,2)*(1.-PLAND)                                      328.   
      PLICE=FDATA(I,J,3)*PLAND                                           329.   
      POCEAN=(1.-PLAND)-POICE                                            330.   
      PEARTH=PLAND-PLICE                                                 331.   
      AT1=AT1+(TX(I,J,1)-273.16)*POCEAN                                  332.   
      BT1=BT1+(TX(I,J,1)-273.16)*PLAND                                   333.   
      CT1=CT1+(TX(I,J,1)-273.16)*POICE                                   334.   
C                                                                        335.   
      SCOVL=0.                                                           336.   
      IF(GDATA(I,J,2).GT.0.) SCOVL=PEARTH                                337.   
      IF(GDATA(I,J,12).GT.0.) SCOVL=SCOVL+PLICE                          338.   
      BSCOV=BSCOV+SCOVL                                                  339.   
      SCOVOI=0.                                                          340.   
      IF(GDATA(I,J,1).GT.0.) SCOVOI=POICE                                341.   
      CSCOV=CSCOV+SCOVOI                                                 342.   
C                                                                        343.   
      SPI=SPI+P(I,J)                                                     344.   
      AIJ(I,J,1)=AIJ(I,J,1)+POICE                                        345.   
      AIJ(I,J,2)=AIJ(I,J,2)+(SCOVOI+SCOVL)                               346.   
      AIJ(I,J,3)=AIJ(I,J,3)+(GDATA(I,J,1)*POICE+GDATA(I,J,2)*PLAND)      347.   
      TS=TX(I,J,1)*((P(I,J)+PTOP)/(SIG(1)*P(I,J)+PTOP))**RBBYG           348.   
C     AIJ(I,J,8)=AIJ(I,J,8)+((P(I,J)+PTOP)*(1.+BBYG*FDATA(I,J,1)/TS)     349.   
C    *  **GBYRB-1000.)                                                   350.   
      AIJ(I,J,29)=AIJ(I,J,29)+(PLICE+PEARTH*GDATA(I,J,6)/                351.   
     *  (GDATA(I,J,5)+GDATA(I,J,6)+1.E-20))                              352.   
C     AIJ(I,J,33)=AIJ(I,J,33)+(TS-273.16)                                353.   
      AIJ(I,J,38)=AIJ(I,J,38)+((P(I,J)+PTOP)*(1.+BBYG*FDATA(I,J,1)/      354.   
     *  BLDATA(I,J,2))**GBYRB-1000.)                                     355.   
  120 CONTINUE                                                           356.   
      AJ(J,22)=AJ(J,22)+AT1                                              357.   
      BJ(J,22)=BJ(J,22)+BT1                                              358.   
      CJ(J,22)=CJ(J,22)+CT1                                              359.   

      if(NOCLM)then
      BJ(J,31)=BJ(J,31)+BSCOV                                            360.   
      CJ(J,31)=CJ(J,31)+CSCOV                                            361.   
      endif

      APJ(J,1)=APJ(J,1)+SPI                                              362.   
C**** GEOPOTENTIALS CALCULATED FOR EACH LAYER                            363.   
  150 DO 160 I=1,IMAX                                                    364.   
      P1=SIG(1)*P(I,J)+PTOP                                              365.   
      PUP=SIG(2)*P(I,J)+PTOP                                             366.   
      IF(ABS(TX(I,J,2)-TX(I,J,1)).LT.EPSLON) GO TO 152                   367.   
      BBYGV=DLOG(TX(I,J,1)/TX(I,J,2))/(RGAS*DLOG(P1/PUP))                368.   
      PHI(I,J,1)=FDATA(I,J,1)+TX(I,J,1)                                  369.   
     *  *(((P(I,J)+PTOP)/P1)**(RGAS*BBYGV)-1.)/BBYGV                     370.   
      PHI(I,J,2)=PHI(I,J,1)+(TX(I,J,1)-TX(I,J,2))/BBYGV                  371.   
      GO TO 154                                                          372.   
  152 PHI(I,J,1)=FDATA(I,J,1)+RGAS*TX(I,J,1)*DLOG((P(I,J)+PTOP)/P1)      373.   
      PHI(I,J,2)=PHI(I,J,1)+RGAS*.5*(TX(I,J,1)+TX(I,J,2))*DLOG(P1/PUP)   374.   
  154 DO 160 L=3,LM                                                      375.   
      PDN=PUP                                                            376.   
      PUP=SIG(L)*P(I,J)+PTOP                                             377.   
      IF(ABS(TX(I,J,L)-TX(I,J,L-1)).LT.EPSLON) GO TO 156                 378.   
      BBYGV=DLOG(TX(I,J,L-1)/TX(I,J,L))/(RGAS*DLOG(PDN/PUP))             379.   
      PHI(I,J,L)=PHI(I,J,L-1)+(TX(I,J,L-1)-TX(I,J,L))/BBYGV              380.   
      GO TO 160                                                          381.   
  156 PHI(I,J,L)=PHI(I,J,L-1)+RGAS*.5*(TX(I,J,L-1)+TX(I,J,L))            382.   
     *  *DLOG(PDN/PUP)                                                   383.   
  160 CONTINUE                                                           384.   
C**** CALCULATE GEOPOTENTIAL HEIGHTS AT SPECIFIC MILLIBAR LEVELS         389.   
C**** NUMBERS ACCUMULATED FOR EACH ODD LEVEL                             414.   
      DO 230 L=1,LM                                                      415.   
      ATX=0.                                                             416.   
      BTX=0.                                                             417.   
      CTX=0.                                                             418.   
      TPI=0.                                                             419.   
      AQ=0.                                                              420.   
      BQ=0.                                                              421.   
      CQ=0.                                                              422.   
      PHIPI=0.                                                           423.   
      QPI=0.                                                             424.   
      THPI=0.                                                            425.   
      RHPI=0.                                                            426.   
      DO 220 I=1,IMAX                                                    427.   
      JR=J
      PLAND=FDATA(I,J,2)                                                 429.   
      POICE=ODATA(I,J,2)*(1.-PLAND)                                      430.   
      POCEAN=(1.-PLAND)-POICE                                            431.   
      SP=P(I,J)                                                          432.   
      ATX=ATX+(TX(I,J,L)-273.16)*POCEAN                                  433.   
      BTX=BTX+(TX(I,J,L)-273.16)*PLAND                                   434.   
      CTX=CTX+(TX(I,J,L)-273.16)*POICE                                   435.   
      AQ=AQ+Q(I,J,L)*SP*POCEAN                                           436.   
      BQ=BQ+Q(I,J,L)*SP*PLAND                                            437.   
      CQ=CQ+Q(I,J,L)*SP*POICE                                            438.   
C                                                                        439.   
C                                                                        440.   
      TPI=TPI+(TX(I,J,L)-273.16)*SP                                      441.   
      PHIPI=PHIPI+PHI(I,J,L)*SP                                          442.   
      QPI=QPI+Q(I,J,L)*SP                                                443.   
      THPI=THPI+T(I,J,L)*SP                                              444.   
      if(TX(I,J,L).gt.273.)then
      QLH=LHE                                                            445.   
      else
      QLH=LHS
      endif
      QSATL=QSAT(TX(I,J,L),SIG(L)*SP+PTOP,QLH)                           446.   
      IF(QSATL.GT.1.) QSATL=1.                                           447.   
      RHPI=RHPI+Q(I,J,L)*SP/QSATL                                        448.   
c     if(J.eq.7.and.(L.eq.6.or.L.eq.2))then
c     print *,' from DIAGA'
c     print *,' J=',J,' L=',L
c     print *,TX(I,J,L),Q(I,J,L)
c     call SAT(TX(I,J,L),SIG(L)*SP+PTOP,QSN,DN)
c     print *,Q(I,J,L)/QSATL,Q(I,J,L)/QSN
c     endif
  220 CONTINUE                                                           449.   
      AJ(J,21)=AJ(J,21)+ATX*DSIG(L)                                      450.   
      BJ(J,21)=BJ(J,21)+BTX*DSIG(L)                                      451.   
      CJ(J,21)=CJ(J,21)+CTX*DSIG(L)                                      452.   
      AJ(J,63)=AJ(J,63)+AQ*DSIG(L)                                       453.   
      BJ(J,63)=BJ(J,63)+BQ*DSIG(L)                                       454.   
      CJ(J,63)=CJ(J,63)+CQ*DSIG(L)                                       455.   
      AJL(J,L,1)=AJL(J,L,1)+TPI                                          456.   
      AJL(J,L,2)=AJL(J,L,2)+PHIPI                                        457.   
      AJL(J,L,3)=AJL(J,L,3)+QPI                                          458.   
      AJL(J,L,17)=AJL(J,L,17)+THPI                                       459.   
      AJL(J,L,18)=AJL(J,L,18)+RHPI                                       460.   
      AJL(J,L,54)=AJL(J,L,54)+TPRIM2(J,L)*PK(1,J,L)**2.                  460.1  
  230 CONTINUE                                                           461.   
C****                                                                    462.   
C**** NORTHWARD GRADIENT OF TEMPERATURE: TROPOSPHERIC AND STRATOSPHERIC  463.   
C****                                                                    464.   
C**** MEAN TROPOSPHERIC NORTHWARD TEMPERATURE GRADIENT                   465.   
      DO 340 L=1,LTM                                                     466.   
      ADTDL=0.                                                           467.   
      BDTDL=0.                                                           468.   
      CDTDL=0.                                                           469.   
      DO 335 I=1,IMAX                                                    470.   
      PLAND=FDATA(I,J,2)                                                 471.   
      POICE=ODATA(I,J,2)*(1.-PLAND)                                      472.   
      POCEAN=(1.-PLAND)-POICE                                            473.   
      ADTDL=ADTDL+(TX(I,JP1,L)-TX(I,JM1,L))*POCEAN                       474.   
      BDTDL=BDTDL+(TX(I,JP1,L)-TX(I,JM1,L))*PLAND                        475.   
      CDTDL=CDTDL+(TX(I,JP1,L)-TX(I,JM1,L))*POICE                        476.   
  335 CONTINUE                                                           477.   
      IF(.NOT.POLE) GO TO 338                                            478.   
      DO 336 I=2,IM                                                      479.   
      ADTDL=ADTDL+(TX(I,JP1,L)-TX(I,JM1,L))*POCEAN                       480.   
      BDTDL=BDTDL+(TX(I,JP1,L)-TX(I,JM1,L))*PLAND                        481.   
      CDTDL=CDTDL+(TX(I,JP1,L)-TX(I,JM1,L))*POICE                        482.   
  336 CONTINUE                                                           483.   
      ADTDL=ADTDL*2./FIM                                                 484.   
      BDTDL=BDTDL*2./FIM                                                 485.   
      CDTDL=CDTDL*2./FIM                                                 486.   
  338 AJ(J,36)=AJ(J,36)+ADTDL*DSIG(L)                                    487.   
      BJ(J,36)=BJ(J,36)+BDTDL*DSIG(L)                                    488.   
  340 CJ(J,36)=CJ(J,36)+CDTDL*DSIG(L)                                    489.   
C**** MEAN STRATOSPHERIC NORTHWARD TEMPERATURE GRADIENT                  490.   
      IF (LS1.GT.LM) GO TO 380                                           491.   
      DO 370 L=LS1,LM                                                    492.   
      ADTDL=0.                                                           493.   
      BDTDL=0.                                                           494.   
      CDTDL=0.                                                           495.   
      DO 350 I=1,IM                                                      496.   
      PLAND=FDATA(I,J,2)                                                 497.   
      POICE=ODATA(I,J,2)*(1.-PLAND)                                      498.   
      POCEAN=(1.-PLAND)-POICE                                            499.   
      ADTDL=ADTDL+(TX(I,JP1,L)-TX(I,JM1,L))*POCEAN                       500.   
      BDTDL=BDTDL+(TX(I,JP1,L)-TX(I,JM1,L))*PLAND                        501.   
      CDTDL=CDTDL+(TX(I,JP1,L)-TX(I,JM1,L))*POICE                        502.   
  350 CONTINUE                                                           503.   
      IF(.NOT.POLE) GO TO 360                                            504.   
      ADTDL=ADTDL*2./FIM                                                 505.   
      BDTDL=BDTDL*2./FIM                                                 506.   
      CDTDL=CDTDL*2./FIM                                                 507.   
  360 AJ(J,35)=AJ(J,35)+ADTDL*DSIG(L)                                    508.   
      BJ(J,35)=BJ(J,35)+BDTDL*DSIG(L)                                    509.   
  370 CJ(J,35)=CJ(J,35)+CDTDL*DSIG(L)                                    510.   
  380 CONTINUE                                                           511.   
C****                                                                    512.   
C**** STATIC STABILITIES: TROPOSPHERIC AND STRATOSPHERIC                 513.   
C****                                                                    514.   
C**** OLD TROPOSPHERIC STATIC STABILITY                                  515.   
      ASS=0.                                                             516.   
      BSS=0.                                                             517.   
      CSS=0.                                                             518.   
      DO 390 I=1,IMAX                                                    519.   
      JR=J
      PLAND=FDATA(I,J,2)                                                 521.   
      POICE=ODATA(I,J,2)*(1.-PLAND)                                      522.   
      POCEAN=(1.-PLAND)-POICE                                            523.   
      SS=(T(I,J,LTM)-T(I,J,1))/(PHI(I,J,LTM)-PHI(I,J,1)+1.E-12)          524.   
      ASS=ASS+SS*POCEAN                                                  525.   
      BSS=BSS+SS*PLAND                                                   526.   
      CSS=CSS+SS*POICE                                                   527.   
C                                                                        528.   
  390 AIJ(I,J,31)=AIJ(I,J,31)+SS                                         529.   
      AJ(J,25)=AJ(J,25)+ASS                                              530.   
      BJ(J,25)=BJ(J,25)+BSS                                              531.   
      CJ(J,25)=CJ(J,25)+CSS                                              532.   
C**** OLD STRATOSPHERIC STATIC STABILITY                                 533.   
      ASS=0.                                                             534.   
      BSS=0.                                                             535.   
      CSS=0.                                                             536.   
      DO 440 I=1,IMAX                                                    537.   
      JR=J
      PLAND=FDATA(I,J,2)                                                 539.   
      POICE=ODATA(I,J,2)*(1.-PLAND)                                      540.   
      POCEAN=(1.-PLAND)-POICE                                            541.   
      SS=(T(I,J,LM)-T(I,J,LTM))/((PHI(I,J,LM)-PHI(I,J,LTM))+1.E-5)       542.   
      ASS=ASS+SS*POCEAN                                                  543.   
      BSS=BSS+SS*PLAND                                                   544.   
      CSS=CSS+SS*POICE                                                   545.   
C                                                                        546.   
  440 CONTINUE                                                           547.   
      AJ(J,24)=AJ(J,24)+ASS                                              548.   
      BJ(J,24)=BJ(J,24)+BSS                                              549.   
      CJ(J,24)=CJ(J,24)+CSS                                              550.   
C****                                                                    551.   
C**** NUMBERS ACCUMULATED FOR THE RADIATION EQUILIBRIUM LAYERS           552.   
C****                                                                    553.   
      DO 470 LR=1,3                                                      554.   
      TRI(LR)=0.                                                         555.   
      DO 460 I=1,IMAX                                                    556.   
  460 TRI(LR)=TRI(LR)+RQT(I,J,LR)                                        557.   
  470 ASJL(J,LR,1)=ASJL(J,LR,1)+(TRI(LR)-273.16*IMAX)                    558.   
      PHIRI=0.                                                           559.   
      DO 480 I=1,IMAX                                                    560.   
  480 PHIRI=PHIRI+(PHI(I,J,LM)+RGAS*.5*(TX(I,J,LM)+RQT(I,J,1))           561.   
     *  *DLOG((SIG(LM)*P(I,J)+PTOP)/PRQ1))                               562.   
      ASJL(J,1,2)=ASJL(J,1,2)+PHIRI                                      563.   
      PHIRI=PHIRI+RGAS*.5*(TRI(1)+TRI(2))*DLNP12                         564.   
      ASJL(J,2,2)=ASJL(J,2,2)+PHIRI                                      565.   
      PHIRI=PHIRI+RGAS*.5*(TRI(2)+TRI(3))*DLNP23                         566.   
      ASJL(J,3,2)=ASJL(J,3,2)+PHIRI                                      567.   
  490 CONTINUE                                                           568.   
C****                                                                    569.   
C**** OUTSIDE J LOOP FOR NON-POLAR PRIMARY GRID ROWS                     570.   
C****                                                                    571.   
      DO 550 J=2,JMM1                                                    572.   
      DXYPJ=DXYP(J)                                                      573.   
      ROSSX=DYP(J)/(DXYPJ*SINP(J))                                       574.   
      ELX=1./SINP(J)                                                     575.   
      PTORQX=COSP(J)*DYP(J)/DXYPJ                                        576.   
C**** NUMBERS ACCUMULATED EXCEPT AT THE POLES                            577.   
      ARIT=0.                                                            578.   
      BRIT=0.                                                            579.   
      CRIT=0.                                                            580.   
      AUTMAX=0.                                                          581.   
      BUTMAX=0.                                                          582.   
      CUTMAX=0.                                                          583.   
      ALT=0.                                                             584.   
      BLT=0.                                                             585.   
      CLT=0.                                                             586.   
      ARIS=0.                                                            587.   
      BRIS=0.                                                            588.   
      CRIS=0.                                                            589.   
      AUSMAX=0.                                                          590.   
      BUSMAX=0.                                                          591.   
      CUSMAX=0.                                                          592.   
      ALS=0.                                                             593.   
      BLS=0.                                                             594.   
      CLS=0.                                                             595.   
      IM1=IM                                                             596.   
      DO 520 I=1,IM                                                      597.   
      JR=J
      PLAND=FDATA(I,J,2)                                                 599.   
      POICE=ODATA(I,J,2)*(1.-PLAND)                                      600.   
      POCEAN=(1.-PLAND)-POICE                                            601.   
      SP=P(I,J)                                                          602.   
C**** NUMBERS ACCUMULATED OVER THE TROPOSPHERE FOR ODD LEVELS            603.   
      UMAX=0.                                                            604.   
      UL=0.                                                              605.   
      DO 500 L=1,LTM                                                     606.   
      UA(L)=U(IM1,J+1,L)+U(I,J+1,L)+U(IM1,J,L)+U(I,J,L)                  607.   
      UAMAX=ABS(UA(L))                                                   608.   
      IF(UAMAX.GT.UMAX) UMAX=UAMAX                                       609.   
      UL=UL+UA(L)*DSIG(L)                                                610.   
  500 CONTINUE                                                           611.   
      DTH=T(I,J,LTM)-T(I,J,1)                                            612.   
      IF(DTH.LT.0.) DTH=0.                                               613.   
      DLNTH=DLOG(T(I,J,LTM)/T(I,J,1))                                    613.5  
      IF (DLNTH.LT.0.) DLNTH=0.                                          613.6  
      DLNP=DLOG((SIG(1)*SP+PTOP)/(SIG(LTM)*SP+PTOP))                     614.   
      DU=UA(LTM)-UA(1)                                                   615.   
      DV=(V(IM1,J+1,LTM)+V(I,J+1,LTM)+V(IM1,J,LTM)+V(I,J,LTM))-          616.   
     *  (V(IM1,J+1,1)+V(I,J+1,1)+V(IM1,J,1)+V(I,J,1))                    617.   
      RI=DTH*DLNP/(DU*DU+DV*DV+1.E-6)                                    618.   
      ARIT=ARIT+RI*POCEAN                                                619.   
      BRIT=BRIT+RI*PLAND                                                 620.   
      CRIT=CRIT+RI*POICE                                                 621.   
C                                                                        622.   
      AUTMAX=AUTMAX+UMAX*POCEAN                                          623.   
      BUTMAX=BUTMAX+UMAX*PLAND                                           624.   
      CUTMAX=CUTMAX+UMAX*POICE                                           625.   
C                                                                        626.   
      EL=SQRT(DLNTH/DLNP)                                                627.   
      ALT=ALT+EL*POCEAN                                                  628.   
      BLT=BLT+EL*PLAND                                                   629.   
      CLT=CLT+EL*POICE                                                   630.   
C                                                                        631.   
      UMAX=0.                                                            632.   
      IF (LS1.GT.LM) GO TO 512 !  NEEDED FOR RUNS WITHOUT A STRATOSPHERE 633.   
      DO 510 L=LS1,LM                                                    634.   
      UA(L)=U(IM1,J+1,L)+U(I,J+1,L)+U(IM1,J,L)+U(I,J,L)                  635.   
      UAMAX=ABS(UA(L))                                                   636.   
      IF(UAMAX.GT.UMAX) UMAX=UAMAX                                       637.   
  510 UL=UL+UA(L)*DSIG(L)                                                638.   
  512 CONTINUE                                                           639.   
      DTH=T(I,J,LM)-T(I,J,LTM)                                           640.   
      IF(DTH.LT.0.) DTH=0.                                               640.5  
      DLNTH=DLOG(T(I,J,LM)/T(I,J,LTM))                                   641.   
      IF (DLNTH.LT.0.) DLNTH=0.                                          641.5  
      DLNP=DLOG((SIG(LTM)*SP+PTOP)/(SIG(LM)*SP+PTOP))                    642.   
      DU=UA(LM)-UA(LTM)                                                  643.   
      DV=(V(IM1,J+1,LM)+V(I,J+1,LM)+V(IM1,J,LM)+V(I,J,LM))-              644.   
     *  (V(IM1,J+1,LTM)+V(I,J+1,LTM)+V(IM1,J,LTM)+V(I,J,LTM))            645.   
      RI=DTH*DLNP/(DU*DU+DV*DV+1.E-6)                                    646.   
      ARIS=ARIS+RI*POCEAN                                                647.   
      BRIS=BRIS+RI*PLAND                                                 648.   
      CRIS=CRIS+RI*POICE                                                 649.   
C                                                                        650.   
      AUSMAX=AUSMAX+UMAX*POCEAN                                          651.   
      BUSMAX=BUSMAX+UMAX*PLAND                                           652.   
      CUSMAX=CUSMAX+UMAX*POICE                                           653.   
C                                                                        654.   
      EL=SQRT(DLNTH/(DLNP))                                              655.   
      ALS=ALS+EL*POCEAN                                                  656.   
      BLS=BLS+EL*PLAND                                                   657.   
      CLS=CLS+EL*POICE                                                   658.   
C                                                                        659.   
  520 IM1=I                                                              660.   
!     AJ(J,27)=AJ(J,27)+ARIT                                             661.   
!     BJ(J,27)=BJ(J,27)+BRIT                                             662.   
!     CJ(J,27)=CJ(J,27)+CRIT                                             663.   
c     AJ(J,29)=AJ(J,29)+AUTMAX*ROSSX                                     664.   
c     BJ(J,29)=BJ(J,29)+BUTMAX*ROSSX                                     665.   
c     CJ(J,29)=CJ(J,29)+CUTMAX*ROSSX                                     666.   
c     AJ(J,38)=AJ(J,38)+ALT*ELX                                          667.   
c     BJ(J,38)=BJ(J,38)+BLT*ELX                                          668.   
c     CJ(J,38)=CJ(J,38)+CLT*ELX                                          669.   
c     AJ(J,26)=AJ(J,26)+ARIS                                             670.   
c     BJ(J,26)=BJ(J,26)+BRIS                                             671.   
c     CJ(J,26)=CJ(J,26)+CRIS                                             672.   
c     AJ(J,28)=AJ(J,28)+AUSMAX*ROSSX                                     673.   
c     BJ(J,28)=BJ(J,28)+BUSMAX*ROSSX                                     674.   
c     CJ(J,28)=CJ(J,28)+CUSMAX*ROSSX                                     675.   
c     AJ(J,37)=AJ(J,37)+ALS*ELX                                          676.   
c     BJ(J,37)=BJ(J,37)+BLS*ELX                                          677.   
c     CJ(J,37)=CJ(J,37)+CLS*ELX                                          678.   
  550 CONTINUE                                                           679.   
C****                                                                    680.   
C**** MEAN TROPOSPHERIC LAPSE RATES:  MOIST CONVECTIVE, ACTUAL,          681.   
C****    DRY ADIABATIC                                                   682.   
C****                                                                    683.   
      QLH=LHE                                                            684.   
      X=RGAS*LHE*LHE/(SHA*461.5)                                         685.   
      DO 570 J=1,JM                                                      686.   
      IMAX=IM                                                            687.   
      IF(J.EQ.1 .OR. J.EQ.JM) IMAX=IM                                    688.   
      DO 570 L=1,LTM                                                     689.   
      AGAMM=0.                                                           690.   
      BGAMM=0.                                                           691.   
      CGAMM=0.                                                           692.   
      AGAMX=0.                                                           693.   
      BGAMX=0.                                                           694.   
      CGAMX=0.                                                           695.   
      DO 560 I=1,IMAX                                                    696.   
      PLAND=FDATA(I,J,2)                                                 697.   
      POICE=ODATA(I,J,2)*(1.-PLAND)                                      698.   
      POCEAN=(1.-PLAND)-POICE                                            699.   
      PRT=(SIG(L)*P(I,J)+PTOP)*RGAS*TX(I,J,L)                            700.   
      ESEPS=QSAT(TX(I,J,L),1.,QLH)                                       701.   
      GAM=(PRT+LHE*ESEPS)/(PRT+X*ESEPS/TX(I,J,L))                        702.   
      AGAMM=AGAMM+GAM*POCEAN                                             703.   
      BGAMM=BGAMM+GAM*PLAND                                              704.   
      CGAMM=CGAMM+GAM*POICE                                              705.   
      IF(L.EQ.1) GO TO 560                                               706.   
      GAM=(TX(I,J,L-1)-TX(I,J,L))/(PHI(I,J,L)-PHI(I,J,L-1))              707.   
      AGAMX=AGAMX+GAM*POCEAN                                             708.   
      BGAMX=BGAMX+GAM*PLAND                                              709.   
      CGAMX=CGAMX+GAM*POICE                                              710.   
  560 CONTINUE                                                           711.   
      AJ(J,65)=AJ(J,65)+DSIG(L)*AGAMM                                    712.   
      BJ(J,65)=BJ(J,65)+DSIG(L)*BGAMM                                    713.   
      CJ(J,65)=CJ(J,65)+DSIG(L)*CGAMM                                    714.   
      IF(L.EQ.1) GO TO 570                                               715.   
      AJ(J,64)=AJ(J,64)+DSIGO(L-1)*AGAMX                                 716.   
      BJ(J,64)=BJ(J,64)+DSIGO(L-1)*BGAMX                                 717.   
      CJ(J,64)=CJ(J,64)+DSIGO(L-1)*CGAMX                                 718.   
  570 CONTINUE                                                           719.   
C**** DRY ADIABATIC LAPSE RATE                                           720.   
      GAMD=.0098                                                         721.   
      DO 600 J=2,JMM1                                                    722.   
      X=SINP(J)*GRAV/(COSP(J)*RGAS*2.*DLAT)                              723.   
      AGAMC=0.                                                           724.   
      BGAMC=0.                                                           725.   
      CGAMC=0.                                                           726.   
      DO 590 I=1,IM                                                      727.   
      PLAND=FDATA(I,J,2)                                                 728.   
      POICE=ODATA(I,J,2)*(1.-PLAND)                                      729.   
      POCEAN=(1.-PLAND)-POICE                                            730.   
      DT2=0.                                                             731.   
      T2=0.                                                              732.   
      DO 580 L=1,LTM                                                     733.   
      DT2=DT2+DSIG(L)*(TX(I,J+1,L)-TX(I,J-1,L))                          734.   
  580 T2=T2+DSIG(L)*TX(I,J,L)                                            735.   
      GAM=GAMD+X*DT2/T2                                                  736.   
      AGAMC=AGAMC+GAM*POCEAN                                             737.   
      BGAMC=BGAMC+GAM*PLAND                                              738.   
  590 CGAMC=CGAMC+GAM*POICE                                              739.   
      AJ(J,66)=AJ(J,66)+AGAMC                                            740.   
      BJ(J,66)=BJ(J,66)+BGAMC                                            741.   
  600 CJ(J,66)=CJ(J,66)+CGAMC                                            742.   
C****                                                                    743.   
C**** MOMENTUM, KINETIC ENERGY, NORTHWARD TRANSPORTS, ANGULAR MOMENTUM   744.   
C****                                                                    745.   
      DO 640 J=2,JM                                                      746.   
      P4I=0.                                                             747.   
      I=IM                                                               748.   
      DO 610 IP1=1,IM                                                    749.   
      P4=P(I,J-1)+P(IP1,J-1)+P(I,J)+P(IP1,J)                             750.   
      P4I=P4I+P4                                                         751.   
      AIJ(I,J,8)=AIJ(I,J,8)+P4                                           752.   
      AIJ(I,J,39)=AIJ(I,J,39)+U(I,J,JET)                                 753.   
      AIJ(I,J,40)=AIJ(I,J,40)+V(I,J,JET)                                 754.   
  610 I=IP1                                                              755.   
      APJ(J,2)=APJ(J,2)+P4I                                              756.   
      DO 640 L=1,LM                                                      757.   
      PU4I=0.                                                            758.   
      PV4I=0.                                                            759.   
      PWW4I=0.                                                           760.   
      PT16I=0.                                                           761.   
      PTV16I=0.                                                          762.   
      PZ16I=0.                                                           763.   
      PZV16I=0.                                                          764.   
      PQ16I=0.                                                           765.   
      PQV16I=0.                                                          766.   
      PWWV4I=0.                                                          767.   
      PUV4I=0.                                                           768.   
      I=IM                                                               769.   
      DO 620 IP1=1,IM                                                    770.   
      P4=P(I,J-1)+P(IP1,J-1)+P(I,J)+P(IP1,J)                             771.   
      PU4I=PU4I+P4*U(I,J,L)                                              772.   
      PV4I=PV4I+P4*V(I,J,L)                                              773.   
      PWW4I=PWW4I+P4*(U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))               774.   
      PWWV4I=PWWV4I+P4*(U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))*V(I,J,L)    775.   
      PUV4I=PUV4I+P4*U(I,J,L)*V(I,J,L)                                   776.   
      T4=TX(I,J-1,L)+TX(IP1,J-1,L)+TX(I,J,L)+TX(IP1,J,L)                 777.   
      PT16I=PT16I+P4*T4                                                  778.   
      PTV16I=PTV16I+P4*T4*V(I,J,L)                                       779.   
      Z4=PHI(I,J-1,L)+PHI(IP1,J-1,L)+PHI(I,J,L)+PHI(IP1,J,L)             780.   
      PZ16I=PZ16I+P4*Z4                                                  781.   
      PZV16I=PZV16I+P4*Z4*V(I,J,L)                                       782.   
      Q4=Q(I,J-1,L)+Q(IP1,J-1,L)+Q(I,J,L)+Q(IP1,J,L)                     783.   
      PQ16I=PQ16I+P4*Q4                                                  784.   
      PQV16I=PQV16I+P4*Q4*V(I,J,L)                                       785.   
      AIJ(I,J,20)=AIJ(I,J,20)+P4*(SHA*T4+Z4)*V(I,J,L)*DSIG(L)*DXV(J)     786.   
         AIJL(I,J,L,1)=AIJL(I,J,L,1)+P4*U(I,J,L)                         787.   
         AIJL(I,J,L,2)=AIJL(I,J,L,2)+P4*V(I,J,L)                         788.   
         AIJL(I,J,L,3)=AIJL(I,J,L,3)+P4*(SHA*T4+Z4)                      789.   
  620 I=IP1                                                              790.   
C C C C C C C C C C C C C C C                                            790.5  
      IF(COEK.EQ.0.) GO TO 609                                           790.51 
      FLUXT=VTH(J,L)*P4*4.*P1000                                         790.52 
C **********
c     FLUXSH=VQ(J,L)*P4*4.                                               790.53 
c
c       VQ transport of Q by eddies
c	VQHD transport of Q by horizontal diffusion
c
      FLUXSH=(VQ(J,L)+VQHD(J,L))*P4*4.                                   790.53 
c
C **********
      PTV16I=PTV16I+FLUXT*FIM                                            790.54 
      PQV16I=PQV16I+FLUXSH*FIM                                           790.55 
      PUV4I=PUV4I+.5*(VU(J,L)+VU(J-1,L))*P4                              790.56 
 609  CONTINUE                                                           790.57 
      AJL(J,L,4)=AJL(J,L,4)+PU4I                                         791.   
      AJL(J,L,5)=AJL(J,L,5)+PV4I                                         792.   
      AJL(J,L,14)=AJL(J,L,14)+VV(J,L)                                    793.   
      AJL(J,L,15)=AJL(J,L,15)+PWW4I                                      794.   
      AJL(J,L,20)=AJL(J,L,20)+PT16I*PV4I/P4I                             795.   
      AJL(J,L,21)=AJL(J,L,21)+PTV16I                                     796.   
      AJL(J,L,22)=AJL(J,L,22)+PZ16I*PV4I/P4I                             797.   
      AJL(J,L,23)=AJL(J,L,23)+PZV16I                                     798.   
      AJL(J,L,24)=AJL(J,L,24)+PQ16I*PV4I/P4I                             799.   
      AJL(J,L,25)=AJL(J,L,25)+PQV16I                                     800.   
c     AJL(J,L,26) is used to calculate change of KIN. EN.
c        due to eddy transport
c     AJL(J,L,26)=AJL(J,L,26)+PWW4I*PV4I/P4I                             801.   
c
      AJL(J,L,27)=AJL(J,L,27)+PWWV4I                                     802.   
      AJL(J,L,48)=AJL(J,L,48)+PU4I*PV4I/P4I                              803.   
      AJL(J,L,49)=AJL(J,L,49)+PUV4I                                      804.   
      AJL(J,L,55)=AJL(J,L,55)+VQHD(J,L)*P4*4.*FIM
  640 CONTINUE                                                           805.   
C****                                                                    806.   
      DO 655 J=1,JM                                                      807.   
      IMAX=IM                                                            808.   
      IF (J.EQ.1 .OR. J.EQ.JM) IMAX=1                                    809.   
      SPI=0.                                                             810.   
      DO 645 I=1,IMAX                                                    811.   
      SPI=SPI+P(I,J)                                                     812.   
  645 CONTINUE                                                           813.   
C****                                                                    814.   
C**** EVEN LEVEL GEOPOTENTIALS, VERTICAL WINDS AND VERTICAL TRANSPORTS   815.   
C****                                                                    816.   
      DO 655 L=1,LMM1                                                    817.   
      SDI=0.                                                             818.   
      WI=0.                                                              819.   
      PZI=0.                                                             820.   
      SDZI=0.                                                            821.   
      PDSE2I=0.                                                          822.   
      SDDS2I=0.                                                          823.   
      PQ2I=0.                                                            824.   
      SDQ2I=0.                                                           825.   
      DO 650 I=1,IMAX                                                    826.   
      JR=J
      PLAND=FDATA(I,J,2)                                                 828.   
      SDI=SDI+SD(I,J,L)                                                  829.   
      PE=SIGE(L+1)*P(I,J)+PTOP                                           830.   
      PKE=EXPBYK(PE)                                                     831.   
      THETA=THBAR(T(I,J,L+1),T(I,J,L))                                   832.   
      W(I,J,L)=SD(I,J,L)*THETA*PKE/PE                                    833.   
      PHIE(I,J,L)=PHI(I,J,L)+SHA*THETA*(PK(I,J,L)-PKE)                   834.   
      WI=WI+W(I,J,L)                                                     835.   
      PZI=PZI+PHIE(I,J,L)*P(I,J)                                         836.   
      SDZI=SDZI+PHIE(I,J,L)*SD(I,J,L)                                    837.   
      PDSE2I=PDSE2I+(SHA*(TX(I,J,L)+TX(I,J,L+1))+2.*PHIE(I,J,L))*P(I,J)  838.   
      SDDS2I=SDDS2I+(SHA*(TX(I,J,L)+TX(I,J,L+1))+2.*PHIE(I,J,L))*        839.   
     *  SD(I,J,L)                                                        840.   
      PQ2I=PQ2I+(Q(I,J,L)*Q(I,J,L+1)/(Q(I,J,L)+Q(I,J,L+1)+1.E-20))*      841.   
     *   P(I,J)                                                          842.   
      SDQ2I=SDQ2I+(Q(I,J,L)*Q(I,J,L+1)/(Q(I,J,L)+Q(I,J,L+1)+             843.   
     *  1.E-20))*SD(I,J,L)                                               844.   
  650 CONTINUE                                                           845.   
C C C C C C C C C C C C C C C                                            845.5  
      IF(COEK.EQ.0.) GO TO 289                                           845.51 
      FLUXT=WTH(J,L)*PA(J) *DXYP(J)*SHA*2.*P1000                         845.52 
      FLUXS=WQ(J,L)*PA(J)*DXYP(J)*.5                                     845.53 
      SDDS2I=SDDS2I+FLUXT*CMAX                                           845.54 
      SDQ2I=SDQ2I+FLUXS*CMAX                                             845.55 
 289  CONTINUE                                                           845.56 
      SDMEAN(J,L)=SDI*BYIM                                               846.   
      AJL(J,L,6)=AJL(J,L,6)+WI                                           847.   
      AJL(J,L,34)=AJL(J,L,34)+(SDZI-PZI*SDI/SPI)                         848.   
      AJL(J,L,30)=AJL(J,L,30)+PDSE2I*SDI/SPI                             849.   
      AJL(J,L,31)=AJL(J,L,31)+SDDS2I                                     850.   
      AJL(J,L,32)=AJL(J,L,32)+PQ2I*SDI/SPI                               851.   
      AJL(J,L,33)=AJL(J,L,33)+SDQ2I                                      852.   
  655 CONTINUE                                                           853.   
C****                                                                    854.   
C**** VERTICAL TRANSPORT OF KINETIC ENERGY AND ANGULAR MOMENTUM          855.   
C****                                                                    856.   
C**** FILL IN AND/OR DOUBLE SD AND SDMEAN AT THE POLES                   857.   
      DO 657 L=1,LMM1                                                    858.   
      SDMEAN(1,L)=2.*FIM*SDMEAN(1,L)                                     859.   
      SDMEAN(JM,L)=2.*FIM*SDMEAN(JM,L)                                   860.   
      SDSP=2.*SD(1,1,L)                                                  861.   
      SDNP=2.*SD(1,JM,L)                                                 862.   
      DO 657 I=1,IM                                                      863.   
      SD(I,1,L)=SDSP                                                     864.   
  657 SD(I,JM,L)=SDNP                                                    865.   
      DO 670 J=2,JM                                                      866.   
      AMA=RADIUS*OMEGA*COSV(J)                                           867.   
      DO 670 L=1,LMM1                                                    868.   
C     TKEM=0.                                                            869.   
      TKET=0.                                                            870.   
      UM=0.                                                              871.   
      UT=0.                                                              872.   
      I=IM                                                               873.   
      DO 660 IP1=1,IM                                                    874.   
      SDU=SD(I,J,L)+SD(IP1,J,L)+SD(I,J-1,L)+SD(IP1,J-1,L)                875.   
      UE=U(I,J,L)+U(I,J,L+1)                                             876.   
      TKE=UE*UE+(V(I,J,L)+V(I,J,L+1))*(V(I,J,L)+V(I,J,L+1))              877.   
C     TKEM=TKEM+TKE                                                      878.   
      TKET=TKET+TKE*SDU                                                  879.   
      UM=UM+UE                                                           880.   
      UT=UT+UE*SDU                                                       881.   
  660 I=IP1                                                              882.   
      AJL(J,L,35)=AJL(J,L,35)+TKET                                       883.   
      AJL(J,L,36)=AJL(J,L,36)+8.*WU(J,L)                                 884.   
      AJL(J,L,37)=AJL(J,L,37)+(UT+4*AMA*FIM*(SDMEAN(J,L)+SDMEAN(J-1,L))) 885.   
  670 CONTINUE                                                           886.   
C****                                                                    887.   
C**** AVAILABLE POTENTIAL ENERGY                                         888.   
C****                                                                    889.   
C**** SET UP FOR CALCULATION                                             890.   
      DO 710 L=1,LM                                                      891.   
  710 GMEAN(L)=0.                                                        892.   
      DO 740 J=1,JM                                                      893.   
      IMAX=IM                                                            894.   
      IF (J.EQ.1 .OR. J.EQ.JM) IMAX=1                                    895.   
      DO 720 I=1,IMAX                                                    896.   
  720 SQRTP(I)=SQRT(P(I,J))                                              897.   
C**** GMEAN CALCULATED FOR EACH LAYER, THJL, THSQJL ARRAYS FILLED        898.   
      DO 730 L=1,LM                                                      899.   
      LDN=LDNA(L)                                                        900.   
      LUP=LUPA(L)                                                        901.   
      THJL(J,L)=0.                                                       902.   
      THSQJL(J,L)=0.                                                     903.   
      DO 730 I=1,IMAX                                                    904.   
      THJL(J,L)=THJL(J,L)+T(I,J,L)*SQRTP(I)                              905.   
      THSQJL(J,L)=THSQJL(J,L)+T(I,J,L)*T(I,J,L)*P(I,J)                   906.   
  730 GMEAN(L)=GMEAN(L)+(SIG(L)*P(I,J)+PTOP)*(T(I,J,LUP)-T(I,J,LDN))*    907.   
     *  DXYP(J)/(P(I,J)*PK(I,J,L))                                       908.   
  740 CONTINUE                                                           909.   
C**** CALCULATE APE                                                      910.   
      DO 760 L=1,LM                                                      911.   
      LP1=LUPA(L)                                                        912.   
      LM1=LDNA(L)                                                        913.   
      THJL(1,L)=THJL(1,L)*FIM                                            914.   
      THJL(JM,L)=THJL(JM,L)*FIM                                          915.   
      THSQJL(1,L)=THSQJL(1,L)*FIM                                        916.   
      THSQJL(JM,L)=THSQJL(JM,L)*FIM                                      917.   
      THGM=0.                                                            918.   
      DO 750 J=1,JM                                                      919.   
  750 THGM=THGM+THJL(J,L)*DXYP(J)                                        920.   
      THGM=THGM/AREAG                                                    921.   
      GMEANL=GMEAN(L)/((SIG(LM1)-SIG(LP1))*AREAG)                        922.   
      DO 760 J=1,JM                                                      923.   
  760 AJL(J,L,16)=AJL(J,L,16)+(THSQJL(J,L)-2.*THJL(J,L)*THGM+THGM*THGM*  924.   
     *  FIM)/GMEANL                                                      925.   
C****                                                                    926.   
C**** OMEGA'*ALPHA' ;  BAROCLINIC EKE GENERATION                         927.   
C****                                                                    928.   
c     IF(JM.NE.24.OR.IM.EQ.1) GO TO 850                                  978.   
  850 CONTINUE                                                          1022.   
C****                                                                   1023.   
C**** ELIASSEN PALM FLUX                                                1024.   
C****                                                                   1025.   
C**** NORTHWARD TRANSPORT                                               1026.   
      DO 868 J=2,JM                                                     1027.   
      BYDXYV=1./DXYV(J)                                                 1028.   
      I=IM                                                              1029.   
      DO 862 IP1=1,IM                                                   1030.   
      PDA(I)=.5*((P(I,J)+P(IP1,J))*DXYS(J)+(P(I,J-1)+P(IP1,J-1))*       1031.   
     *  DXYN(J-1))                                                      1032.   
      PSEC(I)=PDA(I)*BYDXYV                                             1033.   
  862 I=IP1                                                             1034.   
      DO 868 L=1,LM                                                     1035.   
      DUDP=0.                                                           1036.   
      DTHDP=0.                                                          1037.   
      UMN=0.                                                            1038.   
      THMN=0.                                                           1039.   
      LDN=LDNA(L)                                                       1040.   
      LUP=LUPA(L)                                                       1041.   
      I=IM                                                              1042.   
      DO 864 IP1=1,IM                                                   1043.   
      DUDP=DUDP+U(I,J,LUP)-U(I,J,LDN)                                   1044.   
      DTHDP=DTHDP+T(I,J,LUP)+T(I,J-1,LUP)-T(I,J,LDN)-T(I,J-1,LDN)       1045.   
      UMN=UMN+U(I,J,L)                                                  1046.   
      THMN=THMN+T(I,J,L)+T(I,J-1,L)                                     1047.   
      THSEC(I)=T(I,J,L)+T(IP1,J,L)+T(I,J-1,L)+T(IP1,J-1,L)              1048.   
  864 I=IP1                                                             1049.   
      UMN=UMN*BYIM                                                      1050.   
      THMN=2.*THMN/FIM                                                  1051.   
      FPHI=0.                                                           1052.   
      SMALL=.0002*FIM*T(1,J,L)                                          1053.   
      IF (DTHDP.LT.SMALL) DTHDP=SMALL                                   1055.   
      DO 866 I=1,IM                                                     1056.   
  866 FPHI=FPHI+PSEC(I)*(-VU(J,L)+(VTH(J,L)+VTH(J+1,L))*P1000           1057.   
     *  *DUDP/DTHDP)                                                    1058.   
  868 AJL(J,L,41)=AJL(J,L,41)+FPHI                                      1059.   
C**** VERTICAL TRANSPORT                                                1060.   
      DO 878 J=2,JMM1                                                   1061.   
      DO 878 L=1,LMM1                                                   1062.   
      THMN=0.                                                           1063.   
      SDMN=0.                                                           1064.   
      DTHDP=0.                                                          1065.   
      DO 872 I=1,IM                                                     1066.   
      DTHDP=DTHDP+T(I,J,L+1)-T(I,J,L)                                   1067.   
      THMN=THMN+T(I,J,L+1)+T(I,J,L)                                     1068.   
  872 SDMN=SDMN+SD(I,J,L)                                               1069.   
      SMALL=.0001*FIM*T(1,J,L+1)                                        1070.   
      IF (DTHDP.LT.SMALL) DTHDP=SMALL                                   1072.   
      THMN=THMN/FIM                                                     1073.   
      SDMN=SDMN/FIM                                                     1074.   
      DUDX=0.                                                           1075.   
      PVTHP=0.                                                          1076.   
      SDPU=0.                                                           1077.   
      IM1=IM                                                            1078.   
      DO 874 I=1,IM                                                     1079.   
      DUDX=DUDX+DXV(J+1)*(U(I,J+1,L)+U(I,J+1,L+1))-DXV(J)*              1080.   
     *   (U(I,J,L)+U(I,J,L+1))                                          1081.   
      UPE=U(IM1,J,L)+U(IM1,J+1,L)+U(I,J,L)+U(I,J+1,L)+                  1082.   
     *    U(IM1,J,L+1)+U(IM1,J+1,L+1)+U(I,J,L+1)+U(I,J+1,L+1)           1083.   
      VPE=V(IM1,J,L)+V(IM1,J+1,L)+V(I,J,L)+V(I,J+1,L)+                  1084.   
     *    V(IM1,J,L+1)+V(IM1,J+1,L+1)+V(I,J,L+1)+V(I,J+1,L+1)           1085.   
      PVTHP=PVTHP+P(I,J)*8.*(VTH(J,L)+VTH(J+1,L))*P1000                 1086.   
      SDPU=SDPU+8.*WU(J,L)                                              1087.   
  874 IM1=I                                                             1088.   
      AJL(J,L,44)=AJL(J,L,44)+(.5*FIM*F(J)-.25*DUDX)*DSIGO(L)*PVTHP     1089.   
     *   /DTHDP + SDPU                                                  1090.   
  878 CONTINUE                                                          1091.   
C****                                                                   1092.   
C**** POTENTIAL VORTICITY                                               1093.   
C****                                                                   1094.   
C****                                                                   1133.   
C**** LAGRANGIAN MEAN STEAM FUNCTION                                    1134.   
C****                                                                   1135.   
C**** ACCUMULATE TIME USED IN DIAGA                                     1169.   
      CALL CLOCKS (MEND)                                                1170.   
      MINC=MBEGIN-MEND                                                  1171.   
      MDIAG=MDIAG+MINC                                                  1172.   
      MDYN=MDYN-MINC                                                    1173.   
C                                                                       1174.   
      RETURN                                                            1175.   
  997 FORMAT (' DIAGNOSTICS ACCUMULATED ',12I4,15X,2I7)                 1176.   
  999 FORMAT (' DTHETA/DP IS TOO SMALL AT J=',I4,' L=',I4,2F15.6)       1177.   
      END                                                               1178.   
      SUBROUTINE DIAG1(NOCLM)                                           1501.   
C****                                                                   1502.   
C**** THIS SUBROUTINE PRODUCES AREA WEIGHTED STATISTICS OF              1503.   
C****                                                                   1504.   
C   K   N                                                               1505.   
C****                                                                   1506.   
C***1   1  SOLAR RADIATION INCIDENT ON PLANET (W/M**2)                  1507.   
C****                                                                   1508.   
C**1A   2/1  PLANETARY ALBEDO (10**-2)                                  1509.   
C**1B  72/1  PLANETARY ALBEDO VISUAL (10**-2)                           1510.   
C**1C  73/1  PLANETARY ALBEDO NEAR IR (10**-2)                          1511.   
C**1D   6/5  GROUND ALBEDO (10**-2)                                     1512.   
C**1E  74/1  GROUND ALBEDO VISUAL (10**-2)                              1513.   
C**1F  75/1  GROUND ALBEDO NEAR IR (10**-2)                             1514.   
C**1G  76/1  ATMOSPHERIC ALBEDO VISUAL (10**-2)                         1515.   
C**1H  77/1  ATMOSPHERIC ALBEDO NEAR IR (10**-2)                        1516.   
C**1I  78/1  ATMOSPHERIC ABSORPTION VISUAL (10**-2)                     1517.   
C**1J  79/1  ATMOSPHERIC ABSORPTION NEAR IR (10**-2)                    1518.   
C****                                                                   1519.   
C***2   2  SOLAR RADIATION ABSORBED BY PLANET (W/M**2)                  1520.   
C***3   3  SOLAR RADIATION ABSORBED BELOW PTOP (W/M**2)                 1521.   
C***4   4  SOLAR RADIATION ABSORBED BY ATMOSPHERE (W/M**2)              1522.   
C***5   5  SOLAR RADIATION INCIDENT ON GROUND (W/M**2)                  1523.   
C***6   6  SOLAR RADIATION ABSORBED BY GROUND (W/M**2)                  1524.   
C***7   7  THERMAL RADIATION EMITTED BY PLANET (W/M**2)                 1525.   
C***8   8  THERMAL RADIATION AT PTOP (W/M**2)                           1526.   
C***9   9  THERMAL RADIATION EMITTED BY GROUND (W/M**2)                 1527.   
C**10  67  THERMAL RADIATION INCIDENT ON GROUND (W/M**2)                1528.   
C****                                                                   1529.   
C**11  55  BRIGHTNESS TEMPERATURE THROUGH WINDOW REGION (K-273.16)      1530.   
C****  10  NET RADIATION ABSORBED BY PLANET (W/M**2)                    1531.   
C****  11  NET RADIATION ABSORBED BELOW PTOP (W/M**2)                   1532.   
C****  12  NET RADIATION ABSORBED BY GROUND (W/M**2)                    1533.   
C****  13  SENSIBLE HEAT FLUX INTO THE GROUND (W/M**2)                  1534.   
C****  14  EVAPORATION HEAT FLUX INTO THE GROUND (W/M**2)               1535.   
C****  39  PRECIPITATION HEAT FLUX INTO THE GROUND (W/M**2)             1536.   
C****  40  HEAT RUNOFF FROM FIRST GROUND LAYER (W/M**2)                 1537.   
C****  44  NET HEATING AT Z0 (W/M**2)                                   1538.   
C****  42  CONDUCTION AT -Z1 (W/M**2)                                   1539.   
C****                                                                   1540.   
C**21  41  HEAT WATER DUFFUSION AT -Z1 (W/M**2)                         1541.   
C****  16  NET HEATING AT -Z1 (W/M**2)                                  1542.   
C****  43  HEAT RUNOFF FROM SECOND GROUND LAYER (W/M**2)                1543.   
C****  15  CONDUCTION AT -Z1-Z2 (W/M**2)                                1544.   
C****  56  NET HEATING AT -Z1-Z2 (W/M**2)                               1545.   
C****  18  MEAN TEMPERATURE OF FIRST GROUND LAYER (.1 K-273.16)         1546.   
C****  17  MEAN TEMPERATURE OF SECOND GROUND LAYER (.1 K-273.16)        1547.   
C****  23  SURFACE AIR TEMPERATURE (.1 K-273.16)                        1548.   
C****  22  FIRST LAYER AIR TEMPERATURE (.1 K-273.16)                    1549.   
C****  21  COMPOSITE AIR TEMPERATURE (.1 K-273.16)                      1550.   
C****                                                                   1551.   
C**31  35  STRATO TEMPERATURE CHANGE PER DEGREE LATITUDE (10**-2 K)     1552.   
C****  36  TROPO TEMPERATURE CHANGE PER DEGREE LATITUDE (10**-2 K)      1553.   
C****  24  STRATOSPHERIC STATIC STABILITY (10**-3 K/M)                  1554.   
C****  25  TROPOSPHERIC STATIC STABILITY (10**-3 K/M)                   1555.   
C****  26  TEMPERATURE at 2 M hight                                     1556.   
C****  27  Difference T2m from control run                              1557.   
C****  28  STRATOSPHERIC ROSSBY NUMBER (1)                              1558.   
C****  29  TROPOSPHERIC ROSSBY NUMBER (1)                               1559.   
C****  37  L IN THE STRATOSPHERE (10**5 M)                              1560.   
C****  38  L IN THE TROPOSPHERE (10**5 M)                               1561.   
C****                                                                   1562.   
C**41  64  GAM  (10**-3 K/M)                                            1563.   
C****  65  GAMM  (10**-3 K/M)                                           1564.   
C****  66  GAMC  (10**-3 K/M)                                           1565.   
C****  57  INTEGRATED SUPER-SATURATION CLOUD COVER (10**-2)             1566.   
C****  58  INTEGRATED MOIST CONVECTIVE CLOUD COVER (10**-2)             1567.   
C****  59  INTEGRATED TOTAL CLOUD COVER (10**-2)                        1568.   
C****  60  MOIST CONVECTIVE CLOUD DEPTH (100 N)                         1569.   
C****  61  SUPER SATURATION PRECIPITATION (KG/M**2/86400 S)             1570.   
C****  62  MOIST CONVECTIVE PRECIPITATION (KG/M**2/86400 S)             1571.   
C****  20  PRECIPITATION (KG/M**2/86400 S)                              1572.   
C****                                                                   1573.   
C**51  19  EVAPORATION (KG/M**2/86400 S)                                1574.   
C****  63  WATER CONTENT OF ATMOSPHERE (KG/M**2)                        1575.   
C****  45  WATER DIFFUSION AT -Z1 (KG/M**2/86400 S)                     1576.   
C****  54  WATER RUNOFF FOR FIRST GROUND LAYER (KG/M**2/86400 S)        1577.   
C****  46  NET WATER INTO THE FIRST GROUND LAYER (KG/M**2/86400 S)      1578.   
C****  47  WATER RUNOFF FOR SECOND GROUND LAYER (KG/M**2/86400 S)       1579.   
C****  48  NET WATER INTO THE SECOND GROUND LAYER (KG/M**2/86400 S)     1580.   
C****  49  WATER CONTAINED IN FIRST GROUND LAYER (KG/M**2)              1581.   
C****  50  ICE CONTAINED IN FIRST GROUND LAYER (KG/M**2)                1582.   
C****  51  WATER CONTAINED IN SECOND GROUND LAYER (KG/M**2)             1583.   
C****                                                                   1584.   
C**61  52  ICE CONTAINED IN SECOND GROUND LAYER (KG/M**2)               1585.   
C****  53  SNOW DEPTH (KG/M**2)                                         1586.   
C****  31  SNOW COVER (10**-2)                                          1587.   
C****  30  OCEAN ICE COVER (10**-2)                                     1588.   
C****                                                                   1589.   
#include "BD2G04.COM"                                                   1590.   
      LOGICAL NOCLM
      COMMON U,V,T,P,Q                                                  1591.   
      DIMENSION ABCJ(JM0,80,3)                                          1594.   
      EQUIVALENCE (AJ(1,1),ABCJ(1,1,1))                                 1595.   
C **** CLEAR SKY
      common/clrsk/CLEAR(JM0),NCLR(JM0),AJCLR(JM0,12),BJCLR(JM0,12),
     * CJCLR(JM0,12)
      integer CLEAR
C     AJCLR
C            1  SW INC AT P0            RD (AJ(1))
C            2  SW ABS BELOW P0         RD (AJ(2))
C            3  SW ABS BELOW P1         RD (AJ(3))
C            4  SW ABS AT Z0            RD (AJ(6))
C            5  SW INC AT Z0            RD (AJ(5))
C            6  LW INC AT Z0            RD (AJ(67))
C            7  NET LW AT Z0            SF (AJ(9))
C            8  NET LW AT P0            RD (AJ(7))
C            9  NET LW AT P1            RD (AJ(8))
C           10  NET RAD AT P0           DG (AJ(10))
C           11  NET RAD AT P1           DG (AJ(11))
C           12  NET RAD AT Z0           DG (AJ(12))
      dimension ABCJCL(JM0,12,3)
      EQUIVALENCE (AJCLR(1,1),ABCJCL(1,1,1))
C **** CLEAR SKY
C
      DIMENSION JLAT(46),SAREA(46),SPOCEN(46),SPOICE(46),SPLAND(46),    1596.   
     *  S1(46),MLAT(46),FLAT(46),MHEM(2),FHEM(2),WTA(4),WTB(4),WTC(4),  1597.   
     *  INDEX(70),INNUM(10),INDEN(10),IA(72),SCALE(72)                  1598.   
C****                                                                   1599.   
      CHARACTER*16 TERAIN(5)/'    (GLOBAL)','     (LAND)','    (OCEAN)',1600.   
     *  ' (OCEAN ICE)','   (REGIONS)'/                                  1601.   
      CHARACTER*16 TITLE(72),TITLE1(36),TITLE2(36),TITLEA(10)           1602.   
      EQUIVALENCE (TITLE(1),TITLE1(1)),(TITLE(37),TITLE2(1))            1603.   
                         DATA TITLE1/                                   1604.   
     1  ' INC SW(WT/M**2)', '0SW ABS BELOW P0', ' SW ABS BELOW P1',     1605.   
     4  ' SW ABS BY ATMOS', ' SW INC ON Z0   ', ' SW ABS AT Z0   ',     1606.   
     7  '0NET LW AT P0   ', ' NET LW AT P1   ', ' NET LW AT Z0  ',      1607.   
     O  '0NET RAD AT P0  ', ' NET RAD AT P1  ', ' NET RAD AT Z0  ',     1608.   
     3  '0SENSBL HEAT FLX', ' EVAPOR HEAT FLX', ' CONDC AT -Z1-Z2',     1609.   
     6  ' NET HEAT AT -Z1', ' TG2 (.1 C)     ', '0TG1 (.1 C)     ',     1610.   
     9  ' EVAPOR (MM/DAY)', ' PRECIP (MM/DAY)', ' T AIR (.1 C)   ',     1611.   
     2  ' T1 (.1 C)      ', '0T SURF (.1 C)  ',                         1612.   
!    *  '1STAT STB(STRAT)',' STAT STB(TROPO)','0RICH NUM(STRAT)',       1613.   
!    *  ' RICH NUM(TROPO)',' ROSS NUM(STRAT)',' ROSS NUM(TROPO)',       1614.   
     *  '1STAT STB(STRAT)',' STAT STB(TROPO)','0TH2M (DEGREE)',         1613.   
     *  ' RICH NUM(SURF )',' ROSS NUM(STRAT)',' ROSS NUM(TROPO)',       1614.   
     *  ' OCEAN ICE COVER','0SNOW COVER',' TAU L ',' TAU F',            1615.   
     4  ' SURF. FRICTION ', '0DT/DLAT(STRAT) ', ' DT/DLAT(TROPO) '/     1616.   
                         DATA TITLE2/                                   1617.   
     7  ' L(STRAT)(10**5)', ' L(TROP) (10**5)', ' PRECIP HEAT FLX',     1618.   
     O  ' HEAT RUNOFF G1 ', ' HT WTR DIFS -Z1', '0CONDUCTN AT -Z1',     1619.   
     3  '0HEAT RUNOFF G2 ', ' NET HEAT AT Z0 ', '0WTR DIFS AT -Z1',     1620.   
     6  ' NET WTR INTO G1', '0WATER RUNOFF G2', ' NET WTR INTO G2',     1621.   
     9  '0WATER IN G1    ', ' ICE IN G1      ', ' WATER IN G2    ',     1622.   
     2  ' ICE IN G2      ', ' SNOW DEPTH     ', ' WATER RUNOFF G1',     1623.   
     5  ' LW WINDOW BTEMP', ' NET HEAT -Z1-Z2', '0TOT SUP SAT CLD',     1624.   
     8  ' TOT MST CNV CLD', ' TOTAL CLD COVER', ' MC CLD DPTH(MB)',     1625.   
     *  '0SS PRECIP(MM/D)', ' MC PRECIP(MM/D)', ' H2O OF ATM (MM)',     1626.   
     4  '0GAM(K/KM)      ', ' GAMM(K/KM)     ', ' GAMC(K/KM)     ',     1627.   
     *  ' LW INC ON Z0',5*' '/                                          1628.   
      DATA TITLEA/' PLANETARY ALBDO',' PLAN ALB VISUAL',                1629.   
     *  ' PLAN ALB NEARIR', ' SURFACE G ALBDO', ' SURF ALB VISUAL',     1630.   
     *  ' SURF ALB NEARIR', '0ATMO ALB VISUAL', ' ATMO ALB NEARIR',     1631.   
     *  ' ATMO ABS VISUAL', ' ATMO ABS NEARIR'/                         1632.   
      CHARACTER*16 TITCLR(12)
      data TITCLR/
     1  ' INC SW(WT/M**2)', '0SW ABS BELOW P0', ' SW ABS BELOW P1',
     4  ' SW ABS AT Z0   ', ' SW INC ON Z0   ', ' LW INC ON Z0',
     7  '0NET LW AT Z0   ', ' NET LW AT P0   ', ' NET LW AT P1  ',
     *  '0NET RAD AT P0  ', ' NET RAD AT P1  ', ' NET RAD AT Z0  '/
C****                                                                   1633.   
      DATA WTA/1.,0.,1.,0./, WTB/1.,1.,0.,0./, WTC/1.,0.,0.,1./         1634.   
      DATA INDEX/1,2,3,4,5,6,7,8,9,67,  55,10,11,12,13,14,39,40,44,42,  1635.   
     *  41,16,43,15,56,18,17,23,22,21,  35,36,24,25,26,27,28,29,37,38,  1636.   
     *  64,65,66,57,58,59,60,61,62,20,  19,63,45,54,46,47,48,49,50,51,  1637.   
     *  52,53,31,30,32,33,34,3*0/                                       1638.   
      DATA INNUM/2,72,73,6,74,75,76,77,78,79/, INDEN/3*1,5,6*1/         1639.   
C**** IA: 1 CONDENSATION, 2 RADIATION, 3 SURFACE, 4 DIAGA, 0 UNUSED     1640.   
      DATA IA/6*2,  2,2,1,2,2,1,  6*1,  1,1,4,4,3,4,  5*4,1,  6*4,      1641.   
c    *  4,4,4*1,  6*1,  6*1,  2,1,4*2,  1,1,4*4,  2,5*0/                1642.   
     *  1,1,4*1,  6*1,  6*1,  2,1,4*2,  1,1,4*4,  2,5*0/              
      DATA SCALE/6*1.,  6*1.,  4*1.,2*10.,  2*1.,4*10.,  6*100.,        1643.   
     *  6*100.,  6*1.,  6*1.,  6*1.,  2*1.,3*100.,1.,  6*1.,  6*1./     1644.   
c     DATA IFIRST/1/                                                    1645.   
      DATA INQTAB/1,51,3,4,5,6,7,8,48,49,
     &14,15,16,21,23,25,26,88,32,33,
     &34,36,39,77,78,79,11,24,89,52,
     &92,102,12/
c    &101,99,30/
      DATA IFIRST/1/
c     print *,' DIAG1',' IFIRST=',IFIRST
c     print*,'IDACC',IDACC
c     print *,' DT=',DT,' NDYN=',NDYN
c     print *,' INQTAB'
c     print *,INQTAB 
      IF(IFIRST.NE.1) GO TO 100                                         1646.   
      print *,INQTAB 
      IA(32)=1
      IA(33)=1
      IA(34)=1
c     print*,'IDACC',IDACC
c     print *,' IA'
c     do i=1,72
c      print *,i,IA(i)
c     enddo
      IFIRST=0                                                          1647.   
C**** INITIALIZE CERTAIN QUANTITIES  (KD1M LE 69)                       1648.   
      KD1M=67                                                           1649.   
      INC=1+JMM1/24                                                     1650.   
      DTSRCE=DT*NDYN                                                    1651.   
      DTCNDS=DT*NCNDS                                                   1652.   
      JMHALF=JM/2                                                       1653.   
      DO 10 JR=1,46                                                     1654.   
   10 SAREA(JR)=0.                                                      1655.   
      DO 30 J=1,JM                                                      1656.   
      S1(J)=IM                                                          1657.   
      SPLAND(J)=0.                                                      1658.   
      DO 20 I=1,IM                                                      1659.   
      SPLAND(J)=SPLAND(J)+FDATA(I,J,2)                                  1660.   
      JR=J
   20 SAREA(JR)=SAREA(JR)+DXYP(J)                                       1662.   
   30 JLAT(J)=INT(LAT(J)*360./TWOPI+100.5)-100                          1663.   
      S1(1)=1.                                                          1664.   
      S1(JM)=1.                                                         1665.   
      SPLAND(1)=FDATA(1,1,2)                                            1666.   
      SPLAND(JM)=FDATA(1,JM,2)                                          1667.   
      SCALE(9)=1./DTSRCE                                                1668.   
      SCALE(12)=1./DTSRCE                                               1669.   
      SCALE(13)=1./DTSRCE                                               1670.   
      SCALE(14)=1./DTSRCE                                               1671.   
      SCALE(15)=1./DTSRCE                                               1672.   
      SCALE(16)=1./DTSRCE                                               1673.   
      SCALE(19)=SDAY/DTSRCE                                             1674.   
      SCALE(20)=100.*SDAY/(DTCNDS*GRAV)                                 1675.   
      SCALE(24)=1.E3*GRAV*EXPBYK(1000.)                                 1676.   
      SCALE(25)=SCALE(24)                                               1677.   
      SCALE(26)=16.*RGAS                                                1678.   
      SCALE(27)=16.*RGAS                                                1679.   
      SCALE(28)=.25/(2.*OMEGA)                                          1680.   
      SCALE(29)=.25/(2.*OMEGA)                                          1681.   
      if(NOCLM)then
      SCALE(31)=SCALE(31)/IDACC(3)
      endif
      SCALE(35)=.5E2*JMM1/((SIGE(LS1)-SIGE(LM+1)+1.E-12)*180.)          1682.   
      SCALE(36)=.5E2*JMM1/((SIGE(1)-SIGE(LS1))*180.)                    1683.   
c     SCALE(37)=1.E-5*SQRT(RGAS)/(2.*OMEGA)                             1684.   
c     SCALE(37)=0.5
      SCALE(37)=1.0/float(NSURF)
      SCALE(28)=SCALE(37)                                               1685.   
c     38 Total stress TAUL+ROOF STRESS
      SCALE(38)=1.0
c     For T2m
      SCALE(26)=SCALE(23)                                               1685.   
      IA(26)=IA(23)
c     For T2m
c     For dT2m
      SCALE(27)=SCALE(26)                                               1685.   
      IA(27)=IA(26)
c     For dT2m
c     For RIGS
!     SCALE(27)=SCALE(37)                                               1685.   
!     IA(27)=IA(37)
c     For RIGS
      SCALE(39)=1./DTSRCE                                               1686.   
      SCALE(40)=1./DTSRCE                                               1687.   
      SCALE(41)=1./DTSRCE                                               1688.   
      SCALE(42)=1./DTSRCE                                               1689.   
      SCALE(43)=1./DTSRCE                                               1690.   
      SCALE(44)=1./DTSRCE                                               1691.   
      SCALE(45)=SDAY/DTSRCE                                             1692.   
      if(.not.NOCLM)then
c      43 is SWUP at Z0 from CLM
       SCALE(43)=1.0
c      45 is LWUP at Z0 from CLM
       SCALE(45)=1.0
c      16 is vetclm  from CLM
       SCALE(16)=SDAY
c      41 is sevclm  from CLM
       SCALE(41)=SDAY
c      42 is cevclm  from CLM
       SCALE(42)=SDAY
      endif
      SCALE(46)=SDAY/DTSRCE 
      SCALE(47)=SDAY/DTSRCE                                             1693.   
      SCALE(48)=SDAY/DTSRCE                                             1694.   
      SCALE(54)=SDAY/DTSRCE                                             1695.   
      SCALE(56)=1./DTSRCE                                               1696.   
      SCALE(61)=SCALE(20)                                               1697.   
      SCALE(62)=SCALE(20)                                               1698.   
      SCALE(63)=100./GRAV                                               1699.   
      SCALE(64)=1.E3*GRAV/(SIG(1)-SIG(LTM))                             1700.   
      SCALE(65)=1.E3*.0098/(SIGE(1)-SIGE(LTM+1))                        1701.   
      SCALE(66)=1.E3                                                    1702.   
c ============== 020996
c Andrei want me to do this
      if(NSURF.eq.2)then
      SCALE(32)=5.0
      SCALE(33)=5.0
c 0.5*10
      else
      SCALE(32)=10.0
      SCALE(33)=10.0
      endif
c======================
      SCALE(34)=10.
C**** CALCULATE THE DERIVED QUANTITIES                                  1703.   
  100 BYA1=1./(IDACC(1)+1.E-20)                                         1704.   
c     print *,' DIAG1 100'
      A2BYA1=FLOAT(IDACC(2))/FLOAT(IDACC(1))                            1705.   
      A1BYA2=IDACC(1)/(IDACC(2)+1.E-20)                                 1706.   
c     print *,' DTSRCE=',DTSRCE,' A1BYA2=',A1BYA2,' A2BYA1=',A2BYA1
      DO 210 J=1,JM                                                     1722.   
c     print *,' 210 J=',J
      SPOICE(J)=CJ(J,30)*BYA1                                           1723.   
      SPOCEN(J)=S1(J)-SPLAND(J)-SPOICE(J)                               1724.   
c     AJ(J,17)=AJ(J,18)                                                 1725.   
      AJ(J,60)=IDACC(2)*SPOCEN(J)*AJ(J,80)/(AJ(J,58)+1.E-20)            1726.   
      BJ(J,60)=IDACC(2)*SPLAND(J)*BJ(J,80)/(BJ(J,58)+1.E-20)            1727.   
c     if(CJ(J,58).gt.1e-10)then
      CJ(J,60)=IDACC(2)*SPOICE(J)*CJ(J,80)/(CJ(J,58)+1.E-20)            1728.   
c     else
c     CJ(J,60)=0.
c     endif
      DO 210 M=1,3                                                      1729.   
c      print *,'      210 M=',M
      ABCJ(J,4,M)=ABCJ(J,2,M)-ABCJ(J,6,M)                               1730.   
c     ABCJ(J,7,M)=ABCJ(J,70,M)+A2BYA1*ABCJ(J,9,M)/DTSRCE                1731.   
c     ABCJ(J,8,M)=ABCJ(J,71,M)+A2BYA1*ABCJ(J,9,M)/DTSRCE                1732.   
      ABCJ(J,10,M)=ABCJ(J,2,M)+ABCJ(J,7,M)                              1733.   
      ABCJ(J,11,M)=ABCJ(J,3,M)+ABCJ(J,8,M)                              1734.   
      ABCJ(J,12,M)=A1BYA2*ABCJ(J,6,M)*DTSRCE+ABCJ(J,9,M)                1735.   
C  CLEAR SKY
      ABCJCL(J,10,M)=ABCJCL(J,2,M)+ABCJCL(J,8,M)  
      ABCJCL(J,11,M)=ABCJCL(J,3,M)+ABCJCL(J,9,M)  
CCC   ABCJCL(J,7,M)=ABCJCL(J,7,M)/(A1BYA2*DTSRCE)
      ABCJCL(J,12,M)=ABCJCL(J,4,M)+ABCJCL(J,7,M)
C  CLEAR SKY 
      if(NOCLM) then
      ABCJ(J,16,M)=ABCJ(J,41,M)+ABCJ(J,42,M)                            1736.   
      endif
      ABCJ(J,20,M)=ABCJ(J,61,M)+ABCJ(J,62,M)                            1737.   
      if(NOCLM)then
      ABCJ(J,44,M)=ABCJ(J,12,M)+ABCJ(J,13,M)+ABCJ(J,14,M)               1738.   
     *  +ABCJ(J,39,M)-ABCJ(J,40,M)                                      1739.   
      else
      ABCJ(J,44,M)=ABCJ(J,12,M)+ABCJ(J,13,M)+ABCJ(J,14,M)               1738.   
      endif
c     ABCJ(J,46,M)=ABCJ(J,20,M)*SCALE(20)-(ABCJ(J,19,M)+ABCJ(J,45,M)    1740.   
c    *  +ABCJ(J,54,M))*SCALE(19)                                        1741.   
      ABCJ(J,48,M)=ABCJ(J,45,M)-ABCJ(J,47,M)                            1742.   
      ABCJ(J,56,M)=ABCJ(J,15,M)+ABCJ(J,43,M)                            1743.   
  210 CONTINUE                                                          1744.   
      IHOUR0=TOFDY0+.5                                                  1745.   
      IHOUR=TOFDAY+.5                                                   1746.   
      TAUDIF=TAU-TAU0                                                   1747.   
C****                                                                   1748.   
C**** LOOP OVER SURFACE TYPES: GLOBAL, LAND, OCEAN, AND OCEAN ICE       1749.   
C****                                                                   1750.   
      DO 500 M=1,4                                                      1751.   
c     print *,' do 500 M=',M
c     WRITE (6,901) XLABEL                                              1752.   
c     WRITE (6,902) TERAIN(M),IDAY0,IHOUR0,JDATE0,JMNTH0,JYEAR0,TAU0,   1753.   
c    *  IDAY,IHOUR,JDATE,JMONTH,JYEAR,TAU,TAUDIF                        1754.   
c     WRITE (6,903) (JLAT(JM+INC-J),J=INC,JM,INC)                       1755.   
c     WRITE (6,905)                                                     1756.   
      DO 490 K=1,KD1M                                                   1757.   
c     print *,' do 490 K=',K
      N=INDEX(K)                                                        1758.   
      IACC=IDACC(IA(N))                                                 1759.   
      if(K.eq.-16.and.M.eq.2) then
        print *,' DD2'
	print *,' N=',N,' IACC=',IACC,' SCALE',SCALE(N)
	print *,TITLE(n),' M=',M
        print *,' BJ(J,N)=',BJ(J,N)
      endif
      GSUM=0.                                                           1760.   
      GWT=0.                                                            1761.   
      DO 320 JHEMI=1,2                                                  1762.   
c     if(N.eq.60) then
cprint *,' N=',N,' JHEMI=',JHEMI
c     endif
      HSUM=0.                                                           1763.   
      HWT=0.                                                            1764.   
      DO 310 JH=1,JMHALF                                                1765.   
      J=(JHEMI-1)*JMHALF+JH                                             1766.   
      QJ=(AJ(J,N)*WTA(M)+BJ(J,N)*WTB(M)+CJ(J,N)*WTC(M))*SCALE(N)        1767.   
      if(N.eq.-60) then
	print *,AJ(J,N),BJ(J,N),CJ(J,N)
	print *,CJ(J,80),CJ(J,58)
      endif
      WTJ=(SPOCEN (J)*WTA(M)+SPLAND(J)*WTB(M)+SPOICE(J)*WTC(M))*IACC    1768.   
      FLAT(J)=QJ/(WTJ+1.E-20)                                           1769.   
      MLAT(J)=INT(FLAT(J)+10000.5)-10000                                1770.   
      HSUM=HSUM+QJ*DXYP(J)                                              1771.   
  310 HWT=HWT+WTJ*DXYP(J)                                               1772.   
      if(K.eq.-16.and.M.eq.2) then
        print *,' FLAT(J)=',FLAT(J)
      endif
      if(N.eq.-60) then
	print *,' HSUM=',HSUM
	print *,' HWT=',HWT
      endif
      FHEM(JHEMI)=HSUM/(HWT+1.E-20)                                     1773.   
      GSUM=GSUM+HSUM                                                    1774.   
  320 GWT=GWT+HWT                                                       1775.   
      if(N.eq.-60) then
	print *,' GSUM=',GSUM
	print *,' GWT=',GWT
      endif
      FGLOB=GSUM/(GWT+1.E-20)                                           1776.   
         IF(M.EQ.1) CALL KEYD1 (N,FGLOB,FHEM(2))                        1777.   
         DO 323 J=1,JM                                                  1778.   
  323    GBUDG(J,K,M)=FLAT(J)                                           1779.   
         GBUDG(JM+1,K,M)=FHEM(1)                                        1780.   
         GBUDG(JM+2,K,M)=FHEM(2)                                        1781.   
         GBUDG(JM+3,K,M)=FGLOB                                          1782.   
      if(K.eq.-16.and.M.eq.2) then
        print *,' GSUM=',GSUM,' GWT=',GWT
        print *,' FGLOB=',FGLOB
      endif
c     GO TO (350,350,350,350,350,350,  350,350,350,350,350,350,         1783.   
c    *       350,350,350,350,350,350,  340,340,350,350,350,350,         1784.   
c    *       340,350,350,340,340,350,  350,350,350,350,350,350,         1785.   
c    *       350,350,350,350,350,350,  350,350,340,340,340,340,         1786.   
c    *       350,350,350,350,350,340,  350,350,350,350,350,350,         1787.   
c    *       340,340,350,340,340,340,  350,340,350,350,350,350),N       1788.   
c 340 WRITE (6,906) TITLE(N),FGLOB,FHEM(2),FHEM(1),                     1789.   
c    *  (FLAT(JM+INC-J),J=INC,JM,INC)                                   1790.   
c     GO TO 490                                                         1791.   
c 350 WRITE (6,907) TITLE(N),FGLOB,FHEM(2),FHEM(1),                     1792.   
c    *  (MLAT(JM+INC-J),J=INC,JM,INC)                                   1793.   
      if(K.NE.14) go to 491
      DO 570 KCL=1,12
      N=KCL
      GSUM=0.                                                         
      GWT=0.                                                         
      DO 520 JHEMI=1,2                                              
      HSUM=0.                                                     
      HWT=0.                                                    
      DO 510 JH=1,JMHALF                                     
      J=(JHEMI-1)*JMHALF+JH                              
      IACC=NCLR(J)
      IACC=IDACC(IA(1))
      QJ=(AJCLR(J,N)*WTA(M)+BJCLR(J,N)*WTB(M)+CJCLR(J,N)*WTC(M))   
      WTJ=(SPOCEN (J)*WTA(M)+SPLAND(J)*WTB(M)+SPOICE(J)*WTC(M))*IACC
      FLAT(J)=QJ/(WTJ+1.E-20)                                    
      MLAT(J)=INT(FLAT(J)+10000.5)-10000                     
      HSUM=HSUM+QJ*DXYP(J)                                
  510 HWT=HWT+WTJ*DXYP(J)                               
      FHEM(JHEMI)=HSUM/(HWT+1.E-20)                    
      GSUM=GSUM+HSUM                                 
  520 GWT=GWT+HWT                                  
cprint *,' GSUM=',GSUM
cprint *,' GWT=',GWT
      FGLOB=GSUM/(GWT+1.E-20)                     
         DO 453 J=1,JM     
  453    GBUDG(J,KCL+KD1M+10,M)=FLAT(J) 
         GBUDG(JM+1,KCL+KD1M+10,M)=FHEM(1) 
         GBUDG(JM+2,KCL+KD1M+10,M)=FHEM(2)
         GBUDG(JM+3,KCL+KD1M+10,M)=FGLOB 
c     WRITE (6,907) TITCLR(N),FGLOB,FHEM(2),FHEM(1),
c    *  (MLAT(JM+INC-J),J=INC,JM,INC)
  570 CONTINUE
  491 IF(N.NE.1) GO TO 490                                              1794.   
C**** CALCULATE AND PRINT ALBEDOS                                       1795.   
  400 DO 430 KA=1,10                                                    1796.   
cprint *,' KA=',KA
      NN=INNUM(KA)                                                      1797.   
      ND=INDEN(KA)                                                      1798.   
      AMULT=1.                                                          1799.   
      IF(KA.LE.1.OR.KA.EQ.4) AMULT=-1.                                  1800.   
      GSUM=0.                                                           1801.   
      GSUM2=0.                                                          1802.   
      DO 420 JHEMI=1,2                                                  1803.   
      HSUM=0.                                                           1804.   
      HSUM2=0.                                                          1805.   
      DO 410 JH=1,JMHALF                                                1806.   
      J=(JHEMI-1)*JMHALF+JH                                             1807.   
      QNUM=AJ(J,NN)*WTA(M)+BJ(J,NN)*WTB(M)+CJ(J,NN)*WTC(M)              1808.   
      QDEN=AJ(J,ND)*WTA(M)+BJ(J,ND)*WTB(M)+CJ(J,ND)*WTC(M)              1809.   
      FLAT(J)=AMULT*(100.*     QNUM/(QDEN    +1.E-20)-50.)+50.          1810.   
      MLAT(J)=FLAT(J)+.5                                                1811.   
      HSUM=HSUM+QNUM*DXYP(J)                                            1812.   
  410 HSUM2=HSUM2+QDEN*DXYP(J)                                          1813.   
      FHEM(JHEMI)=50.5+AMULT*(100.*HSUM/(HSUM2+1.E-20)-50.)             1814.   
      GSUM=GSUM+HSUM                                                    1815.   
  420 GSUM2=GSUM2+HSUM2                                                 1816.   
cprint *,' GSUM=',GSUM
cprint *,' GSUM2=',GSUM2
      FGLOB=50.5+AMULT*(100.*GSUM/(GSUM2+1.E-20)-50.)                   1817.   
         IF(M.EQ.1.AND.KA.EQ.1) CALL KEYD1A (FGLOB)                     1818.   
         DO 423 J=1,JM                                                  1819.   
  423    GBUDG(J,KA+KD1M,M)=FLAT(J)                                     1820.   
         GBUDG(JM+1,KA+KD1M,M)=FHEM(1)                                  1821.   
         GBUDG(JM+2,KA+KD1M,M)=FHEM(2)                                  1822.   
         GBUDG(JM+3,KA+KD1M,M)=FGLOB                                    1823.   
c     WRITE (6,907) TITLEA(KA),FGLOB,FHEM(2),FHEM(1),                   1824.   
c    *  (MLAT(JM+INC-J),J=INC,JM,INC)                                   1825.   
  430 CONTINUE                                                          1826.   
  490 CONTINUE                                                          1827.   
c     WRITE (6,903) (JLAT(JM+INC-J),J=INC,JM,INC)                       1828.   
c     WRITE (6,905)                                                     1829.   
c     DO 495 LSKIP=1,20                                                 1830.   
c 495 WRITE (6,920)                                                     1831.   
  500 CONTINUE                                                          1832.   
C****                                                                   1833.   
C**** PRODUCE REGIONAL STATISTICS                                       1834.   
C****                                                                   1835.   
      RETURN                                                            1876.   
C****                                                                   1877.   
  901 FORMAT ('1',33A4)                                                 1878.   
  902 FORMAT ('0** BUDGETS',A16,' **    DAY',I5,', HR',I2,' (',I2,A5,   1879.   
     *  I4,')',F8.0,'    TO    DAY',I5,', HR',I2,' (',I2,A5,I4,')',     1880.   
     *  F8.0,'    DIF',F5.0,' HR')                                      1881.   
  903 FORMAT ('0',131('-')/20X,'G     NH    SH  ',24I4)                 1882.   
  904 FORMAT (A16,3I6,2X,24I4)                                          1883.   
  905 FORMAT (1X,131('-'))                                              1884.   
  906 FORMAT (A16,3F6.1,2X,24F4.1)                                      1885.   
  907 FORMAT (A16,3F6.1,2X,24I4)                                        1886.   
  908 FORMAT ('0',17X,'WEST MID- EAST SOU. GRN- MID- NOR. WEST SIBR SOU.1887.   
     * CHNA IND. AUS. NOR. SOU. AFR. AFR. AMZN NOR. MID- NOR. WEST EAST'1888.   
     * /18X,'U.S. U.S. U.S. CNDA LAND EUR. RUSS SIBR PLAT CHNA DSRT DSRT1889.   
     * DSRT SHRA SHRA SAHL RAIN RAIN ATL. ATL. PAC. PAC. PAC. '/1X,     1890.   
     *  131('-'))                                                       1891.   
  909 FORMAT (A16,1X,23I5)                                              1892.   
  910 FORMAT (A16,1X,23F5.1)                                            1893.   
  920 FORMAT (1X)                                                       1894.   
      END                                                               1895.   
      BLOCK DATA a1                                                     2001.   
C****                                                                   2002.   
C**** TITLES FOR SUBROUTINE DIAG2                                       2003.   
C****                                                                   2004.   
      COMMON/D2TTL/TITLE1,TITLE2,TITLE3,                                2005.   
     *  TITLE4,TITLE5,TITLE6,TITLE7,TITLE8,TITLE9,TITLEA,TITLEB,TITLEC  2006.   
     *  ,TITLEN
C    *  ,LINECT,JMHALF,INC,IHOUR0,IHOUR,TAUDIF                          2007.   
                CHARACTER*64 TITLE1(13)/                                2008.   
C****                                                              1-13 2009.   
     1'TEMPERATURE (DEGREES CENTIGRADE)',                               2010.   
     *'HEIGHT (HUNDREDS OF METERS)',                                    2011.   
     3'SPECIFIC HUMIDITY (10**-5 KG H2O/KG AIR)',                       2012.   
     *'RELATIVE HUMIDITY (PERCENT)',                                    2013.   
     *'ZONAL WIND (U COMPONENT) (TENTHS OF METERS/SECOND)',             2014.   
     6'MERIDIONAL WIND (V COMPONENT) (HUNDREDTHS OF METERS/SECOND)',    2015.   
     *'STREAM FUNCTION (10**9 KILOGRAMS/SECOND)',                       2017.   
     8'VERTICAL VELOCITY (10**-5 METERS/SECOND)',                       2018.   
     9'BAROCLINIC EDDY KINETIC ENERGY GEN. (10**-1 WATTS/M**2/SIGMA)',  2019.   
     *'VERTICAL MASS EXCHANGE FROM MOIST CONVECTION (10**9 KG/SECOND)', 2021.   
     *'SOLAR RADIATION HEATING RATE (HUNDREDTHS OF DEGREES KELVIN/DAY)',2023.   
     *'THERMAL RADIATION COOLING RATE (HUNDREDTHS OF DEGREES K/DAY)',   2025.   
     *'TOTAL RADIATION COOLING RATE (10**13 WATTS/UNIT SIGMA)'/         2027.   
                CHARACTER*64 TITLE2(8)/                                 2028.   
C****                                                             14-21 2029.   
     4'HEATING BY LARGE SCALE CONDENSATION (10**13 WATTS/UNIT/SIGMA)',  2030.   
     5'HEATING BY DRY CONVECTION (10**13 WATTS/UNIT SIGMA)',            2032.   
     6'    HEATING BY MOIST CONVECTION (10**13 WATTS/UNIT SIGMA)',      2033.   
     7'STANDING EDDY KINETIC ENERGY (10**4 JOULES/M**2/UNIT SIGMA)',    2035.   
     8'EDDY MERIDIONAL WIND VARIANCE (METER**2/SEC**2) ',               2037.   
     9'TOTAL KINETIC ENERGY (10**4 JOULES/M**2/UNIT SIGMA)',            2038.   
     O'AVAILABLE POTENTIAL ENERGY (10**5 JOULES/M**2/UNIT SIGMA)',      2039.   
     1'POTENTIAL TEMPERATURE (DEGREES KELVIN)'/                         2041.   
                CHARACTER*64 TITLE3(7)/                                 2042.   
C****                                                             22-28 2043.   
     2'NOR. TRANS. OF DRY STAT. ENERGY BY STAND. EDDIES(10**14 W/DSIG)',2044.   
     3'NORTH. TRANS. OF DRY STATIC ENERGY BY EDDIES (10**14 W/DSIG)',   2046.   
     4'TOTAL NORTH. TRANSPORT OF DRY STATIC ENERGY (10**15 W/DSIG)',    2048.   
     5'NORTHWARD TRANSPORT OF LATENT HEAT BY EDDIES (10**13 W/DSIG)',   2050.   
     6'TOTAL NORTHWARD TRANSPORT OF LATENT HEAT (10**14 W/UNIT SIG)',   2052.   
     7'NORTH. TRANSPORT OF STATIC ENERGY BY EDDIES (10**14 W/DSIGMA)',  2054.   
     8'TOTAL NORTHWARD TRANSPORT OF STATIC ENERGY (10**15 W/DSIGMA)'/   2056.   
                CHARACTER*64 TITLE4(5)/                                 2058.   
C****                                                             29-33 2059.   
     9'NORTH. TRANSPORT OF KINETIC ENERGY BY EDDIES (10**12 W/DSIG)',   2060.   
     O'TOTAL NORTHWARD TRANSPORT OF KINETIC ENERGY (10**12 W/DSIG)',    2062.   
     1'NORTH. TRANS. OF ANG. MOMENTUM BY STAND. EDDIES (10**18 J/DSIG)',2064.   
     2'NORTH. TRANS. OF ANG. MOMENTUM BY EDDIES (10**18 J/DSIGMA)',     2066.   
     3'TOTAL NORTHWARD TRANSPORT OF ANG. MOMENTUM (10**19 J/DSIG)'/     2068.   
                CHARACTER*64 TITLE5(6)/                                 2070.   
C****                                                             34-39 2071.   
     4'VERT. TRANS. OF DRY STATIC ENERGY BY EDDIES (10*12 WATTS)',      2072.   
     5'TOTAL LARGE SCALE VERT. TRANS. OF DRY STAT. ENER.(10**14 WATTS)',2074.   
     6'VERTICAL TRANSPORT OF LATENT HEAT BY EDDIES (10**12 WATTS)',     2076.   
     7'TOTAL LARGE SCALE VERT. TRANS. OF LATENT HEAT (10**13 WATTS)',   2078.   
     8'VERTICAL TRANSPORT OF STATIC ENERGY BY EDDIES (10**13 WATTS)',   2080.   
     9'TOTAL LARGE SCALE VERT. TRANS. OF STATIC ENERGY (10**14 W)'/     2082.   
                CHARACTER*64 TITLE6(4)/                                 2084.   
C****                                                             40-43 2085.   
     O'VERTICAL TRANSPORT OF KINETIC ENERGY BY EDDIES (10**11 WATTS)',  2086.   
     1'TOTAL LARGE SCALE VERT. TRANS. OF KINETIC ENERGY (10**11 WATTS)',2088.   
     2'VERT. TRANS. OF ANG. MOMENTUM BY EDDIES (10**16 JOULES)',        2090.   
     3'TOTAL LARGE SCALE VERT. TRANS. OF ANG. MOMENTUM (10**18 JOULES)'/2091.   
                 CHARACTER*64 TITLE7(9)/                                2093.   
C****                                                             44-52 2094.   
     4'CHANGE OF ANG. MOMENTUM BY DRY CONVEC (10**18 J/UNIT SIGMA)',    2095.   
     5'CHANGE OF ANG. MOMENTUM BY MOIST CONV (10**18 J/UNIT SIGMA)',    2097.   
     6'CHANGE OF ANG. MOMENTUM BY DIFFUSION (10**18 J/UNIT SIGMA)',     2099.   
C    7'U WIND AVERAGED OVER I=5-9 (TENTHS OF METERS/SECOND)',           2101.   
     7'NORTHWARD ELIASSEN-PALM FLUX (10**17 JOULES/UNIT SIGMA)',        2102.   
c    8,'V WIND AVERAGED OVER I=5-9 (TENTHS OF METERS/SECOND)',          2103.   
c    9'VERTICAL VELOCITY FOR I=5-9 (10**-5 METERS/SECOND)',             2104.   
     8'SHORTWAVE RADIATION FLUX (W/M**2)',
     9'LONGWAVE RADIATION FLUX (W/M**2)',
C    O'U WIND AVERAGED OVER I=35-3 (TENTHS OF METERS/SECOND)',          2105.   
     O'VERTICAL ELIASSEN-PALM FLUX (10**17 JOULES)',                    2106.   
c    1'V WIND AVERAGED OVER I=35-3 (TENTHS OF METERS/SECOND)',          2107.   
c    2'VERTICAL VELOCITY FOR I=35-3 (10**-5 METERS/SECOND)'/            2108.   
     1'SHORTWAVE RADIATION FLUX CLEAR SKY (W/M**2)',
     2'LONGWAVE RADIATION FLUX CLEAR SKY (W/M**2)'/
                CHARACTER*64 TITLE8(8)/                                 2109.   
C****                                                             53-60 2110.   
     3'POTENTIAL VORTICITY (10**-6 K/(MB-S))',                          2111.   
     4'NORTHWARD TRANSPORT OF Q-G POT. VORTICITY  (10**18 J/DSIG)',     2112.   
     5'P-K BY PRESSURE GRADIENT FORCE  (10**-1 W/M**2/UNIT SIGMA)',     2114.   
     6'Q-G POT. VORTICITY CHANGE OVER LATITUDES (10**-12 1/(SEC-M))',   2116.   
     7'LAGRANGIAN MEANSTREAM FUNCTION (10**9 KILOGRAMS/SECOND)',        2118.   
     8'DYNAMIC CONVERGENCE OF EDDY GEOPOTENTIAL (.1 W/M**2/DSIGMA)',    2119.   
     9'REFRACTION INDEX FOR WAVE NUMBER 1  (10**-8 PER M**2)',          2121.   
     O'REFRACTION INDEX FOR WAVE NUMBER 2  (10**-8 PER M**2)'/          2123.   
                CHARACTER*64 TITLE9(12)/                                2125.   
C****                                                             61-72 2126.   
     1'ZONAL WIND (U COMPONENT) FOR J=11-13 (METERS/SECOND)',           2127.   
     2'MERIDIONAL WIND (V COMPONENT) FOR J=11-13 (METERS/SECOND)',      2128.   
     3'VERTICAL VELOCITY FOR J=11-13 (10**-4 METERS/SEDOND)',           2130.   
     4'TEMPERATURE FOR J=11-13 (DEGREES CENTIGRADE)',                   2131.   
     5'RELATIVE HUMIDITY FOR J=11-13 (PERCENT)',                        2132.   
     6'MOIST CONVECTIVE HEATING FOR J=11-13 (10**13 W/UNIT SIGMA)',     2133.   
     7'TOTAL RADIATIVE COOLING FOR J=11-13 (10**13 W/UNIT SIGMA)',      2135.   
     8' ',                                                              2137.   
     9'VERTICAL VELOCITY AT J=19 (10**-4 METERS/SECOND)',               2138.   
     O'TEMPERATURE AT J=19 (DEGREES CENTIGRADE)',                       2139.   
     1'TOTAL RADIATIVE COOLING AT J=19 (10**13 W/UNIT SIGMA)',          2140.   
     2'ZONAL WIND AT J=19 (METERS/SECOND)'/                             2142.   
                CHARACTER*64 TITLEA(11)/                                2143.   
C****                                                             73-83 2144.   
     3'VERTICAL VELOCITY AT J=21 (10**-4 METERS/SECOND)',               2145.   
     4'TEMPERATURE AT J=21 (DEGREES CENTIGRADE)',                       2146.   
     5'TOTAL RADIATIVE COOLING AT J=21 (10**13 W/UNIT SIGMA)',          2147.   
     6'ZONAL WIND AT J=21 (METERS/SECOND)',                             2149.   
     7'TOTAL CLOUD COVER (.1 * %)',                                     2150.   
     8'SUPER SATURATION CLOUD COVER (.1 * %)',                          2151.   
     9'MOIST CONVECTIVE CLOUD COVER (.1 * %)',                          2152.   
     O'AMPLITUDE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 1 (METERS)',    2153.   
     1'AMPLITUDE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 2 (METERS)',    2155.   
     2'AMPLITUDE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 3 (METERS)',    2157.   
     3'AMPLITUDE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 4 (METERS)'/    2159.   
                CHARACTER*64 TITLEB(9)/                                 2161.   
C****                                                             84-92 2162.   
     4'PHASE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 1 (DEG WEST LONG)', 2163.   
     5'PHASE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 2 (DEG WEST LONG)', 2165.   
     6'PHASE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 3 (DEG WEST LONG)', 2167.   
     7'PHASE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 4 (DEG WEST LONG)', 2169.   
     8'NORTH. TRANS. OF SENSIBLE HEAT BY EDDIES (10**14 W/DSIGMA)',     2171.   
     9'TOTAL NORTHWARD TRANSPORT OF SENSIBLE HEAT (10**15 W/DSIGMA)',   2173.   
     O'VERT. TRANS. OF GEOPOTENTIAL ENERGY BY EDDIES (10**12 WATTS)',   2175.   
     1'TOTAL LARGE SCALE VERT. TRANS. OF GEOPOTEN. ENER. (10**14 W)',   2177.   
     2'SUBGRID SCALE TEMPERATURE VARIANCE (DEGREE**2)  '/               2179.   
                CHARACTER*64 TITLEC(6)/                                 2181.   
C****                                                             93-98 2182.   
     3'DYNAMIC CONVERGENCE OF DRY STATIC ENERGY (10**14 W/INIT SIG)',   2183.   
     4'DIVERGENCE OF THE ELIASSEN-PALM FLUX (10**17 J/UNIT SIGMA)',     2185.   
     5'REFRACTION INDEX FOR WAVE NUMBER 3  (10**-8 PER METWR**2)',      2187.   
     6'REFRACTION INDEX FOR WAVE NUMBER 6  (10**-8 PER METER**2)',      2189.   
     9'REFRACTION INDEX FOR WAVE NUMBER 9  (10**-8 PER METER**2)',      2191.   
     8'CHANGE OF PHASE HEATING BY MOIST CONVECTION (10*13 W/DSIG)'/     2193.   
             CHARACTER*64 TITLEN(5)/
C                                                                99-103
     9' TRANSPORT OF LATENT HEAT BY HOR. DIFF.',
     O' GHANGE OF KIN. ENERGY DUE TO EDDIES',
     1' GHANGE OF SPEC. HUM. BY MOIST CONVECTION ',
     2' TEMPERATURE VARIANCE',
     3' '/
      END                                                               2195.   
      SUBROUTINE DIAG2                                                  2201.   
#include "BD2G04.COM"                                                   2202.   
      COMMON/SPEC2/KMT,KINC,COEK                                        2202.1  
      COMMON U,V,T,P,Q                                                  2203.   
      COMMON/WORK2/SENDEG(72,46),CN(2,37),BYP(46),BYPV(46),BYDXYV(46),  2204.   
     *  AX(JM0,36),ARQX(JM0,3),BX(JM0,36),CX(JM0,36),DX(JM0,36),        2205.   
     *  AMPLTD(JM0,8,4),PHASE(JM0,8,4)                                  2206.   
      COMMON/D2COM/JLAT(46,2),WTJ(46,2,2),                              2207.   
     *  LINECT,JMHALF,INC,IHOUR0,IHOUR,I25,TAUDIF                       2208.   
      common/nqt/NQMAPS
      DIMENSION PL(39),PLE(36),ONES(46),BYDSIG(36),BYDPS(3),BYD2SG(36), 2209.   
     *  BYPKS(3),DACOSV(46),BYDXYP(46),DXYPPO(46),ONESPO(46),DXCOSV(46),2210.   
     *  BYDAPO(46),PMB(7),MW(5),ITIT(5)                                 2211.   
      DATA ITIT/59,60,95,96,97/                                         2212.   
      DATA MW/1,2,3,6,9/                                                2213.   
      DATA PMB/999.9,850.,700.,500.,300.,100.,30./                      2214.   
      DATA ONES/46*1./                                                  2215.   
C**** INITIALIZE CERTAIN QUANTITIES                                     2216.   
      NQMAPS=0
      INC=1+JMM1/24                                                     2217.   
      JMHALF=JM/2                                                       2218.   
      BYIM=1./FIM                                                       2219.   
      BY100G=.01/GRAV                                                   2220.   
      SHA=RGAS/KAPA                                                     2221.   
      DTCNDS=NCNDS*DT                                                   2222.   
      P1000K=EXPBYK(1000.)                                              2223.   
      KM=0                                                              2224.   
      DO 5 K=1,7                                                        2225.   
      IF (PTOP.GT.PMB(K)) GO TO 6                                       2226.   
    5 KM=KM+1                                                           2227.   
    6 ELOFIM=.5*TWOPI-TWOPI/FIM                                         2228.   
      DO 20 L=1,LM                                                      2229.   
      LUP=L+1                                                           2230.   
      LDN=L-1                                                           2231.   
      IF(L.EQ.LM) LUP=LM                                                2232.   
      IF(L.EQ.1) LDN=1                                                  2233.   
      BYD2SG(L)=1./(SIG(LUP)-SIG(LDN))                                  2234.   
      BYDSIG(L)=1./DSIG(L)                                              2235.   
   20 PL(L)=SIG(L)*(PSF-PTOP)+PTOP                                      2236.   
      PL(LM+1)=.75*PTOP                                                 2237.   
      PL(LM+2)=.35*PTOP                                                 2238.   
      PL(LM+3)=.1*PTOP                                                  2239.   
      BYDPS(1)=1./(.5*PTOP)                                             2240.   
      BYDPS(2)=1./(.3*PTOP)                                             2241.   
      BYDPS(3)=1./(.2*PTOP)                                             2242.   
      BYPKS(1)=1./(.75*PTOP)**KAPA                                      2243.   
      BYPKS(2)=1./(.35*PTOP)**KAPA                                      2244.   
      BYPKS(3)=1./(.1*PTOP)**KAPA                                       2245.   
      DO 30 L=1,LM                                                      2246.   
   30 PLE(L)=SIGE(L+1)*(PSF-PTOP)+PTOP                                  2247.   
      DO 40 J=1,JM                                                      2248.   
      DXYPPO(J)=DXYP(J)                                                 2249.   
      BYDXYP(J)=1./DXYP(J)                                              2250.   
      BYDAPO(J)=BYDXYP(J)                                               2251.   
      ONESPO(J)=1.                                                      2252.   
      JLAT(J,1)=INT(.5+(J-1.0)*180./JMM1)-90                            2253.   
      JLAT(J,2)=INT(.5+(J-1.5)*180./JMM1)-90                            2254.   
      WTJ(J,1,1)=1.                                                     2255.   
   40 WTJ(J,2,1)=2.*FIM*DXYP(J)/AREAG                                   2256.   
      DXYPPO(JM)=DXYP(JM)*FIM                                           2257.   
      DXYPPO(1)=DXYP(1)*FIM                                             2258.   
      BYDAPO(1)=BYDAPO(1)*FIM                                           2259.   
      BYDAPO(JM)=BYDAPO(JM)*FIM                                         2260.   
      ONESPO(1)=FIM                                                     2261.   
      ONESPO(JM)=FIM                                                    2262.   
      DO 50 J=2,JM                                                      2263.   
      DXCOSV(J)=DXV(J)*COSV(J)                                          2264.   
      DACOSV(J)=DXYV(J)*COSV(J)                                         2265.   
      BYDXYV(J)=1./DXYV(J)                                              2266.   
      WTJ(J,1,2)=1.                                                     2267.   
   50 WTJ(J,2,2)=2.*FIM*DXYV(J)/AREAG                                   2268.   
      WTJ(JMHALF+1,1,2)=.5                                              2269.   
      WTJ(JMHALF+1,2,2)=WTJ(JMHALF+1,2,2)/2.                            2270.   
      IHOUR0=TOFDY0+.5                                                  2271.   
      IHOUR=TOFDAY+.5                                                   2272.   
      TAUDIF=TAU-TAU0                                                   2273.   
      LINECT=65                                                         2274.   
c     WRITE (6,901)                                                     2275.   
      BYIAC3=1./(IDACC(1)+1.E-20)                                       2276.   
      BYIARD=1./(IDACC(2)+1.E-20)                                       2277.   
      BYIADA=1./(IDACC(4)+1.E-20)                                       2278.   
      BYIMDA=BYIADA*BYIM                                                2279.   
      FIMDA=IDACC(4)*FIM                                                2280.   
      SCLRH=100.*IDACC(4)/IDACC(1)
      SCLTV=IDACC(4)/IDACC(1)
      DO 120 J=1,JM                                                     2281.   
      BYP(J)=1./(APJ(J,1)+1.E-20)                                       2282.   
  120 BYPV(J)=1./(APJ(J,2)+1.E-20)                                      2283.   
C****                                                                   2284.   
C**** PROGNOSTIC QUANTITIES                                             2285.   
C****                                                                   2286.   
C**** TEMPERATURE, HEIGHT, SPECIFIC HUMIDITY, AND RELATIVE HUMIDITY     2287.   
      CALL JLMAPS (1,PL,AJL,ONES,BYP,ONES,LM,2,1,                       2288.   
     *  ASJL,BYIMDA,ONESPO,ONES)                                        2289.   
      SCALES=BYIMDA*BY100G                                              2290.   
      CALL JLMAPS (2,PL,AJL(1,1,2),BY100G,BYP,ONES,LM,2,1,              2291.   
     *  ASJL(1,1,2),SCALES,ONESPO,ONES)                                 2292.   
      CALL JLMAP (3,PL,AJL(1,1,3),1.E5,BYP,ONES,LM,2,1)                 2293.   
      CALL JLMAP (100,PL,AJL(1,1,56),1.E5,BYP,ONES,LM,2,1)
c     change of KIN. EN. due to eddy
c     CALL JLMAP (100,PL,AJL(1,1,26),1.E1,BYP,ONES,LM,2,1)
c     change of KIN. EN. due to eddy
      CALL JLMAP (101,PL,AJL(1,1,57),1.E5,BYP,ONES,LM,2,1)
      CALL JLMAP (4,PL,AJL(1,1,18),100.,BYP,ONES,LM,2,1)                2294.   
c     print *,' DIAGA2 AJL=',AJL(12,1,58),AJL(12,1,58)/IDACC(1)
c     print *,BYP(12),AJL(12,1,58)*BYP(12),APJ(12,1),APJ(12,1)/IDACC(1)
c     print *,IDACC(4),APJ(12,1)/IDACC(4)
C===  CALL JLMAP (102,PL,AJL(1,1,58),SCLRH,BYP,ONES,LM,2,1)
      CALL JLMAP (102,PL,AJL(1,1,59),SCLTV,BYP,ONES,LM,2,1)
C**** U WIND, V WIND, AND STREAM FUNCTION                               2295.   
      CALL JLMAP (5,PL,AJL(1,1,4),1.E1,BYPV,ONES,LM,2,2)                2296.   
      CALL JLMAP (6,PL,AJL(1,1,5),1.E2,BYPV,ONES,LM,2,2)                2297.   
      DO 220 J=2,JM                                                     2298.   
      AX(J,1)=AJL(J,1,5)*DSIG(1)                                        2299.   
      DO 220 L=2,LM                                                     2301.   
  220 AX(J,L)=AX(J,L-1)+AJL(J,L,5)*DSIG(L)                              2303.   
      SCALE=25.E-9*BYIADA/GRAV                                          2304.   
      CALL JLMAP (7,PLE,AX,SCALE,DXV,ONES,LM,2,2)                       2305.   
C**** VERTICAL VELOCITY AND MASS FLUX MOIST CONVECTION                  2307.   
      SCALE=-1.E5*BYIADA*RGAS/(FIM*GRAV)                                2308.   
      CALL JLMAP (8,PLE,AJL(1,1,6),SCALE,BYDAPO,ONES,LMM1,2,1)          2309.   
      SCALE=100.E-9*BYIAC3/(GRAV*DTCNDS)                                2310.   
      CALL JLMAP (10,PLE,AJL(1,1,8),SCALE,DXYPPO,ONES,LMM1,1,1)         2311.   
C****                                                                   2312.   
C**** RADIATION, CONDENSATION AND CONVECTION                            2313.   
C****                                                                   2314.   
C**** SOLAR AND THERMAL RADIATION HEATING                               2315.   
      SCALE=100.E-2*GRAV*SDAY*IDACC(4)*BYIARD/SHA                       2316.   
C     SCALE for 100e-2 degree/day
      SCALES=100.E-2*GRAV*SDAY*BYIM*BYIARD/SHA                          2317.   
      CALL JLMAPS (11,PL,AJL(1,1,9),SCALE,BYP,BYDSIG,LM,2,1,            2318.   
     *  ASJL(1,1,3),SCALES,ONESPO,BYDPS)                                2319.   
      SCALES=-SCALES                                                    2320.   
      SCALE=-SCALE                                                      2321.   
      CALL JLMAPS (12,PL,AJL(1,1,10),SCALE,BYP,BYDSIG,LM,2,1,           2322.   
     *  ASJL(1,1,4),SCALES,ONESPO,BYDPS)                                2323.   
      DO 250 J=1,JM                                                     2324.   
      DO 240 LS=1,3                                                     2325.   
  240 ARQX(J,LS)=ASJL(J,LS,3)+ASJL(J,LS,4)                              2326.   
      DO 250 L=1,LM                                                     2327.   
  250 AX(J,L)=AJL(J,L,9)+AJL(J,L,10)                                    2328.   
      SCALE=-1.E-13*BYIARD                                              2329.   
      SCALES=SCALE*(PSF-PTOP)                                           2330.   
      CALL JLMAPS (13,PL,AX,SCALE,DXYPPO,BYDSIG,LM,1,1,                 2331.   
     *  ARQX,SCALES,DXYPPO,BYDPS)                                       2332.   
C**** SOLAR AND THERMAL RADIATION FLUXES

      SCALE=BYIARD                                      
      CALL JLMAP (48,PL,AJL(1,1,42),SCALE,ONES,ONES,LM,2,1)      
      CALL JLMAP (49,PL,AJL(1,1,43),SCALE,ONES,ONES,LM,2,1)   
      CALL JLMAP (51,PL,AJL(1,1,45),SCALE,ONES,ONES,LM,2,1)        
      CALL JLMAP (52,PL,AJL(1,1,46),SCALE,ONES,ONES,LM,2,1) 

C**** TOTAL, SUPER SATURATION, AND CONVECTIVE CLOUD COVER               2333.   
      SCALE=1000.*BYIARD*BYIM                                           2334.   
      CALL JLMAP (77,PL,AJL(1,1,19),SCALE,ONESPO,ONES,LM,2,1)           2335.   
      CALL JLMAP (78,PL,AJL(1,1,28),SCALE,ONESPO,ONES,LM,2,1)           2336.   
      CALL JLMAP (79,PL,AJL(1,1,29),SCALE,ONESPO,ONES,LM,2,1)           2337.   
C**** SUBGRID SCALE TEMPERATURE DEVIATION                               2338.   
      SCALE=1.00*BYIADA                                                 2339.   
      CALL JLMAP (92,PL,AJL(1,1,54),SCALE,ONES,ONES,LM,2,1)             2340.   
C**** HEATING BY LARGE SCALE CONDENSATION AND DRY AND MOIST CONVECTION  2341.   
C     and vert. diff.
C     SCALE for 10**13W/( UNIT SIGMA)
      SCALE=100.E-13*SHA*BYIAC3/(GRAV*DTCNDS)                           2342.   
c     CALL JLMAP (14,PL,AJL(1,1,11),SCALE,DXYPPO,ONES,LM,1,1)           2343.   
c     CALL JLMAP (15,PL,AJL(1,1,12),SCALE,DXYPPO,ONES,LM,1,1)           2344.   
c     CALL JLMAP (16,PL,AJL(1,1,13),SCALE,DXYPPO,ONES,LM,1,1)           2345.   
cc    CALL JLMAP (99,PL,AJL(1,1,55),SCALE,DXYPPO,ONES,LM,1,1)
c     SCALE for 0.01K/DAY
      SCALE=100.*SDAY*IDACC(4)*BYIAC3/DTCNDS 
      CALL JLMAP (14,PL,AJL(1,1,11),SCALE,BYP,ONES,LM,2,1)    
      CALL JLMAP (15,PL,AJL(1,1,12),SCALE,BYP,ONES,LM,2,1)   
      CALL JLMAP (16,PL,AJL(1,1,13),SCALE,BYP,ONES,LM,2,1)  
c     CALL JLMAP (99,PL,AJL(1,1,55),SCALE,BYP,ONES,LM,2,1)
C****                                                                   2347.   
C**** CALCULATIONS FOR STANDING EDDIES                                  2348.   
C****                                                                   2349.   
      IF(SKIPSE.EQ.1.) GO TO 282                                        2350.   
      DO 255 J=2,JM                                                     2351.   
      DO 255 L=1,LM                                                     2352.   
      AX(J,L)=0.                                                        2353.   
      BX(J,L)=0.                                                        2354.   
  255 CX(J,L)=0.                                                        2355.   
      DO 280 J=2,JM                                                     2356.   
      DO 260 I=1,IM                                                     2357.   
  260 SENDEG(I,J)=0.                                                    2358.   
      DO 280 L=1,LM                                                     2359.   
      PU4TI=0.                                                          2360.   
      PV4TI=0.                                                          2361.   
      DE16TI=0.                                                         2362.   
      SKE4I=0.                                                          2363.   
      SNDEGI=0.                                                         2364.   
      SNAM4I=0.                                                         2365.   
      DO 270 I=1,IM                                                     2366.   
      PU4TI=PU4TI+AIJL(I,J,L,1)                                         2367.   
      PV4TI=PV4TI+AIJL(I,J,L,2)                                         2368.   
      DE16TI=DE16TI+AIJL(I,J,L,3)                                       2369.   
      SKE4I=SKE4I+(AIJL(I,J,L,1)*AIJL(I,J,L,1)                          2370.   
     *            +AIJL(I,J,L,2)*AIJL(I,J,L,2))/AIJ(I,J,8)              2371.   
      SNDEGI=SNDEGI+(AIJL(I,J,L,3)*AIJL(I,J,L,2)/AIJ(I,J,8))            2372.   
      SENDEG(I,J)=SENDEG(I,J)                                           2373.   
     *  +DSIG(L)*(AIJL(I,J,L,3)*AIJL(I,J,L,2)/AIJ(I,J,8))               2374.   
      SNAM4I=SNAM4I+AIJL(I,J,L,1)*AIJL(I,J,L,2)/AIJ(I,J,8)              2375.   
  270 CONTINUE                                                          2376.   
      AX(J,L)=SKE4I-(PU4TI*PU4TI+PV4TI*PV4TI)/APJ(J,2)                  2377.   
      SZNDEG=DE16TI*PV4TI/APJ(J,2)                                      2378.   
      BX(J,L)=SNDEGI-SZNDEG                                             2379.   
      SZNDEG=SZNDEG*DSIG(L)/FIM                                         2380.   
      DO 275 I=1,IM                                                     2381.   
  275 SENDEG(I,J)=SENDEG(I,J)-SZNDEG                                    2382.   
  280 CX(J,L)=SNAM4I-PU4TI*PV4TI/APJ(J,2)                               2383.   
C****                                                                   2384.   
C**** ENERGY                                                            2385.   
C****                                                                   2386.   
C**** STANDING EDDY, EDDY, AND TOTAL KINETIC ENERGY                     2387.   
  282 SCALE=12.5E-4*BYIMDA/GRAV                                         2388.   
  285 DO 288 L=1,LM                                                     2391.   
      DO 288 J=2,JM                                                     2392.   
 288  AX(J,L)=AJL(J,L,14)                                               2393.   
c     CALL JLMAP (18,PL,AX,BYIMDA,ONES,ONES,LM,2,2)                     2394.   
C
      DO L=1,LM 
       DO  J=2,JM   
         AX(J,L)=AJL(J,L,15) 
       enddo
      enddo
c        CALL JLMAP (19,PL,AX,BYIMDA,ONES,ONES,LM,2,2)     
         CALL JLMAP (19,PL,AX,SCALE,ONES,ONES,LM,2,2)     
C**** AVAILABLE POTENTIAL ENERGY, POTENTIAL TEMPERATURE AND VORTICITY   2396.   
      SCALE=50.E-5*RGAS*BYIMDA/GRAV                                     2397.   
      CALL JLMAP (20,PL,AJL(1,1,16),SCALE,ONES,ONES,LM,2,1)             2398.   
      DO 298 LR=1,3                                                     2399.   
      DO 298 J=1,JM                                                     2400.   
  298 ARQX(J,LR)=ASJL(J,LR,1)*BYIMDA*ONESPO(J)+273.16                   2401.   
      CALL JLMAPS (21,PL,AJL(1,1,17),P1000K,BYP,ONES,LM,2,1,            2402.   
     *  ARQX,P1000K,ONES,BYPKS)                                         2403.   
C****                                                                   2406.   
C**** NORTHWARD TRANSPORTS                                              2407.   
C****                                                                   2408.   
C**** NORTHWARD TRANSPORT OF SENSIBLE HEAT BY EDDIES                    2409.   
      SCALE=6.25E-14*SHA*BYIADA/GRAV                                    2410.   
      DO 302 L=1,LM                                                     2411.   
      DO 302 J=2,JM                                                     2412.   
  302 AX(J,L)=AJL(J,L,21)-AJL(J,L,20)                                   2413.   
      CALL JLMAP (88,PL,AX,SCALE,DXV,ONES,LM,1,2)                       2414.   
c     Total NORTHWARD TRANSPORT OF SENSIBLE HEAT
	do L=1,LM
	do J=2,Jm
	 AX(J,L)=AJL(J,L,21)
	end do
	end do
	SCALE=6.25E-15*SHA*BYIADA/GRAV 
	 CALL JLMAP (89,PL,AX,SCALE,DXV,ONES,LM,1,2) 
c
C**** NORTHWARD TRANSPORT OF DRY STATIC ENERGY BY STANDING EDDIES,      2415.   
C****   EDDIES, AND TOTAL                                               2416.   
      SCALE=6.25E-14*BYIADA/GRAV                                        2417.   
      IF(SKIPSE.EQ.1.) GO TO 320                                        2418.   
      CALL JLMAP (22,PL,BX,SCALE,DXV,ONES,LM,1,2)                       2419.   
  320 DO 330 L=1,LM                                                     2420.   
      DO 330 J=2,JM                                                     2421.   
      AX(J,L)=SHA*(AJL(J,L,21)-AJL(J,L,20))+(AJL(J,L,23)-AJL(J,L,22))   2422.   
  330 BX(J,L)=SHA*AJL(J,L,21)+AJL(J,L,23)                               2423.   
      CALL JLMAP (23,PL,AX,SCALE,DXV,ONES,LM,1,2)                       2424.   
      SCALE=SCALE*.1                                                    2425.   
      CALL JLMAP (24,PL,BX,SCALE,DXV,ONES,LM,1,2)                       2426.   
C**** NORTHWARD TRANSPORT OF LATENT HEAT BY EDDIES AND TOTAL            2427.   
      DO 336 L=1,LM                                                     2428.   
      DO 336 J=2,JM                                                     2429.   
      DX(J,L)=AJL(J,L,25)-AJL(J,L,24)                                   2430.   
  336 AX(J,L)=AX(J,L)+LHE*DX(J,L)                                       2431.   
      SCALE=6.25E-13*LHE*BYIADA/GRAV                                    2433.   
      CALL JLMAP (25,PL,DX,SCALE,DXV,ONES,LM,1,2)                       2434.   
      SCALE=SCALE*.1                                                    2435.   
      CALL JLMAP (26,PL,AJL(1,1,25),SCALE,DXV,ONES,LM,1,2)              2436.   
C     NORTHWARD TRANSPORT OF LATENT HEAT BY HOR. DIFF.
c
      CALL JLMAP (99,PL,AJL(1,1,55),SCALE,DXV,ONES,LM,1,2)
c
C**** NORTHWARD TRANSPORT OF STATIC ENERGY BY EDDIES AND TOTAL          2437.   
      DO 340 L=1,LM                                                     2437.11 
      DO 340 J=2,JM                                                     2437.12 
  340 DX(J,L)=BX(J,L)+LHE*AJL(J,L,25)                                   2437.13 
      SCALE=6.25E-14*BYIADA/GRAV                                        2438.   
      CALL JLMAP (27,PL,AX,SCALE,DXV,ONES,LM,1,2)                       2439.   
      SCALE=SCALE*.1                                                    2440.   
      CALL JLMAP (28,PL,DX,SCALE,DXV,ONES,LM,1,2)                       2441.   
C**** NORTHWARD TRANSPORT OF KINETIC ENERGY                             2442.   
      SCALE=6.25E-12*BYIADA/GRAV
      CALL JLMAP (30,PL,AJL(1,1,27),SCALE,DXV,ONES,LM,1,2)
C**** NOR. TRANS. OF ANG. MOMENTUM BY STANDING EDDIES, EDDIES AND TOTAL 2445.   
      SCALE=25.E-18*RADIUS*BYIADA/GRAV                                  2446.   
      IF(SKIPSE.EQ.1.) GO TO 350                                        2447.   
      CALL JLMAP (31,PL,CX,SCALE,DXCOSV,ONES,LM,1,2)                    2448.   
  350 DO 360 L=1,LM                                                     2449.   
      DO 360 J=2,JM                                                     2450.   
      CX(J,L)=AJL(J,L,49)-AJL(J,L,48)                                   2451.   
  360 DX(J,L)=AJL(J,L,49)+RADIUS*OMEGA*COSV(J)*AJL(J,L,5)               2452.   
      CALL JLMAP (32,PL,CX,SCALE,DXCOSV,ONES,LM,1,2)                    2453.   
      SCALE=.1*SCALE                                                    2454.   
      CALL JLMAP (33,PL,DX,SCALE,DXCOSV,ONES,LM,1,2)                    2455.   
C**** NOR. TRANSPORT OF QUASI-GEOSTROPHIC POT. VORTICITY BY EDDIES      2456.   
C****                                                                   2501.   
C**** VERTICAL TRANSPORTS                                               2502.   
C****                                                                   2503.   
C**** VERTICAL TRANSPORT OF GEOPOTENTIAL ENERGY BY EDDIES               2504.   
C**** VERTICAL TRANSPORT OF DRY STATIC ENERGY BY EDDIES AND TOTAL       2507.   
      DO 390 L=1,LMM1                                                   2508.   
      DO 390 J=1,JM                                                     2509.   
      AX(J,L)=AJL(J,L,31)-AJL(J,L,30)                                   2510.   
  390 BX(J,L)=AJL(J,L,33)-AJL(J,L,32)                                   2511.   
      SCALE=-50.E-12*BYIADA/GRAV                                        2512.   
      CALL JLMAP (34,PLE,AX,SCALE,ONESPO,ONES,LMM1,1,1)                 2513.   
      SCALE=SCALE*.01                                                   2514.   
      CALL JLMAP (35,PLE,AJL(1,1,31),SCALE,ONESPO,ONES,LMM1,1,1)        2515.   
C**** VERTICAL TRANSPORT OF LATENT HEAT BY EDDIES AND TOTAL             2516.   
      SCALE=-200.E-12*LHE*BYIADA/GRAV                                   2517.   
      CALL JLMAP (36,PLE,BX,SCALE,ONESPO,ONES,LMM1,1,1)                 2518.   
      SCALE=SCALE*.1                                                    2519.   
      CALL JLMAP (37,PLE,AJL(1,1,33),SCALE,ONESPO,ONES,LMM1,1,1)        2520.   
C**** VERTICAL TRANSPORT OF STATIC ENERGY BY EDDIES AND TOTAL           2521.   
      DO 420 L=1,LMM1                                                   2522.   
      DO 420 J=1,JM                                                     2523.   
      AX(J,L)=AX(J,L)+4.*LHE*BX(J,L)                                    2524.   
  420 BX(J,L)=AJL(J,L,31)+4.*LHE*AJL(J,L,33)                            2525.   
      SCALE=-50.E-13*BYIADA/GRAV                                        2526.   
      CALL JLMAP (38,PLE,AX,SCALE,ONESPO,ONES,LMM1,1,1)                 2527.   
      SCALE=SCALE*.1                                                    2528.   
      CALL JLMAP (39,PLE,BX,SCALE,ONESPO,ONES,LMM1,1,1)                 2529.   
C**** VERTICAL TRANSPORT OF KINETIC ENERGY                              2530.   
C**** VERTICAL TRANSPORT OF ANGULAR MOMENTUM BY LARGE SCALE MOTIONS     2533.   
      SCALE=-12.5E-16*RADIUS*BYIADA/GRAV                                2534.   
      CALL JLMAP (42,PLE,AJL(1,1,36),SCALE,COSV,ONES,LMM1,1,2)          2535.   
      SCALE=1.E-2*SCALE                                                 2536.   
      CALL JLMAP (43,PLE,AJL(1,1,37),SCALE,COSV,ONES,LMM1,1,2)          2537.   
C**** VERTICAL TRANSPORT OF ANGULAR MOMENTUM BY SMALL SCALE MOTIONS     2538.   
      SCALE=100.E-18*RADIUS*BYIAC3/(GRAV*DTCNDS)                        2539.   
      CALL JLMAP (44,PL,AJL(1,1,38),SCALE,DACOSV,ONES,LM,1,2)           2540.   
      CALL JLMAP (45,PL,AJL(1,1,39),SCALE,DACOSV,ONES,LM,1,2)           2541.   
C     CALL JLMAP (46,PL,AJL(1,1,40),SCALE,DACOSV,BYDSIG,LM,1,2)         2542.   
         IF(JM.NE.24) GO TO 425                                         2543.   
      IF(IM.EQ.1) GO TO 425                                             2543.5  
C****                                                                   2544.   
C**** MERIDIONAL LUNES                                                  2545.   
C****                                                                   2546.   
C**** U, V AND W VELOCITY FOR I=5-9                                     2547.   
      SCALE=.2E+1*BYIADA                                                2548.   
C     CALL JLMAP (47,PL,AJL(1,1,41),SCALE,ONES,ONES,LM,2,2)             2549.   
c     CALL JLMAP (48,PL,AJL(1,1,42),SCALE,ONES,ONES,LM,2,2)             2550.   
c     SCALE2=-1.E5*BYIADA*RGAS/(5.*GRAV)                                2551.   
c     CALL JLMAP (49,PLE,AJL(1,1,43),SCALE2,BYDXYP,ONES,LMM1,2,1)       2552.   
C**** U, V AND W VELOCITY FOR I=35-3                                    2553.   
C     CALL JLMAP (50,PL,AJL(1,1,44),SCALE,ONES,ONES,LM,2,2)             2554.   
c     CALL JLMAP (51,PL,AJL(1,1,45),SCALE,ONES,ONES,LM,2,2)             2555.   
c     CALL JLMAP (52,PLE,AJL(1,1,46),SCALE2,BYDXYP,ONES,LMM1,2,1)       2556.   
C****                                                                   2557.   
C**** LATITUDINAL ZONES                                                 2558.   
C****                                                                   2559.   
C**** U, V AND W VELOCITY FOR J=11-13 VS. LONGITUDE                     2560.   
      SCALE=BYIADA/3.                                                   2561.   
      CALL ILMAP (61,PL,AIL(1,1,1),SCALE,ONES,LM,2,2)                   2562.   
C     CALL ILMAP (62,PL,AIL(1,1,2),SCALE,ONES,LM,2,2)                   2563.   
      SCALE =-1.E4*BYIADA*RGAS/(GRAV*(DXYP(11)+DXYP(12)+DXYP(13)))      2564.   
      CALL ILMAP (63,PLE,AIL(1,1,3),SCALE,ONES,LMM1,2,1)                2565.   
C**** TEMPERATURE, RELATIVE HUMIDITY, MOIS CONVECTIVE HEATING, AND      2566.   
C****   RADIATIVE COOLING FOR J=11-13 VS. LONGITUDE                     2567.   
      SCALE=BYIADA/3.                                                   2568.   
      CALL ILMAP (64,PL,AIL(1,1,4),SCALE,ONES,LM,2,1)                   2569.   
      SCALE=1.E2*SCALE                                                  2570.   
      CALL ILMAP (65,PL,AIL(1,1,5),SCALE,ONES,LM,2,1)                   2571.   
      SCALE=100.E-13*SHA*BYIAC3/(GRAV*DTCNDS)                           2572.   
      CALL ILMAP (66,PL,AIL(1,1,6),SCALE,ONES,LM,1,1)                   2573.   
      SCALE=-1.E-13*BYIARD                                              2574.   
      CALL ILMAP (67,PL,AIL(1,1,7),SCALE,BYDSIG,LM,1,1)                 2575.   
C**** AT J=19: W VELOCITY, TEMPERATURE, RADIATION, AND U VELOCITY       2576.   
C     SCALE =-1.E4*BYIADA*RGAS/(GRAV* DXYP(19))                         2577.   
C     CALL ILMAP (69,PLE,AIL(1,1,9),SCALE,ONES,LMM1,2,1)                2578.   
      CALL ILMAP (70,PL,AIL(1,1,10),BYIADA,ONES,LM,2,1)                 2579.   
C     SCALE=-1.E-13*BYIARD                                              2580.   
C     CALL ILMAP (71,PL,AIL(1,1,11),SCALE,BYDSIG,LM,1,1)                2581.   
      SCALE=BYIADA/2.                                                   2582.   
      CALL ILMAP (72,PL,AIL(1,1,12),SCALE,ONES,LM,2,2)                  2583.   
C**** AT J=21: W VELOCITY, TEMPERATURE, RADIATION, AND U VELOCITY       2584.   
C     SCALE =-1.E4*BYIADA*RGAS/(GRAV* DXYP(21))                         2585.   
C     CALL ILMAP (73,PLE,AIL(1,1,13),SCALE,ONES,LMM1,2,1)               2586.   
C     CALL ILMAP (74,PL,AIL(1,1,14),BYIADA,ONES,LM,2,1)                 2587.   
C     SCALE=-1.E-13*BYIARD                                              2588.   
C     CALL ILMAP (75,PL,AIL(1,1,15),SCALE,BYDSIG,LM,1,1)                2589.   
C     SCALE=BYIADA/2.                                                   2590.   
C     CALL ILMAP (76,PL,AIL(1,1,16),SCALE,ONES,LM,2,2)                  2591.   
 425  CONTINUE                                                          2591.5  
C****                                                                   2592.   
C**** ELIASSEN-PALM FLUX : NORTHWARD, VERTICAL, DIVERGENCE              2593.   
C****                                                                   2594.   
      SCALE=100.E-17*BYIADA*RADIUS/GRAV                                 2595.   
      DXCVS=DXCOSV(2)                                                   2599.   
      DO 540 J=2,JMM1                                                   2600.   
      BDN=0.                                                            2601.   
      DXCVN=DXCOSV(J+1)                                                 2603.   
      DO 530 L=1,LM                                                     2604.   
      BUP=AJL(J,L,44)*COSP(J)                                           2605.   
      AX(J,L)=AJL(J+1,L,41)*DXCVN-AJL(J,L,41)*DXCVS+                    2606.   
     *   .125*(BUP-BDN)/DSIG(L)                                         2607.   
  530 BDN=BUP                                                           2608.   
  540 DXCVS=DXCVN                                                       2609.   
      DO 550 L=1,LM                                                     2610.   
      AX(1,L)=0.                                                        2611.   
  550 AX(JM,L)=0.                                                       2612.   
      CALL JLMAP(94,PL,AX,SCALE,ONES,ONES,LM,1,1)                       2613.   
C****                                                                   2614.   
C****                                                                   2615.   
C**** D/DY OF QUASI-GEOSTROPHIC POTENTIAL VORTICITY                     2616.   
C****                                                                   2617.   
      IF(KMT.EQ.1) RETURN                                               2617.5  
      AMA=2.*OMEGA/RADIUS                                               2618.   
      PTOPI=PTOP*FIMDA                                                  2619.   
      DO 580 L=1,LM                                                     2620.   
      LUP=L+1                                                           2621.   
      IF (L.EQ.LM) LUP=LM                                               2622.   
      LDN=L-1                                                           2623.   
      IF (L.EQ.1) LDN=1                                                 2624.   
      DO 570 J=2,JMM1                                                   2625.   
      AX(J,L)=F(J)*AJL(J,L,17)/(DXYP(J)*(AJL(J,L,1)*BYP(J)+273.16) *    2626.   
     *  (AJL(J,LUP,17)-AJL(J,LDN,17))+1.E-20)                           2627.   
      BX(J,L)=((AJL(J,LUP,1)*BYP(J)+273.16)/(APJ(J,1)*SIG(LUP)+PTOPI)-  2628.   
     *  (AJL(J,LDN,1)*BYP(J)+273.16)/(APJ(J,1)*SIG(LDN)+PTOPI))*BYP(J)  2629.   
      CX(J,L)=(AJL(J,L,4)*BYPV(J)*DXV(J)-                               2630.   
     *   AJL(J+1,L,4)*BYPV(J+1)*DXV(J+1))/DXYP(J)                       2631.   
  570 CONTINUE                                                          2632.   
      DX(2,L)=0.                                                        2633.   
      DX(JM,L)=0.                                                       2634.   
      DO 580 J=3,JMM1                                                   2635.   
      DX(J,L)=AMA*COSV(J) + (CX(J,L)-CX(J-1,L) +                        2636.   
     *   .125*(AX(J,L)+AX(J-1,L))*APJ(J,2)*(.25*APJ(J,2)*SIG(L)+PTOPI)* 2637.   
     *   (BX(J,L)-BX(J-1,L)))/DYV(3)                                    2638.   
  580 CONTINUE                                                          2639.   
      CALL JLMAP(56,PL,DX,1.E12,ONES,ONES,LM,2,2)                       2640.   
C****                                                                   2641.   
C**** REFRACTION INDICES FOR WAVES 1 AND 2                              2642.   
C****                                                                   2643.   
      DO 590 L=1,LM                                                     2644.   
      AX(2,L)=0.                                                        2645.   
      AX(JM,L)=0.                                                       2646.   
      LUP=L+1                                                           2647.   
      LDN=L-1                                                           2648.   
      IF(L.EQ.LM) LUP=LM                                                2649.   
      IF(L.EQ.1) LDN=1                                                  2650.   
      DO 590 J=3,JMM1                                                   2651.   
      GBYF=GRAV*DXYP(J)/F(J)                                            2652.   
      SQNBYF=-GBYF*GBYF*(SIG(L)+PTOP*BYP(J)*FIMDA)*                     2653.   
     *   (AJL(J,LUP,17)-AJL(J,LDN,17))*BYD2SG(L)/                       2654.   
     *   (RGAS*(AJL(J,L,1)*BYP(J)+273.16)*(AJL(J,L,17)+1.E-20))         2655.   
      BX(J,L)=SQNBYF                                                    2656.   
  590 CX(J,L)=DX(J,L)*APJ(J,2)/AJL(J,L,4)                               2657.   
      DO 605 M=1,5                                                      2658.   
      SQM=MW(M)*MW(M)                                                   2659.   
      DO 600 J=3,JMM1                                                   2660.   
      BYRCOS=1./(RADIUS*RADIUS*COSV(J)*COSV(J))                         2661.   
      DO 600 L=1,LM                                                     2662.   
      BYHSQ=1./(3434.*(AJL(J,L,1)*BYP(J)+273.16)**2)                    2663.   
  600 AX(J,L)=BX(J,L)*(CX(J,L)-SQM*BYRCOS)-BYHSQ                        2664.   
      IT=ITIT(M)                                                        2665.   
  605 CALL JLMAP(IT,PL,AX,1.E8,ONES,ONES,LM,2,2)                        2666.   
C****                                                                   2667.   
C**** FOURIER ANALYSIS OF GEOPOTENTIAL HEIGHTS FOR WAVE NUMBERS 1 TO 4, 2668.   
C****   AMPLITUDE AND PHASE                                             2669.   
C****                                                                   2670.   
c           LSKIPM=LINECT-63                                            2671.   
c           DO 810 LSKIP=1,LSKIPM                                       2672.   
c 810       WRITE (6,920)                                               2673.   
c           LINECT=63                                                   2674.   
      DO 620 K=1,KM                                                     2675.   
      DO 610 N=1,4                                                      2676.   
      AMPLTD(1,K,N)=0.                                                  2677.   
      AMPLTD(JM,K,N)=0.                                                 2678.   
      PHASE(1,K,N)=0.                                                   2679.   
  610 PHASE(JM,K,N)=0.                                                  2680.   
      DO 620 J=2,JMM1                                                   2681.   
      CALL GETAN (AIJ(1,J,8+K),CN)                                      2682.   
      DO 620 N=1,4                                                      2683.   
      AMPLTD(J,K,N)=SQRT(CN(1,N+1)*CN(1,N+1)+CN(2,N+1)*CN(2,N+1))       2684.   
      PHASE(J,K,N)=(ATAN2(CN(2,N+1),CN(1,N+1))-TWOPI)/N+ELOFIM          2685.   
      IF(PHASE(J,K,N).LE.-.5*TWOPI) PHASE(J,K,N)=PHASE(J,K,N)+TWOPI     2686.   
      PHASE(J,K,N)=-PHASE(J,K,N)                                        2687.   
  620 CONTINUE                                                          2688.   
      SCALE=BYIADA/GRAV                                                 2689.   
      DO 630 N=1,4                                                      2690.   
  630 CALL JLMAP (79+N,PMB,AMPLTD(1,1,N),SCALE,ONES,ONES,KM,2,1)        2691.   
      SCALE=360./TWOPI                                                  2692.   
      DO 640 N=1,4                                                      2693.   
  640 CALL JLMAP (83+N,PMB,PHASE(1,1,N),SCALE,ONES,ONES,KM,2,1)         2694.   
c           LSKIPM=64-LINECT                                            2695.   
c           DO 860 LSKIP=1,LSKIPM                                       2696.   
c 860       WRITE (6,920)                                               2697.   
      RETURN                                                            2698.   
  901 FORMAT (///////                                                   2699.   
     *  '010**14 WATTS = .2067 * 10**19 CALORIES/DAY'/                  2700.   
     *  '010**18 JOULES = .864 * 10**30 GM*CM**2/SEC/DAY'/              2701.   
     *  '0ALL NORTHWARD TRANSPORTS ARE PER UNIT SIGMA')                 2702.   
  920 FORMAT (1X)                                                       2703.   
      END                                                               2704.   
      SUBROUTINE JLMAP (NT,PL,AX,SCALE,SCALEJ,SCALEL,LMAX,JWT,J1)       2801.   
#include "BD2G04.COM"                                                   2802.   
      COMMON U,V,T,P,Q                                                  2803.   
      COMMON/WORK4/MLAT(46),FLAT(46),ASUM(46),FHEM(2),HSUM(2)           2806.   
      COMMON/D2COM/JLAT(46,2),WTJ(46,2,2),                              2807.   
     *  LINECT,JMHALF,INC,IHOUR0,IHOUR,I25,TAUDIF                       2808.   
      COMMON/D2TTL/TITLE(1)                                             2808.1  
      common/nqt/NQMAPS
      DIMENSION AX(JM0,*),ARQX(JM0,*)                                   2809.   
      DIMENSION PL(*),SCALEJ(*),SCALEL(*),SCALJR(*),SCALLR(*)           2810.   
      CHARACTER*4 DASH,WORD(4),TITLE*64                                 2810.1  
      DATA DASH/'----'/,WORD/'SUM','MEAN',' ','.1*'/                    2811.   
C****                                                                   2812.   
C**** PRODUCE A LATITUDE BY LAYER TABLE OF THE ARRAY A                  2813.   
C****                                                                   2814.   
   10 LINECT=LINECT+LMAX+7                                              2815.   
c     IF(LINECT.LE.63) GO TO 20                                         2816.   
c           LSKIPM=64-LINECT+LMAX+7                                     2817.   
c           DO 15 LSKIP=1,LSKIPM                                        2818.   
c  15       WRITE (6,920)                                               2819.   
      JY0=JYEAR0-1900                                                   2820.   
      JY=JYEAR-1900                                                     2821.   
c     WRITE (6,907) (XLABEL(K),K=1,27),JDATE0,JMNTH0,JY0,JDATE,JMONTH,JY2822.   
c     LINECT=LMAX+8                                                     2823.   
c  20 WRITE (6,901) TITLE(NT),(DASH,J=J1,JM,INC)                        2824.   
c     WRITE (6,904) WORD(JWT),(JLAT(JM+J1-J,J1),J=J1,JM,INC)            2825.   
c     WRITE (6,905) (DASH,J=J1,JM,INC)                                  2826.   
      J0=J1-1                                                           2827.   
  100 SDSIG=1.-SIGE(LMAX+1)                                             2828.   
      DO 110 J=J1,JM                                                    2829.   
  110 ASUM(J)=0.                                                        2830.   
      HSUM(1)=0.                                                        2831.   
      HSUM(2)=0.                                                        2832.   
      GSUM=0.                                                           2833.   
      SUMFAC=1.                                                         2834.   
      IWORD=3                                                           2835.   
      IF(NT.NE.1.AND.NT.NE.24.AND.NT.NE.26.AND.NT.NE.28.AND.NT.NE.33)   2836.   
     *  GO TO 112                                                       2837.   
      SUMFAC=10.                                                        2838.   
      IWORD=4                                                           2839.   
  112 DO 140 LX=1,LMAX                                                  2840.   
      L=1+LMAX-LX                                                       2841.   
      FGLOB=0.                                                          2842.   
      DO 130 JHEMI=1,2                                                  2843.   
      FHEM(JHEMI)=0.                                                    2844.   
      DO 120 JH=1,JMHALF                                                2845.   
      J=(JHEMI-1)*(JMHALF-J0)+JH+J0                                     2846.   
      FLAT(J)=AX(J,L)*SCALE*SCALEJ(J)*SCALEL(L)                         2847.   
      MLAT(J)=INT(FLAT(J)+10000.5)-10000                                2848.   
  115 ASUM(J)=ASUM(J)+FLAT(J)*DSIG(L)/SDSIG                             2849.   
  120 FHEM(JHEMI)=FHEM(JHEMI)+FLAT(J)*WTJ(J,JWT,J1)                     2850.   
  130 FGLOB=FGLOB+FHEM(JHEMI)/JWT                                       2851.   
      if(NT.eq.-102)then
          print *,' JLMAP NT=',NT
          print *,AX(12,1),SCALE,SCALEJ(12)
          print *,SCALEL(1),SCALEL(LM)
      endif
c     WRITE (6,902) PL(L),FGLOB,FHEM(2),FHEM(1),                        2852.   
c    *  (MLAT(JM+J1-J),J=J1,JM,INC)                                     2853.   
      do 136 INDEXQ=1,NQTAB
         IF(INQTAB(INDEXQ).NE.NT) GO TO 136                             2860.   
	if(L.eq.-1)then
	print *,' INDEXQ=',INDEXQ,INQTAB(INDEXQ)
	print *,' NT=',nt
	print *,TITLE( NT)
	endif
         J1QT(INDEXQ)=J1
         DO 134 J=J1,JM                                                 2861.   
  134    QTABLE(J,L,INDEXQ)=FLAT(J)                                     2862.   
         QTABLE(JM+1,L,INDEXQ)=FHEM(1)                                  2863.   
         QTABLE(JM+2,L,INDEXQ)=FHEM(2)                                  2864.   
         QTABLE(JM+3,L,INDEXQ)=FGLOB                                    2865.   
  136 CONTINUE
         IF(NT.EQ.5) CALL KEYD2J (L,FLAT)                               2854.
         IF(NT.EQ.7) CALL KEYD2S (L,FLAT)                               2855.
      HSUM(1)=HSUM(1)+FHEM(1)*SUMFAC*DSIG(L)/SDSIG                      2866.   
      HSUM(2)=HSUM(2)+FHEM(2)*SUMFAC*DSIG(L)/SDSIG                      2867.   
  140 GSUM=GSUM+FGLOB*SUMFAC*DSIG(L)/SDSIG                              2868.   
      ASUM(JMHALF+1)=ASUM(JMHALF+1)/J1                                  2869.   
      DO 150 J=J1,JM                                                    2870.   
  150 MLAT(J)=INT(ASUM(J)*SUMFAC+10000.5)-10000                         2871.   
c     WRITE (6,905) (DASH,J=J1,JM,INC)                                  2872.   
      IF(NT.GE.80.AND.NT.LE.87) RETURN                                  2873.   
c     WRITE (6,903) WORD(IWORD),GSUM,HSUM(2),HSUM(1),                   2874.   
c    *  (MLAT(JM+J1-J),J=J1,JM,INC)                                     2875.   
      do 146 INDEXQ=1,NQTAB
         IF(INQTAB(INDEXQ).NE.NT) GO TO 146                             2860.
c         print *,' JLMAP NT=',NT
c        J1QT(INDEXQ)=J1
         L=LM0+1
c        print *,' LM=',L
         DO 144 J=J1,JM                                                 2861.
  144    QTABLE(J,L,INDEXQ)=ASUM(J)*SUMFAC
         QTABLE(JM+1,L,INDEXQ)=HSUM(1)                                  2863.
         QTABLE(JM+2,L,INDEXQ)=HSUM(2)                                  2864.
         QTABLE(JM+3,L,INDEXQ)=GSUM                                     2865.
  146 CONTINUE
         IF(NT.EQ.1) CALL KEYD2T (GSUM,ASUM)                            2876.   
         IF(NT.EQ.18) CALL KEYD2K (ASUM)                                2877.   
         IF(NT.GE.22.AND.NT.LE.33) CALL KEYD2N (NT,ASUM,SUMFAC)         2878.   
c     IF(NT.GE.34) RETURN                                               2879.   
      if(NT.le.43.or.NT.eq.77.or.NT.eq.78.or.NT.eq.79
     *   .or.(NT.ge.88.and.NT.le.92))then
         NQMAPS=NQMAPS+1
         if(NQMAPS.gt.57) then
           print *,' NQMAPS GT 57',NQMAPS
           stop
         endif
         INQMAP(NQMAPS)=NT
         DO 180 J=J1,JM             
  180    QMAPS(J,NQMAPS)=ASUM(J)*SUMFAC         
         QMAPS(JM+1,NQMAPS)=HSUM(1)    
         QMAPS(JM+2,NQMAPS)=HSUM(2)    
         QMAPS(JM+3,NQMAPS)=GSUM       
      end if
      RETURN                                                            2885.   
C****                                                                   2886.   
      ENTRY JLMAPS (NT,PL,AX,SCALE,SCALEJ,SCALEL,LMAX,JWT,J1,           2887.   
     *  ARQX,SCALER,SCALJR,SCALLR)                                      2888.   
c     LINECT=LINECT+LMAX+10                                             2889.   
c     IF(LINECT.LE.63) GO TO 200                                        2890.   
c           LSKIPM=64-LINECT+LMAX+10                                    2891.   
c           DO 190 LSKIP=1,LSKIPM                                       2892.   
c 190       WRITE (6,920)                                               2893.   
c     JY0=JYEAR0-1900                                                   2894.   
c     JY=JYEAR-1900                                                     2895.   
c     WRITE (6,907) (XLABEL(K),K=1,27),JDATE0,JMNTH0,JY0,JDATE,JMONTH,JY2896.   
c     LINECT=LMAX+11                                                    2897.   
c 200 J0=J1-1                                                           2898.   
C**** PRODUCE UPPER STRATOSPHERE NUMBERS FIRST                          2899.   
c     WRITE (6,901) TITLE(NT),(DASH,J=J1,JM,INC)                        2900.   
c     WRITE (6,904) WORD(JWT),(JLAT(JM+J1-J,J1),J=J1,JM,INC)            2901.   
c     WRITE (6,905) (DASH,J=J1,JM,INC)                                  2902.   
c     DO 230 LX=1,3                                                     2903.   
c     L=4-LX                                                            2904.   
c     FGLOB=0.                                                          2905.   
c     DO 220 JHEMI=1,2                                                  2906.   
c     FHEM(JHEMI)=0.                                                    2907.   
c     DO 210 JH=1,JMHALF                                                2908.   
c     J=(JHEMI-1)*(JMHALF-J0)+JH-J0                                     2909.   
c     FLATJ=ARQX(J,L)*SCALER*SCALJR(J)*SCALLR(L)                        2910.   
c     MLAT(J)=INT(FLATJ+10000.5)-10000                                  2911.   
c 210 FHEM(JHEMI)=FHEM(JHEMI)+FLATJ*WTJ(J,JWT,J1)                       2912.   
c 220 FGLOB=FGLOB+FHEM(JHEMI)/JWT                                       2913.   
c 230 WRITE (6,902) PL(L+LM),FGLOB,FHEM(2),FHEM(1),                     2914.   
c    *  (MLAT(JM+J1-J),J=J1,JM,INC)                                     2915.   
      GO TO 100                                                         2916.   
  901 FORMAT ('0',30X,A64/1X,30('-'),24A4)                              2917.   
  902 FORMAT (F6.1,3F8.1,1X,24I4)                                       2918.   
  903 FORMAT (A6,3F8.1,1X,24I4)                                         2919.   
  904 FORMAT (' P(MB) ',A4,' G      NH      SH  ',24I4)                 2920.   
  905 FORMAT (1X,30('-'),24A4)                                          2921.   
  907 FORMAT ('1',27A4,I4,1X,A3,I3,' TO',I3,1X,A3,I3)                   2922.   
  920 FORMAT (1X)                                                       2923.   
      END                                                               2924.   
      SUBROUTINE ILMAP (NT,PL,AX,SCALE,SCALEL,LMAX,JWT,ISHIFT)          3001.   
#include "BD2G04.COM"                                                   3002.   
      COMMON U,V,T,P,Q                                                  3003.   
      COMMON/WORK4/MLON(72),ASUM(72)                                    3004.   
      COMMON/D2TTL/TITLE(1)                                             3005.   
      COMMON/D2COM/JLAT(46,2),WTJ(46,2,2),                              3005.1  
     *  LINECT,JMHALF,INC,IHOUR0,IHOUR,I25,TAUDIF                       3006.   
      DIMENSION AX(36,*)                                                3007.   
      DIMENSION PL(*),SCALEL(*)                                         3008.   
      CHARACTER*4 DASH,WORD(2),TITLE*64                                 3008.1  
      DATA DASH/'----'/,WORD/'SUM','MEAN'/                              3009.   
C****                                                                   3010.   
C**** PRODUCE A LONGITUDE BY LAYER TABLE OF THE ARRAY A                 3011.   
C****                                                                   3012.   
      RETURN                                                            3045.   
  901 FORMAT ('0',30X,A64/1X,14('-'),36A3)                              3046.   
  902 FORMAT (F6.1,F8.1,1X,36I3)                                        3047.   
  903 FORMAT (F14.1,1X,36I3)                                            3048.   
  904 FORMAT (' P(MB)',4X,A4,1X,36I3)                                   3049.   
  905 FORMAT (1X,14('-'),36A3)                                          3050.   
  906 FORMAT (' P(MB)',4X,A4,I2,8I3,I4,26I3)                            3051.   
  907 FORMAT ('1',27A4,I4,1X,A3,I3,' TO',I3,1X,A3,I3)                   3052.   
  920 FORMAT (1X)                                                       3053.   
      END                                                               3054.   
      BLOCK DATA a2                                                     3201.   
C****                                                                   3202.   
C**** TITLES FOR SUBROUTINE DIAG7                                       3203.   
C****                                                                   3204.   
      COMMON/D7COM/TITLE1,TITLE2                                        3205.   
                CHARACTER*64 TITLE1(6)/                                 3206.   
C****                                                               1-6 3207.   
     *'WAVE POWER FOR U NEAR 850 MB AND EQUATOR (DAY*(M/S)**2)',        3208.   
     *'WAVE POWER FOR V NEAR 850 MB AND EQUATOR (DAY*(M/S)**2)',        3210.   
     *'WAVE POWER FOR U NEAR 300 MB AND EQUATOR (10 DAY*(M/S)**2)',     3212.   
     *'WAVE POWER FOR V NEAR 300 MB AND EQUATOR (DAY*(M/S)**2)',        3214.   
     *'WAVE POWER FOR U NEAR 50 MB AND EQUATOR (10 DAY*(M/S)**2)',      3216.   
     *'WAVE POWER FOR V NEAR 50 MB AND EQUATOR (DAY*(M/S)**2)'/         3218.   
                CHARACTER*64 TITLE2(6)/                                 3220.   
C****                                                              7-12 3221.   
     *'WAVE POWER FOR PHI AT 922 MB AND 50 DEG NORTH (10**3 DAY*M**2)', 3222.   
     *'WAVE POWER FOR PHI AT 700 MB AND 50 DEG NORTH (10**3 DAY*M**2)', 3224.   
     *'WAVE POWER FOR PHI AT 500 MB AND 50 DEG NORTH (10**3 DAY*M**2)', 3226.   
     *'WAVE POWER FOR PHI AT 300 MB AND 50 DEG NORTH (10**3 DAY*M**2)', 3228.   
     *'WAVE POWER FOR PHI AT 100 MB AND 50 DEG NORTH (10**4 DAY*M**2)', 3230.   
     *'WAVE POWER FOR PHI AT 10 MB AND 50 DEG NORTH (10**4 DAY*M**2)'/  3232.   
      END                                                               3234.   
      SUBROUTINE DIAG7A                                                 3401.   
C****                                                                   3402.   
C**** THIS SUBROUTINE ACCUMULATES A TIME SEQUENCE FOR SELECTED          3403.   
C**** QUANTITIES AND FROM THAT PRINTS A TABLE OF WAVE FREQUENCIES.      3404.   
C****                                                                   3405.   
#include "BD2G04.COM"                                                   3406.   
      COMMON U,V,T,P,Q                                                  3407.   
      COMMON/WORK3/PHI(IM0,JM0,LM0),HTRD(36,6)                          3408.   
c     COMMON/WORK4/CN(2,37),POWER(120),IPOWER(41),FPE(31)               3409.   
      COMMON/WORK4/CN(2,37),POWER(120),FPE(31),IPOWER(41)
      COMMON/D7COM/TITLE                                                3410.   
      CHARACTER*64 TITLE(12)                                            3411.   
      DIMENSION JLKDEX(6),SCALE(12),PMB(6),GHT(6)                       3412.   
      DATA KM,PMB/6,922.,700.,500.,300.,100.,10./                       3413.   
      DATA NMAX/9/,KQMAX/12/,MMAX/12/,NUAMAX/120/,NUBMAX/15/            3414.   
      DATA SCALE/1.,1., .1,1., .1,1., 4*1.E-3,1.E-4,1.E-5/              3415.   
      DATA GHT/500.,2600.,5100.,8500.,15400.,30000./                    3416.   
      DATA IFIRST/1/                                                    3417.   
      RETURN                                                            3463.   
C****                                                                   3464.   
      ENTRY DIAG7P                                                      3465.   
      RETURN                                                            3540.   
C****                                                                   3541.   
  901 FORMAT ('0',30X,A64,8X,'*1/60 (1/DAY)'/'   PERIOD EASTWARD--',    3542.   
     *   35('---')/' N    -2      *-3   -3.3      -4       -5    -6   -73543.   
     *.5  -10-12-15-20-30-60    60 30 20 15 12 10    7.5    6     5     3544.   
     *   4*   VAR ERR'/'   --',40('---'))                               3545.   
  902 FORMAT (I2,41I3,I4,I4)                                            3546.   
  903 FORMAT ('   --',40('---')/(1X,13F10.4))                           3547.   
  907 FORMAT ('1',27A4,I4,1X,A3,I3,' T0',I3,1X,A3,I3)                   3548.   
  911 FORMAT ('0',30X,A64,8X,'*1/60 (1/DAY)'/'   PERIOD EASTWARD--',    3549.   
     *  35('---')/               ' N   *-4       -5    -6   -7.5  -10-123550.   
     *-15-20-30-60    60 30 20 15 12 10    7.5    6     5        4      3551.   
     * 3.3    3*       2    VAR ERR'/'   --',40('---'))                 3552.   
  920 FORMAT(1X)                                                        3553.   
      END                                                               3554.   
      SUBROUTINE MEM (SERIES,ITM,MMAX,NUAMAX,NUBMAX,POWER,FPE,VAR,PNU)  3801.   
      DIMENSION C(1800),S(1800),B1(62),B2(62),A(12),                    3802.   
     *  AA(11),P(13)                                                    3803.   
      DIMENSION SERIES(*),POWER(*),FPE(*)                               3804.   
C**DOUBLE PRECISION                                                     3805.   
c     REAL*8 PI,ARG,PP,POWERX,P,C,S                                     3806.   
      REAL PI,ARG,PP,POWERX,P,C,S
c     COMPLEX*16 CI,CSUM,SS,A,AA,B1,B2,ANOM,ADEN                        3807.   
      COMPLEX CI,CSUM,SS,A,AA,B1,B2,ANOM,ADEN
      COMPLEX SERIES                                                    3808.   
      PI=3.141592653589793D0                                            3809.   
      CI=DCMPLX(0.D0,1.D0)                                              3810.   
      MMAXP1=MMAX+1                                                     3811.   
C**COSINE AND SINE FUNCTION                                             3812.   
      NUMAX=NUAMAX*NUBMAX                                               3813.   
      DO 20 NU=1,NUMAX                                                  3814.   
      ARG=2.0*PI*DFLOAT(NU)/DFLOAT(NUMAX)                               3815.   
      C(NU)=DCOS(ARG)                                                   3816.   
   20 S(NU)=DSIN(ARG)                                                   3817.   
   50 PP=0.0                                                            3818.   
      DO 60 I=1,ITM                                                     3819.   
   60 PP=PP+SERIES(I)*CONJG(SERIES(I))                                  3820.   
      P(1)=PP/FLOAT(ITM)                                                3821.   
      VAR=P(1)                                                          3822.   
      M=1                                                               3823.   
      B1(1)=SERIES(1)                                                   3824.   
      B2(ITM-1)=SERIES(ITM)                                             3825.   
      ITMM1=ITM-1                                                       3826.   
      DO 70 I=2,ITMM1                                                   3827.   
      B1(I)=SERIES(I)                                                   3828.   
   70 B2(I-1)=SERIES(I)                                                 3829.   
      GO TO 80                                                          3830.   
  100 DO 110 I=1,M                                                      3831.   
  110 AA(I)=A(I)                                                        3832.   
      M=M+1                                                             3833.   
      ITMMM=ITM-M                                                       3834.   
      DO 120 I=1,ITMMM                                                  3835.   
      B1(I)=B1(I)-DCONJG(AA(M-1))*B2(I)                                 3836.   
  120 B2(I)=B2(I+1)-AA(M-1)*B1(I+1)                                     3837.   
   80 ANOM=DCMPLX(0.D0,0.D0)                                            3838.   
      ADEN=DCMPLX(0.D0,0.D0)                                            3839.   
      ITMMM=ITM-M                                                       3840.   
      DO 90 I=1,ITMMM                                                   3841.   
      ANOM=ANOM+DCONJG(B1(I))*B2(I)                                     3842.   
   90 ADEN=ADEN+B1(I)*DCONJG(B1(I))+B2(I)*DCONJG(B2(I))                 3843.   
      A(M)=(ANOM+ANOM)/ADEN                                             3844.   
      P(M+1)=P(M)*(1.0-DCONJG(A(M))*A(M))                               3845.   
      IF(M.EQ.1) GO TO 100                                              3846.   
  130 MM1=M-1                                                           3847.   
      DO 140 I=1,MM1                                                    3848.   
  140 A(I)=AA(I)-A(M)*DCONJG(AA(M-I))                                   3849.   
      IF (M.LT.MMAX) GO TO 100                                          3850.   
C**FINAL PREDICTION ERROR                                               3851.   
      DO 150 M=1,MMAXP1                                                 3852.   
  150 FPE(M)=P(M)*FLOAT(ITM+M-1)/FLOAT(ITM-M+1)                         3853.   
      DO 180 NUA=1,NUAMAX                                               3854.   
      POWERX=0.                                                         3855.   
C**FREQUENCY BAND AVERAGE                                               3856.   
      DO 170 NUB=1,NUBMAX                                               3857.   
      NU=NUB+NUA*NUBMAX+(NUMAX-3*NUBMAX-1)/2                            3858.   
      CSUM=1.                                                           3859.   
      DO 160 M=1,MMAX                                                   3860.   
      NUTM=MOD(NU*M-1,NUMAX)+1                                          3861.   
  160 CSUM=CSUM-A(M)*(C(NUTM)-CI*S(NUTM))                               3862.   
  170 POWERX=POWERX+P(MMAXP1)/(CSUM*DCONJG(CSUM))                       3863.   
      POWER(NUA)=.5*POWERX/FLOAT(NUBMAX)                                3864.   
  180 CONTINUE                                                          3865.   
      PNU=0.0                                                           3866.   
      DO 210 L=1,NUAMAX                                                 3867.   
  210 PNU=PNU+POWER(L)                                                  3868.   
      PNU=PNU/(.5*NUAMAX)                                               3869.   
      RETURN                                                            3870.   
      END                                                               3871.   
      BLOCK DATA a3                                                     4001.   
C****                                                                   4002.   
C**** TITLES, LEGENDS AND CHARACTERS FOR DIAG3                          4003.   
C****                                                                   4004.   
!
! --- Chien Wang some time before 080200
!	To make this work with PGF90: original size of common
!		block was wrong
!
      CHARACTER ACHAR*38,BCHAR*23,CCHAR*38,DCHAR*37,ECHAR*38            4091.   
	character*32 TITLE1, TITLE2, TITLE3
	character*40 legnd1, legnd2
      COMMON/D3COM/TITLE1(3,6),TITLE2(3,6),TITLE3(3,4),	                 4005.   
     *  LEGND1(10),LEGND2(11),ACHAR,BCHAR,CCHAR,                        4006.   
     *  DCHAR,ECHAR
!      COMMON/D3COM/TITLE1(3,6),TITLE2(3,6),TITLE3(3,4),                 4005.   
!     *  LEGND1(10),LEGND2(11),ACHAR,BCHAR,CCHAR,                        4006.   
!     *  DCHAR,ECHAR                                                     4007.   
!        CHARACTER*3 TITLE1,TITLE2,TITLE3
!        CHARACTER*40 LEGND1,LEGND2
!
C****                                                                   4008.   
c                         CHARACTER*32 TITLE1/                          4009.   
                          DATA  TITLE1/   
     1   'TOPOGRAPHY (METERS)',                                         4010.   
     *   'LAND COVERAGE   ',                                            4010.5  
     *   'OCEAN ICE COVERAGE',                                          4010.6  
     *   'SNOW COVERAGE   ',                                            4011.   
     *   'SNOW DEPTH (MM H2O)',                                         4011.5  
     *   'LAND ICE AND FROST COVERAGE',                                 4012.   
C                                                                       4013.   
     7   'PRECIPITATION (MM/DAY)',                                      4014.   
     *   'EVAPORATION (MM/DAY)',                                        4014.5  
     *   'SENSIBLE HEAT FLUX (WATTS/M**2)',                             4015.   
     *   'GROUND WETNESS   ',                                           4015.5  
     *   'GROUND RUNOFF (MM/DAY)',                                      4016.   
     *   'GROUND TEMPERATURE (DEGREES C)',                              4017.   
C                                                                       4018.   
     3   'SURFACE CROSS ISOBAR ANGLE (DEG)',                            4019.   
     *   'JET SPEED (METERS/SEC',                                       4019.5  
     *   'SURFACE WIND SPEED (METERS/SEC)',                             4020.   
     *   'SURF. CROSS ISOBAR ADJ. ANGLE',                               4021.   
     *   'JET DIRECTION (CW NOR)',                                      4021.5  
     *   'SURFACE WIND DIRECTION (CW NOR) '/                            4022.   
c                         CHARACTER*32 TITLE2/                          4023.   
                          DATA TITLE2/  
     9   'TOTAL CLOUD COVER',                                           4024.   
     *   'CONVECTIVE CLOUD COVER',                                      4024.5  
     *   'CLOUD TOP PRESSURE (MB)',                                     4025.   
     *   'LOW LEVEL CLOUDINESS',                                        4025.5  
     *   'MIDDLE LEVEL CLOUDINESS',                                     4026.   
     *   'HIGH LEVEL CLOUDINESS',                                       4027.   
C                                                                       4028.   
     5   'NET RAD. OF PLANET (WATTS/M**2)',                             4029.   
     *   'NET RADIATION AT Z0 (WATTS/M**2)',                            4030.   
     *   'BRIGHTNESS TEMP THRU WNDW(DEG C)',                            4030.5  
     *   'PLANETARY ALBEDO',                                            4031.   
     *   'GROUND ALBEDO   ',                                            4031.5  
     *   'VISUAL ALBEDO   ',                                            4032.   
C                                                                       4033.   
     1   'NET THRML RADIATION (WATTS/M**2)',                            4034.   
     *   'NET HEAT AT Z0 (WATTS/M**2)',                                 4035.   
     *   'TROP STATIC STABILITY (DEG K/KM)',                            4035.5  
     *   'TOTAL NT DRY STAT ENR(10**14 WT)',                            4036.   
     *   'NT DRY STAT ENR BY ST ED(E14 WT)',                            4037.   
     *   'NT DRY STAT ENR BY TR ED(E14 WT)'/                            4038.   
c                         CHARACTER*32 TITLE3/                          4039.   
                          DATA TITLE3/ 
     7   '850 MB HEIGHT (METERS-1500)',                                 4040.   
     *   '700 MB HEIGHT (METERS-3000)',                                 4041.   
     *   '500 MB HEIGHT (METERS-5600)',                                 4042.   
     *   '300 MB HEIGHT (METERS-9500)',                                 4043.   
     *   '100 MB HEIGHT (METERS-16400)',                                4044.   
     *   ' 30 MB HEIGHT (METERS-24000)',                                4045.   
C                                                                       4046.   
     3   'THICKNESS TEMPERATURE 1000-850',                              4047.   
     *   'THICKNESS TEMPERATURE 850-700',                               4048.   
     *   'THICKNESS TEMPERATURE 700-500',                               4049.   
     *   'THICKNESS TEMPERATURE 500-300',                               4050.   
     *   'THICKNESS TEMPERATURE 300-100',                               4051.   
     *   'THICKNESS TEMPERATURE 100-30'/                                4052.   
C****                                                                   4063.   
c     CHARACTER*40 LEGND1/                                              4064.   
      DATA LEGND1/  
     *  '0=0,1=5...9=45,A=50...K=100',                                  4065.   
     *  '0=0...9=90,A=100...I=180...R=270',                             4066.   
     *  '1=.5...9=4.5,A=5...Z=17.5,+=MORE',                             4067.   
     *  '1=1...9=9,A=10...Z=35,+=MORE',                                 4068.   
     *  '1=2...9=18,A=20...Z=70,+=MORE',                                4069.   
C                                                                       4070.   
     *  '1=50...9=450,A=500...Z=1750,+=MORE',                           4071.   
     *  '1=100...9=900,A=1000...Z=3500,+=MORE',                         4072.   
     *  ' ',                                                            4073.   
     *  'A=1...Z=26,3=30...9=90,+=100-150,*=MORE',                      4074.   
     *  '0=0,A=.1...Z=2.6,3=3...9=9,+=10-15'/                           4075.   
c     CHARACTER*40 LEGND2/                                              4076.   
      DATA LEGND2/ 
     *  '-=LESS,Z=-78...0=0...9=27,+=MORE',                             4077.   
     *  '-=LESS,Z=-260...0=0...9=90,+=MORE',                            4078.   
     *  '-=LESS,Z=-520...0=0...9=180,+=MORE',                           4079.   
     *  '-=LESS,Z=-1300...0=0...9=450,+=MORE',                          4080.   
     *  '-=LESS,Z=-2600...0=0...9=900,+=MORE',                          4081.   
C                                                                       4082.   
     *  '-=LESS,Z=-3900...0=0...9=1350,+=MORE',                         4083.   
     *  '-=LESS,Z=-5200...0=0...9=1800,+=MORE',                         4084.   
     *  '-=LESS,9=-.9...0=0,A=.1...Z=2.6,+=MORE',                       4085.   
     *  '-=LESS,9=-45...0=0,A=5...K=45...+=MORE',                       4086.   
     *  '-=LESS,9=-90...0=0,A=10...Z=260,+=MORE',                       4087.   
C                                                                       4088.   
     *  '-=LESS,9=-180...A=20...Z=520,+=MORE'/                          4089.   
C****                                                                   4090.   
      CHARACTER ACHAR*38,BCHAR*23,CCHAR*38,DCHAR*37,ECHAR*38            4091.   
      DATA ACHAR/' 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ+'/              4092.   
      DATA BCHAR/' 0123456789ABCDEFGHIJKX'/                             4095.   
      DATA CCHAR/'-9876543210ABCDEFGHIJKLMNOPQRSTUVWXYZ+'/              4097.   
      DATA DCHAR/' 0ABCDEFGHIJKLMNOPQRSTUVWXYZ3456789+*'/               4100.   
      DATA ECHAR/'-ZYXWVUTSRQPONMLKJIHGFEDCBA0123456789+'/              4103.   
      END                                                               4106.   
      SUBROUTINE DIAG3                                                  4201.   
C****                                                                   4202.   
C**** THIS SUBROUTINE PRODUCES LATITUDE BY LONGITUDE MAPS OF            4203.   
C****                                                                   4204.   
C   K  IND                                                        IDACC 4205.   
C****                                                                   4206.   
C***1      TOPOGRAPHY (M)                                               4207.   
C***2      LAND COVERAGE (10**-2)                                       4208.   
C***3   1  OCEAN ICE COVERAGE (10**-2)                                4 4209.   
C****   2  SNOW COVERAGE (10**-2)                                     4 4210.   
C****   3  SNOW DEPTH (KG H2O/M**2)                                     4211.   
C***6  29  LAND ICE AND FROST COVERAGE (PERCENT)                        4212.   
C****                                                                   4213.   
C***7   5  PRECIPITATION (KG/M**2/86400 S)                            1 4214.   
C****   6  EVAPORATION (KG/M**2/86400 S)                              1 4215.   
C***9   4  SENSIBLE HEAT FLUX (WATTS/METER**2)                          4216.   
C**10   7  BETA, GROUND WETNESS (10**-2)                              3 4217.   
C**11  32  GROUND RUNOFF FROM SURFACE (KG/M**2/86400 S)               1 4218.   
C**12  28  FIRST LAYER GROUND TEMPERATURE (K-273.16)                  1 4219.   
C****                                                                   4220.   
C**13  46  ALPHA0, SURFACE CROSS ISOBAR ANGLE (DEG)                   1 4221.   
C**14  39,40  JET SPEED (M/S)                                         4 4222.   
C**15  36,37  SURFACE WIND SPEED (M/S)                                3 4223.   
C**16  34  SURFACE CROSS ISOBAR ADJUSTMENT ANGLE (DEG)                1 4224.   
C**17  39,40  JET DIRECTION (CLOCKWISE FROM NORTH)                    0 4225.   
C**18  36,37  SURFACE WIND DIRECTION (CLOCKWISE FROM NORTH)           0 4226.   
C****                                                                   4227.   
C**19  19  TOTAL CLOUD COVERAGE (PERCENT)                               4228.   
C**20  17  CLOUD COVERAGE FROM MOIST CONVECTION (PERCENT)               4229.   
C**21  18/19   CLOUD TOP PRESSURE (MILLIBARS)                           4230.   
C**22  41  LOW LEVEL CLOUDINESS (PERCENT)                               4231.   
C**23  42  RMIDDLE LEVEL CLOUDINESS (PERCENT)                           4232.   
C**24  43  HIGH LEVEL CLOUDINESS (PERCENT)                              4233.   
C****                                                                   4234.   
C**25  21+24  RADIATION BALANCE OF PLANET (WATTS/METER**2)              4235.   
C**26  22  RADIATION BALANCE OF GROUND (WATTS/METER**2)                 4236.   
C**27  44  BRIGHTNESS TEMPERATURE THROUGH WINDOW REGION (K-273.16)      4237.   
C**28  24/25  PLANETARY ALBEDO (PERCENT)                                4238.   
C**29  26/27  GROUND ALBEDO (PERCENT)                                   4239.   
C**30  45/25  VISUAL ALBEDO (PERCENT)                                   4240.   
C****                                                                   4241.   
C**31  21  NET THERMAL RADIATION  (WATTA/METER**2)                      4242.   
C**32  23  NET HEAT AT GROUND (WATTS/METER**2)                          4243.   
C**33  31  TROPOSPHERIC STATIC STABILITY                                4244.   
C**34  20  TOTAL NORTH. TRANS. OF DRY STATIC ENERGY  (10**14 WATTS)     4245.   
C**35      STAND. EDDY NORTH. TRANS. OF DRY STATIC ENERGY (10**14 WATTS)4246.   
C**36      TRANS. EDDY NORTH. TRANS. OF DRY STATIC ENERGY (10**14 WATTS)4247.   
C****                                                                   4248.   
C**37  10  850 MB GEOPOTENTIAL HEIGHT (METERS-1500)                     4249.   
C****  11  700 MB GEOPOTENTIAL HEIGHT (METERS-3000)                     4250.   
C****  12  500 MB GEOPOTENTIAL HEIGHT (METERS-5600)                     4251.   
C****  13  300 MB GEOPOTENTIAL HEIGHT (METERS-9500)                     4252.   
C****  14  100 MB GEOPOTENTIAL HEIGHT (METERS-16400)                    4253.   
C****  15  30 MB GEOPOTENTIAL HEIGHT (METERS-24000)                     4254.   
C****                                                                   4255.   
C**43   9,10  THICKNESS TEMPERATURE FROM 1000 TO 850 MB (DEGREES CENT.) 4256.   
C****  10,11  THICKNESS TEMPERATURE FROM 850 TO 700 MB (DEGREES CENT.)  4257.   
C****  11,12  THICKNESS TEMPERATURE FROM 700 TO 500 MB (DEGREES CENT.)  4258.   
C****  12,13  THICKNESS TEMPERATURE FROM 500 TO 300 MB (DEGREES CENT.)  4259.   
C****  13,14  THICKNESS TEMPERATURE FROM 300 TO 100 MB (DEGREES CENT.)  4260.   
C****  14,15  THICKNESS TEMPERATURE FROM 100 TO 30 MB (DEGREES CENT.)   4261.   
#include "BD2G04.COM"                                                   4277.   
      COMMON U,V,T,P,Q                                                  4278.   
      COMMON/WORK2/ENDE16(72,46,2),                                     4279.   
     *  FLAT(3),FGLOBE(3),MLAT(3),MGLOBE(3),GNUM(3),GDEN(3)             4280.   
!
! --- Chien Wang 080200
!	to make this peace of code work with PGF90
!
      CHARACTER*32 TITLE*32
      CHARACTER*4  LEGEND
      CHARACTER ACHAR,BCHAR,CCHAR,DCHAR,ECHAR
      COMMON/D3COM/TITLE(03,16),LEGEND(10,21),ACHAR(38),BCHAR(23),
     *  CCHAR(38),DCHAR(37),ECHAR(38)			      
!      COMMON/D3COM/TITLE(03,16),LEGEND(10,21),ACHAR(38),BCHAR(23),      4281.   
!     *  CCHAR(38),DCHAR(37),ECHAR(38)                                   4282.   
!
C**** ACHAR/ ,0,1,...,8,9,A,B,...,Y,Z,+/                                4285.   
C**** BCHAR/ ,0,1,...,8,9,A,B,...,K,X/                                  4286.   
C**** CCHAR/-,9,8,...,1,0,A,B,...,Y,Z,+/                                4287.   
C**** DCHAR/ ,0,A,B,...,Y,Z,3,4,...,8,9,+,*/                            4288.   
C**** ECHAR/-,Z,Y,...,B,A,0,1,...,8,9,+/                                4289.   
      CHARACTER*1 LINE(72,3),LONGTD(36)                                 4290.   
      DIMENSION IND(48),IA(48),ILEG(3,16),SCALE(48),FAC(48),JGRID(48),  4291.   
     *  PMB(7),GHT(7)                                                   4292.   
      DATA LINE/216*' '/,LONGTD/'+',35*' '/                             4293.   
      DATA IND/3*1,2,3,29,    5, 6, 4, 7,32,28,   46,39,36,34,39,36,    4294.   
     *  19,17,18,41,42,43,   21,22,44,24,26,45,   21,23,31,20, 1, 2,    4295.   
     *  10,11,12,13,14,15,    9,10,11,12,13,14/                         4296.   
      DATA IA/0,0,4*4,        1, 1, 1, 3, 1, 1,    3, 4, 3, 3, 0, 0,    4298.   
     *   2, 2, 0, 2, 2, 2,    2, 1, 2, 0, 0, 0,    2, 1, 4, 4, 4, 4,    4299.   
     *  12*4/                                                           4300.   
      DATA ILEG/7,3*1,9,1,   10,10,12, 1,18,11,   19, 5, 3,19, 2, 2,    4301.   
     *   1, 1, 6, 1, 1, 1,   13,20,11, 1, 1, 1,   13,13, 3,20,20,18,    4302.   
     *  12,13,14,15,15,16,   11,11,11,11,11,11/                         4303.   
      DATA SCALE/1.,3*100.,1.,100.,  3*1.,100.,2*1.,  6*1.,             4305.   
     *  2*100.,1.,3*100.,  3*1.,3*100.,  2*1.,2.,15*1./                 4306.   
      DATA FAC/.01,3*.2,1.,.2,  2*10.,.1,.2,10.,.3333333,               4307.   
     *  .2,.5,2.,.2,2*.1,  2*.2,.02,3*.2,  .05,.1,.3333333,3*.2,        4308.   
     *  2*.05,2.,2*.1,10.,  .1,.05,.02,.01,.01,.006666667,  6*.3333333/ 4309.   
      DATA JGRID/19*1,2,15*1,2,2,1,2,2,8*1/                             4311.   
      DATA PMB/1000.,850.,700.,500.,300.,100.,30./                      4312.   
      DATA GHT/0.,1500.,3000.,5600.,9500.,16400.,24000./                4313.   
C**** INITIALIZE CERTAIN QUANTITIES                                     4314.   
      SHA=RGAS/KAPA                                                     4315.   
      INC=1+JMM1/24                                                     4316.   
      ILINE=36*INC                                                      4317.   
      IQ1=1+IM/(4*INC)                                                  4318.   
      LONGTD(IQ1)=LONGTD(1)                                             4319.   
      IQ2=1+IM/(2*INC)                                                  4320.   
      LONGTD(IQ2)=LONGTD(1)                                             4321.   
      IQ3=1+3*IM/(4*INC)                                                4322.   
      LONGTD(IQ3)=LONGTD(1)                                             4323.   
      BYIM=1./FIM                                                       4324.   
      DTSRCE=NDYN*DT                                                    4325.   
      DTCNDS=NCNDS*DT                                                   4326.   
      SCALE(7)=SDAY/DTCNDS                                              4329.   
      SCALE(8)=SDAY/DTSRCE                                              4330.   
      SCALE(9)=1./DTSRCE                                                4331.   
      SCALE(11)=SDAY/DTSRCE                                             4332.   
      SCALE(13)=360./TWOPI                                              4333.   
      SCALE(16)=360./TWOPI                                              4334.   
      SCALE(26)=1./DTSRCE                                               4335.   
      SCALE(32)=1./DTSRCE                                               4336.   
      SCALE(33)=1.E3*GRAV*EXPBYK(1000.)                                 4337.   
      SCALE(34)=6.25E-14/GRAV                                           4338.   
      SCALE(35)=SCALE(34)                                               4339.   
      SCALE(36)=SCALE(34)                                               4340.   
      DO 70 M=37,42                                                     4341.   
   70 SCALE(M)=1./GRAV                                                  4342.   
      DO 80 M=43,48                                                     4343.   
   80 SCALE(M)=1./(RGAS*DLOG(PMB(M-42)/PMB(M-41)))                      4344.   
C****                                                                   4348.   
      IHOUR0=TOFDY0+.5                                                  4349.   
      IHOUR = TOFDAY + .5                                               4350.   
      TAUDIF=TAU-TAU0                                                   4351.   
      BYIADA=1./(IDACC(4)+1.E-20)                                       4352.   
C****                                                                   4366.   
  160 NDIAG3=46                                                         4367.   
      DO 180 N=1,NDIAG3                                                 4368.   
      IF(JGRID(N).EQ.2) GO TO 180                                       4369.   
      DO 170 I=1,IM                                                     4370.   
      AIJ(I,1,N)=AIJ(1,1,N)                                             4371.   
  170 AIJ(I,JM,N)=AIJ(1,JM,N)                                           4372.   
  180 CONTINUE                                                          4373.   
      DO 610 KPAGE=1,10                                                 4374.   
      IF(KPAGE.GE.7) GO TO 690                                          4375.   
c     WRITE (6,901) XLABEL                                              4376.   
c     WRITE (6,902) IDAY0,IHOUR0,JDATE0,JMNTH0,JYEAR0,TAU0,IDAY,IHOUR,  4377.   
c    *  JDATE,JMONTH,JYEAR,TAU,TAUDIF                                   4378.   
      DO 610 KROW=1,2                                                   4379.   
      KR=2*(KPAGE-1)+KROW                                               4380.   
      IF(KR.GT.16) GO TO 690                                            4381.   
c     WRITE (6,903) (TITLE(K,KR),K=1,03)                                4382.   
      DO 200 KCOLMN=1,3                                                 4383.   
      FGLOBE(KCOLMN)=0.                                                 4384.   
      GNUM(KCOLMN)=0.                                                   4385.   
  200 GDEN(KCOLMN)=0.                                                   4386.   
      DO 550 JX=1,JM,INC                                                4387.   
      J=1+JM-JX                                                         4388.   
      DO 510 KCOLMN=1,3                                                 4389.   
      FLATK=0.                                                          4390.   
      K=3*KR+KCOLMN-3                                                   4391.   
      INDEX=IND(K)                                                      4392.   
      BYIACC=1./(IDACC(IA(K))+1.E-20)                                   4393.   
      GO TO (320,340,400,400,440,400, 440,440,460,400,420,460,          4394.   
     *       420,300,300,420,240,240, 400,400,260,400,400,400,          4395.   
     *       220,420,460,260,260,260, 460,460,380,610,610,610,          4396.   
     *       610,610,610,610,610,610, 610,610,610,610,610,610),K        4397.   
C**** SUM OF TWO ARRAYS                                                 4399.   
  220 DO 230 I=1,IM                                                     4400.   
      A=(AIJ(I,J,21)+AIJ(I,J,24))*SCALE(K)*BYIACC                       4401.   
      FLATK=FLATK+A                                                     4402.   
      N=28.5+A*FAC(K)                                                   4403.   
      IF (N.LT.1 ) N=1                                                  4404.   
      IF (N.GT.38) N=38                                                 4405.   
  230 LINE(I,KCOLMN)=ECHAR(N)                                           4406.   
      GO TO 500                                                         4407.   
C**** WIND DIRECTION                                                    4408.   
  240 IF(J.EQ.1) GO TO 500                                              4409.   
      DO 250 I=1,IM                                                     4410.   
      A=360.*ATAN2(AIJ(I,J,INDEX)+1.E-20,AIJ(I,J,INDEX+1)+1.E-20)/TWOPI 4411.   
      FLATK=FLATK+A                                                     4412.   
      N=2.5+A*FAC(K)                                                    4413.   
      IF(N.LT.2) N=N+36                                                 4414.   
  250 LINE(I,KCOLMN)=ACHAR(N)                                           4415.   
      GO TO 500                                                         4416.   
C**** RATIO OF 2 ARRAYS (MAINLY FOR ALBEDO)                             4417.   
  260 FNUM=0.                                                           4418.   
      FDEN=0.                                                           4419.   
      INDEX2=INDEX+1                                                    4420.   
      IF (INDEX.EQ.45) INDEX2=25                                        4421.   
      DO 270 I=1,IM                                                     4422.   
      A=SCALE(K)*AIJ(I,J,INDEX)/(AIJ(I,J,INDEX2)+1.E-20)                4423.   
      IF(INDEX.EQ.24 .OR. INDEX.EQ.26) A=100.-A                         4424.   
      FNUM=FNUM+AIJ(I,J,INDEX)                                          4425.   
      FDEN=FDEN+AIJ(I,J,INDEX2)                                         4426.   
      N=2.5+A*FAC(K)                                                    4427.   
      IF(A*FAC(K).GE.20.) N=23                                          4428.   
      IF(AIJ(I,J,INDEX2).LE.0.) N=1                                     4429.   
  270 LINE(I,KCOLMN)=ACHAR(N)                                           4430.   
      FLAT(KCOLMN)=SCALE(K)*FNUM/(FDEN+1.E-20)                          4431.   
      IF(INDEX.EQ.24 .OR. INDEX.EQ.26) FLAT(KCOLMN)=100.-FLAT(KCOLMN)   4432.   
      MLAT(KCOLMN)=FLAT(KCOLMN)+.5                                      4433.   
      GNUM(KCOLMN)=GNUM(KCOLMN)+FNUM*DXYP(J)                            4434.   
      GDEN(KCOLMN)=GDEN(KCOLMN)+FDEN*DXYP(J)                            4435.   
      IF(J.GT.INC) GO TO 510                                            4436.   
      FGLOBE(KCOLMN)=SCALE(K)*GNUM(KCOLMN)/(GDEN(KCOLMN)+1.E-20)        4437.   
      IF(INDEX.EQ.24.OR.INDEX.EQ.26) FGLOBE(KCOLMN)=100.-FGLOBE(KCOLMN) 4438.   
      FGLOBE(KCOLMN)=FGLOBE(KCOLMN)*AREAG/(FIM*INC)                     4439.   
      GO TO 510                                                         4440.   
C**** STANDING AND TRANSIENT EDDY NORTHWARD TRANSPORTS OF DSE           4441.   
C 280 IF (SKIPSE.EQ.1.) GO TO 510                                       4442.   
C     DO 290 I=1,IM                                                     4443.   
C     A=ENDE16(I,J,INDEX)*SCALE(K)*BYIACC                               4444.   
C     FLATK=FLATK+A                                                     4445.   
C     N=11.5+A*FAC(K)                                                   4446.   
C     IF(N.LT.1) N=1                                                    4447.   
C     IF(N.GT.38) N=38                                                  4448.   
C 290 LINE(I,KCOLMN)=CCHAR(N)                                           4449.   
C     FLAT(KCOLMN)=FLATK                                                4450.   
C     DAREA=DXYV(J)                                                     4451.   
C     GO TO 505                                                         4452.   
C**** MAGNITUDE OF TWO PERPENDICULAR COMPONENTS                         4453.   
  300 IF(J.EQ.1) GO TO 500                                              4454.   
      DO 310 I=1,IM                                                     4455.   
      A=SQRT(AIJ(I,J,INDEX)**2+AIJ(I,J,INDEX+1)**2)*SCALE(K)*BYIACC     4456.   
      FLATK=FLATK+A                                                     4457.   
      N=2.5+A*FAC(K)                                                    4458.   
      IF(N.GT.38) N=38                                                  4459.   
  310 LINE(I,KCOLMN)=ACHAR(N)                                           4460.   
      GO TO 500                                                         4461.   
C**** SURFACE TOPOGRAPHY                                                4462.   
  320 DO 330 I=1,IM                                                     4463.   
      ZS=FDATA(I,J,1)/GRAV                                              4464.   
      FLATK=FLATK+ZS                                                    4465.   
      N=2.5+.01*ZS                                                      4466.   
      IF (ZS.LE.0.) N=1                                                 4467.   
      IF(N.GT.38) N=38                                                  4468.   
  330 LINE(I,KCOLMN)=ACHAR(N)                                           4469.   
      GO TO 500                                                         4470.   
C**** LAND COVERAGE                                                     4471.   
  340 DO 350 I=1,IM                                                     4472.   
      PLAND=FDATA(I,J,2)*100.                                           4473.   
      FLATK=FLATK+PLAND                                                 4474.   
      N=2.5+PLAND*.2                                                    4475.   
      IF(PLAND.LE.0.) N=1                                               4476.   
      IF(PLAND.GE.100.) N=23                                            4477.   
  350 LINE(I,KCOLMN)=BCHAR(N)                                           4478.   
      GO TO 500                                                         4479.   
C**** THICKNESS TEMPERATURES                                            4480.   
C 360 DO 370 I=1,IM                                                     4481.   
C     A=((AIJ(I,J,INDEX+1)-AIJ(I,J,INDEX))*BYIACC                       4482.   
C    *  +(GHT(INDEX-7)-GHT(INDEX-8))*GRAV)*SCALE(K)-273.16              4483.   
C     FLATK=FLATK+A                                                     4484.   
C     N=28.5+A*FAC(K)                                                   4485.   
C     IF(N.LT.1) N=1                                                    4486.   
C     IF(N.GT.38) N=38                                                  4487.   
C 370 LINE(I,KCOLMN)=ECHAR(N)                                           4488.   
C     GO TO 500                                                         4489.   
C**** POSITIVE QUANTITIES UNIFORMLY SCALED                              4490.   
  380 DO 390 I=1,IM                                                     4491.   
      A=AIJ(I,J,INDEX)*SCALE(K)*BYIACC                                  4492.   
      FLATK=FLATK+A                                                     4493.   
      N=2.5+A*FAC(K)                                                    4494.   
      IF(A.EQ.0.) N=1                                                   4495.   
      IF(N.GT.38) N=38                                                  4496.   
  390 LINE(I,KCOLMN)=ACHAR(N)                                           4497.   
      GO TO 500                                                         4498.   
C**** PERCENTAGES                                                       4499.   
  400 DO 410 I=1,IM                                                     4500.   
      A=AIJ(I,J,INDEX)*SCALE(K)*BYIACC                                  4501.   
      FLATK=FLATK+A                                                     4502.   
      N=2.5+A*FAC(K)                                                    4503.   
      IF(A.LE.0.) N=1                                                   4504.   
      IF(A*FAC(K).GE.20.) N=23                                          4505.   
  410 LINE(I,KCOLMN)=BCHAR(N)                                           4506.   
      GO TO 500                                                         4507.   
C**** SIGNED QUANTITIES UNIFORMLY SCALED (LETTERS +, NUMBERS -)         4508.   
  420 DO 430 I=1,IM                                                     4509.   
      A=AIJ(I,J,INDEX)*SCALE(K)*BYIACC                                  4510.   
      FLATK=FLATK+A                                                     4511.   
      N=11.5+A*FAC(K)                                                   4512.   
      IF(N.LT.1) N=1                                                    4513.   
      IF(N.GT.38) N=38                                                  4514.   
  430 LINE(I,KCOLMN)=CCHAR(N)                                           4515.   
      IF(K.EQ.34) FLATK=FLATK*FIM                                       4516.   
      GO TO 500                                                         4517.   
C**** PRECIPITATION AND EVAPORATION                                     4518.   
  440 DO 450 I=1,IM                                                     4519.   
      A=AIJ(I,J,INDEX)*SCALE(K)*BYIACC                                  4520.   
      FLATK=FLATK+A                                                     4521.   
      N=1                                                               4522.   
      IF(A.LE.0.) GO TO 450                                             4523.   
      N=2.5+A*FAC(K)                                                    4524.   
      IF(N.GT.28) N=(N+263)/10                                          4525.   
      IF(N.GT.35) N=(N+180)/6                                           4526.   
      IF(N.GT.37) N=37                                                  4527.   
  450 LINE(I,KCOLMN)=DCHAR(N)                                           4528.   
      GO TO 500                                                         4529.   
C**** SIGNED QUANTITIES UNIFORMLY SCALED (NUMBERS +, LETTERS -)         4530.   
  460 DO 470 I=1,IM                                                     4531.   
      A=AIJ(I,J,INDEX)*SCALE(K)*BYIACC                                  4532.   
      FLATK=FLATK+A                                                     4533.   
      N=28.5+A*FAC(K)                                                   4534.   
      IF (N.LT.1 ) N=1                                                  4535.   
      IF (N.GT.38) N=38                                                 4536.   
  470 LINE(I,KCOLMN)=ECHAR(N)                                           4537.   
      GO TO 500                                                         4538.   
C**** POSITIVE QUANTITIES NON-UNIFORMLY SCALED                          4539.   
C 480 DO 490 I=1,IM                                                     4540.   
C     A=AIJ(I,J,INDEX)*SCALE(K)*BYIACC                                  4541.   
C     FLATK=FLATK+A                                                     4542.   
C     N=2.5+A*FAC(K)                                                    4543.   
C     IF(N.GE.13) N=(N+123)/10                                          4544.   
C     IF(N.GT.38) N=38                                                  4545.   
C 490 LINE(I,KCOLMN)=ACHAR(N)                                           4546.   
  500 FLAT(KCOLMN)=FLATK*BYIM                                           4547.   
      MLAT(KCOLMN)=INT(FLAT(KCOLMN)+10000.5)-10000                      4548.   
      DAREA=DXYP(J)                                                     4549.   
      IF(JGRID(INDEX).EQ.2) DAREA=DXYV(J)                               4550.   
  505 FGLOBE(KCOLMN)=FGLOBE(KCOLMN)+FLAT(KCOLMN)*DAREA                  4551.   
  510 CONTINUE                                                          4552.   
      GO TO 530
      GO TO (524,520, 520,520, 520,520, 521,520, 526,520, 526,524,      4553.   
     *       527,527, 520,520, 527,527, 527,527),KR                     4554.   
  520 WRITE (6,910) (FLAT(KC),(LINE(I,KC),I=1,ILINE,INC),KC=1,3)        4555.   
      GO TO 530                                                         4556.   
  521 WRITE (6,911) (FLAT(KC),(LINE(I,KC),I=1,ILINE,INC),KC=1,2),       4557.   
     *  MLAT(3),(LINE(I,3),I=1,ILINE,INC)                               4558.   
      GO TO 530                                                         4559.   
  524 WRITE (6,914) MLAT(1),(LINE(I,1),I=1,ILINE,INC),                  4560.   
     *  (FLAT(KC),(LINE(I,KC),I=1,ILINE,INC),KC=2,3)                    4561.   
      GO TO 530                                                         4562.   
  526 WRITE (6,916) (MLAT(KC),(LINE(I,KC),I=1,ILINE,INC),KC=1,2),       4563.   
     *  FLAT(3),(LINE(I,3),I=1,ILINE,INC)                               4564.   
      GO TO 530                                                         4565.   
  527 WRITE (6,917) (MLAT(KC),(LINE(I,KC),I=1,ILINE,INC),KC=1,3)        4566.   
 530  CONTINUE                                                          4567.   
 550  CONTINUE                                                          4575.   
      DO 555 KC=1,3                                                     4576.   
      FGLOBE(KC)=FGLOBE(KC)*FIM*INC/AREAG                               4577.   
  555 MGLOBE(KC)=INT(FGLOBE(KC)+10000.5)-10000                          4578.   
      GO TO 600
      GO TO (574,570, 570,570, 570,570, 571,570, 577,570, 576,570,      4579.   
     *       577,577, 570,570, 577,577, 577,577),KR                     4580.   
  570 WRITE (6,910) (FGLOBE(KC),LONGTD,KC=1,3)                          4581.   
      GO TO 610                                                         4582.   
  571 WRITE (6,911) FGLOBE(1),LONGTD,FGLOBE(2),LONGTD,MGLOBE(3),LONGTD  4583.   
      GO TO 600                                                         4584.   
  574 WRITE (6,914) MGLOBE(1),LONGTD,FGLOBE(2),LONGTD,FGLOBE(3),LONGTD  4585.   
      GO TO 600                                                         4586.   
  576 WRITE (6,916) MGLOBE(1),LONGTD,MGLOBE(2),LONGTD,FGLOBE(3),LONGTD  4587.   
      GO TO 600                                                         4588.   
  577 WRITE (6,917) (MGLOBE(KC),LONGTD,KC=1,3)                          4589.   
  600 WRITE (6,909) ((LEGEND(K,ILEG(KCOLMN,KR)),K=1,10),KCOLMN=1,2),    4590.   
     *  (LEGEND(K,ILEG(3,KR)),K=1,9)                                    4590.1  
  610 CONTINUE                                                          4591.   
  690 CONTINUE                                                          4592.   
C****                                                                   4593.   
C**** PRODUCE FULL PAGE I,J MAPS                                        4594.   
C****                                                                   4595.   
c     WRITE(6,901)XLABEL                                                4596.   
c     WRITE(6,902)IDAY0,IHOUR0,JDATE0,JMNTH0,JYEAR0,TAU0,IDAY,IHOUR,    4597.   
c    *  JDATE,JMONTH,JYEAR,TAU,TAUDIF                                   4598.   
      CALL IJMAP (1,AIJ(1,1,38),BYIADA,JM,IO,IM)                        4599.   
      BYIAC3=1./(IDACC(3)+1.E-20)                                       4600.   
c     WRITE(6,901)XLABEL                                                4601.   
c     WRITE(6,902)IDAY0,IHOUR0,JDATE0,JMNTH0,JYEAR0,TAU0,IDAY,IHOUR,    4602.   
c    *  JDATE,JMONTH,JYEAR,TAU,TAUDIF                                   4603.   
      CALL IJMAP (2,AIJ(1,1,35),BYIAC3,JM,IO,IM)                        4604.   
C     WRITE(6,901)XLABEL                                                4605.   
C     WRITE(6,902)IDAY0,IHOUR0,JDATE0,JMNTH0,JYEAR0,TAU0,IDAY,IHOUR,    4606.   
C    *  JDATE,JMONTH,JYEAR,TAU,TAUDIF                                   4607.   
C     CALL IJMAP (4,AIJ(1,1,8),BYIADA,JM,IO,IM)                         4608.   
C     WRITE(6,901)XLABEL                                                4609.   
C     WRITE(6,902)IDAY0,IHOUR0,JDATE0,JMNTH0,JYEAR0,TAU0,IDAY,IHOUR,    4610.   
C    *  JDATE,JMONTH,JYEAR,TAU,TAUDIF                                   4611.   
C     CALL IJMAP (5,AIJ(1,1,33),BYIADA,JM,IO,IM)                        4612.   
      RETURN                                                            4613.   
C****                                                                   4614.   
  901 FORMAT ('1',33A4)                                                 4615.   
  902 FORMAT ('0',16X,'DAY',I5,', HR',I2,' (',I2,A5,I4,')',F8.0,        4616.   
     *      '    TO    DAY',I5,', HR',I2,' (',I2,A5,I4,')',F8.0,        4617.   
     *  '    DIF',F5.0,' HR')                                           4618.   
  903 FORMAT ('0',6X,A32,13X,A32,13X,A32)                               4619.   
  906 FORMAT ('+',6X,36A1,9X,36A1,9X,36A1)                              4620.   
  909 FORMAT (7X,10A4,5X,10A4,5X,9A4)                                   4621.   
  910 FORMAT (1X,F5.1,1X,36A1,F8.1,1X,36A1,F8.1,1X,36A1)                4622.   
  911 FORMAT (1X,F5.1,1X,36A1,F8.1,1X,36A1,I8,1X,36A1)                  4623.   
  914 FORMAT (1X,I5,1X,36A1,F8.1,1X,36A1,F8.1,1X,36A1)                  4624.   
  916 FORMAT (1X,I5,1X,36A1,I8,1X,36A1,F8.1,1X,36A1)                    4625.   
  917 FORMAT (1X,I5,1X,36A1,I8,1X,36A1,I8,1X,36A1)                      4626.   
      END                                                               4627.   
      SUBROUTINE IJMAP (NT,ARRAY,BYIACC,JM,IO,IM)                       4801.   
      DIMENSION C31(36,24),LON(72),LAT(46),ARRAY(IM,JM)                 4802.   
      CHARACTER*1 LINE(3,72),IDX(12),BLANK,TITLE(5)*48                  4803.   
      DATA IDX/'0','1','2','3','4','5','6','7','8','9','-','*'/         4804.   
      DATA BLANK/' '/                                                   4805.   
C     DATA LINE/216*' '/                                                4806.   
      DATA TITLE/                                                       4807.   
C****                                                                   4808.   
C**** THIS SUBROUTINE PRODUCES NUMERICAL LATITUDE BY LONGITUDE MAPS OF  4809.   
C****                                                                   4810.   
     *  'SEA LEVEL PRESSURE (MB-1000)',                                 4811.   
     *  'SURFACE TEMPERATURE (DEGREES C)',                              4812.   
     *  'INSTANTANEOUS 850 MB HEIGHTS (DEKAMETERS-100)',                4813.   
     *  'SEA LEVEL PRESSURE (MB-1000)  (USING T1)',                     4814.   
     *  'SURFACE TEMPERATURE (DEG C)  (LAPSE RATE FROM T1'/             4815.   
      DATA IFIRST/1/                                                    4815.1  
      IF(IFIRST.NE.1) GO TO 455                                         4815.11 
      IFIRST=0                                                          4815.12 
C****                                                                   4815.2  
C**** INITIALIZE CERTAIN QUANTITIES                                     4815.21 
C****                                                                   4815.22 
      KA=2                                                              4815.24 
c     IO=36                                                             4815.241
c     JM=24                                                             4815.242
c     IM=1                                                              4815.243
      print *,' FROM IJMAP JM=',JM,' IM=',IM,' IO=',IO
      BYIM=1./IM                                                        4815.25 
      INC=1+(JM-1)/24                                                   4815.26 
      ISTEP=INC*2                                                       4815.27 
      IE=36*INC                                                         4815.28 
      LON(1)=-180                                                       4815.29 
      LD=360/IO                                                         4815.3  
      DO 400 I=2,IO                                                     4815.31 
  400 LON(I)=LON(I-1)+LD                                                4815.32 
      DO 450 J=1,JM                                                     4815.33 
  450 LAT(JM-J+1)=INT(.5+(J-1.0)*180./(JM-1))-90                        4815.34 
  455 CONTINUE                                                          4815.35 
C****                                                                   4816.   
c     WRITE(6,900) TITLE(NT)                                            4817.   
c     WRITE (6,910) (I,I=1,IE,INC)                                      4818.   
      DO 300 JX=1,JM                                                    4819.   
      FLAT=0.                                                           4820.   
      J=1+JM-JX                                                         4821.   
      DO 250 I=1,IM                                                     4822.   
      A=ARRAY(I,J)*BYIACC                                               4823.   
      FLAT=FLAT+A                                                       4824.   
      IF (A.LT.999.5.OR.A.GE.-99.5) GO TO 140                           4825.   
      DO 100 K=1,3                                                      4826.   
  100 LINE(K,I)=IDX(12)                                                 4827.   
      GO TO 250                                                         4828.   
  140 DO 150 K=1,3                                                      4829.   
  150 LINE(K,I)=BLANK                                                   4830.   
      JA=NINT(A)                                                        4831.   
      IA=IABS(JA)                                                       4832.   
      IF(IA.GT.99) GO TO 210                                            4833.   
      IF(IA-9) 230,230,220                                              4834.   
  210 LINE(1,I)=IDX(IA/100+1)                                           4835.   
      IA=MOD(IA,100)                                                    4836.   
  220 LINE(2,I)=IDX(IA/10+1)                                            4837.   
      IA=MOD(IA,10)                                                     4838.   
  230 LINE(3,I)=IDX(IA+1)                                               4839.   
      IF(JA.GE.0) GO TO 250                                             4840.   
      IF(JA+9) 240,245,245                                              4841.   
  240 LINE(1,I)=IDX(11)                                                 4842.   
      GO TO 250                                                         4843.   
  245 LINE(2,I)=IDX(11)                                                 4844.   
  250 CONTINUE                                                          4845.   
      FLAT=FLAT*BYIM                                                    4846.   
c     WRITE (6,920) LAT(JX),J,((LINE(K,I),K=1,3),I=1,IE,INC),FLAT       4847.   
c 300 IF(JM.LE.24) WRITE (6,940)                                        4856.   
c     WRITE (6,930) (LON(I),I=1,IM,ISTEP)                               4857.   
  300 continue
      RETURN                                                            4874.   
  900 FORMAT('0',45X,A48)                                               4875.   
  910 FORMAT('0LAT  J/I  ',36I3,5X,'MEAN'//)                            4876.   
  920 FORMAT(2I4,3X,108A1,F9.2)                                         4877.   
  925 FORMAT('+',10X,108A1)                                             4878.   
  930 FORMAT('0  LONG ',18I6)                                           4879.   
  940 FORMAT(' ')                                                       4880.   
      END                                                               4881.   
      BLOCK DATA a4                                                     5001.   
C****                                                                   5002.   
C**** TITLES FOR SUBROUTINE DIAG9                                       5003.   
C****                                                                   5004.   
      COMMON/D9COM/TITLE1,TITLE2,TITLE3,TITLE4                          5005.   
      CHARACTER*32 TITLE1(11)/                                          5006.   
     *  ' INSTANTANE AM (10**9 J*S/M**2) ',                             5007.   
     *  ' CHANGE OF AM BY ADVECTION      ',                             5008.   
     *  ' CHANGE OF AM BY CORIOLIS FORCE ',                             5009.   
     *  ' CHANGE OF AM BY ADVEC + COR    ',                             5010.   
     *  ' CHANGE OF AM BY PRESSURE GRAD  ',                             5011.   
     *  ' CHANGE OF AM BY DYNAMICS       ',                             5012.   
     *  ' CHANGE OF AM BY SURFACE FRIC   ',                             5013.   
     *  ' CHANGE OF AM BY STRATOS DRAG   ',                             5014.   
     *  ' CHANGE OF AM BY FILTER         ',                             5015.   
     *  ' CHANGE OF AM BY DAILY RESTOR   ',                             5016.   
     *  ' SUM OF CHANGES (10**2 J/M**2)  '/                             5017.   
      CHARACTER*32 TITLE2(12)/                                          5018.   
     *  '0INSTANTANEOUS KE (10**3 J/M**2)',                             5019.   
     *  ' CHANGE OF KE BY ADVECTION      ',                             5020.   
     *  ' CHANGE OF KE BY CORIOLIS FORCE ',                             5021.   
     *  ' CHANGE OF KE BY ADVEC + COR    ',                             5022.   
     *  ' CHANGE OF KE BY PRESSURE GRAD  ',                             5023.   
     *  ' CHANGE OF KE BY DYNAMICS       ',                             5024.   
     *  ' CHANGE OF KE BY MOIST CONVEC   ',                             5025.   
     *  ' CHANGE OF KE BY SURF + DRY CONV',                             5026.   
     *  ' CHANGE OF KE BY STRATOS DRAG   ',                             5027.   
     *  ' CHANGE OF KE BY FILTER         ',                             5028.   
     *  ' CHANGE OF KE BY DAILY RESTOR   ',                             5029.   
     *  ' SUM OF CHANGES (10**-3 W/M**2) '/                             5030.   
      CHARACTER*32 TITLE3(5)/                                           5031.   
     *  ' INSTANTANEOUS MASS (KG/M**2)   ',                             5032.   
     *  ' CHANGE OF MASS BY DYNAMICS     ',                             5033.   
     *  ' CHANGE OF MASS BY FILTER       ',                             5034.   
     *  ' CHANGE OF MASS BY DAILY RESTOR ',                             5035.   
     *  ' SUM CHANGES (10**-8 KG/S/M**2) '/                             5036.   
      CHARACTER*32 TITLE4(8)/                                           5037.   
     *  '0INSTANTANE TPE (10**5 J/M**2)  ',                             5038.   
     *  ' CHANGE OF TPE BY DYNAMICS      ',                             5039.   
     *  ' CHANGE OF TPE BY CONDENSATION  ',                             5040.   
     *  ' CHANGE OF TPE BY RADIATION     ',                             5041.   
     *  ' CHANGE OF TPE BY SURFACE INTER ',                             5042.   
     *  ' CHANGE OF TPE BY FILTER        ',                             5043.   
     *  ' CHANGE OF TPE BY DAILY RESTOR  ',                             5044.   
     *  ' SUM OF CHANGES (10**-2 W/M**2) '/                             5045.   
      END                                                               5046.   
      SUBROUTINE DIAG9A (M)                                             5201.   
C****                                                                   5202.   
C**** THIS DIAGNOSTIC ROUTINE KEEPS TRACK OF THE CONSERVATION           5203.   
C**** PROPERTIES OF ANGULAR MOMENTUM, KINETIC ENERGY, MASS, AND         5204.   
C**** TOTAL POTENTIAL ENERGY                                            5205.   
C****                                                                   5206.   
#include "BD2G04.COM"                                                   5207.   
      COMMON U,V,T,P,Q                                                  5208.   
      DIMENSION UX(IO0,JM0,1),VX(IO0,JM0,1)                             5209.   
      COMMON/WORK1/PIT(IM0,JM0),SD(IM0,JM0,LM0-1),PK(IM0,JM0,LM0)       5210.   
      COMMON/WORK2/JLATP(46),JLATV(46),SCALE(36),FGLOB(36),FHEM(2,36),  5211.   
     *  MLAT(46,36),MAREA(46)                                           5212.   
      COMMON/WORK4/PI(46),AM(46),RKE(46),RMASS(46),TPE(46)              5213.   
      COMMON/WORK5/DUT(IO0,JM0,LM0),DVT(IO0,JM0,LM0)                    5214.   
      COMMON/D9COM/TITLE(36)                                            5215.   
      INTEGER NAMOFM(8)/1,6,1,1,7,8,9,10/                               5216.   
      INTEGER NKEOFM(8)/1,17,18,1,19,20,21,22/                          5217.   
      INTEGER NMSOFM(8)/1,25,1,1,1,1,26,27/                             5218.   
      INTEGER NPEOFM(8)/1,30,31,32,33,1,34,35/                          5219.   
      CHARACTER*4 HEMIS(2)/' SH ',' NH '/,DASH/'----'/,TITLE*32         5220.   
C****                                                                   5221.   
C**** THE PARAMETER M INDICATES WHEN DIAG9A IS BEING CALLED             5222.   
C**** M=1  INITIALIZE CURRENT A.M., K.E., MASS, AND T.P.E.              5223.   
C****   2  AFTER DYNAMICS                                               5224.   
C****   3  AFTER CONDENSATION                                           5225.   
C****   4  AFTER RADIATION                                              5226.   
C****   5  AFTER SURFACE INTERACTION AND DRY CONVECTION                 5227.   
C****   6  AFTER STRATOSPHERIC DRAG                                     5228.   
C****   7  AFTER FILTER                                                 5229.   
C****   8  AFTER DAILY RESTORATION                                      5230.   
C****                                                                   5231.   
      RETURN                                                            5332.   
C****                                                                   5333.   
C****                                                                   5334.   
      ENTRY DIAG9D (M,DT1,UX,VX)                                        5335.   
      CALL CLOCKS (MBEGIN)                                              5336.   
C****                                                                   5337.   
C**** THE PARAMETER M INDICATES WHEN DIAG9D IS BEING CALLED             5338.   
C**** M=1  AFTER ADVECTION IN DYNAMICS                                  5339.   
C****   2  AFTER CORIOLIS FORCE IN DYNAMICS                             5340.   
C****   3  AFTER PRESSURE GRADIENT FORCE IN DYNAMICS                    5341.   
C****                                                                   5342.   
      RETURN                                                            5390.   
C****                                                                   5391.   
C****                                                                   5392.   
      ENTRY DIAG9P                                                      5393.   
C****                                                                   5394.   
C**** THIS ENTRY PRODUCES TABLES OF CONSERVATION QUANTITIES             5395.   
C****                                                                   5396.   
      NFILTR=NDYN                                                       5396.1  
      DO 720 J=1,JM                                                     5397.   
      JLATP(J)=INT(.5+(J-1.)*180./JMM1)-90                              5398.   
  720 JLATV(J)=INT(.5+(J-1.5)*180./JMM1)-90                             5399.   
C**** CALCULATE SCALEING FACTORS                                        5400.   
      DTSRCE=DT*NDYN                                                    5401.   
      SCALE(1)=100.E-9*RADIUS/GRAV                                      5402.   
      SCALE(2)=100.E-2*RADIUS/(GRAV*IDACC(1)*DTSRCE+1.E-20)             5403.   
      SCALE(3)=SCALE(2)                                                 5404.   
      SCALE(4)=SCALE(2)                                                 5405.   
      SCALE(5)=SCALE(2)                                                 5406.   
      SCALE(6)=100.E-2*RADIUS/(DTSRCE*GRAV*IDACC(1)+1.E-20)             5407.   
      SCALE(7)=SCALE(6)                                                 5408.   
      SCALE(8)=SCALE(6)                                                 5409.   
      SCALE(9)=100.E-2*RADIUS/(NFILTR*DT*GRAV*IDACC(10)+1.E-20)         5410.   
      SCALE(10)=100.E-2*RADIUS/(SDAY*GRAV*(IDAY-IDAY0)+1.E-20)          5411.   
      SCALE(11)=1.                                                      5412.   
      SCALE(12)=25.E-3/GRAV                                             5413.   
      SCALE(13)=100.E3/(DTSRCE*GRAV*IDACC(1)+1.E-20)                    5414.   
      SCALE(14)=SCALE(13)                                               5415.   
      SCALE(15)=SCALE(13)                                               5416.   
      SCALE(16)=SCALE(13)                                               5417.   
      SCALE(17)=25.E3/(DTSRCE*GRAV*IDACC(1)+1.E-20)                     5418.   
      SCALE(18)=SCALE(17)                                               5419.   
      SCALE(19)=SCALE(17)                                               5420.   
      SCALE(20)=SCALE(17)                                               5421.   
      SCALE(21)=25.E3/(NFILTR*DT*GRAV*IDACC(10)+1.E-20)                 5422.   
      SCALE(22)=25.E3/(SDAY*GRAV*(IDAY-IDAY0)+1.E-20)                   5423.   
      SCALE(23)=1.                                                      5424.   
      SCALE(24)=100.E0/GRAV                                             5425.   
      SCALE(25)=100.E8/(DTSRCE*GRAV*IDACC(1)+1.E-20)                    5426.   
      SCALE(26)=100.E8/(NFILTR*DT*GRAV*IDACC(10)+1.E-20)                5427.   
      SCALE(27)=100.E8/(SDAY*GRAV*(IDAY-IDAY0)+1.E-20)                  5428.   
      SCALE(28)=1.                                                      5429.   
      SCALE(29)=100.E-5/GRAV                                            5430.   
      SCALE(30)=100.E2/(DTSRCE*GRAV*IDACC(1)+1.E-20)                    5431.   
      SCALE(31)=SCALE(30)                                               5432.   
      SCALE(32)=SCALE(30)                                               5433.   
      SCALE(33)=SCALE(30)                                               5434.   
      SCALE(34)=100.E2/(NFILTR*DT*GRAV*IDACC(10)+1.E-20)                5435.   
      SCALE(35)=100.E2/(SDAY*GRAV*(IDAY-IDAY0)+1.E-20)                  5436.   
      SCALE(36)=1.                                                      5437.   
C**** CALCULATE SUMMED QUANTITIES                                       5438.   
      DO 740 J=1,JM                                                     5439.   
      CONSRV(J,4)=CONSRV(J,2)+CONSRV(J,3)                               5440.   
      CONSRV(J,11)=CONSRV(J,6)*SCALE(6)+CONSRV(J,7)*SCALE(7)            5441.   
     *  +CONSRV(J,8)*SCALE(8)+CONSRV(J,9)*SCALE(9)                      5442.   
     *  +CONSRV(J,10)*SCALE(10)                                         5443.   
      CONSRV(J,15)=CONSRV(J,13)+CONSRV(J,14)                            5444.   
      CONSRV(J,23)=CONSRV(J,17)*SCALE(17)+CONSRV(J,18)*SCALE(18)        5445.   
     *  +CONSRV(J,19)*SCALE(19)+CONSRV(J,20)*SCALE(20)                  5446.   
     *  +CONSRV(J,21)*SCALE(21)+CONSRV(J,22)*SCALE(22)                  5447.   
      CONSRV(J,28)=CONSRV(J,25)*SCALE(25)+CONSRV(J,26)*SCALE(26)        5448.   
     *  +CONSRV(J,27)*SCALE(27)                                         5449.   
  740 CONSRV(J,36)=CONSRV(J,30)*SCALE(30)+CONSRV(J,31)*SCALE(31)        5450.   
     *  +CONSRV(J,32)*SCALE(32)+CONSRV(J,33)*SCALE(33)                  5451.   
     *  +CONSRV(J,34)*SCALE(34)+CONSRV(J,35)*SCALE(35)                  5452.   
C**** CALCULATE FINAL ANGULAR MOMENTUM                                  5453.   
      JEQ=1+JM/2                                                        5454.   
      JEQM1=JEQ-1                                                       5455.   
      DO 760 N=1,11                                                     5456.   
      FEQ=CONSRV(JEQ,N)*SCALE(N)*COSV(JEQ)                              5457.   
      FGLOB(N)=FEQ                                                      5458.   
      FHEM(1,N)=.5*FEQ                                                  5459.   
      FHEM(2,N)=.5*FEQ                                                  5460.   
      MLAT(JEQ,N)=INT(FEQ/(FIM*DXYV(JEQ))+1000000.5)-1000000            5461.   
      DO 750 JSH=2,JEQM1                                                5462.   
      JNH=2+JM-JSH                                                      5463.   
      FSH=CONSRV(JSH,N)*SCALE(N)*COSV(JSH)                              5464.   
      FNH=CONSRV(JNH,N)*SCALE(N)*COSV(JNH)                              5465.   
      FGLOB(N)=FGLOB(N)+(FSH+FNH)                                       5466.   
      FHEM(1,N)=FHEM(1,N)+FSH                                           5467.   
      FHEM(2,N)=FHEM(2,N)+FNH                                           5468.   
      MLAT(JSH,N)=INT(FSH/(FIM*DXYV(JSH))+1000000.5)-1000000            5469.   
  750 MLAT(JNH,N)=INT(FNH/(FIM*DXYV(JNH))+1000000.5)-1000000            5470.   
      FGLOB(N)=FGLOB(N)/AREAG                                           5471.   
      FHEM(1,N)=FHEM(1,N)/(.5*AREAG)                                    5472.   
  760 FHEM(2,N)=FHEM(2,N)/(.5*AREAG)                                    5473.   
C**** CALCULATE FINAL KINETIC ENERGY                                    5474.   
      DO 780 N=12,23                                                    5475.   
      FEQ=CONSRV(JEQ,N)*SCALE(N)                                        5476.   
      FGLOB(N)=FEQ                                                      5477.   
      FHEM(1,N)=.5*FEQ                                                  5478.   
      FHEM(2,N)=.5*FEQ                                                  5479.   
      MLAT(JEQ,N)=INT(FEQ/(FIM*DXYV(JEQ))+1000000.5)-1000000            5480.   
      DO 770 JSH=2,JEQM1                                                5481.   
      JNH=2+JM-JSH                                                      5482.   
      FSH=CONSRV(JSH,N)*SCALE(N)                                        5483.   
      FNH=CONSRV(JNH,N)*SCALE(N)                                        5484.   
      FGLOB(N)=FGLOB(N)+(FSH+FNH)                                       5485.   
      FHEM(1,N)=FHEM(1,N)+FSH                                           5486.   
      FHEM(2,N)=FHEM(2,N)+FNH                                           5487.   
      MLAT(JSH,N)=INT(FSH/(FIM*DXYV(JSH))+1000000.5)-1000000            5488.   
  770 MLAT(JNH,N)=INT(FNH/(FIM*DXYV(JNH))+1000000.5)-1000000            5489.   
      FGLOB(N)=FGLOB(N)/AREAG                                           5490.   
      FHEM(1,N)=FHEM(1,N)/(.5*AREAG)                                    5491.   
  780 FHEM(2,N)=FHEM(2,N)/(.5*AREAG)                                    5492.   
C**** CALCUALTE FINAL MASS AND TOTAL POTENTIAL ENERGY                   5493.   
      DO 800 N=24,36                                                    5494.   
      FGLOB(N)=0.                                                       5495.   
      FHEM(1,N)=0.                                                      5496.   
      FHEM(2,N)=0.                                                      5497.   
      DO 790 JSH=1,JEQM1                                                5498.   
      JNH=1+JM-JSH                                                      5499.   
      FSH=CONSRV(JSH,N)*SCALE(N)                                        5500.   
      FNH=CONSRV(JNH,N)*SCALE(N)                                        5501.   
      FGLOB(N)=FGLOB(N)+(FSH+FNH)*DXYP(JSH)                             5502.   
      FHEM(1,N)=FHEM(1,N)+FSH*DXYP(JSH)                                 5503.   
      FHEM(2,N)=FHEM(2,N)+FNH*DXYP(JNH)                                 5504.   
      MLAT(JSH,N)=INT(FSH/FIM+1000000.5)-1000000                        5505.   
  790 MLAT(JNH,N)=INT(FNH/FIM+1000000.5)-1000000                        5506.   
      FGLOB(N)=FGLOB(N)/AREAG                                           5507.   
      FHEM(1,N)=FHEM(1,N)/(.5*AREAG)                                    5508.   
  800 FHEM(2,N)=FHEM(2,N)/(.5*AREAG)                                    5509.   
      AGLOB=1.E-10*AREAG                                                5510.   
      AHEM=1.E-10*(.5*AREAG)                                            5511.   
C**** LOOP OVER HEMISPHERES                                             5512.   
      INC=1+JMM1/24                                                     5513.   
      IHOUR0=TOFDY0+.5                                                  5514.   
      IHOUR=TOFDAY+.5                                                   5515.   
      TAUDIF=TAU-TAU0                                                   5516.   
      DO 870 JHEMIX=1,2                                                 5517.   
      JHEMI=3-JHEMIX                                                    5518.   
c     WRITE (6,901) XLABEL                                              5519.   
c     WRITE (6,902) IDAY0,IHOUR0,JDATE0,JMNTH0,JYEAR0,TAU0,IDAY,IHOUR,  5520.   
c    *  JDATE,JMONTH,JYEAR,TAU,TAUDIF                                   5521.   
      JP1=1+(JHEMI-1)*(JEQ-1)                                           5522.   
      JPM=JHEMI*(JEQ-1)                                                 5523.   
      JV1=2+(JHEMI-1)*(JEQ-2)                                           5524.   
      JVM=JEQ+(JHEMI-1)*(JEQ-2)                                         5525.   
C**** PRODUCE TABLES FOR ANGULAR MOMENTUM AND KINETIC ENERGY            5526.   
c     WRITE (6,903) (DASH,J=JV1,JVM,INC)                                5527.   
c     WRITE (6,904) HEMIS(JHEMI),(JLATV(JV1+JVM-JX),JX=JV1,JVM,INC)     5528.   
c     WRITE (6,903) (DASH,J=JV1,JVM,INC)                                5529.   
c     DO 820 N=1,23                                                     5530.   
c 820 WRITE (6,905) TITLE(N),FGLOB(N),FHEM(JHEMI,N),                    5531.   
c    *  (MLAT(JV1+JVM-JX,N),JX=JV1,JVM,INC)                             5532.   
      DO 830 J=JV1,JVM                                                  5533.   
  830 MAREA(J)=1.E-10*FIM*DXYV(J)+.5                                    5534.   
c     WRITE (6,906) AGLOB,AHEM,(MAREA(JV1+JVM-JX),JX=JV1,JVM,INC)       5535.   
C**** PRODUCE TABLES FOR MASS AND TOTAL POTENTIAL ENERGY                5536.   
c     WRITE (6,907)                                                     5537.   
c     WRITE (6,903) (DASH,J=JP1,JPM,INC)                                5538.   
c     WRITE (6,904) HEMIS(JHEMI),(JLATP(JP1+JPM-JX),JX=JP1,JPM,INC)     5539.   
c     WRITE (6,903) (DASH,J=JP1,JPM,INC)                                5540.   
c     DO 840 N=24,36                                                    5541.   
c 840 WRITE (6,905) TITLE(N),FGLOB(N),FHEM(JHEMI,N),                    5542.   
c    *  (MLAT(JP1+JPM-JX,N),JX=JP1,JPM,INC)                             5543.   
      DO 850 J=JP1,JPM                                                  5544.   
  850 MAREA(J)=1.E-10*FIM*DXYP(J)+.5                                    5545.   
c     WRITE (6,906) AGLOB,AHEM,(MAREA(JP1+JPM-JX),JX=JP1,JPM,INC)       5546.   
c           DO 860 LSKIP=1,10                                           5547.   
c 860       WRITE (6,920)                                               5548.   
  870 CONTINUE                                                          5549.   
      RETURN                                                            5550.   
C****                                                                   5551.   
  901 FORMAT ('1',33A4)                                                 5552.   
  902 FORMAT ('0CONSERVATION QUANTITIES    DAY',I5,', HR',I2,' (',I2,   5553.   
     *  A5,I4,')',F8.0,'    TO    DAY',I5,', HR',I2,' (',I2,A5,I4,')',  5554.   
     *  F8.0,'    DIF',F5.0,' HR'/)                                     5555.   
  903 FORMAT (1X,25('--'),13(A4,'--'))                                  5556.   
  904 FORMAT (35X,'GLOBAL',A7,2X,13I6)                                  5557.   
  905 FORMAT (A32,2F9.2,1X,13I6)                                        5558.   
  906 FORMAT ('0AREA (10**10 M**2)',F22.1,F9.1,1X,13I6)                 5559.   
  907 FORMAT ('0')                                                      5560.   
  920 FORMAT (1X)                                                       5561.   
      END                                                               5562.   
      SUBROUTINE DIAG5A (M25,NDT)                                       6001.   
C****                                                                   6002.   
C**** THIS DIAGNOSTICS ROUTINE PRODUCES A SPECTRAL ANALYSIS OF KINETIC  6003.   
C**** AND AVAILABLE POTENTIAL ENERGIES AND THEIR TRANSFER RATES BY      6004.   
C**** VARIOUS ATMOSPHERIC PROCESSES.                                    6005.   
C****                                                                   6006.   
C**** THE PARAMETER M25 INDICATES WHAT IS STORED IN SPECA(N,M25,KSPHER),6007.   
C**** IT ALSO INDICATES WHEN DIAG5A IS BEING CALLED.                    6008.   
C**** M=1  MEAN STANDING KINETIC ENERGY            BEFORE SOURCES       6009.   
C****   2  MEAN KINETIC ENERGY                     BEFORE DYNAMICS      6010.   
C****   3  MEAN POTENTIAL ENERGY                                        6011.   
C****   4  CONVERSION OF K.E. BY ADVECTION         AFTER ADVECTION      6012.   
C****   5  CONVERSION OF K.E. BY CORIOLIS FORCE    AFTER CORIOLIS TERM  6013.   
C****   6  CONVERSION FROM P.E. INTO K.E.          AFTER PRESS GRAD FORC6014.   
C****   7  CHANGE OF K.E. BY DYNAMICS              AFTER DYNAMICS       6015.   
C****   8  CHANGE OF P.E. BY DYNAMICS                                   6016.   
C****   9  CHANGE OF K.E. BY CONDENSATION          AFTER CONDENSATION   6017.   
C****  10  CHANGE OF P.E. BY CONDENSATION                               6018.   
C****  11  CHANGE OF P.E. BY RADIATION             AFTER RADIATION      6019.   
C****  12  CHANGE OF K.E. BY SURFACE               AFTER SURFACE        6020.   
C****  13  CHANGE OF P.E. BY SURFACE                                    6021.   
C****  14  CHANGE OF K.E. BY FILTER                AFTER FILTER         6022.   
C****  15  CHANGE OF P.E. BY FILTER                                     6023.   
C****  16  CHANGE OF K.E. BY DAILY                 AFTER DAILY          6024.   
C****  17  CHANGE OF P.E. BY DAILY                                      6025.   
C****  18  UNUSED                                                       6026.   
C****  19  LAST KINETIC ENERGY                                          6027.   
C****  20  LAST POTENTIAL ENERGY                                        6028.   
C****                                                                   6029.   
#include "BD2G04.COM"                                                   6030.   
      COMMON U,V,T,P,Q                                                  6031.   
      REAL KE                                                           6032.   
c     REAL*8 TPE,SUMI,SUMT                                              6033.   
      COMMON/WORK1/PIT(IM0,JM0),SD(IM0,JM0,LM0-1),PK(IM0,JM0,LM0)       6034.   
      COMMON/WORK5/DUT(IO0,JM0,LM0),DVT(IO0,JM0,LM0),
     &  FCUV(2,19,JM0,LM0,2),                                           6035.   
     *  FC(2,37),KE(37,8),APE(37,8),VAR(37,4),TPE(2),X(72),             6036.   
     *  SQRTM(72,46),SQRTP(72,46),THJSP(36),THJNP(36),THGM(36),         6037.   
     *  SCALE(20),MN(20),F0(20),FNSUM(20)                               6038.   
      DIMENSION UX(IO0,JM0,*)                                           6039.   
      DIMENSION MTPEOF(20),MAPEOF(8)                                    6040.   
      CHARACTER*8 LATITD(4)/'SOUTHERN','NORTHERN',' EQUATOR','45 NORTH'/6041.   
      CHARACTER*16 SPHERE(2)/'STRATOSPHERE','TROPOSPHERE'/              6042.   
      DATA MTPEOF/0,0,1,0,0,0,0,2,0,3,  4,0,5,0,6,0,7,0,0,8/            6043.   
      DATA MAPEOF/3,8,10,11,13,15,17,20/,IZERO/0/                       6044.   
      NM=1+IM/2                                                         6045.   
      NM8=296                                                           6046.   
      JEQ=1+JM/2                                                        6047.   
      JEQM1=JEQ-1                                                       6048.   
      J45N=2.+.75*JMM1                                                  6049.   
      FIO=IO                                                            6049.5  
      IJL2=IM*JM*LM*2                                                   6050.   
      SHA=RGAS/KAPA                                                     6051.   
      MKE=M25                                                           6052.   
      MAPE=M25                                                          6053.   
C****                                                                   6054.   
C**** KSPHER=1 SOUTHERN STRATOSPHERE       3 NORTHERN STRATOSPHERE      6055.   
C****        2 SOUTHERN TROPOSPHERE        4 NORTHERN TROPOSPHERE       6056.   
C****                                                                   6057.   
C****        5 EQUATORIAL STRATOSPHERE     7 45 DEG NORTH STRATOSPHERE  6058.   
C****        6 EQUATORIAL TROPOSPHERE      8 45 DEG NORTH TROPOSPHERE   6059.   
C****                                                                   6060.   
      GO TO (200,200,810,100,100,  100,200,810,205,810,                 6061.   
     *       296,205,810,205,810,  205,810,810,810,810),M25             6062.   
C****                                                                   6063.   
C**** KINETIC ENERGY                                                    6064.   
C****                                                                   6065.   
C**** TRANSFER RATES FOR KINETIC ENERGY IN THE DYNAMICS                 6066.   
  100 CALL CLOCKS (MBEGIN)                                              6067.   
      DO 110 N=1,NM8                                                    6068.   
  110 KE(N,1)=0.                                                        6069.   
      DO 170 L=1,LM                                                     6070.   
      KSPHER=2                                                          6071.   
      IF (L.GE.LS1) KSPHER=1                                            6072.   
      DO 170 J=2,JM                                                     6073.   
      DO 170 K=IZERO,LM,LM                                              6074.   
      CALL GETAN(DUT(1,J,L+K),FC)                                       6075.   
      DO 120 N=1,NM                                                     6076.   
  120 X(N)=.5*FIM*(FC(1,N)*FCUV(1,N,J,L+K,1)+FC(2,N)*FCUV(2,N,J,L+K,1)) 6077.   
      X(1)=X(1)+X(1)                                                    6078.   
      X(NM)=X(NM)+X(NM)                                                 6079.   
      IF (J.EQ.JEQ) GO TO 150                                           6080.   
      DO 130 N=1,NM                                                     6081.   
  130 KE(N,KSPHER)=KE(N,KSPHER)+X(N)*DSIG(L)                            6082.   
      IF (J.NE.J45N) GO TO 170                                          6083.   
      DO 140 N=1,NM                                                     6084.   
  140 KE(N,KSPHER+4)=KE(N,KSPHER+4)+X(N)*DSIG(L)                        6085.   
      GO TO 170                                                         6086.   
  150 DO 160 N=1,NM                                                     6087.   
      KE(N,KSPHER+4)=KE(N,KSPHER+4)+X(N)*DSIG(L)                        6088.   
      KE(N,KSPHER)=KE(N,KSPHER)+.5D0*X(N)*DSIG(L)                       6089.   
  160 KE(N,KSPHER+2)=KE(N,KSPHER+2)+.5D0*X(N)*DSIG(L)                   6090.   
      IF (K.EQ.LM) KSPHER=KSPHER+2                                      6091.   
  170 CONTINUE                                                          6092.   
      DO 180 KS=1,8                                                     6093.   
      DO 180 N=1,NM                                                     6094.   
  180 SPECA(N,MKE,KS)=SPECA(N,MKE,KS)+KE(N,KS)/NDT                      6095.   
      CALL CLOCKS (MEND)                                                6096.   
      MINC=MBEGIN-MEND                                                  6097.   
      MDIAG=MDIAG+MINC                                                  6098.   
      MDYN=MDYN-MINC                                                    6099.   
      RETURN                                                            6100.   
C**** MASS FOR KINETIC ENERGY                                           6101.   
  200 I=IM                                                              6102.   
      DO 202 J=2,JM                                                     6103.   
      DO 202 IP1=1,IM                                                   6104.   
      SQRTM(I,J)=SQRT(.5*((P(I,J)+P(IP1,J))*DXYS(J)+(P(I,J-1)+          6105.   
     *  P(IP1,J-1))*DXYN(J-1)))                                         6106.   
  202 I=IP1                                                             6107.   
C****                                                                   6108.   
  205 MAPE=MKE+1                                                        6109.   
      DO 206 N=1,NM8                                                    6110.   
  206 KE(N,1)=0.                                                        6111.   
C**** CURRENT KINETIC ENERGY                                            6112.   
      DO 240 L=1,LM                                                     6113.   
      KSPHER=2                                                          6114.   
      IF(L.GE.LS1) KSPHER=1                                             6115.   
      DO 240 J=2,JM                                                     6116.   
      DO 240 K=IZERO,LM,LM                                              6117.   
      DO 210 I=1,IO                                                     6118.   
 210  X(I)=U(1,J,L+K)*SQRTM(1,J)                                        6119.   
c     CALL FRTR (X)                                                     6120.   
      DO 215 N=1,NM                                                     6120.5  
 215  IF(IM.EQ.1) X(N)=X(N)/FIO                                         6120.6  
      IF(J.EQ.JEQ) GO TO 225                                            6121.   
      DO 220 N=1,NM                                                     6122.   
  220 KE(N,KSPHER)=KE(N,KSPHER)+X(N)*DSIG(L)                            6123.   
      IF(J.NE.J45N) GO TO 240                                           6124.   
      DO 222 N=1,NM                                                     6125.   
  222 KE(N,KSPHER+4)=KE(N,KSPHER+4)+X(N)*DSIG(L)                        6126.   
      GO TO 240                                                         6127.   
  225 DO 230 N=1,NM                                                     6128.   
      KE(N,KSPHER+4)=KE(N,KSPHER+4)+X(N)*DSIG(L)                        6129.   
      KE(N,KSPHER)=KE(N,KSPHER)+.5D0*X(N)*DSIG(L)                       6130.   
  230 KE(N,KSPHER+2)=KE(N,KSPHER+2)+.5D0*X(N)*DSIG(L)                   6131.   
      IF(K.EQ.LM) KSPHER=KSPHER+2                                       6132.   
  240 CONTINUE                                                          6133.   
      IF (NDT.EQ.0) GO TO 260                                           6134.   
C**** TRANSFER RATES AS DIFFERENCES OF KINETIC ENERGY                   6135.   
      DO 250 KS=1,8                                                     6136.   
      DO 250 N=1,NM                                                     6137.   
  250 SPECA(N,MKE,KS)=SPECA(N,MKE,KS)+(KE(N,KS)-SPECA(N,19,KS))/NDT     6138.   
  260 DO 270 KS=1,8                                                     6139.   
      DO 270 N=1,NM                                                     6140.   
  270 SPECA(N,19,KS)=KE(N,KS)                                           6141.   
C****                                                                   6142.   
C**** POTENTIAL ENERGY                                                  6143.   
C****                                                                   6144.   
      IF(DOPK.EQ.-1.) GO TO 296                                         6145.   
C**** COMPUTE SQRTP = SQRT(P) AND PK = P**KAPA                          6146.   
      SQRTP1=SQRT(P(1,1))                                               6147.   
      SQRTPM=SQRT(P(1,JM))                                              6148.   
      DO 290 J=2,JMM1                                                   6149.   
      DO 290 I=1,IM                                                     6150.   
  290 SQRTP(I,J)=SQRT(P(I,J))                                           6151.   
      DO 292 I=1,IM                                                     6152.   
      SQRTP(I,1)=SQRTP1                                                 6153.   
  292 SQRTP(I,JM)=SQRTPM                                                6154.   
      IF(DOPK.EQ.0.) GO TO 296                                          6155.   
      DO 294 L=1,LM                                                     6156.   
      DO 294 J=1,JM                                                     6157.   
      DO 294 I=1,IM                                                     6158.   
  294 PK(I,J,L)=EXPBYK(SIG(L)*P(I,J)+PTOP)                              6159.   
  296 DOPK=-1.                                                          6160.   
      DO 298 N=1,NM8                                                    6161.   
  298 APE(N,1)=0.                                                       6162.   
C**** CURRENT AVAILABLE POTENTIAL ENERGY                                6163.   
      LUP=0                                                             6164.   
  300 LUP=LUP+1                                                         6165.   
      THJSP(LUP)=T(1,1,LUP)*SQRTP(1,1)                                  6166.   
      THJNP(LUP)=T(1,JM,LUP)*SQRTP(1,JM)                                6167.   
      THGSUM=FIM*(THJSP(LUP)*DXYP(1)+THJNP(LUP)*DXYP(JM))               6168.   
      DO 320 J=2,JMM1                                                   6169.   
      THJSUM=0.                                                         6170.   
      DO 310 I=1,IM                                                     6171.   
  310 THJSUM=THJSUM+T(I,J,LUP)*SQRTP(I,J)                               6172.   
  320 THGSUM=THGSUM+THJSUM*DXYP(J)                                      6173.   
      THGM(LUP)=THGSUM/AREAG                                            6174.   
      IF(LUP.GE.2) GO TO 350                                            6175.   
      LDN=LUP                                                           6176.   
      L=LUP                                                             6177.   
      GO TO 300                                                         6178.   
  350 DO 360 JHEMI=1,2                                                  6179.   
      DO 360 N=1,NM                                                     6180.   
  360 VAR(N,JHEMI)=0.                                                   6181.   
      VAR(1,1)=.5*(THJSP(L)-THGM(L))**2*DXYP(1)*FIM                     6182.   
      VAR(1,2)=.5*(THJNP(L)-THGM(L))**2*DXYP(JM)*FIM                    6183.   
      GMEAN=((THJSP(LUP)-THJSP(LDN))*DXYP(1)*(SIG(L)*P(1,1)+PTOP)/      6184.   
     *  (SQRTP1*P(1,1)*PK(1,1,L)) + (THJNP(LUP)-THJNP(LDN))*DXYP(JM)*   6185.   
     *  (SIG(L)*P(1,JM)+PTOP)/(SQRTPM*P(1,JM)*PK(1,JM,L)))*FIM          6186.   
      JHEMI=1                                                           6187.   
      DO 388 J=2,JMM1                                                   6188.   
      GMSUM=0.                                                          6189.   
      DO 370 I=1,IO                                                     6190.   
  370 X(I)=T(1,J,L)*SQRTP(1,J)-THGM(L)                                  6191.   
c     INDEX=J+24*(LUP-1)                                                6191.5  
      INDEX=J+JM*(LUP-1)
      GMSUM=FIO*(T(INDEX,1,1)-T(1,J,LDN))*(SIG(L)*P(1,J)+PTOP)/         6192.   
     *  (P(1,J)*PK(1,J,L))                                              6193.   
      GMEAN=GMEAN+GMSUM*DXYP(J)/FIO                                     6194.   
c     CALL FRTR (X)                                                     6195.   
      DO 375 N=1,NM                                                     6195.5  
 375  IF(IM.EQ.1) X(N)=X(N)/FIO                                         6195.6  
      DO 380 N=1,NM                                                     6196.   
  380 VAR(N,JHEMI)=VAR(N,JHEMI)+X(N)*DXYP(J)                            6197.   
      IF(J.NE.JEQ-1) GO TO 384                                          6198.   
      DO 382 N=1,NM                                                     6199.   
  382 VAR(N,3)=X(N)*DXYP(J)                                             6200.   
      JHEMI=2                                                           6201.   
  384 IF(J.NE.J45N-1) GO TO 388                                         6202.   
      DO 386 N=1,NM                                                     6203.   
  386 VAR(N,4)=X(N)*DXYP(J)                                             6204.   
  388 CONTINUE                                                          6205.   
      GMEAN=DSIG(L)*AREAG*(SIG(LDN)-SIG(LUP))/GMEAN                     6206.   
      KS=2                                                              6207.   
      IF(L.GE.LS1) KS=1                                                 6208.   
      DO 400 JHEMI=1,4                                                  6209.   
      DO 390 N=1,NM                                                     6210.   
  390 APE(N,KS)=APE(N,KS)+VAR(N,JHEMI)*GMEAN                            6211.   
  400 KS=KS+2                                                           6212.   
      IF(L.EQ.LM) GO TO 450                                             6213.   
      LDN=L                                                             6214.   
      L=LUP                                                             6215.   
      IF(LUP.LT.LM) GO TO 300                                           6216.   
      GO TO 350                                                         6217.   
C**** CURRENT TOTAL POTENTIAL ENERGY                                    6218.   
  450 DO 480 JHEMI=1,2                                                  6219.   
      JP=1+JMM1*(JHEMI-1)                                               6220.   
      SUMT=0.                                                           6221.   
      DO 455 L=1,LM                                                     6222.   
  455 SUMT=SUMT+T(1,JP,L)*PK(1,JP,L)*DSIG(L)                            6223.   
      TPE(JHEMI)=FIM*DXYP(JP)*(FDATA(1,JP,1)*(P(1,JP)+PTOP)+            6224.   
     *  SUMT*SHA*P(1,JP))                                               6225.   
      DO 480 JH=2,JEQM1                                                 6226.   
      J=JH+(JEQM1-1)*(JHEMI-1)                                          6227.   
      SUMI=0.                                                           6228.   
      DO 470 I=1,IM                                                     6229.   
      SUMT=0.                                                           6230.   
      DO 460 L=1,LM                                                     6231.   
  460 SUMT=SUMT+T(I,J,L)*PK(I,J,L)*DSIG(L)                              6232.   
  470 SUMI=SUMI+FDATA(I,J,1)*(P(I,J)+PTOP)+SUMT*SHA*P(I,J)              6233.   
  480 TPE(JHEMI)=TPE(JHEMI)+SUMI*DXYP(J)                                6234.   
      IF (NDT.EQ.0) GO TO 520                                           6235.   
      MTPE=MTPEOF(MAPE)                                                 6236.   
C**** TRANSFER RATES AS DIFFERENCES FOR POTENTIAL ENERGY                6237.   
      DO 510 KS=1,8                                                     6238.   
      DO 510 N=1,NM                                                     6239.   
  510 SPECA(N,MAPE,KS)=SPECA(N,MAPE,KS)+(APE(N,KS)-SPECA(N,20,KS))/NDT  6240.   
      ATPE(MTPE,1)=ATPE(MTPE,1)+(TPE(1)-ATPE(8,1))/NDT                  6241.   
      ATPE(MTPE,2)=ATPE(MTPE,2)+(TPE(2)-ATPE(8,2))/NDT                  6242.   
  520 DO 530 KS=1,8                                                     6243.   
      DO 530 N=1,NM                                                     6244.   
  530 SPECA(N,20,KS)=APE(N,KS)                                          6245.   
      ATPE(8,1)=TPE(1)                                                  6246.   
      ATPE(8,2)=TPE(2)                                                  6247.   
      CALL CLOCKS (MNOW)                                                6248.   
      MDIAG=MDIAG+MLAST-MNOW                                            6249.   
      MLAST=MNOW                                                        6250.   
      IF(M25.NE.2) RETURN                                               6251.   
C**** ACCUMULATE MEAN KINETIC ENERGY AND MEAN POTENTIAL ENERGY          6252.   
      IDACC(7)=IDACC(7)+1                                               6253.   
      DO 550 KS=1,8                                                     6254.   
      DO 550 N=1,NM                                                     6255.   
      SPECA(N,2,KS)=SPECA(N,2,KS)+KE(N,KS)                              6256.   
  550 SPECA(N,3,KS)=SPECA(N,3,KS)+APE(N,KS)                             6257.   
      ATPE(1,1)=ATPE(1,1)+TPE(1)                                        6258.   
      ATPE(1,2)=ATPE(1,2)+TPE(2)                                        6259.   
      RETURN                                                            6260.   
C****                                                                   6261.   
      ENTRY DIAG5F(UX)                                                  6262.   
C**** FOURIER COEFFICIENTS FOR CURRENT WIND FIELD                       6263.   
C****                                                                   6264.   
      CALL CLOCKS (MBEGIN)                                              6265.   
      DO 590 K=IZERO,LM,LM                                              6266.   
      DO 590 L=1,LM                                                     6267.   
      DO 590 J=2,JM                                                     6268.   
  590 CALL GETAN(UX(1,J,L+K),FCUV(1,1,J,L+K,1))                         6269.   
      IDACC(6)=IDACC(6)+1                                               6270.   
      CALL CLOCKS (MEND)                                                6271.   
      MINC=MBEGIN-MEND                                                  6272.   
      MDIAG=MDIAG+MINC                                                  6273.   
      MDYN=MDYN-MINC                                                    6274.   
      RETURN                                                            6275.   
C****                                                                   6276.   
      ENTRY DIAG5P                                                      6277.   
C**** THIS ENTRY PRINTS THE SPECTRAL ANALYSIS TABLES                    6278.   
C****                                                                   6279.   
      NM=1+IM/2                                                         6280.   
      IF(SKIPSE.GE.1.) GO TO 600                                        6281.   
      JEQ=1+JM/2                                                        6282.   
      J45N=2.+.75*JMM1                                                  6283.   
      FIO=IO                                                            6283.5  
C****                                                                   6284.   
C**** STANDING KINETIC ENERGY                                           6285.   
C****                                                                   6286.   
      DO 710 K=1,8                                                      6287.   
      DO 710 N=1,NM                                                     6288.   
  710 SPECA(N,1,K)=0.                                                   6289.   
      DO 770 L=1,LM                                                     6290.   
      KSPHER=2                                                          6291.   
      IF(L.GE.LS1) KSPHER=1                                             6292.   
      DO 770 J=2,JM                                                     6293.   
      FACTOR=DSIG(L)*FIM*DXYV(J)/APJ(J,2)                               6294.   
      DO 770 K=IZERO,LM,LM                                              6295.   
      DO 720 I=1,IO                                                     6296.   
 720  X(I)=AIJL(1,J,L+K,1)                                              6297.   
c     CALL FRTR (X)                                                     6298.   
      DO 725 N=1,NM                                                     6298.5  
 725  IF(IM.EQ.1) X(N)=X(N)/FIO                                         6298.6  
      IF(J.EQ.JEQ) GO TO 750                                            6299.   
      DO 730 N=1,NM                                                     6300.   
  730 SPECA(N,1,KSPHER)=SPECA(N,1,KSPHER)+X(N)*FACTOR                   6301.   
      IF(J.NE.J45N) GO TO 770                                           6302.   
      DO 740 N=1,NM                                                     6303.   
  740 SPECA(N,1,KSPHER+4)=SPECA(N,1,KSPHER+4)+X(N)*FACTOR               6304.   
      GO TO 770                                                         6305.   
  750 DO 760 N=1,NM                                                     6306.   
      SPECA(N,1,KSPHER+4)=SPECA(N,1,KSPHER+4)+X(N)*FACTOR               6307.   
      SPECA(N,1,KSPHER)=SPECA(N,1,KSPHER)+.5*X(N)*FACTOR                6308.   
  760 SPECA(N,1,KSPHER+2)=SPECA(N,1,KSPHER+2)+.5*X(N)*FACTOR            6309.   
      IF(K.EQ.LM) KSPHER=KSPHER+2                                       6310.   
  770 CONTINUE                                                          6311.   
C****                                                                   6312.   
  600 SCALE(1)=25.E-17/(GRAV*IDACC(4)+1.E-20)                           6313.   
      SCALE(19)=100.E-17/GRAV                                           6314.   
      SCALE(20)=SCALE(19)*RGAS                                          6315.   
      SCALE(2)=SCALE(19)/(IDACC(7)+1.E-20)                              6316.   
      SCALE(3)=SCALE(2)*RGAS                                            6317.   
      SCALE(4)=100.E-12/(GRAV*DT*IDACC(6)+1.E-20)                       6318.   
      SCALE(5)=SCALE(4)                                                 6319.   
      SCALE(6)=SCALE(4)                                                 6320.   
      SCALE(7)=100.E-12/(GRAV*DT*(IDACC(7)+1.E-20))                     6321.   
      SCALE(8)=SCALE(7)*RGAS                                            6322.   
      SCALE(9)=100.E-12/(GRAV*DT*(IDACC(8)+1.E-20))                     6323.   
      SCALE(10)=SCALE(9)*RGAS                                           6324.   
      SCALE(11)=SCALE(10)                                               6325.   
      SCALE(12)=SCALE(9)                                                6326.   
      SCALE(13)=SCALE(10)                                               6327.   
      SCALE(14)=100.E-12/(GRAV*DT*(IDACC(10)+1.E-20))                   6328.   
      SCALE(15)=SCALE(14)*RGAS                                          6329.   
      SCALE(16)=100.E-12/(GRAV*DT*(IDAY-IDAY0+1.E-20))                  6330.   
      SCALE(17)=SCALE(16)*RGAS                                          6331.   
      SCALE(18)=0.                                                      6332.   
      IUNITJ=17                                                         6333.   
      IUNITW=12                                                         6334.   
      IHOUR0=TOFDY0+.5                                                  6335.   
      IHOUR=TOFDAY+.5                                                   6336.   
      DO 690 KPAGE=1,4                                                  6337.   
C**** WRITE HEADINGS                                                    6338.   
c     WRITE (6,901) XLABEL                                              6339.   
c     WRITE (6,902) IDAY0,IHOUR0,JDATE0,JMNTH0,JYEAR0,IDAY,IHOUR,JDATE, 6340.   
c    *  JMONTH,JYEAR,IUNITJ,IUNITW                                      6341.   
      DO 670 KROW=1,2                                                   6342.   
c     IF(JM.GE.25.AND.KROW.EQ.2) WRITE (6,901)                          6343.   
c     WRITE (6,903) LATITD(KPAGE),SPHERE(KROW)                          6344.   
      KSPHER=2*(KPAGE-1)+KROW                                           6345.   
C**** WRITE KINETIC AND AVAILABLE POTENTIAL ENERGY BY WAVE NUMBER       6346.   
C     DO 610 M=1,20                                                     6347.   
C     F0(M)=SPECA(1,M,KSPHER)*SCALE(M)                                  6348.   
C     MN(M)=INT(F0(M)+100000.5)-100000                                  6349.   
C 610 FNSUM(M)=0.                                                       6350.   
C     WRITE (6,904) MN                                                  6351.   
      IF(IM.EQ.1) GO TO 635                                             6351.5  
      DO 630 N=2,NM                                                     6352.   
      KSPHER=2*(KPAGE-1)+KROW                                           6353.   
      DO 620 M=1,20                                                     6354.   
      FNM=SPECA(N,M,KSPHER)*SCALE(M)                                    6355.   
      MN(M)=INT(FNM+100000.5)-100000                                    6356.   
  620 FNSUM(M)=FNSUM(M)+FNM                                             6357.   
      NM1=N-1                                                           6358.   
         IF(KSPHER.LT.8) GO TO 630                                      6359.   
         FN26=SPECA(N,2,6)*SCALE(2)                                     6360.   
         FN28=SPECA(N,2,8)*SCALE(2)                                     6361.   
         CALL KEYD5A (NM1,FN26,FN28)                                    6362.   
c 630 WRITE (6,905) NM1,MN                                              6363.   
  630 continue
 635  CONTINUE                                                          6363.5  
      DO 640 M=1,20                                                     6364.   
  640 MN(M)=INT(FNSUM(M)+100000.5)-100000                               6365.   
c     WRITE (6,906) MN                                                  6366.   
      DO 650 M=1,20                                                     6367.   
  650 MN(M)=INT(FNSUM(M)+F0(M)+100000.5)-100000                         6368.   
c     WRITE (6,907) MN                                                  6369.   
  670 CONTINUE                                                          6370.   
      IF(KPAGE.GE.3) GO TO 690                                          6371.   
C**** WRITE TOTAL POTENTIAL ENERGY                                      6372.   
      DO 680 MTPE=1,8                                                   6373.   
      MAPE=MAPEOF(MTPE)                                                 6374.   
  680 MN(MTPE)=INT(ATPE(MTPE,KPAGE)*SCALE(MAPE)/RGAS+1000000.5)         6375.   
     *  -1000000                                                        6376.   
c     WRITE (6,909) (MN(MTPE),MTPE=1,8)                                 6377.   
      IF(KPAGE.NE.2) GO TO 690                                          6378.   
      DO 685 M=1,20                                                     6379.   
  685 SCALE(M)=SCALE(M)*10.                                             6380.   
      IUNITJ=16                                                         6381.   
      IUNITW=11                                                         6382.   
  690 CONTINUE                                                          6383.   
      RETURN                                                            6384.   
C****                                                                   6385.   
  810 WRITE (6,910) M25                                                 6386.   
      STOP 29                                                           6387.   
  901 FORMAT ('1',33A4)                                                 6388.   
  902 FORMAT ('0** SPECTRAL ANALYSIS **   DAY',I5,', HR',I2,' (',I2,    6389.   
     *  A5,I4,')    TO    DAY',I5,', HR',I2,' (',I2,A5,I4,              6390.   
     *  ')   UNITS 10**',I2,' JOULES AND 10**',I2,' WATTS')             6391.   
  903 FORMAT ('0',50X,A8,A9,A8/                                         6392.   
     *  13X,'MEAN',19X,'DYNAMICS',25X,'SOURCES',16X,'FILTER',8X,        6393.   
     *     'DAILY',16X,'LAST'/                                          6394.   
     *'   N    SKE   KE   APE    KADV  KCOR   P-K  KDYN  PDYN   ',      6395.   
     *     'KCNDS PCNDS   PRAD KSURF PSURF   KFIL  PFIL   KGMP  PGMP',  6396.   
     *     12X,'KE   APE')                                              6397.   
  904 FORMAT ( '0  0',I7,I5,I6,I8,4I6,I8,I6,I7,2I6,I7,I6,I7,2I6,I8,I6/) 6398.   
  905 FORMAT (     I4,I7,I5,I6,I8,4I6,I8,I6,I7,2I6,I7,I6,I7,2I6,I8,I6)  6399.   
  906 FORMAT (' EDDY',I6,I5,I6,I8,4I6,I8,I6,I7,2I6,I7,I6,I7,2I6,I8,I6)  6400.   
  907 FORMAT ('0TOTL',I6,I5,I6,I8,4I6,I8,I6,I7,2I6,I7,I6,I7,2I6,I8,I6)  6401.   
  908 FORMAT ('0')                                                      6402.   
  909 FORMAT (/'0TPE',I18,I32,I14,I7,I12,2I13,I20)                      6403.   
  910 FORMAT ('0INCORRECT VALUE OF M WHEN CALLING DIAG5A.  M=',I5)      6404.   
      END                                                               6405.   
      BLOCK DATA a5                                                     7001.   
C****                                                                   7002.   
C**** TITLES FOR SUBROUTINE DIAG6                                       7003.   
C****                                                                   7004.   
      COMMON/D6COM/TITLE                                                7005.   
      CHARACTER*8 TITLE(50)/                                            7006.   
     *  '0INC SW ',' P ALBD ',' G ALBD ',' ABS ATM',' E CNDS ',         7007.   
     *  '0SRF PRS',' PT 5   ',' PT 4   ',' PT 3   ',' PT 2   ',         7008.   
     *  ' PT 1   ',' TS     ',' TG1    ','0Q 5   ',' Q 4   ',           7009.   
     *  ' Q 3   ',' Q 2   ',' Q 1   ',' QS    ',' QG    ',              7010.   
     *  '0CLD 6  ',' CLD 5  ',' CLD 4  ',' CLD 3  ',' CLD 2  ',         7011.   
     *  ' CLD 1  ',' COVER  ','0SW ON G',' LW AT G',' SNSB HT',         7012.   
     *  ' LAT HT ',' HEAT Z0','0UG*10  ',' VG*10  ',' WG*10  ',         7013.   
     *  ' US*10  ',' VS*10  ',' WS*10  ',' ALPHA0 ','0RIS1*E2',         7014.   
     *  ' RIGS*E2',' CDM*E4 ',' CDH*E4 ',' DGS*10 ',' EDS1*10',         7015.   
     *  '0PPBL   ',' DC FREQ',' LDC*10 ','0PRC*10 ',' EVP*10 '/         7016.   
      END                                                               7017.   
      SUBROUTINE DIAG6                                                  7201.   
C****                                                                   7202.   
C**** THIS SUBROUTINE PRINTS THE DIURNAL CYCLE OF SOME QUANTITIES       7203.   
C****                                                                   7204.   
#include "BD2G04.COM"                                                   7205.   
      COMMON U,V,T,P,Q                                                  7206.   
      DIMENSION SCALE(50),MHOUR(25)                                     7207.   
      COMMON/D6COM/TITLE(50)                                            7208.   
      CHARACTER*8 TITLE                                                 7208.1  
      DATA SCALE/1.,2*100.,2*1.,  5*1.,  3*1.,2*1.E5,  5*1.E5,          7209.   
     *  5*100.,  2*100.,3*1.,  2*1.,3*10.,  3*10.,1.,100.,              7210.   
     *  100.,2*1.E4,2*10.,  1.,100.,10.,2*1./                           7211.   
C****                                                                   7212.   
      IF(IDAY.LE.IDAY0) RETURN                                          7213.   
      DTCNDS=NCNDS*DT                                                   7214.   
      DTSURF=NDYN*DT/NSURF                                              7215.   
      BYIDAC=1./(IDAY-IDAY0)                                            7216.   
      SCALE(5)=100.*RGAS/(KAPA*GRAV*DTCNDS)                             7217.   
      SCALE(28)=1./DTSURF                                               7218.   
      SCALE(29)=1./DTSURF                                               7219.   
      SCALE(30)=1./DTSURF                                               7220.   
      SCALE(31)=1./DTSURF                                               7221.   
      SCALE(32)=1./DTSURF                                               7222.   
      SCALE(39)=360./TWOPI                                              7223.   
      SCALE(49)=100.*100.*SDAY/(DTCNDS*GRAV)                            7224.   
      SCALE(50)=100.*SDAY/DTSURF                                        7225.   
C****                                                                   7226.   
      DO 500 KR=1,4                                                     7227.   
      JY0=JYEAR0-1900                                                   7228.   
      JY=JYEAR-1900                                                     7229.   
c     WRITE (6,901) (XLABEL(K),K=1,27),JDATE0,JMNTH0,JY0,JDATE,JMONTH,JY7230.   
c     WRITE (6,903) NAMD6(KR),IJD6(1,KR),IJD6(2,KR),(I,I=1,24)          7231.   
      DO 500 KQ=1,50                                                    7232.   
      IF(KQ.EQ.48) GO TO 200                                            7233.   
C**** NORMAL QUANTITIES                                                 7234.   
      AVE=0.                                                            7235.   
      DO 120 IH=1,24                                                    7236.   
      AVE=AVE+ADAILY(IH,KQ,KR)                                          7237.   
  120 MHOUR(IH)=INT(ADAILY(IH,KQ,KR)*SCALE(KQ)*BYIDAC+100000.5)-100000  7238.   
      MHOUR(25)=INT(AVE/24.*SCALE(KQ)*BYIDAC+100000.5)-100000           7239.   
      GO TO 500                                                         7240.   
C**** RATIO OF TWO QUANTITIES                                           7241.   
  200 AVEN=0.                                                           7242.   
      AVED=0.                                                           7243.   
      DO 220 IH=1,24                                                    7244.   
      AVEN=AVEN+ADAILY(IH,KQ,KR)                                        7245.   
      AVED=AVED+ADAILY(IH,KQ-1,KR)                                      7246.   
  220 MHOUR(IH)=ADAILY(IH,KQ,KR)*SCALE(KQ)/(ADAILY(IH,KQ-1,KR)+1.E-20)  7247.   
     *  +.5                                                             7248.   
      MHOUR(25)=AVEN*SCALE(KQ)/(AVED+1.E-20)+.5                         7249.   
c 500 WRITE (6,904) TITLE(KQ),MHOUR                                     7250.   
  500 continue
      RETURN                                                            7251.   
C****                                                                   7252.   
  901 FORMAT ('1',27A4,I4,1X,A3,I3,' TO',I3,1X,A3,I3)                   7253.   
  903 FORMAT ('0',A4,I2,',',I2,' ',I2,23I5,'  AVE')                     7254.   
  904 FORMAT (2A4,25I5)                                                 7255.   
      END                                                               7256.   
      SUBROUTINE DIAG4A                                                 8001.   
C****                                                                   8002.   
C**** THIS SUBROUTINE PRODUCES A TIME HISTORY OF ENERGIES               8003.   
C****                                                                   8004.   
#include "BD2G04.COM"                                                   8005.   
      COMMON U,V,T,P,Q                                                  8006.   
      COMMON/WORK1/SUM(20),IK(20)                                       8007.   
      DIMENSION SCALE(20),EHIST(20)                                     8010.   
      IF(IDACC(4).LE.0.OR.IDACC(7).LE.0) RETURN                         8011.   
      JEQ=2.+.5*JMM1                                                    8012.   
      NM=1+IM/2                                                         8013.   
C****                                                                   8014.   
C**** LOAD ENERGIES INTO TIME HISTORY ARRAY                             8015.   
C****                                                                   8016.   
      IDACC5=IDACC(5)+1                                                 8017.   
      IF(SKIPSE.EQ.1.) GO TO 540                                        8018.   
C**** CALCULATE CURRENT SEKE                                            8019.   
      BYIADA=1./IDACC(4)                                                8020.   
      DO 530 L=1,LM                                                     8021.   
      KS=5                                                              8022.   
      IF (L.GE.LS1) KS=15                                               8023.   
      DO 530 J=2,JM                                                     8024.   
      PU4TI=0.                                                          8025.   
      PV4TI=0.                                                          8026.   
      SKE4I=0.                                                          8027.   
      DO 510 I=1,IM                                                     8028.   
      PU4TI=PU4TI+AIJL(I,J,L,1)                                         8029.   
      PV4TI=PV4TI+AIJL(I,J,L,2)                                         8030.   
  510 SKE4I=SKE4I+(AIJL(I,J,L,1)*AIJL(I,J,L,1)                          8031.   
     *            +AIJL(I,J,L,2)*AIJL(I,J,L,2))/AIJ(I,J,8)              8032.   
      SEKE=(SKE4I-(PU4TI*PU4TI+PV4TI*PV4TI)/APJ(J,2))*DXYV(J)*BYIADA    8033.   
      IF(J.EQ.JEQ) GO TO 520                                            8034.   
      ENERGY(KS,IDACC5)=ENERGY(KS,IDACC5)+SEKE*DSIG(L)                  8035.   
      GO TO 530                                                         8036.   
  520 ENERGY(KS,IDACC5)=ENERGY(KS,IDACC5)+.5*SEKE*DSIG(L)               8037.   
      ENERGY(KS+1,IDACC5)=ENERGY(KS+1,IDACC5)+.5*SEKE*DSIG(L)           8038.   
      IF(K.EQ.2) KS=KS+1                                                8039.   
  530 CONTINUE                                                          8040.   
C**** OTHER ENERGIES COME FROM LATEST SPECTRAL ANALYSIS                 8041.   
  540 ENERGY(1,IDACC5)=SPECA(1,19,2)                                    8042.   
      ENERGY(2,IDACC5)=SPECA(1,19,4)                                    8043.   
      ENERGY(7,IDACC5)=SPECA(1,20,2)                                    8044.   
      ENERGY(8,IDACC5)=SPECA(1,20,4)                                    8045.   
      ENERGY(11,IDACC5)=SPECA(1,19,1)                                   8046.   
      ENERGY(12,IDACC5)=SPECA(1,19,3)                                   8047.   
      ENERGY(17,IDACC5)=SPECA(1,20,1)                                   8048.   
      ENERGY(18,IDACC5)=SPECA(1,20,3)                                   8049.   
      IF(IM.EQ.1) GO TO 955                                             8049.5  
      DO 550 N=2,NM                                                     8050.   
      ENERGY(3,IDACC5)=ENERGY(3,IDACC5)+SPECA(N,19,2)                   8051.   
      ENERGY(4,IDACC5)=ENERGY(4,IDACC5)+SPECA(N,19,4)                   8052.   
      ENERGY(9,IDACC5)=ENERGY(9,IDACC5)+SPECA(N,20,2)                   8053.   
      ENERGY(10,IDACC5)=ENERGY(10,IDACC5)+SPECA(N,20,4)                 8054.   
      ENERGY(13,IDACC5)=ENERGY(13,IDACC5)+SPECA(N,19,1)                 8055.   
      ENERGY(14,IDACC5)=ENERGY(14,IDACC5)+SPECA(N,19,3)                 8056.   
      ENERGY(19,IDACC5)=ENERGY(19,IDACC5)+SPECA(N,20,1)                 8057.   
  550 ENERGY(20,IDACC5)=ENERGY(20,IDACC5)+SPECA(N,20,3)                 8058.   
 955  CONTINUE                                                          8058.5  
      IDACC(5)=IDACC5                                                   8059.   
      RETURN                                                            8060.   
C****                                                                   8061.   
      ENTRY DIAG4                                                       8062.   
C**** THIS ENTRY PRODUCES A TIME HISTORY TABLE OF ENERGIES              8063.   
C****                                                                   8064.   
      IDACC5=IDACC(5)                                                   8065.   
      IF(IDACC5.LE.0) RETURN                                            8066.   
      SCALE(1)=100.E-18/GRAV                                            8067.   
      SCALE(2)=SCALE(1)                                                 8068.   
      SCALE(3)=SCALE(1)                                                 8069.   
      SCALE(4)=SCALE(1)                                                 8070.   
      SCALE(5)=.125*SCALE(1)                                            8071.   
      SCALE(6)=SCALE(5)                                                 8072.   
      SCALE(7)=SCALE(1)*RGAS                                            8073.   
      SCALE(8)=SCALE(7)                                                 8074.   
      SCALE(9)=SCALE(7)                                                 8075.   
      SCALE(10)=SCALE(7)                                                8076.   
      SCALE(11)=SCALE(1)                                                8077.   
      SCALE(12)=SCALE(1)                                                8078.   
      SCALE(13)=SCALE(1)                                                8079.   
      SCALE(14)=SCALE(1)                                                8080.   
      SCALE(15)=SCALE(5)                                                8081.   
      SCALE(16)=SCALE(5)                                                8082.   
      SCALE(17)=SCALE(7)                                                8083.   
      SCALE(18)=SCALE(7)                                                8084.   
      SCALE(19)=SCALE(7)                                                8085.   
      SCALE(20)=SCALE(7)                                                8086.   
C****                                                                   8087.   
      IHOUR0=TOFDY0+.5                                                  8088.   
      IHOUR=TOFDAY+.5                                                   8089.   
c     WRITE (6,901) XLABEL                                              8090.   
c     WRITE (6,902) IDAY0,IHOUR0,JDATE0,JMNTH0,JYEAR0,IDAY,IHOUR,JDATE, 8091.   
c    *  JMONTH,JYEAR                                                    8092.   
      DO 110 K=1,20                                                     8093.   
  110 SUM(K)=0.                                                         8094.   
c     WRITE (6,903)                                                     8095.   
      DTAUD4=DT*NDA4/3600.                                              8096.   
      TAUX=TAU0+DT*NCNDS/3600.                                          8097.   
      DO 200 I=1,IDACC5                                                 8098.   
      IDAYX=(TAUX+.001)/24.                                             8099.   
      TOFDYX=TAUX-24.*IDAYX                                             8100.   
      DO 150 K=1,20                                                     8101.   
      IK(K)=ENERGY(K,I)*SCALE(K)+.5                                     8102.   
  150 SUM(K)=SUM(K)+ENERGY(K,I)                                         8103.   
c     WRITE (6,904) IDAYX,TOFDYX,IK                                     8104.   
  200 TAUX=TAUX+DTAUD4                                                  8105.   
      DO 250 K=1,20                                                     8106.   
         EHIST(K)=SUM(K)*SCALE(K)/IDACC5                                8107.   
  250 IK(K)=EHIST(K)+.5                                                 8108.   
c     WRITE (6,905) IK                                                  8109.   
      LSKIPM=54-IDACC5                                                  8110.   
c     DO 260 LSKIP=1,LSKIPM                                             8111.   
c 260 WRITE (6,920)                                                     8112.   
         CALL KEYD4 (IK)                                                8113.   
      RETURN                                                            8114.   
C****                                                                   8115.   
  901 FORMAT ('1',33A4)                                                 8116.   
  902 FORMAT ('0** ENERGY HISTORY **   DAY',I5,', HR',I3,' (',I2,A5,I5, 8117.   
     *  ')    TO    DAY',I5,', HR',I3,' (',I2,A5,I5,                    8118.   
     *  ')    UNITS OF 10**18 JOULES')                                  8119.   
  903 FORMAT ('0',15X,21('-'),' TROPOSPHERE ',22('-'),5X,21('-'),       8120.   
     *  ' STRATOSPHERE ',21('-')/8X,2(11X,'ZKE',8X,'EKE',7X,'SEKE',9X,  8121.   
     *  'ZPE',10X,'EPE')/3X,'DAY  HOUR     SH   NH    SH   NH    SH   NH8122.   
     *    SH    NH     SH    NH      SH   NH    SH   NH    SH   NH     S8123.   
     *H    NH     SH    NH'/1X,132('='))                                8124.   
  904 FORMAT (I6,F6.1,1X,3(I6,I5),2(I7,I6),2X,3(I6,I5),2(I7,I6))        8125.   
  905 FORMAT (1X,132('=')/8X,'MEAN ',3(I6,I5),2(I7,I6),2X,3(I6,I5),     8126.   
     *  2(I7,I6))                                                       8127.   
  920 FORMAT (1X)                                                       8128.   
      END                                                               8129.   
      SUBROUTINE DIAG8(IPFLAG)                                          8601.   
      RETURN                                                            8602.   
      ENTRY ENQJOB                                                      8603.   
      RETURN                                                            8604.   
      END                                                               8605.   
      SUBROUTINE DIAG10(IPFLAG)                                         8801.   
      RETURN                                                            8802.   
      END                                                               8803.   
      SUBROUTINE DIAGKS                                                 9001.   
C****                                                                   9002.   
C**** THIS SUBROUTINE PRODUCES A SUMMARY OF KEY NUMBERS CALCULATED IN   9003.   
C**** OTHER DIAGNOSTIC SUBROUTINES                                      9004.   
C****                                                                   9005.   
C**** CONTENTS OF KEYNR                                                 9006.   
C****  1 MONTH                                                          9007.   
C****  2 TOTAL CLOUD COVER (PERCENT)                                    9008.   
C****  3 SNOW COVER--NORTHERN HEMSIPHERE (PERCENT)                      9009.   
C****  4 ICE COVER--NORTHERN HEMISPHERE (PERCENT)                       9010.   
C****  5 PLANETARY ALBEDO (PERCENT)                                     9011.   
C****  6 SOLAR RADIATION ABSORBED BY ATMOSPHERE (WT/M**2)               9012.   
C****  7 SOLAR RADIATION ABSORBED BY PLANET (WT/M**2)                   9013.   
C****  8 NET HEAT AT GROUND (WT/M**2)                                   9014.   
C****  8 ANGULAR MOMENTUM PER UNIT AREA (10**10 J*SEC/M**2)             9015.   
C****  9 EVAPORATION (.1 MM/DAY)                                        9016.   
C****  9 PRECIPITATION (.1 MM/DAY)                                      9017.   
C**** 10 SENSIBLE HEAT FLUX INTO GROUND (ABS.VALUE)                     9018.   
C**** 11 LATENT HEAT FLUX INTO GROUND (ABS.VALUE)                       9019.   
C**** 12 MEAN GROUND TEMPERATURE (DEGREES K)                            9020.   
C**** 13 MEAN GLOBAL ATMOSPHERIC TEMPERATURE (DEGREES K)                9021.   
C**** 14 MERID. TEMPERATURE GRADIENT (N.HEMISPHERE)                     9022.   
C**** 15 MERID. TEMPERATURE GRADIENT (S.HEMISPHERE)                     9023.   
C**** 16 MEAN TROPOSPHERIC EKE-NORTHERN HEMISPHERE                      9024.   
C**** 17 MEAN TROPOSPHERIC EKE-SOUTHERNN HEMISPHERE                     9025.   
C**** 18 MEAN TROPOSPHERIC ZKE-NORTHERN HEMISPHERE                      9026.   
C**** 19 MEAN TROPOSPHERIC ZKE-SOUTHERN HEMISPHERE                      9027.   
C**** 20 MEAN TROPOSPHERIC EPE-NORTHERN HEMISPHERE                      9028.   
C**** 21 MEAN TROPOSPHERIC ZPE-NORTHERN HEMISPHERE                      9029.   
C**** 22 MEAN EDDY KINETIC ENERGY AT EQUATOR                            9030.   
C**** 23 MAX. MEAN EDDY KINETIC ENERGY IN MID NORTH LATITUDES           9031.   
C**** 24 MAX. ZONAL WIND (U COMPONENT) IN TROPOSPHERE (NH), M/SEC       9032.   
C**** 25 LATITUDE CORRESPONDING TO 24                                   9033.   
C**** 26 MAX. ZONAL WIND (U COMPONENT) IN TROPOSPHERE (SH), M/SEC       9034.   
C**** 27 LATITUDE CORRESPONDING TO 26                                   9035.   
C**** 28-30: 29 IS LARGEST VALUE OF STREAM FUNCTION, POSITIVE OR        9036.   
C****    NEGATIVE; 28 AND 30 ARE THE MAGNITUDES OF THE LARGEST VALUES OF9037.   
C****    OPPOSITE SIGN TO THE NORTH AND SOUTH RESPECTIVELY              9038.   
C**** 31 EKE 'SLOPE' AT EQUATOR, TROPOSPHERE (10**16 JOULES)            9039.   
C**** 32 EKE 'SLOPE' AT 45 DEGREES NORTH, TROPOSPHERE (10**16 JOULES)   9040.   
C****   33-39 REFER TO NORTHERN HEMISPHERE ONLY                         9041.   
C**** 33 MAX.NORTHWARD TRANS. OF DRY STATIC ENERGY BY STANDING EDDIES   9042.   
C**** 34 MAX.NORTHWARD TRANS. OF DRY STATIC ENERGY BY EDDIES            9043.   
C**** 35 MAX. TOTAL NORTH. TRANS. OF DRY STATIC ENERGY                  9044.   
C**** 36 MAX.NORTHWARD TRANS. OF STATIC ENERGY BY EDDIES                9045.   
C**** 37 MAX.TOTAL NORTH. TRANS. OF STATIC ENERGY                       9046.   
C**** 38 LATITUDE CORRESPONDING TO 37                                   9047.   
C**** 39 MAX. NORTH. TRANS. OF ANGULAR MOMENTUM BY STANDING EDDIES      9048.   
C**** 40 MAX. NORTH. TRANS. OF ANGULAR MOMENTUM BY EDDIES               9049.   
C**** 41 MAX. TOTAL NORTH. TRANS. OF ANGULAR MOMENTUM                   9050.   
C**** 42 LATITUDE CORRESPONDING TO 41                                   9051.   
C****                                                                   9052.   
#include "BD2G04.COM"                                                   9053.   
      COMMON U,V,T,P,Q                                                  9054.   
c     7/27/04
c     DIMENSION KEYDS(20)
      DIMENSION KEYDS(42)
      COMMON/WORK4/FKEY(46,36)                                          9057.   
      COMMON/D2COM/JLAT(46,2)                                           9058.   
C     COMMON/KEYS/KEYNR(42,50)                                          9059.   
      DIMENSION ASUM(*),FLAT(*),IK(*)                                   9060.   
      CHARACTER*4 IC,JAN,CKEYNR(42,50)                                  9060.1  
      EQUIVALENCE (CKEYNR,KEYNR)                                        9060.2  
      DATA IC/'IC'/,JAN/'JAN'/                                          9061.   
C****                                                                   9062.   
C**** ENTRIES CALLED FROM DIAG1                                         9063.   
C****                                                                   9064.   
      ENTRY KEYD1 (N,FGLOB,FNH)                                         9065.   
      GO TO (                100,100,100,110,100, 100,100,100,100,115,  9066.   
     *  100,100,120,125,100, 100,100,130,100,135, 100,100,100,100,100,  9067.   
     *  100,100,100,100,140, 145,100,100,100,100, 100,100,100,100,100,  9068.   
     *  100,100,100,150,100, 100,100,100,100,100, 100,100,100,100,100,  9069.   
     *  100,100,100,155),N                                              9070.   
  100 RETURN                                                            9071.   
  110 KEYNR(6,KEYCT)=INT(FGLOB+.5)                                      9072.   
      RETURN                                                            9073.   
  115 KEYNR(7,KEYCT)=INT(FGLOB+.5)                                      9074.   
      RETURN                                                            9075.   
  120 KEYNR(10,KEYCT)=INT(.5-FGLOB)                                     9076.   
      RETURN                                                            9077.   
  125 KEYNR(11,KEYCT)=INT(.5-FGLOB)                                     9078.   
      RETURN                                                            9079.   
  130 KEYNR(12,KEYCT)=INT(.1*FGLOB+.5)                                  9080.   
      RETURN                                                            9081.   
  135 KEYNR(9,KEYCT)=INT(10.*FGLOB+.5)                                  9082.   
      RETURN                                                            9083.   
  140 KEYNR(4,KEYCT)=INT(FNH+.5)                                        9084.   
      RETURN                                                            9085.   
  145 KEYNR(3,KEYCT)=INT(FNH+.5)                                        9086.   
      RETURN                                                            9087.   
  150 KEYNR(8,KEYCT)=INT(FGLOB+100000.5)-100000                         9088.   
      RETURN                                                            9089.   
  155 KEYNR(2,KEYCT)=INT(FGLOB+.5)                                      9090.   
      RETURN                                                            9091.   
C****                                                                   9092.   
      ENTRY KEYD1A (FGLOB)                                              9093.   
      KEYNR(5,KEYCT)=INT(10.*FGLOB+.5)                                  9094.   
      RETURN                                                            9095.   
C****                                                                   9096.   
C**** ENTRIES CALLED FROM DIAG2 VIA JLMAP                               9097.   
C****                                                                   9098.   
      ENTRY KEYD2T (GSUM,ASUM)                                          9099.   
C**** TEMPERATURES                                                      9100.   
      JEQ=2.+.5*JMM1                                                    9101.   
      TEQ=.5*(ASUM(JEQ-1)+ASUM(JEQ))                                    9102.   
      X60=TWOPI/(12.*DLAT)                                              9103.   
      J60=.5+X60                                                        9104.   
      A=DXYP(J60+1)*(X60+.5-J60)                                        9105.   
      TSOU=ASUM(J60+1)*A                                                9106.   
      TNOR=ASUM(JM-J60)*A                                               9107.   
      DO 210 J=1,J60                                                    9108.   
      A=A+DXYP(J)                                                       9109.   
      TSOU=TSOU+ASUM(J)*DXYP(J)                                         9110.   
  210 TNOR=TNOR+ASUM(JM+1-J)*DXYP(J)                                    9111.   
      KEYNR(14,KEYCT)=INT(TEQ-TNOR/A+.5)                                9112.   
      KEYNR(15,KEYCT)=INT(TEQ-TSOU/A+.5)                                9113.   
      KEYNR(13,KEYCT)=INT(.1*GSUM-.5)                                   9114.   
      RETURN                                                            9115.   
C****                                                                   9116.   
      ENTRY KEYD2J (L,FLAT)                                             9117.   
C**** JET STREAMS                                                       9118.   
      IF(L.LT.LM) GO TO 220                                             9119.   
      DO 216 LL=1,LM                                                    9120.   
      IF((PSF-PTOP)*SIG(LL)+PTOP.LT.200.) GO TO 218                     9121.   
  216 CONTINUE                                                          9122.   
  218 LMAX=LL-1                                                         9123.   
  220 IF(L.GT.LMAX) RETURN                                              9124.   
      USLM=-999999.                                                     9125.   
      DO 222 J=3,JEQ                                                    9126.   
      IF(FLAT(J).LT.USLM) GO TO 222                                     9127.   
      USLM=FLAT(J)                                                      9128.   
      JMAX=J                                                            9129.   
  222 CONTINUE                                                          9130.   
      CEPT=.5*(FLAT(JMAX-1)-FLAT(JMAX+1))/                              9131.   
     *  (FLAT(JMAX-1)-2.*FLAT(JMAX)+FLAT(JMAX+1))                       9132.   
      LSLM=INT((JMAX-1.5+CEPT)*DLAT*360/TWOPI+.5)-90                    9133.   
      UNLM=-999999.                                                     9134.   
      DO 224 J=JEQ,JMM1                                                 9135.   
      IF(FLAT(J).LT.UNLM) GO TO 224                                     9136.   
      UNLM=FLAT(J)                                                      9137.   
      JMAX=J                                                            9138.   
  224 CONTINUE                                                          9139.   
      CEPT=.5*(FLAT(JMAX-1)-FLAT(JMAX+1))/                              9140.   
     *  (FLAT(JMAX-1)-2.*FLAT(JMAX)+FLAT(JMAX+1))                       9141.   
      LNLM=INT((JMAX-1.5+CEPT)*DLAT*360/TWOPI+.5)-90                    9142.   
      IF(L.LT.LMAX) GO TO 226                                           9143.   
      USM=USLM                                                          9144.   
      LSM=LSLM                                                          9145.   
      UNM=UNLM                                                          9146.   
      LNM=LNLM                                                          9147.   
      RETURN                                                            9148.   
  226 IF(USLM.LT.USM) GO TO 228                                         9149.   
      USM=USLM                                                          9150.   
      LSM=LSLM                                                          9151.   
  228 IF(UNLM.LT.UNM) GO TO 230                                         9152.   
      UNM=UNLM                                                          9153.   
      LNM=LNLM                                                          9154.   
  230 IF(L.NE.1) RETURN                                                 9155.   
      KEYNR(24,KEYCT)=.1*UNM+.5                                         9156.   
      KEYNR(25,KEYCT)=LNM                                               9157.   
      KEYNR(26,KEYCT)=.1*USM+.5                                         9158.   
      KEYNR(27,KEYCT)=-LSM                                              9159.   
C****                                                                   9160.   
      ENTRY KEYD2S (L,FLAT)                                             9161.   
C**** STREAM FUNCTION                                                   9162.   
      DO 290 J=2,JM                                                     9163.   
  290 FKEY(J,L)=FLAT(J)                                                 9164.   
      IF(L.NE.1) RETURN                                                 9165.   
  300 SAVE=0.                                                           9166.   
      HS=0.                                                             9167.   
      HN=0.                                                             9168.   
      DO 310 K=1,LM                                                     9169.   
      DO 310 I=2,JM                                                     9170.   
      CHECK=ABS(FKEY(I,K))                                              9171.   
      IF(CHECK.LT.SAVE) GO TO 310                                       9172.   
      SAVE=CHECK                                                        9173.   
      INDEX=I                                                           9174.   
      KNDEX=K                                                           9175.   
  310 CONTINUE                                                          9176.   
      SAVE=FKEY(INDEX,KNDEX)                                            9177.   
      ISIGN=1                                                           9178.   
      IF(SAVE.GT.0.0) ISIGN=-1                                          9179.   
      IF(INDEX.LT.4) GO TO 325                                          9180.   
      IEND=INDEX-1                                                      9181.   
      DO 320 K=1,LM                                                     9182.   
      DO 320 I=2,IEND                                                   9183.   
      CHECK=FKEY(I,K)*ISIGN                                             9184.   
  320 IF(CHECK.GT.HS)HS=CHECK                                           9185.   
  325 CONTINUE                                                          9186.   
      IF(INDEX.GT.(JM-2))GO TO 335                                      9187.   
      JSTART=INDEX+1                                                    9188.   
      DO 330 K=1,LM                                                     9189.   
      DO 330 I=JSTART,JM                                                9190.   
      CHECK=FKEY(I,K)*ISIGN                                             9191.   
  330 IF(CHECK.GT.HN)HN=CHECK                                           9192.   
  335 CONTINUE                                                          9193.   
      KEYNR(28,KEYCT)=ABS(HN)+0.5                                       9194.   
      KEYNR(29,KEYCT)=INT(SAVE+10000.5 )-10000                          9195.   
      KEYNR(30,KEYCT)=ABS(HS)+0.5                                       9196.   
      RETURN                                                            9197.   
C****                                                                   9198.   
      ENTRY KEYD2K (ASUM)                                               9199.   
C**** EDDY KINETIC ENERGY                                               9200.   
      KEYNR(22,KEYCT)=INT(ASUM(JEQ)+.5)                                 9201.   
      BIG=-99999.                                                       9202.   
      I35=2.+JMM1*125./180.                                             9203.   
      I70=2.+JMM1*160./180.                                             9204.   
      DO 440 I=I35,I70                                                  9205.   
      IF(ASUM(I).LT.BIG) GO TO 440                                      9206.   
      BIG=ASUM(I)                                                       9207.   
  440 CONTINUE                                                          9208.   
      KEYNR(23,KEYCT)=INT(BIG+.5)                                       9209.   
      RETURN                                                            9210.   
C****                                                                   9211.   
      ENTRY KEYD2N (NT,ASUM,SUMFAC)                                     9212.   
C**** NORTHWARD TRANSPORTS                                              9213.   
  500 BIG=-99999.                                                       9214.   
      JEQP1=JEQ+1                                                       9215.   
      DO 510 I=JEQP1,JM                                                 9216.   
      IF(ASUM(I).LT.BIG) GO TO 510                                      9217.   
      BIG=ASUM(I)                                                       9218.   
      INDEX=I                                                           9219.   
  510 CONTINUE                                                          9220.   
      BIG=BIG*SUMFAC                                                    9221.   
      NTDIF=NT-21                                                       9222.   
      GO TO (392,392,392,390,390,396,394,390,390,400,400,398),NTDIF     9223.   
  390 CONTINUE                                                          9224.   
  392 KEYNR(NT+11,KEYCT)=INT(BIG+.5)                                    9225.   
      RETURN                                                            9226.   
  394 KEYNR(38,KEYCT)=JLAT(INDEX,2)                                     9227.   
  396 KEYNR(NT+9,KEYCT)=INT(BIG+.5)                                     9228.   
      RETURN                                                            9229.   
  398 KEYNR(42,KEYCT)=JLAT(INDEX,2)                                     9230.   
  400 KEYNR(NT+8,KEYCT)=INT(BIG+.5)                                     9231.   
      RETURN                                                            9232.   
C****                                                                   9233.   
C**** ENTRY CALLED FROM DIAG4                                           9234.   
C****                                                                   9235.   
      ENTRY KEYD4 (IK)                                                  9236.   
      KEYNR(16,KEYCT)=(IK(4)+IK(14)+5)/10                               9237.   
      KEYNR(17,KEYCT)=(IK(3)+IK(13)+5)/10                               9238.   
      KEYNR(18,KEYCT)=(IK(2)+IK(12)+5)/10                               9239.   
      KEYNR(19,KEYCT)=(IK(1)+IK(11)+5)/10                               9240.   
      KEYNR(20,KEYCT)=(IK(10)+IK(20)+5)/10                              9241.   
      KEYNR(21,KEYCT)=(IK(8)+IK(18)+5)/10                               9242.   
      RETURN                                                            9243.   
C****                                                                   9244.   
C**** ENTRY CALLED FROM DIAG5                                           9245.   
      ENTRY KEYD5A(NM1,FK1,FK2)                                         9246.   
C**** SPECTRAL ANALYSIS                                                 9247.   
C**** CALCULATES THE 'SLOPE' OF THE MEAN KINETIC ENERGY FOR THE TROPO-  9248.   
C**** SPHERE AT THE EQUATOR AND AT 45 DEGREES NORTH. SLOPE IS DEFINED   9249.   
C**** AS (LNR/LN2)*10 WHERE R IS THE RATIO  OF THE AVERAGE KE IN WAVE   9250.   
C**** NUMBERS 7, 8 AND 9 TO WAVE NUMBERS 11, 12 AND 13.                 9251.   
      NM1M6=NM1-6                                                       9252.   
      GO TO (601,602,602,600,605,606,607),NM1M6                         9253.   
  600 RETURN                                                            9254.   
  601 FEQLO=FK1                                                         9255.   
      F45LO=FK2                                                         9256.   
      RETURN                                                            9257.   
  602 FEQLO=FEQLO+FK1                                                   9258.   
      F45LO=F45LO+FK2                                                   9259.   
      RETURN                                                            9260.   
  605 FEQHI=FK1                                                         9261.   
      F45HI=FK2                                                         9262.   
      RETURN                                                            9263.   
  606 FEQHI=FEQHI+FK1                                                   9264.   
      F45HI=F45HI+FK2                                                   9265.   
      RETURN                                                            9266.   
  607 FEQHI=FEQHI+FK1+1.E-20                                            9267.   
      REQ=FEQLO/FEQHI                                                   9268.   
      KEYNR(31,KEYCT)=10.*DLOG(REQ+1.E-20)/DLOG(1.5)+.5                 9269.   
      F45HI=F45HI+FK2+1.E-20                                            9270.   
      R45=F45LO/F45HI                                                   9271.   
      KEYNR(32,KEYCT)=10.*DLOG(R45+1.E-20)/DLOG(1.5)+.5                 9272.   
      RETURN                                                            9273.   
C****                                                                   9274.   
      ENTRY DIAGKN                                                      9275.   
C**** PRINTS THE TABLE OF KEY NUMBERS                                   9276.   
C****                                                                   9277.   
      IHOUR0=TOFDY0+.5                                                  9278.   
      IHOUR=TOFDAY+.5                                                   9279.   
      TAUDIF=TAU-TAU0                                                   9280.   
      if(KEYCT.gt.50)then
         print *,'1 KEYCT=',KEYCT
         stop
      endif
      CKEYNR(1,KEYCT)=JMNTH0                                            9281.   
      IF(KEYCT.EQ.1) CKEYNR(1,KEYCT)=IC                                 9282.   
      IF(KEYCT.GE.2.AND.CKEYNR(1,KEYCT-1).EQ.JMNTH0) KEYCT=KEYCT-1      9283.   
      if(JYEAR0.ne.JYEAR)then
       JYRPR=JYEAR-KEYCT/12
      WRITE(6,901) XLABEL                                               9284.   
      WRITE(6,910) IDAY0,IHOUR0,JDATE0,JMNTH0,JYEAR0,TAU0,IDAY,IHOUR,   9285.   
     *  JDATE,JMONTH,JYEAR,TAU,TAUDIF                                   9286.   
      WRITE(6,902)                                                      9287.   
      DO 810 I=1,KEYCT                                                  9288.   
c     IF(CKEYNR(1,I).EQ.JAN) WRITE (6,905)                              9289.   
      IF(CKEYNR(1,I).EQ.JAN) THEN
c      print *,JYRPR
       JYRPR=JYRPR+1
      ENDIF
  810 WRITE (6,905) (KEYNR(K,I),K=1,42)                                 9290.   
      WRITE (6,915)                                                     9291.   
      end if
         DO 815 K=1,42                                                  9292.   
  815    KEYDS(K)=KEYNR(K,KEYCT)                                        9293.   
      KEYCT=KEYCT+1                                                     9294.   
      KEYMAX=49                                                         9295.   
      IF(CKEYNR(1,1).NE.IC) KEYMAX=48                                   9296.   
      IF(KEYCT.LE.KEYMAX) RETURN                                        9297.   
C**** ROLL UP KEY NUMBERS 1 YEAR AT A TIME                              9298.   
      DO 820 K=1,36                                                     9299.   
      DO 820 I=1,42                                                     9300.   
  820 KEYNR(I,K)=KEYNR(I,K+KEYMAX-36)                                   9301.   
      DO 880 K=37,50                                                    9302.   
  880 KEYNR(2,K)=0                                                      9303.   
      KEYCT=37                                                          9304.   
      RETURN                                                            9305.   
  901 FORMAT('1',33A4)                                                  9306.   
  902 FORMAT ('0',7X,'NH NH AL AB NT NT PR        T   T-OF-ATM  EKE   ZK9307.   
     *E           EKE   JET-STREAMS STREAM-FN  EKE  NOR-TRAN NOR-TRAN NO9308.   
     *RTH-TRANS'/                                                       9309.   
     *         5X,'CL SN OI BE BY RD HT EC SN LAT OF  GL  GRAD ----- ---9310.   
     *-- EPE ZPE ------ NORTH SOUTH --------- SLOPE DRY-STAT STAT-ENR AN9311.   
     *G MOMENTM'/                                                       9312.   
     *         5X,'CV CV CV DO AT P0 Z0 IP HT  HT GD  OB NH SH NH SH NH 9313.   
     *SH  NH  NH EQ  ML VL LT VL LT NH MAX SH EQ 45 SE ED TL ED TL LT SE9314.   
     * ED TL LT'/)                                                      9315.   
  905 FORMAT (1X,A3,2I3,I2,I4,5I3,I4,I3,I4,6I3,2I4,I3,I4,5I3,I4,13I3)   9316.   
  910 FORMAT ('0',15X,'DAY',I5,', HR',I3,' (',I2,A5,I5,')',F8.1,        9317.   
     *  '   TO   DAY',I5,', HR',I3,' (',I2,A5,I5,')',F8.1,'   DIF',     9318.   
     *  F6.1,' HR',7X,I5,I5)                                            9319.   
  915 FORMAT('0')                                                       9320.   
      END                                                               9321.   
