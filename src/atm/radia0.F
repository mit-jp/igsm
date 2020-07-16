
#include "ctrparam.h"

      SUBROUTINE RADIA0 (IM,JM,CO2,READGHG)                             4501.   
C****                                                                   4502.   
C**** THIS SUBROUTINE SETS THE RADIATION CONTROL PARAMETERS AND         4503.   
C**** CALCULATES AREA WEIGHTED LATITUDES FOR A STANDARD GRID ETC        4504.   
C****                                                                   4505.   

#include "chem_para"
#include "chem_com"
#if ( defined CLM )
#include "CLM.h"
#endif

      REAL LT1,LT2                                                      4506.   
      DIMENSION COSZ(IM,JM),COSZA(IM,JM)                                4507.   
      DIMENSION SINJ(46),COSJ(46),RI(72),SINI(72),COSI(72)              4508.   
      COMMON/WORK5/LT1(72),LT2(72),SLT1(72),SLT2(72),S2LT1(72),S2LT2(72)4509.   
     *  ,DEGLAT(46),DEGLON(72)                                          4510.   
      COMMON/SCJL/SINJ,COSJ
      COMMON/CO2TRND/ALFFOR,CO2TR,YEARGT,CO2IN,INYRAD
     & ,YEARGG,YEARS0,YEARSA,YEARVA
     & ,KTREND
C                                                                       4511.   
C                   RADCOM:      CONTROL/INPUT PARAMETERS               4512.   
C                                                                       4513.   
      COMMON/RADCOM/VADATA(11,4,3),DGLAT(46),DGLON(72),TMINSR,FULGAS(18)4514.   
     A             ,FRACSL,RATQSL,FOGTSL,PTLISO,TLGRAD,TKCICE,FGOLDU(18)4514.5  
     B             ,FLONO3,FRAYLE,FCLDTR,FCLDSR,FALGAE,FMARCL,FEMTRA(6) 4515.   
     C             ,WETTRA,WETSRA,DMOICE,DMLICE,LICETK,NTRCE,FZASRA(6)  4515.5  
     D             ,ID5(5),ITR(4),IMG(2),ILG(2),LAPGAS,KWVCON,NORMS0,NV 4516.   
     E             ,KEEPRH,KEEPAL,ISOSCT,IHGSCT,KFRACC,KGASSR,KAERSR    4516.5  
     F             ,MARCLD,LAYTOP,LMR,LMRP,JMLAT,IMLON,KFORCE,LASTVC    4517.   
C                                                                       4517.5  
C                                BASIC RADCOM INPUT DATA                4518.   
C                                                                       4518.5  
     G             ,PLE(40),HLB(40),TLB(40),TLT(40),TL(40),U0GAS(40,9)  4519.   
     H             ,ULGAS(40,9),TRACER(40,4),RTAU(40),QL(40),RHL(40)    4519.5  
     I             ,POCEAN,PEARTH,POICE,PLICE,AGESN,SNOWE,SNOWOI,SNOWLI 4520.   
     J             ,TGO,TGE,TGOI,TGLI,TS,WS,WEARTH,ZOICE,FSPARE(200)    4520.5  
     K             ,S0,COSZN,PVT(11),BXA(153),SRBXAL(15,2),FRC(5),LUXGAS4520.6  
     L             ,JYEAR,JDAY,JLAT,ILON,MEANAL,KALVIS,ISPARE(25),SGPS  4520.8  
C                                                                       5033.
C                                BASIC RADCOM OUTPUT DATA               5034.
C                                                                       5035.
c    M             ,TRDFLB(40),TRUFLB(40),TRNFLB(40),TRFCRL(40),TRSLCR  5036.
c    N             ,SRDFLB(40),SRUFLB(40),SRNFLB(40),SRFHRL(40),SRSLHR  5037.
c    O             ,SRIVIS,SROVIS,PLAVIS,SRINIR,SRONIR,PLANIR,SRXATM(4) 5038.
c    P             ,SRDVIS,SRUVIS,ALBVIS,SRDNIR,SRUNIR,ALBNIR,FSRNFG(4) 5039.
c    Q             ,SRTVIS,SRRVIS,SRAVIS,SRTNIR,SRRNIR,SRANIR,FTRUFG(4) 5040.
c    R             ,TRDFGW,TRUFGW,TRUFTW,BTEMPW,TRDFSL,TRUFSL,DTRUFG(4) 5041.
c    S             ,TRSLTS,TRSLTG,TRSLWV,TRSLBS,TTRUFG,LBOTCL,LTOPCL    5042.
!     DATA TWOPI/6.283185/,ZERO1/1.E-5/,ZERO2/5.E-3/                    4521.   
!     03/07/06
      DATA TWOPI/6.283185/,ZERO1/1.E-2/,ZERO2/5.E-3/                    4521.   
C**** COMPUTE THE AREA WEIGHTED LATITUDES AND THEIR SINES AND COSINES   4522.   
      JMM1=JM-1                                                         4523.   
      PHIS=-.25*TWOPI                                                   4524.   
      SPHIS=-1.                                                         4525.   
      CPHIS=0.                                                          4526.   
      DO 20 J=1,JMM1                                                    4527.   
      PHIN=(TWOPI/(JMM1+JMM1))*(J-.5*JM)                                4528.   
      SPHIN=SIN(PHIN)                                                   4529.   
      CPHIN=COS(PHIN)                                                   4530.   
      PHIM=(PHIN*SPHIN+CPHIN-PHIS*SPHIS-CPHIS)/(SPHIN-SPHIS)            4531.   
      DEGLAT(J)=(360./TWOPI)*PHIM                                       4532.   
      SINJ(J)=SIN(PHIM)                                                 4533.   
      COSJ(J)=COS(PHIM)                                                 4534.   
      PHIS=PHIN                                                         4535.   
      SPHIS=SPHIN                                                       4536.   
   20 CPHIS=CPHIN                                                       4537.   
      PHIN=.25*TWOPI                                                    4538.   
      SPHIN=1.                                                          4539.   
      CPHIN=0.                                                          4540.   
      PHIM=(PHIN*SPHIN+CPHIN-PHIS*SPHIS-CPHIS)/(SPHIN-SPHIS)            4541.   
      DEGLAT(JM)=(360./TWOPI)*PHIM                                      4542.   
      SINJ(JM)=SIN(PHIM)                                                4543.   
      COSJ(JM)=COS(PHIM)                                                4544.   
#if ( defined CLM35 )
!APS
      DO 21 J=1,JM
      PHIM=-TWOPI*0.25+TWOPI/360.*(J-1)*180./JMM1
!     print *,j,PHIM,PHIM*360./TWOPI
      DEGLAT(J)=(360./TWOPI)*PHIM                                       4532.   
      SINJ(J)=SIN(PHIM)                                                 4533.   
      COSJ(J)=COS(PHIM)                                                 4534.   
   21 CONTINUE
#endif
C**** COMPUTE THE SINES AND COSINES OF LONGITUDE                        4545.   
      DO 40 I=1,IM                                                      4546.   
      RI(I)=(TWOPI/IM)*(I-1)                                            4547.   
      DEGLON(I)=(360./IM)*(I-1)                                         4548.   
      SINI(I)=SIN(RI(I))                                                4549.   
   40 COSI(I)=COS(RI(I))                                                4550.   
C**** MODIFY AND PRINT OUT THE RADIATION CONTROL PARAMETERS             4552.   
      IF (CO2.GT.0.) FULGAS(2)=CO2                                      4553.1  
      DMOICE=10.                                                        4553.11 
C  Convert masking depth over land and ocean ice to meters of water
c     DMOICE=0.01
c     DMLICE=0.01
C  Convert masking depth over land and ocean ice to meters of water
C     FOR THE 'REAL' TRANSIENT EXPERIMENT USE CO2=-FLOAT(KTREND)-AVGAER 4553.15 
C     TO USE THE DEFAULTS EXCEPT FOR CO2 SET CO2.GT.0.                  4553.2  
!     KTREND=-CO2                                                       4553.25 
       print *,' KTREND=',KTREND

!#ifdef PREDICTED_GASES
#if ( defined PREDICTED_GASES  && !defined FIXED_FOR )
        LAPGAS=0
        IMG(2)=9
       print *,' RADIA0 Chem ',1,READGHG,LAPGAS
#endif

#if ( !defined OLD_LAPGAS )
        LAPGAS=0
        IMG(2)=9
#endif
      IF (CO2.GT.0.) GO TO 50                                           4553.3  
      NTRCE=1                                                           4553.35 
      ITR(1)=1                                                          4553.4  
      IF (KTREND.EQ.0) KFORCE=26789                                     4553.45 
      FGOLDU(6)=1.                                                      4553.5  
#if ( defined OLD_LAPGAS )
      LAPGAS=2                                                          4553.55 
#endif
C     FULGAS(8)=0.                                                      4553.6  
C     FULGAS(9)=0.                                                      4553.65 
C     FGOLDU(1)=.005/.012                                               4553.7  
c  50 CALL RCOMP1 (21,0,60)                                             4553.75 
   50 continue
#if ( defined IPCC_FORCING )
      print *,'Run with IPCC FORCING'
#else
#if ( defined ETMINAN_FORCING )
      print *,'Run with ETMINAN FORCING'
#else
      print *,'Run with GISS FORCING'
#endif
#endif
          print *,' RADIA0 JYEAR=',JYEAR,' JDAY=',JDAY
c     IF (KTREND.GT.0) CALL FORSET(1958.,KTREND,1)                      4553.8  
      IF (KTREND.GT.0)then
          TNOW=JYEAR+(JDAY-.5)/365. 
          if(KTREND.GT.20) TNOW=1765.
          print *,' RADIA0 JYEAR=',JYEAR,' JDAY=',JDAY
          if (KTREND.EQ.5)then
!          print *,' TNOW=',TNOW
           print *,' TREF=',TNOW
           CALL FORSET(TNOW,KTREND,1)
          else
!          04/18/2006
           print *,' TREF=',YEARGT
          CALL FORSET(YEARGT,KTREND,1)
         endif
      endif
c   06/20/2005
c        CALL RCOMP1 (521,0,60)
         print *,'Before CALL RCOMP1'
         CALL RCOMP1 (521,0,60,KTREND)
         print *,'After CALL RCOMP1'
c        if(READGHG.lt.0.5) CALL WRITER (1,0)                           4554.   
      RETURN                                                            4555.   
C****                                                                   4556.   
C****                                                                   4557.   
      ENTRY COSZT (IM,JM,SIND,COSD,ROT1,ROT2,COSZ)                      4558.   
C****                                                                   4559.   
C**** THIS ENTRY COMPUTES THE ZENITH ANGLE WEIGHTED BY DAYTIME          4560.   
C**** HOURS FROM ROT1 TO ROT2, GREENWICH MEAN TIME IN RADIANS.  ROT1    4561.   
C**** MUST BE BETWEEN 0 AND 2*PI.  ROT2 MUST BE BETWEEN ROT1 AND        4562.   
C**** ROT1+2*PI.  I=1 MUST LIE ON THE INTERNATIONAL DATE LINE.          4563.   
C****                                                                   4564.   
#if ( defined CLM35 )
! FOR CLM35
      DAYFRAC=ROT1/TWOPI
      SINDD = sin(declin_clm)
      COSD = DSQRT(1.-SINDD*SINDD)                                      8283.
      SIND = SINDD                                                      8284.
!     print *,'FROM COSZT'
!     print *,'declin_clm=',declin_clm
! FOR CLM35
#endif
      DROT=ROT2-ROT1                                                    4565.   
C**** COMPUTE THE SINES AND COSINES OF THE INITIAL AND FINAL GMT'S      4566.   
  100 SR1=SIN(ROT1)                                                     4567.   
      CR1=COS(ROT1)                                                     4568.   
      SR2=SIN(ROT2)                                                     4569.   
      CR2=COS(ROT2)                                                     4570.   
C**** COMPUTE THE INITIAL AND FINAL LOCAL TIMES (MEASURED FROM NOON TO  4571.   
C****   NOON) AND THEIR SINES AND COSINES                               4572.   
      DO 120 I=1,IM                                                     4573.   
      LT1(I)=ROT1+RI(I)                                                 4574.   
      SLT1(I)=SR1*COSI(I)+CR1*SINI(I)                                   4575.   
      LT2(I)=ROT2+RI(I)                                                 4576.   
  120 SLT2(I)=SR2*COSI(I)+CR2*SINI(I)                                   4577.   
C****                                                                   4578.   
C**** CALCULATION FOR POLAR GRID BOXES                                  4579.   
C****                                                                   4580.   
      DO 200 J=1,JM,JMM1                                                4581.   
!     print '(i3,4e14.6)',J,SINJ(J),SIND,COSJ(J),COSD  
!     print '(2e14.6)',ASIN(SINJ(J)),180.*ASIN(SINJ(J))/3.1415
      SJSD=SINJ(J)*SIND                                                 4582.   
      CJCD=COSJ(J)*COSD                                                 4583.   
#if ( defined CLM35 )
      COSZ(1,J)=max(0.0,
     & SJSD-CJCD*cos((float(JDAY)+DAYFRAC)*TWOPI+0.5*TWOPI))
!     print '(4e14.6)',DAYFRAC,float(JDAY)+DAYFRAC,
!    &   cos((float(JDAY)+DAYFRAC)*TWOPI+0.5*TWOPI),COSZ(1,J)
#else
      IF(SJSD+CJCD.LE.ZERO1) GO TO 180                                  4584.   
      IF(SJSD-CJCD.GE.0.) GO TO 160                                     4585.   
C**** AVERAGE COSZ FROM DAWN TO DUSK NEAR THE POLES                     4586.   
      DUSK=ACOS(-SJSD/CJCD)                                             4587.   
      SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD                              4588.   
      DAWN=-DUSK                                                        4589.   
      SDAWN=-SDUSK                                                      4590.   
      COSZ(1,J)=(SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN))/TWOPI             4591.   
      GO TO 200                                                         4592.   
C**** CONSTANT DAYLIGHT NEAR THE POLES                                  4593.   
  160 COSZ(1,J)=SJSD                                                    4594.   
      GO TO 200                                                         4595.   
C**** CONSTANT NIGHTIME NEAR THE POLES                                  4596.   
  180 COSZ(1,J)=0.                                                      4597.   
#endif
  200 CONTINUE                                                          4598.   
C****                                                                   4599.   
C**** LOOP OVER NON-POLAR LATITUDES                                     4600.   
C****                                                                   4601.   
      DO 500 J=2,JMM1                                                   4602.   
!     print '(i3,4e14.6)',J,SINJ(J),SIND,COSJ(J),COSD  
!     print '(2e14.6)',ASIN(SINJ(J)),180.*ASIN(SINJ(J))/3.1415
      SJSD=SINJ(J)*SIND                                                 4603.   
      CJCD=COSJ(J)*COSD                                                 4604.   
#if ( defined CLM35 )
      COSZ(1,J)=max(0.0,
     & SJSD-CJCD*cos((float(JDAY)+DAYFRAC)*TWOPI+0.5*TWOPI))
#else
      IF(SJSD+CJCD.LE.ZERO1) GO TO 460                                  4605.   
      IF(SJSD-CJCD.GE.0.) GO TO 420                                     4606.   
C**** COMPUTE DAWN AND DUSK (AT LOCAL TIME) AND THEIR SINES             4607.   
      DUSK=ACOS(-SJSD/CJCD)                                             4608.   
      SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD                              4609.   
      DAWN=-DUSK                                                        4610.   
      SDAWN=-SDUSK                                                      4611.   
C**** NEITHER CONSTANT DAYTIME NOR CONSTANT NIGHTIME AT THIS LATITUDE,  4612.   
C**** LOOP OVER LONGITUDES                                              4613.   
!     03/07/06 solar radiation
      ZERO2=ZERO1/CJCD
      DO 400 I=1,IM                                                     4614.   
C**** FORCE DUSK TO LIE BETWEEN LT1 AND LT1+2*PI                        4615.   
      IF(DUSK.GT.LT1(I)+ZERO2) GO TO 220                                4616.   
      DUSK=DUSK+TWOPI                                                   4617.   
      DAWN=DAWN+TWOPI                                                   4618.   
  220 IF(DAWN.LT.LT2(I)-ZERO2) GO TO 240                                4619.   
C**** CONTINUOUS NIGHTIME FROM INITIAL TO FINAL TIME                    4620.   
      COSZ(I,J)=0.                                                      4621.   
      GO TO 400                                                         4622.   
  240 IF(DAWN.GE.LT1(I)) GO TO 300                                      4623.   
      IF(DUSK.LT.LT2(I)) GO TO 260                                      4624.   
C**** CONTINUOUS DAYLIGHT FROM INITIAL TIME TO FINAL TIME               4625.   
      COSZ(I,J)=SJSD+CJCD*(SLT2(I)-SLT1(I))/DROT                        4626.   
      GO TO 400                                                         4627.   
  260 IF(DAWN+TWOPI.LT.LT2(I)-ZERO2) GO TO 280                          4628.   
C**** DAYLIGHT AT INITIAL TIME AND NIGHT AT FINAL TIME                  4629.   
      COSZ(I,J)=(SJSD*(DUSK-LT1(I))+CJCD*(SDUSK-SLT1(I)))/DROT          4630.   
      GO TO 400                                                         4631.   
C**** DAYLIGHT AT INITIAL AND FINAL TIMES WITH NIGHTIME IN BETWEEN      4632.   
  280 COSZ(I,J)=(SJSD*(LT2(I)-DAWN-TWOPI+DUSK-LT1(I))+CJCD*             4633.   
     *  (SLT2(I)-SDAWN+SDUSK-SLT1(I)))/DROT                             4634.   
      GO TO 400                                                         4635.   
  300 IF(DUSK.LT.LT2(I)) GO TO 320                                      4636.   
C**** NIGHT AT INITIAL TIME AND DAYLIGHT AT FINAL TIME                  4637.   
      COSZ(I,J)=(SJSD*(LT2(I)-DAWN)+CJCD*(SLT2(I)-SDAWN))/DROT          4638.   
      GO TO 400                                                         4639.   
C**** NIGHTIME AT INITIAL AND FINAL TIMES WITH DAYLIGHT IN BETWEEN      4640.   
  320 COSZ(I,J)=(SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN))/DROT              4641.   
  400 CONTINUE                                                          4642.   
      GO TO 500                                                         4643.   
C**** CONSTANT DAYLIGHT AT THIS LATITUDE                                4644.   
  420 DO 440 I=1,IM                                                     4645.   
  440 COSZ(I,J)=SJSD+CJCD*(SLT2(I)-SLT1(I))/DROT                        4646.   
      GO TO 500                                                         4647.   
C**** CONSTANT NIGHTIME AT THIS LATITUDE                                4648.   
  460 DO 480 I=1,IM                                                     4649.   
  480 COSZ(I,J)=0.                                                      4650.   
#endif
!     print '(4e14.6)',DAYFRAC,float(JDAY)+DAYFRAC,
!    &   cos((float(JDAY)+DAYFRAC)*TWOPI+0.5*TWOPI),COSZ(1,J)
  500 CONTINUE                                                          4651.   
      RETURN                                                            4652.   
C****                                                                   4653.   
C****                                                                   4654.   
      ENTRY COSZS (IM,JM,SIND,COSD,ROT1,ROT2,COSZ,COSZA)                4655.   
C****                                                                   4656.   
C**** THIS ENTRY COMPUTES THE ZENITH ANGLE TWICE, FIRST WEIGHTED BY THE 4657.   
C**** DAYTIME HOURS FROM ROT1 TO ROT2 AND SECONDLY WEIGHTED BY THE      4658.   
C**** INCIDENT SUN LIGHT FROM ROT1 TO ROT2.  COSZT MUST HAVE BEEN       4659.   
C**** CALLED JUST PREVIOUSLY.                                           4660.   
C****                                                                   4661.   
      DROT=ROT2-ROT1                                                    4662.   
C**** COMPUTE THE SINES AND COSINES OF THE INITIAL AND FINAL GMT'S      4663.   
      SR1=SIN(ROT1)                                                     4664.   
      CR1=COS(ROT1)                                                     4665.   
      SR2=SIN(ROT2)                                                     4666.   
      CR2=COS(ROT2)                                                     4667.   
C**** COMPUTE THE INITIAL AND FINAL LOCAL TIMES (MEASURED FROM NOON TO  4668.   
C****   NOON) AND THEIR SINES AND COSINES                               4669.   
      DO 520 I=1,IM                                                     4670.   
      LT1(I)=ROT1+RI(I)                                                 4671.   
      SLT1(I)=SR1*COSI(I)+CR1*SINI(I)                                   4672.   
      CLT1=CR1*COSI(I)-SR1*SINI(I)                                      4673.   
      S2LT1(I)=2.*SLT1(I)*CLT1                                          4674.   
      LT2(I)=ROT2+RI(I)                                                 4675.   
      SLT2(I)=SR2*COSI(I)+CR2*SINI(I)                                   4676.   
      CLT2=CR2*COSI(I)-SR2*SINI(I)                                      4677.   
  520 S2LT2(I)=2.*SLT2(I)*CLT2                                          4678.   
C****                                                                   4679.   
C**** CALCULATION FOR POLAR GRID BOXES                                  4680.   
C****                                                                   4681.   
      DO 600 J=1,JM,JMM1                                                4682.   
      SJSD=SINJ(J)*SIND                                                 4683.   
      CJCD=COSJ(J)*COSD                                                 4684.   
      IF(SJSD+CJCD.LE.ZERO1) GO TO 580                                  4685.   
      IF(SJSD-CJCD.GE.0.) GO TO 560                                     4686.   
C**** AVERAGE COSZ FROM DAWN TO DUSK NEAR THE POLES                     4687.   
      CDUSK=-SJSD/CJCD                                                  4688.   
      DUSK=ACOS(CDUSK)                                                  4689.   
      SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD                              4690.   
      S2DUSK=2.*SDUSK*CDUSK                                             4691.   
      DAWN=-DUSK                                                        4692.   
      SDAWN=-SDUSK                                                      4693.   
      S2DAWN=-S2DUSK                                                    4694.   
      ECOSZ=SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN)                         4695.   
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SDUSK-SDAWN)+                       4696.   
     *  .5*CJCD*(DUSK-DAWN+.5*(S2DUSK-S2DAWN)))                         4697.   
      COSZ(1,J)=ECOSZ/TWOPI                                             4698.   
      COSZA(1,J)=ECOSQZ/ECOSZ                                           4699.   
      GO TO 600                                                         4700.   
C**** CONSTANT DAYLIGHT NEAR THE POLES                                  4701.   
  560 ECOSZ=SJSD*TWOPI                                                  4702.   
      ECOSQZ=SJSD*ECOSZ+.5*CJCD*CJCD*TWOPI                              4703.   
      COSZ(1,J)=ECOSZ/TWOPI                                             4704.   
      COSZA(1,J)=ECOSQZ/ECOSZ                                           4705.   
      GO TO 600                                                         4706.   
C**** CONSTANT NIGHTIME NEAR THE POLES                                  4707.   
  580 COSZ(1,J)=0.                                                      4708.   
      COSZA(1,J)=0.                                                     4709.   
  600 CONTINUE                                                          4710.   
C****                                                                   4711.   
C**** LOOP OVER NON-POLAR LATITUDES                                     4712.   
C****                                                                   4713.   
      DO 900 J=2,JMM1                                                   4714.   
      SJSD=SINJ(J)*SIND                                                 4715.   
      CJCD=COSJ(J)*COSD                                                 4716.   
      IF(SJSD+CJCD.LE.ZERO1) GO TO 860                                  4717.   
      IF(SJSD-CJCD.GE.0.) GO TO 820                                     4718.   
C**** COMPUTE DAWN AND DUSK (AT LOCAL TIME) AND THEIR SINES             4719.   
      CDUSK=-SJSD/CJCD                                                  4720.   
      DUSK=ACOS(CDUSK)                                                  4721.   
      SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD                              4722.   
      S2DUSK=2.*SDUSK*CDUSK                                             4723.   
      DAWN=-DUSK                                                        4724.   
      SDAWN=-SDUSK                                                      4725.   
      S2DAWN=-S2DUSK                                                    4726.   
C**** NEITHER CONSTANT DAYTIME NOR CONSTANT NIGHTIME AT THIS LATITUDE,  4727.   
C**** LOOP OVER LONGITUDES                                              4728.   
!     03/07/06 solar radiation
      ZERO2=ZERO1/CJCD
      DO 800 I=1,IM                                                     4729.   
C**** FORCE DUSK TO LIE BETWEEN LT1 AND LT1+2*PI                        4730.   
      IF(DUSK.GT.LT1(I)+ZERO2) GO TO 620                                4731.   
      DUSK=DUSK+TWOPI                                                   4732.   
      DAWN=DAWN+TWOPI                                                   4733.   
  620 IF(DAWN.LT.LT2(I)-ZERO2) GO TO 640                                4734.   
C**** CONTINUOUS NIGHTIME FROM INITIAL TO FINAL TIME                    4735.   
      COSZ(I,J)=0.                                                      4736.   
      COSZA(I,J)=0.                                                     4737.   
      GO TO 800                                                         4738.   
  640 IF(DAWN.GE.LT1(I)) GO TO 700                                      4739.   
      IF(DUSK.LT.LT2(I)) GO TO 660                                      4740.   
C**** CONTINUOUS DAYLIGHT FROM INITIAL TIME TO FINAL TIME               4741.   
      ECOSZ=SJSD*DROT+CJCD*(SLT2(I)-SLT1(I))                            4742.   
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SLT2(I)-SLT1(I))+                   4743.   
     *  .5*CJCD*(DROT+.5*(S2LT2(I)-S2LT1(I))))                          4744.   
      COSZ(I,J)=ECOSZ/DROT                                              4745.   
      COSZA(I,J)=ECOSQZ/ECOSZ                                           4746.   
      GO TO 800                                                         4747.   
  660 IF(DAWN+TWOPI.LT.LT2(I)-ZERO2) GO TO 680                          4748.   
C**** DAYLIGHT AT INITIAL TIME AND NIGHT AT FINAL TIME                  4749.   
      ECOSZ=SJSD*(DUSK-LT1(I))+CJCD*(SDUSK-SLT1(I))                     4750.   
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SDUSK-SLT1(I))+                     4751.   
     *  .5*CJCD*(DUSK-LT1(I)+.5*(S2DUSK-S2LT1(I))))                     4752.   
      COSZ(I,J)=ECOSZ/DROT                                              4753.   
      COSZA(I,J)=ECOSQZ/ECOSZ                                           4754.   
      GO TO 800                                                         4755.   
C**** DAYLIGHT AT INITIAL AND FINAL TIMES WITH NIGHTIME IN BETWEEN      4756.   
  680 ECOSZ=SJSD*(DROT-DAWN-TWOPI+DUSK)+                                4757.   
     *  CJCD*(SLT2(I)-SDAWN+SDUSK-SLT1(I))                              4758.   
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SDUSK-SLT1(I)+SLT2(I)-SDAWN)+       4759.   
     *  .5*CJCD*(DUSK+DROT-DAWN-TWOPI+                                  4760.   
     *  .5*(S2DUSK-S2LT1(I)+S2LT2(I)-S2DAWN)))                          4761.   
      COSZ(I,J)=ECOSZ/DROT                                              4762.   
      COSZA(I,J)=ECOSQZ/ECOSZ                                           4763.   
      GO TO 800                                                         4764.   
  700 IF(DUSK.LT.LT2(I)) GO TO 720                                      4765.   
C**** NIGHT AT INITIAL TIME AND DAYLIGHT AT FINAL TIME                  4766.   
      ECOSZ=SJSD*(LT2(I)-DAWN)+CJCD*(SLT2(I)-SDAWN)                     4767.   
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SLT2(I)-SDAWN)+                     4768.   
     *  .5*CJCD*(LT2(I)-DAWN+.5*(S2LT2(I)-S2DAWN)))                     4769.   
      COSZ(I,J)=ECOSZ/DROT                                              4770.   
      COSZA(I,J)=ECOSQZ/ECOSZ                                           4771.   
      GO TO 800                                                         4772.   
C**** NIGHTIME AT INITIAL AND FINAL TIMES WITH DAYLIGHT IN BETWEEN      4773.   
  720 ECOSZ=SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN)                         4774.   
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SDUSK-SDAWN)+                       4775.   
     *  .5*CJCD*(DUSK-DAWN+.5*(S2DUSK-S2DAWN)))                         4776.   
      COSZ(I,J)=ECOSZ/DROT                                              4777.   
      COSZA(I,J)=ECOSQZ/ECOSZ                                           4778.   
  800 CONTINUE                                                          4779.   
      GO TO 900                                                         4780.   
C**** CONSTANT DAYLIGHT AT THIS LATITUDE                                4781.   
  820 DO 840 I=1,IM                                                     4782.   
      ECOSZ=SJSD*DROT+CJCD*(SLT2(I)-SLT1(I))                            4783.   
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SLT2(I)-SLT1(I))+                   4784.   
     *  .5*CJCD*(DROT+.5*(S2LT2(I)-S2LT1(I))))                          4785.   
      COSZ(I,J)=ECOSZ/DROT                                              4786.   
  840 COSZA(I,J)=ECOSQZ/ECOSZ                                           4787.   
      GO TO 900                                                         4788.   
C**** CONSTANT NIGHTIME AT THIS LATITUDE                                4789.   
  860 DO 880 I=1,IM                                                     4790.   
      COSZ(I,J)=0.                                                      4791.   
  880 COSZA(I,J)=0.                                                     4792.   
  900 CONTINUE                                                          4793.   
      RETURN                                                            4794.   
      END                                                               4795.   
