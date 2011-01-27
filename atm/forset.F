
#include "ctrparam.h"

      SUBROUTINE FORSET(TREF,KTREND,KWRITE)                             9101.   
C                                                                       9101.1  
C                   RADCOM:      CONTROL/INPUT PARAMETERS               9101.2  
C                                                                       9101.3  
      COMMON/RADCOM/VADATA(11, 4, 3),DLAT(46),DLON(72),TAUMIN,FULGAS(18)9101.4  
     A             ,FRACSL,RATQSL,FOGTSL,PTLISO,TLGRAD,TKCICE,FGOLDH(18)9101.5  
     B             ,FLONO3,FRAYLE,FCLDTR,FCLDSR,FALGAE,FMARCL,FEMTRA(6) 9101.6  
     C             ,WETTRA,WETSRA,DMOICE,DMLICE,LICETK,NTRACE,FZASRA(6) 9101.7  
     D             ,ID5(5),ITR(4),IMG(2),ILG(2),LAPGAS,KWVCON,NORMS0,NV 9101.8  
     E             ,KEEPRH,KEEPAL,ISOSCT,IHGSCT,KGASSR,KAERSR,KFRACC    9101.9  
     F             ,MARCLD,LAYRAD,NL,NLP,JMLAT ,IMLON ,KFORCE,LASTVC    9102.   
C                                                                       9102.1  
C                                BASIC RADCOM INPUT DATA                9102.2  
C                                                                       9102.3  
     G             ,PLB(40),HLB(40),TLB(40),TLT(40),TLM(40),U0GAS(40,9) 9102.4  
     H             ,ULGAS(40,9),TRACER(40,4),CLDTAU(40),SHL(40),RHL(40) 9102.5  
     I             ,POCEAN,PEARTH,POICE,PLICE,AGESN,SNOWE,SNOWOI,SNOWLI 9102.6  
     J             ,TGO,TGE,TGOI,TGLI,TSL,WMAG,WEARTH,ZOICE,FSPARE(200) 9102.7  
     K             ,S0,COSZ,PVT(11),BXA(153),SRBXAL(15,2),FRC(5),LUXGAS 9102.8  
     L             ,JYEAR,JDAY,JLAT,ILON,MEANAL,KALVIS,ISPARE(25),PSIG0 9102.9  
C                                                                       9103.   
C                                BASIC RADCOM OUTPUT DATA               9103.1  
C                                                                       9103.2  
     M             ,TRDFLB(40),TRUFLB(40),TRNFLB(40),TRFCRL(40),TRSLCR  9103.3  
     N             ,SRDFLB(40),SRUFLB(40),SRNFLB(40),SRFHRL(40),SRSLHR  9103.4  
     O             ,SRIVIS,SROVIS,PLAVIS,SRINIR,SRONIR,PLANIR,SRXATM(4) 9103.5  
     P             ,SRDVIS,SRUVIS,ALBVIS,SRDNIR,SRUNIR,ALBNIR,FSRNFG(4) 9103.6  
     Q             ,SRTVIS,SRRVIS,SRAVIS,SRTNIR,SRRNIR,SRANIR,FTRUFG(4) 9103.7  
     R             ,TRDFGW,TRUFGW,TRUFTW,BTEMPW,TRDFSL,TRUFSL,DTRUFG(4) 9103.8  
     S             ,TRSLTS,TRSLTG,TRSLWV,TRSLBS,TTRUFG,LBOTCL,LTOPCL    9103.9  
C                                                                       9104.   
C                                BLOCKD INITIALIZED DEFAULT DATA        9104.1  
C                                                                       9104.2  
      COMMON/BLOCKD/AGOLDH(11, 5),BGOLDH(11, 5),CGOLDH(11, 5)           9104.3  
     T             ,TRAQEX(25,11),TRAQSC(25,11),TRACOS(25,11)           9104.4  
     T             ,TRCQEX(25, 2),TRCQSC(25, 2),TRCCOS(25, 2)           9104.5  
     S             ,SRAQEX( 6,11),SRAQSC( 6,11),SRACOS( 6,11)           9104.6  
     S             ,SRCQEX( 6, 2),SRCQSC( 6, 2),SRCCOS( 6, 2)           9104.7  
     X             ,AOCEAN(25   ),AGSIDV(25, 4),CLDALB(25, 2)           9104.8  
     Y             ,CMANO2(42   ),TRACEG(25,16),PPMV58(9),Z0(9),ZH(9)   9104.9  
     Z             ,ASNALB(15),AOIALB(15),ALIALB(15),NAERO,NGOLDH,NKSR  9105.   
C                                                                       9105.1  
      DIMENSION XNOW(5),XREF(5),XDT0(5),XDAT(5),XRAT(5),KFOR(5)         9105.2  
C                                                                       9105.3  
      parameter(nlat=N_LAT)
      common /ATCO2/atm_co2(nlat)
      LOGICAL wr25,TRANSR,CONTRR,OBSFOR
      common/wrcom/wr25,TRANSR,CONTRR,OBSFOR
      NGAS=5                                                            9105.4  
      print *,'From FORSET TREF=',TREF,' KTREND=',KTREND
      IF(KTREND.EQ.1) CALL ATREND(XREF,TREF,NGAS)                       9105.5  
      IF(KTREND.EQ.2) CALL BTREND(XREF,TREF,NGAS)                       9105.6  
      IF(KTREND.EQ.3) CALL CTREND(XREF,TREF,NGAS)                       9105.7  
      IF(KTREND.EQ.5) CALL GTREND(XREF,TREF,NGAS)
      IF(KTREND.EQ.4) CALL BMTRND(XREF,TREF,NGAS)
      IF(KTREND.EQ.6) CALL BMTRNDMG(XREF,TREF,NGAS)
c     IF(KTREND.GE.21.and.KTREND.LE.29) 
c    &    CALL STBTRND(XREF,TREF,NGAS,KTREND)
      DO 100 I=1,NGAS                                                   9105.8  
      IF(XREF(I).LT.1.E-06) XREF(I)=1.E-06                              9105.9  
 100  KFOR(I)=1                                                         9106.   
         print *,'From FORSET KTREND=',KTREND
         print *,'XREF'
         print *,XREF
!     CFC11 and CFC12 are in ppb to be used in fits
!     DTDX1D DTDX3D DXDT3D
      PPMV58(2)=XREF(1)                                                 9106.1  
      PPMV58(6)=XREF(2)                                                 9106.2  
      PPMV58(7)=XREF(3)                                                 9106.3  
      PPMV58(8)=XREF(4)/1000.0                                          9106.4  
      PPMV58(9)=XREF(5)/1000.0                                          9106.5  
!     CFC11 and CFC12 are in ppm 
C                                                                       9106.6  
      IF(KWRITE.NE.1) GO TO 120                                         9106.7  
      DO 110 I=1,NGAS                                                   9106.8  
 110  XDAT(I)=XREF(I)                                                   9106.9  
      IF(KTREND.EQ.1) WRITE(6,6001)                                     9107.   
      IF(KTREND.EQ.2) WRITE(6,6002)                                     9107.1  
      IF(KTREND.EQ.3) WRITE(6,6003)                                     9107.2  
      IF(KTREND.EQ.4) WRITE(6,6004)
      IF(KTREND.EQ.5) then
         print *,'PPMV58 after GTREND'
         print *,PPMV58
      endif
      IF(KTREND.EQ.6) then
         print *,'PPMV58 after BMTRNDMG'
         print *,PPMV58
      endif
      WRITE(6,6100)                                                     9107.3  
 120  CONTINUE                                                          9107.4  
 6001 FORMAT(1H1,5X,'(OUTPUT DATA)   GCM RAD EXPECTED TEMPERATURES'     9107.5  
     +      ,T55,'PRESENT TREND FORSET INPUT DATA TO GCM'               9107.6  
     +      ,T96,'RATE OF CHANGE/YR OF TRACE GAS AMOUNTS')              9107.7  
 6002 FORMAT(1H1,5X,'(OUTPUT DATA)   GCM RAD EXPECTED TEMPERATURES'     9107.8  
     +      ,T55,'REDUCED TREND FORSET INPUT DATA TO GCM'               9107.9  
     +      ,T96,'RATE OF CHANGE/YR OF TRACE GAS AMOUNTS')              9108.   
 6003 FORMAT(1H1,5X,'(OUTPUT DATA)   GCM RAD EXPECTED TEMPERATURES'     9108.1  
     +      ,T55,'CURTAIL TREND FORSET INPUT DATA TO GCM'               9108.2  
     +      ,T96,'RATE OF CHANGE/YR OF TRACE GAS AMOUNTS')              9108.3  
 6004 FORMAT(1H1,5X,' TREND FROM BOX MODEL')
 6100 FORMAT(6X,6('-'),'(* 3-D)',32('-'),3X,38('-'),3X,38('-')          9108.4  
     +      /1X,'YEAR DTSUM  *DTCO2   DTN2O   DTCH4   DTF11   DTF12'    9108.5  
     +         ,         '   PPMCO2  PPMN20  PPMCH4  PPTF11  PPTF12'    9108.6  
     +         ,         '   RATCO2  RATN2O  RATCH4  RATF11  RATF12')   9108.7  
C                                                                       9108.8  
      RETURN                                                            9108.9  
C                                                                       9109.   
C------------------------------                                         9109.1  
      ENTRY FORGET(TNOW,KTREN,KWRITE)  
C------------------------------                                         9109.3  
C                                                                       9109.4  
      IF(KTREN.EQ.1) CALL ATREND(XNOW,TNOW,NGAS)                        9109.5  
      IF(KTREN.EQ.2) CALL BTREND(XNOW,TNOW,NGAS)                        9109.6  
      IF(KTREN.EQ.3) CALL CTREND(XNOW,TNOW,NGAS)                        9109.7  
      IF(KTREN.EQ.5) CALL GTREND(XNOW,TNOW,NGAS)
      IF(KTREN.EQ.4) CALL BMTRND(XNOW,TNOW,NGAS)
      IF(KTREN.EQ.6) CALL BMTRNDMG(XNOW,TNOW,NGAS)
c     IF(KTREN.GE.21.and.KTREN.LE.29) 
c    &    CALL STBTRND(XNOW,TNOW,NGAS,KTREN)
c     if(KTREN.LE.5)then
      CALL DTDX1D(XNOW,XREF,XDT0,SDT0,KFOR,NGAS)                        9109.8  
      CALL DTDX3D(XNOW,XREF,XDT0,SDT0,KFOR,1)                           9109.9  
      CALL DXDT3D(XNOW,XREF,XDT0,SDT0,KFOR,NGAS)                        9110.   
c     endif
      FULGAS(2)=XNOW(1)/XREF(1)                                         9110.1  
      FULGAS(6)=XNOW(2)/XREF(2)                                         9110.2  
      FULGAS(7)=XNOW(3)/XREF(3)                                         9110.3  
      FULGAS(8)=XNOW(4)/XREF(4)                                         9110.4  
      FULGAS(9)=XNOW(5)/XREF(5)                                         9110.5  
!     IF(KTREND.EQ.6) then
!        print *,'From forget'
!        print *,'XNOW'
!        print *,XNOW
!        print *,'XREF'
!        print *,XREF
!        print *,'FULGAS(8)=',FULGAS(8),FULGAS(8)*PPMV58(8)
!     endif
c        print *,'From forset XNOW(1)=',XNOW(1)
c        print *,'From forset FULGAS(2)=',FULGAS(2),FULGAS(2)*PPMV58(2)
      if(.not.OBSFOR)then
       IF(KWRITE.EQ.1) then
         print *,'From forset XNOW(1)=',XNOW(1)
       ENDIF
       do j=1,nlat
        atm_co2(j)=XNOW(1)
       enddo
      endif
C                                                                       9110.6  
      IF(KWRITE.NE.1) GO TO 220                                         9110.7  
      SDT0=0.0                                                          9110.8  
      DO 210 I=1,NGAS                                                   9110.9  
      SDT0=SDT0+XDT0(I)                                                 9111.   
      XRAT(I)=(XNOW(I)-XDAT(I))/(1.E-10+XDAT(I))                        9111.1  
  210 XDAT(I)=XNOW(I)                                                   9111.2  
      IYEAR=TNOW                                                        9111.3  
      WRITE(6,6200) IYEAR,SDT0,(XDT0(I),I=1,5),(XNOW(I),I=1,5)          9111.4  
     +             ,(XRAT(I),I=1,5)                                     9111.5  
 6200 FORMAT(1X,I4,F6.3,5F8.4,1X,F8.2,4F8.4,1X,5F8.4)                   9111.6  
      NSPACE=IYEAR-(IYEAR/10)*10                                        9111.7  
      IF(NSPACE.EQ.0) WRITE(6,6010)                                     9111.8  
 6010 FORMAT(1H )                                                       9111.9  
 220  CONTINUE                                                          9112.   
C                                                                       9112.1  
      RETURN                                                            9112.2  
      END                                                               9112.3  
