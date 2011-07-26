#include "ctrparam.h"
      SUBROUTINE RADIA                                                  5001.   
C****                                                                   5002.   
C**** THIS SUBROUTINES ADDS THE RADIATION HEATING TO THE TEMPERATURES   5003.   
C****                                                                   5004.   

#include "BD2G04.COM"
#include "chem_para"
#include "chem_com"

      COMMON U,V,T,P,Q                                                  5006.   
      COMMON/ADDALB/BVSURFA,XVSURFA,BNSURFA,XNSURFA
      COMMON/WORK1/CONV(IM0,JM0,LM0),PK(IM0,JM0,LM0),PREC(IM0,JM0),
     &  TPREC(IM0,JM0),                                                 5007.   
     *  COSZ1(IO0,JM0),COSZ2(IO0,JM0),COSZA(IO0,JM0),                   5008.   
     *  TRINCG(IO0,JM0),BTMPW(IO0,JM0),SNFS(IO0,JM0,4),TNFS(IO0,JM0,4), 5009.   
     *  TRHRS(IO0,JM0,3),SRHRS(IO0,JM0,3),ALB(IO0,JM0,9)                5010.   
      COMMON/WORK2/CLDSS(IM0,JM0,LM0),CLDMC(IM0,JM0,LM0),               5011.   
     *  TOTCLD(36)                                                      5012.   
      DIMENSION TRNFP0(JM0),TRNFP1(JM0),ALBJ(JM0,9)
      real ODATA2(JM0,2),GDATA2(JM0,14),BDATA2(JM0,2),FDATA2(JM0,2),
     *     RQT2(JM0,3)
      common/SURRAD/TRSURF(JM0,4),SRSURF(JM0,4)
      common/FORAERSOL/FORSULF,FORBC,FORVOL,FORGHG
      logical FORSULF,FORBC,FORVOL,FORGHG
C     COMMON/WORK4/ IS BEING USED BY THE RADIATION ROUTINES             5013.   
C                                                                       5014.   
C                   RADCOM:      CONTROL/INPUT PARAMETERS               5015.   
C                                                                       5016.   
      COMMON/RADCOM/VADATA(11,4,3),DGLAT(46),DGLON(72),TMINSR,FULGAS(18)5017.   
     A             ,FRACSL,RATQSL,FOGTSL,PTLISO,TLGRAD,TKCICE,FGOLDU(18)5018.   
     B             ,FLONO3,FRAYLE,FCLDTR,FCLDSR,FALGAE,FMARCL,FEMTRA(6) 5019.   
     C             ,WETTRA,WETSRA,DMOICE,DMLICE,LICETK,NTRCE,FZASRA(6)  5020.   
     D             ,ID5(5),ITR(4),IMG(2),ILG(2),LAPGAS,KWVCON,NORMS0,NV 5021.   
     E             ,KEEPRH,KEEPAL,ISOSCT,IHGSCT,KFRACC,KGASSR,KAERSR    5022.   
     F             ,MARCLD,LAYTOP,LMR,LMRP,JMLAT,IMLON,KFORCE,LASTVC    5023.   
C                                                                       5024.   
C                                BASIC RADCOM INPUT DATA                5025.   
C                                                                       5026.   
     G             ,PLE(40),HLB(40),TLB(40),TLT(40),TL(40),U0GAS(40,9)  5027.   
     H             ,ULGAS(40,9),TRACER(40,4),RTAU(40),QL(40),RHL(40)    5028.   
     I             ,POCEAN,PEARTH,POICE,PLICE,AGESN,SNOWE,SNOWOI,SNOWLI 5029.   
     J             ,TGO,TGE,TGOI,TGLI,TS,WS,WEARTH,ZOICE,FSPARE(200)    5030.   
     K             ,S0,COSZ,PVT(11),BXA(153),SRBXAL(15,2),FRC(5),LUXGAS 5031.   
     L             ,JYEARR,JDAYR,JLAT,ILON,MEANAL,KALVIS,ISPARE(25),SGPS5032.   
C                                                                       5033.   
C                                BASIC RADCOM OUTPUT DATA               5034.   
C                                                                       5035.   
     M             ,TRDFLB(40),TRUFLB(40),TRNFLB(40),TRFCRL(40),TRSLCR  5036.   
     N             ,SRDFLB(40),SRUFLB(40),SRNFLB(40),SRFHRL(40),SRSLHR  5037.   
     O             ,SRIVIS,SROVIS,PLAVIS,SRINIR,SRONIR,PLANIR,SRXATM(4) 5038.   
     P             ,SRDVIS,SRUVIS,ALBVIS,SRDNIR,SRUNIR,ALBNIR,FSRNFG(4) 5039.   
     Q             ,SRTVIS,SRRVIS,SRAVIS,SRTNIR,SRRNIR,SRANIR,FTRUFG(4) 5040.   
     R             ,TRDFGW,TRUFGW,TRUFTW,BTEMPW,TRDFSL,TRUFSL,DTRUFG(4) 5041.   
     S             ,TRSLTS,TRSLTG,TRSLWV,TRSLBS,TTRUFG,LBOTCL,LTOPCL    5042.   
      DIMENSION COE(39)                                                 5043.   
      LOGICAL POLE,DC25,HPRNT,WRCLD,CLDFEED                             5044.
#if ( defined OCEAN_3D )
#include "AGRID.h"
#endif
       dimension SWNET(jm0,2),SWIN(jm0,2)

#if ( defined CLM )
#include "CLM.h"
#endif
c

      common /ATCO2/atm_co2(jm0)
      common/conprn/HPRNT
      common/COMCLD/READGHG,PCLOUD,WRCLD,NWRCLD,NWRCL,INYEAR,JNDAY
     &,CFAEROSOL,ALFA,CFBC,cfvolaer
      dimension STAERMN(JM0,12,2000),STAERGL(2000),SAYEARS(2000)
     &  ,JDY(12)
      DATA JDY/31,59,90,120,151,181,212,243,278,304,334,365/
      common/ghstfor/GHSFALB,GHSF,ALBCF,FVOLADD,STRARFOR,S0FOR,CO2FOR,
     &  CO2F,ghostfv(LM0+1),ghostf(LM0+1,JM0),ghflux(LM0+1,JM0)
      logical GHSFALB,GHSF,GSOEQ,STRARFOR,S0FOR,CO2FOR
      common/eqgso/GSOEQ
      common/cldfdb/coefcl(3),CLDFEED
      common/aexpc/AEXP,ISTRT1
      common/ SNOWALB/FRSNALB,III
      common/ S0XR/S0RATE,CFS0X
      common/ BACKGRGHG/GHGBGR(5)
      COMMON/CO2TRND/ALFFOR,CO2TR,YEARGT,CO2IN,INYRAD
     & ,YEARGG,YEARS0,YEARSA,YEARVA
     & ,KTREND
      dimension CLDSSF(JM0,LM0),CLDMCF(JM0,LM0)
     &,BSO4LAND(JM0),BSO4OCEAN(JM0),BSO4TOTAL(JM0)
     &,bcodmn(JM0,LM0,12)
      dimension DSWSRF(jm0),DLWSRF(jm0),DSWVIS(jm0),DSWNIR(jm0)
     &   ,ULWSRF(jm0)
      integer PCLOUD
!     common/TSUR/TSURFC(JM0,0:13),TSURFT(JM0),TSURFD(JM0),DTSURF(JM0)
!    *,cfcld(JM0,3)
#include "TSRF.COM"
      CHARACTER*4 JMNTHF,JMLAST
      DATA JMLAST /'LAST'/
      DATA TF/273.16/,TCIR/258.16/,STBO/.567257E-7/,IFIRST/1/,JDLAST/-9/5045.   
      DATA IRFIRST /1/
C **** CLEAR SKY
      dimension SRHRCL(JM0),TRHRCL(JM0),ALBCL(JM0),SNP1CL(JM0),
     *SNP0CL(JM0),TRINCL(JM0),TRP0CL(JM0),TRP1CL(JM0)
       dimension RTAU0(LM0),TOTCLD0(LM0)
      common/clrsk/CLEAR(JM0),NCLR(JM0),AJCLR(JM0,12),BJCLR(JM0,12),
     * CJCLR(JM0,12)
      integer CLEAR
      real FULCO2,FULCH4,FULN2O,FULF11,FULF12
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
C **** CLEAR SKY
C****                                                                   5046.   
C**** FDATA  2  LAND COVERAGE (1)                                       5047.   
C****        3  RATIO OF LAND ICE COVERAGE TO LAND COVERAGE (1)         5048.   
C****                                                                   5049.   
C**** ODATA  1  OCEAN TEMPERATURE (C)                                   5050.   
C****        2  RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)       5051.   
C****                                                                   5052.   
C**** GDATA  1  OCEAN ICE SNOW AMOUNT (KG/M**2)                         5053.   
C****        2  EARTH SNOW AMOUNT (KG/M**2)                             5054.   
C****        3  OCEAN ICE TEMPERATURE OF FIRST LAYER (C)                5055.   
C****        4  EARTH TEMPERATURE OF FIRST LAYER (C)                    5056.   
C****        5  EARTH WATER OF FIRST LAYER (KG/M**2)                    5057.   
C****        6  EARTH ICE OF FIRST LAYER (KG/M**2)                      5058.   
C****       11  AGE OF SNOW (DAYS)                                      5059.   
C****       12  LAND ICE SNOW AMOUNT (KG/M**2)                          5060.   
C****       13  LAND ICE TEMPERATURE OF FIRST LAYER (C)                 5061.   
C****                                                                   5062.   
C**** BLDATA 1  COMPOSITE SURFACE WIND MAGNITUDE (M/S)                  5063.   
C****        2  COMPOSITE SURFACE AIR TEMPERATURE (K)                   5064.   
C****        5  FREE                                                    5065.   
C****                                                                   5066.   
C**** VDATA  1-8  EARTH RATIOS FOR THE 8 VEGETATION TYPES (1)           5067.   
C****        9  WATER FIELD CAPACITY OF FIRST LAYER (KG/M**2)           5068.   
C****                                                                   5069.   
       dimension GHGBGR0(5)

      IF(MODRD.EQ.0) IDACC(2)=IDACC(2)+1                                5070.   
      IF (IFIRST.NE.1) GO TO 50                                         5071.   
      if(GHSF)then
       do j=1,jm
        ghflux(1,j)=ghostf(1,j)
        do l=2,lm0+1
         ghflux(l,j)= ghflux(l-1,j)+ghostf(l,j)
        enddo
       enddo
       do l=1,lm0+1
         print *,l,ghflux(l,1)
       enddo
      endif
      JDAYR=JNDAY
      JYEARR=INYEAR
c        GHGBGR(1)=GHGBGR0(1)
c        GHGBGR(2)=GHGBGR0(2)
c        GHGBGR(3)=GHGBGR0(3)
c        GHGBGR(4)=GHGBGR0(4)*1000.
c        GHGBGR(5)=GHGBGR0(5)*1000.
c     print *,'Background GHGs'
      print '(5E12.4)',GHGBGR
         nreadcld=0
         nrbyyr=24*365/5
         nrcldmax=20*nrbyyr
c     print *,' CLOUDS for ',nrcldmax/nrbyyr,' years'
      print *,'RADIA'
!     KTREND=-CO2
      JDAY00=-1
      if(CLDFEED)then
c     print *,' radia.f coefcl=',coefcl
      print *,' for low and middle clouds',coefcl(1)
      print *,' for top  clouds',coefcl(2)
      print *,' for MC  clouds',coefcl(3)
      endif
      print *,' READGHG=',READGHG
      print *,' CFAEROSOL=',CFAEROSOL
#if ( defined PREDICTED_BC && !defined FIXED_FOR )
      if(READGHG.eq.0)then
!     For runs with black carbon predicted by chemistry model
       print *,'With black carbon forcing'
       print *,' CFBC=',CFBC
       nrdbc=0
       read(769),bcodmn
      endif
#endif
      do j=1,jm
        BSO4TOTAL(j)=1.0
        BSO4LAND(j)=1.0
        BSO4OCEAN(j)=1.0
      enddo
      print *,' separate caclulations for land and ocean'
      DC25=.TRUE.
!     DC25=.FALSE.
      if(DC25)then
         print *,' with DC'
      else
         print *,' without DC'
         print *,'   subroutine COSZR'
      end if
      if(abs(PCLOUD-3.).gt.1.5.and..NOT.WRCLD)IFIRST=0                  5072.   
      LMP1=LM+1                                                         5072.1  
      DTCNDS=NCNDS*DT                                                   5073.   
C**** SET THE CONTROL PARAMETERS FOR THE RADIATION                      5074.   
      JMLAT=JM                                                          5074.1  
      if(JM.ne.24) then
      DO J=1,JMLAT
        DGLAT(J)=acos(COSP(J))*360./TWOPI
        if(J.le.JMLAT/2)DGLAT(J)=-DGLAT(J)
      END DO
      endif
      IMLON=IO                                                          5074.2  
      LMR=LM+3                                                          5075.   
      COEX=.01*GRAV*KAPA/RGAS                                           5076.   
      PSFMPT=PSF-PTOP                                                   5077.   
      DO 30 L=1,LM                                                      5078.   
      COE(L)=DTCNDS*COEX/DSIG(L)                                        5079.   
   30 PLE(L)=SIGE(L)*(PSF-PTOP)+PTOP                                    5080.   
      PLE(LMP1)=PTOP                                                    5081.   
      PLE(LM+2)=.5*PTOP                                                 5082.   
      PLE(LMR)=.2*PTOP                                                  5083.   
      PLE(LMR+1)=1.E-5                                                  5084.   
      DO 40 LR=LMP1,LMR                                                 5085.   
      COE(LR)=DT*NRAD*COEX/(PLE(LR)-PLE(LR+1))                          5086.   
      QL(LR)=.3E-5                                                      5087.   
   40 RTAU(LR)=0.                                                       5088.   
      DPMICE=10.                                                        5089.   
C     S0X=1.                                                            5089.1  
      KSTREND=S0RATE*100.
      print *,'S0RATE=',S0RATE,' KSTREND=',KSTREND
      print *,'From first call of radia CFS0X=',CFS0X
      S0X0=S0X
      if(S0FOR)then
        print *,' Initial value for solar constant=',S0X0*1367
      endif

#if ( defined FIXED_FOR )
         print *,'Run with fixed forcing'
#endif
        FGOLDU(1)=1.
C
      if(GSOEQ)then
         if(KSTREND.ne.0)then
           print *,'Wrong setting of GSOEQ and KSTREND'
           print *,' GSOEQ=',GSOEQ,' KSTREND=',KSTREND
           stop
         endif
#if ( defined VOL_AER )
           print *,'Wrong setting of GSOEQ and VOL_AER'
           stop
#endif
         print *,
     &   'An equilibrium run for ',YEARGT,' conditions'
         YEARVA=YEARGT
         YEARS0=YEARGT
        call read_staer (NYVADAT,STAERMN,STAERGL,SAYEARS)
        do n=1,NYVADAT
          if(abs(SAYEARS(n)-YEARVA).lt.0.1) then
            FVOL=STAERGL(n)
           go to 768
          endif
        enddo
  768   continue
        RVOL=0.012
        FGOLDU(1)=(RVOL+FVOL)/RVOL
        print *,'STAEROSOL for year ',YEARVA,' STAER=',FVOL,
     &   '   FGOLDU(1)=',FGOLDU(1)
        FGOLDU1=FGOLDU(1)
         if(KSTREND.eq.0)then
         call obssolar(S0X,YEARS0+0.5)
         S0X=S0X/1367.
          print *,'Sconstant for year ',YEARS0,S0X*1367.
         endif
      else
         if(KSTREND.eq.0)then
           print *,' Run with fixed solar constan ',S0X*1367.
         endif
      endif
#if ( defined VOL_AER )
      print *,'With volcanic forcing cfvolaer=',cfvolaer
        FGOLDU0=1.
        print *,'FGOLDU0=',FGOLDU0
       call read_staer (NYVADAT,STAERMN,STAERGL,SAYEARS)
#else
      if(.not.GSOEQ)then
      print *,'No volcanic forcing FGOLDU(1)=',FGOLDU(1)
      endif
#endif
      if(CO2FOR)then
	print *,' Forcing due to ',CO2F,'xCO2 is calculated',
     *   'in stead of cloud forcing'
      endif
      if(S0FOR)then
	print *,' Forcing due to change in solar constant ',
     *   'is calculated in stead of cloud forcing'
      endif
      if(FORBC)then
	print *,' Forcing due to black carbon ',
     *   'is calculated in stead of cloud forcing'
      endif
      if(FORVOL) then
        print *,'VOL AEROSOL FORCING ',
     *   'is calculated in stead of cloud forcing'
       endif
      if(FORGHG) then
        print *,'GHG FORCING ',
     *   'is calculated in stead of cloud forcing'
       endif
      if(STRARFOR)then
        RVOL=0.012
        FVOL=FVOLADD
        FGOLDU1=(RVOL+FVOL)/RVOL
        print *,'Additional STAER=',FVOL,
     & '   FGOLDU1=',FGOLDU1
      endif
      if(GHSFALB)then
        Print *,' Run with ghost albedo forcing'
        print *,'ALBCF=',ALBCF
      endif
      if(GHSF)then
         Print *,' Run with ghost forcing'
         do l=1,lm+1
            ghostfv(l)=0.0
            do j=1,jm
             ghostfv(l)=ghostfv(l)+ghostf(l,j)
            enddo
            ghostfv(l)=ghostfv(l)/float(jm)
        enddo
        print *,' Ghost forcing DFS=',ghostfv(1)
        do l=LM+1,2,-1
         print *,' Forsing at laeyr',l-1,'=',ghostfv(l)
        enddo
c       print *,' 2D ghost forcing'
      endif
      S0X=CFS0X*S0X
!        print *,'Before CALL RADIA0'
      CALL RADIA0 (IO,JM,CO2,READGHG)                                   5090.   
!        print *,'After CALL RADIA0'
         INCHM=NRAD/NDYN                                                5091.   
C**** CLOUD LAYER INDICES USED FOR DIAGNOSTICS                          5092.   
         DO 43 L=1,LM                                                   5093.   
         LLOW=L                                                         5094.   
         IF (.5*(PLE(L+1)+PLE(L+2)).LT.786.) GO TO 44                   5095.   
   43    CONTINUE                                                       5096.   
   44    LMID1=LLOW+1                                                   5097.   
         DO 45 L=LMID1,LM                                               5098.   
         LMID=L                                                         5099.   
         IF (.5*(PLE(L+1)+PLE(L+2)).LT.430.) GO TO 46                   5100.   
   45    CONTINUE                                                       5101.   
   46    LHI1=LMID+1                                                    5102.   
         LHI=LM                                                         5103.   
         IF (LHI1.GT.LHI) LHI=LHI1                                      5104.   
         WRITE (6,47) LLOW,LMID1,LMID,LHI1,LHI                          5105.   
   47    FORMAT (' LOW CLOUDS IN LAYERS 1-',I2,'   MID LEVEL CLOUDS IN',5106.   
     *     ' LAYERS',I3,'-',I2,'   HIGH CLOUDS IN LAYERS',I3,'-',I2)    5107.   
C**** NO RADIATION AVERAGING          IJRA=1   JRA=1   IRA=1            5108.   
C**** RADIATION AVERAGING IN I             2       1       2            5109.   
C**** RADIATION AVERAGING IN I AND J       4       2       2            5110.   
      JRA=(IJRA+2)/3                                                    5111.   
      IRA=IJRA/JRA                                                      5112.   
   50 JALTER=MOD(NSTEP,NRAD*JRA)/NRAD                                   5113.   
      JDAYR=JDAY
      JYEARR=JYEAR
      IALTER=MOD(NSTEP,NRAD*IJRA)/(NRAD*JRA)                            5114.   
      TNOW=JYEAR+(JDAY-.5)/365.                                         5127.1  
!     print *,'RADIA JDAY=',JDAY
      IF (JDAY.NE.JDLAST.AND.KSTREND.NE.0) then
        IF (KSTREND.GT.0) then
           if(abs(CFS0X-1.0).gt.1.e3)then
             print *,'Wrong setting KSTREND=',KSTREND,' CFS0X=',
     &      CFS0X
           endif
           S0X=S0X0*(1.+S0RATE/70.*(TNOW-INYEAR))
c          print *,'From  radia S0X0=',S0X0,' S0X=',S0X
c          print *,'S0X=',S0X,' TNOW=',TNOW
        ELSE
           call obssolar(S0X,TNOW)
           if(JDAY.eq.1.and.MODRD.EQ.0) then
           print *,'From radia TNOW,S0X',TNOW,S0X
           print *,JDAY,JDLAST
           endif
           S0X=S0X/1367.
        ENDIF
      ENDIF
      S0=S0X*1367./RSDIST                                               5115.   
      S0R=S0
!        print *,'S0=',S0
#if ( defined VOL_AER )
          do MNAER=1,12
           if (JDAY.le.JDY(MNAER)) go to 458
          enddo
 458      continue
c         print *,' MNAER=', MNAER,' MONTH=',JMONTH
#endif
!
C**** CALCULATE AVERAGE COSINE OF ZENITH ANGLE FOR CURRENT COMP3 STEP   5116.   
C****   AND RADIATION PERIOD                                            5117.   
      ROT1=TWOPI*TOFDAY/24.                                             5118.   
      if(DC25)then
         ROT2=ROT1+TWOPI*DTCNDS/SDAY                                    5119.
         CALL COSZT (IO,JM,SIND,COSD,ROT1,ROT2,COSZ1)                   5120.
      else
         ROT2=ROT1+TWOPI
         CALL COSZR (IO,JM,SIND,COSD,ROT1,ROT2,COSZ1)
      end if
      if(HPRNT)then
       print *,' radia TAU=',TAU
       print *,' CLDSS'
       print *,(CLDSS(1,7,L),L=1,LM)
       print *,' CLDMC'
       print *,(CLDMC(1,7,L),L=1,LM)
      endif
!print *,' form radia TAU=',TAU,'MODRD=',MODRD
C
      IF(MODRD.NE.0) GO TO 840                                          5121.   
C
      ROT2=ROT1+TWOPI*NRAD*DT/SDAY                                      5122.   
      CALL COSZS (IO,JM,SIND,COSD,ROT1,ROT2,COSZ2,COSZA)                5123.   
C****                                                                   5124.   
C**** COMPUTE EARTH ALBEDOS AND OTHER PARAMETERS FOR BEGINNING OF DAY   5125.   

      TNOW=JYEAR+(JDAY-.5)/365.                                         5127.1  
      KWRITE=0
      if(JMONTH.ne.JMLAST) then
          KWRITE=1
          print *,' FROM radia READGHG=',READGHG
          if(READGHG.eq.2) call tgases(CO2,JMONTH)
c         print *,(zco2(1,j,1)
c    &           *28.97296245/44.0*1.e-3,j=1,jm)
                        !ppb(m) to kg per volume base
          if(READGHG.eq.1) call rtgases(CO2,JMONTH)
#if ( defined PREDICTED_BC && !defined FIXED_FOR )
        if(READGHG.eq.0)then
          nrdbc=nrdbc+1
          do l=1,lm
           do j=1,jm
            bcod(1,j,l)=bcodmn(j,l,nrdbc)
           enddo
          enddo
          if(nrdbc.eq.12)then
            nrdbc=0.
          endif
        endif
#endif
      endif	! JMONTH.ne.JMLAST
      JMLAST=JMONTH
!print *,' form radia TAU=',TAU,'MODRD=',MODRD
C
      IF (JDAY.NE.JDLAST.AND.KTREND.GT.0) THEN
        if(KTREND.NE.7) then
         CALL FORGET(TNOW,KTREND,KWRITE)
        else
         CALL FORGET(YEARGG,KTREND,KWRITE)
        endif
         if(FORGHG.or.CO2FOR) then
         FULCO2=FULGAS(2)
         FULN2O=FULGAS(6)
         FULCH4=FULGAS(7)
         FULF11=FULGAS(8)
         FULF12=FULGAS(9)
         endif
      ENDIF
      IF (JDAY.NE.JDLAST)then
#ifdef PREDICTED_GASES
        call chemglobal(P)
#endif
#if ( defined PREDICTED_BC && !defined FIXED_FOR )
      call sulfr(BSO4LAND,BSO4OCEAN,1986.)
      do j=1,jm
        FLAND=FDATA(1,J,2)
        BSO4TOTAL(j)=BSO4LAND(j)*FLAND+BSO4OCEAN(j)*(1.-FLAND)
      enddo
#endif
       CALL RCOMPT
       if(CLDFEED)then
	DTSURFAV=0.
C
	do j=1,jm
	 DTSURFAV=DTSURFAV+DT2MGL(J)*DXYP(j)
	end do	!j
	 DTSURFAV=DTSURFAV/AREAG
	do j=1,jm
	do k=1,3
	  cfcld(j,k)=1.+coefcl(k)*DTSURFAV
	end do	! k
	end do 	! j
       endif
      ENDIF	! JDAY.NE.JDLAST
      JDLAST=JDAY                                                       5129.   
         IHOUR=1.5+TOFDAY                                               5130.   
      if(CLDFEED)then
       if (KWRITE.eq.1)then
        print *,'cfcld'
        print 9456,cfcld
        print *,' DTSURF'
        print 9456,DT2MGL
        print *,' DTSURFAV=',DTSURFAV
 9456   format(12f6.2)
       endif
C     Save clouds predicted by model
       do  k=1,LM
        do  j=1,JM
          CLDSSF(j,k)=CLDSS(1,j,k)
          CLDMCF(j,k)=CLDMC(1,j,k)
         enddo
        enddo
       do  k=1,LM
        if(k.le.5)then
          k1=1
        else
          k1=2
        endif
        do  j=1,JM
          CLDSS(1,j,k)=cfcld(j,k1)*CLDSS(1,j,k)
          CLDMC(1,j,k)=cfcld(j,3)*CLDMC(1,j,k)
         enddo
        enddo
       endif
C
CB       READING OF CLOUD
      if(abs(PCLOUD-3.).lt.1.5)then
  910   continue
       if(nreadcld.eq.nrcldmax)go to 900
       read(585,END=900)TFDAYF,JDATEF,JMNTHF,CLDSSF,CLDMCF,IRAND
       nreadcld=nreadcld+1
       if(IFIRST.eq.1)then
         print *,' radia.f PCLOUD=',PCLOUD
         if(PCLOUD.eq.2)print *,' FIXED MC and SS CLOUDS'
         if(PCLOUD.eq.4)print *,' FIXED MC CLOUDS ONLY'
         if(PCLOUD.eq.3)print *,' FIXED SS CLOUDS ONLY'
         print *,TOFDAY,JDATE,JMONTH
         print *,TFDAYF,JDATEF,JMNTHF
         print *,' DTCNDS=',DTCNDS/3600.
         print *,' DT*NRAD=',DT*NRAD/3600.
         if(.not.WRCLD)IFIRST=0
       endif
       if(abs(TOFDAY-TFDAYF).gt.1.e-3.or.JDATE.ne.JDATEF.or.
     *       JMONTH.ne.JMNTHF)then
         print *,' RADIA, disagrement in clouds'
         print *,TOFDAY,JDATE,JMONTH
         print *,TFDAYF,JDATEF,JMNTHF
         stop
       endif
       go to 920
  900  rewind 585
       nreadcld=0
         print *,' END OF file85'
         print *,JYEAR
         print *,TOFDAY,JDATE,JMONTH
         print *,' REWIND 85'
       go to 910
  920  continue
c      if (KWRITE.eq.1)then
c       print *,'cfcld'
c       print 9456,cfcld
c       print *,' DTSURF'
c       print 9456,DTSURF
c9456   format(12f6.2)
c      endif
       CALL RINIT (IRAND)
       do 930 k=1,LM
        do 930 j=1,JM
          if(PCLOUD.ne.4)CLDSS(1,j,k)=CLDSSF(j,k)
          if(PCLOUD.ne.3)CLDMC(1,j,k)=CLDMCF(j,k)
  930  continue
      endif
CE     END OF READING OF CLOUD
      if(WRCLD)then
        if(NWRCLD.eq.1)then
           CALL RFINAL(IRAND)
           if(IFIRST.eq.1)print *,' SHORT CLOUDS RECORD'
           write(581)TOFDAY,JDATE,JMONTH,CLDSS,CLDMC,IRAND
        elseif(NWRCLD.eq.2)then
           if(IFIRST.eq.1)print *,' LONG CLOUDS RECORD'
         do 1150 k=1,14
          do 1150 j=1,JM0
            if(k.le.2)then
             ODATA2(j,k)=ODATA(1,j,k)
             BDATA2(j,k)=BLDATA(1,j,k)
             FDATA2(j,k)=FDATA(1,j,k+1)
            endif
            if(k.le.3)RQT2(j,k)=RQT(1,j,k)
            GDATA2(j,k)=GDATA(1,j,k)
 1150   continue
          CALL RFINAL(IRAND)
          write(581)TOFDAY,JDATE,JMONTH,CLDSS,CLDMC,IRAND,
     *              JDAY,JYEAR,T,Q,P,
     *              ODATA2,BDATA2,FDATA2,GDATA2,RQT2
     &             ,CLDSSF,CLDMCF
#if ( defined CLM )
          write(581)asdirclm,asdifclm,aldirclm,aldifclm
#endif
     
        else
          print *,' NWRCLD=',NWRCLD
          stop
        endif
        IFIRST=0
      endif
#if ( defined CPL_CHEM )
!
! --- Chemistry Model Patch 021199
!       retrive cloud coverages
!       for meta model and others
!
        do kchem=1,nlev
        do jchem=1,nlat
        do ichem=1,nlon
          chem_cldss(ichem,jchem,kchem)
     &       = cldss(ichem,jchem,kchem)
          chem_cldmc(ichem,jchem,kchem)
     &       = cldmc(ichem,jchem,kchem)
        end do
        end do
        end do
!
#endif

C****                                                                   5131.   
C**** MAIN J LOOP                                                       5132.   
C****                                                                   5133.   
      DO 600 J=1,JM                                                     5134.   
      IF ((J-1)*(JM-J).NE.0) GO TO 140                                  5135.   
C**** CONDITIONS AT THE POLES                                           5136.   
      POLE=.TRUE.                                                       5137.   
      MODRJ=0                                                           5138.   
      IMAX=1                                                            5139.   
      GO TO 160                                                         5140.   
C**** CONDITIONS AT NON-POLAR POINTS                                    5141.   
  140 POLE=.FALSE.                                                      5142.   
      MODRJ=MOD(J+JALTER,JRA)                                           5143.   
      IMAX=IM                                                           5144.   
  160 XFRADJ=.2+1.2*COSP(J)*COSP(J)                                     5145.   
#if ( defined VOL_AER )
      RVOL=0.012
      iYAER=SAYEARS(1)
      JYEARAER=min(JYEAR-iYAER+1,NYVADAT)
      FVOL=cfvolaer*STAERMN(J,MNAER,JYEARAER)
      FGOLDU(1)=(RVOL+FVOL)/RVOL
      FGOLDU1=FGOLDU(1)
      if (j.eq.-JM/2)then
        print *,'From radia j=',j
        print *,MNAER,JYEAR,JYEARAER
        print *,cfvolaer
        print *,'RVOL=',RVOL,' FVOL=',FVOL
      print *,'FGOLDU1=',FGOLDU1,' FGOLDU0=',FGOLDU0
      endif
#endif

      JLAT=J                                                            5145.1  
      IF(MODRJ.EQ.0) CALL RCOMPJ                                        5146.   
      SWIN(j,1)=0.0
      SWNET(j,1)=0.0
      SWIN(j,2)=0.0
      SWNET(j,2)=0.0
C****                                                                   5147.   
C**** MAIN I LOOP                                                       5148.   
C****                                                                   5149.   
      IM1=IM                                                            5150.   
      DO 500 I=1,IMAX                                                   5151.   
      MODRIJ=MODRJ+MOD(I+IALTER,IRA)                                    5152.   
      IF(POLE) MODRIJ=0                                                 5153.   
         JR=J
C**** DETERMINE FRACTIONS FOR SURFACE TYPES AND COLUMN PRESSURE         5155.   
      PLAND=FDATA(I,J,2)                                                5156.   
      PWATER=1.-PLAND
      POICE=ODATA(I,J,2)*(1.-PLAND)                                     5157.   
      POCEAN=(1.-PLAND)-POICE                                           5158.   
      if(POCEAN.LE.1.E-5)then
         POCEAN=0.
         POICE=PWATER
      endif
      PLICE=FDATA(I,J,3)*PLAND                                          5159.   
      PEARTH=PLAND-PLICE                                                5160.   
      SP=P(I,J)                                                         5161.   
C****                                                                   5162.   
C**** DETERMINE CLOUDS (AND THEIR OPTICAL DEPTHS) SEEN BY RADIATION     5163.   
C****                                                                   5164.   
      X=999999.                                                         5164.1  
c     RANDSS=RANDU(X)                                                   5165.   
c     RANDMC=RANDU(X)                                                   5166.   
      CALL RANDUU(RANDSS,X)
      CALL RANDUU(RANDMC,X)
      if(HPRNT)then
       print *,' radia TAU=',TAU
       print *,' RANDSS=',RANDSS
       print *,' RANDMC=',RANDMC
      endif
C
         CSS=0.                                                         5167.   
         CMC=0.                                                         5168.   
         DEPTH=0.                                                       5169.   
         LTOP=0                                                         5169.1  
      DO 210 L=1,LM                                                     5170.   
      RTAU(L)=0.                                                        5171.   
  210    TOTCLD(L)=0.                                                   5172.   
      DO 240 L=1,LM                                                     5173.   
      IF(CLDSS(I,J,L).LT.RANDSS) GO TO 220                              5174.   
      RTAUSS=.013333*(PTOP-100.+SIG(L)*SP)                              5175.   
      IF(RTAUSS.LT.0.) RTAUSS=0.                                        5176.   
      IF (T(I,J,L)*PK(I,J,L).LT.TCIR) RTAUSS=.3333333                   5177.   
      RTAU(L)=RTAUSS                                                    5178.   
         CSS=1.                                                         5179.   
         AJL(J,L,28)=AJL(J,L,28)+CSS                                    5180.   
         TOTCLD(L)=1.                                                   5181.   
         LTOP=L                                                         5181.1  
  220 IF(CLDMC(I,J,L).LE.RANDMC) GO TO 240                              5182.   
      RTAUMC=DSIG(L)*SP*.08                                             5183.   
      IF(RTAUMC.GT.RTAU(L)) RTAU(L)=RTAUMC                              5184.   
         CMC=1.                                                         5185.   
         AJL(J,L,29)=AJL(J,L,29)+CMC                                    5186.   
         TOTCLD(L)=1.                                                   5187.   
         LTOP=L                                                         5187.1  
         DEPTH=DEPTH+SP*DSIG(L)                                         5188.   
         AJL(J,L,19)=AJL(J,L,19)+TOTCLD(L)                              5189.   
  240    CONTINUE
         AJ(J,57)=AJ(J,57)+CSS*POCEAN                                   5190.   
         BJ(J,57)=BJ(J,57)+CSS*PLAND                                    5191.   
         CJ(J,57)=CJ(J,57)+CSS*POICE                                    5192.   
         DJ(JR,57)=DJ(JR,57)+CSS*DXYP(J)                                5193.   
         AJ(J,58)=AJ(J,58)+CMC*POCEAN                                   5194.   
         BJ(J,58)=BJ(J,58)+CMC*PLAND                                    5195.   
         CJ(J,58)=CJ(J,58)+CMC*POICE                                    5196.   
         DJ(JR,58)=DJ(JR,58)+CMC*DXYP(J)                                5197.   
         AIJ(I,J,17)=AIJ(I,J,17)+CMC                                    5198.   
         AJ(J,80)=AJ(J,80)+DEPTH*POCEAN                                 5199.   
         BJ(J,80)=BJ(J,80)+DEPTH*PLAND                                  5200.   
         CJ(J,80)=CJ(J,80)+DEPTH*POICE                                  5201.   
         DJ(JR,80)=DJ(JR,80)+DEPTH*DXYP(J)                              5202.   
         CLDCV=CMC+CSS-CMC*CSS                                          5203.   
          cldd4tem(j)=cldd4tem(j)+CLDCV
          ncldd4tem(j)=ncldd4tem(j)+1
         AJ(J,59)=AJ(J,59)+CLDCV*POCEAN                                 5204.   
         BJ(J,59)=BJ(J,59)+CLDCV*PLAND                                  5205.   
         CJ(J,59)=CJ(J,59)+CLDCV*POICE                                  5206.   
         DJ(JR,59)=DJ(JR,59)+CLDCV*DXYP(J)                              5207.   
C**** clear sky condinion
         CSS0=CSS
         CMC0=CMC
         do L=1,LM
          RTAU0(L)=RTAU(L)
          TOTCLD0(L)=TOTCLD(L)
         enddo
      if(CMC.le.0.and.CSS.le.0)then      
         CLEAR(J)=1
      else
         CLEAR(J)=0
      endif
       if(STRARFOR.or.CO2FOR.or.S0FOR.or.FORBC.
     &    or.FORVOL.or.FORGHG)then
         CLEAR(J)=0
       endif
C****                                                                   5247.   
C**** SET UP VERTICAL ARRAYS OMITTING THE I AND J INDICES               5248.   
C****                                                                   5249.   
C**** EVEN PRESSURES                                                    5250.   
      DO 340 L=1,LM                                                     5251.   
      PLE(L)=SIGE(L)*SP+PTOP                                            5252.   
C**** TEMPERATURES                                                      5253.   
      TL(L)=T(I,J,L)*PK(I,J,L)                                          5254.   
C**** MOISTURE VARIABLES                                                5255.   
      QL(L)=Q(I,J,L)                                                    5256.   
  340 CONTINUE                                                          5257.   
C****                                                                   5258.   
C**** RADIATION, SOLAR AND THERMAL                                      5259.   
C****                                                                   5260.   

      DO 420 K=1,3                                                      5261.   
  420 TL(LM+K)=RQT(I,J,K)                                               5262.   
      COSZ=COSZA(I,J)                                                   5263.   
      TGO=ODATA(I,J,1)+TF                                               5264.   
      TGOI=GDATA(I,J,3)+TF                                              5265.   
      TGLI=GDATA(I,J,13)+TF                                             5266.   
      TGE=GDATA(I,J,4)+TF                                               5267.   
      TS=BLDATA(I,J,2)                                                  5268.   
      SNOWOI=GDATA(I,J,1)                                               5269.   
      SNOWLI=GDATA(I,J,12)                                              5270.   
      SNOWE=GDATA(I,J,2)                                                5271.   
      AGESN=GDATA(I,J,11)                                               5272.   
c     print *,'From radia J=',J,' TAU=',TAU,' TS=',TS
c     print *,TGO,TGOI,TGLI,TGE
c     print *,SNOWOI,SNOWLI,SNOWE,AGESN
      WEARTH=(GDATA(I,J,5)+GDATA(I,J,6))/(VDATA(I,J,9)+1.E-20)          5273.   
      DO 430 K=1,8                                                      5274.   
  430 PVT(K)=VDATA(I,J,K)                                               5275.   
      WS=BLDATA(I,J,1)                                                  5276.   
      do 439 L=1,LM+1
       SRHR(I,J,L)=0.
       TRHR(I,J,L)=0.
       if(L.le.4)then
        SNFS(I,J,L)=0.
        TNFS(I,J,L)=0.
        if(L.le.3)then
         SRHRS(I,J,L)=0.
         TRHRS(I,J,L)=0.
        endif
       endif
  439 continue
#if ( defined CPL_CHEM )
!
! --- Chemistry Model Patch 020996
!
        do L=1,LM
          solarflux(i,j,L) = 0.0
        enddo
!
#endif

      TRNFP0(J)=0.
      TRNFP1(J)=0.
      TRINCG(I,J)=0.
      BTMPW(I,J)=0.
      SRDAN=0.
      SRNAN=0.
      do 449 K=1,9
       ALB(I,J,K)=0.
       ALBJ(J,K)=0.
  449 continue
      do 499 ii=1,3
      COSZ=COSZA(I,J)
      PLAND=FDATA(I,J,2)        
      PWATER=1.-PLAND
      POICE=ODATA(I,J,2)*(1.-PLAND)             
      POCEAN=(1.-PLAND)-POICE               
      if(POCEAN.LE.1.E-5)then
         POCEAN=0.
         POICE=PWATER
      endif
      PLICE=FDATA(I,J,3)*PLAND            
      PEARTH=PLAND-PLICE
      if(ii.eq.1)then
       BSO4BC=BSO4OCEAN(J)/BSO4TOTAL(J)
       PTYPE=POCEAN
       POICE=0.
       POCEAN=1.
       PLAND=0.
       PEARTH=0.
       PLICE=0.
       TGAL=0.
      else if(ii.eq.3)then
      BSO4BC=BSO4OCEAN(J)/BSO4TOTAL(J)
       PTYPE=POICE
       POICE=1.
       POCEAN=0.
       PLAND=0.
       PEARTH=0.
       PLICE=0.
       TGAL=TGOI
      else
       BSO4BC=BSO4LAND(J)/BSO4TOTAL(J)
       PTYPE=PLAND
       POCEAN=0.            
       POICE=0.
       PWATER=0.   
       PLICE=FDATA(I,J,3)
       PEARTH=1.-PLICE
       TGAL=TGE*PEARTH+TGLI*PLICE
       PLAND=1.
      endif
      if(PTYPE.lt.1.e-10)go to 499
       if(ii.gt.1)then
       if(TGAL.lt.263.)then
        FRSNALB=0.30
       elseif(TGAL.lt.273.)then
        FRSNALB=0.30-0.015*(TGAL-263.)
       else
        FRSNALB=0.15
       endif
      endif	!ii
c     FGOLDU(2)=XFRADJ*(1.-PEARTH)                                      5277.   
c     FGOLDU(3)=XFRADJ*PEARTH                                           5278.   
      FGOLDU(2)=XFRADJ*(1.-PLAND)
      FGOLDU(3)=XFRADJ*PLAND
      if(FORGHG)then
        FULGAS(2)=FULCO2
        FULGAS(6)=FULN2O
        FULGAS(7)=FULCH4
        FULGAS(8)=FULF11
        FULGAS(9)=FULF12
      endif
      if(CO2FOR)then
        if (KTREND.GT.0) THEN
         FULGAS(2)=FULCO2
        else
         FULGAS(2)=CO2
        endif
      endif
      if(FORVOL)then
        FGOLDU(1)=FGOLDU1
      endif
!
#if ( defined IPCC_EMI )
        FULGAS(2)=CO2
!     IF(DMOD(TAU,365.*24.).EQ.0..and.J.eq.-JM/2) then
      IF(DMOD(TAU,365.*24.).EQ.0..and.J.eq.1) then
        print *,'IPCC_EMI CO2 from radia ', CO2
        print *,'FULGAS(2)=',FULGAS(2)
       print *,"tau=",tau,j,ii
       CALL WRITER (1,0)
      endif
#endif

#if ( defined FIXED_FOR )
        FULGAS(2)=1.
#endif
!
      if(S0FOR)then
        S0=S0R
      endif
      CSS=CSS0
      CMC=CMC0
      do L=1,LM
       RTAU(L)=RTAU0(L)
       TOTCLD(L)=TOTCLD0(L)
      enddo

Chemstry Model Patch 092295
c
        PTYPER=max(PTYPE,5.e-3)
        FAERSOL=PLAND/PTYPER*3.0*CFAEROSOL
C   Factor 3 is a "parameterization" for indirec aerosol effect
c	For normal setting CFAEROSOL=1.0 (default)
c	For hight forcing CFAEROSOL=0.5
c	For low forcing CFAEROSOL=2.0
c
c ==========================
#if ( defined PREDICTED_BC && !defined FIXED_FOR )
        FBC=BSO4BC*CFBC
#endif
      ILON=I                                                            5278.1  
      JLAT=J                                                            5278.2  
      if(GHSFALB)then
         BVSURFA=ALBCF
         XVSURFA=ALBCF
         BNSURFA=ALBCF
         XNSURFA=ALBCF
      else
         BVSURFA=0.0
         XVSURFA=0.0
         BNSURFA=0.0
         XNSURFA=0.0
      endif
      if(STRARFOR)then
c       if(j.gt.JM/2)then
        FGOLDU(1)=FGOLDU1
c       else
c       FGOLDU(1)=1.
c       endif
      endif
      if(J.eq.-10)then
        print *,'From radia GHSFALB=',GHSFALB,' ALBCF=',ALBCF
        print *,'From radia BVSURFA=',BVSURFA
      endif
c     if(J.eq.1)then
c      print *,ii,PTYPE,PWATER,POCEAN,POICE,PLAND,PEARTH,PLICE
c      print *,(FGOLDU(if),if=1,5)
c     endif
c         print *,'S0 for real calculation=',S0
      CALL RCOMPX                                                       5279.   
          if (IRFIRST.eq.1.and.READGHG.eq.1)then
            CALL WRITER(12)
            if(ii.ge.2)IRFIRST=0
          endif
      IF(DMOD(TAU,365.*24.).EQ.0..and.J.eq.JM/2) then
!     IF(DMOD(TAU,30.*24.).EQ.0..and.J.eq.1) then
       print *,' tau=',TAU,' J=',J
       CALL WRITER (1,0)
       print *,-TRNFLB(1),TRNFLB(LM+4),TRNFLB(LM+1),TRDFLB(1)
      endif
c     CALL WRITER (13)                                                  5280.   
      SRHR(I,J,1)=SRHR(I,J,1)+SRNFLB(1)*PTYPE
c     GHOST FORCINGS
      if(GHSF)then
 
        do LFR=1,lm+1
        if(LFR.gt.1)then
          TRFCRL(LFR-1)=TRFCRL(LFR-1)-ghostf(LFR,J)
          TRNFLB(LFR)=TRNFLB(LFR)+ghflux(LFR,J) 
        endif
        enddo
      endif
c
      TRHR(I,J,1)=TRHR(I,J,1)+(STBO*(POCEAN*TGO**4+POICE*TGOI**4
     *  +PLICE*TGLI**4+PEARTH*TGE**4)-TRNFLB(1))*PTYPE
      if(GHSF)then
         TRHR(I,J,1)=TRHR(I,J,1)+ghostf(1,J)*PTYPE
         TRDFLB(1)=TRDFLB(1)+ghostf(1,J)
      endif
C *****
      TRSURF(J,ii)=STBO*(POCEAN*TGO**4+POICE*TGOI**4
     *  +PLICE*TGLI**4+PEARTH*TGE**4)-TRNFLB(1)
!     if(TAU.gt.3623.0.and.TAU.lt.3744.0) then
      if(ii.eq.-2) then
       if(J.eq.29) then
!       print *,'From radia TAU=',TAU
!       print *,'J=',j,' ii=',ii
!       print *,POCEAN,POICE,PLICE,PEARTH
!       print *,TGO,TGOI,TGLI,TGE
!       print *,' TRNFLB(1)=',TRNFLB(1)
!       print *,'TRSURF(J,ii)=',TRSURF(J,ii),' TRDFLB(1)=',TRDFLB(1)
!       print *,' TRUFLB(1)=',TRUFLB(1),' SRBO*T**4=',
!    &     STBO*(POCEAN*TGO**4+POICE*TGOI**4
!    *     +PLICE*TGLI**4+PEARTH*TGE**4)
!        print *,'CLD',CLDCV,CMC,CSS
!        if(ii.eq.2)print *,'TGL=',TGE,STBO*TGE**4
!        write(99,*),'From radia ',TAU,j,ii
!        write(99,*), CLDCV,CMC,CSS
         write(129), TAU,TOFDAY,CLDCV,TGE,TRSURF(J,ii),SRDFLB(1)
!        write(99,*), TAU,TRSURF(J,ii),
!    &     STBO*(POCEAN*TGO**4+POICE*TGOI**4
!    *     +PLICE*TGLI**4+PEARTH*TGE**4),
!    &     TRNFLB(1)
       endif
       if(J.eq.30) then
         write(130), TAU,TOFDAY,CLDCV,TGE,TRSURF(J,ii),SRDFLB(1)
       endif
      endif
!     endif
      if(GHSF)then
         TRSURF(J,ii)=TRSURF(J,ii)+ghostf(1,J)
      endif
      SRSURF(J,ii)=SRNFLB(1)
      DO 440 L=1,LM                                                     5284.   
#if ( defined CPL_CHEM )
!
! --- Chemistry Model Patch 120497
!  get solar flux in w/m^2
!
        solarflux(i,j,l) = solarflux(i,j,l)
     &                   + srdflb(l+1)*ptype

c        solarflux(i,j,l) = solarflux(i,j,l)
c     &                   + srnflb(l+1)*ptype
!
#endif
      SRHR(I,J,L+1)=SRHR(I,J,L+1)+SRFHRL(L)*PTYPE
  440 TRHR(I,J,L+1)=TRHR(I,J,L+1)-TRFCRL(L)*PTYPE
      DO 450 LR=1,3                                                     5287.   
      SRHRS(I,J,LR)=SRHRS(I,J,LR)+SRFHRL(LM+LR)*PTYPE
  450 TRHRS(I,J,LR)=TRHRS(I,J,LR)-TRFCRL(LM+LR)*PTYPE
      DO 460 K=1,4                                                      5290.   
      SNFS(I,J,K)=SNFS(I,J,K)+SRNFLB(K+LM)*PTYPE
  460 TNFS(I,J,K)=TNFS(I,J,K)+(TRNFLB(K+LM)-TRNFLB(1))*PTYPE     
         TRNFP0(J)=TRNFP0(J)+TRNFLB(4+LM)*PTYPE
         TRNFP1(J)=TRNFP1(J)+TRNFLB(1+LM)*PTYPE
         TRINCG(I,J)=TRINCG(I,J)+TRDFLB(1)*PTYPE
         BTMPW(I,J)=BTMPW(I,J)+(BTEMPW-TF)*PTYPE
c        ALB(I,J,1)=SRNFLB(1)/(SRDFLB(1)+1.E-20)                        5295.   
         SRDAN=SRDAN+SRDFLB(1)*PTYPE
         SRNAN=SRNAN+SRNFLB(1)*PTYPE
         ALB(I,J,2)=ALB(I,J,2)+PLAVIS*PTYPE 
         ALB(I,J,3)=ALB(I,J,3)+PLANIR*PTYPE 
         ALB(I,J,4)=ALB(I,J,4)+ALBVIS*PTYPE
         ALB(I,J,5)=ALB(I,J,5)+ALBNIR*PTYPE     
         ALB(I,J,6)=ALB(I,J,6)+SRRVIS*PTYPE        
         ALB(I,J,7)=ALB(I,J,7)+SRRNIR*PTYPE      
         ALB(I,J,8)=ALB(I,J,8)+SRAVIS*PTYPE       
         ALB(I,J,9)=ALB(I,J,9)+SRANIR*PTYPE       
         ALB1=SRNFLB(1)/(SRDFLB(1)+1.E-20)
C **********
         ALBJ(J,2)=PLAVIS
         ALBJ(J,3)=PLANIR
         ALBJ(J,4)=ALBVIS
         ALBJ(J,5)=ALBNIR
         ALBJ(J,6)=SRRVIS
         ALBJ(J,7)=SRRNIR
         ALBJ(J,8)=SRAVIS
         ALBJ(J,9)=SRANIR
         ALBJ(J,1)=SRNFLB(1)/(SRDFLB(1)+1.E-20)
C *********
         COSZ=COSZ2(I,J)
         do L=LM,1,-1
           AJL(J,L,43)=AJL(J,L,43)-TRNFLB(L+1)*PTYPE
           AJL(J,L,42)=AJL(J,L,42)+(SRNFLB(L+1)*COSZ)*PTYPE
         enddo
         if(ii.eq.2)then
#if ( defined CLM )
C        for TEM CLM
            DSWSRF(j)=SRDFLB(1)
            DLWSRF(j)=TRSURF(J,2)
            ULWSRF(j)=TRUFLB(1)
            DSWVIS(j)=SRDVIS
            DSWNIR(j)=SRDNIR
C        for TEM CLM
#endif
          PLAND=PTYPE
          BJ(J,1)=BJ(J,1)+(S0*COSZ)*PLAND
          BJ(J,2)=BJ(J,2)+(SRNFLB(4+LM)*COSZ)*PLAND
c         BSNFS1=BSNFS1+(SRNFLB(1+LM)*COSZ)*PLAND
c         BJ(J,5)=BJ(J,5)+(SRNFLB(1)*COSZ/(ALB1+1.E-20))*PLAND
          BJ(J,5)=BJ(J,5)+(SRDFLB(1)*COSZ)*PLAND
          BJ(J,6)=BJ(J,6)+(SRNFLB(1)*COSZ)*PLAND
          BJ(J,55)=BJ(J,55)+(BTEMPW-TF)*PLAND
          BJ(J,67)=BJ(J,67)+TRDFLB(1)*PLAND
          BJ(J,70)=BJ(J,70)-(TRNFLB(4+LM)-TRNFLB(1))*PLAND
          BJ(J,7)=BJ(J,7)-TRNFLB(4+LM)*PLAND
          BJ(J,8)=BJ(J,8)-TRNFLB(1+LM)*PLAND
          BJ(J,3)=BJ(J,3)+(SRNFLB(1+LM)*COSZ)*PLAND
          BJ(J,71)=BJ(J,71)-(TRNFLB(1+LM)-TRNFLB(1))*PLAND
          DO 761 K=2,9
           BJ(J,K+70)=BJ(J,K+70)+(S0*COSZ)*ALBJ(J,K)*PLAND
  761     CONTINUE
          else if(ii.eq.1)then
          POCEAN=PTYPE
          AJ(J,1)=AJ(J,1)+(S0*COSZ)*POCEAN
          AJ(J,2)=AJ(J,2)+(SRNFLB(4+LM)*COSZ)*POCEAN
          AJ(J,5)=AJ(J,5)+(SRDFLB(1)*COSZ)*POCEAN
          AJ(J,6)=AJ(J,6)+(SRNFLB(1)*COSZ)*POCEAN
          AJ(J,55)=AJ(J,55)+(BTEMPW-TF)*POCEAN
          AJ(J,67)=AJ(J,67)+TRDFLB(1)*POCEAN
          AJ(J,70)=AJ(J,70)-(TRNFLB(4+LM)-TRNFLB(1))*POCEAN
          AJ(J,7)=AJ(J,7)-TRNFLB(4+LM)*POCEAN
          AJ(J,8)=AJ(J,8)-TRNFLB(1+LM)*POCEAN
          AJ(J,3)=AJ(J,3)+(SRNFLB(1+LM)*COSZ)*POCEAN
          AJ(J,71)=AJ(J,71)-(TRNFLB(1+LM)-TRNFLB(1))*POCEAN
#if ( defined OCEAN_3D )
          SWIN(j,1)=SRDFLB(1)
          SWNET(j,1)=SRNFLB(1)
#endif
C
          DO  K=2,9
           AJ(J,K+70)=AJ(J,K+70)+(S0*COSZ)*ALBJ(J,K)*POCEAN
          END DO
C  
          else
          POICE=PTYPE
          CJ(J,1)=CJ(J,1)+(S0*COSZ)*POICE 
          CJ(J,2)=CJ(J,2)+(SRNFLB(4+LM)*COSZ)*POICE 
          CJ(J,5)=CJ(J,5)+(SRDFLB(1)*COSZ)*POICE 
          CJ(J,6)=CJ(J,6)+(SRNFLB(1)*COSZ)*POICE 
          CJ(J,55)=CJ(J,55)+(BTEMPW-TF)*POICE 
          CJ(J,67)=CJ(J,67)+TRDFLB(1)*POICE 
          CJ(J,70)=CJ(J,70)-(TRNFLB(4+LM)-TRNFLB(1))*POICE 
          CJ(J,7)=CJ(J,7)-TRNFLB(4+LM)*POICE 
          CJ(J,8)=CJ(J,8)-TRNFLB(1+LM)*POICE 
          CJ(J,3)=CJ(J,3)+(SRNFLB(1+LM)*COSZ)*POICE 
          CJ(J,71)=CJ(J,71)-(TRNFLB(1+LM)-TRNFLB(1))*POICE 
#if ( defined OCEAN_3D )
          SWIN(j,2)=SRDFLB(1)
          SWNET(j,2)=SRNFLB(1)
#endif
C
          DO  K=2,9
           CJ(J,K+70)=CJ(J,K+70)+(S0*COSZ)*ALBJ(J,K)*POICE
          END DO
          endif
      if(CLEAR(J).eq.0)then
       PLAND=FDATA(I,J,2)        
       PWATER=1.-PLAND
       POICE=ODATA(I,J,2)*(1.-PLAND)             
       POCEAN=(1.-PLAND)-POICE               
       if(POCEAN.LE.1.E-5)then
          POCEAN=0.
          POICE=PWATER
       endif
       PLICE=FDATA(I,J,3)*PLAND            
       PEARTH=PLAND-PLICE
       if(ii.eq.1)then
        PTYPE=POCEAN
        POICE=0.
        POCEAN=1.
        PLAND=0.
        PEARTH=0.
        PLICE=0.
       else if(ii.eq.3)then
        PTYPE=POICE
        POICE=1.
        POCEAN=0.
        PLAND=0.
        PEARTH=0.
        PLICE=0.
       else
        PTYPE=PLAND
        POCEAN=0.            
        POICE=0.
        PWATER=0.   
        PLICE=FDATA(I,J,3)
        PEARTH=1.-PLICE
        PLAND=1.
       endif
       COSZ=COSZA(I,J)
       if(STRARFOR)then
         FGOLDU(1)=1.0
       elseif(CO2FOR)then
         FULGAS(2)=CO2F
       elseif(S0FOR)then
          S0=S0X0/CFS0X*1367./RSDIST
          if(JDAY.eq.1.and.MODRD.EQ.0) then
          if(J.eq.24)then
          print *,'For S0FOR ',S0X0,CFS0X,S0X
          print *,'S0 for forcing calculation=',
     &          S0*RSDIST
          endif
          endif
       elseif(FORBC)then
        FBC=0.0
       elseif(FORGHG)then
        CALL FORGET(YEARGT,KTREND,0)
       elseif(FORVOL)then
        FGOLDU(1)=FGOLDU0
       else
        CSS=0.
        CMC=0.
        DEPTH=0.
        LTOP=0.
        do 1210 L=1,LM
          RTAU(L)=0.
          TOTCLD(L)=0.
 1210   continue
       endif
  568  continue
        CALL RCOMPX
       if(FORGHG.or.CO2FOR)then
        IF(DMOD(TAU,365.*24.).EQ.0..and.J.eq.JM/2) then
!       IF(DMOD(TAU,30.*24.).EQ.0..and.J.eq.1) then
         print *,' For CLEAR tau=',TAU,' J=',J
         CALL WRITER (1,0)
         print *,-TRNFLB(1),TRNFLB(LM+4),TRNFLB(LM+1),TRDFLB(1)
        ENDIF
        endif
      endif  ! CLEAR
       SRHRCL(J)=SRNFLB(1)
       TRHRCL(J)=-TRNFLB(1)
       ALBCL(J)=SRNFLB(1)/(SRDFLB(1)+1.e-20)
       SNP1CL(J)=SRNFLB(LM+1)
       SNP0CL(J)=SRNFLB(LM+4)
       TRINCL(J)=TRDFLB(1)
       TRP0CL(J)=TRNFLB(LM+4)
       TRP1CL(J)=TRNFLB(LM+1)
!     if(J.eq.29) then
!     if(TAU.gt.3623.0.and.TAU.lt.3744.0) then
!       print *,' TRHRCL(j)=',TRHRCL(j)
!       print *,' TRINCL(j)=',TRINCL(j)
!     endif
!      endif
       COSZ=COSZ2(I,J)    
       do L=LM,1,-1
        AJL(J,L,46)=AJL(J,L,46)-TRNFLB(L+1)*PTYPE
        AJL(J,L,45)=AJL(J,L,45)+(SRNFLB(L+1)*COSZ)*PTYPE
       enddo
C *******
         NCLR(J)=NCLR(J)+1
         SRHE5=SRHRCL(J)*COSZ/(ALBCL(J)+1.E-20)
         if(ii.eq.1)then
            POCEAN=PTYPE
            AJCLR(J,1)=AJCLR(J,1)+(S0*COSZ)*POCEAN  
            AJCLR(J,2)=AJCLR(J,2)+(SNP0CL(J)*COSZ)*POCEAN   
            AJCLR(J,4)=AJCLR(J,4)+(SRHRCL(J)*COSZ)*POCEAN  
            AJCLR(J,5)=AJCLR(J,5)+SRHE5*POCEAN
            AJCLR(J,6)=AJCLR(J,6)+TRINCL(J)*POCEAN 
            AJCLR(J,8)=AJCLR(J,8)-TRP0CL(J)*POCEAN
            AJCLR(J,9)=AJCLR(J,9)-TRP1CL(J)*POCEAN
            AJCLR(J,3)=AJCLR(J,3)+(SNP1CL(J)*COSZ)*POCEAN
            AJCLR(J,7)=AJCLR(J,7)+TRHRCL(J)*POCEAN
          else if(ii.eq.2)then
            PLAND=PTYPE
            BJCLR(J,1)=BJCLR(J,1)+(S0*COSZ)*PLAND  
            BJCLR(J,2)=BJCLR(J,2)+(SNP0CL(J)*COSZ)*PLAND      
            BJCLR(J,4)=BJCLR(J,4)+(SRHRCL(J)*COSZ)*PLAND    
            BJCLR(J,5)=BJCLR(J,5)+SRHE5*PLAND
            BJCLR(J,6)=BJCLR(J,6)+TRINCL(J)*PLAND    
            BJCLR(J,8)=BJCLR(J,8)-TRP0CL(J)*PLAND
            BJCLR(J,9)=BJCLR(J,9)-TRP1CL(J)*PLAND
            BJCLR(J,3)=BJCLR(J,3)+(SNP1CL(J)*COSZ)*PLAND
            BJCLR(J,7)=BJCLR(J,7)+TRHRCL(J)*PLAND
          else
            POICE=PTYPE
            CJCLR(J,1)=CJCLR(J,1)+(S0*COSZ)*POICE 
            CJCLR(J,2)=CJCLR(J,2)+(SNP0CL(J)*COSZ)*POICE   
            CJCLR(J,4)=CJCLR(J,4)+(SRHRCL(J)*COSZ)*POICE   
            CJCLR(J,5)=CJCLR(J,5)+SRHE5*POICE
            CJCLR(J,6)=CJCLR(J,6)+TRINCL(J)*POICE  
            CJCLR(J,8)=CJCLR(J,8)-TRP0CL(J)*POICE
            CJCLR(J,9)=CJCLR(J,9)-TRP1CL(J)*POICE
            CJCLR(J,3)=CJCLR(J,3)+(SNP1CL(J)*COSZ)*POICE
            CJCLR(J,7)=CJCLR(J,7)+TRHRCL(J)*POICE
          endif
C *********
  499 continue
        ALB(I,J,1)=SRNAN/(SRDAN+1.E-20)
  500 IM1=I                                                             5304.   
C****                                                                   5305.   
C**** END OF MAIN LOOP FOR I INDEX                                      5306.   
C****                                                                   5307.   
  600 CONTINUE                                                          5345.   
C****                                                                   5346.   
C**** END OF MAIN LOOP FOR J INDEX                                      5347.   
C****                                                                   5348.   
C**** ACCUMULATE THE RADIATION DIAGNOSTICS                              5394.   
C****                                                                   5395.   
  700 DO 780 J=1,JM                                                     5396.   
      DXYPJ=DXYP(J)                                                     5397.   
      IMAX=IM                                                           5398.   
      IF(J.EQ.1.OR.J.EQ.JM) IMAX=1                                      5399.   
      DO 720 L=1,LM                                                     5400.   
         ASRHR=0.                                                       5401.   
         ATRHR=0.                                                       5402.   
      DO 710 I=1,IMAX                                                   5403.   
         ASRHR=ASRHR+SRHR(I,J,L+1)*COSZ2(I,J)                           5404.   
  710    ATRHR=ATRHR+TRHR(I,J,L+1)                                      5405.   
         AJL(J,L,9)=AJL(J,L,9)+ASRHR                                    5406.   
  720    AJL(J,L,10)=AJL(J,L,10)+ATRHR                                  5407.   
         ASNFS1=0.                                                      5408.   
         BSNFS1=0.                                                      5409.   
         CSNFS1=0.                                                      5410.   
         ATNFS1=0.                                                      5411.   
         BTNFS1=0.                                                      5412.   
         CTNFS1=0.                                                      5413.   
      DO 770 I=1,IMAX                                                   5414.   
      SP=P(I,J)                                                         5415.   
      COSZ=COSZ2(I,J)                                                   5416.   
      PLAND=FDATA(I,J,2)                                                5417.   
      PWATER=1.-PLAND
      POICE=ODATA(I,J,2)*(1.-PLAND)                                     5418.   
      POCEAN=(1.-PLAND)-POICE                                           5419.   
      if(POCEAN.LE.1.E-5)then
         POCEAN=0.
         POICE=PWATER
      endif
         JR=J
      DO 740 LR=1,3                                                     5421.   
         ASJL(J,LR,3)=ASJL(J,LR,3)+SRHRS(I,J,LR)*COSZ                   5422.   
  740    ASJL(J,LR,4)=ASJL(J,LR,4)+TRHRS(I,J,LR)                        5423.   
         DJ(JR,1)=DJ(JR,1)+(S0*COSZ)*DXYPJ                              5440.   
         DJ(JR,2)=DJ(JR,2)+(SNFS(I,J,4)*COSZ)*DXYPJ                     5444.   
         DJ(JR,3)=DJ(JR,3)+(SNFS(I,J,1)*COSZ)*DXYPJ                     5448.   
         DJ(JR,5)=DJ(JR,5)+(SRHR(I,J,1)*COSZ/(ALB(I,J,1)+1.E-20))*DXYPJ 5452.   
         DJ(JR,6)=DJ(JR,6)+(SRHR(I,J,1)*COSZ)*DXYPJ                     5456.   
         DJ(JR,55)=DJ(JR,55)+BTMPW(I,J)*DXYPJ                           5460.   
         DJ(JR,67)=DJ(JR,67)+TRINCG(I,J)*DXYPJ                          5464.   
         DJ(JR,70)=DJ(JR,70)-TNFS(I,J,4)*DXYPJ                          5468.   
         DJ(JR,71)=DJ(JR,71)-TNFS(I,J,1)*DXYPJ                          5472.   
  770 CONTINUE                                                          5485.   
  780 CONTINUE                                                          5492.   
C****                                                                   5504.   
C**** UPDATE THE TEMPERATURES BY RADIATION                              5505.   
C****                                                                   5506.   
  800 DO 820 J=1,JM                                                     5507.   
      IMAX=IM                                                           5508.   
      IF(J.EQ.1.OR.J.EQ.JM) IMAX=1                                      5509.   
      DO 820 LR=1,3                                                     5510.   
      DO 820 I=1,IMAX                                                   5511.   
 820  RQT(I,J,LR)=RQT(I,J,LR)+(SRHRS(I,J,LR)*COSZ2(I,J)                 5512.   
     *  +TRHRS(I,J,LR))*COE(LR+LM)                                      5513.   
!     JDLAST=JDAY 
  840 DO 860 J=1,JM                                                     5514.   
#if ( defined HR_DATA )
          dswhr(j)=DSWSRF(j)*COSZ1(1,j)
          dlwhr(j)=DLWSRF(j)
#endif
#if ( defined CLM )
          i=1
          dsw4clm(i,j)=DSWSRF(j)*COSZ1(1,j)
          dlw4clm(i,j)=DLWSRF(j)
          swinr4clm(i,j)=DSWNIR(j)*COSZ1(1,j)
          swvis4clm(i,j)=DSWVIS(j)*COSZ1(1,j)
!         if(TAU.gt.3633.0.and.TAU.lt.3744.0.and.J.eq.29) then
          if(J.eq.-29) then
!          print *,TAU,j
!          print *,'DLWSRF(j)=',DLWSRF(j), 
!    &     'dlw4clm(i,j)=',dlw4clm(i,j) 
            write (229,*)TAU,TOFAY,dsw4clm(i,j),dlw4clm(i,j)
          endif
          if(J.eq.-30) then
            write (230,*)TAU,TOFAY,dsw4clm(i,j),dlw4clm(i,j)
          endif
c    For TEM
          swtd4tem(j)=swtd4tem(j)+S0*COSZ1(1,j)
          swsd4tem(j)=swsd4tem(j)+DSWSRF(j)*COSZ1(1,j)
          nradd4tem(j)=nradd4tem(j)+1
#endif
#if ( defined OCEAN_3D )
          solarinc_ocean(J)=solarinc_ocean(J)+SWIN(j,1)*COSZ1(1,j)
          solarnet_ocean(J)=solarnet_ocean(J)+SWNET(j,1)*COSZ1(1,j)
          solarinc_ice(J)=solarinc_ice(J)+SWIN(j,2)*COSZ1(1,j)
          solarnet_ice(J)=solarnet_ice(J)+SWNET(j,2)*COSZ1(1,j)
          navrado(j)=navrado(j)+1
          navrad(j)=navrad(j)+1
#endif
      IMAX=IM                                                           5515.   
      IF(J.EQ.1.OR.J.EQ.JM) IMAX=1                                      5516.   
      if(HPRNT)then
       if(J.eq.1)then
       print *,' radia TAU=',TAU
       print *,' Before 860'
       print *,(T(1,1,L),L=1,LM)
      endif
      endif
      DO 860 L=1,LM                                                     5517.   
      DO 860 I=1,IMAX                                                   5518.   
#if ( defined CPL_CHEM )
!
        coszangle(i,j) = cosz1(i,j)
!
#endif

  860 T(I,J,L)=T(I,J,L)+(SRHR(I,J,L+1)*COSZ1(I,J)+TRHR(I,J,L+1))        5519.   
     *  *COE(L)/(P(I,J)*PK(I,J,L))                                      5520.   
      if(HPRNT)then
       print *,' radia TAU=',TAU
       print *,' after 860'
       print 'B66',T(1,1,1)
       print *,'COSZ1(1,1)',COSZ1(1,1),' COE(1)=',COE(1)
       print *,'SRHR'
       print *,(SRHR(1,1,L+1),L=1,LM)
       print *,'TRHR'
       print *,(TRHR(1,1,L+1),L=1,LM)
       print *,'T'
       print *,(T(1,1,L),L=1,LM)
      endif
c         if(ncallclm.ge.9)stop
      RETURN                                                            5521.   
      END                                                               5522.  
