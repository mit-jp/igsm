c source sokolov users 75004 Aug 15  2006 /home/sokolov/IGSM2/SRC/condse.F
#include "ctrparam.h"

!	==========================================================
!
!	CONDSE.F:  THIS SUBROUTINE ADDS THE CONTRIBUTIONS TO 
!		  TEMPERATURE AND HUMIDITY CAUSED BY CONDENSATION.
!
!	----------------------------------------------------------
!
!	Author of Chemistry Modules:	Chien Wang
!
!	----------------------------------------------------------
! 
!	Revision History:
!	
!	When	Who		What
!	----	----------	-------	
!	073100	Chien Wang	repack based on CliChem3 and add cpp
!     	091901  Chien Wang      make argument of dlog be in r8 format
!     	092001  Chien Wang      add bc and oc
!     	100201  Chien Wang      Eice =0.35
!	062404	Chien Wang	combine bc, oc code with Andrei's
!	==========================================================

      SUBROUTINE CONDSE(mndriver)                                       3001.
C****                                                                   3002.
C**** THIS SUBROUTINE ADDS THE CONTRIBUTIONS TO TEMPERATURE AND         3003.
C**** HUMIDITY CAUSED BY CONDENSATION.                                  3004.
C****                                                                   3005.

#include "BD2G04.COM"                                                   3006.
C
#if ( defined OCEAN_3D || defined ML_2D)
#include "AGRID.h"
cjrs elimated this com file #include "HRD4OCN.COM"
#endif

#if ( defined CLM )
#include "CLM.h"
#endif
c
#if ( defined CPL_CHEM || defined TRACERS)
#include "chem_para"
#include "chem_com"
#endif
#if ( defined TRACERS)
#include "tracer_com.h"
#endif
#if ( defined CPL_CHEM )
!

      dimension xcfc11   (n3d)
      dimension xcfc12   (n3d)
      dimension xxn2o    (n3d)
      dimension xo3      (n3d)
      dimension xco      (n3d)
      dimension xzco2    (n3d)
      dimension xxno     (n3d)
      dimension xxno2    (n3d)
      dimension xxn2o5   (n3d)
      dimension xhno3    (n3d)
      dimension xch4     (n3d)
      dimension xch2o    (n3d)
      dimension xso2     (n3d)
      dimension xh2so4   (n3d)
      dimension xh2o2    (n3d)

      dimension xhfc134a (n3d)
      dimension xpfc     (n3d)
      dimension xsf6     (n3d)

      dimension xbc      (n3d)
      dimension xoc      (n3d)
#endif
#if ( defined CPL_CHEM || defined TRACERS)

      dimension prec_cnv (nlev)
      dimension prec_str (nlev)
!
#endif

      COMMON/SPEC2/KM,KINC,COEK,C3LAND(IO0,JM0),C3OICE(IO0,JM0)         3006.1
     *  ,C3LICE(IO0,JM0),WMGE(IO0,JM0)                                  3006.2
      COMMON U,V,T,P,Q                                                  3007.
      COMMON/WORK1/CONV(IM0,JM0,LM0),PK(IM0,JM0,LM0),PREC(IM0,JM0)
     & ,TPREC(IM0,JM0),                                                 3008.
     *  UC(IM0,JM0,LM0),VC(IM0,JM0,LM0)                                 3009.
      COMMON/WORK2/CLDSS(IM0,JM0,LM0),CLDMC(IM0,JM0,LM0),UCLD(72,9)
     &,VCLD(72,9),                                                      3010.
     *  ID(8),PL(36),PLE(37),PLK(36),TH(36),TL(36),QL(180),             3011.
     *  UL(8,36),UPL(72,36),VPL(72,36),UPUP(72),VPUP(72),               3012.
     *  UUP(8),RA(8),FMXA(36),DSE(36),TCLA(36),TCUP(36),                3013.
     *  X(72),SIGMA1(36),AJ8(36),AJ13(36),AJ50(36)     ,CXCD(36)        3014.
     *  ,DFMX(36),FCD(36),AETA(36),XCD(36),TO(36),QO(36)                3015.
      COMMON/EPARA/VTH(JM0,LM0),WTH(JM0,LM0),VU(JM0,LM0),VV(JM0,LM0),
     &  DQSDT(JM0,LM0)                                                  3015.5
     * ,DWV(JM0),PHIT(JM0,LM0),TPRIM2(JM0,LM0),WU(JM0,LM0),CKS,CKN      3015.51
     * ,WQ(JM0,LM0),VQ(JM0,LM0)                                         3015.52
      common/fixcld/cldssm(JM0,LM0,0:13),cldmcm(JM0,LM0,0:13),
     &  CLDSST(JM0,LM0),
     &  CLDMCT(JM0,LM0)
      DIMENSION XA(1,JM0),XB(1,JM0),CSDATA(JM0,LM0),CMDATA(JM0,LM0)     3015.53
      DIMENSION SHL(180),SHSAT(36),TSAV(36),SIGMA2(36),TX(1,JM0,LM0)    3015.54
     *,QSAV(LM0)
      EQUIVALENCE (SHL(1),QL(1))                                        3015.55
c     DATA CSDATA/                                                      3015.56
c    * 24.2,20.9,49.2,41.5,42.2,52.3,54.4,52.3,44.6,30.7,26.6,24.6,     3015.57
c    * 24.4,26.4,26.0,29.4,35.2,46.3,45.7,36.3,25.4,32.6,38.3,22.3,     3015.58
c    * 10.2,11.1,32.8,11.2,14.6,24.0,12.1,.1,4*0.,                      3015.59
c    * 5*0.,11.9,37.0,34.2,23.3,19.7,21.2,16.5,                         3015.6
c    * 9.5,10.0,25.3,.1,8*0.,                                           3015.61
c    * 9*0.,.2,13.7,14.7,                                               3015.62
c    * 9.2,9.5,31.5,.7,8*0.,                                            3015.63
c    * 7*0.,.2,.6,.1,19.6,12.3,                                         3015.64
c    * 22.8,39.0,53.0,4.9,7.8,29.4,16.8,5*0.,                           3015.65
c    * 5*0.,4.3,36.4,37.6,17.3,17.1,46.7,29.4,                          3015.66
c    * 13.4,18.7,20.7,15.8,25.4,27.9,25.9,3.1,2.1,2.5,1.7,.4,           3015.67
c    * .8,3.1,3.2,3.6,3.6,11.8,29.5,29.7,25.0,18.3,15.7,16.2,           3015.68
c    * 9.9,15.3,20.2,20.3,20.8,19.8,12.7,8.3,12.4,22.8,23.8,19.1,       3015.69
c    * 22.2,26.9,26.5,20.7,15.3,17.2,18.9,19.7,21.4,19.7,12.5,10.2,     3015.7
c    * 3.5,5.4,9.7,10.8,9.9,7.2,3.5,2.1,1.1,1.1,3.2,6.8,                3015.71
c    * 8.3,7.2,5.2,5.7,7.4,8.6,8.2,9.9,11.5,9.2,4.6,4.4,24*0.           3015.72
c    & /
c    & ,242*0.
c    * ,48*0./
c     DATA CMDATA/                                                      3015.73
c    * 12*0.,                                                           3015.74
c    * 12*0.,                                                           3015.75
c    * 3.0,2.6,2.5,4.2,13.7,11.2,6.7,4.3,6.5,8.8,6.4,5.1,               3015.76
c    * 4.5,4.8,7.2,6.3,6.4,6.0,12.1,12.7,13.5,5.6,3.0,3.3,              3015.77
c    * 2.8,1.9,2.5,2.4,5.5,3.9,7.0,4.1,5.3,8.0,6.1,5.0,                 3015.78
c    * 4.5,4.8,7.1,5.5,4.9,4.7,3.6,5.0,15.3,3.9,1.4,2.1,                3015.79
c    * 4.5,3.2,5.1,14.0,20.8,12.6,9.8,3.6,4.0,7.1,5.9,4.7,              3015.8
c    * 4.3,4.8,7.0,4.4,3.6,6.7,13.0,14.6,27.1,10.2,2.0,3.4,             3015.81
c    * 4.7,4.3,8.7,20.1,23.2,24.4,16.4,3.6,4.0,7.1,5.9,4.7,             3015.82
c    * 4.3,4.8,7.0,4.4,3.6,9.3,29.5,19.7,39.4,17.6,3.4,4.8,             3015.83
c    * 0.,0.,0.,6.9,18.5,23.8,15.7,3.9,4.0,7.2,5.9,4.7,                 3015.84
c    * 4.3,4.8,7.0,4.4,3.5,11.2,24.7,19.0,18.5,0.,0.,0.,                3015.85
c    * 7*0.,.1,1.1,6.3,5.5,4.3,                                         3015.86
c    * 4.0,4.9,7.1,2.2,.5,7*0.,48*0.                                    3015.87
c    & /
c    & ,242*0.,48*0./
c    * ,48*0./
      common/COMCLD/READGHG,PCLOUD
      integer PCLOUD
      DIMENSION DSG0(36)                                                3016.
      LOGICAL POLE,SKIPDI,SKIPIF,HPRNT,CONDL                            3017.
     &  ,INIRINST,BARINST,PRNT,HPRNT1
!     common/conprn/HPRNT,JPR,LPR
      common/conprn/HPRNT1,JPR,LPR
      DATA QUP,DSIGUP,CLH/3*0./                                         3018.
      DATA RVAP/461.5/                                                  3019.
      DATA TF/273.16/,TI/233.16/,IFIRST/1/                              3020.
      dimension RHKP(LM0,jm0),RHNEW(JM0)
      QSAT(TM,PR)=.622*EXP(AXCONS+ELHX*BXCONS*(BYTF-1./TM))/PR          3021.
      QSA1(TM,PR)=.622*EXP(AXQSAT-BXQSAT/TM)/PR                         3021.5
      ERFCPI(XX)=.5-XX*(.548-XX*XX*(.139-.0171*XX*XX))                  3022.

!#if ( defined CPL_CHEM )
#if ( defined CPL_CHEM || defined TRACERS)
!
! --- Formula for calculating the Henry's Law Constant
!
      ehenry (AAA,BBB,TM) = AAA*exp(BBB*(1./TM - 0.0033557))

!
! --- Formula for calculating ratio of aqueous to gaseous
!      R = Ha*R*T*L

! 020196

      raq2gas(ehenryx, TM, qqq) = max(0.0,
     &       8.2e-5*ehenryx*TM*qqq)
!
#endif

C****                                                                   3023.
C**** FDATA  2  LAND COVERAGE (1)                                       3024.
C****                                                                   3025.
C**** ODATA  2  RATIO OF LAND ICE COVERAGE TO LAND COVERAGE (1)         3026.
C****                                                                   3027.
C**** GDATA 11  AGE OF SNOW (DAYS)                                      3028.
C****                                                                   3029.
C**** COMPUTE GLOBAL PARAMETERS                                         3030.
      IDACC(1)=IDACC(1)+1                                               3031.
      IF (IFIRST.NE.1) GO TO 50                                         3032.
      SKIPDI=.TRUE.                                                     3031.1
c     SKIPDI=.FALSE.
      SKIPIF=.TRUE.                                                     3031.2
c     SKIPIF=.FALSE.
      INIRINST=.false.
      INIRINST=.true.
      BARINST=.false.
      BARINST=.true.
      JDIFTS=1 
      JDIFTN=JM
      LMDIFT=3
c     LMDIFT=LM
      EDLET=1.
      EDLEQ=1.
      TWOPI=6.283185
      HPRNT=.FALSE.
      print *,' convection before condensation'
      print *,' PCLOUD=',PCLOUD
      print *,' RHNEW is a function of latitude'
      if(SKIPDI)then
      print *,' without vert. diff. for T and Q '
      else
      print *,' vert. diff. for T and Q in',LMDIFT,' layers '      
      print *,' from ',LAT(JDIFTS)*360./TWOPI,' to ',
     * LAT(JDIFTN)*360./TWOPI
      print *,' EDLET=',EDLET,' EDLEQ=',EDLEQ
      endif
      if(SKIPIF)then
      print *,' without vert. diff. for U and V'
      else
      print *,' vert. diff. for U and V in 3 layers '
      endif
      if(INIRINST)then
	print *,' with correction for SYMMETRIC INSTABILITY'
      else
	print *,' without correction for SYMMETRIC INSTABILITY'
      end if
      if(BARINST)then
	print *,' with correction for BAROTROPIC INSTABILITY'
      else
	print *,' without correction for BAROTROPIC INSTABILITY'
      end if
      IFIRST=0                                                          3033.
      NTRACE=0                                                          3033.1
      JDAY00=JDAY-1
      DTCNDS=NCNDS*DT                                                   3034.
      RH0OLD=.80                                                        3034.1
c     RH0OLD=.65
      RH0=0.925	! 2359
!     RH0=0.9 !   630.11
      RH45=0.875	! 2359
!     RH45=0.925	! 470.10 ! 620.11 ! 630.11
#if ( defined CLM35 )
      RH0=0.875	! 360.14
!     RH0=0.90	! 460.15
!     RH0=0.925	! 560.15

      RH45=0.825 ! 360.14
!     RH45=0.85 ! 460.15
!     RH45=0.875 ! 560.15
#endif
!     RH45=0.85	! 480.10
!     RH45=0.80	! 490.10
      print *,' RH0=',RH0,' RH45=',RH45
      RHAV=0.5*(RH0+RH45)
      DRH=0.5*(RH0-RH45)
      do j = 1,jm0

        rhrad = 3.14159*(-90.+4.*(j-1))/180.
        RHNEW(j) =  RHAV+DRH*cos(4.*rhrad)
!       rhrad = 3.14159*(-90.+4.*(j-1)-15.)/180.
!       RHNEW(j) =  RHAV+DRH*sin(6.*rhrad)

        do l=1,3	! Low clouds
#if ( defined CLM35 )
         RHKP(l,j)=0.875*RHNEW(j)  ! 360.14
!        RHKP(l,j)=0.90*RHNEW(j)  ! 660.15 260.15
!        RHKP(l,j)=0.925*RHNEW(j)  ! 960.15
!        RHKP(l,j)=0.912*RHNEW(j)  ! 160.15
!        RHKP(l,j)=0.885*RHNEW(j)  ! 1060.15
!        RHKP(l,j)=0.890*RHNEW(j)  ! 1160.15, 1360.15, 1460.15
!        RHKP(l,j)=0.895*RHNEW(j)  ! 1260.15
#else
         RHKP(l,j)=0.925*RHNEW(j)  ! 2367
#endif
#if ( !defined CLM )
         RHKP(l,j)=0.95*RHNEW(j)  ! 2905.06
         RHKP(l,j)=0.9375*RHNEW(j)  ! 2906.06
         RHKP(l,j)=0.945*RHNEW(j)  ! 2907.06
#endif
        enddo
        do l=4,6	! Middle clouds
!       do l=4,5	! Middle clouds 460.10
         RHKP(l,j)=0.95*RHNEW(j) ! 2367
        enddo
        do l=7,9	! High clouds
!       do l=6,9	! High clouds 460.10
#if ( defined CLM35 )
         RHKP(l,j)=0.995*RHNEW(j) ! 360.14
!        RHKP(l,j)=0.990*RHNEW(j) ! 660.15
!        RHKP(l,j)=0.985*RHNEW(j) ! 760.15
!        RHKP(l,j)=0.975*RHNEW(j) ! 860.15 960.15 260.15 
!        RHKP(l,j)=0.965*RHNEW(j) ! 1360.15
!        RHKP(l,j)=0.955*RHNEW(j) ! 1460.15
#else
         RHKP(l,j)=0.985*RHNEW(j) ! 2358
#endif
#if ( !defined CLM )
         RHKP(l,j)=0.99*RHNEW(j)  ! 2905.06
         RHKP(l,j)=0.995*RHNEW(j)  ! 2906.06 2907.06
#endif
        enddo
      do l=10,LM
       RHKP(l,j)=1.1
      enddo
      enddo
      print *, ' RHNEW=',RHNEW
      print *, ' RHNEW for j=23,34,46'
      print '(3x,3f10.4)',RHNEW(23),RHNEW(34),RHNEW(46)
      print *, ' RHKP/RHNEW '
      do l=lm,1,-1
        print '(i3,2f10.4)',l,SIG(L)*P(1,23)+PTOP,RHKP(l,23)/RHNEW(23)
      enddo
      CSCALE=.6                                                         3034.3
      IQ1=IM/4+1                                                        3035.
      IQ2=IM/2+1                                                        3036.
      SHA=RGAS/KAPA                                                     3037.
      BXCONS=.622/RGAS                                                  3038.
      AXCONS=DLOG(6.1071)                                               3039.
      CLHE=LHE/SHA                                                      3040.
      BYTF=1./TF                                                        3041.
      DTPERD=DTCNDS/SDAY                                                3042.
      AGESNX=1.-DTPERD/50.                                              3043.
C**** PARAMETERS USED FOR CONVECTION                                    3044.
      print *,' RHMAX=',RHMAX
      RVX=0.                                                            3045.
      BX=RHMAX/DTCNDS                                                   3046.
      IMBY2=1                                                           3047.
      NMAX=MIN(IMBY2,17)                                                3048.
      NMIN=MIN(IQ1,7)                                                   3049.
      BYDELN=1./(NMAX+1-NMIN)                                           3050.
      SL1=0.                                                            3051.
      SL4=0.                                                            3052.
      DO 10 N=NMIN,NMAX                                                 3053.
!      ALOGN=DLOG(FLOAT(N))                                              3054.
      ALOGN=LOG(dble(N)) 
      SL1=SL1+ALOGN                                                     3055.
   10 SL4=SL4+ALOGN*ALOGN                                               3056.
      SL4=SL4-SL1*SL1*BYDELN                                            3057.
      SL1=SL1*BYDELN                                                    3058.
      LMCMM1=LMCM-1                                                     3059.
      DSG0(1)=DSIG(1)                                                   3060.
      DO 40 L=1,LMM1                                                    3061.
      WT=1./(L+1)                                                       3062.
   40 DSG0(L+1)=(1.-WT)*DSG0(L)+WT*DSIG(L+1)                            3063.
   50 IF(DOPK.NE.1.) GO TO 58                                           3064.
C**** CALCULATE PK = P**KAPA                                            3065.
      DO 55 J=1,JM                                                      3066.
      DO 55 I=1,IM                                                      3067.
      SP=P(I,J)                                                         3068.
      DO 55 L=1,LM                                                      3069.
      PK(I,J,L)=EXPBYK(SIG(L)*SP+PTOP)                                  3070.
   55 TX(I,J,L)=T(I,J,L)*PK(I,J,L)                                      3070.1
      DOPK=0.                                                           3071.
   58 CONTINUE                                                          3072.
      if(HPRNT)then
        print *,' condse 1'
      print *,' T(J,L)=',T(1,JPR,LPR),' Q(J,L)=',Q(1,JPR,LPR)
      print *,' V(J,L)=',V(1,JPR,LPR),' V(J+1,L)=',V(1,JPR+1,LPR)
      print *,' U(J,L)=',U(1,JPR,LPR),' U(J+1,L)=',U(1,JPR+1,LPR)
      endif
C                                                                       3072.1
C     DO INTERNAL FRICTION FIRST                                        3072.11
C                                                                       3072.12
      IF (SKIPIF) GO TO 66                                              3072.13
      FMU=2.                                                            3072.14
      FCOEF=GRAV*GRAV*FMU*DTCNDS/RGAS                                   3072.15
c     DO 65 J=JDIF,JM-JDIF+2                                            3072.16
      do 65 J=2,JM 
      I=IM                                                              3072.17
      DO 65 IPINC=1,IM                                                  3072.18
      SP=.25*(P(I,J)+P(IPINC,J)+P(I,J-1)+P(IPINC,J-1))                  3072.19
      FCOEF1=FCOEF/(SP*SP)                                              3072.2
      UDN=U(I,J,1)                                                      3072.21
      VDN=V(I,J,1)                                                      3072.22
      TDN=.25*(TX(I,J,1)+TX(IPINC,J,1)+TX(I,J-1,1)+TX(IPINC,J-1,1))     3072.23
c     DO 60 L=2,LM                                                      3072.24
      DO 60 L=2,3
      LM1=L-1                                                           3072.25
      UTP=U(I,J,L)                                                      3072.26
      VUP=V(I,J,L)                                                      3072.27
      TUP=.25*(TX(I,J,L)+TX(IPINC,J,L)+TX(I,J-1,L)+TX(IPINC,J-1,L))     3072.28
      PEUV=SIGE(L)*SP+PTOP                                              3072.29
      RHO=PEUV/(RGAS*.5*(TUP+TDN))                                      3072.3
      TEMP=FCOEF1*(UTP-UDN)*RHO*RHO*RGAS/DSIGO(LM1)                     3072.31
      U(I,J,L)=U(I,J,L)-TEMP/DSIG(L)                                    3072.32
      U(I,J,LM1)=U(I,J,LM1)+TEMP/DSIG(LM1)                              3072.33
      TEMP=FCOEF1*(VUP-VDN)*RHO*RHO*RGAS/DSIGO(LM1)                     3072.34
      V(I,J,L)=V(I,J,L)-TEMP/DSIG(L)                                    3072.35
      V(I,J,LM1)=V(I,J,LM1)+TEMP/DSIG(LM1)                              3072.36
      UDN=UTP                                                           3072.37
      VDN=VUP                                                           3072.38
   60 TDN=TUP                                                           3072.39
   65 I=IPINC                                                           3072.4
   66 CONTINUE                                                          3072.41
      if(HPRNT)then
        print *,' condse 2'
      print *,' J=',JPR,' L=',LPR
      print *,' T(J,L)=',T(1,JPR,LPR),' Q(J,L)=',Q(1,JPR,LPR)
      print *,' V(J,L)=',V(1,JPR,LPR),' V(J+1,L)=',V(1,JPR+1,LPR)
      print *,' U(J,L)=',U(1,JPR,LPR),' U(J+1,L)=',U(1,JPR+1,LPR)
      endif
C                                                                       3072.5
C     PARA. SYMMETRIC INSTABILITY AND BAROTROPIC INSTABILITY            3072.51
C                                                                       3072.52
c     JHALF=JM/2                                                        3072.53
c     JHAM1=JHALF-1                                                     3072.54
c     JHAP3=JHALF+3                                                     3072.55
c     JHAP2=JHALF+2                                                     3072.56
      JVHALF=JM/2+1
      JBAND=4
      if(JM.eq.46)JBAND=8
      JIB=JVHALF-JBAND
      JIE=JVHALF+JBAND-1
      JBB=JIB
      JBE=JIE+1
      DO 168 NITER=1,3                                                  3072.57
      if(HPRNT)then
        print *,' condse 2.1 NITER=',NITER
      print *,' U(J,L)=',U(1,JPR,LPR),' U(J+1,L)=',U(1,JPR+1,LPR)
      endif
	if(INIRINST) then
c     DO 69 J=JHAM1,JHAP2                                               3072.58
      DO 69 J=JIB,JIE
      FTEM=F(J)/DXYP(J)                                                 3072.59
      DO 69 L=1,LM                                                      3072.6
      DUDY=(U(1,J+1,L)*COSV(J+1)-U(1,J,L)*COSV(J))/DYP(J)/COSP(J)       3072.61
      CRI=FTEM*(FTEM-DUDY)                                              3072.62
      IF(CRI.GE.0.) GO TO 69                                            3072.63
      if(HPRNT.and.(J.eq.JPR.or.J.eq.JPR-1))then
       print *,' TAU=',TAU,' J=',J,' L=',L,' NITER=',NITER,' f=',FTEM
       print *,' COSV(J)=',COSV(J),' COSV(J+1)=',COSV(J+1)
       print *,' DYP(J)=',DYP(J),' COSP(J)=',COSP(J)
       print *,' f-dudy=',FTEM-DUDY,' (f-dudy)/f',(FTEM-DUDY)/FTEM
       print *,' U(J,L)=',U(1,J,L),' U(J+1,L)=',U(1,J+1,L)
      endif
      USUM=U(1,J+1,L)*COSV(J+1)+U(1,J,L)*COSV(J)                        3072.64
      U(1,J+1,L)=.5*(FTEM*COSP(J)*DYP(J)+USUM)/COSV(J+1)                3072.65
      U(1,J,L)=(USUM-U(1,J+1,L)*COSV(J+1))/COSV(J)                      3072.66
      if(HPRNT.and.(J.eq.JPR.or.J.eq.JPR-1))then
       print *,' USUM=',USUM
       print *,' U(J,L)=',U(1,J,L),' U(J+1,L)=',U(1,J+1,L)
       print *,' USUMN=',U(1,J+1,L)*COSV(J+1)+U(1,J,L)*COSV(J)
      endif
   69 CONTINUE                                                          3072.67
	end if 	!	INIRINST
	if(BARINST) then
C     BAROTROPIC INSTABILITY                                            3072.68
c     DO 68 J=JHAM1,JHAP3                                               3072.69
      DO 68 J=JBB,JBE
      BETA=(F(J)/DXYP(J)-F(J-1)/DXYP(J-1))/DYV(J)                       3072.7
      DO 68 L=1,LM                                                      3072.73
      PSI=BETA-(U(1,J+1,L)*COSV(J+1)-U(1,J,L)*COSV(J))/                 3072.74
     *  (DYP(J)*DYP(J)*COSP(J))+(U(1,J,L)*COSV(J)-                      3072.75
     *  U(1,J-1,L)*COSV(J-1))/(DYP(J-1)*DYP(J-1)*COSP(J-1))             3072.76
      IF(PSI.GE.0.) GO TO 68                                            3072.77
      if(HPRNT.and.(J.eq.JPR.or.J.eq.JPR-1.or.J.eq.JPR+1))then
      print *,' TAU=',TAU,' J=',J,' L=',L,'NITER=',NITER
      print *,' BETA=',BETA,' PSI=',PSI,' PSI/BETA=',PSI/BETA
      print *,' BETAP1=',BETAP1,' BETAM1=',BETAM1
      print *,' U(J-1,L)=',U(1,J-1,L)
      print *,' U(J,L)=',U(1,J,L),' U(J+1,L)=',U(1,J+1,L)
      endif
      USUM=U(1,J,L)*COSV(J)+U(1,J+1,L)*COSV(J+1)+U(1,J-1,L)*COSV(J-1)   3072.86
      IF ( J.NE.JM/2+1)THEN
      BJJ=1./COSP(J)/DYV(J)
      BJM1=1./COSP(J-1)/DYV(J)
      DJP1=COSV(J)*(BJJ+2.*BJM1)/(COSV(J+1)*(BJJ-BJM1))
      CJP1=(BETA*DYV(J)-USUM*BJM1)/(COSV(J+1)*(BJJ-BJM1))
      DJM1=COSV(J)*(BJM1+2.*BJJ)/(COSV(J-1)*(BJM1-BJJ))
      CJM1=(BETA*DYV(J)-USUM*BJJ)/(COSV(J-1)*(BJM1-BJJ))
      U(1,J,L)=(COSV(J+1)*DJP1*(U(1,J+1,L)-CJP1)+COSV(J)*U(1,J,L)+
     *      COSV(J-1)*DJM1*(U(1,J-1,L)-CJM1))/
     *     (COSV(J+1)*DJP1**2+COSV(J)+COSV(J-1)*DJM1**2)
      U(1,J+1,L)=DJP1*U(1,J,L)+CJP1
      U(1,J-1,L)=DJM1*U(1,J,L)+CJM1
      ELSE
       U(1,J,L)=(USUM-BETA*COSP(J)*DYV(J)**2)/(3.*COSV(J))
       U(1,J+1,L)=1./COSV(J-1)*(USUM-COSV(J)*U(1,J,L)-
     *   COSV(J-1)*(U(1,J-1,L)-U(1,J+1,L)))/
     *   (1.+COSV(J+1)/COSV(J-1))
       U(1,J-1,L)=(USUM-COSV(J)*U(1,J,L)-COSV(J+1)*U(1,J+1,L))/
     *   COSV(J-1)
      ENDIF
      if(HPRNT.and.(J.eq.JPR.or.J.eq.JPR-1.or.J.eq.JPR+1))then
        print *,' TAU=',TAU,' J=',J,' L=',L,'NITER=',NITER
        print *,' USUM=',USUM
        print *,' COSV(J-1)=',COSV(J-1),' FUNM=',FUNM
        print *,' COSP(J)=',COSP(J),' COSP(J-1)=',COSP(J-1)
        print *,' DYV(J)=',DYV(J)
        print *,' U(J-1,L)=',U(1,J-1,L)
        print *,' U(J,L)=',U(1,J,L),' U(J+1,L)=',U(1,J+1,L)
        print *,' USUMN=',U(1,J,L)*COSV(J)+U(1,J+1,L)*COSV(J+1)+
     *    U(1,J-1,L)*COSV(J-1)
      endif
   68 CONTINUE                                                          3072.96
      if(HPRNT)then
        print *,' condse 2.2 NITER=',NITER
      print *,' U(J,L)=',U(1,JPR,LPR),' U(J+1,L)=',U(1,JPR+1,LPR)
      endif
	end if 	!	BARINST
  168 continue
      if(HPRNT)then
        print *,' condse 3'
      print *,' T(J,L)=',T(1,JPR,LPR),' Q(J,L)=',Q(1,JPR,LPR)
      print *,' V(J,L)=',V(1,JPR,LPR),' V(J+1,L)=',V(1,JPR+1,LPR)
      print *,' U(J,L)=',U(1,JPR,LPR),' U(J+1,L)=',U(1,JPR+1,LPR)
      endif
C**** SAVE UC AND VC, AND ZERO OUT CLDSS AND CLDMC                      3073.
   70 DO 75 L=1,LM                                                      3074.
      DO 75 J=1,JM                                                      3075.
      DO 75 I=1,IM                                                      3076.
      UC(I,J,L)=U(I,J,L)                                                3077.
      VC(I,J,L)=V(I,J,L)                                                3078.
      CLDSS(I,J,L)=0.                                                   3079.
   75 CLDMC(I,J,L)=0.                                                   3080.
         IHOUR=1.5+TOFDAY                                               3081.
C****                                                                   3082.
C**** MAIN J LOOP                                                       3083.
C****                                                                   3084.
      DO 810 J=1,JM                                                     3085.
      JHALF=JM/2                                                        3085.5
      COEKD=CKS                                                         3085.6
      IF(J.GT.JHALF) COEKD=CKN                                          3085.7
      IF ((J-1)*(JM-J).NE.0) GO TO 90                                   3086.
C**** CONDITIONS AT THE POLES                                           3087.
      POLE=.TRUE.                                                       3088.
      IMAX=1                                                            3089.
      IF(J.EQ.JM) GO TO 80                                              3090.
      JP=2                                                              3091.
      JVPO=2                                                            3092.
      RAPO=2.*RAPVN(1)                                                  3093.
      RA(1)=RAPO
      GO TO 160                                                         3094.
   80 JP=JMM1                                                           3095.
      JVPO=JM                                                           3096.
      RAPO=2.*RAPVS(JM)                                                 3097.
      RA(1)=RAPO
      GO TO 160                                                         3098.
C**** CONDITIONS AT NON-POLAR POINTS                                    3099.
   90 POLE=.FALSE.                                                      3100.
      JP=J                                                              3101.
      IMAX=IM                                                           3102.
      DO 100 K=1,2                                                      3103.
      RA(K)=RAPVS(J)                                                    3104.
 100  RA(K+2)=RAPVN(J)                                                  3105.
C**** STANDARD DEVIATION FOR TEMPERATURE                                3106.
  160 DO 150 L=1,LM                                                     3107.
      TVAR=0.                                                           3108.
      SUMT=0.                                                           3109.
      PKJ=0.                                                            3110.
      DO 110 I=1,IM                                                     3111.
      PKJ=PKJ+PK(I,JP,L)                                                3112.
  110 SUMT=SUMT+T(I,JP,L)                                               3113.
      IF(KM.EQ.1) GO TO 149                                             3113.5
      DO 120 I=1,IM                                                     3114.
      TDEV=T(I,JP,L)-SUMT/FIM                                           3115.
      X(I)=TDEV                                                         3116.
  120 TVAR=TVAR+TDEV*TDEV                                               3117.
      TVAR=TVAR/FIM                                                     3118.
c     CALL FRTR(X)                                                      3119.
      SL2=0.                                                            3120.
      SL3=0.                                                            3121.
      DO 130 N=NMIN,NMAX                                                3122.
c     ALOGA=DLOG(X(N)+1.E-20)                                           3123.
      ALOGA=LOG(X(N)+1.E-20) 
      SL2=SL2+ALOGA                                                     3124.
      FN=N
  130 SL3=SL3+ALOGA*LOG(FN)
c 130 SL3=SL3+ALOGA*DLOG(FLOAT(N))                                      3125.
      SLOPE=(SL1*SL2-SL3)/SL4                                           3126.
      IF (SLOPE.LT.1.67) SLOPE=1.67                                     3127.
      IF (SLOPE.GT.3.) SLOPE=3.                                         3128.
      SUMXN=0.                                                          3129.
      DO 140 N=1,IMBY2                                                  3130.
  140 SUMXN=SUMXN+X(N)                                                  3131.
      SUMAMK=0.                                                         3132.
      DO 145 N=NMIN,NMAX                                                3133.
  145 SUMAMK=SUMAMK+X(N)*(N**SLOPE)                                     3134.
      SLOPM1=SLOPE-1.                                                   3135.
      XEPE=2.*SUMAMK*BYDELN/((SUMXN+1.E-20)*SLOPM1*(IQ2**SLOPM1))       3136.
 149  SIGMA1(L)=1.4142*SQRT(TPRIM2(JP,L))*PKJ/FIM                       3137.
      SIGMA2(L)=SIGMA1(L)                                               3137.1
!     SIGMA2(L)=2.*PKJ/FIM
  150    CONTINUE                                                       3138.
C****                                                                   3139.
C**** MAIN I LOOP                                                       3140.
C****                                                                   3141.
      IM1=IM                                                            3142.
      DO 700 I=1,IMAX                                                   3143.
         JR=J
C****                                                                   3145.
C**** SET UP VERTICAL ARRAYS, OMITTING THE J AND I SUBSCRIPTS           3146.
C****                                                                   3147.
      PLAND=FDATA(I,J,2)                                                3148.
      PWATER=1.-PLAND
      POICE=ODATA(I,J,2)*(1.-PLAND)                                     3149.
      POCEAN=(1.-PLAND)-POICE                                           3150.
      if(POCEAN.LE.1.E-5)then
         POCEAN=0.
         POICE=PWATER
      endif
!      07/22/2005
      if (pland.lt.1.0)then
       PRLAND=prlnd2total(j,mndriver)
       PROCEAN=(1.-pland*prlnd2total(j,mndriver))
     &      /(1.-pland)
      else
       PRLAND=1.0
       PROCEAN=0.0
      endif
!
C**** PRESSURES, AND PRESSURE TO THE KAPA                               3151.
      SP=P(I,J)                                                         3152.
      DO 170 L=1,LM                                                     3153.
      PL(L)=SIG(L)*SP+PTOP                                              3154.
      PLK(L)=PK(I,J,L)                                                  3155.
C**** TEMPERATURES                                                      3156.
      TH(L)=T(I,J,L)                                                    3157.
      TL(L)=TH(L)*PLK(L)                                                3158.
      QL(L)=Q(I,J,L)                                                    3158.1
      TSAV(L)=TL(L)
      QSAV(L)=QL(L)

#if ( defined CPL_CHEM )
!
      xcfc11(l)=cfc11(i,j,l)
      xcfc12(l)=cfc12(i,j,l)
      xxn2o (l)=xn2o (i,j,l)
      xo3   (l)=o3   (i,j,l)
      xco   (l)=co   (i,j,l)
      xzco2 (l)=zco2 (i,j,l)
      xxno  (l)=xno  (i,j,l)
      xxno2 (l)=xno2 (i,j,l)
      xxn2o5(l)=xn2o5(i,j,l)
      xhno3 (l)=hno3 (i,j,l)
      xch4  (l)=ch4  (i,j,l)
      xch2o (l)=ch2o (i,j,l)
      xso2  (l)=so2  (i,j,l)
      xh2so4(l)=h2so4(i,j,l)
c 062295
      xh2o2 (l)=h2o2 (i,j,l)

        ! === if hfc, pfc, and sf6 are included:
#ifdef INC_3GASES
	! === 032698
      xhfc134a(l) = hfc134a(i,j,l)
      xpfc    (l) = pfc(i,j,l)
      xsf6    (l) = sf6(i,j,l)
	! ===
#endif

      xbc     (l) = bcarbon(i,j,l)
      xoc     (l) = ocarbon(i,j,l)
!
#endif
#if ( defined TRACERS)
         do ntr=1,ntracers
           tracxc(l,ntr)=tracers(i,j,l,ntr)
       enddo
#endif

 170  CONTINUE
      if(HPRNT)then
        print *,' condse after 170 J=',J
      print *,' SP=',SP
      print *,(TH(L),L=1,LM)
      print *,(QL(L),L=1,LM)
      endif
C                                                                       3158.11
C     DO VERTICAL HEAT AND MOISTURE DIFFUSION FIRST                     3158.12
C                                                                       3158.13
      IF (SKIPDI) GO TO 195                                             3158.14
      IF(J.LT.JDIFTS.OR.J.GT.JDIFTN) GO TO 195
      DO 190 LM1=2,LMDIFT                                               3158.15
      L=LM1-1                                                           3158.16
      DTETA=(TH(LM1)-TH(L))*(PLK(LM1)+PLK(L))*.5                        3158.17
      DZUP=SP*DSIG(LM1)*RGAS*TL(LM1)/(PL(LM1)*GRAV)                     3158.18
      DZDN=SP*DSIG(L)*RGAS*TL(L)/(PL(L)*GRAV)                           3158.19
c     EDLE=2.                                                           3158.2
      TEMP=DTCNDS*(DSIG(LM1)+DSIG(L))/(DZUP+DZDN)**2.                   3158.21
      FLE=-2.*EDLET*DTETA*TEMP                                          3158.22
      TL(LM1)=TL(LM1)+FLE/DSIG(LM1)                                     3158.23
      TL(L)=TL(L)-FLE/DSIG(L)                                           3158.24
      TH(LM1)=TL(LM1)/PLK(LM1)                                          3158.25
      TH(L)=TL(L)/PLK(L)                                                3158.26
      DSH=QL(LM1)-QL(L)                                                 3158.27
      ELE=-2.*EDLEQ*DSH*TEMP                                            3158.28
      QL(LM1)=QL(LM1)+ELE/DSIG(LM1)                                     3158.29
      QL(L)=QL(L)-ELE/DSIG(L)                                           3158.3

#if ( defined CPL_CHEM || defined TRACERS)
!
        xxx = -2.0*temp
        xxm1= xxx/dsig(LM1)
        xxL = xxx/dsig(L)
#endif

#if ( defined CPL_CHEM )
!
        ele = (xcfc11(lm1)-xcfc11(l))
        xcfc11(lm1)=xcfc11(lm1)+ele*xxm1
        xcfc11(l)  =xcfc11(l)  -ele*xxL

        ele = (xcfc12(lm1)-xcfc12(l))
        xcfc12(lm1)=xcfc12(lm1)+ele*xxm1
        xcfc12(l)  =xcfc12(l)  -ele*xxL

        ele = (xxn2o (lm1)-xxn2o (l))
        xxn2o (lm1)=xxn2o (lm1)+ele*xxm1
        xxn2o (l)  =xxn2o (l)  -ele*xxL

        ele = (xo3   (lm1)-xo3   (l))
        xo3   (lm1)=xo3   (lm1)+ele*xxm1
        xo3   (l)  =xo3   (l)  -ele*xxL

        ele = (xco   (lm1)-xco   (l))
        xco   (lm1)=xco   (lm1)+ele*xxm1
        xco   (l)  =xco   (l)  -ele*xxL

        ele = (xzco2 (lm1)-xzco2 (l))
        xzco2 (lm1)=xzco2 (lm1)+ele*xxm1
        xzco2 (l)  =xzco2 (l)  -ele*xxL

        ele = (xxno  (lm1)-xxno  (l))
        xxno  (lm1)=xxno  (lm1)+ele*xxm1
        xxno  (l)  =xxno  (l)  -ele*xxL

        ele = (xxno2 (lm1)-xxno2 (l))
        xxno2 (lm1)=xxno2 (lm1)+ele*xxm1
        xxno2 (l)  =xxno2 (l)  -ele*xxL

        ele = (xxn2o5(lm1)-xxn2o5(l))
        xxn2o5(lm1)=xxn2o5(lm1)+ele*xxm1
        xxn2o5(l)  =xxn2o5(l)  -ele*xxL

        ele = (xhno3 (lm1)-xhno3 (l))
        xhno3 (lm1)=xhno3 (lm1)+ele*xxm1
        xhno3 (l)  =xhno3 (l)  -ele*xxL

        ele = (xch4  (lm1)-xch4  (l))
        xch4  (lm1)=xch4  (lm1)+ele*xxm1
        xch4  (l)  =xch4  (l)  -ele*xxL

        ele = (xch2o (lm1)-xch2o (l))
        xch2o (lm1)=xch2o (lm1)+ele*xxm1
        xch2o (l)  =xch2o (l)  -ele*xxL

        ele = (xso2  (lm1)-xso2  (l))
        xso2  (lm1)=xso2  (lm1)+ele*xxm1
        xso2  (l)  =xso2  (l)  -ele*xxL

        ele = (xh2so4(lm1)-xh2so4(l))
        xh2so4(lm1)=xh2so4(lm1)+ele*xxm1
        xh2so4(l)  =xh2so4(l)  -ele*xxL

        ! === if hfc, pfc, and sf6 are included:
#ifdef INC_3GASES
	! === 032698
        ele = (xhfc134a(lm1)-xhfc134a(l))
        xhfc134a(lm1)=xhfc134a(lm1)+ele*xxm1
        xhfc134a(l)  =xhfc134a(l)  -ele*xxL

        ele = (xpfc(lm1)-xpfc(l))
        xpfc(lm1)=xpfc(lm1)+ele*xxm1
        xpfc(l)  =xpfc(l)  -ele*xxL

        ele = (xsf6(lm1)-xsf6(l))
        xsf6(lm1)=xsf6(lm1)+ele*xxm1
        xsf6(l)  =xsf6(l)  -ele*xxL
	! ===
#endif

	ele = (xbc(lm1)-xbc(l))
	xbc(lm1)=xbc(lm1)+ele*xxm1
	xbc(l)  =xbc(l)  -ele*xxL

	ele = (xoc(lm1)-xoc(l))
	xoc(lm1)=xoc(lm1)+ele*xxm1
	xoc(l)  =xoc(l)  -ele*xxL

c 062295
c        ele = (xh2o2(lm1)-xh2o2(l))
c        xh2o2(lm1)=xh2o2(lm1)+ele*xxm1
c        xh2o2(l)  =xh2o2(l)  -ele*xxL

!
#endif
#if ( defined TRACERS)
         do ntr=1,ntracers
         ele =  tracxc(lm1,ntr) -  tracxc(l,ntr)
         tracxc(lm1,ntr) = tracxc(lm1,ntr)+ele*xxm1
         tracxc(l,ntr)=tracxc(l,ntr)-ele*xxL
       enddo
#endif


  190 CONTINUE
c     DO 181 L=1,LM   
c     AJL(J,L,55)=AJL(J,L,55)+(TL(L)-TSAV(L))*SP
c     AJL(J,L,56)=AJL(J,L,56)+(QL(L)-QSAV(L))*SP
c 181 CONTINUE
  195 CONTINUE                                                          3158.31
c     CONDL=.true.
c 824 if(CONDL) go to 871
      DO 180 L=1,LM                                                     3158.32
      TSAV(L)  =TL(L)                                                   3158.5
      QSAV(L)=QL(L)
         AJ13(L)=0.                                                     3159.
         AJ50(L)=0.                                                     3160.
C**** MOISTURE (SPECIFIC HUMIDITY)                                      3161.
      QL(L)=QL(L)                                                       3162.
      XCD(L)=0.                                                         3163.
      DFMX(L)=0.                                                        3164.
      TO(L)=TL(L)                                                       3165.
      QO(L)=QL(L)                                                       3166.
      CXCD(L)=0.                                                        3167.
  180 CONTINUE                                                          3168.
C**** INDICES FOR WINDS                                                 3169.
      ID(1)=I+(J-1)*IM                                                  3172.
      ID(2)=ID(1)+IM*JM*LM                                              3173.
      ID(3)=I+J*IM                                                      3176.
      ID(4)=ID(3)+IM*JM*LM                                              3177.
C**** DETERMINE LATENT HEAT OF EVAPORATION OR SUBLIMATION               3178.
      TPREC(I,J)=TL(1)-TF                                               3179.
      ELHX=LHE                                                          3179.5
      IF (TPREC(I,J ).LT.0.) ELHX=LHS                                   3179.51
      CLH=ELHX/SHA                                                      3179.52
      BXQSAT=ELHX*BXCONS                                                3179.53
      AXQSAT=AXCONS+BXQSAT/TF                                           3179.54
      GAMFAC=CLH*BXQSAT                                                 3179.55
C****                                                                   3180.
C**** CONVECTION AND CLOUDS                                             3181.
C****                                                                   3182.
         HCNDNS=0.                                                      3183.
         CMC=0.                                                         3184.
         DEPTH=0.                                                       3185.
C**** INITIALIZE CONVECTION PARAMETERS                                  3186.
      QSURF=BLDATA(I,J,3)                                               3187.
      DO 225 L=1,LSSM                                                   3188.
         AJ8(L)=0.                                                      3189.
      SHSAT(L)=QSA1(TL(L),PL(L))                                        3189.5
      FMXA(L)=0.                                                        3190.
      IF(POLE) GO TO 222                                                3191.
      DO 220 K=1,4                                                      3192.
  220 UL(K,L)=UC(ID(K),1,L)                                             3193.
      GO TO 225                                                         3194.
  222 DO 223 IPO=1,IM                                                   3195.
      UPL(IPO,L)=UC(IPO,JVPO,L)                                         3196.
  223 VPL(IPO,L)=VC(IPO,JVPO,L)                                         3197.
  225 CONTINUE                                                          3198.
  232 PRCPMC=0.                                                         3199.
      DO 235 L=1,LMCMM1                                                 3201.
      LCOND=L                                                           3202.
      IF (SHSAT(LCOND).LT.QSURF) GO TO 238                              3203.
  235 CONTINUE                                                          3204.
  238 CONTINUE   
      prnt=j.eq.35
      prnt=.false.
      DO 370 LB=LCOND,LMCMM1                                            3205.
      DTCRIT=1.8                               
      SUMTT=0.                                                          3207.
      SUMQT=0.                                                          3208.
      SUMFMX=0.                                                         3209.
      EXPTUP=0.                                                         3210.
      QTCOND=0.                                                         3211.
      FCL=0.                                                            3212.

#if ( defined CPL_CHEM )
!
      sumcfc11=0.0
      sumcfc12=0.0
      sumxn2o =0.0
      sumo3   =0.0
      sumco   =0.0
      sumzco2 =0.0
      sumxno  =0.0
      sumxno2 =0.0
      sumxn2o5=0.0
      sumhno3 =0.0
      sumch4  =0.0
      sumch2o =0.0
      sumso2  =0.0
      sumh2so4=0.0

#ifdef INC_3GASES
! === 032698
      sumhfc134a = 0.0
      sumpfc     = 0.0
      sumsf6     = 0.0
#endif

      sumbc   = 0.0
      sumoc   = 0.0

! 062295
!      sumh2o2 =0.0
!
#endif
#if ( defined TRACERS)
         do ntr=1,ntracers
           tracsum(ntr)=0.0
       enddo
#endif

C**** DIFFERENCES IN STATIC ENERGY AND PRELIMINARY CLOUD TEMPERATURES   3213.
      DSE(LB)=0.                                                        3214.
      DRYSTE=0.                                                         3215.
      PDNK=PLK(LB)                                                      3216.
      SIGT=SIGMA2(LB)                                                   3217.
      TCLA(LB)=0.                                                       3218.
      BYSIGT=1./(SIGT+1.E-10)                                           3219.
      if(prnt)then
        print *,' TAU=',TAU
        print *,' LB=',LB,SIGT,BYSIGT
      endif
      DO 240 L=LB,LMCMM1                                                3220.
      DPHI=(PHIT(J,L)-PHIT(J,LB))/GRAV                                  3221.
      BYTEM=BYSIGT                                                      3222.
C     IF(DPHI.LT..5*(DWV(J)+DWV(J+1))*COEKD) BYTEM=1.E10                3223.
      PUPK=PLK(L+1)                                                     3224.
      THEDGE=THBAR(TH(L+1),TH(L))                                       3225.
      DRYSTE=DRYSTE+(TH(L+1)-THEDGE)*PUPK+(THEDGE-TH(L))*PDNK           3226.
      DSE(L+1)=(DRYSTE+CLH*(SHSAT(L+1)-SHL(LB)))*BYTEM                  3227.
      IF (DSE(L+1).LT.DSE(L)) DSE(L+1)=DSE(L)                           3228.
!     if(prnt)then
!        print 'i4,3f10.4',l,TH(L+1),THEDGE,TH(L)
!        print 'i4,2f10.4',l,DSE(L+1),DRYSTE*BYTEM
!     endif
      TCLA(L+1)=TCLA(L)-TH(L+1)*(PDNK-PUPK)                             3229.
  240 PDNK=PUPK                                                         3230.
      if(prnt)then
        do l=LB,LMCM
         print '(i4,3f10.4)',l,SIG(l),TH(l)*PLK(L),DSE(L)
        enddo
      endif
      L=LMCMM1+2                                                        3231.
  245 L=L-1                                                             3232.
C**** COMPARE STATIC ENERGIES TO FIND CRITICAL TEMPERATURE DEVIATION    3233.
C**** AND RISING MASS (FMX)                                             3234.
      TLOLD=TL(L)                                                       3235.
      SHLOLD=SHL(L)                                                     3236.

#if ( defined CPL_CHEM )
!
      cfc11old=xcfc11(l)
      cfc11cld=xcfc11(lb)

      cfc12old=xcfc12(l)
      cfc12cld=xcfc12(lb)

      xn2oold =xxn2o (l)
      xn2ocld =xxn2o (lb)

      o3old   =xo3   (l)
      o3cld   =xo3   (lb)

      coold   =xco   (l)
      cocld   =xco   (lb)

      zco2old =xzco2 (l)
      zco2cld =xzco2 (lb)

      xnoold  =xxno  (l)
      xnocld  =xxno  (lb)

      xno2old =xxno2  (l)
      xno2cld =xxno2  (lb)

      xn2o5old=xxn2o5(l)
      xn2o5cld=xxn2o5(lb)

      hno3old =xhno3 (l)
      hno3cld =xhno3 (lb)

      ch4old  =xch4  (l)
      ch4cld  =xch4  (lb)

      ch2oold =xch2o (l)
      ch2ocld =xch2o (lb)

      so2old  =xso2  (l)
      so2cld  =xso2  (lb)

      h2so4old=xh2so4(l)
      h2so4cld=xh2so4(lb)

        ! === if hfc, pfc, and sf6 are included:
#ifdef INC_3GASES
	! === 032698
      hfc134aold=xhfc134a(l)
      hfc134acld=xhfc134a(lb)

      pfcold=xpfc(l)
      pfccld=xpfc(lb)

      sf6old=xsf6(l)
      sf6cld=xsf6(lb)
	! ===
#endif

      bcold =xbc(l)
      bccld =xbc(lb)
      
      ocold =xoc(l)
      occld =xoc(lb)

c 062295
c      h2o2old =xh2o2(l)
c      h2o2cld =xh2o2(lb)

!
#endif
#if ( defined TRACERS)
         do ntr=1,ntracers
           tracold(ntr)=tracxc(l,ntr)
           traccld(ntr)=tracxc(lb,ntr)
       enddo
#endif


      DIFFSE=DSE(L)                                                     3237.
      FMX=0.                                                            3238.
      QCOND=0.                                                          3239.
      DSIGDN=DSIG(L)                                                    3240.
      RM=DSIG(LB)/DSIGDN                                                3241.
      BYRM=1./RM                                                        3242.
      CUTOFF=0.                                                         3243.
      IF (RM.GT.1.01) CUTOFF=1.5+BYRM*(1.096*BYRM-2.596)                3244.
C**** CUTOFF RESTRICTS MASS EXCHANGE TO 50% OF THE THINNER LAYER        3245.
      IF (DIFFSE.LT.CUTOFF) DIFFSE=CUTOFF                               3246.
      if(prnt)then
        print *,'L=',L,DSE(l),DTCRIT
      endif
      IF (DTCRIT.LE.DIFFSE+.005) GO TO 269                              3247.
C     FMX =.5-.5*ERF(DIFFSE)-SUMFMX                                     3248.
      FMX=ERFCPI(DIFFSE)-SUMFMX                                         3249.
      DTCRIT=DIFFSE                                                     3250.
      if(prnt)then
        print *,'L=',l,' FMX=',FMX
      endif
C**** DETERMINE CLOUD TEMPERATURE BEFORE CONDENSATION                   3251.
      EXPTDN=EXP(-DTCRIT*DTCRIT)                                        3252.
      DSTEN=.2881*SIGT/FMX*(EXPTDN-EXPTUP)                              3253.
      QWT=0.                                                            3254.
      TWT=1.-QWT                                                        3255.
      DTCL=(TL(LB)-TLOLD)+ TCLA(L) + DSTEN*TWT                          3256.
      SUMTT=SUMTT+FMX*(TL(LB)+DSTEN*TWT)                                3257.
      EXPTUP=EXPTDN                                                     3258.

#if ( defined CPL_CHEM )
!
! --- Accumulated total amount of Tracers:
!
      sumcfc11=sumcfc11+fmx*cfc11cld
      sumcfc12=sumcfc12+fmx*cfc12cld
      sumxn2o =sumxn2o +fmx*xn2ocld
      sumo3   =sumo3   +fmx*o3cld
      sumco   =sumco   +fmx*cocld
      sumzco2 =sumzco2 +fmx*zco2cld
      sumxno  =sumxno  +fmx*xnocld
      sumxno2 =sumxno2 +fmx*xno2cld
      sumxn2o5=sumxn2o5+fmx*xn2o5cld
      sumhno3 =sumhno3 +fmx*hno3cld
      sumch4  =sumch4  +fmx*ch4cld
      sumch2o =sumch2o +fmx*ch2ocld
      sumso2  =sumso2  +fmx*so2cld
      sumh2so4=sumh2so4+fmx*h2so4cld

        ! === if hfc, pfc, and sf6 are included:
#ifdef INC_3GASES
	! === 032698
      sumhfc134a  = sumhfc134a  + fmx*hfc134acld
      sumpfc      = sumpfc      + fmx*pfccld
      sumsf6      = sumsf6      + fmx*sf6cld
	! ===
#endif

      sumbc =sumbc + fmx*bccld
      sumoc =sumoc + fmx*occld

c 062295
c      sumh2o2 =sumh2o2+fmx*h2o2cld
c
!
#endif
#if ( defined TRACERS)
         do ntr=1,ntracers
           tracsum(ntr)=tracsum(ntr)+fmx*traccld(ntr)
       enddo
#endif

C**** FIND CONDENSATION AND CLOUD TEMPERATURE                           3259.
      QCLOUD=SHL(LB)+QWT*DSTEN/CLH                                      3260.
      SUMQT=SUMQT+FMX*QCLOUD                                            3261.
      TCL=TLOLD+DTCL                                                    3262.
      IF (QCLOUD.LE.QSA1(TCL,PL(L))) GO TO 266                          3263.

!#if ( defined CPL_CHEM )
#if ( defined CPL_CHEM || defined TRACERS)
!
      dqtotal = 0.0
!
#endif

      DO 265 N=1,3                                                      3264.
      QSATCL=QSA1(TCL,PL(L))                                            3265.
      GAMA=GAMFAC*QSATCL/(TCL*TCL)                                      3266.
      DQCOND=(QCLOUD-QSATCL)/(1.+GAMA)                                  3267.
      TCL=TCL+CLH*DQCOND                                                3268.

!#if ( defined CPL_CHEM )
#if ( defined CPL_CHEM || defined TRACERS)
!
! --- 062195
!
      dqtotal = dqtotal + dqcond
!
#endif

      QCOND =QCOND +DQCOND                                              3269.
  265 QCLOUD=QCLOUD-DQCOND                                              3270.

#if ( defined CPL_CHEM || defined TRACERS)

	! === convective precipitation
	prec_cnv(l) = dqtotal

#if ( defined CPL_CHEM )
!
! --- 062195:
!  Calculate scavenging of gases by convection
!  -- assume pH is around 5.0
!
! let n(v) & s(VI) disolved almost completely
!     by using a large Henry's Law constant:
!
! 020196:
!
! === Updated Henry's Law constants
!       Chien Wang 102011
!       assume pH = 5.6 for natural rainwater
!       all are actually effective Henry's Law coefficients
!       i.e., 2nd and 3rd coefficients are included
!

      ehenryx = ehenry (1.0e8,0.0,tl(l))


      h2so4cld= h2so4cld
#if ( defined OLD_CHEM )
     &        /(1.0 + raq2gas(1.e10, tl(l), dqtotal) )
#else
     &        /(1.0 + raq2gas(ehenryx, tl(l), dqtotal) )
#endif
      ! 2.1e5*(1 + 15.4*exp(8700/298)/[H+])     ! updated
      ehenryx = ehenry(1.54e1,8700.0,tl(l))
      ehenryx = 2.1e5*(1.0 + ehenryx*10.0**(5.6))
      ehenryx = ehenry (ehenryx,8700.0,tl(l))

      hno3cld = hno3cld
#if ( defined OLD_CHEM )
     &        /(1.0 + raq2gas(1.e10, tl(l), dqtotal) )
#else
     &        /(1.0 + raq2gas(ehenryx, tl(l), dqtotal) )
#endif


#if ( defined OLD_CHEM )
      ehenryx = ehenry (6.3e3,6412.34,tl(l))
#else
      ehenryx = ehenry (3.2e3,6800.00,tl(l))    ! effcoe updated
#endif
      ch2ocld = ch2ocld
     &        /(1.0 + raq2gas(ehenryx, tl(l), dqtotal) )

      ehenryx = ehenry (1.23e3,3120.00,tl(l))
      so2cld  = so2cld
     &        /(1.0 + raq2gas(ehenryx, tl(l), dqtotal) )

!
! === Note: This calculation is not counted!
! Calculate H2O2 also:

#if ( defined OLD_CHEM )
      ehenryx = ehenry (7.45,6620.00,tl(l))
#else
      ehenryx = ehenry (7.45e4,6620.00,tl(l))   ! corrected
#endif
      xh2o2(l)= xh2o2(l)
     &        /(1.0 + raq2gas(ehenryx, tl(l), dqtotal) )

!
! === For radicals apply direct reduction to gaseous phase
!	since convective transport is not involved
! Calculate HO, 062895:

      ehenryx = ehenry (25.0,5280.00,tl(l))
      ho(i,j,l)= ho(i,j,l)
     &         /(1.0 + raq2gas(ehenryx, tl(l), dqtotal) )

! Calculate HO2, 062895:

#if ( defined OLD_CHEM )
      ehenryx = ehenry (1.0e4,6640.00,tl(l))    !2nd reaction = 4.0
#else
      ehenryx = ehenry (2.3e3,6640.00,tl(l))    !2nd reaction = 4.0
#endif
      ho2(i,j,l)= ho2(i,j,l)
     &          /(1.0 + raq2gas(ehenryx, tl(l), dqtotal) )

!
#endif
#endif

      DTCL=DTCL+CLH*QCOND                                               3271.
      FCL=FCL+FMX                                                       3272.
  266 CONTINUE                                                          3273.
  269 SNWFMX=SUMFMX+FMX                                                 3274.
         AJ8(L-1)=AJ8(L-1)+SNWFMX*DSIG(LB)                              3275.
C**** REEVAPORATE WATER CONDENSED IN HIGHER LAYERS                      3276.
  270 IF (SNWFMX.EQ.0.) GO TO 303                                       3277.
      QCONDR=0.                                                         3278.
      IF (SUMFMX.EQ.0.) GO TO 280                                       3279.
      QREEV=QTCOND/SUMFMX                                               3280.
      TFALL=THUP*PLK(L)                                                 3281.
      TDN=TLOLD-CLH*QREEV                                               3282.
      SHDN=SHLOLD+QREEV                                                 3283.
      IF (SHDN.LE.QSA1(TDN,PL(L)))GO TO 280                             3284.
      QCX=SHLOLD                                                        3285.
      TCX=TLOLD                                                         3286.
      QCONDR=QREEV                                                      3287.
      QREEV=0.                                                          3288.
      DO 275 N=1,3                                                      3289.
      QSATCL=QSA1(TCX,PL(L))                                            3290.
      GAMA=GAMFAC*QSATCL/(TCX*TCX)                                      3291.
      DQREEV=(QSATCL-QCX)/(GAMA+1.)                                     3292.
      TCX=TCX-CLH*DQREEV                                                3293.
      QREEV=QREEV+DQREEV                                                3294.
      QCX=QCX+DQREEV                                                    3295.
  275 CONTINUE                                                          3296.
      QCONDR=QCONDR-QREEV                                               3297.
  280 QTCOND=QCOND*FMX+QCONDR*SUMFMX                                    3298.
C**** MIX T,Q,U,TC IN LAYER L                                           3299.
      SHDN=SHLOLD                                                       3300.
      DSH=RM*(FMX*(QCLOUD-SHDN)+SUMFMX*(SHUP+QREEV    -SHDN))           3301.
      SHL(L)=SHDN+DSH                                                   3302.

#if ( defined CPL_CHEM || defined TRACERS)
!
        xrm1   = rm*fmx
        xrm2   = rm*sumfmx
#endif
#if ( defined CPL_CHEM )
      xcfc11(l)=
     &           xrm1*(cfc11cld-cfc11old)
     &          +xrm2*(cfc11up -cfc11old)
     &          +cfc11old

      xcfc12(l)=
     &           xrm1*(cfc12cld-cfc12old)
     &          +xrm2*(cfc12up -cfc12old)
     &          +cfc12old

      xxn2o(l)=
     &           xrm1*(xn2ocld-xn2oold)
     &          +xrm2*(xn2oup -xn2oold)
     &          +xn2oold

      xo3(l)=
     &           xrm1*(o3cld-o3old)
     &          +xrm2*(o3up -o3old)
     &          +o3old

      xco(l)=
     &           xrm1*(cocld-coold)
     &          +xrm2*(coup -coold)
     &          +coold

      xzco2(l)=
     &           xrm1*(zco2cld-zco2old)
     &          +xrm2*(zco2up -zco2old)
     &          +zco2old

      xxno(l)=
     &           xrm1*(xnocld-xnoold)
     &          +xrm2*(xnoup -xnoold)
     &          +xnoold

      xxno2(l)=
     &           xrm1*(xno2cld-xno2old)
     &          +xrm2*(xno2up -xno2old)
     &          +xno2old

      xxn2o5(l)=
     &           xrm1*(xn2o5cld-xn2o5old)
     &          +xrm2*(xn2o5up -xn2o5old)
     &          +xn2o5old

      xhno3(l)=
     &           xrm1*(hno3cld-hno3old)
     &          +xrm2*(hno3up -hno3old)
     &          +hno3old

      xch4(l)=
     &           xrm1*(ch4cld-ch4old)
     &          +xrm2*(ch4up -ch4old)
     &          +ch4old

      xch2o(l)=
     &           xrm1*(ch2ocld-ch2oold)
     &          +xrm2*(ch2oup -ch2oold)
     &          +ch2oold

      xso2(l)=
     &           xrm1*(so2cld-so2old)
     &          +xrm2*(so2up -so2old)
     &          +so2old

      xh2so4(l)=
     &           xrm1*(h2so4cld-h2so4old)
     &          +xrm2*(h2so4up -h2so4old)
     &          +h2so4old

        ! === if hfc, pfc, and sf6 are included:
#ifdef INC_3GASES
	! === 032698:
      xhfc134a(l)=
     &           xrm1*(hfc134acld-hfc134aold)
     &          +xrm2*(hfc134aup -hfc134aold)
     &          +hfc134aold

      xpfc(l)=
     &           xrm1*(pfccld-pfcold)
     &          +xrm2*(pfcup -pfcold)
     &          +pfcold

      xsf6(l)=
     &           xrm1*(sf6cld-sf6old)
     &          +xrm2*(sf6up -sf6old)
     &          +sf6old
	! ===
#endif

      xbc(l)=
     &  	xrm1*(bccld-bcold)
     &         +xrm2*(bcup -bcold)
     &         +bcold

      xoc(l)=
     &  	xrm1*(occld-ocold)
     &         +xrm2*(ocup -ocold)
     &         +ocold

c 062295
c      xh2o2(l)=
c     &          xrm1*(h2o2cld-h2o2old)
c     &         +xrm2*(h2o2up -h2o2old)
c     &         +h2o2old

!
#endif
#if ( defined TRACERS)
         do ntr=1,ntracers
           tracxc(l,ntr)=xrm1*(traccld(ntr)-tracold(ntr))
     &         +xrm2*(tracup(ntr)-tracold(ntr))
     &         +tracold(ntr)
       enddo
#endif

      THDN=TH(L)                                                        3303.
      DTL=RM*(FMX*DTCL+SUMFMX*(TFALL-TLOLD-CLH*QREEV))                  3304.
      TL(L)=TLOLD+DTL                                                   3305.
      TH(L)=TL(L)/PLK(L)                                                3306.
      SHSAT(L)=QSA1(TL(L),PL(L))                                        3307.
      IF(POLE) GO TO 287                                                3308.
      DO 285 K=1,4                                                      3309.
      UDN   =UL(K,L)                                                    3310.
      UL(K,L)=UL(K,L)+RM*RA(K)*(FMX*(UL(K,LB)-UDN)+SUMFMX*(UUP(K)-UDN)) 3311.
  285 UUP(K)=UDN                                                        3312.
      GO TO 290                                                         3313.
  287 CONTINUE
      DO 288 IPOLE=1,IM                                                 3314.
      UPDN=UPL(IPOLE,L)                                                 3315.
      VPDN=VPL(IPOLE,L)                                                 3316.
      UPL(IPOLE,L)=UPL(IPOLE,L)+RM*RA(1)*(FMX*(UPL(IPOLE,LB)-UPDN)+     3317.
     *  SUMFMX*(UPUP(IPOLE)-UPDN))                                      3318.
      VPL(IPOLE,L)=VPL(IPOLE,L)+RM*RA(1)*(FMX*(VPL(IPOLE,LB)-VPDN)+     3319.
     *  SUMFMX*(VPUP(IPOLE)-VPDN))                                      3320.
      UPUP(IPOLE)=UPDN                                                  3321.
  288 VPUP(IPOLE)=VPDN                                                  3322.
  290 IF(NTRACE.EQ.0) GO TO 295                                         3323.
      DO 293 K=1,NTRACE                                                 3324.
      TCDN=SHL(L+K*39)                                                  3325.
      SHL(L+K*39)=TCDN+RM*(FMX*(SHL(LB+K*39)-TCDN)+SUMFMX*              3326.
     *  (TCUP(K)-TCDN))                                                 3327.
  293 TCUP(K)=TCDN                                                      3328.
  295 CONTINUE                                                          3329.
      SUMFMX=SNWFMX                                                     3330.
      FMXA(L)=FMXA(L)+FCL*DSIG(LB)                                      3331.
      CLDMC(I,J,L)=FMXA(L)*BX                                           3332.
      IF (CLDMC(I,J,L).LT.0.) CLDMC(I,J,L)=0.
c     IF (CLDMC(I,J,L).LT.0.005) CLDMC(I,J,L)=0.005
      IF (CLDMC(I,J,L).GT.1.) CLDMC(I,J,L)=1.                           3333.
#if ( defined HR_DATA )
       if(L.le.4)then
        cmcyzhr(L,J)=CLDMC(I,J,L)
       endif
#endif
      THUP=THDN                                                         3334.
      SHUP=SHDN                                                         3335.
      DSIGUP=DSIGDN                                                     3336.

#if ( defined CPL_CHEM )
!
      cfc11up = cfc11old
      cfc12up = cfc12old
      xn2oup  = xn2oold
      o3up    = o3old
      coup    = coold
      zco2up  = zco2old
      xnoup   = xnoold
      xno2up  = xno2old
      xn2o5up = xn2o5old
      hno3up  = hno3old
      ch4up   = ch4old
      ch2oup  = ch2oold
      so2up   = so2old
      h2so4up = h2so4old

        ! === if hfc, pfc, and sf6 are included:
#ifdef INC_3GASES
	! === 032698
      hfc134aup = hfc134aold
      pfcup     = pfcold
      sf6up     = sf6old
	! ===
#endif

      bcup =bcold
      ocup =ocold

! 062295
!      h2o2up =h2o2old
!
#endif
#if ( defined TRACERS)
         do ntr=1,ntracers
           tracup(ntr)=tracold(ntr)
       enddo
#endif


  303 IF (L.GT.LB+1) GO TO 245                                          3337.
      IF (L.EQ.LB) GO TO 355                                            3338.
      L=LB                                                              3339.
      RM=1.                                                             3340.
      FMX=0.                                                            3341.
      FCL=0.                                                            3342.
      DSIGDN=DSIG(LB)                                                   3343.
      TLOLD = (TL(LB)-SUMTT)/(1.-SUMFMX)                                3344.
      SHLOLD=(SHL(LB)-SUMQT)/(1.-SUMFMX)                                3345.

#if ( defined CPL_CHEM || defined TRACERS )
!
      xhaha   = 1./(1.-sumfmx)
#endif
#if ( defined CPL_CHEM )
      cfc11old=(xcfc11(lb)-sumcfc11)*xhaha
      cfc12old=(xcfc12(lb)-sumcfc12)*xhaha
      xn2oold =(xxn2o(lb)-sumxn2o)  *xhaha
      o3old   =(xo3(lb)-sumo3)      *xhaha
      coold   =(xco(lb)-sumco)      *xhaha
      zco2old =(xzco2(lb)-sumzco2)  *xhaha
      xnoold  =(xxno(lb)-sumxno)    *xhaha
      xno2old =(xxno2(lb)-sumxno2)  *xhaha
      xn2o5old=(xxn2o5(lb)-sumxn2o5)*xhaha
      hno3old =(xhno3(lb)-sumhno3)  *xhaha
      ch4old  =(xch4(lb)-sumch4)    *xhaha
      ch2oold =(xch2o(lb)-sumch2o)  *xhaha
      so2old  =(xso2(lb)-sumso2)    *xhaha
      h2so4old=(xh2so4(lb)-sumh2so4)*xhaha

        ! === if hfc, pfc, and sf6 are included:
#ifdef INC_3GASES
	! === 032698
      hfc134aold = (xhfc134a(lb)-sumhfc134a)
     &		  *xhaha
      pfcold     = (xpfc(lb)-sumpfc)
     &		  *xhaha
      sf6old     = (xsf6(lb)-sumsf6)
     &		  *xhaha
	! ===
#endif

      bcold =(xbc(lb)-sumbc)  *xhaha
      ocold =(xoc(lb)-sumoc)  *xhaha

c 062295
c      h2o2old =(xh2o2(lb)-sumh2o2)  *xhaha
c
!
#endif
#if ( defined TRACERS)
         do ntr=1,ntracers
           tracold(ntr)=(tracxc(l,ntr)-tracsum(ntr))*xhaha
       enddo
#endif

      GO TO 270                                                         3346.
  355 CONTINUE                                                          3347.
      PRCPMC=PRCPMC+QTCOND*DSIGDN*SP                                    3348.
  370 CONTINUE                                                          3349.
  216 DO 215 L=1,LTM                                                    3409.
      DTL=TL(L)-TSAV(L)                                                 3410.
         HCNDNS=HCNDNS+DTL*DSIG(L)                                      3411.
         AJL(J,L,13)=AJL(J,L,13)+DTL*SP                                 3412.
         AJL(J,L,57)=AJL(J,L,57)+(QL(L)-QSAV(L))*SP
         IF(J.GE.11.AND.J.LE.13) AIL(I,L,6)=AIL(I,L,6)+DTL*SP*DXYP(J)   3414.
         AJL(J,L,8)=AJL(J,L,8)+AJ8(L)*SP                                3415.
      IF (POLE) GO TO 205                                               3416.
      DO 200 K=1,4                                                      3417.
  200 U(ID(K),1,L)=U(ID(K),1,L)+(UL(K,L)-UC(ID(K),1,L))                 3418.
      GO TO 215                                                         3419.
  205 DO 210 IPO=1,IM                                                   3420.
      U(IPO,JVPO,L)=U(IPO,JVPO,L)+(UPL(IPO,L)-UC(IPO,JVPO,L))           3421.
  210 V(IPO,JVPO,L)=V(IPO,JVPO,L)+(VPL(IPO,L)-VC(IPO,JVPO,L))           3422.
  215 CONTINUE                                                          3423.
      if(HPRNT)then
        print *,' condse 4'
      print *,' T(J,L)=',T(1,JPR,LPR),' Q(J,L)=',Q(1,JPR,LPR)
      print *,' V(J,L)=',V(1,JPR,LPR),' V(J+1,L)=',V(1,JPR+1,LPR)
      print *,' U(J,L)=',U(1,JPR,LPR),' U(J+1,L)=',U(1,JPR+1,LPR)
      endif
      do 873 L=1,LM
      if(PCLOUD.eq.1)then
         CLDMC(I,J,L)=CMDATA(J,L)*.01                                   3471.2
      elseif(PCLOUD.eq.6)then
         CLDMC(I,J,L)=CLDMCT(J,L)
      endif
  873 continue
c     go to 872
  871 CONTINUE
C****                                                                   3424.
C**** LARGE SCALE PRECIPITATION                                         3425.
C****                                                                   3426.
      PRCPSS=0.                                                         3427.
         CSS=0.                                                         3429.
      DQUP=0.                                                           3430.
      ELHXUP=LHE                                                        3431.
      DO 304 LX=1,LM                                                    3432.
      L=LM+1-LX                                                         3433.
      TOLD=TL(L)                                                        3434.
      QOLD=QL(L)                                                        3435.
      ELHX= LHE                                                         3436.
      IF(TOLD.LT.TI) ELHX= LHS                                          3437.
      IF (ELHXUP.EQ.LHS.AND.TOLD.LT.TF) ELHX=LHS                        3438.
      EX=DQUP*DSIGUP/DSIG(L)                                            3439.
      TNEW=TOLD-CLH*EX                                                  3440.
      QNEW=QOLD+EX                                                      3441.
      DQUP=0.                                                           3442.
      QSATL=QSAT(TNEW,PL(L))                                            3443.
      ELHXUP=LHE                                                        3444.
C**** DETERMINE THE CLOUD COVER                                         3445.
CC**  IF (QNEW.LE.1.E-10) GO TO 300                                     3446.
      RHLL=QNEW/QSATL                                                   3446.1
      AJL(J,L,58)=AJL(J,L,58)+RHLL*SP
c     AJL(J,L,59)=AJL(J,L,59)+(RHLL*SP)**2
      IF (QNEW.LE.1.E-10) GO TO 300
        RH0=RHKP(L,j)
      if(HPRNT)then
      if(L.eq.2)then
        print *,' condse CLDSS TAU=',TAU
      print *,TNEW,PL(L),QSATL
      print *,' RHLL=',RHLL,' RH0=',RH0
      endif
      endif
        if(RHLL.gt.RH0)then
         CLDSS(I,J,L)=(RHLL-RH0)/(1.-RH0)  ! 2353.05
!        CLDSS(I,J,L)=((RHLL-RH0)/(1.-RH0) )**2
        else
         CLDSS(I,J,L)=0.
        endif
      if(PCLOUD.eq.5)then
        RH0=RH0OLD
        CLDSS(I,J,L)=CSCALE*(RHLL-RH0)/(1.-RH0)                         3446.2
        IF(PL(L).LT.400.) CLDSS(I,J,L)=.4166667*CLDSS(I,J,L)            3446.21
      endif
      if(PCLOUD.eq.1)then
         CLDSS(I,J,L)=CSDATA(J,L)*.01                                   3471.1
      elseif(PCLOUD.eq.6)then
         CLDSS(I,J,L)=CLDSST(J,L)
      endif
#if ( defined HR_DATA )
       if(L.le.4)then
        pyzhr(L,J)=PL(L)
        tyzhr(L,J)=TL(L)
        rhyzhr(L,J)=RHLL
        cssyzhr(L,J)=CLDSS(I,J,L)
       endif
#endif
      IF(CLDSS(I,J,L).GT.1.) CLDSS(I,J,L)=1.                            3446.3
      IF(CLDSS(I,J,L).LT.0.0) CLDSS(I,J,L)=0.0                          3446.4
c     IF(CLDSS(I,J,L).LT.0.005) CLDSS(I,J,L)=0.005
  300 IF (QNEW.LT.RHNEW(j)*QSATL) GO TO 302                                3455.
      ELHX=LHE                                                          3456.
      IF (TOLD.LT.TF) ELHX=LHS                                          3457.
C     RHNEW=1.                                                          3458.
      CLH=ELHX/SHA                                                      3459.
      GAMFAC=CLH*BXCONS*ELHX                                            3460.

!#if ( defined CPL_CHEM )
#if ( defined CPL_CHEM || defined TRACERS)
!
	dqtotal = 0.0
!
#endif

      DO 301 N=1,3                                                      3461.
      GAMA=GAMFAC*QSATL/(TNEW*TNEW)                                     3462.
      DQ1=(QNEW-QSATL*RHNEW(j))/(1.+GAMA*RHNEW(j))                      3463.

!#if ( defined CPL_CHEM )
#if ( defined CPL_CHEM || defined TRACERS)
!
	dqtotal = dqtotal + dq1
!
#endif

      DQUP=DQUP+DQ1                                                     3464.
      TNEW=TNEW+CLH*DQ1                                                 3465.
      QNEW=QNEW-DQ1                                                     3466.
      QSATL=QSAT(TNEW,PL(L))                                            3467.
  301 CONTINUE

#if ( defined CPL_CHEM || defined TRACERS)

	! === stratform precipitation:
	prec_str(l) = dqtotal

#if ( defined CPL_CHEM )
!
! --- 062195:
!  Calculate scavenging of gases by large-scale
!     precipitation
!  -- assume pH is around 5.0
!
! let n(v) & s(VI) disolved almost completely
!     by using a large Henry's Law constant:
!
!
! === Updated Henry's Law constants
!       Chien Wang 102011
!       assume pH = 5.6 for natural rainwater
!       all are actually effective Henry's Law coefficients
!       i.e., 2nd and 3rd coefficients are included
!

      ehenryx = ehenry (1.0e8,0.0,tl(l))
!
! 020196:

      xh2so4(l)= xh2so4(l)
#if ( defined OLD_CHEM )
     &        /(1.0 + raq2gas(1.e10, tl(l), dqtotal) )
#else
     &        /(1.0 + raq2gas(ehenryx, tl(l), dqtotal) )
#endif
      ! 2.1e5*(1 + 15.4*exp(8700/298)/[H+])     ! updated
      ehenryx = ehenry(1.54e1,8700.0,tl(l))
      ehenryx = 2.1e5*(1.0 + ehenryx*10.0**(5.6))
      ehenryx = ehenry (ehenryx,8700.0,tl(l))

      xhno3(l) = xhno3(l)
#if ( defined OLD_CHEM )
     &        /(1.0 + raq2gas(1.e10, tl(l), dqtotal) )
#else
     &        /(1.0 + raq2gas(ehenryx, tl(l), dqtotal) )
#endif

#if ( defined OLD_CHEM )
      ehenryx = ehenry (6.3e3,6412.34,tl(l))
#else
      ehenryx = ehenry (3.2e3,6800.00,tl(l))    ! effcoe updated
#endif
      xch2o(l) = xch2o(l)
     &         /(1.0 + raq2gas(ehenryx, tl(l), dqtotal) )

      ehenryx  = ehenry (1.23e3,3120.00,tl(l))
      xso2(l)  = xso2(l)
     &         /(1.0 + raq2gas(ehenryx, tl(l), dqtotal) )

! Calculate H2O2 also:

#if ( defined OLD_CHEM )
      ehenryx = ehenry (7.45,6620.00,tl(l))
#else
      ehenryx = ehenry (7.45e4,6620.00,tl(l))   ! corrected
#endif
      xh2o2(l) = xh2o2(l)
     &         /(1.0 + raq2gas(ehenryx, tl(l), dqtotal) )

! Calculate HO, 062895:

      ehenryx = ehenry (25.0,5280.00,tl(l))
      ho(i,j,l)= ho(i,j,l)
     &         /(1.0 + raq2gas(ehenryx, tl(l), dqtotal) )

! Calculate HO2, 062895:

#if ( defined OLD_CHEM )
      ehenryx = ehenry (1.0e4,6640.00,tl(l))    !2nd reaction = 4.0
#else
      ehenryx = ehenry (2.3e3,6640.00,tl(l))    !2nd reaction = 4.0
#endif
      ho2(i,j,l)= ho2(i,j,l)
     &          /(1.0 + raq2gas(ehenryx, tl(l), dqtotal) )

!
#endif
#endif

      DSIGUP=DSIG(L)                                                    3468.
      ELHXUP=ELHX                                                       3469.
  302 TL(L)=TNEW                                                        3470.
      QL(L)=QNEW                                                        3471.
C**** ACCUMULATE SOME DIAGNOSTICS                                       3472.
         HCNDNS=HCNDNS+(TNEW-TOLD)*DSIG(L)                              3473.
  304    AJL(J,L,11)=AJL(J,L,11)+(TNEW-TOLD)*SP                         3474.
      PRCPSS=DQUP*DSIG(1)*SP                                            3475.
c     CONDL=.FALSE.
c     GO TO 824
c 872 CONTINUE
!    07/22/2005 different precipitation over land and ocean
!       PRLAND and PROCEAN are ratios of precip
!        over land and ocean to total precipitation
!
         AJ(J,61)=AJ(J,61)+PRCPSS*POCEAN*PROCEAN                        3476.
         BJ(J,61)=BJ(J,61)+PRCPSS*PLAND*PRLAND                          3477.
         CJ(J,61)=CJ(J,61)+PRCPSS*POICE*PROCEAN                         3478.
         DJ(JR,61)=DJ(JR,61)+PRCPSS*DXYP(J)                             3479.
  305    AJ(J,62)=AJ(J,62)+PRCPMC*POCEAN*PROCEAN                        3480.
         BJ(J,62)=BJ(J,62)+PRCPMC*PLAND*PRLAND                          3481.
         CJ(J,62)=CJ(J,62)+PRCPMC*POICE*PROCEAN                         3482.
         DJ(JR,62)=DJ(JR,62)+PRCPMC*DXYP(J)                             3483.
         DO 390 KR=1,4                                                  3484.
         IF(I.EQ.IJD6(1,KR).AND.J.EQ.IJD6(2,KR)) GO TO 392              3485.
  390    CONTINUE                                                       3486.
      GO TO 400                                                         3487.
  392    ADAILY(IHOUR,5,KR)=ADAILY(IHOUR,5,KR)+HCNDNS*SP                3488.
         ADAILY(IHOUR,49,KR)=ADAILY(IHOUR,49,KR)+PRCPMC+PRCPSS          3489.
  400 PRCP=(PRCPMC+PRCPSS)*100./GRAV                                    3490.
      PREC(I,J)=PRCP                                                    3491.
      IF(TPREC(I,J).GE.0.) PRCP=0.                                      3492.
      GDATA(I,J,11)=(DTPERD+GDATA(I,J,11)*AGESNX)*EXP(-PRCP)            3493.
C**** TOTAL HEATING AND MOISTURE ADJUSTMENT                             3494.
  500 DO 530 L=1,LM                                                     3495.
      T(I,J,L)=TL(L)/PLK(L)                                             3496.

#if ( defined CPL_CHEM )
!
      cfc11(i,j,l)= xcfc11(l)
      cfc12(i,j,l)= xcfc12(l)
      xn2o (i,j,l)= xxn2o (l)
      o3   (i,j,l)= xo3   (l)
      co   (i,j,l)= xco   (l)
      zco2 (i,j,l)= xzco2 (l)
      xno  (i,j,l)= xxno  (l)
      xno2 (i,j,l)= xxno2 (l)
      xn2o5(i,j,l)= xxn2o5(l)
      hno3 (i,j,l)= xhno3 (l)
      ch4  (i,j,l)= xch4  (l)
      ch2o (i,j,l)= xch2o (l)
      so2  (i,j,l)= xso2  (l)
      h2so4(i,j,l)= xh2so4(l)

        ! === if hfc, pfc, and sf6 are included:
#ifdef INC_3GASES
	! === 032698
      hfc134a(i,j,l) = xhfc134a(l)
      pfc(i,j,l)     = xpfc(l)
      sf6(i,j,l)     = xsf6(l)
	! ===
#endif

      bcarbon (i,j,l)= xbc (l)
      ocarbon (i,j,l)= xoc (l)

! 062295
      h2o2 (i,j,l)= xh2o2 (l)
!
#endif
#if ( defined TRACERS)
         do ntr=1,ntracers
           tracers(i,j,l,ntr) = tracxc(l,ntr)
       enddo
#endif


  530 Q(I,J,L)=QL(L)                                                    3497.

!070804
#if ( defined CPL_CHEM || defined TRACERS)
!
	beta = 3600.0*0.15		! dt*correction
!beta = 1.0
        beta = 3600.0*0.05
	t_cnv = max(0.0, prec_cnv(nlev))
	t_str = max(0.0, prec_str(nlev))
	do k=nlev-1,1,-1
	
	   ! === accumulate precipitation
	   !t_cnv = t_cnv + prec_cnv(k)
	   !t_str = t_str + prec_str(k)

	   if ( TX(i,j,k) .le. 273.15 ) then	! Eice =0.35
	     t_cnv = prec_cnv(k)*0.5
	     t_str = prec_str(k)*0.5
	   else
	     t_cnv = prec_cnv(k)
	     t_str = prec_str(k)
	   end if
	   
#if ( defined CPL_CHEM )
	   ! === Wet scavenging by convective precipiation:
	   bcarbon(i,j,k) = bcarbon(i,j,k)
     &			  *(1.0 - 4.4913e-2*t_cnv*beta)
	   if ( bcarbon(i,j,k) .lt. 0.0 ) bcarbon(i,j,k) = 0.0
	   ocarbon(i,j,k) = ocarbon(i,j,k)
     &		  	  *(1.0 - 4.4913e-2*t_cnv*beta)
	   if ( ocarbon(i,j,k) .lt. 0.0 ) ocarbon(i,j,k) = 0.0
	   
	   ! === Wet scavenging by large scale precipitation:
	   bcarbon(i,j,k) = bcarbon(i,j,k)
     &			  *(1.0 - 5.3946e-2*t_str*beta)
	   if ( bcarbon(i,j,k) .lt. 0.0 ) bcarbon(i,j,k) = 0.0
	   ocarbon(i,j,k) = ocarbon(i,j,k)
     &			  *(1.0 - 5.3946e-2*t_str*beta)
	   if ( ocarbon(i,j,k) .lt. 0.0 ) ocarbon(i,j,k) = 0.0
#endif

#if ( defined TRACERS)
! betatr = 1.0
      	 betatr = 3600.0*0.15		! dt*correction
         do ntr=1,ntracers
          if (trtype(ntr).eq.4) then
	   ! === Wet scavenging by convective precipiation:
	   tracers(i,j,k,ntr) = tracers(i,j,k,ntr)
     &			  *(1.0 - 4.4913e-2*t_cnv*betatr)
	   if ( tracers(i,j,k,ntr) .lt. 0.0 ) tracers(i,j,k,ntr) = 0.0
	   ! === Wet scavenging by large scale precipitation:
	   tracers(i,j,k,ntr) = tracers(i,j,k,ntr)
     &			  *(1.0 - 5.3946e-2*t_str*betatr)
	   if ( tracers(i,j,k,ntr) .lt. 0.0 ) tracers(i,j,k,ntr) = 0.0
         endif
       enddo
#endif

	end do

!070804
#endif

  700 IM1=I                                                             3498.
#if ( defined CLM )

      pred4tem(j)=pred4tem(j)+PREC(1,J)
      ewvd4tem(j)=ewvd4tem(j)+QL(1)*P(1,j)*SIG(1)*RVAP/RGAS
      npred4tem(j)=npred4tem(j)+1

c     prhr(j)=PREC(1,J)
c     PRCP=(PRCPMC+PRCPSS)*100./GRAV                                    3490.
      i=1
      pcpl4clm(i,j)=PRCPSS*100./GRAV
      pcpc4clm(i,j)=PRCPMC*100./GRAV
      tpr4clm(i,j)=TPREC(1,J)
#endif
C
#if ( defined OCEAN_3D || defined ML_2D )
      tempr(j)=tempr(j)+TPREC(1,J)
      precip(j)=precip(j)+PREC(1,J)
      if(j.eq.-42)then
       print *,'FROM CONDSE'
       print *,'TPREC=',TPREC(1,J),' PREC=',PREC(1,J)
      endif
      ps4ocean(j)=ps4ocean(j)+(SP+PTOP)
      do l=1,lm
       qyz4ocean(j,l)=qyz4ocean(j,l)+QL(l)
       tyz4ocean(j,l)=tyz4ocean(j,l)+TL(l)
      enddo
#endif
c

C**** END OF MAIN LOOP FOR I INDEX                                      3499.
  810 CONTINUE                                                          3500.
C****                                                                   3501.
C**** END OF MAIN LOOP FOR J INDEX                                      3502.
C****                                                                   3503.
C**** ADD IN CHANGE OF ANG. MOMENTUM BY MOIST CONVECTION FOR DIAGNOSTIC 3504.
         DO 880 L=1,LTM                                                 3505.
         DO 880 J=2,JM                                                  3506.
         DO 880 I=1,IM                                                  3507.
  880    AJL(J,L,39)=AJL(J,L,39)+(U(I,J,L)-UC(I,J,L))*P(I,J)            3508.
      JDAY00=JDAY
      if(HPRNT)then
        print *,' condse 6'
      print *,' T(J,L)=',T(1,JPR,LPR),' Q(J,L)=',Q(1,JPR,LPR)
      print *,' V(J,L)=',V(1,JPR,LPR),' V(J+1,L)=',V(1,JPR+1,LPR)
      print *,' U(J,L)=',U(1,JPR,LPR),' U(J+1,L)=',U(1,JPR+1,LPR)
      endif
C

#if ( defined CPL_CHEM )
!
! --- Chemistry model patch 081795
!  check negative values:
!
      call chemcheck(cfc11)
      call chemcheck(cfc12)
      call chemcheck(xn2o )
      call chemcheck(o3   )
      call chemcheck(co   )
      call chemcheck(zco2 )
      call chemcheck(xno  )
      call chemcheck(xno2 )
      call chemcheck(xn2o5)
      call chemcheck(hno3 )
      call chemcheck(ch4  )
      call chemcheck(ch2o )
      call chemcheck(so2  )
      call chemcheck(h2so4)
      call chemcheck(h2o2 )
      call chemcheck(bcarbon)
      call chemcheck(ocarbon)

        ! === if hfc, pfc, and sf6 are included:
#ifdef INC_3GASES
	! === 032698
      call chemcheck(hfc134a)
      call chemcheck(pfc)
      call chemcheck(sf6)
	! ===
#endif
!
#endif

      RETURN                                                            3509.
      END                                                               3510.
