#include "ctrparam.h"

!       ==========================================================
!
!       EDDYPA_NEW.F: THIS ROUTINE DEALS WITH THE PARAMETERIZATION OF
!                 EDDY TRANSPORTS AND VARIANCES
!       Modified version, includes real averaging
!
!
      SUBROUTINE EDDYPA (UT,VT,TT,PT,QT,U,V,T,P,Q,DTE,EDDFLX)           2181.   
C                                                                       2181.5  
C     THIS ROUTINE DEALS WITH THE PARAMETERIZATION OF                   2182.   
C     EDDY TRANSPORTS AND VARIANCES                                     2182.5  
C                                                                       2183.   
#if ( defined CPL_CHEM || defined TRACERS )
!
#include "chem_para"
#include "chem_com"
!
#endif

#include "BD2G04.COM"

      COMMON/WORK3/PHI(IM0,JM0,LM0)                                     2184.   
      COMMON/WORK4/FD(IM0,JM0)                                          2184.5  
      DIMENSION UT(IM0,JM0,LM0),VT(IM0,JM0,LM0),TT(IM0,JM0,LM0),
     &   PT(IM0,JM0),
     *  QT(IM0,JM0,LM0)                                                 2185.5  
      DIMENSION SINL(JM0),COSL(JM0),                                    2186.   
     *  DXU(JM0),DYU(JM0),DXYU(JM0),PHIS(IO0,JM0),UX(IO0,JM0,2*LM0)     2186.5  
      EQUIVALENCE (SINL(1),SINP(1)),(COSL(1),COSP(1)),                  2187.   
     *  (DXU(1),DXV(1)),(DYU(1),DYV(1)),(DXYU(1),DXYV(1)),              2187.5  
     *  (PHIS,FDATA)                                                    2188.   
c     REAL*4 KAPAP1                                                     2188.5  
      REAL KAPAP1,MOMCOEF
      COMMON/EPARA/VTH(JM0,LM0),WTH(JM0,LM0),VU(JM0,LM0),VV(JM0,LM0)
     & ,DQSDT(JM0,LM0)                                                  2189.   
     * ,DWV(JM0),PHIT(JM0,LM0),TPRIM2(JM0,LM0),WU(JM0,LM0),CKS,CKN      2189.5  
     * ,WQ(JM0,LM0),VQ(JM0,LM0)                                         2190.   
	DIMENSION PHITAV(LM0),TPRIM0(JM0)
      COMMON/SPEC2/KM,KINC,COEK,C3LAND(IO0,JM0),C3OICE(IO0,JM0)         2190.5  
     * ,C3LICE(IO0,JM0),WMGE(IO0,JM0),TSSFC(1,JM0,4)                    2191.   
      DIMENSION TSZM(JM0),THSZM(JM0),FO(JM0),WUDATA(JM0,LM0 )
     & ,VUDATA(JM0,LM0),                                                2191.5  
     &  VSH(JM0,LM0),WSH(JM0,LM0),TEMZ(JM0,LM0),COEKY(JM0)
     & ,DTHDY(JM0,LM0),TEMJL(JM0,LM0+1)                                 2192.   
     * ,DSHDY(JM0,LM0),DSHDP(JM0,LM0+1),DTHDZ(JM0,LM0+1)
     & ,DSHDZ(JM0,LM0+1)                                                2192.5  
     * ,SHSAT(JM0,LM0),DSHDYS(JM0,LM0)                                  2193.   
      DIMENSION DQDYJ(JM0),DQDY(JM0,LM0),UZM(JM0),PSROK(JM0,LM0)
     & ,DQSDT1(JM0,LM0),                                                2193.5  
     *  DUDP(JM0,LM0),TZA(JM0,LM0),DTHDP(JM0,LM0+1),PRESR(JM0,LM0+1)
     & ,PSRK(JM0,LM0+1)                                                 2194.   
     * ,PSI(JM0,LM0),PVOR(JM0,LM0),VQZ(JM0),UNEW(JM0,LM0),UOLD(JM0,LM0)
     & ,FMIX(JM0,LM0)                                                   2194.5  
     * ,SCLH(JM0,LM0),GAMA(JM0),SHJ(JM0),HTJ(JM0),ZSTAR(JM0,LM0)
     & ,BVF(JM0)                                                        2195.   
     * ,PRSRU(JM0,LM0),DTDYJ(JM0),DTDZJ(JM0),DWVJL(JM0,LM0)
     & ,BVFJL(JM0,LM0),                                                 2195.5  
     *  DUDZ(JM0,LM0),DUDZJ(JM0),COEJL(JM0,LM0),DQDY1(JM0,LM0)
     & ,COEJL1(JM0,LM0)                                                 2196.   
     * ,UVAR(JM0),TEMJ(JM0),WTDATA(JM0,LM0),WQDATA(JM0,LM0)
     & ,VTDATA(JM0,LM0),                                                2196.5  
     *  VQDATA(JM0,LM0),VVDATA(JM0,LM0),TTDATA(JM0,LM0),GAMAJ(JM0)
     & ,DTHDZJ(JM0)                                                     2197.   
     * ,BVFM(JM0,LM0),BVFMJ(JM0),DQDZ(JM0,LM0),VLAT(JM0),VVN(9)
     & ,VVS(9)   
     & ,VVNEWA(LM0),FOV(JM0),BETAV(JM0),BETAP(JM0),TEMKGP(JM0)
     & ,TEMKGV(JM0),RDWVERR(JM0),ADWVERR(JM0)
     & ,TEMKV(JM0),TEMKP(JM0),GAMAP(JM0),DWVP(JM0),HTJP(JM0)
     & ,DWVA(JM0),DWVPA(JM0)
      dimension TZN(JM0,LM0),QZN(JM0,LM0),UZN(JM0,LM0),VZN(JM0,LM0)
      common/tprmtg/tprmg(JM0),ntprmg(JM0)
      common/Dscale/DWAV0(JM0)
      logical first,uvtr,wvtr,uttr,wttr,uqtr,wqtr,EDDFLX
     &,OLDUVVERT,HPRNT
      SAVE
      data first/.true./
      DATA VVN/37.6,39.5,57.8,94.8,150.6,236.3,329.5,114.8,43.8/        2261.5  
      DATA VVS/32.7,28.3,33.7,44.9,79.4,151.5,193.6,90.9,15.7/          2262.   
C                                                                       2262.5  
      QSAT(TM,PR,QL)=3.797915*EXP(QL*(7.93252E-6-2.166847E-3/TM))/PR    2263.   
C                                                                       2263.5  
c     HPRNT=TAU.gt.50.0.and.TAU.lt.75.0
      HPRNT=.FALSE.
c     HPRNT=TAU.ge.17520.0.and.TAU.lt.17521.0
c     HPRNT=.TRUE.
      JPR=24
      uvtr=EDDFLX
      wvtr=EDDFLX
      uttr=EDDFLX
      wttr=EDDFLX
      uqtr=EDDFLX
      wqtr=EDDFLX
      if(first)then
      LZ=LM                                                             2264.   
      LZM1=LZ-1                                                         2264.5  
      LZM2=LZ-2                                                         2265.   
      JHALF=JM/2                                                        2265.5  
      JHAFM1=JHALF-1                                                    2266.   
      JHAFP2=JHALF+2                                                    2266.5  
      JMM2=JM-2                                                         2267.   
      IX=IM                                                             2267.5  
      IS=IM                                                             2268.   
      FIS=IS                                                            2268.5  
      CKS=1.                                                            2268.6  
      CKN=1.                                                            2268.7  
      CP=RGAS/KAPA                                                      2269.5  
      SHA=RGAS/KAPA                                                     2270.   
      CONSTA=2.0
      CONSTA=1.5
      CONSTA=1.0
      RCZZ=1.35
c     RCZZ=1.83
      MOMCOEF=1.
c     MOMCOEF=4.*0.5*(0.57+0.71)/CONSTA
      KAPAP1=KAPA+1.                                                    2271.   
      HLAT=LHE                                                          2271.5  
      CLH=HLAT/CP                                                       2272.   
      DQDTX=.622*HLAT/RGAS                                              2272.5  
      P1000=1000.**KAPA                                                 2273.   
      DTDYMN=1.E-6/P1000                                                2273.5  
      DQDYMN=1.E-09                                                     2274.   
      DSCALE=626.*1.e3
c     DSCALE=323.*1.e3
c     DSCALE=1000.*1.e3
      FDR=exp(-DLAT*RADIUS/DSCALE)
      OLDUVVERT=.false.
      DWVCRI=1000.
      DWVCRI=100.
      JINT=6
      IF(JM.EQ.46)JINT=12
c     JINT=JM/2-2
      JEQ=JM/2+1
      JEQP1=JEQ+1   
      JEQM1=JEQ-1  
      JEQM2=JEQ-2 
      LBGN=1     
      ALPHA=1.  
      NVAIT=2
      NVAIT=3
      WMGMIN=3.0
      WMGMIN=0.0
        print *,' Mom. flux param. for dry atmosph.'
        print *,' New vertical avereging'
        print *,' Change in YCUT calculations'
        print *,' with tg difference'
        print *,' CONSTA=',CONSTA
        print *,' RCZZ=',RCZZ
        print *,' DWVCRI=',DWVCRI
        print *,' MRCH=',MRCH
        print *,' FDR=',FDR,' L=',DSCALE
        print *,' JINT=',JINT
        print *,' NVAIT=',NVAIT
        print *,' PSF=',PSF
        print *,' WMGMIN=',WMGMIN
        print *,' With A for momentum=',CONSTA*MOMCOEF
        print *,' UVTR=',uvtr,' WVTR=',wvtr
        print *,' UTTR=',uttr,' WTTR=',wttr
        print *,' UQTR=',uqtr,' WQTR=',wqtr
        if(OLDUVVERT)then
	  print *,' Vertical structure of UV flux from GISS
     &     2D model'
	endif
      DO 3 J=1,JM                                                       2274.5  
  3   FO(J)=F(J)/DXYP(J)                                                2275.   
      DO 4 J=2,JM                                                       2275.5  
         VLAT(J)=(LAT(J)+LAT(J-1))*.5                                   2276.   
         FOV(J)=0.5*(FO(J)+FO(J-1))
         BETAV(J)=(FO(J)-FO(J-1))/DYU(J)
         if(j.NE.JM)BETAP(J)=.5*(FO(J+1)-FO(J-1))/DYU(J)
       TEMY=RADIUS*(.25*TWOPI-ABS(VLAT(J)))
       Xk=1./(RADIUS*COSV(j))
       Yk=TWOPI/(4.*TEMY)
       TEMKGV(j)=SQRT(Xk**2+Yk**2)
  4   CONTINUE
      TPR00=1.0
      TPR90=2.5
      TPR11=0.5*(TPR90+TPR00)
      TPR12=0.5*(TPR90-TPR00)
      print *,'TPR00=',TPR00,' TPR90=',TPR90
      do j=1,jm
         TEMY=RADIUS*(.25*TWOPI-ABS(LAT(J)))
         Xk=1./(RADIUS*COSP(j))
         Yk=TWOPI/(4.*TEMY)
         TEMKGP(j)=SQRT(Xk**2+Yk**2)
         tprimrad = 3.14159*(float(j-1)/22.5)
         TPRIM0(j)=5.5+(4.5*cos(tprimrad)) ! 2323
         TPRIM0(j)=0.0 ! 2324
!        TPRIM0(j)=1.5+(1.0*cos(tprimrad))	! 2327 2332
         TPRIM0(j)=TPR11+(TPR12*cos(tprimrad))	! 2327 2332
!        TPRIM0(j)=2.0 ! 2331
      enddo
c      print *,'LAT'
c      print *,LAT
c      print *,'VLAT'
c      print *,VLAT
c      print *,'COSV'
c      print *,COSV
       print *,'TPRIM0'
       print *,TPRIM0
        first=.false.
      endif
C                                                                       2276.5  
#if ( defined CPL_CHEM || defined TRACERS )
!
      do 1 j=1,nlat
        beta1(j)=0.0
        beta2(j)=0.0
      do 2 k=1,nlev
        beta3(j,k)=0.0
        beta4(j,k)=0.0
        fkt  (j,k)=0.0
2     continue
1     continue
!
#endif

C reassign T,Q,U and V to internal variables

      do l=1,lm
       do j=1,jm
        TZN(j,l)=T(1,j,l)
        QZN(j,l)=Q(1,j,l)
        UZN(j,l)=U(1,j,l)
        VZN(j,l)=V(1,j,l)
       enddo
      enddo
      if(HPRNT)then
       print *,' Eddypa UZN  before  filter'
       print *,'LZM2=',LZM2
         print *,(UZN(13,L),L=1,LM)
         print *,(UZN(14,L),L=1,LM)
       endif
      if(JM.eq.46)then
       CALL FLTR(TZN,JM,LM,1)
       CALL FLTR(QZN,JM,LM,1)
       CALL FLTR(UZN,JM,LM,2)
      endif

      DO 6 J=1,JM                                                       2277.   
      DWV(J)=0.                                                         2277.5  
      DO 6 I=1,IM                                                       2278.   
  6   WMGE(I,J)=WMGMIN                                                  2278.5  
      DO 10 L=1,LM                                                      2279.   
      DO 10 J=1,JM                                                      2279.5  
      VU(J,L)=0.                                                        2280.   
      VQ(J,L)=0.                                                        2280.5  
      VTH(J,L)=0.                                                       2281.   
      WTH(J,L)=0.                                                       2281.5  
      WU(J,L)=0.                                                        2282.   
      WQ(J,L)=0.                                                        2282.5  
      VV(J,L)=0.                                                        2283.   
!
!CAS    Applying latitudinal depedence of convection params
!
!     tprimrad = 3.14159*(float(j-1)/22.5)
!     TPRIMCO=6.0+(4.0*cos(tprimrad))
      TPRIM2(J,L)=TPRIM0(j)/(EXPBYK(SIG(L)*P(1,J)+PTOP))**2. 
!CAS      TPRIM2(J,L)=2./(EXPBYK(SIG(L)*P(1,J)+PTOP))**2.                   2283.5  
!     AJL(J,L,59)=AJL(J,L,59)+2.*PSF
      AJL(J,L,59)=AJL(J,L,59)+TPRIM0(j)*PSF
c     TPRIM2(J,L)=0
      PHIT(J,L)=PHI(1,J,L)                                              2284.   
 10   CONTINUE                                                          2284.5  
      DO L=1,LM
	PHITAV(L)=0.
	CCC=0.
	DO J=1,JM
         PHITAV(L)=PHITAV(L)+PHIT(J,L)*COSP(J)
         CCC=CCC+COSP(J)
	ENDDO
	PHITAV(L)=PHITAV(L)/CCC
      ENDDO
C                                                                       2285.   
C     COMPUTE PRESSURE, P**KAPA, STATIC STABILITY AND SCALE HEIGHT      2285.5  
C                                                                       2286.   
      DO 24 L=1,LZ                                                      2286.5  
      DO 24 J=1,JM                                                      2287.   
      PRESR(J,L+1)= P(1,J)*SIGE(L+1)+PTOP                               2287.5  
      PSRK(J,L+1)=PRESR(J,L+1)**KAPA                                    2288.   
      IF(J.GE.2) PRSRU(J,L)=.5*(P(1,J)+P(1,J-1))*SIG(L)+PTOP            2288.5  
      IF(J.GE.2) PSROK(J,L)=PRSRU(J,L)**KAPA                            2289.   
      TEM=( P(1,J)*SIG(L)+PTOP)**KAPA                                   2289.5  
c24   TZA(J,L)= T(1,J,L)*TEM                                            2290.   
 24   TZA(J,L)= TZN(J,L)*TEM
      DO 23 J=1,JM                                                      2290.5  
      DO 22 L=1,LZM1                                                    2291.   
c     TEM= T(1,J,L+1)- T(1,J,L)                                         2291.5  
      TEM= TZN(J,L+1)- TZN(J,L)
      IF(TEM.LT..1) TEM=.1                                              2292.   
      DTHDP(J,L)=-TEM/(P(1,J)*DSIGO(L))                                 2292.5  
c     DSHDP(J,L)=(Q(1,J,L)-Q(1,J,L+1))/(P(1,J)*DSIGO(L))                2293.   
      DSHDP(J,L)=(QZN(J,L)-QZN(J,L+1))/(P(1,J)*DSIGO(L))
      TEMJL(J,L+1)=-TEM/(SIG(L+1)-SIG(L))*RGAS*                         2293.5  
     *  PSRK(J,L+1)/PRESR(J,L+1)/ P(1,J)                                2294.   
      TEMX=GRAV*PRESR(J,L+1)/(.5*RGAS*(TZA(J,L)+TZA(J,L+1)))            2294.5  

#if ( defined CPL_CHEM || defined TRACERS )
!
      deltap(j,l)=1./(p(1,j)*dsigo(l))
      dp2dz (j,l)=-temx
!
#endif


      DTHDZ(J,L+1)=-DTHDP(J,L)*TEMX                                     2295.   
      DSHDZ(J,L+1)=-DSHDP(J,L)*TEMX                                     2295.5  
 22   CONTINUE                                                          2296.   
      TEMJL(J,1)=(TEMJL(J,2)*(DSIG(1)+DSIG(2))-TEMJL(J,3)*DSIG(1))/     2296.5  
     *  DSIG(2)                                                         2297.   
      IF(TEMJL(J,1).LE..1*TEMJL(J,2)) TEMJL(J,1)=.1*TEMJL(J,2)          2297.5  
      TEMJL(J,LZ+1)=TEMJL(J,LZ)                                         2298.   
      PRESR(J,1)= P(1,J)+PTOP                                           2298.5  
      PSRK(J,1)=PRESR(J,1)**KAPA                                        2299.   
      DTHDZ(J,1)=DTHDZ(J,2)                                             2299.5  
      DSHDZ(J,1)=DSHDZ(J,2)                                             2300.   
      DTHDZ(J,LZ+1)=DTHDZ(J,LZ)                                         2300.5  
      DSHDZ(J,LZ+1)=DSHDZ(J,LZ)                                         2301.   
 23   CONTINUE                                                          2301.5  
      DO 27 L=1,LZ                                                      2302.   
      DO 26 J=2,JMM1                                                    2302.5  
 26   TEMZ(J,L)=0.25*(TEMJL(J,L)*2.+TEMJL(J+1,L)+TEMJL(J-1,L))          2303.   
      TEMJL(1,L)=(2.*TEMJL(1,L)+TEMJL(2,L))/3.                          2303.5  
      TEMJL(JM,L)=(2.*TEMJL(JM,L)+TEMJL(JMM1,L))/3.                     2304.   
      DO 27 J=2,JMM1                                                    2304.5  
      TEMJL(J,L)=TEMZ(J,L)                                              2305.   
 27   CONTINUE                                                          2305.5  
      DO 21 L=1,LZM1                                                    2306.   
      DO 21 J=1,JM                                                      2306.5  
      TEM=2.*PRESR(J,L+1)*GRAV/RGAS/(TZA(J,L)+TZA(J,L+1))               2307.   
      BVFJL(J,L)=TEM*SQRT(TEMJL(J,L+1))                                 2307.5  
C    DRY
c     BVFM(J,L)=BVFJL(J,L)*BVFJL(J,L)+2.*DSHDZ(J,L+1)*CLH*GRAV/         2308.   
c    *          (TZA(J,L)+TZA(J,L+1))                                   2308.5  
      BVFM(J,L)=BVFJL(J,L)*BVFJL(J,L)
c     RVFM is now (Nd)**2 not (Nm)**2
C     DRY
      SCLH(J,L)=RGAS*.5*(TZA(J,L)+TZA(J,L+1))/GRAV                      2309.   
 21   CONTINUE                                                          2309.5  
C                                                                       2310.   
                                                                        2310.5  
C     COMPUTE DEL QS/DEL T                                              2311.   
C                                                                       2311.5  
      DO 30 L=1,LZM1                                                    2312.   
      DO 30 J=1,JM                                                      2312.5  
      TEMZ(J,L)=0.                                                      2313.   
      DQSDT1(J,L)=0.                                                    2313.5  
      DO 25 I=1,IS                                                      2314.   
      PR=(P(I,J)*SIG(L)+PTOP)                                           2314.5  
      TM=TZA(J,L)                                                       2315.   
      TEM=QSAT(TM,PR,HLAT)                                              2315.5  
      SHSAT(J,L)=TEM                                                    2316.   
c     TEM=DQDTX/TM/TM*Q(1,J,L)                                          2316.5  
c     TEM1=DQDTX*Q(1,J,L)/TM/TM                                         2317.   
      TEM=DQDTX/TM/TM*QZN(J,L)
      TEM1=DQDTX*QZN(J,L)/TM/TM
      TEMZ(J,L)=TEMZ(J,L)+TEM                                           2317.5  
      DQSDT1(J,L)=DQSDT1(J,L)+TEM1                                      2318.   
 25   CONTINUE                                                          2318.5  
      DQSDT1(J,L)=DQSDT1(J,L)*CLH/FIS                                   2319.   
 30   DQSDT(J,L)=TEMZ(J,L)/FIS                                          2319.5  
C                                                                       2320.   
C     COMPUTE MERIDIONAL TEMPERATURE AND MOISTURE GRADIENT              2320.5  
C                                                                       2321.   
      DO 3140 L=1,LZ                                                    2321.5  
      DO 3140 J=2,JM                                                    2322.   
c     DTHDY(J,L)=( T(1,J,L)- T(1,J-1,L))/DYU(J)                         2322.5  
      DTHDY(J,L)=( TZN(J,L)- TZN(J-1,L))/DYU(J)
 3140 CONTINUE                                                          2323.   
      DO 3141 L=1,LZ                                                    2323.5  
      DO 3142 J=3,JMM1                                                  2324.   
 3142 TEMZ(J,L)=.25*(2.*DTHDY(J,L)+DTHDY(J+1,L)+DTHDY(J-1,L))           2324.5  
      DTHDY(2,L)=(2.*DTHDY(2,L)+DTHDY(3,L))/3.                          2325.   
      DTHDY(JM,L)=(2.*DTHDY(JM,L)+DTHDY(JMM1,L))/3.                     2325.5  
      DO 3141 J=3,JMM1                                                  2326.   
 3141 DTHDY(J,L)=TEMZ(J,L)                                              2326.5  
      DO 3144 L=1,LZ                                                    2327.   
      DO 3144 J=2,JM                                                    2327.5  
c     DSHDY(J,L)=( Q(1,J,L)- Q(1,J-1,L))/DYU(J)                         2328.   
      DSHDY(J,L)=( QZN(J,L)- QZN(J-1,L))/DYU(J)
 3144 CONTINUE                                                          2328.5  
      DO 3148 L=1,LZ                                                    2329.   
      DO 3148 J=2,JM                                                    2329.5  
      DSHDYS(J,L)=(SHSAT(J,L)-SHSAT(J-1,L))/DYU(J)                      2330.   
 3148 CONTINUE                                                          2330.5  
C     DO 3147 L=1,LZ                                                    2331.   
C     DO 3146 J=3,JMM1                                                  2331.5  
C3146 TEMZ(J,L)=.25*(2.*DSHDY(J,L)+DSHDY(J+1,L)+DSHDY(J-1,L))           2332.   
C     DSHDY(2,L)=(2.*DSHDY(2,L)+DSHDY(3,L))/3.                          2332.5  
C     DSHDY(JM,L)=(2.*DSHDY(JM,L)+DSHDY(JMM1,L))/3.                     2333.   
C     DO 3147 J=3,JMM1                                                  2333.5  
C3147 DSHDY(J,L)=TEMZ(J,L)                                              2334.   
C                                                                       2334.5  
C     COMPUTE WAVE DEPTH                                                2335.   
C                                                                       2335.5  
      DO 3145 J=2,JM                                                    2336.5  
      B45=BETAV(J)
      FCONT=FOV(J)
      DO 3145 L=1,LZM2                                                  2338.   
      TEM1=.5*(BVFJL(J,L)+BVFJL(J-1,L))                                 2338.5  
      TEM=.5*(SCLH(J,L)+SCLH(J-1,L))*TEM1*TEM1*B45*.25*                 2339.   
c    *  (T(1,J,L)+T(1,J,L+1)+T(1,J-1,L+1)+                              2339.5  
     *  (TZN(J,L)+TZN(J,L+1)+TZN(J-1,L+1)+
c    *  T(1,J-1,L))/(FCONT*GRAV*.5*(DTHDY(J,L)+DTHDY(J,L+1))+1.E-20)    2340.   
     *  TZN(J-1,L))/(FCONT*GRAV*.5*(DTHDY(J,L)+DTHDY(J,L+1))+1.E-20)
c     DWVJL(J,L)= .5*(SCLH(J,L)+SCLH(J-1,L))/                           2340.5  
c    *  (1.48*ABS(TEM)+.48)                                             2341.   
      DWVJL(J,L)= .5*(SCLH(J,L)+SCLH(J-1,L))/ 
     &  (2./RCZZ*(1.+ABS(TEM))-1.)
      IF(DWVJL(J,L).LT.DWVCRI) DWVJL(J,L)=DWVCRI                        2341.5  
 3145 CONTINUE                                                          2342.   
C                                                                       2342.5  
C     DETERMINE THE LAT.FOR CUTTING OFF THE WAVE PROPAGATION            2343.   
C                                                                       2343.5  
      DO 3150 J=2,JM                                                    2345.   
      USUM=0.                                                           2345.5  
      DO 3151 L=1,LM                                                    2346.   
c3151 USUM=USUM+U(1,J,L)*DSIG(L)                                        2346.5  
 3151 USUM=USUM+UZN(J,L)*DSIG(L)
 3150 UZM(J)=USUM                                                       2347.   
      JB=JEQ+1
      JE=JB+JINT                                                        2348.   
      YCUTN=0.                                                          2348.5  
      DO 3155 J=JB,JE                                                   2349.   
      JR=JB+JE-J                                                        2349.5  
 3155 IF(UZM(JR).LE.0.) GO TO 3156                                      2350.   
      GO TO 3157                                                        2350.5  
 3156 YCUTN=(UZM(JR)*VLAT(JR+1)-UZM(JR+1)*VLAT(JR))/
     *      (UZM(JR)-UZM(JR+1)-1.E-20)
 3157 CONTINUE
      JE=JEQ-1
      JB=JE-JINT                                                        2352.5  
      YCUTS=0.                                                          2353.   
      DO 3158 JR=JB,JE                                                  2353.5  
 3158 IF(UZM(JR).LE.0.) GO TO 3159                                      2354.   
      GO TO 3160                                                        2354.5  
 3159 YCUTS=(UZM(JR)*VLAT(JR-1)-UZM(JR-1)*VLAT(JR))/
     *      (UZM(JR)-UZM(JR-1)-1.E-20)
 3160 CONTINUE                                                          2356.   
c     print *,'YCUTN=',180./TWOPI*YCUTN,' YCUTS=',180./TWOPI*YCUTS
C                                                                       2356.5  
C                                                                       2356.5
C     NEW WEIGHTING VARIABLES BY WAVE DEPTH        
C                                                                    
c     print *,' DWAV0'
c     print *,DWAV0
      DO J=1,JM
        DWV(J)=DWAV0(J)
      ENDDO
      DO II=1,NVAIT
      DO  J=1,JM                                                
      DO  L=1,LZM1                                               
      TEMZ(J,L)=(TZA(J,L+1)-TZA(J,L))*GRAV/(PHI(1,J,L+1)-PHI(1,J,L))
      ENDDO	! L
      ENDDO	! J
      CALL VWEIGHAV(GRAV,TEMZ,DTDZJ,PHI,PHIS,DWV,1,LZM2,1,JM  ,PTOP, 
     *  SCLH,P,SIGE,DSIGO,CKS,CKN,IO,IM,JM,LM)                      
      DO  J=2,JM                                               
      DO  L=1,LZM2                                            
      TEMZ(J,L)=.5*(DTHDY(J,L)+DTHDY(J,L+1))                     
      ENDDO	! L
      ENDDO	! J
      CALL VWEI1AV(GRAV,TEMZ ,DTDYJ,PHI,PHIS,DWV,1,LZM2,2,JM,PTOP,
     *  SCLH,P,SIGE,DSIGO,CKS,CKN,IO,IM,JM,LM)                   
      DO  L=1,LZM1                                          
      DO  J=1,JM                                           
      TEMZ(J,L)=DTHDZ(J,L+1)                                  
      ENDDO	! L
      ENDDO	! J
      CALL VWEIGHAV(GRAV,TEMZ,DTHDZJ,PHI,PHIS,DWV,1,LZM2,1,JM  ,PTOP, 
     *  SCLH,P,SIGE,DSIGO,CKS,CKN,IO,IM,JM,LM)                       
c     Vertical averaging of N**2
      CALL VWEIGHAV(GRAV,BVFM,BVFMJ,PHI,PHIS,DWV,1,LZM2,1,JM  ,PTOP, 
     *  SCLH,P,SIGE,DSIGO,CKS,CKN,IO,IM,JM,LM)                      
c     Nw=sqrt((N**2)w)
      DO j=1,JM
        BVF(j)=sqrt(BVFMJ(j))
      ENDDO
c     Nw=sqrt((N**2)w)
      CALL VWEIGHAV(GRAV,SCLH,SHJ,   PHI,PHIS,DWV,1,LZM2,1,JM  ,PTOP,  
     *  SCLH,P,SIGE,DSIGO,CKS,CKN,IO,IM,JM,LM)                        
      DO  L=1,LZM1                                               
      DO  J=1,JM                                                
c     TEMZ(J,L)=.5*(T(1,J,L)+T(1,J,L+1))                           
      TEMZ(J,L)=.5*(TZN(J,L)+TZN(J,L+1))
      ENDDO	! L
      ENDDO	! J
      CALL VWEIGHAV(GRAV,TEMZ,THSZM, PHI,PHIS,DWV,1,LZM2,1,JM  ,PTOP, 
     *  SCLH,P,SIGE,DSIGO,CKS,CKN,IO,IM,JM,LM)                       
      DO  L=1,LZM2                                              
      DO  J=2,JM                                               
c     TEMZ(J,L)=.5*(U(1,J,L)+U(1,J,L+1))                          
      TEMZ(J,L)=.5*(UZN(J,L)+UZN(J,L+1))
      ENDDO	! L
      ENDDO	! J
      if(HPRNT)then
       J=JM-1
       print *,'LZ=',LZ
       print *,' Eddypa UZN  before  VWEI1AV'
         print *,(UZN(13,L),L=1,LM)
         print *,(UZN(14,L),L=1,LM)
         print *,DWV(J+1),DWV(J-1),DWV(J)
       endif
      CALL VWEI1AV(GRAV,TEMZ,UZM,PHI,PHIS,DWV,1,LZM2,2,JM,       
     *  PTOP,SCLH,P,SIGE,DSIGO,CKS,CKN,IO,IM,JM,LM)             
      if(HPRNT)then
       print *,' Eddypa UZM after VWEI1AV'
         J=14
         print *,UZM(J+1),UZM(J-1),UZM(J)
       endif
      GAMAD=GRAV/CP                                            
      DO  J=1,JM                                          
      DO  L=1,LZM1                                       
      TEDGE=.5*(TZA(J,L+1)+TZA(J,L))                        
      QSE=QSAT(TEDGE,PRESR(J,L+1),LHE)                     
      TEMZ(J,L)=GAMAD*(1.+LHE*QSE/(RGAS*TEDGE))/          
     *      (1.+.622*LHE*LHE*QSE/(CP*RGAS*TEDGE*TEDGE))  
      ENDDO	! L
      ENDDO	! J
      CALL VWEIGHAV(GRAV,TEMZ,GAMAJ,PHI,PHIS,DWV,1,LZM2,1,JM  ,PTOP,  
     *  SCLH,P,SIGE,DSIGO,CKS,CKN,IO,IM,JM,LM)                       
      DO J=2,JM                                       
      B45=BETAV(J)
      FCONT=FOV(J)
c
c      BVFMJ=(N**2)w
      TEM1=.5*(BVFMJ(J)+BVFMJ(J-1))
      TEM=.5*(SHJ(J)+SHJ(J-1))*TEM1*B45* .5*(THSZM(J)+THSZM(J-1))/
     *  (FCONT*GRAV*DTDYJ(J)+1.E-20)
c
      GAMA(J)=ABS(TEM)                                
      HTJ(J)=.5*(SHJ(J)+SHJ(J-1))/(1.+GAMA(J))       
      IF(HTJ(J).LT.DWVCRI) HTJ(J)=DWVCRI            
c     DWV(J)=.5*(SHJ(J)+SHJ(J-1))/(1.48*GAMA(J)+.48) 
      DWV(J)=.5*(SHJ(J)+SHJ(J-1))/(2./RCZZ*(1.+GAMA(J))-1.)
      IF(DWV(J).LT.DWVCRI) DWV(J)=DWVCRI           
      IF(BVFMJ(J).LT.1.E-5) BVFMJ(J)=1.E-5        
C
       DWV(J)=0.5*(DWV(j)+DWAV0(J))
C
      ENDDO	! J
      IF(BVFMJ(1).LT.1.E-5) BVFMJ(1)=1.E-5      
c     print *,' NEW II=',ii
c     print *,' DWV'
c     print *,DWV
      AERR=0.0
      RERR=0.0
      AERRMAX=0.0
      RERRMAX=0.0
      DO J=2,JM
       ADWVERR(J)=abs(DWV(J)-DWAV0(J))
       RDWVERR(J)=ADWVERR(J)/DWAV0(J)
       AERR=AERR+ADWVERR(J)
       RERR=RERR+RDWVERR(J)
       DWAV0(J)=DWV(J)
       if(ADWVERR(J).gt.AERRMAX)then
            AERRMAX=ADWVERR(J)
            JAMAX=j
       endif
       if(RDWVERR(J).gt.RERRMAX)then
            RERRMAX=RDWVERR(J)
            JRMAX=j
       endif
      ENDDO
c     print *,'Eddy iteration ii=',ii
c     print *,'AERR=',AERR/(jm-1.0),' RERR=',RERR/(jm-1.0)
c     print * ,'MAXAERR=',AERRMAX,' JAMAX=',JAMAX
c     print * ,'MAXRERR=',RERRMAX,' JRMAX=',JRMAX
c     print *,' DWV  D'
c     print *,DWV
c     print *,' HTJ   d'
c     print *,HTJ
c     print *,' SHJ   H'
c     print *,SHJ
      if(RERRMAX.lt.1.e-3) go to 7788
      ENDDO	! iterations
 7788 continue
c     if(ii.gt.20)print *,' NIT=',II
      DO J=1,JM
       DWAV0(J)=DWV(J)
      ENDDO
c     print *,' after  DAAVO calculation'
c     print *,' GAMA '
c     print *,GAMA
c     print *,' DWV'
c     print *,DWV
c     print *,' HTJ'
c     print *,HTJ
C                                                                   
C      DWV (D) and CAMA on P-grid
       do j=2,jmm1
        TEM=SHJ(J)*BVFMJ(J)*BETAP(J)*THSZM(J)/
     &  (FO(J)*GRAV*0.5*(DTDYJ(J)+DTDYJ(J+1)))
        GAMAP(j)=abs(TEM)
        DWVP(J)=SHJ(J)/(2./RCZZ*(1.+GAMAP(J))-1.)
        IF(DWVP(J).LT.DWVCRI) DWVP(J)=DWVCRI           
        HTJP(J)=SHJ(J)/(1.+GAMAP(J))       
        IF(HTJP(J).LT.DWVCRI) HTJP(J)=DWVCRI            
       enddo
c     print *,' GAMAP '
c     print *,(GAMAP(j),j=2,jm-1)
c     print *,' '
c     print *,(0.5*(GAMA(j)+GAMA(J+1)),j=2,jm-1)
c     print *,' DWVP'
c     print *,(DWVP(j),j=2,jm-1)
c     print *,' '
c     print *,(0.5*(DWV(j)+DWV(J+1)),j=2,jm-1)
C       NEW
C                                                                       2387.   
c     Adjust wave length near pole
      DO j=2,JM
      if(j.ne.JEQ)then
      TEM1=.5*(BVF(J)+BVF(J-1))
      TEM=.5*(SHJ(J)+SHJ(J-1))*TEM1/ABS(FOV(J))
      TEMH=TEMKGV(j)*TEM
      TEMKV(J)=sqrt(((1.+GAMA(J))/RCZZ)**2-0.25)
C
       if(TEMKV(J).lt.TEMH)then
         DWVA(J)=.5*(SHJ(J)+SHJ(J-1))/(sqrt(4.*TEMH**2+1.)-1.)
        else
         DWVA(J)=DWV(J)
       endif
      endif
      ENDDO	! J
      do j=2,jmm1
       TEM=BVF(J)*SHJ(J)/abs(FO(J))
       TEMH=TEMKGP(j)*TEM
       TEMKP(J)=sqrt(((1.+GAMAP(J))/RCZZ)**2-0.25)
C
      if(TEMKP(J).lt.TEMH)then
        DWVPA(J)=SHJ(J)/(sqrt(4.*TEMH**2+1.)-1.)
       else
        DWVPA(J)=DWVP(J)
       endif
      enddo
c     print *,' After adjustment'
c     print *,' DWV'
c     print *,DWV
c     print *,' DWVA'
c     print *,DWVA
c     print *,' DWVP'
c     print *,DWVP
c     print *,' DWVPA'
c     print *,DWVPA
c     print *,' HTJ'
c     print *,HTJ
c     stop
C     COMPUTE TEMPERATURE VARIANCES                                     2387.5  
C                                                                       2388.   
      DO 4013 J=2,JMM1                                                  2388.5  
      COEKC=CKS                                                         2389.   
      IF(J.GT.JHALF) COEKC=CKN                                          2389.5  
c     TEM=BVF(J)*.25*(HTJ(J)+HTJ(J+1))*(DTDYJ(J)+DTDYJ(J+1))/           2390.   
      TEM=BVF(J)*.5*HTJP(J)*(DTDYJ(J)+DTDYJ(J+1))/  
     *  (FO(J)+1.E-20)                                                  2390.5  
      IF ((DTDYJ(J)+DTDYJ(J+1))*FO(J).GT.0.) TEM=0.                     2390.6  
c     IF(FO(J).EQ.0.) TEM=0.                                            2391.   
      if(abs(FO(J)).lt.1.e-15) TEM=0.
       PLAND=FDATA(1,J,2)                                               2395.   
       PWATER=1.-PLAND                                                  2395.5  
       PLICE=FDATA(1,J,3)*PLAND                                         2396.   
       PEARTH=PLAND-PLICE                                               2396.5  
       POICE=ODATA(1,J,2)*PWATER                                        2397.   
       POCEAN=PWATER-POICE                                              2397.5  
      if(POCEAN.LE.1.E-5)then
         POCEAN=0.
         POICE=PWATER
      endif
       tland=GDATA(1,J,4)
       toice=GDATA(1,J,3)
       tlice=GDATA(1,J,13)
       tocean=ODATA(1,J,1)
       TSS=POCEAN*tocean+POICE*toice+PLICE*tlice+PEARTH*tland
       TPR1=POCEAN*(tocean-TSS)**2                                
       TPR2=POICE*(toice-TSS)**2
       TPR3=PLICE*(tlice-TSS)**2                                 
       TPR4=PEARTH*(tland-TSS)**2
       tprmg(J)=tprmg(J)+(TPR1+TPR2+TPR3+TPR4)
       ntprmg(J)=ntprmg(J)+1
c     TEMH=.5*(DWV(J)+DWV(J+1))                                         2400.5  
      TEMH=DWVPA(J)
      DO 4013 L=1,LZ                                                    2400.   
      ZPHI=(PHI(1,J,L)-PHIS(1,J))/(GRAV*TEMH+1.E-20)                    2401.5  
C     TEM1=TEM*(1.-.5*ZPHI)                                             2402.   
      TEM1=TEM*(1.-.5*ZPHI)*(1.+DQSDT1(J,L))                            2402.5  
C for analisys
!     TEMAN=TEM*(1.-.5*ZPHI)
C for analisys
      TPRIM2(J,L)=TPRIM2(J,L)+1.*TEM1*TEM1*EXP(-ZPHI/COEKC)/COEKC       2403.
!     2319
!     TPRIM2(J,L)=max(TPRIM2(J,L),1.*TEM1*TEM1*EXP(-ZPHI/COEKC)/COEKC)
      PK2=(EXPBYK(SIG(L)*P(1,J)+PTOP))**2.                              2403.5  
       IF(L.EQ.1) THEN 
        TPRIM2(J,L)=TPRIM2(J,L)+.6065*                                  2404.   
     &   (TPR1+TPR2+TPR3+TPR4)/PK2                                      2404.5  
        AJL(J,L,59)=AJL(J,L,59)+.6065*
     &   (TPR1+TPR2+TPR3+TPR4)*P(1,J)
       ELSEIF(L.EQ.2) THEN 
         TPRIM2(J,L)=TPRIM2(J,L)+.1653*                                 2405.   
     &   (TPR1+TPR2+TPR3+TPR4)/PK2                                      2405.5  
        AJL(J,L,59)=AJL(J,L,59)+.1653*
     &   (TPR1+TPR2+TPR3+TPR4)*P(1,J)
      ENDIF
      AJL(J,L,59)=AJL(J,L,59)+TPRIM2(J,L)
 4013 CONTINUE                                                          2406.   
      tprmg(JM)=0.
      tprmg(1)=0.
      ntprmg(JM)=1.
      ntprmg(1)=1.
CCC  ************
C                                                                       2406.5  
C     COMPUTE WIND VARIANCE AND EDDY TRANSFER COEFFICIENTS              2407.   
C                                                                       2407.5  
      DO 4015 I=1,IM                                                    2408.   
      DO 4015 J=2,JM                                                    2408.5  
      COEKD=CKS                                                         2409.   
      IF(J.GT.JHALF) COEKD=CKN                                          2409.5  
      TEM2=FOV(J)
      TEM3=.5*(THSZM(J)+THSZM(J-1))                                     2410.5  
      TEM=HTJ(J)*DTDYJ(J)*GRAV/(TEM2*TEM3+1.E-20)                       2411.   
      IF (DTDYJ(J)*TEM2.GT.0.) TEM=0.                                   2411.1  
c     IF(TEM2.EQ.0.) TEM=0.                                             2411.5  
      if(abs(TEM2).lt.1.e-15) TEM=0.
C     TEML=3.*1.4142*.5*(BVF(J)+BVF(J-1))*DWV(J)/(TEM2+1.E-20)          2412.   
C     TEML=ABS(TEML)                                                    2412.5  
C     TEMY=RADIUS*(.25*TWOPI-.5*ABS(LAT(J)+LAT(J-1)))                   2413.   
C     RATIO=TEMY/(TEML+1.E-20)                                          2413.5  
C     IF(TEML.GT.TEMY) TEM=TEM*RATIO                                    2414.   
      UVAR(J)=.5*TEM*TEM/COEKD                                          2414.5  
      WMGETEMP=.33*TEM*TEM/COEKD
c     WMGE(I,J)=max(WMGETEMP,WMGMIN)                                     2415.   
      WMGE(I,J)=WMGETEMP+WMGMIN                                          2415.1
 4015 continue
      IF(DMOD(TAU,25.).EQ.0..AND.IABS(MRCH).EQ.1)                       2415.5  
     *  WRITE(6,905)(WMGE(1,J),J=2,JM)                                  2416.   
 905  FORMAT(' WMGE=',24F5.1)                                           2416.5  
C                                                                       2417.   
      DO 2230 J=2,JM                                                    2417.5  
      COEKC=CKS                                                         2418.   
      IF(J.GT.JHALF) COEKC=CKN                                          2418.5  
      TEM=CONSTA*(BVF(J)+BVF(J-1))*ABS(DTDYJ(J))*HTJ(J)*HTJ(J)          2419.   
     *  *GRAV/(THSZM(J)+THSZM(J-1))                                     2419.5  
      COEKY(J)=0.                                                       2420.   
      FVTH=FOV(J)
c     IF(FVTH.NE.0.) COEKY(J)=TEM/(FVTH*FVTH)/COEKC                     2421.   
      if(abs(FVTH).gt.1.e-15)COEKY(J)=TEM/(FVTH*FVTH)/COEKC
      IF(DTDYJ(J)*FVTH.GT.0.) COEKY(J)=0.                               2421.1  
c     TEMH=DWV(J)                                                       2425.5  
      TEMH=DWVA(J)                                                       2425.5  
c
C Old restriction
c
c     TEML=3.*1.4142*(BVF(J)+BVF(J-1))*DWV(J)/((FO(J)+FO(J-1))+1.E-20)  2421.5
c     TEML=ABS(TEML)                                                    2422.
c     TEMY=RADIUS*(.25*TWOPI-.5*ABS(LAT(J)+LAT(J-1)))                   2422.5
c     TEMK=.7071*.5*TWOPI*.25*(BVF(J)+BVF(J-1))*                        2423.
c    *  (SHJ(J)+SHJ(J-1))/(TEMY*ABS(FVTH)+1.E-20)                       2423.5
c     if(TEML.GT.TEMY) then
c        if(TEMK.lt.1.e19)then
c           TEMH=.5*(SHJ(J)+SHJ(J-1))/(SQRT(4.*TEMK*TEMK+1.)-1.)
c        else
c           TEMH=.5*(SHJ(J)+SHJ(J-1))/(2.*TEMK-1.)
c        end if
c     end if
c
C Old restriction
c

      DO 2230 L=1,LZM1                                                  2425.   
      ZSTAR(J,L)=.5*(PHI(1,J,L)+PHI(1,J-1,L)-PHIS(1,J)-PHIS(1,J-1))     2427.   
     *  /(GRAV*TEMH+1.E-20)                                             2427.5  
      EXPFAC=EXP(-ZSTAR(J,L)/COEKC)                                     2428.   
C     TEMV=5.+45.*(5.*SIG(L)*EXP(1.-5.*SIG(L)))**4.                     2428.5  
C    *  +TEMV                                                           2429.5  
      COEJL(J,L)=COEKY(J)*2.                                            2430.   
      COEJL1(J,L)=COEKY(J)*EXPFAC                                       2430.5  
 2230 CONTINUE                                                          2431.   
C                                                                       2431.5  
C     PARA. MERIDIONAL EDDY MOMENTUM TRANSPORT                          2432.   
C                                                                       2432.5  
      DO 35 J=2,JM                                                      2436.   
      BETA=BETAV(J)
      FAVE=FOV(J)
C  DRY
      TEM1=.5*(BVFMJ(J)+BVFMJ(J-1))                                     2437.5
c     Nd instead of Nm     Eq 23 from Part II.
c     See change in calculation of BVFM above
c     BVFMJ=(Nd)**2, while BVF=Nd
C  DRY
      TEM=FAVE*GRAV*DTDYJ(J)*2./(1.  *TEM1*(THSZM(J)+THSZM(J-1))        2438.   
     *  *DWV(J)+1.E-20)                                                 2438.5  
      if(HPRNT)then
       if(J.eq.16)then
       print *,' Eddypa 2350'
         print *,FAVE,DTDYJ(J),TEM1
         print *,THSZM(J),THSZM(J-1),DWV(J)
         print *,j,TEM
       endif
      endif
      IF (J.GT.2.AND.J.LT.JM)                                           2439.   
     *  DUDY2=(UZM(J+1)+UZM(J-1)-2.*UZM(J))/DYP(J)/DYP(J)               2439.5  
     *  -(SINL(J)+SINL(J-1))/(COSL(J)+COSL(J-1))*.5*                    2440.   
     *  (UZM(J+1)-UZM(J-1))/(RADIUS*DYU(J))                             2440.5  
     *  -UZM(J)/(RADIUS*.5*(COSL(J)+COSL(J-1)))**2.                     2441.   
      IF (J.EQ.2) DUDY2=(UZM(3)-2.*UZM(2))/DYP(2)**2.                   2441.5  
     *  -(SINL(2)+SINL(1))/(COSL(2)+COSL(1))*                           2442.   
     *  (UZM(3)-UZM(2))/(RADIUS*DYP(2))                                 2442.5  
     *  -UZM(2)/(RADIUS*.5*(COSL(2)+COSL(1)))**2.                       2443.   
      IF (J.EQ.JM) DUDY2=(UZM(JMM1)-2.*UZM(JM))/DYP(JM)**2.             2443.5  
     *  -(SINL(JMM1)+SINL(JM))/(COSL(JMM1)+COSL(JM))*                   2444.   
     *  (UZM(JM)-UZM(JMM1))/(RADIUS*DYP(JMM1))                          2444.5  
     *  -UZM(JM)/(RADIUS*.5*(COSL(JMM1)+COSL(JM)))**2.                  2445.   
      if(HPRNT)then
       if(J.eq.14)then
       print *,' Eddypa 35 J=',J
         print *,BETAV(J),DUDY2
         print *,UZM(J+1),UZM(J-1),UZM(J)
       endif
      endif
      DO 35 L=1,LZM1                                                    2445.5  
      DQDY(J,L)=  BETAV(J)-DUDY2
C  DRY
c
c35   DQDY1(J,L)=TEM*(1.+.5*(DQSDT1(J,L)+DQSDT1(J-1,L)))+DQDY(J,L)      2446.5
c
c    M=0. Eq 23 from Part II.
c
      DQDY1(J,L)=TEM+DQDY(J,L)
CCCCC  check + or *
c
C  DRY
 35   CONTINUE
      DO 2350 J=2,JM                                                    2447.   
      TEM=0.25*(COSL(J)+COSL(J-1))**2.                                  2447.5  
      DO 2380 L=1,LZM1                                                  2448.   
      TEMZ(J,L)=DQDY(J,L)*TEM                                           2448.5  
      DQDY1(J,L)=COEJL1(J,L)*DQDY1(J,L)*TEM                             2449.   
      if(HPRNT)then
       if(J.eq.14)then
       print *,' Eddypa 2350'
         print *,L,TEM
         print *,DQDY(J,L),COEJL1(J,L),DQDY1(J,L)
       endif
      endif
 2380 CONTINUE                                                          2449.5  
 2350 CONTINUE                                                          2450.   
C**** TAKING VERTICAL AVERAGES                                          2450.5  
      DO 2365 J=2,JM                                                    2451.   
      SUMSIG=0.                                                         2451.5  
      TEMJ(J)=0.                                                        2452.   
      TSZM(J)=0.                                                        2452.5  
      DO 2366 L=LBGN,LZM1                                               2453.   
      SUMSIG=SUMSIG+DSIG(L)                                             2453.5  
      TSZM(J)=TSZM(J)+DQDY1(J,L)*DSIG(L)                                2454.   
 2366  TEMJ(J)=TEMJ(J)+TEMZ(J,L)*DSIG(L)                                2454.5  
      TSZM(J)=TSZM(J)/SUMSIG                                            2455.   
      TEMJ(J) = TEMJ(J)/SUMSIG                                          2455.5  
      TSZM(J)=TSZM(J)*MOMCOEF
 2365 CONTINUE                                                          2456.   
      if(HPRNT)then
       print *,' Eddypa 2365'
         print *,'TEMJ'
         print *,TEMJ
         print *,'TSZM'
         print *,TSZM
      endif
C**** CALCULATE THE ADJUSTMENT COEFFICIENT FOR MOMENTUM BALANCING       2456.5  
      PFRAC=.0                                                          2457.   
      YSCALE=400000./RADIUS                                             2457.5  
      DO 2355 JH=1,2                                                    2458.5  
      JB=2                                                              2459.   
      JE=JM                                                             2459.5  
      IF(JH.EQ.1) JE=JEQ-1                                              2460.   
      IF(JH.EQ.2) JB=JEQ+1                                              2460.5  
      YZERO=YCUTN                                                       2461.   
      IF(JH.EQ.1) YZERO=YCUTS                                           2461.5  
      TEM1=0.                                                           2462.   
      TEM2=0.                                                           2462.5  
      DO 2353 J=JB,JE                                                   2463.   
      FKB=1.                                                            2463.5  
      DELTY=VLAT(J)-YZERO                                               2464.   
c     NLAT=ABS(DELTY/DLAT)+1.                                           2464.5  
      IF(JH.EQ.1.AND.DELTY.GT.0.)                                       2465.   
     *  FKB=EXP(-RADIUS*DELTY/DSCALE)
      IF(JH.EQ.2.AND.DELTY.LT.0.)                                       2466.   
     *  FKB=EXP(RADIUS*DELTY/DSCALE)
      SP=.5*(P(1,J)+P(1,J-1))                                           2467.   
      TEM1=TEM1+TEMJ(J)*SP*FKB                                          2467.5  
      TEM2=TEM2+TSZM(J)*SP                                              2468.   
 2353 CONTINUE                                                          2468.5  
      DO 2354 J=JB,JE                                                   2469.   
      SP=.5*(P(1,J)+P(1,J-1))                                           2469.5  
 2354 TSZM(J)=TSZM(J)-PFRAC*TEM2/(SP*(JE-JB+1))                         2470.   
      COEKB=-TEM2*(1.-PFRAC)/(TEM1+1.E-20)                              2470.5  
C     IF(JH.EQ.1) CKS=1./COEKB                                          2471.   
C     IF(JH.EQ.2) CKN=1./COEKB                                          2471.5  
      IF(DMOD(TAU,25.).EQ.0..AND.IABS(MRCH).EQ.1)                       2472.   
     *  WRITE(6,909) JH,COEKB,TAU                                       2472.5  
      DO 2358 J=JB,JE                                                   2473.   
      FKB=1.                                                            2473.5  
      DELTY=VLAT(J)-YZERO                                               2474.   
c     NLAT=ABS(DELTY/DLAT)+1.                                           2474.5  
      IF(JH.EQ.1.AND.DELTY.GT.0.)                                       2475.   
     *  FKB=EXP(-RADIUS*DELTY/DSCALE)
      IF(JH.EQ.2.AND.DELTY.LT.0.)                                       2476.   
     *  FKB=EXP(RADIUS*DELTY/DSCALE)
      CONKB=COEKB*FKB                                                   2477.   
 2358 TEMJ(J)=DYU(J)*(CONKB*TEMJ(J)+1.000*TSZM(J))                      2477.5  
 2355 CONTINUE                                                          2478.   
 909  FORMAT(10X,'JH =',I5,8X,'KB =',E16.5,5X,'TAU =',F10.2)            2478.5  
C     TOTLJ=JM-3                                                        2479.   
      DO 2451 J=1,JM                                                    2479.5  
 2451 FD(1,J)=PT(1,J)*DXYP(J)                                           2480.   
      FD(1,1)=2.*FD(1,1)                                                2480.5  
      FD(1,JM)=2.*FD(1,JM)                                              2481.   
      DO 2452 J=2,JM                                                    2481.5  
      RFDU=2./(FD(1,J)+FD(1,J-1))                                       2482.   
      DO 2452 L=1,LZM1                                                  2482.5  
 2452 UOLD(J,L)=UT(1,J,L)*RFDU                                          2483.   
c
      if(uvtr)then
c
C**** CALCULATE U'V' CONTRIBUTION IN U-MOMENTUM EQUATION                2483.5  
C     TDISG=1.-SIGE(LZ)                                                 2484.   
C     TEM2=.36-SIGE(LZ)*SIGE(LZ)                                        2484.5  
C     TEM3=.216-SIGE(LZ)**3.                                            2485.   
C     XCONT=TDISG/(TEM2*.3-TEM3/3.)                                     2485.5  
      if(OLDUVVERT)then
C  Old
      VVSUM=0.                                                          2486.   
      DO 2384 L=1,LM                                                    2486.5  
 2384 VVSUM=VVSUM+VVN(L)+VVS(L)                                         2487.   
      FSUM=0.                                                           2487.5  
      DO 2383 L=1,LM                                                    2488.   
 2383 FSUM=FSUM+(VVN(L)+VVS(L))*DSIG(L)/VVSUM                           2488.5  
      else
C New
      VVNEWAV=0.
      ZZD=250.
      CALL VFUNCT(VVD,ZZD)
      DO L=1,LM
	ZZZ=PHITAV(L)/GRAV+250.
	if(L.eq.LM)then
	  ZZU=30000.
	else
	  ZZU=0.5*(PHITAV(L)+PHITAV(L+1))/GRAV+250.
	endif
        CALL VFUNCT(VVNEW,ZZZ)
        CALL VFUNCT(VVU,ZZU)
	VVNEWA(L)=0.25*(VVD+2.*VVNEW+VVU)
        VVNEWAV=VVNEWAV+VVNEWA(L)*DSIG(L)
	VVD=VVU
	ZZD=ZZU
      ENDDO
      if(abs(1.-VVNEWAV).gt.1.e-1)then
         print *,' VVNEWAV=',VVNEWAV
      endif
      endif	! OLDUVVERT
      DO 2385 L=1,LZ                                                    2489.   
      if(OLDUVVERT)then
       TEM=(VVN(L)+VVS(L))/(VVSUM*FSUM)                                 2489.5  
C      IF(SIG(L).LT.0.6) TEM=XCONT*SIG(L)*(0.6-SIG(L))                  2490.   
      else
	TEM=VVNEWA(L)/VVNEWAV
      endif	! OLDUVVERT
C  New Aexp 049.98
      JB=2  
      JE=JMM1
      DO 2360 J=JB,JE  
      VU(J,L)=VU(J-1,L)*P(1,J-1)*(COSL(J-1)/COSL(J))**2.+ 
     *  TEMJ(J)*TEM*.5*(P(1,J)+P(1,J-1))/COSL(J)**2.     
      if(HPRNT)then
c      if(J.eq.JPR-1.and.L.eq.4)then
       if(L.eq.4)then
       print *,' Eddypa 2360'
         print *,VU(J-1,L),TEMJ(J),P(1,J-1)
c        print *,TEMJ(J),TEM,P(1,J)
       endif
      endif
 2360 VU(J,L)=VU(J,L)/P(1,J)                            
C  New
      DO 2370 J=JB,JE                                                   2495.   
      DO 2370 I=1,IX                                                    2495.5  
      if(HPRNT)then
       if(J.eq.JPR-1.and.L.eq.4)then
       print *,' Eddypa 2370'
         print *,DTE,DXP(J),VU(J,L),P(1,J)
         print *,FLUX,UT(I,J+1,L)
       endif
      endif
      FLUX=DTE*VU(J,L)*P (1,J)*DXP(J)                                   2496.   
      UT(I,J+1,L)=UT(I,J+1,L)+FLUX                                      2496.5  
 2370 UT(I,J,L)=UT(I,J,L)-FLUX                                          2497.   
 2385 CONTINUE                                                          2497.5  
c
      if(HPRNT)then
       print *,' Eddypa 2385'
       print *,'UT'
       print *,(UT(1,JPR,L),L=1,LM)
      endif	! PRNT
      endif
c     end uvtr
C                                                                       2498.   
      if(wvtr)then
c
C     PRESCRIBE W'U' BASED ON THE SEMI-SPECTRAL MODEL                   2498.5  
C                                                                       2499.   
C**** ADJUSTMENT OF U TO MAINTAIN TOTAL KINETIC ENERGY                  2499.5  
C                                                                       2500.   
      DO 2454 J=2,JM                                                    2500.5  
      RFDU=2./(FD(1,J)+FD(1,J-1))                                       2501.   
      DO 2454 L=1,LZM1                                                  2501.5  
 2454 UNEW(J,L)=UT(1,J,L)*RFDU                                          2502.   
      DO 2456 L=1,LZM1                                                  2502.5  
      DO 2456 JH=1,2                                                    2503.   
      JB=2                                                              2503.5  
      JE=JM                                                             2504.   
      IF(JH.EQ.1) JE=JEQM1                                              2504.5  
      IF(JH.EQ.2) JB=JEQP1                                              2505.   
      SUM1=0.                                                           2505.5  
      SUM2=0.                                                           2506.   
      DO 2458 J=JB,JE                                                   2506.5  
      TEM=0.5*(FD(1,J)+FD(1,J-1))                                       2507.   
      SUM1=SUM1+TEM*UOLD(J,L)*UOLD(J,L)                                 2507.5  
 2458 SUM2=SUM2+TEM*UNEW(J,L)*UNEW(J,L)                                 2508.   
      DIFF=(SUM2-SUM1)/SUM2                                             2508.5  
      DO 2460 J=JB,JE                                                   2509.   
      TEM=0.5*(FD(1,J)+FD(1,J-1))                                       2509.5  
      TEM1=TEM*UNEW(J,L)*UNEW(J,L)*(1.-DIFF)                            2510.   
      IF(TEM1.LE.0.) GO TO 2465                                         2510.5  
      TEM2=SQRT(TEM1/TEM)                                               2511.   
      IF(UNEW(J,L).GE.0.) UNEW(J,L)=TEM2                                2511.5  
      IF(UNEW(J,L).LT.0.) UNEW(J,L)=-TEM2                               2512.   
      GO TO 2460                                                        2512.5  
 2465 UNEW(J,L)=0.                                                      2513.   
 2460 CONTINUE                                                          2513.5  
 2456 CONTINUE                                                          2514.   
      DO 2462 J=2,JM                                                    2514.5  
      FDU=0.5*(FD(1,J)+FD(1,J-1))                                       2515.   
      DO 2462 L=1,LZM1                                                  2515.5  
 2462 UT(1,J,L)=UNEW(J,L)*FDU                                           2516.   
      if(HPRNT)then
       print *,' Eddypa 2462'
       print *,'UT'
       print *,(UT(1,JPR,L),L=1,LM)
      endif	! PRNT
      endif
c     end wvtr
C                                                                       2516.5  
      if(uvtr)then
C**** CALCULATE OTHER U'V' AND V'V' CONTRIBUTIONS IN MOMENTUM EQUATIONS 2517.   
C                                                                       2517.5  
      DO 2490 L=1,LZM1                                                  2518.   
      DO 2490 J=2,JMM1                                                  2518.5  
      FLUX=DTE*P(1,J)*DXP(J)*.5*(VV(J,L)+VV(J+1,L))                     2519.   
      VT(1,J+1,L)=VT(1,J+1,L)+FLUX                                      2519.5  
 2490 VT(1,J,L)=VT(1,J,L)-FLUX                                          2520.   
      DO 2492 L=1,LZM1                                                  2520.5  
      DO 2494 J=2,JM                                                    2521.   
      TEM=DTE*(DXP(J)-DXP(J-1))*(P(1,J)+P(1,J-1))*.5                    2521.5  
      FLUXU=.5*TEM*(VU(J,L)+VU(J-1,L))                                  2522.   
      FLUXV=TEM*VV(J,L)                                                 2522.5  
      UT(1,J,L)=UT(1,J,L)-FLUXU                                         2523.   
 2494 VT(1,J,L)=VT(1,J,L)+FLUXV                                         2523.5  
 2492 CONTINUE                                                          2524.   
      if(HPRNT)then
       print *,' Eddypa 2492'
       print *,'UT'
       print *,(UT(1,JPR,L),L=1,LM)
      endif	! PRNT
      endif
c     end uvtr
C                                                                       2524.5  
      if(uttr)then
C     PARAMETERIZATION OF MERIDIONAL EDDY HEAT TRANSPORTS               2525.   
C                                                                       2525.5  
      COEKC=1.                                                          2526.   
      DO 2245 L=1,LZ                                                    2526.5  
      DO 2240 J=2,JM                                                    2527.   
      COEKD=CKS                                                         2527.5  
      IF(J.GT.JHALF) COEKD=CKN                                          2528.   
      TEM=.5*(PHI(1,J,L)+PHI(1,J-1,L)-PHIS(1,J)-PHIS(1,J-1))            2528.5  
     *  /GRAV/(DWV(J)+1.E-20)                                           2529.   
      EXPFAC=EXP(-TEM/COEKD)                                            2529.5  
      COEJL(J,L)=COEJL(J,L)*EXPFAC                                      2530.   
      IF(L.EQ.1) COEJL(J,L)=COEJL(J,L)*.3935                            2530.5  
      IF(L.EQ.2) COEJL(J,L)=COEJL(J,L)*.8347                            2531.   
      VTH(J,L)=-.5*COEJL(J,L)*DTDYJ(J)                                  2531.5  
#if ( defined CPL_CHEM || defined TRACERS )
!
! --- Get the Kt
!
!      xxx     =(T(1,j,l)-T(1,j-1,l))/dyv(j)
      xxx     =dthdy(j,l)
      fkt(j,l)=-vth(j,l)/(xxx+1.e-20)
c      fkt(j,l)=dmax1(0.0,fkt(j,l))
c
c 012896: increase midlatitude eddy
c               here 0.0 can be changed to minimum
c               diffussion constant too
c
c       fkt(j,l)=max(2.e6,fkt(j,l)*1.5)
c       fkt(j,l)=max(0.0,fkt(j,l)*1.5)
c
c 111596: climate model change vth to 1.5 time large already:
c
        fkt(j,l)=max(0.0,fkt(j,l))
!
#endif

      FLUXT=DTE*VTH(J,L)*( P(1,J)+ P(1,J-1))*DXU(J)*.5                  2532.   
      TT(1,J,L)=TT(1,J,L)+FLUXT                                         2532.5  
      TT(1,J-1,L)=TT(1,J-1,L)-FLUXT                                     2533.   
 2240 CONTINUE                                                          2533.5  
 2245 CONTINUE                                                          2534.   
      endif
c     end uttr
C
C                                                                       2534.5  
      if(uqtr)then
C     PARAMETERIZATION OF MERIDIONAL EDDY MOISTURE TRANSPORTS           2535.   
C                                                                       2535.5  
      COEKC=1.                                                          2536.   
      DO 2345 L=1,LZ                                                    2536.5  
      DO 2340 J=2,JM                                                    2537.   
      IF(MRCH.EQ.0 ) GO TO 2340                                         2537.5  
      VQ(J,L)=VTH(J,L)*P1000*.5*(DQSDT(J,L)+DQSDT(J-1,L))               2538.   
      IF(DSHDY(J,L)*DTHDY(J,L).LT.0.) VQ(J,L)=0.                        2538.1  
      FLUXSH=VQ(J,L)*DTE*.5*(P(1,J)+P(1,J-1))*DXV(J)                    2538.5  
      IF(FLUXSH.GT..5* QT(1,J-1,L)) FLUXSH=.5* QT(1,J-1,L)              2539.   
      IF(FLUXSH.LT.-.5* QT(1,J,L)) FLUXSH=-.5* QT(1,J,L)                2539.5  
      QT(1,J,L)= QT(1,J,L)+FLUXSH                                       2540.   
      QT(1,J-1,L)= QT(1,J-1,L)-FLUXSH                                   2540.5  
 2340 CONTINUE                                                          2541.   
 2345 CONTINUE                                                          2541.5  
      endif
c     end uqtr
C                                                                       2542.   
      if(wttr) then
C     PARAMETERIZE W'T'                                                 2542.5  
C                                                                       2543.   
      COEF=.25                                                          2543.5  
      DO 4150 J=2,JMM1                                                  2544.   
      COEKD=1.                                                          2544.5  
C     IF(J.GT.JHALF) COEKD=CKN                                          2545.   
      FCONT=ABS(FO(J))                                                  2545.5  
c     TEM2=.5*(HTJ(J)+HTJ(J+1))                                         2546.   
      TEM2=HTJP(J)
      TEM3=ABS(.5*(DTDYJ(J)+DTDYJ(J+1)))                                2546.5  
      TEM=TEM3**3.*TEM2*TEM2*GRAV*GRAV/                                 2547.   
     *  (BVF(J)*THSZM(J)*THSZM(J)*FCONT*FCONT+1.E-20)                   2547.5  
      IF ((DTDYJ(J)+DTDYJ(J+1))*FO(J).GT.0.) TEM=0.                     2547.6  
      RFAC=(GAMAJ(J)/GAMAD)*(DTDZJ(J)+GAMAJ(J))/(DTDZJ(J)+GAMAD)        2551.   
      IF(RFAC.LT.0.) RFAC=0.                                            2551.01 
      RFAC=RFAC+(TEM3*GRAV/(THSZM(J)*FCONT*BVF(J)+1.E-20))**2           2551.1  
      IF(RFAC.LT.0.01) RFAC=0.01                                        2551.5  
      ALAM=.573*SQRT(RFAC)/(1.-.427*RFAC**.302)                         2552.   
c     TEMH=DWV(J)+DWV(J+1)                                              2553.   
      TEMH=DWVPA(J)
      DO 4150 L=1,LZM2                                                  2552.5  
      TEMP=(PHI(1,J,L)+PHI(1,J,L+1)-2.*PHIS(1,J))/                      2554.   
     *  (GRAV*TEMH+1.E-20)                                              2554.5  
      WTH(J,L)=-CONSTA*TEM*EXP(-TEMP/COEKD)*(P(1,J)*SIGE(L+1)+PTOP)     2555.   
     *  /SCLH(J,L)/P(1,J)*TEMP*(1.-.25*TEMP)/COEKD                      2555.5  
      WTH(J,L)=WTH(J,L)*(1.+ALAM/RFAC)/(1.+ALAM)                        2556.   
C     IF(WTH(J,L).GT.0.) WTH(J,L)=0.                                    2556.5  
 4150 CONTINUE                                                          2557.   
      DO 4160 J=2,JMM1                                                  2557.5  
      SP= P(1,J)                                                        2558.   
      DO 4160 I=1,IS                                                    2558.5  
      FLUXDT=0.                                                         2559.   
      DO 4155 L=1,LZM2                                                  2559.5  
      FLUXT=DTE*SP*DXYP(J)*WTH(J,L)                                     2560.   
      TT(I,J,L)=TT(I,J,L)+(FLUXT-FLUXDT)/DSIG(L)                        2560.5  
      FLUXDT=FLUXT                                                      2561.   
 4155 CONTINUE                                                          2561.5  
      TT(I,J,LZM1)=TT(I,J,LZM1)-FLUXDT/DSIG(LZM1)                       2562.   
 4160 CONTINUE                                                          2562.5  
      endif
c     end wttr
C
C                                                                       2563.   
      if(wqtr) then
C     PARAMETERIZE W'Q'                                                 2563.5  
C                                                                       2564.   
      IF(MRCH.ne.0) then
      DO 3170 J=2,JMM1                                                  2564.5  
      SP=P(1,J)                                                         2565.   
      BETA=BETAP(J)
      FCONT=ABS(FO(J))                                                  2565.6  
      TEM3=ABS(.5*(DTDYJ(J)+DTDYJ(J+1)))                                2565.7  
ccc   GAMAP=.5*(GAMA(J)+GAMA(J+1))                                      2566.   
C     FUNK=.55*(1.+GAMAP)*(1.+GAMAP)-.25                                2566.5  
c     FUND=SHJ(J)/(1.48*GAMAP+.48)                                      2567.   
      FUND=SHJ(J)/(2./RCZZ*(1.+GAMAP(j))-1.)
      RFAC=(GAMAJ(J)/GAMAD)*(DTDZJ(J)+GAMAJ(J))/(DTDZJ(J)+GAMAD)        2567.5  
      IF(RFAC.LT.0.) RFAC=0.                                            2567.51 
      RFAC=RFAC+(TEM3*GRAV/(THSZM(J)*FCONT*BVF(J)+1.E-20))**2           2567.6  
c     IF(RFAC.LT.0.01) RFAC=0.01                                        2568.   
      ALAM=.573*SQRT(RFAC)/(1.-.427*RFAC**.302)                         2568.5  
      TEMPB=BETA*SHJ(J)/(FO(J)*GAMAP(j)+1.E-20)                         2569.   
      IF ((DTDYJ(J)+DTDYJ(J+1))*FO(J).GT.0.) TEMPB=0.                   2569.1  

#if ( defined CPL_CHEM || defined TRACERS )
!
! --- Get parameters:
!
      beta1(j)=tempb
      beta2(j)=rfac
!
#endif

      DO 3170 I=1,IS                                                    2569.5  
      FLUXDT=0.                                                         2570.   
      FLUXDS=0.                                                         2570.5  
      DO 3175 L=1,LZM2                                                  2571.   
      ZHEI=.5*(PHI(1,J,L)+PHI(1,J,L+1)-2.*PHIS(1,J))/GRAV               2572.   
      ZOVD=ZHEI/FUND                                                    2572.5  
      TEMDN1=.25*(DSHDY(J,L)+DSHDY(J,L+1)+DSHDY(J+1,L)+DSHDY(J+1,L+1))  2573.   
      IF(J.GT.JHALF.AND.TEMDN1.GT.-DQDYMN) TEMDN1=-DQDYMN               2573.5  
      IF(J.LE.JHALF.AND.TEMDN1.LT.DQDYMN) TEMDN1=DQDYMN                 2574.   
      TERM1=ZOVD+.25*ZOVD*ZOVD*TEMPB*DSHDZ(J,L+1)/TEMDN1                2574.5  
      DSHDZE=CP*(DTDZJ(J)+GAMAJ(J))*(1.-GAMAJ(J)/GAMAD)/HLAT            2575.   
      IF(DSHDZE.LT.0.) DSHDZE=0.                                        2575.1  
      TEMDN2=.25*(DSHDYS(J,L)+DSHDYS(J,L+1)+DSHDYS(J+1,L)+              2575.5  
     *       DSHDYS(J+1,L+1))                                           2576.   
      IF(J.GT.JHALF.AND.TEMDN2.GT.-DQDYMN) TEMDN2=-DQDYMN               2576.5  
      IF(J.LE.JHALF.AND.TEMDN2.LT.DQDYMN) TEMDN2=DQDYMN                 2577.   
      TERM2=ZOVD+.25*ZOVD*ZOVD*TEMPB*DSHDZE/RFAC/TEMDN2                 2577.5  
      TEM=(TEMPB/(1.+ALAM))*(TERM1+ALAM*TERM2/RFAC)                     2578.   
      WQ(J,L)=-.25*(VQ(J,L)+VQ(J,L+1)+VQ(J+1,L)+VQ(J+1,L+1))*TEM*       2578.5  
     *       (P(1,J)*SIGE(L+1)+PTOP)/SCLH(J,L)/P(1,J)                   2579.   

#if ( defined CPL_CHEM || defined TRACERS )
!
! --- Get parameters:
!
      beta3(j,l+1)=zovd
c 112696 added /p(1,j)
      beta4(j,l+1)=(p(1,j)*sige(l+1)+ptop)/sclh(j,l)/p(1,j)
!
#endif

      FLUXS=DTE*SP*DXYP(J)*WQ(J,L)                                      2579.5  
      IF(FLUXS.GT..5* QT(I,J,L+1)*DSIG(L+1)) FLUXS=.5* QT(I,J,L+1)*     2580.   
     *  DSIG(L+1)                                                       2580.5  
      IF(FLUXS.LT.-.5* QT(I,J,L)*DSIG(L)) FLUXS=-.5* QT(I,J,L)*DSIG(L)  2581.   
       QT(I,J,L)= QT(I,J,L)+(FLUXS-FLUXDS)/DSIG(L)                      2581.5  
      FLUXDS=FLUXS                                                      2582.   
 3175 CONTINUE                                                          2582.5  
      QT(I,J,LZM1)= QT(I,J,LZM1)-FLUXDS/DSIG(LZM1) 
 3170 CONTINUE                                                          2583.5  
      endif     ! MRCH
      endif
c     end wqtr
C                                                                       2584.   
c     print *,' TPRIM from eddypa'
c     print *,TPRIM2
c     do j=1,jm
c     do l=1,lm
c     print *,TPRIM2(J,L)*(EXPBYK(SIG(L)*P(1,J)+PTOP))**2. 
c     enddo
c     enddo
    
      RETURN                                                            2584.5  
      END                                                               2585.   
      subroutine FLTR(X,JM,LM,JIN)
      dimension X(JM,LM),Y(JM,LM)
      do l=1,LM
      do j=1,JM
        Y(j,l)=X(j,l)
      enddo
      enddo
      do l=1,LM
      do j=JIN+1,JM-1
        X(j,l)=0.25*(X(j-1,l)+2.0*X(j,l)+X(j+1,l))
      enddo
      X(JIN,l)=(2.0*X(JIN,l)+X(JIN+1,l))/3.0
      X(JM,l)=(2.0*X(JM,l)+X(JM-1,l))/3.0
      enddo
      return
      end
