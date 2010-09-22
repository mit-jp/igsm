
#include "ctrparam.h"

!	==========================================================
!
!	COMP1.F:  NUMERICALLY INTEGRATING THE NON-SOURCE TERMS OF
!		  THE 2-DIMENSIONAL CLIMATE MODEL
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
!
!	==========================================================


      SUBROUTINE COMP1 (UT,VT,TT,PT,QT,U,V,T,P,Q,DT1,NS)                2001.   
C                                                                       2001.5  
C     NUMERICALLY INTEGRATING THE NON-SOURCE TERMS OF                   2002.   
C     THE 2-DIMENSIONAL CLIMATE MODEL                                   2002.5  
C                                                                       2003.   

#if ( defined CPL_CHEM )
!
#include "chem_para"
#include "chem_com"
!
#endif

#include "BD2G04.COM"                                                   2003.5  

      COMMON/WORK1/PIT(IM0,JM0),SD(IM0,JM0,LM0-1)                       2004.   
      COMMON/WORK3/PHI(IM0,JM0,LM0),SPA(IM0,JM0,LM0)                    2004.5  
      COMMON/WORK4/FD(IM0,JM0),FLUXQ(72),DUMMYS(72),DUMMYN(72)          2005.   
      COMMON/WORK5/DUT(IO0,JM0,LM0),DVT(IO0,JM0,LM0)                    2005.5  
      DIMENSION UT(IM0,JM0,LM0),VT(IM0,JM0,LM0),TT(IM0,JM0,LM0)
     &   ,PT(IM0,JM0),                                                  2006.   
     *  QT(IM0,JM0,LM0),PU(IM0,JM0),PV(IM0,JM0),CONV(IM0,JM0,LM0)       2006.5  
      DIMENSION SINL(JM0),COSL(JM0),                                    2007.   
     *  DXU(JM0),DYU(JM0),DXYU(JM0),PHIS(IO0,JM0),UX(IO0,JM0,2*LM0)     2007.5  
      EQUIVALENCE (SINL(1),SINP(1)),(COSL(1),COSP(1)),                  2008.   
     *  (DXU(1),DXV(1)),(DYU(1),DYV(1)),(DXYU(1),DXYV(1)),              2008.5  
     *  (PHIS,FDATA)                                                    2009.   
      EQUIVALENCE (CONV,PIT)                                            2009.5  
c     REAL*4 KAPAP1                                                     2010.   
      REAL KAPAP1
      COMMON/EPARA/VTH(JM0,LM0),WTH(JM0,LM0),VU(JM0,LM0),VV(JM0,LM0)
     &   ,DQSDT(JM0,LM0)                                                2010.5  
     *  ,DWV(JM0),PHIT(JM0,LM0),TPRIM2(JM0,LM0),WU(JM0,LM0),CKS,CKN     2011.   
     *  ,WQ(JM0,LM0),VQ(JM0,LM0),MRCHT                                  2011.5  
      COMMON/INTA/COE1(01,01,01),COE2(01,01,01)                         2012.   
      COMMON/SPEC1/                                                     2012.5  
     *  XA(IM0,JM0),XB(IM0,JM0),YA(IM0,JM0),YB(IM0,JM0),ZA(IM0,JM0)
     &  ,ZB(IM0,JM0)                                                    2013.   
      COMMON/SPEC2/KM,KINC,COEK,C3LAND(IO0,JM0),C3OICE(IO0,JM0)         2013.5  
     *  ,C3LICE(IO0,JM0),WMGE(IO0,JM0),TSSFC(IO0,JM0,4)                 2014.   
      LOGICAL SKIPSI,HPRNT,EDDFLX,FIRST,HPRNT1                          2014.5  
      common/PRNT1/NCOMP
      common/conprn/HPRNT1,JPR,LPR
      data FIRST /.true./
C                                                                       2015.   
      HPRNT=TAU.gt.115.0.and.TAU.lt.125.0
      HPRNT=TAU.lt.125.0
      HPRNT=.false.
      NCOMP=NCOMP+1
      KBGN=KINC+1                                                       2015.5  
      KM2=KM*2-1                                                        2016.   
      KM3=KM2                                                           2016.5  
      INEDDY=1
c     ICHK=(1+(NDYN-1)*2)*5                                             2016.6  
      ICHK=(1+(NDYN-1)*2)*INEDDY
      ICHK1=(1+(NDYN-1)*2)
      IX=IM                                                             2017.   
      IF(KM.EQ.1) IX=1                                                  2017.5  
      IS=IX                                                             2018.   
      FIS=IS                                                            2018.5  
      CP=RGAS/KAPA                                                      2019.   
      SHA=RGAS/KAPA                                                     2019.5  
      KAPAP1=KAPA+1.                                                    2020.   
      JMM2=JM-2                                                         2020.5  
      NTRACE=0                                                          2021.   
      NTLM=NTRACE*LM                                                    2021.5  
      NTP1LM=NTLM+LM                                                    2022.   
      NCOMP3=NDYN                                                       2022.5  
      LMT2=LM+LM                                                        2023.   
      LMM2=LM-2                                                         2023.5  
c     DTE=NDYN*DT*5                                                     2024.   
      DTE=NDYN*DT*INEDDY
      MRCHT=MRCHT+IABS(MRCH)                                            2024.1  
      DTT2=DT1*2.                                                       2024.5  
      DT2=DT1/2.                                                        2025.   
      DT4=DT1/4.                                                        2025.5  
      DT8=DT1/8.                                                        2026.   
      DT12=DT1/12.                                                      2026.5  
      DT24=DT1/24.                                                      2027.   
      DTCP=DT1/CP                                                       2027.5  
      DTDLN=DT1*DLON                                                    2028.   
      DTDLN2=DT2*DLON                                                   2028.5  
      SKIPSI=.TRUE.                                                     2029.   
      FMU=10.E7                                                         2029.5  
      FXCO=DT4/CP                                                       2030.   
      FXCO1=DT2/CP                                                      2030.5  
      HLAT=LHE                                                          2031.   
      CLH=HLAT/CP                                                       2031.5  
      if(FIRST)then
       print *,' EDDY DIFFUSION is CALCULATED EVERY ',INEDDY,' HOURS'
       print *,' DTE=',DTE,DTE/3600.
       print *,' ICHK=',ICHK
       FIRST=.false.
      endif
      DO 5 L=1,LMT2                                                     2032.   
      DO 5 J=1,JM                                                       2032.5  
      DO 5 I=1,IO                                                       2033.   
 5    UX(I,J,L)=U(1,J,L)                                                2033.5  
         IF(MODD5K.LT.MRCH) CALL DIAG5F(UX)                             2034.   
C                                                                       2034.5  
C     SCALE THE PROGNOSTIC VARIABLES                                    2035.   
C                                                                       2035.5  
      DO 50 J=1,JM                                                      2036.   
      DO 50 I=1,IS                                                      2036.5  
 50   FD(I,J)=PT(I,J)*DXYP(J)                                           2037.   

#if ( defined CPL_CHEM )
!
! --- get PAI for chemical advection:
!
      if(mrch.ne.0)goto 6001
      do i=1,n2dh
        p00(i,1)=fd(i,1)
!
! 092595
!
	p4chem0(i,1) = pt(i,1)
      enddo
6001  continue
!
#endif

c     print *,'From comp1 TAU=',TAU,' MRCH=',MRCH
c     print '(12f6.1/11f6.1)',(P(1,j),j=1,jm)
c     print *,' U(jm-1),U(jm),v(jm)'
c     print '(11f6.1)',(U(1,JM-1,l),l=1,lm)
c     print '(11f6.1)',(U(1,JM,l),l=1,lm)
c     print '(11f6.1)',(V(1,JM,l),l=1,lm)
c     print *,' UT(jm-1),UT(jm),VT(jm)'
c     print '(11f6.1)',(UT(1,JM-1,l),l=1,lm)
c     print '(11f6.1)',(UT(1,JM,l),l=1,lm)
c     print '(11f6.1)',(VT(1,JM,l),l=1,lm)
      IF(MRCH.EQ.0) GO TO 58                                            2037.5  
      DO 57 L=1,NTP1LM                                                  2038.   
      DO 57 J=1,JM                                                      2038.5  
      DO 57 I=1,IS                                                      2039.   
 57   QT (I,J,L)=QT (I,J,L)*FD(I,J)                                     2039.5  
 58   CONTINUE                                                          2040.   
      DO 55 L=1,LM                                                      2040.5  
      DO 55 J=1,JM                                                      2041.   
      DO 55 I=1,IX                                                      2041.5  
 55   TT(I,J,L)=TT(I,J,L)*FD(I,J)                                       2042.   
      DO 56 I=1,IX                                                      2042.5  
      FD(I,1)=2.*FD(I,1)                                                2043.   
 56   FD(I,JM)=2.*FD(I,JM)                                              2043.5  
      DO 65 J=2,JM                                                      2044.   
      DO 65 I=1,IX                                                      2044.5  
      FDU=.5*(FD(I,J)+FD(I,J-1))                                        2045.   
      DO 65 L=1,LMT2                                                    2045.5  
 65   UT(I,J,L)=UT(I,J,L)*FDU                                           2046.   
      DO 66 L=1,LMT2                                                    2046.5  
      DO 66 J=2,JM                                                      2047.   
      DO 66 I=1,IO                                                      2047.5  
 66   DUT(I,J,L)=0.                                                     2048.   
      if(HPRNT)then
      print *,' comp1 1 TAU=',TAU,' MRCH=',MRCH,' J=',JPR,' L=',LPR
      print *,' after 66'
      print *,' T(J,L)=',T(1,JPR,LPR),' Q(J,L)=',Q(1,JPR,LPR)
      print *,' TT(J,L)=',TT(1,JPR,LPR),' QT(J,L)=',QT(1,JPR,LPR)
      print *,' V(J,L)=',V(1,JPR,LPR),' V(J,L)=',V(1,JPR,LPR)
      print *,' VT(J,L)=',VT(1,JPR,LPR),' VT(J,L)=',VT(1,JPR,LPR)
      print *,' U(J,L)=',U(1,JPR,LPR),' U(J,L)=',U(1,JPR,LPR)
      print *,' UT(J,L)=',UT(1,JPR,LPR),' UT(J,L)=',UT(1,JPR,LPR)
      endif
C                                                                       2048.5  
C     BEGINNING OF LAYER LOOP                                           2049.   
C                                                                       2049.5  
ccc   L=LM                                                              2050.   
c     print *,' before 5934 MRCH=',MRCH
c     RFDP=2./(P(1,JM-2)*DXYP(JM-2)+P(1,JM-1)*DXYP(JM-1)) 
c     print '(11f6.1)',(UT(1,JM-1,l1)*RFDP,l1=1,lm)
c     RFDP=2./(2.*P(1,JM)*DXYP(JM)+P(1,JM-1)*DXYP(JM-1)) 
c     print '(11f6.1)',(UT(1,JM,l1)*RFDP,l1=1,lm)
c     print '(11f6.1)',(VT(1,JM,l1)*RFDP,l1=1,lm)
      do 5934 l25=1,LM
         L=LM+1-l25
C                                                                       2050.5  
C     COMPUTATION OF MASS FLUX  P,T,PU  PRIMARY GRID ROW                2051.   
C                               PV,U,V  SECONDARY GRID ROW              2051.5  
C                                                                       2052.   
 2150 CONTINUE                                                          2052.5  
      DO 2168 J=2,JM                                                    2053.   
      TEMP=DXU(J)*0.5                                                   2053.5  
      DO 2168 I=1,IX                                                    2054.   
      PV(I,J)=(P(I,J)+P(I,J-1))*V(I,J,L)*TEMP                           2054.5  

#if ( defined CPL_CHEM )
!
! --- get V for chemical advection:
!
      if(mrch.eq.0)then
        pvv(i,j,l)=v(i,j,l)
      endif
!
#endif

2168  continue
      if(HPRNT)then
       print *,' PV L=',L
       print *,PV(1,JPR),P(1,JPR),P(1,JPR-1),V(1,JPR,L),TEMP
      endif
C                                                                       2055.   
C     HORIZONTAL ADVECTION OF MOISTURE                                  2055.5  
C                                                                       2056.   
      IF(MRCH.EQ.0) GO TO 1400                                          2056.5  
      DO 1360 K=L,NTP1LM,LM                                             2057.   
C**** SOUTH-NORTH ADVECTION OF MOISTURE AND TRACE COMPOUNDS             2057.5  
      DO 1330 I=1,IS                                                    2058.   
      FLUXQ(I)=DT2*PV(I,2)*(Q (I,2,K)+Q (I,1,K))                        2058.5  
      IF(FLUXQ(I).GT. QT(I,1,K)) FLUXQ(I)= QT(I,1,K)                    2059.   
      IF(FLUXQ(I).LT.-.5* QT(I,2,K)) FLUXQ(I)=-.5* QT(I,2,K)            2059.5  
 1330  QT(I,1,K)= QT(I,1,K)-FLUXQ(I)                                    2060.   
      DO 1340 J=2,JMM2                                                  2060.5  
      DO 1340 I=1,IS                                                    2061.   
      FLUX=DT2*PV(I,J+1)*( Q(I,J+1,K)+ Q(I,J,K))                        2061.5  
      IF(FLUX.GT..5* QT(I,J,K)) FLUX=.5* QT(I,J,K)                      2062.   
      IF(FLUX.LT.-.5* QT(I,J+1,K)) FLUX=-.5* QT(I,J+1,K)                2062.5  
       QT(I,J,K)= QT(I,J,K)+(FLUXQ(I)-FLUX)                             2063.   
 1340 FLUXQ(I)=FLUX                                                     2063.5  
      DO 1350 I=1,IS                                                    2064.   
      FLUX=DT2*PV(I,JM)*( Q(I,JM,K)+ Q(I,JMM1,K))                       2064.5  
      IF(FLUX.GT..5* QT(I,JMM1,K)) FLUX=.5* QT(I,JMM1,K)                2065.   
      IF(FLUX.LT.- QT(I,JM,K)) FLUX=- QT(I,JM,K)                        2065.5  
       QT(I,JMM1,K)= QT(I,JMM1,K)+(FLUXQ(I)-FLUX)                       2066.   
 1350  QT(I,JM,K)= QT(I,JM,K)+FLUX                                      2066.5  
 1360 CONTINUE                                                          2067.   
 1400 CONTINUE                                                          2067.5  
C                                                                       2068.   
C     HORIZONTAL ADVECTION OF HEAT                                      2068.5  
C                                                                       2069.   
      DO 2215 J=2,JM                                                    2069.5  
      DO 2215 I=1,IX                                                    2070.   
      FLUX=DT2*PV(I,J)                                                  2070.5  
      FLUXT=FLUX*(T(I,J,L)+T(I,J-1,L))                                  2071.   
      TT(I,J,L)=TT(I,J,L)+FLUXT                                         2071.5  
 2215 TT(I,J-1,L)=TT(I,J-1,L)-FLUXT                                     2072.   
C                                                                       2072.5  
C     MERIDIONAL ADVECTION OF U AND V-MOMENTUM                          2073.   
C                                                                       2073.5  
      DO 2336 J=2,JMM1                                                  2074.   
      DO 2336 I=1,IX                                                    2074.5  
      FLUX=DT4*(PV(I,J)+PV(I,J+1))                                      2075.   
      FLUXU=FLUX*(U(I,J,L)+U(I,J+1,L))                                  2075.5  
C*** CONTRIBUTION BY SYMMETRIC INSTABILITY                              2076.   
      IF (MRCH.EQ.0.OR.SKIPSI) GO TO 2335                               2076.5  
      DUDY=(U(I,J+1,L)*COSV(J+1)-U(I,J,L)*COSV(J))/                     2077.   
     *     (DYP(J)*COSP(J))                                             2077.5  
      FTEMP=F(J)/DXYP(J)                                                2078.   
      CRI=FTEMP*(FTEMP-DUDY)                                            2078.5  
      IF (CRI.GE.0.) GO TO 2335                                         2079.   
      FLUXDI=-DT1*FMU*P(I,J)*DUDY*DXP(J)                                2079.5  
      FLUXU=FLUXU+FLUXDI                                                2080.   
 2335 CONTINUE                                                          2080.5  
      UT(I,J+1,L)=UT(I,J+1,L)+FLUXU                                     2081.   
      UT(I,J,L)=UT(I,J,L)-FLUXU                                         2081.5  
      FLUXV=FLUX*(V(I,J,L)+V(I,J+1,L))                                  2082.   
      VT(I,J+1,L)=VT(I,J+1,L)+FLUXV                                     2082.5  
      VT(I,J,L)=VT(I,J,L)-FLUXV                                         2083.   
      if(HPRNT.and.J.eq.JPR.and.L.eq.LPR)then
       print *,' before 2336 c',' J=',J,' L=',L
        print *,'U(I,J+1,L)=',U(I,J+1,L),' U(I,J,L)=',U(I,J,L)
        print *,' DUDY=',DUDY,' F(J)=',F(J),' CRI=',CRI
        print *,'V(I,J,L)=',V(I,J,L),' V(I,J+1,L)=',V(I,J+1,L)
        print *,' PV(I,J)=',PV(I,J),' PV(I,J+1)=',PV(I,J+1)
        print *,' FLUX=',FLUX,' FLUXV=',FLUXV
        print *,' FLUXU=',FLUXU,' FLUXDI=',FLUXDI
      endif
 2336 continue
      if(l25.eq.-LM)then
      print *,' after 2336'
      RFDP=2./(P(1,JM-2)*DXYP(JM-2)+P(1,JM-1)*DXYP(JM-1)) 
      print '(11f6.1)',(UT(1,JM-1,l1)*RFDP,l1=1,lm)
      RFDP=2./(2.*P(1,JM)*DXYP(JM)+P(1,JM-1)*DXYP(JM-1)) 
      print '(11f6.1)',(UT(1,JM,l1)*RFDP,l1=1,lm)
      print '(11f6.1)',(VT(1,JM,l1)*RFDP,l1=1,lm)
      endif
      if(HPRNT)then
      print *,' comp1 1 TAU=',TAU,' MRCH=',MRCH
      print *,' after 2336'
      print *,' T(J,L)=',T(1,JPR,LPR),' Q(J,L)=',Q(1,JPR,LPR)
      print *,' TT(J,L)=',TT(1,JPR,LPR),' QT(J,L)=',QT(1,JPR,LPR)
      print *,' V(J,L)=',V(1,JPR,LPR),' V(J,L)=',V(1,JPR,LPR)
      print *,' VT(J,L)=',VT(1,JPR,LPR),' VT(J,L)=',VT(1,JPR,LPR)
      print *,' U(J,L)=',U(1,JPR,LPR),' U(J,L)=',U(1,JPR,LPR)
      print *,' UT(J,L)=',UT(1,JPR,LPR),' UT(J,L)=',UT(1,JPR,LPR)
      endif
C     DO 2337 J=2,JMM1                                                  2083.1  
C2337 IF(L.EQ.8) WRITE(6,2339) J,MRCH,V(1,J,8),VT(1,J,8),U(1,J,8),      2083.2  
C    *                         UT(1,J,8)                                2083.21 
 2339 FORMAT(1X,'AFTER M ADV---V VT U UT',2I3,4E11.3)                   2083.3  
C                                                                       2083.5  
C     COMPUTE MASS CONVERGENCE                                          2084.   
C                                                                       2084.5  
      DSIGL=DSIG(L)                                                     2085.   
      PVS=PV(1,2)                                                       2085.5  
      PVN=PV(1,JM)                                                      2086.   
C 
      DO 2400 J=2,JMM1                                                  2086.5  
      DO 2400 I=1,IX                                                    2087.   
 2400 CONV(I,J,L)=(PV(I,J)-PV(I,J+1))*DSIGL                             2087.5  
      CONV(1,1,L)=-PVS*DSIGL                                            2088.   
      CONV(1,JM,L)=PVN*DSIGL                                            2088.5  
      if(HPRNT)then
       print *,' PV 3 L=',L,' JPR=',JPR
       print *,PV(1,JPR),PV(1,JPR+1),DSIGL,CONV(1,JPR,L)
       print *,DSIG(1),CONV(1,1,1)
      endif
C
c2409 L=L-1                                                             2089.   
c     IF(L.GE.1) GO TO 2150                                             2089.5  
 5934 continue
      if(HPRNT)then
      print *,' comp1 2 TAU=',TAU,' MRCH=',MRCH
      print *,' T(J,L)=',T(1,JPR,LPR),' Q(J,L)=',Q(1,JPR,LPR)
      print *,' TT(J,L)=',TT(1,JPR,LPR),' QT(J,L)=',QT(1,JPR,LPR)
      print *,' PV(1,J) PV(1,J+1) PIT(1,J)'
      print *,PV(1,JPR),PV(1,JPR+1),PIT(1,JPR)
       print *,' CONV(1,JPR,L)'
       print *,(CONV(1,JPR,L),L=1,LM)
      endif
C                                                                       2090.   
C     END OF LAYER LOOP                                                 2090.5  
C                                                                       2091.   
C     COMPUTE PRESSURE TENDENCY AND SIGMA DOT                           2091.5  
C                                                                       2092.   
C     PIT(I,J)=CONV(I,J,1)                                              2092.5  
      DO 2420 LX=1,LMM1                                                 2093.   
      L=LM+1-LX                                                         2093.5  
      DO 2420 J=1,JM                                                    2094.   
      DO 2420 I=1,IS                                                    2094.5  
 2420 PIT(I,J)=PIT(I,J)+CONV(I,J,L)                                     2095.   
      DO 2430 J=1,JM                                                    2095.5  
      DO 2430 I=1,IS                                                    2096.   
 2430 SD(I,J,LMM1)=CONV(I,J,LM)-DSIG(LM)*PIT(I,J)                       2096.5  
      DO 2440 LX=2,LMM1                                                 2097.   
      L=LM-LX                                                           2097.5  
      DO 2440 J=1,JM                                                    2098.   
      DO 2440 I=1,IS                                                    2098.5  
 2440 SD(I,J,L)=SD(I,J,L+1)+CONV(I,J,L+1)-DSIG(L+1)*PIT(I,J)            2099.   

#if ( defined CPL_CHEM )
!
! --- get Omiga for chemical advection:
!
      if(mrch.ne.0)goto 6009
      do 6008 i=1,nlon
      do 6008 k=1,nlev
      do 6008 j=1,nlat
        pww(i,j,k)=sd(i,j,k)/p00(i,j)
6008  continue

6009  continue
!
#endif

      if(HPRNT)then
       j=jm
       print *,' comp 21 J=',J,' L=',L
       print *,' PIT(1,J)=',PIT(1,J)
       print *,' SD ',(SD(1,J,L),l=1,lm-1)
       print *,' CONV(1,J,L)'
       print *,(CONV(1,J,L),L=1,LM)
      endif
C                                                                       2099.5  
C     COMPUTE THE NEW SURFACE PRESSURE                                  2100.   
C                                                                       2100.5  
      DO 2450 J=1,JM                                                    2101.   
      DO 2450 I=1,IS                                                    2101.5  
      PTNEW=PT(I,J)+DT1*PIT(I,J)/DXYP(J)                                2102.   
      IF(PTNEW.LE.1150..AND.PTNEW.GT.400.) GO TO 2450                   2102.5  
 4563 CONTINUE
      print *,' TAU=',TAU
      WRITE(6,901) I,J,MRCH,P(I,J),PTNEW,(FDATA(1,J,K),K=1,22),         2103.   
     *  (V(1,J,L),L=1,LM),(V(1,J+1,L),L=1,LM),                          2103.1  
     *  (T(1,J,L),L=1,LM),(Q(1,J,L),L=1,LM)                             2103.2  
 901  FORMAT(' PRESSURE DIAGNOSTIC     I,J,MRCH,P,PT=',3I5,2F10.1/,     2103.5  
     * 2(1X,11F10.2/),4(1X,11E11.3/)) 
c    * 2(1X,11F10.2/),4(10X,9E11.3/))                                   2103.6  
      stop
 2450 PT(I,J)=PTNEW                                                     2104.   
C                                                                       2104.5  
C     VERTICAL ADVECTION OF MOMENTUM                                    2105.   
C                                                                       2105.5  
c     print *,' before 2480'
c     RFDP=2./(P(1,JM-2)*DXYP(JM-2)+P(1,JM-1)*DXYP(JM-1)) 
c     print '(11f6.1)',(UT(1,JM-1,l1)*RFDP,l1=1,lm)
c     RFDP=2./(2.*P(1,JM)*DXYP(JM)+P(1,JM-1)*DXYP(JM-1)) 
c     print '(11f6.1)',(UT(1,JM,l1)*RFDP,l1=1,lm)
c     print '(11f6.1)',(VT(1,JM,l1)*RFDP,l1=1,lm)
      DO 2480 L=1,LMM1                                                  2106.   
      LP1=L+1                                                           2106.5  
      WTS=2.                                                            2107.   
      WTN=1.                                                            2107.5  
      DO 2480 J=2,JM                                                    2108.   
      IF(J.EQ.JM) WTN=2.                                                2108.5  
      DO 2470 I=1,IX                                                    2109.   
      SDU=DT4*(SD(I,J,L)*WTN+SD(I,J-1,L)*WTS)                           2109.5  
      SDUDN=SDU/DSIG(L)                                                 2110.   
      SDUUP=SDU/DSIG(LP1)                                               2110.5  
      TEMU=U(I,J,L)+U(I,J,LP1)                                          2111.   
      TEMV=V(I,J,L)+V(I,J,LP1)                                          2111.5  
      UT(I,J,L)=UT(I,J,L)+SDUDN*TEMU                                    2112.   
      UT(I,J,LP1)=UT(I,J,LP1)-SDUUP*TEMU                                2112.5  
      VT(I,J,L)=VT(I,J,L)+SDUDN*TEMV                                    2113.   
 2470 VT(I,J,LP1)=VT(I,J,LP1)-SDUUP*TEMV                                2113.5  
 2480 WTS=1.                                                            2114.   
c     print *,' after 2480'
c     RFDP=2./(P(1,JM-2)*DXYP(JM-2)+P(1,JM-1)*DXYP(JM-1)) 
c     print '(11f6.1)',(UT(1,JM-1,l1)*RFDP,l1=1,lm)
c     RFDP=2./(2.*P(1,JM)*DXYP(JM)+P(1,JM-1)*DXYP(JM-1)) 
c     print '(11f6.1)',(UT(1,JM,l1)*RFDP,l1=1,lm)
c     print '(11f6.1)',(VT(1,JM,l1)*RFDP,l1=1,lm)
      if(HPRNT)then
      print *,' comp1 1 TAU=',TAU,' MRCH=',MRCH
      print *,' after 2480'
      print *,' T(J,L)=',T(1,JPR,LPR),' Q(J,L)=',Q(1,JPR,LPR)
      print *,' TT(J,L)=',TT(1,JPR,LPR),' QT(J,L)=',QT(1,JPR,LPR)
      print *,' V(J,L)=',V(1,JPR,LPR),' V(J,L)=',V(1,JPR,LPR)
      print *,' VT(J,L)=',VT(1,JPR,LPR),' VT(J,L)=',VT(1,JPR,LPR)
      print *,' U(J,L)=',U(1,JPR,LPR),' U(J,L)=',U(1,JPR,LPR)
      print *,' UT(J,L)=',UT(1,JPR,LPR),' UT(J,L)=',UT(1,JPR,LPR)
      endif
C     DO 2481 J=2,JMM1                                                  2114.1  
C2481 WRITE(6,2482) J,MRCH,V(1,J,8),VT(1,J,8),U(1,J,8),UT(1,J,8),       2114.2  
C    *              SD(1,J,8)                                           2114.21 
 2482 FORMAT(1X,'AFTER V ADV---V VT U UT SD',2I3,5E11.3)                2114.3  
         IF(MODD5K.LT.MRCH) CALL DIAG5A(4,MRCH)                         2114.5  
C        IF(MRCH.GT.0) CALL DIAG9D (1,DT1,U,V)                          2115.   
C                                                                       2115.5  
C     VERTICAL ADVECTION OF MOISTURE                                    2116.   
C                                                                       2116.5  
      IF(MRCH.EQ.0) GO TO 1800                                          2117.   
      DO 1730 L=1,LMM1                                                  2117.5  
      LP1=L+1                                                           2118.   
      DO 1710 J=1,JM                                                    2118.5  
      DO 1710 I=1,IS                                                    2119.   
      FLUX=DT1*SD(I,J,L)*2.* Q(I,J,L)* Q(I,J,LP1)/( Q(I,J,L)+           2119.5  
     *   Q(I,J,LP1)+1.E-20)                                             2120.   
      IF(FLUX.GT..5* QT(I,J,LP1)*DSIG(LP1)) FLUX=.5* QT(I,J,LP1)*       2120.5  
     *  DSIG(LP1)                                                       2121.   
      IF(FLUX.LT.-.5* QT(I,J,L)*DSIG(L)) FLUX=-.5* QT(I,J,L)*DSIG(L)    2121.5  
       QT(I,J,L)= QT(I,J,L)+FLUX/DSIG(L)                                2122.   
 1710  QT(I,J,LP1)= QT(I,J,LP1)-FLUX/DSIG(LP1)                          2122.5  
 1730 CONTINUE                                                          2123.   
 1800 CONTINUE                                                          2123.5  
C                                                                       2124.   
C     CORIOLIS FORCE                                                    2124.5  
c     print *,' before 3130'
c     RFDP=2./(P(1,JM-2)*DXYP(JM-2)+P(1,JM-1)*DXYP(JM-1)) 
c     print '(11f6.1)',(UT(1,JM-1,l1)*RFDP,l1=1,lm)
c     RFDP=2./(2.*P(1,JM)*DXYP(JM)+P(1,JM-1)*DXYP(JM-1)) 
c     print '(11f6.1)',(UT(1,JM,l1)*RFDP,l1=1,lm)
c     print '(11f6.1)',(VT(1,JM,l1)*RFDP,l1=1,lm)
C                                                                       2125.   
      DO 3100 I=1,IX                                                    2125.5  
      FD(I,1)=0.                                                        2126.   
 3100 FD(I,JM)=0.                                                       2126.5  
      DO 3130 L=1,LM                                                    2127.   
      DO 3110 J=2,JMM1                                                  2127.5  
      DO 3110 I=1,IX                                                    2128.   
 3110 FD(I,J)=F(J)        +0.5*(U(I,J,L)+U(I,J+1,L))*(DXU(J)-DXU(J+1))  2128.5  
      DO 3115 J=2,JM                                                    2129.   
      DO 3115 I=1,IX                                                    2129.5  
      PU(I,J)=DT4*(P(I,J)+P(I,J-1))*(FD(I,J)+FD(I,J-1))                 2130.   
 3115 CONTINUE                                                          2132.   
      DO 3120 J=2,JM                                                    2132.5  
      DO 3120 I=1,IX                                                    2133.   
      UT(I,J,L)=UT(I,J,L)+     V(I,J,L)*PU(I,J)                         2133.5  
 3120 VT(I,J,L)=VT(I,J,L)-     U(I,J,L)*PU(I,J)                         2134.   
 3130 CONTINUE                                                          2134.5  
c     print *,' after 3130'
c     RFDP=2./(P(1,JM-2)*DXYP(JM-2)+P(1,JM-1)*DXYP(JM-1)) 
c     print '(11f6.1)',(UT(1,JM-1,l1)*RFDP,l1=1,lm)
c     RFDP=2./(2.*P(1,JM)*DXYP(JM)+P(1,JM-1)*DXYP(JM-1)) 
c     print '(11f6.1)',(UT(1,JM,l1)*RFDP,l1=1,lm)
c     print '(11f6.1)',(VT(1,JM,l1)*RFDP,l1=1,lm)
      if(HPRNT)then
      print *,' comp1 1 TAU=',TAU,' MRCH=',MRCH
      print *,' after 3130'
      print *,' T(J,L)=',T(1,JPR,LPR),' Q(J,L)=',Q(1,JPR,LPR)
      print *,' TT(J,L)=',TT(1,JPR,LPR),' QT(J,L)=',QT(1,JPR,LPR)
      print *,' V(J,L)=',V(1,JPR,LPR),' V(J,L)=',V(1,JPR,LPR)
      print *,' VT(J,L)=',VT(1,JPR,LPR),' VT(J,L)=',VT(1,JPR,LPR)
      print *,' U(J,L)=',U(1,JPR,LPR),' U(J,L)=',U(1,JPR,LPR)
      print *,' UT(J,L)=',UT(1,JPR,LPR),' UT(J,L)=',UT(1,JPR,LPR)
      endif
C     DO 3131 J=2,JMM1                                                  2134.6  
C3131 WRITE (6,3132) J,MRCH,V(1,J,8),VT(1,J,8),U(1,J,8),UT(1,J,8)       2134.7  
 3132 FORMAT(1X,'AFTER COR---V VT U UT',2I3,4E11.3)                     2134.8  
         IF(MODD5K.LT.MRCH) CALL DIAG5A(5,MRCH)                         2135.   
C        IF(MRCH.GT.0) CALL DIAG9D (2,DT1,U,V)                          2135.5  
C                                                                       2136.   
C     COMPUTE SPA, PHI AND VERTICAL ADVECTION OF POTENTIAL TEMPERATURE  2136.5  
C                                                                       2137.   
      DO 3070 J=1,JM                                                    2137.5  
      DO 3070 I=1,IX                                                    2138.   
      SUM1=0.                                                           2138.5  
      SUM2=0.                                                           2139.   
      SP=P(I,J)                                                         2139.5  
      PDN=SIG(1)*SP+PTOP                                                2140.   
      PKDN=EXPBYK(PDN)                                                  2140.5  
      FLUXDN=0.                                                         2141.   
      DO 3040 L=1,LMM1                                                  2141.5  
      LP1=L+1                                                           2142.   
      SPA(I,J,L)=RGAS*SP*SIG(L)*T(I,J,L)*PKDN/PDN                       2142.5  
      SUM1=SUM1+SPA(I,J,L)*DSIG(L)                                      2143.   
      PUP=SIG(LP1)*SP+PTOP                                              2143.5  
      PKUP=EXPBYK(PUP)                                                  2144.   
      THETA=0.5*(T(I,J,LP1)+T(I,J,L))                                   2144.5  
      FLUX=DT1*SD(I,J,L)*THETA                                          2145.   
      TT(I,J,L)=TT(I,J,L)+(FLUX-FLUXDN)/DSIG(L)                         2145.5  
      FLUXDN=FLUX                                                       2146.   
      PHI(I,J,LP1)=CP*THETA*(PKDN-PKUP)                                 2146.5  
      PDN=PUP                                                           2147.   
      PKDN=PKUP                                                         2147.5  
 3040 SUM2=SUM2+SIGE(LP1)*PHI(I,J,LP1)                                  2148.   
      SPA(I,J,LM)=SIG(LM)*SP*RGAS*T(I,J,LM)*PKDN/PDN                    2148.5  
      SUM1=SUM1+SPA(I,J,LM)*DSIG(LM)                                    2149.   
      TT(I,J,LM)=TT(I,J,LM)-FLUXDN/DSIG(LM)                             2149.5  
 3050 PHI(I,J,1)=PHIS(I,J)+SUM1-SUM2                                    2150.   
      DO 3070 L=2,LM                                                    2150.5  
 3070 PHI(I,J,L)=PHI(I,J,L)+PHI(I,J,L-1)                                2151.   
      if(HPRNT)then
      print *,' comp1 3 TAU=',TAU,' MRCH=',MRCH
      print *,' T(J,L)=',T(1,JPR,LPR),' Q(J,L)=',Q(1,JPR,LPR)
      print *,' TT(J,L)=',TT(1,JPR,LPR),' QT(J,L)=',QT(1,JPR,LPR)
      endif
C
C                                                                       2151.5  
      DO 3340 L=1,LM                                                    2152.   
C                                                                       2152.5  
C     PRESSURE GRADIENT FORCE                                           2153.   
C                                                                       2153.5  
C     NORTH-SOUTH DERIVATIVE AFFECTS THE V-MOMENTUM                     2154.   
C                                                                       2154.5  
      DO 3236 J=2,JM                                                    2155.   
      DO 3236 I=1,IX                                                    2155.5  
      FLUX=DT2*((P(I,J)+P(I,J-1))*(PHI(I,J,L)-PHI(I,J-1,L))+            2156.   
     *  (SPA(I,J,L)+SPA(I,J-1,L))*(P(I,J)-P(I,J-1)))*DXU(J)             2156.5  
 3236 VT(I,J,L)=VT(I,J,L)-FLUX                                          2157.   
 3340 CONTINUE                                                          2157.5  
      if(HPRNT)then
      print *,' after 3340'
      print *,' comp1 1 TAU=',TAU,' MRCH=',MRCH
      print *,' T(J,L)=',T(1,JPR,LPR),' Q(J,L)=',Q(1,JPR,LPR)
      print *,' TT(J,L)=',TT(1,JPR,LPR),' QT(J,L)=',QT(1,JPR,LPR)
      print *,' V(J,L)=',V(1,JPR,LPR),' V(J,L)=',V(1,JPR,LPR)
      print *,' VT(J,L)=',VT(1,JPR,LPR),' VT(J,L)=',VT(1,JPR,LPR)
      print *,' U(J,L)=',U(1,JPR,LPR),' U(J,L)=',U(1,JPR,LPR)
      print *,' UT(J,L)=',UT(1,JPR,LPR),' UT(J,L)=',UT(1,JPR,LPR)
      endif
C     DO 3341 J=2,JMM1                                                  2157.6  
C3341 WRITE(6,3342) J,MRCH,P(1,J),VT(1,J,8),UT(1,J,8),PHI(1,J,8),       2157.7  
C    *              SPA(1,J,8)                                          2157.71 
 3342 FORMAT(1X,'AFTER P GRAD---P VT UT PHI SPA',2I3,5E11.3)            2157.8  
         IF(MODD5K.LT.MRCH) CALL DIAG5A(6,MRCH)                         2158.   
C        IF(MRCH.GT.0) CALL DIAG9D (3,DT1,U,V)                          2158.5  
C                                                                       2159.   
C     END OF THE PRESSURE GRADIENT FORCE LAYER LOOP                     2159.5  
C                                                                       2160.   
C     CALL EDDY PARA. ROUTINE IF REQUIRED                               2160.5  
C                                                                       2161.   
C     IF (MRCH.NE.0) CALL EDDYPA (UT,VT,TT,PT,QT,U,V,T,P,Q,DT1)         2161.5  
      MODEDY=MOD(MRCHT,ICHK)                                            2161.6  
      MODTVR=MOD(MRCHT,ICHK1)

#if ( defined CPL_CHEM )
!
! --- Get time parameter for eddy calculation:
!
      meddy1=0

      IF (MODEDY.EQ.0.AND.MRCH.NE.0) THEN

        meddy1=1
!
#else
      IF (MODEDY.EQ.0.AND.MRCH.NE.0) THEN
#endif

        EDDFLX=.false.
        IF (MODEDY.eq.0)EDDFLX=.true.
        EDDFLX=.true.
	if(HPRNT)then
       print *,' before eddypa'
      print *,' comp1 1 TAU=',TAU,' MRCH=',MRCH
	print *,' JPR=',jpr,' LPR=',lpr
      print *,' T=',T(1,JPR,LPR),' Q=',Q(1,JPR,LPR)
      print *,' TT=',TT(1,JPR,LPR),' QT=',QT(1,JPR,LPR)
      print *,' V(J,L)=',V(1,JPR,LPR)/FD(1,J)
      print *,' VT(J,L)=',VT(1,JPR,LPR),' VT(J,L)=',VT(1,JPR,LPR)
      print *,' U(J,L)=',U(1,JPR,LPR),' U(J,L)=',U(1,JPR,LPR)
      print *,' UT(J,L)=',UT(1,JPR,LPR),' UT(J,L)=',UT(1,JPR,LPR)
      print *,' U(13,1),U(13,2)'
      print *,' U(J,L)=',U(1,13,1),' U(J,L)=',U(1,13,2)
      endif
!      print *,' call eddypa TAU=',TAU
!      print *,' comp1 MRCHT=',MRCHT,' MRCH=',MRCH
c     RFDU=2./(2.*P(1,JM)*DXYP(JM)+P(1,JM-1)*DXYP(JM-1)) 
c     print *,(VT(1,J,1),J=1,Jm)
        CALL EDDYPA (UT,VT,TT,PT,QT,U,V,T,P,Q,DTE,EDDFLX)
c      print *,' after eddypa'
c     print *,(VT(1,J,1),J=1,Jm)
      if(HPRNT)then
       print *,' after eddypa'
      print *,' comp1 1 TAU=',TAU,' MRCH=',MRCH
      print *,' T=',T(1,JPR,LPR),' Q=',Q(1,JPR,LPR)
      print *,' TT=',TT(1,JPR,LPR),' QT=',QT(1,JPR,LPR)
      print *,' V(J,L)=',V(1,JPR,LPR)/FD(1,J)
      print *,' VT(J,L)=',VT(1,JPR,LPR),' VT(J,L)=',VT(1,JPR,LPR)
      print *,' U(J,L)=',U(1,JPR,LPR),' U(J,L)=',U(1,JPR,LPR)
      print *,' UT(J,L)=',UT(1,JPR,LPR),' UT(J,L)=',UT(1,JPR,LPR)
      endif
      END IF
C
C                                                                       2162.   
C     UNDO SCALING FOR MOISTURE AND TRACE COMPOUNDS                     2162.5  
C                                                                       2163.   
      DO 1910 J=1,JM                                                    2163.5  
      DO 1910 I=1,IM                                                    2164.   
 1910 FD(I,J)=PT(I,J)*DXYP(J)                                           2164.5  

#if ( defined CPL_CHEM )
!
! --- get PAI for chemical advection (mass remapping):
!
      if(mrch.ne.2)goto 6003
      do i=1,n2dh
        p11(i,1)=fd(i,1)
!
! 092595
!
	p4chem1(i,1) = pt(i,1)
      enddo
6003  continue
!
#endif

      IF(MRCH.EQ.0) GO TO 2000                                          2165.   
      DO 1940 L=1,NTP1LM                                                2165.5  
      DO 1928 J=1,JM                                                    2166.   
      DO 1928 I=1,IM                                                    2166.5  
 1928  QT(I,J,L)= QT(I,J,L)/FD(I,J)                                     2167.   
 1940 CONTINUE                                                          2167.5  
 2000 CONTINUE                                                          2168.   
C                                                                       2168.5  
C     UNDO SCALING FOR TEMPERATURE AND WINDS                            2169.   
C                                                                       2169.5  
      DO 3510 L=1,LM                                                    2170.   
      DO 3512 J=1,JM                                                    2170.5  
      DO 3512 I=1,IM                                                    2171.   
 3512 TT(I,J,L)=TT(I,J,L)/FD(I,J)                                       2171.5  
      DO 3513 J=1,JM                                                    2172.   
      DO 3513 I=1,IX                                                    2172.5  
c     IF(TT(I,J,L).LT.20..OR.TT(I,J,L).GT.120.)                         2173.   
      IF(TT(I,J,L).LT.20..OR.TT(I,J,L).GT.130.) then
        WRITE(6,902) I,J,L,MRCH,T(I,J,L),TT(I,J,L),P(I,J),PT(I,J)       2173.5  
        stop
      endif
 902  FORMAT(' POTENTIAL TEMPERATURE DIAGNOSTIC     I,J,L,MRCH,T,TT,P,PT2174.   
     * =',4I4,4F8.1)                                                    2174.5  
 3513 CONTINUE                                                          2175.   
 3510 CONTINUE                                                          2175.5  
      DO 3520 I=1,IM                                                    2176.   
      FD(I,1)=2.*FD(I,1)                                                2176.5  
 3520 FD(I,JM)=2.*FD(I,JM)                                              2177.   
      DO 3542 J=2,JM                                                    2177.5  
      DO 3542 I=1,IM                                                    2178.   
      RFDU=2./(FD(I,J)+FD(I,J-1))                                       2178.5  
      DO 3542 L=1,LMT2                                                  2179.   
 3542 UT(I,J,L)=UT(I,J,L)*RFDU                                          2179.5  
      if(HPRNT)then
      print *,' after 3542'
      print *,' comp1 1 TAU=',TAU,' MRCH=',MRCH
      print *,' T(J,L)=',T(1,JPR,LPR),' Q(J,L)=',Q(1,JPR,LPR)
      print *,' TT(J,L)=',TT(1,JPR,LPR),' QT(J,L)=',QT(1,JPR,LPR)
      print *,' V(J,L)=',V(1,JPR,LPR),' V(J,L)=',V(1,JPR,LPR)
      print *,' VT(J,L)=',VT(1,JPR,LPR),' VT(J,L)=',VT(1,JPR,LPR)
      print *,' U(J,L)=',U(1,JPR,LPR),' U(J,L)=',U(1,JPR,LPR)
      print *,' UT(J,L)=',UT(1,JPR,LPR),' UT(J,L)=',UT(1,JPR,LPR)
      endif
c     print '(12f6.1/11f6.1)',(PT(1,j),j=1,jm)
      RETURN                                                            2180.   
      END                                                               2180.5  
