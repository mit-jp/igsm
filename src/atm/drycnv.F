
#include "ctrparam.h"

!	==========================================================
!
!	DRYCNV.F:  THIS SUBROUTINE MIXES AIR CAUSED BY DRY 
!		     CONVECTION.  SINCE DRY CONVECTION IN THE 
!		     BOUNDARY LAYER IS DONE IN SUBROUTINE SURFCE,
!		     THIS ROUTINE ONLY CHECKS LAYERS 2 TO LM.
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
!	092301	Chien Wang	add bc and oc
!
!	==========================================================


      SUBROUTINE DRYCNV                                                 7501.   

C****                                                                   7502.   
C**** THIS SUBROUTINE MIXES AIR CAUSED BY DRY CONVECTION.  SINCE DRY    7503.   
C**** CONVECTION IN THE BOUNDARY LAYER IS DONE IN SUBROUTINE SURFCE,    7504.   
C**** THIS ROUTINE ONLY CHECKS LAYERS 2 TO LM.                          7505.   
C****                                                                   7506.   

#if ( defined CPL_CHEM || defined TRACERS )
!
#include "chem_para"
#include "chem_com"
!
#endif
#if ( defined TRACERS)
#include "tracer_com.h"
#endif

#include "BD2G04.COM"                                                   7507.   

      COMMON U,V,T,P,Q                                                  7508.   
      COMMON/WORK1/CONV(IM0,JM0,LM0),PK(IM0,JM0,LM0)                    7509.   
      COMMON/WORK2/UT(IM0,JM0,LM0),VT(IM0,JM0,LM0),                     7510.   
     *  RA(8),ID(8),UMS(8)                                              7511.   
      LOGICAL POLE                                                      7512.   
C     DATA RVAP/461.5/                                                  7513.   
      RVX=0.                                                            7514.   
C**** LOAD U,V INTO UT,VT.  UT,VT WILL BE FIXED DURING DRY CONVECTION   7515.   
C****   WHILE U,V WILL BE UPDATED.                                      7516.   
      DO 50 L=1,LM                                                      7517.   
      DO 50 J=2,JM                                                      7518.   
      DO 50 I=1,IM                                                      7519.   
      UT(I,J,L)=U(I,J,L)                                                7520.   
   50 VT(I,J,L)=V(I,J,L)                                                7521.   
C**** OUTSIDE LOOPS OVER J AND I                                        7522.   
      DO 500 J=1,JM                                                     7523.   
      POLE=.FALSE.                                                      7524.   
      IF(J.EQ.1.OR.J.EQ.JM) POLE=.TRUE.                                 7525.   
      IMAX=IM                                                           7526.   
      IF(POLE) IMAX=IM                                                  7527.   
      DO 120 K=1,2                                                      7528.   
      RA(K)=RAPVS(J)                                                    7529.   
 120  RA(K+2)=RAPVN(J)                                                  7530.   
      IM1=IM                                                            7531.   
      DO 500 I=1,IMAX                                                   7532.   
      LMAX=1                                                            7533.   
  130 LMIN=LMAX+1                                                       7534.   
      IF(LMIN.GE.LM) GO TO 500                                          7535.   
      LMAX=LMIN                                                         7536.   
      IF(T(I,J,LMIN)*(1.+Q(I,J,LMIN)*RVX).LE.                           7537.   
     *   T(I,J,LMIN+1)*(1.+Q(I,J,LMIN+1)*RVX)) GO TO 130                7538.   
C**** MIX HEAT AND MOISTURE THROUGHOUT THE UNSTABLE LAYERS              7539.   
      PKMS=PK(I,J,LMIN)*DSIG(LMIN)+PK(I,J,LMIN+1)*DSIG(LMIN+1)          7540.   
      THPKMS=T(I,J,LMIN)*(PK(I,J,LMIN)*DSIG(LMIN))                      7541.   
     *  +T(I,J,LMIN+1)*(PK(I,J,LMIN+1)*DSIG(LMIN+1))                    7542.   
      QMS=Q(I,J,LMIN)*DSIG(LMIN)+Q(I,J,LMIN+1)*DSIG(LMIN+1)             7543.   

#if ( defined CPL_CHEM )
!
! --- 032395
!  sigma of mixing ratios:
!
      cfc11ms=cfc11(i,j,lmin)  *dsig(lmin)
     &       +cfc11(i,j,lmin+1)*dsig(lmin+1)

      cfc12ms=cfc12(i,j,lmin)  *dsig(lmin)
     &       +cfc12(i,j,lmin+1)*dsig(lmin+1)

      xn2oms =xn2o(i,j,lmin)  *dsig(lmin)
     &       +xn2o(i,j,lmin+1)*dsig(lmin+1)

      o3ms   =o3(i,j,lmin)  *dsig(lmin)
     &       +o3(i,j,lmin+1)*dsig(lmin+1)

      coms   =co(i,j,lmin)  *dsig(lmin)
     &       +co(i,j,lmin+1)*dsig(lmin+1)

      zco2ms =zco2(i,j,lmin)  *dsig(lmin)
     &       +zco2(i,j,lmin+1)*dsig(lmin+1)

      xnoms  =xno(i,j,lmin)  *dsig(lmin)
     &       +xno(i,j,lmin+1)*dsig(lmin+1)

      xno2ms =xno2(i,j,lmin)  *dsig(lmin)
     &       +xno2(i,j,lmin+1)*dsig(lmin+1)

      xn2o5ms=xn2o5(i,j,lmin)  *dsig(lmin)
     &       +xn2o5(i,j,lmin+1)*dsig(lmin+1)

      hno3ms =hno3(i,j,lmin)  *dsig(lmin)
     &       +hno3(i,j,lmin+1)*dsig(lmin+1)

      ch4ms  =ch4(i,j,lmin)  *dsig(lmin)
     &       +ch4(i,j,lmin+1)*dsig(lmin+1)

      ch2oms =ch2o(i,j,lmin)  *dsig(lmin)
     &       +ch2o(i,j,lmin+1)*dsig(lmin+1)

      so2ms  =so2(i,j,lmin)  *dsig(lmin)
     &       +so2(i,j,lmin+1)*dsig(lmin+1)

      h2so4ms=h2so4(i,j,lmin)  *dsig(lmin)
     &       +h2so4(i,j,lmin+1)*dsig(lmin+1)

        ! === if hfc, pfc, and sf6 are included:
#ifdef INC_3GASES
	! === 032698
      hfc134ams = hfc134a(i,j,lmin)*dsig(lmin)
     &          + hfc134a(i,j,lmin+1)*dsig(lmin+1)

      pfcms     = pfc(i,j,lmin)*dsig(lmin)
     &          + pfc(i,j,lmin+1)*dsig(lmin+1)

      sf6ms     = sf6(i,j,lmin)*dsig(lmin)
     &          + sf6(i,j,lmin+1)*dsig(lmin+1)
	! ===
#endif

      bcms = bcarbon(i,j,lmin)  *dsig(lmin)
     &     + bcarbon(i,j,lmin+1)*dsig(lmin+1)
      ocms = ocarbon(i,j,lmin)  *dsig(lmin)
     &     + ocarbon(i,j,lmin+1)*dsig(lmin+1)

c 062295
c      h2o2ms =h2o2(i,j,lmin)  *dsig(lmin)
c     &       +h2o2(i,j,lmin+1)*dsig(lmin+1)
!
#endif
#if ( defined TRACERS)
         do ntr=1,ntracers
         tracersms(ntr)=tracers(i,j,lmin,ntr) *dsig(lmin)
     &       +tracers(i,j,lmin+1,ntr)*dsig(lmin+1)
         enddo
#endif

      IF(LMIN+1.GE.LM) GO TO 150                                        7544.   
      TVMS=T(I,J,LMIN)*(1.+Q(I,J,LMIN)*RVX)*(PK(I,J,LMIN)*DSIG(LMIN))   7545.   
     *    +T(I,J,LMIN+1)*(1.+Q(I,J,LMIN+1)*RVX)                         7546.   
     *                                  *(PK(I,J,LMIN+1)*DSIG(LMIN+1))  7547.   
      THETA=TVMS/PKMS                                                   7548.   
      LMINP2=LMIN+2                                                     7549.   
      DO 140 L=LMINP2,LM                                                7550.   
      IF(THETA.LT.T(I,J,L)*(1.+Q(I,J,L)*RVX)) GO TO 160                 7551.   
      PKMS=PKMS+(PK(I,J,L)*DSIG(L))                                     7552.   
      THPKMS=THPKMS+T(I,J,L)*(PK(I,J,L)*DSIG(L))                        7553.   
      QMS=QMS+Q(I,J,L)*DSIG(L)                                          7554.   

#if ( defined CPL_CHEM )
!
! --- sigma of mixing ratios:
!
      cfc11ms=cfc11ms+cfc11(i,j,l)*dsig(l)

      cfc12ms=cfc12ms+cfc12(i,j,l)*dsig(l)

      xn2oms =xn2oms+xn2o(i,j,l)*dsig(l)

      o3ms   =o3ms+o3(i,j,l)*dsig(l)

      coms   =coms+co(i,j,l)*dsig(l)

      zco2ms =zco2ms+zco2(i,j,l)*dsig(l)

      xnoms  =xnoms+xno(i,j,l)*dsig(l)

      xno2ms =xno2ms+xno2(i,j,l)*dsig(l)

      xn2o5ms=xn2o5ms+xn2o5(i,j,l)*dsig(l)

      hno3ms =hno3ms+hno3(i,j,l)*dsig(l)

      ch4ms  =ch4ms+ch4(i,j,l)*dsig(l)

      ch2oms =ch2oms+ch2o(i,j,l)*dsig(l)

      so2ms  =so2ms+so2(i,j,l)*dsig(l)

      h2so4ms=h2so4ms+h2so4(i,j,l)*dsig(l)

        ! === if hfc, pfc, and sf6 are included:
#ifdef INC_3GASES
	! === 032698
      hfc134ams = hfc134ams
     &	        + hfc134a(i,j,l)*dsig(l)

      pfcms     = pfcms
     &	        + pfc(i,j,l)*dsig(l)

      sf6ms     = sf6ms 
     &	        + sf6(i,j,l)*dsig(l)
	! ===
#endif

      bcms = bcms + bcarbon(i,j,l)*dsig(l)
      ocms = ocms + ocarbon(i,j,l)*dsig(l)

c 062295
c      h2o2ms =h2o2ms+h2o2(i,j,l)*dsig(l)

!
#endif
#if ( defined TRACERS)
         do ntr=1,ntracers
         tracersms(ntr)=tracersms(ntr) +tracers(i,j,l,ntr)*dsig(l)
         enddo
#endif

      TVMS=TVMS+T(I,J,L)*(1.+Q(I,J,L)*RVX)*(PK(I,J,L)*DSIG(L))          7555.   
  140 THETA=TVMS/PKMS                                                   7556.   
  150 L=LM+1                                                            7557.   
  160 LMAX=L-1                                                          7558.   
      RDSIGS=1./(SIGE(LMIN)-SIGE(LMAX+1))                               7559.   
      THM=THPKMS/PKMS                                                   7560.   
      QMS=QMS*RDSIGS                                                    7561.   

#if ( defined CPL_CHEM )
!
! --- Get post-transport mixing ratios:
! 
      cfc11ms = cfc11ms*rdsigs

      cfc12ms = cfc12ms*rdsigs

      xn2oms  = xn2oms *rdsigs

      o3ms    = o3ms   *rdsigs

      coms    = coms   *rdsigs

      zco2ms  = zco2ms *rdsigs

      xnoms   = xnoms  *rdsigs

      xno2ms  = xno2ms *rdsigs

      xn2o5ms = xn2o5ms*rdsigs

      hno3ms  = hno3ms *rdsigs

      ch4ms   = ch4ms  *rdsigs

      ch2oms  = ch2oms *rdsigs

      so2ms   = so2ms  *rdsigs

      h2so4ms = h2so4ms*rdsigs

        ! === if hfc, pfc, and sf6 are included:
#ifdef INC_3GASES
	! === 032698
      hfc134ams = hfc134ams*rdsigs

      pfcms     = pfcms*rdsigs

      sf6ms     = sf6ms*rdsigs
	! ===
#endif

      bcms  = bcms*rdsigs
      ocms  = ocms*rdsigs

c 062295
c      h2o2ms  = h2o2ms*rdsigs

!
#endif
#if ( defined TRACERS)
         do ntr=1,ntracers
         tracersms(ntr)=tracersms(ntr) *rdsigs
         enddo
#endif

      DO 180 L=LMIN,LMAX                                                7562.   
         AJL(J,L,12)=AJL(J,L,12)+(THM-T(I,J,L))*PK(I,J,L)*P(I,J)        7563.   
      T(I,J,L)=THM                                                      7564.   
      Q(I,J,L)=QMS                                                      7565.   

#if ( defined CPL_CHEM )
!
! --- Remap mixing ratios:
!
        cfc11(i,j,l)= cfc11ms

        cfc12(i,j,l)= cfc12ms

        xn2o (i,j,l)= xn2oms

        o3   (i,j,l)= o3ms

        co   (i,j,l)= coms

        zco2 (i,j,l)= zco2ms

        xno  (i,j,l)= xnoms

        xno2 (i,j,l)= xno2ms

        xn2o5(i,j,l)= xn2o5ms

        hno3 (i,j,l)= hno3ms

        ch4  (i,j,l)= ch4ms

        ch2o (i,j,l)= ch2oms

        so2  (i,j,l)= so2ms

        h2so4(i,j,l)= h2so4ms

        ! === if hfc, pfc, and sf6 are included:
#ifdef INC_3GASES
	! === 032698
        hfc134a(i,j,l) = hfc134ams

        pfc(i,j,l)     = pfcms

        sf6(i,j,l)     = sf6ms
	! ===
#endif

	bcarbon(i,j,l) = bcms
	ocarbon(i,j,l) = ocms

c 062295
c        h2o2(i,j,l) = h2o2ms
!
#endif     

#if ( defined TRACERS)
         do ntr=1,ntracers
           tracers(i,j,l,ntr) =tracersms(ntr)
         enddo
#endif
180   continue

      IF(POLE) GO TO 300                                                7566.   
C**** MIX MOMENTUM THROUGHOUT UNSTABLE LAYERS AT NON-POLAR GRID BOXES   7567.   
      ID(1)=I+(J-1)*IM                                                  7570.   
      ID(2)=ID(1)+IM*JM*LM                                              7571.   
      ID(3)=I+J*IM                                                      7574.   
      ID(4)=ID(3)+IM*JM*LM                                              7575.   
      DO 240 K=1,4                                                      7576.   
      UMS(K)=0.                                                         7577.   
      DO 220 L=LMIN,LMAX                                                7578.   
  220 UMS(K)=UMS(K)+UT(ID(K),1,L)*DSIG(L)                               7579.   
  240 UMS(K)=UMS(K)*RDSIGS                                              7580.   
      DO 260 L=LMIN,LMAX                                                7581.   
         AJL(J,L,38)=AJL(J,L,38)+(UMS(1)-UT(I,J,L))*.5*                 7582.   
     *     P(I,J)*RA(1)                                                 7583.   
         AJL(J+1,L,38)=AJL(J+1,L,38)+(UMS(3)-                           7584.   
     *     UT(I,J+1,L))*P(I,J)*RA(3)*.5                                 7585.   
      DO 260 K=1,4                                                      7586.   
  260 U(ID(K),1,L)=U(ID(K),1,L)+(UMS(K)-UT(ID(K),1,L))*RA(K)            7587.   
      GO TO 130                                                         7588.   
C**** MIX MOMENTUM THROUGHOUT UNSTABLE LAYERS AT POLAR GRID BOXES       7589.   
  300 JVPO=2                                                            7590.   
      IF(J.EQ.JM) JVPO=JM                                               7591.   
      RAPO=2.*RAPVN(1)                                                  7592.   
      DO 360 IPO=1,IM                                                   7593.   
      UMSPO=0.                                                          7594.   
      VMSPO=0.                                                          7595.   
      DO 320 L=LMIN,LMAX                                                7596.   
      UMSPO=UMSPO+UT(IPO,JVPO,L)*DSIG(L)                                7597.   
  320 VMSPO=VMSPO+VT(IPO,JVPO,L)*DSIG(L)                                7598.   
      UMSPO=UMSPO*RDSIGS                                                7599.   
      VMSPO=VMSPO*RDSIGS                                                7600.   
      DO 340 L=LMIN,LMAX                                                7601.   
      U(IPO,JVPO,L)=U(IPO,JVPO,L)+(UMSPO-UT(IPO,JVPO,L))*RAPO           7602.   
      V(IPO,JVPO,L)=V(IPO,JVPO,L)+(VMSPO-VT(IPO,JVPO,L))*RAPO           7603.   
  340    AJL(JVPO,L,38)=AJL(JVPO,L,38)                                  7604.   
     *  +(UMSPO-UT(IPO,JVPO,L))*P(1,J)*RAPO                             7605.   
  360 CONTINUE                                                          7606.   
      GO TO 130                                                         7607.   
  500 IM1=I                                                             7608.   
      RETURN                                                            7609.   
      END                                                               7610.   
