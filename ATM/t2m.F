
      SUBROUTINE TZM(T,TH,Z,Z0,ZGS,PSMB,TG,TS,RIGS,WS,SHEAT,TGM,LR,EPS)
      PARAMETER(R=287.D0,G=9.8D0,CP=1007.D0,CAR=0.35,
     :RCP=R/CP,RG=R/G,GCP=G/CP)
      LOGICAL FIRST
      DIMENSION AROUGH(20),BROUGH(20),CROUGH(20),DROUGH(20),EROUGH(20)  5823.
      DATA  AROUGH/16.59,13.99,10.4,7.35,5.241,3.926,3.126,2.632,2.319, 5824.
     *2.116,1.982,1.893,1.832,1.788,1.757,1.733,1.714,1.699,1.687,1.677/5825.
      DATA BROUGH/3.245,1.733,0.8481,0.3899,0.1832,0.9026E-1,0.4622E-1, 5826.
     * .241E-1,.1254E-1,.6414E-2,.3199E-2,.1549E-2,.7275E-3,.3319E-3,   5827.
     * .1474E-3,.6392E-4,.2713E-4,.1130E-4,.4630E-5,.1868E-5/           5828.
      DATA  CROUGH/5.111,3.088,1.682,.9239,.5626,.3994,.3282,.3017,.299 5829.
     *,.3114,.3324,.3587,.3881,.4186,.4492,.4792,.5082,.5361,.5627,     5830.
     * .5882/                                                           5831.
      DATA DROUGH/1.24,1.02,0.806,0.682,0.661,0.771,0.797,0.895,0.994,  5832.
     * 1.09,1.18,1.27,1.35,1.43,1.50,1.58,1.65,1.71,1.78,1.84/          5833.
      DATA EROUGH/0.128,0.130,0.141,0.174,0.238,0.330,0.438,0.550,0.660,5834.
     * 0.766,0.866,0.962,1.05,1.14,1.22,1.30,1.37,1.45,1.52,1.58/       5835.
      DATA FIRST/.TRUE./
      PS=PSMB*1.D2
      PSR=PS/R
      ROS=PSR/TG

c     DTH=TS-TG
c     CDN=(CAR/LOG(ZGS/Z0))**2
c     IF(RIGS.LE.0) THEN
C     surface layer has unstable stratification
c        DM=SQRT((1.-AROUGH(LR)*RIGS)*(1.-BROUGH(LR)*RIGS)/
c    *     (1.-CROUGH(LR)*RIGS))
c        DH=1.35*SQRT((1.-DROUGH(LR)*RIGS)/(1.-EROUGH(LR)*RIGS))
c     ELSE
C     surface layer has stable stratification
c        DM=1./(1.+(11.238+89.9*RIGS)*RIGS)
c        DH=1.35/(1.+1.93*RIGS)
c     END IF
c     FM=DM
c     FH=DM*DH
c     HF=-CP*ROS*WS*CDN*FH*DTH
c     TGM1=ROS*CDN*FM*WS**2

c     print *,'FROM TZM'
c     print *,'CDN=',CDN
c     print *,'SHEAT=',SHEAT,' HF=',HF
c     print *,'TGM=',TGM,' TGM1=',TGM1
c     print *,Z,Z0,PS,TG,SHEAT,TGM,LR,EPS
CB  CALCULATION OF T (TEMPERATURE ON HEIGHT Z)
      ZZ=Z/Z0
      SZZ=SQRT(ZZ)
      ZZ1=1.D0/ZZ
      SCD=CAR/LOG(ZZ)
      CD=SCD*SCD
      CDSZZ=CD*SZZ
      HF1=SHEAT/CP
      SRCDTX=SCD*SQRT(ROS*TGM)
      THCOF=HF1/SRCDTX
      RICOF=G*Z*HF1*SRCDTX/(TG*TGM**2)
      RICOF2=RICOF**2
      DX0=10.D0
      IF(SHEAT.LT.0.D0) THEN
        DX0=10.D0
      ELSE
        DX0=-10.D0
      END IF
      DX=DX0
      X1=DX
      CALL DR(X1,RICOF2,DEL1,FM,FH,LR)
      in=0
   10 CONTINUE
      in=in+1
      IF(DEL1.LT.0.D0) THEN
         DX=DX+DX0
         X1=X1+DX
         CALL DR(X1,RICOF2,DEL1,FM,FH,LR)
c       print *,'in=',in,'DEL1=',DEL1
      GO TO 10
      END IF
      X0=X1-DX
      in=0
   20 CONTINUE
      in=in+1
      DLX=0.5D0*(X1-X0)
      IF(ABS(DLX).GT.EPS) THEN
        X=X0+DLX
        CALL DR(X,RICOF2,DEL,FM,FH,LR)
        IF(DEL.LT.0.D0) THEN
           X0=X
        ELSE
           X1=X
        END IF
        GO TO 20
      END IF
      TH=TG-THCOF*SQRT(FM)/FH
      PZ=PS-G*ROS*Z
      T=TH/(1.D0+GCP*Z/TG)
CE
      RETURN
      END
      SUBROUTINE DR(RIGS,RICOF2,DEL,FM,FH,LR)
      DIMENSION AROUGH(20),BROUGH(20),CROUGH(20),DROUGH(20),EROUGH(20)  5823.
      DATA  AROUGH/16.59,13.99,10.4,7.35,5.241,3.926,3.126,2.632,2.319, 5824.
     *2.116,1.982,1.893,1.832,1.788,1.757,1.733,1.714,1.699,1.687,1.677/5825.
      DATA BROUGH/3.245,1.733,0.8481,0.3899,0.1832,0.9026E-1,0.4622E-1, 5826.
     * .241E-1,.1254E-1,.6414E-2,.3199E-2,.1549E-2,.7275E-3,.3319E-3,   5827.
     * .1474E-3,.6392E-4,.2713E-4,.1130E-4,.4630E-5,.1868E-5/           5828.
      DATA  CROUGH/5.111,3.088,1.682,.9239,.5626,.3994,.3282,.3017,.299 5829.
     *,.3114,.3324,.3587,.3881,.4186,.4492,.4792,.5082,.5361,.5627,     5830.
     * .5882/                                                           5831.
      DATA DROUGH/1.24,1.02,0.806,0.682,0.661,0.771,0.797,0.895,0.994,  5832.
     * 1.09,1.18,1.27,1.35,1.43,1.50,1.58,1.65,1.71,1.78,1.84/          5833.
      DATA EROUGH/0.128,0.130,0.141,0.174,0.238,0.330,0.438,0.550,0.660,5834.
     * 0.766,0.866,0.962,1.05,1.14,1.22,1.30,1.37,1.45,1.52,1.58/       5835.
      RI=RIGS
c     print *,'here'
c     print *,RI,LR
c     do ii=1,100
c      print *,'  '
c     enddo
      IF(RIGS.LE.0) THEN
C     surface layer has unstable stratification
         DM=SQRT((1.-AROUGH(LR)*RIGS)*(1.-BROUGH(LR)*RIGS)/
     *     (1.-CROUGH(LR)*RIGS))
         DH=1.35*SQRT((1.-DROUGH(LR)*RIGS)/(1.-EROUGH(LR)*RIGS))
      ELSE
C     surface layer has stable stratification
         DM=1./(1.+(11.238+89.9*RIGS)*RIGS)
         DH=1.35/(1.+1.93*RIGS)
      END IF
      FM=DM
      FH=DM*DH
c     print *,FM,FH
      DEL=RI*RI-RICOF2*FM*FM*FM/(FH*FH)
      RETURN
      END
