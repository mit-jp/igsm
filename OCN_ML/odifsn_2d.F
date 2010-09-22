C $Header: /home/igsm-repository/OCN_ML/odifsn_2d.F,v 1.1 2010/09/03 19:37:33 jscott Exp $
C $Name:  $

#include "ctrparam.h"

c	==========================================================
c
c	ODIFSN.F:  THIS SUBROUTINE CALCULATES THE ANNUAL OCEAN 
c			TEMPERATURE AT THE MAXIMUM MIXED LAYER, 
c			COMPARES THAT TO THE CONTROL RUN'S
c			TEMPERATURE, CALLS SUBROUTINE DIFFUS, 
c			AND REDUCES THE UPPER OCEAN TEMPERATURES 
c			BY THE AMOUNT OF HEAT THAT IS DIFFUSED INTO
c			THE THERMOCLINE
c
c	----------------------------------------------------------
c 
c	Revision History:
c	
c	When	Who		What
c	-----	----------	-------	
c	080100	Chien Wang	repack based on CliChem3 & M24x11,
c				  and add cpp.
c
c	==========================================================


      SUBROUTINE ODIFS_2D
C****                                                                      2.
C**** THIS SUBROUTINE CALCULATES THE ANNUAL OCEAN TEMPERATURE AT THE       3.
C**** MAXIMUM MIXED LAYER, COMPARES THAT TO THE CONTROL RUN'S              4.
C**** TEMPERATURE, CALLS SUBROUTINE DIFFUS, AND REDUCES THE UPPER          5.
C**** OCEAN TEMPERATURES BY THE AMOUNT OF HEAT THAT IS DIFFUSED INTO       6.
C**** THE THERMOCLINE                                                      7.
C****                                                                      8.
CINCL 001                                                                  9.

#include "BD2G04.COM"
#include "AGRID.h"
#include "ODIFF.COM"
#include "BDOCEAN.COM"

c     COMMON/OCN_2D/TG3M2D(IO0,JM0,0:13),RTGO2D(IO0,JM0,LMO),
c    &  DTG32D(IO0,JM0),TG3MA2D(IO0,JM0),TG3MAZ(JM0)
C****                                                                     13.
C**** ACCUMULATE OCEAN TEMPERATURE AT MAXIMUM MIXED LAYER                 14.
C****                                                                     15.
      DIMENSION EDO(IO0,JM0),ADTG3(IO0,JM0),PWATER(IO0,JM0)               16.
     & ,NDAYS(12),TG3MAZ(JM0),EDONEWZ(JM0),EDHZON(JM0)
     & ,ED1(IO0,JM0)
      common/diff/cfdiff,rkv

       common/Garydiff/depthml(jm0),edzcar(jm0),
     & dzg(lmo),dzog(lmo-1),
     & Rco2(jm0,lmo)
       common /Garyvdif/iyearocm,vdfocm,acvdfc,cfocdif,diffcar0
       common /Garyvlog/odifcarbon,ocarcont
       real edzon(jm0)
       logical odifcarbon,HPRNT,ocarcont

      DATA IFIRST/1/                                                      17.
      data NDAYS/31,28,31,30,31,30,31,31,30,31,30,31/
      HPRNT=TAU.ge.17520.0.and.TAU.lt.17545.0
      HPRNT=.false.
      if(HPRNT)then
        print *,' From ODIFS_2D TAU=',TAU,' IFIRST=',IFIRST
      endif
      IF (IFIRST.NE.1) GO TO 90                                           18.
      print *,'ODIFS_2D cfdiff=',cfdiff,' rkv=',rkv
      IFIRST=0                                                            19.
      READ (562) EDO                                                      20.
      READ (563) EDONEWZ
c     EDO - GISS coefficiens base on tritium mixing
c     EDONEWZ - coefficiens based on heat uptake
      do  j=1,JM
        EDONEWZ(j)=cfdiff*EDONEWZ(j)
        do  i=1,IO
         EDO(i,j)=cfdiff*EDO(i,j)
        enddo
      enddo

      READ (652)runt3,TG3M2D
      print *,' TG3 from run ',runt3
c     print *,'  TG3M2D'
c     print 5001,(TG3M2D(1,j,1),j=1,JM)
      DO 50 J=1,JM                                                        22.
      SUM1=0.                                                             23.
      CONT1=0.                                                            24.
      DO 40 I=1,IO                                                        25.
      PWAT=1.-CLAND_OCEAN(I,J)                                                 26.
      IF(PWAT.LE.0.) GO TO 40                                             27.
      CONT1=CONT1+PWAT                                                    28.
      SUM1=SUM1+EDO(I,J)*PWAT                                             29.
   40 CONTINUE                                                            30.
      IF(CONT1.GT.0.) THEN
      edzon(j)= SUM1/CONT1
      ELSE
      edzon(j)= 0.0
      ENDIF
      edhzon(j)= rkv*edzon(j)+(1.-rkv)*edonewz(j)
   50 CONTINUE

c     edzon -zonal average of GISS coefficiens base on tritium mixing
c     edonewz -zonal average of coefficiens based on heat uptake
c     rkv - weight of GISS coefficients in total coefficients
c     edhzon -total zonal average of vertical diffusion coefficiens

      do j=1,jm
       do i=1,io
        if(edzon(j).gt.0.0) then
         ED1(i,j)=edo(i,j)/edzon(j)
         EDO(i,j)=ED1(i,j)*edhzon(j)
        endif
       enddo
       enddo

c     EDO -total vertical diffusion coefficiens for heat

C zonally avegaged EDO 

      print *,'  Vertical diffusion coefficients for heat'
      print *,(edhzon(j),j=1,JM)

      if(odifcarbon)then
!       Kvh=cfdiff*2.5
!       Kvc=diffcar0+cfocdif*Kvh
        acvdfc=diffcar0+cfdiff*cfocdif*2.5
        cfocnew=acvdfc/(cfdiff*2.5)
        do j=1,jm
!       do l=1,lmo
!        edzcar(j,l)=cfocnew*edzon(j)
         edzcar(j)=cfocnew*edhzon(j)
!       enddo
        enddo
!       10/28/06
!       ocarcont = true: continuation of run with ocean carbon model
!       ocarcont = false: run with ocean carbon model from equilibrium
        if (.not.ocarcont) vdfocm=acvdfc
        if(abs(acvdfc-vdfocm).gt.0.0001)then
          print *,'Wrong vertical diffusion coefficeint for carbon'
          print *,' From initial data file=',vdfocm
          print *,' From name.dat file =',acvdfc
          stop
        else
          print *,'Vertical diffusion coefficeint for carbon=',acvdfc
          print *,'diffcar0=',diffcar0,' cfocdif=',cfocdif
          print *,'Kcv/Khv=',cfocnew
          print *,(edzcar(j),j=1,JM)
        endif
      endif	! odifcarbon

c     print *,' DTG3/356'
c     print 5001,(DTG3(1,j)/365.,j=1,JM)
c     print *,' RTGO'
c     print 5001,((RTGO(1,j,k),j=1,JM),k=1,LMO)
 5001 format(24f5.1)
      do J=1,JM
        SUM1=0.                 
        CONT1=0.                 
        do 25 I=1,IO
          PWAT=1.-CLAND_OCEAN(I,J)  
          TG3MA2D(I,J)=0.
          do 26 m=1,12
            TG3MA2D(I,J)=TG3MA2D(I,J)+TG3M2D(I,J,m)/12.
   26     continue
          IF(PWAT.GT.0.) THEN
            CONT1=CONT1+PWAT 
            SUM1=SUM1+TG3MA2D(I,J)*PWAT   
          ENDIF
   25   continue
        IF(CONT1.GT.0.) THEN
          TG3MAZ(j)= SUM1/CONT1
        ENDIF
      enddo
      print *,'JM=',JM,' IO=',IO,'  CONT1=', CONT1
      print *,'  TG3MAZ'
      print 5001,(TG3MAZ(j),j=1,JM)
      return
   90 CONTINUE                                                            39.
C****                                                                     45.
C**** AT THE END OF EACH MONTH, UPDATE THE OCEAN TEMPERATURE              46.
C**** DIFFERENCE AND REPLACE THE MONTHLY SUMMED TEMPERATURE               47.
C****                                                                     48.
      DO 210 J=1,JM                                                       52.
      IMAX=IO                                                             53.
      IF((J.EQ.1).OR.(J.EQ.JM)) IMAX=1                                    54.
      DO 210 I=1,IMAX                                                     55.
      DTG32D(I,J)=ODATA_OCEAN(I,J,5)-TG3MA2D(I,J)
  210 CONTINUE
      if(HPRNT)then
        print *,'TAU=',TAU
        print *,'ODIFS_2D DTG32D'
        print *,(DTG32D(1,J),J=1,JM)
        print *,'ODIFS_2D TG3MA2D'
        print *,(TG3MA2D(1,J),J=1,JM)
        print *,'ODIFS_2D ODATA_OCEAN'
        print *,(ODATA_OCEAN(1,J,5),J=1,JM)
      endif
 945  format(6f10.3)
C****                                                                     59.
C**** DIFFUSE THE OCEAN TEMPERATURE DIFFERENCE OF THE UPPER LAYERS        60.
C**** INTO THE THERMOCLINE AND REDUCE THE UPPER TEMPERATURES BY THE       61.
C**** HEAT THAT IS DIFFUSED DOWNWARD                                      62.
C****                                                                     63.
  300 DO 310 J=1,JM                                                       64.
      IMAX=IO                                                             65.
      IF((J.EQ.1).OR.(J.EQ.JM)) IMAX=1                                    66.
      DO 310 I=1,IMAX                                                     67.
      ADTG3(I,J)=DTG32D(I,J)               
  310 RTGO2D(I,J,1)=ADTG3(I,J)                             
C****                                                                     70.
      DO 320 J=1,JM                                                       71.
      DO 320 I=1,IO                                                       72.
  320 PWATER(I,J)=1.-CLAND_OCEAN(I,J)
C
c      print *,'ODIFS_2D TAU=',TAU,' JDATE=',JDATE,' MONTH=',MONTH
c      print *,' TG3MA2D=',TG3MA2D(1,jm/2),' TO3=',ODATA_OCEAN(1,JM/2,5)
c      print *,' DTG32D=',DTG32D(1,JM/2)
c      print *,' RTGO2D=',(RTGO2D(1,JM/2,l),l=1,LMO)
      CALL DIFFUS2D (IO,JM,SDAY,.5,EDO,Z12O,PWATER,RTGO2D) 
c     print *,'After dffus'
c      print *,' RTGO2D=',(RTGO2D(1,JM/2,l),l=1,LMO)
c      print *,'DTO3=',(RTGO2D(1,JM/2,1)-ADTG3(1,JM/2))
C****                                                                     75.
      DO 340 J=1,JM                                                       76.
         ADT=0.                                                           77.
         CDT=0.                                                           78.
      IMAX=IO                                                             79.
      IF((J.EQ.1).OR.(J.EQ.JM)) IMAX=1                                    80.
      DO 330 I=1,IMAX                                                     81.
      IF(PWATER(I,J).LE.0.) GO TO 330                                     82.
      ODATA_OCEAN(I,J,1)=ODATA_OCEAN(I,J,1)+(RTGO2D(I,J,1)-ADTG3(I,J))         
      ODATA_OCEAN(I,J,4)=ODATA_OCEAN(I,J,4)+(RTGO2D(I,J,1)-ADTG3(I,J))        
      ODATA_OCEAN(I,J,5)=ODATA_OCEAN(I,J,5)+(RTGO2D(I,J,1)-ADTG3(I,J))       
         ADT=ADT+(RTGO2D(I,J,1)-ADTG3(I,J))*Z12O(I,J)*PWATER(I,J)   
     *     *(1.-ODATA_OCEAN(I,J,2))       
         CDT=CDT+(RTGO2D(I,J,1)-ADTG3(I,J))*Z12O(I,J)*PWATER(I,J)    
     *     *ODATA_OCEAN(I,J,2)  
  330 CONTINUE                                                            90.
         AJ(J,68)=AJ(J,68)-ADT                                            91.
         CJ(J,68)=CJ(J,68)-CDT                                            92.
  340 CONTINUE                                                            93.
      RETURN                                                              94.
      END                                                                 95.

#include "ctrparam.h"

      SUBROUTINE DIFFUS2D (IO,JM,DT,ALPHA,ED,Z12O,PWATER,R) 
C****                                                                     97.
C**** THIS SUBROUTINE CALCULATES THE VERTICAL MIXING OF A TRACER, R,      98.
C**** BY DIFFUSION.  LM IS THE NUMBER OF VERTICAL LAYERS.  DT (S) IS      99.
C**** THE TIME STEP.  ALPHA DETERMINES THE TIME SCHEME RANGING FROM      100.
C**** 0 FOR EXPLICIT TO 1 FOR FULLY IMPLICIT.  DZ (M) IS THE DEPTH OF    101.
C**** THE LAYERS, AND DZO (M) IS THE DISTANCE BETWEEN THE CENTERS OF     102.
C**** THE LAYERS.  ED (M**2/S) IS THE DIFFUSION COEFFICIENT BETWEEN      103.
C**** ADJACENT LAYERS.  R IS THE TRACER CONCENTRATION.                   104.
C****                                                                    105.

#include "ODIFF.COM"

      DIMENSION ED(IO,JM),Z12O(IO,JM),PWATER(IO,JM),R(IO,JM,LMO)
      DIMENSION DZ(LMO),DZO(LMO1),AM(LMO),BM(LMO),CM(LMO),DM(LMO)
C**** DEFINE THE VERTICAL LAYERING                                       108.
c     print *,'From diffus ifirst=',ifirst
      data ifirst /1/
      LM=LMO                                                             109.
      LMM1=LM-1                                                          110.
      if(ifirst.eq.1)then
      DZO(1)=10./SQRT(1.7010587)                                         111.
      DZ(2)=10.                                                          112.
      DO 10 L=2,LMM1                                                     113.
      DZO(L)=DZO(L-1)*1.7010587                                          114.
   10 DZ(L+1)=DZ(L)*1.7010587                                            115.
      ZMAX=100.
      do l=1,LM
       zmax=zmax+DZ(l)
      enddo
      print *,'DZ'
      print *,DZ
      print *,'DZO'
      print *,DZO
      print *,' ZMAX_OCEAM=',zmax
      ifirst=0
      endif
C**** LOOP OVER THE HORIZONTAL GRID POINTS                               116.
  100 DO 500 J=1,JM                                                      117.
      IMAX=IO                                                            118.
      IF((J.EQ.1).OR.(J.EQ.JM)) IMAX=1                                   119.
      DO 500 I=1,IMAX                                                    120.
      IF(PWATER(I,J).LE.0.) GO TO 500                                    121.
      DZ(1)=Z12O(I,J)                                                    122.
C**** SET UP TRIDIAGONAL MATRIX ENTRIES AND RIGHT HAND SIDES             123.
      AM(1)=DZ(1)+ALPHA*DT*ED(I,J)/DZO(1)                                124.
      DM(1)=DZ(1)*R(I,J,1)                                               125.
     *  -(1.-ALPHA)*DT*ED(I,J)*(R(I,J,1)-R(I,J,2))/DZO(1)                126.
      L=1                                                                127.
  110 BM(L)=-ALPHA*DT*ED(I,J)/DZO(L)                                     128.
      CM(L+1)=-ALPHA*DT*ED(I,J)/DZO(L)                                   129.
      IF(L.GE.LM-1) GO TO 120                                            130.
      L=L+1                                                              131.
      AM(L)=DZ(L)+ALPHA*DT*(ED(I,J)/DZO(L-1)+ED(I,J)/DZO(L))             132.
      DM(L)=DZ(L)*R(I,J,L)                                               133.
     *  +(1.-ALPHA)*DT*(ED(I,J)*(R(I,J,L-1)-R(I,J,L))/DZO(L-1)           134.
     *  -ED(I,J)*(R(I,J,L)-R(I,J,L+1))/DZO(L))                           135.
      GO TO 110                                                          136.
  120 AM(LM)=DZ(LM)+ALPHA*DT*ED(I,J)/DZO(LM-1)                           137.
      DM(LM)=DZ(LM)*R(I,J,LM)                                            138.
     *  +(1.-ALPHA)*DT*ED(I,J)*(R(I,J,LM-1)-R(I,J,LM))/DZO(LM-1)         139.
C**** ELIMINATE LOWER OFF-DIAGONAL ENTRIES FROM THE MATRIX EQUATION      140.
      DO 210 L=2,LM                                                      141.
      AM(L)=AM(L)-BM(L-1)*CM(L)/AM(L-1)                                  142.
  210 DM(L)=DM(L)-DM(L-1)*CM(L)/AM(L-1)                                  143.
C**** ELIMINATE UPPER OFF-DIAGONAL ENTRIES FROM THE MATRIX EQUATION      144.
      DO 310 LX=2,LM                                                     145.
      L=1+LM-LX                                                          146.
  310 DM(L)=DM(L)-DM(L+1)*BM(L)/AM(L+1)                                  147.
C**** SOLVE THE REDUCED MATRIX EQUATION TO CALCULATE THE NEW R           148.
      DO 410 L=1,LM                                                      149.
  410 R(I,J,L)=DM(L)/AM(L)                                               150.
  500 CONTINUE                                                           151.
      RETURN                                                             152.
      END                                                                153.
