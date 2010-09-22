C $Header: /home/igsm-repository/OCN_ML/ostruc_2d.F,v 1.1 2010/09/03 19:37:33 jscott Exp $
C $Name:  $

#include "ctrparam.h"

      SUBROUTINE OSTRUC_2D                                              8650.   
C****                                                                   8651.   
C**** THIS SUBROUTINE RESTRUCTURES THE OCEAN TEMPERATURE PROFILE        8652.   
C**** WHEN THE MIXED LAYER DEPTHS ARE CHANGED (NORMALLY DONE ONCE       8653.   
C**** A DAY).                                                           8654.   
C**** THE SUBROUTINE ALSO MELTS ICE WHEN TGO > 0 (C).                   8655.   
C****                                                                   8656.   

#include "BD2G04.COM"
#include "AGRID.h"
#include "ODIFF.COM"
#include "BDOCEAN.COM"

      COMMON /HF2ML/ HFML(IO0,JM0)
      DIMENSION OIDPTH(JM0)                                             8659.1  
      LOGICAL FIRST,PRNT
c     REAL LHM
      DATA FIRST /.TRUE./
#include "OCNDATA.COM"
C****                                                                   8662.   
C**** FDATA  2  LAND COVERAGE (1)                                       8663.   
C****                                                                   8664.   
C**** ODATA  1  OCEAN TEMPERATURE OF FIRST LAYER (C)                    8665.   
C****        2  RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)       8666.   
C****        3  OCEAN ICE AMOUNT OF SECOND LAYER (KG/M**2)              8667.   
C****        4  MEAN OCEAN TEMPERATURE OF SECOND LAYER (C)              8668.   
C****        5  OCEAN TEMPERATURE AT BOTTOM OF SECOND LAYER (C)         8669.   
C****                                                                   8670.   
C**** GDATA_ICE  1  OCEAN ICE SNOW AMOUNT (KG/M**2)           
C****            2  OCEAN ICE TEMPERATURE OF FIRST LAYER (C)   
C****            3  OCEAN ICE TEMPERATURE OF SECOND LAYER (C)   
C****                                                                   8674.   
      if(FIRST)then
        ETA=0.25
        print *,' ETA=',ETA
        print *,' IO=',IO,' JM=',JM
        LHM=334000.
        FIRST=.FALSE.
c       print *,' Z12O from OSTRUC'
c       print *,(Z12O(1,j),j=1,JM)
      endif
      ACE1I=Z1I*RHOI                                                    8675.   
      AC2OIM=RHOI*Z2OIM                                                 8675.1  
C****                                                                   8676.   
C**** RESTRUCTURE OCEAN LAYERS                                          8677.   
C****                                                                   8678.   
c     if(JDATE.eq.16.or.JDATE.eq.1)then
c       print *,JMONTH,'  ',JDATE
c       print *,' Z1OOLD from OSTRUC'
c       print *,(Z1OOLD(1,j),j=1,JM)
c       print *,' Z1O from OSTRUC'
c       print *,(Z1O(1,j),j=1,JM)
c      print *,' T1 ocean'
c      print 5001,(ODATA_OCEAN(1,j,1),j=1,JM)
c      print *,' T2 ocean'
c      print 5001,(ODATA_OCEAN(1,j,4),j=1,JM)
c      print *,' T3 ocean'
c      print 5001,(ODATA_OCEAN(1,j,5),j=1,JM)
c5001 format(24f5.1)
c     endif
      DO 200 J=1,JM                                                     8679.   
      IMAX=IO                                                           8680.   
c     IF(J.EQ.1.OR.J.EQ.JM) IMAX=1                                      8681.   
      DO 200 I=1,IMAX                                                   8682.   
      prnt=.false.
      IF(CLAND_OCEAN(I,J).GE.1.) GO TO 200     
      POCEAN=1.-CLAND_OCEAN(I,J)
      if(prnt)then
    
      write (77,*),'       '
      write(77,*),JMONTH,'  ',JDATE
      write (77,*),'Z1OOLD=',Z1OOLD(I,J),' Z12O=',Z12O(I,J)
      write (77,*),'Z1O=',Z1O(I,J)
      write (77,*),'TOAV=',(ODATA_OCEAN(I,J,1)*Z1OOLD(I,J)+
     &  ODATA_OCEAN(I,J,4)*(Z12O(I,J)-Z1OOLD(I,J)))/Z12O(I,J)
      endif
cnew  IF(Z1OOLD(I,J).GE.Z12O(I,J)) GO TO 140                            8684.   
      IF(Z1O(I,J).EQ.Z1OOLD(I,J)) GO TO 200                             8685.   
      WTR1O=RHOW*Z1O(I,J)-ODATA_OCEAN(I,J,2)*(GDATA_ICE(I,J,1)
     &       +ACE1I+ODATA_OCEAN(I,J,3))                                 8686.   
      DWTRO=RHOW*(Z1O(I,J)-Z1OOLD(I,J))                                 8687.   
      WTR2O=RHOW*(Z12O(I,J)-Z1O(I,J))                                   8688.   
      ENRGOIN=(WTR1O-DWTRO)*ODATA_OCEAN(I,J,1)+
     &  (WTR2O+DWTRO)*ODATA_OCEAN(I,J,4)
      ENRG1IN=(WTR1O-DWTRO)*ODATA_OCEAN(I,J,1)
      ENRG2IN=(WTR2O+DWTRO)*ODATA_OCEAN(I,J,4)
cnew
      T2OLD=ODATA_OCEAN(I,J,4)
      T1OLD=ODATA_OCEAN(I,J,1)
      T3OLD=ODATA_OCEAN(I,J,5)
cnew
      IF(DWTRO.GT.0.) GO TO 120                                         8689.   
C**** MIX LAYER DEPTH IS GETTING SHALLOWER                              8690.   
      ODATA_OCEAN(I,J,4)=ODATA_OCEAN(I,J,4)                             8691.   
     &  +((ODATA_OCEAN(I,J,4)-ODATA_OCEAN(I,J,1))*DWTRO/WTR2O+TTRUNC)   8692.   
      if(prnt)then
      write (77,*),' MIX LAYER DEPTH IS GETTING SHALLOWER'
      write (77,*),' ODATA_OCEAN(I,J,4)=', ODATA_OCEAN(I,J,4)
      write (77,*),' ODATA_OCEAN(I,J,1)=', ODATA_OCEAN(I,J,1)
      endif
      GO TO 210                                                         8693.   
C**** MIX LAYER DEPTH IS GETTING DEEPER                                 8694.   
  120 TGAVE=(ODATA_OCEAN(I,J,4)*DWTRO+
     &    (2.*ODATA_OCEAN(I,J,4)-ODATA_OCEAN(I,J,5))*WTR2O)             8695.   
     &  /(WTR2O+DWTRO)                                                  8696.   
cnew  ODATA_OCEAN(I,J,1)=ODATA_OCEAN(I,J,1)+
cnew &((TGAVE-ODATA_OCEAN(I,J,1))*DWTRO/WTR1O                           8697.   
cnew *  +TTRUNC)                                                        8698.   
      IF(Z1O(I,J).GE.Z12O(I,J)) GO TO 140                               8699.   
      ODATA_OCEAN(I,J,4)=ODATA_OCEAN(I,J,4)                             8700.   
     *  +((ODATA_OCEAN(I,J,5)-ODATA_OCEAN(I,J,4))*DWTRO/
     &   (WTR2O+DWTRO)+TTRUNC)                                          8701.   
cnew
      ODATA_OCEAN(I,J,1)=(ODATA_OCEAN(I,J,1)*(WTR1O-DWTRO)+
     & T2OLD*(WTR2O+DWTRO)-ODATA_OCEAN(I,J,4)*WTR2O)/WTR1O 
cnew
      if(prnt)then
      write (77,*),' MIX LAYER DEPTH IS GETTING DEEPER'
      write (77,*),' TGAVE=', TGAVE,' T2OLD=',T2OLD
      write (77,*),' ODATA_OCEAN(I,J,1)=', ODATA_OCEAN(I,J,1)
      write (77,*),' ODATA_OCEAN(I,J,4)=', ODATA_OCEAN(I,J,4)
      endif
      GO TO 210                                                         8702.   
C**** MIXED LAYER DEPTH IS AT ITS MAXIMUM OR TEMP PROFILE IS UNIFORM    8703.   
  140 CONTINUE
      TML=((WTR1O-DWTRO)*ODATA_OCEAN(I,J,1)+
     &  (WTR2O+DWTRO)*ODATA_OCEAN(I,J,4))/WTR1O
      if(prnt)then
      write (77,*),'MIXED LAYER DEPTH IS AT ITS MAXIMUM'
      write (77,*),' TML=', TML
      write (77,*),' ODATA_OCEAN(I,J,1)=', ODATA_OCEAN(I,J,1)
      endif
      ODATA_OCEAN(I,J,1)=TML
      ODATA_OCEAN(I,J,4)=ODATA_OCEAN(I,J,1)   
      ODATA_OCEAN(I,J,5)=ODATA_OCEAN(I,J,1)  
  210 CONTINUE                                                          8706.   
      DENR1=SHW*WTR1O*ODATA_OCEAN(I,J,1)-SHW*T1OLD*(WTR1O-DWTRO)
      DENR2=SHW*WTR2O*ODATA_OCEAN(I,J,4)-SHW*T2OLD*(WTR2O+DWTRO)
      HFML(I,J)=DENR1/(24.*3600.)
      if(prnt)then
c     write (77,*),' DENR1=', DENR1/(24.*3600.)
c    &   ,' DENR2=', DENR2/(24.*3600.)
      write (77,*),' T1OLD=', T1OLD,' T2OLD=', T2OLD
     &  ,' T3OLD=', T3OLD
      write (77,*),' ODATA_OCEAN(I,J,1)=', ODATA_OCEAN(I,J,1)
      write (77,*),' ODATA_OCEAN(I,J,4)=', ODATA_OCEAN(I,J,4)
      write (77,*),' ODATA_OCEAN(I,J,5)=', ODATA_OCEAN(I,J,5)
      write (77,*),'TONAV=',(ODATA_OCEAN(I,J,1)*Z1O(I,J)+
     &  ODATA_OCEAN(I,J,4)*(Z12O(I,J)-Z1O(I,J)))/Z12O(I,J)
      endif
      ENRGOOUT=WTR1O*ODATA_OCEAN(I,J,1)+
     &  WTR2O*ODATA_OCEAN(I,J,4)
      ENRG1OUT=WTR1O*ODATA_OCEAN(I,J,1)
      ENRG2OUT=WTR2O*ODATA_OCEAN(I,J,4)
      DTENRGO=DTENRGO+SHW*(ENRGOOUT-ENRGOIN)*POCEAN*DXYP(J)
      DTENRG1=DTENRG1+SHW*(ENRG1OUT-ENRG1IN)*POCEAN*DXYP(J)
      DTENRG2=DTENRG2+SHW*(ENRG2OUT-ENRG2IN)*POCEAN*DXYP(J)
      QTEM=ODATAT(I,J,1)-T1OLD
      DTO1=ODATA_OCEAN(I,J,1)-T1OLD
  200 CONTINUE                                                          8706.   
C****                                                                   8707.   
C**** REDUCE THE HORIZONTAL EXTENT OF ICE IF OCEAN TEMPERATURE IS WARM  8708.   
C****                                                                   8709.   
      DO 300 J=1,JM                                                     8710.   
      IMAX=IO                                                           8711.   
c     IF(J.EQ.1.OR.J.EQ.JM) IMAX=1                                      8712.   
      DO 300 I=1,IMAX                                                   8713.   
      prnt=.false.
      IF(ODATA_OCEAN(I,J,2).LE.0.) GO TO 300 
      IF(CLAND_OCEAN(I,J).GE.1.) GO TO 300      
      POCEAN=1.-CLAND_OCEAN(I,J)
      if(prnt)then
      write (77,*),' REDUCE THE HORIZONTAL EXTENT OF ICE'
      write (77,*),' ODATA_OCEAN(I,J,1)=', ODATA_OCEAN(I,J,1)
      write (77,*),' ODATA_OCEAN(I,J,2)=', ODATA_OCEAN(I,J,2)
      endif
C**** REDUCE ICE EXTENT IF OCEAN TEMPERATURE IS GREATER THAN 0 DEGREES  8716.   
      IF(ODATA_OCEAN(I,J,1).LE.0.) GO TO 300  
      TGW=ODATA_OCEAN(I,J,1)                 
      ROICE=ODATA_OCEAN(I,J,2)              
      ACE2=ODATA_OCEAN(I,J,3)              
      ACE=GDATA_ICE(I,J,1)+ACE1I+ACE2                                   8720.   
      ENRGI=((GDATA_ICE(I,J,1)+ACE1I)*GDATA_ICE(I,J,2)                  8721.   
     *  +ACE2*GDATA_ICE(I,J,3))*SHI-ACE*LHM                             8722.   
      WTRO=Z1O(I,J)*RHOW                                                8723.   
      WTRW=WTRO-ROICE*ACE                                               8724.   
      ENRGW=WTRW*TGW*SHW                                                8725.   
      ENRGIN=(ROICE*ENRGI+ENRGW)*POCEAN*DXYP(J)
      if(prnt)then
      write (77,*),' ODATA_OCEAN(I,J,3)=', ODATA_OCEAN(I,J,3)
      write (77,*),' SNOW=', GDATA_ICE(I,J,1)
      write (77,*),' TICE1=', GDATA_ICE(I,J,2)
      write (77,*),' TICE2=', GDATA_ICE(I,J,3)
      write (77,*),' ACE=', ACE
      write (77,*),' ROICE*ENRGI=', ROICE*ENRGI,' ENRGW=',ENRGW
      endif
      IF(ROICE*ENRGI+ENRGW.LT.0.) GO TO 230                             8726.   
C**** THE WARM OCEAN MELTS ALL THE SNOW AND ICE                         8727.   
      ODATA_OCEAN(I,J,1)=(ROICE*ENRGI+ENRGW)/(WTRO*SHW)   
      if(prnt)then
      write (77,*),' ODATA_OCEAN(I,J,1)=', ODATA_OCEAN(I,J,1)
      endif
      GO TO 270                                                         8729.   
C**** THE WARM OCEAN COOLS TO 0 DEGREES MELTING SOME SNOW AND ICE       8730.   
  230 ODATA_OCEAN(I,J,1)=0. 
C**** THE FRACTION ETA OF ENERGY IS USED TO REDUCE THE ICE DEPTH        8731.1  
      DACE2=ETA*ENRGW/(ROICE*(GDATA_ICE(I,J,3)*SHI-LHM))                8731.2  
      IF(DACE2.LT.AC2OIM-ACE2) DACE2=AC2OIM-ACE2                        8731.3  
      DENRGI=DACE2*(GDATA_ICE(I,J,3)*SHI-LHM)                           8731.4  
      ROICENEW=(ROICE*ENRGI+ENRGW)/(ENRGI+DENRGI)
      ACE2NEW=ACE2+DACE2
      if(ROICENEW.lt.0.0.or.ACE2NEW.lt.AC2OIM)then
        print *,'From ostruc'
        print *,'J=',j,' ROICENEW=',ROICENEW,' ACE2NEW=',ACE2NEW
        stop
      endif
      ODATA_OCEAN(I,J,2)=(ROICE*ENRGI+ENRGW)/(ENRGI+DENRGI) 
      ODATA_OCEAN(I,J,3)=ACE2+DACE2                        
      if(prnt)then
      write (77,*),' ODATA_OCEAN(I,J,2)=', ODATA_OCEAN(I,J,2)
      write (77,*),' ODATA_OCEAN(I,J,3)=', ODATA_OCEAN(I,J,3)
      endif
      GO TO 290                                                         8733.   
  270 ODATA_OCEAN(I,J,2)=0.                              
      ODATA_OCEAN(I,J,3)=0.                             
      GDATA_ICE(I,J,1)=0.                                               8736.   
      GDATA_ICE(I,J,3)=0.                                               8737.   
      GDATA_ICE(I,J,2)=0.                                               8738.   
      ROICE=ODATA_OCEAN(I,J,2)              

  290 CONTINUE  
      TGW=ODATA_OCEAN(I,J,1)                 
      ROICE=ODATA_OCEAN(I,J,2)              
      ACE2=ODATA_OCEAN(I,J,3)              
      ACE=GDATA_ICE(I,J,1)+ACE1I+ACE2                                   8720.   
      ENRGI=((GDATA_ICE(I,J,1)+ACE1I)*GDATA_ICE(I,J,2)                  8721.   
     *  +ACE2*GDATA_ICE(I,J,3))*SHI-ACE*LHM                             8722.   
      WTRO=Z1O(I,J)*RHOW                                                8723.   
      WTRW=WTRO-ROICE*ACE                                               8724.   
      ENRGW=WTRW*TGW*SHW                                                8725.   
      ENRGOUT=(ROICE*ENRGI+ENRGW)*POCEAN*DXYP(J)
      DTENRG=DTENRG+ENRGOUT-ENRGIN
      if(prnt)then
      write (77,*),' DTENRG=', DTENRG/(24.*3600.)
      write (77,*),' ODATA_OCEAN(I,J,2)=', ODATA_OCEAN(I,J,2)
      write (77,*),' ODATA_OCEAN(I,J,3)=', ODATA_OCEAN(I,J,3)
      endif
  300 CONTINUE                                                          8739.   
C     PRINT OCEAN ICE DEPH                                              8739.1  
c        DO 800 J=1,JM                                                  8739.2  
c 800    OIDPTH(J)=ODATA_OCEAN(1,J,3)/RHOI            
c        IF(DMOD(TAU,720.).EQ.0.) WRITE(6,900) TAU,OIDPTH               8739.4  
c 900 FORMAT (1X,'TAU=',F12.0,' OCN ICE DEPTH OF 2ND LAYER=',/1X,24F5.1)8739.5  
c     if(JDATE.eq.16.or.JDATE.eq.1)then
c      print *,' after OSTRUC'
c      print *,' T1 ocean'
c      print 5001,(ODATA_OCEAN(1,j,1),j=1,JM)
c      print *,' T2 ocean'
c      print 5001,(ODATA_OCEAN(1,j,4),j=1,JM)
c      print *,' T3 ocean'
c      print 5001,(ODATA_OCEAN(1,j,5),j=1,JM)
c     endif
      RETURN                                                            8740.   
      END                                                               8741.   
