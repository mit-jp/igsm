
#include "ctrparam.h"

!	==========================================================
!
!	ML_OCEAN.F:  Lots of utility functions.
!
!	----------------------------------------------------------
! 
!	Revision History:
!	
!	When	Who		What
!	----	----------	-------	
!
!	==========================================================


      SUBROUTINE ML_OCEAN(DTOCN,DTATM,KN) 
C****                                                                   1002.
C**** THIS SUBROUTINE 
C****                                                                   1005.

#include "BD2G04.COM"
#include "AGRID.h"
#include "ODIFF.COM"
#include "BDOCEAN.COM"

c     COMMON /QREL/ DFLDSST,RELAX
      COMMON /HF2ML/ HFML(IO0,JM0)
      COMMON/OCN/TG3M(IM0,JM0,12),RTGO(IM0,JM0,LMO),STG3(IM0,JM0),
     &  DTG3(IM0,JM0)
      DIMENSION AMONTH(12),JDOFM(13),JFDOFM(13)                          1009.
      CHARACTER*4 AMONTH                                                1009.1
      dimension ndmonth(12)
      dimension hfluxoan(jm0),hfluxian(jm0),aroceanan(jm0),
     &  aroicenan(jm0)
      character * 80 dirdat1,dirdat2
      character * 4 anmonth
      logical PRNT,FIRST
c     REAL LHM
      data ifirst /1/
      DATA AMONTH/'JAN','FEB','MAR','APR','MAY','JUNE','JULY','AUG',    1010.
     *  'SEP','OCT','NOV','DEC'/                                        1011.
      DATA JDOFM/0,31,59,90,120,151,181,212,243,273,304,334,365/        1012.
      DATA JFDOFM/1,32,60,91,121,152,182,213,244,274,305,335,1/         1012.1
      DATA ndmonth/31,28,31,30,31,30,31,31,30,31,30,31/
      DATA FIRST /.TRUE./
#include "OCNDATA.COM"

C****                                                                   1016.
      if(ifirst.eq.1)then
         TAUO=0.0
         JDAYLAST=-1
         KOCEAN=KN
         DTOCEAN=DTOCN
         STBO=.5672573E-7
         LHM=334000.
         print *,'KOCEAN=',KOCEAN
         print *,'DTOCEAN=',DTOCEAN
         print *,'DTATM=',DTATM
         print *,'With change in ostruc'
         SPERDAY=24.*3600.
         cfevap=SPERDAY/2.499e6
         NSTPPDAY=SPERDAY/DTOCEAN

         DO 310 MONTH=1,JMPERY   
         IF(JDAY.LE.JDOFM(MONTH+1)) GO TO 320  
  310    CONTINUE
  320    CONTINUE

        call init_ocean

      JDAYLAST=JDAY
      HFMLAV=0.0
      QFGAV=0.0
      NQFGAV=0
      QFTAV=0.0
      QFIAV=0.0
      HFAV=0.0
      HFOAV=0.0
      HFIAV=0.0
      EVIAV=0.0
      PRIAV=0.0
      TPRIAV=0.0
      F2DTAV=0.0
      OTDTAV=0.0
      ENRGMEL=0.0
      POICEAV=0.0
      PWATERAV=0.0
      PWTRAV1=0.0
      ENRGIN=0.0
      ENRGI0=0.0
      ENRGT0=0.0
      ENRGAV=0.0
      SSTAV=0.0
      DOTEMAV=0.0
      DTENRG=0.0
      DTENRGO=0.0
      DTENRG1=0.0
      DTENRG2=0.0
      SO=0.0
      PWAV=0.0
      ACE1I=Z1I*RHOI
      DO J=1,JM0
      DO I=1,IO0
       POCEAN=1.-CLAND_OCEAN(i,j)
       SO=SO+POCEAN*DXYP(J)
      if(POCEAN.gt.0.0)then
      TGW=ODATA_OCEAN(I,J,1)
      ROICE=ODATA_OCEAN(I,J,2)
      PWATER=POCEAN*(1.-ROICE)
      ACE2=ODATA_OCEAN(I,J,3)
      ACE=GDATA_ICE(I,J,1)+ACE1I+ACE2                                   8720.
      ENRGI=((GDATA_ICE(I,J,1)+ACE1I)*GDATA_ICE(I,J,2)
     *  +ACE2*GDATA_ICE(I,J,3))*SHI-ACE*LHM   
      WTRO=Z1O(I,J)*RHOW                                                8723.
      WTRW=WTRO-ROICE*ACE                                               8724.
      ENRGW=WTRW*TGW*SHW                                                8725.
      ENRGI0=ENRGI0+ROICE*ENRGI*POCEAN*DXYP(J)
      ENRGTN=ENRGTN+ENRGW*POCEAN*DXYP(J)
      ENRGIN=ENRGIN+(ROICE*ENRGI+ENRGW)*POCEAN*DXYP(J)
      PWAV=PWAV+PWATER*DXYP(J)
      SSTAV=SSTAV+ODATA_OCEAN(I,J,1)*PWATER*DXYP(J)
      DOTEMAV=DOTEMAV+(ODATA_OCEAN(I,J,1)*Z1O(I,J)+ODATA_OCEAN(I,J,4)
     & *(Z12O(I,J)-Z1O(I,J)))/Z12O(I,J)*POCEAN*DXYP(J)
!     print *,i,j,ODATA_OCEAN(I,J,2),ODATA_OCEAN(I,J,1)
      endif
      ENDDO
      ENDDO
       print *,' ENRGIN=',ENRGIN/SO
       print *,' SSTAV=',SSTAV/PWAV
       print *,' DOTEMAV=',DOTEMAV/SO
!      stop
      SSTAV=0.0
      DOTEMAV=0.0
         DO J=1,JM
         DO I=1,IO
          do k=1,5
           ODATA_OCEANAV(I,J,k)=0.0
           if(k.le.3)then
           GDATA_ICEAV(I,J,k)=0.0
           endif
           enddo
           QFLUXTAV(I,J)=0.0
           HF2MLAV(I,J)=0.0
         ENDDO
         ENDDO
        ifirst=0
         print *,'RELAX=',RELAX
        if(KOCEAN.eq.1)then
        if(RELAX)then
          call qflux_rel
        else
          call qflux
        endif
c         if(TRRUN)then
c          print *,'Before ODIFS_2D in mlocean'
c          call  ODIFS_2D
c         endif
        endif
        return
      endif	! ifirst
c
        if(FIRST)then
         print *,'First start of ML model at TOFDAY=',
     &    TOFDAY,' JDAY=',JDAY
         print *,'MONTH=',MONTH,' JYEAR=',JYEAR
         print *,'JMONTH=',JMONTH,' IDAY=',IDAY
         FIRST=.false.
        endif
        if(JDAY.ne.JDAYLAST)then
         print *,'ML started at TOFDAY=',TOFDAY,' JDAY=',JDAY
         DO 210 MONTH=1,JMPERY   
         IF(JDAY.LE.JDOFM(MONTH+1)) GO TO 220  
  210    CONTINUE
  220    CONTINUE
        endif
        if(JDAY.eq.JFDOFM(1).and.
     &     TOFDAY+0.1.LT.DTOCEAN/3600.)then
c        Last call of ocean for a given year, TOFDAY=0, JDAY=1
         PWTRAV1=PWTRAV1/(365.*NSTPPDAY)
         POICEAV1=SO-PWTRAV1
C        POICEAV and PWATERAV are multiplied by DXYP(J) once at each
C        latitude, SO and PWTRAV1 at each grid point
         POCEANAV=POICEAV+PWATERAV
         HFMLAV=HFMLAV/(SO*365.)
         QFGAV=QFGAV/(365.)
         QFTAV=QFTAV/(365.)
         QFIAV=QFIAV/(365.)
         DTENRG=DTENRG/(SO*24.*3600.)
         DTENRGO=DTENRGO/(SO*365.*24.*3600.)
         DTENRG1=DTENRG1/(SO*365.*24.*3600.)
         DTENRG2=DTENRG2/(SO*365.*24.*3600.)
         ENRGMEL=ENRGMEL/(POICEAV*DTOCEAN)
         F2DTAV=F2DTAV/(POICEAV*DTOCEAN)
         OTDTAV=OTDTAV/(PWATERAV*DTOCEAN)
         HFAV=HFAV/(365.*NSTPPDAY*SO*DTOCEAN)
         HFOAV=HFOAV/(365.*NSTPPDAY*PWTRAV1*DTOCEAN)
         HFIAV=HFIAV/(365.*NSTPPDAY*POICEAV1*DTOCEAN)
         EVIAV=EVIAV/(365.*NSTPPDAY*POICEAV1*DTOCEAN)
         PRIAV=PRIAV/(365.*NSTPPDAY*POICEAV1*DTOCEAN)
         TPRIAV=TPRIAV/(365.*NSTPPDAY*POICEAV1*DTOCEAN)
         HFAV1=(HFOAV*PWTRAV1+(F2DTAV+ENRGMEL)*POICEAV1)/SO
      DO J=1,JM0
      DO I=1,IO0
       POCEAN=1.-CLAND_OCEAN(i,j)
      if(POCEAN.gt.0.0)then
      TGW=ODATA_OCEAN(I,J,1)
      ROICE=ODATA_OCEAN(I,J,2)
      PWATER=POCEAN*(1.-ROICE)
      ACE2=ODATA_OCEAN(I,J,3)
      ACE=GDATA_ICE(I,J,1)+ACE1I+ACE2                                   8720.
      ENRGI=((GDATA_ICE(I,J,1)+ACE1I)*(GDATA_ICE(I,J,2))
     *  +ACE2*(GDATA_ICE(I,J,3)))*SHI-ACE*LHM   
      WTRO=Z1O(I,J)*RHOW                                                8723.
      WTRW=WTRO-ROICE*ACE                                               8724.
      ENRGW=WTRW*TGW*SHW                                                8725.
      ENRGIAV=ENRGI0-ROICE*ENRGI*POCEAN*DXYP(J)
      ENRGTAV=ENRGT0-ENRGW*POCEAN*DXYP(J)
      ENRGAV=ENRGAV+(ROICE*ENRGI+ENRGW)*POCEAN*DXYP(J)
      SSTAV=SSTAV+ODATA_OCEAN(I,J,1)*PWATER*DXYP(J)
      DOTEMAV=DOTEMAV+(ODATA_OCEAN(I,J,1)*Z1O(I,J)+ODATA_OCEAN(I,J,4)
     & *(Z12O(I,J)-Z1O(I,J)))/Z12O(I,J)*POCEAN*DXYP(J)
      endif
      ENDDO
      ENDDO
c     ENRGAV=ENRGAV/(SO*365.*24.*3600.)
         write(78,*),'     '
         write(78,*),' DAY=',JDAY,' JYEAR=',JYEAR
         write(78,*),' TOFDAY=',TOFDAY
         write(78,*),' QFGAV=',QFGAV/SO,' NQFGAV=',NQFGAV
c        write(78,*),' OTDTAV=',OTDTAV
         write(78,*),' HFAV=',HFAV,' HFMLAV=',HFMLAV
c        write(78,*),' SO=',SO,' PWTRAV1=',PWTRAV1
c        write(78,*),' POICEAV1=',POICEAV1
         write(78,*),' QFTAV=',QFTAV/SO,' QFIAV=',QFIAV/SO
         write(78,*),' ENRGMEL=',ENRGMEL,' DTENRG=',DTENRG
         write(78,*),' DTENRGO=',DTENRGO
         write(78,*),' DTENRG1=',DTENRG1,' DTENRG2=',DTENRG2
         write(78,*),' F2DTAV=',F2DTAV
c        write(78,*),' HFAV=',HFAV,' HFOAV=',HFOAV
c        write(78,*),' HFIAV=',HFIAV,' HFAV1=',HFAV1
         write(78,*),' EVIAV=',EVIAV,' PRIAV=',PRIAV
         write(78,*),' TPRIAV=',TPRIAV
         write(78,*),' dENRGAV=',(ENRGAV-ENRGIN)/(SO*365.*24.*3600.)
         write(78,*),' SSTAV=',SSTAV/PWTRAV1,' DOTEMAV=',DOTEMAV/SO
c        backspace (78)
c        read(78,*),JDATE,TAU,QFGAV1
         HFMLAV=0.0
         QFGAV=0.0
         NQFGAV=0
         QFTAV=0.0
         QFIAV=0.0
         HFAV=0.0
         HFOAV=0.0
         HFIAV=0.0
         EVIAV=0.0
         PRIAV=0.0
         TPRIAV=0.0
         F2DTAV=0.0
         OTDTAV=0.0
         ENRGMEL=0.0
         DTENRG=0.0
         POICEAV=0.0
         PWATERAV=0.0
         PWTRAV1=0.0
         ENRGIN=ENRGAV
         ENRGAV=0.0
         SSTAV=0.0
         DOTEMAV=0.0
         DTENRG1=0.0
         DTENRG2=0.0
         endif

        do mn=1,12
c       if(JDAY.eq.JDOFM(mn+1).and.
c    &     TOFDAY+0.1.GT.24.-DTOCEAN/3600.)then
        if(JDAY.eq.JFDOFM(mn+1).and.
     &     TOFDAY+0.1.LT.DTOCEAN/3600.)then
         print *,'ML write averages'
         nmn=mn
c        if(mn.eq.1)nmn=12
         print *,AMONTH(nmn)
         DO J=1,JM
         DO I=1,IO
           do k=1,5
           ODATA_OCEANAV(I,J,k)=ODATA_OCEANAV(I,J,k)
     &           /(NSTPPDAY*ndmonth(nmn))
           if(k.le.3)then
           GDATA_ICEAV(I,J,k)=GDATA_ICEAV(I,J,k)
     &          /(NSTPPDAY*ndmonth(nmn))
           endif
           enddo
           QFLUXTAV(I,J)=QFLUXTAV(I,J)
     &        /(NSTPPDAY*ndmonth(nmn))
           HF2MLAV(I,J)=HF2MLAV(I,J)
     &        /(NSTPPDAY*ndmonth(nmn))
         ENDDO
         ENDDO
         if(WRAVER)then
          write(647),jday,jdate,amonth(nmn),
     &      ODATA_OCEANAV,GDATA_ICEAV,QFLUXTAV,HF2MLAV
          nwro=nwro+1
          rewind 648
          write (648) AEXP,nwro
          rewind 648
         endif
         DO J=1,JM
         DO I=1,IO
          do k=1,5
           ODATA_OCEANAV(I,J,k)=0.0
           if(k.le.3)then
           GDATA_ICEAV(I,J,k)=0.0
           endif
           enddo
           QFLUXTAV(I,J)=0.0
           HF2MLAV(I,J)=0.0
         ENDDO
         ENDDO
         go to 556
        endif
        enddo
  556   continue
      if(JDAY.ne.JDAYLAST)then
        if(KOCEAN.eq.1)then
         if(RELAX)then
          call qflux_rel
         else
          call qflux
         endif
        endif
      endif
        call flux_atm2ocean
        call GROUND_OCEAN
         DO J=1,JM
         DO I=1,IO
           QFLUXTAV(I,J)=QFLUXTAV(I,J)+QFLUXT(I,J)
           HF2MLAV(I,J)=HF2MLAV(I,J)+HFML(I,J)
           do k=1,5
           ODATA_OCEANAV(I,J,k)=ODATA_OCEANAV(I,J,k)+
     &           ODATA_OCEAN(I,J,k)
           if(k.le.3)then
           GDATA_ICEAV(I,J,k)=GDATA_ICEAV(I,J,k)
     &          +GDATA_ICE(I,J,k)
           endif
           enddo
         ENDDO
         ENDDO
      if(JDAY.ne.JDAYLAST)then
c        DO J=1,JM
c        DO I=1,IO
c          QFLUXTAV(I,J)=QFLUXTAV(I,J)+QFLUXT(I,J)
c        ENDDO
c        ENDDO
        call daily_for_ocean
        if(KOCEAN.eq.1)then
          call OSTRUC_2D
          if(TRRUN)then
           call  ODIFS_2D
           do j=1,JM
           DTG3(1,J)=0.0
           do l=1,lmo
             RTGO(1,j,l)=0.0
             CONT1=0.0
             do i=1,IO
               PLAND=CLAND_OCEAN(I,J)
               PWATER=1.-PLAND     
               IF(PWATER.GT.0.) THEN
                 CONT1=CONT1+PWATER             
                 RTGO(1,j,l)=RTGO(1,j,l)+PWATER*RTGO2D(I,J,L)
                 if (l.eq.1)DTG3(1,J)=DTG3(1,J)+DTG32D(I,J)*PWATER
               ENDIF
             enddo
             if(CONT1.GT.0.0)then
                RTGO(1,j,l)=RTGO(1,j,l)/CONT1
               if (l.eq.1)DTG3(1,J)=DTG3(1,J)/CONT1
            endif
          enddo
          enddo
c         print *,'mlocean DTG3'
c         print *,DTG3
          endif	! TRRUN
        endif	! KOCEAN
        JDAYLAST=JDAY
      endif	! JDAYLAST
C****                                                                   1049.
      DO 255 J=1,JM                                                     1088.5
      SUM1=0.                                                           1088.51
      SUM2=0.                                                           1088.511
      SUM3=0.                                                           1088.512
      CONT1=0.                                                          1088.52
      DO 256 I=1,IO                                                     1088.53
      PLAND=CLAND_OCEAN(I,J)         
      POICE= ODATA_OCEAN(I,J,2)*(1.-PLAND)    
      PWATER=1.-PLAND                                                   1088.57
      PWATER=1.-PLAND-POICE
      IF(PWATER.LE.0.) GO TO 256                                        1088.58
      CONT1=CONT1+PWATER                                                1088.59
      SUM1=SUM1+PWATER*ODATA_OCEAN(I,J,1)   
      SUM2=SUM2+PWATER*ODATA_OCEAN(I,J,4)  
      SUM3=SUM3+PWATER*ODATA_OCEAN(I,J,5) 
 256  CONTINUE                                                          1088.61
      IF(CONT1.EQ.0.) GO TO 255                                         1088.62
c     IF (J.EQ.1.OR.J.EQ.JM) GO TO 255                                  1088.63
      SUM1=SUM1/CONT1                                                   1088.64
      SUM2=SUM2/CONT1                                                   1088.642
      SUM3=SUM3/CONT1                                                   1088.643
      osst2(J)=SUM2  
      osst3(J)=SUM3 
      osst(J)=SUM1
 255  CONTINUE                                                          1088.67
c     print*,' FROM ML_OCEAN SST'
c     print *,osst
      DO 257 J=1,JM                                                     1088.68
      SUM1=0.                                                           1088.69
      SUM2=0.                                                           1088.7
      SUM3=0.                                                           1088.7
      SUM4=0.                                                           1088.7
      SUM5=0.                                                           1088.7
      CONT1=0.                                                          1088.71
      DO 254 I=1,IO                                                     1088.72
      POICE=ODATA_OCEAN(I,J,2)*(1.-CLAND_OCEAN(I,J))      
      SUM1=SUM1+POICE                                                   1088.74
      SUM2=SUM2+POICE*ODATA_OCEAN(I,J,3)           
      SUM3=SUM3+POICE*GDATA_ICE(I,J,2)                                  1088.75
      SUM4=SUM4+POICE*GDATA_ICE(I,J,1)                                  1088.75
      SUM5=SUM5+POICE*GDATA_ICE(I,J,3)                                  1088.75
 254  CONT1=CONT1+(1.-CLAND_OCEAN(I,J))   
      IF(SUM1.LE.0.) GO TO 425                                          1088.77
      SUM2=SUM2/SUM1                                                    1088.78
      aoice(J)=SUM2 
      tice1(J)=SUM3/SUM1
      tice2(J)=SUM5/SUM1
      snowice(J)=SUM4/SUM1
 425  CONTINUE                                                          1088.81
      IF(CONT1.LE.0.) GO TO 257                                         1088.82
      RATIO=SUM1/CONT1                                                  1088.83
      foice(J)=RATIO   
 257  CONTINUE                                                          1088.86
c        print *,'after 257',jday,tau
c        print *,'GDATA_ICE(20,36,2)',GDATA_ICE(20,36,2)
c     print *,tice1
c        print *,'GDATA_ICE(20,36,3)',GDATA_ICE(20,36,3)
c     print *,tice2
c     stop
      TAUO=TAUO+DTOCEAN/3600.0
c       print *,'ML done at TOFDAY=',TOFDAY,' JDAY=',JDAY
c       print *,' TOFDAY=',TOFDAY,' JDAY=',JDAY
      RETURN                                                            1108.5
C****                                                                   1109.
  901 FORMAT ('0PRESSURE ADDED IN GMP IS',F10.6/)                       1114.
  902 FORMAT ('0MEAN SURFACE PRESSURE OF THE ATMOSPHERE IS',F10.4)      1115.
  910 FORMAT('1',33A4/)                                                 1116.
  915 FORMAT (47X,'DAY',I5,', HR',I3,' (',I2,A5,I5,')',F8.1)            1117.
  920 FORMAT('1')                                                       1118.
      END                                                               1119.
