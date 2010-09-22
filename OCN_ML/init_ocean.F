C $Header: /home/igsm-repository/OCN_ML/init_ocean.F,v 1.1 2010/09/03 19:37:33 jscott Exp $
C $Name:  $

#include "ctrparam.h"

C Initialize Mixed Layer Ocean model

      subroutine init_ocean

#include "BD2G04.COM"   
#include "AGRID.h"
#include "ODIFF.COM"
#include "BDOCEAN.COM"
   
c     common/diff/cfdiff,rkv,EDOZON(JM0)
      common/diff/cfdiff,rkv
      dimension ndmonth(12),amonth(12)
      dimension hfluxoan(jm0),hfluxian(jm0),aroceanan(jm0),
     &  aroicenan(jm0),z1ozon(jm0),roice(jm0),z2ozon(jm0)
      character *4 amonth,jmonthin
      character *4 anmonth
      character * 80 dirdat2,qffile,zmmaxfl,zmavefl,
     & kvnewfl,
     &povocean,aveocean,flto3,nwravo
      namelist /INPUTON/ RELAX,DFLDSST,dirdat2,qffile,zmmaxfl,zmavefl,
     & kvnewfl,
     &povocean,aveocean,flto3,TRRUN,cdiff,nwravo,WRAVER,CTRUN
      DATA ndmonth/31,28,31,30,31,30,31,31,30,31,30,31/
      DATA AMONTH/'JAN','FEB','MAR','APR','MAY','JUNE','JULY'
     &  ,'AUG','SEP','OCT','NOV','DEC'/

      WRAVER=.false.
      CTRUN=.false.
c     print *,'cfdiff from init_ocean=',cfdiff
      open(614,file='name_ocean.dat',status='OLD')
      read(614,INPUTON)

      print *,'nhours=',nhours
C     OCEAN DEPTH FILES
c        dirdat2='/home/sokolov/cc6/01release/INIT/46lat/'
         id2=index(dirdat2," ")
         CLOSE (525)
         open( unit=525,file=dirdat2(1:id2-1)//zmmaxfl,
     *    status='OLD',form='unformatted')
         open( unit=625,file=dirdat2(1:id2-1)//zmavefl,
     *    status='OLD',form='unformatted')
      READ (625) Z1OAV
      CLOSE (625)
      REWIND 525
      READ (525) Z12O                   

C     Q-flux file
      if(.not.RELAX)then
        open(649,file=qffile,
     &  status='OLD',form='unformatted')
      endif
      if(TRRUN)then
         print *, ' Run with heat diffusion into deep ocean'
        open(652,file=flto3,
     &  status='OLD',form='unformatted')
        open(563,file=kvnewfl,
     &  status='OLD',form='unformatted')
      endif
C     POV_OCEAN
      open (658,
     & file=povocean,
     *    status='OLD',form='unformatted')
      read (658),AEXPO
      read (658),tofdayin,jdayin,jdatein,jmonthin,jyearin
      print *,'Read data for ocean model'
      if(ISTRT1.eq.1)then
      if(abs(AEXPO-AEXP).gt.0.1)then
       print *,'Disagreement between AEXP and AEXPO'
       print *,'AEXP=',AEXP,' AEXPO=',AEXPO
       stop
      endif
       print *,'Restart of run=',AEXP
      else
        print *,'Initial data from run=',AEXP
      endif
      print *,tofdainy,jdayin,jdatein,jmonthin,jyearin
      read (658)ODATA_OCEAN,GDATA_ICE,Z1OOLD,OLDICEM
!     print*,'ODATA_OCEAN(i,8,2)'
!     print *,(ODATA_OCEAN(i,8,2),i=1,io0)
!     print*,'ODATA_OCEAN(i,9,2)'
!     print *,(ODATA_OCEAN(i,9,2),i=1,io0)
        DO  J=1,JM    
        DO  I=1,IO          
         CLAND_OCEAN(i,j)=CLAND4OCEAN(i,j)
        ENDDO
        ENDDO

      if(TRRUN)then
       call ODIFS_2D
       if(ISTRT1.eq.1.or.CTRUN)then
        read (658),RTGO2D
       else
        DO  J=1,JM    
        DO  I=1,IO          
        DO  l=1,LMO          
           RTGO2D(I,J,L)=0.0
        ENDDO
        ENDDO
        ENDDO
       endif
      endif
C     AVE_OCEAN
      if(WRAVER)then
      if(ISTRT1.eq.2)then
       open (647,file=aveocean,
     &    status='old',form='unformatted')
       read(647),ANEXP
       open (648,file=nwravo,
     &    status='old',form='unformatted')
        read(648),ANEXP1,nwro
        if(abs(ANEXP-ANEXP1).gt.0.01)then
         print *,'Wrong aveocean and nwravo files'
         print *,ANEXP,ANEXP1
        endif
      else
       open (647,
     & file=aveocean,
     *    status='new',form='unformatted')
       write(647),AEXP
       open (648,file=nwravo,
     &    status='new',form='unformatted')
       nwro=0
      endif
      endif

      DO  J=1,JM    
      DO  I=1,IO          
       if(ODATA_OCEAN(i,j,2).lt.0.0)ODATA_OCEAN(I,J,3)=0.0
       Z1O(I,J)=Z1OOLD(I,J)
       QFLUXT(I,J)=0.0
      ENDDO
      ENDDO

        call daily_for_ocean



      DO 355 J=1,JM                                                     1088.5
      SUM1=0.                                                           1088.51
      SUM2=0.                                                           1088.511
      SUM3=0.                                                           1088.512
      SUM4=0.                                                           1088.512
      SUM5=0.                                                           1088.512
      CONT1=0.                                                          1088.52
      DO 356 I=1,IO                                                     1088.53
      PLAND=CLAND_OCEAN(I,J)
      POICE= ODATA_OCEAN(I,J,2)*(1.-PLAND)                              1088.55
      PWATER=1.-PLAND                                                   1088.57
      PWATER=1.-PLAND-POICE
      IF(PWATER.LE.0.) GO TO 356                                        1088.58
      CONT1=CONT1+PWATER                                                1088.59
      SUM1=SUM1+PWATER*ODATA_OCEAN(I,J,1)                               1088.6
      SUM2=SUM2+PWATER*ODATA_OCEAN(I,J,4)                               1088.601
      SUM3=SUM3+PWATER*ODATA_OCEAN(I,J,5)                               1088.602
      SUM4=SUM4+Z1OOLD(I,J)*PWATER
      SUM5=SUM5+Z12O(I,J)*PWATER
 356  CONTINUE                                                          1088.61
      IF(CONT1.EQ.0.) GO TO 355                                         1088.62
c     IF (J.EQ.1.OR.J.EQ.JM) GO TO 355                                  1088.63
      SUM1=SUM1/CONT1                                                   1088.64
      SUM2=SUM2/CONT1                                                   1088.642
      SUM3=SUM3/CONT1                                                   1088.643
      osst2(J)=SUM2
      osst3(J)=SUM3
      osst(J)=SUM1
      z1ozon(J)=SUM4/CONT1
      z2ozon(J)=SUM5/CONT1
 355  CONTINUE
      DO 357 J=1,JM                                                     1088.68
      SUM1=0.                                                           1088.69
      SUM2=0.                                                           1088.7
      SUM3=0.                                                           1088.7
      SUM4=0.                                                           1088.7
      SUM5=0.                                                           1088.7
      SUM6=0.                                                           1088.7
      CONT1=0.                                                          1088.71
      DO 354 I=1,IO                                                     1088.72
c     if(ODATA_OCEAN(I,J,2).lt.0.05)then
c        ODATA_OCEAN(I,J,2)=0.0
c        ODATA_OCEAN(I,J,3)=0.0
c        GDATA_ICE(I,J,2)=0.0
c        GDATA_ICE(I,J,1)=0.0
c        GDATA_ICE(I,J,3)=0.0
c     endif
      POICE=ODATA_OCEAN(I,J,2)*(1.-CLAND_OCEAN(I,J))
      SUM1=SUM1+POICE                                                   1088.74
      SUM2=SUM2+POICE*ODATA_OCEAN(I,J,3)                                1088.75
      SUM3=SUM3+POICE*GDATA_ICE(I,J,2)                                  1088.75
      SUM4=SUM4+POICE*GDATA_ICE(I,J,1)                                  1088.75
      SUM5=SUM5+POICE*GDATA_ICE(I,J,3)                                  1088.75
      SUM6=SUM6+POICE*ODATA_OCEAN(I,J,2)
 354  CONT1=CONT1+(1.-CLAND_OCEAN(I,J))
      IF(SUM1.LE.0.) GO TO 435                                          1088.77
      SUM2=SUM2/SUM1                                                    1088.78
      aoice(J)=SUM2
      tice1(J)=SUM3/SUM1
      tice2(J)=SUM5/SUM1
      snowice(J)=SUM4/SUM1
      roice(J)=SUM6/SUM1
 435  CONTINUE                                                          1088.81
      IF(CONT1.LE.0.) GO TO 357                                         1088.82
      RATIO=SUM1/CONT1                                                  1088.83
      foice(J)=RATIO
 357  CONTINUE

         do j=1,jm0
         do i=1,IO0
           ODATA(i,j,1)=osst(j)
           ODATA(i,j,3)=aoice(j)
           ODATA(i,j,2)=foice(j)
           GDATA(i,j,1)=snowice(j)
           GDATA(i,j,3)=tice1(j)
           GDATA(i,j,7)=tice2(j)
         enddo
         enddo


       print *,'IN CONDITION'
       print *,'INDAY=',jdayin
       print *,'SST'
!      print 5001,osst
       print *,osst
       print *,'SST_2'
!      print 5001,osst2
       print *,osst2
       print *,'SST_3'
!      print 5001,osst3
       print *,osst3
        print *,'DICE '
        print *,(aoice(j)/916.,j=1,jm)
        print *,'RICE '
        print *,roice
        print *,'FICE '
        print *,foice
        print *,'TICE '
        print *,tice1
        print *,'TICE_2 '
        print *,tice2
       print *,'SNOW'
       print *,snowice
       print *,'Z1OOLD'
       print *,z1ozon
       print *,'Z12O'
       print *,z2ozon
         ifirst= 0
         print *,'GDATA_ICE(20,36,2)',GDATA_ICE(20,36,2)
      return
 5001 format(24f5.1)
       end
