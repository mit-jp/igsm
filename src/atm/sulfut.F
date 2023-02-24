
#include "ctrparam.h"

!	==========================================================
!
!	SULFUT.F:  Subroutine for setup prescribed sulfate
!			aerosol forcing in the format of
!			surface albedo.
!
!	----------------------------------------------------------
! 
!	Revision History:
!	
!	When	Who		What
!	-----	----------	-------	
!	080100	Chien Wang	repack based on CliChem3 & M24x11,
!				  and add cpp.
!	081100	Chien/Andrei	open cpp control for all case.
!
!	==========================================================

      subroutine sulfr(BSO4LAND,BSO4OCEAN,TNOW)

#include "BD2G04.COM"

C     in original sulf1986.4x5.1986.dat file BSO4LAND and BSO4OCEAN
c     are SO2 loading per m^2 over land and ocean 
C     in  sulf1986.4x5.1986.new.dat file BSO4LAND and BSO4OCEAN
c     are SO2 loading per m^2 over land and ocean multiplied by 
c     land and ocean fraction respectively
      dimension BSO4LAND(JM0),BSO4OCEAN(JM0),
     & BSO4LD86(JM0),BSO4ON86(JM0),
     & IYSULF(2000),CF86(2000),YSULF(2000)
      character * 256 sulf1986,sulf2050,sulfamp,SO2_EM
      character * 80 title
      common /sulfdata/sulf1986,sulf2050,sulfamp,SO2_EM
       logical first
       data first /.true./

       if(first) then
       JM=JM0
       open (575,file=sulf1986,status='old')
       open (577,file=SO2_EM,status='old')
       read (575,'(E13.5)'),(BSO4LD86(J),j=1,JM)
       read (575,'(E13.5)'),(BSO4ON86(J),j=1,JM)
       read (577,'(A80)'),title
       print *,'From sulfut'
       print *,title
       read (577,'(A80)'),title
       print *,title
       read (577,'(A80)'),title
       print *,title
       do ny=1,2000
        read (577,*,end=500),IYSULF(ny),CF86(ny)
        if(IYSULF(ny).eq.1986)cfref=CF86(ny)
!       print *,IYSULF(ny),CF86(ny)
       enddo
 500   continue
       NYS=ny-1
       close (575)
       close (577)
       print *,'NYS=',NYS
       print *,'Sulfate data for years ',IYSULF(1),IYSULF(NYS)
       print *,' BSO4LD86'
       print ('(6F7.4)'),(BSO4LD86(J)*1.e6,j=1,JM)
       print *,' BSO4ON86'
       print ('(6F7.4)'),(BSO4ON86(J)*1.e6,j=1,JM)
       do n=1,NYS
        YSULF(n)=IYSULF(n)+0.5
        CF86(n)=CF86(n)/cfref
       enddo
       first=.false.
       end if
       do n=1,NYS-1
        if(TNOW.gt.YSULF(n).and.TNOW.le.YSULF(n+1))go to 100
       enddo
       print *,' Wrong TNOW TNOW=',TNOW
       stop
  100  continue
        x=(YSULF(n+1)-TNOW)/(YSULF(n+1)-YSULF(n))
        cf=x*CF86(n)+(1.-x)*CF86(n+1)
        do j=1,jm
         BSO4LAND(J)=cf*BSO4LD86(j)
         BSO4OCEAN(J)=cf*BSO4ON86(j)
        enddo
!       print *,' TNOW=',TNOW
!       print *,' n=',n,' YSULF(n)=',YSULF(n),
!    &   ' YSULF(n+1)=',YSULF(n+1)
!      print *,' BSO4LAND'
!      print ('(6F7.4)'),(BSO4LAND(J)*1.e6,j=1,JM)
!      print *,' BSO4OCEAN'
!      print ('(6F7.4)'),(BSO4OCEAN(J)*1.e6,j=1,JM)
       return
       end

      subroutine bchist(bcodland,bcodocean,iybc,nmonth,ndata)
#if ( !defined CPL_CHEM )
#if ( defined HISTORICAL_BC)

#include "BD2G04.COM"

       real *4 bclandm(JM0,LM0,0:13,0:150),
     &   bcoceanm(JM0,LM0,0:13,0:150)
     &   ,bctotalm(JM0,LM0,0:13,0:150)
      real *4 bcarlann(JM,LM,0:13),bcarwann(JM,LM,0:13),
     & bcartann(JM,LM,0:13)
       dimension bcland(JM0,LM0),bcocean(JM0,LM0),bctotal(JM0,LM0)
       dimension bcodland(JM0,LM0),bcodocean(JM0,LM0),bcodtotal(JM0,LM0)
       dimension NDAYMN(12),iyears(0:150)
      character * 256 bc_land,bc_ocean, bc_total
      character * 80 title
       common /bc_datah/bc_land,bc_ocean, bc_total
       logical first
       data NDAYMN /31,28,31,30,31,30,31,31,30,31,30,31/
       data first /.true./

       if(first) then
       JM=JM0
       print *,'Data for black carbon'
       print *,bc_land
!      qex_bc = 8.0e-6/dxyp(j) ! normal
       qex_bc = 8.0e-6 ! normal
       qex_bc = 8.0e-6*1.e9 ! normal
       bcland0 =0.0
       bcocean0=0.0
       iyear0=1860
       open(769,file=bc_land,
     &  status='old',form='unformatted')
       open(768,file=bc_ocean,
     &  status='old',form='unformatted')
!      read(769),iyears
!      read(769),bclandm
!      read(768),iyears
!      read(768),bcoceanm
      read (769),nyt
      read (768),nyt
!     print *,'nyt=',nyt
      do ni=0,nyt
        read (769),iyears(ni)
        read (769),bcarlann
        read (768),iyears(ni)
        read (768),bcarwann
       do imonth=0,13
        do j=1,jm
         do l=1,lm
          bclandm(j,l,imonth,ni)=bcarlann(j,l,imonth)
          bcoceanm(j,l,imonth,ni)=bcarwann(j,l,imonth)
         enddo  ! j
        enddo   ! l
       enddo    ! imonth
       enddo    ! ni
!      print *,iyears
        close (769)
        close (768)
       print *, 'from ',iyears(0),' to ', iyears(150)
        first=.false.
      endif
!     print *,'From bcload ', iybc,nmonth,ndata
      iyr=iybc-1860
      if(ndata.le.16) then
       x=(16-ndata)/31.
       do l=1,LM0
       do j=1,JM0
         bcland(j,l)=bclandm(J,L,nmonth-1,iyr)*x
     &   +bclandm(J,L,nmonth,iyr)*(1-x)
         bcocean(j,l)=bcoceanm(J,L,nmonth-1,iyr)*x
     &   +bcoceanm(J,L,nmonth,iyr)*(1-x)
       enddo
       enddo
      else
       x=1-(ndata-16)/31.
       do l=1,LM0
       do j=1,JM0
         bcland(j,l)=bclandm(J,L,nmonth,iyr)*x
     &   +bclandm(J,L,nmonth+1,iyr)*(1-x)
         bcocean(j,l)=bcoceanm(J,L,nmonth,iyr)*x
     &   +bcoceanm(J,L,nmonth+1,iyr)*(1-x)
       enddo
       enddo
      endif
       do j=1,JM0
        bcodland(j,LM)=bcland(j,LM)*qex_bc
        bcodocean(j,LM)=bcocean(j,LM)*qex_bc
       enddo
       do l=LM0-1,1,-1
       do j=1,JM0
        bcodland(j,L)=bcodland(j,L+1)+bcland(j,L)*qex_bc
        bcodocean(j,L)=bcodocean(j,L+1)+bcocean(j,L)*qex_bc
       enddo
       enddo
!      if (ndata.eq.15) then
!     print *,'From bcload ', iybc,nmonth,ndata
!      print *,' BCOD'
!      do l=1,LM0
!       print *,l,bcodland(23,L),bcodocean(23,L)
!      enddo
!     endif
#endif
#endif
       return
       end
