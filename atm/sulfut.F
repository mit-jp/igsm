
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
      character * 120 sulf1986,sulf2050,sulfamp,SO2_EM
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
c       print *,' TNOW=',TNOW
c       print *,' n=',n,' YSULF(n)=',YSULF(n),
c    &   ' YSULF(n+1)=',YSULF(n+1)
c      print *,' BSO4LAND'
c      print ('(6F7.4)'),(BSO4LAND(J)*1.e6,j=1,JM)
c      print *,' BSO4OCEAN'
c      print ('(6F7.4)'),(BSO4OCEAN(J)*1.e6,j=1,JM)
       return
       end
