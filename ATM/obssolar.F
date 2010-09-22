
#include "ctrparam.h"

!	==========================================================
!
!	OBSSOLAR.F:  This file contains subroitine obssolar.f
!		 which reads observed solar constant
!
!	----------------------------------------------------------
! 
!	Revision History:
!	
!	When	Who		What
!	----	----------	-------	
!	070302	
!	
!
!	==========================================================

      subroutine obssolar(S0const,YEAR1)
      parameter (nyd=2000,nyd1=nyd+1)
C
C
C     For simulations with observed solar constant using equivalent
C     S0conts  from file S0C_data
C     DATA are annual means, but they are
C       used as a data for the middle of the coresponding year
C
C
      dimension  solar(nyd),solar1(nyd1),yearS0(nyd)
      common /solardata/S0C_data
      character * 120 S0C_data
      logical first
      data first /.true./
      if(first)then
         print *,'  Solar constant from observation '
         open (unit=561,file=S0C_data,
     *     status='OLD',form='formatted')
          do 50 i=1,nyd
c        read (561,501,end=500),ighgyr(i),tmp,xco2(i)
            read (561,*,end=500),yearS0(i),solar(i)
c           print *,yearS0(i),solar(i)
   50    continue
  500    continue
         close (561)
         ny=i-1
         print *,'Last year for solardata is ',yearS0(ny)
         ny1=ny+1
         ny2=ny+2
  501    format(2f13.2)
         first=.false.
      endif
      year=year1
c        print *,' year=',year
      if(year.lt.yearS0(1))year=yearS0(1)
      if(year.gt.yearS0(ny))year=yearS0(ny)
      i=year-yearS0(1)+1
      dl=year-yearS0(i)
c     print *,' i=',i,' dl=',dl
c     print *,yearS0(i),yearS0(i+1)
c     print *,solar(i),solar(i+1)
      S0const=solar(i)*(1.-dl)+solar(i+1)*dl
c          print *,' year=',year,' S0x=',S0const
      return
      end
