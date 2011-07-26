
#include "ctrparam.h"

!	==========================================================
!
!	FORCEDOZONE.F:  Routine to get ozone deviations from 
!			specified file and interpolate
!			to the appropriate data
!
!	----------------------------------------------------------
! 
!	Revision History:
!	
!	When	Who		What
!	----	----------	-------	
!	073100	Chien Wang	repack based on CliChem3 & M24x11,
!				  and add cpp.
!
!       110100                 Added data for Jan. 2000, 
!                              dimension of o3f was changed from nm to nm+1
!
!       011503  Chris Forest   Converted to use GISS 1850-2050 ozone data
!                              uses nlev = 14 to adjust ozone profile in three 
!                              layers above 10 hPa top of dynamics levels
!                              (layer edges at 5, 2, and 0.00001 hPa)
!                                 
!	==========================================================

      subroutine getforcedozone( jyear, jday )

      parameter (nm=1728)
c	parameter (nm =252)

#include "chem_para"
#include "chem_com"

      common /o3data/o3_data
      character * 120 o3_data
C      
c     Routine to get ozone deviations from specified file and interpolate
c       to the appropriate date.
C     
C     Author: Chris E Forest. 
C     Added Mar 18, 1998
C
C      tropmass = 28.97296245*1.e-3*0.8/P0
C      pxxx = dpl(l)*part
C      
C      ULGAS(L,3)=o3(ILON,jyyy,l)/48.0
C     &     *pxxx*tropmass
C      
C      o3 = ppb(m)
C      
C      48 = mol weight of o3
C      
C      ULGAS = cm^3 (STP)/cm^2

      logical first
      integer year0, nlevp3 ,yrmax
!     real date, o3f(nlat,nlev+3,nm), o3dummy ! 03/20/2006
      real date, o3f(nlat,nlev+3,nm+12), o3dummy
c      real date, o3f(24,11,nm+1), o3dummy
      character*43 filen
      DIMENSION JDOFM(13)
      data first / .true. /
      DATA JDOFM /0,31,59,90,120,151,181,212,243,273,304,334,365/
      data year0 / 1859 /

c   year0 = 1974 for o3hadley.dat
c   year0 = 1979 for o3giss.MR.dat
c   year0 = 1850 for o3.1850_2050.giss.46x11.dat
c   year0 = 1859 for o3.1859_2002.giss.46x11.dat
c
C      filen = '/data11/ceforest/forcing/dummyozone.dat'
c     filen = '/data11/ceforest/forcing/o3giss.MR.dat'
C      filen = '/data11/ceforest/forcing/o3hadley.dat'

      if (first) then
         print *,'GetForcedOzone: reading data for 14 layers'
          print *,'nlat=',nlat,' nlev=',nlev
         open(unit=500, file=o3_data,status='old')
c
c            GISS: nm = 252 , HadCM2: nm=265, 1850-2050: nm=2412
c
c        nlev = nlev + 3 to include 3 layers above 10 hPa.
c

         nlevp3 = nlev + 3

         do i=1,nm
            do k=1, nlevp3
               do j = 1,nlat
                  read(500,'(E12.4)') o3dummy
                  o3f(j,k,i) = o3dummy*1.e9       ! convert to ppb(m)
               enddo
            enddo
         enddo
         close(500)
C     Added data for Jan. 2000 to allow simulation until Jan 1 2000
c     NOT necessary with 1850-2050 GISS ozone data
c            do k=1, nlev
c               do j = 1,nlat
c                  o3f(j,k,nm+1) = o3f(j,k,nm)
c               enddo
c            enddo
! 03/20/2006 Data for anaditional year are added
!            really only data for January are used
           do m=1, 12
             do k=1, nlev
                do j = 1,nlat
                   o3f(j,k,nm+m) = o3f(j,k,nm-12+m)
                enddo
             enddo
           enddo
         first=.false.
      nyears=nm/12
      yrmax=year0+nyears-1
      print *,'Ozone data for ',nyears,' years '
      print *,'from ',year0,' to ',yrmax
      endif

      
!     iy = jyear - year0 ! 03/17/06
      jyear1=min(jyear,yrmax)
      iy = jyear1 - year0 
!     print *,'Ozone for year ',jyear1,iy
      if (iy.lt.0) then
         do i=1,nlat
            do j=1,nlevp3
               o3dev(1,i,j) = 0.
            enddo
         enddo
      else
C---  o3interp(o3f, o3dev, date, year0)
       
         jdd = jday-15

         if (iy .eq. 0 .and. jdd.le.0) then
            do i=1,nlat
               do j=1,nlevp3
                  dl = ( 31. + float(jdd) )/31.
                  o3dev(1,i,j) = o3f(i,j,1)*dl
               enddo
            enddo
         else
            im =12
            do i=1,12 
               if (jdd.gt.jdofm(i) .and. jdd.le.jdofm(i+1)) then
                  im = i
               endif
            enddo
            if (iy .ge. 1 .and. im.eq.12 .and. jdd.le.0 )  iy = iy-1 
C               if (jdd.le.0) iy = iy-1
C            endif
           
            dl = ( float(jdd) - jdofm(im))/(jdofm(im+1) - jdofm(im)) 
            if (im.eq.12.and.jdd.le.0) dl = ( 31. + float(jdd) )/31.
               
            imm = im + 12*iy
!            print *,jyear1,jday,iy,im,imm
            
            if(imm.lt.1.or.imm.ge.nm+1)then
               print *,' error in o3interp im=',imm
               stop 25
            endif
            
            do j=1,nlat
               do k=1,nlevp3
                  o3dev(1,j,k) = o3f(j,k,imm)*(1.-dl)+o3f(j,k,imm+1)*dl
               enddo
            enddo
C            print *,jyear,jday,iy,im,dl,o3dev(0,0,9)
         endif 
      endif
      
c     print *,'GetForcedOzone: JYEAR, JDAY, IY=', jyear,jday, iy 
c     if (iy.ge.0) print *,'GetForcedOzone: imm=', imm

      return
      end


