
#include "ctrparam.h"

      subroutine  obsco2(iyearr,nm)
#include "BD2G04.COM" 
#include "chem_para"
      parameter (nyears=2000)
      character * 256 oco2file,co2rfile
      dimension co2rat1990(nlat,12),co2rat1860(nlat,12),
     * co2rat2005(nlat,12),co2(nlat,12)
      dimension iyearobs(nyears),anco2obs(nyears)
      common /ATCO2/atm_co2(nlat),oco2file,co2rfile
      data ifirst /1/
      if(ifirst.eq.1) then
        print *,'From obsco2 nyears=', nyears
        print *,oco2file
c       open( unit=645,file=oco2file,
c    &  status='old',form='formatted')
        open( unit=561,file=oco2file,
     &  status='old',form='formatted')
        open( unit=655,file=co2rfile,
     &  status='old',form='formatted')
        do n=1,nyears
c         read(645,'(i4,f8.2)') iyearobs(n),anco2obs(n)
          read(561,end=500,*) iyearobs(n),anco2obs(n)
!         print *, iyearobs(n),anco2obs(n)
        enddo
  500   continue
        close (561)
        print *,'OBSCO2 from ',iyearobs(1),' to ',iyearobs(n-1)
         do j=1,nlat
          do m=1,12
           co2rat1860(j,m)=1.0
          enddo
         enddo
         read(655,'(i5,f7.2)'),iiyear1,cco21
         do j=1,nlat
          read(655,'(i5,12f7.2)'),ii,(co2rat1990(j,m),m=1,12)
         enddo
         print *,iiyear1,cco21
         read(655,'(i5,f7.2)'),iiyear2,cco22
         do j=1,nlat
          read(655,'(i5,12f7.2)'),ii,(co2rat2005(j,m),m=1,12)
         enddo
         print *,iiyear2,cco22
!        dy=1991.-1860.
         ifirst=0
      endif
C    End of initialization
c
      nyear=iyearr-iyearobs(1)+1
      co2an=anco2obs(nyear)
c     print *,'From obsco2 year=',iyearr,' nm=',nm,' co2an=',co2an
      print *,'From obsco2 year,nm,co2an ',iyearr,nm,co2an

      if(iyearr.le.1860)then
         do  j=1,nlat
          atm_co2(j)=anco2obs(nyear)
         enddo
      else
       if(iyearr.le.iiyear1)then
         dy=iiyear1+1-1860.
         x=(iiyear1+1-iyearr)/dy
         do  j=1,nlat
          atm_co2(j)=co2an*(x*co2rat1860(j,nm)+(1.-x)*co2rat1990(j,nm))
         enddo
       else
        if(iyearr.le.iiyear2)then
         dy=iiyear2+1-iyear1
         x=(iiyear2+1-iyearr)/dy
         do  j=1,nlat
          atm_co2(j)=co2an*(x*co2rat1990(j,nm)+(1.-x)*co2rat2005(j,nm))
         enddo
        else
         do  j=1,nlat
          atm_co2(j)=co2an*co2rat2005(j,nm)
         enddo
        endif
       endif
      endif
       if(nm.eq.1)then
        print *,'ATM_CO2 from OBSCO2'
        print '(12f7.1,/,2(11f7.1,/,),12f7.1)',atm_co2
       endif
      return
      end
