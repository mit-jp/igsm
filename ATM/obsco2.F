C $Header: /home/igsm-repository/ATM/obsco2.F,v 1.1 2010/08/30 15:36:35 jscott Exp $
C $Name:  $

#include "ctrparam.h"

      subroutine  obsco2(iyearr,nm)
#include "BD2G04.COM" 
#include "chem_para"
      parameter (nyears=337)
      character * 120 oco2file,co2rfile
      dimension co2rat1990(nlat,12),co2rat1860(nlat,12),
     * co2(nlat,12)
      dimension iyearobs(nyears),anco2obs(nyears)
      common /ATCO2/atm_co2(nlat),oco2file,co2rfile
      data ifirst /1/
      if(ifirst.eq.1) then
c       open( unit=645,file=oco2file,
c    &  status='old',form='formatted')
        open( unit=561,file=oco2file,
     &  status='old',form='formatted')
        open( unit=655,file=co2rfile,
     &  status='old',form='formatted')
        do n=1,nyears
c         read(645,'(i4,f8.2)') iyearobs(n),anco2obs(n)
          read(561,end=500,*) iyearobs(n),anco2obs(n)
        enddo
  500   continue
        close (561)
        print *,'OBSCO2 from ',iyearobs(1),' to ',iyearobs(n-1)
        read(655,'(i5,f7.2)'),iiyear,cco2
         do j=1,nlat
          read(655,'(i5,12f7.2)'),ii,(co2rat1990(j,m),m=1,12)
          do m=1,12
           co2rat1860(j,m)=1.0
          enddo
         enddo
         dy=1991.-1860.
         ifirst=0
      print *,'JM from obsco2=',jm
      print *,iiyear,cco2
      endif
C    End of initialization
c
      nyear=iyearr-iyearobs(1)+1
      co2an=anco2obs(nyear)
c     print *,'From obsco2 year=',iyearr,' nm=',nm,' co2an=',co2an
      print *,'From obsco2 year,nm,co2an ',iyearr,nm,co2an

      if(iyearr.gt.1860)then
       x=(1991.-iyearr)/dy
         do 30 j=1,nlat
          atm_co2(j)=co2an*(x*co2rat1860(j,nm)+(1.-x)*co2rat1990(j,nm))
   30    continue
       else
         do  j=1,nlat
          atm_co2(j)=anco2obs(nyear)
         enddo
       endif
       if(nm.eq.1)then
        print '12f7.1,/,2(11f7.1,/,),12f7.1',atm_co2
       endif
      return
      end
