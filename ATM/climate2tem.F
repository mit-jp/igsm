C $Header: /home/igsm-repository/ATM/climate2tem.F,v 1.1 2010/08/30 15:36:35 jscott Exp $
C $Name:  $

#include "ctrparam.h"
! 	This subroutine is called one a month by driver
!	to normalize monthly climate data for TEM
	subroutine climate2tem(month,ndays)
#include "BD2G04.COM"
#include "TEM.h"        
        print *,'climate4tem'
	do j=1,jm0
c         co24tem(j)=co24tem(j)/ndays
	  temp4tem(j)=temp4tem(j)/ndays
	  sws4tem(j)=sws4tem(j)/ndays
	  pre4tem(j)=pre4tem(j)/ndays
#if ( !defined CPL_CHEM )
	  do n=1,no3
	   o34tem(n,j)=obso3(n,j,month)
	  enddo
#endif
	enddo
        print *,'Month=',month,' CO2(23) for TEM=',co24tem(jm0/2)
c       print *,'O3 for TEM J=23'
c       print *,(o34tem(n,23),n=1,no3)
#if ( defined DATA4TEM )
        write(937),month,ndays,co24tem,temp4tem,pre4tem,sws4tem,
     &		o34tem,dtem4tem
c       write(937),month,ndays,co24tem,
c    &            o34tem,
c    &            temp4tem,
c    &            dtem4tem,
c    &            sws4tem,
c    &            pre4tem,
c    &            strmdur,
c    &            qstrm,
c    &            aet,
c    &            sh2o1m,
c    &            sh2o2m,
c    &            daytsoil,
c    &            daysho,
c    &            rsh2o
#endif
	return
	end
