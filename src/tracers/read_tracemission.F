#include "ctrparam.h"

        subroutine read_tracemission (coef)
!       ===================================
#include "chem_para"
#include "BD2G04.COM"
#include "tracer_com.h"

        logical first
        data first /.true./
        save totaltremis
        if(first) then
          totaltremis=0.0
          first=.false.
        endif
         
        tracemis = 0.0
        totaltremis =0.0

        if (JYEAR.eq.IYREMIS) then
!        print *,'Using TREMIS for ',JYEAR,JDAY
         i=1
         do n=1,ntracers
           do j=1,nlat
            do k=1,nlev_accri
             tracemis (i,j,k,n)=trac_emissions(j,k,JDAY)
             if(n.eq.1) then
             totaltremis=totaltremis+trac_emissions(j,k,JDAY)
             endif
            enddo
          enddo
         enddo
          tracemis = coef*tracemis
!        print *,'TRAC EMIS'
!        do k=1,nlev_accri
!            print *,tracemis (1,39,k,1)
!         enddo
         print *,'TRAC EMIS=',JDAY,totaltremis

        endif
	return
	end
