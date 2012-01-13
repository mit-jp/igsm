
#include "ctrparam.h"

!	============================================================
!
!	TRACNIT.F:	Subroutine for initialization of tracers
!				in MIT Global Chemistry Model
!
!	------------------------------------------------------------
!
!
	subroutine tracinit(ISTRT1)
!	=================================

#include "chem_para"
#include "chem_com"
#include "tracer_com.h"
       data trtype /1,2,0,3,4/
       data treftime /5.0,21.0,0.0,0.0,0.0/
c----------------------------------------------
c   Set  navtrac = 0 and start monthly average
c       calculation:
c
      character * 120 trac_emiss_data,trac_monthly,trac_global
      namelist /TRSERS/ IYREMIS,trac_emiss_data
     &   ,trac_monthly,trac_global

        open(535,file='tracemis.dat')
       print *,' CALL Traser Emissions'
        print *,' nlev_accri=',nlev_accri
        read (535,NML=TRSERS)
        close (535)
        open(664,file=trac_emiss_data,
     &        form='unformatted',
     &        status='old')


        read(664)p_accri
        read(664)trac_emissions
        close (664)
        print *,trac_emiss_data
        print *,'p_accri'
        print *,p_accri

       navtrac = 0
       monthtrac =1 

       tracmonth    = 0.0

        open(269,file='trac_monthly',status='new',
     &        form='unformatted')
        open(276,file='trac_global',status='new',
     &   form='formatted')
c
	if(ISTRT1.eq.0)then
         tracers = 0.0
	endif

c ===== To renew a run:
c
	if(ISTRT1.eq.1)then
 	 read(279)tracers
	 close(279)
	endif
       return
        end

