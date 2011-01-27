! $Id$
	!\begin{verbatim}
	program ProfilerEx
	Use ESMF_TimerMod
	Use ESMF_AppMod
	Use ESMF_MachineMod
	implicit none
	integer :: iter
	integer :: rc
	external do_work
	integer :: ierr
	integer :: node, process, thread
	type(ESMF_App) :: app
	type (ESMF_Log) :: log
	type (ESMF_Log) :: tlog
	character(60) :: str
	logical :: test

#include "ESMF_conf.h"

	app = ESMF_AppNew()

	log = ESMF_LogNew("logfile", ESMF_LOGSTATE_VERBOSE, 1)
	tlog = ESMF_LogNew("timelog", ESMF_LOGSTATE_VERBOSE, 1)


	call ESMF_LogPrint(log, ESMF_LOGLEVEL_INFO,
     &	"Hello from ESMF_ProfilerEx")

	
	call ESMF_TimerInit(ESMF_ALL_TIMERS, ESMF_USRSYS, 1, 
     &	ESMF_WALL, 1, 0)

	call ESMF_TimerStart("Test1", rc)
	test = (rc .eq. ESMF_SUCCESS)
	str = "ESMF_TimerInit: Return Ok"
	call ESMF_ErrorTest(test, str)

!$OMP PARALLEL DO PRIVATE(iter,node,process,thread)
	do iter=1,10
	   call do_work(log, iter)
	   call ESMF_MachinePInfo(node, process, thread)
	   write(*,FMT='(a,I4, a, I4, a, I4, a, I4)') "node=", node,
     &	"proc=", process, "thread=", thread, "iter:", iter
	end do
!$OMP END PARALLEL DO

	call ESMF_TimerStop("Test1", rc)
	test = (rc .eq. ESMF_SUCCESS)
	str = "ESMF_TimerInit: Return Ok"
	call ESMF_ErrorTest(test, str)

! Try default print
	call ESMF_TimerPrint(ESMF_ALL_TIMERS)

	call ESMF_TimerPrint(ESMF_ALL_TIMERS, tlog)


	call ESMF_AppDelete(app)

	end program

	subroutine do_work(log, iter)
	Use ESMF_TimerMod
	Use ESMF_LogMod
	implicit none
	type(ESMF_Log) :: log
	integer :: iter
	integer :: i, j, k
	real(8) :: x, y
	character*7 name
	
	write(name, '(a4,i3.3)') 'iter', iter

	call ESMF_TimerStart(name)
	print *, "Started timer:", name
	x = 0.25
	y = 0.25
	do i = 1, 5*1024
	   x = x * x + y
	enddo
	
	call ESMF_LogPrint(log, ESMF_LOGLEVEL_INFO,
     &	"In iteration:%d\n", iter)
	print *, "ESMF_LOGLEVEL_INFO = ", ESMF_LOGLEVEL_INFO
	write(*,*) 'In iter:', name

	call ESMF_TimerStop(name)
	print *, "Stopped timer:", name


	
	end subroutine !\end{verbatim}
