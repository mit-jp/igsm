! $Id$
	!\begin{verbatim}
	program LogEx
#include "ESMF_conf.h"
	Use ESMF_AppMod
	Use ESMF_LogMod
	implicit none


	integer :: ivar
	integer :: i
	integer :: ierr
	real(8) :: dvar
	real(4) :: fvar
	type(ESMF_Log) :: log
	type(ESMF_Log) :: testlog
        type(ESMF_Log) :: opttest
	type(ESMF_App) :: app

	app = ESMF_AppNew()

	log = ESMF_LogNew("logfile", ESMF_LOGSTATE_VERBOSE, 1)
	testlog = ESMF_LogNew("testlog",
     &	ESMF_LOGSTATE_VERBOSE, 0)
        opttest = ESMF_LogNew(
     &	"optlog", ESMF_LOGSTATE_VERBOSE, 0)

! AAA guarantes the sort comes out right
	call ESMF_LogPrint(testlog, ESMF_LOGLEVEL_INFO,
     &	"AAA: Hello Logfile")



	ivar = 2
	dvar = 3.1459
	fvar = 0.145987

	call ESMF_LogPrint(testlog, ESMF_LOGLEVEL_ERROR,
     &	"Assorted Vars:%s%08d,dvar=%10.2lf,fvar=%10.4f\n",
     &	"i=", ivar, dvar, fvar)

	call ESMF_LogFlush(testlog)
	
	call ESMF_LogPrint(log, ESMF_LOGLEVEL_INFO,
     &	"Hello World Logger\n")
	call ESMF_LogPrint(log, ESMF_LOGLEVEL_INFO, 
     &  "%s %08d, dvar=%lf, fvar=%f\n",
     &	"i=", i, dvar, fvar)

!$OMP PARALLEL DO PRIVATE (i)
	do i=1, 10
		call ESMF_LogPrint(log, ESMF_LOGLEVEL_INFO,
     &	"Hello from do loop, iter= %d\n", i)

	call ESMF_LogPrint(testlog, ESMF_LOGLEVEL_INFO, 
     &	"Called from OMP Parallel do\n")
	enddo
!$OMP END PARALLEL DO

! Now test the option setting
        call ESMF_LogPrint(
     &	opttest,
     &  ESMF_LOGLEVEL_INFO,
     &	 "This should display defaults\n")

        call ESMF_LogSetConfig(opttest,
     &	ESMF_LOGCONFIG_PRINTTIME, 0,
     &	ESMF_NULL)
        call ESMF_LogPrint(
     &	opttest,
     &  ESMF_LOGLEVEL_INFO,
     &	 "This should not display time\n")

        call ESMF_LogSetConfig(opttest,
     &	ESMF_LOGCONFIG_PRINTPID, 0,
     &	ESMF_NULL)
        call ESMF_LogPrint(
     &	opttest,
     &  ESMF_LOGLEVEL_INFO,
     &	 "This should not display time or pid\n")

        call ESMF_LogSetConfig(opttest,
     &	ESMF_LOGCONFIG_PRINTTID, 0,
     &	ESMF_NULL)
        call ESMF_LogPrint(
     &	opttest,
     &  ESMF_LOGLEVEL_INFO,
     &	 "This should not display time or pid or tid\n")

        call ESMF_LogSetConfig(opttest,
     &	ESMF_LOGCONFIG_PRINTTIME, 1,
     &	ESMF_LOGCONFIG_PRINTPID, 1,
     &	ESMF_LOGCONFIG_PRINTTID, 1,
     &	ESMF_NULL)
        call ESMF_LogPrint(opttest, ESMF_LOGLEVEL_INFO, 
     &	"This should display defaults\n")
	call ESMF_LogDelete(testlog)
	call ESMF_LogDelete(log)
        call ESMF_LogDelete(opttest)

	call ESMF_AppDelete(app)



	end program !\end{verbatim}
