# $Id$ 

ALL: 

CFLAGS    =
FFLAGS    =
SOURCEC	  = 
SOURCEF	  = 

SOURCEH	  = 
OBJSC	  =  
OBJSF	  = 
MANSEC	  = 
LIBBASE	  = 
DIRS	  = 
LOCDIR	  = src/Infrastructure/Timer/examples/F90

TESTEXAMPLES_3  = \
	ESMF_TimerEx RUN_ESMF_TimerEx ESMF_TimerEx.rm

TESTEXAMPLES_9  = \
	ESMF_TimerEx RUN_ESMF_TimerExUni ESMF_TimerEx.rm

EXAMPLESF       = ESMF_TimerEx.F

include ${ESMF_DIR}/build/${ESMF_ARCH}/base

ESMF_LogDir:
	mkdir -p ${ESMC_TESTDIR}

ESMF_TimerEx: ESMF_TimerEx.o  chkopts
	-${FLINKER} -o ESMF_TimerEx ESMF_TimerEx.o ${ESMC_TIME_LIB} -lesmf ${MPI_LIB} ${MP_LIB} ${THREAD_LIB} 
	${RM} -f ESMF_TimerEx.o 

RUN_ESMF_TimerEx:
	${MPIRUN} -np 4 ESMF_TimerEx > ${ESMC_TESTDIR}/ESMF_TimerEx.log
	${RM} -f logfile.* stdlog.* timelog.*

RUN_ESMF_TimerExUni:
	${MPIRUN} -np 1 ESMF_TimerEx > ${ESMC_TESTDIR}/ESMF_TimerEx.log
	${RM} -f logfile.* stdlog.* timelog.*












