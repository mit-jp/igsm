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
LOCDIR	  = src/LowLevelUtils/Timer/examples/C

TESTEXAMPLES_1  = \
	ESMC_TimerEx.Alice RUN_ESMC_TimerEx ESMC_TimerEx.rm

TESTEXAMPLES_4  = \
	ESMC_TimerEx.Alice RUN_ESMC_TimerExUNI ESMC_TimerEx.rm

EXAMPLESC       = ESMC_TimerEx.F

include ${ESMF_DIR}/build/${ESMF_ARCH}/base

MFMLOGDIR = ${ESMC_TESTDIR}

ESMC_TimerEx: ESMC_TimerEx.o  chkopts
	-${CLINKER} -o ESMC_TimerEx ESMC_TimerEx.o ${ESMC_TIME_LIB} -lesmf ${MPI_LIB} ${MP_LIB} ${THREAD_LIB}
	${RM} -f ESMC_TimerEx.o

RUN_ESMC_TimerEx:
	-@rm -f testlog stdlog.* logfile.* ESMC_TimerEx.tmps
	${MPIRUN} -np 4 ESMC_TimerEx > ${MFMLOGDIR}/ESMC_TimerEx.log
	@rm -f testlog stdlog.* logfile.* 

RUN_ESMC_TimerExUNI:
	-@rm -f testlog stdlog.* logfile.* ESMC_TimerEx.tmps
	${MPIRUN} -np 1 ESMC_TimerEx > ${MFMLOGDIR}/ESMC_TimerEx.log
	@rm -f testlog stdlog.* logfile.* 

