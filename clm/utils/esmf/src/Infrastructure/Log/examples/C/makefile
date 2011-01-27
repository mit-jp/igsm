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
LOCDIR	  = src/LowLevelUtils/Log/examples/C

TESTEXAMPLES_1  = \
	ESMC_LogEx.Alice RUN_ESMC_LogEx ESMC_LogEx.rm

TESTEXAMPLES_4  = \
	ESMC_LogEx.Alice RUN_ESMC_LogExUNI ESMC_LogEx.rm

EXAMPLESC       = ESMC_LogEx.F

include ${ESMF_DIR}/build/${ESMF_ARCH}/base

MFMLOGDIR = ${ESMC_TESTDIR}

ESMC_LogEx: ESMC_LogEx.o  chkopts
	-${CLINKER} -o ESMC_LogEx ESMC_LogEx.o ${ESMC_TIME_LIB} -lesmf ${MPI_LIB} ${MP_LIB} ${THREAD_LIB} ${PCL_LIB}
	${RM} -f ESMC_LogEx.o

RUN_ESMC_LogEx:
	-@rm -f testlog stdlog.* logfile.* ESMC_LogEx.tmps
	${MPIRUN} -np 4 ESMC_LogEx > /dev/null
	@sort -t, -k2 testlog | cut -d, -f4- > ESMC_LogEx.tmp
	@if (diff ESMC_LogEx.tmp output/ESMC_LogEx.out) then \
	echo "PASS Logger test" > ${MFMLOGDIR}/ESMC_LogEx.log; else \
	echo "FAIL Logger test" > ${MFMLOGDIR}/ESMC_LogEx.log; fi
	@rm -f testlog stdlog.* logfile.* ESMC_LogEx.tmp optlog

RUN_ESMC_LogExUNI:
	-@rm -f testlog stdlog.* logfile.* ESMC_LogEx.tmps
	${MPIRUN} -np 1 ESMC_LogEx > /dev/null
	@sort -t, -k2 testlog | cut -d, -f4- > ESMC_LogEx.tmp
	@if (diff ESMC_LogEx.tmp output/ESMC_LogEx.outuni) then \
	echo "PASS Logger test" > ${MFMLOGDIR}/ESMC_LogEx.log; else \
	echo "FAIL Logger test" > ${MFMLOGDIR}/ESMC_LogEx.log; fi
	@rm -f testlog stdlog.* logfile.* ESMC_LogEx.tmp optlog

