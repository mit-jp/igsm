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
LOCDIR	  = src/Infrastructure/Log/examples/F90

TESTEXAMPLES_3  = \
	ESMF_LogEx RUN_ESMF_LogEx ESMF_LogEx.rm

TESTEXAMPLES_9  = \
	ESMF_LogEx RUN_ESMF_LogExUNI ESMF_LogEx.rm

EXAMPLESF       = ESMF_LogEx.F

include ${ESMF_DIR}/build/${ESMF_ARCH}/base

ESMF_LogDir:
	mkdir -p ${ESMC_TESTDIR}

ESMF_LogEx: ESMF_LogEx.o  chkopts
	-${FLINKER} -o ESMF_LogEx ESMF_LogEx.o ${ESMC_TIME_LIB} -lesmf ${MPI_LIB} ${MP_LIB} ${THREAD_LIB} ${PCL_LIB}
	${RM} -f ESMF_LogEx.o

RUN_ESMF_LogEx:
	-@rm -f testlog stdlog.* logfile.* ESMF_LogEx.tmps
	${MPIRUN} -np 4 ESMF_LogEx > /dev/null
	@sort -t, -k2 testlog | cut -d, -f4- > ESMF_LogEx.tmp
	@if (diff ESMF_LogEx.tmp ${ESMF_DIR}/${LOCDIR}/output/ESMF_LogEx.out) then \
	echo "PASS Logger test" > ${ESMC_TESTDIR}/ESMF_LogEx.log; else \
	echo "FAIL Logger test" > ${ESMC_TESTDIR}/ESMF_LogEx.log; fi
	@rm -f testlog stdlog.* logfile.* ESMF_LogEx.tmp optlog

RUN_ESMF_LogExUNI:
	-@rm -f testlog stdlog.* logfile.* ESMF_LogEx.tmps
	${MPIRUN} -np 1 ESMF_LogEx > /dev/null
	@sort -t, -k2 testlog | cut -d, -f4- > ESMF_LogEx.tmp
	@if (diff ESMF_LogEx.tmp ${ESMF_DIR}/${LOCDIR}/output/ESMF_LogEx.outuni) then \
	echo "PASS Logger test" > ${ESMC_TESTDIR}/ESMF_LogEx.log; else \
	echo "FAIL Logger test" > ${ESMC_TESTDIR}/ESMF_LogEx.log; fi
	@rm -f testlog stdlog.* logfile.* ESMF_LogEx.tmp optlog

