# $Id$ 
ALL: lib 

CFLAGS    = -Mlarge_arrays
FFLAGS    = -Mlarge_arrays
SOURCEC	  = ESMC_Log.c ESMC_LogF.c
SOURCEF	  = ESMF_LogMod.F
SOURCEH	  = ESMC_Log.h
OBJSC     = $(addsuffix .o, $(basename $(SOURCEC)))
OBJSF     = $(addsuffix .o, $(basename $(SOURCEF)))
TEXFILES  = $(addsuffix .tex, $(basename $(SOURCEF)))
MANSEC	  = Timeobj
LIBBASE	  = libesmf
DIRS	  = examples doc
LOCDIR	  = src/Infrastructure/Log/
CLEANFILES = $(LIBNAME) $(addprefix doc/, $(addsuffix .tex, $(basename $(SOURCEF))))

include ${ESMF_DIR}/build/${ESMF_ARCH}/base






