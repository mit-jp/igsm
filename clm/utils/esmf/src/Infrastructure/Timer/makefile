# $Id$ 
ALL: lib 

CFLAGS    = -Mlarge_arrays
FFLAGS    = -Mlarge_arrays
SOURCEC	  = ESMC_Timer.c ESMC_TimerF.c
SOURCEF	  = ESMF_TimerMod.F
SOURCEH	  = ESMC_Timer.h
OBJSC     = $(addsuffix .o, $(basename $(SOURCEC)))
OBJSF     = $(addsuffix .o, $(basename $(SOURCEF)))
TEXFILES  = $(addsuffix .tex, $(basename $(SOURCEF)))
MANSEC	  = Timeobj
LIBBASE	  = libesmf
DIRS	  = examples doc
LOCDIR	  = src/Infrastructure/Timer/
CLEANFILES = $(LIBNAME) $(addprefix doc/, $(addsuffix .tex, $(basename $(SOURCEF))))

include ${ESMF_DIR}/build/${ESMF_ARCH}/base






