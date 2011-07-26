# $Id$ 
ALL: lib 

CFLAGS    = -Mlarge_arrays
FFLAGS    = -Mlarge_arrays
SOURCEC	  = ESMC_Machine.c ESMC_MachineF.c
SOURCEF	  = ESMF_MachineMod.F
SOURCEH	  = ESMC_Machine.h
OBJSC     = $(addsuffix .o, $(basename $(SOURCEC)))
OBJSF     = $(addsuffix .o, $(basename $(SOURCEF)))
TEXFILES  = $(addsuffix .tex, $(basename $(SOURCEF)))
MANSEC	  = Timeobj
LIBBASE	  = libesmf
DIRS	  = doc
LOCDIR	  = src/Infrastructure/Machine/
CLEANFILES = $(LIBNAME) $(addprefix doc/, $(addsuffix .tex, $(basename $(SOURCEF))))

include ${ESMF_DIR}/build/${ESMF_ARCH}/base






