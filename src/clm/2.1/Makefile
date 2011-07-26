#-----------------------------------------------------------------------
# Makefile for building CLM on various platforms.
#
# This makefile assumes the existence the file:
#
# Filepath      The directories, listed one per line, that contain the source
#               code required to build CLM.  This list is used the set the
#               VPATH variable which is used by GNU make to search for
#               dependencies (after it looks in the directory from which
#               it was invoked).  This list of directories, prepended with ./,
#               is also used to construct the list of search directories 
#               used by the preprocessor (as specified by -I command-line options).
#
# The following macros may be set in the user's environment:
#
# ROOTDIR       Root directory for the CLM source distribution.  If ROOTDIR is not
#               set the makefile searches for it in a file called "Rootdir".
# EXENAME       Name to call the executable.  Default: clm
# MODEL_EXEDIR  Directory to build the executable in.  Default: ./
# INC_NETCDF    Directory location of the NetCDF include files.  Default: /usr/local/include
# LIB_NETCDF    Directory location of the NetCDF library.  Default: /usr/local/lib
# INC_MPI       Directory location of the MPI include files.  Default: /usr/local/include
# LIB_MPI       Directory location of the MPI library.  Default: /usr/local/lib
# ESMF_ROOT     Root directory of ESMF source tree.  Default: $(ROOTDIR)/models/utils/esmf
# ESMF_BLD      Upper level build directory for ESMF
#               The ESMF .o files are built in $(ESMF_BLD)/obj
#               The ESMF .m files are built in arch dependent subdirs of $(ESMF_BLD)/mod
#               The ESMF library is built in architecture and optimization
#               dependent subdirectories of $(ESMF_BLD)/lib
#               Default: ./esmf
# DEBUG         Set to TRUE to turn on compiler debugging options.  Default: FALSE 
# SPMD          Whether to build in SPMD mode or not. 	[values TRUE FALSE]
#               Default: read from ./misc.h
#               N.B. To ensure consistency between the Makefile and misc.h the SPMD
#                    macro should not be set.  This forces it to be read from misc.h.
# SMP           Set to TRUE to enable building in SMP mode (uses OpenMP).  
#               Currently implemented for IBM, SGI, linux-pgf90. (default is TRUE on IBM and 
#               linux-pgf90, and depends on SPMD setting on SGI).
# USER_FC       Allow user to override the default Fortran compiler specified in Makefile.
# USER_CC       Allow user to override the default C compiler specified in Makefile (linux only).
# USER_FFLAGS   Additional compiler flags that the user wishes to set.
# NO_SWITCH     On Compaq if the hardware switch is not available
#               set this env variable to "TRUE".
#
# Note: The ESMF library is included in the CLM distribution in 
#       $ROOTDIR/src/utils/esmf and is built using this makefile.  
#------------------------------------------------------------------------

# Set up special characters
null  :=
space := $(null) $(null)
comma := $(null),$(null)

# Determine distribution root directory
ROOTDIR ?= $(if $(realpath Rootdir),$(shell cat Rootdir),.)
ROOTDIR := $(abspath $(ROOTDIR))

# Check for the NetCDF library and include directories 
LIB_NETCDF ?= /usr/local/lib
INC_NETCDF ?= /usr/local/include

# Check for the MPI library and include directories 
LIB_MPI ?= /usr/local/lib
INC_MPI ?= /usr/local/include

# Build the ESMF library
ESMF_ROOT ?= $(ROOTDIR)/utils/esmf
ESMF_BLD ?= $(shell (mkdir ./esmf 2>/dev/null;cd ./esmf;pwd))
ESMF_BOPT := O

ifeq ($(DEBUG),TRUE)
  ESMF_BOPT := g
endif

ESMF_MOD := $(ESMF_BLD)/mod/mod$(ESMF_BOPT)
ESMF_LIB := $(ESMF_BLD)/lib/lib$(ESMF_BOPT)

# Check for directory in which to put executable
MODEL_EXEDIR ?= $(abspath .)

# Check for name of executable
EXENAME ?= clm

# Check if SPMD is defined in "misc.h"
# Ensure that it is defined and not just "undef SPMD" set in file
ifndef SPMD
  SPMDSET := $(if $(realpath misc.h),$(shell /bin/grep SPMD misc.h))
  ifneq ($(findstring define,$(SPMDSET)),)
    SPMD := TRUE
  else
    SPMD := FALSE
  endif
endif

# Load dependency search path.
dirs := . $(if $(realpath Filepath),$(shell cat Filepath))

# Set cpp search path, include netcdf
cpp_dirs := $(dirs) $(INC_NETCDF) $(INC_MPI)
cpp_path := $(foreach dir,$(cpp_dirs),-I$(dir)) # format for command line

# Expand any tildes in directory names. Change spaces to colons.
VPATH := $(foreach dir,$(cpp_dirs),$(wildcard $(dir)))
VPATH := $(ESMF_LIB)/$(ESMF_ARCH):$(subst $(space),:,$(VPATH))

#------------------------------------------------------------------------
# Primary target: build the model
#------------------------------------------------------------------------
TARGET  := $(MODEL_EXEDIR)/$(EXENAME)
CPPDEF  := -DNO_SHR_VMATH -DHIDE_SHR_MSG 

all: $(TARGET) libs
libs: libclm.a libesmf.a

# 'depends.mk' contains a list of make rules for the files in Srcfiles
depends.mk: Filepath Srcfiles
	./depends.pl $^ >$@

# 'Filepath' contains a list of directories to search for the CLM source
Filepath:
	./premake.csh
	# copy source files appropriate for TEM437 coupling
	cp $(wildcard main/TEM437/*.F90) .
	sed 's@$$ROOTDIR@'$(ROOTDIR)'@' Filepath.in >Filepath

# 'Srcfiles' contains a list of all .F, .c and .F90 files
Srcfiles: Filepath
	find . `cat $^` -maxdepth 1 \( -name "*.[Fc]" -o -name "*.F90" \) \
	-printf "%f\n" | sort | uniq >$@

SOURCES := $(shell cat Srcfiles)
OBJS    := $(addsuffix .o, $(basename $(SOURCES)))

$(OBJS): depends.mk preproc.h lnd.stdin

# create the file 'preproc.h'
LSMLON := 1
LSMLAT := 46
preproc.h:
	sed -e 's/$$LSMLON/'$(LSMLON)'/' -e 's/$$LSMLAT/'$(LSMLAT)'/' \
	preproc.h.in >preproc.h

# create the file 'lnd.stdin'
CSMDATA := $(ROOTDIR)/inputdata
lnd.stdin:
	sed -e 's@$$CSMDATA@'$(CSMDATA)'@' lnd.stdin.in >lnd.stdin

$(TARGET): $(OBJS)
	$(FC) -o $@ $^ $(LFLAGCC) -L$(LIB_NETCDF) -lnetcdf -lnetcdff \
	-L$(ESMF_LIB)/$(ESMF_ARCH) -lesmf $(LDFLAGS)

libclm.a: libclm.a($(filter-out program_isgm.o,$(OBJS)))

debug: $(OBJS)
	echo "FFLAGS: $(FFLAGS)"
	echo "LDFLAGS: $(LDFLAGS)"
	echo "OBJS: $(OBJS)"

.SUFFIXES:
.SUFFIXES: .F .F90 .c .o

#***********************************************************************
#********** Architecture-specific flags and rules***********************
#***********************************************************************

# Determine platform 
UNAMES := $(shell uname -s)

#------------------------------------------------------------------------
# Linux
#------------------------------------------------------------------------

ifeq ($(UNAMES),Linux)

ESMF_ARCH = linux

# if the compiler is not specified, choose PGI pgf90

ifndef USER_FC
  FC := pgf90
else
  FC := $(USER_FC)
endif

# if USER_CC is set, use it.  Otherwise use pgcc if pgf90 is the Fortran compiler,
# and cc if not

ifndef USER_CC
  ifeq ($(FC),pgf90)
    CC := pgcc
  else
    CC := cc
  endif
else
  CC := $(USER_CC)
endif

# Figure out the ESMF architecture.
ifeq ($(FC),pgf90)
  ifeq ($(CC),pgcc)
    ESMF_ARCH = linux_pgi
  else
    ESMF_ARCH = linux_gnupgf90
  endif
endif

# if not using pgcc, define USE_GCC to disable threading in the timing library

CFLAGS = $(cpp_path) -DLINUX $(CPPDEF)
ifeq ($(CC),pgcc)
#  CFLAGS += -fast -mcmodel=medium
  CFLAGS += -fast -Mlarge_arrays
else
  CFLAGS += -DUSE_GCC
endif

# pgf90

ifeq ($(FC),pgf90)

  ifndef SMP
    SMP := TRUE
  endif

  FFLAGS    := $(cpp_path) -C -r8 -i4 $(CPPDEF) -Mlarge_arrays \
               -Mrecursive -Mdalign -Mextend -DLINUX -DPGF90 \
               -I$(ESMF_MOD)/$(ESMF_ARCH)
  FREEFLAGS := -Mfree

  ifeq ($(DEBUG),TRUE)
    FFLAGS      += -g -Ktrap=fp -Mbounds
    SPEC_FFLAGS := $(FFLAGS) 
  else
    SPEC_FFLAGS := $(FFLAGS) 
    FFLAGS += -fast
  endif

endif

# lf95
#
# Note that as of lf95 version 6.1 threading does NOT work because of
# ridiculously small per thread stacksize limits.
#
#  -CcdRR8 is an undocumented flag which promotes only vars declared "real", not "real*8"
# --trace produces a call traceback on abort
# --trap causes code to stop on divide by zero or overflow exceptions
# --pca prevents overwriting constant arguments
# --staticlink enables transportability of the executable
# --chk for basic compiler checking (a,e,s,u,x)
# --chkglobal for global checking

ifeq ($(FC),lf95)
  ESMF_ARCH = linux_lf95
  FFLAGS   = $(cpp_path) -CcdRR8 $(CPPDEF) --trace --trap --wide -DLINUX -I$(ESMF_MOD)/$(ESMF_ARCH)
  SPEC_FFLAGS := $(FFLAGS)
  F90FLAGS = $(FFLAGS)
  LDFLAGS  = --staticlink
  ifeq ($(DEBUG),TRUE)
    FFLAGS += -g --chk --pca
  else
    FFLAGS  += -O
  endif
endif

# FFC f90

ifeq ($(FC),f90)
  FFLAGS       = -CcdRR8 -Cpp -X9 -fw -Am -Wa,-W $(cpp_path) $(CPPDEF) -DLINUX \
                 -I$(ESMF_MOD)/$(ESMF_ARCH)
  LDFLAGS      = 
  ifeq ($(DEBUG),TRUE)
    FFLAGS += -g -H aseu
  else
    FFLAGS += -Kfast,eval,fastlib,auto
  endif
endif

# Flags common to all compilers

ifeq ($(SPMD),TRUE)
  FFLAGS  += -I$(INC_MPI) 
  LDFLAGS += -L$(LIB_MPI) -lmpich
else
  FFLAGS  += -DHIDE_MPI 
endif

ifeq ($(FC),lf95)
# do not compile with debug flags under lahey
wrap_nf.o: wrap_nf.F90
	$(FC) -c $(SPEC_FFLAGS) $<
wrap_mpi.o: wrap_mpi.F90
	$(FC) -c $(SPEC_FFLAGS) $<
endif

#  LFLAGCC= -Mscalarsse -lstd -lC -lm -lpgc -lgcc -lc -lstdc++ -I/usr/pkg/x86-64/pgi616/linux86/6.1/include/CC -I/usr/pkg/x86-64/pgi616/linux86/6.1/include -I/usr/lib/gcc-lib/x86_64-redhat-linux/3.2.3/include -I/usr/include -L/usr/lib/gcc-lib/x86_64-redhat-linux/3.2.3 -L/usr/pkg/x86-64/pgi616/linux86/6.1/lib -L/usr/pkg/x86-64/pgi616/linux86/6.1/include/CC -L/usr/pkg/x86-64/pgi616/linux86/6.1/include
  LFLAGCC= -Mscalarsse -lstd -lC -lm -lpgc -lgcc -lc -lstdc++  -I/usr/lib/gcc-lib/x86_64-redhat-linux/3.2.3/include -I/usr/include -L/usr/lib/gcc-lib/x86_64-redhat-linux/3.2.3 
#  LFLAGCC= -Mscalarsse -lstd -lC -lm -lpgc -lgcc -lc -lstdc++ 

.F90.o:
	$(FC) -c $(FREEFLAGS) $(FFLAGS) $<
.F.o:
	$(FC) -c $(FFLAGS) $<
.c.o:
	$(CC) -c $(CFLAGS) $<

endif


#------------------------------------------------------------------------
# Targets/rules that depend on architecture specific variables.
#------------------------------------------------------------------------

# The library is not made to be built in parallel, so we must null passed
# options and specify only one job is to run.
libesmf.a:
	cd $(ESMF_ROOT) && \
	$(MAKE) -j 1 BOPT=$(ESMF_BOPT) ESMF_BUILD=$(ESMF_BLD) ESMF_DIR=$(ESMF_ROOT) ESMF_ARCH=$(ESMF_ARCH)
	cp $(ESMF_LIB)/$(ESMF_ARCH)/$@ ./

time_manager.o : libesmf.a

# Add user defined compiler flags if set, and replace FC if USER option set.
FFLAGS  += $(USER_FFLAGS)
ifdef USER_FC
  FC := $(USER_FC)
endif

.PHONY: clean realclean
clean:
	rm -rf esmf libesmf.a libclm.a
# Files generated by premake.sh
	rm -f Filepath *.F90 misc.h preproc.h lnd.stdin
	rm -f depends.mk Srcfiles *.o *.mod *.stb $(TARGET)

realclean: clean

ifneq ($(MAKECMDGOALS),clean)
  -include depends.mk
endif

