###################################################################
##  A makefile to link the module libraries                      ##
##			---------------------			 ##
## 	  MIT Joint Program for Global Change	                 ##
###################################################################

INCDIR = ./inc

# --- for PGI/LINUX
FC       = pgf90
FFLAGS   = -I$(INCDIR) -r8 -i4 -Mdalign -Msave
CXX      = pgCC
CXXFLAGS = -Mscalarsse -fpic

TARGET = igsm22

# atm must be last for link to succeed?
MODULES = tem ocm ml chem meta atm

LIBS = $(patsubst %,lib%.a,$(MODULES))
clm_LIBS = -L$(realpath ./clm) -lclm -lesmf -L$(LIB_NETCDF) -lnetcdf -lnetcdff

atm_SRCS = $(wildcard atm/*.F)
atm_OBJS = $(patsubst %.F,%.o,$(atm_SRCS))

chem_SRCS = $(wildcard chem/*.F)
chem_OBJS = $(patsubst %.F,%.o,$(chem_SRCS))

# FIXME: not all of the .F90 files in emiprep are used. Why?
emiprep_SRCS = emiprep/eppanew_mod.F90 emiprep/testeppanew.F90
emiprep_OBJS = $(patsubst %.F90,%.o,$(emiprep_SRCS))

meta_SRCS = $(wildcard meta/*.F)
meta_OBJS = $(patsubst %.F,%.o,$(meta_SRCS))

ml_SRCS = $(wildcard ocn_ml/*.F)
ml_OBJS = $(patsubst %.F,%.o,$(ml_SRCS))

ocm_SRCS = $(wildcard ocm/*.F)
ocm_OBJS = $(patsubst %.F,%.o,$(ocm_SRCS))

tem_SRCS = $(wildcard tem/*.cpp)
tem_OBJS = $(patsubst %.cpp,%.o,$(tem_SRCS))

OBJECTS = $(atm_OBJS) $(chem_OBJS) $(clm_OBJS) $(emiprep_OBJS) $(meta_OBJS) $(ml_OBJS) $(ocm_OBJS) $(tem_OBJS)

all: $(TARGET) prep

$(TARGET): $(MODULES) clm
	$(FC) -fastsse -o $@ $(clm_LIBS) $(LIBS)

atm: libatm.a($(atm_OBJS))
chem: libchem.a($(chem_OBJS))
meta: libmeta.a($(meta_OBJS))
ml: libml.a($(ml_OBJS))
ocm: libocm.a($(ocm_OBJS))
tem: libtem.a($(tem_OBJS))

clm: clm/libclm.a clm/libesmf.a

clm/%.a:
	cd clm && make libs

#prep: FFLAGS = -I./emiprep -r8 -i4 -Mdalign -Msave
prep: $(emiprep_OBJS)
	$(FC) -fastsse -o $@ $^

clean:
	rm -f $(LIBS) $(OBJECTS) $(TARGET) prep
	cd clm && make clean

# Implicit rule for .F90 files. Rules for .c, .F etc. are defined by default.
%.o: %.F90
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $< -o $@
