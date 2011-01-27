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
MODULES = tem ml chem meta atm
LIBS = $(patsubst %,lib%.a,$(MODULES))
clm_LIBS = -L$(realpath ./clm) -lclm -lesmf -L$(LIB_NETCDF) -lnetcdf -lnetcdff

atm_SRCS = $(wildcard atm/*.F)
atm_OBJS = $(patsubst %.F,%.o,$(atm_SRCS))

chem_SRCS = $(wildcard chem/*.F)
chem_OBJS = $(patsubst %.F,%.o,$(chem_SRCS))

meta_SRCS = $(wildcard meta/*.F)
meta_OBJS = $(patsubst %.F,%.o,$(meta_SRCS))

ml_SRCS = $(wildcard ocn_ml/*.F)
ml_OBJS = $(patsubst %.F,%.o,$(ml_SRCS))

tem_SRCS = $(wildcard tem/*.cpp)
tem_OBJS = $(patsubst %.cpp,%.o,$(tem_SRCS))

OBJECTS = $(atm_OBJS) $(chem_OBJS) $(meta_OBJS) $(ml_OBJS) $(tem_OBJS)

all: $(MODULES) clm
	$(FC) -fastsse -o $(TARGET) $(clm_LIBS) $(LIBS)

atm: libatm.a($(atm_OBJS))
chem: libchem.a($(chem_OBJS))
meta: libmeta.a($(meta_OBJS))
ml: libml.a($(ml_OBJS))
tem: libtem.a($(tem_OBJS))

clm: clm/libclm.a clm/libesmf.a

clm/%.a:
	cd clm && make libs

clean:
	rm -f $(LIBS) $(OBJECTS) $(TARGET)
	cd clm && make clean

