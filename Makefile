###################################################################
##  A makefile to link the module libraries                      ##
##			---------------------			 ##
## 	  MIT Joint Program for Global Change	                 ##
###################################################################

# directory for header files
INCDIR = ./inc

# compilers and flags for PGI/LINUX
FC       = pgf90
FFLAGS   = -I$(INCDIR) -r8 -i4
CXX      = pgCC
CXXFLAGS = -Mscalarsse -fpic

# lists of source files for IGSM modules, not including CLM
atm_SRCS = $(wildcard atm/*.F)
chem_SRCS = $(wildcard chem/*.F)
meta_SRCS = $(wildcard meta/*.F)
ml_SRCS = $(wildcard ocn_ml/*.F)
ocm_SRCS = $(wildcard ocm/*.F)
tem_SRCS = $(wildcard tem/*.cpp)

# source files for the emissions preprocessor
emiprep_SRCS = $(wildcard emiprep/*.F90)

# targets for IGSM. The syntax 'libfoo.a(foo.o)' compiles foo.F and updates it
# in the library libfoo.o
igsm_OBJECTS = libatm.a($(atm_SRCS:.F=.o)) \
	libchem.a($(chem_SRCS:.F=.o)) \
	clm/libclm.a \
	clm/libesmf.a \
	libmeta.a($(meta_SRCS:.F=.o)) \
	libml.a($(ml_SRCS:.F=.o)) \
	libocm.a($(ocm_SRCS:.F=.o)) \
	libtem.a($(tem_SRCS:.cpp=.o))

# compile flags for specific parts of the model
libatm.a(%.o): FFLAGS += -Mdalign -Msave
libchem.a(%.o): FFLAGS += -Mdalign -fast
libml.a(%.o): FFLAGS += -Mdalign -Msave
libocm.a(%.o): FFLAGS += -Msave

# libraries to link for IGSM — order is important
igsm_LIBS = -L. -lml -latm -locm -Lclm -lclm -lesmf -ltem -lchem -lmeta \
	-L$(LIB_NETCDF) -lnetcdff -lnetcdf -lnetcdf_c++ \
	-lstd -lC -lm -lpgc -lgcc -lc -lstdc++

# excutables to build
all: igsm22 prep

igsm22: $(igsm_OBJECTS)
	# link the resulting libraries
	$(FC) -fastsse -o $@ $(igsm_LIBS)

clm/%.a:
	# trigger the makefile in the CLM subdirectory
	cd clm && make libs

prep: $(emiprep_SRCS:.F90=.o)
	$(FC) -fastsse -o $@ $^

clean:
	rm -f igsm22 prep *.a
	cd clm && make clean

# Pattern rule for .F90 files. Rules for .c, .F etc. are defined by default.
%.o: %.F90
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $< -o $@

