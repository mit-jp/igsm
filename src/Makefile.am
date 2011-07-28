SUBDIRS = atm chem $(CLM_DIR) ocm ocn_ml meta $(TEM_DIR)
DIST_SUBDIRS = atm chem clm/2.1 clm/3.5 ocm ocn_ml meta tem/4.3.8 tem/4.4b

bin_PROGRAMS = igsm22

# There is no top-level source file for the IGSM; set _SOURCES blank so automake
# doesn't go looking for "igsm22.c".
igsm22_SOURCES =

# Use the PGI Fortran, rather than C/C++, linker.
igsm22_LINK = $(FC) -fastsse -o $@

# Libraries to be linked. The order of these is important and should not be
# disturbted: first the mulitlayer ocean, atmosphere and ocean carbon…
igsm22_LDADD = ocn_ml/libml.a \
	atm/libatm.a \
	ocm/libocm.a

# …then the appropriate version of the Common Land Model, as required…
# FIXME: check that -lclm -lesmf is the correct set for CLM 3.5
if WITH_CLM
igsm22_LDADD += -L$(CLM_DIR) -lclm -lesmf
endif

# …then the appropriate version of the Terrestrial Ecosystem Model, as required…
if WITH_TEM
igsm22_LDADD += -L$(TEM_DIR) -ltem
endif

# …finally the atmospheric chemistry, the urban airshed metamodel, the NetCDF
# libraries, and standard C, C++ and Fortran libraries.
# FIXME: determine whether all three NetCDF libraries are necessary
# FIXME: determine which libraries from the last line are redundant
igsm22_LDADD += chem/libchem.a \
	meta/libmeta.a \
	-lnetcdff -lnetcdf -lnetcdf_c++ \
	$(FCLIBS) -lstd -lC -lgcc -lc -lstdc++
