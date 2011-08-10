SUBDIRS = atm \
	chem \
	$(CLM_DIR) \
	inc \
	ocm \
	ocn_ml \
	meta \
	$(TEM_DIR)

DIST_SUBDIRS = atm \
	chem \
	clm/2.1 clm/3.5 \
	inc \
	ocm \
	ocn_ml \
	meta \
	tem/4.3.8 tem/4.4b

AM_CPPFLAGS = -I$(top_srcdir)/src/inc

if IGSM_BIN

bin_PROGRAMS = igsm22

# The driver routines are only built for the igsm22 executable.
igsm22_SOURCES = driver/atm2d_init_fixed.F \
	driver/atm2d_init_vars.F \
	driver/driver_IGSM22.F \
	driver/forward_step_atm2d.F

igsm22_FFLAGS = -Mdalign -Msave

# Use the PGI Fortran, rather than C/C++, linker.
igsm22_LINK = $(FC) -fastsse -o $@

# Libraries to be linked. The order of these is important and should not be
# disturbed: 
#  1. Multilayer ocean, atmosphere, ocean carbon
#  2. Appropriate version of the Common Land Model (as required)
#  3. Appropriate version of the Terrestrial Ecosystem Model (as required)
#  4. Atmospheric chemistry, urban airshed metamodel
#  5. NetCDF libraries
#  6. Standard C, C++ and Fortran libraries
# FIXME: determine whether all three NetCDF libraries are necessary
# FIXME: determine which libraries from the last line are redundant
igsm22_LDADD = ocn_ml/libml.a \
	atm/libatm.a \
	ocm/libocm.a \
	$(CLM_LDADD) \
	$(TEM_LDADD) \
	chem/libchem.a \
	meta/libmeta.a \
	$(NETCDF_LDADD) \
	$(FCLIBS) -lstd -lC -lgcc -lc -lstdc++

endif
