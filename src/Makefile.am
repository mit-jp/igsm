SUBDIRS = atm inc $(MODULE_DIRS)

DIST_SUBDIRS = atm \
	chem \
	clm \
	inc \
	ocm \
	ocn_3d \
	ocn_ml \
	meta \
	tem \
	tracers

AM_CPPFLAGS = -I$(top_srcdir)/src/inc

bin_PROGRAMS = $(BINARIES)

EXTRA_PROGRAMS = igsm22 igsm23_driver igsm23_ocean

LDADD = $(MODULE_LDADD) $(FCLIBS) -lhdf5_hl -lhdf5 -lpgmp -lnuma -lpthread \
	-lpgf90 -lpgf90_rpm1 -lpgf902 -lpgf90rtl -lpgftnrtl -lnspgc -lpgc -lrt \
	-lm -lstd -lC -lgcc -lgcc_s -lc -lstdc++

if IGSM_2D

# no additional sources for the 22 executable
igsm22_SOURCES =

# Use the PGI Fortran, rather than C/C++, linker.
igsm22_LINK = $(FC) -fastsse -o $@

# Libraries to be linked.
# FIXME: determine which libraries from the last set are redundant
#igsm22_LDADD = $(MODULE_LDADD) $(FCLIBS) -lstd -lC -lgcc -lc -lstdc++

else !IGSM_2D

# no additional sources for the 23 executable
igsm23_driver_SOURCES =
igsm23_ocean_SOURCES =

# Use the MPI Fortran, rather than C/C++, linker
igsm23_driver_LINK = $(MPI_FC) -fastsse -o $@
igsm23_ocean_LINK = $(MPI_FC) -fastsse -o $@

# libraries to be linked
#igsm23_driver_LDADD = $(MODULE_LDADD) $(FCLIBS) -lstd -lC -lgcc -lc -lstdc++
#igsm23_ocean_LDADD = $(MODULE_LDADD) $(FCLIBS) -lstd -lC -lgcc -lc -lstdc++

endif !IGSM_2D
