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

EXTRA_PROGRAMS = igsm22 igsm22ocm igsm23_driver igsm23_ocean

# -L$LIB_NETCDF must be inserted *before* $(FCLIBS). pgf90 reports -L/usr/lib64
# as a library directory, but it contains a gcc-compiled NetCDF that cannot be
# linked against PGI-compiled object files.
LDADD = $(MODULE_LDADD) \
	-L$(LIB_NETCDF) \
	$(FCLIBS) \
	-lhdf5_hl -lhdf5 -lpgmp -lnuma -lpthread -lpgf90 -lpgf90_rpm1 -lpgf902 \
	-lpgf90rtl -lpgftnrtl -lnspgc -lpgc -lrt -lm @IGSM_LDADD_EXTRA@ -lgcc \
	-lgcc_s -lc -lstdc++

if IGSM_2D

# no additional sources for the 22 executable
igsm22_SOURCES =

# Use the PGI Fortran, rather than C/C++, linker.
igsm22_LINK = $(FC) -fastsse -Mnorpath -o $@

igsm22ocm_SOURCES = driver_ocm/ocm.F \
	driver_ocm/init_carb_data.F

igsm22ocm_LINK = $(FC) -fastsse -Mnorpath -o $@

else !IGSM_2D

# no additional sources for the 23 executable
igsm23_driver_SOURCES =
igsm23_ocean_SOURCES =

# Use the MPI Fortran, rather than C/C++, linker
igsm23_driver_LINK = $(MPI_FC) -fastsse -o $@
igsm23_ocean_LINK = $(MPI_FC) -fastsse -o $@

endif !IGSM_2D
