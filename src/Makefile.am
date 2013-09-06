SUBDIRS = atm inc $(MODULE_DIRS)

DIST_SUBDIRS = atm \
	chem \
	clm/2.1 \
	clm/3.5 \
	inc \
	ocm \
	ocn_3d \
	ocn_ml \
	meta \
	tem/4.3.8 \
	tem/4.4c \
	tracers

AM_CPPFLAGS = -I$(top_srcdir)/src/inc

bin_PROGRAMS = $(BINARIES)

EXTRA_PROGRAMS = igsm22 igsm23

LDADD = $(MODULE_LDADD) $(FCLIBS) -lstd -lC -lgcc -lgcc_s -lc -lstdc++

if IGSM_2D

# The driver routines are only built for the igsm22 executable.
igsm22_SOURCES = driver/atm2d_init_fixed.F \
	driver/atm2d_init_vars.F \
	driver/driver_IGSM22.F \
	driver/atm2d_finish.F \
	driver/forward_step_atm2d.F

igsm22_FFLAGS = -Mdalign -Msave

# Use the PGI Fortran, rather than C/C++, linker.
igsm22_LINK = $(FC) -fastsse -o $@

# Libraries to be linked.
# FIXME: determine which libraries from the last set are redundant
#igsm22_LDADD = $(MODULE_LDADD) $(FCLIBS) -lstd -lC -lgcc -lc -lstdc++

else !IGSM_2D

# no additional sources for the 23 executable
igsm23_SOURCES =

# Use the MPI Fortran, rather than C/C++, linker
igsm23_LINK = $(MPI_FC) -fastsse -o $@

# libraries to be linked
#igsm23_LDADD = $(MODULE_LDADD) $(FCLIBS) -lstd -lC -lgcc -lc -lstdc++

endif !IGSM_2D
