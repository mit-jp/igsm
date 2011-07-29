# For the meanings of macros AC_ and AM_, see the documentation for autoconf and
# automake:
# - http://www.gnu.org/software/autoconf/manual/autoconf.html
# - http://www.gnu.org/software/automake/manual/automake.html
AC_PREREQ([2.66])
AC_INIT(igsm, [1.0], [globalchange@mit.edu])
AM_INIT_AUTOMAKE
#AC_CONFIG_SRCDIR([src/clm/2.1/misc.h])
AC_CONFIG_HEADERS([config.h])


### Package options
# This section defines the command-line options accepted by 'configure'. For
# usage instructions, run './configure --help'.

# Package option for CLM version. De
AC_ARG_WITH([clm],
  [AS_HELP_STRING([--with-clm],
    [select version of the Common Land Model: 2.1 (default) or 3.5.])],
  [CLM_VERSION=$withval],
  [CLM_VERSION=2.1])
if test "$CLM_VERSION" != no; then
  CLM_DIR="clm/$CLM_VERSION"
fi
AM_CONDITIONAL(WITH_CLM, [test "$CLM_VERSION" != no])
AC_SUBST(CLM_DIR)
AC_SUBST(CLM_VERSION)

# Package option for TEM version
AC_ARG_WITH([clm],
  [AS_HELP_STRING([--with-tem],
    [select version of the Terrestrial Ecosystem Model: 4.3.8 (default) or 4.4b.])],
  [TEM_VERSION=$withval],
  [TEM_VERSION=4.3.8])
if test "$TEM_VERSION" != no; then
  TEM_DIR="tem/$TEM_VERSION"
fi
AM_CONDITIONAL(WITH_TEM, [test "$TEM_VERSION" != no])
AC_SUBST(TEM_DIR)
AC_SUBST(TEM_VERSION)


### Build environment
# TODO: Figure out how to invoke 'module load' when the user runs 'configure'.
# Currently configuration will fail if the user hasn't loaded the pgi and netcdf
# modules.

# C++ compiler
AC_PROG_CXX(pgCC)
# C compiler
AC_PROG_CC(pgcc)
# Fortran 77 compiler
AC_PROG_F77
# TODO: the following three checks may not be necessary. Remove if possible.
AC_F77_DUMMY_MAIN
AC_F77_LIBRARY_LDFLAGS
AC_F77_MAIN
# Fortran 90+ compiler
AC_PROG_FC(pgf90)
# FIXME: this is cheating. The IGSM modules atm, ocm, ocn_ml and meta contain
# Fortran files with the extension '.F', which typically denotes Fortran 77, but
# must be compiled with the pgf90 compiler. It's not clear which variant of
# Fortran these are.
F77=$FC
# TODO: the following three checks may not be necessary. Remove if possible.
AC_FC_DUMMY_MAIN
AC_FC_LIBRARY_LDFLAGS
AC_FC_MAIN
AC_PROG_RANLIB
AC_PROG_MAKE_SET

# Checks for header files.
AC_CHECK_HEADERS([netcdf.h fcntl.h malloc.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_HEADER_STDC
AC_C_INLINE

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([gettimeofday memset modf pow strcasecmp strchr])

# Custom flags
# Include paths to NetCDF when linking. LIB_NETCDF is defined by 'module load
# netcdf'.
LDFLAGS="$LDFLAGS -L${LIB_NETCDF}"
FFLAGS="$FFLAGS -r8 -i4"
# FIXME: find a way to avoid using pwd here (use an autoconf variable)
CPPFLAGS="$CPPFLAGS -I`pwd`/src/inc"
CXXFLAGS="$CXXFLAGS -Mscalarsse -fpic"

### Defines.
# The current effect of these lines is that config.h will contain the same set
# of #defines as ctrparam.h.
# TODO: choose dynamically or define command-line options, as appropriate.
AC_DEFINE(N_LAT, 46, [Number of latitudinal grid points])
AC_DEFINE(N_LEV, 11, [Number of vertical pressure layers])
AC_DEFINE(N_LONO, 72, [TODO: document this])
AC_DEFINE(N_LEVOCEAN, 11, [Number of vertical layers in diffusive ocean model (either 11 or 12)])

# Only one of the following two lines
AC_DEFINE(ML_2D, 1, [Two-dimensional ocean])
#AC_DEFINE(OCEAN_3D, 1, [Three-dimensional ocean])

#AC_DEFINE(CLM, 1, [Coupling CLM])
#AC_DEFINE(CPL_TEM, 1, [Coupling chemistry])
#AC_DEFINE(CPL_NEM, 1, [Coupling natural emission model])
# TODO: logic to reflect the following comment from ctrparam.h: "In the current
# model version if you choose to undefine CPL_CHEM, you have to also choose
# CPL_NEM, CPL_META, and CPL_OCEANCO2. In addition, you should manage to set
# data reading correctly if you want use PREDICTED_GASES and PREDICTED_AEROSOLS.
#AC_DEFINE(CPL_CHEM, 1, [Coupling chemistry])
AC_DEFINE(NYR_CHEM, 1, [Number of year of chemistry model integration (124 for period from 1977 to 2100)])

# Note: with CPL_CHEM undef'd, the following 5 options might not work properly:
#AC_DEFINE(CPL_META, 1, [Coupling parameterized urban airshed model (meta). ])
#AC_DEFINE(CPL_OCEANCO2, 1, [Coupling ocean CO2 model])
#AC_DEFINE(PREDICTED_GASES, 1, [Use predicted concentrations of chemical species by chemistry model to calculate radiative forcing (ifghgpredict)])
#AC_DEFINE(PREDICTED_AEROSOL, 1, [Use predicted aerosol concentrations by chemistry model to calculate radiative forcing (ifaerpredict)])
#AC_DEFINE(PREDICTED_BC, 1, [Use predicted black carbon concentrations by chemistry model to calculate radiative forcing])

#AC_DEFINE(INC_3GASES, 1, [Including "3 gases"])
#AC_DEFINE(O3_RAD, 1, [Use predicted ozone in radiation module])
#AC_DEFINE(O3_FORCED, 1, [Use forced ozone scenario (ifo3forced). If CPL_CHEM is def'd, this parameter is excluded])
#AC_DEFINE(SVI_ALBEDO, 1, [Modify albedo to simulate sulfate aerosol forcing (ifsulfalb). If CPL_CHEM is def'd, this parameter is excluded])

# TODO: document this
#AC_DEFINE(VOL_AER)


### Output Makefiles
AC_CONFIG_FILES([Makefile
	src/Makefile
	src/atm/Makefile
	src/chem/Makefile
	src/meta/Makefile
	src/ocm/Makefile
	src/ocn_ml/Makefile
	src/tem/4.3.8/Makefile
	src/tem/4.4b/Makefile
	util/Makefile
	util/analysis/Makefile
	util/eppa/Makefile])
AC_OUTPUT
