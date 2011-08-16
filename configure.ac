# For the meanings of macros AC_*, AH_*, and AM_, see the documentation for
# autoconf and automake:
# - http://www.gnu.org/software/autoconf/manual/autoconf.html
# - http://www.gnu.org/software/automake/manual/automake.html

### Initialization
AC_PREREQ([2.66])
AC_INIT(igsm, [1.0], [globalchange@mit.edu])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/atm/atmosphere.F])
# Generate ctrparam.h
AC_CONFIG_HEADERS([src/inc/ctrparam.h])
# Keep a comment at the bottom of ctrparam.h(.in) indicating where the user can
# set custom symbols
AH_BOTTOM([/* Add custom IGSM preprocessor #defines below this comment. Some known (but
   undocumented) symbols are: ACCRI, DATA4TEM, FIXED_FOR, FIXED_PCO2, GHS_ALB,
   HR_DATA, IPCC_EMI, OLD_LAPGAS, and ORBITAL_FOR.                            */])


### Package options
# This section defines the command-line options accepted by 'configure'. For
# usage instructions, run './configure --help'.

# Package option for model configurations
AC_ARG_WITH([setup],
  [AS_HELP_STRING([--with-setup],
    [model version to compile: either 'climate' (implies --without-tem), \
    'historical' or 'fullchem' (default)])],
  [if test $withval != climate -a $withval != historical -a $withval != fullchem; then
     AC_MSG_ERROR([invalid value '$withval' for --with-setup; try --help])
   else
     model=$withval
   fi],
  [model=fullchem])

# Package option for ocean carbon
AC_ARG_ENABLE([ocncarbon],
  [AS_HELP_STRING([--enable-ocncarbon],
    [enable the ocean carbon model (default for --with-setup=historical and \
    fullchem)])],
  [ocncarbon=$enableval],
  [if test $model != climate; then
     ocncarbon=yes
   else
     ocncarbon=no
   fi])
# TODO: prevent building src/ocm when ocncarbon != yes

# Package option for the urban metamodel
AC_ARG_ENABLE([meta],
  [AS_HELP_STRING([--enable-meta],
    [enable the urban metamodel (default for --with-setup=fullchem)])],
  [meta=$enableval],
  [if test $model = fullchem; then
     meta=yes
   else
     meta=no
   fi])
# TODO: prevent building src/meta when meta != yes

# Package option for 2D or 3D ocean
AC_ARG_ENABLE([mitgcm],
  [AS_HELP_STRING([--enable-mitgcm],
    [instead of the 2D ocean, produce libraries suitable for linking with the \
    MITgcm three-dimensional ocean model, and don't produce an igsm22 \
    executable])],
  [ocean=3D],
  [ocean=2D])
# See src/Makefile.am
AM_CONDITIONAL([IGSM_BIN], [test $ocean = 2D])

# Package option for CLM version
AC_ARG_WITH([clm],
  [AS_HELP_STRING([--with-clm],
    [select version of the Common Land Model: ARG=2.1 (default) or 3.5; \
     --without-clm implies --without-tem])],
  [CLM_VERSION=$withval],
  [CLM_VERSION=2.1])
if test $CLM_VERSION != no; then
  # Define for ctrparam.h
  AC_DEFINE(CLM, 1, [Coupling CLM])
  # Compiler and linker flags: pass to src/Makefile.am
  AC_SUBST(CLM_DIR, [clm/$CLM_VERSION])
  AC_CONFIG_LINKS([src/$CLM_DIR:$srcdir/src/$CLM_DIR])
  CLM_LDADD="-L$CLM_DIR -lclm"
  if test $CLM_VERSION = 2.1; then
    CLM_LDADD+=" -lesmf"
  fi
  AC_SUBST(CLM_LDADD)
fi

# Package option for TEM version
AC_ARG_WITH([tem],
  [AS_HELP_STRING([--with-tem],
    [select version of the Terrestrial Ecosystem Model: ARG=4.3.8 (default) or \
    4.4b; requires --with-clm])],
  [TEM_VERSION=$withval],
  [if test $model = climate -o $CLM_VERSION = no; then
     TEM_VERSION=no
   else
     TEM_VERSION=4.3.8
   fi])
if test $CLM_VERSION = no -a $TEM_VERSION != no; then
  AC_MSG_ERROR([cannot give both --with-tem and --without-clm])
elif test $TEM_VERSION != no; then
  # Defines for ctrparam.h
  AC_DEFINE(CPL_TEM, 1, [Coupling chemistry])
  AC_DEFINE(CPL_NEM, 1, [Coupling natural emission model])
  # Compiler and linker flags: pass to src/Makefile.am
  AC_SUBST(TEM_DIR, [tem/$TEM_VERSION])
  AC_SUBST(TEM_LDADD, ["-L$TEM_DIR -ltem"])
fi

# Package option for user source overrides
AC_ARG_WITH([usr-src],
  [AS_HELP_STRING([--with-usr-src],
    [ARG is the path to a directory containing replacements for IGSM source \
    files, or 'src/usr' if not given])],
  [usr_srcdir=`eval "echo $withval"`],
  [usr_srcdir=no])
if test "$usr_srcdir" = yes; then
  # user gave --with-usr-src but no directory; assume it's src/usr
  usr_srcdir="$srcdir/src/usr"
elif test "$usr_srcdir" = no -o ! -d "$usr_srcdir"; then
  # something besides --without-usr-src was given, but doesn't exist
  AC_MSG_NOTICE([user source directory '$usr_srcdir' does not exist])
  usr_srcdir=""
fi


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
NETCDF_LDADD="-lnetcdff -lnetcdf -lnetcdf_c++"
AC_SUBST(NETCDF_LDADD)
# Include paths to NetCDF when linking. The variable LIB_NETCDF should be
# defined by 'module load netcdf' or in some other manner.
LDFLAGS+=" -L${LIB_NETCDF}"
FFLAGS+=" -r8 -i4"
CXXFLAGS+=" -Mscalarsse -fpic"


### User source overrides
# only process the overrides if --with-usr-src was provided and a valid path was
# given
if test -n "$usr_srcdir"; then
  usr_links=""
  # consider each of the user's files in $usr_srcdir as a candidate
  for f in `find $usr_srcdir -type f`; do
    # look for a file elsewhere in src/ with the same name
    bn=`basename $f`
    orig=`find $srcdir/src -path $usr_srcdir -prune -o \( -name "$bn*" -print0 \)`
    if test -n "$orig"; then
      # a file was found 
      # path fragment for the target
      opath=`dirname $orig | sed "s@$srcdir/src/@@"`
      # add it to the list of files to be symlinked
      usr_links+="src/$opath/$bn:$f"
    fi
  done
  # let configure create the actual symlinks
  AC_CONFIG_LINKS([$usr_links])
fi


### Defines.
# Common defines for all model configurations
AC_DEFINE(N_LAT, 46, [Number of latitudinal grid points])
AC_DEFINE(N_LEV, 11, [Number of vertical pressure layers])
AC_DEFINE(N_LON0, 72, [TODO: document this symbol])

# Defines depending on 2D vs. 3D ocean
if test $ocean = 2D; then
  AC_DEFINE(ML_2D, 1, [Two-dimensional ocean])
  AC_DEFINE(N_LEVOCEAN, 11, [Number of vertical layers in diffusive ocean model
    (either 11 or 12)])
else
  AC_DEFINE(OCEAN_3D, 1, [Three-dimensional ocean])
fi

# Defines that depend on whether chemistry is included
if test $model = fullchem; then
  AC_DEFINE(CPL_CHEM, 1, [Coupling chemistry])
  # NYR_CHEM, PREDICTED_* always accompany CPL_CHEM
  AC_DEFINE(NYR_CHEM, 1, [Number of year of chemistry model integration (124 for
    period from 1977 to 2100)])
  AC_DEFINE(PREDICTED_GASES, 1, [Use predicted concentrations of chemical
    species by chemistry model to calculate radiative forcing (ifghgpredict)])
  AC_DEFINE(PREDICTED_AEROSOL, 1, [Use predicted aerosol concentrations by
    chemistry model to calculate radiative forcing (ifaerpredict)])
  AC_DEFINE(PREDICTED_BC, 1, [Use predicted black carbon concentrations by
    chemistry model to calculate radiative forcing])
  # more options for chemistry runs only
  AC_DEFINE(INC_3GASES, 1, [Including "3 gases"])
  AC_DEFINE(O3_RAD, 1, [Use predicted ozone in radiation module])
  if test $meta = yes; then
    AC_DEFINE(CPL_META, 1, [Coupling parameterized urban airshed model (meta)])
  fi
  if test $ocncarbon = yes; then
    AC_DEFINE(CPL_OCEANCO2, 1, [Coupling ocean CO2 model])
  fi
fi

# Defines for historical mode only
if test $model = historical; then
  AC_DEFINE(O3_FORCED, 1, [Use forced ozone scenario (ifo3forced)])
  AC_DEFINE(SVI_ALBEDO, 1, [Modify albedo to simulate sulfate aerosol forcing
    (ifsulfalb)])
fi

if test $model = historical -o $model = fullchem; then
  AC_DEFINE(VOL_AER, 1, [TODO: document this symbol])
fi

# Undocumented defines
if false; then
  AC_DEFINE(ACCRI, 1, [TODO: document this symbol])
  AC_DEFINE(DATA4TEM, 1, [TODO: document this symbol])
  AC_DEFINE(FIXED_FOR, 1, [TODO: document this symbol])
  AC_DEFINE(FIXED_PCO2, 1, [TODO: document this symbol])
  AC_DEFINE(GHS_ALB, 1, [TODO: document this symbol])
  AC_DEFINE(HR_DATA, 1, [TODO: document this symbol])
  AC_DEFINE(IPCC_EMI, 1, [TODO: document this symbol])
  AC_DEFINE(OLD_LAPGAS, 1, [TODO: document this symbol])
  AC_DEFINE(ORBITAL_FOR, 1, [TODO: document this symbol])
fi


### Output Makefiles
AC_CONFIG_FILES([Makefile
	src/Makefile
	src/atm/Makefile
	src/chem/Makefile
	src/inc/Makefile
	src/meta/Makefile
	src/ocm/Makefile
	src/ocn_ml/Makefile
	src/tem/4.3.8/Makefile
	src/tem/4.4b/Makefile
	util/Makefile
	util/analysis/Makefile
	util/eppa/Makefile])
AC_OUTPUT
