AC_PREREQ([2.66])
AC_INIT(igsm, [1.0], [globalchange@mit.edu])
AM_INIT_AUTOMAKE
#AC_CONFIG_SRCDIR([src/clm/2.1/misc.h])
AC_CONFIG_HEADERS([config.h])

# Package option for CLM version
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

# Checks for programs.
AC_PROG_CXX(pgCC)
AC_PROG_CC(pgcc)
AC_PROG_F77
AC_F77_DUMMY_MAIN
AC_F77_LIBRARY_LDFLAGS
AC_F77_MAIN
AC_PROG_FC(pgf90)
# FIXME: this is cheating
F77=$FC
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
LDFLAGS="$LDFLAGS -L${LIB_NETCDF}"
FFLAGS="$FFLAGS -r8 -i4"
CXXFLAGS="$CXXFLAGS -Mscalarsse -fpic"
CPPFLAGS="$CPPFLAGS -I`pwd`/src/inc"

AC_CONFIG_FILES([Makefile
	src/Makefile
	src/atm/Makefile
	src/chem/Makefile
	src/meta/Makefile
	src/ocm/Makefile
	src/ocn_ml/Makefile
	src/tem/4.3.8/Makefile
	src/tem/4.4b/Makefile])
AC_OUTPUT
