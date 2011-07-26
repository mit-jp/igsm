/* $Id$ */

#include "ESMC.h"

#include "ESMF_conf.h"

#include "ESMC_App.h"

#include <stdarg.h>

#ifdef ESMC_HAVE_FORTRAN_UNDERSCORE
#define FORTRANUNDERSCORE
#endif


#ifdef MPI_BUILD_PROFILING
#ifdef FORTRANCAPS
#define esmc_appnew_ PESMC_LOGAPPNEW
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmc_appnew_ pesmc_appnew__
#elif !defined(FORTRANUNDERSCORE)
#define esmc_appnew_ pesmc_appnew
#else

#define esmc_appnew_ pesmc_appnew_
#endif
#else
#ifdef FORTRANCAPS
#define esmc_appnew_ ESMC_LOGAPPNEW
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmc_appnew_ esmc_appnew__
#elif !defined(FORTRANUNDERSCORE)
#define esmc_appnew_ esmc_appnew
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef FORTRANCAPS
#define esmc_appdelete_ PESMC_LOGAPPDELETE
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmc_appdelete_ pesmc_appdelete__
#elif !defined(FORTRANUNDERSCORE)
#define esmc_appdelete_ pesmc_appdelete
#else

#define esmc_appdelete_ pesmc_appdelete_
#endif
#else
#ifdef FORTRANCAPS
#define esmc_appdelete_ ESMC_LOGAPPDELETE
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmc_appdelete_ esmc_appdelete__
#elif !defined(FORTRANUNDERSCORE)
#define esmc_appdelete_ esmc_appdelete
#endif
#endif

void esmc_appnew_(ESMC_App *app, int *rc)
{
  *rc = ESMC_AppNew(app);
}

void esmc_appdelete_(ESMC_App *app, int *rc)
{
  *rc = ESMC_AppDelete(*app);
}
