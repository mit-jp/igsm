/* $Id$ */

#include "ESMC.h"

#include "ESMC_Log.h"

#include "ESMC_Fortran.h"

#include <stdarg.h>
#include <malloc.h>

#ifdef ESMC_HAVE_FORTRAN_UNDERSCORE
#define FORTRANUNDERSCORE
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef FORTRANCAPS
#define esmc_logprint_ PESMC_LOGPRINT
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmc_logprint_ pesmc_logprint__
#elif !defined(FORTRANUNDERSCORE)
#define esmc_logprint_ pesmc_logprint
#else

#define esmc_logprint_ pesmc_logprint_
#endif
#else
#ifdef FORTRANCAPS
#define esmc_logprint_ ESMC_LOGPRINT
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmc_logprint_ esmc_logprint__
#elif !defined(FORTRANUNDERSCORE)
#define esmc_logprint_ esmc_logprint
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef FORTRANCAPS
#define esmc_lognew_ PESMC_LOGNEW
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmc_lognew_ pesmc_lognew__
#elif !defined(FORTRANUNDERSCORE)
#define esmc_lognew_ pesmc_lognew
#else

#define esmc_lognew_ pesmc_lognew_
#endif
#else
#ifdef FORTRANCAPS
#define esmc_lognew_ ESMC_LOGNEW
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmc_lognew_ esmc_lognew__
#elif !defined(FORTRANUNDERSCORE)
#define esmc_lognew_ esmc_lognew
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef FORTRANCAPS
#define esmc_logdelete_ PESMC_LOGDELETE
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmc_logdelete_ pesmc_logdelete__
#elif !defined(FORTRANUNDERSCORE)
#define esmc_logdelete_ pesmc_logdelete
#else

#define esmc_logdelete_ pesmc_logdelete_
#endif
#else
#ifdef FORTRANCAPS
#define esmc_logdelete_ ESMC_LOGDELETE
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmc_logdelete_ esmc_logdelete__
#elif !defined(FORTRANUNDERSCORE)
#define esmc_logdelete_ esmc_logdelete
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef FORTRANCAPS
#define esmf_logprint_ PESMF_LOGPRINT
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmf_logprint_ pesmf_logprint__
#elif !defined(FORTRANUNDERSCORE)
#define esmf_logprint_ pesmf_logprint
#else

#define esmf_logprint_ pesmf_logprint_
#endif
#else
#ifdef FORTRANCAPS
#define esmf_logprint_ ESMF_LOGPRINT
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmf_logprint_ esmf_logprint__
#elif !defined(FORTRANUNDERSCORE)
#define esmf_logprint_ esmf_logprint
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef FORTRANCAPS
#define esmf_logsetconfig_ PESMF_LOGSETCONFIG
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmf_logsetconfig_ pesmf_logsetconfig__
#elif !defined(FORTRANUNDERSCORE)
#define esmf_logsetconfig_ pesmf_logsetconfig
#else

#define esmf_logsetconfig_ pesmf_logsetconfig_
#endif
#else
#ifdef FORTRANCAPS
#define esmf_logsetconfig_ ESMF_LOGSETCONFIG
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmf_logsetconfig_ esmf_logsetconfig__
#elif !defined(FORTRANUNDERSCORE)
#define esmf_logsetconfig_ esmf_logsetconfig
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef FORTRANCAPS
#define esmc_logsetstate_ PESMC_LOGSETSTATE
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmc_logsetstate_ pesmc_logsetstate__
#elif !defined(FORTRANUNDERSCORE)
#define esmc_logsetstate_ pesmc_logsetstate
#else

#define esmc_logsetstate_ pesmc_logsetstate_
#endif
#else
#ifdef FORTRANCAPS
#define esmc_logsetstate_ ESMC_LOGSETSTATE
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmc_logsetstate_ esmc_logsetstate__
#elif !defined(FORTRANUNDERSCORE)
#define esmc_logsetstate_ esmc_logsetstate
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef FORTRANCAPS
#define esmc_logflush_ PMFLOGFLUSH
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmc_logflush_ pesmc_logflush__
#elif !defined(FORTRANUNDERSCORE)
#define esmc_logflush_ pesmc_logflush
#else

#define esmc_logflush_ pesmc_logflush_
#endif
#else
#ifdef FORTRANCAPS
#define esmc_logflush_ ESMC_LOGFLUSH
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmc_logflush_ esmc_logflush__
#elif !defined(FORTRANUNDERSCORE)
#define esmc_logflush_ esmc_logflush
#endif
#endif

void esmc_lognew_(ESMC_Log *log, char *logname ESMC_MIXED_LEN(namelen), ESMC_LogState *logstate,
		int *labelio, int *rc ESMC_END_LEN(namelen))
{
  char *t;

  ESMC_FIXCHAR(logname, namelen, t);
  *rc = ESMC_LogNew(log, t, *logstate, *labelio);
  ESMC_FREECHAR(logname, t);
}

void esmc_logdelete_(ESMC_Log *log, int *rc)
{
  *rc = ESMC_LogDelete(*log);
}

void esmc_logprint_(ESMC_Log *log, ESMC_LogLevel *priority, char *format, ...)
{
  va_list len_args;
  va_list args;
  
  va_start(len_args, format);
  va_start(args, format);

  /* We are somewhat restricted in who we can call by
     the stdargs */
  ESMC_LogPrintf(*log, *priority, format, len_args, args);
  
  va_end(args);
  va_end(len_args);
}

/* This can not be in the F90 module due to the stdargs */
void esmf_logprint_(ESMC_Log *log, ESMC_LogLevel *priority, char *format, ...)
{
  va_list len_args;
  va_list args;
  
  va_start(len_args, format);
  va_start(args, format);

  /* We are somewhat restricted in who we can call by
     the stdargs */
  ESMC_LogPrintf(*log, *priority, format, len_args, args);
  
  va_end(args);
  va_end(len_args);
}

/* This can not be in the F90 module due to the stdargs */
void esmf_logsetconfig_(ESMC_Log *log, int *option1, ...)
{
  va_list args;
  
  va_start(args, option1);

  /* We are somewhat restricted in who we can call by
     the stdargs */
  ESMC_LogSetConfigV(*log, *option1, args);
  
  va_end(args);
}

void esmc_logsetstate_(ESMC_Log *log, ESMC_LogState *logstate, int *rc)
{
  *rc = ESMC_LogSetState(*log, *logstate);
}

void esmc_logflush_(ESMC_Log *log, int *rc)
{
  *rc = ESMC_LogFlush(*log);
}

