/* $Id$ */

#include "ESMC.h"
#include "ESMC_Timer.h"
#include "ESMC_Log.h"

#include "ESMC_Fortran.h"
#include <malloc.h>

#ifdef ESMC_HAVE_FORTRAN_UNDERSCORE
#define FORTRANUNDERSCORE
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef FORTRANCAPS
#define esmc_timerinit_ PESMC_TIMERINIT
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmc_timerinit_ pesmc_timerinit__
#elif !defined(FORTRANUNDERSCORE)
#define esmc_timerinit_ pesmc_timerinit
#else

#define esmc_timerinit_ pesmc_timerinit_
#endif
#else
#ifdef FORTRANCAPS
#define esmc_timerinit_ ESMC_TIMERINIT
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmc_timerinit_ esmc_timerinit__
#elif !defined(FORTRANUNDERSCORE)
#define esmc_timerinit_ esmc_timerinit
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef FORTRANCAPS
#define esmf_timerinit_ PESMC_TIMERINIT
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmf_timerinit_ pesmf_timerinit__
#elif !defined(FORTRANUNDERSCORE)
#define esmf_timerinit_ pesmf_timerinit
#else

#define esmf_timerinit_ pesmf_timerinit_
#endif
#else
#ifdef FORTRANCAPS
#define esmf_timerinit_ ESMF_TIMERINIT
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmf_timerinit_ esmf_timerinit__
#elif !defined(FORTRANUNDERSCORE)
#define esmf_timerinit_ esmf_timerinit
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef FORTRANCAPS
#define esmc_timerstart_ PESMF_TIMERSTART
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmc_timerstart_ pesmc_timerstart__
#elif !defined(FORTRANUNDERSCORE)
#define esmc_timerstart_ pesmc_timerstart
#else

#define esmc_timerstart_ pesmc_timerstart_
#endif
#else
#ifdef FORTRANCAPS
#define esmc_timerstart_ ESMC_TIMERSTART
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmc_timerstart_ esmc_timerstart__
#elif !defined(FORTRANUNDERSCORE)
#define esmc_timerstart_ esmc_timerstart
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef FORTRANCAPS
#define esmc_timerstop_ PESMC_TIMERSTOP
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmc_timerstop_ pesmc_timerstop__
#elif !defined(FORTRANUNDERSCORE)
#define esmc_timerstop_ pesmc_timerstop
#else

#define esmc_timerstop_ pesmc_timerstop_
#endif
#else
#ifdef FORTRANCAPS
#define esmc_timerstop_ ESMC_TIMERSTOP
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmc_timerstop_ esmc_timerstop__
#elif !defined(FORTRANUNDERSCORE)
#define esmc_timerstop_ esmc_timerstop
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef FORTRANCAPS
#define esmc_timerstamp_ PESMC_TIMERSTAMP
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmc_timerstamp_ pesmc_timerstamp__
#elif !defined(FORTRANUNDERSCORE)
#define esmc_timerstamp_ pesmc_timerstamp
#else

#define esmc_timerstamp_ pesmc_timerstamp_
#endif
#else
#ifdef FORTRANCAPS
#define esmc_timerstamp_ ESMC_TIMERSTAMP
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmc_timerstamp_ esmc_timerstamp__
#elif !defined(FORTRANUNDERSCORE)
#define esmc_timerstamp_ esmc_timerstamp
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef FORTRANCAPS
#define esmc_timerprint_ PESMC_TIMERPRINT
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmc_timerprint_ pesmc_timerprint__
#elif !defined(FORTRANUNDERSCORE)
#define esmc_timerprint_ pesmc_timerprint
#else

#define esmc_timerprint_ pesmc_timerprint_
#endif
#else
#ifdef FORTRANCAPS
#define esmc_timerprint_ ESMC_TIMERPRINT
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define esmc_timerprint_ esmc_timerprint__
#elif !defined(FORTRANUNDERSCORE)
#define esmc_timerprint_ esmc_timerprint
#endif
#endif


#include <stdarg.h>


void esmc_timerinit_(char *name, int *option1, ...)
{
  /* This function uses var args, so it will have to deal with
     the string length as it picks off args */
  va_list args;

  va_start(args, option1);

  ESMC_TimerInitV(name, *option1, args);
  
  va_end(args);
}

void esmf_timerinit_(char *name, int *option1, ...)
{
  /* This function uses var args, so it will have to deal with
     the string length as it picks off args */
  va_list args;

  va_start(args, option1);

  ESMC_TimerInitV(name, *option1, args);
  
  va_end(args);
}

void esmc_timerstart_(char *name ESMC_MIXED_LEN(len), int *rc ESMC_END_LEN(len))
{
  char *t;

  ESMC_FIXCHAR(name, len, t);
  *rc = ESMC_TimerStart(t);

  ESMC_FREECHAR(name, t);
}

void esmc_timerstop_(char *name ESMC_MIXED_LEN(len), int *rc ESMC_END_LEN(len))
{
  char *t;

  ESMC_FIXCHAR(name, len, t);
  *rc = ESMC_TimerStop(t);
  ESMC_FREECHAR(name, t);
}

void esmc_timerstamp_(double *wall, double *user, double *sys, int *rc)
{
  *rc = ESMC_TimerStamp(wall, user, sys);
}

void esmc_timerprint_(char *name ESMC_MIXED_LEN(len), ESMC_Log *log, int *rc ESMC_END_LEN(len))
{
  char *t;

  ESMC_FIXCHAR(name, len, t);
  *rc = ESMC_TimerPrint(t, *log);
  ESMC_FREECHAR(name, t);
}

