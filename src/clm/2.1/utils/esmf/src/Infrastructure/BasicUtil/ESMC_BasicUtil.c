/* $Id$ */
#include "ESMF_conf.h"

#include "ESMC_BasicUtil.h"

#include "ESMC_Error.h"

#ifdef ESMC_HAVE_MPI
#include <mpi.h> 
#endif


#ifdef ESMC_HAVE_PTHREADS
#include <pthread.h>
#endif

#ifdef ESMC_HAVE_OMP_THREADS
#include <omp.h> 
#endif


#ifdef ESMC_HAVE_PTHREADS
static pthread_mutex_t t_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif

#ifdef ESMC_HAVE_OMP_THREADS
static omp_lock_t lock;
#endif


/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_BasicUtilLockMutex"
int ESMC_BasicUtilLockMutex ()
{
#ifdef ESMC_HAVE_OMP_THREADS
  omp_set_lock (&lock);
#elif defined (ESMC_HAVE_PTHREADS)
  if (pthread_mutex_lock (&t_mutex) != 0)
    {
      ESMC_ERRA(ESMC_ERR_SYS, 0, "pthread_mutex_lock failure\n");
    }
#endif
  return ESMC_SUCCESS;
}

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_BasicUtilUnlockMutex"
int ESMC_BasicUtilUnlockMutex ()
{
#ifdef ESMC_HAVE_OMP_THREADS
  omp_unset_lock (&lock);
#elif defined(ESMC_HAVE_PTHREADS)
  if (pthread_mutex_unlock (&t_mutex) != 0)
    {
      ESMC_ERRA(ESMC_ERR_SYS, 0, "pthread_mutex_unlock failure\n");
    }
#endif
  return ESMC_SUCCESS;
}


/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_BasicUtilInit"
int ESMC_BasicUtilInit()
{
  static int initialized = 0;

  if (initialized)
    return ESMC_SUCCESS;

  initialized = 1;

#ifdef ESMC_HAVE_OMP_THREADS
  omp_init_lock(&lock);
#endif

  return ESMC_SUCCESS;
}
