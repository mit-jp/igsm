/* $Id$ */
#include <stdlib.h>
#include <stdio.h>

#include "ESMC.h"

#include "ESMF_conf.h"

#include "ESMC_App.h"
#include "ESMC_BasicUtil.h"
#include "ESMC_Machine.h"
#include "ESMC_Log.h"

#ifdef ESMC_HAVE_MPI
#include <mpi.h>
#endif

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_AppNew"
int ESMC_AppNew(ESMC_App *app)
{
  int argc;
  int ret;

  *app = (ESMC_App) malloc(sizeof(ESMC_AppClass));

  if ((ret = ESMC_AppConstruct(*app)) != ESMC_SUCCESS)
    {
      free(*app);
      return ret;
    }
  
  return ESMC_SUCCESS;
}

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_AppConstruct"
int ESMC_AppConstruct(ESMC_App app)
{
  static int initialized = 0;
  int ret;
  int mpi_initialized;
  int argc;
  char *argv[] = {"program", NULL};

  if (initialized)
    {
      ESMC_ERRA(ESMC_ERR_BUSY, 0, "There can be only one App.");
    }
  else
    initialized = 1;

/* MPI Init must be called before logger, since logger makes MPI
   calls */
#ifdef ESMC_HAVE_MPI

  MPI_Initialized(&mpi_initialized);
  
  if (!mpi_initialized)
    MPI_Init(&argc, (char ***) &argv);
#endif

  /* Initialize BasicUtils */
  if ((ret = ESMC_BasicUtilInit()) != ESMC_SUCCESS)
    return ret;
  
  
  /* Initialize Machine Model */
  if ((ret = ESMC_MachineNew(&app->machine)) != ESMC_SUCCESS)
    return ret;
  
  if ((ret = ESMC_LogNew(&app->logSTD, "stdlog", ESMC_LOGSTATE_VERBOSE, 1))
      != ESMC_SUCCESS)
    {
      ESMC_MachineDelete(app->machine);
      return ret;
    }

  /* Set up the Timer's stdlog output */
  if ((ret = ESMC_TimerSetSTDLog(app->logSTD)) != ESMC_SUCCESS)
    {
      ESMC_MachineDelete(app->machine);
      ESMC_LogDelete(app->logSTD);
      return ret;
    }

  return ESMC_SUCCESS;
}

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_AppDelete"
int ESMC_AppDelete(ESMC_App app)
{
  int ret;

  if ((ret = ESMC_AppDestruct(app)) != ESMC_SUCCESS)
    {
      return ret;
    }
  free(app);

  return ESMC_SUCCESS;
}

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_AppDestruct"
int ESMC_AppDestruct(ESMC_App app)
{
#ifdef ESMC_HAVE_MPI
  MPI_Finalize();
#endif

  ESMC_MachineDelete(app->machine);

  ESMC_LogDelete(app->logSTD);
  
  return ESMC_SUCCESS;
}
