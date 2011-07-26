/* $Id$ */

#include "ESMC.h"
#ifdef  ESMC_HAVE_OMP_THREADS
#include <omp.h>
#endif


int main(int argc, char *argv[])
{
  ESMC_App app;
  ESMC_Log log, testlog, opttest;
  int i;
  
  ESMC_AppNew(&app);
  
  ESMC_LogNew(&log, "logfile", ESMC_LOGSTATE_VERBOSE, 1);
  ESMC_LogNew(&testlog, "testlog", ESMC_LOGSTATE_VERBOSE, 0);
  ESMC_LogNew(&opttest, "optlog", ESMC_LOGSTATE_VERBOSE, 0);

  ESMC_LogPrint(testlog, ESMC_LOGLEVEL_INFO, "AAA: Hello TestLogger\n");

#pragma omp parallel for private(i)
  for (i = 0; i < 10; i++)
    {
      ESMC_LogPrint(testlog, ESMC_LOGLEVEL_INFO, "In OMP for Iteration\n");
    }

  /* Test option setting */
  ESMC_LogPrint(opttest, ESMC_LOGLEVEL_INFO, "This should display defaults\n");

  ESMC_LogSetConfig(opttest, 
	ESMC_LOGCONFIG_PRINTTIME, 0, 
        NULL);
  ESMC_LogPrint(opttest, ESMC_LOGLEVEL_INFO, "This should not display time\n");

  ESMC_LogSetConfig(opttest, 
	ESMC_LOGCONFIG_PRINTTID, 0, 
        NULL);
  ESMC_LogPrint(opttest, ESMC_LOGLEVEL_INFO, "This should not display time or tid\n");

  ESMC_LogSetConfig(opttest, 
	ESMC_LOGCONFIG_PRINTPID, 0, 
        NULL);
  ESMC_LogPrint(opttest, ESMC_LOGLEVEL_INFO, "This should not display time or tid or pid\n");

  ESMC_LogSetConfig(opttest,
	ESMC_LOGCONFIG_PRINTTIME, 1,
	ESMC_LOGCONFIG_PRINTTID, 1,
	ESMC_LOGCONFIG_PRINTPID, 1,
	NULL);
  ESMC_LogPrint(opttest, ESMC_LOGLEVEL_INFO, "This should display defaults\n");

  
  ESMC_LogDelete(opttest);
  ESMC_LogDelete(testlog);
  ESMC_LogDelete(log);
  ESMC_AppDelete(app);

  return 0;
}
