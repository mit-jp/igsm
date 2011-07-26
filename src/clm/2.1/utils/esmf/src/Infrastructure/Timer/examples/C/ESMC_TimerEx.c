/* $Id$ */

#include "ESMC.h"



int do_work(ESMC_Log log)
{
  int i, j, k;
  double x, y;


  x = y = 0.25;

  for (i = 0; i < 5*1024; i++)
    {
      x = x * x + y;
      ESMC_LogPrint(log, ESMC_LOGLEVEL_INFO, "x=%f\n", x);
    }

  return 0;
}

int main(int argc, char *argv[])
{
  ESMC_App app;
  ESMC_Log log;
  int rc;

  ESMC_AppNew(&app);
  ESMC_LogNew(&log, "timelog", ESMC_LOGSTATE_VERBOSE, 0);

  rc = ESMC_TimerInit(ESMC_ALL_TIMERS,
		     ESMC_USRSYS, 1,
		     ESMC_WALL, 1,
		     0
		     );

  ESMC_ERROR_TEST((rc == ESMC_SUCCESS), "ESMC_TimerInit: Initialize timers");

  rc = ESMC_TimerStart("DavidN");
  ESMC_ERROR_TEST((rc == ESMC_SUCCESS), "ESMC_TimerStart: start a timer");
  
  do_work(log);
  rc = ESMC_TimerStop("DavidN");
  ESMC_ERROR_TEST((rc == ESMC_SUCCESS), "ESMC_TimerStop: stop a timer");
  
  rc = ESMC_TimerPrint(ESMC_ALL_TIMERS, 0);
  ESMC_ERROR_TEST((rc == ESMC_SUCCESS), "ESMC_TimerPrint: Print a timer");
  
  ESMC_LogDelete(log);
  ESMC_AppDelete(app);

  return 0;
}
