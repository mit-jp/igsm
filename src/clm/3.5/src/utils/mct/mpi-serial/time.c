
#include <sys/time.h>
#include <time.h>


#include "mpiP.h"


double CLMMPI_Wtime(void);



double FORT_NAME( clmmpi_wtime, CLMMPI_WTIME )(void)
{
  return(CLMMPI_Wtime());
}



double CLMMPI_Wtime(void)
{
  struct timeval tv;

  if (gettimeofday(&tv,0))
    {
      fprintf(stderr,"MPI_Wtime: error calling gettimeofday()\n");
      abort();
    }


  return((double)(tv.tv_sec) + (double)(tv.tv_usec)/1e6) ;
}



