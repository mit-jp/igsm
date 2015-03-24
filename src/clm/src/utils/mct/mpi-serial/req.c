#include "mpiP.h"


/*
 * COMPLETION
 */



FORT_NAME( clmmpi_test , CLMMPI_TEST)(int *request, int *flag, int *status,
				int *ierror)
{
  *ierror=CLMMPI_Test( (void *)request ,flag,(CLMMPI_Status *)status);
}



int CLMMPI_Test(CLMMPI_Request *request, int *flag, CLMMPI_Status *status)
{
  Req *req;

  if (*request==CLMMPI_REQUEST_NULL)
    {
      status->CLMMPI_TAG= CLMMPI_ANY_TAG;
      status->CLMMPI_SOURCE= CLMMPI_ANY_SOURCE;
      *flag=1;
      return(CLMMPI_SUCCESS);
    }


  req=clmmpi_handle_to_ptr(*request);

  *flag=req->complete;

  if (*flag)
    {
      status->CLMMPI_SOURCE= 0;
      status->CLMMPI_TAG= req->tag;

      clmmpi_free_handle(*request);
      *request=CLMMPI_REQUEST_NULL;
    }

  return(CLMMPI_SUCCESS);
}



FORT_NAME( clmmpi_wait , CLMMPI_WAIT )(int *request, int *status, int *ierror)
{
  *ierror=CLMMPI_Wait( (void *)request, (CLMMPI_Status *)status );
}



int CLMMPI_Wait(CLMMPI_Request *request, CLMMPI_Status *status)
{
  int flag;

  CLMMPI_Test(request,&flag,status);

  if (!flag)
    {
      fprintf(stderr,"MPI_Wait: request not complete, deadlock\n");
      abort();
    }

  return(CLMMPI_SUCCESS);
}


/*********/


FORT_NAME( clmmpi_waitany , CLMMPI_WAITANY )(int *count, int *requests,
				       int *index, int *status, int *ierror)
{

  *ierror=CLMMPI_Waitany(*count, (void *)requests,index,(CLMMPI_Status *)status);
}



int CLMMPI_Waitany(int count, CLMMPI_Request *array_of_requests,
		int *index, CLMMPI_Status *status)
{
  int i;
  int flag;

  for (i=0; i<count; i++)
    {
      CLMMPI_Test(&array_of_requests[i],&flag,status);
      
      if (flag)
	return(CLMMPI_SUCCESS);
    }

  /* none are completed */

  fprintf(stderr,"MPI_Waitany: no requests complete, deadlock\n");
  abort();
}


/*********/


FORT_NAME( clmmpi_waitall , CLMMPI_WAITALL )(int *count, int *array_of_requests,
				       int *array_of_statuses, int *ierror)
{
  *ierror=CLMMPI_Waitall(*count, (void *)array_of_requests,
		      (CLMMPI_Status *)array_of_statuses);

}



int CLMMPI_Waitall(int count, CLMMPI_Request *array_of_requests,
		CLMMPI_Status *array_of_statuses)
{
  int i;
  int flag;

  for (i=0; i<count; i++)
    {
      CLMMPI_Test(&array_of_requests[i],&flag,&array_of_statuses[i]);

      if (!flag)
	{
	  fprintf(stderr,"MPI_Waitall: request not complete, deadlock\n");
	  abort();
	}
    }

  return(CLMMPI_SUCCESS);
}



