
#include "mpiP.h"



/*
 * Communicators
 *
 */



CLMMPI_Comm clmmpi_comm_new(void)
{
  CLMMPI_Comm chandle;
  Comm *cptr;
  static int num=0;

  clmmpi_alloc_handle(&chandle,(void **) &cptr);

  cptr->sendlist=AP_list_new();
  cptr->recvlist=AP_list_new();

  cptr->num=num++;

  return(chandle);
}


/*********/


FORT_NAME( clmmpi_comm_free , CLMMPI_COMM_FREE )(int *comm, int *ierror)
{
  *ierror=CLMMPI_Comm_free(comm);
}


/*
 * CLMMPI_Comm_free()
 *
 * Note: will NOT free any pending CLMMPI_Request handles
 * that are allocated...  correct user code should have
 * already done a Wait or Test to free them.
 *
 */


int CLMMPI_Comm_free(CLMMPI_Comm *comm)
{
  pList sendlist, recvlist;
  int size;
  Comm *mycomm;

  mycomm=clmmpi_handle_to_ptr(*comm);   /* (Comm *)(*comm) */

  sendlist=mycomm->sendlist;
  recvlist=mycomm->recvlist;

  size=AP_list_size(sendlist);
  if (size!=0)
    fprintf(stderr,"MPI_Comm_free: warning: %d pending send reqs\n",
	    size);
  AP_list_free(sendlist);


  size=AP_list_size(recvlist);
  if (size!=0)
    fprintf(stderr,"MPI_Comm_free: warning: %d pending receive reqs\n",
	    size);
  AP_list_free(recvlist);

  clmmpi_free_handle(*comm);            /* free(mycomm); */
  *comm=CLMMPI_COMM_NULL;

  return(CLMMPI_SUCCESS);
}


/*********/



FORT_NAME( clmmpi_comm_size , CLMMPI_COMM_SIZE )(int *comm, int *size, int *ierror)
{
  *ierror=CLMMPI_Comm_size(*comm, size);
}



int CLMMPI_Comm_size(CLMMPI_Comm comm, int *size)
{
  *size=1;

  return(CLMMPI_SUCCESS);
}


/*********/


FORT_NAME( clmmpi_comm_rank , CLMMPI_COMM_RANK )(int *comm, int *rank, int *ierror)
{
  *ierror=CLMMPI_Comm_rank( *comm, rank);
}


int CLMMPI_Comm_rank(CLMMPI_Comm comm, int *rank)
{
  *rank=0;

  return(CLMMPI_SUCCESS);
}



/*********/


FORT_NAME( clmmpi_comm_dup , CLMMPI_COMM_DUP )(int *comm, int *newcomm, int *ierror)
{

  *ierror=CLMMPI_Comm_dup( *comm, newcomm);

}


int CLMMPI_Comm_dup(CLMMPI_Comm comm, CLMMPI_Comm *newcomm)
{
  *newcomm= clmmpi_comm_new();

#ifdef INFO
  fflush(stdout);
  fprintf(stderr,"MPI_Comm_dup: new comm handle=%d\n",*newcomm);
#endif

  return(CLMMPI_SUCCESS);
}


/*********/


int FORT_NAME( clmmpi_comm_create, CLMMPI_COMM_CREATE)
     (int *comm, int *group, int *newcomm, int *ierror)
{
  *ierror=CLMMPI_Comm_create(*comm,*group,newcomm);
}



int CLMMPI_Comm_create(CLMMPI_Comm comm, CLMMPI_Group group, CLMMPI_Comm *newcomm)
{
  if (group==CLMMPI_GROUP_NULL || group==CLMMPI_GROUP_EMPTY)
    *newcomm= CLMMPI_COMM_NULL;
  else
    *newcomm=clmmpi_comm_new();

  return(CLMMPI_SUCCESS);
}



/*********/


FORT_NAME( clmmpi_comm_split, CLMMPI_COMM_SPLIT )
     (int *comm, int *color, int *key, int *newcomm, int *ierror)
{
  *ierror=CLMMPI_Comm_split(*comm,*color,*key,newcomm);

}



int CLMMPI_Comm_split(CLMMPI_Comm comm, int color, int key, CLMMPI_Comm *newcomm)
{
  if (color==CLMMPI_UNDEFINED)
    *newcomm=CLMMPI_COMM_NULL;
  else
    *newcomm= clmmpi_comm_new();

  return(CLMMPI_SUCCESS);
}


/*********/


FORT_NAME( clmmpi_comm_group, CLMMPI_COMM_GROUP )
     (int *comm, int *group, int *ierror)
{
  *ierror= CLMMPI_Comm_group(*comm, group);
}



int CLMMPI_Comm_group(CLMMPI_Comm comm, CLMMPI_Group *group)
{
  if (comm==CLMMPI_COMM_NULL)
    *group= CLMMPI_GROUP_NULL;
  else
    *group= CLMMPI_GROUP_ONE;

  return(CLMMPI_SUCCESS);
}    


