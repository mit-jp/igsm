
#include "mpiP.h"



/*
 * COLLECTIVE
 */


FORT_NAME( clmmpi_barrier , CLMMPI_BARRIER )(int *comm, int *ierror)
{
  *ierror=CLMMPI_Barrier( *comm );
}


int CLMMPI_Barrier(CLMMPI_Comm comm )
{
  return(CLMMPI_SUCCESS);
}


/*********/


FORT_NAME( clmmpi_bcast , CLMMPI_BCAST )(void *buffer, int *count, int *datatype,
				   int *root, int *comm, int *ierror )
{
  *ierror=CLMMPI_Bcast(buffer, *count, *datatype, *root, *comm);
}



int CLMMPI_Bcast(void* buffer, int count, CLMMPI_Datatype datatype,
	      int root, CLMMPI_Comm comm )
{
  if (root!=0)
    {
      fprintf(stderr,"MPI_Bcast: bad root = %d\n",root);
      abort();
    }


  return(CLMMPI_SUCCESS);
}


/*********/


FORT_NAME( clmmpi_gather , CLMMPI_GATHER )
                       (void *sendbuf, int *sendcount, int *sendtype,
			void *recvbuf, int *recvcount, int *recvtype,
			int *root, int *comm, int *ierror)
{
  *ierror=CLMMPI_Gather( sendbuf, *sendcount, *sendtype,
		      recvbuf, *recvcount, *recvtype,
		      *root, *comm);
}


int CLMMPI_Gather(void* sendbuf, int sendcount, CLMMPI_Datatype sendtype,
	       void* recvbuf, int recvcount, CLMMPI_Datatype recvtype,
	       int root, CLMMPI_Comm comm)
{
  if (root!=0)
    {
      fprintf(stderr,"MPI_Gather: bad root = %d\n",root);
      abort();
    }

  memcpy(recvbuf,sendbuf,sendcount*sendtype);

  return(CLMMPI_SUCCESS);
}

/*********/



FORT_NAME( clmmpi_gatherv , CLMMPI_GATHERV )
                        ( void *sendbuf, int *sendcount, int *sendtype,
			  void *recvbuf, int *recvcounts, int *displs,
			  int *recvtype, int *root, int *comm, int *ierror)
{
  *ierror=CLMMPI_Gatherv( sendbuf, *sendcount, *sendtype,
		       recvbuf, recvcounts, displs,
		       *recvtype, *root, *comm);
}


int CLMMPI_Gatherv(void* sendbuf, int sendcount, CLMMPI_Datatype sendtype, 
		void* recvbuf, int *recvcounts, int *displs,
		CLMMPI_Datatype recvtype, int root, CLMMPI_Comm comm)
{
  int offset;

  if (root!=0)
    {
      fprintf(stderr,"MPI_Gatherv: bad root = %d\n",root);
      abort();
    }

  offset=displs[0]*recvtype;
  memcpy( (char *)recvbuf+offset, sendbuf, recvcounts[0] * recvtype);

  return(CLMMPI_SUCCESS);
}



/*********/


FORT_NAME( clmmpi_allgather , CLMMPI_ALLGATHER )
                          ( void *sendbuf, int *sendcount, int *sendtype,
			    void *recvbuf, int *recvcount, int *recvtype,
			    int *comm, int *ierror)
{
  *ierror=CLMMPI_Allgather( sendbuf, *sendcount, *sendtype,
			 recvbuf, *recvcount, *recvtype,
			 *comm );
}


int CLMMPI_Allgather(void* sendbuf, int sendcount, CLMMPI_Datatype sendtype,
		  void* recvbuf, int recvcount, CLMMPI_Datatype recvtype,
		  CLMMPI_Comm comm)
{

  memcpy(recvbuf,sendbuf,sendcount * sendtype);

  return(CLMMPI_SUCCESS);

}


/*********/


FORT_NAME( clmmpi_allgatherv , CLMMPI_ALLGATHERV )
                          ( void *sendbuf, int *sendcount, int *sendtype,
			    void *recvbuf, int *recvcounts, int *displs,
                            int *recvtype, int *comm, int *ierror)
{
  *ierror=CLMMPI_Allgatherv( sendbuf, *sendcount, *sendtype,
			  recvbuf, recvcounts, displs,
                          *recvtype, *comm );
}


int CLMMPI_Allgatherv(void* sendbuf, int sendcount, CLMMPI_Datatype sendtype,
		   void* recvbuf, int *recvcounts, int *displs,
                   CLMMPI_Datatype recvtype, CLMMPI_Comm comm)
{
  int offset;

  offset=displs[0]*recvtype;
  memcpy( (char *)recvbuf+offset, sendbuf, recvcounts[0] * recvtype);

  return(CLMMPI_SUCCESS);
}


/*********/


FORT_NAME( clmmpi_scatterv , CLMMPI_SCATTERV )
                         ( void *sendbuf, int *sendcounts, int *displs,
			   int *sendtype, void *recvbuf, int *recvcount,
			   int *recvtype, int *root, int *comm, int *ierror)
{
  *ierror=CLMMPI_Scatterv(sendbuf, sendcounts, displs,
		       *sendtype, recvbuf, *recvcount,
		       *recvtype, *root, *comm);
}
		       
		       

int CLMMPI_Scatterv(void* sendbuf, int *sendcounts, int *displs, 
		 CLMMPI_Datatype sendtype, void* recvbuf, int recvcount, 
		 CLMMPI_Datatype recvtype, int root, CLMMPI_Comm comm)
{
  int offset;

  if (root!=0)
    {
      fprintf(stderr,"MPI_Scatterv: bad root = %d\n",root);
      abort();
    }

  offset=displs[0]*sendtype;
  memcpy(recvbuf,(char *)sendbuf+offset,sendcounts[0] * sendtype);
  
  return(CLMMPI_SUCCESS);
}



/*********/


FORT_NAME( clmmpi_reduce , CLMMPI_REDUCE )
                       ( void *sendbuf, void *recvbuf, int *count,
			 int *datatype, int *op, int *root, int *comm,
			 int *ierror)
{
  *ierror=CLMMPI_Reduce(sendbuf, recvbuf, *count,
		     *datatype, *op, *root, *comm);
}



int CLMMPI_Reduce(void* sendbuf, void* recvbuf, int count, 
	       CLMMPI_Datatype datatype, CLMMPI_Op op, int root, CLMMPI_Comm comm)

{
  if (root!=0)
    {
      fprintf(stderr,"MPI_Reduce: bad root = %d\n",root);
      abort();
    }

  memcpy(recvbuf,sendbuf,count * datatype);

  return(CLMMPI_SUCCESS);
}


/*********/


FORT_NAME( clmmpi_allreduce , CLMMPI_ALLREDUCE )
                          ( void *sendbuf, void *recvbuf, int *count,
			    int *datatype, int *op, int *comm, int *ierror)
{
  *ierror=CLMMPI_Allreduce(sendbuf, recvbuf, *count,
			*datatype, *op, *comm);

}


int CLMMPI_Allreduce(void* sendbuf, void* recvbuf, int count, 
		  CLMMPI_Datatype datatype, CLMMPI_Op op, CLMMPI_Comm comm)
{

  memcpy(recvbuf,sendbuf,count * datatype);

  return(CLMMPI_SUCCESS);

}


/*********/


FORT_NAME( clmmpi_alltoall , CLMMPI_ALLTOALL )
                        ( void *sendbuf, int *sendcount, int *sendtype,
			  void *recvbuf, int *recvcount, int *recvtype,
                          int *comm, int *ierror )
{
  *ierror=CLMMPI_Alltoall(sendbuf, *sendcount, *sendtype,
		       recvbuf, *recvcount, *recvtype,
		       *comm);
}


int CLMMPI_Alltoall(void *sendbuf, int sendcount, CLMMPI_Datatype sendtype,
		 void *recvbuf, int recvcount, CLMMPI_Datatype recvtype,
		 CLMMPI_Comm comm)
{

  memcpy(recvbuf,sendbuf,sendcount * sendtype);

  return(CLMMPI_SUCCESS);
}


/*********/


FORT_NAME( clmmpi_alltoallv , CLMMPI_ALLTOALLV )
           ( void *sendbuf, int *sendcounts, int *sdispls, int *sendtype,
	     void *recvbuf, int *recvcounts, int *rdispls, int *recvtype,
             int *comm, int *ierror )
{

  *ierror=CLMMPI_Alltoallv(sendbuf, sendcounts, sdispls, *sendtype,
			recvbuf, recvcounts, rdispls, *recvtype,
			*comm);

}

int CLMMPI_Alltoallv(void *sendbuf, int *sendcounts,
		  int *sdispls, CLMMPI_Datatype sendtype,
                  void *recvbuf, int *recvcounts,
		  int *rdispls, CLMMPI_Datatype recvtype,
                  CLMMPI_Comm comm) 

{
  int send_offset;
  int recv_offset;

  send_offset=sdispls[0]*sendtype;
  recv_offset=rdispls[0]*recvtype;


  memcpy( (char *)recvbuf+recv_offset, (char *)sendbuf+send_offset,
	  sendcounts[0] * sendtype);


  return(CLMMPI_SUCCESS);
}
