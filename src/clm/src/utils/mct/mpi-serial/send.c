
#include "mpiP.h"


/*
 * SENDING
 *
 */



static int clmmpi_match_recv(void *r, void *tag)
{
  return( ((Req *)r)->tag == CLMMPI_ANY_TAG ||
	  ((Req *)r)->tag == *((int *)tag) );
}


/*
 *
 */



FORT_NAME( clmmpi_isend , CLMMPI_ISEND )(void *buf, int *count, int *datatype,
   int *dest, int *tag, int *comm, int *req, int *ierror)
{

  *ierror=CLMMPI_Isend(buf,*count,*datatype,*dest,*tag,
		    *comm, (void *)req);

}



int CLMMPI_Isend(void *buf, int count, CLMMPI_Datatype datatype,
	      int dest, int tag, CLMMPI_Comm comm, CLMMPI_Request *request) 
{
  pListitem match;
  Comm *mycomm;
  Req *rreq, *sreq;

  mycomm=clmmpi_handle_to_ptr(comm);   /* (Comm *)comm; */

#ifdef INFO
  fflush(stdout);
  fprintf(stderr,"MPI_Isend: Comm=%d  tag=%d  count=%d type=%d\n",
	 mycomm->num,tag,count,datatype);
#endif

  if (dest!=0)
    {
      fprintf(stderr,"MPI_Isend: send to %d\n",dest);
      abort();
    }

  clmmpi_alloc_handle(request,(void **) &sreq);


  if ( match=AP_list_search_func(mycomm->recvlist,clmmpi_match_recv,&tag) )
    {
      rreq=(Req *)AP_listitem_data(match);
      AP_list_delete_item(mycomm->recvlist,match);

      memcpy(rreq->buf,buf,count * datatype);
      rreq->complete=1;
      rreq->tag=tag;                    /* in case rreq->tag was CLMMPI_ANY_TAG */

      sreq->complete=1;

#ifdef DEBUG
      printf("Completion(send) value=%d tag=%d\n",
	     *((int *)buf),rreq->tag);
#endif

      return(CLMMPI_SUCCESS);
    }

  sreq->buf=buf;
  sreq->tag=tag;
  sreq->complete=0;
  sreq->listitem=AP_list_append(mycomm->sendlist,sreq);

#ifdef INFO
  print_list(mycomm->sendlist,"sendlist for comm ",mycomm->num);
#endif

  return(CLMMPI_SUCCESS);
}


/*********/


FORT_NAME(clmmpi_send, CLMMPI_SEND) ( void *buf, int *count, int *datatype,
 		                int *dest, int *tag, int *comm, int *ierror)
{
  *ierror=CLMMPI_Send(buf, *count, *datatype, *dest, *tag, *comm);
}



int CLMMPI_Send(void* buf, int count, CLMMPI_Datatype datatype,
	     int dest, int tag, CLMMPI_Comm comm)
{
  CLMMPI_Request request;
  CLMMPI_Status status;
  int flag;

#ifdef INFO
  fflush(stdout);
  fprintf(stderr,"MPI_Send: ");
#endif

  CLMMPI_Isend(buf,count,datatype,dest,tag,comm,&request);
  CLMMPI_Wait(&request,&status);


  return(CLMMPI_SUCCESS);
}



