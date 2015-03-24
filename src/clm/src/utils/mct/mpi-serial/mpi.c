

#include "mpiP.h"


/****************************************************************************/

static int initialized=0;


/****************************************************************************/


/*
 * INIT/FINALIZE
 *
 */



FORT_NAME( clmmpi_init_fort , CLMMPI_INIT_FORT)
                          (int *f_CLMMPI_COMM_WORLD,
                           int *f_CLMMPI_ANY_SOURCE, int *f_CLMMPI_ANY_TAG,
                           int *f_CLMMPI_COMM_NULL, int *f_CLMMPI_REQUEST_NULL,
			   int *f_CLMMPI_GROUP_NULL, int *f_CLMMPI_GROUP_EMPTY,
			   int *f_CLMMPI_UNDEFINED,
                           int *f_CLMMPI_MAX_ERROR_STRING, 
                           int *f_CLMMPI_MAX_PROCESSOR_NAME, 
                           int *f_CLMMPI_STATUS_SIZE, 
                           int *f_CLMMPI_SOURCE, int *f_CLMMPI_TAG, int *f_CLMMPI_ERROR,
			   int *f_status,
			   int *fsource, int *ftag, int *ferror,
                           int *f_CLMMPI_INTEGER, void *fint1, void *fint2,
                           int *f_CLMMPI_LOGICAL, void *flog1, void *flog2,
                           int *f_CLMMPI_REAL, void *freal1, void *freal2,
                           int *f_CLMMPI_DOUBLE_PRECISION,
			   void *fdub1, void *fdub2,
			   int *f_CLMMPI_COMPLEX, void *fcomp1, void *fcomp2,
                           int *ierror)
{
  int err;
  int size;
  int offset;

  *ierror=CLMMPI_Init(NULL,NULL);

  err=0;

  /*
   * These 3 macros compare things from clmmpif.h (as passed in by the f_
   * arguments) to the values in C (from #including clmmpi.h).
   *
   * Unfortunately, this kind of thing is done most easily in a nasty
   * looking macto.
   *
   */


  /*
   * verify_eq
   *   compare value of constants in C and fortran
   *   i.e. compare *f_<name> to <name>
   */

#define verify_eq(name)  \
  if (*f_##name != name) \
    { fprintf(stderr,"mpi-serial: clmmpi_init_fort: %s not consistant " \
                     "between clmmpif.h (%d) and clmmpi.h (%d)\n",\
                     #name,*f_##name,name); \
      err=1; }

#define verify_eq_warn(name)  \
  if (*f_##name != name) \
    { fprintf(stderr,"mpi-serial: clmmpi_init_fort: warning: %s not consistant " \
                     "between clmmpif.h (%d) and clmmpi.h (%d)\n",\
                     #name,*f_##name,name); \
    }


  /*
   * verify_size
   *   verify that the type name in fortran has the correct
   *   value (i.e. the size of that data type).
   *   Determine size by subtracting the pointer values of two
   *   consecutive array locations.
   */

#define verify_size(name,p1,p2) \
  if ( (size=((char *)(p2) - (char *)(p1))) != *f_##name ) \
    { fprintf(stderr,"mpi-serial: clmmpi_init_fort: clmmpif.h %s (%d) " \
                     "does not match actual fortran size (%d)\n", \
                     #name,*f_##name,size); \
      err=1; }

  /*
   * verify_field
   *   check the struct member offsets for CLMMPI_Status vs. the
   *   fortan integer array offsets.  E.g. the location of
   *   status->MPI_SOURCE should be the same as STATUS(CLMMPI_SOURCE)
   */

#define verify_field(name) \
  { offset= (char *)&((CLMMPI_Status *)f_status)->name - (char *)f_status; \
    if ( offset != (*f_##name-1)*sizeof(int) ) \
    { fprintf(stderr,"mpi-serial: clmmpi_init_fort: clmmpif.h %s (%d) (%d bytes) " \
                     "is inconsistant w/offset in CLMMPI_Status (%d bytes)\n", \
                    #name,*f_##name,(*f_##name-1)*sizeof(int),offset); \
      err=1; }}



  verify_eq(CLMMPI_COMM_WORLD);
  verify_eq(CLMMPI_ANY_SOURCE);
  verify_eq(CLMMPI_ANY_TAG);
  verify_eq(CLMMPI_COMM_NULL);
  verify_eq(CLMMPI_REQUEST_NULL);
  verify_eq(CLMMPI_GROUP_NULL);
  verify_eq(CLMMPI_GROUP_EMPTY);
  verify_eq(CLMMPI_UNDEFINED);
  verify_eq(CLMMPI_MAX_ERROR_STRING);
  verify_eq(CLMMPI_MAX_PROCESSOR_NAME);

  verify_eq(CLMMPI_STATUS_SIZE);
  verify_field(CLMMPI_SOURCE);
  verify_field(CLMMPI_TAG);
  verify_field(CLMMPI_ERROR);

  verify_eq(CLMMPI_INTEGER);
  verify_size(CLMMPI_INTEGER,fint1,fint2);

  verify_size(CLMMPI_LOGICAL,flog1,flog2);

  verify_eq_warn(CLMMPI_REAL);
  verify_size(CLMMPI_REAL,freal1,freal2);

  verify_eq(CLMMPI_DOUBLE_PRECISION);
  verify_size(CLMMPI_DOUBLE_PRECISION,fdub1,fdub2);

  verify_size(CLMMPI_COMPLEX,fcomp1,fcomp2);

  if (err)
    abort();
}



int CLMMPI_Init(int *argc, char **argv[]) 
{
  CLMMPI_Comm my_comm_world;

  my_comm_world=clmmpi_comm_new();

  if (my_comm_world != CLMMPI_COMM_WORLD)
    {
      fprintf(stderr,"MPI_Init: conflicting CLMMPI_COMM_WORLD\n");
      abort();
    }

  initialized=1;
  return(CLMMPI_SUCCESS);
}


/*********/


FORT_NAME( clmmpi_finalize, CLMMPI_FINALIZE )(int *ierror)
{
  *ierror=CLMMPI_Finalize();
}


/*
 * CLMMPI_Finalize()
 *
 * this library doesn't support re-initializing CLMMPI, so
 * the finalize will just leave everythign as it is...
 *
 */


int CLMMPI_Finalize(void)
{
  initialized=0;

  clmmpi_destroy_handles();

  return(CLMMPI_SUCCESS);
}


/*********/


FORT_NAME( clmmpi_abort , CLMMPI_ABORT )(int *comm, int *errorcode, int *ierror)
{
  *ierror=CLMMPI_Abort( *comm, *errorcode);
}



int CLMMPI_Abort(CLMMPI_Comm comm, int errorcode)
{
  fprintf(stderr,"MPI_Abort: error code = %d\n",errorcode);
  exit(errorcode);
}


/*********/



FORT_NAME( clmmpi_error_string , CLMMPI_ERROR_STRING)
                             (int *errorcode, char *string,
			      int *resultlen, int *ierror)
{
  *ierror=CLMMPI_Error_string(*errorcode, string, resultlen);
}


int CLMMPI_Error_string(int errorcode, char *string, int *resultlen)
{
  sprintf(string,"MPI Error: code %d\n",errorcode);
  *resultlen=strlen(string);

  return(CLMMPI_SUCCESS);
}


/*********/


FORT_NAME( clmmpi_get_processor_name , CLMMPI_GET_PROCESSOR_NAME )
                          (char *name, int *resultlen, int *ierror)
{
  *ierror=CLMMPI_Get_processor_name(name,resultlen);
}


int CLMMPI_Get_processor_name(char *name, int *resultlen)
{
  int ret;

  ret=gethostname(name,CLMMPI_MAX_PROCESSOR_NAME);

  if (ret!=0)
    strncpy(name,"unknown host name",CLMMPI_MAX_PROCESSOR_NAME);


  name[CLMMPI_MAX_PROCESSOR_NAME-1]='\0';  /* make sure NULL terminated */
  *resultlen=strlen(name);

  return(CLMMPI_SUCCESS);
}


/*********/


FORT_NAME( clmmpi_initialized , CLMMPI_INITIALIZED )(int *flag, int *ierror)
{
  *ierror=CLMMPI_Initialized(flag);
}


int CLMMPI_Initialized(int *flag)
{
  *flag= initialized;

  return(CLMMPI_SUCCESS);
}



