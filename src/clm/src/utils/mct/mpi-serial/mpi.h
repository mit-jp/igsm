
typedef int CLMMPI_Comm;
typedef int CLMMPI_Request;

typedef int CLMMPI_Datatype;
typedef int CLMMPI_Op;

#define CLMMPI_COMM_WORLD (1)
#define CLMMPI_COMM_NULL (0)      /* handle 0 maps to NULL */


typedef int CLMMPI_Group;

/* CLMMPI_GROUP_EMPTY and CLMMPI_GROUP_NULL must not conflict with CLMMPI_GROUP_ONE */
#define CLMMPI_GROUP_EMPTY (-1)  
#define CLMMPI_GROUP_NULL  (0)


#define CLMMPI_SUCCESS   (0)

#define CLMMPI_UNDEFINED (-1)     /* value for "color" in e.g. comm_split */


/* The type's value is its size in bytes */

#define CLMMPI_BYTE           (sizeof(char))
#define CLMMPI_CHAR           (sizeof(char))
#define CLMMPI_UNSIGNED_CHAR  (sizeof(unsigned char))
#define CLMMPI_SHORT          (sizeof(short))
#define CLMMPI_UNSIGNED_SHORT (sizeof(unsigned short))
#define CLMMPI_INT            (sizeof(int))
#define CLMMPI_UNSIGNED       (sizeof(unsigned))
#define CLMMPI_LONG           (sizeof(long))
#define CLMMPI_UNSIGNED_LONG  (sizeof(unsigned long))
#define CLMMPI_FLOAT          (sizeof(float))
#define CLMMPI_DOUBLE         (sizeof(double))
#define CLMMPI_LONG_DOUBLE    (sizeof(long double))

/* types for MINLOC and MAXLOC */

#define CLMMPI_FLOAT_INT        (sizeof(struct{float a; int b;}))
#define CLMMPI_DOUBLE_INT       (sizeof(struct{double a; int b;}))
#define CLMMPI_LONG_INT         (sizeof(struct{long a; int b;}))
#define CLMMPI_2INT             (sizeof(struct{int a; int b;}))
#define CLMMPI_SHORT_INT        (sizeof (struct{short a; int b;}))
#define CLMMPI_LONG_DOUBLE_INT  (sizeof (struct{long double a; int b;}))


#define CLMMPI_ANY_TAG (-1)
#define CLMMPI_ANY_SOURCE (-1)

#define CLMMPI_REQUEST_NULL (0)

#define CLMMPI_MAX_ERROR_STRING (128)
#define CLMMPI_MAX_PROCESSOR_NAME (128)


/*
 * CLMMPI_Status
 *
 * definition must be compatible with the clmmpif.h values for
 * CLMMPI_STATUS_SIZE, CLMMPI_SOURCE, CLMMPI_TAG, and CLMMPI_ERROR.
 *
 * Note: The type used for CLMMPI_Status_int must be chosen to match
 * Fortran INTEGER.
 *
 */

typedef int CLMMPI_Status_int;

typedef struct                  /* Fortran: INTEGER status(CLMMPI_STATUS_SIZE) */
{
  CLMMPI_Status_int CLMMPI_SOURCE;    /* Fortran: status(CLMMPI_SOURCE) */
  CLMMPI_Status_int CLMMPI_TAG;       /* Fortran: status(CLMMPI_TAG) */
  CLMMPI_Status_int CLMMPI_ERROR;     /* Fortran: status(CLMMPI_ERROR) */

} CLMMPI_Status;



/* These are provided for Fortran... */

#define CLMMPI_INTEGER           CLMMPI_INT
#define CLMMPI_REAL              CLMMPI_FLOAT
#define CLMMPI_DOUBLE_PRECISION  CLMMPI_DOUBLE

#define CLMMPI_STATUS_SIZE       (sizeof(CLMMPI_Status) / sizeof(int))


/**********************************************************
 *
 * Note: if you need to regenerate the prototypes below,
 * you can use 'protify.awk' and paste the output here.
 *
 */                                                      


extern int CLMMPI_Comm_free(CLMMPI_Comm *comm);
extern int CLMMPI_Comm_size(CLMMPI_Comm comm, int *size);
extern int CLMMPI_Comm_rank(CLMMPI_Comm comm, int *rank);
extern int CLMMPI_Comm_dup(CLMMPI_Comm comm, CLMMPI_Comm *newcomm);
extern int CLMMPI_Comm_create(CLMMPI_Comm comm, CLMMPI_Group group, CLMMPI_Comm *newcomm);
extern int CLMMPI_Comm_split(CLMMPI_Comm comm, int color, int key, CLMMPI_Comm *newcomm);
extern int CLMMPI_Comm_group(CLMMPI_Comm comm, CLMMPI_Group *group);

extern int CLMMPI_Group_incl(CLMMPI_Group group, int n, int *ranks, CLMMPI_Group *newgroup);
extern int CLMMPI_Group_free(CLMMPI_Group *group);


extern int CLMMPI_Init(int *argc, char **argv[]) ;
extern int CLMMPI_Finalize(void);
extern int CLMMPI_Abort(CLMMPI_Comm comm, int errorcode);
extern int CLMMPI_Error_string(int errorcode, char *string, int *resultlen);
extern int CLMMPI_Get_processor_name(char *name, int *resultlen);
extern int CLMMPI_Initialized(int *flag);

extern int CLMMPI_Irecv(void *buf, int count, CLMMPI_Datatype datatype,
                     int source, int tag, CLMMPI_Comm comm, CLMMPI_Request *request);
extern int CLMMPI_Recv(void *buf, int count, CLMMPI_Datatype datatype, int source,
                    int tag, CLMMPI_Comm comm, CLMMPI_Status *status);
extern int CLMMPI_Isend(void *buf, int count, CLMMPI_Datatype datatype,
                     int dest, int tag, CLMMPI_Comm comm, CLMMPI_Request *request) ;
extern int CLMMPI_Send(void* buf, int count, CLMMPI_Datatype datatype,
                    int dest, int tag, CLMMPI_Comm comm);

extern int CLMMPI_Test(CLMMPI_Request *request, int *flag, CLMMPI_Status *status);
extern int CLMMPI_Wait(CLMMPI_Request *request, CLMMPI_Status *status);
extern int CLMMPI_Waitany(int count, CLMMPI_Request *array_of_requests,
                       int *index, CLMMPI_Status *status);
extern int CLMMPI_Waitall(int count, CLMMPI_Request *array_of_requests,
                       CLMMPI_Status *array_of_statuses);

extern int CLMMPI_Barrier(CLMMPI_Comm comm );
extern int CLMMPI_Bcast(void* buffer, int count, CLMMPI_Datatype datatype,
                     int root, CLMMPI_Comm comm );
extern int CLMMPI_Gather(void* sendbuf, int sendcount, CLMMPI_Datatype sendtype,
                      void* recvbuf, int recvcount, CLMMPI_Datatype recvtype,
                      int root, CLMMPI_Comm comm);
extern int CLMMPI_Gatherv(void* sendbuf, int sendcount, CLMMPI_Datatype sendtype, 
                       void* recvbuf, int *recvcounts, int *displs,
                       CLMMPI_Datatype recvtype, int root, CLMMPI_Comm comm);
extern int CLMMPI_Allgather(void* sendbuf, int sendcount, CLMMPI_Datatype sendtype,
                         void* recvbuf, int recvcount, CLMMPI_Datatype recvtype,
                         CLMMPI_Comm comm);
extern int CLMMPI_Allgatherv(void* sendbuf, int sendcount, CLMMPI_Datatype sendtype,
                          void* recvbuf, int *recvcounts, int *displs,
                          CLMMPI_Datatype recvtype, CLMMPI_Comm comm);
extern int CLMMPI_Scatterv(void* sendbuf, int *sendcounts, int *displs, 
                        CLMMPI_Datatype sendtype, void* recvbuf, int recvcount, 
                        CLMMPI_Datatype recvtype, int root, CLMMPI_Comm comm);
extern int CLMMPI_Reduce(void* sendbuf, void* recvbuf, int count, 
                      CLMMPI_Datatype datatype, CLMMPI_Op op, int root, CLMMPI_Comm comm);
extern int CLMMPI_Allreduce(void* sendbuf, void* recvbuf, int count, 
                         CLMMPI_Datatype datatype, CLMMPI_Op op, CLMMPI_Comm comm);


extern double CLMMPI_Wtime(void);

extern int CLMMPI_Alltoall(void *sendbuf, int sendcount, CLMMPI_Datatype sendtype,
                        void *recvbuf, int recvcount, CLMMPI_Datatype recvtype,
                        CLMMPI_Comm comm);

extern int CLMMPI_Alltoallv(void *sendbuf, int *sendcounts,
                         int *sdispls, CLMMPI_Datatype sendtype,
                         void *recvbuf, int *recvcounts,
                         int *rdispls, CLMMPI_Datatype recvtype,
                         CLMMPI_Comm comm) ;
