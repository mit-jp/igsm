
/*
 * Private .h file for CLMMPI
 */


#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#include "listops.h"
#include "mpi.h"


/*
 * Fortran name mangling
 *
 * the configure.ac specifies these
 *
 * cpp does not have the ability to change the case
 * of the argument, so the invocation of the macro
 * has to be give both e.g. FORT_NAME(hello,HELLO)
 * and maps to "hello_", "hello", and "HELLO" repectively.
 *
 * IMPORTANT NOTE:
 * In the case of FORTRAN_GNUF2C (e.g. g95), the rule is this:
 *    name does not contain an underscore -> append *one* underscore
 *    name contains an underscore -> append *two* underscore
 * Since all the clmmpi-serial names exported to fortran start with "mpi_",
 * we only support the latter.
 *
 */


#if   defined(FORTRAN_UNDERSCORE_)
#define FORT_NAME(lower,upper) lower##_
#elif   defined(FORTRANUNDERSCORE)
#define FORT_NAME(lower,upper) lower##_
#elif   defined(FORTRAN_GNUF2C)
#define FORT_NAME(lower,upper) lower##__
#elif defined(FORTRAN_SAME)
#define FORT_NAME(lower,upper) lower
#elif defined(FORTRAN_CAPS_)
#define FORT_NAME(lower,upper) upper
#else
#define FORT_NAME(lower,upper) FORTRAN_MANGLE_ERROR
#error "Unrecognized Fortran-mangle type"
#endif


/*
 * CLMMPI_GROUP_ONE must not conflict with CLMMPI_GROUP_NULL or
 * CLMMPI_GROUP_EMPTY
 */

#define CLMMPI_GROUP_ONE  (1)


/****************************************************************************/


typedef struct
{
  pList sendlist;
  pList recvlist;

  int num;

} Comm;



typedef struct
{
  pListitem listitem;        /* to allow Req to be removed from list */

  int *buf;
  int tag;
  int complete;

} Req;




/****************************************************************************/


extern void *clmmpi_malloc(int size);
extern void clmmpi_free(void *ptr);

extern CLMMPI_Comm clmmpi_comm_new(void);

extern void clmmpi_destroy_handles(void);
extern void clmmpi_alloc_handle(int *handle, void **data);
extern void *clmmpi_handle_to_ptr(int handle);
extern void clmmpi_free_handle(int handle);

