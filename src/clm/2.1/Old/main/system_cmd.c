#include <misc.h>
#include <cfort.h>
#include <stdlib.h>

#if ( defined FORTRANCAPS )
#define system_cmd SYSTEM_CMD
#elif ( defined FORTRANUNDERSCORE )
#define system_cmd system_cmd_
#elif ( defined FORTRANDOUBLEUNDERSCORE )
#define system_cmd system_cmd__
#endif

int system_cmd ( const char *text)
/*

	Wrapper to "C" system command. Allows
	all platforms to have a consistent interface
	that includes a error return code.

*/
{
  int err;   /* Error return code */

  err = system(text);
  return( err );
}

                 
