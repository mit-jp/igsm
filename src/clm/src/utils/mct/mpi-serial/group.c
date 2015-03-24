
#include "mpiP.h"


/*********/


FORT_NAME( clmmpi_group_incl, CLMMPI_GROUP_INCL )
     (int *group, int *n, int *ranks, int *newgroup, int *ierror)
{
  *ierror= CLMMPI_Group_incl(*group, *n, ranks, newgroup);
}


int CLMMPI_Group_incl(CLMMPI_Group group, int n, int *ranks, CLMMPI_Group *newgroup)
{

  if (group==CLMMPI_GROUP_NULL || group==CLMMPI_GROUP_EMPTY)
    *newgroup=group;
  else
    if (n==1 && ranks[0]==0)
      *newgroup=CLMMPI_GROUP_ONE;
    else
      *newgroup=CLMMPI_GROUP_NULL;


  return(CLMMPI_SUCCESS);
}


/*********/


FORT_NAME( clmmpi_group_free, CLMMPI_GROUP_FREE )(int *group, int *ierror)
{
  *ierror= CLMMPI_Group_free(group);
}


int CLMMPI_Group_free(CLMMPI_Group *group)
{
  *group= CLMMPI_GROUP_NULL;

  return(CLMMPI_SUCCESS);
}
