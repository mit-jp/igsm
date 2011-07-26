/* **************************************************************
*****************************************************************
MITATMS437.CPP - object describes physical characteristics of the
	         atmosphere

Modifications:

20060127 - Created by DWK 

****************************************************************
************************************************************* */

 
#include "mitatms437.h"

MITatms43::MITatms43() : Atmosphere43()
{
  int dday;
  
  for( dday = 0; dday < MAXMDAYS; ++dday )
  {
    dayTair[dday] = MISSING;
    rainDuration[dday] = MISSING;
    rainIntensity[dday] = MISSING;
  }
  
};

