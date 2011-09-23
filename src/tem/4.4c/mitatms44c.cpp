/* **************************************************************
*****************************************************************
MITATMS44C.CPP - object describes physical characteristics of the
	         atmosphere

Modifications:

20060127 - Created by DWK 
20080130 - DWK changed include from mitatms437.h to mitatms44a.h
20080130 - DWK changed MITatms43:: to MITatms44::
20080130 - DWK changed Atmosphere43() to Atmosphere44()
20110707 - DWK changed include from mitatms44a.h to mitatms44c.h

****************************************************************
************************************************************* */

 
#include "mitatms44c.h"

MITatms44::MITatms44() : Atmosphere44()
{
  int dday;
  
  for( dday = 0; dday < MAXMDAYS; ++dday )
  {
    dayTair[dday] = MISSING;
    rainDuration[dday] = MISSING;
    rainIntensity[dday] = MISSING;
  }
  
};

