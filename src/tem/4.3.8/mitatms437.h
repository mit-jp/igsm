/* **************************************************************
*****************************************************************
MITATMS437.H - object describes physical characteristics of the
	       atmosphere

Modifications:

20060127 - Created by DWK
                                                                                                                                           
*****************************************************************
************************************************************** */

#ifndef MITATMS437_H
#define MITATMS437_H

// Global constants

#include "mittemconsts43.hpp"

// MITatms43 inherits Atmosphere43

#include "atms437.h"

class MITatms43 : public Atmosphere43
{

  public:

     MITatms43();
  
/* *************************************************************
		 Public Variables
************************************************************* */

     // Daily mean air temperature (degrees C)
//     double dayTair[CYCLE][MAXMDAYS];
     double dayTair[MAXMDAYS];

     double nox;

     // Duration of current rain event (hrs)
//     double rainDuration[CYCLE][MAXMDAYS];
     double rainDuration[MAXMDAYS];

     // Intensity of rain during current event (mm/hr)
//     double rainIntensity[CYCLE][MAXMDAYS];
     double rainIntensity[MAXMDAYS];

};

#endif
