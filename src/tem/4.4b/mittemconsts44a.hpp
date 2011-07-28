// Global constants

#include "temigsmconsts44a.hpp"


// Maximum number of vegetation types in CLM
#ifndef CLMMXNVEG_CONST
#define CLMMXNVEG_CONST
// const int CLMMXNVEG = 19;
   const int CLMMXNVEG = 35;
#endif


// Maximum number of soil layers used in CLM
#ifndef CLMNLAYERS_CONST
#define CLMNLAYERS_CONST
  const int CLMNLAYERS = 10;
#endif


// Maximum number of hours per day
#ifndef MAXDAYHRS_CONST
#define MAXDAYHRS_CONST
  const int MAXDAYHRS = 24;    
#endif 

#ifndef MAX3HRS_CONST
#define MAX3HRS_CONST
  const int MAX3HRS = 8;
#endif 

#ifndef MAXHRS_CONST
#define MAXHRS_CONST
  const int MAXHRS = 24;
#endif

// Maximum number of days per month
#ifndef MAXMDAYS_CONST
#define MAXMDAYS_CONST
  const int MAXMDAYS = 31;
#endif

#ifndef MAXNGRD_CONST
#define MAXNGRD_CONST
  const int MAXNGRD = 46; 
#endif


// maximum number of 0.5 degree latitudinal bands

#ifndef MAXNLAT05_CONST
#define MAXNLAT05_CONST
  const int MAXNLAT05 = 360; 
#endif

#ifndef MXDAYHRS_CONST
#define MXDAYHRS_CONST
  const int MXDAYHRS = 24;
#endif

// maximum number of MIT-IGSM latitudinal bands

#ifndef MXMITNLAT_CONST
#define MXMITNLAT_CONST
  const int MXMITNLAT = 46;   
#endif

#ifndef MXMTHDAYS_CONST
#define MXMTHDAYS_CONST
  const int MXMTHDAYS = 31;   
#endif

//#ifndef MXNCHRTS_CONST
//#define MXNCHRTS_CONST
//  const int MXNCHRTS = 35;   
//#endif

// Maximum number of soil layers
#ifndef MXNLAYERS_CONST
#define MXNLAYERS_CONST
  const int MXNLAYERS = 6;
#endif

// Maximum number of soil layers
#ifndef NLAYERS_CONST
#define NLAYERS_CONST
  const int NLAYERS = 6;
#endif

// Maximum number of soil layers (i.e. levels) used by NEM 
//   for simulating N2O fluxes
#ifndef NLVL_CONST
#define NLVL_CONST
  const int NLVL = 10;
#endif
