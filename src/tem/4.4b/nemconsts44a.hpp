// Global constants

#include "mittemconsts44a.hpp"

// Maximum number of days per year
#ifndef MAXYRDAYS_CONST
#define MAXYRDAYS_CONST
  const int MAXYRDAYS = 365;   
#endif 

// Maximum number of zonal latitudinal bands 
#ifndef MAXZNBND_CONST
#define MAXZNBND_CONST
const int MAXZNBND = 24;      
#endif

#ifndef NELM_CONST
#define NELM_CONST
const int NELM = 5;  // Number of elements for CH4 emissions
#endif

#ifndef NLAY_CONST
#define NLAY_CONST
const int NLAY = 5;  // Number of soil layers
#endif

#ifndef NST_CONST
#define NST_CONST
const int NST = 12;  // Number of soil textures
#endif

// Ratio of Microbial Biomass/Active organic Carbon. 
#ifndef RBO_CONST
#define RBO_CONST
const double RBO = 0.02;
#endif

// Ratio of C/N in microbial biomass.
#ifndef RCNB_CONST
#define RCNB_CONST
const double RCNB = 8;
#endif

// Ratio of C/N in humads.
#ifndef RCNH_CONST
#define RCNH_CONST
const double RCNH = 8;
#endif

// Percentage of labile microbial biomass.
#ifndef SRB_CONST
#define SRB_CONST
const double SRB = 0.9;
#endif
