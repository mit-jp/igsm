/* **************************************************************
*****************************************************************
TIGSMATMS44C.CPP - object describes physical characteristics of 
                   the atmosphere

Modifications:

20060126 - DWK created by modifying atms50b5.cpp
20060126 - DWK changed include from atms50b5.h to atms437.h
20060126 - DWK changed Atmosphere50:: to class Atmosphere43::
20080130 - DWK changed include from atms437.h to tigsmatms44a.h
20080130 - DWK changed Atmosphere43:: to Atmosphere44::
20110707 - DWK changed include from tigsmatms44a.h to
           tigsmatms44c.h

****************************************************************
************************************************************* */

#include<cmath>
  using std::exp;

#include "tigsmatms44c.h"


/* *************************************************************
************************************************************* */
  

Atmosphere44::Atmosphere44()
{
  // Number of days per month
  
  ndays[0] = 31;
  ndays[1] = 28;
  ndays[2] = 31;
  ndays[3] = 30;
  ndays[4] = 31;
  ndays[5] = 30;
  ndays[6] = 31;
  ndays[7] = 31;
  ndays[8] = 30;
  ndays[9] = 31;
  ndays[10] = 30;
  ndays[11] = 31;
	
};

/* *************************************************************
************************************************************* */

/* *************************************************************
************************************************************* */

void Atmosphere44::petjh( const double& nirr, 
                          const double& tair,
                          const int& pdm )
{

  double f;
  double rt;

  f = ((9.0/5.0) * tair) + 32.0;

  rt = nirr * 0.016742;

  pet = ((0.014*f) - 0.37) * rt * ndays[pdm];

  if( pet <= ZERO ) { pet = 0.001; }

};

/* *************************************************************
************************************************************* */

/* *************************************************************
************************************************************* */

void Atmosphere44::precsplt( const double& prec, 
                             const double& tair,
                             double& rain, 
                             double& snowfall )
{


/* *************************************************************
	Willmott's assumptions on snow/rain split:
************************************************************** */

  if( tair >= -1.0 )
  {
    rain = prec;
    
    snowfall = ZERO;
  }
  else
  {
    rain = ZERO;
    
    snowfall = prec;
  }

};

/* *************************************************************
************************************************************* */

/* *************************************************************
************************************************************* */

void Atmosphere44::resetMonthlyFluxes( void )
{
  // Reset monthly fluxes to zero

//  pet = ZERO;

};
/* *************************************************************
************************************************************* */

/* *************************************************************
************************************************************* */

void Atmosphere44::resetYrFluxes( void )
{
  // Reset annual fluxes and summary variables to zero	

  yrrain = ZERO;

  yrsnowfall = ZERO;

  yrpet = ZERO;

};
                   

