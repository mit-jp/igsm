/* *************************************************************
****************************************************************
MITSOIL44A.CPP - object describing general characteristics of 
                 soil

Modifications:

20070128 - Created by DWK by modifying mitsoil437.cpp
20070128 - DWK changed include from mitsoil437.h to mitsoil44.h
20070128 - DWK changed MITsoil43:: to MITsoil44::
20070128 - DWK changed inheritance from Tsoil43 to Tsoil44 
20080130 - DWK changed include from mitsoil44.h to mitsoil44a.h
                       
****************************************************************
************************************************************* */

#include<fstream>

  using std::ifstream;
  using std::ofstream;

#include<cmath>

  using std::exp;
  using std::fabs;
//  using std::modf;
//  using std::pow;

#include<string>
  
  using std::string;


#include "mitsoil44a.h"

/* *************************************************************
************************************************************* */

MITsoil44::MITsoil44( void ) : Tsoil44()
{

//  text  = -99;
//  wsoil = -99;

  // Initialize array with number of days per month
  
  mdays[0] = 31;
  mdays[1] = 28;
  mdays[2] = 31;
  mdays[3] = 30;
  mdays[4] = 31;
  mdays[5] = 30;
  mdays[6] = 31;
  mdays[7] = 31;
  mdays[8] = 30;
  mdays[9] = 31;
  mdays[10] = 30;
  mdays[11] = 31;
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITsoil44::resetMonthlyTraceGasFluxes( void )
{
  ch4consump = ZERO;
  ch4emissions = ZERO;
  ch4flux = ZERO;

  co2dnflux = ZERO;
  co2nflux = ZERO;
  
  n2odnflux = ZERO;
  n2onflux = ZERO;
  n2oflux = ZERO;

  n2flux = ZERO;
	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITsoil44::resetYrTraceGasFluxes( void )
{
  yrch4csmp = ZERO;
  yrch4ems = ZERO;
  yrch4flx = ZERO;

  yrco2dnflx = ZERO;
  yrco2nflx = ZERO;
  
  yrn2odnflx = ZERO;
  yrn2onflx = ZERO;
  yrn2oflx = ZERO;

  yrn2flx = ZERO;
	
};

