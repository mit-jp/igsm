/* *************************************************************
****************************************************************
MITSOIL437.CPP - object describing general characteristics of 
                 soil

Modifications:

20060128 - Created by DWK 
                       
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


#include "mitsoil437.h"

/* *************************************************************
************************************************************* */

MITsoil43::MITsoil43( void ) : Tsoil43()
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

void MITsoil43::resetMonthlyTraceGasFluxes( void )
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

void MITsoil43::resetYrTraceGasFluxes( void )
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

