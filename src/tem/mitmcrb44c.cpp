/* **************************************************************
*****************************************************************
MITMCRB44C.CPP - object describing characteristics of soil 
                 microbes used in the Terrestrial Ecosystem Model 
                 (TEM)

Modifications:

20060128 - DWK created by modifying tmcrb437.cpp
20060128 - DWK changed include from tmcrb437.h to mitmcrb437.h
20080130 - DWK changed include from mitmcrb437.h to mitmcrb44a.h
20080130 - DWK changed Tmicrobe43:: to Tmicrobe44::
20081030 - DWK changed ProcessXML43 to ProcessXML44
20110707 - DWK changed include from mitmcrb44a.h to mitmcrb44c.h
                                  
*****************************************************************
************************************************************** */

#include<cstdio>

  using std::printf;

#include<iostream>

  using std::cin;
  using std::cout;
  using std::ios;
  using std::cerr;
  using std::endl;

#include<fstream>

  using std::ifstream;
  using std::ofstream;

#include<cstdlib>

  using std::exit;
  using std::atof;
  using std:: atoi;

#include<cmath>

  using std::exp;
  using std::log;
  using std::pow;

#include<vector>

  using std::vector;
    
#include<string>
  
  using std::string;


#include "mitmcrb44c.h"

/* **************************************************************
************************************************************** */

Tmicrobe44::Tmicrobe44() : ProcessXML44()
{

  // Initialize variables to MISSING values

  dq10 = MISSING;
  
  nuptake = MISSING;
  netnmin = MISSING;

  rh = MISSING;

  yrnmin = MISSING;
  yrnuptake = MISSING;

  yrrh = MISSING;

};

/* **************************************************************
************************************************************** */

void Tmicrobe44::getvegecd( ofstream& rflog1 )
{

  string ecd;

  cout << "Enter name of the file with microbe parameter ";
  cout << "values (.ECD)" << endl;
  cout << "dependent upon vegetation : " << endl;
  
  cin >> ecd;

  rflog1 << "Enter name of the file with microbe parameter ";
  rflog1 << "values (.ECD)" << endl;
  rflog1 << "dependent upon vegetation: " << ecd << endl << endl;

  getvegecd( ecd );

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void Tmicrobe44::getvegecd( const string& ecd )
{
  ifstream infile;
  int dcmnt;
  int comtype;


  infile.open( ecd.c_str(), ios::in );

  if( !infile )
  {
    cerr << endl << "Cannot open " << ecd;
    cerr << " for microbe ECD input" << endl;
    
    exit( -1 );
  }

  getXMLrootNode( infile, "microbeECD" );

  for( dcmnt = 1; dcmnt < MAXCMNT; ++dcmnt )
  {
    comtype = getXMLcommunityNode( infile, "microbeECD" );

    if( comtype >= MAXCMNT )
    {
      cerr << endl << "comtype is >= MAXCMNT" << endl;
      cerr << "comtype cannot be greater than " << (MAXCMNT-1);
      cerr << " in microbeECD" << endl;
      exit( -1 );
    }

    rhq10[comtype] = getXMLcmntArrayDouble( infile, 
                                            "microbeECD",
                                            "rhq10", 
                                            comtype );

    kn2[comtype] = getXMLcmntArrayDouble( infile, 
                                          "microbeECD",
                                          "kn2", 
                                          comtype );

    moistmin[comtype] = getXMLcmntArrayDouble( infile,
                                               "microbeECD",
                                               "moistmin",
                                               comtype );

    moistopt[comtype] = getXMLcmntArrayDouble( infile,
                                               "microbeECD",
                                               "moistopt",
                                               comtype );

    moistmax[comtype] = getXMLcmntArrayDouble( infile,
                                               "microbeECD",
                                               "moistmax",
                                               comtype );

    endXMLcommunityNode( infile );
  }

  if( dcmnt < MAXCMNT )
  {
    cerr << endl << " Parameters found for only " << dcmnt;
    cerr << " community types out of a maximum of ";
    cerr << (MAXCMNT-1) << " types in microbeECD" << endl;
    
    exit( -1 );
  }

  infile.close();

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

double Tmicrobe44::nminxclm( const int& pdcmnt,
                             const double& soilh2o,
                             const double& soilorgc,
                             const double& soilorgn,
                             const double& availn,
                             const double& ksoil )
{
  double gmin;
  double immb;
  double nmin;
  double tcnsoil;

  tcnsoil = cnsoil[pdcmnt];

  nuptake = ZERO;
  nmin = ZERO;
  if ( soilorgc > ZERO && soilorgn > ZERO )
  {
    immb  = (availn * ksoil) / soilh2o;
    immb /= (kn2[pdcmnt] + immb);
    
    nuptake = nup * immb * decay * rh;
    
    gmin = (soilorgn / soilorgc) * rh;
    
    nmin   = gmin - nuptake;

    if ( nmin >= ZERO ) 
    { 
      nmin *= (soilorgn/soilorgc) * tcnsoil; 
    }
    else { nmin *= (soilorgc/soilorgn) / tcnsoil; }
  }

  return nmin;

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void Tmicrobe44::resetEcds( const int& pcmnt, 
                            const double& psiplusc )
{
  kdc = (kda[pcmnt] * psiplusc) + kdb[pcmnt];

  if( kdc < ZERO ) { kdc = ZERO; }


  nup = (nupa[pcmnt] * psiplusc) + nupb[pcmnt];

  if( nup < ZERO ) { nup = ZERO; }
  

  // Determine the "decay" parameter

  decay = 0.26299 
          + (1.14757 * propftos[pcmnt])
          - (0.42956 * pow( propftos[pcmnt], 2.0 ));

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void Tmicrobe44::resetMonthlyFluxes( void )
{
  // Reset monthly fluxes to zero
  
  // Carbon fluxes
  
  rh = ZERO;

  // Nitrogen fluxes
  
  nuptake = ZERO;
  netnmin = ZERO;
	
};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void Tmicrobe44::resetYrFluxes( void )
{
  // Reset annual fluxes and summary variables to zero

  // Carbon fluxes
  
  yrrh = ZERO;

  // Nitrogen fluxes
  
  yrnuptake = ZERO;
  yrnmin = ZERO;
  
};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

double Tmicrobe44::rhxclm( const double& soilorgc,
                           const double& dq10,
                           const double& moist )
{

  return kd * soilorgc * moist * dq10;

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void Tmicrobe44::setDQ10( const int& pdcmnt, 
                          const double& tair, 
                          const double& snowpack )
{
  // dq10: effect of temperature on decomposition 

  if( snowpack > ZERO ) { dq10 = 1.0; }
  else
  {
    // Use air temperature to calculate dq10
    	
    dq10 = pow( rhq10[pdcmnt], tair / 10.0 );
  }
	
};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

double Tmicrobe44::setRHMOIST( const int& pdcmnt,
                               const double& pcfldcap,
                               const double& vsm,
                               const int& moistlim )
{
  
/* rhxh2o: effect of moisture on decomposition */

  double rhxh2o;
  double vfc;
  
  if( 0 == moistlim )
  {
    vfc = pcfldcap * 0.01;

    rhxh2o = (vfc - moistmin[pdcmnt]) 
             * (vfc - moistmax[pdcmnt]);

    rhxh2o /= rhxh2o - pow( (vfc - moistopt[pdcmnt]), 2.0 );
  }
  else
  {
    rhxh2o = (vsm - moistmin[pdcmnt]) 
             * (vsm - moistmax[pdcmnt]);

    rhxh2o /= rhxh2o - pow( (vsm - moistopt[pdcmnt]), 2.0 );
  }

  if( rhxh2o < ZERO ) { rhxh2o = ZERO; }

  return rhxh2o;
  
};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void Tmicrobe44::showecd( const int& pdcmnt )
{

  cout << endl << "   MICROBIAL PARAMETERS INFLUENCED BY CLIMATE";
  cout << endl << endl;
  printf( "          KN2 = %6.4lf\n", kn2[pdcmnt] );
  printf( "        RHQ10 = %6.2lf\n", rhq10[pdcmnt] );
  printf( "     MOISTMIN = %8.6lf\n", moistmin[pdcmnt] );
  printf( "     MOISTOPT = %8.6lf\n", moistopt[pdcmnt] );
  printf( "     MOISTMAX = %8.6lf\n", moistmax[pdcmnt] );
  printf( "       CNSOIL = %5.2lf\n", cnsoil[pdcmnt] );

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void Tmicrobe44::updateDynamics( const int& pdcmnt,
                                 const double& pcfldcap,
                                 const double& soilorgc,
                                 const double& soilorgn,
                                 const double& soilh2o,
                                 const double& vsm,
                                 const double& availn,
                                 const int& moistlim,
                                 const int& tillflag,
                                 const double& tillfactor,
                                 const double& ksoil )
{
  double rhmoist;

  // rhmoist: effect of moisture on decomposition

  rhmoist = setRHMOIST( pdcmnt, pcfldcap, vsm, moistlim );

  rh = rhxclm( soilorgc, dq10, rhmoist );

  if( rh < ZERO ) { rh = ZERO; }


  // Adjust decomposition rate due to tillage effects
  
  if( 1 == tillflag )
  {
    rh *= tillfactor;
  }


  // Determine Net N Mineralization (microbe.netnmin)

  netnmin = nminxclm( pdcmnt,
                      soilh2o,
                      soilorgc,
                      soilorgn,
                      availn,
                      ksoil );
	
};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

double Tmicrobe44::yrkd( const int& nfeed,
                         const double& yrltrc,
                         const double& yrltrn,
                         const int& pdcmnt )
{

  double yrkd;

  if( yrltrn <= 0.000000000000001 ) { return  yrkd = ZERO; }
  if( yrltrc < ZERO )
  {
    cout << "YRLTRC is < 0.0 in microbe.yrkd()" << endl;

    exit( -1 );
  }
  if( 0 == nfeed ) { yrkd = kdc; }
  else
  {
    yrkd = kdc * pow( (yrltrc/yrltrn), -0.784 )
           / pow( lcclnc[pdcmnt], -0.784 );
  }

  return yrkd;
};

