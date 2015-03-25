/* *************************************************************
****************************************************************
TIGSMSOIL44A.CPP - object describing general characteristics of 
                   soil

Modifications:

20060106 - DWK created by modifying tsoil432b.cpp 
20080130 - DWK changed include from tsoil44.h to tigsmsoil44a.h
20080130 - DWK changed ProcessXML43 to ProcessXML44
20110707 - DWK changed include from tsoil44a.h to tigsmsoil44c.h
                       
****************************************************************
************************************************************* */

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
  using std::pow;
  
#include<string>
  
  using std::string;


#include "tigsmsoil44c.h"

/* *************************************************************
************************************************************* */

Tsoil44::Tsoil44( void ) : ProcessXML44() 
{

  text  = -99;
  wsoil = -99;

  pctsand = MISSING;
  pctsilt = MISSING;
  pctclay = MISSING;
  psiplusc = MISSING;


  awcapmm = MISSING;
  fldcap = MISSING;
  wiltpt = MISSING;
  totpor = MISSING;

  snowpack = MISSING;

  avlh2o = MISSING;
  moist = MISSING;
  pcfc = MISSING;
  pctp = MISSING;
  vsm = MISSING;

  rgrndh2o = MISSING;
  sgrndh2o = MISSING;

  snowinf = MISSING;
  rperc = MISSING;
  sperc = MISSING;
  rrun = MISSING;
  srun = MISSING;
  h2oyld = MISSING;
    
  org.carbon = MISSING;
  org.nitrogen = MISSING;

  availn = MISSING;

  yrsnowpack = MISSING;

  yravlh2o = MISSING;
  yrsmoist = MISSING;
  yrpctp = MISSING;
  yrvsm = MISSING;
  meanvsm = MISSING;

  yrrgrndh2o = MISSING;
  yrsgrndh2o = MISSING;

  yrsnowinf = MISSING;
  yrrperc = MISSING;
  yrsperc = MISSING;
  yrrrun = MISSING;
  yrsrun = MISSING;
  yrh2oyld = MISSING;

  yrorgc = MISSING;
  yrorgn = MISSING;
  yrc2n = MISSING;

  yravln = MISSING;

  ninput = MISSING;
  yrnin = MISSING;

  nlost = MISSING;
  yrnlost = MISSING;

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

void Tsoil44::getecd( ofstream& rflog1 )
{

  string ecd;

  cout << "Enter name of the soil (.ECD) data file with parameter values: ";
  cout << endl;
  
  cin >> ecd;

  rflog1 << "Enter name of the soil (.ECD) data file with parameter values: ";
  rflog1 << ecd << endl;

  getecd( ecd );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil44::getecd( const string& ecd )
{
  ifstream infile;


  infile.open( ecd.c_str(), ios::in );

  if( !infile )
  {
    cerr << endl << "Cannot open " << ecd;
    cerr << " for soil ECD input" << endl;

    exit( -1 );
  }

  getXMLrootNode( infile, "soilECD" );

  pctpora = getXMLdouble( infile, "soilECD", "pctpora" );
  pctporb = getXMLdouble( infile, "soilECD", "pctporb" );

  fldcapa = getXMLdouble( infile, "soilECD", "fldcapa" );
  fldcapb = getXMLdouble( infile, "soilECD", "fldcapb" );

  wiltpta = getXMLdouble( infile, "soilECD", "wiltpta" );
  wiltptb = getXMLdouble( infile, "soilECD", "wiltptb" );

  infile.close();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil44::getrootz( ofstream& rflog1 )
{

  string ecd;

  cout << "Enter name of the data file containing the rooting depths:";
  cout << endl;
  cout << "               (e.g., ROOTZVEG.ECD)" << endl;
  
  cin >> ecd;

  rflog1 << "Enter name of the data file containing the rooting depths:";
  rflog1 << endl;
  rflog1 << "               (e.g., ROOTZVEG.ECD)" << endl;
  rflog1 << ecd << endl;

  getrootz( ecd );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil44::getrootz( const string& ecd )
{
  ifstream infile;
  int dcmnt;
  int comtype;


  infile.open( ecd.c_str(), ios::in );

  if( !infile )
  {
    cerr << endl;
    cerr << "Cannot open " << ecd << " for root ECD input";
    cerr << endl;

    exit( -1 );
  }

  getXMLrootNode( infile, "rootzECD" );

  for( dcmnt = 1; dcmnt < MAXCMNT; ++dcmnt )
  {

    comtype = getXMLcommunityNode( infile, "rootzECD" );

    if( comtype >= MAXCMNT )
    {
      cerr << endl << "comtype is >= MAXCMNT" << endl;
      cerr << "comtype cannot be greater than " << (MAXCMNT-1);
      cerr << " in leafECD" << endl;
      exit( -1 );
    }

    rootza[comtype] = getXMLcmntArrayDouble( infile,
                                             "rootzECD",
                                             "rootza",
                                             comtype );

    rootzb[comtype] = getXMLcmntArrayDouble( infile,
                                             "rootzECD",
                                             "rootzb",
                                             comtype );

    rootzc[comtype] = getXMLcmntArrayDouble( infile,
                                             "rootzECD",
                                             "rootzc",
                                             comtype );

    minrootz[comtype] = getXMLcmntArrayDouble( infile,
                                               "rootzECD",
                                               "minrootz",
                                               comtype );

    endXMLcommunityNode( infile );
  }

  if ( dcmnt < MAXCMNT )
  {
    cerr << endl << " Parameters found for only " << dcmnt;
    cerr << " community types out of a maximum of ";
    cerr << (MAXCMNT-1) << " types in rootzECD" << endl;

    exit( -1 );
  }

  infile.close();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil44::lake( const double& tair,
                    const double& prec,
                    double& rain,
                    double& snowfall,
       		    const double& pet,
                    double& eet )
{

  rgrndh2o = ZERO;
 
  sperc = ZERO;
 
  snowpack = ZERO;
 
  sgrndh2o = ZERO;
 
  moist = ZERO;

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

  eet = pet;

  h2oyld = prec - pet;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil44::percol( const double& rain, const double& avlh2o )
{

  double extra;

  double recharge;

  sperc = ZERO;

  rperc = ZERO;

  recharge = rain + snowinf;

  if( recharge <= ZERO ) { recharge = 0.001; }

  if( (avlh2o + rain + snowinf - eet) > awcapmm )
  {
    extra = rain + snowinf + avlh2o - awcapmm - eet;

    sperc = snowinf * extra / recharge;

    rperc = rain * extra / recharge;
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil44::resetMonthlyFluxes( void )
{
  // Reset monthly fluxes to zero
 
  // Nitrogen fluxes
  
  ninput = ZERO;

  nlost = ZERO;

  // Water fluxes

// Comment out next two lines in MITTEM as these values come
//   from CLM rather than TEM
  
//  ineet = ZERO;

//  eet = ZERO;

  rperc = ZERO;

  sperc = ZERO;

  rrun = ZERO;

  srun = ZERO;

  snowinf = ZERO;

  h2oyld = ZERO;
    	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil44::resetYrFluxes( void )
{
  // Reset annual fluxes and summary variables to zero
  
  // Annual carbon storage

  yrorgc = ZERO;

  // Annual nitrogen storage

  yrorgn = ZERO;

  yravln = ZERO;

  // Annual water storage

  yravlh2o = ZERO;

  yrsmoist = ZERO;

  yrvsm = ZERO;

  yrpctp = ZERO;

  yrsnowpack = ZERO;

  yrrgrndh2o = ZERO;

  yrsgrndh2o = ZERO;


  // Annual nitrogen fluxes

  yrnin = ZERO;

  yrnlost = ZERO;

  // Annual water fluxes

//  yrineet = ZERO;

  yreet = ZERO;

  yrrperc = ZERO;

  yrsperc = ZERO;

  yrrrun = ZERO;

  yrsrun = ZERO;

  yrsnowinf = ZERO;

  yrh2oyld = ZERO;
  	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Tsoil44::rrunoff( const double& rgrndh2o )
{
  double rrunof;

  rrunof = 0.5 * (rgrndh2o + rperc);

  return rrunof;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil44::setKH2O( const double& vsm, 
                       const int& moistlim )
{
  double vfc;
  
  if( 0 == moistlim ) 
  { 
    vfc = pcfldcap * 0.01;

    kh2o = pow( vfc, 3.0 ); 
  }
  else 
  { 
    if( vsm > 1.0 ) { kh2o = 1.0; }
    else { kh2o = pow( vsm, 3.0 ); }
  }

	
};

/* *************************************************************
************************************************************* */

/* *************************************************************
************************************************************* */

void Tsoil44::showecd( void )
{

  cout << endl << "                   SOIL CHARACTERISTICS OF SITE";
  cout << endl << endl;

  printf( "PSAND    = %5.2lf      PSILT = %5.2lf      PCLAY = %5.2lf\n",
          pctsand,
          pctsilt,
          pctclay );

  printf( "POROSITY = %5.2lf   PCFLDCAP = %5.2lf   PCWILTPT = %5.2lf\n",
          pctpor,
          pcfldcap,
          pcwiltpt );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Tsoil44::snowmelt( const double& elev,
                          const double& tair,
                          const double& prevtair,
                          const double& psnowpack )
{

  double snowflux = ZERO;

  if( tair >= -1.0 )
  {
    if( elev <= 500.0 ) { snowflux = psnowpack;}
    else
    {
      if( prevtair < -1.0 ) { snowflux = 0.5 * psnowpack; }
      else { snowflux = psnowpack; }
    }
  }

  if( snowflux < ZERO ) { snowflux = ZERO; }

  return snowflux;
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Tsoil44::srunoff( const double& elev,
                         const double& tair,
                         const double& prevtair,
                         const double& prev2tair,
                         const double& sgrndh2o )
{

  double srunof = ZERO;

  if( tair >= -1.0 )
  {
    if( prevtair < -1.0 )
    {
      srunof = 0.1 * (sgrndh2o + sperc);
    }
    else
    {
      if( prev2tair < -1.0 )
      {
	if( elev <= 500.0 )
        {
          srunof = 0.5 * (sgrndh2o + sperc);
        }
	else { srunof = 0.25 * (sgrndh2o + sperc); }
      }
      else { srunof = 0.5 * (sgrndh2o + sperc); }
    }
  }

  return srunof;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil44::updateHydrology( const double& elev,
                               const double& tair,
                               const double& prevtair,
                               const double& prev2tair,
                               const double& rain,
                               const double& pet,
                               const double& sh2o,
                               const double& rgrndh2o,
                               const double& sgrndh2o,
                               const int& irrgflag,
                               double& irrigate,
                               const int& pdm )
{
  double xrain;

  // Determine available water (avlh2o)
    
  avlh2o = sh2o - wiltpt; 

  if( avlh2o < ZERO )
  {
    avlh2o = ZERO;
  }


  // Determine initial eet (ineet)
  
//  ineet = xeet( rain, pet, avlh2o, pdm );
  eet = xeet( rain, pet, avlh2o, pdm );


//  if( ineet > (rain + snowinf + avlh2o) )
//  {
//    ineet = rain + snowinf + avlh2o;
//  }

  if( eet < ZERO ) { eet = ZERO; }

  if( eet > (rain + snowinf + avlh2o) )
  {
    eet = rain + snowinf + avlh2o;
  }

//  if( ineet < ZERO ) { ineet = ZERO; }
  if( eet < ZERO ) { eet = ZERO; }

//  if( 1 == irrgflag && ineet < pet )
  if( 1 == irrgflag && eet < pet )
   {
    // If irrigated, add just enough water to overcome 
    //   moisture limitations

//    irrigate = pet - ineet;
    irrigate = pet - eet;
    
    xrain = rain + irrigate;


    // Determine adjusted eet when irrigated

    eet = xeet( xrain, pet, avlh2o, pdm );
  }
  else
  {
    irrigate = ZERO;

    xrain = rain;

//    eet = ineet;
  }


  // Determine monthly percent total porosity (pctp)

  pctp = (100.0 * sh2o)  / totpor;


  // Determine volumetric soil moisture (vsm)
  
  vsm = sh2o / (rootz * 1000.0);
 
  if( vsm <= ZERO ) { vsm = 0.001; }


  // Determine percolation of rain water (rperc) and snow melt 
  //   water (sperc) through the soil profile
  
  percol( xrain, avlh2o );

  if( (avlh2o + snowinf + rain + irrigate - eet - rperc - sperc)
       < ZERO )
  {
    eet = avlh2o + snowinf + rain + irrigate - rperc - sperc;
  }


  // Determine runoff derived from rain (soil.rrun) and/or
  //   snow (soil.srun)

  rrun = rrunoff( rgrndh2o );

  srun = srunoff( elev,
                  tair,
                  prevtair,
                  prev2tair,
                  sgrndh2o );
	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil44::updateNLosses( const int& pdcmnt,
                             const double& h2oloss,
                             const double& availn, 
                             const double& soilh2o )
{
  if( soilh2o > ZERO )
  {
    nlost = availn / soilh2o;
  
    nlost *= (h2oloss + (rootz * 1000.0)) 
             / (rootz * 1000.0);
  
    nlost *= nloss[pdcmnt];
  }
  else { nlost = ZERO; }	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil44::updateRootZ( const int& pdcmnt )
{
  rootz = (rootza[pdcmnt] * pow( psiplusc, 2.0 ))
          + (rootzb[pdcmnt] * psiplusc)
          + rootzc[pdcmnt];

  if( rootz < minrootz[pdcmnt] ) { rootz = minrootz[pdcmnt]; }

  pctpor = (pctpora * psiplusc) + pctporb;

  pcfldcap = (fldcapa * psiplusc) + fldcapb;

  pcwiltpt = (wiltpta * psiplusc) + wiltptb;

  totpor  = rootz * pctpor * 10.0;

  fldcap  = rootz * pcfldcap * 10.0;

  wiltpt  = rootz * pcwiltpt * 10.0;

  awcapmm = fldcap - wiltpt;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

double Tsoil44::xeet( const double& rain, 
                      const double& pet,
                      const double& avlh2o,  
                      const int& pdm )
{

  const double edpar = 5.0;

  double aet;
  double def;
  double dsm;
  double ep;
  double gm;
  double prob;
  double rbar;

  if( (rain+snowinf) >= pet )
  {
    aet = pet;
  }
  else
  {
    gm = (1.0 - exp( -edpar * avlh2o / awcapmm) ) 
         / (1.0 - exp( -edpar ));

    ep = pet / ndays[pdm];
    def = ep + awcapmm - avlh2o;
    prob = 1.0 - exp( -0.005*(rain + snowinf) );

    if( prob != ZERO ) 
    { 
      rbar = (rain + snowinf) / (ndays[pdm] * prob); 
    }
    else { rbar = ZERO; }

    if ( rbar != ZERO )
    {
      dsm = rbar*prob*(gm + ((1.0-gm) * exp( -ep/rbar )) 
            - exp( -def/rbar )) 
            - (ep*gm);
    }
    else { dsm = -ep*gm; }

    dsm *= ndays[pdm];

    aet = rain + snowinf - dsm;
    
    if( aet > pet ) { aet = pet; }
  }

  return aet;
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil44::xtext( const int& pdcmnt,
                     const double& pctsilt,
                     const double& pctclay )
{

  totpor = fldcap = wiltpt = MISSING;

  awcapmm =  MISSING;

  psiplusc = (pctsilt + pctclay) * 0.01;

  if( psiplusc < 0.01 ) { psiplusc = 0.01; }

  rootz = (rootza[pdcmnt] * pow( psiplusc, 2.0 ))
          + (rootzb[pdcmnt] * psiplusc)
          + rootzc[pdcmnt];

  if( rootz < minrootz[pdcmnt] ) { rootz = minrootz[pdcmnt]; }

  pctpor = (pctpora * psiplusc) + pctporb;

  pcfldcap = (fldcapa * psiplusc) + fldcapb;

  pcwiltpt = (wiltpta * psiplusc) + wiltptb;

  totpor  = rootz * pctpor * 10.0;

  fldcap  = rootz * pcfldcap * 10.0;

  wiltpt  = rootz * pcwiltpt * 10.0;

  awcapmm = fldcap - wiltpt;

};

