/* *************************************************************
****************************************************************
TIGSMHUMANACT44C.CPP - describes human disturbances to natural 
                       ecosystems

Modifications:

20060423 - DWK created by modifying humnact437.cpp
20060423 - DWK changed include humnact437.h to humnact44.h
20060423 - DWK changed Humanact43:: to Humanact44::
20060423 - DWK added resetMonthlyDisturbFluxes()
20060423 - DWK deleted const int& pdm from function call to 
           resetMonthlyFluxes()
20060423 - DWK added setFireNDEP()
20080130 - DWK changed include from humnact44.h to 
           tigsmhumnact44a.h
20080130 - DWK changed ProcessXML43() to ProcessXML44()
20080826 - DWK changed include from tigsmhumnact44a.h to
           tigsmhumnact44b.h
20110706 - DWK changed include from tigsmhumnact44b.h to
           tigsmhumnact44c.h
20110819 - DWK added slashpar, vconvert, prod10par, prod100par,
                                             
****************************************************************
************************************************************* */

#include<cstdio>

  using std::printf;

#include<iostream>

  using std::cout;
  using std::cin;
  using std::ios;
  using std::cerr;
  using std::endl;

#include<fstream>

  using std::ifstream;
  using std::ofstream;

#include<cstdlib>

  using std::exit;
  using std::atof;
  using std::atoi;

#include<cmath>

  using std::exp;
  using std::pow;
    
#include<string>
  
  using std::string;

#include "tigsmhumnact44c.h"


/* *************************************************************
************************************************************* */

Humanact44::Humanact44() : ProcessXML44()
{

  int i;
  int dm;
  c2n = 54.29;
//  cfall = 0.20;
//  nfall = 0.20;

  state = 0;
  prvstate = 0;

  massbalflg = 1;
  fertflag = 0;
  irrgflag = 0;
  frostflag = 0;
  
  mez = CROPVEG - 1;

  productYear = 0;

  prevPROD1.carbon = ZERO;
  prevPROD1.nitrogen = ZERO;

  prevPROD10.carbon = ZERO;
  prevPROD10.nitrogen = ZERO;

  prevPROD100.carbon = ZERO;
  prevPROD100.nitrogen = ZERO;

  prevCropResidue.carbon = ZERO;
  prevCropResidue.nitrogen = ZERO;

  cropResidue.carbon = ZERO;
  cropResidue.nitrogen = ZERO;

  for( dm = 0; dm < CYCLE; ++dm )
  {  
    initPROD1[dm].carbon = ZERO;
    initPROD1[dm].nitrogen = ZERO;
  }

  for( i = 0; i < 10; ++i )
  {
    initPROD10[i].carbon = ZERO;
    initPROD10[i].nitrogen = ZERO;
  }
  
  for( i = 0; i < 100; ++i )
  {
    initPROD100[i].carbon = ZERO;
    initPROD100[i].nitrogen = ZERO;
  }

  forage.carbon = ZERO;
  forage.nitrogen = ZERO;

  manure.carbon = ZERO;
  manure.nitrogen = ZERO;

  animalresp = ZERO;

  urine = ZERO;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact44::conversion( const int& pdcmnt,
                             const double& vegc,
                             const double& vstrn,
                             const double& vston,
                             const double& soilc,
                             const double& soiln )
{
  // Assume annual conversion fluxes are partitioned equally 
  //   across each of 12 months - determine monthly fluxes
  
  slash.carbon = slashpar[pdcmnt] * vegc / (double) CYCLE;

  slash.nitrogen = slashpar[pdcmnt]  * (vstrn + vston) / (double) CYCLE;

  vconvrtflx.carbon = (vconvert[pdcmnt] * vegc) / (double) CYCLE;

  sconvrtflx.carbon =  (sconvert[pdcmnt] * soilc)/ (double) CYCLE;

  convrtflx.carbon = vconvrtflx.carbon + sconvrtflx.carbon;

  vconvrtflx.nitrogen = ((1.0 - nvretconv[pdcmnt]) 
                        * vconvert[pdcmnt] 
                        * (vstrn + vston))
                        / (double) CYCLE;
                                
  sconvrtflx.nitrogen = ((1.0 - nsretconv[pdcmnt]) 
                        * sconvert[pdcmnt] * soiln)
                        / (double) CYCLE;
                                
  convrtflx.nitrogen = vconvrtflx.nitrogen + sconvrtflx.nitrogen;
                               
  nvretent = (nvretconv[pdcmnt] * vconvert[pdcmnt] * (vstrn + vston)) 
              / (double) CYCLE;

  nsretent = (nsretconv[pdcmnt] * sconvert[pdcmnt] * soiln) 
              / (double) CYCLE;
  
  nretent = nvretent + nsretent;
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact44::createWoodProducts( const int& pdcmnt,
                                     const int& pdyr,
                                     const double& vegc,
                                     const double& vstrn,
                                     const double& vston )
{
  int i;

  formPROD10.carbon  = prod10par[pdcmnt] * vegc;

  formPROD10.nitrogen  = prod10par[pdcmnt] * (vstrn + vston);

  formPROD100.carbon = prod100par[pdcmnt] * vegc;

  formPROD100.nitrogen = prod100par[pdcmnt] * (vstrn + vston);

  if( pdyr < productYear )
  {
    productYear += pdyr;
  }
  else { productYear = pdyr; }

  i = productYear%10;

  initPROD10[i].carbon = formPROD10.carbon;

  initPROD10[i].nitrogen = formPROD10.nitrogen;

  i = productYear%100;

  initPROD100[i].carbon = formPROD100.carbon;

  initPROD100[i].nitrogen = formPROD100.nitrogen;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact44::decayProducts( void )
{

  int i;
  int dm;
  double yrtrash;

  
  PROD1decay.carbon  = ZERO;

  PROD1decay.nitrogen  = ZERO;

  for( dm = 0; dm < CYCLE; ++dm )
  {
    PROD1decay.carbon += initPROD1[dm].carbon / (double) CYCLE;
                        
    PROD1decay.nitrogen += initPROD1[dm].nitrogen 
                          / (double) CYCLE;
  }

  
  PROD10decay.carbon  = ZERO;

  PROD10decay.nitrogen  = ZERO;

  for( i = 0; i < 10; ++i )
  {
    yrtrash = initPROD10[i].carbon * 0.10 / (double) CYCLE;

    PROD10decay.carbon += yrtrash;

    yrtrash = initPROD10[i].nitrogen * 0.10 / (double) CYCLE;

    PROD10decay.nitrogen += yrtrash;
  }


  PROD100decay.carbon = ZERO;

  PROD100decay.nitrogen = ZERO;
  
  for( i = 0; i < 100; ++i )
  {
    yrtrash = initPROD100[i].carbon * 0.01 / (double) CYCLE;

    PROD100decay.carbon += yrtrash;

    yrtrash = initPROD100[i].nitrogen * 0.01 / (double) CYCLE;

    PROD100decay.nitrogen += yrtrash;
  }

  TOTPRODdecay.carbon = PROD1decay.carbon
                        + PROD10decay.carbon
                        + PROD100decay.carbon;

  TOTPRODdecay.nitrogen = PROD1decay.nitrogen
                          + PROD10decay.nitrogen
                          + PROD100decay.nitrogen;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact44::frostDamage( const double& vegc,
                              const double& vstrn, 
                              const double& vston )
{ 
  frostflag = 1;
  
  // Assume all crop biomass becomes stubble
  
  stubble.carbon = vegc;

  stubble.nitrogen = (vstrn+vston);

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact44::getecd( ofstream& rflog1 )
{
  string ecd;

  cout << "Enter name of the data file (.ECD) with agricultural ";
  cout << "parameter values: " << endl;

  cin >> ecd;

  rflog1 << "Enter name of the data file (.ECD) with agricultural ";
  rflog1 << "parameter values: " << ecd << endl;

  getecd( ecd );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact44::getecd( const string& ecd )
{
  int comtype;

  int dcmnt;

  ifstream infile;

  infile.open( ecd.c_str(), ios::in );

  if( !infile )
  {
    cerr << "\nCannot open " << ecd;
    cerr << " for agriculture ECD input" << endl;
    
    exit( -1 );
  }
  
  getXMLrootNode( infile, "agECD" );

  for( dcmnt = 1; dcmnt < MAXCMNT; ++dcmnt )
  {
    comtype = getXMLcommunityNode( infile, "agECD" );

    if( comtype >= MAXCMNT )
    {
      cerr << endl << "comtype is >= MAXCMNT" << endl;
      cerr << "comtype cannot be greater than " << (MAXCMNT-1) << endl;
      cerr << " in agECD" << endl;
      
      exit( -1 );
    }
    

    slashpar[comtype] = getXMLcmntArrayDouble( infile,
                                               "agECD",
                                               "slashpar",
                                               comtype );

    vconvert[comtype] = getXMLcmntArrayDouble( infile,
                                               "agECD",
                                               "vconvert",
                                               comtype );

    prod10par[comtype] = getXMLcmntArrayDouble( infile,
                                                "agECD",
                                                "prod10par",
                                                comtype );

    prod100par[comtype] = getXMLcmntArrayDouble( infile,
                                                 "agECD",
                                                 "prod100par",
                                                 comtype );
 
    sconvert[comtype] = getXMLcmntArrayDouble( infile,
                                               "agECD",
                                               "sconvert",
                                               comtype );

    nvretconv[comtype] = getXMLcmntArrayDouble( infile,
                                                "agECD",
                                                "nvretconv",
                                                comtype );

    nsretconv[comtype] = getXMLcmntArrayDouble( infile,
                                                "agECD",
                                                "nsretconv",
                                                comtype );

    vrespar[comtype] = getXMLcmntArrayDouble( infile,
                                              "agECD",
                                              "vrespar",
                                              comtype );

    tillfactor[comtype] = getXMLcmntArrayDouble( infile,
                                                 "agECD",
                                                 "tillfactor",
                                                 comtype );

    harvstC[comtype] = getXMLcmntArrayDouble( infile,
                                              "agECD",
                                              "harvstC",
                                              comtype );

    harvstN[comtype] = getXMLcmntArrayDouble( infile,
                                              "agECD",
                                              "harvstN",
                                              comtype );

    residueC[comtype] = getXMLcmntArrayDouble( infile,
                                               "agECD",
                                               "residueC",
                                               comtype );

    residueN[comtype] = getXMLcmntArrayDouble( infile,
                                               "agECD",
                                               "residueN",
                                               comtype );

    cropseedC[comtype] = getXMLcmntArrayDouble( infile,
                                                "agECD",
                                                "cropseedC",
                                                comtype );

    cropseedSTRN[comtype] = getXMLcmntArrayDouble( infile,
                                                   "agECD",
                                                   "cropseedSTRN",
                                                   comtype );

    cropseedSTON[comtype] = getXMLcmntArrayDouble( infile,
                                               "agECD",
                                               "cropseedSTON",
                                               comtype );
    endXMLcommunityNode( infile );
  }

  if( dcmnt < MAXCMNT )
  {
    cerr << endl << " Parameters found for only " << dcmnt;
    cerr << " community types out of a maximum of ";
    cerr << (MAXCMNT-1) << " types in agECD" << endl;
    
    exit( -1 );
  }

  infile.close();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact44::grazing( const double& vegc, 
                          const double& strn )
{
  // Assume 5 percent of vegetation biomass is consumed each 
  //   month

  forage.carbon = vegc * 0.05;

  forage.nitrogen = strn * 0.05;

  // Assume 83 percent of forage carbon is respired and 50
  //   percent of nitrogen is mineralized in livestock
  //   metabolism

  animalresp = forage.carbon * 0.83;

  urine = forage.nitrogen * 0.50;

  // Assume remaineder of forage is returned to soil as manure

  manure.carbon = forage.carbon * 0.17;

  manure.nitrogen = forage.nitrogen * 0.50;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact44::harvest( const int& pdm,
                          const double& vegc,
                          const double& vstrn, 
                          const double& vston )
{
//  int kdm;
//  int idm;

  double initresidueC;
  double initresidueN;

//  double residueCflux;
//  double residueNflux;
//  double residueNretent;

  // Determine crop production

  cropprod.carbon = vegc * harvstC[cmnt];

  if( cropprod.carbon < ZERO )
  {
    cropprod.carbon = ZERO;
  }

  cropprod.nitrogen = (vstrn+vston) * harvstN[cmnt];

  if( cropprod.nitrogen < ZERO )
  {
    cropprod.nitrogen = ZERO;
  }

  initPROD1[pdm].carbon  = cropprod.carbon;

  initPROD1[pdm].nitrogen  = cropprod.nitrogen;
  
  // Determine amount of carbon and nitrogen left
  // in crop residue

   initresidueC = vegc - cropprod.carbon;

  if( initresidueC < ZERO )
  {
    initresidueC = ZERO;
  }

  initresidueN = (vstrn+vston) - cropprod.nitrogen;

  if( initresidueN < ZERO )
  {
    initresidueN = ZERO;
  }

  // Determine amount of carbon and nitrogen in
  // crop residue that will be lost from ecosystem

  formCropResidue.carbon = initresidueC
                           * residueC[cmnt];

  formCropResidue.nitrogen = initresidueN
                             * residueN[cmnt];


  // Determine amount of stubble carbon and nitrogen
  // added to soils

  stubble.carbon = initresidueC
                   * ( 1.000000 - residueC[cmnt]);

  stubble.nitrogen = initresidueN
                     * ( 1.000000 - residueN[cmnt]);

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact44::resetMonthlyDisturbFluxes( void )
{
  // Initialize disturbance carbon fluxes to zero

  convrtflx.carbon = ZERO;

  vconvrtflx.carbon = ZERO;

  sconvrtflx.carbon = ZERO;

  slash.carbon = ZERO;

  // Initialize nitrogen fluxes during conversion to zero

  convrtflx.nitrogen = ZERO;

  vconvrtflx.nitrogen = ZERO;

  sconvrtflx.nitrogen = ZERO;

  slash.nitrogen = ZERO;

  nretent = ZERO;

  nvretent = ZERO;

  nsretent = ZERO;

  cropResidueFlux.carbon = ZERO;

  cropResidueFlux.nitrogen = ZERO;  
	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact44::resetMonthlyFluxes( void )
{
  // Reset monthly fluxes to zero
  
  fertn = ZERO;
  
  irrigate = ZERO;

  // Initialize stubble to zero
  
  stubble.carbon = ZERO;

  stubble.nitrogen = ZERO;

  cflux = ZERO;

  
  // Initialize carbon fluxes related to formation and 
  //   decomposition of products and crop residue to zero

  cropprod.carbon = ZERO;
  cropprod.nitrogen = ZERO;

  formCropResidue.carbon = ZERO;
  formCropResidue.nitrogen = ZERO;

  PROD1decay.carbon = ZERO;
  PROD1decay.nitrogen = ZERO;

  PROD10decay.carbon = ZERO;
  PROD10decay.nitrogen = ZERO;
 
  PROD100decay.carbon = ZERO;
  PROD100decay.nitrogen = ZERO;

  formTOTPROD.carbon = ZERO;
  formTOTPROD.nitrogen = ZERO;

  TOTPRODdecay.carbon = ZERO;
  TOTPRODdecay.nitrogen = ZERO;

	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact44::resetPROD( void )
{
  int dm;

  int i;

  prevPROD1.carbon = ZERO;
  prevPROD1.nitrogen = ZERO;

  prevPROD10.carbon = ZERO;
  prevPROD10.nitrogen = ZERO;

  prevPROD100.carbon = ZERO;
  prevPROD100.nitrogen = ZERO;

  prevCropResidue.carbon = ZERO;
  prevCropResidue.nitrogen = ZERO;

  cropResidue.carbon = ZERO;
  cropResidue.nitrogen = ZERO;


  for( dm = 0; dm < CYCLE; ++dm )
  {
    initPROD1[dm].carbon = ZERO;
    initPROD1[dm].nitrogen = ZERO;
  }

  for( i = 0; i < 10; ++i )
  {
    initPROD10[i].carbon = ZERO;
    initPROD10[i].nitrogen = ZERO;
  }

  for( i = 0; i < 100; ++i )
  {
    initPROD100[i].carbon = ZERO;
    initPROD100[i].nitrogen = ZERO;
  }

  PROD1.carbon = ZERO;
  PROD1.nitrogen = ZERO;

  PROD10.carbon = ZERO;
  PROD10.nitrogen = ZERO;

  PROD100.carbon = ZERO;
  PROD100.nitrogen = ZERO;

  TOTPROD.carbon = ZERO;
  TOTPROD.nitrogen = ZERO;

  PROD1decay.carbon = ZERO;
  PROD1decay.nitrogen = ZERO;

  formPROD10.carbon  = ZERO;
  formPROD10.nitrogen  = ZERO;

  PROD10decay.carbon  = ZERO;
  PROD10decay.nitrogen  = ZERO;

  formPROD100.carbon = ZERO;
  formPROD100.nitrogen = ZERO;

  PROD100decay.carbon = ZERO;
  PROD100decay.nitrogen = ZERO;

  formTOTPROD.carbon = ZERO;
  formTOTPROD.nitrogen = ZERO;

  TOTPRODdecay.carbon = ZERO;
  TOTPRODdecay.nitrogen = ZERO;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact44::resetYrFluxes( void )
{
  // Reset annual fluxes and summary variables to zero

  yrfertn = ZERO;

  yrirrig = ZERO;

  yrstubC = ZERO;

  yrstubN = ZERO;

  // Annual carbon fluxes from agricultural conversion

  yrconvrtC = ZERO;

  yrvconvrtC = ZERO;

  yrsconvrtC = ZERO;

  yrslashC = ZERO;

  yrcflux = ZERO;

 // Annual nitrogen fluxes from agricultural conversion

  yrconvrtN = ZERO;

  yrvconvrtN = ZERO;

  yrsconvrtN = ZERO;

  yrslashN = ZERO;

  yrnrent = ZERO;

  yrnvrent = ZERO;

  yrnsrent = ZERO;

  // Annual carbon and nitrogen fluxes in the formation of 
  //   agricultural products

  yrformPROD1C   = ZERO;
  yrformPROD1N   = ZERO;


 // Annual carbon and nitrogen flux from crop residue formation

  yrformResidueC = ZERO;
  yrformResidueN = ZERO;

  // Annual carbon and nitrogen fluxes in the decomposition of 
  //   agricultural products

  yrdecayPROD1C   = ZERO;
  yrdecayPROD1N   = ZERO;

 // Annual carbon and nitrogen fluxes from burning crop residue

  yrfluxResidueC = ZERO;
  yrfluxResidueN = ZERO;

  // Annual carbon and nitrogen fluxes resulting from 
  //   formation of 10-year wood products

  yrformPROD10C  = ZERO;
  yrformPROD10N  = ZERO;

  // Annual carbon and nitrogen fluxes resulting from 
  //   decomposition of 10-year wood products

  yrdecayPROD10C  = ZERO;
  yrdecayPROD10N  = ZERO;

  // Annual carbon and nitrogen fluxes resulting from 
  //   formation of 100-year wood products

  yrformPROD100C = ZERO;
  yrformPROD100N = ZERO;

  // Annual carbon and nitrogen fluxes resulting from 
  //   decomposition of 100-year wood products

  yrdecayPROD100C = ZERO;
  yrdecayPROD100N = ZERO;

  // Annual carbon and nitrogen fluxes resulting from 
  //   formation of all agricultural and wood products

  yrformTOTPRODC = ZERO;
  yrformTOTPRODN = ZERO;

  // Annual carbon and nitrogen fluxes resulting from 
  //   decomposition of all agricultural and wood products

  yrdecayTOTPRODC = ZERO;
  yrdecayTOTPRODN = ZERO;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact44::setAgricFlags( ofstream& rflog1 )
{

  cout << "Are agricultural soils tilled?" << endl;
  cout << "Enter 0 for no:" << endl;
  cout << "Enter 1 for yes:" << endl;
  
  cin >> tillflag;

  rflog1 << "Are agricultural soils tilled?" << endl;
  rflog1 << "Enter 0 for no:" << endl;
  rflog1 << "Enter 1 for yes:" << endl;
  rflog1 << tillflag << endl << endl;

};

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */

void Humanact44::setFireNDEP( void )
{
  double frftemp;

  double newndep;

  // Determine firendep from total nitrogen volatilized during 
  //   fire for reintroduction to the soil (Note: 
  //   vconvrtflx.nitrogen and sconvrtflx.nitrogen are monthly
  //   fluxes that occur over a time period of a year
  
  newndep = vconvrtflx.nitrogen + sconvrtflx.nitrogen;

  
  // Determine the potential number of years until next fire 
  //   event (i.e. Fire Return Interval or FRI)
  
  frftemp = FRI;

  if( frftemp > MAXFRI ) { frftemp = MAXFRI; }

  
  // Assume that (1/FRI) of newndep is deposited each month
  
  firendep = newndep / (frftemp-1);
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact44::setNoCrops( const int& pdm )
{
  cropprod.carbon = ZERO;
  cropprod.nitrogen = ZERO;

  formCropResidue.carbon = ZERO;
  formCropResidue.nitrogen = ZERO;

  if( frostflag != 1 )
  {
    stubble.carbon = ZERO;
    stubble.nitrogen = ZERO;
  }
  
  frostflag = 0;
  
  initPROD1[pdm].carbon = ZERO;
  initPROD1[pdm].nitrogen = ZERO;
  
  	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact44::setNoGrazing( void )
{
  forage.carbon = ZERO;
  forage.nitrogen = ZERO;

  manure.carbon = ZERO;
  manure.nitrogen = ZERO;

  animalresp = ZERO;
 
  urine = ZERO;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact44::setNoWoodProducts( const int& pdyr )
{
  int i;
  	
  formPROD10.carbon = ZERO;
  formPROD10.nitrogen = ZERO;

  formPROD100.carbon = ZERO;
  formPROD100.nitrogen = ZERO;
  
  if( pdyr < productYear )
  {
    productYear += pdyr;
  }
  else { productYear = pdyr; }

  i = productYear%10;

  initPROD10[i].carbon = ZERO;
  initPROD10[i].nitrogen = ZERO;

  i = productYear%100;

  initPROD100[i].carbon = ZERO;
  initPROD100[i].nitrogen = ZERO;
	
};


/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact44::updateCropResidue( void )
{

  cropResidue.carbon = prevCropResidue.carbon
                       + formCropResidue.carbon
                       - cropResidueFlux.carbon;

  cropResidue.nitrogen = prevCropResidue.nitrogen
                         + formCropResidue.nitrogen
                         - cropResidueFlux.nitrogen;

};
/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact44::updateCropResidueFluxes( void )
{
  int kdm;

  cropResidueFlux.carbon = ZERO;
  cropResidueFlux.nitrogen = ZERO;

  if( 1 == state )
  {
    for( kdm = 0; kdm < CYCLE; ++kdm )
    {
      if( harvstC[cmnt] > ZERO )
      {
        cropResidueFlux.carbon += initPROD1[kdm].carbon
                                  * (1.000000 - harvstC[cmnt])
                                  / harvstC[cmnt]
                                  * residueC[cmnt]
                                  / (double) CYCLE;
      }
      else { cropResidueFlux.carbon = ZERO; }

      if( harvstN[cmnt] > ZERO )
      {     
        cropResidueFlux.nitrogen += initPROD1[kdm].nitrogen
                                    * (1.000000 - harvstN[cmnt])
                                    / harvstN[cmnt]
                                    * residueN[cmnt]
                                    / (double) CYCLE;
      }
      else { cropResidueFlux.nitrogen = ZERO; }
    }
  }

};
/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact44::updateProducts( void )
{
  int i;

  // Carbon in Products

  PROD1.carbon = prevPROD1.carbon
                 + cropprod.carbon 
                 - PROD1decay.carbon;

  
  if( PROD1.carbon < 0.1 )
  {
    PROD1.carbon = ZERO;

    for(i = 0; i < CYCLE; ++i )
    {
      initPROD1[i].carbon = ZERO;
      initPROD1[i].nitrogen = ZERO;
    }
  }


  PROD10.carbon = prevPROD10.carbon
                  + formPROD10.carbon
                  - PROD10decay.carbon;
  
  if( PROD10.carbon < 0.1 )
  {
    PROD10.carbon = ZERO;

    for(i = 0; i < 10; ++i )
    {
      initPROD10[i].carbon = ZERO;
      initPROD10[i].nitrogen = ZERO;
    }
  }

  PROD100.carbon = prevPROD100.carbon
                   + formPROD100.carbon
                   - PROD100decay.carbon;
  
  if( PROD100.carbon < 0.1 )
  {
    PROD100.carbon = ZERO;

    for(i = 0; i < 100; ++i )
    {
      initPROD100[i].carbon = ZERO;
      initPROD100[i].nitrogen = ZERO;
    }
  }

  TOTPROD.carbon = PROD1.carbon
                   + PROD10.carbon
                   + PROD100.carbon;

  // Nitrogen in Products

  PROD1.nitrogen = prevPROD1.nitrogen
                   + cropprod.nitrogen
                   - PROD1decay.nitrogen;
  
  if( PROD1.nitrogen < 0.000001 )
  {
    PROD1.nitrogen = ZERO;
  }

  PROD10.nitrogen = prevPROD10.nitrogen
                    + formPROD10.nitrogen
                    - PROD10decay.nitrogen;
  
  if( PROD10.nitrogen < 0.000001 )
  {
    PROD10.nitrogen = ZERO;
  }

  PROD100.nitrogen = prevPROD100.nitrogen
                     + formPROD100.nitrogen
                     - PROD100decay.nitrogen;
  
  if( PROD100.nitrogen < 0.000001 )
  {
    PROD100.nitrogen = ZERO;
  }

  TOTPROD.nitrogen = PROD1.nitrogen
                     + PROD10.nitrogen
                     + PROD100.nitrogen;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact44::updateTotalProductFormation( void )
{

  // Carbon in Total Product Formation

  formTOTPROD.carbon = cropprod.carbon
                       + formPROD10.carbon
                       + formPROD100.carbon;
  
  // Nitrogen in Total Product Formation

  formTOTPROD.nitrogen = cropprod.nitrogen
                         + formPROD10.nitrogen
                         + formPROD100.nitrogen;
  
};

