/* *************************************************************
****************************************************************
HUMANACT438.CPP - describes human disturbances to natural 
                     ecosystems

Modifications:

20060126 - DWK created by modifying humnact50b5.cpp
20060126 - DWK changed include from humnact50b5.h to 
           humnact437.h
20060126 - DWK changed Humanact50:: to Humanact43::
20060126 - DWK changed inheritance of ProcessXML50 to
           ProcessXML43 in Humanact43()
                      
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

#include "humnact438.h"


/* *************************************************************
************************************************************* */

Humanact43::Humanact43() : ProcessXML43()
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

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact43::conversion( const int& pdcmnt,
                             const double& vegc,
                             const double& vstrn,
                             const double& vston,
                             const double& soilc,
                             const double& soiln )
{
  // Assume annual conversion fluxes are partitioned equally 
  //   across each of 12 months - determine monthly fluxes
  
  slash.carbon = slashpar * vegc / (double) CYCLE;
  slash.nitrogen = slashpar  * (vstrn + vston) / (double) CYCLE;

  vconvrtflx.carbon = (vconvert * vegc) / (double) CYCLE;
  sconvrtflx.carbon =  (sconvert * soilc)/ (double) CYCLE;

  convrtflx.carbon = vconvrtflx.carbon + sconvrtflx.carbon;

  vconvrtflx.nitrogen = ((1.0 - nvretconv[pdcmnt]) 
                        * vconvert 
                        * (vstrn + vston))
                        / (double) CYCLE;
                                
  sconvrtflx.nitrogen = ((1.0 - nsretconv[pdcmnt]) 
                        * sconvert * soiln)
                        / (double) CYCLE;
                                
  convrtflx.nitrogen = vconvrtflx.nitrogen + sconvrtflx.nitrogen;
                               
  nvretent = (nvretconv[pdcmnt] * vconvert * (vstrn + vston)) 
              / (double) CYCLE;

  nsretent = (nsretconv[pdcmnt] * sconvert * soiln) 
              / (double) CYCLE;
  
  nretent = nvretent + nsretent;
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact43::createWoodProducts( const int& pdyr,
                                     const double& vegc,
                                     const double& vstrn,
                                     const double& vston )
{
  int i;

  formPROD10.carbon  = prod10par * vegc;
  formPROD10.nitrogen  = prod10par * (vstrn + vston);
  formPROD100.carbon = prod100par * vegc;
  formPROD100.nitrogen = prod100par * (vstrn + vston);

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

void Humanact43::decayProducts( void )
{

  int i;
  int dm;
  double yrtrash;

  
  PROD1decay.carbon  = ZERO;
  PROD1decay.nitrogen  = ZERO;

  for( dm = 0; dm < CYCLE; ++dm )
  {
    PROD1decay.carbon = initPROD1[dm].carbon / (double) CYCLE;
                        
    PROD1decay.nitrogen = initPROD1[dm].nitrogen 
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

void Humanact43::frostDamage( const double& vegc,
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

void Humanact43::getecd( ofstream& rflog1 )
{
  string ecd;

  cout << "Enter name of the data file (.ECD) with agricultural ";
  cout << "parameter values: " << endl;
  cin >> ecd;

  rflog1 << "Enter name of the soil data file (.ECD) with agricultural ";
  rflog1 << "parameter values: " << ecd << endl;

  getecd( ecd );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact43::getecd( const string& ecd )
{
  ifstream infile;
  int dcmnt;
  int comtype;

  infile.open( ecd.c_str(), ios::in );

  if( !infile )
  {
    cerr << "\nCannot open " << ecd;\
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
    

    nvretconv[comtype] = getXMLcmntArrayDouble( infile,
                                                "agECD",
                                                "nvretconv",
                                                comtype );

    nsretconv[comtype] = getXMLcmntArrayDouble( infile,
                                                "agECD",
                                                "nsretconv",
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

void Humanact43::harvest( const int& pdm,
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

  if ( cropprod.carbon < ZERO )
  {
    cropprod.carbon = ZERO;
  }

  cropprod.nitrogen = (vstrn+vston) * harvstN[cmnt];

  if ( cropprod.nitrogen < ZERO )
  {
    cropprod.nitrogen = ZERO;
  }

  initPROD1[pdm].carbon  = cropprod.carbon;
  initPROD1[pdm].nitrogen  = cropprod.nitrogen;
  
  // Determine amount of carbon and nitrogen left
  // in crop residue

   initresidueC = vegc - cropprod.carbon;

  if ( initresidueC < ZERO )
  {
    initresidueC = ZERO;
  }

  initresidueN = (vstrn+vston) - cropprod.nitrogen;

  if ( initresidueN < ZERO )
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

void Humanact43::resetMonthlyFluxes( const int& pdm )
{
  // Reset monthly fluxes to zero
  
  fertn = ZERO;
  
  irrigate = ZERO;

  // Initialize stubble to zero
  
  stubble.carbon = ZERO;
  stubble.nitrogen = ZERO;

  cflux = ZERO;

  if( 0 == pdm )
  {
    // Initialize carbon fluxes during conversion to zero

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
  }
  
  // Initialize carbon fluxes related to formation and 
  //   decomposition of products and crop residue to zero

  cropprod.carbon = ZERO;
  cropprod.nitrogen = ZERO;

  formCropResidue.carbon = ZERO;
  formCropResidue.nitrogen = ZERO;

  PROD1decay.carbon = ZERO;
  PROD1decay.nitrogen = ZERO;

  formPROD10.carbon = ZERO;
  formPROD10.nitrogen = ZERO;

  PROD10decay.carbon = ZERO;
  PROD10decay.nitrogen = ZERO;

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

void Humanact43::resetPROD( void )
{
  int i;
  int dm;

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

void Humanact43::resetYrFluxes( void )
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

void Humanact43::setAgricFlags( ofstream& rflog1 )
{

  cout << "Are agricultural soils tilled?" << endl;
  cout << "Enter 0 for no:" << endl;
  cout << "Enter 1 for yes:" << endl;
  
  cin >> tillflag;

  rflog1 << "Are agricultural soils tilled?" << endl;
  rflog1 << "Enter 0 for no:" << endl;
  rflog1 << "Enter 1 for yes:" << endl;
  rflog1 << tillflag << endl << endl;

  cout << "Do you want to irrigate to avoid moisture limitation of crop productivity?" << endl;
  cout << "Enter 0 for no:" << endl;
  cout << "Enter 1 for yes:" << endl;
  
  cin >> irrgflag;

  rflog1 << "Do you want to irrigate to avoid moisture limitation of crop productivity?" << endl;
  rflog1 << "Enter 0 for no:" << endl;
  rflog1 << "Enter 1 for yes:" << endl;
  rflog1 << irrgflag << endl << endl;

  cout << "Do you want to fertilize to avoid nitrogen limitation of crop productivity?" << endl;
  cout << "Enter 0 for no:" << endl;
  cout << "Enter 1 for yes:" << endl;
  
  cin >> fertflag;

  rflog1 << "Do you want to fertilize to avoid nitrogen limitation of crop productivity?" << endl;
  rflog1 << "Enter 0 for no:" << endl;
  rflog1 << "Enter 1 for yes:" << endl;
  rflog1 << fertflag << endl << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact43::setNoCrops( const int& pdm )
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

void Humanact43::setNoWoodProducts( const int& pdyr )
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

void Humanact43::updateCropResidue( void )
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

void Humanact43::updateCropResidueFluxes( void )
{
  int kdm;

  cropResidueFlux.carbon = ZERO;
  cropResidueFlux.nitrogen = ZERO;

  if ( 1 == state )
  {
    for( kdm = 0; kdm < CYCLE; ++kdm )
    {
      cropResidueFlux.carbon += initPROD1[kdm].carbon
                                * (1.000000 - harvstC[cmnt])
                                / harvstC[cmnt]
                                * residueC[cmnt]
                                / (double) CYCLE;
      
      cropResidueFlux.nitrogen += initPROD1[kdm].nitrogen
                                  * (1.000000 - harvstN[cmnt])
                                  / harvstN[cmnt]
                                  * residueN[cmnt]
                                  / (double) CYCLE;
    }
  }

};
/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact43::updateProducts( void )
{
  int i;

  // Carbon in Products

  PROD1.carbon = prevPROD1.carbon 
                 + cropprod.carbon 
                 - PROD1decay.carbon;
  
  if( PROD1.carbon < 0.1 )
  {
    PROD1.carbon = ZERO;
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

void Humanact43::updateTotalProductFormation( void )
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

