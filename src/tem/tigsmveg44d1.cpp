/* *************************************************************
****************************************************************
TIGSMVEG44D1.CPP - object describing characteristics of vegetation 
                  used in the Terrestrial Ecosystem Model (TEM)

Modifications:

20060126 - DWK created by modifying tveg50b5.cpp 
20060126 - DWK changed include from tveg50b5.h to tveg437.h
20060126 - DWK changed Tveg50:: to Tveg43::
20060126 - DWK changed inheritance from ProcessXML50 to 
           ProcessXML43 in Tveg43()
20060126 - DWK deleted setThawPercent()
20060126 - DWK added double aot40 and int o3flag to 
           updateDynamics()
20060126 - DWK added double gppxio3() and double gppxo3()
20060126 - DWK added double fozone to double nupxclm()
20070201 - DWK changed include tveg438.h to tveg44.h
20070201 - DWK changed Tveg43:: to Tveg44::
20070201 - DWK and TWC added const int& agperennial to 
           updateDynamics()
20070201 - DWK added public function setRESPQ10LaRS()
20080130 - DWK changed include from tveg44.h to tigsmveg44a.h
20080130 - DWK changed ProcessXML43:: ProcessXML44::
20080826 - DWK changed calculation of cmax based on soil texture
           to cmax based on kc and soil texture in resetECDs()
20080826 - DWK changed include from tigsmveg44a.h to
           tigsmveg44b.h
20110707 - DWK changed include from tigsmveg44b.h to
           tigsmveg44c.h
20150429 - DWK changed include tigsmveg44c.h to tigsmveg44d1.h
           
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

#include "tigsmveg44d1.h"


Tveg44::Tveg44() : ProcessXML44()
{
  int dcmnt;

  plant.carbon = MISSING;
  plant.nitrogen = MISSING;

  strctrl.carbon = MISSING;
  strctrl.nitrogen = MISSING;

  labile.carbon = MISSING;
  labile.nitrogen = MISSING;

  ltrfal.carbon = MISSING;
  ltrfal.nitrogen = MISSING;

  nmobil = MISSING;
  nresorb = MISSING;

  inuptake = MISSING;

  nuptake = MISSING;

  suptake = MISSING;
  luptake = MISSING;

  ingpp = MISSING;
  gpp = MISSING;

  innpp = MISSING;
  npp = MISSING;

  rm = MISSING;
  rg = MISSING;
  gpr = MISSING;

  unnormleaf = MISSING;
  leaf = MISSING;

  lai = MISSING;
  fpc = MISSING;


  inprodcn = MISSING;

  yrcarbon = MISSING;
  yrnitrogen = MISSING;

  yrstructn = MISSING;
  yrc2n = MISSING;

  yrstoren = MISSING;

  yrltrfalc = MISSING;
  yrltrfaln = MISSING;

  yringpp = MISSING;
  yrgpp = MISSING;

  yrinnpp = MISSING;
  yrnpp = MISSING;

  yrgpr = MISSING;
  yrrmaint = MISSING;
  yrrgrowth = MISSING;

  yrinnup = MISSING;

  yrnup = MISSING;

  yrsup = MISSING;
  yrlup = MISSING;

  inprodcn = MISSING;

  yrnmobil = MISSING;
  yrnrsorb = MISSING;

  yrunleaf = MISSING;
  yrleaf = MISSING;

  yrfpc = MISSING;
  alleaf = MISSING;
  foliage = MISSING;

  leafyrs = 10;

  for( dcmnt = 0; dcmnt < MAXCMNT; ++dcmnt )
  {
    c2na[dcmnt] = MISSING;
    c2nb[dcmnt] = MISSING;
    c2nmin[dcmnt] = MISSING;
    cnmin[dcmnt] = MISSING;
    initcneven[dcmnt] = MISSING;

    aleaf[dcmnt] = MISSING;
    bleaf[dcmnt] = MISSING;
    cleaf[dcmnt] = MISSING;
    initleafmx[dcmnt] = MISSING;
    minleaf[dcmnt] = MISSING;
    unleaf12[dcmnt] = MISSING;

    cov[dcmnt] = MISSING;
    fpcmax[dcmnt] = MISSING;
    sla[dcmnt] = MISSING;

    kleafc[dcmnt] = MISSING;
    leafmxc[dcmnt] = MISSING;

    cmaxcut[dcmnt] = MISSING;
    cmax1a[dcmnt] = MISSING;
    cmax1b[dcmnt] = MISSING;
    cmax2a[dcmnt] = MISSING;
    cmax2b[dcmnt] = MISSING;

    kc[dcmnt] = MISSING;
    ki[dcmnt] = MISSING;

    tmax[dcmnt] = MISSING;
    tmin[dcmnt] = MISSING;
    toptmin[dcmnt] = MISSING;
    toptmax[dcmnt] = MISSING;

    gva[dcmnt] = MISSING;

    kra[dcmnt] = MISSING;
    krb[dcmnt] = MISSING;

    raq10a0[dcmnt] = MISSING;
    raq10a1[dcmnt] = MISSING;
    raq10a2[dcmnt] = MISSING;
    raq10a3[dcmnt] = MISSING;

    alpha[dcmnt] = MISSING;
    beta[dcmnt] = MISSING;
    gamma[dcmnt] = MISSING;
    qref[dcmnt] = MISSING;
    tref[dcmnt] = MISSING;

    nmaxcut[dcmnt] = MISSING;
    nmax1a[dcmnt] = MISSING;
    nmax1b[dcmnt] = MISSING;
    nmax2a[dcmnt] = MISSING;
    nmax2b[dcmnt] = MISSING;

    kn1[dcmnt] = MISSING;

    cfall[dcmnt] = MISSING;
    nfall[dcmnt] = MISSING;

    labncon[dcmnt] = MISSING;
  }

  cneven = MISSING;
  c2n = MISSING;
  dc2n = MISSING;
  adjc2n = MISSING;

  newleafmx = MISSING;
  prvleafmx = MISSING;

  cmax = MISSING;

  topt = MISSING;
  newtopt = MISSING;

  kr = MISSING;

  nmax = MISSING;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg44::boundTOPT( const int& pcmnt )
{
  if( topt > toptmax[pcmnt] ) {	topt = toptmax[pcmnt]; }

  if( topt < toptmin[pcmnt] ) { topt = toptmin[pcmnt]; }
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Tveg44::deltaleaf( const int& pdcmnt,
                          const double& eet,
                          const double& prveetmx,
                          const double& prvleaf )
{
  double maxeet;

  double normeet;

  double unnormleaf;

  if( prveetmx < 0.1 ) { maxeet = 1.0; }
  else { maxeet = prveetmx; }
  
  normeet = eet / maxeet;

  unnormleaf = (aleaf[pdcmnt] * normeet)
               + (bleaf[pdcmnt] * prvleaf)
               + cleaf[pdcmnt];

  if( unnormleaf < (0.5 * minleaf[pdcmnt]) )
  {
    unnormleaf = 0.5 * minleaf[pdcmnt];
  }

  if( unnormleaf > 1.000000000 )
  {
    unnormleaf = 1.000000000;
  }

  return unnormleaf;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg44::getecd( ofstream& rflog1 )
{

  string ecd;

  cout << "Enter name of the file with the vegetation";
  cout << " parameter values (.ECD):";
  cout << endl;

  cin >> ecd;

  rflog1 << "Enter name of the file with the vegetation";
  rflog1 << " parameter values (.ECD):";
  rflog1 << ecd << endl << endl;

  getecd( ecd );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg44::getecd( const string& ecd )
{
  int comtype;
  int dcmnt;
  ifstream infile;


  infile.open( ecd.c_str(), ios::in );

  if( !infile )
  {
    cerr << endl << "Cannot open " << ecd;
    cerr << " for vegetation ECD input" << endl;
    
    exit( -1 );
  }

  getXMLrootNode( infile, "vegECD" );

  for( dcmnt = 1; dcmnt < MAXCMNT; ++dcmnt )
  {
    comtype = getXMLcommunityNode( infile, "vegECD" );

    if( comtype >= MAXCMNT )
    {
      cerr << endl << "comtype is >= MAXCMNT" << endl;
      cerr << "comtype cannot be greater than ";
      cerr << (MAXCMNT-1) << endl;
      cerr << " in vegECD" << endl;
      
      exit( -1 );
    }

    kc[comtype]= getXMLcmntArrayDouble( infile,
                                        "vegECD",
                                        "kc",
                                        comtype );

    ki[comtype]= getXMLcmntArrayDouble( infile,
                                        "vegECD",
                                        "ki",
                                        comtype );

    gva[comtype]= getXMLcmntArrayDouble( infile,
                                         "vegECD",
                                         "gva",
                                         comtype );

    tmin[comtype]= getXMLcmntArrayDouble( infile,
                                          "vegECD",
                                          "tmin",
                                          comtype );

    toptmin[comtype]= getXMLcmntArrayDouble( infile,
                                             "vegECD",
                                             "toptmin",
                                             comtype );

    toptmax[comtype]= getXMLcmntArrayDouble( infile,
                                             "vegECD",
                                             "toptmax",
                                             comtype );

    tmax[comtype]= getXMLcmntArrayDouble( infile,
                                          "vegECD",
                                          "tmax",
                                          comtype );

    raq10a0[comtype]= getXMLcmntArrayDouble( infile,
                                             "vegECD",
                                             "raq10a0",
                                             comtype );

    raq10a1[comtype]= getXMLcmntArrayDouble( infile,
                                             "vegECD",
                                             "raq10a1",
                                             comtype );

    raq10a2[comtype]= getXMLcmntArrayDouble( infile,
                                             "vegECD",
                                             "raq10a2",
                                             comtype );

    raq10a3[comtype]= getXMLcmntArrayDouble( infile,
                                             "vegECD",
                                             "raq10a3",
                                             comtype );

//    tref[comtype]= getXMLcmntArrayDouble( infile,
//                                          "vegECD",
//                                          "tref",
//                                          comtype );

//    qref[comtype]= getXMLcmntArrayDouble( infile,
//                                          "vegECD",
//                                          "qref",
//                                          comtype );

//    alpha[comtype]= getXMLcmntArrayDouble( infile,
//                                           "vegECD",
//                                           "alpha",
//                                           comtype );

//    beta[comtype]= getXMLcmntArrayDouble( infile,
//                                          "vegECD",
//                                          "beta",
//                                          comtype );

//    gamma[comtype]= getXMLcmntArrayDouble( infile,
//                                           "vegECD",
//                                           "gamma",
//                                           comtype );

    kn1[comtype]= getXMLcmntArrayDouble( infile,
                                         "vegECD",
                                         "kn1",
                                         comtype );

    labncon[comtype]= getXMLcmntArrayDouble( infile,
                                             "vegECD",
                                             "labncon",
                                             comtype );

    leafmxc[comtype]= getXMLcmntArrayDouble( infile,
                                             "vegECD",
                                             "leafmxc",
                                             comtype );

    kleafc[comtype]= getXMLcmntArrayDouble( infile,
                                            "vegECD",
                                            "kleafc",
                                            comtype );

    sla[comtype]= getXMLcmntArrayDouble( infile,
                                         "vegECD",
                                         "sla",
                                         comtype );

    cov[comtype]= getXMLcmntArrayDouble( infile,
                                         "vegECD",
                                         "cov",
                                         comtype );

    fpcmax[comtype]= getXMLcmntArrayDouble( infile,
                                            "vegECD",
                                            "fpcmax",
                                            comtype );

    endXMLcommunityNode( infile );
  }

  if( dcmnt < MAXCMNT )
  {
    cerr << endl << " Parameters found for only " << dcmnt;
    cerr << " community types out of a maximum of ";
    cerr << (MAXCMNT-1) << " types in vegECD" << endl;
    
    exit( -1 );
  }

  infile.close();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg44::getleafecd( ofstream& rflog1 )
{

  string ecd;

  cout << "Enter name of the file with leaf parameter";
  cout << " values (.ECD):" << endl;

  cin >> ecd;

  rflog1 << "Enter name of the file with leaf parameter";
  rflog1 << " values (.ECD):";
  rflog1 << ecd << endl << endl;

  getleafecd( ecd );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg44::getleafecd( const string& ecd )
{
  int comtype;
  int dcmnt;
  ifstream infile;


  infile.open( ecd.c_str(), ios::in );

  if( !infile )
  {
    cerr << endl << "Cannot open " << ecd;
    cerr << " for leaf ECD input" << endl;
    
    exit( -1 );
  }

  getXMLrootNode( infile, "leafECD" );

  for( dcmnt = 1; dcmnt < MAXCMNT; ++dcmnt )
  {
    comtype = getXMLcommunityNode( infile, "leafECD" );

    if( comtype >= MAXCMNT )
    {
      cerr << endl << "comtype is >= MAXCMNT" << endl;
      cerr << "comtype cannot be greater than " << (MAXCMNT-1);
      cerr << " in leafECD" << endl;
      
      exit( -1 );
    }

    minleaf[comtype]= getXMLcmntArrayDouble( infile,
                                             "leafECD",
                                             "minleaf",
                                             comtype );

    aleaf[comtype]= getXMLcmntArrayDouble( infile,
                                           "leafECD",
                                           "aleaf",
                                           comtype );

    bleaf[comtype]= getXMLcmntArrayDouble( infile,
                                           "leafECD",
                                           "bleaf",
                                           comtype );

    cleaf[comtype]= getXMLcmntArrayDouble( infile,
                                           "leafECD",
                                           "cleaf",
                                           comtype );

    endXMLcommunityNode( infile );
  }

  if( dcmnt < MAXCMNT )
  {
    cerr << endl << " Parameters found for only " << dcmnt;
    cerr << " community types out of a maximum of ";
    cerr << (MAXCMNT-1) << " types in leafECD" << endl;
    
    exit( -1 );
  }

  infile.close();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Tveg44::gppxclm( const int& pdcmnt,
                        const double& co2,
                        const double& par,
                        const double& temp,
                        const double& gv,
                        const double& leaf,
                        const double& foliage )
{

  double gpp;


/* *************************************************************
   gpp:    gpp as influenced by carbon dioxide (co2), moisture
	   (gv), phenology (leaf), photosynthetically active
	   radiation (par), air temperature (temp) and freeze-
	   thaw index (thawpercent)

************************************************************* */

  gpp  = co2 * gv;

  gpp *= cmax * foliage / (kc[pdcmnt] + gpp);

  gpp *= leaf * par / (ki[pdcmnt] + par);

  gpp *= temp;


  return gpp;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Tveg44::gppxio3( const double& fozone,
                        const double& eetpet )
{
  double findozone;

  if( ZERO == fozone )
  {
    findozone = 1.0;
  }
  else
  {
    findozone = (1.0-1.0/fozone) * eetpet + (1.0/fozone);
  }

//  findozone = ((1.0/fozone)-1) * exp(-5.0*fprevozone) + 1.0;

  return findozone;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Tveg44::gppxo3( const int& pdcmnt,
                       const double& gpp,
                       const double& d40,
                       const double& eetpet )
{
  double fozone;

  fozone = 1.0 - (o3para[pdcmnt] * pow( 10.0, -6 )
           * (o3parb[pdcmnt] + o3parc[pdcmnt] * gpp) * d40);

  fozone = 1 - eetpet + (eetpet*fozone);

  fozone = fozone + fprevozone - 1.0;


  // Keep fozone between 0.0 and 1.0
  
  if( fozone <= ZERO )
  {
    fozone = ZERO;
  }
  
  if( fozone >= 1.0 )
  { 
    fozone = 1.000000;
  }


  return fozone;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg44::leafinit( ofstream& rflog1 )
{

  int lfswtch;

  cout << "Enter number of years for model run:  " << endl;

  cin >> leafyrs;

  cout << "Do you have a file containing the phenology parameters?:";
  cout << endl;
  cout << "Enter 0 for no:" << endl;
  cout << "Enter 1 for yes:" << endl;

  cin >> lfswtch;

  rflog1 << "Enter number of years for model run:  ";
  rflog1 << leafyrs << endl;
  rflog1 << "Do you have a file containing the phenology parameters?:";
  rflog1 << endl;
  rflog1 << "Enter 0 for no:" << endl;
  rflog1 << "Enter 1 for yes:" << endl;
  rflog1 << lfswtch << endl;

  if( 0 == lfswtch )
  {
    cout << "Enter regression coefficient for relative EET";
    cout << " (i.e. 'a'):  ";

    cin >> aleaf[0];

    cout << "Enter regression coefficient for LAI(t-1)";
    cout << " (i.e. 'b'):";

    cin >> bleaf[0];

    cout << "Enter regression intercept (i.e. 'c'):  ";

    cin >> cleaf[0];

    rflog1 << "Enter regression coefficient for relative EET";
    rflog1 << " (i.e. 'a'): ";
    rflog1 <<  aleaf[0] << endl;
    rflog1 << "Enter regression coefficient for LAI(t-1)";
    rflog1 << " (i.e. 'b'):  ";
    rflog1 << bleaf[0] << endl;
    rflog1 << "Enter regression intercept (i.e. 'c'):  ";
    rflog1 << cleaf[0] << endl;

    minleaf[0] = 2.0;
    
    while( minleaf[0] >= 1.0 )
    {
      cout << "Enter minimum LAI (must be less than 1.0):";

      cin >> minleaf[0];

      rflog1 << "Enter minimum LAI (must be less than 1.0): ";
      rflog1 << minleaf[0];
      rflog1 << endl << endl;
    }
  }
  else { getleafecd( rflog1 ); }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Tveg44::nupxclm( const int& pdcmnt,
                        const double& soilh2o,
                        const double& availn,
                        const double& respq10,
                        const double& ksoil,
                        const double& foliage,
                        const double& fozone )
{

/* *************************************************************
   nuptake:  uptake of nitrogen by plants as influenced by
	     available nitrogen concentration (availn), moisture
	     (ksoil), ozone, and air temperature (respq10)
************************************************************* */

  double nuptake;

  nuptake  = (availn * ksoil) / soilh2o;

  nuptake *= nmax * foliage / (kn1[pdcmnt] + nuptake);

  nuptake *= respq10;

  nuptake *= fozone;  

  return nuptake;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void Tveg44::resetEcds( const int& pcmnt,
                        const double& psiplusc )
{
  // Initialize TEM parameters dependent upon a grid cell's
  //   soil texture

//  if( psiplusc <= cmaxcut[pcmnt] )
//  {
//    cmax = (cmax1a[pcmnt] * psiplusc) + cmax1b[pcmnt];
//  }
//  else
//  {
//    cmax = (cmax2a[pcmnt] * psiplusc) + cmax2b[pcmnt];
//  }

  // Use the next two lines instead for uncertainty analyses 
  //   with kc

  cmax = ((cmax1a[pcmnt] * psiplusc) + cmax1b[pcmnt]) * kc[pcmnt];

  cmax += (cmax2a[pcmnt] * psiplusc) + cmax2b[pcmnt];

  if( cmax < ZERO ) { cmax = ZERO; }

//  cout << "CMAX = " << cmax;
//  cout << " PCMNT = " << pcmnt;
//  cout << " PSIPLUSC = " << psiplusc << endl;
  
  if( psiplusc <= nmaxcut[pcmnt] )
  {
    nmax = (nmax1a[pcmnt] * psiplusc) + nmax1b[pcmnt];
  }
  else
  {
    nmax = (nmax2a[pcmnt] * psiplusc) + nmax2b[pcmnt];
  }

  if( nmax < ZERO ) { nmax = ZERO; }
  
  
};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void Tveg44::resetMonthlyFluxes( void )
{  
  leaf = ZERO;
  lai = ZERO;
  fpc = ZERO;

  // Reset monthly fluxes to zero
  
  ingpp = ZERO;
  gpp = ZERO;
  fozone = ZERO;
  findozone = ZERO;

  innpp = ZERO;
  npp = ZERO;
  
  gpr = ZERO;
  rm = ZERO;
  rg = ZERO;
  
  ltrfal.carbon = ZERO;

  inuptake = ZERO;

  nuptake = ZERO;

  suptake = ZERO;
  luptake = ZERO;
  nmobil = ZERO;
  nresorb = ZERO;

  ltrfal.nitrogen = ZERO;
	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg44::resetNEWTOPT( const int& pcmnt, 
                           const double& tair, 
                           const double& unnrmleaf )
{
  if( ZERO == aleaf[pcmnt] 
      && ZERO == bleaf[pcmnt]
      && 1.0 == cleaf[pcmnt] )
  {
    if( tair > newtopt ) { newtopt = tair; }
  }
  else
  {
    if( unnrmleaf > newleafmx )
    {
      newleafmx = unnrmleaf;

      newtopt = tair;
    }
  }
	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg44::resetYrFluxes( void )
{
  // Reset annual fluxes and summary variables to zero

  yrcarbon = ZERO;

  yrstructn = ZERO;

  yrstoren = ZERO;

  yrnitrogen = ZERO;

  yrc2n = ZERO;

  // Phenology

  yrunleaf = ZERO;

  yrleaf = ZERO;

  yrlai = ZERO;

  yrfpc = ZERO;

  // Carbon fluxes
  
  yringpp = ZERO;

  yrgpp = ZERO;

  yrinnpp = ZERO;

  yrnpp = ZERO;

  yrgpr = ZERO;

  yrrmaint = ZERO;

  yrrgrowth = ZERO;

  yrltrfalc = ZERO;

  // Nitrogen fluxes
  
  yrinnup = ZERO;

  yrnup = ZERO;

  yrsup = ZERO;

  yrlup = ZERO;

  yrnmobil = ZERO;

  yrnrsorb = ZERO;

  yrltrfaln = ZERO;

};
/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Tveg44::rmxclm( const int& pdcmnt,
                       const double& vegc,
                       const double& respq10 )
{

/* *************************************************************
   rm: plant maintenance respiration as influenced by plant
       biomass (vegc) and air temperature (respq10)
************************************************************* */

  double rmaint;

  kr = exp( (kra[pdcmnt]*vegc) + krb[pdcmnt] );

  rmaint  = kr * vegc;

  rmaint *= respq10;

  return rmaint;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Tveg44::setGV( const double& eet,
                      const double& pet,
                      const int& moistlim )
{
  double tstgv;

/* tstgv: effect of moisture on primary productivity */

  if( 0 == moistlim )
  {
    tstgv = 1.0;
  }
  else
  {
    if( ZERO == pet  ) { tstgv = ZERO; }
    else if( eet / pet <= 0.1 )
    {
      tstgv = (-10.0 * pow( (eet / pet), 2.0 ))
              + (2.9 * (eet / pet));
      
      if( tstgv < ZERO ) { tstgv = ZERO; }
    }
    else
    {
      tstgv = 0.1 + (0.9 * eet / pet);
    }
  }

  return tstgv;
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg44::setRESPQ10( const int& pdcmnt, const double& tair )
{
/* *************************************************************
 respq10: effect of temperature on plant respiration
************************************************************* */

  double raq10;

  raq10 = raq10a0[pdcmnt] + (raq10a1[pdcmnt] * tair)
          + (raq10a2[pdcmnt]*pow( tair, 2.0 ))
          + (raq10a3[pdcmnt]*pow( tair, 3.0 ));
  
  respq10 = pow( raq10, tair/10.0 );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg44::setRESPQ10LaRS( const int& pdcmnt, 
                             const double& tair )
{
/* *************************************************************
 respq10: effect of temperature on plant respiration

 Amthor, J. (in review)
************************************************************* */

  double raq10;

  raq10 = qref[pdcmnt]
          * exp( -1.0 * alpha[pdcmnt] * (tair - tref[pdcmnt]) );

  respq10 = pow( raq10, (tair - tref[pdcmnt])/10.0 )
            / ( 1.0 + exp( beta[pdcmnt] - tair )
                + exp( tair - gamma[pdcmnt] ));

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg44::setTEMP( const int& pdcmnt, const double& tair )
{
  /* temp: effect of temperature on primary productivity */

  if( tair <= tmin[pdcmnt] || tair >= tmax[pdcmnt] )
  {
    temp = ZERO;
  }
  else
  {
    if( tair >= topt && tair <= toptmax[pdcmnt] )
    {
      temp = 1.0;
    }
    else
    {
      if( tair > tmin[pdcmnt] && tair < topt )
      {
	temp = (tair - tmin[pdcmnt])
               * (tair - tmax[pdcmnt])
               / ((tair - tmin[pdcmnt])
               * (tair - tmax[pdcmnt])
               - pow( (tair - topt), 2.0 ));
      }
      else
      {
	temp = (tair - tmin[pdcmnt])
               * (tair - tmax[pdcmnt])
               / ((tair - tmin[pdcmnt])
               * (tair - tmax[pdcmnt])
               - pow( (tair - toptmax[pdcmnt]), 2.0 ));
      }
    }
  }
	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg44::showecd( const int& pdcmnt )
{

  cout << endl << "             VEGETATION PARAMETERS INFLUENCED BY CLIMATE";
  cout << endl << endl;

  printf( "  KI = %6.2lf  KC = %6.2lf  KN1 = %6.4lf  GVA = %8.4lf\n\n",
          ki[pdcmnt],
          kc[pdcmnt],
          kn1[pdcmnt],
          gva[pdcmnt] );

  printf( "   TMIN = %5.1lf    TOPTMIN = %5.1lf   TOPTMAX = %5.1lf   TMAX = %5.1lf\n\n",
          tmin[pdcmnt],
          toptmin[pdcmnt],
          toptmax[pdcmnt],
          tmax[pdcmnt] );

  printf( " RAQ10A0 = %7.5lf  RAQ10A1 = %7.5lf  RAQ10A2 = %7.5lf  RAQ10A3 = %7.5lf\n",
          raq10a0[pdcmnt],
          raq10a1[pdcmnt],
          raq10a2[pdcmnt],
          raq10a3[pdcmnt] );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg44::showleaf( const int& pdcmnt )
{

  cout << endl << "         PARAMETERS FOR THE LEAF PHENOLOGY MODEL";
  cout << endl << endl;

  printf( "     ALEAF = %7.5lf     BLEAFC = %8.5lf        CLEAF = %8.5lf\n",
          aleaf[pdcmnt],
          bleaf[pdcmnt],
          cleaf[pdcmnt] );

  printf( "   MINLEAF = %4.2lf       MAXLEAF = %7.4lf      UNLEAF12 = %7.4lf\n",
          minleaf[pdcmnt],
          initleafmx[pdcmnt],
          unleaf12[pdcmnt] );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg44::updateC2N( const int& pdcmnt,
                        const double& yreet,
                        const double& yrpet,
                        const double& currentco2,
                        const double& initco2 )
{
  if( yrpet != ZERO )
  {
    c2n = c2nb[pdcmnt] + c2na[pdcmnt]*(yreet/yrpet);
  }
  else { c2n = c2nb[pdcmnt]; }

  if( yreet > yrpet ) 
  {
  cout  << "c2n=" << c2n << " yreet=" << yreet << " yrpet=" << yrpet << endl;
  }

  if( c2n < c2nmin[pdcmnt] ) { c2n = c2nmin[pdcmnt]; }
  
  adjc2n = 1.0 + (dc2n * (currentco2 - initco2));

  c2n *= adjc2n;

/*cout  << "c2n=" << c2n << " yreet=" << yreet << " yrpet=" << yrpet << endl;
  cout  << "dc2n=" << dc2n << " currentco2=" << currentco2 << " initco2=" << initco2 << endl; */

  cneven = initcneven[pdcmnt] * adjc2n;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg44::updateDynamics( const int& pdcmnt,
                             const double& co2,
                             const double& aot40,
                             const double& ninput,
                             const double& par,
                             const double& pet,
                             const double& prevmaxpet,
                             const double& eet,
                             const double& prevmaxeet,
                             const double& vegc,
                             const double& structn,
                             const double& labilen,
                             const double& soilh2o,
                             const double& availn,
                             const int& moistlim,
                             const int& nfeed,
                             const int& o3flag,
                             const int& agstate,
                             const int& agperennial,
                             const int& fertflag,
                             const double& ksoil,
                             const double& netnmin,
                             double& agfertn )
{	
  double eetpet;

  double gv;
  
  // Determine stage of seasonal phenology

  if( 0 == moistlim )
  {
    unnormleaf = deltaleaf( pdcmnt,
                            pet,
                            prevmaxpet,
                            prevunrmleaf );
  }
  else
  {
    unnormleaf = deltaleaf( pdcmnt,
                            eet,
                            prevmaxeet,
                            prevunrmleaf );                           
  }

  if( prvleafmx <= ZERO ) { leaf = ZERO; }
  else 
  { 
    if( prvleafmx > 1.000000000 )
    {
      prvleafmx = 1.000000000;
    }

    leaf = unnormleaf / prvleafmx; 
  }

  if( leaf < minleaf[pdcmnt] )
  {
    leaf = minleaf[pdcmnt];
  }

  if( leaf > 1.0 ) { leaf = 1.0; }

  alleaf = leafmxc[pdcmnt]
           /(1.0 + kleafc[pdcmnt] * exp( cov[pdcmnt] * vegc ));

  lai = sla[pdcmnt] * alleaf * leaf;

  foliage = alleaf / leafmxc[pdcmnt];

  fpc = 1.0 - exp( -0.5 * lai );


  // gv: effect of moisture on primary productivity

  gv = setGV( eet, pet, moistlim );


  // Determine Gross Primary Production if nitrogen is not
  //   limiting (veg.ingpp)

  ingpp = gppxclm( pdcmnt,
                   co2,
                   par,
                   temp,
                   gv,
                   leaf,
                   foliage );


  if( ingpp < ZERO ) { ingpp = ZERO; }
 

  /* ozone: effect of surface ozone on gpp */

  if( 0 == o3flag )
  {
    fozone = 1.000000;
    fprevozone = 1.000000;
    findozone = 1.000000;
  }
  else
  {
    if( ZERO == pet )
    {
      eetpet = ZERO;
    }
    else
    {
      eetpet = eet / pet;
    }

    fozone = gppxo3( pdcmnt,
                     ingpp,
                     aot40,
                     eetpet );

    ingpp *= fozone;

//    findozone = gppxio3( fozone, eetpet );
//    ingpp *= findozone;
  }


  if( ingpp < ZERO ) { ingpp = ZERO; }


  // Determine Maintenance Respiration of vegetation (veg.rm)

  rm = rmxclm( pdcmnt, vegc, respq10 );

  if( rm < ZERO ) { rm = ZERO; }


  // Determine Net Primary Production if nitrogen is not
  //   limiting (veg.innpp) and Growth Respiration (veg.rg)

  innpp = ingpp - rm;

  rg = ZERO;

  if( innpp > ZERO )
  {
    rg  = 0.2 * innpp;

    innpp *= 0.8;
  }


  // Determine monthly loss of carbon and nitrogen in litterfall

  ltrfal.carbon = cfall[pdcmnt] * vegc;

  if( ltrfal.carbon < ZERO )
  {
    ltrfal.carbon = ZERO;
  }

  ltrfal.nitrogen = nfall[pdcmnt] * structn;
  
  if( ltrfal.nitrogen < ZERO )
  {
    ltrfal.nitrogen = ZERO;
  }


  // Determine nitrogen uptake by vegetation if carbon is not
  //   limiting (veg.inuptake)

  inuptake = nupxclm( pdcmnt,
                      soilh2o,
                      availn,
                      respq10,
                      ksoil,
                      foliage,
                      fozone );



  if( inuptake > availn + netnmin + ninput )
  {
    inuptake = availn + netnmin + ninput;
  }

  if( inuptake < ZERO )
  {
    inuptake = ZERO;
  }


  // Determine how interactions between carbon and nitrogen
  //   availability influence primary production, litterfall
  //   and nitrogen uptake

  gpp = ingpp;

  npp = innpp;


  // Assume CROP NPP is always positive or zero

  if( 1 == agstate && 0 == agperennial )
  {
     if( npp < ZERO ) { npp = ZERO; }
  }
 

  nuptake = inuptake;

  suptake = nuptake;
  luptake = ZERO;

  nmobil = ZERO;
  nresorb = ZERO;


  // Nitrogen feedback of GPP ( 1 == nfeed)

  if( 1 == nfeed )
  {
    // Determine nitrogen resorption (veg.nresorb)
    
    if( ltrfal.nitrogen
         <= ltrfal.carbon / cneven )
    {
      nresorb = ltrfal.carbon / cneven - ltrfal.nitrogen;
    }
    else
    {
      nresorb = ZERO;

      ltrfal.nitrogen = ltrfal.carbon / cneven;
    }

    if( vegc > ZERO )
    {
      nresorb *= (structn / vegc) * c2n;
    }


    // Determine if monthly primary production is carbon-limited  
    //   or nitrogen-limited
    
    if( ZERO == (nuptake + labilen) )
    {
      inprodcn = innpp / 0.000001;
    }
    else 
    {
      inprodcn = innpp / (nuptake + labilen);
    }
    
 
    // If primary production is nitrogen-limited, 
    //   (i.e., veg.inprodcn > veg.cneven) revaluate NPP, RG and 
    //   GPP based on nitrogen availability
     
    if( inprodcn > cneven )
    {
      if( 1 == fertflag )
      {
        // Assume nitrogen is never limiting if fertilized
        // Also assume that fertilized crops are based solely 
        //   on a nitrate economy
          
        nuptake = (npp / cneven) - labilen;

        if( nuptake > inuptake )
        {
          agfertn = nuptake - inuptake;
        }
      }
      else
      {
        npp = cneven * (nuptake + labilen);

        if( npp < ZERO ) { npp = ZERO; }

        rg = 0.25 * npp;

        gpp = npp + rg + rm;

        if( gpp < ZERO ) { gpp = ZERO; }

        nmobil = labilen;
      }
    }


    // If primary production is carbon-limited, 
    //   (i.e., veg.inprodcn < veg.cneven) revaluate nitrogen
    //   uptake by vegetation based on carbon availability

    if( inprodcn <= cneven )
    {
      nuptake = nuptake
                * (inprodcn - cnmin[pdcmnt])
                * (inprodcn - 2 * cneven + cnmin[pdcmnt]);

      nuptake /= ((inprodcn - cnmin[pdcmnt])
                 * (inprodcn - 2 * cneven + cnmin[pdcmnt]))
                 - pow( inprodcn - cneven, 2.0 );


      if( nuptake < ZERO )
      {
        nuptake = ZERO;
      }


   // Drawdown N from labile N pool first to satisfy N
      //   requirement before taking up any nitrogen from the 
      //   soil
      
      if( labilen >= (npp / cneven) )
      {
        nmobil = npp / cneven;
	
	      if( nmobil < ZERO && vegc > ZERO )
        {
	        nmobil *= (structn / vegc) * c2n;
	      }
	
	      suptake = ZERO;
      }
      else
      {
        nmobil = labilen;
	
	      suptake = (npp / cneven) - nmobil;

	      if( suptake < ZERO ) 
 	      { 
	        suptake = ZERO; 
	      }
	
	      if( suptake > nuptake )
        {
          suptake = nuptake;
        }
      }


      // Determine vegetation nitrogen uptake for the labile
      //   N pool
      
      if( (labilen + nuptake - suptake + nresorb - nmobil) 
           < (labncon[pdcmnt] * (structn + suptake 
           - ltrfal.nitrogen - nresorb + nmobil)) )
      {
        luptake = nuptake - suptake;
      }
      else
      {
        luptake = (labncon[pdcmnt] 
                  * (structn + suptake - ltrfal.nitrogen
                  - nresorb + nmobil))  
                  - (labilen + nresorb - nmobil);

	      if( luptake < ZERO ) { luptake = ZERO; }

	      nuptake = suptake + luptake;
      }
    }
  }


  // Determine Gross Plant Respiration (GPR)

  gpr = rm + rg;

};
