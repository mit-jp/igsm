/* *************************************************************
****************************************************************
MDMSOILIGSM44C.CPP - object describing general characteristics
                     of soil for the Methane Dynamics Model (MDM)

Modifications:

20060125 - DWK created by combining and modifying 
           methanediff.cpp, methaneplant.cpp and methaneebu.cpp 
20060125 - DWK added public function void oxygenC() from 
           methaneoxi.cpp
20060125 - DWK added public function double RedoxP() from
           methanepro.cpp 
20060422 - DWK renamed interpolateSM() to be interpolateLayers() 
           and deleted interpolateST()
20080130 - DWK changed include from mdmsoil51.h to 
           mdmsoiligsm44a.h
20080130 - DWK changed MDMsoil51:: to MDMsoil44::
20080130 - DWK changed ProcessXML43() to ProcessXML44()
20110707 - DWK changed include from mdmsoiligsm44a.h to 
           mdmsoiligsm44c.h
                                   
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

  using std::ceil;
  using std::exp;
  using std::pow;
  
#include<string>
  
  using std::string;


#include "mdmsoiligsm44c.h"

/* *************************************************************
************************************************************* */

MDMsoil44::MDMsoil44() : ProcessXML44()
{


  pctsand = MISSING;
  pctsilt = MISSING;
  pctclay = MISSING;


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

double MDMsoil44::CH4ebullition( const double& ch4con )
{
// times 24.0 to consider the conentration at the static state, 
//   then compared to threshold concentration

  // 500 u M for totatly vegetated site
  double const cmin = 500.0; 

  const double ke = 1.0; // unit h-1

  // percent of vegetated area of the site
  double const punveg = 0.01; 

  // Threshold level of methane when ebullition occurs
  double ch4thresh;

  double ch4ebuflx;  // umol liter-1 hr-1

  double fsmzt;
  

// See Walter et al., 2001

  ch4thresh = cmin * ( 1 + punveg * 0.01);

  if( (ch4con * 24.0) >= ch4thresh ) { fsmzt = 1.0; }
  else { fsmzt = ZERO; }

 
  ch4ebuflx = ke * fsmzt * (ch4con - (ch4thresh/24.0));

//  printf(" %3.2f", ch4ebuflx);

  return ch4ebuflx;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double MDMsoil44::CH4plantTransport( const double& tveg, 
                                     const double& prootd,
                                     const double& pz, 
                                     const double& psoilt,
                                     const double& ch4con )
{

// removing ch4 rate from soils by plants

  const double kp  = 0.01; //  h-1

  double ch4plflxr; // umol liter-1 hr-1

  double froott;
  double fgrowt;
  
  froott = froot( prootd, pz );
  
  fgrowt = fgrow( psoilt );
  
  ch4plflxr = kp * tveg * froott * fgrowt * ch4con;

  return ch4plflxr;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double MDMsoil44::dryCH4diffusion( const double& ddf, 
                                   const double& lowb )

{

// calculate the CH4 concentration profile and fluxes

  const int MAXN = 200; //maximum nodes
 
  // number of nodes, which is related to low boundary (lowb)
  int nodes; 
 
  const double x = 1.0; // depth step 10mm or 1cm
 
  // atmospheric ch4 concentration (uM)
// const double aoc = 0.076; 

  // atmospheric ch4 concentration (u mol cm-3)
  const double aoc = 0.076 * 0.001; 
 
  double a[MAXN];
  double b[MAXN];
  double c[MAXN];
  double d[MAXN];
  double u[MAXN];
 
  double k[MAXN]; // conductivity
  double co[MAXN]; //concentration
  double z[MAXN]; //depth
  double df[MAXN];//diffusion

  int i;
  double ch4flx;

  nodes = (int) (ceil( lowb / 10.0 ));

  z[0] = ZERO;
  for( i = 0; i < (nodes+1); ++i )
  {
    df[i] = ddf;
    z[i+1] = z[i] + x;
  }

// k[0] = k_zero;
  k[0] = ddf / 1.0;
  co[0] = aoc;

  for( i = 1; i < (nodes+1); ++i )
  {
    u[i] = ZERO;
    
    if( i < (nodes+1) ) { k[i] = df[i] / (z[i+1] - z[i]); }
    else { k[i] = ZERO; }
    
    a[i+1] = -k[i];
    b[i] = k[i-1] + k[i];
    c[i] = -k[i];
    d[i] = u[i];
  }

  d[1] = d[1] + k[0] * co[0];
  
  for( i = 1; i < nodes; ++i )
  {
    c[i] = c[i] / b[i];
    d[i] = d[i] / b[i];
    b[i+1] = b[i+1] - a[i+1]*c[i];
    d[i+1] = d[i+1] - a[i+1]*d[i];
  }

  co[nodes]= d[nodes] / b[nodes];

  for( i = nodes; i > 0; --i )
  {
    co[i] = d[i] - c[i] * co[i+1];
    ch4con[i] = co[i] / 0.001;  // converse to uM
  }

  // surface flux of methane, positive indicates flux to soil

  ch4flx = - df[1] * (co[1] - co[0]) / 1.0;  // u mol cm-2 s-1

  ch4flx *= pow( 10.0, 4.0 ); // u mol m-2 s-1;
  ch4flx *= pow( 10.0, -6.0 ); // mol m-2 s-1;
  ch4flx *= 16; // g m-2 s-1;
  ch4flx *= 8.64 * pow( 10.0, 4.0 ); // g m-2 day-1;
  ch4flx *= pow( 10.0, 3.0 ); // mg m-2 day-1;

  return ch4flx;
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double MDMsoil44::dryCH4uptake( const double& ddf, 
                                const double& lowb )
{

  const int MAXN = 200; //maximum nodes
 
  // number of nodes, which is related to low boundary (lowb)
  int nodes; 
  
  const double x = 1.0; // depth step 10mm or 1cm
  
  // atmospheric ch4 concentration (uM)
// const double aoc = 0.076; 

  // atmospheric ch4 concentration (u mol cm-3)
  const double aoc = 0.076 * 0.001; 

  double a[MAXN];
  double b[MAXN];
  double c[MAXN];
  double d[MAXN];
  double u[MAXN];
  
  double k[MAXN]; // conductivity
  double co[MAXN]; //concentration
  double z[MAXN]; //depth
  double df[MAXN];//diffusion

  int i;
  double ch4flx;

  nodes = (int) (ceil( lowb / 10.0 ));

  z[0] = 0.0;
  for( i = 0; i < nodes; ++i )
  {
    df[i] = ddf;
    z[i+1] = z[i] + x;
  }

// k[0] = k_zero;
  k[0] = ddf / 1.0;
  co[0] = aoc;

  for( i = 1; i < (nodes+1); ++i )
  {
    u[i] =  - ch4oxirate1[i] 
            * pow( 10.0, -3.0 ) / (3.6 * pow( 10.0, 3.0 ));  // u mol cm-3 s-1
            
    if( i < (nodes+1) ) { k[i] = df[i] / (z[i+1] - z[i]); }
    else { k[i] = ZERO; }
    
    a[i+1] = -k[i];
    b[i] = k[i-1] + k[i];
    c[i] = -k[i];
    d[i] = u[i];
  }

  d[1] = d[1] + k[0] * co[0];
  
  for( i = 1; i < nodes; ++i )
  {
    c[i] = c[i] / b[i];
    d[i] = d[i] / b[i];
    b[i+1] = b[i+1] - a[i+1]*c[i];
    d[i+1] = d[i+1] - a[i+1]*d[i];
  }

  co[nodes] = d[nodes] / b[nodes];

  for( i = nodes; i > 0; --i )
  {
    co[i] = d[i] - c[i] * co[i+1];
    ch4con[i] = co[i] / 0.001;  // converse to uM
  }

  // surface flux of methane, positive indicates flux to soil

  ch4flx = - df[1] * (co[1] - co[0]) / 1.0;  // u mol cm-2 s-1

  ch4flx *= pow( 10.0, 4.0 ); // u mol m-2 s-1;
  ch4flx *= pow( 10.0, -6.0 ); // mol m-2 s-1;
  ch4flx *= 16; // g m-2 s-1;
  ch4flx *= 8.64 * pow( 10.0, 4.0 ); // g m-2 day-1;
  ch4flx *= pow( 10.0, 3.0 ); // mg m-2 day-1;

  return ch4flx;
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double MDMsoil44::fgrow( const double& soilt )
{

//calculate fgrow
  
  const double amin = ZERO;
  const double a = 4.0;

  double fgrowt;
  double tgr;
  double amax;
  double tmat;

  amax = amin + a;
  if( soilt < 5.0 ) { tgr = 2.0; }
  else { tgr = 7.0; }
 
  tmat = tgr + 10.0;

 if( soilt < tgr ) { fgrowt = amin; }

 if( tmat < soilt ) { fgrowt = amax; }

 if( (soilt <= tmat) & (soilt >= tgr) ) 
 {
   fgrowt = amin 
            + a 
            * (1 - pow( (tmat - soilt) / (tmat - tgr), 2.0 ));
 }
 
 return fgrowt;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double MDMsoil44::froot( const double& rootd, const double& z )
{

//calculate froot

  double froott;
  double rtz;

  // Convert rooting depth from meters to millimeters

  rtz = rootd * 1000.0;

  if( z <= rtz ) { froott = 2 * (1 - z / rtz); }
  else { froott = ZERO; }

  return froott;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MDMsoil44::getecd( ofstream& rflog1 ) 
{
  string ecd;

  cout << "Enter name of file (.ECD) with the MDM soil ";
  cout << "parameter values:" << endl;
  
  cin >> ecd;

  rflog1 << "Enter name of file (.ECD) with the MDM soil ";
  rflog1 << "parameter values:" << endl;

  getecd( ecd );
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MDMsoil44::getecd( const string& ecd )
{

  int comtype;
  int dcmnt;
  ifstream infile;


  infile.open( ecd.c_str(), ios::in );

  if( !infile )
  {
    cerr << endl << "Cannot open " << ecd;
    cerr << " for data input" << endl;
    exit( -1 );
  }

  getXMLrootNode( infile, "mdmsoilECD" );

  for( dcmnt = 1; dcmnt < MAXCMNT; ++dcmnt )
  {
    comtype = getXMLcommunityNode( infile, "mdmsoilECD" );

    if( comtype >= MAXCMNT )
    {
      cerr << endl << "comtype is >= MAXCMNT" << endl;
      cerr << "comtype cannot be greater than ";
      cerr << (MAXCMNT-1) << endl;
      cerr << " in mdmsoilECD" << endl;
      
      exit( -1 );
    }

    PA[comtype]= getXMLcmntArrayDouble( infile,
                                        "mdmsoilECD",
                                         "pa",
                                         comtype );

    endXMLcommunityNode( infile );
  }

  if( dcmnt < MAXCMNT )
  {
    cerr << endl << " Parameters found for only " << dcmnt;
    cerr << " community types out of a maximum of ";
    cerr << (MAXCMNT-1) << " types in mdmsoilECD" << endl;
    
    exit( -1 );
  }

  infile.close();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double MDMsoil44::interpolateLayers( const double layerVar[10], 
                                     const double layerThick[10],
                                     const int& nlayers,
                                     const double& depth )
{

// Interpolate moisture to obtain every 1 cm depth soil 
//   moistures with Lagrange method. Numerical method of 
//   computation and FORTRAN lanuage, Xianyi Zheng, 1986
//   temperatures with Lagrange method

// Numerical method of computation and FORTRAN lanuage, 
//   Xianyi Zheng, 1986

  double f;
  double p;
  int i;
  int j;
  double x0[10];
  double y0[10];

  x0[0] = layerThick[0];
  y0[0] = layerVar[0];
  
  for( i = 1; i < nlayers; ++i )
  {
    x0[i] = x0[i-1] + layerThick[i];
    y0[i] = layerVar[i];
  }


  f = ZERO;
  for( i = 0; i < nlayers; ++i )
  {
    p = 1.0;
    for( j = 0; j < nlayers; ++j )
    {
      if( i != j ) 
      { 
        p = p * (depth - x0[j]) / (x0[i] - x0[j]); 
      }
    }

    f = f + p * y0[i];
  }

  return f;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double MDMsoil44::RedoxP( const double& watertable, 
                          const double& depthz, 
                          const double& wfpsl,
                          const int& pcmnt, 
                          const double& lowb )
{

// calculate the redox potential

// depthz -- depth z
// wfpsl -- the water filled pore space -- pctp

  // daily variation of redox potential (mvday-1) 
  //   Zhang et al., set up as 100  
// const double CR = 100.0; 

  // daily variation of redox potential (mvday-1)
  const double CR = 200.0; 

  // area of the cross section of a typical fine root (cm2) 
  const double FRD = 0.0013;  

  // McClaugherty et al., 1982 estimates fine root length 
  //   per m2 ground area = 1300/130 (mm-2)
  // the above translate to root length density 
  //   (cm root cm-3 soil) by assuming the rooting depth 100cm;
  // RLDL = 0.001 cm root cm-3 soil
  const double RLDL = 0.001;


  double AL;
 
  double ehl;


  
  
  double wfps;  // 0 - 1.0  as percentage

  wfps = wfpsl /100.0;

  AL = FRD * PA[pcmnt] * RLDL;

  if ( (depthz <= lowb) && (depthz > watertable) )
  {
    ehl = CR * (AL - 1.0);
  }
  else
  {
    ehl = CR * (AL + 1.0 - wfps);
  }

  return ehl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double MDMsoil44::setDiffusivityFactor( const double& sand, 
                                        const double& silt, 
                                        const double& clay, 
                                        int satyorn )
{

// See Walter et al., 2001
  
  // diffusion coefficient of CH4 in the air and unsaturated 
  //   soil (0.2 cm-2 s-1)
  const double d1 = 0.2;
  
  // diffusion coefficient of CH4 in the saturated soil and 
  //   water (cm-2 s-1) 
  const double d2 = 0.2 * pow( 10.0, -4.0 ); 
  
  // a tortuosity coefficient, suggesting that the apparent 
  //   path is about two-thirds the length of the real average 
  //   path of diffusion in the soil
  const double dd = 0.66;
  
  // relative volume of coarse pores in sandy soils  
  const double pvsand = 0.45; 
  
  const double pvsilt = 0.20; // silt soils
  const double pvclay = 0.14; // clay soils

  // rate of methane oxidation rate at depth z and time t
  double diffu; 
  
  double ftxt;

  ftxt = (sand * pvsand + silt * pvsilt + clay * pvclay) 
         * 0.01;

  if( 0 == satyorn ) { diffu = d1 * dd * ftxt; }
  else { diffu = d2 * dd * ftxt; }
  
  return diffu;

};

/* *************************************************************
************************************************************* */


/***************************************************************
************************************************************* */

void MDMsoil44::setOxygenConc( const double& afp, 
                               const double& lowb )

{
// calculate the oxygen concentration

  // atmospheric oxygen concentration (280 gm-3)
  const double aoc = 280.0; 

  // soil surface O2 uptake rate
  const double azs = 0.0005; 

  // constant b, Sallam et al., 1984
  const double cb = 0.90; 
 
  // constant m, Sallam et al., 1984
  const double cm = 2.36; 

  // binary diffusion coefficients for O2 (m2s-1) 
  //   (Campbell, 1977)
  const double d0 = 0.0000177; 

  // the atmospheric boundary layer conductance (k(0))
  const double k_zero = 0.01; 

  //maximum nodes
  const int MAXN = 200; 
 
  // number of nodes, which is related to low boundary (lowb)
  int nodes; 
 
  // depth step 10mm or 1cm
  const double x = 0.001; 
 

  double a[MAXN];
  double b[MAXN];
  double c[MAXN];
  double co[MAXN]; //concentration
  double d[MAXN]; 
  double df[MAXN];//diffusion
  int i;
  double k[MAXN]; // conductivity
  double u[MAXN];
  double z[MAXN]; //depth

  nodes = (int) (ceil( lowb / 10.0 ));

  z[0] = ZERO;

  for( i = 0; i < (nodes+1); ++i )
  {
   z[i+1] = z[i] + x;
  }

  k[0] = k_zero;
  co[0] = aoc;

  for( i = 1; i < (nodes+1); ++i )
  {
//   z[i+1] = z[i] + x;
    df[i] = d0 * cb * pow( afp, cm );
    
    u[i] = azs * exp(-z[i] / 0.3) * (z[i+1] - z[i-1])/2.0;
   
    if( i < (nodes+1) ) { k[i] = df[i] / (z[i+1] - z[i]); }
    
    else { k[i] = ZERO; }
   
    a[i+1] = -k[i];
    b[i] = k[i-1] + k[i];
    c[i] = -k[i];
    d[i] = u[i];
  }

  d[1] = d[1] + k[0] * co[0];
  
  for( i = 1; i < nodes; ++i )
  {
    c[i] = c[i] / b[i];
    d[i] = d[i] / b[i];
    b[i+1] = b[i+1] - a[i+1]*c[i];
    d[i+1] = d[i+1] - a[i+1]*d[i];
  }

  co[nodes] = d[nodes] / b[nodes];

  for( i = nodes; i > 0; --i )
  {
    co[i] = d[i] - c[i] * co[i+1];
    
    // covert gm-3 to M;  (32 * pow(10,3) g / pow(10,6)cm = 1 M
//    co[i] = co[i] / (32.0 * pow(10,3));   
//    co[i] = co[i] * pow(10,6); //converse M to micro M (uM)
    oxygenc[i] = co[i];
  }
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MDMsoil44::setUnsaturatedZoneMoisture( const double& ph2otable )
{
  // deterine the water content above water table, 
  //   similar to Granberg et al., 1999
  
  // assume water content linearly decrease with depth
  //   the soil moiture is evaluated for every 10 mm, 
  //   i.e. every 1 cm

  // porosity of soil; = 0.9 cm-3cm-3 (Frolking 1996)
  const double phi = 0.9; 

  // minimum water content at the unstaturated zone (= 0.25)
  const double thetamin = 0.42; 

  // minimum depth with the minimum water content 
  //   ( = 0.1 m = 100mm)
  const double zmin = 100; 


  // the gradient in the linearly decreasing interval
  double az; 
  
  int i;
  
  double q1;
  
  // thickness of saturated soil above 300 mm
  //int satThick;

  double thetas;

  // thickness of unsaturated soil above 300 mm
  int unsatThick;
  
  double z;

  
  az = (phi - thetamin) / zmin;

  thetas = max( thetamin, (phi - az * ph2otable) );
  
//  if( ph2otable <= ZERO )
//  {
//    for( i = 0; i < MXMDMWETNLAY; ++i )
//    {
//      unsatthetaWL[i] = 1.0; 
//    }
//  }
//  else
//  {
    unsatThick = (int) (min( ceil( ph2otable/10.0 ), 
                             MXMDMWETNLAY ));
       
    for( i = 0; i < unsatThick; ++i )
    {
      z = 10 * i;
      q1 = pow( (z / ph2otable), 2.0 );
    
      unsatthetaWL[i] = min( phi, 
                             (thetas + (phi - thetas) * q1) );
    }
//  }
	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double MDMsoil44::wetCH4diffusion( const double& lowb, 
                                   const double& watertable )
{

// for anoxic and oxic soils

  const int MAXN = 200; //maximum nodes
 
 // number of nodes, which is related to low boundary (lowb)
  int nodes; 
 
  const double x = 1.0; // depth step 10mm or 1cm
 
 // atmospheric ch4 concentration (uM)
// const double aoc = 0.076; 

  // atmospheric ch4 concentration (u mol cm-3)
  const double aoc = 0.076 * 0.001; 

  double a[MAXN];
  double b[MAXN];
  double c[MAXN];
  double d[MAXN]; 
  double u[MAXN];
  
  double k[MAXN]; // conductivity
  double co[MAXN]; //concentration
  double z[MAXN]; //depth
  double df[MAXN];//diffusion

  int i;
  double ch4flx;

  nodes = (int) (ceil( lowb / 10.0 ));

  z[0] = ZERO;
  
  for( i = 0; i < nodes; ++i )
  {
    df[i] = fdfs[i];
    z[i+1] = z[i] + x;
  }

// k[0] = k_zero;
  k[0] = df[0] / 1.0;
  co[0] = aoc;

  for( i = 1; i < (nodes+1); ++i )
  {
    if( i < (watertable /10.0) ) { u[i] = ZERO; }
    else 
    {
      u[i] = ch4ratesat[i] 
             * pow( 10.0, -3.0 ) / (3.6 * pow( 10.0, 3.0 )); // u mol cm-3 s-1
    }
    
    if( i < (nodes+1) ) { k[i] = df[i] / (z[i+1] - z[i]); }
    else { k[i] = ZERO; }
    
    a[i+1] = -k[i];
    b[i] = k[i-1] + k[i];
    c[i] = -k[i];
    d[i] = u[i];
  }

  d[1] = d[1] + k[0] * co[0];
  
  for( i = 1; i < nodes; ++i )
  {
    c[i] = c[i] / b[i];
    d[i] = d[i] / b[i];
    b[i+1] = b[i+1] - a[i+1]*c[i];
    d[i+1] = d[i+1] - a[i+1]*d[i];
  }

  co[nodes]= d[nodes] / b[nodes];

  for( i = nodes; i > 0; --i )
  {
    co[i] = d[i] - c[i] * co[i+1];
    ch4con_b[i] = co[i] / 0.001;  // converse to uM
  }

  // surface flux of methane, positive indicates flux to soil

  ch4flx = - df[1] * (co[1] - co[0]) / 1.0;  // u mol cm-2 s-1

  ch4flx *= pow( 10.0, 4.0 ); // u mol m-2 s-1;
  ch4flx *= pow( 10.0, -6.0 ); // mol m-2 s-1;
  ch4flx *= 16; // g m-2 s-1;
  ch4flx *= 8.64 * pow( 10.0, 4.0 ); // g m-2 day-1;
  ch4flx *= pow( 10.0, 3.0 ); // mg m-2 day-1;

  return ch4flx;
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double MDMsoil44::wetCH4flux( const double& lowb, 
                              const double& watertable, 
                              const double& soilt )
{

// for anoxic and oxic soils, eventual flux and concentration 
//   profile

  const int MAXN = 200; //maximum nodes
 
  // number of nodes, which is related to low boundary (lowb)
  int nodes; 
 
  // depth step 10mm or 1cm
  const double x = 1.0; 
 
  // atmospheric ch4 concentration (uM)
// const double aoc = 0.076; 
  
  // atmospheric ch4 concentration (u mol cm-3)
  const double aoc = 0.076 * 0.001; 

  double a[MAXN];
  double b[MAXN];
  double c[MAXN];
  double co[MAXN]; //concentration
  double d[MAXN];
  double df[MAXN];//diffusion
  double k[MAXN]; // conductivity
  double u[MAXN];
  double z[MAXN]; //depth

  int i;
  double ch4flx;

  nodes = (int) (ceil( lowb / 10.0 ));

  z[0] = ZERO;
  for( i = 0; i < nodes; ++i )
  {
//   if (soilt < 0.0) { ch4ratesat[i] = 0.0001;}
    df[i] = fdfs[i];
//   if ( (df[i] < pow( 10.0, -3.0 )) || (df[i] > 1.0)) df[i] = 0.01;

//  printf(" %5.3f %5.3f", ch4ratesat[i], ch4oxirate1[i]);
//  if (ch4ratesat[i] > 5.0) ch4ratesat[i] = .010;
//   ch4oxirate1[i] =0.01;

  //printf(" %5.3f", ch4ratesat[i]);

    z[i+1] = z[i] + x;
  }


// k[0] = k_zero;
  k[0] = df[0] / 1.0;
  co[0] = aoc;

  for( i = 1; i < (nodes+1); ++i )
  {
    if( i < (watertable /10.0) )
    {
      u[i]=  - ch4oxirate1[i] 
             * pow( 10.0, -3.0 ) 
             / (3.6 * pow( 10.0, 3.0 ));  // u mol cm-3 s-1;
    }
    else 
    {
      u[i] =  ch4ratesat[i] 
              * pow( 10.0, -3.0 ) 
              / (3.6 * pow( 10.0, 3.0 ));  // u mol cm-3 s-1
    }
    
    if( i < (nodes+1) ) { k[i] = df[i] / (z[i+1] - z[i]); }
    else { k[i] = ZERO; }
    
    a[i+1] = -k[i];
    b[i] = k[i-1] + k[i];
    c[i] = -k[i];
    d[i] = u[i];
  }

  d[1] = d[1] + k[0] * co[0];
  
  for( i = 1; i < nodes; ++i )
  {
    c[i] = c[i] / b[i];
    d[i] = d[i] / b[i];
    b[i+1] = b[i+1] - a[i+1]*c[i];
    d[i+1] = d[i+1] - a[i+1]*d[i];
  }

  co[nodes]= d[nodes] / b[nodes];

  for( i = nodes; i > 0; --i )
  {
    co[i] = d[i] - c[i] * co[i+1];
    ch4con[i] = co[i] / 0.001;  // converse to uM
  }

  // surface flux of methane, positive indicates flux to soil

  ch4flx = - df[1] * (co[1] - co[0]) / 1.0;  // u mol cm-2 s-1

  ch4flx *= pow( 10.0, 4.0 ); // u mol m-2 s-1;
  ch4flx *= pow( 10.0, -6.0 ); // mol m-2 s-1;
  ch4flx *= 16; // g m-2 s-1;
  ch4flx *= 8.64 * pow( 10.0, 4.0 ); // g m-2 day-1;
  ch4flx *= pow( 10.0, 3.0 ); // mg m-2 day-1; 
 
  if( soilt < ZERO ) { ch4flx = 0.001;}
  
  return ch4flx;
  
};





