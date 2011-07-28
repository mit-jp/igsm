/* *************************************************************
****************************************************************
MITTEM44B.CPP - Terrestrial Ecosystem Model Version 4.4
****************************************************************

Modifications:

20060128 - DWK created by modifying ttem437.cpp
20060128 - DWK changed include from ttem437.h to mittem437.h
20060128 - DWK changed TTEM43:: to MITTEM43::
20060128 - DWK added I_CH4EMS, I_CH4CSMP, I_CH4FLX, I_CO2NFLX,  
           I_CO2DNFLX, I_NOFLX, I_N2OFLX, I_N2ONFLX, I_N2ODNFLX
           and I_N2FLX to MITTEM43()
20060129 - DWK changed functions calls to pstate[I_SM} to
           soil.getMOIST() in cropDynamics() and 
           natvegDynamics()            
20060129 - DWK changed functions calls to pstate[I_VSM} to
           soil.getVSM() in cropDynamics() and natvegDynamics()            
20060129 - DWK deleted water variables from boundcon(), delta(),
           ECDsetODEstate(), massbal()
20060129 - DWK added const int equil to initrun()
20060129 - DWK added NEM to stepmonth()
20060202 - DWK added I_NIRR, I_PAR, I_TAIR, I_PREC, I_RAIN,   
           I_SNWFAL, I_CO2, I_AOT40 to MITTEM43(), 
20070129 - DWK changed include from mittem438.h to mittem44.h           
20070129 - DWK changed MITTEM43:: to MITTEM44::
20080130 - DWK changed include from mittem44.h to mittem44a.h
20080826 - DWK changed include from mittem44a.h to mittem44b.h
                                                                                                                                                
****************************************************************
************************************************************** */

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

#ifdef BORLAND_CPP
  #include<stdlib>
#else 
  #include<cstdlib>
#endif

  using std::exit;
  using std::atof;
  using std::atoi;

#include<cmath>

  using std::exp;
  using std::fabs;
  using std::modf;
  using std::pow;

#include<vector>

  using std::vector;
    
#include<string>
  
  using std::string;

// #define CALIBRATE_TEM

#ifdef CALIBRATE_TEM
  #include<conio.h>
  #include<cctype>
    using std::toupper;
#endif

#include "mittem44b.h"

/* **************************************************************
************************************************************** */

// Initialization of static members

int MITTEM44::avlnflag = 0;
int MITTEM44::nfeed = 0;
int MITTEM44::rheqflag = 0;
int MITTEM44::moistlim = 0;
int MITTEM44::o3flag = 0;
int MITTEM44::initbase = 0;
int MITTEM44::baseline = 0;
int MITTEM44::intflag = 0;

int MITTEM44::maxnrun = 0;
int MITTEM44::equil = 0;
int MITTEM44::runsize = 0;
int MITTEM44::maxyears = 0;
int MITTEM44::strteq = 0;
int MITTEM44::endeq = 0;
int MITTEM44::startyr = 0;
int MITTEM44::endyr = 0;
int MITTEM44::diffyr = 0;
int MITTEM44::wrtyr = 0;

double MITTEM44::ctol = 1.0;
double MITTEM44::ntol = 0.02;
double MITTEM44::wtol = 0.01;

// Initialization of adaptive integrator variables

double MITTEM44::inittol = 0.01;
int MITTEM44::maxit = 20;
long MITTEM44::maxitmon = 2000;

double MITTEM44::a1 =   0.115740741;

double   MITTEM44::a3 =   0.548927875;
double  MITTEM44::a31 =   0.09375;
double  MITTEM44::a32 =   0.28125;

double   MITTEM44::a4 =   0.535331384;
double  MITTEM44::a41 =   0.879380974;
double  MITTEM44::a42 =  -3.277196177;
double  MITTEM44::a43 =   3.320892126;

double   MITTEM44::a5 =  -0.20;
double  MITTEM44::a51 =   2.032407407;
double  MITTEM44::a52 =  -8.0;
double  MITTEM44::a53 =   7.173489279;
double  MITTEM44::a54 =  -0.2058966866;

double   MITTEM44::b1 =   0.118518519;
double   MITTEM44::b3 =   0.518986355;
double   MITTEM44::b4 =   0.50613149;
double   MITTEM44::b5 =  -0.18;
double   MITTEM44::b6 =   0.036363636;
double  MITTEM44::b61 =  -0.296296296;
double  MITTEM44::b62 =   2.0;
double  MITTEM44::b63 =  -1.381676413;
double  MITTEM44::b64 =   0.45297271;
double  MITTEM44::b65 =  -0.275;


/* **************************************************************
************************************************************** */

MITTEM44::MITTEM44() : predstr( NUMTEM )
{
  tol = inittol;
  syint = 1;  
  totyr = -99;
  rheqflag = 0;


#ifdef CALIBRATE_TEM
  soilfile = "tsoil43d.ecd";
  rootfile = "trotz43d.ecd";
  vegfile = "tveg43d.ecd";
  leaffile = "tleaf43d.ecd";
  mcrvfile = "tmcrv43d.ecd";
  agfile = "ag437.ecd";

  rheqflag = 0;

  sey[0] = GET_LAI;
  sey[1] = GET_NPP;
  sey[2] = GET_VNUP;

  swy[0] = GET_RAIN;
  swy[1] = GET_SNWINF;
  swy[2] = GET_VSM;
  swy[3] = GET_PET;
  swy[4] = GET_EET;
#endif

 
// Identify potential output variables from TEM

// Ecosystem carbon pools determined by the integrator**********

  // vegetation carbon
  predstr.at( I_VEGC ) = "VEGC";  
  
  // reactive soil organic carbon      
  predstr.at( I_SOLC ) = "SOILORGC";   
  
// Ecosystem nitrogen pools determined by the integrator********

  // vegetation structural nitrogen
  predstr.at( I_STRN ) = "VSTRUCTN";

  // vegetation labile nitrogen
  predstr.at( I_STON ) = "VSTOREN";

  // reactive soil organic nitrogen
  predstr.at( I_SOLN ) = "SOILORGN";  

  // soil available nitrogen
  predstr.at( I_AVLN ) = "AVAILN";

/// Phenology variables determined by the integrator************

  // un-normalized relative phenology variable
  predstr.at( I_UNRMLF ) = "UNRMLEAF";

  // normalized relative phenology variable (0 - 1.0)
  predstr.at( I_LEAF ) = "LEAF";

  // leaf area index
  predstr.at( I_LAI ) = "LAI";    
  
  // foliar projected cover
  predstr.at( I_FPC ) = "FPC";    


// Carbon fluxes for ecosystems ********************************

  // GPP not limited by nutrient availability
  predstr.at( I_INGPP ) = "VEGINGPP";

  // gross primary production
  predstr.at( I_GPP ) = "GPP";
  
  // Direct ozone effects 
  predstr.at( I_FOZONE ) = "FOZONE";
  
  // Indirect ozone effects 
  predstr.at( I_FINDOZONE ) = "FINDOZON";     

 // NPP not limited by nutrient availability
  predstr.at( I_INNPP ) = "VEGINNPP";

  // net primary production
  predstr.at( I_NPP ) = "NPP";      

  // gross plant respiration
  predstr.at( I_GPR ) = "GPR";      

  // vegetation maintenance respiration
  predstr.at( I_RVMNT ) = "RVMAINT";

  // vegetation growth respiration
  predstr.at( I_RVGRW ) = "RVGRWTH";

  // litterfall carbon
  predstr.at( I_LTRC ) = "LTRC";   
  
  // heterotrophic respiration
  predstr.at( I_RH ) = "RH";       

// Nitrogen fluxes for ecosystems determined by the integrator

  // total nitrogen inputs into ecosystem
  predstr.at( I_NINP ) = "NINPUT";

  // nitrogen fertilization
  predstr.at( I_AGFRTN ) = "AGFERTN";

  // VEGNUP not limited by carbon availability
  predstr.at( I_INNUP ) = "VEGINNUP";

  // nitrogen uptake by vegetation
  predstr.at( I_VNUP ) = "VEGNUP";

  // vegetation nitrogen uptake for structural components
  predstr.at( I_VSUP ) = "VEGSUP";

  // vegetation nitrogen uptake for labile components
  predstr.at( I_VLUP ) = "VEGLUP";

  // nitrogen mobilization by vegetation
  predstr.at( I_VNMBL ) = "VNMOBIL";

  // nitrogen resorption by vegetation
  predstr.at( I_VNRSRB ) = "VNRESORB";

  // litterfall nitrogen from vegetation
  predstr.at( I_LTRN ) = "LTRN";

  // total nitrogen immobilization
  predstr.at( I_MNUP ) = "MICRONUP";

  // net nitrogen mineralization
  predstr.at( I_NMIN ) = "NETNMIN";

  // Total nitrogen losses from ecosystems
  predstr.at( I_NLST ) = "NLOST";

// Water fluxes determined by the integrator********************

  // Irrigation
  predstr.at( I_AGIRRIG ) = "IRRIGATE"; 

  // Initial estimated evapotranspiration
  predstr.at( I_INEET ) = "INEET";

  // estimated evapotranspiration
  predstr.at( I_EET ) = "EET";


// Other ecosystem carbon pools ********************************

  // total carbon pool found in ecosystem excluding products
  predstr.at( I_TOTEC ) = "TOTEC";

  // total carbon (including products if present)
  predstr.at( I_TOTC ) = "TOTALC";     

// Other ecosystem nitrogen pools ******************************

  // total nitrogen stored in vegetation
  predstr.at( I_VEGN ) = "VEGN";

// Other ecosystem water pools ******************************

  predstr.at( I_SNWPCK ) = "SNOWPACK";  // snowpack


// Ecosystem water pools determined by the integrator***********

  // available soil moisture
  predstr.at( I_AVLW ) = "AVAILW";

  // soil moisture
  predstr.at( I_SM ) = "SMOIST";       

  // volumetric soil moisture
  predstr.at( I_VSM ) = "VSM";      

  // soil moisture expressed as percent total porosity
  predstr.at( I_PCTP ) = "PCTP";


// Other carbon fluxes for ecosystems **************************

  // net ecosystem production
  predstr.at( I_NEP ) = "NEP";     

  // net carbon exchange of ecosystem with atmosphere
  predstr.at( I_NCE ) = "NCE";

  // Methane emissions from soils
  predstr.at( I_CH4EMS ) = "CH4EMISS";
  
  // Methane consumption by soils
  predstr.at( I_CH4CSMP )= "CH4CSMP";    
  
  // Methane flux (CH4EMISS - CH4CSMP)
  predstr.at( I_CH4FLX ) = "CH4FLUX";     

  // Carbon dioxide flux from NEM nitrification subroutine
  predstr.at( I_CO2NFLX ) = "CO2NFLUX";     

  // Carbon dioxide flux from NEM denitrification subroutine
  predstr.at( I_CO2DNFLX ) = "CO2DNFLX";     

// Other nitrogen fluxes for ecosystems ************************

  // Net nitric oxide flux from soil
  predstr.at( I_NOFLX ) = "NOFLUX";
  
  // Net nitrous oxide flux from soil 
  predstr.at( I_N2OFLX ) = "N2OFLUX";

  // Nitrous oxide flux resulting from nitrification 
  predstr.at( I_N2ONFLX ) = "N2ONFLUX";     

  // Nitrous oxide flux resulting from denitrification
  predstr.at( I_N2ODNFLX ) = "N2ODNFLX";     

  // Net dinitrogen flux from soil
  predstr.at( I_N2FLX ) = "N2FLUX";


// Other water fluxes ******************************************

  // potential evapotranspiration
  predstr.at( I_PET ) = "PET";


// Carbon stocks in products ***********************************

  // carbon in agricultural products
  predstr.at( I_AGPRDC ) = "AGPRODC";

  // carbon pool of products that decompose in 10 years
  predstr.at( I_PROD10C ) = "PROD10C";

  // carbon pool of products that decompose in 100 years
  predstr.at( I_PROD100C ) = "PROD100C";

  // carbon in all product pools
  predstr.at( I_TOTPRDC ) = "TOTPRODC";


// Carbon stocks in crop residue and stubble********************

  // carbon in crop residue
  predstr.at( I_RESIDC ) = "RESIDC";

  // stubble carbon
  predstr.at( I_AGSTUBC ) = "CRPSTUBC";   


// Nitrogen stocks in products *********************************

  // nitrogen in agricultural products
  predstr.at( I_AGPRDN ) = "AGPRODN";

  // nitrogen pool of products that decompose in 10 years
  predstr.at( I_PROD10N ) = "PROD10N";

  // nitrogen pool of products that decompose in 100 years
  predstr.at( I_PROD100N ) = "PROD100N";

  // nitrogen in all product pools
  predstr.at( I_TOTPRDN ) = "TOTPRODN";


// Nitrogen stocks in crop residue and stubble******************

  // nitrogen in crop residue
  predstr.at( I_RESIDN ) = "RESIDN";

  // stubble nitrogen
  predstr.at( I_AGSTUBN ) = "CRPSTUBN";


// Carbon fluxes associated with agricultural conversion *******

  // carbon loss from the ecosystem during conversion
  predstr.at( I_CNVRTC ) = "CONVERTC";

  // carbon loss from vegetation during conversion
  predstr.at( I_VCNVRTC ) = "VCONVRTC";

  // carbon loss from soils during conversion
  predstr.at( I_SCNVRTC ) = "SCONVRTC";

  // carbon associated with slash left after conversion
  predstr.at( I_SLASHC ) = "SLASHC";

  // carbon flux from ecosystem (NEP+CONVERTC)
  predstr.at( I_CFLX ) = "CFLUX";


// Nitrogen fluxes associated with agricultural conversion *****

  // nitrogen loss from the ecosystem during conversion
  predstr.at( I_CNVRTN ) = "CONVERTN";

  // nitrogen loss from vegetation during conversion
  predstr.at( I_VCNVRTN ) = "VCONVRTN";

  // nitrogen loss from soils during conversion
  predstr.at( I_SCNVRTN ) = "SCONVRTN";

  // nitrogen associated with slash left after conversion
  predstr.at( I_SLASHN ) = "SLASHN";

  // Total organic N mineralized and retained in ecosystem
  //   after disturbance
  predstr.at( I_NRETNT ) = "NRETENT";

  // Vegetation N mineralized and retained in ecosystem
  //   after disturbance
  predstr.at( I_NVRTNT ) = "NVRETENT";

  // Soil organic N mineralized and retained in ecosystem
  //   after disturbance
  predstr.at( I_NSRTNT ) = "NSRETENT";


// Carbon and nitrogen fluxes to/from products *****************

  // carbon loss to formation of agricultural products
  predstr.at( I_AGFPRDC ) = "AGFPRODC";

  // nitrogen loss to formation of agricultural products
  predstr.at( I_AGPRDN ) = "AGFPRODN";

  // carbon loss to crop residue
  predstr.at( I_FRESIDC ) = "FRESIDC";

  // nitrogen loss to crop residue
  predstr.at( I_FRESIDN ) = "FRESIDN";

  // carbon loss to resulting from decomposition of agricultural
  //   products
  predstr.at( I_AGPRDFC ) = "AGPRODFC";

  // nitrogen loss resulting from decomposition of agricultural
  //   products
  predstr.at( I_AGPRDFN ) = "AGPRODFN";

  // carbon loss from crop residue
  predstr.at( I_RESIDFC ) = "RESIDFC";

  // nitrogen loss from crop residue
  predstr.at( I_RESIDFN ) = "RESIDFN";

  // carbon loss to formation of products that decompose in
  //  10 years
  predstr.at( I_PRDF10C ) = "PRDF10C";

  // nitrogen loss to formation of products that decompose in
  //   10 years
  predstr.at( I_PRDF10N ) = "PRDF10N";

  // carbon loss resulting from decomposition of PROD10C
  predstr.at( I_PRD10FC ) = "PRD10FC";

  // nitrogen loss resulting from decomposition of PROD10N
  predstr.at( I_PRD10FN ) = "PRD10FN";

  // carbon loss to formation of products that decompose in
  //  100 years
  predstr.at( I_PRDF100C ) = "PRDF100C";

  // nitrogen loss to formation of products that decompose in
  //   100 years
  predstr.at( I_PRDF100N ) = "PRDF100N";

  // carbon loss resulting from decomposition of PROD100C
  predstr.at( I_PRD100FC ) = "PRD100FC";

  // nitrogen loss resulting from decomposition of PROD100N
  predstr.at( I_PRD100FN ) = "PRD100FN";

  // carbon loss to the formation of all products
  predstr.at( I_TOTFPRDC ) = "TOTFPRDC";

  // nitrogen loss to the formation of all products
  predstr.at( I_TOTFPRDN ) = "TOTFPRDN";

  // carbon loss resulting from decomposition of all products
  predstr.at( I_TOTPRDFC ) = "TOTPRDFC";

  // nitrogen loss resulting from decomposition of all products
  predstr.at( I_TOTPRDFN ) = "TOTPRDFN";


// Agro-Ecosystem carbon and nitrogen pools *********************

  // crop carbon
  predstr.at( I_CROPC ) = "CROPC";       

  // carbon in natural vegetation
  predstr.at( I_NATVEGC ) = "NATVEGC";

  // crop nitrogen
  predstr.at( I_CROPN ) = "CROPN";       

  // nitrogen in natural vegetation
  predstr.at( I_NATVEGN ) = "NATVEGN";

  // crop structural N
  predstr.at( I_CSTRN ) = "CROPSTRN";    

  // structural N in natural vegetation
  predstr.at( I_NATSTRN ) = "NATSTRN";

  // crop labile N
  predstr.at( I_CSTON ) = "CROPSTON";    

  // labile N stored in natural vegetation
  predstr.at( I_NATSTON ) = "NATSTON";


// Crop phenology **********************************************

  // unnormalized leaf in crops
  predstr.at( I_CROPULF ) = "CRPUNMLF";

  // unnormalized leaf in natural vegetation
  predstr.at( I_NATULF ) = "NATUNMLF";

  // leaf of crops
  predstr.at ( I_CROPLEAF ) = "CROPLEAF";

  // leaf of natural vegetation
  predstr.at( I_NATLEAF ) = "NATLEAF";

  // leaf area index (LAI) of crops
  predstr.at( I_CROPLAI ) = "CROPLAI";

  // leaf area index (LAI) of natural vegetation
  predstr.at( I_NATLAI ) = "NATLAI";

  // foliar projected cover (FPC) of crops
  predstr.at( I_CROPFPC ) = "CROPFPC";

  // foliar projected cover (FPC) of natural vegetation
  predstr.at( I_NATFPC ) = "NATFPC";


// Additional carbon fluxes for agro-ecosystems *****************

  // GPP of crops not limited by nutrient availability
  predstr.at( I_AGINGPP ) = "CRPINGPP";

  // GPP of natural vegetation not limited by 
  //   nutrient availability
  predstr.at( I_NATINGPP ) = "NATINGPP";

  // gross primary production (GPP) of crops
  predstr.at( I_AGGPP ) = "CROPGPP";

  // gross primary production of natural vegetation
  predstr.at( I_NATGPP ) = "NATGPP";

  // NPP of crops not limited by nutrient availability
  predstr.at( I_AGINNPP ) = "CRPINNPP";

  // NPP of natural vegetation not limited by 
  //   nutrient availability
  predstr.at( I_NATINNPP ) = "NATINNPP";

  // net primary production (NPP) of crops
  predstr.at( I_AGNPP ) = "CROPNPP";

  // net primary production (NPP) of natural vegetation
  predstr.at( I_NATNPP ) = "NATNPP";

  // gross plant respiration of crops
  predstr.at( I_AGGPR ) = "CROPGPR";

  // gross plant respiration of natural vegetation
  predstr.at( I_NATGPR ) = "NATGPR";

  // maintenance respiration of crop plants
  predstr.at( I_AGRVMNT ) = "CRPRMNT";

  // maintenance respiration of natural vegetation
  predstr.at( I_NATRVMNT ) = "NATRVMNT";

  // growth respiration of crop plants
  predstr.at( I_AGRVGRW ) = "CRPRGRW";

  // growth respiration of natural vegetation
  predstr.at( I_NATRVGRW ) = "NATRGRW";

  // litterfall carbon from crops
  predstr.at( I_AGLTRC ) = "CROPLTRC";

  // litterfall carbon from natural vegetation
  predstr.at( I_NATLTRC ) = "NATLTRC";

  // Additional nitrogen fluxes for agro-ecosystems ************

  // nitrogen uptake by crops not limited by carbon availability
  predstr.at( I_AGINNUP ) = "CRPINNUP";

  // nitrogen uptake by natural vegetation not limited by carbon
  //   availability
  predstr.at( I_NATINNUP ) = "NATINNUP";

  // nitrogen uptake by crops
  predstr.at( I_AGVNUP ) = "CROPNUP";

  // nitrogen uptake by natural vegetation
  predstr.at( I_NATVNUP ) = "NATVNUP";

  // nitrogen uptake for structural components of crops
  predstr.at( I_AGVSUP ) = "CROPSUP";

  // nitrogen uptake for structural components of natural
  //  vegetation
  predstr.at( I_NATVSUP ) = "NATVSUP";

  // nitrogen uptake for labile components of crops
  predstr.at( I_AGVLUP ) = "CROPLUP";

  // nitrogen uptake for labile components of natural vegetation
  predstr.at( I_NATVLUP ) = "NATVLUP";

  // nitrogen mobilization by crops
  predstr.at( I_AGVNMBL ) = "CRPNMOBL";

  // nitrogen mobilization by natural vegetation
  predstr.at( I_NATVNMBL ) = "NATVNMBL";

  // nitrogen resorption by crops
  predstr.at( I_AGVNRSRB ) = "CRPNRSRB";

  // nitrogen resorption by natural vegetation
  predstr.at( I_NVNRSRB ) = "NVNRSRB";

  // litterfall nitrogen from crops
  predstr.at( I_AGLTRN ) = "CROPLTRN";

  //litterfall nitrogen from natural vegetation
  predstr.at( I_NATLTRN ) = "NATLTRN";


  // Climate variables 
  predstr.at(I_NIRR ) = "NIRR";
  
  predstr.at( I_PAR ) = "PAR";
  
  predstr.at( I_TAIR ) = "TAIR"; 
  
  predstr.at( I_PREC ) = "PREC"; 

  predstr.at( I_SRFRUN ) = "SURFRUN"; 
  
  predstr.at( I_DRAIN ) = "DRAINAGE";
  
  predstr.at( I_CO2 ) = "ATMSCO2"; 
  
  predstr.at( I_AOT40 ) = "AOT40"; 

  dbugflg = 0;

};

/* **************************************************************
************************* Functions *****************************
************************************************************** */


/* *************************************************************
************************************************************* */

int MITTEM44::adapt( const int& numeq, 
                     double pstate[], 
                     const double& ptol,
                     const int& pdm )
{

  int i;
  double ipart;
  double fpart;
  double time = ZERO;
  double dt = 1.0;
  int mflag = 0;
  long nintmon = 0;
  double oldstate[NUMEQ];

  blackhol = 0;

  while( time != 1.0 )
  {
    test = REJECT;
    if( 1 == syint )
    {
      while( test != ACCEPT )
      {
	      rkf( numeq, pstate, dt, pdm );
	
	      test = boundcon( dum4, error, ptol );

        #ifdef CALIBRATE_TEM
	        if( test != ACCEPT )
          {
            // Display ODE errors to DOS screen

            pcdisplayODEerr( test, pstate );
          }
        #endif
	
	      if( dt <= pow( 0.5, maxit ) )
        {
	        test = ACCEPT;

	        mflag = 1;
        
          if( 0 == nintmon )
          {
            for( i = 0; i < numeq; ++i ) 
            { 
              oldstate[i] = pstate[i]; 
            }
          }
	
	  ++nintmon;
	}

        if( ACCEPT == test )
        {
          for( i = 0; i < numeq; ++i ) { pstate[i] = dum4[i]; }
        
          time += dt;

          #ifdef CALIBRATE_TEM
            // Display time updates to the DOS screen

            pcdisplayDT( time, dt );
          #endif

          fpart = modf( (0.01 + (time/(2.0*dt))),&ipart );
        
          if( fpart < 0.1 && dt < 1.0) { dt *= 2.0; }
        }
        else { dt *= 0.500; }

        if( nintmon == maxitmon )
        {
          time = 1.0;
          blackhol = 1;
        
          for( i = 0; i < numeq; ++i ) { pstate[i] = oldstate[i]; }
        }
      }
    }    /* end rkf integrator */
  }      /* end time while */

  return mflag;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITTEM44::askODE( ofstream& rflog1 )
{


/* **************************************************************
	      Parameters for Adaptive Integrator
************************************************************** */

  cout << endl << "Enter the proportional tolerance for the integrator: ";

  cin >> inittol;

  rflog1 << endl;
  rflog1 << "Enter the proportional tolerance for the integrator: ";
  rflog1 << inittol << endl;

  cout << "Enter the maximum number of iterations in the integrator: ";

  cin >> maxit;

  rflog1 << endl;
  rflog1 << "Enter the maximum number of iterations in the integrator: ";
  rflog1 << maxit << endl;

  cout << "Enter the maximum number of times in a month that the" << endl;
  cout << "integrator can reach the maximum number of iterations: ";

  cin >> maxitmon;

  rflog1 << endl;
  rflog1 << "Enter the maximum number of times in a month that the" << endl;
  rflog1 << "integrator can reach the maximum number of iterations: ";
  rflog1 << maxitmon << endl;

};

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */

int MITTEM44::boundcon( double ptstate[],
                        double err[],
                        const double& ptol )
{

  int test = ACCEPT;

// Check carbon and nitrogen state variables

  if( err[I_VEGC] > fabs( ptol * ptstate[I_VEGC] ) )
  {
    return test = temkey( I_VEGC )+1;
  }

  if( err[I_SOLC] > fabs( ptol * ptstate[I_SOLC] ) )
  {
    return test = temkey( I_SOLC )+1;
  }

  if( 1 == nfeed
      && err[I_STRN] > fabs( ptol * ptstate[I_STRN] ) )
  {
    return test = temkey( I_STRN )+1;
  }

  if( 1 == nfeed
      && err[I_STON] > fabs( ptol * ptstate[I_STON] ) )
  {
    return test = temkey( I_STON )+1;
  }

  if( 1 == nfeed
      && err[I_SOLN] > fabs( ptol * ptstate[I_SOLN] ) )
  {
    return test = temkey( I_SOLN )+1;
  }

  if( 1 == nfeed
      && err[I_AVLN] > fabs( ptol * ptstate[I_AVLN] ) )
  {
    return test = temkey( I_AVLN )+1;
  }

  if( err[I_GPP] > fabs( ptol * ptstate[I_GPP] ) )
  {
    return test = temkey( I_GPP )+1;
  }

  if( err[I_NPP] > fabs( ptol * ptstate[I_NPP] ) )
  {
    return test = temkey( I_NPP )+1;
  }

  if( err[I_RVMNT] > fabs( ptol * ptstate[I_RVMNT] ) )
  {
    return test = temkey( I_RVMNT )+1;
  }

  if( 1 == nfeed
      && err[I_VNUP] > fabs( ptol * ptstate[I_VNUP] ) )
  {
    return test = temkey( I_VNUP )+1;
  }

  if( 1 == nfeed
      && err[I_VSUP] > fabs( ptol * ptstate[I_VSUP] ) )
  {
    return test = temkey( I_VSUP )+1;
  }

  if( 1 == nfeed
      && err[I_VNMBL] > fabs( ptol * ptstate[I_VNMBL] ) )
  {
    return test = temkey( I_VNMBL )+1;
  }

  return test;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITTEM44::cropDynamics( double pstate[] )
{
  int perennial = 0;
                           
  // Assume no agricultural N fertilization (ag.fertn) and no
  //  nutrients resulting from agricultural conversion

  ag.fertn = ZERO;
  soil.setNINPUT( ag.getNRETENT() ); 

  soil.setEET( soil.getINEET() );

  soil.setVSM( (soil.getMOIST() / (soil.getROOTZ() * 1000.0)) );
  
  soil.setKH2O( soil.getVSM(), moistlim );


  // Note: Microbes are assumed to be acting on "old" carbon
  //   (i.e. natural vegetation - veg.cmnt) rather than 
  //   "new" carbon associated with crops (i.e. ag.cmnt)

  microbe.updateDynamics( veg.cmnt,
                          soil.getPCTFLDCAP(),
                          pstate[I_SOLC],
                          pstate[I_SOLN],
                          soil.getMOIST(),
                          soil.getVSM(),
                          pstate[I_AVLN],
                          moistlim,
                          ag.tillflag,
                          ag.getTILLFACTOR( veg.cmnt ),
                          soil.getKH2O() );

  if( ag.getGROWDD() >= GDDSEED )
  {
    if( 0 == moistlim )
    {
      if( ag.getCROPPRVPETMX() < atms.getPET() )
      {
        atms.setPRVPETMX( atms.getPET() );
      }
      else
      {
        atms.setPRVPETMX( ag.getCROPPRVPETMX() );
      }
    }
    else
    {
      if( ag.getCROPPRVEETMX() < soil.getEET() )
      {
        soil.setPRVEETMX( soil.getEET() );
      }
      else
      {
        soil.setPRVEETMX( ag.getCROPPRVEETMX() );
      }
    }    

    veg.updateDynamics( ag.cmnt,
                        atms.getCO2(),
                        atms.getAOT40(),
                        soil.getNINPUT(),
                        atms.getPAR(),
                        atms.getPET(),
                        atms.getPRVPETMX(),
                        soil.getEET(),
                        soil.getPRVEETMX(),
                        pstate[I_VEGC],
                        pstate[I_STRN],
                        pstate[I_STON],
                        soil.getMOIST(),
                        pstate[I_AVLN],
                        moistlim,
                        nfeed,
                        o3flag,
                        ag.state,
                        perennial,
                        ag.fertflag,
                        soil.getKH2O(),
                        microbe.getNETNMIN(),
                        ag.fertn ); 
  }
  else
  {
    // No crop plants exist - set all monthly fluxes to zero
    
    veg.resetMonthlyFluxes();
  }

  soil.setNINPUT( soil.getNINPUT() + ag.fertn );
  

  // Determine nitrogen losses from ecosystem

  if( 1 == avlnflag )
  {
    soil.updateNLosses( veg.cmnt,
                        (soil.getSURFRUN() + soil.getDRAINAGE()),        
                        pstate[I_AVLN], 
                        soil.getMOIST() );

    if( soil.getNLOST() > pstate[I_AVLN] - veg.getNUPTAKE()
         + microbe.getNETNMIN() + soil.getNINPUT() )
    {
      soil.setNLOST( (pstate[I_AVLN] 
                     - veg.getNUPTAKE() 
                     + microbe.getNETNMIN()
                     + soil.getNINPUT()) );
    }

    if( soil.getNLOST() < ZERO )
    {
      soil.setNLOST( ZERO );
      
      microbe.setNETNMIN( (soil.getNLOST() 
                          + veg.getNUPTAKE() 
                          - soil.getNINPUT()
                          - pstate[I_AVLN]) );
    }                        
  }
  else
  {
    // Do not allow changes in available nitrogen in soil
    //   (used during calibration process)

    soil.setNLOST( (soil.getNINPUT()
                   + microbe.getNETNMIN()
                       - veg.getNUPTAKE()) );
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITTEM44::delta( const int& pdm,
	                    double pstate[],
                      double pdstate[] )
{
  switch( ag.state )
  {
    case 1:  cropDynamics( pstate ); break;
    case 2:  pastureDynamics( pdm, pstate ); break;
    case 3:  urbanDynamics( pdm, pstate ); break;
    default: natvegDynamics( pdm, pstate );
  }


  // Describe monthly changes to carbon pools and fluxes for ODE 
  //   state variables (i.e., pdstate)

  // Carbon pools in ecosystems

  pdstate[I_VEGC] = veg.getGPP() 
                    - veg.getGPR() 
                    - veg.getLTRFALC();

  pdstate[I_SOLC] = veg.getLTRFALC()
                    + ag.getSLASHC()
                    - ag.getSCONVRTFLXC()
                    - microbe.getRH();


  // Nitrogen pools in ecosystems

  pdstate[I_STRN] = veg.getSUPTAKE()
                    - veg.getLTRFALN()
                    - veg.getNRESORB()
                    + veg.getNMOBIL();

  pdstate[I_STON] = veg.getLUPTAKE() 
                    + veg.getNRESORB() 
                    - veg.getNMOBIL();

  pdstate[I_SOLN] = veg.getLTRFALN()
                    + ag.getSLASHN()
                    - ag.getSCONVRTFLXN()
                    - ag.getNSRETENT()
                    - microbe.getNETNMIN();

  pdstate[I_AVLN] = soil.getNINPUT()
                    + microbe.getNETNMIN()
                    - veg.getNUPTAKE()
                    - soil.getNLOST();


  // Phenology

  pdstate[I_UNRMLF] = veg.getUNNORMLEAF();

  pdstate[I_LEAF] = veg.getLEAF();

  pdstate[I_LAI] = veg.getLAI();

  pdstate[I_FPC] = veg.getFPC();

  // Carbon fluxes in ecosystems

  pdstate[I_INGPP] = veg.getINGPP();

  pdstate[I_GPP] = veg.getGPP();
  
  pdstate[I_FOZONE] = veg.getFOZONE();

  pdstate[I_FINDOZONE] = veg.getFINDOZONE();

  pdstate[I_INNPP] = veg.getINNPP();

  pdstate[I_NPP] = veg.getNPP();

  pdstate[I_GPR] = veg.getGPR();

  pdstate[I_RVMNT] = veg.getRMAINT();

  pdstate[I_RVGRW] = veg.getRGRWTH();

  pdstate[I_LTRC] = veg.getLTRFALC();
 
  pdstate[I_RH] = microbe.getRH();
  

  // Nitrogen fluxes in ecosystems

  pdstate[I_NINP] = soil.getNINPUT();

  pdstate[I_AGFRTN] = ag.fertn;

  pdstate[I_INNUP] = veg.getINUPTAKE();

  pdstate[I_VNUP] = veg.getNUPTAKE();

  pdstate[I_VSUP] = veg.getSUPTAKE();
  pdstate[I_VLUP] = veg.getLUPTAKE();
  pdstate[I_VNMBL] = veg.getNMOBIL();
  pdstate[I_VNRSRB] = veg.getNRESORB();

  pdstate[I_LTRN] = veg.getLTRFALN();

  pdstate[I_MNUP] = microbe.getNUPTAKE();

  pdstate[I_NMIN] = microbe.getNETNMIN();

  pdstate[I_NLST] = soil.getNLOST();
  
//  pdstate[I_AGIRRIG] = ag.getIRRIGATE();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

#ifdef CALIBRATE_TEM

void MITTEM44::displayOptionalEflx( const seykey& s )
{
  switch( s )
  {
    case GET_LEAF:    cout << "   LEAF "; break;
    case GET_LAI:     cout << "    LAI "; break;
    case GET_FPC:     cout << "    FPC "; break;

    case GET_INGPP:    cout << "  INGPP "; break;
    case GET_GPP:      cout << "   GPP  "; break;
    case GET_INNPP:    cout << "  INNPP "; break;
    case GET_NPP:      cout << "   NPP  "; break;
    case GET_GPR:      cout << "    RA  "; break;
    case GET_RVMNT:    cout << "  RVMNT "; break;
    case GET_RVGRW:    cout << "  RVGRW "; break;
    case GET_LTRC:     cout << "   LTRC "; break;
    case GET_AGSTUBC:  cout << "AGSTUBC "; break;
    case GET_RH:       cout << "    RH  "; break;
    case GET_NEP:      cout << "   NEP  "; break;

    case GET_D40:      cout << " AOT40  "; break;
    case GET_FOZONE:   cout << " FOZONE "; break;

    case GET_NINP:     cout << " NINPUT "; break;
    case GET_AGFRTN:   cout << "AGFERTN "; break;
    case GET_INNUP:    cout << "  INNUP "; break;
    case GET_VNUP:     cout << " UPTAKE "; break;
    case GET_VSUP:     cout << " SUPTAK "; break;
    case GET_VLUP:     cout << " LUPTAK "; break;
    case GET_VNMBL:    cout << " NMOBIL "; break;
    case GET_VNRSRB:   cout << " NRTRAN "; break;
    case GET_LTRN:     cout << "   LTRN "; break;
    case GET_AGSTUBN:  cout << "AGSTUBN "; break;
    case GET_MNUP:     cout << " NIMMOB "; break;
    case GET_NMIN:     cout << "NETNMIN "; break;
    case GET_NLST:     cout << "  NLOST "; break;

    case GET_CNVRTC:   cout << " CNVRTC "; break;
    case GET_VCNVRTC:  cout << "VCNVRTC "; break;
    case GET_SCNVRTC:  cout << "SCNVRTC "; break;
    case GET_SLASHC:   cout << " SLASHC "; break;
    case GET_CFLX:     cout << "  CFLUX "; break;
    case GET_NCE:     cout << "    NCE  "; break;

    case GET_CNVRTN:   cout << " CNVRTN "; break;
    case GET_VCNVRTN:  cout << "VCNVRTN "; break;
    case GET_SCNVRTN:  cout << "SCNVRTN "; break;
    case GET_SLASHN:   cout << " SLASHN "; break;
    case GET_NRETNT:   cout << " NRETNT "; break;
    case GET_NVRTNT:   cout << " NVRTNT "; break;
    case GET_NSRTNT:   cout << " NSRTNT "; break;

    case GET_AGPRDC:   cout << " AGPRODC "; break;
    case GET_PROD10C:  cout << " PROD10C "; break;
    case GET_PROD100C: cout << "PROD100C "; break;
    case GET_RESIDC:   cout << "  RESIDC "; break;

    case GET_AGPRDN:   cout << " AGPRODN "; break;
    case GET_PROD10N:  cout << " PROD10N "; break;
    case GET_PROD100N: cout << "PROD100N "; break;
    case GET_RESIDN:   cout << "  RESIDN "; break;

    case GET_AGFPRDC:  cout << " AGFPRDC "; break;
    case GET_PRDF10C:  cout << " PRDF10C "; break;
    case GET_PRDF100C: cout << "PRDF100C "; break;
    case GET_FRESIDC:  cout << " FRESIDC "; break;
    case GET_AGPRDFC:  cout << " AGPRDFC "; break;
    case GET_PRD10FC:  cout << " PRD10FC "; break;
    case GET_PRD100FC: cout << "PRD100FC "; break;
    case GET_TOTPRDFC: cout << "TOTPRDFC "; break;
    case GET_RESIDFC:  cout << " RESIDFC "; break;

    case GET_AGFPRDN:  cout << " AGFPRDN "; break;
    case GET_PRDF10N:  cout << " PRDF10N "; break;
    case GET_PRDF100N: cout << "PRDF100N "; break;
    case GET_FRESIDN:  cout << " FRESIDN "; break;
    case GET_AGPRDFN:  cout << " AGPRDFN "; break;
    case GET_PRD10FN:  cout << " PRD10FN "; break;
    case GET_PRD100FN: cout << "PRD100FN "; break;
    case GET_TOTPRDFN: cout << "TOTPRDFN "; break;
    case GET_RESIDFN:  cout << " RESIDFN "; break;

    case GET_L2SN:    cout << "   LCON "; break;
  }

};

#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

#ifdef CALIBRATE_TEM

void MITTEM44::displayOptionalWflx( const swykey& s )
{
  switch( s )
  {
    case GET_SH2O:    cout << " SMOIST"; break;
    case GET_PCTP:    cout << "  PCTP "; break;
    case GET_VSM:     cout << "  VSM  "; break;

    case GET_RAIN:    cout << "  RAIN "; break;
    case GET_SNWFAL:  cout << " SNWFAL"; break;
    case GET_SNWINF:  cout << " SNWINF"; break;
    case GET_AGIRRIG: cout << " IRRIG "; break;
    case GET_PET:     cout << "  PET  "; break;
    case GET_INEET:   cout << "  INEET"; break;
    case GET_EET:     cout << "  EET  "; break;
    case GET_RPERC:   cout << "  RPERC"; break;
    case GET_SPERC:   cout << "  SPERC"; break;
    case GET_RRUN:    cout << "  RRUN "; break;
    case GET_SRUN:    cout << "  SRUN "; break;
    case GET_WYLD:    cout << "  WYLD "; break;
  }

};

#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int MITTEM44::ecdqc( const int& dcmnt )
{

  int qc = ACCEPT;

  if( vegca[dcmnt] <= -9999.99 ) { return qc = 101; }
  if( vegcb[dcmnt] <= -9999.99 ) { return qc = 102; }

  if( strna[dcmnt] <= -9999.99 ) { return qc = 103; }
  if( strnb[dcmnt] <= -9999.99 ) { return qc = 104; }

  if( stona[dcmnt] <= -9999.99 ) { return qc = 105; }
  if( stonb[dcmnt] <= -9999.99 ) { return qc = 106; }

  if( solca[dcmnt] <= -9999.99 ) { return qc = 107; }
  if( solcb[dcmnt] <= -9999.99 ) { return qc = 108; }

  if( solna[dcmnt] <= -9999.99 ) { return qc = 109; }
  if( solnb[dcmnt] <= -9999.99 ) { return qc = 110; }

  if( avlna[dcmnt] <= -9999.99 ) { return qc = 111; }
  if( avlnb[dcmnt] <= -9999.99 ) { return qc = 112; }

  if( veg.getUNLEAF12( dcmnt ) <= -9.99 ) { return qc = 113; }
  if( veg.getINITLEAFMX( dcmnt ) <= -9.99 ) { return qc = 114; }

  if( veg.getCMAXCUT( dcmnt ) <= -99.99 ) { return qc = 115; }
  if( veg.getCMAX1A( dcmnt ) <= -9999.99 ) { return qc = 116; }
  if( veg.getCMAX1B( dcmnt ) <= -9999.99 ) { return qc = 117; }
  if( veg.getCMAX2A( dcmnt ) <= -9999.99 ) { return qc = 118; }
  if( veg.getCMAX2B( dcmnt ) <= -9999.99 ) { return qc = 119; }
  if( veg.getCFALL( dcmnt ) <= -99.99 ) { return qc = 120; }

  if( veg.getKRA( dcmnt ) <= -99.99 ) { return qc = 121; }
  if( veg.getKRB( dcmnt ) <= -99.99 ) { return qc = 122; }

  if( microbe.getKDA( dcmnt ) <= -99.99 ) { return qc = 123; }
  if( microbe.getKDB( dcmnt ) <= -99.99 ) { return qc = 124; }

  if( microbe.getLCCLNC( dcmnt ) <= -99.99 ) { return qc = 125; }
  if( microbe.getPROPFTOS( dcmnt ) <= -99.99 ) { return qc = 126; }

  if( veg.getNMAXCUT( dcmnt ) <= -99.99 ) { return qc = 127; }
  if( veg.getNMAX1A( dcmnt ) <= -9999.99 ) { return qc = 128; }
  if( veg.getNMAX1B( dcmnt ) <= -9999.99 ) { return qc = 129; }
  if( veg.getNMAX2A( dcmnt ) <= -9999.99 ) { return qc = 130; }
  if( veg.getNMAX2B( dcmnt ) <= -9999.99 ) { return qc = 131; }

  if( veg.getNFALL( dcmnt ) <= -99.99 ) { return qc = 132; }

  if( microbe.getNFIXPAR( dcmnt ) <= -99.99 ) { return qc = 133; }

  if( microbe.getNUPA( dcmnt ) <= -9999.99 ) { return qc = 134; }
  if( microbe.getNUPB( dcmnt ) <= -9999.99 ) { return qc = 135; }

  if( soil.getNLOSS( dcmnt ) <= -99.99 ) { return qc = 136; }
  
  if( veg.getINITCNEVEN( dcmnt ) <= -9999.99 ) { return qc = 137; }
  if( veg.getCNMIN( dcmnt ) <= -9999.99 ) { return qc = 138; }
  if( veg.getC2NA( dcmnt ) <= -9999.99 ) { return qc = 139; }
  if( veg.getC2NB( dcmnt ) <= -9999.99 ) { return qc = 140; }
  if( veg.getC2NMIN( dcmnt ) <= -9999.99 ) { return qc = 141; }

  if( microbe.getCNSOIL( dcmnt ) <= -9999.99 ) { return qc = 142; }

  if( veg.getO3PARA( dcmnt ) <= -99.99 ) { return qc = 193; }
  if( veg.getO3PARB( dcmnt ) <= -99.99 ) { return qc = 194; }
  if( veg.getO3PARC( dcmnt ) <= -99.99 ) { return qc = 195; }

  return qc;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void MITTEM44::ECDsetODEstate( const int& pdcmnt,
                             const double& psiplusc )
{
  // Initialize the NUMEQ state variables used in the 
  //   ODE integrator from ECD and DAT files 
    
  y[I_VEGC] = vegca[pdcmnt] * psiplusc + vegcb[pdcmnt];
  
  if( y[I_VEGC] < ZERO ) { y[I_VEGC] = ZERO; }


  y[I_SOLC] = solca[pdcmnt] * psiplusc + solcb[pdcmnt];
  
  if( y[I_SOLC] < ZERO ) { y[I_SOLC] = ZERO; }


  y[I_STRN] = strna[pdcmnt] * psiplusc + strnb[pdcmnt];
  
  if( y[I_STRN] < ZERO ) { y[I_STRN] = ZERO; }
    

  y[I_STON] = stona[pdcmnt] * psiplusc + stonb[pdcmnt];

  if( y[I_STON] < ZERO ) { y[I_STON] = ZERO; }
    

  y[I_SOLN] = solna[pdcmnt] * psiplusc + solnb[pdcmnt];
    
  if( y[I_SOLN] < ZERO ) { y[I_SOLN] = ZERO; }
        
   
  y[I_AVLN] = avlna[pdcmnt] * psiplusc + avlnb[pdcmnt];
  
  if( y[I_AVLN] < ZERO ) { y[I_AVLN] = ZERO; }
   

  // Initialize all phenology and flux states to zero
  
  resetODEflux();
  
  // Reinitialize phenology state variables
  
  y[I_UNRMLF] = veg.getUNLEAF12( pdcmnt );

  y[I_LEAF] = veg.getINITLEAFMX( pdcmnt );
    
};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void MITTEM44::getcropecd( const int& dv, const string& agecd )
{
  string agversion;
  string agsitename;
  string agdeveloper;
  string agsitecol;
  string agsiterow;
  string agupdated;

  string agdescription;

  fecd[dv].open( agecd.c_str(), ios::in );

  if( !fecd[dv] )
  {
    cerr << endl << "Cannot open " << agecd;
    cerr << " for ag siteECD input" << endl;
    
    exit( -1 );
  }

  veg.getXMLsiteRootNode( fecd[dv],
                          "siteECD",
                          agversion,
                          agsitename,
                          agsitecol,
                          agsiterow,
                          agdeveloper,
                          agupdated );

  ag.cmnt = veg.getXMLsiteCommunityNode( fecd[dv],
                                         "siteECD",
                                         agdescription );

  vegca[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "vegca",
                                              ag.cmnt );

  vegcb[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "vegcb",
                                              ag.cmnt );

  strna[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "strna",
                                              ag.cmnt );

  strnb[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "strnb",
                                              ag.cmnt );

  solca[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "solca",
                                              ag.cmnt );

  solcb[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "solcb",
                                              ag.cmnt );

  solna[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "solna",
                                              ag.cmnt );

  solnb[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "solnb",
                                              ag.cmnt );

  avlna[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "avlna",
                                              ag.cmnt );

  avlnb[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "avlnb",
                                              ag.cmnt );

  stona[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "stona",
                                              ag.cmnt );

  stonb[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "stonb",
                                              ag.cmnt );

  veg.setUNLEAF12( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "unleaf12",
                                              ag.cmnt ),
                   ag.cmnt );

  veg.setINITLEAFMX( veg.getXMLcmntArrayDouble( fecd[dv],
                                                "siteECD",
                                                "initleafmx",
                                                 ag.cmnt ),
                     ag.cmnt );


  veg.setCMAXCUT( veg.getXMLcmntArrayDouble( fecd[dv],
                                             "siteECD",
                                             "vegcmaxcut",
                                             ag.cmnt ),
                  ag.cmnt );

  veg.setCMAX1A( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegcmax1a",
                                            ag.cmnt ),
                 ag.cmnt );

  veg.setCMAX1B( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegcmax1b",
                                            ag.cmnt ),
                 ag.cmnt );

  veg.setCMAX2A( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegcmax2a",
                                            ag.cmnt ),
                 ag.cmnt );

  veg.setCMAX2B( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegcmax2b",
                                            ag.cmnt ),
                 ag.cmnt );

  veg.setCFALL( veg.getXMLcmntArrayDouble( fecd[dv],
                                           "siteECD",
                                           "vegcfall",
                                           ag.cmnt ),
                ag.cmnt );

  veg.setKRA( veg.getXMLcmntArrayDouble( fecd[dv],
                                         "siteECD",
                                         "vegkra",
                                         ag.cmnt ),
              ag.cmnt );

  veg.setKRB( veg.getXMLcmntArrayDouble( fecd[dv],
                                         "siteECD",
                                         "vegkrb",
                                         ag.cmnt ),
              ag.cmnt );

  microbe.setKDA( veg.getXMLcmntArrayDouble( fecd[dv],
                                             "siteECD",
                                             "microbekda",
                                             ag.cmnt ),
                  ag.cmnt );

  microbe.setKDB( veg.getXMLcmntArrayDouble( fecd[dv],
                                             "siteECD",
                                             "microbekdb",
                                             ag.cmnt ),
                  ag.cmnt );

  microbe.setLCCLNC( veg.getXMLcmntArrayDouble( fecd[dv],
                                                "siteECD",
                                                "microbelcclnc",
                                                ag.cmnt ),
                     ag.cmnt );

  microbe.setPROPFTOS( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "microbepropftos",
                                                  ag.cmnt ),
                       ag.cmnt );

  veg.setNMAXCUT( veg.getXMLcmntArrayDouble( fecd[dv],
                                             "siteECD",
                                             "vegnmaxcut",
                                             ag.cmnt ),
                  ag.cmnt );

  veg.setNMAX1A( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegnmax1a",
                                            ag.cmnt ),
                 ag.cmnt );

  veg.setNMAX1B( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegnmax1b",
                                            ag.cmnt ),
                 ag.cmnt );

  veg.setNMAX2A( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegnmax2a",
                                            ag.cmnt ),
                 ag.cmnt );

  veg.setNMAX2B( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegnmax2b",
                                            ag.cmnt ),
                 ag.cmnt );

  veg.setNFALL( veg.getXMLcmntArrayDouble( fecd[dv],
                                           "siteECD",
                                           "vegnfall",
                                           ag.cmnt ),
                ag.cmnt );

  microbe.setNUPA( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "microbenupa",
                                              ag.cmnt ),
                   ag.cmnt );

  microbe.setNUPB( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "microbenupb",
                                              ag.cmnt ),
                   ag.cmnt );

  soil.setNLOSS( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "soilnloss",
                                            ag.cmnt ),
                 ag.cmnt );

  microbe.setNFIXPAR( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "microbenfixpar",
                                                  ag.cmnt ),
                      ag.cmnt );

  veg.setINITCNEVEN( veg.getXMLcmntArrayDouble( fecd[dv],
                                                "siteECD",
                                                "veginitcneven",
                                                ag.cmnt ),
                     ag.cmnt );

  veg.setCNMIN( veg.getXMLcmntArrayDouble( fecd[dv],
                                           "siteECD",
                                           "vegcnmin",
                                           ag.cmnt ),
                ag.cmnt );

  veg.setC2NA( veg.getXMLcmntArrayDouble( fecd[dv],
                                          "siteECD",
                                          "vegc2na",
                                          ag.cmnt ),
               ag.cmnt );

  veg.setC2NB( veg.getXMLcmntArrayDouble( fecd[dv],
                                          "siteECD",
                                          "vegc2nb",
                                          ag.cmnt ),
               ag.cmnt );

  veg.setC2NMIN( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegc2nmin",
                                            ag.cmnt ),
                 ag.cmnt );

  microbe.setCNSOIL( veg.getXMLcmntArrayDouble( fecd[dv],
                                                "siteECD",
                                                "microbecnsoil",
                                                ag.cmnt ),
                     ag.cmnt );

  veg.setO3PARA( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "o3para",
                                            ag.cmnt ),
                 ag.cmnt );

  veg.setO3PARB( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "o3parb",
                                            ag.cmnt ),
                 ag.cmnt );

  veg.setO3PARC( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "o3parc",
                                            ag.cmnt ),
                 ag.cmnt );

  veg.endXMLcommunityNode( fecd[dv] );

  fecd[dv].close();

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void MITTEM44::getenviron( const int& pdm )
{
  // Determine monthly potential evapotranspiration

  atms.petjh( atms.getNIRR(), atms.getTAIR(), pdm );

  if( atms.getPET() < soil.getINEET() )
  {
    atms.setPET( soil.getINEET() );
  }


  if( soil.getSNOWPACK() < ZERO )
  {
    soil.setSNOWPACK( ZERO );
  }
  
};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

#ifdef CALIBRATE_TEM

double MITTEM44::getOptionalEflx( const int& optflx )
{

  double outflux;

  switch( optflx )
  {
    case GET_LEAF:    outflux = y[I_LEAF]; break;
    case GET_LAI:     outflux = y[I_LAI]; break;
    case GET_FPC:     outflux = y[I_FPC]; break;

    case GET_INGPP:    outflux = y[I_INGPP]; break;
    case GET_GPP:      outflux = y[I_GPP]; break;
    case GET_INNPP:    outflux = y[I_INNPP]; break;
    case GET_NPP:      outflux = y[I_NPP]; break;
    case GET_GPR:      outflux = y[I_GPR]; break;
    case GET_RVMNT:    outflux = y[I_RVMNT]; break;
    case GET_RVGRW:    outflux = y[I_RVGRW]; break;
    case GET_LTRC:     outflux = y[I_LTRC]; break;
    case GET_RH:       outflux = y[I_RH]; break;
    case GET_NEP:      outflux = nep; break;

    case GET_D40:      outflux = atms.getAOT40(); break;
    case GET_FOZONE:   outflux = y[I_FOZONE]; break;

    case GET_NINP:     outflux = y[I_NINP]; break;
    case GET_AGFRTN:   outflux = y[I_AGFRTN]; break;
    case GET_INNUP:    outflux = y[I_INNUP]; break;
    case GET_VNUP:     outflux = y[I_VNUP]; break;
    case GET_VSUP:     outflux = y[I_VSUP]; break;
    case GET_VLUP:     outflux = y[I_VLUP]; break;
    case GET_VNMBL:    outflux = y[I_VNMBL]; break;
    case GET_VNRSRB:   outflux = y[I_VNRSRB]; break;
    case GET_LTRN:     outflux = y[I_LTRN]; break;
    case GET_AGSTUBN:  outflux = ag.getSTUBBLEN(); break;
    case GET_MNUP:     outflux = y[I_MNUP]; break;
    case GET_NMIN:     outflux = y[I_NMIN]; break;
    case GET_NLST:     outflux = y[I_NLST]; break;
    
    case GET_CNVRTC:   outflux = ag.getCONVRTFLXC();  break;
    case GET_VCNVRTC:  outflux = ag.getVCONVRTFLXC();  break;
    case GET_SCNVRTC:  outflux = ag.getSCONVRTFLXC();  break;
    case GET_SLASHC:   outflux = ag.getSLASHC();  break;
    case GET_CFLX:     outflux = ag.getCFLUX();  break;
    case GET_NCE:      outflux = nce;  break;

    case GET_CNVRTN:   outflux = ag.getCONVRTFLXN();  break;
    case GET_VCNVRTN:  outflux = ag.getVCONVRTFLXN();  break;
    case GET_SCNVRTN:  outflux = ag.getSCONVRTFLXN();  break;
    case GET_SLASHN:   outflux = ag.getSLASHN();  break;
    case GET_NRETNT:   outflux = ag.getNRETENT();  break;
    case GET_NVRTNT:   outflux = ag.getNVRETENT();  break;
    case GET_NSRTNT:   outflux = ag.getNSRETENT();  break;

    case GET_AGSTUBC:  outflux = ag.getSTUBBLEC(); break;
    case GET_RESIDC:   outflux = ag.getCROPRESIDUEC();  break;

    case GET_AGPRDC:   outflux = ag.getPROD1C();  break;
    case GET_PROD10C:  outflux = ag.getPROD10C();  break;
    case GET_PROD100C: outflux = ag.getPROD100C();  break;

    case GET_AGPRDN:   outflux = ag.getPROD1N();  break;
    case GET_PROD10N:  outflux = ag.getPROD10N();  break;
    case GET_PROD100N: outflux = ag.getPROD100N();  break;
    case GET_RESIDN:   outflux = ag.getCROPRESIDUEN();  break;

    case GET_FRESIDC:  outflux = ag.getFORMCROPRESIDUEC();  break;

    case GET_AGFPRDC:  outflux = ag.getCROPPRODC();  break;
    case GET_PRDF10C:  outflux = ag.getFORMPROD10C();  break;
    case GET_PRDF100C: outflux = ag.getFORMPROD100C();  break;
    case GET_AGPRDFC:  outflux = ag.getPROD1DECAYC();  break;
    case GET_PRD10FC:  outflux = ag.getPROD10DECAYC();  break;
    case GET_PRD100FC: outflux = ag.getPROD100DECAYC();  break;
    case GET_TOTPRDFC: outflux = ag.getTOTPRODDECAYC();  break;
    case GET_RESIDFC:  outflux = ag.getCROPRESIDUEFLXC();  break;

    case GET_AGFPRDN:  outflux = ag.getCROPPRODN();  break;
    case GET_PRDF10N:  outflux = ag.getFORMPROD10N();  break;
    case GET_PRDF100N: outflux = ag.getFORMPROD100N();  break;
    case GET_FRESIDN:  outflux = ag.getFORMCROPRESIDUEN();  break;
    case GET_AGPRDFN:  outflux = ag.getPROD1DECAYN();  break;
    case GET_PRD10FN:  outflux = ag.getPROD10DECAYN();  break;
    case GET_PRD100FN: outflux = ag.getPROD100DECAYN();  break;
    case GET_TOTPRDFN: outflux = ag.getTOTPRODDECAYN();  break;
    case GET_RESIDFN:  outflux = ag.getCROPRESIDUEFLXN();  break;

    case GET_L2SN:     if ( y[I_VEGC] != ZERO )
                       {
                         outflux = y[I_STON]/y[I_STRN];
                       }
                       else { outflux = MISSING; } break;

    default:           outflux = MISSING;
  }

  return outflux;

};

#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

#ifdef CALIBRATE_TEM

double MITTEM44::getOptionalWflx( const int& optflx )
{

  double outflux;

  switch( optflx )
  {
    case GET_SH2O:    outflux = soil.getMOIST(); break;
    case GET_PCTP:    outflux = soil.getPCTP(); break;
    case GET_VSM:     outflux = 100.0*soil.getVSM(); break;

    case GET_RAIN:    outflux = atms.getRAIN(); break;
    case GET_SNWFAL:  outflux = atms.getSNOWFALL(); break;
    case GET_SNWINF:  outflux = soil.getSNOWINF(); break;
    case GET_AGIRRIG: outflux = ZERO; break;
    case GET_PET:     outflux = atms.getPET(); break;
    case GET_INEET:   outflux = soil.getINEET(); break;
    case GET_EET:     outflux = soil.getEET(); break;
    case GET_RPERC:   outflux = ZERO; break;
    case GET_SPERC:   outflux = ZERO; break;
    case GET_RRUN:    outflux = ZERO; break;
    case GET_SRUN:    outflux = ZERO; break;
    case GET_WYLD:    outflux = soil.getH2OYLD(); break;
    default:          outflux = MISSING;
  }

  return outflux;

};

#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITTEM44::getsitecd( const int& numcmnt, ofstream& rflog1 )
{

  int dv;
  string ecd;

  cout << "Enter name of the site (.ECD) data file with the";
  cout << " parameter values  for cmnt" << endl;

  rflog1 << "Enter name of the site (.ECD) data file with the";
  rflog1 << " parameter values cmnt" << endl;

  for( dv = 0; dv < numcmnt; ++dv )
  {
    cout << (dv+1) << ": ";

    cin >> ecd;

    rflog1 << (dv+1) << ": " << ecd << endl;

    getsitecd( dv, ecd );
  }

  rflog1 << endl;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void MITTEM44::getsitecd( const int& dv, const string&  ecd )
{
  string kstring;

  fecd[dv].open( ecd.c_str(), ios::in );

  if( !fecd[dv] )
  {
    cerr << endl;
    cerr << "Cannot open " << ecd << " for site ECD input";
    cerr << endl;
    
    exit( -1 );
  }

  veg.getXMLsiteRootNode( fecd[dv],
                          "siteECD",
                          version,
                          sitename,
                          kstring,
                          kstring,
                          developer,
                          kstring );

  veg.cmnt = veg.getXMLsiteCommunityNode( fecd[dv],
                                          "siteECD",
                                          description );

  vegca[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "vegca",
                                               veg.cmnt );

  vegcb[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "vegcb",
                                               veg.cmnt );

  strna[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "strna",
                                               veg.cmnt );

  strnb[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "strnb",
                                               veg.cmnt );

  solca[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "solca",
                                               veg.cmnt );

  solcb[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "solcb",
                                               veg.cmnt );

  solna[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "solna",
                                               veg.cmnt );

  solnb[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "solnb",
                                               veg.cmnt );

  avlna[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "avlna",
                                               veg.cmnt );

  avlnb[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "avlnb",
                                               veg.cmnt );

  stona[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "stona",
                                               veg.cmnt );

  stonb[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "stonb",
                                               veg.cmnt );

  veg.setUNLEAF12( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "unleaf12",
                                              veg.cmnt ),
                   veg.cmnt );

  veg.setINITLEAFMX( veg.getXMLcmntArrayDouble( fecd[dv],
                                                "siteECD",
                                                "initleafmx",
                                                veg.cmnt ),
                     veg.cmnt );

  veg.setCMAXCUT( veg.getXMLcmntArrayDouble( fecd[dv],
                                             "siteECD",
                                             "vegcmaxcut",
                                             veg.cmnt ),
                  veg.cmnt );

  veg.setCMAX1A( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegcmax1a",
                                            veg.cmnt ),
                 veg.cmnt );

  veg.setCMAX1B( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegcmax1b",
                                            veg.cmnt ),
                 veg.cmnt );

  veg.setCMAX2A( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegcmax2a",
                                            veg.cmnt ),
                 veg.cmnt );

  veg.setCMAX2B( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegcmax2b",
                                            veg.cmnt ),
                 veg.cmnt );

  veg.setCFALL( veg.getXMLcmntArrayDouble( fecd[dv],
                                           "siteECD",
                                           "vegcfall",
                                           veg.cmnt ),
                veg.cmnt );

  veg.setKRA( veg.getXMLcmntArrayDouble( fecd[dv],
                                         "siteECD",
                                         "vegkra",
                                         veg.cmnt ),
              veg.cmnt );

  veg.setKRB( veg.getXMLcmntArrayDouble( fecd[dv],
                                         "siteECD",
                                         "vegkrb",
                                         veg.cmnt ),
              veg.cmnt );

  microbe.setKDA( veg.getXMLcmntArrayDouble( fecd[dv],
                                             "siteECD",
                                             "microbekda",
                                             veg.cmnt ),
                  veg.cmnt );

  microbe.setKDB( veg.getXMLcmntArrayDouble( fecd[dv],
                                             "siteECD",
                                             "microbekdb",
                                             veg.cmnt ),
                  veg.cmnt );

  microbe.setLCCLNC( veg.getXMLcmntArrayDouble( fecd[dv],
                                                "siteECD",
                                                "microbelcclnc",
                                                veg.cmnt ),
                     veg.cmnt );

  microbe.setPROPFTOS( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "microbepropftos",
                                                  veg.cmnt ),
                       veg.cmnt );

  veg.setNMAXCUT( veg.getXMLcmntArrayDouble( fecd[dv],
                                             "siteECD",
                                             "vegnmaxcut",
                                             veg.cmnt ),
                  veg.cmnt );

  veg.setNMAX1A( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegnmax1a",
                                            veg.cmnt ),
                 veg.cmnt );

  veg.setNMAX1B( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegnmax1b",
                                            veg.cmnt ),
                 veg.cmnt );

  veg.setNMAX2A( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegnmax2a",
                                            veg.cmnt ),
                 veg.cmnt );

  veg.setNMAX2B( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegnmax2b",
                                            veg.cmnt ),
                 veg.cmnt );

  veg.setNFALL( veg.getXMLcmntArrayDouble( fecd[dv],
                                           "siteECD",
                                           "vegnfall",
                                           veg.cmnt ),
                veg.cmnt );

  microbe.setNUPA( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "microbenupa",
                                              veg.cmnt ),
                   veg.cmnt );

  microbe.setNUPB( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "microbenupb",
                                              veg.cmnt ),
                   veg.cmnt );

  soil.setNLOSS( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "soilnloss",
                                            veg.cmnt ),
                 veg.cmnt );

  microbe.setNFIXPAR( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "microbenfixpar",
                                                  veg.cmnt ),
                      veg.cmnt );

  veg.setINITCNEVEN( veg.getXMLcmntArrayDouble( fecd[dv],
                                                "siteECD",
                                                "veginitcneven",
                                                veg.cmnt ),
                     veg.cmnt );

  veg.setCNMIN( veg.getXMLcmntArrayDouble( fecd[dv],
                                           "siteECD",
                                           "vegcnmin",
                                           veg.cmnt ),
                veg.cmnt );

  veg.setC2NA( veg.getXMLcmntArrayDouble( fecd[dv],
                                          "siteECD",
                                          "vegc2na",
                                          veg.cmnt ),
               veg.cmnt );

  veg.setC2NB( veg.getXMLcmntArrayDouble( fecd[dv],
                                          "siteECD",
                                          "vegc2nb",
                                          veg.cmnt ),
               veg.cmnt );

  veg.setC2NMIN( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegc2nmin",
                                            veg.cmnt ),
                 veg.cmnt );

  microbe.setCNSOIL( veg.getXMLcmntArrayDouble( fecd[dv],
                                                "siteECD",
                                                "microbecnsoil",
                                                veg.cmnt ),
                     veg.cmnt );

  veg.setO3PARA( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "o3para",
                                            veg.cmnt ),
                 veg.cmnt );

  veg.setO3PARB( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "o3parb",
                                            veg.cmnt ),
                 veg.cmnt );

  veg.setO3PARC( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "o3parc",
                                            veg.cmnt ),
                 veg.cmnt );

  veg.endXMLcommunityNode( fecd[dv] );

  fecd[dv].close();

};

/* *************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void MITTEM44::initrun( ofstream& rflog1 )
{
  
  avlnflag = nfeed = rheqflag = 0;

/* **************************************************************
		  Run Model with Nitrogen Limitation?
************************************************************** */

  #ifdef CALIBRATE_TEM
    cout << "Do you want to start calibration allowing available ";
    cout << "N to fluctuate?  ";

    if( 'Y' == toupper( getch() ) )
    {
      avlnflag = 1;

      cout << "Y" << endl;
    }
    else
    {
      avlnflag = 0;

      cout << "N" << endl;
    }
  #else
    cout << endl << "Do you want to allow available N to fluctuate?";
    cout << endl;
    cout << "  Enter 0 for No" << endl;
    cout << "  Enter 1 for Yes: ";

    cin >> avlnflag;
  
    rflog1 << endl << "Do you want to allow available N to fluctuate?";
    rflog1 << endl;
    rflog1 << "  Enter 0 for No" << endl;
    rflog1 << "  Enter 1 for Yes: " << endl;
    rflog1 << "avlnflag = " << avlnflag << endl << endl;
  #endif
  
  baseline = initbase = 0;
  
  #ifdef CALIBRATE_TEM
    cout << "Do you want to start calibration with N feedback ";
    cout << "on GPP?  ";
    
    if( 'Y' == toupper( getch() ) )
    {
      nfeed = 1;

      cout << "Y" << endl;
      cout << "Do you want to solve for baseline soil nitrogen?  ";

      if ( 'Y' == toupper( getch() ) )
      {
         baseline = 1;

         cout << "Y" << endl;
      }
      else { cout << "N" << endl; }
    }
    else
    {
       nfeed = 0;

       cout << "N" << endl;
    }
  #else
    cout << endl << "Do you want nitrogen feedback on GPP?" << endl;
    cout << "  Enter 0 for No" << endl;
    cout << "  Enter 1 for Yes: ";

    cin >> nfeed;

    rflog1 << endl << "Do you want nitrogen feedback on GPP?";
    rflog1 << endl;
    rflog1 << "  Enter 0 for No" << endl;
    rflog1 << "  Enter 1 for Yes: " << endl;
    rflog1 << "nfeed = " << nfeed << endl << endl;

    if( 1 == nfeed )
    {
      cout << endl;
      cout << "Do you want to solve for baseline soil nitrogen?";
      cout<< endl;
      cout << "  Enter 0 for No" << endl;
      cout << "  Enter 1 for Yes: ";

      cin >> initbase;

      baseline = initbase;

      rflog1 << endl;
      rflog1 << "Do you want to solve for baseline soil nitrogen?";
      rflog1 << endl;
      rflog1 << "  Enter 0 for No" << endl;
      rflog1 << "  Enter 1 for Yes: " << endl;
      rflog1 << "baseline = " << baseline << endl << endl;
    }
  #endif
  
    
/* **************************************************************
	     Run Model with Moisture Limitation?
************************************************************** */

  moistlim = 0;
  
  #ifdef CALIBRATE_TEM
    cout << "Do you want to start calibration with moisture ";
    cout << "limitation?  ";
  
    if( 'Y' == toupper( getch() ) )
    {
      moistlim = 1;

      cout << "Y" << endl;
    }
    else
    {
      moistlim = 0;

      cout << "N" << endl;
    }
  #else
    cout << endl;
    cout << "Do you want to run the model with moisture limitation?";
    cout << endl;
    cout << "  Enter 0 for No" << endl;
    cout << "  Enter 1 for Yes: ";

    cin >> moistlim;

    rflog1 << endl;
    rflog1 << "Do you want to run the model with moisture limitation?";
    rflog1 << endl;
    rflog1 << "  Enter 0 for No" << endl;
    rflog1 << "  Enter 1 for Yes: " << endl;
    rflog1 << "moistlim = " << moistlim << endl << endl;
  #endif
  

/****************************************************************
			 Run Model with Ozone?
************************************************************** */

  o3flag = 0;
  
  #ifdef CALIBRATE_TEM
    cout << "Do you want to start calibration with ozone?  ";
    
    if( 'Y' == toupper(getch()) )
    {
      o3flag = 1;

      cout << "Y" << endl;
    }
    else
    {
      o3flag = 0;

      cout << "N" << endl;
    }  
  #else
    cout << endl;
    cout << "Do you want to run the model with ozone?" << endl;
    cout << "  Enter 0 for No" << endl;
    cout << "  Enter 1 for Yes: ";

    cin >> o3flag;

    rflog1 << endl;
    rflog1 << "Do you want to run the model with ozone?" << endl;
    rflog1 << "  Enter 0 for No" << endl;
    rflog1 << "  Enter 1 for Yes: " << endl;
    rflog1 << "o3flag = " << o3flag << endl << endl;
  #endif


/* ***************************************************************
	       Details for Steady State Conditions
************************************************************** */

  #ifdef CALIBRATE_TEM
    double co2level;
    double dc2n;

    equil = 0;
    adapttol = 0;
    intbomb = 0;
    tolbomb = 0;

    cout << "Do you want the model to stop at steady state ";
    cout << "conditions?  ";
    
    if( 'Y' == toupper( getch() ) )
    {
      equil = 1;
      
      cout << "Y" << endl;
      cout << endl;
      cout << "How many years do you want to wait before checking ";
      cout << "equilibrium conditions?  ";

      cin >> strteq;

      strteq *= 12;

      if( 0 == nfeed )
      {
        cout << "Do you want decomposition to come into ";
        cout << "equilibrium?  ";
       
        if( 'Y' == toupper(getch()) )
        {
	        rheqflag = 1;
	        
	        cout << "Y" << endl;
        }
        else { cout << "N" << endl; }
      }

      cout << "Enter the absolute tolerance for the water cycle?  ";

      cin >> wtol;

      cout << "Enter the absolute tolerance for the carbon cycle?  ";

      cin >> ctol;

      if( 1 == nfeed )
      {
        rheqflag = 1;
        
        cout << "Enter the absolute tolerance for the nitrogen ";
        cout << "cycle?  ";

        cin >> ntol;
      }

      cout << "Do you want to have the integrator tolerance ";
      cout << "adjusted automatically?  ";
      
      if( 'Y' == toupper( getch() ) )
      {
        cout << "Y" << endl;
        
        adapttol = 1;
        
        askODE( rflog1 );
        
        cout << "Enter the maximum number of years for the model ";
        cout << "to run:  ";

        cin >> maxyears;

        cout << "Enter the maximum number of attempts to reach a ";
        cout << "solution:  ";

        cin >> maxnrun;
      }
      else { cout << "N" << endl; }
    }
    else { cout << "N" << endl << endl; }

    if ( 0 == adapttol ) { askODE( rflog1 ); }

    cout << endl << "Enter the level of carbon dioxide in ppmv: ";
    
    cin >> co2level;
    
    atms.setCO2LEVEL( co2level );

    cout << endl << endl << endl;
    cout << "Enter the factor for changing C:N per ppmv of ";
    cout << "enhanced CO2: " << endl;
    cout << "               (Enter 0.0 for no change)" << endl;
    
    cin >> dc2n;
    
    veg.setDC2N( dc2n );
    
    cout << endl;
  #else
    maxyears = 0;

    maxnrun = 0;

    cout << endl;
    cout << "How many years do you want to wait before checking";
    cout << " equilibrium conditions? ";

    cin >> strteq;

    rflog1 << endl;
    rflog1 << "How many years do you want to wait before checking";
    rflog1 << " equilibrium conditions? ";
    rflog1 << endl;
    rflog1 << "strteq = " << strteq << endl << endl;

    cout << endl;
    cout << "Enter the maximum number of years for the model to run:";

    cin >> maxyears;

    rflog1 << endl;
    rflog1 << "Enter the maximum number of years for the model to run: ";
    rflog1 << endl;
    rflog1 << "maxyears = " << maxyears << endl << endl;

    runsize = maxyears;

    cout << endl;
    cout << "Enter the maximum number of attempts to reach a solution: ";

    cin >> maxnrun;

    rflog1 << endl;
    rflog1 << "Enter the maximum number of attempts to reach a solution: ";
    rflog1 << endl;
    rflog1 << "maxnrun = " << maxnrun << endl << endl;

    if( 0 == nfeed )
    {
      cout << endl;
      cout << "Do you want decomposition to come into equilibrium? ";
      cout << "  Enter 0 for No" << endl;
      cout << "  Enter 1 for Yes: ";

      cin >> rheqflag;

      rflog1 << endl;
      rflog1 << "Do you want decomposition to come into equilibrium? ";
      rflog1 << endl;
      rflog1 << "  Enter 0 for No" << endl;
      rflog1 << "  Enter 1 for Yes: " << endl;
      rflog1 << "rheqflag = " << rheqflag << endl << endl;
    }
   
    wtol = 1000.0;

    cout << endl;
    cout << "What absolute tolerance do you want to use for";
    cout << " checking equilibrium";
    cout << endl;
    cout << "of the water cycle? ";

    cin >> wtol;

    rflog1 << endl;
    rflog1 << "What absolute tolerance do you want to use for";
    rflog1 << " checking equilibrium";
    rflog1 << endl;
    rflog1 << "of the water cycle? wtol = " << wtol;
    rflog1 << endl << endl;

    ctol = 1000.0;

    cout << endl;
    cout << "What absolute tolerance do you want to use for";
    cout << " checking equilibrium";
    cout << endl;
    cout << "of the carbon cycle? ";

    cin >> ctol;

    rflog1 << endl;
    rflog1 << "What absolute tolerance do you want to use for";
    rflog1 << " checking equilibrium";
    rflog1 << endl;
    rflog1 << "of the carbon cycle?" << endl;
    rflog1 << "ctol = " << ctol << endl << endl;

    ntol = 1000.0;

    if( 1 == nfeed )
    {
      rheqflag = 1;

      cout << endl;
      cout << "What absolute tolerance do you want to use for";
      cout << " checking equilibrium";
      cout << endl;
      cout << "of the nitrogen cycle? ";

      cin >> ntol;

      rflog1 << endl;
      rflog1 << "What absolute tolerance do you want to use for";
      rflog1 << " checking equilibrium";
      rflog1 << endl;
      rflog1 << "of the nitrogen cycle?" << endl;
      rflog1 << "ntol = " << ntol << endl << endl;
    }

    if( 0 == equil )
    {

      cout << endl << endl;
      cout << "What year do you want to start collecting output data? ";

      cin >> startyr;

      rflog1 << endl << endl;
      rflog1 << "What year do you want to start collecting output data? ";
      rflog1 << "startyr = " << startyr << endl;

      cout << endl << endl;
      cout << "What year do you want to stop collecting output data? ";

      cin >> endyr;

      rflog1 << endl << endl;
      rflog1 << "What year do you want to stop collecting output data? ";
      rflog1 << "endyr = " << endyr << endl;

      cout << "How often (x years) should data be collected";
      cout << " after the initial year? ";

      cin >> diffyr;

      rflog1 << "How often (x years) should data be collected";
      rflog1 << " after the initial year? ";
      rflog1 << "diffyr = " << diffyr << endl;
    }
  #endif
  
};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void MITTEM44::massbal( double y[NUMEQ], 
                      double prevy[MAXSTATE] )
{
 
/***************** Carbon Cycle Balances **********************/

  if( y[I_INNPP] < y[I_NPP] ) { y[I_INNPP] = y[I_NPP]; }

  if( y[I_INGPP] < y[I_GPP] ) { y[I_INGPP] = y[I_GPP]; }

  if( y[I_GPR] != y[I_GPP] - y[I_NPP] )
  {
    y[I_GPR] = y[I_GPP] - y[I_NPP];
  }

  if( y[I_GPR] != y[I_RVMNT] + y[I_RVGRW] )
  {
    y[I_RVGRW] = y[I_GPR] - y[I_RVMNT];
  }


  if( y[I_VEGC] - prevy[I_VEGC]
      != y[I_NPP] - y[I_LTRC] )
  {
    y[I_LTRC] = y[I_NPP] - y[I_VEGC] + prevy[I_VEGC];
  }


  if( y[I_SOLC] - prevy[I_SOLC] != y[I_LTRC] + ag.getSLASHC() 
       - ag.getSCONVRTFLXC() - y[I_RH] )
  {
    y[I_RH] = y[I_LTRC] 
              + ag.getSLASHC() 
              - ag.getSCONVRTFLXC()
              - y[I_SOLC] 
              + prevy[I_SOLC];
  }
 

  /*********************Nitrogen Cycle Balances**********************/

  if( y[I_VNUP] < ZERO ) { y[I_VNUP] = ZERO; }

  if( y[I_INNUP] < y[I_VNUP] ) { y[I_INNUP] = y[I_VNUP]; }

  if( y[I_VSUP] < ZERO ) { y[I_VSUP] = ZERO; }

  if( y[I_VSUP] > y[I_VNUP] ) { y[I_VSUP] = y[I_VNUP]; }

  if( y[I_VLUP] != y[I_VNUP] - y[I_VSUP] )
  {
    y[I_VLUP] = y[I_VNUP] - y[I_VSUP];
  }


  // DWK modified the following conditions for checking the mass
  //   balance on STON on 0020401

  if( y[I_STON] - prevy[I_STON]
       != y[I_VLUP] + y[I_VNRSRB] - y[I_VNMBL] )
  {
    y[I_VNRSRB] = y[I_STON] 
                  - prevy[I_STON] 
                  + y[I_VNMBL]
                  - y[I_VLUP];
  }

   if( y[I_STRN] - prevy[I_STRN] != y[I_VSUP] - y[I_LTRN]
       - y[I_VNRSRB] + y[I_VNMBL] )
  {
    y[I_LTRN] = y[I_VSUP] 
                - y[I_VNRSRB] 
                + y[I_VNMBL] 
                - y[I_STRN] 
                + prevy[I_STRN];
  }

  if( y[I_SOLN] - prevy[I_SOLN] != y[I_LTRN] + ag.getSLASHN()
      - ag.getSCONVRTFLXN() - ag.getNSRETENT() - y[I_NMIN] )
  {
    y[I_NMIN] = y[I_LTRN] 
                + ag.getSLASHN() 
                - ag.getSCONVRTFLXN()
                - ag.getNSRETENT() 
                - y[I_SOLN]  
                + prevy[I_SOLN];
  }
  

  if( y[I_AGFRTN] < ZERO ) { y[I_AGFRTN] = ZERO; }

  if ( y[I_NINP] != y[I_AGFRTN] + ag.getNRETENT() + ag.getFIRENDEP() )
  {
    y[I_NINP] = y[I_AGFRTN] + ag.getNRETENT() + ag.getFIRENDEP();
  }

  if ( y[I_NINP] < ZERO ) { y[I_NINP] = ZERO; }

  if ( y[I_NLST] < ZERO ) { y[I_NLST] = ZERO; }


  if( y[I_AVLN] - prevy[I_AVLN]
      != y[I_NINP] + y[I_NMIN] - y[I_VNUP] - y[I_NLST] )
  {
    if( y[I_NINP] + y[I_NMIN] - y[I_VNUP] - y[I_NLST] 
          + prevy[I_AVLN] > ZERO )
    {	
      y[I_NLST] =  y[I_NINP] 
                   + y[I_NMIN] 
                   - y[I_VNUP] 
                   - y[I_AVLN]  
                   + prevy[I_AVLN];
    }
    else
    {	
      y[I_NINP] =  y[I_NLST] 
                   - y[I_NMIN] 
                   + y[I_VNUP] 
                   + y[I_AVLN]  
                   - prevy[I_AVLN];
    }   
  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

int MITTEM44::monthlyTransient( const int& pdyr,
                                const int& pdm,
                                const double& ptol )
{
  endeq = 0;
  initFlag = 1;

  // Reset to appropriate rooting depth and texture-dependent
  //   vegetation parameters for each cohort during each month
  //   (Note: microbe and soil parameters are the same between
  //   disturbed [ag.cmnt] and natural [veg.cmnt] land cover

  if( 0 == ag.state 
      || (ag.state > 0 && disturbflag > 0 && pdm < (disturbmonth-1)) )
  {
    soil.updateRootZ( veg.cmnt );

    veg.resetEcds( veg.cmnt, soil.getPSIPLUSC() );
  }
  else
  {
    soil.updateRootZ( ag.cmnt );

    veg.resetEcds( ag.cmnt, soil.getPSIPLUSC() );
  }

  // Reset texture-dependent parameters for microbes
  //   in current cohort

  microbe.resetEcds( veg.cmnt, soil.getPSIPLUSC() );


  // Update monthly carbon, nitrogen and water dynamics of
  //   terrestrial ecosystem
   
  stepmonth( pdyr, pdm, intflag, ptol );


  if( totyr == startyr ) { wrtyr = 0; }

  if( (CYCLE-1) == pdm )
  {
    if( totyr > startyr ) { ++wrtyr; }
  }

  return wrtyr;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITTEM44::natvegDynamics( const int& pdm, double pstate[] )
{
  int agstate = 0;
  int perennial = 1;
  int fertflag = 0;
  int tillflag = 0;                       
  
  // Assume no agricultural N fertilization (ag.fertn) and no
  //  nutrients resulting from agricultural conversion

  ag.fertn = ZERO;

  // Add fireNDEP to ecosystem if fire occurred within the past \
  //   FRI years.  Otherwise do not add any nitrogen to ecosyste
 
  if( firemnthcnt > 0 ) 
  { 
    soil.setNINPUT( ag.getNRETENT() + ag.getFIRENDEP() ); 
  }
  else { soil.setNINPUT( ZERO ); } 

  soil.setEET( soil.getINEET() );

  soil.setVSM( (soil.getMOIST() / (soil.getROOTZ() * 1000.0)) );

  soil.setKH2O( soil.getVSM(), moistlim );


  microbe.updateDynamics( veg.cmnt,
                          soil.getPCTFLDCAP(),
                          pstate[I_SOLC],
                          pstate[I_SOLN],
                          soil.getMOIST(),
                          soil.getVSM(),
                          pstate[I_AVLN],
                          moistlim,
                          tillflag,
                          ag.getTILLFACTOR( veg.cmnt ),
                          soil.getKH2O() );

  if( disturbflag > 1 && pdm == (disturbmonth-1) )
  {
    // Set monthly vegetation fluxes to zero
    	
    veg.resetMonthlyFluxes();  	
  } 
  else
  {
    veg.updateDynamics( veg.cmnt,
                        atms.getCO2(),
                        atms.getAOT40(),
                        soil.getNINPUT(),
                        atms.getPAR(),
                        atms.getPET(),
                        atms.getPRVPETMX(),
                        soil.getEET(),
                        soil.getPRVEETMX(),
                        pstate[I_VEGC],
                        pstate[I_STRN],
                        pstate[I_STON],
                        soil.getMOIST(),
                        pstate[I_AVLN],
                        moistlim,
                        nfeed,
                        o3flag,
                        agstate,
                        perennial,
                        fertflag,
                        soil.getKH2O(),
                        microbe.getNETNMIN(),
                        ag.fertn );
  }                    

  // Determine nitrogen losses from ecosystem

  if( 1 == avlnflag )
  {
    soil.updateNLosses( veg.cmnt,
                        (soil.getSURFRUN() + soil.getDRAINAGE()),
                        pstate[I_AVLN], 
                        soil.getMOIST() );

    if( soil.getNLOST() > pstate[I_AVLN] - veg.getNUPTAKE()
         + microbe.getNETNMIN() + soil.getNINPUT() )
    {
      soil.setNLOST( (pstate[I_AVLN] 
                     - veg.getNUPTAKE() 
                     + microbe.getNETNMIN()
                     + soil.getNINPUT()) );
    }
    if( soil.getNLOST() < ZERO )
    {
      soil.setNLOST( ZERO );
      
      microbe.setNETNMIN( (soil.getNLOST() 
                          + veg.getNUPTAKE() 
                          - soil.getNINPUT()
                          - pstate[I_AVLN]) );
    }                        
  }
  else
  {
    // Do not allow changes in available nitrogen in soil
    //   (used during calibration process)

    soil.setNLOST( (soil.getNINPUT()
                   + microbe.getNETNMIN()
                       - veg.getNUPTAKE()) );
  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void MITTEM44::pastureDynamics( const int& pdm, double pstate[] )
{
  int fertflag = 0;
  
  int perennial = 1;
  
  int tillflag = 0;                       
  
  // Initialize fertilizer input to zero
  
  ag.fertn = ZERO;
 
  
  // Add fireNDEP to ecosystem if fire occurred within the past \
  //   FRI years.  Otherwise do not add any nitrogen to ecosyste
 
  if( firemnthcnt > 0 ) 
  { 
    soil.setNINPUT( ag.getNRETENT() + ag.getFIRENDEP() ); 
  }
  else { soil.setNINPUT( ag.getNRETENT()  ); } 

  soil.setEET( soil.getINEET() );

  soil.setVSM( (soil.getMOIST() / (soil.getROOTZ() * 1000.0)) );

  soil.setKH2O( soil.getVSM(), moistlim );


  microbe.updateDynamics( veg.cmnt,
                          soil.getPCTFLDCAP(),
                          pstate[I_SOLC],
                          pstate[I_SOLN],
                          soil.getMOIST(),
                          soil.getVSM(),
                          pstate[I_AVLN],
                          moistlim,
                          tillflag,
                          ag.getTILLFACTOR( veg.cmnt ),
                          soil.getKH2O() );

  if( disturbflag > 1 && pdm == (disturbmonth-1) )
  {
    // Set monthly vegetation fluxes to zero
    	
    veg.resetMonthlyFluxes();  	
  } 
  else
  {
    veg.updateDynamics( ag.cmnt,
                        atms.getCO2(),
                        atms.getAOT40(),
                        soil.getNINPUT(),
                        atms.getPAR(),
                        atms.getPET(),
                        atms.getPRVPETMX(),
                        soil.getEET(),
                        soil.getPRVEETMX(),
                        pstate[I_VEGC],
                        pstate[I_STRN],
                        pstate[I_STON],
                        soil.getMOIST(),
                        pstate[I_AVLN],
                        moistlim,
                        nfeed,
                        o3flag,
                        ag.state,
                        perennial,
                        fertflag,
                        soil.getKH2O(),
                        microbe.getNETNMIN(),
                        ag.fertn );
  }                    

  // Update NINPUT to account for any fertilizer application
  
//  soil.setNINPUT( soil.getNINPUT() + ag.fertn );
  
  // Determine nitrogen losses from ecosystem

  if( 1 == avlnflag )
  {
    soil.updateNLosses( veg.cmnt,
                        (soil.getSURFRUN() + soil.getDRAINAGE()),
                        pstate[I_AVLN], 
                        soil.getMOIST() );

    if( soil.getNLOST() > pstate[I_AVLN] - veg.getNUPTAKE()
         + microbe.getNETNMIN() + soil.getNINPUT() )
    {
      soil.setNLOST( (pstate[I_AVLN] 
                     - veg.getNUPTAKE() 
                     + microbe.getNETNMIN()
                     + soil.getNINPUT()) );
    }
    if( soil.getNLOST() < ZERO )
    {
      soil.setNLOST( ZERO );
      
      microbe.setNETNMIN( (soil.getNLOST() 
                          + veg.getNUPTAKE() 
                          - soil.getNINPUT()
                          - pstate[I_AVLN]) );
    }                        
  }
  else
  {
    // Do not allow changes in available nitrogen in soil
    //   (used during calibration process)

    soil.setNLOST( (soil.getNINPUT()
                   + microbe.getNETNMIN()
                       - veg.getNUPTAKE()) );
  }

};


/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

#ifdef CALIBRATE_TEM

void MITTEM44::pcdisplayDT( const double& tottime,
                             const double& deltat )
{

  window( 1,15,39,15 );
  gotoxy( 1,1 );
  printf( "TIME = %10.8lf   DT = %10.8lf    ",
          tottime,
          deltat );

};

#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

#ifdef CALIBRATE_TEM

void MITTEM44::pcdisplayMonth( const int& dyr, const int& pdm )
{
  char limitC[2];

  double outflux1;
  double outflux2;
  double outflux3;
  double outflux4;
  double outflux5;
  double outflux6;
  double outflux7;
  double outflux8;

  if( 0 == topwind )
  {
    window( 1,1,80,24 );
    clrscr();
    switch( calwind )
    {

      case 1:  window( 1,1,43,1 );
               cout << "    TIME   AVAILW   RGRNDW  SNOWPCK SRGRNDW";
               break;

      case 2:  window( 1,1,48,1 );
               cout << " TIME    VEG C    STR N   SOIL C   SOIL N  AVALN LABILN";
               break;
    }

    topwind = 1;
  }

  // Assign values of optional variables to outflux? for later 
  //   screen display

  outflux1 = getOptionalEflx( sey[0] );

  outflux2 = getOptionalEflx( sey[1] );

  outflux3 = getOptionalEflx( sey[2] );

  outflux4 = getOptionalWflx( swy[0] );

  outflux5 = getOptionalWflx( swy[1] );

  outflux6 = getOptionalWflx( swy[2] );

  outflux7 = getOptionalWflx( swy[3] );

  outflux8 = getOptionalWflx( swy[4] );

  window( 1,2,80,13 );
  gotoxy( 1,1 );
  delline();
  gotoxy( 1,12 );

  // Display monthly values for selected C and N pools and fluxes

  switch( calwind )
  {
    case 1: printf( "%4d-%2d %9.2lf %8.2lf %8.2lf %7.2lf %6.1lf %6.1lf %6.1lf %6.1lf %7.2lf",
                    dyr,
                    (pdm+1),
                    soil.getAVLH2O(),
                    ZERO,
                    soil.getSNOWPACK(),
                    ZERO,
                    outflux4,
                    outflux5,
                    outflux6,
                    outflux7,
                    outflux8 );
            break;
            
    case 2: // Productivity is nitrogen limited

            if( (y[I_INNPP] > y[I_NPP])
                && (y[I_INNUP] == y[I_VNUP]) )
            {
              strcpy( limitC, "N" );
            }

            // Productivity is either carbon or nitrogen limited

            if( (y[I_INNPP] == y[I_NPP])
                && (y[I_INNUP] == y[I_VNUP]) )
            {
              strcpy( limitC, "E" );
            }

            // Productivity is carbon limited (climate)

            if( (y[I_INNPP] == y[I_NPP])
                && (y[I_INNUP] > y[I_VNUP]) )
            {
              strcpy( limitC, "C" );
            }

            // Productivity is limited by both carbon and nitrogen

            if( (y[I_INNPP] > y[I_NPP])
                && (y[I_INNUP] > y[I_VNUP]) )
            {
              strcpy( limitC, "B" );
            }

            // Unknown limits on productivity

            if( y[I_INNPP] < y[I_NPP]
                || y[I_INNUP] < y[I_VNUP] )
            {
              strcpy( limitC, "?" );
            }

            printf( "%3d-%2d %8.2lf %7.2lf %9.2lf %7.2lf %6.3lf %6.3lf %7.3lf %6.2lf %s %7.3lf",
                    dyr,
                    (pdm+1),
                    y[I_VEGC],
                    y[I_STRN],
                    y[I_SOLC],
                    y[I_SOLN],
                    y[I_AVLN],
                    y[I_STON],
                    outflux1,
                    outflux2,
                    limitC,
                    outflux3 );

            break;
  }

  window( 1,14,80,14 );
  gotoxy( 1,1 );
  delline();
//  printf("                                                                         ");

};

#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

#ifdef CALIBRATE_TEM

void MITTEM44::pcdisplayODEerr( const int& test, 
                                 double pstate[] )
{
  int errtest;

  window( 40,15,80,15 );
  gotoxy( 1,1 );

  if( test != ACCEPT )
  {
    errtest = test - 1;

    switch( errtest )
    {
      case   I_VEGC: printf( "   VEGC = %8.2lf  Error = %11.8lf  ",
                             pstate[I_VEGC],
                             error[I_VEGC] );
                     break;

      case   I_SOLC: printf( "  SOILC = %8.2lf  Error = %11.8lf  ",
                             pstate[I_SOLC],
                             error[I_SOLC] );
                     break;

      case   I_STRN: printf( "   STRN = %8.2lf  Error = %11.8lf  ",
                             pstate[I_STRN],
                             error[I_STRN] );
                     break;

      case   I_STON: printf( "LABILEN = %8.3lf  Error = %11.8lf  ",
                             pstate[I_STON],
                             error[I_STON] );
                     break;

      case   I_SOLN: printf( "  SOILN = %8.2lf  Error = %11.8lf  ",
                             pstate[I_SOLN],
                             error[I_SOLN] );
                     break;

      case   I_AVLN: printf( "   AVLN = %8.5lf  Error = %11.8lf  ",
                             pstate[I_AVLN],
                             error[I_AVLN] );
                     break;

      case  I_INGPP: printf( "INITGPP = %8.1lf  Error = %11.8lf  ",
                             pstate[I_INGPP],
                             error[I_INGPP] );
                     break;

      case    I_GPP: printf( "    GPP = %8.1lf  Error = %11.8lf  ",
                             pstate[I_GPP],
                             error[I_GPP] );
                     break;

      case  I_INNPP: printf( "INITNPP = %8.1lf  Error = %11.8lf  ",
                             pstate[I_INNPP],
                             error[I_INNPP] );
                     break;

      case    I_NPP: printf( "    NPP = %8.1lf  Error = %11.8lf  ",
                             pstate[I_NPP],
                             error[I_NPP] );
                     break;

      case    I_GPR: printf( "     RA = %8.1lf  Error = %11.8lf  ",
                             pstate[I_GPR],
                             error[I_GPR] );
                     break;

      case  I_RVMNT: printf( "     RM = %8.1lf  Error = %11.8lf  ",
                             pstate[I_RVMNT],
                             error[I_RVMNT] );
                     break;

      case  I_RVGRW: printf( "     RG = %8.1lf  Error = %11.8lf  ",
                             pstate[I_RVGRW],
                             error[I_RVGRW] );
                     break;

      case   I_LTRC: printf( "   LTRC = %8.1lf  Error = %11.8lf  ",
                             pstate[I_LTRC],
                             error[I_LTRC] );
                     break;

      case     I_RH: printf( "     RH = %8.1lf  Error = %11.8lf  ",
                             pstate[I_RH],
                             error[I_RH] );
                     break;

      case I_AGFRTN: printf( "AGFERTN = %8.3lf  Error = %11.8lf  ",
                              pstate[I_AGFRTN],
                              error[I_AGFRTN] );
                     break;

      case  I_INNUP: printf( "INUPTAK = %8.3lf  Error = %11.8lf  ",
                             pstate[I_INNUP],
                             error[I_INNUP] );
                     break;

      case   I_VNUP: printf( "NUPTAKE = %8.3lf  Error = %11.8lf  ",
                             pstate[I_VNUP],
                             error[I_VNUP] );
                     break;

      case   I_VSUP: printf( "SUPTAKE = %8.3lf  Error = %11.8lf  ",
                             pstate[I_VSUP],
                             error[I_VSUP] );
                     break;

      case   I_VLUP: printf( "LUPTAKE = %8.3lf  Error = %11.8lf  ",
                             pstate[I_VLUP],
                             error[I_VLUP] );
                     break;

      case  I_VNMBL: printf( " NMOBIL = %8.3lf  Error = %11.8lf  ",
                             pstate[I_VNMBL],
                             error[I_VNMBL] );
                     break;

      case   I_LTRN: printf( "   LTRN = %8.3lf  Error = %11.8lf  ",
                             pstate[I_LTRN],
                             error[I_LTRN] );
                     break;

      case   I_MNUP: printf( "MCRONUP = %8.3lf  Error = %11.8lf  ",
                             pstate[I_MNUP],
                             error[I_MNUP] );
                     break;

      case   I_NMIN: printf( "   NMIN = %8.3lf  Error = %11.8lf  ",
                             pstate[I_NMIN],
                             error[I_NMIN] );
                     break;

      case I_UNRMLF: printf( " UNLEAF = %8.3lf  Error = %11.8lf  ",
                             pstate[I_UNRMLF],
                             error[I_UNRMLF] );
                     break;

      case   I_LEAF: printf( "   LEAF = %8.3lf  Error = %11.8lf  ",
                             pstate[I_LEAF],
                             error[I_LEAF] );
                     break;
    }
  }

};

#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void MITTEM44::resetMonthlyELMNTFluxes( void )
{

  // Reset all monthly fluxes to zero

  atms.resetMonthlyFluxes();
  
  veg.resetMonthlyFluxes();
  
  microbe.resetMonthlyFluxes();
  
  soil.resetMonthlyFluxes();
  
  soil.resetMonthlyTraceGasFluxes();   // For NEM and MDM
    
  ag.resetMonthlyFluxes();
  
  // Carbon fluxes 

  nep = ZERO;

  nce = ZERO;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void MITTEM44::resetODEflux( void )
{
  int i;

  for ( i = MAXSTATE; i < NUMEQ; ++i ) { y[i] = ZERO; }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void MITTEM44::resetYrFluxes( void )
{

  atms.resetYrFluxes();
  
  veg.resetYrFluxes();
  
  microbe.resetYrFluxes();
  
  soil.resetYrFluxes();

  soil.resetYrTraceGasFluxes();   // For NEM and MDM
  
  ag.resetYrFluxes();
  
  yrtotalc = ZERO;

  // Annual carbon fluxes

  yrnep = ZERO;
  yrnce = ZERO;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void MITTEM44::rkf( const int& numeq, 
                    double pstate[], 
                    double& pdt,
                   const int& pdm )
{

  int i;
  double ptdt = ZERO;

  for( i = 0; i < numeq; ++i )
  {
    dum4[i] = dum5[i] = pstate[i];

    yprime[i] = rk45[i] = error[i] = ZERO;
  }

  ptdt = pdt * 0.25;

  delta( pdm, dum4, f11 );
  
  step( numeq, yprime, f11, yprime, a1 );
  
  step( numeq, rk45, f11, rk45, b1 );
  
  step( numeq, dum4, f11, ydum, ptdt );
  
  delta( pdm, ydum, f2 );
  
  for( i = 0; i < numeq; ++i ) 
  {
    f13[i] = a31*f11[i] + a32*f2[i];
  }
  
  step( numeq, dum4, f13, ydum, pdt );
  
  delta( pdm, ydum, f3 );
  
  step( numeq, yprime, f3, yprime, a3 );
  
  step( numeq, rk45, f3, rk45, b3 );
  
  for( i = 0; i < numeq; ++i )
  {
    f14[i] = a41*f11[i] + a42*f2[i] + a43*f3[i];
  }
  
  step( numeq, dum4, f14, ydum, pdt );
  
  delta( pdm, ydum, f4 );
  
  step( numeq, yprime, f4, yprime, a4 );
  
  step( numeq, rk45, f4, rk45, b4 );
  
  for( i = 0; i < numeq; ++i )
  {
    f15[i] = a51*f11[i] + a52*f2[i] + a53*f3[i] + a54*f4[i];
  }
  
  step( numeq, dum4, f15, ydum, pdt );
  
  delta( pdm, ydum, f5 );
  
  step( numeq, yprime, f5, yprime, a5 );
  
  step( numeq, rk45, f5, rk45, b5 );
  
  for( i = 0; i < numeq; ++i )
  {
    f16[i] = b61*f11[i] + b62*f2[i] + b63*f3[i] + b64*f4[i] 
             + b65*f5[i];
  }
  
  step( numeq, dum4, f16, ydum, pdt );
  
  delta( pdm, ydum, f6 );
  
  step( numeq, rk45, f6, rk45, b6 );
  
  step( numeq, dum4, yprime, dum4, pdt );
  
  step( numeq, dum5, rk45, dum5, pdt );
  
  for ( i = 0; i < numeq; ++i )
  {
    error[i] = fabs( dum4[i] - dum5[i] );
  }

};

/***************************************************************
 ***************************************************************/


/* *************************************************************
************************************************************** */

void MITTEM44::setELMNTecd( const int& pdcmnt,
                            const double& psiplusc )
{
  // Initialize TEM parameters dependent upon a grid cell's
  //   soil texture

  veg.resetEcds( pdcmnt, psiplusc );
  
  microbe.resetEcds( pdcmnt, psiplusc );
  	  
};

/* *************************************************************
************************************************************** */

/* *************************************************************
************************************************************** */

void MITTEM44::setEquilC2N( const int& pdcmnt, 
                            const double& co2 )
{

  atms.yrpet = 1.0;

  soil.yreet = 1.0;

  // Determine vegetation C/N parameter as a function of
  //   vegetation type, annual PET, and annual EET (annual EET
  //   initially equal to yrpet)

  veg.updateC2N( pdcmnt,
                 soil.yreet,
                 atms.yrpet,
                 co2,
                 atms.getINITCO2() );

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void MITTEM44::setEquilEvap( const double& nirr,
                             const double& tair,
                             const int& pdm )
{

  // Determine initial values for atms.prvpetmx, atms.prveetmx,
  //   and veg.topt

  atms.petjh( nirr, tair, pdm );

  if( 0 == pdm )
  {
    atms.setPRVPETMX( atms.getPET() );

    soil.setPRVEETMX( atms.getPET() );

    veg.setTOPT( tair );
  }
  else
  {
    if( atms.getPET() > atms.getPRVPETMX() )
    {
      atms.setPRVPETMX( atms.getPET() );

      soil.setPRVEETMX( atms.getPET() );

      veg.setTOPT( tair );
    }
  }
  

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void MITTEM44::setPrevState( void )
{
  for( int i = 0; i < MAXSTATE; ++i ) { prevy[i] = y[i]; }

};

/* *************************************************************
************************************************************* */


/***************************************************************
 ***************************************************************/

void MITTEM44::step( const int& numeq, 
                   double pstate[], 
                   double pdstate[], 
                   double ptstate[],
	           double& pdt )
{

  for( int i = 0; i < numeq; ++i )
  {
    ptstate[i] = pstate[i] + (pdt * pdstate[i]);
  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

int MITTEM44::stepmonth( const int& pdyr,
                         const int& pdm,
                         int& intflag,
                         const double& ptol )
{
  // For NEM and MDM *******************************************
  
  const double gCH4togC = 12.0 / 16.0;  

  int cultIntensity;

  int dday;

  int dlyr; 

  int iwet;

// *************************************************************

  int mintflag;

  // For NEM and MDM *******************************************

  double pctinundated;
 
  double soilLayerSH2O[10];

  //double soilLayerTSOIL[10];
  
  double totalSoilC;

  // Water table depth in millimeters
  double watertableZ;

// *************************************************************
 
  // Reset all monthly fluxes to zero
  
  resetMonthlyELMNTFluxes();


  // Reset ODE state variables representing monthly
  //   fluxes to zero
  
  resetODEflux();

  // If 12 months have passed since disturbance, reset
  //   all immediate fluxes associated with disturbance 
  //   (i.e., slash, convrtflx, nretent) to zero
  
  if( CYCLE == distmnthcnt )
  {
    distmnthcnt = 0;

    ag.resetMonthlyDisturbFluxes();
  } 

  // Count the number of months since disturbance
  
  if( distmnthcnt > 0 ) { ++distmnthcnt; }


  // If (FRI times 12 months have passed since a fire
  //   disturbance, reset ag.firendep to zero
  
  if( (ag.getFRI() * 12) == firemnthcnt )
  {
    firemnthcnt = 0;

    ag.setFIRENDEP( ZERO );
  }

  // Count the number of months since a fire disturbance
  
  if( firemnthcnt > 0 ) { ++firemnthcnt; };
  
  if( 0 == ag.state && 2 == distmnthcnt )
  {
     // Establish natural vegetation the month
     //   following disturbance if cohort is not
     //   in agriculture

      prevy[I_VEGC] = y[I_VEGC] = ag.getNATSEEDC();

      prevy[I_STRN] = y[I_STRN] = ag.getNATSEEDSTRN();

      prevy[I_STON] = y[I_STON] = ag.getNATSEEDSTON();  
  }

  if( 0 == pdm )
  {
    if( 0 == pdyr )
    {
      microbe.setKD( microbe.getKDC() );

      ag.setKD( microbe.getKD() );

      ag.setNATSOIL( y[I_SOLC] );
    }
    else
    {
      if( 0 == ag.state && 0 == ag.prvstate )
      {
        microbe.setKD( microbe.yrkd( nfeed,
                                     veg.yrltrc,
                                     veg.yrltrn,
                                     veg.cmnt ) );

        ag.setKD( microbe.getKD() );

        ag.setNATSOIL( y[I_SOLC] );
      }
      else
      {
        if( y[I_SOLC] < ag.getNATSOIL() )
        {
          microbe.setKD( (ag.getKD() 
                         * y[I_SOLC] 
                         / ag.getNATSOIL()) );
        }
        else { microbe.setKD( ag.getKD() ); }
      }
    }


    // Reset annual fluxes to zero

    resetYrFluxes();
  }

  // Implement disturbance effects

  if( disturbflag > 0 && pdm == (disturbmonth-1) )
  {
    distmnthcnt = 1;

    // Save proportion of vegetation for "seed" for regrowth
    // after abandonment

    ag.setNATSEEDC( (y[I_VEGC] * ag.getVRESPAR()) );

    ag.setNATSEEDSTRN( (y[I_STRN] * ag.getVRESPAR()) );

    ag.setNATSEEDSTON( (y[I_STON] * ag.getVRESPAR()) );


    // Create PROD10 and PROD100 from conversion to agriculture

    ag.createWoodProducts( pdyr,
                           y[I_VEGC],
                           y[I_STRN],
                           y[I_STON] );

      
    // Determine carbon and nitrogen lost during the first year
    //   due to conversion of natural lands to agriculture and/or
    //   use of fuelwood
      
    ag.conversion( veg.cmnt,
                   y[I_VEGC],
                   y[I_STRN],
                   y[I_STON],
                   y[I_SOLC],
                   y[I_SOLN] );


    // Establish crops

    prevy[I_VEGC] = y[I_VEGC] = ZERO;

    prevy[I_STRN] = y[I_STRN] = ZERO;

    prevy[I_STON] = y[I_STON] = ZERO;

    if( 3 == disturbflag )   // fire disturbance = 3
    {
      ag.setFireNDEP();

      firemnthcnt = 1;
    }
      
    if( 1 == disturbflag )
    {
      if( 1 == ag.state )
      {
        // Update rooting depth to be appropriate to crops

        soil.updateRootZ( ag.cmnt );


        // Update soil texture-dependent vegetation parameters
        //   for crops

        veg.resetEcds( ag.cmnt, soil.getPSIPLUSC() );

        // Update other adaptive parameters

        atms.yrpet = 1.0;

        soil.yreet = 1.0;
      
        veg.updateC2N( ag.cmnt,
                       soil.yreet,
                       atms.yrpet,
                       atms.getPREVCO2(),
                       atms.getINITCO2() );

        veg.setPRVLEAFMX( ag.getCROPPRVLEAFMX() );
      
        veg.setTOPT( ag.getCROPTOPT() );
      
        atms.setPRVPETMX( ag.getCROPPRVPETMX() );
      
        soil.setPRVEETMX( ag.getCROPPRVEETMX() );
      
        ag.setPRVCROPNPP( ZERO );
      }
      
      if( 2 == ag.state )
      {
        // Update rooting depth to be appropriate to pasture

        soil.updateRootZ( ag.cmnt );


        // Update soil texture-dependent vegetation parameters
        //   for pasture

        veg.resetEcds( ag.cmnt, soil.getPSIPLUSC() );

        // Update other adaptive parameters

        atms.yrpet = 1.0;
 
        soil.yreet = 1.0;
      
        veg.updateC2N( ag.cmnt,
                       soil.yreet,
                       atms.yrpet,
                       atms.getPREVCO2(),
                       atms.getINITCO2() );

        veg.setPRVLEAFMX( ag.getCROPPRVLEAFMX() );
     
        veg.setTOPT( ag.getCROPTOPT() );
      
        atms.setPRVPETMX( ag.getCROPPRVPETMX() );
      
        soil.setPRVEETMX( ag.getCROPPRVEETMX() );
      
        ag.setPRVCROPNPP( ZERO );
      }      

      if( 3 == ag.state )
      {
        // Update rooting depth to be appropriate to urban areas

        soil.updateRootZ( ag.cmnt );


        // Update soil texture-dependent vegetation parameters
        //   for urban areas

        veg.resetEcds( ag.cmnt, soil.getPSIPLUSC() );

        // Update other adaptive parameters

        atms.yrpet = 1.0;
 
        soil.yreet = 1.0;
      
        veg.updateC2N( ag.cmnt,
                       soil.yreet,
                       atms.yrpet,
                       atms.getPREVCO2(),
                       atms.getINITCO2() );

        veg.setPRVLEAFMX( ag.getCROPPRVLEAFMX() );
     
        veg.setTOPT( ag.getCROPTOPT() );
      
        atms.setPRVPETMX( ag.getCROPPRVPETMX() );
      
        soil.setPRVEETMX( ag.getCROPPRVEETMX() );
      
        ag.setPRVCROPNPP( ZERO );
      }      
    }
  }     
//  else { ag.setNoWoodProducts( pdyr ); }

  // Revert to natural vegetation after cropland abandonment

  if( 0 == ag.state && ag.prvstate > 0 )
  {
    // Update rooting depth to be appropriate to natural vegetation

    soil.updateRootZ( veg.cmnt );

    // Establish natural vegetation

    prevy[I_VEGC] = y[I_VEGC] = ag.getNATSEEDC();

    prevy[I_STRN] = y[I_STRN] = ag.getNATSEEDSTRN();

    prevy[I_STON] = y[I_STON] = ag.getNATSEEDSTON();

    // Update soil texture-dependent vegetation parameters
    //   for natural vegetation

    veg.resetEcds( veg.cmnt, soil.getPSIPLUSC() );


    // Update other adaptive parameters

    atms.yrpet = ag.getNATYRPET();

    soil.yreet = ag.getNATYREET();

    veg.updateC2N( veg.cmnt,
                   soil.yreet,
                   atms.yrpet,
                   atms.getPREVCO2(),
                   atms.getINITCO2() );

    veg.setPRVLEAFMX( ag.getNATPRVLEAFMX() );

    veg.setTOPT( ag.getNATTOPT() );

    atms.setPRVPETMX( ag.getNATPRVPETMX() );

    soil.setPRVEETMX( ag.getNATPRVEETMX() );

    ag.setPRVCROPNPP( ZERO );
  }


  // set ozone efect from previous month for initial month

  if( 0 == pdyr && 0 == pdm ) { veg.setFPREVOZONE( 1.0 ); }


  // Get environmental conditions for month "dm"

  getenviron( pdm );


  // Determine effect of air temperature on GPP (i.e. temp)
  //   and GPR (i.e. respq10)
  
  if( 0 == ag.state )
  {
    veg.setTEMP( veg.cmnt, atms.getTAIR() );

    veg.setRESPQ10( veg.cmnt, atms.getTAIR() );
//    veg.setRESPQ10LaRS( veg.cmnt, atms.getTAIR() );
  }
  else
  {
    veg.setTEMP( ag.cmnt, atms.getTAIR() );

    veg.setRESPQ10( ag.cmnt, atms.getTAIR() );
//    veg.setRESPQ10LaRS( ag.cmnt, atms.getTAIR() );
  }

  // Determine effect of temperature on decomposition 
  //   (i.e. dq10)
 
  //cout << "SWE in MITTEM = " << soil.getSNOWPACK() << endl;

  microbe.setDQ10( veg.cmnt, 
                   atms.getTAIR(), 
                   soil.getSNOWPACK() );

//  microbe.setDQ10LaRS( veg.cmnt, 
//                       atms.getTAIR(), 
//                       soil.getSNOWPACK() );

  // Update growing degree days (GDD) 

  if( atms.getTAIR() >= GDDMIN )
  {
    ag.setGROWDD( ag.getGROWDD() 
                  + ((atms.getTAIR() - GDDMIN) 
                  * atms.getNDAYS( pdm )) );
  }
  else
  {
    // If "cold snap" (i.e. TAIR < GDDMIN) hits after crops 
    //   begin to grow, crops are assumed to die and resulting
    //   detritus is added to soil organic carbon and nitrogen
       	
    if( 1 == ag.state 
        && ag.getGROWDD() > GDDSEED )
    {
      ag.frostDamage( y[I_VEGC], y[I_STRN], y[I_STON] );

      y[I_VEGC] = ZERO;

      prevy[I_VEGC] = ZERO;

      y[I_STRN] = ZERO;

      prevy[I_STRN] = ZERO;

      y[I_STON] = ZERO;

      prevy[I_STON] = ZERO;

      y[I_SOLC] += ag.getSTUBBLEC();

      prevy[I_SOLC] = y[I_SOLC];

      y[I_SOLN] += ag.getSTUBBLEN();

      prevy[I_SOLN] = y[I_SOLN];
    }

    ag.setGROWDD( ZERO );
  }

//  if ( 191 == pdyr ) { dbugflg = 1; }


  // Run TEM for a monthly time step

  mintflag = adapt( NUMEQ, y, ptol, pdm );


  if( 1 == mintflag ) { intflag = 1; }

  if ( blackhol != 0 )
  {
    if( 0 == initFlag ) { qualcon[0] = 10; }
    else { qualcon[pdyr] += 10; }
  }


  // Check mass balance

  massbal( y, prevy );

  // Determine veg.plant.nitrogen
  
  veg.setVEGN( (y[I_STRN] + y[I_STON]) );
  

  // Determine Net Ecosystem Production (nep)

  nep = y[I_NPP] - y[I_RH];


  // Determine fluxes from crop residues

  ag.updateCropResidueFluxes();


  // Determine fluxes from decay of products

  ag.decayProducts();


  // Harvest crops after a specific number of growing degree 
  //   days; reset growing degree days to zero if crops were
  // harvested this month

  if( 1 == ag.state && ag.getGROWDD() >= GDDHARVST )
  {
    ag.harvest( pdm, y[I_VEGC], y[I_STRN], y[I_STON] );

    y[I_VEGC] = ZERO;

    y[I_STRN] = ZERO;

    y[I_STON] = ZERO;

    y[I_SOLC] += ag.getSTUBBLEC();

    y[I_SOLN] += ag.getSTUBBLEN();

    ag.setGROWDD( ZERO );
  }
  else { ag.setNoCrops( pdm ); }


  // Determine crop residue
  
  ag.updateCropResidue();

  
  // Determine standing stock of products

  ag.updateProducts();


  // Determine total amount of products formed
  //   this month
  
  ag.updateTotalProductFormation();

  
  // Graze veg biomass every month if cohort is in pasture

  if( 2 == ag.state )
  {
    ag.grazing(  y[I_VEGC], y[I_STRN] );

    y[I_VEGC] -= ag.getFORAGEC();

    y[I_STRN] -= ag.getFORAGEN();

    y[I_SOLC] += ag.getMANUREC();

    y[I_SOLN] += ag.getMANUREN();

    y[I_AVLN] += ag.getURINE();   
  }
  else { ag.setNoGrazing(); }

  // Determine CFLUX from ecosystem from NEP plus
  //   fluxes from burning associated with agricultural
  //   conversion or 
  
  ag.setCFLUX( (nep 
                - ag.getCONVRTFLXC()
                - ag.getCROPRESIDUEFLXC() 
                - ag.getANIMALRESP()) );


  // Determine Net Carbon Exchange (NCE) with atmosphere
  //   (CFLUX plus decay of agricultural and wood products)
  
  nce = ag.getCFLUX() - ag.getTOTPRODDECAYC();


  // Determine carbon storage in ecosystem
  
  ag.setTOTEC( (y[I_VEGC] 
               + y[I_SOLC]) );


  // Determine total carbon in ecosystems plus
  //   agricultural and wood products
                 
  totalc = ag.getTOTEC() + ag.getTOTPRODC(); 
 
  // Update total loss of nitrogen from ecosystem for flux 
  //   associated with crop residue
    
  soil.setNLOST( (soil.getNLOST()
                  + ag.getCONVRTFLXN()
                  + ag.getCROPRESIDUEFLXN()) );


  if( 0 == initFlag )
  {
    if( 0 == pdm )
    {
      // Determine total organic carbon storage for the 
      //   top 1 m of the soil profile based on Schlesinger (1997)

      switch( veg.getCURRENTVEG() )
      {
        case 1: totalSoilC = 11800.0; break;
        case 2: totalSoilC = 14900.0; break;
        case 3: totalSoilC = 14900.0; break;
        case 4: totalSoilC = 10400.0; break;
        case 5: totalSoilC = 11800.0; break;
        case 6: totalSoilC =  6900.0; break;
        case 7: totalSoilC = 11800.0; break;
        case 8: totalSoilC = 14900.0; break;
        case 9: totalSoilC =  6900.0; break;
        case 10: totalSoilC =  5600.0; break;
        case 11: totalSoilC = 21600.0; break;
        case 12: totalSoilC = 21600.0; break;
        case 13: totalSoilC = 19200.0; break;
        case 14: totalSoilC =  3700.0; break;
        case 15: totalSoilC = 12700.0; break;
        case 16: totalSoilC = 12700.0; break;
        case 17: totalSoilC = 68600.0; break;
        case 18: totalSoilC = 68600.0; break;
        case 19: totalSoilC = 68600.0; break;
        case 20: totalSoilC = 68600.0; break;
        case 21: totalSoilC = 68600.0; break;
        case 22: totalSoilC = 68600.0; break;
        case 23: totalSoilC = 68600.0; break;
        case 24: totalSoilC = 68600.0; break;
        case 25: totalSoilC = 68600.0; break;
        case 26: totalSoilC = 10400.0; break;
        case 27: totalSoilC =  3700.0; break;
        case 28: totalSoilC = 11800.0; break;
        case 29: totalSoilC = 19200.0; break;
        case 32: totalSoilC = 12700.0; break;
        case 33: totalSoilC = 12700.0; break;
        default: totalSoilC = ZERO;
      }

      microbe.nem.setNONSOLC( (totalSoilC - y[I_SOLC]) );

      if( microbe.nem.getNONSOLC() < ZERO )
      { 
        microbe.nem.setNONSOLC( ZERO );
      }
      //cout << "NONSOLC in mittem = " << microbe.nem.getNONSOLC() << endl;
    }
  }


  if( 1 == microbe.nem.startflag
      && 1 == ch4flag )
  {
    // Determine Methane Flux

    

    // Use MDM to estimate CH4 fluxes for arctic or boreal 
    //   ecosystems, otherwise, use NEM
     
    if( 1 == mdmflag
        && (2 == veg.getCURRENTVEG() 
        || 3 == veg.getCURRENTVEG()
        || 8 == veg.getCURRENTVEG()
        || 11 == veg.getCURRENTVEG()
        || 12 == veg.getCURRENTVEG()
        || 21 == veg.getCURRENTVEG()
        || 22 == veg.getCURRENTVEG()) )
    {
      switch( veg.getCURRENTVEG() )
      {
        case 19: iwet = 1;
                 pctinundated = 100.000;
                 break;

        case 20: iwet = 1;
                 pctinundated = 100.000;
                 break;

        case 21: iwet = 1;
                 pctinundated = 100.000;
                 break;

        case 22: iwet = 1;
                 pctinundated = 100.000;
                 break;
      
        default: iwet = 0;
                 pctinundated = ZERO;
      }

      if( 15 == veg.getCURRENTVEG() 
          || 16 == veg.getCURRENTVEG() )
      {
        cultIntensity = 5;
      }
      else { cultIntensity = 1; } 
      
      soil.setCH4EMISS( ZERO );

      soil.setCH4CONSUMP( ZERO );

      soil.setCH4FLUX( ZERO );
      
      for( dday = 0; dday < (int) atms.ndays[pdm]; ++dday )
      {                 
         microbe.mdm.soil.setTSOIL( ((microbe.nem.dayTemp[0][dday] * (1.75/20.0))
                                    + (microbe.nem.dayTemp[1][dday] * (2.75/20.0))
                                    + (microbe.nem.dayTemp[2][dday] * (4.5/20.0))
                                    + (microbe.nem.dayTemp[3][dday] * (7.5/20.0))
                                    + (microbe.nem.dayTemp[4][dday] * (3.5/20.0))) );

        for( dlyr = 0; dlyr < (CLMNLAYERS-1); ++dlyr )
        {
//          soilLayerTSOIL[dlyr] = microbe.nem.dayTemp[dlyr][dday];

          soilLayerSH2O[dlyr] = microbe.nem.dayMoist[dlyr][dday];
        }

        for( dlyr = 0; dlyr < MXMDMNLAY; ++dlyr )
        {
          microbe.mdm.soil.setINTERSOILT( microbe.mdm.soil.getTSOIL(),
                                          dlyr );

//          microbe.mdm.soil.setINTERSOILT( microbe.mdm.soil.interpolateLayers( soilLayerTSOIL, 
//                                                                             soil.layerThick,
//                                                                             (CLMNLAYERS-1),
//                                                                             (dlyr+1) ), 
//                                          dlyr );

          microbe.mdm.soil.setINTERSH2O( microbe.mdm.soil.interpolateLayers( soilLayerSH2O, 
                                                                             soil.layerThick,
                                                                             (CLMNLAYERS-1),
                                                                             (dlyr+1) ),
                                         dlyr );
        }

        
        if( 0 == iwet ){ watertableZ = 500.0; }
        else { watertableZ = 0.00; }

        soil.setPH( 7.5 );

        microbe.mdm.stepday( veg.cmnt,
                             cultIntensity,
                             iwet,
                             (pctinundated * 100.0),
                             watertableZ,
                             soil.getPCTSAND(),
                             soil.getPCTSILT(),
                             soil.getPCTCLAY(),
                             soil.getPCTFLDCAP(),
                             soil.getROOTZ(),
                             soil.getPH(),
                             veg.getNPP(),
                             microbe.mdm.soil.getTSOIL() );

        // Aggregate daily CH4 fluxes from MDM into monthly 
        //   fluxes
        
        soil.setCH4FLUX( (-1.0
                         * 0.001 
                         * (soil.getCH4FLUX() 
                         + (microbe.mdm.soil.getCH4TOT() 
                         * gCH4togC))) ); 
      }
    }
    else
    {   
      soil.setCH4EMISS( microbe.nem.setMonthlyCH4emissions( veg.getCURRENTVEG(),
                                                            pdm,
                                                            atms.getTAIR(),
                                                            atms.getPREC(),
                                                            soil.getEET(),
                                                            lat ) );

      soil.setCH4EMISS( soil.getCH4EMISS() * gCH4togC );

       soil.setCH4FLUX( soil.getCH4EMISS() );
    }
  } 
  else 
  { 
    soil.setCH4EMISS( ZERO ); 

    soil.setCH4FLUX( soil.getCH4EMISS() );
  }
  

  if( 1 == n2oflag
      && (1 == endeq || 1 == initFlag) 
      && y[I_SOLC] > ZERO
      && veg.getCURRENTVEG() > 0 
      && veg.getCURRENTVEG() < 30 )
  {

    if( 1 == microbe.nem.startflag )
    {
      totalSoilC = microbe.nem.getNONSOLC() + y[I_SOLC];

      //cout << " totalSoilC = " << totalSoilC << endl;
      //cout << "Present Month: " << pdm << endl;
      //cout << " NonSOLC = " << microbe.nem.getNONSOLC() << endl;
      //cout << " yISOLC = " << y[I_SOLC] << endl;

      microbe.nem.setPH( 7.00 );

      //cout << " TOPPOR: " << endl;
      //cout << microbe.nem.getTOPPOR() << endl;
      //cout << " TOPDENS: " << endl;
      //cout << microbe.nem.getTOPDENS() << endl;
      //cout << " TOPKSAT: " << endl;
      //cout << microbe.nem.getTOPKSAT() << endl;
      //cout << " dayTemp: " << endl;
      //cout << microbe.nem.dayTemp << endl;
      //cout << " dayMoist: " << endl;
      //cout << microbe.nem.dayMoist << endl;
      //cout << " HourMoist: " << endl;
      //cout << microbe.nem.hourMoist << endl;
      //cout << " dayTair: " << endl;
      //cout << atms.dayTair << endl;
      //cout << " RainDur: " << endl;
      //cout << atms.rainDuration << endl;
      //cout << " RainInt: " << endl;
      //cout << atms.rainIntensity << endl;

      microbe.nem.stepmonth( pdm,
                             veg.getCURRENTVEG(),
                             veg.cmnt,
                             soil.getPCTCLAY(),
                             microbe.nem.getTOPPOR(),
                             microbe.nem.getTOPDENS(),
                             microbe.nem.getPH(),
                             microbe.nem.getTOPKSAT(),
                             totalSoilC,
                             microbe.nem.dayTemp,
                             microbe.nem.dayMoist,
                             microbe.nem.hourMoist,
                             atms.dayTair,
                             atms.rainDuration,
                             atms.rainIntensity );

      soil.setCO2NFLUX( (microbe.nem.getECO2N() * 0.1) );
      
      soil.setCO2DNFLUX( (microbe.nem.getECO2DN() * 0.1) );

      soil.setN2ONFLUX( (microbe.nem.getEN2ON() * 0.0001) );
      
      soil.setN2ODNFLUX( (microbe.nem.getEN2ODN() * 0.0001) );
        
      soil.setN2OFLUX( ((microbe.nem.getEN2ON() 
                         + microbe.nem.getEN2ODN()) 
                         * 0.0001) );

      soil.setN2FLUX( (microbe.nem.getEN2DN() * 0.0001) );

      soil.setCH4CONSUMP( (microbe.nem.getECH4N() * 0.1) );
    }
    else 
    {
      soil.setCO2NFLUX( ZERO );
      
      soil.setCO2DNFLUX( ZERO );
 
      soil.setN2ONFLUX( ZERO );
      
      soil.setN2ODNFLUX( ZERO );
      
      soil.setN2OFLUX( ZERO );

      soil.setN2FLUX( ZERO );

      soil.setCH4CONSUMP( ZERO );
    }
  }  
  else 
  {
    soil.setCO2NFLUX( ZERO );
      
    soil.setCO2DNFLUX( ZERO );
 
    soil.setN2ONFLUX( ZERO );
      
    soil.setN2ODNFLUX( ZERO );
      
    soil.setN2OFLUX( ZERO );

    soil.setN2FLUX( ZERO );

    soil.setCH4CONSUMP( ZERO );
  }
  
  y[I_SOLC] -= soil.getCH4FLUX();


  // Update ANNUAL carbon, nitrogen and water pools and fluxes 
  //   from integrator results

  updateYearSummary();


  #ifdef CALIBRATE_TEM
    // Display monthly results to DOS screen
    pcdisplayMonth( pdyr, pdm );
  #endif


  if( 1 == ag.state && ag.getPRVCROPNPP() < y[I_NPP] )
  {
     ag.setPRVCROPNPP( y[I_NPP] );
  }
  else { ag.setPRVCROPNPP( ZERO ); }

  if( 1 == ag.state && ag.getPRVCROPNPP() > y[I_NPP] )
  {
    y[I_UNRMLF] = ZERO;
  }

  // Reset growing degree days to zero if crops were
  // harvested this month

  if( 1 == ag.state && ag.getGROWDD() >= GDDHARVST )
  {
    ag.setGROWDD( ZERO );
  }

  if( atms.getTAIR() < GDDMIN ) { ag.setGROWDD( ZERO ); }


  if( pdyr > 0 || pdm > 0 )
  {
    veg.setFPREVOZONE( (1.0 - 0.5 * (1.0 - y[I_FOZONE])) );
  }

  if( y[I_NPP] <= ZERO ) { veg.setFPREVOZONE( 1.0 ); }


  // Update atms.prevco2 for next month

  atms.setPREVCO2( atms.getCO2() );

  
  // Update atms.prevtair and atms.prev2tair for next month

  atms.setPREV2TAIR( atms.getPREVTAIR() );

  atms.setPREVTAIR( atms.getTAIR() );


  // Update previous unnormalized relative leaf area for
  //   next month
  
  veg.setPREVUNRMLEAF( y[I_UNRMLF] );
  
  // Update ag.prevPROD1, ag.prevPROD10 and ag.prevPROD100
  // for next month

  ag.setPREVPROD1C( ag.getPROD1C() );
  ag.setPREVPROD1N( ag.getPROD1N() );

  ag.setPREVPROD10C( ag.getPROD10C() );
  ag.setPREVPROD10N( ag.getPROD10N() );

  ag.setPREVPROD100C( ag.getPROD100C() );
  ag.setPREVPROD100N( ag.getPROD100N() );

  // Update ag.prevCropResidue for next month

  ag.setPREVCROPRESIDUEC( ag.getCROPRESIDUEC() );
  ag.setPREVCROPRESIDUEN( ag.getCROPRESIDUEN() );

  //  Update maximum EET, maximum PET, GPP optimum temperature
  //    (veg.topt), and maximum leaf cover (veg.prvleafmx) for
  //    the current year

  if( 0 == pdm )
  {
    soil.setEETMX( soil.getEET() );

    atms.setPETMX( atms.getPET() );

    veg.setNEWTOPT( atms.getTAIR() );

    veg.setNEWLEAFMX( y[I_UNRMLF] );
  }
  else
  {
    if( soil.getEET() > soil.getEETMX() )
    {
      soil.setEETMX( soil.getEET() );
    }

    if( atms.getPET() > atms.getPETMX() )
    {
      atms.setPETMX( atms.getPET() );
    }

    if( 0 == ag.state )
    {
      veg.resetNEWTOPT( veg.cmnt, 
                        atms.getTAIR(), 
                        y[I_UNRMLF] );
    }
    else 
    {
      veg.resetNEWTOPT( ag.cmnt, 
                        atms.getTAIR(), 
                        y[I_UNRMLF] );
    }    
  }


  // Save state of all the ODE state variables 
  //   representing pools to allow checking
  //   of mass balance
  
  setPrevState();


  // Update annual parameters for next year
  
  if( (CYCLE-1) == pdm )
  {
    soil.setPRVEETMX( soil.getEETMX() );

    atms.setPRVPETMX( atms.getPETMX() );

    veg.setTOPT( veg.getNEWTOPT() );

    veg.setPRVLEAFMX( veg.getNEWLEAFMX() );

    // Update optimum temperature parameters for GPP

    if( 0 == ag.state )
    {
      veg.boundTOPT( veg.cmnt );
      
    // Update adaptive parameters

      ag.setNATYRPET( atms.yrpet );

      ag.setNATYREET( soil.yreet );

      ag.setNATPRVPETMX( atms.getPRVPETMX() );

      ag.setNATPRVEETMX( soil.getPRVEETMX() );

      ag.setNATTOPT( veg.getTOPT() );

      ag.setNATPRVLEAFMX( veg.getPRVLEAFMX() );

      // Determine vegetation C/N parameter as a function
      // of vegetation type, annual PET, annual EET,
      // CO2 concentration

      veg.updateC2N( veg.cmnt,
                     soil.yreet,
                     atms.yrpet,
                     atms.getPREVCO2(),
                     atms.getINITCO2() );


    }
    else
    {
      veg.boundTOPT( ag.cmnt );

      ag.setCROPPRVPETMX( atms.getPRVPETMX() );

      ag.setCROPPRVEETMX( soil.getPRVEETMX() );

      ag.setCROPTOPT( veg.getTOPT() );

      ag.setCROPPRVLEAFMX( veg.getPRVLEAFMX() );

     // Determine vegetation C/N parameter as a function of
     //   vegetation type, annual PET, annual EET, CO2
     //   concentration

      veg.updateC2N( ag.cmnt,
                     soil.yreet,
                     atms.yrpet,
                     atms.getPREVCO2(),
                     atms.getINITCO2() );

    }

    // Update next year ag.prvstate with current year ag.state

    ag.prvstate = ag.state;

    veg.yrcarbon  /= 12.0;

    soil.yrorgc /= 12.0;

    yrtotalc /= 12.0;

    veg.yrnitrogen  /= 12.0;

    veg.yrstructn /= 12.0;

    if( veg.yrstructn != ZERO )
    {
      veg.yrc2n  = veg.yrcarbon / veg.yrstructn;
    }

    veg.yrstoren /= 12.0;

    soil.yrorgn /= 12.0;

    if( soil.yrorgn != ZERO )
    {
      soil.yrc2n = soil.yrorgc / soil.yrorgn;
    }

    soil.yravln  /= 12.0;

    soil.yravlh2o /= 12.0;

    soil.yrsmoist /= 12.0;

    soil.yrvsm /= 12.0;

    soil.yrpctp /= 12.0;

    soil.yrsnowpack /= 12.0;

    soil.yrrgrndh2o /= 12.0;

    soil.yrsgrndh2o /= 12.0;

    veg.yrunleaf /= 12.0;

    veg.yrleaf /= 12.0;

    veg.yrlai /= 12.0;

    veg.yrfpc /= 12.0;

    if( 1 == baseline )
    {
      soil.yrnin = ZERO;

      soil.yrnlost = ZERO;
      
      if( y[I_SOLC]/microbe.getCNSOIL( veg.cmnt ) >= y[I_SOLN] )
      {
        soil.yrnin = (y[I_SOLC] 
                     / microbe.getCNSOIL( veg.cmnt )) 
                     - y[I_SOLN];
      }
      else
      {
        soil.yrnlost = y[I_SOLN] 
                       - (y[I_SOLC]
                       / microbe.getCNSOIL( veg.cmnt ));
      }

      y[I_SOLN] = y[I_SOLC] / microbe.getCNSOIL( veg.cmnt );
    }
                       
    if( endeq > 0 ) 
    { 
      ++endeq; 
    }
  }

  return endeq;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int MITTEM44::testEquilibrium( void )
{

  if( 0 == nfeed && 0 == rheqflag
      && (ctol >= fabs( veg.yrnpp - veg.yrltrc )) )
  {
    return 1;
  }

  if( 0 == nfeed && 1 == rheqflag
      && (ctol >= fabs( yrnep ))
      && (ctol >= fabs( veg.yrnpp - veg.yrltrc ))
      && (ctol >= fabs( veg.yrltrc - microbe.yrrh )) )
  {
    return 1;
  }

  if( 1 == nfeed && 1 == rheqflag
      && (ntol >= fabs( soil.yrnin - soil.yrnlost ))
      && (ntol >= fabs( veg.yrnup - veg.yrltrn ))                         
      && (ntol >= fabs( veg.yrnup - microbe.yrnmin ))
      && (ntol >= fabs( veg.yrltrn - microbe.yrnmin ))
       
      && (ctol >= fabs( yrnep ))
      && (ctol >= fabs( veg.yrnpp - veg.yrltrc ))
      && (ctol >= fabs( veg.yrltrc - microbe.yrrh )) )
  {
      return 1;
  }

  return 0;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void MITTEM44::updateYearSummary( void )
{

  // Update sum of annual carbon storage in ecosystems

  veg.yrcarbon  += y[I_VEGC];

  soil.yrorgc += y[I_SOLC];

  yrtotalc += totalc;

  // Update sum of annual nitrogen storage in ecosystems

  veg.yrstructn += y[I_STRN];

  veg.yrstoren += y[I_STON];

  soil.yrorgn += y[I_SOLN];

  soil.yravln  += y[I_AVLN];

  veg.yrnitrogen  += veg.getVEGN();

  // Update sum of annual phenology in natural ecosystems

  veg.yrunleaf += y[I_UNRMLF];

  veg.yrleaf += y[I_LEAF];

  veg.yrlai += y[I_LAI];

  veg.yrfpc += y[I_FPC];

  // Update sum of annual carbon fluxes in ecosystems

  veg.yringpp += y[I_INGPP];

  veg.yrgpp   += y[I_GPP];

  veg.yrinnpp += y[I_INNPP];

  veg.yrnpp   += y[I_NPP];

  veg.yrgpr   += y[I_GPR];

  veg.yrrmaint += y[I_RVMNT];

  veg.yrrgrowth += y[I_RVGRW];

  veg.yrltrc  += y[I_LTRC];

  microbe.yrrh += y[I_RH];

  yrnep += nep;

  yrnce += nce;

  soil.yrch4csmp += soil.getCH4CONSUMP();

  soil.yrch4ems += soil.getCH4EMISS();

  soil.yrch4flx += soil.getCH4FLUX();

  soil.yrco2dnflx += soil.getCO2DNFLUX();

  soil.yrco2nflx += soil.getCO2NFLUX();
  
  
  // Update sum of annual nitrogen fluxes in ecosystems

  soil.yrnin   += y[I_NINP];

  ag.yrfertn += y[I_AGFRTN];

  veg.yrinnup += y[I_INNUP];

  veg.yrnup   += y[I_VNUP];

  veg.yrsup    += y[I_VSUP];

  veg.yrlup    += y[I_VLUP];

  veg.yrnmobil += y[I_VNMBL];

  veg.yrnrsorb += y[I_VNRSRB];

  veg.yrltrn  += y[I_LTRN];

  microbe.yrnuptake += y[I_MNUP];

  microbe.yrnmin  += y[I_NMIN];

  soil.yrnlost += y[I_NLST];

  soil.yrn2odnflx += soil.getN2ODNFLUX();

  soil.yrn2onflx += soil.getN2ONFLUX();

  soil.yrn2oflx += soil.getN2OFLUX();

  soil.yrn2flx += soil.getN2FLUX();


   // Update sum of annual water fluxes in ecosystems

//  ag.yrirrig += ag.getIRRIG();

  soil.yrineet += soil.getINEET();

  soil.yreet += soil.getEET();

  atms.yrpet += atms.getPET();


  ag.yrstubC += ag.getSTUBBLEC();

  ag.yrstubN += ag.getSTUBBLEN();

 // Update sum of annual carbon and nitrogen fluxes from
 //   agricultural conversion

  ag.yrconvrtC += ag.getCONVRTFLXC();

  ag.yrvconvrtC += ag.getVCONVRTFLXC();

  ag.yrsconvrtC += ag.getSCONVRTFLXC();

  ag.yrslashC += ag.getSLASHC();

  ag.yrcflux += ag.getCFLUX();
  
  ag.yrconvrtN += ag.getCONVRTFLXN();

  ag.yrvconvrtN += ag.getVCONVRTFLXN();

  ag.yrsconvrtN += ag.getSCONVRTFLXN();

  ag.yrslashN += ag.getSLASHN();

  ag.yrnrent += ag.getNRETENT();

  ag.yrnvrent += ag.getNVRETENT();

  ag.yrnsrent += ag.getNSRETENT();
  
  ag.yrformResidueC += ag.getFORMCROPRESIDUEC();
  ag.yrformResidueN += ag.getFORMCROPRESIDUEN();

  ag.yrfluxResidueC += ag.getCROPRESIDUEFLXC();
  ag.yrfluxResidueN += ag.getCROPRESIDUEFLXN();


 // Update sum of annual carbon and nitrogen fluxes from
 //   products

  ag.yrformPROD1C += ag.getCROPPRODC();
  ag.yrformPROD1N += ag.getCROPPRODN();

  ag.yrdecayPROD1C += ag.getPROD1DECAYC();
  ag.yrdecayPROD1N += ag.getPROD1DECAYN();

  ag.yrformPROD10C += ag.getFORMPROD10C();
  ag.yrformPROD10N += ag.getFORMPROD10N();

  ag.yrdecayPROD10C += ag.getPROD10DECAYC();
  ag.yrdecayPROD10N += ag.getPROD10DECAYN();

  ag.yrformPROD100C += ag.getFORMPROD100C();
  ag.yrformPROD100N += ag.getFORMPROD100N();

  ag.yrdecayPROD100C += ag.getPROD100DECAYC();
  ag.yrdecayPROD100N += ag.getPROD100DECAYN();

  ag.yrformTOTPRODC += ag.getFORMTOTPRODC();
  ag.yrformTOTPRODN += ag.getFORMTOTPRODN();

  ag.yrdecayTOTPRODC += ag.getTOTPRODDECAYC();
  ag.yrdecayTOTPRODN += ag.getTOTPRODDECAYN();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITTEM44::urbanDynamics( const int& pdm, double pstate[] )
{
  int fertflag = 0;
  int perennial = 1;                       
  int tillflag = 0;                       
  
  // Initialize fertilizer input to zero
  
  ag.fertn = ZERO;

  // Add N from land conversion to pasture
  
  soil.setNINPUT( ag.getNRETENT() );
  
  
  // Add fireNDEP to ecosystem if fire occurred within the past \
  //   FRI years.  Otherwise do not add any nitrogen to ecosyste
 
  if( firemnthcnt > 0 ) 
  { 
    soil.setNINPUT( ag.getNRETENT() + ag.getFIRENDEP() ); 
  }
  else { soil.setNINPUT( ZERO ); } 

  soil.setEET( soil.getINEET() );

  soil.setVSM( (soil.getMOIST() / (soil.getROOTZ() * 1000.0)) );

  soil.setKH2O( soil.getVSM(), moistlim );


  microbe.updateDynamics( veg.cmnt,
                          soil.getPCTFLDCAP(),
                          pstate[I_SOLC],
                          pstate[I_SOLN],
                          soil.getMOIST(),
                          soil.getVSM(),
                          pstate[I_AVLN],
                          moistlim,
                          tillflag,
                          ag.getTILLFACTOR( veg.cmnt ),
                          soil.getKH2O() );

  if( disturbflag > 1 && pdm == (disturbmonth-1) )
  {
    // Set monthly vegetation fluxes to zero
    	
    veg.resetMonthlyFluxes();  	
  } 
  else
  {
    veg.updateDynamics( ag.cmnt,
                        atms.getCO2(),
                        atms.getAOT40(),
                        soil.getNINPUT(),
                        atms.getPAR(),
                        atms.getPET(),
                        atms.getPRVPETMX(),
                        soil.getEET(),
                        soil.getPRVEETMX(),
                        pstate[I_VEGC],
                        pstate[I_STRN],
                        pstate[I_STON],
                        soil.getMOIST(),
                        pstate[I_AVLN],
                        moistlim,
                        nfeed,
                        o3flag,
                        ag.state,
                        perennial,
                        fertflag,
                        soil.getKH2O(),
                        microbe.getNETNMIN(),
                        ag.fertn );
  }                    

  // Update NINPUT to account for any fertilizer application
  
//  soil.setNINPUT( soil.getNINPUT() + ag.fertn );
  
  // Determine nitrogen losses from ecosystem

  if( 1 == avlnflag )
  {
    soil.updateNLosses( veg.cmnt,
                        (soil.getSURFRUN() + soil.getDRAINAGE()),
                        pstate[I_AVLN], 
                        soil.getMOIST() );

    if( soil.getNLOST() > pstate[I_AVLN] - veg.getNUPTAKE()
         + microbe.getNETNMIN() + soil.getNINPUT() )
    {
      soil.setNLOST( (pstate[I_AVLN] 
                     - veg.getNUPTAKE() 
                     + microbe.getNETNMIN()
                     + soil.getNINPUT()) );
    }
    if( soil.getNLOST() < ZERO )
    {
      soil.setNLOST( ZERO );
      
      microbe.setNETNMIN( (soil.getNLOST() 
                          + veg.getNUPTAKE() 
                          - soil.getNINPUT()
                          - pstate[I_AVLN]) );
    }                        
  }
  else
  {
    // Do not allow changes in available nitrogen in soil
    //   (used during calibration process)

    soil.setNLOST( (soil.getNINPUT()
                   + microbe.getNETNMIN()
                       - veg.getNUPTAKE()) );
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITTEM44::writesitecd( ofstream& fout, const int& dcmnt )
{
  fout << "<?xml version = \"1.0\"?>" << endl << endl;

  fout << "<siteECD version = \"" << version << "\"" << endl;
  fout << "  site = \"" << sitename << "\"" << endl;
  fout << "  longitude = \"" << sitecol << "\"" << endl;
  fout << "  latitude = \"" << siterow << "\"" << endl;
  fout << "  developedBy = \"" << developer << "\"" << endl;
  fout << "  updated = \"" << updated << "\">" << endl;
  fout << endl;
  fout << "  <community type = \"" << veg.cmnt << "\"" << endl;
  fout << "    description = \"" << description << "\">" << endl;
  fout << endl;
  fout << "    <vegca>" << vegca[dcmnt] << "</vegca>" << endl;
  fout << "    <vegcb>" << vegcb[dcmnt] << "</vegcb>" << endl;
  fout << "    <strna>" << strna[dcmnt] << "</strna>" << endl;
  fout << "    <strnb>" << strnb[dcmnt] << "</strnb>" << endl;
  fout << "    <solca>" << solca[dcmnt] << "</solca>" << endl;
  fout << "    <solcb>" << solcb[dcmnt] << "</solcb>" << endl;
  fout << "    <solna>" << solna[dcmnt] << "</solna>" << endl;
  fout << "    <solnb>" << solnb[dcmnt] << "</solnb>" << endl;
  fout << "    <avlna>" << avlna[dcmnt] << "</avlna>" << endl;
  fout << "    <avlnb>" << avlnb[dcmnt] << "</avlnb>" << endl;
  fout << "    <stona>" << stona[dcmnt] << "</stona>" << endl;
  fout << "    <stonb>" << stonb[dcmnt] << "</stonb>" << endl;
  fout << "    <unleaf12>" << veg.getUNLEAF12( dcmnt );
  fout << "</unleaf12>" << endl;
  fout << "    <initleafmx>" << veg.getINITLEAFMX( dcmnt );
  fout << "</initleafmx>" << endl;
  fout << endl;
  fout << "    <vegcmaxcut>" << veg.getCMAXCUT( dcmnt );
  fout << "</vegcmaxcut>" << endl;
  fout << "    <vegcmax1a>" << veg.getCMAX1A( dcmnt );
  fout << "</vegcmax1a>" << endl;
  fout << "    <vegcmax1b>" << veg.getCMAX1B( dcmnt );
  fout << "</vegcmax1b>" << endl;
  fout << "    <vegcmax2a>" << veg.getCMAX2A( dcmnt );
  fout << "</vegcmax2a>" << endl;
  fout << "    <vegcmax2b>" << veg.getCMAX2B( dcmnt );
  fout << "</vegcmax2b>" << endl;
  fout << "    <vegcfall>" << veg.getCFALL( dcmnt );
  fout << "</vegcfall>" << endl;
  fout << "    <vegkra>" << veg.getKRA( dcmnt );
  fout << "</vegkra>" << endl;
  fout << "    <vegkrb>" << veg.getKRB( dcmnt );
  fout << "</vegkrb>" << endl;
  fout << "    <microbekda>" << microbe.getKDA( dcmnt );
  fout << "</microbekda>" << endl;
  fout << "    <microbekdb>" << microbe.getKDB( dcmnt );
  fout << "</microbekdb>" << endl;
  fout << "    <microbelcclnc>" << microbe.getLCCLNC( dcmnt );
  fout << "</microbelcclnc>" << endl;
  fout << "    <microbepropftos>" << microbe.getPROPFTOS( dcmnt );
  fout << "</microbepropftos>" << endl;
  fout << "    <vegnmaxcut>" << veg.getNMAXCUT( dcmnt );
  fout << "</vegnmaxcut>" << endl;
  fout << "    <vegnmax1a>" << veg.getNMAX1A( dcmnt );
  fout << "</vegnmax1a>" << endl;
  fout << "    <vegnmax1b>" << veg.getNMAX1B( dcmnt );
  fout << "</vegnmax1b>" << endl;
  fout << "    <vegnmax2a>" << veg.getNMAX2A( dcmnt );
  fout << "</vegnmax2a>" << endl;
  fout << "    <vegnmax2b>" << veg.getNMAX2B( dcmnt );
  fout << "</vegnmax2b>" << endl;
  fout << "    <vegnfall>" << veg.getNFALL( dcmnt );
  fout << "</vegnfall>" << endl;
  fout << "    <microbenupa>" << microbe.getNUPA( dcmnt );
  fout << "</microbenupa>" << endl;
  fout << "    <microbenupb>" << microbe.getNUPB( dcmnt );
  fout << "</microbenupb>" << endl;
  fout << "    <soilnloss>" << soil.getNLOSS( dcmnt );
  fout << "</soilnloss>" << endl;
  fout << "    <microbenfixpar>" << microbe.getNFIXPAR( dcmnt );
  fout << "</microbenfixpar>" << endl;
  fout << "    <veginitcneven>" << veg.getINITCNEVEN( dcmnt );
  fout << "</veginitcneven>" << endl;
  fout << "    <vegcnmin>" << veg.getCNMIN( dcmnt );
  fout << "</vegcnmin>" << endl;
  fout << "    <vegc2na>" << veg.getC2NA( dcmnt );
  fout << "</vegc2na>" << endl;
  fout << "    <vegc2nb>" << veg.getC2NB( dcmnt );
  fout << "</vegc2nb>" << endl;
  fout << "    <vegc2nmin>" << veg.getC2NMIN( dcmnt );
  fout << "</vegc2nmin>" << endl;
  fout << "    <microbecnsoil>" << microbe.getCNSOIL( dcmnt );
  fout << "</microbecnsoil>" << endl;
  fout << "    <o3para>" << veg.getO3PARA( dcmnt );
  fout << "</o3para>" << endl;
  fout << "    <o3parb>" << veg.getO3PARB( dcmnt );
  fout << "</o3parb>" << endl;
  fout << "    <o3parc>" << veg.getO3PARC( dcmnt );
  fout << "</o3parc>" << endl;

  fout << endl;
  fout << "  </community>" << endl;
  fout << "</siteECD>" << endl;


};








