/* *************************************************************
****************************************************************
MITTEM44C.H - Terrestrial Ecosystem Model Version 4.4
****************************************************************

Modifications:

20060128 - DWK created by modifying ttem437.h
20060128 - DWK changed include from temconsts43.hpp to
           newconsts.hpp
20060128 - DWK changed include from tsoil437.h to mitsoil437.h
20060128 - DWK changed include from tmcrb437.h to mitmcrb437.h
20060128 - DWK changed class TTEM43 to class MITTEM43
20060128 - DWK moved I_AVLW, I_SM, I_VSM, I_PCTP, I_RGRW,     
           I_SGRW, to a later position in enum temkey and 
           changed value of MAXWSTAT, NUMWEQ and NUMTEM in
           file temconsts.43.hpp to indicate that water pools
           and fluxes were no longer calculated in the RFK
           integrator
20080128 - DWK added I_CH4EMS, I_CH4CSMP, I_CH4FLX, I_CO2NFLX,  
           I_CO2DNFLX, I_NOFLX, I_N2OFLX, I_N2ONFLX, I_N2ODNFLX
           and I_N2FLX to enum temkey 
20060128 - DWK added const int equil to function call of 
           initrun()
20060128 - DWK added public int ch4flag and int n2oflag
20060202 - DWK added I_NIRR, I_PAR, I_TAIR, I_PREC, I_RAIN,   
           I_SNWFAL, I_CO2, I_AOT40 to temkey
20060306 - DWK added pubic int mdmflag
20070129 - DWK changed include from mitsoil437.h to mitsoil44.h
20070129 - DWK changed include from humnact438.h to humnact44.h
20070129 - DWK changed class MITTEM43 to class MITTEM44
20070201 - DWK changed include from tveg438.h to tveg44.h
20070201 - DWK changed public Tveg43 veg to Tveg44 veg
20070201 - DWK changed include from mitmcrb437.h to mitmcrb44.h
20070201 - DWK changed public Tmicrobe43 microbe to
           Tmicrobe44 microbe
20080130 - DWK changed include from nemconsts.hpp to
           nemconsts44a.hpp
20080130 - DWK changed include from mitatms437.h to mitatms44a.h
20080130 - DWK changed include from tveg44.h to tigsmveg44a.h
20081030 - DWK changed include from mitsoil44.h to mitsoil44a.h
20080130 - DWK changed include from mitmcrb44.h to mitmcrb44a.h 
20080130 - DWK changed include from humnact44.h to 
           tigsmhumnact44a.h
20080130 - DWK changed public MITatms43 atms to MITatms44 atms
20080826 - DWK changed include from tigsmveg44a.h to 
           tigsmveg44b.h
20080826 - DWK changed include from tigsmhumnact44a.h to 
           tigsmhumnact44b.h
20110706 - DWK changed include from tigsmhumnact44b.h to 
           tigsmhumnact44c.h
20110707 - DWK changed include from tigsmveg44b.h to 
           tigsmveg44c.h
20110707 - DWK changed include from nemconsts44a.hpp to
           nemconsts44c.hpp
20110707 - DWK changed include from mitatms44a.h to 
           mitatms44c.h
20110707 - DWK changed include from mitsoil44a.h to mitsoil44c.h
20110707 - DWK changed include from mitmcrb44a.h to mitmcrb44c.h
20140604 - DWK changed access to "yr" variables from public to
           private
                                                                   
****************************************************************

References:

VERSION 4.1

Tian, H., J.M. Melillo, D.W. Kicklighter, A.D. McGuire and J.
  Helfrich.  1999. The sensitvity of terrestrial carbon storage to
  historical climate variability and atmospheric CO2 in the United
  States.  Tellus 51B: 414-452.

VERSION 4.2

McGuire, A.D., S. Sitch, J.S. Clein, R. Dargaville, G. Esser, J. Foley,
  M. Heimann, F. Joos, J. Kaplan, D.W. Kicklighter, R.A. Meier, J.M.
  Melillo, B. Moore III, I.C. Prentice, N. Ramankutty, T. Reichenau,
  A. Schloss, H. Tian, L.J. Williams and U. Wittenberg (2001) Carbon
  balance of the terrestrial biosphere in the twentieth century:
  analyses of CO2, climate and land use effects with four process-
  based ecosystem models.  Global Biogeochemical Cycles 15: 183-206.

Tian, H., J.M. Melillo, D.W. Kicklighter, S. Pan, J. Liu, A.D. McGuire
  and B. Moore III (2003) Regional carbon dynamics in monsoon Asia
  and its implications for the global carbon cycle. Global and
  Planetary Change 37: 201-217.

VERSION 4.3

Felzer, B., D. Kicklighter, J. Melillo, C. Wang, Q. Zhuang, and 
  R. Prinn (2004) Effects of ozone on net primary production and 
  carbon sequestration in the conterminous United States using a 
  biogeochemistry model. Tellus 56B: 230-248.

Runge - Kutta - Fehlberg (RKF) adaptive integrator:

Cheney, W., and D. Kincaid.  1985.  Numerical mathematics and
  computing.  Second edition.  Brooks/ Col Publishing Co.  Monterey,
  CA.  pp. 325-328.

****************************************************************
************************************************************** */

#ifndef MITTEM44C_H
#define MITTEM44C_H
 
#include "nemconsts44c.hpp"


#ifdef CALIBRATE_TEM
  // Additional global constants

  const int WSY = 5;
  const int ESY = 3;
#endif


// Objects describing basic components of the ecosystem

#include "bioms423.hpp"   // MITTTEM44 uses Biomass class


// Objects describing the structure of the ecosystem

#include "mitatms44c.h"  // MITTEM44 uses MITatms44 class
#include "tigsmveg44c.h" // MITTEM44 uses Tveg44 class
#include "mitsoil44c.h"  // MITTEM44 uses MITsoil44 class
#include "mitmcrb44c.h"  // MITTEM44 uses Tmicrobe44 class

// Objects describing the effects of human activities on the 
//   ecosystem

#include "tigsmhumnact44c.h" // MITTEM44 uses Humnact44 class


class MITTEM44
{

  public:

     MITTEM44();

     enum temkey
     {
       I_VEGC,     I_SOLC,          
       
       I_STRN,     I_STON,     I_SOLN,     I_AVLN,
                        
       I_UNRMLF,   I_LEAF,     I_LAI,      I_FPC,

       I_INGPP,    I_GPP,      I_FOZONE,   I_FINDOZONE, I_INNPP,    
       I_NPP,      I_GPR,      I_RVMNT,    I_RVGRW,     I_LTRC,     
       I_RH,                    

       I_NINP,     I_AGFRTN,   I_INNUP,    I_VNUP,      I_VSUP,     
       I_VLUP,     I_VNMBL,    I_VNRSRB,   I_LTRN,      I_MNUP,     
       I_NMIN,     I_NLST,
       
       I_AGIRRIG,  I_INEET,    I_EET,          
                 

       I_TOTEC,    I_TOTC, 
       
       I_VEGN,        
       
       I_SNWPCK,   I_AVLW,     I_SM,       I_VSM,       I_PCTP,      
       

       
       I_NEP,      I_NCE,      I_CH4EMS,   I_CH4CSMP,   I_CH4FLX,   
       I_CO2NFLX,  I_CO2DNFLX,   
       
       I_NOFLX,    I_N2OFLX,   I_N2ONFLX,  I_N2ODNFLX, I_N2FLX, 
       
       I_PET,      
       
       I_AGPRDC,   I_PROD10C,  I_PROD100C, I_TOTPRDC,  
       
       I_RESIDC,   I_AGSTUBC,     
       
       I_AGPRDN,   I_PROD10N,  I_PROD100N, I_TOTPRDN,
       
       I_RESIDN,   I_AGSTUBN,
       
       I_CNVRTC,   I_VCNVRTC,  I_SCNVRTC,  I_SLASHC,    I_CFLX,     
       
       I_CNVRTN,   I_VCNVRTN,  I_SCNVRTN,  I_SLASHN,    I_NRETNT,
       I_NVRTNT,   I_NSRTNT,

       I_AGFPRDC,  I_AGFPRDN,  I_FRESIDC,  I_FRESIDN,   I_AGPRDFC,
       I_AGPRDFN,  I_RESIDFC,  I_RESIDFN,  
       
       I_PRDF10C,  I_PRDF10N,  I_PRD10FC,  I_PRD10FN,   I_PRDF100C,  
       I_PRDF100N, I_PRD100FC, I_PRD100FN, I_TOTFPRDC,  I_TOTFPRDN, 
       I_TOTPRDFC, I_TOTPRDFN,

       I_CROPC,    I_NATVEGC,  I_CROPN,    I_NATVEGN,   I_CSTRN,
       I_NATSTRN,  I_CSTON,    I_NATSTON,

       I_CROPULF,  I_NATULF,   I_CROPLEAF, I_NATLEAF,   I_CROPLAI,
       I_NATLAI,   I_CROPFPC,  I_NATFPC,

       I_AGINGPP,  I_NATINGPP, I_AGGPP,    I_NATGPP,    I_AGINNPP,
       I_NATINNPP, I_AGNPP,    I_NATNPP,   I_AGGPR,     I_NATGPR,
       I_AGRVMNT,  I_NATRVMNT, I_AGRVGRW,  I_NATRVGRW,  I_AGLTRC,
       I_NATLTRC,

       I_AGINNUP,  I_NATINNUP, I_AGVNUP,   I_NATVNUP,   I_AGVSUP,
       I_NATVSUP,  I_AGVLUP,   I_NATVLUP,  I_AGVNMBL,   I_NATVNMBL,
       I_AGVNRSRB, I_NVNRSRB,  I_AGLTRN,   I_NATLTRN,
       
       I_NIRR,     I_PAR,      I_TAIR,     I_PREC,      I_SRFRUN,   
       I_DRAIN,    I_CO2,      I_AOT40
     };

     #ifdef CALIBRATE_TEM    
       enum seykey { NOEKEY,       GET_LEAF,    GET_LAI,     
                     GET_FPC,      

                     GET_INGPP,    GET_GPP,     GET_D40,
                     GET_FOZONE,   GET_INNPP,   GET_NPP,      
                     GET_GPR,     GET_RVMNT,    GET_RVGRW,    
                     GET_LTRC,    GET_AGSTUBC,  GET_RH,       
                     GET_NEP,     

                     GET_NINP,     GET_AGFRTN,  GET_INNUP,   
                     GET_VNUP,     GET_VSUP,    GET_VLUP,     
                     GET_VNMBL,    GET_VNRSRB,   GET_LTRN,     
                     GET_AGSTUBN,  GET_MNUP,    GET_NMIN, 
                     GET_NLST,    

                     GET_CNVRTC,   GET_VCNVRTC, GET_SCNVRTC,
                     GET_SLASHC,   GET_CFLX,    GET_NCE,

                     GET_CNVRTN,   GET_VCNVRTN, GET_SCNVRTN,
                     GET_SLASHN,   GET_NRETNT,  GET_NVRTNT,
                     GET_NSRTNT,

                     GET_AGPRDC,   GET_PROD10C, GET_PROD100C,
                     GET_RESIDC,

                     GET_AGPRDN,   GET_PROD10N, GET_PROD100N,
                     GET_RESIDN,

                     GET_AGFPRDC,  GET_PRDF10C, GET_PRDF100C,
                     GET_FRESIDC,

                     GET_AGPRDFC,  GET_PRD10FC, GET_PRD100FC,
                     GET_TOTPRDFC, GET_RESIDFC,

                     GET_AGFPRDN,  GET_PRDF10N, GET_PRDF100N,
                     GET_FRESIDN,
                     GET_AGPRDFN,  GET_PRD10FN, GET_PRD100FN,
                     GET_TOTPRDFN, GET_RESIDFN,

                     GET_L2SN };


       enum swykey { NOWKEY,    GET_SH2O,   GET_PCTP,   GET_VSM,

                     GET_RAIN,  GET_SNWFAL, GET_SNWINF, GET_AGIRRIG,
                     GET_PET,   GET_INEET,  GET_EET,    GET_RPERC,
                     GET_SPERC, GET_RRUN,  GET_SRUN,    GET_WYLD };

     #endif

/* **************************************************************
			Public Functions
************************************************************** */

     void askODE( ofstream& rflog1 );


     #ifdef CALIBRATE_TEM
       void displayOptionalEflx( const seykey& s );
       void displayOptionalWflx( const swykey& s );
     #endif

    
     int ecdqc( const int& dcmnt );

     void ECDsetODEstate( const int& pdcmnt,
                          const double& psiplusc );
 
     void getco2( void );

     void getcropecd( const int& dv, const string& agecd );


     #ifdef CALIBRATE_TEM
       double getOptionalEflx( const int& optflx );

       double getOptionalWflx( const int& optflx );
     #endif


     void getsitecd( const int& numcmnt, ofstream& rflog1 );

     void getsitecd( const int& dv, const string& ecd );
     
     void initrun( ofstream& rflog1 );

     int monthlyTransient( const int& outyr,
                           const int& pdm,
                           const double& outtol );


     #ifdef CALIBRATE_TEM
       inline seykey& next( seykey& s )
       {
         return s = (GET_L2SN == s) ? GET_LEAF : seykey( s+1 );
       }

       inline swykey& next( swykey& s )
       {
         return s = (GET_WYLD == s) ? GET_SH2O : swykey( s+1 );
       }

       inline seykey& prev( seykey& s )
       {
         return s = (GET_LEAF == s) ? GET_L2SN : seykey( s-1 );
       }

       inline swykey& prev( swykey& s )
       {
         return s = (GET_SH2O == s) ? GET_WYLD : swykey( s-1 );
       }
     #endif


     void resetMonthlyELMNTFluxes( void );

     void resetYrFluxes( void );

     void setELMNTecd( const int& pdcmnt,
                       const double& psiplusc );

     void setEquilC2N( const int& pdcmnt,
                       const double& co2 );

     void setEquilEvap( const double& nirr,
                        const double& tair,
                        const int& pdm );

     void setPrevState( void );
         
     int stepmonth( const int& dyr,
                    const int& pdm,
                    int& intflag,
                    const double& tol );

     int testEquilibrium( void );

     void updateYearSummary( void );
     
     void writesitecd( ofstream& fout, const int& dcmnt );


     // "Get" and "Set" private variables and parameters   
     
     // avlna **************************************************
     
     inline double getAVLNA( const int& pcmnt )
     {
       return avlna[pcmnt];
     }

     inline void setAVLNA( const double& pavlna, 
                           const int& pcmnt )
     {
       avlna[pcmnt] = pavlna;
     }


     // avlnb **************************************************
     
     inline double getAVLNB( const int& pcmnt )
     {
       return avlnb[pcmnt];
     }

     inline void setAVLNB( const double& pavlnb, 
                           const int& pcmnt )
     {
       avlnb[pcmnt] = pavlnb;
     }


    // lat ****************************************************
     
     inline float getLAT( void ) { return lat; }

     inline void setLAT( const float& plat )
     {
       lat = plat;
     }

     // nce ****************************************************
     
     inline double getNCE( void ) { return nce; }


     // nep ****************************************************
          
     inline double getNEP( void ) { return nep; }


     // prevy **************************************************
     
     inline double getPREVY( const int& i ) { return prevy[i]; }

     inline void setPREVY( const double& pprevy, const int& i )
     {
       prevy[i] = pprevy;
     }
     

     // solca **************************************************
     
     inline double getSOLCA( const int& pcmnt )
     {
       return solca[pcmnt];
     }

     inline void setSOLCA( const double& psolca, 
                           const int& pcmnt )
     {
       solca[pcmnt] = psolca;
     }


     // solcb **************************************************
     
     inline double getSOLCB( const int& pcmnt )
     {
       return solcb[pcmnt];
     }

     inline void setSOLCB( const double& psolcb, 
                           const int& pcmnt )
     {
       solcb[pcmnt] = psolcb;
     }

 
     // solna **************************************************
     
     inline double getSOLNA( const int& pcmnt )
     {
       return solna[pcmnt];
     }

     inline void setSOLNA( const double& psolna, 
                           const int& pcmnt )
     {
       solna[pcmnt] = psolna;
     }


     // solnb **************************************************
     
     inline double getSOLNB( const int& pcmnt )
     {
       return solnb[pcmnt];
     }

     inline void setSOLNB( const double& psolnb, 
                           const int& pcmnt )
     {
       solnb[pcmnt] = psolnb;
     }


     // stona **************************************************
     
     inline double getSTONA( const int& pcmnt )
     {
       return stona[pcmnt];
     }

     inline void setSTONA( const double& pstona, 
                           const int& pcmnt )
     {
       stona[pcmnt] = pstona;
     }


     // stonb **************************************************
     
     inline double getSTONB( const int& pcmnt )
     {
       return stonb[pcmnt];
     }

     inline void setSTONB( const double& pstonb, 
                           const int& pcmnt )
     {
       stonb[pcmnt] = pstonb;
     }


     // strna **************************************************
     
     inline double getSTRNA( const int& pcmnt )
     {
       return strna[pcmnt];
     }

     inline void setSTRNA( const double& pstrna, 
                           const int& pcmnt )
     {
       strna[pcmnt] = pstrna;
     }


     // strnb **************************************************
     
     inline double getSTRNB( const int& pcmnt )
     {
       return strnb[pcmnt];
     }

     inline void setSTRNB( const double& pstrnb, 
                           const int& pcmnt )
     {
       strnb[pcmnt] = pstrnb;
     }


     // totalc *************************************************
     
     inline double getTOTALC( void ) { return totalc; }


     // vegca **************************************************
     
     inline double getVEGCA( const int& pcmnt )
     {
       return vegca[pcmnt];
     }

     inline void setVEGCA( const double& pvegca, 
                           const int& pcmnt )
     {
       vegca[pcmnt] = pvegca;
     }


     // vegcb **************************************************
     
     inline double getVEGCB( const int& pcmnt )
     {
       return vegcb[pcmnt];
     }

     inline void setVEGCB( const double& pvegcb, 
                           const int& pcmnt )
     {
       vegcb[pcmnt] = pvegcb;
     }


     // y ******************************************************
          
     inline double getY( const int& i ) { return y[i]; }

     inline void setY( const double& py, const int& i )
     {
       y[i] = py;  
     }

     // yrnce ************************************************
     
     inline double getYRNCE( void )
     {
       return yrnce;
     }

     inline void setYRNCE( const double& pyrnce )
     {
       yrnce = pyrnce;
     }

     inline void updateYRNCE( const double& pnce )
     {
       yrnce += pnce;
     }

    // yrnep ************************************************
     
     inline double getYRNEP( void )
     {
       return yrnep;
     }

     inline void setYRNEP( const double& pyrnep )
     {
       yrnep = pyrnep;
     }

     inline void updateYRNEP( const double& pnep )
     {
       yrnep += pnep;
     }

    // yrtotalc ************************************************
     
     inline double getYRTOTALC( void )
     {
       return yrtotalc;
     }

     inline void setYRTOTALC( const double& pyrtotalc )
     {
       yrtotalc = pyrtotalc;
     }

     inline void updateYRTOTALC( const double& ptotalc )
     {
       yrtotalc += ptotalc;
     }

                          
/* **************************************************************
			 Public Variables
************************************************************** */

     // Input Ecd files with calibration data
     
     #ifdef CALIBRATE_TEM
       string soilfile;
       string rootfile;
       string vegfile;
       string leaffile;
       string mcrvfile;
       string agfile;
     #endif

     #ifdef CALIBRATE_TEM
       int adapttol;
     #endif

    
     Humanact44 ag;

     MITatms44 atms;

     static int avlnflag;
     
     static int baseline;


     #ifdef CALIBRATE_TEM
       int calwind;
     #endif


     // Flag for NEM and MDM: 
     //   0 = no CH4 fluxes calculated
     //   1 = calculate CH4 fluxes
     int ch4flag;     
     
     static double ctol;

     int dbug;

     static int diffyr;

     // Counter of months since disturbance
     int distmnthcnt;
     
     // Flag to indicate that a disturbance has occurred
     //   in a particular year.  Different values represent
     //   different types of disturbances:
     //   0 = no disturbance
     //   1 = sustained disturbed state 
     //       (e.g., conversion to crop, pasture or urban areas)
     //   2 = timber harvest
     //   3 = fire
     int disturbflag;
     
     // Month in which disturbance occurred
     // (e.g. 1 = January, ...., 12 = December)
     int disturbmonth;   

     double elev;              // elevation (m)

     static int endeq;

     static int endyr;

     static int equil;

     // index for vegetation type (ez = veg.temveg-1)
     int ez;

     // Counter of months since fire disturbance
     int firemnthcnt;

     static int initbase;

     int initFlag;
     
     static double inittol;


     #ifdef CALIBRATE_TEM
       int intbomb;
     #endif


     static int intflag;

     static int maxnrun;

     static int maxyears;

     // Flag to indicate whether or not to run the Methane 
     //   Dynamics Module ( 0 = no, 1 = yes)
     int mdmflag;
     
     Tmicrobe44 microbe;

     static int moistlim;

     // Flag for NEM: 
     //   0 = no N2O fluxes calculated
     //   1 = calculate N2O fluxes
     int n2oflag;

     int nattempt;

     static int nfeed;

     static double ntol;

     static int o3flag;
     
     int predflag;

     vector<string> predstr;

     int qualcon[MAXRTIME];

     int retry;

     static int rheqflag;

     static int runsize;


     #ifdef CALIBRATE_TEM
       seykey sey[ESY];
       swykey swy[WSY];
     #endif

    
     MITsoil44 soil;

     static int startyr;

     static int strteq;

     double tol;


     #ifdef CALIBRATE_TEM
       int tolbomb;   
       int topwind;
     #endif


     int totyr;

     Tveg44 veg;

     static int wrtyr;

     static double wtol;


     // Site ECD variables

     string version;
     string sitename;
     string developer;
     float sitecol;
     float siterow;
     long updated;
     string description;


/* **************************************************************
		     Protected Variables
************************************************************** */

  protected:

     int dbugflg;


  private:

/* **************************************************************
		 Private Functions
************************************************************** */

     int adapt( const int& numeq, 
                double pstate[], 
                const double& ptol,
                const int& pdm );

     int boundcon( double ptstate[],
                   double err[],
                   const double& ptol );

     void delta( const int& dm,
                 double pstate[],
                 double pdstate[] );

     void cropDynamics( double pstate[] );

     void getenviron( const int& pdm );

     void massbal( double y[NUMEQ], 
                   double prevy[MAXSTATE] );

     void natvegDynamics( const int& dm, double pstate[] );

     void pastureDynamics( const int& dm, double pstate[] );


     #ifdef CALIBRATE_TEM
       void pcdisplayDT( const double& tottime,
                         const double& deltat );

       void pcdisplayMonth( const int& dyr,
                            const int& dm );

       void pcdisplayODEerr( const int& test, 
                             double pstate[] );
     #endif


     void resetODEflux( void );

     void rkf( const int& numeq, 
               double pstate[], 
               double& pdt, 
               const int& pdm );
                 
     void step( const int& numeq, 
                double pstate[], 
                double pdstate[], 
                double ptstate[], 
                double& pdt );

     void urbanDynamics( const int& dm, double pstate[] );


/* **************************************************************
		     Private Variables
************************************************************** */

     ifstream fecd[MAXCMNT];

     // Latitude
     float lat;

     // Net Carbon Exchange (NCE)
     double nce;      // (g C / (sq. meter * month))

     // Net Ecosystem Production
     double nep;      // (g C / (sq. meter * month))

     // Values of ODE state variables for carbon, nitrogen and
     //   water pools during the previous month 
     double prevy[MAXSTATE];

     // Total carbon storage (veg.plant + soil.org)
     double totalc;   // (g C/sq. meter)

     // Values of ODE state variables for current month
     double y[NUMEQ];

     double yrnce;

     double yrnep;    // (g C / (sq. meter * year))
    
     double yrtotalc;

    
     // Adaptive integrator variables
     
     int blackhol;
     static int maxit;
     static long maxitmon;

     int syint;
     int test;

     double dum4[NUMEQ];
     double error[NUMEQ];

     double dum5[NUMEQ];
     double dumy[NUMEQ];

     double ydum[NUMEQ];
     double dy[NUMEQ];
     double yprime[NUMEQ];
     double rk45[NUMEQ];

     double f11[NUMEQ];
     double f2[NUMEQ];
     double f13[NUMEQ];
     double f3[NUMEQ];
     double f14[NUMEQ];
     double f4[NUMEQ];
     double f15[NUMEQ];
     double f5[NUMEQ];
     double f16[NUMEQ];
     double f6[NUMEQ];


     static double  a1;

     static double  a3;
     static double a31;
     static double a32;

     static double  a4;
     static double a41;
     static double a42;
     static double a43;

     static double  a5;
     static double a51;
     static double a52;
     static double a53;
     static double a54;

     static double  b1;
     static double  b3;
     static double  b4;
     static double  b5;

     static double  b6;
     static double b61;
     static double b62;
     static double b63;
     static double b64;
     static double b65;

     double dummy;


/* **************************************************************
			Private Parameters
************************************************************** */

     double vegca[MAXCMNT];
     double vegcb[MAXCMNT];

     double strna[MAXCMNT];
     double strnb[MAXCMNT];

     double stona[MAXCMNT];
     double stonb[MAXCMNT];

     double solca[MAXCMNT];
     double solcb[MAXCMNT];

     double solna[MAXCMNT];
     double solnb[MAXCMNT];

     double avlna[MAXCMNT];
     double avlnb[MAXCMNT];

 };

#endif
