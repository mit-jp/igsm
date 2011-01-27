/* **************************************************************
*****************************************************************
MITELM438.H - Determines the interface among MITCLM, LCLUC and 
                TEM modules within the MIT IGSM

Modifications:

20060129 - DWK created by modifying telm437.h
20060129 - DWK changed include from telmntcohort437.hpp to
           mitelmntcohort437.hpp 
20060129 - DWK changed include from tclm437.h to mittclm437.h
20060129 - DWK changed include from ttem437.h to mittem437.h
20060129 - DWK deleted include tclmdat437.h
20060129 - DWK changed include from tsoldat437.h to 
           mittsoldat437.h
20060129 - DWK changed Class TEMelmnt43 to Class MITelmnt43
20060129 - DWK added public function void loadMITCLM()
20060129 - DWK added public function void writecflx()
20060129 - DWK changed public TEMclm43 clm to MITclm43 clm
20060127 - DWK changed public ElmntCohort43 cohort[MAXCHRTS] to 
           MITElmntCohort43 cohort[MAXCHRTS]
20060129 - DWK added public ofstream ofstateA, ofstream ofstateB,
           string temstateAfname and string temstateBfname
20060129 - DWK changed public double output[NUMTEM][CYCLE]
           to double output[NUMTEM][MAXCHRTS][CYCLE]
20060129 - DWK changed public TTEM43 tem to MITTEM43 tem
20060129 - DWK added const long& subarea, const double& pet,
           const double& eet and const double& sh2o to function
           call of temgisqc()
20060129 - DWK deleted const double& avtair from function call
           of temgisqc()
20060202 - DWK deleted public functions atmswritemiss() 
           and atmswritepred()
20060320 - DWK added public double initDRAINAGE[MAXCHRTS][CYCLE],
           double initSNOWPACK[MAXCHRTS][CYCLE] and
           double initSURFRUN[MAXCHRTS][CYCLE]

*****************************************************************
************************************************************** */

#ifndef MITELM438_H
#define MITELM438_H

const int TQCZEROFLAG = 31;

//Modules representing climate and TEM

#include "mitclm437.h"   // MITelmnt43 uses the MITclm43 class

// MITelmnt43 uses the MITElmntcohort43 
#include "mitelmntcohort437.hpp"  

#include "tlcluc437.h" // MITelmnt43 uses the TEMlcluc43 class
#include "mittem438.h"   // MITelmnt43 uses the MITTEM43 class

// Modules describing the interface with spatially explicit 
//   data sets

// #include "tclmdat437.h"  //MITelmnt43 uses the Clmdata43 class

//MITelmnt43 uses the MITSoildata43 class
#include "mitsoldat437.h"  

//MITelmnt43 uses the MITSoilLayerdata43 class
#include "mitsolyrdat437.h"  

#include "telvdat437.h"  //MITelmnt43 uses the Elevdata43 class
#include "ttemdat437.h"  //MITelmnt43 uses the Temdata43 class


class MITelmnt43
{

  public:

     MITelmnt43();


/* *************************************************************
		 Public Function Declarations
************************************************************* */

     int coregerr( ofstream& rflog1,
                   const string& varname1,
                   const float& col1,
		   const float& row1,
                   const string& varname2,
                   const float& col2,
		   const float& row2 );

     int equilibrateTEM( const int& pchrt,
                         const double& ptol );

     void getTEMCohortState( const int& pichrt );

     void initializeCohortTEMState( const int& pichrt );

//     void loadMITCLM( ifstream& ifile, MITdata43& clm );

     void outputTEMmonth( const int& pchrt,
                          const int& pdm );

     void readCohortState( ifstream& ifstate,
                           const int& pichrt );
     
     void saveTEMCohortState( const int& pichrt );

     void setCohortTEMState( const MITElmntCohort43& firstchrt,
                             MITElmntCohort43& targetchrt );
     
     int setGIStopography( ofstream& flog1,
                           int& ftlerr,
                           FILE* fstxt,
                           FILE* fslayer,
                           FILE*felev );
          
     void setTEMequilState( ofstream& rflog1,
                            const int& equil,
                            const int& totsptime,
                            const int& pichrt );

     void setTEMmiss( const int& pdyr,
                      const int& equil,
                      const int& totsptime,
                      const int& pichrt );

     void temwritepred( ofstream fout[NUMTEM],
                        const vector<string>& predname,
                        const int& pdyr,
                        const int& pichrt,
                        const int& ntempred );

     void updateTEMmonth( const int& equil,
                          const int& totsptime,
                          const int& outyr,
                          const int& pdm,
                          const int& pichrt );

     void writeCohortState( ofstream& ofstate,
                            const int& pichrt );
                            
     void writecflx( ofstream& ofile, 
                     string varname, 
                     const int& dyr,
                     const int& year, 
                     MITdata43 cflx[MAXRTIME][MXMITNLAT] );


/* *************************************************************
		 Public Variables
************************************************************* */

     int atmstotyr[MAXRTIME];

     // Mean annual air temperature
     double avetair;
     
     long carea;

     double climate[NUMATMS][CYCLE];
     

     MITElmntCohort43 cohort[MAXCHRTS];
     
     float col;

     string contnent;
     
//     double elev;
     
     int fatalerr;

     double initAET[MAXCHRTS][CYCLE];

     double initDRAINAGE[MAXCHRTS][CYCLE];

     double initSH2O[MAXCHRTS][CYCLE];

     double initSNOWPACK[MAXCHRTS][CYCLE];

     double initSURFRUN[MAXCHRTS][CYCLE];

     double lat;

     TEMlcluc43 lcluc;

     double lon;

     int lonlatflag;
     
     // Maximum number of cohorts in an element
     //   during the current year
     int maxcohorts;
     
//     int mez;
     
     MITclm43 mitclm;

     // Maximum monthly air temperature during the year
     double mxtair;

     // Maximum number of "natural" cohorts in an element
     int natcohorts;

     // Number of climate predicted variables
     int natmspred;
     
     // Number of TEM predicted variables
     int ntempred;

     ofstream ofstateA; // Save TEMstate every even "outyr"
     ofstream ofstateB; // Save TEMstate every odd "outyr"

//     vector<string> predstr;

     double output[NUMTEM][MAXCHRTS][CYCLE];
     
     int outyr;

     // Maximum number of cohorts in an element
     //   during the previous year
     int prvmxcohrts;

     string region;
     
     float row;

     long subarea;

     MITTEM43 tem;

     string temstateAfname;
     string temstateBfname;

     int totpred;

     int ttotyr[MAXRTIME];

     int wrtyr;

     // Annual precipitation during the year
     double yrprec;
     
     int year;

/* *************************************************************
		 Private Function Declarations
************************************************************* */

  private:

     int temgisqc( const long& subarea,
                   const double& pctsilt,
                   const double& pctclay,
                   const int& cmnt,
                   const double& elev,
                   const double& nirr,
                   const double& par,
		   const double& tair,
                   const double& mxtair,
                   const double& yrprec,
                   const double& prec,
                   const double& eet,
                   const double& sh2o,
                   const double& spack,
                   const double& srfrun,
                   const double& drain,
                   const double& co2,
                   const double& aot40 );

     int transqc( int& maxyears,
                  int& totyr,
                  double plantc[CYCLE] );

};

#endif
