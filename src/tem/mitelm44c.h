/* **************************************************************
*****************************************************************
MITELM44C.H - Determines the interface among MITCLM, LCLUC and 
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
20070129 - DWK changed include mitelmntcohort437.hpp to
           mitelmntcohort44.hpp
20070129 - DWK changed include tlcluc437.h to tlcluc44.h
20070129 - DWK changed include mittem438.h to mittem44.h
20070129 - DWK changed class MITelmnt43 to class MITelmnt44
20080130 - DWK changed include from mitclm437.h to miclm44a.h
20080130 - DWK changed include from mitelmntcohort44.hpp to
           mitelmntcohort44a.hpp
20080130 - DWK changed include from tlcluc44.h to tigsmlcluc44a.h
20080130 - DWK changed include from mittem44.h to mittem44a.h
20080130 - DWK changed include from mitsoldat437.h to 
           mitsoldat44a.h
20080130 - DWK changed include from mitsolyrdat437.h to
           mitsolyrdat44a.h
20080130 - DWK changed include from ttemdat437.h to 
           tigsmtemdat44a.h
20080130 - DWK changed public MITclm43 mitclm to MITclm44 mitclm
20080826 - DWK changed public string region to 
           string region[MAXCHRTS]
20080827 - DWK changed include from mittem44a.h to mittem44b.h
20110705 - DWK changed include from mitelmntcohort44a.hpp to
           mitelmntcohort44c.hpp
20110706 - DWK added public function landuseChange()
20110707 - DWK changed include from tlcluc44a.h to tigsmlcluc44c.h
20110707 - DWK changed include from mittem44b.h to mittem44c.h
20110707 - DWK changed include from mitclm44a.h to mitclm44c.h
20110707 - DWK changed include from tigsmtemdat44a.h to 
           tigsmtemdat44c.h
20110707 - DWK changed include from mitsoldat44a.h to 
           mitsoldat44c.h
20110707 - DWK changed include from mitsolyrdat44a.h to
           mitsolyrdat44c.h
20110817 - DWK added additional temgisqc()
20110826 - DWK added public double initPET[MAXCHRTS][CYCLE]
20110911 - DWK added private functions createCohortProducts() 
           and updateChangedCohort() 
             
*****************************************************************
************************************************************** */

#ifndef MITELM44C_H
#define MITELM44C_H

const int TQCZEROFLAG = 31;

//Modules representing climate and TEM

#include "mitclm44c.h"   // MITelmnt44 uses the MITclm44 class

// MITelmnt43 uses the MITElmntcohort44 
#include "mitelmntcohort44c.hpp"  

#include "tigsmlcluc44c.h" // MITelmnt44 uses the TEMlcluc44 class
#include "mittem44c.h"     // MITelmnt44 uses the MITTEM44 class

// Modules describing the interface with spatially explicit 
//   data sets

// #include "tigsmclmdat44a.h"  //MITelmnt44 uses the Clmdata44 class

//MITelmnt44 uses the MITSoildata44 class
#include "mitsoldat44c.h"  

//MITelmnt44 uses the MITSoilLayerdata44 class
#include "mitsolyrdat44c.h"  

#include "telvdat437.h"  //MITelmnt44 uses the Elevdata43 class
#include "tigsmtemdat44c.h"  //MITelmnt44 uses the Temdata44 class


class MITelmnt44
{

  public:

     MITelmnt44();


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

     void landUseChange( const int& pdyr );
	
//     void loadMITCLM( ifstream& ifile, MITdata43& clm );

     void outputTEMmonth( const int& pchrt,
                          const int& pdm );

     void readCohortState( ifstream& ifstate,
                           const int& pichrt );
     
     void saveTEMCohortState( const int& pichrt );

     void setCohortTEMState( const MITElmntCohort44& firstchrt,
                             MITElmntCohort44& targetchrt );
     
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
                     MITdata44 cflx[MAXRTIME][MXMITNLAT] );


/* *************************************************************
		 Public Variables
************************************************************* */

     int atmstotyr[MAXRTIME];

     // Mean annual air temperature
     double avetair;
     
     long carea;

     double climate[NUMATMS][CYCLE];
     

     MITElmntCohort44 cohort[MAXCHRTS];
     
     float col;

     string contnent;

     int dbug;
     
//     double elev;
     
     int fatalerr;

     double initAET[MAXCHRTS][CYCLE];

     double initDRAINAGE[MAXCHRTS][CYCLE];

     double initPET[MAXCHRTS][CYCLE];

     double initSH2O[MAXCHRTS][CYCLE];

     double initSNOWPACK[MAXCHRTS][CYCLE];

     double initSURFRUN[MAXCHRTS][CYCLE];

     double lat;

     TEMlcluc44 lcluc;

     double lon;

     int lonlatflag;
     
     // Maximum number of cohorts in an element
     //   during the current year
     int maxcohorts;
     
//     int mez;
     
     MITclm44 mitclm;

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

     string region[MAXCHRTS];
     
     float row;

     long subarea;

     MITTEM44 tem;

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

     void createCohortProducts( const int& pichrt,
                                const double& pproparea,
                                const Biomass& pformPROD10,
                                const Biomass& pformPROD100,
                                const Biomass& pslash,
                                const Biomass& pvconvrtflx,
                                const Biomass& psconvrtflx,
                                const double& pnvretent,
                                const double& pnsretent );

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
                   const double& pet,
                   const double& eet,
                   const double& sh2o,
                   const double& spack,
                   const double& srfrun,
                   const double& drain,
                   const double& co2,
                   const double& aot40 );

     int temgisqc( const long& subarea,
                   const double& pctsilt,
                   const double& pctclay,
                   const int& cmnt,
                   const double& elev,
                   const double& nirr,
                   const double& par,
		   const double& tair,
                   const double& prec,
                   const double& pet,
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

     void updateChangedCohort( const int& pichrt,
                               const long& pinitarea,
                               const double& ppropConvertedArea,
                               const double& ppropAbandonedArea,
                               Tveg44& pabandonedVeg,
                               MITsoil44& pconvertedSoil,
                               MITsoil44& pabandonedSoil,
                               const double& pconvertedNONSOLC,
                               const double& pabandonedNONSOLC,
                               double pconvertedANH4IN[NLVL],
                               double pabandonedANH4IN[NLVL],
                               double pconvertedANO3IN[NLVL],
                               double pabandonedANO3IN[NLVL],
                               double pconvertedDPHUMIN[NLVL],
                               double pabandonedDPHUMIN[NLVL],
                               double pconvertedOCIN[NLVL],
                               double pabandonedOCIN[NLVL],
                               double pconvertedRCLIN[NLVL],
                               double pabandonedRCLIN[NLVL],
                               double pconvertedRCRIN[NLVL],
                               double pabandonedRCRIN[NLVL],
                               double pconvertedRCVLIN[NLVL],
                               double pabandonedRCVLIN[NLVL] );


};

#endif
