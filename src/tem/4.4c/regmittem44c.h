/* *************************************************************
****************************************************************
REGMITTEM44C.H - Extrapolation version of the Terrestrial 
                  Ecosystem Model Version 4.4 used with 
                  interactive version of TEM (IMITTEM44.CXX)
****************************************************************

Modifications:

20060130 - DWK created by modifying regmittem436a.h
20060130 - DWK added include nemconsts.hpp
20060130 - DWK deleted includes tclmlist436.h, tlulclist436.h,
           and mittemlist436.h
20060130 - DWK changed include from tclmdat436.h to tclmdat437.h           
20060130 - DWK changed include from elmnt436.h to elmnt437.h
20060130 - DWK changed include from latdat436.h to latdat437.h
20060130 - DWK changed include from mitelm436a.h to mitelm437a.h
20060130 - DWK updated public function askpred with string 
           vectors rather than character arrays 
20060130 - DWK replaced the character arrays with strings in 
           public function coregTime()
20060130 - DWK deleted public function createMITCLMregion()
20060130 - DWK commented out public function
20060130 - DWK deleted ClmList43 gridclm, LULCList gridlulc,
           and MITTEMList436 gridtem
20060130 - DWK replaced public char predmap[MAXPRED][9] with
           vector<string> clmpredmap( NUMATMS ) and
           vector<string> tempredmap( NUMTEM )
20060130 - DWK added public int spinoutyrs
20060130 - DWK deleted include regmittem436a.cpp from bottom of
           file
20060201 - DWK deleted public functions initializeLCLUCregion()
20060201 - DWK deleted public vector<string> 
           clmpredmap( NUMATMS )
20060201 - DWK deleted public function setClmPred()
20060202 - DWK deleted updateMITTEMregionYr()
20060302 - DWK added public FILE* fslayer
20060420 - DWK added compiler directive STANDALONE_TEM
20070129 - DWK changed include mitelm438.h to mitelm44.h 
20080130 - DWK changed include from mitelm44.h to mitelm44a.h
20080130 - DWK changed incldue from netconsts.hpp to
           nemconsts44a.hpp
20080130 - DWK changed include from tclmdat437.h to 
           tigsmclmdat44a.h
20080826 - DWK changed include from mitelm44a.h to mitelm44b.h
20110705 - DWK changed include from mitelm44b.h to mitelm44c.h
20110707 - DWK changed include from nemconsts44a.hpp to 
           nemconsts44c.hpp
20110707 - DWK changed include from mitdata44a.h to mitdata44c.h
20110707 - DWK changed include from tigsmclmdat44a.h to 
           tigsmclmdat44c.h
20110707 - DWK changed include from climate4tem44a.hpp to 
           climate4tem44c.hpp
20110707 - DWK changed include from mitclmdata44a.h to 
           mitclmdata44c.h
20110707 - DWK changed include from mitdayclmdata44a.h to 
           mitdayclmdata44c.h
20110707 - DWK changed include from mitdaychrtdata44a.h to 
           mitdaychrtdata44c.h
20110707 - DWK changed include from mithrclmdata44a.h to 
           mithrclmdata44c.h
20110707 - DWK changed include from o3data44a.h to o3data44c.h
20110715 - DWK added public functions getGriddedLULCparameters() 
           and setMITLULCgrid()
20110715 - DWK added include landcover4tem44c.hpp
20110715 - DWK changed access to updateLCLUCregion() to occur 
           only in the stand-alone version of TEM
20110715 - DWK added include landareadat44c.h and lulcchrt44c.h
                                                                                              
****************************************************************
************************************************************* */

#ifndef REGMITTEM44C_H
#define REGMITTEM44C_H

// NOTE: If running TEM "offline", make sure the compiler 
//         directive STANDALONE_TEM is DEFINED below.  If  
//         coupling TEM to the IGSM, make sure STANDALONE_TEM is 
//         NOT DEFINED (i.e. comment out next line)

#include "preproc.h"

#include "nemconsts44c.hpp"

#include "climate4tem44c.hpp"

#include "landcover4tem44c.hpp"

#include "elmnt437.h"        // Elmnt43 Class

#include "latdat437.h"       // Latdata43 class

// MITclm43 uses the MITdata43 class
#include "mitdata44c.h" 

#ifdef STANDALONE_TEM 
  #include "tigsmclmdat44c.h"      // Clmdata43 class
#endif

#ifdef STANDALONE_TEM 
  // MITclm43 uses the MITCLMdata43 class
  #include "mitclmdata44c.h" 
#endif

#ifdef STANDALONE_TEM 
  // MITclm43 uses the MITDayChrtData43 class
  #include "mitdaychrtdata44c.h" 
#endif

#ifdef STANDALONE_TEM 
  // MITclm43 uses the MITDayClmData43 class
  #include "mitdayclmdata44c.h" 
#endif

#ifdef STANDALONE_TEM 
  // MITclm43 uses the MITHrClmData43 class
  #include "mithrclmdata44c.h" 
#endif

#ifdef STANDALONE_TEM 
  // MITclm43 uses the O3data class
  #include "o3data44c.h"   
#endif

#ifdef STANDALONE_TEM 
  #include "landareadat44c.h"   
#endif

#ifdef STANDALONE_TEM 
  #include "lulcchrtdat44c.h"   
#endif


#include "mitelm44c.h"       // MITelmnt Class

class RegMITTEM
{

   public:

     RegMITTEM();
     ~RegMITTEM();


/* *************************************************************
		 Function Declarations
************************************************************* */

     int askpred( const vector<string>& pvarname, 
                  const int& posspred, 
                  vector<string>& predmap );
   
     int coregTime( const string& varname1, 
                    const int& year1, 
                    const string& varname2, 
                    const int& year2 );

     void getGriddedLULCparameters (  const int& pigrd,
                                      const int& pichrt,
                                      const int& pdyr );
     void initElmntBatch( void );
       
     void initializeMITTEMregion( const int& pdm );
          
     void setGridParameters( void );
    
     void setInTEMState( void );
    
     void setMITCLMgrid( const int& outyr,
                         const int& pdm, 
                         const long& igrd,
                         const int& ichrt );

     void setMITgridO3( const int& pdm,
                        const long& igrd );

     void setMITLULCgrid( const int& pdyr,
                          const long& igrd );

     void setOutTEMState( void );
    
     void setRunMode( void );          
    
     void setRunTime( void );
    
     void setTEMPred( void );
    
     long skipRecords( void );
    
     void starttem( void );
    
     #ifdef STANDALONE_TEM 
       void updateLCLUCregion( const int& pdyr,
                               FILE* fnumchrts,
                               FILE* flandarea,
                               FILE* flulc );
     #endif

     #ifdef STANDALONE_TEM 
       void updateMITCLMregion( const int& pdyr,
                                const int& pdm,
                                ifstream& ifnirr, 
                                ifstream& iftair,
                                ifstream& ifdayTsoilL1,
                                ifstream& ifdayTsoilL2,
                                ifstream& ifdayTsoilL3,
                                ifstream& ifdayTsoilL4,
                                ifstream& ifdayTsoilL5,
                                ifstream& ifdayTsoilL6,
                                ifstream& ifdayTsoilL7,
                                ifstream& ifdayTsoilL8,
                                ifstream& ifdayTsoilL9,
                                ifstream& ifdayTair, 
                                ifstream& ifrainDur,
                                ifstream& ifrainInt,
                                ifstream& ifdaySML1,
                                ifstream& ifdaySML2,
                                ifstream& ifdaySML3,
                                ifstream& ifdaySML4,
                                ifstream& ifdaySML5,
                                ifstream& ifdaySML6,
                                ifstream& ifdaySML7,
                                ifstream& ifdaySML8,
                                ifstream& ifdaySML9,
                                ifstream& ifhrSML1,
                                ifstream& ifhrSML2,
                                ifstream& ifhrSML3,
                                ifstream& ifhrSML4,
                                ifstream& ifhrSML5,
                                ifstream& ifhrSML6,
                                ifstream& ifprec,
                                ifstream& ifeet,
                                ifstream& ifsh2o1,
                                ifstream& ifsh2o2,
                                ifstream& ifspack,
                                ifstream& ifsrfrun,
                                ifstream& ifdrain,
                                ifstream& ifco2,
                                ifstream& ifo3 );
     #endif
  
     void updateMITTEMregion( const int& pdyr, 
                              const int& pdm );
    

// *************************************************************

     Elmnt43 elmnt;

     int end1;

     int equil;

     int fatalerr;

     FILE* felev;

     ofstream flog1;

     FILE* flonlat;
  
     FILE* fslayer;
     
     FILE* fstxt;

     ofstream ftempred[NUMTEM];

     int glob_count;

     int icount;

     ifstream ifstate;  // Use TEMstate from a specified year

     int istateflag;

     int istateyear;

     string ilandareafname;

     long mxnumgrid;
   
     int numspin;

     ofstream ofstate;  // Save TEMstate for a specified year
     int ostateflag;
     int ostateyear;

     int RTIME;

     int spinflag;
     int spinoutfg;
     int spinoutyrs;
     int spintime;

     MITelmnt44 telmnt[MXMITNLAT];

     int temflag;

     vector<string> tempredmap;

     int totsptime;

     int transtime;

     int yrostateflag;


     #ifdef STANDALONE_TEM 
       // Hourly soil moisture across layers:
       //   0 - 0 to 1.75 cm
       //   1 - 1.75 to 4.50 cm
       //   2 - 4.50 to 9.00 cm
       //   3 - 9.00 to 16.5 cm
       //   4 - 16.5 to 29.0 cm
       //   5 - 29.0 to 49.4 cm 
       MITHrClmData44 thrSoilMoist[NLAYERS][MAXMDAYS][MAXDAYHRS];
     #endif
};

#endif

