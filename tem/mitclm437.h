/* *************************************************************
****************************************************************
MITCLM437.H - physical characteristics of the atmosphere as
                described by the MIT 2D L-O GCM

Modifications:

20060128 - DWK created by modifying mitclm436a.h 
20060128 - DWK added include mittemconsts43.hpp
20060128 - DWK changed include from tclm436.h to tclm437.h
20060128 - DWK changed include from mitdata436.h to 
           mitdata437.h
20060128 - DWK changed include from mitclmdata436.h to 
           mitclmdata437.h
20060128 - DWK changed include from mitdaychrtdata436.h to
           mitdaychrtdata437.h
20060128 - DWK changed include from mitdayclmdata436.h to
           mitdayclmdata437.h
20060128 - DWK changed include from mithrclmdata436.h to
           mithrclmdata437.h
20060128 - DWK changed include from tclmdat436.h to
           tclmdat437.h
20060128 - DWK changed include from o3data436.h to o3data437.h
20060201 - DWK renamed public function mkMITd40() to be
           mkMITaot40()
20060217 - DWK added public string idaySMendL7, 
           string idaySMfnameL7,
           string idaySMendL8, string idaySMfnameL8,
           string idaySMendL9, string idaySMfnameL9,
           string idayTsoilendL7, string idayTsoilfnameL7,
           string idayTsoilendL8, string idayTsoilfnameL8,
           string idayTsoilendL9, and string idayTsoilfnameL9
20060320 - DWK added public string ispackfname, 
           string ispackend,
           string isrfrunfname, string isrfrunend,
           string idrainfname, string idrainend
20060320 - DWK added public initDRAINAGE(), initSNOWPACK() and
           initSURFRUN()
20060420 - DWK added compiler directive STANDALONE_TEM
                      
****************************************************************
************************************************************* */

#ifndef MITCLM437_H
#define MITCLM437_H


// NOTE: If running TEM "offline", make sure the compiler 
//         directive STANDALONE_TEM is DEFINED below.  If  
//         coupling TEM to the IGSM, make sure STANDALONE_TEM is 
//         NOT DEFINED (i.e. comment out next line)

#include "preproc.h"

// Get global constants

#include "mittemconsts43.hpp"

#include "tclm437.h"  // MITclm43 inherits TEMclm43 class


// MITclm43 uses the MITdata43 class
#include "mitdata437.h" 

#ifdef STANDALONE_TEM 
  // MITclm43 uses the MITCLMdata43 class
  #include "mitclmdata437.h" 
#endif

#ifdef STANDALONE_TEM 
  // MITclm43 uses the MITDayChrtData43 class
  #include "mitdaychrtdata437.h" 
#endif

#ifdef STANDALONE_TEM  
  // MITclm43 uses the MITDayClmData43 class
  #include "mitdayclmdata437.h" 
#endif

#ifdef STANDALONE_TEM 
  // MITclm43 uses the MITHrClmData43 class
  #include "mithrclmdata437.h" 
#endif

#ifdef STANDALONE_TEM 
  // MITclm43 uses the Clmdata43 class
  #include "tclmdat437.h" 
#endif

#ifdef STANDALONE_TEM 
  // MITclm43 uses the O3data class
  #include "o3data437.h"   
#endif

class MITclm43 : public TEMclm43 
{

  public:

     MITclm43();

     enum flagkey { OFF, ON };

     enum deltakey { NONE, RATIO, DIFF };

/* *************************************************************
		 Public Functions
************************************************************* */

     void aggregFlux2D( const int& year,
                        const float& lat, 
                        const double& subarea, 
                        const double& influx,
                        MITdata43& outflux );


     #ifdef STANDALONE_TEM 
       void initDaySoilMoist( ofstream& rflog1 );
     #endif

     #ifdef STANDALONE_TEM 
       void initDayTair( ofstream& rflog1 );
     #endif
   
     #ifdef STANDALONE_TEM 
       void initDayTsoil( ofstream& rflog1 );
     #endif
   
     #ifdef STANDALONE_TEM 
       void initDRAINAGE( ofstream& rflog1 );
     #endif

     #ifdef STANDALONE_TEM 
       void initEET( ofstream& rflog1 );
     #endif

     #ifdef STANDALONE_TEM 
       void initHourSoilMoist( ofstream& rflog1 );
     #endif

     void initMITCflux( ofstream& rflog1 );

     void initMITCH4flux( ofstream& rflog1 );

     void initMITN2Oflux( ofstream& rflog1 );

     #ifdef STANDALONE_TEM
       void initMITSolarRad( ofstream& rflog1 );
     #endif

     #ifdef STANDALONE_TEM 
       void initPET( ofstream& rflog1 );
     #endif

     #ifdef STANDALONE_TEM 
       void initRainDuration( ofstream& rflog1 );
     #endif

     #ifdef STANDALONE_TEM 
       void initRainIntensity( ofstream& rflog1 );
     #endif

     #ifdef STANDALONE_TEM 
       void initSNOWPACK( ofstream& rflog1 );
     #endif

     #ifdef STANDALONE_TEM 
       void initSoilH2O( ofstream& rflog1 );
     #endif

     #ifdef STANDALONE_TEM 
       void initSURFRUN( ofstream& rflog1 );
     #endif

     void mkMITaot40( const double ozone3hrs[MAX3HRS],
		      const int& pdm );

     void resetCFLUX( const int& outmon );

     void resetCH4FLUX( const int& outmon );

     void resetN2OFLUX( const int& outmon );

     #ifdef STANDALONE_TEM 
       void setLatBandClm( ifstream& ifile, MITCLMdata43& clm );
     #endif

     void setLatBandClm( ifstream& ifile, MITdata43& clm );

     #ifdef STANDALONE_TEM 
       void setLatBandClm( ifstream& ifile,
                            const int& outmon, 
                            MITDayChrtData43 clm[MAXMDAYS] );
     #endif

     #ifdef STANDALONE_TEM 
       void setLatBandClm( ifstream& ifile, 
                           const int& outmon,
                           MITDayClmData43 clm[MAXMDAYS] );
     #endif

     #ifdef STANDALONE_TEM 
       void setLatBando3( ifstream& ifile, O3data43& clm );
     #endif

     #ifdef STANDALONE_TEM 
       void setLatBandSoil( ifstream& ifile,
                            const int& outmon, 
                            MITDayChrtData43 soilenv[MAXMDAYS] );
     #endif

     #ifdef STANDALONE_TEM 
       void setLatBandSoil( ifstream& ifile,
                            const int& outmon, 
                            MITHrClmData43 soilenv[MAXMDAYS][MAXDAYHRS] );
     #endif

     void setMITCH4Flags( ofstream& rflog1 );

     void setMITN2OFlags( ofstream& rflog1 );

     void setMITSolarRadFlags( ofstream& rflog1,
                               const int& requil );

     void timecheck( ofstream& rflog1 );


/* *************************************************************
		 Public Variables
************************************************************* */

     int begin2end;
     int beginyr;

     int ch4flag;
     
     int cldsyr;
     int co2yr;

     int endyr;

     int firstyear;

     FILE* feet;
     FILE* fnirr;
     FILE* fo3;
     FILE* fpar;
     FILE* fpet;
     FILE* fprec;
     FILE* fsh2o1;
     FILE* fsh2o2;
     FILE* ftair;

     int girryr;

     string icflxend;
     string icflxfname;

     string ich4flxend;
     string ich4flxfname;

     #ifdef STANDALONE_TEM 
       string idayTairend;
       string idayTairfname;
     #endif

     #ifdef STANDALONE_TEM 
       string idaySMendL1;
       string idaySMfnameL1;
     #endif

     #ifdef STANDALONE_TEM 
       string idaySMendL2;
       string idaySMfnameL2;
     #endif

     #ifdef STANDALONE_TEM  
       string idaySMendL3;
       string idaySMfnameL3;
     #endif

     #ifdef STANDALONE_TEM 
       string idaySMendL4;
       string idaySMfnameL4;
     #endif

     #ifdef STANDALONE_TEM 
       string idaySMendL5;
       string idaySMfnameL5;
     #endif

     #ifdef STANDALONE_TEM 
       string idaySMendL6;
       string idaySMfnameL6;
     #endif

     #ifdef STANDALONE_TEM 
       string idaySMendL7;
       string idaySMfnameL7;
     #endif

     #ifdef STANDALONE_TEM 
       string idaySMendL8;
       string idaySMfnameL8;
     #endif

     #ifdef STANDALONE_TEM 
       string idaySMendL9;
       string idaySMfnameL9;
     #endif

     #ifdef STANDALONE_TEM 
       string idayTsoilendL1;
       string idayTsoilfnameL1;
     #endif

     #ifdef STANDALONE_TEM 
       string idayTsoilendL2;
       string idayTsoilfnameL2;
     #endif

     #ifdef STANDALONE_TEM 
       string idayTsoilendL3;
       string idayTsoilfnameL3;
     #endif

     #ifdef STANDALONE_TEM 
       string idayTsoilendL4;
       string idayTsoilfnameL4;
     #endif

     #ifdef STANDALONE_TEM 
       string idayTsoilendL5;
       string idayTsoilfnameL5;
     #endif

     #ifdef STANDALONE_TEM 
       string idayTsoilendL6;
       string idayTsoilfnameL6;
     #endif

     #ifdef STANDALONE_TEM 
       string idayTsoilendL7;
       string idayTsoilfnameL7;
     #endif

     #ifdef STANDALONE_TEM 
       string idayTsoilendL8;
       string idayTsoilfnameL8;
     #endif

     #ifdef STANDALONE_TEM 
       string idayTsoilendL9;
       string idayTsoilfnameL9;
     #endif

     #ifdef STANDALONE_TEM 
       string idrainend;
       string idrainfname;
     #endif

     #ifdef STANDALONE_TEM 
       string ieetend;
       string ieetfname;
     #endif

     #ifdef STANDALONE_TEM 
       string ihrSMendL1;
       string ihrSMfnameL1;
     #endif

     #ifdef STANDALONE_TEM 
       string ihrSMendL2;
       string ihrSMfnameL2;
     #endif

     #ifdef STANDALONE_TEM 
       string ihrSMendL3;
       string ihrSMfnameL3;
     #endif

     #ifdef STANDALONE_TEM 
       string ihrSMendL4;
       string ihrSMfnameL4;
     #endif

     #ifdef STANDALONE_TEM 
       string ihrSMendL5;
       string ihrSMfnameL5;
     #endif

     #ifdef STANDALONE_TEM 
       string ihrSMendL6;
       string ihrSMfnameL6;
     #endif

     string in2oflxend;
     string in2oflxfname;

     #ifdef STANDALONE_TEM 
       string inoxend;
       string inoxfname;
     #endif

     #ifdef STANDALONE_TEM 
       string ipetend;
       string ipetfname;
     #endif
   
     #ifdef STANDALONE_TEM 
       string irainDurend;
       string irainDurfname;
     #endif

     #ifdef STANDALONE_TEM 
       string irainIntend;
       string irainIntfname;
     #endif
   
     #ifdef STANDALONE_TEM 
       string ish2o1end;
       string ish2o1fname;
     #endif

     #ifdef STANDALONE_TEM 
       string ish2o2end;
       string ish2o2fname;
     #endif
   
     #ifdef STANDALONE_TEM 
       string ispackend;
       string ispackfname;
     #endif

     #ifdef STANDALONE_TEM 
       string isrfrunend;
       string isrfrunfname;
     #endif

     int lastyear;

//     double lato3[MAX3HRS];

     int mdays[CYCLE];     

     int midyr;

     int n2oflag;
     
     int nirryr;

     int o3yr;

     int paryr;
     int precyr;
     
     int tairyr;

     // holds data from 2D MIT-IGSM latitudinal bands

     MITdata43 tcflx2D;
     MITdata43 tch4flx2D;

     MITdata43 tn2oflx2D;
         
     int yrmn;


  private:

/* *************************************************************
		 Private Functions
************************************************************* */

     void mkHourlyOzone( const double ozone3hrs[MAX3HRS] );

/* *************************************************************
		 Private Variables
************************************************************* */

     // Daily mean air temperature (degrees C)
     double dayTair[MAXMDAYS];

     double final[MAXHRS];

     // An array of the number of 0.5 degree latitudinal bands 
     //   for each MIT 2-D L-O GCM latitudinal band
     int latnum[MXMITNLAT];
     
     double nox;
     
//     double o3;

     // Duration of current rain event (hrs)
     double rainDuration[MAXMDAYS];

     // Intensity of rain during current event (mm/hr)
     double rainIntensity[MAXMDAYS];
};

#endif
