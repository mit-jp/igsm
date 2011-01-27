/***************************************************************
****************************************************************
METHANOGEN51.H - Describes methanogen dynamics in 
                   terrestrial ecosystems

Created by Q. Zhuang, 19/Feb/2003

Modifications:

20060105 - DWK changed class CH4DMPRO to be 
           class Methanogen51
20060105 - DWK changed MethanePR() to produceCH4()  
20060105 - DWK changed char ecd[80] to be string ecd in
           getecdch4pro() and getch4pro() function calls
20060105 - DWK changed public functions EffectOM(),  
           EffectOMD(), EffectPH(), EffectRX(), and 
           EffectST() to private functions
20060105 - DWK added include methanogen51.cpp to bottom
           of file
20060105 - DWK changed getecdch4pro() to be getecd()
20060105 - DWK deleted getch4pro()
20060207 - DWK deleted public double oxi_c[MAXCMNT]
20060207 - DWK added inheritance of ProcessXML50
                       
****************************************************************
************************************************************* */

#ifndef METHANOGEN51_H
#define METHANOGEN51_H

#include "mdmconsts51.hpp"

// Methanogen51 inherits ProcessXML43 Class
#include "tprocessXML437.h" 

class Methanogen51 : public ProcessXML43
{

   public:

     Methanogen51();
     
/* *************************************************************
			Public Functions
************************************************************* */
    
     void    getecd( ofstream& rflog1 );

     void    getecd( const string& ecd );

     double  produceCH4( const int& pcmnt,
                         const double& vegNPP,  
                         const double& depthz, 
                         const double& rootd, 
                         const double& soilt,   
                         const double& soilph,  
                         const double& ehl );
     
     // "Get" and "Set" private variables and parameters    

     // kc *****************************************************

     inline double getKC( const int& pcmnt ) 
     { 
       return kc[pcmnt]; 
     };

     inline void setKC( const double& pkc, const int& pcmnt ) 
     { 
       kc[pcmnt] = pkc; 
     };


     // lowb ***************************************************

     inline double getLOWB( const int& pcmnt ) 
     { 
       return lowb[pcmnt]; 
     };

     inline void setLOWB( const double& plowb, const int& pcmnt ) 
     { 
       lowb[pcmnt] = plowb; 
     };


     // maxfresh ***********************************************

     inline double getMAXFRESH( const int& pcmnt ) 
     { 
       return maxfresh[pcmnt]; 
     };

     inline void setMAXFRESH( const double& pmaxfresh, 
                              const int& pcmnt ) 
     { 
       maxfresh[pcmnt] = pmaxfresh; 
     };


     // mgo ****************************************************

     inline double getMGO( const int& pcmnt ) 
     { 
       return mgo[pcmnt]; 
     };

     inline void setMGO( const double& pmgo, const int& pcmnt ) 
     { 
       mgo[pcmnt] = pmgo; 
     };


     // pmethaneq **********************************************

     inline double getPMETHANEQ( const int& pcmnt ) 
     { 
       return pmethaneq[pcmnt]; 
     };

     inline void setPMETHANEQ( const double& ppmethaneq, 
                               const int& pcmnt ) 
     { 
       pmethaneq[pcmnt] = ppmethaneq; 
     };


     // proref *************************************************

     inline double getPROREF( const int& pcmnt ) 
     { 
       return proref[pcmnt]; 
     };

     inline void setPROREF( const double& pproref, 
                            const int& pcmnt ) 
     { 
       proref[pcmnt] = pproref; 
     };
     
     
   private:

/* *************************************************************
			Private Functions
************************************************************* */

     double EffectOM( const int& pcmnt,
                      const double& vegNPP );

     double EffectOMD( const int& pcmnt, 
                       const double& depthz, 
                       const double& rootd );
                      
     double EffectPH( const double& soilph );
    
     double EffectRX( const double& ehl );

     double EffectST( const int& pcmnt,
                      const double& soilt );


/* *************************************************************
			Private Variables and Parameters
************************************************************* */

     double annpp;

     double ehlt;

     // lowb -- lower boundary -- parameter say, 1000mm, 1m
     double lowb[MAXCMNT];

     // maxfresh is parameter for the site, annual max value 
     //   (gc m-2 yr-1)
     double maxfresh[MAXCMNT];

     // maximum methane production rate, uMh-1; 
     //   assume 0.3 umh-1 for Alaska tundra (Walter et al., )
     double mgo[MAXCMNT];
     
     //Kc -- half saturation constant of methanogenic acetate 
     //  consumption 0.1 mM, Segers and kengen, 1998
     double kc[MAXCMNT];
     
//     double oxi_c[MAXCMNT];

     double pmethaneq[MAXCMNT];

     double proref[MAXCMNT];

//     double suborg[CYCLE][MXMNDAYS];
     double suborg;
   
};

#endif

