/* **************************************************************
*****************************************************************
TIGSMATMS44A.H - object describes physical characteristics of the
	               atmosphere

Modifications:

20060126 - DWK created by modifying atms50b5.h
20060126 - DWK changed temconsts51.hpp to temconsts43.hpp
20060126 - DWK changed class Atmosphere50 to class Atmosphere43
20060126 - DWK added public inline double getAOT40()
20060126 - DWK added private double aot40
20080130 - DWK changed include from temconsts43.hpp to 
           temigsmconsts44a.hpp
20080130 - DWK changed class Atmosphere43 to class Atmosphere44
                                                                                                                                                      
*****************************************************************
************************************************************** */

// Class uses global constants CYCLE and MAXRTIME

#ifndef ATMS44A_H
#define ATMS44A_H

#include "temigsmconsts44a.hpp"


class Atmosphere44
{

  public:

     Atmosphere44();

/* *************************************************************
		 Public Functions
************************************************************* */

     /* Determine Potential Evapotranspiration based on 
        algorithms from Jensen M. E. and H. R. Haise (1963) 
        Estimating evapotranspiration from solar radiation.  
        Journal of the Irrigation and Drainage Division 4: 
        14-41.  */

     void petjh( const double& nirr,
                 const double& tair,
                 const int& pdm );

     // Determine in precipitation occurs as rain or snow

     void precsplt( const double& prec,
                    const double& tair,
                    double& rain,
                    double& snowfall );

     void resetMonthlyFluxes( void );
     
     void resetYrFluxes( void );


     // "Get" and "Set" private variables

     // aot40 **************************************************
     
     inline double getAOT40( void ) { return aot40; }

     inline void setAOT40( const double& paot40 ) 
     { 
       aot40 = paot40; 
     }


     // avetair ************************************************
     
     inline double getAVETAIR( void ) { return avetair; }

     inline void setAVETAIR( const double& pavetair ) 
     { 
       avetair = pavetair; 
     }

     // clds ***************************************************
     
     inline double getCLDS( void ) { return clds; }

     inline void setCLDS( const double& pclds ) 
     { 
       clds = pclds; 
     }


     // co2 ****************************************************
     
     inline double getCO2( void ) { return co2; }

     inline void setCO2( const double& pco2 ) { co2 = pco2; }

     
     // co2level ***********************************************
     
     inline double getCO2LEVEL( void ) { return co2level; }

     inline void setCO2LEVEL( const double& pco2level ) 
     { 
       co2level = pco2level; 
     }

 
     // girr ***************************************************
     
     inline double getGIRR( void ) { return girr; }

     inline void setGIRR( const double& pgirr ) 
     { 
       girr = pgirr; 
     }


     // initco2 ************************************************
     
     inline double getINITCO2( void ) { return initco2; }

     inline void setINITCO2( const double& pinitco2 ) 
     { 
       initco2 = pinitco2; 
     }


     // mxtair *************************************************
     
     inline double getMXTAIR( void ) { return mxtair; }

     inline void setMXTAIR( const double& pmxtair ) 
     { 
       mxtair = pmxtair; 
     }


     // ndays[] ************************************************
     
     inline double getNDAYS( const int& pdm ) 
     { 
       return ndays[pdm]; 
     }


     // nirr ***************************************************
     
     inline double getNIRR( void ) { return nirr; }

     inline void setNIRR( const double& pnirr ) 
     { 
       nirr = pnirr; 
     }


     // par ****************************************************
     
     inline double getPAR( void ) { return par; }

     inline void setPAR( const double& ppar ) { par = ppar; }


     // pet ****************************************************
     
     inline double getPET( void ) { return pet; }

     inline void setPET( const double& ppet ) { pet = ppet; }


     // petmx **************************************************
     
     inline double getPETMX( void ) { return petmx; }

     inline void setPETMX( const double& ppetmx ) 
     { 
       petmx = ppetmx; 
     }


     // prec ***************************************************
     
     inline double getPREC( void ) { return prec; }

     inline void setPREC( const double& pprec ) 
     { 
       prec = pprec; 
     }


     // prev2tair **********************************************
     
     inline double getPREV2TAIR( void ) { return prev2tair; }

     inline void setPREV2TAIR( const double& pprev2tair ) 
     { 
       prev2tair = pprev2tair; 
     }


     // prevco2 ************************************************
     
     inline double getPREVCO2( void ) { return prevco2; }

     inline void setPREVCO2( const double& pprevco2 ) 
     { 
       prevco2 = pprevco2; 
     }


     // prevtair ***********************************************
     
     inline double getPREVTAIR( void ) { return prevtair; }

     inline void setPREVTAIR( const double& pprevtair ) 
     { 
       prevtair = pprevtair; 
     }


     // prvpetmx ***********************************************
     
     inline double getPRVPETMX( void ) { return prvpetmx; }

     inline void setPRVPETMX( const double& pprvpetmx ) 
     { 
       prvpetmx = pprvpetmx; 
     }


     // rain ***************************************************
     
     inline double getRAIN( void ) { return rain; }

     inline void setRAIN( const double& prain ) 
     { 
       rain = prain; 
     }


     // snowfall ***********************************************
     
     inline double getSNOWFALL( void ) { return snowfall; }

     inline void setSNOWFALL( const double& psnowfall ) 
     { 
       snowfall = psnowfall; 
     }


     // tair ***************************************************
     
     inline double getTAIR( void ) { return tair; }
        
     inline void setTAIR( const double& ptair ) 
     { 
       tair = ptair; 
     }

  



/* **************************************************************
		 Public Variables
************************************************************** */

     // Number of days per month     
     double ndays[CYCLE];   

     // Annual potential evapotranspiration (mm / year)
     double yrpet;

     // Annual total precipitation (mm / year)
     double yrprec;

     // Annual sum of rainfall (mm / year)
     double yrrain;

     // Annual snow ( mm / year)    
     double yrsnowfall;         
  
   
   private:
   
/* **************************************************************
		 Private Variables
************************************************************** */

     // Atmospheric ozone AOT 40 index (i.e. ozone 
     //   concentrations above 40 ppb-hr)
     double aot40;

     // Mean annual air temperature (degrees C)
     double avetair;

     // Cloudiness (%)    
     double clds;   
     
     // Atmospheric carbon dioxide concentration (ppmv)
     double co2; 

     // Constant atmospheric CO2 concentration (ppmv) used to 
     //   calibrate TEM or equilibrate TEM at beginning of 
     //   extrapolation     
     double co2level;

     // Gross Irradiance (cal/(sq. cm * day))
     double girr;

     // initial CO2 concentration (ppmv)     
     double initco2;     

     // Maximum monthly air temperature (degrees C)
     double mxtair;

     // Net Irradiance   (cal/(sq. cm * day))
     double nirr;

     // Photosynthetically Active Radiation  (cal/(sq.cm * day))
     double par;

     // Monthly potential evapotranspiration (mm / month)
     double pet;

     // Maximum PET of current year (mm / month)
     double petmx;

     // Total Precipitation (mm / month)
     double  prec;
     
     // Previous Month's Atmospheric CO2 Concentration (ppmv)
     double prevco2;

     // Previous Month's Air Temperature (degrees C)
     double prevtair;

     // Previous 2 Month's Air Temperature (degrees C)
     double prev2tair;

     // Maximum PET of previous year (mm / month)
     double prvpetmx;

     // Rainfall (mm / month)     
     double  rain;       

     // Snowfall (mm / month)     
     double  snowfall;   

     // Surface Air Temperature (degrees C)
     double tair;

};

#endif
