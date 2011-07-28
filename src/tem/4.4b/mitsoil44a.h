/* **************************************************************
*****************************************************************
MITSOIL44A.H - object describing characteristics of soil used by
     	       the Terrestrial Ecosystem Model (TEM)

Modifications:

20070128 - DWK created by modifying mitsoil437.h 
20070128 - DWK changed include tsoil437.h to tsoil44.h
20070128 - DWK changed class MITsoil43 to class MITsoil44
20070128 - DWK changed inheritance from Tsoil43 to Tsoil44
20080130 - DWK changed include from mittemconsts43.hpp to
           mittemconsts44a.hpp
20080130 - DWK changed include from tsoil44.h to tigsmsoil44a.h
                                 
*****************************************************************
************************************************************** */

#ifndef MITSOIL44_H
#define MITSOIL44_H

#include "mittemconsts44a.hpp"

#include "tigsmsoil44a.h"

class MITsoil44 : public Tsoil44
{

  public:

     MITsoil44( void );

/* **************************************************************
		 Public Functions
************************************************************** */

     void resetMonthlyTraceGasFluxes( void );

     void resetYrTraceGasFluxes( void );
     
     double setDailyCH4emissions( const float& lat, 
                                  const int& dday, 
                                  const double& dayH2OtableZ );

     double setMonthlyCH4emissions( const int& vegtype,
                                    const int& outmon,
                                    const double& tair,
                                    const double& prec,
                                    const double& pet,
                                    const float& lat );
     

     // "Get" and "Set" private variables and parameters   
     
     // ch4consump *********************************************
     
     inline double getCH4CONSUMP( void ) { return ch4consump; }

     inline void setCH4CONSUMP( const double& pch4consump ) 
     { 
       ch4consump = pch4consump; 
     }
     

     // ch4emissions *******************************************
     
     inline double getCH4EMISS( void ) { return ch4emissions; }

     inline void setCH4EMISS( const double& pch4emissions ) 
     { 
       ch4emissions = pch4emissions; 
     }


     // ch4flux ************************************************
     
     inline double getCH4FLUX( void ) { return ch4flux; }

     inline void setCH4FLUX( const double& pch4flux ) 
     { 
       ch4flux = pch4flux; 
     }


     // co2dnflux **********************************************
     
     inline double getCO2DNFLUX( void ) { return co2dnflux; }

     inline void setCO2DNFLUX( const double& pco2dnflux ) 
     { 
       co2dnflux = pco2dnflux; 
     }


     // co2nflux ***********************************************
     
     inline double getCO2NFLUX( void ) { return co2nflux; }

     inline void setCO2NFLUX( const double& pco2nflux ) 
     { 
       co2nflux = pco2nflux; 
     }


     // n2flux *************************************************
     
     inline double getN2FLUX( void ) { return n2flux; }

     inline void setN2FLUX( const double& pn2flux ) 
     {
       n2flux = pn2flux; 
     }


     // n2odnflux **********************************************
     
     inline double getN2ODNFLUX( void ) { return n2odnflux; }

     inline void setN2ODNFLUX( const double& pn2odnflux ) 
     { 
       n2odnflux = pn2odnflux; 
     }


     // n2oflux ************************************************
     
     inline double getN2OFLUX( void ) { return n2oflux; }

     inline void setN2OFLUX( const double& pn2oflux ) 
     { 
       n2oflux = pn2oflux; 
     }


     // n2onflux ***********************************************
     
     inline double getN2ONFLUX( void ) { return n2onflux; }

     inline void setN2ONFLUX( const double& pn2onflux ) 
     { 
       n2onflux = pn2onflux; 
     }


     // pH *****************************************************
     
     inline double getPH( void ) { return pH; }

     inline void setPH( const double& pph ) { pH = pph; }


/* **************************************************************
		 Public Variables
************************************************************** */
          
     double daych4flx[MAXMDAYS];
          
     //Soil density
     double density[CLMNLAYERS];
     
     // Saturated hydraulic conductivity
     double Ksat[CLMNLAYERS];
     
     // Thickness of soil layers
     //   0 -   0    to 1.75 cm
     //   1 -   1.75 to 4.50 cm
     //   2 -   4.50 to 9.00 cm
     //   3 -   9.00 to 16.5 cm
     //   4 -  16.5  to 29.0 cm
     //   5 -  29.0  to 49.4 cm 
     //   6 -  49.4  to 82.9 cm
     //   7 -  83.0  to 138.3 cm
     //   8 - 138.4  to 229.6 cm
     //   9 - 229.7  to 343.3 cm
     double layerThick[CLMNLAYERS];

     // mean annual volumetric soil moisture (%)
     double meanvsm;       
          
     // Soil porosity as fraction of total soil volume
     double porosity[CLMNLAYERS];


     double yrch4csmp;

     double yrch4ems;

     double yrch4flx;

     double yrco2dnflx;

     double yrco2nflx;

     double yrn2flx;

     double yrn2odnflx;

     double yrn2oflx;

     double yrn2onflx;
     
 
  private:
  
/* *************************************************************
		 Private Variables
************************************************************* */

     double ch4consump;

     double ch4emissions;

     double ch4flux;

     double co2dnflux;    

     double co2nflux;

     int mdays[CYCLE];

     // N2 fluxes from soil (g N m^-2 mo^-1)
     double n2flux;

     // N2O fluxes from soil (g N m^-2 mo^-1)
     double n2odnflux;

     double n2oflux;

     double n2onflux;
    
     // Soil pH
     double pH;

};

#endif
