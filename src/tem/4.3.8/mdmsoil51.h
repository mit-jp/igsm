/* *************************************************************
****************************************************************
MDMSOIL51.H - object describing characteristics of soil used by
	      the Methane Dynamics Model (MDM)

Modifications:

20060125 - DWK created by combining and modifying 
           methanediff.hpp, methaneplant.hpp and methaneebu.hpp 
20060125 - DWK added public function void oxygenC() from 
           methaneoxi.hpp
20060125 - DWK added public function double RedoxP() from
           methanepro.hpp 
20060125 - DWK added private double oxygenc[200] from 
           methaneoxi.hpp
20060207 - DWK changed all public variables to be private 
           variables     
20060207 - DWK added private double ph
20060422 - DWK renamed public function interpolateSM() to
           interpolateLayers() and deleted public function
           interpolateST()
                                            
****************************************************************
************************************************************* */

#ifndef MDMSOIL51_H
#define MDMSOIL51_H

// Global constants for MDM
#include "mdmconsts51.hpp"

// Methanogen51 inherits ProcessXML43 Class
#include "tprocessXML437.h" 

#define min( a, b ) ((a) <= (b) ? (a) : (b))
#define max( a, b ) ((a) >= (b) ? (a) : (b))


class MDMsoil51 : public ProcessXML43
{

  public:
    
     MDMsoil51();

/* *************************************************************
		 Public Functions
************************************************************* */


     // CH4ebullition() formerly CH4ebuR()
     
     double CH4ebullition( const double& ch4con );


     // CH4plantTransport() formerly CH4lpR()
     
     double CH4plantTransport( const double& tveg, 
                               const double& prootd,
                               const double& pz, 
                               const double& psoilt,
                               const double& ch4con );

     // dryCH4diffusion() formerly CH4flxandCon()

     double  dryCH4diffusion( const double& ddf, 
                              const double& lowb );

     // dryCH4uptake() formerly ACH4flxandCon()

     double  dryCH4uptake( const double& ddf, 
                           const double& lowb );

     void    getecd( ofstream& rflog1 );

     void    getecd( const string& ecd );

     double interpolateLayers( const double layerVar[10], 
                               const double layerThick[10],
                               const int& nlayers,
                               const double& depth );

     void setOxygenConc( const double& afp, 
                         const double& lowb );

     double RedoxP( const double& watertable, 
                    const double& depthz, 
                    const double& wfpsl,
                    const int& pcmnt, 
                    const double& lowb );

     double  setDiffusivityFactor( const double& sand, 
                                   const double& silt, 
                                   const double& clay,  
                                   int satyorn );

     void setUnsaturatedZoneMoisture( const double& ph2otable );
       
     // wetCH4diffusion() formerly BCH4flxandCon()

     double  wetCH4diffusion( const double& lowb, 
                              const double& watertable );

     // wetCH4flux() formerly CCH4flxandCon()

     double  wetCH4flux( const double& lowb, 
                         const double& watertable, 
                         const double& soilt );


     // "Get" and "Set" private variables and parameters    

     // ch4con *************************************************
     
     inline double getCH4CON( const int& dlyr ) 
     { 
       return ch4con[dlyr]; };

     inline void setCH4CON( const double& pch4con, 
                            const int& dlyr ) 
     { 
       ch4con[dlyr] = pch4con; 
     }


     // ch4con_b ***********************************************
     
     inline double getCH4CON_B( const int& dlyr ) 
     { 
       return ch4con_b[dlyr]; };

     inline void setCH4CON_B( const double& pch4con_b, 
                              const int& dlyr ) 
     { 
       ch4con_b[dlyr] = pch4con_b; 
     }


     // ch4cons **************************************************
     
     inline double getCH4CONS( void ) { return ch4cons; };

     inline void setCH4CONS( const double& pch4cons ) 
     { 
       ch4cons = pch4cons; 
     }


     // ch4ebuflx ***********************************************
     
     inline double getCH4EBUFLX( const int& dlyr ) 
     { 
       return ch4ebuflx[dlyr]; };

     inline void setCH4EBUFLX( const double& pch4ebuflx, 
                               const int& dlyr ) 
     { 
       ch4ebuflx[dlyr] = pch4ebuflx; 
     }


     // ch4emis **************************************************
     
     inline double getCH4EMIS( void ) { return ch4emis; };

     inline void setCH4EMIS( const double& pch4emis ) 
     { 
       ch4emis = pch4emis; 
     }


     // ch4flx **************************************************
     
     inline double getCH4FLX( void ) { return ch4flx; };

     inline void setCH4FLX( const double& pch4flx ) 
     { 
       ch4flx = pch4flx; 
     }


     // ch4oxirate1 ********************************************
     
     inline double getCH4OXIRATE1( const int& dlyr ) 
     { 
       return ch4oxirate1[dlyr]; };

     inline void setCH4OXIRATE1( const double& pch4oxirate1, 
                                 const int& dlyr ) 
     { 
       ch4oxirate1[dlyr] = pch4oxirate1; 
     }


     // ch4plantr ********************************************
     
     inline double getCH4PLANTR( const int& dlyr ) 
     { 
       return ch4plantr[dlyr]; };

     inline void setCH4PLANTR( const double& pch4plantr, 
                               const int& dlyr ) 
     { 
       ch4plantr[dlyr] = pch4plantr; 
     }


     // ch4ratesat ******************************************
     
     inline double getCH4RATESAT( const int& dlyr ) 
     { 
       return ch4ratesat[dlyr]; };

     inline void setCH4RATESAT( const double& pch4ratesat, 
                                const int& dlyr ) 
     { 
       ch4ratesat[dlyr] = pch4ratesat; 
     }


     // ch4pltflx **************************************************
     
     inline double getCH4PLTFLX( void ) { return ch4pltflx; };

     inline void setCH4PLTFLX( const double& pch4pltflx ) 
     { 
       ch4pltflx = pch4pltflx; 
     }


     // ch4rate ************************************************
     
     inline double getCH4RATE( const int& dlyr ) 
     { 
       return ch4rate[dlyr]; };

     inline void setCH4RATE( const double& pch4rate, 
                             const int& dlyr ) 
     { 
       ch4rate[dlyr] = pch4rate; 
     }


     // ch4tot **************************************************
     
     inline double getCH4TOT( void ) { return ch4tot; };

     inline void setCH4TOT( const double& pch4tot ) 
     { 
       ch4tot = pch4tot; 
     }


     // ebuflx *************************************************
     
     inline double getEBUFLX( void ) { return ebuflx; };

     inline void setEBUFLX( const double& pebuflx ) 
     { 
       ebuflx = pebuflx; 
     }


     // ehlt *************************************************
     
     inline double getEHLT( void ) { return ehlt; };

     inline void setEHLT( const double& pehlt ) 
     { 
       ehlt = pehlt; 
     }


     // fdf ****************************************************
     
     inline double getFDF( void ) { return fdf; };

     inline void setFDF( const double& pfdf ) 
     { 
       fdf = pfdf; 
     }


     // fdfs ***************************************************
     
     inline double getFDFS( const int& dlyr ) 
     { 
       return fdfs[dlyr]; };

     inline void setFDFS( const double& pfdfs, 
                          const int& dlyr ) 
     { 
       fdfs[dlyr] = pfdfs; 
     }


     // intersoilh2o *******************************************
     
     inline double getINTERSH2O( const int& dlyr ) 
     { 
       return intersoilh2o[dlyr]; };

     inline void setINTERSH2O( const double& pintersoilh2o, 
                               const int& dlyr ) 
     { 
       intersoilh2o[dlyr] = pintersoilh2o; 
     }


     // intersoilt *********************************************
     
     inline double getINTERSOILT( const int& dlyr ) 
     { 
       return intersoilt[dlyr]; };

     inline void setINTERSOILT( const double& pintersoilt, 
                                const int& dlyr ) 
     { 
       intersoilt[dlyr] = pintersoilt; 
     }


     // oxygenc *********************************************
     
     inline double getOXYGENC( const int& dlyr ) 
     { 
       return oxygenc[dlyr]; };

     inline void setOXYGENC( const double& poxygenc, 
                                const int& dlyr ) 
     { 
       oxygenc[dlyr] = poxygenc; 
     }


     // PA *****************************************************

     inline double getPA( const int& pcmnt ) 
     { 
       return PA[pcmnt]; 
     };

     inline void setPA( const double& pPA, const int& pcmnt ) 
     { 
       PA[pcmnt] = pPA; 
     };


     // sat ****************************************************
     
     inline int getSAT( void ) { return sat; };

     inline void setSAT( const int& psat ) 
     { 
       sat = psat; 
     }

     // tsoil **************************************************
     
     inline double getTSOIL( void ) { return tsoil; };

     inline void setTSOIL( const double& ptsoil ) 
     { 
       tsoil = ptsoil; 
     }


     // tveg ****************************************************
     
     inline double getTVEG( void ) { return tveg; };

     inline void setTVEG( const double& ptveg ) 
     { 
       tveg = ptveg; 
     }


     // unsatthetaWL *******************************************
     
     inline double getUNSATTHETAWL( const int& dlyr ) 
     { 
       return unsatthetaWL[dlyr]; };

     inline void setUNSATTHETAWL( const double& punsatthetaWL, 
                                  const int& dlyr ) 
     { 
       unsatthetaWL[dlyr] = punsatthetaWL; 
     }

    
   private:
   
/* *************************************************************
		 Private Functions
************************************************************* */

     double fgrow( const double& soilt );

     double froot( const double& rootd, const double& z );



/* **************************************************************
		 Private Variables and Parameters
************************************************************** */
     
     double ch4con[MXMDMNLAY]; // CH4 concetration for diffusion
     
     double ch4con_b[MXMDMNLAY]; // CH4 concetration for diffusion

     // Daily efflux of methane via three transport mechanisms
//     double ch4cons[CYCLE][MXMNDAYS]; 
     double ch4cons; 

     // Ebullition flux of methane from each 1-cm soil layer
     double ch4ebuflx[MXMDMNLAY];

     // Daily efflux of methane via three transport mechanisms
//     double ch4emis[CYCLE][MXMNDAYS]; 
     double ch4emis; 

     double ch4flx;
      
     double ch4oxirate1[MXMDMNLAY];

     double ch4plantr[MXMDMNLAY];

//     double ch4pltflx[CYCLE][MXMNDAYS];
     double ch4pltflx;

     // nodes to store CH4
     double ch4rate[MXMDMNLAY];  

     // for deposit saturation zone ch4 concentration
     double ch4ratesat[MXMDMNLAY];  
  
//     double ch4tot[CYCLE][MXMNDAYS];
     double ch4tot;

     // Daily ebullition flux of methane from soil surface
//     double ebuflx[CYCLE][MXMNDAYS];
     double ebuflx;
    
     double ehlt;

     double fdf;

     double fdfs[MXMDMNLAY];

     // Interpolated soil moistures
     double intersoilh2o[MXMDMNLAY];

     // Interpolated soil temperatures
     double intersoilt[MXMDMNLAY];

     double ndays[CYCLE];

     // to store the oxygen concentration (every 10 mm)
     double oxygenc[MXMDMNLAY]; 

    // PA between 0 to 1.0:
    //   grasses and sedges are good transport = 1.0,
    //   tree = 0.5. 
    //   mosses = 0.0,

     double PA[MAXCMNT];

     // Percent clay in soil texture
     double pctclay;        

     //field capacity (%soil volume)
     double pcfldcap;         

     // Percent sand in soil texture
     double pctsand;        
     
     // Percent silt in soil texture
     double pctsilt;

      // pH value
     double ph;              

     //effective rooting depth (m)
     double rootz;            

     // Saturation index 
     //   sat = 0 for unsaturated soils 
     //   sat = 1 for saturated soils
     int sat;


     // Mean soil temperature (degrees C )of top 20 cm of soil 
     //   profile
     double tsoil;

     double tveg;
     
     double unsatthetaWL[MXMDMNLAY];
     
};

#endif
