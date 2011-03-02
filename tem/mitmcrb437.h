/* **************************************************************
*****************************************************************
MITMCRB437.H - object describing characteristics of soil 
                microbes used in the Terrestrial Ecosystem Model 
                (TEM) within the MIT IGSM

Modifications:

20060128 - DWK created by modifying tmcrb437.h
20060128 - DWK changed include from temconsts43.hpp to
           nemconsts.hpp
20060128 - DWK added include nem437.h
20060128 - DWK added public NEM nem
20060301 - DWK added include mdm51.h
20060301 - DWK added public MDM51 mdm
                                
*****************************************************************
************************************************************** */


#ifndef MITTMCRB437_H
#define MITTMCRB437_H

// Tmicrobe43 uses the global constants NUMMSAC and MAXCMNT
#include "nemconsts.hpp"

#include "tprocessXML437.h"

#include "bioms423.hpp"

#include "nem437.h"

 #include "mdm51.h"

class Tmicrobe43: public ProcessXML43
{

   public:

     Tmicrobe43();

/* **************************************************************
		 Public Functions
************************************************************** */

     void getvegecd( ofstream& rflog1 );

     void getvegecd( const string& ecd );

     void resetEcds( const int& pcmnt, const double& psiplusc );
      
     void resetMonthlyFluxes( void );

     void resetYrFluxes( void );

     void setDQ10( const int& pdcmnt, 
                   const double& tair, 
                   const double& snowpack );

     double setRHMOIST( const int& pdcmnt,
                        const double& pcfldcap,
                        const double& vsm,
                        const int& moistlim );

     void showecd( const int& pdcmnt );

     void updateDynamics( const int& pcmnt,
                          const double& pcfldcap,
                          const double& soilorgc,
                          const double& soilorgn,
                          const double& soilh2o,
                          const double& vsm,
                          const double& availn,
                          const int& moistlim,
                          const int& tillflag,
                          const double& tillfactor,
                          const double& ksoil );
     
     double yrkd( const int& nfeed,
                  const double& yrltrc,
                  const double& yrltrn,
                  const int& pdcmnt );


     // "Get" and "Set" private variables and parameters
     
     // cnsoil *************************************************
     
     inline double getCNSOIL( const int& pcmnt ) 
     { 
       return cnsoil[pcmnt]; 
     }

     inline void setCNSOIL( const double& pcnsoil, 
                            const int& pcmnt ) 
     { 
       cnsoil[pcmnt] = pcnsoil; 
     }


     // decay **************************************************
     
     inline double getDECAY( void ) { return decay; }


     // kd *****************************************************
     
     inline double getKD( void ) { return kd; }

     inline void setKD( const double& pkd ) { kd = pkd; }


     // kda ****************************************************
     
     inline double getKDA( const int& pcmnt ) 
     { 
       return kda[pcmnt]; 
     }

     inline void setKDA( const double& pkda, 
                         const int& pcmnt ) 
     { 
       kda[pcmnt] = pkda; 
     }


     // kdb ****************************************************
     
     inline double getKDB( const int& pcmnt ) 
     { 
       return kdb[pcmnt]; 
     }

     inline void setKDB( const double& pkdb, 
                         const int& pcmnt ) 
     { 
       kdb[pcmnt] = pkdb; 
     }


     // kdc ****************************************************
     
     inline double getKDC( void ) { return kdc; }

     inline void setKDC( const double& pkdc ) { kdc = pkdc; }


     // kn2 ****************************************************
     
     inline double getKN2( const int& pcmnt ) 
     { 
       return kn2[pcmnt]; 
     }

     inline void setKN2( const double& pkn2, 
                         const int& pcmnt ) 
     { 
       kn2[pcmnt] = pkn2; 
     }


     // lcclnc *************************************************
     
     inline double getLCCLNC( const int& pcmnt ) 
     { 
       return lcclnc[pcmnt]; 
     }

     inline void setLCCLNC( const double& plcclnc, 
                            const int& pcmnt ) 
     { 
       lcclnc[pcmnt] = plcclnc; 
     }


     // moistmax ***********************************************
     
     inline double getMOISTMAX( const int& pcmnt ) 
     { 
       return moistmax[pcmnt]; 
     }

     inline void setMOISTMAX( const double& pmoistmax, 
                              const int& pcmnt ) 
     { 
       moistmax[pcmnt] = pmoistmax; 
     }


     // moistmin ***********************************************
     
     inline double getMOISTMIN( const int& pcmnt ) 
     { 
       return moistmin[pcmnt]; 
     }

     inline void setMOISTMIN( const double& pmoistmin, 
                              const int& pcmnt ) 
     { 
       moistmin[pcmnt] = pmoistmin; 
     }


     // moistopt ***********************************************
     
     inline double getMOISTOPT( const int& pcmnt ) 
     { 
       return moistopt[pcmnt]; 
     }

     inline void setMOISTOPT( const double& pmoistopt, 
                              const int& pcmnt ) 
     { 
       moistopt[pcmnt] = pmoistopt; 
     }


     // netnmin ************************************************
     
     inline double getNETNMIN( void ) { return netnmin; }

     inline void setNETNMIN( const double& pnetnmin ) 
     { 
       netnmin = pnetnmin; 
     }


     // nfixpar ************************************************
     
     inline double getNFIXPAR( const int& pcmnt ) 
     { 
       return nfixpar[pcmnt]; 
     }

     inline void setNFIXPAR( const double& pnfixpar, 
                             const int& pcmnt ) 
     { 
       nfixpar[pcmnt] = pnfixpar; 
     }


     // nup ****************************************************
     
     inline double getNUP( void ) { return nup; }

     inline void setNUP( const double& pnup ) 
     { 
       nup = pnup; 
     }


     // nupa ***********************************************
     
     inline double getNUPA( const int& pcmnt ) 
     { 
       return nupa[pcmnt]; 
     }

     inline void setNUPA( const double& pnupa, 
                          const int& pcmnt ) 
     { 
       nupa[pcmnt] = pnupa; 
     }


     // nupb ***************************************************
     
     inline double getNUPB( const int& pcmnt ) 
     { 
       return nupb[pcmnt]; 
     }

     inline void setNUPB( const double& pnupb, 
                          const int& pcmnt ) 
     { 
       nupb[pcmnt] = pnupb; 
     }


     // nuptake ************************************************
     
     inline double getNUPTAKE( void ) { return nuptake; }


     // propftos ***********************************************
     
     inline double getPROPFTOS( const int& pcmnt ) 
     { 
       return propftos[pcmnt]; 
     }

     inline void setPROPFTOS( const double& ppropftos, 
                              const int& pcmnt ) 
     { 
       propftos[pcmnt] = ppropftos; 
     }


     // rh *****************************************************
     
     inline double getRH( void ) { return rh; }
    

     // rhq10 **************************************************
     
     inline double getRHQ10( const int& pcmnt ) 
     { 
       return rhq10[pcmnt]; 
     }

     inline void setRHQ10( const double& prhq10, 
                           const int& pcmnt ) 
     { 
       rhq10[pcmnt] = prhq10; 
     }


/* **************************************************************
		 Public Variables
************************************************************** */

     MDM51 mdm;
     
     NEM nem;

     // Annual sum of netnmin
     double yrnmin;         // (g N / (sq. meter * year))

     // Annual sum of nuptake
     double yrnuptake;      // (g N / (sq. meter * year))
     
     // Annual sum of rh
     double yrrh;       // (g C / (sq. meter * year))



   private:

/* **************************************************************
		 Private Functions
************************************************************** */

     double nminxclm( const int& pdcmnt,
                      const double& soilh2o,
                      const double& soilorgc,
                      const double& soilorgn,
                      const double& availn,
                      const double& ksoil );

     double rhxclm( const double& soilorgc,
                    const double& dq10,
                    const double& moist );

   
/* **************************************************************
		 Private Variables
************************************************************** */
     
     // Effect of temperature on decomposition
     double dq10;

     // Net nitrogen mineralization
     double netnmin; // (g N / (sq. meter * month))

     // Total nitrogen uptake or "immobilzation" by microbes
     double nuptake;  // (g N / (sq. meter * month))

     // Heterotrophic respiration
     double rh;  // (g C / (sq. meter * month))


/* *************************************************************
		 Private Parameters
************************************************************* */

     double cnsoil[MAXCMNT];

     // Parameter representing the quality of soil organic matter

     double decay;

     // Biome-specific decomposition parameters for function rhxclm

     double kd;
     double kda[MAXCMNT];
     double kdb[MAXCMNT];
     double kdc;

     double kdin[NUMMSAC];     // kd values read in from file
     double kdsave[NUMMSAC];   // kd values saved to a file

     // Biome-specific half saturation parameter for function
     //   nminxclm describing the effect of available nitrogen
     //   on microbial nitrogen uptake

     double kn2[MAXCMNT];


     double lcclnc[MAXCMNT];

     // Biome-specific parameters describing the influence of 
     //   soil moisture on decomposition (i.e., moist)

     double moistmin[MAXCMNT];
     double moistopt[MAXCMNT];
     double moistmax[MAXCMNT];

    // N fixation parameter

     double nfixpar[MAXCMNT];


     double nup;
     double nupa[MAXCMNT];
     double nupb[MAXCMNT];
     
     double propftos[MAXCMNT];

     double rhq10[MAXCMNT];

};

#endif
