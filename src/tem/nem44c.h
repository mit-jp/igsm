/* *************************************************************
NEM44a.H - MIT Natural Ecosystems Model based on the DNDC model
             of Li et al. (1992)

20060128 - Created by DWK by modifying nem2.h
20060128 - DWK added include nemconsts.hpp
20060128 - DWK deleted include nem2.cpp from bottom of file
20060202 - DWK changed public double 
           dayMoist[NLAYERS][CYCLE][MAXMDAYS] to
           double dayMoist[NLAYERS][MAXMDAYS]
20060202 - DWK changed public double 
           dayTemp[NLAYERS][CYCLE][MAXMDAYS] to
           double dayTemp[NLAYERS][MAXMDAYS]
20060202 - DWK changed public double 
           hourMoist[NLAYERS][MAXMDAYS][CYCLE][MAXDAYHRS] to
           double hourMoist[NLAYERS][MAXMDAYS][MAXDAYHRS]
20060202 - DWK deleted const int dm from function call to 
           setDailyCH4emissions()
20060302 - DWK changed public 
           double dayMoist[NLAYERS][MAXMDAYS] to
           double dayMoist[NLAYERS+3][MAXMDAYS]
20060302 - DWK changed public 
           double dayTemp[NLAYERS][MAXMDAYS] to
           double dayTemp[NLAYERS+3][MAXMDAYS]
20060302 - DWK added public double soilyrthick[NLAYERS+3]
20060327 - DWK added public getecd()
20060327 - DWK added private double initDntrf[MAXCMNT]
20060327 - DWK added inheritance of ProcessXML43
20080130 - DWK changed include from nemconsts.hpp to 
           nemconsts544a.h
20080130 - DWK changed include from tprocessXML437.h to 
           tigsmprocessXML44a.h
20080130 - DWK changed ProcessXML43 to ProcessXML44
20110707 - DWK changed include from nemconsts44a.hpp to
           nemconsts44c.hpp
20110707 - DWK changed include from tigsmprocessXML44a.h to
           tigsmprocessXML44c.h
                                                         
************************************************************* */

/* *************************************************************
References:

Li, C., S. Frolking and T.A. Frolking.  1992.  A model of 
  nitrous oxide evolution from soil driven by rainfall events:
  1. Model structure and sensitivity.  Journal of Geophysical
  Research 97(D9): 9759-9776.

Liu, Y.  1996.  Modeling the Emissions of Nitrous Oxide (N2O) 
  and Methane (CH4) from the Terrestrial Biosphere to the 
  Atmosphere.  MIT Joint Program on the Science and Policy of
  Global Change Report No. 10.  August 1996. 219pp.
  
************************************************************* */

#ifndef NEM_H
#define NEM_H

#include "nemconsts44c.hpp"
#include "tigsmprocessXML44c.h"

class NEM : public ProcessXML44
{
  public: 
  
     NEM();

/* *************************************************************
                       Public Functions
************************************************************* */

     void getecd( ofstream& rflog1 );

     void getecd( const string& ecd );


     double setMonthlyCH4emissions( const int& vegtype,
                                    const int& outmon,
                                    const double& tair,
                                    const double& prec,
                                    const double& pet,
                                    const float& lat );

     void stepmonth( const int& dm,
                     const int& IVEG,
                     const int& pcmnt,
                     const double& clay,
                     const double& totpor,
                     const double& density,
                     const double& pH,
                     const double& Ksat,
                     const double& CIN,
                     double AVTEMP[CLMNLAYERS][MAXMDAYS],
                     double AVWTC[CLMNLAYERS][MAXMDAYS],
                     double HRWTC[NLAYERS][MAXMDAYS][MAXDAYHRS],
                     double TEM[MAXMDAYS],
                     double RAINDUR[MAXMDAYS],
                     double RAININTE[MAXMDAYS] );
               
                     
     // "Get" and "Set" private variables and parameters   

     // ANH4IN *************************************************
          
     inline double getANH4IN( const int& dlvl ) 
     { 
       return ANH4IN[dlvl]; 
     }

     inline void setANH4IN( const double& panh4in, 
                            const int& dlvl )
     {
       ANH4IN[dlvl] = panh4in;  
     }

     
     // ANO3IN *************************************************
          
     inline double getANO3IN( const int& dlvl ) 
     { 
       return ANO3IN[dlvl]; 
     }

     inline void setANO3IN( const double& pano3in, 
                            const int& dlvl )
     {
       ANO3IN[dlvl] = pano3in;  
     }


     // dayMoist ***********************************************
     
//     inline double getDAYMOIST( const int& dlyr, 
//                                const int& dday ) 
//     { 
//       return dayMoist[dlyr][dday]; 
//     }

//     inline void setDAYMOIST( const double& pdaymoist,
//                              const int& dlyr, 
//                              const int& dday ) 
//     { 
//       dayMoist[dlyr][dday] = pdaymoist; 
//     }


     // dayTemp ************************************************
     
//     inline double getDAYTEMP( const int& dlyr, 
//                               const int& dday ) 
//     { 
//       return dayTemp[dlyr][dday]; 
//     }

//     inline void setDAYTEMP( const double& pdaytemp,
//                             const int& dlyr, 
//                             const int& dday ) 
//     { 
//       dayTemp[dlyr][dday] = pdaytemp; 
//     }

     // DPHUMIN ************************************************
          
     inline double getDPHUMIN( const int& dlvl ) 
     { 
       return DPHUMIN[dlvl]; 
     }

     inline void setDPHUMIN( const double& pdphumin, 
                            const int& dlvl )
     {
       DPHUMIN[dlvl] = pdphumin;  
     }

     // ECH4N ************************************************
     
     inline double getECH4N( void ) { return ECH4N; }

     inline void setECH4N( const double& pech4n ) 
     { 
       ECH4N = pech4n; 
     }


     // ECO2DN *************************************************
     
     inline double getECO2DN( void ) { return ECO2DN; }

     inline void setECO2DN( const double& peco2dn ) 
     { 
       ECO2DN = peco2dn; 
     }


     // ECO2N *************************************************
     
     inline double getECO2N( void ) { return ECO2N; }

     inline void setECO2N( const double& peco2n ) 
     { 
       ECO2N = peco2n; 
     }


     // EN2DN *************************************************
     
     inline double getEN2DN( void ) { return EN2DN; }

     inline void setEN2DN( const double& pen2dn ) 
     { 
       EN2DN = pen2dn; 
     }


     // EN2ODN *************************************************
     
     inline double getEN2ODN( void ) { return EN2ODN; }

     inline void setEN2ODN( const double& pen2odn ) 
     { 
       EN2ODN = pen2odn; 
     }


     // EN2ON *************************************************
     
     inline double getEN2ON( void ) { return EN2ON; }

     inline void setEN2ON( const double& pen2on ) 
     { 
       EN2ON = pen2on; 
     }


     // hourMoist ***********************************************
     
//     inline double getHOURMOIST( const int& dlyr, 
//                                 const int& dday,
//                                 const int& dhr ) 
//     { 
//       return hourMoist[dlyr][dday][dhr]; 
//     }

//     inline void setHOURMOIST( const double& phrmoist,
//                               const int& dlyr, 
//                               const int& dday,
//                               const int& dhr ) 
//     { 
//       hourMoist[dlyr][dday][dhr] = phrmoist; 
//     }


     // initDntrf **********************************************
     
     inline double getINITDNTRF( const int& pcmnt ) 
     { 
       return initDntrf[pcmnt]; 
     }

     inline void setINITDNTRF( const double& pidntrf, 
                               const int& pcmnt ) 
     { 
       initDntrf[pcmnt] = pidntrf; 
     }


     // nonreactiveSoilC ***************************************
     
     inline double getNONSOLC( void ) 
     { 
       return nonreactiveSoilC; 
     }

     inline void setNONSOLC( const double& pnsolc ) 
     { 
       nonreactiveSoilC = pnsolc; 
     }


     // OCIN ***************************************************
          
     inline double getOCIN( const int& dlvl ) 
     { 
       return OCIN[dlvl]; 
     }

     inline void setOCIN( const double& pocin, 
                            const int& dlvl )
     {
       OCIN[dlvl] = pocin;  
     }


     // pH******************************************************
     
     inline double getPH( void ) 
     { 
       return pH; 
     }

     inline void setPH( const double& pph ) 
     { 
       pH = pph; 
     }


     // RCLIN **************************************************
          
     inline double getRCLIN( const int& dlvl ) 
     { 
       return RCLIN[dlvl]; 
     }

     inline void setRCLIN( const double& prclin, 
                            const int& dlvl )
     {
       RCLIN[dlvl] = prclin;  
     }


     // RCRIN **************************************************
          
     inline double getRCRIN( const int& dlvl ) 
     { 
       return RCRIN[dlvl]; 
     }

     inline void setRCRIN( const double& prcrin, 
                            const int& dlvl )
     {
       RCRIN[dlvl] = prcrin;  
     }


     // RCVLIN *************************************************
          
     inline double getRCVLIN( const int& dlvl ) 
     { 
       return RCVLIN[dlvl]; 
     }

     inline void setRCVLIN( const double& prcvlin, 
                            const int& dlvl )
     {
       RCVLIN[dlvl] = prcvlin;  
     }


     // topdens ************************************************
     
     inline double getTOPDENS( void ) { return topdens; }

     inline void setTOPDENS( const double& ptopdens ) 
     { 
       topdens = ptopdens; 
     }


     // topKsat ************************************************
     
     inline double getTOPKSAT( void ) { return topKsat; }

     inline void setTOPKSAT( const double& ptopksat ) 
     { 
       topKsat = ptopksat; 
     }


     // toppor *************************************************
     
     inline double getTOPPOR( void ) { return toppor; }

     inline void setTOPPOR( const double& ptoppor ) 
     { 
       toppor = ptoppor; 
     }

     
/* *************************************************************
                       Public Variables
************************************************************* */
     
     int dbugflg;

     double daych4flx[MAXMDAYS];

     // Daily soil moisture across layers:
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
     double dayMoist[CLMNLAYERS][MAXMDAYS];

     // Daily soil temperatures across layers:
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
     double dayTemp[CLMNLAYERS][MAXMDAYS];     

     // Flag to indicate first simulated day
     bool firstday;

     // Hourly soil moisture across layers:
     //   0 -   0    to 1.75 cm
     //   1 -   1.75 to 4.50 cm
     //   2 -   4.50 to 9.00 cm
     //   3 -   9.00 to 16.5 cm
     //   4 -  16.5  to 29.0 cm
     //   5 -  29.0  to 49.4 cm 
     double hourMoist[NLAYERS][MAXMDAYS][MAXDAYHRS];

     
     //Flag to indicate when to start running NEM
     int startflag;


  private:     

/* *************************************************************
                       Private Functions
************************************************************* */

     double adsorbAmmonium( const double& clay,
                            const double& oldNH4,
                            const double& ocini );

     void CARBONTR( const double& carbon, 
                    const int& vegSoilDistrib, 
                    const double& layerThick, 
                    double& topr, 
                    int& toph, 
                    double& chdp, 
                    double& ttoo, 
                    double& ratio );

     double ch4xprec( const double& prec, const double& pet );

     double ch4xtair( const double& tair );

     void decomposeHumads( const double& clayc, 
                           const double& rfm,
                           const double& ddrf );

     void decomposeMicrobialBiomass( const double& rfm,
                                     const double& ddrf );

     void decomposeResidues( const double& rfm,
                             const double& ddrf );

     void denitrification( const int& pcmnt,
                           const int& dlayer,
                           const double& wfps );

     void makeAmmonia( const int& dryday,
                       const int& dlayer,
                       const double& oldNH4 );

     void nitrification( const int& dlayer,
                         const double& atmp, 
                         const double& moist );

     void resetInitialMonthlyConditions( void );
     
     double setDailyCH4emissions( const float& lat,
                                  const int& dday, 
                                  const double& dayH2OtableZ );

     void stepday( const int& DY, 
                   const int& DAY,
                   const double& clay,
                   const double& CLAYC,
                   const double& toppor,
                   double AVTEMP[CLMNLAYERS][MAXMDAYS],
                   double AVWTC[CLMNLAYERS][MAXMDAYS] );

     void stephour( const int& dm,
                    int& dndy,
                    const int& dhr,
                    const int& pcmnt,
                    const int& rainrun,
                    const double& clay,
                    const double& toppor,
                    double HRWTC[NLAYERS][MAXMDAYS][MAXDAYHRS],
                    const double& tt,
                    const double& ttt,
                    const double& lf,
                    const double& tm,
                    const int& maxlayers );

     void uptakeMethane( const int& dy,
                         const double& wrc,
                         const double& wtoc,
                         const double& whmus, 
                         double AVTEMP[CLMNLAYERS][MAXMDAYS] );

/* *************************************************************
                       Private Variables
************************************************************* */

     // Leached NO3 at bottom of profile in a rain event
     //   (KgN/ha/rain event)
     double AAAA;

     double AEN2O[NLVL];

     double ANH4IN[NLVL];
     double ANO3IN[NLVL];

     double AVERMOIS[NLVL];

     // New denitrifier biomass in soil layers at time T
     double B[NLVL+1];

     double BACN[NLVL+1];

     double BACNN;

     // Leached NO2 at bottom of profile in a rain event
     //   (KgN/ha/rain event)
     double BBBB;

     double BX;

     // New soluble C content in soil layers at time T
     double C[NLVL+1];

     double CH4N[MAXMDAYS];
     double CH4UP;

     // Adsorption factor of organic C by clay.
     double CLAYC;

     double CNO3;

     double CO2DN[MAXMDAYS];
     double CO2N[MAXMDAYS];

     // Residue C in KgC/ha/layer.
     double CR[NLVL];

     double CRB[NLVL];

     // Labile bio C,KgC/ha/layer.
     double CRB1[NLVL];
  
     // Carbon in Labile Microbial Biomass in KgC/ha/layer
     double CRB1X;

     double CRB2[NLVL];

    // Carbon in Resistant Microbial Biomass in KgC/ha/layer
     double CRB2X;

     // Totl microbial biomass C.
     double CRBX;

     double CRH[NLVL];

     // Carbon in Labile Humads in KgC/ha/layer
     double CRHL;
  
     // Carbon in Resistant Humads in KgC/ha/layer
     double CRHR;

     // Total carbon in humads (labile + resistant) 
     //   in KgC/ha/layer
     double CRHX;

     double CX;

     // Carbon from decomposed microbial biomass that becomes 
     //   available soluble C
     double DCBAVAI;

     // Carbon from decomposed humads that becomes 
     //   available soluble C
     double DCHAVAI;

     double DCO2[NLVL+1];

     double DCO2T;

     // Carbon from decomposed humads that is recycled 
     //   to micrbial biomass  
     double DHBC;

     // Nitrogen from decomposed humads that is recycled 
     //   to micrbial biomass  
     double DHBN;

     double DHUM[NLVL+1];

     // Denitrifier population in saturated soil layers
     double DPB[NLVL+1];

     // Humus C in KgC/ha/layer.
     double DPHUM[NLVL];

     double DPHUMIN[NLVL];
       
     // Carbon from decomposed humads that is transferred
     //   to humus
     double DPHUMX;

     // Organic C in saturated soil layers
     double DPOC[NLVL];

     // Soluble C in saturated soil layers
     double DPSC[NLVL+1];  

     // CO2 production from decomposition of residue
     double DRCA;

     // Carbon from decomposed residue that is 
     //   transferred to microbial biomass 
     double DRCB;

     // Nitrogen from decomposed residue that is 
     //   transferred to microbial biomass 
     double DRNB;

     // Time step (1 hour) 
     double DT;

     double DYN2EMS;
     
     double DYN2OEMS;

     double ECH4N;

     double ECO2DN;

     double ECO2N;

     // Accumulated N2O emitted in a rain event
     double EE;

     double EMSN2[NLVL+1];
     
     double EMSN2O[NLVL+1];

     double EN2DN;

     double EN2ODN;

     double EN2ON;

     // CO2 production from decomposition of microbial biomass
     double FBCO2;

     // Carbon from decomposed microbial biomass that is 
     //   transferred to humads 
     double FBHC;

     // Nitrogen from decomposed microbial biomass that is 
     //   transferred to humads 
     double FBHN;
     
     double FCO2[NLVL];

     // Accumulated N2 emitted in a rain event
     double FF;

     // CO2 production from decomposition of humads
     double FHCO2;

     // Accumulated N consumed by denitri synthesis
     double FIN;

     int firstn2oyr;
     
     // Total N consumption by denitrifiers
     double FIXN[NLVL+1];

     // Definition of thickness of a layer in m, based on rain 
     //   water infiltration speed; 2 is calibrated factor for 
     //   soil percolate.
     double H;

     // Accumulated humus in profile,KgC/ha
     double HUM;

     // Accumulated CO2 emissions in profile
     double II;

     // Water-filled pore space when denitrification begins to occur
     double initDntrf[MAXCMNT];

     // Half-saturation value of soluble carbon (Kg C /layer/ha)
     double KC;

     // Half-saturation value of N-oxides (Kg N /layer/ha)
     double KN;

     double LENO2[NLVL+1];
     
     double LENO3[NLVL+1];

     // CH4 uptake by soil for a yr,KgC/ha.
     double LICH4UP;

     // Mass of soil in a layer for a ha,Kg/ha.
     double M;

     // Number of days per month
     int mdays[CYCLE];
     
     // Cumulative number of Julian days by the end of the month
     int mjdays[CYCLE];
     
     // Total available nitrogen
     double N;

     // N2 content in soil layers at time T
     double N2[NLVL+1];

     double N2DN[MAXMDAYS];

     // N2O content in soil layers at time T
     double N2O[NLVL+1];
     
     double N2ODN[MAXMDAYS];
    
     double N2ON[MAXMDAYS];

     double N2OX;

     double N2X;

     double NH3[NLVL];

       // Initial NH4 in soil after addition of N-fertilizers
     double NH40;

      // New NH4 pool after residue decomposition
     double NH41;

      // New NH4 pool after microbial biomass decomposition
     double NH43;

     double NH44;

     // New NH4 pool after NH3 loss
     double NH443;

     // New NH4 pool after N2O loss
     double NH45;

     // New NH4+ after nitrification
     double NH46;

     // Initial NH4+ in KgN/ha/layer.
     double NH4INI[NLVL];

     // Nitrification rate in KgN/ha/layer
     double NITRIFY[NLVL+1];

     double NO2[NLVL];

     // NO2 in soil layers at a certain hour
     double NO2P[NLVL+1];
   
     double NO2X;

     double NO3[NLVL];

     // Initial NO3 in soil after addition of N-fertilizers
     double NO30;

      // New NO3 pool after residue decomposition
     double NO31;

     // New NO3 pool after N2O loss (i.e. no change)
     double NO35;

     // New NO3 pool after nitrification
     double NO36;

     // Initial NO3- in KgN/ha/layer.
     double NO3INI[NLVL];

     // NO3 in soil layers at a certain hour
     double NO3P[NLVL+1];

     // NO3 in soil layers at a certain hour
     double NO3X;

     // "Non-reactive soil organic carbon in the soil column
     double nonreactiveSoilC;

     // Microbial biomass N after biomass decomposition.
     double NRB;

     // Total nitrogen in humads 
     double NRH;

     // Nitrogen in labile humads.
     double NRHL;
  
     // Nitrogen in resistant humads.
     double NRHR;

     double OCIN[NLVL];

     // Active organic C(bio+humads),KgC/ha/lyr.
     double OCINI[NLVL];

     // Soil pH
     double pH;

     // Effect of pH on NO3 denitrifiers.
     double PHK1;

     // Effect of pHh on N2O denitrifiers.
     double PHK2;

     // Number of soil layer.
     int Q;

     // Labile residue C in KgC/ha/layer.
     double RCL[NLVL];

     double RCLIN[NLVL];
    
     // Carbon in labile residue in Kg C/ha/layer.
     double RCLX;
  
     // Resistant residue C,KgC/ha/layer.
     double RCR[NLVL];

     double RCRIN[NLVL];
  
     // Carbon in resistant residue in Kg C/ha/layer.
     double RCRX;

     // Very labile residue Cin Kg C/ha/layer.
     double RCVL[NLVL];

     double RCVLIN[NLVL];

     // Carbon in very Labile residue in Kg C/ha/layer.
     double RCVLX;

     double SNH4[NLVL];

     double TCAVAI[NLVL];

     // Effect of temp on denitrifiers
     double TE;

     // Parameter controling denitrifier death & growth.
     double TEM1;

     double TMPCO2DN[NLVL];
     double TMPN2DN[NLVL];
     double TMPN2ODN[NLVL];

     double TOC[NLVL];

     double TON[NLVL];

     // Effective density of the top 30 cm of soil
     double topdens;

     // Effect saturated hydraulic conductivity of the top 30 cm
     //   of soil
     double topKsat;

     // Effective porosity of the top 30 cm of soil
     double toppor;

     double TT;

     // Accumulated N uptake for the dry period,KgN/ha.
     double TUPTN;
  
     // Accumulated NH3 volatilization for the dry period,
     //   Kg N / ha.
     double TVNH3;  
     
     double TZW;

     // Air-occupied volume in a layer in m^3/ha.
     double V;

     // Ammonia volatilization rate in 1/day.
     double VOLNH3;

     // Water volume in layer L,m^3/ha.
     double W;

     // Accumulated nitri N2O emission for the dry period,KgN/ha
     double WAEN2O;

     // Accumulated nitri CO2 emission for the dry period,KgC/ha
     double WFCO2; 

     double WNITRIFY;

     double ZB[NLVL+1];

     double ZC[NLVL+1];

     double ZN2[NLVL+1];

     double ZN2O[NLVL+1];

/* *************************************************************
                       Parameters
************************************************************* */

     // Soil water parameter (range from 4 to 11.4), no unit
     double BETA[NST+1];

     //   TOPR and TOPH are parameters for the description of soil organic matter
     //   profile destribution. Agricultural land is different from natural land.
     //   Setup the profile of organic matter in the soil which is related to 
     //   TOPR and TOPH.
     int TOPH;

};

#endif
