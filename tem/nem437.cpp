/* *************************************************************
NEM437.CPP - MIT Natural Emissions Model

20060128 - DWK created by modifying nem2.cpp
20060128 - DWK added include nem437.h and standard includes
20060327 - DWK added private double initDntrf[MAXCMNT] to
           denitrification()
20060327 - DWK added inheritance of ProcessXML43 to NEM()

************************************************************* */

#include<iostream>

  using std::cerr;
  using std::cin;
  using std::cout;
  using std::endl;
  using std::ios;

#include<fstream>

  using std::ifstream;
  using std::ofstream;

#include<cstdlib>

  using std::exit;
  using std::atof;
  using std::atoi;
  using std::atol;
  
#include<string>
  
  using std::string;

#include<cmath>

  using std::exp;
  using std::fabs;
  using std::log;
  using std::pow;

#include "nem437.h"


/* *************************************************************
************************************************************* */

NEM::NEM() : ProcessXML43()
{
  TZW = ZERO;

  // Number of days per month
  
  mdays[0] = 31;
  mdays[1] = 28;
  mdays[2] = 31;
  mdays[3] = 30;
  mdays[4] = 31;
  mdays[5] = 30;
  mdays[6] = 31;
  mdays[7] = 31;
  mdays[8] = 30;
  mdays[9] = 31;
  mdays[10] = 30;
  mdays[11] = 31;

  mjdays[0] = 31;
  mjdays[1] = 59;
  mjdays[2] = 90;
  mjdays[3] = 120;
  mjdays[4] = 151;
  mjdays[5] = 181;
  mjdays[6] = 212;
  mjdays[7] = 243;
  mjdays[8] = 273;
  mjdays[9] = 304;
  mjdays[10] = 334;
  mjdays[11] = 365;
  

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double NEM::adsorbAmmonium( const double& clay,
                            const double& oldNH4, 
                            const double& ocini )
{
  double ADSORB;
  double EXNH4;
  double newNH4;
  double NH441;
  double PEXNH4;
	
  NH441 = oldNH4; // formerly NH44
  
  NH441 *= (1000.0 / M);   // In gN/Kg soil, units shift.
  
  if( NH441 < 0.14 ) { PEXNH4 = 0.9; }
  else
  {
    PEXNH4 = 0.413 - 0.466*log( NH441 )/2.3026;
    
    if( PEXNH4 < ZERO ) { PEXNH4 = ZERO; }  // added on 3/7/95.
  }
        
  ADSORB = PEXNH4 * ((clay*0.01)/0.63) * (ocini/M/0.04);

  if ( ADSORB > 1.0 ) { ADSORB = 1.0; }
  
  EXNH4 = NH441 * ADSORB;   //Adsorbed NH4+, gN/Kg soil.
  
  newNH4 = NH441 - EXNH4;

//  EXNH4 = EXNH4 * M / 1000.0;
        
  // New NH4+ aftr adsorption in KgN/ha/layer.

  newNH4 *= M / 1000.0;  


  return newNH4;
   
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void NEM::CARBONTR( const double& carbon, 
                    const int& vegSoilDistrib, 
                    const double& layerThick, 
                    double& topr, 
                    int& toph, 
                    double& chdp, 
                    double& ttoo, 
                    double& ratio )
{
	
/* *************************************************************
Transfer 1m soil carbon content to carbon/soil ratio used in the 
  DNDCPT() only 
************************************************************* */

  double A;
  double AA;
  double BB;
  double CTTOO;
  int K;
  int NL;
  double SUM;
  double ttoo1;
	    
  NL = (int)(1.0 / layerThick);
      
  if( 1 == vegSoilDistrib )
  {
    A = 3.41378;
    topr = 2.38867;
    toph = (int)(0.20 / layerThick) + 1;
    chdp = 0.26537;
    SUM = ZERO;
    
    for( K = (toph+1); K <= NL; K++ )
    {
      SUM += (1.0/ pow( topr, (K*layerThick / chdp) ));
    }
    
    AA = carbon * A / 100.0 * layerThick * 100.0;
    
    BB = (NL-toph) / pow( topr, (toph*layerThick / chdp) ) - SUM;
    
    CTTOO = (carbon + AA*BB) / (double) NL;
    
    ttoo = CTTOO / (1.5 * 1000.0 * layerThick);
    
    ttoo1 = AA / (1.5 * 1000.0 * layerThick);
    
    ratio = ttoo1 / ttoo;
  }
  else if( 2 == vegSoilDistrib ) // 'agriculture'
  {
    ttoo = 0.018;
    topr = 2.0;                
    toph = (int) (0.15 / layerThick) + 1;
    chdp = 0.1;
    ratio = ZERO;
  }
  else
  {
//    A = 1.35226
    topr = 1.22568;     // 'grass' or 'trop'
    toph = 0;
    chdp = 0.35772;
//    CTTOO = carbon * A / 100.0
//            * pow( topr, (-H/CHDP) ) 
//            *layerThick * 100.0 !alternative method.
    SUM = ZERO;
    
    for( K = (toph+1); K <= NL; ++K )
    {
      SUM += (1.0/ pow( topr, ((K-toph)*layerThick/chdp) ));
    }
    
    CTTOO = carbon / ((double)(toph) + SUM);
    ttoo = CTTOO / (1.5 * 1000.0 * layerThick);
    ratio = 0.0;
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double NEM::ch4xprec( const double& prec, const double& pet )
{
  double pfactor;
  
  if( pet > ZERO ) 
  { 
    pfactor = prec / pet; 
  }
  else
  {
    pfactor = 1.000000;
  }
  
  if ( pfactor > 1.0 ) { pfactor = 1.000000; }
  
  return pfactor;  
      
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double NEM::ch4xtair( const double& tair )
{
  double tfactor;
  double blo;
  double bb;
  double bm;
  
  if( tair < -1.0 ) { tfactor = ZERO; }
  else
  {
    if( tair < 15.0 ) { bm = -8.0; } // -8.0 = 15.0 - 23.0
    else if( tair < 30.0 ) { bm = tair - 23.0; }
    else { bm = 7.0; }  // 7.0 = 30.0 - 23.0
    
    blo = exp( 0.334 * bm );  
    bb = blo / (1.0 + blo);
    tfactor = bb / 0.6271996;
  }
  
  return tfactor;  
      
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void NEM::decomposeHumads( const double& clayc, 
                           const double& rfm,
                           const double& ddrf )
{
  // If RHRL>0.0,part of resistant humads transfer 
  //   to labile.
  const double RHRL = 0.00;

  double DCHL;
  double DCHR;
  double DLHC;
  double DLHN;
  double DNHL;
  double DNHR;
  double FHNH4;


  // Total decomposition of humads
  double FHC;

  // Specific decomp rate of resist humads,1/day.
  double HRH;
         
   // Specific decomp rate of labile humads,1/day.
  double KRH;


  // Determine carbon transfer from resistant humads to 
  //   labile humads

  DLHC = CRHR * RHRL;  

        
  // Update carbon pools in resistant and labile humads 
        
  CRHR -= DLHC;  
  CRHL += DLHC;

  
  // Determine texture-dependent decomposition rates
  
  HRH = 0.006 * clayc;
  KRH = 0.16 * clayc;        


  // Determine change in resistant and labile humads
  //   carbon pools resulting from decomposition
  	       
  DCHR = rfm * HRH * ddrf * CRHR;
  DCHL = rfm * KRH * ddrf * CRHL;       
  FHC = DCHR + DCHL;  

 
  // Determine fate of decomposed humad carbon 
        
  DHBC = FHC * 0.2;   // recycled to microbial biomass 
  FHCO2 = FHC * 0.4;  // lost as CO2         
  DPHUMX = FHC * 0.4; // transferred to humus

       
  // Available C by decomp of humads.
        
  DCHAVAI = DHBC;     


  // Update humad carbon pools
  
  CRHR -= DCHR;
  CRHL -= DCHL;
  CRHX = CRHR + CRHL;

        
  // Determine nitrogen transfer from resistant humads to 
  //   labile humads
        
  DLHN = DLHC / RCNH;   

        
  // Update nitrogen pools in resistant and labile humads 

  NRHR -= DLHN;
  NRHL += DLHN;

  // Determine change in resistant and labile humads
  //   nitrogen pools resulting from decomposition

  DNHR = DCHR / RCNH;
  DNHL = DCHL / RCNH;

  // Determine fate of decomposed humad nitrogen
        
  DHBN = DHBC / RCNH;   // recycled to microbial biomass 
  FHNH4 = FHCO2 / RCNH; // mineralized to NH4+

  // Update humad nitrogen pools

  NRHR -= DNHR;
  NRHL -= DNHL;
        
  if( NRHR < ZERO ) { NRHR = ZERO; }
  if( NRHL < ZERO ) { NRHL = ZERO; }
        
  NRH = NRHR + NRHL;
        
  // New NH4+ pool after humads decomposition.
        
  NH44 = NH43 + FHNH4;  
  if( NH44 < ZERO ) { NH44 = ZERO; }

        
	
};
/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void NEM::decomposeMicrobialBiomass( const double& rfm,
                                     const double& ddrf )
{
  // Percentage of decomp biomass recycling to new biomass.
  const double EFFAC = 0.2;

  // Percentage of decomposed biomass transferred to humads.
  const double EFFNO = 0.2;

  // Specific rate of decomposition of resistant biomass in 
  //   1/day
  const double HRB = 0.04;
	
  // Specific rate of decomposition of labile biomass in 1/day
  const double KRB = 0.33;

  double DCRB;
  double DCRB1;
  double DCRB2;

  // FBBC is amount of decomposed microbial biomass that is 
  //   recycled to form new microbial biomass
  double FBBC;

  double FBNH4;

	
  // Determine change in resistant and labile microbial biomass
  //   carbon pools resulting from decomposition

  DCRB1 = rfm * KRB * ddrf * CRB1X;
  DCRB2 = rfm * HRB * ddrf * CRB2X;
  DCRB = DCRB1 + DCRB2;

  // Determine fate of decomposed microbial biomass carbon

  FBBC = DCRB * EFFAC;   // recycled to new microbial biomass
  FBHC = DCRB * EFFNO;   // transferred to humads.        
  FBCO2 = DCRB * (1.0 - EFFAC - EFFNO);  // lost as CO2 

  DCBAVAI = FBBC;  // available soluble carbon

  // Update microbial biomass carbon pools: CRB1X is for labile,
  //   CRB2X is for resistant.
        
  CRB1X += (FBBC*SRB - DCRB1);       
  CRB2X += (FBBC*(1.0-SRB) - DCRB2);
  CRBX = CRB1X + CRB2X;              

       
  // Determine fate of decomposed microbial biomass nitrogen
        
  FBHN = FBHC / RCNB;    // transferred to humads
  FBNH4 = FBCO2 / RCNB;  // mineralized to NH4+
        
  // Update microbial biomass nitrogen pool

  NRB = NRB - FBHN - FBNH4;
        
  // New NH4+ pool after bio decomposition
        
  NH43 = NH41 + FBNH4;  
	
};
/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void NEM::decomposeResidues( const double& rfm,
                             const double& ddrf )
{        
  // Ratio of Formed Biomass/CO2 produced in residue 
  //   decomposition
  const double EFFRB = 0.67;

  // Specific rate of CO2 production in labile residue.
  const double KRCL = 0.074;

  // Specific rate of CO2 production in resistant residue.
  const double KRCR = 0.02;

  // Specific rate of CO2 production in very labile residue.
  const double KRCVL = 0.074;

  // Ratio of C/N in labile residue.
  const double RCNRL = 20;
  
  // Ratio of C/N in resistant residue.
  const double RCNRR = 20;
  
  // Ratio of C/N in very labile residue.
  const double RCNRVL = 2.35;


  // The ratio of available C/available N. 
  double CCNN;

  double DRC;

  double DRCL; 
  double DRCLA;
  double DRCLB;

  double DRCR;
  double DRCRA;
  double DRCRB;

  double DRCVL;
  double DRCVLA;
  double DRCVLB;

  // Actual NH4+formed by residue decomposition (KgN/ha/layer)
  double DRNA;

  // The available N limiting factor on decomposition
  double MMNN;

  double PTAN;


  // Determine change in very labile, labile and resistant 
  //   residue carbon pools resulting from decomposition
  
  DRCVL = rfm * KRCVL * ddrf * RCVLX;  
  
  DRCL = rfm * KRCL * ddrf * RCLX;     
  
  DRCR = rfm * KRCR * ddrf * RCRX;    
  
  DRC = DRCVL + DRCL + DRCR;         

  
  PTAN = (DRCVL/RCNRVL + DRCL/RCNRL + DRCR/RCNRR)
         * (1.0 + EFFRB);
               
  if( ZERO == (PTAN+NO30+NH40) ) { MMNN = 0.2; }
  else
  {
    CCNN = DRC * (1.0+EFFRB) / (PTAN+NO30+NH40);
    
    if( CCNN < 9.0 ) { CCNN = 9.0; }
    
    MMNN = 0.2 + 7.2/CCNN;
  }

  // Determine fate of decomposed residue carbon
       
  DRCVLA = MMNN * DRCVL; // lost as CO2 from very labile residue       
  
  DRCLA = MMNN * DRCL;   // lost as CO2 from labil residue       
  
  DRCRA = MMNN * DRCR;   // lost as CO2 from reistant residue       
  
  DRCA = DRCVLA + DRCLA + DRCRA;  // total lost as CO2 

  // C Transferred to microbial biomass from very labile residue

  DRCVLB = DRCVLA * EFFRB; 

  // C Transferred to microbial biomass from labile residue

  DRCLB = DRCLA * EFFRB;

  // C Transferred to microbial biomass from resistant residue

  DRCRB = DRCRA * EFFRB;

  // Total C transferred to microbial biomass from residue

  DRCB = DRCVLB + DRCLB + DRCRB;

  // Update residue carbon pools
  
  RCVLX = RCVLX - DRCVLA - DRCVLB;
  
  RCLX = RCLX - DRCLA - DRCLB;
  
  RCRX = RCRX - DRCRA - DRCRB;
  
  if( RCVLX < ZERO ) { RCVLX = ZERO; }
  
  if( RCLX < ZERO ) { RCLX = ZERO; }
  
  if( RCRX < ZERO ) { RCRX = ZERO; }

  // Determine fate of decomposed residue nitrogen

  DRNA = DRCVLA * 1.67/RCNRVL + DRCLA * 1.67/RCNRL
         + DRCRA * 1.67/RCNRR; // mineralized to NH4+ 
        
  DRNB = DRCB / RCNB;  // transferred to microbial biomass

  // New NH4 and NO3 pools after residue decomposition
  
  if( (ZERO == NH40) && (ZERO == NO30) ) { NH40 = 1.0; }        
       
  NH41 = NH40 + DRNA - DRNB * (NH40/(NH40+NO30));  
  
  NO31 = NO30 - DRNB * (NO30/(NH40+NO30));
  
  if( NO31 < ZERO )
  {
    NO31 = ZERO;                       // New NO3- pool.
    NH41 = NH40 + DRNA - (DRNB-NO30);  // New NH4+ pool.
  }
	
};
/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void NEM::denitrification( const int& pcmnt,
                           const int& dlayer,
                           const double& wfps )
{
  // Maintenance coefficient on carbon (kg C/kg/h)
  const double MC = 0.0076;

  // Maintenance coefficient on N2O (kg N/kg/h)
  const double MN2O = 0.0792;
  
  // Maintenance coefficient on NO2 (kg N/kg/h)
  const double MNO2 = 0.0349;
  
  // Maintenance coefficient on NO3 (kg N/kg/h)
  const double MNO3 = 0.09;

  // N/C ratio of denitrifiers.
  const double R = 0.29;

  // Maximum growth rates of N2O denitrifier (percent/hour)
  const double UMN2O = 0.34;
  
  // Maximum growth rates of NO2 denitrifier (percent/hour)
  const double UMNO2 = 0.67;
  
  // Maximum growth rates of NO3 denitrifier (percent/hour)
  const double UMNO3 = 0.67;

  //   YMC, YMNO3, YMNO2, YMN2O----Maximum growth yield, in Kg denitrifierC
  //   /Kg soluble C, Kg NO3--N, KgNO2--N or KgN2O-N.
  // Maximum growth rate on soluble C (kg C/ kg C)
  const double YMC = 0.503;
  
  // Maximum growth rate on N2O (kg C/ kg N)
  const double YMN2O = 0.151;

  // Maximum growth rate on NO2 (kg C/ kg N)
  const double YMNO2 = 0.428;

  // Maximum growth rate on NO3 (kg C/ kg N)
  const double YMNO3 = 0.401;


  double CCK;

  // N2O consumed by denitrifiers
  double DAN2O;
  
  // NO2 consumed by denitrifiers
  double DANO2;
  
  // NO3 consumed by denitrifiers
  double DANO3;

  double DB;
  double DBD;
  double DBG;

  // Total consumption of soluble C
  double DC;

  double DCO2X;

  double DDN2;

  // Net increase in N2O
  double DDN2O;
  
  // Net increase in NO2
  double DDNO2;
  
  // Total amount of NO3 lost
  double DDNO3;

  double DHUMX;

  // N2O transferred to N2
  double DN2O;
  
  // NO2 transferred to N2O
  double DNO2;

  // NO3 converted to NO2 by nitrate denitrifiers
  double DNO3;

  double N2K;
  double N2OK;
  double N3K;

  //Total growth rate of denitrifiers, in percent/hr.
  double U;

  double UN2O;
  double UNO2;
  double UNO3;


  if( wfps > initDntrf[pcmnt] )
  {           
    //cout << "NEM: SoilMois = " << wfps << endl;
    //cout << "NEM: initDntrf = " << initDntrf[pcmnt] << endl;
    // Soluble C is shared by 3 denitrifiers.
            
    CCK = CX / (KC+CX);  
            
    N3K = NO3X / (KN+NO3X); 
    UNO3 = UMNO3 * CCK * N3K;   // Nitrate denitrifier
    N2K = NO2X / (KN+NO2X);
    UNO2 = UMNO2 * CCK * N2K;   // Nitrite denitrifier
    N2OK = N2OX / (KN+N2OX);
    UN2O = UMN2O * CCK * N2OK;  // N2O denitrifier                        
    U = UNO3 + UNO2 + UN2O;     
                                    
    DBG = U * BX * DT;  //Amount of denitrifier growth.
    DBD = MC * YMC * BX * DT; // Amount of death.
            
    if( TEM1 < ZERO ) { DBG = ZERO; }
    
    if( TEM1 < ZERO ) { DBD = ZERO; } 
    
    DB = DBG - DBD;  // Net increase of denitrifiers
    N = NO3X + NO2X + N2OX;  //Total available N.
             
    if( ZERO == N ) { N = 1.0; }
            
   // DNO3:NO3- converted to NO2- by nitrate denitrifier.
            
    DNO3 = (UNO3/YMNO3 + MNO3*NO3X/N) * PHK1 * BX * TE * DT;
                         
    if( ZERO == CX ) { DNO3 = ZERO; }
            
    if( ZERO == (NO3X+NO2X+N2OX) ) { DANO3 = ZERO; }
    else
    {
      // NO3- consumed by denitrifier synthesis.	
              
      DANO3 = R * DBG * (NO3X/(NO3X+NO2X+N2OX));                                
    }
            
    // Total amount of NO3- lost.
            
    DDNO3 = DNO3 + DANO3;   
            
    // New NO3- in layer L at hr T.
            
    NO3P[dlayer] = NO3X - DDNO3; 

    if( NO3P[dlayer] < ZERO ) { DNO3 = NO3X-DANO3; }
    
    if( NO3P[dlayer] < ZERO ) { NO3P[dlayer] = ZERO; }
            
    // NO2- transferred to N2O.
            
    DNO2 = (UNO2/YMNO2+MNO2*NO2X/N) * BX * TE * DT; 
                                                           
    if( ZERO == CX ) { DNO2 = ZERO; }
            
    if( ZERO == (NO3X+NO2X+N2OX) ) { DANO2 = ZERO; }
    else
    {
      // NO2- consumed by denitrifier synthesis.	
              
      DANO2 = R * DBG * (NO2X/(NO3X+NO2X+N2OX)); 
    }
            
    // Net increase in NO2-.
            
    DDNO2 = DNO3 - DNO2 - DANO2;  
            
    // New NO2- content in layer L at T.
            
    NO2P[dlayer] = NO2X + DDNO2;  
            
    if( NO2P[dlayer] < ZERO ) { DNO2 = DNO3 + NO2X - DANO2; }
    
    if( NO2P[dlayer] < ZERO ) { NO2P[dlayer] = ZERO; }
            
    // N2O transfred to N2.
            
    DN2O = (UN2O/YMN2O+MN2O*N2OX/N) * PHK2 * BX * TE * DT; 
                                                                
    if( ZERO == CX ) { DN2O = ZERO; }
            
    if( 0 == (NO3X+NO2X+N2OX) ) { DAN2O = ZERO; }
    else
    {
      // N2O consumed by denitrifier synthesis.	
              
      DAN2O = R * DBG * (N2OX/(NO3X+NO2X+N2OX));  
    }
            
    // Net increase in N2O.
            
    DDN2O = DNO2 - DN2O - DAN2O;
            
    // New N2O content in layer L at T.
              
    N2O[dlayer]= N2OX + DDN2O;  
            
    if(N2O[dlayer] < ZERO ) { DN2O = DNO2 + N2OX - DAN2O; }
    
    if( N2O[dlayer] < ZERO ) { N2O[dlayer] = ZERO; }
            
    DDN2 = DN2O;
    N2[dlayer] = N2X + DDN2;  // Net N2 content in layer L at T.
    DC = (U/YMC+MC) * BX * DT;  // Total consuptn of soluble C.
            
    if( ZERO == CX ) { DC = ZERO; }
            
    // Production of CO2 in layer L at T.
            
    DCO2X = DC - DBG;  
            
    // Humus formed by dead denitrifier.
            
    DHUMX = DBD; 

            
    // New denitrifier biomass in layer L at T. 
            
    B[dlayer] = BX + DB; 
            
    // New soluble C content in layer L at T.
            
    C[dlayer] = CX - DC; 
            

    // Totl N consumption by denitrification synthesis
            
    FIXN[dlayer] = DANO3 + DANO2 + DAN2O;  
            
    DCO2[dlayer] = DCO2X;
    DHUM[dlayer] = DHUMX;
  }     
  else
  {
    DCO2[dlayer] = ZERO;
    DHUM[dlayer] = ZERO;
    FIXN[dlayer] = ZERO;
  }
	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void NEM::getecd( ofstream& rflog1 ) 
{
  string ecd;

  cout << "Enter name of file (.ECD) with the NEM N2O ";
  cout << "parameter values:" << endl;
  
  cin >> ecd;

  rflog1 << "Enter name of file (.ECD) with the NEM N2O ";
  rflog1 << "parameter values:" << endl;

  getecd( ecd );
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void NEM::getecd( const string& ecd )
{

  int comtype;
  int dcmnt;
  ifstream infile;


  infile.open( ecd.c_str(), ios::in );

  if( !infile )
  {
    cerr << endl << "Cannot open " << ecd;
    cerr << " for data input" << endl;
    exit( -1 );
  }

  getXMLrootNode( infile, "nemECD" );

  for( dcmnt = 1; dcmnt < MAXCMNT; ++dcmnt )
  {
    comtype = getXMLcommunityNode( infile, "nemECD" );

    if( comtype >= MAXCMNT )
    {
      cerr << endl << "comtype is >= MAXCMNT" << endl;
      cerr << "comtype cannot be greater than ";
      cerr << (MAXCMNT-1) << endl;
      cerr << " in nemECD" << endl;
      
      exit( -1 );
    }

    initDntrf[comtype]= getXMLcmntArrayDouble( infile,
                                               "nemECD",
                                               "initDntrf",
                                               comtype );

    endXMLcommunityNode( infile );
  }

  if( dcmnt < MAXCMNT )
  {
    cerr << endl << " Parameters found for only " << dcmnt;
    cerr << " community types out of a maximum of ";
    cerr << (MAXCMNT-1) << " types in nemECD" << endl;
    
    exit( -1 );
  }

  infile.close();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void NEM::makeAmmonia( const int& dryday,
                       const int& dlayer,
                       const double& oldNH4 )
{
  double DNH3;
  double NH30;
  double NH31;
  double NH3X;
  double NH442M;
	
  if( ZERO == oldNH4 )
  {
    NH3X = ZERO;
    NH443 = ZERO;
  }
  else
  {
    // Convert NH4+ from KgN/ha/lyr to molN/l H2O
          	
    NH442M = oldNH4 / W / 14.0; 
          
    // Desolved NH3 produced fro NH4+.
          
    NH3X = 0.2 * NH442M;     
          
    // Convt to KgN/ha/layer.
          
    NH3X *= W * 14.0;    
          
    if( 0 == dryday ) { NH30 = ZERO; }
    else { NH30 = NH3[dlayer]; }
  }
        
  DNH3 = NH3X - NH30;
  NH443 = oldNH4 - DNH3;
  
  if( NH443 < ZERO ) { NH443 = ZERO; }
        
//  NH3X /= 14.0 / (W*1000.0);
  NH3X /= (14.0 * (W*1000.0));
  
  VOLNH3 = NH3X * 2.0 * pow( (0.025/3.1416), 0.5 );
  
  VOLNH3 *= 14.0 * 100000.0;

  NH3X *= 14.0 * (W*1000.0);
  
  NH31 = NH3X - VOLNH3;
  
  if( NH31 < ZERO )
  {
    VOLNH3 = NH3X;
    NH31 = ZERO;
  }
        
  NH3[dlayer] = NH31;
	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void NEM::nitrification( const int& dlayer,
                         const double& atmp, 
                         const double& moist )
{
  // Standard rate of nitrification at 35C.
  const double K35 = 25.0;  
  
  double KA;
  double NNO;
  double RM;
  
  // Effect of temperature 

  if( (atmp >= ZERO) && (atmp < 10.0) )
  { 
    KA = (0.0105*atmp + 0.00095 * pow( atmp, 2.0 ))* K35;
  }
        
  if( (atmp >= 10.0) && (atmp < 35.0) )
  {
    KA = (0.032*atmp -0.12) * K35;
  }
        
  if( (atmp >= 35.0) && (atmp <= 45.0) )
  {
    KA = (-0.1*atmp + 4.5) * K35;
  }

  // Effect of moisture.
        
  if( (moist >= ZERO) && (moist <= 0.9) ) 
  { 
    RM = 1.111 * moist; 
  }
        
  if( (moist > 0.9) && (moist <= 1.0) )
  {  
    RM = 10.0 - 10.0 * moist;
  }
        
  // NO3- produced by nitrification.
        
  NNO = NH45 * (1.0 - exp( -KA )) * RM; 
        
  // New NH4+ after nitrification,KgN/ha/layer.
        
  NH46 = NH45 - NNO;
          
  if( NH46 < ZERO ) 
  {
    NNO = NH45;
    NH46 = ZERO;
  }
        
  // New NO3- after nitrification.
        
  NO36 = NO35 + NNO;  
        
  // Nitrification rate in KgN/ha/layer.
        
  NITRIFY[dlayer] = NNO; 

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void NEM::resetInitialMonthlyConditions( void )
{
  int dlvl;
  int dday;
  
  AAAA = ZERO;
  
  BACNN = ZERO;
  BBBB = ZERO;
  BX = ZERO;

  CH4UP = ZERO;
  CLAYC = ZERO;
  CNO3 = ZERO;
  CRB1X = ZERO;
  CRB2X = ZERO;
  CRBX = ZERO;
  CRHL = ZERO;
  CRHR = ZERO;
  CRHX = ZERO;
  CX = ZERO;

  DCBAVAI = ZERO;
  DCHAVAI = ZERO;
  DCO2T = ZERO;
  DHBC = ZERO;
  DHBN = ZERO;
  DPHUMX = ZERO;
  DRCA = ZERO;
  DRCB = ZERO;
  DRNB = ZERO;
  DT = 1.0;
  DYN2EMS = ZERO;     
  DYN2OEMS = ZERO;

  EE = ZERO;

  FBCO2 = ZERO;
  FBHC = ZERO;
  FBHN = ZERO;
  FF = ZERO;
  FHCO2 = ZERO;
  FIN = ZERO;

  H = ZERO;
  HUM = ZERO;

  II = ZERO;

  KC = ZERO;
  KN = ZERO;

  LICH4UP = ZERO;

  M = ZERO;
    
  N = ZERO;
  N2OX = ZERO;
  N2X = ZERO;

  NH40 = ZERO;
  NH41 = ZERO;
  NH43 = ZERO;
  NH44 = ZERO;
  NH443 = ZERO;
  NH45 = ZERO;
  NH46 = ZERO;

  NO2X = ZERO;
  NO30 = ZERO;
  NO31 = ZERO;
  NO35 = ZERO;
  NO36 = ZERO;
  NO3X = ZERO;

  NRB = ZERO;
  NRH = ZERO;
  NRHL = ZERO;
  NRHR = ZERO;

  PHK1 = ZERO;
  PHK2 = ZERO;

  Q = 0;

  RCLX = ZERO;   
  RCRX = ZERO;
  RCVLX = ZERO;

  TE = ZERO;
  TEM1 = ZERO;
  TT = ZERO;
  TUPTN = ZERO; 
  TVNH3 = ZERO;  
  TZW = ZERO;

  V = ZERO;
  VOLNH3 = ZERO;

  W = ZERO;
  WAEN2O = ZERO;
  WFCO2 = ZERO; 
  WNITRIFY = ZERO;
  
  for( dlvl = 0; dlvl < (NLVL+1); ++dlvl )
  { 
    if( dlvl < NLVL )
    {
      AEN2O[dlvl] = ZERO;
      AVERMOIS[dlvl] = ZERO;

      CR[dlvl] = ZERO;
      CRB[dlvl] = ZERO;
      CRB1[dlvl] = ZERO;
      CRB2[dlvl] = ZERO;
      CRH[dlvl] = ZERO;

      DPHUM[dlvl] = ZERO;
      DPOC[dlvl] = ZERO;

      FCO2[dlvl] = ZERO;

      NH3[dlvl] = ZERO;
      NH4INI[dlvl] = ZERO;
      NO2[dlvl] = ZERO;
      NO3[dlvl] = ZERO;
      NO3INI[dlvl] = ZERO;

      OCINI[dlvl] = ZERO;

      RCL[dlvl] = ZERO;
      RCR[dlvl] = ZERO;
      RCVL[dlvl] = ZERO;

      SNH4[dlvl] = ZERO;

      TCAVAI[dlvl] = ZERO;

      TMPCO2DN[dlvl] = ZERO;
      TMPN2DN[dlvl] = ZERO;
      TMPN2ODN[dlvl] = ZERO;

      TOC[dlvl] = ZERO;
      TON[dlvl] = ZERO;
    }
    
    B[dlvl] = ZERO;
    BACN[dlvl] = ZERO;

    C[dlvl] = ZERO;

    DCO2[dlvl] = ZERO;
    DHUM[dlvl] = ZERO;
    DPB[dlvl] = ZERO;
    DPSC[dlvl] = ZERO;  

    EMSN2[dlvl] = ZERO;   
    EMSN2O[dlvl] = ZERO;

    FIXN[dlvl] = ZERO;

    LENO2[dlvl] = ZERO;     
    LENO3[dlvl] = ZERO;

    N2[dlvl] = ZERO;
    N2O[dlvl] = ZERO;
    NITRIFY[dlvl] = ZERO;
    NO2P[dlvl] = ZERO;
    NO3P[dlvl] = ZERO;

    ZB[dlvl] = ZERO;
    ZC[dlvl] = ZERO;
    ZN2[dlvl] = ZERO;
    ZN2O[dlvl] = ZERO;
  }

  for( dday = 0; dday < MAXMDAYS; ++dday )
  { 
    CH4N[dday] = ZERO;
    CO2DN[dday] = ZERO;
    CO2N[dday] = ZERO;

    N2DN[dday] = ZERO;
    N2ODN[dday] = ZERO;
    N2ON[dday] = ZERO;
  }     

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double NEM::setDailyCH4emissions( const float& lat,
                                  const int& dday, 
                                  const double& dayH2OtableZ )
{
  double daych4;
//  int dlayer;

  // Following assignment added by DWK on 20040816
//  double tzw = -1.0;  
  double tzw;  

  float RLAT;
  
  // Following assignment added by DWK on 20040816
  double h2oTableDepth = dayH2OtableZ; 

/*  for ( dlayer = 0; dlayer < NLAYERS; dlayer++ )
  {
    if ( dayMoist[dlayer][dm][dday] > 0.99 ) 
    {
      tzw = dayTemp[dlayer][dm][dday];
      
      switch ( dlayer )
      {
        case 0:  h2oTableDepth = 0.000000; break;
        case 1:  h2oTableDepth = 1.75; break;
        case 2:  h2oTableDepth = 4.50; break;
        case 3:  h2oTableDepth = 9.00; break;
        case 4:  h2oTableDepth = 16.50; break;
        default: h2oTableDepth = dayH2OtableZ;
      }
      
      dlayer = NLAYERS;  // break from loop
    }
  } 
*/

  if( h2oTableDepth != 9.00 ) 
  {
    h2oTableDepth = 9.00;
  }

  tzw = dayTemp[3][dday];
  RLAT = fabs( lat );
  
//  if( tzw > -1.0 && dayTemp[0][dday] > -1.0 )
  if( tzw > ZERO && dayTemp[0][dday] > ZERO )
  { 
//    DAYCH4 = exp( 3.345 + 0.069*TZW - 0.06*DAYZ0 );
  
    if( RLAT > 60.0 )
    { 
      daych4 = 96.0 
               * exp( 0.121 * (tzw-6.45006) 
                     - 0.0599 * (h2oTableDepth-7.10796) );
    }                 
    else if( RLAT >= 45.0 && RLAT <= 60.0 ) 
    {
      daych4 = 87.0 
               * exp( 0.121 * (tzw-14.51815)
                      -0.0599 * (h2oTableDepth-10.15860) );
    }                  
    else if( RLAT < 45.0 ) 
    {
      daych4 = 135.0
               * exp( 0.121 * (tzw-25.36087)
                    - 0.0599 * (h2oTableDepth-9.80447) );
    }
  }
  else { daych4 = ZERO; }

  return daych4;

	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double NEM::setMonthlyCH4emissions( const int& vegtype,
                                    const int& outmon,
                                    const double& tair,
                                    const double& prec,
                                    const double& pet,
                                    const float& lat )
{
  int dday;
  
  double ch4flx;
  double ratiot;
  double ratiop;
  
  // Assume that water table depth is no greater than 
  //   30 cm in wetlands
  
  double H2OtableZ = 30.0;
  
  switch( vegtype )
  {
    case 11: if( tair > -1.0 )
             {
               ch4flx = 8.5 * (double) mdays[outmon]; 
             }
             else { ch4flx = ZERO; }
             break;

    case 12: if( tair > -1.0 )
             {
               ch4flx = 8.5 * (double) mdays[outmon]; 
             }
             else { ch4flx = ZERO; }
             break;
             
    case 17: ratiot = ch4xtair( tair );
             ratiop = ch4xprec( prec, pet ); 
             ch4flx = 100.0 * ratiot * ratiop 
                      * (double) mdays[outmon];
             break;
             
    case 18: ratiot = ch4xtair( tair );
             ratiop = ch4xprec( prec, pet ); 
             ch4flx = 202.0 * ratiot * ratiop 
                      * (double) mdays[outmon];
             break;
             
    case 19: ratiot = ch4xtair( tair );
             ratiop = ch4xprec( prec, pet ); 
             ch4flx = 100.0 * ratiot * ratiop 
                      * (double) mdays[outmon];
             break;
             
    case 20: ratiot = ch4xtair( tair );
             ratiop = ch4xprec( prec, pet ); 
             ch4flx = 202.0 * ratiot * ratiop 
                      * (double) mdays[outmon];
             break;
             
    case 21: ch4flx = ZERO;
             for( dday = 0; dday < mdays[outmon]; ++dday )
             {
               ch4flx += setDailyCH4emissions( lat,
                                               dday, 
                                               H2OtableZ );
             }
             break;

    case 22: ch4flx = ZERO;
             for( dday = 0; dday < mdays[outmon]; ++dday )
             {
               ch4flx += setDailyCH4emissions( lat,
                                               dday, 
                                               H2OtableZ );
             }
             break;

    case 26: ratiot = ch4xtair( tair );
             ratiop = ch4xprec( prec, pet ); 
             ch4flx = 182.0 * ratiot * ratiop 
                      * (double) mdays[outmon];
             break;

    case 27: ratiot = ch4xtair( tair );
             ratiop = ch4xprec( prec, pet ); 
             ch4flx = 182.0 * ratiot * ratiop 
                      * (double) mdays[outmon];
             break;

    case 28: ratiot = ch4xtair( tair );
             ratiop = ch4xprec( prec, pet ); 
             ch4flx = 182.0 * ratiot * ratiop 
                      * (double) mdays[outmon];
             break;

    case 29: ratiot = ch4xtair( tair );
             ratiop = ch4xprec( prec, pet ); 
             ch4flx = 182.0 * ratiot * ratiop 
                      * (double) mdays[outmon];
             break;
             
    default: ch4flx = ZERO;    	
  }
  
  // Change mg CH4 to g CH4
  
  return (ch4flx * 0.001);
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void NEM::stepday( const int& DY, 
                   const int& DAY,
                   const double& clay,
                   const double& CLAYC,
                   const double& toppor,
                   double AVTEMP[CLMNLAYERS][MAXMDAYS],
                   double AVWTC[CLMNLAYERS][MAXMDAYS] )
{

  // Constants
  
  //  Decomposition rate factor. When tillage occurs, 
  //    DRF is increased.
  const double DRF = 0.02;

  // Percentage of labile humads in total humads.
  const double SRH = 0.16;
     

  double AEN2OX;
  
  double ATEMP;
              
  // Daily nitri N2O emission in profile,KgN/ha/day.
  double DAYAEN2O;
  
  double DDRF;
  
  float DL;
      
  // Daily N uptake in profile,KgN/ha/day.
  double DTUPTN;
  
  // Daily NH3 votalization in profile,KgN/ha/day.
  double DTVNH3;
  
  // Daily nitri CO2 emission in profile,KgC/ha/day.
  double DWFCO2;
         
//  int IJD;
    
  int L;
          
  double MOIS;
  
  double NH442;
  double NH443A;
  double NH444;
  double NH47;

  double NO33;
  double NO34;
  double NO37;

  // Proportions of a 10 cm soil layer
  double plyr1;
  double plyr2;
  double plyr3;
  double plyr4;

  double RFM;
  double RFMM;
  double RFMT;
    
  // NH4+,KgN/ha/layer.
  double TNH4;
  
  // NO3-,KgN/ha/layer.
  double TNO3;
  
  // Accumulated N uptake for the dry period,KgN/ha.
  double TUPTN;
      
  double UPTNL;
           
  // Biomass C in profile in KgC/ha/day.
  double WCRB;
  
  // Humads C in profile,KgC/ha/day.
  double WCRH;
  
  // Humus C in profile,KgC/ha/day.  
  double WHMUS;
  
  // NH4+ mount in profile,KgN/ha/day.
  double WNH4;
    
  // NO3- mount in profile,KgN/ha/day.
  double WNO3;

  double WRC;
    
  // Labile residue C in profile,KgC/ha/day.
  double WRCL;
  
  // Resistant residue C in profile,KgC/ha/day.
  double WRCR;
  
  // Very labile residue C in profile,KgC/ha/day.
  double WRCVL;
  
  // Soluble C in profile,KgC/ha/day.
  double WTCAVAI;
  
  // Active C(bio+humads) in profile,KgC/ha/day.
  double WTOC;
  
  // Active N(bio+humads) in profile,KgC/ha/day.
  double WTON;
      
  WCRB = ZERO;  
  WRCVL = ZERO; 
  WRCL = ZERO;  
  WRCR = ZERO;  
  WTOC = ZERO;  
  WTCAVAI = ZERO;  
  WCRH = ZERO;  
  WTON = ZERO;  
  WHMUS = ZERO; 
  WNH4 = ZERO;  
  WNO3 = ZERO;  
  DAYAEN2O = ZERO; 
  DWFCO2 = ZERO; 
  DTUPTN = ZERO; 
  DTVNH3 = ZERO; 
  
  for( L = 0; L < Q; ++L )
  {
    if( L < 1 ) { DDRF = DRF; }
    
    if( (L >= 1) && (L < 4) ) { DDRF = DRF / 2.0; }
    
    if( L >= 4 ) { DDRF = DRF / 4.0; }

    // set layer temp and moist normalized to 10 layers.
    
    DL = (float) (L+1); 
    
    if( (DL*H) <= 0.1 )
    {
      if( (DL*H) <= 0.0175 )
      {
        ATEMP = AVTEMP[0][DY];
        MOIS = AVWTC[0][DY];
      }
            
      if( ((DL*H) > 0.0175 ) && ((DL*H) <= 0.045) ) 
      {
        plyr1 = 0.0175 / (DL*H);
        plyr2 = 1.0 - plyr1;

        ATEMP = (plyr1 * AVTEMP[0][DY]) 
                + (plyr2 * AVTEMP[1][DY]);

        MOIS = (plyr1 * AVWTC[0][DY]) 
               + (plyr2 * AVWTC[1][DY]);
      }
       
      if( ((DL*H) > 0.045) && ((DL*H) <= 0.090) ) 
      {
        plyr1 = 0.0175 / (DL*H);
        plyr2 = (0.045 - 0.0175) / (DL*H);
        plyr3 = 1.0 - plyr1 - plyr2;

        ATEMP = (plyr1 * AVTEMP[0][DY]) 
                + (plyr2 * AVTEMP[1][DY])
                + (plyr3 * AVTEMP[2][DY]);

        MOIS = (plyr1 * AVWTC[0][DY]) 
               + (plyr2 * AVWTC[1][DY])
               + (plyr3 * AVWTC[2][DY]);
      }

      if( ((DL*H) > 0.090) ) 
      {
        plyr1 = 0.0175 / (DL*H);
        plyr2 = (0.045 - 0.0175) / (DL*H);
        plyr3 = (0.090 - 0.045) / (DL*H);
        plyr4 = 1.0 - plyr1 - plyr2 - plyr3;

        ATEMP = (plyr1 * AVTEMP[0][DY]) 
                + (plyr2 * AVTEMP[1][DY])
                + (plyr3 * AVTEMP[2][DY])
                + (plyr4 * AVTEMP[3][DY]);

        MOIS = (plyr1 * AVWTC[0][DY]) 
               + (plyr2 * AVWTC[1][DY])
               + (plyr3 * AVWTC[2][DY])
               + (plyr4 * AVWTC[3][DY]);
      }

    }
        
    if( ((DL*H) > 0.1) && ((DL*H) <= 0.2) ) 
    {
      if( (DL*H) <= 0.165 ) 
      {
        ATEMP = AVTEMP[3][DY];
        MOIS = AVWTC[3][DY];
      }        
      else 
      {
        plyr2 = ((DL*H) - 0.165) / ((DL*H) - 0.1);
        plyr1 = 1.0 - plyr2;

        ATEMP = (plyr1 * AVTEMP[3][DY])
                + (plyr2 * AVTEMP[4][DY]);

        MOIS = (plyr1 * AVWTC[3][DY])
               + (plyr2 * AVWTC[4][DY]);
      }
    }

    if( ((DL*H) > 0.2) && ((DL*H) <= 0.3) ) 
    {
      if( (DL*H) <= 0.290 ) 
      {
        ATEMP = AVTEMP[4][DY];
        MOIS = AVWTC[4][DY];
      }
      else
      {
        plyr2 = ((DL*H) - 0.290) / ((DL*H) - 0.2);
        plyr1 = 1.0 - plyr2;

        ATEMP = (plyr1 * AVTEMP[4][DY])
                + (plyr2 * AVTEMP[5][DY]);

        MOIS = (plyr1 * AVWTC[4][DY])
               + (plyr2 * AVWTC[5][DY]);
      }

    }

    if( (DL*H) > 0.3 ) 
    {
      ATEMP = AVTEMP[5][DY];
      MOIS = AVWTC[5][DY];
    }

    // Convert MOIS from volumetric soil moisture to 
    //   proportion of porosity

    MOIS /= toppor;

    // Ensure that MOIS is never above 100% porosity

    if( MOIS > 1.0 ) { MOIS = 1.000000; }

    W = V * MOIS;  

    RCVLX = RCVL[L];
    RCLX = RCL[L];
    RCRX = RCR[L];

    if( 0 == DAY ) 
    {
      // set begin year conditions
         
      CRB1X = OCINI[L] * RBO * SRB; 
      CRB2X = OCINI[L] * RBO * (1 - SRB); 
      CRBX = CRB1X + CRB2X;   
      CRHX = OCINI[L] - CRBX; 
      CRHL = CRHX * SRH;      
      CRHR = CRHX * (1 - SRH); 
      NRHR = CRHR / RCNH;     
      NRHL = CRHL / RCNH;     
      NRB = CRBX / RCNB;      
      TNO3 = NO3INI[L];     
      TNH4 = NH4INI[L];      
    }
    else
    {
      CRHX = CRH[L];
      CRHL = CRHX * SRH;
      CRHR = CRHX * (1 - SRH);
      NRHR = CRHR / RCNH;
      NRHL = CRHL / RCNH;
      CRB1X = CRB1[L];
      CRB2X = CRB2[L];
      CRBX = CRB1X + CRB2X;
      NRB = CRBX / RCNB;
      TNO3 = NO3[L];
      TNH4 = SNH4[L];
    }


    // Calculating temperature and moisture reduction factors.
        
    if( MOIS < 0.1 ) { RFMM = MOIS*0.1; }
    
    if( (MOIS >= 0.1) && (MOIS <= 0.6) )
    { 
      RFMM = (MOIS-0.1) * 1.96 + 0.02;
    }
        
    if( (MOIS >=0.6) && (MOIS <= 0.8) )
    { 
      RFMM = 1.0 - (MOIS - 0.6) * 2.5;
    }
        
    if( MOIS > 0.8 ) { RFMM=0.5 - (MOIS - 0.8) * 0.5; }
    
    if( ATEMP < ZERO ) { ATEMP = ZERO; }
    
    if( ATEMP < 30.0 ) { RFMT = ATEMP * 0.06; }
    
    if( (ATEMP >= 30.0) && (ATEMP <= 40.0) ) { RFMT = 1.8; }
    
    if( ATEMP > 40.0 ) 
    { 
      RFMT = 1.8 - (ATEMP - 40.0) * 0.04; 
    }
        
    RFM = RFMM * RFMT;

    // NH40 and NO30 are supposed to be the addition of 
    //   N-fertilizers.
        
    NH40 = TNH4;    
    NO30 = TNO3;    


    decomposeResidues( RFM, DDRF );
                
    // Update carbon and nitrogen pools in resistant
    //   microbial biomass after residue decomposition
        
    CRB1X += (DRCB * SRB);        
    CRB2X += (DRCB * (1.0 - SRB));  
    NRB += DRNB; 
    

    decomposeMicrobialBiomass( RFM, DDRF );
                 
    // Update carbon and nitrogen pools in resistant humads
    //   after microbial biomass decomposition
        

    CRHR += FBHC;           
    NRHR += FBHN;   


    // New NO3- pool after humads decomposition 
    //   (i.e. No change)

    NO33 = NO31;


    decomposeHumads( CLAYC, RFM, DDRF );
                     

    CRB1X += (DHBC * SRB);
    CRB2X += (DHBC * (1.0 - SRB));
    CRBX = CRB1X + CRB2X;
    NRB += DHBN;

        
    // New NO3- pool after humads decomposition 
    //   (i.e. No change)
        
    NO34 = NO33;        
    
    if( NO34 < ZERO ) { NO34 = ZERO; }
        
    // Update NH4 pool after adsorption of NH4+ by clay.

    NH442 = adsorbAmmonium( clay, NH44, OCINI[L] );
        

    // NH3 production 

    makeAmmonia( DAY, L, NH442 );

        
    NH443A = (NH443 / M) * 1000000000.0;

    AEN2OX = 0.0002 * NH443A * ((0.54 + 0.51*ATEMP) /15.8);
        
    // N2O production.
        
    AEN2OX *= (M / 1000000000.0); 

    NH444 = NH443 - AEN2OX;  // New NH4+ in KgN/Ha/layer.
    NH45 = NH444;
    NO35 = NO34;

    // Nitrification

    nitrification( L, ATEMP, MOIS );

       
    // Specify C & N pools in each layer on DAY, KgC 
    //   or N/ha/layer
    
    CRB1[L] = CRB1X;
    CRB2[L] = CRB2X;
    
    CRB[L] = CRB1[L] + CRB2[L]; // Biomass C.
    
    CRH[L] = CRHX;              // Humads C.
    RCVL[L] = RCVLX;
    RCL[L] = RCLX;
    RCR[L] = RCRX;
    
    CR[L] = RCVL[L] + RCL[L] + RCR[L];  // Residues C.
    
    TOC[L] = CRB[L] + CRH[L];           // Active organic C.
    
    TCAVAI[L] = DCBAVAI + DCHAVAI;      // Soluble C.
    
    AEN2O[L] = AEN2OX;                  // N2O produced.
    
    FCO2[L] = DRCA + FBCO2 + FHCO2;     // CO2 produced.
    
    TON[L] = NRB + NRH;                 // Active organic N.
    
    DPHUM[L] = DPHUM[L] + DPHUMX;       // New humads.

    // Total nitri rate in profile for the whole year,KgN/ha.
        
    WNITRIFY += NITRIFY[L]; 
                                            
    WCRB += CRB[L];         // Total biomass in profile.
    WRCVL += RCVL[L];
    WRCL += RCL[L];
    WRCR += RCR[L];         // Residue C in profile.
    WTOC += TOC[L];         // Active organic C in profile.
    WTCAVAI += TCAVAI[L];   // Soluble C in profile.
    WCRH += CRH[L];         // Humads N in profile.
    WTON += TON[L];         // Active organic N in profile.
    WHMUS += DPHUM[L];      // Humus C in profile.

    DAYAEN2O += AEN2O[L];   // Daily N2O prdcd in profile.
    DWFCO2 += FCO2[L];      // Daily CO2 prdcd in profile.
        
    // Daily NH3 emsion from profile,KgN/ha/day
        
    DTVNH3 += VOLNH3; 

    // Accumled CO2 prdcd in profile.
        
    WFCO2 += FCO2[L];       
        
    // Accumulated N2O prdcd in profile.
    //   Accumulated is for the dry period.
        
    WAEN2O += AEN2O[L];     
        
    // Accumlted NH3 emissions fro profile.
    
    TVNH3 = VOLNH3;   

    // model for uptake calc needs developed

    NO37 = NO36; 
    NH47 = NH46; 
    UPTNL = ZERO;    // No uptake below the deepth 30cm.

    DTUPTN += UPTNL;    // Daily.
    TUPTN += UPTNL;     //Accumulated uptake from all layers.
    NO3[L] = NO37;
    SNH4[L] = NH47;
    WNO3 += NO3[L];   // Daily NO3- in profile.
    WNH4 += SNH4[L];  // Daily NH4+ in profile.
  }

  // Total residue C in profile on day DAY.
      
  WRC = WRCVL + WRCL + WRCR;   

  uptakeMethane( DY, WRC, WTOC, WHMUS, AVTEMP );
      
  CH4N[DY] = CH4UP;             // Specify output,KgC/ha/day.
  N2ON[DY] = DAYAEN2O * 1000.0; // Specify output,gN/ha/day.
  CO2N[DY] = DWFCO2;            // Specify output,KgC/ha/day.

  return;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void NEM::stephour( const int& dm,
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
                    const int& maxlayers )
{
  double AD;

  float DL;

  int L;
  int LLL;

  int NOUTBD;

  double PA;
  double PN2; 
  double PN2O;

  // Proportions of a 10 cm soil layer
  double plyr1;
  double plyr2;
  double plyr3;
  double plyr4;

  int T;
        
  T = dhr + (rainrun * 24 ) + 1;
        
  NOUTBD = dndy;
  if( dndy >= mdays[dm] ) { dndy = mdays[dm] - 1; }
        
  for( L = 0; L < maxlayers; ++L )
  {
    DL = (float) (L+1);

    if( (DL*H) <= 0.1 )
    {
      if( (DL*H) <= 0.0175 )
      {
        AVERMOIS[L] = HRWTC[0][dndy][dhr];
      }
            
      if( ((DL*H) > 0.0175 ) && ((DL*H) <= 0.045) ) 
      {
        plyr1 = 0.0175 / (DL*H);
        plyr2 = 1.0 - plyr1;

        AVERMOIS[L] = (plyr1 * HRWTC[0][dndy][dhr])
                      + (plyr2 * HRWTC[1][dndy][dhr]);
      }
       
      if( ((DL*H) > 0.045) && ((DL*H) <= 0.090) ) 
      {
        plyr1 = 0.0175 / (DL*H);
        plyr2 = (0.045 - 0.0175) / (DL*H);
        plyr3 = 1.0 - plyr1 - plyr2;

        AVERMOIS[L] = (plyr1 * HRWTC[0][dndy][dhr])
                      + (plyr2 * HRWTC[1][dndy][dhr])
                      + (plyr3 * HRWTC[2][dndy][dhr]);
      }

      if( ((DL*H) > 0.090) ) 
      {
        plyr1 = 0.0175 / (DL*H);
        plyr2 = (0.045 - 0.0175) / (DL*H);
        plyr3 = (0.090 - 0.045) / (DL*H);
        plyr4 = 1.0 - plyr1 - plyr2 - plyr3;

        AVERMOIS[L] = (plyr1 * HRWTC[0][dndy][dhr])
                      + (plyr2 * HRWTC[1][dndy][dhr])
                      + (plyr3 * HRWTC[2][dndy][dhr])
                      + (plyr4 * HRWTC[3][dndy][dhr]);
      }
    }
        
    if( ((DL*H) > 0.1) && ((DL*H) <= 0.2) ) 
    {
      if( (DL*H) <= 0.165 ) 
      {
        AVERMOIS[L] = HRWTC[3][dndy][dhr];
      }        
      else 
      {
        plyr2 = ((DL*H) - 0.165) / ((DL*H) - 0.1);
        plyr1 = 1.0 - plyr2;

        AVERMOIS[L] = (plyr1 * HRWTC[3][dndy][dhr])
                      + (plyr2 * HRWTC[4][dndy][dhr]);
      }
    }

    if( ((DL*H) > 0.2) && ((DL*H) <= 0.3) ) 
    {
      if ( (DL*H) <= 0.290 ) 
      {
        AVERMOIS[L] = HRWTC[4][dndy][dhr];
      }
      else
      {
        plyr2 = ((DL*H) - 0.290) / ((DL*H) - 0.2);
        plyr1 = 1.0 - plyr2;

        AVERMOIS[L] = (plyr1 * HRWTC[4][dndy][dhr])
                      + (plyr2 * HRWTC[5][dndy][dhr]);
      }
    }

    if( (DL*H) > 0.3 ) 
    {
      AVERMOIS[L] = HRWTC[5][dndy][dhr];
    }

    // Convert soil moisture for volumetric soil moisture to
    // proportion of total porosity

    AVERMOIS[L] /= toppor;

    if( AVERMOIS[L] > 1.0 ) { AVERMOIS[L] = 1.000000; }
  }
        
  dndy = NOUTBD;

  //CAS if ( T < maxlayers ) { LLL = T; }
  //CAS else { LLL = maxlayers; }
  LLL = maxlayers;

  // each layer calc for denitrification.

  for( L = 0; L < LLL; ++L )
  {
    // Define soluble C(KgC/layer/ha) and denitrifier 
    //   in layer L and hr T.

    if( T == (L+1) ) 
    {
      CX = DPSC[L];
      BX = DPB[L];
    }
    else
    {
      CX = ZC[L];
      BX = ZB[L];
    }

    // Define desolved NO3- in layer L at hr T. 

    if( T <= (int) tt ) // During raining hour.
    {
      if ( 0 == L ) 
      {
        // Desolved NO3- content in layer L at hr T.
              
        NO3X = CNO3 + BACN[L] 
               * pow( (1.0-lf), (T-1) ) * lf;        
      }
      else
      {
        NO3X = LENO3[L-1] + BACN[L] 
               * pow( (1.0-lf), (T-(L+1)) ) * lf;
      }
            

      // Total disolved NO3-.
            
      BACNN += BACN[L]* pow( (1.0-lf), (T-(L+1)) ) * lf;  
    }
    else
    {
      // NO3- is the same as the end of last hr.
            	
      NO3X = LENO3[L]; 
    }

    // Define NO2- in layer L at hour T.

    if( 0 == L ) { NO2X = ZERO; }
    else { NO2X = LENO2[L-1]; }
          
    if( T > (int) tt ) { NO2X = LENO2[L]; }

    // Define N2O and N2 in layer L at hour T.

    if( T == (L+1) ) 
    {
      N2OX = ZERO;
      N2X = ZERO;
    }
    else
    {
      N2OX = ZN2O[L];
      N2X = ZN2[L];
    }


    // Denitrify.

    denitrification( pcmnt, L, AVERMOIS[L] );


    // Calculate diffusion of N2O and N2 from soil.

      //cout << "Rain duration for denitrification = " << tt << endl;
    if( T > tt ) 
    {
      //cout << "Diffusion of N2O occurring at T = " << T << endl;
      //cout << "Soil level = " << L << endl;
      //cout << "Averaged soil moisture  = " << AVERMOIS[L] << endl;

      PA = 1.0 - AVERMOIS[L];
      
      AD = 2.0 * (clay * 0.01) / 0.63;
          
      PN2O = (0.0006 + 0.0013 * AD) 
             +(0.013- 0.005 * AD) * PA 
             * pow( 2.0, (tm/20.0) );
                   
      PN2 = 0.017 + (0.025 - 0.0013 * AD) * PA 
            * pow( 2.0, (tm/20.0) );
                  
      EMSN2O[L] = N2O[L] * PN2O;
      EMSN2[L] = N2[L] * PN2;
      N2O[L] -= EMSN2O[L];
      N2[L] -= EMSN2[L];
    }
    

    if( B[L] < ZERO ) { B[L] = ZERO; }
    
    if( C[L] < ZERO ) { C[L] = ZERO; }
    
    if( EMSN2O[L] < ZERO ) { EMSN2O[L] = ZERO; }
    
    if( EMSN2[L] < ZERO ) { EMSN2[L] = ZERO; }
    
    if( DCO2[L] < ZERO ) { DCO2[L] = ZERO; }
    
    if( DHUM[L] <  ZERO ) { DHUM[L] = ZERO; }
    
    if( NO3P[L] < ZERO ) { NO3P[L] = ZERO; }
    
    if( NO2P[L] < ZERO ) { NO2P[L] = ZERO; }
    
    if( N2O[L] < ZERO ) { N2O[L] = ZERO; }
    
    if( N2[L] < ZERO ) { N2[L] = ZERO; }


    // Calculate total pools or fluxes in profile.

    HUM += DHUM[L];   // Accumlted humus in profile,KgC/ha.
    II += DCO2[L];    // Accumlted CO2 emission in profile.
    EE += EMSN2O[L];  // Accumulated N2O emitted in a rain event.
    FF += EMSN2[L];   // Accumulated N2 emitted in a rain event.
    FIN += FIXN[L];   // Accumltd N consumt by denitri synthesis.
          
    if( (L == (maxlayers-1) ) && (T <= ttt) ) { AAAA += NO3P[L]; }
    
    if( (L == (maxlayers-1) ) && (T <= ttt) ) { BBBB += NO2P[L]; } 
          
    DCO2T += DCO2[L];  // Daily CO2 emission.
    DYN2OEMS += EMSN2O[L];  // Daily.
    DYN2EMS += EMSN2[L];     // Daily.
    ZN2O[L] = N2O[L];
    ZN2[L] = N2[L];
    ZB[L] = B[L];
    ZC[L] = C[L];
  }

  for( L = 0; L < LLL; ++L )
  {          
    LENO3[L] = NO3P[L];   
    LENO2[L] = NO2P[L];   
  }
	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void NEM::stepmonth( const int& dm,
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
                     double RAININTE[MAXMDAYS] )

{
// formerly subroutine DNDCPT(KM,NJD,NBEG,NEND)

  // Constants
  
  // Initial NO3- in surface soil in mg N/Kg.
  //const double CM10N11 = 100.0;
  //CAS: Value changed to reduce Jan. shock effect
  const double CM10N11 = 1.0;

  // Tillage reduces proportion of denitrifiers in total microbial
  //   biomass. It's reduction factor.
  const double FD = 0.05;

  // Half-saturation values in Monod equation, Kg C/M^3 water.
  const double KCI = 0.017;

  // Half-saturation values in Monod equation, Kg N/M^3 water.
  const double KNI = 0.083;
    
  // NO3- concentration in rain water in mg/m^3.
  const double RNO3 = ZERO;

  //  Initial NH4+ in surface soil in mg N/Kg.
  //const double SNH41 = 5.0;
  //CAS: Value changed to reduce Jan. shock effect
  const double SNH41 = 0.5;

  // Total thickness of soil profile in m.
  const double TH = 0.3; 

               
  double CHDP;
        
  // Total CO2 emission from all rain events    
  double CSLCO2;
  
  // Total N2 emission from all rain events
  double CSLN2;

  // Total N2O emission from all rain events,KgN/ha.
  double CSLN2O;
      
//  int DAY;
   
  double DAYHUMUS;
   
//  int DAYS;

  int DDAY;
               
  int DHOUR;
      
  int DNDY;
      
  int DY;
      
  double FACTOR;
      
  // Effect of applying anhydrous ammonia on soil pH and soluble C
  double FERENH;
        
  double FINALSN;
          
  //  Total No3- input from rain fall of this rain event,KgN/ha.
  double IIII;

  // Totl NO3- input from all rain events,KgN/ha.
  double IIIII;
  
//  int IJD;
  
  double ININO3;
  double INIOC;
  double INITB;
  double INITSC;
//  double INTOTALN;

  int K;
  
  double kgCin;

//  int KJN;
  
//  int KRDP1;
  
  int L;
  
  // Percent of disolved NO3- by rain water flux in each layer.
  double LF;
  
  // CO2 nitri emission for a yr,KgC/ha.
  double LICO2;
  
  // Totl nitri N2O emission for dry periods.
  double LIN2O;
  
  // Totl NH3 emission for dry periods,KgN/ha.
  double LINH3;
  
  // Totl N uptake for dry periods,KgN/ha.
  double LIUPTS;

  int LN;
          
//  int NADD;
  
  int NBEG;
  int NEND;
  
  int NEOM;
  
//  int NJD;
     
  int NREP;
       
//  int RAINNUM;

//  int RAINRUN;

  double RATIO;
  
  int REC;
    
  // Index to indicate the distribution of organic carbon in the 
  //   top 1 meter of soil
  int soilCdist;
  
//  int T;

  // Temperature at surface soil
  double TM;
      
  double TOPR;
//  int TRN;
 
  // Initial organic C in surface soil in Kg C/Kg.
  double TTOO;
  
  int TTT;
       
  double VS;
  
          
  resetInitialMonthlyConditions();

  // Determine soil characteristics 
  
  CLAYC = log( 0.14 / (clay * 0.01) ) / 2.3026 + 1.0;   

  
  H = 4.0 * 0.005 / totpor;   
  
  M = density * 1000.0 * H * 10000.0;  
  
  V = H * 10000.0 * totpor;   

  
  // Determine number of soil layers
  
  Q = (int)(TH / H);          

  LN = 0;

  NBEG = 0;
  NEND = mdays[dm] - 1;
 
  TEM1 = TEM[NBEG];    

  if( 0 == dm )  // First month 
  {     
    switch ( IVEG )
    {
      case  0: soilCdist = 1;
      case  1: soilCdist = 1;
      case  2: soilCdist = 1;
      case  3: soilCdist = 1;
      case  5: soilCdist = 1;
      case  6: soilCdist = 1;
      case  7: soilCdist = 1;
      case  8: soilCdist = 1;
      case 11: soilCdist = 1;
      case 12: soilCdist = 1;
      case 15: soilCdist = 2;
      default: soilCdist = 3;
    }
     	
    // *************Initialize each layer**********
  
    // Convert soil organic carbon from g C/m^2 to kg C/m^2

//    kgCin = CIN * 0.001 * 2.2;
    kgCin = CIN * 0.001;

    CARBONTR( kgCin, soilCdist, H, TOPR, TOPH, CHDP, TTOO, RATIO );
    
    // Determine amount of organic residue in top layer of soil    
    CR[0] = TTOO * (2.0/5.0) * M;           
    RCVL[0] = 0.0;                
    RCL[0] = 0.0;               
    RCR[0] = CR[0];             
    
    // Determine amount of active organic C(bio+humads) in 
    //   top layer of soil

    OCINI[0] = TTOO * (2.0/5.0) * M;   

    //if (firstn2oyr == 0)
    //{
    
    // Determine initial NO3- in top layer of soil.
    
    //NO3INI[0] = CM10N11 * M * 0.000001;
    
    // Determine initial NH4+ in top layer of soil.
     
    //NH4INI[0] = SNH41 * M * 0.000001;   

    //firstn2oyr = 1;
    //}
    
    // Determine amount of Humus C in top layer of soil.
    
    DPHUM[0] = TTOO/5.0 * M;         
    
    for( L = 0; L < Q; ++L )
    {
      // Initialize carbon stocks in other soil layers
      
      if( (L+1) <= TOPH ) 
      {
        RCVL[L] = RCVL[0];
        RCL[L] = RCL[0];
        RCR[L] = RCR[0];
        CR[L] = RCVL[L] + RCL[L] + RCR[L];
        OCINI[L] = OCINI[0];
        NO3INI[L] = NO3INI[0];
        NH4INI[L] = NH4INI[0];
        DPHUM[L] = DPHUM[0];
        LENO2[L] = ZERO;
        LENO3[L] = ZERO;
        NO3P[L] = ZERO;
        NO2P[L] = ZERO;
      }
      else
      {
        if( 11 == IVEG
            || 12 == IVEG 
            || 6 == IVEG )
        {
          FACTOR = (1.0 - RATIO * (1.0 / pow( TOPR, (TOPH*H/CHDP) ) 
                   - 1.0/ pow( TOPR, ((L+1)*H/CHDP) )));  
     
          if( FACTOR < ZERO ) { FACTOR = ZERO; }              
         
          RCVL[L] = RCVL[0] * FACTOR;                        
          RCL[L] = RCL[0] * FACTOR;                          
          RCR[L] = RCR[0] * FACTOR;                          
          CR[L] = RCVL[L] + RCL[L] + RCR[L];                   
          OCINI[L] = OCINI[0] * FACTOR;                      
          NO3INI[L] = NO3INI[0] * FACTOR;                    
          NH4INI[L] = NH4INI[0] * FACTOR;                    
          DPHUM[L] = DPHUM[0] * FACTOR;
        }                      
        else
        {
          RCVL[L] = RCVL[0] / pow( TOPR, (((L+1)-TOPH)*H/CHDP) );
          RCL[L] = RCL[0] / pow( TOPR, (((L+1)-TOPH)*H/CHDP) );
          RCR[L] = RCR[0] / pow( TOPR, (((L+1)-TOPH)*H/CHDP) );
          CR[L] = RCVL[L] + RCL[L] + RCR[L];
          OCINI[L] = OCINI[0] / pow( TOPR, (((L+1)-TOPH)*H/CHDP) );
          NO3INI[L] = NO3INI[0] / pow( TOPR, (((L+1)-TOPH)*H/CHDP) );
          NH4INI[L] = NH4INI[0] / pow( TOPR, (((L+1)-TOPH)*H/CHDP) );
          DPHUM[L] = DPHUM[0] / pow( TOPR, (((L+1)-TOPH)*H/CHDP) );
        }                                           
      }
    }    
//    if(firstn2oyr == 1)
//    {
//      for( L = 0; L < Q; ++L )  
//      {
      // Initialize carbon and nitrogen stocks for soil layers
      //   during the rest of the year 
            
//        NO3INI[L] = ANO3IN[L];
//        NH4INI[L] = ANH4IN[L];
//      }
//    }                  

  }
  else
  {
    for( L = 0; L < Q; ++L )  
    {
      // Initialize carbon and nitrogen stocks for soil layers
      //   during the rest of the year 
            
      RCVL[L] = RCVLIN[L];              
      RCL[L] = RCLIN[L];                
      RCR[L] = RCRIN[L];                
      
      CR[L] = RCVL[L] + RCL[L]+ RCR[L];
      
      OCINI[L] = OCIN[L];
      NO3INI[L] = ANO3IN[L];
      NH4INI[L] = ANH4IN[L];
      DPHUM[L] = DPHUMIN[L];
    }
  }                  


  
  // *************Initialize each layer**********

  for( K = 0; K < NLVL; ++K )
  {
    TMPCO2DN[K] = ZERO;
    TMPN2ODN[K] = ZERO;
    TMPN2DN[K] = ZERO;
  }

  for( K = 0; K < mdays[dm]; ++K )
  {
    N2ON[K] = ZERO;
    N2ODN[K] = ZERO;
    N2DN[K] = ZERO;
    CO2N[K] = ZERO;
    CO2DN[K] = ZERO;
  }

  // Initialize carbon and nitrogen fluxes to zero
  
  LICH4UP = ZERO;   
  LICO2 = ZERO;     
  LIUPTS = ZERO;    
  LINH3 = ZERO;     
  LIN2O = ZERO;     
  WNITRIFY = ZERO;  
  IIIII = ZERO;       
  CSLN2O = ZERO;    
  CSLN2 = ZERO;     
  CSLCO2 = ZERO;    

  DY = 0;
  DDAY = 0;
  DNDY = 0;

  // Max. number for repeating denitri calculation
    
  NREP = 10;         

  //cout << "NEM: Month = " << dm << endl;
  //cout << "NEM: Veg. Type = " << IVEG << endl;

  for( DY = 0; DY <= NEND; ++DY )
  {
    if( DY >= 1 ) 
    {
      // Determine average daily soil temperature of the 
      //   top 20 cm

      TEM1 = (0.0175/0.200) * AVTEMP[0][DY] 
             + (0.0275/0.200) * AVTEMP[1][DY]
             + (0.0450/0.200) * AVTEMP[2][DY]
             + (0.075/0.200) * AVTEMP[3][DY]
             + (0.035/0.200) * AVTEMP[4][DY];
      
  //cout << "NEM: Day = " << DY << endl;
  //cout << "NEM: Avg. Soil Temp. = " << TEM1 << endl;

      // Initializing for second plus rain events.
      for( L = 0; L < Q; ++L )
      {
        CR[L] = RCVL[L] + RCL[L] + RCR[L];
        
        OCINI[L] = TOC[L]; 
        
        NO3INI[L] = NO3[L] + NO2[L];
        
        NH4INI[L] = SNH4[L];
      }
    }


    // Initialize accumulated carbon and nitrogen fluxes to zero
    
    TUPTN = ZERO;        
    TVNH3 = ZERO;      
    WAEN2O = ZERO;     
    WFCO2 = ZERO;  
   
    stepday( DY, 
             DDAY, 
             clay, 
             CLAYC, 
             toppor,
             AVTEMP, 
             AVWTC );

    //cout << "NEM: Raindur = " << RAINDUR[DY] << endl;
    //cout << "NEM: QStorm = " << RAININTE[DY] << endl;

    if( RAINDUR[DY] > ZERO && RAININTE[DY] > ZERO )
    {
      //  This block guarantees previous rain event does not
      //  affect current rain event except for NO3- and NO2- above layer LN.
    
      for( L = 0; L < Q; ++L )
      {
        B[L] = ZERO;
        C[L] = ZERO;
        N2O[L] = ZERO;
        N2[L] = ZERO;
        EMSN2O[L] = ZERO;
        EMSN2[L] = ZERO;
        DCO2[L] = ZERO;
        DHUM[L] = ZERO;
        FIXN[L] = ZERO;
        ZB[L] = ZERO;
        ZC[L] = ZERO;
        ZN2O[L] = ZERO;
        ZN2[L] = ZERO;
        LENO3[L] = ZERO;
        LENO2[L] = ZERO;
        NO3P[L] = ZERO;
        NO2P[L] = ZERO;
      }

      TT = RAINDUR[DY];

      DNDY = DY; 
      DDAY = 0; 
      NEOM = 0;

      //cout << "NEM: TT = " << TT << endl;
      //cout << "NEM: Duration of Storm = " << RAINDUR[DY] << endl;
      //cout << "NEM: Day of dry period = " << DY << endl;
      //cout << "NEM: Intensity of Storm = " << RAININTE[DY] << endl;


      // Total N uptake for all the Dry periods.
    
      LIUPTS += TUPTN;   
      LINH3 += TVNH3;
      LIN2O += WAEN2O;
      LICO2 += WFCO2;

      DAYHUMUS = ZERO;
      FINALSN = ZERO;
    
      for( L = 0; L < Q; ++L )
      {
        DAYHUMUS += DPHUM[L];
        FINALSN += NO3[L] + SNH4[L];
      }

      // rain event denitrification submodel.

      if( (2.0*RAININTE[DY]/totpor) > Ksat )
      {
        LF = (100.0*Ksat)/(100.0*Ksat + totpor);
      }
      else
      {
        LF = (200.0*RAININTE[DY] / totpor)
              /(200.0*RAININTE[DY] / totpor+totpor);
      }

      DT = 1.0;    // Time step 1 hr.
      VS = H *10000.0;
      W = VS * totpor;
      KC = KCI * W;   // KgC/layer/ha.
      KN = KNI * W;   // KgN/layer/ha.
      CNO3 = RNO3 * W / 1000.0;  // KgN/layer/ha.

      TTT = (int) TT; 

      //cout << "NEM: TTT = " << TTT << endl;

      PHK1 = 7.14 * (pH - 3.8) / 22.8;
      PHK2 = 7.22 * (pH - 4.4) / 18.8;
      
      if( PHK1 < ZERO ) { PHK1 = ZERO; }
      if( PHK2 < ZERO ) { PHK2 = ZERO; }

      if ( RAINDUR[DY] > 24.0 ) 
      { 
        IIII = RAININTE[DY] * 24.0 * 10.0 * RNO3; 
      }
      else {IIII = RAININTE[DY] * RAINDUR[DY] * 10.0 * RNO3; }

      // Define temp at surface soil for denitrification.
//      TM = (AVTEMP[0][dm][DY] + AVTEMP[1][dm][DY]) / 2.0;  

      TM = (0.0175/0.200) * AVTEMP[0][DY] 
           + (0.0275/0.200) * AVTEMP[1][DY]
           + (0.0450/0.200) * AVTEMP[2][DY]
           + (0.075/0.200) * AVTEMP[3][DY]
           + (0.035/0.200) * AVTEMP[4][DY];
    
      // Effect of temp on denitrifiers.
    
      TE = pow( 2.0, ((TM-22.5)/10.0) );  
    
      if ( TM > 75.0 ) { TE = ZERO; }

      EE = ZERO;
      FF = ZERO;
      II = ZERO;
      AAAA = ZERO;
      BBBB = ZERO;
      HUM = ZERO;
      INIOC = ZERO;
      INITSC = ZERO;
      INITB = ZERO;
      ININO3 = ZERO;
      FIN = ZERO;
      BACNN = ZERO;

      FERENH = 1.0;
    
      //CAS if ( TTT < Q ) { LN =(int) TT; }
      //CAS else
      //CAS {
        // Define number of soil layers involved in 
        //   denitrification.
      	
        //CASLN = Q;    
      //CAS }
    	
      LN = Q;    
    
      if( 0 == LN ) { LN = 1; }

      for( L = 0; L < Q; ++L )
      {
        if( (L+1) <= LN )
        {
          // Organic C in saturated layer L.
       
          DPOC[L] = TOC[L];   
        
          // Soluble C in saturated layer L.
        
          DPSC[L] = TCAVAI[L] * FERENH;  
        
          // Denitrifier population in sat layer L. 
         
          DPB[L] = TOC[L] * RBO * FD;  
        
          // NO3- in saturated layer L.
        
          BACN[L] = NO3[L];   
        }  
        else
        {
          DPOC[L] = ZERO;
          DPSC[L] = ZERO;
          DPB[L] = ZERO;
          BACN[L] = ZERO;
        }
      
        INIOC += DPOC[L];
        INITSC += DPSC[L];
        INITB += DPB[L];
      
        // Initial total profile active substrates.
      
        ININO3 += BACN[L];  
      }
    }
    else 
    { 
      TT = 0; 
      //cout << "NEM: TT = " << TT << endl;
    }
        
   
//CAS
//CAS	'DNDY' is now used as conditional for whether denitrification can occur
//CAS   subsequent to a rain event
//CAS
    if( DNDY > 0 && NEOM < NREP )
    {
      //cout << "NEM: Denitrification day of month = " << DNDY << endl;
      //cout << "NEM: # of days for denitrification= " << NEOM << endl;
     
      DCO2T = ZERO;
      DYN2OEMS = ZERO;
      DYN2EMS = ZERO;

      for( DHOUR = 0; DHOUR < MAXDAYHRS; ++DHOUR )
      {
        stephour( dm,
                  DNDY,
                  DHOUR,
                  pcmnt,
                  NEOM,
                  clay,
                  toppor,
                  HRWTC,
                  TT,
                  TTT,
                  LF,
                  TM,
                  Q );                  
      }


      CO2DN[DNDY] = DCO2T;              // Specify output,KgC/ha/day.
      N2ODN[DNDY] = DYN2OEMS * 1000.0;  // Specify output,gN/ha/day.
      N2DN[DNDY] = DYN2EMS * 1000.0;    // Specify output,gN/ha/day.
    
      ++NEOM;

      if( 10 == NEOM ) { DNDY = 0; }
      else { ++DNDY; }
    }

 
    IIIII += IIII;   // Totl NO3- input fro all rain events,KgN/ha.
    CSLN2O += EE;    // Totl N2O emission fro all rain events,KgN/ha.
    CSLN2 += FF;     // Totl N2 emission from all rain events.
    CSLCO2 += II;    // Totl CO2 emission fro all rain events.

    for( REC = 0; REC < LN; REC++ )
    {
      NO3[REC] = NO3P[REC];
      NO2[REC] = NO2P[REC];
    }

    ++DDAY;
  }


  // Save data for the use in future month
  for( L = 0; L < Q; ++L )           
  {
    RCVLIN[L] = RCVL[L];              
    RCLIN[L] = RCL[L];                 
    RCRIN[L] = RCR[L];                
    OCIN[L] = TOC[L];
    ANO3IN[L] = NO3[L] + NO2[L];
    ANH4IN[L] = SNH4[L];
    DPHUMIN[L] = DPHUM[L];
  }

  // Sum up the monthly output for /dndc/ common statement
  EN2ON = ZERO;
  EN2ODN = ZERO;
  EN2DN = ZERO;
  ECO2N = ZERO;
  ECO2DN = ZERO;
  ECH4N = ZERO;
//  NADD = 10;     // Number for adding missed denitri emission

 
  for( K = 0; K < mdays[dm]; ++K )
  {
    ECO2N += CO2N[K];   
    ECH4N += CH4N[K];
    EN2ON += N2ON[K];
    
    ECO2DN += CO2DN[K];
    EN2ODN += N2ODN[K]; 
    EN2DN += N2DN[K]; 
  }


//  INTOTALN = TOTALN + IIIII + WNITRIFY;      // N balance equ.
//  OUTTOTAL = LIUPTS + CSLN2O + LIN2O + CSLN2 + LINH3;
//  ALEACHNO3 = INTOTALN - OUTTOTAL - FINALSN;
//  if( ALEACHNO3 < 0.0 ) { ALEACHNO3 = 0.0; }

  return;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void NEM::uptakeMethane( const int& dy,
                         const double& wrc,
                         const double& wtoc,
                         const double& whmus,
                         double AVTEMP[CLMNLAYERS][MAXMDAYS] )
{
  double TCH4UP;

//  TCH4UP = (AVTEMP[0][dy] + AVTEMP[1][dy]) / 2.0;

  TCH4UP = (0.0175/0.200) * AVTEMP[0][dy] 
           + (0.0275/0.200) * AVTEMP[1][dy]
           + (0.0450/0.200) * AVTEMP[2][dy]
           + (0.075/0.200) * AVTEMP[3][dy]
           + (0.035/0.200) * AVTEMP[4][dy];
      
  // Daily CH4 uptake by soil.
      
  CH4UP = 0.0012 * TCH4UP 
          * (36.26 * log( wrc+wtoc+whmus)/2.3026 - 162.57) 
          / 25.0; 
     
  if( TCH4UP < ZERO ) { CH4UP = ZERO; }
  
  if( CH4UP < ZERO ) { CH4UP = ZERO; }
  
  LICH4UP += CH4UP;  // CH4 uptake in a year,KgC/ha.

};
