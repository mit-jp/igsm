/* **************************************************************
*****************************************************************
MITELM44C.CPP - Runs TEM for a single grid cell

Modifications:

20060129 - DWK created by modifying telm437.cpp
20060129 - DWK changed include from telm437.h to mitelm437a.h
20060129 - DWK changed TEMelmnt43:: to MITelmnt43::
20060202 - DWK added I_NIRR, I_PAR, I_TAIR, I_PREC, I_RAIN,   
           I_SNWFAL, I_CO2, I_AOT40 to outputTEMmonth() and 
           temwritepred()
20060202 - DWK deleted atmswritemiss() and atmswritepred()
20060422 - DWK added MDMnpp to initializeCohortTEMState(),
           readCohortState(), saveTEMCohortState(), 
           setCohortTEMState(), and writeCohortState()
20070129 - DWK changed include mitelm438.h to mitelm44.h
20070129 - DWK changed MITelmnt43:: to MITelmnt44::
20080130 - DWK changed include from mitelm44.h to mitelm44a.h
20080827 - DWK changed include from mitelm44a.h to mitelm44b.h
20110705 - DWK changed include from mitelm44b.h to mitelm44c.h
20110911 - DWK added functions createCohortProducts() and 
           updateChangedCohort() 
20140519 - DWK added "yr" variables to getCohortTEMstate(),
           initializeCohortTEMState(), readCohortState(), 
           saveTEMCohortState(), setCohortTEMState(), and 
           writeCohortState()
20140609 - DWK added  readBinaryTEMCohortState() and
           writeBinaryCohortState()
                                                                 
****************************************************************
************************************************************* */

#include<cstdio>

  using std::fscanf;
  using std::FILE;
  using std::printf;

#include<iostream>

  using std::cout;
  using std::cin;
  using std::ios;
  using std::cerr;
  using std::endl;

#include<fstream>

  using std::ifstream;
  using std::ofstream;

#include<iomanip>

  using std::setprecision;
  using std::setw;

#include<cstdlib>

  using std::exit;
  using std::atof;
  using std::atoi;

#include<cmath>

  using std::exp;
  using std::fabs;
  using std::pow;

#include<vector>

  using std::vector;
      
#include<string>
  
  using std::string;


#include "mitelm44c.h"

/* *************************************************************
************************************************************* */

MITelmnt44::MITelmnt44()
{

  col = MISSING;
  row = MISSING;
  carea = -999;
  subarea = -999;
  fatalerr = 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */
void MITelmnt44::createCohortProducts( const int& pichrt,
                                       const double& pproparea,
                                       const Biomass& pformPROD10,
                                       const Biomass& pformPROD100,
                                       const Biomass& pslash,
                                       const Biomass& pvconvrtflx,
                                       const Biomass& psconvrtflx,
                                       const double& pnvretent,
                                       const double& pnsretent )
{
  int i10;
  int i100;

  i10 = tem.ag.getPRODUCTYEAR() % 10;
  i100 = tem.ag.getPRODUCTYEAR() % 100;

  if( cohort[pichrt].chrtarea > 0 )
  {
    cohort[pichrt].formPROD10.carbon = pformPROD10.carbon
                                       * pproparea
                                       / ((double) cohort[pichrt].chrtarea);

    cohort[pichrt].formPROD10.nitrogen = pformPROD10.nitrogen
                                         * pproparea
                                         / ((double) cohort[pichrt].chrtarea);

    cohort[pichrt].formPROD100.carbon = pformPROD100.carbon
                                        * pproparea
                                        / ((double) cohort[pichrt].chrtarea);

    cohort[pichrt].formPROD100.nitrogen = pformPROD100.nitrogen
                                          * pproparea
                                          / ((double) cohort[pichrt].chrtarea);

    cohort[pichrt].slash.carbon = pslash.carbon
                                  * pproparea
                                  / (((double) cohort[pichrt].chrtarea)
                                  * (double) CYCLE);

    cohort[pichrt].slash.nitrogen = pslash.nitrogen
                                    * pproparea
                                    / (((double) cohort[pichrt].chrtarea)
                                    * (double) CYCLE);

    cohort[pichrt].vconvrtflx.carbon = pvconvrtflx.carbon
                                       * pproparea
                                       / (((double) cohort[pichrt].chrtarea)
                                       * (double) CYCLE);

    if( cohort[pichrt].vconvrtflx.carbon < ZERO )
    {
      cout << " lat = " << row;
      cout << " ichrt = " << pichrt;
      cout << " vconvrtc = " << cohort[pichrt].vconvrtflx.carbon;
      cout << endl;

      exit( -1 );
    }

    cohort[pichrt].sconvrtflx.carbon = psconvrtflx.carbon
                                       * pproparea
                                       / (((double) cohort[pichrt].chrtarea)
                                       * (double) CYCLE);

    if( cohort[pichrt].sconvrtflx.carbon < ZERO )
    {
      cout << " lat = " << row;
      cout << " ichrt = " << pichrt;
      cout << " sconvrtc = " << cohort[pichrt].sconvrtflx.carbon;
      cout << endl;

      exit( -1 );
    }

    cohort[pichrt].vconvrtflx.nitrogen = pvconvrtflx.nitrogen
                                         * pproparea
                                         / (((double) cohort[pichrt].chrtarea)
                                         * (double) CYCLE);
                                
    cohort[pichrt].sconvrtflx.nitrogen = psconvrtflx.nitrogen
                                         * pproparea
                                         / (((double) cohort[pichrt].chrtarea)
                                         * (double) CYCLE);

    cohort[pichrt].nvretent = pnvretent
                              * pproparea
                              / (((double) cohort[pichrt].chrtarea)
                              * (double) CYCLE);

    cohort[pichrt].nsretent = pnsretent
                              * pproparea
                              / (((double) cohort[pichrt].chrtarea)
                              * (double) CYCLE);
  }
 
  cohort[pichrt].initPROD10[i10].carbon = cohort[pichrt].formPROD10.carbon;

  if( cohort[pichrt].prevPROD10.carbon < ZERO )
  {
    cohort[pichrt].prevPROD10.carbon = ZERO;
  }     


  cohort[pichrt].initPROD10[i10].nitrogen = cohort[pichrt].formPROD10.nitrogen;

  if( cohort[pichrt].prevPROD10.nitrogen < ZERO )
  {
    cohort[pichrt].prevPROD10.nitrogen = ZERO;
  }     


  cohort[pichrt].initPROD100[i100].carbon =  cohort[pichrt].formPROD100.carbon;

  if( cohort[pichrt].prevPROD100.carbon < ZERO )
  {
    cohort[pichrt].prevPROD100.carbon = ZERO;
  }     
                                               

  cohort[pichrt].initPROD100[i100].nitrogen = cohort[pichrt].formPROD100.nitrogen;

  if( cohort[pichrt].prevPROD100.nitrogen < ZERO )
  {
    cohort[pichrt].prevPROD100.nitrogen = ZERO;
  }     
 
  cohort[pichrt].convrtflx.carbon = cohort[pichrt].vconvrtflx.carbon
                                    + cohort[pichrt].sconvrtflx.carbon;

  cohort[pichrt].convrtflx.nitrogen = cohort[pichrt].vconvrtflx.nitrogen
                                      + cohort[pichrt].sconvrtflx.nitrogen;

  cohort[pichrt].nretent = cohort[pichrt].nvretent
                           + cohort[pichrt].nsretent;

};

/* *************************************************************
************************************************************* */

int MITelmnt44::coregerr( ofstream& rflog1,
                          const string& varname1,
                          const float& col1,
                          const float& row1,
                          const string& varname2,
                          const float& col2,
                          const float& row2 )
{

  int fatalerr = 0;

  if( col1 != col2 || row1 != row2 )
  {
    fatalerr = 1;

    cout << "ERROR:  " << varname1 << " data and ";
    cout << varname2 << "data are not coregistered." << endl;
    cout << "COL = " << col1 << " and ROW = " << row1;
    cout << " in " << varname1 << " data" << endl;
    cout << "COL = " << col2 << " and ROW = " << row2;
    cout << " in " << varname2 << " data" << endl;

    rflog1 << "ERROR:  " << varname1 << " data and ";
    rflog1 << varname2 << "data are not coregistered." << endl;
    rflog1 << "COL = " << col1 << " and ROW = " << row1;
    rflog1 << " in " << varname1 << " data" << endl;
    rflog1 << "COL = " << col2 << " and ROW = " << row2;
    rflog1 << " in " << varname2 << " data" << endl;
  }

  return fatalerr;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

int MITelmnt44::equilibrateTEM( const int& pchrt,
                                const double& ptol )
{
  // Run TEM until steady state conditions occur (equilibrium)

  int dyr = 0;
  int dm;

/*  if( 14 == (pchrt+1) )
  {
    tem.dbug = 1;
    tem.veg.dbug = 1;
  }
  else
  {
    tem.dbug = 0;
    tem.veg.dbug = 0;
  }
*/ 

  // Initialize standing stocks of carbon and nitrogen from 
  //   calibration ("ECD") data and fluxes for integrator

  tem.ECDsetODEstate( tem.veg.cmnt, tem.soil.getPSIPLUSC() );


  // Set previous value of TEM ODE state variables to the 
  //   current values of TEM ODE state variables for initial
  //   conditions
 
  tem.setPrevState();


  // Initialize all disturbance fluxes to zero

  tem.ag.resetMonthlyDisturbFluxes();

  tem.ag.setFIRENDEP( ZERO );


  // Initialize all agricultural and wood product pools to 
  //   zero

  tem.ag.resetPROD();
  
  tem.totyr = 0;
  tem.endeq = 0;
  tem.intflag = 0;
  tem.initFlag = 0;
  tem.microbe.nem.firstday = 0;

  // Initialize tem.atms.prevtair for water balance (WBM) and 
  //   soil thermal (STM) calculations
  
  tem.atms.setPREVTAIR( tem.atms.getTAIR() );


  // Initialize tem.atms.prev2tair for water balance  
  //   calculations (WBM)

  tem.atms.setPREV2TAIR( tem.atms.getTAIR() );

  
  // Initialize tem.veg.prvleafmx and tem.veg.prevunrmleaf
  //   for phenology calculations
  
  tem.veg.setPRVLEAFMX( tem.veg.getINITLEAFMX( tem.veg.cmnt ) );
  
  tem.veg.setPREVUNRMLEAF( tem.veg.getUNLEAF12( tem.veg.cmnt ) );


  while( (dyr < tem.runsize) && (tem.endeq < 2) )
  {
    for( dm = 0; dm < CYCLE; ++dm )
    {
      // Pass climate data for particular month to TEM

      tem.atms.setNIRR( climate[mitclm.I_NIRR][dm] );

//      cout << climate[mitclm.I_NIRR][dm] << endl;

      tem.atms.setPAR( climate[mitclm.I_PAR][dm] );

//      cout << climate[mitclm._NIRR][dm] << endl;

      tem.atms.setTAIR( climate[mitclm.I_TAIR][dm] );

//      cout << climate[mitclm.I_NIRR][dm] << endl;

      tem.atms.setPET( initPET[pchrt][dm] );

//      cout << climate[mitclm.I_NIRR][dm] << endl;

      tem.soil.setEET( initAET[pchrt][dm] );

//      cout << climate[mitclm.I_NIRR][dm] << endl;

      tem.soil.setMOIST( initSH2O[pchrt][dm] );

//      cout << climate[mitclm.I_NIRR][dm] << endl;

      tem.soil.setSNOWPACK( initSNOWPACK[pchrt][dm] );

//      cout << climate[mitclm.I_NIRR][dm] << endl;

      tem.soil.setSURFRUN( initSURFRUN[pchrt][dm] );

//      cout << climate[mitclm.I_NIRR][dm] << endl;

      tem.soil.setDRAINAGE( initDRAINAGE[pchrt][dm] );

//      cout << climate[mitclm.I_NIRR][dm] << endl;

      tem.atms.setCO2( climate[mitclm.I_CO2][dm] );

//      cout << climate[mitclm.I_NIRR][dm] << endl;

      tem.atms.setAOT40( climate[mitclm.I_AOT40][dm] );

//      cout << climate[mitclm.I_AOT40][dm] << endl;


      tem.endeq = tem.stepmonth( dyr, 
                                 dm, 
                                 tem.intflag, 
                                 ptol );


      // Save TEM output to telmnt[0].output
      
      outputTEMmonth( pchrt, dm );
    }

//    exit( -1 );

    ++dyr;

    ++tem.totyr;


// Check to see if steady state conditions have been reached.

    if( dyr >= tem.strteq && 0 == tem.endeq )
    {
      tem.endeq = tem.testEquilibrium();
    }
  }

  if( tem.totyr >= tem.runsize && tem.endeq < 2 ) 
  { 
    tem.nattempt += 1; 

    tem.initFlag = 0;
  }
  else 
  { 
    tem.initFlag = 1; 
  }

  return tem.nattempt;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void MITelmnt44::getTEMCohortState( const int& pichrt )
{
  int dlyr; // NEM
  int dm;
  int i;
  
// cohort[pichrt].srcCohort = cohort[pichrt].srcCohort

// cohort[pichrt].standage = cohort[pichrt].standage

// cohort[pichrt].chrtarea = cohort[pichrt].chrtarea

// cohort[pichrt].prevchrtarea = cohort[pichrt].prevchrtarea

  tem.veg.setPOTVEG( cohort[pichrt].potveg );

  tem.veg.setCURRENTVEG( cohort[pichrt].currentveg );

  tem.veg.setSUBTYPE( cohort[pichrt].subtype );

  tem.veg.cmnt = cohort[pichrt].cmnt;

  for( i = 0; i < MAXSTATE; ++i )
  {
    tem.setY( cohort[pichrt].y[i], i );

    tem.setPREVY( cohort[pichrt].prevy[i], i );
  }

  tem.ag.cmnt = cohort[pichrt].agcmnt;

  tem.ag.setGROWDD( cohort[pichrt].aggrowdd );

  tem.ag.setKD( cohort[pichrt].agkd );

  tem.ag.prvstate = cohort[pichrt].agprvstate;

  tem.ag.state = cohort[pichrt].agstate;

  tem.veg.setC2N( cohort[pichrt].c2n );
  
  tem.veg.setCNEVEN( cohort[pichrt].cneven );

  tem.ag.setCONVRTFLXC( cohort[pichrt].convrtflx.carbon );
  tem.ag.setCONVRTFLXN( cohort[pichrt].convrtflx.nitrogen );

  tem.ag.setCROPPRVEETMX( cohort[pichrt].cropprveetmx );

  tem.ag.setCROPPRVLEAFMX( cohort[pichrt].cropprvleafmx );

  tem.ag.setCROPPRVPETMX( cohort[pichrt].cropprvpetmx );

  tem.ag.setCROPRESIDUEC( cohort[pichrt].cropResidue.carbon );
  tem.ag.setCROPRESIDUEN( cohort[pichrt].cropResidue.nitrogen );

  tem.ag.setCROPTOPT( cohort[pichrt].croptopt );

  tem.distmnthcnt = cohort[pichrt].distmnthcnt;
  
  tem.disturbflag = cohort[pichrt].disturbflag;

  tem.disturbmonth = cohort[pichrt].disturbmonth;
  
  tem.soil.setEETMX( cohort[pichrt].eetmx );

  tem.ag.fertflag = cohort[pichrt].fertflag;                              

  tem.firemnthcnt = cohort[pichrt].firemnthcnt; 
  
  tem.ag.setFIRENDEP( cohort[pichrt].firendep );

  tem.ag.setFORMPROD10C( cohort[pichrt].formPROD10.carbon );
  tem.ag.setFORMPROD10N( cohort[pichrt].formPROD10.nitrogen );

  tem.ag.setFORMPROD100C( cohort[pichrt].formPROD100.carbon );
  tem.ag.setFORMPROD100N( cohort[pichrt].formPROD100.nitrogen );

  tem.veg.setFPREVOZONE( cohort[pichrt].fprevozone );
  
  tem.ag.setFRI( cohort[pichrt].FRI );

  for( dm = 0; dm < CYCLE; ++dm )
  {
    tem.ag.setINITPROD1C( cohort[pichrt].initPROD1[dm].carbon,
                          dm );

    tem.ag.setINITPROD1N( cohort[pichrt].initPROD1[dm].nitrogen,
                          dm ); 
  }
  
  for( i = 0; i < 10; ++i )
  {
    tem.ag.setINITPROD10C( cohort[pichrt].initPROD10[i].carbon, 
                           i );

    tem.ag.setINITPROD10N( cohort[pichrt].initPROD10[i].nitrogen, 
                           i );
  }
    
  for( i = 0; i < 100; ++i )
  {
    tem.ag.setINITPROD100C( cohort[pichrt].initPROD100[i].carbon, 
                            i );

    tem.ag.setINITPROD100N( cohort[pichrt].initPROD100[i].nitrogen, 
                            i );
  }

  tem.ag.irrgflag = cohort[pichrt].irrgflag;                              
  
  tem.microbe.setKD( cohort[pichrt].kd );

// cohort[pichrt].MDMnpp = cohort[pichrt].MDMnpp

  tem.ag.setNATPRVEETMX( cohort[pichrt].natprveetmx );

  tem.ag.setNATPRVLEAFMX( cohort[pichrt].natprvleafmx );

  tem.ag.setNATPRVPETMX( cohort[pichrt].natprvpetmx );

  tem.ag.setNATSEEDC( cohort[pichrt].natseedC );

  tem.ag.setNATSEEDSTRN( cohort[pichrt].natseedSTRN );

  tem.ag.setNATSEEDSTON( cohort[pichrt].natseedSTON );

  tem.ag.setNATSOIL( cohort[pichrt].natsoil );

  tem.ag.setNATTOPT( cohort[pichrt].nattopt );

  tem.ag.setNATYREET( cohort[pichrt].natyreet );

  tem.ag.setNATYRPET( cohort[pichrt].natyrpet );


  // *************************** NEM ***************************
  
  for( dlyr = 0; dlyr < NLVL; ++dlyr )
  {
    tem.microbe.nem.setANH4IN( cohort[pichrt].NEManh4in[dlyr], 
                               dlyr ); 
  
    tem.microbe.nem.setANO3IN( cohort[pichrt].NEMano3in[dlyr],
                               dlyr );

    tem.microbe.nem.setDPHUMIN( cohort[pichrt].NEMdphumin[dlyr],
                                dlyr ); 

    tem.microbe.nem.setOCIN( cohort[pichrt].NEMocin[dlyr],
                             dlyr );

    tem.microbe.nem.setRCLIN( cohort[pichrt].NEMrclin[dlyr],
                              dlyr ); 
  
    tem.microbe.nem.setRCRIN( cohort[pichrt].NEMrcrin[dlyr],
                              dlyr );
  
    tem.microbe.nem.setRCVLIN( cohort[pichrt].NEMrcvlin[dlyr],
                               dlyr );
  }
  
  tem.microbe.nem.setNONSOLC( cohort[pichrt].NEMnsolc );

//  tem.microbe.nem.setTOPDENS( cohort[pichrt].NEMtopdens );
  
//  tem.microbe.nem.setTOPKSAT( cohort[pichrt].NEMtopksat );
  
//  tem.microbe.nem.setTOPPOR( cohort[pichrt].NEMtoppor );

  // ***********************************************************

  
  tem.veg.setNEWLEAFMX( cohort[pichrt].newleafmx );

  tem.veg.setNEWTOPT( cohort[pichrt].newtopt );

  tem.ag.setNRETENT( cohort[pichrt].nretent );

  tem.ag.setNSRETENT( cohort[pichrt].nsretent );

  tem.ag.setNVRETENT( cohort[pichrt].nvretent );

  tem.atms.setPETMX( cohort[pichrt].petmx );

  tem.atms.setPREV2TAIR( cohort[pichrt].prev2tair );

  tem.atms.setPREVCO2( cohort[pichrt].prevco2 );

  tem.ag.setPREVCROPRESIDUEC( cohort[pichrt].prevCropResidue.carbon );
  tem.ag.setPREVCROPRESIDUEN( cohort[pichrt].prevCropResidue.nitrogen );

  tem.ag.setPREVPROD1C( cohort[pichrt].prevPROD1.carbon );
  tem.ag.setPREVPROD1N( cohort[pichrt].prevPROD1.nitrogen );

  tem.ag.setPREVPROD10C( cohort[pichrt].prevPROD10.carbon );
  tem.ag.setPREVPROD10N( cohort[pichrt].prevPROD10.nitrogen );

  tem.ag.setPREVPROD100C( cohort[pichrt].prevPROD100.carbon );
  tem.ag.setPREVPROD100N( cohort[pichrt].prevPROD100.nitrogen );

//  tem.soil.setPREVSPACK( cohort[pichrt].prevspack );

  tem.atms.setPREVTAIR( cohort[pichrt].prevtair );

  tem.veg.setPREVUNRMLEAF( cohort[pichrt].prevunrmleaf );

//  tem.ag.setPROD10PAR( cohort[pichrt].prod10par ); 

//  tem.ag.setPROD100PAR( cohort[pichrt].prod100par ); 

  tem.ag.setPRODUCTYEAR( cohort[pichrt].productYear );
  
  tem.ag.setPRVCROPNPP( cohort[pichrt].prvcropnpp );

  tem.soil.setPRVEETMX( cohort[pichrt].prveetmx );

  tem.veg.setPRVLEAFMX( cohort[pichrt].prvleafmx );

  tem.atms.setPRVPETMX( cohort[pichrt].prvpetmx );

// cohort[pichrt].qc = cohort[pichrt].qc

//  tem.ag.setSCONVERT( cohort[pichrt].sconvert ); 
  
  tem.ag.setSCONVRTFLXC( cohort[pichrt].sconvrtflx.carbon );
  tem.ag.setSCONVRTFLXN( cohort[pichrt].sconvrtflx.nitrogen );

  tem.ag.setSLASHC( cohort[pichrt].slash.carbon );
  tem.ag.setSLASHN( cohort[pichrt].slash.nitrogen );

//  tem.ag.setSLASHPAR( cohort[pichrt].slashpar ); 
   
  tem.ag.tillflag = cohort[pichrt].tillflag;                           

  tem.veg.setTOPT( cohort[pichrt].topt );

// cohort[pichrt].tqc = cohort[pichrt].tqc 

//  tem.ag.setVCONVERT( cohort[pichrt].vconvert ); 

  tem.ag.setVCONVRTFLXC( cohort[pichrt].vconvrtflx.carbon );
  tem.ag.setVCONVRTFLXN( cohort[pichrt].vconvrtflx.nitrogen );

//  tem.ag.setVRESPAR( cohort[pichrt].vrespar ); 


  tem.ag.setYRSTUBC( cohort[pichrt].yragstubC );
  
  tem.ag.setYRSTUBN( cohort[pichrt].yragstubN );

  tem.ag.setYRCFLUX( cohort[pichrt].yrcflux );
  
  tem.soil.setYRCH4CSMP( cohort[pichrt].yrCH4csmp );
  
  tem.soil.setYRCH4EMISSION( cohort[pichrt].yrCH4ems );
  
  tem.soil.setYRCH4FLUX( cohort[pichrt].yrCH4flx );
  
  tem.soil.setYRCO2DENTRFLUX( cohort[pichrt].yrCO2dnflx );
  
  tem.soil.setYRCO2NTRFLUX( cohort[pichrt].yrCO2nflx );
  
  tem.ag.setYRCONVRTC( cohort[pichrt].yrconvrtC );
  
  tem.ag.setYRCONVRTN( cohort[pichrt].yrconvrtN );

  tem.ag.setYRDECAYPROD1C( cohort[pichrt].yrdecayPROD1C );
  
  tem.ag.setYRDECAYPROD10C( cohort[pichrt].yrdecayPROD10C );
  
  tem.ag.setYRDECAYPROD100C( cohort[pichrt].yrdecayPROD100C );

  tem.ag.setYRDECAYPROD1N( cohort[pichrt].yrdecayPROD1N );
  
  tem.ag.setYRDECAYPROD10N( cohort[pichrt].yrdecayPROD10N );
  
  tem.ag.setYRDECAYPROD100N( cohort[pichrt].yrdecayPROD100N );
  
  tem.ag.setYRDECAYTOTPRODC( cohort[pichrt].yrdecayTOTPRODC );

  tem.ag.setYRDECAYTOTPRODN( cohort[pichrt].yrdecayTOTPRODN );
 
  tem.soil.setYREET( cohort[pichrt].yreet );
  
  tem.ag.setYRFERTN( cohort[pichrt].yrfertn );

  tem.ag.setYRFLUXRESIDUEC( cohort[pichrt].yrfluxResidueC );
  
  tem.ag.setYRFLUXRESIDUEN( cohort[pichrt].yrfluxResidueN );
  
  tem.ag.setYRFORMPROD1C( cohort[pichrt].yrformPROD1C );
  
  tem.ag.setYRFORMPROD10C( cohort[pichrt].yrformPROD10C );
  
  tem.ag.setYRFORMPROD100C( cohort[pichrt].yrformPROD100C );

  tem.ag.setYRFORMPROD1N( cohort[pichrt].yrformPROD1N );
  
  tem.ag.setYRFORMPROD10N( cohort[pichrt].yrformPROD10N );
  
  tem.ag.setYRFORMPROD100N( cohort[pichrt].yrformPROD100N );
  
  tem.ag.setYRFORMRESIDUEC( cohort[pichrt].yrformResidueC );
  
  tem.ag.setYRFORMRESIDUEN( cohort[pichrt].yrformResidueN );

  tem.ag.setYRFORMTOTPRODC( cohort[pichrt].yrformTOTPRODC );

  tem.ag.setYRFORMTOTPRODN( cohort[pichrt].yrformTOTPRODN );

  tem.veg.setYRFPC( cohort[pichrt].yrfpc );
  
  tem.veg.setYRGPP( cohort[pichrt].yrgpp );
  
  tem.veg.setYRGPR( cohort[pichrt].yrgpr );
  
//  tem.soil.setYRH2OYIELD( cohort[pichrt].yrH2Oyld );
  
  tem.microbe.setYRNUPTAKE( cohort[pichrt].yrimmob );
  
//  tem.soil.setYRINEET( cohort[pichrt].yrineet );

  tem.veg.setYRINGPP( cohort[pichrt].yringpp );
  
  tem.veg.setYRINNPP( cohort[pichrt].yrinnpp );
  
//  tem.ag.setYRIRRIG( cohort[pichrt].yrirrig );
  
  tem.veg.setYRLAI( cohort[pichrt].yrlai );
    
  tem.veg.setYRLEAF( cohort[pichrt].yrleaf );

  tem.veg.setYRLTRFALC( cohort[pichrt].yrltrfalC );

  tem.veg.setYRLTRFALN( cohort[pichrt].yrltrfalN );  	
 
  tem.soil.setYRN2FLUX( cohort[pichrt].yrN2flx );
  
  tem.soil.setYRN2ODENTRFLUX( cohort[pichrt].yrN2Odnflx );
  
  tem.soil.setYRN2OFLUX( cohort[pichrt].yrN2Oflx );
  
  tem.soil.setYRN2ONTRFLUX( cohort[pichrt].yrN2Onflx );
  
  tem.setYRNCE( cohort[pichrt].yrnce );

  tem.setYRNEP( cohort[pichrt].yrnep );

  tem.soil.setYRNINPUT( cohort[pichrt].yrninput );

  tem.soil.setYRNLOST( cohort[pichrt].yrnlost );
  
  tem.microbe.setYRNMIN( cohort[pichrt].yrnmin );

  tem.veg.setYRNPP( cohort[pichrt].yrnpp );
  
  tem.ag.setYRNRENT( cohort[pichrt].yrnrent ); 
  
  tem.ag.setYRNSRENT( cohort[pichrt].yrnsrent );
  
  tem.ag.setYRNVRENT( cohort[pichrt].yrnvrent );
  
  tem.atms.setYRPET( cohort[pichrt].yrpet );
  
//  tem.atms.setYRRGRNDH2)( cohort[pichrt].yrrgrndH2O );

  tem.microbe.setYRRH( cohort[pichrt].yrrh );
   
//  tem.atms.setYRRAIN( cohort[pichrt].yrrain );
  
//  tem.soil.setYRRPERC( cohort[pichrt].yrrperc );
  
//  tem.soil.setYRRRUN( cohort[pichrt].yrrrun );

  tem.ag.setYRSCONVRTC( cohort[pichrt].yrsconvrtC );

  tem.ag.setYRCONVRTN( cohort[pichrt].yrsconvrtN );
  
//  tem.soil.setYRSGRNDH2O( cohort[pichrt].yrsgrndH2O );

  tem.ag.setYRSLASHC( cohort[pichrt].yrslashC );

  tem.ag.setYRSLASHN( cohort[pichrt].yrslashN );
  
//  tem.atms.setYRSNOWFALL( cohort[pichrt].yrsnowfall );

//  tem.soil.setYRSNOWINF( cohort[pichrt].yrsnowinf );

//  tem.soil.setYRSNOWPACK( cohort[pichrt].yrsnowpack );

//  tem.soil.setYRAVLH2O( cohort[pichrt].yrsoilavlH2O );
  
  tem.soil.setYRAVLN( cohort[pichrt].yrsoilavlN );

//  tem.soil.setYRC2N( cohort[pichrt].yrsoilC2N );
   
//  tem.soil.setYRSMOIST( cohort[pichrt].yrsoilmoist );
  
  tem.soil.setYRORGC( cohort[pichrt].yrsoilorgC );
  
  tem.soil.setYRORGN( cohort[pichrt].yrsoilorgN );

//  tem.soil.setYRPCTP( cohort[pichrt].yrsoilpctp );
  
//  tem.soil.setYRVSM( cohort[pichrt].yrsoilvsm );

//  tem.soil.setYRSPERC( cohort[pichrt].yrsperc );
  
//  tem.soil.setSRUN( cohort[pichrt].yrsrun );

  tem.setYRTOTALC( cohort[pichrt].yrtotalC );
  
  tem.veg.setYRUNLEAF( cohort[pichrt].yrunleaf );

  tem.ag.setYRVCONVRTC( cohort[pichrt].yrvconvrtC );
  
  tem.ag.setYRCONVRTN( cohort[pichrt].yrvconvrtN );
 
  tem.veg.setYRVEGC( cohort[pichrt].yrvegC );

//  tem.veg.setYRVEGC2N( cohort[pichrt].yrvegC2N );
    
  tem.veg.setYRINNUP( cohort[pichrt].yrveginnup );

 tem.veg.setYRLUP( cohort[pichrt].yrveglup );

  tem.veg.setYRVEGN( cohort[pichrt].yrvegN );

  tem.veg.setYRNMOBIL( cohort[pichrt].yrvegnmobil );

  tem.veg.setYRNRSORB( cohort[pichrt].yrvegnrsorb );

  tem.veg.setYRNUPTAKE( cohort[pichrt].yrvegnup );

  tem.veg.setYRRGROWTH( cohort[pichrt].yrvegrgrowth );

  tem.veg.setYRRMAINT( cohort[pichrt].yrvegrmaint );

  tem.veg.setYRSTOREN( cohort[pichrt].yrvegSTON );
  
  tem.veg.setYRSTRUCTN( cohort[pichrt].yrvegSTRN );

  tem.veg.setYRSUP( cohort[pichrt].yrvegsup );
	
};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void MITelmnt44::initializeCohortTEMState( const int& pichrt )
{
  int dlyr; // NEM
  int dm;
  int i;

// cohort[pichrt].srcCohort initialized by LCLUC

// cohort[pichrt].standage initialized by LCLUC

// cohort[pichrt].chrtarea initialized by LCLUC

  cohort[pichrt].prevchrtarea = cohort[pichrt].chrtarea;

//  cohort[pichrt].potveg initialized by LCLUC

//  cohort[pichrt].currentveg initialized by LCLUC

//  cohort[pichrt].subtype initialized by LCLUC

// cohort[pichrt].cmnt initialized by LCLUC

  for( i = 0; i < MAXSTATE; ++i )
  {
    cohort[pichrt].y[i]  = ZERO;
    cohort[pichrt].prevy[i]  = ZERO;
  }

//  cohort[pichrt].agcmnt initialized by LCLUC

  cohort[pichrt].aggrowdd  = ZERO;

  cohort[pichrt].agkd  = ZERO;

//  cohort[pichrt].agprvstate initialized by LCLUC

//  cohort[pichrt].agstate initialized by LCLUC

  cohort[pichrt].c2n  = ZERO;
  
  cohort[pichrt].cneven  = ZERO;
  
  cohort[pichrt].convrtflx.carbon  = ZERO;
  cohort[pichrt].convrtflx.nitrogen  = ZERO;

  cohort[pichrt].cropprveetmx  = ZERO;

  cohort[pichrt].cropprvleafmx  = ZERO;

  cohort[pichrt].cropprvpetmx  = ZERO;

  cohort[pichrt].cropResidue.carbon  = ZERO;
  cohort[pichrt].cropResidue.nitrogen  = ZERO;

  cohort[pichrt].croptopt  = ZERO;

  cohort[pichrt].distmnthcnt = 0;

//  cohort[pichrt].disturbflag initialized by LCLUC
  cohort[pichrt].disturbflag = 0;
  
//  cohort[pichrt].disturbmonth initialized by LCLUC 
  cohort[pichrt].disturbmonth = 0;
  
  cohort[pichrt].eetmx  = ZERO;

//  cohort[pichrt].fertflag initialized by LCLUC

  cohort[pichrt].firemnthcnt = 0;
  
  cohort[pichrt].firendep  = ZERO;
                                                      
  cohort[pichrt].formPROD10.carbon  = ZERO;
  cohort[pichrt].formPROD10.nitrogen  = ZERO;
  
  cohort[pichrt].formPROD100.carbon  = ZERO;
  cohort[pichrt].formPROD100.nitrogen  = ZERO;

  cohort[pichrt].fprevozone  = ZERO;
     
//  cohort[pichrt].FRI initialized by LCLUC

  for( dm = 0; dm < CYCLE; ++dm )
  {
    cohort[pichrt].initPROD1[dm].carbon = ZERO;
    cohort[pichrt].initPROD1[dm].nitrogen = ZERO;
  }
  
  for( i = 0; i < 10; ++i )
  {
    cohort[pichrt].initPROD10[i].carbon = ZERO; 
    cohort[pichrt].initPROD10[i].nitrogen = ZERO;
  }
    
  for( i = 0; i < 100; ++i )
  {
    cohort[pichrt].initPROD100[i].carbon = ZERO;
    cohort[pichrt].initPROD100[i].nitrogen = ZERO;
  }                            

//  cohort[pichrt].irrgflag initialized by LCLUC
  
  cohort[pichrt].kd  = ZERO;

  cohort[pichrt].MDMnpp  = ZERO;     // MDM

  cohort[pichrt].natprveetmx  = ZERO;

  cohort[pichrt].natprvleafmx  = ZERO;

  cohort[pichrt].natprvpetmx  = ZERO;

  cohort[pichrt].natseedC  = ZERO;

  cohort[pichrt].natseedSTRN  = ZERO;

  cohort[pichrt].natseedSTON  = ZERO;

  cohort[pichrt].natsoil  = ZERO;

  cohort[pichrt].nattopt  = ZERO;

  cohort[pichrt].natyreet  = ZERO;

  cohort[pichrt].natyrpet  = ZERO;


  // ********************* NEM *********************************
  
  for( dlyr = 0; dlyr < NLVL; ++dlyr )
  {
    cohort[pichrt].NEManh4in[dlyr]  = ZERO; 
  
    cohort[pichrt].NEMano3in[dlyr]  = ZERO;

    cohort[pichrt].NEMdphumin[dlyr]  = ZERO; 

    cohort[pichrt].NEMocin[dlyr]  = ZERO;

    cohort[pichrt].NEMrclin[dlyr]  = ZERO; 
  
    cohort[pichrt].NEMrcrin[dlyr]  = ZERO;
  
    cohort[pichrt].NEMrcvlin[dlyr]  = ZERO;
  }
  
  cohort[pichrt].NEMnsolc  = ZERO;

  // cohort[pichrt].NEMtopdens prescribed from soil layer file
  
  // cohort[pichrt].NEMtopksat prescribed from soil layer file
  
  // cohort[pichrt].NEMtoppor prescribed from soil layer file

  // ***********************************************************


  cohort[pichrt].newleafmx  = ZERO;

  cohort[pichrt].newtopt  = ZERO;

  cohort[pichrt].nretent  = ZERO;

  cohort[pichrt].nsretent  = ZERO;

  cohort[pichrt].nvretent  = ZERO;

  cohort[pichrt].petmx  = ZERO;

  cohort[pichrt].prev2tair  = ZERO;

  cohort[pichrt].prevco2  = ZERO;

  cohort[pichrt].prevCropResidue.carbon  = ZERO;
  cohort[pichrt].prevCropResidue.nitrogen  = ZERO;
  
  cohort[pichrt].prevPROD1.carbon  = ZERO; 
  cohort[pichrt].prevPROD1.nitrogen  = ZERO;
  
  cohort[pichrt].prevPROD10.carbon  = ZERO; 
  cohort[pichrt].prevPROD10.nitrogen  = ZERO;
  
  cohort[pichrt].prevPROD100.carbon  = ZERO;
  cohort[pichrt].prevPROD100.nitrogen  = ZERO;

//  cohort[pichrt].prevspack  = ZERO;

  cohort[pichrt].prevtair  = ZERO;
  
  cohort[pichrt].prevunrmleaf  = ZERO;

// cohort[pichrt].prod10par initialized by LCLUC

// cohort[pichrt].prod100par initialized by LCLUC

  cohort[pichrt].productYear = -99;
  
  cohort[pichrt].prvcropnpp  = ZERO;

  cohort[pichrt].prveetmx  = ZERO;

  cohort[pichrt].prvleafmx  = ZERO;

  cohort[pichrt].prvpetmx  = ZERO;

  cohort[pichrt].qc = 0;

// cohort[pichrt].sconvert initialized by LCLUC

  cohort[pichrt].sconvrtflx.carbon  = ZERO;
  cohort[pichrt].sconvrtflx.nitrogen  = ZERO;

  cohort[pichrt].slash.carbon  = ZERO;
  cohort[pichrt].slash.nitrogen  = ZERO;                       

// cohort[pichrt].slashpar initialized by LCLUC

// cohort[pichrt].tillflag initialized by LCLUC

  cohort[pichrt].topt  = ZERO;

  cohort[pichrt].tqc = 0;

// cohort[pichrt].vconvert initialized by LCLUC

  cohort[pichrt].vconvrtflx.carbon  = ZERO;
  cohort[pichrt].vconvrtflx.nitrogen  = ZERO; 

// cohort[pichrt].vrespar initialized by LCLUC

  cohort[pichrt].yragstubC = ZERO;
  
  cohort[pichrt].yragstubN = ZERO;
  
  cohort[pichrt].yrcflux = ZERO;
  
  cohort[pichrt].yrCH4csmp = ZERO;
  
  cohort[pichrt].yrCH4ems = ZERO;
  
  cohort[pichrt].yrCH4flx = ZERO;
  
  cohort[pichrt].yrCO2dnflx = ZERO;
  
  cohort[pichrt].yrCO2nflx = ZERO;
  
  cohort[pichrt].yrconvrtC = ZERO;
  
  cohort[pichrt].yrconvrtN = ZERO;
  
  cohort[pichrt].yrdecayPROD1C = ZERO;
  
  cohort[pichrt].yrdecayPROD10C = ZERO;
  
  cohort[pichrt].yrdecayPROD100C = ZERO;

  cohort[pichrt].yrdecayPROD1N = ZERO;
  
  cohort[pichrt].yrdecayPROD10N = ZERO;
  
  cohort[pichrt].yrdecayPROD100N = ZERO;
  
  cohort[pichrt].yrdecayTOTPRODC = ZERO;

  cohort[pichrt].yrdecayTOTPRODN = ZERO;

  cohort[pichrt].yreet = ZERO;
  
  cohort[pichrt].yrfertn = ZERO;
  
  cohort[pichrt].yrfluxResidueC = ZERO;
  
  cohort[pichrt].yrfluxResidueN = ZERO;

  cohort[pichrt].yrformPROD1C = ZERO;
  
  cohort[pichrt].yrformPROD10C = ZERO;
  
  cohort[pichrt].yrformPROD100C = ZERO;

  cohort[pichrt].yrformPROD1N = ZERO;
  
  cohort[pichrt].yrformPROD10N = ZERO;
  
  cohort[pichrt].yrformPROD100N = ZERO;
  
  cohort[pichrt].yrformResidueC = ZERO;
  
  cohort[pichrt].yrformResidueN = ZERO;

  cohort[pichrt].yrformTOTPRODC = ZERO;

  cohort[pichrt].yrformTOTPRODN = ZERO;

  cohort[pichrt].yrfpc = ZERO;

  cohort[pichrt].yrgpp = ZERO;
  
  cohort[pichrt].yrgpr = ZERO;

//   cohort[pichrt].yrh2oyld = ZERO;
  
  cohort[pichrt].yrimmob = ZERO;
  
//  cohort[pichrt].yrineet = ZERO;
  
  cohort[pichrt].yringpp = ZERO;
  
  cohort[pichrt].yrinnpp = ZERO;
    
//  cohort[pichrt].yrirrig = ZERO;

  cohort[pichrt].yrlai = ZERO;

  cohort[pichrt].yrleaf = ZERO;

  cohort[pichrt].yrltrfalC = ZERO;

  cohort[pichrt].yrltrfalN = ZERO;  	
  
  cohort[pichrt].yrN2flx = ZERO;
  
  cohort[pichrt].yrN2Odnflx = ZERO;
  
  cohort[pichrt].yrN2Oflx = ZERO;
  
  cohort[pichrt].yrN2Onflx = ZERO;
  
  cohort[pichrt].yrnce = ZERO;
  
  cohort[pichrt].yrnep = ZERO;
  
  cohort[pichrt].yrninput = ZERO;
  
  cohort[pichrt].yrnlost = ZERO;
  
  cohort[pichrt].yrnmin = ZERO;
    
  cohort[pichrt].yrnpp = ZERO;
  
  cohort[pichrt].yrnrent = ZERO;
    
  cohort[pichrt].yrnsrent = ZERO;

  cohort[pichrt].yrnvrent = ZERO;
  
  cohort[pichrt].yrpet = ZERO;
    
//  cohort[pichrt].yrrgrndH2O = ZERO; 

  cohort[pichrt].yrrh = ZERO;

//  cohort[pichrt].yrrain = ZERO;
  
//  cohort[pichrt].yrrperc = ZERO;
  
//  cohort[pichrt].yrrrun = ZERO;
    
  cohort[pichrt].yrsconvrtC = ZERO;

  cohort[pichrt].yrsconvrtN = ZERO;

//  cohort[pichrt].yrsgrndH2O = ZERO;
  
  cohort[pichrt].yrslashC = ZERO;

  cohort[pichrt].yrslashN = ZERO;
  
//  cohort[pichrt].yrsnowfall = ZERO;

//  cohort[pichrt].yrsnowinf = ZERO;

//  cohort[pichrt].yrsnowpack = ZERO;

//  cohort[pichrt].yrsoilavlH2O = ZERO;

  cohort[pichrt].yrsoilavlN = ZERO;
  
//  cohort[pichrt].yrsoilC2N = ZERO;

//  cohort[pichrt].yrsoilmoist = ZERO;

  cohort[pichrt].yrsoilorgC = ZERO;
  
  cohort[pichrt].yrsoilorgN = ZERO;

//  cohort[pichrt].yrsoilpctp = ZERO;

//  cohort[pichrt].yrsoilvsm = ZERO;

//  cohort[pichrt].yrsperc = ZERO;
  
//  cohort[pichrt].yrsrun = ZERO;
  
  cohort[pichrt].yrtotalC = ZERO;
  
  cohort[pichrt].yrunleaf = ZERO;

  cohort[pichrt].yrvconvrtC = ZERO;
  
  cohort[pichrt].yrvconvrtN = ZERO;

  cohort[pichrt].yrvegC = ZERO;

// cohort[pichrt].yrvegC2N = ZERO;

  cohort[pichrt].yrveginnup = ZERO;

  cohort[pichrt].yrveglup = ZERO;
    
  cohort[pichrt].yrvegN = ZERO;
  
  cohort[pichrt].yrvegnmobil = ZERO;

  cohort[pichrt].yrvegnrsorb = ZERO;
  
  cohort[pichrt].yrvegnup = ZERO;

  cohort[pichrt].yrvegrgrowth = ZERO;

  cohort[pichrt].yrvegrmaint = ZERO;

  cohort[pichrt].yrvegSTON = ZERO;
  
  cohort[pichrt].yrvegSTRN = ZERO;

  cohort[pichrt].yrvegsup = ZERO;
	
};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void MITelmnt44::landUseChange( const int& pdyr )
{
  double abandonedANH4IN[NLVL];

  double abandonedANO3IN[NLVL];

  double abandonedDPHUMIN[NLVL];

  double abandonedNONSOLC = ZERO;

  double abandonedOCIN[NLVL];

  double abandonedRCLIN[NLVL];

  double abandonedRCRIN[NLVL];

  double abandonedRCVLIN[NLVL]; 

  long adjArea;

  int changedArea;

  double convertedANH4IN[NLVL];

  double convertedANO3IN[NLVL];

  double convertedDPHUMIN[NLVL];

  double convertedNONSOLC = ZERO;

  double convertedOCIN[NLVL];

  double convertedRCLIN[NLVL];

  double convertedRCRIN[NLVL];

  double convertedRCVLIN[NLVL];

  int dcmnt;

  long diffarea;

  long diffBiofuelArea = 0;
  
  long diffCropArea = 0;

  long diffPastureArea = 0;

  int dlyr;

  int dm;


  // Monthly formation of 10-year wood products
  Biomass formPROD10;

  // Monthly formation of 100-year wood products
  Biomass formPROD100;
 
 
  int i10;
  
  int i100;
  
  int ichrt;

  long initarea;  
 

  // Nitrogen from soil organic matter retained after 
  //   conversion
  double nsretent = ZERO;

//  long negadjArea;

  long negdiffarea;
  
  // Nitrogen from vegetation biomass retained after 
  //   conversion
  double nvretent = ZERO;
  
  // Proportion of abandoned area added to cohort
  double propAbandonedArea;

  // Proportion of " extra" biofuels area abandoned
  double propBiofuelArea; 
  
  // Proportion of converted area added to cohort
  double propConvertedArea;

  // Proportion of " extra" crop area abandoned
  double propCropArea; 

  // Proportion of "abandoned" area converted to
  //   other managed land uses

  double propManagedAbandoned; 

  // Proportion of " extra" pasture area abandoned
  double propPastureArea; 

  // Residual soil available N after conversion
//  double residualAvailN;
  
  // Residual soil organic matter after conversion
//  Biomass residualSoilOrg;

  // Soil organic matter lost during conversion to agriculture
  Biomass sconvrtflx;

  // Stock of slash resulting from conversion to agriculture
  Biomass slash;

  // Sum of area gained in managed cohorts (15, 16, 33)
  long sumManagedGain = 0;

  // Sum of area lost in managed cohorts (15, 16, 33)
  long sumManagedLoss = 0;

  // Sum of area gained in natural cohorts
  long sumNaturalGain = 0;

  // Sum of area lost in natural cohorts
  long sumNaturalLoss = 0;

  // Sum of converted area not accounted by net changes in area of Managed Lands
  long xtraManagedLoss = 0;

  // Vegetation biomass lost during conversion to agriculture
  Biomass vconvrtflx;
    
  lcluc.abandonedVeg.setSTRUCTC( ZERO );
  lcluc.abandonedVeg.setSTRUCTN( ZERO );
  lcluc.abandonedVeg.setLABILEN( ZERO );
      
  lcluc.abandonedSoil.setORGC( ZERO );
  lcluc.abandonedSoil.setORGN( ZERO );  	  	
  lcluc.abandonedSoil.setAVLN( ZERO );

  
  lcluc.convertedVeg.setSTRUCTC( ZERO );
  lcluc.convertedVeg.setSTRUCTN( ZERO );
  lcluc.convertedVeg.setLABILEN( ZERO );
  
  lcluc.convertedSoil.setORGC( ZERO );
  lcluc.convertedSoil.setORGN( ZERO );  	  	
  lcluc.convertedSoil.setAVLN( ZERO );

  for( dlyr = 0; dlyr < NLVL; ++dlyr )
  {
    abandonedANH4IN[dlyr] = ZERO;
    abandonedANO3IN[dlyr] = ZERO;
    abandonedDPHUMIN[dlyr] = ZERO;
    abandonedOCIN[dlyr] = ZERO;
    abandonedRCLIN[dlyr] = ZERO;
    abandonedRCRIN[dlyr] = ZERO;
    abandonedRCVLIN[dlyr] = ZERO; 

    convertedANH4IN[dlyr] = ZERO;
    convertedANO3IN[dlyr] = ZERO;
    convertedDPHUMIN[dlyr] = ZERO;
    convertedOCIN[dlyr] = ZERO;
    convertedRCLIN[dlyr] = ZERO;
    convertedRCRIN[dlyr] = ZERO;
    convertedRCVLIN[dlyr] = ZERO;
  }
  
  formPROD10.carbon = ZERO;
  formPROD10.nitrogen = ZERO;
  formPROD100.carbon = ZERO;
  formPROD100.nitrogen = ZERO;

  vconvrtflx.carbon = ZERO;
  vconvrtflx.nitrogen = ZERO;

  sconvrtflx.carbon = ZERO;
  sconvrtflx.nitrogen = ZERO;

  slash.carbon = ZERO;
  slash.nitrogen = ZERO;

  // Food Crops

  cohort[CROPVEG].formPROD10.carbon = ZERO;
  cohort[CROPVEG].formPROD10.nitrogen = ZERO;

  cohort[CROPVEG].formPROD100.carbon = ZERO;
  cohort[CROPVEG].formPROD100.nitrogen = ZERO;
  
  cohort[CROPVEG].slash.carbon = ZERO;
  cohort[CROPVEG].slash.nitrogen = ZERO;

  cohort[CROPVEG].vconvrtflx.carbon = ZERO;
  cohort[CROPVEG].sconvrtflx.carbon = ZERO;
  cohort[CROPVEG].convrtflx.carbon = ZERO;

  cohort[CROPVEG].vconvrtflx.nitrogen = ZERO;
  cohort[CROPVEG].sconvrtflx.nitrogen = ZERO;
  cohort[CROPVEG].convrtflx.nitrogen = ZERO;

  cohort[CROPVEG].nvretent = ZERO;
  cohort[CROPVEG].nsretent = ZERO;
  cohort[CROPVEG].nretent = ZERO;

  // Biofuels

  cohort[BIOFUELS].formPROD10.carbon = ZERO;
  cohort[BIOFUELS].formPROD10.nitrogen = ZERO;

  cohort[BIOFUELS].formPROD100.carbon = ZERO;
  cohort[BIOFUELS].formPROD100.nitrogen = ZERO;
  
  cohort[BIOFUELS].slash.carbon = ZERO;
  cohort[BIOFUELS].slash.nitrogen = ZERO;

  cohort[BIOFUELS].vconvrtflx.carbon = ZERO;
  cohort[BIOFUELS].sconvrtflx.carbon = ZERO;
  cohort[BIOFUELS].convrtflx.carbon = ZERO;

  cohort[BIOFUELS].vconvrtflx.nitrogen = ZERO;
  cohort[BIOFUELS].sconvrtflx.nitrogen = ZERO;
  cohort[BIOFUELS].convrtflx.nitrogen = ZERO;

  cohort[BIOFUELS].nvretent = ZERO;
  cohort[BIOFUELS].nsretent = ZERO;
  cohort[BIOFUELS].nretent = ZERO;

  // Pastures

  cohort[PASTURE].formPROD10.carbon = ZERO;
  cohort[PASTURE].formPROD10.nitrogen = ZERO;

  cohort[PASTURE].formPROD100.carbon = ZERO;
  cohort[PASTURE].formPROD100.nitrogen = ZERO;
  
  cohort[PASTURE].slash.carbon = ZERO;
  cohort[PASTURE].slash.nitrogen = ZERO;

  cohort[PASTURE].vconvrtflx.carbon = ZERO;
  cohort[PASTURE].sconvrtflx.carbon = ZERO;
  cohort[PASTURE].convrtflx.carbon = ZERO;

  cohort[PASTURE].vconvrtflx.nitrogen = ZERO;
  cohort[PASTURE].sconvrtflx.nitrogen = ZERO;
  cohort[PASTURE].convrtflx.nitrogen = ZERO;

  cohort[PASTURE].nvretent = ZERO;
  cohort[PASTURE].nsretent = ZERO;
  cohort[PASTURE].nretent = ZERO;


  // Check for changes in area of cohorts

  changedArea = 0;
  ichrt = 0;

//  cout << " pdyr = " << pdyr;
//  cout << " lat = " << row;
//  cout << endl;

  while( 0 == changedArea && ichrt < maxcohorts )
  {
    diffarea = cohort[ichrt].chrtarea - cohort[ichrt].prevchrtarea;

//    cout << " ichrt = " << ichrt;
//    cout << " chrtarea = " << cohort[ichrt].chrtarea;
//    cout << " prevchrtarea = " << cohort[ichrt].prevchrtarea;
//    cout << " diffarea = " << diffarea;
//    cout << endl;

    if( diffarea != 0 ) { changedArea = 1; }

    ++ichrt;
  }
   
//  cout << " changedArea = " << changedArea;
//  cout << endl;

  // Aggregate carbon and nitrogen in "changed" area into
  //   common pools to represent natural areas that have been
  //   converted to managed lands, or managed lands that have
  //   been abandoned if area changes among cohorts

  if( 1 == changedArea )
  {
    for( ichrt = 0; ichrt < maxcohorts; ++ichrt )
    {
      diffarea = 0;

      if( BAREGRND != cohort[ichrt].currentveg
          && GLACIERS != cohort[ichrt].currentveg
          && LAKES != cohort[ichrt].currentveg )
      {
        diffarea = cohort[ichrt].chrtarea - cohort[ichrt].prevchrtarea;

        switch( cohort[ichrt].currentveg )
        {
          case CROPVEG:  diffCropArea = diffarea; break;
          case BIOFUELS: diffBiofuelArea = diffarea; break;
          case PASTURE:  diffPastureArea = diffarea; break;
        }

        if( diffarea >= 0 )
        {
          if( CROPVEG == cohort[ichrt].currentveg
              || BIOFUELS == cohort[ichrt].currentveg
              || PASTURE == cohort[ichrt].currentveg )
          {
            sumManagedGain += diffarea;

          }
          else
          {
            sumNaturalGain += diffarea;
          }
        }
        else
        {
          if( CROPVEG == cohort[ichrt].currentveg
              || BIOFUELS == cohort[ichrt].currentveg
              || PASTURE == cohort[ichrt].currentveg )
          {
            sumManagedLoss += (diffarea * -1);

          }
          else
          {
            sumNaturalLoss += (diffarea * -1);
          }
        }
      }
	
      if( diffarea < 0 )
      {     
        negdiffarea = diffarea * -1;

        if( CROPVEG != cohort[ichrt].currentveg
  	    && BIOFUELS != cohort[ichrt].currentveg
  	    && PASTURE != cohort[ichrt].currentveg )
        {
          // Aggregate carbon and nitrogen from converted natural areas
          //   into common pools
        
  	  lcluc.convertedVeg.setSTRUCTC( lcluc.convertedVeg.getSTRUCTC() 
                                         + (cohort[ichrt].y[tem.I_VEGC]
  	                                 * negdiffarea) );
  	  	                                     
  	  lcluc.convertedVeg.setSTRUCTN( lcluc.convertedVeg.getSTRUCTN()
                                         + (cohort[ichrt].y[tem.I_STRN]
  	                                 * negdiffarea) );
  	  	                                     
  	  lcluc.convertedVeg.setLABILEN( lcluc.convertedVeg.getLABILEN()
                                         + (cohort[ichrt].y[tem.I_STON]
  	                                 * negdiffarea) );
  	  	
  	  lcluc.convertedSoil.setORGC( lcluc.convertedSoil.getORGC()
                                       + (cohort[ichrt].y[tem.I_SOLC]
  	                               * negdiffarea) );
  	  	                                     
  	  lcluc.convertedSoil.setORGN( lcluc.convertedSoil.getORGN()
                                       + (cohort[ichrt].y[tem.I_SOLN]
  	                               * negdiffarea) );
  	  	                                       	  	
          lcluc.convertedSoil.setAVLN( lcluc.convertedSoil.getAVLN()
                                       + (cohort[ichrt].y[tem.I_AVLN]
  	    	                     * negdiffarea) );

          convertedNONSOLC += (cohort[ichrt].NEMnsolc
                              * negdiffarea);

          for( dlyr = 0; dlyr < NLVL; ++dlyr )
          {
            convertedANH4IN[dlyr] += (cohort[ichrt].NEManh4in[dlyr]
                                     * negdiffarea);

            convertedANO3IN[dlyr] += (cohort[ichrt].NEMano3in[dlyr]
                                     * negdiffarea);

            convertedDPHUMIN[dlyr] += (cohort[ichrt].NEMdphumin[dlyr]
                                      * negdiffarea);

            convertedOCIN[dlyr] += (cohort[ichrt].NEMocin[dlyr]
                                   * negdiffarea);

            convertedRCLIN[dlyr] += (cohort[ichrt].NEMrclin[dlyr]
                                    * negdiffarea);

            convertedRCRIN[dlyr] += (cohort[ichrt].NEMrcrin[dlyr]
                                    * negdiffarea);

            convertedRCVLIN[dlyr] += (cohort[ichrt].NEMrcvlin[dlyr]
                                     * negdiffarea); 
          }


          dcmnt = cohort[ichrt].cmnt;

          // Create PROD10 and PROD100 from conversion to agriculture

          formPROD10.carbon += (cohort[ichrt].y[tem.I_VEGC]
  	                        * negdiffarea 
  	                        * tem.ag.getPROD10PAR( dcmnt ));

          formPROD10.nitrogen += ((cohort[ichrt].y[tem.I_STRN]
                                   + cohort[ichrt].y[tem.I_STON])
  	                           * negdiffarea 
  	                           * tem.ag.getPROD10PAR( dcmnt ));

          formPROD100.carbon += (cohort[ichrt].y[tem.I_VEGC]
  	                         * negdiffarea 
  	                         * tem.ag.getPROD100PAR( dcmnt ));

          formPROD100.nitrogen += ((cohort[ichrt].y[tem.I_STRN]
                                   + cohort[ichrt].y[tem.I_STON])
  	                           * negdiffarea 
  	                           * tem.ag.getPROD100PAR( dcmnt ));

          // Determine conversion fluxes

          slash.carbon += (cohort[ichrt].y[tem.I_VEGC]
  	                  * negdiffarea
  	                  * tem.ag.getSLASHPAR( dcmnt ));

          slash.nitrogen += ((cohort[ichrt].y[tem.I_STRN]
                            + cohort[ichrt].y[tem.I_STON])
  	                    * negdiffarea 
                            * tem.ag.getSLASHPAR( dcmnt ));

          vconvrtflx.carbon += (cohort[ichrt].y[tem.I_VEGC]
  	                       * negdiffarea
  	                       * tem.ag.getVCONVERT( dcmnt ));

          sconvrtflx.carbon += (cohort[ichrt].y[tem.I_SOLC]
  	                       * negdiffarea
  	                       * tem.ag.getSCONVERT( dcmnt ));


          vconvrtflx.nitrogen += ((cohort[ichrt].y[tem.I_STRN]
                                 + cohort[ichrt].y[tem.I_STON])
  	                         * negdiffarea 
  	                         * (1.0 - tem.ag.getNVRETCONV( dcmnt )) 
                                 * tem.ag.getVCONVERT( dcmnt ));
                                
          sconvrtflx.nitrogen += (cohort[ichrt].y[tem.I_SOLN]
  	                         * negdiffarea
  	                         * (1.0 - tem.ag.getNSRETCONV( dcmnt )) 
                                 * tem.ag.getSCONVERT( dcmnt ));
                                
                               
          nvretent += ((cohort[ichrt].y[tem.I_STRN]
                      + cohort[ichrt].y[tem.I_STON])
  	              * negdiffarea 
  	              * tem.ag.getNVRETCONV( dcmnt ) 
  	              * tem.ag.getVCONVERT( dcmnt )); 

          nsretent += (cohort[ichrt].y[tem.I_SOLN]
  	              * negdiffarea
  	              * tem.ag.getNSRETCONV( dcmnt ) 
  	              * tem.ag.getSCONVERT( dcmnt ));
        }
        else
        {
          // Aggregate carbon and nitrogen from abandoned managed areas
          //   into common pools

  	  lcluc.abandonedVeg.setSTRUCTC( lcluc.abandonedVeg.getSTRUCTC()
                                         + (cohort[ichrt].y[tem.I_VEGC]
  	                                 * negdiffarea) );
  	  	                                     
  	  lcluc.abandonedVeg.setSTRUCTN( lcluc.abandonedVeg.getSTRUCTN()
                                         + (cohort[ichrt].y[tem.I_STRN]
  	                                 * negdiffarea) );
  	  	                                     
  	  lcluc.abandonedVeg.setLABILEN( lcluc.abandonedVeg.getLABILEN()
                                         + (cohort[ichrt].y[tem.I_STON]
  	                                 * negdiffarea) );
  	  	
  	  lcluc.abandonedSoil.setORGC( lcluc.abandonedSoil.getORGC()
                                       + (cohort[ichrt].y[tem.I_SOLC]
  	                               * negdiffarea) );
  	  	                                     
  	  lcluc.abandonedSoil.setORGN( lcluc.abandonedSoil.getORGN()
                                       + (cohort[ichrt].y[tem.I_SOLN]
  	     	                       * negdiffarea) );
  	  	                                       	  	
          lcluc.abandonedSoil.setAVLN( lcluc.abandonedSoil.getAVLN()
                                       + (cohort[ichrt].y[tem.I_AVLN]
  	    	                       * negdiffarea) );  	  	

          abandonedNONSOLC += (cohort[ichrt].NEMnsolc
                               * negdiffarea);

          for( dlyr = 0; dlyr < NLVL; ++dlyr )
          {
            abandonedANH4IN[dlyr] += (cohort[ichrt].NEManh4in[dlyr]
                                     * negdiffarea);

            abandonedANO3IN[dlyr] += (cohort[ichrt].NEMano3in[dlyr]
                                     * negdiffarea);

            abandonedDPHUMIN[dlyr] += (cohort[ichrt].NEMdphumin[dlyr]
                                      * negdiffarea);

            abandonedOCIN[dlyr] += (cohort[ichrt].NEMocin[dlyr]
                                   * negdiffarea);

            abandonedRCLIN[dlyr] += (cohort[ichrt].NEMrclin[dlyr]
                                    * negdiffarea);

            abandonedRCRIN[dlyr] += (cohort[ichrt].NEMrcrin[dlyr]
                                    * negdiffarea);

            abandonedRCVLIN[dlyr] += (cohort[ichrt].NEMrcvlin[dlyr]
                                     * negdiffarea); 
          }
  	} 
      }
    } 	

    lcluc.convertedVeg.setPLANTN( lcluc.convertedVeg.getSTRUCTN()
                                  + lcluc.convertedVeg.getLABILEN() ); 

    lcluc.abandonedVeg.setPLANTN( lcluc.abandonedVeg.getSTRUCTN()
                                  + lcluc.abandonedVeg.getLABILEN() ); 

       

    if( sumManagedLoss < sumNaturalGain )
    {
      // Net changes in area do not adequately describe the amount of land
      //   abandoned or converted.  Adjust areas to account for the "extra"
      //   area not considered by net changes in area. 

      xtraManagedLoss = sumNaturalGain - sumManagedLoss;

      sumManagedLoss = sumNaturalGain;


      // Assume "extra" abandoned areas come from managed lands that have
      //   had a net gain in area

      if( diffCropArea > 0 && sumManagedGain > 0 )
      {
        propCropArea = (double) diffCropArea / (double) sumManagedGain;
      }
      else
      {
        propCropArea = ZERO;
      }

      if( diffBiofuelArea > 0 && sumManagedGain > 0 )
      {
        propBiofuelArea = (double) diffBiofuelArea / (double) sumManagedGain;
      }
      else
      {
        propBiofuelArea = ZERO;
      }

      if( diffPastureArea > 0 && sumManagedGain > 0 )
      {
        propPastureArea = (double) diffPastureArea / (double) sumManagedGain;
      }
      else
      {
        propPastureArea = ZERO;
      }


      // Assign carbon and nitrogen from "extra" abandoned managed areas
      //   into common pools and update carbon and nitrogen pools in
      //   affected managed areas

      lcluc.abandonedVeg.setSTRUCTC( lcluc.abandonedVeg.getSTRUCTC()
                                     + (cohort[CROPVEG].y[tem.I_VEGC]
  	                               * propCropArea
                                       * (double) xtraManagedLoss) 
                                     + (cohort[BIOFUELS].y[tem.I_VEGC]
  	                               * propBiofuelArea
                                       * (double) xtraManagedLoss)  
                                     + (cohort[PASTURE].y[tem.I_VEGC]
  	                               * propPastureArea
                                       * (double) xtraManagedLoss) ); 
  	  	                                     
      if( cohort[CROPVEG].prevchrtarea > 0 )
      {
        cohort[CROPVEG].y[tem.I_VEGC] -= ((cohort[CROPVEG].y[tem.I_VEGC]
  	                                 * propCropArea
                                         * (double) xtraManagedLoss)
                                         / (double) cohort[CROPVEG].prevchrtarea);
      }

      if( cohort[BIOFUELS].prevchrtarea > 0 )
      {
        cohort[BIOFUELS].y[tem.I_VEGC] -= ((cohort[BIOFUELS].y[tem.I_VEGC]
  	                                 * propBiofuelArea
                                         * (double) xtraManagedLoss)
                                         / (double) cohort[BIOFUELS].prevchrtarea);
      }

      if( cohort[PASTURE].prevchrtarea > 0 )
      {
        cohort[PASTURE].y[tem.I_VEGC] -= ((cohort[PASTURE].y[tem.I_VEGC]
  	                                 * propPastureArea
                                         * (double) xtraManagedLoss)
                                         / (double) cohort[PASTURE].prevchrtarea);
      }

      lcluc.abandonedVeg.setSTRUCTN( lcluc.abandonedVeg.getSTRUCTN()
                                     + (cohort[CROPVEG].y[tem.I_STRN]
  	                               * propCropArea
                                       * (double) xtraManagedLoss) 
                                     + (cohort[BIOFUELS].y[tem.I_STRN]
  	                               * propBiofuelArea
                                       * (double) xtraManagedLoss)  
                                     + (cohort[PASTURE].y[tem.I_STRN]
  	                               * propPastureArea
                                       * (double) xtraManagedLoss) );
  	  	                                     
      if( cohort[CROPVEG].prevchrtarea > 0 )
      {
        cohort[CROPVEG].y[tem.I_STRN] -= ((cohort[CROPVEG].y[tem.I_STRN]
  	                                 * propCropArea
                                         * (double) xtraManagedLoss)
                                         / (double) cohort[CROPVEG].prevchrtarea);
      }

      if( cohort[BIOFUELS].prevchrtarea > 0 )
      {
        cohort[BIOFUELS].y[tem.I_STRN] -= ((cohort[BIOFUELS].y[tem.I_STRN]
  	                                 * propBiofuelArea
                                         * (double) xtraManagedLoss)
                                         / (double) cohort[BIOFUELS].prevchrtarea);
      }

      if( cohort[PASTURE].prevchrtarea > 0 )
      {
        cohort[PASTURE].y[tem.I_STRN] -= ((cohort[PASTURE].y[tem.I_STRN]
  	                                 * propPastureArea
                                         * (double) xtraManagedLoss)
                                         / (double) cohort[PASTURE].prevchrtarea);
      }

      lcluc.abandonedVeg.setLABILEN( lcluc.abandonedVeg.getLABILEN()
                                     + (cohort[CROPVEG].y[tem.I_STON]
  	                               * propCropArea
                                       * (double) xtraManagedLoss)
                                     + (cohort[BIOFUELS].y[tem.I_STON]
  	                               * propBiofuelArea
                                       * (double) xtraManagedLoss)  
                                     + (cohort[PASTURE].y[tem.I_STON]
  	                               * propPastureArea
                                       * (double) xtraManagedLoss) );
  	  	
      if( cohort[CROPVEG].prevchrtarea > 0 )
      {
        cohort[CROPVEG].y[tem.I_STON] -= ((cohort[CROPVEG].y[tem.I_STON]
  	                                 * propCropArea
                                         * (double) xtraManagedLoss)
                                         / (double) cohort[CROPVEG].prevchrtarea);
      }

      if( cohort[BIOFUELS].prevchrtarea > 0 )
      {
        cohort[BIOFUELS].y[tem.I_STON] -= ((cohort[BIOFUELS].y[tem.I_STON]
  	                                 * propBiofuelArea
                                         * (double) xtraManagedLoss)
                                         / (double) cohort[BIOFUELS].prevchrtarea);
      }

      if( cohort[PASTURE].prevchrtarea > 0 )
      {
        cohort[PASTURE].y[tem.I_STON] -= ((cohort[PASTURE].y[tem.I_STON]
  	                                 * propPastureArea
                                         * (double) xtraManagedLoss)
                                         / (double) cohort[PASTURE].prevchrtarea);
      }

      lcluc.abandonedSoil.setORGC( lcluc.abandonedSoil.getORGC()
                                   + (cohort[CROPVEG].y[tem.I_SOLC]
  	                             * propCropArea
                                     * (double) xtraManagedLoss) 
                                   + (cohort[BIOFUELS].y[tem.I_SOLC]
  	                             * propBiofuelArea
                                     * (double) xtraManagedLoss)  
                                   + (cohort[PASTURE].y[tem.I_SOLC]
  	                             * propPastureArea
                                     * (double) xtraManagedLoss) );
  	  	                                     
      if( cohort[CROPVEG].prevchrtarea > 0 )
      {
        cohort[CROPVEG].y[tem.I_SOLC] -= ((cohort[CROPVEG].y[tem.I_SOLC]
  	                                 * propCropArea
                                         * (double) xtraManagedLoss)
                                         / (double) cohort[CROPVEG].prevchrtarea);
      }

      if( cohort[BIOFUELS].prevchrtarea > 0 )
      {
        cohort[BIOFUELS].y[tem.I_SOLC] -= ((cohort[BIOFUELS].y[tem.I_SOLC]
  	                                 * propBiofuelArea
                                         * (double) xtraManagedLoss)
                                         / (double) cohort[BIOFUELS].prevchrtarea);
      }

      if( cohort[PASTURE].prevchrtarea > 0 )
      {
        cohort[PASTURE].y[tem.I_SOLC] -= ((cohort[PASTURE].y[tem.I_SOLC]
  	                                 * propPastureArea
                                         * (double) xtraManagedLoss)
                                         / (double) cohort[PASTURE].prevchrtarea);
      }

      lcluc.abandonedSoil.setORGN( lcluc.abandonedSoil.getORGN()
                                   + (cohort[CROPVEG].y[tem.I_SOLN]
  	                             * propCropArea
                                     * (double) xtraManagedLoss) 
                                   + (cohort[BIOFUELS].y[tem.I_SOLN]
  	                             * propBiofuelArea
                                     * (double) xtraManagedLoss)  
                                   + (cohort[PASTURE].y[tem.I_SOLN]
  	                             * propPastureArea
                                     * (double) xtraManagedLoss) );
  	  	                                       	  	
      if( cohort[CROPVEG].prevchrtarea > 0 )
      {
        cohort[CROPVEG].y[tem.I_SOLN] -= ((cohort[CROPVEG].y[tem.I_SOLN]
  	                                 * propCropArea
                                         * (double) xtraManagedLoss)
                                         / (double) cohort[CROPVEG].prevchrtarea);
      }

      if( cohort[BIOFUELS].prevchrtarea > 0 )
      {
        cohort[BIOFUELS].y[tem.I_SOLN] -= ((cohort[BIOFUELS].y[tem.I_SOLN]
  	                                 * propBiofuelArea
                                         * (double) xtraManagedLoss)
                                         / (double) cohort[BIOFUELS].prevchrtarea);
      }

      if( cohort[PASTURE].prevchrtarea > 0 )
      {
        cohort[PASTURE].y[tem.I_SOLN] -= ((cohort[PASTURE].y[tem.I_SOLN]
  	                                 * propPastureArea
                                         * (double) xtraManagedLoss)
                                         / (double) cohort[PASTURE].prevchrtarea);
      }

      lcluc.abandonedSoil.setAVLN( lcluc.abandonedSoil.getAVLN()
                                   + (cohort[CROPVEG].y[tem.I_AVLN]
  	                             * propCropArea
                                     * (double) xtraManagedLoss)  
                                   + (cohort[BIOFUELS].y[tem.I_AVLN]
  	                             * propBiofuelArea
                                     * (double) xtraManagedLoss)  
                                   + (cohort[PASTURE].y[tem.I_AVLN]
  	                             * propPastureArea
                                     * (double) xtraManagedLoss) );  	  	

      if( cohort[CROPVEG].prevchrtarea > 0 )
      {
        cohort[CROPVEG].y[tem.I_AVLN] -= ((cohort[CROPVEG].y[tem.I_AVLN]
  	                                 * propCropArea
                                         * (double) xtraManagedLoss)
                                         / (double) cohort[CROPVEG].prevchrtarea);
      }

      if( cohort[BIOFUELS].prevchrtarea > 0 )
      {
        cohort[BIOFUELS].y[tem.I_AVLN] -= ((cohort[BIOFUELS].y[tem.I_AVLN]
  	                                 * propBiofuelArea
                                         * (double) xtraManagedLoss)
                                         / (double) cohort[BIOFUELS].prevchrtarea);
      }

      if( cohort[PASTURE].prevchrtarea > 0 )
      {
        cohort[PASTURE].y[tem.I_AVLN] -= ((cohort[PASTURE].y[tem.I_AVLN]
  	                                 * propPastureArea
                                         * (double) xtraManagedLoss)
                                         / (double) cohort[PASTURE].prevchrtarea);
      }

      abandonedNONSOLC += ((cohort[CROPVEG].NEMnsolc 
  	                     * propCropArea
                             * (double) xtraManagedLoss) 
                           + (cohort[BIOFUELS].NEMnsolc 
  	                     * propBiofuelArea
                             * (double) xtraManagedLoss)   
                           + (cohort[PASTURE].NEMnsolc 
  	                     * propPastureArea
                             * (double) xtraManagedLoss));

      if( cohort[CROPVEG].prevchrtarea > 0 )
      {
        cohort[CROPVEG].NEMnsolc -= ((cohort[CROPVEG].NEMnsolc
  	                           * propCropArea
                                   * (double) xtraManagedLoss)
                                   / (double) cohort[CROPVEG].prevchrtarea);
      }

      if( cohort[BIOFUELS].prevchrtarea > 0 )
      {
        cohort[BIOFUELS].NEMnsolc -= ((cohort[BIOFUELS].NEMnsolc
  	                            * propBiofuelArea
                                    * (double) xtraManagedLoss)
                                    / (double) cohort[BIOFUELS].prevchrtarea);
      }

      if( cohort[PASTURE].prevchrtarea > 0 )
      {
        cohort[PASTURE].NEMnsolc -= ((cohort[PASTURE].NEMnsolc
  	                           * propPastureArea
                                   * (double) xtraManagedLoss)
                                   / (double) cohort[PASTURE].prevchrtarea);
      }

      for( dlyr = 0; dlyr < NLVL; ++dlyr )
      {
        abandonedANH4IN[dlyr] += ((cohort[CROPVEG].NEManh4in[dlyr]
  	                           * propCropArea
                                   * (double) xtraManagedLoss) 
                                 + (cohort[BIOFUELS].NEManh4in[dlyr]
  	                           * propBiofuelArea
                                   * (double) xtraManagedLoss) 
                                 + (cohort[PASTURE].NEManh4in[dlyr]
  	                           * propPastureArea
                                   * (double) xtraManagedLoss));

        if( cohort[CROPVEG].prevchrtarea > 0 )
        {
          cohort[CROPVEG].NEManh4in[dlyr] -= ((cohort[CROPVEG].NEManh4in[dlyr]
  	                                     * propCropArea
                                             * (double) xtraManagedLoss)
                                             / (double) cohort[CROPVEG].prevchrtarea);
        }

        if( cohort[BIOFUELS].prevchrtarea > 0 )
        {
          cohort[BIOFUELS].NEManh4in[dlyr] -= ((cohort[BIOFUELS].NEManh4in[dlyr]
  	                                      * propBiofuelArea
                                              * (double) xtraManagedLoss)
                                              / (double) cohort[BIOFUELS].prevchrtarea);
        }

        if( cohort[PASTURE].prevchrtarea > 0 )
        {
          cohort[PASTURE].NEManh4in[dlyr] -= ((cohort[PASTURE].NEManh4in[dlyr]
  	                                   * propPastureArea
                                           * (double) xtraManagedLoss)
                                           / (double) cohort[PASTURE].prevchrtarea);
        }

        abandonedANO3IN[dlyr] += ((cohort[CROPVEG].NEMano3in[dlyr]
  	                           * propCropArea
                                   * (double) xtraManagedLoss)
                                 + (cohort[BIOFUELS].NEMano3in[dlyr]
  	                           * propBiofuelArea
                                   * (double) xtraManagedLoss) 
                                 + (cohort[PASTURE].NEMano3in[dlyr]
  	                           * propPastureArea
                                   * (double) xtraManagedLoss));

        if( cohort[CROPVEG].prevchrtarea > 0 )
        {
          cohort[CROPVEG].NEMano3in[dlyr] -= ((cohort[CROPVEG].NEMano3in[dlyr]
  	                                     * propCropArea
                                             * (double) xtraManagedLoss)
                                             / (double) cohort[CROPVEG].prevchrtarea);
        }

        if( cohort[BIOFUELS].prevchrtarea > 0 )
        {
          cohort[BIOFUELS].NEMano3in[dlyr] -= ((cohort[BIOFUELS].NEMano3in[dlyr]
  	                                      * propBiofuelArea
                                              * (double) xtraManagedLoss)
                                              / (double) cohort[BIOFUELS].prevchrtarea);
        }

        if( cohort[PASTURE].prevchrtarea > 0 )
        {
          cohort[PASTURE].NEMano3in[dlyr] -= ((cohort[PASTURE].NEMano3in[dlyr]
  	                                     * propPastureArea
                                             * (double) xtraManagedLoss)
                                             / (double) cohort[PASTURE].prevchrtarea);
        }

        abandonedDPHUMIN[dlyr] += ((cohort[CROPVEG].NEMdphumin[dlyr]
  	                            * propCropArea
                                    * (double) xtraManagedLoss) 
                                  + (cohort[BIOFUELS].NEMdphumin[dlyr]
  	                            * propBiofuelArea
                                    * (double) xtraManagedLoss) 
                                  + (cohort[PASTURE].NEMdphumin[dlyr]
  	                            * propPastureArea
                                    * (double) xtraManagedLoss));

        if( cohort[CROPVEG].prevchrtarea > 0 )
        {
          cohort[CROPVEG].NEMdphumin[dlyr] -= ((cohort[CROPVEG].NEMdphumin[dlyr]
  	                                      * propCropArea
                                              * (double) xtraManagedLoss)
                                              / (double) cohort[CROPVEG].prevchrtarea);
        }

        if( cohort[BIOFUELS].prevchrtarea > 0 )
        {
          cohort[BIOFUELS].NEMdphumin[dlyr] -= ((cohort[BIOFUELS].NEMdphumin[dlyr]
  	                                       * propBiofuelArea
                                               * (double) xtraManagedLoss)
                                               / (double) cohort[BIOFUELS].prevchrtarea);
        }

        if( cohort[PASTURE].prevchrtarea > 0 )
        {
          cohort[PASTURE].NEMdphumin[dlyr] -= ((cohort[PASTURE].NEMdphumin[dlyr]
  	                                      * propPastureArea
                                              * (double) xtraManagedLoss)
                                              / (double) cohort[PASTURE].prevchrtarea);
        }

        abandonedOCIN[dlyr] += ((cohort[CROPVEG].NEMocin[dlyr]
  	                         * propCropArea
                                 * (double) xtraManagedLoss) 
                               + (cohort[BIOFUELS].NEMocin[dlyr]
  	                         * propBiofuelArea
                                 * (double) xtraManagedLoss) 
                               + (cohort[PASTURE].NEMocin[dlyr]
  	                         * propPastureArea
                                 * (double) xtraManagedLoss));

        if( cohort[CROPVEG].prevchrtarea > 0 )
        {
          cohort[CROPVEG].NEMocin[dlyr] -= ((cohort[CROPVEG].NEMocin[dlyr]
  	                                   * propCropArea
                                           * (double) xtraManagedLoss)
                                           / (double) cohort[CROPVEG].prevchrtarea);
        }

        if( cohort[BIOFUELS].prevchrtarea > 0 )
        {
          cohort[BIOFUELS].NEMocin[dlyr] -= ((cohort[BIOFUELS].NEMocin[dlyr]
  	                                    * propBiofuelArea
                                            * (double) xtraManagedLoss)
                                            / (double) cohort[BIOFUELS].prevchrtarea);
        }

        if( cohort[PASTURE].prevchrtarea > 0 )
        {
          cohort[PASTURE].NEMocin[dlyr] -= ((cohort[PASTURE].NEMocin[dlyr]
  	                                   * propPastureArea
                                           * (double) xtraManagedLoss)
                                           / (double) cohort[PASTURE].prevchrtarea);
        }

        abandonedRCLIN[dlyr] += ((cohort[CROPVEG].NEMrclin[dlyr]
  	                          * propCropArea
                                  * (double) xtraManagedLoss) 
                                + (cohort[BIOFUELS].NEMrclin[dlyr]
  	                          * propBiofuelArea
                                  * (double) xtraManagedLoss) 
                                + (cohort[PASTURE].NEMrclin[dlyr]
  	                          * propPastureArea
                                  * (double) xtraManagedLoss));

        if( cohort[CROPVEG].prevchrtarea > 0 )
        {
          cohort[CROPVEG].NEMrclin[dlyr] -= ((cohort[CROPVEG].NEMrclin[dlyr]
  	                                    * propCropArea
                                            * (double) xtraManagedLoss)
                                            / (double) cohort[CROPVEG].prevchrtarea);
        }

        if( cohort[BIOFUELS].prevchrtarea > 0 )
        {
          cohort[BIOFUELS].NEMrclin[dlyr] -= ((cohort[BIOFUELS].NEMrclin[dlyr]
  	                                     * propBiofuelArea
                                             * (double) xtraManagedLoss)
                                             / (double) cohort[BIOFUELS].prevchrtarea);
        }

        if( cohort[PASTURE].prevchrtarea > 0 )
        {
          cohort[PASTURE].NEMrclin[dlyr] -= ((cohort[PASTURE].NEMrclin[dlyr]
  	                                    * propPastureArea
                                            * (double) xtraManagedLoss)
                                            / (double) cohort[PASTURE].prevchrtarea);
        }

        abandonedRCRIN[dlyr] += ((cohort[CROPVEG].NEMrcrin[dlyr]
  	                          * propCropArea
                                  * (double) xtraManagedLoss) 
                                + (cohort[BIOFUELS].NEMrcrin[dlyr]
  	                          * propBiofuelArea
                                  * (double) xtraManagedLoss) 
                                + (cohort[PASTURE].NEMrcrin[dlyr]
  	                          * propPastureArea
                                  * (double) xtraManagedLoss));

        if( cohort[CROPVEG].prevchrtarea > 0 )
        {
          cohort[CROPVEG].NEMrcrin[dlyr] -= ((cohort[CROPVEG].NEMrcrin[dlyr]
  	                                    * propCropArea
                                            * (double) xtraManagedLoss)
                                            / (double) cohort[CROPVEG].prevchrtarea);
        }

        if( cohort[BIOFUELS].prevchrtarea > 0 )
        {
          cohort[BIOFUELS].NEMrcrin[dlyr] -= ((cohort[BIOFUELS].NEMrcrin[dlyr]
  	                                     * propBiofuelArea
                                             * (double) xtraManagedLoss)
                                             / (double) cohort[BIOFUELS].prevchrtarea);
        }

        if( cohort[PASTURE].prevchrtarea > 0 )
        {
          cohort[PASTURE].NEMrcrin[dlyr] -= ((cohort[PASTURE].NEMrcrin[dlyr]
  	                                    * propPastureArea
                                            * (double) xtraManagedLoss)
                                            / (double) cohort[PASTURE].prevchrtarea);
        }

        abandonedRCVLIN[dlyr] += ((cohort[CROPVEG].NEMrcvlin[dlyr]
  	                           * propCropArea
                                   * (double) xtraManagedLoss) 
                                 + (cohort[BIOFUELS].NEMrcvlin[dlyr]
  	                           * propBiofuelArea
                                   * (double) xtraManagedLoss) 
                                 + (cohort[PASTURE].NEMrcvlin[dlyr]
  	                           * propPastureArea
                                   * (double) xtraManagedLoss)); 

        if( cohort[CROPVEG].prevchrtarea > 0 )
        {
          cohort[CROPVEG].NEMrcvlin[dlyr] -= ((cohort[CROPVEG].NEMrcvlin[dlyr]
  	                                     * propCropArea
                                             * (double) xtraManagedLoss)
                                             / (double) cohort[CROPVEG].prevchrtarea);
        }

        if( cohort[BIOFUELS].prevchrtarea > 0 )
        {
          cohort[BIOFUELS].NEMrcvlin[dlyr] -= ((cohort[BIOFUELS].NEMrcvlin[dlyr]
  	                                      * propBiofuelArea
                                              * (double) xtraManagedLoss)
                                              / (double) cohort[BIOFUELS].prevchrtarea);
        }

        if( cohort[PASTURE].prevchrtarea > 0 )
        {
          cohort[PASTURE].NEMrcvlin[dlyr] -= ((cohort[PASTURE].NEMrcvlin[dlyr]
  	                                     * propPastureArea
                                             * (double) xtraManagedLoss)
                                             / (double) cohort[PASTURE].prevchrtarea);
        }
      }
    }  

    // Reassign product pools to other managed areas when remaining area of a managed
    //   cohort is abandoned

    if( 0 == cohort[CROPVEG].chrtarea && 0 != cohort[CROPVEG].prevchrtarea )
    { 
      if( 0 != cohort[BIOFUELS].chrtarea && 0 != cohort[BIOFUELS].prevchrtarea )
      {
        // Agricultural products

        for( dm = 0; dm < CYCLE; ++dm )
        {
          cohort[BIOFUELS].initPROD1[dm].carbon = (((double) cohort[CROPVEG].prevchrtarea
                                                  *  cohort[CROPVEG].initPROD1[dm].carbon)
                                                  +  ((double) cohort[BIOFUELS].prevchrtarea
                                                  *  cohort[BIOFUELS].initPROD1[dm].carbon))
                                                     / (double) cohort[BIOFUELS].prevchrtarea;     

          cohort[BIOFUELS].initPROD1[dm].nitrogen = (((double) cohort[CROPVEG].prevchrtarea 
                                                    *  cohort[CROPVEG].initPROD1[dm].nitrogen)
                                                    + ((double) cohort[BIOFUELS].prevchrtarea 
                                                    *  cohort[BIOFUELS].initPROD1[dm].nitrogen))
                                                    / (double) cohort[BIOFUELS].prevchrtarea;   

          cohort[CROPVEG].initPROD1[dm].carbon = ZERO;
          cohort[CROPVEG].initPROD1[dm].nitrogen = ZERO;
        }

        cohort[BIOFUELS].prevPROD1.carbon = (((double) cohort[CROPVEG].prevchrtarea
                                            * cohort[CROPVEG].prevPROD1.carbon)
                                            + ((double) cohort[BIOFUELS].prevchrtarea
                                            * cohort[BIOFUELS].prevPROD1.carbon)) 
                                            / (double) cohort[BIOFUELS].prevchrtarea;     

        cohort[BIOFUELS].prevPROD1.nitrogen = (((double) cohort[CROPVEG].prevchrtarea
                                              * cohort[CROPVEG].prevPROD1.nitrogen)
                                              + ((double) cohort[BIOFUELS].prevchrtarea
                                              * cohort[BIOFUELS].prevPROD1.nitrogen))
                                              / (double) cohort[BIOFUELS].prevchrtarea;     

        cohort[CROPVEG].prevPROD1.carbon = ZERO;
        cohort[CROPVEG].prevPROD1.nitrogen = ZERO;

        // 10-year Woody Products

        for( i10 = 0; i10 < 10; ++i10 )
        {
          cohort[BIOFUELS].initPROD10[i10].carbon = (((double) cohort[CROPVEG].prevchrtarea
                                                    * cohort[CROPVEG].initPROD10[i10].carbon)
                                                    + ((double) cohort[BIOFUELS].prevchrtarea
                                                    * cohort[BIOFUELS].initPROD10[i10].carbon))
                                                     / (double) cohort[BIOFUELS].prevchrtarea;     

          cohort[BIOFUELS].initPROD10[i10].nitrogen = (((double) cohort[CROPVEG].prevchrtarea 
                                                      * cohort[CROPVEG].initPROD10[i10].nitrogen)
                                                      + ((double) cohort[BIOFUELS].prevchrtarea 
                                                      * cohort[BIOFUELS].initPROD10[i10].nitrogen))
                                                       / (double) cohort[BIOFUELS].prevchrtarea;   

          cohort[CROPVEG].initPROD10[i10].carbon = ZERO;
          cohort[CROPVEG].initPROD10[i10].nitrogen = ZERO;
        }

        cohort[BIOFUELS].prevPROD10.carbon = (((double) cohort[CROPVEG].prevchrtarea 
                                             * cohort[CROPVEG].prevPROD10.carbon)
                                             + ((double) cohort[BIOFUELS].prevchrtarea 
                                             * cohort[BIOFUELS].prevPROD10.carbon))
                                             / (double) cohort[BIOFUELS].prevchrtarea; 

        cohort[BIOFUELS].prevPROD10.nitrogen = (((double) cohort[CROPVEG].prevchrtarea 
                                               * cohort[CROPVEG].prevPROD10.nitrogen)
                                               + ((double) cohort[BIOFUELS].prevchrtarea 
                                               * cohort[BIOFUELS].prevPROD10.nitrogen))
                                               / (double) cohort[BIOFUELS].prevchrtarea; 

        cohort[CROPVEG].prevPROD10.carbon = ZERO;
        cohort[CROPVEG].prevPROD10.nitrogen = ZERO;

        // 100-year Woody Products

        for( i100 = 0; i100 < 100; ++i100 )
        {
          cohort[BIOFUELS].initPROD100[i100].carbon = (((double) cohort[CROPVEG].prevchrtarea
                                                      * cohort[CROPVEG].initPROD100[i100].carbon)
                                                      + ((double) cohort[BIOFUELS].prevchrtarea
                                                      * cohort[BIOFUELS].initPROD100[i100].carbon)) 
                                                      / (double) cohort[BIOFUELS].prevchrtarea;     

          cohort[BIOFUELS].initPROD100[i100].nitrogen = (((double) cohort[CROPVEG].prevchrtarea 
                                                        * cohort[CROPVEG].initPROD100[i100].nitrogen)
                                                        + ((double) cohort[BIOFUELS].prevchrtarea 
                                                        * cohort[BIOFUELS].initPROD100[i100].nitrogen)) 
                                                        /(double) cohort[BIOFUELS].prevchrtarea;   

          cohort[CROPVEG].initPROD100[i100].carbon = ZERO;
          cohort[CROPVEG].initPROD100[i100].nitrogen = ZERO;
        }

        cohort[BIOFUELS].prevPROD100.carbon = (((double) cohort[CROPVEG].prevchrtarea 
                                              * cohort[CROPVEG].prevPROD100.carbon)
                                              + ((double) cohort[BIOFUELS].prevchrtarea 
                                              * cohort[BIOFUELS].prevPROD100.carbon))
                                              / (double) cohort[BIOFUELS].prevchrtarea;     

        cohort[BIOFUELS].prevPROD100.nitrogen = (((double) cohort[CROPVEG].prevchrtarea
                                                * cohort[CROPVEG].prevPROD100.nitrogen)
                                                + ((double) cohort[BIOFUELS].prevchrtarea
                                                * cohort[BIOFUELS].prevPROD100.nitrogen))    
                                                / (double) cohort[BIOFUELS].prevchrtarea;     

        cohort[CROPVEG].prevPROD100.carbon = ZERO;
        cohort[CROPVEG].prevPROD100.nitrogen = ZERO;
      }
      else 
      {
        if( 0 != cohort[PASTURE].chrtarea && 0 != cohort[PASTURE].prevchrtarea )
        {
          // Agricultural products

          for( dm = 0; dm < CYCLE; ++dm )
          {
            cohort[PASTURE].initPROD1[dm].carbon = (((double) cohort[CROPVEG].prevchrtarea
                                                   *  cohort[CROPVEG].initPROD1[dm].carbon)
                                                   +  ((double) cohort[PASTURE].prevchrtarea
                                                   *  cohort[PASTURE].initPROD1[dm].carbon))
                                                   / (double) cohort[PASTURE].prevchrtarea;     

            cohort[PASTURE].initPROD1[dm].nitrogen = (((double) cohort[CROPVEG].prevchrtarea 
                                                     *  cohort[CROPVEG].initPROD1[dm].nitrogen)
                                                     + ((double) cohort[PASTURE].prevchrtarea 
                                                     *  cohort[PASTURE].initPROD1[dm].nitrogen))
                                                     / (double) cohort[PASTURE].prevchrtarea;   

            cohort[CROPVEG].initPROD1[dm].carbon = ZERO;
            cohort[CROPVEG].initPROD1[dm].nitrogen = ZERO;
          }

          cohort[PASTURE].prevPROD1.carbon = (((double) cohort[CROPVEG].prevchrtarea
                                             * cohort[CROPVEG].prevPROD1.carbon)
                                             + ((double) cohort[PASTURE].prevchrtarea
                                             * cohort[PASTURE].prevPROD1.carbon)) 
                                             / (double) cohort[PASTURE].prevchrtarea;     

          cohort[PASTURE].prevPROD1.nitrogen = (((double) cohort[CROPVEG].prevchrtarea
                                               * cohort[CROPVEG].prevPROD1.nitrogen)
                                               + ((double) cohort[PASTURE].prevchrtarea
                                               * cohort[PASTURE].prevPROD1.nitrogen))
                                               / (double) cohort[PASTURE].prevchrtarea;     

          cohort[CROPVEG].prevPROD1.carbon = ZERO;
          cohort[CROPVEG].prevPROD1.nitrogen = ZERO;

          // 10-year Woody Products

          for( i10 = 0; i10 < 10; ++i10 )
          {
            cohort[PASTURE].initPROD10[i10].carbon = (((double) cohort[CROPVEG].prevchrtarea
                                                     * cohort[CROPVEG].initPROD10[i10].carbon)
                                                     + ((double) cohort[PASTURE].prevchrtarea
                                                     * cohort[PASTURE].initPROD10[i10].carbon))
                                                     / (double) cohort[PASTURE].prevchrtarea;     

            cohort[PASTURE].initPROD10[i10].nitrogen = (((double) cohort[CROPVEG].prevchrtarea 
                                                       * cohort[CROPVEG].initPROD10[i10].nitrogen)
                                                       + ((double) cohort[PASTURE].prevchrtarea 
                                                       * cohort[PASTURE].initPROD10[i10].nitrogen))
                                                       / (double) cohort[PASTURE].prevchrtarea;   

            cohort[CROPVEG].initPROD10[i10].carbon = ZERO;
            cohort[CROPVEG].initPROD10[i10].nitrogen = ZERO;
          }

          cohort[PASTURE].prevPROD10.carbon = (((double) cohort[CROPVEG].prevchrtarea 
                                              * cohort[CROPVEG].prevPROD10.carbon)
                                              + ((double) cohort[PASTURE].prevchrtarea 
                                              * cohort[PASTURE].prevPROD10.carbon))
                                              / (double) cohort[PASTURE].prevchrtarea; 

          cohort[PASTURE].prevPROD10.nitrogen = (((double) cohort[CROPVEG].prevchrtarea 
                                                * cohort[CROPVEG].prevPROD10.nitrogen)
                                                + ((double) cohort[PASTURE].prevchrtarea 
                                                * cohort[PASTURE].prevPROD10.nitrogen))
                                                / (double) cohort[PASTURE].prevchrtarea; 

          cohort[CROPVEG].prevPROD10.carbon = ZERO;
          cohort[CROPVEG].prevPROD10.nitrogen = ZERO;

          // 100-year Woody Products

          for( i100 = 0; i100 < 100; ++i100 )
          {
            cohort[PASTURE].initPROD100[i100].carbon = (((double) cohort[CROPVEG].prevchrtarea
                                                       * cohort[CROPVEG].initPROD100[i100].carbon)
                                                       + ((double) cohort[PASTURE].prevchrtarea
                                                       * cohort[PASTURE].initPROD100[i100].carbon)) 
                                                       / (double) cohort[PASTURE].prevchrtarea;     

            cohort[PASTURE].initPROD100[i100].nitrogen = (((double) cohort[CROPVEG].prevchrtarea 
                                                          * cohort[CROPVEG].initPROD100[i100].nitrogen)
                                                          + ((double) cohort[PASTURE].prevchrtarea 
                                                          * cohort[PASTURE].initPROD100[i100].nitrogen)) 
                                                          /(double) cohort[PASTURE].prevchrtarea;   

            cohort[CROPVEG].initPROD100[i100].carbon = ZERO;
            cohort[CROPVEG].initPROD100[i100].nitrogen = ZERO;
          }

          cohort[PASTURE].prevPROD100.carbon = (((double) cohort[CROPVEG].prevchrtarea 
                                               * cohort[CROPVEG].prevPROD100.carbon)
                                               + ((double) cohort[PASTURE].prevchrtarea 
                                               * cohort[PASTURE].prevPROD100.carbon))
                                               / (double) cohort[PASTURE].prevchrtarea;     

          cohort[PASTURE].prevPROD100.nitrogen = (((double) cohort[CROPVEG].prevchrtarea
                                                 * cohort[CROPVEG].prevPROD100.nitrogen)
                                                 + ((double) cohort[PASTURE].prevchrtarea
                                                 * cohort[PASTURE].prevPROD100.nitrogen))    
                                                 / (double) cohort[PASTURE].prevchrtarea;     

          cohort[CROPVEG].prevPROD100.carbon = ZERO;
          cohort[CROPVEG].prevPROD100.nitrogen = ZERO;
        } 
      }
    }

    if( 0 == cohort[BIOFUELS].chrtarea && 0 != cohort[BIOFUELS].prevchrtarea )
    { 
      if( 0 != cohort[PASTURE].chrtarea && 0 != cohort[PASTURE].prevchrtarea )
      {
        // Agricultural products

        for( dm = 0; dm < CYCLE; ++dm )
        {
          cohort[PASTURE].initPROD1[dm].carbon = (((double) cohort[BIOFUELS].prevchrtarea
                                                 *  cohort[BIOFUELS].initPROD1[dm].carbon)
                                                 +  ((double) cohort[PASTURE].prevchrtarea
                                                 *  cohort[PASTURE].initPROD1[dm].carbon))
                                                 / (double) cohort[PASTURE].prevchrtarea;     

          cohort[PASTURE].initPROD1[dm].nitrogen = (((double) cohort[BIOFUELS].prevchrtarea 
                                                   *  cohort[BIOFUELS].initPROD1[dm].nitrogen)
                                                   + ((double) cohort[PASTURE].prevchrtarea 
                                                   *  cohort[PASTURE].initPROD1[dm].nitrogen))
                                                   / (double) cohort[PASTURE].prevchrtarea;   

          cohort[BIOFUELS].initPROD1[dm].carbon = ZERO;
          cohort[BIOFUELS].initPROD1[dm].nitrogen = ZERO;
        }

        cohort[PASTURE].prevPROD1.carbon = (((double) cohort[BIOFUELS].prevchrtarea
                                           * cohort[BIOFUELS].prevPROD1.carbon)
                                           + ((double) cohort[PASTURE].prevchrtarea
                                           * cohort[PASTURE].prevPROD1.carbon)) 
                                           / (double) cohort[PASTURE].prevchrtarea;     

        cohort[PASTURE].prevPROD1.nitrogen = (((double) cohort[BIOFUELS].prevchrtarea
                                             * cohort[BIOFUELS].prevPROD1.nitrogen)
                                             + ((double) cohort[PASTURE].prevchrtarea
                                             * cohort[PASTURE].prevPROD1.nitrogen))
                                             / (double) cohort[PASTURE].prevchrtarea;     

        cohort[BIOFUELS].prevPROD1.carbon = ZERO;
        cohort[BIOFUELS].prevPROD1.nitrogen = ZERO;

        // 10-year Woody Products

        for( i10 = 0; i10 < 10; ++i10 )
        {
          cohort[PASTURE].initPROD10[i10].carbon = (((double) cohort[BIOFUELS].prevchrtarea
                                                   * cohort[BIOFUELS].initPROD10[i10].carbon)
                                                   + ((double) cohort[PASTURE].prevchrtarea
                                                   * cohort[PASTURE].initPROD10[i10].carbon))
                                                   / (double) cohort[PASTURE].prevchrtarea;     

          cohort[PASTURE].initPROD10[i10].nitrogen = (((double) cohort[BIOFUELS].prevchrtarea 
                                                     * cohort[BIOFUELS].initPROD10[i10].nitrogen)
                                                     + ((double) cohort[PASTURE].prevchrtarea 
                                                     * cohort[PASTURE].initPROD10[i10].nitrogen))
                                                     / (double) cohort[PASTURE].prevchrtarea;   

          cohort[BIOFUELS].initPROD10[i10].carbon = ZERO;
          cohort[BIOFUELS].initPROD10[i10].nitrogen = ZERO;
        }

        cohort[PASTURE].prevPROD10.carbon = (((double) cohort[BIOFUELS].prevchrtarea 
                                            * cohort[BIOFUELS].prevPROD10.carbon)
                                            + ((double) cohort[PASTURE].prevchrtarea 
                                            * cohort[PASTURE].prevPROD10.carbon))
                                            / (double) cohort[PASTURE].prevchrtarea; 

        cohort[PASTURE].prevPROD10.nitrogen = (((double) cohort[BIOFUELS].prevchrtarea 
                                              * cohort[BIOFUELS].prevPROD10.nitrogen)
                                              + ((double) cohort[PASTURE].prevchrtarea 
                                              * cohort[PASTURE].prevPROD10.nitrogen))
                                              / (double) cohort[PASTURE].prevchrtarea; 

        cohort[BIOFUELS].prevPROD10.carbon = ZERO;
        cohort[BIOFUELS].prevPROD10.nitrogen = ZERO;

        // 100-year Woody Products

        for( i100 = 0; i100 < 100; ++i100 )
        {
          cohort[PASTURE].initPROD100[i100].carbon = (((double) cohort[BIOFUELS].prevchrtarea
                                                     * cohort[BIOFUELS].initPROD100[i100].carbon)
                                                     + ((double) cohort[PASTURE].prevchrtarea
                                                     * cohort[PASTURE].initPROD100[i100].carbon)) 
                                                     / (double) cohort[PASTURE].prevchrtarea;     

          cohort[PASTURE].initPROD100[i100].nitrogen = (((double) cohort[BIOFUELS].prevchrtarea 
                                                        * cohort[BIOFUELS].initPROD100[i100].nitrogen)
                                                        + ((double) cohort[PASTURE].prevchrtarea 
                                                        * cohort[PASTURE].initPROD100[i100].nitrogen)) 
                                                        /(double) cohort[PASTURE].prevchrtarea;   

          cohort[BIOFUELS].initPROD100[i100].carbon = ZERO;
          cohort[BIOFUELS].initPROD100[i100].nitrogen = ZERO;
        }

        cohort[PASTURE].prevPROD100.carbon = (((double) cohort[BIOFUELS].prevchrtarea 
                                             * cohort[BIOFUELS].prevPROD100.carbon)
                                             + ((double) cohort[PASTURE].prevchrtarea 
                                             * cohort[PASTURE].prevPROD100.carbon))
                                             / (double) cohort[PASTURE].prevchrtarea;     

        cohort[PASTURE].prevPROD100.nitrogen = (((double) cohort[BIOFUELS].prevchrtarea
                                               * cohort[BIOFUELS].prevPROD100.nitrogen)
                                               + ((double) cohort[PASTURE].prevchrtarea
                                               * cohort[PASTURE].prevPROD100.nitrogen))    
                                               / (double) cohort[PASTURE].prevchrtarea;     

        cohort[BIOFUELS].prevPROD100.carbon = ZERO;
        cohort[BIOFUELS].prevPROD100.nitrogen = ZERO;
      } 
      else
      {
        if( 0 != cohort[CROPVEG].chrtarea && 0 != cohort[CROPVEG].prevchrtarea )
        {
          // Agricultural products

          for( dm = 0; dm < CYCLE; ++dm )
          {
            cohort[CROPVEG].initPROD1[dm].carbon = (((double) cohort[BIOFUELS].prevchrtarea
                                                   *  cohort[BIOFUELS].initPROD1[dm].carbon)
                                                   +  ((double) cohort[CROPVEG].prevchrtarea
                                                   *  cohort[CROPVEG].initPROD1[dm].carbon))
                                                   / (double) cohort[CROPVEG].prevchrtarea;     

            cohort[CROPVEG].initPROD1[dm].nitrogen = (((double) cohort[BIOFUELS].prevchrtarea 
                                                     *  cohort[BIOFUELS].initPROD1[dm].nitrogen)
                                                     + ((double) cohort[CROPVEG].prevchrtarea 
                                                     *  cohort[CROPVEG].initPROD1[dm].nitrogen))
                                                     / (double) cohort[CROPVEG].prevchrtarea;   

            cohort[BIOFUELS].initPROD1[dm].carbon = ZERO;
            cohort[BIOFUELS].initPROD1[dm].nitrogen = ZERO;
          }

          cohort[CROPVEG].prevPROD1.carbon = (((double) cohort[BIOFUELS].prevchrtarea
                                             * cohort[BIOFUELS].prevPROD1.carbon)
                                             + ((double) cohort[CROPVEG].prevchrtarea
                                             * cohort[CROPVEG].prevPROD1.carbon)) 
                                             / (double) cohort[CROPVEG].prevchrtarea;     

          cohort[CROPVEG].prevPROD1.nitrogen = (((double) cohort[BIOFUELS].prevchrtarea
                                               * cohort[BIOFUELS].prevPROD1.nitrogen)
                                               + ((double) cohort[CROPVEG].prevchrtarea
                                               * cohort[CROPVEG].prevPROD1.nitrogen))
                                               / (double) cohort[CROPVEG].prevchrtarea;     
  
          cohort[BIOFUELS].prevPROD1.carbon = ZERO;
          cohort[BIOFUELS].prevPROD1.nitrogen = ZERO;

        // 10-year Woody Products

        for( i10 = 0; i10 < 10; ++i10 )
        {
          cohort[CROPVEG].initPROD10[i10].carbon = (((double) cohort[BIOFUELS].prevchrtarea
                                                   * cohort[BIOFUELS].initPROD10[i10].carbon)
                                                   + ((double) cohort[CROPVEG].prevchrtarea
                                                   * cohort[CROPVEG].initPROD10[i10].carbon))
                                                   / (double) cohort[CROPVEG].prevchrtarea;     

          cohort[CROPVEG].initPROD10[i10].nitrogen = (((double) cohort[BIOFUELS].prevchrtarea 
                                                     * cohort[BIOFUELS].initPROD10[i10].nitrogen)
                                                     + ((double) cohort[CROPVEG].prevchrtarea 
                                                     * cohort[CROPVEG].initPROD10[i10].nitrogen))
                                                     / (double) cohort[CROPVEG].prevchrtarea;   

          cohort[BIOFUELS].initPROD10[i10].carbon = ZERO;
          cohort[BIOFUELS].initPROD10[i10].nitrogen = ZERO;
        }

        cohort[CROPVEG].prevPROD10.carbon = (((double) cohort[BIOFUELS].prevchrtarea 
                                            * cohort[BIOFUELS].prevPROD10.carbon)
                                            + ((double) cohort[CROPVEG].prevchrtarea 
                                            * cohort[CROPVEG].prevPROD10.carbon))
                                            / (double) cohort[CROPVEG].prevchrtarea; 

        cohort[CROPVEG].prevPROD10.nitrogen = (((double) cohort[BIOFUELS].prevchrtarea 
                                              * cohort[BIOFUELS].prevPROD10.nitrogen)
                                              + ((double) cohort[CROPVEG].prevchrtarea 
                                              * cohort[CROPVEG].prevPROD10.nitrogen))
                                              / (double) cohort[CROPVEG].prevchrtarea; 

        cohort[BIOFUELS].prevPROD10.carbon = ZERO;
        cohort[BIOFUELS].prevPROD10.nitrogen = ZERO;

        // 100-year Woody Products

        for( i100 = 0; i100 < 100; ++i100 )
        {
          cohort[CROPVEG].initPROD100[i100].carbon = (((double) cohort[BIOFUELS].prevchrtarea
                                                     * cohort[BIOFUELS].initPROD100[i100].carbon)
                                                     + ((double) cohort[CROPVEG].prevchrtarea
                                                     * cohort[CROPVEG].initPROD100[i100].carbon)) 
                                                     / (double) cohort[CROPVEG].prevchrtarea;     

          cohort[CROPVEG].initPROD100[i100].nitrogen = (((double) cohort[BIOFUELS].prevchrtarea 
                                                        * cohort[BIOFUELS].initPROD100[i100].nitrogen)
                                                        + ((double) cohort[CROPVEG].prevchrtarea 
                                                        * cohort[CROPVEG].initPROD100[i100].nitrogen)) 
                                                        /(double) cohort[CROPVEG].prevchrtarea;   

          cohort[BIOFUELS].initPROD100[i100].carbon = ZERO;
          cohort[BIOFUELS].initPROD100[i100].nitrogen = ZERO;
        }

        cohort[CROPVEG].prevPROD100.carbon = (((double) cohort[BIOFUELS].prevchrtarea 
                                             * cohort[BIOFUELS].prevPROD100.carbon)
                                             + ((double) cohort[CROPVEG].prevchrtarea 
                                             * cohort[CROPVEG].prevPROD100.carbon))
                                             / (double) cohort[CROPVEG].prevchrtarea;     

        cohort[CROPVEG].prevPROD100.nitrogen = (((double) cohort[BIOFUELS].prevchrtarea
                                               * cohort[BIOFUELS].prevPROD100.nitrogen)
                                               + ((double) cohort[CROPVEG].prevchrtarea
                                               * cohort[CROPVEG].prevPROD100.nitrogen))    
                                               / (double) cohort[CROPVEG].prevchrtarea;     

        cohort[BIOFUELS].prevPROD100.carbon = ZERO;
        cohort[BIOFUELS].prevPROD100.nitrogen = ZERO;

        }
      }
    }

    if( 0 == cohort[PASTURE].chrtarea && 0 != cohort[PASTURE].prevchrtarea )
    {
      if( 0 != cohort[CROPVEG].chrtarea && 0 != cohort[CROPVEG].prevchrtarea )
      {
        // 10-year Woody Products

        for( i10 = 0; i10 < 10; ++i10 )
        {
          cohort[CROPVEG].initPROD10[i10].carbon = (((double) cohort[PASTURE].prevchrtarea
                                                   * cohort[PASTURE].initPROD10[i10].carbon)
                                                   + ((double) cohort[CROPVEG].prevchrtarea
                                                   * cohort[CROPVEG].initPROD10[i10].carbon))
                                                   / (double) cohort[CROPVEG].prevchrtarea;     

          cohort[CROPVEG].initPROD10[i10].nitrogen = (((double) cohort[PASTURE].prevchrtarea 
                                                     * cohort[PASTURE].initPROD10[i10].nitrogen)
                                                     + ((double) cohort[CROPVEG].prevchrtarea 
                                                     * cohort[CROPVEG].initPROD10[i10].nitrogen))
                                                     / (double) cohort[CROPVEG].prevchrtarea;   

          cohort[PASTURE].initPROD10[i10].carbon = ZERO;
          cohort[PASTURE].initPROD10[i10].nitrogen = ZERO;
        }

        cohort[CROPVEG].prevPROD10.carbon = (((double) cohort[PASTURE].prevchrtarea 
                                            * cohort[PASTURE].prevPROD10.carbon)
                                            + ((double) cohort[CROPVEG].prevchrtarea 
                                            * cohort[CROPVEG].prevPROD10.carbon))
                                            / (double) cohort[CROPVEG].prevchrtarea; 

        cohort[CROPVEG].prevPROD10.nitrogen = (((double) cohort[PASTURE].prevchrtarea 
                                              * cohort[PASTURE].prevPROD10.nitrogen)
                                              + ((double) cohort[CROPVEG].prevchrtarea 
                                              * cohort[CROPVEG].prevPROD10.nitrogen))
                                              / (double) cohort[CROPVEG].prevchrtarea; 

        cohort[PASTURE].prevPROD10.carbon = ZERO;
        cohort[PASTURE].prevPROD10.nitrogen = ZERO;

        // 100-year Woody Products

        for( i100 = 0; i100 < 100; ++i100 )
        {
          cohort[CROPVEG].initPROD100[i100].carbon = (((double) cohort[PASTURE].prevchrtarea
                                                     * cohort[PASTURE].initPROD100[i100].carbon)
                                                     + ((double) cohort[CROPVEG].prevchrtarea
                                                     * cohort[CROPVEG].initPROD100[i100].carbon)) 
                                                     / (double) cohort[CROPVEG].prevchrtarea;     

          cohort[CROPVEG].initPROD100[i100].nitrogen = (((double) cohort[PASTURE].prevchrtarea 
                                                       * cohort[PASTURE].initPROD100[i100].nitrogen)
                                                       + ((double) cohort[CROPVEG].prevchrtarea 
                                                       * cohort[CROPVEG].initPROD100[i100].nitrogen)) 
                                                       /(double) cohort[CROPVEG].prevchrtarea;   

          cohort[PASTURE].initPROD100[i100].carbon = ZERO;
          cohort[PASTURE].initPROD100[i100].nitrogen = ZERO;
        }

        cohort[CROPVEG].prevPROD100.carbon = (((double) cohort[PASTURE].prevchrtarea 
                                             * cohort[PASTURE].prevPROD100.carbon)
                                             + ((double) cohort[CROPVEG].prevchrtarea 
                                             * cohort[CROPVEG].prevPROD100.carbon))
                                             / (double) cohort[CROPVEG].prevchrtarea;     

        cohort[CROPVEG].prevPROD100.nitrogen = (((double) cohort[PASTURE].prevchrtarea
                                               * cohort[PASTURE].prevPROD100.nitrogen)
                                               + ((double) cohort[CROPVEG].prevchrtarea
                                               * cohort[CROPVEG].prevPROD100.nitrogen))    
                                               / (double) cohort[CROPVEG].prevchrtarea;     

        cohort[PASTURE].prevPROD100.carbon = ZERO;
        cohort[PASTURE].prevPROD100.nitrogen = ZERO;
      }
      else 
      {
        if( 0 != cohort[BIOFUELS].chrtarea && 0 != cohort[BIOFUELS].prevchrtarea )
        {
          // 10-year Woody Products

          for( i10 = 0; i10 < 10; ++i10 )
          {
            cohort[BIOFUELS].initPROD10[i10].carbon = (((double) cohort[PASTURE].prevchrtarea
                                                      * cohort[PASTURE].initPROD10[i10].carbon)
                                                      + ((double) cohort[BIOFUELS].prevchrtarea
                                                      * cohort[BIOFUELS].initPROD10[i10].carbon))
                                                      / (double) cohort[BIOFUELS].prevchrtarea;     

            cohort[BIOFUELS].initPROD10[i10].nitrogen = (((double) cohort[PASTURE].prevchrtarea 
                                                        * cohort[PASTURE].initPROD10[i10].nitrogen)
                                                        + ((double) cohort[BIOFUELS].prevchrtarea 
                                                        * cohort[BIOFUELS].initPROD10[i10].nitrogen))
                                                        / (double) cohort[BIOFUELS].prevchrtarea;   

            cohort[PASTURE].initPROD10[i10].carbon = ZERO;
            cohort[PASTURE].initPROD10[i10].nitrogen = ZERO;
          }

          cohort[BIOFUELS].prevPROD10.carbon = (((double) cohort[PASTURE].prevchrtarea 
                                               * cohort[PASTURE].prevPROD10.carbon)
                                               + ((double) cohort[BIOFUELS].prevchrtarea 
                                               * cohort[BIOFUELS].prevPROD10.carbon))
                                               / (double) cohort[BIOFUELS].prevchrtarea; 

          cohort[BIOFUELS].prevPROD10.nitrogen = (((double) cohort[PASTURE].prevchrtarea 
                                                 * cohort[PASTURE].prevPROD10.nitrogen)
                                                 + ((double) cohort[BIOFUELS].prevchrtarea 
                                                 * cohort[BIOFUELS].prevPROD10.nitrogen))
                                                 / (double) cohort[BIOFUELS].prevchrtarea; 

          cohort[PASTURE].prevPROD10.carbon = ZERO;
          cohort[PASTURE].prevPROD10.nitrogen = ZERO;

          // 100-year Woody Products

          for( i100 = 0; i100 < 100; ++i100 )
          {
            cohort[BIOFUELS].initPROD100[i100].carbon = (((double) cohort[PASTURE].prevchrtarea
                                                        * cohort[PASTURE].initPROD100[i100].carbon)
                                                        + ((double) cohort[BIOFUELS].prevchrtarea
                                                        * cohort[BIOFUELS].initPROD100[i100].carbon)) 
                                                        / (double) cohort[BIOFUELS].prevchrtarea;     

            cohort[BIOFUELS].initPROD100[i100].nitrogen = (((double) cohort[PASTURE].prevchrtarea 
                                                          * cohort[PASTURE].initPROD100[i100].nitrogen)
                                                          + ((double) cohort[BIOFUELS].prevchrtarea 
                                                          * cohort[BIOFUELS].initPROD100[i100].nitrogen)) 
                                                          /(double) cohort[BIOFUELS].prevchrtarea;   

            cohort[PASTURE].initPROD100[i100].carbon = ZERO;
            cohort[PASTURE].initPROD100[i100].nitrogen = ZERO;
          }

          cohort[BIOFUELS].prevPROD100.carbon = (((double) cohort[PASTURE].prevchrtarea 
                                               * cohort[PASTURE].prevPROD100.carbon)
                                               + ((double) cohort[BIOFUELS].prevchrtarea 
                                               * cohort[BIOFUELS].prevPROD100.carbon))
                                               / (double) cohort[BIOFUELS].prevchrtarea;     

          cohort[BIOFUELS].prevPROD100.nitrogen = (((double) cohort[PASTURE].prevchrtarea
                                                 * cohort[PASTURE].prevPROD100.nitrogen)
                                                 + ((double) cohort[BIOFUELS].prevchrtarea
                                                 * cohort[BIOFUELS].prevPROD100.nitrogen))    
                                                 / (double) cohort[BIOFUELS].prevchrtarea;     

          cohort[PASTURE].prevPROD100.carbon = ZERO;
          cohort[PASTURE].prevPROD100.nitrogen = ZERO;
        } 
      }
    }

 
    // Readjust product pools to account for changes in area of crops and pastures
 
    if( 0 != cohort[CROPVEG].chrtarea && 0 != cohort[CROPVEG].prevchrtarea )
    { 
      // Agricultural products

      for( dm = 0; dm < CYCLE; ++dm )
      {
        cohort[CROPVEG].initPROD1[dm].carbon *= ((double) cohort[CROPVEG].prevchrtarea 
                                                   / (double) cohort[CROPVEG].chrtarea);     

        cohort[CROPVEG].initPROD1[dm].nitrogen *= ((double) cohort[CROPVEG].prevchrtarea 
                                                     / (double) cohort[CROPVEG].chrtarea);   
      }

      cohort[CROPVEG].prevPROD1.carbon *= ((double) cohort[CROPVEG].prevchrtarea 
                                                   / (double) cohort[CROPVEG].chrtarea);     

      cohort[CROPVEG].prevPROD1.nitrogen *= ((double) cohort[CROPVEG].prevchrtarea 
                                                   / (double) cohort[CROPVEG].chrtarea);     

      // 10-year Woody Products

      for( i10 = 0; i10 < 10; ++i10 )
      {
        cohort[CROPVEG].initPROD10[i10].carbon *= ((double) cohort[CROPVEG].prevchrtarea 
                                                   / (double) cohort[CROPVEG].chrtarea);     

        cohort[CROPVEG].initPROD10[i10].nitrogen *= ((double) cohort[CROPVEG].prevchrtarea 
                                                     / (double) cohort[CROPVEG].chrtarea);   
      }

      cohort[CROPVEG].prevPROD10.carbon *= ((double) cohort[CROPVEG].prevchrtarea 
                                           / (double) cohort[CROPVEG].chrtarea); 

      cohort[CROPVEG].prevPROD10.nitrogen *= ((double) cohort[CROPVEG].prevchrtarea 
                                             / (double) cohort[CROPVEG].chrtarea); 

      // 100-year Woody Products

      for( i100 = 0; i100 < 100; ++i100 )
      {
        cohort[CROPVEG].initPROD100[i100].carbon *= ((double) cohort[CROPVEG].prevchrtarea 
                                                     / (double) cohort[CROPVEG].chrtarea);     

        cohort[CROPVEG].initPROD100[i100].nitrogen *= ((double) cohort[CROPVEG].prevchrtarea 
                                                       /(double) cohort[CROPVEG].chrtarea);   
      }

      cohort[CROPVEG].prevPROD100.carbon *= ((double) cohort[CROPVEG].prevchrtarea 
                                            / (double) cohort[CROPVEG].chrtarea);     

      cohort[CROPVEG].prevPROD100.nitrogen *= ((double) cohort[CROPVEG].prevchrtarea    
                                              / (double) cohort[CROPVEG].chrtarea);     
    }

    if( 0 != cohort[BIOFUELS].chrtarea && 0 != cohort[BIOFUELS].prevchrtarea )
    { 
      // Agricultural products

      for( dm = 0; dm < CYCLE; ++dm )
      {
        cohort[BIOFUELS].initPROD1[dm].carbon *= ((double) cohort[BIOFUELS].prevchrtarea 
                                                   / (double) cohort[BIOFUELS].chrtarea);     

        cohort[BIOFUELS].initPROD1[dm].nitrogen *= ((double) cohort[BIOFUELS].prevchrtarea 
                                                     / (double) cohort[BIOFUELS].chrtarea);   
      }

      cohort[BIOFUELS].prevPROD1.carbon *= ((double) cohort[BIOFUELS].prevchrtarea 
                                                   / (double) cohort[BIOFUELS].chrtarea);     

      cohort[BIOFUELS].prevPROD1.nitrogen *= ((double) cohort[BIOFUELS].prevchrtarea 
                                                   / (double) cohort[BIOFUELS].chrtarea);     

      // 10-year Woody Products

      for( i10 = 0; i10 < 10; ++i10 )
      {
        cohort[BIOFUELS].initPROD10[i10].carbon *= ((double) cohort[BIOFUELS].prevchrtarea 
                                                    / (double) cohort[BIOFUELS].chrtarea);     

        cohort[BIOFUELS].initPROD10[i10].nitrogen *= ((double) cohort[BIOFUELS].prevchrtarea 
                                                      / (double) cohort[BIOFUELS].chrtarea);   
      }

      cohort[BIOFUELS].prevPROD10.carbon *= ((double) cohort[BIOFUELS].prevchrtarea 
                                           / (double) cohort[BIOFUELS].chrtarea); 

      cohort[BIOFUELS].prevPROD10.nitrogen *= ((double) cohort[BIOFUELS].prevchrtarea 
                                             / (double) cohort[BIOFUELS].chrtarea); 

      // 100-year Woody Products

      for( i100 = 0; i100 < 100; ++i100 )
      {
        cohort[BIOFUELS].initPROD100[i100].carbon *= ((double) cohort[BIOFUELS].prevchrtarea 
                                                      / (double) cohort[BIOFUELS].chrtarea);     

        cohort[BIOFUELS].initPROD100[i100].nitrogen *= ((double) cohort[BIOFUELS].prevchrtarea 
                                                        / (double) cohort[BIOFUELS].chrtarea);   
      }

      cohort[BIOFUELS].prevPROD100.carbon *= ((double) cohort[BIOFUELS].prevchrtarea 
                                            / (double) cohort[BIOFUELS].chrtarea);     

      cohort[BIOFUELS].prevPROD100.nitrogen *= ((double) cohort[BIOFUELS].prevchrtarea    
                                              / (double) cohort[BIOFUELS].chrtarea);     
    }


    if (0 != cohort[PASTURE].chrtarea && 0 != cohort[PASTURE].prevchrtarea )
    {
      // 10-year Woody Products

      for( i10 = 0; i10 < 10; ++i10 )
      {
        cohort[PASTURE].initPROD10[i10].carbon *= ((double) cohort[PASTURE].prevchrtarea 
                                                   / (double) cohort[PASTURE].chrtarea);     

        cohort[PASTURE].initPROD10[i10].nitrogen *= ((double) cohort[PASTURE].prevchrtarea 
                                                     / (double) cohort[PASTURE].chrtarea);   
      }

      cohort[PASTURE].prevPROD10.carbon *= ((double) cohort[PASTURE].prevchrtarea 
                                           / (double) cohort[PASTURE].chrtarea); 

      cohort[PASTURE].prevPROD10.nitrogen *= ((double) cohort[PASTURE].prevchrtarea 
                                             / (double) cohort[PASTURE].chrtarea); 

      // 100-year Woody Products

      for( i100 = 0; i100 < 100; ++i100 )
      {
        cohort[PASTURE].initPROD100[i100].carbon *= ((double) cohort[PASTURE].prevchrtarea 
                                                     / (double) cohort[PASTURE].chrtarea);     

        cohort[PASTURE].initPROD100[i100].nitrogen *= ( (double) cohort[PASTURE].prevchrtarea 
                                                       / (double) cohort[PASTURE].chrtarea);   
      }

      cohort[PASTURE].prevPROD100.carbon *= ((double) cohort[PASTURE].prevchrtarea 
                                            / (double) cohort[PASTURE].chrtarea);     

      cohort[PASTURE].prevPROD100.nitrogen *= ((double) cohort[PASTURE].prevchrtarea    
                                              / (double) cohort[PASTURE].chrtarea);     
    }


    // Assign current year's products and conversion fluxes to crops and pasture

    if( sumManagedGain > 0 )
    {
      if( pdyr < tem.ag.getPRODUCTYEAR() )
      {
        tem.ag.setPRODUCTYEAR( tem.ag.getPRODUCTYEAR() + pdyr );
      }
      else { tem.ag.setPRODUCTYEAR( pdyr ); }


      // Land conversion to food crops
     
      if( diffCropArea > 0 )
      {
        if( sumNaturalLoss > 0 )
        {
          propConvertedArea = (double) diffCropArea 
                              / (double) sumManagedGain;
        }
        else
        {
          propConvertedArea = ZERO;
        }


        createCohortProducts( CROPVEG,
                              propConvertedArea,
                              formPROD10,
                              formPROD100,
                              slash,
                              vconvrtflx,
                              sconvrtflx,
                              nvretent,
                              nsretent );

        if( cohort[CROPVEG].aggrowdd < ZERO )
        {
          cohort[CROPVEG].aggrowdd = ZERO;
        }
      }

      // Land conversion to biofuel 
 
      if( diffBiofuelArea > 0 )
      {
        if( sumNaturalLoss > 0 )
        {
          propConvertedArea = (double) diffBiofuelArea 
                              / (double) sumManagedGain;
        }
        else
        {
          propConvertedArea = ZERO;
        }

      
        createCohortProducts( BIOFUELS,
                              propConvertedArea,
                              formPROD10,
                              formPROD100,
                              slash,
                              vconvrtflx,
                              sconvrtflx,
                              nvretent,
                              nsretent );

        if( cohort[BIOFUELS].aggrowdd < ZERO )
        {
          cohort[BIOFUELS].aggrowdd = ZERO;
        }
      }

      // Land conversion to pasture
 
      if( diffPastureArea > 0 )
      {
        if( sumNaturalLoss > 0 )
        {
          propConvertedArea = (double) diffPastureArea 
                              / (double) sumManagedGain;
        }
        else
        {
          propConvertedArea = ZERO;
        }


        createCohortProducts( PASTURE,
                              propConvertedArea,
                              formPROD10,
                              formPROD100,
                              slash,
                              vconvrtflx,
                              sconvrtflx,
                              nvretent,
                              nsretent );

      }
    }
  

    // Update carbon and nitrogen stocks of cohorts to reflect
    //   changes in carbon and nitrogen storage resulting from 
    //   land-use change

    for( ichrt = 0; ichrt < maxcohorts; ++ichrt )
    {
      if( BAREGRND != cohort[ichrt].currentveg
          && GLACIERS != cohort[ichrt].currentveg
          && LAKES != cohort[ichrt].currentveg )
      {
        diffarea = cohort[ichrt].chrtarea - cohort[ichrt].prevchrtarea;
        initarea = cohort[ichrt].prevchrtarea; 

  	
        if( diffarea > 0 )
        {
          // Update food crop cohort

  	  if( CROPVEG == cohort[ichrt].currentveg )
  	  {
            if( sumNaturalLoss > 0 
                && sumManagedGain > 0 )
            {
              propConvertedArea = (double) diffCropArea 
                                  / (double) sumManagedGain;
            }
            else
            {
              propConvertedArea = ZERO;
            }

            if( sumManagedLoss > 0
                && sumManagedGain > 0 )
            {
              propManagedAbandoned = (double) (sumManagedLoss - sumNaturalGain)
                                     / (double) sumManagedLoss;
 
              if( propManagedAbandoned < ZERO )
              {
                propManagedAbandoned = ZERO;
              }

              propAbandonedArea = ((double) diffCropArea
                                  / (double) sumManagedGain)
                                  * propManagedAbandoned;
            }
            else
            {
              propAbandonedArea = ZERO;
            }
          }

          // Update Biofuel cohort

          else if( BIOFUELS == cohort[ichrt].currentveg )
          {
            if( sumNaturalLoss > 0 
                && sumManagedGain > 0 )
            {
              propConvertedArea = (double) diffBiofuelArea 
                                  / (double) sumManagedGain;
            }
            else
            {
              propConvertedArea = ZERO;
            }

            if( sumManagedLoss > 0
                && sumManagedGain > 0 )
            {
              propManagedAbandoned = (double) (sumManagedLoss - sumNaturalGain)
                                     / (double) sumManagedLoss;
 
              if( propManagedAbandoned < ZERO )
              {
                propManagedAbandoned = ZERO;
              }

              propAbandonedArea = ((double) diffBiofuelArea
                                  / (double) sumManagedGain)
                                  * propManagedAbandoned;
            }
            else
            {
              propAbandonedArea = ZERO;
            }
          }
       
          // Update Pasture cohort

          else if( PASTURE == cohort[ichrt].currentveg )
          {
            if( sumNaturalLoss > 0 
                && sumManagedGain > 0 )
            {
              propConvertedArea = (double) diffPastureArea 
                                  / (double) sumManagedGain;
            }
            else
            {
              propConvertedArea = ZERO;
            }

            if( sumManagedLoss > 0
                && sumManagedGain > 0 )
            {
              propManagedAbandoned = (double) (sumManagedLoss - sumNaturalGain)
                                     / (double) sumManagedLoss;
 
              if( propManagedAbandoned < ZERO )
              {
                propManagedAbandoned = ZERO;
              }

              propAbandonedArea = ((double) diffPastureArea
                                  / (double) sumManagedGain)
                                  * propManagedAbandoned;
            }
            else
            {
              propAbandonedArea = ZERO;
            }
          }

          // Update natural cohorts

          else  
          {
            propConvertedArea = ZERO;

//            propAbandonedArea = (double) diffarea
//                                / (double) sumNaturalGain;

            propAbandonedArea = (double) diffarea
                                / (double) sumManagedLoss;
          }


          updateChangedCohort( ichrt,
                               initarea,
                               propConvertedArea,
                               propAbandonedArea,
                               lcluc.abandonedVeg,
                               lcluc.convertedSoil,
                               lcluc.abandonedSoil,
                               convertedNONSOLC,
                               abandonedNONSOLC,
                               convertedANH4IN,
                               abandonedANH4IN,
                               convertedANO3IN,
                               abandonedANO3IN,
                               convertedDPHUMIN,
                               abandonedDPHUMIN,
                               convertedOCIN,
                               abandonedOCIN,
                               convertedRCLIN,
                               abandonedRCLIN,
                               convertedRCRIN,
                               abandonedRCRIN,
                               convertedRCVLIN,
                               abandonedRCVLIN );                               


          if( cohort[ichrt].y[tem.I_VEGC] < ZERO )
          {
            cohort[ichrt].y[tem.I_VEGC] = ZERO;
          }

          if( cohort[ichrt].y[tem.I_STRN] < ZERO )
          {
            cohort[ichrt].y[tem.I_STRN] = ZERO;
          }

          if( cohort[ichrt].y[tem.I_STON] < ZERO )
          { 
            cohort[ichrt].y[tem.I_STON] = ZERO; 
          }

          if( cohort[ichrt].y[tem.I_SOLC] < ZERO )
          {
            cohort[ichrt].y[tem.I_SOLC] = ZERO; 
          }

          if( cohort[ichrt].y[tem.I_SOLN] < ZERO )
          {
            cohort[ichrt].y[tem.I_SOLN] = ZERO; 
          }

          if( cohort[ichrt].y[tem.I_AVLN] < ZERO )
          {
            cohort[ichrt].y[tem.I_AVLN] = ZERO; 
          }

          for( dlyr = 0; dlyr < NLVL; ++dlyr )
          {
            if( cohort[ichrt].NEManh4in[dlyr] < ZERO )
            {
              cohort[ichrt].NEManh4in[dlyr] = ZERO;
            }
        
            if( cohort[ichrt].NEMano3in[dlyr] < ZERO )
            {                            
              cohort[ichrt].NEMano3in[dlyr] = ZERO;
            }

            if( cohort[ichrt].NEMdphumin[dlyr] < ZERO )
            {                              
              cohort[ichrt].NEMdphumin[dlyr] = ZERO;
            }

            if( cohort[ichrt].NEMocin[dlyr] < ZERO )
            {                                    
              cohort[ichrt].NEMocin[dlyr] = ZERO;
            }

            if( cohort[ichrt].NEMrclin[dlyr] < ZERO )
            {
              cohort[ichrt].NEMrclin[dlyr] = ZERO;
            }

            if( cohort[ichrt].NEMrcrin[dlyr] < ZERO )
            {
              cohort[ichrt].NEMrcrin[dlyr] = ZERO;
            }

            if( cohort[ichrt].NEMrcvlin[dlyr] < ZERO )
            {
              cohort[ichrt].NEMrcvlin[dlyr] = ZERO;
            }
          }

          cohort[ichrt].prevy[tem.I_VEGC] = cohort[ichrt].y[tem.I_VEGC];

          cohort[ichrt].prevy[tem.I_STRN] = cohort[ichrt].y[tem.I_STRN];

          cohort[ichrt].prevy[tem.I_STON] = cohort[ichrt].y[tem.I_STON];
    
          cohort[ichrt].prevy[tem.I_SOLC] = cohort[ichrt].y[tem.I_SOLC];
 
          cohort[ichrt].prevy[tem.I_SOLN] = cohort[ichrt].y[tem.I_SOLN];

          cohort[ichrt].prevy[tem.I_AVLN] = cohort[ichrt].y[tem.I_AVLN];
        }
      }
    } 
  }

//  if( 5 == pdyr )
//  {
//    exit( -1 );
//  }
    
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITelmnt44::outputTEMmonth( const int& pchrt,
                                 const int& pdm )
{

  // Ecosystem carbon pools determined in integrator

  output[tem.I_VEGC][pchrt][pdm] = tem.getY( tem.I_VEGC );

  output[tem.I_SOLC][pchrt][pdm] = tem.getY( tem.I_SOLC );

  // Ecosystem nitrogen pools determined in integrator

  output[tem.I_STRN][pchrt][pdm] = tem.getY( tem.I_STRN );

  output[tem.I_STON][pchrt][pdm] = tem.getY( tem.I_STON );

  output[tem.I_SOLN][pchrt][pdm] = tem.getY( tem.I_SOLN );
      
  output[tem.I_AVLN][pchrt][pdm] = tem.getY( tem.I_AVLN );

 
  // Ecosystem water pools determined in integrator

  output[tem.I_AVLW][pchrt][pdm] = tem.soil.getAVLH2O();

  output[tem.I_SM][pchrt][pdm] = tem.soil.getMOIST();

  output[tem.I_VSM][pchrt][pdm] = tem.soil.getVSM();

  output[tem.I_PCTP][pchrt][pdm] = tem.soil.getPCTP();


  // Monthly phenology determined in integrator

  output[tem.I_UNRMLF][pchrt][pdm] = tem.getY( tem.I_UNRMLF );

  output[tem.I_LEAF][pchrt][pdm] = tem.getY( tem.I_LEAF );

  output[tem.I_LAI][pchrt][pdm] = tem.getY( tem.I_LAI );

  output[tem.I_FPC][pchrt][pdm] = tem.getY( tem.I_FPC );

  // Monthly carbon fluxes in ecosystems determined in integrator

  output[tem.I_INGPP][pchrt][pdm] = tem.getY( tem.I_INGPP );

  output[tem.I_GPP][pchrt][pdm] = tem.getY( tem.I_GPP );

  output[tem.I_FOZONE][pchrt][pdm] = tem.getY( tem.I_FOZONE );

  output[tem.I_FINDOZONE][pchrt][pdm] = tem.getY( tem.I_FINDOZONE );

  output[tem.I_INNPP][pchrt][pdm] = tem.getY( tem.I_INNPP );

  output[tem.I_NPP][pchrt][pdm] = tem.getY( tem.I_NPP );

  output[tem.I_GPR][pchrt][pdm] = tem.getY( tem.I_GPR );

  output[tem.I_RVMNT][pchrt][pdm] = tem.getY( tem.I_RVMNT );

  output[tem.I_RVGRW][pchrt][pdm] = tem.getY( tem.I_RVGRW );

  output[tem.I_LTRC][pchrt][pdm] = tem.getY( tem.I_LTRC );

  output[tem.I_RH][pchrt][pdm] = tem.getY( tem.I_RH );


  // Monthly nitrogen fluxes in ecosystems determined in 
  //   integrator

  output[tem.I_NINP][pchrt][pdm] = tem.getY( tem.I_NINP );

  output[tem.I_AGFRTN][pchrt][pdm] = tem.getY( tem.I_AGFRTN );

  output[tem.I_INNUP][pchrt][pdm] = tem.getY( tem.I_INNUP );

  output[tem.I_VNUP][pchrt][pdm] = tem.getY( tem.I_VNUP );

  output[tem.I_VSUP][pchrt][pdm] = tem.getY( tem.I_VSUP );

  output[tem.I_VLUP][pchrt][pdm] = tem.getY( tem.I_VLUP );

  output[tem.I_VNMBL][pchrt][pdm] = tem.getY( tem.I_VNMBL );

  output[tem.I_VNRSRB][pchrt][pdm] = tem.getY( tem.I_VNRSRB );

  output[tem.I_LTRN][pchrt][pdm] = tem.getY( tem.I_LTRN );

  output[tem.I_MNUP][pchrt][pdm] = tem.getY( tem.I_MNUP );

  output[tem.I_NMIN][pchrt][pdm] = tem.getY( tem.I_NMIN );

  output[tem.I_NLST][pchrt][pdm] = tem.getY( tem.I_NLST );


  // Monthly water fluxes in ecosystems

  output[tem.I_AGIRRIG][pchrt][pdm] = ZERO;

//  output[tem.I_INEET][pchrt][pdm] = tem.soil.getINEET();

  output[tem.I_EET][pchrt][pdm] = tem.soil.getEET();


  // Other ecosystem carbon pools

  output[tem.I_TOTEC][pchrt][pdm] = tem.ag.getTOTEC();

  output[tem.I_TOTC][pchrt][pdm] = tem.getTOTALC();

  // Other ecosystem nitrogen pools

  output[tem.I_VEGN][pchrt][pdm] = tem.veg.getVEGN();


  // Other ecosystem water pools
  
  output[tem.I_SNWPCK][pchrt][pdm] = tem.soil.getSNOWPACK();

  // Other monthly carbon fluxes in ecosystems

  output[tem.I_NEP][pchrt][pdm] = tem.getNEP();

  output[tem.I_NCE][pchrt][pdm] = tem.getNCE();


  // ******************* NEM or MDM ****************************
  
  output[tem.I_CH4EMS][pchrt][pdm] = tem.soil.getCH4EMISS();

  output[tem.I_CH4CSMP][pchrt][pdm] = tem.soil.getCH4CONSUMP();

  output[tem.I_CH4FLX][pchrt][pdm] = tem.soil.getCH4FLUX();


  // ********************* NEM *********************************
  
  output[tem.I_CO2NFLX][pchrt][pdm] = tem.soil.getCO2NFLUX();

  output[tem.I_CO2DNFLX][pchrt][pdm] = tem.soil.getCO2DNFLUX();

  output[tem.I_N2OFLX][pchrt][pdm] = tem.soil.getN2OFLUX();

  output[tem.I_N2ONFLX][pchrt][pdm] = tem.soil.getN2ONFLUX();

  output[tem.I_N2ODNFLX][pchrt][pdm] = tem.soil.getN2ODNFLUX();

  output[tem.I_N2FLX][pchrt][pdm] = tem.soil.getN2FLUX();

  
  // Other monthly water fluxes in ecosystems

  output[tem.I_PET][pchrt][pdm] = tem.atms.getPET();

  
  // Carbon in Human product pools

  output[tem.I_AGPRDC][pchrt][pdm] = tem.ag.getPROD1C();

  output[tem.I_PROD10C][pchrt][pdm] = tem.ag.getPROD10C();

  output[tem.I_PROD100C][pchrt][pdm] = tem.ag.getPROD100C();

  output[tem.I_TOTPRDC][pchrt][pdm] = tem.ag.getTOTPRODC();

  // Carbon in crop residue pool

  output[tem.I_RESIDC][pchrt][pdm] = tem.ag.getCROPRESIDUEC();

  output[tem.I_AGSTUBC][pchrt][pdm] = tem.ag.getSTUBBLEC();

  // Nitrogen in Human product pools

  output[tem.I_AGPRDN][pchrt][pdm] = tem.ag.getPROD1N();

  output[tem.I_PROD10N][pchrt][pdm] = tem.ag.getPROD10N();

  output[tem.I_PROD100N][pchrt][pdm] = tem.ag.getPROD100N();

  output[tem.I_TOTPRDN][pchrt][pdm] = tem.ag.getTOTPRODN();

  // Nitrogen in crop residue pool

  output[tem.I_RESIDN][pchrt][pdm] = tem.ag.getCROPRESIDUEN();

  output[tem.I_AGSTUBN][pchrt][pdm] = tem.ag.getSTUBBLEN();

  // Monthly carbon fluxes associated with
  //  agricultural conversion

  output[tem.I_CNVRTC][pchrt][pdm] = tem.ag.getCONVRTFLXC();

  output[tem.I_VCNVRTC][pchrt][pdm] = tem.ag.getVCONVRTFLXC();

  output[tem.I_SCNVRTC][pchrt][pdm] = tem.ag.getSCONVRTFLXC();

  output[tem.I_SLASHC][pchrt][pdm] = tem.ag.getSLASHC();

  output[tem.I_CFLX][pchrt][pdm] = tem.ag.getCFLUX();

  // Monthly nitrogen fluxes associated with
  //  agricultural conversion

  output[tem.I_CNVRTN][pchrt][pdm] = tem.ag.getCONVRTFLXN();

  output[tem.I_VCNVRTN][pchrt][pdm] = tem.ag.getVCONVRTFLXN();

  output[tem.I_SCNVRTN][pchrt][pdm] = tem.ag.getSCONVRTFLXN();

  output[tem.I_SLASHN][pchrt][pdm] = tem.ag.getSLASHN();

  output[tem.I_NRETNT][pchrt][pdm] = tem.ag.getNRETENT();

  output[tem.I_NVRTNT][pchrt][pdm] = tem.ag.getNVRETENT();

  output[tem.I_NSRTNT][pchrt][pdm] = tem.ag.getNSRETENT();

  // Monthly carbon and nitrogen fluxes from agricultural
  //   ecosystems

  output[tem.I_AGFPRDC][pchrt][pdm] = tem.ag.getCROPPRODC();
  output[tem.I_AGFPRDN][pchrt][pdm] = tem.ag.getCROPPRODN();

  output[tem.I_FRESIDC][pchrt][pdm] = tem.ag.getFORMCROPRESIDUEC();
  output[tem.I_FRESIDN][pchrt][pdm] = tem.ag.getFORMCROPRESIDUEN();

  output[tem.I_AGPRDFC][pchrt][pdm] = tem.ag.getPROD1DECAYC();
  output[tem.I_AGPRDFN][pchrt][pdm] = tem.ag.getPROD1DECAYN();

  output[tem.I_RESIDFC][pchrt][pdm] = tem.ag.getCROPRESIDUEFLXC();
  output[tem.I_RESIDFN][pchrt][pdm] = tem.ag.getCROPRESIDUEFLXN();


  // Monthly carbon and nitrogen fluxes from products

  output[tem.I_PRDF10C][pchrt][pdm] = tem.ag.getFORMPROD10C();
  output[tem.I_PRDF10N][pchrt][pdm] = tem.ag.getFORMPROD10N();
    	
  output[tem.I_PRD10FC][pchrt][pdm] = tem.ag.getPROD10DECAYC();
  output[tem.I_PRD10FN][pchrt][pdm] = tem.ag.getPROD10DECAYN();

  output[tem.I_PRDF100C][pchrt][pdm] = tem.ag.getFORMPROD100C();
  output[tem.I_PRDF100N][pchrt][pdm] = tem.ag.getFORMPROD100N();
  	
  output[tem.I_PRD100FC][pchrt][pdm] = tem.ag.getPROD100DECAYC();
  output[tem.I_PRD100FN][pchrt][pdm] = tem.ag.getPROD100DECAYN();

  output[tem.I_TOTFPRDC][pchrt][pdm] = tem.ag.getFORMTOTPRODC();
  output[tem.I_TOTFPRDN][pchrt][pdm] = tem.ag.getFORMTOTPRODN();
  	
  output[tem.I_TOTPRDFC][pchrt][pdm] = tem.ag.getTOTPRODDECAYC();
  output[tem.I_TOTPRDFN][pchrt][pdm] = tem.ag.getTOTPRODDECAYN();

  //  Output agricultural area-specific vs natural area-specific
  //    results
  
  if( 1 == tem.ag.state )
  {
    output[tem.I_CROPC][pchrt][pdm] = tem.getY( tem.I_VEGC );
    output[tem.I_NATVEGC][pchrt][pdm] = ZERO;

    output[tem.I_CROPN][pchrt][pdm] = tem.veg.getVEGN();
    output[tem.I_NATVEGN][pchrt][pdm] = ZERO;

    output[tem.I_CSTRN][pchrt][pdm] = tem.getY( tem.I_STRN );
    output[tem.I_NATSTRN][pchrt][pdm] = ZERO;

    output[tem.I_CSTON][pchrt][pdm] = tem.getY( tem.I_STON );
    output[tem.I_NATSTON][pchrt][pdm] = ZERO;

    output[tem.I_CROPULF][pchrt][pdm] = tem.getY( tem.I_UNRMLF );
    output[tem.I_NATULF][pchrt][pdm] = ZERO;

    output[tem.I_CROPLEAF][pchrt][pdm] = tem.getY( tem.I_LEAF );
    output[tem.I_NATLEAF][pchrt][pdm] = ZERO;

    output[tem.I_CROPLAI][pchrt][pdm] = tem.getY( tem.I_LAI );
    output[tem.I_NATLAI][pchrt][pdm] = ZERO;

    output[tem.I_CROPFPC][pchrt][pdm] = tem.getY( tem.I_FPC );
    output[tem.I_NATFPC][pchrt][pdm] = ZERO;

    output[tem.I_AGINGPP][pchrt][pdm] = tem.getY( tem.I_INGPP );
    output[tem.I_NATINGPP][pchrt][pdm] = ZERO;

    output[tem.I_AGGPP][pchrt][pdm] = tem.getY( tem.I_GPP );
    output[tem.I_NATGPP][pchrt][pdm] = ZERO;

    output[tem.I_AGINNPP][pchrt][pdm] = tem.getY( tem.I_INNPP );
    output[tem.I_NATINNPP][pchrt][pdm] = ZERO;

    output[tem.I_AGNPP][pchrt][pdm] = tem.getY( tem.I_NPP );
    output[tem.I_NATNPP][pchrt][pdm] = ZERO;

    output[tem.I_AGGPR][pchrt][pdm] = tem.getY( tem.I_GPR );
    output[tem.I_NATGPR][pchrt][pdm] = ZERO;

    output[tem.I_AGRVMNT][pchrt][pdm] = tem.getY( tem.I_RVMNT );
    output[tem.I_NATRVMNT][pchrt][pdm] = ZERO;

    output[tem.I_AGRVGRW][pchrt][pdm] = tem.getY( tem.I_RVGRW );
    output[tem.I_NATRVGRW][pchrt][pdm] = ZERO;

    output[tem.I_AGLTRC][pchrt][pdm] = tem.getY( tem.I_LTRC );
    output[tem.I_NATLTRC][pchrt][pdm] = ZERO;

    output[tem.I_AGINNUP][pchrt][pdm] = tem.getY( tem.I_INNUP );
    output[tem.I_NATINNUP][pchrt][pdm] = ZERO;

    output[tem.I_AGVNUP][pchrt][pdm] = tem.getY( tem.I_VNUP );
    output[tem.I_NATVNUP][pchrt][pdm] = ZERO;

    output[tem.I_AGVSUP][pchrt][pdm] = tem.getY( tem.I_VSUP );
    output[tem.I_NATVSUP][pchrt][pdm] = ZERO;

    output[tem.I_AGVLUP][pchrt][pdm] = tem.getY( tem.I_VLUP );
    output[tem.I_NATVLUP][pchrt][pdm] = ZERO;

    output[tem.I_AGVNMBL][pchrt][pdm] = tem.getY( tem.I_VNMBL );
    output[tem.I_NATVNMBL][pchrt][pdm] = ZERO;

    output[tem.I_AGVNRSRB][pchrt][pdm] = tem.getY( tem.I_VNRSRB );
    output[tem.I_NVNRSRB][pchrt][pdm] = ZERO;

    output[tem.I_AGLTRN][pchrt][pdm] = tem.getY( tem.I_LTRN );
    output[tem.I_NATLTRN][pchrt][pdm] = ZERO;
  }
  else
  {
    output[tem.I_CROPC][pchrt][pdm] = ZERO;
    output[tem.I_NATVEGC][pchrt][pdm] = tem.getY( tem.I_VEGC );

    output[tem.I_CROPN][pchrt][pdm] = ZERO;
    output[tem.I_NATVEGN][pchrt][pdm] = tem.veg.getVEGN();

    output[tem.I_CSTRN][pchrt][pdm] = ZERO;
    output[tem.I_NATSTRN][pchrt][pdm] = tem.getY( tem.I_STRN );

    output[tem.I_CSTON][pchrt][pdm] = ZERO;
    output[tem.I_NATSTON][pchrt][pdm] = tem.getY( tem.I_STON );

    output[tem.I_CROPULF][pchrt][pdm] = ZERO;
    output[tem.I_NATULF][pchrt][pdm] = tem.getY( tem.I_UNRMLF );

    output[tem.I_CROPLEAF][pchrt][pdm] = ZERO;
    output[tem.I_NATLEAF][pchrt][pdm] = tem.getY( tem.I_LEAF );

    output[tem.I_CROPLAI][pchrt][pdm] = ZERO;
    output[tem.I_NATLAI][pchrt][pdm] = tem.getY( tem.I_LAI );

    output[tem.I_CROPFPC][pchrt][pdm] = ZERO;
    output[tem.I_NATFPC][pchrt][pdm] = tem.getY( tem.I_FPC );

    output[tem.I_AGINGPP][pchrt][pdm] = ZERO;
    output[tem.I_NATINGPP][pchrt][pdm] = tem.getY( tem.I_INGPP );

    output[tem.I_AGGPP][pchrt][pdm] = ZERO;
    output[tem.I_NATGPP][pchrt][pdm] = tem.getY( tem.I_GPP );

    output[tem.I_AGINNPP][pchrt][pdm] = ZERO;
    output[tem.I_NATINNPP][pchrt][pdm] = tem.getY( tem.I_INNPP );

    output[tem.I_AGNPP][pchrt][pdm] = ZERO;
    output[tem.I_NATNPP][pchrt][pdm] = tem.getY( tem.I_NPP );

    output[tem.I_AGGPR][pchrt][pdm] = ZERO;
    output[tem.I_NATGPR][pchrt][pdm] = tem.getY( tem.I_GPR );

    output[tem.I_AGRVMNT][pchrt][pdm] = ZERO;
    output[tem.I_NATRVMNT][pchrt][pdm] = tem.getY( tem.I_RVMNT );

    output[tem.I_AGRVGRW][pchrt][pdm] = ZERO;
    output[tem.I_NATRVGRW][pchrt][pdm] = tem.getY( tem.I_RVGRW );

    output[tem.I_AGLTRC][pchrt][pdm] = ZERO;
    output[tem.I_NATLTRC][pchrt][pdm] = tem.getY( tem.I_RVGRW );

    output[tem.I_AGINNUP][pchrt][pdm] = ZERO;
    output[tem.I_NATINNUP][pchrt][pdm] = tem.getY( tem.I_INNUP );

    output[tem.I_AGVNUP][pchrt][pdm] = ZERO;
    output[tem.I_NATVNUP][pchrt][pdm] = tem.getY( tem.I_VNUP );

    output[tem.I_AGVSUP][pchrt][pdm] = ZERO;
    output[tem.I_NATVSUP][pchrt][pdm] = tem.getY( tem.I_VSUP );

    output[tem.I_AGVLUP][pchrt][pdm] = ZERO;
    output[tem.I_NATVLUP][pchrt][pdm] = tem.getY( tem.I_VLUP );

    output[tem.I_AGVNMBL][pchrt][pdm] = ZERO;
    output[tem.I_NATVNMBL][pchrt][pdm] = tem.getY( tem.I_VNMBL );

    output[tem.I_AGVNRSRB][pchrt][pdm] = ZERO;
    output[tem.I_NVNRSRB][pchrt][pdm] = tem.getY( tem.I_VNRSRB );

    output[tem.I_AGLTRN][pchrt][pdm] = ZERO;
    output[tem.I_NATLTRN][pchrt][pdm] = tem.getY( tem.I_LTRN );
  }


  // Climate variables
  
  output[tem.I_NIRR][pchrt][pdm] = tem.atms.getPAR(); 

  output[tem.I_PAR][pchrt][pdm] = tem.atms.getPAR(); 

  output[tem.I_TAIR][pchrt][pdm] = tem.atms.getTAIR(); 

  output[tem.I_PREC][pchrt][pdm] = tem.atms.getPREC(); 

  output[tem.I_SRFRUN][pchrt][pdm] = tem.soil.getSURFRUN(); 

  output[tem.I_DRAIN][pchrt][pdm] = tem.soil.getDRAINAGE();   

  output[tem.I_CO2][pchrt][pdm] = tem.atms.getCO2();

  output[tem.I_AOT40][pchrt][pdm] = tem.atms.getAOT40();  

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void MITelmnt44::readBinaryCohortState( ifstream& ifstate,
                                        const int& pichrt )
{
  
  int ichrt;

  float tstcol;
  
  float tstrow;
  
  ifstate.read( (char *)(&tstcol), sizeof( tstcol ) );
  ifstate.read( (char *)(&tstrow), sizeof( tstrow ) );

  if( (tstcol != col) || (tstrow != row) )
  {
    cout << "TEMSTATE is not co-registered with other data sets!" << endl;
    cout << "TEMSTATE Lon = " << tstcol << " Lat = " << tstrow << endl;
    cout << "Other Lon = " << col << " Lat = " << row << endl;

    exit( -1 );
  }

  ifstate.read( (char *)(&ichrt), sizeof( ichrt ) ); 
  ifstate.read( (char *)(&cohort[pichrt]), sizeof( MITElmntCohort44 ) );

  if( ifstate.fail() )
  {
  	cout << "Problem with reading in TEMSTATE at " << endl;
    cout << "Lon = " << tstcol << " Lat = " << tstrow;
    cout << " Cohort = " << ichrt << endl;
  
    exit( -1 );
  }
  	 	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void MITelmnt44::readCohortState( ifstream& ifstate, 
                                  const int& pichrt )
{
  int i;
  int dlyr;
  int dm;
  double dumdouble;
  int dumint;
  
  float tstcol;
  float tstrow;

  ifstate >> tstcol;    // Longitude of element
  ifstate >> tstrow;    // Latitude of element

  ifstate.seekg( 0, ios::cur );

  if( (tstcol != col) || (tstrow != row) )
  {
    cout << "TEMSTATE is not co-registered with other data sets!" << endl;
    cout << "TEMSTATE Lon = " << tstcol << " Lat = " << tstrow << endl;
    cout << "Other Lon = " << col << " Lat = " << row << endl;

    exit( -1 );
  }

  ifstate >> dumint;   // ichrt+1
  
  ifstate >> cohort[pichrt].srcCohort;
  
  ifstate >> cohort[pichrt].standage;
  
  ifstate >> cohort[pichrt].chrtarea;
  
  ifstate >> cohort[pichrt].prevchrtarea;

  ifstate >> cohort[pichrt].potveg;
  
  ifstate >> cohort[pichrt].currentveg;
  
  ifstate >> cohort[pichrt].subtype;
  
  ifstate >> cohort[pichrt].cmnt;

 	
  for( i = 0; i < MAXSTATE; ++i )
  {
    ifstate >> cohort[pichrt].y[i];
    ifstate >> cohort[pichrt].prevy[i];
  }

  ifstate >> cohort[pichrt].agcmnt;

  ifstate >> cohort[pichrt].aggrowdd;

  ifstate >> cohort[pichrt].agkd;

  ifstate >> dumint;

//  ifstate >> cohort[pichrt].agprvstate;

  ifstate >> dumint;

//  ifstate >> cohort[pichrt].agstate;

  ifstate >> cohort[pichrt].c2n;
  
  ifstate >> cohort[pichrt].cneven;

  ifstate >> dumdouble;

//  ifstate >> cohort[pichrt].convrtflx.carbon;

  ifstate >> dumdouble;

//  ifstate >> cohort[pichrt].convrtflx.nitrogen;

  ifstate >> cohort[pichrt].cropprveetmx;

  ifstate >> cohort[pichrt].cropprvleafmx;

  ifstate >> cohort[pichrt].cropprvpetmx;

  ifstate >> cohort[pichrt].cropResidue.carbon;
  ifstate >> cohort[pichrt].cropResidue.nitrogen;

  ifstate >> cohort[pichrt].croptopt;

  ifstate >> cohort[pichrt].distmnthcnt;

  ifstate >> cohort[pichrt].disturbflag;
  
  ifstate >> cohort[pichrt].disturbmonth;
  
  ifstate >> cohort[pichrt].eetmx;

  ifstate >> dumint;
//  ifstate >> cohort[pichrt].fertflag;                              

  ifstate >> cohort[pichrt].firemnthcnt;
  
  ifstate >> cohort[pichrt].firendep;
  
  ifstate >> cohort[pichrt].formPROD10.carbon;
  ifstate >> cohort[pichrt].formPROD10.nitrogen;

  ifstate >> cohort[pichrt].formPROD100.carbon;
  ifstate >> cohort[pichrt].formPROD100.nitrogen;

  ifstate >> cohort[pichrt].fprevozone;

  ifstate >> cohort[pichrt].FRI;

  for( dm = 0; dm < CYCLE; ++dm )
  {  
    ifstate >> cohort[pichrt].initPROD1[dm].carbon;
    ifstate >> cohort[pichrt].initPROD1[dm].nitrogen;
  }
  
  for( i = 0; i < 10; ++i )
  {
    ifstate >> cohort[pichrt].initPROD10[i].carbon;
    ifstate >> cohort[pichrt].initPROD10[i].nitrogen;
  }
    
  for( i = 0; i < 100; ++i )
  {
    ifstate >> cohort[pichrt].initPROD100[i].carbon;
    ifstate >> cohort[pichrt].initPROD100[i].nitrogen;
  }

  ifstate >> cohort[pichrt].irrgflag;                              
  
  ifstate >> cohort[pichrt].kd;

  ifstate >> cohort[pichrt].MDMnpp;       // MDM

  ifstate >> cohort[pichrt].natprveetmx;

  ifstate >> cohort[pichrt].natprvleafmx;

  ifstate >> cohort[pichrt].natprvpetmx;

  ifstate >> cohort[pichrt].natseedC;

  ifstate >> cohort[pichrt].natseedSTRN;

  ifstate >> cohort[pichrt].natseedSTON;

  ifstate >> cohort[pichrt].natsoil;

  ifstate >> cohort[pichrt].nattopt;

  ifstate >> cohort[pichrt].natyreet;

  ifstate >> cohort[pichrt].natyrpet;


  // ************************* NEM *****************************
  
  for( dlyr = 0; dlyr < NLVL; ++dlyr )
  {
    ifstate >> cohort[pichrt].NEManh4in[dlyr]; 
  
    ifstate >> cohort[pichrt].NEMano3in[dlyr];

    ifstate >> cohort[pichrt].NEMdphumin[dlyr]; 

    ifstate >> cohort[pichrt].NEMocin[dlyr];

    ifstate >> cohort[pichrt].NEMrclin[dlyr]; 
  
    ifstate >> cohort[pichrt].NEMrcrin[dlyr];
  
    ifstate >> cohort[pichrt].NEMrcvlin[dlyr];
  }

  ifstate >> cohort[pichrt].NEMnsolc;
  
//  ifstate >> cohort[pichrt].NEMtopdens;
  
//  ifstate >> cohort[pichrt].NEMtopksat;
  
//  ifstate >> cohort[pichrt].NEMtoppor;

  // ***********************************************************
  
  
  ifstate >> cohort[pichrt].newleafmx;

  ifstate >> cohort[pichrt].newtopt;

  ifstate >> dumdouble;
//  ifstate >> cohort[pichrt].nretent;

  ifstate >> dumdouble;
//  ifstate >> cohort[pichrt].nsretent;

  ifstate >> dumdouble;
//  ifstate >> cohort[pichrt].nvretent;

  ifstate >> cohort[pichrt].petmx;

  ifstate >> cohort[pichrt].prev2tair;

  ifstate >> cohort[pichrt].prevco2;

  ifstate >> cohort[pichrt].prevCropResidue.carbon;
  ifstate >> cohort[pichrt].prevCropResidue.nitrogen;

  ifstate >> cohort[pichrt].prevPROD1.carbon;
  ifstate >> cohort[pichrt].prevPROD1.nitrogen;

  ifstate >> cohort[pichrt].prevPROD10.carbon;
  ifstate >> cohort[pichrt].prevPROD10.nitrogen;

  ifstate >> cohort[pichrt].prevPROD100.carbon;
  ifstate >> cohort[pichrt].prevPROD100.nitrogen;

//  ifstate >> cohort[pichrt].prevspack;

  ifstate >> cohort[pichrt].prevtair;

  ifstate >> cohort[pichrt].prevunrmleaf;
  
//  ifstate >> dumdouble;
//  ifstate >> cohort[pichrt].prod10par; 

//  ifstate >> dumdouble;
//  ifstate >> cohort[pichrt].prod100par; 

  ifstate >> cohort[pichrt].productYear;

  ifstate >> cohort[pichrt].prvcropnpp;

  ifstate >> cohort[pichrt].prveetmx;

  ifstate >> cohort[pichrt].prvleafmx;

  ifstate >> cohort[pichrt].prvpetmx;

  ifstate >> cohort[pichrt].qc;

//  ifstate >> dumdouble;
//  ifstate >> cohort[pichrt].sconvert; 
  
  ifstate >> dumdouble;
//  ifstate >> cohort[pichrt].sconvrtflx.carbon;
  ifstate >> dumdouble;
//  ifstate >> cohort[pichrt].sconvrtflx.nitrogen;

  ifstate >> dumdouble;
//  ifstate >> cohort[pichrt].slash.carbon;

  ifstate >> dumdouble;
//  ifstate >> cohort[pichrt].slash.nitrogen;

//  ifstate >> dumdouble;
//  ifstate >> cohort[pichrt].slashpar; 

  ifstate >> cohort[pichrt].tillflag;                           

  ifstate >> cohort[pichrt].topt;

  ifstate >> cohort[pichrt].tqc;

//  ifstate >> dumdouble;
//  ifstate >> cohort[pichrt].vconvert; 

  ifstate >> dumdouble;
//  ifstate >> cohort[pichrt].vconvrtflx.carbon;
  ifstate >> dumdouble;
//  ifstate >> cohort[pichrt].vconvrtflx.nitrogen;

//  ifstate >> dumdouble;
//  ifstate >> cohort[pichrt].vrespar; 

  ifstate >> cohort[pichrt].yragstubC;
  
  ifstate >> cohort[pichrt].yragstubN;

  ifstate >> cohort[pichrt].yrcflux;
  
  ifstate >> cohort[pichrt].yrCH4csmp;
  
  ifstate >> cohort[pichrt].yrCH4ems;
  
  ifstate >> cohort[pichrt].yrCH4flx;
  
  ifstate >> cohort[pichrt].yrCO2dnflx;
  
  ifstate >> cohort[pichrt].yrCO2nflx;
  
  ifstate >> cohort[pichrt].yrconvrtC;
  
  ifstate >> cohort[pichrt].yrconvrtN;
  
  ifstate >> cohort[pichrt].yrdecayPROD1C;
  
  ifstate >> cohort[pichrt].yrdecayPROD10C;
  
  ifstate >> cohort[pichrt].yrdecayPROD100C;

  ifstate >> cohort[pichrt].yrdecayPROD1N;
  
  ifstate >> cohort[pichrt].yrdecayPROD10N;
  
  ifstate >> cohort[pichrt].yrdecayPROD100N;
  
  ifstate >> cohort[pichrt].yrdecayTOTPRODC;

  ifstate >> cohort[pichrt].yrdecayTOTPRODN;

  ifstate >> cohort[pichrt].yreet;
  
  ifstate >> cohort[pichrt].yrfertn;
  
  ifstate >> cohort[pichrt].yrfluxResidueC;
  
  ifstate >> cohort[pichrt].yrfluxResidueN;

  ifstate >> cohort[pichrt].yrformPROD1C;
  
  ifstate >> cohort[pichrt].yrformPROD10C;
  
  ifstate >> cohort[pichrt].yrformPROD100C;

  ifstate >> cohort[pichrt].yrformPROD1N;
  
  ifstate >> cohort[pichrt].yrformPROD10N;
  
  ifstate >> cohort[pichrt].yrformPROD100N;
  
  ifstate >> cohort[pichrt].yrformResidueC;
  
  ifstate >> cohort[pichrt].yrformResidueN;

  ifstate >> cohort[pichrt].yrformTOTPRODC;

  ifstate >> cohort[pichrt].yrformTOTPRODN;

  ifstate >> cohort[pichrt].yrfpc;

  ifstate >> cohort[pichrt].yrgpp;
  
  ifstate >> cohort[pichrt].yrgpr;

//  ifstate >> cohort[pichrt].yrh2oyld;
  
  ifstate >> cohort[pichrt].yrimmob;
  
//  ifstate >> cohort[pichrt].yrineet;
  
  ifstate >> cohort[pichrt].yringpp;
  
  ifstate >> cohort[pichrt].yrinnpp;

//  ifstate >> cohort[pichrt].yrirrig;
    
  ifstate >> cohort[pichrt].yrlai;

  ifstate >> cohort[pichrt].yrleaf;

  ifstate >> cohort[pichrt].yrltrfalC;

  ifstate >> cohort[pichrt].yrltrfalN;  	
  
  ifstate >> cohort[pichrt].yrN2flx;
  
  ifstate >> cohort[pichrt].yrN2Odnflx;
  
  ifstate >> cohort[pichrt].yrN2Oflx;
  
  ifstate >> cohort[pichrt].yrN2Onflx;
  
  ifstate >> cohort[pichrt].yrnce;
  
  ifstate >> cohort[pichrt].yrnep;
  
  ifstate >> cohort[pichrt].yrninput;
  
  ifstate >> cohort[pichrt].yrnlost;
  
  ifstate >> cohort[pichrt].yrnmin;
    
  ifstate >> cohort[pichrt].yrnpp;
  
  ifstate >> cohort[pichrt].yrnrent;
    
  ifstate >> cohort[pichrt].yrnsrent;

  ifstate >> cohort[pichrt].yrnvrent;
  
  ifstate >> cohort[pichrt].yrpet;
  
//  ifstate >> cohort[pichrt].yrrgrndH2O;  
  
  ifstate >> cohort[pichrt].yrrh;
    
//  ifstate >> cohort[pichrt].yrrain;
  
//  ifstate >> cohort[pichrt].yrrperc;
  
//  ifstate >> cohort[pichrt].yrrrun;

  ifstate >> cohort[pichrt].yrsconvrtC;

  ifstate >> cohort[pichrt].yrsconvrtN;

//  ifstate >> cohort[pichrt].yrsgrndh2o;
  
  ifstate >> cohort[pichrt].yrslashC;

  ifstate >> cohort[pichrt].yrslashN;
  
//  ifstate >> cohort[pichrt].yrsnowfall;

//  ifstate >> cohort[pichrt].yrsnowinf;

//  ifstate >> cohort[pichrt].yrsnowpack;

//  ifstate >> cohort[pichrt].yrsoilavlH2O;

  ifstate >> cohort[pichrt].yrsoilavlN;
  
//  ifstate >> cohort[pichrt].yrsoilC2N; 
  
//  ifstate >> cohort[pichrt].yrsoilmoist;

  ifstate >> cohort[pichrt].yrsoilorgC;
  
  ifstate >> cohort[pichrt].yrsoilorgN;
  
//  ifstate >> cohort[pichrt].yrsoilpctp;

//  ifstate >> cohort[pichrt].yrsoilvsm;

//  ifstate >> cohort[pichrt].yrsperc;
  
//  ifstate >> cohort[pichrt].yrsrun;

  ifstate >> cohort[pichrt].yrtotalC;
  
  ifstate >> cohort[pichrt].yrunleaf;

  ifstate >> cohort[pichrt].yrvconvrtC;
  
  ifstate >> cohort[pichrt].yrvconvrtN;

  ifstate >> cohort[pichrt].yrvegC;

//  ifstate >> cohort[pichrt].yrC2N;

  ifstate >> cohort[pichrt].yrveginnup;
    
  ifstate >> cohort[pichrt].yrveglup;

  ifstate >> cohort[pichrt].yrvegN;
  
  ifstate >> cohort[pichrt].yrvegnmobil;

  ifstate >> cohort[pichrt].yrvegnrsorb;

  ifstate >> cohort[pichrt].yrvegnup;

  ifstate >> cohort[pichrt].yrvegrgrowth;

  ifstate >> cohort[pichrt].yrvegrmaint;

  ifstate >> cohort[pichrt].yrvegSTON;
  
  ifstate >> cohort[pichrt].yrvegSTRN;

  ifstate >> cohort[pichrt].yrvegsup;

  ifstate.seekg( 0, ios::cur );
    	
};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */
     
void MITelmnt44::saveTEMCohortState( const int& pichrt )
{
  int dlyr;   // NEM
  int dm;
  int i;
 
 //  cohort[pichrt].srcCohort = cohort[pichrt].srcCohort

//  cohort[pichrt].standage = cohort[pichrt].standage

//  cohort[pichrt].chrtarea = cohort[pichrt].chrtarea
 
//  cohort[pichrt].prevchrtarea = cohort[pichrt].prevchrtarea
  
  cohort[pichrt].potveg = tem.veg.getPOTVEG();

  cohort[pichrt].currentveg = tem.veg.getCURRENTVEG();

  cohort[pichrt].subtype = tem.veg.getSUBTYPE();

  cohort[pichrt].cmnt = tem.veg.cmnt;

  for( i = 0; i < MAXSTATE; ++i )
  {
    cohort[pichrt].y[i] = tem.getY( i );
    cohort[pichrt].prevy[i] = tem.getPREVY( i );
  }

  cohort[pichrt].agcmnt = tem.ag.cmnt;

  cohort[pichrt].aggrowdd = tem.ag.getGROWDD();

  cohort[pichrt].agkd = tem.ag.getKD();

  cohort[pichrt].agprvstate = tem.ag.prvstate;

  cohort[pichrt].agstate = tem.ag.state;

  cohort[pichrt].c2n = tem.veg.getC2N();
  
  cohort[pichrt].cneven = tem.veg.getCNEVEN();

  cohort[pichrt].convrtflx.carbon = tem.ag.getCONVRTFLXC();
  cohort[pichrt].convrtflx.nitrogen = tem.ag.getCONVRTFLXN();

  cohort[pichrt].cropprveetmx = tem.ag.getCROPPRVEETMX();

  cohort[pichrt].cropprvleafmx = tem.ag.getCROPPRVLEAFMX();

  cohort[pichrt].cropprvpetmx = tem.ag.getCROPPRVPETMX();

  cohort[pichrt].cropResidue.carbon = tem.ag.getCROPRESIDUEC();
  cohort[pichrt].cropResidue.nitrogen = tem.ag.getCROPRESIDUEN();

  cohort[pichrt].croptopt = tem.ag.getCROPTOPT();

  cohort[pichrt].distmnthcnt = tem.distmnthcnt;

  cohort[pichrt].disturbflag = tem.disturbflag;
  
  cohort[pichrt].disturbmonth = tem.disturbmonth;
  
  cohort[pichrt].eetmx = tem.soil.getEETMX();

  cohort[pichrt].fertflag = tem.ag.fertflag;                              

  cohort[pichrt].firemnthcnt = tem.firemnthcnt;
  
  cohort[pichrt].firendep = tem.ag.getFIRENDEP();
  
  cohort[pichrt].formPROD10.carbon = tem.ag.getFORMPROD10C();
  cohort[pichrt].formPROD10.nitrogen = tem.ag.getFORMPROD10N();

  cohort[pichrt].formPROD100.carbon = tem.ag.getFORMPROD100C();
  cohort[pichrt].formPROD100.nitrogen = tem.ag.getFORMPROD100N();

  cohort[pichrt].fprevozone = tem.veg.getFPREVOZONE();

  cohort[pichrt].FRI = tem.ag.getFRI();

  for( dm = 0; dm < CYCLE; ++dm )
  {  
    cohort[pichrt].initPROD1[dm].carbon = tem.ag.getINITPROD1C( dm );
    cohort[pichrt].initPROD1[dm].nitrogen = tem.ag.getINITPROD1N( dm );
  }
  
  for( i = 0; i < 10; ++i )
  {
    cohort[pichrt].initPROD10[i].carbon = tem.ag.getINITPROD10C( i );
    cohort[pichrt].initPROD10[i].nitrogen = tem.ag.getINITPROD10N( i );
  }
    
  for( i = 0; i < 100; ++i )
  {
    cohort[pichrt].initPROD100[i].carbon = tem.ag.getINITPROD100C( i );
    cohort[pichrt].initPROD100[i].nitrogen = tem.ag.getINITPROD100N( i );
  }

  cohort[pichrt].irrgflag = tem.ag.irrgflag;                              
  
  cohort[pichrt].kd = tem.microbe.getKD();

  cohort[pichrt].MDMnpp = tem.getY( tem.I_NPP );       // MDM

  cohort[pichrt].natprveetmx = tem.ag.getNATPRVEETMX();

  cohort[pichrt].natprvleafmx = tem.ag.getNATPRVLEAFMX();

  cohort[pichrt].natprvpetmx = tem.ag.getNATPRVPETMX();

  cohort[pichrt].natseedC = tem.ag.getNATSEEDC();

  cohort[pichrt].natseedSTRN = tem.ag.getNATSEEDSTRN();

  cohort[pichrt].natseedSTON = tem.ag.getNATSEEDSTON();

  cohort[pichrt].natsoil = tem.ag.getNATSOIL();

  cohort[pichrt].nattopt = tem.ag.getNATTOPT();

  cohort[pichrt].natyreet = tem.ag.getNATYREET();

  cohort[pichrt].natyrpet = tem.ag.getNATYRPET();


  // ************************** NEM ****************************
  
  for( dlyr = 0; dlyr < NLVL; ++dlyr )
  {
    cohort[pichrt].NEManh4in[dlyr] = tem.microbe.nem.getANH4IN( dlyr ); 
  
    cohort[pichrt].NEMano3in[dlyr] = tem.microbe.nem.getANO3IN( dlyr );

    cohort[pichrt].NEMdphumin[dlyr] = tem.microbe.nem.getDPHUMIN( dlyr ); 

    cohort[pichrt].NEMocin[dlyr] = tem.microbe.nem.getOCIN( dlyr );

    cohort[pichrt].NEMrclin[dlyr] = tem.microbe.nem.getRCLIN( dlyr ); 
  
    cohort[pichrt].NEMrcrin[dlyr] = tem.microbe.nem.getRCRIN( dlyr );
  
    cohort[pichrt].NEMrcvlin[dlyr] = tem.microbe.nem.getRCVLIN( dlyr );
  }
  
  cohort[pichrt].NEMnsolc = tem.microbe.nem.getNONSOLC();

//  cohort[pichrt].NEMtopdens = tem.microbe.nem.getTOPDENS();
  
//  cohort[pichrt].NEMtopksat = tem.microbe.nem.getTOPKSAT();
  
//  cohort[pichrt].NEMtoppor = tem.microbe.nem.getTOPPOR();

  // ***********************************************************
  
  cohort[pichrt].newleafmx = tem.veg.getNEWLEAFMX();

  cohort[pichrt].newtopt = tem.veg.getNEWTOPT();

  cohort[pichrt].nretent = tem.ag.getNRETENT();

  cohort[pichrt].nsretent = tem.ag.getNSRETENT();

  cohort[pichrt].nvretent = tem.ag.getNVRETENT();

  cohort[pichrt].petmx = tem.atms.getPETMX();

  cohort[pichrt].prev2tair = tem.atms.getPREV2TAIR();

  cohort[pichrt].prevco2 = tem.atms.getPREVCO2();

  cohort[pichrt].prevCropResidue.carbon = tem.ag.getPREVCROPRESIDUEC();
  cohort[pichrt].prevCropResidue.nitrogen = tem.ag.getPREVCROPRESIDUEN();

  cohort[pichrt].prevPROD1.carbon = tem.ag.getPREVPROD1C();
  cohort[pichrt].prevPROD1.nitrogen = tem.ag.getPREVPROD1N();

  cohort[pichrt].prevPROD10.carbon = tem.ag.getPREVPROD10C();
  cohort[pichrt].prevPROD10.nitrogen = tem.ag.getPREVPROD10N();

  cohort[pichrt].prevPROD100.carbon = tem.ag.getPREVPROD100C();
  cohort[pichrt].prevPROD100.nitrogen = tem.ag.getPREVPROD100N();

//  cohort[pichrt].prevspack = tem.soil.getPREVSPACK();

  cohort[pichrt].prevtair = tem.atms.getPREVTAIR();

  cohort[pichrt].prevunrmleaf = tem.veg.getPREVUNRMLEAF();
  
//  cohort[pichrt].prod10par = tem.ag.getPROD10PAR(); 

//  cohort[pichrt].prod100par = tem.ag.getPROD100PAR(); 

  cohort[pichrt].productYear = tem.ag.getPRODUCTYEAR();
  
  cohort[pichrt].prvcropnpp = tem.ag.getPRVCROPNPP();

  cohort[pichrt].prveetmx = tem.soil.getPRVEETMX();

  cohort[pichrt].prvleafmx = tem.veg.getPRVLEAFMX();

  cohort[pichrt].prvpetmx = tem.atms.getPRVPETMX();

// cohort[pichrt].qc determined in setTEMequilState() 

//  cohort[pichrt].sconvert = tem.ag.getSCONVERT(); 
  
  cohort[pichrt].sconvrtflx.carbon = tem.ag.getSCONVRTFLXC();
  cohort[pichrt].sconvrtflx.nitrogen = tem.ag.getSCONVRTFLXN();

  cohort[pichrt].slash.carbon = tem.ag.getSLASHC();
  cohort[pichrt].slash.nitrogen = tem.ag.getSLASHN();

//  cohort[pichrt].slashpar = tem.ag.getSLASHPAR(); 
 
  cohort[pichrt].tillflag = tem.ag.tillflag;                           

  cohort[pichrt].topt = tem.veg.getTOPT();

// cohort[pichrt].tqc determined in setTEMequilState()

//  cohort[pichrt].vconvert = tem.ag.getVCONVERT(); 

  cohort[pichrt].vconvrtflx.carbon = tem.ag.getVCONVRTFLXC();
  cohort[pichrt].vconvrtflx.nitrogen = tem.ag.getVCONVRTFLXN();

//  cohort[pichrt].vrespar = tem.ag.getVRESPAR(); 

  cohort[pichrt].yragstubC = tem.ag.getYRSTUBC();
  
  cohort[pichrt].yragstubN = tem.ag.getYRSTUBN();

  cohort[pichrt].yrcflux = tem.ag.getYRCFLUX();
  
  cohort[pichrt].yrCH4csmp = tem.soil.getYRCH4CSMP();
  
  cohort[pichrt].yrCH4ems = tem.soil.getYRCH4EMISSION();
  
  cohort[pichrt].yrCH4flx = tem.soil.getYRCH4FLUX();
  
  cohort[pichrt].yrCO2dnflx = tem.soil.getYRCO2DENTRFLUX();
  
  cohort[pichrt].yrCO2nflx = tem.soil.getYRCO2NTRFLUX();
  
  cohort[pichrt].yrconvrtC = tem.ag.getYRCONVRTC();
  
  cohort[pichrt].yrconvrtN = tem.ag.getYRCONVRTN();
  
  cohort[pichrt].yrdecayPROD1C = tem.ag.getYRDECAYPROD1C();
  
  cohort[pichrt].yrdecayPROD10C = tem.ag.getYRDECAYPROD10C();
  
  cohort[pichrt].yrdecayPROD100C = tem.ag.getYRDECAYPROD100C();

  cohort[pichrt].yrdecayPROD1N = tem.ag.getYRDECAYPROD1N();
  
  cohort[pichrt].yrdecayPROD10N = tem.ag.getYRDECAYPROD10N();
  
  cohort[pichrt].yrdecayPROD100N = tem.ag.getYRDECAYPROD100N();
  
  cohort[pichrt].yrdecayTOTPRODC = tem.ag.getYRDECAYTOTPRODC();

  cohort[pichrt].yrdecayTOTPRODN = tem.ag.getYRDECAYTOTPRODN();

  cohort[pichrt].yreet = tem.soil.getYREET();
  
  cohort[pichrt].yrfertn = tem.ag.getYRFERTN();
  
  cohort[pichrt].yrfluxResidueC = tem.ag.getYRFLUXRESIDUEC();
  
  cohort[pichrt].yrfluxResidueN = tem.ag.getYRFLUXRESIDUEN();

  cohort[pichrt].yrformPROD1C = tem.ag.getYRFORMPROD1C();
  
  cohort[pichrt].yrformPROD10C = tem.ag.getYRFORMPROD10C();
  
  cohort[pichrt].yrformPROD100C = tem.ag.getYRFORMPROD100C();

  cohort[pichrt].yrformPROD1N = tem.ag.getYRFORMPROD1N();
  
  cohort[pichrt].yrformPROD10N = tem.ag.getYRFORMPROD10N();
  
  cohort[pichrt].yrformPROD100N = tem.ag.getYRFORMPROD100N();
  
  cohort[pichrt].yrformResidueC = tem.ag.getYRFORMRESIDUEC();
  
  cohort[pichrt].yrformResidueN = tem.ag.getYRFORMRESIDUEN();

  cohort[pichrt].yrformTOTPRODC = tem.ag.getYRFORMTOTPRODC();

  cohort[pichrt].yrformTOTPRODN = tem.ag.getYRFORMTOTPRODN();

  cohort[pichrt].yrfpc = tem.veg.getYRFPC();

  cohort[pichrt].yrgpp = tem.veg.getYRGPP();
  
  cohort[pichrt].yrgpr = tem.veg.getYRGPR();
  
//  cohort[pichrt].yrh2oyld = tem.soil.getYRH2OYIELD();
  
  cohort[pichrt].yrimmob = tem.microbe.getYRNUPTAKE();
  
//  cohort[pichrt].yrineet = tem.soil.getYRINEET();
  
  cohort[pichrt].yringpp = tem.veg.getYRINGPP();
  
  cohort[pichrt].yrinnpp = tem.veg.getYRINNPP();
      
//  cohort[pichrt].yrirrig = tem.ag.getYRIRRIG();

  cohort[pichrt].yrlai = tem.veg.getYRLAI();
    
  cohort[pichrt].yrleaf = tem.veg.getYRLEAF();

  cohort[pichrt].yrltrfalC = tem.veg.getYRLTRFALC();

  cohort[pichrt].yrltrfalN = tem.veg.getYRLTRFALN();  	
  
  cohort[pichrt].yrN2flx = tem.soil.getYRN2FLUX();
  
  cohort[pichrt].yrN2Odnflx = tem.soil.getYRN2ODENTRFLUX();
  
  cohort[pichrt].yrN2Oflx = tem.soil.getYRN2OFLUX();
  
  cohort[pichrt].yrN2Onflx = tem.soil.getYRN2ONTRFLUX();
  
  cohort[pichrt].yrnce = tem.getYRNCE();
  
  cohort[pichrt].yrnep = tem.getYRNEP();
  
  cohort[pichrt].yrninput = tem.soil.getYRNINPUT();
  
  cohort[pichrt].yrnlost = tem.soil.getYRNLOST();
  
  cohort[pichrt].yrnmin = tem.microbe.getYRNMIN();
   
  cohort[pichrt].yrnpp = tem.veg.getYRNPP();
  
  cohort[pichrt].yrnrent = tem.ag.getYRNRENT(); 
  
  cohort[pichrt].yrnsrent = tem.ag.getYRNSRENT();

  cohort[pichrt].yrnvrent = tem.ag.getYRNVRENT();
  
  cohort[pichrt].yrpet = tem.atms.getYRPET();
  
//  cohort[pichrt].yrrgrndH2O = tem.soil.getYRRGRNDH2O(); 
  
  cohort[pichrt].yrrh = tem.microbe.getYRRH();
    
//  cohort[pichrt].yrrain = tem.atms.getYRRAIN();
  
//  cohort[pichrt].yrrperc = tem.soil.getYRRPERC();
  
//  cohort[pichrt].yrrrun = tem.soil.getYRRRUN();

  cohort[pichrt].yrsconvrtC = tem.ag.getYRSCONVRTC();

  cohort[pichrt].yrsconvrtN = tem.ag.getYRSCONVRTN();
  
//  cohort[pichrt].yrsgrndH2O = tem.soil.getYRSGRNDH2O();

  cohort[pichrt].yrslashC = tem.ag.getYRSLASHC();

  cohort[pichrt].yrslashN = tem.ag.getYRSLASHN();
  
//  cohort[pichrt].yrsnowfall = tem.atms.getYRSNOWFALL();

//  cohort[pichrt].yrsnowinf = tem.soil.getYRSNOWINF();

//  cohort[pichrt].yrsnowpack = tem.soil.getYRSNOWPACK();

//  cohort[pichrt].yrsoilavlH2O = tem.soil.getYRAVLH2O();

  cohort[pichrt].yrsoilavlN = tem.soil.getYRAVLN();

//  cohort[pichrt].yrsoilC2N = tem.soil.getYRC2N();

//  cohort[pichrt].yrsoilmoist = tem.soil.getYRSMOIST();
  
  cohort[pichrt].yrsoilorgC = tem.soil.getYRORGC();
  
  cohort[pichrt].yrsoilorgN = tem.soil.getYRORGN();
  
//  cohort[pichrt].yrsoilpctp = tem.soil.getYRPCTP();

//  cohort[pichrt].yrsoilvsm = tem.soil.getYRVSM();

//  cohort[pichrt].yrsperc = tem.soil.getYRSPERC();
  
//  cohort[pichrt].yrsrun = tem.soil.getYRSRUN();

  cohort[pichrt].yrtotalC = tem.getYRTOTALC();
  
  cohort[pichrt].yrunleaf = tem.veg.getYRUNLEAF();

  cohort[pichrt].yrvconvrtC = tem.ag.getYRVCONVRTC();
  
  cohort[pichrt].yrvconvrtN = tem.ag.getYRVCONVRTN();

  cohort[pichrt].yrvegC = tem.veg.getYRVEGC();
    
//  cohort[pichrt].yrvegC2N = tem.veg.getYRC2N();

  cohort[pichrt].yrveginnup = tem.veg.getYRINNUP();

  cohort[pichrt].yrveglup = tem.veg.getYRLUP();

  cohort[pichrt].yrvegN = tem.veg.getYRVEGN();

  cohort[pichrt].yrvegnmobil = tem.veg.getYRNMOBIL();

  cohort[pichrt].yrvegnrsorb = tem.veg.getYRNRSORB();

  cohort[pichrt].yrvegnup = tem.veg.getYRNUPTAKE();
  
  cohort[pichrt].yrvegrgrowth = tem.veg.getYRRGROWTH();

  cohort[pichrt].yrvegrmaint = tem.veg.getYRRMAINT();

  cohort[pichrt].yrvegSTON = tem.veg.getYRSTOREN();
  
  cohort[pichrt].yrvegSTRN = tem.veg.getYRSTRUCTN();

  cohort[pichrt].yrvegsup = tem.veg.getYRSUP();

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void MITelmnt44::setCohortTEMState( const MITElmntCohort44& firstchrt,
                                    MITElmntCohort44& targetchrt )
{
  int dlyr;  // NEM
  int dm;
  int i;
  
// targetchrt.srcCohort initialized elsewhere (see xtem606.updateTTEMgridcell())
  
// targetchrt.standage initialized by LCLUC
  
// targetchrt.chrtarea initialized by LCLUC

// targetchrt.prevchrtarea initialized by LCLUC
  
// targetchrt.potveg initialized by LCLUC
  
// targetchrt.currentveg initialized by LCLUC

// targetchrt.subtype initialized by LCLUC

// targetchrt.cmnt initialized by LCLUC

  
  for( i = 0; i < MAXSTATE; ++i )
  {
    targetchrt.y[i] = firstchrt.y[i];
    targetchrt.prevy[i] = firstchrt.prevy[i];
  }

// targetchrt.agcmnt initialized by LCLUC

  targetchrt.aggrowdd = firstchrt.aggrowdd;

  targetchrt.agkd = firstchrt.agkd;

// targetchrt.agprvstate initialized by LCLUC

// targetchrt.agstate initialized by LCLUC

  targetchrt.c2n = firstchrt.c2n;
  
  targetchrt.cneven = firstchrt.cneven;

  targetchrt.convrtflx.carbon = firstchrt.convrtflx.carbon;
  targetchrt.convrtflx.nitrogen = firstchrt.convrtflx.nitrogen;

  targetchrt.cropprveetmx = firstchrt.cropprveetmx;

  targetchrt.cropprvleafmx = firstchrt.cropprvleafmx;

  targetchrt.cropprvpetmx = firstchrt.cropprvpetmx;

  targetchrt.cropResidue.carbon = firstchrt.cropResidue.carbon;
  targetchrt.cropResidue.nitrogen = firstchrt.cropResidue.nitrogen;

  targetchrt.croptopt = firstchrt.croptopt;

  targetchrt.distmnthcnt = firstchrt.distmnthcnt;

// targetchrt.disturbflag initialized by LCLUC

// targetchrt.disturbmonth initialized by LCLUC

  targetchrt.eetmx = firstchrt.eetmx;

// targetchrt.fertflag initialized by LCLUC

  targetchrt.firemnthcnt = firstchrt.firemnthcnt;

  targetchrt.firendep = firstchrt.firendep;

  targetchrt.formPROD10.carbon = firstchrt.formPROD10.carbon;
  targetchrt.formPROD10.nitrogen = firstchrt.formPROD10.nitrogen;

  targetchrt.formPROD100.carbon = firstchrt.formPROD100.carbon;
  targetchrt.formPROD100.nitrogen = firstchrt.formPROD100.nitrogen;

  targetchrt.fprevozone = firstchrt.fprevozone;

// targetchrt.FRI initialized by LCLUC

  for( dm = 0; dm < CYCLE; ++dm )
  {  
    targetchrt.initPROD1[dm].carbon = firstchrt.initPROD1[dm].carbon;
    targetchrt.initPROD1[dm].nitrogen = firstchrt.initPROD1[dm].nitrogen;
  }
  
  for( i = 0; i < 10; ++i )
  {
    targetchrt.initPROD10[i].carbon = firstchrt.initPROD10[i].carbon;
    targetchrt.initPROD10[i].nitrogen = firstchrt.initPROD10[i].nitrogen;
  }
    
  for( i = 0; i < 100; ++i )
  {
    targetchrt.initPROD100[i].carbon = firstchrt.initPROD100[i].carbon;
    targetchrt.initPROD100[i].nitrogen = firstchrt.initPROD100[i].nitrogen;
  }
 
// targetchrt.irrgflag initialized by LCLUC

  targetchrt.kd = firstchrt.kd;

  targetchrt.MDMnpp = firstchrt.MDMnpp;     // MDM

  targetchrt.natprveetmx = firstchrt.natprveetmx;

  targetchrt.natprvleafmx = firstchrt.natprvleafmx;

  targetchrt.natprvpetmx = firstchrt.natprvpetmx;

  targetchrt.natseedC = firstchrt.natseedC;

  targetchrt.natseedSTRN = firstchrt.natseedSTRN;

  targetchrt.natseedSTON = firstchrt.natseedSTON;

  targetchrt.natsoil = firstchrt.natsoil;

  targetchrt.nattopt = firstchrt.nattopt;

  targetchrt.natyreet = firstchrt.natyreet;

  targetchrt.natyrpet = firstchrt.natyrpet;


  // *************************** NEM ***************************
  
  for( dlyr = 0; dlyr < NLVL; ++dlyr )
  {
    targetchrt.NEManh4in[dlyr] = firstchrt.NEManh4in[dlyr]; 
  
    targetchrt.NEMano3in[dlyr] = firstchrt.NEMano3in[dlyr];

    targetchrt.NEMdphumin[dlyr] = firstchrt.NEMdphumin[dlyr]; 

    targetchrt.NEMocin[dlyr] = firstchrt.NEMocin[dlyr];

    targetchrt.NEMrclin[dlyr] = firstchrt.NEMrclin[dlyr]; 
  
    targetchrt.NEMrcrin[dlyr] = firstchrt.NEMrcrin[dlyr];
  
    targetchrt.NEMrcvlin[dlyr] = firstchrt.NEMrcvlin[dlyr];
  }
  
  targetchrt.NEMnsolc = firstchrt.NEMnsolc;

//  targetchrt.NEMtopdens = firstchrt.NEMtopdens;
  
//  targetchrt.NEMtopksat = firstchrt.NEMtopksat;
  
//  targetchrt.NEMtoppor = firstchrt.NEMtoppor;

  // ***********************************************************
  
  
  targetchrt.newleafmx =  firstchrt.newleafmx;

  targetchrt.newtopt = firstchrt.newtopt;

  targetchrt.nretent = firstchrt.nretent;

  targetchrt.nsretent = firstchrt.nsretent;

  targetchrt.nvretent = firstchrt.nvretent;

  targetchrt.petmx = firstchrt.petmx;

  targetchrt.prev2tair = firstchrt.prev2tair;

  targetchrt.prevco2 = firstchrt.prevco2;

  targetchrt.prevCropResidue.carbon = firstchrt.prevCropResidue.carbon;
  targetchrt.prevCropResidue.nitrogen = firstchrt.prevCropResidue.nitrogen;

  targetchrt.prevPROD1.carbon = firstchrt.prevPROD1.carbon;
  targetchrt.prevPROD1.nitrogen = firstchrt.prevPROD1.nitrogen;

  targetchrt.prevPROD10.carbon = firstchrt.prevPROD10.carbon;
  targetchrt.prevPROD10.nitrogen = firstchrt.prevPROD10.nitrogen;

  targetchrt.prevPROD100.carbon = firstchrt.prevPROD100.carbon;
  targetchrt.prevPROD100.nitrogen = firstchrt.prevPROD100.nitrogen;

//  targetchrt.prevspack = firstchrt.prevspack;

  targetchrt.prevtair = firstchrt.prevtair;
  
  targetchrt.prevunrmleaf = firstchrt.prevunrmleaf;

// targetchrt.prod10par initialized by LCLUC
  
// targetchrt.prod100par initialized by LCLUC

  targetchrt.productYear = firstchrt.productYear;
  
  targetchrt.prvcropnpp = firstchrt.prvcropnpp;

  targetchrt.prveetmx = firstchrt.prveetmx;

  targetchrt.prvleafmx = firstchrt.prvleafmx;

  targetchrt.prvpetmx = firstchrt.prvpetmx;

  targetchrt.qc = firstchrt.qc;

// targetchrt.sconvert initialized by LCLUC

  targetchrt.sconvrtflx.carbon = firstchrt.sconvrtflx.carbon;
  targetchrt.sconvrtflx.nitrogen = firstchrt.sconvrtflx.nitrogen;

  targetchrt.slash.carbon = firstchrt.slash.carbon;
  targetchrt.slash.nitrogen = firstchrt.slash.nitrogen;

// targetchrt.slashpar initialized by LCLUC

// targetchrt.tillflag initialized by LCLUC

  targetchrt.topt = firstchrt.topt;

  targetchrt.tqc = firstchrt.tqc;

// targetchrt.vconvert initialized by LCLUC

  targetchrt.vconvrtflx.carbon = firstchrt.vconvrtflx.carbon;
  targetchrt.vconvrtflx.nitrogen = firstchrt.vconvrtflx.nitrogen;

// targetchrt.vrespar initialized by LCLUC

  targetchrt.yragstubC = firstchrt.yragstubC;
  
  targetchrt.yragstubN = firstchrt.yragstubN;

  targetchrt.yrcflux = firstchrt.yrcflux;
  
  targetchrt.yrCH4csmp = firstchrt.yrCH4csmp;
  
  targetchrt.yrCH4ems = firstchrt.yrCH4ems;
  
  targetchrt.yrCH4flx = firstchrt.yrCH4flx;
  
  targetchrt.yrCO2dnflx = firstchrt.yrCO2dnflx;
  
  targetchrt.yrCO2nflx = firstchrt.yrCO2nflx;
  
  targetchrt.yrconvrtC = firstchrt.yrconvrtC;
  
  targetchrt.yrconvrtN = firstchrt.yrconvrtN;
  
  targetchrt.yrdecayPROD1C = firstchrt.yrdecayPROD1C;
  
  targetchrt.yrdecayPROD10C = firstchrt.yrdecayPROD10C;
  
  targetchrt.yrdecayPROD100C = firstchrt.yrdecayPROD100C;

  targetchrt.yrdecayPROD1N = firstchrt.yrdecayPROD1N;
  
  targetchrt.yrdecayPROD10N = firstchrt.yrdecayPROD10N;
  
  targetchrt.yrdecayPROD100N = firstchrt.yrdecayPROD100N;
  
  targetchrt.yrdecayTOTPRODC = firstchrt.yrdecayTOTPRODC;

  targetchrt.yrdecayTOTPRODN = firstchrt.yrdecayTOTPRODN;

  targetchrt.yreet = firstchrt.yreet;
  
  targetchrt.yrfertn = firstchrt.yrfertn;
  
  targetchrt.yrfluxResidueC = firstchrt.yrfluxResidueC;
  
  targetchrt.yrfluxResidueN = firstchrt.yrfluxResidueN;

  targetchrt.yrformPROD1C = firstchrt.yrformPROD1C;
  
  targetchrt.yrformPROD10C = firstchrt.yrformPROD10C;
  
  targetchrt.yrformPROD100C = firstchrt.yrformPROD100C;

  targetchrt.yrformPROD1N = firstchrt.yrformPROD1N;
  
  targetchrt.yrformPROD10N = firstchrt.yrformPROD10N;
  
  targetchrt.yrformPROD100N = firstchrt.yrformPROD100N;
  
  targetchrt.yrformResidueC = firstchrt.yrformResidueC;
  
  targetchrt.yrformResidueN = firstchrt.yrformResidueN;

  targetchrt.yrformTOTPRODC = firstchrt.yrformTOTPRODC;

  targetchrt.yrformTOTPRODN = firstchrt.yrformTOTPRODN;

  targetchrt.yrfpc = firstchrt.yrfpc;

  targetchrt.yrgpp = firstchrt.yrgpp;
  
  targetchrt.yrgpr = firstchrt.yrgpr;

//  targetchrt.yrh2oyld = firstchrt.yrh2oyld;
  
  targetchrt.yrimmob = firstchrt.yrimmob;
  
//  targetchrt.yrineet = firstchrt.yrineet;
  
  targetchrt.yringpp = firstchrt.yringpp;
  
  targetchrt.yrinnpp = firstchrt.yrinnpp;
   
//  targetchrt.yrirrig = firstchrt.yrirrig;

  targetchrt.yrlai = firstchrt.yrlai;

  targetchrt.yrleaf = firstchrt.yrleaf;

  targetchrt.yrltrfalC = firstchrt.yrltrfalC;

  targetchrt.yrltrfalN = firstchrt.yrltrfalN;  	
  
  targetchrt.yrN2flx = firstchrt.yrN2flx;
  
  targetchrt.yrN2Odnflx = firstchrt.yrN2Odnflx;
  
  targetchrt.yrN2Oflx = firstchrt.yrN2Oflx;
  
  targetchrt.yrN2Onflx = firstchrt.yrN2Onflx;
  
  targetchrt.yrnce = firstchrt.yrnce;
  
  targetchrt.yrnep = firstchrt.yrnep;
  
  targetchrt.yrninput = firstchrt.yrninput;
  
  targetchrt.yrnlost = firstchrt.yrnlost;
  
  targetchrt.yrnmin = firstchrt.yrnmin;
   
  targetchrt.yrnpp = firstchrt.yrnpp;
  
  targetchrt.yrnrent = firstchrt.yrnrent;
    
  targetchrt.yrnsrent = firstchrt.yrnsrent;

  targetchrt.yrnvrent = firstchrt.yrnvrent;
  
  targetchrt.yrpet = firstchrt.yrpet;
    
//  targetchrt.yrrgrndh2o = firstchrt.yrrgrndH2O; 

  targetchrt.yrrh = firstchrt.yrrh;
    
//  targetchrt.yrrain = firstchrt.yrrain;
  
//  targetchrt.yrrperc = firstchrt.yrrperc;
  
//  targetchrt.yrrrun = firstchrt.yrrrun;

  targetchrt.yrsconvrtC = firstchrt.yrsconvrtC;

  targetchrt.yrsconvrtN = firstchrt.yrsconvrtN;
  
//  targetchrt.yrsgrndh2o = firstchrt.yrsgrndH2O;

  targetchrt.yrslashC = firstchrt.yrslashC;

  targetchrt.yrslashN = firstchrt.yrslashN;
  
//  targetchrt.yrsnowfall = firstchrt.yrsnowfall;

//  targetchrt.yrsnowinf = firstchrt.yrsnowinf;

//  targetchrt.yrsnowpack = firstchrt.yrsnowpack;

//  targetchrt.yrsoilavlH2O = firstchrt.yrsoilavlH2O;

  targetchrt.yrsoilavlN = firstchrt.yrsoilavlN;
  
//  targetchrt.yrsoilC2N = firstchrt.yrsoilC2N;

//  targetchrt.yrsoilmoist = firstchrt.yrsoilmoist;

  targetchrt.yrsoilorgC = firstchrt.yrsoilorgC;
  
  targetchrt.yrsoilorgN = firstchrt.yrsoilorgN;
    
//  targetchrt.yrsoilpctp = firstchrt.yrsoilpctp;

//  targetchrt.yrsoilvsm = firstchrt.yrsoilvsm;

//  targetchrt.yrsperc = firstchrt.yrsperc;
  
//  targetchrt.yrsrun = firstchrt.yrsrun;

  targetchrt.yrtotalC = firstchrt.yrtotalC;
  
  targetchrt.yrunleaf = firstchrt.yrunleaf;

  targetchrt.yrvconvrtC = firstchrt.yrvconvrtC;
  
  targetchrt.yrvconvrtN = firstchrt.yrvconvrtN;

  targetchrt.yrvegC = firstchrt.yrvegC;

//  targetchrt.yrvegC2N = firstchrt.yrvegC2N;

  targetchrt.yrveginnup = firstchrt.yrveginnup;
    
  targetchrt.yrveglup = firstchrt.yrveglup;

  targetchrt.yrvegN = firstchrt.yrvegN;
  
  targetchrt.yrvegnmobil = firstchrt.yrvegnmobil;

  targetchrt.yrvegnrsorb = firstchrt.yrvegnrsorb;

  targetchrt.yrvegnup = firstchrt.yrvegnup;

  targetchrt.yrvegrgrowth = firstchrt.yrvegrgrowth;

  targetchrt.yrvegrmaint = firstchrt.yrvegrmaint;

  targetchrt.yrvegSTON = firstchrt.yrvegSTON;
  
  targetchrt.yrvegSTRN = firstchrt.yrvegSTRN;

  targetchrt.yrvegsup = firstchrt.yrvegsup;

};

/* *************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

int MITelmnt44::setGIStopography( ofstream& rflog1,
                                  int& ftlerr,
                                  FILE* fstxt,
                                  FILE* fslayer,
                                  FILE* felev )
{

  int dlyr;   // NEM
  int gisend;

  MITSoildata44 fao;

  MITSoilLayerdata44 clmsoil;

  Elevdata43 elv;

  gisend = fao.getdel( fstxt );
 
  if( -1 == gisend )
  {
    rflog1 << "Ran out of Soil texture data" << endl << endl;
 
    exit( -1 );
  }
  
  ftlerr = coregerr( rflog1,
                     "Climate",
                     col,
                     row,
                     "TEXTURE",
                     fao.col,
                     fao.row );

  tem.soil.setPCTSAND( fao.pctsand );

  tem.soil.setPCTSILT( fao.pctsilt );

  tem.soil.setPCTCLAY( fao.pctclay );

  tem.soil.setWSOIL( fao.wsoil );

  tem.soil.setPH( fao.pH );
  

  if( 1 == tem.ch4flag || 1 == tem.n2oflag )
  {
    for ( dlyr = 0; dlyr < CLMNLAYERS; ++dlyr )
    {
      gisend = clmsoil.getdel( fslayer );

      if( -1 == gisend )
      {
        rflog1 << "Ran out of Soil layer data" << endl << endl;
 
        exit( -1 );
      }
    
      tem.soil.layerThick[dlyr] = clmsoil.thickness;
      
      tem.soil.porosity[dlyr] = clmsoil.porosity;
      
      tem.soil.density[dlyr] = (1.0 - tem.soil.porosity[dlyr]) * 2.7000; 
    
      tem.soil.Ksat[dlyr] = clmsoil.Ksat;
    }


    // Use the "cohort" version if soil porosity, soil density 
    //   and saturated hydraulic conductivity vary with
    //   vegetation types.  Otherwise, just read in these data
    //   directly into the associated NEM variables
    
    // Determine effective soil porosity of the top 30 cm
  
//    cohort[0].NEMtoppor = (1.75/29.0) * tem.soil.porosity[0]
//                          + (2.75/29.0) * tem.soil.porosity[1]
//                          + (4.5/29.0) * tem.soil.porosity[2]
//                          + (7.5/29.0) * tem.soil.porosity[3]
//                          + (12.5/29.0) * tem.soil.porosity[4];

    // Determine effective soil density of the top 30 cm
  
//    cohort[0].NEMtopdens = (1.75/29.0) * tem.soil.density[0]
//                           + (2.75/29.0) * tem.soil.density[1]
//                           + (4.5/29.0) * tem.soil.density[2]
//                           + (7.5/29.0) * tem.soil.density[3]
//                           + (12.5/29.0) * tem.soil.density[4];

    // Determine effective saturated soil hydraulic conductivity 
    //   of the top 30 cm
  
//    cohort[0].NEMtopksat = (1.75/29.0) * tem.soil.Ksat[0]
//                           + (2.75/29.0) * tem.soil.Ksat[1]
//                           + (4.5/29.0) * tem.soil.Ksat[2]
//                           + (7.5/29.0) * tem.soil.Ksat[3]
//                           + (12.5/29.0) * tem.soil.Ksat[4];

    // Determine effective soil porosity of the top 30 cm

    tem.microbe.nem.setTOPPOR( ((1.75/29.0) * tem.soil.porosity[0]
                               + (2.75/29.0) * tem.soil.porosity[1]
                               + (4.5/29.0) * tem.soil.porosity[2]
                               + (7.5/29.0) * tem.soil.porosity[3]
                               + (12.5/29.0) * tem.soil.porosity[4]) );

    // Determine effective soil density of the top 30 cm

    tem.microbe.nem.setTOPDENS( ((1.75/29.0) * tem.soil.density[0]
                                + (2.75/29.0) * tem.soil.density[1]
                                + (4.5/29.0) * tem.soil.density[2]
                                + (7.5/29.0) * tem.soil.density[3]
                                + (12.5/29.0) * tem.soil.density[4]) );

    // Determine effective saturated soil hydraulic conductivity 
    //   of the top 30 cm

    tem.microbe.nem.setTOPKSAT( ((1.75/29.0) * tem.soil.Ksat[0]
                                + (2.75/29.0) * tem.soil.Ksat[1]
                                + (4.5/29.0) * tem.soil.Ksat[2]
                                + (7.5/29.0) * tem.soil.Ksat[3]
                                + (12.5/29.0) * tem.soil.Ksat[4]) );
 
  }
  else 
  {
//    cohort[0].NEMtoppor = ZERO;

//    cohort[0].NEMtopdens = ZERO;

//    cohort[0].NEMtopksat = ZERO;

    tem.microbe.nem.setTOPPOR( ZERO ); 

    tem.microbe.nem.setTOPDENS( ZERO ); 

    tem.microbe.nem.setTOPKSAT( ZERO ); 
  }
  
  gisend = elv.getdel( felev );
  
  if( gisend == -1 )
  {
    rflog1 << "Ran out of Elevation data" << endl << endl;
 
    exit( -1 );
  }
  
  ftlerr = coregerr( rflog1,
                     "Climate",
                     col,
                     row,
                     "ELEV",
                     elv.col,
                     elv.row );

  tem.elev = elv.elev;

  return gisend;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void MITelmnt44::setTEMequilState( ofstream& rflog1,
                                   const int& equil,
                                   const int& totsptime,
                                   const int& pichrt )
{
  const int dyr = 0;
  
  int dm;

//  double totalSoilC;
//  double outn2oflx;

  // Set all TEM state-related variables in cohort to MISSING
  //   (i.e. start with "clean slate" in cohort)
   
  initializeCohortTEMState( pichrt );

  // Assign cohort data to TEM (i.e. transfer information from 
  //   the land cover/land use module to TEM and start with 
  //   "clean slate" for TEM cohort data)
   
  getTEMCohortState( pichrt );
  
  
  cohort[pichrt].qc = ACCEPT;

  cohort[pichrt].tqc = ACCEPT;

  tem.qualcon[dyr] = 0;

  tem.totyr = 0;

//  tem.ag.fert1950flag = 0;

  tem.microbe.nem.startflag = 0;
  
  tem.atms.setMXTAIR( mxtair );
  
  tem.atms.setYRPREC( yrprec );


  for( dm = 0; dm < CYCLE; ++dm )
  {                                
    // Pass climate data for particular month to TEM

    tem.atms.setNIRR( climate[mitclm.I_NIRR][dm] );

    tem.atms.setPAR( climate[mitclm.I_PAR][dm] );

    tem.atms.setTAIR( climate[mitclm.I_TAIR][dm] );

    tem.atms.setPREC( climate[mitclm.I_PREC][dm] );

    tem.atms.setPET( initPET[pichrt][dm] );

    tem.soil.setEET( initAET[pichrt][dm] );

    tem.soil.setMOIST( initSH2O[pichrt][dm] );

    tem.soil.setSNOWPACK( initSNOWPACK[pichrt][dm] );

    tem.soil.setSURFRUN( initSURFRUN[pichrt][dm] );

    tem.soil.setDRAINAGE( initDRAINAGE[pichrt][dm] );

    tem.atms.setCO2( climate[mitclm.I_CO2][dm] );

    tem.atms.setAOT40( climate[mitclm.I_AOT40][dm] );


    // Check TEM input for valid data 

    cohort[pichrt].qc = temgisqc( cohort[pichrt].chrtarea,
                                  tem.soil.getPCTSILT(),
                                  tem.soil.getPCTCLAY(),
                                  tem.veg.cmnt,
                                  tem.elev,
                                  tem.atms.getNIRR(),
                                  tem.atms.getPAR(),
                                  tem.atms.getTAIR(),
                                  tem.atms.getMXTAIR(),
                                  tem.atms.getYRPREC(),
                                  tem.atms.getPREC(),
                                  tem.atms.getPET(),
                                  tem.soil.getEET(),
                                  tem.soil.getMOIST(),
                                  tem.soil.getSNOWPACK(),
                                  tem.soil.getSURFRUN(),
                                  tem.soil.getDRAINAGE(),
                                  tem.atms.getCO2(),
                                  tem.atms.getAOT40() );


    if( cohort[pichrt].qc != ACCEPT ) 
    { 
      rflog1 << "temgisqc = " << cohort[pichrt].qc;
      rflog1 << " for cohort " << (pichrt+1);
      rflog1 << " during month " << (dm+1) << endl;
      break; 
    }
     
    // Determine initial values for tem.atms.prvpetmx, 
    //   tem.atms.prveetmx and and tem.veg.topt based on
    //   long-term mean climate

    tem.setEquilEvap( tem.atms.getNIRR(), 
                      tem.atms.getTAIR(), 
                      dm );
  }

//  exit( -1 );

   // Determine soil properties of element based on 
  //   soil texture

  tem.soil.xtext( tem.veg.cmnt,
                  tem.soil.getPCTSILT(),
                  tem.soil.getPCTCLAY() );


  // Assume potential vegetation when determining 
  //   equilibrium conditions

  tem.ag.state = 0;
  tem.ag.prvstate = 0;

  tem.ag.tillflag = 0;
  tem.ag.fertflag = 0;
  tem.ag.irrgflag = 0;

  tem.disturbflag = 0;
  tem.distmnthcnt = 0;

  // Initialize all fluxes and disturbance-related variables to zero

  tem.ag.setNATSEEDC( ZERO );
  tem.ag.setNATSEEDSTRN( ZERO );
  tem.ag.setNATSEEDSTON( ZERO );
  tem.ag.setCROPPRVLEAFMX( ZERO );
  tem.ag.setCROPTOPT( ZERO );
  tem.ag.setCROPPRVPETMX( ZERO );
  tem.ag.setCROPPRVEETMX( ZERO );
  tem.ag.setGROWDD( ZERO );
  tem.ag.setPRVCROPNPP( ZERO );

  tem.ag.resetMonthlyDisturbFluxes();
  tem.resetMonthlyELMNTFluxes();
  tem.resetYrFluxes();

   
  // Check TEM parameters for specific vegetation types

  if( ACCEPT == cohort[pichrt].qc ) 
  { 
    cohort[pichrt].qc = tem.ecdqc( tem.veg.cmnt ); 

    if( cohort[pichrt].qc != ACCEPT )
    {
      // Note: If a TEM parameter is invalid, 
      //   cohort[pichrt].qc will have a value greater than 
      //   100
      	
      rflog1 << "temecdqc = " << cohort[pichrt].qc << endl;
    }
  }

  if( cohort[pichrt].qc != ACCEPT )
  {
    // If environmental conditions are too extreme for the 
    //   existence of vegetation (e.g., no precipitation or 
    //   constant freezing air temperatures), assign zero to 
    //   all TEM variables if the plant community is anything
    //   besides ice and open water; and all TEM parameters 
    //   are valid (i.e. cohort[pichrt].qc < 100 )

    // Set tqc flag to assign zero to all TEM variables 
    //   during simulation
    	
    cohort[pichrt].tqc = TQCZEROFLAG; 


    // Set missing values to telmnt[0].output

    setTEMmiss( dyr,
                equil,
                totsptime,
                pichrt  );
  }
  else // "cohort[pichrt].qc == ACCEPT"
  {

/* *************************************************************
                   Start Equilibrium Conditions
************************************************************* */

    // Determine soil properties of element based on 
    //   soil texture
    
//    tem.soil.xtext( tem.veg.cmnt, 
//                    tem.soil.getPCTSILT(), 
//                    tem.soil.getPCTCLAY() );

    
    // Initialize tem.atms.prevco2

    tem.atms.setPREVCO2( tem.atms.getCO2LEVEL() );


    // Initialize TEM parameters based on element's 
    //   (i.e. grid cell) vegetation type, soil texture
    //   and atmospheric CO2 concentration
      
    tem.setELMNTecd( tem.veg.cmnt, tem.soil.getPSIPLUSC() );

    tem.setEquilC2N( tem.veg.cmnt, 
                     tem.atms.getPREVCO2() );

          
    // "While" loop to allow adaptive integrator tolerance 
    //   (i.e. tem.tol) to be reduced if chaotic behavior 
    //   occurs


    // Try up to "tem.maxnrun" times to equilibrate TEM.  If
    //   TEM does not equilibrate within "tem.runsize" 
    //   iterations, decrease tem.tol by an order of magnitude
    //   and try again
     
    tem.nattempt = 0;

    tem.tol = tem.inittol;

    tem.baseline = tem.initbase;

    tem.initFlag = 0;
        
    while( tem.nattempt < tem.maxnrun 
           && 0 == tem.initFlag )
    {               
      tem.nattempt = equilibrateTEM( pichrt, tem.tol );
	                               
      if( tem.nattempt < tem.maxnrun 
          && 0 == tem.initFlag ) 
      { 
      	tem.tol /= 10.0; 
      }
    }

    // Update summary variables for initial agricultural 
    //   state of cohort at end of equilibrium portion 
    //   of the TEM simulation

    tem.ag.setNATSEEDC( ZERO );

    tem.ag.setNATSEEDSTRN( ZERO );

    tem.ag.setNATSEEDSTON( ZERO );

    tem.ag.setCROPPRVLEAFMX( 1.0 );

    tem.ag.setCROPTOPT( tem.veg.getTOPT() );

    tem.ag.setCROPPRVPETMX( tem.atms.getPRVPETMX() );

    tem.ag.setCROPPRVEETMX( tem.soil.getPRVEETMX() );

    tem.ag.setPRVCROPNPP( ZERO );


    // Save quality control information about the simulation 
    //   conditions when the equilibrium portion ended
    //   (i.e. did the carbon and nitrogen fluxes really come 
    //         to equilibrium or was the run terminated after
    //         running chaotically for a specified maximum
    //         number of years?) 

    tem.qualcon[dyr] += (tem.nattempt + 1);
    
    // If simulation is part of a transient simulation, reset
    //   tem.totyr to represent an actual year rather than 
    //   the number of iterations required to reach equilibrum
    
    if( 0 == equil )
    {
      tem.totyr = tem.startyr - totsptime - 1;
 
      ttotyr[dyr] = tem.totyr;
      
      cohort[pichrt].tqc = transqc( tem.maxyears, 
	                                  tem.totyr, 
	                                  output[tem.I_VEGC][pichrt] );
    }
    else { ttotyr[dyr] = tem.totyr; }
  } // End of "cohort.qc == ACCEPT"

  // Save TEM state of cohort to telmnt[0].cohort

  saveTEMCohortState( pichrt );

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void MITelmnt44::setTEMmiss( const int& pdyr,
                             const int& equil,
                             const int& totsptime,
                             const int& pichrt )
{
  int dm;
  int i;
  
  if( 0 == equil )
  {
    ttotyr[pdyr] = tem.startyr 
                   - totsptime - 1 
                   + (pdyr * tem.diffyr);
  }
  else
  {
    ttotyr[pdyr] = -999;
  }

  tem.totyr = ttotyr[pdyr];

  if( TQCZEROFLAG == cohort[pichrt].tqc )
  {
    if( 1 == equil ) { ttotyr[pdyr] = 1; }

    // Assign zero to all TEM state variables
      
    for( i = 0; i < MAXSTATE; ++i )
    {
      tem.setY( ZERO, i );
      tem.setPREVY(ZERO, i );
    }

    for( i = MAXSTATE; i < NUMEQ; ++i )
    {
      tem.setY( ZERO, i );
    }
      
    // Assign zero to all TEM output variables
      
    for(i = 0; i < NUMTEM; ++i )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        output[i][pichrt][dm] = ZERO;
      }
    }
  }
  else
  {
    // Assign missing values to grid cells that are covered by ice or open
    // water, or where TEM did not converge on a solution

    for( i = 0; i < MAXSTATE; ++i )
    {
      tem.setY( MISSING, i );
      tem.setPREVY( MISSING, i );
    }

    for( i = MAXSTATE; i < NUMEQ; ++i )
    {
      tem.setY( MISSING, i );
    }
      
    for(i = 0; i < NUMTEM; ++i )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        output[i][pichrt][dm] = MISSING;
      }
    }
  }

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

int MITelmnt44::temgisqc( const long& subarea,
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
                          const double& aot40 )


{
  int qc;

  qc = ACCEPT;

  if( subarea < 1 ) { return qc = 1; }

  if( pctsilt < ZERO ) { return qc = 2; }

  if( pctclay < ZERO ) { return qc = 3; }

  if( cmnt < 2 || cmnt > NUMVEG ) { return qc = 4; }

  if( elev <= -999.0 ) { return qc = 5;}

  if( nirr <= -1.0 ) { return qc = 6; }

  if( par <= -1.0 ) { return qc = 7; }

  if( tair <= -99.0 ) { return qc = 8; }

  if( mxtair < -1.0 ) { return qc = 9; }

  if( yrprec <= ZERO ) { return qc = 10; }

  if( prec <= -1.0 ) { return qc = 11; }

  if( pet <= -1.0 ) { return qc = 12; }

  if( eet <= -1.0 ) { return qc = 13; }

  if( sh2o <= -1.0 ) { return qc = 14; }

  if( spack <= -1.0 ) { return qc = 15; }

  if( srfrun <= -1.0 ) { return qc = 16; }

  if( drain <= -1.0 ) { return qc = 17; }

  if( co2 <= -1.0 ) { return qc = 18; }

  if( aot40 <= -1.0 ) { return qc = 19; }

  return qc;

};

/* *************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

int MITelmnt44::temgisqc( const long& subarea,
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
                          const double& aot40 )


{
  int qc;

  qc = ACCEPT;

  if( subarea < 1 ) { return qc = 1; }

  if( pctsilt < ZERO ) { return qc = 2; }

  if( pctclay < ZERO ) { return qc = 3; }

  if( cmnt < 2 || cmnt > NUMVEG ) { return qc = 4; }

  if( elev <= -999.0 ) { return qc = 5;}

  if( nirr <= -1.0 ) { return qc = 6; }

  if( par <= -1.0 ) { return qc = 7; }

  if( tair <= -99.0 ) { return qc = 8; }

  if( prec <= -1.0 ) { return qc = 11; }

  if( pet <= -1.0 ) { return qc = 12; }

  if( eet <= -1.0 ) { return qc = 13; }

  if( sh2o <= -1.0 ) { return qc = 14; }

  if( spack <= -1.0 ) { return qc = 15; }

  if( srfrun <= -1.0 ) { return qc = 16; }

  if( drain <= -1.0 ) { return qc = 17; }

  if( co2 <= -1.0 ) { return qc = 18; }

  if( aot40 <= -1.0 ) { return qc = 19; }

  return qc;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void MITelmnt44::temwritepred( ofstream fout[NUMTEM],
                               const vector<string>& predname,
                               const int& pdyr,
                               const int& pichrt,
                               const int& ntempred )
{
  // Covert cal/cm2/day to W/m2 (4.186 Joules / calorie)
  const double  cal2Watts = 0.4845;

  // Units conversion from grams to milligrams
  const double GRAMS2MG = 1000.0;
   
  // Units conversion from proportion to percent
  const double PROP2PCT = 100.0;
  
  int i;
  int dm;
  Temdata44 tempred;


  for( i = 0; i < ntempred; ++i )
  {
    // ************** Carbon stocks in ecosystems  *************


    if( predname.at( i ) == tem.predstr.at( tem.I_VEGC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_VEGC][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_SOLC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_SOLC][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_TOTEC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_TOTEC][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_TOTC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_TOTC][pichrt][dm];
      }
    }


    // *************** Nitrogen stocks in ecosystems ***********

    else if( predname.at( i ) == tem.predstr.at( tem.I_STRN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_STRN][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_STON ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_STON][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_SOLN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_SOLN][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AVLN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AVLN][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_VEGN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_VEGN][pichrt][dm];
      }
    }


    // *****************Water stocks in ecosystems *************

    else if( predname.at( i ) == tem.predstr.at( tem.I_AVLW ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AVLW][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_SM ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_SM][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_VSM ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_VSM][pichrt][dm] * PROP2PCT; 
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PCTP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PCTP][pichrt][dm] * PROP2PCT;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_SNWPCK ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_SNWPCK][pichrt][dm];
      }
    }


   // ******************** Phenology ***************************


    else if( predname.at( i ) == tem.predstr.at( tem.I_UNRMLF ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_UNRMLF][pichrt][dm] * PROP2PCT;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_LEAF ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_LEAF][pichrt][dm] * PROP2PCT;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_LAI ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_LAI][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_FPC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_FPC][pichrt][dm] * PROP2PCT;
      }
    }


    // *************** Carbon fluxes in ecosystems *************


    else if( predname.at( i ) == tem.predstr.at( tem.I_INGPP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_INGPP][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_GPP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_GPP][pichrt][dm];
      }
    }

    // *********************** Ozone Effects *******************

    else if( predname.at( i ) == tem.predstr.at( tem.I_FOZONE ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_FOZONE][pichrt][dm] * PROP2PCT;  

      }
    }
    else if( predname.at( i ) == tem.predstr.at( tem.I_FINDOZONE ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_FINDOZONE][pichrt][dm] * PROP2PCT;  

      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_INNPP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_INNPP][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NPP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NPP][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_GPR ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_GPR][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_RVMNT ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_RVMNT][pichrt][dm];  
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_RVGRW ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_RVGRW][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_LTRC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_LTRC][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_RH ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_RH][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NEP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NEP][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NCE ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NCE][pichrt][dm];
      }
    }

    // ****************** NEM and MDM Methane Fluxes ***********

    else if( predname.at( i ) == tem.predstr.at( tem.I_CH4EMS ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CH4EMS][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_CH4CSMP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CH4CSMP][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_CH4FLX ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CH4FLX][pichrt][dm] * GRAMS2MG;
      }
    }


    // ********************* NEM CO2 Fluxes ********************

    else if( predname.at( i ) == tem.predstr.at( tem.I_CO2NFLX ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CO2NFLX][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_CO2DNFLX ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CO2DNFLX][pichrt][dm] * GRAMS2MG;
      }
    }

    // ************** Nitrogen fluxes in ecosystems ************


    else if( predname.at( i ) == tem.predstr.at( tem.I_NINP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NINP][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGFRTN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGFRTN][pichrt][dm] * GRAMS2MG; 
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_INNUP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_INNUP][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_VNUP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_VNUP][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_VSUP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_VSUP][pichrt][dm] * GRAMS2MG;  
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_VLUP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_VLUP][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_VNMBL ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_VNMBL][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_VNRSRB ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_VNRSRB][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_LTRN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_LTRN][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_MNUP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_MNUP][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NMIN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NMIN][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NLST ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NLST][pichrt][dm] * GRAMS2MG;
      }
    }


    // ***************** NEM NO, N2O and N2 Fluxes *************

    else if( predname.at( i ) == tem.predstr.at( tem.I_NOFLX ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NOFLX][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_N2OFLX ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_N2OFLX][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_N2ONFLX ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_N2ONFLX][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_N2ODNFLX ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_N2ODNFLX][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_N2FLX ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_N2FLX][pichrt][dm] * GRAMS2MG;
      }
    }

    // *****************Water fluxes in ecosystems *************

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGIRRIG ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGIRRIG][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_INEET ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_INEET][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_EET ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_EET][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PET ) )
    {
      for ( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PET][pichrt][dm];
      }
    }


// ************** Carbon stocks in products ********************


    else if( predname.at( i ) == tem.predstr.at( tem.I_AGPRDC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGPRDC][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PROD10C ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PROD10C][pichrt][dm];      
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PROD100C ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PROD100C][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_TOTPRDC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_TOTPRDC][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_RESIDC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_RESIDC][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGSTUBC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGSTUBC][pichrt][dm];
      }
    }

    // ************** Nitrogen stocks in products **************


    else if( predname.at( i ) == tem.predstr.at( tem.I_AGPRDN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGPRDN][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PROD10N ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PROD10N][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PROD100N ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PROD100N][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_TOTPRDN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_TOTPRDN][pichrt][dm] * GRAMS2MG; 
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_RESIDN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_RESIDN][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGSTUBN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGSTUBN][pichrt][dm];  
      }
    }


    // *** Carbon fluxes during agricultural conversion ********


    else if( predname.at( i ) == tem.predstr.at( tem.I_CNVRTC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CNVRTC][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_VCNVRTC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_VCNVRTC][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_SCNVRTC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_SCNVRTC][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_SLASHC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_SLASHC][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_CFLX ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CFLX][pichrt][dm];
      }
    }


    // *** Nitrogen fluxes during agricultural conversion ******


    else if( predname.at( i ) == tem.predstr.at( tem.I_CNVRTN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CNVRTN][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_VCNVRTN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_VCNVRTN][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_SCNVRTN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_SCNVRTN][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_SLASHN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
      	 tempred.mon[dm] = output[tem.I_SLASHN][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NRETNT ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NRETNT][pichrt][dm]; 
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NVRTNT ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NVRTNT][pichrt][dm]; 
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NSRTNT ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NSRTNT][pichrt][dm];
      }
    }


    // ************** Carbon fluxes to/from products ***********

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGFPRDC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGFPRDC][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PRDF10C ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PRDF10C][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PRDF100C ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PRDF100C][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_TOTFPRDC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_TOTFPRDC][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_FRESIDC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_FRESIDC][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGPRDFC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGPRDFC][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PRD10FC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PRD10FC][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PRD100FC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PRD100FC][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_TOTPRDFC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_TOTPRDFC][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_RESIDFC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_RESIDFC][pichrt][dm];
      }
    }

    // ************** Nitrogen fluxes to/from products *********


    else if( predname.at( i ) == tem.predstr.at( tem.I_AGFPRDN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGFPRDN][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PRDF10N ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PRDF10N][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PRDF100N ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PRDF100N][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_TOTFPRDN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_TOTFPRDN][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_FRESIDN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_FRESIDN][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGPRDFN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGPRDFN][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PRD10FN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PRD10FN][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PRD100FN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PRD100FN][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_TOTPRDFN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_TOTPRDFN][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_RESIDFN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_RESIDFN][pichrt][dm];
      }
    }

    // ************** Carbon stocks in crops   *****************


    else if( predname.at( i ) == tem.predstr.at( tem.I_CROPC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CROPC][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATVEGC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATVEGC][pichrt][dm];
      }
    }


    // ************** Nitrogen stocks in crops *****************


    else if( predname.at( i ) == tem.predstr.at( tem.I_CROPN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CROPN][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATVEGN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATVEGN][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_CSTRN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CSTRN][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATSTRN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATSTRN][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_CSTON ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CSTON][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATSTON ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATSTON][pichrt][dm];
      }
    }

    // ******************** Crop Phenology *********************


    else if( predname.at( i ) == tem.predstr.at( tem.I_CROPULF ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CROPULF][pichrt][dm] * PROP2PCT;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATULF ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATULF][pichrt][dm] * PROP2PCT;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_CROPLEAF ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CROPLEAF][pichrt][dm] * PROP2PCT;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATLEAF ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATLEAF][pichrt][dm] * PROP2PCT;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_CROPLAI ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CROPLAI][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATLAI ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATLAI][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_CROPFPC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CROPFPC][pichrt][dm] * PROP2PCT;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATFPC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATFPC][pichrt][dm] * PROP2PCT;
      }
    }

    // ************** Carbon fluxes in croplands ***************

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGINGPP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGINGPP][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATINGPP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATINGPP][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGGPP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGGPP][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATGPP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATGPP][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGINNPP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGINNPP][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATINNPP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATINNPP][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGNPP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGNPP][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATNPP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATNPP][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGGPR ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGGPR][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATGPR ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATGPR][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGRVMNT ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGRVMNT][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATRVMNT ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATRVMNT][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGRVGRW ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGRVGRW][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATRVGRW ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATRVGRW][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGLTRC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGLTRC][pichrt][dm]; 
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATLTRC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATLTRC][pichrt][dm];
      }
    }

    // ************** Nitrogen fluxes in croplands *************


    else if( predname.at( i ) == tem.predstr.at( tem.I_AGINNUP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGINNUP][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATINNUP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATINNUP][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGVNUP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGVNUP][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATVNUP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATVNUP][pichrt][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGVSUP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGVSUP][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATVSUP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATVSUP][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGVLUP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGVLUP][pichrt][dm];  
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATVLUP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATVLUP][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGVNMBL ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGVNMBL][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATVNMBL ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATVNMBL][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGVNRSRB ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGVNRSRB][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NVNRSRB ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NVNRSRB][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGLTRN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGLTRN][pichrt][dm] * GRAMS2MG; 
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATLTRN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATLTRN][pichrt][dm] * GRAMS2MG;
      }
    }

    
    // ******* IGSM Climate Output in TEM output format ********
    
    else if( predname.at( i ) == tem.predstr.at( tem.I_NIRR ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NIRR][pichrt][dm] * cal2Watts; 
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PAR ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PAR][pichrt][dm] * cal2Watts; 
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_TAIR ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_TAIR][pichrt][dm]; 
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PREC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PREC][pichrt][dm]; 
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_SRFRUN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_SRFRUN][pichrt][dm]; 
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_DRAIN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_DRAIN][pichrt][dm];   
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_CO2 ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CO2][pichrt][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AOT40 ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AOT40][pichrt][dm];  
      }
    }

    else
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = MISSING;
      }
    }

    // Write output data to files

    if( predname.at( i ) == tem.predstr.at( tem.I_VSM ) 
        || predname.at( i ) == tem.predstr.at( tem.I_PCTP ) 
        || predname.at( i ) == tem.predstr.at( tem.I_LEAF ) )
    {
      tempred.poutdel( fout[i],
                       col,
                       row,
                       predname.at( i ),
                       (pichrt+1),
                       cohort[pichrt].standage,
                       tem.veg.getPOTVEG(),
                       tem.veg.getCURRENTVEG(),
                       tem.veg.getSUBTYPE(),
                       tem.veg.cmnt,
                       (PROP2PCT * tem.soil.getPSIPLUSC()),
                       tem.qualcon[pdyr],
                       carea,
                       cohort[pichrt].chrtarea,
                       ttotyr[pdyr],
                       tempred.mon,
                       region[pichrt] );
    }
    else
    {
      tempred.outdel( fout[i],
                      col,
                      row,
                      predname.at( i ),
                      (pichrt+1),
                      cohort[pichrt].standage,
                      tem.veg.getPOTVEG(),
                      tem.veg.getCURRENTVEG(),
                      tem.veg.getSUBTYPE(),
                      tem.veg.cmnt,
                      (PROP2PCT * tem.soil.getPSIPLUSC()),
                      tem.qualcon[pdyr],
                      carea,
                      cohort[pichrt].chrtarea,
                      ttotyr[pdyr],
                      tempred.mon,
                      region[pichrt] );
    }
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int MITelmnt44::transqc( int& maxyears,
                         int& totyr,
                         double plantc[CYCLE] )
{

  int dm;
  int qc;
  double sumcarbon = ZERO;
  qc = ACCEPT;

  if( totyr < 0 || totyr >= maxyears ) { return qc = 30; }

  for( dm = 0; dm < CYCLE; ++dm ) { sumcarbon += plantc[dm]; }

  if( sumcarbon <= 0.1 ) { return qc = TQCZEROFLAG; }

  return qc;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */
void MITelmnt44::updateChangedCohort( const int& pichrt,
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
                                      double pabandonedRCVLIN[NLVL] )                               

{
  int dlyr;


  cohort[pichrt].y[tem.I_VEGC] = ((cohort[pichrt].y[tem.I_VEGC]  
                                 * (double) pinitarea) 
                                 + (pabandonedVeg.getSTRUCTC() 
                                 * ppropAbandonedArea))
                                 / (double) cohort[pichrt].chrtarea;

  cohort[pichrt].y[tem.I_STRN] = ((cohort[pichrt].y[tem.I_STRN] 
                                 * (double) pinitarea)
                                 + (pabandonedVeg.getSTRUCTN() 
                                 * ppropAbandonedArea))
                                 / (double) cohort[pichrt].chrtarea;

  cohort[pichrt].y[tem.I_STON] = ((cohort[pichrt].y[tem.I_STON] 
                                 * (double) pinitarea)
                                 + (pabandonedVeg.getLABILEN() 
                                 * ppropAbandonedArea))
                                 / (double) cohort[pichrt].chrtarea;
    
  cohort[pichrt].y[tem.I_SOLC] = ((cohort[pichrt].y[tem.I_SOLC] 
                                 * (double) pinitarea)
                                 + (pconvertedSoil.getORGC() 
                                 * ppropConvertedArea)
                                 + (pabandonedSoil.getORGC() 
                                 * ppropAbandonedArea))
                                 / (double) cohort[pichrt].chrtarea;

  cohort[pichrt].y[tem.I_SOLN] = ((cohort[pichrt].y[tem.I_SOLN] 
                                 * (double) pinitarea)
                                 + (pconvertedSoil.getORGN() 
                                 * ppropConvertedArea)
                                 + (pabandonedSoil.getORGN() 
                                 * ppropAbandonedArea))
                                 / (double) cohort[pichrt].chrtarea;

  cohort[pichrt].y[tem.I_AVLN] = ((cohort[pichrt].y[tem.I_AVLN] 
                                 * (double) pinitarea)
                                 + (pconvertedSoil.getAVLN() 
                                 * ppropConvertedArea)
                                 + (pabandonedSoil.getAVLN() 
                                 * ppropAbandonedArea))
                                 / (double) cohort[pichrt].chrtarea;

  cohort[pichrt].NEMnsolc = ((cohort[pichrt].NEMnsolc 
                            * (double) pinitarea)
                            + (pconvertedNONSOLC 
                            * ppropConvertedArea)
                            + (pabandonedNONSOLC 
                            * ppropAbandonedArea))
                            / (double) cohort[pichrt].chrtarea;

  for( dlyr = 0; dlyr < NLVL; ++dlyr )
  {
    cohort[pichrt].NEManh4in[dlyr] = ((cohort[pichrt].NEManh4in[dlyr] 
                                     * (double) pinitarea)
                                     + (pconvertedANH4IN[dlyr] 
                                     * ppropConvertedArea)
                                     + (pabandonedANH4IN[dlyr] 
                                     * ppropAbandonedArea))
                                     / (double) cohort[pichrt].chrtarea;
                                    
    cohort[pichrt].NEMano3in[dlyr] = ((cohort[pichrt].NEMano3in[dlyr] 
                                     * (double) pinitarea) 
                                     + (pconvertedANO3IN[dlyr] 
                                     * ppropConvertedArea)
                                     + (pabandonedANO3IN[dlyr] 
                                     * ppropAbandonedArea))
                                     / (double) cohort[pichrt].chrtarea;
                             
    cohort[pichrt].NEMdphumin[dlyr] = ((cohort[pichrt].NEMdphumin[dlyr] 
                                      * (double) pinitarea)
                                      + (pconvertedDPHUMIN[dlyr] 
                                      * ppropConvertedArea)
                                      + (pabandonedDPHUMIN[dlyr] 
                                      * ppropAbandonedArea))
                                      / (double) cohort[pichrt].chrtarea;
                                    
    cohort[pichrt].NEMocin[dlyr] = ((cohort[pichrt].NEMocin[dlyr] 
                                    * (double) pinitarea)
                                   + (pconvertedOCIN[dlyr] 
                                   * ppropConvertedArea)
                                   + (pabandonedOCIN[dlyr] 
                                   * ppropAbandonedArea))
                                   / (double) cohort[pichrt].chrtarea;

    cohort[pichrt].NEMrclin[dlyr] = ((cohort[pichrt].NEMrclin[dlyr] 
                                    * (double) pinitarea)
                                    + (pconvertedRCLIN[dlyr] 
                                    * ppropConvertedArea)
                                    + (pabandonedRCLIN[dlyr] 
                                    * ppropAbandonedArea))
                                    / (double) cohort[pichrt].chrtarea;

    cohort[pichrt].NEMrcrin[dlyr] = ((cohort[pichrt].NEMrcrin[dlyr] 
                                    * (double) pinitarea)
                                    + (pconvertedRCRIN[dlyr] 
                                    * ppropConvertedArea)
                                    + (pabandonedRCRIN[dlyr] 
                                    * ppropAbandonedArea))
                                    / (double) cohort[pichrt].chrtarea;

    cohort[pichrt].NEMrcvlin[dlyr] = ((cohort[pichrt].NEMrcvlin[dlyr] 
                                     * (double) pinitarea)
                                     + (pconvertedRCVLIN[dlyr] 
                                     * ppropConvertedArea)
                                     + (pabandonedRCVLIN[dlyr] 
                                     * ppropAbandonedArea))
                                     / (double) cohort[pichrt].chrtarea;
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITelmnt44::updateTEMmonth( const int& equil,
                                 const int& totsptime,
                                 const int& pdyr,
                                 const int& pdm,
                                 const int& pichrt )
{  

  // Pass cohort characteristics information to TEM

  getTEMCohortState( pichrt );

/*  if( 1 == pdyr && 30 == (pichrt+1) && -2.0 == row )
  {
    tem.dbug = 1;
    tem.veg.dbug = 1;

    cout << " pdyr = " << pdyr;
    cout << " pdm = " << pdm;
    cout << " pichrt = " << pichrt;
    cout << endl;
  }
  else 
  { 
    tem.dbug = 0; 
    tem.veg.dbug = 0;
  }
*/

  // Assign land cover data to cohort

//  tem.ag.setGROWDD( cohort[pichrt].aggrowdd );

//  tem.ag.setSLASHC( cohort[pichrt].slash.carbon );
//  tem.ag.setSLASHN( cohort[pichrt].slash.nitrogen );

//  tem.ag.setVCONVRTFLXC( cohort[pichrt].vconvrtflx.carbon );
//  tem.ag.setVCONVRTFLXN( cohort[pichrt].vconvrtflx.nitrogen );

//  tem.ag.setSCONVRTFLXC( cohort[pichrt].sconvrtflx.carbon );
//  tem.ag.setSCONVRTFLXN( cohort[pichrt].sconvrtflx.nitrogen );

//  tem.ag.setCONVRTFLXC( cohort[pichrt].vconvrtflx.carbon
//                        + cohort[pichrt].sconvrtflx.carbon );

//  tem.ag.setCONVRTFLXN( cohort[pichrt].vconvrtflx.nitrogen
//                        + cohort[pichrt].sconvrtflx.nitrogen );

//  tem.ag.setNVRETENT( cohort[pichrt].nvretent );
//  tem.ag.setNSRETENT( cohort[pichrt].nsretent );
//  tem.ag.setNRETENT( cohort[pichrt].nvretent
//                     + cohort[pichrt].nsretent );


  if( 1 == cohort[pichrt].agstate 
      && 0 == cohort[pichrt].prevchrtarea 
      && cohort[pichrt].chrtarea > 0 )
  {
    tem.ag.setCROPPRVLEAFMX( 1.0 );

    tem.ag.setCROPTOPT( tem.veg.getTOPT() );

    tem.ag.setCROPPRVPETMX( tem.atms.getPRVPETMX() );

    tem.ag.setCROPPRVEETMX( tem.soil.getPRVEETMX() );

    tem.ag.setPRVCROPNPP( ZERO );
  }


  // Reset fluxes estimated by TEM to zero

  tem.resetMonthlyELMNTFluxes();

  if( pdm > 0 )
  {
    tem.ag.setFORMPROD10C( ZERO );
    tem.ag.setFORMPROD10N( ZERO );

    tem.ag.setFORMPROD100C( ZERO );
    tem.ag.setFORMPROD100N( ZERO );
  }


  cohort[pichrt].tqc = ACCEPT; 


  // Check TEM input for valid data 

  cohort[pichrt].qc = temgisqc( cohort[pichrt].chrtarea,
                                tem.soil.getPCTSILT(),
                                tem.soil.getPCTCLAY(),
                                tem.veg.cmnt,
                                tem.elev,
                                tem.atms.getNIRR(),
                                tem.atms.getPAR(),
                                tem.atms.getTAIR(),
                                tem.atms.getPREC(),
                                tem.atms.getPET(),
                                tem.soil.getEET(),
                                tem.soil.getMOIST(),
                                tem.soil.getSNOWPACK(),
                                tem.soil.getSURFRUN(),
                                tem.soil.getDRAINAGE(),
                                tem.atms.getCO2(),
                                tem.atms.getAOT40() );


  // Check TEM parameters for specific vegetation types

  if( ACCEPT == cohort[pichrt].qc ) 
  { 
    cohort[pichrt].qc = tem.ecdqc( tem.veg.cmnt ); 
  }


  if( ACCEPT == cohort[pichrt].qc )
  {
    tem.baseline = 0;

    tem.wrtyr = -99;

    tem.totyr = tem.startyr 
                - totsptime - 1 
                + (pdyr * tem.diffyr);


    // Allow optimum N fertilization of crops for land cover 16

     if( 16 == tem.veg.getCURRENTVEG() )
     {
       tem.ag.fertflag = 1;
     }
     else
     {
       tem.ag.fertflag = 0;
     }


    // Run the Terrestrial Ecosystem Model (TEM) under 
    //   transient conditions

    wrtyr = tem.monthlyTransient( pdyr, 
                                  pdm, 
                                  tem.tol );
   

    // Save TEM output to telmnt[0].output

    outputTEMmonth( pichrt, pdm );

    ttotyr[pdyr] = tem.totyr;
  } // End of qc == ACCEPT and tqc = ACCEPT
  else
  {
    // Set tqc flag to assign zero to all TEM variables 
    //   during simulation
      	
    cohort[pichrt].tqc = TQCZEROFLAG; 


    // Set missing values to telmnt[0].output

    setTEMmiss( pdyr,
                equil,
                totsptime,
                pichrt );
  }  


  // Save TEM state for cohort
  
  saveTEMCohortState( pichrt );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void MITelmnt44::writeBinaryCohortState( ofstream& ofstate,
                                         const int& pichrt )
{
  int ichrt = pichrt;
  
  ofstate.write( (char *)(&col), sizeof( col ) );
  ofstate.write( (char *)(&row), sizeof( row ) );
  ofstate.write( (char *)(&ichrt), sizeof( ichrt ) ); 
  ofstate.write( (char *)(&cohort[pichrt]), sizeof( MITElmntCohort44 ) );
  	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void MITelmnt44::writeCohortState( ofstream& ofstate,
                                   const int& pichrt )
{
  int i;
  int dlyr;
  int dm;

  ofstate << col << " ";    

  ofstate << row << " ";

  ofstate << (pichrt+1) << " ";
  
  ofstate << cohort[pichrt].srcCohort << " ";

  ofstate << cohort[pichrt].standage << " ";

  ofstate << cohort[pichrt].chrtarea << " ";

  ofstate << cohort[pichrt].prevchrtarea << " ";

  ofstate << cohort[pichrt].potveg << " ";

  ofstate << cohort[pichrt].currentveg << " ";

  ofstate << cohort[pichrt].subtype << " ";

  ofstate << cohort[pichrt].cmnt << " ";

 	
  for( i = 0; i < MAXSTATE; ++i )
  {
    ofstate << cohort[pichrt].y[i] << " ";
    ofstate << cohort[pichrt].prevy[i] << " ";
  }

  ofstate << cohort[pichrt].agcmnt << " ";

  ofstate << cohort[pichrt].aggrowdd << " ";

  ofstate << cohort[pichrt].agkd << " ";

  ofstate << cohort[pichrt].agprvstate << " ";

  ofstate << cohort[pichrt].agstate << " ";

  ofstate << cohort[pichrt].c2n << " ";
  
  ofstate << cohort[pichrt].cneven << " ";

  ofstate << cohort[pichrt].convrtflx.carbon << " ";
  ofstate << cohort[pichrt].convrtflx.nitrogen << " ";

  ofstate << cohort[pichrt].cropprveetmx << " ";

  ofstate << cohort[pichrt].cropprvleafmx << " ";

  ofstate << cohort[pichrt].cropprvpetmx << " ";

  ofstate << cohort[pichrt].cropResidue.carbon << " ";
  ofstate << cohort[pichrt].cropResidue.nitrogen << " ";

  ofstate << cohort[pichrt].croptopt << " ";

  ofstate << cohort[pichrt].distmnthcnt << " ";

  ofstate << cohort[pichrt].disturbflag << " ";
  
  ofstate << cohort[pichrt].disturbmonth << " ";
  
  ofstate << cohort[pichrt].eetmx << " ";

  ofstate << cohort[pichrt].fertflag << " ";                              

  ofstate << cohort[pichrt].firemnthcnt << " ";                              

  ofstate << cohort[pichrt].firendep << " ";                              

  ofstate << cohort[pichrt].formPROD10.carbon << " ";
  ofstate << cohort[pichrt].formPROD10.nitrogen << " ";

  ofstate << cohort[pichrt].formPROD100.carbon << " ";
  ofstate << cohort[pichrt].formPROD100.nitrogen << " ";

  ofstate << cohort[pichrt].fprevozone << " ";

  ofstate << cohort[pichrt].FRI << " ";

  for( dm = 0; dm < CYCLE; ++dm )
  {  
    ofstate << cohort[pichrt].initPROD1[dm].carbon << " ";
    ofstate << cohort[pichrt].initPROD1[dm].nitrogen << " ";
  }
  
  for( i = 0; i < 10; ++i )
  {
    ofstate << cohort[pichrt].initPROD10[i].carbon << " ";
    ofstate << cohort[pichrt].initPROD10[i].nitrogen << " ";
  }
    
  for( i = 0; i < 100; ++i )
  {
    ofstate << cohort[pichrt].initPROD100[i].carbon << " ";
    ofstate << cohort[pichrt].initPROD100[i].nitrogen << " ";
  }

  ofstate << cohort[pichrt].irrgflag << " ";                              
  
  ofstate << cohort[pichrt].kd << " ";

  ofstate << cohort[pichrt].MDMnpp << " ";

  ofstate << cohort[pichrt].natprveetmx << " ";

  ofstate << cohort[pichrt].natprvleafmx << " ";

  ofstate << cohort[pichrt].natprvpetmx << " ";

  ofstate << cohort[pichrt].natseedC << " ";

  ofstate << cohort[pichrt].natseedSTRN << " ";

  ofstate << cohort[pichrt].natseedSTON << " ";

  ofstate << cohort[pichrt].natsoil << " ";

  ofstate << cohort[pichrt].nattopt << " ";

  ofstate << cohort[pichrt].natyreet << " ";

  ofstate << cohort[pichrt].natyrpet << " ";


  // ********************* NEM *********************************
  
  for( dlyr = 0; dlyr < NLVL; ++dlyr )
  {
    ofstate <<  cohort[pichrt].NEManh4in[dlyr] << " "; 
  
    ofstate << cohort[pichrt].NEMano3in[dlyr] << " ";

    ofstate << cohort[pichrt].NEMdphumin[dlyr] << " "; 

    ofstate << cohort[pichrt].NEMocin[dlyr] << " ";

    ofstate << cohort[pichrt].NEMrclin[dlyr] << " "; 
  
    ofstate << cohort[pichrt].NEMrcrin[dlyr] << " ";
  
    ofstate << cohort[pichrt].NEMrcvlin[dlyr] << " ";
  }
  
  ofstate << cohort[pichrt].NEMnsolc << " ";

  //ofstate << cohort[pichrt].NEMtopdens << " ";
  
  //ofstate << cohort[pichrt].NEMtopksat << " ";
  
  //ofstate << cohort[pichrt].NEMtoppor << " ";

  // ***********************************************************
  
  
  ofstate << cohort[pichrt].newleafmx << " ";

  ofstate << cohort[pichrt].newtopt << " ";

  ofstate << cohort[pichrt].nretent << " ";

  ofstate << cohort[pichrt].nsretent << " ";

  ofstate << cohort[pichrt].nvretent << " ";

  ofstate << cohort[pichrt].petmx << " ";

  ofstate << cohort[pichrt].prev2tair << " ";

  ofstate << cohort[pichrt].prevco2 << " ";

  ofstate << cohort[pichrt].prevCropResidue.carbon << " ";
  ofstate << cohort[pichrt].prevCropResidue.nitrogen << " ";

  ofstate << cohort[pichrt].prevPROD1.carbon << " ";
  ofstate << cohort[pichrt].prevPROD1.nitrogen << " ";

  ofstate << cohort[pichrt].prevPROD10.carbon << " ";
  ofstate << cohort[pichrt].prevPROD10.nitrogen << " ";

  ofstate << cohort[pichrt].prevPROD100.carbon << " ";
  ofstate << cohort[pichrt].prevPROD100.nitrogen << " ";

//  ofstate << cohort[pichrt].prevspack << " ";

  ofstate << cohort[pichrt].prevtair << " ";

  ofstate << cohort[pichrt].prevunrmleaf << " ";
  
//  ofstate << cohort[pichrt].prod10par << " "; 

//  ofstate << cohort[pichrt].prod100par << " "; 

  ofstate << cohort[pichrt].productYear << " ";
  
  ofstate << cohort[pichrt].prvcropnpp << " ";

  ofstate << cohort[pichrt].prveetmx << " ";

  ofstate << cohort[pichrt].prvleafmx << " ";

  ofstate << cohort[pichrt].prvpetmx << " ";

  ofstate << cohort[pichrt].qc << " ";

//  ofstate << cohort[pichrt].sconvert << " "; 
  
  ofstate << cohort[pichrt].sconvrtflx.carbon << " ";
  ofstate << cohort[pichrt].sconvrtflx.nitrogen << " ";
  
  ofstate << cohort[pichrt].slash.carbon << " ";
  ofstate << cohort[pichrt].slash.nitrogen << " ";

//  ofstate << cohort[pichrt].slashpar << " "; 
   
  ofstate << cohort[pichrt].tillflag << " ";                           

  ofstate << cohort[pichrt].topt << " ";

  ofstate << cohort[pichrt].tqc << " ";

//  ofstate << cohort[pichrt].vconvert << " "; 

  ofstate << cohort[pichrt].vconvrtflx.carbon << " ";
  ofstate << cohort[pichrt].vconvrtflx.nitrogen << " ";

//  ofstate << cohort[pichrt].vrespar << " "; 

  ofstate << cohort[pichrt].yragstubC << " ";
  
  ofstate << cohort[pichrt].yragstubN << " ";
    
  ofstate << cohort[pichrt].yrcflux << " ";
  
  ofstate << cohort[pichrt].yrCH4csmp << " ";
  
  ofstate << cohort[pichrt].yrCH4ems << " ";
  
  ofstate << cohort[pichrt].yrCH4flx << " ";
  
  ofstate << cohort[pichrt].yrCO2dnflx << " ";
  
  ofstate << cohort[pichrt].yrCO2nflx << " ";
  
  ofstate << cohort[pichrt].yrconvrtC << " ";
  
  ofstate << cohort[pichrt].yrconvrtN << " ";
  
  ofstate << cohort[pichrt].yrdecayPROD1C << " ";
  
  ofstate << cohort[pichrt].yrdecayPROD10C << " ";
  
  ofstate << cohort[pichrt].yrdecayPROD100C << " ";

  ofstate << cohort[pichrt].yrdecayPROD1N << " ";
  
  ofstate << cohort[pichrt].yrdecayPROD10N << " ";
  
  ofstate << cohort[pichrt].yrdecayPROD100N << " ";
  
  ofstate << cohort[pichrt].yrdecayTOTPRODC << " ";

  ofstate << cohort[pichrt].yrdecayTOTPRODN << " ";

  ofstate << cohort[pichrt].yreet << " ";
  
  ofstate << cohort[pichrt].yrfertn << " ";
  
  ofstate << cohort[pichrt].yrfluxResidueC << " ";
  
  ofstate << cohort[pichrt].yrfluxResidueN << " ";

  ofstate << cohort[pichrt].yrformPROD1C << " ";
  
  ofstate << cohort[pichrt].yrformPROD10C << " ";
  
  ofstate << cohort[pichrt].yrformPROD100C << " ";

  ofstate << cohort[pichrt].yrformPROD1N << " ";
  
  ofstate << cohort[pichrt].yrformPROD10N << " ";
  
  ofstate << cohort[pichrt].yrformPROD100N << " ";
  
  ofstate << cohort[pichrt].yrformResidueC << " ";
  
  ofstate << cohort[pichrt].yrformResidueN << " ";

  ofstate << cohort[pichrt].yrformTOTPRODC << " ";

  ofstate << cohort[pichrt].yrformTOTPRODN << " ";

  ofstate << cohort[pichrt].yrfpc << " ";

  ofstate << cohort[pichrt].yrgpp << " ";
  
  ofstate << cohort[pichrt].yrgpr << " ";
  
//  ofstate << cohort[pichrt].yrH2Oyld << " ";

  ofstate << cohort[pichrt].yrimmob << " ";
  
//  ofstate << cohort[pichrt].yrineet << " ";
  
  ofstate << cohort[pichrt].yringpp << " ";
  
  ofstate << cohort[pichrt].yrinnpp << " ";
  
//  ofstate << cohort[pichrt].yrirrig << " ";

  ofstate << cohort[pichrt].yrlai << " ";

  ofstate << cohort[pichrt].yrleaf << " ";

  ofstate << cohort[pichrt].yrltrfalC << " ";

  ofstate << cohort[pichrt].yrltrfalN << " ";  	

  ofstate << cohort[pichrt].yrN2flx << " ";
  
  ofstate << cohort[pichrt].yrN2Odnflx << " ";
  
  ofstate << cohort[pichrt].yrN2Oflx << " ";
  
  ofstate << cohort[pichrt].yrN2Onflx << " ";
  
  ofstate << cohort[pichrt].yrnce << " ";
  
  ofstate << cohort[pichrt].yrnep << " ";
  
  ofstate << cohort[pichrt].yrninput << " ";
  
  ofstate << cohort[pichrt].yrnlost << " ";
  
  ofstate << cohort[pichrt].yrnmin << " ";
    
  ofstate << cohort[pichrt].yrnpp << " ";

  ofstate << cohort[pichrt].yrnrent << " ";
  
  ofstate << cohort[pichrt].yrnsrent << " ";

  ofstate << cohort[pichrt].yrnvrent << " ";
  
  ofstate << cohort[pichrt].yrpet << " ";
  
//  ofstate << cohort[pichrt].yrrgrndh2o << " ";
 
  ofstate << cohort[pichrt].yrrh << " ";
  
//  ofstate << cohort[pichrt].yrrain << " ";
  
//  ofstate << cohort[pichrt].yrrperc << " ";

//  ofstate << cohort[pichrt].yrrrun << " ";

  ofstate << cohort[pichrt].yrsconvrtC << " ";

  ofstate << cohort[pichrt].yrsconvrtN << " ";
  
//  ofstate << cohort[pichrt].yrsgrndh2o << " ";

  ofstate << cohort[pichrt].yrslashC << " ";

  ofstate << cohort[pichrt].yrslashN << " ";

//  ofstate << cohort[pichrt].yrsnowfall << " ";

//  ofstate << cohort[pichrt].yrsnowinf << " ";

//  ofstate << cohort[pichrt].yrsnowpack << " ";

//  ofstate << cohort[pichrt].yrsoilavlH2O << " ";
  
  ofstate << cohort[pichrt].yrsoilavlN << " ";

//  ofstate << cohort[pichrt].yrsoilC2N << " ";

//  ofstate << cohort[pichrt].yrsoilmoist << " ";
  
  ofstate << cohort[pichrt].yrsoilorgC << " ";
  
  ofstate << cohort[pichrt].yrsoilorgN << " ";

//  ofstate << cohort[pichrt].yrsoilpctp << " ";

//  ofstate << cohort[pichrt].yrsoilvsm << " ";

//  ofstate << cohort[pichrt].yrsperc << " ";

//  ofstate << cohort[pichrt].yrsrun << " ";

  ofstate << cohort[pichrt].yrtotalC << " ";
  
  ofstate << cohort[pichrt].yrunleaf << " ";

  ofstate << cohort[pichrt].yrvconvrtC << " ";
  
  ofstate << cohort[pichrt].yrvconvrtN << " ";

  ofstate << cohort[pichrt].yrvegC << " ";

//  ofstate << cohort[pichrt].yrvegC2N << " ";
    
  ofstate << cohort[pichrt].yrveginnup << " ";

  ofstate << cohort[pichrt].yrveglup << " ";
  
  ofstate << cohort[pichrt].yrvegN << " ";
  
  ofstate << cohort[pichrt].yrvegnmobil << " ";

  ofstate << cohort[pichrt].yrvegnrsorb << " ";
    
  ofstate << cohort[pichrt].yrvegnup << " ";

  ofstate << cohort[pichrt].yrvegrgrowth << " ";

  ofstate << cohort[pichrt].yrvegrmaint << " ";

  ofstate << cohort[pichrt].yrvegSTON << " ";
  
  ofstate << cohort[pichrt].yrvegSTRN << " ";

  ofstate << cohort[pichrt].yrvegsup << " ";

  ofstate << endl;
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITelmnt44::writecflx( ofstream& ofile, 
                            string varname, 
                            const int& dyr,
                            const int& year, 
                            MITdata44 cflx[MAXRTIME][MXMITNLAT] )
{

  int dm;

  int dmlat;

  string mname;

  ofile.setf( ios::fixed, ios::floatfield );
  ofile.setf( ios::showpoint );
  ofile.precision(0);


  for( dm = 0; dm < CYCLE; ++dm )
  {
    switch( dm )
    {
      case 0:  mname = "JAN"; break;
      case 1:  mname = "FEB"; break;
      case 2:  mname = "MAR"; break;
      case 3:  mname = "APR"; break;
      case 4:  mname = "MAY"; break;
      case 5:  mname = "JUN"; break;
      case 6:  mname = "JUL"; break;
      case 7:  mname = "AUG"; break;
      case 8:  mname = "SEP"; break;
      case 9:  mname = "OCT"; break;
      case 10: mname = "NOV"; break;
      case 11: mname = "DEC"; break;
    }
    ofile << year << "  ";

    ofile << mname << "  ";

    ofile << varname << endl;;

    for( dmlat = 0; dmlat < (MXMITNLAT/2); ++dmlat )
    {
      ofile << "  " << setprecision( 2 ) << setw( 8 );
      
      ofile << cflx[dyr][dm].latband[dmlat];
    }
    
    ofile << endl;

    for( dmlat = (MXMITNLAT/2); dmlat < MXMITNLAT; ++dmlat )
    {
      ofile << "  " << setprecision( 2 ) << setw( 8 );
      
      ofile << cflx[dyr][dm].latband[dmlat];
    }
    
    ofile << endl;
  }

};
