/* *************************************************************
REGMITTEM44C.CPP - Extrapolation version of the Terrestrial 
                     Ecosystem Model Version 4.4 used with 
                     interactive version of TEM (IMITTEM44.CXX)
****************************************************************


Modifications:

20060130 - DWK created by modifying regmittem436a.cpp
20060130 - DWK updated code in askpred(), setClmPred() and 
           setTEMPred()
20060201 - DWK added include regmittem437.h and standard 
           includes
20060201 - DWK deleted initializeLCLUCregion()
20060202 - DWK deleted setClmPred()
20060202 - DWK added inheritance of tempredmap( NUMTEM ) to
           RegMITTEM()
20060202 - DWK deleted updateMITTEMregionYr()
20060420 - DWK added compiler directive STANDALONE_TEM
20070129 - DWK changed include regmittem438.h to regmittem44.h
20080130 - DWK changed include from regmittem44.h to
           regmittem44a.h
20080827 - DWK changed include from regmitt44a.h to regmittem44b.h
20110705 - DWK added code to limit RTIME to MAXRTIME in 
           setRunTime()
20110705 - DWK changed include from regmitt44b.h to regmittem44c.h
20110705 - DWK changed updateLCLUCregion() to update MXCOHRTS file
           only at beginning of simulation
20110714 - DWK added relationship between igsmveg = 16 and
           clmveg = 15 in setMITCLM()
20110714 - DWK changed relationship between igsmveg = 33 and 
           clmveg = 15 to clmveg = 16 in setMITCLM()
20111019 - DWK added PET to updateMITCLMregion()
20111020 - DWK deleted mapping of igsmveg to clmveg in 
           setMITCLM() as CLM now provides data for all IGSMVEG
           types
                        
****************************************************************
************************************************************* */

// When using the Borland compiler DEFINE the compiler directive 
//   BORLAND_CPP below. Otherwise, DEFINE the compiler 
//   directive ANSI_CPP below. 

//#define BORLAND_CPP
#define ANSI_CPP

// NOTE: If running TEM "offline", make sure the compiler 
//         directive STANDALONE_TEM is DEFINED below.  If  
//         coupling TEM to the IGSM, make sure STANDALONE_TEM is 
//         NOT DEFINED (i.e. comment out next line)

#include "preproc.h"

#include<cstdio>

  using std::fopen;
  using std::fscanf;
  using std::FILE;

#include<iostream>

  using std::cerr;
  using std::cin;
  using std::cout;
  using std::ios;
  using std::endl;

#include<fstream>
  
  using std::ifstream;
  using std::ofstream;

#include<iomanip>

  using std::setprecision;

#include<cstdlib>

  using std::exit;

#include<vector>
  
  using std::vector;

#include<cmath>

#include<cctype>

  using std::toupper;

#include<string>

  using std::string;

#include<sstream>

  using std::ostringstream;

#ifdef ANSI_CPP

  #include<ctime>

  using std::time_t;
  using std::ctime;

#endif

#ifdef BORLAND_CPP

  #include<time>

  using std::time_t;
  using std::ctime;

#endif


#include "regmittem44c.h"

/* *************************************************************
************************************************************* */

RegMITTEM::RegMITTEM() : tempredmap( NUMTEM ) { };  

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

RegMITTEM::~RegMITTEM() // Destructor
{

  int i;

// Finished processing all elements - close open files

  if( fatalerr != 0 )
  {
    if( elmnt.grdcnt != -99 && elmnt.count <= elmnt.grdcnt )
    {
      cout << "FATAL ERROR! Program Terminated" << endl;
    }
    flog1 << "FATAL ERROR! Program Terminated" << endl;
  }
 
  if( 0 == telmnt[0].lonlatflag ) { fclose( flonlat ); }
  
   
  if( 1 == temflag )
  {
    fclose( fstxt );
    fclose( felev );
    
    for( i = 0; i < telmnt[0].ntempred; ++i ) 
    { 
      ftempred[i].close(); 
    }
    
    if( istateflag != 0 ) { ifstate.close(); }
    if( ostateflag != 0 ) { ofstate.close(); }
  }


  if( 1 == telmnt[0].tem.ch4flag 
      || 1 == telmnt[0].tem.n2oflag )
  {
    fclose( fslayer );
  }

  
  cout << "Closed all files!" << endl << endl;
  flog1 << "Closed all files!" << endl << endl;

  flog1.close();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int RegMITTEM::askpred( const vector<string>& pvarname, 
                        const int& posspred,
                        vector<string>& predmap ) 
{
  const int MAXCOLUMNS = 7;
  int count = 0;

  int numpred;
  
  int i;
  int j;
  int cnt;
  int length;
  
  cout << endl << endl;
  cout << "           POSSIBLE OUTPUT VARIABLES:";
  cout << endl << endl;

  vector<string>::const_iterator dn;
  
  for( dn = pvarname.begin(); dn != pvarname.end(); ++dn )
  {
    cout << std::setw(10) << *dn << " ";
    
    ++count;
  
    if( 0 == count%MAXCOLUMNS ) { cout << endl; }
  }

  cout << endl << endl;

  flog1 << endl << endl;
  flog1 << "           POSSIBLE OUTPUT VARIABLES:";
  flog1 << endl << endl;

  count = 0;
  
  for( dn = pvarname.begin(); dn != pvarname.end(); ++dn )
  {
    flog1 << std::setw(10) << *dn << " ";
    
    ++count;
    
    if( 0 == count%MAXCOLUMNS ) { flog1 << endl; }
  }

  flog1 << endl << endl;

  cout << endl << endl;
  cout << "How many variables are to be mapped (max ";
  cout << posspred << ") in output files?  ";

  cin >> numpred;

  cout << numpred << endl;

  flog1 << endl << endl;
  flog1 << "How many variables are to be mapped (max ";
  flog1 << posspred << ") in output files?";
  flog1 << numpred << endl << endl;

  cout << "Please enter output variable: " << endl;
  flog1 << "Please enter output variable: " << endl;

  for( i = 0; i < numpred; ++i )
  {
    cnt = i + 1;
    cout << cnt << " ";

    cin >> predmap.at( i );

    length = predmap.at( i ).length();
    
    for( j = 0; j < length; ++j ) 
    { 
      predmap.at( i ).at( j ) = toupper( predmap.at( i ).at( j ) ); 
    }
    
    flog1 << cnt << " " << predmap.at( i ) << endl;
  }

  return numpred;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int RegMITTEM::coregTime( const string& varname1, 
                          const int& year1, 
                          const string& varname2, 
                          const int& year2 )
{

  int fatalerr = 0;

  if( year1 != year2 )
  {
    fatalerr = 1;

    cout << "ERROR:  " << varname1 << " data and ";
    cout << varname2 << "data are not coregistered." << endl;
    cout << "Year in " << varname1 << " is " << year1;
    cout << " Year in " << varname2 << " is " << year2;
    cout << endl << endl;

    flog1 << "ERROR:  " << varname1 << " data and ";
    flog1 << varname2 << "data are not coregistered." << endl;
    flog1 << "Year in " << varname1 << " is " << year1;
    flog1 << " Year in " << varname2 << " is " << year2;
    flog1 << endl << endl;

  }

  return fatalerr;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void RegMITTEM::getGriddedLULCparameters (  const int& pigrd,
                                            const int& pichrt,
                                           const int& pdyr )
{                                 
  telmnt[pigrd].cohort[pichrt].srcCohort = (pichrt+1);
  telmnt[pigrd].cohort[pichrt].standage = MAXSTANDAGE;
  telmnt[pigrd].cohort[pichrt].potveg = pichrt;
  telmnt[pigrd].cohort[pichrt].currentveg = pichrt;
  telmnt[pigrd].cohort[pichrt].subtype = pichrt;

  telmnt[pigrd].cohort[pichrt].tillflag = 0;
  telmnt[pigrd].cohort[pichrt].irrgflag = 0;
  telmnt[pigrd].cohort[pichrt].disturbflag = 0;
  telmnt[pigrd].cohort[pichrt].disturbmonth = 0;
  telmnt[pigrd].cohort[pichrt].FRI = MAXFRI;
  telmnt[pigrd].region[pichrt] = "GLOBE";

  telmnt[pigrd].cohort[pichrt].cmnt = telmnt[pigrd].lcluc.getCommunityType( pichrt );

  telmnt[pigrd].cohort[pichrt].agcmnt = telmnt[pigrd].cohort[pichrt].cmnt;

  switch( pichrt )
  {
    case 0:  telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.50; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.50;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.00; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.00; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 1:  telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.33; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.40;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.20; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.07; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 2:  telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.33; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.40;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.20; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.07; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 3:  telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.33; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.40;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.20; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.07; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 4:  telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.33; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.40;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.27; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.00; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 5:  telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.33; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.40;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.20; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.07; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 6:  telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.50; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.40;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.10; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.00; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 7:  telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.33; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.40;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.20; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.07; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 8:  telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.33; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.40;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.20; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.07; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 9:  telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.50; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.40;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.10; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.00; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 10: telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.50; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.40;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.10; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.00; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 11: telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.50; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.50;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.00; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.00; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 12: telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.50; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.50;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.00; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.00; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 13: telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.50; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.50;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.00; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.00; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 14: telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.50; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.50;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.00; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.00; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 15: telmnt[pigrd].cohort[pichrt].agstate = 1;
             telmnt[pigrd].cohort[pichrt].agprvstate = 1;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.00; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.00;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.00; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.00; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.00; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 16: telmnt[pigrd].cohort[pichrt].agstate = 1;
             telmnt[pigrd].cohort[pichrt].agprvstate = 1;
             telmnt[pigrd].cohort[pichrt].fertflag = 1;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.00; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.00;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.00; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.00; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.00; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 17: telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.33; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.40;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.27; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.00; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 18: telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.50; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.50;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.00; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.00; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 19: telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.33; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.40;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.20; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.07; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 20: telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.50; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.50;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.00; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.00; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 21: telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.33; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.40;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.20; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.07; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 22: telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.50; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.50;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.00; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.00; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 23: telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.33; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.40;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.27; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.00; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 24: telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.50; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.50;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.00; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.00; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 25: telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.50; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.50;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.00; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.00; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 26: telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.33; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.40;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.27; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.00; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 27: telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.50; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.50;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.00; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.00; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 28: telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.33; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.40;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.20; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.07; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 29: telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.50; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.50;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.00; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.00; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.01; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    case 33: telmnt[pigrd].cohort[pichrt].agstate = 2;
             telmnt[pigrd].cohort[pichrt].agprvstate = 2;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.00; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.00;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.00; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.00; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.00; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
             break;
    default: telmnt[pigrd].cohort[pichrt].agstate = 0;
             telmnt[pigrd].cohort[pichrt].agprvstate = 0;
             telmnt[pigrd].cohort[pichrt].fertflag = 0;
             telmnt[pigrd].cohort[pichrt].slashpar = 0.00; 
             telmnt[pigrd].cohort[pichrt].vconvert = 0.00;
             telmnt[pigrd].cohort[pichrt].prod10par = 0.00; 
             telmnt[pigrd].cohort[pichrt].prod100par = 0.00; 
             telmnt[pigrd].cohort[pichrt].vrespar = 0.00; 
             telmnt[pigrd].cohort[pichrt].sconvert = 0.00;
  }

    
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void RegMITTEM::initElmntBatch( void )
{
  int dv;
  int igrd;

  for( igrd = 1; igrd < MAXNGRD; ++igrd )
  {
    if( 1 == temflag )
    {
      telmnt[igrd].ntempred = telmnt[0].ntempred;

      telmnt[igrd].tem.startyr = telmnt[0].tem.startyr;

      telmnt[igrd].lonlatflag = telmnt[0].lonlatflag;

      telmnt[igrd].tem.ag.tlulcflag = telmnt[0].tem.ag.tlulcflag;

      telmnt[igrd].tem.ch4flag = telmnt[0].tem.ch4flag;
      telmnt[igrd].tem.mdmflag = telmnt[0].tem.mdmflag;
      telmnt[igrd].tem.n2oflag = telmnt[0].tem.n2oflag;

      // Initialize CO2 for all grid cells

      telmnt[igrd].tem.atms.setINITCO2( telmnt[0].tem.atms.getINITCO2() );
      telmnt[igrd].tem.atms.setCO2LEVEL( telmnt[0].tem.atms.getCO2LEVEL() );

      telmnt[igrd].tem.veg.setDC2N( telmnt[0].tem.veg.getDC2N() );

      // Soil texture specific TEM parameters

      telmnt[igrd].tem.soil.setPCTPORA( telmnt[0].tem.soil.getPCTPORA() );
      telmnt[igrd].tem.soil.setPCTPORB( telmnt[0].tem.soil.getPCTPORB() );
      telmnt[igrd].tem.soil.setFLDCAPA( telmnt[0].tem.soil.getFLDCAPA() );
      telmnt[igrd].tem.soil.setFLDCAPB( telmnt[0].tem.soil.getFLDCAPB() );
      telmnt[igrd].tem.soil.setWILTPTA( telmnt[0].tem.soil.getWILTPTA() );
      telmnt[igrd].tem.soil.setWILTPTB( telmnt[0].tem.soil.getWILTPTB() );


      for( dv = 1; dv < MAXCMNT; ++dv )
      {
  	// Vegetation-specific rooting depth parameters
  	
  	telmnt[igrd].tem.soil.setROOTZA( telmnt[0].tem.soil.getROOTZA( dv ),
  	                                 dv );
  	                                 
	telmnt[igrd].tem.soil.setROOTZB( telmnt[0].tem.soil.getROOTZB( dv ),
	                                 dv );
	                                 
	telmnt[igrd].tem.soil.setROOTZC( telmnt[0].tem.soil.getROOTZC( dv ),
	                                 dv );
	                                 
	telmnt[igrd].tem.soil.setMINROOTZ( telmnt[0].tem.soil.getMINROOTZ( dv ),
	                                   dv );
	
	// Vegetation-specific leaf phenology parameters
	 
	telmnt[igrd].tem.veg.setMINLEAF( telmnt[0].tem.veg.getMINLEAF( dv ),
	                                 dv );
	                                 
	telmnt[igrd].tem.veg.setALEAF( telmnt[0].tem.veg.getALEAF( dv ),
	                               dv );
	                               
	telmnt[igrd].tem.veg.setBLEAF( telmnt[0].tem.veg.getBLEAF( dv ),
	                               dv );
	                               
	telmnt[igrd].tem.veg.setCLEAF( telmnt[0].tem.veg.getCLEAF( dv ),
	                               dv );

        // Site-specific TEM parameters

	telmnt[igrd].tem.setVEGCA( telmnt[0].tem.getVEGCA( dv ),
	                           dv );
	                           
	telmnt[igrd].tem.setVEGCB( telmnt[0].tem.getVEGCB( dv ),
	                           dv );
	                           
	telmnt[igrd].tem.setSTRNA( telmnt[0].tem.getSTRNA( dv ),
	                           dv );
	                           
	telmnt[igrd].tem.setSTRNB( telmnt[0].tem.getSTRNB( dv ),
	                           dv );
	                           
	telmnt[igrd].tem.setSOLCA( telmnt[0].tem.getSOLCA( dv ),
	                           dv );
	                           
	telmnt[igrd].tem.setSOLCB( telmnt[0].tem.getSOLCB( dv ),
	                           dv );
	                           
	telmnt[igrd].tem.setSOLNA( telmnt[0].tem.getSOLNA( dv ),
	                           dv );
	                           
	telmnt[igrd].tem.setSOLNB( telmnt[0].tem.getSOLNB( dv ),
	                           dv );
	                           
	telmnt[igrd].tem.setAVLNA( telmnt[0].tem.getAVLNA( dv ),
	                           dv );
	                           
	telmnt[igrd].tem.setAVLNB( telmnt[0].tem.getAVLNB( dv ),
	                           dv );
	                           
	telmnt[igrd].tem.setSTONA( telmnt[0].tem.getSTONA( dv ),
	                           dv );
	                           
	telmnt[igrd].tem.setSTONB( telmnt[0].tem.getSTONB( dv ),
	                           dv );
	                           
	telmnt[igrd].tem.veg.setUNLEAF12( telmnt[0].tem.veg.getUNLEAF12( dv ),
	                                  dv );
	                                  
	telmnt[igrd].tem.veg.setINITLEAFMX( telmnt[0].tem.veg.getINITLEAFMX( dv ),
	                                    dv );
	                                    
	telmnt[igrd].tem.veg.setCMAXCUT( telmnt[0].tem.veg.getCMAXCUT( dv ),
	                                 dv );
	                                 
	telmnt[igrd].tem.veg.setCMAX1A( telmnt[0].tem.veg.getCMAX1A( dv ),
	                                dv );
	                                
	telmnt[igrd].tem.veg.setCMAX1B( telmnt[0].tem.veg.getCMAX1B( dv ),
	                                dv );
	                                
	telmnt[igrd].tem.veg.setCMAX2A( telmnt[0].tem.veg.getCMAX2A( dv ),
	                                dv );
	                                
	telmnt[igrd].tem.veg.setCMAX2B( telmnt[0].tem.veg.getCMAX2B( dv ),
	                                dv );
	                                
	telmnt[igrd].tem.veg.setCFALL( telmnt[0].tem.veg.getCFALL( dv ),
	                               dv );
	                               
	telmnt[igrd].tem.veg.setKRA( telmnt[0].tem.veg.getKRA( dv ),
	                             dv );
	                             
	telmnt[igrd].tem.veg.setKRB( telmnt[0].tem.veg.getKRB( dv ),
	                             dv );
	                             
	telmnt[igrd].tem.microbe.setKDA( telmnt[0].tem.microbe.getKDA( dv ),
	                                 dv );
	telmnt[igrd].tem.microbe.setKDB( telmnt[0].tem.microbe.getKDB( dv ),
	                                 dv );
	                                 
	telmnt[igrd].tem.microbe.setLCCLNC( telmnt[0].tem.microbe.getLCCLNC( dv ),
	                                    dv );
	                                    
	telmnt[igrd].tem.microbe.setPROPFTOS( telmnt[0].tem.microbe.getPROPFTOS( dv ),
	                                      dv );
	                                      
	telmnt[igrd].tem.veg.setNMAXCUT( telmnt[0].tem.veg.getNMAXCUT( dv ),
	                                 dv );
	                                 
	telmnt[igrd].tem.veg.setNMAX1A( telmnt[0].tem.veg.getNMAX1A( dv ),
	                                dv );
	                                
	telmnt[igrd].tem.veg.setNMAX1B( telmnt[0].tem.veg.getNMAX1B( dv ),
	                                dv );
	                                
	telmnt[igrd].tem.veg.setNMAX2A( telmnt[0].tem.veg.getNMAX2A( dv ),
	                                dv );
	                                
	telmnt[igrd].tem.veg.setNMAX2B( telmnt[0].tem.veg.getNMAX2B( dv ),
	                                dv );
	                                
	telmnt[igrd].tem.veg.setNFALL( telmnt[0].tem.veg.getNFALL( dv ),
	                               dv );
	                               
	telmnt[igrd].tem.microbe.setNUPA( telmnt[0].tem.microbe.getNUPA( dv ),
	                                  dv );
	                                  
	telmnt[igrd].tem.microbe.setNUPB( telmnt[0].tem.microbe.getNUPB( dv ),
	                                  dv );
	                                  
	telmnt[igrd].tem.soil.setNLOSS( telmnt[0].tem.soil.getNLOSS( dv ),
	                                dv );
	                                
	telmnt[igrd].tem.microbe.setNFIXPAR( telmnt[0].tem.microbe.getNFIXPAR( dv ),
	                                     dv );
	                                     
	telmnt[igrd].tem.veg.setINITCNEVEN( telmnt[0].tem.veg.getINITCNEVEN( dv ),
	                                    dv );
	                                    
	telmnt[igrd].tem.veg.setCNMIN( telmnt[0].tem.veg.getCNMIN( dv ),
	                               dv );
	                               
	telmnt[igrd].tem.veg.setC2NA( telmnt[0].tem.veg.getC2NA( dv ),
	                              dv );
	                              
	telmnt[igrd].tem.veg.setC2NB( telmnt[0].tem.veg.getC2NB( dv ),
	                              dv );
	                              
	telmnt[igrd].tem.veg.setC2NMIN( telmnt[0].tem.veg.getC2NMIN( dv ),
	                                dv );
	                                
	telmnt[igrd].tem.microbe.setCNSOIL( telmnt[0].tem.microbe.getCNSOIL( dv ),
	                                    dv );
	                                    
	telmnt[igrd].tem.veg.setO3PARA( telmnt[0].tem.veg.getO3PARA( dv ),
	                                dv );
	                                
	telmnt[igrd].tem.veg.setO3PARB( telmnt[0].tem.veg.getO3PARB( dv ),
	                                dv );
	                                
	telmnt[igrd].tem.veg.setO3PARC( telmnt[0].tem.veg.getO3PARC( dv ),
	                                dv );

        // Vegetation-specific TEM parameters

	telmnt[igrd].tem.veg.setKC( telmnt[0].tem.veg.getKC( dv ),
	                            dv );
	                            
	telmnt[igrd].tem.veg.setKI( telmnt[0].tem.veg.getKI( dv ),
	                            dv );
	                            
	telmnt[igrd].tem.veg.setGVA( telmnt[0].tem.veg.getGVA( dv ),
	                             dv );
	                             
	telmnt[igrd].tem.veg.setTMIN( telmnt[0].tem.veg.getTMIN( dv ),
	                              dv );
	                              
	telmnt[igrd].tem.veg.setTOPTMIN( telmnt[0].tem.veg.getTOPTMIN( dv ),
	                                 dv );
	                                 
	telmnt[igrd].tem.veg.setTOPTMAX( telmnt[0].tem.veg.getTOPTMAX( dv ),
	                                 dv );
	                                 
	telmnt[igrd].tem.veg.setTMAX( telmnt[0].tem.veg.getTMAX( dv),
	                              dv );
	                              
	telmnt[igrd].tem.veg.setRAQ10A0( telmnt[0].tem.veg.getRAQ10A0( dv ),
	                                 dv );
	                                 
	telmnt[igrd].tem.veg.setRAQ10A1( telmnt[0].tem.veg.getRAQ10A1( dv ),
	                                 dv );
	                                 
	telmnt[igrd].tem.veg.setRAQ10A2( telmnt[0].tem.veg.getRAQ10A2( dv ),
	                                 dv );
	                                 
	telmnt[igrd].tem.veg.setRAQ10A3( telmnt[0].tem.veg.getRAQ10A3( dv ),
	                                 dv );
	                                 
	telmnt[igrd].tem.veg.setKN1( telmnt[0].tem.veg.getKN1( dv ),
	                             dv );
	                             
	telmnt[igrd].tem.veg.setLABNCON( telmnt[0].tem.veg.getLABNCON( dv ),
	                                 dv );

        telmnt[igrd].tem.veg.setLEAFMXC( telmnt[0].tem.veg.getLEAFMXC( dv ),
                                        dv );
                                        
	telmnt[igrd].tem.veg.setKLEAFC( telmnt[0].tem.veg.getKLEAFC( dv ),
	                                dv );
	                                
	telmnt[igrd].tem.veg.setSLA( telmnt[0].tem.veg.getSLA( dv ),
	                             dv );
	                             
	telmnt[igrd].tem.veg.setCOV( telmnt[0].tem.veg.getCOV( dv ),
	                             dv );
	                             
	telmnt[igrd].tem.veg.setFPCMAX( telmnt[0].tem.veg.getFPCMAX( dv ),
	                                dv );

        // Vegetation-specific TEM parameters for microbes

	telmnt[igrd].tem.microbe.setRHQ10( telmnt[0].tem.microbe.getRHQ10( dv ),
	                                   dv );
	                                   
	telmnt[igrd].tem.microbe.setKN2( telmnt[0].tem.microbe.getKN2( dv ),
	                                 dv );
	                                 
	telmnt[igrd].tem.microbe.setMOISTMIN( telmnt[0].tem.microbe.getMOISTMIN( dv ),
	                                      dv );
	                                      
	telmnt[igrd].tem.microbe.setMOISTOPT( telmnt[0].tem.microbe.getMOISTOPT( dv ),
	                                      dv );
	                                      
	telmnt[igrd].tem.microbe.setMOISTMAX( telmnt[0].tem.microbe.getMOISTMAX( dv ),
	                                      dv );

        // Vegetation-specific TEM parameters for agricultural systems

	                                 
	telmnt[igrd].tem.ag.setSLASHPAR( telmnt[0].tem.ag.getSLASHPAR( dv ),
	                                 dv );

	telmnt[igrd].tem.ag.setVCONVERT( telmnt[0].tem.ag.getVCONVERT( dv ),
	                                 dv );

	telmnt[igrd].tem.ag.setPROD10PAR( telmnt[0].tem.ag.getPROD10PAR( dv ),
	                                  dv );

	telmnt[igrd].tem.ag.setPROD100PAR( telmnt[0].tem.ag.getPROD100PAR( dv ),
	                                   dv );

	telmnt[igrd].tem.ag.setSCONVERT( telmnt[0].tem.ag.getSCONVERT( dv ),
	                                 dv );

	telmnt[igrd].tem.ag.setNVRETCONV( telmnt[0].tem.ag.getNVRETCONV( dv ),
	                                  dv );
	                                  
	telmnt[igrd].tem.ag.setNSRETCONV( telmnt[0].tem.ag.getNSRETCONV( dv ),
	                                  dv );
	                                  
	telmnt[igrd].tem.ag.setVRESPAR( telmnt[0].tem.ag.getVRESPAR( dv ),
	                                dv );

	telmnt[igrd].tem.ag.setTILLFACTOR( telmnt[0].tem.ag.getTILLFACTOR( dv ),
	                                   dv );
	                                   
	telmnt[igrd].tem.ag.setHARVSTC( telmnt[0].tem.ag.getHARVSTC( dv ),
	                                dv );
	                                
	telmnt[igrd].tem.ag.setHARVSTN( telmnt[0].tem.ag.getHARVSTN( dv ),
	                                dv );
	                                
	telmnt[igrd].tem.ag.setRESIDUEC( telmnt[0].tem.ag.getRESIDUEC( dv ),
	                                 dv );
	                                 
	telmnt[igrd].tem.ag.setRESIDUEN( telmnt[0].tem.ag.getRESIDUEN( dv ),
	                                 dv );
	                                 
	telmnt[igrd].tem.ag.setCROPSEEDC( telmnt[0].tem.ag.getCROPSEEDC( dv ),
	                                  dv );
	                                  
	telmnt[igrd].tem.ag.setCROPSEEDSTRN( telmnt[0].tem.ag.getCROPSEEDSTRN( dv ),
	                                     dv );
	                                     
	telmnt[igrd].tem.ag.setCROPSEEDSTON( telmnt[0].tem.ag.getCROPSEEDSTON( dv ),
	                                     dv );


        if( 1 == telmnt[0].tem.mdmflag )
        {
          // Initialize MDM dryMethanotroph parameters for all grid cells
        
          telmnt[igrd].tem.microbe.mdm.dryMethanotroph.setOMAX( telmnt[0].tem.microbe.mdm.dryMethanotroph.getOMAX( dv ),
                                                                dv );
                                      
          telmnt[igrd].tem.microbe.mdm.dryMethanotroph.setKC( telmnt[0].tem.microbe.mdm.dryMethanotroph.getKC( dv ),
                                                              dv );
                                                    
          telmnt[igrd].tem.microbe.mdm.dryMethanotroph.setOCH4Q10( telmnt[0].tem.microbe.mdm.dryMethanotroph.getOCH4Q10( dv ),
                                                                   dv );
                                                         
          telmnt[igrd].tem.microbe.mdm.dryMethanotroph.setKO( telmnt[0].tem.microbe.mdm.dryMethanotroph.getKO( dv ),
                                                              dv );
                                                    
          telmnt[igrd].tem.microbe.mdm.dryMethanotroph.setOXI_C( telmnt[0].tem.microbe.mdm.dryMethanotroph.getOXI_C( dv ),
                                                                 dv );
                                                       
          telmnt[igrd].tem.microbe.mdm.dryMethanotroph.setAFP( telmnt[0].tem.microbe.mdm.dryMethanotroph.getAFP( dv ),
                                                               dv );
                                                     
          telmnt[igrd].tem.microbe.mdm.dryMethanotroph.setMVMAX( telmnt[0].tem.microbe.mdm.dryMethanotroph.getMVMAX( dv ),
                                                                 dv );
                                                       
          telmnt[igrd].tem.microbe.mdm.dryMethanotroph.setMVMIN( telmnt[0].tem.microbe.mdm.dryMethanotroph.getMVMIN( dv ),
                                                                 dv );
                                                       
          telmnt[igrd].tem.microbe.mdm.dryMethanotroph.setMVOPT( telmnt[0].tem.microbe.mdm.dryMethanotroph.getMVOPT( dv ),
                                                                 dv );
                                                       
          telmnt[igrd].tem.microbe.mdm.dryMethanotroph.setOXIREF( telmnt[0].tem.microbe.mdm.dryMethanotroph.getOXIREF( dv ),
                                                                  dv );


          // Initialize MDM wetMethanotroph parameters for all grid cells
        
          telmnt[igrd].tem.microbe.mdm.wetMethanotroph.setOMAX( telmnt[0].tem.microbe.mdm.wetMethanotroph.getOMAX( dv ),
                                                                dv );
                                      
          telmnt[igrd].tem.microbe.mdm.wetMethanotroph.setKC( telmnt[0].tem.microbe.mdm.wetMethanotroph.getKC( dv ),
                                                              dv );
                                                    
          telmnt[igrd].tem.microbe.mdm.wetMethanotroph.setOCH4Q10( telmnt[0].tem.microbe.mdm.wetMethanotroph.getOCH4Q10( dv ),
                                                                   dv );
                                                         
          telmnt[igrd].tem.microbe.mdm.wetMethanotroph.setKO( telmnt[0].tem.microbe.mdm.wetMethanotroph.getKO( dv ),
                                                              dv );
                                                    
          telmnt[igrd].tem.microbe.mdm.wetMethanotroph.setOXI_C( telmnt[0].tem.microbe.mdm.wetMethanotroph.getOXI_C( dv ),
                                                                 dv );
                                                       
          telmnt[igrd].tem.microbe.mdm.wetMethanotroph.setAFP( telmnt[0].tem.microbe.mdm.wetMethanotroph.getAFP( dv ),
                                                               dv );
                                                     
          telmnt[igrd].tem.microbe.mdm.wetMethanotroph.setMVMAX( telmnt[0].tem.microbe.mdm.wetMethanotroph.getMVMAX( dv ),
                                                                 dv );
                                                       
          telmnt[igrd].tem.microbe.mdm.wetMethanotroph.setMVMIN( telmnt[0].tem.microbe.mdm.wetMethanotroph.getMVMIN( dv ),
                                                                 dv );
                                                       
          telmnt[igrd].tem.microbe.mdm.wetMethanotroph.setMVOPT( telmnt[0].tem.microbe.mdm.wetMethanotroph.getMVOPT( dv ),
                                                                 dv );
                                                       
          telmnt[igrd].tem.microbe.mdm.wetMethanotroph.setOXIREF( telmnt[0].tem.microbe.mdm.wetMethanotroph.getOXIREF( dv ),
                                                                  dv );


          // Initialize MDM methanogen parameters for all grid cells
          
          telmnt[igrd].tem.microbe.mdm.methanogen.setMGO( telmnt[0].tem.microbe.mdm.methanogen.getMGO( dv ),
                                                          dv );
                                                          
          telmnt[igrd].tem.microbe.mdm.methanogen.setKC( telmnt[0].tem.microbe.mdm.methanogen.getKC( dv ),                                               
                                                         dv );
                                                         
          telmnt[igrd].tem.microbe.mdm.methanogen.setPMETHANEQ( telmnt[0].tem.microbe.mdm.methanogen.getPMETHANEQ( dv ),
                                                                dv );
                                                               
          telmnt[igrd].tem.microbe.mdm.methanogen.setMAXFRESH( telmnt[0].tem.microbe.mdm.methanogen.getMAXFRESH( dv ),
                                                               dv ); 
                                                               
          telmnt[igrd].tem.microbe.mdm.methanogen.setLOWB( telmnt[0].tem.microbe.mdm.methanogen.getLOWB( dv ),
                                                           dv );
                                                           
          telmnt[igrd].tem.microbe.mdm.methanogen.setPROREF( telmnt[0].tem.microbe.mdm.methanogen.getPROREF( dv ),
                                                             dv );


          // Initialize MDM soil parameters for all grid cells

          telmnt[igrd].tem.microbe.mdm.soil.setPA( telmnt[0].tem.microbe.mdm.soil.getPA( dv ), 
                                                   dv );
        } 

        if( 1 == telmnt[0].tem.n2oflag )
        {
          // Initialize NEM N2O parameters for all grid cells

          telmnt[igrd].tem.microbe.nem.setINITDNTRF( telmnt[0].tem.microbe.nem.getINITDNTRF( dv ),
                                                     dv );
        }
      }
    }
  }

};

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */

void RegMITTEM::initializeMITTEMregion( const int& pdm )
{

  int ichrt;
  long igrd;
  
  time_t timer;
 
  int kdm;

  int dyr = 0;

//  const double gCH4togC = 12.0 / 16.0;

  int cultIntensity;

  int dday;

  int dlyr; 

  int iwet;

  double outch4flx;
   
  double pctinundated;
 
  double soilLayerSH2O[10];

  //double soilLayerTSOIL[10];
  
  double totalSoilC;

  // Water table depth in millimeters
  double watertableZ;

  // Load monthly climate data into climate, initAET, 
  //   initSH2O, initSNOWPACK, initSURFRUN and initDRAINAGE 
  //   arrays.  After gathering 12 months of data run
  //   telmnt[igrd].setTEMequilState() to initialize TEM
  //   for each grid element
  
  for( igrd = 0; igrd < mxnumgrid; ++igrd )
  {

    // Pass land cover information to telmnt[igrd]

    telmnt[igrd].col = 0.0;  
    telmnt[igrd].row = -90.0 + (4.0 * igrd);
    telmnt[igrd].maxcohorts = MAXCHRTS;
    telmnt[igrd].natcohorts = 1;

    telmnt[igrd].carea = (long) round( landcover4tem_.elmentArea[igrd] 
                         * landcover4tem_.landFrac[igrd] );

    telmnt[igrd].contnent = "GLOBE";

    // Determine maximum number of cohorts from previous year
      
    telmnt[igrd].prvmxcohrts = MAXCHRTS;


    if( (0 == istateflag && 0 == pdm) 
        || (istateflag > 0 && (CYCLE-1) == pdm) )
    {
      // Get spatially explicit soil texture and elevation
    
      telmnt[igrd].setGIStopography( flog1, 
                                     fatalerr, 
                                     fstxt,
                                     fslayer, 
                                     felev );

      // Copy parameters from first grid cell 9i.e. telmnt[0] to 
      //   all grid cells in study region

      initElmntBatch();
    }
   
    for( ichrt = 0; ichrt < telmnt[igrd].maxcohorts; ++ichrt )
    {
      // Determine land cover characteristics of cohort

      getGriddedLULCparameters ( igrd, ichrt, dyr );

      telmnt[igrd].cohort[ichrt].cmnt = telmnt[0].lcluc.getCommunityType( telmnt[igrd].cohort[ichrt].currentveg );

      // Determine soil properties based on soil texture
          
      telmnt[igrd].tem.soil.xtext( telmnt[igrd].cohort[ichrt].cmnt, 
                                   telmnt[igrd].tem.soil.getPCTSILT(), 
                                   telmnt[igrd].tem.soil.getPCTCLAY() );

      setMITCLMgrid( dyr, pdm, igrd, ichrt );

       
      if( 0 == ichrt )
      {
        telmnt[igrd].climate[telmnt[igrd].mitclm.I_NIRR][pdm] = telmnt[igrd].tem.atms.getNIRR();

        telmnt[igrd].climate[telmnt[igrd].mitclm.I_PAR][pdm] = telmnt[igrd].tem.atms.getPAR();

        telmnt[igrd].climate[telmnt[igrd].mitclm.I_TAIR][pdm] = telmnt[igrd].tem.atms.getTAIR();

        telmnt[igrd].climate[telmnt[igrd].mitclm.I_PREC][pdm] = telmnt[igrd].tem.atms.getPREC();

        telmnt[igrd].climate[telmnt[igrd].mitclm.I_CO2][pdm] = telmnt[igrd].tem.atms.getCO2();

        telmnt[igrd].climate[telmnt[igrd].mitclm.I_AOT40][pdm] = telmnt[igrd].tem.atms.getAOT40();         
      }
      
      telmnt[igrd].initPET[ichrt][pdm] = telmnt[igrd].tem.atms.getPET();

      telmnt[igrd].initAET[ichrt][pdm] = telmnt[igrd].tem.soil.getINEET();

      telmnt[igrd].initSH2O[ichrt][pdm] = telmnt[igrd].tem.soil.getMOIST();

      telmnt[igrd].initSNOWPACK[ichrt][pdm] = telmnt[igrd].tem.soil.getSNOWPACK();

      telmnt[igrd].initSURFRUN[ichrt][pdm] = telmnt[igrd].tem.soil.getSURFRUN();

      telmnt[igrd].initDRAINAGE[ichrt][pdm] = telmnt[igrd].tem.soil.getDRAINAGE();
    }
  }

  if( (CYCLE-1) == pdm )
  {    
    if( istateflag > 0 )
    {
      telmnt[0].tem.totyr = -99;

      cout << endl;
      flog1 << endl << endl;

      elmnt.show( flog1, 
                  MISSING, 
                  MISSING,
                  telmnt[0].tem.totyr, 
                  telmnt[0].tem.inittol );


// *************************************************************
// Skip to the desired record in the GIS data sets

//  grdcnt = skipRecords();

//*********************************************************** */


      if ( 1 == yrostateflag )
      {
        telmnt[0].ofstateA.open( telmnt[0].temstateAfname.c_str(), 
                                 ios::out );
      }


      // Reset aggregated CFLUX (by MIT-IGSM latitudinal band) 
      //   to zero

      telmnt[0].mitclm.resetCFLUX( pdm );

      // Reset aggregated CH4FLUX (by MIT-IGSM latitudinal band) 
      //   to zero

      telmnt[0].mitclm.resetCH4FLUX( pdm );

      // Reset aggregated N2OFLUX (by MIT-IGSM latitudinal band) 
      //   to zero

      telmnt[0].mitclm.resetN2OFLUX( pdm );


      igrd = 0;  // reinitialize grid cell count for current year

      while ( igrd < mxnumgrid && 0 == fatalerr ) // Grid cell loop
      {
        telmnt[igrd].tem.totyr = -99;
    
        telmnt[igrd].tem.setLAT( telmnt[igrd].row );
   
/* *************************************************************
		BEGIN VEGETATION MOSAIC LOOP
************************************************************* */

        for( ichrt = 0; ichrt < telmnt[igrd].maxcohorts; ++ichrt )
        {
          // Reset monthly fluxes estimated by TEM to zero

          telmnt[igrd].tem.resetMonthlyELMNTFluxes();

          telmnt[igrd].tem.ag.setFORMPROD10C( ZERO );
          telmnt[igrd].tem.ag.setFORMPROD10N( ZERO );

          telmnt[igrd].tem.ag.setFORMPROD100C( ZERO );
          telmnt[igrd].tem.ag.setFORMPROD100N( ZERO );
  

          // Read in initial TEM state determined in a previous 
          //  TEM simulation to telmnt[igrd].cohort

          telmnt[igrd].readCohortState( ifstate, ichrt );

          
          landcover4tem_.initialCohortArea[igrd][ichrt] = telmnt[igrd].cohort[ichrt].chrtarea;

          if( telmnt[igrd].cohort[ichrt].eetmx < telmnt[igrd].tem.soil.getEET() )
          {
            telmnt[igrd].cohort[ichrt].eetmx = telmnt[igrd].tem.soil.getEET();
            telmnt[igrd].cohort[ichrt].prveetmx = telmnt[igrd].cohort[ichrt].eetmx;
          }

          if( telmnt[igrd].tem.atms.getPET() < telmnt[igrd].tem.soil.getEET() )
          {
            telmnt[igrd].tem.atms.setPET( telmnt[igrd].tem.soil.getEET() );
          }

          if( telmnt[igrd].cohort[ichrt].petmx < telmnt[igrd].tem.atms.getPET() )
          {
            telmnt[igrd].cohort[ichrt].petmx = telmnt[igrd].tem.atms.getPET();
            telmnt[igrd].cohort[ichrt].prvpetmx = telmnt[igrd].cohort[ichrt].petmx;
          }


          // Determine soil properties based on soil texture
        
          telmnt[igrd].tem.soil.xtext( telmnt[igrd].tem.veg.cmnt, 
                                       telmnt[igrd].tem.soil.getPCTSILT(), 
                                       telmnt[igrd].tem.soil.getPCTCLAY() );


          if( (CROPVEG == ichrt || BIOFUELS == ichrt || PASTURE == ichrt) )
          {   
            telmnt[igrd].tem.microbe.resetEcds( telmnt[igrd].tem.veg.cmnt, 
                                                telmnt[igrd].tem.soil.getPSIPLUSC() );   

            telmnt[igrd].cohort[ichrt].kd = telmnt[igrd].tem.microbe.getKDC();

            telmnt[igrd].cohort[ichrt].agkd = telmnt[igrd].cohort[ichrt].kd;

            if( telmnt[igrd].cohort[ichrt].natsoil <= ZERO )
            {
              telmnt[igrd].cohort[ichrt].natsoil = telmnt[igrd].cohort[ichrt].y[telmnt[igrd].tem.I_SOLC];
            }

            if( telmnt[igrd].cohort[ichrt].fprevozone < ZERO ) 
            {
              telmnt[igrd].cohort[ichrt].fprevozone = 1.0;
            }
          }


          // Pass telmnt[igrd].cohort information to TEM
            
          telmnt[igrd].getTEMCohortState( ichrt );


          // Set texture dependent parameters for TEM
          
          telmnt[igrd].tem.setELMNTecd( telmnt[igrd].tem.veg.cmnt,
                                        telmnt[igrd].tem.soil.getPSIPLUSC() );
                                        

          // Run NEM for December of year previous to startyear,
          //   if CH4 and N2O fluxes are desired
          
          telmnt[igrd].tem.microbe.nem.startflag = 1;

          if ( 1 == telmnt[igrd].tem.microbe.nem.startflag 
               && 1 == telmnt[igrd].tem.ch4flag )
          {
            // Determine Methane Flux

            // Use MDM to estimate CH4 fluxes for arctic or boreal 
            //   ecosystems, otherwise, use NEM
     
            if( 1 == telmnt[igrd].tem.mdmflag
                && (2 == telmnt[igrd].tem.veg.getCURRENTVEG() 
                || 3 == telmnt[igrd].tem.veg.getCURRENTVEG()
                || 8 == telmnt[igrd].tem.veg.getCURRENTVEG()
                || 11 == telmnt[igrd].tem.veg.getCURRENTVEG()
                || 12 == telmnt[igrd].tem.veg.getCURRENTVEG()
                || 21 == telmnt[igrd].tem.veg.getCURRENTVEG()
                || 22 == telmnt[igrd].tem.veg.getCURRENTVEG()) )
            {
              switch( telmnt[igrd].tem.veg.getCURRENTVEG() )
              {
                case 19: iwet = 1;
                         pctinundated = 100.000;
                         break;

                case 20: iwet = 1;
                         pctinundated = 100.000;
                         break;

                case 21: iwet = 1;
                         pctinundated = 100.000;
                         break;

                case 22: iwet = 1;
                         pctinundated = 100.000;
                         break;
      
                default: iwet = 0;
                         pctinundated = ZERO;
              }

              if( 15 == telmnt[igrd].tem.veg.getCURRENTVEG() 
                  || 16 == telmnt[igrd].tem.veg.getCURRENTVEG() )
              {
                cultIntensity = 5;
              }
              else { cultIntensity = 1; } 
      
              telmnt[igrd].tem.soil.setCH4EMISS( ZERO );
              telmnt[igrd].tem.soil.setCH4CONSUMP( ZERO );
              telmnt[igrd].tem.soil.setCH4FLUX( ZERO );
      
              for( dday = 0; dday < (int) telmnt[igrd].tem.atms.ndays[pdm]; ++dday )
              {                 
                telmnt[igrd].tem.microbe.mdm.soil.setTSOIL( ((telmnt[igrd].tem.microbe.nem.dayTemp[0][dday] * (1.75/20.0))
                                           + (telmnt[igrd].tem.microbe.nem.dayTemp[1][dday] * (2.75/20.0))
                                           + (telmnt[igrd].tem.microbe.nem.dayTemp[2][dday] * (4.5/20.0))
                                           + (telmnt[igrd].tem.microbe.nem.dayTemp[3][dday] * (7.5/20.0))
                                           + (telmnt[igrd].tem.microbe.nem.dayTemp[4][dday] * (3.5/20.0))) );
 
                for( dlyr = 0; dlyr < (CLMNLAYERS-1); ++dlyr )
                {
//                  soilLayerTSOIL[dlyr] = telmnt[igrd].tem.microbe.nem.dayTemp[dlyr][dday];
                  soilLayerSH2O[dlyr] = telmnt[igrd].tem.microbe.nem.dayMoist[dlyr][dday];
                }

                for( dlyr = 0; dlyr < MXMDMNLAY; ++dlyr )
                {
                  telmnt[igrd].tem.microbe.mdm.soil.setINTERSOILT( telmnt[igrd].tem.microbe.mdm.soil.getTSOIL(),
                                                                   dlyr );

//                  telmnt[igrd].tem.microbe.mdm.soil.setINTERSOILT( telmnt[igrd].tem.microbe.mdm.soil.interpolateLayers( soilLayerTSOIL, 
//                                                                                                                        telmnt[igrd].tem.soil.layerThick,
//                                                                                                                        (CLMNLAYERS-1),
//                                                                 (dlyr+1) ), 
//                                                                   dlyr );

                  telmnt[igrd].tem.microbe.mdm.soil.setINTERSH2O( telmnt[igrd].tem.microbe.mdm.soil.interpolateLayers( soilLayerSH2O, 
                                                                                                                       telmnt[igrd].tem.soil.layerThick,
                                                                                                                       (CLMNLAYERS-1),
                                                                                                                       (dlyr+1) ),
                                                                                                                       dlyr );
                }

        
                if( 0 == iwet ){ watertableZ = 500.0; }
                else { watertableZ = 0.00; }

                telmnt[igrd].tem.soil.setPH( 7.5 );

                telmnt[igrd].tem.microbe.mdm.stepday( telmnt[igrd].tem.veg.cmnt,
                                                      cultIntensity,
                                                      iwet,
                                                      (pctinundated * 100.0),
                                                      watertableZ,
                                                      telmnt[igrd].tem.soil.getPCTSAND(),
                                                      telmnt[igrd].tem.soil.getPCTSILT(),
                                                      telmnt[igrd].tem.soil.getPCTCLAY(),
                                                      telmnt[igrd].tem.soil.getPCTFLDCAP(),
                                                      telmnt[igrd].tem.soil.getROOTZ(),
                                                      telmnt[igrd].tem.soil.getPH(),
                                                      telmnt[igrd].cohort[ichrt].MDMnpp,
                                                      telmnt[igrd].tem.microbe.mdm.soil.getTSOIL() );

                // Aggregate daily CH4 fluxes from MDM into monthly 
                //   fluxes
         
                telmnt[igrd].tem.soil.setCH4FLUX( (-1.0
                                                   * 0.001 
                                                   * (telmnt[igrd].tem.soil.getCH4FLUX() 
                                                   + telmnt[igrd].tem.microbe.mdm.soil.getCH4TOT())) ); 
              }
            }
            else
            {   
              telmnt[igrd].tem.soil.setCH4EMISS( telmnt[igrd].tem.microbe.nem.setMonthlyCH4emissions( telmnt[igrd].tem.veg.getCURRENTVEG(),
                                                                                                      pdm,
                                                                                                      telmnt[igrd].tem.atms.getTAIR(),
                                                                                                      telmnt[igrd].tem.atms.getPREC(),
                                                                                                      telmnt[igrd].tem.soil.getEET(),
                                                                                                      telmnt[igrd].tem.getLAT() ) );

              telmnt[igrd].tem.soil.setCH4EMISS( telmnt[igrd].tem.soil.getCH4EMISS() );

              telmnt[igrd].tem.soil.setCH4FLUX( telmnt[igrd].tem.soil.getCH4EMISS() );
            }

            // Agggregate methane fluxes for latitudinal band
            
            outch4flx = telmnt[igrd].tem.soil.getCH4FLUX() 
                        * 1000000000.0;

            telmnt[igrd].mitclm.aggregFlux2D( telmnt[igrd].tem.totyr, 
                                              telmnt[igrd].row, 
                                              (double) telmnt[igrd].cohort[ichrt].chrtarea, 
                                              outch4flx,
                                              telmnt[0].mitclm.tch4flx2D );

          } 
  

          if( 1 == telmnt[igrd].tem.n2oflag
              && telmnt[igrd].tem.getY(telmnt[igrd].tem.I_SOLC) > ZERO
              && telmnt[igrd].tem.veg.getCURRENTVEG() > 0 
              && telmnt[igrd].tem.veg.getCURRENTVEG() < 30 )
          {
            if ( 1 == telmnt[igrd].tem.microbe.nem.startflag )
            {
            if (telmnt[igrd].tem.microbe.nem.getNONSOLC() < ZERO )
            {
            telmnt[igrd].tem.microbe.nem.setNONSOLC( ZERO );
            }

              totalSoilC = telmnt[igrd].tem.microbe.nem.getNONSOLC() + telmnt[igrd].tem.getY(telmnt[igrd].tem.I_SOLC);

              //cout << " NONSOLC before NEM in regmittem = " << telmnt[igrd].tem.microbe.nem.getNONSOLC() << endl;
              //cout << " TOPPOR: " << endl;
              //cout << telmnt[igrd].tem.microbe.nem.getTOPPOR() << endl;
              //cout << " TOPDENS: " << endl;
              //cout << telmnt[igrd].tem.microbe.nem.getTOPDENS() << endl;
              //cout << " TOPKSAT: " << endl;
              //cout << telmnt[igrd].tem.microbe.nem.getTOPKSAT() << endl;
              //cout << " dayTemp: " << endl;
              //cout << microbe.nem.dayTemp << endl;
              //cout << " dayMoist: " << endl;
              //cout << microbe.nem.dayMoist << endl;
              //cout << " HourMoist: " << endl;
              //cout << microbe.nem.hourMoist << endl;
              //cout << " dayTair: " << endl;
              //cout << atms.dayTair << endl;
              //cout << " RainDur: " << endl;
              //cout << atms.rainDuration << endl;
              //cout << " RainInt: " << endl;
              //cout << atms.rainIntensity << endl;
        
              telmnt[igrd].tem.microbe.nem.setPH( 7.00 );

              setMITCLMgrid( 0 , (CYCLE-1) , igrd, ichrt );

              telmnt[igrd].tem.microbe.nem.stepmonth( pdm,
                                                      telmnt[igrd].tem.veg.getCURRENTVEG(),
                                                      telmnt[igrd].tem.veg.cmnt,
                                                      telmnt[igrd].tem.soil.getPCTCLAY(),
                                                      telmnt[igrd].tem.microbe.nem.getTOPPOR(),
                                                      telmnt[igrd].tem.microbe.nem.getTOPDENS(),
                                                      telmnt[igrd].tem.microbe.nem.getPH(),
                                                      telmnt[igrd].tem.microbe.nem.getTOPKSAT(),
                                                      totalSoilC,
                                                      telmnt[igrd].tem.microbe.nem.dayTemp,
                                                      telmnt[igrd].tem.microbe.nem.dayMoist,
                                                      telmnt[igrd].tem.microbe.nem.hourMoist,
                                                      telmnt[igrd].tem.atms.dayTair,
                                                      telmnt[igrd].tem.atms.rainDuration,
                                                      telmnt[igrd].tem.atms.rainIntensity );

              telmnt[igrd].tem.soil.setN2OFLUX( ((telmnt[igrd].tem.microbe.nem.getEN2ON() 
                                                + telmnt[igrd].tem.microbe.nem.getEN2ODN()) 
                                                * 0.0001) );

              //cout << "REGMITTEM: N2OFLUX = " << telmnt[igrd].tem.microbe.nem.getEN2ON() << endl;
              //cout << "REGMITTEM: DN2OFLUX = " << telmnt[igrd].tem.microbe.nem.getEN2ODN() << endl;
              //cout << "REGMITTEM: CAREA = " << telmnt[igrd].cohort[ichrt].chrtarea << endl;

              // Agggregate methane fluxes for latitudinal band

              telmnt[igrd].mitclm.aggregFlux2D( telmnt[igrd].tem.totyr, 
                                                telmnt[igrd].row, 
                                                telmnt[igrd].cohort[ichrt].chrtarea,
                                               (telmnt[igrd].tem.soil.getN2OFLUX() * 1000000000.0),
                                               telmnt[0].mitclm.tn2oflx2D );

            }
          }  
        }

        elmnt.show( flog1, 
                    telmnt[igrd].col, 
                    telmnt[igrd].row,
                    telmnt[igrd].tem.totyr, 
                    telmnt[igrd].tem.tol );

        ++igrd;
      }
    }
    else
    {      	    
      for( igrd = 0; igrd < mxnumgrid; ++igrd )
      {
        for( ichrt = 0; ichrt < telmnt[igrd].maxcohorts; ++ichrt )
        {
          // Pass lulcdat information to telmnt[0].cohort
    
          telmnt[igrd].cohort[ichrt].chrtarea = (long) (telmnt[igrd].carea
                                                * landcover4tem_.fracLandArea[igrd][ichrt]);  

          telmnt[igrd].cohort[ichrt].prevchrtarea = telmnt[igrd].cohort[ichrt].chrtarea;
  
          telmnt[igrd].cohort[ichrt].cmnt = telmnt[0].lcluc.getCommunityType( telmnt[igrd].cohort[ichrt].currentveg );

          telmnt[igrd].setTEMequilState( flog1,
                                         equil,
                                         totsptime,
                                         ichrt );

          if( telmnt[igrd].tem.intflag > 0 )
          {
            if( elmnt.count < elmnt.grdcnt )
            {
              cout << "Integration terminated before attaining ";
              cout << "tolerance level" << endl;
            }

            flog1 << "Integration terminated before attaining ";
            flog1 << "tolerance level" << endl;

            telmnt[igrd].tem.intflag = 0;
          }
      
          // Write out telmnt[0].cohort to output file for
          //   potential use in a future TEM simulation
     
          if( 1 == ostateflag )
          {
            telmnt[igrd].writeCohortState( ofstate, ichrt );
          }

          if( 1 == yrostateflag )
          {
            telmnt[igrd].writeCohortState( telmnt[0].ofstateA, 
                                           ichrt );
          }
      
          // Write selected TEM variables from telmnt[0].output to 
          //   outfile files
  
          if( 1 == spinoutfg || 1 == equil )
          {
            telmnt[igrd].temwritepred( ftempred, 
                                       tempredmap, 
                                       dyr, 
                                       ichrt,
                                       telmnt[0].ntempred );
          }
        } // End of cohort loop
 

        elmnt.show( flog1, 
                    telmnt[igrd].col, 
                    telmnt[igrd].row,
                    telmnt[igrd].tem.totyr, 
                    telmnt[igrd].tem.tol );
      } // End of grid cell loop
    }
  }

  if ( 1 == yrostateflag ) { telmnt[0].ofstateA.close(); }

  timer = time( NULL );
  flog1 << "Finished year " << telmnt[0].tem.totyr;
  flog1 << " at " << ctime( &timer );
  
};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void RegMITTEM::setGridParameters( void )
{
  string ifilename;
  
  
  cout << "How many elements are in the gridded data sets?";

  cin >> mxnumgrid;

  flog1 << endl << endl;
  flog1 << "How many elements are in the gridded data sets? ";
  flog1 << mxnumgrid << endl << endl;

  cout << "How do you locate your grid cells?" << endl;
  cout << "Enter 0 for column/row:" << endl;
  cout << "Enter 1 for longitude/latitude: ";

  cin >> telmnt[0].lonlatflag;

  flog1 << "How do you locate your grid cells?" << endl;
  flog1 << "Enter 0 for column/row:" << endl;
  flog1 << "Enter 1 for longitude/latitude: ";
  flog1 << telmnt[0].lonlatflag << endl << endl;

  if ( 0 == telmnt[0].lonlatflag )
  {
    cout << "Please enter the name of the file containing the";
    cout << " latitude data: " << endl;
    cout << "               (e.g., LAT.GIS) " << endl;

    cin >> ifilename;

    flog1 << "Please enter the name of the file containing the";
    flog1 << " latitude data: " << endl;
    flog1 << "               (e.g., LAT.GIS) " << endl;
    flog1 << ifilename << endl << endl;

    flonlat = fopen( ifilename.c_str(), "r" );

    if ( !flonlat )
    {
      cerr << "\nCannot open " << ifilename;
      cerr << " for data input" << endl;
      exit( -1 );
    }
  }

 // Use all records in GIS data sets (i.e. desired coverage)?

  elmnt.ask( flog1 );
	
};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void RegMITTEM::setInTEMState( void )
{

  string ifilename;
  
	
  istateflag = 0;
  istateyear = -99;

  if( 0 == equil )
  {
    cout << endl;
    cout << "Do you want to use spatially explicit data for";
    cout << " initial conditions? " << endl;

    cout << "  Enter 0 for NO:" << endl;

    cout << "  Enter 1 if spatially explicit data represents the";
    cout << " state of TEM at equilibrium conditions:" << endl;

    cout << "  Enter 2 if spatially explicit data represents the";
    cout << " state of TEM at the end of a specific year:" << endl;

    cin >> istateflag;

    flog1 << endl;
    flog1 << "Do you want to use spatially explicit data for";
    flog1 << " intial conditions? " << endl;

    flog1 << "  Enter 0 for NO:" << endl;

    flog1 << "  Enter 1 if spatially explicit data represents the";
    flog1 << " state of TEM at equilibrium conditions:" << endl;

    flog1 << "  Enter 2 if spatially explicit data represents the";
    flog1 << " state of TEM at the end of a specific year:" << endl;

    flog1 << "istateflag = " << istateflag << endl << endl;

    if( istateflag > 0 )
    {
      cout << "Please enter the name of the file containing the";
      cout << " initial TEM state data: " << endl;
      cout << "               (e.g., TEMINIT.GIS) " << endl;

      cin >> ifilename;

      flog1 << "Please enter the name of the file containing the";
      flog1 << " initial TEM state data: " << endl;
      flog1 << "               (e.g., TEMINIT.GIS) " << endl;
      flog1 << ifilename << endl << endl;

      ifstate.open( ifilename.c_str(), ios::in );

      if( !ifstate )
      {
        cout << endl;
        cout << "Cannot open " << ifilename << " for data input";
        cout << endl;

        flog1 << endl;
        flog1 << "Cannot open " << ifilename << " for data input";
        flog1 << endl;
      
        exit( -1 );
      }

      if( 2 == istateflag )
      {
        cout << "Please enter the year that you wish to use as";
        cout << " the initial TEM state: ";

        cin >> istateyear;

        flog1 << "Please enter the year that you wish to use as";
        flog1 << " the initial TEM state: ";
        flog1 << istateyear << endl << endl;
      }
    }
  }
  
};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void RegMITTEM::setMITCLMgrid( const int& outyr, 
                               const int& pdm,
                               const long& igrd,
                               const int& ichrt )
{
 
  // Correction to convert degrees Kelvin to degrees Celsius
  static const double KELVIN0 = -273.15;
  
  int dday;
  int dhr;
  int dlyr;

  int igsmveg;
    
  if ( 0 == ichrt )
  {
    // Determine climate based on MIT-IGSM latitudinal bands

    if ( 0 == pdm )
    {
      telmnt[igrd].mxtair = -999.9;
      telmnt[igrd].yrprec = ZERO;
    }   

    
    // Get monthly surface solar radiation
    
    telmnt[igrd].tem.atms.setNIRR( climate4tem_.swrs[igrd] );


    // Get monthly photosynthetically active radiation (PAR)
    
    telmnt[igrd].tem.atms.setPAR( telmnt[igrd].tem.atms.getNIRR() * 0.5 );
    
    
    // Get monthly air temperature
    
    telmnt[igrd].tem.atms.setTAIR( (climate4tem_.temp[igrd] 
                                    + KELVIN0) );

                                  
    // Determine maximum monthly air temperature for the year
    
    if ( telmnt[igrd].mxtair < telmnt[igrd].tem.atms.getTAIR() )
    {
      telmnt[igrd].mxtair = telmnt[igrd].tem.atms.getTAIR();
    }


    // Get monthly precipitation
 
    telmnt[igrd].tem.atms.setPREC( (climate4tem_.pre[igrd]
                                   * telmnt[0].mitclm.ndays[pdm]) );

    
    // Determine annual precipitation

    telmnt[igrd].yrprec += telmnt[igrd].tem.atms.getPREC();

    if ( 1 == telmnt[0].mitclm.n2oflag )
    {
      for ( dday = 0; dday < telmnt[0].mitclm.mdays[pdm]; ++dday )
      {
     	telmnt[igrd].tem.atms.dayTair[dday] 
     	     = climate4tem_.daytemp[igrd][dday];
      	           	     
      	telmnt[igrd].tem.atms.dayTair[dday] += KELVIN0;

      	telmnt[igrd].tem.atms.rainDuration[dday] 
      	     = climate4tem_.strmdur[igrd][dday];

      	telmnt[igrd].tem.atms.rainIntensity[dday] 
      	     = climate4tem_.qstrm[igrd][dday];      	
      }
    }


    // Get monthly atmospheric CO2 
        
    telmnt[igrd].tem.atms.setCO2( climate4tem_.co2[igrd] );
    
    // Get monthly AOT40 ozone values
      
    setMITgridO3( pdm, igrd );
  }    
  

  igsmveg = telmnt[igrd].cohort[ichrt].currentveg;

  // Determine hydrology inputs for each cohort within a 
  //   latitudinal band 

  telmnt[igrd].tem.atms.setPET( (climate4tem_.pet[igrd][igsmveg]
                                * telmnt[igrd].mitclm.ndays[pdm]) );

  if ( telmnt[igrd].tem.atms.getPET() < ZERO )
  {
    telmnt[igrd].tem.atms.setPET( ZERO );
  }    


  telmnt[igrd].tem.soil.setINEET( (climate4tem_.aet[igrd][igsmveg]
                                  * telmnt[igrd].mitclm.ndays[pdm]) );
           
  if ( telmnt[igrd].tem.soil.getINEET() < ZERO )
  {
    telmnt[igrd].tem.soil.setINEET( ZERO );
  }    
  

  // Adjust soil moistures to TEM rooting depths

  if ( telmnt[igrd].tem.soil.getROOTZ() <= 1.0 )
  {
    telmnt[igrd].tem.soil.setMOIST( (climate4tem_.sh2o1m[igrd][igsmveg]
                                    * telmnt[igrd].tem.soil.getROOTZ()) ); 
  }
  else if ( telmnt[igrd].tem.soil.getROOTZ() <= 2.0 )
  {
    telmnt[igrd].tem.soil.setMOIST( (climate4tem_.sh2o1m[igrd][igsmveg]
                                    * (2.0 - telmnt[igrd].tem.soil.getROOTZ()))
                                    + (climate4tem_.sh2o2m[igrd][igsmveg]
                                    * (telmnt[igrd].tem.soil.getROOTZ() - 1.0)) ); 
  }
  else 
  {        
    telmnt[igrd].tem.soil.setMOIST( (climate4tem_.sh2o2m[igrd][igsmveg]
                                    * (telmnt[igrd].tem.soil.getROOTZ() / 2.0)) );
  }

  
  if ( telmnt[igrd].tem.soil.getMOIST() < ZERO )
  {
    telmnt[igrd].tem.soil.setMOIST( ZERO );
  }    


  // Set soil moisture of wetlands to field capacity

  if ( 17 == igsmveg
       || 18 == igsmveg
       || 19 == igsmveg
       || 20 == igsmveg
       || 21 == igsmveg
       || 22 == igsmveg
       || 23 == igsmveg
       || 24 == igsmveg
       || 25 == igsmveg )
  {
    telmnt[igrd].tem.soil.setMOIST( telmnt[igrd].tem.soil.getFLDCAP() );
  }    
 
  if ( 1 == telmnt[0].tem.microbe.nem.startflag
       && (1 == telmnt[0].mitclm.ch4flag 
       || 1 == telmnt[0].mitclm.n2oflag) ) 
  { 
    for ( dday = 0; dday < telmnt[0].mitclm.mdays[pdm]; ++dday )
    {
      for ( dlyr = 0; dlyr < (CLMNLAYERS-1); ++dlyr )
      {
        telmnt[igrd].tem.microbe.nem.dayTemp[dlyr][dday] = climate4tem_.daytsoil[dlyr][igrd][igsmveg][dday];

        telmnt[igrd].tem.microbe.nem.dayTemp[dlyr][dday] += KELVIN0; 
       
        telmnt[igrd].tem.microbe.nem.dayMoist[dlyr][dday] = climate4tem_.daysh2o[dlyr][igrd][igsmveg][dday];
      }
    }
  }

   
  telmnt[igrd].tem.soil.setSNOWPACK( climate4tem_.swe[igrd][igsmveg] );

  telmnt[igrd].tem.soil.setSURFRUN( (climate4tem_.sfr[igrd][igsmveg]
                                    * telmnt[igrd].mitclm.ndays[pdm]) );

  telmnt[igrd].tem.soil.setDRAINAGE( (climate4tem_.drn[igrd][igsmveg]
                                     * telmnt[igrd].mitclm.ndays[pdm]) );


  if ( 1 == telmnt[0].tem.microbe.nem.startflag
       && 1 == telmnt[0].mitclm.n2oflag )
  {
    for ( dday = 0; dday < telmnt[0].mitclm.mdays[pdm]; ++dday )
    {
      for ( dlyr = 0; dlyr < NLAYERS; ++dlyr )
      {
        for ( dhr = 0; dhr < MAXDAYHRS; ++dhr )
        {
          telmnt[igrd].tem.microbe.nem.hourMoist[dlyr][dday][dhr] = climate4tem_.hrsh2o[dlyr][igrd][igsmveg][dday][dhr];
        }
      }                                       
    }
  }
      
  
  // Set year

  telmnt[igrd].mitclm.year = telmnt[0].mitclm.startyr 
                          - totsptime - 1 + outyr;
  telmnt[igrd].tem.totyr = telmnt[0].mitclm.year;

  telmnt[igrd].tem.atms.setMXTAIR( telmnt[igrd].mxtair );
  telmnt[igrd].tem.atms.yrprec = telmnt[igrd].yrprec;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void RegMITTEM::setMITgridO3( const int& pdm,
                              const long& igrd ) 
{
  if( 1 == telmnt[0].mitclm.to3flag )
  {    
    // Convert 3-hourly ozone latitudinal band data from 
    //    MIT climate module into AOT40 values
    
    telmnt[igrd].mitclm.mkMITaot40( climate4tem_.o3[igrd],
                                    pdm );

    if( telmnt[igrd].mitclm.getAOT40() <= ZERO )
    {
      telmnt[igrd].tem.atms.setAOT40( ZERO );
    }
    else
    {
      telmnt[igrd].tem.atms.setAOT40( telmnt[igrd].mitclm.getAOT40() );
    } 
  }
  else
  {
    telmnt[0].tem.atms.setAOT40( ZERO );
  }    
	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void RegMITTEM::setMITLULCgrid( const int& pdyr,
                                const long& pigrd )
{
  int ichrt;


  // Update prevchrtarea with area of cohort from 
  //  the previous year
        
  for( ichrt = 0; ichrt < MAXCHRTS; ++ichrt )
  {
    telmnt[pigrd].cohort[ichrt].prevchrtarea = telmnt[pigrd].cohort[ichrt].chrtarea;
  }
  

  for( ichrt = 0; ichrt < MAXCHRTS; ++ichrt )
  {       
    // Pass lulcdat information to telmnt[0].cohort
    
    telmnt[pigrd].cohort[ichrt].chrtarea = (long) round( telmnt[pigrd].carea
                                               * landcover4tem_.fracLandArea[pigrd][ichrt] );  
  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void RegMITTEM::setOutTEMState( void )
{
  string ifilename;

  ostringstream tempfname;
  
	
  cout << "Do you want to save the state of TEM as a spatially";
  cout << " explicit data set each year?";
  cout << endl;
  cout << "  Enter 0 for NO:" << endl;
  cout << "  Enter 1 for YES:" << endl;

  cin >> yrostateflag;

  flog1 << "Do you want to save the state of TEM as a spatially";
  flog1 << " explicit data set each year?";
  flog1 << endl;
  flog1 << yrostateflag << endl << endl;

  if( 1 == yrostateflag )
  {
    // Identify path to files temstateA and temstateB

    cout << "Enter the pathname to the two files that save ";
    cout << "the 'state' of TEM every other year:" << endl;

    cin >> ifilename;

    flog1 << "Enter the pathname to the two files that save ";
    flog1 << "the 'state' of TEM every other year:" << endl;
    flog1 << ifilename << endl << endl;

    tempfname.str( "" );
    tempfname << ifilename << "temstateA";
    telmnt[0].temstateAfname = tempfname.str();
             
    tempfname.str( "" );
    tempfname << ifilename << "temstateB";
    telmnt[0].temstateBfname = tempfname.str();
  }

  ostateflag = 0;
  ostateyear = -99;

  cout << endl;
  cout << "Do you want to save the state of TEM as a spatially";
  cout << " explicit data set for a specified year? " << endl;

  cout << "  Enter 0 for NO:" << endl;

  cout << "  Enter 1 for saving the state of TEM at equilibrium";
  cout << " conditions:" << endl;

  cout << "  Enter 2 for saving the state of TEM at the end of";
  cout << " a specific year:" << endl;

  cin >> ostateflag;

  flog1 << endl;
  flog1 << "Do you want to save the state of TEM as a spatially";
  flog1 << " explicit data set for a specified year? " << endl;
  flog1 << "  Enter 0 for NO:" << endl;
  flog1 << "  Enter 1 for saving the state of TEM at equilibrium";
  flog1 << " conditions:" << endl;
  flog1 << "  Enter 2 for saving the state of TEM at the end of";
  flog1 << " a specific year:" << endl;
  flog1 << ostateflag << endl << endl;

  if( ostateflag != 0 )
  {
    cout << "Please enter the name of the file to contain";
    cout << " the 'state' of TEM: " << endl;
    cout << "               (e.g., TEMSTATE.GIS) " << endl;

    cin >> ifilename;

    flog1 << "Please enter the name of the file to contain";
    flog1 << " the 'state' of TEM: " << endl;
    flog1 << "               (e.g., TEMSTATE.GIS) " << endl;
    flog1 << ifilename << endl << endl;

    ofstate.open( ifilename.c_str(), ios::out );

    if( 2 == ostateflag )
    {
      cout << "Please enter the year that you wish to save";
      cout << " the TEM state: ";

      cin >> ostateyear;

      flog1 << "Please enter the year that you wish to save";
      flog1 << " the TEM state: ";
      flog1 << ostateyear << endl << endl;
    }
  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void RegMITTEM::setRunMode( void )
{

  // Assign a log file for the simulation

  flog1.open( "mittem44.log" );

/* *************************************************************
  Run equilibrium simulation or transient simulation ?
************************************************************* */

  cout << endl;
  cout << "Do you want to run the model only for";
  cout << " steady state conditions ? " << endl;

  cout << " Enter 0 for transient simulation" << endl;
  cout << " Enter 1 for steady state simulation" << endl;

  cin >> equil;

  flog1 << endl;
  flog1 << "Do you want to run the model only for";
  flog1 << " steady state conditions ? " << endl;
  flog1 << " Enter 0 for transient simulation" << endl;
  flog1 << " Enter 1 for steady state simulation" << endl;
  flog1 << "equil = " << equil << endl << endl;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void RegMITTEM::setRunTime( void )
{
	
  RTIME = 1;

  if ( 0 == equil )
  {
    cout << "What is the first year of the transient climate ";
    cout << "data?" << endl;

    cin >> telmnt[0].mitclm.startyr;

    flog1 << "What is the first year of the transient climate ";
    flog1 << "data?" << endl;
    flog1 << telmnt[0].mitclm.startyr << endl << endl;

    telmnt[0].tem.startyr = telmnt[0].mitclm.startyr;

    flog1 << telmnt[0].mitclm.startyr << endl << endl;


    spinflag = 0;
    numspin = 0;
    spintime = 0;
    totsptime = 0;
    
    if( 0 == istateflag )
    {	
      // Start transient TEM run from equilibrium conditions 
      //   (i.e. spinflag == 0) or with a "spin up" period to 
      //   remove potential artifacts associated with changing 
      //   from equilibrium conditions to transient conditions 
      //   from model results

      cout << "Do you want to start the transient with a spin up ";
      cout << "period? " << endl;
      cout << "Enter 0 for no:" << endl;
      cout << "Enter 1 for yes: ";

      cin  >> spinflag;

      flog1 << "Do you want to start the transient with a spin up ";
      flog1 << "period? " << endl;
      flog1 << "Enter 0 for no:" << endl;
      flog1 << "Enter 1 for yes: " << endl;
      flog1 << "spinflag = " << spinflag << endl << endl;

      if( 1 == spinflag )
      {
        // Specify conditions for initializing TEM with a transient 
        //   "spin up" period for a grid cell

        cout << "How many spins do you want in the spin up ";
        cout << "period? ";
      
        cin >> numspin;
      
        flog1 << "How many spins do you want in the spin up ";
        flog1 << "period? " << endl;
        flog1 << "numspin = " << numspin << endl << endl;

        cout << "How many years per spin? ";

        cin >> spintime;

        flog1 << "How many years per spin? " << endl;
        flog1 << "spintime = " << spintime << endl << endl;

        totsptime = spintime * numspin;
      
        flog1 << "totsptime = " << totsptime << endl << endl;
      
        RTIME += totsptime;
      }
    }
    
    // Specify conditions for the "non-spin up" part of the 
    //   transient TEM run

    cout << endl;
    cout << "How many years do you run for transient ";
    cout << "simulations ? " << endl;
    
    cin >> transtime;

    flog1 << endl;
    flog1 << "How many years do you run for transient ";
    flog1 << "simulations ? " << endl;
    flog1 << "transtime = " << transtime << endl << endl;

    RTIME += transtime;
    flog1 << "RTIME = " << RTIME << endl << endl;
    
    if( RTIME > MAXRTIME )
    {
      cout << "RTIME = " << RTIME << endl;
    	cout << "RTIME is more than MAXRTIME ";
    	cout << "MAXRTIME = " << MAXRTIME << endl;
    	cout << "Increase MAXRTIME to conduct longer simulations";
    	cout << endl;

    	flog1 << "RTIME is more than MAXRTIME ";
    	flog1 << "(MAXRTIME = " << MAXRTIME << endl;
    	flog1 << "Increase MAXRTIME to conduct longer simulations";
    	flog1 << endl;
    	
    	exit( -1 );
    }
  }
  else 
  { 
    totsptime = RTIME; 
    spintime = 1;
  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void RegMITTEM::setTEMPred( void )
{

  int i;
  string tempredfile;

  telmnt[0].tem.predflag = 0;
  
  cout << endl;
  cout << "Do you wish spatially explicit output data from TEM?";
  cout << endl;
  cout << "  Enter 0 for no" << endl;
  cout << "  Enter 1 for yes: ";

  cin >> telmnt[0].tem.predflag;

  flog1 << endl;
  flog1 << "Do you wish spatially explicit output data from TEM?";
  flog1 << endl;
  flog1 << "  Enter 0 for no" << endl;
  flog1 << "  Enter 1 for yes: " << endl;
  flog1 << "telmnt[0].tem.predflag = " << telmnt[0].tem.predflag;
  flog1 << endl << endl;

  if( 1 == telmnt[0].tem.predflag )
  {
    telmnt[0].ntempred = askpred( telmnt[0].tem.predstr, 
                                  NUMTEM,
                                  tempredmap );
                                  

    for( i = 0; i < telmnt[0].ntempred; ++i )
    {
      cout << endl;
      cout << "Enter the name of the OUTPUT file to contain ";
      cout << tempredmap[i] << ":  ";
      
      cin >> tempredfile;
      
      flog1 << "Enter the name of the OUTPUT file to contain ";
      flog1 << tempredmap[i] << ":  " << tempredfile << endl;
      
      ftempred[i].open( tempredfile.c_str(), ios::out );
    }

    spinoutfg = 0;
    spinoutyrs = 0;
    
    if( 1 == spinflag )
    {
      cout << "Do you want to save TEM output from ";
      cout << "the spin-up period?" << endl;
      cout << "Enter 0 for no:" << endl;
      cout << "Enter 1 for all spin-up years: ";
      cout << "Enter 2 for some spinup years: ";
      
      cin >> spinoutfg;

      flog1 << "Do you want to save TEM output from ";
      flog1 << "the spin-up period?" << endl;
      flog1 << "Enter 0 for no:" << endl;
      flog1 << "Enter 1 for all spin-up years: ";
      flog1 << "Enter 2 for some spinup years: ";
      flog1 << "spinoutfg = " << spinoutfg << endl << endl;

      if( 2 == spinoutfg )
      {
        cout << "How many years of spin-up do you want to ";
        cout << "save in the TEM output?" << endl;
        
        cin >> spinoutyrs;
        
        flog1 << "How many years of spin-up do you want to ";
        flog1 << "save in the TEM output?" << endl;     
        flog1 << "spinoutyrs = " << spinoutyrs << endl << endl;
      }
    }
  }
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

long RegMITTEM::skipRecords( void )
{

  long elmntcnt = 0;
  string contnent;

  end1 = 1;
  fatalerr = 0;

  if( 0 == elmnt.strtflag )
  {

    for( elmntcnt = 0; elmntcnt < elmnt.numskip; ++elmntcnt )
    {       
      if( 1 == temflag )
      {
        telmnt[0].setGIStopography( flog1, 
                                    fatalerr, 
                                    fstxt,
                                    fslayer, 
                                    felev );
      }      
    }
  }
  
  return elmntcnt;  // Number of grid cells skipped!

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void RegMITTEM::starttem( void )
{
  double dc2n;
  string ifilename;  
  int numcmnt;
    
  cout << "Do you want to run the terrestrial ecosystem model ";
  cout << "(TEM)?" << endl;
  cout << "Enter 0 for No" << endl;
  cout << "Enter 1 for Yes: " << endl;

  cin >> temflag;

  flog1 << "Do you want to run the terrestrial ecosystem model ";
  flog1 << "(TEM)?" << endl;
  flog1 << "  Enter 0 for No" << endl;
  flog1 << "  Enter 1 for Yes" << endl;
  flog1 << "temflag = " << temflag << endl << endl;

  telmnt[0].ntempred = 0;

  if ( 1 == temflag )
  {
    // Hand off climate startyr and flags to TEM
    
    telmnt[0].tem.startyr = telmnt[0].mitclm.startyr;
    telmnt[0].tem.ch4flag = telmnt[0].mitclm.ch4flag;
    telmnt[0].tem.n2oflag = telmnt[0].mitclm.n2oflag;
    telmnt[0].tem.atms.setINITCO2( telmnt[0].mitclm.getINITCO2() );
    telmnt[0].tem.atms.setCO2LEVEL( telmnt[0].mitclm.getCO2LEVEL() );

    // Hand off tlulcflag from LCLUC module to TEM
    telmnt[0].tem.ag.tlulcflag = telmnt[0].lcluc.tlulcflag;
   
    telmnt[0].tem.initrun( flog1 );
    telmnt[0].tem.askODE( flog1 );


    // Get soil texture dependent parameters

    telmnt[0].tem.soil.getecd( flog1 );

    // Get vegetation type dependent parameters

    telmnt[0].tem.soil.getrootz( flog1 );
    telmnt[0].tem.veg.getecd( flog1 );
    telmnt[0].tem.veg.getleafecd( flog1 );
    telmnt[0].tem.microbe.getvegecd( flog1 );

    // Get parameters associated with human disturbance 
    //   activities

    if ( 1 == telmnt[0].tem.ag.tlulcflag )
    {
      telmnt[0].tem.ag.getecd( flog1 );
      telmnt[0].tem.ag.setAgricFlags( flog1 );
    }


    if( 1 == telmnt[0].tem.ch4flag 
        && 1 == telmnt[0].tem.mdmflag )
    {
      telmnt[0].tem.microbe.mdm.methanogen.getecd( flog1 );
      telmnt[0].tem.microbe.mdm.dryMethanotroph.getecd( flog1 );
      telmnt[0].tem.microbe.mdm.wetMethanotroph.getecd( flog1 );
      telmnt[0].tem.microbe.mdm.soil.getecd( flog1 );
    }

    if( 1 == telmnt[0].tem.n2oflag )
    {
      telmnt[0].tem.microbe.nem.getecd( flog1 );
    }

    // Get calibration site specific parameters

    cout << "Please enter the number of community types with ";
    cout << "calibration data:";

    cin >> numcmnt;

    flog1 << endl << endl << "Please enter the number of ";
    flog1 << "community types with calibration data:";
    flog1 << numcmnt << endl;

    telmnt[0].tem.getsitecd( numcmnt, flog1 );


    cout << "Please enter the name of the file containing the ";
    cout << "soil texture data:" << endl;
    cout << "               (e.g., TEXTURE.GIS) " << endl;
    
    cin >> ifilename;
    
    flog1 << "Please enter the name of the file containing the ";
    flog1 << "soil texture data:";
    flog1 << endl;
    flog1 << "               (e.g., TEXTURE.GIS) " << endl;
    flog1 << ifilename << endl << endl;

    fstxt = fopen( ifilename.c_str(), "r" );

    if ( !fstxt )
    {
      flog1 << endl << "Cannot open " << ifilename;
      flog1 << " for soil texture data input" << endl;
      
      exit( -1 );
    }

    if( 1 == telmnt[0].tem.ch4flag 
        || 1 == telmnt[0].tem.n2oflag )
    {
      cout << "Please enter the name of the file containing ";
      cout << "the soil layer characteristics data:" << endl;
      cout << "               (e.g., SOILLAYER.GIS) " << endl;
    
      cin >> ifilename;

      flog1 << "Please enter the name of the file containing ";
      flog1 << "the soil layer characteristics data:" << endl;
      flog1 << "               (e.g., SOILLAYER.GIS) " << endl;
      flog1 << ifilename << endl << endl;

      fslayer = fopen( ifilename.c_str(), "r" );

      if ( !fslayer )
      {
        flog1 << endl << "Cannot open " << ifilename;
        flog1 << " for soil layer data input" << endl;
      
        exit( -1 );
      }
    }

    cout << "Please enter the name of the file containing the ";
    cout << "elevation data: " << endl;
    cout << "               (e.g., ELEV.GIS) " << endl;
    
    cin >> ifilename;
    
    flog1 << "Please enter the name of the file containing the ";
    flog1 << "elevation data: " << endl;
    flog1 << "               (e.g., ELEV.GIS) " << endl;
    flog1 << ifilename << endl << endl;
    
    felev = fopen( ifilename.c_str(), "r" );

    if ( !felev )
    {
      flog1 << "\nCannot open " << ifilename;
      flog1 << " for data input" << endl;
      exit( -1 );
    }

    cout << endl << endl;
    cout << "Enter the factor for changing C:N per ppmv of ";
    cout << "enhanced CO2:" << endl;
    cout << "                     (Enter 0.0 for no change): ";
    cout << endl;
 
    cin >> dc2n;

    telmnt[0].tem.veg.setDC2N( dc2n );
    
    flog1 << endl;
    flog1 << "Enter the factor for changing C:N per ppmv of ";
    flog1 << "enhanced CO2:" << endl;
    flog1 << "                     (Enter 0.0 for no change): ";
    flog1 << endl;
    flog1 << "telmnt[0].tem.veg.dc2n = ";
    flog1 << telmnt[0].tem.veg.getDC2N() << endl << endl;

    setTEMPred();
 

    // Identify files to receive TEM States from particular 
    //   years in current simulation
  
    setOutTEMState();

  } // End of "temflag == 1"
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

#ifdef STANDALONE_TEM
void RegMITTEM::updateMITCLMregion( const int& pdyr,
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
                                    ifstream& ifpet,
                                    ifstream& ifeet,
                                    ifstream& ifsh2o1,
                                    ifstream& ifsh2o2,
                                    ifstream& ifspack,
                                    ifstream& ifsrfrun,
                                    ifstream& ifdrain,
                                    ifstream& ifco2,
                                    ifstream& ifo3 )
{
  int dlyr;
  int dgrd;
  int dveg;

  int dday;
  int dhr;
  int dt;

  MITdata44 tco22D;

  // Daily soil moisture across layers:
  //   0 - 0 to 1.75 cm
  //   1 - 1.75 to 4.50 cm
  //   2 - 4.50 to 9.00 cm
  //   3 - 9.00 to 16.5 cm
  //   4 - 16.5 to 29.0 cm
  //   5 - 29.0 to 49.4 cm 
  MITDayChrtData44 tdaySoilMoist[CLMNLAYERS-1][MAXMDAYS];

  MITDayClmData44 tdayTair[MAXMDAYS];

  // Daily soil temperatures across layers:
  //   0 - 0 to 1.75 cm
  //   1 - 1.75 to 4.50 cm
  //   2 - 4.50 to 9.00 cm
  //   3 - 9.00 to 16.5 cm
  //   4 - 16.5 to 29.0 cm
  //   5 - 29.0 to 49.4 cm      
  MITDayChrtData44 tdayTsoil[CLMNLAYERS-1][MAXMDAYS];     

  MITCLMdata44 tdrain2D;

  MITCLMdata44 teet2D;

  MITdata44 tnirr2D;

  O3data44 to32D;

  MITCLMdata44 tpet2D;

  MITdata44 tprec2D;

  // Daily rain duration of current event
  MITDayClmData44 trainDuration[MAXMDAYS];

  // Daily rain intensity of current event
  MITDayClmData44 trainIntensity[MAXMDAYS];

  MITCLMdata44 tsrfrun2D;

  // Soil moisture (mm) in top 1 meter
  MITCLMdata44 tsh2o1m2D;

  // Soil moisture (mm) in top 2 meters
  MITCLMdata44 tsh2o2m2D;
          
  MITCLMdata44 tspack2D;

  MITdata44 ttair2D;


  // Monthly net irradiance or surface solar radiation *********
       
  telmnt[0].mitclm.setLatBandClm( ifnirr, 
                                  tnirr2D );

  for( dgrd = 0; dgrd < MAXNGRD; ++dgrd )
  {
    telmnt[dgrd].year = tnirr2D.year; 
  }


  // Monthly air temperature ***********************************
  
  telmnt[0].mitclm.setLatBandClm( iftair, 
                                  ttair2D );

  fatalerr = coregTime( "NIRR", 
                        telmnt[0].year, 
                        "TAIR", 
                        ttair2D.year );

  if( fatalerr != 0 ) { exit( -1 ); }

  if( 1 == telmnt[0].tem.microbe.nem.startflag
       && (1 == telmnt[0].mitclm.ch4flag 
       || 1 == telmnt[0].mitclm.n2oflag) )
  {
    // Daily soil temperature of layer 1 ***********************
    
    telmnt[0].mitclm.setLatBandSoil( ifdayTsoilL1,
                                     pdm, 
                                     tdayTsoil[0] );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "DAYTSOIL1", 
                          tdayTsoil[0][0].year );

    if( fatalerr != 0 ) { exit( -1 ); }


    // Daily soil temperature of layer 2 ***********************

    telmnt[0].mitclm.setLatBandSoil( ifdayTsoilL2,
                                     pdm, 
                                     tdayTsoil[1] );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "DAYTSOIL2", 
                          tdayTsoil[1][0].year );

    if( fatalerr != 0 ) { exit( -1 ); }


    // Daily soil temperature of layer 3 ***********************

    telmnt[0].mitclm.setLatBandSoil( ifdayTsoilL3,
                                     pdm, 
                                     tdayTsoil[2] );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "DAYTSOIL3", 
                          tdayTsoil[2][0].year );

    if( fatalerr != 0 ) { exit( -1 ); }


    // Daily soil temperature of layer 4 ***********************

    telmnt[0].mitclm.setLatBandSoil( ifdayTsoilL4,
                                     pdm, 
                                     tdayTsoil[3] );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "DAYTSOIL4", 
                          tdayTsoil[3][0].year );

    if( fatalerr != 0 ) { exit( -1 ); }


    // Daily soil temperature of layer 5 ***********************

    telmnt[0].mitclm.setLatBandSoil( ifdayTsoilL5,
                                     pdm, 
                                     tdayTsoil[4] );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "DAYTSOIL5", 
                          tdayTsoil[4][0].year );

    if( fatalerr != 0 ) { exit( -1 ); }


    // Daily soil temperature of layer 6 ***********************

    telmnt[0].mitclm.setLatBandSoil( ifdayTsoilL6,
                                     pdm, 
                                     tdayTsoil[5] );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "DAYTSOIL6", 
                          tdayTsoil[5][0].year );

    if( fatalerr != 0 ) { exit( -1 ); }


    // Daily soil temperature of layer 7 ***********************

    telmnt[0].mitclm.setLatBandSoil( ifdayTsoilL7,
                                     pdm, 
                                     tdayTsoil[6] );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "DAYTSOIL7", 
                          tdayTsoil[6][0].year );

    if( fatalerr != 0 ) { exit( -1 ); }


    // Daily soil temperature of layer 8 ***********************

    telmnt[0].mitclm.setLatBandSoil( ifdayTsoilL8,
                                     pdm, 
                                     tdayTsoil[7] );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "DAYTSOIL8", 
                          tdayTsoil[7][0].year );

    if( fatalerr != 0 ) { exit( -1 ); }


    // Daily soil temperature of layer 9 ***********************

    telmnt[0].mitclm.setLatBandSoil( ifdayTsoilL9,
                                     pdm, 
                                     tdayTsoil[8] );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "DAYTSOIL9", 
                          tdayTsoil[8][0].year );

    if( fatalerr != 0 ) { exit( -1 ); }


    // Daily soil moisture of layer 1 **************************

    telmnt[0].mitclm.setLatBandSoil( ifdaySML1,
                                     pdm, 
                                     tdaySoilMoist[0] );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "DAYSH2O1", 
                          tdaySoilMoist[0][0].year );

    if( fatalerr != 0 ) { exit( -1 ); }


    // Daily soil moisture of layer 2 **************************

    telmnt[0].mitclm.setLatBandSoil( ifdaySML2,
                                     pdm, 
                                     tdaySoilMoist[1] );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "DAYSH2O2", 
                          tdaySoilMoist[1][0].year );

    if( fatalerr != 0 ) { exit( -1 ); }


    // Daily soil moisture of layer 3 **************************

    telmnt[0].mitclm.setLatBandSoil( ifdaySML3,
                                     pdm, 
                                     tdaySoilMoist[2] );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "DAYSH2O3", 
                          tdaySoilMoist[2][0].year );

    if( fatalerr != 0 ) { exit( -1 ); }


    // Daily soil moisture of layer 4 **************************

    telmnt[0].mitclm.setLatBandSoil( ifdaySML4,
                                     pdm, 
                                     tdaySoilMoist[3] );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "DAYSH2O4", 
                          tdaySoilMoist[3][0].year );

    if( fatalerr != 0 ) { exit( -1 ); }


    // Daily soil moisture of layer 5 **************************

    telmnt[0].mitclm.setLatBandSoil( ifdaySML5,
                                     pdm, 
                                     tdaySoilMoist[4] );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "DAYSH2O5", 
                          tdaySoilMoist[4][0].year );

    if( fatalerr != 0 ) { exit( -1 ); }


    // Daily soil moisture of layer 6 **************************

    telmnt[0].mitclm.setLatBandSoil( ifdaySML6,
                                     pdm, 
                                     tdaySoilMoist[5] );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "DAYSH2O6", 
                          tdaySoilMoist[5][0].year );

    if( fatalerr != 0 ) { exit( -1 ); }


    // Daily soil moisture of layer 7 **************************

    telmnt[0].mitclm.setLatBandSoil( ifdaySML7,
                                     pdm, 
                                     tdaySoilMoist[6] );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "DAYSH2O7", 
                          tdaySoilMoist[6][0].year );

    if( fatalerr != 0 ) { exit( -1 ); }


    // Daily soil moisture of layer 8 **************************

    telmnt[0].mitclm.setLatBandSoil( ifdaySML8,
                                     pdm, 
                                     tdaySoilMoist[7] );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "DAYSH2O8", 
                          tdaySoilMoist[7][0].year );

    if( fatalerr != 0 ) { exit( -1 ); }


    // Daily soil moisture of layer 9 **************************

    telmnt[0].mitclm.setLatBandSoil( ifdaySML9,
                                     pdm, 
                                     tdaySoilMoist[8] );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "DAYSH2O9", 
                          tdaySoilMoist[8][0].year );

    if( fatalerr != 0 ) { exit( -1 ); }


    // Assign daily soil profile data from input files to 
    //   climate4tem_ variables for later use with TEM
    
    for( dlyr = 0; dlyr < (CLMNLAYERS-1); ++dlyr )
    {
      for( dgrd = 0; dgrd < MAXNGRD; ++dgrd )
      {
        for( dveg = 0; dveg < CLMMXNVEG; ++dveg )
        {
          for( dday = 0; dday < telmnt[0].mitclm.mdays[pdm]; ++dday )
          {
            climate4tem_.daytsoil[dlyr][dgrd][dveg][dday] = 
                     tdayTsoil[dlyr][dday].latband[dgrd][dveg];

            climate4tem_.daysh2o[dlyr][dgrd][dveg][dday] = 
                     tdaySoilMoist[dlyr][dday].latband[dgrd][dveg];
          }
        }
      }
    }
  }
    
  if ( 1 == telmnt[0].tem.microbe.nem.startflag
       && 1 == telmnt[0].mitclm.n2oflag )
  {
    // Daily air temperature ***********************************
    
    telmnt[0].mitclm.setLatBandClm( ifdayTair, 
                                    pdm,
                                    tdayTair );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "DAYTAIR", 
                          tdayTair[0].year );

    if( fatalerr != 0 ) { exit( -1 ); }


    // Daily rain duration *************************************

    telmnt[0].mitclm.setLatBandClm( ifrainDur,
                                    pdm, 
                                    trainDuration );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "RAINDUR", 
                          trainDuration[0].year );

    if( fatalerr != 0 ) { exit( -1 ); }


    // Daily rain intensity *************************************

    telmnt[0].mitclm.setLatBandClm( ifrainInt,
                                    pdm, 
                                    trainIntensity );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "RAININT", 
                          trainIntensity[0].year );

    if( fatalerr != 0 ) { exit( -1 ); }


    // Hourly soil moisture of layer 1 *************************
    
    telmnt[0].mitclm.setLatBandSoil( ifhrSML1,
                                     pdm, 
                                     thrSoilMoist[0] );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "HRSH2O1", 
                          thrSoilMoist[0][0][0].year );

    if( fatalerr != 0 ) { exit( -1 ); }


    // Hourly soil moisture of layer 2 *************************

    telmnt[0].mitclm.setLatBandSoil( ifhrSML2,
                                     pdm, 
                                     thrSoilMoist[1] );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "HRSH2O2", 
                          thrSoilMoist[1][0][0].year );

    if( fatalerr != 0 ) { exit( -1 ); }


    // Hourly soil moisture of layer 3 *************************

    telmnt[0].mitclm.setLatBandSoil( ifhrSML3,
                                     pdm, 
                                     thrSoilMoist[2] );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "HRSH2O3", 
                          thrSoilMoist[2][0][0].year );

    if( fatalerr != 0 ) { exit( -1 ); }


    // Hourly soil moisture of layer 4 *************************

    telmnt[0].mitclm.setLatBandSoil( ifhrSML4,
                                     pdm, 
                                     thrSoilMoist[3] );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "HRSH2O4", 
                          thrSoilMoist[3][0][0].year );

    if( fatalerr != 0 ) { exit( -1 ); }


    // Hourly soil moisture of layer 5 *************************

    telmnt[0].mitclm.setLatBandSoil( ifhrSML5,
                                     pdm, 
                                     thrSoilMoist[4] );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "HRSH2O5", 
                          thrSoilMoist[4][0][0].year );

    if( fatalerr != 0 ) { exit( -1 ); }


    // Hourly soil moisture of layer 6 *************************

    telmnt[0].mitclm.setLatBandSoil( ifhrSML6,
                                     pdm, 
                                     thrSoilMoist[5] );

    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "HRSH2O6", 
                          thrSoilMoist[5][0][0].year );

    if( fatalerr != 0 ) { exit( -1 ); }

    // Assign daily data from input files to climate4tem_ 
    //   variables for later use with TEM

    for( dgrd = 0; dgrd < MAXNGRD; ++dgrd )
    {
      for( dday = 0; dday < telmnt[0].mitclm.mdays[pdm]; ++dday )
      {
        climate4tem_.daytemp[dgrd][dday] = tdayTair[dday].latband[dgrd];
        climate4tem_.strmdur[dgrd][dday] = trainDuration[dday].latband[dgrd];
        climate4tem_.qstrm[dgrd][dday] = trainIntensity[dday].latband[dgrd];
      }
    }

    // Assign hourly soil profile data from input files to 
    //   climate4tem_ variables for later use with TEM

    for( dlyr = 0; dlyr < MXNLAYERS; ++dlyr )
    {
      for( dgrd = 0; dgrd < MAXNGRD; ++dgrd )
      {
        for( dveg = 0; dveg < CLMMXNVEG; ++dveg )
        {
          for( dday = 0; dday < telmnt[0].mitclm.mdays[pdm]; ++dday )
          {
            for( dhr = 0; dhr < MXDAYHRS; ++dhr )
            {
              climate4tem_.hrsh2o[dlyr][dgrd][dveg][dday][dhr] = 
                      thrSoilMoist[dlyr][dday][dhr].latband[dgrd][dveg];
            }
          }
        }
      }
    }
  }
    

  // Monthly precipitation *************************************
  
  telmnt[0].mitclm.setLatBandClm( ifprec, 
                                  tprec2D );

  fatalerr = coregTime( "NIRR", 
                        telmnt[0].year, 
                        "PREC", 
                        tprec2D.year );

  if( fatalerr != 0 ) { exit( -1 ); }


  // Monthly potential evapotranspiration *************************
  
  telmnt[0].mitclm.setLatBandClm( ifpet, 
                                  tpet2D );

  fatalerr = coregTime( "NIRR", 
                        telmnt[0].year, 
                        "PET", 
                        tpet2D.year );

  if( fatalerr != 0 ) { exit( -1 ); }


  // Monthly actual evapotranspiration *************************
  
  telmnt[0].mitclm.setLatBandClm( ifeet, 
                                  teet2D );

  fatalerr = coregTime( "NIRR", 
                        telmnt[0].year, 
                        "AET", 
                        teet2D.year );

  if( fatalerr != 0 ) { exit( -1 ); }


  // Monthly snowpack ******************************************
  
  telmnt[0].mitclm.setLatBandClm( ifspack, 
                                  tspack2D );

  fatalerr = coregTime( "NIRR", 
                        telmnt[0].year, 
                        "SNOWPACK", 
                        tspack2D.year );

  if( fatalerr != 0 ) { exit( -1 ); }


  // Monthly surface runoff ************************************
  
  telmnt[0].mitclm.setLatBandClm( ifsrfrun, 
                                  tsrfrun2D );

  fatalerr = coregTime( "NIRR", 
                        telmnt[0].year, 
                        "SURFRUN", 
                        tsrfrun2D.year );

  if( fatalerr != 0 ) { exit( -1 ); }


  // Monthly drainage ******************************************
  
  telmnt[0].mitclm.setLatBandClm( ifdrain, 
                                  tdrain2D );

  fatalerr = coregTime( "NIRR", 
                        telmnt[0].year, 
                        "DRAINAGE", 
                        tdrain2D.year );

  if( fatalerr != 0 ) { exit( -1 ); }


  // Monthly soil moisture for top 1 meter of soil profile *****
  
  telmnt[0].mitclm.setLatBandClm( ifsh2o1, 
                                  tsh2o1m2D );

  fatalerr = coregTime( "NIRR", 
                        telmnt[0].year, 
                        "1MSH2O", 
                        tsh2o1m2D.year );

  if( fatalerr != 0 ) { exit( -1 ); }


  // Monthly soil moisture for top 2 meters of soil profile ****

  telmnt[0].mitclm.setLatBandClm( ifsh2o2, 
                                  tsh2o2m2D );

  fatalerr = coregTime( "NIRR", 
                        telmnt[0].year, 
                        "2MSH2O", 
                        tsh2o2m2D.year );

  if( fatalerr != 0 ) { exit( -1 ); }


  // Monthly atmospheric CO2 concentration**********************
  
  telmnt[0].mitclm.setLatBandClm( ifco2, 
                                  tco22D );

  if( 0 == equil )
  {
    if( pdyr > totsptime )
    {
      fatalerr = coregTime( "NIRR", 
                            telmnt[0].year, 
                            "CO2", 
                            tco22D.year );

      if( fatalerr != 0 ) { exit( -1 ); }
    }
  }

  // 3-hour atmospheric ozone concentration *******
  
  telmnt[0].mitclm.setLatBando3( ifo3, 
                                 to32D );

  if( pdyr > totsptime )
  {
    fatalerr = coregTime( "NIRR", 
                          telmnt[0].year, 
                          "AOT40", 
                          to32D.year );

    if( fatalerr != 0 ) { exit( -1 ); }
  }


  // Assign monthly data from input files to climate4tem_ 
  //   variables for later use with TEM

  for( dgrd = 0; dgrd < MAXNGRD; ++dgrd )
  {
    climate4tem_.swrs[dgrd] = tnirr2D.latband[dgrd];
    climate4tem_.temp[dgrd] = ttair2D.latband[dgrd];
    climate4tem_.pre[dgrd] = tprec2D.latband[dgrd];
    climate4tem_.co2[dgrd] = tco22D.latband[dgrd];
  }

  // Assign daily data for different land cover cohorts from 
  //   input files to climate4tem_ variables for later use 
  //   with TEM
  
  for( dgrd = 0; dgrd < MAXNGRD; ++dgrd )
  {
    for( dveg = 0; dveg < CLMMXNVEG; ++dveg )
    {
      climate4tem_.pet[dgrd][dveg] = tpet2D.latband[dgrd][dveg];
      climate4tem_.aet[dgrd][dveg] = teet2D.latband[dgrd][dveg];
      climate4tem_.sh2o1m[dgrd][dveg] = tsh2o1m2D.latband[dgrd][dveg];
      climate4tem_.sh2o2m[dgrd][dveg] = tsh2o2m2D.latband[dgrd][dveg];
      climate4tem_.swe[dgrd][dveg] = tspack2D.latband[dgrd][dveg];
      climate4tem_.sfr[dgrd][dveg] = tsrfrun2D.latband[dgrd][dveg];
      climate4tem_.drn[dgrd][dveg] = tdrain2D.latband[dgrd][dveg];
    }
  }


  // Assign 3-hour ozone data from input files to climate4tem_ 
  //   variable for later use with TEM

  for( dgrd = 0; dgrd < MAXNGRD; ++dgrd )
  {
    for( dt = 0; dt < MAX3HRS; ++dt )
    {
      climate4tem_.o3[dgrd][dt] = to32D.latbandhr[dgrd][dt]; 
    }
  }

};
#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

#ifdef STANDALONE_TEM
void RegMITTEM::updateLCLUCregion( const int& pdyr,
                                   FILE* fnumchrts,
                                   FILE* flandarea,
                                   FILE* flulc )
{
  int gisend;
  int ichrt;
  long igrd;

  LandAreadata44 landareadat;
  LulcCohortdata44 lulcdat[MAXCHRTS];  
  MaxCohortdata43 mxcohrtdat;

  fatalerr = 0;
  igrd = 0;

  if( 0 == pdyr )
  {
    while ( igrd < MAXNGRD && 0 == fatalerr )
    {    

      // Get the total number of cohorts in the grid cell

      gisend = mxcohrtdat.getdel( fnumchrts );

      if( -1 == gisend )
      {
        cout << "Ran out of Number of Cohorts data";
        cout << endl << endl;
        flog1 << "Ran out of Number of Cohorts data";
        flog1 << endl << endl;

        exit( -1 );
      }

      telmnt[igrd].col = mxcohrtdat.col;
      telmnt[igrd].row = mxcohrtdat.row;


      // Get land area information for the grid cell

      gisend = landareadat.getdel( flandarea );

      if( -1 == gisend ) 
      { 
        cout << "Ran out of Land area data";
        cout << endl << endl;

        flog1 << "Ran out of Land area data";
        flog1 << endl << endl;
        
        exit( -1 );
      } 

      // Check data for spatial coregistration errors

      fatalerr = telmnt[0].coregerr( flog1, 
                                     "MXCOHRTS", 
                                     telmnt[igrd].col,  
                                     telmnt[igrd].row, 
                                     "LANDAREA", 
                                     landareadat.col, 
                                     landareadat.row );

      if( fatalerr != 0 ) { exit( -1 ); }


      landcover4tem_.elmentArea[igrd] = landareadat.elemntArea;
      landcover4tem_.landFrac[igrd] = landareadat.landFrac;

      ++igrd;

    }
  }
 
  igrd = 0;

  while ( igrd < MAXNGRD && 0 == fatalerr )
  {
    // Get fraction of land area of each total number of cohorts in the grid cell

    if( (0 == telmnt[0].lcluc.tlulcflag && 0 == pdyr)
        || 1 == telmnt[0].lcluc.tlulcflag  )
    {
      if( pdyr > 0 )
      {         	              
        // Update prevchrtarea with area of cohort from 
        //  the previous year
         
        for( ichrt = 0; ichrt < telmnt[igrd].maxcohorts; ++ichrt )
        {
          telmnt[igrd].cohort[ichrt].prevchrtarea = telmnt[igrd].cohort[ichrt].chrtarea;
        }
      }

      // Get current land use/land cover cohort data for the grid cell 

      for( ichrt = 0; ichrt < mxcohrtdat.total; ++ichrt )
      {
        gisend = lulcdat[ichrt].getdel( flulc );

        if( -1 == gisend ) 
        { 
          flog1 << "Ran out of Land cover/land use data";
          flog1 << endl << endl;
            
          exit( -1 );
        }

        // Check data for spatial coregistration errors

        fatalerr = telmnt[0].coregerr( flog1, 
                                       "LANDAREA", 
                                       telmnt[igrd].col,  
                                       telmnt[igrd].row, 
                                       "LULC", 
                                       lulcdat[ichrt].col, 
                                       lulcdat[ichrt].row );

        if( fatalerr != 0 ) { exit( -1 ); }

        landcover4tem_.fracLandArea[igrd][ichrt] = lulcdat[ichrt].fracLandArea;
      }      
    }  
  
    ++igrd;
  }

};
#endif

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */

void RegMITTEM::updateMITTEMregion( const int& pdyr, 
                                    const int& pdm )
{
  const double gCtogCH4 = 16.0 / 12.0;

  int ichrt;
  long igrd;
  
  double outch4flx;

  int tchrt;

  time_t timer;


  icount = 0;

// *************************************************************
// Skip to the desired record in the GIS data sets

//  grdcnt = skipRecords();

//************************************************************ */


  // Save TEM state of December for every-other-year to a 
  //   "temstate" file
   
  if ( (CYCLE-1) == pdm )
  {
    if ( 1 == yrostateflag )
    {
      if ( 0 == pdyr%2 )
      {
        telmnt[0].ofstateA.open( telmnt[0].temstateAfname.c_str(), 
                                 ios::out );
      }
      else
      {
        telmnt[0].ofstateB.open( telmnt[0].temstateBfname.c_str(), 
                                 ios::out );
      }
    }
  }

  // Reset aggregated CFLUX (by MIT-IGSM latitudinal band) 
  //   to zero

  telmnt[0].mitclm.resetCFLUX( pdm );

  // Reset aggregated CH4FLUX (by MIT-IGSM latitudinal band) 
  //   to zero

  telmnt[0].mitclm.resetCH4FLUX( pdm );

  // Reset aggregated N2OFLUX (by MIT-IGSM latitudinal band) 
  //   to zero

  telmnt[0].mitclm.resetN2OFLUX( pdm );


  igrd = 0;  // reinitialize grid cell count for current year
  while ( igrd < mxnumgrid && 0 == fatalerr ) // Grid cell loop
  {
    telmnt[igrd].tem.setLAT( telmnt[igrd].row );


/* *************************************************************
           INITIALIZE TEM STATE FOR NEW COHORTS
************************************************************* */
  
    if( telmnt[igrd].maxcohorts > telmnt[igrd].prvmxcohrts )
    {
      for( ichrt = telmnt[igrd].prvmxcohrts; 
           ichrt < telmnt[igrd].maxcohorts; 
           ++ichrt )
      {     

        tchrt = telmnt[igrd].cohort[ichrt].srcCohort - 1;

        telmnt[igrd].setCohortTEMState( telmnt[igrd].cohort[tchrt],
                                        telmnt[igrd].cohort[ichrt] );
      }
    }

/* *************************************************************
           IMPLEMENT LAND-USE CHANGE
************************************************************* */
 

    if( pdyr > 0 && 0 == pdm )
    {
      // Determine area changes in land cover due to land-use 
      //   change
 
      setMITLULCgrid( pdyr, igrd );


      // Determine carbon and nitrogen fluxes due to land 
      //   conversions and readjust carbon and nitrogen pools 
      //   among managed and natural land cover cohorts based
      //   on area changes determined above
      
      telmnt[igrd].landUseChange( pdyr );
    }
    
/* *************************************************************
                     UPDATE TEM FOR GRID CELL
************************************************************* */

    
/* *************************************************************
		BEGIN VEGETATION MOSAIC LOOP
************************************************************* */

    for( ichrt = 0; ichrt < telmnt[igrd].maxcohorts; ++ichrt )
    {
      telmnt[igrd].tem.veg.cmnt = telmnt[igrd].cohort[ichrt].cmnt;


      // Redetermine soil characteristics (instead of saving it 
      //   to memory)
    
      telmnt[igrd].tem.soil.xtext( telmnt[igrd].tem.veg.cmnt, 
                                   telmnt[igrd].tem.soil.getPCTSILT(), 
                                   telmnt[igrd].tem.soil.getPCTCLAY() );

      
      // Obtain gridded climate from MIT-IGSM latitudinal bands

      setMITCLMgrid( pdyr, pdm, igrd, ichrt );

      if( istateflag > 0 && 1 == pdyr )
      {
        if( telmnt[igrd].cohort[ichrt].eetmx < telmnt[igrd].tem.soil.getEET() )
        {
          telmnt[igrd].cohort[ichrt].eetmx = telmnt[igrd].tem.soil.getEET();
        }

        if( telmnt[igrd].tem.atms.getPET() < telmnt[igrd].tem.soil.getEET() )
        {
          telmnt[igrd].tem.atms.setPET( telmnt[igrd].tem.soil.getEET() );
        }

        if( telmnt[igrd].cohort[ichrt].petmx < telmnt[igrd].tem.atms.getPET() )
        {
          telmnt[igrd].cohort[ichrt].petmx = telmnt[igrd].tem.atms.getPET();
        }

        if( CROPVEG == ichrt || BIOFUELS == ichrt || PASTURE == ichrt )
        { 
          // Initialize previous monthly CO2 to current monthly CO2 if 
          //   managed lands have not been initialized yet
 
          if( telmnt[igrd].tem.atms.getPREVCO2() < ZERO )
          {
            telmnt[igrd].tem.atms.setPREVCO2( telmnt[igrd].tem.atms.getCO2() );
          }
 

          // Determine vegetation C/N parameter as a function of
          //   vegetation type, annual PET, annual EET, CO2
          //   concentration

          telmnt[igrd].tem.veg.updateC2N( telmnt[igrd].cohort[ichrt].cmnt,
                                          telmnt[igrd].tem.soil.yreet,
                                          telmnt[igrd].tem.atms.yrpet,
                                          telmnt[igrd].tem.atms.getPREVCO2(),
                                          telmnt[igrd].tem.atms.getINITCO2() );

          telmnt[igrd].cohort[ichrt].cneven = telmnt[igrd].tem.veg.getCNEVEN();

/*          if( PASTURE == ichrt && telmnt[igrd].cohort[ichrt].chrtarea > 0 )
          {
            cout << " lat = " << telmnt[igrd].row << endl;
            cout << " cmnt = " << telmnt[igrd].cohort[ichrt].cmnt << endl;
            cout << " yreet = " << telmnt[igrd].tem.soil.yreet << endl;
            cout << " yrpet = " << telmnt[igrd].tem.atms.yrpet << endl;
            cout << " prevco2 = " << telmnt[igrd].tem.atms.getPREVCO2() << endl;
            cout << " initco2 = " << telmnt[igrd].tem.atms.getINITCO2() << endl;
            cout << " cneven = " <<  telmnt[igrd].cohort[ichrt].cneven << endl;
            exit( -1 );           
          }
*/
          telmnt[igrd].cohort[ichrt].productYear = telmnt[igrd].tem.ag.getPRODUCTYEAR();
        }
      }

      telmnt[igrd].updateTEMmonth( equil, 
                                   totsptime, 
                                   pdyr, 
                                   pdm,
                                   ichrt );

      if( (telmnt[igrd].tem.veg.getCURRENTVEG() > 0
          && telmnt[igrd].tem.veg.getCURRENTVEG() < 30)
          || 33 == telmnt[igrd].tem.veg.getCURRENTVEG() )
      { 
        telmnt[igrd].mitclm.aggregFlux2D( telmnt[igrd].tem.totyr,
                                          telmnt[igrd].row,
                                          (double) telmnt[igrd].cohort[ichrt].chrtarea,
                                          telmnt[igrd].tem.getNCE(),
                                          telmnt[0].mitclm.tcflx2D );
      }


      // Determine aggregated monthly methane fluxes from grid cell
      //   used by the MIT IGSM (i.e. 4 degree latitudinal bands)
      
      if( 1 == telmnt[igrd].tem.ch4flag 
          && ((telmnt[igrd].tem.veg.getCURRENTVEG() > 0
          && telmnt[igrd].tem.veg.getCURRENTVEG() < 30)
          || 33 == telmnt[igrd].tem.veg.getCURRENTVEG()) )
      { 
        outch4flx = telmnt[igrd].tem.soil.getCH4FLUX() 
                    * gCtogCH4
                    * 1000000000.0;

        telmnt[igrd].mitclm.aggregFlux2D( telmnt[igrd].tem.totyr, 
                                          telmnt[igrd].row, 
                                          (double) telmnt[igrd].cohort[ichrt].chrtarea, 
                                          outch4flx,
                                          telmnt[0].mitclm.tch4flx2D );
      } 


      // Determine aggregated monthly nitrous fluxes from grid cell
      //   used by the MIT IGSM (i.e. 4 degree latitudinal bands)
      //cout << " NONSOLC before aggregFlux in UpdateMITRegion = " << telmnt[igrd].cohort[ichrt].NEMnsolc << endl;

      if ( 1 == telmnt[igrd].tem.n2oflag
           && telmnt[igrd].tem.getY( telmnt[igrd].tem.I_SOLC ) > ZERO 
           && ((telmnt[igrd].tem.veg.getCURRENTVEG() > 0
           && telmnt[igrd].tem.veg.getCURRENTVEG() < 30)
           || 33 == telmnt[igrd].tem.veg.getCURRENTVEG()) )
      {
        telmnt[igrd].mitclm.aggregFlux2D( telmnt[igrd].tem.totyr, 
                                          telmnt[igrd].row, 
                                          (double) telmnt[igrd].cohort[ichrt].chrtarea, 
                                          (telmnt[igrd].tem.soil.getN2OFLUX() * 1000000000.0),
                                          telmnt[0].mitclm.tn2oflx2D );
      }
      
      if( pdm == (CYCLE-1) )
      {
        if( 2 == ostateflag && telmnt[igrd].tem.totyr == ostateyear )
        {
          telmnt[igrd].writeCohortState( ofstate, ichrt );
        }

        if( 1 == yrostateflag )
        {
          if ( 0 == pdyr%2 ) 
          { 
            telmnt[igrd].writeCohortState( telmnt[0].ofstateA, 
                                           ichrt );
          }
          else
          { 
            telmnt[igrd].writeCohortState( telmnt[0].ofstateB, 
                                           ichrt );
          }
          
        }

        if ( (1 == spinoutfg && telmnt[igrd].tem.totyr < telmnt[0].tem.startyr)
             || (2 == spinoutfg 
             && telmnt[igrd].tem.totyr >= (telmnt[0].tem.startyr-spinoutyrs))
             || (telmnt[igrd].tem.totyr >= telmnt[0].tem.startyr 
             && telmnt[igrd].tem.totyr <= telmnt[0].tem.endyr)
             && 0 == (telmnt[0].wrtyr%telmnt[0].tem.diffyr) )
        {


          // Output TEM transient results for specified years to files

          telmnt[igrd].temwritepred( ftempred, 
                                     tempredmap, 
                                     pdyr, 
                                     ichrt,
                                     telmnt[0].ntempred );

        }
      }
    } // End of ichrt loop
 

    if( (0 == istateflag && totsptime == pdyr && (CYCLE-1) == pdm)
        || istateflag > 0 )
    {
      telmnt[igrd].tem.microbe.nem.startflag = 1;
    }


    if( pdyr < 3 && (CYCLE-1) == pdm )
    {
      elmnt.show( flog1, 
                  telmnt[igrd].col, 
                  telmnt[igrd].row,
                  telmnt[igrd].tem.totyr, 
                  telmnt[igrd].tem.tol );
    }
   
    ++igrd;
  } // End of igrd loop


//  if( 1 == pdyr ) { exit( -1 ); }


  if( (CYCLE-1) == pdm )
  {
    if ( 1 == yrostateflag ) 
    {
      if ( 0 == pdyr%2 ) { telmnt[0].ofstateA.close(); }
      else { telmnt[0].ofstateB.close(); }
    }

    timer = time( NULL );
    flog1 << "Finished year " << telmnt[0].tem.totyr;
    flog1 << " at " << ctime( &timer );
  }
  
};

