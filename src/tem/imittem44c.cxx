/* *************************************************************
****************************************************************
IMITTEM44C.CXX - Interactive version of Terrestrial Ecosystem 
                  Model Version 4.4 for use with the MIT IGSM 
****************************************************************


Modifications:

20060202 - DWK created by modifying imittem436.cpp 
20060202 - DWK changed include from regmittem436.h to 
           regmittem437.h
20060202 - DWK replaced xtem.updateMITTEMregionYr() with 
           xtem.updateMITTEMregion()
20060204 - DWK deleted void assignInitMITCO2TEM()
20070129 - DWK changed include regmittem438.h to
           regmittem44.h
20080130 - DWK changed include from regmittem44.h to 
           regmittem44a.h
20080131 - DWK added xtem.telmnt[0].lcluc.initPotvegMaxCohorts()
           to initLCLUC()
20080131 - DWK added xtem.telmnt[0].lcluc.initPotvegCohorts()
           to initLCLUC()
20080826 - DWK changed include from regmittem44a.h to 
           regmittem44b.h
20110705 - DWK modified initLCLUC() and updateLCLUC() to read in 
           xtem.telmnt[0].lcluc.initPotvegMaxCohorts() only once
           for number of maximum number of cohorts, but still 
           allowed the LULCCOHRT file to vary over time
20110705 - DWK changed include from regmittem44b.h to 
           regmittem44c.h
20111019 - DWK added PET to assignMITCLM2TEM()
                                                                                            
****************************************************************
************************************************************* */

// NOTE:  To help avoid duplication of code, MITTEM uses 
//        compiler directives for different applications of the 
//        model or to handle compiler-specific quirks.  Below, 
//        we identify these compiler directives and indicate
//        where in the code they are defined.

// NOTE:  The Borland C++ Builder compiler stores time variables
//        in a different library <time> than the g++ or
//        Portland Group C++ complier <ctime>.  When using the 
//        Borland compiler DEFINE the compiler directive 
//        BORLAND_CPP below. Otherwise, DEFINE the compiler 
//        directive ANSI_CPP below.  In addition, these compiler 
//        directives need to be consistently defined in
//        regmittem437.cpp and celmnt437.cpp  

//#define BORLAND_CPP
#define ANSI_CPP


// NOTE: imittem44a.cxx is the outer shell of the "offline" 
//       version of MITTEM.  Make sure the compiler directive 
//       STANDALONE_TEM is DEFINED below and in regmittem437.h,
//       regmittem437.cpp, climate4tem437.hpp, mitclm437.h 
//       and mitclm437.cpp

#define STANDALONE_TEM


// NOTE: mittem44.cpp contains DOS WINDOWS-specific code used  
//       to calibrate TEM.  This code is not used when 
//       extrapolating the model with spatially explicit input
//       data sets.  To avoid using this code, make sure the 
//       compiler directive CALIBRATE_TEM in mittem437.cpp is 
//       be NOT DEFINED (i.e. commented out)



#include<cstdio>

  using std::fopen;
  using std::fclose;
  using std::printf;
  using std::sprintf;
  using std::fscanf;
  using std::FILE;

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

  using std::setiosflags;
  using std::setw;
  using std::setprecision;

#include<cstdlib>

  using std::exit;
  
#include<cmath>

#include<vector>

  using std::vector;
  
#include<cctype>

  using std::toupper;

#include<cstring>
  
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


#include "regmittem44c.h"    // RegMITTEM Class

void assignMITCLM2TEM( const int& pdyr, 
                       const int& pdm, 
                       RegMITTEM& xtem );
                       
void assignTEM2MITCLM( const int& pdyr, 
                       const int& pdm,
                       RegMITTEM& xtem );

void initLCLUC( RegMITTEM& xtem );
void initMITCLM( RegMITTEM& xtem );
void initTEM( RegMITTEM& xtem );

void updateLCLUC( int& pdyr, RegMITTEM& xtem );


RegMITTEM xtem;


/* *************************************************************
**********************START MAIN PROGRAM************************
************************************************************* */

int main()
{
  int dm;
  int dyr;


  // Run model in steady state mode or dynamic mode? 

  xtem.setRunMode();


  // Specify if TEM is to be started from calibration data
  //   or from spatially explicit conditions determined from
  //   a previous run

  xtem.setInTEMState();


  // Determine time frame of simulation
  
  xtem.setRunTime();

  
  // Specify number of grid cells in the study region and 
  //   determine if spatially explicit data is referenced by 
  //   column and row or by longitude and latitude
  
  xtem.setGridParameters();


  // Initialize climate module for study region

  initMITCLM( xtem );


  // Initialize land cover module for study region

  initLCLUC( xtem );

  
  // Initialize TEM for study region
  
  initTEM( xtem );


  // Begin simulation of transient climate and terrestrial 
  //   ecosystem response

  if( 0 == xtem.equil )
  {
    for( dyr = 1; dyr < xtem.RTIME; ++dyr )
    {
      // Run land cover module or read in land cover data from 
      //   file to update land cover characteristics for study 
      //   region

      updateLCLUC( dyr, xtem ); 


      // Run TEM for study region 
    
      for ( dm = 0; dm < CYCLE; ++dm )
      {
        // Run climate module or read climate data from file to
        //   update climate for study region

        assignMITCLM2TEM( dyr, dm, xtem );


        // Run TEM for study region to obtain monthly estimates
       
        xtem.updateMITTEMregion( dyr, dm );
      

        if ( 2 == xtem.istateflag 
             || (xtem.istateflag < 2 && dyr > xtem.totsptime) ) 
        {  

          // Assign NEP (CFLUX) determined by TEM to a climate 
          //   model variable

          assignTEM2MITCLM( dyr, dm, xtem );
        }
      }
    }
  }
  
  if ( 0 == xtem.fatalerr )
  {
    cout << "Extrapolation successfully completed -";
    cout << " Congratulations!" << endl;

    xtem.flog1 << "Extrapolation successfully completed -";
    xtem.flog1 << " Congratulations!" << endl;
  }

  return 0;

};

/* *************************************************************
******************** End of Main ******************************* 
************************************************************* */

/* *************************************************************
************************************************************** */

void assignMITCLM2TEM( const int& pdyr, 
                       const int& pdm,
                       RegMITTEM& xtem )
{
  //int tdm;

  ostringstream tempfname;
  string nirrname;
  string tairname;
  string dayTsoilname;
  string dayTairname;
  string rainDurname;
  string rainIntname;
  string daySMname;
  string hrSMname;
  string precname;
  string petname;
  string eetname;
  string sh2o1name;
  string sh2o2name;
  string spackname;
  string srfrunname;
  string drainname;
  string co2name;
  string o3name; 

  static ifstream ifmitnirr;
  static ifstream ifmittair;

  static ifstream ifmitdayTsoilL1;
  static ifstream ifmitdayTsoilL2;
  static ifstream ifmitdayTsoilL3;
  static ifstream ifmitdayTsoilL4;
  static ifstream ifmitdayTsoilL5;
  static ifstream ifmitdayTsoilL6;
  static ifstream ifmitdayTsoilL7;
  static ifstream ifmitdayTsoilL8;
  static ifstream ifmitdayTsoilL9;

  static ifstream ifmitdayTair;
  static ifstream ifmitrainDur;
  static ifstream ifmitrainInt; 

  static ifstream ifmitdaySML1; 
  static ifstream ifmitdaySML2; 
  static ifstream ifmitdaySML3; 
  static ifstream ifmitdaySML4; 
  static ifstream ifmitdaySML5; 
  static ifstream ifmitdaySML6; 
  static ifstream ifmitdaySML7; 
  static ifstream ifmitdaySML8; 
  static ifstream ifmitdaySML9; 

  static ifstream ifmithrSML1;   
  static ifstream ifmithrSML2;   
  static ifstream ifmithrSML3;   
  static ifstream ifmithrSML4;   
  static ifstream ifmithrSML5;   
  static ifstream ifmithrSML6;   

  static ifstream ifmitprec;
  static ifstream ifmitpet;
  static ifstream ifmiteet;
  static ifstream ifmitsh2o1;
  static ifstream ifmitsh2o2;

  static ifstream ifmitspack;
  static ifstream ifmitsrfrun;
  static ifstream ifmitdrain;

  static ifstream ifmitco2;
  static ifstream ifmito3;

 
  int co2tstyr;
  int tstyr;

  if ( 0 == pdm )
  {
    // Read in MIT-IGSM climate data sets for the transient 
    //   spin-up portion of TEM initialization


    if( xtem.istateflag < 2 && pdyr < (xtem.totsptime+1) )
    {
      tstyr = (pdyr-1)%xtem.spintime;
    }
    else
    {
      // Read in MIT-IGSM transient climate data sets for the TEM
      // transient simulation

      if( xtem.istateflag < 2 ) 
      { 
        tstyr = pdyr - xtem.totsptime - 1; 
      }
      else { tstyr = pdyr - 1; }
    }
 
    // Open MIT cloudiness or solar radiation file

    if( 1 == xtem.telmnt[0].mitclm.tcldsflag )
    {
      tempfname.str( "" );
      tempfname << xtem.telmnt[0].mitclm.inirrfname
                << (xtem.telmnt[0].mitclm.startyr+tstyr)
                << xtem.telmnt[0].mitclm.inirrend;
        
      nirrname = tempfname.str();

      ifmitnirr.open( nirrname.c_str(), ios::in );

      if( !ifmitnirr )
      {
        xtem.flog1 << endl;
        xtem.flog1 << "Cannot open " << nirrname;
        xtem.flog1 << " for NIRR data input" << endl << endl;
          
        exit( -1 );
      }
    }
    else 
    { 
      ifmitnirr.open( xtem.telmnt[0].mitclm.inirrfname.c_str(), 
                      ios::in );

      if( !ifmitnirr )
      {
        xtem.flog1 << endl << "Cannot open ";
        xtem.flog1 << xtem.telmnt[0].mitclm.inirrfname;
        xtem.flog1 << " for NIRR data input" << endl << endl;
          
        exit( -1 );
      }
    }


    // Open MIT air temperature file
    if( 1 == xtem.telmnt[0].mitclm.ttairflag )
    {
      tempfname.str( "" );
      tempfname << xtem.telmnt[0].mitclm.itairfname
                << (xtem.telmnt[0].mitclm.startyr+tstyr)
                << xtem.telmnt[0].mitclm.itairend;
      
      tairname = tempfname.str();
      ifmittair.open( tairname.c_str(), ios::in );

      if( !ifmittair )
      {
        xtem.flog1 << endl;
        xtem.flog1 << "Cannot open " << tairname;
        xtem.flog1 << " for TAIR data input" << endl << endl;

        exit( -1 );
      }
    }
    else
    {
      ifmittair.open( xtem.telmnt[0].mitclm.itairfname.c_str(), 
                      ios::in );

      if( !ifmittair )
      {
        xtem.flog1 << endl << "Cannot open ";
        xtem.flog1 << xtem.telmnt[0].mitclm.itairfname;
        xtem.flog1 << " for TAIR data input" << endl << endl;

        exit( -1 );
      }
    }

    if( 1 == xtem.telmnt[0].mitclm.ch4flag
        || 1 == xtem.telmnt[0].mitclm.n2oflag )
    {
      if( 1 == xtem.telmnt[0].mitclm.ttairflag )
      {
        // Open MIT daily layer 1 soil temperature file
        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.idayTsoilfnameL1
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                  << xtem.telmnt[0].mitclm.idayTsoilendL1;
      
        dayTsoilname = tempfname.str();

        ifmitdayTsoilL1.open( dayTsoilname.c_str(), ios::in );

        if( !ifmitdayTsoilL1 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << dayTsoilname;
          xtem.flog1 << " for layer 1 DAYTSOIL data input";
          xtem.flog1 << endl << endl;

          exit( -1 );
        }
      

        // Open MIT daily layer 2 soil temperature file
        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.idayTsoilfnameL2
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                  << xtem.telmnt[0].mitclm.idayTsoilendL2;
      
        dayTsoilname = tempfname.str();

        ifmitdayTsoilL2.open( dayTsoilname.c_str(), ios::in );

        if( !ifmitdayTsoilL2 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << dayTsoilname;
          xtem.flog1 << " for layer 2 DAYTSOIL data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }
 

        // Open MIT daily layer 3 soil temperature file
        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.idayTsoilfnameL3
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                << xtem.telmnt[0].mitclm.idayTsoilendL3;
      
        dayTsoilname = tempfname.str();

        ifmitdayTsoilL3.open( dayTsoilname.c_str(), ios::in );

        if( !ifmitdayTsoilL3 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << dayTsoilname;
          xtem.flog1 << " for layer 3 DAYTSOIL data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }


        // Open MIT daily layer 4 soil temperature file
        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.idayTsoilfnameL4
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                  << xtem.telmnt[0].mitclm.idayTsoilendL4;
      
        dayTsoilname = tempfname.str();

        ifmitdayTsoilL4.open( dayTsoilname.c_str(), ios::in );

        if( !ifmitdayTsoilL4 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << dayTsoilname;
          xtem.flog1 << " for layer 4 DAYTSOIL data input";
          xtem.flog1 << endl << endl;
       
          exit( -1 );
        }


        // Open MIT daily layer 5 soil temperature file
        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.idayTsoilfnameL5
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                  << xtem.telmnt[0].mitclm.idayTsoilendL5;
    
        dayTsoilname = tempfname.str();

        ifmitdayTsoilL5.open( dayTsoilname.c_str(), ios::in );

        if( !ifmitdayTsoilL5 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << dayTsoilname;
          xtem.flog1 << " for layer 5 DAYTSOIL data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }


        // Open MIT daily layer 6 soil temperature file
        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.idayTsoilfnameL6
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                  << xtem.telmnt[0].mitclm.idayTsoilendL6;
      
        dayTsoilname = tempfname.str();

        ifmitdayTsoilL6.open( dayTsoilname.c_str(), ios::in );

        if( !ifmitdayTsoilL6 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << dayTsoilname;
          xtem.flog1 << " for layer 6 DAYTSOIL data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }


        // Open MIT daily layer 7 soil temperature file
        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.idayTsoilfnameL7
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                  << xtem.telmnt[0].mitclm.idayTsoilendL7;
      
        dayTsoilname = tempfname.str();

        ifmitdayTsoilL7.open( dayTsoilname.c_str(), ios::in );

        if( !ifmitdayTsoilL7 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << dayTsoilname;
          xtem.flog1 << " for layer 7 DAYTSOIL data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }


        // Open MIT daily layer 8 soil temperature file
        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.idayTsoilfnameL8
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                  << xtem.telmnt[0].mitclm.idayTsoilendL8;
      
        dayTsoilname = tempfname.str();

        ifmitdayTsoilL8.open( dayTsoilname.c_str(), ios::in );

        if( !ifmitdayTsoilL8 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << dayTsoilname;
          xtem.flog1 << " for layer 8 DAYTSOIL data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }


        // Open MIT daily layer 9 soil temperature file
        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.idayTsoilfnameL9
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                  << xtem.telmnt[0].mitclm.idayTsoilendL9;
      
        dayTsoilname = tempfname.str();

        ifmitdayTsoilL9.open( dayTsoilname.c_str(), ios::in );

        if( !ifmitdayTsoilL9 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << dayTsoilname;
          xtem.flog1 << " for layer 9 DAYTSOIL data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }
      }
      else
      {
        // Open MIT daily layer 1 soil temperature file
        ifmitdayTsoilL1.open( xtem.telmnt[0].mitclm.idayTsoilfnameL1.c_str(), 
                              ios::in );

        if( !ifmitdayTsoilL1 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.idayTsoilfnameL1;
          xtem.flog1 << " for layer 1 DAYTSOIL data input";
          xtem.flog1 << endl << endl;

          exit( -1 );
        }


        // Open MIT daily layer 2 soil temperature file
        ifmitdayTsoilL2.open( xtem.telmnt[0].mitclm.idayTsoilfnameL2.c_str(), 
                              ios::in );

        if( !ifmitdayTsoilL2 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.idayTsoilfnameL2;
          xtem.flog1 << " for layer 2 DAYTSOIL data input";
          xtem.flog1 << endl << endl;

          exit( -1 );
        }


        // Open MIT daily layer 3 soil temperature file
        ifmitdayTsoilL3.open( xtem.telmnt[0].mitclm.idayTsoilfnameL3.c_str(), 
                              ios::in );

        if( !ifmitdayTsoilL3 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.idayTsoilfnameL3;
          xtem.flog1 << " for layer 3 DAYTSOIL data input";
          xtem.flog1 << endl << endl;

          exit( -1 );
        }


        // Open MIT daily layer 4 soil temperature file
        ifmitdayTsoilL4.open( xtem.telmnt[0].mitclm.idayTsoilfnameL4.c_str(), 
                              ios::in );

        if( !ifmitdayTsoilL4 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.idayTsoilfnameL4;
          xtem.flog1 << " for layer 4 DAYTSOIL data input";
          xtem.flog1 << endl << endl;

          exit( -1 );
        }


        // Open MIT daily layer 5 soil temperature file
        ifmitdayTsoilL5.open( xtem.telmnt[0].mitclm.idayTsoilfnameL5.c_str(), 
                              ios::in );

        if( !ifmitdayTsoilL5 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.idayTsoilfnameL5;
          xtem.flog1 << " for layer 5 DAYTSOIL data input";
          xtem.flog1 << endl << endl;

          exit( -1 );
        }


        // Open MIT daily layer 6 soil temperature file
        ifmitdayTsoilL6.open( xtem.telmnt[0].mitclm.idayTsoilfnameL6.c_str(), 
                              ios::in );

        if( !ifmitdayTsoilL6 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.idayTsoilfnameL6;
          xtem.flog1 << " for layer 6 DAYTSOIL data input";
          xtem.flog1 << endl << endl;

          exit( -1 );
        }


        // Open MIT daily layer 7 soil temperature file
        ifmitdayTsoilL7.open( xtem.telmnt[0].mitclm.idayTsoilfnameL7.c_str(), 
                              ios::in );

        if( !ifmitdayTsoilL7 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.idayTsoilfnameL7;
          xtem.flog1 << " for layer 7 DAYTSOIL data input";
          xtem.flog1 << endl << endl;

          exit( -1 );
        }


        // Open MIT daily layer 8 soil temperature file
        ifmitdayTsoilL8.open( xtem.telmnt[0].mitclm.idayTsoilfnameL8.c_str(), 
                              ios::in );

        if( !ifmitdayTsoilL8 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.idayTsoilfnameL8;
          xtem.flog1 << " for layer 8 DAYTSOIL data input";
          xtem.flog1 << endl << endl;

          exit( -1 );
        }


        // Open MIT daily layer 9 soil temperature file
        ifmitdayTsoilL9.open( xtem.telmnt[0].mitclm.idayTsoilfnameL9.c_str(), 
                              ios::in );

        if( !ifmitdayTsoilL9 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.idayTsoilfnameL9;
          xtem.flog1 << " for layer 9 DAYTSOIL data input";
          xtem.flog1 << endl << endl;

          exit( -1 );
        }
      }

      if( 1 == xtem.telmnt[0].mitclm.tprecflag )
      {
        // Open MIT daily layer 1 soil moisture file
        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.idaySMfnameL1
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                  << xtem.telmnt[0].mitclm.idaySMendL1;
     
        daySMname = tempfname.str();

        ifmitdaySML1.open( daySMname.c_str(), ios::in );

        if( !ifmitdaySML1 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << daySMname;
          xtem.flog1 << " for daily layer 1 soil moisture data input";
          xtem.flog1 << endl << endl;
       
          exit( -1 );
        }

        // Open MIT daily layer 2 soil moisture file
        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.idaySMfnameL2
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                  << xtem.telmnt[0].mitclm.idaySMendL2;
        
        daySMname = tempfname.str();

        ifmitdaySML2.open( daySMname.c_str(), ios::in );

        if( !ifmitdaySML2 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << daySMname;
          xtem.flog1 << " for daily layer 2 soil moisture data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }

        // Open MIT daily layer 3 soil moisture file
        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.idaySMfnameL3
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                  << xtem.telmnt[0].mitclm.idaySMendL3;
      
        daySMname = tempfname.str();

        ifmitdaySML3.open( daySMname.c_str(), ios::in );

        if( !ifmitdaySML3 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << daySMname;
          xtem.flog1 << " for daily layer 3 soil moisture data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }

        // Open MIT daily layer 4 soil moisture file
        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.idaySMfnameL4
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                  << xtem.telmnt[0].mitclm.idaySMendL4;
      
        daySMname = tempfname.str();

        ifmitdaySML4.open( daySMname.c_str(), ios::in );

        if( !ifmitdaySML4 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << daySMname;
          xtem.flog1 << " for daily layer 4 soil moisture data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }

        // Open MIT daily layer 5 soil moisture file
        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.idaySMfnameL5
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                  << xtem.telmnt[0].mitclm.idaySMendL5;
      
        daySMname = tempfname.str();

        ifmitdaySML5.open( daySMname.c_str(), ios::in );

        if( !ifmitdaySML5 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << daySMname;
          xtem.flog1 << " for daily layer 5 soil moisture data input";
          xtem.flog1 << endl << endl;
      
          exit( -1 );
        }

        // Open MIT daily layer 6 soil moisture file
        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.idaySMfnameL6
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                  << xtem.telmnt[0].mitclm.idaySMendL6;
      
        daySMname = tempfname.str();

        ifmitdaySML6.open( daySMname.c_str(), ios::in );

        if( !ifmitdaySML6 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << daySMname;
          xtem.flog1 << " for daily layer 6 soil moisture data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }


        // Open MIT daily layer 7 soil moisture file
        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.idaySMfnameL7
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                  << xtem.telmnt[0].mitclm.idaySMendL7;
      
        daySMname = tempfname.str();

        ifmitdaySML7.open( daySMname.c_str(), ios::in );

        if( !ifmitdaySML7 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << daySMname;
          xtem.flog1 << " for daily layer 7 soil moisture data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }


        // Open MIT daily layer 8 soil moisture file

        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.idaySMfnameL8
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                  << xtem.telmnt[0].mitclm.idaySMendL8;
      
        daySMname = tempfname.str();

        ifmitdaySML8.open( daySMname.c_str(), ios::in );

        if( !ifmitdaySML8 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << daySMname;
          xtem.flog1 << " for daily layer 8 soil moisture data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }


        // Open MIT daily layer 9 soil moisture file

        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.idaySMfnameL9
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                  << xtem.telmnt[0].mitclm.idaySMendL9;
      
        daySMname = tempfname.str();

        ifmitdaySML9.open( daySMname.c_str(), ios::in );

        if( !ifmitdaySML9 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << daySMname;
          xtem.flog1 << " for daily layer 9 soil moisture data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }
      }
      else
      {
        ifmitdaySML1.open( xtem.telmnt[0].mitclm.idaySMfnameL1.c_str(), 
                           ios::in );

        if( !ifmitdaySML1 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.idaySMfnameL1;
          xtem.flog1 << " for daily layer 1 soil moisture data input";
          xtem.flog1 << endl << endl;
       
          exit( -1 );
        }


        ifmitdaySML2.open( xtem.telmnt[0].mitclm.idaySMfnameL2.c_str(), 
                           ios::in );

        if( !ifmitdaySML2 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.idaySMfnameL2;
          xtem.flog1 << " for daily layer 2 soil moisture data input";
          xtem.flog1 << endl << endl;
       
          exit( -1 );
        }


        ifmitdaySML3.open( xtem.telmnt[0].mitclm.idaySMfnameL3.c_str(), 
                           ios::in );

        if( !ifmitdaySML3 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.idaySMfnameL3;
          xtem.flog1 << " for daily layer 3 soil moisture data input";
          xtem.flog1 << endl << endl;
       
          exit( -1 );
        }


        ifmitdaySML4.open( xtem.telmnt[0].mitclm.idaySMfnameL4.c_str(), 
                           ios::in );

        if( !ifmitdaySML4 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.idaySMfnameL4;
          xtem.flog1 << " for daily layer 4 soil moisture data input";
          xtem.flog1 << endl << endl;
       
          exit( -1 );
        }


        ifmitdaySML5.open( xtem.telmnt[0].mitclm.idaySMfnameL5.c_str(), 
                           ios::in );

        if( !ifmitdaySML5 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.idaySMfnameL5;
          xtem.flog1 << " for daily layer 5 soil moisture data input";
          xtem.flog1 << endl << endl;
       
          exit( -1 );
        }


        ifmitdaySML6.open( xtem.telmnt[0].mitclm.idaySMfnameL6.c_str(), 
                           ios::in );

        if( !ifmitdaySML6 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.idaySMfnameL6;
          xtem.flog1 << " for daily layer 6 soil moisture data input";
          xtem.flog1 << endl << endl;
       
          exit( -1 );
        }


        ifmitdaySML7.open( xtem.telmnt[0].mitclm.idaySMfnameL7.c_str(), 
                           ios::in );

        if( !ifmitdaySML7 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.idaySMfnameL7;
          xtem.flog1 << " for daily layer 7 soil moisture data input";
          xtem.flog1 << endl << endl;
       
          exit( -1 );
        }


        ifmitdaySML8.open( xtem.telmnt[0].mitclm.idaySMfnameL8.c_str(), 
                           ios::in );

        if( !ifmitdaySML8 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.idaySMfnameL8;
          xtem.flog1 << " for daily layer 8 soil moisture data input";
          xtem.flog1 << endl << endl;
       
          exit( -1 );
        }


        ifmitdaySML9.open( xtem.telmnt[0].mitclm.idaySMfnameL9.c_str(), 
                           ios::in );

        if( !ifmitdaySML9 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.idaySMfnameL9;
          xtem.flog1 << " for daily layer 9 soil moisture data input";
          xtem.flog1 << endl << endl;
       
          exit( -1 );
        }
      }
    }


    if( 1 == xtem.telmnt[0].mitclm.n2oflag )
    {
      if( 1 == xtem.telmnt[0].mitclm.ttairflag )
      {
        // Open MIT daily air temperature file
        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.idayTairfname
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                  << xtem.telmnt[0].mitclm.idayTairend;
      
        dayTairname = tempfname.str();

        ifmitdayTair.open( dayTairname.c_str(), ios::in );

        if( !ifmitdayTair )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << dayTairname;
          xtem.flog1 << " for DAYTAIR data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }
      }
      else
      {
        ifmitdayTair.open( xtem.telmnt[0].mitclm.idayTairfname.c_str(), 
                           ios::in );

        if( !ifmitdayTair )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.idayTairfname;
          xtem.flog1 << " for DAYTAIR data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }
      }

      if( 1 == xtem.telmnt[0].mitclm.tprecflag )
      {
        // Open MIT daily rain duration file
        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.irainDurfname
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                  << xtem.telmnt[0].mitclm.irainDurend;
     
        rainDurname = tempfname.str();

        ifmitrainDur.open( rainDurname.c_str(), ios::in );

        if( !ifmitrainDur )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << rainDurname;
          xtem.flog1 << " for rain duration data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }

        // Open MIT daily rain intensity file
        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.irainIntfname
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                  << xtem.telmnt[0].mitclm.irainIntend;
      
        rainIntname = tempfname.str();

        ifmitrainInt.open( rainIntname.c_str(), ios::in );

        if( !ifmitrainInt )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << rainIntname;
          xtem.flog1 << " for rain intensity data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }

        // Open MIT hourly layer 1 soil moisture file
        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.ihrSMfnameL1
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                  << xtem.telmnt[0].mitclm.ihrSMendL1;
      
        hrSMname = tempfname.str();

        ifmithrSML1.open( hrSMname.c_str(), ios::in );

        if( !ifmithrSML1 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << hrSMname;
          xtem.flog1 << " for hourly layer 1 soil moisture data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }

        // Open MIT hourly layer 2 soil moisture file
        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.ihrSMfnameL2
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                  << xtem.telmnt[0].mitclm.ihrSMendL2;
      
        hrSMname = tempfname.str();

        ifmithrSML2.open( hrSMname.c_str(), ios::in );

        if( !ifmithrSML2 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << hrSMname;
          xtem.flog1 << " for hourly layer 2 soil moisture data input";
          xtem.flog1 << endl << endl;
       
          exit( -1 );
        }

        // Open MIT hourly layer 3 soil moisture file
        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.ihrSMfnameL3
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                  << xtem.telmnt[0].mitclm.ihrSMendL3;
      
        hrSMname = tempfname.str();

        ifmithrSML3.open( hrSMname.c_str(), ios::in );

        if( !ifmithrSML3 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << hrSMname;
          xtem.flog1 << " for hourly layer 3 soil moisture data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }

        // Open MIT hourly layer 4 soil moisture file
        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.ihrSMfnameL4
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                  << xtem.telmnt[0].mitclm.ihrSMendL4;
      
        hrSMname = tempfname.str();

        ifmithrSML4.open( hrSMname.c_str(), ios::in );

        if( !ifmithrSML4 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << hrSMname;
          xtem.flog1 << " for hourly layer 4 soil moisture data input";
          xtem.flog1 << endl << endl;
      
          exit( -1 );
        }

        // Open MIT hourly layer 5 soil moisture file
        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.ihrSMfnameL5
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                  << xtem.telmnt[0].mitclm.ihrSMendL5;
     
        hrSMname = tempfname.str();

        ifmithrSML5.open( hrSMname.c_str(), ios::in );

        if( !ifmithrSML5 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << hrSMname;
          xtem.flog1 << " for hourly layer 5 soil moisture data input";
          xtem.flog1 << endl << endl;
       
          exit( -1 );
        }

        // Open MIT hourly layer 6 soil moisture file
        tempfname.str( "" );
        tempfname << xtem.telmnt[0].mitclm.ihrSMfnameL6
                  << (xtem.telmnt[0].mitclm.startyr+tstyr)
                  << xtem.telmnt[0].mitclm.ihrSMendL6;
     
        hrSMname = tempfname.str();

        ifmithrSML6.open( hrSMname.c_str(), ios::in );

        if( !ifmithrSML6 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open " << hrSMname;
          xtem.flog1 << " for hourly layer 6 soil moisture data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }
      }
      else
      {
        ifmitrainDur.open( xtem.telmnt[0].mitclm.irainDurfname.c_str(), 
                           ios::in );

        if( !ifmitrainDur )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.irainDurfname;
          xtem.flog1 << " for rain duration data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }


        ifmitrainInt.open( xtem.telmnt[0].mitclm.irainIntfname.c_str(), 
                           ios::in );

        if( !ifmitrainInt )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.irainIntfname;
          xtem.flog1 << " for rain intensity data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }


        ifmithrSML1.open( xtem.telmnt[0].mitclm.ihrSMfnameL1.c_str(), 
                          ios::in );

        if( !ifmithrSML1 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.ihrSMfnameL1;
          xtem.flog1 << " for hourly layer 1 soil moisture data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }


        ifmithrSML2.open( xtem.telmnt[0].mitclm.ihrSMfnameL2.c_str(), 
                          ios::in );

        if( !ifmithrSML2 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.ihrSMfnameL2;
          xtem.flog1 << " for hourly layer 2 soil moisture data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }


        ifmithrSML3.open( xtem.telmnt[0].mitclm.ihrSMfnameL3.c_str(), 
                          ios::in );

        if( !ifmithrSML3 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.ihrSMfnameL3;
          xtem.flog1 << " for hourly layer 3 soil moisture data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }


        ifmithrSML4.open( xtem.telmnt[0].mitclm.ihrSMfnameL4.c_str(), 
                          ios::in );

        if( !ifmithrSML4 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.ihrSMfnameL4;
          xtem.flog1 << " for hourly layer 4 soil moisture data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }


        ifmithrSML5.open( xtem.telmnt[0].mitclm.ihrSMfnameL5.c_str(), 
                          ios::in );

        if( !ifmithrSML5 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.ihrSMfnameL5;
          xtem.flog1 << " for hourly layer 5 soil moisture data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }


        ifmithrSML6.open( xtem.telmnt[0].mitclm.ihrSMfnameL6.c_str(), 
                          ios::in );

        if( !ifmithrSML6 )
        {
          xtem.flog1 << endl;
          xtem.flog1 << "Cannot open ";
          xtem.flog1 << xtem.telmnt[0].mitclm.ihrSMfnameL6;
          xtem.flog1 << " for hourly layer 6 soil moisture data input";
          xtem.flog1 << endl << endl;
        
          exit( -1 );
        }
      }
    }
    
    
    if( 1 == xtem.telmnt[0].mitclm.tprecflag )
    {
      // Open MIT estimated precipitation file
  
      tempfname.str( "" );
      tempfname << xtem.telmnt[0].mitclm.iprecfname
                << (xtem.telmnt[0].mitclm.startyr+tstyr)
                << xtem.telmnt[0].mitclm.iprecend;
      
      precname = tempfname.str();

      ifmitprec.open( precname.c_str(), ios::in );

      if( !ifmitprec )
      {
        xtem.flog1 << endl << "Cannot open " << precname;
        xtem.flog1 << " for PREC data input" << endl << endl;
        
        exit( -1 );
      }


      // Open MIT potential evapotranspiration file
  
      tempfname.str( "" );
      tempfname << xtem.telmnt[0].mitclm.ipetfname
                << (xtem.telmnt[0].mitclm.startyr+tstyr)
                << xtem.telmnt[0].mitclm.ipetend;
      
      petname = tempfname.str();

      ifmitpet.open( petname.c_str(), ios::in );

      if( !ifmitpet )
      {
        xtem.flog1 << endl << "Cannot open " << petname;
        xtem.flog1 << " for PET data input" << endl << endl;
        
        exit( -1 );
      }

      // Open MIT estimated actual evapotranspiration file
  
      tempfname.str( "" );
      tempfname << xtem.telmnt[0].mitclm.ieetfname
                << (xtem.telmnt[0].mitclm.startyr+tstyr)
                << xtem.telmnt[0].mitclm.ieetend;
      
      eetname = tempfname.str();

      ifmiteet.open( eetname.c_str(), ios::in );

      if( !ifmiteet )
      {
        xtem.flog1 << endl << "Cannot open " << eetname;
        xtem.flog1 << " for EET data input" << endl << endl;
        
        exit( -1 );
      }

      // Open MIT monthly soil moisture file for the top 1 meter
  
      tempfname.str( "" );
      tempfname << xtem.telmnt[0].mitclm.ish2o1fname
                << (xtem.telmnt[0].mitclm.startyr+tstyr)
                << xtem.telmnt[0].mitclm.ish2o1end;
     
      sh2o1name = tempfname.str();

      ifmitsh2o1.open( sh2o1name.c_str(), ios::in );

      if( !ifmitsh2o1 )
      {
        xtem.flog1 << endl << "Cannot open " << sh2o1name;
        xtem.flog1 << " for SOILH2O data input for top  1 meter";
        xtem.flog1 << endl << endl;
       
        exit( -1 );
      }

      // Open MIT monthly soil moisture file for the top 2 meters
  
      tempfname.str( "" );
      tempfname << xtem.telmnt[0].mitclm.ish2o2fname
                << (xtem.telmnt[0].mitclm.startyr+tstyr)
                << xtem.telmnt[0].mitclm.ish2o2end;
      
      sh2o2name = tempfname.str();

      ifmitsh2o2.open( sh2o2name.c_str(), ios::in );

      if( !ifmitsh2o2 )
      {
        xtem.flog1 << endl << "Cannot open " << sh2o2name;
        xtem.flog1 << " for SOILH2O data input for top  2 meters";
        xtem.flog1 << endl << endl;
        
        exit( -1 );
      }

      // Open MIT monthly snowpack
  
      tempfname.str( "" );
      tempfname << xtem.telmnt[0].mitclm.ispackfname
                << (xtem.telmnt[0].mitclm.startyr+tstyr)
                << xtem.telmnt[0].mitclm.ispackend;
      
      spackname = tempfname.str();

      ifmitspack.open( spackname.c_str(), ios::in );

      if( !ifmitspack )
      {
        xtem.flog1 << endl << "Cannot open " << spackname;
        xtem.flog1 << " for SNOWPACK data";
        xtem.flog1 << endl << endl;
        
        exit( -1 );
      }

      // Open MIT monthly surface runoff file
  
      tempfname.str( "" );
      tempfname << xtem.telmnt[0].mitclm.isrfrunfname
                << (xtem.telmnt[0].mitclm.startyr+tstyr)
                << xtem.telmnt[0].mitclm.isrfrunend;
      
      srfrunname = tempfname.str();

      ifmitsrfrun.open( srfrunname.c_str(), ios::in );

      if( !ifmitsrfrun )
      {
        xtem.flog1 << endl << "Cannot open " << srfrunname;
        xtem.flog1 << " for SURFRUN data";
        xtem.flog1 << endl << endl;
        
        exit( -1 );
      }

      // Open MIT monthly drainage file
  
      tempfname.str( "" );
      tempfname << xtem.telmnt[0].mitclm.idrainfname
                << (xtem.telmnt[0].mitclm.startyr+tstyr)
                << xtem.telmnt[0].mitclm.idrainend;
      
      drainname = tempfname.str();

      ifmitdrain.open( drainname.c_str(), ios::in );

      if( !ifmitdrain )
      {
        xtem.flog1 << endl << "Cannot open " << drainname;
        xtem.flog1 << " for DRAINAGE data";
        xtem.flog1 << endl << endl;
        
        exit( -1 );
      }
    }
    else
    {
      ifmitprec.open( xtem.telmnt[0].mitclm.iprecfname.c_str(), 
                      ios::in );

      if( !ifmitprec )
      {
        xtem.flog1 << endl << "Cannot open ";
        xtem.flog1 << xtem.telmnt[0].mitclm.iprecfname;
        xtem.flog1 << " for PREC data input" << endl << endl;
        
        exit( -1 );
      }

      ifmitpet.open( xtem.telmnt[0].mitclm.ipetfname.c_str(), 
                     ios::in );

      if( !ifmitpet )
      {
        xtem.flog1 << endl << "Cannot open ";
        xtem.flog1 << xtem.telmnt[0].mitclm.ipetfname;
        xtem.flog1 << " for PET data input" << endl << endl;
        
        exit( -1 );
      }

      ifmiteet.open( xtem.telmnt[0].mitclm.ieetfname.c_str(), 
                     ios::in );

      if( !ifmiteet )
      {
        xtem.flog1 << endl << "Cannot open ";
        xtem.flog1 << xtem.telmnt[0].mitclm.ieetfname;
        xtem.flog1 << " for EET data input" << endl << endl;
        
        exit( -1 );
      }

      ifmitsh2o1.open( xtem.telmnt[0].mitclm.ish2o1fname.c_str(), 
                       ios::in );

      if( !ifmitsh2o1 )
      {
        xtem.flog1 << endl << "Cannot open ";
        xtem.flog1 << xtem.telmnt[0].mitclm.ish2o1fname;
        xtem.flog1 << " for SOILH2O data input of the top 1 meter";
        xtem.flog1 << endl << endl;
       
        exit( -1 );
      }

      ifmitsh2o2.open( xtem.telmnt[0].mitclm.ish2o2fname.c_str(), 
                       ios::in );

      if( !ifmitsh2o2 )
      {
        xtem.flog1 << endl << "Cannot open ";
        xtem.flog1 << xtem.telmnt[0].mitclm.ish2o2fname;
        xtem.flog1 << " for SOILH2O data input of the top 2 meters";
        xtem.flog1 << endl << endl;
        
        exit( -1 );
      }

      ifmitspack.open( xtem.telmnt[0].mitclm.ispackfname.c_str(), 
                       ios::in );

      if( !ifmitspack )
      {
        xtem.flog1 << endl << "Cannot open ";
        xtem.flog1 << xtem.telmnt[0].mitclm.ispackfname;
        xtem.flog1 << " for SNOWPACK data input";
        xtem.flog1 << endl << endl;
        
        exit( -1 );
      }

      ifmitsrfrun.open( xtem.telmnt[0].mitclm.isrfrunfname.c_str(), 
                        ios::in );

      if( !ifmitsrfrun )
      {
        xtem.flog1 << endl << "Cannot open ";
        xtem.flog1 << xtem.telmnt[0].mitclm.isrfrunfname;
        xtem.flog1 << " for SURFRUN data input";
        xtem.flog1 << endl << endl;
        
        exit( -1 );
      }

      ifmitdrain.open( xtem.telmnt[0].mitclm.idrainfname.c_str(), 
                       ios::in );

      if( !ifmitdrain )
      {
        xtem.flog1 << endl << "Cannot open ";
        xtem.flog1 << xtem.telmnt[0].mitclm.idrainfname;
        xtem.flog1 << " for DRAINAGE data input";
        xtem.flog1 << endl << endl;
        
        exit( -1 );
      }
    }
    
    // Open MIT ozone file

    if( 1 == xtem.telmnt[0].mitclm.to3flag )
    {
      if( xtem.istateflag < 2 && pdyr < (xtem.totsptime+1) )
      {
        co2tstyr = 0;
      }
      else { co2tstyr = tstyr; }

      tempfname.str( "" );
      tempfname << xtem.telmnt[0].mitclm.io3fname
                << (xtem.telmnt[0].mitclm.startyr+co2tstyr)
                << xtem.telmnt[0].mitclm.io3end;
      
      o3name = tempfname.str();

      ifmito3.open( o3name.c_str(), ios::in );
      
      if( !ifmito3 )
      {
        xtem.flog1 << endl << "Cannot open " << o3name;
        xtem.flog1 << " for O3 data input" << endl << endl;
        
        exit( -1 );
      }
    }
    else
    {
      ifmito3.open( xtem.telmnt[0].mitclm.io3fname.c_str(), 
                    ios::in );

      if( !ifmito3 )
      {
        xtem.flog1 << endl << "Cannot open "; 
        xtem.flog1<< xtem.telmnt[0].mitclm.io3fname;
        xtem.flog1 << " for O3 data input" << endl << endl;
        
        exit( -1 );
      }
    }

    // Open MIT atmospheric CO2 file

    if( xtem.telmnt[0].mitclm.tco2flag > 0 )
    {
      if( xtem.istateflag < 2 && pdyr < (xtem.totsptime+1) )
      {
        co2tstyr = 0;
      }
      else { co2tstyr = tstyr; }
    
      tempfname.str( "" );
      tempfname << xtem.telmnt[0].mitclm.ico2fname
                << (xtem.telmnt[0].mitclm.startyr+co2tstyr)
                << xtem.telmnt[0].mitclm.ico2end;
      
      co2name = tempfname.str();

      ifmitco2.open( co2name.c_str(), ios::in );
      
      if( !ifmitco2 )
      {
        xtem.flog1 << endl << "Cannot open " << co2name;
        xtem.flog1 << " for CO2 data input" << endl << endl;
        
        exit( -1 );
      }
    }
    else
    {
      ifmitco2.open( xtem.telmnt[0].mitclm.ico2fname.c_str(), 
                     ios::in );

      if( !ifmitco2 )
      {
        xtem.flog1 << endl << "Cannot open ";
        xtem.flog1 << xtem.telmnt[0].mitclm.ico2fname;
        xtem.flog1 << " for CO2 data input" << endl << endl;
        
        exit( -1 );
      }
    }
  }
  
 
  xtem.updateMITCLMregion( pdyr,
                           pdm, 
                           ifmitnirr, 
                           ifmittair, 
                           ifmitdayTsoilL1,
                           ifmitdayTsoilL2,
                           ifmitdayTsoilL3,
                           ifmitdayTsoilL4,
                           ifmitdayTsoilL5,
                           ifmitdayTsoilL6,
                           ifmitdayTsoilL7,
                           ifmitdayTsoilL8,
                           ifmitdayTsoilL9,
                           ifmitdayTair,
                           ifmitrainDur,
                           ifmitrainInt,
                           ifmitdaySML1,
                           ifmitdaySML2,
                           ifmitdaySML3,
                           ifmitdaySML4,
                           ifmitdaySML5,
                           ifmitdaySML6,
                           ifmitdaySML7,
                           ifmitdaySML8,
                           ifmitdaySML9,
                           ifmithrSML1,
                           ifmithrSML2,
                           ifmithrSML3,
                           ifmithrSML4,
                           ifmithrSML5,
                           ifmithrSML6,
                           ifmitprec,
                           ifmitpet,
                           ifmiteet,
                           ifmitsh2o1,
                           ifmitsh2o2,
                           ifmitspack,
                           ifmitsrfrun,
                           ifmitdrain,
                           ifmitco2,
                           ifmito3 );

 
  if( (CYCLE-1) == pdm )
  {  
    ifmitnirr.close();
 
    ifmittair.close();

    if( 1 == xtem.telmnt[0].mitclm.ch4flag 
        || 1 == xtem.telmnt[0].mitclm.n2oflag )
    {   
      ifmitdayTsoilL1.close(); 
      ifmitdayTsoilL2.close(); 
      ifmitdayTsoilL3.close(); 
      ifmitdayTsoilL4.close(); 
      ifmitdayTsoilL5.close(); 
      ifmitdayTsoilL6.close(); 
      ifmitdayTsoilL7.close(); 
      ifmitdayTsoilL8.close(); 
      ifmitdayTsoilL9.close(); 

      ifmitdaySML1.close();  
      ifmitdaySML2.close();  
      ifmitdaySML3.close();  
      ifmitdaySML4.close();  
      ifmitdaySML5.close();  
      ifmitdaySML6.close();  
      ifmitdaySML7.close();  
      ifmitdaySML8.close();  
      ifmitdaySML9.close();  
    }

    if( 1 == xtem.telmnt[0].mitclm.n2oflag )
    { 
      ifmitdayTair.close();  
      ifmitrainDur.close(); 
      ifmitrainInt.close(); 

      ifmithrSML1.close();  
      ifmithrSML2.close();  
      ifmithrSML3.close();  
      ifmithrSML4.close();  
      ifmithrSML5.close();  
      ifmithrSML6.close();  
    }

    ifmitprec.close();
  
    ifmitpet.close();
  
    ifmiteet.close();
  
    ifmitsh2o1.close();
    ifmitsh2o2.close();
  
    ifmitspack.close();
 
    ifmitsrfrun.close();

    ifmitdrain.close();

    ifmitco2.close();
  
    ifmito3.close();
  }
  
};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void assignTEM2MITCLM( const int& pdyr, 
                       const int& pdm, 
                       RegMITTEM& xtem )
{
  int nepyr;
  ostringstream tempfname;
  
  string cflxname;
  string ch4flxname;
  string n2oflxname;

  static ofstream ftcflx;
  static ofstream ftch4flx;
  static ofstream ftn2oflx;

  if( 0 == pdm )
  {
    if( xtem.istateflag < 2 )
    {
      nepyr = pdyr - xtem.totsptime - 1;
    }
    else
    {
      nepyr = pdyr - 1;
    }

    tempfname.str( "" );
    tempfname << xtem.telmnt[0].mitclm.icflxfname
              << (xtem.telmnt[0].mitclm.startyr+nepyr)
              << xtem.telmnt[0].mitclm.icflxend;
  
    cflxname = tempfname.str();

    ftcflx.open( cflxname.c_str(), ios::app );


    if( 1 == xtem.telmnt[0].tem.ch4flag )
    {
      tempfname.str( "" );
      tempfname << xtem.telmnt[0].mitclm.ich4flxfname
                << (xtem.telmnt[0].mitclm.startyr+nepyr)
                << xtem.telmnt[0].mitclm.ich4flxend;
    
      ch4flxname = tempfname.str();

      ftch4flx.open( ch4flxname.c_str(), ios::app );
    }


    if( 1 == xtem.telmnt[0].tem.n2oflag )
    {
      tempfname.str( "" );
      tempfname << xtem.telmnt[0].mitclm.in2oflxfname
                << (xtem.telmnt[0].mitclm.startyr+nepyr)
                << xtem.telmnt[0].mitclm.in2oflxend;
    
      n2oflxname = tempfname.str();

      ftn2oflx.open( n2oflxname.c_str(), ios::app );
    }
  }


  xtem.telmnt[0].mitclm.tcflx2D.writeclm( ftcflx, 
                                          "CFLUX", 
                                          pdm,
                                          xtem.telmnt[0].mitclm.tcflx2D );


  if( 1 == xtem.telmnt[0].tem.ch4flag )
  {
    xtem.telmnt[0].mitclm.tch4flx2D.writeclm( ftch4flx, 
                                              "CH4FLUX", 
                                              pdm,
                                              xtem.telmnt[0].mitclm.tch4flx2D );
  }


  if( 1 == xtem.telmnt[0].tem.n2oflag )
  {
    xtem.telmnt[0].mitclm.tn2oflx2D.writeclm( ftn2oflx, 
                                              "N2OFLUX", 
                                              pdm,
                                              xtem.telmnt[0].mitclm.tn2oflx2D );
  }

  
  if( (CYCLE-1) == pdm )
  {
    ftcflx.close();

    if( 1 == xtem.telmnt[0].tem.ch4flag )
    {
      ftch4flx.close();
    }

    if( 1 == xtem.telmnt[0].tem.n2oflag )
    {
      ftn2oflx.close();
    }
  }
  
};

/* ***************************************************************
*************************************************************** */


/* ***************************************************************
*************************************************************** */

void initLCLUC( RegMITTEM& xtem )
{
  int dyr;
  int ichrt;
  
  // Hand 0ff startyr from clm
  
  xtem.telmnt[0].lcluc.startyr = xtem.telmnt[0].mitclm.startyr;


  // Get community type parameterizations
  
  xtem.telmnt[0].lcluc.getvtype( xtem.flog1 );
  

  // Set flags land cover /land-use change data

  xtem.telmnt[0].lcluc.setLCLUCFlags( xtem.flog1, xtem.equil );


  // Set filename for spatially explicit data for the maximum
  //   number of cohorts of land cover cohorts in a grid cell

  xtem.telmnt[0].lcluc.initPotvegMaxCohorts( xtem.flog1 );


  // Set filename for spatially explicit data for the land area
  //   in a grid cell
 
  cout << "Please enter the file name containing the land area ";
  cout << "of grid cells in the region:";
  cout << endl;
  cout << "        (e.g., LANDAREA) " << endl;

  cin >> xtem.ilandareafname;
  
  xtem.flog1 << "Please enter the file name containing the land area ";
  xtem.flog1 << "of grid cells in the region:";
  xtem.flog1 << endl;
  xtem.flog1 << "        (e.g., LANDAREA) " << endl  ;

  xtem.flog1 << xtem.ilandareafname << endl << endl;

  if( 0 == xtem.equil )
  {
    xtem.telmnt[0].lcluc.initCohorts( xtem.flog1 );
  }
  else
  {
    xtem.telmnt[0].lcluc.initPotvegCohorts( xtem.flog1 );
  }

  dyr = 0;
  updateLCLUC( dyr, xtem ); 
     
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void initMITCLM( RegMITTEM& xtem )
{
  double co2level;
  int dm;
  int dyr;
  double initco2;
   
  // Set filenames for spatially explicit cloudiness or 
  //   solar radiation data
  
  xtem.telmnt[0].mitclm.setMITSolarRadFlags( xtem.flog1, xtem.equil );
  xtem.telmnt[0].mitclm.initMITSolarRad( xtem.flog1 ); 

  
 // Set filenames for spatially explicit air temperature data

  xtem.telmnt[0].mitclm.setTairFlags( xtem.flog1, xtem.equil );
  xtem.telmnt[0].mitclm.initTair( xtem.flog1 );


 // Set flag for constant or transient hydrology

  xtem.telmnt[0].mitclm.setPrecFlags( xtem.flog1, xtem.equil );
  xtem.telmnt[0].mitclm.initPrec( xtem.flog1 ); 
 
 // Set filenames for spatially explicit potential 
 //   evapotranspiration data
  
  xtem.telmnt[0].mitclm.initPET( xtem.flog1 );
  
  // Set filenames for spatially explicit estimated actual 
  //   evapotranspiration data
  
  xtem.telmnt[0].mitclm.initEET( xtem.flog1 );

  // Set filenames for spatially explicit soil moisture data
  
  xtem.telmnt[0].mitclm.initSoilH2O( xtem.flog1 );

  // Set filenames for spatially explicit snowpack data

  xtem.telmnt[0].mitclm.initSNOWPACK( xtem.flog1 );

  // Set filenames for spatially explicit surface runoff data

  xtem.telmnt[0].mitclm.initSURFRUN( xtem.flog1 );

  // Set filenames for spatially explicit drainage data

  xtem.telmnt[0].mitclm.initDRAINAGE( xtem.flog1 );

  // Set filenames for spatially explicit ozone data

  xtem.telmnt[0].mitclm.setO3Flags( xtem.flog1, xtem.equil );
  xtem.telmnt[0].mitclm.initO3( xtem.flog1 );

 	
  // Determine initial atmospheric CO2 concentration for initial 
  //  equilibrium portion of the simulation
  
  cout << endl << endl;
  cout << "Enter the initial concentration of carbon";
  cout << " dioxide in ppmv: ";

  cin >> initco2;

  xtem.telmnt[0].mitclm.setINITCO2( initco2 );

  xtem.flog1 << endl << endl;
  xtem.flog1 << "Enter the initial concentration of carbon";
  xtem.flog1 << " dioxide in ppmv: ";
  xtem.flog1 << xtem.telmnt[0].mitclm.getINITCO2() << endl << endl;

  cout << "Enter the final equilibrium concentration of";
  cout << " carbon dioxide in ppmv: ";

  cin >> co2level;

  xtem.telmnt[0].mitclm.setCO2LEVEL( co2level );

  xtem.flog1 << "Enter the final equilibrium concentration of";
  xtem.flog1 << " carbon dioxide in ppmv: ";
  xtem.flog1 << xtem.telmnt[0].mitclm.getCO2LEVEL() << endl << endl;
  
  
  // Get annual atmospheric CO2 concentration data for transient 
  //   portion of the simulation
  
  xtem.telmnt[0].mitclm.setCO2Flags( xtem.flog1, xtem.equil );
  xtem.telmnt[0].mitclm.initCO2( xtem.flog1 );
  
  if( 0 == xtem.telmnt[0].mitclm.tco2flag )
  {
    cout << "Please enter the file name containing ";
    cout << "the annual atmospheric CO2 data" << endl;
    cout << "               (e.g., CO2) " << endl;

    cin >> xtem.telmnt[0].mitclm.ico2fname;

    xtem.flog1 << "Please enter the file name containing ";
    xtem.flog1 << "the annual transient atmospheric CO2 ";
    xtem.flog1 << "data: " << endl;
    xtem.flog1 << "               (e.g., CO2) " << endl;
    xtem.flog1 << xtem.telmnt[0].mitclm.ico2fname << endl << endl;
  }


  // Set file names for CFLUX output from TEM

  xtem.telmnt[0].mitclm.initMITCflux( xtem.flog1 );

  // Set file names for CH4FLUX output from TEM

  xtem.telmnt[0].mitclm.setMITCH4Flags( xtem.flog1 );

  xtem.telmnt[0].tem.mdmflag = 0;
  
  if( 1 == xtem.telmnt[0].mitclm.ch4flag )
  {
    cout << "Do you want to run the Methane Dynamics Model?:";
    cout << endl;
    cout << "Enter 0 for No (run NEM-CH4 instead):" << endl;
    cout << "Enter 1 for Yes: ";
    
    cin >> xtem.telmnt[0].tem.mdmflag;

    xtem.flog1 << "Do you want to run the Methane Dynamics Model?:";
    xtem.flog1 << endl;
    xtem.flog1 << "Enter 0 for No (run NEM-CH4 instead):" << endl;
    xtem.flog1 << "Enter 1 for Yes: ";    
    xtem.flog1 << xtem.telmnt[0].tem.mdmflag << endl << endl;
  }

  xtem.telmnt[0].mitclm.initMITCH4flux( xtem.flog1 );

  
  // Set file names for N2OFLUX output from TEM

  xtem.telmnt[0].mitclm.setMITN2OFlags( xtem.flog1 );
  xtem.telmnt[0].mitclm.initMITN2Oflux( xtem.flog1 );

  if ( 1 == xtem.telmnt[0].mitclm.tcldsflag
       && 1 == xtem.telmnt[0].mitclm.ttairflag 
       && 1 == xtem.telmnt[0].mitclm.tprecflag )
  {
    dyr = 0;
    dm = 0;

    assignMITCLM2TEM( dyr, dm, xtem );
  }
  else
  {
    dyr = 0;
    dm = 0;

    assignMITCLM2TEM( dyr, dm, xtem );


//    cout << "Source code has not yet been evaluated";
//    cout << " for this set of conditions" << endl;

//    xtem.flog1 << "Source code has not yet been evaluated";
//    xtem.flog1 << " for this set of conditions" << endl;

//    exit( 0 );
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void initTEM( RegMITTEM& xtem )
{
  int dm;
  int dyr; 
  
  // Initialize flags in TEM to specify operating conditions 
  //  for the simulation

  xtem.starttem();

  // Initialize 1nitClimate, initAET and initSoilH2O arrays with
  //   12 months of data from the IGSM and run TEM to equilibrium 
  //   conditions using the baseline climate if starting from 
  //   calibration data (i.e. 0 == xtem.istateflag ); or initialize 
  //   TEM with December climate if starting with a temstate file

  dyr = 0;

  if( 0 == xtem.istateflag ) 
  { 

    dm = 0;
    xtem.initializeMITTEMregion( dm );

    for( dm = 1; dm < CYCLE; ++dm )
    {
      assignMITCLM2TEM( dyr, dm, xtem );
      xtem.initializeMITTEMregion( dm );
    }     
  }
  else
  {
    // Read in 12 months of data from input data sets

    for( dm = 1; dm < CYCLE; ++dm )
    {
      assignMITCLM2TEM( dyr, dm, xtem );
    }


    // Use December climate for initialization

    dm = CYCLE-1;
    xtem.initializeMITTEMregion( dm );
  }
        	
};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void updateLCLUC( int& pdyr, RegMITTEM& xtem )
{  
  FILE* iflulc;

  FILE* iflandarea;

  FILE* ifnumchrts;

  string lulcname;

  string mxcohrtname;
  
  ostringstream tempfname;

  int tstyr;

  if( 0 == pdyr )
  {
    // Open maximum cohort file 

    ifnumchrts = fopen( xtem.telmnt[0].lcluc.ipotmxcohrtfname.c_str(), "r" );

    if( !ifnumchrts )
    {
      xtem.flog1 << endl << "Cannot open " << xtem.telmnt[0].lcluc.ipotmxcohrtfname;
      xtem.flog1 << " for MXCOHRTS data input";
      xtem.flog1 << endl << endl;

      exit( -1 );
    }

    // Open land area file

    iflandarea = fopen( xtem.ilandareafname.c_str(), "r" );

    if( !iflandarea )
    {
      xtem.flog1 << endl << "Cannot open " << xtem.ilandareafname;
      xtem.flog1 << " for LANDAREA data input";
      xtem.flog1 << endl << endl;

      exit( -1 );
    }
  }

 
  if( 0 == xtem.telmnt[0].lcluc.tlulcflag )
//      || (0 == xtem.istateflag && 0 == pdyr) )
  {
    // Open land use/land cover cohort file

    iflulc = fopen( xtem.telmnt[0].lcluc.ipotlulcfname.c_str(), "r" );

    if( !iflulc )
    {
      xtem.flog1 << endl << "Cannot open " << xtem.telmnt[0].lcluc.ipotlulcfname;
      xtem.flog1 << " for POTVEG LULCCHRT data input";
      xtem.flog1 << endl << endl;

      exit( -1 );
    }  	
  }
  else	  
  {   
    if( 2 == xtem.istateflag && 0 == pdyr ) { tstyr = -1; }
    else if( xtem.istateflag < 2 && pdyr < (xtem.totsptime+1) )
    {
      tstyr = 0;
    }
    else
    {
      if( xtem.istateflag < 2 )
      {
        tstyr = pdyr - xtem.totsptime - 1;
      }
      else { tstyr = pdyr - 1; }
    }


    // Open the same LULC data file (i.e., the first year of data) 
    //   for spin-up  
      
    tempfname.str( "" );
    
    tempfname << xtem.telmnt[0].lcluc.ilulcfname
              << (xtem.telmnt[0].mitclm.startyr+tstyr) 
              << xtem.telmnt[0].lcluc.ilulcend;
      
    lulcname = tempfname.str();

    iflulc = fopen( lulcname.c_str(), "r" );

    if( !iflulc )
    {
      xtem.flog1 << endl << "Cannot open " << lulcname;
      xtem.flog1 << " for data input" << endl << endl;

      exit( -1 );
    }
  }
  
  xtem.updateLCLUCregion( pdyr, ifnumchrts, iflandarea, iflulc );

//  fclose( ifnumchrts );
   
//  fclose( iflandarea );

//  fclose( iflulc );

};
