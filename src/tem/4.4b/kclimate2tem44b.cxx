/* This file contains four subroutines which are called by the driver:
   tem_init, tem, tem2climate and climate2tem
*/

#include<cstdio>
 
  using std::fclose;
  using std::FILE;
  using std::fopen;
  using std::fscanf;
  using std::printf;
  using std::sprintf;
 
#include<iostream>

  using std::cerr;
  using std::cin;
  using std::cout;
  using std::endl;
  using std::ios;
 
#include<fstream>

  using std::ifstream;
  using std::ofstream;

#include<iomanip>

  using std::setiosflags;
  using std::setprecision;
  using std::setw;


#include<cstdlib>

  using std::exit;

#include<cmath>


#include<vector>

  using std::vector;

#include<cctype>        // Use with ANSII C++

  using std::toupper;

#include<cstring>       // Use with ANSII C++

#include<string>

  using std::string;

#include<sstream>

  using std::ostringstream;

#include<ctime>

  using std::ctime;
  using std::time_t;

#include"nemconsts44a.hpp"

#include "regmittem44b.h"

extern "C" {

//  #include "climate4tem437.hpp"
//  extern "C" struct
  //{
  //double co2[MAXNGRD];     // Monthly atmospheric CO2
  //double o3[MAXNGRD][MAX3HRS]; // Monthly atmospheric O3
  //double temp[MAXNGRD];    // Monthly air temperature
  //double daytemp[MAXNGRD][MXMTHDAYS]; //Daily air temperature
  //double swrs[MAXNGRD];    // Monthly short-wave radiation at surface
  //double pre[MAXNGRD];     // Monthly precipitation
  //double strmdur[MAXNGRD][MXMTHDAYS];  // Daily duration of storms
  //double qstrm[MAXNGRD][MXMTHDAYS];    // Daily intensity of storms
  //double aet[MAXNGRD][MXNCHRTS];     // Monthly actual evapotranspiration
  //double sh2o1m[MAXNGRD][MXNCHRTS];  // Monthly soil moisture in top 1 m
  //double sh2o2m[MAXNGRD][MXNCHRTS];  // Monthly soil moisture in top 2 m?
  //double daytsoil[MXNLAYERS][MAXNGRD][MXNCHRTS][MXMTHDAYS];  // Daily soil temperature
  //double daysh2o[MXNLAYERS][MAXNGRD][MXNCHRTS][MXMTHDAYS];   // Daily soil moisture
  //double hrsh2o[MXNLAYERS][MAXNGRD][MXNCHRTS][MXMTHDAYS][MXDAYHRS]; // Hourly soil moisture
  //int storms[MAXNGRD][MXMTHDAYS];  // Daily number of storms in a month
  //int strmdry[MAXNGRD][MXMTHDAYS]; // Number of dry days per month
  //} climate4tem_;

 extern "C" struct { 
    double nep4chem[MAXNGRD];
    double ch44chem[MAXNGRD];
    double n2o4chem[MAXNGRD];
  } upt4chem_;

  void updatelcluc_( const int& pdyr );

RegMITTEM xtem;

/* *************************************************************
************************************************************* */



/* *************************************************************
************************************************************* */

  extern void tem2climate_( const int& outyr, const int& dm )
  {
    // This subroutine "averages" data from TEM grid to climate 
    //   model grid

    int dmlat;
  
    for( dmlat=0; dmlat <  MAXNGRD; ++dmlat ) 
    { 
      upt4chem_.nep4chem[dmlat] = 0.0;
      upt4chem_.ch44chem[dmlat] = 0.0;
      upt4chem_.n2o4chem[dmlat] = 0.0;
    }

    for( dmlat=0; dmlat <  MAXNGRD; ++dmlat )
    {
      upt4chem_.nep4chem[dmlat] = upt4chem_.nep4chem[dmlat]
                                  + xtem.telmnt[0].mitclm.tcflx2D.latband[dmlat];
      upt4chem_.ch44chem[dmlat] = upt4chem_.ch44chem[dmlat]
                                  + xtem.telmnt[0].mitclm.tch4flx2D.latband[dmlat];
      upt4chem_.n2o4chem[dmlat] = upt4chem_.n2o4chem[dmlat]
                                  + xtem.telmnt[0].mitclm.tn2oflx2D.latband[dmlat];
    }

   //cout << "NEP4CHEM from TEM=" << upt4chem_.nep4chem[23]<< endl ;

   //cout << "NEP4CHEM from TEM=" <<  xtem.telmnt[0].mitclm.tcflx2D.latband[23] << endl ;

  }

/* *************************************************************
************************************************************* */

  extern void tem_( int& year, int& mn )
  {
    // Extrapolates TEM across the globe for a single year
  
    cout << " Call TEM  " << year << " " <<  mn << endl;
    
    xtem.updateMITTEMregion( year,mn );
  }


/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

  extern void temclm_init_(void)
  {
    double co2level;
    double initco2;
    
    // Initialize TEM flags, parameters, etc.
    
    cout << "From  temclm_init before xtem.initrun " << endl;
    
//  xtem.flog1.open( "mittem4.log",ios::out );

    xtem.setRunMode();

    xtem.setInTEMState();

    xtem.setRunTime();

    xtem.setGridParameters();

    xtem.telmnt[0].mitclm.setMITSolarRadFlags( xtem.flog1, 
                                               xtem.equil );

    xtem.telmnt[0].mitclm.setTairFlags( xtem.flog1, xtem.equil );

    xtem.telmnt[0].mitclm.setPrecFlags( xtem.flog1, xtem.equil );

    xtem.telmnt[0].mitclm.setO3Flags( xtem.flog1, xtem.equil );

    // Determine initial atmospheric CO2 concentration for 
    //   initial equilibrium portion of the simulation

    cout << endl << endl;
    cout << "Enter the initial concentration of carbon ";
    cout << "dioxide in ppmv: ";
  
    cin >> initco2;
    cout << initco2 << endl;
  
    xtem.telmnt[0].mitclm.setINITCO2( initco2 );

    xtem.flog1 << endl << endl;
    xtem.flog1 << "Enter the initial concentration of carbon ";
    xtem.flog1 << "dioxide in ppmv: ";
    xtem.flog1 << xtem.telmnt[0].mitclm.getINITCO2();
    xtem.flog1 << endl << endl;

    cout << "Enter the final equilibrium concentration of ";
    cout << "carbon dioxide in ppmv: " << endl;
  
    cin >> co2level;
    cout << co2level << endl;
  
    xtem.telmnt[0].mitclm.setCO2LEVEL( co2level );
    
    xtem.flog1 << "Enter the final equilibrium concentration ";
    xtem.flog1 << "of carbon dioxide in ppmv: ";
    xtem.flog1 << xtem.telmnt[0].mitclm.getCO2LEVEL();
    xtem.flog1 << endl << endl;

//  cout << "From  temclm_init CO2 concentration " << xtem.telmnt[0].mitclm.co2level << endl;

    xtem.telmnt[0].mitclm.setCO2Flags( xtem.flog1, xtem.equil );
    
    xtem.telmnt[0].mitclm.setMITCH4Flags( xtem.flog1 );
  
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

    xtem.telmnt[0].mitclm.setMITN2OFlags( xtem.flog1 );

  }


/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */
  
  extern void lcluc_init_( void )
  {
    int dyr;

    // Hand off startyr from mitclm

    xtem.telmnt[0].lcluc.startyr = xtem.telmnt[0].mitclm.startyr;

    // Get community type parameterizations

    xtem.telmnt[0].lcluc.getvtype( xtem.flog1 );

    xtem.telmnt[0].lcluc.setLCLUCFlags( xtem.flog1,
                                        xtem.equil );

    // Set filename for spatially explicit number of cohorts data

    xtem.telmnt[0].lcluc.initMaxCohorts( xtem.flog1 );


    // Set filename for spatially explicit land use/ land cover 
    //   cohort data

    xtem.telmnt[0].lcluc.initCohorts( xtem.flog1 );

    dyr = 0;
    updatelcluc_( dyr );

  }


/* *************************************************************
************************************************************* */

  extern void tem_cleanup_( void )
  {
    // Close files and free up memory used by TEM
    
    //  xtem.cleanup();

  }

/* *************************************************************
************************************************************* */

  void updatelcluc_( const int& pdyr )
  {
    string lulcname;
    string mxcohrtname;
    ostringstream tempfname;

    FILE* ifnumchrts;
    FILE* iflulc;
  
    int tstyr;
  
    if( 1 == xtem.telmnt[0].lcluc.tlulcflag )
    {

      if( 0 == pdyr ) { tstyr = 0; }
      else if( xtem.istateflag < 2 && pdyr < (xtem.totsptime+1) )
      {
        tstyr = 1;
      }
      else
      {
        if( xtem.istateflag < 2 )
        {
          tstyr = pdyr - xtem.totsptime;
        }
        else { tstyr = pdyr; }
      }

      // Read in land cover/land use change data sets for the 
      //   transient spin-up portion of TEM initialization


      // Open maximum cohort file

      tempfname.str( "" ); 
      tempfname << xtem.telmnt[0].lcluc.imxcohrtfname
                << (xtem.telmnt[0].mitclm.startyr+tstyr) 
                << xtem.telmnt[0].lcluc.imxcohrtend;

      mxcohrtname = tempfname.str();

      ifnumchrts = fopen( mxcohrtname.c_str(), "r" );

      if( !ifnumchrts )
      {
        xtem.flog1 << endl << "Cannot open " << mxcohrtname;
        xtem.flog1 << " for MXCOHRTS data input" << endl << endl;

        exit( -1 );
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
        xtem.flog1 << "\nCannot open " << lulcname;
        xtem.flog1 << " for data input" << endl << endl;
      
        exit( -1 );
      }
    }
    else
    {
      // Open maximum cohort file

      tempfname.str( "" );
      tempfname << xtem.telmnt[0].lcluc.imxcohrtfname
                << xtem.telmnt[0].lcluc.imxcohrtend;

      mxcohrtname = tempfname.str();

      ifnumchrts = fopen( mxcohrtname.c_str(), "r" );

      if( !ifnumchrts )
      {
        xtem.flog1 << endl << "Cannot open " << mxcohrtname;
        xtem.flog1 << " for MXCOHRTS data input" << endl << endl;

        exit( -1 );
      }

      // Open land use/land cover cohort file

      tempfname.str( "" );
      tempfname << xtem.telmnt[0].lcluc.ilulcfname
                << xtem.telmnt[0].lcluc.ilulcend;

      lulcname = tempfname.str();

      iflulc = fopen( lulcname.c_str(), "r" );

      if( !iflulc )
      {
        xtem.flog1 << endl << "Cannot open " << lulcname;
        xtem.flog1 << " for LULCCHRT data input" << endl << endl;

        exit( -1 );
      }
    }

    xtem.updateLCLUCregion( pdyr, ifnumchrts, iflulc );
    
    fclose( ifnumchrts );
    fclose( iflulc );

  }

/* *************************************************************
************************************************************* */

  extern void tem_init_()
  {
    int dm;
  
    xtem.starttem();
  
    if( 0 == xtem.istateflag ) { dm = 0; }
    else { dm = CYCLE-1; }
  
    xtem.initializeMITTEMregion( dm );
  
    cout << "From  tem_init " << endl;

  }

}     // For extern "C"
