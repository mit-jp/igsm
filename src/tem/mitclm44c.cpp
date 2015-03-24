/* *************************************************************
****************************************************************
MITCLM44A.CPP - physical characteristics of the atmosphere
                  as described by the MIT 2D L-O GCM

20060128 - DWK created by modifying mitclm436a.cpp
20060128 - DWK added include mitclm437a.h and standard includes
20060201 - DWK renamed mkMITd40() to be mkMITaot40()
20060320 - DWK added initDRAINAGE(), initSNOWPACK() and
           initSURFRUN()
20060420 - DWK added compiler directive STANDALONE_TEM
20080130 - DWK changed include from mitclm437.h to mitclm44a.h
20080130 - DWK changed MITclm43:: to MITclm44::
20080130 - DWK changed TEMclm43() to TEMclm44()
20110707 - DWK changed include mitclm44a.h to mitclm44c.h
          
****************************************************************
************************************************************* */


// NOTE: If running TEM "offline", make sure the compiler 
//         directive STANDALONE_TEM is DEFINED below.  If  
//         coupling TEM to the IGSM, make sure STANDALONE_TEM is 
//         NOT DEFINED (i.e. comment out next line)

#include "preproc.h"

#include<cstdio>

  using std::fscanf;
  using std::FILE;

#include<iostream>

  using std::cout;
  using std::cin;
  using std::ios;
  using std::endl;

#include<fstream>

  using std::ifstream;
  using std::ofstream;

#ifdef BORLAND_CPP
  #include<stdlib>
#else 
  #include<cstdlib>
#endif

  using std::exit;

#include<iomanip>

  using std::setprecision;
  using std::setw;

#include<vector>

  using std::vector;

#include<string>

  using std::string;


#include "mitclm44c.h"

/* *************************************************************
************************************************************* */

MITclm44::MITclm44() : TEMclm44()
{
  int dday;
  int dm;
  int dmlat;
//  int dmlon;

  // Assign the number of 0.5 degree latitudinal bands for 
  //   each MIT 2-D L-O GCM latitudinal band

  latnum[0] = 4;   // 90.0 degrees S to 88.0 degrees S
  
  // 88.0 degrees S to 88.0 degrees N by 4 degree 
  //   latitudinal band increments
  for( dmlat = 1; dmlat < (MXMITNLAT-1); ++dmlat )
  {
    latnum[dmlat] = 8;   
  }

  latnum[MXMITNLAT-1] = 4; // 88.0 degrees N to 90.0 degrees N
  
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



  for( dm = 0; dm < CYCLE; ++dm ) 
  {  
    resetCFLUX( dm );  
    resetCH4FLUX( dm );  
    resetN2OFLUX( dm );  
  }

  for( dday = 0; dday < MAXMDAYS; ++dday )
  {
    dayTair[dday] = MISSING;
    rainDuration[dday] = MISSING;
    rainIntensity[dday] = MISSING;
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITclm44::aggregFlux2D( const int& year,
                             const float& lat, 
                             const double& subarea, 
                             const double& influx,
                             MITdata44& outflux )
{

  int dmlat = 0;

  float pos = -90.0;

  double tstflx;

  if( influx <= MISSING ) { tstflx = ZERO; }
  else { tstflx = influx; }

  while( lat > pos ) 
  {
    pos += latnum[dmlat] * 0.5;

    ++dmlat;
  }

  if( lat < pos ) { --dmlat; }

  outflux.year = year;

  outflux.latband[dmlat] += (tstflx * subarea * 0.000001);

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

#ifdef STANDALONE_TEM
void MITclm44::initDaySoilMoist( ofstream& rflog1 )
{

  if( 1 == ch4flag 
      || 1 == n2oflag )
  {
    if( 1 == tprecflag )
    {
      cout << "Please enter the first part of the file name ";
      cout << "containing the daily soil moisture data";
      cout << "for soil layer 1";
      cout << endl;
      cout << "               (e.g., DAYSM) " << endl;
    
      cin >> idaySMfnameL1;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> idaySMendL1;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the daily soil moisture data";
      rflog1 << "for soil layer 1";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYSM) " << endl;
      rflog1 << idaySMfnameL1 << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << idaySMendL1 << endl << endl;


      cout << "Please enter the first part of the file name ";
      cout << "containing the daily soil moisture data";
      cout << "for soil layer 2";
      cout << endl;
      cout << "               (e.g., DAYSM) " << endl;
    
      cin >> idaySMfnameL2;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> idaySMendL2;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the daily soil moisture data";
      rflog1 << "for soil layer 2";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYSM) " << endl;
      rflog1 << idaySMfnameL2 << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << idaySMendL2 << endl << endl;


      cout << "Please enter the first part of the file name ";
      cout << "containing the daily soil moisture data";
      cout << "for soil layer 3";
      cout << endl;
      cout << "               (e.g., DAYSM) " << endl;
    
      cin >> idaySMfnameL3;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> idaySMendL3;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the daily soil moisture data";
      rflog1 << "for soil layer 3";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYSM) " << endl;
      rflog1 << idaySMfnameL3 << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << idaySMendL3 << endl << endl;


      cout << "Please enter the first part of the file name ";
      cout << "containing the daily soil moisture data";
      cout << "for soil layer 4";
      cout << endl;
      cout << "               (e.g., DAYSM) " << endl;
    
      cin >> idaySMfnameL4;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> idaySMendL4;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the daily soil moisture data";
      rflog1 << "for soil layer 4";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYSM) " << endl;
      rflog1 << idaySMfnameL4 << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << idaySMendL4 << endl << endl;


      cout << "Please enter the first part of the file name ";
      cout << "containing the daily soil moisture data";
      cout << "for soil layer 5";
      cout << endl;
      cout << "               (e.g., DAYSM) " << endl;
    
      cin >> idaySMfnameL5;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> idaySMendL5;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the daily soil moisture data";
      rflog1 << "for soil layer 5";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYSM) " << endl;
      rflog1 << idaySMfnameL5 << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << idaySMendL5 << endl << endl;


      cout << "Please enter the first part of the file name ";
      cout << "containing the daily soil moisture data";
      cout << "for soil layer 6";
      cout << endl;
      cout << "               (e.g., DAYSM) " << endl;
    
      cin >> idaySMfnameL6;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> idaySMendL6;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the daily soil moisture data";
      rflog1 << "for soil layer 6";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYSM) " << endl;
      rflog1 << idaySMfnameL6 << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << idaySMendL6 << endl << endl;


      cout << "Please enter the first part of the file name ";
      cout << "containing the daily soil moisture data";
      cout << "for soil layer 7";
      cout << endl;
      cout << "               (e.g., DAYSM) " << endl;
    
      cin >> idaySMfnameL7;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> idaySMendL7;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the daily soil moisture data";
      rflog1 << "for soil layer 7";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYSM) " << endl;
      rflog1 << idaySMfnameL7 << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << idaySMendL7 << endl << endl;


      cout << "Please enter the first part of the file name ";
      cout << "containing the daily soil moisture data";
      cout << "for soil layer 8";
      cout << endl;
      cout << "               (e.g., DAYSM) " << endl;
    
      cin >> idaySMfnameL8;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> idaySMendL8;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the daily soil moisture data";
      rflog1 << "for soil layer 8";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYSM) " << endl;
      rflog1 << idaySMfnameL8 << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << idaySMendL8 << endl << endl;


      cout << "Please enter the first part of the file name ";
      cout << "containing the daily soil moisture data";
      cout << "for soil layer 9";
      cout << endl;
      cout << "               (e.g., DAYSM) " << endl;
    
      cin >> idaySMfnameL9;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> idaySMendL9;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the daily soil moisture data";
      rflog1 << "for soil layer 9";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYSM) " << endl;
      rflog1 << idaySMfnameL9 << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << idaySMendL9 << endl << endl;
    }
    else
    {
      cout << "Please enter the file name ";
      cout << "containing the daily soil moisture data";
      cout << "for soil layer 1";
      cout << endl;
      cout << "               (e.g., DAYSM) " << endl;
    
      cin >> idaySMfnameL1;

      rflog1 << "Please enter the file name ";
      rflog1 << "containing the daily soil moisture data";
      rflog1 << "for soil layer 1";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYSM) " << endl;
      rflog1 << idaySMfnameL1 << endl << endl;


      cout << "Please enter the file name ";
      cout << "containing the daily soil moisture data";
      cout << "for soil layer 2";
      cout << endl;
      cout << "               (e.g., DAYSM) " << endl;
    
      cin >> idaySMfnameL2;

      rflog1 << "Please enter the file name ";
      rflog1 << "containing the daily soil moisture data";
      rflog1 << "for soil layer 2";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYSM) " << endl;
      rflog1 << idaySMfnameL2 << endl << endl;


      cout << "Please enter the file name ";
      cout << "containing the daily soil moisture data";
      cout << "for soil layer 3";
      cout << endl;
      cout << "               (e.g., DAYSM) " << endl;
    
      cin >> idaySMfnameL3;

      rflog1 << "Please enter the file name ";
      rflog1 << "containing the daily soil moisture data";
      rflog1 << "for soil layer 3";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYSM) " << endl;
      rflog1 << idaySMfnameL3 << endl << endl;


      cout << "Please enter the file name ";
      cout << "containing the daily soil moisture data";
      cout << "for soil layer 4";
      cout << endl;
      cout << "               (e.g., DAYSM) " << endl;
    
      cin >> idaySMfnameL4;

      rflog1 << "Please enter the file name ";
      rflog1 << "containing the daily soil moisture data";
      rflog1 << "for soil layer 4";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYSM) " << endl;
      rflog1 << idaySMfnameL4 << endl << endl;


      cout << "Please enter the file name ";
      cout << "containing the daily soil moisture data";
      cout << "for soil layer 5";
      cout << endl;
      cout << "               (e.g., DAYSM) " << endl;
    
      cin >> idaySMfnameL5;

      rflog1 << "Please enter the file name ";
      rflog1 << "containing the daily soil moisture data";
      rflog1 << "for soil layer 5";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYSM) " << endl;
      rflog1 << idaySMfnameL5 << endl << endl;


      cout << "Please enter the file name ";
      cout << "containing the daily soil moisture data";
      cout << "for soil layer 6";
      cout << endl;
      cout << "               (e.g., DAYSM) " << endl;
    
      cin >> idaySMfnameL6;

      rflog1 << "Please enter the file name ";
      rflog1 << "containing the daily soil moisture data";
      rflog1 << "for soil layer 6";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYSM) " << endl;
      rflog1 << idaySMfnameL6 << endl << endl;


      cout << "Please enter the file name ";
      cout << "containing the daily soil moisture data";
      cout << "for soil layer 7";
      cout << endl;
      cout << "               (e.g., DAYSM) " << endl;
    
      cin >> idaySMfnameL7;

      rflog1 << "Please enter the file name ";
      rflog1 << "containing the daily soil moisture data";
      rflog1 << "for soil layer 7";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYSM) " << endl;
      rflog1 << idaySMfnameL7 << endl << endl;


      cout << "Please enter the file name ";
      cout << "containing the daily soil moisture data";
      cout << "for soil layer 8";
      cout << endl;
      cout << "               (e.g., DAYSM) " << endl;
    
      cin >> idaySMfnameL8;

      rflog1 << "Please enter the file name ";
      rflog1 << "containing the daily soil moisture data";
      rflog1 << "for soil layer 8";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYSM) " << endl;
      rflog1 << idaySMfnameL8 << endl << endl;


      cout << "Please enter the file name ";
      cout << "containing the daily soil moisture data";
      cout << "for soil layer 9";
      cout << endl;
      cout << "               (e.g., DAYSM) " << endl;
    
      cin >> idaySMfnameL9;

      rflog1 << "Please enter the file name ";
      rflog1 << "containing the daily soil moisture data";
      rflog1 << "for soil layer 9";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYSM) " << endl;
      rflog1 << idaySMfnameL9 << endl << endl;
    }
  }

};
#endif

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

#ifdef STANDALONE_TEM
void MITclm44::initDayTair( ofstream& rflog1 )
{

  if( 1 == n2oflag )
  {
    if( 1 == ttairflag )
    {
      cout << "Please enter the first part of the file name ";
      cout << "containing the daily air temperature data";
      cout << endl;
      cout << "               (e.g., DAYTAIR) " << endl;
    
      cin >> idayTairfname;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> idayTairend;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the daily air temperature data";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYTAIR) " << endl;
      rflog1 << idayTairfname << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << idayTairend << endl << endl;
    }
    else
    {
      cout << "Please enter the file name ";
      cout << "containing the daily air temperature data";
      cout << endl;
      cout << "               (e.g., DAYTAIR) " << endl;
    
      cin >> idayTairfname;

      rflog1 << "Please enter the file name ";
      rflog1 << "containing the daily air temperature data";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYTAIR) " << endl;
      rflog1 << idayTairfname << endl << endl;
    }
  }


};
#endif

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

#ifdef STANDALONE_TEM
void MITclm44::initDayTsoil( ofstream& rflog1 )
{

  if( 1 == ch4flag 
      || 1 == n2oflag )
  {
    if( 1 == ttairflag )
    {
      cout << "Please enter the first part of the file name ";
      cout << "containing the daily soil temperature data ";
      cout << "for soil layer 1";
      cout << endl;
      cout << "               (e.g., DAYTSOIL) " << endl;
    
      cin >> idayTsoilfnameL1;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> idayTsoilendL1;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the daily soil temperature data";
      rflog1 << "for soil layer 1";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYTSOIL) " << endl;
      rflog1 << idayTsoilfnameL1 << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << idayTsoilendL1 << endl << endl;


      cout << "Please enter the first part of the file name ";
      cout << "containing the daily soil temperature data";
      cout << "for soil layer 2";
      cout << endl;
      cout << "               (e.g., DAYTSOIL) " << endl;
    
      cin >> idayTsoilfnameL2;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> idayTsoilendL2;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the daily soil temperature data";
      rflog1 << "for soil layer 2";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYTSOIL) " << endl;
      rflog1 << idayTsoilfnameL2 << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << idayTsoilendL2 << endl << endl;


      cout << "Please enter the first part of the file name ";
      cout << "containing the daily soil temperature data";
      cout << "for soil layer 3";
      cout << endl;
      cout << "               (e.g., DAYTSOIL) " << endl;
    
      cin >> idayTsoilfnameL3;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> idayTsoilendL3;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the daily soil temperature data";
      rflog1 << "for soil layer 3";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYTSOIL) " << endl;
      rflog1 << idayTsoilfnameL3 << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << idayTsoilendL3 << endl << endl;


      cout << "Please enter the first part of the file name ";
      cout << "containing the daily soil temperature data";
      cout << "for soil layer 4";
      cout << endl;
      cout << "               (e.g., DAYTSOIL) " << endl;
    
      cin >> idayTsoilfnameL4;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> idayTsoilendL4;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the daily soil temperature data";
      rflog1 << "for soil layer 4";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYTSOIL) " << endl;
      rflog1 << idayTsoilfnameL4 << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << idayTsoilendL4 << endl << endl;


      cout << "Please enter the first part of the file name ";
      cout << "containing the daily soil temperature data";
      cout << "for soil layer 5";
      cout << endl;
      cout << "               (e.g., DAYTSOIL) " << endl;
    
      cin >> idayTsoilfnameL5;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> idayTsoilendL5;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the daily soil temperature data";
      rflog1 << "for soil layer 5";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYTSOIL) " << endl;
      rflog1 << idayTsoilfnameL5 << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << idayTsoilendL5 << endl << endl;


      cout << "Please enter the first part of the file name ";
      cout << "containing the daily soil temperature data";
      cout << "for soil layer 6";
      cout << endl;
      cout << "               (e.g., DAYTSOIL) " << endl;
    
      cin >> idayTsoilfnameL6;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> idayTsoilendL6;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the daily soil temperature data";
      rflog1 << "for soil layer 6";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYTSOIL) " << endl;
      rflog1 << idayTsoilfnameL6 << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << idayTsoilendL6 << endl << endl;


      cout << "Please enter the first part of the file name ";
      cout << "containing the daily soil temperature data";
      cout << "for soil layer 7";
      cout << endl;
      cout << "               (e.g., DAYTSOIL) " << endl;
    
      cin >> idayTsoilfnameL7;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> idayTsoilendL7;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the daily soil temperature data";
      rflog1 << "for soil layer 7";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYTSOIL) " << endl;
      rflog1 << idayTsoilfnameL7 << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << idayTsoilendL7 << endl << endl;


      cout << "Please enter the first part of the file name ";
      cout << "containing the daily soil temperature data";
      cout << "for soil layer 8";
      cout << endl;
      cout << "               (e.g., DAYTSOIL) " << endl;
    
      cin >> idayTsoilfnameL8;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> idayTsoilendL8;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the daily soil temperature data";
      rflog1 << "for soil layer 8";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYTSOIL) " << endl;
      rflog1 << idayTsoilfnameL8 << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << idayTsoilendL8 << endl << endl;


      cout << "Please enter the first part of the file name ";
      cout << "containing the daily soil temperature data";
      cout << "for soil layer 9";
      cout << endl;
      cout << "               (e.g., DAYTSOIL) " << endl;
    
      cin >> idayTsoilfnameL9;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> idayTsoilendL9;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the daily soil temperature data";
      rflog1 << "for soil layer 9";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYTSOIL) " << endl;
      rflog1 << idayTsoilfnameL9 << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
     rflog1 << idayTsoilendL9 << endl << endl;
    }
    else
    {
      cout << "Please enter the file name ";
      cout << "containing the daily soil temperature data ";
      cout << "for soil layer 1";
      cout << endl;
      cout << "               (e.g., DAYTSOIL) " << endl;
    
      cin >> idayTsoilfnameL1;

      rflog1 << "Please enter the file name ";
      rflog1 << "containing the daily soil temperature data";
      rflog1 << "for soil layer 1";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYTSOIL) " << endl;
      rflog1 << idayTsoilfnameL1 << endl << endl;


      cout << "Please enter the file name ";
      cout << "containing the daily soil temperature data";
      cout << "for soil layer 2";
      cout << endl;
      cout << "               (e.g., DAYTSOIL) " << endl;
    
      cin >> idayTsoilfnameL2;

      rflog1 << "Please enter the file name ";
      rflog1 << "containing the daily soil temperature data";
      rflog1 << "for soil layer 2";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYTSOIL) " << endl;
      rflog1 << idayTsoilfnameL2 << endl << endl;


      cout << "Please enter the file name ";
      cout << "containing the daily soil temperature data";
      cout << "for soil layer 3";
      cout << endl;
      cout << "               (e.g., DAYTSOIL) " << endl;
    
      cin >> idayTsoilfnameL3;

      rflog1 << "Please enter the file name ";
      rflog1 << "containing the daily soil temperature data";
      rflog1 << "for soil layer 3";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYTSOIL) " << endl;
      rflog1 << idayTsoilfnameL3 << endl << endl;


      cout << "Please enter the file name ";
      cout << "containing the daily soil temperature data";
      cout << "for soil layer 4";
      cout << endl;
      cout << "               (e.g., DAYTSOIL) " << endl;
    
      cin >> idayTsoilfnameL4;

      rflog1 << "Please enter the file name ";
      rflog1 << "containing the daily soil temperature data";
      rflog1 << "for soil layer 4";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYTSOIL) " << endl;
      rflog1 << idayTsoilfnameL4 << endl << endl;


      cout << "Please enter the file name ";
      cout << "containing the daily soil temperature data";
      cout << "for soil layer 5";
      cout << endl;
      cout << "               (e.g., DAYTSOIL) " << endl;
    
      cin >> idayTsoilfnameL5;

      rflog1 << "Please enter the file name ";
      rflog1 << "containing the daily soil temperature data";
      rflog1 << "for soil layer 5";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYTSOIL) " << endl;
      rflog1 << idayTsoilfnameL5 << endl << endl;


      cout << "Please enter the file name ";
      cout << "containing the daily soil temperature data";
      cout << "for soil layer 6";
      cout << endl;
      cout << "               (e.g., DAYTSOIL) " << endl;
    
      cin >> idayTsoilfnameL6;

      rflog1 << "Please enter the file name ";
      rflog1 << "containing the daily soil temperature data";
      rflog1 << "for soil layer 6";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYTSOIL) " << endl;
      rflog1 << idayTsoilfnameL6 << endl << endl;


      cout << "Please enter the file name ";
      cout << "containing the daily soil temperature data";
      cout << "for soil layer 7";
      cout << endl;
      cout << "               (e.g., DAYTSOIL) " << endl;
    
      cin >> idayTsoilfnameL7;

      rflog1 << "Please enter the file name ";
      rflog1 << "containing the daily soil temperature data";
      rflog1 << "for soil layer 7";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYTSOIL) " << endl;
      rflog1 << idayTsoilfnameL7 << endl << endl;


      cout << "Please enter the file name ";
      cout << "containing the daily soil temperature data";
      cout << "for soil layer 8";
      cout << endl;
      cout << "               (e.g., DAYTSOIL) " << endl;
    
      cin >> idayTsoilfnameL8;


      rflog1 << "Please enter the file name ";
      rflog1 << "containing the daily soil temperature data";
      rflog1 << "for soil layer 8";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYTSOIL) " << endl;
      rflog1 << idayTsoilfnameL8 << endl << endl;


      cout << "Please enter the file name ";
      cout << "containing the daily soil temperature data";
      cout << "for soil layer 9";
      cout << endl;
      cout << "               (e.g., DAYTSOIL) " << endl;
    
      cin >> idayTsoilfnameL9;

      rflog1 << "Please enter the file name ";
      rflog1 << "containing the daily soil temperature data";
      rflog1 << "for soil layer 9";
      rflog1 << endl;
      rflog1 << "               (e.g., DAYTSOIL) " << endl;
      rflog1 << idayTsoilfnameL9 << endl << endl;
    }
  }

};
#endif

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

#ifdef STANDALONE_TEM
void MITclm44::initDRAINAGE ( ofstream& rflog1 )
{

  if( 1 == tprecflag )
  {
    cout << "Please enter the first part of the file name ";
    cout << "containing the monthly drainage data: " << endl;
    cout << "            (e.g., DRAINAGE) " << endl;
    
    cin >> idrainfname;

    cout << "Please enter the file extension (include the '.'): ";
    
    cin >> idrainend;

    rflog1 << "Please enter the first part of the file name ";
    rflog1 << "containing the monthly drainage data: " << endl;
    rflog1 << "          (e.g., DRAINAGE) " << endl;
    rflog1 << idrainfname << endl << endl;
    rflog1 << "Please enter the file extension (include the '.'): ";
    rflog1 << idrainend << endl << endl;
  }
  else
  {
    cout << "Please enter the file name containing the ";
    cout << "long-term mean monthly drainage data: ";
    cout << endl;
    cout << "             (e.g., DRAINAGE) " << endl;

    cin >> idrainfname;

    rflog1 << "Please enter the file name containing the ";
    rflog1 << "long-term mean monthly drainage data: ";
    rflog1 << endl;
    rflog1 << "           (e.g., DRAINAGE) " << endl;
    rflog1 << idrainfname << endl << endl;
  }

};
#endif

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

#ifdef STANDALONE_TEM
void MITclm44::initEET ( ofstream& rflog1 )
{

  if( 1 == tprecflag )
  {
    cout << "Please enter the first part of the file name ";
    cout << "containing the monthly estimated ";
    cout << "evapotranspiration data: " << endl;
    cout << "               (e.g., EET) " << endl;
    
    cin >> ieetfname;

    cout << "Please enter the file extension (include the '.'): ";
    
    cin >> ieetend;

    rflog1 << "Please enter the first part of the file name ";
    rflog1 << "containing the monthly estimated ";
    rflog1 << "evapotranspiration data: " << endl;
    rflog1 << "               (e.g., EET) " << endl;
    rflog1 << ieetfname << endl << endl;
    rflog1 << "Please enter the file extension (include the '.'): ";
    rflog1 << ieetend << endl << endl;
  }
  else
  {
    cout << "Please enter the file name containing the ";
    cout << "long-term mean monthly estimated ";
    cout << "evapotranspiration data: ";
    cout << endl;
    cout << "               (e.g., EET) " << endl;

    cin >> ieetfname;

    rflog1 << "Please enter the file name containing the ";
    rflog1 << "long-term mean monthly estimated ";
    rflog1 << "evapotranspiration data: ";
    rflog1 << endl;
    rflog1 << "               (e.g., EET) " << endl;
    rflog1 << ieetfname << endl << endl;
  }

};
#endif

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

#ifdef STANDALONE_TEM
void MITclm44::initHourSoilMoist( ofstream& rflog1 )
{

  if( 1 == n2oflag )
  {
    if( 1 == tprecflag )
    {
      cout << "Please enter the first part of the file name ";
      cout << "containing the hourly soil moisture data";
      cout << "for soil layer 1";
      cout << endl;
      cout << "               (e.g., HRSM) " << endl;
    
      cin >> ihrSMfnameL1;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> ihrSMendL1;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the hourly soil moisture data";
      rflog1 << "for soil layer 1";
      rflog1 << endl;
      rflog1 << "               (e.g., HRSM) " << endl;
      rflog1 << ihrSMfnameL1 << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << ihrSMendL1 << endl << endl;

      cout << "Please enter the first part of the file name ";
      cout << "containing the hourly soil moisture data";
      cout << "for soil layer 2";
      cout << endl;
      cout << "               (e.g., HRSM) " << endl;
    
      cin >> ihrSMfnameL2;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> ihrSMendL2;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the hourly soil moisture data";
      rflog1 << "for soil layer 2";
      rflog1 << endl;
      rflog1 << "               (e.g., HRSM) " << endl;
      rflog1 << ihrSMfnameL2 << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << ihrSMendL2 << endl << endl;

      cout << "Please enter the first part of the file name ";
      cout << "containing the hourly soil moisture data";
      cout << "for soil layer 3";
      cout << endl;
      cout << "               (e.g., HRSM) " << endl;
    
      cin >> ihrSMfnameL3;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> ihrSMendL3;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the hourly soil moisture data";
      rflog1 << "for soil layer 3";
      rflog1 << endl;
      rflog1 << "               (e.g., HRSM) " << endl;
      rflog1 << ihrSMfnameL3 << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << ihrSMendL3 << endl << endl;

      cout << "Please enter the first part of the file name ";
      cout << "containing the hourly soil moisture data";
      cout << "for soil layer 4";
      cout << endl;
      cout << "               (e.g., HRSM) " << endl;
    
      cin >> ihrSMfnameL4;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> ihrSMendL4;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the hourly soil moisture data";
      rflog1 << "for soil layer 4";
      rflog1 << endl;
      rflog1 << "               (e.g., HRSM) " << endl;
      rflog1 << ihrSMfnameL4 << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << ihrSMendL4 << endl << endl;

      cout << "Please enter the first part of the file name ";
      cout << "containing the hourly soil moisture data";
      cout << "for soil layer 5";
      cout << endl;
      cout << "               (e.g., HRSM) " << endl;
    
      cin >> ihrSMfnameL5;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> ihrSMendL5;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the hourly soil moisture data";
      rflog1 << "for soil layer 5";
      rflog1 << endl;
      rflog1 << "               (e.g., HRSM) " << endl;
      rflog1 << ihrSMfnameL5 << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << ihrSMendL5 << endl << endl;

      cout << "Please enter the first part of the file name ";
      cout << "containing the hourly soil moisture data";
      cout << "for soil layer 6";
      cout << endl;
      cout << "               (e.g., HRSM) " << endl;
    
      cin >> ihrSMfnameL6;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> ihrSMendL6;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the hourly soil moisture data";
      rflog1 << "for soil layer 6";
      rflog1 << endl;
      rflog1 << "               (e.g., HRSM) " << endl;
      rflog1 << ihrSMfnameL6 << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << ihrSMendL6 << endl << endl;
    }
    else
    {
      cout << "Please enter the file name ";
      cout << "containing the hourly soil moisture data";
      cout << "for soil layer 1";
      cout << endl;
      cout << "               (e.g., HRSM) " << endl;
    
      cin >> ihrSMfnameL1;

      rflog1 << "Please enter the file name ";
      rflog1 << "containing the hourly soil moisture data";
      rflog1 << "for soil layer 1";
      rflog1 << endl;
      rflog1 << "               (e.g., HRSM) " << endl;
      rflog1 << ihrSMfnameL1 << endl << endl;

      cout << "Please enter the file name ";
      cout << "containing the hourly soil moisture data";
      cout << "for soil layer 2";
      cout << endl;
      cout << "               (e.g., HRSM) " << endl;
    
      cin >> ihrSMfnameL2;

      rflog1 << "Please enter the file name ";
      rflog1 << "containing the hourly soil moisture data";
      rflog1 << "for soil layer 2";
      rflog1 << endl;
      rflog1 << "               (e.g., HRSM) " << endl;
      rflog1 << ihrSMfnameL2 << endl << endl;
 
      cout << "Please enter the file name ";
      cout << "containing the hourly soil moisture data";
      cout << "for soil layer 3";
      cout << endl;
      cout << "               (e.g., HRSM) " << endl;
    
      cin >> ihrSMfnameL3;

      rflog1 << "Please enter the file name ";
      rflog1 << "containing the hourly soil moisture data";
      rflog1 << "for soil layer 3";
      rflog1 << endl;
      rflog1 << "               (e.g., HRSM) " << endl;
      rflog1 << ihrSMfnameL3 << endl << endl;

      cout << "Please enter the file name ";
      cout << "containing the hourly soil moisture data";
      cout << "for soil layer 4";
      cout << endl;
      cout << "               (e.g., HRSM) " << endl;
    
      cin >> ihrSMfnameL4;

      rflog1 << "Please enter the file name ";
      rflog1 << "containing the hourly soil moisture data";
      rflog1 << "for soil layer 4";
      rflog1 << endl;
      rflog1 << "               (e.g., HRSM) " << endl;
      rflog1 << ihrSMfnameL4 << endl << endl;

      cout << "Please enter the file name ";
      cout << "containing the hourly soil moisture data";
      cout << "for soil layer 5";
      cout << endl;
      cout << "               (e.g., HRSM) " << endl;
    
      cin >> ihrSMfnameL5;

      rflog1 << "Please enter the file name ";
      rflog1 << "containing the hourly soil moisture data";
      rflog1 << "for soil layer 5";
      rflog1 << endl;
      rflog1 << "               (e.g., HRSM) " << endl;
      rflog1 << ihrSMfnameL5 << endl << endl;

      cout << "Please enter the file name ";
      cout << "containing the hourly soil moisture data";
      cout << "for soil layer 6";
      cout << endl;
      cout << "               (e.g., HRSM) " << endl;
    
      cin >> ihrSMfnameL6;

      rflog1 << "Please enter the file name ";
      rflog1 << "containing the hourly soil moisture data";
      rflog1 << "for soil layer 6";
      rflog1 << endl;
      rflog1 << "               (e.g., HRSM) " << endl;
      rflog1 << ihrSMfnameL6 << endl << endl;
    }
  }

};
#endif

/* **************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void MITclm44::initMITCflux( ofstream& rflog1 )
{

  cout << "Please enter the first part of the file name";
  cout << " containing the TEM monthly net carbon flux data: ";
  cout << endl;
  cout << "               (e.g., NEP, CFLX) " << endl;

  cin >> icflxfname;

  cout << "Please enter the file extension (include the '.'): ";

  cin >> icflxend;

  rflog1 << "Please enter the first part of the file name";
  rflog1 << " containing the TEM monthly net carbon flux data: ";
  rflog1 << endl;
  rflog1 << "               (e.g., NEP, CFLX) " << endl;
  rflog1 << icflxfname << endl << endl;
  rflog1 << "Please enter the file extension (include the '.'): ";
  rflog1 << icflxend << endl << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITclm44::initMITCH4flux( ofstream& rflog1 )
{

  if( 1 == ch4flag )
  {  
    cout << "Please enter the first part of the file name";
    cout << " containing the TEM monthly methane flux data: ";
    cout << endl;
    cout << "               (e.g., CH4FLX) " << endl;

    cin >> ich4flxfname;

    cout << "Please enter the file extension (include the '.'): ";

    cin >> ich4flxend;

    rflog1 << "Please enter the first part of the file name";
    rflog1 << " containing the TEM monthly methane flux data: ";
    rflog1 << endl;
    rflog1 << "               (e.g., CH4FLX) " << endl;
    rflog1 << ich4flxfname << endl << endl;
    rflog1 << "Please enter the file extension (include the '.'): ";
    rflog1 << ich4flxend << endl << endl;

    #ifdef STANDALONE_TEM
      initDayTsoil( rflog1 );

      initDaySoilMoist( rflog1 );
    #endif
  }

    
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITclm44::initMITN2Oflux( ofstream& rflog1 )
{
  if( 1 == n2oflag )
  {  
    cout << "Please enter the first part of the file name";
    cout << " containing the TEM monthly nitrous oxide flux data: ";
    cout << endl;
    cout << "               (e.g., N2OFLX) " << endl;

    cin >> in2oflxfname;

    cout << "Please enter the file extension (include the '.'): ";

    cin >> in2oflxend;

    rflog1 << "Please enter the first part of the file name";
    rflog1 << " containing the TEM monthly nitrous oxide flux data: ";
    rflog1 << endl;
    rflog1 << "               (e.g., N2OFLX) " << endl;
    rflog1 << in2oflxfname << endl << endl;
    rflog1 << "Please enter the file extension (include the '.'): ";
    rflog1 << in2oflxend << endl << endl;

    #ifdef STANDALONE_TEM

      if ( 0 == ch4flag ) 
      {
        initDayTsoil( rflog1 );

        initDaySoilMoist( rflog1 );
      }
    
      initDayTair( rflog1 );

      initHourSoilMoist( rflog1 );

      initRainDuration( rflog1 );

      initRainIntensity( rflog1 );    // end comment here
    #endif
  } 

};

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */

#ifdef STANDALONE_TEM
void MITclm44::initMITSolarRad( ofstream& rflog1 )
{

  if( 1 == tcldsflag )
  {
    cout << "Please enter the first part of the file name";
    cout << " containing the monthly surface solar";
    cout << " radiation data: " << endl;
    cout << "               (e.g., NIRR) " << endl;

    rflog1 << "Please enter the first part of the file name";
    rflog1 << " containing the monthly surface solar";
    rflog1 << " radiation data: " << endl;
    rflog1 << "               (e.g., NIRR) " << endl;

    cin >> inirrfname;

    rflog1 << inirrfname << endl << endl;

    cout << "Please enter the file extension (include the '.'): ";
        
    cin >> inirrend;

    rflog1 << "Please enter the file extension (include the '.'): ";
    rflog1 << inirrend << endl << endl; 
  }
  else
  {
    cout << "Please enter the file name containing the";
    cout << " long-term mean monthly surface solar";
    cout << " radiation data: " << endl;
    cout << "               (e.g., NIRR) " << endl;

    rflog1 << "Please enter the file name containing the";
    rflog1 << " long-term monthly surface solar";
    rflog1 << " radiation data: " << endl;
    rflog1 << "               (e.g., NIRR) " << endl;

    cin >> inirrfname;
    
    rflog1 << inirrfname << endl << endl;
  }

};
#endif

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

#ifdef STANDALONE_TEM
void MITclm44::initPET( ofstream& rflog1 )
{

  if( 1 == tprecflag )
  {
    cout << "Please enter the first part of the file name ";
    cout << "containing the monthly potential ";
    cout << "evapotranspiration data: " << endl;
    cout << "               (e.g., PET) " << endl;
    
    cin >> ipetfname;

    cout << "Please enter the file extension (include the '.'): ";
    
    cin >> ipetend;

    rflog1 << "Please enter the first part of the file name ";
    rflog1 << "containing the monthly potential ";
    rflog1 << "evapotranspiration data: " << endl;
    rflog1 << "               (e.g., PET) " << endl;
    rflog1 << ipetfname << endl << endl;
    rflog1 << "Please enter the file extension (include the '.'): ";
    rflog1 << ipetend << endl << endl;
  }
  else
  {
    cout << "Please enter the file name containing the ";
    cout << "long-term mean monthly potential ";
    cout << "evapotranspiration data: ";
    cout << endl;
    cout << "               (e.g., PET) " << endl;

    cin >> ipetfname;

    rflog1 << "Please enter the file name containing the ";
    rflog1 << "long-term mean monthly potential ";
    rflog1 << "evapotranspiration data: ";
    rflog1 << endl;
    rflog1 << "               (e.g., PET) " << endl;
    rflog1 << ipetfname << endl << endl;
  }

};
#endif

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

#ifdef STANDALONE_TEM
void MITclm44::initRainDuration( ofstream& rflog1 )
{

  if( 1 == n2oflag )
  {
    if( 1 == tprecflag )
    {
      cout << "Please enter the first part of the file name ";
      cout << "containing the daily rain duration data";
      cout << endl;
      cout << "               (e.g., RAINEVNTS) " << endl;
    
      cin >> irainDurfname;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> irainDurend;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the daily rain duration data";
      rflog1 << endl;
      rflog1 << "               (e.g., RAINEVNTS) " << endl;
      rflog1 << irainDurfname << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << irainDurend << endl << endl;
    }
    else
    {
      cout << "Please enter the file name ";
      cout << "containing the daily rain duration data";
      cout << endl;
      cout << "               (e.g., RAINEVNTS) " << endl;
    
      cin >> irainDurfname;

      rflog1 << "Please enter the file name ";
      rflog1 << "containing the daily rain duration data";
      rflog1 << endl;
      rflog1 << "               (e.g., RAINEVNTS) " << endl;
      rflog1 << irainDurfname << endl << endl;
    }
  }

};
#endif

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

#ifdef STANDALONE_TEM
void MITclm44::initRainIntensity( ofstream& rflog1 )
{

  if( 1 == n2oflag )
  {
    if( 1 == tprecflag )
    {
      cout << "Please enter the first part of the file name ";
      cout << "containing the daily rain intensity data";
      cout << endl;
      cout << "               (e.g., RAININT) " << endl;
    
      cin >> irainIntfname;

      cout << "Please enter the file extension (include the '.'): ";
    
      cin >> irainIntend;

      rflog1 << "Please enter the first part of the file name ";
      rflog1 << "containing the daily rain intensity data";
      rflog1 << endl;
      rflog1 << "               (e.g., RAININT) " << endl;
      rflog1 << irainIntfname << endl << endl;
      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << irainIntend << endl << endl;
    }
    else
    {
      cout << "Please enter the file name ";
      cout << "containing the daily rain intensity data";
      cout << endl;
      cout << "               (e.g., RAININT) " << endl;
    
      cin >> irainIntfname;

      rflog1 << "Please enter the file name ";
      rflog1 << "containing the daily rain intensity data";
      rflog1 << endl;
      rflog1 << "               (e.g., RAININT) " << endl;
      rflog1 << irainIntfname << endl << endl;
    }

  }

};
#endif

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

#ifdef STANDALONE_TEM
void MITclm44::initSNOWPACK ( ofstream& rflog1 )
{

  if( 1 == tprecflag )
  {
    cout << "Please enter the first part of the file name ";
    cout << "containing the monthly snowpack data: " << endl;
    cout << "            (e.g., SNOWPACK) " << endl;
    
    cin >> ispackfname;

    cout << "Please enter the file extension (include the '.'): ";
    
    cin >> ispackend;

    rflog1 << "Please enter the first part of the file name ";
    rflog1 << "containing the monthly snowpack data: " << endl;
    rflog1 << "           (e.g., SNOWPACK) " << endl;
    rflog1 << ispackfname << endl << endl;
    rflog1 << "Please enter the file extension (include the '.'): ";
    rflog1 << ispackend << endl << endl;
  }
  else
  {
    cout << "Please enter the file name containing the ";
    cout << "long-term mean monthly snowpack data: ";
    cout << endl;
    cout << "           (e.g., SNOWPACK) " << endl;

    cin >> ispackfname;

    rflog1 << "Please enter the file name containing the ";
    rflog1 << "long-term mean monthly snowpack data: ";
    rflog1 << endl;
    rflog1 << "         (e.g., SNOWPACK) " << endl;
    rflog1 << ispackfname << endl << endl;
  }

};
#endif

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

#ifdef STANDALONE_TEM
void MITclm44::initSoilH2O( ofstream& rflog1 )
{

  if( 1 == tprecflag )
  {
    cout << "Please enter the first part of the file name ";
    cout << "containing the monthly soil moisture data: " << endl;
    cout << "for the top 1 meter: " << endl; 
    cout << "            (e.g., SOILH2O1M) " << endl;
    
    cin >> ish2o1fname;

    cout << "Please enter the file extension (include the '.'): ";
    
    cin >> ish2o1end;

    rflog1 << "Please enter the first part of the file name ";
    rflog1 << "containing the monthly soil moisture data: " << endl;
    rflog1 << "for the top 1 meter: " << endl; 
    rflog1 << "            (e.g., SOILH2O1M) " << endl;
    rflog1 << ish2o1fname << endl << endl;
    rflog1 << "Please enter the file extension (include the '.'): ";
    rflog1 << ish2o1end << endl << endl;

    cout << "Please enter the first part of the file name ";
    cout << "containing the monthly soil moisture data: " << endl;
    cout << "for the top 2 meters: " << endl; 
    cout << "            (e.g., SOILH2O2M) " << endl;
    
    cin >> ish2o2fname;

    cout << "Please enter the file extension (include the '.'): ";
    
    cin >> ish2o2end;

    rflog1 << "Please enter the first part of the file name ";
    rflog1 << "containing the monthly soil moisture data: " << endl;
    rflog1 << "for the top 2 meter: " << endl; 
    rflog1 << "            (e.g., SOILH2O2M) " << endl;
    rflog1 << ish2o2fname << endl << endl;
    rflog1 << "Please enter the file extension (include the '.'): ";
    rflog1 << ish2o2end << endl << endl;
  }
  else
  {
    cout << "Please enter the file name containing the ";
    cout << "long-term mean monthly soil moisture data: ";
    cout << endl;
    cout << "for the top 1 meter: " << endl; 
    cout << "            (e.g., SOILH2O1M) " << endl;

    cin >> ish2o1fname;

    rflog1 << "Please enter the file name containing the ";
    rflog1 << "long-term mean monthly soil moisture data: ";
    rflog1 << endl;
    rflog1 << "for the top 1 meter: " << endl; 
    rflog1 << "            (e.g., SOILH2O1M) " << endl;
    rflog1 << ish2o1fname << endl << endl;

    cout << "Please enter the file name containing the ";
    cout << "long-term mean monthly soil moisture data: ";
    cout << endl;
    cout << "for the top 2 meters: " << endl; 
    cout << "            (e.g., SOILH2O2M) " << endl;

    cin >> ish2o2fname;

    rflog1 << "Please enter the file name containing the ";
    rflog1 << "long-term mean monthly soil moisture data: ";
    rflog1 << endl;
    rflog1 << "for the top 2 meter: " << endl; 
    rflog1 << "            (e.g., SOILH2O2M) " << endl;
    rflog1 << ish2o2fname << endl << endl;
  }

};
#endif

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */

#ifdef STANDALONE_TEM
void MITclm44::initSURFRUN ( ofstream& rflog1 )
{

  if( 1 == tprecflag )
  {
    cout << "Please enter the first part of the file name ";
    cout << "containing the monthly surface runoff ";
    cout << "data: " << endl;
    cout << "            (e.g., SURFRUN) " << endl;
    
    cin >> isrfrunfname;

    cout << "Please enter the file extension (include the '.'): ";
    
    cin >> isrfrunend;

    rflog1 << "Please enter the first part of the file name ";
    rflog1 << "containing the monthly surface runoff ";
    rflog1 << "data: " << endl;
    rflog1 << "          (e.g., SURFRUN) " << endl;
    rflog1 << isrfrunfname << endl << endl;
    rflog1 << "Please enter the file extension (include the '.'): ";
    rflog1 << isrfrunend << endl << endl;
  }
  else
  {
    cout << "Please enter the file name containing the ";
    cout << "long-term mean monthly surface runoff ";
    cout << "data: ";
    cout << endl;
    cout << "            (e.g., SURFRUN) " << endl;

    cin >> isrfrunfname;

    rflog1 << "Please enter the file name containing the ";
    rflog1 << "long-term mean monthly surface runoff ";
    rflog1 << "data: ";
    rflog1 << endl;
    rflog1 << "          (e.g., SURFRUN) " << endl;
    rflog1 << isrfrunfname << endl << endl;
  }

};
#endif

/* **************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void MITclm44::mkHourlyOzone( const double ozone3hrs[MAX3HRS] )
{
  int hr;
  
  int t = 0;  

  for( hr = 0; hr < 24; ++hr )
  {
    if( 0 == ((hr+3) % 3) ) 
    {
      if( 21 == hr ) 
      {
        final[hr] = ozone3hrs[t];
        
	      final[hr+1] = ozone3hrs[t] 
	                    + (ozone3hrs[0] - ozone3hrs[t]) / 3.0;
	
	      final[hr+2] = ozone3hrs[t] 
	                    + 2 * (ozone3hrs[0] - ozone3hrs[t]) / 3.0;
      }
      else
      {
	      final[hr] = ozone3hrs[t];
	
	      final[hr+1] = ozone3hrs[t] 
	                    + (ozone3hrs[t+1] - ozone3hrs[t]) / 3.0;
	
	      final[hr+2] = ozone3hrs[t] 
	                    + 2 * (ozone3hrs[t+1] - ozone3hrs[t]) / 3.0;
	
	      ++t;
      }
    }
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITclm44::mkMITaot40( const double ozone3hrs[MAX3HRS], 
		           const int& pdm )
{
  double correct;

  double dum;

  int hr;

  double sumtmp;

//  correction factor based on ratio of MATCH results to Chien's model

  correct = 1.2;

  sumtmp = ZERO;


  // Interpolate 3-hourly ozone data to hourly data
  //   (i.e. values for final[MAXHR] )
    
  mkHourlyOzone( ozone3hrs );
  
  
  for( hr = 6; hr < 19; ++hr ) 
  {
    dum = final[hr];
    
    if( dum < 10 ) { dum = 10.0; }
    
    if( dum >= 40.0 ) 
    {
      sumtmp += (dum - 40.0);
    }
  }

  setAOT40( sumtmp * ndays[pdm] * correct ) ;


};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITclm44::resetCFLUX( const int& outmon )
{

  int dmlat;

  tcflx2D.year = -99;

  tcflx2D.mon = outmon;
  
  for( dmlat = 0; dmlat < MXMITNLAT; dmlat++ )
  {
    tcflx2D.latband[dmlat] = ZERO;
  }


};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITclm44::resetCH4FLUX( const int& outmon )
{

  int dmlat;

  tch4flx2D.year = -99;

  tch4flx2D.mon = outmon;
  
  for( dmlat = 0; dmlat < MXMITNLAT; dmlat++ )
  {
    tch4flx2D.latband[dmlat] = ZERO;
  }


};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITclm44::resetN2OFLUX( const int& outmon )
{

  int dmlat;

  tn2oflx2D.year = -99;

  tn2oflx2D.mon = outmon;
  
  for( dmlat = 0; dmlat < MXMITNLAT; dmlat++ )
  {
    tn2oflx2D.latband[dmlat] = ZERO;
  }


};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

#ifdef STANDALONE_TEM
void MITclm44::setLatBandClm( ifstream& ifile, 
                              MITCLMdata44& clm )
{

  clm.get( ifile );

  ifile.seekg( 0, ios::cur );

};
#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITclm44::setLatBandClm( ifstream& ifile, MITdata44& clm )
{

  clm.get( ifile );

  ifile.seekg( 0, ios::cur );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

#ifdef STANDALONE_TEM
void MITclm44::setLatBandClm( ifstream& ifile, 
                              const int& outmon,
                              MITDayChrtData44 clm[MAXMDAYS] )
{
  int dday;
  
  for( dday = 0; dday < mdays[outmon]; dday++ )
  {
    clm[dday].get( ifile );
  }

};
#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

#ifdef STANDALONE_TEM
void MITclm44::setLatBandClm( ifstream& ifile, 
                              const int& outmon,
                              MITDayClmData44 clm[MAXMDAYS] )
{
  int dday;
  
  for( dday = 0; dday < mdays[outmon]; dday++ )
  {
    clm[dday].get( ifile );
  }

};
#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

#ifdef STANDALONE_TEM
void MITclm44::setLatBando3( ifstream& ifile, O3data44& clm )
{

  clm.get( ifile );

  ifile.seekg( 0, ios::cur );

};
#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

#ifdef STANDALONE_TEM
void MITclm44::setLatBandSoil( ifstream& ifile,
                               const int& outmon, 
                               MITDayChrtData44 soilenv[MAXMDAYS] )
{
  int dday;
  
  for( dday = 0; dday < mdays[outmon]; dday++ )
  {
    soilenv[dday].get( ifile );
  }
};
#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

#ifdef STANDALONE_TEM
void MITclm44::setLatBandSoil( ifstream& ifile,
                               const int& outmon, 
                               MITHrClmData44 soilenv[MAXMDAYS][MAXDAYHRS] )
{
  int dday;
  int dhr;

  for( dday = 0; dday < mdays[outmon]; dday++ )
  {
    for (dhr = 0; dhr < MAXDAYHRS; dhr++ )
    {	
      soilenv[dday][dhr].get( ifile );
    }
  }

};
#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void MITclm44::setMITCH4Flags( ofstream& rflog1 )
{
  cout << "Do you want to generate CH4 fluxes?:" << endl;
  cout << "Enter 0 for No:" << endl;
  cout << "Enter 1 for Yes: ";
    
  cin >> ch4flag;

  rflog1 << "Do you want to generate CH4 fluxes?:" << endl;
  rflog1 << "Enter 0 for No:" << endl;
  rflog1 << "Enter 1 for Yes: ";
  rflog1 << "ch4flag = " << ch4flag << endl << endl;
  	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void MITclm44::setMITN2OFlags( ofstream& rflog1 )
{
  cout << "Do you want to generate N2O fluxes?:" << endl;
  cout << "Enter 0 for No:" << endl;
  cout << "Enter 1 for Yes: ";
    
  cin >> n2oflag;

  rflog1 << "Do you want to generate N2O fluxes?:" << endl;
  rflog1 << "Enter 0 for No:" << endl;
  rflog1 << "Enter 1 for Yes: ";
  rflog1 << "n2oflag = " << n2oflag << endl << endl;
  	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void MITclm44::setMITSolarRadFlags( ofstream& rflog1,
                                    const int& requil )
{

  tcldsflag = 0;

  if( 0 == requil )
  {
    cout << "Do you have transient solar radiation data?:";
    cout << endl;

    rflog1 << endl;
    rflog1 << "Do you have transient solar radiation data?:";
    rflog1 << endl;
  
    cout << "Enter 0 for no:" << endl;
    cout << "Enter 1 for yes: ";

    cin >> tcldsflag;

    rflog1 << "Enter 0 for no:" << endl;
    rflog1 << "Enter 1 for yes: " << endl;
    rflog1 << "tcldsflag = " << tcldsflag << endl << endl;
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITclm44::timecheck( ofstream& rflog1 )
{

/* *************************************************************
  Look for temporal co-registration problems among climate data
************************************************************* */

  // Cloudiness and gross irradiance spatially explicit data sets

  if( 0 == cldflag ) { cldsyr = nirryr; }

  if( 1 == tcldsflag && 0 == sradflag  
       && 1 == cldflag  && cldsyr != girryr ) 
  { 
    cout << "Cloud data is for year " << cldsyr << endl; 
    cout << "GIRR data is for year " << girryr << endl;
    rflog1 << "Cloud data is for year " << cldsyr << endl;
    rflog1 << "GIRR data is for year " << girryr << endl;
    
    exit( -1 ); 
  }

  // Net and gross irradiance spatially explicit data sets

  if( 1 == tcldsflag && 0 == sradflag 
      && 0 == cldflag && cldsyr != girryr ) 
  { 
    cout << "NIRR data is for year " << cldsyr << endl;
    cout << "GIRR data is for year " << girryr << endl;
    rflog1 << "Cloud data is for year " << cldsyr << endl;
    rflog1 << "GIRR data is for year " << girryr << endl;
    
    exit( -1 ); 
  }

  // NIRR and photosyntheically active radiation (PAR) spatially 
  //   explicit data sets

  if( 1 == tcldsflag && 1 == parflag && cldsyr != paryr ) 
  { 
    cout << "NIRR data is for year " << cldsyr << endl;
    cout << "PAR data is for year " << paryr << endl;
    rflog1 << "NIRR data is for year " << cldsyr << endl;
    rflog1 << "PAR data is for year " << paryr << endl;
    
    exit( -1 ); 
  }

  // Cloudiness and air temperature spatially explicit data sets

  if( 1 == tcldsflag && 1 == ttairflag && cldsyr != tairyr ) 
  { 
    if( 1 == cldflag ) 
    { 
      cout << "Cloud data is for year " << cldsyr << endl; 
    }
    else { cout << "NIRR data is for year " << cldsyr << endl; }
    
    cout << "Tair data is for year " << tairyr << endl;
    rflog1 << "Cloud data is for year " << cldsyr << endl;
    rflog1 << "Tair data is for year " << tairyr << endl;
    
    exit( -1 ); 
  }

  // Air temperature and preipitation spatially explicit data sets

  if( 1 == ttairflag && 1 == tprecflag && tairyr != precyr ) 
  { 
    cout << "Tair data is for year " << tairyr << endl;
    cout << "Prec data is for year " << precyr << endl;
    rflog1 << "Tair data is for year " << tairyr << endl;
    rflog1 << "Prec data is for year " << precyr << endl;
    
    exit( -1 ); 
  }

  // Cloudiness and precipitation spatially explicit data sets

  if( 1 == tcldsflag && 1 == tprecflag && cldsyr != precyr ) 
  { 
    if( 1 == cldflag ) 
    { 
      cout << "Cloud data is for year " << cldsyr << endl; 
    }
    else { cout << "NIRR data is for year " << cldsyr << endl; }
    
    cout << "Prec data is for year " << precyr << endl;
    rflog1 << "Cloud data is for year " << cldsyr << endl;
    rflog1 << "Prec data is for year " << precyr << endl;
    
    exit( -1 ); 
  }
    
  // Cloudiness and atmospheric CO2

  if( 1 == tcldsflag && 1 == tco2flag && cldsyr != co2yr ) 
  { 
    if( 1 == cldflag ) 
    { 
      cout << "Cloud data is for year " << cldsyr << endl; 
    }
    else { cout << "NIRR data is for year " << cldsyr << endl; }
    
    cout << "CO2 data is for year " << co2yr << endl;
    rflog1 << "Cloud data is for year " << cldsyr << endl;
    rflog1 << "CO2 data is for year " << co2yr << endl;
    
    exit( -1 ); 
  }

  // Air temperature and atmospheric CO2
    
  if( 1 == ttairflag && 1 == tco2flag && tairyr != co2yr ) 
  { 
    cout << "Tair data is for year " << tairyr << endl;
    cout << "CO2 data is for year " << co2yr << endl;
    rflog1 << "Tair data is for year " << tairyr << endl;
    rflog1 << "CO2 data is for year " << co2yr << endl;
    
    exit( -1 ); 
  }

  // Precipitation and atmospheric CO2
    
  if( 1 == tprecflag && 1 == tco2flag && precyr != co2yr ) 
  { 
    cout << "Prec data is for year " << precyr << endl;
    cout << "CO2 data is for year " << co2yr << endl;
    rflog1 << "Prec data is for year " << precyr << endl;
    rflog1 << "CO2 data is for year " << co2yr << endl;
    
    exit( -1 ); 
  }

  // Cloudiness and atmospheric O3

  if( 1 == tcldsflag && 1 == to3flag && cldsyr != o3yr ) 
  { 
    cout << "Cloud data is for year " << cldsyr << endl;
    cout << "O3 data is for year " << o3yr << endl;
    rflog1 << "Cloud data is for year " << cldsyr << endl;
    rflog1 << "O3 data is for year " << o3yr << endl;
    
    exit( -1 ); 
  }

  // Air temperature and atmospheric O3
    
  if( 1 == ttairflag && 1 == to3flag && tairyr != o3yr ) 
  { 
    cout << "Tair data is for year " << tairyr << endl;
    cout << "O3 data is for year " << o3yr << endl;
    rflog1 << "Tair data is for year " << tairyr << endl;
    rflog1 << "O3 data is for year " << o3yr << endl;
    
    exit( -1 ); 
  }

  // Precipitation and atmospheric O3
    
  if( 1 == tprecflag && 1 == to3flag && precyr != o3yr ) 
  { 
    cout << "Prec data is for year " << precyr << endl;
    cout << "O3 data is for year " << o3yr << endl;
    rflog1 << "Prec data is for year " << precyr << endl;
    rflog1 << "O3 data is for year " << o3yr << endl;
    
    exit( -1 ); 
  }

};

