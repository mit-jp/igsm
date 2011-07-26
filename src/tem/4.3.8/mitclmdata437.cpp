/* *************************************************************
MITCLMDATA437.CPP - object to read and write CO2, cloudiness, 
                    temperature and precipitation data from the 
                   2-D MIT L-O climate model

Modifications:

20060128 - DWK created by modifying mitclmdata436.cpp 
20060128 - DWK added include mitclmdata437.h and standard 
           includes
20060128 - DWK changed char mname[4] to string mname in 
           get() and writeclm()
                      
****************************************************************
************************************************************* */

#include<cstdio>

  using std::fscanf;
  using std::FILE;

#include<iostream>

  using std::ios;
  using std::endl;

#include<fstream>

  using std::ifstream;
  using std::ofstream;

#include<iomanip>

  using std::setprecision;
  using std::setw;

#include<string>

  using std::string;
  

#include "mitclmdata437.h"


/* *************************************************************
************************************************************* */

MITCLMdata43::MITCLMdata43( void ) { };


/* *************************************************************
************************************************************* */

void MITCLMdata43::get( ifstream& ifile )
{

  int dmlat;
  int ichrt;
  string mname;

  ifile >> year;
  ifile >> mname;
  
  if ( "JAN" == mname ) { mon = 0; }
  if ( "FEB" == mname ) { mon = 1; }
  if ( "MAR" == mname ) { mon = 2; }
  if ( "APR" == mname ) { mon = 3; }
  if ( "MAY" == mname ) { mon = 4; }
  if ( "JUN" == mname ) { mon = 5; }
  if ( "JUL" == mname ) { mon = 6; }
  if ( "AUG" == mname ) { mon = 7; }
  if ( "SEP" == mname ) { mon = 8; }
  if ( "OCT" == mname ) { mon = 9; }
  if ( "NOV" == mname ) { mon = 10; }
  if ( "DEC" == mname ) { mon = 11; }
  
  ifile >> varname;

  for ( dmlat = 0; dmlat < MXMITNLAT; ++dmlat )
  {
    for ( ichrt = 0; ichrt < 19; ++ichrt )
    {
      ifile >> latband[dmlat][ichrt];
    }
  }

  ifile.seekg( 0, ios::cur );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITCLMdata43::writeclm( ofstream& ofile, 
                             string varname, 
                             const int& dm, 
                             const MITCLMdata43& clm )
{

  int dmlat;
  int ichrt;
  string mname;

  ofile.setf( ios::fixed,ios::floatfield );
  ofile.setf( ios::showpoint );
  ofile.precision( 0 );

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
  ofile << clm.year << "  ";
  ofile << mname << "  ";
  ofile << varname << endl;
  
  for ( dmlat = 0; dmlat < MXMITNLAT; ++dmlat )
  {
    for ( ichrt = 0; ichrt < 19; ++ichrt )
    {
      ofile << "  " << setprecision(2) << setw(8) << clm.latband[dmlat][ichrt];
    }
    
    ofile << endl;
  }

};
