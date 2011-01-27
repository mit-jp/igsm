/* *************************************************************
MITDAYCLMDATA437.CPP - object to read and write daily CO2,  
                       temperature, precipitation and solar 
                       radiation data from the MIT-IGSM model

Modifications:

20060128 - DWK created by modifying mitdayclmdata436.cpp 
20060128 - DWK added include mitclmdata437.h and standard 
           includes
20060128 - DWK changed char mname[4] to string mname in 
           writeclm()
           
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
  

#include "mitdayclmdata437.h"

/* *************************************************************
************************************************************* */

MITDayClmData43::MITDayClmData43( void ) { };


/* *************************************************************
************************************************************* */

void MITDayClmData43::get( ifstream& ifile )
{

  int dmlat;

  ifile >> year;  
  ifile >> day;
  ifile >> varname;

  for ( dmlat = 0; dmlat < MXMITNLAT; ++dmlat )
  {
    ifile >> latband[dmlat];
  }

  ifile.seekg( 0, ios::cur );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITDayClmData43::writeclm( ofstream& ofile, 
                                string varname, 
                                const int& dm, 
                                const MITDayClmData43& clm )
{

  int dmlat;
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
  ofile << clm.day << "  ";
  ofile << varname << endl;

  for ( dmlat = 0; dmlat < (MXMITNLAT/2); ++dmlat )
  {
    ofile << "  " << setprecision( 2 ) << setw( 8 ) << clm.latband[dmlat];
  }
  ofile << endl;

  for ( dmlat = (MXMITNLAT/2); dmlat < MXMITNLAT; ++dmlat )
  {
    ofile << "  " << setprecision( 2 ) << setw( 8 ) << clm.latband[dmlat];
  }
  ofile << endl;

};
