/* *************************************************************
MITHRCLMDATA44A.CPP - object to read and write hourly 
                      temperature, precipitation and solar 
                      radiation files from the MIT-IGSM climate 
                      model

Modifications:

20060128 - DWK created by modifying mithrclmdata436.cpp 
20060128 - DWK added include mithrclmdata437.h and standard 
           includes
20060128 - DWK changed char mname[4] to string mname in 
           writeclm()
20060128 - DWK changed char hrname[4] to string hrname in get()
20060128 - DWK changed char dayname1[9] to string dayname2 in
           get()
20060128 - DWK changed char dayname2[9] to string dayname2 in 
           get()
20080130 - DWK changed include from mithrclmdata437.h to
           mithrclmdata44a.h
20080130 - DWK changed MITHrClmData43:: to MITHrClmData44::
                       
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
  

#include "mithrclmdata44a.h"

/* *************************************************************
************************************************************* */

MITHrClmData44::MITHrClmData44( void ) { };


/* *************************************************************
************************************************************* */

void MITHrClmData44::get( ifstream& ifile )
{
  string dayname1;

  string dayname2;

  int dmlat;

  int hour24;

  string hrname;

  int ichrt;
  
  ifile >> year;

  ifile >> hour24;

  ifile >> hrname;

  ifile >> dayname1;

  ifile >> dayname2;

  ifile >> day;

  ifile >> varname;
  
  for( dmlat = 0; dmlat < MXMITNLAT; ++dmlat )
  {
    for( ichrt = 0; ichrt < 19; ++ichrt )
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

void MITHrClmData44::writeclm( ofstream& ofile, 
                               string varname, 
                               const int& dm, 
                               const MITHrClmData44& clm )
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

  ofile << clm.day << " ";

  ofile << clm.hour << " ";

  ofile << varname << endl;
  
  for( dmlat = 0; dmlat < MXMITNLAT; ++dmlat )
  {
    for( ichrt = 0; ichrt < 19; ++ichrt )
    {
      ofile << "  " << setprecision( 2 ) << setw( 8 ) << clm.latband[dmlat][ichrt];
    }
    
    ofile << endl;
  }

};
