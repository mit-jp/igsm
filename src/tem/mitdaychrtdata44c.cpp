/* *************************************************************
MITDAYCHRTDATA44C.CPP - object to read and write daily CO2,  
                        temperature, precipitation and solar 
                        radiation data from the MIT-IGSM model

Modifications:

20060128 - DWK created by modifying mitdaychrtdata436.cpp 
20060128 - DWK added include mitdaychrtdata437.h and standard 
           includes
20060128 - DWK changed char mname[4] to string mname in 
           writeclm()
20080130 - DWK changed include from mitdaychrtdata437.h to
           mitdaychrtdata44a.h
20080130 - DWK changed MITDayChrtData43:: to MITDayChrtData44::
20110707 - DWK changed include from mitdaychrtdata44a.h to
           mitdaychrtdata44c.h
                       
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
  

#include "mitdaychrtdata44c.h"

/* *************************************************************
************************************************************* */

MITDayChrtData44::MITDayChrtData44( void ) { };


/* *************************************************************
************************************************************* */

void MITDayChrtData44::get( ifstream& ifile )
{
  int dmlat;

  int iveg;

  ifile >> year;  

  ifile >> day;

  ifile >> varname;

  for( dmlat = 0; dmlat < MXMITNLAT; ++dmlat )
  {
    for( iveg = 0; iveg < CLMMXNVEG; ++iveg )
    {
      ifile >> latband[dmlat][iveg];
    }
  }

  ifile.seekg( 0, ios::cur );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITDayChrtData44::writeclm( ofstream& ofile, 
                                 string varname, 
                                 const int& dm, 
                                 const MITDayChrtData44& clm )
{
  int dmlat;

  int iveg;

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
  
  for( dmlat = 0; dmlat < MXMITNLAT; ++dmlat )
  {
    for( iveg = 0; iveg < CLMMXNVEG; ++iveg )
    {
      ofile << "  " << setprecision( 2 ) << setw( 8 ) << clm.latband[dmlat][iveg];
    }
    
    ofile << endl;
  }

};
