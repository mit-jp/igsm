/* *************************************************************
MITDAYCHRTDATA437.CPP - object to read and write daily CO2,  
                        temperature, precipitation and solar 
                        radiation data from the MIT-IGSM model

Modifications:

20060128 - DWK created by modifying mitdaychrtdata436.cpp 
20060128 - DWK added include mitdaychrtdata437.h and standard 
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
  

#include "mitdaychrtdata437.h"

/* *************************************************************
************************************************************* */

MITDayChrtData43::MITDayChrtData43( void ) { };


/* *************************************************************
************************************************************* */

void MITDayChrtData43::get( ifstream& ifile )
{

  int dmlat;
  int ichrt;

  ifile >> year;  
  ifile >> day;
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

void MITDayChrtData43::writeclm( ofstream& ofile, 
                                 string varname, 
                                 const int& dm, 
                                 const MITDayChrtData43& clm )
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
  ofile << clm.day << "  ";
  ofile << varname << endl;
  
  for ( dmlat = 0; dmlat < MXMITNLAT; ++dmlat )
  {
    for ( ichrt = 0; ichrt < 19; ++ichrt )
    {
      ofile << "  " << setprecision( 2 ) << setw( 8 ) << clm.latband[dmlat][ichrt];
    }
    
    ofile << endl;
  }

};
