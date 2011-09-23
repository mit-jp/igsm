/* *************************************************************
O3DATA44C.CPP - object to read and write O3 data from the 
                  2-D MIT L-O climate model

20060128 - created by DWK by modifying o3data436.cpp 
20060128 - DWK added include o3data437.h and standard includes
20080130 - DWK changed include from o3data437.h to o3data44a.h
20080130 - DWK change O3data43:: to O3data44::
20110707 - DWK changed include from o3data44a.h to o3data44c.h

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
  

#include "o3data44c.h"



/* *************************************************************
************************************************************* */

O3data44::O3data44( void ) { };

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void O3data44::get( ifstream& ifile )
{
  int dmhr;

  int dmlat;

  ifile >> year;
 
  ifile >> mon;
 
  ifile >> varname;
 
  for( dmlat = 0; dmlat < MXMITNLAT; ++dmlat )
  {
    for( dmhr = 0; dmhr < MAX3HRS; ++dmhr )
    {
      ifile >> latbandhr[dmlat][dmhr];
    }
  }

  ifile.seekg( 0, ios::cur );

};
