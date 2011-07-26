/* *************************************************************
O3DATA4376.CPP - object to read and write O3 data from the 
                  2-D MIT L-O climate model

20060128 - created by DWK by modifying o3data436.cpp 
20060128 - DWK added include o3data437.h and standard includes

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
  

#include "o3data437.h"



/* *************************************************************
************************************************************* */

O3data43::O3data43( void ) { };

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void O3data43::get( ifstream& ifile )
{

  int dmlat;
  int dmhr;

  ifile >> year;
  ifile >> mon;
  ifile >> varname;
 
  for ( dmlat = 0; dmlat < MXMITNLAT; ++dmlat )
  {
    for ( dmhr = 0; dmhr < MAX3HRS; ++dmhr )
    {
      ifile >> latbandhr[dmlat][dmhr];
    }
  }

  ifile.seekg( 0, ios::cur );

};
