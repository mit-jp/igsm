/* **************************************************************
LATDAT437.CPP - object to read and write the structure of 
                   latitude and longitude data from/to files

Modifications:

20060113 - DWK created by modifying latdat425.cpp
20060113 - DWK added include latdat437.h and standard includes
20060113 - DWK changed Latdata:: to Latdata43::
20060113 - DWK changed char varname[9] to string varname
20060113 - DWK changed char contnent[9] to string contnent

************************************************************** */

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

#include<string>

  using std::string;
  

#include "latdat437.h"

Latdata43::Latdata43( void )
{

  latend = 1;
  lagpos = -99;
  curpos = 0;

};

/* **************************************************************
                    Public Functions
************************************************************** */

/* *************************************************************
************************************************************* */

int Latdata43::get( ifstream& infile )
{
  lagpos = infile.tellg();

  infile >> col;

  infile >> row;

  infile >> varname;

  infile >> lat;

  infile >> lon;

  infile >> contnent;

  infile.seekg( 0, ios::cur );

  curpos = infile.tellg();

  if( curpos < (lagpos + 10) ) { latend = -1; }

  return latend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Latdata43::getdel( FILE* infile )
{
  char tmpvarname[80];

  char tmpcontnent[80];
  
  latend = fscanf( infile,
                   "%f,%f, %s ,%lf,%lf, %s",
                   &col,
                   &row,
                   tmpvarname,
                   &lat,
                   &lon,
                   tmpcontnent );

  varname = tmpvarname;

  contnent = tmpcontnent;
  
  return latend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Latdata43::out( ofstream& ofile, 
                     const float& col, 
                     const float& row, 
                     const string& varname,
                     const double& lat, 
                     const double& lon, 
                     const string& contnent )
{
  ofile.setf( ios::fixed,ios::floatfield );
  ofile.setf( ios::showpoint );
  ofile.precision( 1 );

  ofile << col << ' ';

  ofile << row << ' ';

  ofile << varname << ' ';

  ofile << setprecision( 7 ) << lat << ' ';

  ofile << lon << ' ';

  ofile << contnent;

  ofile << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Latdata43::outdel( ofstream& ofile, 
                        const float& col, 
                        const float& row, 
                        const string& varname,
                        const double& lat, 
                        const double& lon, 
                        const string& contnent )
{
  ofile.setf( ios::fixed,ios::floatfield );
  ofile.setf( ios::showpoint );
  ofile.precision( 1 );

  ofile << col << ",";

  ofile << row << ", ";

  ofile << varname << " ,";

  ofile << setprecision( 7 ) << lat << ",";

  ofile << lon << ", ";

  ofile << contnent;

  ofile << endl;

};

