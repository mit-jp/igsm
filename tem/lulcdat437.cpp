/* *************************************************************
LULCDAT437.CPP - object to read and write the structure of 
                   land use/land cover data from/to files

Modifications:

20060113 - DWK created by modifying lulcdat425.cpp
20060113 - DWK changed Lulcdata:: to Lulcdata43::
20060113 - DWK added int currentveg and int potveg in functions
20060113 - DWK renamed long carea as long chrtarea in functions
20060113 - DWK renamed string contnent to string region 
           in functions
20060113 - DWK added int subtype to functions
20060113 - DWK added include lulcdat437.h and standard includes
20060218 - DWK added int isrccohort to functions

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

#include<string>

  using std::string;
  

#include "lulcdat437.h"


Lulcdata43::Lulcdata43( void )
{

  lulcend = 1;
  lagpos = -99;
  curpos = 0;

};

/* **************************************************************
                    Public Functions
************************************************************** */

/* *************************************************************
************************************************************* */

int Lulcdata43::get( ifstream& infile )
{

  lagpos = infile.tellg();

  infile >> col;
  infile >> row;
  infile >> varname;
  infile >> icohort;
  infile >> isrccohort;
  infile >> chrtarea;
  infile >> year;
  infile >> potveg;
  infile >> currentveg;
  infile >> subtype;
  infile >> FRF;
  infile >> agstate;
  infile >> agprevstate;
  infile >> pctag;
  infile >> tillflag;
  infile >> fertflag;
  infile >> irrgflag;
  infile >> RAP;
  infile >> slashpar;
  infile >> vconvert;
  infile >> prod10par;
  infile >> prod100par;
  infile >> vrespar;
  infile >> sconvert;
  infile >> region;

  infile.seekg( 0, ios::cur );
  curpos = infile.tellg();

  if( curpos < (lagpos + 10) ) { lulcend = -1; }

  return lulcend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Lulcdata43::getdel( FILE* infile )
{
  char tmpvarname[40];
  char tmpregion[80];

  lulcend = fscanf( infile, 
                    "%f,%f, %s ,%d,%d,%ld,%d,%d,%d,%d,%d,%d,%d,%lf,%d,%d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf, %s",
                    &col, 
                    &row, 
                    tmpvarname,
                    &icohort,
                    &isrccohort, 
                    &chrtarea, 
                    &year,
                    &potveg,
                    &currentveg,
                    &subtype,
                    &FRF, 
                    &agstate, 
                    &agprevstate, 
                    &pctag,
                    &tillflag,
                    &fertflag,
                    &irrgflag, 
                    &RAP, 
                    &slashpar, 
                    &vconvert,
                    &prod10par, 
                    &prod100par, 
                    &vrespar, 
                    &sconvert,
                    tmpregion );

  varname = tmpvarname;
  region = tmpregion;
  
  return lulcend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Lulcdata43::out( ofstream& ofile,  
                      const float& col, 
                      const float& row, 
                      const string& varname,
                      const int& icohort,
                      const int& isrccohort, 
                      const long& chrtarea, 
                      const int& year,
                      const int& potveg,
                      const int& currentveg, 
                      const int& subtype,
                      const int& FRF, 
                      const int& agstate, 
                      const int& agprevstate, 
                      const double& pctag, 
                      const int& tillflag, 
                      const int& fertflag, 
                      const int& irrgflag, 
                      const double& RAP, 
                      const double& slashpar, 
                      const double& vconvert,
                      const double& prod10par, 
                      const double& prod100par, 
                      const double& vrespar, 
                      const double& sconvert,
                      const string& region )
{

  ofile.setf( ios::fixed,ios::floatfield );
  ofile.setf( ios::showpoint );
  ofile.precision( 1 );

  ofile << col << ' ';
  ofile << row << ' ';
  ofile << varname << ' ';
  ofile << setprecision( 0 ) << icohort << ' ';
  ofile << setprecision( 0 ) << isrccohort << ' ';
  ofile << setprecision( 0 ) << chrtarea << ' ';
  ofile << year << ' ';
  ofile << setprecision( 0 ) << potveg << ' ';
  ofile << setprecision( 0 ) << currentveg << ' ';
  ofile << setprecision( 0 ) << subtype << ' ';
  ofile << FRF << ' ';
  ofile << agstate << ' ';
  ofile << agprevstate << ' ';
  ofile << setprecision( 3 ) << pctag << ' ';
  ofile << setprecision( 0 ) << tillflag << ' ';
  ofile << setprecision( 0 ) << fertflag << ' ';
  ofile << setprecision( 0 ) << irrgflag << ' ';
  ofile << setprecision( 2 ) << RAP << ' ';
  ofile << setprecision( 3 ) << slashpar << ' ';
  ofile << vconvert << ' ';
  ofile << prod10par << ' ';
  ofile << prod100par << ' ';
  ofile << vrespar << ' ';
  ofile << sconvert << ' ';
  ofile << region;
  ofile << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Lulcdata43::outdel( ofstream& ofile, 
                         const float& col, 
                         const float& row, 
                         const string& varname,
                         const int& icohort,
                         const int& isrccohort, 
                         const long& chrtarea, 
                         const int& year,
                         const int& potveg,
                         const int& currentveg,
                         const int& subtype, 
                         const int& FRF, 
                         const int& agstate, 
                         const int& agprevstate, 
                         const double& pctag, 
                         const int& tillflag, 
                         const int& fertflag, 
                         const int& irrgflag, 
                         const double& RAP, 
                         const double& slashpar, 
                         const double& vconvert,
                         const double& prod10par, 
                         const double& prod100par, 
                         const double& vrespar, 
                         const double& sconvert,
                         const string& region )
{

  ofile.setf( ios::fixed,ios::floatfield );
  ofile.setf( ios::showpoint );
  ofile.precision( 1 );

  ofile << col << ",";
  ofile << row << ", ";
  ofile << varname << " ,";
  ofile << icohort << ",";
  ofile << isrccohort << ",";
  ofile << chrtarea << ",";
  ofile << year << ",";
  ofile << potveg << ",";
  ofile << currentveg << ",";
  ofile << subtype << ",";
  ofile << FRF << ",";
  ofile << agstate << ",";
  ofile << agprevstate << ",";
  ofile << setprecision( 3 ) << pctag << ",";
  ofile << tillflag << ",";
  ofile << fertflag << ",";
  ofile << irrgflag << ",";
  ofile << setprecision( 2 ) << RAP << ",";
  ofile << setprecision( 3 ) << slashpar << ",";
  ofile << setprecision( 3 ) << vconvert << ",";
  ofile << setprecision( 3 ) << prod10par << ",";
  ofile << setprecision( 3 ) << prod100par << ",";
  ofile << setprecision( 3 ) << vrespar << ",";
  ofile << setprecision( 3 ) << sconvert << ", ";
  ofile << region;
  ofile << endl;

};

