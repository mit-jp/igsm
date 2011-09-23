/* *************************************************************
TELVDAT437.CPP - object to read and write the structure of 
                  elevation data from/to files used by the Water 
                  Balance Model

Modifications:

20060114 - DWK created by modifying telvdat425.cpp
20060114 - DWK added include telvdat437.h and standard includes
20060114 - DWK changed Elevdata:: to Elevdata43::

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

#include<string>

  using std::string;
  

#include "telvdat437.h"



Elevdata43::Elevdata43( void )
{

  elvend = 1;
  lagpos = -99;
  curpos = 0;

};

/* *************************************************************
************************************************************* */

int Elevdata43::get( ifstream& infile )
{

  lagpos = infile.tellg();

  infile >> col;
  infile >> row;
  infile >> varname;
  infile >> carea;
  infile >> elev;
  infile >> contnent;

  infile.seekg( 0, ios::cur );
  
  curpos = infile.tellg();

  if( curpos < (lagpos + 10) ) { elvend = -1; }

  return elvend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Elevdata43::getdel( FILE* infile )
{
  char tmpvarname[80];
  char tmpcontnent[80];
  
  elvend = fscanf( infile,
                   "%f,%f, %s ,%ld,%lf, %s",
                   &col,
                   &row,
                   tmpvarname,
                   &carea,
                   &elev,
                   tmpcontnent );

  varname = tmpvarname;
  contnent = tmpcontnent;

  return elvend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Elevdata43::out( ofstream& ofile, 
                      const float& col, 
                      const float& row, 
                      const string& varname,
                      const long& carea, 
                      const double& elev, 
                      const string& contnent )
{

  ofile.setf( ios::fixed,ios::floatfield );
  ofile.setf( ios::showpoint );
  ofile.precision( 1 );

  ofile << col << ' ';
  ofile << row << ' ';
  ofile << varname << ' ';
  ofile << setprecision( 0 ) << carea << ' ';
  ofile << setprecision( 1 ) << elev << ' ';
  ofile << contnent;
  ofile << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Elevdata43::outdel( ofstream& ofile, 
                         const float& col, 
                         const float& row, 
                         const string& varname,
                         const long& carea, 
                         const double& elev, 
                         const string& contnent )
{

  ofile.setf( ios::fixed,ios::floatfield );
  ofile.setf( ios::showpoint );
  ofile.precision( 1 );

  ofile << col << ",";
  ofile << row << ", ";
  ofile << varname << " ,";
  ofile << setprecision( 0 ) << carea << ",";
  ofile << setprecision( 1 ) << elev << ", ";
  ofile << contnent;
  ofile << endl;

};

