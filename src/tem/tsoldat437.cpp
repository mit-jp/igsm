/* *************************************************************
TSOLDAT437.CPP - object to read and write the structure of soil
                    texture data from/to files used by the  
                    the Terrestrial Ecosystem Model (TEM)

Modifications:

20060114 - DWK created by modifying tsoldat425.cpp
20060114 - DWK added include tsoldat437.h and standard includes
20060114 - DWK changed Soildata:: to Soildata43::

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
  

#include "tsoldat437.h"



Soildata43::Soildata43( void )
{

  soilend = 1;
  lagpos = -99;
  curpos = 0;

};

/* *************************************************************
************************************************************* */

int Soildata43::get( ifstream& infile )
{

  lagpos = infile.tellg();

  infile >> col;
  infile >> row;
  infile >> varname;
  infile >> carea;
  infile >> pctsand;
  infile >> pctsilt;
  infile >> pctclay;
  infile >> wsoil;
  infile  >> source;
  infile >> contnent;

  infile.seekg( 0, ios::cur );
  
  curpos = infile.tellg();

  if( curpos < (lagpos + 10) ) { soilend = -1; }

  return soilend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Soildata43::getdel( FILE* infile )
{
  char tmpvarname[80];
  char tmpsource[80];
  char tmpcontnent[80];
  
  soilend = fscanf( infile,
                    "%f,%f, %s ,%ld,%lf,%lf,%lf,%d, %s , %s ",
                    &col,
                    &row,
                    tmpvarname,
                    &carea,
                    &pctsand,
                    &pctsilt,
                    &pctclay,
                    &wsoil,
                    tmpsource,
                    tmpcontnent );
  
  varname = tmpvarname;
  source = tmpsource;
  contnent = tmpcontnent;
  
  return soilend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Soildata43::out( ofstream& ofile, 
                      const float& col, 
                      const float& row, 
                      const string& varname,
                      const long& carea, 
                      const double& pctsand, 
                      const double& pctsilt, 
                      const double& pctclay,
                      const int& wsoil, 
                      const string& source, 
                      const string& contnent )
{

   ofile.setf( ios::fixed,ios::floatfield );
   ofile.setf( ios::showpoint );
   ofile.precision( 1 );

  ofile << col << ' ';
  ofile << row << ' ';
  ofile << varname << ' ';
  ofile << setprecision( 0 ) << carea << ' ';
  ofile << setprecision( 2 ) << pctsand << ' ';
  ofile << pctsilt << ' ';
  ofile << pctclay << ' ';
  ofile << setprecision( 0 ) << wsoil << ' ';
  ofile << source << ' ';
  ofile << contnent;
  ofile << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Soildata43::outdel( ofstream& ofile, 
                         const float& col, 
                         const float& row, 
                         const string& varname,
                         const long& carea, 
                         const double& pctsand, 
                         const double& pctsilt, 
                         const double& pctclay,
                         const int& wsoil, 
                         const string& source, 
                         const string& contnent )
{

   ofile.setf( ios::fixed,ios::floatfield );
   ofile.setf( ios::showpoint );
   ofile.precision( 1 );

  ofile << col << ",";
  ofile << row << ", ";
  ofile << varname << " ,";
  ofile << setprecision( 0 ) << carea << ",";
  ofile << setprecision( 2 ) << pctsand << ",";
  ofile << pctsilt << ",";
  ofile << pctclay << ",";
  ofile << setprecision( 0 ) << wsoil << ", ";
  ofile << source << " , ";
  ofile << contnent;
  ofile << endl;

};

