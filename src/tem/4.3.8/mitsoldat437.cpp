/* *************************************************************
MITSOLDAT437.CPP - object to read and write the structure of 
                   soil characteristics data from/to files used 
                   by the Terrestrial Ecosystem Model (TEM)

Modifications:

20060129 - DWK created by modifying mitsoldat436.cpp
20060129 - DWK added include mitsoldat437.h and standard 
           includes
20060129 - DWK changed public char region[MAXREGNAME] to 
           string region in functions
20060129 - DWK changed public char char source[9] to 
           string source in functions
20060129 - DWK changed public char char char varname[MAXVARNAME] 
           to string varname in functions
20060306 - DWK public deleted porosity[] amd Ksat[] from 
           functions
           
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

#include "mitsoldat437.h"

/* *************************************************************
************************************************************* */

MITSoildata43::MITSoildata43( void )
{

  soilend = 1;
  lagpos = -99;
  curpos = 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int MITSoildata43::get( ifstream& infile )
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
  infile >> pH;
    infile  >> source;
  infile >> region;

  infile.seekg( 0, ios::cur );
  curpos = infile.tellg();

  if ( curpos < (lagpos + 10) ) { soilend = -1; }

  return soilend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int MITSoildata43::getdel( FILE* infile )
{
  char tmpvarname[40];
  char tmpsource[40];
  char tmpregion[40];
	

  soilend = fscanf( infile,
                    "%f,%f, %s ,%ld,%lf,%lf,%lf,%d,%lf, %s , %s ",
                    &col,
                    &row,
                    tmpvarname,
                    &carea,
                    &pctsand,
                    &pctsilt,
                    &pctclay,
                    &wsoil,
                    &pH,
                    tmpsource,
                    tmpregion );

  varname = tmpvarname;
  source = tmpsource;
  region = tmpregion;
  
  return soilend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITSoildata43::out( ofstream& ofile, 
                         const float& col, 
                         const float& row, 
                         const string& varname,
                         const long& carea, 
                         const double& pctsand, 
                         const double& pctsilt, 
                         const double& pctclay,
                         const int& wsoil, 
                         const double& pH,
                         const string& source, 
                         const string& region )
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
  ofile << setprecision( 2 ) << pH << ' ';
  ofile << source << ' ';
  ofile << region;
  ofile << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITSoildata43::outdel( ofstream& ofile, 
                            const float& col, 
                            const float& row, 
                            const string& varname,
                            const long& carea, 
                            const double& pctsand, 
                            const double& pctsilt, 
                            const double& pctclay,
                            const int& wsoil,
                            const double& pH, 
                            const string& source, 
                            const string& region )
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
  ofile << setprecision( 0 ) << wsoil << ",";
  ofile << setprecision( 2 ) << pH << ",";
  ofile << source << " , ";
  ofile << region;
  ofile << endl;

};

