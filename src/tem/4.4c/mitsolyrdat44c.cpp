/* *************************************************************
MITSOLYRDAT44C.CPP - object to read and write the structure of 
                     soil layer characteristics data from/to 
                     files used by the Terrestrial Ecosystem 
                     Model (TEM)

Modifications:

20060302 - DWK created by modifying mitsoldat437.cpp
20080130 - DWK changed include from mitsolyrdat437.h to
           mitsolyrdat44a.h
20080130 - DWK changed MITSoilLayerdata43:: to 
           MITSoilLayerdata44::
20110707 - DWK changed include from mitsolyrdat44a.h to 
           mitsolyrdat44c.h
                       
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

#include "mitsolyrdat44c.h"

/* *************************************************************
************************************************************* */

MITSoilLayerdata44::MITSoilLayerdata44( void )
{

  soilend = 1;
  lagpos = -99;
  curpos = 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int MITSoilLayerdata44::get( ifstream& infile )
{
  
  lagpos = infile.tellg();

  infile >> col;
  infile >> row;
  infile >> varname;
  infile >> layer;  
  infile >> thickness;
  infile >> porosity;
  infile >> Ksat;
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

int MITSoilLayerdata44::getdel( FILE* infile )
{
  char tmpvarname[40];
  char tmpsource[40];
  char tmpregion[40];
	

  soilend = fscanf( infile,
                    "%f,%f, %s ,%d,%lf,%lf,%lf, %s , %s ",
                    &col,
                    &row,
                    tmpvarname,
                    &layer,
                    &thickness,
                    &porosity,
                    &Ksat,
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

void MITSoilLayerdata44::out( ofstream& ofile, 
                         const float& col, 
                         const float& row, 
                         const string& varname,
                         const int& layer, 
                         const double& thickness,
                         const double& porosity,
                         const double& Ksat, 
                         const string& source, 
                         const string& region )
{
  ofile.setf( ios::fixed,ios::floatfield );
  ofile.setf( ios::showpoint );
  ofile.precision( 1 );

  ofile << col << ' ';
  ofile << row << ' ';
  ofile << varname << ' ';
  ofile << layer << ' ';
  ofile << setprecision( 2 ) << thickness << ' ';
  ofile << setprecision( 3 ) << porosity << ' ';
  ofile << setprecision( 6 ) << Ksat << ' ';

  ofile << source << ' ';
  ofile << region;
  ofile << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MITSoilLayerdata44::outdel( ofstream& ofile, 
                                 const float& col, 
                                 const float& row, 
                                 const string& varname,
                                 const int& layer, 
                                 const double& thickness, 
                                 const double& porosity,
                                 const double& Ksat, 
                                 const string& source, 
                                 const string& region )
{

  ofile.setf( ios::fixed,ios::floatfield );
  ofile.setf( ios::showpoint );
  ofile.precision( 1 );

  ofile << col << ",";
  ofile << row << ", ";
  ofile << varname << " ,";
  ofile << layer << ",";

  ofile << setprecision( 2 ) << thickness << ",";
  ofile << setprecision( 3 ) << porosity << ' ';
  ofile << setprecision( 6 ) << Ksat << ' ';


  ofile << source << " , ";
  ofile << region;
  ofile << endl;

};

