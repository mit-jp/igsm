/* *************************************************************
LANDAREADAT44C.CPP - object to read and write the structure 
                     of area and land fraction of a grid cell 
                     from/to files used by the 
                     Terrestrial Ecosystem Model

Modifications:

20110715 - DWK created by modifying tcohortdat437.cpp
20110715 - DWK changed include from tmaxcohortdat437.h to
           landareadat44c.h
20110715 - DWK changed MaxCohortdata43:: to LandAreadata44::


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

#include "landareadat44c.h" 

/* *************************************************************
************************************************************* */

LandAreadata44::LandAreadata44( void )
{

  chrtend = 1;
  lagpos = -99;
  curpos = 0;

};

/* *************************************************************
************************************************************* */

int LandAreadata44::get( ifstream& infile )
{

  lagpos = infile.tellg();

  infile >> col;
  infile >> row;
  infile >> varname;
  infile >> elemntArea;
  infile >> landFrac;
  infile >> contnent;

  infile.seekg( 0, ios::cur );
  
  curpos = infile.tellg();

  if( curpos < (lagpos + 10) ) { chrtend = -1; }

  return chrtend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int LandAreadata44::getdel( FILE* infile )
{
  char tmpvarname[40];
  char tmpcontnent[40];
  
  chrtend = fscanf( infile,"%f,%f, %s ,%f,%f, %s",
                   &col, 
                   &row, 
                   tmpvarname, 
                   &elemntArea,
                   &landFrac,  
                   tmpcontnent );

  varname = tmpvarname;
  contnent = tmpcontnent;
  
  return chrtend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void LandAreadata44::out( ofstream& ofile, 
                          const float& col, 
                          const float& row, 
                          const string& varname, 
                          const float& elemntArea,
                          const float& landFrac,  
                          const string& contnent )
{

  ofile.setf( ios::fixed,ios::floatfield );
  ofile.setf( ios::showpoint );
  ofile.precision( 1 );

  ofile << col << ' ';
  ofile << row << ' ';
  ofile << varname << ' ';
  ofile << setprecision( 8 ) << elemntArea << ' ';
  ofile << setprecision( 8 ) << landFrac << ' ';
  ofile << contnent;
  ofile << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void LandAreadata44::outdel( ofstream& ofile, 
                             const float& col, 
                             const float& row, 
                             const string& varname, 
                             const float& elmentArea,  
                             const float& landFrac,
                             const string& contnent )
{

  ofile.setf( ios::fixed,ios::floatfield );
  ofile.setf( ios::showpoint );
  ofile.precision( 1 );

  ofile << col << ",";
  ofile << row << ", ";
  ofile << varname << " ,";
  ofile << setprecision( 8 ) << elemntArea << ",";
  ofile << setprecision( 8 ) << landFrac << ", ";
  ofile << contnent;
  ofile << endl;

};

