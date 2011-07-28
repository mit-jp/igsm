/* *************************************************************
TMAXCOHORTDAT437.CPP - object to read and write the structure 
                          of maximum number of cohort grid cell 
                          data from/to files used by the 
                          Terrestrial Ecosystem Model

Modifications:

20060114 - DWK created by modifying tcohortdat433.cpp
20060114 - DWK changed Cohortdata43:: to MaxCohortdata43::
20060114 - DWK added int year to functions  
20060114 - DWK added int natchrts to functions
20060114 - DWK addedin include tmaxcohortdat437.h and standard
           includes

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

#include "tmaxcohortdat437.h" 

/* *************************************************************
************************************************************* */

MaxCohortdata43::MaxCohortdata43( void )
{

  chrtend = 1;
  lagpos = -99;
  curpos = 0;

};

/* *************************************************************
************************************************************* */

int MaxCohortdata43::get( ifstream& infile )
{

  lagpos = infile.tellg();

  infile >> col;
  infile >> row;
  infile >> varname;
  infile >> carea;
  infile >> year;
  infile >> total;
  infile >> natchrts;
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

int MaxCohortdata43::getdel( FILE* infile )
{
  char tmpvarname[40];
  char tmpcontnent[40];
  
  chrtend = fscanf( infile,"%f,%f, %s ,%ld,%d,%d,%d, %s",
                   &col, 
                   &row, 
                   tmpvarname, 
                   &carea,
                   &year, 
                   &total,
                   &natchrts, 
                   tmpcontnent );

  varname = tmpvarname;
  contnent = tmpcontnent;
  
  return chrtend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MaxCohortdata43::out( ofstream& ofile, 
                           const float& col, 
                           const float& row, 
                           const string& varname, 
                           const long& carea,
                           const int& year,  
                           const int& total,
                           const int& natchrts, 
                           const string& contnent )
{

  ofile.setf( ios::fixed,ios::floatfield );
  ofile.setf( ios::showpoint );
  ofile.precision( 1 );

  ofile << col << ' ';
  ofile << row << ' ';
  ofile << varname << ' ';
  ofile << setprecision( 0 ) << carea << ' ';
  ofile << setprecision( 0 ) << year << ' ';
  ofile << setprecision( 0 ) << total << ' ';
  ofile << setprecision( 0 ) << natchrts << ' ';
  ofile << contnent;
  ofile << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void MaxCohortdata43::outdel( ofstream& ofile, 
                              const float& col, 
                              const float& row, 
                              const string& varname, 
                              const long& carea,  
                              const int& year,
                              const int& total,
                              const int& natchrts, 
                              const string& contnent )
{

  ofile.setf( ios::fixed,ios::floatfield );
  ofile.setf( ios::showpoint );
  ofile.precision( 1 );

  ofile << col << ",";
  ofile << row << ", ";
  ofile << varname << " ,";
  ofile << setprecision( 0 ) << carea << ",";
  ofile << setprecision( 0 ) << year << ",";
  ofile << setprecision( 0 ) << total << ",";
  ofile << setprecision( 0 ) << natchrts << ", ";
  ofile << contnent;
  ofile << endl;

};

