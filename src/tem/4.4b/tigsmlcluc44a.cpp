/* **************************************************************
*****************************************************************
TIGSMLCLUC44A.CPP - determines both potential and "actual" land
                    cover, and identifies land uses within a grid
                    cell

20060422 - DWK created by modifying tlcluc44.cpp
20060422 - DWK changed include tlcluc437.h to tlcluc44.h
20060422 - DWK changed TEMlcluc43:: to TEMlcluc44::
20070426 - DWK changed include tlcluc44.h to tlcluc44a.h
20080130 - DWK changed include from tlcluc44a.h to tigsmlcluc44a.h
20080130 - DWK changed Biome43() to Biome44()
20080131 - DWK added initPotvegCohorts()
20080131 - DWK added initPotvegMaxCohorts()
           
*****************************************************************
************************************************************** */

#include<cstdio>

  using std::fscanf;
  using std::FILE;

#include<iostream>

  using std::cin;
  using std::cout;
  using std::ios;
  using std::endl;

#include<fstream>

  using std::ifstream;
  using std::ofstream;

#include<iomanip>

  using std::setprecision;

#include<string>

  using std::string;


#include "tigsmlcluc44a.h"

/* *************************************************************
************************************************************* */

TEMlcluc44::TEMlcluc44() : Biome44()
{

  cohorts.col = -999.9;
  cohorts.row = -999.9;
  cohorts.total = -99;
  cohorts.natchrts = -99;
  
  lulc.year = -999;
  lulc.agstate = -99;
  lulc.agprevstate = -99;
  lulc.standage = -99;

  maxtype = -999;

  agcmnt = -999;
  cmnt = -999;

  lastyr = -999;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int TEMlcluc44::getCohort( FILE* flulc )
{

  int gisend = lulc.getdel( flulc );
  
  if( -1 == gisend ) { return gisend; }

  return gisend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int TEMlcluc44::getNumberOfCohorts( FILE* fnchrts )
{

  int gisend = cohorts.getdel( fnchrts );
  
  if( -1 == gisend ) { return gisend; }

  return gisend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TEMlcluc44::initCohorts( ofstream& rflog1 )
{

  cout << "Please enter the first part of the file name ";
  cout << "containing the cohort data: " << endl;
  cout << "               (e.g., COHORT) " << endl;

  cin >> ilulcfname;

  cout << "Please enter the file extension (include the '.'): ";

  cin >> ilulcend;

  rflog1 << "Please enter the first part of the file name ";
  rflog1 << "containing the cohort data: " << endl;
  rflog1 << "               (e.g., COHORT) " << endl;
  rflog1 << ilulcfname << endl << endl;
  rflog1 << "Please enter the file extension (include the '.'): ";
  rflog1 << endl;
  rflog1 << ilulcend << endl << endl;

  if( 1 == tlulcflag )
  {
    cout << "Please enter the last year for which you have ";
    cout << "cohort data: " << endl;

    cin >> lastyr;

    rflog1 << "Please enter the last year for which you have ";
    rflog1 << "cohort data: " << endl;
    rflog1 << lastyr << endl << endl;
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TEMlcluc44::initMaxCohorts( ofstream& rflog1 )
{
 
  cout << "Please enter the first part of the file name ";
  cout << "containing the maximun number of cohort data";
  cout << endl;
  cout << "        (e.g., MXCOHRTS) " << endl;

  cin >> imxcohrtfname;
  
  cout << "Please enter the file extension (include the '.'): ";

  cin >> imxcohrtend;

  rflog1 << "Please enter the first part of the file name ";
  rflog1 << "containing the maximun number of cohort data";
  rflog1 << endl;
  rflog1 << "        (e.g., MXCOHRTS) " << endl;
  rflog1 << imxcohrtfname << endl << endl;
  rflog1 << "Please enter the file extension (include the '.'): ";
  rflog1 << imxcohrtend << endl << endl;

  return;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TEMlcluc44::initPotvegCohorts( ofstream& rflog1 )
{

  cout << "Please enter the file name containing the cohort ";
  cout << "data for potential vegetation: " << endl;
  cout << "               (e.g., POTCOHORT) " << endl;

  cin >> ipotlulcfname;

  rflog1 << "Please enter the file name containing the cohort ";
  rflog1 << "data for potential vegetation: " << endl;
  rflog1 << "               (e.g., POTCOHORT) " << endl;
  rflog1 << ipotlulcfname << endl << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TEMlcluc44::initPotvegMaxCohorts( ofstream& rflog1 )
{
 
  cout << "Please enter the file name containing the maximum ";
  cout << "number of cohorts data for potential vegetation:";
  cout << endl;
  cout << "        (e.g., POTMXCOHRTS) " << endl;

  cin >> ipotmxcohrtfname;
  
  rflog1 << "Please enter the file name containing the maximum ";
  rflog1 << "number of cohorts data for potential vegetation:";
  rflog1 << endl;
  rflog1 << "        (e.g., POTMXCOHRTS) " << endl;
  rflog1 << ipotmxcohrtfname << endl << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TEMlcluc44::setLCLUCFlags( ofstream& rflog1, 
                                const int& requil )
{

  tlulcflag = 0;

  if( 0 == requil )
  {
    cout << "Do you have transient land use data?:" << endl;
    cout << "Enter 0 for No:" << endl;
    cout << "Enter 1 for Yes: ";

    cin >> tlulcflag;

    rflog1 << "Do you have transient land use data?:" << endl;
    rflog1 << "Enter 0 for No:" << endl;
    rflog1 << "Enter 1 for Yes: ";
    rflog1 << tlulcflag << endl << endl;
  }

};

