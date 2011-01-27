/* **************************************************************
*****************************************************************
TLCLUC437.CPP - determines both potential and "actual" land
                 cover, and identifies land uses within a grid
                 cell

20060116 - DWK created by modifying tlcluc431b.cpp
20060116 - DWK added getNumberOfCohorts(), getVegtype(),  
           initCohorts(), initMaxCohorts()
20060116 - DWK added cohorts.col, cohorts.row, cohorts.total, 
           and vegtype.temveg to TEMlcluc43()
20060116 - DWK deleted potveg.col, potveg.row and potveg.temveg 
           from TEMlcluc43()
20060116 - DWK deleted getPotentialVeg() and  initPotentialVeg()
20060116 - DWK renamed getLandUse() as getCohort()
20060116 - DWK added include tlcluc437.h and standard includes
           
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


#include "tlcluc437.h"

/* *************************************************************
************************************************************* */

TEMlcluc43::TEMlcluc43() : Biome43()
{

  cohorts.col = -999.9;
  cohorts.row = -999.9;
  cohorts.total = -99;
  cohorts.natchrts = -99;
  
  lulc.year = -999;
  lulc.agstate = -99;
  lulc.agprevstate = -99;
  lulc.RAP = -999.9;

  maxtype = -999;

  agcmnt = -999;
  cmnt = -999;

  lastyr = -999;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int TEMlcluc43::getCohort( FILE* flulc )
{

  int gisend = lulc.getdel( flulc );
  
  if( -1 == gisend ) { return gisend; }

  return gisend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int TEMlcluc43::getNumberOfCohorts( FILE* fnchrts )
{

  int gisend = cohorts.getdel( fnchrts );
  
  if( -1 == gisend ) { return gisend; }

  return gisend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TEMlcluc43::initCohorts( ofstream& rflog1 )
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

void TEMlcluc43::initMaxCohorts( ofstream& rflog1 )
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

void TEMlcluc43::setLCLUCFlags( ofstream& rflog1, 
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

