/* *************************************************************
****************************************************************
TIGSMBIOME44A.CPP - object describing general characteristics of
                    vegetation mosaic used in the Terrestrial
	                  Ecosystem Model (TEM)

Modifications:

20060126 - DWK created by modifying tbiome50b5.cpp
20060126 - DWK changed include tbiome50b5.h to tbiome437.h 
20060126 - DWK changed Biome50:: to Biome43::
20060126 - DWK changed inheritance from ProcessXML50 to 
           ProcessXML43 in Biome43()
20080130 - DWK changed include from tbiome437.h to 
           tigsmbiome44a.h
20080130 - DWK changed Biome43:: to Biome44::
20080130 - DWK changed ProcessXML43 to ProcessXML44
20080130 - DWK first introduced TEM constants FIRSTVEG,
           MISSCOMMUN and MISSVEG to code
                      
****************************************************************
************************************************************* */

#include<iostream>

  using std::cin;
  using std::cout;
  using std::ios;
  using std::cerr;
  using std::endl;

#include<fstream>

  using std::ifstream;
  using std::ofstream;

#include<cstdlib>

  using std::exit;
  using std::atof;
  using std::atoi;
  
#include<string>
  
  using std::string;

#include<sstream>

  using std::ostringstream;


#include "tigsmbiome44a.h"


/* *************************************************************
************************************************************* */

Biome44::Biome44( void ) : ProcessXML44()
{
  temveg = -99;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Biome44::getCommunityType( const int& tveg )
{
  int communtype;
  
  
  if( tveg < FIRSTVEG || tveg >= NUMVEG )
  {
    communtype = MISSCOMMUN;
  }
  else { communtype = subtype[tveg][0]; }

  return communtype;	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Biome44::getVegMosaic( const int& tveg )
{
  int maxtype;
  
  if( tveg < FIRSTVEG || tveg >= NUMVEG )
  {
    maxtype = 1;
  }
  else { maxtype = numtype[tveg]; }

  return maxtype;
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Biome44::getVegSubarea( const int& tveg,
                               const int& dtype,
                               const int& carea )
{
  double sarea;

  if( tveg < FIRSTVEG || tveg >= NUMVEG )
  {
    sarea = (double) carea;
  }
  else { sarea = (double) carea * pcttype[tveg][dtype] * 0.01; }

  return sarea;
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Biome44::getVegSubtype( const int& tveg, const int& dtype )
{
  int vegtype;

  if( tveg < FIRSTVEG || tveg >= NUMVEG )
  {
    vegtype = MISSVEG;
  }
  else { vegtype = subtype[tveg][dtype]; }

  return vegtype;
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Biome44::getvtype( ofstream& rflog1 )
{

  string ecd;

  cout << endl;
  cout << "Enter name of the file prescribing vegetation mosaics (.ECD):";
  cout << endl;
  
  cin >> ecd;

  rflog1 << endl;
  rflog1 << "Enter name of the file prescribing vegetation mosaics (.ECD):";
  rflog1 << endl << ecd << endl << endl;

  getvtype( ecd );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Biome44::getvtype( const string& ecd )
{
  int dv;

  int dtype;

  int ez;

  ifstream infile;

  ostringstream tempString;

  infile.open( ecd.c_str(), ios::in );

  if( !infile )
  {
    cerr << endl << "Cannot open " << ecd;
    cerr << " for community ECD input" << endl;
    exit( -1 );
  }

  getXMLrootNode( infile, "communityECD" );


  for( dv = FIRSTVEG; dv < NUMVEG; ++dv )
  {

    ez = getXMLtemvegNode( infile, "communityECD" );

    if( ez > NUMVEG )
    {
      cerr << endl << "TEMVEG type " << ez << endl;
      cerr << " cannot be greater than " << NUMVEG;
      cerr << " in communityECD" << endl;
      exit( -1 );
    }


    numtype[ez] = getXMLtvegArrayInt( infile,
                                      "communityECD",
                                      "numtype",
                                      ez );

    for( dtype = 0; dtype < NUMMSAC; ++dtype )
    {
      tempString.str( "" );
      tempString << "subtype" << (dtype+1);
      subtype[ez][dtype] = getXMLtvegArrayInt( infile,
                                               "communityECD",
                                               tempString.str(),
                                               ez );

      tempString.str( "" );
      tempString << "pcttype" << (dtype+1);
      pcttype[ez][dtype] = getXMLtvegArrayDouble( infile,
                                                  "communityECD",
                                                  tempString.str(),
                                                  ez );
    }

    endXMLtvegNode( infile );
  }


  if( dv < NUMVEG )
  {
    cerr << endl << " Parameters found for only " << (dv - FIRSTVEG +1);
    cerr << " community types out of a maximum of ";
    cerr << NUMVEG << " types in communityECD" << endl;

    exit( -1 );
  }

  infile.close();

};


