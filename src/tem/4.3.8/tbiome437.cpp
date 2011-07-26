/* *************************************************************
****************************************************************
TBIOME437.CPP - object describing general characteristics of
                   vegetation mosaic used in the Terrestrial
	           Ecosystem Model (TEM)

Modifications:

20060126 - DWK created by modifying tbiome50b5.cpp
20060126 - DWK changed include tbiome50b5.h to tbiome437.h 
20060126 - DWK changed Biome50:: to Biome43::
20060126 - DWK changed inheritance from ProcessXML50 to 
           ProcessXML43 in Biome43()
           
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


#include "tbiome437.h"


/* *************************************************************
************************************************************* */

Biome43::Biome43( void ) : ProcessXML43()
{
  temveg = -99;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Biome43::getCommunityType( const int& tveg )
{
  int mez;
  int communtype;
  
  mez = tveg - 1;
  if( mez < 0 || mez >= NUMVEG )
  {
    communtype = 1;
  }
  else { communtype = subtype[mez][0]; }

  return communtype;	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Biome43::getVegMosaic( const int& tveg )
{
  int mez;
  int maxtype;

  mez = tveg - 1;
  if( mez < 0 || mez >= NUMVEG )
  {
    maxtype = 1;
  }
  else { maxtype = numtype[mez]; }

  return maxtype;
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Biome43::getVegSubarea( const int& tveg,
                               const int& dtype,
                               const int& carea )
{
  int mez;
  double sarea;

  mez = tveg - 1;
  if( mez < 0 || mez >= NUMVEG )
  {
    sarea = (double) carea;
  }
  else { sarea = (double) carea * pcttype[mez][dtype] * 0.01; }

  return sarea;
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Biome43::getVegSubtype( const int& tveg, const int& dtype )
{
  int mez;
  int vegtype;

  mez = tveg - 1;
  if( mez < 0 || mez >= NUMVEG )
  {
    vegtype = 1;
  }
  else { vegtype = subtype[mez][dtype]; }

  return vegtype;
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Biome43::getvtype( ofstream& rflog1 )
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

void Biome43::getvtype( const string& ecd )
{
  ifstream infile;
  int dv;
  int dtype;
  int vegtype;
  int ez;

  ostringstream tempString;

  infile.open( ecd.c_str(), ios::in );

  if( !infile )
  {
    cerr << endl << "Cannot open " << ecd;
    cerr << " for community ECD input" << endl;
    exit( -1 );
  }

  getXMLrootNode( infile, "communityECD" );


  for( dv = 0; dv < NUMVEG; ++dv )
  {

    vegtype = getXMLtemvegNode( infile, "communityECD" );

    if( vegtype > NUMVEG )
    {
      cerr << endl << "TEMVEG type " << vegtype << endl;
      cerr << " cannot be greater than " << NUMVEG;
      cerr << " in communityECD" << endl;
      exit( -1 );
    }

    ez = vegtype - 1;

    numtype[ez] = getXMLtvegArrayInt( infile,
                                      "communityECD",
                                      "numtype",
                                      vegtype );

    for( dtype = 0; dtype < NUMMSAC; ++dtype )
    {
      tempString.str( "" );
      tempString << "subtype" << (dtype+1);
      subtype[ez][dtype] = getXMLtvegArrayInt( infile,
                                               "communityECD",
                                               tempString.str(),
                                               vegtype );

      tempString.str( "" );
      tempString << "pcttype" << (dtype+1);
      pcttype[ez][dtype] = getXMLtvegArrayDouble( infile,
                                                  "communityECD",
                                                  tempString.str(),
                                                  vegtype );
    }

    endXMLtvegNode( infile );
  }


  if( dv < NUMVEG )
  {
    cerr << endl << " Parameters found for only " << dv;
    cerr << " community types out of a maximum of ";
    cerr << (NUMVEG-1) << " types in communityECD" << endl;
    exit( -1 );
  }

  infile.close();

};


