/* *************************************************************
METHANOGENIGSM44D1.CPP - Describes methanogen dynamics in 
                        terrestrial ecosystems

Created by Q. Zhuang, 19/Feb/2003

Modifications:

20060105 - DWK changed CH4DMPRO:: to be Methanogen51::
20060105 - DWK changed MethanePR() to produceCH4()  
20060105 - DWK changed char ecd[80] to be string ecd in
           getecdch4pro() and getch4pro()
20060105 - DWK changed getch4pro() to be getecd()
20060105 - DWK deleted getecdch4pro()
20060207 - DWK added inheritance of ProcessXML50() to 
           Methanogen51()
20060207 - DWK added XML commands to getecd()
20080130 - DWK changed include from methanogen51.h to
           methanogenigsm44a.h
20080130 - DWK changed Methanogen51:: to Methanogen44::
20080130 - DWK changed ProcessXML43() to ProcessXML44()
20110707 - DWK changed include from methanogenigsm44a.h to
           methanogenigsm44c.h
20150428 - DWK changed include from methanogenigsm44c.h to
           methanogenigsm44d1.h
           
************************************************************* */

#include<iostream>

  using std::cout;
  using std::cin;
  using std::ios;
  using std::cerr;
  using std::endl;

#include<fstream>

  using std::ifstream;
  using std::ofstream;

#include<cstdlib>

  using std::exit;

#include<cmath>

  using std::exp;
  using std::pow;

#include<string>
  using std::string;
  
#include "methanogenigsm44d1.h"

/* *************************************************************
************************************************************* */

Methanogen44::Methanogen44() : ProcessXML44()
{


};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Methanogen44::EffectOM( const int& pcmnt,
                               const double& vegNPP )
{
  // Effects of organic matter on methanogenesis.  Based on
  //   equation 2 in Walter and Heimann (2000)

  // vegNPP is either prescribed or simulated by TEM

  double fsom;
 
  // fresh organic matter from the vegetation
  double freshorg; 


  if( vegNPP <= ZERO ) { freshorg = 0.01 * maxfresh[pcmnt]; }
  else { freshorg = maxfresh[pcmnt] + vegNPP; }

  fsom = freshorg / maxfresh[pcmnt];

 return fsom;
 
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Methanogen44::EffectOMD( const int& pcmnt,
                                const double& depthz, 
                                const double& rootd )
{

  // Effects of organic matter distribution in the soil on 
  //   methanogenesis.  Based on equations 3 and 4 in
  // Walter and Heimann (2000)

  // the index of organic matter distribution as function of the 
  //   rooting depth and upper and lower boundary  
  double fcdis; 
  
  double rtz;

  // for vegetated area, vegetation type (2-34)
  // hydm_xtext() produce rooting depth, as meter, 
  // need to transfer to mm = 1000* rootz

  rtz = rootd * 1000.0;  // convert m to mm for rooting depth

  if( pcmnt != 0 )
  {
    if( (depthz < rtz) && (depthz > ZERO) ) { fcdis = 1.0; }
    if( (depthz > rtz) && (depthz < lowb[pcmnt]) )
    {
      fcdis = exp( -1.0 * (depthz - rtz) / 10.0 );
    }
  }
  else // unvegetated area
  {
    fcdis = 0.875 * exp( -1.0 * depthz / 10.0 );
  }


 return fcdis;
}

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Methanogen44::EffectPH( const double& soilph )
{
// Effects of soil PH to methanogenesis
// need reading in soil pH, similar to elevation

  // const double phmin = 5.5;
  const double phmin = 3.0;
  const double phmax = 9.0;
  const double phopt = 7.5;

  double fph;
  double v1;
  double v2; 
  double v3;

  v1 = soilph - phmin;
  v2 = soilph - phmax;
  v3 = soilph - phopt;

  fph = v1 * v2 / (v1 * (v2 - pow( v3, 2.0 )));

  return fph;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Methanogen44::EffectRX( const double& ehl )
{

// Effects of soil redox potential to methanogenesis

  double frxp;

  if( ehl <= -200.0 ) { frxp = 1.0; }
  else
  {
    if( (ehl > -200) && (ehl <= -100.0) ) 
    {
      frxp = -0.01 * ehl - 1;
    }
    
    if( ehl >= -100.0 ) { frxp = ZERO; }
  }

 return frxp;
 
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Methanogen44::EffectST( const int& pcmnt,
                               const double& soilt )
{

// Effects of soil temperature to methanogenesis

  double fmst;
  double tt;

//  double dynamicQ10;
//  double tempt1;
//  double tempt2;
//  double tempt;  // in unit K

  // reference soil temperature for methanogenesis, 7
//  double tref = 10.0; 

  // ch4q10 is a coefficient with a constant, 
  //   (0.6 as used in Walter et al., )

//  fmst = exp( soilt ) 
//         * pow( pmethaneq[pcmnt], (soilt - tref) / 10.0 );
  
  tt = (soilt - proref[pcmnt]) / 10.0;

  // The following 6 lines to implement the dynamics of Q10 
  //   approach, [Gedney and Cox]
  
  // Oke, T. R. Boundary Layer Climates
//  tempt1 = proref + 273.2;   
//  tempt1 = pow( tempt1, 2.0 );
//  tempt2 = soilt + 273.2;
//  tempt2 = pow( tempt2, 2.0 );
//  tempt = tempt1 / tempt2;
//  dynamicQ10 = pow( pmethaneq[pcmnt], tempt );

//  fmst = pow( dynamicQ10, tt );

  fmst = pow( pmethaneq[pcmnt], tt );

//  printf( " %3.2f %3.2f %3.2f %3.2f %3.2f\n", 
//          pmethaneq[pcmnt], 
//          tt, 
//          soilt, 
//          proref[pcmnt], 
//          fmst );

  return fmst;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Methanogen44::getecd( ofstream& rflog1 ) 
{
  string ecd;

  cout << "Enter name of file (.ECD) with the methanogen ";
  cout << "parameter values:" << endl;
  
  cin >> ecd;

  rflog1 << "Enter name of file (.ECD) with the methanogen ";
  rflog1 << "parameter values:" << endl;

  getecd( ecd );
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Methanogen44::getecd( const string& ecd )
{

  int comtype;
  int dcmnt;
  ifstream infile;


  infile.open( ecd.c_str(), ios::in );

  if( !infile )
  {
    cerr << endl << "Cannot open " << ecd;
    cerr << " for data input" << endl;
    exit( -1 );
  }

  getXMLrootNode( infile, "methanogenECD" );

  for( dcmnt = 1; dcmnt < MAXCMNT; ++dcmnt )
  {
    comtype = getXMLcommunityNode( infile, "methanogenECD" );

    if( comtype >= MAXCMNT )
    {
      cerr << endl << "comtype is >= MAXCMNT" << endl;
      cerr << "comtype cannot be greater than ";
      cerr << (MAXCMNT-1) << endl;
      cerr << " in methanogenECD" << endl;
      
      exit( -1 );
    }

    mgo[comtype]= getXMLcmntArrayDouble( infile,
                                         "methanogenECD",
                                         "mgo",
                                         comtype );

    kc[comtype]= getXMLcmntArrayDouble( infile,
                                        "methanogenECD",
                                        "kc",
                                        comtype );

    pmethaneq[comtype]= getXMLcmntArrayDouble( infile,
                                               "methanogenECD",
                                               "pmethaneq",
                                               comtype );

    maxfresh[comtype]= getXMLcmntArrayDouble( infile,
                                              "methanogenECD",
                                              "maxfresh",
                                              comtype );

    lowb[comtype]= getXMLcmntArrayDouble( infile,
                                          "methanogenECD",
                                          "lowb",
                                          comtype );

    proref[comtype]= getXMLcmntArrayDouble( infile,
                                            "methanogenECD",
                                            "proref",
                                            comtype );

    endXMLcommunityNode( infile );
  }

  if( dcmnt < MAXCMNT )
  {
    cerr << endl << " Parameters found for only " << dcmnt;
    cerr << " community types out of a maximum of ";
    cerr << (MAXCMNT-1) << " types in methanogenECD" << endl;
    
    exit( -1 );
  }

  infile.close();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Methanogen44::produceCH4( const int& pcmnt, 
                                 const double& vegNPP,  
                                 const double& depthz, 
                                 const double& rootd,   
                                 const double& soilt, 
                                 const double& soilph,  
                                 const double& ehl )
{

// Methane production rate occurred between the water table and 
//   the lower boundary the soil column is devided as 5 cm depth 
//   interval ( formerly MethanePR())

  // rate of methane production at depth z and time t
  double methanepr; 

  double fsom;
  double fcdis;
  double fmst;
  double fph;
  double frx;
  

  // Effects of the amount of organic matter on methanogenesis
  
  fsom = EffectOM( pcmnt, vegNPP );
  
  
  // Effects of the distribution of organic matter in the soil  
  //   on methanogenesis
  
  fcdis =  EffectOMD( pcmnt, depthz, rootd ); 
  
  
  // Effects of soil temperature to methanogenesis
  
  fmst = EffectST( pcmnt, soilt );
  
  
  // Effects of soil pH on methanogenesis
  
  fph = EffectPH( soilph );
  
  
  // Effects of soil redox potential to methanogenesis
  
  frx = EffectRX( ehl );
  
  methanepr = mgo[pcmnt] * fsom * fcdis * fmst * fph * frx;


  return methanepr;

};



