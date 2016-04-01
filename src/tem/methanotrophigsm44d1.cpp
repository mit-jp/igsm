/* *************************************************************
METHANOTROPHIGSM44C.CPP - Describes activity of methanotrophs
                          (i.e. methane oxidizers)
                     
Created by Q. Zhuang, 11/March/2003

Modifications:

20050316 - DWK changed CH4DMOXI:: to Methanotroph51::
20050316 - DWK renamed MethaneOR() to consumeCH4()
20050316 - DWK changed char ecd[80] to string ecd in 
           getch4oxi() and getecdch4oxi()
20060105 - DWK changed getch4oxi() to be getecd()
20060105 - DWK deleted getecdch4oxi()
20060207 - DWK added double oxi_c[MAXCMNT] to getecd()
20060207 - DWK added inheritance of ProcessXML50() to 
           Methanotroph51()
20060207 - DWK added XML commands to getecd()
20080130 - DWK changed include from methanotroph51.h to
           methanotrophigsm44a.h
20080130 - DWK changed Methanotroph51:: to Methanotroph44::
20080130 - DWK changed ProcessXML43() to ProcessXML44()
20110707 - DWK changed include from methanotrophigsm44a.h to
           methanotrophigsm44c.h
20150428 - DWK changed include from methanotrophigsm44c.h to
           methanotrophigsm44d1.h
           
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
  
#include "methanotrophigsm44d1.h"

/* *************************************************************
************************************************************* */

Methanotroph44::Methanotroph44() : ProcessXML44()
{


};

/* *************************************************************
************************************************************* */


/***************************************************************
************************************************************* */

double Methanotroph44::consumeCH4( const int& pcmnt,
                                   const double& ch4con, 
                                   const double& oxyc, 
                                   const double& soilt,
                                   const double& soilm, 
                                   const double& ehl, 
                                   const int& icul )
{
// Methane oxidation occurred between the soil surface and the 
//   water table the soil column is devided as 1 cm depth 
//   interval (formerly known as MethaneOR())

  double fmc; 
  double foxy; 
  double fst;
  double fsm; 
  double frx; 
  double fct;
  
  // rate of methane oxidation rate at depth z and time t
  double methaneor; 

  // Effects of methane concentration on methanotrophic reaction
  
  fmc = EffectMC( pcmnt, ch4con );
  
  
  // Effects of oxygen concentration to oxidation
  
  foxy = EffectOXY( pcmnt, oxyc );
  
  
  // Effects of soil temperature on methanogenesis

  fst = EffectST( pcmnt, soilt );
  
  
  // Effects of soil moisture on methanogenesis

  fsm = EffectSM( pcmnt, soilm );
  
  
  // Effects of soil redox potential on methanotrophy

  frx = EffectRX( ehl );
  
  
  //effects of agricultural intensification on methanotrophy

  fct = EffectCT( icul );

  methaneor = omax[pcmnt] * fmc * foxy * fst * fsm * frx * fct;

  return methaneor;
  
};

/***************************************************************
************************************************************* */


/***************************************************************
************************************************************* */

double Methanotroph44::EffectCT( const int& icul )
{

//effects of agricultural intensification
//According to Ridgwell et al., 1999

  double fct;
  double temp;
  // ICUL is the fractional intensity of cultivation 
  //   (Matthews, 1983)
  
  switch( icul )
  {
    case 1:  temp = ZERO; break;
    case 2:  temp = 0.2; break;
    case 3:  temp = 0.5; break;
    case 4:  temp = 0.75; break;
    case 5:  temp = 1.0; break;
    default: temp = ZERO; 
  }
 
  fct = 1.0 - 0.75 * temp;

  return fct;
  
};

/***************************************************************
************************************************************* */


/***************************************************************
************************************************************* */

double Methanotroph44::EffectMC( const int& pcmnt,
                                 const double& ch4con )
{

// Effects of methane concentration on methanotrophic reaction

  double fmc;
  
  // ch4con is soil methane concentration at depth z and time t

 fmc = ch4con / (kc[pcmnt] + ch4con);

 return fmc;
 
};

/***************************************************************
************************************************************* */


/***************************************************************
************************************************************* */

double Methanotroph44::EffectOXY( const int& pcmnt,
                                  const double& oxyc )
{

// Effects of oxygen concentration to oxidation

  double foxy;
  double kk;
    
  // oxyc denotes the oxygen concentration in the soil
 
  kk = ko[pcmnt] * 32 / pow( 10.0, 3.0 );
  foxy = oxyc / (kk + oxyc);

 return foxy;
 
};

/***************************************************************
************************************************************* */


/***************************************************************
************************************************************* */

double Methanotroph44::EffectRX( const double& ehl )
{

// Effects of soil redox potential on methanotrophy
// Follow Zhang et al., 2001

  double frx;

  if( (ehl <= -100.0) && (ehl >= -200) )  
  {
    frx = 0.0075 * ehl + 1.5;
  }
  
  if( (ehl > -100.0) && (ehl <= 200.0) )  
  {
    frx = 0.00083 * ehl + 0.833;
  }
  
  if( (ehl > 200) && (ehl <= 600) ) { frx = 1.0; }

 return frx;
 
};

/***************************************************************
************************************************************* */


/***************************************************************
************************************************************* */

double Methanotroph44::EffectSM( const int& pcmnt,
                                 const double& soilm )
{

// Similar to TEM for the heterotrophic respiration Tian et al.,

  double fsm;
  double temp1;
  double temp2;
  
  temp1 = (soilm - mvmin[pcmnt]) * ( soilm - mvmax[pcmnt]);
  temp2 = pow( (soilm - mvopt[pcmnt]), 2.0 );
  fsm = temp1 / (temp1 - temp2);

  return fsm;

};

/***************************************************************
************************************************************* */


/***************************************************************
************************************************************* */

double Methanotroph44::EffectST( const int& pcmnt,
                                 const double& soilt )
{

// Effects of soil temperature on methanotrophy

  double fst;

// reference soil temperature for methanotrophy,
// double tref = 23.0; 

// TSOIL(z,t) is the soil temperature at depth z and time t.

// tref is the reference soil temperature (oC) (e.g. 23oC).

  fst = pow( och4q10[pcmnt], (soilt - oxiref[pcmnt]) / 10.0 );
  
  // to limit the oxidation 20 Oct 2003
  if( soilt < ZERO ) { fst = ZERO; } 

  return fst;

};

/***************************************************************
************************************************************* */



/***************************************************************
************************************************************* */

void Methanotroph44::getecd( ofstream& rflog1 )
{

  string ecd;

  cout << "Enter name of the file (.ECD) with methanotroph ";
  cout << "parameter values:" << endl;
  
  cin >> ecd;

  rflog1 << "Enter name of the file (.ECD) with methanotroph ";
  rflog1 << "parameter values:";
  rflog1 << ecd << endl << endl;

  getecd( ecd );

};

/***************************************************************
************************************************************* */


/***************************************************************
************************************************************* */

void Methanotroph44::getecd( const string& ecd )
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

  getXMLrootNode( infile, "methanotrophECD" );

  for( dcmnt = 1; dcmnt < MAXCMNT; ++dcmnt )
  {
    comtype = getXMLcommunityNode( infile, "methanotrophECD" );

    if( comtype >= MAXCMNT )
    {
      cerr << endl << "comtype is >= MAXCMNT" << endl;
      cerr << "comtype cannot be greater than ";
      cerr << (MAXCMNT-1) << endl;
      cerr << " in methanotrophECD" << endl;
      
      exit( -1 );
    }

    omax[comtype]= getXMLcmntArrayDouble( infile,
                                          "methanotrophECD",
                                          "omax",
                                          comtype );

    kc[comtype]= getXMLcmntArrayDouble( infile,
                                        "methanotrophECD",
                                        "kc",
                                        comtype );

    och4q10[comtype]= getXMLcmntArrayDouble( infile,
                                             "methanotrophECD",
                                             "och4q10",
                                             comtype );

    ko[comtype]= getXMLcmntArrayDouble( infile,
                                        "methanotrophECD",
                                        "ko",
                                        comtype );

    oxi_c[comtype]= getXMLcmntArrayDouble( infile,
                                           "methanotrophECD",
                                           "oxi_c",
                                           comtype );

    afp[comtype]= getXMLcmntArrayDouble( infile,
                                         "methanotrophECD",
                                         "afp",
                                         comtype );

    mvmax[comtype]= getXMLcmntArrayDouble( infile,
                                           "methanotrophECD",
                                           "mvmax",
                                           comtype );

    mvmin[comtype]= getXMLcmntArrayDouble( infile,
                                           "methanotrophECD",
                                           "mvmin",
                                           comtype );

    mvopt[comtype]= getXMLcmntArrayDouble( infile,
                                           "methanotrophECD",
                                           "mvopt",
                                           comtype );

    oxiref[comtype]= getXMLcmntArrayDouble( infile,
                                            "methanotrophECD",
                                            "oxiref",
                                            comtype );

    endXMLcommunityNode( infile );
  }

  if( dcmnt < MAXCMNT )
  {
    cerr << endl << " Parameters found for only " << dcmnt;
    cerr << " community types out of a maximum of ";
    cerr << (MAXCMNT-1) << " types in methanotrophECD" << endl;
    
    exit( -1 );
  }
  
  infile.close();

};


