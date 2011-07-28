/* *************************************************************
****************************************************************
TIGSMLCLUC44A.H - determines both potential and "actual" land
                  cover, and identifies land uses within a grid
                  cell

20060422 - DWK created by modifying tlcluc437.h
20060422 - DWK changed include from lulcdat437.h to lulcdat44.h
20060422 - DWK changed public Lulcdata43 lulc to Lulcdata44 lulc
20060422 - DWK deleted public int backcastflag int backcastyears
           and FILE* frandom
20070426 - DWK changed lulcdat44.h to lulcdat44a.h
20080130 - DWK changed include from tbiome437.h to 
           tigsmbiome44a.h
20080130 - DWK changed Biome43 to Biome44
20080131 - DWK added public function initPotvegCohorts()
20080131 - DWK added public function initPotvegMaxCohorts()
20080131 - DWK added public string ipotlulcfname
20080131 - DWK added public string ipotmxcohrtfname
            
****************************************************************
************************************************************* */

#ifndef TIGSMLCLUC44A_H
#define TIGSMLCLUC44A_H

const int MXFRI = 2000;

// TEMlcluc43 uses the MaxCohortdata43 class
#include "tmaxcohortdat437.h" 

// TEMlcluc43 uses the LULCdata44 class
#include "lulcdat44a.h"    

// TEMlcluc43 inherits the Biome43 class
#include "tigsmbiome44a.h"    

class TEMlcluc44 : public Biome44
{

  public:

     TEMlcluc44();

 /* ************************************************************
		 Public Functions
************************************************************* */

     int getCohort( FILE* flulc );

     int getNumberOfCohorts( FILE* fnchrts );

     void initCohorts( ofstream& rflog1 );

     void initMaxCohorts( ofstream& rflog1 );

     void initPotvegCohorts( ofstream& rflog1 );

     void initPotvegMaxCohorts( ofstream& rflog1 );

     void setLCLUCFlags( ofstream& rflog1, const int& requil );


/* *************************************************************
		 Public Variables
************************************************************* */

     int agcmnt;

     MaxCohortdata43 cohorts;

     int currentveg;
     
     int endyr;
                 
     FILE* ifnumcohrts;
     
     FILE* iflulc;

     string ilulcend;
     string ilulcfname;

     string imxcohrtend;
     string imxcohrtfname;

     string ipotlulcfname;

     string ipotmxcohrtfname;
     
     int lastyr;
     
     Lulcdata44 lulc;

     int maxtype;

     int potveg;
 
     int startyr;
     
     // portion of carea in a particular cohort, km^2
     long subarea;

     int tlulcflag;

};

#endif
