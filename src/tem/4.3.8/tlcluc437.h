/* *************************************************************
****************************************************************
TLCLUC437.H - determines both potential and "actual" land
               cover, and identifies land uses within a grid
               cell

20060116 - DWK created by modifying tlcluc431b.h
20060116 - DWK changed include from tbiome431a.h to tbiome437.h
20060116 - DWK changed char ilulcfname[MAXFNAME] and
           char ilulcend[25] to string ilulcfname and
           string ilulcend, respectively
20060116 - DWK deleted include tlcluc431b.cpp at bottom of file
20060116 - DWK added global const int MXFRF = 2000;
20060116 - DWK added include tmaxcohortdat437.h
20060116 - DWK changed include from tvegdat425.h to tvegdat437.h
20060116 - DWK changed include from lulcdat425.h to lulcdat437.h
20060116 - DWK added public functions getNumberOfCohorts(), 
           getVegtype() initCohorts() and initMaxCohorts() 
20060116 - DWK added public int backcastflag, nt backcastyears, 
           Cohortdata43 cohorts, int endyr, FILE* frandom, 
           int FRF, FILE* ifcohrts, FILE* ifvtype, 
           string ivegfname, string ivegend, double prod10par,
           double prod100par, double sconvert, double slashpar,
           int startyr, long subarea, double vconvert, 
           Vegdata43 vegtype and double vrespar
20060116 - DWK deleted FILE* ifpotveg and Vegdata potveg  
20060116 - DWK deleted public function getPotentialVeg() and 
           initPotentialVeg()
20060116 - DWK renamed getLandUse() as getCohort()
20060116 - DWK added public string imxcohrtend and 
           string imxcohrtfname
20060116 - DWK deleted public Vegdata43 vegtype 
                                                                    
****************************************************************
************************************************************* */

#ifndef TLCLUC437_H
#define TLCLUC437_H

const int MXFRF = 2000;

// TEMlcluc43 uses the MaxCohortdata43 class
#include "tmaxcohortdat437.h" 

// TEMlcluc43 uses the LULCdata43 class
#include "lulcdat437.h"    

// TEMlcluc43 inherits the Biome43 class
#include "tbiome437.h"    

class TEMlcluc43 : public Biome43
{

  public:

     TEMlcluc43();

 /* ************************************************************
		 Public Functions
************************************************************* */

     int getCohort( FILE* flulc );
     int getNumberOfCohorts( FILE* fnchrts );

     void initCohorts( ofstream& rflog1 );
     void initMaxCohorts( ofstream& rflog1 );

     void setLCLUCFlags( ofstream& rflog1, const int& requil );


/* *************************************************************
		 Public Variables
************************************************************* */

     int agcmnt;

     // Flag to indicate if the transient simulation will be 
     //   preceeded with recurring disturbance
     
     int backcastflag;

     // number of years needed for backcasting a cohort
          
     int backcastyears;

     MaxCohortdata43 cohorts;

     int currentveg;
     
     int endyr;

     FILE* frandom;
              
     int FRF;  // Fire Return Frequency
     
     FILE* ifnumcohrts;
     FILE* iflulc;

     FILE* ifvtype;

     string ilulcend;
     string ilulcfname;

     string imxcohrtend;
     string imxcohrtfname;
     
     string ivegfname;
     string ivegend;

     int lastyr;
     
     Lulcdata43 lulc;

     int maxtype;

//     Vegdata potveg;
     int potveg;
 
     double prod10par;
     double prod100par;
     
     double sconvert;
     double slashpar;

     int startyr;
     
     // portion of carea in a particular cohort, km^2
     long subarea;

     int tlulcflag;

     double vconvert;

     double vrespar;

};

#endif
