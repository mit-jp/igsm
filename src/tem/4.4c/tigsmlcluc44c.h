/* *************************************************************
****************************************************************
TIGSMLCLUC44C.H - determines both potential and "actual" land
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
20110705 - DWK added include tigsmveg44a.h
20110705 - DWK added include mitsoil44a.h
20110705 - DWK added public Tveg44 abandonedVeg and 
           Tveg44 convertedVeg
20110705 - DWK added public MITsoil44 abandonedSoil and
           MITsoil44 convertedSoil
20110707 - DWK changed include from tigsmveg44a.h to 
           tigsmveg44c.h
20110707 - DWK changed include mitsoil44a.h to mitsoil44c.h
20110707 - DWK changed include tigsmbiome44a.h to 
           tigsmbiome44c.h
                      
****************************************************************
************************************************************* */

#ifndef TIGSMLCLUC44C_H
#define TIGSMLCLUC44C_H

const int MXFRI = 2000;

// TEMlcluc44 uses the MaxCohortdata43 class
#include "tmaxcohortdat437.h" 

// TEMlcluc44 uses the LULCdata44 class
#include "lulcdat44a.h"    

// TEMlcluc44 inherits the Biome44 class
#include "tigsmbiome44c.h"    

// TEMlcluc44 uses the Tveg44 class
#include "tigsmveg44c.h"

// TEMlcluc44 uses the MITsoil44 class
#include "mitsoil44c.h"

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

     MITsoil44 abandonedSoil;

     Tveg44 abandonedVeg;
     
     int agcmnt;

     MaxCohortdata43 cohorts;

     MITsoil44 convertedSoil;

     Tveg44 convertedVeg;

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
