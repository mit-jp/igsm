/* *************************************************************
O3DATA44C.H - object to read and write the ozone files from the 
              MIT L-O climate model

20060128 - created by DWK by modifying o3data436.h 
20060128 - DWK added include mittemconsts43.hpp
20060128 - DWK changed public char varname[MAXVARNAME] to
           string varname
20060128 - DWK deleted include o3data436.cpp at bottom of file
20080130 - DWK changed include from mittemconst43.hpp to
           mittemconsts44a.hpp
20080130 - DWK changed class O3data43 to class O3data44
20110707 - DWK changed include from mittemconsts44a.hpp to
           mittemconsts44c.hpp
20150429 - DWK changed include from mittemconsts44c.hpp to
           temigsmconstants.hpp
                      
************************************************************* */

#ifndef O3DATA44D_H
#define O3DATA44D_H

#include "temigsmconstants.hpp"

class O3data44 
{

  public:

    O3data44( void );

/* *************************************************************
                Public Functions
************************************************************* */

    void get( ifstream& file );

/* *************************************************************
                Public variables
************************************************************* */

    // data in latitude bands from the 2-D MIT L-O climate model
    double latbandhr[MXMITNLAT][MAX3HRS]; 

    int mon; // month
    
    string varname;

    int year; // year of transient data
    
};

#endif

