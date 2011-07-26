/* *************************************************************
O3DATA437.H - object to read and write the ozone files from the 
              MIT L-O climate model

20060128 - created by DWK by modifying o3data436.h 
20060128 - DWK added include mittemconsts43.hpp
20060128 - DWK changed public char varname[MAXVARNAME] to
           string varname
20060128 - DWK deleted include o3data436.cpp at bottom of file
           
************************************************************* */

#ifndef O3DATA437_H
#define O3DATA437_H

#include "mittemconsts43.hpp"

class O3data43 
{

  public:

    O3data43( void );

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

