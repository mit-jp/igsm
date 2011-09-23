/* *************************************************************
MITDAYCLMDATA44C.H - object to read and write daily temperature,
                     precipitation and cloudiness files from the  
                     2-D MIT L-O climate model

Modifications:

20060128 - DWK created by modifying mitclmdata436.h
20060128 - DWK added include mittemconsts43.hpp
20060128 - DWK changed public char varname[MAXVARNAME] to
           string varname
20060128 - DWK deleted include mitclmdata436.cpp at bottom of 
           file
20080130 - DWK changed include from mittemconsts43.hpp to
           mittemconsts44a.hpp
20080130 - DWK changed class MITDayClmData43 to class
           MITDayClmData44
20110707 - DWK changed include from mittemconsts44a.hpp to
           mittemconsts44c.hpp
                                             
************************************************************* */

#ifndef MITDAYCLMDATA44A_H
#define MITDAYCLMDATA44A_H

#include "mittemconsts44c.hpp"

class MITDayClmData44
{

  public:

    MITDayClmData44( void );

/* *************************************************************
                Public Functions
************************************************************* */

    void get( ifstream& file );
    
    void writeclm( ofstream& ofile, 
                   string varname, 
                   const int& dm, 
                   const MITDayClmData44& clm );

/* *************************************************************
                Public variables
************************************************************* */

    int day;

    // data in latitude bands from the 2-D MIT L-O climate model
    double latband[MXMITNLAT]; 

    int mon; // month

    string varname;

    int year; // year of transient data
    
};

#endif

