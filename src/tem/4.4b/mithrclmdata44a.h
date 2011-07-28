/* *************************************************************
MITHRCLMDATA44A.H - object to read and write hourly temperature,
                     precipitation and solar radiation files 
                     from the MIT-IGSM climate model

Modifications:

20060128 - DWK created by modifying mithrclmdata436.h
20060128 - DWK added include mittemconsts43.hpp
20060128 - DWK changed public char varname[MAXVARNAME] to
           string varname
20060128 - DWK deleted include mithrclmdata436.cpp at bottom of 
           file
20080130 - DWK changed include from mittemconsts43.hpp to
           mittemconsts44a.hpp
20080130 - DWK changed class MITHrClmData43 to class
           MITHrClmData44
                                             
************************************************************* */

#ifndef MITHRCLMDATA44A_H
#define MITHRCLMDATA44A_H

#include "mittemconsts44a.hpp"


class MITHrClmData44
{

  public:

    MITHrClmData44( void );

/* *************************************************************
                Public Functions
************************************************************* */

    void get( ifstream& file );
    
    void writeclm( ofstream& ofile, 
                   string varname, 
                   const int& dm, 
                   const MITHrClmData44& clm );

/* *************************************************************
                Public variables
************************************************************* */

    int day;
 
    int hour;
 
    // data in latitude bands from the 2-D MIT L-O climate model
    double latband[MXMITNLAT][MAXCHRTS]; 

    int mon; // month
    
    string varname;
    

    int year; // year of transient data

};

#endif

