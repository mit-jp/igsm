/* *************************************************************
MITHRCLMDATA437.H - object to read and write hourly temperature,
                     precipitation and solar radiation files 
                     from the MIT-IGSM climate model

Modifications:

20060128 - DWK created by modifying mithrclmdata436.h
20060128 - DWK added include mittemconsts43.hpp
20060128 - DWK changed public char varname[MAXVARNAME] to
           string varname
20060128 - DWK deleted include mithrclmdata436.cpp at bottom of 
           file
                       
************************************************************* */

#ifndef MITHRCLMDATA437_H
#define MITHRCLMDATA437_H

#include "mittemconsts43.hpp"


class MITHrClmData43 
{

  public:

    MITHrClmData43( void );

/* *************************************************************
                Public Functions
************************************************************* */

    void get( ifstream& file );
    
    void writeclm( ofstream& ofile, 
                   string varname, 
                   const int& dm, 
                   const MITHrClmData43& clm );

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

