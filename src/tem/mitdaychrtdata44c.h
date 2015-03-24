/* *************************************************************
MITDAYCHRTDATA44C.H - object to read and write daily 
                      temperature, precipitation and cloudiness 
                      files from the 2-D MIT L-O climate model

Modifications:

20060128 - DWK created by modifying mitdaychrtdata437.h
20060128 - DWK added include mittemconsts43.hpp
20060128 - DWK changed public char varname[MAXVARNAME] to
           string varname
20060128 - DWK deleted include mitdaychrtdata436.cpp at bottom 
           of file
20080130 - DWK changes include from mittemconsts43.hpp to
           mittemconsts44a.hpp
20080130 - DWK changes class MITDayChrtData43 to class
           MITDayChrtData44
20110707 - DWK changed included from mittemconsts44a.hpp to
           mittemconsts44c.hpp                                            
************************************************************* */

#ifndef MITDAYCHRTDATA44C_H
#define MITDAYCHRTDATA44C_H

#include "mittemconsts44c.hpp"

class MITDayChrtData44 
{

  public:

    MITDayChrtData44( void );

/* *************************************************************
                Public Functions
************************************************************* */

    void get( ifstream& file );
    
    void writeclm( ofstream& ofile, 
                   string varname, 
                   const int& dm, 
                   const MITDayChrtData44& clm );

/* *************************************************************
                Public variables
************************************************************* */

    int day;

    // data in latitude bands from the 2-D MIT L-O climate model
    double latband[MXMITNLAT][CLMMXNVEG];

    int mon; // month

    string varname;

    int year; // year of transient data
    
};

#endif

