/* *************************************************************
MITCLMDATA437.H - object to read and write the temperature,
                  precipitation and cloudiness files from the  
                  2-D MIT L-O climate model

Modifications:

20060128 - DWK created by modifying mitdata436.h
20060128 - DWK added include mittemconsts43.hpp
20060128 - DWK changed public char varname[MAXVARNAME] to
           string varname
20060128 - DWK deleted include mitclmdata436.cpp at bottom 
           of file
                       
************************************************************* */

#ifndef MITCLMDATA437_H
#define MITCLMDATA437_H

#include "mittemconsts43.hpp"

class MITCLMdata43 
{

  public:

    MITCLMdata43( void );

/* *************************************************************
                Public Functions
************************************************************* */

    void get( ifstream& file );
    
    void writeclm( ofstream& ofile, 
                   string varname, 
                   const int& dm, 
                   const MITCLMdata43& clm);

/* *************************************************************
                Public variables
************************************************************* */

    // data in latitude bands from the 2-D MIT L-O climate model
    double latband[MXMITNLAT][MAXCHRTS]; 
    
    int mon; // month
    
    string varname;

    int year; // year of transient data
    
};

#endif

