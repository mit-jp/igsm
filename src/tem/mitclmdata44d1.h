/* *************************************************************
MITCLMDATA44D1.H - object to read and write the temperature,
                  precipitation and cloudiness files from the  
                  2-D MIT L-O climate model

Modifications:

20060128 - DWK created by modifying mitdata436.h
20060128 - DWK added include mittemconsts43.hpp
20060128 - DWK changed public char varname[MAXVARNAME] to
           string varname
20060128 - DWK deleted include mitclmdata436.cpp at bottom 
           of file
20080130 - DWK changed include from mittemconsts43.hpp to
           mittemconsts44a.hpp
20080130 - DWK changed class MITCLMdata43 to class
           MITCLMdata44
20110707 - DWK changed include from mittemconsts44a.hpp to
           mittemconsts44c.hpp
20150428 - DWK changed include from mittemconsts44c.hpp to
           temigsmconstants.hpp
                                             
************************************************************* */

#ifndef MITCLMDATA44D_H
#define MITCLMDATA44D_H

#include "temigsmconstants.hpp"

class MITCLMdata44 
{

  public:

    MITCLMdata44( void );

/* *************************************************************
                Public Functions
************************************************************* */

    void get( ifstream& file );
    
    void writeclm( ofstream& ofile, 
                   string varname, 
                   const int& dm, 
                   const MITCLMdata44& clm);

/* *************************************************************
                Public variables
************************************************************* */

    // data in latitude bands from the 2-D MIT L-O climate model
    double latband[MXMITNLAT][CLMMXNVEG]; 
    
    int mon; // month
    
    string varname;

    int year; // year of transient data
    
};

#endif

