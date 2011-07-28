/* **************************************************************
LATDAT437.H - object to read and write the structure of latitude
                 and longitude data from/to files

Modifications:

20060113 - DWK created by modifying latdat425.h
20060113 - DWK changed class Latdata to class Latdata43
20060113 - DWK changed public char varname[9] to string varname
20060113 - DWK changed public char contnent9] to string contnent
20060113 - DWK deleted include latdat425.cpp from bottom of file

************************************************************** */

#ifndef LATDAT437_H
#define LATDAT437_H

class Latdata43
{

  public:

     Latdata43( void );

/* **************************************************************
                      Public Functions
************************************************************** */

// read data structure.
     int get( ifstream& infile );
     
     int getdel( FILE* infile );

//write data structure.

     void out( ofstream& ofile, 
               const float& col, 
               const float& row, 
               const string& varname, 
               const double& lat, 
               const double& lon, 
               const string& contnent );

     void outdel( ofstream& ofile, 
                  const float& col, 
                  const float& row, 
                  const string& varname, 
                  const double& lat, 
                  const double& lon, 
                  const string& contnent );


/* **************************************************************
                     Public Variables
************************************************************** */

     float col;         // column of grid cell
     float row;         // row or of grid cell
     string varname;    // "LATITUDE?"
     double lat;        // latitude of grid cell (degrees)
     double lon;        // longitude of grid cell (degrees)
     string contnent;   // name of continent containing grid cell


  private:

/* **************************************************************
                      Private Variables
************************************************************** */

     int latend;
     long curpos;
     long lagpos;

};

#endif

