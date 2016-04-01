/* *************************************************************
TIGSMCLMDAT44C.H - object to read and write the structure of the
                   climate data from files used by the 
                   Terrestrial Ecosystem Model (TEM)

20060114 - DWK created by modifying tclmdat425.h
20060114 - DWK changed class Clmdata to class Clmdata43
20060114 - DWK changed long year to int year in
           out(), outdel(), pctout() and pctoutdel()
20060114 - DWK changed char varname[9] to string varname in
           out(), outdel(), pctout() and pctoutdel()
20060114 - DWK changed char contnent[9] to string contnent in
           out(), outdel(), pctout() and pctoutdel()
20060114 - DWK changed public variable char varname[9] to
           string varname
20060114 - DWK changed public variable char contnent[9] to
           string contnent
20060114 - DWK deleted tclmdat425.cpp at bottom of file
20080130 - DWK changed include from temconsts43.hpp to
           temigsmconsts44a.hpp
20080130 - DWK changed class Clmdata43 to class Clmdata44
20110707 - DWK changed include from temigsmconsts44a.hpp to
           temigsmconsts44c.hpp
20150429 - DWK changed include temigsmconsts44c.hpp to
           temigsmconstants.hpp
                                
************************************************************* */

#ifndef TIGSMCLMDAT44D_H
#define TIGSMCLMDAT44D_H

#include "temigsmconstants.hpp"

class Clmdata44 
{
  
  public:
    
     Clmdata44( void );

/* *************************************************************
		      Public Functions
************************************************************* */

     // read data structure.
     
     int get( ifstream& ifile );
     int getdel( FILE* ifile );
     
     //write data structure.
     
     void out( ofstream& ofile, 
               const float& col, 
               const float& row, 
               const string& varname, 
               const int& carea, 
               const int& year, 
               double mon[CYCLE], 
               const string& contnent );
               
     void outdel( ofstream& ofile, 
                  const float& col, 
                  const float& row, 
                  const string& varname, 
                  const int& carea, 
                  const int& year, 
                  double mon[CYCLE], 
                  const string& contnent );
                  
     void pctout( ofstream& ofile, 
                  const float& col, 
                  const float& row, 
                  const string& varname, 
                  const int& carea, 
                  const int& year, 
                  double mon[CYCLE], 
                  const string& contnent );
                  
     void poutdel( ofstream& ofile, 
                   const float& col, 
                   const float& row, 
                   const string& varname, 
                   const int& carea, 
                   const int& year, 
                   double mon[CYCLE], 
                   const string& contnent );


/* *************************************************************
		     Public Variables
************************************************************* */

     // column or longitude of grid cell (degrees)	  
     float col;          

     // row or latitude of grid cell (degrees)
     float row;          
 
     // climate variable name
     string varname;     

     // area covered by grid cell (sq. km)
     int carea;          

      // date (year) of data
     long year;         

     // annual sum of monthly data for grid cell
     double total;       

      // maximum monthly value for grid cell
     double max;        

     // mean annual value for grid cell
     double ave;         

     // minimum monthly value for grid cell
     double min;         

     // monthly values for the grid cell
     double mon[CYCLE];  

      // name of continent containing grid cell
     string contnent;   


  private:

/* *************************************************************
		      Private Variables
************************************************************* */

     int clmend;
     long curpos;
     long lagpos;

};

#endif
