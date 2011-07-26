/* *************************************************************
MITSOLDAT437.H - object to read and write the structure of soil
                 characteristrics data from/to files used by the  
               Terrestrial Ecosystem Model (TEM)

Modifications:

20060129 - DWK created by modifying mitsoldat436.h
20060129 - DWK added include nemconsts.hpp
20060129 - DWK changed public char region[MAXREGNAME] to 
           string region
20060129 - DWK changed public char char source[9] to 
           string source
20060129 - DWK changed public char char char varname[MAXVARNAME] 
           to string varname
20060129 - DWK deleted include mitsoldat436.cpp from bottom of
           the file
20060306 - DWK public deleted porosity[] amd Ksat[]                     

************************************************************* */

#ifndef MITSOLDAT437_H
#define MITSOLDAT437_H

#include "nemconsts.hpp"

class MITSoildata43 
{
  
  public:
	 
     MITSoildata43( void );

/* *************************************************************
		    Public Functions
************************************************************* */
     
     int get( ifstream& infile );
     
     int getdel( FILE* infile );
     
     void out( ofstream& ofile, 
               const float& col, 
               const float& row, 
               const string& varname, 
               const long& carea, 
               const double& pctsand, 
               const double& pctsilt, 
               const double& pctclay, 
               const int& wsoil,
               const double& pH,
               const string& source, 
               const string& region );
     
     void outdel( ofstream& ofile, 
                  const float& col, 
                  const float& row, 
                  const string& varname, 
                  const long& carea, 
                  const double& pctsand, 
                  const double& pctsilt, 
                  const double& pctclay, 
                  const int& wsoil,
                  const double& pH, 
                  const string& source, 
                  const string& region );

/* *************************************************************
		     Public Variables
************************************************************* */   
          
     // area covered by grid cell (sq. km)     
     long carea;

     // column or longitude of grid cell (degrees)	 
     float col;           
    
     // percent clay of grid cell's soil texture
     double pctclay;      
     
     // percent sand of grid cell's soil texture           
     double pctsand;      
     
     // percent silt of grid cell's soil texture
     double pctsilt;      
     
     // Soil pH
     double pH;
    
     // name of continent containing grid cell      
     string region;    

     // row or latitude of grid cell (degrees)
     float row;           

     // reference to data source
     string source;

     // "TEXTURE"
     string varname;

     // wetland soil type designation (categorical data)
     int wsoil;           


  private:

/* *************************************************************
		      Private Variables
************************************************************* */

     int soilend;
     long curpos;
     long lagpos;

};

#endif

