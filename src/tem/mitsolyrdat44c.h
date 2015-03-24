/* *************************************************************
MITSOLYRDAT44C.H - object to read and write the structure of 
                   soil characteristrics data from/to files used 
                   by the Terrestrial Ecosystem Model (TEM)

Modifications:

20060302 - DWK created by modifying mitsoldat437.h
20080130 - DWK changed include from netconsts.hpp to
           nemconsts44a.hpp
20080130 - DWK changed class MITSoilLayerdata43 to 
           class MITSoilLayerdata44
20110707 - DWK changed include from nemconsts44a.hpp to 
           nemconsts44c.hpp
                                            
************************************************************* */

#ifndef MITSOLYRDAT44C_H
#define MITSOLYRDAT44C_H

#include "nemconsts44c.hpp"

class MITSoilLayerdata44 
{
  
  public:
	 
     MITSoilLayerdata44( void );

/* *************************************************************
		    Public Functions
************************************************************* */
     
     int get( ifstream& infile );
     
     int getdel( FILE* infile );
     
     void out( ofstream& ofile, 
               const float& col, 
               const float& row, 
               const string& varname,
               const int& layer,
               const double& thickness, 
               const double& porosity,
               const double& Ksat, 
               const string& source, 
               const string& region );
     
     void outdel( ofstream& ofile, 
                  const float& col, 
                  const float& row, 
                  const string& varname, 
                  const int& layer,
                  const double& thickness, 
                  const double& porosity,
                  const double& Ksat, 
                  const string& source, 
                  const string& region );


/* *************************************************************
		     Public Variables
************************************************************* */   
          
     // column or longitude of grid cell (degrees)	 
     float col;           

     // Saturated hydraulic conductivity for soil layers
     double Ksat;
     
     // Index of soil layer starting from soil surface
     int layer;
     
     // Porosity of soil layers as a fraction of soil volume
     double porosity;    

     // Thickness of soil layer
     double thickness;
          
     // name of continent containing grid cell      
     string region;    

     // row or latitude of grid cell (degrees)
     float row;           

     // reference to data source
     string source;

     // "SOILAYER"
     string varname;


  private:

/* *************************************************************
		      Private Variables
************************************************************* */

     int soilend;
     long curpos;
     long lagpos;

};

#endif

