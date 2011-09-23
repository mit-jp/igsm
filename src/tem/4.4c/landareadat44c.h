/* *************************************************************
LANDAREADAT44C.H - object to read and write the structure of 
                       grid cell data representing land area 
                       from/to files used by the 
                       Terrestrial Ecosystem Model

Modifications:

20110715 - DWK created by modifying tcohortdat437.h
20110715 - DWK changed class MaxCohortdata43 to 
           class LandAreadata44 
20110715 - DWK deleted int year, int total, and int natchrts
20110715 - DWK changed long carea to float carea
20110715 - DWK added float landfrac
                                 
****************************************************************
************************************************************* */

#ifndef LANDAREADAT44_H
#define LANDAREADAT44_H

class LandAreadata44 
{
   
  public:
          
     LandAreadata44( void );

/* *************************************************************
                      Public Functions
************************************************************* */

// read data structure.
     int get( ifstream& infile );

     int getdel( FILE* infile );

//write data structure.
     void out( ofstream& ofile, 
               const float& col, 
               const float& row, 
               const string& varname, 
               const float& elemntArea,
               const float& landFrac,  
               const string& contnent );
               
     void outdel( ofstream& ofile, 
                  const float& col, 
                  const float& row, 
                  const string& varname, 
                  const float& elementArea,
                  const float& landFrac,  
                  const string& contnent );

          
/* *************************************************************
                     Public Variables
************************************************************* */

     // area covered by grid cell (sq. km)
     float elemntArea;          

     // column or longitude of grid cell (degrees)
     float col;          
     
     // name of continent containing grid cell      
     string contnent;   

     // fraction of grid cell covered by land
     float landFrac;

     // row or latitude of grid cell (degrees)
     float row;               
         
     
     // "LANDAREA"
     string varname;    

  private:

/* *************************************************************
                      Private Variables
************************************************************* */

     int chrtend;
     long curpos;
     long lagpos;

};

#endif

