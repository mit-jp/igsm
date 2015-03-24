/* *************************************************************
LULCCHRTDAT44C.H - object to read and write the structure of land 
                    use / land cover cohort data from/to files

20110715 - DWK created by modifying lulcdat44a.h
20110715 - DWK changed class Lulcdata44 to 
           class LulcCohortdata44
20110715 - DWK deleted int isrccohort, int standage, 
           long chrtarea, int potveg, int subtype, int agstate, 
           int agprevstate, int tillflag, int fertflag, 
           int irrgflag, int disturbflag, int disturbmonth, 
           int FRI, double slashpar, double vconvert, 
           double prod10par, double prod100par, double vrespar, 
           and double sconvert 
20110715 - DWK added float fracLandArea
                               
************************************************************* */

#ifndef LULCCHRTDAT44C_H
#define LULCCHRTDAT44C_H

class LulcCohortdata44 
{

  public:

     LulcCohortdata44( void );

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
               const int& year,
               const int& icohort,
               const int& currentveg,
               const float& fracLandArea,
               const string& region );
     
     void outdel( ofstream& ofile, 
                  const float& col, 
                  const float& row, 
                  const string& varname,
                  const int& year,
                  const int& icohort,
                  const int& currentveg,
                  const float& fracLandArea,
                  const string& region );


/* *************************************************************
                     Public Variables
************************************************************* */


     // column or longitude of grid cell (degrees)
     float col;          
     
     // Index of current vegetation type
     int currentveg;

     // fraction of land area covered by land cover cohort 
     //   within a grid cell
     float fracLandArea; 

     // Cohort index
     int icohort;

     // name of region containing cohort
     string region;   
    
     // row or latitude of grid cell (degrees)
     float row;          
   
     // varname = "LULCCHRT" for this data set
     string varname;   

     
     // Year data represents
     int year;


  private:

/* *************************************************************
                      Private Variables
************************************************************* */

     int lulcend;
     long curpos;
     long lagpos;

};

#endif

