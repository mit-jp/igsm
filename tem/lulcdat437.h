/* *************************************************************
LULCDAT437.H - object to read and write the structure of land 
                  use / land cover cohort data from/to files

20060113 - DWK created by modifying lulcdat425.h
20060113 - DWK changed class Lulcdata to class Lulcdata43
20060113 - DWK added public int currentveg and int potveg
20060113 - DWK renamed public long carea as long chrtarea
20060113 - DWK renamed public string contnent to string region
20060113 - DWK added public int subtype
20060113 - DWK added public int agprevstate
20060113 - DWK added public int icohort
20060113 - DWK deleted standard string include and 
           lulcdat425.cpp at bottom of file   
20060218 - DWK added public int isrccohort
                   
************************************************************* */

#ifndef LULCDAT437_H
#define LULCDAT437_H

class Lulcdata43 
{

  public:

     Lulcdata43( void );

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
               const int& icohort,
               const int& isrccohort, 
               const long& chrtarea, 
               const int& year,
               const int& potveg,
               const int& currentveg, 
               const int& subtype,
               const int& FRF, 
               const int& agstate, 
               const int& agprevstate, 
               const double& pctag,
               const int& tillflag, 
               const int& fertflag, 
               const int& irrgflag, 
               const double& RAP, 
               const double& slashpar, 
               const double& vconvert, 
               const double& prod10par,
               const double& prod100par, 
               const double& vrespar, 
               const double& sconvert, 
               const string& region );
     
     void outdel( ofstream& ofile, 
                  const float& col, 
                  const float& row, 
                  const string& varname,
                  const int& icohort,
                  const int& isrccohort, 
                  const long& chrtarea, 
                  const int& year,
                  const int& potveg,
                  const int& currentveg, 
                  const int& subtype,
                  const int& FRF, 
                  const int& agstate, 
                  const int& agprevstate, 
                  const double& pctag, 
                  const int& tillflag, 
                  const int& fertflag, 
                  const int& irrgflag, 
                  const double& RAP, 
                  const double& slashpar, 
                  const double& vconvert, 
                  const double& prod10par,
                  const double& prod100par, 
                  const double& vrespar, 
                  const double& sconvert, 
                  const string& region );


/* *************************************************************
                     Public Variables
************************************************************* */

     // flag whether or not grid cell was cultivated the previous
     //   year
     
     int agprevstate;        

     // flag whether or not grid cell is cultivated
     int agstate;        

     // area of a cohort within a grid cell
     long chrtarea; 

     // column or longitude of grid cell (degrees)
     float col;          
                        
     // Index of current vegetation type
     int currentveg;
          
     // Flag to indicate whether cropland in grid cell 
     //   experiences optimum fertilizer application 
     int fertflag;
     
     // Fire Return Frequency
     int FRF;            

     // Cohort index
     int icohort;

     // Flag to indicate whether cropland in grid cell 
     //   experiences optimum irrigation 
     int irrgflag;

     // Index of cohort that was source of current cohort
     int isrccohort;
          
     // percent of grid cell covered with agriculture
     double pctag;       
     
     // Index of potential vegetation biome type
     int potveg;

     // Proportion of vegetation carbon extracted for 
     // 10-year wood products
     double prod10par;
     
     // Proportion of vegetation carbon extracted for 
     // 100-year wood products
     double prod100par;
     
     // relative agricultural production
     double RAP;        
 
     // name of region containing cohort
     string region;   
    
     // row or latitude of grid cell (degrees)
     float row;          
   
     // Proportion of soil organic carbon lost 
     // during a disturbance
     double sconvert;

     // Proportion of vegetation left as slash when
     // ecosystem is disturbed
     double slashpar;

     // Index for vegetation subtype
     int subtype;
     
     // Flag to indicate whether cropland in grid cell 
     //   is tilled ( = 1 ) or not ( = 0 ) 
     int tillflag;
     
     // varname = "AGRICUL" for this data set
     string varname;   

     // Proportion of vegetation carbon lost 
     // during a disturbance
     double vconvert;
     
     double vrespar;
     
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

