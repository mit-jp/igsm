/* *************************************************************
ELMNT437.H - Contains functions to manage elements used in GIS
             (Note: Elements are usually grid cells)
             
Modifications:

20060113 - DWK created by modifying elmnt50b5.h
20060113 - DWK change class Elmnt50 to class Elmnt43
                      
************************************************************** */

#ifndef ELMNT437_H
#define ELMNT437_H

class Elmnt43 
{
   
   public:

     Elmnt43( void );

/* **************************************************************
                 Public Functions
************************************************************** */

     void ask( ofstream& rflog1 );
     
     int coregerr( ofstream& rflog1, 
                   const string& varname1, 
                   const float& col1, 
                   const float& row1, 
                   const string& varname2, 
                   const float& col2, 
                   const float& row2 );
     
     void show( ofstream& rflog1, 
                const float& col, 
                const float& row );
     
     void show( ofstream& rflog1, 
                const float& col, 
                const float& row, 
                const long& totyr, 
                const double& tol );

/* **************************************************************
                 Public Variables
************************************************************** */
         
     // Column or longitude of element
     float col;

     // Area of element
     long carea;

     // Continent location of element
     string contnent;

     // Count of elements in a region
     long count;

     // Mean Elevation of element
     double elev;

     int end;

     long grdcnt;

     long numskip;

     // Row or latitude of element
     float row;

     int stopflag;

     int strtflag;

     int  totyr;


   private:

/* **************************************************************
                 Private Variables
************************************************************** */

     int endflag;
     long numgrids;
};

#endif

