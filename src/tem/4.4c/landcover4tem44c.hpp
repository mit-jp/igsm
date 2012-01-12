/* *************************************************************
LANDCOVER4TEM44C.HPP  - Land cover variables shared between TEM 
                        and the MIT IGSM
************************************************************* */

#ifndef IGSMLULC4TEM_H
#define IGSMLULC4TEM_H


// NOTE: If running TEM "offline", make sure the compiler 
//         directive STANDALONE_TEM is DEFINED below.  If  
//         coupling TEM to the IGSM, make sure STANDALONE_TEM is 
//         NOT DEFINED (i.e. comment out next line)

#include "preproc.h"

// Input Variables to be obtained from other IGSM modules

#include "mittemconsts44c.hpp"

#ifdef STANDALONE_TEM
  static struct
  {
    // Total area of elemnt (e.g. grid cell) that includes area
    //   of all land and water within the element (sq. km)
    //   for a region of MAXNGRD elements
    double elmentArea[MAXNGRD];     

    // Fraction of land area that is covered by a land cover
    //   cohort for MAXCHRTS cohorts within an element for a 
    //   region of MAXNGRD elements
    double fracLandArea[MAXNGRD][CLMMXNVEG]; 
  
    // Fraction of total area in an element that is covered by
    //   land for a region of MAXNGRD elements
    double landFrac[MAXNGRD]; 

    // Initial area of vegetation cohort
    long initialCohortArea[MAXNGRD][CLMMXNVEG]; 
   
  } landcover4tem_; 

#else
  extern "C" struct
  {
    // Total area of elemnt (e.g. grid cell) that includes area
    //   of all land and water within the element (sq. km)
    //   for a region of MAXNGRD elements
    double elmentArea[MAXNGRD];     

    // Fraction of land area that is covered by a land cover
    //   cohort for MAXCHRTS cohorts within an element for a 
    //   region of MAXNGRD elements
    double fracLandArea[MAXNGRD][CLMMXNVEG]; 
  
    // Fraction of total area in an element that is covered by
    //   land for a region of MAXNGRD elements
    double landFrac[MAXNGRD]; 

    // Initial area of vegetation cohort
    long initialCohortArea[MAXNGRD][CLMMXNVEG]; 

  } landcover4tem_; 
#endif
#endif
