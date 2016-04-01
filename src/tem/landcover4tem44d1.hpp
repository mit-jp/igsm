/* *************************************************************
LANDCOVER4TEM44D1.HPP  - Land cover variables shared between TEM 
                        and the MIT IGSM
                        
Modifications:
20140828 - DWK changed double landFrac to long landArea
20140828 - DWK changed double fracLandArea[MAXNGRD][CLMMXNVEG] 
           to long cohortArea[MAXNGRD][CLMMXNVEG] 
20150428 - DWK changed include preproc.h to preproc.hpp
20150428 - DWK changed include from mittemconsts44c.hpp to
           temigsmconstants.hpp
20150428 - DWK added pctsand[], pctsilt[], pctclay[], pH[],
           Ksat[], pporosity[], thickness[] 

************************************************************* */

#ifndef IGSMLULC4TEM_H
#define IGSMLULC4TEM_H


// NOTE: If running TEM "offline", make sure the compiler 
//         directive STANDALONE_TEM is DEFINED below.  If  
//         coupling TEM to the IGSM, make sure STANDALONE_TEM is 
//         NOT DEFINED (i.e. comment out next line)

#include "preproc.hpp"

// Input Variables to be obtained from other IGSM modules

#include "temigsmconstants.hpp"

#ifdef STANDALONE_TEM
  static struct
  {
    // Total area of elemnt (e.g. grid cell) that includes area
    //   of all land and water within the element (sq. km)
    //   for a region of MAXNGRD elements
    double elmentArea[MAXNGRD];     

    // Land area that is covered by a land cover
    //   cohort for MAXCHRTS cohorts within an element for a 
    //   region of MAXNGRD elements
    long cohortArea[MAXNGRD][CLMMXNVEG]; 
  
    // Total area in an element that is covered by
    //   land for a region of MAXNGRD elements
    long landArea[MAXNGRD]; 

    // Initial area of vegetation cohort
    long initialCohortArea[MAXNGRD][CLMMXNVEG]; 
   
   // Percent clay of grid cell's soil texture
   double pctclay[MAXNGRD];
   
   // Percent sand of grid cell's soil texture
   double pctsand[MAXNGRD];
   
   // Percent silt of grid cell's soil texture
   double pctsilt[MAXNGRD];
   
   // Soil pH
   double pH[MAXNGRD];
   
   // Saturated hydraulic conductivity for soil layers
   double Ksat[MAXNGRD][CLMNLAYERS];
   
   // Index of soil layer starting from soil surface
   //  int layer[MAXNGRD][CLMNLAYERS];
   
   // Porosity of soil layers as a fraction of soil volume
   double porosity[MAXNGRD][CLMNLAYERS];
   
   // Thickness of soil layer
   double thickness[MAXNGRD][CLMNLAYERS];
   
  } landcover4tem_; 

#else
  extern "C" struct
  {
    // Total area of elemnt (e.g. grid cell) that includes area
    //   of all land and water within the element (sq. km)
    //   for a region of MAXNGRD elements
    double elmentArea[MAXNGRD];     

    // Land area that is covered by a land cover
    //   cohort for MAXCHRTS cohorts within an element for a 
    //   region of MAXNGRD elements
    long cohortArea[MAXNGRD][CLMMXNVEG]; 
    // Fraction of total area in an element that is covered by
  
    //   land for a region of MAXNGRD elements
    long landArea[MAXNGRD]; 

    // Initial area of vegetation cohort
    long initialCohortArea[MAXNGRD][CLMMXNVEG]; 

   // Percent clay of grid cell's soil texture
   double pctclay[MAXNGRD];
   
   // Percent sand of grid cell's soil texture
   double pctsand[MAXNGRD];
   
   // Percent silt of grid cell's soil texture
   double pctsilt[MAXNGRD];
   
   // Soil pH
   double pH[MAXNGRD];
   
   // Saturated hydraulic conductivity for soil layers
   double Ksat[MAXNGRD][CLMNLAYERS];
   
   // Index of soil layer starting from soil surface
   //  int layer[MAXNGRD][CLMNLAYERS];
   
   // Porosity of soil layers as a fraction of soil volume
   double porosity[MAXNGRD][CLMNLAYERS];
   
   // Thickness of soil layer
   double thickness[MAXNGRD][CLMNLAYERS];

  } landcover4tem_; 
#endif
#endif
