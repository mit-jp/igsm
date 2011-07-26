/* *************************************************************
CLIMATE4TEM.HPP  - Variables shared between TEM and MIT IGSM
************************************************************* */

#ifndef IGSMCLM4TEM_H
#define IGSMCLM4TEM_H


// NOTE: If running TEM "offline", make sure the compiler 
//         directive STANDALONE_TEM is DEFINED below.  If  
//         coupling TEM to the IGSM, make sure STANDALONE_TEM is 
//         NOT DEFINED (i.e. comment out next line)

#include "preproc.h"

// Input Variables to be obtained from other IGSM modules

#include "nemconsts.hpp"

#ifdef STANDALONE_TEM
  static struct
  {
    // Monthly atmospheric CO2
    double co2[MAXNGRD];     

    // Monthly atmospheric O3
    double o3[MAXNGRD][MAX3HRS]; 
  
    // Monthly air temperature
    double temp[MAXNGRD];    
  
    //Daily air temperature
    double daytemp[MAXNGRD][MXMTHDAYS]; 
  
    // Monthly short-wave radiation at surface
    double swrs[MAXNGRD];    
  
    // Monthly precipitation
    double pre[MAXNGRD];    
  
    // Daily duration of storms
    double strmdur[MAXNGRD][MXMTHDAYS];  
  
    // Daily intensity of storms
    double qstrm[MAXNGRD][MXMTHDAYS];    
  
    // Monthly actual evapotranspiration
    double aet[MAXNGRD][MXNCHRTS];     
  
    // Monthly soil moisture in top 1 m
    double sh2o1m[MAXNGRD][MXNCHRTS];  
  
    // Monthly soil moisture in top 2 m
    double sh2o2m[MAXNGRD][MXNCHRTS];  
  
    // Monthly snow water equivalent
    double swe[MAXNGRD][MXNCHRTS];

    // Monthly surface runoff (mm per day)
    double sfr[MAXNGRD][MXNCHRTS];

    // Monthly drainage (mm per day)
    double drn[MAXNGRD][MXNCHRTS];
  
    // Daily soil temperature
    double daytsoil[CLMNLAYERS][MAXNGRD][MXNCHRTS][MXMTHDAYS];  
  
    // Daily soil moisture
    double daysh2o[CLMNLAYERS][MAXNGRD][MXNCHRTS][MXMTHDAYS];   
  
    // Hourly soil moisture
    double hrsh2o[MXNLAYERS][MAXNGRD][MXNCHRTS][MXMTHDAYS][MXDAYHRS]; 

  } climate4tem_; 
#else
  extern "C" struct
  {
    // Monthly atmospheric CO2
    double co2[MAXNGRD];     

    // Monthly atmospheric O3
    double o3[MAXNGRD][MAX3HRS]; 
  
    // Monthly air temperature
    double temp[MAXNGRD];    
  
    //Daily air temperature
    double daytemp[MAXNGRD][MXMTHDAYS]; 
  
    // Monthly short-wave radiation at surface
    double swrs[MAXNGRD];    
  
    // Monthly precipitation
    double pre[MAXNGRD];    
  
    // Daily duration of storms
    double strmdur[MAXNGRD][MXMTHDAYS];  
  
    // Daily intensity of storms
    double qstrm[MAXNGRD][MXMTHDAYS];    
  
    // Monthly actual evapotranspiration
    double aet[MAXNGRD][MXNCHRTS];     
  
    // Monthly soil moisture in top 1 m
    double sh2o1m[MAXNGRD][MXNCHRTS];  
  
    // Monthly soil moisture in top 2 m
    double sh2o2m[MAXNGRD][MXNCHRTS];  
  
    // Monthly snow water equivalent
    double swe[MAXNGRD][MXNCHRTS];

    // Monthly surface runoff (mm per day)
    double sfr[MAXNGRD][MXNCHRTS];

    // Monthly drainage (mm per day)
    double drn[MAXNGRD][MXNCHRTS];
  
    // Daily soil temperature
    double daytsoil[CLMNLAYERS][MAXNGRD][MXNCHRTS][MXMTHDAYS];  
  
    // Daily soil moisture
    double daysh2o[CLMNLAYERS][MAXNGRD][MXNCHRTS][MXMTHDAYS];   
  
    // Hourly soil moisture
    double hrsh2o[MXNLAYERS][MAXNGRD][MXNCHRTS][MXMTHDAYS][MXDAYHRS]; 

  } climate4tem_; 
#endif
#endif

