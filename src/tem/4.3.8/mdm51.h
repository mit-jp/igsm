/* *************************************************************
****************************************************************
MDM51.H - Methane Dynamics Model 


Modifications:

20060104 - DWK created by modifying dattem423e1.hpp
20060104 - DWK changed class TTEM to class MDM51

****************************************************************

References:

Zhuang, Q., J. M. Melillo, D. W. Kicklighter, R. G. Prinn, A. D. 
  McGuire, P. A. Steudler, B. S. Felzer and S. Hu (2004) Methane 
  fluxes between terrestrial ecosystems and the atmosphere at 
  northern high latitudes during the past century: a 
  retrospective analysis with a process-based biogeochemistry 
  model. Global Biogeochemical Cycles, 18, GB3010, 
  doi:10.1029/2004GB002239.

****************************************************************
************************************************************* */

#ifndef MDM51_H
#define MDM51_H

// global constants for MDM
#include "mdmconsts51.hpp"

// for methane production
#include "methanogen51.h" // MDM uses Methanogen class

// for methane oxidation
#include "methanotroph51.h" // MDM uses Methanotroph class

// for diffusion processes
#include "mdmsoil51.h"


class MDM51 
{

  public:

     MDM51();

     enum mdmkey { I_CH4CON,   I_CH4EMI,   I_TOTFLX };


/* **************************************************************
			Public Functions
************************************************************** */

     void stepday( const int& pcmnt,
                   const int& picult,
                   const int& pwet,
                   const double& pfrin,
                   const double& ph2otable,
                   const double& ppctsand,
                   const double& ppctsilt,
                   const double& ppctclay,
                   const double& ppcfldcap,
                   const double& prootz,
                   const double& pph,
                   const double& psuborg,
                   const double& ptsoil );

                   
/* *************************************************************
			 Public Variables
************************************************************* */

     // Methanotrophs in upland areas
     Methanotroph51 dryMethanotroph;

      // The value is the percent area of the cell which is
      //   inundated, multiplied by 100.
     double frin;            

     // cultivation intensity
     double icult;           

     Methanogen51 methanogen;
     
     vector<string> predstr;

     MDMsoil51 soil;

     double wet;             //wetland (type 1-5)

     // Methanotrophs in wetland areas above the water table
     Methanotroph51 wetMethanotroph;
          
};

#endif
