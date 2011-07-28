/* *************************************************************
****************************************************************
MDMIGSM44A.H - Methane Dynamics Model 


Modifications:

20060104 - DWK created by modifying dattem423e1.hpp
20060104 - DWK changed class TTEM to class MDM51
20080130 - DWK changed include from mdmconsts51.hpp to
           mdmigsmconsts44a.hpp
20080130 - DWK changed include from methanotroph51.h to
           methanotrophigsm44a.h
20080130 - DWK changed include from methanotroph51.h to 
           methanotrophigsm44a.h
20080130 - DWK changed include from mdmsoil51.h to 
           mdmsoiligsm44a.h
20080130 - DWK changed class MDM51 to class MDM44
20080130 - DWK changed public Methanotroph51 dryMethanotroph to
           Methanotroph44 dryMethanotroph
20080130 - DWK changed public Methanogen51 methanogen to
           Methanogen44 methanogen 
20080130 - DWK changed public MDMsoil51 soil to MDMsoil44 soil
20080130 - DWK changed public Methanotroph51 wetMethanotroph to
           Methanotroph44 wetMethanotroph
                                                                 
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

#ifndef MDM44A_H
#define MDM44A_H

// global constants for MDM
#include "mdmigsmconsts44a.hpp"

// for methane production
#include "methanogenigsm44a.h" // MDM uses Methanogen class

// for methane oxidation
#include "methanotrophigsm44a.h" // MDM uses Methanotroph class

// for diffusion processes
#include "mdmsoiligsm44a.h"


class MDM44 
{

  public:

     MDM44();

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
     Methanotroph44 dryMethanotroph;

      // The value is the percent area of the cell which is
      //   inundated, multiplied by 100.
     double frin;            

     // cultivation intensity
     double icult;           

     Methanogen44 methanogen;
     
     vector<string> predstr;

     MDMsoil44 soil;

     double wet;             //wetland (type 1-5)

     // Methanotrophs in wetland areas above the water table
     Methanotroph44 wetMethanotroph;
          
};

#endif
