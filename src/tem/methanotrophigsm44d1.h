/* *************************************************************
****************************************************************
METHANOTROPHIGSM44D1.H 

Created by Q. Zhuang, 19/Feb/2003

Modifications:

20060105 - DWK changed class CH4DMOXI to class Methanotroph51
20060105 - DWK renamed MethaneOR() to consumeCH4()
20060105 - DWK changed char ecd[80] to string ecd in 
           getch4oxi() and getecdch4oxi() function calls
20060105 - DWK added include methanotroph51.cpp to bottom of 
           file
20060105 - DWK changed getecdch4oxi() to be getecd()
20060105 - DWK deleted getch4oxi()
20060207 - DWK added double oxi_c[MAXCMNT]
20060207 - DWK added include tprocessXML51.h
20060207 - DWK added inheritance of ProcessXML50
20080130 - DWK changed include from mdmconsts51.hpp to
           mdmigsmconsts44a.hpp
20080130 - DWK changed include from tprocessXML437.h to
           tigsmprocessXML44a.h
20080130 - DWK changed class Methanotroph51 to class
           Methanotroph44
20080130 - DWK changed ProcessXML43 to ProcessXML44
20110707 - DWK changed include from mdmigsmconsts44a.hpp to
           mdmigsmconsts44c.hpp
20110707 - DWK changed include from tigsmprocessXML44a.h to
           tigsmprocessXML44c.h
20150428 - DWK changed include from mdmigsmconsts44c.hpp to
           temigsmconstants.hpp
20150428 - DWK changed include from tigsmprocessXML44c.h 
           to tigsmprocessXML44d1.h
                                                       
****************************************************************
************************************************************* */

#ifndef METHANOTROPH44D_H
#define METHANOTROPH44D_H

#include "temigsmconstants.hpp"

// Methanogen44 inherits ProcessXML44 Class
#include "tigsmprocessXML44d1.h" 

class Methanotroph44 : public ProcessXML44
{

   public:

     Methanotroph44();
     
/* *************************************************************
			Public Functions
************************************************************* */

     double consumeCH4( const int& pcmnt,
                        const double& ch4con, 
                        const double& oxyc, 
                        const double& soilt,
                        const double& soilm, 
                        const double& ehl, 
                        const int& icul );

     void    getecd( ofstream& rflog1 );

     void    getecd( const string& ecd );
                
                                             
     // "Get" and "Set" private variables and parameters    

     // afp ***************************************************

     inline double getAFP( const int& pcmnt ) 
     { 
       return afp[pcmnt]; 
     }
     
     inline void setAFP( const double& pafp, const int& pcmnt ) 
     { 
       afp[pcmnt] = pafp; 
     }


     // consp **************************************************
     
     inline double getCONSP( void ) 
     { 
       return consp; 
     };

     inline void setCONSP( const double& pconsp ) 
     { 
       consp = pconsp; 
     }


     // kc *****************************************************

     inline double getKC( const int& pcmnt ) 
     { 
       return kc[pcmnt]; 
     }
     
     inline void setKC( const double& pkc, const int& pcmnt ) 
     { 
       kc[pcmnt] = pkc; 
     }


     // k0 *****************************************************

     inline double getKO( const int& pcmnt ) 
     { 
       return ko[pcmnt]; 
     }
     
     inline void setKO( const double& pko, const int& pcmnt ) 
     { 
       ko[pcmnt] = pko; 
     }


     // mvmax **************************************************

     inline double getMVMAX( const int& pcmnt ) 
     { 
       return mvmax[pcmnt]; 
     }
     
     inline void setMVMAX( const double& pmvmax, 
                           const int& pcmnt ) 
     { 
       mvmax[pcmnt] = pmvmax; 
     }


     // mvmin **************************************************

     inline double getMVMIN( const int& pcmnt ) 
     { 
       return mvmin[pcmnt]; 
     }
     
     inline void setMVMIN( const double& pmvmin, 
                           const int& pcmnt ) 
     { 
       mvmin[pcmnt] = pmvmin; 
     }


     // mvopt **************************************************

     inline double getMVOPT( const int& pcmnt ) 
     { 
       return mvopt[pcmnt]; 
     }
     
     inline void setMVOPT( const double& pmvopt, 
                           const int& pcmnt ) 
     { 
       mvopt[pcmnt] = pmvopt; 
     }


     // och4q10 ************************************************

     inline double getOCH4Q10( const int& pcmnt ) 
     { 
       return och4q10[pcmnt]; 
     }
     
     inline void setOCH4Q10( const double& poch4q10, 
                             const int& pcmnt ) 
     { 
       och4q10[pcmnt] = poch4q10; 
     }


     // omax ***************************************************

     inline double getOMAX( const int& pcmnt ) 
     { 
       return omax[pcmnt]; 
     }
     
     inline void setOMAX( const double& pomax, const int& pcmnt ) 
     { 
       omax[pcmnt] = pomax; 
     }


     // oxi_c **************************************************

     inline double getOXI_C( const int& pcmnt ) 
     { 
       return oxi_c[pcmnt]; 
     }
     
     inline void setOXI_C( const double& poxi_c, 
                           const int& pcmnt ) 
     { 
       oxi_c[pcmnt] = poxi_c; 
     }


     // oxiref *************************************************

     inline double getOXIREF( const int& pcmnt ) 
     { 
       return oxiref[pcmnt]; 
     }
     
     inline void setOXIREF( const double& poxiref, 
                            const int& pcmnt ) 
     { 
       oxiref[pcmnt] = poxiref; 
     }


   private:

/* *************************************************************
			Private Functions
************************************************************* */

     double EffectCT( const int& icul );

     double EffectMC( const int& pcmnt,
                      const double& ch4con );

     double EffectOXY( const int& pcmnt,
                       const double& oxyc );
                       
     double EffectRX( const double& ehl );

     double EffectSM( const int& pcmnt,
                      const double& soilm );

     double EffectST( const int& pcmnt,
                      const double& soilt );                                              


/* *************************************************************
			Private Variables and Parameters
************************************************************* */

  double afp[MAXCMNT];

  double amax[MAXCMNT];

  double consp; // hourly integration

//  double dailyconsp; // daily consumption

  // Kc is Michaelis-Menten coefficient. (Walter et al., 2000; 
  //   Bender and Conrad, 1992)
  // Walter et al. 2000 set kc as 5uM
  double kc[MAXCMNT];

  // Ko is Michaelis-Menten coefficient, typical value ranged 
  //   from 37 to 200 uM as concentration
  double ko[MAXCMNT];
  
  // volumetric soil moisture (cm-3/cm-3, %)
  double mvmax[MAXCMNT];
  double mvmin[MAXCMNT];
  double mvopt[MAXCMNT];
  
  //och4q10  is a coefficient with a constant, Walter et al., 
  //  set it as 2.0; the values is between 1.4 and 2.1 according 
  //  to Dunfield et al., 1993; Knoblauch, 1994
  double och4q10[MAXCMNT];

  // maximum methane oxidation rate, uMh-1; 5-50 
  //   [Dunfield et al., 1993; Knoblauch, 1994; Krumholz et al., 
  //   1995; Moore and Dalva, 1997; Sundh et al., 1994; Watson 
  //   et al., 1997] 
  // Walter et al. 2000 set omax as 20 uMh-1
  double omax[MAXCMNT];

  double oxi_c[MAXCMNT];

  // reference temperature (oC)  
  double oxiref[MAXCMNT]; 
  
};

#endif
