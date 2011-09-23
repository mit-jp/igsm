/* **************************************************************
*****************************************************************
TIGSMVEG44C.H -  Vegetation characteristics used in the 
                 Terrestrial Ecosystem Model (TEM). Based on 
                 tveg44.h

Modifications:

20060126 - DWK created by modifying tveg50b5.h 
20060126 - DWK changed include from tprocessXML51.h to
           tprocessXML437.h
20060126 - DWK changed class Tveg50 to class Tveg43
20060126 - DWK changed inheritance from ProcessXML50 to 
           ProcessXML43
20060126 - DWK deleted public function void setThawPercent()
20060126 - DWK added const double& aot40 and const int& o3flag 
           to public function void updateDynamics()
20060126 - DWK added public functions inline double 
           getFINDOZONE(), inline double getFOZONE(), inline 
           double getFPREVOZONE(), inline double getO3PARA(),
           inline double getO3PARB() and inline double 
           getO3PARC()
20060126 - DWK added private functions double gppxio3() and
           double gppxo3()
20060126 - DWK added const double& fozone to function call of
           double nupxclm()
20060126 - DWK added private double findozone, double fozone,
           double fprevozone, double o3para[MAXCMNT];
           double o3parb[MAXCMNT] and double o3parc[MAXCMNT]
20060221 - DWK deleted private function rq10()
20070201 - DWK changed class Tveg43 to class Tveg44
20070201 - DWK and TWC added const int& agperennial to 
           updateDynamics()
20070201 - DWK added public function setRESPQ10LaRS()
20070201 - DWK added private double tref[MAXCMNT],
           qref[MAXCMNT], alpha[MAXCMNT[MAXCMNT], beta[MAXCMNT],
           and gamma[MAXCMNT]
20080130 - DWK changed include from temconsts43.hpp to
           temigsmconsts44a.hpp
20080130 - DWK changed include from tprocessXML437.h to
           tigsmprocessXML44a.h
20080130 - DWK changed ProcessXML43 to ProcessXML44
20080826 - DWK modified corresponding file tigsmveg44b.cp
20110707 - DWK added inline functions, setLABILEN(), 
           setSTRUCTC(), setSTRUCTN(), getPLANTN(), and
           setPLANTN() 
20110707 - DWK changed include temigsmconsts44a.hpp to
           temigsmconsts44c.hpp
20110707 - DWK changed include from tigsmprocessXML44a.h to
           tigsmprocessXML44c.h
                                                        
****************************************************************
************************************************************* */

#ifndef TIGSMVEG44C_H
#define TIGSMVEG44C_H

// Tveg44 also uses the global constants NUMVEG
#include "temigsmconsts44c.hpp"

// Tveg44 uses Biomass Class
#include "bioms423.hpp"   

// Tveg44 inherits ProcessXML44 Class
#include "tigsmprocessXML44c.h" 

class Tveg44 : public ProcessXML44
{

  public:

     Tveg44();

/* **************************************************************
		 Public Functions
************************************************************** */

     // Limit parameter topt to a biome-specific range of 
     //   air temperatures
     
     void boundTOPT( const int& pcmnt );

     void getecd( ofstream& rflog1 );

     void getecd( const string& ecd );

     void getleafecd( ofstream& rflog1 );

     void getleafecd( const string& ecd );

     void leafinit( ofstream& rflog1 );

     void resetEcds( const int& pcmnt, const double& psiplusc );

     void resetMonthlyFluxes( void );

     void resetNEWTOPT( const int& pcmnt, 
                        const double& tair, 
                        const double& unnrmleaf );

     void resetYrFluxes( void );

     double setGV( const double& eet,
                   const double& pet,
                   const int& moistlim );

     void setRESPQ10( const int& pdcmnt, const double& tair );

     void setRESPQ10LaRS( const int& pdcmnt, const double& tair );

     void setTEMP( const int& pdcmnt, const double& tair );

     void   showecd( const int& pdcmnt );
     
     void   showleaf( const int& pdcmnt );

     void   updateC2N( const int& pdcmnt,
                       const double& yreet,
                       const double& yrpet,
                       const double& currentco2,
                       const double& initco2 );

     void updateDynamics( const int& pdcmnt,
                          const double& co2,
                          const double& aot40,
                          const double& ninput,
                          const double& par,
                          const double& pet,
                          const double& prevmaxpet,
                          const double& eet,
                          const double& prevmaxeet,
                          const double& vegc,
                          const double& structn,
                          const double& labilen,
                          const double& soilh2o,
                          const double& availn,
                          const int& moistlim,
                          const int& nfeed,
                          const int& o3flag,
                          const int& agstate,
                          const int& agperennial,
                          const int& fertflag,
                          const double& ksoil,
                          const double& netnmin,
                          double& agfertn );


     // "Get" and "Set" private variables and parameters
     
     // adjc2n *************************************************

     inline double getADJC2N( void ) { return adjc2n; } 

     inline void setADJC2N( const double& padjc2n ) 
     { 
       adjc2n = padjc2n; 
     };


     // aleaf **************************************************
     
     inline double getALEAF( const int& pcmnt ) 
     { 
       return aleaf[pcmnt]; 
     };

     inline void setALEAF( const double& paleaf, 
                           const int& pcmnt ) 
     { 
       aleaf[pcmnt] = paleaf; 
     };

 
     // alpha *****************************************************

     inline double getALPHA( const int& pcmnt )
     {
       return alpha[pcmnt];
     };

     inline void setALPHA( const double& palpha, const int& pcmnt )
     {
       alpha[pcmnt] = palpha;
     };


     // beta *****************************************************

     inline double getBETA( const int& pcmnt )
     {
       return beta[pcmnt];
     };

     inline void setBETA( const double& pbeta, const int& pcmnt )
     {
       beta[pcmnt] = pbeta;
     };


     // bleaf **************************************************
     
     inline double getBLEAF( const int& pcmnt ) 
     { 
       return bleaf[pcmnt]; 
     };

     inline void setBLEAF( const double& pbleaf, 
                           const int& pcmnt ) 
     { 
       bleaf[pcmnt] = pbleaf; 
     };

 
     // c2n ****************************************************
     
     inline double getC2N( void ) { return c2n; };
     
     inline void setC2N( const double& pc2n ) { c2n = pc2n; };


     // c2na ***************************************************
     
     inline double getC2NA( const int& pcmnt ) 
     {
       return c2na[pcmnt]; 
     };
     
     inline void setC2NA( const double& pc2na, 
                          const int& pcmnt ) 
     { 
       c2na[pcmnt] = pc2na; 
     };


     // c2nb ***************************************************
     
     inline double getC2NB( const int& pcmnt ) 
     { 
       return c2nb[pcmnt]; 
     };
     
     inline void setC2NB( const double& pc2nb, 
                          const int& pcmnt ) 
     { 
       c2nb[pcmnt] = pc2nb; 
     };


     // c2nmin *************************************************
     
     inline double getC2NMIN( const int& pcmnt ) 
     { 
       return c2nmin[pcmnt]; 
     };
     
     inline void setC2NMIN( const double& pc2nmin, 
                            const int& pcmnt ) 
     { 
       c2nmin[pcmnt] = pc2nmin; 
     };


     // cfall **************************************************
     
     inline double getCFALL( const int& pcmnt ) 
     { 
       return cfall[pcmnt]; 
     };

     inline void setCFALL( const double& pcfall, 
                           const int& pcmnt ) 
     { 
       cfall[pcmnt] = pcfall; 
     };


     // cleaf **************************************************
     
     inline double getCLEAF( const int& pcmnt ) 
     { 
       return cleaf[pcmnt]; 
     };

     inline void setCLEAF( const double& pcleaf, 
                           const int& pcmnt ) 
     { 
       cleaf[pcmnt] = pcleaf; 
     };

 
     // cmax ***************************************************
     
     inline double getCMAX( void ) { return cmax; };

     inline void setCMAX( const double& pcmax ) 
     { 
       cmax = pcmax; 
     };


     // cmaxcut ************************************************
     
     inline double getCMAXCUT( const int& pcmnt ) 
     { 
       return cmaxcut[pcmnt]; 
     };

     inline void setCMAXCUT( const double& pcmaxcut, 
                             const int& pcmnt ) 
     { 
       cmaxcut[pcmnt] = pcmaxcut; 
     };


     // cmax1a *************************************************
     
     inline double getCMAX1A( const int& pcmnt ) 
     { 
       return cmax1a[pcmnt]; 
     };

     inline void setCMAX1A( const double& pcmax1a, 
                            const int& pcmnt ) 
     { 
       cmax1a[pcmnt] = pcmax1a; 
     };


     // cmax1b *************************************************
     
     inline double getCMAX1B( const int& pcmnt ) 
     { 
       return cmax1b[pcmnt]; 
     };

     inline void setCMAX1B( const double& pcmax1b, 
                            const int& pcmnt ) 
     { 
       cmax1b[pcmnt] = pcmax1b; 
     };


     // cmax2a *************************************************
     
     inline double getCMAX2A( const int& pcmnt ) 
     { 
       return cmax2a[pcmnt]; 
     };

     inline void setCMAX2A( const double& pcmax2a, 
                            const int& pcmnt ) 
     { 
       cmax2a[pcmnt] = pcmax2a; 
     };


     // cmax2b *************************************************
     
     inline double getCMAX2B( const int& pcmnt ) 
     { 
       return cmax2b[pcmnt]; 
     };

     inline void setCMAX2B( const double& pcmax2b, 
                            const int& pcmnt ) 
     { 
       cmax2b[pcmnt] = pcmax2b; 
     };


     // cneven *************************************************
     
     inline double getCNEVEN( void ) { return cneven; };

     inline void setCNEVEN( const double& pcneven ) 
     { 
       cneven = pcneven; 
     };


     // cnmin **************************************************
     
     inline double getCNMIN( const int& pcmnt ) 
     { 
       return cnmin[pcmnt]; 
     };

     inline void setCNMIN( const double& pcnmin, 
                           const int& pcmnt ) 
     { 
       cnmin[pcmnt] = pcnmin; 
     };


     // cov ****************************************************
     
     inline double getCOV( const int& pcmnt ) 
     { 
       return cov[pcmnt]; 
     };

     inline void setCOV( const double& pcov, const int& pcmnt ) 
     { 
       cov[pcmnt] = pcov; 
     };


     // currentveg *********************************************
     
     inline int getCURRENTVEG( void ) { return currentveg; };

     inline void setCURRENTVEG( const int& ptveg ) 
     { 
       currentveg = ptveg; 
     };
     
     
     // dc2n ***************************************************
     
     inline double getDC2N( void ) { return dc2n; };

     inline void setDC2N( const double& pdc2n ) 
     { 
       dc2n = pdc2n; 
     };
     
     
     // findozone **********************************************
     
     inline double getFINDOZONE( void ) { return findozone; };


     // foliage *************************************************
     
     inline double getFOLIAGE( void ) { return foliage; };


     // fozone *************************************************
     
     inline double getFOZONE( void ) { return fozone; };


     // fpc ****************************************************
     
     inline double getFPC( void ) { return fpc; };


     // fpcmax *************************************************
     
     inline double getFPCMAX( const int& pcmnt ) 
     { 
       return fpcmax[pcmnt]; 
     };

     inline void setFPCMAX( const double& pfpcmx, 
                            const int& pcmnt ) 
     { 
       fpcmax[pcmnt] = pfpcmx; 
     };


     // fprevozone *********************************************
     
     inline double getFPREVOZONE( void ) { return fprevozone; };

     inline void setFPREVOZONE( const double& pfprevozone ) 
     { 
       fprevozone = pfprevozone; 
     };


     // gamma *****************************************************

     inline double getGAMMA( const int& pcmnt )
     {
       return gamma[pcmnt];
     };

     inline void setGAMMA( const double& pgamma, const int& pcmnt )
     {
       gamma[pcmnt] = pgamma;
     };


     // gpp ****************************************************
     
     inline double getGPP( void ) { return gpp; };


     // gpr ****************************************************
     
     inline double getGPR( void ) { return gpr; };


     // gva ****************************************************
     
     inline double getGVA( const int& pcmnt ) 
     { 
       return gva[pcmnt]; 
     };

     inline void setGVA( const double& pgva, 
                         const int& pcmnt ) 
     { 
       gva[pcmnt] = pgva; 
     };


     // ingpp **************************************************
     
     inline double getINGPP( void ) { return ingpp; };


     // initcneven *********************************************
     
     inline double getINITCNEVEN( const int& pcmnt ) 
     { 
       return initcneven[pcmnt]; 
     };

     inline void setINITCNEVEN( const double& pincneven, 
                                const int& pcmnt ) 
     { 
       initcneven[pcmnt] = pincneven; 
     };


     // initleafmx *********************************************
     
     inline double getINITLEAFMX( const int& pcmnt ) 
     { 
      return initleafmx[pcmnt]; 
     };

     inline void setINITLEAFMX( const double& pinleafmx, 
                                const int& pcmnt ) 
     { 
       initleafmx[pcmnt] = pinleafmx; 
     };


     // innpp **************************************************
     
     inline double getINNPP( void ) { return innpp; };


     // inuptake ***********************************************
     
     inline double getINUPTAKE( void ) 
     { 
       return inuptake; 
     };

     
     // kc *****************************************************
     
     inline double getKC( const int& pcmnt ) 
     { 
       return kc[pcmnt]; 
     };

     inline void setKC( const double& pkc, const int& pcmnt ) 
     { 
       kc[pcmnt] = pkc; 
     };


     // ki *****************************************************
     
     inline double getKI( const int& pcmnt ) 
     { 
       return ki[pcmnt]; 
     };

     inline void setKI( const double& pki, const int& pcmnt ) 
     { 
       ki[pcmnt] = pki; 
     };


     // kn1 ****************************************************
     
     inline double getKN1( const int& pcmnt ) 
     { 
       return kn1[pcmnt]; 
     };

     inline void setKN1( const double& pkn1, const int& pcmnt ) 
     { 
       kn1[pcmnt] = pkn1; 
     };


     // kra ****************************************************
     
     inline double getKRA( const int& pcmnt ) 
     { 
       return kra[pcmnt]; 
     };

     inline void setKRA( const double& pkra, const int& pcmnt ) 
     { 
       kra[pcmnt] = pkra; 
     };


     // krb ****************************************************
     
     inline double getKRB( const int& pcmnt ) 
     { 
       return krb[pcmnt]; 
     };

     inline void setKRB( const double& pkrb, const int& pcmnt ) 
     { 
       krb[pcmnt] = pkrb; 
     };


     // kleaf **************************************************
     
     inline double getKLEAFC( const int& pcmnt ) 
     { 
       return kleafc[pcmnt]; 
     };

     inline void setKLEAFC( const double& pkleafc, 
                            const int& pcmnt ) 
     { 
       kleafc[pcmnt] = pkleafc; 
     };


     // labile.nitrogen ****************************************
     
     inline double getLABILEN( void ) 
     { 
       return labile.nitrogen; 
     };

     inline void setLABILEN( const double& plabilen ) 
     { 
       labile.nitrogen = plabilen; 
     };

     // labncon ************************************************
     
     inline double getLABNCON( const int& pcmnt ) 
     { 
       return labncon[pcmnt]; 
     };

     inline void setLABNCON( const double& plabncon, 
                             const int& pcmnt ) 
     { 
       labncon[pcmnt] = plabncon; 
     };


     // lai ****************************************************
     
     inline double getLAI( void ) { return lai; };


     // leaf ***************************************************
     
     inline double getLEAF( void ) { return leaf; };
     

     // leafmxc ************************************************
     
     inline double getLEAFMXC( const int& pcmnt ) 
     { 
       return leafmxc[pcmnt]; 
     };

     inline void setLEAFMXC( const double& pleafmxc, 
                             const int& pcmnt ) 
     { 
       leafmxc[pcmnt] = pleafmxc; 
     };


     // ltrfal.carbon ******************************************  
     
     inline double getLTRFALC( void ) { return ltrfal.carbon; };

     
     // ltrfal.nitrogen ****************************************
     
     inline double getLTRFALN( void ) 
     { 
       return ltrfal.nitrogen; 
     };


     // luptake ************************************************
     
     inline double getLUPTAKE( void ) { return luptake; };


     // minleaf ************************************************
     
     inline double getMINLEAF( const int& pcmnt ) 
     { 
       return minleaf[pcmnt]; 
     };

     inline void setMINLEAF( const double& pminleaf, 
                             const int& pcmnt ) 
     { 
       minleaf[pcmnt] = pminleaf; 
     };

 
     // newleafmx **********************************************
     
     inline double getNEWLEAFMX( void ) { return newleafmx; };

     inline void setNEWLEAFMX( const double& pnewleafmx ) 
     { 
       newleafmx = pnewleafmx; 
     };


     // newtopt ************************************************
     
     inline double getNEWTOPT( void ) { return newtopt; };

     inline void setNEWTOPT( const double& pnewtopt ) 
     { 
       newtopt = pnewtopt; 
     };


     // nfall **************************************************
     
     inline double getNFALL( const int& pcmnt ) 
     { 
       return nfall[pcmnt]; 
     };

     inline void setNFALL( const double& pnfall, 
                           const int& pcmnt ) 
     { 
       nfall[pcmnt] = pnfall; 
     };

 
     // nmax ***************************************************     
     
     inline double getNMAX( void ) { return nmax; };

     inline void setNMAX( const double& pnmax ) 
     { 
       nmax = pnmax; 
     };


     // nmaxcut ************************************************
     
     inline double getNMAXCUT( const int& pcmnt ) 
     { 
       return nmaxcut[pcmnt]; 
     };

     inline void setNMAXCUT( const double& pnmaxcut, 
                             const int& pcmnt ) 
     { 
       nmaxcut[pcmnt] = pnmaxcut; 
     };


     // nmax1a *************************************************
     
     inline double getNMAX1A( const int& pcmnt ) 
     { 
       return nmax1a[pcmnt]; 
     };

     inline void setNMAX1A( const double& pnmax1a, 
                            const int& pcmnt ) 
     { 
       nmax1a[pcmnt] = pnmax1a; 
     };


     // nmax1b *************************************************
     
     inline double getNMAX1B( const int& pcmnt ) 
     { 
       return nmax1b[pcmnt]; 
     };

     inline void setNMAX1B( const double& pnmax1b, 
                            const int& pcmnt ) 
     { 
       nmax1b[pcmnt] = pnmax1b; 
     };


     // nmax2a *************************************************
     
     inline double getNMAX2A( const int& pcmnt ) 
     { 
       return nmax2a[pcmnt]; 
     };

     inline void setNMAX2A( const double& pnmax2a, 
                            const int& pcmnt ) 
     { 
       nmax2a[pcmnt] = pnmax2a; 
     };


     // nmax2b *************************************************
     
     inline double getNMAX2B( const int& pcmnt ) 
     { 
       return nmax2b[pcmnt]; 
     };

     inline void setNMAX2B( const double& pnmax2b, 
                            const int& pcmnt ) 
     { 
       nmax2b[pcmnt] = pnmax2b; 
     };


     // nmobil *************************************************
     
     inline double getNMOBIL( void ) { return nmobil; };


     // npp ****************************************************
     
     inline double getNPP( void ) { return npp; };


     // nresorb ************************************************
     
     inline double getNRESORB( void ) { return nresorb; };


     // nuptake ************************************************
     
     inline double getNUPTAKE( void ) { return nuptake; };

     
     // o3para *************************************************
     
     inline double getO3PARA( const int& pcmnt ) 
     { 
       return o3para[pcmnt]; 
     };

     inline void setO3PARA( const double& po3para, 
                            const int& pcmnt ) 
     { 
       o3para[pcmnt] = po3para; 
     };


     // o3parb *************************************************
     
     inline double getO3PARB( const int& pcmnt ) 
     { 
       return o3parb[pcmnt]; 
     };

     inline void setO3PARB( const double& po3parb, 
                            const int& pcmnt ) 
     { 
       o3parb[pcmnt] = po3parb; 
     };


     // o3parc *************************************************
     
     inline double getO3PARC( const int& pcmnt ) 
     { 
       return o3parc[pcmnt]; 
     };

     inline void setO3PARC( const double& po3parc, 
                            const int& pcmnt ) 
     { 
       o3parc[pcmnt] = po3parc; 
     };


     // plant.nitrogen ****************************************
     
     inline double getPLANTN( void ) 
     { 
       return plant.nitrogen; 
     };

     inline void setPLANTN( const double& pplantn ) 
     { 
       plant.nitrogen = pplantn; 
     };

     // potveg *************************************************
     
     inline int getPOTVEG( void ) { return potveg; };

     inline void setPOTVEG( const int& ptveg ) 
     { 
       potveg = ptveg; 
     };


     // prevunrmleaf *******************************************
     
     inline double getPREVUNRMLEAF( void ) 
     { 
       return prevunrmleaf; 
     };

     inline void setPREVUNRMLEAF( const double& pprevunrmleaf ) 
     { 
       prevunrmleaf = pprevunrmleaf; 
     };


     // prvleafmx **********************************************
     
     inline double getPRVLEAFMX( void ) { return prvleafmx; };

     inline void setPRVLEAFMX( const double& pprvleafmx ) 
     { 
       prvleafmx = pprvleafmx; 
     };


     // qref *****************************************************

     inline double getQREF( const int& pcmnt )
     {
       return qref[pcmnt];
     };

     inline void setQREF( const double& pqref, const int& pcmnt )
     {
       qref[pcmnt] = pqref;
     };


     // raq10a0 ************************************************
     
     inline double getRAQ10A0( const int& pcmnt ) 
     { 
       return raq10a0[pcmnt]; 
     };

     inline void setRAQ10A0( const double& praq10a0, 
                             const int& pcmnt ) 
     { 
       raq10a0[pcmnt] = praq10a0; 
     };


     // raq10a1 ************************************************
     
     inline double getRAQ10A1( const int& pcmnt ) 
     { 
       return raq10a1[pcmnt]; 
     };

     inline void setRAQ10A1( const double& praq10a1, 
                             const int& pcmnt ) 
     { 
       raq10a1[pcmnt] = praq10a1; 
     };


     // raq10a2 ************************************************
     
     inline double getRAQ10A2( const int& pcmnt ) 
     { 
       return raq10a2[pcmnt]; 
     };

     inline void setRAQ10A2( const double& praq10a2, 
                             const int& pcmnt ) 
     { 
       raq10a2[pcmnt] = praq10a2; 
     };


     // raq10a3 ************************************************
     
     inline double getRAQ10A3( const int& pcmnt ) 
     { 
       return raq10a3[pcmnt]; 
     };

     inline void setRAQ10A3( const double& praq10a3, 
                             const int& pcmnt ) 
     { 
       raq10a3[pcmnt] = praq10a3; 
     };


     // rg *****************************************************
     
     inline double getRGRWTH( void ) { return rg; };


     // rm *****************************************************
     
     inline double getRMAINT( void ) { return rm; };


     // sla ****************************************************
     
     inline double getSLA( const int& pcmnt ) 
     { 
       return sla[pcmnt]; 
     };

     inline void setSLA( const double& psla, const int& pcmnt ) 
     { 
       sla[pcmnt] = psla; 
     };


     // strctrl.carbon ***************************************
          
     inline double getSTRUCTC( void ) 
     { 
       return strctrl.carbon; 
     };

     inline void setSTRUCTC( const double& pstructc ) 
     { 
       strctrl.carbon = pstructc; 
     };

     // strctrl.nitrogen ***************************************
          
     inline double getSTRUCTN( void ) 
     { 
       return strctrl.nitrogen; 
     };

     inline void setSTRUCTN( const double& pstructn ) 
     { 
       strctrl.nitrogen = pstructn; 
     };

     // subtype ************************************************
          
     inline int getSUBTYPE( void ) { return subtype; };

     inline void setSUBTYPE( const int& psubtype ) 
     { 
       subtype = psubtype; 
     };


     // suptake ************************************************
          
     inline double getSUPTAKE( void ) { return suptake; };


     // tmax ***************************************************
     
     inline double getTMAX( const int& pcmnt ) 
     { 
       return tmax[pcmnt]; 
     };

     inline void setTMAX( const double& ptmax, 
                          const int& pcmnt ) 
     { 
       tmax[pcmnt] = ptmax; 
     };


     // tmin ***************************************************
     
     inline double getTMIN( const int& pcmnt ) 
     { 
       return tmin[pcmnt]; 
     };

     inline void setTMIN( const double& ptmin, 
                          const int& pcmnt ) 
     { 
       tmin[pcmnt] = ptmin; 
     };


     // topt ***************************************************
     
     inline double getTOPT( void ) { return topt; };

     inline void setTOPT( const double& ptopt ) 
     { 
       topt = ptopt; 
     };


     // toptmax ************************************************
     
     inline double getTOPTMAX( const int& pcmnt ) 
     { 
       return toptmax[pcmnt]; 
     };

     inline void setTOPTMAX( const double& ptoptmax, 
                             const int& pcmnt ) 
     { 
       toptmax[pcmnt] = ptoptmax; 
     };


     // toptmin ************************************************
     
     inline double getTOPTMIN( const int& pcmnt ) 
     { 
       return toptmin[pcmnt]; 
     };

     inline void setTOPTMIN( const double& ptoptmin, 
                             const int& pcmnt ) 
     { 
       toptmin[pcmnt] = ptoptmin; 
     };


     // tref *****************************************************

     inline double getTREF( const int& pcmnt )
     {
       return tref[pcmnt];
     };

     inline void setTREF( const double& ptref, const int& pcmnt )
     {
       tref[pcmnt] = ptref;
     };


     // unleaf12 ***********************************************
     
     inline double getUNLEAF12( const int& pcmnt ) 
     { 
       return unleaf12[pcmnt]; 
     };

     inline void setUNLEAF12( const double& punleaf12, 
                              const int& pcmnt ) 
     { 
       unleaf12[pcmnt] = punleaf12; 
     };


     // unnormleaf *********************************************
     
     inline double getUNNORMLEAF( void ) { return unnormleaf; };


     // plant.carbon *******************************************
     
     inline double getVEGC( void ) { return plant.carbon; };


     // plant.nitrogen *****************************************
      
     inline double getVEGN( void ) { return plant.nitrogen; };

     inline void setVEGN( const double& pvegn )
     {
       plant.nitrogen = pvegn;
     };




/* *************************************************************
		 Public Variables
************************************************************* */


     // Index for community type
     int cmnt;


     // ratio of yrcarbon to yrnitrogen
     double yrc2n;            

     // Annual sum of plant.carbon
     double yrcarbon;          

     // Sum of monthly FPC
     double yrfpc;        

     // annual sum of monthly GPP
     double yrgpp;             

      // Annual GPR
     double yrgpr;            

     // Annual sum of ingpp
     double yringpp;           

     // Annual sum of innpp
     double yrinnpp;           

     // Annual sum of innup
     double yrinnup;           

     double yrinpr;

      // Sum of monthly LAI
     double yrlai;            

     // mean annual normalized leaf phenology
     double yrleaf;      

     // Annual sum of ltrfal.carbon
     double yrltrc;      
     
     // Annual sum of ltrfal.nitrogen      
     double yrltrn;            

     // Annual sum of luptake
     double yrlup;             

     // Annual sum of plant.nitrogen
     double yrnitrogen;        

     // Annual sum of nmobil
     double yrnmobil;          

     // Annual sum of npp
     double yrnpp;             

     // Annual sum of nresorb
     double yrnrsorb;          

     // Annual sum of nuptake
     double yrnup;             

     double yrprod;

     // Annual sum of rg
     double yrrgrowth;

     // Annual sum of rm
     double yrrmaint;

     // Annual sum of labile.nitrogen
     double yrstoren;          
     
      // Annual sum of strctrl.nitrogen
     double yrstructn;       
     
     // Annual sum of suptake
     double yrsup;             

     // Annual sum of unnormleaf
     double yrunleaf;          


  private:

/* **************************************************************
		 Private Functions
************************************************************** */

     double deltaleaf( const int& pdcmnt,
                       const double& eet,
                       const double& prveetmx,
                       const double& prvleaf );

     double gppxclm( const int& pdcmnt,
                     const double& co2,
                     const double& par,
                     const double& temp,
                     const double& gv,
                     const double& leaf,
                     const double& foliage );

     double gppxio3( const double& fozone, 
                     const double& eetpet );

     double gppxo3( const int& pdcmnt,
                    const double& gpp,
                    const double& d40,
                    const double& eetpet );

     double nupxclm( const int& pdcmnt,
                     const double& soilh2o,
                     const double& availn,
                     const double& respq10,
                     const double& ksoil,
                     const double& foliage,
                     const double& fozone );

     double rmxclm( const int& pdcmnt,
                    const double& vegc,
                    const double& respq10 );


/* **************************************************************
		 Private Variables
************************************************************** */
     
     double alleaf;

     // Index for current vegetation type
     int currentveg;

     // Multiplier of indirect ozone effects on GPP
     double findozone;
     
     double foliage;

     // Multiplier of direct ozone effects on GPP
     double fozone;

     // Monthly foliar projective cover
     double fpc;        

     // Previous month's value of fozone
     double fprevozone;

     // Monthly gross primary productivity (GPP)
     double gpp;        

     // Monthly gross plant respiration (rm + rg)
     double gpr;        

      // Initial monthly gross primary productivity
     double ingpp;     

     // Initial net primary productivity
     double innpp;      

     // Initial C/N of biomass production
     double inprodcn;          

     // Initial N uptake by plants
     double inuptake;   

     // Labile plant biomass
     Biomass labile;    

     // Monthly leaf area index
     double lai;        

     // monthly normalized leaf phenology
     double leaf; 

     // Number of annual iterations for determining monthly
     //   phenology
     int leafyrs;
     
     // Monthly litterfall
     Biomass ltrfal;    

     // Monthly N uptake by plants for labile N
     double luptake;    

     // Updated maximum leaf for current year 
     double newleafmx;

     // Updated optimum air temperature for current year
     double newtopt;

     // Monthly N mobilization by plants
     double nmobil;     

     // Monthly net primary productivity (NPP)
     double npp;        

     // Monthly N resorption by plants
     double nresorb;    

     // Monthly N uptake by plants
     double nuptake;    

     // whole plant biomass (structural + labile)
     Biomass plant;     

     // Index for potential vegetation biome type
     int potveg;

     // Unnormalized relative leaf area of previous month
     double prevunrmleaf;

     // Maximum relative leaf area of previous year
     double prvleafmx;

     // Effect of air temperature on plant respiration
     double respq10;

     // Monthly growth respiration
     double rg;         

     // Monthly maintenance respiration
     double rm;         

     // Structural plant biomass
     Biomass strctrl;   

     // Index for vegetation subtype
     int subtype;
     
     // Monthly N uptake by plants for structural N
     double suptake;    

     // Effect of air temperature on GPP
     double temp;

     // Monthly unnormalized leaf phenology
     double unnormleaf; 


/* *************************************************************
		 Private Parameters
************************************************************* */

     // Biome-specific plant respiration parameters

     double alpha[MAXCMNT];
     double beta[MAXCMNT];
     double gamma[MAXCMNT];
     double qref[MAXCMNT];
     double tref[MAXCMNT];

     // Biome-specific vegetation C/N parameters

     double adjc2n;
     double c2n;
     double c2na[MAXCMNT];
     double c2nb[MAXCMNT];
     double c2nmin[MAXCMNT];
     double cnmin[MAXCMNT];
     double dc2n;

     double cneven;
     double initcneven[MAXCMNT];

     double cfall[MAXCMNT];  // proportion of vegetation carbon

     // Biome-specific carbon uptake parameters for function gppxclm

     double cmax;
     double cmaxcut[MAXCMNT];
     double cmax1a[MAXCMNT];
     double cmax1b[MAXCMNT];
     double cmax2a[MAXCMNT];
     double cmax2b[MAXCMNT];

     // Biome-specific respiration parameters for function rmxclm

     double kr;
     double kra[MAXCMNT];
     double krb[MAXCMNT];

     // Biome-specific phenology parameters

     double aleaf[MAXCMNT];
     double bleaf[MAXCMNT];
     double cleaf[MAXCMNT];
     double initleafmx[MAXCMNT];
     double minleaf[MAXCMNT];
     double unleaf12[MAXCMNT];


     // Biome-specific foliage projection cover parameters

     double cov[MAXCMNT];
     double fpcmax[MAXCMNT];
     double sla[MAXCMNT];

     // Biome-specific parameter to describe the sensitivity of GPP
     //   to evapotranspiration
     
     double gva[MAXCMNT];

     // Biome-specific half saturation parameter for function 
     //   gppxclm describing the effects of solar atmospheric 
     //   carbon dioxide concentration on GPP

     double kc[MAXCMNT];

     // Biome-specific half saturation parameter for function 
     //   gppxclm describing the effects of photosybtheically 
     //   active radiation on GPP

     double ki[MAXCMNT];

     double kn1[MAXCMNT];

     // Biome-specific allocation parameters

     double kleafc[MAXCMNT];
     double leafmxc[MAXCMNT];

     double labncon[MAXCMNT];

     double nfall[MAXCMNT];  // proportion of vegetation nitrogen


     // Biome-specific nitrogen uptake parameters for function nupxclm

     double nmax;
     double nmaxcut[MAXCMNT];
     double nmax1a[MAXCMNT];
     double nmax1b[MAXCMNT];
     double nmax2a[MAXCMNT];
     double nmax2b[MAXCMNT];


     // Biome-specific ozone parameters

     double o3para[MAXCMNT];
     double o3parb[MAXCMNT];
     double o3parc[MAXCMNT];


     // Biome-specific parameters for function rq10 to describe the
     //   effect of temperature on plant respiration

     double raq10a0[MAXCMNT];
     double raq10a1[MAXCMNT];
     double raq10a2[MAXCMNT];
     double raq10a3[MAXCMNT];

     // Biome-specific proportion of gpr lost as root respiration

     double rroot[MAXCMNT];

    // Element-specific optimum temperature for GPP

     double tmax[MAXCMNT];
     double tmin[MAXCMNT];
     double topt;
     double toptmax[MAXCMNT];
     double toptmin[MAXCMNT];

};

#endif

