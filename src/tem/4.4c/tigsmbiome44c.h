/* *************************************************************
****************************************************************
TIGSMBIOME44C.H - object describing general characteristics of
                  vegetation mosaic used in the Terrestrial 
                  Ecosystem Model (TEM)

Modifications:

20060114 - DWK created by modifying tbiome50b5.h
20060114 - DWK changed include from temconsts51.hpp to 
           temconsts43.hpp
20060114 - DWK changed include from tprocessXML50b5.h to
           tprocessXML437.h
20060114 - DWK changed inheritance from ProcessXML50 to 
           ProcessXML43
20080130 - DWK changed include from temconst43.hpp to
           temigsmconsts44a.hpp
20080130 - DWK changed include from tprocessXML437.h to
           tigsmprocessCXML44a.h
20080130 - DWK changed class Biome43 to class Biome44
20080130 - DWK changed ProcessXML43 to ProcessXML44
20110707 - DWK changed include from temigsmconsts44a.hpp
           to temigsmconst44c.hpptigsmprocessCXML44a.h
20110707 - DWK changed include from tigsmprocessCXML44a.h to
           tigsmprocessCXML44c.h
                                                       
****************************************************************
************************************************************* */

#ifndef TIGSMBIOME44C_H
#define TIGSMBIOME44C_H

#include "temigsmconsts44c.hpp"
#include "tigsmprocessXML44c.h"


class Biome44 : public ProcessXML44
{

  public:

     Biome44( void );

 /* ************************************************************
		 Public Functions
************************************************************* */

     int   getCommunityType( const int& tveg );
     int    getVegMosaic( const int& tveg );
     double getVegSubarea( const int& tveg,
                           const int& dtype,
                           const int& carea );
      int   getVegSubtype( const int& tveg, const int& dtype );
     void   getvtype( ofstream& rflog1 );
     void   getvtype( const string& ecd );

/* *************************************************************
		 Public Variables
************************************************************* */

     // biome type or ecozone (categorical data)
     int temveg;

     // vegetation community type (categorical data)
     int cmnt;

     // number of community types in a vegetation mosaic
     int numtype[NUMVEG];

     // community types of a vegetation mosaic
     int subtype[NUMVEG][NUMMSAC];

     // percent coverage of a community type in a vegetation
     //   mosaic
     double pcttype[NUMVEG][NUMMSAC];

     //Description of a vegetation community type

     string cmnt_name;

     // Area covered by a vegetation community type

     double subarea;

};

#endif

