/* *************************************************************
****************************************************************
TBIOME437.H - object describing general characteristics of
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
                                            
****************************************************************
************************************************************* */

#ifndef TBIOME437_H
#define TBIOME437_H

#include "temconsts43.hpp"
#include "tprocessXML437.h"


class Biome43 : public ProcessXML43
{

  public:

     Biome43( void );

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

