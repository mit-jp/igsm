; $Id: find_cells_by_country.pro,v 1.2 2007/07/27 18:57:29 bmy Exp $
;-----------------------------------------------------------------------
;+
; NAME:
;        FIND_CELLS_BY_COUNTRY
;
; PURPOSE:
;        Returns an index array which can be used to determine 
;        which CTM grid boxes lie within a given country.
;
; CATEGORY:
;        GAMAP Utilities
;
; CALLING SEQUENCE:
;        RESULT = FIND_CELLS_BY_COUNTRY( COUNTRYID, GRIDINFO [, Keywords] )
;
; INPUTS:
;        COUNTRYID -> Name of the country you are interested in, OR
;             the country ID number ID from "countries.table".  
;
;        GRIDINFO -> Structure from CTM_GRID which defines the output
;             CTM model grid in which you wish to find a given country.
;
; KEYWORD PARAMETERS:
;        /MAXIMIZE -> Set this switch to search for all grid boxes that
;             contain any portion of the specified country.  The
;             default is to determine the specified country by the
;             center of the grid box.
;
;        /INDEX -> Set this switch to return RESULT as an 1-D index 
;             vector (i.e. similar to output from the WHERE command). 
;
; OUTPUTS:
;        RESULT -> Integer index array for OUTGRID.  Grid boxes where 
;             RESULT[I,J] = 1 designate the desired country.  If the
;             /INDEX flag is set then RESULT will be an 1-D index
;             vector (i.e. similar to output from the WHERE command). 
;       
; SUBROUTINES:
;        Internal Subroutines:
;        ============================
;        GET_COUNTRY_NUMBER (function)
;
;        External Subroutines Required:
;        ===========================================
;        CHKSTRU  (function)    CTM_GRID (function)
;        CTM_TYPE (function)    DATATYPE (function)
;        STRBREAK (function) 
;
; REQUIREMENTS:
;        None
;
; NOTES:
;        (1) Requires the following input files:
;            (a) countries.table                
;            (b) countries.generic.025x025.gif
;
;        (2) The search algorithm is brute-force (i.e. FOR loops
;            over lat & lon).  Maybe in a later version we can 
;            optimize this w/ IDL array notation.
;
; MODIFICATION HISTORY:
;  tmf & bmy, 01 Jul 2006: GAMAP VERSION 2.05
;                          - Initial version
;  bmy & phs, 27 Jul 2007: GAMAP VERSION 2.10
;                          - Use FILE_WHICH to find countries.table
;                          - Use FILE_WHICH to find countries.gif 
;
;-
; Copyright (C) 2006-2007, Martin Schultz,
; Bob Yantosca and Philippe Le Sager, Harvard University
; This software is provided as is without any warranty whatsoever. 
; It may be freely used, copied or distributed for non-commercial 
; purposes. This copyright notice must be kept with any copy of 
; this software. If this software shall be used commercially or 
; sold as part of a larger package, please contact the author.
; Bugs and comments should be directed to bmy@io.as.harvard.edu
; or phs@io.harvard.edu with subject "IDL routine find_cells_by_country"
;-----------------------------------------------------------------------


function Get_Country_Number, CountryID
   
   ;====================================================================
   ; Internal function GET_COUNTRY_NUMBER reads the lookup table
   ; "countries.table" and returns the country number that corresponds
   ; to a given country name. (bmy, 10/4/06)
   ;====================================================================

   ; Initialize
   Line        = ''
   CountryName = StrUpCase( StrTrim( CountryID, 2 ) )
   CountryNum  = -1

   ; Use FILE_WHICH to look for the "countries.table" file, and
   ; failing that, in the directories specified in the !PATH variable
   FileName = File_Which( 'countries.table', /Include_Current_Dir )
   FileName = Expand_Path( FileName )

   ; Open the countries file
   Open_File, FileName, Ilun, /Get_LUN
  
   ; Skip header line
   ReadF, Ilun, Line

   ; Loop thru the file
   while ( not EOF( Ilun ) ) do begin

      ; Read each line
      ReadF, Ilun, Line
      
      ; Split the line by tabs
      Result = StrBreak( Line, Byte(9B) )

      ; Full Country name is in the 6th column of the file
      ThisCountry = StrUpCase( StrTrim( Result[5], 2 ) )

      ; Return the country number if there's a match
      if ( ThisCountry eq CountryName ) then begin
         CountryNum = Long( Result[0] )
         goto, Quit1
      endif

      ; Undefine stuff for next pass 
      Undefine, Result
    
   endwhile

Quit1:

   ; Close file
   Close,    Ilun
   Free_LUN, Ilun

   ; Return 
   return, CountryNum
end

;------------------------------------------------------------------------------

function Find_Cells_By_Country, CountryID,   GridInfo, $
                                Index=Index, Maximize=Maximize, _EXTRA=e

   ;====================================================================
   ; Initialization
   ;====================================================================

   ; External functions 
   FORWARD_FUNCTION ChkStru, CTM_Type, CTM_Grid, DataType, StrBreak

   ; Fields for GRIDINFO test
   Test = [ 'IMX', 'JMX', 'XMID', 'YMID', 'XEDGE', 'YEDGE' ]

   ; Error check inputs
   Index    = Keyword_Set( Index    )
   Maximize = Keyword_Set( Maximize )
   if ( N_Elements( CountryID ) ne 1  ) then Message, 'COUNTRYID not passed!'
   if ( not ChkStru( GridInfo, Test ) ) then Message, 'Invalid GRIDINFO!'

   ; Get the country number (match to name if necessary)
   if ( DataType( CountryID, /Name ) eq 'STRING' )      $
      then CountryNum = Get_Country_Number( CountryID ) $
      else CountryNum = Long( CountryID )

   ; Exit if country isn't found
   if ( CountryNum lt 0 ) then Message, 'Invalid COUNTRYID!'

   ; Input grid (0.25 x 0.25)
   FineType                = CTM_Type( 'GENERIC', Resolution=[0.25,0.25] )
   FineGrid                = CTM_Grid( FineType )

   ; Define some 0.25 x 0.25 parameters 
   IMIN                    = FineGrid.XEdge[0]
   IMAX                    = FineGrid.XEdge[FineGrid.IMX]
   JMIN                    = FineGrid.YEdge[1]
   IMX                     = FineGrid.IMX
   JMX                     = FineGrid.JMX

   ; Get the TEMPWORLD array from a GIF file
   FileName = File_Which( 'countries.025x025.gif', /Include_Current_Dir )
   FileName = Expand_Path( FileName )   
   Read_GIF, FileName, TempWorld

   ; Patch a strip of 20 grid longitude to the east of WORLD matrix 
   ; for cross-dateline grids
   TempWorld               = Long( Temporary( Tempworld ) )
   Temp                    = Reform( TempWorld[0L:29L,*] )
   World                   = LonArr( IMX + 30L, JMX )
   World[0L:IMX-1L,*]      = Reform( TempWorld )
   World[IMX:IMX+30L-1L,*] = Reform( Temp )

   ; Array for output (on model grid)
   CountryFlag             = IntArr( GridInfo.IMX, GridInfo.JMX )

   ;====================================================================
   ; Determine which output grid boxes are w/in our desired country
   ;====================================================================

   if ( Maximize eq 1 ) then begin

      ;-----------------------------------------------------------------
      ; If /MAXIMIZE is set, then consider a match if any part of the
      ; desired country lies w/in the output grid box (I,J)
      ;-----------------------------------------------------------------

      ; Loop over output grid boxes
      For J = 0L, GridInfo.JMX-1L do begin
      For I = 0L, GridInfo.IMX-1L do begin

         ; Lon edges of the output grid box
         LonEdge1 = GridInfo.XEdge[I]
         LonEdge2 = GridInfo.XEdge[I+1]

         ; Lat edges of the output grid box
         LatEdge1 = GridInfo.YEdge[J]
         LatEdge2 = GridInfo.YEdge[J+1]

         ; If the western most grid edge is westward of IMIN, add 360. 
         If ( LonEdge1 lt IMIN ) then begin
            LonEdge1 = LonEdge1 + 360.0
            LonEdge2 = LonEdge2 + 360.0
         endif

         ; Find the 0.25 x 0.25 longitude edges that 
         ; correspond to the output grid box edges
         II1 = Floor( ( LonEdge1 - IMIN ) / FineGrid.DI )
         II2 = Floor( ( LonEdge2 - IMIN ) / FineGrid.DI )

         ; Find the 0.25 x 0.25 lat edge corresponding to LATEDGE1
         if ( LatEdge1 lt JMIN )                                     $
            then JJ1 = 0                                             $
            else JJ1 = 1 + Floor( ( LatEdge1 - JMIN ) / FineGrid.DJ )

         ; Find the 0.25 x 0.25 lat edge corresponding to LATEDGE2
         if ( LatEdge2 lt JMIN )                                     $
            then JJ2 = 0                                             $
            else JJ2 = 1 + Floor( ( LatEdge2 - JMIN ) / FineGrid.DJ )

         ; Fix: since IDL starts at zero, the max value 
         ; of JJ1 or JJ2 should be FineGrid.JMX-1.
         if ( JJ1 eq FineGrid.JMX ) then JJ1 = FineGrid.JMX - 1
         if ( JJ2 eq FineGrid.JMX ) then JJ2 = FineGrid.JMX - 1      

         ; COVERAGE is the subset of WORLD between the 0.25 x 0.25
         ; lon and lat edges II1:II2 and JJ1:JJ2   
         Coverage = Reform( World[ II1:II2, JJ1:JJ2 ] )

         ; Test if any part of COVERAGE lies w/in our desired country
         Ind = Where( Coverage[*,*] eq CountryNum )
         if ( ind[0] ge 0 ) then CountryFlag[I,J] = 1

         ; Undefine stuff for next pass
         Undefine, Coverage
         UnDefine, Ind

      endfor
      endfor

   endif else begin

      ;-----------------------------------------------------------------
      ; Otherwise, consider a match if the center of the output grid
      ; box (I,J) lies w/in the desired country.  This is the default.
      ;-----------------------------------------------------------------

      ; Loop over output grid boxes
      for J = 0L, GridInfo.JMX-1L do begin  
      for I = 0L, GridInfo.IMX-1L do begin

         ; Lon and lat of the center of the model grid
         Lon = GridInfo.XMid[I]
         Lat = GridInfo.YMid[J]

         ; Get corresponding lon index of 0.25 x 0.25 grid
         II = Floor( ( Lon - IMIN ) / FineGrid.DI )

         ; Get corresponding lat index of 0.25 x 0.25 grid
         if ( Lat lt JMIN )                                    $
            then JJ = 0                                        $
            else JJ = 1 + Floor( ( Lat - JMIN ) / FineGrid.DJ )

         ; Test if the output grid box is w/in our desired country
         if ( World[II,JJ] eq CountryNum ) then CountryFlag[I,J] = 1

      endfor
      endfor

   endelse

   ;====================================================================
   ; Cleanup and quit
   ;====================================================================
Quit:

   ; Return as 1-D index vector or 2-D array
   if ( Index ) then begin
      Ind = Where( CountryFlag gt 0 ) 
      return, Ind
   endif else begin
      return, CountryFlag
   endelse
end


