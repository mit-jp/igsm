; $Id: ctm_grid.pro,v 1.3 2007/08/16 14:43:13 bmy Exp $
;-----------------------------------------------------------------------
;+
; NAME:
;        CTM_GRID  (function)
;
; PURPOSE:
;        Set-up of the model grid for various models. While the
;        horizontal grid information can be computed from a few
;        basic parameters, the vertical structure is hardwired
;        for each model seperately. Currently, the following models
;        are supported: GEOS1, GEOS_STRAT, GEOS-3, GEOS-4, GEOS-5,
;        GISS_II, GISS_II_PRIME, FSU, FVCCM, and MATCH.
;       
; CATEGORY
;        GAMAP Utilities, GAMAP Models & Grids, Structures
;
; CALLING SEQUENCE:
;        RESULT = CTM_GRID( MTYPE [, Keywords ] )
;
; INPUTS:
;        MTYPE --> a MODELINFO structure as returned from function 
;             CTM_TYPE.  This structure  must contain the following tags: 
;             resolution (a 2 element vector with DI, DJ), halfpolar, 
;             center180, name, nlayers, ptop, psurf, hybrid.)
; 
; OUTPUT:
;        RESULT -> A structure, scalor, or vector variable depending 
;             on the output keywords below.
;
; KEYWORD PARAMETERS:
;        PSURF -> specify surface pressure in mb. Overrides 
;            value passed in through modelinfo structure.
;
;        The following keywords define output options. If none of these
;        is set, a structure is returned that contains all parameters.
;	   IMX   (int   ) -> Maximum I (longitude) dimension  [alt: NLON]
;          JMX   (int   ) -> Maximum J (latitude ) dimension  [alt: NLAT]
;          DI    (flt   ) -> Delta-I interval between grid box edges
;          DJ    (flt   ) -> Delta-J interval between grid box edges
;          YEDGE (fltarr) -> Array of latitude  edges
;          YMID  (fltarr) -> Array of latitude  centers
;          XEDGE (fltarr) -> Array of longitude edges   
;          XMID  (fltarr) -> Array of longitude centers
;
;          LMX   (int)    -> Maximum altitude level (layers)  [alt: NLAYERS]
;          SIGMID (fltarr)-> Array of sigma center values
;          SIGEDGE (..)   -> Array of sigma edges
;          ETAMID  (..)   -> Array of ETA center values
;          ETAEDGE (..)   -> Array of ETA edge values
;          PMID    (..)   -> Array of approx. pressure values for sigma centers
;          PEDGE   (..)   -> dto. for sigma edges
;          ZMID    (..)   -> Array of approx. mean altitudes for sigma centers
;          ZEDGE   (..)   -> dto. for sigma edges
;
;        /NO_VERTICAL --> do not return vertical grid info in structure
;             In this case the MTYPE structure only needs to contain
;             resolution, halfpolar and center180. This keyword is ignored if 
;             a specific (vertical) output option is requested.
;
; SUBROUTINES:
;        External Subroutines Required:
;        ===============================
;        USSA_ALT (function)
;        GETSIGMA (function)
;        GETETA   (function)
;
; REQUIREMENTS:
;        Best if used with function CTM_TYPE.  Also requires functions
;        GETSIGMA and GETETA for definition of vertical layers.
;
; NOTES:
;        This routine is not very efficient in that it always computes 
;        all the available information. But since it will not be
;        called too often and does not handle large amounts of data, 
;        we won't worry about computational efficiency here.
;
; EXAMPLE:
;        MTYPE = CTM_TYPE( 'GEOS1' )
;        MGRID = CTM_GRID( MTYPE )
;
;        ; This will define the following structure (help,mgrid,/stru)
;
;        ** Structure <10323808>, 17 tags, length=1624, refs=1:
;           IMX             INT             72
;           JMX             INT             46
;           DI              FLOAT           5.00000
;           DJ              FLOAT           4.00000
;           XEDGE           FLOAT     Array[73]     (-177.500, -172.500, ...)
;           YEDGE           FLOAT     Array[47]     (-90, -88, -84, ...)
;           XMID            FLOAT     Array[72]     (-180.0, -175.0, ...)
;           YMID            FLOAT     Array[46]     (-89, -86, -82, ...)
;           LMX             INT             20
;           SIGMID          FLOAT     Array[20]     (0.993936, 0.971301, ...)
;           SIGEDGE         FLOAT     Array[21]     (1.00000,  0.987871, ...)
;           ETAMID          FLOAT     Array[20]     (all zeroes)
;           ETAEDGE         FLOAT     Array[21]     (all zeroes)
;           PMID            FLOAT     Array[20]     (980.082, 957.990, ...)
;           PEDGE           FLOAT     Array[21]     (986.000, 974.162, ...)
;           ZEDGE           FLOAT     Array[21]
;           ZMID            FLOAT     Array[20]
;
;        ; Or, with the use of output keywords:
;             PRINT, CTM_GRID( CTM_TYPE( 'GISS_II' ), /PMID )
;
;        ; IDL will print
;             986.000      935.897      855.733      721.458      551.109 
;             390.781      255.503      150.287      70.1236      10.0000
;
;        ; A more awkward example (see yourself):
;             HELP, CTM_GRID( { RESOLUTION : [3.6,3.0],           $
;                               HALFPOLAR  : 0,                   $
;                               CENTER180  : 0 }, /NO_VERT), /STRU
;
; MODIFICATION HISTORY:
;        bmy, 19 Aug 1997: VERSION 1.00
;        bmy, 24 Sep 1997: VERSION 1.01
;        mgs, 26 Feb 1998: Version 2.00  - rewritten as a function
;        mgs, 27 Feb 1998: - added vertical information
;        mgs, 02 Mar 1998: - better defined interface with CTM_MODEL_TYPE
;        bmy, 07 Apr 1998: - Renamed 
;        mgs, 24 Apr 1998: - changed structure to named structure
;        mgs, 04 May 1998: - changed back because of conflicting sizes
;        mgs, 07 May 1998: - Renamed to CTM_GRID
;                          - x coordinates now start with -182.5 for 
;                            center180 grids
;        bmy, 19 Jun 1998: - now uses USSA_ALT to compute altitudes
;                            from pressure coordinates
;                          - fixed some comments
;                          - added FORWARD_FUNCTION statement
;        mgs, 30 Jun 1999: - added PSURF keyword
;        bmy, 27 Jul 1999: GAMAP VERSION 1.42
;                          - now can compute pressure levels and
;                            edges for hybrid sigma-pressure grids
;                          - a few cosmetic changes
;        bmy, 03 Jan 2000: - more cosmetic changes
;        bmy, 20 Apr 2000: GAMAP VERSION 1.45
;                          - now returns correct YMID values for FSU grid
;        bmy, 15 Sep 2000: GAMAP VERSION 1.46
;                          - fixed bug for computing YMID for grids
;                            w/o halfpolar boxes.  This also fixes the
;                            previous problem w/ the FSU grid.
;        bmy, 03 Jul 2001: GAMAP VERSION 1.48
;                          - If MTYPE.NLAYERS = 0, then create a
;                            return structure w/o vertical level info
;        bmy, 06 Nov 2001: GAMAP VERSION 1.49
;                          - added ETAMID, ETAEDGE keywords
;                          - added ETAMID, ETAEDGE tags to return structure
;                          - now calls GETETA to return ETA coordinates
;                            for hybrid models (such as GEOS-4/fvDAS)
;                          - updated comments
;        bmy, 18 Oct 2002: GAMAP VERSION 1.52
;                          - deleted obsolete commented-out code
;        bmy, 04 Nov 2003: GAMAP VERSION 2.01
;                          - Use STRPOS to test for GEOS4 or 
;                            GEOS4_30L model names
;                          - Now treat GISS_II_PRIME 23-layer model
;                            as a hybrid grid instead of using the
;                            obsolete "fake" formulation.
;        bmy, 28 Jun 2006: GAMAP VERSION 2.05
;                          - bug fix for multi-level GENERIC grids
;  bmy & phs, 16 Aug 2007: GAMAP VERSION 2.10
;                          - now compute mXEDGE, mXMID, mYEDGE, mYMID
;                            as double precision, and cast back to
;                            float, so that we get correct values for
;                            high-res grids like 0.5 x 0.666,
;                          - cosmetic changes
;                          - Special handling for GEOS-5 
;                          - Now use USAGE, ROUTINE_NAME() instead of
;                            function USE_CTM_GRID to display options
;                          - Now return IMX, JMX as type LONG
;
;-
; Copyright (C) 1999-2007, Martin Schultz,
; Bob Yantosca and Philippe Le Sager, Harvard University
; This software is provided as is without any warranty whatsoever. 
; It may be freely used, copied or distributed for non-commercial 
; purposes. This copyright notice must be kept with any copy of 
; this software. If this software shall be used commercially or 
; sold as part of a larger package, please contact the author.
; Bugs and comments should be directed to bmy@io.as.harvard.edu
; or phs@io.harvard.edu with subject "IDL routine ctm_grid"
;-----------------------------------------------------------------------


function CTM_GRID, MTYPE,                                             $
                   PSURF=PSURF,     IMX=IMX,         JMX=JMX,         $
                   NLON=NLON,       NLAT=NLAT,       DI=DI,           $
                   DJ=DJ,           XEDGE=XEDGE,     YEDGE=YEDGE,     $
                   XMID=XMID,       YMID=YMID,       LMX=LMX,         $
                   NLAYERS=NLAYERS, SIGMID=SIGMID,   SIGEDGE=SIGEDGE, $
                   ETAMID=ETAMID,   ETAEDGE=ETAEDGE, PMID=PMID,       $
                   PEDGE=PEDGE,     ZMID=ZMID,       ZEDGE=ZEDGE,     $
                   NO_VERTICAL=NO_VERTICAL,          HELP=HELP
 
   ; Pass external functions
   FORWARD_FUNCTION USSA_Alt, GetSigma, GetEta

   ;====================================================================
   ;  Error checking
   ;====================================================================
   on_error, 2

   if (keyword_set(HELP)) then begin
      usage, Routine_Name()
      return,-1
   endif

   if (n_params() lt 1) then begin
      print, 'Too few parameters...'
      Usage, Routine_Name()
      return,-1
   endif

   ; check if MTYPE has the correct structure type 
   fields = ['resolution','halfpolar','center180']

   notvert = keyword_set(NO_VERTICAL) and $
      not(keyword_set(LMX)    or keyword_set(SIGEDGE) or $
          keyword_set(SIGMID) or keyword_set(PEDGE) or $
          keyword_set(PMID)   or keyword_set(ZEDGE) or $
          keyword_set(ZMID) )
   
   if (not notvert) then $
      fields = [ fields, 'name','nlayers','ptop','psurf' ]

   if (not chkstru(MTYPE,fields,/VERBOSE)) then begin
      Usage, Routine_Name()
      return,-1
   endif

   ; extra check: see if resolution is 2 element vector
   if (n_elements(MTYPE.resolution) ne 2) then begin
      print,'CTM_GRID: ** MODEL resolution invalid !! **'
      return,-1
   endif

   ; If there are no vertical layers specified in MTYPE, then do 
   ; not include vertical layers in the return structure (bmy, 7/3/01)
   if ( MType.NLayers eq 0 ) then NotVert = 1

   ;====================================================================
   ; Extract resolution info from structure
   ;====================================================================

   mDI = MTYPE.resolution(0)    ; longitude interval
   mDJ = MTYPE.resolution(1)    ; latitude interval

   ;====================================================================
   ; Compute number of grid boxes from resolution - 
   ; take care of polar boxes
   ;
   ; NOTE: ATTENTION needed if model is added that does not have 
   ;       half size polar boxes!
   ;====================================================================

   ;---------------------------------------------------
   ; Prior to 8/16/07:
   ; Make IMX, JMX variables longwords (bmy, 8/16/07)
   ;mIMX  =  fix(360./mDI)
   ;mJMX  =  fix(180./mDJ) + MTYPE.halfpolar
   ;---------------------------------------------------
   mIMX = Long( 360. / mDI )
   mJMX = Long( 180. / mDJ ) + MTYPE.halfpolar

   ;====================================================================
   ; Compute longitude and latitude edge and center vectors
   ; NOTE: Need to compute in DOUBLE PRECISION and cast back to FLOAT
   ;====================================================================
   mXEDGE       = DIndGen( mIMX+1 )*mDI - 180d0 - ( mDI/2d0 * MTYPE.Center180 )
   mXEDGE       = float( mXEDGE )

   mXMID        = mXEDGE - mDI/2d0
   mXMID        = mXMID(1:*)
   mXMID        = Float( mXMID )

   mYEDGE       = DIndGen( mJMX+1 )*mDJ - 90d0 - ( mDJ/2d0 * MTYPE.HalfPolar )
   mYEDGE(0)    = -90.
   mYEDGE(mJMX) =  90.
   mYEDGE       = Float( mYEDGE )
  
   mYMID        = DIndGen( mJMX )*mDJ - 90d0
   mYMID        = Float( mYMID )

   ; Bug fix for grids w/o halfpolar boxes (bmy, 9/15/00)
   if ( MTYPE.HalfPolar ) then begin
      mYMID(0)      = ( mYEDGE(0) + mYEDGE(1) ) / 2.
      mYMID(mJMX-1) = - mYMID(0)
   endif else begin
      ; For grids w/o halfpolar boxes, add 1/2 of the 
      ; box size to MYMID as computed above (bmy, 9/15/00)
      mYMID = mYMID + ( mYEDGE(1) - mYEDGE(0) ) / 2.
   endelse
   
   ; numerical fix
   ind = where(abs(mXEDGE) lt 1.e-5)
   if (ind(0) ge 0) then mXEDGE(ind) = 0.
   ind = where(abs(mXMID) lt 1.e-5)
   if (ind(0) ge 0) then mXMID(ind) = 0.
   ind = where(abs(mYEDGE) lt 1.e-5)
   if (ind(0) ge 0) then mYEDGE(ind) = 0.
   ind = where(abs(mYMID) lt 1.e-5)
   if (ind(0) ge 0) then mYMID(ind) = 0.

   ;====================================================================
   ; Get vertical structure information
   ;====================================================================

   ; no vertical information requested
   if ( NotVert ) then goto, output  

   ; create copy of NLAYERS to pass to getsigma
   NLayers = MTYPE.nlayers

   ; keyword psurf overrides grid definition structure
   if ( N_Elements( PSURF ) ne 1 ) then PSURF = MTYPE.psurf
   
   ; For MOPITT grid only: Always set PSURF = 1000 hPa, so that the 
   ; pressure levels will be computed correctly (clh, bmy, 10/15/02)
   if ( StrTrim( MType.Family, 2 ) eq 'MOPITT' ) then Psurf = 1000.0

   ; Get sigma/eta and pressure coordinates
   if ( MType.Hybrid ) then begin

      ;=================================================================
      ; HYBRID GRIDS (e.g. GEOS-4, GEOS-5, GISS_II_PRIME 23L) 
      ; We must use the ETA coordinate for vertical levels
      ;=================================================================

      if ( StrPos( MType.Name, 'GEOS5' ) ge 0 ) then begin
      
         ;--------------------------------------------------------------
         ; Special handling for GEOS-5: Read pressure edges (EPRESS)
         ; and centers (CPRESS) from disk.  Then compute the eta edges
         ; (EETA) and centers (CETA) accordingly. (bmy, 5/22/07)
         ;--------------------------------------------------------------
      
         ; Get GEOS-5 pressure edges & centers 
         Get_Geos5_Press, EPRESS, CPRESS,  $
            NLayers=NLayers, Res=MTYPE.Resolution, PSurf=PSurf
            
         ; Compute ETA on box edges & centers from
         EETA = ( EPRESS - MTYPE.Ptop ) / ( PSurf - MTYPE.Ptop )
         CETA = ( CPRESS - MTYPE.Ptop ) / ( PSurf - MTYPE.Ptop )

         ; GEOS-5 is a hybrid grid, set sigma levels to zero
         ESIG = FltArr( N_Elements( EPRESS ) )
         CSIG = FltArr( N_Elements( CPRESS ) )

      endif else begin

         ;--------------------------------------------------------------
         ; Other hybrid grids (e.g. GEOS-4): First get the ETA edges
         ; (EETA) and centers (CETA).  Then use those, along with the
         ; surface pressure PSURF to compute the pressure edges 
         ; (EPRESS) and centers (CPRESS).
         ;--------------------------------------------------------------

         ; Get ETA coordinate at grid box centers
         CETA   = GetEta( MType.Name, NLayers, PSurf=PSurf, /Center )

         ; check if number of layers is as expected
         if ( NLayers ne MTYPE.nlayers ) $
            then Message, '** number of layers inconsistent ! **'

         ; Get ETA coordinates at grid box edges
         EETA   = GetEta( MTYPE.name, nlayers, PSurf=PSurf, /Edges )

         ; Construct pressures at box centers and edges
         ; Pressure = ( ETA * ( Psurf - PTOP ) ) + PTOP
         CPRESS = CETA * ( psurf - MTYPE.ptop ) + MTYPE.ptop
         EPRESS = EETA * ( psurf - MTYPE.ptop ) + MTYPE.ptop
      
         ; Set ESIG, CSIG to zero, we don't have sigma coordinates
         ESIG   = FltArr( N_Elements( EETA ) )
         CSIG   = FltArr( N_Elements( CETA ) )

      endelse

   endif else begin

      ;=================================================================
      ; PURE SIGMA GRIDS (e.g. GEOS-1, GEOS-STRAT, GEOS-3)
      ; We must use the SIGMA coordinate for vertical levels
      ;=================================================================      

      ; For multi-level generic grids, skip to output (bmy, 6/28/06)
      if ( MType.Name eq 'GENERIC' ) then begin
         CSIG   = 0.0
         ESIG   = 0.0
         CPRESS = 0.0
         EPRESS = 0.0
         EETA   = 0.0
         CETA   = 0.0
         mZEDGE = 0.0
         mZMID  = 0.0
         goto, output
      endif

      ; Get SIGMA centers
      CSIG = GetSigma( MTYPE.name, NLayers, /CENTER )

      ; check if number of layers is as expected
      if ( nlayers ne MTYPE.nlayers ) $
         then Message, '** number of layers inconsistent ! **'
   
      ; Get SIGMA edges
      ESIG = GetSigma( MTYPE.name, NLayers, /EDGES )
   
      ; compute approx. pressure levels
      ;     sigma = (P - Ptop) / (Psurf - Ptop)
   
      ; Pressure levels for pure sigma-level models
      ; This includes GEOS-1, GEOS-STRAT, GEOS-3, GISS_II, FSU, etc
      CPRESS = CSIG * ( psurf - MTYPE.ptop ) + MTYPE.ptop
      EPRESS = ESIG * ( psurf - MTYPE.ptop ) + MTYPE.ptop

      ; Set EETA and CETA to zero for pure sigma models
      EETA = FltArr( N_Elements( ESIG ) )
      CETA = FltArr( N_Elements( CSIG ) )

   endelse

   ; Call function USSA_ALT to convert Edge and center 
   ; pressures to altitudes (bmy, 6/19/98) 
   mZEDGE = USSA_Alt( EPRESS )
   mZMID  = USSA_Alt( CPRESS )

   ;====================================================================
   ; handle output request : if an output keyword is set, return the 
   ; respective value, else create a structure with all the information
   ;====================================================================

output:
   ; only alternative name for IMX
   if ( keyword_set( NLON    ) ) then return, mIMX 
   if ( keyword_set( NLAT    ) ) then return, mJMX
   if ( keyword_set( IMX     ) ) then return, mIMX
   if ( keyword_set( JMX     ) ) then return, mJMX
   if ( keyword_set( DI      ) ) then return, mDI
   if ( keyword_set( DJ      ) ) then return, mDJ
   if ( keyword_set( XEDGE   ) ) then return, mXEDGE
   if ( keyword_set( YEDGE   ) ) then return, mYEDGE
   if ( keyword_set( XMID    ) ) then return, mXMID
   if ( keyword_set( YMID    ) ) then return, mYMID
   
   if ( keyword_set( LMX     ) ) then return, MTYPE.nlayers
   if ( keyword_set( SIGEDGE ) ) then return, ESIG
   if ( keyword_set( SIGMID  ) ) then return, CSIG
   if ( keyword_set( ETAEDGE ) ) then return, EETA
   if ( keyword_set( ETAMID  ) ) then return, CETA
   if ( keyword_set( PEDGE   ) ) then return, EPRESS
   if ( keyword_set( PMID    ) ) then return, CPRESS
   if ( keyword_set( ZEDGE   ) ) then return, mZEDGE
   if ( keyword_set( ZMID    ) ) then return, mZMID

   result = { IMX:mIMX,     JMX:mJMX,     DI:mDI,     DJ:mDJ,     $
              XEDGE:mXEDGE, YEDGE:mYEDGE, XMID:mXMID, YMID:mYMID }

   ; Added ETAMID, ETAEDGE tags (bmy, 11/6/01)
   if (not notvert) then $
      result = CREATE_STRUCT(result,                   $
                             'LMX',     MTYPE.nlayers, $
                             'SIGMID',  CSIG,          $
                             'SIGEDGE', ESIG,          $
                             'ETAMID',  CETA,          $
                             'ETAEDGE', EETA,          $
                             'PMID',    CPRESS,        $
                             'PEDGE',   EPRESS,        $
                             'ZEDGE',   mZEDGE,        $
                             'ZMID',    mZMID )

   return,result


end



