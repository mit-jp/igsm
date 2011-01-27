#include <misc.h>
#include <preproc.h>

module clm_varpar

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: clm_varpar
!
! !DESCRIPTION: 
! Module containing clm2 parameters
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! Define land surface 2-d grid. This sets the model resolution according
! to cpp directives LSMLON and LSMLAT in preproc.h. 
!
  integer, parameter :: lsmlon = LSMLON  !maximum number of longitude points on lsm grid
  integer, parameter :: lsmlat = LSMLAT  !number of latitude points on lsm grid
!
! Define maximum number of PFT patches per grid cell and set
! patch number for urban, lake, wetland, and glacier patches
!
  integer, parameter :: maxpatch_pft = 16                !maximum number of PFT subgrid patches per grid cell
  integer, parameter :: npatch_urban = maxpatch_pft + 1 !urban   patch number: 1 to maxpatch
  integer, parameter :: npatch_lake  = npatch_urban + 1 !lake    patch number: 1 to maxpatch
  integer, parameter :: npatch_wet   = npatch_lake  + 1 !wetland patch number: 1 to maxpatch
  integer, parameter :: npatch_gla   = npatch_wet   + 1 !glacier patch number: 1 to maxpatch
  integer, parameter :: maxpatch     = npatch_gla       !maximum number of subgrid patches per grid cell
!
! Define number of levels
!
  integer, parameter :: nlevsoi     =  10   !number of soil layers
  integer, parameter :: nlevlak     =  10   !number of lake layers
  integer, parameter :: nlevsno     =   5   !maximum number of snow layers
!
! Define miscellaneous parameters
!
  integer, parameter :: numwat      =   5   !number of water types (soil, ice, 2 lakes, wetland)
  integer, parameter :: numpft      =  16   !number of plant types
  integer, parameter :: npftpar     =  32   !number of pft parameters (in LPJ) (DGVM only)
  integer, parameter :: numcol      =   8   !number of soil color types
  integer, parameter :: numrad      =   2   !number of solar radiation bands: vis, nir
  integer, parameter :: ndst        =   4   !number of dust size classes (BGC only)
  integer, parameter :: dst_src_nbr =   3   !number of size distns in src soil (BGC only)
  integer, parameter :: sz_nbr      = 200   !number of sub-grid bins in large bin of dust size distribution (BGC only)
  integer, parameter :: nvoc        =   5   !number of voc categories (BGC only)
!
! Define parameters for RTM river routing model
!
  integer, parameter :: rtmlon = 720  !number of rtm longitudes
  integer, parameter :: rtmlat = 360  !number of rtm latitudes
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!----------------------------------------------------------------------- 

end module clm_varpar

