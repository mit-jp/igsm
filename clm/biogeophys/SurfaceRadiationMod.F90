#include <misc.h>
#include <preproc.h>

module SurfaceRadiationMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: SurfaceRadiationMod
! 
! !DESCRIPTION: 
! Calculate solar fluxes absorbed by vegetation and ground surface
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: SurfaceRadiation ! Solar fluxes absorbed by veg and ground surface
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!----------------------------------------------------------------------- 

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SurfaceRadiation
!
! !INTERFACE:
  subroutine SurfaceRadiation (c)
!
! !DESCRIPTION: 
! Solar fluxes absorbed by vegetation and ground surface
! Note possible problem when land is on different grid than atmosphere.
! Land may have sun above the horizon (coszen > 0) but atmosphere may
! have sun below the horizon (forc_solad = 0 and forc_solai = 0). This is okay
! because all fluxes (absorbed, reflected, transmitted) are multiplied
! by the incoming flux and all will equal zero.
! Atmosphere may have sun above horizon (forc_solad > 0 and forc_solai > 0) but
! land may have sun below horizon. This is okay because fabd, fabi,
! ftdd, ftid, and ftii all equal zero so that sabv=sabg=fsa=0. Also,
! albd and albi equal one so that fsr=forc_solad+forc_solai. In other words, all
! the radiation is reflected. NDVI should equal zero in this case.
! However, the way the code is currently implemented this is only true
! if (forc_solad+forc_solai)|vis = (forc_solad+forc_solai)|nir.
! Output variables are parsun,parsha,sabv,sabg,fsa,fsr,ndvi
!
! !USES:
    use clmtype
!
! !ARGUMENTS:
    implicit none
    type (column_type),target,intent(inout):: c		!column derived type
!
! !CALLED FROM:
! subroutine Biogeophysics1 in module Biogeophysics1Mod
! subroutine BiogeophysicsLake in module BiogeophysicsLakeMod
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! 2/18/02, Peter Thornton: Migrated to new data structures. Added a pft loop.
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in scalars
!
    real(r8),pointer:: fsun   !sunlit fraction of canopy
    real(r8),pointer:: elai   !one-sided leaf area index with burying by snow
    real(r8),pointer:: esai   !one-sided stem area index with burying by snow
!
! local pointers to original implicit out scalars
!
    real(r8),pointer:: laisun !sunlit leaf area
    real(r8),pointer:: laisha !shaded leaf area
    real(r8),pointer:: sabg   !solar radiation absorbed by ground (W/m**2)
    real(r8),pointer:: sabv   !solar radiation absorbed by vegetation (W/m**2)
    real(r8),pointer:: fsa    !solar radiation absorbed (total) (W/m**2)
    real(r8),pointer:: parsun !average absorbed PAR for sunlit leaves (W/m**2)
    real(r8),pointer:: parsha !average absorbed PAR for shaded leaves (W/m**2)
    real(r8),pointer:: fsr    !solar radiation reflected (W/m**2)
    real(r8),pointer:: ndvi   !Normalized Difference Vegetation Index
!
! local pointers to original implicit in arrays
!
    real(r8),dimension(:),pointer:: forc_solad   !direct beam radiation (W/m**2)
    real(r8),dimension(:),pointer:: forc_solai   !diffuse radiation (W/m**2)
    real(r8),dimension(:),pointer:: fabd         !flux absorbed by veg per unit direct flux
    real(r8),dimension(:),pointer:: fabi         !flux absorbed by veg per unit diffuse flux
    real(r8),dimension(:),pointer:: ftdd         !down direct flux below veg per unit dir flx
    real(r8),dimension(:),pointer:: ftid         !down diffuse flux below veg per unit dir flx
    real(r8),dimension(:),pointer:: ftii         !down diffuse flux below veg per unit dif flx
    real(r8),dimension(:),pointer:: albgrd       !ground albedo (direct)
    real(r8),dimension(:),pointer:: albgri       !ground albedo (diffuse)
    real(r8),dimension(:),pointer:: albd         !surface albedo (direct)
    real(r8),dimension(:),pointer:: albi         !surface albedo (diffuse)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer:: ib              ! waveband number (1=vis, 2=nir)
    integer::  nband          ! number of solar radiation waveband classes
    real(r8):: abs            ! absorbed solar radiation (W/m**2)
    real(r8):: rnir           ! reflected solar radiation [nir] (W/m**2)
    real(r8):: rvis           ! reflected solar radiation [vis] (W/m**2)
    real(r8):: laifra         ! leaf area fraction of canopy
    real(r8):: trd            ! transmitted solar radiation: direct (W/m**2)
    real(r8):: tri            ! transmitted solar radiation: diffuse (W/m**2)
    real(r8):: cad(numrad)    ! direct beam absorbed by canopy (W/m**2)
    real(r8):: cai(numrad)    ! diffuse radiation absorbed by canopy (W/m**2)
    real(r8):: fsha           ! shaded fraction of canopy
    real(r8):: vai            ! total leaf area index + stem area index, one sided
    real(r8):: mpe            ! prevents overflow for division by zero
    integer :: pi             ! pft index

    ! local pointers to derived subtypes
    type(atm2lnd_flux_type),pointer::   a2lf
    type(column_pstate_type),pointer::  cps
    type(pft_type),pointer::            p
    type(pft_pstate_type),pointer::     pps
    type(pft_eflux_type),pointer::      pef
!-----------------------------------------------------------------------

    ! assign local pointers to derived subtypes (column-level)
    cps => c%cps
    a2lf => c%a2lf

    ! assign local pointers to derived type array members (column level)
    forc_solad => a2lf%forc_solad
    forc_solai => a2lf%forc_solai
    albgrd => cps%albgrd
    albgri => cps%albgri

    ! Begin pft loop
    mpe   = 1.e-06
    nband = numrad

    do pi=1,cps%npfts
       ! Assign local pointers to derived subtypes (pft-level)
       p => c%p(pi)
       pps => p%pps
       pef => p%pef

       ! Apssign local pointers to derived type scalar members (pft-level)
       fsun => pps%fsun
       elai => pps%elai
       esai => pps%esai
       laisun => pps%laisun
       laisha => pps%laisha
       ndvi => pps%ndvi
       sabg => pef%sabg
       sabv => pef%sabv
       fsa => pef%fsa
       fsr => pef%fsr
       parsun => pef%parsun
       parsha => pef%parsha

       ! Assign local pointers to derived type array members (pft level)
       fabd => pps%fabd
       fabi => pps%fabi
       ftdd => pps%ftdd
       ftid => pps%ftid
       ftii => pps%ftii
       albd => pps%albd
       albi => pps%albi

       ! Begin pft-level calculations
       fsha = 1.-fsun
       laisun = elai*fsun
       laisha = elai*fsha
       vai = elai+ esai

       ! Zero summed solar fluxes
       sabg = 0.
       sabv = 0.
       fsa  = 0.

       ! Loop over nband wavebands
       do ib = 1, nband
          ! Absorbed by canopy
          cad(ib)  = forc_solad(ib)*fabd(ib)
          cai(ib)  = forc_solai(ib)*fabi(ib)
          sabv = sabv + cad(ib) + cai(ib)
          fsa  = fsa  + cad(ib) + cai(ib)

          ! Transmitted = solar fluxes incident on ground
          trd = forc_solad(ib)*ftdd(ib)
          tri = forc_solad(ib)*ftid(ib) + forc_solai(ib)*ftii(ib)

          ! Solar radiation absorbed by ground surface
          abs = trd*(1.-albgrd(ib)) + tri*(1.-albgri(ib))
          sabg = sabg + abs
          fsa  = fsa  + abs

       end do

       ! Partition visible canopy absorption to sunlit and shaded fractions
       ! to get average absorbed par for sunlit and shaded leaves
       laifra = elai / max(vai,mpe)
       if (fsun > 0.) then
          parsun = (cad(1) + cai(1)) * laifra
          parsha = 0._r8
       else
          parsun = 0._r8
          parsha = 0._r8
       endif

       ! NDVI and reflected solar radiation
       rvis = albd(1)*forc_solad(1) + albi(1)*forc_solai(1)
       rnir = albd(2)*forc_solad(2) + albi(2)*forc_solai(2)
       fsr = rvis + rnir
       ndvi = (rnir-rvis) / max(rnir+rvis,mpe)

    end do ! (end of pfts loop)

    return
  end subroutine SurfaceRadiation

end module SurfaceRadiationMod
