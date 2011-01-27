#include <misc.h>
#include <preproc.h>

module Hydrology2Mod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: Hydrology2Mod
! 
! !DESCRIPTION: 
! Calculation of soil/snow hydrology
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: Hydrology2        ! Calcultes soil/snow hydrology
!
! !REVISION HISTORY:
! 2/28/02 Peter Thornton: Migrated to new data structures.
!
!EOP
!----------------------------------------------------------------------- 

  
contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Hydrology2
!
! !INTERFACE:
  subroutine Hydrology2 (c)
!
! !DESCRIPTION: 
! This is the main subroutine to execute the calculation of soil/snow 
! hydrology
! Calling sequence is: 
!  Hydrology2:                       surface hydrology driver
!    -> SnowWater:                   change of snow mass and snow water
!                                    onto soil
!    -> SurfaceRunoff:               surface runoff 
!    -> Infiltration:                infiltration into surface soil layer
!    -> SoilWater:                   soil water movement between layers
!          -> Tridiagonal            tridiagonal matrix solution
!    -> Drainage:                    subsurface runoff  
!    -> SnowCompaction:              compaction of snow layers
!    -> CombineSnowLayers:           combine snow layers that are thinner
!                                     than minimum
!    -> DivideSnowLayers:            subdivide snow layers that are thicker
!                                     than maximum
!    -> WetIceHydrology:             calculate hydrology for wetland and land
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use clm_varcon, only : denh2o, denice, istice, istwet, istsoil, spval
    use clm_varpar, only : nlevsoi, nlevsno
    use SnowHydrologyMod, only : SnowCompaction, CombineSnowLayers, DivideSnowLayers, SnowWater
    use SoilHydrologyMod, only : Infiltration, SoilWater, Drainage, SurfaceRunoff
    use WetIceHydrologyMod, only : WetIceHydrology
!
! !ARGUMENTS:
    implicit none
    type (column_type), target, intent(inout) :: c !column derived type
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
!
    integer , pointer :: itypwat !water type
    integer , pointer :: snl     !number of snow layers
    real(r8), pointer :: h2ocan  !canopy water (mm H2O)
    real(r8), pointer :: h2osno  !snow water (mm H2O)
!
! local pointers to implicit out scalars
!
    real(r8), pointer :: endwb   !water mass end of the time step
    real(r8), pointer :: snowage !non dimensional snow age [-]
    real(r8), pointer :: wf      !soil water as frac. of whc for top 0.5 m
    real(r8), pointer :: t_snow  !vertically averaged snow temperature
    real(r8), pointer :: t_grnd  !ground temperature (Kelvin)
    real(r8), pointer :: snowice !average snow ice lens
    real(r8), pointer :: snowliq !average snow liquid water
!
! local pointers to implicit in arrays
!
    real(r8), dimension(:), pointer :: watsat !volumetric soil water at saturation (porosity)
    real(r8), dimension(:), pointer :: sucsat !minimum soil suction (mm)
    real(r8), dimension(:), pointer :: bsw    !Clapp and Hornberger "b"
    real(r8), dimension(:), pointer :: z      !layer depth  (m)
!
! local pointers to implicit inout arrays
!
    real(r8), dimension(:), pointer :: dz     !layer thickness depth (m)
    real(r8), dimension(:), pointer :: zi     !interface depth (m)
!
! local pointers to implicit out arrays
!
    real(r8), dimension(:), pointer :: t_soisno     !soil temperature (Kelvin)
    real(r8), dimension(:), pointer :: h2osoi_ice   !ice lens (kg/m2)
    real(r8), dimension(:), pointer :: h2osoi_liq   !liquid water (kg/m2)
    real(r8), dimension(:), pointer :: h2osoi_vol   !volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: j                  !do loop index
    real(r8) :: zwice              !the sum of ice mass of soil (kg/m2)
    real(r8) :: vol_liq(1:nlevsoi) !partial volume of liquid water in layer
    real(r8) :: s(1:nlevsoi)       !wetness of soil (including ice)
    real(r8) :: zwt                !water table depth
    real(r8) :: fcov               !fractional area with water table at surface
    real(r8) :: dwat(1:nlevsoi)    !change in soil water
    real(r8) :: hk(1:nlevsoi)      !hydraulic conductivity (mm h2o/s)
    real(r8) :: dhkdw(1:nlevsoi)   !d(hk)/d(vol_liq)
    type(column_pstate_type), pointer ::  cps   !local pointers to derived subtypes
    type(column_estate_type), pointer ::  ces   !local pointers to derived subtypes
    type(column_wstate_type), pointer ::  cws   !local pointers to derived subtypes
    type(water_balance_type), pointer ::  cwbal !local pointers to derived subtypes 
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes
    cps => c%cps
    ces => c%ces
    cws => c%cws
    cwbal => c%cwbal

    ! Assign local pointers to derived type scalar members
    itypwat => cps%lps%itypwat
    endwb => cwbal%endwb
    snl => cps%snl
    snowage => cps%snowage
    t_snow => ces%t_snow
    t_grnd => ces%t_grnd
    h2ocan => cws%pws_a%h2ocan
    h2osno => cws%h2osno
    wf => cps%wf
    snowice => cws%snowice
    snowliq => cws%snowliq

    ! Assign local pointers to derived type array members
    watsat => cps%lps%watsat
    sucsat => cps%lps%sucsat
    bsw => cps%lps%bsw
    z => cps%z
    dz => cps%dz
    zi => cps%zi
    t_soisno => ces%t_soisno
    h2osoi_ice => cws%h2osoi_ice
    h2osoi_liq => cws%h2osoi_liq
    h2osoi_vol => cws%h2osoi_vol

    ! Determine the change of snow mass and the snow water onto soil
    call SnowWater (c)

    ! Determine soil hydrology
    if (itypwat == istsoil) then
       call SurfaceRunoff (c, zwice, vol_liq, s, zwt, fcov)
       call Infiltration (c)
       call SoilWater (c, vol_liq, dwat, hk, dhkdw)
       call Drainage (c, zwice, vol_liq, s, zwt, &
            fcov, hk, dhkdw, dwat)
    endif

    if (snl < 0) then
       ! Natural compaction and metamorphosis.
       call SnowCompaction (c)

       ! Combine thin snow elements
       call CombineSnowLayers (c)

       ! Divide thick snow elements
       call DivideSnowLayers (c)

       ! Set empty snow layers to zero
       if (snl > -nlevsno) then
          snowage = 0.
          do j = -nlevsno+1, snl
             h2osoi_ice(j) = 0.
             h2osoi_liq(j) = 0.
             t_soisno(j) = 0.
             dz(j) = 0.
             z(j) = 0.
             zi(j-1) = 0.
          enddo
       endif

    endif

    ! Vertically average t_soisno and sum of h2osoi_liq and h2osoi_ice 
    ! over all snow layers for the given patch - these will be written out
    ! to the history tape
    if (snl < 0) then
       t_snow  = 0.
       snowice = 0.
       snowliq = 0.
       do j = snl+1, 0
          t_snow  = t_snow  + t_soisno(j)
          snowice = snowice + h2osoi_ice(j)
          snowliq = snowliq + h2osoi_liq(j)
       end do
       t_snow = t_snow/abs(snl)
    else
       t_snow  = spval
       snowice = spval
       snowliq = spval
    endif

    ! Update ground temperature
    t_grnd = t_soisno(snl+1)

    ! Determine volumetric soil water
    do j = 1,nlevsoi
       h2osoi_vol(j) = h2osoi_liq(j)/(dz(j)*denh2o) &
            + h2osoi_ice(j)/(dz(j)*denice)
    end do

    ! Determine ending water balance
    endwb=h2ocan+h2osno
    do j = 1, nlevsoi
       endwb = endwb + h2osoi_ice(j) + h2osoi_liq(j)
    enddo

    ! Determine wetland and land ice hydrology (must be placed here 
    ! since need snow updated from CombineSnowLayers)
    if (itypwat==istwet .or. itypwat==istice) call WetIceHydrology (c)

  end subroutine Hydrology2

end module Hydrology2Mod



