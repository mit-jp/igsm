#include <misc.h>
#include <preproc.h>

module SoilHydrologyMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: SoilHydrologyMod
! 
! !DESCRIPTION: 
! Calculate soil hydrology
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: SurfaceRunoff  ! Calculate surface runoff
  public :: Infiltration   ! Calculate infiltration into surface soil layer
  public :: SoilWater      ! Calculate soil hydrology
  public :: Drainage       ! Calculate subsurface drainage
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
! !IROUTINE: SurfaceRunoff
!
! !INTERFACE:
  subroutine SurfaceRunoff (c, zwice, vol_liq, s, zwt, &
                            fcov ) 
!
! !DESCRIPTION: 
! Calculate surface runoff
! The original code was provide by Robert E. Dickinson based on 
! following clues:  exponential decrease of Ksat, a water table 
! level determination level including highland and lowland levels 
! and fractional area of wetland (water table above the surface). 
! Runoff is parameterized from the lowlands in terms of precip 
! incident on wet areas and a base flow, where these are estimated 
! using ideas from TOPMODEL.
! The original scheme was modified by Z.-L. Yang and G.-Y. Niu,
! o  using a new method to determine water table depth and
!    the fractional wet area (fcov)
! o  computing runoff (surface and subsurface) from this
!    fraction and the remaining fraction (i.e. 1-fcov)
! o  for the 1-fcov part, using BATS1e method to compute
!    surface and subsurface runoff.
! The original code on soil moisture and runoff were provided by 
! R. E. Dickinson in July 1996.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use globals
    use clm_varcon, only : denice, denh2o, wimp
    use clm_varpar, only : nlevsoi
!
! !ARGUMENTS:
    implicit none
    type (column_type),target,intent(inout):: c		!column derived type
    real(r8), intent(out) :: zwice              ! the sum of ice mass of soil (kg/m2)
    real(r8), intent(out) :: vol_liq(1:nlevsoi) ! partial volume of liquid water in layer
    real(r8), intent(out) :: s(1:nlevsoi)       ! wetness of soil (including ice)
    real(r8), intent(out) :: zwt                ! water table depth
    real(r8), intent(out) :: fcov               ! fractional area with water table at surface
!
! !CALLED FROM:
! subroutine BiogeophysicsLake in module BiogeophysicsLakeMod 
! subroutine Hydrology2 in module Hydrology2Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 12 November 1999:  Z.-L. Yang and G.-Y. Niu
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! 2/26/02, Peter Thornton: Migrated to new data structures.
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
    real(r8),pointer:: wtfact        !Fraction of model area with high water table
    real(r8),pointer:: qflx_top_soil !net water input into soil from top (mm/s)
!
! local pointers to original implicit out arguments
!
    real(r8),pointer:: qflx_surf     !surface runoff (mm H2O /s)
!
! local pointers to original implicit in arrays
!
    real(r8),dimension(:),pointer:: watsat !volumetric soil water at saturation (porosity)
    real(r8),dimension(:),pointer:: dz     !layer depth (m)
    real(r8),dimension(:),pointer:: zi     !interface level below a "z" level (m)
    real(r8),dimension(:),pointer:: h2osoi_ice   !ice lens (kg/m2) 
    real(r8),dimension(:),pointer:: h2osoi_liq   !liquid water (kg/m2)
!
! local pointers to original implicit out arrays
!
    real(r8),dimension(:),pointer:: eff_porosity !effective porosity = porosity - vol_ice
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer:: i                   ! do loop index
    real(r8):: vol_ice(1:nlevsoi) ! partial volume of ice lens in layer
    real(r8):: zmean              ! The surface soil layers contributing to runoff
    real(r8):: wmean              ! The averaged soil wetness in surface soil layers
    real(r8):: fz                 ! coefficient for water table depth
    type(column_pstate_type),pointer::  cps !local pointers to derived subtypes
    type(column_wstate_type),pointer::  cws !local pointers to derived subtypes
    type(column_wflux_type),pointer::   cwf !local pointers to derived subtypes
!-----------------------------------------------------------------------

    ! assign local pointers to derived subtypes
    cps => c%cps
    cws => c%cws
    cwf => c%cwf

    ! assign local pointers to derived type scalar members
    wtfact => cps%lps%gps%wtfact
    qflx_top_soil => cwf%qflx_top_soil
    qflx_surf => cwf%qflx_surf

    ! assign local pointers to derived type array members
    watsat => cps%lps%watsat
    dz => cps%dz
    zi => cps%zi
    h2osoi_ice => cws%h2osoi_ice
    h2osoi_liq => cws%h2osoi_liq
    eff_porosity => cps%eff_porosity

    ! Porosity of soil, partial volume of ice and liquid
    zwice = 0.
    do i = 1,nlevsoi
       zwice = zwice + h2osoi_ice(i)
       vol_ice(i) = min(watsat(i), h2osoi_ice(i)/(dz(i)*denice))
       eff_porosity(i) = watsat(i)-vol_ice(i)
       vol_liq(i) = min(eff_porosity(i), h2osoi_liq(i)/(dz(i)*denh2o))
    enddo

    ! Calculate wetness of soil
    do i = 1,nlevsoi
       s(i) = min(1._r8,(vol_ice(i)+vol_liq(i))/watsat(i))
    end do

    ! Determine water table 
    wmean = 0.
    fz    = 1.0
    do i  = 1, nlevsoi
       wmean = wmean + s(i)*dz(i)
    enddo
    zwt = fz * (zi(nlevsoi) - wmean)

    ! Saturation fraction
    fcov = wtfact*min(1._r8,exp(-zwt))

    ! Currently no overland flow parameterization in code is considered
    ! qflx_surf = 0.   Zong-Liang Yang & G.-Y. Niu                         
    !*modified surface runoff according to the concept of TOPMODEL   
    wmean = 0.
    zmean = 0.
    do i = 1, 3
       zmean = zmean + dz(i)
       wmean = wmean + s(i) * dz(i)
    enddo
    wmean = wmean / zmean

    ! If top soil layer is impermeable then all qflx_top_soil goes to surface runoff
    if (eff_porosity(1) < wimp) then
       qflx_surf =  max(0._r8, fcov*qflx_top_soil) + &
            max(0._r8, (1.-fcov)*qflx_top_soil)
    else
       qflx_surf =  max(0._r8, fcov*qflx_top_soil) + &
            max(0._r8, (1.-fcov)*min(1._r8,wmean**4)*qflx_top_soil)
    endif

  end subroutine SurfaceRunoff

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Infiltration
!
! !INTERFACE:
  subroutine Infiltration (c) 
!
! !DESCRIPTION: 
! Calculate infiltration into surface soil layer (minus the evaporation)
! The original code on soil moisture and runoff were provided by 
! R. E. Dickinson in July 1996.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
!
! !ARGUMENTS:
    implicit none
    type (column_type),target,intent(inout):: c		!column derived type
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 12 November 1999:  Z.-L. Yang and G.-Y. Niu
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! 2/27/02, Peter Thornton: Migrated to new data structures.
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
    integer,pointer:: snl            !number of snow layers
    real(r8),pointer:: qflx_top_soil !net water input into soil from top (mm/s)
    real(r8),pointer:: qflx_surf     !surface runoff (mm H2O /s)
    real(r8),pointer:: qflx_evap_grnd!ground surface evaporation rate (mm H2O/s) [+]
!
! local pointers to original implicit out arguments
!
    real(r8),pointer:: qflx_infl     !infiltration (mm H2O /s)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    type(column_pstate_type),pointer::  cps !local pointer to derived subtypes
    type(column_wflux_type),pointer::   cwf !local pointer to derived subtypes
!-----------------------------------------------------------------------

    ! assign local pointers to derived subtypes
    cps => c%cps
    cwf => c%cwf

    ! assign local pointers to derived type scalar members
    snl => cps%snl
    qflx_top_soil => cwf%qflx_top_soil
    qflx_surf => cwf%qflx_surf
    qflx_evap_grnd => cwf%pwf_a%qflx_evap_grnd
    qflx_infl => cwf%qflx_infl

    ! Infiltration into surface soil layer (minus the evaporation)
    if (snl+1 >= 1) then
       qflx_infl = qflx_top_soil - qflx_surf - qflx_evap_grnd
    else
       qflx_infl = qflx_top_soil - qflx_surf
    endif

  end subroutine Infiltration

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SoilWater
!
! !INTERFACE:
  subroutine SoilWater (c, vol_liq, dwat, hk, dhkdw)
!
! !DESCRIPTION: 
! Soil hydrology
! Soil moisture is predicted from a 10-layer model (as with soil 
! temperature), in which the vertical soil moisture transport is governed
! by infiltration, runoff, gradient diffusion, gravity, and root 
! extraction through canopy transpiration.  The net water applied to the
! surface layer is the snowmelt plus precipitation plus the throughfall 
! of canopy dew minus surface runoff and evaporation. 
!
! The vertical water flow in an unsaturated porous media is described by
! Darcy's law, and the hydraulic conductivity and the soil negative 
! potential vary with soil water content and soil texture based on the work 
! of Clapp and Hornberger (1978) and Cosby et al. (1984). The equation is
! integrated over the layer thickness, in which the time rate of change in
! water mass must equal the net flow across the bounding interface, plus the
! rate of internal source or sink. The terms of water flow across the layer
! interfaces are linearly expanded by using first-order Taylor expansion.  
! The equations result in a tridiagonal system equation. 
!
! Note: length units here are all millimeter 
! (in temperature subroutine uses same soil layer 
! structure required but lengths are m)
!
! Richards equation:
!
! d wat      d     d wat d psi
! ----- = - -- [ k(----- ----- - 1) ] + S
!   dt      dz       dz  d wat
!
! where: wat = volume of water per volume of soil (mm**3/mm**3)
! psi = soil matrix potential (mm)
! dt  = time step (s)
! z   = depth (mm)
! dz  = thickness (mm)
! qin = inflow at top (mm h2o /s) 
! qout= outflow at bottom (mm h2o /s)
! s   = source/sink flux (mm h2o /s) 
! k   = hydraulic conductivity (mm h2o /s)
!
!                       d qin                  d qin
! qin[n+1] = qin[n] +  --------  d wat(j-1) + --------- d wat(j)
!                       d wat(j-1)             d wat(j)
!                ==================|================= 
!                                  < qin 
!
!                 d wat(j)/dt * dz = qin[n+1] - qout[n+1] + S(j) 
!
!                                  > qout
!                ==================|================= 
!                        d qout               d qout
! qout[n+1] = qout[n] + --------- d wat(j) + --------- d wat(j+1)
!                        d wat(j)             d wat(j+1)
!
!
! Solution: linearize k and psi about d wat and use tridiagonal 
! system of equations to solve for d wat, 
! where for layer j
!
!
! r_j = a_j [d wat_j-1] + b_j [d wat_j] + c_j [d wat_j+1]
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use globals, only: dtime
    use clm_varcon, only: wimp
    use clm_varpar, only : nlevsoi
    use shr_const_mod, only : SHR_CONST_TKFRZ, SHR_CONST_LATICE, SHR_CONST_G
    use TridiagonalMod, only : Tridiagonal
!
! !ARGUMENTS:
    implicit none
    type (column_type),target,intent(inout):: c		!column derived type
    real(r8), intent(in) :: vol_liq(1 : nlevsoi) ! soil water per unit volume [mm/mm]
    real(r8), intent(out) :: dwat (1 : nlevsoi)  ! change of soil water [m3/m3]
    real(r8), intent(out) :: hk   (1 : nlevsoi)  ! hydraulic conductivity [mm h2o/s]
    real(r8), intent(out) :: dhkdw(1 : nlevsoi)  ! d(hk)/d(vol_liq)
!
! !CALLED FROM:
! subroutine Hydrology2 in module Hydrology2Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! 2/27/02, Peter Thornton: Migrated to new data structures. Includes
! treatment of multiple PFTs on a single soil column.
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
    real(r8),pointer:: smpmin     !restriction for min of soil potential (mm)
    real(r8),pointer:: qflx_infl     !infiltration (mm H2O /s)
    real(r8),pointer:: qflx_tran_veg !vegetation transpiration (mm H2O/s) (+ = to atm)
!
! local pointers to original implicit in arrays
!
    real(r8),dimension(:),pointer:: eff_porosity!effective porosity = porosity - vol_ice
    real(r8),dimension(:),pointer:: watsat    !volumetric soil water at saturation (porosity)
    real(r8),dimension(:),pointer:: hksat     !hydraulic conductivity at saturation (mm H2O /s)
    real(r8),dimension(:),pointer:: bsw       !Clapp and Hornberger "b"
    real(r8),dimension(:),pointer:: sucsat    !minimum soil suction (mm)
    real(r8),dimension(:),pointer:: t_soisno  !soil temperature (Kelvin)
!
! local pointers to original implicit inout arrays
!
    real(r8),dimension(:),pointer:: h2osoi_liq   !liquid water (kg/m2)
!
! local pointers to original implicit out arrays
!
    real(r8),dimension(:),pointer:: rootr_column !effective fraction of roots in each soil layer
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer:: j                  ! do loop indices 
    real(r8):: amx(1:nlevsoi)    ! "a" left off diagonal of tridiagonal matrix
    real(r8):: bmx(1:nlevsoi)    ! "b" diagonal column for tridiagonal matrix
    real(r8):: cmx(1:nlevsoi)    ! "c" right off diagonal tridiagonal matrix
    real(r8):: z (1 : nlevsoi)   ! layer depth [mm]
    real(r8):: dz(1 : nlevsoi)   ! layer thickness [mm]
    real(r8):: den               ! used in calculating qin, qout
    real(r8):: dqidw0            ! d(qin)/d(vol_liq(i-1))
    real(r8):: dqidw1            ! d(qin)/d(vol_liq(i))
    real(r8):: dqodw1            ! d(qout)/d(vol_liq(i))
    real(r8):: dqodw2            ! d(qout)/d(vol_liq(i+1))
    real(r8):: dsmpdw(1:nlevsoi) ! d(smp)/d(vol_liq)
    real(r8):: num               ! used in calculating qin, qout
    real(r8):: qin               ! flux of water into soil layer [mm h2o/s]
    real(r8):: qout              ! flux of water out of soil layer [mm h2o/s]
    real(r8):: rmx(1:nlevsoi)    ! "r" forcing term of tridiagonal matrix
    real(r8):: s_node            ! soil wetness
    real(r8):: s1                ! "s" at interface of layer
    real(r8):: s2                ! k*s**(2b+2)
    real(r8):: smp(1:nlevsoi)    ! soil matrix potential [mm]
    real(r8):: sdamp             ! extrapolates soiwat dependence of evaporation
    integer :: pi                !pft index
    real(r8):: pft_wt            !pft weight for averaging
    real(r8):: pft_qflx_tran_veg !vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8):: temp1 
    real(r8):: temp2 ! accumulator for rootr weighting
    type(column_pstate_type),pointer::  cps  !local pointers to derived subtypes
    type(column_estate_type),pointer::  ces  !local pointers to derived subtypes
    type(column_wstate_type),pointer::  cws  !local pointers to derived subtypes
    type(column_wflux_type),pointer::   cwf  !local pointers to derived subtypes
    type(pft_type),pointer::            p    !local pointers to derived subtypes
    type(pft_pstate_type),pointer::     pps  !local pointers to derived subtypes 
    type(pft_wflux_type),pointer::      pwf  !local pointers to derived subtypes 
!-----------------------------------------------------------------------

    ! assign local pointers to derived subtypes
    cps => c%cps
    ces => c%ces
    cws => c%cws
    cwf => c%cwf

    ! assign local pointers to derived type scalar members
    smpmin => cps%lps%smpmin
    qflx_infl => cwf%qflx_infl
    qflx_tran_veg => cwf%pwf_a%qflx_tran_veg

    ! assign local pointers to derived type array members
    watsat => cps%lps%watsat
    hksat => cps%lps%hksat
    bsw => cps%lps%bsw
    sucsat => cps%lps%sucsat
    eff_porosity => cps%eff_porosity
    rootr_column => cps%rootr_column
    t_soisno => ces%t_soisno
    h2osoi_liq => cws%h2osoi_liq

    ! Because the depths in this routine are in mm, using local
    ! variable arrays instead of pointers
    do j = 1, nlevsoi
       z(j) = cps%z(j)*1.e3
       dz(j) = cps%dz(j)*1.e3
    enddo

    ! First step is to calculate the column-level effective rooting
    ! fraction in each soil layer. This is done outside the usual
    ! PFT-to-column averaging routines because it is not a simple
    ! weighted average of the PFT level rootr arrays. Instead, the
    ! weighting depends on both the per-unit-area transpiration 
    ! of the PFT and the PFTs area relative to all PFTs.
    do j = 1, nlevsoi
       rootr_column(j) = 0.
    end do

    temp2 = 0.
    do pi=1,cps%npfts
       p => c%p(pi)
       pps => p%pps
       pwf => p%pwf

       temp1 = pwf%qflx_tran_veg * pps%wt
       do j=1,nlevsoi
          rootr_column(j) = rootr_column(j) + pps%rootr(j) * temp1
       end do
       temp2 = temp2 + temp1
    end do

    ! normalize column-level rootr
    if (temp2 /= 0.) then
       do j=1,nlevsoi
          rootr_column(j) = rootr_column(j)/temp2
       end do
    endif

    ! set initial values
    sdamp = 0.

    ! Set hydraulic conductivity to zero if effective porosity 5% in any of 
    ! two neighbor layers or liquid content (theta) less than 0.001
    do j = 1, nlevsoi
       if (      (eff_porosity(j) < wimp) &
            .or. (eff_porosity(min(nlevsoi,j+1)) < wimp) &
            .or. (vol_liq(j) <= 1.e-3))then
          hk(j) = 0.
          dhkdw(j) = 0.
       else
          s1 = 0.5*(vol_liq(j)+vol_liq(min(nlevsoi,j+1))) / &
               (0.5*(watsat(j)+watsat(min(nlevsoi,j+1))))
          s2 = hksat(j)*s1**(2.*bsw(j)+2.)
          hk(j) = s1*s2
          dhkdw(j) = (2.*bsw(j)+3.)*s2*0.5/watsat(j)
          if(j == nlevsoi) dhkdw(j) = dhkdw(j) * 2.
       endif
    enddo

    ! Evaluate hydraulic conductivity, soil matric potential,
    ! d(smp)/d(vol_liq), and d(hk)/d(vol_liq).
    do j = 1, nlevsoi
       if (t_soisno(j)>SHR_CONST_TKFRZ) then

          s_node = max(vol_liq(j)/watsat(j),0.01_r8)
          s_node = min(1.0_r8,s_node)
          smp(j) = -sucsat(j)*s_node**(-bsw(j))
          smp(j) = max(smpmin, smp(j))        ! Limit soil suction
          dsmpdw(j) = -bsw(j)*smp(j)/(s_node*watsat(j))

       else
          ! When ice is present, the matric potential is only related to temperature
          ! by (Fuchs et al., 1978: Soil Sci. Soc. Amer. J. 42(3):379-385)
          ! Unit 1 Joule = 1 (kg m2/s2), J/kg /(m/s2) ==> m ==> 1e3 mm 
          smp(j) = 1.e3 * SHR_CONST_LATICE/SHR_CONST_G*(t_soisno(j)-SHR_CONST_TKFRZ)/t_soisno(j)
          smp(j) = max(smpmin, smp(j))        ! Limit soil suction
          dsmpdw(j) = 0.

       endif
    enddo

    ! Set up r, a, b, and c vectors for tridiagonal solution
    ! Node j=1
    j      = 1
    qin    = qflx_infl
    den    = (z(j+1)-z(j))
    num    = (smp(j+1)-smp(j)) - den
    qout   = -hk(j)*num/den
    dqodw1 = -(-hk(j)*dsmpdw(j)   + num*dhkdw(j))/den
    dqodw2 = -( hk(j)*dsmpdw(j+1) + num*dhkdw(j))/den
    rmx(j) =  qin - qout - qflx_tran_veg*rootr_column(j)
    amx(j) =  0.
    bmx(j) =  dz(j)*(sdamp+1./dtime) + dqodw1
    cmx(j) =  dqodw2

    ! Nodes j=2 to j=nlevsoi-1
    do j = 2, nlevsoi - 1
       den    = (z(j) - z(j-1))
       num    = (smp(j)-smp(j-1)) - den
       qin    = -hk(j-1)*num/den
       dqidw0 = -(-hk(j-1)*dsmpdw(j-1) + num*dhkdw(j-1))/den
       dqidw1 = -( hk(j-1)*dsmpdw(j)   + num*dhkdw(j-1))/den
       den    = (z(j+1)-z(j))
       num    = (smp(j+1)-smp(j)) - den
       qout   = -hk(j)*num/den
       dqodw1 = -(-hk(j)*dsmpdw(j)   + num*dhkdw(j))/den
       dqodw2 = -( hk(j)*dsmpdw(j+1) + num*dhkdw(j))/den
       rmx(j) =  qin - qout - qflx_tran_veg*rootr_column(j)
       amx(j) = -dqidw0
       bmx(j) =  dz(j)/dtime - dqidw1 + dqodw1
       cmx(j) =  dqodw2
    enddo

    ! Node j=nlevsoi
    j      = nlevsoi
    den    = (z(j) - z(j-1))
    num    = (smp(j)-smp(j-1)) - den
    qin    = -hk(j-1)*num/den
    dqidw0 = -(-hk(j-1)*dsmpdw(j-1) + num*dhkdw(j-1))/den
    dqidw1 = -( hk(j-1)*dsmpdw(j)   + num*dhkdw(j-1))/den
    qout   =  hk(j)
    dqodw1 =  dhkdw(j)
    rmx(j) =  qin - qout - qflx_tran_veg*rootr_column(j)
    amx(j) = -dqidw0
    bmx(j) =  dz(j)/dtime - dqidw1 + dqodw1
    cmx(j) =  0.

    ! Solve for dwat
    call Tridiagonal (nlevsoi, amx, bmx, cmx, rmx, &
                      dwat)

    ! Renew the mass of liquid water
    do j= 1,nlevsoi
       h2osoi_liq(j) = h2osoi_liq(j) + dwat(j)*dz(j)
    enddo

  end subroutine SoilWater

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Drainage
!
! !INTERFACE:
  subroutine Drainage (c,  zwice, vol_liq, s, zwt, &
       fcov, hk, dhkdw, dwat) 
!
! !DESCRIPTION: 
! Calculate subsurface drainage
! The original code was provide by Robert E. Dickinson based on 
! following clues:  exponential decrease of Ksat, a water table 
! level determination level including highland and lowland levels 
! and fractional area of wetland (water table above the surface). 
! Runoff is parameterized from the lowlands in terms of precip 
! incident on wet areas and a base flow, where these are estimated 
! using ideas from TOPMODEL.
! The original scheme was modified by Z.-L. Yang and G.-Y. Niu,
! *  using a new method to determine water table depth and
!    the fractional wet area (fcov)
! *  computing runoff (surface and subsurface) from this
!    fraction and the remaining fraction (i.e. 1-fcov)
! *  for the 1-fcov part, using BATS1e method to compute
!    surface and subsurface runoff.
! The original code on soil moisture and runoff were provided by 
! R. E. Dickinson in July 1996.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use globals, only: dtime
    use clm_varcon, only: pondmx
    use clm_varpar, only : nlevsoi
!
! !ARGUMENTS:
    implicit none
    type (column_type),target,intent(inout):: c		!column derived type
    real(r8), intent(in) :: vol_liq(1:nlevsoi) ! partial volume of liquid water in layer
    real(r8), intent(in) :: zwice              ! the sum of ice mass of soil (kg/m2)
    real(r8), intent(in) :: s(1:nlevsoi)       ! wetness of soil (including ice)
    real(r8), intent(in) :: zwt                ! water table depth
    real(r8), intent(in) :: fcov               ! fractional area with water table at surface
    real(r8), intent(in) :: hk(1:nlevsoi)      ! hydraulic conductivity (mm h2o/s)
    real(r8), intent(in) :: dhkdw(1:nlevsoi)   ! d(hk)/d(vol_liq)
    real(r8), intent(in) :: dwat(1:nlevsoi)    ! change in soil water
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 12 November 1999:  Z.-L. Yang and G.-Y. Niu
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in scalars
!
    integer,pointer:: snl            !number of snow layers
    real(r8),pointer:: qflx_snowcap  !excess precipitation due to snow capping (mm H2O /s) [+]
    real(r8),pointer:: qflx_dew_grnd !ground surface dew formation (mm H2O /s) [+]
    real(r8),pointer:: qflx_dew_snow !surface dew added to snow pack (mm H2O /s) [+]
    real(r8),pointer:: qflx_sub_snow !sublimation rate from snow pack (mm H2O /s) [+]
!
! local pointers to original implicit out scalars
!
    real(r8),pointer:: qflx_drain    !sub-surface runoff (mm H2O /s)
    real(r8),pointer:: qflx_qrgwl    !qflx_surf at glaciers, wetlands, lakes
    real(r8),pointer:: eflx_impsoil  !implicit evaporation for soil temperature equation
!
! local pointers to original implicit in arrays
!
    real(r8),dimension(:),pointer:: dz           !layer depth (m)
    real(r8),dimension(:),pointer:: bsw          !Clapp and Hornberger "b"
    real(r8),dimension(:),pointer:: eff_porosity !effective porosity = porosity - vol_ice
!
! local pointers to original implicit inout arrays
!
    real(r8),dimension(:),pointer:: h2osoi_ice   !ice lens (kg/m2)
    real(r8),dimension(:),pointer:: h2osoi_liq   !liquid water (kg/m2)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer i                ! do loop index
    real(r8) xs              ! excess soil water above saturation
    real(r8) dzmm(1:nlevsoi) ! layer thickness (mm)
    real(r8) watmin          ! minimum soil moisture
    real(r8) hksum           ! summation of hydraulic cond for layers 6->9
    real(r8) zsat            ! hydraulic conductivity weighted soil thickness
    real(r8) wsat            ! hydraulic conductivity weighted soil wetness
    real(r8) qflx_drain_wet  ! subsurface runoff from "wet" part (mm h2o/s)
    real(r8) qflx_drain_dry  ! subsurface runoff from "dry" part (mm h2o/s)
    real(r8) dzksum          ! hydraulic conductivity weighted soil thickness

    ! local pointers to derived subtypes
    type(column_pstate_type),pointer::  cps
    type(column_wstate_type),pointer::  cws
    type(column_wflux_type),pointer::   cwf
    type(column_eflux_type),pointer::   cef
!-----------------------------------------------------------------------

    ! assign local pointers to derived subtypes
    cps => c%cps
    cws => c%cws
    cwf => c%cwf
    cef => c%cef

    ! assign local pointers to derived type scalar members
    snl => cps%snl
    qflx_snowcap => cwf%pwf_a%qflx_snowcap
    qflx_dew_grnd => cwf%pwf_a%qflx_dew_grnd
    qflx_dew_snow => cwf%pwf_a%qflx_dew_snow
    qflx_sub_snow => cwf%pwf_a%qflx_sub_snow
    qflx_drain => cwf%qflx_drain
    qflx_qrgwl => cwf%qflx_qrgwl
    eflx_impsoil => cef%eflx_impsoil

    ! assign local pointers to derived type array members
    dz => cps%dz
    bsw => cps%lps%bsw
    eff_porosity => cps%eff_porosity
    h2osoi_liq => cws%h2osoi_liq
    h2osoi_ice => cws%h2osoi_ice

    ! begin calculation
    ! Streamflow and total runoff

    ! Convert layer thicknesses from m to mm
    do i = 1,nlevsoi
       dzmm(i) = dz(i)*1.e3
    enddo

    ! The amount of streamflow is assumed maintained by flow from the 
    ! lowland water table with different levels contributing according to 
    ! their thickness and saturated hydraulic conductivity, i.e. a given 
    ! layer below the water table interface loses water at a rate per unit 
    ! depth given by qflx_drain*hk/(sum over all layers below this water table 
    ! of hk*dz). Because this is a slow smooth process, and not strongly 
    ! coupled to water in any one layer, it should remain stable for 
    ! explicit time differencing. Hence, for simplicity it is removed
    ! explicitly prior to the main soil water calculation.
    ! Another assumption: no subsurface runoff for ice mixed soil 
    ! Zong-Liang Yang & G.-Y. Niu                                         
    qflx_drain = 0.                      ! subsurface runoff
    qflx_drain_wet = 0.                      ! subsurface runoff
    qflx_drain_dry = 0.                      ! subsurface runoff

    hksum = 0.
    do i = 6,nlevsoi-1
       hksum = hksum + hk(i)
    enddo
    if (zwice <= 0. .AND. hksum > 0.) then
       zsat = 0.
       wsat = 0.
       dzksum = 0.
       do i = 6,nlevsoi-1
          zsat = zsat + dz(i)*hk(i)
          wsat = wsat + s(i)*dz(i)*hk(i)
          dzksum  = dzksum   + hk(i)*dz(i)
       enddo
       wsat = wsat / zsat

       qflx_drain_dry = (1.-fcov)*4.e-2* wsat ** (2.*bsw(1)+3.)  ! mm/s
       qflx_drain_wet = fcov * 1.e-5 * exp(-zwt)                     ! mm/s
       qflx_drain = qflx_drain_dry + qflx_drain_wet

       do i = 6, nlevsoi-1
          h2osoi_liq(i) = h2osoi_liq(i) &
               - dtime*qflx_drain*dz(i)*hk(i)/dzksum
       enddo
    endif

    ! Limit h2osoi_liq to be greater than or equal to watmin. 
    ! Get water needed to bring h2osoi_liq equal watmin from lower layer. 
    watmin = 0.0
    do i = 1, nlevsoi-1
       if (h2osoi_liq(i) < 0.) then
          xs = watmin-h2osoi_liq(i)
       else
          xs = 0.
       end if
       h2osoi_liq(i  ) = h2osoi_liq(i  ) + xs
       h2osoi_liq(i+1) = h2osoi_liq(i+1) - xs
    end do
    i = nlevsoi
    if (h2osoi_liq(i) < watmin) then
       xs = watmin-h2osoi_liq(i)
    else
       xs = 0.
    end if
    h2osoi_liq(i) = h2osoi_liq(i) + xs
    qflx_drain = qflx_drain - xs/dtime

    ! Determine water in excess of saturation
    xs = max(0._r8, h2osoi_liq(1)-(pondmx+eff_porosity(1)*dzmm(1)))
    if (xs > 0.) h2osoi_liq(1) = pondmx+eff_porosity(1)*dzmm(1)

    do i = 2,nlevsoi
       xs = xs + max(h2osoi_liq(i)-eff_porosity(i)*dzmm(i), 0._r8)  ! [mm]
       h2osoi_liq(i) = min(eff_porosity(i)*dzmm(i), h2osoi_liq(i))
    enddo

    ! Sub-surface runoff and drainage 
    qflx_drain = qflx_drain + xs/dtime  &
         + hk(nlevsoi) + dhkdw(nlevsoi)*dwat(nlevsoi) ! [mm/s]

    ! Set imbalance for snow capping
    qflx_qrgwl = qflx_snowcap

    ! Implicit evaporation term is now zero
    eflx_impsoil = 0.

    ! Renew the ice and liquid mass due to condensation
    if (snl+1 >= 1) then
       h2osoi_liq(1) = h2osoi_liq(1) + qflx_dew_grnd*dtime
       h2osoi_ice(1) = h2osoi_ice(1) + (qflx_dew_snow-qflx_sub_snow)* &
            dtime
    endif

  end subroutine Drainage

end module SoilHydrologyMod
