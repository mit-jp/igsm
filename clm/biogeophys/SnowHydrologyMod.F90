#include <misc.h>
#include <preproc.h>

module SnowHydrologyMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: SnowHydrologyMod
! 
! !DESCRIPTION: 
! Calculate snow hydrology
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: SnowWater         ! Change of snow mass and the snow water onto soil
  public :: SnowCompaction    ! Change in snow layer thickness due to compaction 
  public :: CombineSnowLayers ! Combine snow layers less than a min thickness
  public :: DivideSnowLayers  ! Subdivide snow layers if they exceed maximum thickness
!
! !PRIVATE MEMBER FUNCTIONS
  private :: Combo            ! Returns the combined variables: dz, t, wliq, wice. 
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
! !IROUTINE: SnowWater
!
! !INTERFACE:
  subroutine SnowWater (c)
!
! !DESCRIPTION: 
! Evaluate the change of snow mass and the snow water onto soil.
! Water flow within snow is computed by an explicit and non-physical 
! based scheme, which permits a part of liquid water over the holding 
! capacity (a tentative value is used, i.e. equal to 0.033*porosity) to
! percolate into the underlying layer.  Except for cases where the 
! porosity of one of the two neighboring layers is less than 0.05, zero 
! flow is assumed. The water flow out of the bottom of the snow pack will 
! participate as the input of the soil water and runoff. 
!
! !USES:
    use clmtype
    use globals, only: dtime
    use clm_varcon, only : denh2o, denice, wimp, ssi
!
! !ARGUMENTS:
    implicit none
    type (column_type),target,intent(inout):: c		!column derived type
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! 15 November 2000: Mariana Vertenstein
! 2/26/02, Peter Thornton: Migrated to new data structures.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer,pointer:: snl               !number of snow layers
    logical,pointer:: do_capsnow        !true => do snow capping
    real(r8),pointer:: qflx_snomelt     !snow melt (mm H2O /s)
    real(r8),pointer:: qflx_rain_grnd   !rain on ground after interception (mm H2O/s) [+]
    real(r8),pointer:: qflx_sub_snow    !sublimation rate from snow pack (mm H2O /s) [+]
    real(r8),pointer:: qflx_evap_grnd   !ground surface evaporation rate (mm H2O/s) [+]
    real(r8),pointer:: qflx_dew_snow    !surface dew added to snow pack (mm H2O /s) [+]
    real(r8),pointer:: qflx_dew_grnd    !ground surface dew formation (mm H2O /s) [+]
!
! local pointers to implicit out arguments
!
    real(r8),pointer:: qflx_top_soil    !net water input into soil from top (mm/s)
!
! local pointers to implicit in arrays
!
    real(r8),dimension(:),pointer:: dz  !layer depth (m)
!
! local pointers to implicit inout arrays
!
    real(r8),dimension(:),pointer:: h2osoi_ice   !ice lens (kg/m2) 
    real(r8),dimension(:),pointer:: h2osoi_liq   !liquid water (kg/m2)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  j                            !do loop/array indices
    real(r8) qin                          !water flow into the elmement (mm/s)
    real(r8) qout                         !water flow out of the elmement (mm/s)
    real(r8) wgdif                        !ice mass after minus sublimation
    real(r8) vol_liq(c%cps%snl+1 : 0)     !partial volume of liquid water in layer
    real(r8) vol_ice(c%cps%snl+1 : 0)     !partial volume of ice lens in layer
    real(r8) eff_porosity(c%cps%snl+1: 0) !effective porosity = porosity - vol_ice
    real(r8) qout_snowb                   !rate of water out of snow bottom (mm/s)
    type(column_pstate_type),pointer::  cps   !local pointers to derived subtypes
    type(column_wstate_type),pointer::  cws   !local pointers to derived subtypes
    type(column_wflux_type),pointer::   cwf   !local pointers to derived subtypes
    type(pft_wflux_type),pointer::      pwf_a !column-average over all pfts
!-----------------------------------------------------------------------

    ! assign local pointers to derived subtypes
    cps => c%cps
    cws => c%cws
    cwf => c%cwf
    pwf_a => cwf%pwf_a

    ! assign local pointers to derived type scalar members
    snl => cps%snl
    do_capsnow => cps%do_capsnow
    qflx_snomelt => cwf%qflx_snomelt
    qflx_rain_grnd => pwf_a%qflx_rain_grnd
    qflx_sub_snow => pwf_a%qflx_sub_snow
    qflx_evap_grnd => pwf_a%qflx_evap_grnd
    qflx_dew_snow => pwf_a%qflx_dew_snow
    qflx_dew_grnd => pwf_a%qflx_dew_grnd
    qflx_top_soil => cwf%qflx_top_soil

    ! assign local pointers to derived type array members
    dz => cps%dz
    h2osoi_ice => cws%h2osoi_ice
    h2osoi_liq => cws%h2osoi_liq

    ! begin calculation
    if (snl+1 >=1) then
       qflx_top_soil = qflx_rain_grnd + qflx_snomelt

    else
       ! Renew the mass of ice lens (h2osoi_ice) and liquid (h2osoi_liq) in the
       ! surface snow layer, resulting from sublimation (frost) / evaporation (condense)
       if (do_capsnow) then
          wgdif = h2osoi_ice(snl+1) - qflx_sub_snow*dtime
          h2osoi_ice(snl+1) = wgdif
          if (wgdif < 0.) then
             h2osoi_ice(snl+1) = 0.
             h2osoi_liq(snl+1) = h2osoi_liq(snl+1) + wgdif
          endif
          h2osoi_liq(snl+1) = h2osoi_liq(snl+1) - qflx_evap_grnd*dtime
          h2osoi_liq(snl+1) = max(0._r8, h2osoi_liq(snl+1))
       else
          wgdif = h2osoi_ice(snl+1) + (qflx_dew_snow - qflx_sub_snow)*dtime
          h2osoi_ice(snl+1) = wgdif
          if (wgdif < 0.) then
             h2osoi_ice(snl+1) = 0.
             h2osoi_liq(snl+1) = h2osoi_liq(snl+1) + wgdif
          endif
          h2osoi_liq(snl+1) = h2osoi_liq(snl+1) +  &
               (qflx_rain_grnd + qflx_dew_grnd - qflx_evap_grnd)*dtime
          h2osoi_liq(snl+1) = max(0._r8, h2osoi_liq(snl+1))
       endif

       ! Porosity and partial volume
       do j = snl+1, 0
          vol_ice(j) = min(1._r8, h2osoi_ice(j)/(dz(j)*denice))
          eff_porosity(j) = 1. - vol_ice(j)
          vol_liq(j) = min(eff_porosity(j),h2osoi_liq(j)/(dz(j)*denh2o))
       enddo

       ! Capillary forces within snow are usually two or more orders of magnitude
       ! less than those of gravity. Only gravity terms are considered. 
       ! the genernal expression for water flow is "K * ss**3", however, 
       ! no effective parameterization for "K".  Thus, a very simple consideration 
       ! (not physically based) is introduced: 
       ! when the liquid water of layer exceeds the layer's holding 
       ! capacity, the excess meltwater adds to the underlying neighbor layer.
       qin = 0.
       do j= snl+1, 0
          h2osoi_liq(j) = h2osoi_liq(j) + qin
          if (j <= -1) then
             ! No runoff over snow surface, just ponding on surface
             if (eff_porosity(j)<wimp .OR. eff_porosity(j+1)<wimp) then
                qout = 0.
             else
                qout = max(0._r8,(vol_liq(j)-ssi*eff_porosity(j))*dz(j))
                qout = min(qout,(1.-vol_ice(j+1)-vol_liq(j+1))*dz(j+1))
             endif
          else
             qout = max(0._r8,(vol_liq(j)-ssi*eff_porosity(j))*dz(j))
          endif
          qout = qout*1000.
          h2osoi_liq(j) = h2osoi_liq(j) - qout
          qin = qout
       enddo

       qout_snowb = qout/dtime
       qflx_top_soil = qout_snowb

    endif

  end subroutine SnowWater

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SnowCompaction
!
! !INTERFACE:
  subroutine SnowCompaction (c)
!
! !DESCRIPTION: 
! Determine the change in snow layer thickness due to compaction and 
! settling. 
! Three metamorphisms of changing snow characteristics are implemented, 
! i.e., destructive, overburden, and melt. The treatments of the former 
! two are from SNTHERM.89 and SNTHERM.99 (1991, 1999). The contribution 
! due to melt metamorphism is simply taken as a ratio of snow ice 
! fraction after the melting versus before the melting. 
!
! !USES:
    use clmtype
    use globals, only: dtime
    use clm_varcon, only : denice, denh2o, tfrz
!
! !ARGUMENTS:
    implicit none
    type (column_type),target,intent(inout):: c		!column derived type
!
! !CALLED FROM:
! subroutine Hydrology2 in module Hydrology2Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! 2/28/02, Peter Thornton: Migrated to new data structures
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
!
    integer,pointer:: snl   !number of snow layers
!
! local pointers to implicit in arrays
!
    integer,dimension(:),pointer:: imelt        !flag for melting (=1), freezing (=2), Not=0
    real(r8),dimension(:),pointer:: frac_iceold  !fraction of ice relative to the tot water
    real(r8),dimension(:),pointer:: t_soisno     !soil temperature (Kelvin)
    real(r8),dimension(:),pointer:: h2osoi_ice   !ice lens (kg/m2)
    real(r8),dimension(:),pointer:: h2osoi_liq   !liquid water (kg/m2)
!
! local pointers to implicit inout arrays
!
    real(r8),dimension(:),pointer:: dz           !layer depth (m)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer i       ! do loop index
    real(r8) c2     ! = 21e-3 [m3/kg]
    real(r8) c3     ! = 2.777e-6 [1/s]
    real(r8) c4     ! = 0.04 [1/K]
    real(r8) c5     ! = 2.0
    real(r8) dm     ! Upper Limit on Destructive Metamorphism Compaction [kg/m3]
    real(r8) eta0   ! The Viscosity Coefficient Eta0 [kg-s/m2]
    real(r8) burden ! pressure of overlying snow [kg/m2]
    real(r8) ddz1   ! Rate of settling of snowpack due to destructive metamorphism.
    real(r8) ddz2   ! Rate of compaction of snowpack due to overburden.
    real(r8) ddz3   ! Rate of compaction of snowpack due to melt [1/s]
    real(r8) dexpf  ! expf=exp(-c4*(273.15-t_soisno)).
    real(r8) fi     ! Fraction of ice relative to the total water content at current
                    ! time step
    real(r8) td     ! t_soisno - tfrz [K]
    real(r8) pdzdtc ! Nodal rate of change in fractional-thickness due to compaction
                    ! [fraction/s]
    real(r8) void   ! void (1 - vol_ice - vol_liq)
    real(r8) wx     ! water mass (ice+liquid) [kg/m2]
    real(r8) bi     ! partial density of ice [kg/m3]
    type(column_pstate_type),pointer::  cps !local pointers to derived subtypes
    type(column_estate_type),pointer::  ces !local pointers to derived subtypes
    type(column_wstate_type),pointer::  cws !local pointers to derived subtypes

    data c2,c3,c4,c5/23.e-3, 2.777e-6, 0.04, 2.0/
    data dm/100./
    data eta0/9.e5/
!-----------------------------------------------------------------------

    ! assign local pointers to derived subtypes
    cps => c%cps
    ces => c%ces
    cws => c%cws

    ! assign local pointers to derived type scalar members
    snl => cps%snl

    ! assign local pointers to derived type array members
    dz => cps%dz
    imelt => cps%imelt
    frac_iceold => cps%frac_iceold
    t_soisno => ces%t_soisno
    h2osoi_ice => cws%h2osoi_ice
    h2osoi_liq => cws%h2osoi_liq

    ! Begin calculation
    burden = 0.0

    do i = snl+1, 0

       wx = h2osoi_ice(i) + h2osoi_liq(i)
       void=1.-(h2osoi_ice(i)/denice+h2osoi_liq(i)/denh2o)/dz(i)

       ! Disallow compaction for water saturated node and lower ice lens node.
       if(void <= 0.001 .or. h2osoi_ice(i) <= .1)then
          burden = burden+wx
          cycle
       endif

       bi = h2osoi_ice(i) / dz(i)
       fi = h2osoi_ice(i) / wx
       td = tfrz-t_soisno(i)
       dexpf = exp(-c4*td)

       ! Settling as a result of destructive metamorphism
       ddz1 = -c3*dexpf
       if(bi > dm) ddz1 = ddz1*exp(-46.0e-3*(bi-dm))

       ! Liquid water term
       if(h2osoi_liq(i) > 0.01*dz(i)) ddz1=ddz1*c5

       ! Compaction due to overburden
       ddz2 = -burden*exp(-0.08*td - c2*bi)/eta0

       ! Compaction occurring during melt
       if(imelt(i) == 1)then
          ddz3 = - 1./dtime * max(0._r8,(frac_iceold(i) - fi)/frac_iceold(i))
       else
          ddz3 = 0.
       endif

       ! Time rate of fractional change in dz (units of s-1)
       pdzdtc = ddz1+ddz2+ddz3

       ! The change in dz due to compaction
       dz(i) = dz(i)*(1.+pdzdtc*dtime)

       ! Pressure of overlying snow
       burden = burden+wx

    enddo

  end subroutine SnowCompaction

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CombineSnowLayers
!
! !INTERFACE:
  subroutine CombineSnowLayers (c) 
!
! !DESCRIPTION: 
! Combine snow layers that are less than a minimum thickness or mass
! If the snow element thickness or mass is less than a prescribed minimum,
! then it is combined with a neighboring element.  The subroutine 
! clm_combo.f90 then executes the combination of mass and energy.
!
! !USES:
    use clmtype
    use clm_varcon, only : istsoil
!
! !ARGUMENTS:
    implicit none
    type (column_type),target,intent(inout):: c		!column derived type
!
! !CALLED FROM:
! subroutine Hydrology2 in module Hydrology2Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! 2/28/02, Peter Thornton: Migrated to new data structures.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
!
    integer,pointer:: itypwat  !water type
!
! local pointers to implicit inout scalars
!
    integer ,pointer:: snl     !number of snow layers
    real(r8),pointer:: h2osno  !snow water (mm H2O)
    real(r8),pointer:: snowdp  !snow height (m)
!
! local pointers to implicit inout arrays
!
    real(r8),dimension(:),pointer:: dz           !layer depth (m)
    real(r8),dimension(:),pointer:: zi           !interface level below a "z" level (m)
    real(r8),dimension(:),pointer:: t_soisno     !soil temperature (Kelvin)
    real(r8),dimension(:),pointer:: h2osoi_ice   !ice lens (kg/m2)
    real(r8),dimension(:),pointer:: h2osoi_liq   !liquid water (kg/m2)
!
! local pointers to implicit out arrays
!
    real(r8),dimension(:),pointer:: z            !layer thickness (m)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer i         ! do loop index 
    integer j         ! node index
    integer k         ! do loop index
    integer l         ! node index
    integer msn_old   ! number of snow layer 1 (top) to msn0 (bottom)
    integer mssi      ! node index
    integer neibor    ! adjacent node selected for combination
    real(r8) zwice    ! total ice mass in snow
    real(r8) zwliq    ! total liquid water in snow
    type(column_pstate_type),pointer::  cps !local pointers to derived subtypes
    type(column_estate_type),pointer::  ces !local pointers to derived subtypes 
    type(column_wstate_type),pointer::  cws !local pointers to derived subtypes 
    real(r8) dzmin(5) ! minimum of snow layer 1 (top) to msn0 (bottom)
    data dzmin /0.010, 0.015, 0.025, 0.055, 0.115/
!-----------------------------------------------------------------------

    ! assign local pointers to derived subtypes
    cps => c%cps
    ces => c%ces
    cws => c%cws

    ! assign local pointers to derived type scalar members
    itypwat => cps%lps%itypwat
    snl => cps%snl
    snowdp => cps%snowdp
    h2osno => cws%h2osno

    ! assign local pointers to derived type array members
    dz => cps%dz
    zi => cps%zi
    z => cps%z
    t_soisno => ces%t_soisno
    h2osoi_ice => cws%h2osoi_ice
    h2osoi_liq => cws%h2osoi_liq

    ! Check the mass of ice lens of snow, when the total is less than a small value,
    ! combine it with the underlying neighbor.
    msn_old = snl
    do j = msn_old+1, 0
       if(h2osoi_ice(j) <= .1)then

          if (itypwat == istsoil) then
             h2osoi_liq(j+1) = h2osoi_liq(j+1) + h2osoi_liq(j)
             h2osoi_ice(j+1) = h2osoi_ice(j+1) + h2osoi_ice(j)
          else if (itypwat /= istsoil .and. j /= 0) then
             h2osoi_liq(j+1) = h2osoi_liq(j+1) + h2osoi_liq(j)
             h2osoi_ice(j+1) = h2osoi_ice(j+1) + h2osoi_ice(j)
          endif

          ! shift all elements above this down one.
          if(j > snl+1 .AND. snl < -1)then
             do i =  j, snl+2, -1
                t_soisno(i) = t_soisno(i-1)
                h2osoi_liq(i) = h2osoi_liq(i-1)
                h2osoi_ice(i) = h2osoi_ice(i-1)
                dz(i) = dz(i-1)
             enddo
          endif
          snl = snl + 1
       endif
    enddo

    if(snl == 0)then
       h2osno = 0.
       snowdp = 0.
       return
    else
       h2osno = 0.
       snowdp = 0.
       zwice = 0.
       zwliq = 0.
       do j = snl + 1, 0
          h2osno = h2osno + h2osoi_ice(j) + h2osoi_liq(j)
          snowdp = snowdp + dz(j)
          zwice = zwice + h2osoi_ice(j)
          zwliq = zwliq + h2osoi_liq(j)
       enddo
    endif

    ! Check the snow depth
    if(snowdp < 0.01)then       !!! all snow gone
       snl = 0
       h2osno = zwice
       if(h2osno <= 0.) snowdp = 0.

       ! The liquid water assumed ponding on soil surface.
       if (itypwat == istsoil) then
          h2osoi_liq(1) = h2osoi_liq(1) + zwliq
       endif

       return
    else                        !!! snow layers combined

       ! Two or more layers 
       if(snl < -1)then
          msn_old = snl
          mssi = 1
          do i = msn_old+1, 0

             ! If top node is removed, combine with bottom neighbor.
             if(dz(i) < dzmin(mssi))then
                if(i == snl+1)then
                   neibor = i + 1

                   ! If the bottom neighbor is not snow, combine with the top neighbor.
                else if(i == 0)then
                   neibor = i - 1


                   ! If none of the above special cases apply, combine with the thinnest neighbor
                else
                   neibor = i + 1
                   if((dz(i-1)+dz(i)) < (dz(i+1)+dz(i))) neibor = i-1
                endif

                ! Node l and j are combined and stored as node j.
                if(neibor > i)then
                   j = neibor
                   l = i
                else
                   j = i
                   l = neibor
                endif

                call Combo (dz(j), h2osoi_liq(j), h2osoi_ice(j), &
                     t_soisno(j), dz(l), h2osoi_liq(l), &
                     h2osoi_ice(l), t_soisno(l) )

                ! Now shift all elements above this down one.
                if(j-1 > snl+1) then
                   do k = j-1, snl+2, -1
                      t_soisno(k) = t_soisno(k-1)
                      h2osoi_ice(k) = h2osoi_ice(k-1)
                      h2osoi_liq(k) = h2osoi_liq(k-1)
                      dz(k) = dz(k-1)
                   enddo
                endif

                snl = snl + 1

                if(snl >= -1) EXIT

                ! The layer thickness is greater than the prescribed minimum value                      
             else
                mssi = mssi + 1 
             endif
          enddo

       endif

       ! Reset the node depth and the depth of layer interface
       do k = 0, snl+1, -1
          z(k) = zi(k) - 0.5*dz(k)
          zi(k-1) = zi(k) - dz(k)
       enddo

    endif                       !!! snow layers combined

  end subroutine CombineSnowLayers

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: DivideSnowLayers 
!
! !INTERFACE:
  subroutine DivideSnowLayers (c)
!
! !DESCRIPTION: 
! Subdivides snow layers if they exceed their prescribed maximum thickness.
!
! !USES:
    use clmtype
!
! !ARGUMENTS:
    implicit none
    type (column_type),target,intent(inout):: c		!column derived type
!
! !CALLED FROM:
! subroutine Hydrology2 in module Hydrology2Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! 2/28/02, Peter Thornton: Migrated to new data structures.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit inout scalars
!
    integer,pointer:: snl
!
! local pointers to implicit inout arrays
!
    real(r8),dimension(:),pointer:: dz
    real(r8),dimension(:),pointer:: zi
    real(r8),dimension(:),pointer:: t_soisno
    real(r8),dimension(:),pointer:: h2osoi_ice
    real(r8),dimension(:),pointer:: h2osoi_liq
!
! local pointers to implicit out arrays
!
    real(r8),dimension(:),pointer:: z
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer k         ! number of do looping
    integer msno      ! number of snow layer 1 (top) to msno (bottom)
    real(r8) drr      ! thickness of the combined [m] 
    real(r8) dzsno(5) ! Snow layer thickness [m] 
    real(r8) swice(5) ! Partial volume of ice [m3/m3]
    real(r8) swliq(5) ! Partial volume of liquid water [m3/m3]
    real(r8) tsno(5)  ! Nodel temperature [K]
    real(r8) zwice    !
    real(r8) zwliq    !
    real(r8) propor   !
    type(column_pstate_type),pointer::  cps  !local pointers to derived subtypes
    type(column_estate_type),pointer::  ces  !local pointers to derived subtypes
    type(column_wstate_type),pointer::  cws  !local pointers to derived subtypes
!-----------------------------------------------------------------------

    ! assign local pointers to derived subtypes
    cps => c%cps
    ces => c%ces
    cws => c%cws

    ! assign local pointers to derived type scalar members
    snl => cps%snl

    ! assign local pointers to derived type array members
    dz => cps%dz
    zi => cps%zi
    z => cps%z
    t_soisno => ces%t_soisno
    h2osoi_ice => cws%h2osoi_ice
    h2osoi_liq => cws%h2osoi_liq

    msno = abs(snl)
    do k = 1, msno
       dzsno(k) = dz  (k + snl)
       swice(k) = h2osoi_ice(k + snl)
       swliq(k) = h2osoi_liq(k + snl)
       tsno(k)  = t_soisno (k + snl)
    enddo

    if (msno == 1) then
       if (dzsno(1) > 0.03) then
          msno = 2

          ! Specified a new snow layer
          dzsno(1) = dzsno(1)/2.
          swice(1) = swice(1)/2.
          swliq(1) = swliq(1)/2.
          dzsno(2) = dzsno(1)
          swice(2) = swice(1)
          swliq(2) = swliq(1)
          tsno(2)  = tsno(1)
       endif
    endif

    if (msno > 1) then
       if (dzsno(1) > 0.02) then
          drr = dzsno(1) - 0.02
          propor = drr/dzsno(1)
          zwice = propor*swice(1)
          zwliq = propor*swliq(1)
          propor = 0.02/dzsno(1)
          swice(1) = propor*swice(1)
          swliq(1) = propor*swliq(1)
          dzsno(1) = 0.02

          call Combo (dzsno(2), swliq(2), swice(2), tsno(2), drr, &
               zwliq, zwice, tsno(1))

          if (msno <= 2 .AND. dzsno(2) > 0.07) then

             ! Subdivided a new layer
             msno = 3
             dzsno(2) = dzsno(2)/2.
             swice(2) = swice(2)/2.
             swliq(2) = swliq(2)/2.
             dzsno(3) = dzsno(2)
             swice(3) = swice(2)
             swliq(3) = swliq(2)
             tsno(3)  = tsno(2)
          endif
       endif
    endif

    if (msno > 2) then
       if (dzsno(2) > 0.05) then
          drr = dzsno(2) - 0.05
          propor = drr/dzsno(2)
          zwice = propor*swice(2)
          zwliq = propor*swliq(2)
          propor = 0.05/dzsno(2)
          swice(2) = propor*swice(2)
          swliq(2) = propor*swliq(2)
          dzsno(2) = 0.05

          call Combo (dzsno(3), swliq(3), swice(3), tsno(3), drr, &
               zwliq, zwice, tsno(2))

          if (msno <= 3 .AND. dzsno(3) > 0.18) then

             ! Subdivided a new layer
             msno =  4
             dzsno(3) = dzsno(3)/2.
             swice(3) = swice(3)/2.
             swliq(3) = swliq(3)/2.
             dzsno(4) = dzsno(3)
             swice(4) = swice(3)
             swliq(4) = swliq(3)
             tsno(4)  = tsno(3)
          endif
       endif
    endif

    if (msno > 3) then
       if (dzsno(3) > 0.11) then
          drr = dzsno(3) - 0.11
          propor = drr/dzsno(3)
          zwice = propor*swice(3)
          zwliq = propor*swliq(3)
          propor = 0.11/dzsno(3)
          swice(3) = propor*swice(3)
          swliq(3) = propor*swliq(3)
          dzsno(3) = 0.11

          call Combo (dzsno(4), swliq(4), swice(4), tsno(4), drr, &
               zwliq, zwice, tsno(3))

          if (msno <= 4 .AND. dzsno(4) > 0.41) then

             ! Subdivided a new layer
             msno = 5
             dzsno(4) = dzsno(4)/2.
             swice(4) = swice(4)/2.
             swliq(4) = swliq(4)/2.
             dzsno(5) = dzsno(4)
             swice(5) = swice(4)
             swliq(5) = swliq(4)
             tsno(5)  = tsno(4)
          endif
       endif
    endif

    if (msno > 4) then
       if (dzsno(4) > 0.23) then
          drr = dzsno(4) - 0.23
          propor = drr/dzsno(4)
          zwice = propor*swice(4)
          zwliq = propor*swliq(4)
          propor = 0.23/dzsno(4)
          swice(4) = propor*swice(4)
          swliq(4) = propor*swliq(4)
          dzsno(4) = 0.23

          call Combo (dzsno(5), swliq(5), swice(5), tsno(5), drr, &
               zwliq, zwice, tsno(4))

       endif
    endif

    snl = - msno

    do k = snl+1, 0
       dz(k)   = dzsno(k - snl)
       h2osoi_ice(k) = swice(k - snl)
       h2osoi_liq(k) = swliq(k - snl)
       t_soisno(k)  = tsno (k - snl)
    enddo

    do k = 0, snl+1, -1
       z(k)    = zi(k) - 0.5*dz(k)
       zi(k-1) = zi(k) - dz(k)
    enddo

  end subroutine DivideSnowLayers

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Combo
!
! !INTERFACE:
  subroutine Combo (dz,  wliq,  wice, t, dz2, &
                    wliq2, wice2, t2) 
!
! !DESCRIPTION: 
! Combines two elements and returns the following combined
! variables: dz, t, wliq, wice. 
! The combined temperature is based on the equation:
! the sum of the enthalpies of the two elements =  
! that of the combined element.
!
! !USES:
    use clm_varcon,  only : cpice, cpliq, tfrz, hfus         
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: dz2     ! nodal thickness of 2 elements being combined [m]
    real(r8), intent(in) :: wliq2   ! liquid water of element 2 [kg/m2]
    real(r8), intent(in) :: wice2   ! ice of element 2 [kg/m2]
    real(r8), intent(in) :: t2      ! nodal temperature of element 2 [K]
    real(r8), intent(inout) :: dz   ! nodal thickness of 1 elements being combined [m]
    real(r8), intent(inout) :: wliq ! liquid water of element 1
    real(r8), intent(inout) :: wice ! ice of element 1 [kg/m2]
    real(r8), intent(inout) :: t    ! nodel temperature of elment 1 [K]
!
! !CALLED FROM:
! subroutine CombineSnowLayers in this module 
! subroutine DivideSnowLayers in this module
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!EOP
!
! !LOCAL VARIABLES:
!
    real(r8) dzc   ! Total thickness of nodes 1 and 2 (dzc=dz+dz2).
    real(r8) wliqc ! Combined liquid water [kg/m2]
    real(r8) wicec ! Combined ice [kg/m2]
    real(r8) tc    ! Combined node temperature [K]
    real(r8) h     ! enthalpy of element 1 [J/m2]
    real(r8) h2    ! enthalpy of element 2 [J/m2]
    real(r8) hc    ! temporary
!-----------------------------------------------------------------------

    dzc = dz+dz2
    wicec = (wice+wice2)
    wliqc = (wliq+wliq2)
    h =(cpice*wice+cpliq*wliq)* &
         (t-tfrz)+hfus*wliq
    h2=(cpice*wice2+cpliq*wliq2)* &
         (t2-tfrz)+hfus*wliq2

    hc = h + h2
    if(hc < 0.)then
       tc = tfrz + hc/(cpice* &
            wicec+cpliq*wliqc)
    else if(hc.le.hfus*wliqc)then
       tc = tfrz
    else
       tc = tfrz + (hc - hfus*wliqc)/ &
            (cpice*wicec+cpliq*wliqc)
    endif

    dz = dzc
    wice = wicec 
    wliq = wliqc
    t = tc

  end subroutine Combo

end module SnowHydrologyMod








