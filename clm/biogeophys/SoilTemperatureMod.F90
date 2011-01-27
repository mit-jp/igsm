#include <misc.h>
#include <preproc.h>

module SoilTemperatureMod 

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: SoilTemperatureMod 
! 
! !DESCRIPTION: 
! Calculates snow and soil temperatures including phase change
!
! !USES:
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: SoilTemperature  ! Snow and soil temperatures including phase change
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: SoilThermProp   ! Set therm conductivities and heat cap of snow/soil layers
  private :: PhaseChange     ! Calculation of the phase change within snow and soil layers
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
! !IROUTINE: SoilTemperature
!
! !INTERFACE:
  subroutine SoilTemperature (c, xmf, fact)
!
! !DESCRIPTION: 
! Snow and soil temperatures including phase change
! o The volumetric heat capacity is calculated as a linear combination 
!   in terms of the volumetric fraction of the constituent phases. 
! o The thermal conductivity of soil is computed from 
!   the algorithm of Johansen (as reported by Farouki 1981), and the 
!   conductivity of snow is from the formulation used in
!   SNTHERM (Jordan 1991).
! o Boundary conditions:  
!   F = Rnet - Hg - LEg (top),  F= 0 (base of the soil column).
! o Soil / snow temperature is predicted from heat conduction 
!   in 10 soil layers and up to 5 snow layers. 
!   The thermal conductivities at the interfaces between two 
!   neighboring layers (j, j+1) are derived from an assumption that 
!   the flux across the interface is equal to that from the node j 
!   to the interface and the flux from the interface to the node j+1. 
!   The equation is solved using the Crank-Nicholson method and 
!   results in a tridiagonal system equation.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use globals, only : dtime
    use clm_varcon, only : sb, capr, cnfac
    use clm_varpar, only : nlevsoi
    use TridiagonalMod, only : Tridiagonal
!
! !ARGUMENTS:
    implicit none
    type (column_type),target,intent(inout):: c		!column derived type
    real(r8), intent(out) :: xmf  ! total latent heat of phase change of ground water
    real(r8), intent(out) :: fact(c%cps%snl+1 : nlevsoi)  ! used in computing tridiagonal matrix
!
! !CALLED FROM:
! subroutine Biogeophysics2 in module Biogeophysics2Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! 12/19/01, Peter Thornton
! Changed references for tg to t_grnd, for consistency with the 
! rest of the code (tg eliminated as redundant)
! 2/14/02, Peter Thornton: Migrated to new data structures. Added pft loop
! in calculation of net ground heat flux.
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
    real(r8), pointer :: forc_lwrad      !downward infrared (longwave) radiation (W/m**2)
    integer , pointer :: snl             !number of snow layers
    real(r8), pointer :: htvp            !latent heat of vapor of water (or sublimation) [j/kg]
    real(r8), pointer :: emg             !ground emissivity
    real(r8), pointer :: cgrnd           !deriv. of soil energy flux wrt to soil temp [w/m2/k]
    real(r8), pointer :: dlrad           !downward longwave radiation blow the canopy [W/m2]
    real(r8), pointer :: sabg            !solar radiation absorbed by ground (W/m**2)
    integer , pointer :: frac_veg_nosno  !fraction of vegetation not covered by snow (0 OR 1 now) [-] (new)
    real(r8), pointer :: eflx_sh_grnd    !sensible heat flux from ground (W/m**2) [+ to atm]
    real(r8), pointer :: qflx_evap_soi   !soil evaporation (mm H2O/s) (+ = to atm)
!
! local pointers to  original implicit inout arguments
!
    real(r8), pointer :: t_grnd          !ground surface temperature [K]

! local pointers to original implicit out arguments
! these two are new variables added to clmtype at the pft level
!
    real(r8), pointer :: eflx_gnet       !net ground heat flux into the surface (W/m**2)
    real(r8), pointer :: dgnetdT         !temperature derivative of ground net heat flux  
!
! local pointers to original implicit in arrays
!
    real(r8), dimension(:), pointer:: zi        !interface level below a "z" level (m)
    real(r8), dimension(:), pointer:: dz        !layer depth (m)
    real(r8), dimension(:), pointer:: z         !layer thickness (m)
    real(r8), dimension(:), pointer:: t_soisno  !soil temperature (Kelvin)
    real(r8), dimension(:), pointer:: tssbef    !soil/snow temperature before update
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer i,j                          !do loop indices
    real(r8) at(c%cps%snl+1 : nlevsoi)   !"a" vector for tridiagonal matrix
    real(r8) bt(c%cps%snl+1 : nlevsoi)   !"b" vector for tridiagonal matrix
    real(r8) ct(c%cps%snl+1 : nlevsoi)   !"c" vector for tridiagonal matrix
    real(r8) rt(c%cps%snl+1 : nlevsoi)   !"r" vector for tridiagonal solution
    real(r8) cv(c%cps%snl+1 : nlevsoi)   !heat capacity [J/(m2 K)]
    real(r8) tk(c%cps%snl+1 : nlevsoi)   !thermal conductivity [W/(m K)]
    real(r8) fn  (c%cps%snl+1 : nlevsoi) !heat diffusion through the layer interface [W/m2]
    real(r8) fn1 (c%cps%snl+1 : nlevsoi) !heat diffusion through the layer interface [W/m2]
    real(r8) dzm                         !used in computing tridiagonal matrix
    real(r8) dzp                         !used in computing tridiagonal matrix
    real(r8) hs                          !net energy flux into the surface (w/m2)
    real(r8) brr(c%cps%snl+1 : nlevsoi)  !temporary set
    real(r8) dhsdT                       !d(hs)/dT
    integer :: pi                        !pft index
    real(r8):: temp1                     !temporary variable
    real(r8):: temp2                     !temporary variable

    ! local pointers to derived subtypes
    type(atm2lnd_flux_type) , pointer :: a2lf
    type(column_pstate_type), pointer :: cps
    type(column_estate_type), pointer :: ces
    type(column_eflux_type) , pointer :: cef
    type(pft_type)          , pointer :: p
    type(pft_pstate_type)   , pointer :: pps
    type(pft_eflux_type)    , pointer :: pef
    type(pft_wflux_type)    , pointer :: pwf
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes (column-level)
    a2lf => c%a2lf
    cps => c%cps
    ces => c%ces
    cef => c%cef

    ! Assign local pointers to derived subtypes components (column-level)
    forc_lwrad => a2lf%forc_lwrad
    snl => cps%snl
    htvp => cps%htvp
    emg => cps%emg
    t_grnd => ces%t_grnd
    zi => cps%zi
    dz => cps%dz
    z => cps%z
    t_soisno => ces%t_soisno
    tssbef => ces%tssbef

    ! Ground surface and soil temperatures 

    ! Thermal conductivity and Heat capacity
    call SoilThermProp (tk, cv, c)

    ! Net ground heat flux into the surface and its temperature derivative
    ! Added a pfts loop here to get the average of hs and dhsdT over 
    ! all PFTs on the column.
    ! precalculate the terms that do not depend on PFT...
    temp1 = emg*sb*t_grnd**4
    temp2 = 4.*emg * sb * t_grnd**3

    ! begin pft loop
    hs = 0.
    dhsdT = 0.
    do pi=1,cps%npfts
       ! assign local pointers to derived subtypes (pft-level)
       p => c%p(pi)
       pps => p%pps
       pef => p%pef
       pwf => p%pwf

       ! assign local pointers to derived subtype components (pft-level)
       frac_veg_nosno => pps%frac_veg_nosno
       cgrnd => pef%cgrnd
       dlrad => pef%dlrad
       sabg => pef%sabg
       eflx_sh_grnd => pef%eflx_sh_grnd
       qflx_evap_soi => pwf%qflx_evap_soi
       eflx_gnet => pef%eflx_gnet 
       dgnetdT => pef%dgnetdT 

       ! calculate the averages
       eflx_gnet = sabg + dlrad + (1-frac_veg_nosno)*emg*forc_lwrad - temp1 &
            - (eflx_sh_grnd+qflx_evap_soi*htvp)
       dgnetdT = - cgrnd - temp2
       hs = hs + eflx_gnet * pps%wt
       dhsdT = dhsdT + dgnetdT * pps%wt

    end do ! (end of pfts loop)

    j = snl+1
    fact(j) = dtime / cv(j) * dz(j) / (0.5*(z(j)-zi(j-1)+capr*(z(j+1)-zi(j-1))))

    do j = snl+1 + 1, nlevsoi
       fact(j) = dtime/cv(j)
    enddo

    do j = snl+1, nlevsoi - 1
       fn(j) = tk(j)*(t_soisno(j+1)-t_soisno(j))/(z(j+1)-z(j))
    enddo
    fn(nlevsoi) = 0.

    ! Set up vector r and vectors a, b, c that define tridiagonal 
    ! matrix and solve system
    j     = snl+1
    dzp   = z(j+1)-z(j)
    at(j) = 0.
    bt(j) = 1+(1.-cnfac)*fact(j)*tk(j)/dzp-fact(j)*dhsdT
    ct(j) =  -(1.-cnfac)*fact(j)*tk(j)/dzp
    rt(j) = t_soisno(j) +  fact(j)*( hs - dhsdT*t_soisno(j) + cnfac*fn(j) )

    do j    = snl+1 + 1, nlevsoi - 1
       dzm   = (z(j)-z(j-1))
       dzp   = (z(j+1)-z(j))

       at(j) =   - (1.-cnfac)*fact(j)* tk(j-1)/dzm
       bt(j) = 1.+ (1.-cnfac)*fact(j)*(tk(j)/dzp + tk(j-1)/dzm)
       ct(j) =   - (1.-cnfac)*fact(j)* tk(j)/dzp

       rt(j) = t_soisno(j) + cnfac*fact(j)*( fn(j) - fn(j-1) )
    enddo

    j     =  nlevsoi
    dzm   = (z(j)-z(j-1))
    at(j) =   - (1.-cnfac)*fact(j)*tk(j-1)/dzm
    bt(j) = 1.+ (1.-cnfac)*fact(j)*tk(j-1)/dzm
    ct(j) = 0.
    rt(j) = t_soisno(j) - cnfac*fact(j)*fn(j-1)

    i = size(at)
    call Tridiagonal (i, at, bt, ct, rt, t_soisno(snl+1:nlevsoi))

    ! Melting or Freezing 
    do j = snl+1, nlevsoi - 1
       fn1(j) = tk(j)*(t_soisno(j+1)-t_soisno(j))/(z(j+1)-z(j))
    enddo
    fn1(nlevsoi) = 0.

    j = snl+1
    brr(j) = cnfac*fn(j) + (1.-cnfac)*fn1(j)

    do j = snl+1 + 1, nlevsoi
       brr(j) = cnfac*(fn(j)-fn(j-1)) + (1.-cnfac)*(fn1(j)-fn1(j-1))
    enddo

    call PhaseChange (fact(snl+1:nlevsoi), brr(snl+1:nlevsoi), hs, dhsdT, &
         tssbef(snl+1:nlevsoi), xmf, c)

    t_grnd = t_soisno(snl+1)

  end subroutine SoilTemperature

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SoilThermProp
!
! !INTERFACE:
  subroutine SoilThermProp (tk, cv, c)
!
! !DESCRIPTION: 
! Calculation of thermal conductivities and heat capacities of 
! snow/soil layers
! (1) The volumetric heat capacity is calculated as a linear combination 
!     in terms of the volumetric fraction of the constituent phases. 
!
! (2) The thermal conductivity of soil is computed from the algorithm of
!     Johansen (as reported by Farouki 1981), and of snow is from the
!     formulation used in SNTHERM (Jordan 1991).
! The thermal conductivities at the interfaces between two neighboring 
! layers (j, j+1) are derived from an assumption that the flux across 
! the interface is equal to that from the node j to the interface and the 
! flux from the interface to the node j+1. 
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use clm_varcon, only : denh2o, denice, tfrz,   tkwat, tkice, tkair, &
         cpice,  cpliq,  istice, istwet
    use clm_varpar, only : nlevsoi
!
! !ARGUMENTS:
    implicit none
    type (column_type),target,intent(inout):: c		!column derived type
    real(r8), intent(out) :: cv(c%cps%snl+1:nlevsoi) ! heat capacity [J/(m2 K)]
    real(r8), intent(out) :: tk(c%cps%snl+1:nlevsoi) ! thermal conductivity [W/(m K)]
!
! !CALLED FROM:
! subroutine SoilTemperature in this module
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! 2/13/02, Peter Thornton: migrated to new data structures
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in scalars
!
    integer , pointer :: itypwat    !water type
    integer , pointer :: snl        !number of snow layers
    real(r8), pointer :: h2osno     !snow water (mm H2O) 
!
! local pointers to original implicit in arrays
!
    real(r8), dimension(:), pointer:: watsat       !volumetric soil water at saturation (porosity)
    real(r8), dimension(:), pointer:: tksatu       !thermal conductivity, saturated soil [W/m-K]
    real(r8), dimension(:), pointer:: tkmg         !thermal conductivity, soil minerals  [W/m-K]
    real(r8), dimension(:), pointer:: tkdry        !thermal conductivity, dry soil (W/m/Kelvin)
    real(r8), dimension(:), pointer:: csol         !heat capacity, soil solids (J/m**3/Kelvin)
    real(r8), dimension(:), pointer:: dz           !layer depth (m)
    real(r8), dimension(:), pointer:: zi           !interface level below a "z" level (m)
    real(r8), dimension(:), pointer:: z            !layer thickness (m)
    real(r8), dimension(:), pointer:: t_soisno     !soil temperature (Kelvin)
    real(r8), dimension(:), pointer:: h2osoi_liq   !liquid water (kg/m2)
    real(r8), dimension(:), pointer:: h2osoi_ice   !ice lens (kg/m2)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer i                       ! do loop index
    real(r8) bw                     ! partial density of water (ice + liquid)
    real(r8) dksat                  ! thermal conductivity for saturated soil (j/(k s m))
    real(r8) dke                    ! kersten number
    real(r8) fl                     ! fraction of liquid or unfrozen water to total water
    real(r8) satw                   ! relative total water content of soil.
    real(r8) thk(c%cps%snl+1:nlevsoi) ! thermal conductivity of layer

    ! local pointers to derived subtypes
    type(landunit_pstate_type),pointer :: lps
    type(column_pstate_type)  ,pointer :: cps
    type(column_wstate_type)  ,pointer :: cws
    type(column_estate_type)  ,pointer :: ces
!-----------------------------------------------------------------------

    ! Assign local pointers to derived types (column-level)
    lps => c%cps%lps
    cps => c%cps
    cws => c%cws
    ces => c%ces

    ! Assign local pointers to derived type components (column-level)
    itypwat => lps%itypwat
    snl => cps%snl
    h2osno => cws%h2osno

    ! Assign local pointers to implicit in arrays
    watsat => lps%watsat
    tksatu => lps%tksatu
    tkmg => lps%tkmg
    tkdry => lps%tkdry
    csol => lps%csol
    dz => cps%dz
    zi => cps%zi
    z => cps%z
    t_soisno => ces%t_soisno
    h2osoi_liq =>  cws%h2osoi_liq
    h2osoi_ice =>  cws%h2osoi_ice

    ! Thermal conductivity of soil from Farouki (1981)
    do i = 1, nlevsoi
       if (itypwat/=istwet .AND. itypwat/=istice) then
          satw = (h2osoi_liq(i)/denh2o+  &
               h2osoi_ice(i)/denice)/(dz(i)*watsat(i))
          satw = min(1._r8, satw)
          if (satw > .1e-6) then
             fl = h2osoi_liq(i)/(h2osoi_ice(i)+h2osoi_liq(i))
             if(t_soisno(i) >= tfrz) then       ! Unfrozen soil
                dke = max(0._r8, log10(satw) + 1.0)
                dksat = tksatu(i)
             else                               ! Frozen soil
                dke = satw
                dksat = tkmg(i)*0.249**(fl*watsat(i))*2.29**watsat(i)
             endif
             thk(i) = dke*dksat + (1.-dke)*tkdry(i)
          else
             thk(i) = tkdry(i)
          endif
       else
          thk(i) = tkwat
          if (t_soisno(i) < tfrz) thk(i) = tkice
       endif
    enddo

    ! Thermal conductivity of snow, which from Jordan (1991) pp. 18
    if(snl+1 < 1)then
       do i = snl+1, 0
          bw = (h2osoi_ice(i)+h2osoi_liq(i))/dz(i)
          thk(i) = tkair + (7.75e-5 *bw + 1.105e-6*bw*bw) &
               *(tkice-tkair)
       enddo
    endif

    ! Thermal conductivity at the layer interface
    do i = snl+1, nlevsoi-1
       tk(i) = thk(i)*thk(i+1)*(z(i+1)-z(i)) &
            /(thk(i)*(z(i+1)-zi(i))+thk(i+1)*(zi(i)-z(i)))
    enddo
    tk(nlevsoi) = 0.

    ! Soil heat capacity, from de Vires (1963)
    do i = 1, nlevsoi
       if (itypwat/=istwet .AND. itypwat/=istice) then
          cv(i) = csol(i)*(1-watsat(i))*dz(i) +   &
               (h2osoi_ice(i)*cpice + h2osoi_liq(i)*cpliq)
       else
          cv(i) = (h2osoi_ice(i)*cpice + h2osoi_liq(i)*cpliq)
       endif
    enddo
    if (snl+1 == 1 .AND. h2osno > 0.) cv(1) = cv(1) + cpice*h2osno

    ! Snow heat capacity
    if (snl+1 < 1)then
       do i = snl+1, 0
          cv(i) = cpliq*h2osoi_liq(i) + cpice*h2osoi_ice(i)
       enddo
    endif

  end subroutine SoilThermProp

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: PhaseChange
!
! !INTERFACE:
  subroutine PhaseChange (fact,   brr, hs, dhsdT, &   
                          tssbef, xmf, c) 
!
! !DESCRIPTION: 
! Calculation of the phase change within snow and soil layers:
! (1) Check the conditions for which the phase change may take place, 
!     i.e., the layer temperature is great than the freezing point 
!     and the ice mass is not equal to zero (i.e. melting), 
!     or the layer temperature is less than the freezing point 
!     and the liquid water mass is not equal to zero (i.e. freezing).
! (2) Assess the rate of phase change from the energy excess (or deficit) 
!     after setting the layer temperature to freezing point.
! (3) Re-adjust the ice and liquid mass, and the layer temperature
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use globals, only: dtime
    use clm_varcon, only : tfrz, hfus
    use clm_varpar, only : nlevsoi
!
! !ARGUMENTS:
    implicit none
    type (column_type),target,intent(inout):: c	      !column derived type
    real(r8), intent(in) :: tssbef(c%cps%snl+1:nlevsoi) !temperature at previous time step [K]
    real(r8), intent(in) :: brr   (c%cps%snl+1:nlevsoi) !
    real(r8), intent(in) :: fact  (c%cps%snl+1:nlevsoi) !temporary variables
    real(r8), intent(in) :: hs                          !net ground heat flux into the surface
    real(r8), intent(in) :: dhsdT                       !temperature derivative of "hs"
    real(r8), intent(out):: xmf                         !total latent heat of phase change
!
! !CALLED FROM:
! subroutine SoilTemperature in this module
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! 2/14/02, Peter Thornton: Migrated to new data structures.
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in scalars
!
    integer , pointer :: snl          !number of snow layers
    real(r8), pointer :: h2osno       !snow water (mm H2O)
!
! local pointers to original implicit inout scalars
!
    real(r8), pointer :: snowdp       !snow height (m)
!
! local pointers tooriginal implicit out scalars
!
    real(r8), pointer :: qflx_snomelt !snow melt (mm H2O /s)
    real(r8), pointer :: eflx_snomelt !snow melt heat flux (W/m**2)
!
! local pointers to original implicit in arrays
!
    real(r8), dimension(:), pointer :: h2osoi_liq  !liquid water (kg/m2) (new)
    real(r8), dimension(:), pointer :: h2osoi_ice  !ice lens (kg/m2) (new)
!
! local pointers to original implicit inout arrays
!
    real(r8), dimension(:), pointer :: t_soisno !soil temperature (Kelvin)
!
! local pointers to original implicit out arrays
!
    integer, dimension(:), pointer ::  imelt !flag for melting (=1), freezing (=2), Not=0 (new)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  j                         !do loop index
    real(r8) hm(c%cps%snl+1:nlevsoi)   !energy residual [W/m2]
    real(r8) xm(c%cps%snl+1:nlevsoi)   !melting or freezing within a time step [kg/m2]
    real(r8) heatr                     !energy residual or loss after melting or freezing
    real(r8) temp1                     !temporary variables [kg/m2]
    real(r8), dimension(c%cps%snl+1:nlevsoi) :: wmass0, wice0, wliq0
    real(r8)  propor,tinc

    ! local pointers to derived subtypes
    type(column_pstate_type), pointer :: cps
    type(column_wstate_type), pointer :: cws
    type(column_estate_type), pointer :: ces
    type(column_wflux_type) , pointer :: cwf
    type(column_eflux_type) , pointer :: cef
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes (column-level)
    cps => c%cps
    cws => c%cws
    ces => c%ces
    cwf => c%cwf
    cef => c%cef

    ! Assign local pointers to derived subtypes components (column-level)
    snl => cps%snl
    h2osno => cws%h2osno
    snowdp => cps%snowdp
    qflx_snomelt => cwf%qflx_snomelt  
    eflx_snomelt => cef%eflx_snomelt  

    ! Assign local pointers to implicit in, inout, and out arrays
    h2osoi_liq => cws%h2osoi_liq
    h2osoi_ice => cws%h2osoi_ice
    imelt => cps%imelt
    t_soisno => ces%t_soisno

    ! Initialization 
    qflx_snomelt = 0.
    xmf = 0.
    do j = snl+1, nlevsoi
       imelt(j) = 0
       hm(j) = 0.
       xm(j) = 0.
       wice0(j) = h2osoi_ice(j)
       wliq0(j) = h2osoi_liq(j)
       wmass0(j) = h2osoi_ice(j) + h2osoi_liq(j)
    enddo

    ! Melting identification
    ! If ice exists above melt point, melt some to liquid.
    do j = snl+1, nlevsoi
       if (h2osoi_ice(j) > 0. .AND. t_soisno(j) > tfrz) then
          imelt(j) = 1
          t_soisno(j) = tfrz
       endif
       ! Freezing identification
       ! If liquid exists below melt point, freeze some to ice.
       if (h2osoi_liq(j) > 0. .AND. t_soisno(j) < tfrz) then
          imelt(j) = 2
          t_soisno(j) = tfrz
       endif
    enddo

    ! If snow exists, but its thickness is less than the critical value (0.01 m)
    if (snl+1 == 1 .AND. h2osno > 0.) then
       if (t_soisno(1) > tfrz) then
          imelt(1) = 1
          t_soisno(1) = tfrz
       endif
    endif

    ! Calculate the energy surplus and loss for melting and freezing
    do j = snl+1, nlevsoi
       if (imelt(j) > 0) then
          tinc = t_soisno(j)-tssbef(j)
          if (j > snl+1) then
             hm(j) = brr(j) - tinc/fact(j)
          else
             hm(j) = hs + dhsdT*tinc + brr(j) - tinc/fact(j)
          endif
       endif
    enddo

    ! These two errors were checked carefully.  They result from the 
    ! computed error of "Tridiagonal-Matrix" in subroutine "thermal".
    do j = snl+1, nlevsoi
       if (imelt(j) == 1 .AND. hm(j) < 0.) then
          hm(j) = 0.
          imelt(j) = 0
       endif

       if (imelt(j) == 2 .AND. hm(j) > 0.) then
          hm(j) = 0.
          imelt(j) = 0
       endif
    enddo

    ! The rate of melting and freezing
    do j = snl+1, nlevsoi

       if (imelt(j) > 0 .and. abs(hm(j)) > .0) then
          xm(j) = hm(j)*dtime/hfus                        ! kg/m2

          ! If snow exists, but its thickness is less than the critical value
          ! (1 cm). Note: more work is needed to determine how to tune the
          ! snow depth for this case
          if (j == 1) then
             if ((snl+1 == 1) .AND. (h2osno > 0.) .AND. (xm(j) > 0.)) then
                temp1 = h2osno                                        ! kg/m2
                h2osno = max(0._r8,temp1-xm(j))
                propor = h2osno/temp1
                snowdp = propor * snowdp
                heatr = hm(j) - hfus*(temp1-h2osno)/dtime         ! W/m2
                if (heatr > 0.) then
                   xm(j) = heatr*dtime/hfus                           ! kg/m2
                   hm(j) = heatr                                          ! W/m2
                else
                   xm(j) = 0.
                   hm(j) = 0.
                endif
                qflx_snomelt = max(0._r8,(temp1-h2osno))/dtime   ! kg/(m2 s)
                xmf = hfus*qflx_snomelt
             endif
          endif

          heatr = 0.
          if (xm(j) > 0.) then
             h2osoi_ice(j) = max(0._r8, wice0(j)-xm(j))
             heatr = hm(j) - hfus*(wice0(j)-h2osoi_ice(j))/dtime
          else if (xm(j) < 0.) then
             h2osoi_ice(j) = min(wmass0(j), wice0(j)-xm(j))
             heatr = hm(j) - hfus*(wice0(j)-h2osoi_ice(j))/dtime
          endif

          h2osoi_liq(j) = max(0._r8,wmass0(j)-h2osoi_ice(j))

          if (abs(heatr) > 0.) then
             if (j > snl+1) then
                t_soisno(j) = t_soisno(j) + fact(j)*heatr
             else
                t_soisno(j) = t_soisno(j) + fact(j)*heatr/(1.-fact(j)*dhsdT)
             endif
             if (h2osoi_liq(j)*h2osoi_ice(j)>0.) t_soisno(j) = tfrz
          endif

          xmf = xmf + hfus * (wice0(j)-h2osoi_ice(j))/dtime

          if (imelt(j) == 1 .AND. j < 1) then
             qflx_snomelt = qflx_snomelt + max(0._r8,(wice0(j)-h2osoi_ice(j)))/dtime
          endif

       endif

    enddo

    ! Needed for history file output
    eflx_snomelt = qflx_snomelt * hfus 

  end subroutine PhaseChange

end module SoilTemperatureMod 
