#include <misc.h>
#include <preproc.h>

module SurfaceAlbedoMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: SurfaceAlbedoMod
! 
! !DESCRIPTION: 
! Performs surface albedo calculations
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: SurfaceAlbedo  ! Surface albedo and two-stream fluxes
  public :: SnowAge        ! Update snow age
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: SnowAlbedo    ! Determine snow albedos
  private :: SoilAlbedo    ! Determine ground surface albedo 
  private :: TwoStream     ! Two-stream fluxes for canopy radiative transfer 
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
! !IROUTINE: SurfaceAlbedo
!
! !INTERFACE:
  subroutine SurfaceAlbedo (c, caldayp1, eccen, obliqr, lambm0, mvelpp)
!
! !DESCRIPTION: 
! Surface albedo and two-stream fluxes
! Surface albedos. Also fluxes (per unit incoming direct and diffuse
! radiation) reflected, transmitted, and absorbed by vegetation. 
! Also sunlit fraction of the canopy. 
! The calling sequence is:
! -> SurfaceAlbedo:     albedos for next time step
!    -> SnowAlbedo:   snow albedos: direct beam
!    -> SnowAlbedo:   snow albedos: diffuse
!    -> SoilAlbedo:   soil/lake/glacier/wetland albedos
!    -> TwoStream:    absorbed, reflected, transmitted solar fluxes (vis dir)
!    -> TwoStream:    absorbed, reflected, transmitted solar fluxes (vis dif)
!    -> TwoStream:    absorbed, reflected, transmitted solar fluxes (nir dir)
!    -> TwoStream:    absorbed, reflected, transmitted solar fluxes (nir dif)
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use shr_orb_mod
!
! !ARGUMENTS:
    implicit none
    type (column_type),target,intent(inout):: c		!column derived type
    real(r8), intent(in) :: caldayp1 ! calendar day at Greenwich (1.00, ..., 365.99)
    real(r8), intent(in) :: eccen    ! Earth's orbital eccentricity
    real(r8), intent(in) :: obliqr   ! Earth's obliquity in radians
    real(r8), intent(in) :: lambm0   ! Mean longitude of perihelion at the vernal equinox
                                     ! (radians)
    real(r8), intent(in) :: mvelpp   ! Earth's moving vernal equinox long. of perihelion
!
! !CALLED FROM:
! subroutine lpjreset1 in module DGVMMod (only applicable when cpp token DGVM is defined)
! subroutine driver   
! subroutine iniTimeVar
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! 2/1/02, Peter Thornton: Migrate to new data structures
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
!
    real(r8), pointer :: lat    !latitude (radians)
    real(r8), pointer :: lon    !longitude (radians)
    real(r8), pointer :: elai   !one-sided leaf area index with burying by snow
    real(r8), pointer :: esai   !one-sided stem area index with burying by snow
    real(r8), pointer :: h2osno !snow water (mm H2O)
    real(r8), pointer :: snowage   !non dimensional snow age [-]
!
! local pointers toimplicit out scalars
!
    real(r8), pointer :: fsun   !sunlit fraction of canopy
!
! local pointers toimplicit in arrays
!
    real(r8), dimension(:), pointer :: rhol    !leaf reflectance: 1=vis, 2=nir
    real(r8), dimension(:), pointer :: rhos    !stem reflectance: 1=vis, 2=nir
    real(r8), dimension(:), pointer :: taul    !leaf transmittance: 1=vis, 2=nir
    real(r8), dimension(:), pointer :: taus    !stem transmittance: 1=vis, 2=nir 
!
! local pointers toimplicit out arrays
!
    real(r8), dimension(:), pointer :: albgrd  !ground albedo (direct)
    real(r8), dimension(:), pointer :: albgri  !ground albedo (diffuse)
    real(r8), dimension(:), pointer :: albd    !surface albedo (direct)
    real(r8), dimension(:), pointer :: albi    !surface albedo (diffuse)
    real(r8), dimension(:), pointer :: fabd    !flux absorbed by veg per unit direct flux
    real(r8), dimension(:), pointer :: fabi    !flux absorbed by veg per unit diffuse flux
    real(r8), dimension(:), pointer :: ftdd    !down direct flux below veg per unit dir flx
    real(r8), dimension(:), pointer :: ftid    !down diffuse flux below veg per unit dir flx
    real(r8), dimension(:), pointer :: ftii    !down diffuse flux below veg per unit dif flx
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: ib              !band index
    integer  :: ic              !0=unit incoming direct; 1=unit incoming diffuse
    real(r8) :: wl              !fraction of LAI+SAI that is LAI
    real(r8) :: ws              !fraction of LAI+SAI that is SAI
    real(r8) :: mpe = 1.e-06    !prevents overflow for division by zero
    real(r8) :: vai             !elai+esai
    real(r8) :: rho(numrad)     !leaf/stem refl weighted by fraction LAI and SAI
    real(r8) :: tau(numrad)     !leaf/stem tran weighted by fraction LAI and SAI
    real(r8) :: ftdi(numrad)    !down direct flux below veg per unit dif flux = 0
    real(r8) :: albsnd(numrad)  !snow albedo (direct)
    real(r8) :: albsni(numrad)  !snow albedo (diffuse)
    real(r8) :: gdir            !aver projected leaf/stem area in solar direction
    real(r8) :: ext             !optical depth direct beam per unit LAI+SAI
    real(r8) :: delta           !solar declination angle in radians
    real(r8) :: eccf            !earth orbit eccentricity factor
    real(r8) :: coszen          !cosine solar zenith angle for next time step
    integer  :: pi              !pft index

    ! local pointers to derived subtypes
    type(column_pstate_type), pointer :: cps
    type(column_wstate_type), pointer :: cws
    type(pft_type)          , pointer :: p
    type(pft_pstate_type)   , pointer :: pps
    type(pft_epc_type)      , pointer :: pepc
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes (column-level)
    cps => c%cps
    cws => c%cws

    ! Assign local pointers to derived subtypes components (column-level)
    lat => cps%lps%gps%lat
    lon => cps%lps%gps%lon
    h2osno => cws%h2osno
    snowage => cps%snowage

    ! Assign local pointers to column-level implicit out arrays
    albgrd => cps%albgrd
    albgri => cps%albgri

    ! Perform column-level calculations

    ! Note: call to shr_orb_decl can be moved into driver, and just 
    ! pass delta into SurfaceAlbedo() as an argument. This changes only
    ! with timestep, not with gridcell.
    ! Solar declination  for next time step

    call shr_orb_decl (caldayp1, eccen, mvelpp, lambm0, obliqr, &
         delta, eccf)

    ! Actually, coszen only needs to be calculated once per gridcell, so
    ! save a lot of calls if done inside gridcell loop in driver, and
    ! just pass in coszen to SurfaceAlbedo()                  
    ! Cosine solar zenith angle for next time step

    coszen = shr_orb_cosz (caldayp1, lat, lon, delta)

    ! Initialize output because solar radiation only done if coszen > 0
    ! Note: not using the local array pointers in this pft loop because it
    ! is so short...

    do ib = 1, numrad
       albgrd(ib) = 0._r8
       albgri(ib) = 0._r8
       do pi = 1, cps%npfts
          pps => c%p(pi)%pps
          pps%albd(ib)   = 1.
          pps%albi(ib)   = 1.
          pps%fabd(ib)   = 0._r8
          pps%fabi(ib)   = 0._r8
          pps%ftdd(ib)   = 0._r8
          pps%ftid(ib)   = 0._r8
          pps%ftii(ib)   = 0._r8
          if (ib==1) pps%fsun = 0.
       end do
    end do

    ! Return if coszen is not positive 

    if (coszen <= 0._r8) RETURN 

    ! Snow albedos: only if h2osno > 0

    if (h2osno > 0._r8 ) then
       ic=0; call SnowAlbedo (coszen, snowage, ic, albsnd)
       ic=1; call SnowAlbedo (coszen, snowage, ic, albsni)
    else
       albsnd(:) = 0._r8
       albsni(:) = 0._r8
    endif

    ! Ground surface albedos

    call SoilAlbedo (c, coszen, albsnd, albsni) 

    ! Begin pft loop

    do pi = 1, cps%npfts 

       ! Assign local pointers to derived subtypes
       p => c%p(pi)
       pps => p%pps
       pepc => p%pepc

       ! Assign local pointers to input arrays
       rhol => pepc%rhol
       rhos => pepc%rhos
       taul => pepc%taul
       taus => pepc%taus

       ! Assign local pointers to input scalars
       elai => pps%elai
       esai => pps%esai  
       fsun => pps%fsun 

       ! Assign local pointers to output arrays
       albd => pps%albd
       albi => pps%albi
       fabd => pps%fabd
       fabi => pps%fabi
       ftdd => pps%ftdd
       ftid => pps%ftid
       ftii => pps%ftii

       vai = elai + esai
       if (vai /= 0.) then  ! vegetated patch

          ! Weight reflectance/transmittance by lai and sai

          wl = elai / max( vai,mpe )
          ws = esai / max( vai,mpe )
          do ib = 1, numrad
             rho(ib) = max( rhol(ib)*wl + rhos(ib)*ws, mpe )
             tau(ib) = max( taul(ib)*wl + taus(ib)*ws, mpe )
          end do

          ! Calculate surface albedos and flux

          call TwoStream (p, coszen, vai, rho, tau, albgrd, albgri, &
               fabd, albd, ftdd, ftid, &
               fabi, albi, ftdi, ftii, gdir)

          ! Sunlit fraction of canopy. Set fsun = 0 if fsun < 0.01.

          ext = gdir/coszen * sqrt(1.-rho(1)-tau(1))
          fsun = (1.-exp(-ext*vai)) / max(ext*vai,mpe)
          ext = fsun                                       !temporary fsun
          if (ext < 0.01) then
             wl = 0._r8                                        !temporary fsun
          else
             wl = ext                                          !temporary fsun
          end if
          fsun = wl

       else     ! non-vegetated patch

          do ib = 1,numrad
             fabd(ib) = 0.
             fabi(ib) = 0.
             ftdd(ib) = 1.
             ftid(ib) = 0.
             ftii(ib) = 1.
             albd(ib) = albgrd(ib)
             albi(ib) = albgri(ib)
             fsun     = 0.
          end do

       endif

    end do ! end of pft loop

    return
  end subroutine SurfaceAlbedo

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SnowAlbedo
!
! !INTERFACE:
  subroutine SnowAlbedo (coszen, snowage, ind, alb)
!
! !DESCRIPTION: 
! Determine snow albedos
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: coszen  ! cosine solar zenith angle for next time step
    real(r8), intent(in) :: snowage ! non dimensional snow age [-]
    integer , intent(in) :: ind     ! 0=direct beam, 1=diffuse radiation
    real(r8), intent(out):: alb(2)  ! snow albedo by waveband (assume 2 wavebands)
!
! !CALLED FROM:
! subroutine SurfaceAlbedo in this module
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! 2/5/02, Peter Thornton: Migrated to new data structures. Eliminated
! reference to derived types in this subroutine, and made consistent use
! of the assumption that numrad = 2, with index values: 1=visible,2=NIR
!
!EOP
!
! !LOCAL VARIABLES:
!
! variables and constants for snow albedo calculation
    integer  :: ib          !waveband class
    real(r8) :: snal0 = 0.95 !vis albedo of new snow for sza<60
    real(r8) :: snal1 = 0.65 !nir albedo of new snow for sza<60
    real(r8) :: cons  = 0.2  !constant for visible snow albedo calculation [-]
    real(r8) :: conn  = 0.5  !constant for nir snow albedo calculation [-]
    real(r8) :: sl    = 2.0  !factor that helps control alb zenith dependence [-]
    real(r8) :: age          !factor to reduce visible snow alb due to snow age [-]
    real(r8) :: albs         !temporary vis snow albedo
    real(r8) :: albl         !temporary nir snow albedo
    real(r8) :: cff          !snow alb correction factor for zenith angle > 60 [-]
    real(r8) :: czf          !solar zenith correction for new snow albedo [-]
!-----------------------------------------------------------------------

    ! this code assumes that numrad = 2 , with the following
    ! index values: 1 = visible, 2 = NIR

    ! Zero albedos
    alb(1) = 0._r8
    alb(2) = 0._r8

    ! CLM Albedo for snow cover.
    ! Snow albedo depends on snow-age, zenith angle, and thickness of snow,
    ! age gives reduction of visible radiation

    ! Correction for snow age
    age = 1.-1./(1.+snowage)
    albs = snal0*(1.-cons*age)
    albl = snal1*(1.-conn*age)

    if (ind == 0) then
       ! czf corrects albedo of new snow for solar zenith
       cff    = ((1.+1./sl)/(1.+max(0.001_r8,coszen)*2.*sl )- 1./sl)
       cff    = max(cff,0._r8)
       czf    = 0.4*cff*(1.-albs)
       albs = albs+czf
       czf    = 0.4*cff*(1.-albl)
       albl = albl+czf
    endif

    alb(1) = albs
    alb(2) = albl

    return
  end subroutine SnowAlbedo

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SoilAlbedo
!
! !INTERFACE:
  subroutine SoilAlbedo (c, coszen, albsnd, albsni)
!
! !DESCRIPTION: 
! Determine ground surface albedo, accounting for snow
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use clm_varpar, only : numrad
    use clm_varcon, only : albsat, albdry, alblak, albice, tfrz, istice, istsoil
!
! !ARGUMENTS:
    implicit none
    type (column_type), target, intent(inout) :: c  !column derived type
    real(r8), intent(in) :: coszen          !cos solar zenith angle next time step
    real(r8), intent(in) :: albsnd(numrad)  !snow albedo (direct)
    real(r8), intent(in) :: albsni(numrad)  !snow albedo (diffuse)
!
! !CALLED FROM:
! subroutine SurfaceAlbedo in this module
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! 2/5/02, Peter Thornton: Migrated to new data structures.
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in scalars
!
    integer , pointer :: itypwat    !water type
    integer , pointer :: isoicol   !soil color class
    real(r8), pointer :: t_grnd    !ground temperature (Kelvin)
    real(r8), pointer :: frac_sno  !fraction of ground covered by snow (0 to 1)
!
! local pointers to original implicit in arrays
!
    real(r8), dimension(:), pointer:: h2osoi_vol !volumetric soil water [m3/m3]
!
! local pointers to original implicit out arrays
!
    real(r8), dimension(:), pointer:: albgrd  !ground albedo (direct)
    real(r8), dimension(:), pointer:: albgri  !ground albedo (diffuse)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: nband =numrad !number of solar radiation waveband classes
    integer  :: ib            !waveband number (1=vis, 2=nir)
    real(r8) :: inc           !soil water correction factor for soil albedo
    real(r8) :: albsod        !soil albedo (direct)
    real(r8) :: albsoi        !soil albedo (diffuse)

    ! local pointers to derived subtypes
    type(column_pstate_type), pointer :: cps
    type(column_estate_type), pointer :: ces
    type(column_wstate_type), pointer :: cws
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes (column-level)
    cps => c%cps
    ces => c%ces
    cws => c%cws

    ! Assign local pointers to derived subtypes components (column-level)
    itypwat => cps%lps%itypwat
    isoicol => cps%lps%isoicol
    t_grnd => ces%t_grnd
    frac_sno => cps%frac_sno

    ! Assign local pointers to implicit in arrays
    h2osoi_vol => cws%h2osoi_vol

    ! Assign local pointers to implicit out arrays
    albgrd => cps%albgrd
    albgri => cps%albgri

    do ib = 1, nband
       if (itypwat == istsoil)  then               !soil
          inc    = max(0.11-0.40*h2osoi_vol(1), 0._r8)
          albsod = min(albsat(isoicol,ib)+inc, albdry(isoicol,ib))
          albsoi = albsod
       else if (itypwat == istice)  then           !land ice
          albsod = albice(ib)
          albsoi = albsod
       else if (t_grnd > tfrz) then                !unfrozen lake, wetland
          albsod = 0.05/(max(0.001_r8,coszen) + 0.15)
          albsoi = albsod
       else                                        !frozen lake, wetland
          albsod = alblak(ib)
          albsoi = albsod
       end if

       albgrd(ib) = albsod*(1.-frac_sno) + albsnd(ib)*frac_sno
       albgri(ib) = albsoi*(1.-frac_sno) + albsni(ib)*frac_sno

    end do

    return
  end subroutine SoilAlbedo

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: TwoStream
!
! !INTERFACE:
  subroutine TwoStream  (p, coszen, vai, rho, tau, albgrd, albgri, &
       fabd, albd, ftdd, ftid, &
       fabi, albi, ftdi, ftii, gdir)
!
! !DESCRIPTION: 
! Two-stream fluxes for canopy radiative transfer
! Use two-stream approximation of Dickinson (1983) Adv Geophysics
! 25:305-353 and Sellers (1985) Int J Remote Sensing 6:1335-1372
! to calculate fluxes absorbed by vegetation, reflected by vegetation,
! and transmitted through vegetation for unit incoming direct or diffuse 
! flux given an underlying surface with known albedo.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use clm_varpar, only : numrad
    use clm_varcon, only : omegas, tfrz, betads, betais
!
! !ARGUMENTS:
    implicit none
    type (pft_type), target, intent(inout) :: p	 !pft derived type
    real(r8), intent(in)  :: coszen          !cosine solar zenith angle for next time step
    real(r8), intent(in)  :: vai             !elai+esai
    real(r8), intent(in)  :: rho(numrad)     !leaf/stem refl weighted by fraction LAI and SAI
    real(r8), intent(in)  :: tau(numrad)     !leaf/stem tran weighted by fraction LAI and SAI
    real(r8), intent(in)  :: albgrd(numrad)  !ground albedo (direct)
    real(r8), intent(in)  :: albgri(numrad)  !ground albedo (diffuse)
    real(r8), intent(out) :: albd(numrad)    !surface albedo (direct)
    real(r8), intent(out) :: albi(numrad)    !surface albedo (diffuse)
    real(r8), intent(out) :: fabd(numrad)    !flux absorbed by veg per unit direct flux
    real(r8), intent(out) :: fabi(numrad)    !flux absorbed by veg per unit diffuse flux
    real(r8), intent(out) :: ftdd(numrad)    !down direct flux below veg per unit dir flx
    real(r8), intent(out) :: ftdi(numrad)    !down direct flux below veg per unit dif flux 
    real(r8), intent(out) :: ftid(numrad)    !down diffuse flux below veg per unit dir flx
    real(r8), intent(out) :: ftii(numrad)    !down diffuse flux below veg per unit dif flx
    real(r8), intent(out) :: gdir            !aver projected leaf/stem area in solar direction
!
! !CALLED FROM:
! subroutine SurfaceAlbedo in this module
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! Modified for speedup: Mariana Vertenstein, 8/26/02
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
!
    real(r8), pointer :: xl     !leaf/stem orientation index
    real(r8), pointer :: t_veg  !vegetation temperature (Kelvin)
    real(r8), pointer :: fwet   !fraction of canopy that is wet (0 to 1)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  i,j      ! array indices
    integer  ic       !0=unit incoming direct; 1=unit incoming diffuse
    integer  ib       !waveband number
    real(r8) cosz     ! 0.001 <= coszen <= 1.000
    real(r8) asu      ! single scattering albedo
    real(r8) chil     ! -0.4 <= xl <= 0.6
    real(r8) twostext ! optical depth of direct beam per unit leaf area
    real(r8) avmu     ! average diffuse optical depth
    real(r8) omega    ! fraction of intercepted radiation that is scattered
    real(r8) omegal   ! omega for leaves
    real(r8) betai    ! upscatter parameter for diffuse radiation
    real(r8) betail   ! betai for leaves
    real(r8) betad    ! upscatter parameter for direct beam radiation
    real(r8) betadl   ! betad for leaves

    ! local temporaries
    real(r8) tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9
    real(r8) p1,p2,p3,p4,s1,s2,u1,u2,u3
    real(r8) b,c1,d,d1,d2,f,h,h1,h2,h3,h4,h5,h6,h7,h8,h9,h10
    real(r8) phi1,phi2,sigma
    real(r8) temp0,temp1,temp2

    ! local pointers to derived subtypes
    type(pft_pstate_type)   , pointer :: pps
    type(pft_epc_type)      , pointer :: pepc
    type(pft_estate_type)   , pointer :: pes
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes (pft-level)

    pps => p%pps
    pepc => p%pepc
    pes => p%pes

    ! Assign local pointers to derived subtypes components (pft-level)

    fwet => pps%fwet
    xl => pepc%xl
    t_veg => pes%t_veg

    ! Calculate two-stream parameters omega, betad, betai, avmu, gdir, twostext.
    ! Omega, betad, betai are adjusted for snow. Values for omega*betad 
    ! and omega*betai are calculated and then divided by the new omega
    ! because the product omega*betai, omega*betad is used in solution. 
    ! Also, the transmittances and reflectances (tau, rho) are linear 
    ! weights of leaf and stem values.

    cosz = max(0.001_r8, coszen)
    chil = min( max(xl, -0.4_r8), 0.6_r8 )
    if (abs(chil) <= 0.01) chil = 0.01
    phi1 = 0.5 - 0.633*chil - 0.330*chil*chil
    phi2 = 0.877 * (1.-2.*phi1)
    gdir = phi1 + phi2*cosz
    twostext = gdir/cosz
    avmu = ( 1. - phi1/phi2 * log((phi1+phi2)/phi1) ) / phi2
    temp0 = gdir + phi2*cosz
    temp1 = phi1*cosz
    temp2 = ( 1. - temp1/temp0 * log((temp1+temp0)/temp1) )

    do ib = 1, numrad

       omegal = rho(ib) + tau(ib)
       asu = 0.5*omegal*gdir/temp0 *temp2 
       betadl = (1.+avmu*twostext)/(omegal*avmu*twostext)*asu
       betail = 0.5 * ((rho(ib)+tau(ib)) + (rho(ib)-tau(ib)) * ((1.+chil)/2.)**2) / omegal

       ! Adjust omega, betad, and betai for intercepted snow

       if (t_veg > tfrz) then                             !no snow
          tmp0 = omegal           
          tmp1 = betadl 
          tmp2 = betail  
       else
          tmp0 =   (1.-fwet)*omegal        + fwet*omegas(ib)           
          tmp1 = ( (1.-fwet)*omegal*betadl + fwet*omegas(ib)*betads ) / tmp0
          tmp2 = ( (1.-fwet)*omegal*betail + fwet*omegas(ib)*betais ) / tmp0
       end if
       omega = tmp0           
       betad = tmp1 
       betai = tmp2  

       ! Absorbed, reflected, transmitted fluxes per unit incoming radiation

       b = 1. - omega + omega*betai
       c1 = omega*betai
       tmp0 = avmu*twostext
       d = tmp0 * omega*betad
       f = tmp0 * omega*(1.-betad)
       tmp1 = b*b - c1*c1
       h = sqrt(tmp1) / avmu
       sigma = tmp0*tmp0 - tmp1
       p1 = b + avmu*h
       p2 = b - avmu*h
       p3 = b + tmp0
       p4 = b - tmp0
       s1 = exp(-h*vai)
       s2 = exp(-twostext*vai)

       ! Determine fluxes for vegetated patch for unit incoming direct 
       ! Loop over incoming direct and incoming diffuse
       ! 0=unit incoming direct; 1=unit incoming diffuse

       do ic = 0,1

          if (ic == 0) then
             u1 = b - c1/albgrd(ib)
             u2 = b - c1*albgrd(ib)
             u3 = f + c1*albgrd(ib)
          else
             u1 = b - c1/albgri(ib)
             u2 = b - c1*albgri(ib)
             u3 = f + c1*albgri(ib)
          end if

          tmp2 = u1 - avmu*h
          tmp3 = u1 + avmu*h
          d1 = p1*tmp2/s1 - p2*tmp3*s1
          tmp4 = u2 + avmu*h
          tmp5 = u2 - avmu*h
          d2 = tmp4/s1 - tmp5*s1
          h1 = -d*p4 - c1*f
          tmp6 = d - h1*p3/sigma
          tmp7 = ( d - c1 - h1/sigma*(u1+tmp0) ) * s2
          h2 = ( tmp6*tmp2/s1 - p2*tmp7 ) / d1
          h3 = - ( tmp6*tmp3*s1 - p1*tmp7 ) / d1
          h4 = -f*p3 - c1*d
          tmp8 = h4/sigma
          tmp9 = ( u3 - tmp8*(u2-tmp0) ) * s2
          h5 = - ( tmp8*tmp4/s1 + tmp9 ) / d2
          h6 = ( tmp8*tmp5*s1 + tmp9 ) / d2
          h7 = (c1*tmp2) / (d1*s1)
          h8 = (-c1*tmp3*s1) / d1
          h9 = tmp4 / (d2*s1)
          h10 = (-tmp5*s1) / d2

          if (ic == 0) then

             ! Downward direct and diffuse fluxes below vegetation
             ftdd(ib) = s2
             ftid(ib) = h4*s2/sigma + h5*s1 + h6/s1

             ! Flux reflected by vegetation
             albd(ib) = h1/sigma + h2 + h3

             ! Flux absorbed by vegetation
             fabd(ib) = 1. - albd(ib) - (1.-albgrd(ib))*ftdd(ib) - (1.-albgri(ib))*ftid(ib)

          else

             ! Downward direct and diffuse fluxes below vegetation
             ftdi(ib) = 0._r8
             ftii(ib) = h9*s1 + h10/s1

             ! Flux reflected by vegetation
             albi(ib) = h7 + h8

             ! Flux absorbed by vegetation
             fabi(ib) = 1. - albi(ib) - (1.-albgrd(ib))*ftdi(ib) - (1.-albgri(ib))*ftii(ib)

          end if

       end do   ! end of ic loop

    end do   ! end of ib loop

    return
  end subroutine TwoStream

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SnowAge
!
! !INTERFACE:
  subroutine SnowAge (c)
!
! !DESCRIPTION: 
! Updates snow age Based on BATS code.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use globals, only: dtime
    use clm_varcon, only : tfrz
!
! !ARGUMENTS:
    implicit none
    type (column_type),target,intent(inout):: c		!column derived type
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Original Code:  Robert Dickinson
! 15 September 1999: Yongjiu Dai; Integration of code into CLM
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! 3/4/02, Peter Thornton: Migrated to new data structures.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
!
    real(r8),pointer:: t_grnd     !ground temperature (Kelvin)
    real(r8),pointer:: h2osno     !snow water (mm H2O)
    real(r8),pointer:: h2osno_old !snow mass for previous time step (kg/m2)
!
! local pointers to implicit inout scalars
!
    real(r8),pointer:: snow_age    !non dimensional snow age [-]
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    real(r8) age1 ! snow aging factor due to crystal growth [-]
    real(r8) age2 ! snow aging factor due to surface growth [-]
    real(r8) age3 ! snow aging factor due to accum of other particles [-]
    real(r8) arg  ! temporary variable used in snow age calculation [-]
    real(r8) arg2 ! temporary variable used in snow age calculation [-]
    real(r8) dela ! temporary variable used in snow age calculation [-]
    real(r8) dels ! temporary variable used in snow age calculation [-]
    real(r8) sge  ! temporary variable used in snow age calculation [-]

    ! local pointers to derived subtypes
    type(column_pstate_type),pointer::  cps
    type(column_estate_type),pointer::  ces
    type(column_wstate_type),pointer::  cws
!-----------------------------------------------------------------------

    ! assign local pointers to derived subtypes
    cps => c%cps
    ces => c%ces
    cws => c%cws

    ! assign local pointers to derived type scalar members
    snow_age => cps%snowage
    t_grnd => ces%t_grnd
    h2osno => cws%h2osno
    h2osno_old => cws%h2osno_old

    ! Begin calculation

    if (h2osno <= 0.) then

       snow_age = 0.

    else if (h2osno > 800.) then       ! Over Antarctica

       snow_age = 0.

    else                               ! Away from Antarctica

       age3  = 0.3
       arg   = 5.e3*(1./tfrz-1./t_grnd)
       arg2  = min(0._r8,10.*arg)
       age2  = exp(arg2)
       age1  = exp(arg)
       dela  = 1.e-6*dtime*(age1+age2+age3)
       dels  = 0.1*max(0.0_r8, h2osno-h2osno_old)
       sge   = (snow_age+dela)*(1.0-dels)
       snow_age   = max(0.0_r8,sge)

    endif

  end subroutine SnowAge

end module SurfaceAlbedoMod
