#include <misc.h>
#include <preproc.h>

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: iniTimeVar
!
! !INTERFACE:
subroutine iniTimeVar (readini, eccen, obliqr, lambm0 , mvelpp)
!
! !DESCRIPTION: 
! Initializes the following time varying variables:
! water      : h2osno, h2ocan, h2osoi_liq, h2osoi_ice, h2osoi_vol
! snow       : snowdp, snowage, snl, dz, z, zi 
! temperature: t_soisno, t_veg, t_grnd
! The variable, h2osoi_vol, is needed by clm_soilalb -this is not needed on 
! restart since it is computed before the soil albedo computation is called.
! The remaining variables are initialized by calls to ecosystem dynamics
! and albedo subroutines. 
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  use clmpoint
  use globals
  use clm_varpar, only : nlevsoi, nlevsno, nlevlak
  use clm_varcon, only : bdsno, istice, istwet, istsoil, denice, denh2o, spval, zlnd 
  use shr_const_mod, only : SHR_CONST_TKFRZ
  use inicFileMod, only : inicrd
  use time_manager, only : get_nstep, get_curr_calday    
  use spmdMod, only : masterproc
  use EcosystemDynMod, only : EcosystemDyn
  use FracWetMod, only : FracWet
  use SurfaceAlbedoMod, only : SurfaceAlbedo
!
! !ARGUMENTS:
  implicit none
  logical , intent(in) :: readini  !true if read in initial data set
  real(r8), intent(in) :: eccen    !Earth's orbital eccentricity
  real(r8), intent(in) :: obliqr   !Earth's obliquity in radians
  real(r8), intent(in) :: lambm0   !Mean longitude of perihelion at the vernal equinox (radians)
  real(r8), intent(in) :: mvelpp   !Earth's moving vernal equinox long. of perihelion + pi (radians)
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
  character(len=16) :: initype               !type of initial dataset
  integer  :: j,ci,pi                        !loop indices
  logical  :: doalb                          !true => albedo time step
  real(r8) :: calday                         !calendar day
  integer  :: ier                            !MPI return code
  type(gridcell_type)       , pointer :: g   !local pointer to derived subtype
  type(landunit_type)       , pointer :: l   !local pointer to derived subtype
  type(column_type)         , pointer :: c   !local pointer to derived subtype
  type(pft_type)            , pointer :: p   !local pointer to derived subtype
  type(gridcell_pstate_type), pointer :: gps !local pointer to derived subtype
  type(landunit_pstate_type), pointer :: lps !local pointer to derived subtype
  type(column_pstate_type)  , pointer :: cps !local pointer to derived subtype	
  type(column_estate_type)  , pointer :: ces !local pointer to derived subtype	
  type(column_wstate_type)  , pointer :: cws !local pointer to derived subtype	
  type(pft_pstate_type)     , pointer :: pps !local pointer to derived subtype	
  type(pft_wstate_type)     , pointer :: pws !local pointer to derived subtype	
  type(pft_estate_type)     , pointer :: pes !local pointer to derived subtype	
!-----------------------------------------------------------------------

  ! Initialize water and temperature based on:
  ! readini = true : read initial data set -- requires netCDF codes
  ! readini = false: arbitrary initialization

  if (readini) then

     if ( masterproc ) write (6,*) 'Reading initial data '

     ! Open and read data from netCDF initial dataset

     call inicrd () 

     ! Determine volumetric soil water

     do ci = cols1d%beg,cols1d%end
        c => cpoint%col(ci)%c
        cws => c%cws
        cps => c%cps
        do j = 1,nlevsoi
           cws%h2osoi_vol(j) = cws%h2osoi_liq(j)/(cps%dz(j)*denh2o) &
                + cws%h2osoi_ice(j)/(cps%dz(j)*denice)
        end do
     end do

  else

     if ( masterproc ) write (6,*) 'Setting initial data to non-spun up values'

     ! ========================================================================
     ! Set snow water 
     ! ========================================================================

     ! NOTE: h2ocan, h2osno, snowdp and snowage has valid values everywhere

     do ci = cols1d%beg,cols1d%end
        c => cpoint%col(ci)%c
        cws => c%cws
        cps => c%cps
        lps => cps%lps

        ! canopy water (pft level and column averaged)
        do pi = 1,cps%npfts 
           p => c%p(pi)
           pws => p%pws
           pws%h2ocan = 0._r8
        end do
        cws%pws_a%h2ocan = 0._r8

        ! snow water
        if (lps%itypwat == istice) then
           cws%h2osno = 1000._r8
        else
           cws%h2osno = 0._r8
        endif

        ! snow depth
        cps%snowdp  = cws%h2osno/bdsno

        ! snow age
        cps%snowage = 0.                
     end do

     ! ========================================================================
     ! Set snow layer number, depth and thickiness 
     ! ========================================================================

     call snowdp2lev ()

     ! ========================================================================
     ! Set snow/soil temperature
     ! ========================================================================

     ! NOTE: 
     ! t_soisno only has valid values over non-lake
     ! t_lake   only has valid values over lake
     ! t_grnd has valid values over all land
     ! t_veg  has valid values over all land

     do ci = cols1d%beg,cols1d%end
        c => cpoint%col(ci)%c
        cps => c%cps
        ces => c%ces
        lps => cps%lps

        ces%t_soisno(-nlevsno+1:nlevsoi) = spval
        ces%t_lake(1:nlevlak) = spval
        if (.not. lps%lakpoi) then  !not lake
           ces%t_soisno(-nlevsno+1:0) = spval
           if (cps%snl < 0) then    !snow layer temperatures
              do j = cps%snl+1, 0
                 ces%t_soisno(j) = 250._r8
              enddo
           endif
           if (lps%itypwat == istice) then
              do j = 1, nlevsoi
                 ces%t_soisno(j) = 250._r8
              end do
           else if (lps%itypwat == istwet) then
              do j = 1, nlevsoi
                 ces%t_soisno(j) = 277._r8
              end do
           else
              do j = 1, nlevsoi
                 ces%t_soisno(j) = 283._r8
              end do
           endif
           ces%t_grnd = ces%t_soisno(cps%snl+1)
        else                           !lake
           ces%t_lake(1:nlevlak) = 277._r8
           ces%t_grnd = ces%t_lake(1)
        endif
        do pi = 1,cps%npfts 
           p => c%p(pi)
           pes => p%pes
           pes%t_veg = 283._r8
        end do   ! end of pi loop 
     end do   !end of ci loop

     ! ========================================================================
     ! Set snow/soil ice and liquid mass
     ! ========================================================================

     ! volumetric water is set first and liquid content and ice lens are obtained
     ! NOTE: h2osoi_vol, h2osoi_liq and h2osoi_ice only have valid values over soil

     do ci = cols1d%beg,cols1d%end
        c => cpoint%col(ci)%c
        cws => c%cws
        cps => c%cps
        ces => c%ces
        lps => cps%lps

        cws%h2osoi_vol(         1:nlevsoi) = spval
        cws%h2osoi_liq(-nlevsno+1:nlevsoi) = spval
        cws%h2osoi_ice(-nlevsno+1:nlevsoi) = spval

        if (.not. lps%lakpoi) then  !not lake
           ! volumetric water
           do j = 1,nlevsoi
              if (lps%itypwat == istsoil) then 
                 cws%h2osoi_vol(j) = 0.3_r8
              else
                 cws%h2osoi_vol(j) = 1.0_r8
              endif
              cws%h2osoi_vol(j) = min(cws%h2osoi_vol(j),lps%watsat(j))
           end do

           ! liquid water and ice (note that dz is in meters below)
           if (cps%snl < 0) then    !snow 
              do j = cps%snl+1, 0
                 cws%h2osoi_ice(j) = cps%dz(j)*250.
                 cws%h2osoi_liq(j) = 0._r8
              enddo
           endif
           do j = 1, nlevsoi           !soil layers
              if (ces%t_soisno(j) <= SHR_CONST_TKFRZ) then
                 cws%h2osoi_ice(j)  = cps%dz(j)*denice*cws%h2osoi_vol(j) 
                 cws%h2osoi_liq(j) = 0._r8
              else
                 cws%h2osoi_ice(j) = 0._r8
                 cws%h2osoi_liq(j) = cps%dz(j)*denh2o*cws%h2osoi_vol(j)
              endif
           enddo
        endif
     end do    !end of ci loop

  end if  ! end of arbitrary initialization if-block

  ! ========================================================================
  ! Remaining variables are initialized by calls to ecosystem dynamics and
  ! albedo subroutines. 
  ! Note: elai, esai, frac_veg_nosno_alb are computed in Ecosysdyn and needed
  ! by Fwet and SurfaceAlbedo
  ! frac_veg_nosno is needed by fwet
  ! fwet is needed in routine clm_twostream (called by SurfaceAlbedo)
  ! frac_sno is needed by SoilAlbedo (called by SurfaceAlbedo)
  ! ========================================================================

  nstep = get_nstep()
  calday = get_curr_calday()
  doalb = .true.

  do ci = cols1d%beg,cols1d%end
     c => cpoint%col(ci)%c
     cps => c%cps

     ! Determine variaglesneeded for albedo subroutines
     call EcosystemDyn(c, doalb, .false.)	

     ! Compute frac_sno
     cps%frac_sno = cps%snowdp/(10.*zlnd + cps%snowdp)  

     ! Compute fwet
     cps%pps_a%fwet = 0.0
     do pi=1, cps%npfts
        p => c%p(pi)
        pps => p%pps
        pps%frac_veg_nosno = pps%frac_veg_nosno_alb  
        call FracWet(p)
        cps%pps_a%fwet = cps%pps_a%fwet + (pps%fwet * pps%wt)
     end do

     call SurfaceAlbedo (c, calday, eccen, obliqr, lambm0, mvelpp)
  end do

  return
end subroutine iniTimeVar
