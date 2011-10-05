#include <misc.h>
#include <preproc.h>

module accFldsMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: accFldsMod
!
! !DESCRIPTION:
! This module contains subroutines that initialize, update and extract
! the user-specified fields over user-defined intervals. Each interval
! and accumulation type is unique to each field processed.
! Subroutine [initAccumFlds] defines the fields to be processed
! and the type of accumulation. Subroutine [updateAccumFlds] does
! the actual accumulation for a given field. Fields are accumulated
! by calls to subroutine [update_accum_field]. To accumulate a field,
! it must first be defined in subroutine [initAccumFlds] and then
! accumulated by calls to [updateAccumFlds].
! Four types of accumulations are possible:
!   o average over time interval
!   o running mean over time interval
!   o running accumulation over time interval
! Time average fields are only valid at the end of the averaging interval.
! Running means are valid once the length of the simulation exceeds the
! averaging interval. Accumulated fields are continuously accumulated.
! The trigger value "-99999." resets the accumulation to zero.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use abortutils,   only: endrun
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: initAccFlds     ! Initialization accumulator fields
  public :: initAccClmtype  ! Initialize clmtype variables obtained from accum fields
  public :: updateAccFlds   ! Update accumulator fields
!
! !REVISION HISTORY:
! Created by M. Vertenstein 03/2003
!
!EOP

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initAccFlds()
!
! !INTERFACE:
  subroutine initAccFlds()
!
! !DESCRIPTION:
! Initializes accumulator and sets up array of accumulated fields
!
! !USES:
    use accumulMod   , only : init_accum_field, print_accum_fields
    use clm_time_manager , only : get_step_size
    use shr_const_mod, only : SHR_CONST_CDAY, SHR_CONST_TKFRZ
    use nanMod       , only : bigint
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY::
! Created by M. Vertenstein 03/2003
!
!EOP
!
! LOCAL VARIABLES:
!
    integer :: dtime                     !time step size
    integer, parameter :: not_used = bigint
!------------------------------------------------------------------------

    ! Hourly average of 2m temperature.

    dtime = get_step_size()
    call init_accum_field(name='TREFAV', units='K', &
         desc='average over an hour of 2-m temperature', &
         accum_type='timeavg', accum_period=nint(3600._r8/dtime), &
         subgrid_type='pft', numlev=1, init_value=0._r8)

! Running monthly/daily means and hourly outputs for TEM
		 
#if (defined COUP_TEM)

    call init_accum_field(name='EPOT30', units='mm/s', &
         desc='Daily mean Potential Evapotranspiration', &
         accum_type='timeavg', accum_period=-1, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field(name='EVAP30', units='mm/s', &
         desc='Daily mean of Evapotranspiration', &
         accum_type='timeavg', accum_period=-1, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field(name='DRAI30', units='mm/s', &
         desc='Daily mean of Drainage', &
         accum_type='timeavg', accum_period=-1, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field(name='SURF30', units='mm/s', &
         desc='Daily mean of Surface Runoff', &
         accum_type='timeavg', accum_period=-1, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field(name='SWE30', units='mm', &
         desc='Daily Snow Water Equivalent', &
         accum_type='timeavg', accum_period=-1, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field(name='S1M30', units='m', &
         desc='Daily mean of 1 meter soil-water', &
         accum_type='timeavg', accum_period=-1, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field(name='S2M30', units='m', &
         desc='Daily mean of 2 meter soil-water', &
         accum_type='timeavg', accum_period=-1, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field (name='TSOIDY', units='K', &
         desc='Daily mean of soil temperature', &
         accum_type='timeavg', accum_period=-1, &
         subgrid_type='pft', numlev=10,init_value=SHR_CONST_TKFRZ+20._r8, &
         type2d='levsoi')

    call init_accum_field (name='HSOIDY', units='K', &
         desc='Daily mean of soil temperature', &
         accum_type='timeavg', accum_period=-1, &
         subgrid_type='pft', numlev=10,init_value=0._r8, &
         type2d='levsoi')

    call init_accum_field (name='PRCPDY', units='mm/s', &
         desc='Daily mean of precipitation', &
         accum_type='timeavg', accum_period=-1, &
         subgrid_type='pft', numlev=1,init_value=0._r8)

    call init_accum_field (name='STRMDY', units='s', &
         desc='Daily mean of Storm Duration', &
         accum_type='timeavg', accum_period=-1, &
         subgrid_type='pft', numlev=1,init_value=0._r8)

    !call init_accum_field(name='HSOIHR', units='m', &
    !     desc='Hourly average of soil water', &
    !     accum_type='timeavg', accum_period=nint(3600._r8/dtime), &
    !     subgrid_type='pft', numlev=10, init_value=0._r8, &
    !     type2d='levsoi')

#endif

#if (defined DGVM)
    ! 30-day average of 2m temperature.

    call init_accum_field (name='TDA', units='K', &
         desc='30-day average of 2-m temperature', &
         accum_type='timeavg', accum_period=-30, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! The following are running means.
    ! The accumulation period is set to 10 days for a 10-day running mean.

    call init_accum_field (name='T10', units='K', &
         desc='10-day running mean of 2-m temperature', &
         accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1,init_value=SHR_CONST_TKFRZ+20._r8)

    call init_accum_field (name='FNPSN10', units='UMOL/M2S', &
         desc='10-day running mean net cpy photosynth', &
         accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field (name='PREC365', units='MM H2O/S', &
         desc='365-day running mean of total precipitation', &
         accum_type='runmean', accum_period=-365, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! The following are accumulated fields.
    ! These types of fields are accumulated until a trigger value resets
    ! the accumulation to zero (see subroutine update_accum_field).
    ! Hence, [accper] is not valid.

    call init_accum_field (name='AGDD0', units='K', &
         desc='growing degree-days base 0C', &
         accum_type='runaccum', accum_period=not_used, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field (name='AGDD5', units='K', &
         desc='growing degree-days base -5C', &
         accum_type='runaccum', accum_period=not_used, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field (name='AGDDTW', units='K', &
         desc='growing degree-days base twmax', &
         accum_type='runaccum', accum_period=not_used, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field (name='AGDD', units='K', &
         desc='growing degree-days base 5C', &
         accum_type='runaccum', accum_period=not_used,  &
         subgrid_type='pft', numlev=1, init_value=0._r8)
#endif

    ! Print output of accumulated fields

    call print_accum_fields()

  end subroutine initAccFlds

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: updateAccFlds
!
! !INTERFACE:
  subroutine updateAccFlds()
!
! !DESCRIPTION:
! Update and/or extract accumulated fields
!
! !USES:
    use clmtype
    use clm_atmlnd   , only : clm_a2l, atm_a2l
    use decompMod    , only : get_proc_bounds
    use clm_varcon   , only : spval, istwet, istsoil, istdlak, isturb, istice
    use pftvarcon    , only : pftpar
    use shr_const_mod, only : SHR_CONST_CDAY, SHR_CONST_TKFRZ
    use clm_time_manager , only : get_step_size, get_nstep, is_end_curr_day, get_curr_date, get_prev_date
    use accumulMod   , only : update_accum_field, extract_accum_field
    use clm_varpar, only : nlevsoi

!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by M. Vertenstein 03/2003
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: ptype(:)            ! pft vegetation
    integer , pointer :: ltype(:)            ! land type for grid
    integer , pointer :: clandunit(:)        ! landunit indice for column
    integer , pointer :: plandunit(:)        ! landunit indice for pft
    integer , pointer :: pgridcell(:)        ! index into gridcell level quantities
    integer , pointer :: cgridcell(:)        ! index into gridcell level quantities
    integer , pointer :: pcolumn(:)          ! index into column level quantities
    real(r8), pointer :: forc_t(:)           ! atmospheric temperature (Kelvin)
    real(r8), pointer :: forc_rain(:)        ! rain rate [mm/s]
    real(r8), pointer :: forc_snow(:)        ! snow rate [mm/s]
    real(r8), pointer :: forc_strm_dur(:)    ! storm duration (s) 
    real(r8), pointer :: frmf(:)             ! leaf maintenance respiration  (umol CO2 /m**2 /s)
    real(r8), pointer :: fpsn(:)             ! photosynthesis (umol CO2 /m**2 /s)
    real(r8), pointer :: t_ref2m(:)          ! 2 m height surface air temperature (Kelvin)
    real(r8), pointer :: qflx_evap_tot(:)    ! Total Evapotranspiration (mm/s)
    real(r8), pointer :: qflx_surf(:)        ! Surface Runoff (mm/s)
    real(r8), pointer :: qflx_drain(:)       ! Drainage (mm/s)
    real(r8), pointer :: h2osno(:)           ! Snow Water Equivalent (mm)
    real(r8), pointer :: h2osoi_liq(:,:)     ! Liquid Soil Water (mm)
    real(r8), pointer :: h2osoi_ice(:,:)     ! Frozen Soil Water (mm)
    real(r8), pointer :: t_soisno(:,:)       ! Soil Temperature (K)
	
#if (defined COUP_TEM)
include 'TEM.inc'
    real(r8), pointer :: qflx_evap_pet1(:)   ! Potential Evapotranspiration (mm/s)
    real(r8), pointer :: qflx_evap_pet2(:)   ! Potential Evapotranspiration (mm/s)
    real(r8), pointer :: s1msat(:)           ! Saturated 1 meter soil moisture (mm)
    real(r8), pointer :: s2msat(:)           ! Saturated 2 meter soil moisture (mm)
    real(r8), pointer :: epot30r(:,:)        ! Array for 30-day average potential evap. (mm/s)
    real(r8), pointer :: evap30r(:,:)        ! Array for 30-day average evapotranspiration (mm/s)
    real(r8), pointer :: surf30r(:,:)        ! Array for 30-day average surface runoff (mm/s)
    real(r8), pointer :: drai30r(:,:)        ! Array for 30-day average drainage (mm/s)
    real(r8), pointer :: swe30r(:,:)         ! Array for 30-day average snow water-equivalent depth (mm)
    real(r8), pointer :: s1m30r(:,:)         ! Array for 30-day average 1 meter soil moisture (mm)
    real(r8), pointer :: s2m30r(:,:)         ! Array for 30-day average 2 meter soil moisture (mm)
    real(r8), pointer :: evap30(:)           ! 30-day average evapotranspiration (mm/s)
    real(r8), pointer :: epot30(:)           ! 30-day average potential evap. (mm/s)
    real(r8), pointer :: surf30(:)           ! 30-day average surface runoff (mm/s)
    real(r8), pointer :: drai30(:)           ! 30-day average drainage (mm/s)
    real(r8), pointer :: swe30(:)            ! 30-day average snow water-equivalent depth (mm)
    real(r8), pointer :: s1m30(:)            ! 30-day average 1 meter soil moisture (mm)
    real(r8), pointer :: s2m30(:)            ! 30-day average 2 meter soil moisture (mm)
    real(r8), pointer :: dyh_soi(:,:)        ! Daily mean soil-water profile 
    real(r8), pointer :: dyt_soi(:,:)        ! Daily mean soil-temperature profile (K)
    real(r8), pointer :: hrh_soi(:,:)        ! Hourly soil-water profile
    real(r8), pointer :: watsat(:,:)         ! Porosity
    real(r8), pointer :: dz(:,:)             ! Soil Layer Thickness
    real(r8), pointer :: dypcp(:)            ! Daily mean preciptiation (mm/s)
    real(r8), pointer :: dystrm(:)           ! Daily storm duration (s)
    real(r8), pointer :: dystrm_inst(:)      ! daily max of storm duration
    real(r8), pointer :: wtgcell(:)          ! weight (relative to land/gridcell)
    real(r8) :: wetlndf(mxmsaics,igsmlat)    !Wetland type fractional area of zonal band
    data ((wetlndf(i,j), i=1,17), j=1,46) / 782*0.0 /
    data ((wetlndf(i,j), i=31,35), j=1,46) / 230*0.0 /
    data (wetlndf(18,i), i=1,46) / 15*0.0, &
          0.11319,0.520604,0.703485,0.351305,0.402655,0.700123,0.585064,0.489737,0.722638, &
          0.385134,0.171036,0.253563,0.866345,0.820354,0.0989037,0.008995,15*0.0 /
    data (wetlndf(19,i), i=1,46) / 15*0.0, &
	  0.0292921,0.170999,0.195895,0.194403,0.36631,0.0982119,0.104799,0.111005,0.0993542, &
	  0.326655,0.7263,0.639274,0.133655,0.166269,0.441134,0.0316716, 15*0.0 /
    data (wetlndf(20,i), i=1,46) / 9*0.0,5*0.592545,0.418213,0.129242,0.0,0.0217687,3*0.0, &
          0.00111177,0.00819991,3*0.0,0.0288033,0.0,0.0133779,0.0747448,0.189465,0.145859, &
          0.238796,0.0887471,0.112111,0.00375326,0.0120095, 9*0.0 /
    data (wetlndf(21,i), i=1,46) / 9*0.0,5*0.407455,0.261997,0.230551,0.308398,0.03036,0.386048,4*0.0, &
          0.0405736,0.0,0.0,0.0050863,0.0,0.0,0.00609636,0.187267,0.28152,0.04693,0.87239,0.274892,0.213517, &
	  0.0307741,0.0,0.0431368,0.0149039,6*0.0 /
    data (wetlndf(22,i), i=1,46) / 22*0.0,0.219633,0.0604931,5*0.0,0.218677,0.231012,0.3662,0.360184,0.0388628, &
	  0.587504,0.657288,0.687025,0.606983,0.623981,0.718816,0.0717533,0.0717533,4*0.0 /
    data (wetlndf(23,i), i=1,46) / 35*0.0,0.101793,0.270191,0.393017,0.332882,0.26628,0.928247,0.928247,4*0.0 /
    data (wetlndf(24,i), i=1,46) / 17*0.0,0.0484911,0.0682437,0.0,0.0,0.0709459,0.0758356,0.0769407, &
          0.24811,0.0420167,0.0732732,0.0,0.0,0.0735224,16*0.0 /
    data (wetlndf(25,i), i=1,46) / 29*0.0,0.0869223,0.330944,0.0909473,0.354089,0.0,0.0254925,11*0.0 /
    data (wetlndf(26,i), i=1,46) / 15*0.0,0.341153,19*0.0,0.0236483,10*0.0 /
    data (wetlndf(27,i), i=1,46) / 15*0.0,0.0469222,3*0.0,0.129817,0.0997384,0.188806,0.0538077,0.0, &
          0.0260633,0.0433985,20*0.0 /
    data (wetlndf(28,i), i=1,46) / 15*0.0,0.109649,3*0.0,0.101218,0.101927,0.0492723,0.0,0.0,0.0140386, &
          5*0.0,0.00393006,15*0.0 /
    data (wetlndf(29,i), i=1,46) / 14*0.0,0.159895,15*0.0,0.0103222,0.115474,14*0.0 /
    data (wetlndf(30,i), i=1,46) / 14*0.0,0.159895,7*0.0,0.0417816,2*0.0,0.0172484,4*0.0,0.00639219,15*0.0 /

#endif

!
! local pointers to implicit out arguments
!
    real(r8), pointer :: t_ref2m_min(:)      ! daily minimum of average 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_max(:)      ! daily maximum of average 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_min_inst(:) ! instantaneous daily min of average 2 m height surface air temp (K)
    real(r8), pointer :: t_ref2m_max_inst(:) ! instantaneous daily max of average 2 m height surface air temp (K)
    real(r8), pointer :: t10(:)              ! 10-day running mean of the 2 m temperature (K)
    real(r8), pointer :: t_mo(:)             ! 30-day average temperature (Kelvin)
    real(r8), pointer :: t_mo_min(:)         ! annual min of t_mo (Kelvin)
    real(r8), pointer :: fnpsn10(:)          ! 10-day running mean net photosynthesis
    real(r8), pointer :: prec365(:)          ! 365-day running mean of tot. precipitation
    real(r8), pointer :: agdd0(:)            ! accumulated growing degree days above 0 deg C
    real(r8), pointer :: agdd5(:)            ! accumulated growing degree days above -5
    real(r8), pointer :: agddtw(:)           ! accumulated growing degree days above twmax
    real(r8), pointer :: agdd(:)             ! accumulated growing degree days above 5
!
!EOP
!
! OTHER LOCAL VARIABLES:
    integer :: g,l,c,p,i,j,k,m,lev,hr    ! indices
    integer :: itypveg,iveg,ii           ! vegetation type indices
    integer :: itypcol                   ! landunit type
    integer :: dtime                     ! timestep size [seconds]
    integer :: nstep                     ! timestep number
    integer :: year,yearp                ! year (0, ...) for nstep
    integer :: month,monthp              ! month (1, ..., 12) for nstep
    integer :: day,dayp                  ! day of month (1, ..., 31) for nstep
    integer :: secs,secsp                ! seconds into current date for nstep
    logical :: end_cd                    ! temporary for is_end_curr_day() value
    integer :: ier                       ! error status
    integer :: begp, endp                ! per-proc beginning and ending pft indices
    integer :: begc, endc                ! per-proc beginning and ending column indices
    integer :: begl, endl                ! per-proc beginning and ending landunit indices
    integer :: begg, endg                ! per-proc gridcell ending gridcell indices
    real(r8), pointer :: rbufslp(:)      ! temporary single level - pft level
    real(r8), pointer :: rbufmlp(:,:)    ! temporary multiple level - pft level
    real(r8), pointer :: rbufslc(:)      ! temporary single level - column level
!------------------------------------------------------------------------

    ! Determine necessary indices

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! Assign local pointers to derived subtypes components (gridcell-level)

    forc_t => clm_a2l%forc_t
    forc_rain => clm_a2l%forc_rain
    forc_snow => clm_a2l%forc_snow

    ! Assign local pointers to derived subtypes components (pft-level)

    ptype            => clm3%g%l%c%p%itype
    ltype            => clm3%g%l%itype
    pgridcell        => clm3%g%l%c%p%gridcell
    cgridcell        => clm3%g%l%c%gridcell
    pcolumn          => clm3%g%l%c%p%column
    plandunit        => clm3%g%l%c%p%landunit
    clandunit        => clm3%g%l%c%landunit
    t_ref2m          => clm3%g%l%c%p%pes%t_ref2m
    fpsn             => clm3%g%l%c%p%pcf%fpsn
    frmf             => clm3%g%l%c%p%pcf%frmf
    t_ref2m_max_inst => clm3%g%l%c%p%pes%t_ref2m_max_inst
    t_ref2m_min_inst => clm3%g%l%c%p%pes%t_ref2m_min_inst
    t_ref2m_max      => clm3%g%l%c%p%pes%t_ref2m_max
    t_ref2m_min      => clm3%g%l%c%p%pes%t_ref2m_min
    t_mo             => clm3%g%l%c%p%pdgvs%t_mo
    t_mo_min         => clm3%g%l%c%p%pdgvs%t_mo_min
    t10              => clm3%g%l%c%p%pdgvs%t10
    fnpsn10          => clm3%g%l%c%p%pdgvs%fnpsn10
    prec365          => clm3%g%l%c%p%pdgvs%prec365
    agdd0            => clm3%g%l%c%p%pdgvs%agdd0
    agdd5            => clm3%g%l%c%p%pdgvs%agdd5
    agddtw           => clm3%g%l%c%p%pdgvs%agddtw
    agdd             => clm3%g%l%c%p%pdgvs%agdd

#if (defined COUP_TEM)	
    epot30r          => clm3%g%l%c%p%pps%epot30r
    evap30r          => clm3%g%l%c%p%pps%evap30r
    drai30r          => clm3%g%l%c%p%pps%drai30r
    surf30r          => clm3%g%l%c%p%pps%surf30r
    swe30r           => clm3%g%l%c%p%pps%swe30r
    s1m30r           => clm3%g%l%c%p%pps%s1m30r
    s2m30r           => clm3%g%l%c%p%pps%s2m30r
    epot30           => clm3%g%l%c%p%pps%epot30
    evap30           => clm3%g%l%c%p%pps%evap30
    drai30           => clm3%g%l%c%p%pps%drai30
    surf30           => clm3%g%l%c%p%pps%surf30
    swe30            => clm3%g%l%c%p%pps%swe30
    s1m30            => clm3%g%l%c%p%pps%s1m30
    s2m30            => clm3%g%l%c%p%pps%s2m30
    dyt_soi          => clm3%g%l%c%p%pps%dyt_soi
    dyh_soi          => clm3%g%l%c%p%pps%dyh_soi
    hrh_soi          => clm3%g%l%c%p%pps%hrh_soi
    watsat           => clm3%g%l%c%cps%watsat
    dz               => clm3%g%l%c%cps%dz
    dypcp            => clm3%g%l%c%p%pps%dypcp
    dystrm           => clm3%g%l%c%p%pps%dystrm
    dystrm_inst      => clm3%g%l%c%p%pps%dystrm_inst
    qflx_evap_pet1   => clm3%g%l%c%p%pwf%qflx_evap_pet1
    qflx_evap_pet2   => clm3%g%l%c%p%pwf%qflx_evap_pet2
    qflx_evap_tot    => clm3%g%l%c%p%pwf%qflx_evap_tot
    qflx_surf        => clm3%g%l%c%cwf%qflx_surf
    qflx_drain       => clm3%g%l%c%cwf%qflx_drain
    h2osno           => clm3%g%l%c%cws%h2osno
    h2osoi_liq       => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_ice       => clm3%g%l%c%cws%h2osoi_ice
    t_soisno         => clm3%g%l%c%ces%t_soisno
#if (defined STOCHASTIC)
    forc_strm_dur    => atm_a2l%forc_strm_dur
#endif
    wtgcell          => clm3%g%l%c%p%wtgcell
#endif

    ! Determine calendar information

    dtime = get_step_size()
    nstep = get_nstep()
    call get_curr_date (year, month, day, secs)

    ! Don't do any accumulation if nstep is zero
    ! (only applies to coupled or cam mode)

    if (nstep == 0) return

    ! NOTE: currently only single level pft fields are used below
    ! Variables are declared above that should make it easy to incorporate
    ! multi-level or single-level fields of any subgrid type

    ! Allocate needed dynamic memory for single level pft field

    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(6,*)'update_accum_hist allocation error for rbuf1dp'
       call endrun
    endif

    ! Accumulate and extract TREFAV - hourly average 2m air temperature
    ! Used to compute maximum and minimum of hourly averaged 2m reference
    ! temperature over a day. Note that "spval" is returned by the call to
    ! accext if the time step does not correspond to the end of an
    ! accumulation interval. First, initialize the necessary values for
    ! an initial run at the first time step the accumulator is called

    call update_accum_field  ('TREFAV', t_ref2m, nstep)
    call extract_accum_field ('TREFAV', rbufslp, nstep)

    end_cd = is_end_curr_day()
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (rbufslp(p) /= spval) then
          t_ref2m_max_inst(p) = max(rbufslp(p), t_ref2m_max_inst(p))
          t_ref2m_min_inst(p) = min(rbufslp(p), t_ref2m_min_inst(p))
       endif
       if (end_cd) then
          t_ref2m_max(p) = t_ref2m_max_inst(p)
          t_ref2m_min(p) = t_ref2m_min_inst(p)
          t_ref2m_max_inst(p) = -spval
          t_ref2m_min_inst(p) =  spval
       else if (secs == int(dtime)) then
          t_ref2m_max(p) = spval
          t_ref2m_min(p) = spval
       endif
    end do

!CAS
!CAS    Code augments for TEM/NEM coupling
!CAS
#if (defined COUP_TEM)
    hr = (secs/dtime) + 1
    allocate(rbufmlp(begp:endp,1:nlevsoi), stat=ier)
    allocate(rbufslc(begc:endc), stat=ier)
    allocate(s1msat(begc:endc))
    allocate(s2msat(begc:endc))

!CAS	NOTE: In situations where multiple storms would occur during the day, this 
!CAS			procedure will take the daily mean. 
!CAS
!CAS    For remainder of variables - use the accum/average routines accordingly
#if (defined STOCHASTIC)
    do p = begp,endp
       g = pgridcell(p)
       rbufslp(p) = forc_strm_dur(g)/float(dtime)
    end do
    call update_accum_field  ('STRMDY', rbufslp, nstep)
    call extract_accum_field ('STRMDY', dystrm, nstep)
#endif

!dir$ concurrent
!cdir nodep

!dir$ concurrent
!cdir nodep

    do p = begp,endp
       rbufslp(p) = qflx_evap_pet1(p)
    end do
    call update_accum_field  ('EPOT30', rbufslp, nstep)
    call extract_accum_field ('EPOT30', epot30, nstep)

    do p = begp,endp
       rbufslp(p) = qflx_evap_tot(p)
    end do
    call update_accum_field  ('EVAP30', rbufslp, nstep)
    call extract_accum_field ('EVAP30', evap30, nstep)

    do p = begp,endp
       c = pcolumn(p)
       rbufslp(p) = qflx_surf(c)
    end do
    call update_accum_field  ('SURF30', rbufslp, nstep)
    call extract_accum_field ('SURF30', surf30, nstep)

    do p = begp,endp
       c = pcolumn(p)
       rbufslp(p) = qflx_drain(c)
    end do
    call update_accum_field  ('DRAI30', rbufslp, nstep)
    call extract_accum_field ('DRAI30', drai30, nstep)

    do p = begp,endp
       c = pcolumn(p)
       rbufslp(p) = h2osno(c)
    end do
    call update_accum_field  ('SWE30', rbufslp, nstep)
    call extract_accum_field ('SWE30', swe30, nstep)

    do p = begp,endp
       rbufslp(p)=0.0
       c = pcolumn(p)
       s1msat(c)=0.0
       do lev=1,7
         rbufslp(p) = rbufslp(p)+h2osoi_liq(c,lev)
         s1msat(c) = s1msat(c)+(1000.0*watsat(c,lev)*dz(c,lev))
       enddo
!CAS
!CAS    Scale last layer to estimate soil water in top 1 meter
!CAS
       s1msat(c) = s1msat(c)+((0.1711/0.5539)*1000.0*watsat(c,lev)*dz(c,lev))
       rbufslp(p) = rbufslp(p)+((0.1711/0.5539)*h2osoi_liq(c,8))
    end do
    call update_accum_field  ('S1M30', rbufslp, nstep)
    call extract_accum_field ('S1M30', s1m30, nstep)

    do p = begp,endp
       rbufslp(p)=0.0
       c = pcolumn(p)
       s2msat(c)=0.0
       do lev=1,8
         rbufslp(p) = rbufslp(p)+h2osoi_liq(c,lev)
         s2msat(c) = s2msat(c)+(1000.0*watsat(c,lev)*dz(c,lev))
       enddo
!CAS
!CAS    Scale last layer to estimate soil water in top 2 meter
!CAS
       s2msat(c) = s2msat(c)+((0.61716/0.91329)*1000.0*watsat(c,lev)*dz(c,lev))
       rbufslp(p) = rbufslp(p)+((0.61716/0.91329)*h2osoi_liq(c,9))
    end do
    call update_accum_field  ('S2M30', rbufslp, nstep)
    call extract_accum_field ('S2M30', s2m30, nstep)
	
! CAS
! CAS	Check that daily timeaveraging is completed for 30-day running mean
! CAS	If value of daily 'timeavg' for evap is /= bogus value
! CAS	then assume all variables are ready for running averaging
! CAS

    m = (endp - begp)/2
    if (evap30(m) /= spval) then
      do i=30,2,-1
         epot30r(i,:) = epot30r(i-1,:)
         evap30r(i,:) = evap30r(i-1,:)
         surf30r(i,:) = surf30r(i-1,:)
         drai30r(i,:) = drai30r(i-1,:)
         swe30r(i,:) = swe30r(i-1,:)
         s1m30r(i,:) = s1m30r(i-1,:)
         s2m30r(i,:) = s2m30r(i-1,:)
      enddo
      epot30r(1,:) = epot30
      evap30r(1,:) = evap30
      surf30r(1,:) = surf30
      drai30r(1,:) = drai30
      swe30r(1,:) = swe30
      s1m30r(1,:) = s1m30
      s2m30r(1,:) = s2m30
    endif

!    call get_prev_date (yearp, monthp, dayp, secsp)
!    if ( day /= dayp ) then
!      write (6,*) epot30
!    endif

!CAS Converting to volumetric soil-water for daily and hourly profiles

    do p = begp,endp
       c = pcolumn(p)
       do lev=1,nlevsoi
         rbufmlp(p,lev) = (h2osoi_liq(c,lev)+h2osoi_ice(c,lev))/(1000.0*dz(c,lev))
         hrh_soi(p,lev) = rbufmlp(p,lev)
       enddo
    end do
    call update_accum_field  ('HSOIDY', rbufmlp, nstep)
    call extract_accum_field ('HSOIDY', dyh_soi, nstep)

    do p = begp,endp
      c = pcolumn(p)
      do lev=1,nlevsoi
        rbufmlp(p,lev) = t_soisno(c,lev)
      enddo
    end do
    call update_accum_field  ('TSOIDY', rbufmlp, nstep)
    call extract_accum_field ('TSOIDY', dyt_soi, nstep)

    do p = begp,endp
       g = pgridcell(p)
       rbufslp(p) = forc_rain(g) + forc_snow(g)
    end do
    call update_accum_field  ('PRCPDY', rbufslp, nstep)
    call extract_accum_field ('PRCPDY', dypcp, nstep)

    mitpet = 0.0
    mitaet = 0.0
    mitsfr = 0.0
    mitdrn = 0.0
    mitswe = 0.0
    mitsh2o1m = 0.0
    mitsh2o2m = 0.0
    do p = begp,endp
       g = pgridcell(p)
       c = pcolumn(p)
       l = plandunit(p)
       itypcol = ltype(l)
       if(itypcol==2) then   ! LAND ICE CASE
         iveg=31
       elseif(itypcol==3.or.itypcol==4) then  ! DEEP/SHALLOW LAKE CASE
         iveg=32
       elseif(itypcol==6) then  ! URBAN CASE
         iveg=35
       else
         iveg = ptype(p)+1
       endif
!
!	 Adjust/assign grid indice according to IGSM zonal land bands
!        Note:  This procedure also accounts for IGSMVEG types having a 
!               spectrum of wetland types - yet CLM makes no distinction at this
!               time.
!
       if (g > 7 ) g = g + 2
!
!CAS
!CAS Fill Vegetated subgrid elements
!CAS
       if (ltype(l)==istsoil) then
         do i=1,30
           mitpet(iveg,g) = mitpet(iveg,g) + (86400*epot30r(i,p)/30.)
           mitaet(iveg,g) = mitaet(iveg,g) + (86400*evap30r(i,p)/30.)
           mitsfr(iveg,g) = mitsfr(iveg,g) + (86400*surf30r(i,p)/30.)
           mitdrn(iveg,g) = mitdrn(iveg,g) + (86400*drai30r(i,p)/30.)
           mitsh2o1m(iveg,g) = mitsh2o1m(iveg,g) + (s1m30r(i,p)/30.)
           mitsh2o2m(iveg,g) = mitsh2o2m(iveg,g) + (s2m30r(i,p)/30.)
    	   mitswe(iveg,g) = mitswe(iveg,g)+(swe30r(i,p)/30.)
!CAS
!CAS    Fill in Wetland values (NOTE OFFSET WITH INDICE = IGSMVEG+1)
!CAS    NOTE: IN THIS CONSTRUCTION WETLANDS NOT ASSUMED TO BE SATURATED
!CAS
           mitsh2o1m(18:26,g) = s1msat(c)
           mitsh2o2m(18:26,g) = s2msat(c)
    	   mitswe(18:26,g) = mitswe(iveg,g)
         enddo
!CAS
!CAS    Fill-in vegetation area fraction
!CAS

         vegfrac(iveg,g) = wtgcell(p)

!CAS
!CAS    Similar to check for monthly data
!CAS    if daily value for t_soi is not bogus value, then all daily data 
!CAS    is ready to be written to TEM/NEM common block
!CAS
         if (dyt_soi(p,1) /= spval) then
           do i = 1,17; mitdaytsoil(day,i,g,:) = dyt_soi(p,:); enddo
           do i = 27,35; mitdaytsoil(day,i,g,:) = dyt_soi(p,:); enddo
         endif
         if (dyh_soi(p,1) /= spval) then
           do i = 1,17; mitdaysh2o(day,i,g,:) = dyh_soi(p,:); enddo
           do i = 27,35; mitdaysh2o(day,i,g,:) = dyh_soi(p,:); enddo
         endif
!
!      Note: Calculations in this loop assume that timestep is an hour!
!
         do k = 1,6
          if (hrh_soi(p,k) /= spval) then
            do i = 1,17; mithrsh2o(hr,day,i,g,k) = hrh_soi(p,k); enddo
            do i = 27,35; mithrsh2o(hr,day,i,g,k) = hrh_soi(p,k); enddo
          endif
         enddo
!CAS
!CAS    Fill Wetland grid sub elements
!CAS
       elseif (ltype(l) == istwet) then
         if (dyt_soi(p,1) /= spval) then
           do i = 18,26; mitdaytsoil(day,i,g,:) = dyt_soi(p,:); enddo
         endif
         if (dyh_soi(p,1) /= spval) then
           do i = 18,26; mitdaysh2o(day,i,g,:) = dyh_soi(p,:); enddo
         endif
!
!      Note: Calculations in this loop assume that timestep is an hour!
!
         do k = 1,6
          if (hrh_soi(p,k) /= spval) then
            do i = 18,26; mithrsh2o(hr,day,i,g,k) = hrh_soi(p,k); enddo
          endif
         enddo
!
!	Note: Calculations in this loop assume that timestep is an hour!
!
#if (defined STOCHASTIC)
         if (dypcp(p) /= spval) mitqstrm(day,g) = dypcp(p)*dtime/10.0
         if (dystrm(p) /= spval .and. dystrm(p) /= 0.0) mitstrmdur(day,g) = dystrm(p)
#endif

!CAS
!CAS    Given that CLMVEG only consisders one type of Wetland
!CAS    Fill-in all IGSMVEG wetland varieties (IGSMVEGs 18 thru 30)
!CAS    and wetland subtypes partitioned according to 0.5x0.5 degree 
!CAS    IGSMVEG fractions (wetlndf, prescribed in DATA statement above)
!CAS
!CAS:AR5 - THIS IS NOW DONE IN CLIMATE2TEM TO ENSURE CONSISTENCY WITH TEM WETLAND
!CAS:AR5   INSTEAD - ALL WETLAND AREA IS PUT INTO VEGFRAC(18,*)
          vegfrac(18,g) = wtgcell(p)
!         do i = 18,30
!           vegfrac(i,g) = wetlndf(i,g)*wtgcell(p)
!         enddo
!CAS
!CAS    Fill land-ice grid sub elements
!CAS
       elseif (ltype(l) == istice) then
         vegfrac(31,g) = wtgcell(p)
!CAS
!CAS    Fill Lake grid sub elements
!CAS
       elseif (ltype(l) == istdlak) then
         vegfrac(32,g) = wtgcell(p)
!CAS
!CAS    Fill urban grid sub elements
!CAS
       elseif (ltype(l) == isturb) then
         vegfrac(35,g) = wtgcell(p)
       endif
    enddo
!CAS
!CAS	check that all vegfrac latitudes sum to 1
!CAS
!    write (6,*) 'VEGFRAC by Latitude:'
!    do l = 1,igsmlat
!      write (6,*) l,sum(vegfrac(1:35,l))
!    enddo
    do c = begc,endc
      l = clandunit(c)
      g = cgridcell(c)
      if (g > 7 ) g = g + 2
!        !write (6,*) 'WETLAND CASE @ grid:',g,'column',c
      do iveg = 18,26
!
!   To account for CLM wetland evaporation difficiency (uses low soil evaporation)
!   Set wetland evaporation to maximum evaporation from PFT at lat.
!
        mitaet(iveg,g) = 0.0
        mitpet(iveg,g) = 0.0
!        mitsfr(iveg,g) = 0.0
!        mitdrn(iveg,g) = 0.0
        do i = 1,17
          mitaet(iveg,g) = max(mitaet(i,g),mitaet(iveg,g))
          mitpet(iveg,g) = max(mitpet(i,g),mitpet(iveg,g))
!          mitsfr(iveg,g) = max(mitsfr(i,g),mitsfr(iveg,g))
!          mitdrn(iveg,g) = max(mitdrn(i,g),mitdrn(iveg,g))
        enddo
      enddo
    enddo

! CAS   Setting IGSMVEG=32 (rice paddies) according to CLMVEG=16 
! CAS   Setting IGSMVEG=33 (pastures) according to CLMVEG=16
! CAS   NOTE: ARRAY INDICE OFFSET (IGSMVEG+1 for mit*** arrays)
    do g = 1,igsmlat
      if (vegfrac(17,g) /= 0.0 ) then
        iveg = 17
        mitpet(33,g) = mitpet(iveg,g)
        mitaet(33,g) = mitaet(iveg,g)
        mitsfr(33,g) = mitsfr(iveg,g)
        mitdrn(33,g) = mitdrn(iveg,g)
        mitsh2o1m(33,g) = mitsh2o1m(iveg,g) 
        mitsh2o2m(33,g) = mitsh2o2m(iveg,g)
        !vegfrac(33,g) = vegfrac(iveg,g)

!CAS:AR5
!CAS:AR5  Null for AR5:  No rice paddies considered
!CAS:AR5
        vegfrac(33,g) = 0.0

        mitpet(34,g) = mitpet(iveg,g)
        mitaet(34,g) = mitaet(iveg,g) 
        mitsfr(34,g) = mitsfr(iveg,g)
        mitdrn(34,g) = mitdrn(iveg,g)
        mitsh2o1m(34,g) = mitsh2o1m(iveg,g) 
        mitsh2o2m(34,g) = mitsh2o2m(iveg,g)
        vegfrac(34,g) = vegfrac(iveg,g)

! CAS:AR5
! CAS:   Corresponding area is set to zero because IGSMVEG/TEM is using this array element
! CAS:   for fertilized crop area
! CAS:AR5
        vegfrac(iveg,g) = 0.0  
      endif
! CAS
! CAS   Populate floodplain values
! CAS   Setting IGSMVEG=26 (floodplains-tree-tropical) according to CLMVEG=4
! CAS   NOTE: ARRAY INDICE OFFSET (IGSMVEG+1 for mit*** arrays)
! CAS
      if (vegfrac(5,g) /= 0.0 ) then
        iveg = 5
      elseif (vegfrac(7,g) /= 0.0) then
        iveg = 7
      else
        iveg = 18
      endif
      mitpet(27,g) = mitpet(iveg,g) 
      mitaet(27,g) = mitaet(iveg,g) 
      mitsfr(27,g) = mitsfr(iveg,g)
      mitdrn(27,g) = mitdrn(iveg,g)
      mitsh2o1m(27,g) = mitsh2o1m(iveg,g) 
      mitsh2o2m(27,g) = mitsh2o2m(iveg,g)
! CAS
! CAS   Setting IGSMVEG=27 (floodplains-notree-tropical) according to CLMVEG=14
! CAS   If no cover for CLMVEG = 14 use CLMVEG = 13
! CAS   NOTE: ARRAY INDICE OFFSET (IGSMVEG+1 for mit*** arrays)
! CAS
      if (vegfrac(15,g) /= 0.0 ) then
        iveg = 15 
      elseif (vegfrac(14,g) /= 0.0) then
        iveg = 14
      else
        iveg = 18
      endif
      mitpet(28,g) = mitpet(iveg,g) 
      mitaet(28,g) = mitaet(iveg,g) 
      mitsfr(28,g) = mitsfr(iveg,g)
      mitdrn(28,g) = mitdrn(iveg,g)
      mitsh2o1m(28,g) = mitsh2o1m(iveg,g) 
      mitsh2o2m(28,g) = mitsh2o2m(iveg,g)
! CAS
! CAS   Setting IGSMVEG=28 (floodplains-tree-temperate) according to CLMVEG=7
! CAS   If no cover for CLMVEG = 7 use CLMVEG = 5
! CAS   NOTE: ARRAY INDICE OFFSET (IGSMVEG+1 for mit*** arrays)
! CAS
      if (vegfrac(8,g) /= 0.0 ) then
        iveg = 8
      elseif (vegfrac(6,g) /= 0.0) then
        iveg = 6
      else
        iveg = 18
      endif
      mitpet(29,g) = mitpet(iveg,g) 
      mitaet(29,g) = mitaet(iveg,g) 
      mitsfr(29,g) = mitsfr(iveg,g)
      mitdrn(29,g) = mitdrn(iveg,g)
      mitsh2o1m(29,g) = mitsh2o1m(iveg,g) 
      mitsh2o2m(29,g) = mitsh2o2m(iveg,g)
! CAS
! CAS   Setting IGSMVEG=29 (floodplains-notree-temperate) according to CLMVEG=13
! CAS   If no cover for CLMVEG = 13 use CLMVEG = 14
! CAS   NOTE: ARRAY INDICE OFFSET (IGSMVEG+1 for mit*** arrays)
! CAS
      if (vegfrac(14,g) /= 0.0 ) then
        iveg = 14 
      elseif (vegfrac(14,g) /= 0.0) then
        iveg = 15
      else
        iveg = 18
      endif
      mitpet(30,g) = mitpet(iveg,g) 
      mitaet(30,g) = mitaet(iveg,g) 
      mitsfr(30,g) = mitsfr(iveg,g)
      mitdrn(30,g) = mitdrn(iveg,g)
      mitsh2o1m(30,g) = mitsh2o1m(iveg,g) 
      mitsh2o2m(30,g) = mitsh2o2m(iveg,g)
    enddo !igsmlat
    where ( mitpet < mitaet ) mitpet=mitaet
!
!    call get_prev_date (yearp, monthp, dayp, secsp)
!    if ( month /= monthp ) then
!       write (6,*) 'Monthly Evapotranspiration Output to TEM:'
!       do j = 1,46
!         write(6,*) (mitaet(i,j),i=1,35)
!       enddo
!       write (6,*) 'Monthly Drainage to TEM:'
!       do j = 1,46
!         write(6,*) (mitdrn(i,j),i=1,35)
!       enddo
!       write (6,*) 'Monthly Potential Evapotranspiration Output to TEM:'
!       do j = 1,46
!         write(6,*) (mitpet(i,j),i=1,35)
!       enddo
!       write (6,*) 'Monthly SWE Output to TEM:'
!       do j = 1,46
!         write(6,*) (mitswe(i,j),i=1,35)
!       enddo
!       write (6,*) 'Monthly 1m Soil Moisture Output to TEM:'
!       do j = 1,46
!         write(6,*) (mitsh2o1m(i,j),i=1,35)
!       enddo
!       write (6,*) 'Monthly 2m Soil Moisture Output to TEM:'
!       do j = 1,46
!         write(6,*) (mitsh2o2m(i,j),i=1,35)
!       enddo
!       write (6,*) 'Daily Top Soil Layer Temperature:'
!       do j = 1,46
!         write(6,*) (mitdaytsoil(day,i,j,1),i=1,35)
!       enddo
!       write (6,*) 'Daily Top Soil Layer Moisture:'
!       do j = 1,46
!         write(6,*) (mitdaysh2o(day,i,j,1),i=1,35)
!       enddo
!       write (6,*) 'Hourly Top Soil Layer Moisture:'
!       do j = 1,46
!         write(6,*) ((mithrsh2o(hr,i,1,j,k),i=1,5),k=1,4)
!       enddo
!       write (6,*) 'Storm Duration (Hours):'
!       do j = 1,46
!         write(6,*) mitstrmdur(day,j)
!       enddo
!    endif
!
!  Deallocate buffers
!
    deallocate(rbufmlp)
    deallocate(rbufslc)
    deallocate(s1msat)
    deallocate(s2msat)
#endif   ! (COUP_TEM)

#if (defined DGVM)
    ! Accumulate and extract TDA
    ! (accumulates TBOT as 30-day average)
    ! Also determine t_mo_min

!dir$ concurrent
!cdir nodep
    do p = begp,endp
       g = pgridcell(p)
       rbufslp(p) = forc_t(g)
    end do
    call update_accum_field  ('TDA', rbufslp, nstep)
    call extract_accum_field ('TDA', rbufslp, nstep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       t_mo(p) = rbufslp(p)
       t_mo_min(p) = min(t_mo_min(p), rbufslp(p))
    end do

    ! Accumulate and extract T10
    !(acumulates TSA as 10-day running mean)

    call update_accum_field  ('T10', t_ref2m, nstep)
    call extract_accum_field ('T10', t10, nstep)

    ! Accumulate and extract FNPSN10
    !(accumulates fpsn-frmf as 10-day running mean)

!dir$ concurrent
!cdir nodep
    do p = begp,endp
       rbufslp(p) = fpsn(p) - frmf(p)
    end do
    call update_accum_field  ('FNPSN10', rbufslp, nstep)
    call extract_accum_field ('FNPSN10', fnpsn10, nstep)

    ! Accumulate and extract PREC365
    ! (accumulates total precipitation as 365-day running mean)

!dir$ concurrent
!cdir nodep
    do p = begp,endp
       g = pgridcell(p)
       rbufslp(p) = forc_rain(g) + forc_snow(g)
    end do
    call update_accum_field  ('PREC365', rbufslp, nstep)
    call extract_accum_field ('PREC365', prec365, nstep)

    ! Accumulate growing degree days based on 10-day running mean temperature.
    ! Accumulate GDD above 0C and -5C using extracted t10 from accumulated variable.
    ! The trigger to reset the accumulated values to zero is -99999.
    ! agddtw is currently reset at the end of each year in subr. lpj

    ! Accumulate and extract AGDDO

!dir$ concurrent
!cdir nodep
    do p = begp,endp
       rbufslp(p) = (t10(p) - SHR_CONST_TKFRZ) * dtime / SHR_CONST_CDAY
       if (rbufslp(p) < 0._r8) rbufslp(p) = -99999._r8
    end do
    call update_accum_field  ('AGDD0', rbufslp, nstep)
    call extract_accum_field ('AGDD0', agdd0, nstep)

    ! Accumulate and extract AGDD5

!dir$ concurrent
!cdir nodep
    do p = begp,endp
       rbufslp(p) = (t10(p) - (SHR_CONST_TKFRZ - 5.0_r8))*dtime / SHR_CONST_CDAY
       if (rbufslp(p) < 0._r8) rbufslp(p) = -99999._r8
    end do
    call update_accum_field  ('AGDD5', rbufslp, nstep)
    call extract_accum_field ('AGDD5', agdd5, nstep)

    ! Accumulate and extract AGDDTW

!dir$ concurrent
!cdir nodep
    do p = begp,endp
       itypveg = itype(p)
       rbufslp(p) = max(0.0_r8, (t10(p) - (SHR_CONST_TKFRZ+pftpar(itypveg,31))) &
            * dtime/SHR_CONST_CDAY)
    end do
    call update_accum_field  ('AGDDTW', rbufslp, nstep)
    call extract_accum_field ('AGDDTW', agddtw, nstep)

    ! Accumulate and extract AGDD

!dir$ concurrent
!cdir nodep
    do p = begp,endp
       rbufslp(p) = max(0.0_r8, (t_ref2m(p) - (SHR_CONST_TKFRZ + 5.0_r8)) &
            * dtime/SHR_CONST_CDAY)
    end do
    call update_accum_field  ('AGDD', rbufslp, nstep)
    call extract_accum_field ('AGDD', agdd, nstep)
#endif

    ! Deallocate dynamic memory

    deallocate(rbufslp)

  end subroutine updateAccFlds

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: initAccClmtype
!
! !INTERFACE:
  subroutine initAccClmtype
!
! !DESCRIPTION:
! Initialize clmtype variables that are associated with
! time accumulated fields. This routine is called in an initial run
! at nstep=0 for cam and csm mode and at nstep=1 for offline mode.
! This routine is also always called for a restart run and
! therefore must be called after the restart file is read in
! and the accumulated fields are obtained.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use decompMod   , only : get_proc_bounds, get_proc_global
    use accumulMod  , only : extract_accum_field
    use clm_time_manager, only : get_nstep
    use clm_varctl  , only : nsrest
    use clm_varcon  , only : spval
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
! !LOCAL VARIABLES:
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: t_ref2m_min(:)      ! daily minimum of average 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_max(:)      ! daily maximum of average 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_min_inst(:) ! instantaneous daily min of average 2 m height surface air temp (K)
    real(r8), pointer :: t_ref2m_max_inst(:) ! instantaneous daily max of average 2 m height surface air temp (K)
    real(r8), pointer :: t10(:)              ! 10-day running mean of the 2 m temperature (K)
    real(r8), pointer :: t_mo(:)             ! 30-day average temperature (Kelvin)
    real(r8), pointer :: fnpsn10(:)          ! 10-day running mean net photosynthesis
    real(r8), pointer :: prec365(:)          ! 365-day running mean of tot. precipitation
    real(r8), pointer :: agdd0(:)            ! accumulated growing degree days above 0 deg C
    real(r8), pointer :: agdd5(:)            ! accumulated growing degree days above -5
    real(r8), pointer :: agddtw(:)           ! accumulated growing degree days above twmax
    real(r8), pointer :: agdd(:)             ! accumulated growing degree days above 5
    real(r8), pointer :: evap30(:)           ! 30-day average evapotranspiration (mm/day)
    real(r8), pointer :: evap30c(:)          ! 30-day average evapotranspiration (mm/day)
    real(r8), pointer :: surf30(:)           ! 30-day average evapotranspiration (mm/day)
    real(r8), pointer :: drai30(:)           ! 30-day average evapotranspiration (mm/day)
    real(r8), pointer :: swe30(:)            ! 30-day average evapotranspiration (mm/day)
!
! !LOCAL VARIABLES:
!
!EOP
!
! !OTHER LOCAL VARIABLES:
    integer :: p,c          ! indices
    integer :: nstep        ! time step
    integer :: ier          ! error status
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    real(r8), pointer :: rbufslp(:)  ! temporary
    real(r8), pointer :: rbufslc(:)  ! temporary
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (pft-level)

    t_ref2m_max_inst => clm3%g%l%c%p%pes%t_ref2m_max_inst
    t_ref2m_min_inst => clm3%g%l%c%p%pes%t_ref2m_min_inst
    t_ref2m_max      => clm3%g%l%c%p%pes%t_ref2m_max
    t_ref2m_min      => clm3%g%l%c%p%pes%t_ref2m_min
    t10              => clm3%g%l%c%p%pdgvs%t10
    t_mo             => clm3%g%l%c%p%pdgvs%t_mo
    fnpsn10          => clm3%g%l%c%p%pdgvs%fnpsn10
    prec365          => clm3%g%l%c%p%pdgvs%prec365
    agdd0            => clm3%g%l%c%p%pdgvs%agdd0
    agdd5            => clm3%g%l%c%p%pdgvs%agdd5
    agddtw           => clm3%g%l%c%p%pdgvs%agddtw
    agdd             => clm3%g%l%c%p%pdgvs%agdd
    evap30           => clm3%g%l%c%p%pps%evap30
    evap30c          => clm3%g%l%c%cwf%evap30c
    drai30           => clm3%g%l%c%p%pps%drai30
    surf30           => clm3%g%l%c%p%pps%surf30
    swe30            => clm3%g%l%c%p%pps%swe30

    ! Determine necessary indices

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! Determine time step

    nstep = get_nstep()

    ! Initialize 2m ref temperature max and min values

    if (nsrest == 0) then
!dir$ concurrent
!cdir nodep
       do p = begp,endp
          t_ref2m_max(p) = spval
          t_ref2m_min(p) = spval
          t_ref2m_max_inst(p) = -spval
          t_ref2m_min_inst(p) =  spval
       end do
    end if

#if (defined COUP_TEM)

    ! Allocate needed dynamic memory for single level pft field

    allocate(rbufslp(begp:endp), stat=ier)
    allocate(rbufslc(begc:endc), stat=ier)
    if (ier/=0) then
       write(6,*)'update_accum_hist allocation error for rbufslp'
       call endrun
    endif

    ! Initialize clmtype variables that are to be time accumulated

    call extract_accum_field ('EVAP30', rbufslp, nstep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       evap30(p) = rbufslp(p)
    end do

!    call extract_accum_field ('EVAP30C', rbufslc, nstep)
!    do c = begc,endc
!       evap30c(c) = rbufslc(c)
!    end do

    call extract_accum_field ('SURF30', rbufslp, nstep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       surf30(p) = rbufslp(p)
    end do

    call extract_accum_field ('DRAI30', rbufslp, nstep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       drai30(p) = rbufslp(p)
    end do

    call extract_accum_field ('SWE30', rbufslp, nstep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       swe30(p) = rbufslp(p)
    end do
#endif

#if (defined DGVM)

    ! Allocate needed dynamic memory for single level pft field

    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(6,*)'update_accum_hist allocation error for rbufslp'
       call endrun
    endif

    ! Initialize clmtype variables that are to be time accumulated

    call extract_accum_field ('T10', rbufslp, nstep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       t10(p) = rbufslp(p)
    end do

    call extract_accum_field ('TDA', rbufslp, nstep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       t_mo(p) = rbufslp(p)
    end do

    call extract_accum_field ('AGDD0', rbufslp, nstep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       agdd0(p) = rbufslp(p)
    end do

    call extract_accum_field ('AGDD5', rbufslp, nstep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       agdd5(p) = rbufslp(p)
    end do

    call extract_accum_field ('FNPSN10', rbufslp, nstep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       fnpsn10(p) = rbufslp(p)
    end do

    call extract_accum_field ('PREC365', rbufslp, nstep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       prec365(p) = rbufslp(p)
    end do

    call extract_accum_field ('AGDDTW', rbufslp, nstep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       agddtw(p) = rbufslp(p)
    end do

    call extract_accum_field ('AGDD', rbufslp, nstep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       agdd(p) = rbufslp(p)
    end do

    deallocate(rbufslp)
#endif

  end subroutine initAccClmtype

end module accFldsMod
