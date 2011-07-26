#include <misc.h>
#include <preproc.h>

module histFileMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: histFileMod
! 
! !DESCRIPTION: 
! Module containing methods to for clm2 history file handling
! 
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
#if (defined COUP_TEM)
    use clm_varpar, only : lsmlon, lsmlat, nlevsoi
#endif
  implicit none
  save
  private
  include "netcdf.inc"
!
! !PUBLIC TYPES
!
! Constants
!
  integer , public, parameter :: max_tapes = 6          ! max number of history tapes
  integer , public, parameter :: max_flds = 1000        ! max number of history fields 
  integer , public, parameter :: max_namlen = 32        ! maximum number of characters for field name
  real(r8), public, parameter :: spval = 1.e36          ! fill value
!
! Counters
!
  integer , public :: ntapes = 0          ! index of max history file requested 
!
! Namelist
!
  integer :: ni                          ! implicit index below
  logical, public :: &
       hist_empty_htapes  = .false.      ! namelist: flag indicates no default history fields
  integer, public :: &
       hist_ndens(max_tapes) = 2         ! namelist: output density of netcdf history files
  integer, public :: &
       hist_mfilt(max_tapes) = 30        ! namelist: number of time samples per tape
  logical, public :: &
       hist_dov2xy(max_tapes) = (/.true.,(.true.,ni=2,max_tapes)/) ! namelist: true=> do grid averaging
#if (defined COUP_TEM)
  logical, public :: ready4tem = .false.
#endif
  integer, public :: &
       hist_nhtfrq(max_tapes) = (/0, (-24, ni=2,max_tapes)/)        ! namelist: history write freq(0=monthly)
  character(len=1), public :: &
       hist_avgflag_pertape(max_tapes) = (/(' ',ni=1,max_tapes)/)   ! namelist: per tape averaging flag
  character(len=max_namlen), public :: &
       hist_type1d_pertape(max_tapes)  = (/(' ',ni=1,max_tapes)/)   ! namelist: per tape type1d

  character(len=max_namlen+2), public :: &
       hist_fincl1(max_flds) = (/(' ',ni=1,max_flds)/) ! namelist: list of fields to add 
  character(len=max_namlen+2), public :: &
       hist_fincl2(max_flds) = (/(' ',ni=1,max_flds)/) ! namelist: list of fields to add  
  character(len=max_namlen+2), public :: &
       hist_fincl3(max_flds) = (/(' ',ni=1,max_flds)/) ! namelist: list of fields to add  
  character(len=max_namlen+2), public :: &
       hist_fincl4(max_flds) = (/(' ',ni=1,max_flds)/) ! namelist: list of fields to add  
  character(len=max_namlen+2), public :: &
       hist_fincl5(max_flds) = (/(' ',ni=1,max_flds)/) ! namelist: list of fields to add  
  character(len=max_namlen+2), public :: &
       hist_fincl6(max_flds) = (/(' ',ni=1,max_flds)/) ! namelist: list of fields to add  
  character(len=max_namlen+2), public :: &
       fincl(max_flds,max_tapes)    ! namelist-equivalence list of fields to add 

  character(len=max_namlen), public :: &
       hist_fexcl1(max_flds) = (/(' ',ni=1,max_flds)/) ! namelist: list of fields to remove 
  character(len=max_namlen), public :: &
       hist_fexcl2(max_flds) = (/(' ',ni=1,max_flds)/) ! namelist: list of fields to remove 
  character(len=max_namlen), public :: &
       hist_fexcl3(max_flds) = (/(' ',ni=1,max_flds)/) ! namelist: list of fields to remove 
  character(len=max_namlen), public :: &
       hist_fexcl4(max_flds) = (/(' ',ni=1,max_flds)/) ! namelist: list of fields to remove 
  character(len=max_namlen), public :: &
       hist_fexcl5(max_flds) = (/(' ',ni=1,max_flds)/) ! namelist: list of fields to remove 
  character(len=max_namlen), public :: &
       hist_fexcl6(max_flds) = (/(' ',ni=1,max_flds)/) ! namelist: list of fields to remove
  character(len=max_namlen), public :: &
       fexcl(max_flds,max_tapes)    ! namelist-equivalence list of fields to remove
!
! Equivalence used to satisfy namelist on a wide variety of platforms
! NOTE: It is *ASSUMED* that max_tapes is 6
!  
  equivalence (hist_fincl1,fincl(1,1))
  equivalence (hist_fincl2,fincl(1,2))
  equivalence (hist_fincl3,fincl(1,3))
  equivalence (hist_fincl4,fincl(1,4))
  equivalence (hist_fincl5,fincl(1,5))
  equivalence (hist_fincl6,fincl(1,6))

  equivalence (hist_fexcl1,fexcl(1,1))
  equivalence (hist_fexcl2,fexcl(1,2))
  equivalence (hist_fexcl3,fexcl(1,3))
  equivalence (hist_fexcl4,fexcl(1,4))
  equivalence (hist_fexcl5,fexcl(1,5))
  equivalence (hist_fexcl6,fexcl(1,6))
!
! Restart
!
  logical, public :: if_writrest    ! true=> write restart file now  
!
! !PUBLIC MEMBER FUNCTIONS
  public  :: masterlist_build ! Build master field list of all possible history file fields 
  public  :: htapes_build     ! Initialize history file handler for initial or continue run
  public  :: htapes_wrapup    ! Write and/or dispose history tape(s) 
  public  :: update_hbuf      ! Updates history buffer 
  public  :: restart_history  ! Read/write history file restart data
!
! !REVISION HISTORY
! Created by Mariana Vertenstein
!
!EOP
!
! PRIVATE MEMBER FUNCTIONS
  private :: masterlist_addfld         ! Add a field to the master field list   
  private :: masterlist_change_active  ! Define default contents of history files
  private :: masterlist_make_active    ! Add a field to the given history file default "on" list 
  private :: masterlist_change_timeavg ! Override default history tape contents for specific tape
  private :: htapes_fieldlist          ! Define the contents of each history file based on namelist  
  private :: htape_addfld              ! Add a field to the active list for a history tape
  private :: htape_create              ! Define contents of history file t
  private :: htape_timeconst           ! Write time constant values to primary history tape    
  private :: hfields_normalize         ! Normalize history file fields by number of accumulations 
  private :: hfields_zero              ! Zero out accumulation and hsitory buffers for a tape
  private :: hfields_write             ! Write a variable to a history tape 
  private :: list_index                ! Find index of field in exclude list
  private :: getname                   ! Retrieve name portion of input "inname"
  private :: set_hist_filename         ! Determine history dataset filenames    

! PRIVATE TYPES
! Constants
!
  integer, parameter :: max_chars = 128        ! max chars for char variables
  integer, parameter :: not_valid = -999       ! value for not defined hpindex
!
! Derived types
!
  type field_info
     character(len=max_namlen) :: name         ! field name
     character(len=max_chars) :: long_name     ! long name
     character(len=max_chars) :: units         ! units
     integer :: numlev                         ! vertical dimension 
     integer :: beg1d                          ! on-node 1d start index
     integer :: end1d                          ! on-node 1d end index
     integer :: num1d                          ! total number of 1d [gridcells,landunits,columns, or patches] (overall all nodes)
     character(len=8) :: type1d                ! 1d vertical type ("gridcell","landunit","column" or "pft")
     integer :: hpindex                        ! history pointer index selected out of possible list
     integer :: hpindices(4)                   ! possible history pointer indices for different 1d choices
  end type field_info
  
  type master_entry
     type (field_info)  :: field               ! field information
     logical            :: actflag(max_tapes)  ! active/inactive flag
     character(len=1)   :: avgflag(max_tapes)  ! time averaging flag ("X","A","M" or "I",)
  end type master_entry
  
  type history_entry
     type (field_info) :: field                ! field information
     character(len=1)  :: avgflag              ! time averaging flag
     real(r8), pointer :: hbuf(:,:)            ! history buffer (dimensions: dim1d x numlev)
     integer , pointer :: nacs(:,:)            ! accumulation counter (dimensions: dim1d x numlev)
  end type history_entry                       
                                               
  type history_tape                            
     integer  :: nflds                         ! number of active fields on tape
     integer  :: ntimes                        ! current number of time samples on tape
     integer  :: mfilt                         ! maximum number of time samples per tape      
     integer  :: nhtfrq                        ! number of time samples per tape
     integer  :: ncprec                        ! netcdf output precision
     logical  :: dov2xy                        ! true => do xy average for all fields
     logical  :: is_endhist                    ! true => current time step is end of history interval 
     real(r8) :: begtime                       ! time at beginning of history averaging interval
     type (history_entry) :: hlist(max_flds)   ! array of active history tape entries
  end type history_tape                          
!                                                 
! Master list: an array of master_entry entities
!
  type (master_entry) :: masterlist(max_flds)  ! master field list
!  
! History tape: an array of history_tape entities (only active fields)
!  
  type (history_tape) :: tape(max_tapes)       ! array history tapes
!  
! Namelist input
!  
!
! Counters
!
  integer :: nfmaster = 0                        ! number of fields in master field list
!  
! Other variables
!  
  character(len=16) :: host                      ! host name
  character(len= 8) :: logname                   ! user name
  character(len=max_chars) :: locfnh(max_tapes)  ! local history file names
  logical :: htapes_defined = .false.            ! flag indicates history contents have been defined
!
! NetCDF dimension Id's
!  
  integer :: lon_dimid                      ! longitude dimension id
  integer :: lat_dimid                      ! latitude dimension id
  integer :: levsoi_dimid                   ! soil layer dimension id
  integer :: time_dimid                     ! time dimension id
  integer :: hist_interval_dimid            ! time bounds dimension id
  integer :: strlen_dimid                   ! string dimension id
  integer :: gridcell_dimid                 ! 1d grid dimension id
  integer :: landunit_dimid                 ! 1d landunit dimension id
  integer :: column_dimid                   ! 1d column dimension id
  integer :: pft_dimid                      ! 1d pft dimension id

! NetCDF dimension time invariant Id's

  integer :: nfid(max_tapes)             ! file ids
  integer :: lonvar_id(max_tapes)        ! full grid longitude coordinate variable ids
  integer :: latvar_id(max_tapes)        ! full grid latitude  coordinate variable ids
  integer :: levvar_id(max_tapes)        ! soil level coordinate variable ids
  integer :: longxy_id(max_tapes)        ! 2d longitudes (longxy) ids 
  integer :: latixy_id(max_tapes)        ! 2d latitudes (latixy) ids
  integer :: area_id(max_tapes)          ! 2d area (area) ids
  integer :: landfrac_id(max_tapes)      ! 2d land fraction ids
  integer :: numlon_id(max_tapes)        ! number of longitudes at each latitude ids
  integer :: landmask_id(max_tapes)      ! 2d land/ocean mask (landmask) ids
#if (defined OFFLINE)
  integer :: edgen_id(max_tapes)         ! northern edge of grid (lsmedge(1)) ids
  integer :: edgee_id(max_tapes)         ! eastern  edge of grid (lsmedge(2)) ids
  integer :: edges_id(max_tapes)         ! southern edge of grid (lsmedge(3)) ids 
  integer :: edgew_id(max_tapes)         ! western  edge of grid (lsmedge(4)) ids 
#endif

! NetCDF dimension time variant Id's

  integer :: varid(max_flds,max_tapes)   ! variable ids
  integer :: mcdate_id(max_tapes)        ! current date, yyyymmdd format (mcdate) ids
  integer :: mcsec_id(max_tapes)         ! current seconds in day (mcsec) ids
  integer :: mdcur_id(max_tapes)         ! current day (from base day) (mdcur) ids
  integer :: mscur_id(max_tapes)         ! current seconds of current day (mdcur) ids
  integer :: nstep_id(max_tapes)         ! current nstep ids 
  integer :: time_var_id(max_tapes)      ! time coordinate variable ids
  integer :: time_bounds_id(max_tapes)   ! time bounds ids
  integer :: date_written_id(max_tapes)  ! current system date ids
  integer :: time_written_id(max_tapes)  ! current system time ids
  
  integer :: zsoi_id(max_tapes)          ! time invariant soil depth id
  integer :: dzsoi_id(max_tapes)         ! time invariant soil thickness id
  integer :: watsat_id(max_tapes)        ! time invariant soil porosity id
  integer :: hksat_id(max_tapes)         ! time invariant sat. hydraulic conductivity
  integer :: sucsat_id(max_tapes)        ! time invariant soil maxtrix potential id
  integer :: bsw_id(max_tapes)           ! time invariant soil water retention slope id
#if (defined COUP_TEM) 
    integer, parameter :: mxmsaics = 35
    integer, parameter :: mxmthdys = 31
    integer, parameter :: mxdayhrs = 24
    integer, parameter :: max3hrs = 8
    integer, parameter :: mxnlayers = 6
    integer, parameter :: totlayers = 10
    real(r8) :: mitco2(lsmlat)
    real(r8) :: mito3(max3hrs,lsmlat)
    real(r8) :: mittemp(lsmlat)
    real(r8) :: mitdaytemp(mxmthdys,lsmlat)
    real(r8) :: mitswrs(lsmlat)
    real(r8) :: mitpre(lsmlat)
    real(r8) :: mitstrmdur(mxmthdys,lsmlat)
    real(r8) :: mitqstrm(mxmthdys,lsmlat)
    real(r8) :: mitaet(mxmsaics,lsmlat)
    real(r8) :: mitsh2o1m(mxmsaics,lsmlat)
    real(r8) :: mitsh2o2m(mxmsaics,lsmlat)
    real(r8) :: mitswe(mxmsaics,lsmlat)
    real(r8) :: mitsfr(mxmsaics,lsmlat)
    real(r8) :: mitdrn(mxmsaics,lsmlat)
    real(r8) :: mitdaytsoil(mxmthdys,mxmsaics,lsmlat,totlayers)
    real(r8) :: mitdaysh2o(mxmthdys,mxmsaics,lsmlat,totlayers)
    real(r8) :: mithrsh2o(mxdayhrs,mxmthdys,mxmsaics,lsmlat,mxnlayers)
    common/climate4tem/mitco2,mito3,mittemp,mitdaytemp,mitswrs,mitpre,mitstrmdur,&
                       mitqstrm,mitaet,mitsh2o1m,mitsh2o2m,mitswe,mitsfr,mitdrn,&
                       mitdaytsoil,mitdaysh2o,&
                       mithrsh2o
#endif
!----------------------------------------------------------------------- 
  
contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: masterlist_build
!
! !INTERFACE:
  subroutine masterlist_build ()
!
! !DESCRIPTION: 
! Build master field list of all possible fields in a history file.  Each field has 
! associated with it a "long_name" netcdf attribute that describes what the field is, 
! and a "units" attribute. A subroutine is called to add each field to the masterlist.
! The array hpindices has the form (gridcell_index, landunit_index, column_index,
! pft_index). The value of type1d massed to masterlist_add_fld determines which of
! the hpindices array values is used to determine the 1d type (gridcell, landunit,
! column or pft) of the output history buffer (i.e. it determines beg1d and end1d
! of the history buffer field)
!
! !USES
    use clmpoint
    use clm_varpar, only : nlevsoi
    use spmdMod, only : masterproc
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY
! Created by Mariana Vertenstein 
!
!EOP
!
! LOCAL VARIABLES
    integer :: nf            ! masterlist field counter
    integer :: hpindices(4)  ! pointer indices into clmtype derived types
!-----------------------------------------------------------------------

    ! Snow properties (will be vertically averaged over the snow profile)
    
    hpindices = (/-1, -1, ic_cps_snowdp, not_valid/)
    call masterlist_addfld (fname='SNOWDP', type1d='column', units='m', numlev=1, &
         avgflag='A', long_name='snow height', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cps_frac_sno, not_valid/)
    call masterlist_addfld (fname='SNOWCOV', type1d='column', units='%', numlev=1, &
         avgflag='A', long_name='snow cover', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cps_snowage, not_valid/)
    call masterlist_addfld (fname='SNOWAGE', type1d='column', units='unitless', numlev=1, &
         avgflag='A', long_name='snow age', hpindices=hpindices)

    ! Temperatures
    
    hpindices = (/-1, -1, ic_ces_pes_a_t_ref2m, ip_pes_t_ref2m/)
    call masterlist_addfld (fname='TSA', type1d='column', units='K', numlev=1, &
         avgflag='A', long_name='2m air temperature', hpindices=hpindices)        
                                                                       
    hpindices = (/-1, -1, ic_ces_pes_a_t_af, ip_pes_t_af/)
    call masterlist_addfld (fname='TCAS', type1d='column', units='K', numlev=1, &
         avgflag='A', long_name='canopy air temperature', hpindices=hpindices)        
                                                                       
    hpindices = (/-1, -1, ic_ces_pes_a_t_veg, ip_pes_t_veg/)
    call masterlist_addfld (fname='TV', type1d='column', units='K', numlev=1, &
         avgflag='A', long_name='vegetation temperature', hpindices=hpindices)    
                                                                       
    hpindices = (/-1, -1, ic_ces_t_grnd, not_valid/)
    call masterlist_addfld (fname='TG', type1d='column', units='K', numlev=1, &
         avgflag='A', long_name='ground temperature', hpindices=hpindices)        
                                                                       
    hpindices = (/-1, -1, ic_ces_t_snow, not_valid/)
    call masterlist_addfld (fname='TSNOW', type1d='column', units='K', numlev=1, &
         avgflag='A', long_name='snow temperature', hpindices=hpindices)          
                                                                       
    hpindices = (/-1, -1, ic_ces_t_soisno, not_valid/)
    call masterlist_addfld (fname='TSOI', type1d='column', units='K', numlev=nlevsoi, &
         avgflag='A', long_name='soil temperature', hpindices=hpindices)    
                                                                       
    hpindices = (/-1, -1, ic_ces_t_lake, not_valid/)
    call masterlist_addfld (fname='TLAKE', type1d='column', units='K', numlev=nlevsoi, &
         avgflag='A', long_name='lake temperature', hpindices=hpindices)

    ! Surface radiation                                          

    hpindices = (/-1, -1, ic_cef_pef_a_fsa, ip_pef_fsa/)
    call masterlist_addfld (fname='FSA', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='absorbed solar radiation', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cef_pef_a_parsun, ip_pef_parsun/)
    call masterlist_addfld (fname='PARSUN', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='absorbed PAR direct', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cef_pef_a_fsr, ip_pef_fsr/)
    call masterlist_addfld (fname='FSR', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='reflected solar radiation', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cps_pps_a_ndvi, ip_pps_ndvi/)
    call masterlist_addfld (fname='NDVI', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='surface ndvi', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cef_pef_a_eflx_lwrad_net, ip_pef_eflx_lwrad_net/)
    call masterlist_addfld (fname='FIRA', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='net infrared (longwave) radiation', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cef_pef_a_eflx_lwrad_out, ip_pef_eflx_lwrad_out/)
    call masterlist_addfld (fname='FIRE', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='emitted infrared (longwave) radiation', hpindices=hpindices)

    ! Surface energy fluxes                                      

    hpindices = (/-1, -1, ic_cef_pef_a_eflx_lh_vegt, ip_pef_eflx_lh_vegt/)
    call masterlist_addfld (fname='FCTR', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='canopy transpiration', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cef_pef_a_eflx_lh_vege, ip_pef_eflx_lh_vege/)
    call masterlist_addfld (fname='FCEV', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='canopy evaporation', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cef_pef_a_eflx_lh_grnd, ip_pef_eflx_lh_grnd/)
    call masterlist_addfld (fname='FGEV', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='ground evaporation', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cef_pef_a_eflx_sh_tot, ip_pef_eflx_sh_tot/)
    call masterlist_addfld (fname='FSH', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='sensible heat', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cef_pef_a_eflx_soil_grnd, ip_pef_eflx_soil_grnd/)
    call masterlist_addfld (fname='FGR', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='heat flux into soil', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cef_eflx_snomelt, not_valid/)
    call masterlist_addfld (fname='FSM', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='snow melt heat flux', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cmf_pmf_a_taux, ip_pmf_taux/)
    call masterlist_addfld (fname='TAUX', type1d='column', units='kg/m/s^2', numlev=1, &
         avgflag='A', long_name='zonal surface stress', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cmf_pmf_a_tauy, ip_pmf_tauy/)
    call masterlist_addfld (fname='TAUY', type1d='column', units='kg/m/s^2', numlev=1, &
         avgflag='A', long_name='meridional surface stress', hpindices=hpindices)

    ! Vegetation phenology

    hpindices = (/-1, -1, ic_cps_pps_a_elai, ip_pps_elai/)
    call masterlist_addfld (fname='ELAI', type1d='column', units='m^2/m^2', &
         numlev=1, avgflag='A', long_name='exposed one-sided leaf area index', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cps_pps_a_esai, ip_pps_esai/)
    call masterlist_addfld (fname='ESAI', type1d='column', units='m^2/m^2', &
         numlev=1, avgflag='A', long_name='exposed one-sided stem area index', hpindices=hpindices)

    ! Canopy physiology                                          

    hpindices = (/-1, -1, ic_cps_pps_a_rssun, ip_pps_rssun/)
    call masterlist_addfld (fname='RSSUN', type1d='column', units='s/m', numlev=1, &
         avgflag='M', long_name='sunlit leaf stomatal resistance', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cps_pps_a_rssha, ip_pps_rssha/)
    call masterlist_addfld (fname='RSSHA', type1d='column', units='s/m', numlev=1, &
         avgflag='M', long_name='shaded leaf stomatal resistance', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cps_pps_a_btran, ip_pps_btran/)
    call masterlist_addfld (fname='BTRAN', type1d='column', units='s/m', numlev=1, &
         avgflag='A', long_name='transpiration beta factor', hpindices=hpindices)

    ! Hydrology

    hpindices = (/-1, -1, ic_cws_h2osno, not_valid/)
    call masterlist_addfld (fname='H2OSNO', type1d='column', units='mm', numlev=1, &
         avgflag='A', long_name='snow depth (liquid water)', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cws_pws_a_h2ocan, ip_pws_h2ocan/)
    call masterlist_addfld (fname='H2OCAN', type1d='column', units='mm', numlev=1, &
         avgflag='A', long_name='intercepted water', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cws_h2osoi_vol, not_valid/)
    call masterlist_addfld (fname='H2OSOI', type1d='column', units='mm3/mm3', numlev=nlevsoi, &
         avgflag='A', long_name='volumetric soil water', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cws_h2osoi_liq, not_valid/)
    call masterlist_addfld (fname='SOILLIQ', type1d='column', units='kg/m2', numlev=nlevsoi, &
         avgflag='A', long_name='soil liquid water', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cws_h2osoi_ice, not_valid/)
    call masterlist_addfld (fname='SOILICE', type1d='column', units='kg/m2', numlev=nlevsoi, &
         avgflag='A', long_name='soil liquid ice', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cws_snowliq, not_valid/)
    call masterlist_addfld (fname='SNOWLIQ', type1d='column', units='kg/m2', numlev=1, &
         avgflag='A', long_name='snow liquid water', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cws_snowice, not_valid/)
    call masterlist_addfld (fname='SNOWICE', type1d='column', units='kg/m2', numlev=1, &
         avgflag='A', long_name='snow liquid ice', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cwf_qflx_infl, not_valid/)
    call masterlist_addfld (fname='QINFL', type1d='column', units='mm/s', numlev=1, &
         avgflag='A', long_name='infiltration', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cwf_qflx_surf, not_valid/)
    call masterlist_addfld (fname='QOVER', type1d='column', units='mm/s', numlev=1, &
         avgflag='A', long_name='surface runoff', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cwf_qflx_qrgwl, not_valid/)
    call masterlist_addfld (fname='QRGWL', type1d='column', units='mm/s', numlev=1, &
         avgflag='A', long_name='surface runoff at glaciers, wetlands, lakes', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cwf_qflx_drain, not_valid/)
    call masterlist_addfld (fname='QDRAI', type1d='column', units='mm/s', numlev=1, &
         avgflag='A', long_name='sub-surface drainage', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cwf_pwf_a_qflx_prec_intr, ip_pwf_qflx_prec_intr/)
    call masterlist_addfld (fname='QINTR', type1d='column', units='mm/s', numlev=1, &
         avgflag='A', long_name='interception', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cwf_pwf_a_qflx_prec_grnd, ip_pwf_qflx_prec_grnd/)
    call masterlist_addfld (fname='QDRIP', type1d='column', units='mm/s', numlev=1, &
         avgflag='A', long_name='throughfall', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cwf_qflx_snomelt, not_valid/)
    call masterlist_addfld (fname='QMELT', type1d='column', units='mm/s', numlev=1, &
         avgflag='A', long_name='snow melt', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cwf_pwf_a_qflx_evap_soi, ip_pwf_qflx_evap_soi/)
    call masterlist_addfld (fname='QSOIL', type1d='column', units='mm/s', numlev=1, &
         avgflag='A', long_name='ground evaporation', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cwf_pwf_a_qflx_evap_can, ip_pwf_qflx_evap_can/)
    call masterlist_addfld (fname='QVEGE', type1d='column', units='mm/s', numlev=1, &
         avgflag='A', long_name='canopy evaporation', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cwf_pwf_a_qflx_tran_veg, ip_pwf_qflx_tran_veg/)
    call masterlist_addfld (fname='QVEGT', type1d='column', units='mm/s', numlev=1, &
         avgflag='A', long_name='canopy transpiration', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cwf_pwf_a_qflx_efpot, ip_pwf_qflx_efpot/)
    call masterlist_addfld (fname='PET', type1d='column', units='mm/s', numlev=1, &
         avgflag='A', long_name='potential transpiration', hpindices=hpindices)

#if (defined RTM)
    ! RTM River Routing

    hpindices = (/-1, -1, ic_cwf_qchan2, not_valid/)
    call masterlist_addfld (fname='QCHANR', type1d='column', units='m3/s', numlev=1, &
         avgflag='A', long_name='RTM river flow (maximum subgrid flow)', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cwf_qchocn2, not_valid/)
    call masterlist_addfld (fname='QCHOCNR', type1d='column', units='m3/s', numlev=1, &
         avgflag='A', long_name='RTM river discharge into ocean', hpindices=hpindices)
#endif

    ! Water and energy balance checks

    hpindices = (/-1, -1, ic_cebal_errsoi, not_valid/)
    call masterlist_addfld (fname='ERRSOI', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='soil/lake energy conservation error', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cebal_errseb, not_valid/)
    call masterlist_addfld (fname='ERRSEB', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='surface energy conservation error', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cebal_errsol, not_valid/)
    call masterlist_addfld (fname='ERRSOL', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='solar radiation conservation error', hpindices=hpindices)

    hpindices = (/-1, -1, ic_cwbal_errh2o, not_valid/)
    call masterlist_addfld (fname='ERRH2O', type1d='column', units='mm', numlev=1, &
         avgflag='A', long_name='total water conservation error', hpindices=hpindices)

    ! Atmospheric forcing                                

    hpindices = (/ig_a2lf_forc_rain, not_valid, not_valid, not_valid/)
    call masterlist_addfld (fname='RAIN', type1d='gridcell', units='mm/s', numlev=1, &
         avgflag='A', long_name='atmospheric rain', hpindices=hpindices)

    hpindices = (/ig_a2lf_forc_snow, not_valid, not_valid, not_valid/)
    call masterlist_addfld (fname='SNOW', type1d='gridcell', units='mm/s', numlev=1, &
         avgflag='A', long_name='atmospheric snow', hpindices=hpindices)

    hpindices = (/ig_a2lf_forc_strm_int, not_valid, not_valid, not_valid/)
    call masterlist_addfld (fname='QSTRM', type1d='gridcell', units='mm/s', numlev=1, &
         avgflag='A', long_name='Storm Intensity', hpindices=hpindices)

#if (defined STOCHASTIC)
    hpindices = (/ig_a2lf_forc_strm_dur, not_valid, not_valid, not_valid/)
    call masterlist_addfld (fname='STRMDUR', type1d='gridcell', units='s', numlev=1, &
         avgflag='A', long_name='Storm Duration', hpindices=hpindices)

    hpindices = (/ig_a2lf_forc_strm_dry, not_valid, not_valid, not_valid/)
    call masterlist_addfld (fname='STRMDRY', type1d='gridcell', units='s', numlev=1, &
         avgflag='I', long_name='Inter-Storm Period', hpindices=hpindices)

    hpindices = (/ig_a2lf_forc_strms, not_valid, not_valid, not_valid/)
    call masterlist_addfld (fname='STORMS', type1d='gridcell', units='# of events', numlev=1, &
         avgflag='X', long_name='# of storms', hpindices=hpindices)
#endif

    hpindices = (/ig_a2ls_forc_t, not_valid, not_valid, not_valid/)
    call masterlist_addfld (fname='TBOT', type1d='gridcell', units='K', numlev=1, &
         avgflag='A', long_name='atmospheric air temperature', hpindices=hpindices)

    hpindices = (/ig_a2ls_forc_th, not_valid, not_valid, not_valid/)
    call masterlist_addfld (fname='THBOT', type1d='gridcell', units='K', numlev=1, &
         avgflag='A', long_name='atmospheric air potential temperature', hpindices=hpindices)

    hpindices = (/ig_a2ls_forc_wind, not_valid, not_valid, not_valid/)
    call masterlist_addfld (fname='WIND', type1d='gridcell', units='m/s', numlev=1, &
         avgflag='A', long_name='atmospheric wind velocity magnitude', hpindices=hpindices)

    hpindices = (/ig_a2ls_forc_q, not_valid, not_valid, not_valid/) 
    call masterlist_addfld (fname='QBOT', type1d='gridcell', units='kg/kg', numlev=1, &
         avgflag='A', long_name='atmospheric specific humidity', hpindices=hpindices)

    hpindices = (/ig_a2ls_forc_hgt, not_valid, not_valid, not_valid/) 
    call masterlist_addfld (fname='ZBOT', type1d='gridcell', units='m', numlev=1, &
         avgflag='A', long_name='atmospheric reference height', hpindices=hpindices)

    hpindices = (/ig_a2lf_forc_lwrad, not_valid, not_valid, not_valid/) 
    call masterlist_addfld (fname='FLDS', type1d='gridcell', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='atmospheric longwave radiation', hpindices=hpindices)

    hpindices = (/ig_a2lf_forc_solar,  not_valid, not_valid, not_valid/) 
    call masterlist_addfld (fname='FSDS', type1d='gridcell', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='atmospheric incident solar radiation', hpindices=hpindices)

    ! Print master field list

    if (masterproc) then
       write(6,*)'MASTERLIST_BUILD: number of master fields = ',nfmaster
       write(6,*)' ******* MASTER FIELD LIST *******'
       do nf = 1,nfmaster
          write(6,9000)nf, masterlist(nf)%field%name, masterlist(nf)%field%units
9000      format (i5,1x,a8,1x,a8)
       end do
    end if
    
  end subroutine masterlist_build
  
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: masterlist_addfld
!
! !INTERFACE:
  subroutine masterlist_addfld (fname, type1d, units, numlev, avgflag, &
       long_name, hpindices)
!
! !DESCRIPTION: 
! Add a field to the master field list. Put input arguments of 
! field name, units, number of levels, averaging flag, and long name 
! into a type entry in the global master field list (masterlist).
!
! !USES
    use clmtype, only : grid1d, land1d, cols1d, pfts1d
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: fname      ! field name
    character(len=*), intent(in) :: type1d     ! 1d vertical type ("gridcell","landunit","column" or "pft")
    character(len=*), intent(in) :: units      ! units of field
    character(len=1), intent(in) :: avgflag    ! time averaging flag
    character(len=*), intent(in) :: long_name  ! long name of field
    integer, intent(in) :: numlev              ! number of vertical levels (dimension and loop)
    integer, intent(in) :: hpindices(4)        ! possible list of clm pointer indices
!
! !REVISION HISTORY
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES
    integer :: n     ! loop index
    integer :: f     ! masterlist index  
    integer :: hpindex ! appropriate index from list of hpindices
!-----------------------------------------------------------------------

    ! Ensure that new field is not all blanks

    if (fname == ' ') then
       write(6,*)'MASTERLIST_ADDFLD: blank field name not allowed'
    end if
    
    ! Ensure that new field doesn't already exist
    
    do n = 1,nfmaster
       if (masterlist(n)%field%name == fname) then
          write(6,*)'MASTERLIST_ADDFLD:', fname, ' already on list'
          call endrun ()
       end if
    end do
    
    ! Increase number of fields on master field list

    nfmaster = nfmaster + 1  
    f = nfmaster
    
    ! Check number of fields in master list against maximum number for master list
    
    if (nfmaster > max_flds) then
       write(6,*)'MASTERLIST_ADDFLD: Too many fields for primary history file -- max_flds,nfmaster=', &
            max_flds, nfmaster
       call endrun
    end if
    
    ! Add field to master list
    
    masterlist(f)%field%name      = fname
    masterlist(f)%field%long_name = long_name
    masterlist(f)%field%units     = units
    masterlist(f)%field%numlev    = numlev
    masterlist(f)%field%hpindices(:) = hpindices(:)
    masterlist(f)%field%type1d    = type1d

    select case (type1d)
    case ('gridcell') 
       hpindex = hpindices(1)
       if (hpindex == -1) then
          write(6,*)' a field was chosen incorrectly for gridcell 1d output'
	  call endrun
       endif	
       masterlist(f)%field%hpindex = hpindex
       masterlist(f)%field%beg1d = grid1d%beg
       masterlist(f)%field%end1d = grid1d%end
       masterlist(f)%field%num1d = grid1d%num
    case ('landunit') 
       hpindex = hpindices(2)
       if (hpindex == -1) then
          write(6,*)' a field was chosen incorrectly for landunit 1d output'
	  call endrun
       endif	
       masterlist(f)%field%hpindex = hpindex
       masterlist(f)%field%beg1d = land1d%beg
       masterlist(f)%field%end1d = land1d%end
       masterlist(f)%field%num1d = land1d%num
    case ('column')
       hpindex = hpindices(3)
       if (hpindex == -1) then
          write(6,*)' a field was chosen incorrectly for landunit 1d output'
	  call endrun
       endif	
       masterlist(f)%field%hpindex = hpindex
       masterlist(f)%field%beg1d = cols1d%beg
       masterlist(f)%field%end1d = cols1d%end
       masterlist(f)%field%num1d = cols1d%num
    case ('pft')
       hpindex = hpindices(4)
       if (hpindex == -1) then
          write(6,*)' a field was chosen incorrectly for pft 1d output'
	  call endrun
       endif	
       masterlist(f)%field%hpindex = hpindex
       masterlist(f)%field%beg1d = pfts1d%beg
       masterlist(f)%field%end1d = pfts1d%end
       masterlist(f)%field%num1d = pfts1d%num
    case default
       write(6,*)'MASTERLIST_ADDFLD: unknown 1d vertical type=',type1d
       call endrun ()
    end select
    
    ! The following two fields are used only in master field list, 
    ! NOT in the runtime active field list
    ! ALL FIELDS IN THE MASTER LIST ARE INITIALIZED WITH THE ACTIVE
    ! FLAG SET TO FALSE 
    
    masterlist(f)%avgflag(:) = avgflag
    masterlist(f)%actflag(:) = .false.
    
    return
  end subroutine masterlist_addfld

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: htapes_build
!
! !INTERFACE:
  subroutine htapes_build ()
!
! !DESCRIPTION: 
! Initialize history file for initial or continuation run.
! For example, on an initial run, this routine initializes "ntapes"
! history files.  On a restart run, this routine only initializes history 
! files declared beyond what existed on the previous run.  
! Files which already existed on the previous run have 
! already been initialized (i.e. named and opened) in routine restart_history.
! Loop over tapes and fields per tape setting appropriate variables and
! calling appropriate routines
!
! !USES
    use spmdMod, only : masterproc
    use time_manager, only: get_curr_date, get_curr_time
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES
    integer :: i                   ! index
    integer :: t, f                ! tape, field indices
    integer :: day, sec            ! day and seconds from base date
    integer :: numlev              ! number of vertical levels (dimension and loop)
    integer :: beg1d,end1d         ! 1d beginning and ending indices
!-----------------------------------------------------------------------

    if (masterproc) then
       write(6,*) 'Initializing clm2 history files'
       write(6,'(72a1)') ("-",i=1,60)
    endif

    ! Get users logname and machine hostname
     
    if ( masterproc )then
       logname = ' '
       call getenv ('LOGNAME',logname)
       if (logname=='        ') then
          write(6,*)'PARSE_NAMELIST: Cannot find LOGNAME environment variable'
          call endrun
       end if
       host = ' '
       call getenv ('HOST',host)
    end if
    
    ! Set default history contents for all tapes
    
    call masterlist_change_active ()
    
    ! Override averaging flag for all fields on a particular tape 
    ! if namelist input so specifies
    
    do t=1,max_tapes
       if (hist_avgflag_pertape(t) /= ' ') then
          call masterlist_change_timeavg (t)
       end if
    end do
    
    ! Define field list information for all history files.  
    ! Update ntapes to reflect number of active history files 
    ! (note, branch runs can have additional auxiliary history files
    ! declared).
    
    call htapes_fieldlist ()
    
    ! Determine elapased time since reference date
    
    call get_curr_time(day, sec)  
    
    ! Set number of time samples in each history file and 
    ! time of a beginning of current averaging interval.
    ! (note - the following entries will be overwritten by history restart)
    ! Determine if xy averaging is done for all fields on the tape, etc
    ! Note - with netcdf, only 1 (nf_double) and 2 (nf_float) are allowed

    do t=1,ntapes
       tape(t)%ntimes = 0            
       tape(t)%begtime = day + sec/86400._r8
       tape(t)%dov2xy = hist_dov2xy(t)
       tape(t)%nhtfrq = hist_nhtfrq(t)
       if (hist_nhtfrq(t) == 0) hist_mfilt(t) = 1
       tape(t)%mfilt = hist_mfilt(t)
       if (hist_ndens(t) == 1) then
          tape(t)%ncprec = nf_double
       else
          tape(t)%ncprec = nf_float
       endif
    end do
    
    ! Alloccate and initialize history buffer and related info
    
    do t = 1,ntapes
       do f = 1,tape(t)%nflds
          numlev = tape(t)%hlist(f)%field%numlev
          beg1d = tape(t)%hlist(f)%field%beg1d
          end1d = tape(t)%hlist(f)%field%end1d
          allocate (tape(t)%hlist(f)%hbuf(beg1d:end1d,numlev))
          tape(t)%hlist(f)%hbuf(:,:) = 0._r8
          allocate (tape(t)%hlist(f)%nacs(beg1d:end1d,numlev))
          tape(t)%hlist(f)%nacs(:,:) = 0
       end do
    end do
    
    if (masterproc) then
       write(6,*) 'Successfully initialized clm2 history files'
       write(6,'(72a1)') ("-",i=1,60)
    endif
    
    return
  end subroutine htapes_build

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: masterlist_change_active
!
! !INTERFACE:
  subroutine masterlist_change_active ()
!
! !DESCRIPTION: 
! Define default contents of history files. Call masterlist_make_active 
! for each field. Arguments are field name, tape index, and
! modification to averaging flag (blank means no modification)
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

    ! first tape (monthly output by default) (no default second tape)
    ! if all fields are kept at their default values then time avgeraging 
    ! flag does not have to be entered as a variable

    call masterlist_make_active (name='SNOWDP  ', tape_index=1)
    call masterlist_make_active (name='SNOWCOV ', tape_index=1)
    call masterlist_make_active (name='SNOWAGE ', tape_index=1)
    call masterlist_make_active (name='TSA     ', tape_index=1) 
    call masterlist_make_active (name='TCAS    ', tape_index=1) 
    call masterlist_make_active (name='TV      ', tape_index=1) 
    call masterlist_make_active (name='TG      ', tape_index=1) 
    call masterlist_make_active (name='TSNOW   ', tape_index=1) 
    call masterlist_make_active (name='TSOI    ', tape_index=1) 
    call masterlist_make_active (name='TLAKE   ', tape_index=1) 
    call masterlist_make_active (name='FSA     ', tape_index=1) 
    call masterlist_make_active (name='PARSUN  ', tape_index=1) 
    call masterlist_make_active (name='FSR     ', tape_index=1) 
    call masterlist_make_active (name='NDVI    ', tape_index=1) 
    call masterlist_make_active (name='FIRA    ', tape_index=1) 
    call masterlist_make_active (name='FIRE    ', tape_index=1) 
    call masterlist_make_active (name='FCTR    ', tape_index=1) 
    call masterlist_make_active (name='FCEV    ', tape_index=1) 
    call masterlist_make_active (name='FGEV    ', tape_index=1) 
    call masterlist_make_active (name='FSH     ', tape_index=1) 
    call masterlist_make_active (name='FGR     ', tape_index=1) 
    call masterlist_make_active (name='FSM     ', tape_index=1) 
    call masterlist_make_active (name='TAUX    ', tape_index=1) 
    call masterlist_make_active (name='TAUY    ', tape_index=1) 
    call masterlist_make_active (name='ELAI    ', tape_index=1) 
    call masterlist_make_active (name='ESAI    ', tape_index=1) 
    call masterlist_make_active (name='RSSUN   ', tape_index=1) 
    call masterlist_make_active (name='RSSHA   ', tape_index=1) 
    call masterlist_make_active (name='BTRAN   ', tape_index=1) 
    call masterlist_make_active (name='H2OSOI  ', tape_index=1) 
    call masterlist_make_active (name='H2OSNO  ', tape_index=1) 
    call masterlist_make_active (name='H2OCAN  ', tape_index=1) 
    call masterlist_make_active (name='SOILLIQ ', tape_index=1) 
    call masterlist_make_active (name='SOILICE ', tape_index=1) 
    call masterlist_make_active (name='SNOWLIQ ', tape_index=1) 
    call masterlist_make_active (name='SNOWICE ', tape_index=1) 
    call masterlist_make_active (name='QINFL   ', tape_index=1) 
    call masterlist_make_active (name='QOVER   ', tape_index=1) 
    call masterlist_make_active (name='QRGWL   ', tape_index=1) 
    call masterlist_make_active (name='QDRAI   ', tape_index=1) 
    call masterlist_make_active (name='QINTR   ', tape_index=1) 
    call masterlist_make_active (name='QDRIP   ', tape_index=1) 
    call masterlist_make_active (name='QMELT   ', tape_index=1) 
    call masterlist_make_active (name='QSOIL   ', tape_index=1) 
    call masterlist_make_active (name='QVEGE   ', tape_index=1) 
    call masterlist_make_active (name='QVEGT   ', tape_index=1) 
    call masterlist_make_active (name='PET     ', tape_index=1) 
#if (defined RTM)
    call masterlist_make_active (name='QCHANR  ', tape_index=1) 
    call masterlist_make_active (name='QCHOCNR ', tape_index=1) 
#endif
    call masterlist_make_active (name='ERRSOI  ', tape_index=1) 
    call masterlist_make_active (name='ERRSEB  ', tape_index=1) 
    call masterlist_make_active (name='ERRSOL  ', tape_index=1) 
    call masterlist_make_active (name='ERRH2O  ', tape_index=1) 
    call masterlist_make_active (name='RAIN    ', tape_index=1) 
    call masterlist_make_active (name='SNOW    ', tape_index=1) 
    call masterlist_make_active (name='QSTRM   ', tape_index=1) 
#if (defined STOCHASTIC)
    call masterlist_make_active (name='STRMDUR ', tape_index=1) 
    call masterlist_make_active (name='STRMDRY ', tape_index=1) 
    call masterlist_make_active (name='STORMS ', tape_index=1) 
#endif
    call masterlist_make_active (name='TBOT    ', tape_index=1) 
    call masterlist_make_active (name='THBOT   ', tape_index=1) 
    call masterlist_make_active (name='WIND    ', tape_index=1) 
    call masterlist_make_active (name='QBOT    ', tape_index=1) 
    call masterlist_make_active (name='ZBOT    ', tape_index=1) 
    call masterlist_make_active (name='FLDS    ', tape_index=1) 
    call masterlist_make_active (name='FSDS    ', tape_index=1) 

    return
  end subroutine masterlist_change_active
  
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: masterlist_make_active
!
! !INTERFACE:
  subroutine masterlist_make_active (name, tape_index, avgflag)
!
! !DESCRIPTION: 
! Add a field to the default "on" list for a given history file
! Also change the default time averaging flag if requested
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: name          ! field name
    integer, intent(in) :: tape_index             ! history tape index
    character(len=1), intent(in), optional :: avgflag  ! time averaging flag
!
! !REVISION HISTORY
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES
    integer :: f            ! field index
    logical :: found        ! flag indicates field found in masterlist
!-----------------------------------------------------------------------

    ! Check validity of input arguments

    if (tape_index > max_tapes) then
       write(6,*)'MASTERLIST_MAKE_ACTIVE: tape index=', tape_index, ' is too big'
       call endrun
    end if
    
    if (present(avgflag)) then
       if ( avgflag /= ' ' .and. &
            avgflag /= 'A' .and. avgflag /= 'I' .and. &
            avgflag /= 'X' .and. avgflag /= 'M') then
          write(6,*)'MASTERLIST_MAKE_ACTIVE: unknown averaging flag=', avgflag
          call endrun
       endif
    end if
    
    ! Look through master list for input field name.  
    ! When found, set active flag for that tape to true.  
    ! Also reset averaging flag if told to use other than default.
    
    found = .false.
    do f = 1,nfmaster
       if (trim(name) == trim(masterlist(f)%field%name)) then
          masterlist(f)%actflag(tape_index) = .true.
          if (present(avgflag)) then
             if (avgflag/= ' ') masterlist(f)%avgflag(tape_index) = avgflag
          end if
          found = .true.
          exit
       end if
    end do
    if (.not. found) then
       write(6,*)'MASTERLIST_MAKE_ACTIVE: field=', name, ' not found'
       call endrun
    end if
    
    return
  end subroutine masterlist_make_active
    
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: masterlist_change_timeavg
!
! !INTERFACE:
  subroutine masterlist_change_timeavg (t)
!
! !DESCRIPTION: 
! Override default history tape contents for a specific tape
! Copy the flag into the master field list
!
! !USES
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: t         ! history tape index
!
! !REVISION HISTORY
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES
    integer :: f                     ! field index
    character(len=1) :: avgflag      ! lcl equiv of hist_avgflag_pertape(t) 
!-----------------------------------------------------------------------

    avgflag = hist_avgflag_pertape(t)
    
    do f = 1,nfmaster
       select case (avgflag)
       case ('A')
          masterlist(f)%avgflag(t) = avgflag
       case ('I')
          masterlist(f)%avgflag(t) = avgflag
       case ('X')
          masterlist(f)%avgflag(t) = avgflag
       case ('M')
          masterlist(f)%avgflag(t) = avgflag
       case default
          write(6,*)'masterlist_change_timeavg: unknown avgflag=',avgflag
          call endrun ()
       end select
    end do
    
  end subroutine masterlist_change_timeavg
         
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: htapes_fieldlist
!
! !INTERFACE:
  subroutine htapes_fieldlist ()
!
! !DESCRIPTION: 
! Define the contents of each history file based on namelist
! input for initial or branch run, and restart data if a restart run.
! Use arrays fincl and fexcl to modify default history tape contents.
! Then sort the result alphanumerically.
!
! !USES
    use spmdMod, only : masterproc
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES
    integer :: t, f                ! tape, field indices
    integer :: ff                  ! index into include, exclude and fprec list
    character(len=max_namlen) :: name       ! field name portion of fincl (i.e. no avgflag separator)
    character(len=max_namlen) :: mastername ! name from masterlist field
    character(len=1)  :: avgflag    ! averaging flag
    character(len=1)  :: prec_acc   ! history buffer precision flag
    character(len=1)  :: prec_wrt   ! history buffer write precision flag
    type (history_entry) :: tmp     ! temporary used for swapping
!-----------------------------------------------------------------------

    ! First ensure contents of fincl and fexcl are valid names

    do t = 1,max_tapes
       f = 1
       do while (f < max_flds .and. fincl(f,t) /= ' ')
          name = getname (fincl(f,t))
          do ff = 1,nfmaster
             mastername = masterlist(ff)%field%name
             if (name == mastername) exit
          end do
          if (name /= mastername) then
             write(6,*)'HTAPES_FIELDLIST: ', trim(name), ' in fincl(', f, ') ',&
                  'for history tape ',t,' not found'
             call endrun
          end if
          f = f + 1
       end do
       
       f = 1
       do while (f < max_flds .and. fexcl(f,t) /= ' ')
          do ff = 1,nfmaster
             mastername = masterlist(ff)%field%name
             if (fexcl(f,t) == mastername) exit
          end do
          if (fexcl(f,t) /= mastername) then
             write(6,*)'HTAPES_FIELDLIST: ', fexcl(f,t), ' in fexcl(', f, ') ', &
                  'for history tape ',t,' not found'
             call endrun
          end if
          f = f + 1
       end do
    end do
    
    tape(:)%nflds = 0
    do t = 1,max_tapes

       ! Loop through the masterlist set of field names and determine if any of those
       ! are in the FINCL or FEXCL arrays
       ! The call to list_index determines the index in the FINCL or FEXCL arrays
       ! that the masterlist field corresponds to
       ! Add the field to the tape if specified via namelist (FINCL[1-max_tapes]), 
       ! or if it is on by default and was not excluded via namelist (FEXCL[1-max_tapes]).
       
       do f = 1,nfmaster
          mastername = masterlist(f)%field%name
          call list_index (fincl(1,t), mastername, ff)
          
          if (ff > 0) then
             
             ! if field is in include list, ff > 0 and htape_addfld
             ! will not be called for field

             avgflag = getflag (fincl(ff,t))
             call htape_addfld (t, f, avgflag)
             
          else if (.not. hist_empty_htapes) then
             
             ! find index of field in exclude list

             call list_index (fexcl(1,t), mastername, ff)

             ! if field is in exclude list, ff > 0 and htape_addfld
             ! will not be called for field
             ! if field is not in exclude list, ff =0 and htape_addfld
             ! will be called for field (note that htape_addfld will be
             ! called below only if field is not in exclude list OR in
             ! include list

             if (ff == 0 .and. masterlist(f)%actflag(t)) then
                call htape_addfld (t, f, ' ')
             end if
             
          end if
       end do
       
       ! Specification of tape contents now complete.  
       ! Sort each list of active entries 
       
       do f = tape(t)%nflds-1,1,-1
          do ff = 1,f
             if (tape(t)%hlist(ff)%field%name > tape(t)%hlist(ff+1)%field%name) then
                
                tmp = tape(t)%hlist(ff)
                tape(t)%hlist(ff  ) = tape(t)%hlist(ff+1)
                tape(t)%hlist(ff+1) = tmp
                
             else if (tape(t)%hlist(ff  )%field%name == tape(t)%hlist(ff+1)%field%name) then
                
                write(6,*)'HTAPES_FIELDLIST: Duplicate field ', tape(t)%hlist(ff  )%field%name
                write(6,*)'t,ff,name=',t,ff,tape(t)%hlist(ff  )%field%name
                call endrun
                
             end if
          end do
       end do
       
       if (masterproc) then
          if (tape(t)%nflds > 0) then
             write(6,*)'HTAPES_FIELDLIST: Included fields tape ',t,'=',tape(t)%nflds
          end if
          do f = 1,tape(t)%nflds
             write(6,*) f,' ',tape(t)%hlist(f)%field%name, &
                  tape(t)%hlist(f)%field%numlev,' ',tape(t)%hlist(f)%avgflag
          end do
       end if
    end do
    
    ! Determine total number of active history tapes
    
    ntapes = 0
    do t = max_tapes,1,-1
       if (tape(t)%nflds > 0) then
          ntapes = t
          exit
       end if
    end do
    
    ! Ensure there are no "holes" in tape specification, i.e. empty tapes.
    ! Enabling holes should not be difficult if necessary.
    
    do t = 1,ntapes
       if (tape(t)%nflds  ==  0) then
          write(6,*)'HTAPES_FIELDLIST: Tape ',t,' is empty'
          call endrun
       end if
    end do
    
    ! Check that the number of history files declared does not exceed
    ! the maximum allowed.
    
    if (ntapes > max_tapes) then
       write(6,*) 'HTAPES_FIELDLIST: Too many history files declared, max=',max_tapes
       write(6,*)'To increase, change parameter max_tapes.'
       call endrun
    end if

    ! Change 1d output per tape output flag if requested - only for history
    ! tapes greater than 1 where 2d xy averaging is not enabled

    do t = 2,ntapes
       if (hist_type1d_pertape(t) /= ' ' .and. (.not. hist_dov2xy(t))) then
          select case (trim(hist_type1d_pertape(t)))
          case ('PFTS','COLS', 'LAND', 'GRID')
             write(6,*)'history tape ',t,' will have 1d type output of ',hist_type1d_pertape(t)
          case default
             write(6,*)'HTAPES_FIELDLIST: unknown namelist type1d per tape=',hist_type1d_pertape(t)
             call endrun ()
          end select
       end if
    end do

    if (masterproc) then
       write(6,*) 'There will be a total of ',ntapes,' history tapes'
       do t=1,ntapes
          write(6,*)
          if (hist_nhtfrq(t) == 0) then
             write(6,*)'History tape ',t,' write frequency is MONTHLY'
          else
             write(6,*)'History tape ',t,' write frequency = ',hist_nhtfrq(t)
          endif
          if (hist_dov2xy(t)) then 
             write(6,*)'All fields on history tape ',t,' are grid averaged'
          else
             write(6,*)'All fields on history tape ',t,' are not grid averaged'
          end if
          write(6,*)'Number of time samples on history tape ',t,' is ',hist_mfilt(t)
          write(6,*)'Output precision on history tape ',t,'=',hist_ndens(t)
          write(6,*)
       end do
    end if
    
    ! Set flag indicating h-tape contents are now defined (needed by masterlist_addfld)
    
    htapes_defined = .true.      
    
    return
  end subroutine htapes_fieldlist
  
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: htape_addfld
!
! !INTERFACE:
  subroutine htape_addfld (t, f, avgflag)
!
! !DESCRIPTION: 
! Add a field to the active list for a history tape. Copy the data from 
! the master field list to the active list for the tape.
!
! !USES
    use clmtype, only : grid1d, land1d, cols1d, pfts1d
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: t                 ! history tape index
    integer, intent(in) :: f                 ! field index from master field list
    character(len=1), intent(in) :: avgflag  ! time averaging flag
!
! !REVISION HISTORY
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES
    integer :: n                        ! field index on defined tape
    character(len=4) :: hist_type1d     ! 1d type ('GRID','LAND','COLS','PFTS')
    integer :: index                    ! clmtype pointer index
!-----------------------------------------------------------------------

    ! Ensure that it is not to late to add a field to the history tape

    if (htapes_defined) then
       write(6,*)'HTAPE_ADDFLD: Attempt to add field ',masterlist(f)%field%name,&
            ' after history files set'
       call endrun ()
    end if
    
    tape(t)%nflds = tape(t)%nflds + 1
    n = tape(t)%nflds
    
    ! Copy field information
    
    tape(t)%hlist(n)%field = masterlist(f)%field
    
    ! Reset 1dtype field and field pointer based on namelist hist_type1d_pertape
    ! Only applies to tapes other than primary and when 2d xy averaging is not active

    if (t > 1 .and. hist_type1d_pertape(t) /= ' ' .and. (.not. hist_dov2xy(t))) then
       select case (trim(hist_type1d_pertape(t)))
       case('GRID')
          index = tape(t)%hlist(n)%field%hpindices(1)
          if (index == -1) then
             write(6,*)'HTAPE_ADDFLD: the field ',tape(t)%hlist(n)%field%name, &
                  ' was chosen incorrectly for gridcell 1d output'
             call endrun
          endif
          tape(t)%hlist(n)%field%hpindex = index
          tape(t)%hlist(n)%field%beg1d = grid1d%beg
          tape(t)%hlist(n)%field%end1d = grid1d%end
          tape(t)%hlist(n)%field%num1d = grid1d%num
          tape(t)%hlist(n)%field%type1d = 'gridcell'
       case('LAND')
          index = tape(t)%hlist(n)%field%hpindices(2)
          if (index == -1) then
             write(6,*)'HTAPE_ADDFLD: the field ',tape(t)%hlist(n)%field%name, &
                  ' was chosen incorrectly for landunit 1d output'
             call endrun
          endif
          tape(t)%hlist(n)%field%hpindex = index
          tape(t)%hlist(n)%field%beg1d = land1d%beg
          tape(t)%hlist(n)%field%end1d = land1d%end
          tape(t)%hlist(n)%field%num1d = land1d%num
          tape(t)%hlist(n)%field%type1d = 'landunit'
       case('COLS')
          index = tape(t)%hlist(n)%field%hpindices(3)
          if (index == -1) then
             write(6,*)'HTAPE_ADDFLD: the field ',tape(t)%hlist(n)%field%name, &
                  ' was chosen incorrectly for column 1d output'
             call endrun
          endif
          tape(t)%hlist(n)%field%hpindex = index
          tape(t)%hlist(n)%field%beg1d = cols1d%beg
          tape(t)%hlist(n)%field%end1d = cols1d%end
          tape(t)%hlist(n)%field%num1d = cols1d%num
          tape(t)%hlist(n)%field%type1d = 'column'
       case ('PFTS')
          index = tape(t)%hlist(n)%field%hpindices(4)
          if (index == -1) then
             write(6,*)'HTAPE_ADDFLD: the field ',tape(t)%hlist(n)%field%name, &
                  ' was chosen incorrectly for pft 1d output'
             call endrun
          endif
          tape(t)%hlist(n)%field%hpindex = index
          tape(t)%hlist(n)%field%beg1d = pfts1d%beg
          tape(t)%hlist(n)%field%end1d = pfts1d%end
          tape(t)%hlist(n)%field%num1d = pfts1d%num
          tape(t)%hlist(n)%field%type1d = 'pft'
       case default
          write(6,*)'HTAPE_ADDFLD: unknown namelist input hist_type1d_pertape=', &
               hist_type1d_pertape(t)
          call endrun
       end select
    endif

    ! Override the default averaging (masterlist) averaging flag if non-blank

    select case (avgflag)
    case (' ')
       tape(t)%hlist(n)%avgflag = masterlist(f)%avgflag(t)
    case ('A','I','X','M') 
       tape(t)%hlist(n)%avgflag = avgflag
    case default
       write(6,*)'HTAPE_ADDFLD: unknown avgflag=', avgflag
       call endrun
    end select
    
#if (defined DEBUG)
    write(6,*)'HTAPE_ADDFLD: field ',tape(t)%hlist(n)%field%name,' added as field number ',n,' on tape ',t
    write(6,*)'units      =',tape(t)%hlist(n)%field%units
    write(6,*)'numlev     =',tape(t)%hlist(n)%field%numlev
    write(6,*)'avgflag    =',tape(t)%hlist(n)%avgflag
#endif
    
    return
  end subroutine htape_addfld

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: update_hbuf
!
! !INTERFACE:
  subroutine update_hbuf ()
!
! !DESCRIPTION: 
! Accumulate (or take min, max, etc. as appropriate) input field
! into its history buffer for appropriate tapes
!
! !USES
    use clmpoint
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES
    integer :: k                   ! 1d index
    integer :: lev                 ! level index 
    integer :: t                   ! tape index
    integer :: f                   ! field index
    integer :: hpindex             ! history pointer index   
    integer :: beg1d,end1d         ! 1d beginning and ending indices
    integer :: numlev              ! number of vertical levels (dimension and loop)
    character(len=1)  :: avgflag   ! time averaging flag
    real(r8), pointer :: hbuf(:,:) ! history buffer
    integer , pointer :: nacs(:,:) ! accumulation counter
!-----------------------------------------------------------------------
    
    do t = 1,ntapes
       
       do f = 1,tape(t)%nflds
       
          hpindex = tape(t)%hlist(f)%field%hpindex

          ! perform time averaging

          avgflag = tape(t)%hlist(f)%avgflag
          beg1d   = tape(t)%hlist(f)%field%beg1d
          end1d   = tape(t)%hlist(f)%field%end1d
          numlev  = tape(t)%hlist(f)%field%numlev
          nacs   => tape(t)%hlist(f)%nacs
          hbuf   => tape(t)%hlist(f)%hbuf
       
          ! if field is not defined for index type, 
          ! just set to spval and continue to next field

          if (hpindex == not_valid) then
             hbuf(:,:) = spval
             cycle
          endif
       
          select case (avgflag)
          case ('I') ! Instantaneous
             
             if (numlev == 1) then
!$OMP PARALLEL DO PRIVATE (k)
                do k = beg1d,end1d
                   if (clmpointer(hpindex)%val(k)%rp /= spval) then 
                      hbuf(k,1) = clmpointer(hpindex)%val(k)%rp
                   else
                      hbuf(k,1) = spval
                   endif
                   nacs(k,1) = 1
                end do
             else 
                do lev = 1,numlev
!$OMP PARALLEL DO PRIVATE (k)
                   do k = beg1d,end1d
                      if (clmpointer(hpindex)%val(k)%rap(lev) /= spval) then 
                         hbuf(k,lev) = clmpointer(hpindex)%val(k)%rap(lev)
                      else
                         hbuf(k,lev) = spval
                      endif
                      nacs(k,lev) = 1
                   end do
                end do
             endif
             
          case ('A') ! Time average
             
             if (numlev == 1) then
!$OMP PARALLEL DO PRIVATE (k)
                do k = beg1d,end1d
                   if (clmpointer(hpindex)%val(k)%rp /= spval) then 
                      if (nacs(k,1) == 0) hbuf(k,1) = 0.
                      hbuf(k,1) = hbuf(k,1) + clmpointer(hpindex)%val(k)%rp
                      nacs(k,1) = nacs(k,1) + 1
                   else
                      if (nacs(k,1) == 0) hbuf(k,1) = spval
                   endif
                end do
             else 
                do lev = 1,numlev
!$OMP PARALLEL DO PRIVATE (k)
                   do k = beg1d,end1d
                      if (clmpointer(hpindex)%val(k)%rap(lev) /= spval) then 
                         if (nacs(k,lev) == 0) hbuf(k,lev) = 0.
                         hbuf(k,lev) = hbuf(k,lev) + clmpointer(hpindex)%val(k)%rap(lev)
                         nacs(k,lev) = nacs(k,lev) + 1
                      else
                         if (nacs(k,lev) == 0) hbuf(k,lev) = spval
                      endif
                   end do
                end do
             endif
             
          case ('X') ! Maximum over time
                
             if (numlev == 1) then
!$OMP PARALLEL DO PRIVATE (k)
                do k = beg1d,end1d
                   if (clmpointer(hpindex)%val(k)%rp /= spval) then 
                      if (nacs(k,1) == 0) hbuf(k,1) = -1.e50
                      hbuf(k,1) = max( hbuf(k,1), clmpointer(hpindex)%val(k)%rp )
                   else
                      hbuf(k,1) = spval
                   endif
                   nacs(k,1) = 1
                end do
             else 
                do lev = 1,numlev
!$OMP PARALLEL DO PRIVATE (k)
                   do k = beg1d,end1d
                      if (clmpointer(hpindex)%val(k)%rap(lev) /= spval) then 
                         if (nacs(k,lev) == 0) hbuf(k,lev) = -1.e50
                         hbuf(k,lev) = max( hbuf(k,lev), clmpointer(hpindex)%val(k)%rap(lev) )
                      else
                         hbuf(k,lev) = spval
                      end if
                      nacs(k,lev) = 1
                   end do
                end do
             endif
             
          case ('M') ! Minimum over time
                
             if (numlev == 1) then
!$OMP PARALLEL DO PRIVATE (k)
                do k = beg1d,end1d
                   if (clmpointer(hpindex)%val(k)%rp /= spval) then 
                      if (nacs(k,1) == 0) hbuf(k,1) = +1.e50
                      hbuf(k,1) = min( hbuf(k,1), clmpointer(hpindex)%val(k)%rp )
                   else
                      hbuf(k,1) = spval
                   endif
                   nacs(k,1) = 1
                end do
             else
                do lev = 1,numlev
!$OMP PARALLEL DO PRIVATE (k)
                   do k = beg1d,end1d
                      if (clmpointer(hpindex)%val(k)%rap(lev) /= spval) then 
                         if (nacs(k,lev) == 0) hbuf(k,lev)= +1.e50
                         hbuf(k,lev) = min( hbuf(k,lev), clmpointer(hpindex)%val(k)%rap(lev) )
                      else
                         hbuf(k,lev) = spval
                      endif
                      nacs(k,lev) = 1
                   end do
                end do
             endif
             
          case default
             
             write(6,*)'UPDATE_HBUF: invalid time averaging flag ', avgflag
             call endrun ()
             
          end select

       end do   ! end of field loop 
    end do   ! end of tape loop
          
    return
  end subroutine update_hbuf
     
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: hfields_normalize
!
! !INTERFACE:
  subroutine hfields_normalize (t)
!
! !DESCRIPTION: 
! Normalize fields on a history file by the number of accumulations
! Loop over fields on the tape.  Need averaging flag and number of
! accumulations to perform normalization.
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: t       ! tape index
!
! !REVISION HISTORY
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES
    integer :: f                   ! field index
    integer :: k                   ! 1d index
    integer :: lev                 ! level index
    integer :: beg1d,end1d         ! 1d beginning and ending indices 
    integer :: numlev              ! number of history buffer levels 
    character(len=1)  :: avgflag   ! averaging flag
    real(r8), pointer :: hbuf(:,:) ! history buffer
    integer , pointer :: nacs(:,:) ! accumulation counter
!-----------------------------------------------------------------------
    
    ! Normalize by number of accumulations for time averaged case

    do f = 1,tape(t)%nflds
       avgflag = tape(t)%hlist(f)%avgflag
       beg1d   = tape(t)%hlist(f)%field%beg1d
       end1d   = tape(t)%hlist(f)%field%end1d
       numlev  = tape(t)%hlist(f)%field%numlev
       nacs   => tape(t)%hlist(f)%nacs
       hbuf   => tape(t)%hlist(f)%hbuf

       do lev = 1, numlev
          do k = beg1d, end1d
             if (avgflag == 'A' .and. nacs(k,lev) /= 0) then
                hbuf(k,lev) = hbuf(k,lev) / float(nacs(k,lev))
             end if
          end do
       end do
    end do
        
    return
  end subroutine hfields_normalize

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: hfields_zero
!
! !INTERFACE:
  subroutine hfields_zero (t)
!
! !DESCRIPTION: 
! Zero out accumulation and history buffers for a given history tape.
! Loop through fields on the tape.
!
! !USES
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: t     ! tape index
!
! !REVISION HISTORY
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES
    integer :: f                 ! field index
!-----------------------------------------------------------------------

    do f = 1,tape(t)%nflds
       tape(t)%hlist(f)%hbuf(:,:) = 0.
       tape(t)%hlist(f)%nacs(:,:) = 0
    end do
    
  end subroutine hfields_zero

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: htape_create
!
! !INTERFACE:
  subroutine htape_create (t)
!
! !DESCRIPTION: 
! Define contents of history file t. Issue the required netcdf 
! wrapper calls to define the history file contents.
!
! !USES
    use clm_varpar, only : lsmlon, lsmlat, nlevsoi
    use clmtype, only : grid1d, land1d, cols1d, pfts1d
    use clm_varctl, only : caseid, ctitle, frivinp_rtm, fsurdat, finidat, fpftcon
    use clm_varsur, only : fullgrid, offline_rdgrid, longxy, latixy, area, &
         landfrac, landmask, numlon, lsmedge
    use clm_varcon, only : zsoi
    use fileutils, only : get_filename
    use time_manager, only : get_ref_date	
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: t   ! tape index
!
! !REVISION HISTORY
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES
    integer :: f                   ! field index
    integer :: yr,mon,day,nbsec    ! year,month,day,seconds components of a date
    integer :: hours,minutes,secs  ! hours,minutes,seconds of hh:mm:ss
    integer :: numlev              ! number of vertical levels (dimension and loop)
    integer :: omode               ! returned mode from netCDF call
    integer :: ncprec              ! output netCDF write precision
    integer :: dim1_id(1)          ! netCDF dimension id for 1-d variables
    integer :: dim2_id(2)          ! netCDF dimension id for 2-d variables
    integer :: dim3_id(3)          ! netCDF dimension id for 3-d variables
    integer :: dim4_id(4)          ! netCDF dimension id for 4-d variables
    integer :: ret                 ! netCDF error status   

    integer :: grid1d_lon_id(max_tapes)    ! netCDF variable id
    integer :: grid1d_lat_id(max_tapes)    ! netCDF variable id
    integer :: grid1d_ixy_id(max_tapes)    ! netCDF variable id
    integer :: grid1d_jxy_id(max_tapes)    ! netCDF variable id

    integer :: land1d_lon_id(max_tapes)    ! netCDF variable id
    integer :: land1d_lat_id(max_tapes)    ! netCDF variable id
    integer :: land1d_ixy_id(max_tapes)    ! netCDF variable id
    integer :: land1d_jxy_id(max_tapes)    ! netCDF variable id
    integer :: land1d_gindex_id(max_tapes) ! netCDF variable id
    integer :: land1d_wtxy_id(max_tapes)   ! netCDF variable id
    integer :: land1d_itypwat_id(max_tapes)! netCDF variable id

    integer :: cols1d_lon_id(max_tapes)    ! netCDF variable id
    integer :: cols1d_lat_id(max_tapes)    ! netCDF variable id
    integer :: cols1d_ixy_id(max_tapes)    ! netCDF variable id
    integer :: cols1d_jxy_id(max_tapes)    ! netCDF variable id
    integer :: cols1d_gindex_id(max_tapes) ! netCDF variable id
    integer :: cols1d_lindex_id(max_tapes) ! netCDF variable id
    integer :: cols1d_wtxy_id(max_tapes)   ! netCDF variable id
    integer :: cols1d_wtlnd_id(max_tapes)  ! netCDF variable id
    integer :: cols1d_itypwat_id(max_tapes)! netCDF variable id

    integer :: pfts1d_lon_id(max_tapes)    ! netCDF variable id
    integer :: pfts1d_lat_id(max_tapes)    ! netCDF variable id
    integer :: pfts1d_ixy_id(max_tapes)    ! netCDF variable id
    integer :: pfts1d_jxy_id(max_tapes)    ! netCDF variable id
    integer :: pfts1d_gindex_id(max_tapes) ! netCDF variable id
    integer :: pfts1d_lindex_id(max_tapes) ! netCDF variable id
    integer :: pfts1d_cindex_id(max_tapes) ! netCDF variable id
    integer :: pfts1d_wtxy_id(max_tapes)   ! netCDF variable id
    integer :: pfts1d_wtlnd_id(max_tapes)  ! netCDF variable id
    integer :: pfts1d_wtcol_id(max_tapes)  ! netCDF variable id
    integer :: pfts1d_itypwat_id(max_tapes)! netCDF variable id
    integer :: pfts1d_itypveg_id(max_tapes)! netCDF variable id

    character(len=  8) :: curdate  ! current date
    character(len=  8) :: curtime  ! current time 
    character(len= 10) :: basedate ! base date (yyyymmdd)
    character(len=  8) :: basesec  ! base seconds
    character(len=256) :: name     ! name of attribute
    character(len=256) :: unit     ! units of attribute
    character(len=256) :: str      ! global attribute string 
    character(len=  1) :: avgflag  ! time averaging flag 
    real(r8) :: lonvar(lsmlon)     ! only used for full grid 
    real(r8) :: latvar(lsmlat)     ! only used for full grid
!-----------------------------------------------------------------------

    ! Create new netCDF file. File will be in define mode
    
    write(6,*)'Opening netcdf htape ', trim(locfnh(t))
    call wrap_create (trim(locfnh(t)), nf_clobber, nfid(t))
    
    ret = nf_set_fill (nfid(t), nf_nofill, omode)
    
    ! define output write precsion for tape
    
    ncprec = tape(t)%ncprec
    
    ! Create global attributes. Attributes are used to store information
    ! about the data set. Global attributes are information about the
    ! data set as a whole, as opposed to a single variable

    call wrap_put_att_text (nfid(t), NF_GLOBAL, 'Conventions', 'CF-1.0')
    
    call getdatetime(curdate, curtime)
    str = 'created on ' // curdate // ' ' // curtime
    call wrap_put_att_text(nfid(t), NF_GLOBAL, 'history', trim(str))
    
    call getenv ('LOGNAME', str)
    call wrap_put_att_text (nfid(t), NF_GLOBAL, 'logname',trim(str))
    
    call getenv ('HOST', str)
    call wrap_put_att_text (nfid(t), NF_GLOBAL, 'host', trim(str))
    
    str = 'Community Land Model: CLM2'
    call wrap_put_att_text (nfid(t), NF_GLOBAL, 'source', trim(str))
     
    str = '$Name$'
    call wrap_put_att_text (nfid(t), NF_GLOBAL, 'version', trim(str))
     
    str = '$Id$'
    call wrap_put_att_text (nfid(t), NF_GLOBAL, 'revision_id', trim(str))

    str = ctitle 
    call wrap_put_att_text (nfid(t), NF_GLOBAL, 'case_title', trim(str))

    str = caseid
    call wrap_put_att_text (nfid(t), NF_GLOBAL, 'case_id', trim(str))

    if (fsurdat == ' ') then
       str = 'created at run time'
    else
       str = get_filename(fsurdat)
    endif
    call wrap_put_att_text(nfid(t), NF_GLOBAL, 'Surface_dataset', trim(str))

    if (finidat == ' ') then
       str = 'arbitrary initialization'
    else
       str = get_filename(finidat)
    endif
    call wrap_put_att_text(nfid(t), NF_GLOBAL, 'Initial_conditions_dataset', trim(str))

    str = get_filename(fpftcon)
    call wrap_put_att_text(nfid(t), NF_GLOBAL, 'PFT_physiological_constants_dataset', trim(str))

    if (frivinp_rtm /= ' ') then
       str = get_filename(frivinp_rtm)
       call wrap_put_att_text(nfid(t), NF_GLOBAL, 'RTM_input_datset', trim(str))
    endif
     
    ! Define dimensions. Array dimensions are referenced by an
    ! associated dimenision id: e.g., lon_id -> lon.
    ! Time is an unlimited dimension.
    ! Character string is treated as an array of characters. 
     
     if (.not. tape(t)%dov2xy) then
        call wrap_def_dim (nfid(t), 'gridcell', grid1d%num, gridcell_dimid)
        call wrap_def_dim (nfid(t), 'landunit', land1d%num, landunit_dimid)
        call wrap_def_dim (nfid(t), 'column', cols1d%num, column_dimid)
        call wrap_def_dim (nfid(t), 'pft', pfts1d%num, pft_dimid)
     end if
     call wrap_def_dim (nfid(t), 'lon', lsmlon, lon_dimid)
     call wrap_def_dim (nfid(t), 'lat', lsmlat, lat_dimid)
     call wrap_def_dim (nfid(t), 'levsoi', nlevsoi, levsoi_dimid)
     call wrap_def_dim (nfid(t), 'time', nf_unlimited, time_dimid)
     call wrap_def_dim (nfid(t), 'hist_interval', 2, hist_interval_dimid)
     call wrap_def_dim (nfid(t), 'string_length', 8, strlen_dimid)

    ! Define time-independent grid variables 
    ! Coordinate variables (including time)

     if (fullgrid) then
        dim1_id(1) = lon_dimid
        call wrap_def_var (nfid(t), 'lon' , ncprec, 1, dim1_id, lonvar_id(t))
        call wrap_put_att_text (nfid(t), lonvar_id(t), 'long_name',&
             'coordinate longitude')
        call wrap_put_att_text (nfid(t), lonvar_id(t), 'units'    ,&
             'degrees_east')
        
        dim1_id(1) = lat_dimid
        call wrap_def_var (nfid(t), 'lat' , ncprec, 1, dim1_id, latvar_id(t))
        call wrap_put_att_text (nfid(t), latvar_id(t), 'long_name',&
             'coordinate latitude')
        call wrap_put_att_text (nfid(t), latvar_id(t), 'units'    ,&
             'degrees_north')
     endif
     
     dim1_id(1) = levsoi_dimid
     call wrap_def_var (nfid(t), 'levsoi' , ncprec, 1, dim1_id, levvar_id(t))
     call wrap_put_att_text (nfid(t), levvar_id(t), 'long_name',&
          'coordinate soil levels')
     call wrap_put_att_text (nfid(t), levvar_id(t), 'units'    ,&
          'm')
     
     dim1_id(1) = time_dimid
     call get_ref_date(yr, mon, day, nbsec)
     hours   = nbsec / 3600
     minutes = (nbsec - hours*3600) / 60
     secs    = (nbsec - hours*3600 - minutes*60)
     write(basedate,80) yr,mon,day
80   format(i4.4,'-',i2.2,'-',i2.2)
     write(basesec ,90) hours, minutes, secs
90   format(i2.2,':',i2.2,':',i2.2)
     unit = 'days since ' // basedate // " " // basesec
     call wrap_def_var (nfid(t), 'time', ncprec, 1, dim1_id, time_var_id(t))
     call wrap_put_att_text (nfid(t), time_var_id(t), 'long_name','time')
     call wrap_put_att_text (nfid(t), time_var_id(t), 'units'    ,unit)
     call wrap_put_att_text (nfid(t), time_var_id(t), 'calendar' ,'noleap')
     
#if (defined OFFLINE)
     if (.not. offline_rdgrid) then
        call wrap_def_var (nfid(t), 'edgen', ncprec, 0, 0, edgen_id(t))
        call wrap_put_att_text (nfid(t), edgen_id(t), 'long_name',&
             'northern edge of surface grid')
        call wrap_put_att_text (nfid(t), edgen_id(t), 'units'    ,&
             'degrees_north')
       
        call wrap_def_var (nfid(t), 'edgee', ncprec, 0, 0, edgee_id(t))
        call wrap_put_att_text (nfid(t), edgee_id(t), 'long_name',&
             'eastern edge of surface grid')
        call wrap_put_att_text (nfid(t), edgee_id(t), 'units'    ,&
             'degrees_east')
        
        call wrap_def_var (nfid(t), 'edges', ncprec, 0, 0, edges_id(t))
        call wrap_put_att_text (nfid(t), edges_id(t), 'long_name',&
             'southern edge of surface grid')
        call wrap_put_att_text (nfid(t), edges_id(t), 'units'    ,&
             'degrees_north')
        
        call wrap_def_var (nfid(t), 'edgew', ncprec, 0, 0, edgew_id(t))
        call wrap_put_att_text (nfid(t), edgew_id(t) , 'long_name',&
             'western edge of surface grid')
        call wrap_put_att_text (nfid(t), edgew_id(t) , 'units'    ,&
             'degrees_east')
     endif
#endif

     ! Longitude, latitude, surface type: real (lsmlon x lsmlat)

     dim2_id(1) = lon_dimid
     dim2_id(2) = lat_dimid
     
     if (fullgrid) then
        name = 'longitude'
        unit = 'degrees_east'
        call wrap_def_var (nfid(t), 'longxy' , ncprec, 2, dim2_id, longxy_id(t))
     else
        name = 'rlongitude'
        unit = 'degrees_east'
        call wrap_def_var (nfid(t), 'rlongxy', ncprec, 2, dim2_id, longxy_id(t))
     endif
     call wrap_put_att_text (nfid(t), longxy_id(t), 'long_name',name)
     call wrap_put_att_text (nfid(t), longxy_id(t), 'units'    ,unit)
     
     call wrap_def_var (nfid(t), 'latixy', ncprec, 2, dim2_id, latixy_id(t))
     call wrap_put_att_text (nfid(t), latixy_id(t), 'long_name',&
          'latitude')
     call wrap_put_att_text (nfid(t), latixy_id(t), 'units'    ,&
          'degrees_north')
     
     call wrap_def_var (nfid(t), 'area', ncprec, 2, dim2_id, area_id(t))
     call wrap_put_att_text (nfid(t), area_id(t), 'long_name',&
          'grid cell areas')
     call wrap_put_att_text (nfid(t), area_id(t), 'units'    ,&
          'km^2')
     
     call wrap_def_var (nfid(t), 'landfrac', ncprec, 2, dim2_id, landfrac_id(t))
     call wrap_put_att_text (nfid(t), landfrac_id(t), 'long_name',&
          'land fraction')
     
     ! Number of longitudes per latitude (reduced grid only)
     
     dim1_id(1) = lat_dimid
     call wrap_def_var (nfid(t), 'numlon', nf_int, 1, dim1_id, numlon_id(t))
     call wrap_put_att_text (nfid(t), numlon_id(t), 'long_name', &
          'number of longitudes at each latitude')
     
     ! Surface type
     
     call wrap_def_var (nfid(t), 'landmask', nf_int, 2, dim2_id, landmask_id(t))
     call wrap_put_att_text (nfid(t), landmask_id(t),'long_name',&
          'land/ocean mask (0.=ocean and 1.=land)')
     
     ! Define time-independent 1d-indices surrounding datatypes

     if (.not. tape(t)%dov2xy) then

        ! 1d gridcell info

        dim1_id(1) = gridcell_dimid

        call wrap_def_var(nfid(t), 'grid1d_lon', ncprec, 1, dim1_id, grid1d_lon_id(t))
        call wrap_put_att_text (nfid(t), grid1d_lon_id(t), 'long_name', &
             '1d gridcell longitude')
        call wrap_put_att_text (nfid(t), grid1d_lon_id(t), 'units',&
             'degrees_east')

        call wrap_def_var(nfid(t), 'grid1d_lat', ncprec, 1, dim1_id, grid1d_lat_id(t))
        call wrap_put_att_text (nfid(t), grid1d_lat_id(t), 'long_name', &
             '1d gridcell latitude')
        call wrap_put_att_text (nfid(t), grid1d_lat_id(t), 'units',&
             'degrees_north')

        call wrap_def_var(nfid(t), 'grid1d_ixy', nf_int, 1, dim1_id, grid1d_ixy_id(t))
        call wrap_put_att_text (nfid(t), grid1d_ixy_id(t), 'long_name', &
             '2d longitude index of corresponding gridcell')

        call wrap_def_var(nfid(t), 'grid1d_jxy', nf_int, 1, dim1_id, grid1d_jxy_id(t))
        call wrap_put_att_text (nfid(t), grid1d_jxy_id(t), 'long_name', &
             '2d latitude index of corresponding gridcell')

        ! 1d landunit info

        dim1_id(1) = landunit_dimid

        call wrap_def_var(nfid(t), 'land1d_lon', ncprec, 1, dim1_id, land1d_lon_id(t))
        call wrap_put_att_text (nfid(t), land1d_lon_id(t), 'long_name', &
             '1d landunit longitude')
        call wrap_put_att_text (nfid(t), land1d_lon_id(t), 'units', &
             'degrees_east')

        call wrap_def_var(nfid(t), 'land1d_lat', ncprec, 1, dim1_id, land1d_lat_id(t))
        call wrap_put_att_text (nfid(t), land1d_lat_id(t), 'long_name', &
             '1d landunit latitude')
        call wrap_put_att_text (nfid(t), land1d_lat_id(t), 'units',&
             'degrees_north')

        call wrap_def_var(nfid(t), 'land1d_ixy', nf_int, 1, dim1_id, land1d_ixy_id(t))
        call wrap_put_att_text (nfid(t), land1d_ixy_id(t), 'long_name', &
             '2d longitude index of corresponding landunit')

        call wrap_def_var(nfid(t), 'land1d_jxy', nf_int, 1, dim1_id, land1d_jxy_id(t))
        call wrap_put_att_text (nfid(t), land1d_jxy_id(t), 'long_name', &
             '2d latitude index of corresponding landunit')

        call wrap_def_var(nfid(t), 'land1d_gi', nf_int, 1, dim1_id, land1d_gindex_id(t))
        call wrap_put_att_text (nfid(t), land1d_gindex_id(t) , 'long_name', &
             '1d grid index of corresponding landunit')

        call wrap_def_var(nfid(t), 'land1d_wtxy', ncprec, 1, dim1_id, land1d_wtxy_id(t))
        call wrap_put_att_text (nfid(t), land1d_wtxy_id(t) , 'long_name', &
             'landunit weight relative to corresponding gridcell')

        call wrap_def_var(nfid(t), 'land1d_itypwat', nf_int, 1, dim1_id, land1d_itypwat_id(t))
        call wrap_put_att_text (nfid(t), land1d_itypwat_id(t) , 'long_name', &
             'landunit water type')

        ! 1d column info

        dim1_id(1) = column_dimid

        call wrap_def_var(nfid(t), 'cols1d_lon', ncprec, 1, dim1_id, cols1d_lon_id(t))
        call wrap_put_att_text (nfid(t), cols1d_lon_id(t), 'long_name', &
             '1d column longitude')
        call wrap_put_att_text (nfid(t), cols1d_lon_id(t), 'units', &
             'degrees_east')

        call wrap_def_var(nfid(t), 'cols1d_lat', ncprec, 1, dim1_id, cols1d_lat_id(t))
        call wrap_put_att_text (nfid(t), cols1d_lat_id(t), 'long_name', &
             '1d column latitude')
        call wrap_put_att_text (nfid(t), cols1d_lat_id(t), 'units', &
             'degrees_north')

        call wrap_def_var(nfid(t), 'cols1d_ixy', nf_int, 1, dim1_id, cols1d_ixy_id(t))
        call wrap_put_att_text (nfid(t), cols1d_ixy_id(t), 'long_name', &
             '2d longitude index of corresponding landunit')

        call wrap_def_var(nfid(t), 'cols1d_jxy', nf_int, 1, dim1_id, cols1d_jxy_id(t))
        call wrap_put_att_text (nfid(t), cols1d_jxy_id(t), 'long_name', &
             '2d latitude index of corresponding landunit')

        call wrap_def_var(nfid(t), 'cols1d_gi ', nf_int, 1, dim1_id, cols1d_gindex_id(t))
        call wrap_put_att_text (nfid(t), cols1d_gindex_id(t) , 'long_name', &
             '1d grid index of corresponding column')

        call wrap_def_var(nfid(t), 'cols1d_li ', nf_int, 1, dim1_id, cols1d_lindex_id(t))
        call wrap_put_att_text (nfid(t), cols1d_lindex_id(t) , 'long_name', &
             '1d landunit index of corresponding column')

        call wrap_def_var(nfid(t), 'cols1d_wtxy', ncprec, 1, dim1_id, cols1d_wtxy_id(t))
        call wrap_put_att_text (nfid(t), cols1d_wtxy_id(t) , 'long_name', &
             'column weight relative to corresponding gridcell')

        call wrap_def_var(nfid(t), 'cols1d_wtlnd', ncprec, 1, dim1_id, cols1d_wtlnd_id(t))
        call wrap_put_att_text (nfid(t), cols1d_wtlnd_id(t) , 'long_name', &
             'column weight relative to corresponding landunit')

        call wrap_def_var(nfid(t), 'cols1d_itypwat', nf_int, 1, dim1_id, cols1d_itypwat_id(t))
        call wrap_put_att_text (nfid(t), cols1d_itypwat_id(t) , 'long_name', &
             'column water type')

        ! 1d pft info

        dim1_id(1) = pft_dimid

        call wrap_def_var(nfid(t), 'pfts1d_lon', ncprec, 1, dim1_id, pfts1d_lon_id(t))
        call wrap_put_att_text (nfid(t), pfts1d_lon_id(t), 'long_name', &
             '1d pft longitude')
        call wrap_put_att_text (nfid(t), pfts1d_lon_id(t), 'units', &
             'degrees_east')

        call wrap_def_var(nfid(t), 'pfts1d_lat', ncprec, 1, dim1_id, pfts1d_lat_id(t))
        call wrap_put_att_text (nfid(t), pfts1d_lat_id(t), 'long_name', &
             '1d pft latitude')
        call wrap_put_att_text (nfid(t), pfts1d_lat_id(t), 'units', &
             'degrees_north')

        call wrap_def_var(nfid(t), 'pfts1d_ixy', nf_int, 1, dim1_id, pfts1d_ixy_id(t))
        call wrap_put_att_text (nfid(t), pfts1d_ixy_id(t), 'long_name', &
             '2d longitude index of corresponding landunit')

        call wrap_def_var(nfid(t), 'pfts1d_jxy', nf_int, 1, dim1_id, pfts1d_jxy_id(t))
        call wrap_put_att_text (nfid(t), pfts1d_jxy_id(t), 'long_name', &
             '2d latitude index of corresponding landunit')

        call wrap_def_var(nfid(t), 'pfts1d_gi', nf_int, 1, dim1_id, pfts1d_gindex_id(t))
        call wrap_put_att_text (nfid(t), pfts1d_gindex_id(t) , 'long_name', &
             '1d grid index of corresponding pft')

        call wrap_def_var(nfid(t), 'pfts1d_li', nf_int, 1, dim1_id, pfts1d_lindex_id(t))
        call wrap_put_att_text (nfid(t), pfts1d_lindex_id(t) , 'long_name', &
             '1d landunit index of corresponding pft')

        call wrap_def_var(nfid(t), 'pfts1d_ci ', nf_int, 1, dim1_id, pfts1d_cindex_id(t))
        call wrap_put_att_text (nfid(t), pfts1d_cindex_id(t) , 'long_name', &
             '1d column index of corresponding pft')

        call wrap_def_var(nfid(t), 'pfts1d_wtxy', ncprec, 1, dim1_id, pfts1d_wtxy_id(t))
        call wrap_put_att_text (nfid(t), pfts1d_wtxy_id(t), 'long_name', &
             'pft weight relative to corresponding gridcell')

        call wrap_def_var(nfid(t), 'pfts1d_wtlnd', ncprec, 1, dim1_id, pfts1d_wtlnd_id(t))
        call wrap_put_att_text (nfid(t), pfts1d_wtlnd_id(t), 'long_name', &
             'pft weight relative to corresponding landunit')

        call wrap_def_var(nfid(t), 'pfts1d_wtcol', ncprec, 1, dim1_id, pfts1d_wtcol_id(t))
        call wrap_put_att_text (nfid(t), pfts1d_wtcol_id(t), 'long_name', &
             'pft weight relative to corresponding column')

        call wrap_def_var(nfid(t), 'pfts1d_itypwat', nf_int, 1, dim1_id, pfts1d_itypwat_id(t))
        call wrap_put_att_text (nfid(t), pfts1d_itypwat_id(t) , 'long_name', &
             'pft water type')

        call wrap_def_var(nfid(t), 'pfts1d_itypveg', nf_int, 1, dim1_id, pfts1d_itypveg_id(t))
        call wrap_put_att_text (nfid(t), pfts1d_itypveg_id(t) , 'long_name', &
             'pft vegetation type')

     endif

     ! Define time-constant soil variables for primary tape only
     
     if (t == 1) then
        if (tape(t)%dov2xy) then
           dim3_id(1) = lon_dimid         
           dim3_id(2) = lat_dimid         
           dim3_id(3) = levsoi_dimid
           call wrap_def_var(nfid(t), 'ZSOI'  , ncprec, 3, dim3_id, zsoi_id(t))
           call wrap_def_var(nfid(t), 'DZSOI' , ncprec, 3, dim3_id, dzsoi_id(t))
           call wrap_def_var(nfid(t), 'WATSAT', ncprec, 3, dim3_id, watsat_id(t))
           call wrap_def_var(nfid(t), 'HKSAT' , ncprec, 3, dim3_id, hksat_id(t))
           call wrap_def_var(nfid(t), 'SUCSAT', ncprec, 3, dim3_id, sucsat_id(t))
           call wrap_def_var(nfid(t), 'BSW'   , ncprec, 3, dim3_id, bsw_id(t))
        else
           dim2_id(1) = column_dimid
           dim2_id(2) = levsoi_dimid
           call wrap_def_var(nfid(t), 'ZSOI'  , ncprec, 2, dim2_id, zsoi_id(t))
           call wrap_def_var(nfid(t), 'DZSOI' , ncprec, 2, dim2_id, dzsoi_id(t))
           call wrap_def_var(nfid(t), 'WATSAT', ncprec, 2, dim2_id, watsat_id(t))
           call wrap_def_var(nfid(t), 'HKSAT' , ncprec, 2, dim2_id, hksat_id(t))
           call wrap_def_var(nfid(t), 'SUCSAT', ncprec, 2, dim2_id, sucsat_id(t))
           call wrap_def_var(nfid(t), 'BSW'   , ncprec, 2, dim2_id, bsw_id(t))
        endif
        call wrap_put_att_text (nfid(t), zsoi_id(t), 'long_name', &
             'soil depth')
        call wrap_put_att_text (nfid(t), zsoi_id(t), 'units', &
             'm')
        call wrap_put_att_realx(nfid(t), zsoi_id(t), 'missing_value',&
             ncprec, 1 ,spval)
        call wrap_put_att_realx(nfid(t), zsoi_id(t), '_FillValue', &
             ncprec, 1, spval)

        call wrap_put_att_text (nfid(t), dzsoi_id(t), 'long_name', &
             'soil thickness')
        call wrap_put_att_text (nfid(t), dzsoi_id(t), 'units', &
             'm')
        call wrap_put_att_realx(nfid(t), dzsoi_id(t), 'missing_value',&
             ncprec, 1 ,spval)
        call wrap_put_att_realx(nfid(t), dzsoi_id(t), '_FillValue', &
             ncprec, 1, spval)

        call wrap_put_att_text (nfid(t), watsat_id(t), 'long_name', &
             'saturated soil water content (porosity)')
        call wrap_put_att_text (nfid(t), watsat_id(t), 'units', &
             'mm3/mm3 ')
        call wrap_put_att_realx(nfid(t), watsat_id(t), 'missing_value',&
             ncprec, 1 ,spval)
        call wrap_put_att_realx(nfid(t), watsat_id(t), '_FillValue', &
             ncprec, 1, spval)

        call wrap_put_att_text (nfid(t), hksat_id(t), 'long_name', &
             'saturated hydraulic conductivity')
        call wrap_put_att_text (nfid(t), hksat_id(t), 'units', &
             'mm H2O /s ')
        call wrap_put_att_realx(nfid(t), hksat_id(t), 'missing_value',&
             ncprec, 1 ,spval)
        call wrap_put_att_realx(nfid(t), hksat_id(t), '_FillValue', &
             ncprec, 1, spval)

        call wrap_put_att_text (nfid(t), sucsat_id(t), 'long_name', &
             'saturated soil matric potential')
        call wrap_put_att_text (nfid(t), sucsat_id(t), 'units', &
             'mm')
        call wrap_put_att_realx(nfid(t), sucsat_id(t), 'missing_value',&
             ncprec, 1 ,spval)
        call wrap_put_att_realx(nfid(t), sucsat_id(t), '_FillValue', &
             ncprec, 1, spval)

        call wrap_put_att_text (nfid(t), bsw_id(t), 'long_name', &
             'slope of soil water retention curve')
        call wrap_put_att_text (nfid(t), bsw_id(t), 'units', &
             'unitless')
        call wrap_put_att_realx(nfid(t), bsw_id(t), 'missing_value',&
             ncprec, 1 ,spval)
        call wrap_put_att_realx(nfid(t), bsw_id(t), '_FillValue', &
             ncprec, 1, spval)
     endif

     ! --------------------------------------------------------------------
     ! Define time-dependent variables: time information
     ! --------------------------------------------------------------------
     
     ! Date and time
     
     dim1_id(1) = time_dimid
     
     call wrap_def_var (nfid(t) , 'mcdate', nf_int, 1, dim1_id  , mcdate_id(t))
     call wrap_put_att_text (nfid(t), mcdate_id(t), 'long_name', &
          'current date (YYYYMMDD)')
     
     call wrap_def_var (nfid(t) , 'mcsec' , nf_int, 1, dim1_id , mcsec_id(t))
     call wrap_put_att_text (nfid(t), mcsec_id(t), 'long_name', &
          'current seconds of current date')
     call wrap_put_att_text (nfid(t), mcsec_id(t), 'units'    , &
          's')
     
     call wrap_def_var (nfid(t) , 'mdcur' , nf_int, 1, dim1_id , mdcur_id(t))
     call wrap_put_att_text (nfid(t), mdcur_id(t), 'long_name', &
          'current day (from base day)')
     
     call wrap_def_var (nfid(t) , 'mscur' , nf_int, 1, dim1_id , mscur_id(t))
     call wrap_put_att_text (nfid(t), mscur_id(t), 'long_name', &
          'current seconds of current day')
     
     call wrap_def_var (nfid(t) , 'nstep' , nf_int, 1, dim1_id , nstep_id(t))
     call wrap_put_att_text (nfid(t), nstep_id(t), 'long_name', &
          'time step')
     
     dim2_id(1) = hist_interval_dimid
     dim2_id(2) = time_dimid
     call wrap_def_var (nfid(t), 'time_bounds', nf_double, 2, dim2_id, time_bounds_id(t))
     call wrap_put_att_text (nfid(t), time_bounds_id(t), 'long_name', &
          'history time interval endpoints')

     dim2_id(1) = strlen_dimid
     dim2_id(2) = time_dimid
     call wrap_def_var (nfid(t), 'date_written', nf_char, 2, dim2_id, date_written_id(t))
     call wrap_def_var (nfid(t), 'time_written', nf_char, 2, dim2_id, time_written_id(t))

     ! Create variables and attributes for field list

     do f = 1,tape(t)%nflds

        numlev = tape(t)%hlist(f)%field%numlev
        
        if (numlev == 1) then
           if (tape(t)%dov2xy) then
              dim3_id(1) = lon_dimid         
              dim3_id(2) = lat_dimid         
              dim3_id(3) = time_dimid         
              call wrap_def_var(nfid(t), tape(t)%hlist(f)%field%name, ncprec, 3, dim3_id, varid(f,t))
           else
              select case (tape(t)%hlist(f)%field%type1d)
              case ('gridcell') 
                 dim2_id(1) = gridcell_dimid
              case ('landunit') 
                 dim2_id(1) = landunit_dimid
              case ('column')
                 dim2_id(1) = column_dimid
              case ('pft')
                 dim2_id(1) = pft_dimid
              case default
                 write(6,*)'HTAPE_CREATE: unknown 1d type=',tape(t)%hlist(f)%field%type1d
                 call endrun ()
              end select
              dim2_id(2) = time_dimid                                              
              call wrap_def_var(nfid(t), tape(t)%hlist(f)%field%name, ncprec, 2, dim2_id, varid(f,t))
           endif
        else 
           if (tape(t)%dov2xy) then
              dim4_id(1) = lon_dimid      
              dim4_id(2) = lat_dimid   
              dim4_id(3) = levsoi_dimid   
              dim4_id(4) = time_dimid      
              call wrap_def_var(nfid(t), tape(t)%hlist(f)%field%name, ncprec, 4, dim4_id, varid(f,t))
           else
              select case (tape(t)%hlist(f)%field%type1d)
              case ('gridcell') 
                 dim3_id(1) = gridcell_dimid 
              case ('landunit') 
                 dim3_id(1) = landunit_dimid
              case ('column')
                 dim3_id(1) = column_dimid
              case ('pft')
                 dim3_id(1) = pft_dimid
              case default
                 write(6,*)'HTAPE_CREATE: unknown 1d type=',tape(t)%hlist(f)%field%type1d
                 call endrun ()
              end select
              dim3_id(2) = levsoi_dimid                                             
              dim3_id(3) = time_dimid                                             
              call wrap_def_var(nfid(t), tape(t)%hlist(f)%field%name, ncprec, 3, dim3_id, varid(f,t))
           endif
        endif

        str = tape(t)%hlist(f)%field%units
        if ( str(1:1) /= ' ' ) then
           call wrap_put_att_text (nfid(t), varid(f,t), 'units', str)
        end if

        str = tape(t)%hlist(f)%field%long_name
        call wrap_put_att_text (nfid(t), varid(f,t), 'long_name', str)

        avgflag = tape(t)%hlist(f)%avgflag
        select case (avgflag)
        case ('A')
           str = 'mean'
        case ('I')
           str = 'instantaneous'
        case ('X')
           str = 'maximum'
        case ('M')
           str = 'minimum'
        case default
           write(6,*)' HTAPE_CREATE: unknown time averaging flag (avgflag)=',avgflag
           call endrun ()
        end select
        call wrap_put_att_text (nfid(t), varid(f,t), 'cell_method', 'time: '//str)

        call wrap_put_att_realx(nfid(t), varid(f,t), '_FillValue', ncprec, 1, spval)

        call wrap_put_att_realx(nfid(t), varid(f,t), 'missing_value', ncprec,1 ,spval)

     end do

     ! End of define mode 

     ret = nf_enddef (nfid(t))
     write(6,*)'HTAPE_CREATE: Successfully defined netcdf history file ',t

     ! Write out time-invariant grid variables

#if (defined OFFLINE)
     if (.not. offline_rdgrid) then
        call wrap_put_var_realx (nfid(t), edgen_id(t), lsmedge(1))
        call wrap_put_var_realx (nfid(t), edgee_id(t), lsmedge(2))
        call wrap_put_var_realx (nfid(t), edges_id(t), lsmedge(3))
        call wrap_put_var_realx (nfid(t), edgew_id(t), lsmedge(4))
     endif
#endif
     if (fullgrid) then
        lonvar(1:lsmlon) = longxy(1:lsmlon,1)
        call wrap_put_var_realx (nfid(t), lonvar_id(t), lonvar)
        latvar(1:lsmlat) = latixy(1,1:lsmlat)
        call wrap_put_var_realx (nfid(t), latvar_id(t), latvar)
     endif
     call wrap_put_var_realx (nfid(t), levvar_id(t), zsoi)
     call wrap_put_var_realx (nfid(t), longxy_id(t), longxy)
     call wrap_put_var_realx (nfid(t), latixy_id(t), latixy)
     call wrap_put_var_realx (nfid(t), area_id(t), area) 
     call wrap_put_var_realx (nfid(t), landfrac_id(t), landfrac) 
     call wrap_put_var_int   (nfid(t), landmask_id(t), landmask)
     call wrap_put_var_int   (nfid(t), numlon_id(t), numlon)
     
     ! Output time-independent 1d-indices surrounding datatypes FILL

     if (.not. tape(t)%dov2xy) then
        call wrap_put_var_realx(nfid(t), grid1d_lon_id(t), grid1d%londeg(:))
        call wrap_put_var_realx(nfid(t), grid1d_lat_id(t), grid1d%latdeg(:))
        call wrap_put_var_int  (nfid(t), grid1d_ixy_id(t), grid1d%ixy(:))
        call wrap_put_var_int  (nfid(t), grid1d_jxy_id(t), grid1d%jxy(:))

        call wrap_put_var_realx(nfid(t), land1d_lon_id(t), land1d%londeg(:))
        call wrap_put_var_realx(nfid(t), land1d_lat_id(t), land1d%latdeg(:))
        call wrap_put_var_int  (nfid(t), land1d_ixy_id(t), land1d%ixy(:))
        call wrap_put_var_int  (nfid(t), land1d_jxy_id(t), land1d%jxy(:))
        call wrap_put_var_int  (nfid(t), land1d_gindex_id(t), land1d%gindex(:))
        call wrap_put_var_realx(nfid(t), land1d_wtxy_id(t), land1d%wtxy(:))
        call wrap_put_var_int  (nfid(t), land1d_itypwat_id(t), land1d%itypwat(:))

        call wrap_put_var_realx(nfid(t), cols1d_lon_id(t), cols1d%londeg(:))
        call wrap_put_var_realx(nfid(t), cols1d_lat_id(t), cols1d%latdeg(:))
        call wrap_put_var_int  (nfid(t), cols1d_ixy_id(t), cols1d%ixy(:))
        call wrap_put_var_int  (nfid(t), cols1d_jxy_id(t), cols1d%jxy(:))
        call wrap_put_var_int  (nfid(t), cols1d_gindex_id(t), cols1d%gindex(:))
        call wrap_put_var_int  (nfid(t), cols1d_lindex_id(t), cols1d%lindex(:))
        call wrap_put_var_realx(nfid(t), cols1d_wtxy_id(t), cols1d%wtxy(:))
        call wrap_put_var_realx(nfid(t), cols1d_wtlnd_id(t), cols1d%wtlnd(:))
        call wrap_put_var_int  (nfid(t), cols1d_itypwat_id(t), cols1d%itypwat(:))

        call wrap_put_var_realx(nfid(t), pfts1d_lon_id(t), pfts1d%londeg(:))
        call wrap_put_var_realx(nfid(t), pfts1d_lat_id(t), pfts1d%latdeg(:))
        call wrap_put_var_int  (nfid(t), pfts1d_ixy_id(t), pfts1d%ixy(:))
        call wrap_put_var_int  (nfid(t), pfts1d_jxy_id(t), pfts1d%jxy(:))
        call wrap_put_var_int  (nfid(t), pfts1d_gindex_id(t), pfts1d%gindex(:))
        call wrap_put_var_int  (nfid(t), pfts1d_lindex_id(t), pfts1d%lindex(:))
        call wrap_put_var_int  (nfid(t), pfts1d_cindex_id(t), pfts1d%cindex(:))
        call wrap_put_var_realx(nfid(t), pfts1d_wtxy_id(t), pfts1d%wtxy(:))
        call wrap_put_var_realx(nfid(t), pfts1d_wtlnd_id(t), pfts1d%wtlnd(:))
        call wrap_put_var_realx(nfid(t), pfts1d_wtcol_id(t), pfts1d%wtcol(:))
        call wrap_put_var_int  (nfid(t), pfts1d_itypwat_id(t), pfts1d%itypwat(:))
        call wrap_put_var_int  (nfid(t), pfts1d_itypveg_id(t), pfts1d%itypveg(:))
     endif

  end subroutine htape_create

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: htape_timeconst
!
! !INTERFACE:
  subroutine htape_timeconst (t)
!
! !DESCRIPTION: 
! Write time constant values to primary history tape.
! Issue the required netcdf wrapper calls to define the history file contents.
!
! !USES
    use clmpoint
    use clm_varpar, only : lsmlon, lsmlat, nlevsoi
    use clmtype, only : cols1d
    use mapxy, only : vec2xy
#if ( defined SPMD )
    use spmdMod, only : masterproc, gather_data_to_master
#else
    use spmdMod, only : masterproc
#endif
    use system_messages
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: t   ! tape index
!
! !REVISION HISTORY
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES
    integer :: ci,k,f,lev,ifld           ! indices
    integer :: dim1_id(1)                ! netCDF dimension id for 1-d variables
    integer :: dim2_id(2)                ! netCDF dimension id for 2-d variables
    integer :: dim3_id(3)                ! netCDF dimension id for 3-d variables
    integer :: dim4_id(4)                ! netCDF dimension id for 4-d variables
    integer :: vid(6),idx(6)             ! temporaries for time-invariant output
    integer :: start(4)                  ! starting indices for netcdf field
    integer :: count(4)                  ! count values for netcdf field
    integer :: num1d,beg1d,end1d         ! 1d size, beginning and ending indices
    integer :: ier                       ! error status 
    real(r8), pointer :: rloc1d(:)       ! temporary
    real(r8), pointer :: rglob1d(:)      ! temporary
    real(r8), pointer :: rglob2d(:,:)    ! temporary
    type(column_type), pointer :: c      !local pointer to derived subtype
    real(r8) :: mlfxy(lsmlon,lsmlat,nlevsoi)! grid-average multi-level soil field
!-----------------------------------------------------------------------

    ! Output time constant history variables
    
    beg1d = cols1d%beg
    end1d = cols1d%end
    num1d = cols1d%num
    
    if (tape(t)%dov2xy) then
       start(1) = 1
       count(1) = lsmlon
       start(2) = 1
       count(2) = lsmlat
       start(3) = 1
       count(3) = nlevsoi
    else
       start(1) = 1
       count(1) = num1d
       start(2) = 1
       count(2) = nlevsoi
    endif
    
    
    allocate (rloc1d(beg1d:end1d), stat=ier)
    call allocation_err(ier, 'htape_timeconst', 'rloc1d', end1d-beg1d+1)
    if (masterproc) then
       allocate(rglob1d(num1d)) 
       call allocation_err(ier, 'htape_timeconst', 'rglob1d', num1d)
       allocate(rglob2d(num1d,nlevsoi))
       call allocation_err(ier, 'htape_timeconst', 'rglob1d', num1d*nlevsoi)
    endif

    vid(1) = zsoi_id(t)
    vid(2) = dzsoi_id(t)
    vid(3) = watsat_id(t)
    vid(4) = sucsat_id(t)
    vid(5) = bsw_id(t)
    vid(6) = hksat_id(t)

    do ifld = 1,6
       do lev = 1,nlevsoi
          do ci = beg1d,end1d
             c => cpoint%col(ci)%c
             if (c%cps%lps%lakpoi) then
                rloc1d(ci) = spval
             else
                if (ifld ==1) rloc1d(ci) = c%cps%z(lev)
                if (ifld ==2) rloc1d(ci) = c%cps%dz(lev)
                if (ifld ==3) rloc1d(ci) = c%cps%lps%watsat(lev)
                if (ifld ==4) rloc1d(ci) = c%cps%lps%sucsat(lev)
                if (ifld ==5) rloc1d(ci) = c%cps%lps%bsw(lev)
                if (ifld ==6) rloc1d(ci) = c%cps%lps%hksat(lev)
             endif
          end do
#if (defined SPMD)
          call gather_data_to_master(rloc1d, rglob1d, clmlevel=cols1d%name)
#else
          rglob1d(:) = rloc1d(:)
#endif
          if (masterproc) then
             rglob2d(:,lev) = rglob1d(:)
          endif
       end do
       if (masterproc) then
          if (tape(t)%dov2xy) then
             do lev = 1,nlevsoi
                rglob1d(:) = rglob2d(:,lev)
                call vec2xy (cols1d%name, rglob1d, spval, mlfxy(1,1,lev))
             end do
             call wrap_put_vara_realx (nfid(t), vid(ifld), start, count, mlfxy)
          else
             call wrap_put_vara_realx (nfid(t), vid(ifld), start, count, rglob2d)
          endif
       endif
    end do

    deallocate (rloc1d)
    if (masterproc) then
       deallocate(rglob2d)
       deallocate(rglob1d)
    endif

    return
  end subroutine htape_timeconst 

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: hfields_write
!
! !INTERFACE:
  subroutine hfields_write (t)
!
! !DESCRIPTION: 
! Write history tape. If SPMD, first gather the data to the master processor.  
! Issue the netcdf call to write the variable.
!
! !USES
    use clmpoint
    use clm_varpar, only : lsmlon, lsmlat, nlevsoi
#if (defined COUP_TEM)
!    use clmtype, only : grid1d,cols1d,pfts1d,pps
    use clmtype
    use time_manager, only : get_curr_date, get_nstep
#endif
    use mapxy, only : vec2xy
#if ( defined SPMD )
    use spmdMod, only : masterproc, gather_data_to_master
#else
    use spmdMod, only : masterproc
#endif
    use shr_sys_mod, only : shr_sys_flush
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: t                ! tape index
!
! !REVISION HISTORY
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES
    integer :: f                            ! field index 
    integer :: k                            ! 1d index
    integer :: lev                          ! level index
    integer :: size                         ! allocation size    
    integer :: start(4)                     ! array of starting field indices for netcdf field
    integer :: count(4)                     ! array of count values for netcdf field
    integer :: numlev                       ! number of vertical levels 
    integer :: num1d                        ! 1d size 
    integer :: beg1d,end1d                  ! 1d beginning and ending indices
    character(len=8)  :: type1d             ! 1d type
    real(r8), pointer :: rloc1d(:)          ! temporary
    real(r8), pointer :: rglob1d(:)         ! temporary
    real(r8), pointer :: rglob2d(:,:)       ! temporary
    real(r8), pointer :: temploc2d(:,:)     ! temporary
    real(r8), pointer :: tempglob2d(:,:)    ! temporary
    real(r8) :: slfxy(lsmlon,lsmlat)        ! grid-average single-level field
    real(r8) :: mlfxy(lsmlon,lsmlat,nlevsoi)! grid-average multi-level soil field
    integer  :: hpindex                     ! temporary - testing
    integer  :: iii                         !jrs hack
#if (defined COUP_TEM) 
    real(r8), pointer :: rqvege1d(:)        ! temporary
    real(r8), pointer :: rqvegt1d(:)        ! temporary
    real(r8), pointer :: rqover1d(:)        ! temporary
    real(r8), pointer :: rqdrai1d(:)        ! temporary
    real(r8), pointer :: rqsoil1d(:)        ! temporary
    real(r8), pointer :: rqstrm1d(:)        ! temporary
    real(r8), pointer :: rstrmdur1d(:)      ! temporary
    real(r8), pointer :: rstorms1d(:)       ! temporary
    real(r8), pointer :: rstrmdry1d(:)      ! temporary
    real(r8), pointer :: rsice2d(:,:)       ! temporary
    real(r8), pointer :: rsliq2d(:,:)       ! temporary
    real(r8), pointer :: rsnice2d(:)       ! temporary
    real(r8), pointer :: rsnliq2d(:)       ! temporary
    real(r8), pointer :: rstemp2d(:,:)      ! temporary
    integer :: i,ii,j,jj,ip,pp              !indices
    integer :: mcsec,hr,day,mon,yr          !time indices
    integer :: nstep
    integer, parameter :: mxmsaics = 35
    integer, parameter :: mxmthdys = 31
    integer, parameter :: mxdayhrs = 24
    integer, parameter :: max3hrs = 8
    integer, parameter :: mxnlayers = 6
    integer, parameter :: totlayers = 10
    integer, save :: xstrmlast(lsmlat)
!    integer, save :: xstrmadj(lsmlat)
    integer :: exstrm
    real(r8) :: mitco2(lsmlat)
    real(r8) :: mito3(max3hrs,lsmlat)
    real(r8) :: mittemp(lsmlat)
    real(r8) :: mitdaytemp(mxmthdys,lsmlat)
    real(r8) :: mitswrs(lsmlat)
    real(r8) :: mitpre(lsmlat)
    real(r8) :: mitstrmdur(mxmthdys,lsmlat)
    real(r8) :: mitqstrm(mxmthdys,lsmlat),laststrm
    real(r8) :: mitaet(mxmsaics,lsmlat)
    real(r8) :: mitsh2o1m(mxmsaics,lsmlat)
    real(r8) :: mitsh2o2m(mxmsaics,lsmlat)
    real(r8) :: mitswe(mxmsaics,lsmlat)
    real(r8) :: mitsfr(mxmsaics,lsmlat)
    real(r8) :: mitdrn(mxmsaics,lsmlat)
    real(r8) :: mitdaytsoil(mxmthdys,mxmsaics,lsmlat,totlayers)
    real(r8) :: mitdaysh2o(mxmthdys,mxmsaics,lsmlat,totlayers)
    real(r8) :: mithrsh2o(mxdayhrs,mxmthdys,mxmsaics,lsmlat,mxnlayers)
    integer :: mitstorm(mxmthdys,lsmlat)
    common/climate4tem/mitco2,mito3,mittemp,mitdaytemp,mitswrs,mitpre,mitstrmdur,&
                       mitqstrm,mitaet,mitsh2o1m,mitsh2o2m,mitswe,mitsfr,mitdrn,&
                       mitdaytsoil,mitdaysh2o,&
                       mithrsh2o
    logical :: firstday
    logical :: firsthr
    logical :: dodry
    data firstday/.true./
    data firsthr/.true./
    character(len=14) :: fldlst(40)
    data fldlst( 1) /'pfts1d_lat    '/
    data fldlst( 2) /'cols1d_typwat '/
    data fldlst( 3) /'pfts1d_typveg '/
    data fldlst( 4) /'grid1d_lat    '/
    data fldlst(10) /'TBOT          '/
    data fldlst(11) /'QVEGE         '/
    data fldlst(12) /'QVEGT         '/
    data fldlst(13) /'QSOIL         '/
    data fldlst(14) /'SOILICE       '/
    data fldlst(15) /'SOILLIQ       '/
    data fldlst(16) /'SNOWLIQ       '/
    data fldlst(17) /'SNOWICE       '/
    data fldlst(18) /'QOVER         '/
    data fldlst(19) /'QDRAI         '/
    data fldlst(21) /'QSTRM         '/
    data fldlst(22) /'STRMDUR       '/
    data fldlst(23) /'STORMS        '/
    data fldlst(24) /'STRMDRY       '/
    data fldlst(25) /'TSOI          '/
    data fldlst(26) /'H2OSOI        '/
    data fldlst(27) /'TSA           '/
    data fldlst(30) /'H2OSOI        '/
    real :: lats(lsmlat)                !atm latitude points
    data lats /-90,-86,-82,-78,-74,-70,-66,-62,-58,-54,-50,-46,-42,-38,-34,-30,-26,-22,-18,-14,-10,-6,-2, &
               2,6,10,14,18,22,26,30,34,38,42,46,50,54,58,62,66,70,74,78,82,86,90/
    integer :: mdays(12),ndays
    integer :: prvmon,last
    data mdays /31,28,31,30,31,30,31,31,30,31,30,31/
#endif
#if ( defined COUP_URBAN)
    real :: tscreen(lsmlon,lsmlat,mxmsaics)      ! Screen temperature (for each PFT) to Urban airshed model
    real :: dlypcp(lsmlon,lsmlat,mxmsaics)       ! Daily precipitation (for each PFT) to Urbran airshed model
    real(r8), pointer :: rst2m2d(:)              ! temporary pointer for screen temperature
    real(r8), pointer :: pcp2pft                 !PCP Distribution factor across PFTs
    type(pft_type)             , pointer :: p    ! local pointers to derived subtypes
    type(column_type)          , pointer :: c    ! local pointers to derived subtypes
    type(pft_pstate_type)      , pointer :: pps  ! local pointers to derived subtypes
    common/chem_meta_gls/tscreen,dlypcp
#endif
!-----------------------------------------------------------------------

#if (defined COUP_TEM)
   if (t == 2) then
!
!	Monthly PFT Data for TEM
!
    mitsh2o1m=0.0
    mitsh2o2m=0.0
    mitswe=0.0
    do f=1,tape(t)%nflds
!      write(6,*)'HFIELD_SAVE: saving ',tape(t)%hlist(f)%field%name,' for monthly TEM data'
      if (tape(t)%hlist(f)%field%name == fldlst(11)) rqvege1d => tape(t)%hlist(f)%hbuf(:,1)
      if (tape(t)%hlist(f)%field%name == fldlst(12)) rqvegt1d => tape(t)%hlist(f)%hbuf(:,1)
      if (tape(t)%hlist(f)%field%name == fldlst(13)) rqsoil1d => tape(t)%hlist(f)%hbuf(:,1)
      if (tape(t)%hlist(f)%field%name == fldlst(14)) rsliq2d => tape(t)%hlist(f)%hbuf
      if (tape(t)%hlist(f)%field%name == fldlst(15)) rsice2d => tape(t)%hlist(f)%hbuf
      if (tape(t)%hlist(f)%field%name == fldlst(16)) rsnliq2d => tape(t)%hlist(f)%hbuf(:,1)
      if (tape(t)%hlist(f)%field%name == fldlst(17)) rsnice2d => tape(t)%hlist(f)%hbuf(:,1)
      if (tape(t)%hlist(f)%field%name == fldlst(18)) rqover1d => tape(t)%hlist(f)%hbuf(:,1)
      if (tape(t)%hlist(f)%field%name == fldlst(19)) rqdrai1d => tape(t)%hlist(f)%hbuf(:,1)
    enddo
    num1d = tape(1)%hlist(1)%field%num1d
    do j=1,lsmlat
      do jj=1,num1d
          !write(6,*) 'Latitude of PFT ',pfts1d%latdeg(jj)
       if (pfts1d%latdeg(jj)==lats(j)) then
         if(cols1d%itypwat(jj)==2) then
           ip=31
         elseif(cols1d%itypwat(jj)==5) then
           ip=18
         elseif(cols1d%itypwat(jj)==3) then
           ip=32
         elseif(cols1d%itypwat(jj)==1) then
           do pp=0,16
             if (pfts1d%itypveg(jj)==pp) ip=pp+1
           enddo
         else
           write(6,*) 'Land cover type not found for latitude: ',lats(j)
           write(6,*) 'Land cover type = ',cols1d%itypwat(jj)
         endif
!
!       Fill in data for latitude and land cover type/PFT
!
         mitaet(ip,j)=(rqvege1d(jj)+rqvegt1d(jj)+rqsoil1d(jj))*86400
         mitsfr(ip,j)=(rqover1d(jj))*86400
         mitdrn(ip,j)=(rqdrai1d(jj))*86400
         mitswe(ip,j)=rsnliq2d(jj)+rsnice2d(jj)
         if ( mitswe(ip,j) > 1e+35 ) mitswe(ip,j) = 0.0
         do lev=1,8
           mitsh2o1m(ip,j)=mitsh2o1m(ip,j)+rsliq2d(jj,lev)+rsice2d(jj,lev)
         enddo
         mitsh2o2m(ip,j)=mitsh2o1m(ip,j)+rsliq2d(jj,9)+rsice2d(jj,9)
       endif
     enddo !jj
!     do jj=1,grid1d%num
!       if (x1d(4,jj)==lats(j)) then
!          mittemp(j)=rtbot1d(jj)
!       endif
!     enddo
    enddo !j
!10  format(19f10.3)
!    write(6,*) 'AET'
!    do j=1,lsmlat
!      write(6,10)(mitaet(ip,j),ip=1,19)
!    enddo
!    write(6,*) 'SM1M'
!    do j=1,lsmlat
!      write(6,10) (mitsh2o1m(ip,j),ip=1,19)
!    enddo
!    write(6,*) 'SM2M'
!    do j=1,lsmlat
!      write(6,10)(mitsh2o2m(ip,j),ip=1,19)
!    enddo
!
!	Daily PFT Data NEM (CH4 and N2O)
!
   elseif (t == 3) then
!
!	Daily PFT Data NEM (CH4 and N2O)
!
    num1d = tape(1)%hlist(1)%field%num1d
    do f=1,tape(t)%nflds
!      write(6,*)'HFIELD_SAVE: saving ',tape(t)%hlist(f)%field%name,' for Daily TEM data'
      if (tape(t)%hlist(f)%field%name == fldlst(21)) rqstrm1d => tape(t)%hlist(f)%hbuf(:,1)
      if (tape(t)%hlist(f)%field%name == fldlst(22)) rstrmdur1d => tape(t)%hlist(f)%hbuf(:,1)
      if (tape(t)%hlist(f)%field%name == fldlst(23)) rstorms1d => tape(t)%hlist(f)%hbuf(:,1)
      if (tape(t)%hlist(f)%field%name == fldlst(25)) rstemp2d => tape(t)%hlist(f)%hbuf
      if (tape(t)%hlist(f)%field%name == fldlst(26)) rsliq2d => tape(t)%hlist(f)%hbuf
#if ( defined COUP_URBAN)
      if (tape(t)%hlist(f)%field%name == fldlst(27)) rst2m2d => tape(t)%hlist(f)%hbuf(:,1)
#endif
    enddo

    ! Set calendar for current time step

    call get_curr_date (yr, mon, day, mcsec) 
!    write (6,*) 'mon= ',mon,' day= ',day,'sec= ',mcsec
!
!	Offset day for the daily processing
!
    prvmon=mon-1
    if (prvmon == 0) prvmon=12
     day = day - 1
    if (day < 1) then
      day = mdays(prvmon)
    endif
!    write(6,*)'HFIELD_SAVE: saving at day ',day
    if (firstday) then     
       xstrmlast=0
       mitqstrm=0.
       mitstrmdur=0.
       mitstorm=0
       firstday=.false.
     endif
!
!	At beginning of each month, adjust counters so that
!	number of storms for NEM is equal to zero
!

    do j=1,lsmlat
      do jj=1,num1d
       if (pfts1d%latdeg(jj)==lats(j)) then
         if(cols1d%itypwat(jj)==2) then
           ip=31
         elseif(cols1d%itypwat(jj)==5) then
           ip=18
         elseif(cols1d%itypwat(jj)==3) then
           ip=32
         elseif(cols1d%itypwat(jj)==1) then
           do pp=0,16
             if (pfts1d%itypveg(jj)==pp) ip=pp+1
           enddo
         else
           write(6,*) 'Land cover type not found for latitude: ',lats(j)
           write(6,*) 'Land cover type = ',cols1d%itypwat(jj)
         endif
!
!       Fill in data for latitude and land cover type/PFT
!
#if ( defined COUP_URBAN)
         tscreen(1,j,ip) = rst2m2d(jj)
         p => c%p(jj)
         pps => p%pps
         pcp2pft => pps%pcp2pft
         dlypcp(1,j,ip) = (rqstrm1d(jj)*24*3600)*pcp2pft
#endif
         do lev=1,totlayers
           mitdaytsoil(day,ip,j,lev)=rstemp2d(jj,lev)
           mitdaysh2o(day,ip,j,lev)=rsliq2d(jj,lev)
         enddo
       endif
      enddo !jj
    do jj=1,grid1d%num
       if (grid1d%latdeg(jj)==lats(j)) then
!
! 	lump multiple storms that occur in a day
!	for NEM, which cannot distinguish them.
!
        exstrm = rstorms1d(jj)-xstrmlast(jj)-1
        if ( exstrm > 0 ) then
!          write(6,*) lats(jj),rstorms1d(jj)-xstrmadj(jj),xstrmlast(jj),exstrm
!
!	increase daily average storm duration to reflect lumped multiple events
!	Total lumped storm duration cannot be more than a day
!
          rstrmdur1d(jj) = min((exstrm+1)*rstrmdur1d(jj),3600*24.)
        endif
        xstrmlast(jj) = rstorms1d(jj)
!
!       Convert precipitation statistics to default units for NEM
!
         mitqstrm(day,j)=(rqstrm1d(jj)*3600)/10.
         mitstrmdur(day,j)=rstrmdur1d(jj)/3600
       endif
    enddo !jj
   enddo

!
!	Modifying daily precip statistics for N2O model
!
    where (mitqstrm > 1.0e+20) mitqstrm=0.
    where (mitstrmdur > 1.0e+20) mitstrmdur=0.
    where (mitstorm > 1.0e+20) mitstorm=0.
!
!	Minimum value of storm duration is one-hour (timestep of model)
!
    where ( mitstrmdur < 1. .and. mitqstrm > 0.) mitstrmdur = 1.
!
!    write(6,*) 'day',day,'storms',mitstorm(day,23)
!    write(6,*) 'day',day,'qstorm',mitqstrm(day,23)
!    write(6,*) 'day',day,'strmdur',mitstrmdur(day,23)
!    do j=1,lsmlat
!      j = 23
!      write(6,*) 'LAT=',lats(j)
!      write(6,*) 'day',day,'storms',(mitstorm(i,j),i=1,31)
!      write(6,*) 'day',day,'qstorm',(mitqstrm(i,j),i=1,31)
!      write(6,*) 'day',day,'strmdur',(mitstrmdur(i,j),i=1,31)
!      j = 33
!      write(6,*) 'LAT=',lats(j)
!      write(6,*) 'day',day,'storms',(mitstorm(i,j),i=1,31)
!      write(6,*) 'day',day,'qstorm',(mitqstrm(i,j),i=1,31)
!      write(6,*) 'day',day,'strmdur',(mitstrmdur(i,j),i=1,31)
!    enddo
!
!	Processing Daily Storm Statistics at end of month
!	for NEM Compatibility
!
   if (day > 1) then 
!	
!     Keep storm that lasts < day confined to a calendar day
!     (no midnight carry-over of storm seen in NEM)
!
     do j=1,lsmlat
       laststrm=mitqstrm(day-1,j)
       if (mitqstrm(day,j) == laststrm) then
          if (mitstrmdur(day,j) > 0. .and. mitstrmdur(day,j) <= 24.) then
            mitstrmdur(day,j) = 0.
            mitqstrm(day,j) = 0.
          endif
       endif
    enddo !j
!      j = 23
!      write(6,*) 'REPROCESSED STORM STATISTICS'
!      write(6,*) '============================'
!      write(6,*) 'LAT=',lats(j)
!      write(6,*) 'day',day,'strmdry',(mitstrmdry(i,j),i=1,31)
!      write(6,*) 'day',day,'storms',(mitstorm(i,j),i=1,31)
!      write(6,*) 'day',day,'qstorm',(mitqstrm(i,j),i=1,31)
!      write(6,*) 'day',day,'strmdur',(mitstrmdur(i,j),i=1,31)
!      j = 33
!      write(6,*) 'LAT=',lats(j)
!      write(6,*) 'day',day,'strmdry',(mitstrmdry(i,j),i=1,31)
!      write(6,*) 'day',day,'storms',(mitstorm(i,j),i=1,31)
!      write(6,*) 'day',day,'qstorm',(mitqstrm(i,j),i=1,31)
!      write(6,*) 'day',day,'strmdur',(mitstrmdur(i,j),i=1,31)
   endif
   elseif (t == 4) then
!
!	Hourly PFT Data for TEM
!
    num1d = tape(1)%hlist(1)%field%num1d
    do f=1,tape(t)%nflds
      !write(6,*)'HFIELD_SAVE: saving ',tape(t)%hlist(f)%field%name,' for hourly TEM data'
      if (tape(t)%hlist(f)%field%name == fldlst(30)) rsliq2d => tape(t)%hlist(f)%hbuf
    enddo

    ! Set calendar for current time step

    call get_curr_date (yr, mon, day, mcsec) 
    hr = (mcsec/3600)+1
    !write(6,*)'HFIELD_SAVE: saving at hour',hr,' day ',day

    do j=1,lsmlat
      do jj=1,num1d
       if (pfts1d%latdeg(jj)==lats(j)) then
         if(cols1d%itypwat(jj)==2) then
           ip=31
         elseif(cols1d%itypwat(jj)==5) then
           ip=18
         elseif(cols1d%itypwat(jj)==3) then
           ip=32
         elseif(cols1d%itypwat(jj)==1) then
           do pp=0,16
             if (pfts1d%itypveg(jj)==pp) ip=pp+1
           enddo
         else
           write(6,*) 'Land cover type not found for latitude: ',lats(j)
           write(6,*) 'Land cover type = ',cols1d%itypwat(jj)
         endif
!
!       Fill in data for latitude and land cover type/PFT
!
         do lev=1,mxnlayers
           mithrsh2o(hr,day,ip,j,lev)=rsliq2d(jj,lev)
         enddo
       endif
     enddo !jj
!    if (j==33) write(6,*)mithrsh2o(hr,day,:,j,3)
    enddo !j
    if (firsthr) then
      write (6,*) 'Assigning first hour of soil water profile'
      mithrsh2o(hr-1,day,:,:,:)=mithrsh2o(hr,day,:,:,:)
      firsthr=.false.
    endif
   else
#endif
    ! Start loop over fields

    do f=1,tape(t)%nflds

#if (defined DEBUG)
       write(6,*)'HFIELD_WRITE: writing ',tape(t)%hlist(f)%field%name
#endif

       numlev = tape(t)%hlist(f)%field%numlev
       type1d = tape(t)%hlist(f)%field%type1d
       beg1d = tape(t)%hlist(f)%field%beg1d
       end1d = tape(t)%hlist(f)%field%end1d
       num1d = tape(t)%hlist(f)%field%num1d
       
       if (numlev == 1) then
          if (tape(t)%dov2xy) then
             start(1) = 1
             count(1) = lsmlon
             start(2) = 1
             count(2) = lsmlat
             start(3) = tape(t)%ntimes
             count(3) = 1
          else
             start(1) = 1
             count(1) = num1d
             start(2) = tape(t)%ntimes
             count(2) = 1
          endif
       else if (numlev > 1) then
          if (tape(t)%dov2xy) then
             start(1) = 1
             count(1) = lsmlon
             start(2) = 1
             count(2) = lsmlat
             start(3) = 1
             count(3) = numlev
             start(4) = tape(t)%ntimes
             count(4) = 1
          else
             start(1) = 1
             count(1) = num1d
             start(2) = 1
             count(2) = numlev
             start(3) = tape(t)%ntimes
             count(3) = 1
          endif
       else
          write(6,*)'HFIELD_WRITE: bad numlev=', numlev
          call endrun
       end if
        
#if (defined DEBUG)
       write(6,*)'HFIELD_WRITE: writing time indx ', tape(t)%ntimes, &
            ' field id', varid(f,t), &
            ' field name ', tape(t)%hlist(f)%field%name
       write(6,*)'start=',start
       write(6,*)'count=',count
       call print_memusage ()
#endif
     
       ! Write single-level fields
        
       if (numlev == 1) then
#if (defined SPMD)
          allocate (rloc1d(beg1d:end1d))
          do k = beg1d,end1d
             rloc1d(k) = tape(t)%hlist(f)%hbuf(k,1)
          end do
          if (masterproc) allocate(rglob1d(num1d)) 
          call gather_data_to_master(rloc1d, rglob1d, clmlevel=type1d)
#else
          rglob1d => tape(t)%hlist(f)%hbuf(:,1)
#endif
          if (masterproc) then
             if (tape(t)%dov2xy) then
                call vec2xy (type1d, rglob1d, spval, slfxy)
!start jrs hack to fix netcdf write garbage crash
                if (varid(f,t) == 79) then
                  do iii = 1,9
                     slfxy(1,iii)=0.
                  end do
                  do iii = 44,46
                     slfxy(1,iii)=0.
                  end do
                endif
!end jrs hack
                call wrap_put_vara_realx (nfid(t), varid(f,t), start, count, slfxy)
             else
                call wrap_put_vara_realx (nfid(t), varid(f,t), start, count, rglob1d)
             endif
          endif
#if (defined SPMD)
          deallocate(rloc1d)
          if (masterproc) deallocate(rglob1d)
#endif
       endif
        
       ! Write multi-level fields
       
       if (numlev > 1) then
#if (defined SPMD)
          allocate(temploc2d(numlev,beg1d:end1d))
          if (masterproc) then
             allocate(tempglob2d(numlev,num1d))
             allocate(rglob2d(num1d,numlev))
          endif
          do lev = 1,numlev
             do k = beg1d, end1d
                temploc2d(lev,k) = tape(t)%hlist(f)%hbuf(k,lev)
             end do
          end do
          call gather_data_to_master(temploc2d, tempglob2d, clmlevel=type1d)
          if (masterproc) then
             do lev = 1,numlev
                do k = 1,num1d
                   rglob2d(k,lev) = tempglob2d(lev,k)
                end do
             end do
          endif
#else
          rglob2d => tape(t)%hlist(f)%hbuf
#endif
          if (masterproc) then
             if (tape(t)%dov2xy) then
                do lev = 1,numlev
                   rglob1d => rglob2d(:,lev)
                   call vec2xy (type1d, rglob1d, spval, mlfxy(1,1,lev))
                end do
                call wrap_put_vara_realx (nfid(t), varid(f,t), start, count, mlfxy)
             else
                call wrap_put_vara_realx (nfid(t), varid(f,t), start, count, rglob2d)
             endif
          endif
#if (defined SPMD)
          deallocate(temploc2d)
          if (masterproc) then
             deallocate(tempglob2d)
             deallocate(rglob2d)
          endif
#endif
        end if

     ! End of loop over fields

     end do   
        
#if (defined COUP_TEM)
     endif
#endif
     return
  end subroutine hfields_write

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: htapes_wrapup
!
! !INTERFACE:
  subroutine htapes_wrapup ()
!
! !DESCRIPTION: 
! Write and/or dispose history tape(s)
! Determine if next time step is beginning of history interval and if so:
!   increment the current time sample counter, open a new history file 
!   if needed (i.e., when ntim = 1), write history data to current history 
!   file, reset field accumulation counters to zero
! If primary history file is full or at the last time step of the simulation,
!   write restart dataset and close and dispose all history fiels 
! If history file is full or at the last time step of the simulation:
!   close history file and dispose to mass store (only if file is open)
!   reset time sample counter to zero if file is full
! Daily-averaged data for the first day in September are written on 
!   date = 00/09/02 with mscur = 0
! Daily-averaged data for the first day in month mm are written on 
!   date = yyyy/mm/02 with mscur = 0
! Daily-averaged data for the 30th day (last day in September) are written on 
!   date = 0000/10/01 mscur = 0
! Daily-averaged data for the last day in month mm are written on 
!   date = yyyy/mm+1/01 with mscur = 0
!
! !USES
    use clm_varctl, only : archive_dir, mss_wpass, mss_irt
    use fileutils, only : set_filename, putfil
    use spmdMod, only : masterproc
    use time_manager, only : get_nstep, get_curr_date, get_curr_time, get_prev_date
    use shr_const_mod, only: SHR_CONST_CDAY
    use shr_sys_mod, only : shr_sys_flush
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES
    integer :: t                      !tape index 
    integer :: f                      !field index
    integer :: ier                    !error code
    integer :: nstep                  !current step
    integer :: day                    !current day (1 -> 31)
    integer :: mon                    !current month (1 -> 12)
    integer :: yr                     !current year (0 -> ...)
    integer :: mcsec                  !seconds of current date
    integer :: mdcur                  !current day
    integer :: mscur                  !seconds of current day 
    integer :: mcdate                 !current date 
    integer :: daym1                  !nstep-1 day (1 -> 31)
    integer :: monm1                  !nstep-1 month (1 -> 12)
    integer :: yrm1                   !nstep-1 year (0 -> ...)
    integer :: mcsecm1                !nstep-1 time of day [seconds]
    integer :: start1                 !start 1d value
    integer :: count1                 !count 1d value
    integer :: start2(2)              !start 2d values 
    integer :: count2(2)              !count 2d values 
    character(len=256) :: rem_dir     !remote (archive) directory
    character(len=256) :: rem_fn      !remote (archive) filename
    character(len=  8) :: cdate       !system date
    character(len=  8) :: ctime       !system time
    real(r8):: time                   !current time
    real(r8):: timedata(2)            !time interval boundaries
    logical :: if_stop                !true => last time step of run 
    logical :: if_disphist(max_tapes) !true => save and dispose history file 
    logical :: if_remvhist(max_tapes) !true => remove local history file after dispose
!-----------------------------------------------------------------------

    ! get current step

    nstep = get_nstep()

    ! Set calendar for current time step

    call get_curr_date (yr, mon, day, mcsec) 
    mcdate = yr*10000 + mon*100 + day

    ! Set calendar for current for previous time step
    
    call get_prev_date (yrm1, monm1, daym1, mcsecm1)
    
    ! Set elapased time since reference date
    
    call get_curr_time(mdcur, mscur)  
    
    ! Loop over active history tapes, create new history files if necessary 
    ! and write data to history files if end of history interval.
    
    do t = 1, ntapes
       
       ! Skip nstep=0 if monthly average
       
       if (nstep==0 .and. tape(t)%nhtfrq==0) cycle
       
       ! Determine if end of history interval
       
       tape(t)%is_endhist = .false.
       if (tape(t)%nhtfrq==0) then   !monthly average
          if (mon /= monm1) tape(t)%is_endhist = .true.
#if (defined COUP_TEM)
          if (mon /= monm1 .and. t == 2) ready4tem = .true.
#endif 
       else 
          if (mod(nstep,tape(t)%nhtfrq) == 0) tape(t)%is_endhist = .true.
       end if

       ! If end of history interval

       if (tape(t)%is_endhist) then

          ! Normalize history buffer if time averaged

          call hfields_normalize(t)

          ! Increment current time sample counter.  

          tape(t)%ntimes = tape(t)%ntimes + 1

          ! Create history file if appropriate and build time comment

          ! If first time sample, generate unique history file name, open file
          ! and write time constant fields to file

#if (defined COUP_TEM)
         if (t == 1 .or. t > 4) then
#endif
          if (tape(t)%ntimes == 1) then
             if (masterproc) then
                locfnh(t) = set_hist_filename (hist_freq=tape(t)%nhtfrq, hist_file=t)
                write(6,*)'HTAPES_WRAPUP: Creating history file ', trim(locfnh(t)), &
                     ' at nstep = ',get_nstep()
                call htape_create (t)
             endif

             ! Write time constant history variables to primary tape
             if (t==1) call htape_timeconst(t)
          endif

          ! Write time information to history file

          if (masterproc) then

             start1 = tape(t)%ntimes
             count1 = 1

             call get_curr_time (mdcur, mscur)  
             call get_curr_date (yr, mon, day, mcsec) 
             mcdate = yr*10000 + mon*100 + day

             call wrap_put_vara_int (nfid(t), mcdate_id(t), start1, count1, mcdate)
             call wrap_put_vara_int (nfid(t), mcsec_id(t), start1, count1, mcsec)
             call wrap_put_vara_int (nfid(t), mdcur_id(t), start1, count1, mdcur)
             call wrap_put_vara_int (nfid(t), mscur_id(t), start1, count1, mscur)
             call wrap_put_vara_int (nfid(t), nstep_id(t), start1, count1, nstep)

             time = mdcur + mscur/86400._r8
             call wrap_put_vara_realx (nfid(t), time_var_id(t), start1, count1, time)

             start2(1) = 1
             start2(2) = tape(t)%ntimes
             count2(1) = 2
             count2(2) = 1
             timedata(1) = tape(t)%begtime
             timedata(2) = time
             call wrap_put_vara_realx (nfid(t), time_bounds_id(t), start2, count2, timedata)
             
             start2(1) = 1
             start2(2) = tape(t)%ntimes
             count2(1) = 8
             count2(2) = 1
             call getdatetime (cdate, ctime)
             call wrap_put_vara_text (nfid(t), date_written_id(t), start2, count2, cdate)
             call wrap_put_vara_text (nfid(t), time_written_id(t), start2, count2, ctime)

             write(6,*)
             write(6,*)'HTAPES_WRAPUP: Writing current time sample to local history file ', &
                  trim(locfnh(t)),' at nstep = ',get_nstep(), &
                  ' for history time interval beginning at ', tape(t)%begtime, &
                  ' and ending at ',time
             write(6,*)
             call shr_sys_flush(6)

             ! Update beginning time of next interval

             tape(t)%begtime = time  

          endif

#if (defined COUP_TEM)
         endif 
#endif
          ! Write history time samples

          call hfields_write(t)

          ! Zero necessary history buffers

          call hfields_zero(t)

       end if

    end do  ! end loop over history tapes

    ! Determine if file needs to be closed and disposed
    
    call do_close_dispose (ntapes, tape(:)%ntimes, tape(:)%mfilt, if_stop, if_disphist, if_remvhist)
    
    ! Close open history file and dispose to mass store
    ! Auxilary files may have been closed and saved off without being full,
    ! must reopen the files

    do t = 1, ntapes
#if (defined COUP_TEM)
     if (t == 1 .or. t > 4) then
#endif
       if (if_disphist(t)) then
          if (masterproc) then
             if (tape(t)%ntimes /= 0) then
                write(6,*)
                write(6,*) 'HTAPES_WRAPUP: Closing local history file ',&
                     trim(locfnh(t)),' at nstep = ', get_nstep()
                write(6,*)
                call wrap_close(nfid(t))
                rem_dir = trim(archive_dir) // '/hist/'
                rem_fn = set_filename(rem_dir, locfnh(t))
                call putfil (locfnh(t), rem_fn, mss_wpass, mss_irt, if_remvhist(t))
                if (.not.if_stop .and. (tape(t)%ntimes/=tape(t)%mfilt)) then
                   call wrap_open (trim(locfnh(t)), NF_WRITE, nfid(t))
                end if
             else
                write(6,*)'HTAPES_WRAPUP: history tape ',t,': no open file to close'
             end if
          endif
       endif
#if (defined COUP_TEM)
      endif 
#endif
    end do

    ! Determine if time to write restarts before resetting number of
    ! time samples below

#if (defined OFFLINE) || (defined COUP_CAM)
    if_writrest = .false.
    if (tape(1)%ntimes == tape(1)%mfilt) if_writrest = .true.
#endif

    ! Reset number of time samples to zero if file is full

    do t = 1, ntapes
#if (defined COUP_TEM)
     if (t == 1 .or. t > 4) then
#endif
       if (if_disphist(t) .and. tape(t)%ntimes==tape(t)%mfilt) then
          tape(t)%ntimes = 0               
       end if
#if (defined COUP_TEM)
      endif 
#endif
    end do

  end subroutine htapes_wrapup

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restart_history
!
! !INTERFACE:
  subroutine restart_history (nio, flag)
!
! !DESCRIPTION: 
! Read/write history file restart data.
! If the current history file(s) are not full, file(s) are opened 
! so that subsequent time samples are added until the file is full. 
! A new history file is used on a branch run. 
!
! !USES
     use iobinary
     use clmtype, only : grid1d, land1d, cols1d, pfts1d
     use clm_varctl, only : archive_dir, nsrest, mss_irt
     use fileutils, only : set_filename, getfil
#if (defined SPMD)
     use spmdMod, only : masterproc, mpicom, MPI_REAL8, MPI_INTEGER, MPI_CHARACTER
#else
     use spmdMod, only : masterproc
#endif
!
! !ARGUMENTS:
     implicit none
     integer, intent(in) :: nio               !restart unit 
     character(len=*), intent(in) :: flag     !'read' or 'write'
!
! !REVISION HISTORY
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES
     integer :: k                             !1d index
     integer :: lev                           !level index  
     integer :: t                             !tape index
     integer :: f                             !field index
     integer :: ier                           !error status  
     character(len=256) :: fnameh(max_tapes)  !full name of history file
     character(len=256) :: rem_dir            !remote (archive) directory
     character(len=256) :: filename           !generic filename
     character(len=  8) :: type1d             !1d type
     integer :: num1d,beg1d,end1d             !1d size, beginning and ending indices  
     integer :: numlev                        !number of vertical levels 
     real(r8), pointer  :: hbuf(:,:)          !history buffer
     integer , pointer  :: nacs(:,:)          !accumulation counter
     integer , pointer, dimension(:)   :: ibuf1d  !temporary
     real(r8), pointer, dimension(:)   :: rbuf1d  !temporary
     integer , pointer, dimension(:,:) :: ibuf2d  !temporary
     real(r8), pointer, dimension(:,:) :: rbuf2d  !temporary
!-----------------------------------------------------------------------

     ! If branch run, initialize file times and return

     if (flag == 'read') then
        if (nsrest == 3) then
           do t = 1,ntapes
              tape(t)%ntimes = 0
           end do
           RETURN
        end if
     endif
     
     ! Read/write history file data only for restart run (not for branch run)

     if (masterproc) then

        if (flag == 'write') then
           write (nio, iostat=ier) ntapes, varid, fincl, fexcl
        else if (flag == 'read' ) then
           read (nio, iostat=ier) ntapes, varid, fincl, fexcl
        endif
        if (ier /= 0) then
           write (6,*)' RESTART_HISTORY read/write error 1',ier,' on i/o unit = ',nio
           call endrun
        end if

        do t = 1,ntapes
           if (flag == 'write') then
              write (nio, iostat=ier) &
                   tape(t)%nflds, tape(t)%ntimes, tape(t)%nhtfrq, tape(t)%mfilt, &
                   tape(t)%ncprec, tape(t)%begtime, tape(t)%is_endhist, tape(t)%dov2xy
           else if (flag == 'read') then
              read (nio, iostat=ier) &
                   tape(t)%nflds, tape(t)%ntimes, tape(t)%nhtfrq, tape(t)%mfilt, &
                   tape(t)%ncprec, tape(t)%begtime, tape(t)%is_endhist, tape(t)%dov2xy
           endif
           if (ier /= 0) then
              write (6,*)' RESTART_HISTORY read/write error 2',ier,' on i/o unit = ',nio
              call endrun
           end if
           do f=1,tape(t)%nflds
              if (flag == 'write') then
                 write (nio,iostat=ier) tape(t)%hlist(f)%field%name,        &
                                        tape(t)%hlist(f)%field%long_name,   &
                                        tape(t)%hlist(f)%field%units,       &
                                        tape(t)%hlist(f)%field%numlev,      &
                                        tape(t)%hlist(f)%field%hpindex,     &
                                        tape(t)%hlist(f)%field%hpindices,   &
                                        tape(t)%hlist(f)%field%type1d,      &
                                        tape(t)%hlist(f)%avgflag            
              else if (flag == 'read') then
                 read  (nio,iostat=ier) tape(t)%hlist(f)%field%name,        &
                                        tape(t)%hlist(f)%field%long_name,   &
                                        tape(t)%hlist(f)%field%units,       &
                                        tape(t)%hlist(f)%field%numlev,      &
                                        tape(t)%hlist(f)%field%hpindex,     &
                                        tape(t)%hlist(f)%field%hpindices,   &
                                        tape(t)%hlist(f)%field%type1d,      &
                                        tape(t)%hlist(f)%avgflag            
              end if
              if (ier /= 0) then
                 write(6,*)'HIST RESART read/write error 3: End or error condition writing ', &
                      'history restart field ',f,' from tape ',t
                 call endrun
              end if
           end do
        end do

     endif  ! end of if-masterproc block

#if ( defined SPMD )
     if (flag == 'read') then

        ! Broadcast history information from masterprocessor

        call mpi_bcast (ntapes, 1, mpi_integer, 0, mpicom, ier)
        do t = 1,ntapes
           call mpi_bcast (tape(t)%nflds, 1, MPI_INTEGER, 0, mpicom, ier)
           call mpi_bcast (tape(t)%ntimes, 1, MPI_INTEGER, 0, mpicom, ier)
           call mpi_bcast (tape(t)%nhtfrq, 1, MPI_INTEGER ,0,mpicom, ier)
           call mpi_bcast (tape(t)%mfilt, 1, MPI_INTEGER, 0, mpicom, ier)
           call mpi_bcast (tape(t)%ncprec, 1, MPI_INTEGER, 0, mpicom, ier)
           call mpi_bcast (tape(t)%is_endhist, 1, MPI_LOGICAL, 0, mpicom, ier)
           call mpi_bcast (tape(t)%dov2xy, 1, MPI_LOGICAL, 0, mpicom, ier)
           call mpi_bcast (tape(t)%begtime, 1, MPI_REAL8, 0, mpicom, ier)
           do f = 1,tape(t)%nflds
              call mpi_bcast (tape(t)%hlist(f)%field%name, 8, MPI_CHARACTER, 0, mpicom, ier)
              call mpi_bcast (tape(t)%hlist(f)%field%units, max_chars, MPI_CHARACTER, 0, mpicom, ier)
              call mpi_bcast (tape(t)%hlist(f)%field%numlev, 1, MPI_INTEGER, 0, mpicom, ier)
              call mpi_bcast (tape(t)%hlist(f)%field%hpindex, 1, MPI_INTEGER, 0, mpicom, ier)
              call mpi_bcast (tape(t)%hlist(f)%field%hpindices, 4, MPI_INTEGER, 0, mpicom, ier)
              call mpi_bcast (tape(t)%hlist(f)%field%type1d, 8, MPI_CHARACTER, 0, mpicom, ier)
              call mpi_bcast (tape(t)%hlist(f)%avgflag, 1, MPI_CHARACTER, 0, mpicom, ier)
           end do
        end do
     endif
#endif

     ! Allocate memory for history buffers - read only

     if (flag == 'read') then
        do t = 1,ntapes
           do f = 1,tape(t)%nflds
              numlev = tape(t)%hlist(f)%field%numlev
              type1d = tape(t)%hlist(f)%field%type1d
              select case (type1d)
              case ('gridcell') 
                 num1d = grid1d%num 
                 beg1d = grid1d%beg
                 end1d = grid1d%end 
              case ('landunit') 
                 num1d = land1d%num 
                 beg1d = land1d%beg
                 end1d = land1d%end
              case ('column')
                 num1d = cols1d%num 
                 beg1d = cols1d%beg
                 end1d = cols1d%end
              case ('pft')
                 num1d = pfts1d%num 
                 beg1d = pfts1d%beg
                 end1d = pfts1d%end
              case default
                 write(6,*)'RESTART_HISTORY: read unknown 1d vertical type=',type1d    
                 call endrun ()
              end select
              tape(t)%hlist(f)%field%num1d = num1d
              tape(t)%hlist(f)%field%beg1d = beg1d
              tape(t)%hlist(f)%field%end1d = end1d

              allocate (tape(t)%hlist(f)%hbuf(beg1d:end1d,numlev))
              tape(t)%hlist(f)%hbuf(:,:) = 0._r8
              allocate (tape(t)%hlist(f)%nacs(beg1d:end1d,numlev))
              tape(t)%hlist(f)%nacs(:,:) = 0
           end do
        end do
     endif

     ! Loop over tapes - only read/write accumulators and counters if needed
     ! (if not end of history interval)
     
     do t = 1,ntapes
        if (.not. tape(t)%is_endhist) then

           if (masterproc) then
              if (flag == 'write') then
                 write(6,*)'RESTART_HISTORY: Writing history restart information for tape ',t
              else
                 write(6,*)'RESTART_HISTORY: Reading history restart information for tape ',t
              endif
           endif
              
           do f=1,tape(t)%nflds
              
              type1d = tape(t)%hlist(f)%field%type1d
              numlev = tape(t)%hlist(f)%field%numlev
              num1d =  tape(t)%hlist(f)%field%num1d
              beg1d =  tape(t)%hlist(f)%field%beg1d
              end1d =  tape(t)%hlist(f)%field%end1d
              nacs  => tape(t)%hlist(f)%nacs
              hbuf  => tape(t)%hlist(f)%hbuf

              if (numlev == 1) then
                 allocate(ibuf1d(num1d))
                 allocate(rbuf1d(num1d))
                 if (flag == 'read' ) then
                    call readin (nio, ibuf1d, clmlevel=type1d)
                    call readin (nio, rbuf1d, clmlevel=type1d)
                    do k = beg1d,end1d
                       nacs(k,1) = ibuf1d(k)
                       hbuf(k,1) = rbuf1d(k)
                    end do
                 else if (flag == 'write') then
                    do k = beg1d,end1d
                       ibuf1d(k) = nacs(k,1)
                       rbuf1d(k) = hbuf(k,1)
                    end do
                    call wrtout (nio, ibuf1d, clmlevel=type1d)
                    call wrtout (nio, rbuf1d, clmlevel=type1d)
                 endif
                 deallocate(ibuf1d)
                 deallocate(rbuf1d)
              else
                 allocate(ibuf2d(numlev,num1d))
                 allocate(rbuf2d(numlev,num1d))
                 if (flag == 'read') then
                    call readin (nio, ibuf2d, clmlevel=type1d)
                    call readin (nio, rbuf2d, clmlevel=type1d)
                    do lev = 1,numlev
                       do k = beg1d,end1d
                          nacs(k,lev) = ibuf2d(lev,k)
                          hbuf(k,lev) = rbuf2d(lev,k)
                       end do
                    end do
                 else if (flag == 'write') then
                    do lev = 1,numlev
                       do k = beg1d,end1d
                          ibuf2d(lev,k) = nacs(k,lev)
                          rbuf2d(lev,k) = hbuf(k,lev)
                       end do
                    end do
                    call wrtout (nio, ibuf2d, clmlevel=type1d)
                    call wrtout (nio, rbuf2d, clmlevel=type1d)
                 endif
                 deallocate(ibuf2d)
                 deallocate(rbuf2d)
              endif
              
           end do   ! end of fields do-loop
        endif   ! end of if-end-of-history-interval block
     end do   ! end of tape do-loop

     ! Read names of history files. If history file is not full,
     ! open history file and obtain time dependent netCDF variable id's

     if (flag == 'read') then
        if (masterproc) then
           do t = 1,ntapes
              read (nio) fnameh(t)
              if (tape(t)%ntimes /= 0) then
                 call getfil (fnameh(t), locfnh(t), 0)
                 call wrap_open (locfnh(t), nf_write, nfid(t))
                 call wrap_inq_varid (nfid(t), 'mcdate', mcdate_id(t))
                 call wrap_inq_varid (nfid(t), 'mcsec' , mcsec_id(t))
                 call wrap_inq_varid (nfid(t), 'mdcur' , mdcur_id(t))
                 call wrap_inq_varid (nfid(t), 'mscur' , mscur_id(t))
                 call wrap_inq_varid (nfid(t), 'nstep' , nstep_id(t))
                 call wrap_inq_varid (nfid(t), 'time'  , time_var_id(t))
                 call wrap_inq_varid (nfid(t), 'time_bounds' , time_bounds_id(t))
                 call wrap_inq_varid (nfid(t), 'date_written', date_written_id(t))
                 call wrap_inq_varid (nfid(t), 'time_written', time_written_id(t))
              end if
           end do
        endif
     endif   
     
     ! Write name of current history file(s) 
     
     if (flag == 'write') then
        if (masterproc) then
           do t = 1,ntapes  
              if (mss_irt == 0) then
                 filename = locfnh(t)
              else
                 rem_dir = trim(archive_dir) // '/hist/'
                 filename = set_filename(rem_dir, locfnh(t))
              endif
              write(nio) filename
           end do
        endif
     endif

     return
   end subroutine restart_history

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: getname
!
! !INTERFACE:
   character(len=max_namlen) function getname (inname)
!
! !DESCRIPTION: 
! Retrieve name portion of inname. If an averaging flag separater character is 
! present (":") in inname, lop it off
!
! !ARGUMENTS:
   implicit none
   character(len=*), intent(in) :: inname
!
! !REVISION HISTORY
! Created by Jim Rosinski
!
!EOP
!
! LOCAL VARIABLES
   integer :: length
   integer :: i
!-----------------------------------------------------------------------

   length = len (inname)
   
   if (length < max_namlen .or. length > max_namlen+2) then
      write(6,*) 'getname: bad length=',length
      call endrun
   end if
   
   getname = ' '
   do i = 1,max_namlen
      if (inname(i:i) == ':') exit
      getname(i:i) = inname(i:i)
   end do
   
   return
   end function getname

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: getflag
!
! !INTERFACE:
   character(len=1) function getflag (inname)
!
! !DESCRIPTION: 
! Retrieve flag portion of inname. If an averaging flag separater character 
! is present (":") in inname, return the character after it as the flag
!
! !ARGUMENTS:
   implicit none
   character(len=*) inname   ! character string
!
! !REVISION HISTORY
! Created by Jim Rosinski
!
!EOP
!
! LOCAL VARIABLES
   integer :: length         ! length of inname
   integer :: i              ! loop index
!-----------------------------------------------------------------------

   length = len (inname)

   if (length < max_namlen .or. length > max_namlen+2) then
      write(6,*) 'getname: bad length=',length
      call endrun
   end if
   
   getflag = ' '
   do i = 1,length
      if (inname(i:i) == ':') then
         getflag = inname(i+1:i+1)
         exit
      end if
   end do
   
   return
   end function getflag

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: list_index
!
! !INTERFACE:
   subroutine list_index (list, name, index)
!
! !DESCRIPTION: 
!
! !USES
!
! !ARGUMENTS:
     implicit none
     character(len=*), intent(in) :: list(max_flds)  ! input list of names, possibly ":" delimited
     character(len=max_namlen), intent(in) :: name   ! name to be searched for
     integer, intent(out) :: index                   ! index of "name" in "list"
!
! !REVISION HISTORY
! Created by Jim Rosinski
!
!EOP
!
! LOCAL VARIABLES
     character(len=max_namlen) :: listname           ! input name with ":" stripped off.
     integer f                                       ! field index
!-----------------------------------------------------------------------
     
     ! Only list items
     
     index = 0
     do f=1,max_flds  
        listname = getname (list(f))
        if (listname == ' ') exit
        if (listname == name) then
           index = f
           exit
        end if
     end do
     
     return
   end subroutine list_index

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_hist_filename
!
! !INTERFACE:
  character(len=256) function set_hist_filename (hist_freq, hist_file)
!
! !DESCRIPTION: 
! Determine history dataset filenames
!
! !USES
    use clm_varctl, only : caseid
    use time_manager, only : get_curr_date, get_prev_date
!
! !ARGUMENTS:
   implicit none
   integer, intent(in)  :: hist_freq   !history file frequency
   integer, intent(in)  :: hist_file   !history file index 
!
! !REVISION HISTORY
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES
   character(len=256) :: cdate       !date char string
   character(len=  1) :: hist_index  !p,1 or 2 (currently)
   integer :: day                    !day (1 -> 31)
   integer :: mon                    !month (1 -> 12)
   integer :: yr                     !year (0 -> ...)
   integer :: sec                    !seconds into current day
!-----------------------------------------------------------------------

   if (hist_freq == 0 ) then   !monthly
      call get_prev_date (yr, mon, day, sec) 
      write(cdate,'(i4.4,"-",i2.2)') yr,mon
   else                        !other 
      call get_curr_date (yr, mon, day, sec) 
      write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,sec
   endif
   write(hist_index,'(i1.1)') hist_file - 1
   set_hist_filename = "./"//trim(caseid)//".clm2.h"//hist_index//"."//&
        trim(cdate)//".nc"
   
  end function set_hist_filename

end module histFileMod


