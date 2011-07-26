#include <misc.h>
#include <preproc.h>

module atmdrvMod

#if (defined OFFLINE)

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: atmdrvMod
!
! !DESCRIPTION: 
! Read and generate atmospheric grid data at model resolution
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_const_mod, only : SHR_CONST_TKFRZ, SHR_CONST_PSTD
  use clm_varpar, only : lsmlon, lsmlat
#if (defined PCP2PFT)
  use pcp2pftMod, only : pcp2land
#endif
#if (defined SPMD)
  use spmdMod, only : masterproc, mpicom, MPI_REAL8, MPI_INTEGER
#else
  use spmdMod, only : masterproc
#endif
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: atmdrv       ! read atmospheric data
  public :: atm_getgrid  ! read atmospheric grid

!                              
! !REVISION HISTORY:
! Created by Gordon Bonan, Sam Levis and Mariana Vertenstein
!
!EOP
!
! PRIVATE MEMBER FUNCTIONS:
  private :: atm_openfile  ! open atmospheric forcing netCDF file
  private :: atm_readdata  ! read atmospheric forcing data
  private :: interpa2si    ! initialize atm->land model interpolation
  private :: interpa2s     ! area average fields from atmosphere to surface grid
!
! PRIVATE TYPES:
  private
!
! logical variables for file manipuation
!
  logical :: open_data=.true.             !true => open data file (first tstep of the run or month)
  logical :: allocated_data=.false.       !true => allocate dynamic data
!
! atmospheric grid data
!
#if (defined COUP_MIT2D)
  real(r8) :: psmit (lsmlon,lsmlat)
  real(r8) :: pcplmit (lsmlon,lsmlat)
  real(r8) :: pcpcmit (lsmlon,lsmlat)
  real(r8) :: tprmit (lsmlon,lsmlat)
  real(r8) :: tslmit (lsmlon,lsmlat)
  real(r8) :: qsmit (lsmlon,lsmlat)
  real(r8) :: wsmit (lsmlon,lsmlat)
  real(r8) :: usmit (lsmlon,lsmlat)
  real(r8) :: vsmit (lsmlon,lsmlat)
  real(r8) :: dswmit (lsmlon,lsmlat)
  real(r8) :: dlwmit (lsmlon,lsmlat)
  real(r8) :: pco2mit (lsmlon,lsmlat)
  real(r8) :: swnirmit (lsmlon,lsmlat)
  real(r8) :: swparmit (lsmlon,lsmlat)
  real(r8) :: tglbtrnd 
#if (defined PCPTRND)
  common/mit2din/psmit,pcplmit,pcpcmit,tprmit,tslmit,qsmit,wsmit,usmit,vsmit,dswmit,dlwmit,pco2mit,swnirmit,swparmit,tglbtrnd
#else
  common/mit2din/psmit,pcplmit,pcpcmit,tprmit,tslmit,qsmit,wsmit,usmit,vsmit,dswmit,dlwmit,pco2mit,swnirmit,swparmit
#endif
#endif
  integer  :: atmlon                      !number of atm longitudes
  integer  :: atmlat                      !number of atm latitudes
  real(r8) :: edge_a(4)                   !N,E,S,W edges of atm grid
  integer , allocatable :: numlon_a(:)    !number of lon points at each lat
  real(r8), allocatable :: latixy_a(:,:)  !latitude of grid cell (degrees)
  real(r8), allocatable :: longxy_a(:,:)  !longitude of grid cell (degrees)
!
! atmospheric forcing variables on atmospheric grid
!
  real(r8), allocatable :: x(:,:,:)            !temp. array in which atm data is stored
  real(r8), allocatable :: forc_txy_a (:,:)    !atm bottom level temperature (Kelvin)
  real(r8), allocatable :: forc_uxy_a (:,:)    !atm bottom level zonal wind (m/s)
  real(r8), allocatable :: forc_vxy_a (:,:)    !atm bottom level meridional wind (m/s)
  real(r8), allocatable :: forc_qxy_a (:,:)    !atm bottom level specific humidity (kg/kg)
  real(r8), allocatable :: zgcmxy_a (:,:)      !atm bottom level height above surface (m)
  real(r8), allocatable :: prcxy_a  (:,:)      !convective precipitation rate (mm H2O/s)
  real(r8), allocatable :: prlxy_a  (:,:)      !large-scale precipitation rate (mm H2O/s)
  real(r8), allocatable :: flwdsxy_a(:,:)      !downward longwave rad onto surface (W/m**2)
  real(r8), allocatable :: fswdsxy_a(:,:)      !downward longwave rad onto surface (W/m**2)
  real(r8), allocatable :: forc_solsxy_a (:,:) !vis direct beam solar rad onto srf (W/m**2)
  real(r8), allocatable :: forc_sollxy_a (:,:) !nir direct beam solar rad onto srf (W/m**2)
  real(r8), allocatable :: forc_solsdxy_a(:,:) !vis diffuse solar rad onto srf (W/m**2)
  real(r8), allocatable :: forc_solldxy_a(:,:) !nir diffuse solar rad onto srf(W/m**2)
  real(r8), allocatable :: forc_pbotxy_a (:,:) !atm bottom level pressure (Pa)
  real(r8), allocatable :: forc_psrfxy_a (:,:) !atm surface pressure (Pa)
!
! atmospheric forcing variables on land model grid
!
  real(r8) :: forc_txy (lsmlon,lsmlat)         !atm bottom level temperature (Kelvin)
  real(r8) :: forc_uxy (lsmlon,lsmlat)         !atm bottom level zonal wind (m/s)
  real(r8) :: forc_vxy (lsmlon,lsmlat)         !atm bottom level meridional wind (m/s)
  real(r8) :: forc_qxy (lsmlon,lsmlat)         !atm bottom level specific humidity (kg/kg)
  real(r8) :: zgcmxy (lsmlon,lsmlat)           !atm bottom level height above surface (m)
  real(r8) :: prcxy  (lsmlon,lsmlat)           !convective precipitation rate (mm H2O/s)
  real(r8) :: prlxy  (lsmlon,lsmlat)           !large-scale precipitation rate (mm H2O/s)
  real(r8) :: flwdsxy(lsmlon,lsmlat)           !downward longwave rad onto surface (W/m**2)
  real(r8) :: forc_solsxy (lsmlon,lsmlat)      !vis direct beam solar rad onto srf (W/m**2)
  real(r8) :: forc_sollxy (lsmlon,lsmlat)      !nir direct beam solar rad onto srf (W/m**2)
  real(r8) :: forc_solsdxy(lsmlon,lsmlat)      !vis diffuse solar rad onto srf (W/m**2)
  real(r8) :: forc_solldxy(lsmlon,lsmlat)      !nir diffuse solar rad onto srf(W/m**2)
  real(r8) :: forc_pbotxy (lsmlon,lsmlat)      !atm bottom level pressure (Pa)
  real(r8) :: forc_psrfxy (lsmlon,lsmlat)      !atm surface pressure (Pa)
!CAS
!CAS	Stochastic treatment of precipitation added for MIT2D model implementation/coupling
!CAS
!CAS    C. Adam Schlosser, June, 2004
!CAS
#if (defined STOCHASTIC)
    real,save :: t_dry(lsmlon,lsmlat)
    real,save :: t_storm(lsmlon,lsmlat)
    real,save :: pcpc_resid(lsmlon,lsmlat)
    real,save :: pcpl_resid(lsmlon,lsmlat)
    real,save :: prc_poiss(lsmlon,lsmlat)
    real,save :: prl_poiss(lsmlon,lsmlat)
    integer,save :: storms(lsmlon,lsmlat)
    integer,save :: storm_dur(lsmlon,lsmlat)
    integer,save :: dtcumu(lsmlon,lsmlat)
    integer :: ist,jst
    integer,parameter :: ngrids=lsmlon*lsmlat
    data storms/ngrids*0/
#endif
!
! atmosphere grid to land model surface grid mapping for each land grid cell:
!
  integer, parameter :: mxovr =10          !maximum number of overlapping cells
  integer :: novr_a2s(lsmlon,lsmlat)       !number    of overlapping atm cells
  integer :: iovr_a2s(lsmlon,lsmlat,mxovr) !lon index of overlapping atm cells
  integer :: jovr_a2s(lsmlon,lsmlat,mxovr) !lat index of overlapping atm cells
  real(r8):: wovr_a2s(lsmlon,lsmlat,mxovr) !weight    of overlapping atm cells
!
! file netCDF id's
!
  integer :: ncid                !netCDF dataset id
  integer :: nvar                !number of variables in the data file
  integer :: nlon                !number of atm longitude points
  integer :: nlat                !number of atm latitude points
  integer :: ntim                !number of atm time slices per data file
  character(len=8) :: varnam(99) !variable names of atm. fields
!----------------------------------------------------------------------- 

contains

!=======================================================================

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: atmdrv
!
! !INTERFACE:
  subroutine atmdrv(nstep)
!
! !DESCRIPTION: 
! This code reads in atmospheric fields from an input file and generates 
! the required atmospheric forcing. These data files have [atmmin] minute 
! average data for each month. Input data files are named in month-year 
! format (e.g., 09-0001 contains 240 3-hour time slices of data, 30*8, for 
! September of year one). The model will cycle through however many full 
! years of data are available [pyr]. At least one full year of data is
! necessary for cycling. The model may start on any base date, as long as
! this date corresponds to an existing data file. A run need not be an
! exact multiple of a year.
!
! ============================
! Possible atmospheric fields:
! ============================
! Name     Description                              Required/Optional
! -----------------------------------------------------------------------------
! TBOT     temperature (K)                          Required
! WIND     wind:sqrt(u**2+v**2) (m/s)               Required
! QBOT     specific humidity (kg/kg)                Required
! Tdew     dewpoint temperature (K)                 Alternative to Q
! RH       relative humidity (percent)              Alternative to Q
! ZBOT     reference height (m)                     optional 
! PSRF     surface pressure (Pa)                    optional
! FSDS     total incident solar radiation (W/m**2)  Required
! FSDSdir  direct incident solar radiation (W/m**2) optional (replaces FSDS)
! FSDSdif  diffuse incident solar rad (W/m**2)      optional (replaces FSDS)
! FLDS     incident longwave radiation (W/m**2)     optional
! PRECTmms total precipitation (mm H2O / sec)       Required
! PRECCmms convective precipitation (mm H2O / sec)  optional (replaces PRECT)
! PRECLmms large-scale precipitation (mm H2O / sec) optional (replaces PRECT)
! 
! ============
! Data format:
! ============
! Data format is netCDF with dimensions longitude x latitude
! for each time slice and field. Variable names can be as in above list
! or can be reset to desired names using [fldlst] in code below.
!
! ===============
! Namelist input:
! ===============
! character*256 offline_atmdir = directory for input atm data files (can be Mass Store)
!
! !USES:
    use nanMod
    use clmtype
    use clm_varctl, only : offline_atmdir
    use clm_varcon, only : rair, cpair, po2, pco2, tcrit, spval
#if (defined STOCHASTIC)
    use clm_varsur, only : exp_tstorm, exp_tdry, trnd_tdry
#endif
    use time_manager, only : get_step_size, get_curr_calday, get_curr_date
    use fileutils, only : getfil
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    integer, intent(in) :: nstep    !current time step
!
! !REVISION HISTORY:
! Created by Sam Levis
!
!EOP
!
! LOCAL VARIABLES:
#if (defined STOCHASTIC)
    real(r8):: trn_tdry               !Transient Value of Expected Storm Interval
#endif
    integer :: i,j,k,gi,gf            !indices
    integer :: itimlast               !last time index used in atmrd
    real(r8):: calday                 !calendar day at Greenwich (1.00 -> 365.99)
    integer :: kda                    !day (1 -> 31)
    integer :: kmo                    !month (1 -> 12)
    integer :: kyr                    !year (0 -> ...)
    integer :: ksec                   !current seconds of current date (0 -> 86400)
    integer :: mcdate                 !current date in integer format [yyyymmdd]
    integer :: mcsec                  !current time of day [seconds]
    integer :: dtime                  !time step size  
    integer :: minpday = 1440         !minutes per day
    integer :: secpmin = 60           !seconds per minute
    integer, SAVE :: itim             !time index used in atmrd
    integer, SAVE :: atmmin           !temporal resolution of atm data (in minutes)
    character(len=256), SAVE :: locfn !full file name in case atmdir is in MSS
    type(atm2lnd_state_type), pointer :: a2ls !local pointer to derived subtype
    type(atm2lnd_flux_type) , pointer :: a2lf !local pointer to derived subtype
!------------------------------------------------------------------------


    ! -----------------------------------------------------------------
    ! Open netcdf file and read data every [atmmin] minutes
    ! -----------------------------------------------------------------

    ! Calendar information for current [nstep]

    dtime = get_step_size()
    calday = get_curr_calday()
    call get_curr_date(kyr, kmo, kda, mcsec)
    mcdate = kyr*10000 + kmo*100 + kda

    ! If 1st tstep of the run or of the month, then enter next if-block
    ! Rest flag to open file and set flag to read data

    if (open_data) then
#if (!defined COUP_MIT2D)
    locfn = " "
    call atm_openfile (kda, kmo, kyr, locfn, itim, atmmin)
#else
    locfn = trim(offline_atmdir)
    atmmin = dtime/60.
    itim=0
#endif
    endif

    ! Calculate time index

!    write (6,*) 'ATMDRV: ',itim
!    write (6,*) 'ATMDRV: ',atmmin
    itimlast = itim
    itim = 1 + ((kda - 1)*minpday + (mcsec - dtime)/secpmin)/atmmin
    if (dtime == int(atmmin*secpmin)) then   
       if (kda == 1 .and. mcsec == 0) itim = itimlast+1
    else
       if (kda == 1 .and. mcsec == 0) itim = itimlast
    endif

    ! Determine if new data is to be read

    if (open_data .or. mod(nstep-1,atmmin*secpmin/dtime) ==0 ) then

       ! Read data for current time slice

!       if (masterproc) then
!          write (6,*)
!          write (6,'(72a1)') ("-",i=1,60)
!          write (6,*)'nstep= ',nstep,' date= ',mcdate,' sec= ',mcsec
!          if ( len_trim(locfn) > 0 )then
!             write (6,*)'ATM: attempting to read data from ',trim(locfn)
!          end if
!          write (6,'(72a1)') ("-",i=1,60)
!          write (6,*)
!       endif
#if (defined COUP_MIT2D)
       forc_txy=tslmit
       forc_uxy=usmit
       forc_vxy=vsmit
       forc_qxy=qsmit
       prcxy=pcpcmit/3600.
       prlxy=pcplmit/3600.
       flwdsxy=dlwmit
       forc_solsdxy=0.3*swparmit
       forc_solsxy=0.7*swparmit
       forc_solldxy=0.3*swnirmit
       forc_sollxy=0.7*swnirmit
       forc_pbotxy=psmit*100
       forc_psrfxy=psmit*100
       zgcmxy(:,:) = 30.
#else
       call atm_readdata (locfn, kmo, itim)

       ! Map 2d atmospheric fields from atmospheric grid to land model surface grid. 
       ! Area-average absolute value of winds (i.e., regardless of
       ! direction) since land model cares about magnitude not direction.
       ! Then need to adjust resultant stresses for direction of wind.

       call interpa2s (forc_txy_a    , forc_txy    , zgcmxy_a      , zgcmxy      , &
                       forc_uxy_a    , forc_uxy    , forc_vxy_a    , forc_vxy    , &
                       forc_qxy_a    , forc_qxy    , prcxy_a       , prcxy       , &
                       prlxy_a       , prlxy       , flwdsxy_a     , flwdsxy     , &
                       forc_solsxy_a , forc_solsxy , forc_sollxy_a , forc_sollxy , &
                       forc_solsdxy_a, forc_solsdxy, forc_solldxy_a, forc_solldxy, &
                       forc_pbotxy_a , forc_pbotxy , forc_psrfxy_a , forc_psrfxy )

#endif
#if (defined PCP2PFT)
!
!       Adjust Zonal Precipitation to Overland Precipitation
!
!CAS       prlxy = pcp2land*prlxy
!CAS       prcxy = pcp2land*prcxy
!CAS
!CAS	This adjustment has been moved to the IGSM driver.F routine - 06/01/05
!CAS
#endif
!CAS
!CAS	Stochastic treatment of precipitation added for MIT2D model implementation/coupling
!CAS
!CAS    C. Adam Schlosser, June, 2004
!CAS
!CAS    print value of global temperature anomaly for PCP Freq. Trend
!CAS
#if (defined STOCHASTIC)
       if (kda == 1 .and. mcsec == 0) print *, 'Global Land Temperature Trend for Precip. Freq.:', tglbtrnd

          do j=1,lsmlat
          do i=1,lsmlon
!
!	Check for beginning of month, if so, reset storm count
!

       if (kda == 1 .and. mcsec == 0) then
          storms = 0
       endif

!
!	Scale the expected value of tdry according to global temperature trend 
!	and normalized trend in tdry
!
!	NOTE: Formula assumes that normalized trend is in units of %/deg. C or K
!
#if (defined PCPTRND)

          trn_tdry = exp_tdry(i,j,kmo)*(1+(trnd_tdry(i,j,kmo)*tglbtrnd/100.))

!          write (6,*) i,j,kmo
!          write (6,*) trn_tdry,trnd_tdry(i,j,kmo)
!          write (6,*) prcxy(i,j),prlxy(i,j)
!
          call stoch(prcxy(i,j),prlxy(i,j),prc_poiss(i,j),prl_poiss(i,j), &
                     storm_dur(i,j),exp_tstorm(i,j,kmo),trn_tdry,&
                     t_storm(i,j),t_dry(i,j),dtcumu(i,j),storms(i,j),pcpc_resid(i,j),pcpl_resid(i,j))

#else

          call stoch(prcxy(i,j),prlxy(i,j),prc_poiss(i,j),prl_poiss(i,j), &
                     storm_dur(i,j),exp_tstorm(i,j,kmo),exp_tdry(i,j,kmo), &
                     t_storm(i,j),t_dry(i,j),dtcumu(i,j),storms(i,j),pcpc_resid(i,j),pcpl_resid(i,j))

#endif

!            write (6,*)'ATM: stochastic pcp',lsmlat,lsmlon
!            write (6,*) prcxy(i,j),prlxy(i,j)
!            write (6,*) prc_poiss(i,j),prl_poiss(i,j)
!            write (6,*) storm_dur(i,j),t_storm(i,j),t_dry(i,j),dtcumu(i,j)
!            write (6,*) pcpc_resid(i,j),pcpl_resid(i,j)
          enddo
          enddo
#endif

       ! Map atmospheric fields to force land model: [lsmlon] x [lsmlat] grid ->
       ! [numland] vector of land points -> [numpatch] vector of subgrid patches

!$OMP PARALLEL DO PRIVATE (gi,i,j,a2ls,a2lf)
       do gi = 1,clm%mps%ngridcells
          i = clm%g(gi)%gps%ixy
          j = clm%g(gi)%gps%jxy
          a2ls => clm%g(gi)%a2ls
          a2lf => clm%g(gi)%a2lf

          !States
          a2ls%forc_t = forc_txy(i,j)
          a2ls%forc_u = forc_uxy(i,j)
          a2ls%forc_v = forc_vxy(i,j)
          a2ls%forc_wind = sqrt(forc_uxy(i,j)**2 + forc_vxy(i,j)**2)
          a2ls%forc_q = forc_qxy(i,j)
          a2ls%forc_hgt = zgcmxy(i,j)
          a2ls%forc_hgt_u = zgcmxy(i,j) !observational height of wind [m] 
          a2ls%forc_hgt_t = zgcmxy(i,j) !observational height of temp [m] 
          a2ls%forc_hgt_q = zgcmxy(i,j) !observational height of humidity [m] 
          a2ls%forc_pbot = forc_pbotxy(i,j)
          a2ls%forc_psrf = forc_psrfxy(i,j)
          a2ls%forc_th  = a2ls%forc_t * (a2ls%forc_psrf/a2ls%forc_pbot)**(rair/cpair)
          a2ls%forc_vp  = a2ls%forc_q*a2ls%forc_pbot / (0.622+0.378*a2ls%forc_q)
          a2ls%forc_rho = (a2ls%forc_pbot-0.378*a2ls%forc_vp) / (rair*a2ls%forc_t)
          a2ls%forc_co2 = pco2*a2ls%forc_pbot
          a2ls%forc_o2  = po2*a2ls%forc_pbot

          !Fluxes
          a2lf%forc_lwrad = flwdsxy(i,j)
          a2lf%forc_solad(1) = forc_solsxy(i,j)
          a2lf%forc_solad(2) = forc_sollxy(i,j)
          a2lf%forc_solai(1) = forc_solsdxy(i,j)
          a2lf%forc_solai(2) = forc_solldxy(i,j)
          a2lf%forc_solar = forc_solsxy(i,j) + forc_sollxy(i,j) + forc_solsdxy(i,j) + forc_solldxy(i,j)

#if (defined STOCHASTIC)
          if (prc_poiss(i,j)+prl_poiss(i,j)==0) then
            a2lf%forc_strm_int = spval
            a2lf%forc_strm_dur = spval
          else
            a2lf%forc_strm_int = prc_poiss(i,j)+prl_poiss(i,j)
            a2lf%forc_strm_dur = storm_dur(i,j)*dtime
          endif
          a2lf%forc_strm_dry = t_dry(i,j)
          a2lf%forc_strms = storms(i,j)
#endif
  

          ! Set upper limit of air temperature for snowfall at 275.65K.
          ! This cut-off was selected based on Fig. 1, Plate 3-1, of Snow
          ! Hydrology (1956).

          if ( (prcxy(i,j)+prlxy(i,j)) > 0.) then
             if (a2ls%forc_t > (SHR_CONST_TKFRZ + tcrit)) then
                a2ls%itypprc = 1
                a2lf%forc_rain = prcxy(i,j) + prlxy(i,j)
                a2lf%forc_snow = 0.
             else
                a2ls%itypprc = 2
                a2lf%forc_rain = 0.
                a2lf%forc_snow = prcxy(i,j) + prlxy(i,j)
             endif
          else
             a2ls%itypprc = 0
             a2lf%forc_rain = 0.
             a2lf%forc_snow = 0
          endif
       end do

    end if

    ! Reset open_data 

    if (open_data) then
       open_data = .false.    !reset to false 
    elseif (kda == 1 .and. mcsec == 0) then
       open_data = .true.     !for next time step
    endif

    return
  end subroutine atmdrv

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: atm_getgrid
!
! !INTERFACE:
  subroutine atm_getgrid()
!
! !DESCRIPTION: 
! Read atmospheric grid 
!
! !USES:
    use nanMod
    use clm_varctl, only : offline_atmdir
    use clm_varcon, only : rair, cpair, po2, pco2
    use clm_varsur, only : numlon, longxy, latixy, lsmedge
    use fileutils, only : getfil
    use time_manager, only : get_curr_date
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer :: kda                !day (1 -> 31)
    integer :: kmo                !month (1 -> 12)
    integer :: kyr                !year (0 -> ...)
    integer :: ksec               !current seconds of current date (0 -> 86400)
    integer :: mcsec              !current time of day [seconds]
    character(len=  7) :: ext     !month-year extension, e.g., 01-0005
    character(len=256) :: filenam !full file name, atmdir + ext
    character(len=256) :: locfn   !full file name in case atmdir is in MSS
    logical :: lexist             !true => file exists, used when looking for a file
    integer :: dimid              !netCDF dimension id
    integer :: varid              !netCDF variable id
    integer :: ier                !error status
!------------------------------------------------------------------------

    ! ----------------------------------------------------------------------
    ! Read offline grid data and allocate dynamic memory
    ! ----------------------------------------------------------------------

    if (masterproc) then

       ! Build [month]-[year] extension for file name to be read
       ! append extension to path name to get full file name
#if (!defined COUP_MIT2D)
       call get_curr_date(kyr, kmo, kda, mcsec)
       write (ext,'(i4.4,"-",i2.2)') kyr,kmo
       filenam = trim(offline_atmdir) // '/' // ext // '.nc' 
       call getfil(filenam, locfn, 1)
#else
       locfn = trim(offline_atmdir) 
       write(6,*) 'ATM_GETGRID: atm datafile - ',locfn
#endif
       inquire (file = locfn, exist = lexist)
       if (.not. lexist) then
          write(6,*) 'ATM_GETGRID error: could not find initial atm datafile' 
          call endrun
       endif

       ! Open netCDF data file and get lengths of lat,lon,time dimensions

       call wrap_open (locfn, nf_nowrite, ncid)
       
       call wrap_inq_dimid  (ncid, 'lon', dimid)
       call wrap_inq_dimlen (ncid, dimid, atmlon)
       
       call wrap_inq_dimid  (ncid, 'lat', dimid)
       call wrap_inq_dimlen (ncid, dimid, atmlat)
       
       call wrap_inq_dimid  (ncid, 'time', dimid)
       call wrap_inq_dimlen (ncid, dimid, ntim)
       if (ntim == 0) then
          write (6,*) 'ATM_GETGRID error: zero input time slices'
          call endrun
       end if
       
    endif  !end of if-masterproc block
#if (defined SPMD) 
    call mpi_bcast (atmlon, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (atmlat, 1, MPI_INTEGER, 0, mpicom, ier)
#endif

    ! Allocate space for dynamic variables

    if (.not. allocated_data) then
       if (masterproc) write(6,*)' ATM_GETRID: allocating dynamic space'
       allocate (numlon_a(atmlat))      
       allocate (latixy_a(atmlon,atmlat))    
       allocate (longxy_a(atmlon,atmlat))    
       allocate (forc_txy_a(atmlon,atmlat))   
       allocate (forc_uxy_a(atmlon,atmlat))   
       allocate (forc_vxy_a(atmlon,atmlat))   
       allocate (forc_qxy_a(atmlon,atmlat))   
       allocate (zgcmxy_a(atmlon,atmlat))   
       allocate (prcxy_a(atmlon,atmlat))   
       allocate (prlxy_a(atmlon,atmlat))   
       allocate (flwdsxy_a(atmlon,atmlat))   
       allocate (forc_solsxy_a(atmlon,atmlat))   
       allocate (forc_sollxy_a(atmlon,atmlat))   
       allocate (forc_solsdxy_a(atmlon,atmlat))   
       allocate (forc_solldxy_a(atmlon,atmlat))   
       allocate (forc_pbotxy_a(atmlon,atmlat))   
       allocate (forc_psrfxy_a(atmlon,atmlat))   
       allocate (x(atmlon,atmlat,14))
       allocated_data = .true.
    endif

    ! Extract atmospheric data grid information and close file

    if (masterproc) then
       
       call wrap_inq_varid(ncid, 'EDGEN', varid)
       call wrap_get_var_realx(ncid, varid, edge_a(1))
       
       call wrap_inq_varid(ncid, 'EDGEE', varid)
       call wrap_get_var_realx(ncid, varid, edge_a(2))
       
       call wrap_inq_varid(ncid, 'EDGES', varid)
       call wrap_get_var_realx(ncid, varid, edge_a(3))
       
       call wrap_inq_varid(ncid, 'EDGEW', varid)
       call wrap_get_var_realx(ncid, varid, edge_a(4))
       
       call wrap_inq_varid(ncid, 'LONGXY', varid)
       call wrap_get_var_realx(ncid, varid, longxy_a)
       
       call wrap_inq_varid(ncid,'LATIXY', varid)
       call wrap_get_var_realx(ncid, varid, latixy_a)
       
       call wrap_close (ncid)
       write (6,*) 'ATM_GETGRID: closing data for ',trim(locfn)
       
    endif !end of if-masterproc block

#if (defined SPMD) 
    call mpi_bcast (edge_a  , size(edge_a)  , MPI_REAL8 , 0, mpicom, ier)
    call mpi_bcast (longxy_a, size(longxy_a), MPI_REAL8 , 0, mpicom, ier)
    call mpi_bcast (latixy_a, size(latixy_a), MPI_REAL8 , 0, mpicom, ier)
#endif

    return
  end subroutine atm_getgrid

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: atm_openfile
!
! !INTERFACE:
  subroutine atm_openfile (kda, kmo, kyr, locfn, itim, atmmin)
!
! !DESCRIPTION: 
! Open atmospheric forcing netCDF file
!
! !USES:
    use clm_varctl, only : offline_atmdir
    use fileutils, only : getfil 
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc' 
    integer, intent(in)  :: kda              !day (1 -> 31)
    integer, intent(in)  :: kmo              !month (1 -> 12)
    integer, intent(in)  :: kyr              !year (0 -> ...)
    character(len=*), intent(inout) :: locfn !history file to open and read
    integer, intent(out) :: itim             !time index used in atmrd
    integer, intent(out) :: atmmin           !temporal resolution of atm data (in minutes)
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer :: i,j,k,n                !do loop indices
    integer :: dimid                  !netCDF dimension id
    integer :: status                 !netCDF error status
    integer :: cyc                    !current cycle
    integer :: nyr                    !year extension to data file name, e.g., 05
    integer :: pyr = 0                !complete years of data since basedate
    integer :: nmo = 0                !number of months past basedate
    integer :: mmo = -1               !number of months past end of data 
    integer :: ier                    !error status
    character(len=256) :: filenam     !full file name, atmdir + ext
    character(len=256) :: locfnlast   !full file name saved from last call
    character(len=  7) :: ext         !month-year extension, e.g., 01-0005
    logical :: lexist                 !true => file exists, used when looking for a file
    integer :: ndaypm(12) =      &    !number of days per month
    	     (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)              
    integer :: minpday = 1440         !minutes per day
!------------------------------------------------------------------------

    if (masterproc) then
#if (!defined COUP_MIT2D)

       ! Build [month]-[year] extension for file name to be read
       ! append extension to path name to get full file name
       ! first save old locfn, then obtain the file if it exists, and
       ! finally check if file to be opened exists (e.g. if are at the
       ! end of the set of files to be cycled, then the next file will
       ! not exists and will have to rewind to beginning of data set
       
       write (ext,'(i4.4,"-",i2.2)') kyr,kmo
       filenam = trim(offline_atmdir) //'/'// ext // '.nc' 
       locfnlast = locfn
       call getfil(filenam, locfn, 1)
       inquire (file = locfn, exist = lexist)

       ! If file exists...
       !   makes sure that if run was restarted within the same month
       !   that the original 'initial' run was started, then nmo will 
       !   not count that month twice
       ! If file doesn't exist...
       !   rewind to beginning of data set (by repeating the process of 
       !   creating extension, appending to path name and looking for 
       !   the file
       
       if (lexist) then
          if (locfn /= locfnlast) nmo = nmo + 1 !months past base date
          pyr = nmo / 12      !years past base date 
       elseif (nmo == 0) then
          write(6,*) 'ATM_OPENFILE error: could not find initial atm datafile' 
          call endrun
       else
          mmo = mmo + 1       !(months-1) past end of data
          cyc = mmo / (pyr*12) !cycles since end of data
          nyr = kyr - pyr * (cyc + 1) !rewind to beginning of dataset
          write (ext,'(i4.4,"-",i2.2)') kyr,kmo
          filenam = trim(offline_atmdir) //'/'// ext // '.nc' 
          call getfil(filenam, locfn, 1)
          inquire (file = locfn, exist = lexist)
          if (.not. lexist) then
             if (nmo < 12) then
                write(6,*) 'ATM_OPENFILE error: You must supply at least a'
                write(6,*) 'year of input data if you wish to'
                write(6,*) 'cycle through it more than once'
             else
                write(6,*) 'ATM_OPENFILE error: Not finding ',locfn
             end if
             call endrun
          end if
       end if
!
       ! Open netCDF data file and get lengths of lat,lon,time dimensions
       ! Do this only at the first timestep of the run or of the month
#else
       locfn = trim(offline_atmdir) 
       write(6,*) 'ATM_GETGRID: atm datafile - ',locfn
#endif

       call wrap_open (locfn, nf_nowrite, ncid)

       call wrap_inq_dimid  (ncid, 'lon', dimid)
       call wrap_inq_dimlen (ncid, dimid, nlon)

       call wrap_inq_dimid  (ncid, 'lat', dimid)
       call wrap_inq_dimlen (ncid, dimid, nlat)

       call wrap_inq_dimid  (ncid, 'time', dimid)
       call wrap_inq_dimlen (ncid, dimid, ntim)

       if (ntim == 0) then
          write (6,*) 'ATM_OPENFILE error: zero input time slices'
          call endrun
       end if
       if (nlon /= atmlon) then
          write (6,*) 'ATM_OPENFILE error: nlon = ',nlon, &
               ' in data file not equal to atmlon = ',atmlon,' first read in'
          call endrun
       end if
       if (nlat /= atmlat) then
          write (6,*) 'ATMRD error: nlat = ',nlat, &
               ' in data file not equal to atmlat = ',atmlat,' first read in'
          call endrun
       end if

       ! Get variable names 

       status = nf_inq_nvars(ncid, nvar)
       if (status /= nf_noerr) then
          write (6,*) ' ATM_OPENFILE netcdf error = ',nf_strerror(status)
          call endrun
       end if
       do i = 1, nvar
          call wrap_inq_varname (ncid, i, varnam(i))
       end do

    endif     !end of if-masterproc block

    ! Calculate temporal resolution of the dataset in minutes and
    ! reset time index

#if (defined SPMD) 
    call mpi_bcast (ntim, 1, MPI_INTEGER, 0, mpicom, ier)
#endif
    atmmin = ndaypm(kmo) * minpday / ntim
    itim = 0                 
    return
  end subroutine atm_openfile

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: atm_readdata
!
! !INTERFACE:
  subroutine atm_readdata (fname, kmo, itim)
!
! !DESCRIPTION: 
! read atmospheric forcing single level fields from netCDF file
! this is done at every subsequent call to atmrd
! this includes a second call at the 1st tstep of the run or of the month
! Use nf_get_vara_real to read data for variable referenced by
! variable id = [i]
! nf_get_vara_real (ncid, i, beg, len, varval)
! integer beg - a vector of integers specifying the index in the
!               variable where the first data value is read.
!               must be dimensioned same as the variable's dimension
!               and a starting value must be given for each dimension
! integer len - a vector of integers specifying the number of data
!               values, for each of the variable's dimensions, to read.
!               must be dimensioned same as the variable's dimension
!               and a length must be given for each dimension
!
! !USES:
    use clm_varcon  , only : sb
    use fileutils   , only : getfil
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    character(len=*), intent(in) :: fname           !history file to open and read
    integer, intent(in)  :: kmo, itim               !current month and time index
!
! !REVISION HISTORY:
! Created by Sam Levis
!
!EOP
!
! LOCAL VARIABLES:
    integer :: i,j,k,n                  !do loop indices
    integer :: ier                      !error status
    integer :: varid                    !netCDF variable id
    integer :: status                   !netCDF error status
    integer :: beg3d(3)                 !netCDF 3-d start index (where to read first value)
    integer :: len3d(3)                 !netCDF 3-d count index (number of values to read)
!
! atm input field names 
!
    character(len=8) :: fldlst(14)      !name of possible atm fields in input file
    data fldlst( 1) /'TBOT    '/
    data fldlst( 2) /'WIND    '/
    data fldlst( 3) /'QBOT    '/ 
    data fldlst( 4) /'Tdew    '/ 
    data fldlst( 5) /'RH      '/ 
    data fldlst( 6) /'ZBOT    '/
    data fldlst( 7) /'PSRF    '/
    data fldlst( 8) /'FSDS    '/
    data fldlst( 9) /'FSDSdir '/
    data fldlst(10) /'FSDSdif '/
    data fldlst(11) /'FLDS    '/
    data fldlst(12) /'PRECTmms'/
    data fldlst(13) /'PRECCmms'/
    data fldlst(14) /'PRECLmms'/

    real(r8) ea                    !atmospheric emissivity

    logical atmread_err

    ! use polynomials to calculate saturation vapor pressure and derivative with
    ! respect to temperature: over water when t > 0 c and over ice when t <= 0 c
    ! required to convert relative humidity to specific humidity

    real(r8) esatw                 !saturation vapor pressure over water (Pa)
    real(r8) esati                 !saturation vapor pressure over ice (Pa)
    real(r8) e                     !vapor pressure (Pa)
    real(r8) qsat                  !saturation specific humidity (kg/kg)
    real(r8) a0,a1,a2,a3,a4,a5,a6  !coefficients for esat over water
    real(r8) b0,b1,b2,b3,b4,b5,b6  !coefficients for esat over ice
    real(r8) tdc, t                !Kelvins to Celcius function and its input

    parameter (a0=6.107799961    , a1=4.436518521e-01, &
               a2=1.428945805e-02, a3=2.650648471e-04, &
               a4=3.031240396e-06, a5=2.034080948e-08, &
               a6=6.136820929e-11)

    parameter (b0=6.109177956    , b1=5.034698970e-01, &
               b2=1.886013408e-02, b3=4.176223716e-04, &
               b4=5.824720280e-06, b5=4.838803174e-08, &
               b6=1.838826904e-10)
!
! function declarations
!
    tdc(t) = min( 50., max(-50.,(t-SHR_CONST_TKFRZ)) )
    esatw(t) = 100.*(a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*a6))))))
    esati(t) = 100.*(b0+t*(b1+t*(b2+t*(b3+t*(b4+t*(b5+t*b6))))))
!------------------------------------------------------------------------


    ! Read single level fields

    if (masterproc) then

       ! initialize fields to the flag value

       x(:,:,:) = -1.         

       ! read input data single-level fields

       beg3d(1) = 1     ;  len3d(1) = atmlon
       beg3d(2) = 1     ;  len3d(2) = atmlat
       beg3d(3) = itim  ;  len3d(3) = 1
       do k = 1, 14
          do n = 1, nvar
             if (varnam(n) == fldlst(k)) then
!          write (6,*) '---------------------------------------'
!          write (6,*) 'ATMRD: getting data for ', varnam(n), 'with ID: ',n
!          write (6,*) '---------------------------------------'
!          write (6,*)
                call wrap_get_vara_realx (ncid,n,beg3d,len3d,x(1,1,k))
             end if
          end do              !end loop of fields in input file
       end do                 !end loop of fields expected in input file

       ! Close file at the end of the month
       ! NOTE: as written will not close file if run ends mid-month

       if (itim == ntim) then
          call wrap_close (ncid)
          write (6,*) '---------------------------------------'
          write (6,*) 'ATMRD: closing data for ',trim(fname)
          write (6,*) '---------------------------------------'
          write (6,*)
       end if

    endif	!end of if-masterproc block

#if (defined SPMD) 
    call mpi_bcast (x, size(x), MPI_REAL8, 0, mpicom, ier)
#endif

    ! ----------------------------------------------------------------------
    ! Determine 2d atmospheric fields
    ! Follow order in fldlst(14) to determine what was read and what was not 
    ! ----------------------------------------------------------------------

    ! Loop over atmospheric longitudes and latitudes

    atmread_err = .false.
!$OMP PARALLEL DO PRIVATE (i,j,e,ea,qsat)
    do j = 1, atmlat
       do i = 1, atmlon

          ! FORC_TXY
          if (nint(x(i,j,1)) == -1) then
             write(6,*)'ATM error: TBOT has not been read by atmrd'
             atmread_err = .true.
          else if (x(i,j,1) < 50. ) then
             write(6,*) 'ATM error: TBOT appears to be in deg C'
             write(6,*) 'Converting to Kelvins now'
             write (6,*) '---------------------------------------'
             write (6,*) 'ATMRD: ',fldlst(1), ' = ',x(i,j,1)
             write (6,*) '---------------------------------------'
             write (6,*)
             forc_txy_a(i,j) = x(i,j,1) + SHR_CONST_TKFRZ
          else
             forc_txy_a(i,j) = x(i,j,1)
          end if

          ! FORC_UXY, FORC_VXY
          if (nint(x(i,j,2)) == -1) then
             write(6,*)'ATM error: WIND has not been read by atmrd'
             atmread_err = .true.
          else
             forc_uxy_a(i,j) = x(i,j,2) / sqrt(2.)
             forc_vxy_a(i,j) = x(i,j,2) / sqrt(2.)
          end if

          ! FORC_PSRFXY, FORC_PBOTXY
          if (nint(x(i,j,7)) == -1) then
             forc_psrfxy_a(i,j) = SHR_CONST_PSTD
          else
             forc_psrfxy_a(i,j) = x(i,j,7)
          end if
          forc_pbotxy_a(i,j)  = forc_psrfxy_a(i,j)

          !FORC_QXY
          if (nint(x(i,j,3)) == -1) then
             if (nint(x(i,j,4)) == -1) then
                if (nint(x(i,j,5)) == -1) then
                   write(6,*)'ATM error: Humidity has not been'
                   write(6,*)'read by atmrd'
                   atmread_err = .true.
                else          !using RH as %
                   if (forc_txy_a(i,j) > SHR_CONST_TKFRZ) then
                      e = x(i,j,5)/100. * esatw(tdc(forc_txy_a(i,j)))
                   else
                      e = x(i,j,5)/100. * esati(tdc(forc_txy_a(i,j)))
                   end if
                end if
                forc_qxy_a(i,j) = 0.622*e / (forc_pbotxy_a(i,j) - 0.378*e)
             else             !using Tdew
                if (x(i,j,4) < 50.) then
                   write(6,*)'ATM warning: Tdew appears to be in'
                   write(6,*)'deg C, so converting to Kelvin'
                   x(i,j,4) = x(i,j,4) + SHR_CONST_TKFRZ
                end if
                if (x(i,j,4) > forc_txy_a(i,j)) then
                   write(6,*)'ATM warning: Dewpt temp > temp!'
                end if
                if (x(i,j,4) > SHR_CONST_TKFRZ) then
                   e = esatw(tdc(x(i,j,4)))
                else
                   e = esati(tdc(x(i,j,4)))
                end if
                forc_qxy_a(i,j) = 0.622*e / (forc_pbotxy_a(i,j) - 0.378*e)
             end if
          else                !using QBOT in kg/kg
             if (forc_txy_a(i,j) > SHR_CONST_TKFRZ) then
                e = esatw(tdc(forc_txy_a(i,j)))
             else
                e = esati(tdc(forc_txy_a(i,j)))
             end if
             qsat = 0.622*e / (forc_pbotxy_a(i,j) - 0.378*e)
             if (qsat < x(i,j,3)) then
                forc_qxy_a(i,j) = qsat
!                 write(6,*)'ATM warning: qsat < q!'
             else
                forc_qxy_a(i,j) = x(i,j,3)
             end if
          end if

          ! ZGCMXY
          if (nint(x(i,j,6)) == -1) then
             zgcmxy_a(i,j) = 30.
          else
             zgcmxy_a(i,j) = x(i,j,6)
          end if

          ! FORC_SOLSXY, FORC_SOLLXY, FORC_SOLSDXY, FORC_SOLLDXY

          if (nint(x(i,j,9))==-1.or.nint(x(i,j,10))==-1) then
             if (nint(x(i,j,8)) /= -1) then
                forc_solsxy_a(i,j)  = 0.7 * (0.5 * x(i,j,8))
                forc_sollxy_a(i,j)  = forc_solsxy_a(i,j)
                forc_solsdxy_a(i,j) = 0.3 * (0.5 * x(i,j,8))
                forc_solldxy_a(i,j) = forc_solsdxy_a(i,j)
             else
                write(6,*)'ATM error: neither FSDSdir/dif nor'
                write(6,*)'       FSDS have been read in by atmrd'
                atmread_err = .true.
             end if
          else
             forc_solsxy_a(i,j)  = 0.5 * x(i,j,9)
             forc_sollxy_a(i,j)  = forc_solsxy_a(i,j)
             forc_solsdxy_a(i,j) = 0.5 * x(i,j,10)
             forc_solldxy_a(i,j) = forc_solsdxy_a(i,j)
          end if

          ! PRCXY, PRLXY

          if (nint(x(i,j,13))==-1.or.nint(x(i,j,14))==-1) then
             if (nint(x(i,j,12)).ne.-1) then
                prcxy_a(i,j) = 0.1 * x(i,j,12)
                prlxy_a(i,j) = 0.9 * x(i,j,12)
             else
                write(6,*)'ATM error: neither PRECC/L nor PRECT'
                write(6,*)'           have been read in by atmrd'
                atmread_err = .true.
             end if
          else
             prcxy_a(i,j) = x(i,j,13)
             prlxy_a(i,j) = x(i,j,14)
          end if

          ! FLWDSXY

          if (nint(x(i,j,11)) == -1) then
             e = forc_psrfxy_a(i,j) * forc_qxy_a(i,j) / (0.622 + 0.378 * forc_qxy_a(i,j))
             ea = 0.70 + 5.95e-05 * 0.01*e * exp(1500.0/forc_txy_a(i,j))             
             flwdsxy_a(i,j) = ea * sb * forc_txy_a(i,j)**4
          else
             flwdsxy_a(i,j) = x(i,j,11) 
          end if

       end do                 !end loop of latitudes
    end do                    !end loop of longitudes
!$OMP END PARALLEL DO

    if (atmread_err) then
       write(6,*) 'atm_readdata: error reading atm data'
       call endrun
    end if

    return
  end subroutine atm_readdata

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: interpa2si
!
! !INTERFACE:
  subroutine interpa2si ()
!
! !DESCRIPTION: 
! Initialize variables for atm->land model surface interp
!
! !USES:
    use clm_varsur, only : numlon, longxy, latixy, lsmedge, lonw, lats, area
    use areaMod  
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Gordon Bonan
!
!EOP
!
! LOCAL VARIABLES:
    integer i,j,k                        !indices
    real(r8), allocatable :: lon_a(:,:)  !atm grid longitude cell edges
    real(r8), allocatable :: lat_a(:)    !atm grid latitude cell edges
    real(r8), allocatable :: area_a(:,:) !atm grid grid cell areas
    real(r8), allocatable :: mask_a(:,:) !dummy field: atm grid mask
    real(r8), allocatable :: mask_s(:,:) !dummy field: land model grid mask
!------------------------------------------------------------------------

    ! Dynamically allocate memory

    allocate (lon_a(atmlon+1,atmlat))
    allocate (lat_a(atmlat+1))
    allocate (area_a(atmlon,atmlat))
    allocate (mask_a(atmlon,atmlat))
    allocate (mask_s(lsmlon,lsmlat))

    if ( masterproc )then
       write (6,*) 'Attempting to initialize atm->land model grid interpolation .....'
       write (6,*) 'Initializing atm -> srf interpolation .....'
    end if

    ! --------------------------------------------------------------------
    ! Map from atmosphere grid to surface grid
    ! --------------------------------------------------------------------

    ! determine numlon for atmosphere grid 

    numlon_a(:) = 0
    do j = 1, atmlat
       do i = 1, atmlon
          if (longxy_a(i,j) /= 1.e36) numlon_a(j) = numlon_a(j) + 1
       end do
    end do

    ! [mask_a] = 1 means all grid cells on atm grid, regardless of whether
    ! land or ocean, will contribute to surface grid.

    do j = 1, atmlat
       do i = 1, numlon_a(j)
          mask_a(i,j) = 1.
       end do
    end do

    ! [mask_s] = 1 means all the surface grid is land. Used as dummy
    ! variable so code will not abort with false, non-valid error check

    do j = 1, lsmlat
       do i = 1, numlon(j)
          mask_s(i,j) = 1.
       end do
    end do

    ! For each surface grid cell: get lat [jovr_a2s] and lon [iovr_a2s] indices 
    ! and weights [wovr_a2s] of overlapping atm grid cells 

    call celledge (atmlat    , atmlon    , numlon_a  , longxy_a  , & 
                   latixy_a  , edge_a(1) , edge_a(2) , edge_a(3) , & 
                   edge_a(4) , lat_a     , lon_a     )

    call cellarea (atmlat    , atmlon    , numlon_a  , lat_a     , &
                   lon_a     , edge_a(1) , edge_a(2) , edge_a(3) , &
                   edge_a(4) , area_a    )

    call areaini (atmlon, atmlat, numlon_a, lon_a, lat_a, area_a, mask_a, &
                  lsmlon, lsmlat, numlon  , lonw , lats , area  , mask_s, &
                  mxovr , novr_a2s, iovr_a2s, jovr_a2s, wovr_a2s )

    deallocate (lon_a)
    deallocate (lat_a)
    deallocate (area_a)
    deallocate (mask_a)
    deallocate (mask_s)

    if ( masterproc )then
       write (6,*) 'Successfully made atm -> srf interpolation'
       write (6,*) 'Successfully initialized area-averaging interpolation'
       write (6,*)
    end if

    return
  end subroutine interpa2si

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: interpa2s
!
! !INTERFACE:
  subroutine interpa2s (forc_t_a  , forc_t_s  , zgcm_a  , zgcm_s  , &
                        forc_u_a  , forc_u_s  , forc_v_a  , forc_v_s  , &
                        forc_q_a  , forc_q_s  , prc_a   , prc_s   , &
                        prl_a   , prl_s   , flwds_a , flwds_s , &
                        forc_sols_a  , forc_sols_s  , forc_soll_a  , forc_soll_s  , &
                        forc_solsd_a , forc_solsd_s , forc_solld_a , forc_solld_s , &
                        forc_pbot_a  , forc_pbot_s  , forc_psrf_a  , forc_psrf_s  )
!
! !DESCRIPTION: 
! Area average fields from atmosphere grid to surface grid
!
! !USES:
    use clm_varsur, only : numlon, longxy, latixy, lsmedge
    use areaMod  
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  ::  forc_t_a(atmlon,atmlat)     !atm bottom level temperature (Kelvin)
    real(r8), intent(in)  ::  zgcm_a(atmlon,atmlat)       !atm bottom level height above surface (m)
    real(r8), intent(in)  ::  forc_u_a(atmlon,atmlat)     !atm bottom level zonal wind (m/s)
    real(r8), intent(in)  ::  forc_v_a(atmlon,atmlat)     !atm bottom level meridional wind (m/s)
    real(r8), intent(in)  ::  forc_q_a(atmlon,atmlat)     !atm bottom level specific humidity (kg/kg)
    real(r8), intent(in)  ::  prc_a(atmlon,atmlat)        !convective precipitation rate (mm H2O/s)
    real(r8), intent(in)  ::  prl_a(atmlon,atmlat)        !large-scale precipitation rate (mm H2O/s)
    real(r8), intent(in)  ::  flwds_a(atmlon,atmlat)      !downward longwave rad onto surface (W/m**2)
    real(r8), intent(in)  ::  forc_sols_a(atmlon,atmlat)  !vis direct beam solar rad onto srf (W/m**2)
    real(r8), intent(in)  ::  forc_soll_a(atmlon,atmlat)  !nir direct beam solar rad onto srf (W/m**2)
    real(r8), intent(in)  ::  forc_solsd_a(atmlon,atmlat) !vis diffuse solar rad onto srf (W/m**2)
    real(r8), intent(in)  ::  forc_solld_a(atmlon,atmlat) !nir diffuse solar rad onto srf(W/m**2)
    real(r8), intent(in)  ::  forc_pbot_a(atmlon,atmlat)  !atm bottom level pressure (Pa)
    real(r8), intent(in)  ::  forc_psrf_a(atmlon,atmlat)  !atm surface pressure (Pa)

    real(r8), intent(out) ::  forc_t_s(lsmlon,lsmlat)     !atm bottom level temperature (Kelvin)
    real(r8), intent(out) ::  zgcm_s(lsmlon,lsmlat)       !atm bottom level height above surface (m)
    real(r8), intent(out) ::  forc_u_s(lsmlon,lsmlat)     !atm bottom level zonal wind (m/s)
    real(r8), intent(out) ::  forc_v_s(lsmlon,lsmlat)     !atm bottom level meridional wind (m/s)
    real(r8), intent(out) ::  forc_q_s(lsmlon,lsmlat)     !atm bottom level specific humidity (kg/kg)
    real(r8), intent(out) ::  prc_s(lsmlon,lsmlat)        !convective precipitation rate (mm H2O/s)
    real(r8), intent(out) ::  prl_s(lsmlon,lsmlat)        !large-scale precipitation rate (mm H2O/s)
    real(r8), intent(out) ::  flwds_s(lsmlon,lsmlat)      !downward longwave rad onto surface (W/m**2)
    real(r8), intent(out) ::  forc_sols_s(lsmlon,lsmlat)  !vis direct beam solar rad onto srf (W/m**2)
    real(r8), intent(out) ::  forc_soll_s(lsmlon,lsmlat)  !nir direct beam solar rad onto srf (W/m**2)
    real(r8), intent(out) ::  forc_solsd_s(lsmlon,lsmlat) !vis diffuse solar rad onto srf (W/m**2)
    real(r8), intent(out) ::  forc_solld_s(lsmlon,lsmlat) !nir diffuse solar rad onto srf(W/m**2)
    real(r8), intent(out) ::  forc_pbot_s(lsmlon,lsmlat)  !atm bottom level pressure (Pa)
    real(r8), intent(out) ::  forc_psrf_s(lsmlon,lsmlat)  !atm surface pressure (Pa)
!
! !REVISION HISTORY:
! Created by Gordon Bonan
!
!EOP
!
! LOCAL VARIABLES:
    integer  :: i,j                  !longitude,latitude loop indices
    real(r8) :: forc_u(atmlon,atmlat)  !dummy wind (u)
    real(r8) :: forc_v(atmlon,atmlat)  !dummy wind (v)
    logical  :: initinterp = .false. !interpolation initialization flag
!------------------------------------------------------------------------

    ! Initialize

    if (.not. initinterp) then
       call interpa2si
       initinterp = .true.
    endif

    ! area-average absolute value of winds (i.e., regardless of
    ! direction) since land model cares about magnitude not direction.
    ! then need to adjust resultant stresses for direction of wind. 
    
!$OMP PARALLEL DO PRIVATE (j,i)
    do j = 1, atmlat
       do i = 1, numlon_a(j)
          forc_u(i,j) = abs(forc_u_a(i,j))
          forc_v(i,j) = abs(forc_v_a(i,j))
       end do
    end do
!$OMP END PARALLEL DO

    call areaave (atmlat   , atmlon   , numlon_a , forc_t_a  , &
                  lsmlat   , lsmlon   , numlon   , forc_t_s  , &
                  iovr_a2s , jovr_a2s , wovr_a2s , mxovr   )

    call areaave (atmlat   , atmlon   , numlon_a , zgcm_a  , &
                  lsmlat   , lsmlon   , numlon   , zgcm_s  , &
                  iovr_a2s , jovr_a2s , wovr_a2s , mxovr   )

    call areaave (atmlat   , atmlon   , numlon_a , forc_u    , &
                  lsmlat   , lsmlon   , numlon   , forc_u_s  , &
                  iovr_a2s , jovr_a2s , wovr_a2s , mxovr   )

    call areaave (atmlat   , atmlon   , numlon_a , forc_v    , &
                  lsmlat   , lsmlon   , numlon   , forc_v_s  , &
                  iovr_a2s , jovr_a2s , wovr_a2s , mxovr   )

    call areaave (atmlat   , atmlon   , numlon_a , forc_q_a  , &
                  lsmlat   , lsmlon   , numlon   , forc_q_s  , &
                  iovr_a2s , jovr_a2s , wovr_a2s , mxovr   )

    call areaave (atmlat   , atmlon   , numlon_a , forc_pbot_a  , &
                  lsmlat   , lsmlon   , numlon   , forc_pbot_s  , &
                  iovr_a2s , jovr_a2s , wovr_a2s , mxovr   )

    call areaave (atmlat   , atmlon   , numlon_a , forc_psrf_a  , &
                  lsmlat   , lsmlon   , numlon   , forc_psrf_s  , &
                  iovr_a2s , jovr_a2s , wovr_a2s , mxovr   )

    call areaave (atmlat   , atmlon   , numlon_a , prc_a   , &
                  lsmlat   , lsmlon   , numlon   , prc_s   , &
                  iovr_a2s , jovr_a2s , wovr_a2s , mxovr   )

    call areaave (atmlat   , atmlon   , numlon_a , prl_a   , &
                  lsmlat   , lsmlon   , numlon   , prl_s   , &
                  iovr_a2s , jovr_a2s , wovr_a2s , mxovr   )

    call areaave (atmlat   , atmlon   , numlon_a , flwds_a , &
                  lsmlat   , lsmlon   , numlon   , flwds_s , &
                  iovr_a2s , jovr_a2s , wovr_a2s , mxovr   )

    call areaave (atmlat   , atmlon   , numlon_a , forc_sols_a  , &
                  lsmlat   , lsmlon   , numlon   , forc_sols_s  , &
                  iovr_a2s , jovr_a2s , wovr_a2s , mxovr   )

    call areaave (atmlat   , atmlon   , numlon_a , forc_soll_a  , &
                  lsmlat   , lsmlon   , numlon   , forc_soll_s  , &
                  iovr_a2s , jovr_a2s , wovr_a2s , mxovr   )

    call areaave (atmlat   , atmlon   , numlon_a , forc_solsd_a , &
                  lsmlat   , lsmlon   , numlon   , forc_solsd_s , &
                  iovr_a2s , jovr_a2s , wovr_a2s , mxovr   )

    call areaave (atmlat   , atmlon   , numlon_a , forc_solld_a , &
                  lsmlat   , lsmlon   , numlon   , forc_solld_s , &
                  iovr_a2s , jovr_a2s , wovr_a2s , mxovr   )

    return
  end subroutine interpa2s

!=======================================================================

#endif

end module atmdrvMod

