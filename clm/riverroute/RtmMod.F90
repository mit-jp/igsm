#include <misc.h>
#include <preproc.h>

module RtmMod

#if (defined RTM) 

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: 
! 
! !DESCRIPTION: 
! River Routing Model (U. of Texas River Transport Model) 
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varpar, only : lsmlon, lsmlat, rtmlon, rtmlat 
  use shr_sys_mod, only : shr_sys_flush
!
! !PUBLIC TYPES:
  implicit none
  save
!
! RTM grid info
!
  real(r8), public, allocatable :: latixy_r(:,:)      !rtm latitudes  of grid cells (degrees)       
  real(r8), public, allocatable :: longxy_r(:,:)      !rtm longitudes of grid cells (degrees)       
  real(r8), public :: area_r(rtmlon,rtmlat)           !rtm gridcell area (km^2)          
  integer , public :: mask_r(rtmlon,rtmlat)           !rtm landmask (land=1,ocean=0)
!
! RTM runoff for coupled communication
!
  integer , public, allocatable :: ocnrof_iindx(:)    !rtm longitude index of ocean runoff point
  integer , public, allocatable :: ocnrof_jindx(:)    !rtm latitude index of ocean runoff point
  real(r8), public, allocatable :: ocnrof_vec(:)      !rtm runoff vector (1/2 deg grid, kg/m^2/s)
!
! RTM history file variables
!
  real(r8), public, pointer :: qchan2(:)              !global river (channel) flow (m**3 H2O /s)
  real(r8), public, pointer :: qchocn2(:)             !global river (channel) flow into ocean (m**3/s)
!
! time averaging for rtm calculation
!
  real(r8), public, pointer, dimension(:) :: totrunin_ave   !time averaged vector of input fluxes
  real(r8), public, pointer, dimension(:) :: prec_ave       !time averaged vector of precipitation
  real(r8), public, pointer, dimension(:) :: evap_ave       !time averaged vector of evaporation
  real(r8), public :: delt_rtm                              !rtm time step
  integer , public :: ncount_rtm                            !number of time samples to average over
!
! Fluxes
!
  real(r8), public  :: volr(rtmlon,rtmlat)            !water volume in cell (m^3)
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: Rtmgridini   ! Initialize RTM grid and land mask
  public :: Rtmlandini   ! Initialize RTM-land interpolation weights
  public :: Rtmfluxini   ! Initialize RTM fluxout 
  public :: Rtmriverflux ! Interface with RTM river routing model
  public :: restart_rtm  ! Read/write RTM restart data
!
! !PRIVATE MEMBER FUNCTIONS:
!  private :: RTM         ! River routing model (based on U. Texas code)
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!EOP
!
! !PRIVATE TYPES:
!
! RTM grid info
!
  private
  integer  :: numlon_r(rtmlat)               !number of lon points at each lat
  real(r8), dimension(4) :: rtmedge = (/ 90., 180., -90., -180. /)  !N,E,S,W edges of rtm grid
!
! land model to RTM mapping. for each rtm grid cell:
!
  integer  :: mxovr_s2r                      !max number of overlapping cells
  integer  :: novr_s2r(rtmlon,rtmlat)        !number    of overlapping lsm cells
  integer , allocatable :: iovr_s2r(:,:,:)   !lon index of overlapping lsm cells
  integer , allocatable :: jovr_s2r(:,:,:)   !lat index of overlapping lsm cells
  real(r8), allocatable :: wovr_s2r(:,:,:)   !weight    of overlapping lsm cells
!
! inputs for rtm calculation
!
  real(r8), pointer, dimension(:) :: totrunin      !input runoff
  real(r8), pointer, dimension(:) :: prec          !input precipitation
  real(r8), pointer, dimension(:) :: evap          !input evaporation
  real(r8), pointer, dimension(:) :: totrunin_glob !temporary for MPI gather
  real(r8), pointer, dimension(:) :: prec_glob     !temporary for MPI gather
  real(r8), pointer, dimension(:) :: evap_glob     !temporary for MPI gather
!
! fluxes
!
  integer  :: rdirc(0:rtmlon+1,0:rtmlat+1)   !rtm river flow direction (0-8)
  real(r8) :: fluxout(0:rtmlon+1,0:rtmlat+1) !water flux out of cell (m^3/s)
  real(r8) :: ddist(rtmlon,rtmlat)           !downstream distance (m)
  real(r8) :: rivarea(rtmlon,rtmlat)         !cell area (m^2)
!
  real(r8), allocatable :: latsh(:)          !southern edge of cells at rtm grid     
  real(r8), allocatable :: lonwh(:,:)        !western  edge of cells at rtm grid     
!
! inputs to RTM at 1/2 degree resolution
!
  real(r8) :: totrunin_r(rtmlon,rtmlat)      !surface runoff (mm/s)
!
! outputs returned from RTM at 1/2 degree resolution
!
  real(r8) :: flxlnd_r(rtmlon,rtmlat)        !river flux (m**3/s)
  real(r8) :: flxocn_r(rtmlon,rtmlat)        !river flux to the ocean (m**3/s)
  real(r8) :: dvolrdt_r(rtmlon,rtmlat)       !change in storage (mm/s)
  real(r8) :: volrtm(rtmlon,rtmlat)          !change in storage (m**3/s)
  real(r8) :: runrtm(rtmlon,rtmlat)          !input runoff on rtm grid (m**3/s)
!
! RTM water flux into cell
!
  real(r8) :: sfluxin(rtmlon,rtmlat)         !water flux into cell (m3/s)
!
! global averaging
!
  character(len=*),parameter :: F40="('(diag) ',a17,'    date  ', &
  &     '   prec        evap        runoff(lnd)   runoff(rtm) dvoldt(rtm) runoff-ocn(rtm)  (m^3/sec)')"
  character(len=*),parameter :: F41="('(diag) ',a17,'   nstep  ', &
  &     '   prec        evap        runoff(lnd)   runoff(rtm) dvoldt(rtm) runoff-ocn(rtm)  (m^3/sec)')"
  character(len=*),parameter :: F21="('(diag) ',a17,' ----------------------', &
  &     7('----------'))"
  character(len=*),parameter :: F22="('(diag) ',a17,i8,6(d13.4))"
!
  real(r8) :: prec_global            !total precipitation (m^3/sec) 
  real(r8) :: evap_global            !total evaporation (m^3/sec)
  real(r8) :: runlnd_global          !total input runoff on land grid (m^3/sec)
  real(r8) :: runrtm_global          !total input runoff on rtm grid (m^3/sec)
  real(r8) :: ocnrtm_global          !total ocean runoff on rtm grid (m^3/sec)
  real(r8) :: volrtm_global          !total change in storage on rtm (m^3/sec)
  integer  :: ncount_global          !global counter 
  integer  :: yrold                  !old year
  real(r8), parameter :: effvel = 0.35    !effective velocity (m/s)

!Temporaries
#if (defined SPMD)
  real(r8), pointer, dimension(:) :: spmdtemp
#endif 
!----------------------------------------------------------------------- 

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtmgridini
!
! !INTERFACE:
  subroutine Rtmgridini
!
! !DESCRIPTION: 
! Initialize RTM grid and land mask (U. of Texas River Transport Model)
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8	
    use spmdMod, only : masterproc
    use areaMod, only : celledge, cellarea
    use clm_varctl, only : frivinp_rtm
    use clm_varcon, only : re
    use shr_const_mod, only: SHR_CONST_PI
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: ioff(0:8) = (/0,0,1,1,1,0,-1,-1,-1/) !calc dist as in hydra
    integer  :: joff(0:8) = (/0,1,1,0,-1,-1,-1,0,1/) !of grid cell down stream
    integer  :: i,j,k,n                       !loop indices
    integer  :: i2,j2                         !downstream i and j
    real(r8) :: deg2rad                       !pi/180
    real(r8) :: dx                            !lon dist. between grid cells (m)
    real(r8) :: dy                            !lat dist. between grid cells (m)
    real(r8) :: dist(rtmlon,rtmlat)           !dist. of the grid cell down stream (m) 
    real(r8) :: tempg(rtmlon,rtmlat)          !temporary buffer
    integer  :: tempgp(0:rtmlon+1,0:rtmlat+1) !temporary buffer  
!-----------------------------------------------------------------------
  
    if (masterproc) then
       
       ! Useful constants and initial values
       
       write(6,*)'Columns in RTM = ',rtmlon
       write(6,*)'Rows in RTM    = ',rtmlat
       
       allocate(latixy_r(rtmlon,rtmlat))
       allocate(longxy_r(rtmlon,rtmlat))
       allocate(latsh(rtmlat+1))        !southern edge of cells at rtm grid     
       allocate(lonwh(rtmlon+1,rtmlat)) !western  edge of cells at rtm grid     
       
       deg2rad = SHR_CONST_PI / 180.
       volr(:,:) = 0.
       fluxout(:,:) = 0.
       flxocn_r(:,:) = 0.
       flxlnd_r(:,:) = 0.

       ! Open and read input data (river direction file)
       ! rtm operates from south to north and from the dateline
       
       open (1,file=frivinp_rtm)
       write(6,*)'opened river direction data'
       do j = 1,rtmlat
          numlon_r(j) = 0
          do i = 1,rtmlon
             read(1,*) latixy_r(i,j),longxy_r(i,j),tempg(i,j)
             if (longxy_r(i,j) /= 1.e36) numlon_r(j) = numlon_r(j) + 1
             tempgp(i,j) = nint(tempg(i,j))
          enddo
       enddo
       close(1)
       write(6,*)'closed river direction data'
       write(6,*)
       
       ! Determine RTM celledges, areas and interpolation masks
       
       call celledge (rtmlat    , rtmlon    , numlon_r  , longxy_r  , &
                      latixy_r  , rtmedge(1), rtmedge(2), rtmedge(3), &
                      rtmedge(4), latsh     , lonwh     )

       call cellarea (rtmlat    , rtmlon    , numlon_r  , latsh     , lonwh , &
                      rtmedge(1), rtmedge(2), rtmedge(3), rtmedge(4), area_r) 

       ! Determine rtm mask, downstream distance and area
       
       ! determine rtm ocn/land mask
       
       do i=1,rtmlon
          tempgp(i,0)        = tempgp(mod(i+rtmlon/2-1,rtmlon)+1,1)
          tempgp(i,rtmlat+1) = tempgp(mod(i+rtmlon/2-1,rtmlon)+1,rtmlat)
          if (tempgp(i,0)        /= 0) tempgp(i,0)        = mod(tempgp(i,0)       +4-1,8)+1
          if (tempgp(i,rtmlat+1) /= 0) tempgp(i,rtmlat+1) = mod(tempgp(i,rtmlat+1)+4-1,8)+1
       enddo
       do j=0,rtmlat+1
          tempgp(0,j) =tempgp(rtmlon,j)
          tempgp(rtmlon+1,j)=tempgp(1,j)
       enddo
       
       do j=0,rtmlat+1
          do i=0,rtmlon+1
             rdirc(i,j)=tempgp(i,j)
          enddo
       enddo
     
       do j=1,rtmlat
          do i=1,rtmlon
             if (rdirc(i,j) == 0) then
                mask_r(i,j) = 0
             else
                mask_r(i,j) = 1
             end if
          enddo
       enddo
       
       ! determine downstream distance - instead of reading a distance file 
       ! calculate the downstream distance as in hydra
       
       do j=1,rtmlat
          do i=1,rtmlon
             i2 = i + ioff(tempgp(i,j))
             j2 = j + joff(tempgp(i,j))
             if (i2 == 0) i2 = 2                 !avoids i2 out of bounds in the following
             if (i2 == rtmlon+1) i2 = rtmlon-1   !avoids i2 out of bounds in the following  
             dy = deg2rad * abs(latixy_r(i,j)-latixy_r(i2,j2)) * re*1000.
             dx = deg2rad * abs(longxy_r(i,j)-longxy_r(i2,j2)) * re*1000. &
                  *0.5*(cos(latixy_r(i,j)*deg2rad)+cos(latixy_r(i2,j2)*deg2rad))
             dist(i,j) = sqrt(dx*dx + dy*dy)
             ddist(i,j) = dist(i,j)
             rivarea(i,j)=1.e6 * area_r(i,j)     !convert into m**2
          enddo
       enddo
       
    endif  ! end of if-masterproc block

    return
  end subroutine Rtmgridini

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtmlandini
!
! !INTERFACE:
  subroutine Rtmlandini
!
! !DESCRIPTION: 
! Initialize RTM-land interpolation weights 
! (U. of Texas River Transport Model)
! and variables related to runoff time averaging
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8	
    use clmpoint
    use spmdMod, only : masterproc
    use areaMod, only : areaini, mkmxovr            
    use clm_varsur, only : numlon, area, lats, lonw, landmask
    use time_manager, only : get_curr_date
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: i,j,k,n                    !loop indices
    integer  :: is,js                      !land model grid indices
    real(r8) :: maskone_s(lsmlon,lsmlat)   !dummy field: see below                 
    real(r8) :: maskone_r(rtmlon,rtmlat)   !dummy field: see below                 
    integer  :: ocnrof_mask(rtmlon,rtmlat) !rtm mask for ocean points with possible nonzero runoff
    integer  :: ocnrof_num                 !number of valid ocean points with possible nonzero runoff
    integer  :: yrnew                      !year (0, ...)
    integer  :: mon                        !month (1, ..., 12) 
    integer  :: day                        !day of month (1, ..., 31)
    integer  :: ncsec                      !seconds of current date
    integer  :: begc,endc,numc             !column 1d indices
!-----------------------------------------------------------------------

    if (masterproc) then

       ! --------------------------------------------------------------------
       ! The following section allows RTM and land model to coexist at different
       ! horizontal resolutions
       ! --------------------------------------------------------------------

       write(6,*)
       write(6,*) 'Initializing area-averaging interpolation for RTM.....'

       ! To find fraction of each land model grid cell that is land based on rtm grid.
       ! For this purpose, want all rtm grid cells to contribute to grid cell 
       ! average on land model grid, i.e., all cells used regardless of whether land 
       ! or ocean. Do this by setting [maskone_s] = 1 

       ! [maskone_s] = 1 means all grid cells on land model grid, regardless of whether
       ! land or ocean, will contribute to rtm grid.

       do j = 1, lsmlat
          do i = 1, numlon(j)
             maskone_s(i,j) = 1.
          end do
       end do

       ! [maskone_r] = 1 means all the rtm grid is land. Used as dummy
       ! variable so code will not abort with false, non-valid error check

       do j = 1, rtmlat
          do i = 1, numlon_r(j)
             maskone_r(i,j) = 1.
          end do
       end do

       ! --------------------------------------------------------------------
       ! Map weights from land model grid to rtm grid
       ! --------------------------------------------------------------------

       write(6,*) 'Initializing land model -> rtm interpolation .....'

       ! For each rtm grid cell: get lat [jovr_s2r] and lon [iovr_s2r] indices 
       ! and weights [wovr_s2r] of overlapping atm grid cells 

       call mkmxovr (lsmlon, lsmlat, numlon  , lonw , lats , &
            rtmlon, rtmlat, numlon_r, lonwh, latsh, &
            mxovr_s2r     , novr_s2r)
       
       allocate(iovr_s2r(rtmlon,rtmlat,mxovr_s2r))
       allocate(jovr_s2r(rtmlon,rtmlat,mxovr_s2r))
       allocate(wovr_s2r(rtmlon,rtmlat,mxovr_s2r))

       call areaini (lsmlon   , lsmlat  , numlon  , lonw , lats , area  , maskone_s,  &
            rtmlon   , rtmlat  , numlon_r, lonwh, latsh, area_r, maskone_r,  &
            mxovr_s2r, novr_s2r, iovr_s2r, jovr_s2r, wovr_s2r)

       write(6,*) 'Successfully made land model -> rtm interpolation'
       write(6,*)

#if (defined COUP_CSM)

       ! --------------------------------------------------------------------
       ! Determine which ocean cells might have runoff values. 
       ! --------------------------------------------------------------------

       ! First loop over all ocean points and determine which are at the 
       ! end of rivers by examining if any neighboring points are land and 
       ! if that land neighbor points into this ocean point. Next loop over all
       ! ocean points and determine which overlap with at least one land cell.

       ocnrof_num = 0
       ocnrof_mask(:,:) = 0
       do j=1,rtmlat
          do i=1,rtmlon
             if (mask_r(i,j) == 0) then
                if (rdirc(i  ,j-1)==1) ocnrof_mask(i,j) = 1
                if (rdirc(i-1,j-1)==2) ocnrof_mask(i,j) = 1
                if (rdirc(i-1,j  )==3) ocnrof_mask(i,j) = 1
                if (rdirc(i-1,j+1)==4) ocnrof_mask(i,j) = 1
                if (rdirc(i  ,j+1)==5) ocnrof_mask(i,j) = 1
                if (rdirc(i+1,j+1)==6) ocnrof_mask(i,j) = 1
                if (rdirc(i+1,j  )==7) ocnrof_mask(i,j) = 1
                if (rdirc(i+1,j-1)==8) ocnrof_mask(i,j) = 1
                if (ocnrof_mask(i,j) == 0) then
                   do n=1,novr_s2r(i,j)
                      is = iovr_s2r(i,j,n)
                      js = jovr_s2r(i,j,n)
                      if (landmask(is,js)==1 .and. wovr_s2r(i,j,n)>0.) then
                         ocnrof_mask(i,j) = 1
                      end if
                   end do
                endif
             endif
             if (ocnrof_mask(i,j) == 1) ocnrof_num = ocnrof_num +1
          enddo
       enddo

       ! allocate ocean runoff vector and indices and determine indices
       ! need to reset ocnrof_num to 0 and do the counting again because need to first
       ! first count to allocate vector and must now count to actually determine indices

       allocate(ocnrof_vec  (ocnrof_num)) 
       allocate(ocnrof_iindx(ocnrof_num)) 
       allocate(ocnrof_jindx(ocnrof_num))

       ocnrof_num = 0
       do j=1,rtmlat
          do i=1,rtmlon
             if (ocnrof_mask(i,j) == 1) then
                ocnrof_num = ocnrof_num + 1
                ocnrof_iindx(ocnrof_num) = i
                ocnrof_jindx(ocnrof_num) = j
                ocnrof_vec(ocnrof_num) = 0.
             endif
          end do
       enddo

#endif

       ! Deallocate memory for rtm grid  - needed to be done here because
       ! rtm grid information had to be sent to coupler between calls to 
       ! Rtmgridini and Rtmlandini

       deallocate(latixy_r)
       deallocate(longxy_r)
       deallocate(latsh)
       deallocate(lonwh) 

    endif ! end of if-masterproc block

    ! Initialize rtm time averaging variables
    ! Upon restart, the following variables will get new values 
    ! from the restart file - these values are only valid for
    ! initial runs

    ncount_rtm    = 0

    ncount_global = 0

    prec_global   = 0.  
    evap_global   = 0.
    runlnd_global = 0.
    runrtm_global = 0.
    volrtm_global = 0.
    ocnrtm_global = 0.

    call get_curr_date(yrold, mon, day, ncsec)

    begc = cols1d%beg
    endc = cols1d%end
    numc = cols1d%num

    allocate (totrunin_ave(numc)); totrunin_ave(:) = 0.
    allocate (prec_ave(numc)); prec_ave(:) = 0.
    allocate (evap_ave(numc)); evap_ave(:) = 0.

    allocate (totrunin(numc)); totrunin(:) = 0.
    allocate (prec(numc)); prec(:) = 0.
    allocate (evap(numc)); evap(:) = 0.

#if (defined SPMD)
    if (masterproc) then
       allocate (totrunin_glob(numc)); totrunin_glob(:) = 0.
       allocate (prec_glob(numc)); prec_glob(:) = 0.
       allocate (evap_glob(numc)); evap_glob(:) = 0.
    endif
#endif

    allocate (qchan2(numc)); qchan2(:) = 0.
    allocate (qchocn2(numc)); qchocn2(:) = 0.
#if (defined SPMD)
    allocate (spmdtemp(numc)); spmdtemp(:) = 0.
#endif

    do k = begc,endc
       clmpointer(ic_cwf_qchan2)%val(k)%rp = 0.
       clmpointer(ic_cwf_qchocn2)%val(k)%rp = 0.
    end do
    
    return
  end subroutine Rtmlandini

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtmfluxini()
!
! !INTERFACE:
  subroutine Rtmfluxini()
!
! !DESCRIPTION: 
! Initialize RTM fluxout for case of initial run when initial data is 
! read in. For restart run, RTM fluxout is read from restart dataset.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8	
    use time_manager, only : get_step_size
    use clm_varctl, only : rtm_nsteps
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j	 !indices
    integer :: delt_rtm    !delt for rtm
!-----------------------------------------------------------------------

    delt_rtm = rtm_nsteps*get_step_size()   
    do j=1,rtmlat
       do i=1,rtmlon
          if (mask_r(i,j)==1) then
             fluxout(i,j) = volr(i,j) * effvel/ddist(i,j)
             fluxout(i,j) = min(fluxout(i,j), volr(i,j) / delt_rtm)
          else
             fluxout(i,j) = 0.  
          endif
       enddo
    enddo

    return
  end subroutine Rtmfluxini

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtmriverflux
!
! !INTERFACE:
  subroutine Rtmriverflux ()
!
! !DESCRIPTION: 
! Interface with RTM river routing model
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8	
    use clmtype
    use clmpoint
    use clm_varpar, only : lsmlon, lsmlat
    use clm_varsur, only : numlon, area, landfrac
    use clm_varctl, only : rtm_nsteps
#if (defined SPMD)
    use spmdMod, only : masterproc, scatter_data_from_master, gather_data_to_master, iam
#else
    use spmdMod, only : masterproc
#endif
    use time_manager, only : get_step_size, get_curr_date, get_nstep
    use mapxy, only : vec2xy, xy2vec
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!EOP
!
! !LOCAL VARIABLES:
!
! misc variables
!
    integer  :: gi,li,ci,pi                   !indices 
    integer  :: io,jo,ir,jr,is,js             !mapping indices
    integer  :: k,n,i,j	                      !indices
    integer  :: index                         !index into 1d column array 
    real(r8) :: wt                            !weight
    type(column_type)       , pointer :: c    !local pointer to derived subtype
    type(column_wflux_type) , pointer :: cwf  !local pointer to derived subtype
    type(atm2lnd_flux_type) , pointer :: a2lf !local pointer to derived subtype
!
! inputs to RTM at land model resolution
!
    real(r8) :: totruninxy(lsmlon,lsmlat)     !surface runoff (mm H2O /s)
    real(r8) :: precxy(lsmlon,lsmlat)         !precipitation (mm H2O /s)
    real(r8) :: evapxy(lsmlon,lsmlat)         !evaporation (mm H2O /s)
!
! outputs returned from RTM converted to land model resolution 
!
    real(r8) :: flxout_s(lsmlon,lsmlat)       !river flow (m**3)
    real(r8) :: flxocn_s(lsmlon,lsmlat)       !flow into ocean (m**3)
!
! global balance
!
    integer  :: yrnew         !year (0, ...)
    integer  :: mon           !month (1, ..., 12) 
    integer  :: day           !day of month (1, ..., 31)
    integer  :: ncsec         !seconds of current date  
    integer  :: ncdate        !current date   
    real(r8) :: prec_sum      !total precipitation (m^3/sec) 
    real(r8) :: evap_sum      !total evaporation (m^3/sec)
    real(r8) :: runlnd_sum    !total input runoff on land grid (m^3/sec)
    real(r8) :: runrtm_sum    !total input runoff on rtm grid (m^3/sec)
    real(r8) :: ocnrtm_sum    !total ocean runoff on rtm grid (m^3/sec)
    real(r8) :: volrtm_sum    !total change in storage on rtm (m^3/sec)
    integer  :: nstep         !time step index

    integer :: begc,endc,numc                 !column 1d indices
    character(len=16) :: namec                !column 1d name 
    integer :: ier                            !error status
!-----------------------------------------------------------------------

#if (defined TIMING_BARRIERS)
    call t_startf ('sync_clmrtm')
    call mpibarrier (mpicom)
    call t_stopf ('sync_clmrtm')
#endif

    ! --------------------------------------------------------------------
    ! RTM inputs 
    ! --------------------------------------------------------------------

    ! Make gridded representation of runoff from subgrid data
    ! total surface runoff = surface runoff on soils 
    ! + runoff on glaciers, wetlands, and lakes (P-E) 

!$OMP PARALLEL DO PRIVATE (gi,a2lf,li,ci,c,cwf,index)
    do gi = 1, clm%mps%ngridcells
       a2lf => clm%g(gi)%a2lf
       do li = 1, clm%g(gi)%gps%nlandunits
          do ci = 1, clm%g(gi)%l(li)%lps%ncolumns
             c => clm%g(gi)%l(li)%c(ci)
             cwf => c%cwf
             index = c%cps%index1d
             totrunin(index) = cwf%qflx_surf + cwf%qflx_qrgwl + cwf%qflx_drain
             prec(index) = a2lf%forc_rain + a2lf%forc_snow
             evap(index) = cwf%pwf_a%qflx_evap_tot
          end do
       end do
    end do

    ! --------------------------------------------------------------------
    ! Average fluxes for RTM calculation if appropriate
    ! --------------------------------------------------------------------

    namec = cols1d%name
    begc = cols1d%beg
    endc = cols1d%end
    numc = cols1d%num

    if (rtm_nsteps <= 1) then

       ! RTM averaging is not done

#if (defined SPMD)
       call gather_data_to_master(totrunin, totrunin_glob, clmlevel=namec)
       call gather_data_to_master(prec, prec_glob, clmlevel=namec)
       call gather_data_to_master(evap, evap_glob, clmlevel=namec)
#else
       totrunin_glob => totrunin
       prec_glob => prec
       evap_glob => evap
#endif           
       if (masterproc) then
          call vec2xy (namec, totrunin_glob, 0._r8, totruninxy)
          call vec2xy (namec, prec_glob, 0._r8, precxy)
          call vec2xy (namec, evap_glob, 0._r8, evapxy)
          delt_rtm = get_step_size()
       endif

    else

       ! RTM averaging is done only done by master processor - however
       ! all SPMD processe will continue below

       do k = begc,endc
          totrunin_ave(k) = totrunin_ave(k) + totrunin(k)
          prec_ave(k) = prec_ave(k) + prec(k)
          evap_ave(k) = evap_ave(k) + evap(k)
       end do
       if (masterproc) then
          ncount_rtm = ncount_rtm + 1     
       endif
       nstep = get_nstep()

       if ((mod(nstep,rtm_nsteps)==0) .and. (nstep>1)) then
#if (defined SPMD)
          call gather_data_to_master(totrunin_ave, totrunin_glob, clmlevel=namec)
          call gather_data_to_master(prec_ave, prec_glob, clmlevel=namec)
          call gather_data_to_master(evap_ave, evap_glob, clmlevel=namec)
#else
          totrunin_glob => totrunin_ave
          prec_glob => prec_ave
          evap_glob => evap_ave
#endif           
          if (masterproc) then
             do k = 1,numc
                totrunin_glob(k) = totrunin_glob(k)/ncount_rtm
                prec_glob(k) = prec_glob(k)/ncount_rtm
                evap_glob(k) = evap_glob(k)/ncount_rtm
             end do
             call vec2xy(namec, totrunin_glob, 0._r8, totruninxy)
             call vec2xy(namec, prec_glob, 0._r8, precxy)
             call vec2xy(namec, evap_glob, 0._r8, evapxy)
             delt_rtm = ncount_rtm*get_step_size()   !compute delt for rtm
             ncount_rtm = 0                          !reset counter to 0
          endif
          do k = begc,endc
             totrunin_ave(k) = 0.                    !reset averager
             prec_ave(k) = 0.
             evap_ave(k) = 0.
          end do
       else
!$OMP PARALLEL DO PRIVATE (k)
          do k = begc,endc
             clmpointer(ic_cwf_qchan2)%val(k)%rp = qchan2(k)
             clmpointer(ic_cwf_qchocn2)%val(k)%rp = qchocn2(k)
          end do
          RETURN
       endif

    endif

    ! --------------------------------------------------------------------
    ! RTM runoff - only master processor performs runoff computation
    ! --------------------------------------------------------------------

    if (masterproc) then

       ! Map from land model grid to RTM grid (intepolate to 1/2 degree resolution)

!$OMP PARALLEL DO PRIVATE (jr,ir,n,is,js,wt)
       do jr = 1,rtmlat
          do ir = 1,numlon_r(jr)
             totrunin_r(ir,jr) = 0.
             do n = 1, novr_s2r(ir,jr)
                if (wovr_s2r(ir,jr,n) > 0.) then
                   is = iovr_s2r(ir,jr,n)
                   js = jovr_s2r(ir,jr,n)
                   wt = wovr_s2r(ir,jr,n)
                   totrunin_r(ir,jr) = totrunin_r(ir,jr) + wt*totruninxy(is,js)*landfrac(is,js)
                end if
             end do
          end do
       end do

       ! Determine fluxes on 1/2 degree grid

       call Rtm 

#if (defined COUP_CSM)

       ! Determine ocean runoff vector to send to coupler 

       do n = 1,size(ocnrof_vec)
          i = ocnrof_iindx(n)
          j = ocnrof_jindx(n)
          ocnrof_vec(n) = flxocn_r(i,j)/(area_r(i,j)*1000.) ! units of kg/m^2/s
       end do

#endif

       ! Determine ocean runoff and total runoff on land model grid and
       ! compute global input runoff on rtm grid, global ocean runoff on
       ! rtm grid and global change in storage on rtm grid

       do js = 1, lsmlat
          do is =1, numlon(js)
             flxout_s(is,js)  = 0.
             flxocn_s(is,js)  = 0.
          end do
       end do
       do jr = 1,rtmlat
          do ir = 1,numlon_r(jr)
             do n = 1, novr_s2r(ir,jr)
                if (wovr_s2r(ir,jr,n) > 0.) then
                   is = iovr_s2r(ir,jr,n)
                   js = jovr_s2r(ir,jr,n)
                   wt = wovr_s2r(ir,jr,n)
                   flxocn_s(is,js)  = flxocn_s(is,js) + wt*flxocn_r(ir,jr)
                   flxout_s(is,js)  = max(flxout_s(is,js), wt*flxlnd_r(ir,jr))
                end if
             end do
             runrtm(ir,jr) = totrunin_r(ir,jr)*1000.*area_r(ir,jr)
             volrtm(ir,jr) = dvolrdt_r(ir,jr)*1000.*area_r(ir,jr)
          end do
       end do

       ! Determine global quantities and increment global counter

!$OMP PARALLEL DO PRIVATE (js,is)
       do js = 1,lsmlat
          do is = 1,numlon(js)
             precxy(is,js) = precxy(is,js)*area(is,js)*1000.*landfrac(is,js)
             evapxy(is,js) = evapxy(is,js)*area(is,js)*1000.*landfrac(is,js)
             totruninxy(is,js) = totruninxy(is,js)*area(is,js)*1000.*landfrac(is,js)
          end do
       end do

       prec_sum   = sum(precxy)
       evap_sum   = sum(evapxy)
       runlnd_sum = sum(totruninxy)
       runrtm_sum = sum(runrtm)
       volrtm_sum = sum(volrtm)
       ocnrtm_sum = sum(flxocn_r)

       prec_global   = prec_global   + prec_sum
       evap_global   = evap_global   + evap_sum
       runlnd_global = runlnd_global + runlnd_sum
       runrtm_global = runrtm_global + runrtm_sum
       volrtm_global = volrtm_global + volrtm_sum
       ocnrtm_global = ocnrtm_global + ocnrtm_sum

       ncount_global = ncount_global + 1

       ! Print out diagnostics if appropriate

!       write(6,*)
!       write(6,F41)'water inst   '
!       write(6,F21)'water inst   '
!       write(6,F22)'water inst   ',get_nstep(), prec_sum, evap_sum, &
!            runlnd_sum, runrtm_sum, volrtm_sum, ocnrtm_sum
!       write(6,*)

       call get_curr_date(yrnew, mon, day, ncsec)
       ncdate = yrnew*10000 + mon*100 + day
       if (yrnew /= yrold) then
          prec_global   = prec_global/ncount_global
          evap_global   = evap_global/ncount_global
          runlnd_global = runlnd_global/ncount_global
          runrtm_global = runrtm_global/ncount_global
          volrtm_global = volrtm_global/ncount_global
          ocnrtm_global = ocnrtm_global/ncount_global
          ncount_global = 0
          write(6,*)
          write(6,F40)'water tavg   '
          write(6,F21)'water tavg   '
          write(6,F22)'water tavg   ',ncdate, prec_global, evap_global,&
               runlnd_global, runrtm_global, volrtm_global, ocnrtm_global
          write(6,*)
       endif
       yrold = yrnew

       ! Convert gridded output to vector form 

       call xy2vec (namec, flxocn_s, qchocn2)
       call xy2vec (namec, flxout_s, qchan2)

    endif  !end of if-masterproc block

#if (defined SPMD)
    if (masterproc) spmdtemp(:) = qchocn2(:)
    call scatter_data_from_master(qchocn2, spmdtemp, clmlevel=namec)
    if (masterproc) spmdtemp(:) = qchocn2(:)
    call scatter_data_from_master(qchan2, qchan2, clmlevel=namec)
#endif
    do k = begc,endc
       clmpointer(ic_cwf_qchan2)%val(k)%rp = qchan2(k)
       clmpointer(ic_cwf_qchocn2)%val(k)%rp = qchocn2(k)
    end do

  end subroutine Rtmriverflux

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtm
!
! !INTERFACE:
  subroutine Rtm 
!
! !DESCRIPTION: 
! River routing model (based on U. Texas code)
! input is totrunin_r
! input/output is fluxout, volr
! outputs are dvolrdt_r, flxocn_r, flxlnd_r
!
! !USES:
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine Rtmriverflux in this module
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: i, j                        !loop indices
    real(r8) :: buffern(rtmlon)             !temp buffer
    real(r8) :: buffers(rtmlon)             !temp buffer
    real(r8) :: dvolrdt                     !change in storage (m3/s)
    real(r8) :: sumdvolr(rtmlat)            !global sum (m3/s)
    real(r8) :: sumrunof(rtmlat)            !global sum (m3/s)
    real(r8) :: sumdvolr_tot                !global sum (m3/s)
    real(r8) :: sumrunof_tot                !global sum (m3/s)
!-----------------------------------------------------------------------

    ! Determine fluxout at extended points and at southern and northern outer lats

    fluxout(rtmlon+1,rtmlat+1) = fluxout(1,1)
    fluxout(rtmlon+1,0)        = fluxout(1,rtmlat)
    fluxout(0,0)               = fluxout(rtmlon,rtmlat)
    fluxout(0,rtmlat+1)        = fluxout(rtmlon,1)

    do i=1,rtmlon
       fluxout(i,0)        = fluxout(i,rtmlat)
       fluxout(i,rtmlat+1) = fluxout(i,1)
       buffern(i)          = fluxout(i,rtmlat)
       buffers(i)          = fluxout(i,1)
    enddo
    do j=1,rtmlat
       fluxout(0,j)        = fluxout(rtmlon,j)
       fluxout(rtmlon+1,j) = fluxout(1,j)
    enddo
    do i=0,rtmlon+1
       fluxout(i,0)        = buffern(mod(i+rtmlon/2-1,rtmlon)+1)
       fluxout(i,rtmlat+1) = buffers(mod(i+rtmlon/2-1,rtmlon)+1)
    enddo

    ! Determine cell-to-cell transport - calculate sfluxin

!$OMP PARALLEL DO PRIVATE (I,J)
    do j=1,rtmlat
       do i=1,rtmlon
          sfluxin(i,j) = 0.
          if (rdirc(i  ,j-1)==1) sfluxin(i,j) = sfluxin(i  ,j) + fluxout(i  ,j-1)
          if (rdirc(i-1,j-1)==2) sfluxin(i,j) = sfluxin(i  ,j) + fluxout(i-1,j-1)
          if (rdirc(i-1,j  )==3) sfluxin(i,j) = sfluxin(i  ,j) + fluxout(i-1,j  )
          if (rdirc(i-1,j+1)==4) sfluxin(i,j) = sfluxin(i  ,j) + fluxout(i-1,j+1)
          if (rdirc(i  ,j+1)==5) sfluxin(i,j) = sfluxin(i  ,j) + fluxout(i  ,j+1)
          if (rdirc(i+1,j+1)==6) sfluxin(i,j) = sfluxin(i  ,j) + fluxout(i+1,j+1)
          if (rdirc(i+1,j  )==7) sfluxin(i,j) = sfluxin(i  ,j) + fluxout(i+1,j  )
          if (rdirc(i+1,j-1)==8) sfluxin(i,j) = sfluxin(i  ,j) + fluxout(i+1,j-1)
       enddo
    enddo

    ! Loops above and below must remain separate because fluxout is updated below

    sumdvolr(:) = 0.
    sumrunof(:) = 0.
!$OMP PARALLEL DO PRIVATE (I,J,DVOLRDT)
    do j=1,rtmlat
       do i=1,rtmlon

          ! calculate change in cell storage volume change units for 
          !totrunin from kg/m2s==mm/s -> m3/s

          dvolrdt = sfluxin(i,j) - fluxout(i,j) + 0.001*totrunin_r(i,j)*rivarea(i,j)

          ! calculate flux out of a cell:
          ! land: do not permit change in cell storage volume greater than volume present
          ! make up for the difference with an extraction term (eg from aquifers)
          ! ocean: do not permit negative change in cell storage volume,
          ! because at ocean points cell storage volume equals zero
          ! water balance check (in mm/s), convert runinxy from mm/s to m/s (* 1.e-3) 
          ! and land model area from km**2 to m**2 (* 1.e6)

          if (mask_r(i,j) == 1) then         ! land points
             volr(i,j)    = volr(i,j) + dvolrdt*delt_rtm
             fluxout(i,j) = volr(i,j) * effvel/ddist(i,j)
             fluxout(i,j) = min(fluxout(i,j), volr(i,j) / delt_rtm)
             flxlnd_r(i,j) = fluxout(i,j)
             flxocn_r(i,j) = 0.
          else                               ! ocean points
             flxlnd_r(i,j) = 0.
             flxocn_r(i,j) = dvolrdt
          endif
          sumdvolr(j) = sumdvolr(j) + dvolrdt
          sumrunof(j) = sumrunof(j) + totrunin_r(i,j)*1000.*area_r(i,j)
          dvolrdt_r(i,j) = 1000.*dvolrdt/rivarea(i,j)

       enddo
    enddo

    ! Global water balance calculation and error check

    sumdvolr_tot = 0.
    sumrunof_tot = 0.
    do j = 1,rtmlat
       sumdvolr_tot = sumdvolr_tot + sumdvolr(j)
       sumrunof_tot = sumrunof_tot + sumrunof(j)
    end do
    if (abs((sumdvolr_tot-sumrunof_tot)/sumrunof_tot) > 0.01) then
       write(6,*) 'RTM Error: sumdvolr= ',sumdvolr_tot,&
            ' not equal to sumrunof= ',sumrunof_tot
       call endrun
    end if

    return
  end subroutine Rtm

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restart_rtm
!
! !INTERFACE:
  subroutine restart_rtm (nio, flag)
!
! !DESCRIPTION: 
! Read/write RTM restart data
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8	          
    use clmpoint
    use clmtype, only : cols1d
#if (defined SPMD)
    use spmdMod, only : masterproc, scatter_data_from_master, gather_data_to_master
#else
    use spmdMod, only : masterproc
#endif
    use shr_sys_mod, only : shr_sys_flush
!
! !ARGUMENTS: 
    implicit none
    integer, intent(in) :: nio             !restart unit 
    character(len=*), intent(in) :: flag   !'read' or 'write'
!
! !CALLED FROM:
! subroutine restart in module restFileMod
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: k                             !indices 
    integer :: begc,endc,numc                !column 1d indices
    character (len=16) :: namec              !column 1d name 
!-----------------------------------------------------------------------

    begc = cols1d%beg
    endc = cols1d%end
    numc = cols1d%num
    namec = cols1d%name

    if (flag == 'read') then
       if (masterproc) then
          read (nio) volr
          read (nio) fluxout
          read (nio) ncount_rtm, ncount_global, yrold
          read (nio) prec_global, evap_global, runlnd_global, &
               runrtm_global, volrtm_global, ocnrtm_global
          read (nio) totrunin_ave, prec_ave, evap_ave, qchan2, qchocn2
#if (defined COUP_CSM)
          read (nio) ocnrof_vec     
#endif
       endif
#if (defined SPMD)
       call scatter_data_from_master (spmdtemp, totrunin_ave, clmlevel=namec)
       totrunin_ave(begc:endc) = spmdtemp(begc:endc)
       call scatter_data_from_master (spmdtemp, prec_ave, clmlevel=namec)
       prec_ave(begc:endc) = spmdtemp(begc:endc)
       call scatter_data_from_master (spmdtemp, evap_ave, clmlevel=namec)
       evap_ave(begc:endc) = spmdtemp(begc:endc)
       call scatter_data_from_master (spmdtemp, qchocn2, clmlevel=namec)
       qchocn2(begc:endc) = spmdtemp(begc:endc)
       call scatter_data_from_master (spmdtemp, qchan2, clmlevel=namec)
       qchan2(begc:endc) = spmdtemp(begc:endc)
#endif
       do k = begc,endc
          clmpointer(ic_cwf_qchan2)%val(k)%rp = qchan2(k)
          clmpointer(ic_cwf_qchocn2)%val(k)%rp = qchocn2(k)
       end do
    endif

    if (flag == 'write') then
#if (defined SPMD)
       spmdtemp(begc:endc) = totrunin_ave(begc:endc)
       call gather_data_to_master (spmdtemp, totrunin_ave, clmlevel=namec)
       spmdtemp(begc:endc) = prec_ave(begc:endc)
       call gather_data_to_master (spmdtemp, prec_ave, clmlevel=namec)
       spmdtemp(begc:endc) = evap_ave(begc:endc)
       call gather_data_to_master (spmdtemp, evap_ave, clmlevel=namec)
       spmdtemp(begc:endc) = qchocn2(begc:endc)
       call gather_data_to_master (spmdtemp, qchocn2, clmlevel=namec)
       spmdtemp(begc:endc) = qchan2(begc:endc)
       call gather_data_to_master (spmdtemp, qchan2, clmlevel=namec)
#endif
       if (masterproc) then
          write (nio) volr
          write (nio) fluxout
          write (nio) ncount_rtm, ncount_global, yrold
          write (nio) prec_global, evap_global, runlnd_global, &
               runrtm_global, volrtm_global, ocnrtm_global
          write (nio) totrunin_ave, prec_ave, evap_ave, qchan2, qchocn2
#if (defined COUP_CSM)
          write (nio) ocnrof_vec    
#endif
       endif
    endif

    return
  end subroutine restart_rtm

#endif

end module RtmMod


