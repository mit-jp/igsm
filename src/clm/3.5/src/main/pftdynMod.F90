#include <misc.h>
#include <preproc.h>

module pftdynMod

!---------------------------------------------------------------------------
!BOP
!
! !MODULE: pftdynMod
!
! !USES:
  use spmdMod
  use clmtype
  use decompMod   , only : gsmap_lnd_gdc2glo,perm_lnd_gdc2glo
  use decompMod   , only : get_proc_bounds
  use ncdio
  use clm_varsur  , only : pctspec
  use clm_varpar  , only : max_pft_per_col
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun
!
! !DESCRIPTION:
! Determine pft weights at current year using dynamic landuse datasets.
! ASSUMES that only have one dynamic landuse dataset.
!
! !PUBLIC TYPES:
  implicit none
  private
  save
  public :: pftdyn_init
  public :: pftdyn_interp
  public :: pftdyn_wbal_init
  public :: pftdyn_wbal
!
! !REVISION HISTORY:
! Created by Gordon Bonan, Sam Levis and Mariana Vertenstein
! Updated by Erwan Monier (07/15/2011)
! The daily interpolation is changed to an annual interpolation (for missing
! years in the data file) since PFTs do not change on a daily basis, but most
! likely during winter after harvest. If simulation year is before the data
! year range, then use first year of data file, if simulation year is past
! the data year range, then use last year of data file.
! And if rampYear_dynpft > 0, then set pft to ramped year.
!
!EOP
!
! ! PRIVATE TYPES
  real(r8), parameter :: days_per_year = 365._r8
  integer , pointer   :: yearspft(:)
  real(r8), pointer   :: wtpft1(:,:)   
  real(r8), pointer   :: wtpft2(:,:)   
  real(r8), pointer   :: wtcol_old(:)
  integer :: nt1
  integer :: nt2
  integer :: ncid
!---------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftdyn_init
!
! !INTERFACE:
  subroutine pftdyn_init()
!
! !DESCRIPTION:
! Initialize dynamic landuse dataset (position it to the right time samples to 
! that bound the initial model date
!
! !USES:
    use clm_time_manager, only : get_curr_date
    use clm_varctl  , only : fpftdyn, rampYear_dynpft
    use clm_varpar  , only : lsmlon, lsmlat, numpft
    use fileutils   , only : getfil
    use spmdGathScatMod, only : gather_data_to_master
    use clm_varcon  , only : istsoil
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: i,j,m,n,g,p,l                   ! indices
    integer  :: ntimes                          ! number of input time samples
    real(r8) :: sumpct                          ! sum for error check
    integer  :: varid                           ! netcdf ids
    integer  :: year                            ! year (0, ...) for nstep+1
    integer  :: mon                             ! month (1, ..., 12) for nstep+1
    integer  :: day                             ! day of month (1, ..., 31) for nstep+1
    integer  :: sec                             ! seconds into current date for nstep+1
    integer  :: ier                             ! error status
    logical  :: found                           ! true => input dataset bounding dates found
    integer  :: begg,endg                       ! beg/end indices for land gridcells
    integer  :: begl,endl                       ! beg/end indices for land landunits
    integer  :: begc,endc                       ! beg/end indices for land columns
    integer  :: begp,endp                       ! beg/end indices for land pfts
    real(r8), pointer :: pctgla(:)          ! percent of gcell is glacier
    real(r8), pointer :: pctlak(:)          ! percent of gcell is lake
    real(r8), pointer :: pctwet(:)          ! percent of gcell is wetland
    real(r8), pointer :: pcturb(:)          ! percent of gcell is urbanized
    type(gridcell_type), pointer :: gptr        ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr         ! pointer to landunit derived subtype
    type(pft_type)     , pointer :: pptr         ! pointer to pft derived subtype
    character(len=256) :: locfn                 ! local file name
    character(len= 32) :: subname='pftdyn_init' ! subroutine name
 !-----------------------------------------------------------------------

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    allocate(pctgla(begg:endg),pctlak(begg:endg))
    allocate(pctwet(begg:endg),pcturb(begg:endg))

    ! Set pointers into derived type

    gptr => clm3%g
    lptr => clm3%g%l
    pptr => clm3%g%l%c%p

    ! pctspec must be saved between time samples
    ! position to first time sample - assume that first time sample must match
    ! starting date
    ! check consistency -  special landunits, grid, frac and mask
    ! only do this once

    ! read data PCT_PFT corresponding to correct year

    allocate(wtpft1(begg:endg,0:numpft), wtpft2(begg:endg,0:numpft), stat=ier)
    if (ier /= 0) then
       write(6,*)'pctpft_dyn_init allocation error for wtpft1, wtpft2'
       call endrun()
    end if

    allocate(wtcol_old(begp:endp),stat=ier)
    if (ier /= 0) then
       write(6,*)'pctpft_dyn_init allocation error for wtcol_old'
       call endrun()
    end if

! Initialize pptr%wtcol_old as pptr%wtcol. If reading from the initial condition or restart file
! this will save the pptr%wtcol from these files into pptr%wtcol_old and
! pptr%wtcol will be updated to the current weight from the surface data file in
! the pftdyn_interp subroutine

    do p = begp,endp
       l = pptr%landunit(p)
       if (lptr%itype(l) == istsoil) then
          wtcol_old(p)   = pptr%wtcol(p)
       end if
    end do

    if (masterproc) then

       ! Obtain file

       write (6,*) 'Attempting to read pft dynamic landuse data .....'
       call getfil (fpftdyn, locfn, 0)
       call check_ret(nf_open(locfn, 0, ncid), subname)

       ! Obtain pft years from dynamic landuse file
! Modified the name of variable read (time exists for other variables in that
! surface data file, so we use year) Erwan Monier (07/14/2011)

       call check_ret(nf_inq_dimid(ncid, 'year', varid), subname)
       call check_ret(nf_inq_dimlen(ncid, varid, ntimes), subname)

       ! Consistency check

       call check_dim(ncid, 'lsmpft', numpft+1)

    endif

    call clmmpi_bcast(ntimes,1,CLMMPI_INTEGER,0,clmmpicom,ier)

    allocate (yearspft(ntimes), stat=ier)
    if (ier /= 0) then
       write(6,*)'pctpft_dyn_init allocation error for yearspft'; call endrun()
    end if

    if (masterproc) then
       call check_ret(nf_inq_varid(ncid, 'year', varid), subname)
       call check_ret(nf_get_var_int(ncid, varid, yearspft), subname)
    endif

    call clmmpi_bcast(yearspft,ntimes,CLMMPI_INTEGER,0,clmmpicom,ier)

    ! Determine if current date spans the years
    ! If current year is less than first dynamic PFT timeseries year,
    ! then use the first year from dynamic pft file for both nt1 and nt2,
    ! forcing constant weights until the model year enters the dynamic
    ! pft dataset timeseries range.
    ! If current year is equal to or greater than the last dynamic pft
    ! timeseries year, then use the last year for both nt1 and nt2, 
    ! forcing constant weights for the remainder of the simulation.
    ! This mechanism permits the introduction of a dynamic pft period in the
    ! middle
    ! of a simulation, with constant weights before and after the dynamic
    ! period.

    if (rampYear_dynpft /= 0) then
       year = rampYear_dynpft
    else
       call get_curr_date(year, mon, day, sec)
    endif

    if (year < yearspft(1)) then
       nt1 = 1
       nt2 = 1
       found = .true.
    else if (year >= yearspft(ntimes)) then
       nt1 = ntimes
       nt2 = ntimes
       found = .true.
    else
       found = .false.
       do n = 1,ntimes-1 
          if (year == yearspft(n)) then
             nt1 = n
             nt2 = n
             found = .true.
          end if   
       end do
! If the year is not found in the pft data file, then find years that bracket
! the current
       if (.not. found) then
          nt2 = 1
          do while(yearspft(nt2) < year) 
             nt2 = nt2 + 1
          end do
          nt1 = nt2 - 1
       end if
    end if

    deallocate(pctgla,pctlak,pctwet,pcturb)

  end subroutine pftdyn_init

  
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftdyn_interp
!
! !INTERFACE:
  subroutine pftdyn_interp()
!
! !DESCRIPTION:
! Time interpolate dynamic landuse data to get pft weights for model time
!
! !USES:
    use clm_time_manager, only : get_curr_date, get_curr_calday
    use clm_varcon  , only : istsoil
    use clm_varpar  , only : numpft, lsmlon, lsmlat
    use clm_varctl  , only : rampYear_dynpft
!
! !ARGUMENTS:
    implicit none
!
!EOP
!
! !LOCAL VARIABLES:
    logical  :: update_pft       ! true => the pft needs to be updated to correct year
    logical  :: found            ! true => current year found in pft input dataset
    logical  :: first            ! true => first time step of initial, restart or branch run
    integer  :: i,j,m,p,l,g,n,c    ! indices
    integer  :: ntimes           ! number of input time samples
    integer  :: year             ! year (0, ...) for nstep+1
    integer  :: mon              ! month (1, ..., 12) for nstep+1
    integer  :: day              ! day of month (1, ..., 31) for nstep+1
    integer  :: sec              ! seconds into current date for nstep+1
    real(r8) :: cday             ! current calendar day (1.0 = 0Z on Jan 1)
    integer  :: ier              ! error status
    real(r8) :: deltay,fact              ! time interpolation factors
    integer  :: begg,endg                       ! beg/end indices for land gridcells
    integer  :: begl,endl                       ! beg/end indices for land landunits
    integer  :: begc,endc                       ! beg/end indices for land columns
    integer  :: begp,endp                       ! beg/end indices for land pfts
    real(r8), pointer :: wtpfttot1(:)           ! summation of pft weights for renormalization
    real(r8), pointer :: wtpfttot2(:)           ! summation of pft weights for renormalization
    type(gridcell_type), pointer :: gptr         ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr         ! pointer to landunit derived subtype
    type(pft_type)     , pointer :: pptr         ! pointer to pft derived subtype
    character(len=32) :: subname='pftdyn_interp' ! subroutine name
!-----------------------------------------------------------------------

    data first /.true./

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    ! Set pointers into derived type

    gptr => clm3%g
    lptr => clm3%g%l
    pptr => clm3%g%l%c%p

    ! Interpolat pctpft to current year - output in pctpft_intp
    ! Map interpolated pctpft to subgrid weights
    ! assumes that maxpatch_pft = numpft + 1, that each landunit has only 1 column, 
    ! SCAM and DGVM have not been defined, and that create_croplandunit = .false.

    ! Get current date or use rampYear_dynpft
    update_pft = .false.

    if (rampYear_dynpft /= 0) then
       year = rampYear_dynpft
       ! No need to update the pft since it is ramped at year rampYear_dynpft 
       ! and the weigths have already been set in surfrdMod
    else
       call get_curr_date(year, mon, day, sec)

    ! Update the year to read when the current year of the simulation is greater
    ! than the year at which the pft is read from the input data file

       if (nt1 == nt2 .and. year > yearspft(nt1)) then
          update_pft = .true.
       end if

       if (nt1 /= nt2 .and. year == yearspft(nt2)) then
          update_pft = .true.
       end if

    end if

! Update the pft if new year or if first time step of the run (whether initial,
! restart or branch run) in order to rectify the weights read from the initial
! or restart file.
    if (update_pft .or. first) then
       if (year >= yearspft(size(yearspft))) then
          nt1 = size(yearspft)
          nt2 = size(yearspft)
          found = .true.
       else
          found = .false.
          do n = 1,size(yearspft)-1
             if (year == yearspft(n)) then
                nt1 = n
                nt2 = n
                found = .true.
             end if   
          end do
! If the year is not found in the pft data file, then find years that bracket
! the current year of the simulation
          if (.not. found) then
             nt2 = 1
             do while(yearspft(nt2) < year) 
                nt2 = nt2 + 1
             end do
             nt1 = nt2 - 1
          end if
       end if
    ! Now that the times have been updated, update the pctpft (if the current year
    ! is found in the pft data file, then the wtpft1 and wtpft2 will be equal)

       call pftdyn_getdata(nt1, wtpft1, begg,endg,0,numpft)
       call pftdyn_getdata(nt2, wtpft2, begg,endg,0,numpft)

       do m = 0,numpft
          do g = begg,endg
             wtpft1(g,m) = wtpft1(g,m)/100._r8
             wtpft2(g,m) = wtpft2(g,m)/100._r8
          end do
       end do

       allocate(wtpfttot1(begc:endc),wtpfttot2(begc:endc))
       wtpfttot1(:) = 0._r8
       wtpfttot2(:) = 0._r8

!dir$ concurrent
!cdir nodep
       do p = begp,endp
          c = pptr%column(p)
          g = pptr%gridcell(p)
          l = pptr%landunit(p)
          if (lptr%itype(l) == istsoil) then
             m = pptr%itype(p)
             wtpfttot1(c) = wtpfttot1(c)+pptr%wtgcell(p)
             if (nt1 == nt2) then
! Then the current year exists in the pft input data file
! No need to interpolate
                pptr%wtgcell(p)   = wtpft1(g,m)
             else
! We interpolate the pft between years of missing pft data
                deltay = yearspft(nt2) - yearspft(nt1)
                fact = (yearspft(nt2) - year)/deltay
                pptr%wtgcell(p)   = wtpft2(g,m) + fact*(wtpft1(g,m)-wtpft2(g,m))
             end if
             pptr%wtlunit(p)   = pptr%wtgcell(p) / lptr%wtgcell(l)
             pptr%wtcol(p)     = pptr%wtgcell(p) / lptr%wtgcell(l)
             wtpfttot2(c) = wtpfttot2(c)+pptr%wtgcell(p)
          end if
       end do

       deallocate(wtpfttot1,wtpfttot2)

! Now that we have updated the pfts, we need to call the column-level
! water mass-balance correction subroutine
       call pftdyn_wbal()

! After the water mass-balance correction, we update pptr%wtcol_old
       do p = begp,endp
          l = pptr%landunit(p)
          if (lptr%itype(l) == istsoil) then
             wtcol_old(p)   = pptr%wtcol(p)
          end if
       end do


! End of the update
    first = .false.
    end if

  end subroutine pftdyn_interp

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftdyn_getdata
!
! !INTERFACE:
  subroutine pftdyn_getdata(ntime, pctpft, lb1,ub1,lb2,ub2)
!
! !DESCRIPTION:
! Obtain dynamic landuse data (pctpft) and make sure that
! percentage of PFTs sum to 100% cover for vegetated landunit
!
! !USES:
    use clm_varpar  , only : numpft, lsmlon, lsmlat
    use pftvarcon   , only : crop
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    integer , intent(in)  :: ntime
    integer , intent(in)  :: lb1,ub1,lb2,ub2
    real(r8), intent(out) :: pctpft(lb1:ub1,lb2:ub2)
!
!EOP
!
! !LOCAL VARIABLES:
! changes by Erwan start here
    integer  :: i,j,m,n,nl
! changes by Erwan end here
    integer  :: begg,endg         
! changes by Erwan start here
    real(r8) :: sumpct                     ! temporary
! changes by Erwan end here
    integer  :: start(4), count(4)                ! input sizes
    real(r8),pointer :: arrayl(:)                 ! temporary array
    character(len=32) :: subname='pftdyn_getdata' ! subroutine name
!-----------------------------------------------------------------------
    
    call get_proc_bounds(begg,endg)

    allocate(arrayl(begg:endg))
    do n = 0,numpft
       start(1) = 1
       count(1) = lsmlon
       start(2) = 1
       count(2) = lsmlat
       start(3) = n+1      ! dataset is 1:numpft+1, not 0:numpft
       count(3) = 1
       start(4) = ntime 
       count(4) = 1
       call ncd_iolocal(ncid, 'PCT_PFT', 'read', arrayl, begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, start, count)
       pctpft(begg:endg,n) = arrayl(begg:endg)
    enddo
    deallocate(arrayl)

! Changes by Erwan start here
! Scale the pctpft to ensure that the sum of the pctpft is equal to 100-pctspec
! (should be the case in the data file read, but make sure anyway)

    do nl = begg,endg
       if (pctspec(nl) < 100._r8) then
          sumpct = 0._r8
          do m = 0,numpft
             sumpct = sumpct + pctpft(nl,m)
          end do
          if (abs(sumpct+pctspec(nl)-100._r8) > 1.0e-6) then
             do m = 0,numpft
                pctpft(nl,m) = pctpft(nl,m) * (100._r8-pctspec(nl))/sumpct
             end do
          end if
       end if
    end do
! Changes by Erwan end here
    
  end subroutine pftdyn_getdata

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftdyn_wbal_init
!
! !INTERFACE:
  subroutine pftdyn_wbal_init()
!
! !DESCRIPTION:
! initialize the column-level mass-balance correction term.
! Called in every timestep.
!
! !USES:
!
! !ARGUMENTS:
    implicit none
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: begp, endp    ! proc beginning and ending pft indices
    integer  :: begc, endc    ! proc beginning and ending column indices
    integer  :: begl, endl    ! proc beginning and ending landunit indices
    integer  :: begg, endg    ! proc beginning and ending gridcell indices
    integer  :: c             ! indices
    type(column_type),   pointer :: cptr         ! pointer to column derived subtype
!-----------------------------------------------------------------------

    ! Set pointers into derived type

    cptr => clm3%g%l%c

    ! Get relevant sizes

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! set column-level canopy water mass balance correction flux
    ! term to 0 at the beginning of every timestep
    
!dir$ concurrent
!cdir nodep
    do c = begc,endc
       cptr%cwf%h2ocan_loss(c) = 0._r8
    end do
    
  end subroutine pftdyn_wbal_init

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftdyn_wbal
!
! !INTERFACE:
  subroutine pftdyn_wbal()
!
! !DESCRIPTION:
! modify pft-level state and flux variables to maintain water balance with
! dynamic pft-weights.
!
! !USES:
    use clm_varcon  , only : istsoil
    use clm_time_manager, only : get_step_size
!
! !ARGUMENTS:
    implicit none
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: begp, endp    ! proc beginning and ending pft indices
    integer  :: begc, endc    ! proc beginning and ending column indices
    integer  :: begl, endl    ! proc beginning and ending landunit indices
    integer  :: begg, endg    ! proc beginning and ending gridcell indices
    integer  :: pi,p,c,l,g    ! indices
    integer  :: ier           ! error code
    real(r8) :: dwt           ! change in pft weight (relative to column)
    real(r8) :: dtime         ! land model time step (sec)
    real(r8) :: init_h2ocan   ! initial canopy water mass
    real(r8) :: new_h2ocan    ! canopy water mass after weight shift
    real(r8), allocatable :: loss_h2ocan(:) ! canopy water mass loss due to weight shift
    type(landunit_type), pointer :: lptr         ! pointer to landunit derived subtype
    type(column_type),   pointer :: cptr         ! pointer to column derived subtype
    type(pft_type)   ,   pointer :: pptr         ! pointer to pft derived subtype
    character(len=32) :: subname='pftdyn_wbal' ! subroutine name
!-----------------------------------------------------------------------

    ! Set pointers into derived type

    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    ! Get relevant sizes

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! Allocate loss_h2ocan
    allocate(loss_h2ocan(begp:endp), stat=ier)
    if (ier /= 0) then
          write(6,*)subname,' allocation error for loss_h2ocan'; call endrun()
    end if

    ! Get time step

    dtime = get_step_size()

    ! set column-level canopy water mass balance correction flux
    ! term to 0 at the beginning of every weight-shifting timestep

!dir$ concurrent
!cdir nodep
    do c = begc,endc
       cptr%cwf%h2ocan_loss(c) = 0._r8
    end do

!dir$ concurrent
!cdir nodep

    do p = begp,endp
       l = pptr%landunit(p)
       loss_h2ocan(p) = 0._r8

       if (lptr%itype(l) == istsoil) then

          ! calculate the change in weight for the timestep
          dwt = pptr%wtcol(p)-wtcol_old(p)

          if (dwt > 0._r8) then
          
             ! if the pft gained weight, then the 
             ! initial canopy water state is redistributed over the
             ! new (larger) area, conserving mass.

             pptr%pws%h2ocan(p) = pptr%pws%h2ocan(p) * (wtcol_old(p)/pptr%wtcol(p))
          
          else if (dwt < 0._r8) then
          
             ! if the pft lost weight on the timestep, then the canopy water
             ! mass associated with the lost weight is directed to a 
             ! column-level flux term that gets added to the precip flux
             ! for every pft calculation in Hydrology1()
             
             init_h2ocan = pptr%pws%h2ocan(p) * wtcol_old(p)
             loss_h2ocan(p) = pptr%pws%h2ocan(p) * (-dwt)
             new_h2ocan = init_h2ocan - loss_h2ocan(p)
             if (abs(new_h2ocan) < 1e-8_r8) then
                new_h2ocan = 0._r8
                loss_h2ocan(p) = init_h2ocan
             end if
             if (pptr%wtcol(p) /= 0._r8) then  
                pptr%pws%h2ocan(p) = new_h2ocan/pptr%wtcol(p)
             else
                pptr%pws%h2ocan(p) = 0._r8
                loss_h2ocan(p) = init_h2ocan
             end if 
          

          end if
          
       end if
       
    end do

!dir$ nointerchange
    do pi = 1,max_pft_per_col
!dir$ concurrent
!cdir nodep
       do c = begc,endc
          if (pi <= cptr%npfts(c)) then
             p = cptr%pfti(c) + pi - 1
             cptr%cwf%h2ocan_loss(c) = cptr%cwf%h2ocan_loss(c) + loss_h2ocan(p)/dtime
          end if
       end do
    end do

    ! Deallocate loss_h2ocan
    deallocate(loss_h2ocan)
    
  end subroutine pftdyn_wbal

end module pftdynMod
