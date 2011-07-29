#include <misc.h>
#include <preproc.h>

#if (defined 2dato3dl) 

module setcxyMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: setcxyMod
!
! !DESCRIPTION: 
! Module containing routines to read precipitation and air-temperature distribution
! coefficients across latitude bands, for each month.
!
! This module and routines have been deliberately written (though not exclusively)
! for the 3D implementation of CLM within the MIT Integrated Global System Model
! (with a zonal atmosphere)
!
! Author: C. Adam Schlosser 5/2009
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use abortutils  , only : endrun
  use clmtype
  use clm_varpar  , only : lsmlon, lsmlat
!
! !PUBLIC TYPES
  implicit none
  save
  public
!
! !PUBLIC MEMBER FUNCTIONS
  public :: readmonthlycxy  !read monthly data a month
!
! !REVISION HISTORY
!
!EOP
!
! PRIVATE MEMBER FUNCTIONS
!
! PRIVATE TYPES
!  integer, private :: Months1 = -999  !saved month index
!----------------------------------------------------------------------- 

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: readMonthlyCxy
!
! !INTERFACE:
  subroutine readmonthlycxy(fcxy,monthp)
!
! !DESCRIPTION: 
! Read monthly Precipitation and Temperature Distribution data
!
! !USES
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype	
    use clm_varctl, only : fsurdat
    use clm_varpar, only : nlevsoi, numpft, &
                           maxpatch_pft, maxpatch_cft, maxpatch, &
                           npatch_urban, npatch_lake, npatch_wet, npatch_glacier
    use clm_varsur, only : wtxy, vegxy
    use fileutils, only : getfil
    use clm_varsur  , only : pctspec
    use ncdio
    use spmdMod              
    use clm_varctl,   only : fsurdat,scmlat, scmlon, single_column 
    use domainMod   , only : domain_type
    use decompMod    , only : get_proc_bounds, ldecomp
    use clm_varcon, only : istsoil

#if (defined SPMD)
    use spmdMod, only : masterproc, clmmpicom, CLMMPI_REAL8
#else
    use spmdMod, only : masterproc
#endif
    use clm_time_manager, only : get_nstep
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: fcxy    !file with precipitation distr. data
    integer, intent(in) :: monthp           !months to be interpolated (1 to 12)
!
! !REVISION HISTORY
! Created by Adam Schlosser
!
!EOP
!
! LOCAL VARIABLES
  character(len=256) :: locfn       !local file name
  integer :: i,j,k,m,n,pi,pt        !indices
  integer :: jn,js,ci,clst          !indices
  integer :: ncid,dimid,varid       !input netCDF id's
  integer :: beg4d(4),len4d(4)      !netCDF variable edges
  integer :: beg3d(3),len3d(3)      !netCDF variable edges
  integer :: ntim                   !number of input data time samples
  integer :: nlon_i                 !number of input data longitudes
  integer :: nlat_i                 !number of input data latitudes
  integer :: nlev_i                 !number of input data latitudes
  integer :: npft_i                 !number of input data pft types
  integer :: ier                    !error code 
  integer :: iwtyp                  !water type (soil,lake,wetland,glacier)
  real :: cxysatgrid(lsmlon,lsmlat,12)
  real :: cxypcpgrid(lsmlon,lsmlat,12)
  type(pft_type) , pointer :: p                    !pointer to derived subtype
  type(column_type)         , pointer :: c   !local pointer to derived subtype
  type(pft_pstate_type), pointer :: pps            !local pointers to derived subtypes
  type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
  type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
  type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
  type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
  integer :: latmin,latmax                         !beginning and ending latitudes
  integer :: begp, endp      ! per-proc beginning and ending pft indices
  integer :: begc, endc      ! per-proc beginning and ending column indices
  integer :: begl, endl      ! per-proc beginning and ending landunit indices
  integer :: begg, endg      ! per-proc gridcell ending gridcell indices
  integer , pointer :: itype(:)          ! landunit type
  integer :: iland           ! Land type indice
  character(len=32) :: subname = 'SETCXY:'     ! subroutine name


!-----------------------------------------------------------------------

  ! Set pointers into derived type

  gptr => clm3%g
  lptr => clm3%g%l
  cptr => clm3%g%l%c
  pptr => clm3%g%l%c%p
  itype => clm3%g%l%itype

  ! Determine processor bounds

  call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! ----------------------------------------------------------------------
    ! Open monthly vegetation file
    ! Read data and convert from [lsmlon] x [lsmlat] grid to patch data 
    ! ----------------------------------------------------------------------

    if (masterproc) then

       write (6,*) subname,'Attempting to read 2D atmos. to 3D land coefficients...'
	   
! ----------------------------------------------------------------------	   
!
!	CAS: 05/09:  For now - all netcdf read and file access statements commented out
!				 but kept as placeholders for eventual reformatting of input data.
! ----------------------------------------------------------------------	   

!       call getfil (fsurdat, locfn, 0)
!       call check_ret(nf_open(locfn, 0, ncid), subname)

!       call check_ret(nf_inq_dimid (ncid, 'lsmlon', dimid), subname)
!       call check_ret(nf_inq_dimlen(ncid, dimid, nlon_i), subname)

!       call check_ret(nf_inq_dimid (ncid, 'lsmlat', dimid), subname)
!       call check_ret(nf_inq_dimlen(ncid, dimid, nlat_i), subname)

!       call check_ret(nf_inq_varid(ncid, 'PFT', varid), subname)
!       call check_ret(nf_get_vara_int(ncid, varid, beg3d, len3d, pctpft), subname)

!       beg3d(1) = 1         ; len3d(1) = lsmlon
!       beg3d(2) = 1         ; len3d(2) = lsmlat
!       beg3d(3) = 1         ; len3d(3) = maxpatch_pft
 
!       call check_ret(nf_inq_varid(ncid, 'PCT_PFT', varid), subname)
!       call check_ret(nf_get_vara_real(ncid, varid, beg3d, len3d, pctpft), subname) 
!       call check_ret(nf_close(ncid), subname)

!	write (6,*) 'Configuring PCP and SAT distribution across the global land grid'
!       write (6,*) 'Attempting to read monthly Cxy data .....'
!       write (6,*) 'nstep = ',get_nstep(),' for month = ', monthp

!       call getfil (fpcp2pft, locfn, 0)
!       call check_ret(nf_open(locfn, 0, ncid), subname)

!       call check_ret(nf_inq_dimid (ncid, 'lsmlon', dimid), subname)
!       call check_ret(nf_inq_dimlen(ncid, dimid, nlon_i), subname)

!       if (nlon_i /= lsmlon) then
!          write(6,*)'ReadMonthlyPcp2Pft: parameter lsmlon= ', &
!               lsmlon,'does not equal input nlat_i= ',nlon_i
!          call endrun
!       endif

!       call check_ret(nf_inq_dimid (ncid, 'lsmlat', dimid), subname)
!       call check_ret(nf_inq_dimlen(ncid, dimid, nlat_i), subname)

!       if (nlat_i /= lsmlat) then
!          write(6,*)'ReadMonthlyPcp2Pft: parameter lsmlat= ', &
!               lsmlat,'does not equal input nlat_i= ',nlat_i
!          call endrun
!       endif

!       call check_ret(nf_inq_dimid (ncid, 'lsmpft', dimid), subname)
!       call check_ret(nf_inq_dimlen(ncid, dimid, npft_i), subname)

!       if (npft_i /= maxpatch_pft) then
!          write(6,*)'ReadMonthlyPcp2Pft: parameter maxpatch_pft = ', &
!               maxpatch_pft,'does not equal input npft_i= ',npft_i
!          call endrun
!       endif

!       call check_ret(nf_inq_dimid (ncid, 'time', dimid), subname)
!       call check_ret(nf_inq_dimlen(ncid, dimid, ntim), subname)

!       beg4d(1) = 1         ; len4d(1) = nlon_i
!       beg4d(2) = 1         ; len4d(2) = nlat_i
!       beg4d(3) = 1         ; len4d(3) = npft_i
!       beg4d(4) = monthp    ; len4d(4) = 1

!       beg3d(1) = 1         ; len3d(1) = nlon_i
!       beg3d(2) = 1         ; len3d(2) = nlat_i
!       beg3d(3) = monthp    ; len3d(3) = 1

       !write (6,*) 'ReadPCP2PFT: maxpatch = ', maxpatch_pft, ' NPFT_I: ',npft_i

!       call check_ret(nf_inq_varid(ncid, 'PCP2LND', varid), subname)
!       call check_ret(nf_get_vara_real(ncid, varid, beg3d, len3d, pcp2land), subname)

!       call check_ret(nf_inq_varid(ncid, 'PCP2PFT', varid), subname)
!       call check_ret(nf_get_vara_real(ncid, varid, beg4d, len4d, pcp2pftxy ), subname)

    fl4pcp = trim(igsmdir) // 'data4tem' // flyr
    print *,'Opening data files: '
    print *, fl4clm
    print *, 'Opening CLM forcing:',fl4clm
    call opnfil(fl4clm,20,'u')


	do m = 1,monthp
          read (icxytunit,*) cxysatgrd
	  read (icxypunit,*) cxypcpgrd
        enddo

        write (6,*) 'setcxyMod: PCP Array'
        write (6,*) cxypcpgrd
        write (6,*) 'setcxyMod: SAT Array'
        write (6,*) cxysatgrd

!		Map coefficients from the lat/lon grid to the appropiate clm vector elements
!	    NOTE: THIS CODE ASSUMES THAT THE CLM AND CXY GRIDS ARE ALREADY COMPATIBLE!!!

       !write (6,*) 'SetCxy - Beginning and ending pft indices: ',begp,endp
       do pi = begp,endp
          i = ldecomp%gdc2i(pptr%gridcell(pi))
          j = ldecomp%gdc2j(pptr%gridcell(pi))
!
!  CAS:	CXYSAT AND CXYPCP ARE DECLARED IN CLMTYPE.F90
!
          pptr%pps%cxysat(pi) = cxysatgrd(i,j)
          pptr%pps%cxypcp(pi) = cxypcpgrd(i,j)			 
       end do

!      call check_ret(nf_close(ncid), subname)

       write (6,*) subname,'Successfully read monthly Cxy distr. data for'
       write (6,*) 'month ', monthp
       write (6,*)

    endif ! end of if-masterproc if block

#if ( defined SPMD )
    ! pass surface data to all processors
    call clmmpi_bcast (p%pps%cxysat, size(p%pps%cxysat), CLMMPI_REAL8, 0, clmmpicom, ier)
    call clmmpi_bcast (p%pps%cxypcp, size(p%pps%cxypcp), CLMMPI_REAL8, 0, clmmpicom, ier)
#endif

    return
  end subroutine readMonthlyPcp2Pft

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: setcxy
!
! !INTERFACE:
!  subroutine setcxy
!
! !DESCRIPTION: 
! Dynamically allocate memory and set to signaling NaN
!
! !USES
!    use nanMod
!
! !ARGUMENTS:
!    implicit none
!
! !REVISION HISTORY
!
!EOP
!-----------------------------------------------------------------------
!
!
!
end module setcxyMod

#else

!The following is only here since empty file won't compile
subroutine setcxy 
  write(6,*) 'SETCXY: this routine should not be called'
  return
end subroutine setcxy2pft

#endif

