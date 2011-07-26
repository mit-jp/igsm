#include <misc.h>
#include <preproc.h>

module mvegFileMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: mvegFileMod
!
! !DESCRIPTION: 
! Module containing routines to read monthly vegetation data for two 
! months, interpolate monthly vegetation data and allocate and 
! initialize dynamic memory
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
!
! !PUBLIC TYPES
  implicit none
  save
  real(r8), public :: timwt(2)              !time weights for month 1 and month 2
  real(r8), public, allocatable :: mlai1(:) !lai for interpolation (month 1)
  real(r8), public, allocatable :: mlai2(:) !lai for interpolation (month 2)
  real(r8), public, allocatable :: msai1(:) !sai for interpolation (month 1)
  real(r8), public, allocatable :: msai2(:) !sai for interpolation (month 2)
  real(r8), public, allocatable :: mhvt1(:) !top vegetation height for interpolation (month 1)
  real(r8), public, allocatable :: mhvt2(:) !top vegetation height for interpolation (month 2)
  real(r8), public, allocatable :: mhvb1(:) !bottom vegetation height for interpolation(month 1)
  real(r8), public, allocatable :: mhvb2(:) !bottom vegetation height for interpolation(month 2)
  save
!
! !PUBLIC MEMBER FUNCTIONS
  public  :: interpMonthlyVeg        !interpolate monthly vegetation data
  public  :: monthveg_ini            !allocate and initialize dynamic memory
!
! !REVISION HISTORY
!
!EOP
!
! PRIVATE MEMBER FUNCTIONS
  private :: readMonthlyVegetation   !read monthly vegetation data for two months
!
! PRIVATE TYPES
  integer, private  :: InterpMonths1 = -999   !saved month index
!----------------------------------------------------------------------- 

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: interpMonthlyVeg
!
! !INTERFACE:
  subroutine interpMonthlyVeg (fveg, kmo, kda)
!
! !DESCRIPTION: 
! Determine if 2 new months of data are to be read
!
! !USES
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: fveg  !file with monthly vegetation data
    integer, intent(in) :: kmo             !month (1, ..., 12)
    integer, intent(in) :: kda             !day of month (1, ..., 31)
!
! !REVISION HISTORY
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES
    real(r8):: t           !a fraction: kda/ndaypm
    integer :: it(2)       !month 1 and month 2 (step 1)
    integer :: months(2)   !months to be interpolated (1 to 12)
    integer, dimension(12) :: ndaypm= &
         (/31,28,31,30,31,30,31,31,30,31,30,31/) !days per month
!-----------------------------------------------------------------------

    t = (kda-0.5) / ndaypm(kmo)
    it(1) = t + 0.5
    it(2) = it(1) + 1
    months(1) = kmo + it(1) - 1
    months(2) = kmo + it(2) - 1
    if (months(1) <  1) months(1) = 12
    if (months(2) > 12) months(2) = 1
    timwt(1) = (it(1)+0.5) - t
    timwt(2) = 1.-timwt(1)

    if (InterpMonths1 /= months(1)) then
       call readMonthlyVegetation (fveg, kmo, kda, months)
       InterpMonths1 = months(1)
    endif

    return
  end subroutine interpMonthlyVeg

!=======================================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: readMonthlyVegetation
!
! !INTERFACE:
  subroutine readMonthlyVegetation (fveg, kmo, kda, months)
!
! !DESCRIPTION: 
! Read monthly vegetation data for two consec. months
!
! !USES
    use clmtype	
    use clm_varpar, only : lsmlon, lsmlat, maxpatch_pft, maxpatch 
    use fileutils, only : getfil
#if (defined SPMD)
    use spmdMod, only : masterproc, mpicom, MPI_REAL8
#else
    use spmdMod, only : masterproc
#endif
    use time_manager, only : get_nstep
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: fveg    !file with monthly vegetation data
    integer, intent(in) :: kmo              !month (1, ..., 12)
    integer, intent(in) :: kda              !day of month (1, ..., 31)
    integer, intent(in) :: months(2)        !months to be interpolated (1 to 12)
!
! !REVISION HISTORY
! Created by Sam Levis
!
!EOP
!
! LOCAL VARIABLES
    character(len=256) :: locfn       !local file name
    integer :: i,j,k,m,pi             !indices
    integer :: ncid,dimid,varid       !input netCDF id's
    integer :: beg4d(4),len4d(4)      !netCDF variable edges
    integer :: ntim                   !number of input data time samples
    integer :: nlon_i                 !number of input data longitudes
    integer :: nlat_i                 !number of input data latitudes
    integer :: npft_i                 !number of input data pft types
    integer :: ier                    !error code 
    real(r8):: laixy(lsmlon,lsmlat,maxpatch_pft) !lai read from input files
    real(r8):: saixy(lsmlon,lsmlat,maxpatch_pft) !sai read from input files
    real(r8):: hvtxy(lsmlon,lsmlat,maxpatch_pft) !top vegetation height
    real(r8):: hvbxy(lsmlon,lsmlat,maxpatch_pft) !bottom vegetation height
!-----------------------------------------------------------------------

    ! ----------------------------------------------------------------------
    ! Open monthly vegetation file
    ! Read data and convert from [lsmlon] x [lsmlat] grid to patch data 
    ! ----------------------------------------------------------------------

    if (masterproc) then

       write (6,*) 'Attempting to read monthly vegetation data .....'
       write (6,*) 'nstep = ',get_nstep(),' month = ',kmo,' day = ',kda

       call getfil (fveg, locfn, 0)
       call wrap_open(locfn, 0, ncid)

       call wrap_inq_dimid  (ncid, 'lsmlon', dimid)
       call wrap_inq_dimlen (ncid, dimid, nlon_i)
       if (nlon_i /= lsmlon) then
          write(6,*)'ReadMonthlyVegetation: parameter lsmlon= ', &
               lsmlon,'does not equal input nlat_i= ',nlon_i
          call endrun
       endif

       call wrap_inq_dimid  (ncid, 'lsmlat', dimid)
       call wrap_inq_dimlen (ncid, dimid, nlat_i)
       if (nlat_i /= lsmlat) then
          write(6,*)'ReadMonthlyVegetation: parameter lsmlat= ', &
               lsmlat,'does not equal input nlat_i= ',nlat_i
          call endrun
       endif

       call wrap_inq_dimid  (ncid, 'lsmpft', dimid)
       call wrap_inq_dimlen (ncid, dimid, npft_i)
       if (npft_i /= maxpatch_pft) then
          write(6,*)'ReadMonthlyVegetation: parameter maxpatch_pft = ', &
               maxpatch_pft,'does not equal input npft_i= ',npft_i
          call endrun
       endif

       call wrap_inq_dimid  (ncid, 'time', dimid)
       call wrap_inq_dimlen (ncid, dimid, ntim)

       do k=1,2   !loop over months and read vegetated data

          beg4d(1) = 1         ; len4d(1) = nlon_i
          beg4d(2) = 1         ; len4d(2) = nlat_i
          beg4d(3) = 1         ; len4d(3) = npft_i
          beg4d(4) = months(k) ; len4d(4) = 1

          call wrap_inq_varid (ncid, 'MONTHLY_LAI', varid)
          call wrap_get_vara_realx (ncid, varid, beg4d, len4d, laixy )

          call wrap_inq_varid (ncid, 'MONTHLY_SAI', varid)
          call wrap_get_vara_realx (ncid, varid, beg4d, len4d, saixy )

          call wrap_inq_varid (ncid, 'MONTHLY_HEIGHT_TOP', varid)
          call wrap_get_vara_realx (ncid, varid, beg4d, len4d, hvtxy)

          call wrap_inq_varid (ncid, 'MONTHLY_HEIGHT_BOT', varid)
          call wrap_get_vara_realx (ncid, varid, beg4d, len4d, hvbxy)

          ! store data directly in clmtype structure
          ! only vegetated pfts have nonzero values

          do pi = 1, pfts1d%num
             i = pfts1d%ixy(pi)
             j = pfts1d%jxy(pi)
             m = pfts1d%mxy(pi)
             if (m <= maxpatch_pft) then ! vegetated patch
                if (k == 1) then
                   mlai1(pi) = laixy(i,j,m)
                   msai1(pi) = saixy(i,j,m)
                   mhvt1(pi) = hvtxy(i,j,m)
                   mhvb1(pi) = hvbxy(i,j,m)
                else !if (k == 2)
                   mlai2(pi) = laixy(i,j,m)
                   msai2(pi) = saixy(i,j,m)
                   mhvt2(pi) = hvtxy(i,j,m)
                   mhvb2(pi) = hvbxy(i,j,m)
                end if
             else                        ! non-vegetated patch
                if (k == 1) then
                   mlai1(pi) = 0. 
                   mlai2(pi) = 0.
                   msai1(pi) = 0. 
                   msai2(pi) = 0.
                else !if (k == 2)
                   mhvt1(pi) = 0. 
                   mhvt2(pi) = 0.
                   mhvb1(pi) = 0. 
                   mhvb2(pi) = 0.
                endif
             end if
          end do

       end do   ! end of loop over months

       call wrap_close(ncid)

       write (6,*) 'Successfully read monthly vegetation data for'
       write (6,*) 'month ', months(1), ' and month ', months(2)
       write (6,*)

    endif ! end of if-masterproc if block


#if ( defined SPMD )
    ! pass surface data to all processors

    call mpi_bcast (mlai1, size(mlai1), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (mlai2, size(mlai2), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (msai1, size(msai1), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (msai2, size(msai2), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (mhvt1, size(mhvt1), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (mhvt2, size(mhvt2), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (mhvb1, size(mhvb1), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (mhvb2, size(mhvb2), MPI_REAL8, 0, mpicom, ier)
#endif

    return
  end subroutine readMonthlyVegetation

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: monthveg_ini
!
! !INTERFACE:
  subroutine monthveg_ini ()
!
! !DESCRIPTION: 
! Dynamically allocate memory and set to signaling NaN
!
! !USES
    use nanMod
    use clmtype, only : pfts1d
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY
!
!EOP
!-----------------------------------------------------------------------

    allocate (mlai1(pfts1d%num))
    allocate (mlai2(pfts1d%num))
    allocate (msai1(pfts1d%num))
    allocate (msai2(pfts1d%num))
    allocate (mhvt1(pfts1d%num))
    allocate (mhvt2(pfts1d%num))
    allocate (mhvb1(pfts1d%num))
    allocate (mhvb2(pfts1d%num))

    mlai1(:) = nan
    mlai2(:) = nan
    msai1(:) = nan
    msai2(:) = nan
    mhvt1(:) = nan
    mhvt2(:) = nan
    mhvb1(:) = nan
    mhvb2(:) = nan

    return
  end subroutine monthveg_ini

end module mvegFileMod
