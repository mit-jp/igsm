#include <misc.h>
#include <preproc.h>

#if (defined PCP2PFT) 

module pcp2pftMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: pcp2pftMod
!
! !DESCRIPTION: 
! Module containing routines to read precipitation distribution
! parameters across the PFT types, for each month.
!
! This module and routines have been deliberately written (though not exclusively)
! for the zonal implementation of CLM with the MIT Integrated Global System Model
!
! Author: C. Adam Schlosser 3/2005
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
!
! !PUBLIC TYPES
  implicit none
  save
  public
  real :: pcp2land(lsmlon,lsmlat)
!
! !PUBLIC MEMBER FUNCTIONS
  public :: readmonthlypcp2pft  !read monthly precipitation data a month
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
! !IROUTINE: readMonthlyPcp2Pft
!
! !INTERFACE:
  subroutine readMonthlyPcp2Pft (fpcp2pft,monthp)
!
! !DESCRIPTION: 
! Read monthly Precipitation Distribution data
!
! !USES
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype	
    use clm_varctl, only : fsurdat
    use clm_varpar, only : lsmlon, lsmlat, maxpatch_pft, maxpatch 
    use fileutils, only : getfil
    use clmpoint
#if (defined SPMD)
    use spmdMod, only : masterproc, mpicom, MPI_REAL8
#else
    use spmdMod, only : masterproc
#endif
    use time_manager, only : get_nstep
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: fpcp2pft    !file with precipitation distr. data
    integer, intent(in) :: monthp           !months to be interpolated (1 to 12)
!
! !REVISION HISTORY
! Created by Sam Levis
!
!EOP
!
! LOCAL VARIABLES
    character(len=256) :: locfn       !local file name
    integer :: i,j,k,m,pi,pt          !indices
    integer :: jn,js                  !indices
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
    real :: pcp2pftxy(lsmlon,lsmlat,maxpatch_pft)    !read from input files
    real :: srrgtn,srrgts                            !Surrogate values of pcp2pft ratio 
    real :: denom                     !denominator term
    real(r8):: sumpcp2pft(lsmlon,lsmlat)             !read from input files
    integer :: pft(lsmlon,lsmlat,maxpatch_pft)       !array for pctpft
    real(r8):: pctpft(lsmlon,lsmlat,maxpatch_pft)    !array for pft values
    real(r8):: sumpctpft(lsmlon,lsmlat)              !array for summed pctpft
    type(pft_type) , pointer :: p                    !pointer to derived subtype
    type(pft_pstate_type), pointer :: pps            !local pointers to derived subtypes
    integer :: begp, endp                            !beginning and ending pft indices
    integer :: latmin,latmax                         !beginning and ending latitudes

!-----------------------------------------------------------------------

    ! ----------------------------------------------------------------------
    ! Open monthly vegetation file
    ! Read data and convert from [lsmlon] x [lsmlat] grid to patch data 
    ! ----------------------------------------------------------------------

    if (masterproc) then

       write (6,*) 'Attempting to read surface boundary data .....'

       call getfil (fsurdat, locfn, 0)
       call wrap_open(locfn, 0, ncid)

       call wrap_inq_dimid  (ncid, 'lsmlon', dimid)
       call wrap_inq_dimlen (ncid, dimid, nlon_i)

       call wrap_inq_dimid  (ncid, 'lsmlat', dimid)
       call wrap_inq_dimlen (ncid, dimid, nlat_i)

       call wrap_inq_dimid  (ncid, 'nlevsoi', dimid)
       call wrap_inq_dimlen (ncid, dimid, nlev_i)

       call wrap_inq_dimid  (ncid, 'lsmpft', dimid)
       call wrap_inq_dimlen (ncid, dimid, npft_i)

       call wrap_inq_varid (ncid, 'PFT', varid)
       call wrap_get_var_int (ncid, varid, pft)

       call wrap_inq_varid (ncid, 'PCT_PFT', varid)
       call wrap_get_var_realx (ncid, varid, pctpft)

       write (6,*) 'Successfully read surface boundary data'
       write (6,*)

       call wrap_close(ncid)

       write (6,*) 'Configuring PCP distribution across PFTs'

       write (6,*) 'Attempting to read monthly pcp2pft data .....'
       write (6,*) 'nstep = ',get_nstep(),' for month = ', monthp

       call getfil (fpcp2pft, locfn, 0)
       call wrap_open(locfn, 0, ncid)

       call wrap_inq_dimid  (ncid, 'lsmlon', dimid)
       call wrap_inq_dimlen (ncid, dimid, nlon_i)
       if (nlon_i /= lsmlon) then
          write(6,*)'ReadMonthlyPcp2Pft: parameter lsmlon= ', &
               lsmlon,'does not equal input nlat_i= ',nlon_i
          call endrun
       endif

       call wrap_inq_dimid  (ncid, 'lsmlat', dimid)
       call wrap_inq_dimlen (ncid, dimid, nlat_i)
       if (nlat_i /= lsmlat) then
          write(6,*)'ReadMonthlyPcp2Pft: parameter lsmlat= ', &
               lsmlat,'does not equal input nlat_i= ',nlat_i
          call endrun
       endif

       call wrap_inq_dimid  (ncid, 'lsmpft', dimid)
       call wrap_inq_dimlen (ncid, dimid, npft_i)
       if (npft_i /= maxpatch_pft) then
          write(6,*)'ReadMonthlyPcp2Pft: parameter maxpatch_pft = ', &
               maxpatch_pft,'does not equal input npft_i= ',npft_i
          call endrun
       endif

       call wrap_inq_dimid  (ncid, 'time', dimid)
       call wrap_inq_dimlen (ncid, dimid, ntim)

          beg4d(1) = 1         ; len4d(1) = nlon_i
          beg4d(2) = 1         ; len4d(2) = nlat_i
          beg4d(3) = 1         ; len4d(3) = npft_i
          beg4d(4) = monthp    ; len4d(4) = 1

          beg3d(1) = 1         ; len3d(1) = nlon_i
          beg3d(2) = 1         ; len3d(2) = nlat_i
          beg3d(3) = monthp    ; len3d(3) = 1

!CAS          call wrap_inq_varid (ncid, 'PCP2LND', varid)
!CAS          call wrap_get_vara_realx (ncid, varid, beg3d, len3d, pcp2land )

          call wrap_inq_varid (ncid, 'PCP2PFT', varid)
          call wrap_get_vara_realx (ncid, varid, beg4d, len4d, pcp2pftxy )

          ! store data directly in clmtype structure
          ! only vegetated pfts have nonzero values

          sumpcp2pft = 0.
          sumpctpft = 0.
          begp = pfts1d%beg
          endp = pfts1d%end
          !write (6,*) 'ReadMonthlyPcp2Pft - Beginning and ending pft indices: ',begp,endp
          do pi = begp,endp
             p => ppoint%pft(pi)%p
             i = pfts1d%ixy(pi)
             j = pfts1d%jxy(pi)
             m = pfts1d%itypveg(pi)
             iwtyp = pfts1d%itypwat(pi)
             pt = m+1
!CAS
!CAS	Loop through pft array to find consistent mapping between pcp2pft and pctpft
!CAS
             do k = 1, maxpatch_pft
               if ( pft(i,j,k) == m ) exit
             enddo
             if (m <= maxpatch_pft .and. iwtyp == 1 ) then ! vegetated patch
                if ( pcp2pftxy(i,j,pt) == 0 ) then
                 write(6,*)'ReadMonthlyPcp2Pft: Non-zero ratio not found ', &
                       'at i = ',i,' and j = ',j, 'for PFT = ',pt,'  ...  ',&
                       '% area coverage',pfts1d%wtlnd(pi), &
                       'Stopping run in CLM'
                 call endrun
                endif
                !write(6,*)'ReadMonthlyPcp2Pft: pcp2pft = ',pcp2pftxy(i,j,pt), &
                !      'at i = ',i,' and j = ',j, 'for PFT = ',m,'  ...  ',&
                !      '% area coverage',pctpft(i,j,k),' PCTPFT value of PFT: ',pft(i,j,k)
                !write(6,*)'ReadMonthlyPcp2Pft: WATER TYPE = ',iwtyp
                !p%pps%pcp2pft = pcp2pftxy(i,j,pt)
                sumpcp2pft(i,j) = sumpcp2pft(i,j) + pcp2pftxy(i,j,pt)
                sumpctpft(i,j) = sumpctpft(i,j) + pctpft(i,j,k)
             else                        ! non-vegetated patch
                p%pps%pcp2pft = 1.
             end if
          end do
          do pi = begp,endp
            p => ppoint%pft(pi)%p
            i = pfts1d%ixy(pi)
            j = pfts1d%jxy(pi)
            m = pfts1d%itypveg(pi)
            iwtyp = pfts1d%itypwat(pi)
            pt = m + 1
!CAS
!CAS	Loop through pft array to find consistent mapping between pcp2pft and pctpft
!CAS
             do k = 1, maxpatch_pft
               if ( pft(i,j,k) == m ) exit
             enddo
              if ( m <= maxpatch_pft .and. iwtyp == 1 ) then ! vegetated patch
               if (pctpft(i,j,k) /= 0 ) then
                   denom = (sumpcp2pft(i,j)*pctpft(i,j,k))/sumpctpft(i,j)
                   p%pps%pcp2pft = pcp2pftxy(i,j,pt)/denom
                   !write(6,*)'ReadMonthlyPcp2Pft: at i = ',i,' j = ',j,' k =',k
                   !write(6,*)'for PFT = ',m,' and check ',pft(i,j,k)
                   !write(6,*)'SUMPCP2PFT @ column',pi,' = ',sumpcp2pft(i,j)
                   !write(6,*)'PCP2PFT @ column',pi,' = ',pcp2pftxy(i,j,pt)
                   !write(6,*)'PCTPFT @ column',pi,' = ',pctpft(i,j,k)
                   !write(6,*)'SCALE PCP2PFT @ column',pi,' = ', p%pps%pcp2pft
               else
                 write(6,*)'ReadMonthlyPcp2Pft: PCTPFT = 0 ', &
                       'at i = ',i,' and j = ',j, 'for PFT = ',m,'.',&
                       'Should be non-zero: Stopping CLM'
                 call endrun
               endif
              else
!CAS                 write(6,*)'ReadMonthlyPcp2Pft: Non-vegetated land patch:',&
!CAS                       'at i = ',i,' and j = ',j, 'for PFT = ',pt,'.',&
!CAS                       'Setting PCP2PFT distr. ratio to 1'
                 p%pps%pcp2pft = 1.
              endif
          end do

       call wrap_close(ncid)

       write (6,*) 'Successfully read monthly precipitation distr. data for'
       write (6,*) 'month ', monthp
       write (6,*)

    endif ! end of if-masterproc if block


#if ( defined SPMD )
    ! pass surface data to all processors

    call mpi_bcast (p%pps%pcp2pft, size(p%pps%pcp2pft), MPI_REAL8, 0, mpicom, ier)
#endif

    return
  end subroutine readMonthlyPcp2Pft

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: monthpcp2pft_ini
!
! !INTERFACE:
!  subroutine monthpcp2pft_ini ()
!
! !DESCRIPTION: 
! Dynamically allocate memory and set to signaling NaN
!
! !USES
!    use nanMod
!    use clmtype, only : pfts1d
!
! !ARGUMENTS:
!    implicit none
!
! !REVISION HISTORY
!
!EOP
!-----------------------------------------------------------------------
!
!    allocate (pcp2pft(pfts1d%num))
!
!    pcp2pft(:) = nan
!
!    return
!  end subroutine monthpcp2pft_ini
!
end module pcp2pftMod

#else

!The following is only here since empty file won't compile
subroutine pcp2pft 
  write(6,*) 'PCP2PFT: this routine should not be called'
  return
end subroutine pcp2pft

#endif

