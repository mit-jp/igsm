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
  use abortutils  , only : endrun
  use clmtype
  use clm_varpar  , only : lsmlon, lsmlat
!
! !PUBLIC TYPES
  implicit none
  save
  public
!  real :: pcp2land(lsmlon,lsmlat)
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
    character(len=*), intent(in) :: fpcp2pft            !file with precipitation distr. data
    integer, intent(in) :: monthp                       !months to be interpolated (1 to 12)
!
! !REVISION HISTORY
! Created by Sam Levis
!
!EOP
!
! LOCAL VARIABLES
  character(len=256) :: locfn                           !local file name
  integer :: i,j,k,m,n,pi,pt                            !indices
  integer :: jn,js,ci,clst                              !indices
  integer :: ncid,dimid,varid                           !input netCDF id's
  integer :: beg4d(4),len4d(4)                          !netCDF variable edges
  integer :: beg3d(3),len3d(3)                          !netCDF variable edges
  integer :: ntim                                       !number of input data time samples
  integer :: nlon_i                                     !number of input data longitudes
  integer :: nlat_i                                     !number of input data latitudes
  integer :: nlev_i                                     !number of input data latitudes
  integer :: npft_i                                     !number of input data pft types
  integer :: ier                                        !error code 
  integer :: iwtyp                                      !water type (soil,lake,wetland,glacier)
  real :: pcp2pftxy(lsmlon,lsmlat,numpft)               !read from input files
  real :: srrgtn,srrgts                                 !Surrogate values of pcp2pft ratio 
  real :: denom                                         !denominator term
  real(r8):: sumpcp2pft(lsmlon,lsmlat)                  !read from input files
  integer :: pft(lsmlon,lsmlat,maxpatch_pft)            !array for pctpft
  real :: pctpft(lsmlon,lsmlat,maxpatch_pft)            !array for pft values
  real(r8), pointer :: wtcol(:)
  real(r8):: sumpctpft(lsmlon,lsmlat)                   !array for summed pctpft
  type(pft_type) , pointer :: p                         !pointer to derived subtype
  type(column_type)         , pointer :: c              !local pointer to derived subtype
  type(pft_pstate_type), pointer :: pps                 !local pointers to derived subtypes
  type(gridcell_type), pointer :: gptr                  ! pointer to gridcell derived subtype
  type(landunit_type), pointer :: lptr                  ! pointer to landunit derived subtype
  type(column_type)  , pointer :: cptr                  ! pointer to column derived subtype
  type(pft_type)     , pointer :: pptr                  ! pointer to pft derived subtype
  integer :: latmin,latmax                              !beginning and ending latitudes
  integer :: begp, endp                                 ! per-proc beginning and ending pft indices
  integer :: begc, endc                                 ! per-proc beginning and ending column indices
  integer :: begl, endl                                 ! per-proc beginning and ending landunit indices
  integer :: begg, endg                                 ! per-proc gridcell ending gridcell indices
  integer , pointer :: itype(:)                         ! landunit type
  integer :: iland                                      ! Land type indice
  character(len=32) :: subname = 'pcp2pft'              ! subroutine name

!-----------------------------------------------------------------------

  ! Set pointers into derived type

  gptr => clm3%g
  lptr => clm3%g%l
  cptr => clm3%g%l%c
  pptr => clm3%g%l%c%p
  itype => clm3%g%l%itype
  wtcol   => clm3%g%l%c%p%wtcol

  ! Determine processor bounds

  call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! ----------------------------------------------------------------------
    ! Open monthly vegetation file
    ! Read data and convert from [lsmlon] x [lsmlat] grid to patch data 
    ! ----------------------------------------------------------------------

    if (masterproc) then

       write (6,*) 'Attempting to read surface boundary data .....'

       call getfil (fsurdat, locfn, 0)
       call check_ret(nf_open(locfn, 0, ncid), subname)

!       call check_ret(nf_inq_dimid (ncid, 'lsmlon', dimid), subname)
!       call check_ret(nf_inq_dimlen(ncid, dimid, nlon_i), subname)

!       call check_ret(nf_inq_dimid (ncid, 'lsmlat', dimid), subname)
!       call check_ret(nf_inq_dimlen(ncid, dimid, nlat_i), subname)

!       call check_ret(nf_inq_dimid (ncid, 'nlevsoi', dimid), subname)
!       call check_ret(nf_inq_dimlen(ncid, dimid, nlev_i), subname)

       call check_ret(nf_inq_dimid (ncid, 'lsmpft', dimid), subname)
       call check_ret(nf_inq_dimlen(ncid, dimid, npft_i), subname)

!       call check_ret(nf_inq_varid(ncid, 'PFT', varid), subname)
!       call check_ret(nf_get_vara_int(ncid, varid, beg3d, len3d, pctpft), subname)

       beg3d(1) = 1         ; len3d(1) = lsmlon
       beg3d(2) = 1         ; len3d(2) = lsmlat
       beg3d(3) = 1         ; len3d(3) = maxpatch_pft

       !write (6,*) 'ReadPCP2PFT: maxpatch = ', maxpatch_pft, ' NPFT_I: ',npft_i
 
       call check_ret(nf_inq_varid(ncid, 'PCT_PFT', varid), subname)
       call check_ret(nf_get_vara_real(ncid, varid, beg3d, len3d, pctpft), subname)
 
       !write (6,*) 'Successfully read surface boundary data'
       !write (6,*) pctpft

       call check_ret(nf_close(ncid), subname)

       write (6,*) 'Configuring PCP distribution across PFTs'

       write (6,*) 'Attempting to read monthly pcp2pft data .....'
       write (6,*) 'nstep = ',get_nstep(),' for month = ', monthp

       call getfil (fpcp2pft, locfn, 0)
       call check_ret(nf_open(locfn, 0, ncid), subname)

       call check_ret(nf_inq_dimid (ncid, 'lsmlon', dimid), subname)
       call check_ret(nf_inq_dimlen(ncid, dimid, nlon_i), subname)

       if (nlon_i /= lsmlon) then
          write(6,*)'ReadMonthlyPcp2Pft: parameter lsmlon= ', &
               lsmlon,'does not equal input nlat_i= ',nlon_i
          call endrun
       endif

       call check_ret(nf_inq_dimid (ncid, 'lsmlat', dimid), subname)
       call check_ret(nf_inq_dimlen(ncid, dimid, nlat_i), subname)

       if (nlat_i /= lsmlat) then
          write(6,*)'ReadMonthlyPcp2Pft: parameter lsmlat= ', &
               lsmlat,'does not equal input nlat_i= ',nlat_i
          call endrun
       endif

       call check_ret(nf_inq_dimid (ncid, 'lsmpft', dimid), subname)
       call check_ret(nf_inq_dimlen(ncid, dimid, npft_i), subname)

!       if (npft_i /= numpft) then
!          write(6,*)'ReadMonthlyPcp2Pft: parameter numpft = ', &
!               numpft,'does not equal input npft_i= ',npft_i
!          call endrun
!       endif

       call check_ret(nf_inq_dimid (ncid, 'time', dimid), subname)
       call check_ret(nf_inq_dimlen(ncid, dimid, ntim), subname)

       beg4d(1) = 1         ; len4d(1) = nlon_i
       beg4d(2) = 1         ; len4d(2) = nlat_i
       beg4d(3) = 1         ; len4d(3) = npft_i
       beg4d(4) = monthp    ; len4d(4) = 1

       beg3d(1) = 1         ; len3d(1) = nlon_i
       beg3d(2) = 1         ; len3d(2) = nlat_i
       beg3d(3) = monthp    ; len3d(3) = 1

       !write (6,*) 'ReadPCP2PFT: maxpatch = ', maxpatch_pft, ' NPFT_I: ',npft_i

!       call check_ret(nf_inq_varid(ncid, 'PCP2LND', varid), subname)
!       call check_ret(nf_get_vara_real(ncid, varid, beg3d, len3d, pcp2land), subname)

       call check_ret(nf_inq_varid(ncid, 'PCP2PFT', varid), subname)
       call check_ret(nf_get_vara_real(ncid, varid, beg4d, len4d, pcp2pftxy ), subname)

       !write (6,*) 'ReadMonthlyPcp2Pft: PCP2PFT Array'
       !write (6,*) pcp2pftxy

       ! store data directly in clmtype structure
       ! only vegetated pfts have nonzero values

       sumpcp2pft = 0.
       sumpctpft = 0.
       !write (6,*) 'ReadMonthlyPcp2Pft - Beginning and ending pft indices: ',begp,endp
          do pi = begp,endp
             i = ldecomp%gdc2i(pptr%gridcell(pi))
             j = ldecomp%gdc2j(pptr%gridcell(pi))
             m = pptr%itype(pi)
             iland = clm3%g%l%c%p%landunit(pi)
             iwtyp = clm3%g%l%itype(iland)
             pt = m + 1
             !write (6,*) 'For PFT INDICE: ',pi,' PFT TYPE = ', m
             !write (6,*) 'I = ',i,'  J = ',j,'  ITYP = ',iwtyp, '  ISTSOIL = ',istsoil
!CAS
!CAS	Loop through pft array to find consistent mapping between pcp2pft and pctpft
!CAS
             if (m <= numpft .and. iwtyp == istsoil ) then ! vegetated patch
               if ( pcp2pftxy(i,j,pt) == 0 .and. pt <= 15 ) then
                 pcp2pftxy(i,j,pt) = wtcol(pi)
                 !write(6,*)'ReadMonthlyPcp2Pft: Non-zero ratio not found ', &
                     !'at i = ',i,' and j = ',j, 'for PFT = ',pt,'  ...  ',&
                     !'PCP2PFT Values for latitude:',&
                     !pcp2pftxy(i,j,:),&
                     !'Configuring pcp2pft ratio to 1:'
                     !pcp2pftxy(i,j,pt) = pctpft(i,j,pt)/100.
               endif
               !write (6,*) 'pcp2pftxy is now set to: ', pcp2pftxy(i,j,pt)
               !write(6,*)'ReadMonthlyPcp2Pft: pcp2pft = ',pcp2pftxy(i,j,pt), &
               !      'at i = ',i,' and j = ',j, 'for PFT = ',m,'  ...  ',&
               !      '% area coverage',pctpft(i,j,k),' PCTPFT value of PFT: ',pft(i,j,k)
               !write(6,*)'ReadMonthlyPcp2Pft: WATER TYPE = ',iwtyp
!CAS
!CAS	AR5 LAND USE
!CAS	Addiing provision for crop and pasture types (pft16 and pft17)
!CAS
               if (pt >= 16) then
	         sumpcp2pft(i,j) = sumpcp2pft(i,j) + wtcol(pi)
       	       else
                 sumpcp2pft(i,j) = sumpcp2pft(i,j) + pcp2pftxy(i,j,pt)
	       endif
             else                        ! non-vegetated patch
                pptr%pps%pcp2pft(pi) = 1.
             end if
          end do
          !write (6,*) 'ReadMonthlyPcp2Pft:  SUMPCTPFT'
          !write (6,*) sumpctpft
          !write (6,*) 'ReadMonthlyPcp2Pft:  SUMPCP2PFT'
          !write (6,*) sumpcp2pft
          srrgtn = 0.
          clst = clm3%g%l%c%p%column(begp)
          do pi = begp,endp
            i = ldecomp%gdc2i(pptr%gridcell(pi))
            j = ldecomp%gdc2j(pptr%gridcell(pi))
            m = pptr%itype(pi)
            iland = clm3%g%l%c%p%landunit(pi)
            ci = clm3%g%l%c%p%column(pi)
            if ( ci /= clst ) then
              !if (srrgtn < 0.99 .or. srrgtn > 1.01 ) then
                !print *,'PCP2PFT: Aggregate weighting for latitude band not equal to 1.'
                !print *,'PCP2PFT: Aggregate weighting = ',srrgtn
                !print *,'PCP2PFT: at columnn = ',clst
              !endif
              srrgtn = 0.
            endif
            clst = ci
            n = clm3%g%l%c%npfts(ci)
            iwtyp = clm3%g%l%itype(iland)
            pt = m + 1
!CAS
!CAS	Loop through pft array to find consistent mapping between pcp2pft and pctpft
!CAS
            if ( m <= numpft .and. iwtyp == istsoil ) then ! vegetated patch
              !if (pctpft(i,j,pt) /= 0 .and. pt <= 15) then
              if (pt <= 15) then
                 !denom = (sumpcp2pft(i,j)*pctpft(i,j,pt))/sumpctpft(i,j)
		 denom = (sumpcp2pft(i,j)*wtcol(pi))
                 pptr%pps%pcp2pft(pi) = pcp2pftxy(i,j,pt)/denom
                 !write(6,*)'ReadMonthlyPcp2Pft: at i = ',i,' j = ',j,' pt =',pt
                 !write(6,*)'for PFT = ',m,' and check ',pft(i,j,k)
                 !write(6,*)'SUMPCP2PFT @ column',pi,' = ',sumpcp2pft(i,j)
                 !write(6,*)'PCP2PFT @ column',pi,' = ',pcp2pftxy(i,j,pt)
                 !write(6,*)'PCTPFT @ column',pi,' = ',pctpft(i,j,k)
                 !write(6,*)'PCP2PFT @ column',pi,' = ', pptr%pps%pcp2pft(pi)
              elseif (pt >= 16) then
                 !write (6,*) 'Pasture and/or Crop @ latitude ',j,' for PFT',pt,' Setting pcp2pft=1'
		 denom = (sumpcp2pft(i,j)*wtcol(pi))
                 pptr%pps%pcp2pft(pi) = wtcol(pi)/denom
	      else
                 write(6,*)'ReadMonthlyPcp2Pft: PCTPFT = 0 ', &
                    'at i = ',i,' and j = ',j, 'for PFT = ',pt,'.',&
                    'Should be non-zero: Stopping CLM'
                    call endrun
              endif
              srrgtn = srrgtn + (pptr%pps%pcp2pft(pi)*wtcol(pi))
            else
!CAS              write(6,*)'ReadMonthlyPcp2Pft: Non-vegetated land patch:',&
!CAS                       'at i = ',i,' and j = ',j, 'for PFT = ',pt,'.',&
!CAS                       'Setting PCP2PFT distr. ratio to 1'
              pptr%pps%pcp2pft(pi) = 1.
            endif
          end do

       call check_ret(nf_close(ncid), subname)

       write (6,*) 'Successfully read monthly precipitation distr. data for'
       write (6,*) 'month ', monthp
       write (6,*)

    endif ! end of if-masterproc if block


#if ( defined SPMD )
    ! pass surface data to all processors

    call clmmpi_bcast (p%pps%pcp2pft, size(p%pps%pcp2pft), CLMMPI_REAL8, 0, clmmpicom, ier)
#endif

    return
  end subroutine readMonthlyPcp2Pft

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: pcp2pftnoluc
!
! !INTERFACE:
!  subroutine pcp2pftnoluc ()
!
! !DESCRIPTION: 
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

