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
  subroutine readMonthlyPcp2Pft(fpcp2pft,monthp)
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
! Created by Adam Schlosser
!
!EOP
!
! LOCAL VARIABLES
  character(len=256) :: locfn                           !local file name
  integer, parameter :: npcp = 16                       ! # of pcp2pft types
  integer :: i,j,k,m,n,pi,pt                            !indices
  integer :: jn,js,ci                                   !indices
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
  !real :: pcp2pftxy(lsmlon,lsmlat,npcp,12)              ! PCP2PFT input array
  real :: pcp2pftxy(lsmlon,lsmlat,npcp)                 ! PCP2PFT input array
  real :: denom                                         !denominator term
  real(r8):: sumpcp2pft(lsmlon,lsmlat)                  !PCP2PFT accumulator
  real(r8):: sumpcp2pft1(lsmlon,lsmlat)                 !Buffer for PCP2PFT accumulator
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
  logical :: pcp2pftchk                                 ! pcp2pft logical check

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

       beg3d(1) = 1         ; len3d(1) = lsmlon
       beg3d(2) = 1         ; len3d(2) = lsmlat
       beg3d(3) = 1         ; len3d(3) = maxpatch_pft

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

       !if ( maxpatch_pft /= 17 ) then
       !  write (6,*) 'ReadPCP2PFT: maxpatch = ', maxpatch_pft, ' NPFT_I: ',npft_i,' numpft = ', numpft
       !  write (6,*) 'ReadPCP2PFT: maxpatch_pft must be equal to 17 in this IGSM configuration'
       !  call endrun
       !endif

       call check_ret(nf_inq_varid(ncid, 'PCP2PFT', varid), subname)
       call check_ret(nf_get_vara_real(ncid, varid, beg4d, len4d, pcp2pftxy ), subname)

       !write (6,*) 'ReadMonthlyPcp2Pft: PCP2PFT Array'
       !write (6,*) pcp2pftxy

       ! store data directly in clmtype structure
       ! only vegetated pfts have nonzero values

       sumpcp2pft = 0.0
       sumpctpft = 0.0
       pptr%pps%pcp2pft = 0.0
       !write (6,*) 'ReadMonthlyPcp2Pft - Beginning and ending pft indices: ',begp,endp
        do pi = begp,endp
          i = ldecomp%gdc2i(pptr%gridcell(pi))
          j = ldecomp%gdc2j(pptr%gridcell(pi))
          m = pptr%itype(pi)
          iland = clm3%g%l%c%p%landunit(pi)
          ci = clm3%g%l%c%p%column(pi)
          n = clm3%g%l%c%npfts(ci)
          iwtyp = clm3%g%l%itype(iland)
          pt = m + 1
!
!	Loop through pft array to find consistent mapping between pcp2pft and pctpft
!
          if ( iwtyp == istsoil ) then ! natural vegetated patch
            if (pt <= 15) then
              pptr%pps%pcp2pft(pi) = pcp2pftxy(i,j,pt)
              if ( wtcol(pi) > 0. .and. pcp2pftxy(i,j,pt) == 0. ) then
                write (6,*) 'PCP2PFT: WARNING - NO PCP2PFT VALUE AVAILABLE FOR VEGETATION COVER:'
                write (6,*) 'PCP2PFT: WTCOL',wtcol(pi),' @ J = ',j,' PFT = ',pt
                write (6,*) 'PCP2PFT: PCP2PFTXY',pcp2pftxy(i,j,pt)
                call endrun
              endif
              sumpcp2pft(i,j) = sumpcp2pft(i,j) + (pptr%pps%pcp2pft(pi)*wtcol(pi))
            else
!
!		Sum up area weight of pasture/crop for grid/latitude
!
              sumpctpft(i,j) = sumpctpft(i,j) + wtcol(pi)
            endif
          else
!                 write(6,*)'ReadMonthlyPcp2Pft: Non-vegetated land patch:',&
!                          'at i = ',i,' and j = ',j, 'for PFT = ',pt,'.',&
!                          'Setting PCP2PFT distr. ratio to 1'
            pptr%pps%pcp2pft(pi) = 1.
          endif
        end do
!
        do pi = begp,endp
          i = ldecomp%gdc2i(pptr%gridcell(pi))
          j = ldecomp%gdc2j(pptr%gridcell(pi))
          m = pptr%itype(pi)
          ci = clm3%g%l%c%p%column(pi)
          iwtyp = clm3%g%l%itype(iland)
          pt = m + 1
!
!	   Adjust crop/pasture pcp2pft to satisfy water balance:
!          This is achieved by divying residual precipitation (due to land use
!          from vegetated areas) over to crop/pasture by relative area.
!
!          NOTE: This translation of precipitation is only done for latitudes in
!                crop area exceeds 1e-5 in relative (vegetated) area.
!
          if (pt >= 16) then
               if (sumpctpft(i,j) > 1.0e-5 .and. sumpcp2pft(i,j) >= 0.0) then
                 pptr%pps%pcp2pft(pi) = (1-sumpcp2pft(i,j))/sumpctpft(i,j)
               else
                 pptr%pps%pcp2pft(pi) = 0.0
               endif
               if (sumpctpft(i,j) > 1.0e-5 .and. (sumpcp2pft(i,j)-1) > 1.0e-6 ) then
                 write (6,*) 'PCP2PFT: Insufficient PCP for Pasture/Crop @ latitude ',j,' for PFT ',pt
                 write (6,*) 'RESIDUAL PCP = ',1-sumpcp2pft(i,j)
                 write (6,*) 'CROP/PASTURE FRACTION = ',sumpctpft(i,j)
                 write (6,*) 'PCP2PFT = ',pptr%pps%pcp2pft(pi)
                 call endrun
               endif
          endif
        end do
        sumpcp2pft = 0.0
        do pi = begp,endp
          i = ldecomp%gdc2i(pptr%gridcell(pi))
          j = ldecomp%gdc2j(pptr%gridcell(pi))
          m = pptr%itype(pi)
          iland = clm3%g%l%c%p%landunit(pi)
          iwtyp = clm3%g%l%itype(iland)
          pt = m + 1
!
!	Loop through pft array to integrate precipitation over ALL vegetation (natural and crops)
!
          if ( iwtyp == istsoil ) then ! vegetated patch
            sumpcp2pft(i,j) = sumpcp2pft(i,j) + (pptr%pps%pcp2pft(pi)*wtcol(pi))
!
!       If CROP/PASTURE: CHECK THAT RESULTING PCP2PFT COEFF DOES NOT
!       EXCEED ANY OF THE PCP2PFT VALUES OF THE NATURAL VEGETATION
!
            if ( pt > 15 ) then
              pcp2pftchk = .false.
              do n = 1,15
                if (pptr%pps%pcp2pft(pi) < pcp2pftxy(i,j,n)) pcp2pftchk = .true.
              enddo
              if (.not.pcp2pftchk) then
                  write (6,*) 'PCP2PFT: Crop/Pasture value too high'
                  write (6,*) 'PCP2PFT(',pt,') =  ',pptr%pps%pcp2pft(pi)
                  write (6,*) 'PCP2PFT(',n,') =  ',pcp2pftxy(i,j,n)
                  call endrun
              endif
            endif
          endif
        end do
!
!       The following block of code checks for whether the PCP2PFT consistency
!       holds - that is, the sum of the product of PCP2PFT and WTCOL across a
!       latitude must equal 1 for potential vegetation.  If not, the PCP2PFT values
!       are rescaled accordingly (and not the wtcol values).
!
        do pi = begp,endp
          i = ldecomp%gdc2i(pptr%gridcell(pi))
          j = ldecomp%gdc2j(pptr%gridcell(pi))
          m = pptr%itype(pi)
          iland = clm3%g%l%c%p%landunit(pi)
          iwtyp = clm3%g%l%itype(iland)
          pt = m + 1
          if (sumpcp2pft(i,j) /= 1.0 .and. iwtyp == istsoil) then
            pptr%pps%pcp2pft(pi) = pptr%pps%pcp2pft(pi)/sumpcp2pft(i,j)
          endif
        enddo   

       !write (6,*) 'PCP2PFTMTH: SUMPCP2PFT'
       !write (6,*) sumpcp2pft

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
! !IROUTINE: pcp2pftinit
!
! !INTERFACE:
  subroutine pcp2pftinit(fpcp2pft)
!
! !DESCRIPTION: 
!
! !USES
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
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: fpcp2pft            !file with precipitation distr. data
!
! !REVISION HISTORY
! Created by Adam Schlosser
!
!EOP
!
! LOCAL VARIABLES
  character(len=256) :: locfn                           !local file name
  integer, parameter :: npcp = 16                       ! # of pcp2pft types
  integer :: i,j,k,m,n,pi,pt                            !indices
  integer :: jn,js,ci                                   !indices
  integer :: ncid,dimid,varid                           !input netCDF id's
  integer :: beg4d(4),len4d(4)                          !netCDF variable edges
  integer :: ntim                                       !number of input data time samples
  integer :: nlon_i                                     !number of input data longitudes
  integer :: nlat_i                                     !number of input data latitudes
  integer :: nlev_i                                     !number of input data latitudes
  integer :: npft_i                                     !number of input data pft types
  integer :: ier                                        !error code 
  integer :: iwtyp                                      !water type (soil,lake,wetland,glacier)
  real :: denom                                         !denominator term
  real(r8):: sumpcp2pft(lsmlon,lsmlat)                  !PCP2PFT accumulator
  integer :: pft(lsmlon,lsmlat,maxpatch_pft)            !array for pctpft
  real :: pctpft(lsmlon,lsmlat,maxpatch_pft)            !array for pft values
  real(r8), pointer :: wtcol(:)
  real(r8):: sumpctpft(lsmlon,lsmlat)                   !array for summed pctpft
  real :: pcp2pftxy(lsmlon,lsmlat,npcp,12)              ! PCP2PFT input array
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
  character(len=3) :: mthl(12)
  data (mthl(i),i=1,12) / 'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC' /
  
!
!EOP
!-----------------------------------------------------------------------
!
  ! Set pointers into derived type

  gptr => clm3%g
  lptr => clm3%g%l
  cptr => clm3%g%l%c
  pptr => clm3%g%l%c%p
  itype => clm3%g%l%itype
  wtcol => clm3%g%l%c%p%wtcol

  ! Determine processor bounds

  call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! ----------------------------------------------------------------------
    ! Open monthly vegetation file
    ! Read data and convert from [lsmlon] x [lsmlat] grid to patch data 
    ! ----------------------------------------------------------------------

     write (6,*) 'Initializing PCP distribution across PFTs'

     write (6,*) 'Attempting to read monthly pcp2pft data .....'

     call getfil (fpcp2pft, locfn, 0)
     call check_ret(nf_open(locfn, nf_write, ncid), subname)

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

     call check_ret(nf_inq_dimid (ncid, 'time', dimid), subname)
     call check_ret(nf_inq_dimlen(ncid, dimid, ntim), subname)

     write (6,*) 'PCP2PFTINIT: NTIM ',ntim

     beg4d(1) = 1         ; len4d(1) = nlon_i
     beg4d(2) = 1         ; len4d(2) = nlat_i
     beg4d(3) = 1         ; len4d(3) = npft_i
     beg4d(4) = 1         ; len4d(4) = 12

     !if ( maxpatch_pft /= 17 ) then
     !  write (6,*) 'ReadPCP2PFT: maxpatch = ', maxpatch_pft, ' NPFT_I: ',npft_i,' numpft = ', numpft
     !  write (6,*) 'ReadPCP2PFT: maxpatch_pft must be equal to 17 in this IGSM configuration'
     !  call endrun
     !endif

     call check_ret(nf_inq_varid(ncid, 'PCP2PFT', varid), subname)
     call check_ret(nf_get_vara_real(ncid, varid, beg4d, len4d, pcp2pftxy), subname)

     write (6,*) 'ReadMonthlyPcp2Pft: JAN PCP2PFT Array'
     write (6,*) pcp2pftxy(:,:,:,1)
     write (6,*) 'ReadMonthlyPcp2Pft: JUN PCP2PFT Array'
     write (6,*) pcp2pftxy(:,:,:,6)
     write (6,*) 'ReadMonthlyPcp2Pft: DEC PCP2PFT Array'
     write (6,*) pcp2pftxy(:,:,:,12)

     ! store data directly in clmtype structure
     ! only vegetated pfts have nonzero values

     do n = 1,12
       sumpcp2pft = 0.0
       write (6,*) 'Initial processing of pcp2pft data for: ',mthl(n)
       do pi = begp,endp
         i = ldecomp%gdc2i(pptr%gridcell(pi))
         j = ldecomp%gdc2j(pptr%gridcell(pi))
         m = pptr%itype(pi)
         iland = clm3%g%l%c%p%landunit(pi)
         iwtyp = clm3%g%l%itype(iland)
         pt = m + 1
!CAS
!CAS	Loop through pft array to integrate precipitation over potential vegetation
!CAS
         if ( iwtyp == istsoil .and. pt <= 15) then ! natural vegetated patch
            pptr%pps%pcp2pftmth(pi,n) = pcp2pftxy(i,j,pt,n)
            sumpcp2pft(i,j) = sumpcp2pft(i,j) + (pptr%pps%pcp2pftmth(pi,n)*wtcol(pi))
         endif
       end do
!
!       The following block of code checks for whether the PCP2PFT consistency
!       holds - that is, the sum of the product of PCP2PFT and WTCOL across a
!       latitude must equal 1 for potential vegetation.  If not, the PCP2PFT values
!       are rescaled accordingly (and not the wtcol values).
!
       do pi = begp,endp
         i = ldecomp%gdc2i(pptr%gridcell(pi))
         j = ldecomp%gdc2j(pptr%gridcell(pi))
         m = pptr%itype(pi)
         iland = clm3%g%l%c%p%landunit(pi)
         iwtyp = clm3%g%l%itype(iland)
         pt = m + 1
         if (sumpcp2pft(i,j) /= 1.0 .and. pt <= 15 .and. iwtyp == istsoil) then
           pptr%pps%pcp2pftmth(pi,n) = pptr%pps%pcp2pftmth(pi,n)/sumpcp2pft(i,j)
           pcp2pftxy(i,j,pt,n) = pptr%pps%pcp2pftmth(pi,n)
         endif
       enddo   
       write (6,*) 'PCP2PFTINIT: SUMPCP2PFT'
       write (6,*) sumpcp2pft
     enddo   

     call check_ret(nf_inq_varid(ncid, 'PCP2PFT', varid), subname)
     call check_ret(nf_put_vara_real(ncid, varid, beg4d, len4d, pcp2pftxy), subname)
     call check_ret(nf_close(ncid), subname)

     write (6,*) 'Initial processing of pcp2pft data completed'
     return
  end subroutine pcp2pftinit
end module pcp2pftMod

#else

!The following is only here since empty file won't compile
subroutine pcp2pft 
  write(6,*) 'PCP2PFT: this routine should not be called'
  return
end subroutine pcp2pft

#endif

