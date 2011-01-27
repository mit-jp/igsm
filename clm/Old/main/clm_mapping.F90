#include <misc.h>
#include <preproc.h>

module clm_mapping

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: clm_mapping
!
! !DESCRIPTION: 
! Initializes sub-grid mapping
!                              
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
#if (defined COUP_CAM)
  use clmtype
#endif
  use clm_varpar, only : lsmlon, lsmlat, maxpatch, maxpatch_pft
  use clm_varsur, only : numlon, landmask, area, latixy, longxy
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public clm_map   ! Initialize sub-grid mapping 
  public clm_map1d ! Initialize 1d mapping arrays
!
! !REVISION HISTORY:
! Created by Peter Thornton and Mariana Vertenstein 
!
!EOP
!----------------------------------------------------------------------- 

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_map
!
! !INTERFACE:
  subroutine clm_map (vegxy, wtxy) 
!
! !DESCRIPTION: 
! Initialize sub-grid mapping and allocates space for derived type hierarchy.
!
! !USES
#if (!defined COUP_CAM)
    use clmtype
#endif
    use lnd_grid, only : surface_grid_init, get_clump_info
#if (defined COUP_CAM)
    use lp_coupling, only : lp_coupling_init
#endif 
    use pftvarcon, only : noveg
#if (defined SPMD)
    use spmdMod, only : masterproc, spmd_init_arrays, gather_data_to_master, &
         proc_gridtot, proc_landtot, proc_coltot, proc_pfttot, npes
#else
    use spmdMod, only : masterproc
#endif
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: vegxy(lsmlon,lsmlat,maxpatch) !PFT type 
    real(r8), intent(in) :: wtxy(lsmlon,lsmlat,maxpatch)  !subgrid patch weights
!
! !REVISION HISTORY:
! Created by Peter Thornton and Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: i,j,m,n,gi,li,ci,pi    !indices
    integer  :: proc                   !SPMD processor id
    integer  :: nvegpatch,nspecpatch   !indices for pft and special cover types
    integer  :: nlandunits             !temporary index
    integer  :: ncolumns               !temporary index
    integer  :: npfts                  !temporary index
    integer  :: begg,endg              !beginning and ending 1d grid indices
    integer  :: begl,endl              !beginning and ending 1d land unit indices
    integer  :: begc,endc              !beginning and ending 1d column indices
    integer  :: begp,endp              !beginning and ending 1d pft indices
    integer  :: gindex                 !index into global 1d grid array
    integer  :: lindex                 !index into global 1d land unit array
    integer  :: cindex                 !index into global 1d land column array
    integer  :: pindex                 !index into global 1d land pft array
    real(r8) :: weight, vegpatchwt     !temporaries
    type(gridcell_type)       , pointer :: g   !pointer to derived subtype
    type(landunit_type)       , pointer :: l   !pointer to derived subtype
    type(column_type)         , pointer :: c   !pointer to derived subtype
    type(pft_type)            , pointer :: p   !pointer to derived subtype
    type(model_pstate_type)   , pointer :: mps !pointer to derived subtype
    type(gridcell_pstate_type), pointer :: gps !pointer to derived subtype
    type(column_pstate_type)  , pointer :: cps !pointer to derived subtype
    type(landunit_pstate_type), pointer :: lps !pointer to derived subtype
    type(pft_pstate_type)     , pointer :: pps !pointer to derived subtype
!------------------------------------------------------------------------

#if (defined SPMD)
    ! Initialize per-processor spmd arrays 

    call spmd_init_arrays
#endif

    call surface_grid_init(wtxy)

#if (defined SPMD)
    ! Determine total subgrid components for each process

    do proc=0,npes-1
       call get_clump_info (proc, proc_gridtot(proc), proc_landtot(proc), &
            proc_coltot(proc), proc_pfttot(proc))
    end do
#endif

#if (defined COUP_CAM)
    ! Initialize mapping between the atmosphere physics chunks and
    ! and the land gridcell clumps

    call lp_coupling_init()
#endif


    ! set up shorthand for processor beginning and ending 1d indices

    begg = grid1d%beg
    endg = grid1d%end
    begl = land1d%beg
    endl = land1d%end
    begc = cols1d%beg
    endc = cols1d%end
    begp = pfts1d%beg
    endp = pfts1d%end

    ! --------------------------------------------------------------------
    ! Allocate dynamic memory for the CLM derived type hierarchy
    ! For now, the vegetated patches will all be gathered on a single landunit,
    ! with each vegetated type having its own column on that landunit.  The
    ! special patches (urban, lake, wetland, glacier) each get 
    ! their own landunit having a single column and 1 non-vegetated pfts
    ! --------------------------------------------------------------------

    ! assign local pointer for simpler referencing
    mps => clm%mps
    mps%ngridcells = endg - begg + 1  

    ! Allocate for the array of gridcells to be handled by this process
    ! and for the arrays of area, weight and type for each gridcell.

    allocate(clm%g (mps%ngridcells))
    allocate(clm%ga(mps%ngridcells))
    allocate(clm%gw(mps%ngridcells))
    allocate(clm%gt(mps%ngridcells))

    gi = 0
    do j = 1, lsmlat
       do i = 1, numlon(j)
          if (landmask(i,j)==1) then
             gi = gi + 1                           
             if (gi>=begg .and. gi<=endg) then
                clm%g(gi-begg+1)%gps%ixy = i !longitude index
                clm%g(gi-begg+1)%gps%jxy = j !latitude index       
             end if
          endif
       end do
    end do

    ! set beggining 1d global indices
    ! (subtract 1 so that first index is beggrid for grid index, etc)
    gindex = begg - 1  
    lindex = begl - 1
    cindex = begc - 1
    pindex = begp - 1

    !loop over the land grid
    !assign local pointer for simpler referencing
    !note that nlandunits below refers to the number of landunits in
    !a gridcell, above it referred to the total number of landunits 
    !per processo

    do gi = 1,mps%ngridcells

       ! Set up pointers 
       g => clm%g(gi)					
       gps => g%gps

       !increment 1d global index
       gindex = gindex + 1
       gps%index1d = gindex

       ! Get grid indices
       i = gps%ixy
       j = gps%jxy

       ! Set area, weight, and type information for this gridcell.
       ! For now there is only one type of gridcell, value = 1
       ! Still need to resolve the calculation of area for the gridcell
       ! Note that the area, wt, and type are stored in an array at
       ! the higher level in hierarchy, so that the higher level has information
       ! about all members of the next lower level, and the same values
       ! are stored as scalars in each member of the lower level, so
       ! that each member has information about their own area, weight, and type.
       clm%ga(gi) = area(i,j)
!      clm%gw(gi) = clm%ga(gi)/mps%area  
       clm%gt(gi) = 1
       gps%area  = clm%ga(gi)
       gps%wt    = clm%gw(gi)
       gps%itype = clm%gt(gi)

       ! Assign pointer to higher level
       gps%mps => mps

       ! Count the number of vegetation patches
       ! and find the total weight for the landunit that will hold
       ! all the pfts. This code assumes that wtxy() values sum to
       ! 1.0 for all patches on a land gridcell.
       nvegpatch = 0
       vegpatchwt = 0.
       do m = 1, maxpatch_pft                           
          if (wtxy(i,j,m) > 0.) then
             nvegpatch = nvegpatch + 1
             vegpatchwt = vegpatchwt + wtxy(i,j,m)
          end if
       end do

       ! Count the number of special types (urban, lake, wetland, glacier).
       nspecpatch = 0
       do m = maxpatch_pft+1, maxpatch
          if (wtxy(i,j,m) > 0.) nspecpatch = nspecpatch + 1
       end do

       ! Set the total number of landunits to nspec
       ! Add 1 additional landunit if there are any vegetated patches
       nlandunits = nspecpatch
       if (nvegpatch > 0) nlandunits = nlandunits + 1
       gps%nlandunits = nlandunits

       ! Allocate for the array of landunits on this gridcell
       ! and for the arrays of area, weight and type for each landunit.
       ! pft landunits will have type 0, and the special type
       ! landunits will have type as set by npatch_urban, npatch_lake, etc.
       allocate(g%l(nlandunits))
       allocate(g%la(nlandunits))
       allocate(g%lw(nlandunits))
       allocate(g%lt(nlandunits))

       !*** Loop through landunits for this gridcell.***
       li = 0

       ! The first landunit contains all the vegetated patches
       if (nvegpatch > 0) then
          li = li + 1

          ! Set up pointers
          l => g%l(li)
          lps => l%lps

          ! Set special landunit flag to false
          lps%ifspecial = .false.

          ! Increment 1d global land index
          lindex = lindex + 1
          lps%index1d = lindex

          ! Set area, weight, and type information for this landunit
          ! the weights for each landunit on the grid cell must add
          ! to one when summed over the grid cell
          weight = vegpatchwt
          g%la(li) = clm%ga(gi) * weight
          g%lw(li) = weight
          g%lt(li) = 0
          lps%area  = g%la(li)
          lps%wt    = g%lw(li)
          lps%itype = g%lt(li)

          ! Set grid index and weight (relative to grid cell) 
          lps%ixy = i
          lps%jxy = j
          lps%wtxy = lps%area / gps%area

          ! Assign pointers to higher level
          l%a2ls => g%a2ls
          l%a2lf => g%a2lf
          l%lps%gps => g%gps

          ! Allocate memory for the array of columns on this landunit
          ! For the vegetated landunit, ncolumns = nvegpatch
          ncolumns = nvegpatch

          lps%ncolumns = ncolumns
          allocate(l%c(ncolumns))
          allocate(l%ca(ncolumns))
          allocate(l%cw(ncolumns))
          allocate(l%ct(ncolumns))

          ! Loop through regular (vegetated) patches, assign one column for each
          ! vegetated patch with non-zero weight. The weights for each column on
          ! the vegetated landunit must add to one when summed over the landunit,
          ! so the wtxy() values are taken relative to the total vegpatchwt
          ci = 0
          do m = 1, maxpatch_pft                           
             if (wtxy(i,j,m) > 0.) then
                ci = ci + 1

                ! Set up pointers
                c => l%c(ci) 
                cps => c%cps

                !increment 1d global column index
                cindex = cindex + 1
                cps%index1d = cindex

                ! Set area, weight, and type information for this column
                ! For now all columns have the same type, value = 1
                weight = wtxy(i,j,m) / vegpatchwt
                l%ca(ci) = g%la(li) * weight
                l%cw(ci) = weight
                l%ct(ci) = vegxy(i,j,m)  
                cps%area  = l%ca(ci)
                cps%wt    = l%cw(ci)
                cps%itype = vegxy(i,j,m) 

                ! Set grid index and weight (relative to grid cell) 
                cps%ixy = i
                cps%jxy = j
                cps%wtxy = cps%area / gps%area

                ! Assign pointers to higher level
                c%a2ls => g%a2ls
                c%a2lf => g%a2lf
                c%cps%lps => l%lps

                ! Allocate memory for the array of pfts on this column
                ! and for the arrays of area, weight and type for each pft.
                ! For now each column in the vegetated landunit is allocated
                ! a single pft. Later there will be the option for more than
                ! one pft per column.
                npfts = 1
                cps%npfts = npfts
                allocate(c%p(npfts))
                allocate(c%pa(npfts))
                allocate(c%pw(npfts))
                allocate(c%pt(npfts))

                ! Loop through the pfts for this column.
                ! For now there is only one pft for each column
                ! but later there will be the possibility of multiple
                ! pfts on each column
                do pi = 1, npfts										

                   ! Set up pointers
                   p => c%p(pi)
                   pps => p%pps

                   ! Increment 1d global pft index
                   pindex = pindex + 1
                   pps%index1d = pindex

                   ! Set area, weight (relative to column) and type information for this pft
                   ! For now, a single pft per column, so weight = 1
                   ! pft type comes from the m dimension of wtxy()
                   weight = 1.0/npfts
                   c%pa(pi) = l%ca(ci) * weight
                   c%pw(pi) = weight
                   c%pt(pi) = vegxy(i,j,m)
                   pps%area    = c%pa(pi)
                   pps%wt      = c%pw(pi)
                   pps%itype   = c%pt(pi)
                   pps%itypveg = c%pt(pi)

                   ! Set grid index, weight (relative to grid cell) 
                   ! and m index (needed for laixy, etc. reference)
                   pps%mxy = m
                   pps%ixy = i
                   pps%jxy = j
                   pps%wtxy = pps%area / gps%area

                   ! Assign pointers to higher level
                   p%a2ls => g%a2ls
                   p%a2lf => g%a2lf
                   pps%cps => c%cps

                end do   ! end loop through pfts
             end if   ! end if non-zero weight	
          end do   ! end loop through the possible vegetated patch indices

       end if   ! end if any vegetated patches

       ! Loop through the special landunits (urban, lake, wetland, glacier).
       if (nspecpatch > 0) then
          do m=maxpatch_pft+1, maxpatch
             if (wtxy(i,j,m) > 0.) then
                li = li + 1

                ! Set up pointers
                l => g%l(li)  
                lps => l%lps

                ! Set special landunit flag to true
                lps%ifspecial = .true.

                ! Increment 1d global land index
                lindex = lindex + 1
                lps%index1d = lindex

                ! Set area, weight, and type information for this landunit.
                weight = wtxy(i,j,m)
                g%la(li) = clm%ga(gi) * weight
                g%lw(li) = weight
                g%lt(li) = m
                lps%area = g%la(li)
                lps%wt   = g%lw(li)
                lps%itype= g%lt(li)

                ! Set grid index and weight (relative to grid cell) 
                lps%ixy = i
                lps%jxy = j
                lps%wtxy = lps%area / gps%area

                ! Assign pointers to higher level
                l%a2ls => g%a2ls
                l%a2lf => g%a2lf
                lps%gps => g%gps

                ! Allocate memory for the array of columns on this landunit
                ! For the special landunits, there is only one column 
                ! For now there is only one type of column, value = 1.
                ! Later, the age classes will be implemented on different
                ! columns within the same landunit, so the column type
                ! will correspond to an age class
                ncolumns = 1
                lps%ncolumns = ncolumns
                allocate(l%c(ncolumns))
                allocate(l%ca(ncolumns))
                allocate(l%cw(ncolumns))
                allocate(l%ct(ncolumns))

                ! Loop through columns for this landunit
                ! to set area, weight, and type information for this landunit
                ! We know that there is only one column for the special
                ! landunits, but the loop is included for consistency.
                do ci=1, ncolumns

                   ! Set up pointers
                   c => l%c(ci) 
                   cps => c%cps

                   ! Increment 1d global column index
                   cindex = cindex + 1
                   cps%index1d = cindex

                   ! Set area, weight, and type information for this column
                   ! For now all columns have the same type, value = 1
                   weight = 1.0/ncolumns
                   l%ca(ci) = g%la(li) * weight
                   l%cw(ci) = weight
                   l%ct(ci) = 1
                   cps%area  = l%ca(ci)
                   cps%wt    = l%cw(ci)
                   cps%itype = l%ct(ci)

                   ! Set grid index and weight (relative to grid cell) 
                   cps%ixy = i
                   cps%jxy = j
                   cps%wtxy = cps%area / gps%area

                   ! Assign pointers to higher level
                   c%a2ls => g%a2ls
                   c%a2lf => g%a2lf
                   cps%lps => l%lps

                   ! Set one non-vegetated pft for this column 
                   npfts = 1
                   cps%npfts = npfts
                   allocate(c%p(npfts))
                   allocate(c%pa(npfts))
                   allocate(c%pw(npfts))
                   allocate(c%pt(npfts))

                   ! Set up pointers for non-vegetated pft
                   pi = 1
                   p => c%p(pi) 
                   pps => p%pps

                   ! Increment 1d global pft index
                   pindex = pindex + 1
                   pps%index1d = pindex

                   ! Set area, weight (relative to column), and type information 
                   ! for this non-vegetated pft
                   weight = 1.0/npfts
                   c%pa(pi) = l%ca(ci) * weight
                   c%pw(pi) = weight
                   c%pt(pi) = noveg
                   pps%area    = c%pa(pi)
                   pps%wt      = c%pw(pi)
                   pps%itype   = c%pt(pi)
                   pps%itypveg = c%pt(pi)

                   ! Set grid index, weight (relative to grid cell) and  
                   ! m index (needed for laixy, etc. reference)
                   pps%mxy = m
                   pps%ixy = i
                   pps%jxy = j
                   pps%wtxy = pps%area / gps%area

                   ! Assign pointers to higher level
                   p%a2ls => g%a2ls
                   p%a2lf => g%a2lf
                   pps%cps => c%cps

                end do   ! end loop through ncolumns
             end if   ! end if this patch# has weight
          end do   ! end loop through special landunits
       end if   ! end if any special landunits

    end do  ! end loop grid cells

  end subroutine clm_map

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_map1d
!
! !INTERFACE:
  subroutine clm_map1d()
!
! !DESCRIPTION: 
! Set up 1d array of weights and indices for xy mapping
! Note: if DGVM is defined, weights are updated in DGVM mode 
!
! !USES:
#if (!defined COUP_CAM)
    use clmtype
#endif
    use clmpoint, only : cpoint, ppoint
#if (defined SPMD)
    use spmdMod, only : masterproc, spmd_init_arrays, gather_data_to_master
#else
    use spmdMod, only : masterproc
#endif
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: gi,li,ci,pi            !indices
    integer  :: gindex                 !index into global 1d gridcell array
    integer  :: lindex                 !index into global 1d landunit array
    integer  :: cindex                 !index into global 1d column array
    integer  :: pindex                 !index into global 1d pft array
    integer  :: ixy,jxy                !gridcell indices
    real(r8) :: londeg                 !gridcell longitude (degrees)
    real(r8) :: latdeg                 !gridcell latitude (degrees)
    integer  :: itypwat                !landunit water type 
    integer  :: begc,endc              !beginning and ending 1d column indices
    integer  :: begp,endp              !beginning and ending 1d pft indices
    integer  :: numg                   !total numger of grid cells
    integer  :: numl                   !total numger of land units
    integer  :: numc                   !total numger of columns
    integer  :: nump                   !total numger of pfts
    type(gridcell_type), pointer :: g  !pointer to derived subtype
    type(landunit_type), pointer :: l  !pointer to derived subtype
    type(column_type)  , pointer :: c  !pointer to derived subtype
    type(pft_type)     , pointer :: p  !pointer to derived subtype
#if (defined SPMD)
    real(r8), pointer :: rloc(:)       !temporaries for mpi gather
    integer , pointer :: iloc(:)       !temporaries for mpi gather
    real(r8), pointer :: rglob(:)      !temporaries for mpi gather
    integer , pointer :: iglob(:)      !temporaries for mpi gather
#endif
!-------------------------------------------------------------------------

    grid1d%name = 'gridcell'
    land1d%name = 'landunit'
    cols1d%name = 'column'
    pfts1d%name = 'pft'

    numg = grid1d%num
    numl = land1d%num
    numc = cols1d%num
    nump = pfts1d%num

    allocate (grid1d%wtxy(numg))
    allocate (grid1d%ixy(numg))
    allocate (grid1d%jxy(numg))
    allocate (grid1d%latdeg(numg))
    allocate (grid1d%londeg(numg))

    allocate (land1d%wtxy(numl))
    allocate (land1d%ixy(numl))
    allocate (land1d%jxy(numl))
    allocate (land1d%gindex(numl))
    allocate (land1d%latdeg(numl))
    allocate (land1d%londeg(numl))
    allocate (land1d%itypwat(numl))

    allocate (cols1d%wtxy(numc))
    allocate (cols1d%wtlnd(numc))
    allocate (cols1d%ixy(numc))
    allocate (cols1d%jxy(numc))
    allocate (cols1d%gindex(numc))
    allocate (cols1d%lindex(numc))
    allocate (cols1d%latdeg(numc))
    allocate (cols1d%londeg(numc))
    allocate (cols1d%itypwat(numc))

    allocate (pfts1d%wtxy(nump))
    allocate (pfts1d%wtlnd(nump))
    allocate (pfts1d%wtcol(nump))
    allocate (pfts1d%ixy(nump))
    allocate (pfts1d%jxy(nump))
    allocate (pfts1d%mxy(nump))
    allocate (pfts1d%gindex(nump))
    allocate (pfts1d%lindex(nump))
    allocate (pfts1d%cindex(nump))
    allocate (pfts1d%latdeg(nump))
    allocate (pfts1d%londeg(nump))
    allocate (pfts1d%itypwat(nump))
    allocate (pfts1d%itypveg(nump))

    do gi = 1, clm%mps%ngridcells
       g => clm%g(gi)					
       gindex = g%gps%index1d 
       ixy = g%gps%ixy
       jxy = g%gps%jxy
       latdeg = latixy(ixy,jxy)
       londeg = longxy(ixy,jxy)
       grid1d%wtxy(gindex) = 1._r8
       grid1d%ixy(gindex) = ixy
       grid1d%jxy(gindex) = jxy
       grid1d%latdeg(gindex) = latdeg
       grid1d%londeg(gindex) = londeg

       do li = 1, g%gps%nlandunits
          l => g%l(li)
          lindex = l%lps%index1d 
          itypwat = l%lps%itypwat
          land1d%ixy(lindex) = ixy
          land1d%jxy(lindex) = jxy
          land1d%latdeg(lindex) = latdeg
          land1d%londeg(lindex) = londeg
          land1d%gindex(lindex) = gindex
          land1d%itypwat(lindex) = itypwat
          land1d%wtxy(lindex) = l%lps%area / g%gps%area

          do ci = 1, l%lps%ncolumns
             c => l%c(ci) 
             cindex = c%cps%index1d 
             cols1d%ixy(cindex) = ixy
             cols1d%jxy(cindex) = jxy
             cols1d%latdeg(cindex) = latdeg
             cols1d%londeg(cindex) = londeg
             cols1d%gindex(cindex) = gindex
             cols1d%lindex(cindex) = lindex
             cols1d%itypwat(cindex) = itypwat
             cols1d%wtxy(cindex) = c%cps%area / g%gps%area
             cols1d%wtlnd(cindex) = c%cps%area / l%lps%area

             do pi = 1, c%cps%npfts
                p => c%p(pi)
                pindex = p%pps%index1d 
                pfts1d%ixy(pindex) = ixy
                pfts1d%jxy(pindex) = jxy
                pfts1d%mxy(pindex) = p%pps%mxy
                pfts1d%latdeg(pindex) = latdeg
                pfts1d%londeg(pindex) = londeg
                pfts1d%gindex(pindex) = gindex
                pfts1d%lindex(pindex) = lindex
                pfts1d%cindex(pindex) = cindex
                pfts1d%itypwat(pindex) = itypwat
                pfts1d%itypveg(pindex) = p%pps%itypveg
                pfts1d%wtxy(pindex) = p%pps%area / g%gps%area
                pfts1d%wtlnd(pindex) = p%pps%area / l%lps%area
                pfts1d%wtcol(pindex) = p%pps%area / c%cps%area
             end do
          end do
       end do
    end do

#if (defined SPMD)

    ! masterproc gridcell gather 

    if (masterproc) then
       allocate(rglob(numg))
       allocate(iglob(numg))
    endif

    rloc => grid1d%wtxy
    call gather_data_to_master(rloc, rglob, clmlevel=grid1d%name)
    if (masterproc) grid1d%wtxy(:) = rglob(:)

    iloc => grid1d%ixy
    call gather_data_to_master(iloc, iglob, clmlevel=grid1d%name)
    if (masterproc) grid1d%ixy(:) = iglob(:)

    iloc => grid1d%jxy
    call gather_data_to_master(iloc, iglob, clmlevel=grid1d%name)
    if (masterproc) grid1d%jxy(:) = iglob(:)

    rloc => grid1d%latdeg
    call gather_data_to_master(rloc, rglob, clmlevel=grid1d%name)
    if (masterproc) grid1d%latdeg(:) = rglob(:)

    rloc => grid1d%londeg
    call gather_data_to_master(rloc, rglob, clmlevel=grid1d%name)
    if (masterproc) grid1d%londeg(:) = rglob(:)

    if (masterproc) then
       deallocate(rglob)
       deallocate(iglob)
    endif

    ! masterproc landunit gather

    if (masterproc) then
       allocate(rglob(numl))
       allocate(iglob(numl))
    endif

    rloc => land1d%wtxy
    call gather_data_to_master(rloc, rglob, clmlevel=land1d%name)
    if (masterproc) land1d%wtxy(:) = rglob(:)

    iloc => land1d%ixy
    call gather_data_to_master(iloc, iglob, clmlevel=land1d%name)
    if (masterproc) land1d%ixy(:) = iglob(:)

    iloc => land1d%jxy
    call gather_data_to_master(iloc, iglob, clmlevel=land1d%name)
    if (masterproc) land1d%jxy(:) = iglob(:)

    iloc => land1d%gindex
    call gather_data_to_master(iloc, iglob, clmlevel=land1d%name)
    if (masterproc) land1d%gindex(:) = iglob(:)

    rloc => land1d%latdeg
    call gather_data_to_master(rloc, rglob, clmlevel=land1d%name)
    if (masterproc) land1d%latdeg(:) = rglob(:)

    rloc => land1d%londeg
    call gather_data_to_master(rloc, rglob, clmlevel=land1d%name)
    if (masterproc) land1d%londeg(:) = rglob(:)

    iloc => land1d%itypwat
    call gather_data_to_master(iloc, iglob, clmlevel=land1d%name)
    if (masterproc) land1d%itypwat(:) = iglob(:)

    if (masterproc) then
       deallocate(rglob)
       deallocate(iglob)
    endif

    ! masterproc column gather

    if (masterproc) then
       allocate(rglob(numc))
       allocate(iglob(numc))
    endif

    rloc => cols1d%wtxy
    call gather_data_to_master(rloc, rglob, clmlevel=cols1d%name)
    if (masterproc) cols1d%wtxy(:) = rglob(:)

    rloc => cols1d%wtlnd
    call gather_data_to_master(rloc, rglob, clmlevel=cols1d%name)
    if (masterproc) cols1d%wtlnd(:) = rglob(:)

    iloc => cols1d%ixy
    call gather_data_to_master(iloc, iglob, clmlevel=cols1d%name)
    if (masterproc) cols1d%ixy(:) = iglob(:)

    iloc => cols1d%jxy
    call gather_data_to_master(iloc, iglob, clmlevel=cols1d%name)
    if (masterproc) cols1d%jxy(:) = iglob(:)

    iloc => cols1d%gindex
    call gather_data_to_master(iloc, iglob, clmlevel=cols1d%name)
    if (masterproc) cols1d%gindex(:) = iglob(:)

    iloc => cols1d%lindex
    call gather_data_to_master(iloc, iglob, clmlevel=cols1d%name)
    if (masterproc) cols1d%lindex(:) = iglob(:)

    rloc => cols1d%latdeg
    call gather_data_to_master(rloc, rglob, clmlevel=cols1d%name)
    if (masterproc) cols1d%latdeg(:) = rglob(:)

    rloc => cols1d%londeg
    call gather_data_to_master(rloc, rglob, clmlevel=cols1d%name)
    if (masterproc) cols1d%londeg(:) = rglob(:)

    iloc => cols1d%itypwat
    call gather_data_to_master(iloc, iglob, clmlevel=cols1d%name)
    if (masterproc) cols1d%itypwat(:) = iglob(:)

    if (masterproc) then
       deallocate(rglob)
       deallocate(iglob)
    endif

    ! masterproc pft gather 

    if (masterproc) then
       allocate(rglob(nump))
       allocate(iglob(nump))
    endif

    rloc => pfts1d%wtxy
    call gather_data_to_master(rloc, rglob, clmlevel=pfts1d%name)
    if (masterproc) pfts1d%wtxy(:) = rglob(:)

    rloc => pfts1d%wtlnd
    call gather_data_to_master(rloc, rglob, clmlevel=pfts1d%name)
    if (masterproc) pfts1d%wtlnd(:) = rglob(:)

    rloc => pfts1d%wtcol
    call gather_data_to_master(rloc, rglob, clmlevel=pfts1d%name)
    if (masterproc) pfts1d%wtcol(:) = rglob(:)

    iloc => pfts1d%ixy
    call gather_data_to_master(iloc, iglob, clmlevel=pfts1d%name)
    if (masterproc) pfts1d%ixy(:) = iglob(:)

    iloc => pfts1d%jxy
    call gather_data_to_master(iloc, iglob, clmlevel=pfts1d%name)
    if (masterproc) pfts1d%jxy(:) = iglob(:)

    iloc => pfts1d%mxy
    call gather_data_to_master(iloc, iglob, clmlevel=pfts1d%name)
    if (masterproc) pfts1d%mxy(:) = iglob(:)

    iloc => pfts1d%gindex
    call gather_data_to_master(iloc, iglob, clmlevel=pfts1d%name)
    if (masterproc) pfts1d%gindex(:) = iglob(:)

    iloc => pfts1d%lindex
    call gather_data_to_master(iloc, iglob, clmlevel=pfts1d%name)
    if (masterproc) pfts1d%lindex(:) = iglob(:)

    iloc => pfts1d%cindex
    call gather_data_to_master(iloc, iglob, clmlevel=pfts1d%name)
    if (masterproc) pfts1d%cindex(:) = iglob(:)

    rloc => pfts1d%latdeg
    call gather_data_to_master(rloc, rglob, clmlevel=pfts1d%name)
    if (masterproc) pfts1d%latdeg(:) = rglob(:)

    rloc => pfts1d%londeg
    call gather_data_to_master(rloc, rglob, clmlevel=pfts1d%name)
    if (masterproc) pfts1d%londeg(:) = rglob(:)

    iloc => pfts1d%itypwat
    call gather_data_to_master(iloc, iglob, clmlevel=pfts1d%name)
    if (masterproc) pfts1d%itypwat(:) = iglob(:)

    iloc => pfts1d%itypveg
    call gather_data_to_master(iloc, iglob, clmlevel=pfts1d%name)
    if (masterproc) pfts1d%itypveg(:) = iglob(:)

    if (masterproc) then
       deallocate(rglob)
       deallocate(iglob)
    endif

#endif

    ! When COMPETE is not defined, it is assumed that each column has 
    ! one pft and indices into the 1d column vector must equal indices 
    ! into the 1d pft vector 

    begc = cols1d%beg
    endc = cols1d%end
    begp = pfts1d%beg
    nump = pfts1d%num
    endp = pfts1d%end
    numc = cols1d%num

    if (numc /= nump) then
       write(6,*)'total number of columns must equal total number of pfts'
       write(6,*)'total number of columns = ',numc
       write(6,*)'total number of pfts = ',nump
       call endrun
    endif
    if (begc /= begp) then
       write(6,*)'beginning index of columns must equal beginning index of pfts'
       write(6,*)'beginning index of columns = ',begc
       write(6,*)'beginning index of pfts = ',begp
       call endrun
    endif
    if (endc /= endp) then
       write(6,*)'ending index of columns must equal ending index of pfts'
       write(6,*)'ending index of columns = ',endc
       write(6,*)'ending index of pfts = ',endp
       call endrun
    endif
    do pi = begp,endp
       p => ppoint%pft(pi)%p
       c => cpoint%col(pi)%c
       if (p%pps%index1d /= c%cps%index1d) then
          write(6,*)'COMPETE not defined: column indices must match pft indices'
          write(6,*)'indexp= ',p%pps%index1d,' indexc= ',c%cps%index1d
          call endrun
       end if
    end do

  end subroutine clm_map1d

end module clm_mapping

