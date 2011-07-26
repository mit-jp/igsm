#include <misc.h>
#include <preproc.h>

module lp_coupling

#if (defined COUP_CAM)

!------------------------------------------------------------------------------
!BOP
!
! !MODULE: lp_coupling
!
! !DESCRIPTION:
! Module provides coupling between the atmosphere physics (decomposed into
! chunks) and the land (decomposed into clumps).
!
! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8

#if (defined SPMD)
   use mpishorthand, only: mpir8, mpicom
   use spmd_dyn, only: npes
   use pmgrid, only  : iam
#else
   use spmdMod, only: npes, iam
#endif
   use lnd_grid, only: get_nclumps, get_clump_owner_id, get_clump_ncells_id, &
        get_clump_coord_id, get_clump_gcell_info
   use phys_grid, only: get_chunk_coord_owner_p
!
! !PUBLIC TYPES:
   implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
   public lp_coupling_init             ! initialize clump<-->chunk mapping
   public lp_coupling_finalize         ! destroy clump<-->chunk mapping
   public alltoall_clump_to_chunk_init ! communicate fluxes from lnd to atm
   public alltoall_clump_to_chunk      ! communicate fluxes from lnd to atm
   public alltoall_chunk_to_clump      ! communicate fluxes from atm to lnd
   SAVE
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
! 2002.11.18  Mariana Vertenstein and Forrest Hoffman Updated to clm2.1.
!
!EOP
!
! !PRIVATE TYPES
   private

   type clump2chunk
      integer :: lchunk
      integer :: col
   end type clump2chunk
   type(clump2chunk), dimension(:,:), allocatable, private :: clump2chunks

   type chunk2clump
      integer :: clumpid
      integer :: cell
   end type chunk2clump
   type(chunk2clump), dimension(:,:), allocatable, private :: chunk2clumps

   real(r8), dimension(:), target, allocatable :: lp_sendbuf ! lnd->phys send buf
   real(r8), dimension(:), target, allocatable :: lp_recvbuf ! lnd->phys receive buf
   real(r8), dimension(:), target, allocatable :: pl_sendbuf ! phys->lnd send buf
   real(r8), dimension(:), target, allocatable :: pl_recvbuf ! phys->lnd receive buf

   integer, dimension(:), allocatable :: lp_blkcnts ! l->p send/p->l recv blocks
   integer, dimension(:), allocatable :: lp_sndcnts ! lnd->phys send counts
   integer, dimension(:), allocatable :: lp_rcvcnts ! lnd->phys receive counts
   integer, dimension(:), allocatable :: lp_sdispls ! lnd->phys send dsplsmnt
   integer, dimension(:), allocatable :: lp_rdispls ! lnd->phys receive dsplsmnt

   integer, dimension(:), allocatable :: pl_blkcnts ! p->l send/l->p recv blocks
   integer, dimension(:), allocatable :: pl_sndcnts ! phys->lnd send counts
   integer, dimension(:), allocatable :: pl_rcvcnts ! phys->lnd receive counts
   integer, dimension(:), allocatable :: pl_sdispls ! phys->lnd send dsplsmnt
   integer, dimension(:), allocatable :: pl_rdispls ! phys->lnd receive dsplsmnt

   integer, parameter :: pl_nval = 16        ! phys->lnd flux values
   integer, parameter :: lp_nval = 13        ! lnd->phys flux values

   logical :: lpc_init_flag = .false.        ! set if initialized
!------------------------------------------------------------------------------
! $Id$
! $Author$
!------------------------------------------------------------------------------

   contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lp_coupling_init
!
! !INTERFACE:
   subroutine lp_coupling_init()
!
! !DESCRIPTION:
! This subroutine initializes the mapping between the atmosphere physics chunks
! and the land clumps.  It may (and must) be called repeatedly to re-initialize 
! the mapping if the decomposition of either the atmosphere physics or the land
! changes.  It allocates communication buffers constructs vectors of counts and
! displacements used for subsequent communication between MPI processes.
!
! !ARGUMENTS:
   implicit none
!
! !LOCAL VARIABLES:
   integer :: p, c, g                            ! loop indices
   integer :: nclumps                            ! number of clumps defined
   integer :: ncells                             ! number of clump cells
   integer :: clump_owner                        ! clump owner
   integer, dimension(:), allocatable :: lons    ! clump longitudes
   integer, dimension(:), allocatable :: lats    ! clump latitudes
   integer, dimension(:), allocatable :: lchnks  ! chunk ids
   integer, dimension(:), allocatable :: cols    ! chunk columns
   integer, dimension(:), allocatable :: chunk_owners  ! chunk owners
   integer :: max_gpc = 0                        ! max cells per clump
   integer :: ier                                ! error codes
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
! 2002.11.18  Mariana Vertenstein and Forrest Hoffman Updated to clm2.1.
!
!EOP
!------------------------------------------------------------------------------
! $Id$
! $Author$
!------------------------------------------------------------------------------

   ! If already initialized, then deallocate buffers and re-initialize everything

   call lp_coupling_finalize()

   allocate(lp_blkcnts(0:npes-1), lp_sndcnts(0:npes-1), lp_rcvcnts(0:npes-1), &
            pl_blkcnts(0:npes-1), pl_sndcnts(0:npes-1), pl_rcvcnts(0:npes-1), &
            stat=ier)
   if (ier /= 0) then
      write (6,*) 'lp_coupling_init(): allocation error for lp_blkcnts, ', &
         'lp_sndcnts, pl_blkcnts, and pl_sndcnts'
      call endrun
   end if
   lp_blkcnts(:) = 0
   lp_sndcnts(:) = 0
   lp_rcvcnts(:) = 0
   pl_blkcnts(:) = 0
   pl_sndcnts(:) = 0
   pl_rcvcnts(:) = 0

   ! Determine max_gpc and allocate dynamic memory

   nclumps = get_nclumps()
   do c = 1,nclumps
      ncells = get_clump_ncells_id(c)
      if (ncells > max_gpc) max_gpc = ncells
   end do
   allocate(lons(max_gpc), lats(max_gpc), lchnks(max_gpc), cols(max_gpc), &
        chunk_owners(max_gpc), stat=ier)
   if (ier /= 0) then
      write (6,*) 'lp_coupling_init(): allocation error for local lons, ', &
         'lats, lchnks, cols, and chunk_owners variables'
      call endrun
   end if

   ! Found above already
   ! nclumps = get_nclumps()
   ! lp_blkcnts: for each cam pid, determine the total number of sends that
   ! will be sent from clm
   ! pl_blkcnts: for each clm pid, determine the total number of sends that
   ! will be sent from cam

   do c = 1,nclumps
      clump_owner = get_clump_owner_id(c)
      ncells = get_clump_ncells_id(c)
      call get_clump_coord_id(c, ncells, lons, lats)
      call get_chunk_coord_owner_p(ncells, lons, lats, lchnks, cols, chunk_owners)
      do g = 1,ncells
         if (clump_owner == iam) then
            p = chunk_owners(g)
            lp_blkcnts(p) = lp_blkcnts(p) + 1
         endif
         if (chunk_owners(g) == iam) then
            p = clump_owner
            pl_blkcnts(p) = pl_blkcnts(p) + 1
         endif
      end do
   end do

   allocate(clump2chunks(0:npes-1, 1:maxval(pl_blkcnts)), &
            chunk2clumps(0:npes-1, 1:maxval(lp_blkcnts)), stat=ier)
   if (ier /= 0) then
      write (6,*) 'lp_coupling_init(): allocation error for clump2chunks ', &
         'and chunk2clumps'
      call endrun
   end if
   clump2chunks(:,:)%lchunk = 0
   clump2chunks(:,:)%col = 0
   chunk2clumps(:,:)%clumpid = 0
   chunk2clumps(:,:)%cell = 0


   ! Found above already
   ! nclumps = get_nclumps()

   lp_blkcnts(:) = 0
   pl_blkcnts(:) = 0
   do c = 1,nclumps
      clump_owner = get_clump_owner_id(c)
      ncells = get_clump_ncells_id(c)
      call get_clump_coord_id(c, ncells, lons, lats)
      call get_chunk_coord_owner_p(ncells,lons,lats,lchnks,cols,chunk_owners)
      do g = 1,ncells
         if (clump_owner == iam) then
            p = chunk_owners(g)
            lp_blkcnts(p) = lp_blkcnts(p) + 1
            chunk2clumps(p,lp_blkcnts(p))%clumpid=c
            chunk2clumps(p,lp_blkcnts(p))%cell = g
         end if
         if (chunk_owners(g) == iam) then
            p = clump_owner
            pl_blkcnts(p) = pl_blkcnts(p) + 1
            clump2chunks(p,pl_blkcnts(p))%lchunk = lchnks(g)
            clump2chunks(p,pl_blkcnts(p))%col = cols(g)
         end if
      end do
   end do

   deallocate(lons, lats, lchnks, cols, chunk_owners)

   pl_sndcnts(:) = pl_blkcnts(:) * pl_nval
   lp_rcvcnts(:) = pl_blkcnts(:) * lp_nval
   lp_sndcnts(:) = lp_blkcnts(:) * lp_nval
   pl_rcvcnts(:) = lp_blkcnts(:) * pl_nval

   allocate(pl_sendbuf(0:sum(pl_sndcnts)-1), lp_recvbuf(0:sum(lp_rcvcnts)-1), &
      stat=ier)
   if (ier /= 0) then
      write (6,*) 'lp_coupling_init(): allocation error for pl_sendbuf and ', &
         'lp_recvbuf'
      call endrun
   end if
   allocate(lp_sendbuf(0:sum(lp_sndcnts)-1), pl_recvbuf(0:sum(pl_rcvcnts)-1), &
      stat=ier)
   if (ier /= 0) then
      write (6,*) 'lp_coupling_init(): allocation error for lp_sendbuf and ', &
         'pl_recvbuf'
      call endrun
   end if

   allocate(lp_sdispls(0:npes-1), lp_rdispls(0:npes-1), pl_sdispls(0:npes-1), &
      pl_rdispls(0:npes-1), stat=ier)
   if (ier /= 0) then
      write (6,*) 'lp_coupling_init(): allocation error for lp_sdispls, ', &
         'lp_rdispls, pl_sdispls, and pl_rdispls'
      call endrun
   end if

   lp_sdispls(0) = 0
   lp_rdispls(0) = 0
   pl_sdispls(0) = 0
   pl_rdispls(0) = 0
   do p = 1,npes-1
      lp_sdispls(p) = lp_sdispls(p-1) + lp_sndcnts(p-1)
      lp_rdispls(p) = lp_rdispls(p-1) + lp_rcvcnts(p-1)
      pl_sdispls(p) = pl_sdispls(p-1) + pl_sndcnts(p-1)
      pl_rdispls(p) = pl_rdispls(p-1) + pl_rcvcnts(p-1)
   end do

   lpc_init_flag = .true.

   end subroutine lp_coupling_init

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lp_coupling_finalize
!
! !INTERFACE:
   subroutine lp_coupling_finalize()
!
! !ARGUMENTS:
   implicit none
!
! !DESCRIPTION:
! This subroutine destroys the mapping between the atmsphere physics chunks and
! the land clumps if the lpc\_init\_flag flag is set.  It is called from
! lp\_coupling\_init() to ensure memory is recycled when a new mapping is to be
! created.
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
! 2002.11.18  Mariana Vertenstein and Forrest Hoffman Updated to clm2.1.
!
!EOP
!------------------------------------------------------------------------------
! $Id$
! $Author$
!------------------------------------------------------------------------------

   if (lpc_init_flag) then
      deallocate(clump2chunks, chunk2clumps)
      deallocate(lp_sendbuf, lp_recvbuf, pl_sendbuf, pl_recvbuf)
      deallocate(lp_blkcnts, pl_blkcnts)
      deallocate(lp_sndcnts, lp_rcvcnts, pl_sndcnts, pl_rcvcnts)
      deallocate(lp_sdispls, lp_rdispls, pl_sdispls, pl_rdispls)
      lpc_init_flag = .false.
   endif

   end subroutine lp_coupling_finalize

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: alltoall_clump_to_chunk_init
!
! !INTERFACE:
   subroutine alltoall_clump_to_chunk_init(srfflx2d)
!
! !DESCRIPTION:
! This subroutine performs the initial communication from the land model to the
! atmosphere physics (from clumps to chunks) based on the mapping constructed
! in lp\_coupling\_init().
!
! !USES:
   use ppgrid, only: begchunk, endchunk
   use comsrf, only: srfflx_parm, snowhland
   use clm_varcon, only: sb
   use clmtype
   use clmpoint, only : gpoint
!
! !ARGUMENTS:
   implicit none
   type(srfflx_parm), intent(inout), dimension(begchunk:endchunk) :: srfflx2d
!
! !LOCAL VARIABLES:
   integer  :: p, n, k, m                          ! loop indices
   integer  :: begg, numg                          ! gridcell beginning index and count
   integer  :: lchnk, i                            ! local chunk and column
   real(r8), dimension(:), pointer :: lp_rbufp     ! recv buffer pointer
   integer  :: ier                                 ! returned error code
   integer  :: is,gi,li,ci,pi                  ! indices
   integer  :: indexg                              ! 1d gridcell index 
   type(gridcell_type)       , pointer :: g        ! local pointer to derived subtype
   type(landunit_type)       , pointer :: l        ! local pointer to derived subtype
   type(column_type)         , pointer :: c        ! local pointer to derived subtype
   type(gridcell_pstate_type), pointer :: gps      ! local pointer to derived subtype
   type(landunit_pstate_type), pointer :: lps      ! local pointer to derived subtype
   type(column_pstate_type)  , pointer :: cps      ! local pointer to derived subtype
   type(column_estate_type)  , pointer :: ces      ! local pointer to derived subtype 
   type(column_wstate_type)  , pointer :: cws      ! local pointer to derived subtype
   type(pft_pstate_type)     , pointer :: pps      ! local pointer to derived subtype
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
! 2002.11.18  Mariana Vertenstein and Forrest Hoffman Updated to clm2.1.
!
!EOP
!------------------------------------------------------------------------------
! $Id$
! $Author$
!------------------------------------------------------------------------------

   ! Fill lnd->phys send buffer

   lp_sendbuf(:) = 0.0_r8

!$OMP PARALLEL DO PRIVATE(p,n,gi,g,gps,li,l,lps,ci,c,cps,cws,ces,is,pi,pps)
   do p = 0, npes-1
      do n = 1, lp_blkcnts(p)

         ! Determine clump 1d gridcell index
         call get_clump_gcell_info(chunk2clumps(p,n)%clumpid, &
              chunk2clumps(p,n)%cell, gi)

         ! Loop over all grid cells in clump
         g => gpoint%grd(gi)%g
         gps => g%gps
         do li = 1, gps%nlandunits
            l => g%l(li)
            lps => l%lps
            do ci = 1, lps%ncolumns
               c => l%c(ci)
               cps => c%cps
               cws => c%cws
               ces => c%ces
               is = lp_sdispls(p)+(n-1)*lp_nval + 0
               lp_sendbuf(is) = lp_sendbuf(is) + ces%t_grnd * cps%wtxy       ! tsxy
               is = lp_sdispls(p)+(n-1)*lp_nval + 5
               lp_sendbuf(is) = lp_sendbuf(is) + cws%h2osno/1000. * cps%wtxy ! snow [mm]->[m]
               is = lp_sdispls(p)+(n-1)*lp_nval + 10
               lp_sendbuf(is) = lp_sendbuf(is) + sb*(ces%t_grnd**4) * cps%wtxy 
               do pi = 1, cps%npfts
                  pps => c%p(pi)%pps
                  is = lp_sdispls(p)+(n-1)*lp_nval + 1
                  lp_sendbuf(is) = lp_sendbuf(is) + pps%albd(1) * pps%wtxy   ! asdir 
                  is = lp_sdispls(p)+(n-1)*lp_nval + 2
                  lp_sendbuf(is) = lp_sendbuf(is) + pps%albd(2) * pps%wtxy   ! aldir
                  is = lp_sdispls(p)+(n-1)*lp_nval + 3
                  lp_sendbuf(is) = lp_sendbuf(is) + pps%albi(1) * pps%wtxy   ! asdif
                  is = lp_sdispls(p)+(n-1)*lp_nval + 4
                  lp_sendbuf(is) = lp_sendbuf(is) + pps%albi(2) * pps%wtxy   ! aldif
               end do
            end do
         end do
         is = lp_sdispls(p)+(n-1)*lp_nval + 6
         lp_sendbuf(is) = 1.e36 
         is = lp_sdispls(p)+(n-1)*lp_nval + 7
         lp_sendbuf(is) = 1.e36 
         is = lp_sdispls(p)+(n-1)*lp_nval + 8
         lp_sendbuf(is) = 1.e36 
         is = lp_sdispls(p)+(n-1)*lp_nval + 9
         lp_sendbuf(is) = 1.e36 
         is = lp_sdispls(p)+(n-1)*lp_nval + 11
         lp_sendbuf(is) = 1.e36 
         is = lp_sdispls(p)+(n-1)*lp_nval + 12
         lp_sendbuf(is) = 1.e36 
      end do
   end do

#ifdef SPMD
   call mpi_alltoallv (lp_sendbuf, lp_sndcnts, lp_sdispls, mpir8, &
                       lp_recvbuf, lp_rcvcnts, lp_rdispls, mpir8, mpicom, ier)
   if (ier /= 0) then
      write (6,*) 'alltoall_clump_to_chunk_init(): MPI error ', ier, &
         ' from mpi_alltoallv()'
      call endrun
   endif
   lp_rbufp => lp_recvbuf
#else
   lp_rbufp => lp_sendbuf
#endif

   ! Extract lnd->phys receive buffer

!$OMP PARALLEL DO PRIVATE(p,n,lchnk,i)
   do p = 0,npes-1
      do n = 1,pl_blkcnts(p)
         lchnk = clump2chunks(p,n)%lchunk
         i = clump2chunks(p,n)%col
         srfflx2d(lchnk)%ts(i)    = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+0)
         srfflx2d(lchnk)%asdir(i) = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+1)
         srfflx2d(lchnk)%aldir(i) = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+2)
         srfflx2d(lchnk)%asdif(i) = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+3)
         srfflx2d(lchnk)%aldif(i) = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+4)
         snowhland(i,lchnk)       = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+5)
         srfflx2d(lchnk)%lwup(i)  = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+10)
      end do
   end do

   return
   end subroutine alltoall_clump_to_chunk_init

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: alltoall_clump_to_chunk
!
! !INTERFACE:
   subroutine alltoall_clump_to_chunk(srfflx2d)
!
! !DESCRIPTION:
! This subroutine performs the communication from the land model to the
! atmosphere physics (from clumps to chunks) based on the mapping constructed
! in lp\_coupling\_init().
!
! !USES:
   use ppgrid, only: begchunk, endchunk
   use comsrf, only: srfflx_parm, snowhland
   use constituents, only: pcnst, pnats
   use clmtype
   use clmpoint, only : gpoint
!
! !ARGUMENTS:
   implicit none
   type(srfflx_parm), intent(inout), dimension(begchunk:endchunk) :: srfflx2d
!
! !LOCAL VARIABLES:
   integer  :: p, n, k, m                          ! loop indices
   integer  :: bpatch, npatch                      ! patch id and count
   integer  :: lchnk, i                            ! local chunk and column
   real(r8) :: wt                                  ! patch wt
   real(r8), dimension(:), pointer :: lp_rbufp     ! recv buffer pointer
   integer  :: ier                                 ! returned error code
   integer  :: is                                  ! send buffer index
   integer  :: gi,li,ci,pi                         ! clmtype 1d indices
   integer  :: indexg                              ! clmtype gridcell 1d index 
   type(gridcell_type)       , pointer :: g        ! local pointer to derived subtype
   type(landunit_type)       , pointer :: l        ! local pointer to derived subtype
   type(column_type)         , pointer :: c        ! local pointer to derived subtype
   type(gridcell_pstate_type), pointer :: gps      ! local pointer to derived subtype
   type(landunit_pstate_type), pointer :: lps      ! local pointer to derived subtype
   type(column_pstate_type)  , pointer :: cps      ! local pointer to derived subtype
   type(column_wstate_type)  , pointer :: cws      ! local pointer to derived subtype
   type(pft_pstate_type)     , pointer :: pps      ! local pointer to derived subtype
   type(pft_mflux_type)      , pointer :: pmf      ! local pointer to derived subtype
   type(pft_wflux_type)      , pointer :: pwf      ! local pointer to derived subtype
   type(pft_eflux_type)      , pointer :: pef      ! local pointer to derived subtype 
   type(pft_estate_type)     , pointer :: pes      ! local pointer to derived subtype
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
! 2002.11.18  Mariana Vertenstein and Forrest Hoffman Updated to clm2.1.
!
!EOP
!------------------------------------------------------------------------------
! $Id$
! $Author$
!------------------------------------------------------------------------------

   ! Fill lnd->phys send buffer

   lp_sendbuf(:) = 0.0_r8

!$OMP PARALLEL DO PRIVATE(p,n,gi,g,gps,li,l,lps,ci,c,cps,cws,is,pi,pps,pwf,pef,pes,pmf)
   do p = 0,npes-1
      do n = 1,lp_blkcnts(p)

         ! Determine clump gridcell info
         call get_clump_gcell_info(chunk2clumps(p,n)%clumpid, &
              chunk2clumps(p,n)%cell, gi)

         ! Loop over all grid cells in clump
         g => gpoint%grd(gi)%g
         gps => g%gps
         do li = 1, gps%nlandunits
            l => g%l(li)
            lps => l%lps
            do ci = 1, lps%ncolumns
               c => l%c(ci)
               cps => c%cps
               cws => c%cws
               is = lp_sdispls(p)+(n-1)*lp_nval+5   ! snow [mm]->[m]
               lp_sendbuf(is) = lp_sendbuf(is) +  cws%h2osno/1000. * cps%wtxy
               do pi = 1, cps%npfts
                  pps => c%p(pi)%pps
                  pwf => c%p(pi)%pwf
                  pef => c%p(pi)%pef
                  pes => c%p(pi)%pes
                  pmf => c%p(pi)%pmf
                  is = lp_sdispls(p)+(n-1)*lp_nval+0   ! tsxy
                  lp_sendbuf(is) = lp_sendbuf(is) + pes%t_rad_pft * pps%wtxy 
                  is = lp_sdispls(p)+(n-1)*lp_nval+1   ! asdir
                  lp_sendbuf(is) = lp_sendbuf(is) + pps%albd(1) * pps%wtxy  
                  is = lp_sdispls(p)+(n-1)*lp_nval+2   ! aldir
                  lp_sendbuf(is) = lp_sendbuf(is) + pps%albd(2) * pps%wtxy 
                  is = lp_sdispls(p)+(n-1)*lp_nval+3   ! asdif
                  lp_sendbuf(is) = lp_sendbuf(is) + pps%albi(1) * pps%wtxy 
                  is = lp_sdispls(p)+(n-1)*lp_nval+4   ! aldif
                  lp_sendbuf(is) = lp_sendbuf(is) + pps%albi(2) * pps%wtxy 
                  is = lp_sdispls(p)+(n-1)*lp_nval+6   ! taux
                  lp_sendbuf(is) = lp_sendbuf(is) + pmf%taux * pps%wtxy 
                  is = lp_sdispls(p)+(n-1)*lp_nval+7   ! tauy
                  lp_sendbuf(is) = lp_sendbuf(is) + pmf%tauy * pps%wtxy 
                  is = lp_sdispls(p)+(n-1)*lp_nval+8   ! lhflx
                  lp_sendbuf(is) = lp_sendbuf(is) + pef%eflx_lh_tot * pps%wtxy 
                  is = lp_sdispls(p)+(n-1)*lp_nval+9   ! shflx
                  lp_sendbuf(is) = lp_sendbuf(is) + pef%eflx_sh_tot * pps%wtxy 
                  is = lp_sdispls(p)+(n-1)*lp_nval+10  ! lwrad
                  lp_sendbuf(is) = lp_sendbuf(is) + pef%eflx_lwrad_out * pps%wtxy 
                  is = lp_sdispls(p)+(n-1)*lp_nval+11  ! qflx
                  lp_sendbuf(is) = lp_sendbuf(is) + pwf%qflx_evap_tot * pps%wtxy 
                  is = lp_sdispls(p)+(n-1)*lp_nval+12  ! tref
                  lp_sendbuf(is) = lp_sendbuf(is) + pes%t_ref2m * pps%wtxy 
               end do   ! end pft loop   
            end do   ! end column loop
         end do   ! end landunit loop
      end do   ! end loop over destination processes
   end do   ! end loop over source processes

#ifdef SPMD
   call mpi_alltoallv (lp_sendbuf, lp_sndcnts, lp_sdispls, mpir8, &
                       lp_recvbuf, lp_rcvcnts, lp_rdispls, mpir8, mpicom, ier)
   if (ier /= 0) then
      write (6,*) 'alltoall_clump_to_chunk(): MPI error ', ier, &
         ' from mpi_alltoallv()'
      call endrun
   endif
   lp_rbufp => lp_recvbuf
#else
   lp_rbufp => lp_sendbuf
#endif

   ! Extract lnd->phys receive buffer

!$OMP PARALLEL DO PRIVATE(p,n,lchnk,i,m)
   do p = 0,npes-1
      do n = 1,pl_blkcnts(p)
         lchnk = clump2chunks(p,n)%lchunk
         i = clump2chunks(p,n)%col
         srfflx2d(lchnk)%ts(i)     = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+0)
         srfflx2d(lchnk)%asdir(i)  = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+1)
         srfflx2d(lchnk)%aldir(i)  = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+2)
         srfflx2d(lchnk)%asdif(i)  = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+3)
         srfflx2d(lchnk)%aldif(i)  = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+4)
         snowhland(i,lchnk)        = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+5)
         srfflx2d(lchnk)%wsx(i)    = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+6)
         srfflx2d(lchnk)%wsy(i)    = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+7)
         srfflx2d(lchnk)%lhf(i)    = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+8)
         srfflx2d(lchnk)%shf(i)    = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+9)
         srfflx2d(lchnk)%lwup(i)   = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+10)
         srfflx2d(lchnk)%cflx(i,1) = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+11)

         ! Reset all other constituent surface fluxes to zero over land
         do m = 2,pcnst+pnats
            srfflx2d(lchnk)%cflx(i,m) = 0.0_r8
         end do
         srfflx2d(lchnk)%tref(i)  = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+12)
      end do
   end do

   return
   end subroutine alltoall_clump_to_chunk

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: alltoall_chunk_to_clump
!
! !INTERFACE:
   subroutine alltoall_chunk_to_clump(srf_state)
!
! !DESCRIPTION:
! This subroutine performs communication from the atmosphere physics to the
! land model (from chunks to clumps) based on the mapping constructed in
! lp\_coupling\_init().
!
! !USES:
   use ppgrid, only: begchunk, endchunk
   use comsrf, only: surface_state
   use clmtype
   use clmpoint, only : gpoint
   use clm_varcon, only: rair, po2, pco2
!
! !ARGUMENTS:
   implicit none
   type(surface_state), intent(in), dimension(begchunk:endchunk) :: srf_state
!
! !LOCAL VARIABLES:
   integer  :: p, n, k                             ! loop indices
   integer  :: bpatch, npatch                      ! patch id and count
   integer  :: lchnk, i                            ! local chunk and column
   real(r8) :: forc_rainc, forc_rainl              ! rainxy [mm/s]
   real(r8) :: forc_snowc, forc_snowl              ! snowfxy [mm/s]
   real(r8), dimension(:), pointer :: pl_rbufp     ! recv buffer pointer
   integer  :: ier                                 ! returned error code
   integer  :: gi                                  ! 1d grid index
   type(gridcell_type)     , pointer :: g          ! local pointer to derived subtype
   type(atm2lnd_state_type), pointer :: a2ls       ! local pointer to derived subtype
   type(atm2lnd_flux_type) , pointer :: a2lf       ! local pointer to derived subtype
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
! 2002.11.18  Mariana Vertenstein and Forrest Hoffman Updated to clm2.1.
!
!EOP
!------------------------------------------------------------------------------
! $Id$
! $Author$
!------------------------------------------------------------------------------

   ! Fill phys->lnd send buffer

!$OMP PARALLEL DO PRIVATE(p,n,lchnk,i)
   do p = 0,npes-1
      do n = 1,pl_blkcnts(p)
         lchnk = clump2chunks(p,n)%lchunk
         i = clump2chunks(p,n)%col

         ! Atmoshperic state variable [m]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+0) = srf_state(lchnk)%zbot(i)

         ! Atmoshperic state variable [m/s]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+1) = srf_state(lchnk)%ubot(i)

         ! Atmoshperic state variable [m/s]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+2) = srf_state(lchnk)%vbot(i)

         ! Atmoshperic state variable [K]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+3) = srf_state(lchnk)%thbot(i)

         ! Atmoshperic state variable [kg/kg]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+4) = srf_state(lchnk)%qbot(i)

         ! Atmoshperic state variable [Pa]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+5) = srf_state(lchnk)%pbot(i)

         ! Atmoshperic state variable [K]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+6) = srf_state(lchnk)%tbot(i)

         ! Atmoshperic flux [W/m^2]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+7) = srf_state(lchnk)%flwds(i)

         ! Atmoshperic flux [m/s] -> [mm/s]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+8) = &
              srf_state(lchnk)%precsc(i) * 1000.

         ! Atmoshperic flux [m/s] -> [mm/s]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+9) = &
            srf_state(lchnk)%precsl(i) * 1000.

         ! Atmoshperic flux [m/s] -> [mm/s]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+10) = &
            (srf_state(lchnk)%precc(i) - srf_state(lchnk)%precsc(i)) * 1000.

         ! Atmoshperic flux [m/s] -> [mm/s]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+11) = &
            (srf_state(lchnk)%precl(i) - srf_state(lchnk)%precsl(i)) * 1000.

         ! Atmoshperic flux [W/m^2]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+12) = srf_state(lchnk)%soll(i)

         ! Atmoshperic flux [W/m^2]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+13) = srf_state(lchnk)%sols(i)

         ! Atmoshperic flux [W/m^2]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+14) = srf_state(lchnk)%solld(i)

         ! Atmoshperic flux [W/m^2]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+15) = srf_state(lchnk)%solsd(i)

      end do
   end do

#ifdef SPMD
   call mpi_alltoallv (pl_sendbuf, pl_sndcnts, pl_sdispls, mpir8, &
                       pl_recvbuf, pl_rcvcnts, pl_rdispls, mpir8, mpicom, ier)
   if (ier /= 0) then
      write (6,*) 'alltoall_chunk_to_clump(): MPI error ', ier, &
         ' from mpi_alltoallv()'
      call endrun
   endif
   pl_rbufp => pl_recvbuf
#else
   pl_rbufp => pl_sendbuf
#endif

   ! Extract phys->lnd receive buffer

!$OMP PARALLEL DO PRIVATE(p,n,gi,g,a2ls,a2lf,forc_rainc,forc_rainl,forc_snowc,forc_snowl)
   do p = 0,npes-1
      do n = 1,lp_blkcnts(p)

         ! Determine clump gridcell info
         call get_clump_gcell_info(chunk2clumps(p,n)%clumpid, &
              chunk2clumps(p,n)%cell, gi)

         ! Loop over all grid cells in clump
         g => gpoint%grd(gi)%g
         a2ls => g%a2ls
         a2lf => g%a2lf

         a2ls%forc_hgt      = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+0)
         a2ls%forc_u        = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+1)
         a2ls%forc_v        = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+2)
         a2ls%forc_th       = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+3)
         a2ls%forc_q        = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+4)
         a2ls%forc_pbot     = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+5)
         a2ls%forc_t        = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+6)
         a2lf%forc_lwrad    = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+7)
         forc_snowc         = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+8)
         forc_snowl         = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+9)
         forc_rainc         = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+10)
         forc_rainl         = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+11)
         a2lf%forc_solad(2) = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+12)
         a2lf%forc_solad(1) = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+13)
         a2lf%forc_solai(2) = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+14)
         a2lf%forc_solai(1) = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+15)

         ! Determine derived quantities
         
         a2ls%forc_hgt_u = a2ls%forc_hgt !obs height of wind [m]
         a2ls%forc_hgt_t = a2ls%forc_hgt !obs height of temperature [m]
         a2ls%forc_hgt_q = a2ls%forc_hgt !obs height of humidity [m]
         a2ls%forc_vp    = a2ls%forc_q*a2ls%forc_pbot / (0.622+0.378*a2ls%forc_q)
         a2ls%forc_rho   = (a2ls%forc_pbot-0.378*a2ls%forc_vp) / (rair*a2ls%forc_t)
         a2ls%forc_co2   = pco2*a2ls%forc_pbot
         a2ls%forc_o2    = po2*a2ls%forc_pbot
         a2ls%forc_wind  = sqrt(a2ls%forc_u**2 + a2ls%forc_v**2)
         a2lf%forc_solar = a2lf%forc_solad(1) + a2lf%forc_solai(1) + &
                           a2lf%forc_solad(2) + a2lf%forc_solai(2)

         ! Determine precipitation needed by clm

#ifdef PERGRO
         ! For error growth only, allow rain, not snowfall
         forc_rainc           = forc_rainc + forc_snowc
         forc_rainl           = forc_rainl + forc_snowl
         forc_snowc           = 0.0_r8
         forc_snowl           = 0.0_r8
#endif
         a2lf%forc_rain = forc_rainc + forc_rainl
         a2lf%forc_snow = forc_snowc + forc_snowl
         
         if (a2lf%forc_snow > 0.0_r8  .and. a2lf%forc_rain > 0.0_r8) then
            write(6,*)' ERROR: clm2 cannot currently handle both non-zero rain and snow'
            write(6,*)' grid point at i = ',grid1d%ixy(gi),'and j= ',grid1d%jxy(gi),' has ' 
            write(6,*)' snow= ',a2lf%forc_snow,' rain= ',a2lf%forc_rain
            call endrun
         elseif (a2lf%forc_rain > 0.) then
            a2ls%itypprc = 1
         elseif (a2lf%forc_snow > 0.) then
            a2ls%itypprc = 2
         else
            a2ls%itypprc = 0
         endif
         
      end do   ! end loop over destination processes
   end do   ! end loop over source processes

   end subroutine alltoall_chunk_to_clump

#endif

end module lp_coupling
