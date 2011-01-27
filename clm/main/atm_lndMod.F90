#include <misc.h>
#include <preproc.h>

module atm_lndMod

#if (defined COUP_CAM)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Atm - Land interface module
! 
! Method: 
! If running as part of cam, the land surface model must use the same 
! grid as the cam. The land surface model calculates its own net solar 
! radiation and net longwave radiation at the surface. The net longwave 
! radiation at the surface will differ somewhat from that calculated in 
! the atmospheric model because the atm model will use the upward 
! longwave flux (or radiative temperature) from the previous time
! step whereas the land surface model uses the flux for the current
! time step. The net solar radiation should equal that calculated
! in the atmospheric model. If not, there is a problem in how the models
! are coupled.
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id$
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid, only: plon, plat
  use rgrid, only: nlon
  use ppgrid, only: begchunk, endchunk
  use phys_grid, only: gather_chunk_to_field
  use comsrf, only :snowhland, srfflx_state2d, srfflx_parm2d,srfflx_parm,surface_state,landfrac
  use history, only :  ctitle, inithist, nhtfrq, mfilt
  use filenames, only: caseid
  use shr_const_mod, only: SHR_CONST_PI
  implicit none
  
  private              ! By default make data private
  integer :: landmask(plon,plat) !2d land mask
  
  public atmlnd_ini, atmlnd_drv ! Public interfaces


!===============================================================================
CONTAINS
!===============================================================================

  subroutine atmlnd_ini(srfflx2d)

    use initializeMod, only : initialize           !initialization of clm 
    use lp_coupling, only: alltoall_clump_to_chunk_init
#if ( defined SPMD )
    use mpishorthand
#endif
    use commap    
    use time_manager, only: get_nstep
    use filenames,    only: mss_irt

#include <comsol.h>
#include <comctl.h>

!-----------------------------------------------------------------------
! Initialize land surface model and obtain relevant atmospheric model 
! arrays back from (i.e. albedos, surface temperature and snow cover over land)
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
    integer  :: i,lat               !indices
    integer  :: nstep               !current timestep number
    integer  :: ier                 !returned error code
    real(r8) :: oro_glob(plon,plat) !global oro field
    real(r8) :: lsmlandfrac(plon,plat) !2d fractional land
    real(r8) :: latixy(plon,plat)   !2d latitude  grid (degrees)
    real(r8) :: longxy(plon,plat)   !2d longitude grid (degrees)
    real(r8) :: pi
    type(srfflx_parm), intent(inout), dimension(begchunk:endchunk) :: srfflx2d
!-----------------------------------------------------------------------
! Time management variables.

    nstep = get_nstep()

! Initialize land model

#if (defined SPMD) && (defined TIMING_BARRIERS)
    call t_startf('sync_gather_oro')
    call mpi_barrier(mpicom, ier)
    if (ier /= 0) then
       write (6,*) 'atmlnd_ini(): MPI error ', ier, ' from mpi_barrier()'
       call endrun
    endif
    call t_stopf('sync_gather_oro')
#endif

    call gather_chunk_to_field(1,1,1,plon,landfrac,oro_glob)
#if (defined SPMD)
    call mpi_bcast (oro_glob, size(oro_glob), mpir8, 0, mpicom, ier)
    if (ier /= 0) then
       write (6,*) 'atmlnd_ini(): MPI error ', ier, ' from mpi_bcast()'
       call endrun
    endif
#endif

    pi = SHR_CONST_PI
    longxy(:,:) = 1.e36
    do lat = 1,plat
       do i = 1,nlon(lat)
          longxy(i,lat) = (i-1)*360.0/nlon(lat)
          latixy(i,lat) = (180./pi)*clat(lat)
          if (oro_glob(i,lat) > 0.) then
             landmask(i,lat) = 1
             lsmlandfrac(i,lat) = oro_glob(i,lat)
          else
             landmask(i,lat) = 0
             lsmlandfrac(i,lat) = 0.
          endif
       end do
    end do

! Initialize albedos, surface temperature, upward longwave radiation,
! and snow depth for land points (used for initial run only)

    call  initialize(eccen    , obliqr   , lambm0  , mvelpp  , caseid  , &
                     ctitle   , nsrest   , nstep   , iradsw  , inithist, &
                     nhtfrq(1), mfilt(1) , longxy  , latixy  , nlon    , &
                     landmask , lsmlandfrac , mss_irt)

! For initial run only - get 2d data back from land model (Note that 
! in SPMD case, only masterproc contains valid recv2d data) and 
! split 2d data into appropriate arrays contained in module comsrf. 

    if (nstep == 0) then
#if (defined SPMD) && (defined TIMING_BARRIERS)
       call t_startf('sync_cl2ch_ini')
       call mpi_barrier(mpicom, ier)
       if (ier /= 0) then
          write (6,*) 'atmlnd_ini(): MPI error ', ier, ' from mpi_barrier()'
          call endrun
       endif
       call t_stopf('sync_cl2ch_ini')
#endif
       call t_startf('clump2chunkini')
       call alltoall_clump_to_chunk_init(srfflx2d)
       call t_stopf('clump2chunkini')
    endif
    
    return
  end subroutine atmlnd_ini

!===============================================================================

  subroutine atmlnd_drv (nstep, iradsw, eccen, obliqr, lambm0, mvelpp,&
                         srf_state,srfflx2d)

!-----------------------------------------------------------------------
! Pack data to be sent to land model into a single array. 
! Send data to land model and call land model driver. 
! Receive data back from land model in a single array.
! Unpack this data into component arrays. 
! NOTE: component arrays are contained in module comsrf.
! When coupling to an atmospheric model: solar radiation depends on 
! surface albedos from the previous time step (based on current
! surface conditions and solar zenith angle for next time step).
! Longwave radiation depends on upward longwave flux from previous
! time step.
!-----------------------------------------------------------------------

#if ( defined SPMD )
    use mpishorthand
#endif
    use comsrf, only:surface_state
    use lp_coupling, only: alltoall_chunk_to_clump, alltoall_clump_to_chunk
!---------------------------Arguments----------------------------------- 
    integer , intent(in) :: nstep    !Current time index
    integer , intent(in) :: iradsw   !Iteration frequency for shortwave radiation
    real(r8), intent(in) :: eccen    !Earth's orbital eccentricity
    real(r8), intent(in) :: obliqr   !Earth's obliquity in radians
    real(r8), intent(in) :: lambm0   !Mean longitude of perihelion at the vernal equinox (radians)
    real(r8), intent(in) :: mvelpp   !Earth's moving vernal equinox longitude of perihelion + pi (radians)
   type(srfflx_parm), intent(inout), dimension(begchunk:endchunk) :: srfflx2d
   type(surface_state), intent(inout), dimension(begchunk:endchunk) :: srf_state
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
    integer :: i,lat,m,n,lchnk,ncols !indices
    integer :: ier         !returned error code
    logical doalb          !true if surface albedo calculation time step
!-----------------------------------------------------------------------

! -----------------------------------------------------------------
! Determine doalb
! [doalb] is a logical variable that is true when the next time
! step is a radiation time step. This allows for the fact that
! an atmospheric model may not do the radiative calculations 
! every time step. For example:
!      nstep dorad doalb
!        1     F     F
!        2     F     T
!        3     T     F
!        4     F     F
!        5     F     T
!        6     T     F
! The following expression for doalb is for example only (it is 
! specific to the NCAR CAM). This variable must be calculated
! appropriately for the host atmospheric model
! -----------------------------------------------------------------

    doalb = iradsw==1 .or. (mod(nstep,iradsw)==0 .and. nstep+1/=1)

#if (defined SPMD) && (defined TIMING_BARRIERS)
   call t_startf('sync_chnk2clmp')
   call mpi_barrier(mpicom, ier)
   if (ier /= 0) then
      write (6,*) 'atmlnd_drv(): MPI error ', ier, ' from mpi_barrier()'
      call endrun
   endif
   call t_stopf('sync_chnk2clmp')
#endif

   call t_startf('chunk2clump')
   call alltoall_chunk_to_clump(srf_state)
   call t_stopf('chunk2clump')

! Call land model driver

    call driver (doalb, eccen, obliqr, lambm0, mvelpp)

#if (defined SPMD) && (defined TIMING_BARRIERS)
   call t_startf('sync_clmp2chnk')
   call mpi_barrier(mpicom, ier)
   if (ier /= 0) then
      write (6,*) 'atmlnd_drv(): MPI error ', ier, ' from mpi_barrier()'
      call endrun
   endif
   call t_stopf('sync_clmp2chnk')
#endif

   call t_startf('clump2chunk')
   call alltoall_clump_to_chunk(srfflx2d)
   call t_stopf('clump2chunk')

    return
  end subroutine atmlnd_drv

!===============================================================================

#endif        

end module atm_lndMod
