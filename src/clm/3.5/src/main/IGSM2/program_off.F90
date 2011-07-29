#include <misc.h>
#include <preproc.h>

#if (defined OFFLINE)

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: program_off
!
! !INTERFACE:
! !INTERFACE:
#if  (defined COUP_MIT2D)
subroutine clm4mit2d
#else
PROGRAM program_off
#endif
!
! !DESCRIPTION:
! "off-line" code to mimic coupling to an atmospheric model.
! This program is an "off-line" driver for clm3.
! This code can be used to run the clm3 uncoupled from any atmospheric model.
! The appropriate atmospheric forcing is provided in module [atmdrvMod.F90]
! o If running as an offline driver, the land surface model may use
!   a different grid than the input atmospheric data. The atmospheric
!   data is then interpolated to the land model grid inside the
!   atmospheric driver module [atmdrvMod.F90].
! o If running as part of cam, the land surface model must use the
!   same grid as the cam.
! o If running through the flux coupler, the land surface model grid
!   is interpolated to the atmospheric grid inside the flux coupler
! o To map from the atmospheric grid to the land grid, the atmospheric
!   model must provide latitudes and longitudes (degrees) for each grid
!   point and the North, East, South, and West edges of atmospheric grid.
!   Comparable data for the land grid are provided by the land model.
!   When mapping from land to atm grid, an atm grid cell that is part
!   land and part ocean (as defined by the land surface grid) will have
!   fluxes only based on the land portion.
! o The zenith angle calculation is for the NEXT time step rather
!   than the current time step. Make sure the calendar day is for
!   the NEXT time step. Make sure the calendar day is for Greenwich
!   time (see next comment).
! o The land surface model calculates its own net solar radiation and
!   net longwave radiation at the surface. The net longwave radiation
!   at the surface will differ somewhat from that calculated in the
!   atmospheric model because the atm model will use the upward
!   longwave flux (or radiative temperature) from the previous time
!   step whereas the land surface model uses the flux for the current
!   time step. The net solar radiation should equal that calculated
!   in the atmospheric model. If not, there is a problem in how
!   the models are coupled.
!
! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8, SHR_KIND_CL
  use shr_orb_mod          
  use clm_varorb      , only : eccen, mvelpp, lambm0, obliqr, obliq, &
                               iyear_AD, nmvelp
  use clm_comp        , only : clm_init0, clm_init1, clm_init2, clm_run1, clm_run2
  use clm_time_manager, only : is_last_step, advance_timestep, get_nstep, get_curr_date
  use atmdrvMod       , only : atmdrv, atmdrv_init
  use abortutils      , only : endrun
  use fileutils  , only : relavu, opnfil, getfil, getavu
  use controlMod      , only : control_setNL
  use clm_mct_mod
  use spmdMod  
  use ESMF_Mod
  use perf_mod
!
! !ARGUMENTS:
  implicit none
#if (defined COUP_MIT2D)
  logical,save :: first !First step logical
  data first/.true./
  integer :: yr, mon, day, mcsec
#endif
  include 'IGSM2.inc'
#if (defined COUP_TEM)
  include 'TEM.inc'
#endif
!
! !REVISION HISTORY:
! Author: Gordon Bonan and Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
  integer  :: nstep     ! time step index
  real(r8) :: dtime     ! time step increment (sec)
  integer  :: ier       ! error code

! Orbital information after call to routine shr_orbit_params

  logical  :: log_print    ! true=> print diagnostics
  real(r8) :: eccf         ! earth orbit eccentricity factor
  logical  :: clmmpi_running  ! true => CLMMPI is initialized 
  integer  :: clmmpicom_glob  ! CLMMPI communicator

  character(len=SHR_KIND_CL) :: nlfilename = "lnd.stdin"
!-----------------------------------------------------------------------

  ! -----------------------------------------------------------------
  ! Initialize CLMMPI
  ! -----------------------------------------------------------------
#if (defined COUP_MIT2D)
  if (first) then
#endif

  call clmmpi_initialized (clmmpi_running, ier)
  if (.not. clmmpi_running) call clmmpi_init(ier)
  clmmpicom_glob = CLMMPI_COMM_WORLD
  call spmd_init(clmmpicom_glob)
  call mct_world_init(1,clmmpicom_glob,clmmpicom,comp_id)

  call t_startf('init')

  ! -----------------------------------------------------------------
  ! Initialize ESMF (needed for time-manager)
  ! -----------------------------------------------------------------

  call ESMF_Initialize()

  ! -----------------------------------------------------------------
  ! Initialize timing library, and set full path to namelist
  ! -----------------------------------------------------------------

  call control_setNL( nlfilename )     ! Set namelist
  call t_initf(nlfilename, LogPrint=masterproc, CLMMpicom=clmmpicom, &
               MasterTask=masterproc)

  ! -----------------------------------------------------------------
  ! Initialize Orbital parameters
  ! -----------------------------------------------------------------

  ! obliq, eccen and nmvelp are determined based on value of iyear_AD

  if (masterproc) then
     log_print = .true.
  else
     log_print = .false.
  end if
  iyear_AD = 1950
  obliq    = SHR_ORB_UNDEF_REAL
  eccen    = SHR_ORB_UNDEF_REAL
  nmvelp   = SHR_ORB_UNDEF_REAL
  call shr_orb_params (iyear_AD, eccen, obliq, nmvelp, obliqr, &
                       lambm0, mvelpp, log_print)

  ! -----------------------------------------------------------------
  ! Initialize land model
  ! -----------------------------------------------------------------

  call clm_init0()
  call clm_init1()
  call clm_init2()

  ! -----------------------------------------------------------------
  ! Initialize "external" atmospheric forcing
  ! -----------------------------------------------------------------

  ! Read atmospheric forcing dataset one time to obtain the longitudes
  ! and latitudes of the atmospheric dataset, as well as the edges. When
  ! coupled to atm model, these are input variables. If no
  ! atmospheric data files are provided, model uses dummy atmospheric
  ! forcing and sets atmospheric grid to land grid.
  
  if (masterproc) write (6,*) 'Attempting to set up atmospheric grid '
  call atmdrv_init()
  if (masterproc) write (6,*) 'Successfully set up atmospheric grid '
#if (defined FORC_IGSM)
   call opnfil('data4clm',20,'u')
   call opnfil('data4tem',21,'u')
#endif

  call t_stopf('init')
  
  ! -----------------------------------------------------------------
  ! Time stepping loop
  ! -----------------------------------------------------------------

!  call t_barrierf('barrieri',clmmpicom)
  call t_startf('runtotal')

#if (defined COUP_MIT2D)
  first=.false.
 else
#else
  do
#endif

     ! Current atmospheric state and fluxes for all [atmlon] x [atmlat] points.
     nstep = get_nstep()

#if (defined FORC_IGSM)
     call readIGSM2(nstep)
#endif

     call t_startf('atmdrv')
     call atmdrv(nstep)
     call t_stopf('atmdrv')

     !  call t_barrierf('barrier1b',clmmpicom)
     ! Run

     call clm_run1()

     !  call t_barrierf('barrier2b',clmmpicom)

     call clm_run2()

     !  call t_barrierf('barrierd2p',clmmpicom)
     ! Determine if time to stop

#if (!defined COUP_MIT2D)
     if (is_last_step()) exit
#endif

     ! Increment time step

     call advance_timestep()

#if (!defined COUP_MIT2D)
  end do
#else
   call get_curr_date (yr, mon, day, mcsec)

   endif
   if (is_last_step()) then
#endif
  call t_stopf('runtotal')

  ! -----------------------------------------------------------------
  ! Exit gracefully
  ! -----------------------------------------------------------------

#if (defined BGL)
     call print_stack_size()
#endif

  if (masterproc) then
     write(6,*)'SUCCESFULLY TERMINATING CLM MODEL at nstep= ',get_nstep()
  endif
  call t_prf('timing_all',clmmpicom)
  call t_finalizef()

  ! Finalize ESMF
  call ESMF_Finalize()

#if (defined COUP_MIT2D)
  endif
  return
  end subroutine clm4mit2d
#else
  stop
end program program_off
#endif

#else

!The following is only here since empty file won't compile
subroutine program_off_stub
  write(6,*) 'PROGRAM_OFF: this routine should not be called'
  return
end subroutine program_off_stub

#endif
