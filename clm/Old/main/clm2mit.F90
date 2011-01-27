#include <misc.h>
#include <preproc.h>

module clm2mit

#if (defined COUP_MIT2D)

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: clm2mit
!
! !DESCRIPTION: 
!	
!	This module and subsequent subroutines are adapted from 'clm_csmMod.F90'
!	and used to pass CLM outputs to the MIT 2D atmospheric model
!
!				CAS 09/01/04
!	
! Set of routines that define communication between the
! land model and flux coupler. The order of sends/receives is: 
! 1) receive orbital data from coupler
! 2) send control data (grids and masks) to coupler 
!    land grid does not have valid data, runoff grid does
! 3) receive valid land grid from flux coupler
! 4) send compressed runoff information to flux coupler
! 5) start normal send/recv communication patterm   
! 
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use nanMod
  use clm_varpar                      !parameters
#if (defined SPMD)
  use spmdMod, only : masterproc, mpicom, gather_data_to_master, scatter_data_from_master
#else
  use spmdMod, only : masterproc
#endif	
  use shr_msg_mod                     !csm_share message passing routines
  use shr_sys_mod, only: shr_sys_irtc !csm_share system utility routines
  use system_messages, only : allocation_err !allocation error output
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
!  public :: csm_recvorb     !Receive initial control data
!  public :: csm_sendcontrol !Send first control data and first "invalid" grid
!  public :: csm_recvgrid    !Receive grid and land mask 
!  public :: csm_sendrunoff  !Send valid runoff information
!  public :: mit_sendalb     !Send initial albedos, surface temp and snow data 
!  public :: csm_dosndrcv    !logic for determining if send/recv
  public :: mit_recv        !receive data from flux coupler
  public :: mit_send        !send data to flux coupler 
  public :: mit_send1st     !send data to flux coupler 
  public :: mit_flxave      !flux averaging rougine
!  public :: restart_coupler !restart code 
!  public :: compat_check_spval
!  public :: csm_compat      !checks compatibility of messages send/received
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
!
  private
!
! PRIVATE MEMBER FUNCTIONS:
  private :: mit_allocate         !csm dynamic memory allocation
!
! PRIVATE TYPES:
!
! Buffer information 
!
  integer, parameter :: nibuff = 100    ! cpl ->atm msg, initial
  integer, parameter :: ncbuff_max=2000 ! Max size of character data from cpl

  integer   :: ncbuff             !Size of character data from cpl
  integer   :: ibuffr(nibuff)     !Integer buffer from cpl
  integer   :: ibuffs(nibuff)     !Integer buffer to   cpl
  real(r8)  :: rbuff(nibuff)      !Floating pt buffer from cpl
  character :: cbuff(ncbuff_max)  !Character data recieved from cpl
!
! Timing information 
!
  logical :: csm_timing
  integer :: irtc_w               !rc ticks when waiting for msg
  integer :: irtc_r               !rtc tics when msg recved
  integer :: irtc_s               !rtc tics when msg sent
!
! Send/recv buffers 
!
  integer, parameter :: nsend_csm = 14
  integer, parameter :: nrecv_csm = 16
  integer, parameter :: nsend_mit = 60
  integer, parameter :: nrecv_mit = 16

  real(r8) :: send2d(lsmlon,lsmlat,nsend_mit)  !2d send buffer
  real(r8) :: recv2d(lsmlon,lsmlat,nrecv_mit)  !2d recv buffer
  
  real(r8), pointer   :: recv1d_glob(:,:) 
  real(r8), pointer   :: recv1d_loc(:,:)
  real(r8), pointer   :: send1d_glob(:,:)  
  real(r8), pointer   :: send1d_loc(:,:)  

  logical :: debug_flag   !received from coupler
!
! Flux averaging arrays and counters
!
  integer  :: icnt  !step counter for flux averager 
  integer  :: ncnt  !number of steps over which to average output fluxes 
  real(r8) :: rncnt !reciprocal of ncnt

  real(r8), allocatable :: taux_ave(:)     !averaged array
  real(r8), allocatable :: tauy_ave(:)     !averaged array
  real(r8), allocatable :: lhflx_ave(:)    !averaged array
  real(r8), allocatable :: shflx_ave(:)    !averaged array
  real(r8), allocatable :: lwup_ave(:)     !averaged array
  real(r8), allocatable :: qflx_ave(:)     !averaged array
  real(r8), allocatable :: swabs_ave(:)    !averaged array
!
! When to send/receive messages to coupler and when to make restart and stop
!
  integer, private:: ncpday         !number of send/recv calls per day
  logical, public :: dorecv         !receive data from coupler this step
  logical, public :: dosend         !send data to coupler this step
  logical, public :: csmstop_next   !received stop at eod signal and will stop on next ts
  logical, public :: csmstop_now    !received stop now signal from coupler
  logical, public :: csmrstrt       !restart write signal received from coupler
!
! Indices for send/recv fields
!
  integer, parameter :: irecv_hgt    = 1    !zgcmxy       Atm state m    
  integer, parameter :: irecv_u      = 2    !forc_uxy     Atm state m/s  
  integer, parameter :: irecv_v      = 3    !forc_vxy     Atm state m/s  
  integer, parameter :: irecv_th     = 4    !forc_thxy    Atm state K    
  integer, parameter :: irecv_q      = 5    !forc_qxy     Atm state kg/kg
  integer, parameter :: irecv_pbot   = 6    !ptcmxy       Atm state Pa   
  integer, parameter :: irecv_t      = 7    !forc_txy     Atm state K    
  integer, parameter :: irecv_lwrad  = 8    !flwdsxy      Atm flux  W/m^2
  integer, parameter :: irecv_rainc  = 9    !rainxy       Atm flux  mm/s 
  integer, parameter :: irecv_rainl  = 10   !rainxy       Atm flux  mm/s 
  integer, parameter :: irecv_snowc  = 11   !snowfxy      Atm flux  mm/s 
  integer, parameter :: irecv_snowl  = 12   !snowfxl      Atm flux  mm/s 
  integer, parameter :: irecv_soll   = 13   !forc_sollxy  Atm flux  W/m^2
  integer, parameter :: irecv_sols   = 14   !forc_solsxy  Atm flux  W/m^2
  integer, parameter :: irecv_solld  = 15   !forc_solldxy Atm flux  W/m^2
  integer, parameter :: irecv_solsd  = 16   !forc_solsdxy Atm flux  W/m^2

  integer, parameter :: isend_trad   = 1
  integer, parameter :: isend_asdir  = 2
  integer, parameter :: isend_aldir  = 3
  integer, parameter :: isend_asdif  = 4
  integer, parameter :: isend_aldif  = 5 
  integer, parameter :: isend_sno    = 6 
  integer, parameter :: isend_taux   = 7 
  integer, parameter :: isend_tauy   = 8 
  integer, parameter :: isend_lhflx  = 9 
  integer, parameter :: isend_shflx  = 10
  integer, parameter :: isend_lwup   = 11
  integer, parameter :: isend_qflx   = 12
  integer, parameter :: isend_tref2m = 13
  integer, parameter :: isend_swabs  = 14
  integer, parameter :: isend_sro    = 15
  integer, parameter :: isend_ssr    = 16
  integer, parameter :: isend_glr    = 17
  integer, parameter :: isend_tflx   = 18
  integer, parameter :: isend_tgnd   = 19
  integer, parameter :: isend_snc    = 20
  integer, parameter :: isend_vet    = 21
  integer, parameter :: isend_sev    = 22
  integer, parameter :: isend_cev    = 23
  integer, parameter :: isend_emv    = 24
!
  integer, parameter :: isend_h2ol   = 30
  integer, parameter :: isend_h2oi   = 40
  integer, parameter :: isend_tsoi   = 50

!
! MIT 2D COMMON BLOCK
!
  real(r8) :: lhfclm(lsmlon,lsmlat)
  real(r8) :: shfclm(lsmlon,lsmlat)
  real(r8) :: tauxclm(lsmlon,lsmlat)
  real(r8) :: tauyclm(lsmlon,lsmlat)
  real(r8) :: asdirclm(lsmlon,lsmlat)
  real(r8) :: aldirclm(lsmlon,lsmlat)
  real(r8) :: asdifclm(lsmlon,lsmlat)
  real(r8) :: aldifclm(lsmlon,lsmlat)
  real(r8) :: emvclm(lsmlon,lsmlat)
  real(r8) :: sroclm(lsmlon,lsmlat)
  real(r8) :: ssrclm(lsmlon,lsmlat)
  real(r8) :: glrclm(lsmlon,lsmlat)
  real(r8) :: vetclm(lsmlon,lsmlat)
  real(r8) :: sevclm(lsmlon,lsmlat)
  real(r8) :: cevclm(lsmlon,lsmlat)
  real(r8) :: lwuclm(lsmlon,lsmlat)
  real(r8) :: tref2mclm(lsmlon,lsmlat)
  real(r8) :: tflxclm(lsmlon,lsmlat)
  real(r8) :: tgndclm(lsmlon,lsmlat)
  real(r8) :: snwdclm(lsmlon,lsmlat)
  real(r8) :: snwcclm(lsmlon,lsmlat)
  real(r8) :: tsoiclm(lsmlon,lsmlat,nlevsoi)
  real(r8) :: h2olclm(lsmlon,lsmlat,nlevsoi)
  real(r8) :: h2oiclm(lsmlon,lsmlat,nlevsoi)
  common/clm4mit/lhfclm,shfclm,tauxclm,tauyclm,asdirclm,aldirclm,asdifclm,aldifclm,emvclm,&
         sroclm,ssrclm,glrclm,vetclm,sevclm,cevclm,lwuclm,tref2mclm,tflxclm,tgndclm,&
         snwdclm,snwcclm,tsoiclm,h2oiclm,h2olclm

!
! CCSM timers
!
  logical  :: timer_lnd_sendrecv = .false. !true => timer is on
  logical  :: timer_lnd_recvsend = .false. !true => timer is on
!
! Input dtime
!
  integer, public :: csm_dtime  !value passed to coupler on initialization set to namelist input 
                                !consistency check done later that restart file gives same value
!----------------------------------------------------------------------- 

contains

!===============================================================================
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_allocate
!
! !INTERFACE:
  subroutine mit_allocate ()
!
! !DESCRIPTION: 
! Allocate dynamic memory for csm coupling
!
! !USES:
    use clmtype
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
! 
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                 !error code 
    integer :: numg, begg, endg    !grid 1d indices 
!-----------------------------------------------------------------------

    numg = grid1d%num
    begg = grid1d%beg
    endg = grid1d%end

    if (masterproc) then
       allocate (recv1d_glob(nrecv_csm,numg), STAT=ier) 
       call allocation_err(ier, 'csm_allocate', 'recv1d_glob', numg)
       recv1d_glob(:,:) = nan

       allocate (send1d_glob(nsend_mit,numg), STAT=ier)
       call allocation_err(ier, 'csm_allocate', 'send1d_glob', numg)
       send1d_glob(:,:) = nan
    endif

#if (defined SPMD)
    allocate (send1d_loc(nsend_mit,begg:endg), STAT=ier)
    call allocation_err(ier, 'csm_allocate', 'send1d_loc', endg-begg+1)
    send1d_loc(:,:) = nan

    allocate (recv1d_loc(nrecv_csm,begg:endg), STAT=ier) 
    call allocation_err(ier, 'csm_allocate', 'recv1d_loc', endg-begg+1)
    recv1d_loc(:,:) = nan
#else
    send1d_loc => send1d_glob
    recv1d_loc => recv1d_glob
#endif

    return
  end subroutine mit_allocate
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_recv
! 
! !INTERFACE:
  subroutine mit_recv()
!
! !DESCRIPTION: 
!  Recveive and map data from flux coupler
!
! !USES:
    use clmpoint
    use clm_varcon, only : rair, po2, pco2
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: i,j,n,gi,indexg               !indices 
    real(r8) :: forc_rainc                    !rainxy Atm flux mm/s   
    real(r8) :: forc_rainl                    !rainxy Atm flux mm/s   
    real(r8) :: forc_snowc                    !snowfxy Atm flux  mm/s 
    real(r8) :: forc_snowl                    !snowfxl Atm flux  mm/s 
    integer  :: ier                           !return error code
    type(gridcell_type)     , pointer :: g    !pointer to derived subtype
    type(atm2lnd_state_type), pointer :: a2ls !pointer to derived subtype
    type(atm2lnd_flux_type) , pointer :: a2lf !pointer to derived subtype
!-----------------------------------------------------------------------

    ! Start timers

     if (timer_lnd_sendrecv) then
        call t_stopf ('lnd_sendrecv') ; timer_lnd_sendrecv = .false.
     endif

     call t_startf('lnd_recv')

     ! Receive message from flux coupler

     if (masterproc) then
        ibuffr(:)     = 0            
        recv2d(:,:,:) = 1.e30
        if (csm_timing) irtc_w = shr_sys_irtc()
!        call shr_msg_recv_i (ibuffr, size(ibuffr), SHR_MSG_TID_CPL, SHR_MSG_TAG_C2L)
!        call shr_msg_recv_r (recv2d, size(recv2d), SHR_MSG_TID_CPL, SHR_MSG_TAG_C2L)
        if (csm_timing) then
           irtc_r = shr_sys_irtc()
           write(6,9099) irtc_w,'d->l waiting'
9099       format('[mp timing]  irtc = ',i20,' ',a)
           write(6,9099) irtc_r,'d->l received'
        end if
        
        ! Do global integrals of fluxes if flagged

        if (debug_flag) then	
           write(6,*)
           write(6,100) 'lnd','recv', irecv_lwrad, global_sum(recv2d(1,1,irecv_lwrad),1.e30), ' lwrad'
           write(6,100) 'lnd','recv', irecv_rainc, global_sum(recv2d(1,1,irecv_rainc),1.e30), ' rainc'
           write(6,100) 'lnd','recv', irecv_rainl, global_sum(recv2d(1,1,irecv_rainl),1.e30), ' rainl'
           write(6,100) 'lnd','recv', irecv_snowc, global_sum(recv2d(1,1,irecv_snowc),1.e30), ' snowc'
           write(6,100) 'lnd','recv', irecv_snowl, global_sum(recv2d(1,1,irecv_snowl),1.e30), ' snowl'
           write(6,100) 'lnd','recv', irecv_soll , global_sum(recv2d(1,1,irecv_soll ),1.e30), ' soll '
           write(6,100) 'lnd','recv', irecv_sols , global_sum(recv2d(1,1,irecv_sols ),1.e30), ' sols '
           write(6,100) 'lnd','recv', irecv_solld, global_sum(recv2d(1,1,irecv_solld),1.e30), ' solld'
           write(6,100) 'lnd','recv', irecv_solsd, global_sum(recv2d(1,1,irecv_solsd),1.e30), ' solsd'
100        format('comm_diag ',a3,1x,a4,1x,i3,es26.19,a)
           write(6,*)
        endif
     endif  ! end of if-masteproc

#if (defined SPMD)
     call mpi_bcast (ibuffr, size(ibuffr), MPI_INTEGER, 0, mpicom, ier)    
     if (ier /= MPI_SUCCESS) then
        write(6,*) 'MPI BCAST ERROR: for ibuffr in csm_recv'
        call endrun
     end if
#endif

     ! Stop timer

     call t_stopf('lnd_recv') 

     ! Check if end of run now, if so stop (each processor does this) 

     csmstop_now = .false.
     if (ibuffr(3) /= 0) then 
        csmstop_now = .true.
        if (timer_lnd_recvsend) call t_stopf('lnd_recvsend')
        if (timer_lnd_sendrecv) call t_stopf('lnd_sendrecv')
        write(6,*)'(CSM_RECV) stop now signal from flux coupler'
        write(6,*)'(CSM_RECV) ibuffr(3) = ',ibuffr(3)
        if (masterproc) then
           write(6,9001)
           write(6,9002) ibuffr(4)
           write(6,9003)
9001       format(/////' ===========> Terminating CLM Model')
9002       format(     '      Date: ',i8)
9003       format(/////' <=========== CLM Model Terminated')
        endif
        RETURN
     endif

     ! More timer logic
     
     if (.not. timer_lnd_recvsend) then
        call t_startf('lnd_recvsend') ; timer_lnd_recvsend = .true.
     endif
     
     ! Map 2d received fields on [lsmlon]x[lsmlat] grid to subgrid vectors
     
     if (masterproc) then
        do gi = 1,grid1d%num
           i = grid1d%ixy(gi)
           j = grid1d%jxy(gi)
           recv1d_glob(:nrecv_csm,gi) = recv2d(i,j,:nrecv_csm)
       end do
    endif
    
#if (defined SPMD)
    call scatter_data_from_master (recv1d_loc, recv1d_glob, clmlevel=grid1d%name)
#endif

    ! Split data from coupler into component arrays. 
    ! Note that the precipitation fluxes received  from the coupler 
    ! are in units of kg/s/m^2. To convert these precipitation rates 
    ! in units of mm/sec, one must divide by 1000 kg/m^3 and multiply 
    ! by 1000 mm/m resulting in an overall factor of unity. 
    ! Below the units are therefore given in mm/s.


!$OMP PARALLEL DO PRIVATE (gi,g,a2ls,a2lf,indexg,forc_rainc,forc_rainl,forc_snowc,forc_snowl)
    do gi = 1,clm%mps%ngridcells
       g => clm%g(gi)					
       a2ls => g%a2ls
       a2lf => g%a2lf
       indexg = g%gps%index1d

       a2ls%forc_hgt      = recv1d_loc(irecv_hgt  ,indexg)    !zgcmxy  Atm state m
       a2ls%forc_u        = recv1d_loc(irecv_u    ,indexg)    !forc_uxy  Atm state m/s
       a2ls%forc_v        = recv1d_loc(irecv_v    ,indexg)    !forc_vxy  Atm state m/s
       a2ls%forc_th       = recv1d_loc(irecv_th   ,indexg)    !forc_thxy Atm state K
       a2ls%forc_q        = recv1d_loc(irecv_q    ,indexg)    !forc_qxy  Atm state kg/kg
       a2ls%forc_pbot     = recv1d_loc(irecv_pbot ,indexg)    !ptcmxy  Atm state Pa
       a2ls%forc_t        = recv1d_loc(irecv_t    ,indexg)    !forc_txy  Atm state K
       a2lf%forc_lwrad    = recv1d_loc(irecv_lwrad,indexg)    !flwdsxy Atm flux  W/m^2
       forc_rainc         = recv1d_loc(irecv_rainc,indexg)    !mm/s 
       forc_rainl         = recv1d_loc(irecv_rainl,indexg)    !mm/s 
       forc_snowc         = recv1d_loc(irecv_snowc,indexg)    !mm/s
       forc_snowl         = recv1d_loc(irecv_snowl,indexg)    !mm/s
       a2lf%forc_solad(2) = recv1d_loc(irecv_soll ,indexg)    !forc_sollxy  Atm flux  W/m^2
       a2lf%forc_solad(1) = recv1d_loc(irecv_sols ,indexg)    !forc_solsxy  Atm flux  W/m^2 
       a2lf%forc_solai(2) = recv1d_loc(irecv_solld,indexg)    !forc_solldxy Atm flux  W/m^2
       a2lf%forc_solai(1) = recv1d_loc(irecv_solsd,indexg)    !forc_solsdxy Atm flux  W/m^2

       ! Determine derived quantities

       a2ls%forc_hgt_u = a2ls%forc_hgt    !observational height of wind [m] 
       a2ls%forc_hgt_t = a2ls%forc_hgt    !observational height of temperature [m]  
       a2ls%forc_hgt_q = a2ls%forc_hgt    !observational height of humidity [m]      
       a2ls%forc_vp    = a2ls%forc_q * a2ls%forc_pbot / (0.622 + 0.378 * a2ls%forc_q)
       a2ls%forc_rho   = (a2ls%forc_pbot - 0.378 * a2ls%forc_vp) / (rair * a2ls%forc_t)
       a2ls%forc_co2   = pco2 * a2ls%forc_pbot
       a2ls%forc_o2    = po2 * a2ls%forc_pbot
       a2ls%forc_wind  = sqrt(a2ls%forc_u**2 + a2ls%forc_v**2)
       a2lf%forc_solar = a2lf%forc_solad(1) + a2lf%forc_solai(1) + &
                         a2lf%forc_solad(2) + a2lf%forc_solai(2)

       ! Determine precipitation needed by clm

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
    end do

    return
  end subroutine mit_recv

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mit_send
!
! !INTERFACE:
  subroutine mit_send()
!   	
! !DESCRIPTION: 
! Send data to MIT common block
!
! !USES:
    use clmtype
    use clmpoint
    use clm_varctl, only : csm_doflxave
    use clm_varsur, only : landmask
    use time_manager, only : get_curr_date
    use clm_varpar, only: nlevsoi
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
!
! !LOCAL VARIABLES:
!
! MIT 2D COMMON BLOCK
!
  real(r8) :: lhfclm(lsmlon,lsmlat)
  real(r8) :: shfclm(lsmlon,lsmlat)
  real(r8) :: tauxclm(lsmlon,lsmlat)
  real(r8) :: tauyclm(lsmlon,lsmlat)
  real(r8) :: asdirclm(lsmlon,lsmlat)
  real(r8) :: aldirclm(lsmlon,lsmlat)
  real(r8) :: asdifclm(lsmlon,lsmlat)
  real(r8) :: aldifclm(lsmlon,lsmlat)
!  real(r8) :: emvclm(lsmlon,lsmlat)
  real(r8) :: sroclm(lsmlon,lsmlat)
  real(r8) :: ssrclm(lsmlon,lsmlat)
  real(r8) :: glrclm(lsmlon,lsmlat)
  real(r8) :: vetclm(lsmlon,lsmlat)
  real(r8) :: sevclm(lsmlon,lsmlat)
  real(r8) :: cevclm(lsmlon,lsmlat)
  real(r8) :: lwuclm(lsmlon,lsmlat)
  real(r8) :: tref2mclm(lsmlon,lsmlat)
  real(r8) :: tflxclm(lsmlon,lsmlat)
  real(r8) :: tgndclm(lsmlon,lsmlat)
  real(r8) :: snwdclm(lsmlon,lsmlat)
  real(r8) :: snwcclm(lsmlon,lsmlat)
  real(r8) :: tsoiclm(lsmlon,lsmlat,nlevsoi)
  real(r8) :: h2olclm(lsmlon,lsmlat,nlevsoi)
  real(r8) :: h2oiclm(lsmlon,lsmlat,nlevsoi)
  common/clm4mit/lhfclm,shfclm,tauxclm,tauyclm,asdirclm,aldirclm,asdifclm,aldifclm,&
         sroclm,ssrclm,glrclm,vetclm,sevclm,cevclm,lwuclm,tref2mclm,tflxclm,tgndclm,&
         snwdclm,snwcclm,tsoiclm,h2oiclm,h2olclm
  real*4 :: c2mout(lsmlon,lsmlat)
    logical :: first
    integer :: i,j,n,gi,ci,pi   !indices
    integer :: irec             !indices
    integer :: lev              !soil layer indice
    integer :: yr               !current year 
    integer :: mon              !current month 
    integer :: day              !current day (0, 1, ...)
    integer :: ncsec            !current seconds of current date (0, ..., 86400)
    integer :: ncdate           !current date (yymmdd format) (e.g., 021105)
    type(gridcell_type)       , pointer :: g   !pointer to derived subtype
    type(landunit_type)       , pointer :: l   !pointer to derived subtype
    type(column_type)         , pointer :: c   !pointer to derived subtype
    type(pft_type)            , pointer :: p   !pointer to derived subtype
    type(gridcell_pstate_type), pointer :: gps !pointer to derived subtype
    type(landunit_pstate_type), pointer :: lps !pointer to derived subtype
    type(column_pstate_type)  , pointer :: cps !pointer to derived subtype
    type(column_wstate_type)  , pointer :: cws !pointer to derived subtype
    type(column_estate_type)  , pointer :: ces !pointer to derived subtype
    type(column_eflux_type)   , pointer :: cef !pointer to derived subtype
    type(column_wflux_type)   , pointer :: cwf !pointer to derived subtype
    type(pft_pstate_type)     , pointer :: pps !pointer to derived subtype
    type(pft_mflux_type)      , pointer :: pmf !pointer to derived subtype
    type(pft_wflux_type)      , pointer :: pwf !pointer to derived subtype
    type(pft_eflux_type)      , pointer :: pef !pointer to derived subtype 
    type(pft_estate_type)     , pointer :: pes !pointer to derived subtype
    data irec/1/
    data first/.true./
!-----------------------------------------------------------------------

    ! Determine 1d vector of fields that will be sent to coupler.
    ! Coupler has convention that fluxes are positive downward.
    if (first) then
      call mit_allocate()
      first = .false.
    endif
    send1d_loc(:,:) = 0.

    do ci = cols1d%beg,cols1d%end
       c => cpoint%col(ci)%c
       cps => c%cps
       cws => c%cws
       ces => c%ces
       cwf => c%cwf
       cef => c%cef
       gi = cols1d%gindex(ci)
       send1d_loc(isend_sno,gi) = send1d_loc(isend_sno,gi) + cws%h2osno/1000. * cps%wtxy 
       send1d_loc(isend_snc,gi) = send1d_loc(isend_snc,gi) + cps%frac_sno * cps%wtxy 
       send1d_loc(isend_sro,gi) = send1d_loc(isend_sro,gi) + cwf%qflx_surf * cps%wtxy
       send1d_loc(isend_ssr,gi) = send1d_loc(isend_ssr,gi) + cwf%qflx_drain * cps%wtxy
       send1d_loc(isend_glr,gi) = send1d_loc(isend_glr,gi) + cwf%qflx_qrgwl * cps%wtxy
       send1d_loc(isend_cev,gi) = send1d_loc(isend_cev,gi) + cwf%pwf_a%qflx_evap_can * cps%wtxy
       send1d_loc(isend_sev,gi) = send1d_loc(isend_sev,gi) + cwf%pwf_a%qflx_evap_soi * cps%wtxy
       !send1d_loc(isend_vet,gi) = send1d_loc(isend_vet,gi) + cwf%pwf_a%qflx_evap_tot * cps%wtxy
       send1d_loc(isend_tgnd,gi) = send1d_loc(isend_tgnd,gi) + ces%t_grnd * cps%wtxy
       send1d_loc(isend_lhflx,gi) = send1d_loc(isend_lhflx,gi) - cef%pef_a%eflx_lh_tot * cps%wtxy 
       do lev = 1,nlevsoi
         !CAS
	 !CAS  A check must be made to see if column is a water-type of lake, if so, adjust grid-weighted value accordingly
         !CAS
         if ( cols1d%itypwat(ci) /= 3 ) then
           send1d_loc(isend_h2ol+lev,gi) = send1d_loc(isend_h2ol+lev,gi) + cws%h2osoi_liq(lev) * cps%wtxy
           send1d_loc(isend_h2oi+lev,gi) = send1d_loc(isend_h2oi+lev,gi) + cws%h2osoi_ice(lev) * cps%wtxy
           send1d_loc(isend_tsoi+lev,gi) = send1d_loc(isend_tsoi+lev,gi) + ces%t_soisno(lev) * cps%wtxy
         else
           !print *,'WTXY @ grid = ',gi,' and ITYPE = ',cols1d%itypwat(ci),' : ',cps%wtxy
           send1d_loc(isend_h2ol+lev,gi) = send1d_loc(isend_h2ol+lev,gi)/(1-cps%wtxy)
           send1d_loc(isend_h2oi+lev,gi) = send1d_loc(isend_h2oi+lev,gi)/(1-cps%wtxy)
           send1d_loc(isend_tsoi+lev,gi) = send1d_loc(isend_tsoi+lev,gi)/(1-cps%wtxy)
         endif
         !print *,'H2OSOI @ grid = ',gi,' and lev = ',lev,' : ',cws%h2osoi_liq(lev)
         !print *,'TSOI @ grid = ',gi,' and lev = ',lev,' : ',ces%t_soisno(lev)
       enddo
    end do

    do pi = pfts1d%beg,pfts1d%end
       p => ppoint%pft(pi)%p
       pps => p%pps
       pes => p%pes
       pef => p%pef
       pmf => p%pmf
       pwf => p%pwf
       gi = pfts1d%gindex(pi)
       send1d_loc(isend_trad  ,gi) = send1d_loc(isend_trad ,gi) + pes%t_rad_pft * pps%wtxy 
       send1d_loc(isend_asdir ,gi) = send1d_loc(isend_asdir,gi) + pps%albd(1) * pps%wtxy  
       send1d_loc(isend_aldir ,gi) = send1d_loc(isend_aldir,gi) + pps%albd(2) * pps%wtxy 
       send1d_loc(isend_asdif ,gi) = send1d_loc(isend_asdif,gi) + pps%albi(1) * pps%wtxy 
       send1d_loc(isend_aldif ,gi) = send1d_loc(isend_aldif,gi) + pps%albi(2) * pps%wtxy 
       send1d_loc(isend_emv ,gi) = send1d_loc(isend_emv,gi) + pps%emv * pps%wtxy 
       send1d_loc(isend_tref2m,gi) = send1d_loc(isend_tref2m,gi) + pes%t_ref2m * pps%wtxy 
       send1d_loc(isend_tflx,gi) = send1d_loc(isend_tflx,gi) + pes%t_af * pps%wtxy 
       send1d_loc(isend_taux ,gi) = send1d_loc(isend_taux ,gi) - pmf%taux * pps%wtxy 
       send1d_loc(isend_tauy ,gi) = send1d_loc(isend_tauy ,gi) - pmf%tauy * pps%wtxy 
       send1d_loc(isend_vet,gi) = send1d_loc(isend_vet,gi) + pwf%qflx_tran_veg * pps%wtxy
       send1d_loc(isend_shflx,gi) = send1d_loc(isend_shflx,gi) - pef%eflx_sh_tot * pps%wtxy 
       send1d_loc(isend_lwup ,gi) = send1d_loc(isend_lwup ,gi) - pef%eflx_lwrad_out * pps%wtxy 
    end do

    ! Map fields from 1d-grid vector to 2d-grid
    ! NOTE: snow is sent as zero over non-land because currently 
    ! the ocn and sea-ice send no snow cover to cpl and so the cpl 
    ! sends back zero snow over non-land to  the atm (the atm and 
    ! land grid are currently assumed to be identical)
    
!       do n = 1,nsend_mit
!          where( landmask(:,:) > 0 ) 
!             send2d(:,:,n) = 0.
!          elsewhere
!             send2d(:,:,n) = 1.e30
!          endwhere
!       end do
!       send2d(:,:,isend_sno) = 0.     !reset snow to 0 everywhere
       do gi = 1, grid1d%num
         i = grid1d%ixy(gi)
         j = grid1d%jxy(gi)
         lhfclm(i,j) = send1d_loc(isend_lhflx,gi)
         shfclm(i,j) = send1d_loc(isend_shflx,gi)
         tauxclm(i,j) = send1d_loc(isend_taux,gi)
         tauyclm(i,j) = send1d_loc(isend_tauy,gi)
         asdirclm(i,j) = send1d_loc(isend_asdir,gi)
         aldirclm(i,j) = send1d_loc(isend_aldir,gi)
         asdifclm(i,j) = send1d_loc(isend_asdif,gi)
         aldifclm(i,j) = send1d_loc(isend_aldif,gi)
         emvclm(i,j) = send1d_loc(isend_emv,gi)
         sroclm(i,j) = send1d_loc(isend_sro,gi)
         ssrclm(i,j) = send1d_loc(isend_ssr,gi)
         glrclm(i,j) = send1d_loc(isend_glr,gi)
         vetclm(i,j) = send1d_loc(isend_vet,gi)
         cevclm(i,j) = send1d_loc(isend_cev,gi)
         sevclm(i,j) = send1d_loc(isend_sev,gi)
         lwuclm(i,j) = send1d_loc(isend_lwup,gi)
         tref2mclm(i,j) = send1d_loc(isend_tref2m,gi)
         tflxclm(i,j) = send1d_loc(isend_tflx,gi)
         tgndclm(i,j) = send1d_loc(isend_tgnd,gi)
         snwdclm(i,j) = send1d_loc(isend_sno,gi)
         snwcclm(i,j) = send1d_loc(isend_snc,gi)
         do lev = 1,nlevsoi
           h2olclm(i,j,lev) = send1d_loc(isend_h2ol+lev,gi) 
           h2oiclm(i,j,lev) = send1d_loc(isend_h2oi+lev,gi)
           tsoiclm(i,j,lev) = send1d_loc(isend_tsoi+lev,gi)
         enddo
       end do

!
!	write data for consistency check
!
!       c2mout=lhfclm
!       write (111,rec=irec) c2mout
!       irec=irec+1
!       c2mout=shfclm
!       write (111,rec=irec) c2mout
!       irec=irec+1
!       c2mout=tauxclm
!       write (111,rec=irec) c2mout
!       irec=irec+1
!       c2mout=tauyclm
!       write (111,rec=irec) c2mout
!       irec=irec+1
!       c2mout=asdirclm
!       write (111,rec=irec) c2mout
!       irec=irec+1
!       c2mout=aldirclm
!       write (111,rec=irec) c2mout
!       irec=irec+1
!       c2mout=asdifclm
!       write (111,rec=irec) c2mout
!       irec=irec+1
!       c2mout=aldifclm
!       write (111,rec=irec) c2mout
!       irec=irec+1
!       c2mout=emvclm
!       write (111,rec=irec) c2mout
!       irec=irec+1
!       c2mout=vetclm
!       write (111,rec=irec) c2mout
!       irec=irec+1
!       c2mout=cevclm
!       write (111,rec=irec) c2mout
!       irec=irec+1
!       c2mout=sevclm
!       write (111,rec=irec) c2mout
!       irec=irec+1
!       c2mout=sroclm
!       write (111,rec=irec) c2mout
!       irec=irec+1
!       c2mout=ssrclm
!       write (111,rec=irec) c2mout
!       irec=irec+1
!       c2mout=glrclm
!       write (111,rec=irec) c2mout
!       irec=irec+1
!       c2mout=lwuclm
!       write (111,rec=irec) c2mout
!       irec=irec+1
!       c2mout=tref2mclm
!       write (111,rec=irec) c2mout
!       irec=irec+1
!       c2mout=tflxclm
!       write (111,rec=irec) c2mout
!       irec=irec+1
!       c2mout=tgndclm
!       write (111,rec=irec) c2mout
!       irec=irec+1
!       c2mout=snwdclm
!       write (111,rec=irec) c2mout
!       irec=irec+1
!       c2mout=snwcclm
!       write (111,rec=irec) c2mout
!       irec=irec+1
!       do lev = 1,nlevsoi
!         c2mout=h2olclm(:,:,lev)
!         print *,'H2OSOI TO MIT2D: lev = ',lev
!         print *,c2mout
!         write (111,rec=irec) c2mout
!         irec=irec+1
!       enddo
!       do lev = 1,nlevsoi
!         c2mout=h2oiclm(:,:,lev)
!         print *,'H2OSOI TO MIT2D: lev = ',lev
!         print *,c2mout
!         write (111,rec=irec) c2mout
!         irec=irec+1
!       enddo
!       do lev = 1,nlevsoi
!         c2mout=tsoiclm(:,:,lev)
!         print *,'TSOI TO MIT2D: lev = ',lev
!         print *,c2mout
!         write (111,rec=irec) c2mout
!         irec=irec+1
!       enddo
  
     ! Do global integrals if flag is set

       if (debug_flag) then	
          write(6,*)
          write(6,100) 'lnd','send', isend_taux , global_sum(send2d(1,1,isend_taux ),1.e30), ' taux'
          write(6,100) 'lnd','send', isend_tauy , global_sum(send2d(1,1,isend_tauy ),1.e30), ' tauy'
          write(6,100) 'lnd','send', isend_lhflx, global_sum(send2d(1,1,isend_lhflx),1.e30), ' lhflx'
          write(6,100) 'lnd','send', isend_shflx, global_sum(send2d(1,1,isend_shflx),1.e30), ' shflx'
          write(6,100) 'lnd','send', isend_lwup , global_sum(send2d(1,1,isend_lwup ),1.e30), ' lwup'
          write(6,100) 'lnd','send', isend_qflx , global_sum(send2d(1,1,isend_qflx ),1.e30), ' qflx'
          write(6,100) 'lnd','send', isend_swabs, global_sum(send2d(1,1,isend_swabs),1.e30), ' swabs'
          write(6,*)
100       format('comm_diag ',a3,1x,a4,1x,i3,es26.19,a)
       endif
   return
   end subroutine mit_send

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mit_send1st
!
! !INTERFACE:
  subroutine mit_send1st()
!   	
! !DESCRIPTION: 
! Send data to MIT common block
!
! !USES:
    use clmtype
    use clmpoint
    use clm_varctl, only : csm_doflxave
    use clm_varsur, only : landmask
    use time_manager, only : get_curr_date
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
!
! !LOCAL VARIABLES:
!
! MIT 2D COMMON BLOCK
!
  real(r8) :: lhfclm(lsmlon,lsmlat)
  real(r8) :: shfclm(lsmlon,lsmlat)
  real(r8) :: tauxclm(lsmlon,lsmlat)
  real(r8) :: tauyclm(lsmlon,lsmlat)
  real(r8) :: asdirclm(lsmlon,lsmlat)
  real(r8) :: aldirclm(lsmlon,lsmlat)
  real(r8) :: asdifclm(lsmlon,lsmlat)
  real(r8) :: aldifclm(lsmlon,lsmlat)
  real(r8) :: emvclm(lsmlon,lsmlat)
  real(r8) :: sroclm(lsmlon,lsmlat)
  real(r8) :: ssrclm(lsmlon,lsmlat)
  real(r8) :: glrclm(lsmlon,lsmlat)
  real(r8) :: vetclm(lsmlon,lsmlat)
  real(r8) :: sevclm(lsmlon,lsmlat)
  real(r8) :: cevclm(lsmlon,lsmlat)
  real(r8) :: lwuclm(lsmlon,lsmlat)
  real(r8) :: tref2mclm(lsmlon,lsmlat)
  real(r8) :: tflxclm(lsmlon,lsmlat)
  real(r8) :: tgndclm(lsmlon,lsmlat)
  real(r8) :: snwdclm(lsmlon,lsmlat)
  real(r8) :: snwcclm(lsmlon,lsmlat)
  real(r8) :: tsoiclm(lsmlon,lsmlat,nlevsoi)
  real(r8) :: h2olclm(lsmlon,lsmlat,nlevsoi)
  real(r8) :: h2oiclm(lsmlon,lsmlat,nlevsoi)
  common/clm4mit/lhfclm,shfclm,tauxclm,tauyclm,asdirclm,aldirclm,asdifclm,aldifclm,emvclm,&
         sroclm,ssrclm,glrclm,vetclm,sevclm,cevclm,lwuclm,tref2mclm,tflxclm,tgndclm,&
         snwdclm,snwcclm,tsoiclm,h2oiclm,h2olclm
  real*4 :: c2mout(lsmlon,lsmlat)
    logical :: first
    integer :: i,j,n,gi,ci,pi   !indices
    integer :: irec             !indices
    integer :: lev              !soil layer indice
    integer :: yr               !current year 
    integer :: mon              !current month 
    integer :: day              !current day (0, 1, ...)
    integer :: ncsec            !current seconds of current date (0, ..., 86400)
    integer :: ncdate           !current date (yymmdd format) (e.g., 021105)
    type(gridcell_type)       , pointer :: g   !pointer to derived subtype
    type(landunit_type)       , pointer :: l   !pointer to derived subtype
    type(column_type)         , pointer :: c   !pointer to derived subtype
    type(pft_type)            , pointer :: p   !pointer to derived subtype
    type(gridcell_pstate_type), pointer :: gps !pointer to derived subtype
    type(landunit_pstate_type), pointer :: lps !pointer to derived subtype
    type(column_pstate_type)  , pointer :: cps !pointer to derived subtype
    type(column_wstate_type)  , pointer :: cws !pointer to derived subtype
    type(column_estate_type)  , pointer :: ces !pointer to derived subtype
    type(pft_pstate_type)     , pointer :: pps !pointer to derived subtype
    data irec/1/

    ! Send data to the MIT common block
    ! Allocate csm dynamic memory

    call mit_allocate()
    send1d_loc(:,:) = 0.
    do pi = pfts1d%beg,pfts1d%end
       p => ppoint%pft(pi)%p
       pps => p%pps
       gi = pfts1d%gindex(pi)
       send1d_loc(isend_asdir ,gi) = send1d_loc(isend_asdir,gi) + pps%albd(1) * pps%wtxy  
       send1d_loc(isend_aldir ,gi) = send1d_loc(isend_aldir,gi) + pps%albd(2) * pps%wtxy 
       send1d_loc(isend_asdif ,gi) = send1d_loc(isend_asdif,gi) + pps%albi(1) * pps%wtxy 
       send1d_loc(isend_aldif ,gi) = send1d_loc(isend_aldif,gi) + pps%albi(2) * pps%wtxy 
    end do
    do gi = 1, grid1d%num
      i = grid1d%ixy(gi)
      j = grid1d%jxy(gi)
      asdirclm(i,j) = send1d_loc(isend_asdir,gi)
      aldirclm(i,j) = send1d_loc(isend_aldir,gi)
      asdifclm(i,j) = send1d_loc(isend_asdif,gi)
      aldifclm(i,j) = send1d_loc(isend_aldif,gi)
    end do
    do j = 1,lsmlat
      write(6,*) j,asdirclm(1,j)
      write(6,*) j,aldirclm(1,j)
      write(6,*) j,asdifclm(1,j)
      write(6,*) j,aldifclm(1,j)
    enddo
  return
  end subroutine mit_send1st

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_flxave
!
! !INTERFACE:
  subroutine mit_flxave() 
!   	
! !DESCRIPTION: 
! Average output fluxes for flux coupler
! Add land surface model output fluxes to accumulators every time step.
! When icnt==ncnt, compute the average flux over the time interval.
!
! !USES:
!
    use clmtype
    use clmpoint
    use clm_varctl, only : irad
    use time_manager, only : get_nstep
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: pi             !indices
    integer :: nstep          !model time step
    integer :: begp,endp      !pft 1d indices
    type(pft_type)       , pointer :: p   !pointer to derived subtype
    type(pft_pstate_type), pointer :: pps !pointer to derived subtype
    type(pft_mflux_type) , pointer :: pmf !pointer to derived subtype
    type(pft_wflux_type) , pointer :: pwf !pointer to derived subtype
    type(pft_eflux_type) , pointer :: pef !pointer to derived subtype 
    type(pft_estate_type), pointer :: pes !pointer to derived subtype
!-----------------------------------------------------------------------

    ! Allocate dynamic memory if necessary

    begp = pfts1d%beg
    endp = pfts1d%end

    if (.not. allocated(taux_ave)) then
       allocate (taux_ave(begp:endp)) ; taux_ave(:) = nan 
    endif
    if (.not. allocated(tauy_ave)) then
       allocate (tauy_ave(begp:endp)) ; tauy_ave(:) = nan 
    endif
    if (.not. allocated(lhflx_ave)) then
       allocate (lhflx_ave(begp:endp)); lhflx_ave(:) = nan 
    endif
    if (.not. allocated(shflx_ave)) then
       allocate (shflx_ave(begp:endp)); shflx_ave(:) = nan 
    endif
    if (.not. allocated(lwup_ave)) then
       allocate (lwup_ave(begp:endp)) ; lwup_ave(:) = nan 
    endif
    if (.not. allocated(qflx_ave)) then
       allocate (qflx_ave(begp:endp)) ; qflx_ave(:) = nan 
    endif
    if (.not. allocated(swabs_ave)) then
       allocate (swabs_ave(begp:endp)) ; swabs_ave(:) = nan 
    endif

    ! Determine output flux averaging interval

    nstep = get_nstep()
    if (dorecv) then
       icnt = 1
       if ( nstep==0 ) then
          ncnt = irad + 1
       else
          ncnt = irad
       endif
       rncnt = 1./ncnt
    endif

    if (icnt == 1) then
       
       ! Initial call of averaging interval, copy data to accumulators

       do pi = begp,endp
          p => ppoint%pft(pi)%p
          pps => p%pps
          pwf => p%pwf
          pef => p%pef
          pes => p%pes
          pmf => p%pmf
          taux_ave(pi)  = pmf%taux 
          tauy_ave(pi)  = pmf%tauy 
          lhflx_ave(pi) = pef%eflx_lh_tot
          shflx_ave(pi) = pef%eflx_sh_tot
          lwup_ave(pi)  = pef%eflx_lwrad_out
          qflx_ave(pi)  = pwf%qflx_evap_tot
          swabs_ave(pi) = pef%fsa
       end do

    else if (icnt == ncnt) then

       ! Final call of averaging interval, complete averaging 

       do pi = begp,endp
          p => ppoint%pft(pi)%p
          pps => p%pps
          pwf => p%pwf
          pef => p%pef
          pes => p%pes
          pmf => p%pmf
          taux_ave(pi)  = rncnt * (taux_ave(pi)  + pmf%taux) 
          tauy_ave(pi)  = rncnt * (tauy_ave(pi)  + pmf%tauy) 
          lhflx_ave(pi) = rncnt * (lhflx_ave(pi) + pef%eflx_lh_tot)
          shflx_ave(pi) = rncnt * (shflx_ave(pi) + pef%eflx_sh_tot)
          lwup_ave(pi)  = rncnt * (lwup_ave(pi)  + pef%eflx_lwrad_out)
          qflx_ave(pi)  = rncnt * (qflx_ave(pi)  + pwf%qflx_evap_tot)
          swabs_ave(pi) = rncnt * (swabs_ave(pi) + pef%fsa)
       end do

    else

       ! Intermediate call, add data to accumulators

       do pi = begp,endp
          p => ppoint%pft(pi)%p
          pps => p%pps
          pwf => p%pwf
          pef => p%pef
          pes => p%pes
          pmf => p%pmf
          taux_ave(pi)  = (taux_ave(pi)  + pmf%taux) 
          tauy_ave(pi)  = (tauy_ave(pi)  + pmf%tauy) 
          lhflx_ave(pi) = (lhflx_ave(pi) + pef%eflx_lh_tot)
          shflx_ave(pi) = (shflx_ave(pi) + pef%eflx_sh_tot)
          lwup_ave(pi)  = (lwup_ave(pi)  + pef%eflx_lwrad_out)
          qflx_ave(pi)  = (qflx_ave(pi)  + pwf%qflx_evap_tot)
          swabs_ave(pi) = (swabs_ave(pi) + pef%fsa)
       end do

    end if
    
    ! Increment counter

    icnt = icnt + 1
    
    return
  end subroutine mit_flxave

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: global_sum
!
! !INTERFACE:
  real(r8) function global_sum(flux, spval)
!
! !DESCRIPTION: 
! Performs a global sum of fluxes
!
! !USES:
    use clm_varsur, only : area                 !km^2
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: flux(lsmlon,lsmlat) !W/m2, Kg/m2-s or N/m2
    real(r8), intent(in) :: spval               !points to not include in global sum
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j         !indices
!------------------------------------------------------------------------

    global_sum = 0.
    do j = 1,lsmlat
       do i = 1,lsmlon
          if (flux(i,j) /= spval) then
             global_sum = global_sum + flux(i,j)*area(i,j)*1.e6
          endif
       end do
    end do

    return
  end function global_sum

!===============================================================================

#endif

end module clm2mit
