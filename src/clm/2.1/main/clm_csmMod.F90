#include <misc.h>
#include <preproc.h>

module clm_csmMod

#if (defined COUP_CSM)

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: clm_csmMod
!
! !DESCRIPTION: 
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
  public :: csm_recvorb     !Receive initial control data
  public :: csm_sendcontrol !Send first control data and first "invalid" grid
  public :: csm_recvgrid    !Receive grid and land mask 
  public :: csm_sendrunoff  !Send valid runoff information
  public :: csm_sendalb     !Send initial albedos, surface temp and snow data 
  public :: csm_dosndrcv    !logic for determining if send/recv
  public :: csm_recv        !receive data from flux coupler
  public :: csm_send        !send data to flux coupler 
  public :: csm_flxave      !flux averaging rougine
  public :: restart_coupler !restart code 
  public :: compat_check_spval
  public :: csm_compat      !checks compatibility of messages send/received
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
!
  private
!
! PRIVATE MEMBER FUNCTIONS:
  private :: csm_allocate         !csm dynamic memory allocation
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

  real(r8) :: send2d(lsmlon,lsmlat,nsend_csm)  !2d send buffer
  real(r8) :: recv2d(lsmlon,lsmlat,nrecv_csm)  !2d recv buffer
  
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
! !IROUTINE: csm_recvorb
!
! !INTERFACE: 
  subroutine csm_recvorb (eccen, obliqr, lambm0, mvelpp)
!
! !DESCRIPTION: 
!  Receive initial integer, real and character control data from 
!  the flux coupler.  Parse it into the variables used by the land
!  model. Receive information contains orbital parameters from the 
!  flux coupler.
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(out) :: eccen  !Earth's eccentricity of orbit
    real(r8), intent(out) :: obliqr !Earth's obliquity in radians
    real(r8), intent(out) :: lambm0 !Mean longitude of perihelion at the vernal equinox (radians)
    real(r8), intent(out) :: mvelpp !Earth's moving vernal equinox long of perihelion plus pi (radians)
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: i,j,k                !indices
    integer  :: ierr                 !error code                     
    integer  :: info_time            !T => turn on msg-passing timing
    integer  :: maj_vers             !Major version of message passed from the cpl
    integer  :: min_vers             !Minor version of message passed from the cpl
    real(r8) :: spval                !float-pt buffer special value from coupler
!------------------------------------------------------------------------

    if (masterproc) then

       ! Receive control message from coupler

       ibuffr(:) = 0    

       call shr_msg_recv_i (ibuffr, size(ibuffr), SHR_MSG_TID_CPL, SHR_MSG_TAG_C2LI)

       ierr       = ibuffr( 1)  !error code
       info_time  = ibuffr(11)  !T => turn on msg-passing timing
       maj_vers   = ibuffr(40)  !Coupler message major version
       min_vers   = ibuffr(41)  !Coupler message minor version
       ncbuff     = ibuffr(42)  !Size of character data to receive
       write(6,*) '(CSM_RECVORB): recd d->l initial ibuffr msg_id = ',SHR_MSG_TAG_C2LI

       ! Determine debug flag

       if (ibuffr(12) >= 2) then	
          debug_flag = .true.
       else
          debug_flag = .false.
       endif

       ! Check that the version of the message from the coupler is valid

       call csm_compat(maj_vers, min_vers,SHR_MSG_L_MAJ_V04, SHR_MSG_L_MIN_V00)

       ! Receive orbital parameters from coupler

       rbuff(:) = 0.0

       call shr_msg_recv_r (rbuff, size(rbuff), SHR_MSG_TID_CPL, SHR_MSG_TAG_C2LI)

       spval  = rbuff(1)      !Coupler float-pt special data flag
       eccen  = rbuff(2)      !Earth's eccentricity of orbit
       obliqr = rbuff(3)      !Earth's Obliquity radians
       lambm0 = rbuff(4)      !Earth's Long. of prehelian at v-equinox
       mvelpp = rbuff(5)      !Earth's Moving vernal equinox of orbit + pi
       write(6,*)'(CSM_RECVORB): recd d->l initial real buff msg_id = ',SHR_MSG_TAG_C2LI
 
       ! Check that data received is good data and not the special value
 
       call compat_check_spval(spval, eccen ,'Eccentricity'     )
       call compat_check_spval(spval, obliqr,'Obliquity'        )
       call compat_check_spval(spval, lambm0,'Long of perhelion')
       call compat_check_spval(spval, mvelpp,'Move long of perh')
       
       write(6,*)'(CSM_RECVORB): eccen:  ', eccen
       write(6,*)'(CSM_RECVORB): obliqr: ', obliqr
       write(6,*)'(CSM_RECVORB): lambm0: ', lambm0
       write(6,*)'(CSM_RECVORB): mvelpp: ', mvelpp
       
       ! Receive character data message from coupler and determine if will output csm timing info

       if (ncbuff > 0) then
          call shr_msg_recv_c (cbuff, ncbuff, SHR_MSG_TID_CPL, SHR_MSG_TAG_C2LI)
          write(6,*)'(CSM_RECVORB): recd d->a initial char. buf msg_id = ',SHR_MSG_TAG_C2LI
          write(6,*)'(CSM_RECVORB): Char: ',(cbuff(i), i = 1,ncbuff)
       end if

       if (info_time == 0) then
          csm_timing = .false.
       else
          csm_timing = .true.
       endif

    end if
#if ( defined SPMD )
   call mpi_bcast(spval , 1, MPI_REAL8, 0, mpicom, ierr)
   if (ierr/= MPI_SUCCESS) then
      write(6,*) 'MPI BCAST ERROR: for spval in csm_recvorb'
      call endrun
   end if
   call mpi_bcast(eccen , 1, MPI_REAL8, 0, mpicom, ierr)
   if (ierr /= MPI_SUCCESS) then
      write(6,*) 'MPI BCAST ERROR: for eccen in csm_recvorb'
      call endrun
   end if
   call mpi_bcast(obliqr, 1, MPI_REAL8, 0, mpicom, ierr)
   if (ierr /= MPI_SUCCESS) then
      write(6,*) 'MPI BCAST ERROR: for obliqr in csm_recvorb'
      call endrun
   end if
   call mpi_bcast(lambm0, 1, MPI_REAL8, 0, mpicom, ierr)
   if (ierr /= MPI_SUCCESS) then
      write(6,*) 'MPI BCAST ERROR: for lambm0 in csm_recvorb'
      call endrun
   end if
   call mpi_bcast(mvelpp, 1, MPI_REAL8, 0, mpicom, ierr)
   if (ierr /= MPI_SUCCESS) then
      write(6,*) 'MPI BCAST ERROR: for mvelpp in csm_recvorb'
      call endrun
   end if
#endif

   return
 end subroutine csm_recvorb 

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_sendcontrol
!
! !INTERFACE:
 subroutine csm_sendcontrol(irad)
!
! !DESCRIPTION: 
!  Send first control data to flux coupler and "invalid" grid
!  containing special value data.
!  The coupler treats points where the mask is nonzero as points where 
!  you could possibly do a calculation (in the case of the runoff, this 
!  corresponds to all rtm ocean points). The coupler then defines a "key"
!  as points where the model can give you valid data (in the case of runoff, 
!  this corresponds to points where the land model will give you valid 
!  compressed data points). The key can be 0 where the mask is 1. However,
!  the key cannot be 1 where the mask is 0 unless the data is also zero. 
!  In the case of runoff, the key the coupler builds is time invariant.
!  Send first control data to flux coupler and "invalid" grid
!  containing special value data
!
! !USES:
    use clm_varctl, only : csm_doflxave, nsrest
    use RtmMod, only : area_r, longxy_r, latixy_r, mask_r
    use clm_varcon, only : re	
    use time_manager, only : get_step_size
    use shr_const_mod, only : SHR_CONST_CDAY
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: irad  !frequency of radiation computation 
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
!
! !LOCAL VARIABLES
    real(r8) rtemp_lnd(lsmlon*lsmlat*4) !temporary vector
    integer  itemp_lnd(lsmlon*lsmlat)   !temporary vector 
    real(r8) rtemp_rtm(rtmlon*rtmlat*4) !temporary vector
    real(r8) temp_area(rtmlon,rtmlat)   !temporary area 
    integer  temp_mask(rtmlon,rtmlat)   !temporary mask
    real(r8) dtime                      !step size
!------------------------------------------------------------------------

    if (masterproc) then

       rtemp_lnd(:) = 1.e30
       itemp_lnd(:) = 999
       rtemp_rtm(:) = 1.e30

       ! determine number of send/recv calls steps per day to flux coupler
    
       if (nsrest == 0) then
          dtime = get_step_size()
       else
          dtime = csm_dtime
       endif
       if (csm_doflxave) then
          ncpday = nint(SHR_CONST_CDAY/dtime)/irad
       else
          ncpday = nint(SHR_CONST_CDAY/dtime)
       endif
       
       ! send integer control information 

       ibuffs(:)  = 0                !initialize ibuffs
       ibuffs(7)  = lsmlon           !number of land longitudes
       ibuffs(8)  = lsmlat           !number of land latitudes
       ibuffs(9)  = ncpday           !number of land send/recv calls per day
       ibuffs(34) = 1                !T(or 1) => requests cpl to send valid land domain info back
       ibuffs(36) = 0                !size of compressed runoff vector, if zero then not sending compressed info 
       ibuffs(37) = rtmlon           !number of longitudes in uncompressed 2d runoff array
       ibuffs(38) = rtmlat           !number of latitudes  in uncompressed 2d runoff array

       call shr_msg_send_i (ibuffs ,nibuff, SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)

       ! send "invalid" land model grid and mask data

       call shr_msg_send_r (rtemp_lnd(1:lsmlon*lsmlat  ), lsmlon*lsmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)
       call shr_msg_send_r (rtemp_lnd(1:lsmlon*lsmlat  ), lsmlon*lsmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)
       call shr_msg_send_r (rtemp_lnd(1:lsmlon*lsmlat*4), lsmlon*lsmlat*4, SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)
       call shr_msg_send_r (rtemp_lnd(1:lsmlon*lsmlat*4), lsmlon*lsmlat*4, SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)
       call shr_msg_send_r (rtemp_lnd(1:lsmlon*lsmlat  ), lsmlon*lsmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)
       call shr_msg_send_i (itemp_lnd(1:lsmlon*lsmlat  ), lsmlon*lsmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)

       ! send "valid" RTM grid and mask data 
      
       temp_area(:,:) = area_r(:,:)/(re*re)  !convert from km^2 to radians^2 before sending to coupler
       temp_mask(:,:) = 1 - mask_r(:,:)      !make coupler runoff mask 1 over ocean and 0 over land 

       call shr_msg_send_r (longxy_r , rtmlon*rtmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)
       call shr_msg_send_r (latixy_r , rtmlon*rtmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)
       call shr_msg_send_r (rtemp_rtm, rtmlon*rtmlat*4, SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)
       call shr_msg_send_r (rtemp_rtm, rtmlon*rtmlat*4, SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)
       call shr_msg_send_r (temp_area, rtmlon*rtmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)
       call shr_msg_send_i (temp_mask, rtmlon*rtmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)

       write(6,*)'(CSM_SENDCONTROL): there will be ',ncpday, &
            ' send/recv calls per day from the land model to the flux coupler'
       write(6,*)'(CSM_SENDCONTROL):sent l->d control data msg_id = ',SHR_MSG_TAG_L2CI

    endif

    return
  end subroutine csm_sendcontrol

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_recvgrid
!
! !INTERFACE:
  subroutine csm_recvgrid (cam_longxy, cam_latixy, cam_numlon, cam_landfrac, cam_landmask) 
!
! !DESCRIPTION: 
!  Receive valid land grid and land mask from coupler
!
! !ARGUMENTS:
    implicit none
    integer , intent(out) :: cam_numlon(lsmlat)           !cam number of longitudes 
    real(r8), intent(out) :: cam_longxy(lsmlon,lsmlat)    !cam lon values
    real(r8), intent(out) :: cam_latixy(lsmlon,lsmlat)    !cam lat values 
    real(r8), intent(out) :: cam_landfrac(lsmlon,lsmlat)  !cam fractional land
    integer , intent(out) :: cam_landmask(lsmlon,lsmlat)  !cam land mask
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
!
! !LOCAL VARIABLES:
    integer  i,j                      !loop indices
    real(r8) xe(4,lsmlon,lsmlat)      !coupler land grid edges 
    real(r8) ye(4,lsmlon,lsmlat)      !coupler land grid edges
    real(r8) area_a(lsmlon,lsmlat)    !coupler atm grid areas
    integer  mask_a(lsmlon,lsmlat)    !coupler atm valid grid mask 
!------------------------------------------------------------------------

    if (masterproc) then

       ibuffr(:) = 0
       
       call shr_msg_recv_i (ibuffr       ,nibuff         , SHR_MSG_TID_CPL, SHR_MSG_TAG_C2LI)
       call shr_msg_recv_r (cam_longxy   ,lsmlon*lsmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_C2LI)
       call shr_msg_recv_r (cam_latixy   ,lsmlon*lsmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_C2LI)
       call shr_msg_recv_r (xe           ,lsmlon*lsmlat*4, SHR_MSG_TID_CPL, SHR_MSG_TAG_C2LI)
       call shr_msg_recv_r (ye           ,lsmlon*lsmlat*4, SHR_MSG_TID_CPL, SHR_MSG_TAG_C2LI)
       call shr_msg_recv_r (area_a       ,lsmlon*lsmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_C2LI)
       call shr_msg_recv_r (cam_landfrac ,lsmlon*lsmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_C2LI)
       call shr_msg_recv_i (cam_landmask ,lsmlon*lsmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_C2LI)
       call shr_msg_recv_i (mask_a       ,lsmlon*lsmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_C2LI)
       
       write(6,*)'(CSM_SENDGRID):recd d->l land grid, msg_id= ',SHR_MSG_TAG_C2L
       
       ! USE mask_a to determine number of valid longitudes for each latitude band 
       ! this is the only use for mask_a
       
       cam_numlon(:) = 0
       do j = 1,lsmlat
          do i= 1,lsmlon
             if (mask_a(i,j) /= 0) cam_numlon(j) = cam_numlon(j)+1
          end do
       end do
       
    endif  !end of if-masterproc block
    
    return
  end subroutine csm_recvgrid

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_sendrunoff
!
! !INTERFACE:
  subroutine csm_sendrunoff()
!
! !DESCRIPTION: 
! Send valid runoff information back to flux coupler
!
! !USES:
    use RtmMod, only : ocnrof_iindx, ocnrof_jindx, ocnrof_vec
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
! -----------------------------------------------------------------

    if (masterproc) then

       ! Send integer buffer control info

       ibuffs(7)  = lsmlon           !number of lsm longitudes
       ibuffs(8)  = lsmlat           !number of lsm latitudes
       ibuffs(36) = size(ocnrof_vec) !number of data points in compressed runoff data 
       ibuffs(37) = rtmlon           !number of longitudes in uncompressed 2d runoff array
       ibuffs(38) = rtmlat           !number of latitudes  in uncompressed 2d runoff array

       call shr_msg_send_i(ibuffs       ,nibuff          , SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)

       ! Send runoff vector compression info

       call shr_msg_send_i(ocnrof_iindx, size(ocnrof_vec), SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)
       call shr_msg_send_i(ocnrof_jindx, size(ocnrof_vec), SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)

       write(6,*) '(CSM_SENDROF):sent l->d valid initial runoff info msg_id = ',SHR_MSG_TAG_L2CI

    endif

    return
  end subroutine csm_sendrunoff

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_allocate
!
! !INTERFACE:
  subroutine csm_allocate ()
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

       allocate (send1d_glob(nsend_csm,numg), STAT=ier)
       call allocation_err(ier, 'csm_allocate', 'send1d_glob', numg)
       send1d_glob(:,:) = nan
    endif

#if (defined SPMD)
    allocate (send1d_loc(nsend_csm,begg:endg), STAT=ier)
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
  end subroutine csm_allocate

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_sendalb
!
! !INTERFACE:
  subroutine csm_sendalb()
!
! !DESCRIPTION: 
! Send initial albedos, surface temperature and snow data to the 
! flux coupler
!
! !USES:
    use clmtype
    use clmpoint
    use clm_varsur                    
    use clm_varctl  , only : csm_doflxave, nsrest
    use time_manager, only : get_curr_date, get_prev_date
    use RtmMod      , only : ocnrof_vec
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
    integer  :: i,j,n       !loop indices
    integer  :: gi,li,ci,pi !indices
    real(r8) :: wt          !weight  
    integer  :: yr          !current year 
    integer  :: mon         !current month 
    integer  :: day         !current day (0, 1, ...)
    integer  :: ncsec       !current seconds of current date (0, ..., 86400)
    integer  :: ncdate      !current date (yymmdd format) (e.g., 021105)
    type(gridcell_type)       , pointer :: g   !local pointer to derived subtype
    type(landunit_type)       , pointer :: l   !local pointer to derived subtype
    type(column_type)         , pointer :: c   !local pointer to derived subtype
    type(pft_type)            , pointer :: p   !local pointer to derived subtype
    type(gridcell_pstate_type), pointer :: gps !local pointer to derived subtype
    type(landunit_pstate_type), pointer :: lps !local pointer to derived subtype
    type(column_pstate_type)  , pointer :: cps !local pointer to derived subtype
    type(column_estate_type)  , pointer :: ces !local pointer to derived subtype 
    type(column_wstate_type)  , pointer :: cws !local pointer to derived subtype
    type(pft_pstate_type)     , pointer :: pps !local pointer to derived subtype
! -----------------------------------------------------------------

    ! Allocate csm dynamic memory

    call csm_allocate()
    
    ! Send albedos

    if (nsrest == 0 ) then   !initial run

       ! On initial timestep ONLY: determine 1d vector of states that will be sent
       ! to coupler and map fields from 1d subgrid vector to 2d [lsmlon]x[lsmlat] grid.

       send1d_loc(:,:) = 0.
       do ci = cols1d%beg,cols1d%end
          c => cpoint%col(ci)%c
          cps => c%cps
          cws => c%cws
          ces => c%ces
          gi = cols1d%gindex(ci)
          send1d_loc(isend_trad,gi) = send1d_loc(isend_trad,gi) + ces%t_grnd * cps%wtxy
          send1d_loc(isend_sno ,gi) = send1d_loc(isend_sno,gi) + cws%h2osno/1000. * cps%wtxy 
          do pi = 1, cps%npfts
             p => c%p(pi)
             pps => p%pps
             send1d_loc(isend_asdir,gi) = send1d_loc(isend_asdir,gi) + pps%albd(1) * pps%wtxy  
             send1d_loc(isend_aldir,gi) = send1d_loc(isend_aldir,gi) + pps%albd(2) * pps%wtxy 
             send1d_loc(isend_asdif,gi) = send1d_loc(isend_asdif,gi) + pps%albi(1) * pps%wtxy 
             send1d_loc(isend_aldif,gi) = send1d_loc(isend_aldif,gi) + pps%albi(2) * pps%wtxy 
          end do
       end do

#if (defined SPMD)
       call gather_data_to_master(send1d_loc, send1d_glob, clmlevel=grid1d%name)
#endif

       if (masterproc ) then
          do n=1,nsend_csm
             where (landmask(:,:) > 0) 
                send2d(:,:,n) = 0.
             elsewhere
                send2d(:,:,n) = 1.e30
             end where
          end do
          send2d(:,:,isend_sno) = 0. ! snow initialized to 0 everywhere
          do gi = 1, grid1d%num
             i = grid1d%ixy(gi)
             j = grid1d%jxy(gi)
             send2d(i,j,:nsend_csm) = send1d_glob(:nsend_csm,gi)
          end do
       endif

    else  ! restart run

       ! On a restart run, no meaningful data is sent to the flux coupler - 
       ! this includes ocnrof_vec (which should only contain zero values)
       ! since the runoff code (riverfluxrtm) has not been called yet
       
       if (masterproc) send2d(:,:,:) = 1.e30  ! this will be sent on a restart timestep

    endif  
      
    ! Send data to coupler
    ! Determine time index to send to coupler. Note that for a restart run,
    ! the next time step is nstep+1. But must send current time step to flux couper here.

    if (masterproc) then

       if (nsrest == 0) then
          call get_curr_date (yr, mon, day, ncsec) 
       else
          call get_prev_date (yr, mon, day, ncsec)
       endif
       ncdate = yr*10000 + mon*100 + day

       ibuffs(4) = ncdate  !model date (yyyymmdd)
       ibuffs(5) = ncsec   !elapsed seconds in model date

       call shr_msg_send_i (ibuffs    , size(ibuffs)    , SHR_MSG_TID_CPL, SHR_MSG_TAG_L2C)
       call shr_msg_send_r (send2d    , size(send2d)    , SHR_MSG_TID_CPL, SHR_MSG_TAG_L2C)
       call shr_msg_send_r (ocnrof_vec, size(ocnrof_vec), SHR_MSG_TID_CPL, SHR_MSG_TAG_L2C)

       if (csm_timing) then
          irtc_s = shr_sys_irtc()
          write(6,9099) irtc_s,'l->d sending'
9099      format('[mp timing]  irtc = ',i20,' ',a)
       end if

    endif  ! end of if_masterproc

    return
  end subroutine csm_sendalb

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_dosndrecv
!
! !INTERFACE:
  subroutine csm_dosndrcv (doalb)
!
! !DESCRIPTION: 
! Determine when to send and receive messages to/from the
! flux coupler on this time-step.
! Determine if send/receive information to/from flux coupler
! Send msgs (land state and fluxes) to the flux coupler only when 
! doalb is true (i.e. on time steps before the atm does a solar
! radiation computation). Receive msgs (atm state) from the
! flux coupler only when dorad is true (i.e. on time steps 
! when the atm does a solar radiation computation). 
! The fluxes are then averaged between the send and receive calls. 
! 
! !USES:
    use clm_varctl   , only : csm_doflxave
    use time_manager , only : get_step_size, get_nstep
    use shr_const_mod, only : SHR_CONST_CDAY
!
! !ARGUMENTS:
    implicit none
    logical, intent(in) :: doalb  !true=>next timestep a radiation time step
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: ntspday           !model steps per day
    real(r8) :: dtime             !step size (seconds)   
    integer  :: nstep             !time step 
!-----------------------------------------------------------------------

    ! Determine if send/receive information to/from flux coupler

    nstep = get_nstep()
    if (csm_doflxave) then
       if (nstep == 0) then
          dorecv = .true.
          dosend = .false.
       else if (nstep == 1) then
          dorecv = .false.
          dosend = doalb
       else
          dorecv = dosend
          dosend = doalb
       endif
    else
       if (nstep == 0) then
          dorecv = .true.
          dosend = .false.
       else if (nstep == 1) then
          dorecv = .false.
          dosend = .true.
       else
          dorecv = .true.
          dosend = .true.
       endif
    endif

    ! If at end of day: check if should write restart file or stop 
    ! at next time step. Note, these statements must appear here since 
    ! ibuffr is not received at every time step when flux averaging occurs.

    csmstop_next = .false.
    csmrstrt     = .false.
    dtime        = get_step_size()
    ntspday      = nint(SHR_CONST_CDAY/dtime)
    if (mod(nstep,ntspday) == 0) then
       if (ibuffr(2) /= 0) then  !stop at end of day
          csmstop_next = .true.  !will stop on next time step
          write(6,*)'(CSM_DOSNDRCV) output restart and history files at nstep = ',nstep
       endif
       if (ibuffr(21) /= 0) then !write restart at end of day
          csmrstrt = .true.      !will write restart now
          write(6,*)'(CSM_DOSNDRCV) output restart and history files at nstep = ',nstep
       endif
    endif

    return
  end subroutine csm_dosndrcv

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_recv
! 
! !INTERFACE:
  subroutine csm_recv()
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
        call shr_msg_recv_i (ibuffr, size(ibuffr), SHR_MSG_TID_CPL, SHR_MSG_TAG_C2L)
        call shr_msg_recv_r (recv2d, size(recv2d), SHR_MSG_TID_CPL, SHR_MSG_TAG_C2L)
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
  end subroutine csm_recv

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_send
!
! !INTERFACE:
  subroutine csm_send()
!   	
! !DESCRIPTION: 
! Send data to the flux coupler
!
! !USES:
    use clmtype
    use clmpoint
    use clm_varctl, only : csm_doflxave
    use clm_varsur, only : landmask
    use time_manager, only : get_curr_date
    use RtmMod, only : ocnrof_vec
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
    integer :: i,j,n,gi,ci,pi   !indices
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
    type(pft_pstate_type)     , pointer :: pps !pointer to derived subtype
    type(pft_mflux_type)      , pointer :: pmf !pointer to derived subtype
    type(pft_wflux_type)      , pointer :: pwf !pointer to derived subtype
    type(pft_eflux_type)      , pointer :: pef !pointer to derived subtype 
    type(pft_estate_type)     , pointer :: pes !pointer to derived subtype
!-----------------------------------------------------------------------

    ! Send data to the flux coupler

    if (timer_lnd_recvsend) then
       call t_stopf ('lnd_recvsend') ; timer_lnd_recvsend = .false.
    endif

    ! Start timer

    call t_startf('lnd_send')

    ! Determine 1d vector of fields that will be sent to coupler.
    ! Coupler has convention that fluxes are positive downward.

    send1d_loc(:,:) = 0.

    do ci = cols1d%beg,cols1d%end
       c => cpoint%col(ci)%c
       cps => c%cps
       cws => c%cws
       gi = cols1d%gindex(ci)
       send1d_loc(isend_sno,gi) = send1d_loc(6,gi) + cws%h2osno/1000. * cps%wtxy 
    end do

    do pi = pfts1d%beg,pfts1d%end
       p => ppoint%pft(pi)%p
       pps => p%pps
       pes => p%pes
       gi = pfts1d%gindex(pi)
       send1d_loc(isend_trad  ,gi) = send1d_loc(isend_trad ,gi) + pes%t_rad_pft * pps%wtxy 
       send1d_loc(isend_asdir ,gi) = send1d_loc(isend_asdir,gi) + pps%albd(1) * pps%wtxy  
       send1d_loc(isend_aldir ,gi) = send1d_loc(isend_aldir,gi) + pps%albd(2) * pps%wtxy 
       send1d_loc(isend_asdif ,gi) = send1d_loc(isend_asdif,gi) + pps%albi(1) * pps%wtxy 
       send1d_loc(isend_aldif ,gi) = send1d_loc(isend_aldif,gi) + pps%albi(2) * pps%wtxy 
       send1d_loc(isend_tref2m,gi) = send1d_loc(isend_tref2m,gi) + pes%t_ref2m * pps%wtxy 
       if (csm_doflxave) then
          send1d_loc(isend_taux ,gi) = send1d_loc(isend_taux ,gi) - taux_ave(pi) * pps%wtxy
          send1d_loc(isend_tauy ,gi) = send1d_loc(isend_tauy ,gi) - tauy_ave(pi) * pps%wtxy
          send1d_loc(isend_lhflx,gi) = send1d_loc(isend_lhflx,gi) - lhflx_ave(pi) * pps%wtxy
          send1d_loc(isend_shflx,gi) = send1d_loc(isend_shflx,gi) - shflx_ave(pi) * pps%wtxy
          send1d_loc(isend_lwup ,gi) = send1d_loc(isend_lwup ,gi) - lwup_ave(pi) * pps%wtxy
          send1d_loc(isend_qflx ,gi) = send1d_loc(isend_qflx ,gi) - qflx_ave(pi) * pps%wtxy
          send1d_loc(isend_swabs,gi) = send1d_loc(isend_swabs,gi) - swabs_ave(pi) * pps%wtxy
       else
          pef => p%pef
          pmf => p%pmf
          pwf => p%pwf
          send1d_loc(isend_taux ,gi) = send1d_loc(isend_taux ,gi) - pmf%taux * pps%wtxy 
          send1d_loc(isend_tauy ,gi) = send1d_loc(isend_tauy ,gi) - pmf%tauy * pps%wtxy 
          send1d_loc(isend_lhflx,gi) = send1d_loc(isend_lhflx,gi) - pef%eflx_lh_tot * pps%wtxy 
          send1d_loc(isend_shflx,gi) = send1d_loc(isend_shflx,gi) - pef%eflx_sh_tot * pps%wtxy 
          send1d_loc(isend_lwup ,gi) = send1d_loc(isend_lwup ,gi) - pef%eflx_lwrad_out * pps%wtxy 
          send1d_loc(isend_qflx ,gi) = send1d_loc(isend_qflx ,gi) - pwf%qflx_evap_tot * pps%wtxy 
          send1d_loc(isend_swabs,gi) = send1d_loc(isend_swabs,gi) - pef%fsa * pps%wtxy 
       endif
    end do

#if (defined SPMD)
    call gather_data_to_master(send1d_loc, send1d_glob, clmlevel=grid1d%name)
#endif
    
    ! Map fields from 1d-grid vector to 2d-grid
    ! NOTE: snow is sent as zero over non-land because currently 
    ! the ocn and sea-ice send no snow cover to cpl and so the cpl 
    ! sends back zero snow over non-land to  the atm (the atm and 
    ! land grid are currently assumed to be identical)
    
    if (masterproc) then

       do n = 1,nsend_csm
          where( landmask(:,:) > 0 ) 
             send2d(:,:,n) = 0.
          elsewhere
             send2d(:,:,n) = 1.e30
          endwhere
       end do
       send2d(:,:,isend_sno) = 0.     !reset snow to 0 everywhere
       do gi = 1, grid1d%num
          i = grid1d%ixy(gi)
          j = grid1d%jxy(gi)
          send2d(i,j,:nsend_csm) = send1d_glob(:nsend_csm,gi)
       end do
       
       call get_curr_date (yr, mon, day, ncsec) 
       ncdate = yr*10000 + mon*100 + day

       ibuffs(4)  = ncdate           !model date (yyyymmdd)
       ibuffs(5)  = ncsec            !elapsed seconds in current date
       
       call shr_msg_send_i (ibuffs    , size(ibuffs)    , SHR_MSG_TID_CPL, SHR_MSG_TAG_L2C)
       call shr_msg_send_r (send2d    , size(send2d)    , SHR_MSG_TID_CPL, SHR_MSG_TAG_L2C)
       call shr_msg_send_r (ocnrof_vec, size(ocnrof_vec), SHR_MSG_TID_CPL, SHR_MSG_TAG_L2C)

       if (csm_timing) then
          irtc_s = shr_sys_irtc()
          write(6,9099) irtc_s,'l->d sending'
9099      format('[mp timing]  irtc = ',i20,' ',a)
       end if
       
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

    endif  ! end of if_masterproc

    ! Stop timers

    call t_stopf('lnd_send')

    if (.not. timer_lnd_recvsend) then
       call t_startf('lnd_sendrecv') ; timer_lnd_sendrecv = .true.
    endif

    return
  end subroutine csm_send

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_flxave
!
! !INTERFACE:
  subroutine csm_flxave() 
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
  end subroutine csm_flxave

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: compat_check
!
! !INTERFACE:
  subroutine compat_check_spval( spval, data, string )
!   	
! !DESCRIPTION: 
! Check that the given piece of real data sent from the coupler is valid
! data and not the couplers special data flag.  This ensures that the data
! you expect is actually being sent by the coupler.
!
! !PARAMETERS
    implicit none
    real(r8), intent(in) :: spval
    real(r8), intent(in) :: data
    character(len=*), intent(in) :: string
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
! -----------------------------------------------------------------

    if ( spval == data )then
       write(6,*)'ERROR:(compat_check_spval) msg incompatibility'
       write(6,*)'ERROR: I expect to recieve the data type: ',string
       write(6,*)'from CPL, but all I got was the special data flag'
       write(6,*)'coupler must not be sending this data, you are'
       write(6,*)'running with an incompatable version of the coupler'
       call endrun
    end if
    return
  end subroutine compat_check_spval
  
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: subroutine csm_compat
!
! !INTERFACE:
  subroutine csm_compat(cpl_maj_vers, cpl_min_vers, expect_maj_vers, expect_min_vers)
!   	
! !DESCRIPTION: 
! Checks that the message recieved from the coupler is compatable
! with the type of message that I expect to recieve.  If the minor
! version numbers differ I print a warning message.  If the major
! numbers differ I abort since that means that the change is
! drastic enough that I can't run with the differences.
! Original Author: Erik Kluzek Dec/97
!
! !PARAMETERS
    implicit none
    integer, intent(in) :: cpl_maj_vers    ! major version from coupler initial ibuffr array
    integer, intent(in) :: cpl_min_vers    ! minor version from coupler initial ibuffr array
    integer, intent(in) :: expect_maj_vers ! major version of the coupler I'm expecting
    integer, intent(in) :: expect_min_vers ! minor version of the coupler I'm expecting
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
! -----------------------------------------------------------------

    write(6,*)'(cpl_COMPAT): This is revision: $Revision$'
    write(6,*)'              Tag: $Name$'
    write(6,*)'              of the message compatability interface:'

    if ( cpl_min_vers /= expect_min_vers )then
       write(6,*) 'WARNING(cpl_compat):: Minor version of coupler ', &
            'messages different than expected: '
       write(6,*) 'The version of the coupler being used is: ',&
            cpl_min_vers
       write(6,*) 'The version I expect is:                  ',&
            expect_min_vers
    end if
    
    if ( cpl_maj_vers /= expect_maj_vers )then
       write(6,*) 'ERROR(cpl_compat):: Major version of coupler ', &
            'messages different than expected: '
       write(6,*) 'The version of the coupler being used is: ',&
            cpl_maj_vers
       write(6,*) 'The version I expect is:                  ',&
            expect_maj_vers
       call endrun
    end if
    
    return
  end subroutine csm_compat

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restart_coupler
!
! !INTERFACE:
  subroutine restart_coupler (nio, flag)
!
! !DESCRIPTION: 
!  Read/write restart data needed for running in flux coupled mode
!
! !USES:
    use clm_varctl, only : csm_doflxave
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: nio             !restart unit 
    character(len=*), intent(in) :: flag   !"read" or "write"
!
! !Revision History:
!  02.09.17  Mariana Vertenstein: moved code to be part of ccsm module
!
!EOP
!
! !LOCAL VARIABLES:
    logical :: flxave_res   !flux averaging flag read from restart file
    integer :: ier          !mpi return error code
!-----------------------------------------------------------------------

    if (flag == 'read') then
       if (masterproc) then
          read(nio) flxave_res
          read(nio) dosend
          if ((flxave_res .and. .not.csm_doflxave).or.(.not.flxave_res .and. csm_doflxave)) then
             write(6,*)'(RESTART_COUPLER): flxave value from namelist ',csm_doflxave, &
                  ' must be the same as flxave value from restart dataset ',flxave_res
             call endrun
          endif
          if (flxave_res .and. .not. dosend) then
             write(6,*)'(RESTART_COUPLER): assume that current flux coupled model ', &
                  'with flux averaging must stop on a time step where dosend (doalb) is true'
             call endrun
          end if
       endif
#if ( defined SPMD ) 
       call mpi_bcast(dosend    , 1, MPI_LOGICAL, 0, mpicom, ier)
       if (ier /= MPI_SUCCESS) then
          write(6,*) 'MPI BCAST ERROR: for dosend in restart_coupler'
          call endrun
       end if
       call mpi_bcast(flxave_res, 1, MPI_INTEGER, 0, mpicom, ier)
       if (ier /= MPI_SUCCESS) then
          write(6,*) 'MPI BCAST ERROR: for flxave_res in restart_coupler'
          call endrun
       end if
#endif
    endif

    if (flag == 'write') then
       if (masterproc) then
          write(nio) csm_doflxave
          write(nio) dosend
       endif
    end if

    return
end subroutine restart_coupler

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

end module clm_csmMod
