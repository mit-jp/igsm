#include <misc.h>
#include <preproc.h>

module restFileMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: restFileMod
! 
! !DESCRIPTION: 
! Reads from or  writes to/ the CLM restart file. 
! 
! !USES
  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmdMod, only : masterproc
  use shr_sys_mod, only : shr_sys_flush
!
! !PUBLIC TYPES:
  implicit none
! save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: restart
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: restart_setup               !Setup restart files, perform consistency checks
  private :: restart_time                !Read/write restart time manager data
  private :: restart_biogeophys          !Read/write restart biogeophysics data
  private :: restart_wrapup              !Close restart file and write restart pointer file 
  private :: write_rest_pfile            !Writes restart pointer file
  private :: set_restart_filename        !Sets restart filename

!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
! 
! !PRIVATE TYPES
  private
  character(len=256) :: locfnr           !restart file name 
  integer, parameter :: rest_id = 6      !restart id
!
!----------------------------------------------------------------------- 

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restart
!
! !INTERFACE:
  subroutine restart (flag)
!
! !DESCRIPTION: 
! Read CLM restart file.
!
! !USES:
    use clm_varctl, only : nsrest, rpntdir, rpntfil, nrevsn, & 
                           archive_dir, mss_irt, mss_wpass, &
                           csm_doflxave, caseid
#if (defined RTM)
    use RtmMod, only : restart_rtm
#endif
    use histFileMod, only : restart_history
#if (defined COUP_CSM)
    use clm_csmMod, only : restart_coupler
#endif
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    character(len=*), intent(in) :: flag   !'read' or 'write'
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: nio                         !Fortran unit number
!-----------------------------------------------------------------------

    if (flag /= 'read' .and. flag  /= 'write') then 
       write(6,*)'ERROR: flag must be set to read or write'
       call endrun
    endif

    call restart_setup(nio, flag)

    call restart_time(nio, flag)
          
    call restart_biogeophys(nio, flag)
    
#if (defined COUP_CSM)
    call restart_coupler(nio, flag)
#endif

#if (defined RTM)
    call restart_rtm(nio, flag)   
#endif

    call restart_history(nio, flag)

    call restart_wrapup(nio, flag)

  end subroutine restart

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restart_setup
!
! !INTERFACE:
  subroutine restart_setup (nio, flag)
!
! !DESCRIPTION: 
    ! Setup restart file and perform necessary consistency checks
!
! !USES:
    use fileutils, only : opnfil, getfil, getavu, relavu
    use time_manager, only : get_step_size, get_nstep 
    use clm_varctl, only : nsrest, nrevsn, rpntdir, rpntfil, caseid 
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    integer, intent(out) :: nio             !restart unit 
    character(len=*), intent(in) :: flag   !'read' or 'write'
!
! !CALLED FROM:
! subroutine restart in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer i,n                            !indices
    character(len=256) :: fnamer           !full name of restart file
    character(len= 32) :: casename         !case name read in from restart
!-----------------------------------------------------------------------

    if (masterproc) then

       if (flag == 'write') then
          ! Open restart file for writing
          write (6,*) 'Writout out restart data at nstep = ',get_nstep()
          write (6,'(72a1)') ("-",i=1,60)
          locfnr = set_restart_filename()
          nio = getavu()
          call opnfil (locfnr, nio, 'u')
       else 
          ! Obtain the restart file for reading
          ! For restart runs, the restart pointer file contains the 
          ! full pathname of the restart file. For branch runs, the 
          ! namelist variable [nrevsn] contains the full pathname of the
          ! restart file.  New history files are created for branch runs.
          write (6,*) 'Reading in restart data .....'
          write (6,'(72a1)') ("-",i=1,60)
          if (nsrest == 1) then     !restart
             nio = getavu()
             locfnr = trim(rpntdir) //'/'// trim(rpntfil)
             call opnfil (locfnr, nio, 'f')
             read (nio,'(a80)') fnamer
             call relavu (nio)
          else                      !branch run
             fnamer = nrevsn
          end if
          nio = getavu()
          call getfil (fnamer, locfnr, 0)
          call opnfil (locfnr, nio, 'u')
       endif

       if (flag == 'write') then
          ! Write restart id
          write(nio) rest_id
       else
          ! Check restart id (restart file id must match restart code version id)
          read(nio) n
          if (n == rest_id) then
             write(6,*)'RESTRD: using restart file id ',rest_id
          else
             write(6,*)'RESTRD: input restart file id ',n, &
                  ' does not match required restart id ',rest_id
             call endrun
          endif
       endif
       
       if (flag == 'write') then
          ! Write case name
          write(nio) caseid
       else
          ! Check case name consistency (case name must be different for branch run)
          read (nio) casename
          if (nsrest == 3 .and. casename==caseid) then
             write(6,*) 'RESTRD ERROR: Must change case name on branch run'
             write(6,*) 'Prev case = ',trim(casename),' current case = ',trim(caseid)
             call endrun
          end if
       endif

    endif  !masterproc if-block
          
  end subroutine restart_setup

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restart_time
!
! !INTERFACE:
  subroutine restart_time (nio, flag)
!
! !DESCRIPTION: 
! Read/write time manager information to/from restart
!
! !USES:
    use globals
#if (defined COUP_CSM)
    use controlMod, only : csm_dtime       !dtime from input namelist
#endif
    use time_manager, only : get_step_size, timemgr_write_restart, &
                             timemgr_read_restart, timemgr_restart, get_nstep
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: nio             !restart unit 
    character(len=*), intent(in) :: flag   !'read' or 'write'
!
! !CALLED FROM:
! subroutine restart in this module
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
#if (defined COUP_CAM)
    integer :: clm_nstep                     !nstep from restart file
#endif
!-----------------------------------------------------------------------

#if (defined OFFLINE)

    if (flag == 'write') then
       if (masterproc) call timemgr_write_restart(nio)
    else
       if (masterproc) call timemgr_read_restart(nio)
       call timemgr_restart()
    end if

#elif (defined COUP_CAM)    

    if (masterproc) then
       if (flag == 'write') then
          ! Write out time step

          write(nio) get_nstep()
       else
          ! Read in time step - check that clm restart time step is the
          ! same as cam restart time step. In cam mode, the time manager
          ! restart variables have already been read in by calls to 
          ! routines timemgr_read_restart and timemgr_restart.

          read(nio) clm_nstep
          if ((clm_nstep+1) /= get_nstep()) then
             write(6,*)'RESTART_TIME: incompatibility in clm and cam restart dates'
             write(6,*)'  restart step from cam = ',get_nstep()
             write(6,*)'  restart step from clm = ',clm_nstep + 1
             call endrun
          endif
       endif
    endif

#elif (defined COUP_CSM)

    if (flag == 'write') then
       if (masterproc) call timemgr_write_restart(nio)
    else
       if (masterproc) call timemgr_read_restart(nio)
       call timemgr_restart()
    end if
    if (masterproc) then
       if (csm_dtime /= get_step_size()) then
          write(6,*)'(RESTRD): error '
          write(6,*)'namelist dtime on restart does not match input dtime'
          call endrun
       endif
    endif

#endif

    ! Set model time step. (On an initial run, this is set in routine iniTimeConst.F90) 

    dtime = get_step_size()

  end subroutine restart_time

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restart_biogeophys
!
! !INTERFACE:
  subroutine restart_biogeophys (nio, flag)
!
! !DESCRIPTION: 
! Read/Write biogeophysics information to/from restart file. 
!
! !USES:
    use clmtype
    use clmpoint
    use iobinary
    use clm_varpar, only : nlevsoi, numrad, lsmlon, lsmlat
    use clm_varcon, only : denice, denh2o
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: nio             !restart unit 
    character(len=*), intent(in) :: flag  !'read' or 'write'
!
! !CALLED FROM:
! subroutine restart in this module
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer gi,li,ci,pi,index,j     !indices
    real(r8):: pftsum               !temporary used for pft averaging for columns
    type(gridcell_type)  , pointer :: g     !local pointer to derived subtype
    type(landunit_type)  , pointer :: l     !local pointer to derived subtype
    type(column_type)    , pointer :: c     !local pointer to derived subtype
    type(pft_type)       , pointer :: p     !local pointer to derived subtype
    integer , pointer, dimension(:) :: ibuf1dg
    integer , pointer, dimension(:) :: ibuf1dl
    integer , pointer, dimension(:) :: ibuf1dc
    integer , pointer, dimension(:) :: ibuf1dp
    real(r8), pointer, dimension(:) :: rbuf1dg
    real(r8), pointer, dimension(:) :: rbuf1dl
    real(r8), pointer, dimension(:) :: rbuf1dc
    real(r8), pointer, dimension(:) :: rbuf1dp
    real(r8), pointer, dimension(:,:) :: rbuf2dc
    real(r8), pointer, dimension(:,:) :: rbuf2dp
    integer begc,endc,numc          !column 1d indices
    integer begp,endp,nump          !pft 1d indices
    character(len=16) namec         !column 1d name 
    character(len=16) namep         !pft 1d name 
!-----------------------------------------------------------------------

    ! Set up shorthand for derived type names

    begc = cols1d%beg
    endc = cols1d%end
    numc = cols1d%num
    namec = cols1d%name

    begp = pfts1d%beg
    endp = pfts1d%end
    nump = pfts1d%num
    namep = pfts1d%name

   ! Allocate necessary 1d buffers

    allocate (ibuf1dg(grid1d%num))
    allocate (ibuf1dl(land1d%num))
    allocate (ibuf1dc(cols1d%num))
    allocate (ibuf1dp(pfts1d%num))

    allocate (rbuf1dg(grid1d%num))
    allocate (rbuf1dl(land1d%num))
    allocate (rbuf1dc(cols1d%num))
    allocate (rbuf1dp(pfts1d%num))

    ! Column physical state - snl
    if (flag == 'read') call readin (nio, ibuf1dc, clmlevel=cols1d%name)  
    do ci = begc,endc
       c => cpoint%col(ci)%c
       if (flag == 'read' ) c%cps%snl = ibuf1dc(ci)
       if (flag == 'write') ibuf1dc(ci) = c%cps%snl
    end do
    if (flag == 'write') call wrtout (nio, ibuf1dc, clmlevel=cols1d%name)  

    ! Column physical state - snowdp
    if (flag == 'read') call readin (nio, rbuf1dc, clmlevel=cols1d%name)  
    do ci = begc,endc
       c => cpoint%col(ci)%c
       if (flag == 'read' ) c%cps%snowdp = rbuf1dc(ci)
       if (flag == 'write') rbuf1dc(ci) = c%cps%snowdp  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dc, clmlevel=cols1d%name)  

    ! Column physical state - snowage
    if (flag == 'read') call readin (nio, rbuf1dc, clmlevel=cols1d%name)  
    do ci = begc,endc
       c => cpoint%col(ci)%c
       if (flag == 'read' ) c%cps%snowage = rbuf1dc(ci)
       if (flag == 'write') rbuf1dc(ci) = c%cps%snowage  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dc, clmlevel=cols1d%name)  

    ! Column physical state - frac_sno
    if (flag == 'read') call readin (nio, rbuf1dc, clmlevel=cols1d%name)  
    do ci = begc,endc
       c => cpoint%col(ci)%c
       if (flag == 'read' ) c%cps%frac_sno = rbuf1dc(ci)
       if (flag == 'write') rbuf1dc(ci) = c%cps%frac_sno 
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dc, clmlevel=cols1d%name)  

    ! Column physical state - dz (snow)
    allocate (rbuf2dc(-nlevsno+1:0,numc))
    if (flag == 'read') call readin (nio, rbuf2dc, clmlevel=cols1d%name)  
    do ci = begc,endc
       c => cpoint%col(ci)%c
       do j = -nlevsno+1,0
          if (flag == 'read' ) c%cps%dz(j) = rbuf2dc(j,ci)
          if (flag == 'write') rbuf2dc(j,ci) = c%cps%dz(j) 
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dc, clmlevel=cols1d%name)  
    deallocate(rbuf2dc)

    ! Column physical state - z (snow)
    allocate (rbuf2dc(-nlevsno+1:0,numc))
    if (flag == 'read') call readin (nio, rbuf2dc, clmlevel=cols1d%name)  
    do ci = begc,endc
       c => cpoint%col(ci)%c
       do j = -nlevsno+1,0
          if (flag == 'read' ) c%cps%z(j) = rbuf2dc(j,ci)
          if (flag == 'write') rbuf2dc(j,ci) = c%cps%z(j) 
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dc, clmlevel=cols1d%name)  
    deallocate(rbuf2dc)

    allocate (rbuf2dc(-nlevsno:0,numc))
    if (flag == 'read') call readin (nio, rbuf2dc, clmlevel=cols1d%name)  
    do ci = begc,endc
       c => cpoint%col(ci)%c
       do j = -nlevsno,0
          if (flag == 'read' ) c%cps%zi(j) = rbuf2dc(j,ci)
          if (flag == 'write') rbuf2dc(j,ci) = c%cps%zi(j) 
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dc, clmlevel=cols1d%name)  
    deallocate(rbuf2dc)

    !pft type physical state variable - albd
    allocate (rbuf2dp(numrad,nump))
    if (flag == 'read') call readin (nio, rbuf2dp, clmlevel=pfts1d%name)  
    do pi = begp,endp
       p => ppoint%pft(pi)%p
       do j = 1,numrad
          if (flag == 'read' ) p%pps%albd(j) = rbuf2dp(j,pi)
          if (flag == 'write') rbuf2dp(j,pi) = p%pps%albd(j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dp, clmlevel=pfts1d%name)  
    deallocate(rbuf2dp)

    !pft type physical state variable - albi
    allocate (rbuf2dp(numrad,nump))
    if (flag == 'read') call readin (nio, rbuf2dp, clmlevel=pfts1d%name)  
    do pi = begp,endp
       p => ppoint%pft(pi)%p
       do j = 1,numrad
          if (flag == 'read' ) p%pps%albi(j) = rbuf2dp(j,pi)
          if (flag == 'write') rbuf2dp(j,pi) = p%pps%albi(j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dp, clmlevel=pfts1d%name)  
    deallocate(rbuf2dp)

    !column type physical state variable - albgrd
    allocate (rbuf2dc(numrad,numc))
    if (flag == 'read') call readin (nio, rbuf2dc, clmlevel=cols1d%name)  
    do ci = begc,endc
       c => cpoint%col(ci)%c
       do j = 1,numrad
          if (flag == 'read' ) c%cps%albgrd(j) = rbuf2dc(j,ci)
          if (flag == 'write') rbuf2dc(j,ci) = c%cps%albgrd(j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dc, clmlevel=cols1d%name)  
    deallocate(rbuf2dc)

    !column type physical state variable - albgri
    allocate (rbuf2dc(numrad,numc))
    if (flag == 'read') call readin (nio, rbuf2dc, clmlevel=cols1d%name)  
    do ci = begc,endc
       c => cpoint%col(ci)%c
       do j = 1,numrad
          if (flag == 'read' ) c%cps%albgri(j) = rbuf2dc(j,ci)
          if (flag == 'write') rbuf2dc(j,ci) = c%cps%albgri(j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dc, clmlevel=cols1d%name)  
    deallocate(rbuf2dc)

   ! column water state variable - h2osno
    if (flag == 'read') call readin (nio, rbuf1dc, clmlevel=cols1d%name)  
    do gi = 1,clm%mps%ngridcells
       do li = 1,clm%g(gi)%gps%nlandunits
          do ci = 1,clm%g(gi)%l(li)%lps%ncolumns
             c => clm%g(gi)%l(li)%c(ci)
             index = c%cps%index1d
             if (flag == 'read' ) c%cws%h2osno = rbuf1dc(index)
             if (flag == 'write') rbuf1dc(index) = c%cws%h2osno
          end do
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dc, clmlevel=cols1d%name)  

   ! column water state variable - h2osoi_liq
    allocate (rbuf2dc(-nlevsno+1:nlevsoi,numc))
    if (flag == 'read') call readin (nio, rbuf2dc, clmlevel=cols1d%name)  
    do ci = begc,endc
       c => cpoint%col(ci)%c
       do j = -nlevsno+1,nlevsoi
          if (flag == 'read' ) c%cws%h2osoi_liq(j) = rbuf2dc(j,ci)
          if (flag == 'write') rbuf2dc(j,ci) = c%cws%h2osoi_liq(j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dc, clmlevel=cols1d%name)  
    deallocate(rbuf2dc)

   ! column water state variable - h2osoi_ice
    allocate (rbuf2dc(-nlevsno+1:nlevsoi,numc))
    if (flag == 'read') call readin (nio, rbuf2dc, clmlevel=cols1d%name)  
    do ci = begc,endc
       c => cpoint%col(ci)%c
       do j = -nlevsno+1,nlevsoi
          if (flag == 'read' ) c%cws%h2osoi_ice(j) = rbuf2dc(j,ci)
          if (flag == 'write') rbuf2dc(j,ci) = c%cws%h2osoi_ice(j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dc, clmlevel=cols1d%name)  
    deallocate(rbuf2dc)

   ! column energy state variable - t_grnd
    if (flag == 'read') call readin (nio, rbuf1dc, clmlevel=cols1d%name)  
    do ci = begc,endc
       c => cpoint%col(ci)%c
       if (flag == 'read' ) c%ces%t_grnd = rbuf1dc(ci)
       if (flag == 'write') rbuf1dc(ci) = c%ces%t_grnd  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dc, clmlevel=cols1d%name)  

   ! column energy state variable - t_soisno
    allocate (rbuf2dc(-nlevsno+1:nlevsoi,numc))
    if (flag == 'read') call readin (nio, rbuf2dc, clmlevel=cols1d%name)  
    do ci = begc,endc
       c => cpoint%col(ci)%c
       do j = -nlevsno+1,nlevsoi
          if (flag == 'read' ) c%ces%t_soisno(j) = rbuf2dc(j,ci)
          if (flag == 'write') rbuf2dc(j,ci) = c%ces%t_soisno(j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dc, clmlevel=cols1d%name)  
    deallocate(rbuf2dc)

    !column type energy state variable - t_lake
    allocate (rbuf2dc(1:nlevlak,numc))
    if (flag == 'read') call readin (nio, rbuf2dc, clmlevel=cols1d%name)  
    do ci = begc,endc
       c => cpoint%col(ci)%c
       do j = 1,nlevlak
          if (flag == 'read' ) c%ces%t_lake(j) = rbuf2dc(j,ci)
          if (flag == 'write') rbuf2dc(j,ci) = c%ces%t_lake(j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dc, clmlevel=cols1d%name)  
    deallocate(rbuf2dc)

    ! pft type physical state variable - frac_veg_nosno_alb 
    if (flag == 'read') call readin (nio, ibuf1dp, clmlevel=pfts1d%name)  
    do pi = begp,endp
       p => ppoint%pft(pi)%p
       if (flag == 'read' ) p%pps%frac_veg_nosno_alb = ibuf1dp(pi)
       if (flag == 'write') ibuf1dp(pi) = p%pps%frac_veg_nosno_alb 
    end do
    if (flag == 'write') call wrtout (nio, ibuf1dp, clmlevel=pfts1d%name)

    ! pft type physical state variable - fwet
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=pfts1d%name)  
    do pi = begp,endp
       p => ppoint%pft(pi)%p
       if (flag == 'read' ) p%pps%fwet = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pps%fwet 
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=pfts1d%name)  

    ! pft type physical state variable - tlai
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=pfts1d%name)  
    do pi = begp,endp
       p => ppoint%pft(pi)%p
       if (flag == 'read' ) p%pps%tlai = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pps%tlai 
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=pfts1d%name)  

    ! pft type physical state variable - tsai
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=pfts1d%name)  
    do pi = begp,endp
       p => ppoint%pft(pi)%p
       if (flag == 'read' ) p%pps%tsai = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pps%tsai 
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=pfts1d%name)  

    ! pft type physical state variable - elai
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=pfts1d%name)  
    do pi = begp,endp
       p => ppoint%pft(pi)%p
       if (flag == 'read' ) p%pps%elai = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pps%elai  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=pfts1d%name)  

    ! pft type physical state variable - esai
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=pfts1d%name)  
    do pi = begp,endp
       p => ppoint%pft(pi)%p
       if (flag == 'read' ) p%pps%esai = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi)= p%pps%esai 
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=pfts1d%name)  

    ! pft type physical state variable - fsun
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=pfts1d%name)  
    do pi = begp,endp
       p => ppoint%pft(pi)%p
       if (flag == 'read') then
          p%pps%fsun = rbuf1dp(pi)
       else
          rbuf1dp(pi)= p%pps%fsun
       end if
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=pfts1d%name)  

    ! pft type physical state variable - htop
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=pfts1d%name)  
    do pi = begp,endp
       p => ppoint%pft(pi)%p
       if (flag == 'read' ) p%pps%htop = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi)= p%pps%htop
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=pfts1d%name)  

    ! pft type physical state variable - hbot
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=pfts1d%name)  
    do pi = begp,endp
       p => ppoint%pft(pi)%p
       if (flag == 'read' ) p%pps%hbot = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi)= p%pps%hbot
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=pfts1d%name)  

    ! pft type physical state variable - fabd
    allocate (rbuf2dp(numrad,nump))
    if (flag == 'read') call readin (nio, rbuf2dp, clmlevel=pfts1d%name)  
    do pi = begp,endp
       p => ppoint%pft(pi)%p
       do j = 1,numrad
          if (flag == 'read' ) p%pps%fabd(j) = rbuf2dp(j,pi)
          if (flag == 'write') rbuf2dp(j,pi) = p%pps%fabd(j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dp, clmlevel=pfts1d%name)  
    deallocate(rbuf2dp)

    ! pft type physical state variable - fabi
    allocate (rbuf2dp(numrad,nump))
    if (flag == 'read') call readin (nio, rbuf2dp, clmlevel=pfts1d%name)  
    do pi = begp,endp
       p => ppoint%pft(pi)%p
       do j = 1,numrad
          if (flag == 'read' ) p%pps%fabi(j) = rbuf2dp(j,pi)
          if (flag == 'write') rbuf2dp(j,pi) = p%pps%fabi(j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dp, clmlevel=pfts1d%name)  
    deallocate(rbuf2dp)

    ! pft type physical state variable - ftdd
    allocate (rbuf2dp(numrad,nump))
    if (flag == 'read') call readin (nio, rbuf2dp, clmlevel=pfts1d%name)  
    do pi = begp,endp
       p => ppoint%pft(pi)%p
       do j = 1,numrad
          if (flag == 'read' ) p%pps%ftdd(j) = rbuf2dp(j,pi)
          if (flag == 'write') rbuf2dp(j,pi) = p%pps%ftdd(j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dp, clmlevel=pfts1d%name)  
    deallocate(rbuf2dp)

    ! pft type physical state variable - ftid
    allocate (rbuf2dp(numrad,nump))
    if (flag == 'read') call readin (nio, rbuf2dp, clmlevel=pfts1d%name)  
    do pi = begp,endp
       p => ppoint%pft(pi)%p
       do j = 1,numrad
          if (flag == 'read' ) p%pps%ftid(j) = rbuf2dp(j,pi)
          if (flag == 'write') rbuf2dp(j,pi) = p%pps%ftid(j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dp, clmlevel=pfts1d%name)  
    deallocate(rbuf2dp)

    ! pft type physical state variable - ftii
    allocate (rbuf2dp(numrad,nump))
    if (flag == 'read') call readin (nio, rbuf2dp, clmlevel=pfts1d%name)  
    do pi = begp,endp
       p => ppoint%pft(pi)%p
       do j = 1,numrad
          if (flag == 'read' ) p%pps%ftii(j) = rbuf2dp(j,pi)
          if (flag == 'write') rbuf2dp(j,pi) = p%pps%ftii(j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dp, clmlevel=pfts1d%name)  
    deallocate(rbuf2dp)

    ! pft type energy state variable - t_veg
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=pfts1d%name)  
    do pi = begp,endp
       p => ppoint%pft(pi)%p
       if (flag == 'read' ) p%pes%t_veg = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pes%t_veg
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=pfts1d%name)  

    ! pft type water state variable - h2ocan 
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=pfts1d%name)  
    do pi = begp,endp
       p => ppoint%pft(pi)%p
       if (flag == 'read' ) p%pws%h2ocan = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pws%h2ocan   
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=pfts1d%name)  

    ! For read only:
    ! Determine average over all column pfts for h2ocan, needed by begwb 
    ! computation in routine driver.F90)
    if (flag == 'read' ) then
       do ci = begc,endc
          c => cpoint%col(ci)%c
          pftsum = 0.0
          do pi = 1,c%cps%npfts
             pftsum = pftsum + c%p(pi)%pws%h2ocan * c%pw(pi)
          end do
          c%cws%pws_a%h2ocan = pftsum
       end do
    endif

    ! For read only:
    ! Determine volumetric soil water
    if (flag == 'read' ) then
       do ci = begc,endc
          c => cpoint%col(ci)%c
          do j = 1,nlevsoi
             c%cws%h2osoi_vol(j) = c%cws%h2osoi_liq(j)/(c%cps%dz(j)*denh2o) &
                                 + c%cws%h2osoi_ice(j)/(c%cps%dz(j)*denice)
          end do
       end do
    endif
    
    deallocate (ibuf1dg)
    deallocate (ibuf1dl)
    deallocate (ibuf1dc)
    deallocate (ibuf1dp)
                     
    deallocate (rbuf1dg)
    deallocate (rbuf1dl)
    deallocate (rbuf1dc)
    deallocate (rbuf1dp)

  end subroutine restart_biogeophys

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restart_wrapup
!
! !INTERFACE:
  subroutine restart_wrapup (nio, flag)       
!
! !DESCRIPTION: 
! Close and archive restart file and write restart pointer file if
! in write mode, otherwise just close restart file if in read mode
!
! !USES:
    use clm_varctl, only : mss_irt, mss_wpass, archive_dir
#if (defined COUP_CSM)
    use clm_csmMod, only : csmstop_next
#endif
    use fileutils, only : putfil, relavu, set_filename 
    use time_manager, only : is_last_step
!
! !ARGUMENTS:
    implicit none
    integer, intent(inout) :: nio          !restart unit 
    character(len=*), intent(in) :: flag   !'read' or 'write'
!
! !CALLED FROM:
! subroutine restart in this module
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer i                      !index
    logical :: lremove             !true => remove file after archive
    character(len=256) :: rem_fn   !remote (archive) filename
    character(len=256) :: rem_dir  !remote (archive) directory
!-----------------------------------------------------------------------

   if (masterproc) then
      if (flag == 'write') then
         ! close and archive restart file and write restart pointer file
         lremove = .true.
#if (defined OFFLINE) || (defined COUP_CAM)
         if (is_last_step()) lremove = .false.
#elif (defined COUP_CSM)
         if (csmstop_next) lremove = .false.
#endif
         call relavu (nio)
         write(6,*) 'Successfully wrote local restart file ',trim(locfnr)
         if (mss_irt > 0) then
            rem_dir = trim(archive_dir) // '/rest/'
            rem_fn = set_filename(rem_dir, locfnr)
            call putfil (locfnr, rem_fn, mss_wpass, mss_irt, lremove)
         endif
         call write_rest_pfile( )
         write(6,'(72a1)') ("-",i=1,60)
         write(6,*)
      else
         ! close file
         call relavu (nio)
         write(6,'(72a1)') ("-",i=1,60)
         write(6,*) 'Successfully read restart data for restart run'
         write(6,*)
      endif
   endif
   
   return
  end subroutine restart_wrapup
 
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_rest_pfile
!
! !INTERFACE:
  subroutine write_rest_pfile ( )
!
! !DESCRIPTION: 
! Open restart pointer file. Write names of current restart and 
! history files. If using mass store, these are the mass store
! names except if mss_irt=0 (no mass store files written). Close. 
!
! !USES:
    use clm_varctl, only : rpntdir, mss_irt, archive_dir, rpntfil
    use fileutils, only : set_filename, relavu
    use fileutils, only : getavu, opnfil
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine restart in this module
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: m                           !index
    integer :: nio                         !restart pointer file 
    character(len=256) :: filename         !local file name
    character(len=256) :: rem_dir          !remote directory
!-----------------------------------------------------------------------

    nio = getavu()
    filename= trim(rpntdir) //'/'// trim(rpntfil)
    call opnfil (filename, nio, 'f')

    ! write name of restart file to pointer file 
    if (mss_irt == 0) then
       write(nio,'(a)') locfnr
    else
       rem_dir = trim(archive_dir) // '/rest/'
       write(nio,'(a)') set_filename(rem_dir, locfnr)
    endif

    ! add comments to pointer file of all files that are needed for restart
    write(nio,*)'The following lines list files needed for restart - do not edit'

    ! only write names of open history files
!    do m = 1,nhist
!       if (locfnh(m) /= ' ') then
!          if (mss_irt == 0) then
!             filename = locfnh(m)
!          else
!             rem_dir = trim(archive_dir) // '/hist/'
!             filename = set_filename(rem_dir, locfnh(m))
!          endif
!          write(nio, '(a)') filename
!       end if
!    end do

    call relavu (nio)
    write(6,*)'Successfully wrote local restart pointer file'

  end subroutine write_rest_pfile

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_restart_filename
!
! !INTERFACE:
  character(len=256) function set_restart_filename ()
!
! !DESCRIPTION: 
!
! !USES:
    use clm_varctl  , only : caseid
    use time_manager, only : get_curr_date
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine restart in this module
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=256) :: cdate       !date char string
    integer :: day                    !day (1 -> 31)
    integer :: mon                    !month (1 -> 12)
    integer :: yr                     !year (0 -> ...)
    integer :: sec                    !seconds into current day
!-----------------------------------------------------------------------
#if (defined ONE_RESTART)
    set_restart_filename = "./"//trim(caseid)//".clm2.rst"
#else
    call get_curr_date (yr, mon, day, sec) 
    write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,sec
    set_restart_filename = "./"//trim(caseid)//".clm2.r."//trim(cdate)
#endif

    return
  end function set_restart_filename

end module restFileMod



