      subroutine stoch(pcpc,pcpl,pcpcpoiss,pcplpoiss, &
                       storm,exptr,expta,trout,taout, &
                       dtcum,strmnum,pcpcacc,pcplacc, &
                       pcprnd1,pcprnd2,prnt_stoch)
      use clm_time_manager, only : get_step_size
      implicit none
      real*8 :: pcpl
      real*8 :: pcpc
      integer,parameter :: nranmx=10000
      integer :: storm
      integer,save :: tr
      integer,save :: ta
      integer :: dtcum
      integer :: strmnum
      integer :: irn, iter
      real,save :: dt
      real :: expon
      real*8 :: pcpcacc
      real*8 :: pcplacc
      real :: pcpcpoiss
      real :: pcplpoiss
      real*8 :: trout
      real*8 :: taout
      real :: expta
      real :: exptr
      real*8 :: pcprnd1
      real*8 :: pcprnd2
      logical,save :: start
      logical :: prnt_stoch
      data start/.true./
!
!  This was adapted from the file relfns.f that was used for Milly (WRR,
!  30(7), 1994). 
!
!       initialize constants (flags)
!
     if (start) then
       dt=get_step_size()
!
!  Initialize accumlators
!
       start=.false.
       write (6,*) 'Initializing Stochastic Forcing Routine'
     endif
     tr=int(trout/dt)
     ta=int(taout/dt)
!      write (6,*) dtcum,dt,tr,ta
!      write (6,*) pcprnd1, pcprnd2
!
!	check whether at beginning, within or after pcp Poisson-box event
!
     if (dtcum==0) then
!
!	Pick ta and tr, begin accumulating precip for next Poisson box event
!
!	Calculate the pcp rate of current Poisson event, accumulate event counter
!
       storm=tr
       pcpcacc=0.
       pcplacc=0.
!
!	Round storm duration to the nearest whole number of timesteps
!	If tr < model timestep, set tr = model timestep
!
       tr=0
       expon = -alog(pcprnd1)
       tr = int(((exptr/dt)*expon)+0.5)
       if ( tr < 1 ) tr = 1
!
!	Round inter-storm period to nearest whole number of timesteps
!
       ta=0
       iter = 1
!      do while (ta<=tr)
          expon = -alog(pcprnd2)
          ta = int(((expta/dt)*expon)+0.5)
          if ( ta < 1 ) ta = 1
!      enddo
       taout=ta*dt
       trout=tr*dt
     if (prnt_stoch ) then
         write (6,50) dtcum
        write (6,20) pcprnd1,pcprnd2
 20      format ('Rand1 storm period: ',e20.12,'  Rand2 duration ',e20.12)
         write (6,10) ta,tr
 10      format ('Inter storm period: ',i8,' timesteps;  Storm duration ',i8,' timesteps')
      endif
!       write (6,70) taout,trout
!70      format ('Inter storm period: ',e20.12,'  Storm duration ',e20.12)
     endif
     dtcum=dtcum+1
     if (prnt_stoch ) then
         write (6,50) dtcum
        write (6,30) pcpcacc,pcplacc
     endif
     if (dtcum<=ta) then
        pcpcacc=pcpcacc+pcpc
        pcplacc=pcplacc+pcpl
        pcpcpoiss=0.
        pcplpoiss=0.
     elseif(dtcum <ta+tr) then
       pcpcpoiss=pcpcacc/float(tr)+pcpc
       pcplpoiss=pcplacc/float(tr)+pcpl
     else
       pcpcpoiss=pcpcacc/float(tr)+pcpc
       pcplpoiss=pcplacc/float(tr)+pcpl
        dtcum=0
        strmnum=strmnum+1
     endif
     if (prnt_stoch ) then
         write (6,10) ta,tr
        write (6,60) pcpc,pcpl
 60      format (' IGSM pcpc= ',e20.12,'  IGSM pcpl= ',e20.12)
        write (6,30) pcpcacc,pcplacc
 30      format (' pcpcacc= ',e20.12,'  pcplacc= ',e20.12)
        write (6,40) pcpcpoiss,pcplpoiss
 40      format (' pcpcpoiss= ',e20.12,'  pcplpoiss= ',e20.12)
         write (6,50) dtcum
 50      format ('DTCUM= ',i8,' timesteps;')
     endif
     pcpc=pcpcpoiss
     pcpl=pcplpoiss
     return
     end
