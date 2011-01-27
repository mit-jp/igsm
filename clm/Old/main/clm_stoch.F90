      subroutine stoch(pcpc,pcpl,pcpcpoiss,pcplpoiss, &
                       storm,exptr,expta,trout,taout, &
                       dtcum,strmnum,pcpcacc,pcplacc)
      use time_manager, only : get_step_size
      implicit none
      real*8 :: pcpl
      real*8 :: pcpc
      integer,parameter :: nranmx=10000
      integer,save :: iexpon
      integer :: storm
      integer,save :: tr
      integer,save :: ta
      integer :: dtcum
      integer :: strmnum
      integer :: irn
      real,save :: expon(nranmx)
      real,save :: dt
      real :: pcpcacc
      real :: pcplacc
      real :: pcpcpoiss
      real :: pcplpoiss
      real :: trout
      real :: taout
      real :: expta
      real :: exptr
      logical,save :: start
      data start/.true./
      data tr/0/,ta/0/,iexpon/0/
!
!  This was adapted from the file relfns.f that was used for Milly (WRR,
!  30(7), 1994). 
!
!       initialize constants (flags)
!
      if (start) then
        dt=get_step_size()
        call random(expon,nranmx,iexpon)
        start=.false.
        write (6,*) 'Initializing Stochastic Forcing Routine'
      else
        tr=int(trout/dt)
        ta=int(taout/dt)
      endif
!      write (6,*) dtcum,dt,tr,ta
!
!	check whether at beginning, within or after pcp Poisson-box event
!
     if (dtcum==0) then
!
!	Pick ta and tr, begin accumulating precip for next Poisson box event
!
!	Calculate the pcp rate of current Poisson event, accumulate event counter
!
        pcpcpoiss=pcpcacc/float(tr)
        pcplpoiss=pcplacc/float(tr)
        storm=tr
        pcpcacc=0.
        pcplacc=0.
        iexpon = iexpon + 1
        if (iexpon.gt.nranmx) call random(expon,nranmx,iexpon)
!
!	Round storm duration to the nearest whole number of timesteps
!	If tr < model timestep, set tr = model timestep
!
        tr=0
          tr = int(((exptr/dt)*expon(iexpon))+0.5)
          if ( tr < 1 ) tr = 1
          iexpon = iexpon + 1
          if (iexpon.gt.nranmx) call random(expon,nranmx,iexpon)
!
!	Round inter-storm period to nearest whole number of timesteps
!	Assure that the storm duration is shorter than or equal to inter-storm arrival period
!
        ta=0
        do while (ta<=tr)
          ta = int(((expta/dt)*expon(iexpon))+0.5)
          iexpon = iexpon + 1
          if (iexpon.gt.nranmx) call random(expon,nranmx,iexpon)
        enddo
        taout=ta*dt
        trout=tr*dt
!        write (6,10) ta,tr
!10      format ('Inter storm period: ',i8,' timesteps;  Storm duration ',i8,' timesteps')
      endif
      dtcum=dtcum+1
      if (dtcum==ta) then
        dtcum=0
        strmnum=strmnum+1
      endif
      if (dtcum>storm) then
        pcpcpoiss=0.
        pcplpoiss=0.
      endif
      pcpcacc=pcpcacc+pcpc
      pcplacc=pcplacc+pcpl
      pcpc=pcpcpoiss
      pcpl=pcplpoiss
      return
      end
!===============================================================================
      subroutine random (expon,nranmx,iexpon)
      real expon(nranmx)
       call random_number (harvest=expon)
       do irn = 1, nranmx
        do while (expon(irn)==0.or.expon(irn)==1) 
          expon(irn)=ran(irn)
        enddo
        expon(irn) = -alog(expon(irn))
        if (expon(irn)==0) expon(irn)=-alog(ran(nranmx))
       enddo
       iexpon=0
       return
      end
!===============================================================================
