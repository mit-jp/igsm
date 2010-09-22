C $Header: /home/igsm-repository/CHEM/chemlsodesrc.F,v 1.1 2010/09/02 17:56:19 jscott Exp $
C $Name:  $

c===============================================================c
c	LSODESRC:	Collection of subroutines or functions  c	
c			used in Chien Wang's  photochemistry	c
c			model from ODEPACK			c
c	------------------------------------------------	c
c	Chien Wang						c
c	MIT, Rm.E40-269, Cambridge, MA 02139			c
c								c
c	December 5, 1995					c
c	------------------------------------------------	c
c	Table of Contents					c
c								c
c		lsodenew:	subroutine	-- TC1		c
c		stodenew:	subroutine	-- TC2		c
c		ewset:		subroutine	-- TC3		c
c		solsy:		subroutine	-- TC4		c
c		intdy:		subroutine	-- TC5		c
c		cfode:		subroutine	-- TC6		c
c		sgefa:		subroutine	-- TC7		c
c		sgbfa: 		subroutine	-- TC8		c
c		sgesl:		subroutine	-- TC9		c
c		sgbsl:		subroutine	-- TC10		c
c		sscal:		subroutine	-- TC11		c
c		saxpysmp:	subroutine	-- TC12		c
c		xerrwv:		subroutine	-- TC13		c
c		r1mach:		function	-- TC14		c
c		vnorm:		function	-- TC15		c
c		isamax:		function	-- TC16		c
c		sdot:		function	-- TC17		c
c		prepj64:	subroutine	-- TC18		c
c===============================================================c

!	------------------------------------------------------------
!
!	Revision History:
!	
!	When	Who		What
!	----	----------	-------	
!	092001	Chien Wang	rev. of some old style code
!
!	==========================================================

c===============================================================
c -- TC1
c
      subroutine lsodenew (f, neq, y, t, tout, itol, rtol, atol, itask,
     &            istate, iopt, rwork, lrw, iwork, liw, jac, mf)
c     ==========================================================

c===============================================================c
c   LSODENEW.F:	A simplified version of original lsode.f    	c
c 		for cases where					c
c 		ISTATE = 1 &					c
c 		ITASK  = 1 initially				c
c		IOPT   = 0					c
c		ITOL   = 1					c
c        ------------------------------------------------       c
c								c
c  		Chien Wang					c
c		MIT Joint Program for Science and Policy 	c
c			of Global Change			c
c								c
c		March 20, 1995					c
c===============================================================c

      external f, jac
      integer neq, itol, itask, istate, iopt, lrw, iwork, liw, mf
      real y, t, tout, rtol, atol, rwork
      dimension neq(*), y(*), rtol(*), atol(*), rwork(lrw), iwork(liw)

c-----------------------------------------------------------------------
c this is the march 30, 1987 version of
c lsode.. livermore solver for ordinary differential equations.
c this version is in single precision.
c
c lsode solves the initial value problem for stiff or nonstiff
c systems of first order ode-s,
c     dy/dt = f(t,y) ,  or, in component form,
c     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(neq)) (i = 1,...,neq).
c lsode is a package based on the gear and gearb packages, and on the
c october 23, 1978 version of the tentative odepack user interface
c standard, with minor modifications.
c-----------------------------------------------------------------------
c reference..
c     alan c. hindmarsh,  odepack, a systematized collection of ode
c     solvers, in scientific computing, r. s. stepleman et al. (eds.),
c     north-holland, amsterdam, 1983, pp. 55-64.
c-----------------------------------------------------------------------
c author and contact.. alan c. hindmarsh,
c                      computing and mathematics research div., l-316
c                      lawrence livermore national laboratory
c                      livermore, ca 94550.
c-----------------------------------------------------------------------

      external prepj, solsy
      integer illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     1   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns
      integer icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     1   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, i1, i2, iflag, imxer, kgo, lf0,
     1   leniw, lenrw, lenwm, ml, mord, mu, mxhnl0, mxstp0
      real rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real atoli, ayi, big, ewti, h0, hmax, hmx, rh, rtoli,
     1   tcrit, tdist, tnext, tol, tolsf, tp, size, sum, w0,
     2   r1mach, vnorm
      dimension mord(2)
      logical ihit
c-----------------------------------------------------------------------
c the following internal common block contains
c (a) variables which are local to any subroutine but whose values must
c     be preserved between calls to the routine (own variables), and
c (b) variables which are communicated between subroutines.
c the structure of the block is as follows..  all real variables are
c listed first, followed by all integers.  within each type, the
c variables are grouped with those local to subroutine lsode first,
c then those local to subroutine stode, and finally those used
c for communication.  the block is declared in subroutines
c lsode, intdy, stode, prepj, and solsy.  groups of variables are
c replaced by dummy arrays in the common declarations in routines
c where those variables are not used.
c-----------------------------------------------------------------------
      common /ls0001/ rowns(209),
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     2   illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     3   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
c
c      data  mord(1),mord(2)/12,5/, mxstp0/500/, mxhnl0/10/
      data  mord(1),mord(2)/12,5/, mxstp0/260/, mxhnl0/10/
      data illin/0/, ntrep/0/
c-----------------------------------------------------------------------
c block a.
c this code block is executed on every call.
c it tests istate and itask for legality and branches appropriately.
c if istate .gt. 1 but the flag init shows that initialization has
c not yet been done, an error return occurs.
c if istate = 1 and tout = t, jump to block g and return immediately.
c-----------------------------------------------------------------------

      init = 0
      ntrep = 0
c-----------------------------------------------------------------------
c block b.
c the next code block is executed for the initial call (istate = 1),
c or for a continuation call with parameter changes (istate = 3).
c it contains checking of all inputs and various initializations.
c
c first check legality of the non-optional inputs neq, itol, iopt,
c mf, ml, and mu.
c-----------------------------------------------------------------------

      n = neq(1)
      meth = mf/10
      miter = mf - 10*meth
      if (miter .gt. 3)then
        ml = iwork(1)
        mu = iwork(2)
      endif

c next process and check the optional inputs. --------------------------
      maxord = mord(meth)
      mxstep = mxstp0
      mxhnil = mxhnl0
      h0 = 0.0e0
      hmxi = 0.0e0
      hmin = 0.0e0

c-----------------------------------------------------------------------
c set work array pointers and check lengths lrw and liw.
c pointers to segments of rwork and iwork are named by prefixing l to
c the name of the segment.  e.g., the segment yh starts at rwork(lyh).
c segments of rwork (in order) are denoted  yh, wm, ewt, savf, acor.
c-----------------------------------------------------------------------

      lyh = 21
      if (istate .eq. 1) nyh = n
      lwm = lyh + (maxord + 1)*nyh
      if (miter .eq. 0) lenwm = 0
      if (miter .eq. 1 .or. miter .eq. 2) lenwm = n*n + 2
      if (miter .eq. 3) lenwm = n + 2
      if (miter .ge. 4) lenwm = (2*ml + mu + 1)*n + 2
      lewt = lwm + lenwm
      lsavf = lewt + n
      lacor = lsavf + n
      lenrw = lacor + n - 1
      iwork(17) = lenrw
      liwm = 1
      leniw = 20 + n
      if (miter .eq. 0 .or. miter .eq. 3) leniw = 20
      iwork(18) = leniw

c check rtol and atol for legality. ------------------------------------
      rtoli = rtol(1)
      atoli = atol(1)

c      if(itol.eq.2)then
c      do 70 i = 1,n
c        atoli = atol(i)
c70    continue
c      endif

c-----------------------------------------------------------------------
c block c.
c the next block is for the initial call only (istate = 1).
c it contains all remaining initializations, the initial call to f,
c and the calculation of the initial step size.
c the error weights in ewt are inverted after being loaded.
c-----------------------------------------------------------------------

      uround = r1mach(4)
      tn = t
      jstart = 0
      if (miter .gt. 0) rwork(lwm) = sqrt(uround)
      nhnil = 0
      nst = 0
      nje = 0
      nslast = 0
      hu = 0.0
      nqu = 0
      ccmax = 0.3
      maxcor = 3
      msbp = 20
      mxncf = 10

c initial call to f.  (lf0 points to yh(*,2).) -------------------------
      lf0 = lyh + nyh
      call f (neq, t, y, rwork(lf0))
      nfe = 1

c load the initial value vector in yh. ---------------------------------
      iii = lyh-1
      do 115 i = 1,n
 115    rwork(i+iii) = y(i)

c load and invert the ewt array.  (h is temporarily set to 1.0.) -------
      nq = 1
      h = 1.0
      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))

      iii = lewt-1
      do 120 i = lewt,iii+n
 120    rwork(i) = 1.0/rwork(i)

c-----------------------------------------------------------------------
c the coding below computes the step size, h0, to be attempted on the
c first step, unless the user has supplied a value for this.
c first check that tout - t differs significantly from zero.
c a scalar tolerance quantity tol is computed, as max(rtol(i))
c if this is positive, or max(atol(i)/abs(y(i))) otherwise, adjusted
c so as to be between 100*uround and 1.0e-3.
c then the computed value h0 is given by..
c                                      neq
c   h0**2 = tol / ( w0**-2 + (1/neq) * sum ( f(i)/ywt(i) )**2  )
c                                       1
c where   w0     = max ( abs(t), abs(tout) ),
c         f(i)   = i-th component of initial value of f,
c         ywt(i) = ewt(i)/tol  (a weight for y(i)).
c the sign of h0 is inferred from the initial values of tout and t.
c-----------------------------------------------------------------------

      if (h0 .eq. 0.0e0)then

        tdist = abs(tout - t)
        w0 = max(abs(t),abs(tout))
        tol = rtol(1)

c        if (itol .gt. 2)then
c      do 130 i = 1,n
c 130    tol = max(tol,rtol(i))
c        endif

        if (tol .le. 0.0e0)then
          atoli = atol(1)
      do 150 i = 1,n
c        if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
        ayi = abs(y(i))
        if (ayi .ne. 0.0e0) tol = max(tol,atoli/ayi)
 150    continue
        endif

      tol = max(tol,100.0e0*uround)
      tol = min(tol,0.001e0)
      sum = vnorm (n, rwork(lf0), rwork(lewt))
      sum = 1.0e0/(tol*w0*w0) + tol*sum**2
      h0 = 1.0e0/sqrt(sum)
      h0 = min(h0,tdist)
      h0 = sign(h0,tout-t)

      endif

c adjust h0 if necessary to meet hmax bound. ---------------------------
      rh = abs(h0)*hmxi
      if (rh .gt. 1.0e0) h0 = h0/rh
c load h with h0 and scale yh(*,2) by h0. ------------------------------
      h = h0
      do 190 i = 1,n
 190    rwork(i+lf0-1) = h0*rwork(i+lf0-1)

c-----------------------------------------------------------------------
c block e.
c the next block is normally executed for all calls and contains
c the call to the one-step core integrator stode.
c
c this is a looping point for the integration steps.
c
c first check for too many steps being taken, update ewt (if not at
c start of problem), check for too much accuracy being requested, and
c check for h below the roundoff level in t.
c-----------------------------------------------------------------------
      
      do 270 while ((tn - tout)*h .lt. 0.0)

      tolsf = uround*vnorm (n, rwork(lyh), rwork(lewt))

      if (tolsf .gt. 1.0)then
       tolsf = tolsf*2.0
       go to 580
      endif

      call stodenew (neq, y, rwork(lyh), nyh, 
     &               rwork(lyh), rwork(lewt),
     &   rwork(lsavf), rwork(lacor), rwork(lwm), iwork(liwm),
     &   f, jac, prepj, solsy)

      init = 1

      if ((nst-nslast) .ge. mxstep) go to 580

      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))

      do 260 i = 1,n
        if (rwork(i+lewt-1) .le. 0.0e0) go to 580
 260    rwork(i+lewt-1) = 1.0e0/rwork(i+lewt-1)

270   continue

      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      t = tout

      istate = 2
      illin = 0
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      return
c

c set y vector, t, illin, and optional outputs. ------------------------
 580  do 590 i = 1,n
 590    y(i) = rwork(i+lyh-1)
      t = tn
      illin = 0
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      return

c----------------------- end of subroutine lsode -----------------------
      end


c================================================================
c -- TC2
c
      subroutine stodenew (neq, y, yh, nyh, yh1, ewt, savf, acor,
     &   wm, iwm, f, jac, pjac, slvs)
c     ==========================================================

c---------------------------------------------------------------c
c  STODENEW.F:	A simplified version of STODE.F			c
c			for JSTART >= 0				c
c	--------------------------------------------		c
c								c
c  		Chien Wang					c
c								c
c  Last revised:	March 20, 1995				c
c								c
c---------------------------------------------------------------c

clll. optimize
      external f, jac, pjac, slvs
      integer neq, nyh, iwm
      integer iownd, ialth, ipup, lmax, meo, nqnyh, nslp,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, i1, iredo, iret, j, jb, m, ncf, newq
      real y, yh, yh1, ewt, savf, acor, wm
      real conit, crate, el, elco, hold, rmax, tesco,
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real dcon, ddn, del, delp, dsm, dup, exdn, exsm, exup,
     1   r, rh, rhdn, rhsm, rhup, told, vnorm
      dimension neq(*), y(*), yh(nyh,*), yh1(*), ewt(*), savf(*),
     1   acor(*), wm(*), iwm(*)
      common /ls0001/ conit, crate, el(13), elco(13,12),
     1   hold, rmax, tesco(3,12),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround, iownd(14),
     3   ialth, ipup, lmax, meo, nqnyh, nslp,
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
c-----------------------------------------------------------------------
c stode performs one step of the integration of an initial value
c problem for a system of ordinary differential equations.
c note.. stode is independent of the value of the iteration method
c indicator miter, when this is .ne. 0, and hence is independent
c of the type of chord method used, or the jacobian structure.
c communication with stode is done with the following variables..
c
c neq    = integer array containing problem size in neq(1), and
c          passed as the neq argument in all calls to f and jac.
c y      = an array of length .ge. n used as the y argument in
c          all calls to f and jac.
c yh     = an nyh by lmax array containing the dependent variables
c          and their approximate scaled derivatives, where
c          lmax = maxord + 1.  yh(i,j+1) contains the approximate
c          j-th derivative of y(i), scaled by h**j/factorial(j)
c          (j = 0,1,...,nq).  on entry for the first step, the first
c          two columns of yh must be set from the initial values.
c nyh    = a constant integer .ge. n, the first dimension of yh.
c yh1    = a one-dimensional array occupying the same space as yh.
c ewt    = an array of length n containing multiplicative weights
c          for local error measurements.  local errors in y(i) are
c          compared to 1.0/ewt(i) in various error tests.
c savf   = an array of working storage, of length n.
c          also used for input of yh(*,maxord+2) when jstart = -1
c          and maxord .lt. the current order nq.
c acor   = a work array of length n, used for the accumulated
c          corrections.  on a successful return, acor(i) contains
c          the estimated one-step local error in y(i).
c wm,iwm = real and integer work arrays associated with matrix
c          operations in chord iteration (miter .ne. 0).
c pjac   = name of routine to evaluate and preprocess jacobian matrix
c          and p = i - h*el0*jac, if a chord method is being used.
c slvs   = name of routine to solve linear system in chord iteration.
c ccmax  = maximum relative change in h*el0 before pjac is called.
c h      = the step size to be attempted on the next step.
c          h is altered by the error control algorithm during the
c          problem.  h can be either positive or negative, but its
c          sign must remain constant throughout the problem.
c hmin   = the minimum absolute value of the step size h to be used.
c hmxi   = inverse of the maximum absolute value of h to be used.
c          hmxi = 0.0 is allowed and corresponds to an infinite hmax.
c          hmin and hmxi may be changed at any time, but will not
c          take effect until the next change of h is considered.
c tn     = the independent variable. tn is updated on each step taken.
c jstart = an integer used for input only, with the following
c          values and meanings..
c               0  perform the first step.
c           .gt.0  take a new step continuing from the last.
c              -1  take the next step with a new value of h, maxord,
c                    n, meth, miter, and/or matrix parameters.
c              -2  take the next step with a new value of h,
c                    but with other inputs unchanged.
c          on return, jstart is set to 1 to facilitate continuation.
c kflag  = a completion code with the following meanings..
c               0  the step was succesful.
c              -1  the requested error could not be achieved.
c              -2  corrector convergence could not be achieved.
c              -3  fatal error in pjac or slvs.
c          a return with kflag = -1 or -2 means either
c          abs(h) = hmin or 10 consecutive failures occurred.
c          on a return with kflag negative, the values of tn and
c          the yh array are as of the beginning of the last
c          step, and h is the last step size attempted.
c maxord = the maximum order of integration method to be allowed.
c maxcor = the maximum number of corrector iterations allowed.
c msbp   = maximum number of steps between pjac calls (miter .gt. 0).
c mxncf  = maximum number of convergence failures allowed.
c meth/miter = the method flags.  see description in driver.
c n      = the number of first-order differential equations.
c-----------------------------------------------------------------------
      kflag = 0
      told = tn
      ncf = 0
      ierpj = 0
      iersl = 0
      jcur = 0
      icf = 0
      delp = 0.0e0
c-----------------------------------------------------------------------
c on the first call, the order is set to 1, and other variables are
c initialized.  rmax is the maximum ratio by which h can be increased
c in a single step.  it is initially 1.e4 to compensate for the small
c initial h, but then is normally equal to 10.  if a failure
c occurs (in corrector convergence or error test), rmax is set at 2
c for the next increase.
c
c cfode is called to get all the integration coefficients for the
c current meth.  then the el vector and related constants are reset
c whenever the order nq is changed, or at the start of the problem.
c-----------------------------------------------------------------------
      if(jstart.eq.0)then
        lmax = maxord + 1
        nq = 1
        l = 2
        ialth = 2
        rmax = 10000.0e0
        rc = 0.0e0
        el0 = 1.0e0
        crate = 0.7e0
        hold = h
        meo = meth
        nslp = 0
        ipup = miter
        iret = 3

      call cfode (meth, elco, tesco)

      endif
c-----------------------------------------------------------------------

      if (jstart .gt. 0) go to 200

 150  do 155 i = 1,l
 155    el(i) = elco(i,nq)
      nqnyh = nq*nyh
      rc = rc*el(1)/el0
      el0 = el(1)
      conit = 0.5e0/float(nq+2)

      if(iret.eq.3) go to 200

c-----------------------------------------------------------------------
c if h is being changed, the h ratio rh is checked against
c rmax, hmin, and hmxi, and the yh array rescaled.  ialth is set to
c l = nq + 1 to prevent a change of h for that many steps, unless
c forced by a convergence or error test failure.
c-----------------------------------------------------------------------

 170  rh = max(rh,hmin/abs(h))
      rh = min(rh,rmax)
      rh = rh/max(1.0e0,abs(h)*hmxi*rh)
      r = 1.0e0
      do 180 j = 2,l
        r = r*rh
        do 180 i = 1,n
 180      yh(i,j) = yh(i,j)*r
      h = h*rh
      rc = rc*rh
      ialth = l
      if (iredo .eq. 0) go to 690
c-----------------------------------------------------------------------
c this section computes the predicted values by effectively
c multiplying the yh array by the pascal triangle matrix.
c rc is the ratio of new to old values of the coefficient  h*el(1).
c when rc differs from 1 by more than ccmax, ipup is set to miter
c to force pjac to be called, if a jacobian is involved.
c in any case, pjac is called at least every msbp steps.
c-----------------------------------------------------------------------
 200  if (abs(rc-1.0e0) .gt. ccmax) ipup = miter
      if (nst .ge. nslp+msbp) ipup = miter
      tn = tn + h
      i1 = nqnyh + 1
      do 215 jb = 1,nq
        i1 = i1 - nyh
cdir$ ivdep
        do 210 i = i1,nqnyh
 210      yh1(i) = yh1(i) + yh1(i+nyh)
 215    continue
c-----------------------------------------------------------------------
c up to maxcor corrector iterations are taken.  a convergence test is
c made on the r.m.s. norm of each correction, weighted by the error
c weight vector ewt.  the sum of the corrections is accumulated in the
c vector acor(i).  the yh array is not altered in the corrector loop.
c-----------------------------------------------------------------------
 220  m = 0
      do 230 i = 1,n
 230    y(i) = yh(i,1)
      call f (neq, tn, y, savf)
      nfe = nfe + 1
c-----------------------------------------------------------------------
c if indicated, the matrix p = i - h*el(1)*j is reevaluated and
c preprocessed before starting the corrector iteration.  ipup is set
c to 0 as an indicator that this has been done.
c-----------------------------------------------------------------------
      if (ipup .gt. 0)then

      call pjac (neq, y, yh, nyh, ewt, acor, savf, wm, iwm, f, jac)
      ipup = 0
      rc = 1.0e0
      nslp = nst
      crate = 0.7e0
      if (ierpj .ne. 0) go to 430

      endif

      do 260 i = 1,n
 260    acor(i) = 0.0e0
 270  if (miter .ne. 0) go to 350
c-----------------------------------------------------------------------
c in the case of functional iteration, update y directly from
c the result of the last function evaluation.
c-----------------------------------------------------------------------
      do 290 i = 1,n
        savf(i) = h*savf(i) - yh(i,2)
 290    y(i) = savf(i) - acor(i)
      del = vnorm (n, y, ewt)
      do 300 i = 1,n
        y(i) = yh(i,1) + el(1)*savf(i)
 300    acor(i) = savf(i)
      go to 400
c-----------------------------------------------------------------------
c in the case of the chord method, compute the corrector error,
c and solve the linear system with that as right-hand side and
c p as coefficient matrix.
c-----------------------------------------------------------------------
 350  do 360 i = 1,n
 360    y(i) = h*savf(i) - (yh(i,2) + acor(i))
      call slvs (wm, iwm, y, savf)
      if (iersl .lt. 0) go to 430
      if (iersl .gt. 0) go to 410
      del = vnorm (n, y, ewt)
      do 380 i = 1,n
        acor(i) = acor(i) + y(i)
 380    y(i) = yh(i,1) + el(1)*acor(i)
c-----------------------------------------------------------------------
c test for convergence.  if m.gt.0, an estimate of the convergence
c rate constant is stored in crate, and this is used in the test.
c-----------------------------------------------------------------------
 400  if (m .ne. 0) crate = max(0.2e0*crate,del/delp)
      dcon = del*min(1.0e0,1.5e0*crate)/(tesco(2,nq)*conit)
      if (dcon .le. 1.0e0) go to 450
      m = m + 1
      if (m .eq. maxcor) go to 410
      if (m .ge. 2 .and. del .gt. 2.0e0*delp) go to 410
      delp = del
      call f (neq, tn, y, savf)
      nfe = nfe + 1
      go to 270
c-----------------------------------------------------------------------
c the corrector iteration failed to converge.
c if miter .ne. 0 and the jacobian is out of date, pjac is called for
c the next try.  otherwise the yh array is retracted to its values
c before prediction, and h is reduced, if possible.  if h cannot be
c reduced or mxncf failures have occurred, exit with kflag = -2.
c-----------------------------------------------------------------------
 410  if (miter .eq. 0 .or. jcur .eq. 1) go to 430
      icf = 1
      ipup = miter
      go to 220
 430  icf = 2
      ncf = ncf + 1
      rmax = 2.0e0
      tn = told
      i1 = nqnyh + 1
      do 445 jb = 1,nq
        i1 = i1 - nyh
cdir$ ivdep
        do 440 i = i1,nqnyh
 440      yh1(i) = yh1(i) - yh1(i+nyh)
 445    continue
      if (ierpj .lt. 0 .or. iersl .lt. 0) go to 680
      if (abs(h) .le. hmin*1.00001e0) go to 670
      if (ncf .eq. mxncf) go to 670
      rh = 0.25e0
      ipup = miter
      iredo = 1
      go to 170
c-----------------------------------------------------------------------
c the corrector has converged.  jcur is set to 0
c to signal that the jacobian involved may need updating later.
c the local error test is made and control passes to statement 500
c if it fails.
c-----------------------------------------------------------------------
 450  jcur = 0
      if (m .eq. 0) dsm = del/tesco(2,nq)
      if (m .gt. 0) dsm = vnorm (n, acor, ewt)/tesco(2,nq)
      if (dsm .gt. 1.0e0) go to 500
c-----------------------------------------------------------------------
c after a successful step, update the yh array.
c consider changing h if ialth = 1.  otherwise decrease ialth by 1.
c if ialth is then 1 and nq .lt. maxord, then acor is saved for
c use in a possible order increase on the next step.
c if a change in h is considered, an increase or decrease in order
c by one is considered also.  a change in h is made only if it is by a
c factor of at least 1.1.  if not, ialth is set to 3 to prevent
c testing for that many steps.
c-----------------------------------------------------------------------
      kflag = 0
      iredo = 0
      nst = nst + 1
      hu = h
      nqu = nq
      do 470 j = 1,l
        do 470 i = 1,n
 470      yh(i,j) = yh(i,j) + el(j)*acor(i)
      ialth = ialth - 1
      if (ialth .eq. 0) go to 520
      if (ialth .gt. 1) go to 700
      if (l .eq. lmax) go to 700
      do 490 i = 1,n
 490    yh(i,lmax) = acor(i)
      go to 700
c-----------------------------------------------------------------------
c the error test failed.  kflag keeps track of multiple failures.
c restore tn and the yh array to their previous values, and prepare
c to try the step again.  compute the optimum step size for this or
c one lower order.  after 2 or more failures, h is forced to decrease
c by a factor of 0.2 or less.
c-----------------------------------------------------------------------
 500  kflag = kflag - 1
      tn = told
      i1 = nqnyh + 1
      do 515 jb = 1,nq
        i1 = i1 - nyh
cdir$ ivdep
        do 510 i = i1,nqnyh
 510      yh1(i) = yh1(i) - yh1(i+nyh)
 515    continue
      rmax = 2.0e0
      if (abs(h) .le. hmin*1.00001e0) go to 660
      if (kflag .le. -3) go to 640
      iredo = 2
      rhup = 0.0e0
      go to 540
c-----------------------------------------------------------------------
c regardless of the success or failure of the step, factors
c rhdn, rhsm, and rhup are computed, by which h could be multiplied
c at order nq - 1, order nq, or order nq + 1, respectively.
c in the case of failure, rhup = 0.0 to avoid an order increase.
c the largest of these is determined and the new order chosen
c accordingly.  if the order is to be increased, we compute one
c additional scaled derivative.
c-----------------------------------------------------------------------
 520  rhup = 0.0e0
      if (l .eq. lmax) go to 540
      do 530 i = 1,n
 530    savf(i) = acor(i) - yh(i,lmax)
      dup = vnorm (n, savf, ewt)/tesco(3,nq)
      exup = 1.0e0/float(l+1)
      rhup = 1.0e0/(1.4e0*dup**exup + 0.0000014e0)
 540  exsm = 1.0e0/float(l)
      rhsm = 1.0e0/(1.2e0*dsm**exsm + 0.0000012e0)
      rhdn = 0.0e0
      if (nq .eq. 1) go to 560
      ddn = vnorm (n, yh(1,l), ewt)/tesco(1,nq)
      exdn = 1.0e0/float(nq)
      rhdn = 1.0e0/(1.3e0*ddn**exdn + 0.0000013e0)
 560  if (rhsm .ge. rhup) go to 570
      if (rhup .gt. rhdn) go to 590
      go to 580
 570  if (rhsm .lt. rhdn) go to 580
      newq = nq
      rh = rhsm
      go to 620
 580  newq = nq - 1
      rh = rhdn
      if (kflag .lt. 0 .and. rh .gt. 1.0e0) rh = 1.0e0
      go to 620
 590  newq = l
      rh = rhup
      if (rh .lt. 1.1e0) go to 610
      r = el(l)/float(l)
      do 600 i = 1,n
 600    yh(i,newq+1) = acor(i)*r
      go to 630
 610  ialth = 3
      go to 700
 620  if ((kflag .eq. 0) .and. (rh .lt. 1.1e0)) go to 610
      if (kflag .le. -2) rh = min(rh,0.2e0)
c-----------------------------------------------------------------------
c if there is a change of order, reset nq, l, and the coefficients.
c in any case h is reset according to rh and the yh array is rescaled.
c then exit from 690 if the step was ok, or redo the step otherwise.
c-----------------------------------------------------------------------
      if (newq .eq. nq) go to 170
 630  nq = newq
      l = nq + 1
      iret = 2
      go to 150
c-----------------------------------------------------------------------
c control reaches this section if 3 or more failures have occured.
c if 10 failures have occurred, exit with kflag = -1.
c it is assumed that the derivatives that have accumulated in the
c yh array have errors of the wrong order.  hence the first
c derivative is recomputed, and the order is set to 1.  then
c h is reduced by a factor of 10, and the step is retried,
c until it succeeds or h reaches hmin.
c-----------------------------------------------------------------------
 640  if (kflag .eq. -10) go to 660
      rh = 0.1e0
      rh = max(hmin/abs(h),rh)
      h = h*rh
      do 645 i = 1,n
 645    y(i) = yh(i,1)
      call f (neq, tn, y, savf)
      nfe = nfe + 1
      do 650 i = 1,n
 650    yh(i,2) = h*savf(i)
      ipup = miter
      ialth = 5
      if (nq .eq. 1) go to 200
      nq = 1
      l = 2
      iret = 3
      go to 150
c-----------------------------------------------------------------------
c all returns are made through this section.  h is saved in hold
c to allow the caller to change h on the next step.
c-----------------------------------------------------------------------
 660  kflag = -1
      go to 720
 670  kflag = -2
      go to 720
 680  kflag = -3
      go to 720
 690  rmax = 10.0e0
 700  r = 1.0e0/tesco(2,nqu)
      do 710 i = 1,n
 710    acor(i) = acor(i)*r
 720  hold = h
      jstart = 1
      return
c----------------------- end of subroutine stode -----------------------
      end


c===================================================================
c -- TC3
c
      subroutine ewset (n, itol, rtol, atol, ycur, ewt)
c     =================================================

clll. optimize
c-----------------------------------------------------------------------
c this subroutine sets the error weight vector ewt according to
c     ewt(i) = rtol(i)*abs(ycur(i)) + atol(i),  i = 1,...,n,
c with the subscript on rtol and/or atol possibly replaced by 1 above,
c depending on the value of itol.
c-----------------------------------------------------------------------
c
c  Chien Wang 
c    rewritten 031795
c
      integer n, itol
      integer i
      real rtol, atol, ycur, ewt
      dimension rtol(*), atol(*), ycur(n), ewt(n)
c
      aaa = rtol(1)
      bbb = atol(1)

      do 1 i=1,n
        ewt(i) = aaa*abs(ycur(i)) + bbb
1     continue

      return
       end


c==============================================================
c -- TC4
c
      subroutine solsy (wm, iwm, x, tem)
c     =================================

c
c 	Chien Wang
c	MIT
c	Rewrote 122695
c

clll. optimize
      integer iwm
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, meband, ml, mu
      real wm, x, tem
      real rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real di, hl0, phl0, r
      dimension wm(*), iwm(*), x(*), tem(*)
      common /ls0001/ rowns(209),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     3   iownd(14), iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
c-----------------------------------------------------------------------
c this routine manages the solution of the linear system arising from
c a chord iteration.  it is called if miter .ne. 0.
c if miter is 1 or 2, it calls sgesl to accomplish this.
c if miter = 3 it updates the coefficient h*el0 in the diagonal
c matrix, and then computes the solution.
c if miter is 4 or 5, it calls sgbsl.
c communication with solsy uses the following variables..
c wm    = real work space containing the inverse diagonal matrix if
c         miter = 3 and the lu decomposition of the matrix otherwise.
c         storage of matrix elements starts at wm(3).
c         wm also contains the following matrix-related data..
c         wm(1) = sqrt(uround) (not used here),
c         wm(2) = hl0, the previous value of h*el0, used if miter = 3.
c iwm   = integer work space containing pivot information, starting at
c         iwm(21), if miter is 1, 2, 4, or 5.  iwm also contains band
c         parameters ml = iwm(1) and mu = iwm(2) if miter is 4 or 5.
c x     = the right-hand side vector on input, and the solution vector
c         on output, of length n.
c tem   = vector of work space of length n, not used in this version.
c iersl = output flag (in common).  iersl = 0 if no trouble occurred.
c         iersl = 1 if a singular matrix arose with miter = 3.
c this routine also uses the common variables el0, h, miter, and n.
c-----------------------------------------------------------------------
      iersl = 0

	if(miter.eq.1.or.miter.eq.2)
     &    call sgesl (wm(3), n, n, iwm(21), x, 0)
	if(miter.eq.3)
     &    call solsy2(wm, iwm, x, tem)
	if(miter.eq.4.or.miter.eq.5)then
   	  ml = iwm(1)
          mu = iwm(2)
          meband = 2*ml + mu + 1
          call sgbsl (wm(3), meband, n, 
     &                ml, mu, iwm(21), x, 0)
	endif
	
	return
	 end

c
	subroutine solsy2(wm, iwm, x, tem)
c	==================================
c
c	Chien Wang
c	MIT
c	122695
c
      integer iwm
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, meband, ml, mu
      real wm, x, tem
      real rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real di, hl0, phl0, r
      dimension wm(*), iwm(*), x(*), tem(*)
      common /ls0001/ rowns(209),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     3   iownd(14), iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu

      	  iersl = 0
   	  phl0 = wm(2)
      	  hl0 = h*el0
      	  wm(2) = hl0

          if (hl0 .ne. phl0)then
      	    r = hl0/phl0
            do 320 i = 1,n
        	di = 1.0e0 - r*(1.0e0 - 1.0e0/wm(i+2))
        	if (abs(di) .ne. 0.0e0)then
        	wm(i+2) = 1.0e0/di
	        else
	          iersl = 1
	          goto 360
	   	endif
320	    continue

	  else
    	    do 340 i = 1,n
       	      x(i) = wm(i+2)*x(i)
340	    continue
	  endif

360     return
      	 end

c---- end of subroutine solsy and new solsy2


c==============================================================
c -- TC5
c
      subroutine intdy (t, k, yh, nyh, dky, iflag)
c     ===========================================

clll. optimize
      integer k, nyh, iflag
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, ic, j, jb, jb2, jj, jj1, jp1
      real t, yh, dky
      real rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real c, r, s, tp
      dimension yh(nyh,*), dky(*)
      common /ls0001/ rowns(209),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     3   iownd(14), iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
c-----------------------------------------------------------------------
c intdy computes interpolated values of the k-th derivative of the
c dependent variable vector y, and stores it in dky.  this routine
c is called within the package with k = 0 and t = tout, but may
c also be called by the user for any k up to the current order.
c (see detailed instructions in the usage documentation.)
c-----------------------------------------------------------------------
c the computed values in dky are gotten by interpolation using the
c nordsieck history array yh.  this array corresponds uniquely to a
c vector-valued polynomial of degree nqcur or less, and dky is set
c to the k-th derivative of this polynomial at t.
c the formula for dky is..
c              q
c  dky(i)  =  sum  c(j,k) * (t - tn)**(j-k) * h**(-j) * yh(i,j+1)
c             j=k
c where  c(j,k) = j*(j-1)*...*(j-k+1), q = nqcur, tn = tcur, h = hcur.
c the quantities  nq = nqcur, l = nq+1, n = neq, tn, and h are
c communicated by common.  the above sum is done in reverse order.
c iflag is returned negative if either k or t is out of bounds.
c-----------------------------------------------------------------------
      iflag = 0
      if (k .lt. 0 .or. k .gt. nq) go to 80
      tp = tn - hu -  100.0e0*uround*(tn + hu)
      if ((t-tp)*(t-tn) .gt. 0.0e0) go to 90
c
      s = (t - tn)/h
      ic = 1
      if (k .eq. 0) go to 15
      jj1 = l - k
      do 10 jj = jj1,nq
 10     ic = ic*jj
 15   c = float(ic)
      do 20 i = 1,n
 20     dky(i) = c*yh(i,l)
      if (k .eq. nq) go to 55
      jb2 = nq - k
      do 50 jb = 1,jb2
        j = nq - jb
        jp1 = j + 1
        ic = 1
        if (k .eq. 0) go to 35
        jj1 = jp1 - k
        do 30 jj = jj1,j
 30       ic = ic*jj
 35     c = float(ic)
        do 40 i = 1,n
 40       dky(i) = c*yh(i,jp1) + s*dky(i)
 50     continue
      if (k .eq. 0) return
 55   r = h**(-k)
      do 60 i = 1,n
 60     dky(i) = r*dky(i)
      return
c
 80   call xerrwv(30hintdy--  k (=i1) illegal      ,
     1   30, 51, 0, 1, k, 0, 0, 0.0e0, 0.0e0)
      iflag = -1
      return
 90   call xerrwv(30hintdy--  t (=r1) illegal      ,
     1   30, 52, 0, 0, 0, 0, 1, t, 0.0e0)
      call xerrwv(
     1  60h      t not in interval tcur - hu (= r1) to tcur (=r2)      ,
     1   60, 52, 0, 0, 0, 0, 2, tp, tn)
      iflag = -2
      return
c----------------------- end of subroutine intdy -----------------------
      end


c=============================================================
c -- TC6
c
      subroutine cfode (meth, elco, tesco)
c     ===================================

clll. optimize
      integer meth
      integer i, ib, nq, nqm1, nqp1
      real elco, tesco
      real agamq, fnq, fnqm1, pc, pint, ragq,
     1   rqfac, rq1fac, tsign, xpin
      dimension elco(13,12), tesco(3,12)
c-----------------------------------------------------------------------
c cfode is called by the integrator routine to set coefficients
c needed there.  the coefficients for the current method, as
c given by the value of meth, are set for all orders and saved.
c the maximum order assumed here is 12 if meth = 1 and 5 if meth = 2.
c (a smaller value of the maximum order is also allowed.)
c cfode is called once at the beginning of the problem,
c and is not called again unless and until meth is changed.
c
c the elco array contains the basic method coefficients.
c the coefficients el(i), 1 .le. i .le. nq+1, for the method of
c order nq are stored in elco(i,nq).  they are given by a genetrating
c polynomial, i.e.,
c     l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
c for the implicit adams methods, l(x) is given by
c     dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.
c for the bdf methods, l(x) is given by
c     l(x) = (x+1)*(x+2)* ... *(x+nq)/k,
c where         k = factorial(nq)*(1 + 1/2 + ... + 1/nq).
c
c the tesco array contains test constants used for the
c local error test and the selection of step size and/or order.
c at order nq, tesco(k,nq) is used for the selection of step
c size at order nq - 1 if k = 1, at order nq if k = 2, and at order
c nq + 1 if k = 3.
c-----------------------------------------------------------------------

	if(meth.eq.1) call cfode1(meth, elco, tesco)
	if(meth.eq.2) call cfode2(meth, elco, tesco)

	return
	 end

c
	subroutine cfode1(meth,elco,tesco)
c	==================================

c
c	Chien Wang
c	MIT
c	122695
c

      integer meth
      integer i, ib, nq, nqm1, nqp1
      real elco, tesco
      real agamq, fnq, fnqm1, pc, pint, ragq,
     1   rqfac, rq1fac, tsign, xpin
      dimension elco(13,12), tesco(3,12)
      dimension pc(12)

      elco(1,1) = 1.0e0
      elco(2,1) = 1.0e0
      tesco(1,1) = 0.0e0
      tesco(2,1) = 2.0e0
      tesco(1,2) = 1.0e0
      tesco(3,12) = 0.0e0
      pc(1) = 1.0e0
      rqfac = 1.0e0
      do 140 nq = 2,12
c-----------------------------------------------------------------------
c the pc array will contain the coefficients of the polynomial
c     p(x) = (x+1)*(x+2)*...*(x+nq-1).
c initially, p(x) = 1.
c-----------------------------------------------------------------------
        rq1fac = rqfac
        rqfac = rqfac/float(nq)
        nqm1 = nq - 1
        fnqm1 = float(nqm1)
        nqp1 = nq + 1
c form coefficients of p(x)*(x+nq-1). ----------------------------------
        pc(nq) = 0.0e0
        do 110 ib = 1,nqm1
          i = nqp1 - ib
 110      pc(i) = pc(i-1) + fnqm1*pc(i)
        pc(1) = fnqm1*pc(1)
c compute integral, -1 to 0, of p(x) and x*p(x). -----------------------
        pint = pc(1)
        xpin = pc(1)/2.0e0
        tsign = 1.0e0
        do 120 i = 2,nq
          tsign = -tsign
          pint = pint + tsign*pc(i)/float(i)
 120      xpin = xpin + tsign*pc(i)/float(i+1)
c store coefficients in elco and tesco. --------------------------------
        elco(1,nq) = pint*rq1fac
        elco(2,nq) = 1.0e0
        do 130 i = 2,nq
 130      elco(i+1,nq) = rq1fac*pc(i)/float(i)
        agamq = rqfac*xpin
        ragq = 1.0e0/agamq
        tesco(2,nq) = ragq
        if (nq .lt. 12) tesco(1,nqp1) = ragq*rqfac/float(nqp1)
        tesco(3,nqm1) = ragq
 140    continue
      return
	end

c
	subroutine cfode2(meth, elco, tesco)
c	====================================

c
c	Chien Wang
c	MIT
c	122695
c

      integer meth
      integer i, ib, nq, nqm1, nqp1
      real elco, tesco
      real agamq, fnq, fnqm1, pc, pint, ragq,
     1   rqfac, rq1fac, tsign, xpin
      dimension elco(13,12), tesco(3,12)
      dimension pc(12)
 
      pc(1) = 1.0e0
      rq1fac = 1.0e0
      do 230 nq = 1,5
c-----------------------------------------------------------------------
c the pc array will contain the coefficients of the polynomial
c     p(x) = (x+1)*(x+2)*...*(x+nq).
c initially, p(x) = 1.
c-----------------------------------------------------------------------
        fnq = float(nq)
        nqp1 = nq + 1
c form coefficients of p(x)*(x+nq). ------------------------------------
        pc(nqp1) = 0.0e0
        do 210 ib = 1,nq
          i = nq + 2 - ib
 210      pc(i) = pc(i-1) + fnq*pc(i)
        pc(1) = fnq*pc(1)
c store coefficients in elco and tesco. --------------------------------
        do 220 i = 1,nqp1
 220      elco(i,nq) = pc(i)/pc(2)
        elco(2,nq) = 1.0e0
        tesco(1,nq) = rq1fac
        tesco(2,nq) = float(nqp1)/elco(1,nq)
        tesco(3,nq) = float(nq+2)/elco(1,nq)
        rq1fac = rq1fac/fnq
 230    continue
      return
	end

c---- end of subroutine cfode and new cfode1 & cfode2 

c================================================================
c -- TC7
c
      subroutine sgefa(a,lda,n,ipvt,info)
c     ==================================

clll. optimize
      integer lda,n,ipvt(*),info
      real a(lda,*)
c
c     sgefa factors a real matrix by gaussian elimination.
c
c     sgefa is usually called by sgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for sgeco) = (1 + 9/n)*(time for sgefa) .
c
c     on entry
c
c        a       real(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that sgesl or sgedi will divide by zero
c                     if called.  use  rcond  in sgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 07/14/77 .
c     cleve moler, university of new mexico, argonne national labs.
c
c     subroutines and functions
c
c     blas saxpy,sscal,isamax
c
c     internal variables
c
      real t
      integer isamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = isamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0e0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0e0/a(k,k)
            call sscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call saxpysmp(n-k,t,a(k+1,k),a(k+1,j))
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0e0) info = n
      return
      end

c===============================================================
c -- TC8
c
      subroutine sgbfa(abd,lda,n,ml,mu,ipvt,info)
c     ===========================================

clll. optimize
      integer lda,n,ml,mu,ipvt(*),info
      real abd(lda,*)
c
c     sgbfa factors a real band matrix by elimination.
c
c     sgbfa is usually called by sgbco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c
c     on entry
c
c        abd     real(lda, n)
c                contains the matrix in band storage.  the columns
c                of the matrix are stored in the columns of  abd  and
c                the diagonals of the matrix are stored in rows
c                ml+1 through 2*ml+mu+1 of  abd .
c                see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array  abd .
c                lda must be .ge. 2*ml + mu + 1 .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c                0 .le. ml .lt. n .
c
c        mu      integer
c                number of diagonals above the main diagonal.
c                0 .le. mu .lt. n .
c                more efficient if  ml .le. mu .
c     on return
c
c        abd     an upper triangular matrix in band storage and
c                the multipliers which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that sgbsl will divide by zero if
c                     called.  use  rcond  in sgbco for a reliable
c                     indication of singularity.
c
c     band storage
c
c           if  a  is a band matrix, the following program segment
c           will set up the input.
c
c                   ml = (band width below the diagonal)
c                   mu = (band width above the diagonal)
c                   m = ml + mu + 1
c                   do 20 j = 1, n
c                      i1 = max0(1, j-mu)
c                      i2 = min0(n, j+ml)
c                      do 10 i = i1, i2
c                         k = i - j + m
c                         abd(k,j) = a(i,j)
c                10    continue
c                20 continue
c
c           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
c           in addition, the first  ml  rows in  abd  are used for
c           elements generated during the triangularization.
c           the total number of rows needed in  abd  is  2*ml+mu+1 .
c           the  ml+mu by ml+mu  upper left triangle and the
c           ml by ml  lower right triangle are not referenced.
c
c     linpack. this version dated 07/14/77 .
c     cleve moler, university of new mexico, argonne national labs.
c
c     subroutines and functions
c
c     blas saxpy,sscal,isamax
c     fortran max0,min0
c
c     internal variables
c
      real t
      integer i,isamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1
c
c
      m = ml + mu + 1
      info = 0
c
c     zero initial fill-in columns
c
      j0 = mu + 2
      j1 = min0(n,m) - 1
      if (j1 .lt. j0) go to 30
      do 20 jz = j0, j1
         i0 = m + 1 - jz
         do 10 i = i0, ml
            abd(i,jz) = 0.0e0
   10    continue
   20 continue
   30 continue
      jz = j1
      ju = 0
c
c     gaussian elimination with partial pivoting
c
      nm1 = n - 1
      if (nm1 .lt. 1) go to 140
      do 130 k = 1, nm1
         kp1 = k + 1
c
c        zero next fill-in column
c
         jz = jz + 1
         if (jz .gt. n) go to 60
            if (ml .lt. 1) go to 50
            do 40 i = 1, ml
               abd(i,jz) = 0.0e0
   40       continue
   50       continue
   60    continue
c
c        find l = pivot index
c
         lm = min0(ml,n-k)
         l = isamax(lm+1,abd(m,k),1) + m - 1
         ipvt(k) = l + k - m
c
c        zero pivot implies this column already triangularized
c
         if (abd(l,k) .eq. 0.0e0) go to 110
c
c           interchange if necessary
c
            if (l .eq. m) go to 70
               t = abd(l,k)
               abd(l,k) = abd(m,k)
               abd(m,k) = t
   70       continue
c
c           compute multipliers
c
            t = -1.0e0/abd(m,k)
            call sscal(lm,t,abd(m+1,k),1)
c
c           row elimination with column indexing
c
            ju = min0(max0(ju,mu+ipvt(k)),n)
            mm = m
            if (ju .lt. kp1) go to 100
            do 90 j = kp1, ju
               l = l - 1
               mm = mm - 1
               t = abd(l,j)
               if (l .eq. mm) go to 80
                  abd(l,j) = abd(mm,j)
                  abd(mm,j) = t
   80          continue
               call saxpysmp(lm,t,abd(m+1,k),abd(mm+1,j))
   90       continue
  100       continue
         go to 120
  110    continue
            info = k
  120    continue
  130 continue
  140 continue
      ipvt(n) = n
      if (abd(m,n) .eq. 0.0e0) info = n
      return
      end


c===========================================================
c -- TC9
c
      subroutine sgesl(a,lda,n,ipvt,b,job)
c     ===================================

clll. optimize
      integer lda,n,ipvt(*),job
      real a(lda,*),b(*)
c
c     sgesl solves the real system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by sgeco or sgefa.
c
c     on entry
c
c        a       real(lda, n)
c                the output from sgeco or sgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from sgeco or sgefa.
c
c        b       real(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if sgeco has set rcond .gt. 0.0
c        or sgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call sgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call sgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 07/14/77 .
c     cleve moler, university of new mexico, argonne national labs.
c
c     subroutines and functions
c
c     blas saxpy,sdot
c
c     internal variables
c
      real sdot,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call saxpysmp(n-k,t,a(k+1,k),b(k+1))
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call saxpysmp(k-1,t,a(1,k),b(1))
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = sdot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + sdot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end


c============================================================
c -- TC10
c
      subroutine sgbsl(abd,lda,n,ml,mu,ipvt,b,job)
c     ===========================================

clll. optimize
      integer lda,n,ml,mu,ipvt(*),job
      real abd(lda,*),b(*)
c
c     sgbsl solves the real band system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by sgbco or sgbfa.
c
c     on entry
c
c        abd     real(lda, n)
c                the output from sgbco or sgbfa.
c
c        lda     integer
c                the leading dimension of the array  abd .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c
c        mu      integer
c                number of diagonals above the main diagonal.
c
c        ipvt    integer(n)
c                the pivot vector from sgbco or sgbfa.
c
c        b       real(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b , where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if sgbco has set rcond .gt. 0.0
c        or sgbfa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call sgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call sgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 07/14/77 .
c     cleve moler, university of new mexico, argonne national labs.
c
c     subroutines and functions
c
c     blas saxpy,sdot
c     fortran min0
c
c     internal variables
c
      real sdot,t
      integer k,kb,l,la,lb,lm,m,nm1
c
      m = mu + ml + 1
      nm1 = n - 1
      if (job .ne. 0) go to 60
c
c        job = 0 , solve  a * x = b
c        first solve l*y = b
c
         if (ml .eq. 0) go to 40
            if (nm1 .lt. 1) go to 30
            do 20 k = 1, nm1
               lm = min0(ml,n-k)
               l = ipvt(k)
               t = b(l)
               if (l .eq. k) go to 10
                  b(l) = b(k)
                  b(k) = t
   10          continue
               call saxpysmp(lm,t,abd(m+1,k),b(k+1))
   20       continue
   30       continue
   40    continue
c
c        now solve  u*x = y
c
         do 50 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/abd(m,k)
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = -b(k)
            call saxpysmp(lm,t,abd(la,k),b(lb))
   50    continue
      go to 120
   60 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 70 k = 1, n
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = sdot(lm,abd(la,k),1,b(lb),1)
            b(k) = (b(k) - t)/abd(m,k)
   70    continue
c
c        now solve trans(l)*x = y
c
         if (ml .eq. 0) go to 110
            if (nm1 .lt. 1) go to 100
            do 90 kb = 1, nm1
               k = n - kb
               lm = min0(ml,n-k)
               b(k) = b(k) + sdot(lm,abd(m+1,k),1,b(k+1),1)
               l = ipvt(k)
               if (l .eq. k) go to 80
                  t = b(l)
                  b(l) = b(k)
                  b(k) = t
   80          continue
   90       continue
  100       continue
  110    continue
  120 continue
      return
      end


c========================================================
c -- TC11
c
      subroutine sscal(n,sa,sx,incx)
c     =============================

clll. optimize
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to 1.
c     jack dongarra, linpack, 6/17/77.
c
      real sa,sx(*)
      integer i,incx,m,mp1,n,nincx
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        sx(i) = sa*sx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sx(i) = sa*sx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        sx(i) = sa*sx(i)
        sx(i + 1) = sa*sx(i + 1)
        sx(i + 2) = sa*sx(i + 2)
        sx(i + 3) = sa*sx(i + 3)
        sx(i + 4) = sa*sx(i + 4)
   50 continue
      return
      end


c=====================================================
c -- TC12
c
      subroutine saxpysmp(n,sa,sx,sy)
c     =====================================

c=====================================================
c	saxpysmp:  a simplified version of saxpy with
c		   incx = incy = 1
c	---------------------------------------------
c	Chien Wang
c	MIT
c
c	Creates:	010795
c=====================================================

clll. optimize
c
c     constant times a vector plus a vector.
c     uses unrolled loop for increments equal to one.
c     jack dongarra, linpack, 6/17/77.
c
      real sx(*),sy(*),sa
      integer i,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (sa .eq. 0.0) return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
      m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sy(i) + sa*sx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        sy(i) = sy(i) + sa*sx(i)
        sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
        sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
        sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
   50 continue
      return
      end


c=======================================================================
c -- TC13
c
      subroutine xerrwv (msg, nmes, nerr, level, ni, i1, i2, nr, r1, r2)
c     ==================================================================

      integer msg, nmes, nerr, level, ni, i1, i2, nr,
     1   i, lun, lunit, mesflg, ncpw, nch, nwds
      real r1, r2
      dimension msg(nmes)
c-----------------------------------------------------------------------
c subroutines xerrwv, xsetf, and xsetun, as given here, constitute
c a simplified version of the slatec error handling package.
c written by a. c. hindmarsh at llnl.  version of march 30, 1987.
c
c all arguments are input arguments.
c
c msg    = the message (hollerith literal or integer array).
c nmes   = the length of msg (number of characters).
c nerr   = the error number (not used).
c level  = the error level..
c          0 or 1 means recoverable (control returns to caller).
c          2 means fatal (run is aborted--see note below).
c ni     = number of integers (0, 1, or 2) to be printed with message.
c i1,i2  = integers to be printed, depending on ni.
c nr     = number of reals (0, 1, or 2) to be printed with message.
c r1,r2  = reals to be printed, depending on nr.
c
c note..  this routine is machine-dependent and specialized for use
c in limited context, in the following ways..
c 1. the number of hollerith characters stored per word, denoted
c    by ncpw below, is a data-loaded constant.
c 2. the value of nmes is assumed to be at most 60.
c    (multi-line messages are generated by repeated calls.)
c 3. if level = 2, control passes to the statement   stop
c    to abort the run.  this statement may be machine-dependent.
c 4. r1 and r2 are assumed to be in single precision and are printed
c    in e21.13 format.
c 5. the common block /eh0001/ below is data-loaded (a machine-
c    dependent feature) with default values.
c    this block is needed for proper retention of parameters used by
c    this routine which the user can reset by calling xsetf or xsetun.
c    the variables in this block are as follows..
c       mesflg = print control flag..
c                1 means print all messages (the default).
c                0 means no printing.
c       lunit  = logical unit number for messages.
c                the default is 6 (machine-dependent).
c-----------------------------------------------------------------------
c the following are instructions for installing this routine
c in different machine environments.
c
c to change the default output unit, change the data statement below.
c
c for some systems, the data statement below must be replaced
c by a separate block data subprogram.
c
c for a different number of characters per word, change the
c data statement setting ncpw below, and format 10.  alternatives for
c various computers are shown in comment cards.
c
c for a different run-abort command, change the statement following
c statement 100 at the end.
c-----------------------------------------------------------------------
      common /eh0001/ mesflg, lunit
c
      data mesflg/1/, lunit/6/
c-----------------------------------------------------------------------
c the following data-loaded value of ncpw is valid for the cdc-6600
c and cdc-7600 computers.
c     data ncpw/10/
c the following is valid for the cray-1 computer.
      data ncpw/8/
c the following is valid for the burroughs 6700 and 7800 computers.
c     data ncpw/6/
c the following is valid for the pdp-10 computer.
c     data ncpw/5/
c the following is valid for the vax computer with 4 bytes per integer,
c and for the ibm-360, ibm-370, ibm-303x, and ibm-43xx computers.
c     data ncpw/4/
c the following is valid for the pdp-11, or vax with 2-byte integers.
c     data ncpw/2/
c-----------------------------------------------------------------------
c
      if (mesflg .eq. 0) go to 100
c get logical unit number. ---------------------------------------------
      lun = lunit
c get number of words in message. --------------------------------------
      nch = min0(nmes,60)
      nwds = nch/ncpw
      if (nch .ne. nwds*ncpw) nwds = nwds + 1
c write the message. ---------------------------------------------------
      write (lun, 10) (msg(i),i=1,nwds)
c-----------------------------------------------------------------------
c the following format statement is to have the form
c 10  format(1x,mmann)
c where nn = ncpw and mm is the smallest integer .ge. 60/ncpw.
c the following is valid for ncpw = 10.
c 10  format(1x,6a10)
c the following is valid for ncpw = 8.
  10  format(1x,8a8)
c the following is valid for ncpw = 6.
c 10  format(1x,10a6)
c the following is valid for ncpw = 5.
c 10  format(1x,12a5)
c the following is valid for ncpw = 4.
c 10  format(1x,15a4)
c the following is valid for ncpw = 2.
c 10  format(1x,30a2)
c-----------------------------------------------------------------------
      if (ni .eq. 1) write (lun, 20) i1
 20   format(6x,"in above message,  i1 =",i10)
      if (ni .eq. 2) write (lun, 30) i1,i2
 30   format(6x,"in above message,  i1 =",i10,3x,"i2 =",i10)
      if (nr .eq. 1) write (lun, 40) r1
 40   format(6x,"in above message,  r1 =",e21.13)
      if (nr .eq. 2) write (lun, 50) r1,r2
 50   format(6x,"in above,  r1 =",e21.13,3x,"r2 =",e21.13)
c abort the run if level = 2. ------------------------------------------
 100  if (level .ne. 2) return
      stop
c----------------------- end of subroutine xerrwv ----------------------
      end


c=============================================================
c -- TC14
c
      real function r1mach (idum)
c     ===========================

      integer idum
c-----------------------------------------------------------------------
c this routine computes the unit roundoff of the machine.
c this is defined as the smallest positive machine number
c u such that  1.0 + u .ne. 1.0
c-----------------------------------------------------------------------
      real u, comp
      u = 1.0e0
 10   u = u*0.5e0
      comp = 1.0e0 + u
      if (comp .ne. 1.0e0) go to 10
      r1mach = u*2.0e0
      return
c----------------------- end of function r1mach ------------------------
      end


c========================================================
c -- TC15
c
      real function vnorm (n, v, w)
c     ============================

clll. optimize
c-----------------------------------------------------------------------
c this function routine computes the weighted root-mean-square norm
c of the vector of length n contained in the array v, with weights
c contained in the array w of length n..
c   vnorm = sqrt( (1/n) * sum( v(i)*w(i) )**2 )
c-----------------------------------------------------------------------
      integer n,   i
      real v, w,   sum
      dimension v(n), w(n)
      sum = 0.0e0
      do 10 i = 1,n
 10     sum = sum + (v(i)*w(i))**2
      vnorm = sqrt(sum/float(n))
      return
c----------------------- end of function vnorm -------------------------
      end


c=========================================================
c -- TC16
c
      integer function isamax(n,sx,incx)
c     =================================

clll. optimize
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 6/17/77.
c
      real sx(*),smax
      integer i,incx,ix,n
c
      isamax = 1
      if(n.le.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smax = abs(sx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(abs(sx(ix)).le.smax) go to 5
         isamax = i
         smax = abs(sx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smax = abs(sx(1))
      do 30 i = 2,n
         if(abs(sx(i)).le.smax) go to 30
         isamax = i
         smax = abs(sx(i))
   30 continue
      return
      end


c=======================================================
c -- TC17
c
      real function sdot(n,sx,incx,sy,incy)
c     =====================================

clll. optimize
c
c     forms the dot product of a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 6/17/77.
c
      real sx(*),sy(*),stemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      stemp = 0.0e0
      sdot = 0.0e0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp = stemp + sx(ix)*sy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      sdot = stemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = stemp + sx(i)*sy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) +
     *   sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
   50 continue
   60 sdot = stemp
      return
      end


c================================================================
c -- TC18
c
      subroutine prepj (neq, y, yh, nyh, ewt, ftem, savf, wm, iwm,
     1   f, jac)
c     ===========================================================

clll. optimize
      external f, jac
      integer neq, nyh, iwm
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, i1, i2, ier, ii, j, j1, jj, lenp,
     1   mba, mband, meb1, meband, ml, ml3, mu, np1
      real y, yh, ewt, ftem, savf, wm
      real rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real con, di, fac, hl0, r, r0, srur, yi, yj, yjj,
     1   vnorm
      dimension neq(*), y(*), yh(nyh,*), ewt(*), ftem(*), savf(*),
     1   wm(*), iwm(*)
      common /ls0001/ rowns(209),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     3   iownd(14), iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
c-----------------------------------------------------------------------
c prepj is called by stode to compute and process the matrix
c p = i - h*el(1)*j , where j is an approximation to the jacobian.
c here j is computed by the user-supplied routine jac if
c miter = 1 or 4, or by finite differencing if miter = 2, 3, or 5.
c if miter = 3, a diagonal approximation to j is used.
c j is stored in wm and replaced by p.  if miter .ne. 3, p is then
c subjected to lu decomposition in preparation for later solution
c of linear systems with p as coefficient matrix. this is done
c by sgefa if miter = 1 or 2, and by sgbfa if miter = 4 or 5.
c
c in addition to variables described previously, communication
c with prepj uses the following..
c y     = array containing predicted values on entry.
c ftem  = work array of length n (acor in stode).
c savf  = array containing f evaluated at predicted y.
c wm    = real work space for matrices.  on output it contains the
c         inverse diagonal matrix if miter = 3 and the lu decomposition
c         of p if miter is 1, 2 , 4, or 5.
c         storage of matrix elements starts at wm(3).
c         wm also contains the following matrix-related data..
c         wm(1) = sqrt(uround), used in numerical jacobian increments.
c         wm(2) = h*el0, saved for later use if miter = 3.
c iwm   = integer work space containing pivot information, starting at
c         iwm(21), if miter is 1, 2, 4, or 5.  iwm also contains band
c         parameters ml = iwm(1) and mu = iwm(2) if miter is 4 or 5.
c el0   = el(1) (input).
c ierpj = output error flag,  = 0 if no trouble, .gt. 0 if
c         p matrix found to be singular.
c jcur  = output flag = 1 to indicate that the jacobian matrix
c         (or approximation) is now current.
c this routine also uses the common variables el0, h, tn, uround,
c miter, n, nfe, and nje.
c-----------------------------------------------------------------------
      nje = nje + 1
      ierpj = 0
      jcur = 1
      hl0 = h*el0
      go to (100, 200, 300, 400, 500), miter
c if miter = 1, call jac and multiply by scalar. -----------------------
 100  lenp = n*n
      do 110 i = 1,lenp
 110    wm(i+2) = 0.0e0
      call jac (neq, tn, y, 0, 0, wm(3), n)
      con = -hl0
      do 120 i = 1,lenp
 120    wm(i+2) = wm(i+2)*con
      go to 240
c if miter = 2, make n calls to f to approximate j. --------------------
 200  fac = vnorm (n, savf, ewt)
      r0 = 1000.0e0*abs(h)*uround*float(n)*fac
      if (r0 .eq. 0.0e0) r0 = 1.0e0
      srur = wm(1)
      j1 = 2
      do 230 j = 1,n
        yj = y(j)
        r = max(srur*abs(yj),r0/ewt(j))
        y(j) = y(j) + r
        fac = -hl0/r
        call f (neq, tn, y, ftem)
        do 220 i = 1,n
 220      wm(i+j1) = (ftem(i) - savf(i))*fac
        y(j) = yj
        j1 = j1 + n
 230    continue
      nfe = nfe + n
c add identity matrix. -------------------------------------------------
 240  j = 3
      np1 = n + 1
      do 250 i = 1,n
        wm(j) = wm(j) + 1.0e0
 250    j = j + np1
c do lu decomposition on p. --------------------------------------------
      call sgefa (wm(3), n, n, iwm(21), ier)
      if (ier .ne. 0) ierpj = 1
      return
c if miter = 3, construct a diagonal approximation to j and p. ---------
 300  wm(2) = hl0
      r = el0*0.1e0
      do 310 i = 1,n
 310    y(i) = y(i) + r*(h*savf(i) - yh(i,2))
      call f (neq, tn, y, wm(3))
      nfe = nfe + 1
      do 320 i = 1,n
        r0 = h*savf(i) - yh(i,2)
        di = 0.1e0*r0 - h*(wm(i+2) - savf(i))
        wm(i+2) = 1.0e0
        if (abs(r0) .lt. uround/ewt(i)) go to 320
        if (abs(di) .eq. 0.0e0) go to 330
        wm(i+2) = 0.1e0*r0/di
 320    continue
      return
 330  ierpj = 1
      return
c if miter = 4, call jac and multiply by scalar. -----------------------
 400  ml = iwm(1)
      mu = iwm(2)
      ml3 = ml + 3
      mband = ml + mu + 1
      meband = mband + ml
      lenp = meband*n
      do 410 i = 1,lenp
 410    wm(i+2) = 0.0e0
      call jac (neq, tn, y, ml, mu, wm(ml3), meband)
      con = -hl0
      do 420 i = 1,lenp
 420    wm(i+2) = wm(i+2)*con
      go to 570
c if miter = 5, make mband calls to f to approximate j. ----------------
 500  ml = iwm(1)
      mu = iwm(2)
      mband = ml + mu + 1
      mba = min0(mband,n)
      meband = mband + ml
      meb1 = meband - 1
      srur = wm(1)
      fac = vnorm (n, savf, ewt)
      r0 = 1000.0e0*abs(h)*uround*float(n)*fac
      if (r0 .eq. 0.0e0) r0 = 1.0e0
      do 560 j = 1,mba
        do 530 i = j,n,mband
          yi = y(i)
          r = max(srur*abs(yi),r0/ewt(i))
 530      y(i) = y(i) + r
        call f (neq, tn, y, ftem)
        do 550 jj = j,n,mband
          y(jj) = yh(jj,1)
          yjj = y(jj)
          r = max(srur*abs(yjj),r0/ewt(jj))
          fac = -hl0/r
          i1 = max0(jj-mu,1)
          i2 = min0(jj+ml,n)
          ii = jj*meb1 - ml + 2
          do 540 i = i1,i2
 540        wm(ii+i) = (ftem(i) - savf(i))*fac
 550      continue
 560    continue
      nfe = nfe + mba
c add identity matrix. -------------------------------------------------
 570  ii = mband + 2
      do 580 i = 1,n
        wm(ii) = wm(ii) + 1.0e0
 580    ii = ii + meband
c do lu decomposition of p. --------------------------------------------
      call sgbfa (wm(3), meband, n, ml, mu, iwm(21), ier)
      if (ier .ne. 0) ierpj = 1
      return
c----------------------- end of subroutine prepj -----------------------
      end
