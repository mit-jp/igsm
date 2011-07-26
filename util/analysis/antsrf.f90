program main
  integer YEARF(10),YEARL(10)
  dimension nytf(10)
  namelist /PAR/JM,IM,LM,IM0,YEARF,YEARL,nmf
  open( 19,file='antsrf.name',status='old')
  read(19,PAR)

  LMP1=LM+1
  JMP3=JM+3
  JMP1=JM+1
  nytt=0

  do nn=1,nmf
    print *,YEARL(nn),YEARF(nn)
    nyt=YEARL(nn)-YEARF(nn)
    print *,' nyt=',nyt
    nytf(nn)=nyt
    nytt=nytt+nyt
  enddo
  print *,' nytt=',nytt
  call srftemp(JM,IM,LM,IM0,LMP1,JMP1,JMP3,YEARF,YEARL,nmf,nytt,nytf)
end

subroutine srftemp(JM,IM,LM,IM0,LMP1,JMP1,JMP3,YEARFA,YEARLA,nmf,nyt,nytf)
  PARAMETER (NQTAB=33)
  dimension GBUDG(JMP3,89,4), QMAPS(JMP3,57), QTABLE(JMP3,LMP1,NQTAB), &
    INQTAB(NQTAB), J1QT(NQTAB), INQMAP(57)
  character *256 flname(25),flaver,flout,floutascii,flsen
!  INCLUDE 'COMP.COM'
  character *4 JMONTH,JMNTH0,AMONTH(12),MONTHF,MONTHL
  character *4 YEAR,YEAR0,XLABLE(33)
  character * 3 FOR
  character * 16 SK
!  character * 11 SENS
  equivalence (SK,XLABLE(21))
  character * 9 DATAFR
  character *10 var
  integer YEAR1,YEAR2
  integer YEARF,YEARL
  integer YEARFA(10),YEARLA(10)
  dimension nytf(10)
  dimension fland(JM),fice(JM),fwater(JM),rice(JM),focean(JM),FDATA(IM0,JM,3), &
    ficeav(JM),fwatav(JM)
  dimension yp(JM),cosp(JM),cosv(JM),yv(JM+1), dxyp(JM), fricea(5,nyt), &
    oicean(nyt),ttg(nyt)
  namelist /ave/ datafr, flname, flaver, var, ivr, opf, floutascii
  logical OPF
!  namelist/AVE/DATAFR,flaver,flname,flout,FOR,floutascii,flsen
  data AMONTH/'JAN','FEB','MAR','APR','MAY','JUNE','JULY','AUG','SEP','OCT',&
    'NOV','DEC'/
  dimension ndmonth(12),ndseason(4),ndaver(12)
  DATA ndmonth/31,28,31,30,31,30,31,31,30,31,30,31/
!  DJF MAM JJA SON
  DATA ndseason/90,92,92,91/
  nave=0
  OPF=.false.
  read(19,AVE)
  if(var.eq.'GLOBAL')ntype=1
  if(var.eq.'LAND  ')ntype=2
  if(var.eq.'OCEAN ')ntype=3
  print *,'ivr=',ivr
  print *,'var=',var
  pi=4.*atan(1.)
  iyp=JM
  dyp=pi/float(iyp-1)
  yp(1)=-pi/2.
  cosp(1)=0.
  do j=2,iyp
    yp(j)=yp(j-1)+dyp
    cosp(j)=cos(yp(j))
    cosv(j)=0.5*(cosp(j)+cosp(j-1))
  enddo
  cosv(iyp/2+1)=1.
  dxyp(1)=0.5*cosv(2)*0.5
  dxyp(JM)=0.5*cosv(JM)*0.5
  do j=2,iyp-1
    dxyp(j)=0.5*(cosv(j)+cosv(j+1))
  enddo
  ss=0.0
  do j=1,jm
    ss=ss+dxyp(j)
  enddo
  print *,' nyt=',nyt
!  print *,flsen
!  open (66,file='flsen',status='old',form='formatted')
!  open (66,file=flsen,status='old',form='formatted')
  n0=0
  do nn=1,nmf
    YEARF=YEARFA(nn)
    YEARL=YEARLA(nn)
    print *,YEARF,YEARL 
    nyt1=YEARL-YEARF
    if(nn.gt.1)n0=n0+nytf(nn-1)
    print *,nn,n0,nyt1
    open( unit=45,file=flname(nn), status='old', form='unformatted')
    if (opf) then
      open(unit=55, file=flaver, status='new', form='unformatted')
      open(unit=441, file=floutascii, status='new', form='formatted')
    endif
    print *,flname(nn)
    read(45) ANEXP,TAU,XLABLE
    print *,' ANEXP=',ANEXP
!    print *,XLABLE
    if(nn.eq.1)ANEXP0=ANEXP
    if(ANEXP.ne.ANEXP0)then
      print *,'Disagreement in run numbers'
      print *,'nn=',nn,' ANEXP=',ANEXP,' ANEXP0=',ANEXP0
      if(ANEXP0.eq.115.03.and.ANEXP.eq.2500.03)then
        print *,'continue'
      else
!        stop
      endif
    endif
100 continue
    READ (45,end=800) ANEXP,JDATE,JMONTH,JYEAR,JDATE0,JMNTH0,JYEAR0
    do i=1,12
      if(JMONTH.eq.AMONTH(i))JM1=i
      if(JMNTH0.eq.AMONTH(i))JM0=i
    enddo
    if(DATAFR.eq.'MONTHLY')then
!     nrpy = number of records per year
      nrpy=12
      do ii=1,12
        ndaver(ii)=ndmonth(ii)
      enddo
      if(JYEAR.lt.YEARF) goto 100
      if(JYEAR0.eq.YEARF) backspace(45)
    elseif(DATAFR.eq.'SEASONALY')then
!      nrpy = number of records per year
      nrpy=4
      do ii=1,4
        ndaver(ii)=ndseason(ii)
      enddo
!      print *,' Check condition'
!      print *,JYEAR,YEARF,JM0
      if(JYEAR.lt.YEARF) goto 100
      if(JYEAR.eq.YEARF.and.JM0.eq.12)backspace(45)
    elseif(DATAFR.eq.'ANNUALLY')then
!      nrpy = number of records per year
      nrpy=1
      ndaver(1)=365
!      print *,' Check condition'
!      print *,JYEAR,YEARF,JM0
      if(JYEAR-1.lt.YEARF)go to 100
      if(JYEAR-1.eq.YEARF)backspace(45)
    else
      print *,' DATAFR is wrong'
      print *,' DATAFR =',DATAFR
      stop
    endif
!    if(JM0.eq.JMF.and.JYEAR0.eq.YEARF)backspace(45)
    fnrpy=float(nrpy)
    fnrpy=nrpy
!    print *,' number of records per year',fnrpy
    do ny=1,nyt1
      ttg(n0+ny)=0.0
      si=0.0
      do mn=1,nrpy
        READ (45,end=800) ANEXP,JDATE,JMONTH,JYEAR,JDATE0,JMNTH0,JYEAR0,GBUDG,&
          QMAPS,QTABLE,INQTAB,J1QT,INQMAP
!          print *,'read for averaging begin'
!          print *,JYEAR0,JYEAR
!          print *,JMNTH0,JMONTH
!          print *,'read for averaging end'
        do i=1,12
          if(JMONTH.eq.AMONTH(i))JM1=i
          if(JMNTH0.eq.AMONTH(i))JM0=i
        enddo
        if(ny.eq.1.and.mn.eq.1)then
          print *,' FROM'
          print *,'    ',JMNTH0,JDATE0,JYEAR0,JM0
!          YEAR1=JYEAR0
          if(n0.eq.0) YEAR1=JYEAR0
        end if
        if(ny.eq.nyt1.and.mn.eq.nrpy)then
          print *,'  TO'
          print *,'    ',JMONTH,JDATE,JYEAR,JM1
          YEAR2=JYEAR
        end if
!        ttg(n0+ny)=ttg(n0+ny)+GBUDG(JMP3,ivr,ntype)/fnrpy
        ttg(n0+ny)=ttg(n0+ny)+GBUDG(JMP3,ivr,ntype)*ndaver(mn)/365.
        do j=1,JM
          fri=0.01*GBUDG(j,64,1)
          si=si+fri*dxyp(j)*ndaver(mn)
        enddo
      enddo
!    print *,' ny=',ny
      oicean(n0+ny)=si/ss/365.0
    enddo
    close(45)
  enddo !nn
  YEAR1=YEARFA(1)
  dty=1.
!  print *,'RUN'
!  print *,ANEXP,YEAR1,YEAR2,nyt,dty
  print *,' tgan'
  if (nyt.gt.10)then
    print *,(0.1*ttg(n),n=1,nyt,nyt/10)
  else
    print *,(0.1*ttg(n),n=1,nyt)
  endif
!  print *,' seaice'
!  print *,(oicean(n),n=1,nyt)
!  do n=1,nyt
!    print *,YEAR1+n-1,0.1*ttg(n)
!  enddo
!  YEAR1=1
  YEAR2=YEAR1+nyt
!  ANEXP=2910.03
!  ANEXP=3811.03
  if(.not.OPF) then
    open(unit=55, file=flaver, status='new', form='unformatted')
    open(unit=441, file=floutascii, status='new', form='formatted')
  endif
  print *,ANEXP,YEAR1,YEAR2,nyt,dty
  write (55) ANEXP, YEAR1, YEAR2, nyt, dty
  write (55),ttg, oicean
  do ii=1,nyt
    write (441,'i5,f10.3'),YEAR1+ii-1,0.1*ttg(ii)
  enddo

  goto 900
800 continue
  print *,' end of file'
  print *,'    ',JMNTH0,JDATE0,JYEAR0,JM0
900 continue
  close(55)
  close(441)
120 format(f6.2)
109 format (a3)
110 format (i4)
901 format(f6.2)
end
