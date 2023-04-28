 SUBROUTINE TESTEPPA(SO2EM,grids1,filein1)

 
  USE EPPANEW_MOD
  parameter ( nchemyr=1,IYEARL=2100,IYEARF=2006)
   REAL,  INTENT(OUT) :: SO2EM
   CHARACTER(LEN=255), INTENT(IN) :: GRIDS1
   CHARACTER(LEN=255), INTENT(IN) :: filein1
  integer n_total_urbana (nlon,nlat,nchemyr)
  integer n_urbana (3,nlon,nlat,nchemyr)
      common /EPPAEMISS/edailyf11     (nlon,nlat,nchemyr), &
       edailyf12     (nlon,nlat,nchemyr), &
       edailyn2o     (nlon,nlat,nchemyr),&
       edailyco      (nlon,nlat,nchemyr),&
       edailynox     (nlon,nlat,nchemyr),&
       edailych4     (nlon,nlat,nchemyr),&
       edailyso2     (nlon,nlat,nchemyr),&
       edailyco2     (nlon,nlat,nchemyr),&
       edailyhfc134a (nlon,nlat,nchemyr),&
       edailypfc     (nlon,nlat,nchemyr),&
       edailysf6     (nlon,nlat,nchemyr),&
       edailyuco     (nlon,nlat,nchemyr),&
       edailyunmv    (nlon,nlat,nchemyr),&
       edailyunox    (nlon,nlat,nchemyr),&
       edailyusox    (nlon,nlat,nchemyr),&
       edailybc      (nlon,nlat,nchemyr),&
       edailynh3     (nlon,nlat,nchemyr),&
       edailyoc      (nlon,nlat,nchemyr),&
       edailyubc     (nlon,nlat,nchemyr),&
       edailyunh3    (nlon,nlat,nchemyr),&
       edailyuoc     (nlon,nlat,nchemyr),&
       ehpbl         (nlon,nlat)
    common /chem_meta_para2/n_total_urbana,&
       n_urbana
  character *80 names
  print *,nlat
  print *,nlon
  LYEAR1=(IYEARL-2000)/5+2
  if(LYEAR1.ne.LYEAR) then
   print *,'LYEAR=',LYEAR
   print *,'YEARL=',YEARL
   print *,'LYEAR1=',LYEAR1
   stop
  endif

  open (100, file='edaily.dat',form='unformatted')
! open (101, file='edaily.so2',form='unformatted')
  open (101, file='edaily.so2',form='formatted')
      write(100)IYEARF,IYEARL
 grids=grids1
 print *,grids
 filein=filein1
 print *,filein
! nem=15
  nem=1
  cfemi=1./365.*1.e15
  do iyear=IYEARF,IYEARL
   print *,iyear,nem
   print *,'CALL EPPA_DRIVER'
   
   CALL EPPA_DRIVER(iyear)
!  print *,' SO21997=',SO21997
!  print *,'EMISTEMPCO2'
!  print *,EMISTEMPCO2
    edailyf11(:,:,nem)=EMISTEMPC11*cfemi
    edailyf12(:,:,nem)=EMISTEMPc12*cfemi
    edailyn2o(:,:,nem)=EMISTEMPN2o*cfemi
    edailyco(:,:,nem)=EMISTEMPCO*cfemi
    edailynox(:,:,nem)=EMISTEMPNO2*cfemi
    edailych4(:,:,nem)=EMISTEMPCH4*cfemi
    edailyso2(:,:,nem)=EMISTEMPSO2*cfemi
    edailyco2(:,:,nem)=EMISTEMPCO2*cfemi
    edailyhfc134a(:,:,nem)=EMISTEMPhfc*cfemi
    edailypfc(:,:,nem)=EMISTEMPpfc*cfemi
    edailysf6(:,:,nem)=EMISTEMPsf6*cfemi
    edailyuco(:,:,nem)=EMISTEMPuco*cfemi
    edailyunmv(:,:,nem)=EMISTEMPunmv*cfemi
    edailyunox(:,:,nem)=EMISTEMPuno2*cfemi
    edailyusox(:,:,nem)=EMISTEMPuso2*cfemi
    edailybc(:,:,nem)=EMISTEMPBC*cfemi
    edailynh3(:,:,nem)=EMISTEMPNH3*cfemi
    edailyoc(:,:,nem)=EMISTEMPOC*cfemi
    edailyubc(:,:,nem)=EMISTEMPubc*cfemi
    edailyunh3(:,:,nem)=EMISTEMPunh3*cfemi
    edailyuoc(:,:,nem)=EMISTEMPuOC*cfemi
    n_total_urbana(:,:,nem)=N_TOTAL_URBAN
    n_urbana(:,:,:,nem)=N_URBAN
      write(100)iyear
      write(100)edailyf11, edailyf12, edailyn2o, edailyco, edailynox, &
       edailych4,edailyso2,edailyco2,edailyhfc134a,edailypfc,edailysf6, &
       edailyuco,edailyunmv, edailyunox,edailyusox,edailybc,edailynh3, &
       edailyoc,edailyubc,edailyunh3,edailyuoc,n_total_urbana,n_urbana
  enddo
  close(100)
   rewind(17)
   read (17,*),names
   print *, names
  do i=IYEARF,IYEARL
   read (17,*),iyear,co2,f11,f12,n2o,co,no,ch4,so2,hfc,pfc,sf6,bc,nh3,oc
   print *,i,iyear
   if (iyear.eq.2006) then
    SO21997=so2
    print *,iyear,SO21997
    exit
   endif
  enddo
 SO2EM=SO21997
!      write (101),SO21997
       write (101,*),SO21997
  close(101)
 
 RETURN
 END SUBROUTINE TESTEPPA
