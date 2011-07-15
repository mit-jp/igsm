#include "ctrparam.h
      PROGRAM TESTEPPA

 
       USE EPPANEW_MOD

#include "chem_para"
#include "chem_com"
#include "chem_meta"

       print *,nlat
       print *,nlon
!      open (10, file='edaily.new.data',form='unformatted')
       open (100, file='edaily.dat',form='unformatted')
       open (101, file='edaily.so2',form='formatted')
       nem=15
       cfemi=1./365.*1.e15
! do iyear=1991,2100
       do iyear=1997,2010
        print *,iyear,nem
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
         n_total_urban(:,nem)=N_TOTAL_URBAN
         n_urban(:,:,nem)=N_URBAN
       nem=nem+1
       enddo
!  print *,'edailyco2'
!  print *,(edailyco2(1,i,nem),i=1,46)
!  write (100),EMISTEMPC11,EMISTEMPc12,EMISTEMPN2o,EMISTEMPCO,EMISTEMPNO2, &
!          EMISTEMPCH4,EMISTEMPSO2,EMISTEMPCO2,EMISTEMPhfc,EMISTEMPpfc,EMISTEMPsf6, &
!          EMISTEMPuco,EMISTEMPunmv,EMISTEMPuno2,EMISTEMPuso2, &
!          EMISTEMPBC,EMISTEMPNH3,EMISTEMPOC, &
!          EMISTEMPubc,EMISTEMPunh3,EMISTEMPuOC,n_total_urban,n_urban
      write(100)edailyf11, edailyf12, edailyn2o, edailyco, edailynox, 
     $  edailych4,edailyso2,edailyco2,edailyhfc134a,edailypfc,edailysf6, 
     $  edailyuco,edailyunmv, edailyunox,edailyusox,edailybc,edailynh3, 
     $  edailyoc,edailyubc,edailyunh3,edailyuoc,n_total_urbana,n_urbana
       close(100)
       write (101,*),SO21997
       close(101)
 

      END PROGRAM TESTEPPA
