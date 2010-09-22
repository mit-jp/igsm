
!       ==========================================================
!
!       DRIVER.F:  The main program of the MIT Global Climate and
!                       Biogeochemistry Model IGSM2.2
!
!       ----------------------------------------------------------
!
!       This version calls Climate-Chemistry Model (subroutine atmosphere)
!               and TEM
!
!               MIT Joint Program for the Science and Policy
!                       of Global Change
!
!       ----------------------------------------------------------

C     === Local variables ===
      INTEGER iloop
      INTEGER nCouplePer
      INTEGER myIter, myThid   !for compatibility with MITGCM ocean
      REAL*8  myTime 
      CHARACTER *10 buf
      CHARACTER *10 buf1

      myTime=0.D0      !variables not used here
      myIter=0

      CALL date(buf)
      CALL time(buf1)
      print *,'Run started at ',buf1,'on ',buf

      CALL ATM2D_INIT_FIXED( nCouplePer, myThid)
      CALL ATM2D_INIT_VARS( myThid)

      CALL date(buf)
      CALL time(buf1)
      print *,'Initialization finished at ',buf1,'on ',buf

      DO iloop=1, nCouplePer
        
        CALL FORWARD_STEP_ATM2D(iloop, myTime, myIter, myThid)

      ENDDO !  End of time loop

C      CALL DRIVER_FINISH   ! not doing anything

      CALL date(buf)
      CALL time(buf1)
      print *,'Run ended at ',buf1,' on', buf
     
      STOP
      END
