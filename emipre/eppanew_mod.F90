!------------------------------------------------------------------------------
!          MIT INTEGRATED GLOBAL SYSTEMS MODEL                                 !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: eppanew_mod.F90
!
! !DESCRIPTION: This module eventually will call eppa every 5 years 
!   and puts emissions into arrays useable for the IGSM . 
!\\
!\\
! !INTERFACE: 
!
   MODULE EPPANEW_MOD
! 
! !USES:
!
  
  IMPLICIT NONE
!!PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: EPPA_DRIVER
!
! !PRIVATE MEMBER FUNCTIONS:
!
!
! !REVISION HISTORY:
!  17Jul09: Noelle Selin
!
! !REMARKS:
!  Right now, this only does preprocessing. Later it will include 
!  EPPA call
!
!  Input: a file from EPPA which contains
!  ag and non-ag emissions of:
!  nox, so2, co, voc, bc, oc, amo, ch4, n2o, co2 
!  called chm.put
!  Outputs: EMISTEMP arrays that can be used in the IGSM 
!           One for each species
!  Based in part from edaily.f, emipop122207.f90 by CWang+MSarofim
!
!EOP
!------------------------------------------------------------------------------
! 
!  ! PRIVATE DATA MEMBERS:

      include "size.h"   !! Size parameters      
   !MODULE VARIABLES
      REAL                  :: EMISTEMPCO2(NLON,NLAT)
      REAL                  :: EMISTEMPCO(NLON, NLAT)
      REAL                  :: EMISTEMPSO2(NLON, NLAT)
      REAL                  :: EMISTEMPBC(NLON, NLAT)
      REAL                  :: EMISTEMPOC(NLON, NLAT)
      REAL                  :: EMISTEMPNO2(NLON, NLAT)
      REAL                  :: EMISTEMPNH3(NLON, NLAT)
      REAL                  :: EMISTEMPVOC(NLON, NLAT)
      REAL                  :: EMISTEMPN2o(NLON, NLAT)
      REAL                  :: EMISTEMPCH4(NLON, NLAT)
      REAL                  :: EMISTEMPc12(NLON, NLAT)
      REAL                  :: EMISTEMPC11(NLON, NLAT)
      REAL                  :: EMISTEMPhfc(NLON, NLAT)
      REAL                  :: EMISTEMPpfc(NLON, NLAT)
      REAL                  :: EMISTEMPsf6(NLON, NLAT)
      REAL                  :: SO21997

    !URBANS
      REAL                  :: EMISTEMPuno2(NLON, NLAT)
      REAL                  :: EMISTEMPuso2(NLON, NLAT)
      REAL                  :: EMISTEMPuco(NLON, NLAT)
      REAL                  :: EMISTEMPunmv(NLON, NLAT)
      REAL                  :: EMISTEMPubc(NLON, NLAT)
      REAL                  :: EMISTEMPunh3(NLON, NLAT)
      REAL                  :: EMISTEMPuOC(NLON, NLAT)

      INTEGER                  :: N_URBAN(3, NLON, NLAT)
      INTEGER               :: N_TOTAL_URBAN(NLON, NLAT)

      REAL                  :: REGIONMAP(GLON, GLAT, LREGIONS)
      REAL                  :: CITYINDEX(GLON, GLAT)

!APS
       CHARACTER(LEN=255)           :: GRIDS
       CHARACTER(LEN=255)           :: FILEIN
!APS

       CONTAINS

!------------------------------------------------------------------------------
!          MIT INTEGRATED GLOBAL SYSTEMS MODEL                                 !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  EPPA_DRIVER
!
! !DESCRIPTION: This is the driver routine for EPPA in the IGSM. 
! Should be  called from the coupler for each year.
! Input is year in integer format. 
!\\
!\\
! !INTERFACE: 
!
  SUBROUTINE EPPA_DRIVER(YEAR)
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN) :: year  
!   CHARACTER(LEN=255), INTENT(IN) :: GRIDS  
!
! !REVISION HISTORY: 
!  17Jul09: Noelle Selin
!  20Sep10: Now saves SO2 1997 value for IGSM input per A. Sokolov;
!  automatically reads EPPA4 or EPPA5
! !REMARKS:
! 
!
!EOP
!------------------------------------------------------------------------------
!BOC
   !LOCAL VARIABLES
      LOGICAL, SAVE         :: FIRST_EPPA = .TRUE. 
      CHARACTER(LEN=3),SAVE :: EPPAREG(LREGIONS+1)
      CHARACTER(LEN=3)      :: SPECIES
      REAL                  :: OUTEMISYEAR1(LSPECIES, LREGIONS)
      REAL                  :: OUTEMISYEAR2(LSPECIES, LREGIONS)
      REAL                  :: OUTEMISS(LSPECIES, LREGIONS)
      INTEGER               :: EPPAYEAR
      INTEGER               :: YEARINDEX, Y, L1, L2
      INTEGER,SAVE          :: EPPAYEARS(LYEAR)
      REAL                  :: TESTYEAR, TOTALSF6
      REAL                  :: TEMP2(GLON, GLAT)
      REAL                  :: YR1, YR2, X
      REAL, SAVE            :: INPEMISS(LSPECIES, LREGIONS+1, LYEAR)
      CHARACTER(LEN=5), SAVE:: EPPAFLAG

      ! First call read_regions and figure out if it's EPPA4 or EPPA5
       IF (FIRST_EPPA) THEN
!         CALL EPPA_READ(1997, EPPAREG, OUTEMISYEAR1, EPPAFLAG)
          CALL EPPA_READ( EPPAREG, INPEMISS, EPPAFLAG)
          CALL EPPA_GET( 1, OUTEMISYEAR1,INPEMISS)
          CALL READ_REGIONS(EPPAREG, EPPAFLAG)
          SO21997=SUM(OUTEMISYEAR1(13,:)+OUTEMISYEAR1(20,:))+25.6 !add natural emissions
          !PRINT*, SO21997
           DO Y = 1,LYEAR
             EPPAYEARS(Y)=INPEMISS(1,1,Y)
!            WRITE(6,*), EPPAYEARS(Y)
           ENDDO
       
          FIRST_EPPA = .FALSE.
       ENDIF 
       
      !Pre-1997 is the same for both EPPA4 and EPPA5
       IF (YEAR .ge. 1991 .and. YEAR .lt. 1997) THEN
       ! read EPPA 1997 values
       

!      CALL EPPA_READ(1997, EPPAREG, OUTEMISYEAR1, EPPAFLAG)
       CALL EPPA_GET( 1, OUTEMISYEAR1,INPEMISS)

       OUTEMISS=OUTEMISYEAR1
       
       CALL PRE_EPPA_EMIS(OUTEMISS, YEAR) 
     
       !Return if we aren't at 1991 yet
       
       
       ENDIF
       
       IF (YEAR .ge. 1997) THEN
            
           DO Y = 1,LYEAR
             IF( YEAR.ge.EPPAYEARS(Y)) then
              L1=Y
              YR1=EPPAYEARS(Y)
              IF(YEAR.lt.2100) THEN
               L2=Y+1
               YR2=EPPAYEARS(Y+1)
              ENDIF
             ENDIF
           ENDDO
        
             WRITE(6,*), YEAR,YR1,YR2

           CALL EPPA_GET( L1, OUTEMISYEAR1,INPEMISS)
           IF(YEAR.lt.2100) THEN
             X=(YEAR-YR1)/(YR2-YR1)
             CALL EPPA_GET( L2, OUTEMISYEAR2,INPEMISS)
          !Here, interpolate between the two years
          ! save the resulting interpolation in outemiss
             OUTEMISS=OUTEMISYEAR1*(1.-X)+OUTEMISYEAR2*X
           ELSE
          !At EPPAYEAR=2100, don't try to read 2105! 
             OUTEMISS=OUTEMISYEAR1
           ENDIF

          !Call the program to update the emistemp variables
          CALL EPPA_EMIS(OUTEMISS, YEAR)
             
       
 
       END IF

       call write_output(YEAR)
       !print*, emistempco2
      ! print*, emistempso2
     !  write(7), emistempso2

 END SUBROUTINE EPPA_DRIVER
!EOC
!------------------------------------------------------------------------------
!          MIT INTEGRATED GLOBAL SYSTEMS MODEL                                 !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  EPPA_EMIS
!
! !DESCRIPTION: Does emissions for EPPA YEARS
!\\
!\\
! !INTERFACE:
!
         SUBROUTINE EPPA_EMIS(OUTEMISS, YEAR)
!
! !INPUT PARAMETERS: 
!
       REAL, INTENT(IN) :: OUTEMISS(:,:)
       INTEGER, INTENT(IN) :: YEAR
!
! !OUTPUT PARAMETERS 
!
!
! !REVISION HISTORY: 
!  30 Jul09: Noelle Selin
!
! !REMARKS:
!
!EOP
!------------------------------------------------------------------------------
!BOC
     !Local variables

      REAL                  :: DUMMYEMISS(LSPECIES, LREGIONS)
      dummyemiss=0d0
          !Here, call all the emissions which are gridded according to 
          !Edgar ag and nonag, first the noncarbon species
          !SO2 ag = 13, so2 non-ag =20

          !Need to call NOx first, because it updates the urban cityindex
          CALL EMIS_EPPA_NONC('nox',YEAR,OUTEMISS(12,:),OUTEMISS(19,:), emistempno2)
       !   print*, sum(outemiss(12,:)+outemiss(19,:))
       !   stop
          CALL EMIS_EPPA_NONC('so2',YEAR,OUTEMISS(13,:),OUTEMISS(20,:), emistempso2)
        
        
       
          !NO2 ag=12, no2 non-ag=19

          !bc ag=14, bc non-ag=21 [grid according to sox for now]
          CALL EMIS_EPPA_NONC('bc1',YEAR,OUTEMISS(14,:),OUTEMISS(21,:), emistempbc)
    
           !oc ag=15, oc non-ag=22 [grid according to sox for now]
          CALL EMIS_EPPA_NONC('oc1',YEAR,OUTEMISS(15,:),OUTEMISS(22,:), emistempoc)
       
            !n20 ag=4, n2o non-ag=6
          CALL EMIS_EPPA_NONC('n2o',YEAR, OUTEMISS(4,:),OUTEMISS(6,:), EMISTEMPN2o)
      
          !nmvoc ag=11, nmvoc nonag=18
          CALL EMIS_EPPA_NONC('nmv',YEAR,OUTEMISS(11,:),OUTEMISS(18,:), EMISTEMPVOC)
      
          !NH3 ag = 16, nh3 non-ag=23 [grid by nox for now]
          CALL EMIS_EPPA_NONC('nh3',YEAR,OUTEMISS(16,:),OUTEMISS(23,:), EMISTEMPNH3)
      
         
       CALL EMIS_EPPA_NONC('sf6',YEAR,DUMMYEMISS(9,:),OUTEMISS(8,:), EMISTEMPsf6)
       !   print*, 'total sf6', sum(emistempsf6)

          
       CALL EMIS_EPPA_NONC('pfc',YEAR,DUMMYEMISS(9,:),OUTEMISS(7,:), EMISTEMPpfc)
      
          !Now CO2, CH4, methane together
          CALL EMIS_EPPA_CARBON(OUTEMISS, YEAR, EMISTEMPCO2, EMISTEMPCO, EMISTEMPCH4)
     

 
          ! CFC11, CFC12, HFC, SF6
          !only HFCs need the EPPA values, so otherwise send a zero value
          DUMMYEMISS=0d0
          CALL EMIS_CFCs('c12', YEAR, DUMMYEMISS(9,:))
       !   print*, 'total cfc12', sum(emistempc12)
          CALL EMIS_CFCs('c11', YEAR, DUMMYEMISS(9,:))
        !  print*, 'total cfc11', sum(emistempc11)
          CALL EMIS_CFCs('sf6', YEAR, DUMMYEMISS(9,:))
      !    print*, 'total sf6', sum(emistempsf6)
          ! for HFC134a use eppa
          CALL EMIS_CFCs('134', YEAR, OUTEMISS(9,:))
       !    print*, 'total hfc', sum(emistemphfc)
        
          
      END SUBROUTINE EPPA_EMIS
!EOC
!------------------------------------------------------------------------------
!          MIT INTEGRATED GLOBAL SYSTEMS MODEL                                 !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  PRE_EPPA_EMIS
!
! !DESCRIPTION: Does emissions for 1991-1997
!\\
!\\
! !INTERFACE:
!
         SUBROUTINE PRE_EPPA_EMIS(OUTEMISS, YEAR)
!
! !INPUT PARAMETERS: 
!
       REAL, INTENT(IN) :: OUTEMISS(:,:)
       INTEGER, INTENT(IN) :: YEAR
!
! !OUTPUT PARAMETERS 
!
!
! !REVISION HISTORY: 
!  30 Jul09: Noelle Selin
!
! !REMARKS:
!
!EOP
!------------------------------------------------------------------------------
!BOC
     !Local variables

      REAL                  :: DUMMYEMISS(LSPECIES, LREGIONS)
      REAL                  :: NATTEMP(GLON, GLAT)
      REAL                  :: NATTEMP2(NLON, NLAT)
      REAL                  :: SO2FACTOR
         DUMMYEMISS=0d0
!For NOx, linearly interpolate from 1977=9.8135e4 as NO to 1997
          CALL EMIS_EPPA_NONC('nox',1997,OUTEMISS(12,:),OUTEMISS(19,:), emistempno2)
          EMISTEMPNO2=((9.8135e4)+(sum(emistempno2)-9.8135e4)/20.*(year-1977))*emistempno2/sum(emistempno2)
          
   
       !   print*, 'total nox in', sum(outemiss(12,:)+outemiss(19,:))
! For BC, OC, NH3 use 1997 values
          !bc ag=14, bc non-ag=21 [grid according to sox for now]
          CALL EMIS_EPPA_NONC('bc1',1997, OUTEMISS(14,:),OUTEMISS(21,:), emistempbc)
    
           !oc ag=15, oc non-ag=22 [grid according to sox for now]
          CALL EMIS_EPPA_NONC('oc1',1997,OUTEMISS(15,:),OUTEMISS(22,:), emistempoc)
    
         !NH3 ag = 16, nh3 non-ag=23 [grid by nox for now]
          CALL EMIS_EPPA_NONC('nh3',YEAR,OUTEMISS(16,:),OUTEMISS(23,:), EMISTEMPNH3)
   
 
          CALL EMIS_EPPA_NONC('sf6', YEAR,DUMMYEMISS(9,:),OUTEMISS(8,:), EMISTEMPsf6)
   
           CALL EMIS_EPPA_NONC('pfc',YEAR,DUMMYEMISS(9,:),OUTEMISS(7,:), EMISTEMPpfc)
     
            EMISTEMPPFC=EMISTEMPPFC*(15./SUM(EMISTEMPPFC))

!  FOR CFCs, call with year
          dummyemiss=0d0
          CALL EMIS_CFCs('c12', YEAR, DUMMYEMISS(9,:))
     
          CALL EMIS_CFCs('c11', YEAR, DUMMYEMISS(9,:))
      
          CALL EMIS_CFCs('sf6', YEAR, DUMMYEMISS(9,:))
  
          CALL EMIS_CFCs('134', YEAR, OUTEMISS(9,:))
         
        


! for N20, linearly interpolate from 1987=0
           CALL EMIS_EPPA_NONC('n2o',YEAR, OUTEMISS(4,:),OUTEMISS(6,:), EMISTEMPN2o)
          emistempn2o=(sum(emistempn2o)/9.*(year-1987))*emistempno2/sum(emistempno2)

! for CO2, linearly interpolate from 1977=1.76e7
! need to call EMIS_EPPA_CARBON to subtract out CO, CH4
          CALL EMIS_EPPA_CARBON(OUTEMISS,YEAR, EMISTEMPCO2, EMISTEMPCO, EMISTEMPCH4)
      
          EMISTEMPCO2=((1.76e7)+(sum(emistempco2)-1.76e7)/20.*(year-1977))*emistempco2/sum(emistempco2)

          EMISTEMPCH4=((3.09e5)+(sum(emistempch4)-3.09e5)/20.*(year-1977))*emistempch4/sum(emistempch4)
 
          EMISTEMPCO=((1.33e6)+(sum(emistempco)-1.33e6)/20.*(year-1977))*emistempco/sum(emistempco)


       
! For SO2, exrapolate back using factors (Sokolov)
          !SO2 ag = 13, so2 non-ag =20
          CALL EMIS_EPPA_NONC('so2',YEAR,OUTEMISS(13,:),OUTEMISS(20,:), emistempso2)
          ! put in if statements here that scale emistempso2 by Sokolov factors
              IF (YEAR == 1991) THEN 
             EMISTEMPSO2=EMISTEMPSO2*1.073
             ELSE IF (YEAR == 1992) THEN 
                EMISTEMPSO2=EMISTEMPSO2*1.047
             ELSE IF (YEAR == 1993) THEN 
               EMISTEMPSO2=EMISTEMPSO2*1.031
             ELSE IF (YEAR == 1994) THEN 
               EMISTEMPSO2=EMISTEMPSO2*1.027
             ELSE IF (YEAR == 1995) THEN 
                EMISTEMPSO2=EMISTEMPSO2*1.020
             ELSE IF (YEAR == 1996) THEN 
                EMISTEMPSO2=EMISTEMPSO2*1.011
             END IF

          
      END SUBROUTINE PRE_EPPA_EMIS
!EOC
!------------------------------------------------------------------------------
!          MIT INTEGRATED GLOBAL SYSTEMS MODEL                                 !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  EPPA_READ
!
! !DESCRIPTION: This subroutine reads a chm.put file from EPPA. 
! It updates OUTEMISS (nspecies[incl ag and nonag separately]*nregions)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EPPA_READ( EPPAREG, INPEMISS, EPPAFLAG)
!
!
! !OUTPUT PARAMETERS:
    !character string of EPPA region names
    CHARACTER(LEN=3), INTENT(OUT) :: EPPAREG(LREGIONS+1)
    REAL, INTENT(OUT)         :: INPEMISS(LSPECIES, LREGIONS+1, LYEAR)
    CHARACTER(LEN=5), INTENT(OUT)        :: EPPAFLAG
! !REVISION HISTORY: 
!  17Jul09: Noelle Selin
!
! !REMARKS:
! This code is a bit messy because of the structure of the EPPA file. Any changes
! to the EPPA output file structure will necessitate edits in this code.
!
!EOP
!------------------------------------------------------------------------------
!BOC
!  !LOCAL VARIABLES:
       !species x region+1 x lyears    
       CHARACTER(LEN=255)      :: FILENAME
       REAL                    :: TEMP, RYEAR
       INTEGER                 :: IUNIT, IOS, R, J, Y, L
       INTEGER                 :: YEARINDEX
       CHARACTER(LEN=10)       :: JUNK

       IUNIT=39
       
    

!      FILENAME = 'eppa5chm.put'
       FILENAME = TRIM(FILEIN)

       ! Open the CHM.PUT file here
       OPEN(IUNIT, FILE=FILENAME, STATUS='OLD', IOSTAT=IOS)

       !Error check!
       IF ( IOS /= 0 )  WRITE( 6, '(/,a)' ) 'Error: CHM.PUT not found'

!    Original code below from Marcus
!	read the emissions data for CO2, CH4, N2O, NOx, SOx, CO, NMVOC, 
!       HFC, PFC, SF6,BC,OC
!	agricultural data and non-agric data (ff) -> from EPPA
!	Order of EPPA Regions: 
!       (1)USA  (2)CAN  (3)MEX  (4)JPN  (5)ANZ  (6)EUR (7)EET (8)FSU  (9)ASI  
!       (10) CHN  (11)IND  (12)IDZ (13)AFR  (14)MES  (15)LAM  (16)ROW
!	
!	order of species as written into chm.put from EPPA:
!	(1)CO2 ag, (2)CO2 non-ag, (3)CH4ag, (4)N2O ag, (5)CH4 non-ag,(6)
!       N2O non-ag,(7) PFC,
!	(8) SF6, (9) HFC, (10)CO ag, (11)NMV ag, (12)NO2 ag, (13)SO2 ag,
!       (14)BC ag, (15) OC ag,
!	(16) NH3 ag, (17) CO non-ag, (18)NMV non-ag, (19)NO2 non-ag,
!       (20)SO2 non-ag, 
!	(21)BC non-ag, (22)OC non-ag, (23)NH3 non-ag
!	  
!	units: Tg C for CO2, Tg CH4, Tg N2O, Tg CO, Tg NO2, Tg SO2, 
!       Tg NMVOC, Tg BC, Tg OC, 
!	Tg NH3,	PFC: Mg CF4 equivalents, HFC: Mg HFC-134a equivalents, Mg SF6
        

        ! This will be put into the inpemiss array.CHM.PUT
        ! Inpemiss is species x region+1 x lyears
        ! in the below, loops over j,r,y refer to species, region, year
        ! l loops over blank lines

        ! First, loop over the first 11 rows of eppa blank lines
         READ(IUNIT,*) (JUNK, eppaflag)
         ! print*, eppaflag
        !  stop

	DO L = 2,11
	   READ(IUNIT,*) 
	ENDDO
      
    
        ! Then read the EPPA region text. There are 17 here because the first is the year column title	
      
	READ(IUNIT,*) (EPPAREG(R), R=1,LREGIONS+1)

        ! Then you finally get down to the emissions in the chm.put file. Read them in turn.
     
        ! Reading CO2-AG, loop over years and regions
	DO Y = 1,LYEAR
           
             
		READ(IUNIT,*) (INPEMISS(1,R,Y), R=1,LREGIONS+1)
         
	ENDDO
       	 
        ! Skip 2 lines after CO2_ag
	READ(IUNIT,*)
	READ(IUNIT,*)
       
        ! Reading CO2 non-ag, loop over years and regions
	DO Y = 1,LYEAR
              
		READ(IUNIT,*) (INPEMISS(2,R,Y), R=1,LREGIONS+1)

	ENDDO
  
        ! Skip 2 lines after CO2_nonag
	DO L = 1,5	
		READ(IUNIT,*)
	ENDDO
        
        ! Reading CH4 and N20 ag, loop over years and regions
	DO J = 3,4
	READ(IUNIT,*)
	  DO Y = 1,LYEAR
               
		READ(IUNIT,*) (INPEMISS(J,R,Y), R=1,LREGIONS+1)
          ENDDO
	ENDDO

        ! Skip many lines
	DO L=1,3*(LYEAR+1)+2
		READ(IUNIT,*)
	ENDDO

        !  Read (5)CH4 non-ag,(6) N2O non-ag,(7) PFC, (8) SF6, (9) HFC
        ! loop over years and regions
	DO J = 5,9
	READ(IUNIT,*)
	  DO Y = 1,LYEAR
              
		READ(IUNIT,*) (INPEMISS(J,R,Y), R=1,LREGIONS+1)

	  ENDDO
	ENDDO
	
        !Skip many more lines
	DO L = 1,5
		READ(IUNIT,*)
	ENDDO

        !(10)CO ag, (11)NMV ag, (12)NO2 ag, (13)SO2 ag,
        !(14)BC ag, (15) OC ag,
        !(16) NH3 ag,
        ! loop over years and regions
	DO J = 10,16
	READ(IUNIT,*)
	   DO Y = 1,LYEAR
               
		READ(IUNIT,*) (INPEMISS(J,R,Y), R=1,LREGIONS+1)

	  ENDDO
	ENDDO

        !Skip many more lines
	READ(IUNIT,*)
	READ(IUNIT,*)
 
        ! (17) CO non-ag, (18)NMV non-ag, (19)NO2 non-ag, 
        ! (20)SO2 non-ag, 
        ! (21)BC non-ag, (22)OC non-ag, (23)NH3 non-ag

	DO J = 17,23
        READ(IUNIT,*)
	   DO Y = 1,LYEAR
               
		READ(IUNIT,*) (INPEMISS(J,R,Y), R=1,LREGIONS+1)

	  ENDDO
	ENDDO       
        
       CLOSE( IUNIT )
       
       
       ! Regurn to calling program
       END SUBROUTINE EPPA_READ
!EOC
!BOP
!
! !IROUTINE:  EPPA_GET
!
! !DESCRIPTION: This subroutine pull relevant year from INPEMISS
! It updates OUTEMISS (nspecies[incl ag and nonag separately]*nregions)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EPPA_GET(YEARINDEX, OUTEMISS, INPEMISS)
!
! !INPUT PARAMETERS:
!
    ! year
    INTEGER, INTENT(IN)     :: YEARINDEX
    REAL, INTENT(IN)         :: INPEMISS(LSPECIES, LREGIONS+1, LYEAR)
!
!
! !OUTPUT PARAMETERS:
    !character string of EPPA region names
    REAL, INTENT(OUT)         :: OUTEMISS(LSPECIES, LREGIONS)

!  !LOCAL VARIABLES:
      INTEGER                 ::R, J

        !Pull only the relevant year from INPEMISS
        DO J = 1,LSPECIES
        DO R = 1,LREGIONS
        OUTEMISS(J,R)=INPEMISS(J,R+1,YEARINDEX)
        ENDDO
        ENDDO
       ! Regurn to calling program
       END SUBROUTINE EPPA_GET
     
!------------------------------------------------------------------------------
!          MIT INTEGRATED GLOBAL SYSTEMS MODEL                                 !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  READ_EDGAR_FLAG
!
! !DESCRIPTION: This subroutine reads EDGAR data given species, flag and 
! a map. The EDGAR data is used for the gridding initial map.
!\\
!\\
! !INTERFACE:
!
       SUBROUTINE READ_EDGAR_FLAG(SPEC, FLAG, MAP)
!
! !INPUT PARAMETERS: 
!
!   SPECIES
      CHARACTER(LEN=3), INTENT(IN) :: SPEC
      CHARACTER(LEN=3), INTENT(IN) :: FLAG
!
! !OUTPUT PARAMETERS 
!
      REAL, INTENT(OUT)            :: MAP(GLON,GLAT)
!
! !REVISION HISTORY: 
!  17Jul09: Noelle Selin
!
! !REMARKS:
!  The flag indicates how the species should be gridded. You can modify and add
!  another gridding scheme by adding another flag to the code. 
!  This routine calls READ_EDGAR_FILE to do the actual reading.
!EOP
!------------------------------------------------------------------------------
!BOC
! Local Variables: 
      CHARACTER(LEN=255)           :: FILENAME, FOLDER
      INTEGER                      :: IUNIT
      INTEGER                      :: IOS3,  N, K

      CHARACTER(LEN=3)             :: NONAGLIST(27) 
         CHARACTER(LEN=3)             :: HFCLIST(5) 
      CHARACTER(LEN=3)             :: AGLIST(13), TOTLIST(40), biolist(4)
      
      !Initialize AG_EMIS
      
      IUNIT=3
      
      FOLDER=TRIM(GRIDS)//'edgar/' &
      //'edgar_32ft2000_'//SPEC//'/'
     
   !   if  (spec .eq. 'co2') then
   !   FOLDER='grids/edgar/edgar_3_co2/'
   !   print*, 'testco2'
   !   end if

      !Here, you could put some if statements to say which ag and nonaglist to read depending on
      !what species it is.
      IF (FLAG .eq. 'nag') THEN
      NONAGLIST=(/'b10', 'b20', 'b30', 'b40', 'b51','f10', 'f20', 'f30', 'f40', 'f51', 'f54', &
           'f57', 'f58', 'f60','f70','f80','f90', 'i10', 'i20', 'i30', 'i40', 'i50','i70','w10',&
           'w20','w30','w40'/)

       !Loop over the non-ag categories in EDGAR 
      DO N=1, 27
  
      FILENAME=trim(FOLDER)//nonaglist(n)//'00'//SPEC//'.1x1'
      CALL READ_EDGAR_FILE(FILENAME, MAP)
      
      ENDDO

      ELSE IF (FLAG .eq. 'agr') THEN
     !The below need to be checked against adding categories in the future
     AGLIST=(/'l10','l15','l20','l30','l41','l42','l43','l44','l47','l50','l60','l71','l75'/)
      
   
          !Loop over the ag categories in EDGAR 
      DO K=1,13
   
      !Open the file, and check to see if it exists
      
      FILENAME=trim(FOLDER)//aglist(k)//'00'//SPEC//'.1x1'
      CALL READ_EDGAR_FILE(FILENAME, MAP)
     enddo

      ELSE IF (FLAG .eq. 'hfn') THEN
     
      hfcLIST=(/'h32','h33','h34','h36', 'h38'/)
      
   
          !Loop over the ag categories in EDGAR 
      DO K=1,5
   
      !Open the file, and check to see if it exists
      
      FILENAME=trim(FOLDER)//hfclist(k)//'00'//SPEC//'.1x1'
      CALL READ_EDGAR_FILE(FILENAME, MAP)
    
       
      ENDDO

    !for gridding relative to total emissions
    ELSE IF (FLAG .eq. 'tot') THEN
      totLIST=(/'b10', 'b20', 'b30', 'b40', 'b51','f10', 'f20', 'f30', 'f40', 'f51', 'f54', &
           'f57', 'f58', 'f60','f70','f80','f90', 'i10', 'i20', 'i30', 'i40', 'i50','i70','w10',&
           'w20','w30','w40','l10','l15','l20','l30','l41','l42','l43','l44','l47','l50','l60','l71','l75'/)
       !Loop over the non-ag categories in EDGAR 
      DO N=1,40
   
      FILENAME=trim(FOLDER)//totlist(n)//'00'//SPEC//'.1x1'
      CALL READ_EDGAR_FILE(FILENAME, MAP)
      
      ENDDO
  

   !for gridding relative to biomass
    ELSE IF (FLAG .eq. 'bio') THEN
      bioLIST=(/'l41','l42','l43','l44'/)
       !Loop over the L categories in EDGAR 
      DO N=1,4
   
      FILENAME=trim(FOLDER)//biolist(n)//'00'//SPEC//'.1x1'
      CALL READ_EDGAR_FILE(FILENAME, MAP)
      
      ENDDO
   ENDIF
  
       ! Return to calling program
       END SUBROUTINE READ_EDGAR_FLAG
!EOC
!------------------------------------------------------------------------------
!          MIT INTEGRATED GLOBAL SYSTEMS MODEL                                 !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  READ_EDGAR_FILE
!
! !DESCRIPTION: This subroutine reads the ag and nonag distribution maps 
! from EDGAR data. It's called multiple times from READ_EDGAR_FLAG. 
!\\
!\\
! !INTERFACE:
!
     SUBROUTINE READ_EDGAR_FILE(FILENAME, EMIS)
!
! !INPUT PARAMETERS: 
!
!   FILENAME
      CHARACTER(LEN=255), INTENT(IN)  :: FILENAME
!
! !OUTPUT PARAMETERS 
!
      REAL, INTENT(OUT)            :: EMIS(GLON, GLAT)
!
! !REVISION HISTORY: 
!  17Jul09: Noelle Selin
!
! !REMARKS:
!  The EDGAR format is consistent through all emissions files. If for some reason 
!  that changes, you will need to update this routine. 
!EOP
!------------------------------------------------------------------------------
!BOC
! Local Variables
!  
      INTEGER                      :: IUNIT, LON, LAT, L, I, IOS3
      REAL                        :: VAL
        IUNIT=25
      !zero emis
     ! EMIS=0d0
       ! open the file
       OPEN(IUNIT, FILE=FILENAME, STATUS='OLD', &
            IOSTAT=IOS3)
       

      !If the file doesn't exist, don't read it!
      IF ( IOS3 == 0 )  THEN
   
        
         ! Skip 12 lines
         DO L = 1,12
	   READ(IUNIT,*) 
         ENDDO
      
         ! Read the rest of the file 
         DO 
           READ(IUNIT,*, iostat=i) LON, LAT, VAL
           ! if (val .lt. 0d0) then
            ! print*, 'lon', lon
            ! print*, 'lat', lat
            ! print*, 'val', val
            ! print*, 'filename', filename
            ! stop
            ! end if
           EMIS(LON+181,LAT+91)=EMIS(LON+181,LAT+91)+VAL
           IF(i /= 0)  EXIT
         ENDDO
          
        
     ENDIF
         CLOSE(IUNIT)
         
       ! Return to calling program
       END SUBROUTINE READ_EDGAR_FILE
!EOC
!------------------------------------------------------------------------------
!          MIT INTEGRATED GLOBAL SYSTEMS MODEL                                 !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  EMIS_EPPA_NONC
!
! !DESCRIPTION: Does emissions for noncarbon gases that use EPPA emissions
!\\
!\\
! !INTERFACE:
!
         SUBROUTINE EMIS_EPPA_NONC (SPECIN, YEAR, EMIS_AG, EMIS_NONAG, EMISSPEC)

!
! !INPUT PARAMETERS: 
!
       REAL, INTENT(IN) :: EMIS_AG(:), EMIS_NONAG(:)
       CHARACTER(LEN=3), INTENT(IN) :: SPECIN
       INTEGER, INTENT(IN) :: YEAR
!
! !OUTPUT PARAMETERS 
!
      REAL, INTENT(OUT):: EMISSPEC(NLON, NLAT)
!
! !REVISION HISTORY: 
!  17Jul09: Noelle Selin
!
! !REMARKS:
! Carbon and non-carbon gases are done separately because of the need to carbon balance.
! For simplicity's sake, gridding is done by SO2 for SO2, BC, and OC and by NOX for NH3
! This can easily be changed by adding info to READ_EDGAR_FLAG and additional calls here
!EOP
!------------------------------------------------------------------------------
!BOC
! Local Variables
         REAL             :: AG_EMIS(GLON, GLAT)
         REAL             :: NONAG_EMIS(GLON, GLAT)
         REAL             :: REGTOT, TEMPAG(GLAT,GLON)
         REAL                  :: GLONTEMP2(GLON, GLAT)
         REAL             :: EPPAAG(GLON, GLAT), EPPANONAG(GLON, GLAT)
         REAL             :: EMISTEMP(GLON,GLAT)
         REAL             :: NATTEMP(GLON,GLAT)
         REAL             :: UTEMP(GLON,GLAT)
         REAL             :: THRESHOLD
         INTEGER          :: R, II, JJ
         CHARACTER(LEN=3) :: SPEC
        
         
        
         

         !Here and below, deal with the gridding

         EMISTEMP=0d0
         nattemp=0d0
         ! grid bc, oc as SO2 for now
         IF (specin .eq. 'bc1' .or. specin .eq. 'oc1') then
         spec='so2'
         ! grid nh3 by nox for now
         ELSE IF (specin .eq. 'nh3') THEN
         spec='nox'
         ELSE
         SPEC=SPECIN
         endif
  
         !Call read_edgar_maps
         !CALL READ_EDGAR_MAPS(SPEC, AG_EMIS, NONAG_EMIS)
         CALL READ_EDGAR_FLAG(SPEC, 'agr', AG_EMIS)
         CALL READ_EDGAR_FLAG(SPEC, 'nag', NONAG_EMIS)      
         !  write(101), ag_emis+nonag_emis
         !These functions distribute the EPPA emissions according to the EDGAR grids
         !CAll normalize_maps for AG EMISSIONS 
         CALL NORMALIZE_MAPS (EMIS_AG, AG_EMIS, EPPAAG)
     
        
         !CAll normalize_maps for NONAG EMISSIONS 
         CALL NORMALIZE_MAPS (EMIS_NONAG, NONAG_EMIS, EPPANONAG)
        
         !call add_natural here if > 1997
         IF (YEAR .ge. 1997) THEN
            CALL ADD_NATURAL(SPECIN, NATTEMP)
         ENDIF
        
        

         !Add up AG and NONAG emissions and regrid to IGSM2D; UNIT CONVERT TG to GG
         !add an if statement here for the unit conversion
         ! add natural emissions here also
         EMISTEMP = (EPPAAG + EPPANONAG) * 1000d0 + NATTEMP
    
         IF (SPEC .eq. 'sf6' .and. year .lt. 1997) then
            !SF6 from 0 to  6300 t/yr 1970-1996 and keep same rate
            CALL READ_EDGAR_FLAG('sf6','hfn', glontemp2)
            
            glontemp2= 0.24231*(year-1970)*glontemp2/sum(glontemp2)
         
           emistemp=glontemp2
           end if
           if (spec .eq. 'sf6' .and. year .ge. 1997) then
              emistemp=emistemp/1000.
            end if
         IF (SPEC .EQ. 'pfc') then
             emistemp=emistemp/1000.
             end if
    
         
         !********THIS IS WHAT WILL NEED TO BE CHANGED FOR 'NEW' METAMODEL'***********
         !Now, call the urban emissions before you grid to 2D.
         IF (specin .eq. 'nox') THEN
 

            !Change units to NO
             EMISTEMP=EMISTEMP*30./46.
         !>5 kg N per day per km2 
         ! x 95^2 km/gridsquare *365 days/1yr * 30 g NO/14 g N * 1e-6 Gg/kg
        ! THRESHOLD = 5.*95.*95.*365.*30./14.*1d-6

             !Test threshold=10 kg N to get 200 cities
              THRESHOLD = 10.*95.*95.*365.*30./14.*1d-6
       !  WHERE (EMISTEMP .ge. 35.294)
         WHERE (EMISTEMP .ge. THRESHOLD)
         CITYINDEX=1
         END WHERE
       !  print*, sum(cityindex)
        ! stop
         !N_TOTAL_URBAN is a 2D grid of cityindex
         
         N_TOTAL_URBAN=GRID_2D(CITYINDEX)
         CALL MAKE_NURBANS(YEAR, N_TOTAL_URBAN, N_URBAN)
         
          ENDIF
         
        
         
         ! multiply by cityindex to get urban, then grid2d to get EMISTEMP

         
         UTEMP=EMISTEMP*CITYINDEX
         IF (specin .eq. 'nox') THEN
         EMISTEMPuno2=grid_2D(UTEMP)
        ! write(81), utemp
         ELSE IF (specin .eq. 'so2') THEN
         EMISTEMPuso2 = GRID_2D(UTEMP)
       !  write(82), utemp
         ELSE IF (specin .eq. 'nmv') THEN
       !  write(83), utemp
         EMISTEMPunmv=GRID_2D(UTEMP)
         ELSE IF (specin .eq. 'bc1') THEN
         EMISTEMPubc=GRID_2D(UTEMP)
       !  write(84), utemp
         ELSE IF (specin .eq. 'nh3') THEN
         EMISTEMPunh3=GRID_2D(UTEMP)
       !  write(85), utemp
         ELSE IF (specin .eq. 'oc1') THEN
         EMISTEMPuoc=grid_2D(UTEMP)
       !  write(86), utemp
         END IF

         
         !********END OF WHAT WILL NEED TO BE CHANGED FOR 'NEW' METAMODEL'***********
         ! This calls the 1x1->2D gridding. If you want another grid write another
         ! function and call it here.
         EMISSPEC= GRID_2D(EMISTEMP)
 
       ! Return to calling program
       END SUBROUTINE EMIS_EPPA_NONC
!EOC
!------------------------------------------------------------------------------
!          MIT INTEGRATED GLOBAL SYSTEMS MODEL                                 !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  EMIS_EPPA_CARBON
!
! !DESCRIPTION: Does emissions for carbon gases
!\\
!\\
! !INTERFACE:
!
         SUBROUTINE EMIS_EPPA_CARBON(OUTEMISS, YEAR, EMISCO2, EMISCO, EMISCH4)

!
! !INPUT PARAMETERS: 
!
       REAL, INTENT(IN) :: OUTEMISS(:,:)
       INTEGER, INTENT(IN):: YEAR
!
! !OUTPUT PARAMETERS 
!
      REAL, INTENT(OUT):: EMISCO(NLON, NLAT), EMISCO2(NLON,NLAT), EMISCH4(NLON, NLAT)
!
! !REVISION HISTORY: 
!  17Jul09: Noelle Selin
!
! !REMARKS:
!
!EOP
!------------------------------------------------------------------------------
!BOC
! Local Variables
         REAL             :: AG_EMIS(GLON, GLAT)
         REAL             :: NONAG_EMIS(GLON, GLAT)
         REAL             :: REGTOT, TEMPAG(GLAT,GLON)
         REAL             :: EPPAAG(GLON, GLAT), EPPANONAG(GLON, GLAT)
         REAL             :: TEMPCO(GLON,GLAT), TEMPCO2(GLON,GLAT),TEMPCH4(GLON,GLAT)
         REAL             :: CH4ag(GLON,GLAT)
         REAL             :: NATTEMP(GLON,GLAT), nattemp2(glon, glat), nattemp3(glon, glat)
         INTEGER          :: R, II, JJ
         
         
        NATTEMP=0d0
       
    !CARBON MONOXIDE
         !Initialize ag_emis, nonag_emis
         AG_EMIS=0d0
         NONAG_EMIS=0d0
         !Call read_edgar_maps (CO is 'co1' for character length consistency)
         CALL READ_EDGAR_FLAG('co1', 'agr', AG_EMIS)
         CALL READ_EDGAR_FLAG('co1', 'nag', NONAG_EMIS)
                         
         !These functions distribute the EPPA emissions according to the EDGAR grids
         !CAll normalize_maps for AG EMISSIONS 
         CALL NORMALIZE_MAPS (OUTEMISS(10,:), AG_EMIS, EPPAAG)
         
         !CAll normalize_maps for NONAG EMISSIONS 
         CALL NORMALIZE_MAPS (OUTEMISS(17,:), NONAG_EMIS, EPPANONAG)
         
      
         !Add up AG and NONAG emissions and regrid to IGSM2D; UNIT CONVERT TG to GG
         !add an if statement here for the unit conversion
         ! add natural emissions here also
         TEMPCO = (EPPAAG + EPPANONAG) * 1000d0 
     
     !METHANE
             !Initialize ag_emis, nonag_emis
         AG_EMIS=0d0
         NONAG_EMIS=0d0
         !Call read_edgar_maps
         CALL READ_EDGAR_FLAG('ch4', 'agr', AG_EMIS)
         CALL READ_EDGAR_FLAG('ch4', 'nag', NONAG_EMIS)
                         
         !These functions distribute the EPPA emissions according to the EDGAR grids
         !CAll normalize_maps for AG EMISSIONS 
         CALL NORMALIZE_MAPS (OUTEMISS(3,:), AG_EMIS, EPPAAG)
         
         !CAll normalize_maps for NONAG EMISSIONS 
         CALL NORMALIZE_MAPS (OUTEMISS(5,:), NONAG_EMIS, EPPANONAG)
        
         !Add up AG and NONAG emissions and regrid to IGSM2D; UNIT CONVERT TG to GG
         !add an if statement here for the unit conversion
         ! add natural emissions here also
         TEMPCH4 = (EPPAAG + EPPANONAG) * 1000d0 
         
         !save methane ag separately for subtracting out of CO2
         ch4ag=eppaag*1000d0

      !CARBON DIOXIDE
             !Initialize ag_emis, nonag_emis
         AG_EMIS=0d0
         NONAG_EMIS=0d0
         !Call read_edgar_maps
        
         CALL READ_EDGAR_FLAG('co2','nag', NONAG_EMIS)
          CALL READ_EDGAR_FLAG('co2', 'agr', AG_EMIS)
                         
      !   write(21), ag_emis
     !    write(22), nonag_emis

         !These functions distribute the EPPA emissions according to the EDGAR grids
         !CAll normalize_maps for AG EMISSIONS 
         CALL NORMALIZE_MAPS (OUTEMISS(1,:), AG_EMIS, EPPAAG)
         
         !CAll normalize_maps for NONAG EMISSIONS 
         CALL NORMALIZE_MAPS (OUTEMISS(2,:), NONAG_EMIS, EPPANONAG)
         
        

         !Add up AG and NONAG emissions and regrid to IGSM2D; UNIT CONVERT TG to GG
         !convert units to Gg CO2 from Tg CO2
         !I think this is in Tg C: convert to CO2 here?
         ! add natural emissions here also
         TEMPCO2 = (EPPAAG + EPPANONAG) * 1000d0 * 44./12.

    
        
         
           nattemp=0d0
    
         
          nattemp2=0d0
          nattemp3=0d0
           
        !Adjust CO2 emissions by subtracting carbon emitted as CH4, CO
             CALL ADD_NATURAL('co1', NATTEMP2)
             TEMPCO=TEMPCO+NATTEMP2
              
              ! print*, 'co', sum(tempco)
            CALL ADD_NATURAL('ch4', NATTEMP3)
             TEMPCH4=TEMPCH4+NATTEMP3
                !   print*, 'co2', sum(tempco2)*12./44.
    

        !subtract natural ch4, agric ch4, and total co for carbon balancing
        TEMPCO2 = TEMPCO2 - (nattemp3*44./16.) - (ch4ag*44./16.) - (tempco*44./28.)
            !  print*, sum(tempco2)
           !   stop
       
        ! This calls the 1x1->2D gridding. If you want another grid write another
         ! function and call it here.
         EMISCO= GRID_2D(TEMPCO)
         EMISCO2= GRID_2D(TEMPCO2)
         EMISCH4= GRID_2D(TEMPCH4)
        
         !URBAN CO: NEEDS TO CHANGE FOR NEW METAMODEL *************
        EMISTEMPuco=grid_2d(TEMPCO*CITYINDEX)
     !  write(9), tempco
     !  write(10), tempco2
     !  write(11), tempch4
       ! Return to calling program
       END SUBROUTINE EMIS_EPPA_CARBON
!EOC
!------------------------------------------------------------------------------
!          MIT INTEGRATED GLOBAL SYSTEMS MODEL                                 !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  EMIS_CFCs
!
! !DESCRIPTION: Does emissions for CFCs and other gases
!\\
!\\
! !INTERFACE:
!
       SUBROUTINE EMIS_CFCS(SPEC, YEAR, HFCEMISS)

         !
         ! !INPUT PARAMETERS: 
         !
         INTEGER, INTENT(IN) :: YEAR
         REAL, INTENT(IN)    :: HFCEMISS(:)
         CHARACTER(LEN=3), INTENT(IN) :: SPEC
         !
         ! !OUTPUT PARAMETERS 
         !
         !
         ! !REVISION HISTORY: 
         !  17Jul09: Noelle Selin
         !
         ! !REMARKS:
         !
         !EOP
         !------------------------------------------------------------------------------
         !BOC
         ! Local Variables
         REAL                  :: GLONTEMP(GLON,GLAT)
         REAL                  :: GLONTEMP2(GLON,GLAT)
         REAL                  :: GLONTEMP3(GLON,GLAT)
         REAL                  :: HFCTEST(GLON, GLAT)
         CHARACTER(LEN=255)    :: FOLDER, CFCFILE, FILENAME2
         REAL                  :: TOTALc12, TOTALC11, TOTALHFC, HFCTOTAL(GLON, GLAT)

         glontemp=0d0
         glontemp2=0d0
         glontemp3=0d0
         !add stuff for CFCs -- needs to be scaled for different years
         !Folder where the cfc grids are
      !   folder='/Users/noelle/Documents/MIT07/igsm/newprocess/grids/edgar/'
        folder=TRIM(GRIDS)//'edgar/'

         IF (SPEC .eq. 'c11') THEN
            !CFC11
            !Use edgar emissions of cfc11 and scale according to extrapolation post 2000\
            cfcfile='hcgridtot/i918611.1x1'
            filename2=TRIM(FOLDER)//trim(cfcfile)
            CALL READ_EDGAR_FILE(FILENAME2, GLONTEMP3)
            EMISTEMPc11=GRID_2D(GLONTEMP3)
         !    print*, TRIM(FOLDER)//trim(cfcfile)
          !   print*, FILENAME2
           !  print*, sum(glontemp3)

         !   stop
            if (year .ge. 2010) THEN 
               TOTALC11 = 0d0
                ELSE IF (YEAR .eq. 2009) THEN
               TOTALc11 = 5d0
            ELSE IF (YEAR .eq. 2008) THEN
               TOTALc11 = 5d0
            ELSE IF (YEAR .le. 2007 .and. YEAR .ge. 2001) THEN
               TOTALc11=70-(YEAR-2000)*10+10
            ELSE IF (YEAR .eq. 2000) THEN
               TOTALc11=74.8
            ELSE IF (YEAR .eq. 1999) THEN
               TOTALc11=79.4
            ELSE IF (YEAR .eq. 1998) THEN
               TOTALc11=84.7
            ELSE IF (YEAR .eq. 1997) THEN
               TOTALc11=92.4
              ELSE IF (YEAR .eq. 1996) THEN
               TOTALc11=100.6
            ELSE IF (YEAR .eq. 1995) THEN
               TOTALc11=106.2
            ELSE IF (YEAR .eq. 1994) THEN
               TOTALc11=118.5
            ELSE IF (YEAR .eq. 1993) THEN
               TOTALc11=198.9
            ELSE IF (YEAR .eq. 1992) THEN
               TOTALc11=212.0
            ELSE IF (YEAR .eq. 1991) THEN
               TOTALc11=223.0

            ENDIF

            EMISTEMPc11=TOTALC11*EMISTEMPc11/sum(EMISTEMPc11)



         ELSE IF (SPEC .eq. 'c12') THEN 

            !CFC12
            !Use edgar emissions of cfc12 and scale according to extrapolation post 2000\
            cfcfile='hcgridtot/i918612.1x1'
             filename2=TRIM(FOLDER)//trim(cfcfile)
            CALL READ_EDGAR_FILE(FILENAME2, GLONTEMP)
            EMISTEMPc12=GRID_2D(GLONTEMP)

            if (year .ge. 2010) THEN 
               TOTALC12 = 0d0
             ELSE IF (YEAR .eq. 2009) THEN
               TOTALc12 = 5d0
            ELSE IF (YEAR .eq. 2008) THEN
               TOTALc12 = 5d0
            ELSE IF (YEAR .le. 2007 .and. YEAR .ge. 2004) THEN
               TOTALc12=35-(YEAR-2004)*5
            ELSE IF (YEAR .lt. 2004 .and. YEAR .ge. 2001) THEN
               TOTALc12=95-(YEAR-2001)*20
            ELSE IF (YEAR .eq. 2000) THEN
               TOTALc12=125.0
            ELSE IF (YEAR .eq. 1999) THEN
               TOTALc12=155.2
            ELSE IF (YEAR .eq. 1998) THEN
               TOTALc12=182.0
            ELSE IF (YEAR .eq. 1997) THEN
               TOTALc12=207.7
            ELSE IF (YEAR .eq. 1996) THEN
               TOTALc12=233.1
            ELSE IF (YEAR .eq. 1995) THEN
               TOTALc12=255.5
            ELSE IF (YEAR .eq. 1994) THEN
               TOTALc12=277.0
            ELSE IF (YEAR .eq. 1993) THEN
               TOTALc12=300.5
            ELSE IF (YEAR .eq. 1992) THEN
               TOTALc12=319.1
            ELSE IF (YEAR .eq. 1991) THEN
               TOTALc12=344.6

            ENDIF

            EMISTEMPc12=TOTALC12*EMISTEMPc12/sum(EMISTEMPc12)
            
            glontemp=0d0
         ELSE IF (SPEC .eq. '134') THEN
            !For hfc 134a use grid for cfc11 edgar but eppa emissions
            cfcfile='hcgridtot/i918611.1x1'
             filename2=TRIM(FOLDER)//trim(cfcfile)
            CALL READ_EDGAR_FILE(FILENAME2, HFCTEST)
             !    write(78), hfctest
             !     print*, 'beginHFC************'
            CALL NORMALIZE_MAPS(HFCEMISS,HFCTEST, GLONTEMP2)
             !     print*, hfcemiss
             !     print*, year, sum(glontemp2)
             !     print*, year, sum(hfcemiss)
            IF (year .ge. 1997) then
               
               EMISTEMPhfc=GRID_2D(GLONTEMP2)

            ELSE
                
               HFCTOTAL=glontemp2*(year-1991)/(1997-1991)
           
               EMISTEMPHFC=GRID_2D(HFCTOTAL)
            ENDIF

         ELSE
            RETURN
         ENDIF

         ! Return to calling program
       END SUBROUTINE EMIS_CFCs
!EOC
!-------------------------------if-----------------------------------------------
!          MIT INTEGRATED GLOBAL SYSTEMS MODEL                                 !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  READ_REGIONS
!
! !DESCRIPTION: This subroutine reads the regions for EPPA and stores 
! them in an array
! The array regionmap will have 1 if gridbox is in each region, 0 if not
! regionmap dimensions are longitude x latitude x nregions
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_REGIONS(EPPAREG, EPPAFLAG )
!
! !INPUT PARAMETERS: 
!
    ! EPPA REGION NAMES
    CHARACTER(LEN=3), INTENT(IN) :: EPPAREG(LREGIONS+1)    
    CHARACTER(LEN=5), INTENT(IN) :: EPPAFLAG
!
!
! !OUTPUT PARAMETERS:
    !map of regions with 1, 0 in each gridbox
! !REVISION HISTORY: 
!  17Jul09: Noelle Selin
!
! !REMARKS:
! This program is called only once per IGSM run.
! Right now, the input files are based on EPPA4.
! To change this to EPPA5 only the input files need to be changed.
! Region names are read from chm.put
! EPPAFLAG now indicates EPPA5 or EPPA4: nes (21 jun 2010)
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !Local Variables: 

       INTEGER                 :: R, IOS2, II, JJ
       REAL                    :: TEMPMAP(GLON,GLAT)
       REAL                    :: TEMPUNIT
       INTEGER, PARAMETER      :: UU=2
      
     
       ! Loop over each of the EPPA regions in order
       
       DO R=1, LREGIONS
       ! Open the 1x1 file here
       OPEN(UNIT=UU, FILE=TRIM(GRIDS)//EPPAFLAG//'/eppa.'//EPPAREG(R+1)//'.1x1',  FORM='UNFORMATTED',STATUS='OLD', &
         IOSTAT=IOS2)
       !Error check!
       !write(6,*), eppaflag
       IF ( IOS2 /= 0 )  WRITE( 6, '(/,a)' ) 'Error: Region file not found'
       !WRITE(6,*), 'grids/'//EPPAFLAG//'/eppa.'//EPPAREG(R+1)//'.1x1'
          READ(UU) TEMPMAP
          
          REGIONMAP(:,:,R)=TEMPMAP
       
       ENDDO

       !Return to calling program   
       END SUBROUTINE READ_REGIONS
!EOC
!------------------------------------------------------------------------------
!          MIT INTEGRATED GLOBAL SYSTEMS MODEL                                 !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  NORMALIZE_MAPS
!
! !DESCRIPTION:  !This subroutine takes as input the EDGAR maps of emissions
! It normalizes them so the total emissions in each region is equal to one
! for each of the EPPA regions. e 1 if gridbox is in each region, 0 if not
! regionmap dimensions are longitude x latitude x nregions
!\\
!\\
! !INTERFACE:
!
     SUBROUTINE NORMALIZE_MAPS(REGIONEMIS, EDGARMAP, EPPAMAP)
!
! !INPUT PARAMETERS: 
!
         REAL, INTENT(IN)      :: REGIONEMIS(LREGIONS)
         REAL, INTENT(IN)      :: EDGARMAP(GLON, GLAT)  
!
!
! !OUTPUT PARAMETERS:
          REAL, INTENT(OUT)     :: EPPAMAP(GLON, GLAT)
! !REVISION HISTORY: 
!  17Jul09: Noelle Selin
!
! !REMARKS:
! This program is called only once per IGSM run.
! Right now, the input files are based on EPPA4.
! To change this to EPPA5 only the input files need to be changed.
! Region names are read from chm.put
! EPPAFLAG now reads whether it's EPPA4 or EPPA5 (nes, 21 june 2010)
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !Local Variables: 
 
         REAL                  :: EMISTEMPR(GLON, GLAT)
         REAL                  :: REGTOT
         INTEGER               :: II, JJ, R

       !Loop over regions
       ! Initialize eppamap
        EPPAMAP=0d0
         DO R=1, LREGIONS
         EMISTEMPR=0d0
         REGTOT=0d0
       
         !Make the total EDGAR emissions in that region 1

          !I can't figure out how to make this work w/o a do loop...
         !Find the total emissions for each of the edgar regions
          DO JJ=1,GLON
         DO II=1,GLAT
         REGTOT= REGTOT+REGIONMAP(JJ,II,R) * EDGARMAP(JJ,II)
          ENDDO
         ENDDO
         
   
         !Normalize the emission map (AG_EMIS) for that region. 
         !This means dividing all the nonzero elements by the total
         !I can't figure out any way around a loop here to avoid the 
         !dividing by zero thing
         
         DO JJ=1,GLON
         DO II=1,GLAT
             
         IF (REGIONMAP(JJ,II,R) /= 0D0 .AND. EDGARMAP(JJ,II) /= 0D0) THEN
         
         EMISTEMPR(JJ,II)=EDGARMAP(JJ,II)/REGTOT
         
         ELSE
         EMISTEMPR(JJ,II)=0D0
         ENDIF
          ENDDO
         ENDDO
     !    print*, sum(emistempr)
         !Multiply EMISTEMPR by the regional emissions to get the emissionsmap
         !Add to the pregious eppamap
         EPPAMAP(:,:)=EPPAMAP(:,:)+(EMISTEMPR(:,:)*REGIONEMIS(R))
         
          ENDDO
      !   print*, 'TOTAL COMPARISON'
       !  print*, sum(eppamap)
        ! print*, sum(regionemis)
         !print*, sum(eppamap)-sum(regionemis)
        

       END SUBROUTINE NORMALIZE_MAPS

!EOC
!------------------------------------------------------------------------------
!          MIT INTEGRATED GLOBAL SYSTEMS MODEL                                 !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ADD_NATURAL
!
! !DESCRIPTION:  !This subroutine takes as input name of species
! then, it adds natural emissions to the emistemp file (gridded appropriately)
!\\
!\\
! !INTERFACE:
!
     SUBROUTINE ADD_NATURAL(SPEC, EMISTEMP)
!
! !INPUT PARAMETERS: 
!
         CHARACTER(LEN=3), INTENT(IN)      :: SPEC       
!
!
! !OUTPUT PARAMETERS:
         REAL, INTENT(OUT)      :: EMISTEMP(GLON, GLAT)
! !REVISION HISTORY: 
!  29Jul09: Noelle Selin
!
! !REMARKS:
! This program adds natural emissions appropriate to the grid
! Right now, it assumes that NEM is active. Need to add ifdef statements for
! Non-NEM runs.
! Need to decide how to grid each of these emissions (suggest EDGAR categories)
!EOP
!------------------------------------------------------------------------------
!BOC
! !Local Variables: 
 
         REAL                  :: EMISTEMPR(GLON, GLAT)
         REAL                  :: GRIDNAT(GLON, GLAT)
         REAL                  :: GRIDNAT2(GLON, GLAT)
         REAL                  :: landmap(glon, glat)
         INTEGER               :: R
   EMISTEMP=0d0
   IF (SPEC .eq. 'so2') THEN
   !Add 0.256e5, gridded according to agricultural emissions for now
   CALL READ_EDGAR_FLAG('so2','agr', GRIDNAT)
   EMISTEMP=0.256e5*(GRIDNAT)/sum(GRIDNAT)

   
   ELSE IF (SPEC .eq. 'nox') THEN
   !Add 2.143e4, gridded according to agricultural emissions (as NO2)
   CALL READ_EDGAR_FLAG('nox','agr', GRIDNAT)
   EMISTEMP=2.143e4*(GRIDNAT)/sum(GRIDNAT)*46./30.
 
  
   ELSE IF (SPEC .eq. 'co1') THEN
   !Add 370 Tg gridded according to agricultural emissions
   CALL READ_EDGAR_FLAG('co1','agr', GRIDNAT)
   EMISTEMP=0.370e6*(GRIDNAT)/sum(GRIDNAT)

   ELSE IF (SPEC .eq. 'ch4') THEN
   !Add 0.40e5 Gg gridded according to biomass burning emissions
   CALL READ_EDGAR_FLAG('ch4','bio', GRIDNAT)
!  EMISTEMP=0.40e5*(GRIDNAT)/sum(GRIDNAT)
!CLM35
!  EMISTEMP=0.00e5*(GRIDNAT)/sum(GRIDNAT)
   EMISTEMP=0.1e5*(GRIDNAT)/sum(GRIDNAT)
 

   ELSE IF (SPEC .eq. 'co2') THEN

   CALL READ_EDGAR_FLAG('co2','bio', GRIDNAT)
   !Forest regrowth
   EMISTEMP=-0.6e6*44./12.*(GRIDNAT)/sum(GRIDNAT)
   !Add ocean equally distributed for now
   !EMISTEMPR=(7.34e6/(glon*glat))
   !EMISTEMP=EMISTEMP+EMISTEMPR
   LANDMAP=0d0
   DO R=1, LREGIONS
   LANDMAP(:,:)=LANDMAP(:,:)+REGIONMAP(:,:,R)
   END DO 
   LANDMAP=(LANDMAP-1d0)*-1d0
   EMISTEMPR=(7.34e6/sum(landmap))*landmap
   EMISTEMP=EMISTEMP+EMISTEMPR
     ELSE
   EMISTEMP=EMISTEMP
   ENDIF
   

    END SUBROUTINE ADD_NATURAL
!EOC
!------------------------------------------------------------------------------
!          MIT INTEGRATED GLOBAL SYSTEMS MODEL                                 !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  MAKE_NURBANS
!
! !DESCRIPTION: This makes the nurban and nurbantotal arrays.
!  In the old metamodel these don't change. 
!\\
!\\
! !INTERFACE: 
!
  SUBROUTINE MAKE_NURBANS(year, n_total_urban, n_urban)
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: year
     INTEGER :: N_TOTAL_URBAN(NLON, NLAT)
     INTEGER, INTENT(OUT):: N_URBAN(3,NLON, NLAT)
 
     INTEGER :: J, II
!
! !REVISION HISTORY: 
!  17Jul09: Noelle Selin
!
! !REMARKS:
! 
!
!EOP
!------------------------------------------------------------------------------
!BOC
   !LOCAL VARIABLES
  
!
! --- urban data
!
    IF (YEAR .le. 1997) THEN
	! total
	n_total_urban(1, 1:nlat) = (/ 		        &
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4,	&
	12, 6, 4, 4, 0, 1, 9, 29, 14, 9, 3, 3, 3,	&
	6, 11, 18, 24, 53, 60, 54, 19, 20, 14, 5,	&	
	0, 0, 0, 0, 0, 0, 0, 0 /)
    END IF
      
 	! low    
	n_urban(1,1, 1:nlat)     = (/ 	                &
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,             &
	3, 3, 9, 5, 3, 3, 0, 1, 7, 22, 11, 7, 3, 3, 3,  &
	5, 8, 14, 18, 40, 45, 41, 14, 15, 11, 4, 0, 0,  &
	0, 0, 0, 0, 0, 0 /)

	! medium
	n_urban(2,1, 1:nlat)     = (/			&
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 3, 	&
	1, 1, 1, 0, 0, 2, 6, 3, 2, 0, 0, 0, 1, 3, 4, 	&
	5, 11, 12, 11, 5, 4, 3, 1, 0, 0, 0, 0, 		&
	0, 0, 0, 0 /)

	! high 
	n_urban(3,1, 1:nlat)     = (/  			&
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 			&
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,			&
	0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 		&
	3, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
     
        !This is code that calculates the number of urban regions of n=3 types, given 
        ! the total number of cities for the year
         ii=1
 	  DO j=1,nlat
	    IF (n_total_urban(ii,j).lt.4) THEN
		n_urban(1,ii,j) = n_total_urban(ii,j)
		n_urban(2,ii,j) = 0
		n_urban(3,ii,j) = 0
	    ELSE IF (n_total_urban(ii,j).lt.20) THEN
		n_urban(1,ii,j) = int(n_total_urban(ii,j)*0.75)
		n_urban(2,ii,j) = n_total_urban(ii,j) - n_urban(1,ii,j)
		n_urban(3,ii,j) = 0
	    ELSE
		n_urban(1,ii,j) = int(n_total_urban(ii,j)*0.75)
		n_urban(2,ii,j) = int(n_total_urban(ii,j)*0.20)
		n_urban(3,ii,j) = n_total_urban(ii,j) - n_urban(1,ii,j) - n_urban(2,ii,j)
	    END IF
	  END DO

 END SUBROUTINE MAKE_NURBANS
!EOC


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------ 
!          MIT INTEGRATED GLOBAL SYSTEMS MODEL                                 !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GRID_2D
!
! !DESCRIPTION: This function does a 2D grid
!\\
!\\
! !INTERFACE:
!
   FUNCTION GRID_2D (INGRID) RESULT (OUTGRID)
!
! !INPUT PARAMETERS: 
!
  REAL, INTENT(IN) :: INGRID(GLON, GLAT)
!
! !OUTPUT PARAMETERS:
!
      REAL :: OUTGRID(NLON,NLAT)
!
! !REVISION HISTORY: 
!  17July09: Noelle Selin
!
! !REMARKS:
!  Here, I'm taking a quick and dirty approach and just dividing by 4.
!   When we do the 3D model someone (not me) needs to write a decent
!   regridding program. nes, 30 april 09.
!
!EOP
!------------------------------------------------------------------------------
!BOC   
   INTEGER   :: II, JJ, NII, NJJ, LAT, NEWI, INDEX
  
   !initialize outgrid
   outgrid=0d0
   index=0
   DO II=1, GLAT
   
    NEWI=FLOOR(II/4d0)+1
   
   DO JJ=1,GLON
      OUTGRID(1,NEWI)=OUTGRID(1,NEWI)+INGRID(JJ,II)

   ENDDO
   ENDDO

      END FUNCTION GRID_2D
!EOC


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------ 
!          MIT INTEGRATED GLOBAL SYSTEMS MODEL                                 !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  WRITE_OUTPUT
!
! !DESCRIPTION: This subroutine writes totals for comparison
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE WRITE_OUTPUT(YEAR)
!
! !INPUT PARAMETERS: 
      INTEGER, INTENT(IN)::YEAR
!
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  7Oct09: Noelle Selin
!
! !REMARKS:
!  
!
!EOP
!------------------------------------------------------------------------------
!BOC   
    REAL :: totco2g, cfc11g, cfc12g, totn2og,totcog, totnoxg
    REAL :: totch4g, totso2g, tothfcg, totpfcg, totsf6g
    REAL :: totbc, totnh3, totoc
    LOGICAL, SAVE         :: FIRST_out = .TRUE. 

       IF(FIRST_OUT) Then
!open(17,file='edailytotals',status='new')
	open(17,file='edailytotals')
   !	write(17,*)'Annual Emission'
	write(17,900)
900	format(' Year  CO2:Pg  F11:Gg  F12:Gg  N2O:Tg  CO:Pg &
                  & NO:Tg  CH4:Pg  SO2:Tg  HFC:Tg  PFC:Gg  SF6:Gg &
		                          &  BC:Tg  NH3:Tg   OC:Tg')
        first_out=.false.
       END IF
  
       totco2g=sum(emistempco2)
       cfc11g=sum(emistempc11)
       cfc12g=sum(emistempc12)
       totn2og=sum(emistempn2o)
       totcog=sum(emistempco)
       totnoxg=sum(emistempno2)
       totch4g=sum(emistempch4)
       totso2g=sum(emistempso2)
       tothfcg=sum(emistemphfc)
       totpfcg=sum(emistemppfc)
       totsf6g=sum(emistempsf6)
       totbc=sum(emistempbc)
       totnh3=sum(emistempnh3)
       totoc=sum(emistempoc)


 

	  write(17,901)year,					&
                       totco2g*1.e-6,cfc11g,cfc12g,		&
                       totn2og*1.e-3,totcog*1.e-6,		&
                       totnoxg*1.e-3,totch4g*1.e-6,		&
                       totso2g*1.e-3,tothfcg*1.e-3,		&
                       totpfcg,totsf6g,				&
		       totbc*1.e-3,totnh3*1.e-3,totoc*1.e-3

901	format(i5,3f8.2,2f8.3,f8.2,f8.3,f8.2,3f8.3,3f8.2)


      END subroutine write_output
!EOC

!------------------------------------------------------------------------------
        !End of module
	END MODULE EPPANEW_MOD



!EOC
