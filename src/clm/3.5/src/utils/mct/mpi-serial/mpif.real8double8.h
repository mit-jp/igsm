
!!!
!!! clmmpif.real8double8.h (real is 8 bytes, double precision 8)
!!!
!!!   NOTE: clmmpif.h is copied from either clmmpif.XXXX.h when the
!!!   library is built.  Do not modify clmmpif.h, the changes will 
!!!   get clobbered.
!!!


!
! CLMMPI_COMM_WORLD
!

	INTEGER CLMMPI_COMM_WORLD
        parameter (clmmpi_comm_world=1)


!
! source,tag
!

	integer CLMMPI_ANY_SOURCE, CLMMPI_ANY_TAG
        parameter (clmmpi_any_source=-1, clmmpi_any_tag= -1)


        integer CLMMPI_COMM_NULL, CLMMPI_REQUEST_NULL
        parameter (CLMMPI_COMM_NULL=0, CLMMPI_REQUEST_NULL=0)

        integer CLMMPI_GROUP_NULL, CLMMPI_GROUP_EMPTY
        parameter (CLMMPI_GROUP_NULL=0, CLMMPI_GROUP_EMPTY= -1)

        integer CLMMPI_MAX_ERROR_STRING
        parameter (CLMMPI_MAX_ERROR_STRING=128)

        integer CLMMPI_MAX_PROCESSOR_NAME
        parameter (CLMMPI_MAX_PROCESSOR_NAME=128)


        integer CLMMPI_SUCCESS
        parameter (CLMMPI_SUCCESS=0)

        integer CLMMPI_UNDEFINED
        parameter (CLMMPI_UNDEFINED= -1)


!
! CLMMPI_Status
!
! The values in this section MUST match the struct definition
! in clmmpi.h
!


        INTEGER CLMMPI_STATUS_SIZE
        PARAMETER (CLMMPI_STATUS_SIZE=3)

        INTEGER CLMMPI_SOURCE, CLMMPI_TAG, CLMMPI_ERROR
        PARAMETER(CLMMPI_SOURCE=1, CLMMPI_TAG=2, CLMMPI_ERROR=3)



!
! CLMMPI_Datatype values
!  
! The value is the size of the datatype in bytes.
! Change if necessary for the machine in question.
! (The clmmpi.h file uses sizeof(), so it should be more
! portable).
! 
!


	INTEGER CLMMPI_BYTE
	PARAMETER (CLMMPI_BYTE=1)

	INTEGER CLMMPI_CHARACTER
	PARAMETER (CLMMPI_CHARACTER=1)

	INTEGER CLMMPI_REAL4
	PARAMETER (CLMMPI_REAL4=4)

	INTEGER CLMMPI_REAL8
	PARAMETER (CLMMPI_REAL8=8)

	INTEGER CLMMPI_INTEGER
	PARAMETER (CLMMPI_INTEGER=4)

	INTEGER CLMMPI_LOGICAL
	PARAMETER (CLMMPI_LOGICAL=4)

!!!!!!!
	INTEGER CLMMPI_REAL
	PARAMETER (CLMMPI_REAL=8)

	INTEGER CLMMPI_DOUBLE_PRECISION
	PARAMETER (CLMMPI_DOUBLE_PRECISION=8)
!!!!!!!

	integer CLMMPI_COMPLEX
	parameter (CLMMPI_COMPLEX=2*CLMMPI_REAL)

        integer CLMMPI_2REAL
        parameter (CLMMPI_2REAL=2*CLMMPI_REAL)

        integer CLMMPI_2DOUBLE_PRECISION
        parameter (CLMMPI_2DOUBLE_PRECISION=2*CLMMPI_DOUBLE_PRECISION)

        integer CLMMPI_2INTEGER
        parameter (CLMMPI_2INTEGER=2*CLMMPI_INTEGER)

!
! CLMMPI_Op values
!
! (All are handled as no-op so no value is necessary)
!

        INTEGER CLMMPI_SUM, CLMMPI_MAX, CLMMPI_MIN, CLMMPI_PROD, CLMMPI_LAND, CLMMPI_BAND
        INTEGER CLMMPI_LOR, CLMMPI_BOR, CLMMPI_LXOR, CLMMPI_BXOR, CLMMPI_MINLOC
        INTEGER CLMMPI_MAXLOC
        INTEGER CLMMPI_OP_NULL

!
! CLMMPI_Wtime
!

        DOUBLE PRECISION CLMMPI_WTIME
        EXTERNAL CLMMPI_WTIME
