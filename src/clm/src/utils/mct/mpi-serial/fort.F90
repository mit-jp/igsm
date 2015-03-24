

	subroutine clmmpi_init(ierror)

        implicit none
        include "mpif.h"

	integer fint(2)
	logical flog(2)
	real freal(2)
	double precision fdub(2)
	complex fcomp(2)
	integer status(CLMMPI_STATUS_SIZE)

        integer ierror


        !!
        !! Pass values from clmmpif.h to the C side
        !! to check for consistency clmmpi.h and hardware sizes.
        !!

        call clmmpi_init_fort( CLMMPI_COMM_WORLD, &
	                    CLMMPI_ANY_SOURCE, CLMMPI_ANY_TAG, &
                            CLMMPI_COMM_NULL, CLMMPI_REQUEST_NULL, &
                            CLMMPI_GROUP_NULL, CLMMPI_GROUP_EMPTY, &
                            CLMMPI_UNDEFINED, &
                            CLMMPI_MAX_ERROR_STRING, &
                            CLMMPI_MAX_PROCESSOR_NAME, &
                            CLMMPI_STATUS_SIZE, &
                            CLMMPI_SOURCE, CLMMPI_TAG, CLMMPI_ERROR, &
                            status, status(CLMMPI_SOURCE), &
                            status(CLMMPI_TAG), status(CLMMPI_ERROR), &
                            CLMMPI_INTEGER, fint(1), fint(2), &
                            CLMMPI_LOGICAL, flog(1), flog(2), &
                            CLMMPI_REAL, freal(1), freal(2), &
                            CLMMPI_DOUBLE_PRECISION, fdub(1), fdub(2), &
                            CLMMPI_COMPLEX, fcomp(1), fcomp(2), &
                            IERROR )


        return
	end





