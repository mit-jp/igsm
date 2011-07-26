#include <misc.h>
#include <preproc.h>

module pftvarcon

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: pftvarcon
! 
! !DESCRIPTION: 
! Module containing vegetation constants and method to 
! eads and initialize vegetation (PFT) constants.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varpar
!
! !PUBLIC TYPES:
  implicit none
  save
!
! Vegetation type constants
!
  character(len=40) pftname(0:numpft) !PFT description

  integer ncorn                  !value for corn
  integer nwheat                 !value for wheat
  integer noveg                  !value for not vegetated 
  integer ntree                  !value for last type of tree

  real(r8) dleaf(0:numpft)       !characteristic leaf dimension (m) 
  real(r8) c3psn(0:numpft)       !photosynthetic pathway: 0. = c4, 1. = c3
  real(r8) vcmx25(0:numpft)      !max rate of carboxylation at 25C (umol CO2/m**2/s)
  real(r8) mp(0:numpft)          !slope of conductance-to-photosynthesis relationship
  real(r8) qe25(0:numpft)        !quantum efficiency at 25C (umol CO2 / umol photon)
  real(r8) xl(0:numpft)          !leaf/stem orientation index
  real(r8) rhol(0:numpft,numrad) !leaf reflectance: 1=vis, 2=nir 
  real(r8) rhos(0:numpft,numrad) !stem reflectance: 1=vis, 2=nir 
  real(r8) taul(0:numpft,numrad) !leaf transmittance: 1=vis, 2=nir 
  real(r8) taus(0:numpft,numrad) !stem transmittance: 1=vis, 2=nir 
  real(r8) z0mr(0:numpft)        !ratio of momentum roughness length to canopy top height (-)
  real(r8) displar(0:numpft)     !ratio of displacement height to canopy top height (-)
  real(r8) roota_par(0:numpft)   !CLM rooting distribution parameter [1/m]
  real(r8) rootb_par(0:numpft)   !CLM rooting distribution parameter [1/m]

  real(r8) sla(0:numpft)              !sp. leaf area [m2 leaf g-1 carbon]
  real(r8) pftpar(0:numpft,1:npftpar) !the rest for use with DGVM
  real(r8) lm_sapl(0:numpft)
  real(r8) sm_sapl(0:numpft)
  real(r8) hm_sapl(0:numpft)
  real(r8) rm_sapl(0:numpft)
  logical  tree(0:numpft)
  logical  summergreen(0:numpft)
  logical  raingreen(0:numpft)

  real(r8), parameter :: reinickerp = 1.6 !parameter in allometric equation
  real(r8), parameter :: wooddens = 2.0e5 !wood density (gC/m3)
  real(r8), parameter :: latosa = 8.0e3   !ratio of leaf area to sapwood cross-sectional 
                                          !area (Shinozaki et al 1964a,b)
  real(r8), parameter :: allom1 = 100.0   !parameters in allometric
  real(r8), parameter :: allom2 =  40.0
  real(r8), parameter :: allom3 =   0.5
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: pftconrd ! Read and initialize vegetation (PFT) constants 
!
! !REVISION HISTORY:
! Created by Sam Levis (put into module form by Mariana Vertenstein)
!
!EOP
!----------------------------------------------------------------------- 

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: pftconrd
!
! !INTERFACE:
  subroutine pftconrd
!
! !DESCRIPTION: 
! Read and initialize vegetation (PFT) constants 
!
! !USES:
    use fileutils, only : opnfil, getfil, relavu, getavu
    use clm_varctl, only : fpftcon	
#if (defined SPMD)
    use spmdMod, only : masterproc, mpicom, MPI_REAL8 
#else
    use spmdMod, only : masterproc        
#endif
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! routine initialize in module initializeMod
!
! !REVISION HISTORY:
! Created by Gordon Bonan
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=256) :: locfn !local file name
    integer :: i,n              !loop indices
    integer :: ier              !error code
!-----------------------------------------------------------------------

    ! Set specific vegetation type values

    ncorn  = 15
    nwheat = 16

    ! Set value for last type of tree

    ntree = 8  !value for last type of tree

    ! Set value for non-vegetated

    noveg = 0  !value

    ! Assign unit number to file. Get local file. 
    ! Open file and read PFT's.
    ! Close and release file.

    if (masterproc) then
       write (6,*) 'Attempting to read PFT physiological data .....'
       n = getavu()
       call getfil (fpftcon, locfn, 0)
       call opnfil (locfn, n, 'f')
       do i = 1, numpft
          read (n,*)  pftname(i),              &
                      z0mr(i)   , displar(i), dleaf(i)  , c3psn(i)  , &
                      vcmx25(i) , mp(i)     , qe25(i)   , rhol(i,1) , &
                      rhol(i,2) , rhos(i,1) , rhos(i,2) , taul(i,1) , &
                      taul(i,2) , taus(i,1) , taus(i,2) , xl(i)     , &
                      roota_par(i), rootb_par(i) 
       end do
       call relavu (n)
    endif

    ! Define PFT zero to be bare ground

    pftname(noveg) = 'not_vegetated'
    z0mr(noveg) = 0.
    displar(noveg) = 0.
    dleaf(noveg) = 0.
    c3psn(noveg) = 1.
    vcmx25(noveg) = 0.
    mp(noveg) = 9.
    qe25(noveg) = 0.
    rhol(noveg,1) = 0.
    rhol(noveg,2) = 0.
    rhos(noveg,1) = 0.
    rhos(noveg,2) = 0.
    taul(noveg,1) = 0.
    taul(noveg,2) = 0.
    taus(noveg,1) = 0.
    taus(noveg,2) = 0.
    xl(noveg) = 0.
    roota_par(noveg) = 0.
    rootb_par(noveg) = 0.

#if ( defined SPMD )
    ! pass surface data to all processors

    call mpi_bcast (z0mr, size(z0mr), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (displar, size(displar), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (dleaf, size(dleaf), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (c3psn, size(c3psn), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (vcmx25, size(vcmx25), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (mp, size(mp), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (qe25, size(qe25), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (rhol, size(rhol), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (rhos, size(rhos), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (taul, size(taul), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (taus, size(taus), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (xl, size(xl), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (roota_par, size(roota_par), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (rootb_par, size(rootb_par), MPI_REAL8, 0, mpicom, ier)
#endif

    if (masterproc) then
       write (6,*) 'Successfully read PFT physiological data'
       write (6,*)
    endif

    return 
  end subroutine pftconrd

end module pftvarcon

