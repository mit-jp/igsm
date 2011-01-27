#include <misc.h>
#include <preproc.h>

module mapxy

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: mapxy
! 
! !DESCRIPTION: 
! Utilities to perfrom grid-average from subgrid 1d vector and 
! grid average to 1d vector
!
! !PUBLIC TYPES
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS
  public :: vec2xy ! Perfrom grid-average from subgrid 1d vector   
  public :: xy2vec ! Convert a grid-average field to subgrid vector
!
! !REVISION HISTORY
! Created by Mariana Vertenstein
!
!EOP
!----------------------------------------------------------------------- 

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: vec2xy
!
! !INTERFACE:
  subroutine vec2xy (clmlevel, fldv, fldxyini, fldxy)
!
! !DESCRIPTION: 
! Perfrom grid-average from subgrid 1d vector 
! Subgrid to grid average mapping: average a subgrid input vector 
! [fldv] of length to a 2-d [lsmlon] x [lsmlat] output array [fldxy]
! setting non-land points to [nonland].  Averaging is only done for 
! points that are not equal to "spval".
!
! !USES
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype, only : grid1d, land1d, cols1d, pfts1d
    use clm_varsur, only : landmask
    use clm_varpar, only : lsmlon, lsmlat
    use clm_varcon, only : spval
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: clmlevel      !type of clm 1d array
    real(r8), pointer, dimension(:) :: fldv       !1d vector input
    real(r8), intent(in)  :: fldxyini             !initial value of fldxy
    real(r8), intent(out) :: fldxy(lsmlon,lsmlat) !xy gridded output
!
! !REVISION HISTORY
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES
    integer :: i,j,gi,li,ci,pi      !indices
    real(r8):: sumwt(lsmlon,lsmlat) !sum of wt
!-----------------------------------------------------------------------

    fldxy(:,:) = fldxyini
    sumwt(:,:) = 0.
    
    select case (clmlevel)
    case('gridcell')
       do gi = 1, grid1d%num
          if (fldv(gi) /= spval) then
             i = grid1d%ixy(gi) 
             j = grid1d%jxy(gi) 
             if (sumwt(i,j)==0.) fldxy(i,j) = 0.
             fldxy(i,j) = fldxy(i,j) + grid1d%wtxy(gi)*fldv(gi)
             sumwt(i,j) = sumwt(i,j) + grid1d%wtxy(gi)
          endif
       end do
    case('landunit')
       do li = 1, land1d%num
          if (fldv(li) /= spval) then
             i = land1d%ixy(li) 
             j = land1d%jxy(li) 
             if (sumwt(i,j)==0.) fldxy(i,j) = 0.
             fldxy(i,j) = fldxy(i,j) + land1d%wtxy(li)*fldv(li)
             sumwt(i,j) = sumwt(i,j) + land1d%wtxy(li)
          endif
       end do
    case('column')
       do ci = 1, cols1d%num
          if (fldv(ci) /= spval) then
             i = cols1d%ixy(ci) 
             j = cols1d%jxy(ci) 
             if (sumwt(i,j)==0.) fldxy(i,j) = 0.
             fldxy(i,j) = fldxy(i,j) + cols1d%wtxy(ci)*fldv(ci)
             sumwt(i,j) = sumwt(i,j) + cols1d%wtxy(ci)
          endif
       end do
    case('pft')
       do pi = 1, pfts1d%num
          if (fldv(pi) /= spval) then
             i = pfts1d%ixy(pi)      !longitude index for land point
             j = pfts1d%jxy(pi)      !latitude index for land point
             if (sumwt(i,j)==0.) fldxy(i,j) = 0.
             fldxy(i,j) = fldxy(i,j) + pfts1d%wtxy(pi)*fldv(pi)
             sumwt(i,j) = sumwt(i,j) + pfts1d%wtxy(pi)
          endif
       end do
    case default
       write(6,*) 'VEC2XY: Invalid expansion character: ',trim(clmlevel)
       call endrun
    end select
    
    where (landmask(:,:) == 1 .and. sumwt(:,:) /= 0.)
       fldxy(:,:) = fldxy(:,:)/sumwt(:,:)
    endwhere
    
    return
  end subroutine vec2xy
  
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: xy2vec
!
! !INTERFACE:
  subroutine xy2vec (clmlevel, fldxy, fldv)
!
! !DESCRIPTION: 
! Convert a grid-average field to a subgrid vector
! This code converts a grid-average field [fldxy] dimensioned
! [lsmlon] x [lsmlat] to a 1d gridcell, landunit, column or 
! pft vector [fldv]
!
! !USES
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype, only : grid1d, land1d, cols1d, pfts1d
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: clmlevel !type of clm 1d array
    real(r8), intent(in) :: fldxy(:,:)       !gridded input
    real(r8), pointer, dimension(:) :: fldv  !subgrid vector output
!
! !REVISION HISTORY
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES
    integer :: i,j,gi,li,ci,pi   !indices
    integer :: ni, nf            !beginning and ending 1d indices
!-----------------------------------------------------------------------

    ni = lbound(fldv,dim=1)
    nf = ubound(fldv,dim=1)
    
    select case (clmlevel)
    case('gridcell')
       do gi = ni,nf
          i = grid1d%ixy(gi)
          j = grid1d%jxy(gi)
          fldv(gi) = fldxy(i,j)
       end do
    case('landunit')
       do li = ni,nf
          i = land1d%ixy(li)
          j = land1d%jxy(li)
          fldv(li) = fldxy(i,j)
       end do
    case('column')
       do ci = ni,nf
          i = cols1d%ixy(ci)
          j = cols1d%jxy(ci)
          fldv(ci) = fldxy(i,j)
       end do
    case('pft')
       do pi = ni,nf
          i = pfts1d%ixy(pi)
          j = pfts1d%jxy(pi)
          fldv(pi) = fldxy(i,j)
       end do
    case default
       write(6,*) 'XY2VEC: Invalid expansion character: ',trim(clmlevel)
       call endrun
    end select
    
    return
  end subroutine xy2vec
  
!=======================================================================

end module mapxy

