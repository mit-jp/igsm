!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: snowdp2lev
!
! !INTERFACE:
subroutine snowdp2lev ()
!
! !DESCRIPTION: 
! Create snow layers and interfaces given snow depth.
! Note that cps%zi(0) is set in routine iniTimeConst.
!
! !USES
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  use clm_varpar, only : nlevsno, nlevlak
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES
  integer j,gi,li,ci                           !indices
  type(gridcell_type)       , pointer :: g     !local pointer to derived subtype
  type(landunit_type)       , pointer :: l     !local pointer to derived subtype
  type(column_type)         , pointer :: c     !local pointer to derived subtype
  type(model_pstate_type)   , pointer :: mps   !local pointer to derived subtype
  type(gridcell_pstate_type), pointer :: gps   !local pointer to derived subtype
  type(landunit_pstate_type), pointer :: lps   !local pointer to derived subtype
  type(column_pstate_type)  , pointer :: cps   !local pointer to derived subtype
!-----------------------------------------------------------------------

  ! Assign local pointer for simpler referencing
  mps => clm%mps

  ! Grid level initialization
  do gi = 1,mps%ngridcells
     g => clm%g(gi)					
     gps => g%gps
     do li = 1,gps%nlandunits
        l => g%l(li)
        lps => l%lps
        do ci = 1,lps%ncolumns
           c => l%c(ci)
           cps => c%cps
           cps%dz(-nlevsno+1:0) = 1.e36 
           cps%z (-nlevsno+1:0) = 1.e36 
           cps%zi(-nlevsno:-1)  = 1.e36 
           if (.not. lps%lakpoi) then  !not lake
              if (cps%snowdp < 0.01) then
                 cps%snl = 0
                 cps%dz(-nlevsno+1:0) = 0.
                 cps%z (-nlevsno+1:0) = 0.
                 cps%zi(-nlevsno+0:0) = 0.
              else
                 if ((cps%snowdp >= 0.01) .AND. (cps%snowdp <= 0.03)) then
                    cps%snl = -1
                    cps%dz(0)  = cps%snowdp
                 else if ((cps%snowdp > 0.03) .AND. (cps%snowdp <= 0.04)) then
                    cps%snl = -2
                    cps%dz(-1) = cps%snowdp/2.
                    cps%dz( 0) = cps%dz(-1)
                 else if ((cps%snowdp > 0.04) .AND. (cps%snowdp <= 0.07)) then
                    cps%snl = -2
                    cps%dz(-1) = 0.02
                    cps%dz( 0) = cps%snowdp - cps%dz(-1)
                 else if ((cps%snowdp > 0.07) .AND. (cps%snowdp <= 0.12)) then
                    cps%snl = -3
                    cps%dz(-2) = 0.02
                    cps%dz(-1) = (cps%snowdp - 0.02)/2.
                    cps%dz( 0) = cps%dz(-1)
                 else if ((cps%snowdp > 0.12) .AND. (cps%snowdp <= 0.18)) then
                    cps%snl = -3
                    cps%dz(-2) = 0.02
                    cps%dz(-1) = 0.05
                    cps%dz( 0) = cps%snowdp - cps%dz(-2) - cps%dz(-1)
                 else if ((cps%snowdp > 0.18) .AND. (cps%snowdp <= 0.29)) then
                    cps%snl = -4
                    cps%dz(-3) = 0.02
                    cps%dz(-2) = 0.05
                    cps%dz(-1) = (cps%snowdp - cps%dz(-3) - cps%dz(-2))/2.
                    cps%dz( 0) = cps%dz(-1)
                 else if ((cps%snowdp > 0.29) .AND. (cps%snowdp <= 0.41)) then
                    cps%snl = -4
                    cps%dz(-3) = 0.02
                    cps%dz(-2) = 0.05
                    cps%dz(-1) = 0.11
                    cps%dz( 0) = cps%snowdp - cps%dz(-3) - cps%dz(-2) - cps%dz(-1)
                 else if ((cps%snowdp > 0.41) .AND. (cps%snowdp <= 0.64)) then
                    cps%snl = -5
                    cps%dz(-4) = 0.02
                    cps%dz(-3) = 0.05
                    cps%dz(-2) = 0.11
                    cps%dz(-1) = (cps%snowdp - cps%dz(-4) - cps%dz(-3) - cps%dz(-2))/2.
                    cps%dz( 0) = cps%dz(-1)
                 else if (cps%snowdp > 0.64) then 
                    cps%snl = -5
                    cps%dz(-4) = 0.02
                    cps%dz(-3) = 0.05
                    cps%dz(-2) = 0.11
                    cps%dz(-1) = 0.23
                    cps%dz( 0)=cps%snowdp-cps%dz(-4)-cps%dz(-3)-cps%dz(-2)-cps%dz(-1)
                 endif
                 do j = 0, cps%snl+1, -1
                    cps%z(j)    = cps%zi(j) - 0.5*cps%dz(j)
                    cps%zi(j-1) = cps%zi(j) - cps%dz(j)
                 enddo
              endif
           else   !lake points
              cps%snl = 0
              cps%dz(-nlevsno+1:0) = 0.
              cps%z (-nlevsno+1:0) = 0.
              cps%zi(-nlevsno+0:0) = 0.
           end if  ! end if-lake  block

        end do ! end column level loop
     end do  ! end landunit level loop
  end do  ! end grid level loop

end subroutine snowdp2lev
