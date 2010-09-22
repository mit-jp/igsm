C $Header: /home/igsm-repository/ATM/fltr4tem.F,v 1.1 2010/08/30 15:36:35 jscott Exp $
C $Name:  $


        subroutine fltr4tem(z,iy)
        dimension z(iy),zt(iy)
          do j=2,iy-2
             zt(j)=0.25*(z(j-1)+2.*z(j)+z(j+1))
          enddo
          do j=2,iy-1
           z(j)=zt(j)
          enddo
        return
        end
