C $Header: /home/igsm-repository/INC/PRNT.COM,v 1.1 2010/08/31 19:27:01 jscott Exp $
C $Name:  $

      print *,'J=',JPR,' L=',LPR,' P PT PX'
      print *,P(1,JPR),PT(1,JPR),PX(1,JPR)
      print *,' T(J,L)=',T(1,JPR,LPR),' Q(J,L)=',Q(1,JPR,LPR)
      print *,' TT(J,L)=',TT(1,JPR,LPR),' QT(J,L)=',QT(1,JPR,LPR)
      print *,' TX(J,L)=',TX(1,JPR,LPR),' DOPK=',DOPK
      print *,' V(J,L)=',V(1,JPR,LPR),' V(J+1,L)=',V(1,JPR+1,LPR)
      print *,' VT(J,L)=',VT(1,JPR,LPR),' VT(J+1,L)=',VT(1,JPR+1,LPR)
      print *,' VX(J,L)=',VX(1,JPR,LPR),' VX(J+1,L)=',VX(1,JPR+1,LPR)
      print *,' U(J,L)=',U(1,JPR,LPR),' U(J+1,L)=',U(1,JPR+1,LPR)
      print *,' UT(J,L)=',UT(1,JPR,LPR),' UT(J+1,L)=',UT(1,JPR+1,LPR)
      print *,' UX(J,L)=',UX(1,JPR,LPR),' UX(J+1,L)=',UX(1,JPR+1,LPR)
      print *,'ODATA(1,JPR,1)',ODATA(1,JPR,1)
c       print*,' AJL(59) from main'
c       do l=1,lm
c       do j=1,jm
c       if(j.eq.20.or.j.eq.2)then
c       print '(2i5,f10.2)',(l,j,AJL(j,l,59))
c       endif
c       enddo
c       enddo
