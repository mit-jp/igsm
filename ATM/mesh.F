
#include "ctrparam.h"

!	==========================================================
!
!	MESH.F:  subroutines for setting vertical coordinates 
!			of the model.
!
!	----------------------------------------------------------
!
!	Revision History:
!	
!	When	Who		What
!	----	----------	-------	
!	073100	Chien Wang	repack based on CliChem3 & M24x11
!
!	==========================================================

        SUBROUTINE MESH09(N,DST,DSB,SIGF,SIGH,DSIGH,DSIGF)

	IMPLICIT REAL*8 (A-H,O-Z)
	dimension sigf(n),sigh(0:n),dsigh(n),dsigf(n)
	dimension sig0(9),sige0(0:9)
creal *4  sig0(9),sige0(0:9)
CB for LM=9
      DATA SIG0/.974264e0,.907372e0,.796957e0,.640124e0,.470418e0,
     &  .318899e0,    
     *  .195759e0,.094938e0,.016897e0/                             
      DATA SIGE0/1.e0,.948665e0,.866530e0,.728953e0,.554415e0,.390144e0,   
     *  .251540e0,.143737e0,.061602e0,0.e0/                           
	ptop=10.e0
	psrf=984.e0
        do l=0,9
         sigh(n-l)=sige0(l)
        enddo
        do l=1,9
          sigf(n-l+1)=sig0(l)
        enddo
        print *,' '
        print '(A35,I2)','VERTICAL GRID FOR LM=',N
        print *,' '
        print *,'   L       SIGH       PH        SIGF        PF'
        print *,' '
	do l=0,n
	pp=sigh(l)*(psrf-ptop)+ptop
	if(l.eq.0)then
	 print '(i5,2f10.2)',l,sigh(l),pp
	else
	 pp1=sigf(l)*(psrf-ptop)+ptop
	 print '(i5,4f10.2)',l,sigh(l),pp,sigf(l),pp1
	endif
	enddo
	return 
	end

        SUBROUTINE MESH11(N,DST,DSB,SIGF,SIGH,DSIGH,DELSF)

	IMPLICIT REAL*8 (A-H,O-Z)
	dimension sigf(n),sigh(0:n),dsigh(n),dsigf(n)
	real sig0(11),sige0(0:11),prnew(5)
	data prnew/30.,60.,100.,150.,200./
CB for LM=11
      DATA SIGE0/1.,.948665,.866530,.728953,.554415,.390144,
     *  .251540,.143737,.092402,.0513347,.02053388,0./
      DATA SIG0/.974264,.907372,.796957,.640124,.470418,.318899,
     *  .195759,.112936,.066735,.0308008,.0051335/
CE for LM=11
	do l=0,11
	 sigh(n-l)=sige0(l)
	enddo
	do l=1,11
	  sigf(n-l+1)=sig0(l)
	enddo
	ptop=10.
	psrf=984.
        do l=1,4
          pp=sigh(l)*(psrf-ptop)+ptop
          ppm1=sigh(l-1)*(psrf-ptop)+ptop
          ppf=sqrt(pp*ppm1)
          sigf(l)=(ppf-ptop)/(psrf-ptop)
        enddo
        print *,' '
        print '(A35,I2)','VERTICAL GRID FOR LM=',N
        print *,' '
        print *,'   L       SIGH       PH        SIGF        PF'
        print *,' '
	do l=0,n
	pp=sigh(l)*(psrf-ptop)+ptop
	if(l.eq.0)then
	 print '(i5,f12.5,f10.2)',l,sigh(l),pp
	else
	 pp1=sigf(l)*(psrf-ptop)+ptop
	 print '(i5,2(f12.5,f10.2))',l,sigh(l),pp,sigf(l),pp1
	endif
	enddo
	return 
	end
