subroutine readIGSM2()
include 'IGSM2.inc'
    integer,parameter:: fldsize=14
	real :: latincr
    character(len=8) :: fldlst(fldsize)      !name of possible atm fields in input file
    data fldlst( 1) /'TBOT    '/
    data fldlst( 2) /'WIND    '/
    data fldlst( 3) /'QBOT    '/
    data fldlst( 4) /'Tdew    '/
    data fldlst( 5) /'RH      '/
    data fldlst( 6) /'ZBOT    '/
    data fldlst( 7) /'PSRF    '/
    data fldlst( 8) /'FSDS    '/
    data fldlst( 9) /'FSDSdir '/
    data fldlst(10) /'FSDSdif '/
    data fldlst(11) /'FLDS    '/
    data fldlst(12) /'PRECTmms'/
    data fldlst(13) /'PRECCmms'/
    data fldlst(14) /'PRECLmms'/

     read(20) psmit,pcplmit,pcpcmit,tprmit,tslmit,qsmit, &
          wsmit,usmit,vsmit,dswmit,dlwmit,pco2mit,swnirmit,swparmit
!
!	First - fill in (laboriously) empty latitude bands with data
!
     psmit(1,8)=psmit(1,10)
     pcplmit(1,8)=pcplmit(1,10)
     pcpcmit(1,8)=pcpcmit(1,10)
     tprmit(1,8)=tprmit(1,10)
     tslmit(1,8)=tslmit(1,10)
     qsmit(1,8)=qsmit(1,10)
     wsmit(1,8)=wsmit(1,10)
     usmit(1,8)=usmit(1,10)
     vsmit(1,8)=vsmit(1,10)
     dswmit(1,8)=dswmit(1,10)
     dlwmit(1,8)=dlwmit(1,10)
     pco2mit(1,8)=pco2mit(1,10)
     swnirmit(1,8)=swnirmit(1,10)
     swparmit(1,8)=swparmit(1,10)
     psmit(1,9)=psmit(1,10)
     pcplmit(1,9)=pcplmit(1,10)
     pcpcmit(1,9)=pcpcmit(1,10)
     tprmit(1,9)=tprmit(1,10)
     tslmit(1,9)=tslmit(1,10)
     qsmit(1,9)=qsmit(1,10)
     wsmit(1,9)=wsmit(1,10)
     usmit(1,9)=usmit(1,10)
     vsmit(1,9)=vsmit(1,10)
     dswmit(1,9)=dswmit(1,10)
     dlwmit(1,9)=dlwmit(1,10)
     pco2mit(1,9)=pco2mit(1,10)
     swnirmit(1,9)=swnirmit(1,10)
     swparmit(1,9)=swparmit(1,10)
     psmit(1,45)=psmit(1,44)
     pcplmit(1,45)=pcplmit(1,44)
     pcpcmit(1,45)=pcpcmit(1,44)
     tprmit(1,45)=tprmit(1,44)
     tslmit(1,45)=tslmit(1,44)
     qsmit(1,45)=qsmit(1,44)
     wsmit(1,45)=wsmit(1,44)
     usmit(1,45)=usmit(1,44)
     vsmit(1,45)=vsmit(1,44)
     dswmit(1,45)=dswmit(1,44)
     dlwmit(1,45)=dlwmit(1,44)
     pco2mit(1,45)=pco2mit(1,44)
     swnirmit(1,45)=swnirmit(1,44)
     swparmit(1,45)=swparmit(1,44)
     psmit(1,46)=psmit(1,44)
     pcplmit(1,46)=pcplmit(1,44)
     pcpcmit(1,46)=pcpcmit(1,44)
     tprmit(1,46)=tprmit(1,44)
     tslmit(1,46)=tslmit(1,44)
     qsmit(1,46)=qsmit(1,44)
     wsmit(1,46)=wsmit(1,44)
     usmit(1,46)=usmit(1,44)
     vsmit(1,46)=vsmit(1,44)
     dswmit(1,46)=dswmit(1,44)
     dlwmit(1,46)=dlwmit(1,44)
     pco2mit(1,46)=pco2mit(1,44)
     swnirmit(1,46)=swnirmit(1,44)
     swparmit(1,46)=swparmit(1,44)
	 latincr = 2.
	 do j = 1, datlat
	   j1 = int(j/latincr)+1
           x(:,j,1)=tslmit(1,j1)
	   x(:,j,2)=sqrt((usmit(1,j1)*usmit(1,j1))+(vsmit(1,j1)*vsmit(1,j1)))
	   x(:,j,3)=qsmit(1,j1)
	   x(:,j,7)=psmit(1,j1)*100
	   x(:,j,6) = 100.0
	   x(:,j,9) = 0.7*swparmit(1,j1)+0.7*swnirmit(1,j1)
	   x(:,j,10) = 0.3*swparmit(1,j1)+0.3*swnirmit(1,j1)
	   x(:,j,11) = dlwmit(1,j1)
	   x(:,j,13) = pcpcmit(1,j1)/3600.
	   x(:,j,14) = pcplmit(1,j1)/3600.
	 enddo
     !write (6,*) 'PROGRAM_OFF:'
     !write (6,*) 'Reading IGSM2 forcing data for step: ',nstep
end subroutine readIGSM2
