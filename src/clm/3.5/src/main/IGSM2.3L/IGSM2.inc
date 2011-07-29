!
! atmospheric gridded data
!   
    integer,parameter :: igsmlon=1
    integer,parameter :: igsmlat=46
    double precision ::  psmit (igsmlon,igsmlat)
    double precision ::  pcptmit (igsmlon,igsmlat)
    double precision ::  pcplmit (igsmlon,igsmlat)
    double precision ::  pcpcmit (igsmlon,igsmlat)
    double precision ::  tprmit (igsmlon,igsmlat)
    double precision ::  tslmit (igsmlon,igsmlat)
    double precision ::  qsmit (igsmlon,igsmlat)
    double precision ::  wsmit (igsmlon,igsmlat)
    double precision ::  usmit (igsmlon,igsmlat)
    double precision ::  vsmit (igsmlon,igsmlat)
    double precision ::  dswmit (igsmlon,igsmlat)
    double precision ::  dlwmit (igsmlon,igsmlat)
    double precision ::  pco2mit (igsmlon,igsmlat)
!    double precision ::  po3mit (igsmlon,igsmlat)
    double precision ::  swnirmit (igsmlon,igsmlat)
    double precision ::  swparmit (igsmlon,igsmlat)
if (defined PCPTRND)
    real(r8) :: tglbtrnd
    common/mit2din/psmit,pcplmit,pcpcmit,tprmit,tslmit,qsmit,wsmit,usmit,vsmit,dswmit,dlwmit,pco2mit,swnirmit,swparmit,tglbtrnd
#else
    common/mit2din/psmit,pcplmit,pcpcmit,tprmit,tslmit,qsmit,wsmit,usmit,vsmit,dswmit,dlwmit,pco2mit,swnirmit,swparmit
#endif
