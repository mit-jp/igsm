!       ============================================================
!
!       TRACER_COM:    An include file consists of COMMONs describing
!                       mixing ratios of tracers
!                       for MIT IGSM
      parameter (ntracers=5)
!     nlon, nlat and nlev are setup in chem_para
!
      real tracers(nlon,nlat,nlev,ntracers),
     &  treftime(ntracers),tottracers(ntracers)
      integer  trtype(ntracers)
      common /trac1/ tracers,tottracers,treftime,trtype
!
      real tracmass(ntracers),tracmass1(ntracers),
     &   tracersms(ntracers)
      common /trac2/ tracmass,tracmass1,tracersms
!
      real tracmonth(nlon,nlat,nlev,ntracers),
     &  tracglob(ntracers)
      common /trac3/ tracmonth,tracglob,navtrac,monthntrac
!
      real tracemis(nlon,nlat,nlev_accri,ntracers)
     &  ,trac_emissions(nlat,nlev_accri,365)
      integer  IYREMIS
      common /trac4/ tracemis,trac_emissions,IYREMIS
!
      real tracxc(nlev,ntracers),tracsum(ntracers),
     &  tracold(ntracers), traccld(ntracers),tracup(ntracers)
      common /trac5/ tracxc,tracsum,tracold,traccld,tracup
!
!
!     data trtype /1,2,0,3,4/
!      types of tracers
!      0 - inert tracer (no losses)
!      1 - NOx like tracer with prescibed e-folding time
!      2 -  tropospheric tracer (Ozone)with prescibed e-folding time
!      3 - Black carbon like tracer with dry removal
!      4 - Black carbon like tracer with dry and wet removal
!     data treftime /5.0,21.0,0.0,0.0,0.0/
