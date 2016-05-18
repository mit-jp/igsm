
#include "ctrparam.h"

      program redgary
#include "ODIFF.COM"
      parameter(nlat=N_LAT,nlatp1=nlat+1,nlon=N_LON0)
      character *3 reg,areg(3)
      real* 4 Rco2av0(nlat,lmo),pC_atmav0(nlat),pC_mlaav0(nlat),
c
     &Flux_atm_mlav0(nlat),Hav0(nlat),THETAHav0(nlat),HenryEav0(nlat),
     &DIC_mlav0(nlat),pC_mlav0(nlat)
      real *4 anexp
     & ,tsrfm(nlat,12),wsrfm(nlat,12),asrfm(nlat,12)
     & ,tsrf(nlat),wsrf(nlat),asrf(nlat)
      real * 4  co2rat1990(nlat,12),co2rat1860(nlat,12)
c
      common /Garyav/RCO2av(nlat,lmo),pC_atmav(nlat),pC_mlaav(nlat),
     &Flux_atm_mlav(nlat),Hav(nlat),
     &DIC_mlav(nlat),TCO2_mlav(nlat),AlkCav(nlat),
     &AlkBav(nlat),AlkSav(nlat),AlkWav(nlat),
     &CO2av(nlat),CO3av(nlat)
      real * 4 RCO2av,pC_atmav,pC_mlaav,
     &Flux_atm_mlav,Hav,
     &DIC_mlav,TCO2_mlav,AlkCav,
     &AlkBav,AlkSav,AlkWav,
     &CO2av,CO3av
c
      real  FDATA(N_LON0,N_LAT,3)
       real fricem(nlat,12)
     * ,dxyp(nlat),yp(nlat),cosp(nlat),cosv(nlat),yv(nlatp1)
     & ,dxv(nlat),dyv(nlat),fluxao(nlat),heE(nlat)
     & ,volume_ml(nlat)
      integer Ico2(nlat,lmo)
      common/Garyclim/t_srf(nlat),w_srf(nlat),area_ml(nlat)
     &  ,area_srf(nlat)
      common/Garydiff/depthml(nlat),ed(nlat),dzg(lmo),dzog(lmo-1),
     &Rco2(nlat,lmo),edohd(lmo),zg(lmo),focean(nlat)
      common /Garyflux/pC_atm(nlat),wind_amp
      real  Rco2dec(nlat,lmo),Hdec(nlat),kvc,kvc0,kvh0,kvcmin
      common /Garychem/Hg(nlat)
      namelist /CO2UP/nrun,areg,kvh0,kvcmin,rkckh,restart,lyeartr,
     &  co2eq,nyeareq,iyeareq,
     & ocean_eq,deepo_eq,ocean_tr,deepo_tr,ouptake
     &  ,dirdat2,oco2file,co2rfile,fl_init_dic
     &  ,fl_init_alkt,fl_init_salt
     &  ,fl_ed,fl_zoav
     &  ,eqclimate,trclimate
      character * 128 ocean_eq,deepo_eq,ocean_tr,deepo_tr,ouptake
     &  ,dirdat2,oco2file,co2rfile,fl_init_dic
     &  ,fl_init_alkt,fl_init_salt
     &  ,fl_ed,fl_zoav
     &  ,eqclimate,trclimate
      logical restart
C
      print *,'Start of run with off line OCM model'
      print *,'N_LAT=',N_LAT
      rkckh=0.42
      rkckh=1.375
      rkckh=0.275
      kvcmin=3.00
      open (615,file='ocup.name')
      read (615,CO2UP)
      id2=index(dirdat2," ")
c     print *,'From ocm dirdat2=',dirdat2
      open( unit=526,file=dirdat2(1:id2-1)//'FILE26',
     *    status='OLD',form='unformatted')
      rad=6375000.
      pi=4.*atan(1.)
      ncalldec=0
      dyp=pi/float(nlat-1)
      yp(1)=-pi/2.
      cosp(1)=0.
      do 65 j=2,nlat
        yp(j)=yp(j-1)+dyp
        cosp(j)=cos(yp(j))
        cosv(j)=0.5*(cosp(j)+cosp(j-1))
        dxv(j)=2.*pi*rad*cosv(j)
        dyv(j)=rad*(yp(j)-yp(j-1))
   65 continue
      cosv(nlat/2+1)=1.
      dxyp(1)=0.5*cosv(2)*0.5
      dxyp(nlat)=0.5*cosv(nlat)*0.5
      do 66 j=2,nlat-1
         dxyp(j)=0.5*(cosv(j)+cosv(j+1))
   66 continue
      do j=1,nlat
        dxyp(j)=dxyp(j)*dyp*2.*pi*rad**2
      enddo
      sn=0.
      ss=0.
      do 67 j=1,nlat/2
        sn=sn+dxyp(j)
        ss=ss+dxyp(j+12)
   67 continue
C *******
      atmmass=(ss+sn)*975.*100./9.81
      Chien=1./1000.*1.e6/12.*29/atmmass
      Chien=1.e-15/chien
C *******
      READ (526) FDATA
      do 284 J=1,nlat
        ILAND=0.
        CONT1=0.
        do 285 I=1,nlon
           PLAND=FDATA(I,J,2)
           CONT1=CONT1+PLAND
           ILAND=ILAND+1
  285   continue
        do 286 I=1,nlon
           IF(ILAND.GT.0)FDATA(I,J,2)=CONT1/ILAND
           IF(FDATA(I,J,2).LT.0.05)FDATA(I,J,2)=0.
  286   continue
        focean(J)=1.-FDATA(1,J,2)
  284 continue
c      print *,'FOCEAN'
c      print *,focean
C
cinclude 'DEPTHED'
      open (767,file=fl_ed,form='unformatted')
      open (768,file=fl_zoav,form='unformatted')
      read(767),ED
      read(768),depthml
c     do j=1,nlat
c       depthml(j)=100.
c     enddo
c
      wind_amp=365./12.
      chdif=1.0
!      kvh0 kv for heat
!      rkckh ratio of Kvh to Kvc 
!!  New  kvcmin minimal value of Kvc
      print *,'kvh0=',kvh0,' Kvc/Kvh=',rkckh
!     kvc0=kvh0*rkckh
      kvc0=kvcmin+kvh0*rkckh
      print *,'kvcmin=',kvcmin
      print *,'kvc0=',kvc0
      cfvdif=kvc0/2.5
!     cfvdif miltiplier for GISS diffusion coeficients 
      do j=1,nlat
	ed(j)=ed(j)*cfvdif
        area_ml(j)=dxyp(j)*focean(j)
        volume_ml(j) = area_ml(j)*depthml(j)
      end do    !       j
      print *,'ED from ocm.F'
      print *,ed
c     print *,'ML_VOL'
c     print *,volume_ml
      dzog(1)=10./SQRT(1.7010587)
      dzg(2)=10.
      do l=2,lmo-1
        dzog(L)=dzog(L-1)*1.7010587
        dzg(L+1)=dzg(L)*1.7010587
      end do
      zg(1)=50.
      dzg(1)=100.
      do l=2,lmo
        zg(l)=zg(l-1)+0.5*(dzg(l-1)+dzg(l))
      end do
        print *,'Vertical levels and horizontal diffusion coefficients'
      do l=1,lmo
c       edohd(l)=3.0e4/(zg(l)/zg(1))**0.9
c       edohd(l)=2.5e4/(zg(l)/zg(1))**1.0
        edohd(l)=1.55e4-9.231e3*(atan((zg(l)-300)/50))
        edohd(l)=chdif*edohd(l)
        print *,'l=',l,' Z=',zg(l),' Kh=',edohd(l)
      end do
c     iyeareq=1860
c     co2eq=288.
      print *,' iyeareq=',iyeareq,'  co2eq=',co2eq
	 open(667,file=fl_init_dic,
     &   form='unformatted',status='old')
	 open(668,file=fl_init_alkt,
     &   form='unformatted',status='old')
	 open(669,file=fl_init_salt,
     &   form='unformatted',status='old')
      call init_carb_data(co2eq)
c
C     Do loop for diferent reg
      do nnn=1,nrun
      reg=areg(nnn)
      print *,' nnn=',nnn,' reg=',reg
      if(reg.eq.'equ')then
         open (602,file=ocean_eq,status='new',
     &   form='unformatted')
	 open(685,file=eqclimate,form='unformatted')
	 read(685)tsrfm,wsrfm,asrfm
	 open(603,file=deepo_eq,
     &   form='unformatted',status='new')
	 nyear=nyeareq
         kvc=kvc0
         iyrocu=iyeareq-1
         do j=1,nlat
          Hg(j)=1.0e-8
         enddo
      elseif(reg.eq.'tra')then
         open(686,file=trclimate,form='unformatted',
     &   status='old')
         open(684,file=oco2file,
     &   form='formatted',
     &   status='old')
         if(.not.restart) then
           open (602,file=ocean_tr,status='new',
     &     form='unformatted')
	   open(603,file=deepo_eq,
     &     form='unformatted',status='old')
           rewind 603
           nyear=lyeartr-iyeareq+1
           read(603)iyrocu,kvc
           if(abs(kvc-kvc0).gt.0.1) then
              print *,' Disagreement in Kvc'
              print *,' Kvc=',kvc,' Kvc0=',kvc0
              stop
           endif
           read(603)Hdec
           read(603)Rco2dec
           close (603)
	   open(603,file=deepo_tr,
     &     form='unformatted',status='new')
         else
C      Restart
           open (602,file=ocean_tr,status='old',
     &     form='unformatted')
	   open(603,file=deepo_tr,
     &     form='unformatted',status='old')
           read(603)iyrocu,kvc
           if(abs(kvc-kvc0).gt.0.1) then
              print *,' Disagreement in Kvc'
              print *,' Kvc=',kvc,' Kvc0=',kvc0
              stop
           endif
           read(603)Hdec
           read(603)Rco2dec
           idydata=-iyeareq+iyrocu+1
           nyear=lyeartr-iyrocu+1
           do nd=1,idydata
            do mn=1,12
             read(602)
             read(602)
             read(602)
             read(602)
            enddo
           enddo
         endif
         do j=1,nlat
          Hg(j)=Hdec(j)
         enddo
         print *,' INITIAL RCO2'
         do k=1,lmo
          do j=1,nlat
            Rco2(j,k)=Rco2dec(j,k)
            Ico2(j,k)=Rco2(j,k)*1000
            Rco2dec(j,k)=0.0
          end do
         end do
         print 5003,((Ico2(j,k)
     &   ,j=1,nlat),k=1,LMO)
         open( unit=655,file=co2rfile,
     *   status='old',form='formatted')
         read(655,'(i5,f7.2)'),iiyear,cco2
CCCC
         do j=1,nlat
           read(655,'(i5,12f7.2)'),ii,(co2rat1990(j,m),m=1,12)
           do m=1,12
            co2rat1860(j,m)=1.0
c           print *,co2rat1860(j,m),co2rat1990(j,m)
          enddo
         enddo
         dy=1990.-1860.
      end if
      if(restart)then
        open(676,file=ouptake,form='formatted',status='old')
        do nd=1,idydata
          read(676,'(i4,f11.6)')iyear,co2
        enddo
      else
        open(676,file=ouptake,form='formatted',status='new')
      endif

      fluxavy5=0.
C
C    Start of integration
c     print *,'nyear=',nyear
c     print *,'idydata=',idydata
	if(reg.eq.'tra')then
           do nd=1,iyeareq-1860
            read(684,*)iyear,co2
            print *,'year=',iyear,' co2=',co2
           enddo
         endif
C *********
          if(restart) then
           do nd=1,idydata
            read(684,100)iyear,co2
            do mn=1,12
             read(686,end=500)asrf,tsrf,wsrf,wsrf
            enddo
           enddo
          endif
C *********
      do ny=1,nyear
        fluxavy=0.
        fluxg=0.
cprint *,' ny=',ny
	if(reg.eq.'tra')then
          if (restart)then
           do nd=1,idydata
            read(684,100)iyear,co2
           enddo
          endif
          read(684,*)iyear,co2
          if(ny.eq.1.and.(iyear.ne.iyrocu+1))then
            print *,' Wrong CO2 year'
            print *,' iyear=',iyear,' iyrocu=',iyrocu
            stop
          endif
  100     format(i4,f8.2)
c         print *,iyear,co2
          iyrocu=iyear
	end if
	do m=1,12
c         call zerogary(ncallgary)
          call zerogary
          if(reg.eq.'equ')then
            do j=1,nlat
              t_srf(j)=tsrfm(j,m)+273.16
              w_srf(j)=wsrfm(j,m)
              area_srf(j)=asrfm(j,m)
           end do	!	j
          else
           read(686,end=500)asrf,tsrf,wsrf,wsrf
           do j=1,nlat
             t_srf(j)=tsrf(j)+273.16
             w_srf(j)=wsrf(j)
             area_srf(j)=asrf(j)
           end do  !       j
           if(reg.eq.'tra')then
c           print *,reg,' ',ny,m
             x=(1990.-iyear)/dy
             do j=1,nlat
               pC_atm(j)=co2*1.e-6
               pC_atm(j)=co2*1.e-6
     &          *(x*co2rat1860(j,m)+(1.-x)*co2rat1990(j,m))
             end do
           end if
          end if
          if(ny.eq.1)then
c    print *,' m=',m
c           print *,' tg C'
c           print 5004,(t_srf(j)-273.16,j=1,nlat)
c           print *,' ws'
c           print 5001,w_srf
c           print *,' ar*1.e-13'
c           print 5001,(area_srf(j)*1.e-13,j=1,nlat)
          end if
c         ncallgary=ncallgary+1
!         call carb_mxdlyr_chem(focean)
!         call carb_airsea_flx
          dtco2=3600.*24.*30.
          call carb_chem_ocmip(focean)
          call carb_airsea_flx(dtco2)
          call diffusco2(lmo,nlat,dtco2,0.5,ed,depthml,volume_ml,
     &                dzg,dzog,RCO2)
          call hdocean(RCO2,focean,dxv,dyv,DXYP,DEPTHML,edohd,dtco2)
          call avegary
          if(reg.eq.'equ')then
             if(ny.gt.nyear-50)then
c               call wrgary(ncallgary,fluxao,area_srf)
                call wrgary
                do j=1,nlat
                  Fluxao(j)=Flux_atm_mlav(j)
                end do !        j
                if(m.eq.12)call avedec(ncalldec,Hdec,Rco2dec)
              fluxav=0.
              do j=1,nlat
                 fluxav=fluxav+fluxao(j)
              end do    !       j
              fluxavy=fluxavy+fluxav
              if(m.eq.12)then
              write(676,'(i4,f11.6)'),ny,fluxavy*12.*1.e-15
              endif
       end if
          else 
c    call wrgary(ncallgary,fluxao,area_srf)
	    call wrgary
                do j=1,nlat
                  Fluxao(j)=Flux_atm_mlav(j)
                end do !        j
	    if(reg.eq.'tra')then
              fluxav=0.
              do j=1,nlat
                 fluxav=fluxav+fluxao(j)
              end do    !       j
              fluxavy=fluxavy+fluxav
              if(m.eq.12)then
              write(676,'(i4,f11.6)'),iyear,fluxavy*12.*1.e-15
	      end if
            endif
c           if(ny.gt.nyear-10)then
            if(ny.gt.nyear)then
              if(m.eq.12)then
                 print *,fluxavy*12.*1.e-15
                 fluxavy5=fluxavy5+fluxavy/10.
                 print *,fluxavy5*12.*1.e-15
              endif
	    end if
	    if(ny.eq.nyear.and.m.eq.12) then
               print *,' nyear=',nyear,' m=',m
c              print *,' fluxavy=',fluxavy,fluxavy*chien
c              fluxao is in moles/month
c              fluxav and fluxavy  are in moles/year
               fluxavco2=fluxavy*12.*1.e-15
c              fluxavco2 is in GtC/year (12 is weith of carbon)
                print *,iyear,co2
               print *,' fluxavco2=',fluxavco2,' Gtc/y'
               call avedec(ncalldec,Hdec,Rco2dec)
               call avedec(ncalldec,Hdec,Rco2dec)
	     end if
	  end if
cend if
	end do	!	m
      end do	!	ny
      if(reg.eq.'equ'.or.reg.eq.'tra')then
c       print *,'ncalldec=',ncalldec
c       print *,(Rco2dec(j,1),j=1,nlat)
        do j=1,nlat
         Hdec(j)=Hdec(j)/ncalldec
         do l=1,lmo
           Rco2dec(j,l)=Rco2dec(j,l)/ncalldec
         end do !       l
        end do  !       j
        do l=1,lmo
           Rco2dec(1,l)=Rco2dec(3,l)
           Rco2dec(2,l)=Rco2dec(3,l)
        end do
c       print *,(Rco2dec(j,1),j=1,nlat)
        write(603)iyrocu,kvc
        write(603)Hdec
        write(603)Rco2dec
        ncalldec=0.0
        do j=1,nlat
         Hdec(j)=0.0
         do l=1,lmo
           Rco2dec(j,l)=0.0
         end do !       l
        end do  !       j
	close (103)
      end if
      close (102)
      end do 	! nnn
      print *,'End of run with off line OCM'
      stop 25
  500 print *,' end of file 86'
 5001 format(24f5.2)
 5002 format(24e10.2)
 5003 format(24i5)
 5004 format(24f5.1)
      end
      subroutine avedec(ncalldec,Hdec,Rco2dec)
#include "ODIFF.COM"
!#include "chem_para"
#include "com_carb"
      real  Rco2dec(nlat,lmo),Hdec(nlat)
      ncalldec=ncalldec+1
      do j=1,nlat
         Hdec(j)=Hdec(j)+H(j)
         do l=1,lmo
           Rco2dec(j,l)=Rco2dec(j,l)+rco2(j,l)
         end do !       l
      end do  !       j
      return
      end
