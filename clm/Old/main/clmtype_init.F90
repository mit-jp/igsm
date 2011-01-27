!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clmtype_init
!
! !INTERFACE:
subroutine clmtype_init()
!
! !DESCRIPTION: 
! Initialize clmtype components to signaling nan 
! The following clmtype components should NOT be initialized here
! since they are set in routine clm_map which is called before this
! routine is invoked
! *%area, *%wt, *%wtxy, *%ixy, *%jxy, *%mxy, *%index1d, *%ifspecial
! *%type, *%itypveg, *%ngridcells, *%nlandunits, *%ncolumns, *%npfts
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmpoint
  use nanMod
!
! !ARGUMENTS:
  implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
  integer :: gi,li,ci,pi            !indices
  integer :: ivt                    !vegetation type index
  integer :: begg, endg             !per-process beginning and ending gridcell indices
  integer :: begl, endl             !per-process beginning and ending land unit indices
  integer :: begc, endc             !per-process beginning and ending column indices 
  integer :: begp, endp             !per-process beginning and ending pft indices
  type(gridcell_type), pointer :: g !pointer to derived subtype
  type(landunit_type), pointer :: l !pointer to derived subtype
  type(column_type)  , pointer :: c !pointer to derived subtype
  type(pft_type)     , pointer :: p !pointer to derived subtype
!------------------------------------------------------------------------

  ! Set beginning and ending indices

  begg = grid1d%beg
  endg = grid1d%end
  begl = land1d%beg
  endl = land1d%end
  begc = cols1d%beg
  endc = cols1d%end
  begp = pfts1d%beg
  endp = pfts1d%end

  ! energy balance structures

  clm%mebal%errsoi = nan       
  clm%mebal%errseb = nan       
  clm%mebal%errsol = nan       
  clm%mebal%errlon = nan       
  clm%mebal%acc_errsoi = 0._r8 
  clm%mebal%acc_errseb = 0._r8 
  clm%mebal%acc_errsol = 0._r8 

  do gi = begg,endg
     g => gpoint%grd(gi)%g
     g%gebal%errsoi = nan       
     g%gebal%errseb = nan       
     g%gebal%errsol = nan       
     g%gebal%errlon = nan       
     g%gebal%acc_errsoi = 0._r8 
     g%gebal%acc_errseb = 0._r8 
     g%gebal%acc_errsol = 0._r8 
  end do

  do li = begl,endl
     l => lpoint%lnd(li)%l
     l%lebal%errsoi = nan       
     l%lebal%errseb = nan       
     l%lebal%errsol = nan       
     l%lebal%errlon = nan       
     l%lebal%acc_errsoi = 0._r8 
     l%lebal%acc_errseb = 0._r8 
     l%lebal%acc_errsol = 0._r8 
  end do

  do ci = begc,endc
     c => cpoint%col(ci)%c
     c%cebal%errsoi = nan       
     c%cebal%errseb = nan       
     c%cebal%errsol = nan       
     c%cebal%errlon = nan       
     c%cebal%acc_errsoi = 0._r8 
     c%cebal%acc_errseb = 0._r8 
     c%cebal%acc_errsol = 0._r8 
  end do

  do pi = begp,endp
     p => ppoint%pft(pi)%p
     p%pebal%errsoi = nan       
     p%pebal%errseb = nan       
     p%pebal%errsol = nan       
     p%pebal%errlon = nan       
     p%pebal%acc_errsoi = 0._r8 
     p%pebal%acc_errseb = 0._r8 
     p%pebal%acc_errsol = 0._r8 
  end do

  ! water balance structures

  clm%mwbal%begwb = nan        
  clm%mwbal%endwb = nan        
  clm%mwbal%errh2o = nan       
  clm%mwbal%acc_errh2o = 0._r8 

  do gi = begg,endg
     g => gpoint%grd(gi)%g
     g%gwbal%begwb = nan        
     g%gwbal%endwb = nan        
     g%gwbal%errh2o = nan       
     g%gwbal%acc_errh2o = 0._r8 
  end do

  do li = begl,endl
     l => lpoint%lnd(li)%l
     l%lwbal%begwb = nan        
     l%lwbal%endwb = nan        
     l%lwbal%errh2o = nan       
     l%lwbal%acc_errh2o = 0._r8 
  end do

  do ci = begc,endc
     c => cpoint%col(ci)%c
     c%cwbal%begwb = nan        
     c%cwbal%endwb = nan        
     c%cwbal%errh2o = nan       
     c%cwbal%acc_errh2o = 0._r8 
  end do

  do ci = begc,endc
     c => cpoint%col(ci)%c
     c%cwbal%begwb = nan        
     c%cwbal%endwb = nan        
     c%cwbal%errh2o = nan       
     c%cwbal%acc_errh2o = 0._r8 
  end do

  do pi = begp,endp
     p => ppoint%pft(pi)%p
     p%pwbal%begwb = nan        
     p%pwbal%endwb = nan        
     p%pwbal%errh2o = nan       
     p%pwbal%acc_errh2o = 0._r8 
  end do

  ! pft physical state variables 

  do pi = begp,endp
     p => ppoint%pft(pi)%p
     p%pps%frac_veg_nosno = bigint	   
     p%pps%frac_veg_nosno_alb = bigint  
     p%pps%rsw = nan			
     p%pps%emv = nan			
     p%pps%z0mv = nan		
     p%pps%z0hv = nan		
     p%pps%z0qv = nan		
     p%pps%rootfr(:nlevsoi) = nan
     p%pps%rootr (:nlevsoi) = nan
     p%pps%dewmx = nan		
     p%pps%rssun = nan		
     p%pps%rssha = nan		
     p%pps%laisun = nan		
     p%pps%laisha = nan		
     p%pps%btran = nan		
     p%pps%fsun = nan		
     p%pps%tlai = nan		
     p%pps%tsai = nan		
#if (defined PCP2PFT)
     p%pps%pcp2pft = nan		
#endif
     p%pps%elai = nan		
     p%pps%esai = nan		
     p%pps%igs = nan			
     p%pps%stembio = nan		
     p%pps%rootbio = nan		
     p%pps%fwet = nan		
     p%pps%fdry = nan		
     p%pps%dt_veg = nan		
     p%pps%htop = nan		
     p%pps%hbot = nan		
     p%pps%z0m = nan			
     p%pps%displa = nan		
     p%pps%albd(:numrad) = nan		
     p%pps%albi(:numrad) = nan		
     p%pps%fabd(:numrad) = nan		
     p%pps%fabi(:numrad) = nan		
     p%pps%ftdd(:numrad) = nan		
     p%pps%ftid(:numrad) = nan		
     p%pps%ftii(:numrad) = nan		
     p%pps%ndvi = nan		
     p%pps%u10 = nan			
     p%pps%fv = nan			
     p%pps%ram1 = nan			
  end do

  ! pft ecophysiological constants 

  do ivt = 0,numpft 
     pftcon(ivt)%ncorn = bigint		
     pftcon(ivt)%nwheat = bigint		
     pftcon(ivt)%noveg = bigint	        
     pftcon(ivt)%ntree = bigint		
     pftcon(ivt)%smpmax = -1.5e5           
     pftcon(ivt)%foln = nan	        
     pftcon(ivt)%dleaf = nan		
     pftcon(ivt)%c3psn = nan		
     pftcon(ivt)%vcmx25 = nan		
     pftcon(ivt)%mp = nan			
     pftcon(ivt)%qe25 = nan		
     pftcon(ivt)%xl = nan			
     pftcon(ivt)%rhol(:numrad) = nan        
     pftcon(ivt)%rhos(:numrad) = nan        
     pftcon(ivt)%taul(:numrad) = nan        
     pftcon(ivt)%taus(:numrad) = nan        
     pftcon(ivt)%z0mr = nan		
     pftcon(ivt)%displar = nan		
     pftcon(ivt)%roota_par = nan		
     pftcon(ivt)%rootb_par = nan		
     pftcon(ivt)%sla = nan			
  end do

  ! pft DGVM-specific ecophysiological constants 

  do ivt = 0,numpft 
     dgv_pftcon(ivt)%respcoeff = nan      
     dgv_pftcon(ivt)%flam = nan           
     dgv_pftcon(ivt)%resist = nan         
     dgv_pftcon(ivt)%l_turn = nan         
     dgv_pftcon(ivt)%l_long = nan         
     dgv_pftcon(ivt)%s_turn = nan         
     dgv_pftcon(ivt)%r_turn = nan         
     dgv_pftcon(ivt)%l_cton = nan         
     dgv_pftcon(ivt)%s_cton = nan         
     dgv_pftcon(ivt)%r_cton = nan         
     dgv_pftcon(ivt)%l_morph = nan        
     dgv_pftcon(ivt)%l_phen = nan         
     dgv_pftcon(ivt)%lmtorm  = nan        
     dgv_pftcon(ivt)%crownarea_max = nan  
     dgv_pftcon(ivt)%init_lai = nan       
     dgv_pftcon(ivt)%x  = nan             
     dgv_pftcon(ivt)%tcmin = nan          
     dgv_pftcon(ivt)%tcmax = nan          
     dgv_pftcon(ivt)%gddmin = nan         
     dgv_pftcon(ivt)%twmax = nan          
     dgv_pftcon(ivt)%lm_sapl = nan
     dgv_pftcon(ivt)%sm_sapl = nan
     dgv_pftcon(ivt)%hm_sapl = nan
     dgv_pftcon(ivt)%rm_sapl = nan
     dgv_pftcon(ivt)%tree = .false.
     dgv_pftcon(ivt)%summergreen = .false.
     dgv_pftcon(ivt)%raingreen = .false.
     dgv_pftcon(ivt)%reinickerp = nan
     dgv_pftcon(ivt)%wooddens = nan	 
     dgv_pftcon(ivt)%latosa = nan	    
     dgv_pftcon(ivt)%allom1 = nan	    
     dgv_pftcon(ivt)%allom2 = nan     
     dgv_pftcon(ivt)%allom3 = nan	    
  end do

  ! pft energy state variables 

  do pi = begp,endp
     p => ppoint%pft(pi)%p
     p%pes%t_ref2m = nan         
     p%pes%t_veg = nan           
     p%pes%t_rad_pft = nan	    
  end do


  ! pft water state variables 

  do pi = begp,endp
     p => ppoint%pft(pi)%p
     p%pws%h2ocan = nan  
  end do

  ! pft DGVM state variables 

  do pi = begp,endp
     p => ppoint%pft(pi)%p
     p%pdgvs%agdd0 = nan              
     p%pdgvs%agdd5  = nan             
     p%pdgvs%agddtw = nan             
     p%pdgvs%agdd = nan               
     p%pdgvs%t10 = nan                
     p%pdgvs%t_mo = nan               
     p%pdgvs%t_mo_min = nan           
     p%pdgvs%fnpsn10 = nan            
     p%pdgvs%prec365 = nan            
     p%pdgvs%agdd20 = nan             
     p%pdgvs%tmomin20 = nan           
     p%pdgvs%t10min = nan             
     p%pdgvs%tsoi25 = nan             
     p%pdgvs%annpsn = nan             
     p%pdgvs%annpsnpot = nan          
     p%pdgvs%present = .false.            
     p%pdgvs%dphen = nan              
     p%pdgvs%leafon = nan             
     p%pdgvs%leafof = nan             
     p%pdgvs%nind = nan               
     p%pdgvs%lm_ind = nan             
     p%pdgvs%sm_ind = nan             
     p%pdgvs%hm_ind = nan             
     p%pdgvs%rm_ind = nan             
     p%pdgvs%lai_ind = nan            
     p%pdgvs%fpcinc = nan             
     p%pdgvs%fpcgrid = nan            
     p%pdgvs%crownarea = nan          
     p%pdgvs%bm_inc = nan             
     p%pdgvs%afmicr = nan             
     p%pdgvs%firelength  = nan        
     p%pdgvs%litterag = nan           
     p%pdgvs%litterbg = nan           
     p%pdgvs%cpool_fast = nan         
     p%pdgvs%cpool_slow = nan         
     p%pdgvs%k_fast_ave = nan         
     p%pdgvs%k_slow_ave = nan         
     p%pdgvs%litter_decom_ave = nan   
     p%pdgvs%turnover_ind = nan       
  end do


  ! pft energy flux variables 

  do pi = begp,endp
     p => ppoint%pft(pi)%p
     p%pef%sabg = nan		
     p%pef%sabv = nan		
     p%pef%fsa = nan			
     p%pef%fsr = nan			
     p%pef%parsun = nan      
     p%pef%parsha = nan      
     p%pef%dlrad = nan    
     p%pef%ulrad = nan    
     p%pef%eflx_lh_tot = nan		
     p%pef%eflx_lh_grnd = nan	
     p%pef%eflx_soil_grnd = nan	
     p%pef%eflx_sh_tot = nan 
     p%pef%eflx_sh_grnd = nan   
     p%pef%eflx_sh_veg = nan 
     p%pef%eflx_lh_vege = nan   
     p%pef%eflx_lh_vegt = nan   
     p%pef%cgrnd = nan    
     p%pef%cgrndl = nan      
     p%pef%cgrnds = nan      
     p%pef%eflx_gnet = nan    
     p%pef%dgnetdT = nan      
     p%pef%eflx_lwrad_out = nan	
     p%pef%eflx_lwrad_net = nan	
  end do

  ! pft momentum flux variables 

  do pi = begp,endp
     p => ppoint%pft(pi)%p
     p%pmf%taux = nan  
     p%pmf%tauy = nan  
  end do

  ! pft water flux variables 

  do pi = begp,endp
     p => ppoint%pft(pi)%p
     p%pwf%qflx_prec_intr = nan	
     p%pwf%qflx_prec_grnd = nan	
     p%pwf%qflx_rain_grnd = nan	
     p%pwf%qflx_snow_grnd = nan	
     p%pwf%qflx_snowcap = nan	
     p%pwf%qflx_evap_veg = nan	
     p%pwf%qflx_tran_veg = nan	
     p%pwf%qflx_evap_can = nan       
     p%pwf%qflx_evap_soi = nan	
     p%pwf%qflx_evap_tot = nan	
     p%pwf%qflx_evap_grnd = nan	
     p%pwf%qflx_dew_grnd = nan	
     p%pwf%qflx_sub_snow = nan	
     p%pwf%qflx_dew_snow = nan	
  end do


  ! pft carbon flux variables 

  do pi = begp,endp
     p => ppoint%pft(pi)%p
     p%pcf%psnsun = nan		
     p%pcf%psnsha = nan		
     p%pcf%fpsn = nan		
     p%pcf%frm = nan		
     p%pcf%frmf = nan		
     p%pcf%frms = nan		
     p%pcf%frmr = nan		
     p%pcf%frg = nan		
     p%pcf%dmi = nan		
     p%pcf%fco2 = nan		
     p%pcf%fmicr = nan		
  end do


  ! pft VOC flux variables 

  do pi = begp,endp
     p => ppoint%pft(pi)%p
     p%pvf%vocflx_tot = nan   	
     p%pvf%vocflx(:nvoc) = nan	
     p%pvf%vocflx_1 = nan      	
     p%pvf%vocflx_2 = nan      	
     p%pvf%vocflx_3 = nan      	
     p%pvf%vocflx_4 = nan      	
     p%pvf%vocflx_5 = nan      	
  end do

  ! pft dust flux variables 

  do pi = begp,endp
     p => ppoint%pft(pi)%p
     p%pdf%flx_mss_vrt_dst(:ndst) = nan 
     p%pdf%flx_mss_vrt_dst_tot = nan   
     p%pdf%vlc_trb(:ndst) = nan
     p%pdf%vlc_trb_1 = nan
     p%pdf%vlc_trb_2 = nan
     p%pdf%vlc_trb_3 = nan
     p%pdf%vlc_trb_4 = nan
  end do

  ! column physical state variables 

  do ci = begc,endc
     c => cpoint%col(ci)%c
     c%cps%gwc_thr = nan		
     c%cps%mss_frc_cly_vld = nan	
     c%cps%mbl_bsn_fct = nan		
     c%cps%do_capsnow = .false.
     c%cps%snowdp = nan		
     c%cps%snowage = nan		
     c%cps%frac_sno = nan		
     c%cps%zi(-nlevsno+0:nlevsoi) = nan          
     c%cps%dz(-nlevsno+1:nlevsoi) = nan          
     c%cps%z (-nlevsno+1:nlevsoi) = nan          
     c%cps%frac_iceold(-nlevsno+1:nlevsoi) = nan 
     c%cps%imelt(-nlevsno+1:nlevsoi) = bigint    
     c%cps%eff_porosity(:nlevsoi) = nan           
     c%cps%sfact = nan		
     c%cps%sfactmax = nan		
     c%cps%emg = nan			
     c%cps%z0mg = nan		
     c%cps%z0hg = nan		
     c%cps%z0qg = nan		
     c%cps%htvp = nan		
     c%cps%beta = nan		
     c%cps%zii = nan			
     c%cps%ur = nan			
     c%cps%albgrd(:numrad) = nan 
     c%cps%albgri(:numrad) = nan 
     c%cps%rootr_column(:nlevsoi) = nan		
     c%cps%wf = nan                  
  end do

  ! column energy state variables 

  do ci = begc,endc
     c => cpoint%col(ci)%c
     c%ces%t_grnd = nan 			
     c%ces%dt_grnd	 = nan		        
     c%ces%t_rad_column = nan		
     c%ces%t_soisno(-nlevsno+1:nlevsoi) = nan
     c%ces%t_lake(1:nlevlak)= nan            
     c%ces%tssbef(-nlevsno:nlevsoi) = nan    
     c%ces%t_snow = nan			
     c%ces%thv = nan       			
     c%ces%thm = nan			        
  end do

  ! column water state variables 

  do ci = begc,endc
     c => cpoint%col(ci)%c
     c%cws%h2osno = nan			
     c%cws%h2osoi_liq(-nlevsno+1:nlevsoi) = nan 
     c%cws%h2osoi_ice(-nlevsno+1:nlevsoi) = nan 
     c%cws%h2osoi_vol(:nlevsoi) = nan
     c%cws%h2osno_old = nan			
     c%cws%qg = nan				
     c%cws%dqgdT = nan			
     c%cws%snowice = nan		
     c%cws%snowliq = nan		
  end do

  ! column carbon state variables 

  do ci = begc,endc
     c => cpoint%col(ci)%c
     c%ccs%soilc = nan			
  end do

  ! column energy flux variables 

  do ci = begc,endc
     c => cpoint%col(ci)%c
     c%cef%eflx_snomelt = nan	
     c%cef%eflx_impsoil = nan	
  end do

  ! column water flux variables 

  do ci = begc,endc
     c => cpoint%col(ci)%c
     c%cwf%qflx_infl = nan		
     c%cwf%qflx_surf = nan		
     c%cwf%qflx_drain = nan		
     c%cwf%qflx_top_soil = nan	
     c%cwf%qflx_snomelt = nan	
     c%cwf%qflx_qrgwl = nan		
     c%cwf%qmelt = nan		
     c%cwf%qchan2 = nan              
     c%cwf%qchocn2 = nan             
  end do

  ! land unit physical state variables 

  do li = begl,endl
     l => lpoint%lnd(li)%l
     l%lps%itypwat = bigint		     
     l%lps%isoicol = bigint		     
     l%lps%lakpoi = .false.		     
     l%lps%bsw   (:nlevsoi) = nan          
     l%lps%watsat(:nlevsoi) = nan          
     l%lps%hksat (:nlevsoi) = nan          
     l%lps%sucsat(:nlevsoi) = nan          
     l%lps%csol  (:nlevsoi) = nan          
     l%lps%tkmg  (:nlevsoi) = nan          
     l%lps%tkdry (:nlevsoi) = nan          
     l%lps%tksatu(:nlevsoi) = nan          
     l%lps%smpmin = nan		     
  end do

  ! gridcell gridcell: physical state variables 

  do gi = begg,endg
     g => gpoint%grd(gi)%g
     g%gps%lat = nan		    
     g%gps%lon = nan		    
     g%gps%latdeg = nan	    
     g%gps%londeg = nan	    
     g%gps%wtfact = nan	    
  end do

  ! gridcell: atmosphere -> land state variables 

  do gi = begg,endg
     g => gpoint%grd(gi)%g
     g%a2ls%itypprc = bigint
     g%a2ls%forc_t = nan		
     g%a2ls%forc_u = nan		
     g%a2ls%forc_v = nan		
     g%a2ls%forc_wind = nan       
     g%a2ls%forc_q = nan		
     g%a2ls%forc_hgt = nan	
     g%a2ls%forc_hgt_u = nan	
     g%a2ls%forc_hgt_t = nan	
     g%a2ls%forc_hgt_q = nan	
     g%a2ls%forc_pbot = nan	
     g%a2ls%forc_th = nan		
     g%a2ls%forc_vp = nan		
     g%a2ls%forc_rho = nan	
     g%a2ls%forc_co2 = nan	
     g%a2ls%forc_o2 = nan		
     g%a2ls%forc_psrf = nan	
  end do

  ! gridcell: land -> atmosphere state variables 

  do gi = begg,endg
     g => gpoint%grd(gi)%g
     g%l2as%t_veg = nan	
     g%l2as%t_grnd = nan 	
     g%l2as%t_rad = nan	
     g%l2as%t_ref2m = nan 	
     g%l2as%t_soisno(-nlevsno+1:nlevsoi) = nan
     g%l2as%t_lake(1:nlevlak) = nan
     g%l2as%t_snow = nan
     g%l2as%dt_veg = nan
     g%l2as%dt_grnd = nan
     g%l2as%albd(:numrad) = nan
     g%l2as%albi(numrad) = nan
     g%l2as%albgrd(:numrad) = nan
     g%l2as%albgri(:numrad) = nan
  end do

  ! gridcell: atmosphere -> land flux variables 

  do gi = begg,endg
     g => gpoint%grd(gi)%g
     g%a2lf%forc_lwrad = nan
     g%a2lf%forc_solad(numrad) = nan
     g%a2lf%forc_solai(numrad) = nan
     g%a2lf%forc_solar = nan        
     g%a2lf%forc_rain = nan		
     g%a2lf%forc_snow = nan		
  end do

  ! gridcell: land -> atmosphere flux variables  

  do gi = begg,endg
     g => gpoint%grd(gi)%g
     g%l2af%taux = nan		
     g%l2af%tauy = nan		
     g%l2af%eflx_lwrad_out = nan	
     g%l2af%eflx_lwrad_net = nan	
     g%l2af%eflx_sh_tot = nan	
     g%l2af%eflx_sh_veg = nan	
     g%l2af%eflx_sh_grnd = nan	
     g%l2af%eflx_lh_tot = nan	
     g%l2af%eflx_lh_vege = nan	
     g%l2af%eflx_lh_vegt = nan	
     g%l2af%eflx_lh_grnd = nan	
     g%l2af%eflx_soil_grnd = nan	
     g%l2af%eflx_snomelt = nan	
     g%l2af%eflx_impsoil = nan        	
  end do

  return
end subroutine clmtype_init
