#include <misc.h>
#include <preproc.h>

module clmpoint

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: clmpoint
!
! !DESCRIPTION: 
! Module containing 1d real and integer pointers for clmtype entries
! Use of these "pointer" arrays significantly speeds up model I/O
!
! !USES:
  use clmtype
  use clm_varpar, only : nlevsno, nlevsoi, nlevlak, numrad
!
! !PUBLIC TYPES:
  implicit none
  save
!
! 1d real and integer pointers for clmtype entries
!
  type value_pointer
     real(r8), pointer :: rp
     real(r8), dimension(:), pointer:: rap
  end type value_pointer

  type clm_pointer
     type (value_pointer), dimension(:), pointer :: val
  end type clm_pointer

  integer, parameter :: max_mapflds = 500
  type (clm_pointer) :: clmpointer(max_mapflds)

! 1d gridcell pointer
!
  type grd_single_pointer
     type(gridcell_type), pointer :: g
  end type grd_single_pointer

  type grd_array_pointer
     type (grd_single_pointer), dimension(:), pointer :: grd
  end type grd_array_pointer

  type (grd_array_pointer) :: gpoint
!
! 1d landunit pointer
!
  type lnd_single_pointer
     type(landunit_type), pointer :: l
  end type lnd_single_pointer

  type lnd_array_pointer
     type (lnd_single_pointer), dimension(:), pointer :: lnd
  end type lnd_array_pointer

  type (lnd_array_pointer) :: lpoint
!
! 1d column pointer
!
  type col_single_pointer
     type(column_type), pointer :: c
  end type col_single_pointer

  type col_array_pointer
     type (col_single_pointer), dimension(:), pointer :: col
  end type col_array_pointer

  type (col_array_pointer) :: cpoint
!
! 1d pft pointer
!
  type pft_single_pointer
     type(pft_type), pointer :: p
  end type pft_single_pointer

  type pft_array_pointer
     type (pft_single_pointer), dimension(:), pointer :: pft
  end type pft_array_pointer

  type (pft_array_pointer) :: ppoint
!
! pft indices
!
  integer, parameter :: ip_pps_rsw = 1
  integer, parameter :: ip_pps_emv = ip_pps_rsw + 1 
  integer, parameter :: ip_pps_z0mv = ip_pps_emv + 1 
  integer, parameter :: ip_pps_z0hv = ip_pps_z0mv + 1 
  integer, parameter :: ip_pps_z0qv = ip_pps_z0hv + 1 
  integer, parameter :: ip_pps_dewmx = ip_pps_z0qv + 1 
  integer, parameter :: ip_pps_rssun = ip_pps_dewmx + 1 
  integer, parameter :: ip_pps_rssha = ip_pps_rssun + 1 
  integer, parameter :: ip_pps_laisun = ip_pps_rssha + 1 
  integer, parameter :: ip_pps_laisha = ip_pps_laisun + 1 
  integer, parameter :: ip_pps_btran = ip_pps_laisha + 1 
  integer, parameter :: ip_pps_fsun = ip_pps_btran + 1 
  integer, parameter :: ip_pps_tlai = ip_pps_fsun + 1 
  integer, parameter :: ip_pps_tsai = ip_pps_tlai + 1 
!CAS add PFT Preipitation Distribution
  integer, parameter :: ip_pps_pcp2pft = ip_pps_tsai + 1
!CAS add PFT Preipitation Distribution
  integer, parameter :: ip_pps_elai = ip_pps_pcp2pft + 1
  integer, parameter :: ip_pps_esai = ip_pps_elai + 1 
  integer, parameter :: ip_pps_igs = ip_pps_esai + 1 
  integer, parameter :: ip_pps_stembio = ip_pps_igs + 1 
  integer, parameter :: ip_pps_rootbio = ip_pps_stembio + 1 
  integer, parameter :: ip_pps_fwet = ip_pps_rootbio + 1 
  integer, parameter :: ip_pps_fdry = ip_pps_fwet + 1 
  integer, parameter :: ip_pps_dt_veg = ip_pps_fdry + 1 
  integer, parameter :: ip_pps_htop = ip_pps_dt_veg + 1 
  integer, parameter :: ip_pps_hbot = ip_pps_htop + 1 
  integer, parameter :: ip_pps_z0m = ip_pps_hbot + 1 
  integer, parameter :: ip_pps_displa = ip_pps_z0m + 1 
  integer, parameter :: ip_pps_ndvi = ip_pps_displa + 1 
  integer, parameter :: ip_pps_u10 = ip_pps_ndvi + 1 
  integer, parameter :: ip_pps_fv = ip_pps_u10 + 1  
  integer, parameter :: ip_pps_area = ip_pps_fv + 1
  integer, parameter :: ip_pps_wt = ip_pps_area + 1
  integer, parameter :: ip_pps_wtxy = ip_pps_wt + 1
  integer, parameter :: ip_pps_rootfr = ip_pps_wtxy + 1
  integer, parameter :: ip_pps_rootr = ip_pps_rootfr + 1
  integer, parameter :: ip_pps_albd = ip_pps_rootr + 1 
  integer, parameter :: ip_pps_albi = ip_pps_albd + 1
  integer, parameter :: ip_pps_fabd = ip_pps_albi + 1
  integer, parameter :: ip_pps_fabi = ip_pps_fabd + 1
  integer, parameter :: ip_pps_ftdd = ip_pps_fabi + 1
  integer, parameter :: ip_pps_ftid = ip_pps_ftdd + 1
  integer, parameter :: ip_pps_ftii = ip_pps_ftid + 1
  integer, parameter :: ip_pes_t_ref2m = ip_pps_ftii + 1
  integer, parameter :: ip_pes_t_af = ip_pes_t_ref2m + 1
  integer, parameter :: ip_pes_t_veg = ip_pes_t_af + 1
  integer, parameter :: ip_pes_t_rad_pft = ip_pes_t_veg + 1
  integer, parameter :: ip_pws_h2ocan = ip_pes_t_rad_pft + 1
  integer, parameter :: ip_pef_sabg = ip_pws_h2ocan + 1
  integer, parameter :: ip_pef_sabv = ip_pef_sabg + 1
  integer, parameter :: ip_pef_fsa = ip_pef_sabv + 1
  integer, parameter :: ip_pef_fsr = ip_pef_fsa + 1
  integer, parameter :: ip_pef_parsun = ip_pef_fsr + 1
  integer, parameter :: ip_pef_parsha = ip_pef_parsun + 1
  integer, parameter :: ip_pef_dlrad = ip_pef_parsha + 1
  integer, parameter :: ip_pef_ulrad = ip_pef_dlrad + 1
  integer, parameter :: ip_pef_eflx_lh_tot = ip_pef_ulrad + 1
  integer, parameter :: ip_pef_eflx_lh_grnd = ip_pef_eflx_lh_tot + 1
  integer, parameter :: ip_pef_eflx_soil_grnd = ip_pef_eflx_lh_grnd + 1
  integer, parameter :: ip_pef_eflx_sh_tot = ip_pef_eflx_soil_grnd + 1
  integer, parameter :: ip_pef_eflx_sh_grnd = ip_pef_eflx_sh_tot + 1
  integer, parameter :: ip_pef_eflx_sh_veg = ip_pef_eflx_sh_grnd + 1
  integer, parameter :: ip_pef_eflx_lh_vege = ip_pef_eflx_sh_veg + 1
  integer, parameter :: ip_pef_eflx_lh_vegt = ip_pef_eflx_lh_vege + 1
  integer, parameter :: ip_pef_cgrnd = ip_pef_eflx_lh_vegt + 1
  integer, parameter :: ip_pef_cgrndl = ip_pef_cgrnd + 1
  integer, parameter :: ip_pef_cgrnds = ip_pef_cgrndl + 1
  integer, parameter :: ip_pef_eflx_gnet = ip_pef_cgrnds + 1
  integer, parameter :: ip_pef_dgnetdT = ip_pef_eflx_gnet + 1
  integer, parameter :: ip_pef_eflx_lwrad_out = ip_pef_dgnetdT + 1
  integer, parameter :: ip_pef_eflx_lwrad_net = ip_pef_eflx_lwrad_out + 1
  integer, parameter :: ip_pmf_taux = ip_pef_eflx_lwrad_net + 1
  integer, parameter :: ip_pmf_tauy = ip_pmf_taux + 1
  integer, parameter :: ip_pwf_qflx_prec_intr = ip_pmf_tauy + 1
  integer, parameter :: ip_pwf_qflx_prec_grnd = ip_pwf_qflx_prec_intr + 1
  integer, parameter :: ip_pwf_qflx_rain_grnd = ip_pwf_qflx_prec_grnd + 1
  integer, parameter :: ip_pwf_qflx_snow_grnd = ip_pwf_qflx_rain_grnd + 1
  integer, parameter :: ip_pwf_qflx_snowcap = ip_pwf_qflx_snow_grnd + 1
!CAS add PET, Storm Statistics
  integer, parameter :: ip_pwf_qflx_efpot = ip_pwf_qflx_snowcap + 1
  integer, parameter :: ip_pwf_qflx_strm_int = ip_pwf_qflx_efpot + 1
  integer, parameter :: ip_pwf_qflx_strm_dur = ip_pwf_qflx_strm_int + 1
  integer, parameter :: ip_pwf_qflx_strm_dry = ip_pwf_qflx_strm_dur + 1
  integer, parameter :: ip_pwf_qflx_strms = ip_pwf_qflx_strm_dry + 1
!CAS add PET, Storm Statistics
  integer, parameter :: ip_pwf_qflx_evap_veg = ip_pwf_qflx_strms + 1
  integer, parameter :: ip_pwf_qflx_tran_veg = ip_pwf_qflx_evap_veg + 1
  integer, parameter :: ip_pwf_qflx_evap_soi = ip_pwf_qflx_tran_veg + 1
  integer, parameter :: ip_pwf_qflx_evap_tot = ip_pwf_qflx_evap_soi + 1
  integer, parameter :: ip_pwf_qflx_evap_grnd = ip_pwf_qflx_evap_tot + 1
  integer, parameter :: ip_pwf_qflx_dew_grnd = ip_pwf_qflx_evap_grnd + 1
  integer, parameter :: ip_pwf_qflx_sub_snow = ip_pwf_qflx_dew_grnd + 1
  integer, parameter :: ip_pwf_qflx_dew_snow = ip_pwf_qflx_sub_snow + 1
  integer, parameter :: ip_pwf_qflx_evap_can =  ip_pwf_qflx_dew_snow + 1
  integer, parameter :: ip_pcf_psnsun = ip_pwf_qflx_evap_can + 1
  integer, parameter :: ip_pcf_psnsha = ip_pcf_psnsun + 1
  integer, parameter :: ip_pcf_fpsn = ip_pcf_psnsha + 1
  integer, parameter :: ip_pcf_frm = ip_pcf_fpsn + 1
  integer, parameter :: ip_pcf_frmf = ip_pcf_frm + 1
  integer, parameter :: ip_pcf_frms = ip_pcf_frmf + 1
  integer, parameter :: ip_pcf_frmr = ip_pcf_frms + 1
  integer, parameter :: ip_pcf_frg = ip_pcf_frmr + 1
  integer, parameter :: ip_pcf_dmi = ip_pcf_frg + 1
  integer, parameter :: ip_pcf_fco2 = ip_pcf_dmi + 1
  integer, parameter :: ip_pcf_fmicr = ip_pcf_fco2 + 1
  integer, parameter :: ip_pvf_vocflx_tot = ip_pcf_fmicr + 1
  integer, parameter :: ip_pvf_vocflx_1 = ip_pvf_vocflx_tot + 1
  integer, parameter :: ip_pvf_vocflx_2 = ip_pvf_vocflx_1 + 1
  integer, parameter :: ip_pvf_vocflx_3 = ip_pvf_vocflx_2 + 1
  integer, parameter :: ip_pvf_vocflx_4 = ip_pvf_vocflx_3 + 1
  integer, parameter :: ip_pvf_vocflx_5 = ip_pvf_vocflx_4 + 1
  integer, parameter :: ip_pdf_flx_mss_vrt_dst_tot = ip_pvf_vocflx_5 + 1
  integer, parameter :: ip_pdf_vlc_trb_1 = ip_pdf_flx_mss_vrt_dst_tot + 1
  integer, parameter :: ip_pdf_vlc_trb_2 = ip_pdf_vlc_trb_1 + 1
  integer, parameter :: ip_pdf_vlc_trb_3 = ip_pdf_vlc_trb_2 + 1
  integer, parameter :: ip_pdf_vlc_trb_4 = ip_pdf_vlc_trb_3 + 1
  integer, parameter :: ip_pdgvs_t_mo = ip_pdf_vlc_trb_4 + 1
  integer, parameter :: ip_pdgvs_t_mo_min = ip_pdgvs_t_mo + 1
  integer, parameter :: ip_pdgvs_t10 = ip_pdgvs_t_mo_min + 1
  integer, parameter :: ip_pdgvs_fnpsn10 = ip_pdgvs_t10 + 1
  integer, parameter :: ip_pdgvs_prec365 = ip_pdgvs_fnpsn10 + 1
  integer, parameter :: ip_pdgvs_agdd0 = ip_pdgvs_prec365 + 1
  integer, parameter :: ip_pdgvs_agdd5 = ip_pdgvs_agdd0 + 1
  integer, parameter :: ip_pdgvs_agddtw = ip_pdgvs_agdd5 + 1
  integer, parameter :: ip_pdgvs_agdd = ip_pdgvs_agddtw + 1 
!
! column indices
!
  integer, parameter :: ic_cps_pps_a_rsw = ip_pdgvs_agdd + 1
  integer, parameter :: ic_cps_pps_a_emv = ic_cps_pps_a_rsw + 1
  integer, parameter :: ic_cps_pps_a_z0mv = ic_cps_pps_a_emv + 1
  integer, parameter :: ic_cps_pps_a_z0hv = ic_cps_pps_a_z0mv + 1
  integer, parameter :: ic_cps_pps_a_z0qv = ic_cps_pps_a_z0hv + 1
  integer, parameter :: ic_cps_pps_a_dewmx = ic_cps_pps_a_z0qv + 1
  integer, parameter :: ic_cps_pps_a_rssun = ic_cps_pps_a_dewmx + 1
  integer, parameter :: ic_cps_pps_a_rssha = ic_cps_pps_a_rssun + 1
  integer, parameter :: ic_cps_pps_a_laisun = ic_cps_pps_a_rssha + 1
  integer, parameter :: ic_cps_pps_a_laisha = ic_cps_pps_a_laisun + 1
  integer, parameter :: ic_cps_pps_a_btran = ic_cps_pps_a_laisha + 1
  integer, parameter :: ic_cps_pps_a_fsun = ic_cps_pps_a_btran + 1
  integer, parameter :: ic_cps_pps_a_tlai = ic_cps_pps_a_fsun + 1
  integer, parameter :: ic_cps_pps_a_tsai = ic_cps_pps_a_tlai + 1
!CAS add PFT Preipitation Distribution
  integer, parameter :: ic_cps_pps_a_pcp2pft = ic_cps_pps_a_tsai + 1
!CAS add PFT Preipitation Distribution
  integer, parameter :: ic_cps_pps_a_elai = ic_cps_pps_a_pcp2pft + 1
  integer, parameter :: ic_cps_pps_a_esai = ic_cps_pps_a_elai + 1
  integer, parameter :: ic_cps_pps_a_igs = ic_cps_pps_a_esai + 1
  integer, parameter :: ic_cps_pps_a_stembio = ic_cps_pps_a_igs + 1
  integer, parameter :: ic_cps_pps_a_rootbio = ic_cps_pps_a_stembio + 1
  integer, parameter :: ic_cps_pps_a_fwet = ic_cps_pps_a_rootbio + 1
  integer, parameter :: ic_cps_pps_a_fdry = ic_cps_pps_a_fwet + 1
  integer, parameter :: ic_cps_pps_a_dt_veg = ic_cps_pps_a_fdry + 1
  integer, parameter :: ic_cps_pps_a_htop = ic_cps_pps_a_dt_veg + 1
  integer, parameter :: ic_cps_pps_a_hbot = ic_cps_pps_a_htop + 1
  integer, parameter :: ic_cps_pps_a_z0m = ic_cps_pps_a_hbot + 1
  integer, parameter :: ic_cps_pps_a_displa = ic_cps_pps_a_z0m + 1
  integer, parameter :: ic_cps_pps_a_ndvi = ic_cps_pps_a_displa + 1
  integer, parameter :: ic_cps_pps_a_u10 = ic_cps_pps_a_ndvi + 1
  integer, parameter :: ic_cps_pps_a_fv = ic_cps_pps_a_u10 + 1
  integer, parameter :: ic_cps_pps_a_area = ic_cps_pps_a_fv + 1
  integer, parameter :: ic_cps_pps_a_wt = ic_cps_pps_a_area + 1
  integer, parameter :: ic_cps_pps_a_wtxy = ic_cps_pps_a_wt + 1
  integer, parameter :: ic_cps_vwc_thr = ic_cps_pps_a_wtxy + 1
  integer, parameter :: ic_cps_mss_frc_cly_vld = ic_cps_vwc_thr + 1
  integer, parameter :: ic_cps_mbl_bsn_fct = ic_cps_mss_frc_cly_vld + 1
  integer, parameter :: ic_cps_snowdp = ic_cps_mbl_bsn_fct + 1
  integer, parameter :: ic_cps_snowage = ic_cps_snowdp + 1
  integer, parameter :: ic_cps_frac_sno = ic_cps_snowage + 1
  integer, parameter :: ic_cps_sfact = ic_cps_frac_sno + 1
  integer, parameter :: ic_cps_sfactmax = ic_cps_sfact + 1
  integer, parameter :: ic_cps_emg = ic_cps_sfactmax + 1
  integer, parameter :: ic_cps_z0mg = ic_cps_emg + 1
  integer, parameter :: ic_cps_z0hg = ic_cps_z0mg + 1
  integer, parameter :: ic_cps_z0qg = ic_cps_z0hg + 1
  integer, parameter :: ic_cps_htvp = ic_cps_z0qg + 1
  integer, parameter :: ic_cps_beta = ic_cps_htvp + 1
  integer, parameter :: ic_cps_zii = ic_cps_beta + 1
  integer, parameter :: ic_cps_ur = ic_cps_zii + 1
  integer, parameter :: ic_cps_wf = ic_cps_ur + 1
  integer, parameter :: ic_cps_area = ic_cps_wf + 1
  integer, parameter :: ic_cps_wt = ic_cps_area + 1
  integer, parameter :: ic_cps_wtxy = ic_cps_wt + 1
  integer, parameter :: ic_cps_pps_a_rootfr = ic_cps_wtxy + 1
  integer, parameter :: ic_cps_pps_a_rootr = ic_cps_pps_a_rootfr + 1
  integer, parameter :: ic_cps_frac_iceold = ic_cps_pps_a_rootr + 1
  integer, parameter :: ic_cps_eff_porosity = ic_cps_frac_iceold + 1
  integer, parameter :: ic_cps_rootr_column = ic_cps_eff_porosity + 1
  integer, parameter :: ic_cps_albd = ic_cps_rootr_column + 1
  integer, parameter :: ic_cps_albi = ic_cps_albd + 1
  integer, parameter :: ic_cps_fabd = ic_cps_albi + 1
  integer, parameter :: ic_cps_fabi = ic_cps_fabd + 1
  integer, parameter :: ic_cps_ftdd = ic_cps_fabi + 1
  integer, parameter :: ic_cps_ftid = ic_cps_ftdd + 1
  integer, parameter :: ic_cps_ftii = ic_cps_ftid + 1
  integer, parameter :: ic_cps_albgrd = ic_cps_ftii + 1
  integer, parameter :: ic_cps_albgri = ic_cps_albgrd + 1
  integer, parameter :: ic_ces_pes_a_t_ref2m = ic_cps_albgri + 1
  integer, parameter :: ic_ces_pes_a_t_af = ic_ces_pes_a_t_ref2m + 1
  integer, parameter :: ic_ces_pes_a_t_veg = ic_ces_pes_a_t_af + 1
  integer, parameter :: ic_ces_pes_a_t_rad_pft = ic_ces_pes_a_t_veg + 1
  integer, parameter :: ic_ces_t_grnd = ic_ces_pes_a_t_rad_pft + 1
  integer, parameter :: ic_ces_dt_grnd = ic_ces_t_grnd + 1
  integer, parameter :: ic_ces_t_rad_column = ic_ces_dt_grnd + 1
  integer, parameter :: ic_ces_t_snow = ic_ces_t_rad_column + 1
  integer, parameter :: ic_ces_thv = ic_ces_t_snow + 1
  integer, parameter :: ic_ces_thm = ic_ces_thv + 1
  integer, parameter :: ic_ces_t_lake = ic_ces_thm + 1
  integer, parameter :: ic_ces_t_soisno = ic_ces_t_lake + 1
  integer, parameter :: ic_ces_tssbef = ic_ces_t_soisno + 1
  integer, parameter :: ic_cws_pws_a_h2ocan = ic_ces_tssbef + 1
  integer, parameter :: ic_cws_h2osno = ic_cws_pws_a_h2ocan + 1
  integer, parameter :: ic_cws_qg = ic_cws_h2osno + 1
  integer, parameter :: ic_cws_dqgdT = ic_cws_qg + 1
  integer, parameter :: ic_cws_snowice = ic_cws_dqgdT + 1
  integer, parameter :: ic_cws_snowliq = ic_cws_snowice + 1
  integer, parameter :: ic_cws_h2osoi_liq = ic_cws_snowliq + 1
  integer, parameter :: ic_cws_h2osoi_ice = ic_cws_h2osoi_liq + 1
  integer, parameter :: ic_cws_h2osoi_vol = ic_cws_h2osoi_ice + 1
  integer, parameter :: ic_ccs_soilc = ic_cws_h2osoi_vol + 1
  integer, parameter :: ic_cef_pef_a_sabg = ic_ccs_soilc + 1
  integer, parameter :: ic_cef_pef_a_sabv = ic_cef_pef_a_sabg + 1
  integer, parameter :: ic_cef_pef_a_fsa = ic_cef_pef_a_sabv + 1
  integer, parameter :: ic_cef_pef_a_fsr = ic_cef_pef_a_fsa + 1
  integer, parameter :: ic_cef_pef_a_parsun = ic_cef_pef_a_fsr + 1
  integer, parameter :: ic_cef_pef_a_parsha = ic_cef_pef_a_parsun + 1
  integer, parameter :: ic_cef_pef_a_dlrad = ic_cef_pef_a_parsha + 1
  integer, parameter :: ic_cef_pef_a_ulrad = ic_cef_pef_a_dlrad + 1
  integer, parameter :: ic_cef_pef_a_eflx_lh_tot = ic_cef_pef_a_ulrad + 1
  integer, parameter :: ic_cef_pef_a_eflx_lh_grnd = ic_cef_pef_a_eflx_lh_tot + 1
  integer, parameter :: ic_cef_pef_a_eflx_soil_grnd = ic_cef_pef_a_eflx_lh_grnd + 1
  integer, parameter :: ic_cef_pef_a_eflx_sh_tot = ic_cef_pef_a_eflx_soil_grnd + 1
  integer, parameter :: ic_cef_pef_a_eflx_sh_grnd = ic_cef_pef_a_eflx_sh_tot + 1
  integer, parameter :: ic_cef_pef_a_eflx_sh_veg = ic_cef_pef_a_eflx_sh_grnd + 1
  integer, parameter :: ic_cef_pef_a_eflx_lh_vege = ic_cef_pef_a_eflx_sh_veg + 1
  integer, parameter :: ic_cef_pef_a_eflx_lh_vegt = ic_cef_pef_a_eflx_lh_vege + 1
  integer, parameter :: ic_cef_pef_a_cgrnd = ic_cef_pef_a_eflx_lh_vegt + 1
  integer, parameter :: ic_cef_pef_a_cgrndl = ic_cef_pef_a_cgrnd + 1
  integer, parameter :: ic_cef_pef_a_cgrnds = ic_cef_pef_a_cgrndl + 1
  integer, parameter :: ic_cef_pef_a_eflx_gnet = ic_cef_pef_a_cgrnds + 1
  integer, parameter :: ic_cef_pef_a_dgnetdT = ic_cef_pef_a_eflx_gnet + 1
  integer, parameter :: ic_cef_pef_a_eflx_lwrad_out = ic_cef_pef_a_dgnetdT + 1
  integer, parameter :: ic_cef_pef_a_eflx_lwrad_net = ic_cef_pef_a_eflx_lwrad_out + 1
  integer, parameter :: ic_cef_eflx_snomelt = ic_cef_pef_a_eflx_lwrad_net + 1
  integer, parameter :: ic_cef_eflx_impsoil = ic_cef_eflx_snomelt + 1
  integer, parameter :: ic_cmf_pmf_a_taux = ic_cef_eflx_impsoil + 1
  integer, parameter :: ic_cmf_pmf_a_tauy = ic_cmf_pmf_a_taux + 1
  integer, parameter :: ic_cwf_pwf_a_qflx_prec_intr = ic_cmf_pmf_a_tauy + 1
  integer, parameter :: ic_cwf_pwf_a_qflx_prec_grnd = ic_cwf_pwf_a_qflx_prec_intr + 1
  integer, parameter :: ic_cwf_pwf_a_qflx_rain_grnd = ic_cwf_pwf_a_qflx_prec_grnd + 1
  integer, parameter :: ic_cwf_pwf_a_qflx_snow_grnd = ic_cwf_pwf_a_qflx_rain_grnd + 1
  integer, parameter :: ic_cwf_pwf_a_qflx_snowcap = ic_cwf_pwf_a_qflx_snow_grnd + 1
!CAS add PET and Storm Statistics
  integer, parameter :: ic_cwf_pwf_a_qflx_efpot = ic_cwf_pwf_a_qflx_snowcap + 1
  integer, parameter :: ic_cwf_pwf_a_qflx_strm_int = ic_cwf_pwf_a_qflx_efpot + 1
  integer, parameter :: ic_cwf_pwf_a_qflx_strm_dur = ic_cwf_pwf_a_qflx_strm_int + 1
  integer, parameter :: ic_cwf_pwf_a_qflx_strm_dry = ic_cwf_pwf_a_qflx_strm_dur + 1
  integer, parameter :: ic_cwf_pwf_a_qflx_strms = ic_cwf_pwf_a_qflx_strm_dry + 1
!CAS add PET and Storm Statistics
  integer, parameter :: ic_cwf_pwf_a_qflx_evap_veg = ic_cwf_pwf_a_qflx_strms + 1
  integer, parameter :: ic_cwf_pwf_a_qflx_tran_veg = ic_cwf_pwf_a_qflx_evap_veg + 1
  integer, parameter :: ic_cwf_pwf_a_qflx_evap_soi = ic_cwf_pwf_a_qflx_tran_veg + 1
  integer, parameter :: ic_cwf_pwf_a_qflx_evap_tot = ic_cwf_pwf_a_qflx_evap_soi + 1
  integer, parameter :: ic_cwf_pwf_a_qflx_evap_grnd = ic_cwf_pwf_a_qflx_evap_tot + 1
  integer, parameter :: ic_cwf_pwf_a_qflx_dew_grnd = ic_cwf_pwf_a_qflx_evap_grnd + 1
  integer, parameter :: ic_cwf_pwf_a_qflx_sub_snow = ic_cwf_pwf_a_qflx_dew_grnd + 1
  integer, parameter :: ic_cwf_pwf_a_qflx_dew_snow = ic_cwf_pwf_a_qflx_sub_snow + 1
  integer, parameter :: ic_cwf_pwf_a_qflx_evap_can = ic_cwf_pwf_a_qflx_dew_snow + 1
  integer, parameter :: ic_cwf_qflx_infl = ic_cwf_pwf_a_qflx_evap_can + 1
  integer, parameter :: ic_cwf_qflx_surf = ic_cwf_qflx_infl + 1
  integer, parameter :: ic_cwf_qflx_drain = ic_cwf_qflx_surf + 1
  integer, parameter :: ic_cwf_qflx_top_soil = ic_cwf_qflx_drain + 1
  integer, parameter :: ic_cwf_qflx_snomelt = ic_cwf_qflx_top_soil + 1
  integer, parameter :: ic_cwf_qflx_qrgwl = ic_cwf_qflx_snomelt + 1
  integer, parameter :: ic_cwf_qmelt = ic_cwf_qflx_qrgwl + 1
  integer, parameter :: ic_cwf_qchan2 = ic_cwf_qmelt + 1
  integer, parameter :: ic_cwf_qchocn2 = ic_cwf_qchan2 + 1
  integer, parameter :: ic_ccf_pcf_a_psnsun = ic_cwf_qchocn2 + 1
  integer, parameter :: ic_ccf_pcf_a_psnsha = ic_ccf_pcf_a_psnsun + 1
  integer, parameter :: ic_ccf_pcf_a_fpsn = ic_ccf_pcf_a_psnsha + 1
  integer, parameter :: ic_ccf_pcf_a_frm = ic_ccf_pcf_a_fpsn + 1
  integer, parameter :: ic_ccf_pcf_a_frmf = ic_ccf_pcf_a_frm + 1
  integer, parameter :: ic_ccf_pcf_a_frms = ic_ccf_pcf_a_frmf + 1
  integer, parameter :: ic_ccf_pcf_a_frmr = ic_ccf_pcf_a_frms + 1
  integer, parameter :: ic_ccf_pcf_a_frg = ic_ccf_pcf_a_frmr + 1
  integer, parameter :: ic_ccf_pcf_a_dmi = ic_ccf_pcf_a_frg + 1
  integer, parameter :: ic_ccf_pcf_a_fco2 = ic_ccf_pcf_a_dmi + 1
  integer, parameter :: ic_ccf_pcf_a_fmicr = ic_ccf_pcf_a_fco2 + 1
  integer, parameter :: ic_cebal_errsoi = ic_ccf_pcf_a_fmicr + 1
  integer, parameter :: ic_cebal_errseb = ic_cebal_errsoi + 1
  integer, parameter :: ic_cebal_errsol = ic_cebal_errseb + 1
  integer, parameter :: ic_cwbal_errh2o = ic_cebal_errsol + 1
!
! grid cell indices
!
  integer, parameter :: ig_a2ls_forc_t = ic_cwbal_errh2o + 1
  integer, parameter :: ig_a2ls_forc_u = ig_a2ls_forc_t + 1
  integer, parameter :: ig_a2ls_forc_v = ig_a2ls_forc_u + 1
  integer, parameter :: ig_a2ls_forc_wind = ig_a2ls_forc_v + 1
  integer, parameter :: ig_a2ls_forc_q = ig_a2ls_forc_wind + 1
  integer, parameter :: ig_a2ls_forc_hgt = ig_a2ls_forc_q + 1
  integer, parameter :: ig_a2ls_forc_hgt_u = ig_a2ls_forc_hgt + 1
  integer, parameter :: ig_a2ls_forc_hgt_t = ig_a2ls_forc_hgt_u + 1
  integer, parameter :: ig_a2ls_forc_hgt_q = ig_a2ls_forc_hgt_t + 1
  integer, parameter :: ig_a2ls_forc_pbot = ig_a2ls_forc_hgt_q + 1
  integer, parameter :: ig_a2ls_forc_th = ig_a2ls_forc_pbot + 1
  integer, parameter :: ig_a2ls_forc_vp = ig_a2ls_forc_th + 1
  integer, parameter :: ig_a2ls_forc_rho = ig_a2ls_forc_vp + 1
  integer, parameter :: ig_a2ls_forc_co2 = ig_a2ls_forc_rho + 1
  integer, parameter :: ig_a2ls_forc_o2 = ig_a2ls_forc_co2 + 1
  integer, parameter :: ig_a2ls_forc_psrf = ig_a2ls_forc_o2 + 1
  integer, parameter :: ig_a2lf_forc_lwrad = ig_a2ls_forc_psrf + 1
  integer, parameter :: ig_a2lf_forc_rain = ig_a2lf_forc_lwrad + 1
  integer, parameter :: ig_a2lf_forc_snow = ig_a2lf_forc_rain + 1
!CAS add PET and Storm Statistics
  integer, parameter :: ig_a2lf_forc_strm_int = ig_a2lf_forc_snow + 1
  integer, parameter :: ig_a2lf_forc_strm_dur = ig_a2lf_forc_strm_int + 1
  integer, parameter :: ig_a2lf_forc_strm_dry = ig_a2lf_forc_strm_dur + 1
  integer, parameter :: ig_a2lf_forc_strms = ig_a2lf_forc_strm_dry + 1
!CAS add PET and Storm Statistics
  integer, parameter :: ig_a2lf_forc_solad = ig_a2lf_forc_strms + 1
  integer, parameter :: ig_a2lf_forc_solai = ig_a2lf_forc_solad + 1
  integer, parameter :: ig_a2lf_forc_solar = ig_a2lf_forc_solai + 1
  integer, parameter :: ilast = ig_a2lf_forc_solar
!
! !PUBLIC MEMBER FUNCTIONS:
  public clmpoint_init  ! Initialize pointer arrays and allocate necessary memory
!                              
! !REVISION HISTORY:
! Created by Mariana Vertenstein 
!
!EOP
!----------------------------------------------------------------------- 

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clmpoint_init
!
! !INTERFACE:
  subroutine clmpoint_init()
!
! !DESCRIPTION: 
! Initialize pointer arrays and allocate necessary memory
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
    integer :: gi,li,ci,pi               !indices
    integer :: gindex                    !gridcell index
    integer :: lindex                    !land unit index
    integer :: cindex	 		 !column index
    integer :: pindex			 !pft index
    integer :: begg, endg                !per-process beginning and ending gridcell indices
    integer :: begl, endl                !per-process beginning and ending land unit indices
    integer :: begc, endc                !per-process beginning and ending column indices 
    integer :: begp, endp                !per-process beginning and ending pft indices
    type(gridcell_type), pointer :: g 	 !local pointer to derived subtype
    type(landunit_type), pointer :: l 	 !local pointer to derived subtype
    type(column_type)  , pointer :: c 	 !local pointer to derived subtype
    type(pft_type)     , pointer :: p 	 !local pointer to derived subtype
!------------------------------------------------------------------------

    ! check that ilast is not too large

    if (ilast > max_mapflds) then
       write(6,*)' error: ilast = ',ilast,' greater than max_mapflds= ',max_mapflds
    endif

    ! determine per-process beginning and ending indices

    begp = pfts1d%beg
    endp = pfts1d%end
    begc = cols1d%beg
    endc = cols1d%end
    begl = land1d%beg
    endl = land1d%end
    begg = grid1d%beg
    endg = grid1d%end

    ! allocate gridcell, landunit, column and pft pointers
    
    allocate (gpoint%grd(begg:endg))
    allocate (lpoint%lnd(begl:endl))
    allocate (cpoint%col(begc:endc))
    allocate (ppoint%pft(begp:endp))
    
    ! Set up gridcell, landunit, column and pft pointers
    
    do gi = 1,clm%mps%ngridcells
       g => clm%g(gi)
       gindex = g%gps%index1d
       gpoint%grd(gindex)%g => g
       do li = 1,g%gps%nlandunits
          l => g%l(li)
          lindex = l%lps%index1d
          lpoint%lnd(lindex)%l => l
          do ci = 1,l%lps%ncolumns
             c => l%c(ci)
             cindex = c%cps%index1d
             cpoint%col(cindex)%c => c
             do pi = 1,c%cps%npfts 
                p => c%p(pi)
                pindex = p%pps%index1d
                ppoint%pft(pindex)%p => p
             end do
          end do
       end do
    end do

    ! pft variables allocation

    allocate (clmpointer(ip_pps_rssun)%val(begp:endp))
    allocate (clmpointer(ip_pps_rssha)%val(begp:endp))
    allocate (clmpointer(ip_pps_btran)%val(begp:endp))
#if (defined PCP2PFT)
    allocate (clmpointer(ip_pps_pcp2pft)%val(begp:endp))  !CAS PCP2PFT pointer
#endif
    allocate (clmpointer(ip_pps_elai)%val(begp:endp))
    allocate (clmpointer(ip_pps_esai)%val(begp:endp))
    allocate (clmpointer(ip_pps_ndvi)%val(begp:endp))
    allocate (clmpointer(ip_pps_wtxy)%val(begp:endp))

    allocate (clmpointer(ip_pes_t_ref2m)%val(begp:endp))
    allocate (clmpointer(ip_pes_t_af)%val(begp:endp))
    allocate (clmpointer(ip_pes_t_veg)%val(begp:endp))
    allocate (clmpointer(ip_pws_h2ocan)%val(begp:endp)) 
    allocate (clmpointer(ip_pef_fsa)%val(begp:endp)) 
    allocate (clmpointer(ip_pef_parsun)%val(begp:endp)) 
    allocate (clmpointer(ip_pef_fsr)%val(begp:endp)) 

    allocate (clmpointer(ip_pef_eflx_lh_grnd)%val(begp:endp)) 
    allocate (clmpointer(ip_pef_eflx_soil_grnd)%val(begp:endp)) 
    allocate (clmpointer(ip_pef_eflx_sh_tot)%val(begp:endp)) 
    allocate (clmpointer(ip_pef_eflx_lh_vege)%val(begp:endp)) 
    allocate (clmpointer(ip_pef_eflx_lh_vegt)%val(begp:endp)) 
    allocate (clmpointer(ip_pef_eflx_lwrad_out)%val(begp:endp)) 
    allocate (clmpointer(ip_pef_eflx_lwrad_net)%val(begp:endp)) 

    allocate (clmpointer(ip_pmf_taux)%val(begp:endp)) 
    allocate (clmpointer(ip_pmf_tauy)%val(begp:endp)) 

    allocate (clmpointer(ip_pwf_qflx_prec_intr)%val(begp:endp)) 
    allocate (clmpointer(ip_pwf_qflx_prec_grnd)%val(begp:endp)) 
    allocate (clmpointer(ip_pwf_qflx_efpot)%val(begp:endp)) 
    allocate (clmpointer(ip_pwf_qflx_tran_veg)%val(begp:endp)) 
    allocate (clmpointer(ip_pwf_qflx_evap_soi)%val(begp:endp)) 
    allocate (clmpointer(ip_pwf_qflx_evap_tot)%val(begp:endp)) 
    allocate (clmpointer(ip_pwf_qflx_evap_can)%val(begp:endp)) 

    ! column variables allocation

    allocate (clmpointer(ic_cps_pps_a_rssun)%val(begc:endc)) 
    allocate (clmpointer(ic_cps_pps_a_rssha)%val(begc:endc)) 
    allocate (clmpointer(ic_cps_pps_a_btran)%val(begc:endc)) 
#if (defined PCP2PFT)
    allocate (clmpointer(ic_cps_pps_a_pcp2pft)%val(begc:endc))  !CAS PCP2PFT pointer
#endif
    allocate (clmpointer(ic_cps_pps_a_elai)%val(begc:endc)) 
    allocate (clmpointer(ic_cps_pps_a_esai)%val(begc:endc)) 
    allocate (clmpointer(ic_cps_pps_a_ndvi)%val(begc:endc)) 
    allocate (clmpointer(ic_cps_pps_a_wtxy)%val(begc:endc)) 

    allocate (clmpointer(ic_cps_snowdp)%val(begc:endc)) 
    allocate (clmpointer(ic_cps_frac_sno)%val(begc:endc)) 
    allocate (clmpointer(ic_cps_snowage)%val(begc:endc)) 
    allocate (clmpointer(ic_cps_wtxy)%val(begc:endc)) 

    allocate (clmpointer(ic_ces_pes_a_t_ref2m)%val(begc:endc)) 
    allocate (clmpointer(ic_ces_pes_a_t_af)%val(begc:endc)) 
    allocate (clmpointer(ic_ces_pes_a_t_veg)%val(begc:endc)) 
    allocate (clmpointer(ic_ces_t_grnd)%val(begc:endc)) 
    allocate (clmpointer(ic_ces_t_rad_column)%val(begc:endc)) 
    allocate (clmpointer(ic_ces_t_snow)%val(begc:endc)) 
    allocate (clmpointer(ic_ces_t_lake)%val(begc:endc)) 
    allocate (clmpointer(ic_ces_t_soisno)%val(begc:endc)) 

    allocate (clmpointer(ic_cws_pws_a_h2ocan)%val(begc:endc)) 
    allocate (clmpointer(ic_cws_h2osno)%val(begc:endc)) 
    allocate (clmpointer(ic_cws_snowice)%val(begc:endc)) 
    allocate (clmpointer(ic_cws_snowliq)%val(begc:endc)) 
    allocate (clmpointer(ic_cws_h2osoi_liq)%val(begc:endc)) 
    allocate (clmpointer(ic_cws_h2osoi_ice)%val(begc:endc)) 
    allocate (clmpointer(ic_cws_h2osoi_vol)%val(begc:endc)) 

    allocate (clmpointer(ic_cef_pef_a_fsa)%val(begc:endc)) 
    allocate (clmpointer(ic_cef_pef_a_parsun)%val(begc:endc)) 
    allocate (clmpointer(ic_cef_pef_a_fsr)%val(begc:endc)) 
    allocate (clmpointer(ic_cef_pef_a_eflx_lh_grnd)%val(begc:endc)) 
    allocate (clmpointer(ic_cef_pef_a_eflx_soil_grnd)%val(begc:endc)) 
    allocate (clmpointer(ic_cef_pef_a_eflx_sh_tot)%val(begc:endc)) 
    allocate (clmpointer(ic_cef_pef_a_eflx_lh_vege)%val(begc:endc)) 
    allocate (clmpointer(ic_cef_pef_a_eflx_lh_vegt)%val(begc:endc)) 
    allocate (clmpointer(ic_cef_pef_a_eflx_lwrad_out)%val(begc:endc)) 
    allocate (clmpointer(ic_cef_pef_a_eflx_lwrad_net)%val(begc:endc)) 
    allocate (clmpointer(ic_cef_eflx_snomelt)%val(begc:endc)) 
    allocate (clmpointer(ic_cmf_pmf_a_taux)%val(begc:endc)) 
    allocate (clmpointer(ic_cmf_pmf_a_tauy)%val(begc:endc)) 

    allocate (clmpointer(ic_cwf_pwf_a_qflx_prec_intr)%val(begc:endc)) 
    allocate (clmpointer(ic_cwf_pwf_a_qflx_prec_grnd)%val(begc:endc)) 
    allocate (clmpointer(ic_cwf_pwf_a_qflx_efpot)%val(begc:endc)) 
    allocate (clmpointer(ic_cwf_pwf_a_qflx_tran_veg)%val(begc:endc)) 
    allocate (clmpointer(ic_cwf_pwf_a_qflx_evap_soi)%val(begc:endc)) 
    allocate (clmpointer(ic_cwf_pwf_a_qflx_evap_tot)%val(begc:endc)) 
    allocate (clmpointer(ic_cwf_pwf_a_qflx_evap_can)%val(begc:endc)) 

    allocate (clmpointer(ic_cwf_qflx_infl)%val(begc:endc)) 
    allocate (clmpointer(ic_cwf_qflx_surf)%val(begc:endc)) 
    allocate (clmpointer(ic_cwf_qflx_drain)%val(begc:endc)) 
    allocate (clmpointer(ic_cwf_qflx_snomelt)%val(begc:endc)) 
    allocate (clmpointer(ic_cwf_qflx_qrgwl)%val(begc:endc)) 
    allocate (clmpointer(ic_cwf_qmelt)%val(begc:endc)) 

#if (defined RTM)
    allocate (clmpointer(ic_cwf_qchan2)%val(begc:endc)) 
    allocate (clmpointer(ic_cwf_qchocn2)%val(begc:endc)) 
#endif

    allocate (clmpointer(ic_cebal_errsoi)%val(begc:endc)) 
    allocate (clmpointer(ic_cebal_errseb)%val(begc:endc)) 
    allocate (clmpointer(ic_cebal_errsol)%val(begc:endc)) 
    allocate (clmpointer(ic_cwbal_errh2o)%val(begc:endc)) 

    ! grid cell variable allocation

    allocate (clmpointer(ig_a2ls_forc_t)%val(begg:endg)) 
    allocate (clmpointer(ig_a2ls_forc_u)%val(begg:endg)) 
    allocate (clmpointer(ig_a2ls_forc_v)%val(begg:endg)) 
    allocate (clmpointer(ig_a2ls_forc_wind)%val(begg:endg)) 
    allocate (clmpointer(ig_a2ls_forc_q)%val(begg:endg)) 
    allocate (clmpointer(ig_a2ls_forc_hgt)%val(begg:endg)) 
    allocate (clmpointer(ig_a2ls_forc_hgt_u)%val(begg:endg)) 
    allocate (clmpointer(ig_a2ls_forc_hgt_t)%val(begg:endg)) 
    allocate (clmpointer(ig_a2ls_forc_hgt_q)%val(begg:endg)) 
    allocate (clmpointer(ig_a2ls_forc_pbot)%val(begg:endg)) 
    allocate (clmpointer(ig_a2ls_forc_th)%val(begg:endg)) 
    allocate (clmpointer(ig_a2ls_forc_vp)%val(begg:endg)) 
    allocate (clmpointer(ig_a2ls_forc_rho)%val(begg:endg)) 
    allocate (clmpointer(ig_a2ls_forc_co2)%val(begg:endg)) 
    allocate (clmpointer(ig_a2ls_forc_o2)%val(begg:endg)) 
    allocate (clmpointer(ig_a2ls_forc_psrf)%val(begg:endg)) 
    allocate (clmpointer(ig_a2lf_forc_lwrad)%val(begg:endg)) 
    allocate (clmpointer(ig_a2lf_forc_rain)%val(begg:endg)) 
    allocate (clmpointer(ig_a2lf_forc_snow)%val(begg:endg)) 
    allocate (clmpointer(ig_a2lf_forc_strm_int)%val(begp:endp)) 
    allocate (clmpointer(ig_a2lf_forc_strm_dur)%val(begp:endp)) 
    allocate (clmpointer(ig_a2lf_forc_strm_dry)%val(begp:endp)) 
    allocate (clmpointer(ig_a2lf_forc_strms)%val(begp:endp)) 
    allocate (clmpointer(ig_a2lf_forc_solad)%val(begg:endg))
    allocate (clmpointer(ig_a2lf_forc_solai)%val(begg:endg))
    allocate (clmpointer(ig_a2lf_forc_solar)%val(begg:endg))

    ! Set up clm pointers

    do gi = 1,clm%mps%ngridcells
       g => clm%g(gi)
       gindex = g%gps%index1d
       do li = 1,g%gps%nlandunits
          l => g%l(li)
          lindex = l%lps%index1d
          do ci = 1,l%lps%ncolumns
             c => l%c(ci)
             cindex = c%cps%index1d
             do pi = 1,c%cps%npfts
                p => c%p(pi)
                pindex = p%pps%index1d

                ! pft pps single-level real fields
                clmpointer(ip_pps_rssun)%val(pindex)%rp => p%pps%rssun
                clmpointer(ip_pps_rssha)%val(pindex)%rp => p%pps%rssha
                clmpointer(ip_pps_btran)%val(pindex)%rp => p%pps%btran
#if (defined PCP2PFT)
                clmpointer(ip_pps_pcp2pft)%val(pindex)%rp => p%pps%pcp2pft
#endif
                clmpointer(ip_pps_elai)%val(pindex)%rp => p%pps%elai
                clmpointer(ip_pps_esai)%val(pindex)%rp => p%pps%esai
                clmpointer(ip_pps_ndvi)%val(pindex)%rp => p%pps%ndvi
                clmpointer(ip_pps_wtxy)%val(pindex)%rp => p%pps%wtxy

                ! pft pes single-level  fields
                clmpointer(ip_pes_t_ref2m)%val(pindex)%rp => p%pes%t_ref2m
                clmpointer(ip_pes_t_af)%val(pindex)%rp => p%pes%t_af
                clmpointer(ip_pes_t_veg)%val(pindex)%rp => p%pes%t_veg

                ! pft pws single_level  fields
                clmpointer(ip_pws_h2ocan)%val(pindex)%rp =>  p%pws%h2ocan

                ! pft pef single-level  fields
                clmpointer(ip_pef_fsa)%val(pindex)%rp => p%pef%fsa
                clmpointer(ip_pef_parsun)%val(pindex)%rp => p%pef%parsun
                clmpointer(ip_pef_fsr)%val(pindex)%rp => p%pef%fsr

                clmpointer(ip_pef_eflx_lh_grnd)%val(pindex)%rp => p%pef%eflx_lh_grnd
                clmpointer(ip_pef_eflx_soil_grnd)%val(pindex)%rp => p%pef%eflx_soil_grnd
                clmpointer(ip_pef_eflx_sh_tot)%val(pindex)%rp => p%pef%eflx_sh_tot
                clmpointer(ip_pef_eflx_lh_vege)%val(pindex)%rp => p%pef%eflx_lh_vege
                clmpointer(ip_pef_eflx_lh_vegt)%val(pindex)%rp => p%pef%eflx_lh_vegt
                clmpointer(ip_pef_eflx_lwrad_out)%val(pindex)%rp => p%pef%eflx_lwrad_out
                clmpointer(ip_pef_eflx_lwrad_net)%val(pindex)%rp => p%pef%eflx_lwrad_net

                ! pft pmf single-level  fields
                clmpointer(ip_pmf_taux)%val(pindex)%rp => p%pmf%taux
                clmpointer(ip_pmf_tauy)%val(pindex)%rp => p%pmf%tauy

                ! pft pwf single-level  fields
                clmpointer(ip_pwf_qflx_prec_intr)%val(pindex)%rp => p%pwf%qflx_prec_intr
                clmpointer(ip_pwf_qflx_prec_grnd)%val(pindex)%rp => p%pwf%qflx_prec_grnd
                clmpointer(ip_pwf_qflx_efpot)%val(pindex)%rp => p%pwf%qflx_efpot
                clmpointer(ip_pwf_qflx_tran_veg)%val(pindex)%rp => p%pwf%qflx_tran_veg
                clmpointer(ip_pwf_qflx_evap_soi)%val(pindex)%rp => p%pwf%qflx_evap_soi
                clmpointer(ip_pwf_qflx_evap_tot)%val(pindex)%rp => p%pwf%qflx_evap_tot
                clmpointer(ip_pwf_qflx_evap_can)%val(pindex)%rp => p%pwf%qflx_evap_can

             end do ! end of PFTs loop

             ! column cps single-level fields - pft average fields
             clmpointer(ic_cps_pps_a_rssun)%val(cindex)%rp => c%cps%pps_a%rssun
             clmpointer(ic_cps_pps_a_rssha)%val(cindex)%rp => c%cps%pps_a%rssha
             clmpointer(ic_cps_pps_a_btran)%val(cindex)%rp => c%cps%pps_a%btran
#if (defined PCP2PFT)
             clmpointer(ic_cps_pps_a_pcp2pft)%val(cindex)%rp => c%cps%pps_a%pcp2pft
#endif
             clmpointer(ic_cps_pps_a_elai)%val(cindex)%rp => c%cps%pps_a%elai
             clmpointer(ic_cps_pps_a_esai)%val(cindex)%rp => c%cps%pps_a%esai
             clmpointer(ic_cps_pps_a_ndvi)%val(cindex)%rp => c%cps%pps_a%ndvi
             clmpointer(ic_cps_pps_a_wtxy)%val(cindex)%rp => c%cps%pps_a%wtxy

             ! column cps single-level fields - fields defined at the column level
             clmpointer(ic_cps_snowdp)%val(cindex)%rp => c%cps%snowdp
             clmpointer(ic_cps_frac_sno)%val(cindex)%rp => c%cps%frac_sno
             clmpointer(ic_cps_snowage)%val(cindex)%rp => c%cps%snowage
             clmpointer(ic_cps_wtxy)%val(cindex)%rp => c%cps%wtxy

             ! column ces single-level fields - pft average fields
             clmpointer(ic_ces_pes_a_t_ref2m)%val(cindex)%rp => c%ces%pes_a%t_ref2m
             clmpointer(ic_ces_pes_a_t_af)%val(cindex)%rp => c%ces%pes_a%t_af
             clmpointer(ic_ces_pes_a_t_veg)%val(cindex)%rp => c%ces%pes_a%t_veg

             ! column ces single-level - fields defined at the column level
             clmpointer(ic_ces_t_grnd)%val(cindex)%rp => c%ces%t_grnd
             clmpointer(ic_ces_t_rad_column)%val(cindex)%rp => c%ces%t_rad_column
             clmpointer(ic_ces_t_snow)%val(cindex)%rp => c%ces%t_snow

             ! column ces multi-level (soil) fields
             clmpointer(ic_ces_t_lake)%val(cindex)%rap => c%ces%t_lake
             clmpointer(ic_ces_t_soisno)%val(cindex)%rap => c%ces%t_soisno

             ! column cws single-level fields - pft average fields
             clmpointer(ic_cws_pws_a_h2ocan)%val(cindex)%rp => c%cws%pws_a%h2ocan

             ! column cws single-level fields - fields defined at the column level
             clmpointer(ic_cws_h2osno)%val(cindex)%rp => c%cws%h2osno
             clmpointer(ic_cws_snowice)%val(cindex)%rp => c%cws%snowice
             clmpointer(ic_cws_snowliq)%val(cindex)%rp => c%cws%snowliq

             ! column cws multi-level (soil)  fields
             clmpointer(ic_cws_h2osoi_liq)%val(cindex)%rap => c%cws%h2osoi_liq
             clmpointer(ic_cws_h2osoi_ice)%val(cindex)%rap => c%cws%h2osoi_ice
             clmpointer(ic_cws_h2osoi_vol)%val(cindex)%rap => c%cws%h2osoi_vol

             ! column cef single-level fields - pft average vals
             clmpointer(ic_cef_pef_a_fsa)%val(cindex)%rp => c%cef%pef_a%fsa
             clmpointer(ic_cef_pef_a_parsun)%val(cindex)%rp => c%cef%pef_a%parsun
             clmpointer(ic_cef_pef_a_fsr)%val(cindex)%rp => c%cef%pef_a%fsr
             clmpointer(ic_cef_pef_a_eflx_lh_grnd)%val(cindex)%rp => c%cef%pef_a%eflx_lh_grnd
             clmpointer(ic_cef_pef_a_eflx_soil_grnd)%val(cindex)%rp => c%cef%pef_a%eflx_soil_grnd
             clmpointer(ic_cef_pef_a_eflx_sh_tot)%val(cindex)%rp => c%cef%pef_a%eflx_sh_tot
             clmpointer(ic_cef_pef_a_eflx_lh_vege)%val(cindex)%rp => c%cef%pef_a%eflx_lh_vege
             clmpointer(ic_cef_pef_a_eflx_lh_vegt)%val(cindex)%rp => c%cef%pef_a%eflx_lh_vegt
             clmpointer(ic_cef_pef_a_eflx_lwrad_out)%val(cindex)%rp => c%cef%pef_a%eflx_lwrad_out
             clmpointer(ic_cef_pef_a_eflx_lwrad_net)%val(cindex)%rp => c%cef%pef_a%eflx_lwrad_net

             ! column cef single-level fields - fields defined at the column level
             clmpointer(ic_cef_eflx_snomelt)%val(cindex)%rp => c%cef%eflx_snomelt

             ! column cmf single-level fields - pft average vals
             clmpointer(ic_cmf_pmf_a_taux)%val(cindex)%rp => c%cmf%pmf_a%taux
             clmpointer(ic_cmf_pmf_a_tauy)%val(cindex)%rp => c%cmf%pmf_a%tauy

             ! column cwf single-level fields - pft average vals
             clmpointer(ic_cwf_pwf_a_qflx_prec_intr)%val(cindex)%rp => c%cwf%pwf_a%qflx_prec_intr
             clmpointer(ic_cwf_pwf_a_qflx_prec_grnd)%val(cindex)%rp => c%cwf%pwf_a%qflx_prec_grnd
             clmpointer(ic_cwf_pwf_a_qflx_efpot)%val(cindex)%rp => c%cwf%pwf_a%qflx_efpot
             clmpointer(ic_cwf_pwf_a_qflx_tran_veg)%val(cindex)%rp => c%cwf%pwf_a%qflx_tran_veg
             clmpointer(ic_cwf_pwf_a_qflx_evap_soi)%val(cindex)%rp => c%cwf%pwf_a%qflx_evap_soi
             clmpointer(ic_cwf_pwf_a_qflx_evap_tot)%val(cindex)%rp => c%cwf%pwf_a%qflx_evap_tot
             clmpointer(ic_cwf_pwf_a_qflx_evap_can)%val(cindex)%rp => c%cwf%pwf_a%qflx_evap_can

             ! column cwf single-level fields - fields defined at the column level
             clmpointer(ic_cwf_qflx_infl)%val(cindex)%rp => c%cwf%qflx_infl
             clmpointer(ic_cwf_qflx_surf)%val(cindex)%rp => c%cwf%qflx_surf
             clmpointer(ic_cwf_qflx_drain)%val(cindex)%rp => c%cwf%qflx_drain
             clmpointer(ic_cwf_qflx_snomelt)%val(cindex)%rp => c%cwf%qflx_snomelt
             clmpointer(ic_cwf_qflx_qrgwl)%val(cindex)%rp => c%cwf%qflx_qrgwl
             clmpointer(ic_cwf_qmelt)%val(cindex)%rp => c%cwf%qmelt
#if (defined RTM)
             clmpointer(ic_cwf_qchan2)%val(cindex)%rp => c%cwf%qchan2
             clmpointer(ic_cwf_qchocn2)%val(cindex)%rp => c%cwf%qchocn2
#endif
             clmpointer(ic_cebal_errsoi)%val(cindex)%rp => c%cebal%errsoi
             clmpointer(ic_cebal_errseb)%val(cindex)%rp => c%cebal%errseb
             clmpointer(ic_cebal_errsol)%val(cindex)%rp => c%cebal%errsol
             clmpointer(ic_cwbal_errh2o)%val(cindex)%rp => c%cwbal%errh2o

          end do ! end of columns loop

       end do ! end of landunits loop

       ! gridcell a2ls single-level  fields

       clmpointer(ig_a2ls_forc_t)%val(gindex)%rp => g%a2ls%forc_t
       clmpointer(ig_a2ls_forc_u)%val(gindex)%rp => g%a2ls%forc_u
       clmpointer(ig_a2ls_forc_v)%val(gindex)%rp => g%a2ls%forc_v
       clmpointer(ig_a2ls_forc_wind)%val(gindex)%rp => g%a2ls%forc_wind
       clmpointer(ig_a2ls_forc_q)%val(gindex)%rp => g%a2ls%forc_q
       clmpointer(ig_a2ls_forc_hgt)%val(gindex)%rp => g%a2ls%forc_hgt
       clmpointer(ig_a2ls_forc_hgt_u)%val(gindex)%rp => g%a2ls%forc_hgt_u
       clmpointer(ig_a2ls_forc_hgt_t)%val(gindex)%rp => g%a2ls%forc_hgt_t
       clmpointer(ig_a2ls_forc_hgt_q)%val(gindex)%rp => g%a2ls%forc_hgt_q
       clmpointer(ig_a2ls_forc_pbot)%val(gindex)%rp => g%a2ls%forc_pbot
       clmpointer(ig_a2ls_forc_th)%val(gindex)%rp => g%a2ls%forc_th
       clmpointer(ig_a2ls_forc_vp)%val(gindex)%rp => g%a2ls%forc_vp
       clmpointer(ig_a2ls_forc_rho)%val(gindex)%rp => g%a2ls%forc_rho
       clmpointer(ig_a2ls_forc_co2)%val(gindex)%rp => g%a2ls%forc_co2
       clmpointer(ig_a2ls_forc_o2)%val(gindex)%rp => g%a2ls%forc_o2
       clmpointer(ig_a2ls_forc_psrf)%val(gindex)%rp => g%a2ls%forc_psrf

       ! gridcell a2lf single-level  fields
       
       clmpointer(ig_a2lf_forc_lwrad)%val(gindex)%rp => g%a2lf%forc_lwrad
       clmpointer(ig_a2lf_forc_rain)%val(gindex)%rp => g%a2lf%forc_rain
       clmpointer(ig_a2lf_forc_snow)%val(gindex)%rp => g%a2lf%forc_snow
       clmpointer(ig_a2lf_forc_strm_int)%val(gindex)%rp => g%a2lf%forc_strm_int
       clmpointer(ig_a2lf_forc_strm_dur)%val(gindex)%rp => g%a2lf%forc_strm_dur
       clmpointer(ig_a2lf_forc_strm_dry)%val(gindex)%rp => g%a2lf%forc_strm_dry
       clmpointer(ig_a2lf_forc_strms)%val(gindex)%rp => g%a2lf%forc_strms
       clmpointer(ig_a2lf_forc_solar)%val(gindex)%rp => g%a2lf%forc_solar

       ! gridcell a2lf multi-level (rad)  fields
       
       clmpointer(ig_a2lf_forc_solad)%val(gindex)%rap => g%a2lf%forc_solad
       clmpointer(ig_a2lf_forc_solai)%val(gindex)%rap => g%a2lf%forc_solai

    end do ! end of gridcells loop

    return
  end subroutine clmpoint_init

end module clmpoint
