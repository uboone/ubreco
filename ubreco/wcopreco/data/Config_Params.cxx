#include "Config_Params.h"

namespace wcopreco {
  Config_Params::Config_Params()
  {
  }

  void Config_Params::set_num_channels(int n) {
      _cfg_deconvolver._num_channels = n ;
      _cfg_flashesbeam._num_channels = n ;
      _cfg_hitfinder_beam._num_channels = n ;
      _cfg_opflash._num_channels = n ;
      _cfg_saturation_merger._num_channels = n ;


  }

  void Config_Params::set_tick_width_us(float width) {
      _cfg_deconvolver._tick_width_us = width ;
      _cfg_hitfinder_beam._tick_width_us = width ;
      _cfg_opflash._tick_width_us = width ;
      _cfg_saturation_merger._tick_width_us = width ;


  }

  void Config_Params::set_rebin_frac(int n) {
      _cfg_flashesbeam._rebin_frac = n ;
      _cfg_hitfinder_beam._rebin_frac = n ;
      _cfg_opflash._rebin_frac = n ;

  }

  void Config_Params::set_nbins_beam(int n) {
      _cfg_deconvolver._nbins_beam = n ;
      _cfg_flashesbeam._nbins_beam = n ;
      _cfg_hitfinder_beam._nbins_beam = n ;
  }

  void Config_Params::set_baseline_difference_max(double p) {
      _cfg_deconvolver._baseline_difference_max = p ;
      _cfg_saturation_merger._baseline_difference_max = p ;
  }

  void Config_Params::set_nbins_baseline_search(int n) {
      _cfg_deconvolver._nbins_baseline_search = n ;
      _cfg_saturation_merger._nbins_baseline_search = n ;
  }

  void Config_Params::set_baseline_default(int baseline) {
      _cfg_cophit._baseline_default = baseline ;
      _cfg_saturation_merger._baseline_default = baseline ;
  }

  void Config_Params::set_nbins_cosmic(int n) {
      _cfg_cophit._nbins_cosmic = n ;
  }

  void Config_Params::set_COphit_integral_thresh(double thresh) {
      _cfg_cophit._COphit_integral_thresh = thresh ;
  }

  void Config_Params::set_COphit_baseline_diff_thresh(double thresh) {
      _cfg_cophit._COphit_baseline_diff_thresh = thresh ;
  }

  void Config_Params::set_pe_factor(double factor) {
      _cfg_cophit._pe_factor = factor ;
  }


  void Config_Params::set_Baseline_uncertainty(double unc) {
      _cfg_cophit._Baseline_uncertainty = unc ;
  }

  void Config_Params::set_Baseline_unc_bad_baseline(double unc) {
      _cfg_cophit._Baseline_unc_bad_baseline = unc ;
  }

  void Config_Params::set_cal_integral_p0(double p) {
      _cfg_cophit._cal_integral_p0 = p ;
  }

  void Config_Params::set_cal_integral_p1(double p) {
      _cfg_cophit._cal_integral_p1 = p ;
  }

  void Config_Params::set_cal_integral_p2(double p) {
      _cfg_cophit._cal_integral_p2 = p ;
  }

  void Config_Params::set_cal_integral_p3(double p) {
      _cfg_cophit._cal_integral_p3 = p ;
  }

  void Config_Params::set_cal_integral_p4(double p) {
      _cfg_cophit._cal_integral_p4 = p ;
  }

  void Config_Params::set_cal_integral_p5(double p) {
      _cfg_cophit._cal_integral_p5 = p ;
  }

  void Config_Params::set_channel_status_v(std::vector<bool> is_ch_good_v) {
      _cfg_cophit._channel_status_v = is_ch_good_v ;
  }

  void Config_Params::set_high_freq_p0(double p) {
      _cfg_deconvolver._high_freq_p0 = p ;
  }

  void Config_Params::set_high_freq_p1(double p) {
      _cfg_deconvolver._high_freq_p1 = p ;
  }

  void Config_Params::set_latelight_filter_p0(double p) {
      _cfg_deconvolver._latelight_filter_p0 = p ;
  }

  void Config_Params::set_latelight_filter_p1(double p) {
      _cfg_deconvolver._latelight_filter_p1 = p ;
  }

  void Config_Params::set_latelight_filter_p2(double p) {
      _cfg_deconvolver._latelight_filter_p2 = p ;
  }

  void Config_Params::set_baseline_safety_subtraction(double sub) {
      _cfg_deconvolver._baseline_safety_subtraction = sub ;
  }

  void Config_Params::set_xq(double x) {
      _cfg_deconvolver._xq = x ;
  }

  void Config_Params::set_xq_diff(double diff) {
      _cfg_deconvolver._xq_diff = diff ;
  }

  void Config_Params::set_n_bins_end_wfm(int n) {
      _cfg_deconvolver._n_bins_end_wfm = n ;
  }

  void Config_Params::set_small_content_bump(double bump) {
      _cfg_deconvolver._small_content_bump = bump ;
  }

  void Config_Params::set_bflash_pe_thresh(double thresh) {
      _cfg_flashesbeam._bflash_pe_thresh = thresh ;
  }

  void Config_Params::set_bflash_mult_thresh(double thresh) {
      _cfg_flashesbeam._bflash_mult_thresh = thresh ;
  }

  void Config_Params::set_bflash_bin_diff_p0(int p) {
      _cfg_flashesbeam._bflash_bin_diff_p0 = p ;
  }

  void Config_Params::set_bflash_bin_diff_p1(int p) {
      _cfg_flashesbeam._bflash_bin_diff_p1 = p ;
  }

  void Config_Params::set_bflash_bin_diff_p2(int p) {
      _cfg_flashesbeam._bflash_bin_diff_p2 = p ;
  }

  void Config_Params::set_KS_test_thresh(double thresh) {
      _cfg_flashesbeam._KS_test_thresh = thresh ;
  }

  void Config_Params::set_bflash_bin_start_cushion(int cushion) {
      _cfg_flashesbeam._bflash_bin_start_cushion = cushion ;
  }

  void Config_Params::set_flash_filter_time_thresh(double thresh) {
      _cfg_flashfiltering._flash_filter_time_thresh = thresh ;
  }

  void Config_Params::set_flash_filter_pe_thresh(double thresh) {
      _cfg_flashfiltering._flash_filter_pe_thresh = thresh ;
  }

  void Config_Params::set_do_swap_channels(bool do_swap) {
      _cfg_flashfiltering._do_swap_channels = do_swap ;
  }

  void Config_Params::set_l1_content_thresh(double thresh) {
      _cfg_hitfinder_beam._l1_content_thresh = thresh ;
  }

  void Config_Params::set_frac_G_t2_first(double frac) {
      _cfg_hitfinder_beam._frac_G_t2_first = frac ;
  }

  void Config_Params::set_frac_G_sametime(double frac) {
      _cfg_hitfinder_beam._frac_G_sametime = frac ;
  }

  void Config_Params::set_G_p0(double p) {
      _cfg_hitfinder_beam._G_p0 = p ;
  }

  void Config_Params::set_G_p1(double p) {
      _cfg_hitfinder_beam._G_p1 = p ;
  }

  void Config_Params::set_G_p2(double p) {
      _cfg_hitfinder_beam._G_p2 = p ;
  }
  void Config_Params::set_Lasso_p0(double p) {
      _cfg_hitfinder_beam._Lasso_p0 = p ;
  }

  void Config_Params::set_Lasso_p1(int p) {
      _cfg_hitfinder_beam._Lasso_p1 = p ;
  }

  void Config_Params::set_Lasso_p2(double p) {
      _cfg_hitfinder_beam._Lasso_p2 = p ;
  }
  void Config_Params::set_totPE_v_thresh(double thresh) {
      _cfg_hitfinder_beam._totPE_v_thresh = thresh ;
  }

  void Config_Params::set_mult_v_thresh(double thresh) {
      _cfg_hitfinder_beam._mult_v_thresh = thresh ;
  }

  void Config_Params::set_l1_mult_v_thresh(double thresh) {
      _cfg_hitfinder_beam._l1_mult_v_thresh = thresh ;
  }

  void Config_Params::set_ophit_group_t_diff_max(double max) {
      _cfg_hitfinder_cosmic._ophit_group_t_diff_max = max ;
  }

  void Config_Params::set_PE_err_cosmic(double err) {
      _cfg_opflash._PE_err_cosmic = err ;
  }

  void Config_Params::set_PE_subtract(double sub) {
      _cfg_opflash._PE_subtract = sub ;
  }

  void Config_Params::set_flash_low_time_cushion(int bin) {
      _cfg_opflash._flash_low_time_cushion = bin ;
  }

  void Config_Params::set_flash_high_time_cushion(int bin) {
      _cfg_opflash._flash_high_time_cushion = bin ;
  }

  void Config_Params::set_PE_err_beam(double err) {
      _cfg_opflash._PE_err_beam = err ;
  }

  void Config_Params::set_mult_content_thresh(double thresh) {
      _cfg_opflash._mult_content_thresh = thresh ;
  }

  void Config_Params::set_mult_required(int req) {
      _cfg_opflash._mult_required = req ;
  }

  void Config_Params::set_PE_err_stat_beam(double err) {
      _cfg_opflash._PE_err_stat_beam = err ;
  }

  void Config_Params::set_PE_err_unc_beam(double err) {
      _cfg_opflash._PE_err_unc_beam = err ;
  }

  void Config_Params::set_addl1_pe_thresh(double thresh) {
      _cfg_opflash._addl1_pe_thresh = thresh ;
  }

  void Config_Params::set_addl1_mult_thresh(double thresh) {
      _cfg_opflash._addl1_mult_thresh = thresh ;
  }

  void Config_Params::set_sat_threshold(short thresh) {
      _cfg_saturation_merger._sat_threshold = thresh ;
  }

  void Config_Params::set_cosmic_tick_window(int n) {
      _cfg_saturation_merger._cosmic_tick_window = n ;
  }

  void Config_Params::set_low_bound_baseline_search(double d) {
      _cfg_saturation_merger._low_bound_baseline_search = d ;
  }

  void Config_Params::set_high_bound_baseline_search(double d) {
      _cfg_saturation_merger._high_bound_baseline_search = d ;
  }

  void Config_Params::set_nbins_saturation_threshold(int n) {
      _cfg_saturation_merger._nbins_saturation_threshold = n ;
  }

  void Config_Params::set_scaling_by_channel(std::vector<float> scales_v) {
      _cfg_saturation_merger._scaling_by_channel = scales_v ;
  }

  void Config_Params::set_rc_tau_badch(double tau) {
      _cfg_ub_rc._rc_tau_badch = tau ;
  }

  void Config_Params::set_rc_tau_goodch(double tau) {
      _cfg_ub_rc._rc_tau_goodch = tau ;
  }

  void Config_Params::set_spe_p0(double p) {
      _cfg_ub_spe._spe_p0 = p ;
  }

  void Config_Params::set_spe_p1(double p) {
      _cfg_ub_spe._spe_p1 = p ;
  }

  //function to check that common parameters are all filled the same!
void Config_Params::Check_common_parameters(){
  //num_channels
  if(_cfg_deconvolver._num_channels == _cfg_flashesbeam._num_channels
    && _cfg_deconvolver._num_channels == _cfg_hitfinder_beam._num_channels
    && _cfg_deconvolver._num_channels == _cfg_opflash._num_channels
    && _cfg_deconvolver._num_channels == _cfg_saturation_merger._num_channels){}
  else{
    std::cout << "ERROR!!! DIFFERENT VALUES SET FOR _num_channels: \n";
    std::cout << "Deconvolver " << _cfg_deconvolver._num_channels <<"\n";
    std::cout << "flashesbeam "<< _cfg_flashesbeam._num_channels <<"\n";
    std::cout << "hitfinder " << _cfg_hitfinder_beam._num_channels <<"\n";
    std::cout << "opflash "<< _cfg_opflash._num_channels <<"\n";
    std::cout << "saturation_merger " << _cfg_saturation_merger._num_channels <<"\n";
  }

  //tick width
  if( _cfg_deconvolver._tick_width_us == _cfg_hitfinder_beam._tick_width_us
    && _cfg_deconvolver._tick_width_us == _cfg_opflash._tick_width_us
    && _cfg_deconvolver._tick_width_us == _cfg_saturation_merger._tick_width_us){}
  else{
    std::cout << "ERROR!!! DIFFERENT VALUES SET FOR _tick_width_us: \n";
    std::cout << "Deconvolver "<<_cfg_deconvolver._tick_width_us <<"\n";
    std::cout << "hitfinder " << _cfg_hitfinder_beam._tick_width_us <<"\n";
    std::cout <<  "opflash "<< _cfg_opflash._tick_width_us <<"\n";
    std::cout << "saturation_merger " << _cfg_saturation_merger._tick_width_us <<"\n";
  }

  //_rebin_frac
 if(_cfg_flashesbeam._rebin_frac == _cfg_hitfinder_beam._rebin_frac
   && _cfg_flashesbeam._rebin_frac == _cfg_opflash._rebin_frac){}
 else{
   std::cout << "ERROR!!! DIFFERENT VALUES SET FOR _rebin_frac: \n";
   std::cout << "flashesbeam "<<_cfg_flashesbeam._rebin_frac <<"\n";
   std::cout << "hitfinder_beam " << _cfg_hitfinder_beam._rebin_frac <<"\n";
   std::cout <<  "opflash "<< _cfg_opflash._rebin_frac<<"\n";
 }

 //_nbins_beam
 if(_cfg_flashesbeam._nbins_beam == _cfg_hitfinder_beam._nbins_beam
   && _cfg_flashesbeam._nbins_beam == _cfg_deconvolver._nbins_beam){}
 else{
   std::cout << "ERROR!!! DIFFERENT VALUES SET FOR _nbins_beam: \n";
   std::cout << "flashesbeam "<<_cfg_flashesbeam._nbins_beam <<"\n";
   std::cout << "hitfinder_beam " << _cfg_hitfinder_beam._nbins_beam <<"\n";
   std::cout <<  "deconvolver "<< _cfg_deconvolver._nbins_beam<<"\n";
 }

 //_baseline_difference_max
 if(_cfg_deconvolver._baseline_difference_max == _cfg_saturation_merger._baseline_difference_max){}
 else{
   std::cout << "ERROR!!! DIFFERENT VALUES SET FOR _baseline_difference_max: \n";
   std::cout << "saturation_merger " <<_cfg_saturation_merger._baseline_difference_max  <<"\n";
   std::cout <<  "deconvolver "<< _cfg_deconvolver._baseline_difference_max<<"\n";
 }

 //_nbins_baseline_search
 if(_cfg_deconvolver._nbins_baseline_search  == _cfg_saturation_merger._nbins_baseline_search ){}
 else{
   std::cout << "ERROR!!! DIFFERENT VALUES SET FOR _nbins_baseline_search : \n";
   std::cout << "saturation_merger " <<_cfg_saturation_merger._nbins_baseline_search  <<"\n";
   std::cout <<  "deconvolver "<< _cfg_deconvolver._nbins_baseline_search <<"\n";
 }

 //_baseline_default
 if(_cfg_cophit._baseline_default  == _cfg_saturation_merger._baseline_default ){}
 else{
   std::cout << "ERROR!!! DIFFERENT VALUES SET FOR _baseline_default : \n";
   std::cout << "saturation_merger " <<_cfg_saturation_merger._baseline_default  <<"\n";
   std::cout <<  "deconvolver "<< _cfg_cophit._baseline_default <<"\n";
 }

 /// ch status length
 if(int(_cfg_cophit._channel_status_v.size()) != _cfg_deconvolver._num_channels ){
   std::cout << "ERROR VECTOR OF CHANNEL STATUS ISN'T RIGHT LENGTH!! \n";
   std::cout << "num channels " << _cfg_deconvolver._num_channels << "\n";
   std::cout << "length of channel status " << _cfg_cophit._channel_status_v.size() << "\n";

 }

}

}
