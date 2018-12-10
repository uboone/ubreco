#ifndef CONFIG_PARAMS_H
#define CONFIG_PARAMS_H

//Include other Configs:
#include "Config_COpHit.h"
#include "Config_Deconvolver.h"
#include "Config_FlashesBeam.h"
#include "Config_FlashFiltering.h"
#include "Config_Hitfinder_Beam.h"
#include "Config_Hitfinder_Cosmic.h"
#include "Config_Opflash.h"
#include "Config_Saturation_Merger.h"
#include "Config_UB_rc.h"
#include "Config_UB_spe.h"

//Include basics:
#include <vector>
#include <iostream>
#include <sstream>

namespace wcopreco {
    class Config_Params {
    public:
      Config_Params();
      ~Config_Params() {};
      //Common
      void set_num_channels(int n);
      void set_nbins_beam(int n);
      void set_nbins_cosmic(int n);
      void set_tick_width_us(float width);
      void set_rebin_frac(int frac);
      void set_baseline_difference_max(double dif);
      void set_baseline_default( int baseline);
      //COphit
      void set_COphit_integral_thresh(double thresh);
      void set_COphit_baseline_diff_thresh(double thresh);
      void set_pe_factor(double factor);
      void set_Baseline_uncertainty(double unc);
      void set_Baseline_unc_bad_baseline(double unc);
      void set_cal_integral_p0(double p);
      void set_cal_integral_p1(double p);
      void set_cal_integral_p2(double p);
      void set_cal_integral_p3(double p);
      void set_cal_integral_p4(double p);
      void set_cal_integral_p5(double p);
      void set_channel_status_v(std::vector<bool> is_ch_bad_v );

      //Deconvolver
      void set_high_freq_p0(double p);
      void set_high_freq_p1(double p);
      void set_latelight_filter_p0(double p);
      void set_latelight_filter_p1(double p);
      void set_latelight_filter_p2(double p);
      void set_baseline_safety_subtraction(double sub);
      void set_xq(double x);
      void set_xq_diff(double x);
      void set_n_bins_end_wfm(int n);
      void set_small_content_bump(double bump);
      void set_nbins_baseline_search(int n);
      //Flashesbeam
      void set_bflash_pe_thresh(double thresh);
      void set_bflash_mult_thresh(double thresh);
      void set_bflash_bin_diff_p0(int p);
      void set_bflash_bin_diff_p1(int p);
      void set_bflash_bin_diff_p2(int p);
      void set_KS_test_thresh(double thresh);
      void set_bflash_bin_start_cushion(int cushion);
      //FlashFiltering
      void set_flash_filter_time_thresh(double thresh);
      void set_flash_filter_pe_thresh(double thresh);
      void set_do_swap_channels(bool do_swap);
      //Hitfinder_Beam
      void set_l1_content_thresh(double thresh);
      void set_frac_G_t2_first(double frac);
      void set_frac_G_sametime(double frac);
      void set_G_p0(double p);
      void set_G_p1(double p);
      void set_G_p2(double p);
      void set_Lasso_p0(double p);
      void set_Lasso_p1(int p);
      void set_Lasso_p2(double p);
      void set_totPE_v_thresh(double thresh);
      void set_mult_v_thresh(double thresh);
      void set_l1_mult_v_thresh(double thresh);
      //Hitfinder_cosmic
      void set_ophit_group_t_diff_max(double max);
      //Opflash
      void set_PE_err_cosmic(double err);
      void set_PE_subtract(double sub);
      void set_flash_low_time_cushion(int bin);
      void set_flash_high_time_cushion(int bin);
      void set_PE_err_beam(double err);
      void set_mult_content_thresh(double thresh);
      void set_mult_required(int req);
      void set_PE_err_stat_beam(double err);
      void set_PE_err_unc_beam(double err);
      void set_addl1_pe_thresh(double thresh);
      void set_addl1_mult_thresh(double thresh);
      // Saturation Merger
      void set_sat_threshold(short thresh);
      void set_cosmic_tick_window(int n);
      void set_low_bound_baseline_search(double d);
      void set_high_bound_baseline_search(double d);
      void set_nbins_saturation_threshold(int n);
      void set_scaling_by_channel(std::vector<float> scales_v);
      // UB_rc
      void set_rc_tau_badch(double tau);
      void set_rc_tau_goodch(double tau);
      // UB_spe
      void set_spe_p0(double p);
      void set_spe_p1(double p);


      void _set_cfg_cophit(Config_COpHit cfg){_cfg_cophit = cfg;}
      Config_COpHit _get_cfg_cophit(){return _cfg_cophit;}

      void _set_cfg_deconvolver(Config_Deconvolver cfg){_cfg_deconvolver = cfg;}
      Config_Deconvolver _get_cfg_deconvolver(){return _cfg_deconvolver;}

      void _set_cfg_flashesbeam(Config_FlashesBeam cfg){_cfg_flashesbeam = cfg;}
      Config_FlashesBeam _get_cfg_flashesbeam(){return _cfg_flashesbeam;}

      void _set_cfg_flashfiltering(Config_FlashFiltering cfg){_cfg_flashfiltering = cfg;}
      Config_FlashFiltering _get_cfg_flashfiltering(){return _cfg_flashfiltering;}

      void _set_cfg_hitfinder_beam(Config_Hitfinder_Beam cfg){_cfg_hitfinder_beam = cfg;}
      Config_Hitfinder_Beam _get_cfg_hitfinder_beam(){return _cfg_hitfinder_beam;}

      void _set_cfg_hitfinder_cosmic(Config_Hitfinder_Cosmic cfg){_cfg_hitfinder_cosmic = cfg;}
      Config_Hitfinder_Cosmic _get_cfg_hitfinder_cosmic(){return _cfg_hitfinder_cosmic;}

      void _set_cfg_opflash(Config_Opflash cfg){_cfg_opflash = cfg;}
      Config_Opflash _get_cfg_opflash(){return _cfg_opflash;}

      void _set_cfg_saturation_merger(Config_Saturation_Merger cfg){_cfg_saturation_merger = cfg;}
      Config_Saturation_Merger _get_cfg_saturation_merger(){return _cfg_saturation_merger;}

      void _set_cfg_ub_rc(Config_UB_rc cfg){_cfg_ub_rc = cfg;}
      Config_UB_rc _get_cfg_ub_rc(){return _cfg_ub_rc;}

      void _set_cfg_cfg_ub_spe(Config_UB_spe cfg){_cfg_ub_spe = cfg;}
      Config_UB_spe _get_cfg_ub_spe(){return _cfg_ub_spe;}

      //function to check that common parameters are all filled the same!
      void Check_common_parameters();

    protected:
      Config_COpHit             _cfg_cophit ;
      Config_Deconvolver        _cfg_deconvolver ;
      Config_FlashesBeam        _cfg_flashesbeam ;
      Config_FlashFiltering     _cfg_flashfiltering ;
      Config_Hitfinder_Beam     _cfg_hitfinder_beam ;
      Config_Hitfinder_Cosmic   _cfg_hitfinder_cosmic ;
      Config_Opflash            _cfg_opflash ;
      Config_Saturation_Merger  _cfg_saturation_merger ;
      Config_UB_rc              _cfg_ub_rc ;
      Config_UB_spe             _cfg_ub_spe ;

    };

  }

  #endif
