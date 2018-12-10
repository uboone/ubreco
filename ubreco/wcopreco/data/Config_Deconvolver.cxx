#include "Config_Deconvolver.h"

namespace wcopreco {

  Config_Deconvolver::Config_Deconvolver()
    {
        _num_channels = 32; 
        _nbins_beam = 1500;
        _baseline_difference_max = 8;
        _tick_width_us = .015625;
        _high_freq_p0 = 0.45;
        _high_freq_p1 = 3.07;
        _latelight_filter_p0 = 0.05;
        _latelight_filter_p1 = 0.45;
        _latelight_filter_p2 = 3.07;
        _baseline_safety_subtraction = 100;
        _xq  = 0.5 ;
        _xq_diff = 0.34 ;
        _n_bins_end_wfm = 4 ;
        _small_content_bump = 0.01;
        _nbins_baseline_search = 20;

    }

}
