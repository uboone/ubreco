#include "Config_Hitfinder_Beam.h"

namespace wcopreco {

  Config_Hitfinder_Beam::Config_Hitfinder_Beam()
    {
        _num_channels = 32;
        _nbins_beam = 1500 ;
        _rebin_frac = 6 ; //250 are 1500/6
        _l1_content_thresh = 0.3 ;
        _frac_G_t2_first = 0.75 ;
        _frac_G_sametime = 0.25  ;
        _G_p0 = 6.0 ;
        _G_p1 = 3.0 ;
        _G_p2 = 1.5 ;
        _Lasso_p0 = 5.0 ;
        _Lasso_p1 = 100000;
        _Lasso_p2 = 0.05 ;
        _totPE_v_thresh = 0.2 ;
        _mult_v_thresh  = 1.5 ;
        _l1_mult_v_thresh  = 1.0 ;

        _tick_width_us = .015625;
    }

}
