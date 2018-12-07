#include "Config_COpHit.h"

namespace wcopreco{
  Config_COpHit::Config_COpHit(){
    //defaults
    _nbins_cosmic = 40;
    _COphit_integral_thresh = 20e3 ;
    _COphit_baseline_diff_thresh = 50 ;
    _pe_factor  = 2.0;
    _Baseline_uncertainty  = 0.03 ;
    _Baseline_unc_bad_baseline  = 2.0 ;
    _cal_integral_p0 = 4000;
    _cal_integral_p1 = 1.06241e+01;
    _cal_integral_p2 = 2.01214e-04;
    _cal_integral_p3 = 4.5715824e4;
    _cal_integral_p4 = 8.62296e+00;
    _cal_integral_p5 = 6.76898e-04;
    _baseline_default =2050;
    _channel_status_v = {true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true};
    _channel_status_v[28] = false ;

  }
}
