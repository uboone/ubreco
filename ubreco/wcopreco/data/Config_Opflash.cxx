#include "Config_Opflash.h"

namespace wcopreco{
  Config_Opflash::Config_Opflash(){
    //defaults
     _tick_width_us = .015625;
     _rebin_frac = 6;
     _num_channels = 32;
     _PE_err_cosmic = 6.4;
     _PE_subtract = 0.15;
     _flash_low_time_cushion  = -3;
     _flash_high_time_cushion = +37;
     _PE_err_beam = 0.2;
     _mult_content_thresh = 1.5;
     _mult_required = 3;
     _PE_err_stat_beam = 1.875;
     _PE_err_unc_beam = 0.02;
     _addl1_pe_thresh = 10.0;
     _addl1_mult_thresh = 3.0;

  }
}
