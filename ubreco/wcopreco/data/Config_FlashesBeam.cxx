#include "Config_FlashesBeam.h"

namespace wcopreco{
  Config_FlashesBeam::Config_FlashesBeam(){
    //defaults
    _rebin_frac = 6;
    _nbins_beam = 1500;
     _num_channels = 32;
      _bflash_pe_thresh  = 6.0 ;
     _bflash_mult_thresh = 3.0 ;
     _bflash_bin_diff_p0 = 78 ;  
     _bflash_bin_diff_p1 = 4  ;
     _bflash_bin_diff_p2  = 15 ;
     _KS_test_thresh  = 0.1 ;
     _bflash_bin_start_cushion = 2  ;

  }
}
