#include "Config_FlashFiltering.h"

namespace wcopreco{
  Config_FlashFiltering::Config_FlashFiltering(){
    //defaults
    _flash_filter_time_thresh = 2.4;
    _flash_filter_pe_thresh = 0.7;
    _do_swap_channels  = true;
  }
}
