#include "Config_Saturation_Merger.h"

namespace wcopreco {

  Config_Saturation_Merger::Config_Saturation_Merger()
    {
       _num_channels = 32 ;
       _sat_threshold = 4080 ;
       _baseline_default =2050;
       _baseline_difference_max = 8;
       _cosmic_tick_window = 20;
       _tick_width_us = .015625;
       _low_bound_baseline_search = 1500;
       _high_bound_baseline_search = 2500;
       _nbins_baseline_search = 20;
       _nbins_saturation_threshold = 3;
       _scaling_by_channel =  {10.13, 10.20, 10.13, 10.05, 9.96, 9.95, 10.04, 9.58, 9.42, 9.81, 9.25, 9.61, 9.56, 9.35, 9.99, 9.66, 10.26, 9.82, 10.32, 10.08, 9.77, 9.64, 10.14, 9.74, 9.76, 10.10, 10.48, 9.81, 9.99, 9.79, 9.93, 10.01};
     }

     void Config_Saturation_Merger::_set_scaling_by_channel(std::vector<float> scalings_v) {
       _scaling_by_channel = scalings_v;
       if (int(scalings_v.size()) != _num_channels) {
         std::cout << "Careful, size of scaling vector is not the same as the number of channels currently set (try to set _num_channels first if you think they should be the same).\n";
         std::cout << _num_channels << " Is the number of channels, " << scalings_v.size() << "  Is the size of scalings vector set\n";
       }
     }

}
