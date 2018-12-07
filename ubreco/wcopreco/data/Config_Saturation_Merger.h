#ifndef CONFIG_SATURATION_MERGER_H
#define CONFIG_SATURATION_MERGER_H


#include <vector>
#include <iostream>
#include <sstream>

namespace wcopreco {

  class Config_Saturation_Merger {
  public:
    Config_Saturation_Merger();
    ~Config_Saturation_Merger() {};
     int     _num_channels ; //Number of PMT channels
     short   _sat_threshold; // Value at which a tick is marked saturated
     double  _baseline_default; //Default Assumed Baseline
     double  _baseline_difference_max; //Maximum allowed adjustment to new baseline from original one found
     int     _cosmic_tick_window; //How many ticks in half a cosmic waveform
     float   _tick_width_us;  //Tick width in microseconds
     double  _low_bound_baseline_search; //Minimum baseline allowed to be found
     double  _high_bound_baseline_search; //Maximum baseline allowed to be found
     int     _nbins_baseline_search; //Number of bins at the beginning of wfm to search through to find a modal baseline
     int     _nbins_saturation_threshold; //number of bins required to be saturated in order to mark the wfm for replacement if possible
     std::vector<float> _scaling_by_channel; //The scaling from low gain to high gain by channel.

     void _set_num_channels(int n){_num_channels = n;}
     int _get_num_channels(){return _num_channels;}

     void _set_sat_threshold(short s){_sat_threshold = s;}
     short _get_sat_threshold(){return _sat_threshold;}

     void _set_baseline_default(double b){_baseline_default = b;}
     double _get_baseline_default(){return _baseline_default;}

     void _set_baseline_difference_max(double b){_baseline_difference_max = b;}
     double _get_baseline_difference_max(){return _baseline_difference_max;}

     void _set_cosmic_tick_window(int n){_cosmic_tick_window = n;}
     int _get_cosmic_tick_window(){return _cosmic_tick_window;}

     void _set_tick_width_us(float width){_tick_width_us = width;}
     float _get_tick_width_us(){return _tick_width_us;}

     void _set_low_bound_baseline_search(double bound){_low_bound_baseline_search = bound;}
     double _get_low_bound_baseline_search(){return _low_bound_baseline_search;}

     void _set_high_bound_baseline_search(double bound){_high_bound_baseline_search = bound;}
     double _get_high_bound_baseline_search(){return _high_bound_baseline_search;}

     void _set_nbins_baseline_search(int n){_nbins_baseline_search = n;}
     int _get_nbins_baseline_search(){return _nbins_baseline_search;}

     void _set_nbins_saturation_threshold(int n){_nbins_saturation_threshold = n;}
     int _get_nbins_saturation_threshold(){return _nbins_saturation_threshold;}

     void _set_scaling_by_channel(std::vector<float> scalings_v);
     std::vector<float> _get_scaling_by_channel(){return _scaling_by_channel;}

  protected:


  };

}
#endif
