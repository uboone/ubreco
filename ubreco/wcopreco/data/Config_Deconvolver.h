#ifndef CONFIG_DECONVOLVER_H
#define CONFIG_DECONVOLVER_H

namespace wcopreco {
  class Config_Deconvolver {
  public:
    Config_Deconvolver();
    ~Config_Deconvolver() {};
   int     _num_channels; //Number of channels reading out waveforms
   int     _nbins_beam ; //Nbins in a beam waveform (waveform being deconvolved)
   double  _baseline_difference_max; //Maximum allowed adjustment to new baseline from original one found
   float   _tick_width_us; //Tick width in microseconds
   double  _high_freq_p0;  // HighFreqFilter parameter (see function in Deconvolver.cxx)
   double  _high_freq_p1;  // HighFreqFilter parameter (see function in Deconvolver.cxx)
   double  _latelight_filter_p0;  // LateLightFilter parameter (see function in Deconvolver.cxx)
   double  _latelight_filter_p1;  // LateLightFilter parameter (see function in Deconvolver.cxx)
   double  _latelight_filter_p2;  // LateLightFilter parameter (see function in Deconvolver.cxx)
   double  _baseline_safety_subtraction; // deconv Remove_Baseline_Secondary
   double  _xq; //Middle Quartile Fitting parameter
   double  _xq_diff; //Quartile Fitting upper and lower dist from middle parameter
   int     _n_bins_end_wfm; //Number of bins at end of deconvolved wfm set to zero
   double  _small_content_bump; //Small offset amount added to content in wfm after deconvolution.     -
   int     _nbins_baseline_search; //Number of bins to search in order to determine modal baseline     -

   void _set_num_channels(int n){_num_channels = n;}
   int _get_num_channels(){return _num_channels;}

   void _set_nbins_beam(int n){_nbins_beam = n;}
   int _get_nbins_beam(){return _nbins_beam;}

   void _set_baseline_difference_max(double b){_baseline_difference_max = b;}
   double _get_baseline_difference_max(){return _baseline_difference_max;}

   void _set_tick_width_us(float width){_tick_width_us = width;}
   float _get_tick_width_us(){return _tick_width_us;}

   void _set_high_freq_p0(double p0){_high_freq_p0 = p0;}
   double _get_high_freq_p0(){return _high_freq_p0;}

   void _set_high_freq_p1(double p1){_high_freq_p1 = p1;}
   double _get_high_freq_p1(){return _high_freq_p1;}

   void _set_latelight_filter_p0(double p0){_latelight_filter_p0 = p0;}
   double _get_latelight_filter_p0(){return _latelight_filter_p0;}

   void _set_latelight_filter_p1(double p1){_latelight_filter_p1 = p1;}
   double _get_latelight_filter_p1(){return _latelight_filter_p1;}

   void _set_latelight_filter_p2(double p2){_latelight_filter_p2 = p2;}
   double _get_latelight_filter_p2(){return _latelight_filter_p2;}

   void _set_baseline_safety_subtraction(double adjust){_baseline_safety_subtraction = adjust;}
   double _get_baseline_safety_subtraction(){return _baseline_safety_subtraction;}

   void _set_xq(double value) {_xq = value;}
   double _get_xq() {return _xq;}

   void _set_xq_diff(double value) {_xq_diff = value;}
   double _get_xq_diff() {return _xq_diff;}

   void _set_n_bins_end_wfm(int n){_n_bins_end_wfm = n;}
   int _get_n_bins_end_wfm(){return _n_bins_end_wfm;}

   void _set_small_content_bump(double n){_small_content_bump = n;}
   double _get_small_content_bump(){return _small_content_bump;}

   void _set_nbins_baseline_search(int n){_nbins_baseline_search = n;}
   int _get_nbins_baseline_search(){return _nbins_baseline_search;}


  protected:


  };

}
#endif
