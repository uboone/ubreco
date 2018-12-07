#ifndef CONFIG_HITFINDER_BEAM_H
#define CONFIG_HITFINDER_BEAM_H

namespace wcopreco {
  class Config_Hitfinder_Beam {
  public:
    Config_Hitfinder_Beam();
    ~Config_Hitfinder_Beam() {};
     int     _num_channels; //Number of pmt channels providing waveforms
     int     _nbins_beam; //Original number of bins in a beam waveform
     int     _rebin_frac; //Fraction to rebin beam wfms into. By default 1500bin /6frac = 250 bin final waveform size
     double  _l1_content_thresh; // Minimum needed in bin to add a value to the l1 fitter
     double  _frac_G_t2_first; //Fraction expected that t2 is before t1 in l1 fitter
     double  _frac_G_sametime; //Fraction expected that t2 ==t1 in fitter
     double  _G_p0; //Parameter in G of L1 fitter
     double  _G_p1; //Parameter in G of L1 fitter
     double  _G_p2 ; //Parameter in G of L1 fitter
     double  _Lasso_p0; //Lasso Model Parameter
     int     _Lasso_p1; //Lasso Model Parameter
     double  _Lasso_p2; //Lasso Model Parameter
     double  _totPE_v_thresh; //Minimum value in rebinned content allowed to add to total PE
     double  _mult_v_thresh; //Minimum value in rebinned content allowed in order to add +1 to multiplicity
     double  _l1_mult_v_thresh; //Hitfinder beam

     double _tick_width_us; //Width of original bin in microseconds

     void _set_num_channels(int n){_num_channels = n;}
     int _get_num_channels(){return _num_channels;}

     void _set_nbins_beam(int n){_nbins_beam = n;}
     int _get_nbins_beam(){return _nbins_beam;}

     void _set_rebin_frac(int n){_rebin_frac = n;}
     int _get_rebin_frac(){return _rebin_frac;}

     void _set_l1_content_thresh(double value) {_l1_content_thresh = value;}
     double _get_l1_content_thresh() {return _l1_content_thresh;}

     void _set_frac_G_t2_first(double value) {_frac_G_t2_first = value;}
     double _get_frac_G_t2_first() {return _frac_G_t2_first;}

     void _set_frac_G_sametime(double value) {_frac_G_sametime = value;}
     double _get_frac_G_sametime() {return _frac_G_sametime;}

     void _set_G_p0(double value) {_G_p0 = value;}
     double _get_G_p0() {return _G_p0;}

     void _set_G_p1(double value) {_G_p1 = value;}
     double _get_G_p1() {return _G_p1;}

     void _set_G_p2(double value) {_G_p2 = value;}
     double _get_G_p2() {return _G_p2;}

     void _set_Lasso_p0(double value) {_Lasso_p0 = value;}
     double _get_Lasso_p0() {return _Lasso_p0;}

     void _set_Lasso_p1(int value) {_Lasso_p1 = value;}
     int _get_Lasso_p1() {return _Lasso_p1;}

     void _set_Lasso_p2(double value) {_Lasso_p2 = value;}
     double _get_Lasso_p2() {return _Lasso_p2;}

     void _set_totPE_v_thresh(double value) {_totPE_v_thresh = value;}
     double _get_totPE_v_thresh() {return _totPE_v_thresh;}

     void _set_mult_v_thresh(double value) {_mult_v_thresh = value;}
     double _get_mult_v_thresh() {return _mult_v_thresh;}

     void _set_l1_mult_v_thresh(double value) {_l1_mult_v_thresh = value;}
     double _get_l1_mult_v_thresh() {return _l1_mult_v_thresh;}

     void _set_tick_width_us(double value) {_tick_width_us = value;}
     double _get_tick_width_us() {return _tick_width_us;}


  protected:


  };

}
#endif
