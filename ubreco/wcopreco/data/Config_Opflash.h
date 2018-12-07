#ifndef CONFIG_OPFLASH_H
#define CONFIG_OPFLASH_H

namespace wcopreco{

  class Config_Opflash {
  public:
    Config_Opflash();
    ~Config_Opflash() {};

    //bin width should be derived from detector constants when in larlite/larsoft
    float   _tick_width_us;
    //num beam bins/ rebinsize
    int   _rebin_frac;
    //num of channels in the detector
    int     _num_channels;
    //  pe error for cosmic, default: 11/sqrt(3.)
    double  _PE_err_cosmic;
    //// 250 kHz at 0.6 us - subtract from each bin
    double  _PE_subtract;
    // used in setting the low time
    int     _flash_low_time_cushion;
    //used in setting the high time
    int     _flash_high_time_cushion;
    //base pe error for beam flashes
    double  _PE_err_beam;
    //mult threshold to add to mult
    double  _mult_content_thresh;
    //mult required to save beam
    int     _mult_required;
    //stat error in beam flashes :  7.5 us * random noise ...
    double  _PE_err_stat_beam;
    //used in calcluating pe error for beam
    double  _PE_err_unc_beam;
    //threshold to add l1 pe
    double  _addl1_pe_thresh;
    //threshold to add l1 mult
    double  _addl1_mult_thresh;

    void _set_tick_width_us(float width){_tick_width_us = width;}
    float _get_tick_width_us(){return _tick_width_us;}

    void _set_rebin_frac(float frac){_rebin_frac = frac;}
    float _get_rebin_frac(){return _rebin_frac;}

    void _set_num_channels(int n){_num_channels = n;}
    int _get_num_channels(){return _num_channels;}

    void _set_PE_err_cosmic(double value) {_PE_err_cosmic = value;}
    double _get_PE_err_cosmic() {return _PE_err_cosmic;}

    void _set_PE_subtract(double value) {_PE_subtract = value;}
    double _get_PE_subtract() {return _PE_subtract;}

    void _set_flash_low_time_cushion(int value) {_flash_low_time_cushion = value;}
    int _get_flash_low_time_cushion() {return _flash_low_time_cushion;}

    void _set_flash_high_time_cushion(int value) {_flash_high_time_cushion = value;}
    int _get_flash_high_time_cushion() {return _flash_high_time_cushion;}

    void _set_PE_err_beam(double value) {_PE_err_beam = value;}
    double _get_PE_err_beam() {return _PE_err_beam;}

    void _set_mult_content_thresh(double value) {_mult_content_thresh = value;}
    double _get_mult_content_thresh() {return _mult_content_thresh;}

    void _set_mult_required(int value) {_mult_required = value;}
    int _get_mult_required() {return _mult_required;}

    void _set_PE_err_stat_beam(double value) {_PE_err_stat_beam = value;}
    double _get_PE_err_stat_beam() {return _PE_err_stat_beam;}

    void _set_PE_err_unc_beam(double value) {_PE_err_unc_beam = value;}
    double _get_PE_err_unc_beam() {return _PE_err_unc_beam;}

    void _set_addl1_pe_thresh(double value) {_addl1_pe_thresh = value;}
    double _get_addl1_pe_thresh() {return _addl1_pe_thresh;}

    void _set_addl1_mult_thresh(double value) {_addl1_mult_thresh = value;}
    double _get_addl1_mult_thresh() {return _addl1_mult_thresh;}
  protected:

  };

}

#endif
