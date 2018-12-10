#ifndef CONFIG_FLASHESBEAM_H
#define CONFIG_FLASHESBEAM_H

namespace wcopreco{

  class Config_FlashesBeam {
  public:
    Config_FlashesBeam();
    ~Config_FlashesBeam() {};

    //num beam bins/ rebinsize
    float   _rebin_frac;
    //number of bins in a beam waveform
    int     _nbins_beam;
    //number of _num_channels
    int     _num_channels;
    //min pe threshold for beam flashes
     double  _bflash_pe_thresh;
     //min mult threshold for beam flashes
     double  _bflash_mult_thresh;
     // set of parameters used to threshold the time diff between beam flashes
     int     _bflash_bin_diff_p0;
     int     _bflash_bin_diff_p1;
     int     _bflash_bin_diff_p2;
     //min ks_test value
     double  _KS_test_thresh;
     //start the flash this many bins in
     int     _bflash_bin_start_cushion;


    void _set_rebin_frac(float frac){_rebin_frac = frac;}
    float _get_rebin_frac(){return _rebin_frac;}

    void _set_nbins_beam(int n){_nbins_beam = n;}
    float _get_nbins_beam(){return _nbins_beam;}

    void _set_num_channels(int n){_num_channels = n;}
    int _get_num_channels(){return _num_channels;}

    void _set_bflash_pe_thresh(double value) {_bflash_pe_thresh = value;}
      double _get_bflash_pe_thresh() {return _bflash_pe_thresh;}

      void _set_bflash_mult_thresh(double value) {_bflash_mult_thresh = value;}
      double _get_bflash_mult_thresh() {return _bflash_mult_thresh;}

      void _set_bflash_bin_diff_p0(int value) {_bflash_bin_diff_p0 = value;}
      int _get_bflash_bin_diff_p0() {return _bflash_bin_diff_p0;}

      void _set_bflash_bin_diff_p1(int value) {_bflash_bin_diff_p1 = value;}
      int _get_bflash_bin_diff_p1() {return _bflash_bin_diff_p1;}

      void _set_bflash_bin_diff_p2(int value) {_bflash_bin_diff_p2 = value;}
      int _get_bflash_bin_diff_p2() {return _bflash_bin_diff_p2;}

      void _set_KS_test_thresh(double value) {_KS_test_thresh = value;}
      double _get_KS_test_thresh() {return _KS_test_thresh;}

      void _set_bflash_bin_start_cushion(int value) {_bflash_bin_start_cushion = value;}
      int _get_bflash_bin_start_cushion() {return _bflash_bin_start_cushion;}


  protected:

  };

}

#endif
