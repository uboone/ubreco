#ifndef CONFIG_COPHIT_H
#define CONFIG_COPHIT_H

#include <vector>
namespace wcopreco{

  class Config_COpHit {
  public:
    Config_COpHit();
    ~Config_COpHit() {};

    //num bins for a cosmic flash
   int    _nbins_cosmic;
   //default baseline
   int    _baseline_default;
   //min integral value for good baseline
   double  _COphit_integral_thresh;
   //max diff between baseline and default baseline for good baseline
   double  _COphit_baseline_diff_thresh;
   // taking into account factor of 2 for 0.6 us window
   double  _pe_factor;
   //Percent uncertainty in baseline
   double  _Baseline_uncertainty;
   //percent uncertainty in a bad baseline
   double  _Baseline_unc_bad_baseline;
   //parameters used for calculating integral for bad ch 28
   double  _cal_integral_p0;
   double  _cal_integral_p1;
   double  _cal_integral_p2;
   double  _cal_integral_p3;
   double  _cal_integral_p4;
   double  _cal_integral_p5;
   std::vector<bool>  _channel_status_v ; //Vector of length number of channels, true is good channels, false is bad channel

   void _set_baseline_default(double b){_baseline_default = b;}
   double _get_baseline_default(){return _baseline_default;}

    void _set_nbins_cosmic(int n){_nbins_cosmic = n;}
    float _get_nbins_cosmic(){return _nbins_cosmic;}

    void _set_COphit_integral_thresh(double value) {_COphit_integral_thresh = value;}
    double _get_COphit_integral_thresh() {return _COphit_integral_thresh;}

    void _set_COphit_baseline_diff_thresh(double value) {_COphit_baseline_diff_thresh = value;}
    double _get_COphit_baseline_diff_thresh() {return _COphit_baseline_diff_thresh;}

    void _set_pe_factor(double value) {_pe_factor = value;}
    double _get_pe_factor() {return _pe_factor;}

    void _set_Baseline_uncertainty(double value) {_Baseline_uncertainty = value;}
    double _get_Baseline_uncertainty() {return _Baseline_uncertainty;}

    void _set_Baseline_unc_bad_baseline(double value) {_Baseline_unc_bad_baseline = value;}
    double _get_Baseline_unc_bad_baseline() {return _Baseline_unc_bad_baseline;}

    void _set_cal_integral_p0(double value) {_cal_integral_p0 = value;}
    double _get_cal_integral_p0() {return _cal_integral_p0;}

    void _set_cal_integral_p1(double value) {_cal_integral_p1 = value;}
    double _get_cal_integral_p1() {return _cal_integral_p1;}

    void _set_cal_integral_p2(double value) {_cal_integral_p2 = value;}
    double _get_cal_integral_p2() {return _cal_integral_p2;}

    void _set_cal_integral_p3(double value) {_cal_integral_p3 = value;}
    double _get_cal_integral_p3() {return _cal_integral_p3;}

    void _set_cal_integral_p4(double value) {_cal_integral_p4 = value;}
    double _get_cal_integral_p4() {return _cal_integral_p4;}

    void _set_cal_integral_p5(double value) {_cal_integral_p5 = value;}
    double _get_cal_integral_p5() {return _cal_integral_p5;}



  protected:

  };

}

#endif
