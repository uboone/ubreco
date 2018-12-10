#include "COphit.h"

using namespace wcopreco;

wcopreco::COphit::COphit(int ch_no, OpWaveform *wfm, double time, double gain, double gain_err , const Config_COpHit &config)
  : channel_no(ch_no)
  , time(time)
  , gain(gain)
  , gain_err(gain_err)
  , _cfg(config)
{
  // calculate baseline
  baseline = wfm->at(0);

  // calculate peak and integral
  peak = 0;
  integral = 0;

  for (int i=0; i!=_cfg._nbins_cosmic; i++){
    double content = wfm->at(i) - baseline;
    if (content > peak){
      peak  = content;
    }
    integral += content;
  }

  // calculate PE and its error ...
  good_baseline = false;
  if (integral > _cfg._COphit_integral_thresh || fabs(baseline-_cfg._baseline_default)<_cfg._COphit_baseline_diff_thresh) good_baseline = true;

  // special treatment of FEM channel 28 ...
  if (_cfg._channel_status_v.at(channel_no) == false) integral = cal_integral(peak);
  if (!good_baseline) integral = cal_integral(peak);

  if (good_baseline){
    PE = integral / gain * _cfg._pe_factor;
    // take care of PE ...
    if (PE < 0) PE = 0;

    PE_err = sqrt(pow(PE * gain_err/gain ,2) // gain uncertainties
		  + 2 * PE // statistical term ...
		  + pow(_cfg._nbins_cosmic/sqrt(3.)/gain * 2,2)  // basline (1/sqrt(3.) ADC with 40 time tics ...)
		  + pow(PE*_cfg._Baseline_uncertainty,2) // 3% relative uncertainties (Baseline guess) ...
		  );

  }else{
    PE = integral/gain * 2;
    if (PE < 0) PE = 0;
    PE_err = _cfg._Baseline_unc_bad_baseline* PE;
  }

}

wcopreco::COphit::~COphit(){
}

double wcopreco::COphit::cal_integral(double peak){
  double content;
  if (peak <=_cfg._cal_integral_p0){
    content = _cfg._cal_integral_p1 * peak + _cfg._cal_integral_p2 * pow(peak,2);
  }else{
    content = _cfg._cal_integral_p3 + _cfg._cal_integral_p4 * (peak - _cfg._cal_integral_p0) + _cfg._cal_integral_p5 * pow(peak-_cfg._cal_integral_p0,2);
  }
  return content;
}
