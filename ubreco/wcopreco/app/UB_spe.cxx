#include "UB_spe.h"




namespace wcopreco {

wcopreco::UB_spe::UB_spe(bool mult_flag, float op_gain, const Config_UB_spe &configSPE)
: kernel_fourier ("UB_SPE",mult_flag), _cfgSPE(configSPE)
 {
   gain = op_gain;
 }

std::vector<double> wcopreco::UB_spe::Get_wfm(int nbins, float tick_width_ns)

  {
    std::vector<double> wfm(nbins,0);
    double value=0;
    int size = wfm.size();
    double X;

    for (int i=0; i<size; i++) {
      X = (double(i)+0.5);
      value = 1./_cfgSPE._spe_p0*pow(X/_cfgSPE._spe_p1,4)*exp(-X/_cfgSPE._spe_p1)*gain;
      wfm[i] = value;
    }

    return wfm;
  }

}
