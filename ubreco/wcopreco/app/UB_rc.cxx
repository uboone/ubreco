#include "UB_rc.h"

using namespace wcopreco;

 UB_rc::UB_rc(bool mult_flag, bool is_bad, const Config_UB_rc &configRC)
 : kernel_fourier ("UB_RC",mult_flag), _cfgRC(configRC)
  {
    bad_ch = is_bad;
  }

std::vector<double> wcopreco::UB_rc::Get_wfm(int nbins, float tick_width_ns)
{
  std::vector<double> wfm(nbins,0);
  int size = wfm.size();
  double rc_tau = 0;
  if (bad_ch){
    rc_tau = _cfgRC._rc_tau_badch;
  }
  else {
    rc_tau = _cfgRC._rc_tau_goodch;
  }

  double X;
  for (int i=0; i<size; i++) {
    X = i+.5;

    double content = -1./rc_tau * exp(-X/rc_tau);
    if (i==0) content += 1;
    wfm[i] = content;
  }

  return wfm;

}
