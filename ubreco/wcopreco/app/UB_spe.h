#ifndef UB_SPE_H
#define UB_SPE_H

#include "kernel_fourier.h"
#include "Config_UB_spe.h"

#include "TH1D.h"
#include "TVirtualFFT.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"

#include <vector>
#include <iostream>
#include <sstream>
#include <cmath>

namespace wcopreco {

  class UB_spe : public kernel_fourier {
  public:
     UB_spe(bool mult_flag, float gain, const Config_UB_spe &configSPE);
     virtual ~UB_spe() {};

    std::vector<double> Get_wfm(int nbins, float tick_width_ns);
    // void Get_pow_spec(int nbins, float tick_width_ns, std::vector<double>* mag, std::vector<double>* phase);
    float gain;

  protected:
    Config_UB_spe _cfgSPE;



  };

}

#endif
