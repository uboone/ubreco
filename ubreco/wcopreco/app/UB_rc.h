#ifndef UB_RC_H
#define UB_RC_H

#include "kernel_fourier.h"
#include "Config_UB_rc.h"

#include <vector>
#include <iostream>
#include <sstream>
#include <cmath>
#include <string>

namespace wcopreco {

  class UB_rc : public kernel_fourier {
  public:
     UB_rc(bool mult_flag, bool bad_ch, const Config_UB_rc &configRC);
     virtual ~UB_rc() {};

    std::vector<double> Get_wfm(int nbins, float tick_width_ns);
    //std::vector<double> Get_pow_spec(int nbins, float tick_width_ns, std::vector<double>* mag, std::vector<double>* phase){};

    std::string name;
    // 0 = divide 1 = multiply
    bool bad_ch;


  protected:

    Config_UB_rc _cfgRC;

  };

}
#endif
