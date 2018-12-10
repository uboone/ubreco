#ifndef KERNEL_FOURIER_H
#define KERNEL_FOURIER_H

#include "TH1D.h"
#include "TVirtualFFT.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"

#include <vector>
#include <string>
#include <iostream>
#include <sstream>

namespace wcopreco {

  class kernel_fourier {
  public:
    kernel_fourier(std::string word, bool flag);
    virtual ~kernel_fourier() {};

    virtual std::vector<double> Get_wfm(int nbins, float tick_width_ns) = 0;
      /*
      This function is not declared, and cannot be used without a coded childclass
      version. It should in return a vector of doubles with size
      nbins that represents the waveform.
      */
    virtual void Get_pow_spec(int nbins, float tick_width_ns, std::vector<double>* mag, std::vector<double>* phase);
      /*
      Get_pow_spec is a function that takes in the number of bins, and their tick width in ns, as well as
      two pass by reference vectors of doubles. It then calls Get_wfm, which is child class specific,
      and performs a fourier transform filling the vector<double>s with magnitude and phase.
      */
    std::string name;

    bool mult_div;
    // A boolean representing whether the kernel is multiplied or divided in the
    // deconvolution stage.
    // 0 = divide 1 = multiply
  protected:



  };

}

#endif
