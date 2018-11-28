#ifndef DECONVOLVER_H
#define DECONVOLVER_H

//deconv functions
#include "kernel_fourier.h"
#include "kernel_fourier_container.h"
#include "LassoModel.h"
#include "ElasticNetModel.h"
#include "LinearModel.h"
#include "Opflash.h"
#include "COphit.h"

//data
#include "OpWaveform.h"
#include "OpWaveformCollection.h"
#include "EventOpWaveforms.h"
#include "Config_Deconvolver.h"


//c++ includes
#include <vector>
#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <Eigen/Dense>

//root
#include "TCanvas.h"
#include "TLine.h"

namespace wcopreco{


  class Deconvolver {
  public:
    Deconvolver(OpWaveformCollection &merged_beam, bool with_filters, std::vector<kernel_fourier_container> &input_k_container_v, const Config_Deconvolver &); //OpWaveform op_wfm, kernel_fourier_shape kernel_fourier, noise_remover noise, std::vector<short???> LL_shape
    ~Deconvolver() {};


    const std::vector<kernel_fourier_container> & get_kernel_container_v() {return *kernel_container_v;}
    kernel_fourier_container get_kernel_container_entry(int channel) {return kernel_container_v->at(channel);}
    // void add_kernel_container_entry(kernel_fourier *kernel, int channel =-1);
    // void clear_kernels();
    // void add_kernel_container


    void set_filter_status(bool status) {filter_status = status;}
    OpWaveformCollection Deconvolve_Collection(OpWaveformCollection & merged_beam);
    double HighFreqFilter(double frequency);
    double LateLightFilter(double frequency2);
    void Remove_Baseline_Leading_Edge(OpWaveform &wfm);
    void Remove_Baseline_Secondary(OpWaveform &wfm);
    OpWaveform Deconvolve_One_Wfm(OpWaveform &wfm, const kernel_fourier_container &kernel_container);
    std::pair<double,double> cal_mean_rms(std::vector<double> wfm, int nbin);




    double KS_maxdiff(int n, double *array1, double *array2);

    /*
    The KS_maxdiff function takes two array PDFs of n elements, and transforms
    them into CDFs then it finds the maximum difference between the CDFs and
    returns it.
    */
    //UB_rc Make_UB_rc(int ch);
  protected:
    Config_Deconvolver _cfg;
    int nbins;

    bool filter_status;

    std::vector<float>  op_gain;
    const std::vector<kernel_fourier_container> *kernel_container_v;
    OpWaveformCollection deconvolved_collection;
  };

}
#endif
