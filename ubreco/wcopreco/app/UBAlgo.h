#ifndef UBALGO_H
#define UBALGO_H

//WCOpReco includes
#include "Config_Params.h"
#include "OpWaveform.h"
#include "OpWaveformCollection.h"
#include "EventOpWaveforms.h"
#include "Deconvolver.h"
#include "HitFinder_cosmic.h"
#include "Flashes_cosmic.h"
#include "Flashes_beam.h"
#include "HitFinder_beam.h"
#include "FlashFiltering.h"

//Config includes
#include "Config_Params.h"

//UB specific classes
#include "UBEventWaveform.h"
#include "DataReader.h"
#include "Saturation_Merger.h"
#include "UB_rc.h"
#include "UB_spe.h"

//root includes
#include "TH1S.h"
#include "TF1.h"
#include "TVirtualFFT.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"

//c++ includes
#include <iostream>
#include <sstream>
#include <time.h>

namespace wcopreco  {

  //This is a class to run the WCOpReco code, for uboone
  class UBAlgo {
  public:
    //UBAlgo(const Config_Params &cfg_all);
    UBAlgo();
    ~UBAlgo();

    void Configure(const Config_Params &cfg_all);
    void SaturationCorrection(UBEventWaveform *_UB_Ev_wfm);
    void set_merged_beam(OpWaveformCollection &ext_merged_beam);
    void set_merged_cosmic(OpWaveformCollection &ext_merged_cosmic);
    void Run(std::vector<float> * op_gain,
	     std::vector<float> * op_gainerror,
	     std::vector<wcopreco::kernel_fourier_container> * kernel_container_v);

    OpflashSelection get_flashes_cosmic(){return flashes_cosmic;};
    OpflashSelection get_flashes_beam(){return flashes_beam;};
    OpflashSelection get_flashes(){return flashes;};
    OpWaveformCollection get_merged_beam(){return merged_beam;};
    OpWaveformCollection get_merged_cosmic(){return merged_cosmic;};
    void clear_flashes();

  protected:
    Config_Params _cfg;
    OpflashSelection flashes_cosmic;
    OpflashSelection flashes_beam;
    OpflashSelection flashes;
    OpWaveformCollection merged_beam;
    OpWaveformCollection merged_cosmic;

  };

}
#endif
