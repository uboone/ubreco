#include "OpWaveform.h"
#include "OpWaveformCollection.h"
#include "EventOpWaveforms.h"
#include "UBEventWaveform.h"
#include "DataReader.h"
#include "Deconvolver.h"
#include "Saturation_Merger.h"
#include "HitFinder_cosmic.h"
#include "Flashes_cosmic.h"
#include "Flashes_beam.h"
#include "HitFinder_beam.h"
#include "FlashFiltering.h"
#include "UBAlgo.h"
#include "Config_Params.h"

//root includes
#include "TH1S.h"
#include "TF1.h"
#include "TVirtualFFT.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"

#include <iostream>
#include <sstream>
#include <time.h>

//Construct the vector of kernel containers (one container per channel)
std::vector<wcopreco::kernel_fourier_container> Build_UB_kernels(wcopreco::Config_Params cfg_all, std::vector<float> op_gain){
  std::vector<wcopreco::kernel_fourier_container> kernel_container_v;
  kernel_container_v.resize(cfg_all._get_cfg_deconvolver()._get_num_channels());
  for (int i =0 ; i<cfg_all._get_cfg_deconvolver()._get_num_channels(); i++){

    wcopreco::UB_spe *spe = new wcopreco::UB_spe(true, op_gain.at(i), cfg_all._get_cfg_ub_spe()); //Place UB_spe on heap, so object not deleted
    kernel_container_v.at(i).add_kernel(spe);
    if ( false == cfg_all._get_cfg_cophit()._channel_status_v.at(i) ){
      wcopreco::UB_rc *rc_bad_ch = new wcopreco::UB_rc(true, true, cfg_all._get_cfg_ub_rc());
      kernel_container_v.at(i).add_kernel(rc_bad_ch);
    }
    else{
      wcopreco::UB_rc *rc_good_ch = new wcopreco::UB_rc(true, false, cfg_all._get_cfg_ub_rc());
      kernel_container_v.at(i).add_kernel(rc_good_ch);
    }

  }

  return kernel_container_v;
}

int main(){
  //clock_t t;
  //t=clock();

  //Set the filepath
  std::string file = "./celltreeMC_modified.root";
  std::cout << "\n\nFilepath is set to:   " << file << std::endl;

  //Open the reader, choose event number, create the UBEventWaveform _UB_Ev_wfm
  wcopreco::DataReader reader(&file);
  wcopreco::UBEventWaveform _UB_Ev_wfm;

  wcopreco::Config_Params cfg_all;
  cfg_all.Check_common_parameters();

  //make root file of outputs to test
  TFile output("OurOutput.root", "RECREATE");
  //create TTree
  TTree *OpReco = new TTree("OpReco", "A tree to hold outputs from OpReco");
  //make objects to save
  std::vector<std::vector<double>> Tmerged_beam;
  std::vector<std::vector<double>> Tmerged_cosmic;
  std::vector<double> TtotalPE_all;
  std::vector<double> Tflashtime_all;
  std::vector<double> TtotalPE_cosmic;
  std::vector<double> Tflashtime_cosmic;
  std::vector<double> TtotalPE_beam;
  std::vector<double> Tflashtime_beam;
  std::vector<std::vector<double>> TPEperPMT_all;
  std::vector<std::vector<double>> TPEperPMT_cosmic;
  std::vector<std::vector<double>> TPEperPMT_beam;
  //create branches
  OpReco->Branch("merged_beam", &Tmerged_beam);
  OpReco->Branch("merged_cosmic", &Tmerged_cosmic);
  OpReco->Branch("totalPE_all", &TtotalPE_all);
  OpReco->Branch("flashtime_all", &Tflashtime_all);
  OpReco->Branch("totalPE_cosmic", &TtotalPE_cosmic);
  OpReco->Branch("flashtime_cosmic", &Tflashtime_cosmic);
  OpReco->Branch("totalPE_beam", &TtotalPE_beam);
  OpReco->Branch("flashtime_beam", &Tflashtime_beam);
  OpReco->Branch("PEperPMT_all", &TPEperPMT_all);
  OpReco->Branch("PEperPMT_cosmic", &TPEperPMT_cosmic);
  OpReco->Branch("PEperPMT_beam", &TPEperPMT_beam);
  //loop over all events
  for (int EVENT_NUM = 0; EVENT_NUM < 3; EVENT_NUM ++){
    wcopreco::UBAlgo opreco_run;
    opreco_run.Configure(cfg_all);
    // create UBEventWaveform object
    _UB_Ev_wfm = reader.Reader(EVENT_NUM);
    std::vector<float> op_gain = _UB_Ev_wfm.get_op_gain();
    std::vector<float> op_gainerror = _UB_Ev_wfm.get_op_gainerror();
    //create vector of kernels
    std::vector<wcopreco::kernel_fourier_container> kernel_container_v = Build_UB_kernels(cfg_all, op_gain);
    //saturation correction
    opreco_run.SaturationCorrection(&_UB_Ev_wfm);
    //deconvole beam, find hits, make flashes
    opreco_run.Run(&op_gain, &op_gainerror, &kernel_container_v);
    //get outputs
    wcopreco::OpWaveformCollection merged_beam = opreco_run.get_merged_beam();
    wcopreco::OpWaveformCollection merged_cosmic = opreco_run.get_merged_cosmic();
    wcopreco::OpflashSelection flashes = opreco_run.get_flashes();
    wcopreco::OpflashSelection flashes_cosmic = opreco_run.get_flashes_cosmic();
    wcopreco::OpflashSelection flashes_beam = opreco_run.get_flashes_beam();
    // check size of outputs
    //std::cout << flashes.size() << " :Number of flashes\n";
    //std::cout << flashes_cosmic.size() << " :Number of cosmic flashes\n";
    //std::cout << flashes_beam.size() << " :Number of beam flashes\n";
    //std::cout << flashes_beam[0]->get_total_PE() << " :beam flash[0] totPE\n";
    // fill root branches
    int num_channels = cfg_all._get_cfg_deconvolver()._get_num_channels();
    int num_flashes = flashes.size();
    //int num_flashes_cosmic = flashes_cosmic.size();
    //int num_flashes_beam = flashes_beam.size();
    /*    
    Tmerged_beam.resize(num_channels);
    for (int i = 0; i<num_channels; i++){
      Tmerged_beam.at(i).resize(merged_beam.at(i).size());
      for (int j = 0; j<merged_beam.at(i).size(); j++){
	Tmerged_beam.at(i).at(j) = merged_beam.at(i).at(j);
      }
    }
    Tmerged_cosmic.resize(num_channels);
    for (int i = 0; i<num_channels; i++){
      Tmerged_cosmic.at(i).resize(merged_cosmic.at(i).size());
      for (int j = 0; j<merged_cosmic.at(i).size(); j++){
	Tmerged_cosmic.at(i).at(j) = merged_cosmic.at(i).at(j);
      }
    }
    */
    TtotalPE_all.resize(num_flashes);
    Tflashtime_all.resize(num_flashes);
    for (int i =0; i<num_flashes; i++){
      TtotalPE_all.at(i) = flashes.at(i)->get_total_PE();
      Tflashtime_all.at(i) = flashes.at(i)->get_time();
    }
    /*
    TtotalPE_cosmic.resize(num_flashes_cosmic);
    Tflashtime_cosmic.resize(num_flashes_cosmic);
    for (int i =0; i<num_flashes_cosmic; i++){
      TtotalPE_cosmic.at(i) = flashes_cosmic.at(i)->get_total_PE();
      Tflashtime_cosmic.at(i) = flashes_cosmic.at(i)->get_time();
    }
    TtotalPE_beam.resize(num_flashes_beam);
    Tflashtime_beam.resize(num_flashes_beam);
    for (int i =0; i<num_flashes_beam; i++){
      TtotalPE_beam.at(i) = flashes_beam.at(i)->get_total_PE();
      Tflashtime_beam.at(i) = flashes_beam.at(i)->get_time();
    }
    */
    TPEperPMT_all.resize(num_flashes);
    for (int i = 0; i < num_flashes; i++){
      TPEperPMT_all.at(i).resize(num_channels);
      for (int j = 0; j < num_channels; j++){
	TPEperPMT_all.at(i).at(j) = flashes.at(i)->get_PE(j);
      }
    }
    /*
    TPEperPMT_cosmic.resize(num_flashes_cosmic);
    for (int i = 0; i < num_flashes_cosmic; i++){
      TPEperPMT_cosmic.at(i).resize(num_channels);
      for (int j = 0; j < num_channels; j++){
	TPEperPMT_cosmic.at(i).at(j) = flashes_cosmic.at(i)->get_PE(j);
      }
    }
    TPEperPMT_beam.resize(num_flashes_beam);
    for (int i = 0; i < num_flashes_beam; i++){
      TPEperPMT_beam.at(i).resize(num_channels);
      for (int j = 0; j < num_channels; j++){
	TPEperPMT_beam.at(i).at(j) = flashes_beam.at(i)->get_PE(j);
      }
    }
    */
    // for (int i =0 ; i<flashes.size(); i++) {
    //     std::cout << flashes.at(i)->get_time() << "\n";
    // }
    OpReco->Fill();
    
  }
  output.Write();
  output.Close();
  
  return 0;
}
