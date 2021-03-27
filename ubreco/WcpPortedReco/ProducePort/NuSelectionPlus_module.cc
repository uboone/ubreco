#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "ubobj/WcpPort/NuSelectionSTM.h"

#include <memory>
#include <string>
#include <dirent.h>
#include <iostream>
#include <fstream>

namespace nsm{
  class NuSelectionPlus;
}

class nsm::NuSelectionPlus : public art::EDProducer {
public:
  explicit NuSelectionPlus(fhicl::ParameterSet const & p);
  NuSelectionPlus(NuSelectionPlus const &) = delete;
  NuSelectionPlus(NuSelectionPlus &&) = delete;
  NuSelectionPlus & operator = (NuSelectionPlus const &) = delete;
  NuSelectionPlus & operator = (NuSelectionPlus &&) = delete;

private:
  void produce(art::Event &e) override;

  std::string fInput;
};

nsm::NuSelectionPlus::NuSelectionPlus(fhicl::ParameterSet const & p) : EDProducer{p}
{
  fInput = p.get<std::string>("PortInput");

  produces<std::vector<nsm::NuSelectionSTM> >();
}

void nsm::NuSelectionPlus::produce(art::Event &e){

  auto outputNsSTMVec = std::make_unique< std::vector<nsm::NuSelectionSTM> >();
  
  std::string path(fInput);
  std::cout<<"INPUT FILE NAME: "<<fInput<<std::endl;
  std::ifstream fin;
  fin.open(path.c_str(), std::ios::in);     
  if(!fin.good()) std::cout<<"INPUT FILE NOT FOUND!"<<std::endl;

  std::string runno="";
  int flash_id=-1;
  int cluster_id=-1;
  float flash_time=-1.;
  int event_type=-1;  
  int flag_low_energy=-1; 
  int flag_LM=-1; 
  int flag_TGM=-1;
  int flag_FC=-1; 
  int flag_STM=-1; 
  int flag_full_detector_dead=-1; 
  float cluster_length=-1.; 

  std::string event_runinfo = std::to_string((int)e.run())+"_"+std::to_string((int)e.subRun())+"_"+std::to_string((int)e.event());
 
  while(!fin.eof()){
	fin >> runno >> flash_id >> cluster_id >> flash_time >> event_type >> flag_low_energy >> flag_LM >> flag_TGM >> flag_FC
		>> flag_STM >> flag_full_detector_dead >> cluster_length;
	
	if( runno != event_runinfo ) continue;
	nsm::NuSelectionSTM nstm;
	nstm.SetEventType( event_type );
	nstm.SetLowEnergy( flag_low_energy );
	nstm.SetLM( flag_LM );
	nstm.SetTGM( flag_TGM );
	nstm.SetSTM( flag_STM );
	nstm.SetFullDead( flag_full_detector_dead );
	nstm.SetClusterLength( cluster_length );
        outputNsSTMVec->push_back( nstm );
	std::cout<<"EVENT FOUND: "<<runno<<std::endl;
	// For those events with multiple in-beam activity, this (the earliest one) will introduce an uncertainty! probably != what we ported in Reco1.5
  	break;
  } 

  e.put(std::move(outputNsSTMVec));
  fin.close();
}

DEFINE_ART_MODULE(nsm::NuSelectionPlus)
