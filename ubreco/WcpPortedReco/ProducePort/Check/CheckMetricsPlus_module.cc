#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "ubobj/WcpPort/NuSelectionSTM.h"

#include "TTree.h"
#include "TFile.h"

class CheckMetricsPlus;

class CheckMetricsPlus : public art::EDAnalyzer {
public:
  explicit CheckMetricsPlus(fhicl::ParameterSet const& p);

  CheckMetricsPlus(CheckMetricsPlus const&) = delete;
  CheckMetricsPlus(CheckMetricsPlus&&) = delete;
  CheckMetricsPlus& operator=(CheckMetricsPlus const&) = delete;
  CheckMetricsPlus& operator=(CheckMetricsPlus&&) = delete;

  void reconfigure(fhicl::ParameterSet const& pset);
  void analyze(art::Event const& e) override;

private:

  std::string fSTMLabel;
};


CheckMetricsPlus::CheckMetricsPlus(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
  art::ServiceHandle<art::TFileService> tfs;
  this->reconfigure(p);
}

void CheckMetricsPlus::reconfigure(fhicl::ParameterSet const& pset)
{
  fSTMLabel = pset.get<std::string>("STMLabel");;
}

void CheckMetricsPlus::analyze(art::Event const& e)
{

  std::cout<<" RUN: "<<e.run()<<"\n SUBRUN: "<<e.subRun()<<"\n EVENT: "<<e.id().event()<<std::endl;

  art::Handle<std::vector<nsm::NuSelectionSTM> > match_handle;
  e.getByLabel(fSTMLabel,match_handle);
  std::vector<art::Ptr<nsm::NuSelectionSTM> > match_vec;
  art::fill_ptr_vector(match_vec,match_handle);
  for(size_t i=0; i<match_vec.size(); i++){
	  art::Ptr<nsm::NuSelectionSTM> m = match_vec.at(i);
	  std::cout<<"Event Type: "<<m->GetEventType()
		  <<"\n Flag Low Energy: "<<m->GetLowEnergy()
		  <<"\n Flag Light Mismatch: "<<m->GetLM()
		  <<"\n Flag TGM: "<<m->GetTGM()
		  <<"\n Flag STM: "<<m->GetSTM()
		  <<"\n Flag Full Detector Dead: "<<m->GetFullDead()
		  <<"\n Cluster Length: "<<m->GetClusterLength()
		  <<std::endl;
  }

}

DEFINE_ART_MODULE(CheckMetricsPlus)
