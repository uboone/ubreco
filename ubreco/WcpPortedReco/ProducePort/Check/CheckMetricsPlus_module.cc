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
  fSTMLabel = pset.get<std::string>("STMLabel");
}

void CheckMetricsPlus::analyze(art::Event const& e)
{

  std::cout<<" RUN: "<<e.run()<<"\n SUBRUN: "<<e.subRun()<<"\n EVENT: "<<e.id().event()<<std::endl;

  auto const& match_vec = e.getProduct<std::vector<nsm::NuSelectionSTM>>(fSTMLabel);
  for(nsm::NuSelectionSTM const& m : match_vec) {
          std::cout<<"Event Type: "<<m.GetEventType()
                  <<"\n Flag Low Energy: "<<m.GetLowEnergy()
                  <<"\n Flag Light Mismatch: "<<m.GetLM()
                  <<"\n Flag TGM: "<<m.GetTGM()
                  <<"\n Flag STM: "<<m.GetSTM()
                  <<"\n Flag Full Detector Dead: "<<m.GetFullDead()
                  <<"\n Cluster Length: "<<m.GetClusterLength()
		  <<std::endl;
  }

}

DEFINE_ART_MODULE(CheckMetricsPlus)
