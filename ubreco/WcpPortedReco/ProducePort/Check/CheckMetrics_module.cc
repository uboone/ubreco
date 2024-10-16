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

#include "ubobj/WcpPort/NuSelectionContainment.h"
#include "ubobj/WcpPort/NuSelectionMatch.h"
#include "ubobj/WcpPort/NuSelectionTruth.h"
#include "ubobj/WcpPort/NuSelectionCharge.h"

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

class CheckMetrics;

class CheckMetrics : public art::EDAnalyzer {
public:
  explicit CheckMetrics(fhicl::ParameterSet const& p);

  CheckMetrics(CheckMetrics const&) = delete;
  CheckMetrics(CheckMetrics&&) = delete;
  CheckMetrics& operator=(CheckMetrics const&) = delete;
  CheckMetrics& operator=(CheckMetrics&&) = delete;

  void reconfigure(fhicl::ParameterSet const& pset);
  void analyze(art::Event const& e) override;

private:

  std::string fContainmentLabel;
  std::string fChargeLabel;
  std::string fTruthLabel;
  std::string fMatchLabel;
  bool fMC;
};


CheckMetrics::CheckMetrics(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
  art::ServiceHandle<art::TFileService> tfs;
  this->reconfigure(p);
}

void CheckMetrics::reconfigure(fhicl::ParameterSet const& pset)
{
  fContainmentLabel = pset.get<std::string>("ContainmentLabel");;
  fChargeLabel = pset.get<std::string>("ChargeLabel");;
  fTruthLabel = pset.get<std::string>("TruthLabel");;
  fMatchLabel = pset.get<std::string>("MatchLabel");;
  fMC = pset.get<bool>("MC");
}

void CheckMetrics::analyze(art::Event const& e)
{

  std::cout<<" RUN: "<<e.run()<<"\n SUBRUN: "<<e.subRun()<<"\n EVENT: "<<e.id().event()<<std::endl;

  auto const& containment_vec = e.getProduct<std::vector<nsm::NuSelectionContainment>>(fContainmentLabel);
  std::cout<<"--- containment ---"<<std::endl;
  for(nsm::NuSelectionContainment const& c : containment_vec) {
    std::cout<<" flash found: "<<c.GetFlashFound()
             <<"\n flash time: "<<c.GetFlashTime()
             <<"\n measured PE: "<<c.GetFlashMeasPe()
             <<"\n predicted PE: "<<c.GetFlashPredPe()
             <<"\n match found? "<<c.GetMatchFound()
             <<"\n match type: "<<c.GetMatchType()
             <<"\n Is fully contained? "<<c.GetIsFC()
             <<"\n Is through going muon? "<<c.GetIsTGM()
             <<"\n Is NOT fully contained due to fiducial volume flag? "<<c.GetNotFCFV()
             <<"\n Is NOT fully contained due to signal processing flag? "<<c.GetNotFCSP()
             <<"\n Is NOT fully contained due to dead channel flag? "<<c.GetNotFCDC()
             <<"\n Match charge: "<<c.GetCharge()
             <<"\n Match energy: "<<c.GetEnergy()
	     <<std::endl;
  }

  auto const& charge_vec = e.getProduct<std::vector<nsm::NuSelectionCharge>>(fChargeLabel);
  for(nsm::NuSelectionCharge const& c : charge_vec) {
    std::cout<<" U plane integrated charge: "<<c.GetChargeU()
             <<"\n V plane integrated charge: "<<c.GetChargeV()
             <<"\n Y plane integrated charge: "<<c.GetChargeY()
	     <<std::endl;
  }

  if(fMC==true){

    auto const& truth_vec = e.getProduct<std::vector<nsm::NuSelectionTruth>>(fTruthLabel);
    for(nsm::NuSelectionTruth const& t : truth_vec) {
      std::cout<<" Truth is CC? "<<t.GetIsCC()
               <<"\n Truth is eligible? "<<t.GetIsEligible()
               <<"\n Truth is fully conained? "<<t.GetIsFC()
               <<"\n Truth vertex is in TPC active volume (5,5,10 bounds) ? "<<t.GetIsVtxInside()
               <<"\n Truth neutrino PDG: "<<t.GetNuPdg()
               <<"\n Truth vertex position: ( "<<t.GetVtxX()<<", "<<t.GetVtxY()<<","<<t.GetVtxZ()<<" )"
               <<"\n Truth time: "<<t.GetTime()
               <<"\n Truth neutrino energy: "<<t.GetNuEnergy()
               <<"\n Truth neutrion energy (electrons) deposited inside TPC: "<<t.GetEnergyInside()<<" ("<<t.GetElectronInside()<<")"
	       <<std::endl;
    }

    auto const& match_vec = e.getProduct<std::vector<nsm::NuSelectionMatch>>(fMatchLabel);
    for(nsm::NuSelectionMatch const& m : match_vec) {
      std::cout<<" Match completeness: "<<m.GetCompleteness()
               <<"\n Match completeness energy: "<<m.GetCompletenessEnergy()
               <<"\n Match purity: "<<m.GetPurity()
               <<"\n Match purity XY: "<<m.GetPurityXY()
               <<"\n Match purity XZ: "<<m.GetPurityXZ()
	       <<std::endl;
    }

  }
}

DEFINE_ART_MODULE(CheckMetrics)
