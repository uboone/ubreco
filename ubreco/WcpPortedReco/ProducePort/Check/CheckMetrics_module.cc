#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "ubreco/WcpPortedReco/ProducePort/WCP/NuSelectionContainment.h"
#include "ubreco/WcpPortedReco/ProducePort/WCP/NuSelectionMatch.h"
#include "ubreco/WcpPortedReco/ProducePort/WCP/NuSelectionTruth.h"
#include "ubreco/WcpPortedReco/ProducePort/WCP/NuSelectionCharge.h"

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

  art::Handle<std::vector<nsm::NuSelectionContainment> > containment_handle;
  e.getByLabel(fContainmentLabel,containment_handle);
  std::vector<art::Ptr<nsm::NuSelectionContainment> > containment_vec;
  art::fill_ptr_vector(containment_vec,containment_handle);
  std::cout<<"--- containment ---"<<std::endl;
  for(size_t i=0; i<containment_vec.size(); i++){
    art::Ptr<nsm::NuSelectionContainment> c = containment_vec.at(i);
    std::cout<<" flash found: "<<c->GetFlashFound()
	     <<"\n flash time: "<<c->GetFlashTime()
	     <<"\n measured PE: "<<c->GetFlashMeasPe()
	     <<"\n predicted PE: "<<c->GetFlashPredPe()
	     <<"\n match found? "<<c->GetMatchFound()
	     <<"\n match type: "<<c->GetMatchType()
	     <<"\n Is fully contained? "<<c->GetIsFC()
	     <<"\n Is through going muon? "<<c->GetIsTGM()
	     <<"\n Is NOT fully contained due to fiducial volume flag? "<<c->GetNotFCFV()
	     <<"\n Is NOT fully contained due to signal processing flag? "<<c->GetNotFCSP()
	     <<"\n Is NOT fully contained due to dead channel flag? "<<c->GetNotFCDC()
	     <<std::endl;
  }

  art::Handle<std::vector<nsm::NuSelectionCharge> > charge_handle;
  e.getByLabel(fChargeLabel,charge_handle);
  std::vector<art::Ptr<nsm::NuSelectionCharge> > charge_vec;
  art::fill_ptr_vector(charge_vec,charge_handle);
  for(size_t i=0; i<charge_vec.size(); i++){
    art::Ptr<nsm::NuSelectionCharge> c = charge_vec.at(i);
    std::cout<<" U plane integrated charge: "<<c->GetChargeU()
	     <<"\n V plane integrated charge: "<<c->GetChargeV()
	     <<"\n Y plane integrated charge: "<<c->GetChargeY()
	     <<std::endl;
  }

  if(fMC==true){

    art::Handle<std::vector<nsm::NuSelectionTruth> > truth_handle;
    e.getByLabel(fTruthLabel,truth_handle);
    std::vector<art::Ptr<nsm::NuSelectionTruth> > truth_vec;
    art::fill_ptr_vector(truth_vec,truth_handle);
    for(size_t i=0; i<truth_vec.size(); i++){
      art::Ptr<nsm::NuSelectionTruth> t = truth_vec.at(i);
      std::cout<<" Truth is CC? "<<t->GetIsCC()
	       <<"\n Truth is eligible? "<<t->GetIsEligible()
	       <<"\n Truth is fully conained? "<<t->GetIsFC()
	       <<"\n Truth vertex is in TPC active volume (5,5,10 bounds) ? "<<t->GetIsVtxInside()
	       <<"\n Truth neutrino PDG: "<<t->GetNuPdg()
	       <<"\n Truth vertex position: ( "<<t->GetVtxX()<<", "<<t->GetVtxY()<<","<<t->GetVtxZ()<<" )"
	       <<"\n Truth time: "<<t->GetTime()
	       <<"\n Truth neutrino energy: "<<t->GetNuEnergy()
	       <<"\n Truth neutrion energy (electrons) deposited inside TPC: "<<t->GetEnergyInside()<<" ("<<t->GetElectronInside()<<")"
	       <<std::endl;
    }

    art::Handle<std::vector<nsm::NuSelectionMatch> > match_handle;
    e.getByLabel(fMatchLabel,match_handle);
    std::vector<art::Ptr<nsm::NuSelectionMatch> > match_vec;
    art::fill_ptr_vector(match_vec,match_handle);
    for(size_t i=0; i<match_vec.size(); i++){
      art::Ptr<nsm::NuSelectionMatch> m = match_vec.at(i);
      std::cout<<" Match completeness: "<<m->GetCompleteness()
	       <<"\n Match completeness energy: "<<m->GetCompletenessEnergy()
	       <<"\n Match purity: "<<m->GetPurity()
	       <<"\n Match purity XY: "<<m->GetPurityXY()
	       <<"\n Match purity XZ: "<<m->GetPurityXZ()
	       <<"\n Match charge: "<<m->GetCharge()
	       <<"\n Match energy: "<<m->GetEnergy()
	       <<std::endl;
    }

  }
}

DEFINE_ART_MODULE(CheckMetrics)
