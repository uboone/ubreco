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

#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

class CheckHitsAgainst;


class CheckHitsAgainst : public art::EDAnalyzer {
public:
  explicit CheckHitsAgainst(fhicl::ParameterSet const& p);

  CheckHitsAgainst(CheckHitsAgainst const&) = delete;
  CheckHitsAgainst(CheckHitsAgainst&&) = delete;
  CheckHitsAgainst& operator=(CheckHitsAgainst const&) = delete;
  CheckHitsAgainst& operator=(CheckHitsAgainst&&) = delete;

  void reconfigure(fhicl::ParameterSet const& pset);
  void analyze(art::Event const& e) override;
  void setup_ttree(art::ServiceHandle<art::TFileService> tfs);

private:

  std::string fScLabel;
  std::string fLsHitLabel;
  std::string fHitLabel;
  std::string fWireLabel;
  bool fMC;
  int fTickOffset;

  TH2F *ls_hit_u, *ls_hit_v, *ls_hit_y;
  TH2F *hit_u, *hit_v, *hit_y;
  TH2F *wire_u, *wire_v, *wire_y;
  TH2F *sim_u, *sim_v, *sim_y;

  //double G4RefTime = -4050.0;
  double TriggerOffsetTPC = -1600.0;
  double DefaultTrigTime = 4050.0;
  double TickPeriod = 0.5;
};

CheckHitsAgainst::CheckHitsAgainst(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  ls_hit_u(nullptr),
  ls_hit_v(nullptr),
  ls_hit_y(nullptr),
  hit_u(nullptr),
  hit_v(nullptr),
  hit_y(nullptr),
  wire_u(nullptr),
  wire_v(nullptr),
  wire_y(nullptr),
  sim_u(nullptr),
  sim_v(nullptr),
  sim_y(nullptr)
{
  art::ServiceHandle<art::TFileService> tfs;
  setup_ttree(tfs);
  this->reconfigure(p);
}

void CheckHitsAgainst::reconfigure(fhicl::ParameterSet const& pset)
{
  fScLabel = pset.get<std::string>("ScLabel");;
  fHitLabel = pset.get<std::string>("HitLabel");;
  fLsHitLabel = pset.get<std::string>("LsHitLabel");;
  fWireLabel = pset.get<std::string>("WireLabel");;
  fMC = pset.get<bool>("MC");
  fTickOffset = pset.get<int>("TickOffset");
}

void CheckHitsAgainst::analyze(art::Event const& e)
{
  art::Handle<std::vector<recob::Hit> > hit_handle;
  e.getByLabel(fHitLabel,hit_handle);
  std::vector<art::Ptr<recob::Hit> > hit_vec;
  art::fill_ptr_vector(hit_vec,hit_handle);
  for(size_t i=0; i<hit_vec.size(); i++){
    art::Ptr<recob::Hit> hit = hit_vec.at(i);
    int refChannel = hit->Channel();
    float refTime = hit->PeakTime();
    float refCharge = hit->Integral();
    if(hit->View()==0){
      hit_u->SetBinContent(refChannel+1,(int)refTime,refCharge);
    }
    else if(hit->View()==1){
      hit_v->SetBinContent(refChannel+1-2400,(int)refTime,refCharge);
    }
    else{
      hit_y->SetBinContent(refChannel+1-4800,(int)refTime,refCharge);
    }
  }

  art::Handle<std::vector<recob::Hit> > ls_hit_handle;
  e.getByLabel(fLsHitLabel,ls_hit_handle);
  std::vector<art::Ptr<recob::Hit> > ls_hit_vec;
  art::fill_ptr_vector(ls_hit_vec,ls_hit_handle);
  for(size_t i=0; i<ls_hit_vec.size(); i++){
    art::Ptr<recob::Hit> hit = ls_hit_vec.at(i);
    int refChannel = hit->Channel();
    float refTime = hit->PeakTime();
    float refCharge = hit->Integral();
    if(hit->View()==0){
      ls_hit_u->SetBinContent(refChannel+1,(int)refTime,refCharge);
    }
    else if(hit->View()==1){
      ls_hit_v->SetBinContent(refChannel+1-2400,(int)refTime,refCharge);
    }
    else{
      ls_hit_y->SetBinContent(refChannel+1-4800,(int)refTime,refCharge);
    }
  }

  if(fMC==true){
    auto const& scHandle = e.getValidHandle<std::vector<sim::SimChannel>>(fScLabel);
    auto const& sc_vec(*scHandle);
    for(size_t i=0; i<sc_vec.size(); i++){
      
      sim::SimChannel sc = sc_vec[i];
      auto const timeSlices = sc.TDCIDEMap();
      for(auto timeSlice : timeSlices){
	
	auto tdc = timeSlice.first;
	
	auto const& energyDeposits = timeSlice.second;
	for(auto energyDeposit : energyDeposits){
	  
	  unsigned int tick = (unsigned int) tdc - ((DefaultTrigTime + TriggerOffsetTPC)/TickPeriod);
	  
	  if(sc.Channel()<2400){
	    sim_u->SetBinContent(sc.Channel()+1,tick-fTickOffset,energyDeposit.numElectrons);
	  }
	  if(sc.Channel()>=2400 && sc.Channel()<4800){
	    sim_v->SetBinContent(sc.Channel()-2400+1,tick-fTickOffset,energyDeposit.numElectrons);
	  }
	  if(sc.Channel()>4800){
	    sim_y->SetBinContent(sc.Channel()-4800+1,tick-fTickOffset,energyDeposit.numElectrons);
	  }
	}
      }
    }
  }
  if(fMC==false){
    art::Handle<std::vector<recob::Wire> > wire_handle;
    e.getByLabel(fWireLabel,wire_handle);
    std::vector<art::Ptr<recob::Wire> > wire_vec;
    art::fill_ptr_vector(wire_vec,wire_handle);
    for(size_t i=0; i<wire_vec.size(); i++){
      art::Ptr<recob::Wire> wire = wire_vec.at(i);
      int refChannel = wire->Channel();
      std::vector<float> wf = wire->Signal();
      int nbin = (int)wf.size();
      if(wire->View()==0){
	for(int j=0; j<nbin; j++){
	  wire_u->SetBinContent(refChannel+1,j+1,wf[j]);
	}
      }
      else if(wire->View()==1){
	for(int j=0; j<nbin; j++){
	  wire_v->SetBinContent(refChannel+1-2400,j+1,wf[j]);
	}
      }
      else{
	for(int j=0; j<nbin; j++){
	  wire_y->SetBinContent(refChannel+1-4800,j+1,wf[j]);
	}
      }
    }
  }
}

void CheckHitsAgainst::setup_ttree(art::ServiceHandle<art::TFileService> tfs)
{
  art::TFileDirectory histoDir = tfs->mkdir("Histo");
  ls_hit_u = histoDir.make<TH2F>("ls_hit_u","",2400,0,2400,6400,0,6400);
  ls_hit_v = histoDir.make<TH2F>("ls_hit_v","",2400,2400,4800,6400,0,6400);
  ls_hit_y = histoDir.make<TH2F>("ls_hit_y","",3456,4800,8256,6400,0,6400);
  hit_u = histoDir.make<TH2F>("hit_u","",2400,0,2400,6400,0,6400);
  hit_v = histoDir.make<TH2F>("hit_v","",2400,2400,4800,6400,0,6400);
  hit_y = histoDir.make<TH2F>("hit_y","",3456,4800,8256,6400,0,6400);
  wire_u = histoDir.make<TH2F>("wire_u","",2400,0,2400,6400,0,6400);
  wire_v = histoDir.make<TH2F>("wire_v","",2400,2400,4800,6400,0,6400);
  wire_y = histoDir.make<TH2F>("wire_y","",3456,4800,8256,6400,0,6400);
  sim_u = histoDir.make<TH2F>("sim_u","",2400,0,2400,6400,0,6400);
  sim_v = histoDir.make<TH2F>("sim_v","",2400,2400,4800,6400,0,6400);
  sim_y = histoDir.make<TH2F>("sim_y","",3456,4800,8256,6400,0,6400);
}


DEFINE_ART_MODULE(CheckHitsAgainst)
