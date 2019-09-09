////////////////////////////////////////////////////////////////////////
// Class:       CheckWires
// Plugin Type: analyzer (art v3_01_02)
// File:        CheckWires_module.cc
//
// Generated at Sun Jun 30 14:01:53 2019 by Brooke Russell using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

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
#include "lardataobj/RecoBase/Wire.h"

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

class CheckWires;

class CheckWires : public art::EDAnalyzer {
public:
  explicit CheckWires(fhicl::ParameterSet const& p);

  CheckWires(CheckWires const&) = delete;
  CheckWires(CheckWires&&) = delete;
  CheckWires& operator=(CheckWires const&) = delete;
  CheckWires& operator=(CheckWires&&) = delete;

  void reconfigure(fhicl::ParameterSet const& pset);
  void analyze(art::Event const& e) override;
  void setup_ttree(art::ServiceHandle<art::TFileService> tfs);
private:

  std::string fScLabel;
  std::string fSigprocLabel;
  std::string fWireLabel;
  int fTickOffset;

  TH2F *wire_u, *wire_v, *wire_y;
  TH2F *sp_u, *sp_v, *sp_y;
  TH2F *sim_u, *sim_v, *sim_y;

  double G4RefTime = -4050.0;
  double TriggerOffsetTPC = -1600.0;
  double DefaultTrigTime = 4050.0;
  double TickPeriod = 0.5;
};

CheckWires::CheckWires(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  wire_u(nullptr),
  wire_v(nullptr),
  wire_y(nullptr),
  sp_u(nullptr),
  sp_v(nullptr),
  sp_y(nullptr),
  sim_u(nullptr),
  sim_v(nullptr),
  sim_y(nullptr)
{
  art::ServiceHandle<art::TFileService> tfs;
  setup_ttree(tfs);
  this->reconfigure(p);
}

void CheckWires::reconfigure(fhicl::ParameterSet const& pset)
{
  fScLabel = pset.get<std::string>("ScLabel");;
  fSigprocLabel = pset.get<std::string>("SigprocLabel");
  fWireLabel = pset.get<std::string>("WireLabel");;
  fTickOffset = pset.get<int>("TickOffset");
}

void CheckWires::analyze(art::Event const& e)
{
  art::Handle<std::vector<recob::Wire> > wire_handle;
  e.getByLabel(fWireLabel,wire_handle);
  std::vector<art::Ptr<recob::Wire> > wire_vec;
  art::fill_ptr_vector(wire_vec,wire_handle);
  for(size_t i=0; i<wire_vec.size(); i++){
    art::Ptr<recob::Wire> wire = wire_vec.at(i);
    int refChannel = wire->Channel();
    std::vector<float> wf = wire->Signal();
    int nbin = (int)wf.size();
    if(refChannel<2400){
      for(int j=0; j<nbin; j++){
	wire_u->SetBinContent(refChannel+1,j+1,wf[j]);
      }
    }
    else if(refChannel>=2400 && refChannel<4800){
      for(int j=0; j<nbin; j++){
	wire_v->SetBinContent(refChannel-2400+1,j+1,wf[j]);
      }
    }
    else{
      for(int j=0; j<nbin; j++){
	wire_y->SetBinContent(refChannel-4800+1,j+1,wf[j]);
      }
    }
  }

  art::Handle<std::vector<recob::Wire> > sp_wire_handle;
  e.getByLabel(fSigprocLabel,sp_wire_handle);
  std::vector<art::Ptr<recob::Wire> > sp_wire_vec;
  art::fill_ptr_vector(sp_wire_vec,sp_wire_handle);
  for(size_t i=0; i<sp_wire_vec.size(); i++){
    art::Ptr<recob::Wire> wire = sp_wire_vec.at(i);
    int refChannel = wire->Channel();
    std::vector<float> wf = wire->Signal();
    int nbin = (int)wf.size();
    if(refChannel<2400){
      for(int j=0; j<nbin; j++){
	sp_u->SetBinContent(refChannel+1,j+1,wf[j]);
      }
    }
    else if(refChannel>=2400 && refChannel<4800){
      for(int j=0; j<nbin; j++){
	sp_v->SetBinContent(refChannel-2400+1,j+1,wf[j]);
      }
    }
    else{
      for(int j=0; j<nbin; j++){
	sp_y->SetBinContent(refChannel-4800+1,j+1,wf[j]);
      }
    }
  }

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
	  sim_u->SetBinContent(sc.Channel()+1,tick+1-fTickOffset,energyDeposit.numElectrons);
	}
	if(sc.Channel()>=2400 && sc.Channel()<4800){
	  sim_v->SetBinContent(sc.Channel()-2400+1,tick+1-fTickOffset,energyDeposit.numElectrons);
	}
	if(sc.Channel()>4800){
	  sim_y->SetBinContent(sc.Channel()-4800+1,tick+1-fTickOffset,energyDeposit.numElectrons);
	}
      }
    }
  }
}

void CheckWires::setup_ttree(art::ServiceHandle<art::TFileService> tfs)
{
  art::TFileDirectory histoDir = tfs->mkdir("Histo");
  wire_u = histoDir.make<TH2F>("wire_u","",2400,0,2400,6400,0,6400);
  wire_v = histoDir.make<TH2F>("wire_v","",2400,2400,4800,6400,0,6400);
  wire_y = histoDir.make<TH2F>("wire_y","",3456,4800,8256,6400,0,6400);
  sp_u = histoDir.make<TH2F>("sp_u","",2400,0,2400,6400,0,6400);
  sp_v = histoDir.make<TH2F>("sp_v","",2400,2400,4800,6400,0,6400);
  sp_y = histoDir.make<TH2F>("sp_y","",3456,4800,8256,6400,0,6400);
  sim_u = histoDir.make<TH2F>("sim_u","",2400,0,2400,6400,0,6400);
  sim_v = histoDir.make<TH2F>("sim_v","",2400,2400,4800,6400,0,6400);
  sim_y = histoDir.make<TH2F>("sim_y","",3456,4800,8256,6400,0,6400);
}

DEFINE_ART_MODULE(CheckWires)
