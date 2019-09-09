////////////////////////////////////////////////////////////////////////
// Class:       CheckFlash
// Plugin Type: analyzer (art v3_01_02)
// File:        CheckFlash_module.cc
//
// Generated at Thu Jul  4 11:47:48 2019 by Brooke Russell using cetskelgen
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

#include "lardataobj/RecoBase/OpFlash.h"

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

class CheckFlash;

class CheckFlash : public art::EDAnalyzer {
public:
  explicit CheckFlash(fhicl::ParameterSet const& p);

  CheckFlash(CheckFlash const&) = delete;
  CheckFlash(CheckFlash&&) = delete;
  CheckFlash& operator=(CheckFlash const&) = delete;
  CheckFlash& operator=(CheckFlash&&) = delete;

  void reconfigure(fhicl::ParameterSet const& pset);
  void analyze(art::Event const& e) override;
  void setup_ttree(art::ServiceHandle<art::TFileService> tfs);

private:

  TH1F *ls, *wcp_1, *wcp_2, *wcp_3;

  float fLsLow, fLsHigh;
  std::string fLsFlashLabel;
  std::string fWcpFlashLabel;
};


CheckFlash::CheckFlash(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  ls(nullptr),
  wcp_1(nullptr),
  wcp_2(nullptr),
  wcp_3(nullptr)
{
  art::ServiceHandle<art::TFileService> tfs;
  setup_ttree(tfs);
  this->reconfigure(p);
}

void CheckFlash::reconfigure(fhicl::ParameterSet const& pset)
{
  fLsLow = pset.get<float>("LsLow");
  fLsHigh = pset.get<float>("LsHigh");
  fLsFlashLabel = pset.get<std::string>("LsFlashLabel");
  fWcpFlashLabel = pset.get<std::string>("WcpFlashLabel");
}

void CheckFlash::analyze(art::Event const& e)
{
  art::Handle<std::vector<recob::OpFlash> > ls_handle;
  e.getByLabel(fLsFlashLabel,ls_handle);
  std::vector<art::Ptr<recob::OpFlash> > ls_vec;
  art::fill_ptr_vector(ls_vec,ls_handle);
  for(size_t i=0; i<ls_vec.size(); i++){
    art::Ptr<recob::OpFlash> flash = ls_vec.at(i);
    double time = flash->Time();
    if(time < fLsLow || time > fLsHigh ) continue;
    std::cout<<"LArSoft flash time: "<< time <<std::endl;
    std::vector<double> pe = flash->PEs();
    for(size_t j=0; j<pe.size(); j++){
      ls->SetBinContent(j+1, pe.at(j));
    }
  }

  art::Handle<std::vector<recob::OpFlash> > wcp_handle;
  e.getByLabel(fWcpFlashLabel,wcp_handle);
  std::vector<art::Ptr<recob::OpFlash> > wcp_vec;
  art::fill_ptr_vector(wcp_vec,wcp_handle);
  for(size_t i=0; i<wcp_vec.size(); i++){
    art::Ptr<recob::OpFlash> flash = wcp_vec.at(i);
    double time = flash->Time();
    std::cout<<"WCP flash time: "<< time <<std::endl;
    std::vector<double> pe = flash->PEs();
    for(size_t j=0; j<pe.size(); j++){
      if(i==0){ wcp_1->SetBinContent(j+1, pe.at(j)); }
      if(i==1){ wcp_2->SetBinContent(j+1, pe.at(j)); }
      if(i==2){ wcp_3->SetBinContent(j+1, pe.at(j)); }
    }
  }
}

void CheckFlash::setup_ttree(art::ServiceHandle<art::TFileService> tfs)
{
  art::TFileDirectory histoDir = tfs->mkdir("Histo");
  ls = histoDir.make<TH1F>("ls","",32,0,32);
  wcp_1 = histoDir.make<TH1F>("wcp_1","",32,0,32);
  wcp_2 = histoDir.make<TH1F>("wcp_2","",32,0,32);
  wcp_3 = histoDir.make<TH1F>("wcp_3","",32,0,32);
}

DEFINE_ART_MODULE(CheckFlash)
