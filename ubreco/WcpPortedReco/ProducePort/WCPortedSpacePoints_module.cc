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

#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"

#include <memory>
#include <string>
#include <dirent.h>
#include <iostream>

namespace WCPsp {
  class WCPortedSpacePoints;
}

class WCPsp::WCPortedSpacePoints : public art::EDProducer {
public:
  explicit WCPortedSpacePoints(fhicl::ParameterSet const & p);
  WCPortedSpacePoints(WCPortedSpacePoints const &) = delete;
  WCPortedSpacePoints(WCPortedSpacePoints &&) = delete;
  WCPortedSpacePoints & operator = (WCPortedSpacePoints const &) = delete;
  WCPortedSpacePoints & operator = (WCPortedSpacePoints &&) = delete;

private:
  void produce(art::Event &e) override;

  std::string fInput;
  std::string fTreeName;
  bool fMainCluster;
  std::string fSpacePointLabel;
  short fTickOffset;
};

WCPsp::WCPortedSpacePoints::WCPortedSpacePoints(fhicl::ParameterSet const & p) : EDProducer{p}
{
  fInput = p.get<std::string>("PortInput");
  fTreeName = p.get<std::string>("TreeName");
  fMainCluster = p.get<bool>("MainCluster");
  fSpacePointLabel = p.get<std::string>("SpacePointLabel");
  fTickOffset = p.get<short>("TickOffset");

  produces<std::vector<std::array<float, 4>>>(); // x, y, z, q
}

void WCPsp::WCPortedSpacePoints::produce(art::Event &e){

  auto outputSpacePointVec = std::make_unique<std::vector<std::array<float, 4>>>();

  std::cout << "lhagaman modified, adding WC no-trajectory-fitting neutrino cluster spacepoints here:" << std::endl;

  std::string path(fInput);
  std::string file = (path);
  std::cout<<"INPUT FILE NAME: "<<file<<std::endl;
  TFile *fin = new TFile(file.c_str());
  TTree *tin = (TTree*)fin->Get(fTreeName.c_str());
 
  int run=-1, subrun=-1, event=-1, /*cluster_id=-1,*/ main_flag=-1, time_slice=-1, ch_u=-1, ch_v=-1, ch_w=-1;
  double x=-1., y=-1., z=-1., q=-1., nq=-1;
  tin->SetBranchAddress("run",&run);
  tin->SetBranchAddress("subrun",&subrun);
  tin->SetBranchAddress("event",&event);
  //tin->SetBranchAddress("cluster_id",&cluster_id);
  tin->SetBranchAddress("main_flag",&main_flag);
  tin->SetBranchAddress("time_slice",&time_slice);
  tin->SetBranchAddress("ch_u",&ch_u);
  tin->SetBranchAddress("ch_v",&ch_v);
  tin->SetBranchAddress("ch_w",&ch_w);
  tin->SetBranchAddress("x",&x);
  tin->SetBranchAddress("y",&y);
  tin->SetBranchAddress("z",&z);
  tin->SetBranchAddress("q",&q);
  tin->SetBranchAddress("nq",&nq);

  for(int i=0; i<tin->GetEntries(); i++){
    tin->GetEntry(i);

    if(fMainCluster==true && main_flag!=1) continue;

    if( run!=(int)e.run() || 
        subrun!=(int)e.subRun() || 
        event!=(int)e.id().event() ) continue;

    std::array<float, 4> xyzq = std::array<float, 4>{(float) x, (float) y, (float) z, (float) q};

    outputSpacePointVec->emplace_back(xyzq);
  }

  fin->Close();
  std::cout<<" space point vector size: "<<outputSpacePointVec->size()<<std::endl;
  e.put(std::move(outputSpacePointVec));
}

DEFINE_ART_MODULE(WCPsp::WCPortedSpacePoints)
