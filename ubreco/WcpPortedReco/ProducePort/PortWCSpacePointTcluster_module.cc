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

#include "ubobj/WcpPort/SpacePointStructs.h"

#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"

#include <memory>
#include <string>
#include <dirent.h>
#include <iostream>

namespace WCTclustersp {
  class WCTclusterPortedSpacePoints;
}

class WCTclustersp::WCTclusterPortedSpacePoints : public art::EDProducer {
public:
  explicit WCTclusterPortedSpacePoints(fhicl::ParameterSet const & p);
  WCTclusterPortedSpacePoints(WCTclusterPortedSpacePoints const &) = delete;
  WCTclusterPortedSpacePoints(WCTclusterPortedSpacePoints &&) = delete;
  WCTclusterPortedSpacePoints & operator = (WCTclusterPortedSpacePoints const &) = delete;
  WCTclusterPortedSpacePoints & operator = (WCTclusterPortedSpacePoints &&) = delete;

private:
  void produce(art::Event &e) override;
};

WCTclustersp::WCTclusterPortedSpacePoints::WCTclusterPortedSpacePoints(fhicl::ParameterSet const & p) : EDProducer{p}
{
  produces<std::vector<TclusterSpacePoint>>();
}

void WCTclustersp::WCTclusterPortedSpacePoints::produce(art::Event &e){

  auto outputSpacePointVec = std::make_unique<std::vector<TclusterSpacePoint>>();

  std::cout << "Adding T_cluster (WC no-trajectory-fitting all-clusters) spacepoints here:" << std::endl;

  std::string file = "./WCPwork/nue_" + std::to_string((int) e.run()) + "_" + std::to_string((int) e.subRun()) + "_" + std::to_string((int) e.id().event()) + ".root";
  std::cout << "loading file: " << file << std::endl;

  try {
    TFile *fin = new TFile(file.c_str());
    TTree *tin = (TTree*)fin->Get("T_cluster");

    double x, y, z, q, nq;
    int cluster_id, real_cluster_id, sub_cluster_id, time_slice, ch_u, ch_v, ch_w;
    tin->SetBranchAddress("x", &x);
    tin->SetBranchAddress("y", &y);
    tin->SetBranchAddress("z", &z);
    tin->SetBranchAddress("q", &q);
    tin->SetBranchAddress("nq", &nq);
    tin->SetBranchAddress("cluster_id", &cluster_id);
    tin->SetBranchAddress("real_cluster_id", &real_cluster_id);
    tin->SetBranchAddress("sub_cluster_id", &sub_cluster_id);
    tin->SetBranchAddress("time_slice", &time_slice);
    tin->SetBranchAddress("ch_u", &ch_u);
    tin->SetBranchAddress("ch_v", &ch_v);
    tin->SetBranchAddress("ch_w", &ch_w);

    for(int i=0; i<tin->GetEntries(); i++){
      tin->GetEntry(i);
      TclusterSpacePoint xyzq = TclusterSpacePoint{x, y, z, q, nq, cluster_id, real_cluster_id, sub_cluster_id, time_slice, ch_u, ch_v, ch_w};
      outputSpacePointVec->emplace_back(xyzq);
    }

    fin->Close();
    std::cout << " space point vector size: "<<outputSpacePointVec->size()<<std::endl;
  } catch (std::exception &e) {
    std::cout << "Exception: " << e.what() << std::endl;
    std::cout << "Due to exception, adding 0 spacepoints..." << std::endl;
  }
  
  e.put(std::move(outputSpacePointVec));
}

DEFINE_ART_MODULE(WCTclustersp::WCTclusterPortedSpacePoints)
