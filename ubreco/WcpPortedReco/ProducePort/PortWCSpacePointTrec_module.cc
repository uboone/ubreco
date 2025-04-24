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

namespace WCTrecsp {
  class WCTrecPortedSpacePoints;
}

class WCTrecsp::WCTrecPortedSpacePoints : public art::EDProducer {
public:
  explicit WCTrecPortedSpacePoints(fhicl::ParameterSet const & p);
  WCTrecPortedSpacePoints(WCTrecPortedSpacePoints const &) = delete;
  WCTrecPortedSpacePoints(WCTrecPortedSpacePoints &&) = delete;
  WCTrecPortedSpacePoints & operator = (WCTrecPortedSpacePoints const &) = delete;
  WCTrecPortedSpacePoints & operator = (WCTrecPortedSpacePoints &&) = delete;

private:
  void produce(art::Event &e) override;
};

WCTrecsp::WCTrecPortedSpacePoints::WCTrecPortedSpacePoints(fhicl::ParameterSet const & p) : EDProducer{p}
{
  produces<std::vector<TrecSpacePoint>>();
}

void WCTrecsp::WCTrecPortedSpacePoints::produce(art::Event &e){

  auto outputSpacePointVec = std::make_unique<std::vector<TrecSpacePoint>>();

  std::cout << "Adding T_rec (WC no-trajectory-fitting neutrino cluster) spacepoints here:" << std::endl;

  std::string file = "./WCPwork/nue_" + std::to_string((int) e.run()) + "_" + std::to_string((int) e.subRun()) + "_" + std::to_string((int) e.id().event()) + ".root";
  std::cout << "loading file: " << file << std::endl;

  try {
    TFile *fin = new TFile(file.c_str());
    TTree *tin = (TTree*)fin->Get("T_rec");

    double x, y, z, q, nq;
    int cluster_id, real_cluster_id, sub_cluster_id;
    tin->SetBranchAddress("x", &x);
    tin->SetBranchAddress("y", &y);
    tin->SetBranchAddress("z", &z);
    tin->SetBranchAddress("q", &q);
    tin->SetBranchAddress("nq", &nq);
    tin->SetBranchAddress("cluster_id", &cluster_id);
    tin->SetBranchAddress("real_cluster_id", &real_cluster_id);
    tin->SetBranchAddress("sub_cluster_id", &sub_cluster_id);
    
    for(int i=0; i<tin->GetEntries(); i++){
      tin->GetEntry(i);
      TrecSpacePoint xyzqi = TrecSpacePoint{x, y, z, q, nq, cluster_id, real_cluster_id, sub_cluster_id};
      outputSpacePointVec->emplace_back(xyzqi);
    }

    fin->Close();
    std::cout << " space point vector size: "<<outputSpacePointVec->size()<<std::endl;
  } catch (std::exception &e) {
    std::cout << "Exception: " << e.what() << std::endl;
    std::cout << "Due to exception, adding 0 spacepoints..." << std::endl;
  }
  
  e.put(std::move(outputSpacePointVec));
}

DEFINE_ART_MODULE(WCTrecsp::WCTrecPortedSpacePoints)
