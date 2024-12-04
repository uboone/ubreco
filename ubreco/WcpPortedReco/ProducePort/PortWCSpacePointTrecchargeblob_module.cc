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

#include "ubreco/WcpPortedReco/ProducePort/SimpleSpacePointID.h"

#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"

#include <memory>
#include <string>
#include <dirent.h>
#include <iostream>

namespace WCTrecchargeblobsp {
  class WCTrecchargeblobPortedSpacePoints;
}

class WCTrecchargeblobsp::WCTrecchargeblobPortedSpacePoints : public art::EDProducer {
public:
  explicit WCTrecchargeblobPortedSpacePoints(fhicl::ParameterSet const & p);
  WCTrecchargeblobPortedSpacePoints(WCTrecchargeblobPortedSpacePoints const &) = delete;
  WCTrecchargeblobPortedSpacePoints(WCTrecchargeblobPortedSpacePoints &&) = delete;
  WCTrecchargeblobPortedSpacePoints & operator = (WCTrecchargeblobPortedSpacePoints const &) = delete;
  WCTrecchargeblobPortedSpacePoints & operator = (WCTrecchargeblobPortedSpacePoints &&) = delete;

private:
  void produce(art::Event &e) override;
};

WCTrecchargeblobsp::WCTrecchargeblobPortedSpacePoints::WCTrecchargeblobPortedSpacePoints(fhicl::ParameterSet const & p) : EDProducer{p}
{
  produces<std::vector<SimpleSpacePointID>>();
}

void WCTrecchargeblobsp::WCTrecchargeblobPortedSpacePoints::produce(art::Event &e){

  auto outputSpacePointVec = std::make_unique<std::vector<SimpleSpacePointID>>();

  std::cout << "Adding T_rec_charge_blob (WC trajectory-fitted neutrino-cluster) spacepoints here:" << std::endl;

  std::string file = "./WCPwork/nue_" + std::to_string((int) e.run()) + "_" + std::to_string((int) e.subRun()) + "_" + std::to_string((int) e.id().event()) + ".root";
  std::cout << "loading file: " << file << std::endl;

  try {
    TFile *fin = new TFile(file.c_str());
    TTree *tin = (TTree*)fin->Get("T_rec_charge_blob");

    double x, y, z, q;
    int real_cluster_id;
    tin->SetBranchAddress("x", &x);
    tin->SetBranchAddress("y", &y);
    tin->SetBranchAddress("z", &z);
    tin->SetBranchAddress("q", &q);
    tin->SetBranchAddress("real_cluster_id", &real_cluster_id);

    for(int i=0; i<tin->GetEntries(); i++){
      tin->GetEntry(i);
      SimpleSpacePointID xyzq = SimpleSpacePointID{x, y, z, q, real_cluster_id};
      outputSpacePointVec->emplace_back(xyzq);
    }

    fin->Close();
    std::cout << " space point vector size: "<<outputSpacePointVec->size()<<std::endl;
  } catch (std::exception &e) {
    // print exception
    std::cout << "Exception: " << e.what() << std::endl;
    std::cout << "Due to exception, adding 0 spacepoints..." << std::endl;
  }
  
  e.put(std::move(outputSpacePointVec));
}

DEFINE_ART_MODULE(WCTrecchargeblobsp::WCTrecchargeblobPortedSpacePoints)
