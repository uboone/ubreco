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

#include "ubreco/WcpPortedReco/ProducePort/SpacePointStructs.h"

#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"

#include <memory>
#include <string>
#include <dirent.h>
#include <iostream>

namespace WCTrecchargesp {
  class WCTrecchargePortedSpacePoints;
}

class WCTrecchargesp::WCTrecchargePortedSpacePoints : public art::EDProducer {
public:
  explicit WCTrecchargePortedSpacePoints(fhicl::ParameterSet const & p);
  WCTrecchargePortedSpacePoints(WCTrecchargePortedSpacePoints const &) = delete;
  WCTrecchargePortedSpacePoints(WCTrecchargePortedSpacePoints &&) = delete;
  WCTrecchargePortedSpacePoints & operator = (WCTrecchargePortedSpacePoints const &) = delete;
  WCTrecchargePortedSpacePoints & operator = (WCTrecchargePortedSpacePoints &&) = delete;

private:
  void produce(art::Event &e) override;
};

WCTrecchargesp::WCTrecchargePortedSpacePoints::WCTrecchargePortedSpacePoints(fhicl::ParameterSet const & p) : EDProducer{p}
{
  produces<std::vector<TrecchargeSpacePoint>>();
}

void WCTrecchargesp::WCTrecchargePortedSpacePoints::produce(art::Event &e){

  auto outputSpacePointVec = std::make_unique<std::vector<TrecchargeSpacePoint>>();

  std::cout << "Adding T_rec_charge (WC trajectory-fitted neutrino-cluster) spacepoints here:" << std::endl;

  std::string file = "./WCPwork/nue_" + std::to_string((int) e.run()) + "_" + std::to_string((int) e.subRun()) + "_" + std::to_string((int) e.id().event()) + ".root";
  std::cout << "loading file: " << file << std::endl;

  try {
    TFile *fin = new TFile(file.c_str());
    TTree *tin = (TTree*)fin->Get("T_rec_charge");

    double x, y, z, q, nq, chi2, ndf, pu, pv, pw, pt, reduced_chi2, rr;
    int cluster_id, real_cluster_id, sub_cluster_id, flag_vertex, flag_shower;
    tin->SetBranchAddress("x", &x);
    tin->SetBranchAddress("y", &y);
    tin->SetBranchAddress("z", &z);
    tin->SetBranchAddress("q", &q);
    tin->SetBranchAddress("nq", &nq);
    tin->SetBranchAddress("cluster_id", &cluster_id);
    tin->SetBranchAddress("real_cluster_id", &real_cluster_id);
    tin->SetBranchAddress("sub_cluster_id", &sub_cluster_id);
    tin->SetBranchAddress("chi2", &chi2);
    tin->SetBranchAddress("ndf", &ndf);
    tin->SetBranchAddress("pu", &pu);
    tin->SetBranchAddress("pv", &pv);
    tin->SetBranchAddress("pw", &pw);
    tin->SetBranchAddress("pt", &pt);
    tin->SetBranchAddress("reduced_chi2", &reduced_chi2);
    tin->SetBranchAddress("flag_vertex", &flag_vertex);
    tin->SetBranchAddress("flag_shower", &flag_shower);
    tin->SetBranchAddress("rr", &rr);


    for(int i=0; i<tin->GetEntries(); i++){
      tin->GetEntry(i);
      TrecchargeSpacePoint xyzqi = TrecchargeSpacePoint{x, y, z, q, nq, cluster_id, real_cluster_id, sub_cluster_id, chi2, ndf, pu, pv, pw, pt, reduced_chi2, flag_vertex, flag_shower, rr};
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

DEFINE_ART_MODULE(WCTrecchargesp::WCTrecchargePortedSpacePoints)
