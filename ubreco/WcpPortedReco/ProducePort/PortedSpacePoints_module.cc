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

namespace psp {
  class PortedSpacePoints;
}

class psp::PortedSpacePoints : public art::EDProducer {
public:
  explicit PortedSpacePoints(fhicl::ParameterSet const & p);
  PortedSpacePoints(PortedSpacePoints const &) = delete;
  PortedSpacePoints(PortedSpacePoints &&) = delete;
  PortedSpacePoints & operator = (PortedSpacePoints const &) = delete;
  PortedSpacePoints & operator = (PortedSpacePoints &&) = delete;

private:
  void produce(art::Event &e) override;

  std::string fInput;
  std::string fTreeName;
  bool fMainCluster;
  std::string fSpacePointLabel;
  short fTickOffset;
};

psp::PortedSpacePoints::PortedSpacePoints(fhicl::ParameterSet const & p) : EDProducer{p}
{
  fInput = p.get<std::string>("PortInput");
  fTreeName = p.get<std::string>("TreeName");
  fMainCluster = p.get<bool>("MainCluster");
  fSpacePointLabel = p.get<std::string>("SpacePointLabel");
  fTickOffset = p.get<short>("TickOffset");

  produces<std::vector<recob::SpacePoint> >();
}

void psp::PortedSpacePoints::produce(art::Event &e){

  auto outputSpacePointVec = std::make_unique< std::vector<recob::SpacePoint> >();

  std::cout << "lhagaman modified, adding Pandora spacepoints here:" << std::endl;

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

    int id = -1;
    double xyz[3] = {0., 0., 0.};
    double xyz_err[6] = {0., 0., 0., 0., 0., 0.};
    double chisq = 0.;

    xyz[0]=x; // x alignment not confirmed
    xyz[1]=y;
    xyz[2]=z;
    recob::SpacePoint sp(xyz, xyz_err, chisq, id);
    outputSpacePointVec->emplace_back(sp);
  }

  fin->Close();
  std::cout<<" space point vector size: "<<outputSpacePointVec->size()<<std::endl;
  e.put(std::move(outputSpacePointVec));
}

DEFINE_ART_MODULE(psp::PortedSpacePoints)
