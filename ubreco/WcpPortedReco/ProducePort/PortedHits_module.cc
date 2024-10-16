#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/Hit.h"
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "TTree.h"
#include "TBranch.h"
#include "TH2F.h"
#include "TFile.h"

#include <memory>
#include <string>
#include <dirent.h>
#include <iostream>

namespace ph {
  class PortedHits;
}

class ph::PortedHits : public art::EDProducer {
public:
  explicit PortedHits(fhicl::ParameterSet const & p);
  PortedHits(PortedHits const &) = delete;
  PortedHits(PortedHits &&) = delete;
  PortedHits & operator = (PortedHits const &) = delete;
  PortedHits & operator = (PortedHits &&) = delete;

private:
  void produce(art::Event &e) override;

  std::string fInput;
  std::string fTreeName;
  bool fMainCluster;
  short fTickOffset;
  float fChargeScaling;
  float fRebin;
};

ph::PortedHits::PortedHits(fhicl::ParameterSet const & p) : EDProducer{p}
{
  fInput = p.get<std::string>("PortInput");
  fTreeName = p.get<std::string>("TreeName");
  fMainCluster = p.get<bool>("MainCluster");
  fTickOffset = p.get<short>("TickOffset");
  fChargeScaling = p.get<float>("ChargeScaling");
  fRebin = p.get<float>("Rebin");
  produces<std::vector<recob::Hit> >();
}

void ph::PortedHits::produce(art::Event &e){
  std::unique_ptr<std::vector<recob::Hit> >hit_collection(new std::vector<recob::Hit>);

  std::string path(fInput);

  std::string file = (path);
  std::cout<<"INPUT FILE NAME: "<<file<<std::endl;
  TFile *fin = new TFile(file.c_str());
  TTree *tin = (TTree*)fin->Get(fTreeName.c_str());
  
  int run, subrun, event;
  int channel, start_tick, main_flag;
  float charge, charge_error;
  tin->SetBranchAddress("run",&run);
  tin->SetBranchAddress("subrun",&subrun);
  tin->SetBranchAddress("event",&event);
  tin->SetBranchAddress("channel",&channel);
  tin->SetBranchAddress("start_tick",&start_tick);
  tin->SetBranchAddress("main_flag",&main_flag);
  tin->SetBranchAddress("charge",&charge);
  tin->SetBranchAddress("charge_error",&charge_error);

  //auto const &detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
  art::ServiceHandle<geo::Geometry> geom;
  int cryostat_no=0, tpc_no=0, plane_no=0;
  for(int i=0; i<tin->GetEntries(); i++){
    tin->GetEntry(i);

    if(fMainCluster==true && main_flag!=1) continue;

    if(run!=(int)e.run() || subrun!=(int)e.subRun() || event!=(int)e.id().event()) continue;

    if(channel<2400){ 
      plane_no=0; 
    }
    else if(channel>=2400 && channel<4800){ 
      plane_no=1; 
      channel=channel-2400; 
    }
    else{ 
      plane_no=2; 
      channel=channel-4800; 
    }
    geo::WireID wire(cryostat_no,tpc_no,plane_no,channel);
    raw::ChannelID_t chan = geom->PlaneWireToChannel(wire); // removed -1 

    //raw::TDCtick_t start_tdc = detClocks->TPCTick2TDC(start_tick-fTickOffset);
    //raw::TDCtick_t end_tdc = detClocks->TPCTick2TDC(start_tick+(fRebin-1)-fTickOffset);
	
    hit_collection->push_back(recob::Hit(chan,
					 start_tick-fTickOffset,
					 start_tick+3-fTickOffset,
					 start_tick+(fRebin/2.)-fTickOffset,
					 fRebin/2.,
					 fRebin/2.,
					 charge*fChargeScaling/fRebin, 
					 charge_error*fChargeScaling/fRebin,
					 charge*fChargeScaling,
					 charge*fChargeScaling,
					 charge*fChargeScaling,
					 charge_error*fChargeScaling,
					 -1, // 1
					 0,
					 1,
					 1,
					 geom->View(chan),
					 geom->SignalType(chan),
					 wire));
  }
  
  fin->Close();

  std::cout<<" [threshold] hit collection size: "<<hit_collection->size()<<std::endl;
  e.put(std::move(hit_collection));
}

DEFINE_ART_MODULE(ph::PortedHits)
