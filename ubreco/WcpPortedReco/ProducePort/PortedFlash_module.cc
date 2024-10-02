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

#include "lardataobj/RecoBase/OpFlash.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "ubcore/Geometry/UBOpReadoutMap.h"

#include "TTree.h"
#include "TBranch.h"
#include "TH2F.h"
#include "TFile.h"

#include <memory>
#include <string>
#include <dirent.h>
#include <iostream>

namespace pf {
  class PortedFlash;
}

class pf::PortedFlash : public art::EDProducer {
public:
  explicit PortedFlash(fhicl::ParameterSet const & p);
  PortedFlash(PortedFlash const &) = delete;
  PortedFlash(PortedFlash &&) = delete;
  PortedFlash & operator = (PortedFlash const &) = delete;
  PortedFlash & operator = (PortedFlash &&) = delete;

private:
  void produce(art::Event &e) override;

  std::string fInput;
  std::string fTreeName;
};

pf::PortedFlash::PortedFlash(fhicl::ParameterSet const & p) : EDProducer{p}
{
  fInput = p.get<std::string>("PortInput");
  fTreeName = p.get<std::string>("TreeName");

  produces<std::vector<recob::OpFlash> >();
  produces<std::vector<int> >("flashType");
  produces<std::vector<double> >("flashLowTime");
  produces<std::vector<double> >("flashHighTime");
}

void pf::PortedFlash::produce(art::Event &e){

  auto const& channelMap = ::art::ServiceHandle<geo::WireReadout const>()->Get();

  auto outputFlashVec = std::make_unique< std::vector<recob::OpFlash> >();
  auto outputIntVec = std::make_unique< std::vector<int> >();
  auto outputDoubleVecA = std::make_unique< std::vector<double> >();
  auto outputDoubleVecB = std::make_unique< std::vector<double> >();

  std::string path(fInput);
  std::string file = (path);
  std::cout<<"INPUT FILE NAME: "<<file<<std::endl;
  TFile *fin = new TFile(file.c_str());
  TTree *tin = (TTree*)fin->Get(fTreeName.c_str());
  
  int run=-1, subrun=-1, event=-1, type=-1;
  double totalPE=-1., time=-1., low_time=-1., high_time=-1.;
  std::vector<double> *pe = new std::vector<double>;
  std::vector<double> *pe_err = new std::vector<double>;
  tin->SetBranchAddress("run",&run);
  tin->SetBranchAddress("subrun",&subrun);
  tin->SetBranchAddress("event",&event);
  tin->SetBranchAddress("type",&type);
  tin->SetBranchAddress("totalPE",&totalPE);
  tin->SetBranchAddress("time",&time);
  tin->SetBranchAddress("low_time",&low_time);
  tin->SetBranchAddress("high_time",&high_time);
  tin->SetBranchAddress("pe",&pe);
  tin->SetBranchAddress("pe_err",&pe_err);

  std::cout<<"==========================================="<<std::endl;
  std::cout<<"Run "<<e.run()<<"   Subrun "<< e.subRun()<<"    Event "<<e.id().event()<<std::endl;

  for(int i=0; i<tin->GetEntries(); i++){
    tin->GetEntry(i);
    if( run!=(int)e.run() || 
        subrun!=(int)e.subRun() || 
	event!=(int)e.id().event() ) continue;

    double Ycenter=0., Zcenter=0., Ywidth=0., Zwidth=0.;
    double sumy=0., sumz=0., sumy2=0., sumz2=0.;
    double totalPE=0.;

    std::cout<<pe->size()<<" pe vector size"<<std::endl;

    // may have multiple flashes in beam gate window (pe size > 32)
    // assuming incidence of 3 flashes in beam gate window is incredibly rare...
    size_t peSize = pe->size();

    std::vector<double> tempA(32,0.0);
    for(size_t j=0; j<32; j++){
      tempA.at(j)=pe->at(j);
      auto const PMTxyz = channelMap.OpDetGeoFromOpChannel(j).GetCenter();
      sumy += pe->at(j)*PMTxyz.Y();
      sumy2 += pe->at(j)*PMTxyz.Y()*PMTxyz.Y();
      sumz += pe->at(j)*PMTxyz.Z();
      sumz2 += pe->at(j)*PMTxyz.Z()*PMTxyz.Z();
      totalPE += pe->at(j);
    }
    Ycenter = sumy/totalPE;
    Zcenter = sumz/totalPE;
    if( (sumy2*totalPE - sumy*sumy) > 0. ){ Ywidth = std::sqrt(sumy2*totalPE - sumy*sumy) / totalPE; }
    if( (sumz2*totalPE - sumz*sumz) > 0. ){ Zwidth = std::sqrt(sumz2*totalPE - sumz*sumz) / totalPE; }

    recob::OpFlash flash(time,
			 high_time-low_time,
			 -1.,
			 2, 
			 tempA, //*pe,
			 1,
			 1,
			 1,
			 Ycenter,
			 Ywidth,
			 Zcenter,
			 Zwidth);

    outputFlashVec->emplace_back( flash );
    outputIntVec->push_back( type );
    outputDoubleVecA->push_back( low_time );
    outputDoubleVecB->push_back( high_time );

    if(peSize>32){
      Ycenter=0.; Zcenter=0.; Ywidth=0.; Zwidth=0.;
      sumy=0.; sumz=0.; sumy2=0.; sumz2=0.;
      totalPE=0.;
      std::vector<double> tempB(32,0.0);
      for(size_t j=32; j<64; j++){
	tempB.at(j-32)=pe->at(j);
        auto const PMTxyz = channelMap.OpDetGeoFromOpChannel(j-32).GetCenter();
        sumy += pe->at(j)*PMTxyz.Y();
        sumy2 += pe->at(j)*PMTxyz.Y()*PMTxyz.Y();
        sumz += pe->at(j)*PMTxyz.Z();
        sumz2 += pe->at(j)*PMTxyz.Z()*PMTxyz.Z();
	totalPE += pe->at(j);
      }
      Ycenter = sumy/totalPE;
      Zcenter = sumz/totalPE;
      if( (sumy2*totalPE - sumy*sumy) > 0. ){ Ywidth = std::sqrt(sumy2*totalPE - sumy*sumy) / totalPE; }
      if( (sumz2*totalPE - sumz*sumz) > 0. ){ Zwidth = std::sqrt(sumz2*totalPE - sumz*sumz) / totalPE; }
      recob::OpFlash flash(time,
			   high_time-low_time,
			   -1.,
			   2, 
			   tempB, //*pe,
			   1,
			   1,
			   1,
			   Ycenter,
			   Ywidth,
			   Zcenter,
			   Zwidth);
      outputFlashVec->emplace_back( flash );
      outputIntVec->push_back( type );
      outputDoubleVecA->push_back( low_time );
      outputDoubleVecB->push_back( high_time );
    }
  }

  fin->Close();
  std::cout<<" flash vector size: "<<outputFlashVec->size()<<std::endl;
  e.put(std::move(outputFlashVec));
  e.put(std::move(outputIntVec), "flashType");
  e.put(std::move(outputDoubleVecA), "flashLowTime");
  e.put(std::move(outputDoubleVecB), "flashHighTime");
}

DEFINE_ART_MODULE(pf::PortedFlash)
