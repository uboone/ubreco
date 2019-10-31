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
#include "ubobj/WcpPort/NuSelectionContainment.h"
#include "ubobj/WcpPort/NuSelectionMatch.h"
#include "ubobj/WcpPort/NuSelectionTruth.h"
#include "ubobj/WcpPort/NuSelectionCharge.h"

#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"

#include <memory>
#include <string>
#include <dirent.h>
#include <iostream>

namespace nsm{
  class NuSelectionMetrics;
}

class nsm::NuSelectionMetrics : public art::EDProducer {
public:
  explicit NuSelectionMetrics(fhicl::ParameterSet const & p);
  NuSelectionMetrics(NuSelectionMetrics const &) = delete;
  NuSelectionMetrics(NuSelectionMetrics &&) = delete;
  NuSelectionMetrics & operator = (NuSelectionMetrics const &) = delete;
  NuSelectionMetrics & operator = (NuSelectionMetrics &&) = delete;

private:
  void produce(art::Event &e) override;

  std::string fInput;
  std::string fTreeEvalName;
  std::string fTreeChargeName;
  bool fMainCluster;
  bool fMC;
};

nsm::NuSelectionMetrics::NuSelectionMetrics(fhicl::ParameterSet const & p) : EDProducer{p}
{
  fInput = p.get<std::string>("PortInput");
  fTreeEvalName = p.get<std::string>("TreeEvalName");
  fTreeChargeName = p.get<std::string>("TreeChargeName");
  fMainCluster = p.get<bool>("MainCluster");
  fMC = p.get<bool>("MC");

  produces<std::vector<nsm::NuSelectionContainment> >();
  produces<std::vector<nsm::NuSelectionCharge> >();
  if(fMC==true){
    produces<std::vector<nsm::NuSelectionMatch> >();
    produces<std::vector<nsm::NuSelectionTruth> >();
  }
}

void nsm::NuSelectionMetrics::produce(art::Event &e){

  auto outputNscontainmentVec = std::make_unique< std::vector<nsm::NuSelectionContainment> >();
  auto outputNschargeVec = std::make_unique< std::vector<nsm::NuSelectionCharge> >();
  auto outputNsmatchVec = std::make_unique< std::vector<nsm::NuSelectionMatch> >();
  auto outputNstruthVec = std::make_unique< std::vector<nsm::NuSelectionTruth> >();
  
  std::string path(fInput);
  std::string file = (path);
  std::cout<<"INPUT FILE NAME: "<<file<<std::endl;
  TFile *fin = new TFile(file.c_str());
  TTree *tin = (TTree*)fin->Get(fTreeEvalName.c_str());
  
  int run=-1, subrun=-1, event=-1;

  float truth_nuEnergy=-1., truth_energyInside=-1., truth_electronInside=-1.;
  int truth_nuPdg=-1;
  bool truth_isCC=false, truth_isEligible=false, truth_isFC=false, truth_vtxInside=false;
  float truthVtxX=-1., truthVtxY=-1., truthVtxZ=-1., truth_nuTime=-1.;

  bool flash_found=false;
  float flash_time=-1., flash_measPe=-1., flash_predPe=-1.;

  bool match_found=false;
  float match_completeness=-1., match_completeness_energy=-1., match_purity=-1., match_purity_xy=-1., match_purity_xz=-1., match_charge=-1., match_energy=-1.;
  unsigned int match_type=0;
  bool match_isFC=false, match_isTgm=false, match_notFC_FV=false, match_notFC_SP=false, match_notFC_DC=false;

  tin->SetBranchAddress("run",&run);
  tin->SetBranchAddress("subrun",&subrun);
  tin->SetBranchAddress("event",&event);
  tin->SetBranchAddress("flash_found",&flash_found);
  tin->SetBranchAddress("flash_time",&flash_time);
  tin->SetBranchAddress("flash_measPe",&flash_measPe);
  tin->SetBranchAddress("flash_predPe",&flash_predPe);
  tin->SetBranchAddress("match_found",&match_found);
  tin->SetBranchAddress("match_type",&match_type);
  tin->SetBranchAddress("match_isFC",&match_isFC);
  tin->SetBranchAddress("match_isTgm",&match_isTgm);
  tin->SetBranchAddress("match_notFC_FV",&match_notFC_FV);
  tin->SetBranchAddress("match_notFC_SP",&match_notFC_SP);
  tin->SetBranchAddress("match_notFC_DC",&match_notFC_DC);
  tin->SetBranchAddress("match_charge",&match_charge);
  tin->SetBranchAddress("match_energy",&match_energy);

  if(fMC==true){
    tin->SetBranchAddress("truth_nuEnergy",&truth_nuEnergy);
    tin->SetBranchAddress("truth_energyInside",&truth_energyInside);
    tin->SetBranchAddress("truth_electronInside",&truth_electronInside);
    tin->SetBranchAddress("truth_nuPdg",&truth_nuPdg);
    tin->SetBranchAddress("truth_isCC",&truth_isCC);
    tin->SetBranchAddress("truth_isEligible",&truth_isEligible);
    tin->SetBranchAddress("truth_isFC",&truth_isFC);
    tin->SetBranchAddress("truth_vtxInside",&truth_vtxInside);
    tin->SetBranchAddress("truth_vtxX",&truthVtxX);
    tin->SetBranchAddress("truth_vtxY",&truthVtxY);
    tin->SetBranchAddress("truth_vtxZ",&truthVtxZ);
    tin->SetBranchAddress("truth_nuTime",&truth_nuTime);

    tin->SetBranchAddress("match_completeness",&match_completeness);
    tin->SetBranchAddress("match_completeness_energy",&match_completeness_energy);
    tin->SetBranchAddress("match_purity",&match_purity);
    tin->SetBranchAddress("match_purity_xy",&match_purity_xy);
    tin->SetBranchAddress("match_purity_xz",&match_purity_xz);
  }

  
  for(int i=0; i<tin->GetEntries(); i++){
    tin->GetEntry(i);
    if( run!=(int)e.run() || 
        subrun!=(int)e.subRun() || 
	event!=(int)e.id().event() ) continue;

    nsm::NuSelectionContainment nsc;
    nsc.SetFlashFound( flash_found );
    nsc.SetFlashTime( flash_time );
    nsc.SetFlashMeasPe( flash_measPe );
    nsc.SetFlashPredPe( flash_predPe );
    nsc.SetMatchFound( match_found );
    nsc.SetMatchType( match_type );
    nsc.SetIsFC( match_isFC );
    nsc.SetIsTGM( match_isTgm );
    nsc.SetNotFCFV( match_notFC_FV );
    nsc.SetNotFCSP( match_notFC_SP );
    nsc.SetNotFCDC( match_notFC_DC );
    nsc.SetCharge( match_charge );
    nsc.SetEnergy( match_energy );
    outputNscontainmentVec->push_back( nsc );

    if(fMC==true){
      nsm::NuSelectionMatch nsm;
      nsm.SetCompleteness( match_completeness );
      nsm.SetCompletenessEnergy( match_completeness_energy );
      nsm.SetPurity( match_purity );
      nsm.SetPurityXY( match_purity_xy );
      nsm.SetPurityXZ( match_purity_xz );
      outputNsmatchVec->push_back( nsm );

      nsm::NuSelectionTruth nst;
      nst.SetIsCC( truth_isCC );
      nst.SetIsEligible( truth_isEligible );
      nst.SetIsFC( truth_isFC );
      nst.SetVtxInside( truth_vtxInside );
      nst.SetNuPdg( truth_nuPdg );
      nst.SetVtxX( truthVtxX );
      nst.SetVtxY( truthVtxY );
      nst.SetVtxZ( truthVtxZ );
      nst.SetTime( truth_nuTime );
      nst.SetNuEnergy( truth_nuEnergy );
      nst.SetEnergyInside( truth_energyInside );
      nst.SetElectronInside( truth_electronInside );
      outputNstruthVec->push_back( nst );
    }
  }

  TTree *t = (TTree*)fin->Get(fTreeChargeName.c_str());
  int channel, start_tick, main_flag;
  float charge, charge_error;
  t->SetBranchAddress("run",&run);
  t->SetBranchAddress("subrun",&subrun);
  t->SetBranchAddress("event",&event);
  t->SetBranchAddress("channel",&channel);
  t->SetBranchAddress("start_tick",&start_tick);
  t->SetBranchAddress("main_flag",&main_flag);
  t->SetBranchAddress("charge",&charge);
  t->SetBranchAddress("charge_error",&charge_error);

  float u=0., v=0., y=0.;
  for(int i=0; i<t->GetEntries(); i++){
    t->GetEntry(i);

    if(fMainCluster==true && main_flag!=1) continue;

    if( run!=(int)e.run() || 
        subrun!=(int)e.subRun() || 
	event!=(int)e.id().event() ) continue;

    if(channel<2400){ u += charge; }
    else if(channel>=2400 && channel<4800){ v += charge; }
    else{ y += charge; }

  }

  NuSelectionCharge nscharge;
  nscharge.SetChargeU( u );
  nscharge.SetChargeV( v );
  nscharge.SetChargeY( y );
  outputNschargeVec->push_back( nscharge );

  e.put(std::move(outputNscontainmentVec));
  e.put(std::move(outputNschargeVec));
  if(fMC==true){
    e.put(std::move(outputNsmatchVec));
    e.put(std::move(outputNstruthVec));
  }
  fin->Close();
  delete fin;
}

DEFINE_ART_MODULE(nsm::NuSelectionMetrics)
