////////////////////////////////////////////////////////////////////////
// Class:       WCPselection
// Plugin Type: analyzer (art v3_01_02)
// File:        WCPselection_module.cc
//
// Generated at Fri Oct  4 15:24:42 2019 by Hanyu Wei using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "ubobj/WcpPort/NuSelectionContainment.h"
#include "ubobj/WcpPort/NuSelectionMatch.h"
#include "ubobj/WcpPort/NuSelectionTruth.h"
#include "ubobj/WcpPort/NuSelectionCharge.h"


#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TString.h"

#include <iostream>
#include <fstream>


class WCPselection;


class WCPselection : public art::EDAnalyzer {
public:
  explicit WCPselection(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WCPselection(WCPselection const&) = delete;
  WCPselection(WCPselection&&) = delete;
  WCPselection& operator=(WCPselection const&) = delete;
  WCPselection& operator=(WCPselection&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void endSubRun(art::SubRun const& sr) override;

  // user defined
  void reconfigure(fhicl::ParameterSet const& pset);
  void initOutput();

private:

  // Declare member data here.

  // fcl config
  std::string fContainmentLabel;
  std::string fChargeLabel;
  std::string fTruthLabel;
  std::string fMatchLabel;
  bool fMC;
  
  // output 
  TTree* fTreeEval; 
  Int_t           f_run;
  Int_t           f_subRun;
  Int_t           f_event;
  Bool_t          f_flash_found;
  Float_t         f_flash_time;
  Float_t         f_flash_measPe;
  Float_t         f_flash_predPe;
  Bool_t          f_match_found;
  UInt_t          f_match_type;
  Bool_t          f_match_isFC;
  Bool_t          f_match_isTgm;
  Bool_t          f_match_notFC_FV;
  Bool_t          f_match_notFC_SP;
  Bool_t          f_match_notFC_DC;
  Float_t         f_match_charge; // main flag collection plane charge
  Float_t         f_match_energy; 
  Float_t	  f_match_chargeU; 
  Float_t	  f_match_chargeV; 
  Float_t	  f_match_chargeY; 
  Float_t	  f_match_energyY;
  Bool_t	  f_lightmismatch;
 
  Float_t         f_truth_nuEnergy;
  Float_t         f_truth_energyInside;
  Float_t         f_truth_electronInside;
  Int_t           f_truth_nuPdg;
  Bool_t          f_truth_isCC;
  Bool_t          f_truth_isEligible;
  Bool_t          f_NC_truth_isEligible;
  Bool_t          f_truth_isFC;
  Bool_t          f_truth_vtxInside;
  Float_t         f_truth_vtxX;
  Float_t         f_truth_vtxY;
  Float_t         f_truth_vtxZ;
  Float_t         f_truth_nuTime;
  Float_t         f_match_completeness;
  Float_t         f_match_completeness_energy;
  Float_t         f_match_purity;
  Float_t         f_match_purity_xz;
  Float_t         f_match_purity_xy;
 

  TTree* fTreePot;
  std::string fPOT_inputTag; 
  bool fPOT_counting=false; 
 // beamdata:bnbETOR875
 // EXT?
 // MC [label: geneartor] : [instance: ]

  int frun;
  int fsubRun;  
  double fpot_tor875;
  double fpot_tor875good;
  double fspill_tor875;
  double fspill_tor875good;
   

  // output-- histograms
  // Selection
 
  int fhist_nbins;
  float fhist_binlow;
  float fhist_binup;
 
  float fcut_lowenergy;
  float fcut_completeness;

  // not FC
  TH1F* h1_total;
  TH1F* h1_cosmic;
  TH1F* h1_numuCC_FV;
  TH1F* h1_numuCC_nFV;
  TH1F* h1_numuNC_FV;
  TH1F* h1_numuNC_nFV;
  TH1F* h1_numubar;
  TH1F* h1_nue;
  TH1F* h1_nuebar;
  // FC
  TH1F* h1_FC_total;
  TH1F* h1_FC_cosmic;
  TH1F* h1_FC_numuCC_FV;
  TH1F* h1_FC_numuCC_nFV;
  TH1F* h1_FC_numuNC_FV;
  TH1F* h1_FC_numuNC_nFV;
  TH1F* h1_FC_numubar;
  TH1F* h1_FC_nue;
  TH1F* h1_FC_nuebar;

  // breakdown of inefficiency
  TH1F* h2_nuCC_FV; 
  TH1F* h2_nuCC_FCselected;

  TH1F* h2_flash_found;
  TH1F* h2_match_found;
  TH1F* h2_match_energy;
  TH1F* h2_lightmismatch; 
  TH1F* h2_Tgm;
  TH1F* h2_FC;
  TH1F* h2_cosmic_match;  

  TH1F* h2_nuNC_FV; 
  TH1F* h2_nuNC_FCselected;

  TH1F* h2_NC_flash_found;
  TH1F* h2_NC_match_found;
  TH1F* h2_NC_match_energy;
  TH1F* h2_NC_lightmismatch; 
  TH1F* h2_NC_Tgm;
  TH1F* h2_NC_FC;
  TH1F* h2_NC_cosmic_match;  
};


WCPselection::WCPselection(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  
  // fcl config
  reconfigure(p);

  // T_eval / event
  // Histograms / event
  // T_pot / subrun
  initOutput();

}

void WCPselection::reconfigure(fhicl::ParameterSet const& pset)
{
  std::cout<<"------------ WCPselection::reconfigure ----------"<<std::endl;  

  fContainmentLabel = pset.get<std::string>("ContainmentLabel");
  fChargeLabel = pset.get<std::string>("ChargeLabel");
  fTruthLabel = pset.get<std::string>("TruthLabel");
  fMatchLabel = pset.get<std::string>("MatchLabel");
  fMC = pset.get<bool>("MC"); // overlay and full mc

  fhist_nbins = pset.get<int>("Hist_nbins");
  fhist_binlow = pset.get<float>("Hist_binlow");
  fhist_binup = pset.get<float>("Hist_binup");
  fcut_lowenergy = pset.get<float>("Cut_lowE");
  fcut_completeness = pset.get<float>("Cut_completeness");

  fPOT_inputTag = pset.get<std::string>("POT_inputTag");
  fPOT_counting = pset.get<bool>("POT_counting");
 
}

void WCPselection::initOutput()
{
  std::cout<<"------------ WCPselection::initOutput ----------"<<std::endl;  

  art::ServiceHandle<art::TFileService> tfs;
  fTreeEval = tfs->make<TTree>("T_eval", "T_eval");
  
  fTreeEval->Branch("run", 			&f_run);
  fTreeEval->Branch("subrun", 			&f_subRun);
  fTreeEval->Branch("event", 			&f_event);
  fTreeEval->Branch("flash_found", 		&f_flash_found);
  fTreeEval->Branch("flash_time", 		&f_flash_time);
  fTreeEval->Branch("flash_measPe", 		&f_flash_measPe);
  fTreeEval->Branch("flash_predPe", 		&f_flash_predPe);
  fTreeEval->Branch("match_found", 		&f_match_found);
  fTreeEval->Branch("match_type", 		&f_match_type);
  fTreeEval->Branch("match_isFC", 		&f_match_isFC);
  fTreeEval->Branch("match_isTgm", 		&f_match_isTgm);
  fTreeEval->Branch("match_notFC_FV", 		&f_match_notFC_FV);
  fTreeEval->Branch("match_notFC_SP", 		&f_match_notFC_SP);
  fTreeEval->Branch("match_notFC_DC", 		&f_match_notFC_DC);
  fTreeEval->Branch("match_chargeU", 		&f_match_chargeU);
  fTreeEval->Branch("match_chargeV", 		&f_match_chargeV);
  fTreeEval->Branch("match_chargeY", 		&f_match_chargeY);
  fTreeEval->Branch("match_energyY", 		&f_match_energyY);
  fTreeEval->Branch("light_mismatch", 		&f_lightmismatch);
  fTreeEval->Branch("match_charge", 		&f_match_charge);
  fTreeEval->Branch("match_energy", 		&f_match_energy);
  if( fMC==true ){
  fTreeEval->Branch("truth_nuEnergy", 		&f_truth_nuEnergy);
  fTreeEval->Branch("truth_energyInside", 	&f_truth_energyInside);
  fTreeEval->Branch("truth_electronInside", 	&f_truth_electronInside);
  fTreeEval->Branch("truth_nuPdg", 		&f_truth_nuPdg);
  fTreeEval->Branch("truth_isCC", 		&f_truth_isCC);
  fTreeEval->Branch("truth_isEligible", 	&f_truth_isEligible);
  fTreeEval->Branch("truth_NCisEligible", 	&f_NC_truth_isEligible);
  fTreeEval->Branch("truth_isFC", 		&f_truth_isFC);
  fTreeEval->Branch("truth_vtxInside", 		&f_truth_vtxInside);
  fTreeEval->Branch("truth_vtxX", 		&f_truth_vtxX);
  fTreeEval->Branch("truth_vtxY", 		&f_truth_vtxY);
  fTreeEval->Branch("truth_vtxZ", 		&f_truth_vtxZ);
  fTreeEval->Branch("truth_nuTime", 		&f_truth_nuTime);
  fTreeEval->Branch("match_completeness", 	&f_match_completeness);
  fTreeEval->Branch("match_completeness_energy",&f_match_completeness_energy);
  fTreeEval->Branch("match_purity", 		&f_match_purity);
  fTreeEval->Branch("match_purity_xz", 		&f_match_purity_xz);
  fTreeEval->Branch("match_purity_xy", 		&f_match_purity_xy);
  }
  fTreePot = tfs->make<TTree>("T_pot", "T_pot");
  fTreePot->Branch("runNo", &f_run);
  fTreePot->Branch("subRunNo", &f_subRun);
  fTreePot->Branch("pot_tor875", &fpot_tor875);
  fTreePot->Branch("pot_tor875good", &fpot_tor875good);
  fTreePot->Branch("spill_tor875", &fspill_tor875);
  fTreePot->Branch("spill_tor875good", &fspill_tor875good);
 


  int nbins = fhist_nbins; // 33
  float bin_low = fhist_binlow; // 0
  float bin_up = fhist_binup; // 1980 MeV
 
  TString roostr = "";
  //art::TFileDirectory hist = tfs->mkdir("hist");
  if( fMC==true ){
  roostr="h1_cosmic"; 
  h1_cosmic = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up); 
  roostr="h1_numuCC_FV"; 
  h1_numuCC_FV = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up); 
  roostr="h1_numuCC_nFV"; 
  h1_numuCC_nFV = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up); 
  roostr="h1_numuNC_FV"; 
  h1_numuNC_FV = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up); 
  roostr="h1_numuNC_nFV"; 
  h1_numuNC_nFV = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up); 
  roostr="h1_numubar"; 
  h1_numubar = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up); 
  roostr="h1_nue"; 
  h1_nue = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up); 
  roostr="h1_nuebar"; 
  h1_nuebar = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up); 
 
  roostr="h1_FC_cosmic"; 
  h1_FC_cosmic = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up); 
  roostr="h1_FC_numuCC_FV"; 
  h1_FC_numuCC_FV = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up); 
  roostr="h1_FC_numuCC_nFV"; 
  h1_FC_numuCC_nFV = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up); 
  roostr="h1_FC_numuNC_FV"; 
  h1_FC_numuNC_FV = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up); 
  roostr="h1_FC_numuNC_nFV"; 
  h1_FC_numuNC_nFV = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up); 
  roostr="h1_FC_numubar"; 
  h1_FC_numubar = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up); 
  roostr="h1_FC_nue"; 
  h1_FC_nue = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up); 
  roostr="h1_FC_nuebar"; 
  h1_FC_nuebar = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up); 

  roostr="h2_nuCC_FV"; 
  h2_nuCC_FV = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);
  roostr="h2_nuCC_FCselected";
  h2_nuCC_FCselected = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);
  roostr="h2_flash_found";
  h2_flash_found = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);
  roostr="h2_match_found";
  h2_match_found = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);
  roostr="h2_match_energy";
  h2_match_energy = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);
  roostr="h2_lightmismatch"; 
  h2_lightmismatch = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up); 
  roostr="h2_Tgm";
  h2_Tgm = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);
  roostr="h2_FC";
  h2_FC = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);
  roostr="h2_cosmic_match";  
  h2_cosmic_match = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);  
  
  roostr="h2_nuNC_FV"; 
  h2_nuNC_FV = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);
  roostr="h2_nuNC_FCselected";
  h2_nuNC_FCselected = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);
  roostr="h2_NC_flash_found";
  h2_NC_flash_found = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);
  roostr="h2_NC_match_found";
  h2_NC_match_found = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);
  roostr="h2_NC_match_energy";
  h2_NC_match_energy = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);
  roostr="h2_NC_lightmismatch"; 
  h2_NC_lightmismatch = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up); 
  roostr="h2_NC_Tgm";
  h2_NC_Tgm = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);
  roostr="h2_NC_FC";
  h2_NC_FC = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);
  roostr="h2_NC_cosmic_match";  
  h2_NC_cosmic_match = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);  
  }

  roostr="h1_total";
  h1_total = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);
  roostr="h1_FC_total";
  h1_FC_total = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);

}

void WCPselection::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  // read NuMetrics
	std::cout<<" RUN: "<<e.run()<<"\n SUBRUN: "<<e.subRun()<<"\n EVENT: "<<e.id().event()<<std::endl;
	f_run = e.run();
 	f_subRun = e.subRun();
	f_event = e.id().event();

	art::Handle<std::vector<nsm::NuSelectionContainment> > containment_handle;
	e.getByLabel(fContainmentLabel,containment_handle);
	std::vector<art::Ptr<nsm::NuSelectionContainment> > containment_vec;
	art::fill_ptr_vector(containment_vec,containment_handle);
	std::cout<<"--- NuSelectionContainment ---"<<std::endl;
	if(containment_vec.size()!=1) {
		std::cout<<"WARNING: >1 in-beam matched TPC activity?!" << std::endl;
		return;
	} 
	for(size_t i=0; i<containment_vec.size(); i++){
		art::Ptr<nsm::NuSelectionContainment> c = containment_vec.at(i);
		f_flash_found = c->GetFlashFound();
		f_flash_time = c->GetFlashTime();
		f_flash_measPe = c->GetFlashMeasPe();
		f_flash_predPe = c->GetFlashPredPe();
		f_match_found = c->GetMatchFound();
		f_match_type = c->GetMatchType();
		f_match_isFC = c->GetIsFC();
		f_match_isTgm = c->GetIsTGM();
		f_match_notFC_FV = c->GetNotFCFV();
		f_match_notFC_SP = c->GetNotFCSP();
		f_match_notFC_DC = c->GetNotFCDC();
		f_match_charge = c->GetCharge();
		f_match_energy = c->GetEnergy();
	}

	art::Handle<std::vector<nsm::NuSelectionCharge> > charge_handle;
	e.getByLabel(fChargeLabel,charge_handle);
	std::vector<art::Ptr<nsm::NuSelectionCharge> > charge_vec;
	art::fill_ptr_vector(charge_vec,charge_handle);
	std::cout<<"--- NuSelectionCharge  ---"<<std::endl;
	if(charge_vec.size()!=1) {
		std::cout<<"WARNING: >1 in-beam matched TPC activity?!" << std::endl;
		return;
	} 
	for(size_t i=0; i<charge_vec.size(); i++){
		art::Ptr<nsm::NuSelectionCharge> c = charge_vec.at(i);
		f_match_chargeU = c->GetChargeU();
		f_match_chargeV = c->GetChargeV();
		f_match_chargeY = c->GetChargeY();
	}

	if(fMC==true){

	art::Handle<std::vector<nsm::NuSelectionTruth> > truth_handle;
	e.getByLabel(fTruthLabel,truth_handle);
	std::vector<art::Ptr<nsm::NuSelectionTruth> > truth_vec;
	art::fill_ptr_vector(truth_vec,truth_handle);
	std::cout<<"--- NuSelectionTruth  ---"<<std::endl;
	if(truth_vec.size()!=1) {
		std::cout<<"WARNING: >1 in-beam matched TPC activity?!" << std::endl;
		return;
	} 
	for(size_t i=0; i<truth_vec.size(); i++){
		art::Ptr<nsm::NuSelectionTruth> t = truth_vec.at(i);
		f_truth_nuEnergy = t->GetNuEnergy();
		f_truth_energyInside = t->GetEnergyInside();
		f_truth_electronInside = t->GetElectronInside();
		f_truth_nuPdg = t->GetNuPdg();
		f_truth_isCC = t->GetIsCC();
		f_truth_isEligible = t->GetIsEligible();
		f_truth_isFC = t->GetIsFC();
		f_truth_vtxInside = t->GetIsVtxInside();
		f_truth_vtxX = t->GetVtxX();
		f_truth_vtxY = t->GetVtxY();
		f_truth_vtxZ = t->GetVtxZ();
		f_truth_nuTime = t->GetTime();
	}

	art::Handle<std::vector<nsm::NuSelectionMatch> > match_handle;
	e.getByLabel(fMatchLabel,match_handle);
	std::vector<art::Ptr<nsm::NuSelectionMatch> > match_vec;
	art::fill_ptr_vector(match_vec,match_handle);
	std::cout<<"--- NuSelectionMatch  ---"<<std::endl;
	if(match_vec.size()!=1) {
		std::cout<<"WARNING: >1 in-beam matched TPC activity?!" << std::endl;
		return;
	} 
	for(size_t i=0; i<match_vec.size(); i++){
		art::Ptr<nsm::NuSelectionMatch> m = match_vec.at(i);
		f_match_completeness = m->GetCompleteness();
		f_match_completeness_energy = m->GetCompletenessEnergy();
		f_match_purity = m->GetPurity();
		f_match_purity_xz = m->GetPurityXZ();
		f_match_purity_xy = m->GetPurityXY();
	}
	
	if ( !f_truth_isCC && f_truth_vtxInside && f_truth_nuPdg==14 ) f_NC_truth_isEligible = true;
	else f_NC_truth_isEligible = false;	

	}
	


	if ( fMC==false ) f_truth_isEligible = false;
	if ( fMC==false ) f_NC_truth_isEligible = false;
	f_match_energyY = f_match_chargeY*23.6*1e-6/0.55;
	if( f_match_type&1U || (f_match_type>>1)&1U ) f_lightmismatch = true;
	else f_lightmismatch = false;

	// Do selection & Fill histograms
	if( f_truth_isEligible ) h2_nuCC_FV->Fill(f_truth_energyInside);
	if( f_NC_truth_isEligible ) h2_nuNC_FV->Fill(f_truth_energyInside);
	if( !f_flash_found ) {
		if(f_truth_isEligible) h2_flash_found->Fill(f_truth_energyInside);	
		if(f_NC_truth_isEligible) h2_NC_flash_found->Fill(f_truth_energyInside);	
	}
	else if( !f_match_found) {
		if(f_truth_isEligible) h2_match_found->Fill(f_truth_energyInside);		
		if(f_NC_truth_isEligible) h2_NC_match_found->Fill(f_truth_energyInside);		
	}
	else if( f_match_energyY<fcut_lowenergy ){
		if(f_truth_isEligible) h2_match_energy->Fill(f_truth_energyInside);		
		if(f_NC_truth_isEligible) h2_NC_match_energy->Fill(f_truth_energyInside);		
	}
	else if( f_lightmismatch ){
		if(f_truth_isEligible) h2_lightmismatch->Fill(f_truth_energyInside);		
		if(f_NC_truth_isEligible) h2_NC_lightmismatch->Fill(f_truth_energyInside);		
	}
	else if( f_match_isTgm ){
		if(f_truth_isEligible) h2_Tgm->Fill(f_truth_energyInside);		
		if(f_NC_truth_isEligible) h2_NC_Tgm->Fill(f_truth_energyInside);		
	}
	else if( !f_match_isFC ){
		if(f_truth_isEligible) h2_FC->Fill(f_truth_energyInside);
		if(f_NC_truth_isEligible) h2_NC_FC->Fill(f_truth_energyInside);
		
		h1_total->Fill(f_match_energyY);
		if( fMC && f_match_completeness_energy/f_truth_energyInside>fcut_completeness ){
			if( f_truth_isCC && f_truth_nuPdg==14 && f_truth_vtxInside) h1_numuCC_FV->Fill(f_match_energyY);
			if( f_truth_isCC && f_truth_nuPdg==14 && !f_truth_vtxInside) h1_numuCC_nFV->Fill(f_match_energyY);		
			if( !f_truth_isCC && f_truth_nuPdg==14 && f_truth_vtxInside) h1_numuNC_FV->Fill(f_match_energyY);
			if( !f_truth_isCC && f_truth_nuPdg==14 && !f_truth_vtxInside) h1_numuNC_nFV->Fill(f_match_energyY);		
			if( f_truth_nuPdg==-14) h1_numubar->Fill(f_match_energyY);		
			if( f_truth_nuPdg==12) h1_nue->Fill(f_match_energyY);		
			if( f_truth_nuPdg==-12) h1_nuebar->Fill(f_match_energyY);		
		}
		else if( fMC ){
			h1_cosmic->Fill(f_match_energyY);
		}	
	
	}
	else{
		if( f_match_completeness_energy/f_truth_energyInside<fcut_completeness ){
			if(f_truth_isEligible) h2_cosmic_match->Fill(f_truth_energyInside);
			if(f_NC_truth_isEligible) h2_NC_cosmic_match->Fill(f_truth_energyInside);
		}
		else{
			if(f_truth_isEligible) h2_nuCC_FCselected->Fill(f_truth_energyInside);
			if(f_NC_truth_isEligible) h2_nuNC_FCselected->Fill(f_truth_energyInside);
		}
		
		h1_FC_total->Fill(f_match_energyY);
		if( fMC && f_match_completeness_energy/f_truth_energyInside>fcut_completeness ){
			if( f_truth_isCC && f_truth_nuPdg==14 && f_truth_vtxInside) h1_FC_numuCC_FV->Fill(f_match_energyY);
			if( f_truth_isCC && f_truth_nuPdg==14 && !f_truth_vtxInside) h1_FC_numuCC_nFV->Fill(f_match_energyY);		
			if( !f_truth_isCC && f_truth_nuPdg==14 && f_truth_vtxInside) h1_FC_numuNC_FV->Fill(f_match_energyY);
			if( !f_truth_isCC && f_truth_nuPdg==14 && !f_truth_vtxInside) h1_FC_numuNC_nFV->Fill(f_match_energyY);		
			if( f_truth_nuPdg==-14) h1_FC_numubar->Fill(f_match_energyY);		
			if( f_truth_nuPdg==12) h1_FC_nue->Fill(f_match_energyY);		
			if( f_truth_nuPdg==-12) h1_FC_nuebar->Fill(f_match_energyY);		
		}
		else if( fMC ){
			h1_FC_cosmic->Fill(f_match_energyY);
		}	
	}

	fTreeEval->Fill();
}

void WCPselection::endSubRun(art::SubRun const& sr)
{
  // Implementation of optional member function here.

  // POT counting
  if( fPOT_counting==true ){
	art::Handle<sumdata::POTSummary> pots;
	if(! sr.getByLabel(fPOT_inputTag, pots)){
		std::cout << "WARNING:  no sumdata::POTSummary inputTag " << fPOT_inputTag << std::endl;
		return;
	}
	sumdata::POTSummary const& p1(*pots);
	fpot_tor875 = p1.totpot;
	fpot_tor875good = p1.totgoodpot;
	fspill_tor875 = p1.totspills;
	fspill_tor875good = p1.goodspills;
  }
  fTreePot->Fill(); 
  // check if MC POT is cumulative for each file as MCC8 overlay
  // sumdata::POTSummary used for Overlay and ?MC
  // Zarko's script used for data (bnb, ext)
  // If fPOT_counting is false, we still need run and subrun info
}

DEFINE_ART_MODULE(WCPselection)
