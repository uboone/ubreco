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

#include "art_root_io/TFileService.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h" 
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"

#include "ubobj/WcpPort/NuSelectionContainment.h"
#include "ubobj/WcpPort/NuSelectionMatch.h"
#include "ubobj/WcpPort/NuSelectionTruth.h"
#include "ubobj/WcpPort/NuSelectionCharge.h"
#include "ubobj/WcpPort/NuSelectionSTM.h"
#include "ubobj/WcpPort/NuSelectionBDT.h"


#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TVector3.h"

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
  void resetOutput();
  void ShowerID(int trackId);
  void MuonID(int trackId);
  void save_weights(art::Event const& e);

private:

  // Declare member data here.

  // fcl config
  std::string fContainmentLabel;
  std::string fChargeLabel;
  std::string fTruthLabel;
  std::string fMatchLabel;
  std::string fSTMLabel;
  bool fMC;
  bool f_wirecellPF;
  bool fSaveWeights;
  bool f_BDTvars;

  bool fPFValidation; // switch of particle flow validation
  std::string fPFInputTag; // inputTag -- label:instance:process
  std::string fPFtruthInputTag; // inputTag -- label:instance:process
  float fthreshold_showerKE;
  std::vector<int> fPrimaryID;
  std::vector<int> fMuonID;
  std::vector<int> fShowerID;
  std::map<int, simb::MCParticle> fParticleMap; // map from trackId to particle instance
 
  // output
  /// PF validation 
  /// when fPFValidation is true
  TTree* fPFeval;
  Int_t		f_neutrino_type;
  Float_t	f_reco_nuvtxX;
  Float_t	f_reco_nuvtxY;
  Float_t	f_reco_nuvtxZ;
  Float_t	f_reco_showervtxX; // primary shower [highest energy electron & not pi0 daughters] 
  Float_t	f_reco_showervtxY;
  Float_t	f_reco_showervtxZ;
  Float_t 	f_reco_showerKE;
  Float_t	f_reco_muonvtxX; // 
  Float_t	f_reco_muonvtxY;
  Float_t	f_reco_muonvtxZ;
  Float_t	f_reco_muonMomentum[4];
  Float_t	f_truth_corr_nuvtxX; // truth -(SCE)-> SED -(nu time offset)-> reco [trigger offset O(10) ns ignored]
  Float_t	f_truth_corr_nuvtxY;
  Float_t	f_truth_corr_nuvtxZ;
  Float_t	f_truth_corr_showervtxX; // truth -(SCE)-> SED -(nu time offset)-> reco [trigger offset O(10) ns ignored]
  Float_t	f_truth_corr_showervtxY;
  Float_t	f_truth_corr_showervtxZ;
  Float_t 	f_truth_showerKE;
  Float_t	f_truth_corr_muonvtxX; // primary muon, vtx = nu vtx 
  Float_t	f_truth_corr_muonvtxY;
  Float_t	f_truth_corr_muonvtxZ;
  Float_t	f_truth_muonvtxX; 
  Float_t	f_truth_muonvtxY;
  Float_t	f_truth_muonvtxZ;
  Float_t	f_truth_muonendX; // may not in TPC active 
  Float_t	f_truth_muonendY;
  Float_t	f_truth_muonendZ;
  Float_t	f_truth_muonMomentum[4];
  Int_t		f_truth_nuIntType;
  /// other truth info as follows save in this tree

  /// BDT input vars
  TTree* fBDT;
  Int_t M_tag;
  Float_t M_shower_energy; // MeV
  Float_t M_shower_angle; // degree, relative to beam
  Float_t M_seg_max_length; // cm 
  Float_t M_seg_max_angle; // degree, relative to shower
  Float_t M_acc_forward_length; // cm
  Float_t M_acc_backward_length; // cm
  Int_t M_flag_shwtrunck_outside;
  Int_t E_tag;
  Float_t E_shower_energy; // MeV
  Float_t E_shower_length; // cm
  Float_t E_shower_angle; // degree, relative to beam
  Float_t E_muon_energy; // MeV
  Float_t E_muon_length; // cm
  Int_t E_flag_numuCC;
  Int_t flag_others;
  ///

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
  
  Float_t 	  f_weight_spline;
  Float_t	  f_weight_cv;
 
  Int_t		  f_stm_eventtype;
  Int_t		  f_stm_lowenergy;
  Int_t		  f_stm_LM;
  Int_t		  f_stm_TGM;
  Int_t		  f_stm_STM;
  Int_t		  f_stm_FullDead;
  Float_t	  f_stm_clusterlength;

  TTree* fTreePot;
  std::string fPOT_inputTag; 
  bool fPOT_counting=false; 
 // beamdata:bnbETOR875
 // EXT?
 // MC [label: generator] : [instance: ]

  //int frun;
  //int fsubRun;  
  double fpot_tor875;
  double fpot_tor875good;
  double fspill_tor875;
  double fspill_tor875good;
 

  // output-- histograms
  // Selection
 
  int fhist_nbins;
  float fhist_binlow;
  float fhist_binup;
 
  float fcut_lowenergy; //MeV
  float fcut_completeness;
  float fcut_clusterlength; //cm
  /*
  //obselete
  // not STM
  TH1F* h1_total;
  TH1F* h1_cosmic;
  TH1F* h1_numuCC_FV;
  TH1F* h1_numuCC_nFV;
  TH1F* h1_numuNC_FV;
  TH1F* h1_numuNC_nFV;
  TH1F* h1_numubar;
  TH1F* h1_nue;
  TH1F* h1_nuebar;
  // STM
  TH1F* h1_STM_total;
  TH1F* h1_STM_cosmic;
  TH1F* h1_STM_numuCC_FV;
  TH1F* h1_STM_numuCC_nFV;
  TH1F* h1_STM_numuNC_FV;
  TH1F* h1_STM_numuNC_nFV;
  TH1F* h1_STM_numubar;
  TH1F* h1_STM_nue;
  TH1F* h1_STM_nuebar;

  // breakdown of inefficiency
  TH1F* h2_nuCC_FV; 
  TH1F* h2_nuCC_STMselected;

  TH1F* h2_flash_found;
  TH1F* h2_match_found;
  TH1F* h2_match_energy;
  TH1F* h2_lightmismatch; 
  TH1F* h2_Tgm;
  TH1F* h2_STM;
  TH1F* h2_cosmic_match;  

  TH1F* h2_nuNC_FV; 
  TH1F* h2_nuNC_STMselected;

  TH1F* h2_NC_flash_found;
  TH1F* h2_NC_match_found;
  TH1F* h2_NC_match_energy;
  TH1F* h2_NC_lightmismatch; 
  TH1F* h2_NC_Tgm;
  TH1F* h2_NC_STM;
  TH1F* h2_NC_cosmic_match; 
  */ 
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
  f_wirecellPF = pset.get<bool>("wirecellPF", false);
  f_BDTvars = pset.get<bool>("BDTvars", false);
  fSaveWeights = pset.get<bool>("SaveWeights", false); // GENIE weights
  fSTMLabel = pset.get<std::string>("STMLabel");
  
  fhist_nbins = pset.get<int>("Hist_nbins");
  fhist_binlow = pset.get<float>("Hist_binlow");
  fhist_binup = pset.get<float>("Hist_binup");
  fcut_lowenergy = pset.get<float>("Cut_lowE");
  fcut_completeness = pset.get<float>("Cut_completeness");
  fcut_clusterlength = pset.get<float>("Cut_clusterlength");

  fPOT_inputTag = pset.get<std::string>("POT_inputTag");
  fPOT_counting = pset.get<bool>("POT_counting");

  fPFValidation = pset.get<bool>("PF_validation");
  fPFInputTag = pset.get<std::string>("PF_inputtag");
  fPFtruthInputTag = pset.get<std::string>("PFtruth_inputtag");
  fthreshold_showerKE = pset.get<float>("Threshold_showerKE"); // GeV 
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
  fTreeEval->Branch("stm_eventtype",		&f_stm_eventtype);
  fTreeEval->Branch("stm_lowenergy",		&f_stm_lowenergy);
  fTreeEval->Branch("stm_LM",			&f_stm_LM);
  fTreeEval->Branch("stm_TGM",			&f_stm_TGM);
  fTreeEval->Branch("stm_STM",			&f_stm_STM);
  fTreeEval->Branch("stm_FullDead",		&f_stm_FullDead);
  fTreeEval->Branch("stm_clusterlength",	&f_stm_clusterlength);

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

  fTreeEval->Branch("weight_spline", 		&f_weight_spline); //MicroBooNE GENIE tune on top of weight_CV; weight_spline*weight_cv = weight
  fTreeEval->Branch("weight_cv",		&f_weight_cv); //MicroBooNE MCC9 untuned GENIE v3

  }
  fTreePot = tfs->make<TTree>("T_pot", "T_pot");
  fTreePot->Branch("runNo", &f_run);
  fTreePot->Branch("subRunNo", &f_subRun);
  fTreePot->Branch("pot_tor875", &fpot_tor875);
  fTreePot->Branch("pot_tor875good", &fpot_tor875good);
  fTreePot->Branch("spill_tor875", &fspill_tor875);
  fTreePot->Branch("spill_tor875good", &fspill_tor875good);
 
  /// PF validation 
  fPFeval = tfs->make<TTree>("T_PFeval", "T_PFeval");
  fPFeval->Branch("run", 			&f_run);
  fPFeval->Branch("subrun", 			&f_subRun);
  fPFeval->Branch("event", 			&f_event);
  fPFeval->Branch("neutrino_type", 		&f_neutrino_type);
  fPFeval->Branch("reco_nuvtxX", 		&f_reco_nuvtxX);
  fPFeval->Branch("reco_nuvtxY", 		&f_reco_nuvtxY);
  fPFeval->Branch("reco_nuvtxZ", 		&f_reco_nuvtxZ);
  fPFeval->Branch("reco_showervtxX", 		&f_reco_showervtxX);
  fPFeval->Branch("reco_showervtxY", 		&f_reco_showervtxY);
  fPFeval->Branch("reco_showervtxZ", 		&f_reco_showervtxZ);
  fPFeval->Branch("reco_showerKE", 		&f_reco_showerKE);
  fPFeval->Branch("reco_muonvtxX", 		&f_reco_muonvtxX);
  fPFeval->Branch("reco_muonvtxY", 		&f_reco_muonvtxY);
  fPFeval->Branch("reco_muonvtxZ", 		&f_reco_muonvtxZ);
  fPFeval->Branch("reco_muonMomentum", 		&f_reco_muonMomentum, "reco_muonMomentum[4]/F");
  if( fMC==true ){
  fPFeval->Branch("truth_corr_nuvtxX", 		&f_truth_corr_nuvtxX);
  fPFeval->Branch("truth_corr_nuvtxY", 		&f_truth_corr_nuvtxY);
  fPFeval->Branch("truth_corr_nuvtxZ", 		&f_truth_corr_nuvtxZ);
  fPFeval->Branch("truth_corr_showervtxX", 	&f_truth_corr_showervtxX);
  fPFeval->Branch("truth_corr_showervtxY", 	&f_truth_corr_showervtxY);
  fPFeval->Branch("truth_corr_showervtxZ", 	&f_truth_corr_showervtxZ);
  fPFeval->Branch("truth_showerKE", 		&f_truth_showerKE);
  fPFeval->Branch("truth_corr_muonvtxX",	&f_truth_corr_muonvtxX);
  fPFeval->Branch("truth_corr_muonvtxY", 	&f_truth_corr_muonvtxY);
  fPFeval->Branch("truth_corr_muonvtxZ", 	&f_truth_corr_muonvtxZ);
  fPFeval->Branch("truth_muonvtxX",		&f_truth_muonvtxX);
  fPFeval->Branch("truth_muonvtxY", 		&f_truth_muonvtxY);
  fPFeval->Branch("truth_muonvtxZ", 		&f_truth_muonvtxZ);
  fPFeval->Branch("truth_muonendX",		&f_truth_muonendX);
  fPFeval->Branch("truth_muonendY", 		&f_truth_muonendY);
  fPFeval->Branch("truth_muonendZ", 		&f_truth_muonendZ);
  fPFeval->Branch("truth_muonMomentum", 	&f_truth_muonMomentum, "truth_muonMomentum[4]/F");
  fPFeval->Branch("truth_nuEnergy", 		&f_truth_nuEnergy);
  fPFeval->Branch("truth_energyInside", 	&f_truth_energyInside);
  fPFeval->Branch("truth_electronInside", 	&f_truth_electronInside);
  fPFeval->Branch("truth_nuPdg", 		&f_truth_nuPdg);
  fPFeval->Branch("truth_isCC", 		&f_truth_isCC);
  fPFeval->Branch("truth_vtxX", 		&f_truth_vtxX);
  fPFeval->Branch("truth_vtxY", 		&f_truth_vtxY);
  fPFeval->Branch("truth_vtxZ", 		&f_truth_vtxZ);
  fPFeval->Branch("truth_nuTime", 		&f_truth_nuTime);
  fPFeval->Branch("truth_nuIntType",		&f_truth_nuIntType);
  }

  fBDT = tfs->make<TTree>("T_BDTvars", "T_BDTvars");
  fBDT->Branch("M_tag",				&M_tag); 
  fBDT->Branch("M_shower_energy",		&M_shower_energy); 
  fBDT->Branch("M_shower_angle",		&M_shower_angle);
  fBDT->Branch("M_seg_max_length",		&M_seg_max_length);
  fBDT->Branch("M_seg_max_angle",		&M_seg_max_angle);
  fBDT->Branch("M_acc_forward_length",		&M_acc_forward_length);
  fBDT->Branch("M_acc_backward_length",		&M_acc_backward_length);
  fBDT->Branch("M_flag_shwtrunck_outside",	&M_flag_shwtrunck_outside);
  fBDT->Branch("E_tag",				&E_tag);
  fBDT->Branch("E_shower_energy",		&E_shower_energy);
  fBDT->Branch("E_shower_length",		&E_shower_length);
  fBDT->Branch("E_shower_angle",		&E_shower_angle);
  fBDT->Branch("E_muon_energy",			&E_muon_energy);
  fBDT->Branch("E_muon_length",			&E_muon_length);
  fBDT->Branch("E_flag_numuCC",			&E_flag_numuCC);
  fBDT->Branch("flag_others",			&flag_others); 

 
  /// 

  //obselete
  /*
  int nbins = fhist_nbins; // 33
  float bin_low = fhist_binlow; // 0
  float bin_up = fhist_binup; // 1980 MeV
 
  TString roostr = "";
  //art::TFileDirectory hist = tfs->mkdir("hist");
  if( fMC ){
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
 
  roostr="h1_STM_cosmic"; 
  h1_STM_cosmic = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up); 
  roostr="h1_STM_numuCC_FV"; 
  h1_STM_numuCC_FV = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up); 
  roostr="h1_STM_numuCC_nFV"; 
  h1_STM_numuCC_nFV = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up); 
  roostr="h1_STM_numuNC_FV"; 
  h1_STM_numuNC_FV = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up); 
  roostr="h1_STM_numuNC_nFV"; 
  h1_STM_numuNC_nFV = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up); 
  roostr="h1_STM_numubar"; 
  h1_STM_numubar = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up); 
  roostr="h1_STM_nue"; 
  h1_STM_nue = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up); 
  roostr="h1_STM_nuebar"; 
  h1_STM_nuebar = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up); 

  roostr="h2_nuCC_FV"; 
  h2_nuCC_FV = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);
  roostr="h2_nuCC_STMselected";
  h2_nuCC_STMselected = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);
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
  roostr="h2_STM";
  h2_STM = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);
  roostr="h2_cosmic_match";  
  h2_cosmic_match = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);  
  
  roostr="h2_nuNC_FV"; 
  h2_nuNC_FV = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);
  roostr="h2_nuNC_STMselected";
  h2_nuNC_STMselected = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);
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
  roostr="h2_NC_STM";
  h2_NC_STM = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);
  roostr="h2_NC_cosmic_match";  
  h2_NC_cosmic_match = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);  
  }
  roostr="h1_total";
  h1_total = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);
  roostr="h1_STM_total";
  h1_STM_total = tfs->make<TH1F>(roostr, roostr, nbins, bin_low, bin_up);
  */
}

void WCPselection::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  // reset output at the beginning of each event
  resetOutput();

  // read NuMetrics
	std::cout<<" RUN: "<<e.run()<<"\n SUBRUN: "<<e.subRun()<<"\n EVENT: "<<e.id().event()<<std::endl;
	f_run = e.run();
 	f_subRun = e.subRun();
	f_event = e.id().event();

        auto const& containment_vec = e.getProduct<std::vector<nsm::NuSelectionContainment>>(fContainmentLabel);
	std::cout<<"--- NuSelectionContainment ---"<<std::endl;
	if(containment_vec.size()!=1) {
		std::cout<<"WARNING: >1 in-beam matched TPC activity?!" << std::endl;
		return;
	} 
        for(nsm::NuSelectionContainment const& c : containment_vec) {
                f_flash_found = c.GetFlashFound();
                f_flash_time = c.GetFlashTime();
                f_flash_measPe = c.GetFlashMeasPe();
                f_flash_predPe = c.GetFlashPredPe();
                f_match_found = c.GetMatchFound();
                f_match_type = c.GetMatchType();
                f_match_isFC = c.GetIsFC();
                f_match_isTgm = c.GetIsTGM();
                f_match_notFC_FV = c.GetNotFCFV();
                f_match_notFC_SP = c.GetNotFCSP();
                f_match_notFC_DC = c.GetNotFCDC();
                f_match_charge = c.GetCharge();
                f_match_energy = c.GetEnergy();
	}

        auto const& charge_vec = e.getProduct<std::vector<nsm::NuSelectionCharge>>(fChargeLabel);
	std::cout<<"--- NuSelectionCharge  ---"<<std::endl;
	if(charge_vec.size()!=1) {
		std::cout<<"WARNING: >1 in-beam matched TPC activity?!" << std::endl;
		return;
	} 
        for(nsm::NuSelectionCharge const& c : charge_vec) {
                f_match_chargeU = c.GetChargeU();
                f_match_chargeV = c.GetChargeV();
                f_match_chargeY = c.GetChargeY();
	}

        auto const& stm_vec = e.getProduct<std::vector<nsm::NuSelectionSTM>>(fSTMLabel);
	std::cout<<"--- NuSelectionSTM ---"<<std::endl;
	if(stm_vec.size()>1) {
		std::cout<<"WARNING: >1 in-beam matched TPC activity?!" << std::endl;
		return;
	} 
	if(stm_vec.size()<1) {
		f_stm_eventtype = -1;
		f_stm_lowenergy = -1;
		f_stm_LM = -1;
		f_stm_TGM = -1;
		f_stm_STM = -1;
		f_stm_FullDead = -1;
		f_stm_clusterlength = -1.0;
	} 
        for(nsm::NuSelectionSTM const& s : stm_vec) {
                f_stm_eventtype = s.GetEventType();
                f_stm_lowenergy = s.GetLowEnergy();
                f_stm_LM = s.GetLM();
                f_stm_TGM = s.GetTGM();
                f_stm_STM = s.GetSTM();
                f_stm_FullDead = s.GetFullDead();
                f_stm_clusterlength = s.GetClusterLength();
	}

	if(fMC==true){

        auto const& truth_vec = e.getProduct<std::vector<nsm::NuSelectionTruth>>(fTruthLabel);
	std::cout<<"--- NuSelectionTruth  ---"<<std::endl;
	if(truth_vec.size()!=1) {
		std::cout<<"WARNING: >1 in-beam matched TPC activity?!" << std::endl;
		return;
	} 
        for(nsm::NuSelectionTruth const& t : truth_vec) {
                f_truth_nuEnergy = t.GetNuEnergy();
                f_truth_energyInside = t.GetEnergyInside();
                f_truth_electronInside = t.GetElectronInside();
                f_truth_nuPdg = t.GetNuPdg();
                f_truth_isCC = t.GetIsCC();
                f_truth_isEligible = t.GetIsEligible();
                f_truth_isFC = t.GetIsFC();
                f_truth_vtxInside = t.GetIsVtxInside();
                f_truth_vtxX = t.GetVtxX();
                f_truth_vtxY = t.GetVtxY();
                f_truth_vtxZ = t.GetVtxZ();
                f_truth_nuTime = t.GetTime();
	}

        auto const& match_vec = e.getProduct<std::vector<nsm::NuSelectionMatch>>(fMatchLabel);
	std::cout<<"--- NuSelectionMatch  ---"<<std::endl;
	if(match_vec.size()!=1) {
		std::cout<<"WARNING: >1 in-beam matched TPC activity?!" << std::endl;
		return;
	} 
        for(nsm::NuSelectionMatch const& m : match_vec) {
                f_match_completeness = m.GetCompleteness();
                f_match_completeness_energy = m.GetCompletenessEnergy();
                f_match_purity = m.GetPurity();
                f_match_purity_xz = m.GetPurityXZ();
                f_match_purity_xy = m.GetPurityXY();
	}

	/// save GENIE weights
	if(fSaveWeights) save_weights(e);	
	/// end
	
	if ( !f_truth_isCC && f_truth_vtxInside && f_truth_nuPdg==14 ) f_NC_truth_isEligible = true;
	else f_NC_truth_isEligible = false;	

	}
	
	if ( fMC==false ) f_truth_isEligible = false;
	if ( fMC==false ) f_NC_truth_isEligible = false;
	f_match_energyY = f_match_chargeY*23.6*1e-6/0.55;
	if( f_match_type&1U || (f_match_type>>1)&1U ) f_lightmismatch = true;
	else f_lightmismatch = false;

	/// PF validation starts
	// reco start [nested loop]
      if( f_wirecellPF ){
        auto particleHandle = e.getHandle<std::vector<simb::MCParticle>>(fPFInputTag);
        if (! particleHandle) return;

        for (auto const& particle: *particleHandle) {
                int trkID = particle.TrackId();
                fParticleMap[trkID] = particle;
                if(particle.Mother() == 0){
			if(fPrimaryID.size()<1){ // fill once
                        const TLorentzVector& position = particle.Position(0);
			f_reco_nuvtxX = position.X(); // units: cm inherit from larsoft				
			f_reco_nuvtxY = position.Y(); // units: cm inherit from larsoft				
			f_reco_nuvtxZ = position.Z(); // units: cm inherit from larsoft		
                        f_neutrino_type = particle.StatusCode(); // neutrino type
			}
			fPrimaryID.push_back(trkID);
		}
	}
	for (size_t i=0; i<fPrimaryID.size(); i++){
		//std::cout<<"Primary particle:  "<< fPrimaryID.at(i) <<std::endl;
		MuonID(fPrimaryID.at(i));
		ShowerID(fPrimaryID.at(i)); // nested loop to dump descendant (not one generation) shower IDs with additional constraint
	}
	//muon
	float MuonKE = 0; //temp shower kinetic energy
	for (size_t j=0; j<fMuonID.size(); j++){
		auto const& p = fParticleMap[fMuonID.at(j)];
		const TLorentzVector& pos = p.Position(0);
		const TLorentzVector& momentum = p.Momentum(0);
		float tempKE = momentum.E() - momentum.M();
		if( MuonKE>=tempKE ) continue;
		MuonKE = tempKE;
		f_reco_muonvtxX = pos.X(); // cm	
		f_reco_muonvtxY = pos.Y();	
		f_reco_muonvtxZ = pos.Z();	
		f_reco_muonMomentum[0] = momentum.Px();	
		f_reco_muonMomentum[1] = momentum.Py();	
		f_reco_muonMomentum[2] = momentum.Pz();	
		f_reco_muonMomentum[3] = momentum.E(); // GeV
	}

 	//shower	
	float ShowerKE = 0; //temp shower kinetic energy
	bool ShowerPrimary = false; //temp shower primary particle? 
	for (size_t j=0; j<fShowerID.size(); j++){
		//std::cout<<"Shower: "<< fShowerID.at(j) <<std::endl;
		auto const& p = fParticleMap[fShowerID.at(j)];
		bool IsPrimary = false;
		// find primary shower (electron): primary particle > others; high energy > low energyi
		const TLorentzVector& pos = p.Position(0);
		TVector3 vnu(f_reco_nuvtxX, f_reco_nuvtxY, f_reco_nuvtxZ);
		TVector3 vshower(pos.X(), pos.Y(), pos.Z());
		// here the "primary" may be from a pi0 but connected to reco nu vertex
		if ( p.Mother()==0 || (vnu-vshower).Mag()<1 ) IsPrimary = true;
		const TLorentzVector& momentum = p.Momentum(0);
		float tempKE = momentum.E() - momentum.M();
		if( (ShowerKE>=tempKE || ShowerPrimary) && !IsPrimary ) continue; // logic bug: ignore coincidence equal shower KE particle ordering
		if( IsPrimary && ShowerPrimary && ShowerKE>=tempKE) continue;
		if( IsPrimary ) ShowerPrimary = true; // a primary shower recorded
		ShowerKE = tempKE;
		f_reco_showervtxX = pos.X(); // units: cm inherit from larsoft				
		f_reco_showervtxY = pos.Y(); // units: cm inherit from larsoft				
		f_reco_showervtxZ = pos.Z(); // units: cm inherit from larsoft		
		f_reco_showerKE = ShowerKE;
	}
	//std::cout<<"Primary shower KE: "<< ShowerKE << std::endl;
	// reco end

	if(fMC == true){ 	
	/// truth start
        auto particleHandle2 = e.getHandle<std::vector<simb::MCParticle>>(fPFtruthInputTag);
        if (! particleHandle2) return;

	// Get space charge correction
	auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
	
        for (auto const& particle: *particleHandle2){
                if( particle.Mother() == 0 && (particle.PdgCode() == 11 || particle.PdgCode() == -11) ){
                        const TLorentzVector& position = particle.Position(0);
			//f_truth_corr_showervtxX = position.X(); // units: cm inherit from larsoft				
			//f_truth_corr_showervtxY = position.Y(); // units: cm inherit from larsoft				
			//f_truth_corr_showervtxZ = position.Z(); // units: cm inherit from larsoft		
			auto sce_offset = SCE->GetPosOffsets(geo::Point_t(position.X(), position.Y(), position.Z()));
			f_truth_corr_showervtxX = position.X() - sce_offset.X();
			f_truth_corr_showervtxY = position.Y() + sce_offset.Y();
			f_truth_corr_showervtxZ = position.Z() + sce_offset.Z();
			f_truth_corr_showervtxX = (f_truth_corr_showervtxX + 0.6)*1.101/1.098 + position.T()*1e-3*1.101*0.1; //T: ns; 1.101 mm/us
			std::cout<<"Shower info: "<<position.X() <<", "<<position.Y() <<", "<<position.Z()<<", "<<position.T()<<" ns"<<std::endl;
			std::cout<<"Shower vertex SCE offset: "<<sce_offset.X() <<", "<<sce_offset.Y() <<", "<<sce_offset.Z()<<std::endl;
                        const TLorentzVector& showerMom = particle.Momentum(0);
			f_truth_showerKE = showerMom.E() - showerMom.M();
				
		}
                if( particle.Mother()==0 && (particle.PdgCode() == 13 || particle.PdgCode() == -13) ){
                        const TLorentzVector& position = particle.Position(0);
			f_truth_muonvtxX = position.X();
			f_truth_muonvtxY = position.Y();
			f_truth_muonvtxZ = position.Z();
			auto sce_offset = SCE->GetPosOffsets(geo::Point_t(position.X(), position.Y(), position.Z()));
			f_truth_corr_muonvtxX = position.X() - sce_offset.X();
			f_truth_corr_muonvtxY = position.Y() + sce_offset.Y();
			f_truth_corr_muonvtxZ = position.Z() + sce_offset.Z();
			f_truth_corr_muonvtxX = (f_truth_corr_muonvtxX + 0.6)*1.101/1.098 + position.T()*1e-3*1.101*0.1; //T: ns; 1.101 mm/us
			
                        const TLorentzVector& endposition = particle.EndPosition();
			f_truth_muonendX = endposition.X();
			f_truth_muonendY = endposition.Y();
			f_truth_muonendZ = endposition.Z();

                        const TLorentzVector& momentum = particle.Momentum(0);
			f_truth_muonMomentum[0] = momentum.Px();
			f_truth_muonMomentum[1] = momentum.Py();
			f_truth_muonMomentum[2] = momentum.Pz();
			f_truth_muonMomentum[3] = momentum.E();
		}

	}
	auto sce_offset = SCE->GetPosOffsets(geo::Point_t(f_truth_vtxX, f_truth_vtxY, f_truth_vtxZ));
	f_truth_corr_nuvtxX = f_truth_vtxX - sce_offset.X();
	f_truth_corr_nuvtxY = f_truth_vtxY + sce_offset.Y();
	f_truth_corr_nuvtxZ = f_truth_vtxZ + sce_offset.Z();
	f_truth_corr_nuvtxX = (f_truth_corr_nuvtxX + 0.6)*1.101/1.098 + f_truth_nuTime*1.101*0.1; //nuTime: us; 1.101 mm/us
	std::cout<<"Neutrino info: "<<f_truth_vtxX <<", "<<f_truth_vtxY <<", "<<f_truth_vtxZ<<", "<<f_truth_nuTime<<" us"<<std::endl;
	std::cout<<"Neutrino vertex SCE offset: "<<sce_offset.X() <<", "<<sce_offset.Y() <<", "<<sce_offset.Z()<<std::endl;
	
	// neutrino interaction type. Integer, see MCNeutrino.h for more details.
        auto const& mclist = e.getProduct<std::vector<simb::MCTruth>>("generator");

	if (mclist.size()>0) {
                simb::MCTruth const& mctruth = mclist.front();
                if (mctruth.NeutrinoSet()) {
                        simb::MCNeutrino nu = mctruth.GetNeutrino();
			f_truth_nuIntType = nu.InteractionType();
			// one can access more neutrino GENIE info
		}
	}
	}
	//
	
	/// truth end	
	std::cout<<"Corrected Truth Neutrino vertex: ("<<f_truth_corr_nuvtxX<<", "<<f_truth_corr_nuvtxY<<", "<<f_truth_corr_nuvtxZ<<")"<<"\n";
	std::cout<<"Corrected Truth Shower vertex: ("<<f_truth_corr_showervtxX<<", "<<f_truth_corr_showervtxY<<", "<<f_truth_corr_showervtxZ<<")"<<"\n";
	std::cout<<"Reco neutrino vertex: ("<<f_reco_nuvtxX<<", "<<f_reco_nuvtxY<<", "<<f_reco_nuvtxZ<<")"<<"\n";
	std::cout<<"Reco shower vertex: ("<<f_reco_showervtxX<<", "<<f_reco_showervtxY<<", "<<f_reco_showervtxZ<<")"<<"\n";
	std::cout<<"Reco/True shower KE: "<<f_reco_showerKE<<" / "<<f_truth_showerKE<<"\n";
	/// PF validation ends
      }

      /// BDT input variables
      if(f_BDTvars){
        auto bdthandle = e.getHandle<std::vector<nsm::NuSelectionBDT>>(fPFInputTag);
        if (! bdthandle) return;

        auto const& bdtvec = *bdthandle;
	std::cout<<"--- NuSelectionBDT ---"<<std::endl;
	if(bdtvec.size()>1) {
		std::cout<<"WARNING: >1 set of BDT input variables" << std::endl;
		return;
	} 
	if(bdtvec.size()<1) {
		//  this maybe redundant since the initialization is done in class instance.
		M_tag=-1;
		M_shower_energy=-1; // MeV
		M_shower_angle=-1; // degree, relative to beam
		M_seg_max_length=-1; // cm 
		M_seg_max_angle=-1; // degree, relative to shower
		M_acc_forward_length=-1; // cm
		M_acc_backward_length=-1; // cm
		M_flag_shwtrunck_outside=-1;
		E_tag=-1;
		E_shower_energy=-1; // MeV
		E_shower_length=-1; // cm
		E_shower_angle=-1; // degree, relative to beam
		E_muon_energy=-1; // MeV
		E_muon_length=-1; // cm
		E_flag_numuCC=-1;
		flag_others=-1;
	} 
	for(size_t i=0; i<bdtvec.size(); i++){
/*		art::Ptr<nsm::NuSelectionBDT> bdt = bdtvec.at(i);
		M_tag=bdt->GetCaseMinput().tag;
		M_shower_energy=bdt->GetCaseMinput().shower_energy; // MeV
		M_shower_angle=bdt->GetCaseMinput().shower_angle; // degree, relative to beam
		M_seg_max_length=bdt->GetCaseMinput().seg_max_length; // cm 
		M_seg_max_angle=bdt->GetCaseMinput().seg_max_angle; // degree, relative to shower
		M_acc_forward_length=bdt->GetCaseMinput().acc_forward_length; // cm
		M_acc_backward_length=bdt->GetCaseMinput().acc_backward_length; // cm
		M_flag_shwtrunck_outside=bdt->GetCaseMinput().flag_shwtrunck_outside;
		E_tag=bdt->GetCaseEinput().tag;
		E_shower_energy=bdt->GetCaseEinput().shower_energy; // MeV
		E_shower_length=bdt->GetCaseEinput().shower_length; // cm
		E_shower_angle=bdt->GetCaseEinput().shower_angle; // degree, relative to beam
		E_muon_energy=bdt->GetCaseEinput().muon_energy; // MeV
		E_muon_length=bdt->GetCaseEinput().muon_length; // cm
		E_flag_numuCC=bdt->GetCaseEinput().flag_numuCC;
		flag_others=bdt->GetFlagOthers();
*/
	}
	std::cout<<"BDT input vars check: \n"<<
	"Case M: "<<
	M_tag<<" "<<
	M_shower_energy<<" "<< // MeV
	M_shower_angle<<" "<< // degree, relative to beam
	M_seg_max_length<<" "<< // cm 
	M_seg_max_angle<<" "<< // degree, relative to shower
	M_acc_forward_length<<" "<< // cm
	M_acc_backward_length<<" "<< // cm
	M_flag_shwtrunck_outside<<"\n"
	"Case E: "<<
	E_tag<<" "<<
	E_shower_energy<<" "<< // MeV
	E_shower_length<<" "<< // cm
	E_shower_angle<<" "<< // degree, relative to beam
	E_muon_energy<<" "<< // MeV
	E_muon_length<<" "<< // cm
	E_flag_numuCC<<"\n"
	"Flag others: "<<flag_others<<"\n";
      }   
	
	// [obselete] Do selection & Fill histograms
	/*
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
	else if( f_match_energy<fcut_lowenergy ){
		if(f_truth_isEligible) h2_match_energy->Fill(f_truth_energyInside);		
		if(f_NC_truth_isEligible) h2_NC_match_energy->Fill(f_truth_energyInside);		
	}
	else if( f_lightmismatch || f_stm_LM || f_stm_lowenergy){
		if(f_truth_isEligible) h2_lightmismatch->Fill(f_truth_energyInside);		
		if(f_NC_truth_isEligible) h2_NC_lightmismatch->Fill(f_truth_energyInside);		
	}
	else if( f_match_isTgm || f_stm_TGM ){
		if(f_truth_isEligible) h2_Tgm->Fill(f_truth_energyInside);		
		if(f_NC_truth_isEligible) h2_NC_Tgm->Fill(f_truth_energyInside);		
	}
	else if( f_stm_STM || f_stm_FullDead || f_stm_clusterlength<=fcut_clusterlength ){
		if(f_truth_isEligible) h2_STM->Fill(f_truth_energyInside);
		if(f_NC_truth_isEligible) h2_NC_STM->Fill(f_truth_energyInside);
		
		h1_total->Fill(f_match_energy);
		if( fMC && f_match_completeness_energy/f_truth_energyInside>fcut_completeness ){
			if( f_truth_isCC && f_truth_nuPdg==14 && f_truth_vtxInside) h1_numuCC_FV->Fill(f_match_energy);
			if( f_truth_isCC && f_truth_nuPdg==14 && !f_truth_vtxInside) h1_numuCC_nFV->Fill(f_match_energy);		
			if( !f_truth_isCC && f_truth_nuPdg==14 && f_truth_vtxInside) h1_numuNC_FV->Fill(f_match_energy);
			if( !f_truth_isCC && f_truth_nuPdg==14 && !f_truth_vtxInside) h1_numuNC_nFV->Fill(f_match_energy);		
			if( f_truth_nuPdg==-14) h1_numubar->Fill(f_match_energy);		
			if( f_truth_nuPdg==12) h1_nue->Fill(f_match_energy);		
			if( f_truth_nuPdg==-12) h1_nuebar->Fill(f_match_energy);		
		}
		else if( fMC ){
			h1_cosmic->Fill(f_match_energy);
		}	
	
	}
	else{
		if( f_match_completeness_energy/f_truth_energyInside<fcut_completeness ){
			if(f_truth_isEligible) h2_cosmic_match->Fill(f_truth_energyInside);
			if(f_NC_truth_isEligible) h2_NC_cosmic_match->Fill(f_truth_energyInside);
		}
		else{
			if(f_truth_isEligible) h2_nuCC_STMselected->Fill(f_truth_energyInside);
			if(f_NC_truth_isEligible) h2_nuNC_STMselected->Fill(f_truth_energyInside);
		}
		
		h1_STM_total->Fill(f_match_energy);
		if( fMC && f_match_completeness_energy/f_truth_energyInside>fcut_completeness ){
			if( f_truth_isCC && f_truth_nuPdg==14 && f_truth_vtxInside) h1_STM_numuCC_FV->Fill(f_match_energy);
			if( f_truth_isCC && f_truth_nuPdg==14 && !f_truth_vtxInside) h1_STM_numuCC_nFV->Fill(f_match_energy);		
			if( !f_truth_isCC && f_truth_nuPdg==14 && f_truth_vtxInside) h1_STM_numuNC_FV->Fill(f_match_energy);
			if( !f_truth_isCC && f_truth_nuPdg==14 && !f_truth_vtxInside) h1_STM_numuNC_nFV->Fill(f_match_energy);		
			if( f_truth_nuPdg==-14) h1_STM_numubar->Fill(f_match_energy);		
			if( f_truth_nuPdg==12) h1_STM_nue->Fill(f_match_energy);		
			if( f_truth_nuPdg==-12) h1_STM_nuebar->Fill(f_match_energy);		
		}
		else if( fMC ){
			h1_STM_cosmic->Fill(f_match_energy);
		}	
	}
	/// histogram filling
	*/	
	fTreeEval->Fill();
	fPFeval->Fill();
	fBDT->Fill();
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

void WCPselection::ShowerID(int trackId)
{
  //std::cout<<"TrackID: "<<trackId<<std::endl;
  // nested loop to fill ShowerID into a vector
  // initial trackId should be looped over primary particles
  auto const& p = fParticleMap[trackId];
  // if (p.PdgCode() == 111) return; // not pi0 induced --> deferred to tagger
  const TLorentzVector& mom = p.Momentum(0);
  if (mom.E()-mom.M()<fthreshold_showerKE) return; //energy threshold
  if (p.PdgCode() == 11 || p.PdgCode() == -11) fShowerID.push_back(trackId); // key: fill shower ID
  
  // keep looping
  int Ndaughters =  p.NumberDaughters();
  if (Ndaughters == 0) return;
  for(int i=0; i<Ndaughters; i++){
    //std::cout<<"Daughter ID: "<<p.Daughter(i)<<std::endl;
    ShowerID(p.Daughter(i));
  }
}

void WCPselection::MuonID(int trackId)
{
  auto const& p = fParticleMap[trackId];
  if (p.Mother()==0 && (p.PdgCode() == 13 || p.PdgCode() == -13)) fMuonID.push_back(trackId);
}


void WCPselection::resetOutput()
{
	// live period within each event
	// maybe redundant here
	M_tag=-1;
	M_shower_energy=-1; // MeV
	M_shower_angle=-1; // degree, relative to beam
	M_seg_max_length=-1; // cm 
	M_seg_max_angle=-1; // degree, relative to shower
	M_acc_forward_length=-1; // cm
	M_acc_backward_length=-1; // cm
	M_flag_shwtrunck_outside=-1;
	E_tag=-1;
	E_shower_energy=-1; // MeV
	E_shower_length=-1; // cm
	E_shower_angle=-1; // degree, relative to beam
	E_muon_energy=-1; // MeV
	E_muon_length=-1; // cm
	E_flag_numuCC=-1;
        flag_others=-1;
	
	f_neutrino_type = -1;
	f_reco_nuvtxX = -1;
	f_reco_nuvtxY = -1;
	f_reco_nuvtxZ = -1;
	f_reco_showervtxX = -1; // primary shower [highest energy electron & not pi0 daughters] 
	f_reco_showervtxY = -1;
	f_reco_showervtxZ = -1;
	f_reco_showerKE = -1;
	f_reco_muonvtxX = -1; //
	f_reco_muonvtxY = -1;
	f_reco_muonvtxZ = -1;
	f_reco_muonMomentum[0] = -1;
	f_reco_muonMomentum[1] = -1;
	f_reco_muonMomentum[2] = -1;
	f_reco_muonMomentum[3] = -1;
	f_truth_corr_nuvtxX = -1; // truth -(SCE)-> SED -(nu time offset)-> reco [trigger offset O(10) ns ignored]
	f_truth_corr_nuvtxY = -1;
	f_truth_corr_nuvtxZ = -1;
	f_truth_corr_showervtxX = -1; // truth -(SCE)-> SED -(nu time offset)-> reco [trigger offset O(10) ns ignored]
	f_truth_corr_showervtxY = -1;
	f_truth_corr_showervtxZ = -1;
	f_truth_showerKE = -1;
	f_truth_corr_muonvtxX = -1; //
	f_truth_corr_muonvtxY = -1;
	f_truth_corr_muonvtxZ = -1;
	f_truth_muonvtxX = -1; //
	f_truth_muonvtxY = -1;
	f_truth_muonvtxZ = -1;
	f_truth_muonendX = -1; //
	f_truth_muonendY = -1;
	f_truth_muonendZ = -1;
	f_truth_muonMomentum[0] = -1;
	f_truth_muonMomentum[1] = -1;
	f_truth_muonMomentum[2] = -1;
	f_truth_muonMomentum[3] = -1;
	f_truth_nuIntType = -1;

	fPrimaryID.clear();	
	fShowerID.clear();
	fMuonID.clear();
	fPrimaryID.shrink_to_fit();
	fShowerID.shrink_to_fit();
	fMuonID.shrink_to_fit();
	fParticleMap.clear();

 	f_weight_spline = -1.0;
	f_weight_cv = -1.0;
	
	f_stm_eventtype = -1;
	f_stm_lowenergy = -1;
	f_stm_LM = -1;
	f_stm_TGM = -1;
	f_stm_STM = -1;
	f_stm_FullDead = -1;
	f_stm_clusterlength = -1;
}

void WCPselection::save_weights(art::Event const& e)
{ 
  // Use the EventWeight producer label here
  art::Handle<std::vector<evwgh::MCEventWeight> > weightsHandle;
  e.getByLabel("eventweight", weightsHandle);
  
  // Loop through these objects for each neutrino vertex in the event
  for(size_t i=0; i<weightsHandle->size(); i++){
    const evwgh::MCEventWeight& mc_weights = weightsHandle->at(i);
    // Loop over all of the weights in the MCEventWeight object
    for ( const auto& pair : mc_weights.fWeight ) {
      std::string knob_name = pair.first;
      std::vector<double> weights = pair.second;
      //std::cout<<"Knob name: "<<knob_name<<std::endl; 
      //std::cout<<"Weight size: "<<weights.size()<<std::endl; 

      if( knob_name == "TunedCentralValue_UBGenie"){
          f_weight_cv = weights.at(0);
      }
      if (knob_name == "splines_general_Spline"){
          f_weight_spline = weights.at(0);
      }
    }
  }
  
  //std::cout<<"cv weight: "<<f_weight_cv<<std::endl;
  //std::cout<<"spline weight: "<<f_weight_spline<<std::endl;

}

DEFINE_ART_MODULE(WCPselection)
