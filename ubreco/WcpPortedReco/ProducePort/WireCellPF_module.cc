////////////////////////////////////////////////////////////////////////
// Class:       WireCellPF
// Plugin Type: producer (art v3_01_02)
// File:        WireCellPF_module.cc
//
// Generated at Fri May  8 11:02:15 2020 by Hanyu Wei using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "ubobj/WcpPort/NuSelectionBDT.h"
#include "ubobj/WcpPort/NuSelectionKINE.h"

#include <memory>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

#define MAX_TRACKS 1000

namespace nsm {
  class WireCellPF;
}


class nsm::WireCellPF : public art::EDProducer {
public:
  explicit WireCellPF(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WireCellPF(WireCellPF const&) = delete;
  WireCellPF(WireCellPF&&) = delete;
  WireCellPF& operator=(WireCellPF const&) = delete;
  WireCellPF& operator=(WireCellPF&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  unsigned int NumberOfPF;
  std::string fInput; // input ROOT file for each event
  std::string fInput_prefix; // file name prefix
  std::string fInput_tree;  // particle flow
  std::string fInput_tree2; // nu selection tagger results
  std::string fInput_tree3; // BDT input variables
  std::string fInput_tree4; // KINE input variables

  bool f_ssmBDT;
  bool f_PFport;
  bool f_BDTport;
  bool f_KINEport;
};


nsm::WireCellPF::WireCellPF(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
  , NumberOfPF    (0)
  , fInput        ("")
  , fInput_prefix (p.get<std::string>("PFInput_prefix", "nue")  )
  , fInput_tree   (p.get<std::string>("PFInput_tree", "TMC")    )
  , fInput_tree2  (p.get<std::string>("PFInput_tree2", "T_match"))
  , fInput_tree3  (p.get<std::string>("PFInput_BDT", "T_tagger"))
  , fInput_tree4  (p.get<std::string>("PFInput_KINE", "T_kine"))
  , f_ssmBDT      (p.get<bool>("ssmBDT", false))
  , f_PFport	  (p.get<bool>("PFport", true))
  , f_BDTport	  (p.get<bool>("BDTport", true))
  , f_KINEport	  (p.get<bool>("KINEport", true))
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  MF_LOG_DEBUG("WireCellPF") << "Debug: WireCellPF() begins";
  if(f_PFport) produces< std::vector<simb::MCParticle> >();
  if(f_BDTport) produces< std::vector<nsm::NuSelectionBDT> >();
  if(f_KINEport) produces< std::vector<nsm::NuSelectionKINE> >();
  if (f_PFport) produces< std::vector<double> >("NewPMTInfoPePred");
  if (f_PFport) produces< std::vector<double> >("NewPMTInfoPeMeas");
  if (f_PFport) produces< std::vector<double> >("NewPMTInfoPeMeasErr");
  if (f_PFport) produces< int >("NewPMTInfoTPCClusterID");
  if (f_PFport) produces< int >("NewPMTInfoFlashID");
  if (f_PFport) produces< double >("NewPMTInfoStrength");
  if (f_PFport) produces< int >("NewPMTInfoEventType");
  if (f_PFport) produces< double >("NewPMTInfoKSDistance");
  if (f_PFport) produces< double >("NewPMTInfoChi2");
  if (f_PFport) produces< int >("NewPMTInfoNDF");
  if (f_PFport) produces< double >("NewPMTInfoClusterLength");
  if (f_PFport) produces< int >("NewPMTInfoNeutrinoType");
  if (f_PFport) produces< double >("NewPMTInfoFlashTime");
  MF_LOG_DEBUG("WireCellPF") << "Debug: WireCellPF() ends";

}

void nsm::WireCellPF::produce(art::Event& e)
{
  // Implementation of required member function here.
  bool badinput = false;
  TFile* fin = nullptr;
  std::string event_runinfo = std::to_string((int)e.run())+"_"+std::to_string((int)e.subRun())+"_"+std::to_string((int)e.event());
  fInput = fInput_prefix+"_"+event_runinfo+".root";
  mf::LogInfo("WireCellPF") <<"INPUT FILE NAME: "<< fInput <<"\n";
  std::ifstream fcheck(fInput.c_str());
  if(!fcheck.good()) {
	mf::LogInfo("WireCellPF") <<"INPUT FILE NOT FOUND!"<<std::endl;
	badinput = true;
  }
  else{
	fin = new TFile(fInput.c_str());
	if(fin->IsZombie()){
	  mf::LogError("WireCellPF") <<"INPUT FILE CANNOT OPEN: "<< fInput <<"\n";
	  badinput = true;
	}
  }

std::cout << "f_PFport: " << f_PFport << "\n";

if(f_PFport){
  auto outputPF = std::make_unique< std::vector<simb::MCParticle> >();
  if(badinput){
	e.put(std::move(outputPF));
    std::cout << "badinput, not loading from T_match tree\n";
    // Put empty products for all NewPMTInfo* products when input is bad
    e.put(std::make_unique<std::vector<double>>(), "NewPMTInfoPePred");
    e.put(std::make_unique<std::vector<double>>(), "NewPMTInfoPeMeas");
    e.put(std::make_unique<std::vector<double>>(), "NewPMTInfoPeMeasErr");
    e.put(std::make_unique<int>(-1), "NewPMTInfoTPCClusterID");
    e.put(std::make_unique<int>(-1), "NewPMTInfoFlashID");
    e.put(std::make_unique<double>(-1.0), "NewPMTInfoStrength");
    e.put(std::make_unique<int>(-1), "NewPMTInfoEventType");
    e.put(std::make_unique<double>(-1.0), "NewPMTInfoKSDistance");
    e.put(std::make_unique<double>(-1.0), "NewPMTInfoChi2");
    e.put(std::make_unique<int>(-1), "NewPMTInfoNDF");
    e.put(std::make_unique<double>(-1.0), "NewPMTInfoClusterLength");
    e.put(std::make_unique<int>(-1), "NewPMTInfoNeutrinoType");
    e.put(std::make_unique<double>(-1.0), "NewPMTInfoFlashTime");
	//return;
  }
  else{
  /// various nu taggers saved in an integer:
  // 2nd bit: cosmic (including dirt)
  // 3rd bit: numu CC
  // 4th bit: NC
  // 5th bit: long muon
  // 6th bit: nue CC

  std::cout << "preparing to read T_match tree\n";

  Int_t tpc_cluster_id;
  Int_t flash_id;
  Double_t strength;
  Double_t pe_pred[32];
  Double_t pe_meas[32];
  Double_t pe_meas_err[32];
  Int_t event_type;
  Double_t ks_dis;
  Double_t chi2;
  Int_t ndf;
  Double_t cluster_length;
  Int_t neutrino_type;
  Double_t flash_time;
  
  TTree *tree2 = (TTree*)fin->Get(fInput_tree2.c_str());
  if(tree2) {
    std::cout << "in WireCellPF, tree2 containing T_match found\n";
    tree2->SetBranchStatus("*", 0);

    tree2->SetBranchStatus("tpc_cluster_id", 1);
    tree2->SetBranchAddress("tpc_cluster_id",&tpc_cluster_id);

    tree2->SetBranchStatus("flash_id", 1);
    tree2->SetBranchAddress("flash_id",&flash_id);

    tree2->SetBranchStatus("strength", 1);
    tree2->SetBranchAddress("strength",&strength);

    tree2->SetBranchStatus("pe_pred", 1);
    tree2->SetBranchAddress("pe_pred",&pe_pred);

    tree2->SetBranchStatus("pe_meas", 1);
    tree2->SetBranchAddress("pe_meas",&pe_meas);

    tree2->SetBranchStatus("pe_meas_err", 1);
    tree2->SetBranchAddress("pe_meas_err",&pe_meas_err);

    tree2->SetBranchStatus("event_type", 1);
    tree2->SetBranchAddress("event_type",&event_type);

    tree2->SetBranchStatus("ks_dis", 1);
    tree2->SetBranchAddress("ks_dis",&ks_dis);

    tree2->SetBranchStatus("chi2", 1);
    tree2->SetBranchAddress("chi2",&chi2);

    tree2->SetBranchStatus("ndf", 1);
    tree2->SetBranchAddress("ndf",&ndf);

    tree2->SetBranchStatus("cluster_length", 1);
    tree2->SetBranchAddress("cluster_length",&cluster_length);

    tree2->SetBranchStatus("neutrino_type", 1);
    tree2->SetBranchAddress("neutrino_type",&neutrino_type);

    tree2->SetBranchStatus("flash_time", 1);
    tree2->SetBranchAddress("flash_time",&flash_time);

    for(int i=0; i<tree2->GetEntries(); i++){
      tree2->GetEntry(i);
      if(neutrino_type>1) break; // this should be the in-beam flash match
    }
    std::cout << "Should be in-beam flash match, neutrino_type: " << neutrino_type << "\n";

    std::cout << "T_match/pe_meas[32] values: ";
    for(int i=0; i<32; i++){
      std::cout << pe_meas[i] << " ";
    }
    std::cout << "\n";

    std::cout << "T_match/pe_pred[32] values: ";
    for(int i=0; i<32; i++){
      std::cout << pe_pred[i] << " ";
    }
    std::cout << "\n";

    // put pred_pe and meas_pe into the resulting file with e.put(std::move( ))
    
    // Create vectors for pe_pred and pe_meas data
    auto output_pe_pred = std::make_unique< std::vector<double> >();
    auto output_pe_meas = std::make_unique< std::vector<double> >();
    auto output_pe_meas_err = std::make_unique< std::vector<double> >();

    // Copy the arrays to vectors
    for(int i=0; i<32; i++){
      output_pe_pred->push_back(pe_pred[i]);
      output_pe_meas->push_back(pe_meas[i]);
      output_pe_meas_err->push_back(pe_meas_err[i]);
    }
    
    e.put(std::move(output_pe_pred), "NewPMTInfoPePred");
    e.put(std::move(output_pe_meas), "NewPMTInfoPeMeas");
    e.put(std::move(output_pe_meas_err), "NewPMTInfoPeMeasErr");
    e.put(std::make_unique<int>(tpc_cluster_id), "NewPMTInfoTPCClusterID");
    e.put(std::make_unique<int>(flash_id), "NewPMTInfoFlashID");
    e.put(std::make_unique<double>(strength), "NewPMTInfoStrength");
    e.put(std::make_unique<int>(event_type), "NewPMTInfoEventType");
    e.put(std::make_unique<double>(ks_dis), "NewPMTInfoKSDistance");
    e.put(std::make_unique<double>(chi2), "NewPMTInfoChi2");
    e.put(std::make_unique<int>(ndf), "NewPMTInfoNDF");
    e.put(std::make_unique<double>(cluster_length), "NewPMTInfoClusterLength");
    e.put(std::make_unique<int>(neutrino_type), "NewPMTInfoNeutrinoType");
    e.put(std::make_unique<double>(flash_time), "NewPMTInfoFlashTime");

  } else {
    std::cout << "T_match tree not found, filling with defaults\n";
    e.put(std::make_unique<std::vector<double>>(), "NewPMTInfoPePred");
    e.put(std::make_unique<std::vector<double>>(), "NewPMTInfoPeMeas");
    e.put(std::make_unique<std::vector<double>>(), "NewPMTInfoPeMeasErr");
    e.put(std::make_unique<int>(-1), "NewPMTInfoTPCClusterID");
    e.put(std::make_unique<int>(-1), "NewPMTInfoFlashID");
    e.put(std::make_unique<double>(-1.0), "NewPMTInfoStrength");
    e.put(std::make_unique<int>(-1), "NewPMTInfoEventType");
    e.put(std::make_unique<double>(-1.0), "NewPMTInfoKSDistance");
    e.put(std::make_unique<double>(-1.0), "NewPMTInfoChi2");
    e.put(std::make_unique<int>(-1), "NewPMTInfoNDF");
    e.put(std::make_unique<double>(-1.0), "NewPMTInfoClusterLength");
    e.put(std::make_unique<int>(-1), "NewPMTInfoNeutrinoType");
    e.put(std::make_unique<double>(-1.0), "NewPMTInfoFlashTime");
  }

  TTree *tree = (TTree*)fin->Get(fInput_tree.c_str());
  if(tree) {
  /// Wire-Cell Particle Flow
  int mc_Ntrack;  // number of tracks in MC
  int mc_id[MAX_TRACKS];  // track id; size == mc_Ntrack
  int mc_pdg[MAX_TRACKS];  // track particle pdg; size == mc_Ntrack
  int mc_process[MAX_TRACKS];  // track generation process code; size == mc_Ntrack
  int mc_mother[MAX_TRACKS];  // mother id of this track; size == mc_Ntrack
  int mc_included[MAX_TRACKS];  // whether this particle should be included in the energy calculation; size == mc_Ntrack
  float mc_startXYZT[MAX_TRACKS][4];  // start position of this track; size == mc_Ntrack
  float mc_endXYZT[MAX_TRACKS][4];  // end position of this track; size == mc_Ntrack
  float mc_startMomentum[MAX_TRACKS][4];  // start momentum of this track; size == mc_Ntrack
  float mc_endMomentum[MAX_TRACKS][4];  // end momentum of this track; size == mc_Ntrack
  std::vector<std::vector<int> > *mc_daughters = new std::vector<std::vector<int>>;  // daughters id of this track; vector
  //std::cout<<"Check point 0.1: "<<std::endl;

  tree->SetBranchAddress("mc_Ntrack"       , &mc_Ntrack);
  tree->SetBranchAddress("mc_id"           , &mc_id);
  tree->SetBranchAddress("mc_pdg"          , &mc_pdg);
  tree->SetBranchAddress("mc_process"      , &mc_process);
  tree->SetBranchAddress("mc_mother"       , &mc_mother);
  tree->SetBranchAddress("mc_included"     , &mc_included);
  tree->SetBranchAddress("mc_daughters"    , &mc_daughters);
  tree->SetBranchAddress("mc_startXYZT"    , &mc_startXYZT);     // unit: cm
  tree->SetBranchAddress("mc_endXYZT"      , &mc_endXYZT);       // unit: cm
  tree->SetBranchAddress("mc_startMomentum", &mc_startMomentum); // unit: GeV
  tree->SetBranchAddress("mc_endMomentum"  , &mc_endMomentum);   // unit: GeV
  //std::cout<<"Check point 0.2: "<<std::endl;
  tree->GetEntry(0); // one file one event (entry)
  //std::cout<<"Check point 0.3: "<<std::endl;

  for(int i=0; i<mc_Ntrack; i++){
	//std::cout<<"Check point 1: "<<i<<std::endl;
	// (trackID, pdg, string process, motherID, mass, status)
	// status --> neutrino_type for our reco level taggers
	// mass --> not an actual mass reconstruction since PID tells us the mass already if you believe
	// mass --> container of "mc_included" info to indicate if this particle should be included in the nu energy calculation
	simb::MCParticle p(mc_id[i], mc_pdg[i], "unknown", mc_mother[i], mc_included[i], neutrino_type);
	//std::cout<<"Check point 2: "<<i<<std::endl;
  	int Ndaughter = mc_daughters->at(i).size();
	//std::cout<<"Check point 3: "<<i<<std::endl;
	for(int j=0; j<Ndaughter; j++ ){
		p.AddDaughter(mc_daughters->at(i).at(j));
	}
	/// Currently only two points: start and end
	TLorentzVector startposition(mc_startXYZT[i][0], mc_startXYZT[i][1], mc_startXYZT[i][2], mc_startXYZT[i][3]);
	TLorentzVector endposition(mc_endXYZT[i][0], mc_endXYZT[i][1], mc_endXYZT[i][2], mc_endXYZT[i][3]);
	TLorentzVector startMomentum(mc_startMomentum[i][0], mc_startMomentum[i][1], mc_startMomentum[i][2], mc_startMomentum[i][3]);
	TLorentzVector endMomentum(mc_endMomentum[i][0], mc_endMomentum[i][1], mc_endMomentum[i][2], mc_endMomentum[i][3]);
	p.AddTrajectoryPoint(startposition, startMomentum);
	p.AddTrajectoryPoint(endposition, endMomentum);
	///
  	outputPF->push_back(std::move(p));
  }

  NumberOfPF++;

  }
  else {
    mf::LogError("WireCellPF") <<"TTree "<< fInput_tree <<" not found in file " << fInput <<"\n";
  }

  e.put(std::move(outputPF));
  }
}

float numu_cc_flag;
if(f_BDTport){
  auto outputBDTvars = std::make_unique< std::vector<nsm::NuSelectionBDT> >();
  if(badinput){
	e.put(std::move(outputBDTvars));
	//return;
  }
  else{
  nsm::NuSelectionBDT nsmbdt;

  TTree *tree3 = (TTree*)fin->Get(fInput_tree3.c_str());
  if(tree3){

  /// define variables and set branch address

  float ssm_flag_st_kdar = 0.;
  float ssm_Nsm = 0.;
  float ssm_Nsm_wivtx = 0.;

  float ssm_dq_dx_fwd_1 = 0.;
  float ssm_dq_dx_fwd_2 = 0.;
  float ssm_dq_dx_fwd_3 = 0.;
  float ssm_dq_dx_fwd_4 = 0.;
  float ssm_dq_dx_fwd_5 = 0.;
  float ssm_dq_dx_bck_1 = 0.;
  float ssm_dq_dx_bck_2 = 0.;
  float ssm_dq_dx_bck_3 = 0.;
  float ssm_dq_dx_bck_4 = 0.;
  float ssm_dq_dx_bck_5 = 0.;
  float ssm_d_dq_dx_fwd_12 = 0.;
  float ssm_d_dq_dx_fwd_23 = 0.;
  float ssm_d_dq_dx_fwd_34 = 0.;
  float ssm_d_dq_dx_fwd_45 = 0.;
  float ssm_d_dq_dx_bck_12 = 0.;
  float ssm_d_dq_dx_bck_23 = 0.;
  float ssm_d_dq_dx_bck_34 = 0.;
  float ssm_d_dq_dx_bck_45 = 0.;
  float ssm_max_dq_dx_fwd_3 = 0.;
  float ssm_max_dq_dx_fwd_5 = 0.;
  float ssm_max_dq_dx_bck_3 = 0.;
  float ssm_max_dq_dx_bck_5 = 0.;
  float ssm_max_d_dq_dx_fwd_3 = 0.;
  float ssm_max_d_dq_dx_fwd_5 = 0.;
  float ssm_max_d_dq_dx_bck_3 = 0.;
  float ssm_max_d_dq_dx_bck_5 = 0.;
  float ssm_medium_dq_dx = 0.;
  float ssm_medium_dq_dx_bp = 0.;
      //angluar info
  float ssm_angle_to_z = 0.;
  float ssm_angle_to_target = 0.;
  float ssm_angle_to_absorber = 0.;
  float ssm_angle_to_vertical = 0.;
      //directional info
  float ssm_x_dir = 0.;
  float ssm_y_dir = 0.;
  float ssm_z_dir = 0.;
      //energy info
  float ssm_kine_energy = 0.;
  float ssm_kine_energy_reduced = 0.;
      //general properties
  float ssm_vtx_activity = 0.;
  float ssm_pdg = 0.;
  float ssm_dQ_dx_cut = 0.;
  float ssm_score_mu_fwd = 0.;
  float ssm_score_p_fwd = 0.;
  float ssm_score_e_fwd = 0.;
  float ssm_score_mu_bck = 0.;
  float ssm_score_p_bck = 0.;
  float ssm_score_e_bck = 0.;
  float ssm_score_mu_fwd_bp = 0.;
  float ssm_score_p_fwd_bp = 0.;
  float ssm_score_e_fwd_bp = 0.;
      //track "straighness"
  float ssm_length = 0.;
  float ssm_direct_length = 0.;
  float ssm_length_ratio = 0.;
  float ssm_max_dev = 0.;
    //number of other particles
  float ssm_n_prim_tracks_1 = 0.;
  float ssm_n_prim_tracks_3 = 0.;
  float ssm_n_prim_tracks_5 = 0.;
  float ssm_n_prim_tracks_8 = 0.;
  float ssm_n_prim_tracks_11 = 0.;
  float ssm_n_all_tracks_1 = 0.;
  float ssm_n_all_tracks_3 = 0.;
  float ssm_n_all_tracks_5 = 0.;
  float ssm_n_all_tracks_8 = 0.;
  float ssm_n_all_tracks_11 = 0.;
  float ssm_n_daughter_tracks_1 = 0.;
  float ssm_n_daughter_tracks_3 = 0.;
  float ssm_n_daughter_tracks_5 = 0.;
  float ssm_n_daughter_tracks_8 = 0.;
  float ssm_n_daughter_tracks_11 = 0.;
  float ssm_n_daughter_all_1 = 0.;
  float ssm_n_daughter_all_3 = 0.;
  float ssm_n_daughter_all_5 = 0.;
  float ssm_n_daughter_all_8 = 0.;
  float ssm_n_daughter_all_11 = 0.;
    //properties of leading other primary track
  float ssm_prim_track1_pdg = 0.;
  float ssm_prim_track1_score_mu_fwd = 0.;
  float ssm_prim_track1_score_p_fwd = 0.;
  float ssm_prim_track1_score_e_fwd = 0.;
  float ssm_prim_track1_score_mu_bck = 0.;
  float ssm_prim_track1_score_p_bck = 0.;
  float ssm_prim_track1_score_e_bck = 0.;
  float ssm_prim_track1_length = 0.;
  float ssm_prim_track1_direct_length = 0.;
  float ssm_prim_track1_length_ratio = 0.;
  float ssm_prim_track1_max_dev = 0.;
  float ssm_prim_track1_kine_energy_range = 0.;
  float ssm_prim_track1_kine_energy_range_mu = 0.;
  float ssm_prim_track1_kine_energy_range_p = 0.;
  float ssm_prim_track1_kine_energy_range_e = 0.;
  float ssm_prim_track1_kine_energy_cal = 0.;
  float ssm_prim_track1_medium_dq_dx = 0.;
  float ssm_prim_track1_x_dir = 0.;
  float ssm_prim_track1_y_dir = 0.;
  float ssm_prim_track1_z_dir = 0.;
  float ssm_prim_track1_add_daught_track_counts_1 = 0.;
  float ssm_prim_track1_add_daught_all_counts_1 = 0.;
  float ssm_prim_track1_add_daught_track_counts_5 = 0.;
  float ssm_prim_track1_add_daught_all_counts_5 = 0.;
  float ssm_prim_track1_add_daught_track_counts_11 = 0.;
  float ssm_prim_track1_add_daught_all_counts_11 = 0.;
  //properties of sub-leading other primary track
  float ssm_prim_track2_pdg = 0.;
  float ssm_prim_track2_score_mu_fwd = 0.;
  float ssm_prim_track2_score_p_fwd = 0.;
  float ssm_prim_track2_score_e_fwd = 0.;
  float ssm_prim_track2_score_mu_bck = 0.;
  float ssm_prim_track2_score_p_bck = 0.;
  float ssm_prim_track2_score_e_bck = 0.;
  float ssm_prim_track2_length = 0.;
  float ssm_prim_track2_direct_length = 0.;
  float ssm_prim_track2_length_ratio = 0.;
  float ssm_prim_track2_max_dev = 0.;
  float ssm_prim_track2_kine_energy_range = 0.;
  float ssm_prim_track2_kine_energy_range_mu = 0.;
  float ssm_prim_track2_kine_energy_range_p = 0.;
  float ssm_prim_track2_kine_energy_range_e = 0.;
  float ssm_prim_track2_kine_energy_cal = 0.;
  float ssm_prim_track2_medium_dq_dx = 0.;
  float ssm_prim_track2_x_dir = 0.;
  float ssm_prim_track2_y_dir = 0.;
  float ssm_prim_track2_z_dir = 0.;
  float ssm_prim_track2_add_daught_track_counts_1 = 0.;
  float ssm_prim_track2_add_daught_all_counts_1 = 0.;
  float ssm_prim_track2_add_daught_track_counts_5 = 0.;
  float ssm_prim_track2_add_daught_all_counts_5 = 0.;
  float ssm_prim_track2_add_daught_track_counts_11 = 0.;
  float ssm_prim_track2_add_daught_all_counts_11 = 0.;
    //properties of leading daughter track
  float ssm_daught_track1_pdg = 0.;
  float ssm_daught_track1_score_mu_fwd = 0.;
  float ssm_daught_track1_score_p_fwd = 0.;
  float ssm_daught_track1_score_e_fwd = 0.;
  float ssm_daught_track1_score_mu_bck = 0.;
  float ssm_daught_track1_score_p_bck = 0.;
  float ssm_daught_track1_score_e_bck = 0.;
  float ssm_daught_track1_length = 0.;
  float ssm_daught_track1_direct_length = 0.;
  float ssm_daught_track1_length_ratio = 0.;
  float ssm_daught_track1_max_dev = 0.;
  float ssm_daught_track1_kine_energy_range = 0.;
  float ssm_daught_track1_kine_energy_range_mu = 0.;
  float ssm_daught_track1_kine_energy_range_p = 0.;
  float ssm_daught_track1_kine_energy_range_e = 0.;
  float ssm_daught_track1_kine_energy_cal = 0.;
  float ssm_daught_track1_medium_dq_dx = 0.;
  float ssm_daught_track1_x_dir = 0.;
  float ssm_daught_track1_y_dir = 0.;
  float ssm_daught_track1_z_dir = 0.;
  float ssm_daught_track1_add_daught_track_counts_1 = 0.;
  float ssm_daught_track1_add_daught_all_counts_1 = 0.;
  float ssm_daught_track1_add_daught_track_counts_5 = 0.;
  float ssm_daught_track1_add_daught_all_counts_5 = 0.;
  float ssm_daught_track1_add_daught_track_counts_11 = 0.;
  float ssm_daught_track1_add_daught_all_counts_11 = 0.;
    //properties of sub-leading daughter track
  float ssm_daught_track2_pdg = 0.;
  float ssm_daught_track2_score_mu_fwd = 0.;
  float ssm_daught_track2_score_p_fwd = 0.;
  float ssm_daught_track2_score_e_fwd = 0.;
  float ssm_daught_track2_score_mu_bck = 0.;
  float ssm_daught_track2_score_p_bck = 0.;
  float ssm_daught_track2_score_e_bck = 0.;
  float ssm_daught_track2_length = 0.;
  float ssm_daught_track2_direct_length = 0.;
  float ssm_daught_track2_length_ratio = 0.;
  float ssm_daught_track2_max_dev = 0.;
  float ssm_daught_track2_kine_energy_range = 0.;
  float ssm_daught_track2_kine_energy_range_mu = 0.;
  float ssm_daught_track2_kine_energy_range_p = 0.;
  float ssm_daught_track2_kine_energy_range_e = 0.;
  float ssm_daught_track2_kine_energy_cal = 0.;
  float ssm_daught_track2_medium_dq_dx = 0.;
  float ssm_daught_track2_x_dir = 0.;
  float ssm_daught_track2_y_dir = 0.;
  float ssm_daught_track2_z_dir = 0.;
  float ssm_daught_track2_add_daught_track_counts_1 = 0.;
  float ssm_daught_track2_add_daught_all_counts_1 = 0.;
  float ssm_daught_track2_add_daught_track_counts_5 = 0.;
  float ssm_daught_track2_add_daught_all_counts_5 = 0.;
  float ssm_daught_track2_add_daught_track_counts_11 = 0.;
  float ssm_daught_track2_add_daught_all_counts_11 = 0.;
  //properties of leading other primary shower
  float ssm_prim_shw1_pdg = 0.;
  float ssm_prim_shw1_score_mu_fwd = 0.;
  float ssm_prim_shw1_score_p_fwd = 0.;
  float ssm_prim_shw1_score_e_fwd = 0.;
  float ssm_prim_shw1_score_mu_bck = 0.;
  float ssm_prim_shw1_score_p_bck = 0.;
  float ssm_prim_shw1_score_e_bck = 0.;
  float ssm_prim_shw1_length = 0.;
  float ssm_prim_shw1_direct_length = 0.;
  float ssm_prim_shw1_length_ratio = 0.;
  float ssm_prim_shw1_max_dev = 0.;
  float ssm_prim_shw1_kine_energy_range = 0.;
  float ssm_prim_shw1_kine_energy_range_mu = 0.;
  float ssm_prim_shw1_kine_energy_range_p = 0.;
  float ssm_prim_shw1_kine_energy_range_e = 0.;
  float ssm_prim_shw1_kine_energy_cal = 0.;
  float ssm_prim_shw1_kine_energy_best = 0.;
  float ssm_prim_shw1_medium_dq_dx = 0.;
  float ssm_prim_shw1_x_dir = 0.;
  float ssm_prim_shw1_y_dir = 0.;
  float ssm_prim_shw1_z_dir = 0.;
  float ssm_prim_shw1_add_daught_track_counts_1 = 0.;
  float ssm_prim_shw1_add_daught_all_counts_1 = 0.;
  float ssm_prim_shw1_add_daught_track_counts_5 = 0.;
  float ssm_prim_shw1_add_daught_all_counts_5 = 0.;
  float ssm_prim_shw1_add_daught_track_counts_11 = 0.;
  float ssm_prim_shw1_add_daught_all_counts_11 = 0.;
    //properties of sub-leading other primary shower
  float ssm_prim_shw2_pdg = 0.;
  float ssm_prim_shw2_score_mu_fwd = 0.;
  float ssm_prim_shw2_score_p_fwd = 0.;
  float ssm_prim_shw2_score_e_fwd = 0.;
  float ssm_prim_shw2_score_mu_bck = 0.;
  float ssm_prim_shw2_score_p_bck = 0.;
  float ssm_prim_shw2_score_e_bck = 0.;
  float ssm_prim_shw2_length = 0.;
  float ssm_prim_shw2_direct_length = 0.;
  float ssm_prim_shw2_length_ratio = 0.;
  float ssm_prim_shw2_max_dev = 0.;
  float ssm_prim_shw2_kine_energy_range = 0.;
  float ssm_prim_shw2_kine_energy_range_mu = 0.;
  float ssm_prim_shw2_kine_energy_range_p = 0.;
  float ssm_prim_shw2_kine_energy_range_e = 0.;
  float ssm_prim_shw2_kine_energy_cal = 0.;
  float ssm_prim_shw2_kine_energy_best = 0.;
  float ssm_prim_shw2_medium_dq_dx = 0.;
  float ssm_prim_shw2_x_dir = 0.;
  float ssm_prim_shw2_y_dir = 0.;
  float ssm_prim_shw2_z_dir = 0.;
  float ssm_prim_shw2_add_daught_track_counts_1 = 0.;
  float ssm_prim_shw2_add_daught_all_counts_1 = 0.;
  float ssm_prim_shw2_add_daught_track_counts_5 = 0.;
  float ssm_prim_shw2_add_daught_all_counts_5 = 0.;
  float ssm_prim_shw2_add_daught_track_counts_11 = 0.;
  float ssm_prim_shw2_add_daught_all_counts_11 = 0.;
  //properties of leading daughter shower
  float ssm_daught_shw1_pdg = 0.;
  float ssm_daught_shw1_score_mu_fwd = 0.;
  float ssm_daught_shw1_score_p_fwd = 0.;
  float ssm_daught_shw1_score_e_fwd = 0.;
  float ssm_daught_shw1_score_mu_bck = 0.;
  float ssm_daught_shw1_score_p_bck = 0.;
  float ssm_daught_shw1_score_e_bck = 0.;
  float ssm_daught_shw1_length = 0.;
  float ssm_daught_shw1_direct_length = 0.;
  float ssm_daught_shw1_length_ratio = 0.;
  float ssm_daught_shw1_max_dev = 0.;
  float ssm_daught_shw1_kine_energy_range = 0.;
  float ssm_daught_shw1_kine_energy_range_mu = 0.;
  float ssm_daught_shw1_kine_energy_range_p = 0.;
  float ssm_daught_shw1_kine_energy_range_e = 0.;
  float ssm_daught_shw1_kine_energy_cal = 0.;
  float ssm_daught_shw1_kine_energy_best = 0.;
  float ssm_daught_shw1_medium_dq_dx = 0.;
  float ssm_daught_shw1_x_dir = 0.;
  float ssm_daught_shw1_y_dir = 0.;
  float ssm_daught_shw1_z_dir = 0.;
  float ssm_daught_shw1_add_daught_track_counts_1 = 0.;
  float ssm_daught_shw1_add_daught_all_counts_1 = 0.;
  float ssm_daught_shw1_add_daught_track_counts_5 = 0.;
  float ssm_daught_shw1_add_daught_all_counts_5 = 0.;
  float ssm_daught_shw1_add_daught_track_counts_11 = 0.;
  float ssm_daught_shw1_add_daught_all_counts_11 = 0.;
    //properties of sub-leading daughter shower
  float ssm_daught_shw2_pdg = 0.;
  float ssm_daught_shw2_score_mu_fwd = 0.;
  float ssm_daught_shw2_score_p_fwd = 0.;
  float ssm_daught_shw2_score_e_fwd = 0.;
  float ssm_daught_shw2_score_mu_bck = 0.;
  float ssm_daught_shw2_score_p_bck = 0.;
  float ssm_daught_shw2_score_e_bck = 0.;
  float ssm_daught_shw2_length = 0.;
  float ssm_daught_shw2_direct_length = 0.;
  float ssm_daught_shw2_length_ratio = 0.;
  float ssm_daught_shw2_max_dev = 0.;
  float ssm_daught_shw2_kine_energy_range = 0.;
  float ssm_daught_shw2_kine_energy_range_mu = 0.;
  float ssm_daught_shw2_kine_energy_range_p = 0.;
  float ssm_daught_shw2_kine_energy_range_e = 0.;
  float ssm_daught_shw2_kine_energy_cal = 0.;
  float ssm_daught_shw2_kine_energy_best = 0.;
  float ssm_daught_shw2_medium_dq_dx = 0.;
  float ssm_daught_shw2_x_dir = 0.;
  float ssm_daught_shw2_y_dir = 0.;
  float ssm_daught_shw2_z_dir = 0.;
  float ssm_daught_shw2_add_daught_track_counts_1 = 0.;
  float ssm_daught_shw2_add_daught_all_counts_1 = 0.;
  float ssm_daught_shw2_add_daught_track_counts_5 = 0.;
  float ssm_daught_shw2_add_daught_all_counts_5 = 0.;
  float ssm_daught_shw2_add_daught_track_counts_11 = 0.;
  float ssm_daught_shw2_add_daught_all_counts_11 = 0.;
    //event level properties
  float ssm_nu_angle_z = 0.;
  float ssm_nu_angle_target = 0.;
  float ssm_nu_angle_absorber = 0.;
  float ssm_nu_angle_vertical = 0.;
  float ssm_con_nu_angle_z = 0.;
  float ssm_con_nu_angle_target = 0.;
  float ssm_con_nu_angle_absorber = 0.;
  float ssm_con_nu_angle_vertical = 0.;
  float ssm_prim_nu_angle_z = 0.;
  float ssm_prim_nu_angle_target = 0.;
  float ssm_prim_nu_angle_absorber = 0.;
  float ssm_prim_nu_angle_vertical = 0.;
  float ssm_track_angle_z = 0.;
  float ssm_track_angle_target = 0.;
  float ssm_track_angle_absorber = 0.;
  float ssm_track_angle_vertical = 0.;
  float ssm_vtxX = 0.;
  float ssm_vtxY = 0.;
  float ssm_vtxZ = 0.;
    //off vertex stuff
  float ssm_offvtx_length = 0.;
  float ssm_offvtx_energy = 0.;
  float ssm_n_offvtx_tracks_1 = 0.;
  float ssm_n_offvtx_tracks_3 = 0.;
  float ssm_n_offvtx_tracks_5 = 0.;
  float ssm_n_offvtx_tracks_8 = 0.;
  float ssm_n_offvtx_tracks_11 = 0.;
  float ssm_n_offvtx_showers_1 = 0.;
  float ssm_n_offvtx_showers_3 = 0.;
  float ssm_n_offvtx_showers_5 = 0.;
  float ssm_n_offvtx_showers_8 = 0.;
  float ssm_n_offvtx_showers_11 = 0.;
    //properties of leading off vertex track
  float ssm_offvtx_track1_pdg = 0.;
  float ssm_offvtx_track1_score_mu_fwd = 0.;
  float ssm_offvtx_track1_score_p_fwd = 0.;
  float ssm_offvtx_track1_score_e_fwd = 0.;
  float ssm_offvtx_track1_score_mu_bck = 0.;
  float ssm_offvtx_track1_score_p_bck = 0.;
  float ssm_offvtx_track1_score_e_bck = 0.;
  float ssm_offvtx_track1_length = 0.;
  float ssm_offvtx_track1_direct_length = 0.;
  float ssm_offvtx_track1_max_dev = 0.;
  float ssm_offvtx_track1_kine_energy_range = 0.;
  float ssm_offvtx_track1_kine_energy_range_mu = 0.;
  float ssm_offvtx_track1_kine_energy_range_p = 0.;
  float ssm_offvtx_track1_kine_energy_range_e = 0.;
  float ssm_offvtx_track1_kine_energy_cal = 0.;
  float ssm_offvtx_track1_medium_dq_dx = 0.;
  float ssm_offvtx_track1_x_dir = 0.;
  float ssm_offvtx_track1_y_dir = 0.;
  float ssm_offvtx_track1_z_dir = 0.;
  float ssm_offvtx_track1_dist_mainvtx = 0.;
    //properties of leading off vertex shower
  float ssm_offvtx_shw1_pdg_offvtx = 0.;
  float ssm_offvtx_shw1_score_mu_fwd = 0.;
  float ssm_offvtx_shw1_score_p_fwd = 0.;
  float ssm_offvtx_shw1_score_e_fwd = 0.;
  float ssm_offvtx_shw1_score_mu_bck = 0.;
  float ssm_offvtx_shw1_score_p_bck = 0.;
  float ssm_offvtx_shw1_score_e_bck = 0.;
  float ssm_offvtx_shw1_length = 0.;
  float ssm_offvtx_shw1_direct_length = 0.;
  float ssm_offvtx_shw1_max_dev = 0.;
  float ssm_offvtx_shw1_kine_energy_best = 0.;
  float ssm_offvtx_shw1_kine_energy_range = 0.;
  float ssm_offvtx_shw1_kine_energy_range_mu = 0.;
  float ssm_offvtx_shw1_kine_energy_range_p = 0.;
  float ssm_offvtx_shw1_kine_energy_range_e = 0.;
  float ssm_offvtx_shw1_kine_energy_cal = 0.;
  float ssm_offvtx_shw1_medium_dq_dx = 0.;
  float ssm_offvtx_shw1_x_dir = 0.;
  float ssm_offvtx_shw1_y_dir = 0.;
  float ssm_offvtx_shw1_z_dir = 0.;
  float ssm_offvtx_shw1_dist_mainvtx = 0.;
    // Spacepoints
  int ssmsp_Ntrack = 0;
  std::vector<int> *ssmsp_Nsp= new std::vector<int>;
  int ssmsp_Nsp_tot = 0;
  std::vector<int> *ssmsp_pdg= new std::vector<int>;
  std::vector<int> *ssmsp_id= new std::vector<int>;
  std::vector<int> *ssmsp_mother= new std::vector<int>;
  std::vector<float> *ssmsp_x= new std::vector<float>;
  std::vector<float> *ssmsp_y= new std::vector<float>;
  std::vector<float> *ssmsp_z= new std::vector<float>;
  std::vector<float> *ssmsp_dx= new std::vector<float>;
  std::vector<float> *ssmsp_dQ= new std::vector<float>;
  std::vector<float> *ssmsp_KE= new std::vector<float>;
  std::vector<float> *ssmsp_containing_shower_id= new std::vector<float>;
  std::vector<float> *ssmsp_containing_shower_ke= new std::vector<float>;
  std::vector<float> *ssmsp_containing_shower_flag= new std::vector<float>;
    //Kine vars
  float ssm_kine_reco_Enu = 0.; // ssm_kinetic energy  + additional energy ...
  float ssm_kine_reco_add_energy = 0.;  // mass, binding energy ...
  std::vector<float> *ssm_kine_energy_particle = new std::vector<float>;  // energy of each particle
  std::vector<int> *ssm_kine_energy_info = new std::vector<int>; // what kind of energy reconstruction?
  std::vector<int> *ssm_kine_particle_type = new std::vector<int>;
  std::vector<int> *ssm_kine_energy_included = new std::vector<int>; // included in the neutrino energy calculation?
  float ssm_kine_pio_mass = 0.; // mass
  int ssm_kine_pio_flag = 0; // 0 not filled, 1, with vertex: CCpio, 2 without vertex: NCpi0
  float ssm_kine_pio_vtx_dis = 0.;
  float ssm_kine_pio_energy_1 = 0.;
  float ssm_kine_pio_theta_1 = 0.;
  float ssm_kine_pio_phi_1 = 0.;
  float ssm_kine_pio_dis_1 = 0.;
  float ssm_kine_pio_energy_2 = 0.;
  float ssm_kine_pio_theta_2 = 0.;
  float ssm_kine_pio_phi_2 = 0.;
  float ssm_kine_pio_dis_2 = 0.;
  float ssm_kine_pio_angle = 0.;

  if(f_ssmBDT){

    tree3->SetBranchAddress("ssm_flag_st_kdar",&ssm_flag_st_kdar);
    tree3->SetBranchAddress("ssm_Nsm",&ssm_Nsm);
    tree3->SetBranchAddress("ssm_Nsm_wivtx",&ssm_Nsm_wivtx);

    //only filled if there is one ssm
    //properties of the ssm
    //dq/dx info
      tree3->SetBranchAddress("ssm_dq_dx_fwd_1", &ssm_dq_dx_fwd_1);
      tree3->SetBranchAddress("ssm_dq_dx_fwd_2", &ssm_dq_dx_fwd_2);
      tree3->SetBranchAddress("ssm_dq_dx_fwd_3", &ssm_dq_dx_fwd_3);
      tree3->SetBranchAddress("ssm_dq_dx_fwd_4", &ssm_dq_dx_fwd_4);
      tree3->SetBranchAddress("ssm_dq_dx_fwd_5", &ssm_dq_dx_fwd_5);
      tree3->SetBranchAddress("ssm_dq_dx_bck_1", &ssm_dq_dx_bck_1);
      tree3->SetBranchAddress("ssm_dq_dx_bck_2", &ssm_dq_dx_bck_2);
      tree3->SetBranchAddress("ssm_dq_dx_bck_3", &ssm_dq_dx_bck_3);
      tree3->SetBranchAddress("ssm_dq_dx_bck_4", &ssm_dq_dx_bck_4);
      tree3->SetBranchAddress("ssm_dq_dx_bck_5", &ssm_dq_dx_bck_5);
      tree3->SetBranchAddress("ssm_d_dq_dx_fwd_12", &ssm_d_dq_dx_fwd_12);
      tree3->SetBranchAddress("ssm_d_dq_dx_fwd_23", &ssm_d_dq_dx_fwd_23);
      tree3->SetBranchAddress("ssm_d_dq_dx_fwd_34", &ssm_d_dq_dx_fwd_34);
      tree3->SetBranchAddress("ssm_d_dq_dx_fwd_45", &ssm_d_dq_dx_fwd_45);
      tree3->SetBranchAddress("ssm_d_dq_dx_bck_12", &ssm_d_dq_dx_bck_12);
      tree3->SetBranchAddress("ssm_d_dq_dx_bck_23", &ssm_d_dq_dx_bck_23);
      tree3->SetBranchAddress("ssm_d_dq_dx_bck_34", &ssm_d_dq_dx_bck_34);
      tree3->SetBranchAddress("ssm_d_dq_dx_bck_45", &ssm_d_dq_dx_bck_45);
      tree3->SetBranchAddress("ssm_max_dq_dx_fwd_3", &ssm_max_dq_dx_fwd_3);
      tree3->SetBranchAddress("ssm_max_dq_dx_fwd_5", &ssm_max_dq_dx_fwd_5);
      tree3->SetBranchAddress("ssm_max_dq_dx_bck_3", &ssm_max_dq_dx_bck_3);
      tree3->SetBranchAddress("ssm_max_dq_dx_bck_5", &ssm_max_dq_dx_bck_5);
      tree3->SetBranchAddress("ssm_max_d_dq_dx_fwd_3", &ssm_max_d_dq_dx_fwd_3);
      tree3->SetBranchAddress("ssm_max_d_dq_dx_fwd_5", &ssm_max_d_dq_dx_fwd_5);
      tree3->SetBranchAddress("ssm_max_d_dq_dx_bck_3", &ssm_max_d_dq_dx_bck_3);
      tree3->SetBranchAddress("ssm_max_d_dq_dx_bck_5", &ssm_max_d_dq_dx_bck_5);
      tree3->SetBranchAddress("ssm_medium_dq_dx", &ssm_medium_dq_dx);
      tree3->SetBranchAddress("ssm_medium_dq_dx_bp", &ssm_medium_dq_dx_bp);
    //angluar info
      tree3->SetBranchAddress("ssm_angle_to_z", &ssm_angle_to_z);
      tree3->SetBranchAddress("ssm_angle_to_target", &ssm_angle_to_target);
      tree3->SetBranchAddress("ssm_angle_to_absorber", &ssm_angle_to_absorber);
      tree3->SetBranchAddress("ssm_angle_to_vertical", &ssm_angle_to_vertical);
    //directional info
      tree3->SetBranchAddress("ssm_x_dir", &ssm_x_dir);
      tree3->SetBranchAddress("ssm_y_dir", &ssm_y_dir);
      tree3->SetBranchAddress("ssm_z_dir", &ssm_z_dir);
    //energy info
      tree3->SetBranchAddress("ssm_kine_energy", &ssm_kine_energy);
      tree3->SetBranchAddress("ssm_kine_energy_reduced", &ssm_kine_energy_reduced);
    //general properties
      tree3->SetBranchAddress("ssm_vtx_activity", &ssm_vtx_activity);
      tree3->SetBranchAddress("ssm_pdg", &ssm_pdg);
      tree3->SetBranchAddress("ssm_dQ_dx_cut", &ssm_dQ_dx_cut);
      tree3->SetBranchAddress("ssm_score_mu_fwd", &ssm_score_mu_fwd);
      tree3->SetBranchAddress("ssm_score_p_fwd", &ssm_score_p_fwd);
      tree3->SetBranchAddress("ssm_score_e_fwd", &ssm_score_e_fwd);
      tree3->SetBranchAddress("ssm_score_mu_bck", &ssm_score_mu_bck);
      tree3->SetBranchAddress("ssm_score_p_bck", &ssm_score_p_bck);
      tree3->SetBranchAddress("ssm_score_e_bck", &ssm_score_e_bck);
      tree3->SetBranchAddress("ssm_score_mu_fwd_bp", &ssm_score_mu_fwd_bp);
      tree3->SetBranchAddress("ssm_score_p_fwd_bp", &ssm_score_p_fwd_bp);
      tree3->SetBranchAddress("ssm_score_e_fwd_bp", &ssm_score_e_fwd_bp);
    //track "straighness"
      tree3->SetBranchAddress("ssm_length", &ssm_length);
      tree3->SetBranchAddress("ssm_direct_length", &ssm_direct_length);
      tree3->SetBranchAddress("ssm_length_ratio", &ssm_length_ratio);
      tree3->SetBranchAddress("ssm_max_dev", &ssm_max_dev);
    //number of other particles
      tree3->SetBranchAddress("ssm_n_prim_tracks_1", &ssm_n_prim_tracks_1);
      tree3->SetBranchAddress("ssm_n_prim_tracks_3", &ssm_n_prim_tracks_3);
      tree3->SetBranchAddress("ssm_n_prim_tracks_5", &ssm_n_prim_tracks_5);
      tree3->SetBranchAddress("ssm_n_prim_tracks_8", &ssm_n_prim_tracks_8);
      tree3->SetBranchAddress("ssm_n_prim_tracks_11", &ssm_n_prim_tracks_11);
      tree3->SetBranchAddress("ssm_n_all_tracks_1", &ssm_n_all_tracks_1);
      tree3->SetBranchAddress("ssm_n_all_tracks_3", &ssm_n_all_tracks_3);
      tree3->SetBranchAddress("ssm_n_all_tracks_5", &ssm_n_all_tracks_5);
      tree3->SetBranchAddress("ssm_n_all_tracks_8", &ssm_n_all_tracks_8);
      tree3->SetBranchAddress("ssm_n_all_tracks_11", &ssm_n_all_tracks_11);
      tree3->SetBranchAddress("ssm_n_daughter_tracks_1", &ssm_n_daughter_tracks_1);
      tree3->SetBranchAddress("ssm_n_daughter_tracks_3", &ssm_n_daughter_tracks_3);
      tree3->SetBranchAddress("ssm_n_daughter_tracks_5", &ssm_n_daughter_tracks_5);
      tree3->SetBranchAddress("ssm_n_daughter_tracks_8", &ssm_n_daughter_tracks_8);
      tree3->SetBranchAddress("ssm_n_daughter_tracks_11", &ssm_n_daughter_tracks_11);
      tree3->SetBranchAddress("ssm_n_daughter_all_1", &ssm_n_daughter_all_1);
      tree3->SetBranchAddress("ssm_n_daughter_all_3", &ssm_n_daughter_all_3);
      tree3->SetBranchAddress("ssm_n_daughter_all_5", &ssm_n_daughter_all_5);
      tree3->SetBranchAddress("ssm_n_daughter_all_8", &ssm_n_daughter_all_8);
      tree3->SetBranchAddress("ssm_n_daughter_all_11", &ssm_n_daughter_all_11);
    //properties of leading other primary track
      tree3->SetBranchAddress("ssm_prim_track1_pdg", &ssm_prim_track1_pdg);
      tree3->SetBranchAddress("ssm_prim_track1_score_mu_fwd", &ssm_prim_track1_score_mu_fwd);
      tree3->SetBranchAddress("ssm_prim_track1_score_p_fwd", &ssm_prim_track1_score_p_fwd);
      tree3->SetBranchAddress("ssm_prim_track1_score_e_fwd", &ssm_prim_track1_score_e_fwd);
      tree3->SetBranchAddress("ssm_prim_track1_score_mu_bck", &ssm_prim_track1_score_mu_bck);
      tree3->SetBranchAddress("ssm_prim_track1_score_p_bck", &ssm_prim_track1_score_p_bck);
      tree3->SetBranchAddress("ssm_prim_track1_score_e_bck", &ssm_prim_track1_score_e_bck);
      tree3->SetBranchAddress("ssm_prim_track1_length", &ssm_prim_track1_length);
      tree3->SetBranchAddress("ssm_prim_track1_direct_length", &ssm_prim_track1_direct_length);
      tree3->SetBranchAddress("ssm_prim_track1_length_ratio", &ssm_prim_track1_length_ratio);
      tree3->SetBranchAddress("ssm_prim_track1_max_dev", &ssm_prim_track1_max_dev);
      tree3->SetBranchAddress("ssm_prim_track1_kine_energy_range", &ssm_prim_track1_kine_energy_range);
      tree3->SetBranchAddress("ssm_prim_track1_kine_energy_range_mu", &ssm_prim_track1_kine_energy_range_mu);
      tree3->SetBranchAddress("ssm_prim_track1_kine_energy_range_p", &ssm_prim_track1_kine_energy_range_p);
      tree3->SetBranchAddress("ssm_prim_track1_kine_energy_range_e", &ssm_prim_track1_kine_energy_range_e);
      tree3->SetBranchAddress("ssm_prim_track1_kine_energy_cal", &ssm_prim_track1_kine_energy_cal);
      tree3->SetBranchAddress("ssm_prim_track1_medium_dq_dx", &ssm_prim_track1_medium_dq_dx);
      tree3->SetBranchAddress("ssm_prim_track1_x_dir", &ssm_prim_track1_x_dir);
      tree3->SetBranchAddress("ssm_prim_track1_y_dir", &ssm_prim_track1_y_dir);
      tree3->SetBranchAddress("ssm_prim_track1_z_dir", &ssm_prim_track1_z_dir);
      tree3->SetBranchAddress("ssm_prim_track1_add_daught_track_counts_1", &ssm_prim_track1_add_daught_track_counts_1);
      tree3->SetBranchAddress("ssm_prim_track1_add_daught_all_counts_1", &ssm_prim_track1_add_daught_all_counts_1);
      tree3->SetBranchAddress("ssm_prim_track1_add_daught_track_counts_5", &ssm_prim_track1_add_daught_track_counts_5);
      tree3->SetBranchAddress("ssm_prim_track1_add_daught_all_counts_5", &ssm_prim_track1_add_daught_all_counts_5);
      tree3->SetBranchAddress("ssm_prim_track1_add_daught_track_counts_11", &ssm_prim_track1_add_daught_track_counts_11);
      tree3->SetBranchAddress("ssm_prim_track1_add_daught_all_counts_11", &ssm_prim_track1_add_daught_all_counts_11);
    //properties of sub-leading other primary track
      tree3->SetBranchAddress("ssm_prim_track2_pdg", &ssm_prim_track2_pdg);
      tree3->SetBranchAddress("ssm_prim_track2_score_mu_fwd", &ssm_prim_track2_score_mu_fwd);
      tree3->SetBranchAddress("ssm_prim_track2_score_p_fwd", &ssm_prim_track2_score_p_fwd);
      tree3->SetBranchAddress("ssm_prim_track2_score_e_fwd", &ssm_prim_track2_score_e_fwd);
      tree3->SetBranchAddress("ssm_prim_track2_score_mu_bck", &ssm_prim_track2_score_mu_bck);
      tree3->SetBranchAddress("ssm_prim_track2_score_p_bck", &ssm_prim_track2_score_p_bck);
      tree3->SetBranchAddress("ssm_prim_track2_score_e_bck", &ssm_prim_track2_score_e_bck);
      tree3->SetBranchAddress("ssm_prim_track2_length", &ssm_prim_track2_length);
      tree3->SetBranchAddress("ssm_prim_track2_direct_length", &ssm_prim_track2_direct_length);
      tree3->SetBranchAddress("ssm_prim_track2_length_ratio", &ssm_prim_track2_length_ratio);
      tree3->SetBranchAddress("ssm_prim_track2_max_dev", &ssm_prim_track2_max_dev);
      tree3->SetBranchAddress("ssm_prim_track2_kine_energy_range", &ssm_prim_track2_kine_energy_range);
      tree3->SetBranchAddress("ssm_prim_track2_kine_energy_range_mu", &ssm_prim_track2_kine_energy_range_mu);
      tree3->SetBranchAddress("ssm_prim_track2_kine_energy_range_p", &ssm_prim_track2_kine_energy_range_p);
      tree3->SetBranchAddress("ssm_prim_track2_kine_energy_range_e", &ssm_prim_track2_kine_energy_range_e);
      tree3->SetBranchAddress("ssm_prim_track2_kine_energy_cal", &ssm_prim_track2_kine_energy_cal);
      tree3->SetBranchAddress("ssm_prim_track2_medium_dq_dx", &ssm_prim_track2_medium_dq_dx);
      tree3->SetBranchAddress("ssm_prim_track2_x_dir", &ssm_prim_track2_x_dir);
      tree3->SetBranchAddress("ssm_prim_track2_y_dir", &ssm_prim_track2_y_dir);
      tree3->SetBranchAddress("ssm_prim_track2_z_dir", &ssm_prim_track2_z_dir);
      tree3->SetBranchAddress("ssm_prim_track2_add_daught_track_counts_1", &ssm_prim_track2_add_daught_track_counts_1);
      tree3->SetBranchAddress("ssm_prim_track2_add_daught_all_counts_1", &ssm_prim_track2_add_daught_all_counts_1);
      tree3->SetBranchAddress("ssm_prim_track2_add_daught_track_counts_5", &ssm_prim_track2_add_daught_track_counts_5);
      tree3->SetBranchAddress("ssm_prim_track2_add_daught_all_counts_5", &ssm_prim_track2_add_daught_all_counts_5);
      tree3->SetBranchAddress("ssm_prim_track2_add_daught_track_counts_11", &ssm_prim_track2_add_daught_track_counts_11);
      tree3->SetBranchAddress("ssm_prim_track2_add_daught_all_counts_11", &ssm_prim_track2_add_daught_all_counts_11);
    //properties of leading daughter track
      tree3->SetBranchAddress("ssm_daught_track1_pdg", &ssm_daught_track1_pdg);
      tree3->SetBranchAddress("ssm_daught_track1_score_mu_fwd", &ssm_daught_track1_score_mu_fwd);
      tree3->SetBranchAddress("ssm_daught_track1_score_p_fwd", &ssm_daught_track1_score_p_fwd);
      tree3->SetBranchAddress("ssm_daught_track1_score_e_fwd", &ssm_daught_track1_score_e_fwd);
      tree3->SetBranchAddress("ssm_daught_track1_score_mu_bck", &ssm_daught_track1_score_mu_bck);
      tree3->SetBranchAddress("ssm_daught_track1_score_p_bck", &ssm_daught_track1_score_p_bck);
      tree3->SetBranchAddress("ssm_daught_track1_score_e_bck", &ssm_daught_track1_score_e_bck);
      tree3->SetBranchAddress("ssm_daught_track1_length", &ssm_daught_track1_length);
      tree3->SetBranchAddress("ssm_daught_track1_direct_length", &ssm_daught_track1_direct_length);
      tree3->SetBranchAddress("ssm_daught_track1_length_ratio", &ssm_daught_track1_length_ratio);
      tree3->SetBranchAddress("ssm_daught_track1_max_dev", &ssm_daught_track1_max_dev);
      tree3->SetBranchAddress("ssm_daught_track1_kine_energy_range", &ssm_daught_track1_kine_energy_range);
      tree3->SetBranchAddress("ssm_daught_track1_kine_energy_range_mu", &ssm_daught_track1_kine_energy_range_mu);
      tree3->SetBranchAddress("ssm_daught_track1_kine_energy_range_p", &ssm_daught_track1_kine_energy_range_p);
      tree3->SetBranchAddress("ssm_daught_track1_kine_energy_range_e", &ssm_daught_track1_kine_energy_range_e);
      tree3->SetBranchAddress("ssm_daught_track1_kine_energy_cal", &ssm_daught_track1_kine_energy_cal);
      tree3->SetBranchAddress("ssm_daught_track1_medium_dq_dx", &ssm_daught_track1_medium_dq_dx);
      tree3->SetBranchAddress("ssm_daught_track1_x_dir", &ssm_daught_track1_x_dir);
      tree3->SetBranchAddress("ssm_daught_track1_y_dir", &ssm_daught_track1_y_dir);
      tree3->SetBranchAddress("ssm_daught_track1_z_dir", &ssm_daught_track1_z_dir);
      tree3->SetBranchAddress("ssm_daught_track1_add_daught_track_counts_1", &ssm_daught_track1_add_daught_track_counts_1);
      tree3->SetBranchAddress("ssm_daught_track1_add_daught_all_counts_1", &ssm_daught_track1_add_daught_all_counts_1);
      tree3->SetBranchAddress("ssm_daught_track1_add_daught_track_counts_5", &ssm_daught_track1_add_daught_track_counts_5);
      tree3->SetBranchAddress("ssm_daught_track1_add_daught_all_counts_5", &ssm_daught_track1_add_daught_all_counts_5);
      tree3->SetBranchAddress("ssm_daught_track1_add_daught_track_counts_11", &ssm_daught_track1_add_daught_track_counts_11);
      tree3->SetBranchAddress("ssm_daught_track1_add_daught_all_counts_11", &ssm_daught_track1_add_daught_all_counts_11);
    //properties of sub-leading daughter track
      tree3->SetBranchAddress("ssm_daught_track2_pdg", &ssm_daught_track2_pdg);
      tree3->SetBranchAddress("ssm_daught_track2_score_mu_fwd", &ssm_daught_track2_score_mu_fwd);
      tree3->SetBranchAddress("ssm_daught_track2_score_p_fwd", &ssm_daught_track2_score_p_fwd);
      tree3->SetBranchAddress("ssm_daught_track2_score_e_fwd", &ssm_daught_track2_score_e_fwd);
      tree3->SetBranchAddress("ssm_daught_track2_score_mu_bck", &ssm_daught_track2_score_mu_bck);
      tree3->SetBranchAddress("ssm_daught_track2_score_p_bck", &ssm_daught_track2_score_p_bck);
      tree3->SetBranchAddress("ssm_daught_track2_score_e_bck", &ssm_daught_track2_score_e_bck);
      tree3->SetBranchAddress("ssm_daught_track2_length", &ssm_daught_track2_length);
      tree3->SetBranchAddress("ssm_daught_track2_direct_length", &ssm_daught_track2_direct_length);
      tree3->SetBranchAddress("ssm_daught_track2_length_ratio", &ssm_daught_track2_length_ratio);
      tree3->SetBranchAddress("ssm_daught_track2_max_dev", &ssm_daught_track2_max_dev);
      tree3->SetBranchAddress("ssm_daught_track2_kine_energy_range", &ssm_daught_track2_kine_energy_range);
      tree3->SetBranchAddress("ssm_daught_track2_kine_energy_range_mu", &ssm_daught_track2_kine_energy_range_mu);
      tree3->SetBranchAddress("ssm_daught_track2_kine_energy_range_p", &ssm_daught_track2_kine_energy_range_p);
      tree3->SetBranchAddress("ssm_daught_track2_kine_energy_range_e", &ssm_daught_track2_kine_energy_range_e);
      tree3->SetBranchAddress("ssm_daught_track2_kine_energy_cal", &ssm_daught_track2_kine_energy_cal);
      tree3->SetBranchAddress("ssm_daught_track2_medium_dq_dx", &ssm_daught_track2_medium_dq_dx);
      tree3->SetBranchAddress("ssm_daught_track2_x_dir", &ssm_daught_track2_x_dir);
      tree3->SetBranchAddress("ssm_daught_track2_y_dir", &ssm_daught_track2_y_dir);
      tree3->SetBranchAddress("ssm_daught_track2_z_dir", &ssm_daught_track2_z_dir);
      tree3->SetBranchAddress("ssm_daught_track2_add_daught_track_counts_1", &ssm_daught_track2_add_daught_track_counts_1);
      tree3->SetBranchAddress("ssm_daught_track2_add_daught_all_counts_1", &ssm_daught_track2_add_daught_all_counts_1);
      tree3->SetBranchAddress("ssm_daught_track2_add_daught_track_counts_5", &ssm_daught_track2_add_daught_track_counts_5);
      tree3->SetBranchAddress("ssm_daught_track2_add_daught_all_counts_5", &ssm_daught_track2_add_daught_all_counts_5);
      tree3->SetBranchAddress("ssm_daught_track2_add_daught_track_counts_11", &ssm_daught_track2_add_daught_track_counts_11);
      tree3->SetBranchAddress("ssm_daught_track2_add_daught_all_counts_11", &ssm_daught_track2_add_daught_all_counts_11);
    //properties of leading other primary shower
      tree3->SetBranchAddress("ssm_prim_shw1_pdg", &ssm_prim_shw1_pdg);
      tree3->SetBranchAddress("ssm_prim_shw1_score_mu_fwd", &ssm_prim_shw1_score_mu_fwd);
      tree3->SetBranchAddress("ssm_prim_shw1_score_p_fwd", &ssm_prim_shw1_score_p_fwd);
      tree3->SetBranchAddress("ssm_prim_shw1_score_e_fwd", &ssm_prim_shw1_score_e_fwd);
      tree3->SetBranchAddress("ssm_prim_shw1_score_mu_bck", &ssm_prim_shw1_score_mu_bck);
      tree3->SetBranchAddress("ssm_prim_shw1_score_p_bck", &ssm_prim_shw1_score_p_bck);
      tree3->SetBranchAddress("ssm_prim_shw1_score_e_bck", &ssm_prim_shw1_score_e_bck);
      tree3->SetBranchAddress("ssm_prim_shw1_length", &ssm_prim_shw1_length);
      tree3->SetBranchAddress("ssm_prim_shw1_direct_length", &ssm_prim_shw1_direct_length);
      tree3->SetBranchAddress("ssm_prim_shw1_length_ratio", &ssm_prim_shw1_length_ratio);
      tree3->SetBranchAddress("ssm_prim_shw1_max_dev", &ssm_prim_shw1_max_dev);
      tree3->SetBranchAddress("ssm_prim_shw1_kine_energy_range", &ssm_prim_shw1_kine_energy_range);
      tree3->SetBranchAddress("ssm_prim_shw1_kine_energy_range_mu", &ssm_prim_shw1_kine_energy_range_mu);
      tree3->SetBranchAddress("ssm_prim_shw1_kine_energy_range_p", &ssm_prim_shw1_kine_energy_range_p);
      tree3->SetBranchAddress("ssm_prim_shw1_kine_energy_range_e", &ssm_prim_shw1_kine_energy_range_e);
      tree3->SetBranchAddress("ssm_prim_shw1_kine_energy_cal", &ssm_prim_shw1_kine_energy_cal);
      tree3->SetBranchAddress("ssm_prim_shw1_kine_energy_best", &ssm_prim_shw1_kine_energy_best);
      tree3->SetBranchAddress("ssm_prim_shw1_medium_dq_dx", &ssm_prim_shw1_medium_dq_dx);
      tree3->SetBranchAddress("ssm_prim_shw1_x_dir", &ssm_prim_shw1_x_dir);
      tree3->SetBranchAddress("ssm_prim_shw1_y_dir", &ssm_prim_shw1_y_dir);
      tree3->SetBranchAddress("ssm_prim_shw1_z_dir", &ssm_prim_shw1_z_dir);
      tree3->SetBranchAddress("ssm_prim_shw1_add_daught_track_counts_1", &ssm_prim_shw1_add_daught_track_counts_1);
      tree3->SetBranchAddress("ssm_prim_shw1_add_daught_all_counts_1", &ssm_prim_shw1_add_daught_all_counts_1);
      tree3->SetBranchAddress("ssm_prim_shw1_add_daught_track_counts_5", &ssm_prim_shw1_add_daught_track_counts_5);
      tree3->SetBranchAddress("ssm_prim_shw1_add_daught_all_counts_5", &ssm_prim_shw1_add_daught_all_counts_5);
      tree3->SetBranchAddress("ssm_prim_shw1_add_daught_track_counts_11", &ssm_prim_shw1_add_daught_track_counts_11);
      tree3->SetBranchAddress("ssm_prim_shw1_add_daught_all_counts_11", &ssm_prim_shw1_add_daught_all_counts_11);
    //properties of sub-leading other primary shower
      tree3->SetBranchAddress("ssm_prim_shw2_pdg", &ssm_prim_shw2_pdg);
      tree3->SetBranchAddress("ssm_prim_shw2_score_mu_fwd", &ssm_prim_shw2_score_mu_fwd);
      tree3->SetBranchAddress("ssm_prim_shw2_score_p_fwd", &ssm_prim_shw2_score_p_fwd);
      tree3->SetBranchAddress("ssm_prim_shw2_score_e_fwd", &ssm_prim_shw2_score_e_fwd);
      tree3->SetBranchAddress("ssm_prim_shw2_score_mu_bck", &ssm_prim_shw2_score_mu_bck);
      tree3->SetBranchAddress("ssm_prim_shw2_score_p_bck", &ssm_prim_shw2_score_p_bck);
      tree3->SetBranchAddress("ssm_prim_shw2_score_e_bck", &ssm_prim_shw2_score_e_bck);
      tree3->SetBranchAddress("ssm_prim_shw2_length", &ssm_prim_shw2_length);
      tree3->SetBranchAddress("ssm_prim_shw2_direct_length", &ssm_prim_shw2_direct_length);
      tree3->SetBranchAddress("ssm_prim_shw2_length_ratio", &ssm_prim_shw2_length_ratio);
      tree3->SetBranchAddress("ssm_prim_shw2_max_dev", &ssm_prim_shw2_max_dev);
      tree3->SetBranchAddress("ssm_prim_shw2_kine_energy_range", &ssm_prim_shw2_kine_energy_range);
      tree3->SetBranchAddress("ssm_prim_shw2_kine_energy_range_mu", &ssm_prim_shw2_kine_energy_range_mu);
      tree3->SetBranchAddress("ssm_prim_shw2_kine_energy_range_p", &ssm_prim_shw2_kine_energy_range_p);
      tree3->SetBranchAddress("ssm_prim_shw2_kine_energy_range_e", &ssm_prim_shw2_kine_energy_range_e);
      tree3->SetBranchAddress("ssm_prim_shw2_kine_energy_cal", &ssm_prim_shw2_kine_energy_cal);
      tree3->SetBranchAddress("ssm_prim_shw2_kine_energy_best", &ssm_prim_shw2_kine_energy_best);
      tree3->SetBranchAddress("ssm_prim_shw2_medium_dq_dx", &ssm_prim_shw2_medium_dq_dx);
      tree3->SetBranchAddress("ssm_prim_shw2_x_dir", &ssm_prim_shw2_x_dir);
      tree3->SetBranchAddress("ssm_prim_shw2_y_dir", &ssm_prim_shw2_y_dir);
      tree3->SetBranchAddress("ssm_prim_shw2_z_dir", &ssm_prim_shw2_z_dir);
      tree3->SetBranchAddress("ssm_prim_shw2_add_daught_track_counts_1", &ssm_prim_shw2_add_daught_track_counts_1);
      tree3->SetBranchAddress("ssm_prim_shw2_add_daught_all_counts_1", &ssm_prim_shw2_add_daught_all_counts_1);
      tree3->SetBranchAddress("ssm_prim_shw2_add_daught_track_counts_5", &ssm_prim_shw2_add_daught_track_counts_5);
      tree3->SetBranchAddress("ssm_prim_shw2_add_daught_all_counts_5", &ssm_prim_shw2_add_daught_all_counts_5);
      tree3->SetBranchAddress("ssm_prim_shw2_add_daught_track_counts_11", &ssm_prim_shw2_add_daught_track_counts_11);
      tree3->SetBranchAddress("ssm_prim_shw2_add_daught_all_counts_11", &ssm_prim_shw2_add_daught_all_counts_11);
    //properties of leading daughter shower
      tree3->SetBranchAddress("ssm_daught_shw1_pdg", &ssm_daught_shw1_pdg);
      tree3->SetBranchAddress("ssm_daught_shw1_score_mu_fwd", &ssm_daught_shw1_score_mu_fwd);
      tree3->SetBranchAddress("ssm_daught_shw1_score_p_fwd", &ssm_daught_shw1_score_p_fwd);
      tree3->SetBranchAddress("ssm_daught_shw1_score_e_fwd", &ssm_daught_shw1_score_e_fwd);
      tree3->SetBranchAddress("ssm_daught_shw1_score_mu_bck", &ssm_daught_shw1_score_mu_bck);
      tree3->SetBranchAddress("ssm_daught_shw1_score_p_bck", &ssm_daught_shw1_score_p_bck);
      tree3->SetBranchAddress("ssm_daught_shw1_score_e_bck", &ssm_daught_shw1_score_e_bck);
      tree3->SetBranchAddress("ssm_daught_shw1_length", &ssm_daught_shw1_length);
      tree3->SetBranchAddress("ssm_daught_shw1_direct_length", &ssm_daught_shw1_direct_length);
      tree3->SetBranchAddress("ssm_daught_shw1_length_ratio", &ssm_daught_shw1_length_ratio);
      tree3->SetBranchAddress("ssm_daught_shw1_max_dev", &ssm_daught_shw1_max_dev);
      tree3->SetBranchAddress("ssm_daught_shw1_kine_energy_range", &ssm_daught_shw1_kine_energy_range);
      tree3->SetBranchAddress("ssm_daught_shw1_kine_energy_range_mu", &ssm_daught_shw1_kine_energy_range_mu);
      tree3->SetBranchAddress("ssm_daught_shw1_kine_energy_range_p", &ssm_daught_shw1_kine_energy_range_p);
      tree3->SetBranchAddress("ssm_daught_shw1_kine_energy_range_e", &ssm_daught_shw1_kine_energy_range_e);
      tree3->SetBranchAddress("ssm_daught_shw1_kine_energy_cal", &ssm_daught_shw1_kine_energy_cal);
      tree3->SetBranchAddress("ssm_daught_shw1_kine_energy_best", &ssm_daught_shw1_kine_energy_best);
      tree3->SetBranchAddress("ssm_daught_shw1_medium_dq_dx", &ssm_daught_shw1_medium_dq_dx);
      tree3->SetBranchAddress("ssm_daught_shw1_x_dir", &ssm_daught_shw1_x_dir);
      tree3->SetBranchAddress("ssm_daught_shw1_y_dir", &ssm_daught_shw1_y_dir);
      tree3->SetBranchAddress("ssm_daught_shw1_z_dir", &ssm_daught_shw1_z_dir);
      tree3->SetBranchAddress("ssm_daught_shw1_add_daught_track_counts_1", &ssm_daught_shw1_add_daught_track_counts_1);
      tree3->SetBranchAddress("ssm_daught_shw1_add_daught_all_counts_1", &ssm_daught_shw1_add_daught_all_counts_1);
      tree3->SetBranchAddress("ssm_daught_shw1_add_daught_track_counts_5", &ssm_daught_shw1_add_daught_track_counts_5);
      tree3->SetBranchAddress("ssm_daught_shw1_add_daught_all_counts_5", &ssm_daught_shw1_add_daught_all_counts_5);
      tree3->SetBranchAddress("ssm_daught_shw1_add_daught_track_counts_11", &ssm_daught_shw1_add_daught_track_counts_11);
      tree3->SetBranchAddress("ssm_daught_shw1_add_daught_all_counts_11", &ssm_daught_shw1_add_daught_all_counts_11);
    //properties of sub-leading daughter shower
      tree3->SetBranchAddress("ssm_daught_shw2_pdg", &ssm_daught_shw2_pdg);
      tree3->SetBranchAddress("ssm_daught_shw2_score_mu_fwd", &ssm_daught_shw2_score_mu_fwd);
      tree3->SetBranchAddress("ssm_daught_shw2_score_p_fwd", &ssm_daught_shw2_score_p_fwd);
      tree3->SetBranchAddress("ssm_daught_shw2_score_e_fwd", &ssm_daught_shw2_score_e_fwd);
      tree3->SetBranchAddress("ssm_daught_shw2_score_mu_bck", &ssm_daught_shw2_score_mu_bck);
      tree3->SetBranchAddress("ssm_daught_shw2_score_p_bck", &ssm_daught_shw2_score_p_bck);
      tree3->SetBranchAddress("ssm_daught_shw2_score_e_bck", &ssm_daught_shw2_score_e_bck);
      tree3->SetBranchAddress("ssm_daught_shw2_length", &ssm_daught_shw2_length);
      tree3->SetBranchAddress("ssm_daught_shw2_direct_length", &ssm_daught_shw2_direct_length);
      tree3->SetBranchAddress("ssm_daught_shw2_length_ratio", &ssm_daught_shw2_length_ratio);
      tree3->SetBranchAddress("ssm_daught_shw2_max_dev", &ssm_daught_shw2_max_dev);
      tree3->SetBranchAddress("ssm_daught_shw2_kine_energy_range", &ssm_daught_shw2_kine_energy_range);
      tree3->SetBranchAddress("ssm_daught_shw2_kine_energy_range_mu", &ssm_daught_shw2_kine_energy_range_mu);
      tree3->SetBranchAddress("ssm_daught_shw2_kine_energy_range_p", &ssm_daught_shw2_kine_energy_range_p);
      tree3->SetBranchAddress("ssm_daught_shw2_kine_energy_range_e", &ssm_daught_shw2_kine_energy_range_e);
      tree3->SetBranchAddress("ssm_daught_shw2_kine_energy_cal", &ssm_daught_shw2_kine_energy_cal);
      tree3->SetBranchAddress("ssm_daught_shw2_kine_energy_best", &ssm_daught_shw2_kine_energy_best);
      tree3->SetBranchAddress("ssm_daught_shw2_medium_dq_dx", &ssm_daught_shw2_medium_dq_dx);
      tree3->SetBranchAddress("ssm_daught_shw2_x_dir", &ssm_daught_shw2_x_dir);
      tree3->SetBranchAddress("ssm_daught_shw2_y_dir", &ssm_daught_shw2_y_dir);
      tree3->SetBranchAddress("ssm_daught_shw2_z_dir", &ssm_daught_shw2_z_dir);
      tree3->SetBranchAddress("ssm_daught_shw2_add_daught_track_counts_1", &ssm_daught_shw2_add_daught_track_counts_1);
      tree3->SetBranchAddress("ssm_daught_shw2_add_daught_all_counts_1", &ssm_daught_shw2_add_daught_all_counts_1);
      tree3->SetBranchAddress("ssm_daught_shw2_add_daught_track_counts_5", &ssm_daught_shw2_add_daught_track_counts_5);
      tree3->SetBranchAddress("ssm_daught_shw2_add_daught_all_counts_5", &ssm_daught_shw2_add_daught_all_counts_5);
      tree3->SetBranchAddress("ssm_daught_shw2_add_daught_track_counts_11", &ssm_daught_shw2_add_daught_track_counts_11);
      tree3->SetBranchAddress("ssm_daught_shw2_add_daught_all_counts_11", &ssm_daught_shw2_add_daught_all_counts_11);
    //event level properties
      tree3->SetBranchAddress("ssm_nu_angle_z", &ssm_nu_angle_z);
      tree3->SetBranchAddress("ssm_nu_angle_target", &ssm_nu_angle_target);
      tree3->SetBranchAddress("ssm_nu_angle_absorber", &ssm_nu_angle_absorber);
      tree3->SetBranchAddress("ssm_nu_angle_vertical", &ssm_nu_angle_vertical);
      tree3->SetBranchAddress("ssm_con_nu_angle_z", &ssm_con_nu_angle_z);
      tree3->SetBranchAddress("ssm_con_nu_angle_target", &ssm_con_nu_angle_target);
      tree3->SetBranchAddress("ssm_con_nu_angle_absorber", &ssm_con_nu_angle_absorber);
      tree3->SetBranchAddress("ssm_con_nu_angle_vertical", &ssm_con_nu_angle_vertical);
      tree3->SetBranchAddress("ssm_prim_nu_angle_z", &ssm_prim_nu_angle_z);
      tree3->SetBranchAddress("ssm_prim_nu_angle_target", &ssm_prim_nu_angle_target);
      tree3->SetBranchAddress("ssm_prim_nu_angle_absorber", &ssm_prim_nu_angle_absorber);
      tree3->SetBranchAddress("ssm_prim_nu_angle_vertical", &ssm_prim_nu_angle_vertical);
      tree3->SetBranchAddress("ssm_track_angle_z", &ssm_track_angle_z);
      tree3->SetBranchAddress("ssm_track_angle_target", &ssm_track_angle_target);
      tree3->SetBranchAddress("ssm_track_angle_absorber", &ssm_track_angle_absorber);
      tree3->SetBranchAddress("ssm_track_angle_vertical", &ssm_track_angle_vertical);
      tree3->SetBranchAddress("ssm_vtxX", &ssm_vtxX);
      tree3->SetBranchAddress("ssm_vtxY", &ssm_vtxY);
      tree3->SetBranchAddress("ssm_vtxZ", &ssm_vtxZ);
    //off vertex stuff
      tree3->SetBranchAddress("ssm_offvtx_length",&ssm_offvtx_length);
      tree3->SetBranchAddress("ssm_offvtx_energy",&ssm_offvtx_energy);
      tree3->SetBranchAddress("ssm_n_offvtx_tracks_1",&ssm_n_offvtx_tracks_1);
      tree3->SetBranchAddress("ssm_n_offvtx_tracks_3",&ssm_n_offvtx_tracks_3);
      tree3->SetBranchAddress("ssm_n_offvtx_tracks_5",&ssm_n_offvtx_tracks_5);
      tree3->SetBranchAddress("ssm_n_offvtx_tracks_8",&ssm_n_offvtx_tracks_8);
      tree3->SetBranchAddress("ssm_n_offvtx_tracks_11",&ssm_n_offvtx_tracks_11);
      tree3->SetBranchAddress("ssm_n_offvtx_showers_1",&ssm_n_offvtx_showers_1);
      tree3->SetBranchAddress("ssm_n_offvtx_showers_3",&ssm_n_offvtx_showers_3);
      tree3->SetBranchAddress("ssm_n_offvtx_showers_5",&ssm_n_offvtx_showers_5);
      tree3->SetBranchAddress("ssm_n_offvtx_showers_8",&ssm_n_offvtx_showers_8);
      tree3->SetBranchAddress("ssm_n_offvtx_showers_11",&ssm_n_offvtx_showers_11);
    //properties of leading off vertex track
      tree3->SetBranchAddress("ssm_offvtx_track1_pdg",&ssm_offvtx_track1_pdg);
      tree3->SetBranchAddress("ssm_offvtx_track1_score_mu_fwd",&ssm_offvtx_track1_score_mu_fwd);
      tree3->SetBranchAddress("ssm_offvtx_track1_score_p_fwd",&ssm_offvtx_track1_score_p_fwd);
      tree3->SetBranchAddress("ssm_offvtx_track1_score_e_fwd",&ssm_offvtx_track1_score_e_fwd);
      tree3->SetBranchAddress("ssm_offvtx_track1_score_mu_bck",&ssm_offvtx_track1_score_mu_bck);
      tree3->SetBranchAddress("ssm_offvtx_track1_score_p_bck",&ssm_offvtx_track1_score_p_bck);
      tree3->SetBranchAddress("ssm_offvtx_track1_score_e_bck",&ssm_offvtx_track1_score_e_bck);
      tree3->SetBranchAddress("ssm_offvtx_track1_length",&ssm_offvtx_track1_length);
      tree3->SetBranchAddress("ssm_offvtx_track1_direct_length",&ssm_offvtx_track1_direct_length);
      tree3->SetBranchAddress("ssm_offvtx_track1_max_dev",&ssm_offvtx_track1_max_dev);
      tree3->SetBranchAddress("ssm_offvtx_track1_kine_energy_range",&ssm_offvtx_track1_kine_energy_range);
      tree3->SetBranchAddress("ssm_offvtx_track1_kine_energy_range_mu",&ssm_offvtx_track1_kine_energy_range_mu);
      tree3->SetBranchAddress("ssm_offvtx_track1_kine_energy_range_p",&ssm_offvtx_track1_kine_energy_range_p);
      tree3->SetBranchAddress("ssm_offvtx_track1_kine_energy_range_e",&ssm_offvtx_track1_kine_energy_range_e);
      tree3->SetBranchAddress("ssm_offvtx_track1_kine_energy_cal",&ssm_offvtx_track1_kine_energy_cal);
      tree3->SetBranchAddress("ssm_offvtx_track1_medium_dq_dx",&ssm_offvtx_track1_medium_dq_dx);
      tree3->SetBranchAddress("ssm_offvtx_track1_x_dir",&ssm_offvtx_track1_x_dir);
      tree3->SetBranchAddress("ssm_offvtx_track1_y_dir",&ssm_offvtx_track1_y_dir);
      tree3->SetBranchAddress("ssm_offvtx_track1_z_dir",&ssm_offvtx_track1_z_dir);
      tree3->SetBranchAddress("ssm_offvtx_track1_dist_mainvtx",&ssm_offvtx_track1_dist_mainvtx);
    //properties of leading off vertex shower
      tree3->SetBranchAddress("ssm_offvtx_shw1_pdg_offvtx",&ssm_offvtx_shw1_pdg_offvtx);
      tree3->SetBranchAddress("ssm_offvtx_shw1_score_mu_fwd",&ssm_offvtx_shw1_score_mu_fwd);
      tree3->SetBranchAddress("ssm_offvtx_shw1_score_p_fwd",&ssm_offvtx_shw1_score_p_fwd);
      tree3->SetBranchAddress("ssm_offvtx_shw1_score_e_fwd",&ssm_offvtx_shw1_score_e_fwd);
      tree3->SetBranchAddress("ssm_offvtx_shw1_score_mu_bck",&ssm_offvtx_shw1_score_mu_bck);
      tree3->SetBranchAddress("ssm_offvtx_shw1_score_p_bck",&ssm_offvtx_shw1_score_p_bck);
      tree3->SetBranchAddress("ssm_offvtx_shw1_score_e_bck",&ssm_offvtx_shw1_score_e_bck);
      tree3->SetBranchAddress("ssm_offvtx_shw1_length",&ssm_offvtx_shw1_length);
      tree3->SetBranchAddress("ssm_offvtx_shw1_direct_length",&ssm_offvtx_shw1_direct_length);
      tree3->SetBranchAddress("ssm_offvtx_shw1_max_dev",&ssm_offvtx_shw1_max_dev);
      tree3->SetBranchAddress("ssm_offvtx_shw1_kine_energy_best",&ssm_offvtx_shw1_kine_energy_best);
      tree3->SetBranchAddress("ssm_offvtx_shw1_kine_energy_range",&ssm_offvtx_shw1_kine_energy_range);
      tree3->SetBranchAddress("ssm_offvtx_shw1_kine_energy_range_mu",&ssm_offvtx_shw1_kine_energy_range_mu);
      tree3->SetBranchAddress("ssm_offvtx_shw1_kine_energy_range_p",&ssm_offvtx_shw1_kine_energy_range_p);
      tree3->SetBranchAddress("ssm_offvtx_shw1_kine_energy_range_e",&ssm_offvtx_shw1_kine_energy_range_e);
      tree3->SetBranchAddress("ssm_offvtx_shw1_kine_energy_cal",&ssm_offvtx_shw1_kine_energy_cal);
      tree3->SetBranchAddress("ssm_offvtx_shw1_medium_dq_dx",&ssm_offvtx_shw1_medium_dq_dx);
      tree3->SetBranchAddress("ssm_offvtx_shw1_x_dir",&ssm_offvtx_shw1_x_dir);
      tree3->SetBranchAddress("ssm_offvtx_shw1_y_dir",&ssm_offvtx_shw1_y_dir);
      tree3->SetBranchAddress("ssm_offvtx_shw1_z_dir",&ssm_offvtx_shw1_z_dir);
      tree3->SetBranchAddress("ssm_offvtx_shw1_dist_mainvtx",&ssm_offvtx_shw1_dist_mainvtx);
    // Spacepoints
      tree3->SetBranchAddress("ssmsp_Ntrack", &ssmsp_Ntrack);
      tree3->SetBranchAddress("ssmsp_Nsp", &ssmsp_Nsp);
      tree3->SetBranchAddress("ssmsp_Nsp_tot", &ssmsp_Nsp_tot);
      tree3->SetBranchAddress("ssmsp_pdg", &ssmsp_pdg);
      tree3->SetBranchAddress("ssmsp_id", &ssmsp_id);
      tree3->SetBranchAddress("ssmsp_mother", &ssmsp_mother);
      tree3->SetBranchAddress("ssmsp_x", &ssmsp_x);
      tree3->SetBranchAddress("ssmsp_y", &ssmsp_y);
      tree3->SetBranchAddress("ssmsp_z", &ssmsp_z);
      tree3->SetBranchAddress("ssmsp_dx", &ssmsp_dx);
      tree3->SetBranchAddress("ssmsp_dQ", &ssmsp_dQ);
      tree3->SetBranchAddress("ssmsp_KE", &ssmsp_KE);
      tree3->SetBranchAddress("ssmsp_containing_shower_id", &ssmsp_containing_shower_id);
      tree3->SetBranchAddress("ssmsp_containing_shower_ke", &ssmsp_containing_shower_ke);
      tree3->SetBranchAddress("ssmsp_containing_shower_flag", &ssmsp_containing_shower_flag);
    // Kine vars
      if(f_KINEport){
        tree3->SetBranchAddress("ssm_kine_reco_Enu", &ssm_kine_reco_Enu); // ssm_kinetic energy  + additional energy ...
        tree3->SetBranchAddress("ssm_kine_reco_add_energy", &ssm_kine_reco_add_energy);  // mass, binding energy ...
        tree3->SetBranchAddress("ssm_kine_energy_particle", &ssm_kine_energy_particle);  // energy of each particle
        tree3->SetBranchAddress("ssm_kine_energy_info", &ssm_kine_energy_info); // what kind of energy reconstruction?
        tree3->SetBranchAddress("ssm_kine_particle_type", &ssm_kine_particle_type);
        tree3->SetBranchAddress("ssm_kine_energy_included", &ssm_kine_energy_included); // included in the neutrino energy calculation?
        tree3->SetBranchAddress("ssm_kine_pio_mass", &ssm_kine_pio_mass); // mass
        tree3->SetBranchAddress("ssm_kine_pio_flag", &ssm_kine_pio_flag); // 0 not filled, 1, with vertex: CCpio, 2 without vertex: NCpi0
        tree3->SetBranchAddress("ssm_kine_pio_vtx_dis", &ssm_kine_pio_vtx_dis);
        tree3->SetBranchAddress("ssm_kine_pio_energy_1", &ssm_kine_pio_energy_1);
        tree3->SetBranchAddress("ssm_kine_pio_theta_1", &ssm_kine_pio_theta_1);
        tree3->SetBranchAddress("ssm_kine_pio_phi_1", &ssm_kine_pio_phi_1);
        tree3->SetBranchAddress("ssm_kine_pio_dis_1", &ssm_kine_pio_dis_1);
        tree3->SetBranchAddress("ssm_kine_pio_energy_2", &ssm_kine_pio_energy_2);
        tree3->SetBranchAddress("ssm_kine_pio_theta_2", &ssm_kine_pio_theta_2);
        tree3->SetBranchAddress("ssm_kine_pio_phi_2", &ssm_kine_pio_phi_2);
        tree3->SetBranchAddress("ssm_kine_pio_dis_2", &ssm_kine_pio_dis_2);
        tree3->SetBranchAddress("ssm_kine_pio_angle", &ssm_kine_pio_angle);
      }
    }

  // single photon
  float shw_sp_num_mip_tracks;
  float shw_sp_num_muons;
  float shw_sp_num_pions;
  float shw_sp_num_protons;
  float shw_sp_proton_length_1;
  float shw_sp_proton_dqdx_1;
  float shw_sp_proton_energy_1;
  float shw_sp_proton_length_2;
  float shw_sp_proton_dqdx_2;
  float shw_sp_proton_energy_2;
  float shw_sp_n_good_showers;
  float shw_sp_n_20mev_showers;
  float shw_sp_n_br1_showers;
  float shw_sp_n_br2_showers;
  float shw_sp_n_br3_showers;
  float shw_sp_n_br4_showers;
  float shw_sp_n_20br1_showers;
  std::vector<int> *shw_sp_20mev_showers = new std::vector<int>;
  std::vector<int> *shw_sp_br1_showers = new std::vector<int>;
  std::vector<int> *shw_sp_br2_showers = new std::vector<int>;
  std::vector<int> *shw_sp_br3_showers = new std::vector<int>;
  std::vector<int> *shw_sp_br4_showers = new std::vector<int>;
  float shw_sp_shw_vtx_dis;
  float shw_sp_max_shw_dis;
  tree3->SetBranchAddress("shw_sp_num_mip_tracks",&shw_sp_num_mip_tracks);
  tree3->SetBranchAddress("shw_sp_num_muons",&shw_sp_num_muons);
  tree3->SetBranchAddress("shw_sp_num_pions",&shw_sp_num_pions);
  tree3->SetBranchAddress("shw_sp_num_protons",&shw_sp_num_protons);
  tree3->SetBranchAddress("shw_sp_proton_length_1",&shw_sp_proton_length_1);
  tree3->SetBranchAddress("shw_sp_proton_dqdx_1",&shw_sp_proton_dqdx_1);
  tree3->SetBranchAddress("shw_sp_proton_energy_1",&shw_sp_proton_energy_1);
  tree3->SetBranchAddress("shw_sp_proton_length_2",&shw_sp_proton_length_2);
  tree3->SetBranchAddress("shw_sp_proton_dqdx_2",&shw_sp_proton_dqdx_2);
  tree3->SetBranchAddress("shw_sp_proton_energy_2",&shw_sp_proton_energy_2);
  tree3->SetBranchAddress("shw_sp_n_good_showers",&shw_sp_n_good_showers);
  tree3->SetBranchAddress("shw_sp_n_20mev_showers",&shw_sp_n_20mev_showers);
  tree3->SetBranchAddress("shw_sp_n_br1_showers",&shw_sp_n_br1_showers);
  tree3->SetBranchAddress("shw_sp_n_br2_showers",&shw_sp_n_br2_showers);
  tree3->SetBranchAddress("shw_sp_n_br3_showers",&shw_sp_n_br3_showers);
  tree3->SetBranchAddress("shw_sp_n_br4_showers",&shw_sp_n_br4_showers);
  tree3->SetBranchAddress("shw_sp_n_20br1_showers",&shw_sp_n_20br1_showers);
  tree3->SetBranchAddress("shw_sp_20mev_showers",&shw_sp_20mev_showers);
  tree3->SetBranchAddress("shw_sp_br1_showers",&shw_sp_br1_showers);
  tree3->SetBranchAddress("shw_sp_br2_showers",&shw_sp_br2_showers);
  tree3->SetBranchAddress("shw_sp_br3_showers",&shw_sp_br3_showers);
  tree3->SetBranchAddress("shw_sp_br4_showers",&shw_sp_br4_showers);
  tree3->SetBranchAddress("shw_sp_shw_vtx_dis",&shw_sp_shw_vtx_dis);
  tree3->SetBranchAddress("shw_sp_max_shw_dis",&shw_sp_max_shw_dis);

  float shw_sp_filled;
  float shw_sp_flag;
  float shw_sp_energy;
  float shw_sp_vec_dQ_dx_0;
  float shw_sp_vec_dQ_dx_1;
  float shw_sp_max_dQ_dx_sample;
  float shw_sp_n_below_threshold;
  float shw_sp_n_below_zero;
  float shw_sp_n_lowest;
  float shw_sp_n_highest;
  float shw_sp_lowest_dQ_dx;
  float shw_sp_highest_dQ_dx;
  float shw_sp_medium_dQ_dx;
  float shw_sp_stem_length;
  float shw_sp_length_main;
  float shw_sp_length_total;
  float shw_sp_angle_beam;
  float shw_sp_iso_angle;
  float shw_sp_n_vertex;
  float shw_sp_n_good_tracks;
  float shw_sp_E_indirect_max_energy;
  float shw_sp_flag_all_above;
  float shw_sp_min_dQ_dx_5;
  float shw_sp_n_other_vertex;
  float shw_sp_n_stem_size;
  float shw_sp_flag_stem_trajectory;
  float shw_sp_min_dis;
  tree3->SetBranchAddress("shw_sp_filled",&shw_sp_filled);
  tree3->SetBranchAddress("shw_sp_flag",&shw_sp_flag);
  tree3->SetBranchAddress("shw_sp_energy",&shw_sp_energy);
  tree3->SetBranchAddress("shw_sp_vec_dQ_dx_0",&shw_sp_vec_dQ_dx_0);
  tree3->SetBranchAddress("shw_sp_vec_dQ_dx_1",&shw_sp_vec_dQ_dx_1);
  tree3->SetBranchAddress("shw_sp_max_dQ_dx_sample",&shw_sp_max_dQ_dx_sample);
  tree3->SetBranchAddress("shw_sp_n_below_threshold",&shw_sp_n_below_threshold);
  tree3->SetBranchAddress("shw_sp_n_below_zero",&shw_sp_n_below_zero);
  tree3->SetBranchAddress("shw_sp_n_lowest",&shw_sp_n_lowest);
  tree3->SetBranchAddress("shw_sp_n_highest",&shw_sp_n_highest);
  tree3->SetBranchAddress("shw_sp_lowest_dQ_dx",&shw_sp_lowest_dQ_dx);
  tree3->SetBranchAddress("shw_sp_highest_dQ_dx",&shw_sp_highest_dQ_dx);
  tree3->SetBranchAddress("shw_sp_medium_dQ_dx",&shw_sp_medium_dQ_dx);
  tree3->SetBranchAddress("shw_sp_stem_length",&shw_sp_stem_length);
  tree3->SetBranchAddress("shw_sp_length_main",&shw_sp_length_main);
  tree3->SetBranchAddress("shw_sp_length_total",&shw_sp_length_total);
  tree3->SetBranchAddress("shw_sp_angle_beam",&shw_sp_angle_beam);
  tree3->SetBranchAddress("shw_sp_iso_angle",&shw_sp_iso_angle);
  tree3->SetBranchAddress("shw_sp_n_vertex",&shw_sp_n_vertex);
  tree3->SetBranchAddress("shw_sp_n_good_tracks",&shw_sp_n_good_tracks);
  tree3->SetBranchAddress("shw_sp_E_indirect_max_energy",&shw_sp_E_indirect_max_energy);
  tree3->SetBranchAddress("shw_sp_flag_all_above",&shw_sp_flag_all_above);
  tree3->SetBranchAddress("shw_sp_min_dQ_dx_5",&shw_sp_min_dQ_dx_5);
  tree3->SetBranchAddress("shw_sp_n_other_vertex",&shw_sp_n_other_vertex);
  tree3->SetBranchAddress("shw_sp_n_stem_size",&shw_sp_n_stem_size);
  tree3->SetBranchAddress("shw_sp_flag_stem_trajectory",&shw_sp_flag_stem_trajectory);
  tree3->SetBranchAddress("shw_sp_min_dis",&shw_sp_min_dis);

  float shw_sp_vec_median_dedx;
  float shw_sp_vec_mean_dedx;
  float shw_sp_vec_dQ_dx_2;
  float shw_sp_vec_dQ_dx_3;
  float shw_sp_vec_dQ_dx_4;
  float shw_sp_vec_dQ_dx_5;
  float shw_sp_vec_dQ_dx_6;
  float shw_sp_vec_dQ_dx_7;
  float shw_sp_vec_dQ_dx_8;
  float shw_sp_vec_dQ_dx_9;
  float shw_sp_vec_dQ_dx_10;
  float shw_sp_vec_dQ_dx_11;
  float shw_sp_vec_dQ_dx_12;
  float shw_sp_vec_dQ_dx_13;
  float shw_sp_vec_dQ_dx_14;
  float shw_sp_vec_dQ_dx_15;
  float shw_sp_vec_dQ_dx_16;
  float shw_sp_vec_dQ_dx_17;
  float shw_sp_vec_dQ_dx_18;
  float shw_sp_vec_dQ_dx_19;
  tree3->SetBranchAddress("shw_sp_vec_median_dedx",&shw_sp_vec_median_dedx);
  tree3->SetBranchAddress("shw_sp_vec_mean_dedx",&shw_sp_vec_mean_dedx);
  tree3->SetBranchAddress("shw_sp_vec_dQ_dx_2",&shw_sp_vec_dQ_dx_2);
  tree3->SetBranchAddress("shw_sp_vec_dQ_dx_3",&shw_sp_vec_dQ_dx_3);
  tree3->SetBranchAddress("shw_sp_vec_dQ_dx_4",&shw_sp_vec_dQ_dx_4);
  tree3->SetBranchAddress("shw_sp_vec_dQ_dx_5",&shw_sp_vec_dQ_dx_5);
  tree3->SetBranchAddress("shw_sp_vec_dQ_dx_6",&shw_sp_vec_dQ_dx_6);
  tree3->SetBranchAddress("shw_sp_vec_dQ_dx_7",&shw_sp_vec_dQ_dx_7);
  tree3->SetBranchAddress("shw_sp_vec_dQ_dx_8",&shw_sp_vec_dQ_dx_8);
  tree3->SetBranchAddress("shw_sp_vec_dQ_dx_9",&shw_sp_vec_dQ_dx_9);
  tree3->SetBranchAddress("shw_sp_vec_dQ_dx_10",&shw_sp_vec_dQ_dx_10);
  tree3->SetBranchAddress("shw_sp_vec_dQ_dx_11",&shw_sp_vec_dQ_dx_11);
  tree3->SetBranchAddress("shw_sp_vec_dQ_dx_12",&shw_sp_vec_dQ_dx_12);
  tree3->SetBranchAddress("shw_sp_vec_dQ_dx_13",&shw_sp_vec_dQ_dx_13);
  tree3->SetBranchAddress("shw_sp_vec_dQ_dx_14",&shw_sp_vec_dQ_dx_14);
  tree3->SetBranchAddress("shw_sp_vec_dQ_dx_15",&shw_sp_vec_dQ_dx_15);
  tree3->SetBranchAddress("shw_sp_vec_dQ_dx_16",&shw_sp_vec_dQ_dx_16);
  tree3->SetBranchAddress("shw_sp_vec_dQ_dx_17",&shw_sp_vec_dQ_dx_17);
  tree3->SetBranchAddress("shw_sp_vec_dQ_dx_18",&shw_sp_vec_dQ_dx_18);
  tree3->SetBranchAddress("shw_sp_vec_dQ_dx_19",&shw_sp_vec_dQ_dx_19);

  float shw_sp_pio_filled;
  float shw_sp_pio_flag;
  float shw_sp_pio_mip_id;
  float shw_sp_pio_flag_pio;
  float shw_sp_pio_1_flag;
  float shw_sp_pio_1_mass;
  float shw_sp_pio_1_pio_type;
  float shw_sp_pio_1_energy_1;
  float shw_sp_pio_1_energy_2;
  float shw_sp_pio_1_dis_1;
  float shw_sp_pio_1_dis_2;
  std::vector<float> *shw_sp_pio_2_v_flag = new std::vector<float>;
  std::vector<float> *shw_sp_pio_2_v_dis2 = new std::vector<float>;
  std::vector<float> *shw_sp_pio_2_v_angle2 = new std::vector<float>;
  std::vector<float> *shw_sp_pio_2_v_acc_length = new std::vector<float>;
  tree3->SetBranchAddress("shw_sp_pio_filled",&shw_sp_pio_filled);
  tree3->SetBranchAddress("shw_sp_pio_flag",&shw_sp_pio_flag);
  tree3->SetBranchAddress("shw_sp_pio_mip_id",&shw_sp_pio_mip_id);
  tree3->SetBranchAddress("shw_sp_pio_flag_pio",&shw_sp_pio_flag_pio);
  tree3->SetBranchAddress("shw_sp_pio_1_flag",&shw_sp_pio_1_flag);
  tree3->SetBranchAddress("shw_sp_pio_1_mass",&shw_sp_pio_1_mass);
  tree3->SetBranchAddress("shw_sp_pio_1_pio_type",&shw_sp_pio_1_pio_type);
  tree3->SetBranchAddress("shw_sp_pio_1_energy_1",&shw_sp_pio_1_energy_1);
  tree3->SetBranchAddress("shw_sp_pio_1_energy_2",&shw_sp_pio_1_energy_2);
  tree3->SetBranchAddress("shw_sp_pio_1_dis_1",&shw_sp_pio_1_dis_1);
  tree3->SetBranchAddress("shw_sp_pio_1_dis_2",&shw_sp_pio_1_dis_2);
  tree3->SetBranchAddress("shw_sp_pio_2_v_flag",&shw_sp_pio_2_v_flag);
  tree3->SetBranchAddress("shw_sp_pio_2_v_dis2",&shw_sp_pio_2_v_dis2);
  tree3->SetBranchAddress("shw_sp_pio_2_v_angle2",&shw_sp_pio_2_v_angle2);
  tree3->SetBranchAddress("shw_sp_pio_2_v_acc_length",&shw_sp_pio_2_v_acc_length);

  float shw_sp_lem_flag;
  float shw_sp_lem_shower_total_length;
  float shw_sp_lem_shower_main_length;
  float shw_sp_lem_n_3seg;
  float shw_sp_lem_e_charge;
  float shw_sp_lem_e_dQdx;
  float shw_sp_lem_shower_num_segs;
  float shw_sp_lem_shower_num_main_segs;
  tree3->SetBranchAddress("shw_sp_lem_flag",&shw_sp_lem_flag);
  tree3->SetBranchAddress("shw_sp_lem_shower_total_length",&shw_sp_lem_shower_total_length);
  tree3->SetBranchAddress("shw_sp_lem_shower_main_length",&shw_sp_lem_shower_main_length);
  tree3->SetBranchAddress("shw_sp_lem_n_3seg",&shw_sp_lem_n_3seg);
  tree3->SetBranchAddress("shw_sp_lem_e_charge",&shw_sp_lem_e_charge);
  tree3->SetBranchAddress("shw_sp_lem_e_dQdx",&shw_sp_lem_e_dQdx);
  tree3->SetBranchAddress("shw_sp_lem_shower_num_segs",&shw_sp_lem_shower_num_segs);
  tree3->SetBranchAddress("shw_sp_lem_shower_num_main_segs",&shw_sp_lem_shower_num_main_segs);

  float shw_sp_br_filled;
  float shw_sp_br1_flag;
  float shw_sp_br1_1_flag;
  float shw_sp_br1_1_shower_type;
  float shw_sp_br1_1_vtx_n_segs;
  float shw_sp_br1_1_energy;
  float shw_sp_br1_1_n_segs;
  float shw_sp_br1_1_flag_sg_topology;
  float shw_sp_br1_1_flag_sg_trajectory;
  float shw_sp_br1_1_sg_length;
  float shw_sp_br1_2_flag;
  float shw_sp_br1_2_energy;
  float shw_sp_br1_2_n_connected;
  float shw_sp_br1_2_max_length;
  float shw_sp_br1_2_n_connected_1;
  float shw_sp_br1_2_vtx_n_segs;
  float shw_sp_br1_2_n_shower_segs;
  float shw_sp_br1_2_max_length_ratio;
  float shw_sp_br1_2_shower_length;
  float shw_sp_br1_3_flag;
  float shw_sp_br1_3_energy;
  float shw_sp_br1_3_n_connected_p;
  float shw_sp_br1_3_max_length_p;
  float shw_sp_br1_3_n_shower_segs;
  float shw_sp_br1_3_flag_sg_topology;
  float shw_sp_br1_3_flag_sg_trajectory;
  float shw_sp_br1_3_n_shower_main_segs;
  float shw_sp_br1_3_sg_length;
  tree3->SetBranchAddress("shw_sp_br_filled",&shw_sp_br_filled);
  tree3->SetBranchAddress("shw_sp_br1_flag",&shw_sp_br1_flag);
  tree3->SetBranchAddress("shw_sp_br1_1_flag",&shw_sp_br1_1_flag);
  tree3->SetBranchAddress("shw_sp_br1_1_shower_type",&shw_sp_br1_1_shower_type);
  tree3->SetBranchAddress("shw_sp_br1_1_vtx_n_segs",&shw_sp_br1_1_vtx_n_segs);
  tree3->SetBranchAddress("shw_sp_br1_1_energy",&shw_sp_br1_1_energy);
  tree3->SetBranchAddress("shw_sp_br1_1_n_segs",&shw_sp_br1_1_n_segs);
  tree3->SetBranchAddress("shw_sp_br1_1_flag_sg_topology",&shw_sp_br1_1_flag_sg_topology);
  tree3->SetBranchAddress("shw_sp_br1_1_flag_sg_trajectory",&shw_sp_br1_1_flag_sg_trajectory);
  tree3->SetBranchAddress("shw_sp_br1_1_sg_length",&shw_sp_br1_1_sg_length);
  tree3->SetBranchAddress("shw_sp_br1_2_flag",&shw_sp_br1_2_flag);
  tree3->SetBranchAddress("shw_sp_br1_2_energy",&shw_sp_br1_2_energy);
  tree3->SetBranchAddress("shw_sp_br1_2_n_connected",&shw_sp_br1_2_n_connected);
  tree3->SetBranchAddress("shw_sp_br1_2_max_length",&shw_sp_br1_2_max_length);
  tree3->SetBranchAddress("shw_sp_br1_2_n_connected_1",&shw_sp_br1_2_n_connected_1);
  tree3->SetBranchAddress("shw_sp_br1_2_vtx_n_segs",&shw_sp_br1_2_vtx_n_segs);
  tree3->SetBranchAddress("shw_sp_br1_2_n_shower_segs",&shw_sp_br1_2_n_shower_segs);
  tree3->SetBranchAddress("shw_sp_br1_2_max_length_ratio",&shw_sp_br1_2_max_length_ratio);
  tree3->SetBranchAddress("shw_sp_br1_2_shower_length",&shw_sp_br1_2_shower_length);
  tree3->SetBranchAddress("shw_sp_br1_3_flag",&shw_sp_br1_3_flag);
  tree3->SetBranchAddress("shw_sp_br1_3_energy",&shw_sp_br1_3_energy);
  tree3->SetBranchAddress("shw_sp_br1_3_n_connected_p",&shw_sp_br1_3_n_connected_p);
  tree3->SetBranchAddress("shw_sp_br1_3_max_length_p",&shw_sp_br1_3_max_length_p);
  tree3->SetBranchAddress("shw_sp_br1_3_n_shower_segs",&shw_sp_br1_3_n_shower_segs);
  tree3->SetBranchAddress("shw_sp_br1_3_flag_sg_topology",&shw_sp_br1_3_flag_sg_topology);
  tree3->SetBranchAddress("shw_sp_br1_3_flag_sg_trajectory",&shw_sp_br1_3_flag_sg_trajectory);
  tree3->SetBranchAddress("shw_sp_br1_3_n_shower_main_segs",&shw_sp_br1_3_n_shower_main_segs);
  tree3->SetBranchAddress("shw_sp_br1_3_sg_length",&shw_sp_br1_3_sg_length);


  float shw_sp_br2_flag;
  float shw_sp_br2_flag_single_shower;
  float shw_sp_br2_num_valid_tracks;
  float shw_sp_br2_energy;
  float shw_sp_br2_angle1;
  float shw_sp_br2_angle2;
  float shw_sp_br2_angle;
  float shw_sp_br2_angle3;
  float shw_sp_br2_n_shower_main_segs;
  float shw_sp_br2_max_angle;
  float shw_sp_br2_sg_length;
  float shw_sp_br2_flag_sg_trajectory;
  tree3->SetBranchAddress("shw_sp_br2_flag",&shw_sp_br2_flag);
  tree3->SetBranchAddress("shw_sp_br2_flag_single_shower",&shw_sp_br2_flag_single_shower);
  tree3->SetBranchAddress("shw_sp_br2_num_valid_tracks",&shw_sp_br2_num_valid_tracks);
  tree3->SetBranchAddress("shw_sp_br2_energy",&shw_sp_br2_energy);
  tree3->SetBranchAddress("shw_sp_br2_angle1",&shw_sp_br2_angle1);
  tree3->SetBranchAddress("shw_sp_br2_angle2",&shw_sp_br2_angle2);
  tree3->SetBranchAddress("shw_sp_br2_angle",&shw_sp_br2_angle);
  tree3->SetBranchAddress("shw_sp_br2_angle3",&shw_sp_br2_angle3);
  tree3->SetBranchAddress("shw_sp_br2_n_shower_main_segs",&shw_sp_br2_n_shower_main_segs);
  tree3->SetBranchAddress("shw_sp_br2_max_angle",&shw_sp_br2_max_angle);
  tree3->SetBranchAddress("shw_sp_br2_sg_length",&shw_sp_br2_sg_length);
  tree3->SetBranchAddress("shw_sp_br2_flag_sg_trajectory",&shw_sp_br2_flag_sg_trajectory);


  float shw_sp_br3_flag;
  float shw_sp_br3_1_flag;
  float shw_sp_br3_1_energy;
  float shw_sp_br3_1_n_shower_segments;
  float shw_sp_br3_1_sg_flag_trajectory;
  float shw_sp_br3_1_sg_direct_length;
  float shw_sp_br3_1_sg_length;
  float shw_sp_br3_1_total_main_length;
  float shw_sp_br3_1_total_length;
  float shw_sp_br3_1_iso_angle;
  float shw_sp_br3_1_sg_flag_topology;
  float shw_sp_br3_2_flag;
  float shw_sp_br3_2_n_ele;
  float shw_sp_br3_2_n_other;
  float shw_sp_br3_2_energy;
  float shw_sp_br3_2_total_main_length;
  float shw_sp_br3_2_total_length;
  float shw_sp_br3_2_other_fid;
  std::vector<float> *shw_sp_br3_3_v_flag = new std::vector<float>;
  std::vector<float> *shw_sp_br3_3_v_energy = new std::vector<float>;
  std::vector<float> *shw_sp_br3_3_v_angle = new std::vector<float>;
  std::vector<float> *shw_sp_br3_3_v_dir_length = new std::vector<float>;
  std::vector<float> *shw_sp_br3_3_v_length = new std::vector<float>;
  float shw_sp_br3_4_flag;
  float shw_sp_br3_4_acc_length;
  float shw_sp_br3_4_total_length;
  float shw_sp_br3_4_energy;
  std::vector<float> *shw_sp_br3_5_v_flag = new std::vector<float>;
  std::vector<float> *shw_sp_br3_5_v_dir_length = new std::vector<float>;
  std::vector<float> *shw_sp_br3_5_v_total_length = new std::vector<float>;
  std::vector<float> *shw_sp_br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  std::vector<float> *shw_sp_br3_5_v_n_seg = new std::vector<float>;
  std::vector<float> *shw_sp_br3_5_v_angle = new std::vector<float>;
  std::vector<float> *shw_sp_br3_5_v_sg_length = new std::vector<float>;
  std::vector<float> *shw_sp_br3_5_v_energy = new std::vector<float>;
  std::vector<float> *shw_sp_br3_5_v_n_main_segs = new std::vector<float>;
  std::vector<float> *shw_sp_br3_5_v_n_segs = new std::vector<float>;
  std::vector<float> *shw_sp_br3_5_v_shower_main_length = new std::vector<float>;
  std::vector<float> *shw_sp_br3_5_v_shower_total_length = new std::vector<float>;
  std::vector<float> *shw_sp_br3_6_v_flag = new std::vector<float>;
  std::vector<float> *shw_sp_br3_6_v_angle = new std::vector<float>;
  std::vector<float> *shw_sp_br3_6_v_angle1 = new std::vector<float>;
  std::vector<float> *shw_sp_br3_6_v_flag_shower_trajectory = new std::vector<float>;
  std::vector<float> *shw_sp_br3_6_v_direct_length = new std::vector<float>;
  std::vector<float> *shw_sp_br3_6_v_length = new std::vector<float>;
  std::vector<float> *shw_sp_br3_6_v_n_other_vtx_segs = new std::vector<float>;
  std::vector<float> *shw_sp_br3_6_v_energy = new std::vector<float>;
  float shw_sp_br3_7_flag;
  float shw_sp_br3_7_energy;
  float shw_sp_br3_7_min_angle;
  float shw_sp_br3_7_sg_length;
  float shw_sp_br3_7_shower_main_length;
  float shw_sp_br3_8_flag;
  float shw_sp_br3_8_max_dQ_dx;
  float shw_sp_br3_8_energy;
  float shw_sp_br3_8_n_main_segs;
  float shw_sp_br3_8_shower_main_length;
  float shw_sp_br3_8_shower_length;
  tree3->SetBranchAddress("shw_sp_br3_flag",&shw_sp_br3_flag);
  tree3->SetBranchAddress("shw_sp_br3_1_flag",&shw_sp_br3_1_flag);
  tree3->SetBranchAddress("shw_sp_br3_1_energy",&shw_sp_br3_1_energy);
  tree3->SetBranchAddress("shw_sp_br3_1_n_shower_segments",&shw_sp_br3_1_n_shower_segments);
  tree3->SetBranchAddress("shw_sp_br3_1_sg_flag_trajectory",&shw_sp_br3_1_sg_flag_trajectory);
  tree3->SetBranchAddress("shw_sp_br3_1_sg_direct_length",&shw_sp_br3_1_sg_direct_length);
  tree3->SetBranchAddress("shw_sp_br3_1_sg_length",&shw_sp_br3_1_sg_length);
  tree3->SetBranchAddress("shw_sp_br3_1_total_main_length",&shw_sp_br3_1_total_main_length);
  tree3->SetBranchAddress("shw_sp_br3_1_total_length",&shw_sp_br3_1_total_length);
  tree3->SetBranchAddress("shw_sp_br3_1_iso_angle",&shw_sp_br3_1_iso_angle);
  tree3->SetBranchAddress("shw_sp_br3_1_sg_flag_topology",&shw_sp_br3_1_sg_flag_topology);
  tree3->SetBranchAddress("shw_sp_br3_2_flag",&shw_sp_br3_2_flag);
  tree3->SetBranchAddress("shw_sp_br3_2_n_ele",&shw_sp_br3_2_n_ele);
  tree3->SetBranchAddress("shw_sp_br3_2_n_other",&shw_sp_br3_2_n_other);
  tree3->SetBranchAddress("shw_sp_br3_2_energy",&shw_sp_br3_2_energy);
  tree3->SetBranchAddress("shw_sp_br3_2_total_main_length",&shw_sp_br3_2_total_main_length);
  tree3->SetBranchAddress("shw_sp_br3_2_total_length",&shw_sp_br3_2_total_length);
  tree3->SetBranchAddress("shw_sp_br3_2_other_fid",&shw_sp_br3_2_other_fid);
  tree3->SetBranchAddress("shw_sp_br3_3_v_flag",&shw_sp_br3_3_v_flag);
  tree3->SetBranchAddress("shw_sp_br3_3_v_energy",&shw_sp_br3_3_v_energy);
  tree3->SetBranchAddress("shw_sp_br3_3_v_angle",&shw_sp_br3_3_v_angle);
  tree3->SetBranchAddress("shw_sp_br3_3_v_dir_length",&shw_sp_br3_3_v_dir_length);
  tree3->SetBranchAddress("shw_sp_br3_3_v_length",&shw_sp_br3_3_v_length);
  tree3->SetBranchAddress("shw_sp_br3_4_flag", &shw_sp_br3_4_flag);
  tree3->SetBranchAddress("shw_sp_br3_4_acc_length", &shw_sp_br3_4_acc_length);
  tree3->SetBranchAddress("shw_sp_br3_4_total_length", &shw_sp_br3_4_total_length);
  tree3->SetBranchAddress("shw_sp_br3_4_energy", &shw_sp_br3_4_energy);
  tree3->SetBranchAddress("shw_sp_br3_5_v_flag", &shw_sp_br3_5_v_flag);
  tree3->SetBranchAddress("shw_sp_br3_5_v_dir_length", &shw_sp_br3_5_v_dir_length);
  tree3->SetBranchAddress("shw_sp_br3_5_v_total_length", &shw_sp_br3_5_v_total_length);
  tree3->SetBranchAddress("shw_sp_br3_5_v_flag_avoid_muon_check", &shw_sp_br3_5_v_flag_avoid_muon_check);
  tree3->SetBranchAddress("shw_sp_br3_5_v_n_seg", &shw_sp_br3_5_v_n_seg);
  tree3->SetBranchAddress("shw_sp_br3_5_v_angle", &shw_sp_br3_5_v_angle);
  tree3->SetBranchAddress("shw_sp_br3_5_v_sg_length", &shw_sp_br3_5_v_sg_length);
  tree3->SetBranchAddress("shw_sp_br3_5_v_energy", &shw_sp_br3_5_v_energy);
  tree3->SetBranchAddress("shw_sp_br3_5_v_n_main_segs", &shw_sp_br3_5_v_n_main_segs);
  tree3->SetBranchAddress("shw_sp_br3_5_v_n_segs", &shw_sp_br3_5_v_n_segs);
  tree3->SetBranchAddress("shw_sp_br3_5_v_shower_main_length", &shw_sp_br3_5_v_shower_main_length);
  tree3->SetBranchAddress("shw_sp_br3_5_v_shower_total_length", &shw_sp_br3_5_v_shower_total_length);
  tree3->SetBranchAddress("shw_sp_br3_6_v_flag",&shw_sp_br3_6_v_flag);
  tree3->SetBranchAddress("shw_sp_br3_6_v_angle",&shw_sp_br3_6_v_angle);
  tree3->SetBranchAddress("shw_sp_br3_6_v_angle1",&shw_sp_br3_6_v_angle1);
  tree3->SetBranchAddress("shw_sp_br3_6_v_flag_shower_trajectory",&shw_sp_br3_6_v_flag_shower_trajectory);
  tree3->SetBranchAddress("shw_sp_br3_6_v_direct_length",&shw_sp_br3_6_v_direct_length);
  tree3->SetBranchAddress("shw_sp_br3_6_v_length",&shw_sp_br3_6_v_length);
  tree3->SetBranchAddress("shw_sp_br3_6_v_n_other_vtx_segs",&shw_sp_br3_6_v_n_other_vtx_segs);
  tree3->SetBranchAddress("shw_sp_br3_6_v_energy",&shw_sp_br3_6_v_energy);
  tree3->SetBranchAddress("shw_sp_br3_7_flag",&shw_sp_br3_7_flag);
  tree3->SetBranchAddress("shw_sp_br3_7_energy",&shw_sp_br3_7_energy);
  tree3->SetBranchAddress("shw_sp_br3_7_min_angle",&shw_sp_br3_7_min_angle);
  tree3->SetBranchAddress("shw_sp_br3_7_sg_length",&shw_sp_br3_7_sg_length);
  tree3->SetBranchAddress("shw_sp_br3_7_main_length",&shw_sp_br3_7_shower_main_length);
  tree3->SetBranchAddress("shw_sp_br3_8_flag",&shw_sp_br3_8_flag);
  tree3->SetBranchAddress("shw_sp_br3_8_max_dQ_dx",&shw_sp_br3_8_max_dQ_dx);
  tree3->SetBranchAddress("shw_sp_br3_8_energy",&shw_sp_br3_8_energy);
  tree3->SetBranchAddress("shw_sp_br3_8_n_main_segs",&shw_sp_br3_8_n_main_segs);
  tree3->SetBranchAddress("shw_sp_br3_8_shower_main_length",&shw_sp_br3_8_shower_main_length);
  tree3->SetBranchAddress("shw_sp_br3_8_shower_length",&shw_sp_br3_8_shower_length);


  float shw_sp_br4_flag;
  float shw_sp_br4_1_flag;
  float shw_sp_br4_1_shower_main_length;
  float shw_sp_br4_1_shower_total_length;
  float shw_sp_br4_1_min_dis;
  float shw_sp_br4_1_energy;
  float shw_sp_br4_1_flag_avoid_muon_check;
  float shw_sp_br4_1_n_vtx_segs;
  float shw_sp_br4_1_n_main_segs;
  float shw_sp_br4_2_flag;
  float shw_sp_br4_2_ratio_45;
  float shw_sp_br4_2_ratio_35;
  float shw_sp_br4_2_ratio_25;
  float shw_sp_br4_2_ratio_15;
  float shw_sp_br4_2_energy;
  float shw_sp_br4_2_ratio1_45;
  float shw_sp_br4_2_ratio1_35;
  float shw_sp_br4_2_ratio1_25;
  float shw_sp_br4_2_ratio1_15;
  float shw_sp_br4_2_iso_angle;
  float shw_sp_br4_2_iso_angle1;
  float shw_sp_br4_2_angle;
  tree3->SetBranchAddress("shw_sp_br4_flag", &shw_sp_br4_flag);
  tree3->SetBranchAddress("shw_sp_br4_1_flag", &shw_sp_br4_1_flag);
  tree3->SetBranchAddress("shw_sp_br4_1_shower_main_length", &shw_sp_br4_1_shower_main_length);
  tree3->SetBranchAddress("shw_sp_br4_1_shower_total_length", &shw_sp_br4_1_shower_total_length);
  tree3->SetBranchAddress("shw_sp_br4_1_min_dis", &shw_sp_br4_1_min_dis);
  tree3->SetBranchAddress("shw_sp_br4_1_energy", &shw_sp_br4_1_energy);
  tree3->SetBranchAddress("shw_sp_br4_1_flag_avoid_muon_check", &shw_sp_br4_1_flag_avoid_muon_check);
  tree3->SetBranchAddress("shw_sp_br4_1_n_vtx_segs", &shw_sp_br4_1_n_vtx_segs);
  tree3->SetBranchAddress("shw_sp_br4_1_n_main_segs", &shw_sp_br4_1_n_main_segs);
  tree3->SetBranchAddress("shw_sp_br4_2_flag", &shw_sp_br4_2_flag);
  tree3->SetBranchAddress("shw_sp_br4_2_ratio_45", &shw_sp_br4_2_ratio_45);
  tree3->SetBranchAddress("shw_sp_br4_2_ratio_35", &shw_sp_br4_2_ratio_35);
  tree3->SetBranchAddress("shw_sp_br4_2_ratio_25", &shw_sp_br4_2_ratio_25);
  tree3->SetBranchAddress("shw_sp_br4_2_ratio_15", &shw_sp_br4_2_ratio_15);
  tree3->SetBranchAddress("shw_sp_br4_2_energy",   &shw_sp_br4_2_energy);
  tree3->SetBranchAddress("shw_sp_br4_2_ratio1_45", &shw_sp_br4_2_ratio1_45);
  tree3->SetBranchAddress("shw_sp_br4_2_ratio1_35", &shw_sp_br4_2_ratio1_35);
  tree3->SetBranchAddress("shw_sp_br4_2_ratio1_25", &shw_sp_br4_2_ratio1_25);
  tree3->SetBranchAddress("shw_sp_br4_2_ratio1_15", &shw_sp_br4_2_ratio1_15);
  tree3->SetBranchAddress("shw_sp_br4_2_iso_angle", &shw_sp_br4_2_iso_angle);
  tree3->SetBranchAddress("shw_sp_br4_2_iso_angle1", &shw_sp_br4_2_iso_angle1);
  tree3->SetBranchAddress("shw_sp_br4_2_angle", &shw_sp_br4_2_angle);

  float shw_sp_hol_flag;
  float shw_sp_hol_1_flag;
  float shw_sp_hol_1_n_valid_tracks;
  float shw_sp_hol_1_min_angle;
  float shw_sp_hol_1_energy;
  float shw_sp_hol_1_flag_all_shower;
  float shw_sp_hol_1_min_length;
  float shw_sp_hol_2_flag;
  float shw_sp_hol_2_min_angle;
  float shw_sp_hol_2_medium_dQ_dx;
  float shw_sp_hol_2_ncount;
  float shw_sp_hol_2_energy;
  tree3->SetBranchAddress("shw_sp_hol_flag", &shw_sp_hol_flag);
  tree3->SetBranchAddress("shw_sp_hol_1_flag", &shw_sp_hol_1_flag);
  tree3->SetBranchAddress("shw_sp_hol_1_n_valid_tracks", &shw_sp_hol_1_n_valid_tracks);
  tree3->SetBranchAddress("shw_sp_hol_1_min_angle", &shw_sp_hol_1_min_angle);
  tree3->SetBranchAddress("shw_sp_hol_1_energy", &shw_sp_hol_1_energy);
  tree3->SetBranchAddress("shw_sp_hol_1_flag_all_shower", &shw_sp_hol_1_flag_all_shower);
  tree3->SetBranchAddress("shw_sp_hol_1_min_length", &shw_sp_hol_1_min_length);
  tree3->SetBranchAddress("shw_sp_hol_2_flag", &shw_sp_hol_2_flag);
  tree3->SetBranchAddress("shw_sp_hol_2_min_angle", &shw_sp_hol_2_min_angle);
  tree3->SetBranchAddress("shw_sp_hol_2_medium_dQ_dx", &shw_sp_hol_2_medium_dQ_dx);
  tree3->SetBranchAddress("shw_sp_hol_2_ncount", &shw_sp_hol_2_ncount);
  tree3->SetBranchAddress("shw_sp_hol_2_energy", &shw_sp_hol_2_energy);


  float shw_sp_lol_flag;
  std::vector<float> *shw_sp_lol_1_v_flag = new std::vector<float>;
  std::vector<float> *shw_sp_lol_1_v_energy = new std::vector<float>;
  std::vector<float> *shw_sp_lol_1_v_vtx_n_segs = new std::vector<float>;
  std::vector<float> *shw_sp_lol_1_v_nseg = new std::vector<float>;
  std::vector<float> *shw_sp_lol_1_v_angle = new std::vector<float>;
  std::vector<float> *shw_sp_lol_2_v_flag = new std::vector<float>;
  std::vector<float> *shw_sp_lol_2_v_length = new std::vector<float>;
  std::vector<float> *shw_sp_lol_2_v_angle = new std::vector<float>;
  std::vector<float> *shw_sp_lol_2_v_type = new std::vector<float>;
  std::vector<float> *shw_sp_lol_2_v_vtx_n_segs = new std::vector<float>;
  std::vector<float> *shw_sp_lol_2_v_energy = new std::vector<float>;
  std::vector<float> *shw_sp_lol_2_v_shower_main_length = new std::vector<float>;
  std::vector<float> *shw_sp_lol_2_v_flag_dir_weak = new std::vector<float>;
  float shw_sp_lol_3_flag;
  float shw_sp_lol_3_angle_beam;
  float shw_sp_lol_3_n_valid_tracks;
  float shw_sp_lol_3_min_angle;
  float shw_sp_lol_3_vtx_n_segs;
  float shw_sp_lol_3_energy;
  float shw_sp_lol_3_shower_main_length;
  float shw_sp_lol_3_n_out;
  float shw_sp_lol_3_n_sum;
  tree3->SetBranchAddress("shw_sp_lol_flag",&shw_sp_lol_flag);
  tree3->SetBranchAddress("shw_sp_lol_1_v_flag",&shw_sp_lol_1_v_flag);
  tree3->SetBranchAddress("shw_sp_lol_1_v_energy",&shw_sp_lol_1_v_energy);
  tree3->SetBranchAddress("shw_sp_lol_1_v_vtx_n_segs",&shw_sp_lol_1_v_vtx_n_segs);
  tree3->SetBranchAddress("shw_sp_lol_1_v_nseg",&shw_sp_lol_1_v_nseg);
  tree3->SetBranchAddress("shw_sp_lol_1_v_angle",&shw_sp_lol_1_v_angle);
  tree3->SetBranchAddress("shw_sp_lol_2_v_flag",&shw_sp_lol_2_v_flag);
  tree3->SetBranchAddress("shw_sp_lol_2_v_length",&shw_sp_lol_2_v_length);
  tree3->SetBranchAddress("shw_sp_lol_2_v_angle",&shw_sp_lol_2_v_angle);
  tree3->SetBranchAddress("shw_sp_lol_2_v_type",&shw_sp_lol_2_v_type);
  tree3->SetBranchAddress("shw_sp_lol_2_v_vtx_n_segs",&shw_sp_lol_2_v_vtx_n_segs);
  tree3->SetBranchAddress("shw_sp_lol_2_v_energy",&shw_sp_lol_2_v_energy);
  tree3->SetBranchAddress("shw_sp_lol_2_v_shower_main_length",&shw_sp_lol_2_v_shower_main_length);
  tree3->SetBranchAddress("shw_sp_lol_2_v_flag_dir_weak",&shw_sp_lol_2_v_flag_dir_weak);
  tree3->SetBranchAddress("shw_sp_lol_3_flag",&shw_sp_lol_3_flag);
  tree3->SetBranchAddress("shw_sp_lol_3_angle_beam",&shw_sp_lol_3_angle_beam);
  tree3->SetBranchAddress("shw_sp_lol_3_n_valid_tracks",&shw_sp_lol_3_n_valid_tracks);
  tree3->SetBranchAddress("shw_sp_lol_3_min_angle",&shw_sp_lol_3_min_angle);
  tree3->SetBranchAddress("shw_sp_lol_3_vtx_n_segs",&shw_sp_lol_3_vtx_n_segs);
  tree3->SetBranchAddress("shw_sp_lol_3_energy",&shw_sp_lol_3_energy);
  tree3->SetBranchAddress("shw_sp_lol_3_shower_main_length",&shw_sp_lol_3_shower_main_length);
  tree3->SetBranchAddress("shw_sp_lol_3_n_out",&shw_sp_lol_3_n_out);
  tree3->SetBranchAddress("shw_sp_lol_3_n_sum",&shw_sp_lol_3_n_sum);



	  float cosmic_filled;
	  float cosmic_flag;
	  float cosmic_n_solid_tracks;
	  float cosmic_energy_main_showers;
	  float cosmic_energy_direct_showers;
	  float cosmic_energy_indirect_showers;
	  float cosmic_n_direct_showers;
	  float cosmic_n_indirect_showers;
	  float cosmic_n_main_showers;
	  tree3->SetBranchAddress("cosmic_filled",&cosmic_filled);
	  tree3->SetBranchAddress("cosmic_flag",&cosmic_flag);
	  tree3->SetBranchAddress("cosmic_energy_main_showers",&cosmic_energy_main_showers);
	  tree3->SetBranchAddress("cosmic_energy_indirect_showers",&cosmic_energy_indirect_showers);
	  tree3->SetBranchAddress("cosmic_energy_direct_showers",&cosmic_energy_direct_showers);
	  tree3->SetBranchAddress("cosmic_n_direct_showers",&cosmic_n_direct_showers);
	  tree3->SetBranchAddress("cosmic_n_indirect_showers",&cosmic_n_indirect_showers);
	  tree3->SetBranchAddress("cosmic_n_solid_tracks",&cosmic_n_solid_tracks);
	  tree3->SetBranchAddress("cosmic_n_main_showers",&cosmic_n_main_showers);


	  float gap_filled;
	  float gap_flag;
	  float gap_flag_prolong_u;
	  float gap_flag_prolong_v;
	  float gap_flag_prolong_w;
	  float gap_flag_parallel;
	  float gap_n_points;
	  float gap_n_bad;
	  float gap_energy;
	  float gap_num_valid_tracks;
	  float gap_flag_single_shower;
	  tree3->SetBranchAddress("gap_flag",&gap_flag);
	  tree3->SetBranchAddress("gap_flag_prolong_u",&gap_flag_prolong_u);
	  tree3->SetBranchAddress("gap_flag_prolong_v",&gap_flag_prolong_v);
	  tree3->SetBranchAddress("gap_flag_prolong_w",&gap_flag_prolong_w);
	  tree3->SetBranchAddress("gap_flag_parallel",&gap_flag_parallel);
	  tree3->SetBranchAddress("gap_n_points",&gap_n_points);
	  tree3->SetBranchAddress("gap_n_bad",&gap_n_bad);
	  tree3->SetBranchAddress("gap_energy",&gap_energy);
	  tree3->SetBranchAddress("gap_num_valid_tracks",&gap_num_valid_tracks);
	  tree3->SetBranchAddress("gap_flag_single_shower",&gap_flag_single_shower);
	  tree3->SetBranchAddress("gap_filled",&gap_filled);


	  float mip_quality_filled;
	  float mip_quality_flag;
	  float mip_quality_energy;
	  float mip_quality_overlap;
	  float mip_quality_n_showers;
	  float mip_quality_n_tracks;
	  float mip_quality_flag_inside_pi0;
	  float mip_quality_n_pi0_showers;
	  float mip_quality_shortest_length;
	  float mip_quality_acc_length;
	  float mip_quality_shortest_angle;
	  float mip_quality_flag_proton;
	  tree3->SetBranchAddress("mip_quality_filled",&mip_quality_filled);
	  tree3->SetBranchAddress("mip_quality_flag",&mip_quality_flag);
	  tree3->SetBranchAddress("mip_quality_energy",&mip_quality_energy);
	  tree3->SetBranchAddress("mip_quality_overlap",&mip_quality_overlap);
	  tree3->SetBranchAddress("mip_quality_n_showers",&mip_quality_n_showers);
	  tree3->SetBranchAddress("mip_quality_n_tracks",&mip_quality_n_tracks);
	  tree3->SetBranchAddress("mip_quality_flag_inside_pi0",&mip_quality_flag_inside_pi0);
	  tree3->SetBranchAddress("mip_quality_n_pi0_showers",&mip_quality_n_pi0_showers);
	  tree3->SetBranchAddress("mip_quality_shortest_length",&mip_quality_shortest_length);
	  tree3->SetBranchAddress("mip_quality_acc_length",&mip_quality_acc_length);
	  tree3->SetBranchAddress("mip_quality_shortest_angle",&mip_quality_shortest_angle);
	  tree3->SetBranchAddress("mip_quality_flag_proton",&mip_quality_flag_proton);


	  float mip_filled;
	  float mip_flag;
	  float mip_energy;
	  float mip_n_end_reduction;
	  float mip_n_first_mip;
	  float mip_n_first_non_mip;
    float mip_n_first_non_mip_1;
    float mip_n_first_non_mip_2;
    float mip_vec_dQ_dx_0;
    float mip_vec_dQ_dx_1;
    float mip_max_dQ_dx_sample;
    float mip_n_below_threshold;
    float mip_n_below_zero;
    float mip_n_lowest;
    float mip_n_highest;
    float mip_lowest_dQ_dx;
    float mip_highest_dQ_dx;
    float mip_medium_dQ_dx;
    float mip_stem_length;
    float mip_length_main;
    float mip_length_total;
    float mip_angle_beam;
    float mip_iso_angle;
    float mip_n_vertex;
    float mip_n_good_tracks;
    float mip_E_indirect_max_energy;
    float mip_flag_all_above;
    float mip_min_dQ_dx_5;
    float mip_n_other_vertex;
    float mip_n_stem_size;
    float mip_flag_stem_trajectory;
    float mip_min_dis;
	  tree3->SetBranchAddress("mip_filled",&mip_filled);
	  tree3->SetBranchAddress("mip_flag",&mip_flag);
	  tree3->SetBranchAddress("mip_energy",&mip_energy);
	  tree3->SetBranchAddress("mip_n_end_reduction",&mip_n_end_reduction);
	  tree3->SetBranchAddress("mip_n_first_mip",&mip_n_first_mip);
	  tree3->SetBranchAddress("mip_n_first_non_mip",&mip_n_first_non_mip);
	  tree3->SetBranchAddress("mip_n_first_non_mip_1",&mip_n_first_non_mip_1);
	  tree3->SetBranchAddress("mip_n_first_non_mip_2",&mip_n_first_non_mip_2);
	  tree3->SetBranchAddress("mip_vec_dQ_dx_0",&mip_vec_dQ_dx_0);
	  tree3->SetBranchAddress("mip_vec_dQ_dx_1",&mip_vec_dQ_dx_1);
	  tree3->SetBranchAddress("mip_max_dQ_dx_sample",&mip_max_dQ_dx_sample);
	  tree3->SetBranchAddress("mip_n_below_threshold",&mip_n_below_threshold);
	  tree3->SetBranchAddress("mip_n_below_zero",&mip_n_below_zero);
	  tree3->SetBranchAddress("mip_n_lowest",&mip_n_lowest);
	  tree3->SetBranchAddress("mip_n_highest",&mip_n_highest);
	  tree3->SetBranchAddress("mip_lowest_dQ_dx",&mip_lowest_dQ_dx);
	  tree3->SetBranchAddress("mip_highest_dQ_dx",&mip_highest_dQ_dx);
	  tree3->SetBranchAddress("mip_medium_dQ_dx",&mip_medium_dQ_dx);
	  tree3->SetBranchAddress("mip_stem_length",&mip_stem_length);
	  tree3->SetBranchAddress("mip_length_main",&mip_length_main);
	  tree3->SetBranchAddress("mip_length_total",&mip_length_total);
	  tree3->SetBranchAddress("mip_angle_beam",&mip_angle_beam);
	  tree3->SetBranchAddress("mip_iso_angle",&mip_iso_angle);
	  tree3->SetBranchAddress("mip_n_vertex",&mip_n_vertex);
	  tree3->SetBranchAddress("mip_n_good_tracks",&mip_n_good_tracks);
	  tree3->SetBranchAddress("mip_E_indirect_max_energy",&mip_E_indirect_max_energy);
	  tree3->SetBranchAddress("mip_flag_all_above",&mip_flag_all_above);
	  tree3->SetBranchAddress("mip_min_dQ_dx_5",&mip_min_dQ_dx_5);
	  tree3->SetBranchAddress("mip_n_other_vertex",&mip_n_other_vertex);
	  tree3->SetBranchAddress("mip_n_stem_size",&mip_n_stem_size);
	  tree3->SetBranchAddress("mip_flag_stem_trajectory",&mip_flag_stem_trajectory);
	  tree3->SetBranchAddress("mip_min_dis",&mip_min_dis);


	  float mip_vec_dQ_dx_2;
	  float mip_vec_dQ_dx_3;
	  float mip_vec_dQ_dx_4;
	  float mip_vec_dQ_dx_5;
	  float mip_vec_dQ_dx_6;
	  float mip_vec_dQ_dx_7;
	  float mip_vec_dQ_dx_8;
	  float mip_vec_dQ_dx_9;
	  float mip_vec_dQ_dx_10;
	  float mip_vec_dQ_dx_11;
	  float mip_vec_dQ_dx_12;
	  float mip_vec_dQ_dx_13;
	  float mip_vec_dQ_dx_14;
	  float mip_vec_dQ_dx_15;
	  float mip_vec_dQ_dx_16;
	  float mip_vec_dQ_dx_17;
	  float mip_vec_dQ_dx_18;
	  float mip_vec_dQ_dx_19;
	  tree3->SetBranchAddress("mip_vec_dQ_dx_2",&mip_vec_dQ_dx_2);
	  tree3->SetBranchAddress("mip_vec_dQ_dx_3",&mip_vec_dQ_dx_3);
	  tree3->SetBranchAddress("mip_vec_dQ_dx_4",&mip_vec_dQ_dx_4);
	  tree3->SetBranchAddress("mip_vec_dQ_dx_5",&mip_vec_dQ_dx_5);
	  tree3->SetBranchAddress("mip_vec_dQ_dx_6",&mip_vec_dQ_dx_6);
	  tree3->SetBranchAddress("mip_vec_dQ_dx_7",&mip_vec_dQ_dx_7);
	  tree3->SetBranchAddress("mip_vec_dQ_dx_8",&mip_vec_dQ_dx_8);
	  tree3->SetBranchAddress("mip_vec_dQ_dx_9",&mip_vec_dQ_dx_9);
	  tree3->SetBranchAddress("mip_vec_dQ_dx_10",&mip_vec_dQ_dx_10);
	  tree3->SetBranchAddress("mip_vec_dQ_dx_11",&mip_vec_dQ_dx_11);
	  tree3->SetBranchAddress("mip_vec_dQ_dx_12",&mip_vec_dQ_dx_12);
	  tree3->SetBranchAddress("mip_vec_dQ_dx_13",&mip_vec_dQ_dx_13);
	  tree3->SetBranchAddress("mip_vec_dQ_dx_14",&mip_vec_dQ_dx_14);
	  tree3->SetBranchAddress("mip_vec_dQ_dx_15",&mip_vec_dQ_dx_15);
	  tree3->SetBranchAddress("mip_vec_dQ_dx_16",&mip_vec_dQ_dx_16);
	  tree3->SetBranchAddress("mip_vec_dQ_dx_17",&mip_vec_dQ_dx_17);
	  tree3->SetBranchAddress("mip_vec_dQ_dx_18",&mip_vec_dQ_dx_18);
	  tree3->SetBranchAddress("mip_vec_dQ_dx_19",&mip_vec_dQ_dx_19);


	  float pio_filled;
	  float pio_flag;
	  float pio_mip_id;
	  float pio_flag_pio;
	  float pio_1_flag;
	  float pio_1_mass;
	  float pio_1_pio_type;
	  float pio_1_energy_1;
	  float pio_1_energy_2;
	  float pio_1_dis_1;
	  float pio_1_dis_2;
	  std::vector<float> *pio_2_v_flag = new std::vector<float>;
	  std::vector<float> *pio_2_v_dis2 = new std::vector<float>;
	  std::vector<float> *pio_2_v_angle2 = new std::vector<float>;
	  std::vector<float> *pio_2_v_acc_length = new std::vector<float>;
	  tree3->SetBranchAddress("pio_filled",&pio_filled);
	  tree3->SetBranchAddress("pio_flag",&pio_flag);
	  tree3->SetBranchAddress("pio_mip_id",&pio_mip_id);
	  tree3->SetBranchAddress("pio_flag_pio",&pio_flag_pio);
	  tree3->SetBranchAddress("pio_1_flag",&pio_1_flag);
	  tree3->SetBranchAddress("pio_1_mass",&pio_1_mass);
	  tree3->SetBranchAddress("pio_1_pio_type",&pio_1_pio_type);
	  tree3->SetBranchAddress("pio_1_energy_1",&pio_1_energy_1);
	  tree3->SetBranchAddress("pio_1_energy_2",&pio_1_energy_2);
	  tree3->SetBranchAddress("pio_1_dis_1",&pio_1_dis_1);
	  tree3->SetBranchAddress("pio_1_dis_2",&pio_1_dis_2);
	  tree3->SetBranchAddress("pio_2_v_flag",&pio_2_v_flag);
	  tree3->SetBranchAddress("pio_2_v_dis2",&pio_2_v_dis2);
	  tree3->SetBranchAddress("pio_2_v_angle2",&pio_2_v_angle2);
	  tree3->SetBranchAddress("pio_2_v_acc_length",&pio_2_v_acc_length);


	  float sig_flag;
	  std::vector<float> *sig_1_v_flag= new std::vector<float>;
	  std::vector<float> *sig_1_v_angle = new std::vector<float>;
	  std::vector<float> *sig_1_v_flag_single_shower= new std::vector<float>;
	  std::vector<float> *sig_1_v_energy= new std::vector<float>;
	  std::vector<float> *sig_1_v_energy_1= new std::vector<float>;
	  std::vector<float> *sig_2_v_flag= new std::vector<float>;
	  std::vector<float> *sig_2_v_energy= new std::vector<float>;
	  std::vector<float> *sig_2_v_shower_angle= new std::vector<float>;
	  std::vector<float> *sig_2_v_flag_single_shower= new std::vector<float>;
	  std::vector<float> *sig_2_v_medium_dQ_dx= new std::vector<float>;
	  std::vector<float> *sig_2_v_start_dQ_dx= new std::vector<float>;
	  tree3->SetBranchAddress("sig_flag",&sig_flag);
	  tree3->SetBranchAddress("sig_1_v_flag",&sig_1_v_flag);
	  tree3->SetBranchAddress("sig_1_v_angle",&sig_1_v_angle);
	  tree3->SetBranchAddress("sig_1_v_flag_single_shower",&sig_1_v_flag_single_shower);
	  tree3->SetBranchAddress("sig_1_v_energy",&sig_1_v_energy);
	  tree3->SetBranchAddress("sig_1_v_energy_1",&sig_1_v_energy_1);
	  tree3->SetBranchAddress("sig_2_v_flag",&sig_2_v_flag);
	  tree3->SetBranchAddress("sig_2_v_energy",&sig_2_v_energy);
	  tree3->SetBranchAddress("sig_2_v_shower_angle",&sig_2_v_shower_angle);
	  tree3->SetBranchAddress("sig_2_v_flag_single_shower",&sig_2_v_flag_single_shower);
	  tree3->SetBranchAddress("sig_2_v_medium_dQ_dx",&sig_2_v_medium_dQ_dx);
	  tree3->SetBranchAddress("sig_2_v_start_dQ_dx",&sig_2_v_start_dQ_dx);


	  float mgo_flag;
	  float mgo_energy;
	  float mgo_max_energy;
	  float mgo_total_energy;
	  float mgo_n_showers;
	  float mgo_max_energy_1;
	  float mgo_max_energy_2;
	  float mgo_total_other_energy;
	  float mgo_n_total_showers;
	  float mgo_total_other_energy_1;
	  tree3->SetBranchAddress("mgo_flag",&mgo_flag);
	  tree3->SetBranchAddress("mgo_energy",&mgo_energy);
	  tree3->SetBranchAddress("mgo_max_energy",&mgo_max_energy);
	  tree3->SetBranchAddress("mgo_total_energy",&mgo_total_energy);
	  tree3->SetBranchAddress("mgo_n_showers",&mgo_n_showers);
	  tree3->SetBranchAddress("mgo_max_energy_1",&mgo_max_energy_1);
	  tree3->SetBranchAddress("mgo_max_energy_2",&mgo_max_energy_2);
	  tree3->SetBranchAddress("mgo_total_other_energy",&mgo_total_other_energy);
	  tree3->SetBranchAddress("mgo_n_total_showers",&mgo_n_total_showers);
	  tree3->SetBranchAddress("mgo_total_other_energy_1",&mgo_total_other_energy_1);


	  float mgt_flag;
	  float mgt_flag_single_shower;
	  float mgt_max_energy;
	  float mgt_energy;
	  float mgt_total_other_energy;
	  float mgt_max_energy_1;
	  float mgt_e_indirect_max_energy;
	  float mgt_e_direct_max_energy;
	  float mgt_n_direct_showers;
	  float mgt_e_direct_total_energy;
	  float mgt_flag_indirect_max_pio;
	  float mgt_e_indirect_total_energy;
	  tree3->SetBranchAddress("mgt_flag",&mgt_flag);
	  tree3->SetBranchAddress("mgt_flag_single_shower",&mgt_flag_single_shower);
	  tree3->SetBranchAddress("mgt_max_energy",&mgt_max_energy);
	  tree3->SetBranchAddress("mgt_energy",&mgt_energy);
	  tree3->SetBranchAddress("mgt_total_other_energy",&mgt_total_other_energy);
	  tree3->SetBranchAddress("mgt_max_energy_1",&mgt_max_energy_1);
	  tree3->SetBranchAddress("mgt_e_indirect_max_energy",&mgt_e_indirect_max_energy);
	  tree3->SetBranchAddress("mgt_e_direct_max_energy",&mgt_e_direct_max_energy);
	  tree3->SetBranchAddress("mgt_n_direct_showers",&mgt_n_direct_showers);
	  tree3->SetBranchAddress("mgt_e_direct_total_energy",&mgt_e_direct_total_energy);
	  tree3->SetBranchAddress("mgt_flag_indirect_max_pio",&mgt_flag_indirect_max_pio);
	  tree3->SetBranchAddress("mgt_e_indirect_total_energy",&mgt_e_indirect_total_energy);


	  float stw_flag;
	  float stw_1_flag;
	  float stw_1_energy;
	  float stw_1_dis;
	  float stw_1_dQ_dx;
	  float stw_1_flag_single_shower;
	  float stw_1_n_pi0;
	  float stw_1_num_valid_tracks;
	  std::vector<float> *stw_2_v_flag = new std::vector<float>;
	  std::vector<float> *stw_2_v_medium_dQ_dx = new std::vector<float>;
	  std::vector<float> *stw_2_v_energy = new std::vector<float>;
	  std::vector<float> *stw_2_v_angle = new std::vector<float>;
	  std::vector<float> *stw_2_v_dir_length = new std::vector<float>;
	  std::vector<float> *stw_2_v_max_dQ_dx = new std::vector<float>;
	  std::vector<float> *stw_3_v_flag = new std::vector<float>;
	  std::vector<float> *stw_3_v_angle = new std::vector<float>;
	  std::vector<float> *stw_3_v_dir_length = new std::vector<float>;
	  std::vector<float> *stw_3_v_energy = new std::vector<float>;
	  std::vector<float> *stw_3_v_medium_dQ_dx = new std::vector<float>;
	  std::vector<float> *stw_4_v_flag = new std::vector<float>;
	  std::vector<float> *stw_4_v_angle = new std::vector<float>;
	  std::vector<float> *stw_4_v_dis = new std::vector<float>;
	  std::vector<float> *stw_4_v_energy = new std::vector<float>;
	  tree3->SetBranchAddress("stw_flag", &stw_flag);
	  tree3->SetBranchAddress("stw_1_flag",&stw_1_flag);
	  tree3->SetBranchAddress("stw_1_energy",&stw_1_energy);
	  tree3->SetBranchAddress("stw_1_dis",&stw_1_dis);
	  tree3->SetBranchAddress("stw_1_dQ_dx",&stw_1_dQ_dx);
	  tree3->SetBranchAddress("stw_1_flag_single_shower",&stw_1_flag_single_shower);
	  tree3->SetBranchAddress("stw_1_n_pi0",&stw_1_n_pi0);
	  tree3->SetBranchAddress("stw_1_num_valid_tracks",&stw_1_num_valid_tracks);
	  tree3->SetBranchAddress("stw_2_v_flag", &stw_2_v_flag);
	  tree3->SetBranchAddress("stw_2_v_medium_dQ_dx", &stw_2_v_medium_dQ_dx);
	  tree3->SetBranchAddress("stw_2_v_energy", &stw_2_v_energy);
	  tree3->SetBranchAddress("stw_2_v_angle", &stw_2_v_angle);
	  tree3->SetBranchAddress("stw_2_v_dir_length", &stw_2_v_dir_length);
	  tree3->SetBranchAddress("stw_2_v_max_dQ_dx", &stw_2_v_max_dQ_dx);
	  tree3->SetBranchAddress("stw_3_v_flag",&stw_3_v_flag);
	  tree3->SetBranchAddress("stw_3_v_angle",&stw_3_v_angle);
	  tree3->SetBranchAddress("stw_3_v_dir_length",&stw_3_v_dir_length);
	  tree3->SetBranchAddress("stw_3_v_energy",&stw_3_v_energy);
	  tree3->SetBranchAddress("stw_3_v_medium_dQ_dx",&stw_3_v_medium_dQ_dx);
	  tree3->SetBranchAddress("stw_4_v_flag",&stw_4_v_flag);
	  tree3->SetBranchAddress("stw_4_v_angle",&stw_4_v_angle);
	  tree3->SetBranchAddress("stw_4_v_dis",&stw_4_v_dis);
	  tree3->SetBranchAddress("stw_4_v_energy",&stw_4_v_energy);


	  float spt_flag;
	  float spt_flag_single_shower;
	  float spt_energy;
	  float spt_shower_main_length;
	  float spt_shower_total_length;
	  float spt_angle_beam;
	  float spt_angle_vertical;
	  float spt_max_dQ_dx;
	  float spt_angle_beam_1;
	  float spt_angle_drift;
	  float spt_angle_drift_1;
	  float spt_num_valid_tracks;
	  float spt_n_vtx_segs;
	  float spt_max_length;
	  tree3->SetBranchAddress("spt_flag", &spt_flag);
	  tree3->SetBranchAddress("spt_flag_single_shower", &spt_flag_single_shower);
	  tree3->SetBranchAddress("spt_energy", &spt_energy);
	  tree3->SetBranchAddress("spt_shower_main_length", &spt_shower_main_length);
	  tree3->SetBranchAddress("spt_shower_total_length", &spt_shower_total_length);
	  tree3->SetBranchAddress("spt_angle_beam", &spt_angle_beam);
	  tree3->SetBranchAddress("spt_angle_vertical", &spt_angle_vertical);
	  tree3->SetBranchAddress("spt_max_dQ_dx", &spt_max_dQ_dx);
	  tree3->SetBranchAddress("spt_angle_beam_1", &spt_angle_beam_1);
	  tree3->SetBranchAddress("spt_angle_drift", &spt_angle_drift);
	  tree3->SetBranchAddress("spt_angle_drift_1", &spt_angle_drift_1);
	  tree3->SetBranchAddress("spt_num_valid_tracks", &spt_num_valid_tracks);
	  tree3->SetBranchAddress("spt_n_vtx_segs", &spt_n_vtx_segs);
	  tree3->SetBranchAddress("spt_max_length", &spt_max_length);


	  float stem_len_flag;
	  float stem_len_energy;
	  float stem_len_length;
	  float stem_len_flag_avoid_muon_check;
	  float stem_len_num_daughters;
	  float stem_len_daughter_length;
	  tree3->SetBranchAddress("stem_len_flag", &stem_len_flag);
	  tree3->SetBranchAddress("stem_len_energy", &stem_len_energy);
	  tree3->SetBranchAddress("stem_len_length", &stem_len_length);
	  tree3->SetBranchAddress("stem_len_flag_avoid_muon_check", &stem_len_flag_avoid_muon_check);
	  tree3->SetBranchAddress("stem_len_num_daughters", &stem_len_num_daughters);
	  tree3->SetBranchAddress("stem_len_daughter_length", &stem_len_daughter_length);


	  float lem_flag;
	  float lem_shower_total_length;
	  float lem_shower_main_length;
	  float lem_n_3seg;
	  float lem_e_charge;
	  float lem_e_dQdx;
	  float lem_shower_num_segs;
	  float lem_shower_num_main_segs;
	  tree3->SetBranchAddress("lem_flag",&lem_flag);
	  tree3->SetBranchAddress("lem_shower_total_length",&lem_shower_total_length);
	  tree3->SetBranchAddress("lem_shower_main_length",&lem_shower_main_length);
	  tree3->SetBranchAddress("lem_n_3seg",&lem_n_3seg);
	  tree3->SetBranchAddress("lem_e_charge",&lem_e_charge);
	  tree3->SetBranchAddress("lem_e_dQdx",&lem_e_dQdx);
	  tree3->SetBranchAddress("lem_shower_num_segs",&lem_shower_num_segs);
	  tree3->SetBranchAddress("lem_shower_num_main_segs",&lem_shower_num_main_segs);


	  float brm_flag;
	  float brm_n_mu_segs;
	  float brm_Ep;
	  float brm_energy;
	  float brm_acc_length;
	  float brm_shower_total_length;
	  float brm_connected_length;
	  float brm_n_size;
	  float brm_acc_direct_length;
	  float brm_n_shower_main_segs;
	  float brm_n_mu_main;
	  tree3->SetBranchAddress("brm_flag",&brm_flag);
	  tree3->SetBranchAddress("brm_n_mu_segs",&brm_n_mu_segs);
	  tree3->SetBranchAddress("brm_Ep",&brm_Ep);
	  tree3->SetBranchAddress("brm_energy",&brm_energy);
	  tree3->SetBranchAddress("brm_acc_length",&brm_acc_length);
	  tree3->SetBranchAddress("brm_shower_total_length",&brm_shower_total_length);
	  tree3->SetBranchAddress("brm_connected_length",&brm_connected_length);
	  tree3->SetBranchAddress("brm_n_size",&brm_n_size);
	  tree3->SetBranchAddress("brm_acc_direct_length",&brm_acc_direct_length);
	  tree3->SetBranchAddress("brm_n_shower_main_segs",&brm_n_shower_main_segs);
	  tree3->SetBranchAddress("brm_n_mu_main",&brm_n_mu_main);


	  float cme_flag;
	  float cme_mu_energy;
	  float cme_energy;
	  float cme_mu_length;
	  float cme_length;
	  float cme_angle_beam;
	  tree3->SetBranchAddress("cme_flag",&cme_flag);
	  tree3->SetBranchAddress("cme_mu_energy",&cme_mu_energy);
	  tree3->SetBranchAddress("cme_energy",&cme_energy);
	  tree3->SetBranchAddress("cme_mu_length",&cme_mu_length);
	  tree3->SetBranchAddress("cme_length",&cme_length);
	  tree3->SetBranchAddress("cme_angle_beam",&cme_angle_beam);


	  float anc_flag;
	  float anc_energy;
	  float anc_angle;
	  float anc_max_angle;
	  float anc_max_length;
	  float anc_acc_forward_length;
	  float anc_acc_backward_length;
	  float anc_acc_forward_length1;
	  float anc_shower_main_length;
	  float anc_shower_total_length;
	  float anc_flag_main_outside;
	  tree3->SetBranchAddress("anc_flag",&anc_flag);
	  tree3->SetBranchAddress("anc_energy",&anc_energy);
	  tree3->SetBranchAddress("anc_angle",&anc_angle);
	  tree3->SetBranchAddress("anc_max_angle",&anc_max_angle);
	  tree3->SetBranchAddress("anc_max_length",&anc_max_length);
	  tree3->SetBranchAddress("anc_acc_forward_length",&anc_acc_forward_length);
	  tree3->SetBranchAddress("anc_acc_backward_length",&anc_acc_backward_length);
	  tree3->SetBranchAddress("anc_acc_forward_length1",&anc_acc_forward_length1);
	  tree3->SetBranchAddress("anc_shower_main_length",&anc_shower_main_length);
	  tree3->SetBranchAddress("anc_shower_total_length",&anc_shower_total_length);
	  tree3->SetBranchAddress("anc_flag_main_outside",&anc_flag_main_outside);


	  float stem_dir_filled;
	  float stem_dir_flag;
	  float stem_dir_flag_single_shower;
	  float stem_dir_angle;
	  float stem_dir_energy;
	  float stem_dir_angle1;
	  float stem_dir_angle2;
	  float stem_dir_angle3;
	  float stem_dir_ratio;
	  tree3->SetBranchAddress("stem_dir_filled",&stem_dir_filled);
	  tree3->SetBranchAddress("stem_dir_flag",&stem_dir_flag);
	  tree3->SetBranchAddress("stem_dir_flag_single_shower",&stem_dir_flag_single_shower);
	  tree3->SetBranchAddress("stem_dir_angle",&stem_dir_angle);
	  tree3->SetBranchAddress("stem_dir_energy",&stem_dir_energy);
	  tree3->SetBranchAddress("stem_dir_angle1",&stem_dir_angle1);
	  tree3->SetBranchAddress("stem_dir_angle2",&stem_dir_angle2);
	  tree3->SetBranchAddress("stem_dir_angle3",&stem_dir_angle3);
	  tree3->SetBranchAddress("stem_dir_ratio",&stem_dir_ratio);


	  float vis_flag;
	  float vis_1_filled;
	  float vis_1_flag;
	  float vis_1_n_vtx_segs;
	  float vis_1_energy;
	  float vis_1_num_good_tracks;
	  float vis_1_max_angle;
	  float vis_1_max_shower_angle;
	  float vis_1_tmp_length1;
	  float vis_1_tmp_length2;
	  float vis_1_particle_type;
	  float vis_2_filled;
	  float vis_2_flag;
	  float vis_2_n_vtx_segs;
	  float vis_2_min_angle;
	  float vis_2_min_weak_track;
	  float vis_2_angle_beam;
	  float vis_2_min_angle1;
	  float vis_2_iso_angle1;
	  float vis_2_min_medium_dQ_dx;
	  float vis_2_min_length;
	  float vis_2_sg_length;
	  float vis_2_max_angle;
	  float vis_2_max_weak_track;
	  tree3->SetBranchAddress("vis_flag",&vis_flag);
	  tree3->SetBranchAddress("vis_1_filled",&vis_1_filled);
	  tree3->SetBranchAddress("vis_1_flag",&vis_1_flag);
	  tree3->SetBranchAddress("vis_1_n_vtx_segs",&vis_1_n_vtx_segs);
	  tree3->SetBranchAddress("vis_1_energy",&vis_1_energy);
	  tree3->SetBranchAddress("vis_1_num_good_tracks",&vis_1_num_good_tracks);
	  tree3->SetBranchAddress("vis_1_max_angle",&vis_1_max_angle);
	  tree3->SetBranchAddress("vis_1_max_shower_angle",&vis_1_max_shower_angle);
	  tree3->SetBranchAddress("vis_1_tmp_length1",&vis_1_tmp_length1);
	  tree3->SetBranchAddress("vis_1_tmp_length2",&vis_1_tmp_length2);
	  tree3->SetBranchAddress("vis_1_particle_type",&vis_1_particle_type);
	  tree3->SetBranchAddress("vis_2_filled",&vis_2_filled);
	  tree3->SetBranchAddress("vis_2_flag",&vis_2_flag);
	  tree3->SetBranchAddress("vis_2_n_vtx_segs",&vis_2_n_vtx_segs);
	  tree3->SetBranchAddress("vis_2_min_angle",&vis_2_min_angle);
	  tree3->SetBranchAddress("vis_2_min_weak_track",&vis_2_min_weak_track);
	  tree3->SetBranchAddress("vis_2_angle_beam",&vis_2_angle_beam);
	  tree3->SetBranchAddress("vis_2_min_angle1",&vis_2_min_angle1);
	  tree3->SetBranchAddress("vis_2_iso_angle1",&vis_2_iso_angle1);
	  tree3->SetBranchAddress("vis_2_min_medium_dQ_dx",&vis_2_min_medium_dQ_dx);
	  tree3->SetBranchAddress("vis_2_min_length",&vis_2_min_length);
	  tree3->SetBranchAddress("vis_2_sg_length",&vis_2_sg_length);
	  tree3->SetBranchAddress("vis_2_max_angle",&vis_2_max_angle);
	  tree3->SetBranchAddress("vis_2_max_weak_track",&vis_2_max_weak_track);


	  float br_filled;
	  float br1_flag;
	  float br1_1_flag;
	  float br1_1_shower_type;
	  float br1_1_vtx_n_segs;
	  float br1_1_energy;
	  float br1_1_n_segs;
	  float br1_1_flag_sg_topology;
	  float br1_1_flag_sg_trajectory;
	  float br1_1_sg_length;
	  float br1_2_flag;
	  float br1_2_energy;
	  float br1_2_n_connected;
	  float br1_2_max_length;
	  float br1_2_n_connected_1;
	  float br1_2_vtx_n_segs;
	  float br1_2_n_shower_segs;
	  float br1_2_max_length_ratio;
	  float br1_2_shower_length;
	  float br1_3_flag;
	  float br1_3_energy;
	  float br1_3_n_connected_p;
	  float br1_3_max_length_p;
	  float br1_3_n_shower_segs;
	  float br1_3_flag_sg_topology;
	  float br1_3_flag_sg_trajectory;
	  float br1_3_n_shower_main_segs;
	  float br1_3_sg_length;
	  tree3->SetBranchAddress("br_filled",&br_filled);
	  tree3->SetBranchAddress("br1_flag",&br1_flag);
	  tree3->SetBranchAddress("br1_1_flag",&br1_1_flag);
	  tree3->SetBranchAddress("br1_1_shower_type",&br1_1_shower_type);
	  tree3->SetBranchAddress("br1_1_vtx_n_segs",&br1_1_vtx_n_segs);
	  tree3->SetBranchAddress("br1_1_energy",&br1_1_energy);
	  tree3->SetBranchAddress("br1_1_n_segs",&br1_1_n_segs);
	  tree3->SetBranchAddress("br1_1_flag_sg_topology",&br1_1_flag_sg_topology);
	  tree3->SetBranchAddress("br1_1_flag_sg_trajectory",&br1_1_flag_sg_trajectory);
	  tree3->SetBranchAddress("br1_1_sg_length",&br1_1_sg_length);
	  tree3->SetBranchAddress("br1_2_flag",&br1_2_flag);
	  tree3->SetBranchAddress("br1_2_energy",&br1_2_energy);
	  tree3->SetBranchAddress("br1_2_n_connected",&br1_2_n_connected);
	  tree3->SetBranchAddress("br1_2_max_length",&br1_2_max_length);
	  tree3->SetBranchAddress("br1_2_n_connected_1",&br1_2_n_connected_1);
	  tree3->SetBranchAddress("br1_2_vtx_n_segs",&br1_2_vtx_n_segs);
	  tree3->SetBranchAddress("br1_2_n_shower_segs",&br1_2_n_shower_segs);
	  tree3->SetBranchAddress("br1_2_max_length_ratio",&br1_2_max_length_ratio);
	  tree3->SetBranchAddress("br1_2_shower_length",&br1_2_shower_length);
	  tree3->SetBranchAddress("br1_3_flag",&br1_3_flag);
	  tree3->SetBranchAddress("br1_3_energy",&br1_3_energy);
	  tree3->SetBranchAddress("br1_3_n_connected_p",&br1_3_n_connected_p);
	  tree3->SetBranchAddress("br1_3_max_length_p",&br1_3_max_length_p);
	  tree3->SetBranchAddress("br1_3_n_shower_segs",&br1_3_n_shower_segs);
	  tree3->SetBranchAddress("br1_3_flag_sg_topology",&br1_3_flag_sg_topology);
	  tree3->SetBranchAddress("br1_3_flag_sg_trajectory",&br1_3_flag_sg_trajectory);
	  tree3->SetBranchAddress("br1_3_n_shower_main_segs",&br1_3_n_shower_main_segs);
	  tree3->SetBranchAddress("br1_3_sg_length",&br1_3_sg_length);


	  float br2_flag;
	  float br2_flag_single_shower;
	  float br2_num_valid_tracks;
	  float br2_energy;
	  float br2_angle1;
	  float br2_angle2;
	  float br2_angle;
	  float br2_angle3;
	  float br2_n_shower_main_segs;
	  float br2_max_angle;
	  float br2_sg_length;
	  float br2_flag_sg_trajectory;
	  tree3->SetBranchAddress("br2_flag",&br2_flag);
	  tree3->SetBranchAddress("br2_flag_single_shower",&br2_flag_single_shower);
	  tree3->SetBranchAddress("br2_num_valid_tracks",&br2_num_valid_tracks);
	  tree3->SetBranchAddress("br2_energy",&br2_energy);
	  tree3->SetBranchAddress("br2_angle1",&br2_angle1);
	  tree3->SetBranchAddress("br2_angle2",&br2_angle2);
	  tree3->SetBranchAddress("br2_angle",&br2_angle);
	  tree3->SetBranchAddress("br2_angle3",&br2_angle3);
	  tree3->SetBranchAddress("br2_n_shower_main_segs",&br2_n_shower_main_segs);
	  tree3->SetBranchAddress("br2_max_angle",&br2_max_angle);
	  tree3->SetBranchAddress("br2_sg_length",&br2_sg_length);
	  tree3->SetBranchAddress("br2_flag_sg_trajectory",&br2_flag_sg_trajectory);


	  float br3_flag;
	  float br3_1_flag;
	  float br3_1_energy;
	  float br3_1_n_shower_segments;
	  float br3_1_sg_flag_trajectory;
	  float br3_1_sg_direct_length;
	  float br3_1_sg_length;
	  float br3_1_total_main_length;
	  float br3_1_total_length;
	  float br3_1_iso_angle;
	  float br3_1_sg_flag_topology;
	  float br3_2_flag;
	  float br3_2_n_ele;
	  float br3_2_n_other;
	  float br3_2_energy;
	  float br3_2_total_main_length;
	  float br3_2_total_length;
	  float br3_2_other_fid;
	  std::vector<float> *br3_3_v_flag = new std::vector<float>;
	  std::vector<float> *br3_3_v_energy = new std::vector<float>;
	  std::vector<float> *br3_3_v_angle = new std::vector<float>;
	  std::vector<float> *br3_3_v_dir_length = new std::vector<float>;
	  std::vector<float> *br3_3_v_length = new std::vector<float>;
	  float br3_4_flag;
	  float br3_4_acc_length;
	  float br3_4_total_length;
	  float br3_4_energy;
	  std::vector<float> *br3_5_v_flag = new std::vector<float>;
	  std::vector<float> *br3_5_v_dir_length = new std::vector<float>;
	  std::vector<float> *br3_5_v_total_length = new std::vector<float>;
	  std::vector<float> *br3_5_v_flag_avoid_muon_check = new std::vector<float>;
	  std::vector<float> *br3_5_v_n_seg = new std::vector<float>;
	  std::vector<float> *br3_5_v_angle = new std::vector<float>;
	  std::vector<float> *br3_5_v_sg_length = new std::vector<float>;
	  std::vector<float> *br3_5_v_energy = new std::vector<float>;
	  std::vector<float> *br3_5_v_n_main_segs = new std::vector<float>;
	  std::vector<float> *br3_5_v_n_segs = new std::vector<float>;
	  std::vector<float> *br3_5_v_shower_main_length = new std::vector<float>;
	  std::vector<float> *br3_5_v_shower_total_length = new std::vector<float>;
	  std::vector<float> *br3_6_v_flag = new std::vector<float>;
	  std::vector<float> *br3_6_v_angle = new std::vector<float>;
	  std::vector<float> *br3_6_v_angle1 = new std::vector<float>;
	  std::vector<float> *br3_6_v_flag_shower_trajectory = new std::vector<float>;
	  std::vector<float> *br3_6_v_direct_length = new std::vector<float>;
	  std::vector<float> *br3_6_v_length = new std::vector<float>;
	  std::vector<float> *br3_6_v_n_other_vtx_segs = new std::vector<float>;
	  std::vector<float> *br3_6_v_energy = new std::vector<float>;
	  float br3_7_flag;
	  float br3_7_energy;
	  float br3_7_min_angle;
	  float br3_7_sg_length;
	  float br3_7_shower_main_length;
	  float br3_8_flag;
	  float br3_8_max_dQ_dx;
	  float br3_8_energy;
	  float br3_8_n_main_segs;
	  float br3_8_shower_main_length;
	  float br3_8_shower_length;
	  tree3->SetBranchAddress("br3_flag",&br3_flag);
	  tree3->SetBranchAddress("br3_1_flag",&br3_1_flag);
	  tree3->SetBranchAddress("br3_1_energy",&br3_1_energy);
	  tree3->SetBranchAddress("br3_1_n_shower_segments",&br3_1_n_shower_segments);
	  tree3->SetBranchAddress("br3_1_sg_flag_trajectory",&br3_1_sg_flag_trajectory);
	  tree3->SetBranchAddress("br3_1_sg_direct_length",&br3_1_sg_direct_length);
	  tree3->SetBranchAddress("br3_1_sg_length",&br3_1_sg_length);
	  tree3->SetBranchAddress("br3_1_total_main_length",&br3_1_total_main_length);
	  tree3->SetBranchAddress("br3_1_total_length",&br3_1_total_length);
	  tree3->SetBranchAddress("br3_1_iso_angle",&br3_1_iso_angle);
	  tree3->SetBranchAddress("br3_1_sg_flag_topology",&br3_1_sg_flag_topology);
	  tree3->SetBranchAddress("br3_2_flag",&br3_2_flag);
	  tree3->SetBranchAddress("br3_2_n_ele",&br3_2_n_ele);
	  tree3->SetBranchAddress("br3_2_n_other",&br3_2_n_other);
	  tree3->SetBranchAddress("br3_2_energy",&br3_2_energy);
	  tree3->SetBranchAddress("br3_2_total_main_length",&br3_2_total_main_length);
	  tree3->SetBranchAddress("br3_2_total_length",&br3_2_total_length);
	  tree3->SetBranchAddress("br3_2_other_fid",&br3_2_other_fid);
	  tree3->SetBranchAddress("br3_3_v_flag",&br3_3_v_flag);
	  tree3->SetBranchAddress("br3_3_v_energy",&br3_3_v_energy);
	  tree3->SetBranchAddress("br3_3_v_angle",&br3_3_v_angle);
	  tree3->SetBranchAddress("br3_3_v_dir_length",&br3_3_v_dir_length);
	  tree3->SetBranchAddress("br3_3_v_length",&br3_3_v_length);
	  tree3->SetBranchAddress("br3_4_flag", &br3_4_flag);
	  tree3->SetBranchAddress("br3_4_acc_length", &br3_4_acc_length);
	  tree3->SetBranchAddress("br3_4_total_length", &br3_4_total_length);
	  tree3->SetBranchAddress("br3_4_energy", &br3_4_energy);
	  tree3->SetBranchAddress("br3_5_v_flag", &br3_5_v_flag);
	  tree3->SetBranchAddress("br3_5_v_dir_length", &br3_5_v_dir_length);
	  tree3->SetBranchAddress("br3_5_v_total_length", &br3_5_v_total_length);
	  tree3->SetBranchAddress("br3_5_v_flag_avoid_muon_check", &br3_5_v_flag_avoid_muon_check);
	  tree3->SetBranchAddress("br3_5_v_n_seg", &br3_5_v_n_seg);
	  tree3->SetBranchAddress("br3_5_v_angle", &br3_5_v_angle);
	  tree3->SetBranchAddress("br3_5_v_sg_length", &br3_5_v_sg_length);
	  tree3->SetBranchAddress("br3_5_v_energy", &br3_5_v_energy);
	  tree3->SetBranchAddress("br3_5_v_n_main_segs", &br3_5_v_n_main_segs);
	  tree3->SetBranchAddress("br3_5_v_n_segs", &br3_5_v_n_segs);
	  tree3->SetBranchAddress("br3_5_v_shower_main_length", &br3_5_v_shower_main_length);
	  tree3->SetBranchAddress("br3_5_v_shower_total_length", &br3_5_v_shower_total_length);
	  tree3->SetBranchAddress("br3_6_v_flag",&br3_6_v_flag);
	  tree3->SetBranchAddress("br3_6_v_angle",&br3_6_v_angle);
	  tree3->SetBranchAddress("br3_6_v_angle1",&br3_6_v_angle1);
	  tree3->SetBranchAddress("br3_6_v_flag_shower_trajectory",&br3_6_v_flag_shower_trajectory);
	  tree3->SetBranchAddress("br3_6_v_direct_length",&br3_6_v_direct_length);
	  tree3->SetBranchAddress("br3_6_v_length",&br3_6_v_length);
	  tree3->SetBranchAddress("br3_6_v_n_other_vtx_segs",&br3_6_v_n_other_vtx_segs);
	  tree3->SetBranchAddress("br3_6_v_energy",&br3_6_v_energy);
	  tree3->SetBranchAddress("br3_7_flag",&br3_7_flag);
	  tree3->SetBranchAddress("br3_7_energy",&br3_7_energy);
	  tree3->SetBranchAddress("br3_7_min_angle",&br3_7_min_angle);
	  tree3->SetBranchAddress("br3_7_sg_length",&br3_7_sg_length);
	  tree3->SetBranchAddress("br3_7_main_length",&br3_7_shower_main_length);
	  tree3->SetBranchAddress("br3_8_flag",&br3_8_flag);
	  tree3->SetBranchAddress("br3_8_max_dQ_dx",&br3_8_max_dQ_dx);
	  tree3->SetBranchAddress("br3_8_energy",&br3_8_energy);
	  tree3->SetBranchAddress("br3_8_n_main_segs",&br3_8_n_main_segs);
	  tree3->SetBranchAddress("br3_8_shower_main_length",&br3_8_shower_main_length);
	  tree3->SetBranchAddress("br3_8_shower_length",&br3_8_shower_length);


	  float br4_flag;
	  float br4_1_flag;
	  float br4_1_shower_main_length;
	  float br4_1_shower_total_length;
	  float br4_1_min_dis;
	  float br4_1_energy;
	  float br4_1_flag_avoid_muon_check;
	  float br4_1_n_vtx_segs;
	  float br4_1_n_main_segs;
	  float br4_2_flag;
	  float br4_2_ratio_45;
	  float br4_2_ratio_35;
	  float br4_2_ratio_25;
	  float br4_2_ratio_15;
	  float br4_2_energy;
	  float br4_2_ratio1_45;
	  float br4_2_ratio1_35;
	  float br4_2_ratio1_25;
	  float br4_2_ratio1_15;
	  float br4_2_iso_angle;
	  float br4_2_iso_angle1;
	  float br4_2_angle;
	  tree3->SetBranchAddress("br4_flag", &br4_flag);
	  tree3->SetBranchAddress("br4_1_flag", &br4_1_flag);
	  tree3->SetBranchAddress("br4_1_shower_main_length", &br4_1_shower_main_length);
	  tree3->SetBranchAddress("br4_1_shower_total_length", &br4_1_shower_total_length);
	  tree3->SetBranchAddress("br4_1_min_dis", &br4_1_min_dis);
	  tree3->SetBranchAddress("br4_1_energy", &br4_1_energy);
	  tree3->SetBranchAddress("br4_1_flag_avoid_muon_check", &br4_1_flag_avoid_muon_check);
	  tree3->SetBranchAddress("br4_1_n_vtx_segs", &br4_1_n_vtx_segs);
	  tree3->SetBranchAddress("br4_1_n_main_segs", &br4_1_n_main_segs);
	  tree3->SetBranchAddress("br4_2_flag", &br4_2_flag);
	  tree3->SetBranchAddress("br4_2_ratio_45", &br4_2_ratio_45);
	  tree3->SetBranchAddress("br4_2_ratio_35", &br4_2_ratio_35);
	  tree3->SetBranchAddress("br4_2_ratio_25", &br4_2_ratio_25);
	  tree3->SetBranchAddress("br4_2_ratio_15", &br4_2_ratio_15);
	  tree3->SetBranchAddress("br4_2_energy",   &br4_2_energy);
	  tree3->SetBranchAddress("br4_2_ratio1_45", &br4_2_ratio1_45);
	  tree3->SetBranchAddress("br4_2_ratio1_35", &br4_2_ratio1_35);
	  tree3->SetBranchAddress("br4_2_ratio1_25", &br4_2_ratio1_25);
	  tree3->SetBranchAddress("br4_2_ratio1_15", &br4_2_ratio1_15);
	  tree3->SetBranchAddress("br4_2_iso_angle", &br4_2_iso_angle);
	  tree3->SetBranchAddress("br4_2_iso_angle1", &br4_2_iso_angle1);
	  tree3->SetBranchAddress("br4_2_angle", &br4_2_angle);


	  float tro_flag;
	  std::vector<float> *tro_1_v_flag= new std::vector<float>;
	  std::vector<float> *tro_1_v_particle_type= new std::vector<float>;
	  std::vector<float> *tro_1_v_flag_dir_weak= new std::vector<float>;
	  std::vector<float> *tro_1_v_min_dis = new std::vector<float>;
	  std::vector<float> *tro_1_v_sg1_length = new std::vector<float>;
	  std::vector<float> *tro_1_v_shower_main_length = new std::vector<float>;
	  std::vector<float> *tro_1_v_max_n_vtx_segs= new std::vector<float>;
	  std::vector<float> *tro_1_v_tmp_length = new std::vector<float>;
	  std::vector<float> *tro_1_v_medium_dQ_dx = new std::vector<float>;
	  std::vector<float> *tro_1_v_dQ_dx_cut = new std::vector<float>;
	  std::vector<float> *tro_1_v_flag_shower_topology= new std::vector<float>;
	  std::vector<float> *tro_2_v_flag= new std::vector<float>;
	  std::vector<float> *tro_2_v_energy = new std::vector<float>;
	  std::vector<float> *tro_2_v_stem_length = new std::vector<float>;
	  std::vector<float> *tro_2_v_iso_angle = new std::vector<float>;
	  std::vector<float> *tro_2_v_max_length = new std::vector<float>;
	  std::vector<float> *tro_2_v_angle = new std::vector<float>;
	  float tro_3_flag;
	  float tro_3_stem_length;
	  float tro_3_n_muon_segs;
	  float tro_3_energy;
	  std::vector<float> *tro_4_v_flag= new std::vector<float>;
	  std::vector<float> *tro_4_v_dir2_mag = new std::vector<float>;
	  std::vector<float> *tro_4_v_angle = new std::vector<float>;
	  std::vector<float> *tro_4_v_angle1 = new std::vector<float>;
	  std::vector<float> *tro_4_v_angle2 = new std::vector<float>;
	  std::vector<float> *tro_4_v_length = new std::vector<float>;
	  std::vector<float> *tro_4_v_length1 = new std::vector<float>;
	  std::vector<float> *tro_4_v_medium_dQ_dx = new std::vector<float>;
	  std::vector<float> *tro_4_v_end_dQ_dx = new std::vector<float>;
	  std::vector<float> *tro_4_v_energy = new std::vector<float>;
	  std::vector<float> *tro_4_v_shower_main_length = new std::vector<float>;
	  std::vector<float> *tro_4_v_flag_shower_trajectory= new std::vector<float>;
	  std::vector<float> *tro_5_v_flag = new std::vector<float>;
	  std::vector<float> *tro_5_v_max_angle = new std::vector<float>;
	  std::vector<float> *tro_5_v_min_angle = new std::vector<float>;
	  std::vector<float> *tro_5_v_max_length = new std::vector<float>;
	  std::vector<float> *tro_5_v_iso_angle = new std::vector<float>;
	  std::vector<float> *tro_5_v_n_vtx_segs= new std::vector<float>;
	  std::vector<float> *tro_5_v_min_count= new std::vector<float>;
	  std::vector<float> *tro_5_v_max_count= new std::vector<float>;
	  std::vector<float> *tro_5_v_energy = new std::vector<float>;
	  tree3->SetBranchAddress("tro_flag",&tro_flag);
	  tree3->SetBranchAddress("tro_1_v_flag",&tro_1_v_flag);
	  tree3->SetBranchAddress("tro_1_v_particle_type",&tro_1_v_particle_type);
	  tree3->SetBranchAddress("tro_1_v_flag_dir_weak",&tro_1_v_flag_dir_weak);
	  tree3->SetBranchAddress("tro_1_v_min_dis",&tro_1_v_min_dis);
	  tree3->SetBranchAddress("tro_1_v_sg1_length",&tro_1_v_sg1_length);
	  tree3->SetBranchAddress("tro_1_v_shower_main_length",&tro_1_v_shower_main_length);
	  tree3->SetBranchAddress("tro_1_v_max_n_vtx_segs",&tro_1_v_max_n_vtx_segs);
	  tree3->SetBranchAddress("tro_1_v_tmp_length",&tro_1_v_tmp_length);
	  tree3->SetBranchAddress("tro_1_v_medium_dQ_dx",&tro_1_v_medium_dQ_dx);
	  tree3->SetBranchAddress("tro_1_v_dQ_dx_cut",&tro_1_v_dQ_dx_cut);
	  tree3->SetBranchAddress("tro_1_v_flag_shower_topology",&tro_1_v_flag_shower_topology);
	  tree3->SetBranchAddress("tro_2_v_flag",&tro_2_v_flag);
	  tree3->SetBranchAddress("tro_2_v_energy",&tro_2_v_energy);
	  tree3->SetBranchAddress("tro_2_v_stem_length",&tro_2_v_stem_length);
	  tree3->SetBranchAddress("tro_2_v_iso_angle",&tro_2_v_iso_angle);
	  tree3->SetBranchAddress("tro_2_v_max_length",&tro_2_v_max_length);
	  tree3->SetBranchAddress("tro_2_v_angle",&tro_2_v_angle);
	  tree3->SetBranchAddress("tro_3_flag",&tro_3_flag);
	  tree3->SetBranchAddress("tro_3_stem_length",&tro_3_stem_length);
	  tree3->SetBranchAddress("tro_3_n_muon_segs",&tro_3_n_muon_segs);
	  tree3->SetBranchAddress("tro_3_energy",&tro_3_energy);
	  tree3->SetBranchAddress("tro_4_v_flag",&tro_4_v_flag);
	  tree3->SetBranchAddress("tro_4_v_dir2_mag",&tro_4_v_dir2_mag);
	  tree3->SetBranchAddress("tro_4_v_angle",&tro_4_v_angle);
	  tree3->SetBranchAddress("tro_4_v_angle1",&tro_4_v_angle1);
	  tree3->SetBranchAddress("tro_4_v_angle2",&tro_4_v_angle2);
	  tree3->SetBranchAddress("tro_4_v_length",&tro_4_v_length);
	  tree3->SetBranchAddress("tro_4_v_length1",&tro_4_v_length1);
	  tree3->SetBranchAddress("tro_4_v_medium_dQ_dx",&tro_4_v_medium_dQ_dx);
	  tree3->SetBranchAddress("tro_4_v_end_dQ_dx",&tro_4_v_end_dQ_dx);
	  tree3->SetBranchAddress("tro_4_v_energy",&tro_4_v_energy);
	  tree3->SetBranchAddress("tro_4_v_shower_main_length",&tro_4_v_shower_main_length);
	  tree3->SetBranchAddress("tro_4_v_flag_shower_trajectory",&tro_4_v_flag_shower_trajectory);
	  tree3->SetBranchAddress("tro_5_v_flag",&tro_5_v_flag);
	  tree3->SetBranchAddress("tro_5_v_max_angle",&tro_5_v_max_angle);
	  tree3->SetBranchAddress("tro_5_v_min_angle",&tro_5_v_min_angle);
	  tree3->SetBranchAddress("tro_5_v_max_length",&tro_5_v_max_length);
	  tree3->SetBranchAddress("tro_5_v_iso_angle",&tro_5_v_iso_angle);
	  tree3->SetBranchAddress("tro_5_v_n_vtx_segs",&tro_5_v_n_vtx_segs);
	  tree3->SetBranchAddress("tro_5_v_min_count",&tro_5_v_min_count);
	  tree3->SetBranchAddress("tro_5_v_max_count",&tro_5_v_max_count);
	  tree3->SetBranchAddress("tro_5_v_energy",&tro_5_v_energy);


	  float hol_flag;
	  float hol_1_flag;
	  float hol_1_n_valid_tracks;
	  float hol_1_min_angle;
	  float hol_1_energy;
	  float hol_1_flag_all_shower;
	  float hol_1_min_length;
	  float hol_2_flag;
	  float hol_2_min_angle;
	  float hol_2_medium_dQ_dx;
	  float hol_2_ncount;
	  float hol_2_energy;
	  tree3->SetBranchAddress("hol_flag", &hol_flag);
	  tree3->SetBranchAddress("hol_1_flag", &hol_1_flag);
	  tree3->SetBranchAddress("hol_1_n_valid_tracks", &hol_1_n_valid_tracks);
	  tree3->SetBranchAddress("hol_1_min_angle", &hol_1_min_angle);
	  tree3->SetBranchAddress("hol_1_energy", &hol_1_energy);
	  tree3->SetBranchAddress("hol_1_flag_all_shower", &hol_1_flag_all_shower);
	  tree3->SetBranchAddress("hol_1_min_length", &hol_1_min_length);
	  tree3->SetBranchAddress("hol_2_flag", &hol_2_flag);
	  tree3->SetBranchAddress("hol_2_min_angle", &hol_2_min_angle);
	  tree3->SetBranchAddress("hol_2_medium_dQ_dx", &hol_2_medium_dQ_dx);
	  tree3->SetBranchAddress("hol_2_ncount", &hol_2_ncount);
	  tree3->SetBranchAddress("hol_2_energy", &hol_2_energy);


	  float lol_flag;
	  std::vector<float> *lol_1_v_flag= new std::vector<float>;
	  std::vector<float> *lol_1_v_energy = new std::vector<float>;
	  std::vector<float> *lol_1_v_vtx_n_segs= new std::vector<float>;
	  std::vector<float> *lol_1_v_nseg= new std::vector<float>;
	  std::vector<float> *lol_1_v_angle= new std::vector<float>;
	  std::vector<float> *lol_2_v_flag = new std::vector<float>;
	  std::vector<float> *lol_2_v_length= new std::vector<float>;
	  std::vector<float> *lol_2_v_angle= new std::vector<float>;
	  std::vector<float> *lol_2_v_type= new std::vector<float>;
	  std::vector<float> *lol_2_v_vtx_n_segs= new std::vector<float>;
	  std::vector<float> *lol_2_v_energy= new std::vector<float>;
	  std::vector<float> *lol_2_v_shower_main_length= new std::vector<float>;
	  std::vector<float> *lol_2_v_flag_dir_weak= new std::vector<float>;
	  float lol_3_flag;
	  float lol_3_angle_beam;
	  float lol_3_n_valid_tracks;
	  float lol_3_min_angle;
	  float lol_3_vtx_n_segs;
	  float lol_3_energy;
	  float lol_3_shower_main_length;
	  float lol_3_n_out;
	  float lol_3_n_sum;
	  tree3->SetBranchAddress("lol_flag",&lol_flag);
	  tree3->SetBranchAddress("lol_1_v_flag",&lol_1_v_flag);
	  tree3->SetBranchAddress("lol_1_v_energy",&lol_1_v_energy);
	  tree3->SetBranchAddress("lol_1_v_vtx_n_segs",&lol_1_v_vtx_n_segs);
	  tree3->SetBranchAddress("lol_1_v_nseg",&lol_1_v_nseg);
	  tree3->SetBranchAddress("lol_1_v_angle",&lol_1_v_angle);
	  tree3->SetBranchAddress("lol_2_v_flag",&lol_2_v_flag);
	  tree3->SetBranchAddress("lol_2_v_length",&lol_2_v_length);
	  tree3->SetBranchAddress("lol_2_v_angle",&lol_2_v_angle);
	  tree3->SetBranchAddress("lol_2_v_type",&lol_2_v_type);
	  tree3->SetBranchAddress("lol_2_v_vtx_n_segs",&lol_2_v_vtx_n_segs);
	  tree3->SetBranchAddress("lol_2_v_energy",&lol_2_v_energy);
	  tree3->SetBranchAddress("lol_2_v_shower_main_length",&lol_2_v_shower_main_length);
	  tree3->SetBranchAddress("lol_2_v_flag_dir_weak",&lol_2_v_flag_dir_weak);
	  tree3->SetBranchAddress("lol_3_flag",&lol_3_flag);
	  tree3->SetBranchAddress("lol_3_angle_beam",&lol_3_angle_beam);
	  tree3->SetBranchAddress("lol_3_n_valid_tracks",&lol_3_n_valid_tracks);
	  tree3->SetBranchAddress("lol_3_min_angle",&lol_3_min_angle);
	  tree3->SetBranchAddress("lol_3_vtx_n_segs",&lol_3_vtx_n_segs);
	  tree3->SetBranchAddress("lol_3_energy",&lol_3_energy);
	  tree3->SetBranchAddress("lol_3_shower_main_length",&lol_3_shower_main_length);
	  tree3->SetBranchAddress("lol_3_n_out",&lol_3_n_out);
	  tree3->SetBranchAddress("lol_3_n_sum",&lol_3_n_sum);

	  float cosmict_flag_1; // fiducial volume vertex
	  float cosmict_flag_2;  // single muon
	  float cosmict_flag_3;  // single muon (long)
	  float cosmict_flag_4;  // kinematics muon
	  float cosmict_flag_5; // kinematics muon (long)
	  float cosmict_flag_6; // special ...
	  float cosmict_flag_7;  // muon+ michel
	  float cosmict_flag_8;  // muon + michel + special
	  float cosmict_flag_9;  // this tagger is relevant for nueCC, see "cosmic tagger ones, one case of cosmics ..." (frist one ...)
	  std::vector<float> *cosmict_flag_10 = new std::vector<float>;  // front upstream (dirt)
	  float cosmict_flag;
	  float cosmict_2_filled;
	  float cosmict_2_particle_type;
	  float cosmict_2_n_muon_tracks;
	  float cosmict_2_total_shower_length;
	  float cosmict_2_flag_inside;
	  float cosmict_2_angle_beam;
	  float cosmict_2_flag_dir_weak;
	  float cosmict_2_dQ_dx_end;
	  float cosmict_2_dQ_dx_front;
	  float cosmict_2_theta;
	  float cosmict_2_phi;
	  float cosmict_2_valid_tracks;
	  float cosmict_3_filled;
	  float cosmict_3_flag_inside;
	  float cosmict_3_angle_beam;
	  float cosmict_3_flag_dir_weak;
	  float cosmict_3_dQ_dx_end;
	  float cosmict_3_dQ_dx_front;
	  float cosmict_3_theta;
	  float cosmict_3_phi;
	  float cosmict_3_valid_tracks;
	  float cosmict_4_filled;
	  float cosmict_4_flag_inside;
	  float cosmict_4_angle_beam;
	  float cosmict_4_connected_showers;  // need to be careful about the nueCC ...
	  float cosmict_5_filled;
	  float cosmict_5_flag_inside;
	  float cosmict_5_angle_beam;
	  float cosmict_5_connected_showers;
	  float cosmict_6_filled;
	  float cosmict_6_flag_dir_weak;
	  float cosmict_6_flag_inside;
	  float cosmict_6_angle;
	  float cosmict_7_filled;
	  float cosmict_7_flag_sec;
	  float cosmict_7_n_muon_tracks;
	  float cosmict_7_total_shower_length;
	  float cosmict_7_flag_inside;
	  float cosmict_7_angle_beam;
	  float cosmict_7_flag_dir_weak;
	  float cosmict_7_dQ_dx_end;
	  float cosmict_7_dQ_dx_front;
	  float cosmict_7_theta;
	  float cosmict_7_phi;
	  float cosmict_8_filled;
	  float cosmict_8_flag_out;
	  float cosmict_8_muon_length;
	  float cosmict_8_acc_length;
	  std::vector<float> *cosmict_10_flag_inside= new std::vector<float>;
	  std::vector<float> *cosmict_10_vtx_z= new std::vector<float>;
	  std::vector<float> *cosmict_10_flag_shower= new std::vector<float>;
	  std::vector<float> *cosmict_10_flag_dir_weak= new std::vector<float>;
	  std::vector<float> *cosmict_10_angle_beam= new std::vector<float>;
	  std::vector<float> *cosmict_10_length = new std::vector<float>;
	  tree3->SetBranchAddress("cosmict_flag_1",&cosmict_flag_1);
	  tree3->SetBranchAddress("cosmict_flag_2",&cosmict_flag_2);
	  tree3->SetBranchAddress("cosmict_flag_3",&cosmict_flag_3);
	  tree3->SetBranchAddress("cosmict_flag_4",&cosmict_flag_4);
	  tree3->SetBranchAddress("cosmict_flag_5",&cosmict_flag_5);
	  tree3->SetBranchAddress("cosmict_flag_6",&cosmict_flag_6);
	  tree3->SetBranchAddress("cosmict_flag_7",&cosmict_flag_7);
	  tree3->SetBranchAddress("cosmict_flag_8",&cosmict_flag_8);
	  tree3->SetBranchAddress("cosmict_flag_9",&cosmict_flag_9);
	  tree3->SetBranchAddress("cosmict_flag_10",&cosmict_flag_10);
	  tree3->SetBranchAddress("cosmict_flag",&cosmict_flag);
	  tree3->SetBranchAddress("cosmict_2_filled",&cosmict_2_filled);
	  tree3->SetBranchAddress("cosmict_2_particle_type",&cosmict_2_particle_type);
	  tree3->SetBranchAddress("cosmict_2_n_muon_tracks",&cosmict_2_n_muon_tracks);
	  tree3->SetBranchAddress("cosmict_2_total_shower_length",&cosmict_2_total_shower_length);
	  tree3->SetBranchAddress("cosmict_2_flag_inside",&cosmict_2_flag_inside);
	  tree3->SetBranchAddress("cosmict_2_angle_beam",&cosmict_2_angle_beam);
	  tree3->SetBranchAddress("cosmict_2_flag_dir_weak",&cosmict_2_flag_dir_weak);
	  tree3->SetBranchAddress("cosmict_2_dQ_dx_end",&cosmict_2_dQ_dx_end);
	  tree3->SetBranchAddress("cosmict_2_dQ_dx_front",&cosmict_2_dQ_dx_front);
	  tree3->SetBranchAddress("cosmict_2_theta",&cosmict_2_theta);
	  tree3->SetBranchAddress("cosmict_2_phi",&cosmict_2_phi);
	  tree3->SetBranchAddress("cosmict_2_valid_tracks",&cosmict_2_valid_tracks);
	  tree3->SetBranchAddress("cosmict_3_filled",&cosmict_3_filled);
	  tree3->SetBranchAddress("cosmict_3_flag_inside",&cosmict_3_flag_inside);
	  tree3->SetBranchAddress("cosmict_3_angle_beam",&cosmict_3_angle_beam);
	  tree3->SetBranchAddress("cosmict_3_flag_dir_weak",&cosmict_3_flag_dir_weak);
	  tree3->SetBranchAddress("cosmict_3_dQ_dx_end",&cosmict_3_dQ_dx_end);
	  tree3->SetBranchAddress("cosmict_3_dQ_dx_front",&cosmict_3_dQ_dx_front);
	  tree3->SetBranchAddress("cosmict_3_theta",&cosmict_3_theta);
	  tree3->SetBranchAddress("cosmict_3_phi",&cosmict_3_phi);
	  tree3->SetBranchAddress("cosmict_3_valid_tracks",&cosmict_3_valid_tracks);
	  tree3->SetBranchAddress("cosmict_4_filled",&cosmict_4_filled);
	  tree3->SetBranchAddress("cosmict_4_flag_inside",&cosmict_4_flag_inside);
	  tree3->SetBranchAddress("cosmict_4_angle_beam",&cosmict_4_angle_beam);
	  tree3->SetBranchAddress("cosmict_4_connected_showers",&cosmict_4_connected_showers);
	  tree3->SetBranchAddress("cosmict_5_filled",&cosmict_5_filled);
	  tree3->SetBranchAddress("cosmict_5_flag_inside",&cosmict_5_flag_inside);
	  tree3->SetBranchAddress("cosmict_5_angle_beam",&cosmict_5_angle_beam);
	  tree3->SetBranchAddress("cosmict_5_connected_showers",&cosmict_5_connected_showers);
	  tree3->SetBranchAddress("cosmict_6_filled",&cosmict_6_filled);
	  tree3->SetBranchAddress("cosmict_6_flag_dir_weak",&cosmict_6_flag_dir_weak);
	  tree3->SetBranchAddress("cosmict_6_flag_inside",&cosmict_6_flag_inside);
	  tree3->SetBranchAddress("cosmict_6_angle",&cosmict_6_angle);
	  tree3->SetBranchAddress("cosmict_7_filled",&cosmict_7_filled);
	  tree3->SetBranchAddress("cosmict_7_flag_sec",&cosmict_7_flag_sec);
	  tree3->SetBranchAddress("cosmict_7_n_muon_tracks",&cosmict_7_n_muon_tracks);
	  tree3->SetBranchAddress("cosmict_7_total_shower_length",&cosmict_7_total_shower_length);
	  tree3->SetBranchAddress("cosmict_7_flag_inside",&cosmict_7_flag_inside);
	  tree3->SetBranchAddress("cosmict_7_angle_beam",&cosmict_7_angle_beam);
	  tree3->SetBranchAddress("cosmict_7_flag_dir_weak",&cosmict_7_flag_dir_weak);
	  tree3->SetBranchAddress("cosmict_7_dQ_dx_end",&cosmict_7_dQ_dx_end);
	  tree3->SetBranchAddress("cosmict_7_dQ_dx_front",&cosmict_7_dQ_dx_front);
	  tree3->SetBranchAddress("cosmict_7_theta",&cosmict_7_theta);
	  tree3->SetBranchAddress("cosmict_7_phi",&cosmict_7_phi);
	  tree3->SetBranchAddress("cosmict_8_filled",&cosmict_8_filled);
	  tree3->SetBranchAddress("cosmict_8_flag_out",&cosmict_8_flag_out);
	  tree3->SetBranchAddress("cosmict_8_muon_length",&cosmict_8_muon_length);
	  tree3->SetBranchAddress("cosmict_8_acc_length",&cosmict_8_acc_length);
	  tree3->SetBranchAddress("cosmict_10_flag_inside",&cosmict_10_flag_inside);
	  tree3->SetBranchAddress("cosmict_10_vtx_z",&cosmict_10_vtx_z);
	  tree3->SetBranchAddress("cosmict_10_flag_shower",&cosmict_10_flag_shower);
	  tree3->SetBranchAddress("cosmict_10_flag_dir_weak",&cosmict_10_flag_dir_weak);
	  tree3->SetBranchAddress("cosmict_10_angle_beam",&cosmict_10_angle_beam);
	  tree3->SetBranchAddress("cosmict_10_length",&cosmict_10_length);

	  // numu tagger
	  //float numu_cc_flag;
	  std::vector<float> *numu_cc_flag_1= new std::vector<float>;
	  std::vector<float> *numu_cc_1_particle_type= new std::vector<float>;
	  std::vector<float> *numu_cc_1_length= new std::vector<float>;
	  std::vector<float> *numu_cc_1_medium_dQ_dx= new std::vector<float>;
	  std::vector<float> *numu_cc_1_dQ_dx_cut= new std::vector<float>;
	  std::vector<float> *numu_cc_1_direct_length= new std::vector<float>;
	  std::vector<float> *numu_cc_1_n_daughter_tracks= new std::vector<float>;
	  std::vector<float> *numu_cc_1_n_daughter_all= new std::vector<float>;
	  std::vector<float> *numu_cc_flag_2= new std::vector<float>;
	  std::vector<float> *numu_cc_2_length= new std::vector<float>;
	  std::vector<float> *numu_cc_2_total_length= new std::vector<float>;
	  std::vector<float> *numu_cc_2_n_daughter_tracks= new std::vector<float>;
	  std::vector<float> *numu_cc_2_n_daughter_all = new std::vector<float>;
	  float numu_cc_flag_3;
	  float numu_cc_3_particle_type;
	  float numu_cc_3_max_length;
	  float numu_cc_3_acc_track_length;
	  float numu_cc_3_max_length_all;
	  float numu_cc_3_max_muon_length;
	  float numu_cc_3_n_daughter_tracks;
	  float numu_cc_3_n_daughter_all;
	  tree3->SetBranchAddress("numu_cc_flag",&numu_cc_flag);
	  tree3->SetBranchAddress("numu_cc_flag_1",&numu_cc_flag_1);
	  tree3->SetBranchAddress("numu_cc_1_particle_type",&numu_cc_1_particle_type);
	  tree3->SetBranchAddress("numu_cc_1_length",&numu_cc_1_length);
	  tree3->SetBranchAddress("numu_cc_1_medium_dQ_dx",&numu_cc_1_medium_dQ_dx);
	  tree3->SetBranchAddress("numu_cc_1_dQ_dx_cut",&numu_cc_1_dQ_dx_cut);
	  tree3->SetBranchAddress("numu_cc_1_direct_length",&numu_cc_1_direct_length);
	  tree3->SetBranchAddress("numu_cc_1_n_daughter_tracks",&numu_cc_1_n_daughter_tracks);
	  tree3->SetBranchAddress("numu_cc_1_n_daughter_all",&numu_cc_1_n_daughter_all);
	  tree3->SetBranchAddress("numu_cc_flag_2",&numu_cc_flag_2);
	  tree3->SetBranchAddress("numu_cc_2_length",&numu_cc_2_length);
	  tree3->SetBranchAddress("numu_cc_2_total_length",&numu_cc_2_total_length);
	  tree3->SetBranchAddress("numu_cc_2_n_daughter_tracks",&numu_cc_2_n_daughter_tracks);
	  tree3->SetBranchAddress("numu_cc_2_n_daughter_all",&numu_cc_2_n_daughter_all);
	  tree3->SetBranchAddress("numu_cc_flag_3",&numu_cc_flag_3);
	  tree3->SetBranchAddress("numu_cc_3_particle_type",&numu_cc_3_particle_type);
	  tree3->SetBranchAddress("numu_cc_3_max_length",&numu_cc_3_max_length);
	  tree3->SetBranchAddress("numu_cc_3_track_length",&numu_cc_3_acc_track_length);
	  tree3->SetBranchAddress("numu_cc_3_max_length_all",&numu_cc_3_max_length_all);
	  tree3->SetBranchAddress("numu_cc_3_max_muon_length",&numu_cc_3_max_muon_length);
	  tree3->SetBranchAddress("numu_cc_3_n_daughter_tracks",&numu_cc_3_n_daughter_tracks);
	  tree3->SetBranchAddress("numu_cc_3_n_daughter_all",&numu_cc_3_n_daughter_all);

	  // BDT scores
	  float cosmict_2_4_score;
	  float cosmict_3_5_score;
	  float cosmict_6_score;
	  float cosmict_7_score;
	  float cosmict_8_score;
	  // vector ...
	  float cosmict_10_score;
	  // vector
	  float numu_1_score;
	  float numu_2_score;
	  // scalar
	  float numu_3_score;
	  // total one
	  float cosmict_score;
	  float numu_score;
	  // nue BDTs
	  float mipid_score;
	  float gap_score;
	  float hol_lol_score;
	  float cme_anc_score;
	  float mgo_mgt_score;
	  float br1_score;
	  float br3_score;
	  float br3_3_score;
	  float br3_5_score;
	  float br3_6_score;
	  float stemdir_br2_score;
	  float trimuon_score;
	  float br4_tro_score;
	  float mipquality_score;
	  float pio_1_score;
	  float pio_2_score;
	  float stw_spt_score;
	  float vis_1_score;
	  float vis_2_score;
	  float stw_2_score;
	  float stw_3_score;
	  float stw_4_score;
	  float sig_1_score;
	  float sig_2_score;
	  float lol_1_score;
	  float lol_2_score;
	  float tro_1_score;
	  float tro_2_score;
	  float tro_4_score;
	  float tro_5_score;
	  float nue_score;
	  tree3->SetBranchAddress("cosmict_2_4_score",&cosmict_2_4_score);
	  tree3->SetBranchAddress("cosmict_3_5_score",&cosmict_3_5_score);
	  tree3->SetBranchAddress("cosmict_6_score",&cosmict_6_score);
	  tree3->SetBranchAddress("cosmict_7_score",&cosmict_7_score);
	  tree3->SetBranchAddress("cosmict_8_score",&cosmict_8_score);
	  tree3->SetBranchAddress("cosmict_10_score",&cosmict_10_score);
	  tree3->SetBranchAddress("numu_1_score",&numu_1_score);
	  tree3->SetBranchAddress("numu_2_score",&numu_2_score);
	  tree3->SetBranchAddress("numu_3_score",&numu_3_score);
	  tree3->SetBranchAddress("cosmict_score",&cosmict_score);
	  tree3->SetBranchAddress("numu_score",&numu_score);
	  tree3->SetBranchAddress("mipid_score",&mipid_score);
	  tree3->SetBranchAddress("gap_score",&gap_score);
	  tree3->SetBranchAddress("hol_lol_score",&hol_lol_score);
	  tree3->SetBranchAddress("cme_anc_score",&cme_anc_score);
	  tree3->SetBranchAddress("mgo_mgt_score",&mgo_mgt_score);
	  tree3->SetBranchAddress("br1_score",&br1_score);
	  tree3->SetBranchAddress("br3_score",&br3_score);
	  tree3->SetBranchAddress("br3_3_score",&br3_3_score);
	  tree3->SetBranchAddress("br3_5_score",&br3_5_score);
	  tree3->SetBranchAddress("br3_6_score",&br3_6_score);
	  tree3->SetBranchAddress("stemdir_br2_score",&stemdir_br2_score);
	  tree3->SetBranchAddress("trimuon_score",&trimuon_score);
	  tree3->SetBranchAddress("br4_tro_score",&br4_tro_score);
	  tree3->SetBranchAddress("mipquality_score",&mipquality_score);
	  tree3->SetBranchAddress("pio_1_score",&pio_1_score);
	  tree3->SetBranchAddress("pio_2_score",&pio_2_score);
	  tree3->SetBranchAddress("stw_spt_score",&stw_spt_score);
	  tree3->SetBranchAddress("vis_1_score",&vis_1_score);
	  tree3->SetBranchAddress("vis_2_score",&vis_2_score);
	  tree3->SetBranchAddress("stw_2_score",&stw_2_score);
	  tree3->SetBranchAddress("stw_3_score",&stw_3_score);
	  tree3->SetBranchAddress("stw_4_score",&stw_4_score);
	  tree3->SetBranchAddress("sig_1_score",&sig_1_score);
	  tree3->SetBranchAddress("sig_2_score",&sig_2_score);
	  tree3->SetBranchAddress("lol_1_score",&lol_1_score);
	  tree3->SetBranchAddress("lol_2_score",&lol_2_score);
	  tree3->SetBranchAddress("tro_1_score",&tro_1_score);
	  tree3->SetBranchAddress("tro_2_score",&tro_2_score);
	  tree3->SetBranchAddress("tro_4_score",&tro_4_score);
	  tree3->SetBranchAddress("tro_5_score",&tro_5_score);
	  tree3->SetBranchAddress("nue_score",&nue_score);



  /// Read and assign values
  tree3->GetEntry(0); // rare case: multiple in-beam matched activity

  nsm::NuSelectionBDT::stkdar _stkdar_init = {
          ssm_flag_st_kdar,
          ssm_Nsm,
          ssm_Nsm_wivtx,
          ssm_dq_dx_fwd_1,
          ssm_dq_dx_fwd_2,
          ssm_dq_dx_fwd_3,
          ssm_dq_dx_fwd_4,
          ssm_dq_dx_fwd_5,
          ssm_dq_dx_bck_1,
          ssm_dq_dx_bck_2,
          ssm_dq_dx_bck_3,
          ssm_dq_dx_bck_4,
          ssm_dq_dx_bck_5,
          ssm_d_dq_dx_fwd_12,
          ssm_d_dq_dx_fwd_23,
          ssm_d_dq_dx_fwd_34,
          ssm_d_dq_dx_fwd_45,
          ssm_d_dq_dx_bck_12,
          ssm_d_dq_dx_bck_23,
          ssm_d_dq_dx_bck_34,
          ssm_d_dq_dx_bck_45,
          ssm_max_dq_dx_fwd_3,
          ssm_max_dq_dx_fwd_5,
          ssm_max_dq_dx_bck_3,
          ssm_max_dq_dx_bck_5,
          ssm_max_d_dq_dx_fwd_3,
          ssm_max_d_dq_dx_fwd_5,
          ssm_max_d_dq_dx_bck_3,
          ssm_max_d_dq_dx_bck_5,
          ssm_medium_dq_dx,
          ssm_medium_dq_dx_bp,
          ssm_angle_to_z,
          ssm_angle_to_target,
          ssm_angle_to_absorber,
          ssm_angle_to_vertical,
          ssm_x_dir,
          ssm_y_dir,
          ssm_z_dir,
          ssm_kine_energy,
          ssm_kine_energy_reduced,
          ssm_vtx_activity,
          ssm_pdg,
          ssm_dQ_dx_cut,
          ssm_score_mu_fwd,
          ssm_score_p_fwd,
          ssm_score_e_fwd,
          ssm_score_mu_bck,
          ssm_score_p_bck,
          ssm_score_e_bck,
          ssm_score_mu_fwd_bp,
          ssm_score_p_fwd_bp,
          ssm_score_e_fwd_bp,
          ssm_length,
          ssm_direct_length,
          ssm_length_ratio,
          ssm_max_dev,
          ssm_n_prim_tracks_1,
          ssm_n_prim_tracks_3,
          ssm_n_prim_tracks_5,
          ssm_n_prim_tracks_8,
          ssm_n_prim_tracks_11,
          ssm_n_all_tracks_1,
          ssm_n_all_tracks_3,
          ssm_n_all_tracks_5,
          ssm_n_all_tracks_8,
          ssm_n_all_tracks_11,
          ssm_n_daughter_tracks_1,
          ssm_n_daughter_tracks_3,
          ssm_n_daughter_tracks_5,
          ssm_n_daughter_tracks_8,
          ssm_n_daughter_tracks_11,
          ssm_n_daughter_all_1,
          ssm_n_daughter_all_3,
          ssm_n_daughter_all_5,
          ssm_n_daughter_all_8,
          ssm_n_daughter_all_11,
          ssm_prim_track1_pdg,
          ssm_prim_track1_score_mu_fwd,
          ssm_prim_track1_score_p_fwd,
          ssm_prim_track1_score_e_fwd,
          ssm_prim_track1_score_mu_bck,
          ssm_prim_track1_score_p_bck,
          ssm_prim_track1_score_e_bck,
          ssm_prim_track1_length,
          ssm_prim_track1_direct_length,
          ssm_prim_track1_length_ratio,
          ssm_prim_track1_max_dev,
          ssm_prim_track1_kine_energy_range,
          ssm_prim_track1_kine_energy_range_mu,
          ssm_prim_track1_kine_energy_range_p,
          ssm_prim_track1_kine_energy_range_e,
          ssm_prim_track1_kine_energy_cal,
          ssm_prim_track1_medium_dq_dx,
          ssm_prim_track1_x_dir,
          ssm_prim_track1_y_dir,
          ssm_prim_track1_z_dir,
          ssm_prim_track1_add_daught_track_counts_1,
          ssm_prim_track1_add_daught_all_counts_1,
          ssm_prim_track1_add_daught_track_counts_5,
          ssm_prim_track1_add_daught_all_counts_5,
          ssm_prim_track1_add_daught_track_counts_11,
          ssm_prim_track1_add_daught_all_counts_11,
	  ssm_prim_track2_pdg,
          ssm_prim_track2_score_mu_fwd,
          ssm_prim_track2_score_p_fwd,
          ssm_prim_track2_score_e_fwd,
          ssm_prim_track2_score_mu_bck,
          ssm_prim_track2_score_p_bck,
          ssm_prim_track2_score_e_bck,
          ssm_prim_track2_length,
          ssm_prim_track2_direct_length,
          ssm_prim_track2_length_ratio,
          ssm_prim_track2_max_dev,
          ssm_prim_track2_kine_energy_range,
          ssm_prim_track2_kine_energy_range_mu,
          ssm_prim_track2_kine_energy_range_p,
          ssm_prim_track2_kine_energy_range_e,
          ssm_prim_track2_kine_energy_cal,
          ssm_prim_track2_medium_dq_dx,
          ssm_prim_track2_x_dir,
          ssm_prim_track2_y_dir,
          ssm_prim_track2_z_dir,
          ssm_prim_track2_add_daught_track_counts_1,
          ssm_prim_track2_add_daught_all_counts_1,
          ssm_prim_track2_add_daught_track_counts_5,
          ssm_prim_track2_add_daught_all_counts_5,
          ssm_prim_track2_add_daught_track_counts_11,
          ssm_prim_track2_add_daught_all_counts_11,
          ssm_daught_track1_pdg,
          ssm_daught_track1_score_mu_fwd,
          ssm_daught_track1_score_p_fwd,
          ssm_daught_track1_score_e_fwd,
          ssm_daught_track1_score_mu_bck,
          ssm_daught_track1_score_p_bck,
          ssm_daught_track1_score_e_bck,
          ssm_daught_track1_length,
          ssm_daught_track1_direct_length,
          ssm_daught_track1_length_ratio,
          ssm_daught_track1_max_dev,
          ssm_daught_track1_kine_energy_range,
          ssm_daught_track1_kine_energy_range_mu,
          ssm_daught_track1_kine_energy_range_p,
          ssm_daught_track1_kine_energy_range_e,
          ssm_daught_track1_kine_energy_cal,
          ssm_daught_track1_medium_dq_dx,
          ssm_daught_track1_x_dir,
          ssm_daught_track1_y_dir,
          ssm_daught_track1_z_dir,
          ssm_daught_track1_add_daught_track_counts_1,
          ssm_daught_track1_add_daught_all_counts_1,
          ssm_daught_track1_add_daught_track_counts_5,
          ssm_daught_track1_add_daught_all_counts_5,
          ssm_daught_track1_add_daught_track_counts_11,
          ssm_daught_track1_add_daught_all_counts_11,
          ssm_daught_track2_pdg,
          ssm_daught_track2_score_mu_fwd,
          ssm_daught_track2_score_p_fwd,
          ssm_daught_track2_score_e_fwd,
          ssm_daught_track2_score_mu_bck,
          ssm_daught_track2_score_p_bck,
          ssm_daught_track2_score_e_bck,
          ssm_daught_track2_length,
          ssm_daught_track2_direct_length,
          ssm_daught_track2_length_ratio,
          ssm_daught_track2_max_dev,
          ssm_daught_track2_kine_energy_range,
          ssm_daught_track2_kine_energy_range_mu,
          ssm_daught_track2_kine_energy_range_p,
          ssm_daught_track2_kine_energy_range_e,
          ssm_daught_track2_kine_energy_cal,
          ssm_daught_track2_medium_dq_dx,
          ssm_daught_track2_x_dir,
          ssm_daught_track2_y_dir,
          ssm_daught_track2_z_dir,
          ssm_daught_track2_add_daught_track_counts_1,
          ssm_daught_track2_add_daught_all_counts_1,
          ssm_daught_track2_add_daught_track_counts_5,
          ssm_daught_track2_add_daught_all_counts_5,
          ssm_daught_track2_add_daught_track_counts_11,
          ssm_daught_track2_add_daught_all_counts_11,
          ssm_prim_shw1_pdg,
          ssm_prim_shw1_score_mu_fwd,
          ssm_prim_shw1_score_p_fwd,
          ssm_prim_shw1_score_e_fwd,
          ssm_prim_shw1_score_mu_bck,
          ssm_prim_shw1_score_p_bck,
          ssm_prim_shw1_score_e_bck,
          ssm_prim_shw1_length,
          ssm_prim_shw1_direct_length,
          ssm_prim_shw1_length_ratio,
          ssm_prim_shw1_max_dev,
          ssm_prim_shw1_kine_energy_range,
          ssm_prim_shw1_kine_energy_range_mu,
          ssm_prim_shw1_kine_energy_range_p,
          ssm_prim_shw1_kine_energy_range_e,
          ssm_prim_shw1_kine_energy_cal,
          ssm_prim_shw1_kine_energy_best,
          ssm_prim_shw1_medium_dq_dx,
          ssm_prim_shw1_x_dir,
          ssm_prim_shw1_y_dir,
          ssm_prim_shw1_z_dir,
          ssm_prim_shw1_add_daught_track_counts_1,
          ssm_prim_shw1_add_daught_all_counts_1,
          ssm_prim_shw1_add_daught_track_counts_5,
          ssm_prim_shw1_add_daught_all_counts_5,
          ssm_prim_shw1_add_daught_track_counts_11,
          ssm_prim_shw1_add_daught_all_counts_11,
          ssm_prim_shw2_pdg,
          ssm_prim_shw2_score_mu_fwd,
          ssm_prim_shw2_score_p_fwd,
          ssm_prim_shw2_score_e_fwd,
          ssm_prim_shw2_score_mu_bck,
          ssm_prim_shw2_score_p_bck,
          ssm_prim_shw2_score_e_bck,
          ssm_prim_shw2_length,
          ssm_prim_shw2_direct_length,
          ssm_prim_shw2_length_ratio,
          ssm_prim_shw2_max_dev,
          ssm_prim_shw2_kine_energy_range,
          ssm_prim_shw2_kine_energy_range_mu,
          ssm_prim_shw2_kine_energy_range_p,
          ssm_prim_shw2_kine_energy_range_e,
          ssm_prim_shw2_kine_energy_cal,
          ssm_prim_shw2_kine_energy_best,
          ssm_prim_shw2_medium_dq_dx,
          ssm_prim_shw2_x_dir,
          ssm_prim_shw2_y_dir,
          ssm_prim_shw2_z_dir,
          ssm_prim_shw2_add_daught_track_counts_1,
          ssm_prim_shw2_add_daught_all_counts_1,
          ssm_prim_shw2_add_daught_track_counts_5,
          ssm_prim_shw2_add_daught_all_counts_5,
          ssm_prim_shw2_add_daught_track_counts_11,
          ssm_prim_shw2_add_daught_all_counts_11,
          ssm_daught_shw1_pdg,
          ssm_daught_shw1_score_mu_fwd,
          ssm_daught_shw1_score_p_fwd,
          ssm_daught_shw1_score_e_fwd,
          ssm_daught_shw1_score_mu_bck,
          ssm_daught_shw1_score_p_bck,
          ssm_daught_shw1_score_e_bck,
          ssm_daught_shw1_length,
          ssm_daught_shw1_direct_length,
          ssm_daught_shw1_length_ratio,
          ssm_daught_shw1_max_dev,
          ssm_daught_shw1_kine_energy_range,
          ssm_daught_shw1_kine_energy_range_mu,
          ssm_daught_shw1_kine_energy_range_p,
          ssm_daught_shw1_kine_energy_range_e,
          ssm_daught_shw1_kine_energy_cal,
          ssm_daught_shw1_kine_energy_best,
          ssm_daught_shw1_medium_dq_dx,
          ssm_daught_shw1_x_dir,
          ssm_daught_shw1_y_dir,
          ssm_daught_shw1_z_dir,
          ssm_daught_shw1_add_daught_track_counts_1,
          ssm_daught_shw1_add_daught_all_counts_1,
          ssm_daught_shw1_add_daught_track_counts_5,
          ssm_daught_shw1_add_daught_all_counts_5,
          ssm_daught_shw1_add_daught_track_counts_11,
          ssm_daught_shw1_add_daught_all_counts_11,
          ssm_daught_shw2_pdg,
          ssm_daught_shw2_score_mu_fwd,
          ssm_daught_shw2_score_p_fwd,
          ssm_daught_shw2_score_e_fwd,
          ssm_daught_shw2_score_mu_bck,
          ssm_daught_shw2_score_p_bck,
          ssm_daught_shw2_score_e_bck,
          ssm_daught_shw2_length,
          ssm_daught_shw2_direct_length,
          ssm_daught_shw2_length_ratio,
          ssm_daught_shw2_max_dev,
          ssm_daught_shw2_kine_energy_range,
          ssm_daught_shw2_kine_energy_range_mu,
          ssm_daught_shw2_kine_energy_range_p,
          ssm_daught_shw2_kine_energy_range_e,
          ssm_daught_shw2_kine_energy_cal,
          ssm_daught_shw2_kine_energy_best,
          ssm_daught_shw2_medium_dq_dx,
          ssm_daught_shw2_x_dir,
          ssm_daught_shw2_y_dir,
          ssm_daught_shw2_z_dir,
          ssm_daught_shw2_add_daught_track_counts_1,
          ssm_daught_shw2_add_daught_all_counts_1,
          ssm_daught_shw2_add_daught_track_counts_5,
          ssm_daught_shw2_add_daught_all_counts_5,
          ssm_daught_shw2_add_daught_track_counts_11,
          ssm_daught_shw2_add_daught_all_counts_11,
          ssm_nu_angle_z,
          ssm_nu_angle_target,
          ssm_nu_angle_absorber,
          ssm_nu_angle_vertical,
          ssm_con_nu_angle_z,
          ssm_con_nu_angle_target,
          ssm_con_nu_angle_absorber,
          ssm_con_nu_angle_vertical,
          ssm_prim_nu_angle_z,
          ssm_prim_nu_angle_target,
          ssm_prim_nu_angle_absorber,
          ssm_prim_nu_angle_vertical,
          ssm_track_angle_z,
          ssm_track_angle_target,
          ssm_track_angle_absorber,
          ssm_track_angle_vertical,
          ssm_vtxX,
	  ssm_vtxY,
	  ssm_vtxZ,
	  ssm_offvtx_length,
          ssm_offvtx_energy,
	  ssm_n_offvtx_tracks_1,
          ssm_n_offvtx_tracks_3,
          ssm_n_offvtx_tracks_5,
          ssm_n_offvtx_tracks_8,
          ssm_n_offvtx_tracks_11,
          ssm_n_offvtx_showers_1,
          ssm_n_offvtx_showers_3,
          ssm_n_offvtx_showers_5,
          ssm_n_offvtx_showers_8,
          ssm_n_offvtx_showers_11,
          ssm_offvtx_track1_pdg,
          ssm_offvtx_track1_score_mu_fwd,
          ssm_offvtx_track1_score_p_fwd,
          ssm_offvtx_track1_score_e_fwd,
          ssm_offvtx_track1_score_mu_bck,
          ssm_offvtx_track1_score_p_bck,
          ssm_offvtx_track1_score_e_bck,
          ssm_offvtx_track1_length,
          ssm_offvtx_track1_direct_length,
          ssm_offvtx_track1_max_dev,
          ssm_offvtx_track1_kine_energy_range,
          ssm_offvtx_track1_kine_energy_range_mu,
          ssm_offvtx_track1_kine_energy_range_p,
          ssm_offvtx_track1_kine_energy_range_e,
          ssm_offvtx_track1_kine_energy_cal,
          ssm_offvtx_track1_medium_dq_dx,
          ssm_offvtx_track1_x_dir,
          ssm_offvtx_track1_y_dir,
          ssm_offvtx_track1_z_dir,
          ssm_offvtx_track1_dist_mainvtx,
          ssm_offvtx_shw1_pdg_offvtx,
          ssm_offvtx_shw1_score_mu_fwd,
          ssm_offvtx_shw1_score_p_fwd,
          ssm_offvtx_shw1_score_e_fwd,
          ssm_offvtx_shw1_score_mu_bck,
          ssm_offvtx_shw1_score_p_bck,
          ssm_offvtx_shw1_score_e_bck,
          ssm_offvtx_shw1_length,
          ssm_offvtx_shw1_direct_length,
          ssm_offvtx_shw1_max_dev,
          ssm_offvtx_shw1_kine_energy_best,
          ssm_offvtx_shw1_kine_energy_range,
          ssm_offvtx_shw1_kine_energy_range_mu,
          ssm_offvtx_shw1_kine_energy_range_p,
          ssm_offvtx_shw1_kine_energy_range_e,
          ssm_offvtx_shw1_kine_energy_cal,
          ssm_offvtx_shw1_medium_dq_dx,
          ssm_offvtx_shw1_x_dir,
          ssm_offvtx_shw1_y_dir,
          ssm_offvtx_shw1_z_dir,
	  ssm_offvtx_shw1_dist_mainvtx,
          ssmsp_Ntrack,
          ssmsp_Nsp,
          ssmsp_Nsp_tot,
          ssmsp_pdg,
          ssmsp_id,
          ssmsp_mother,
          ssmsp_x,
          ssmsp_y,
          ssmsp_z,
          ssmsp_dx,
          ssmsp_dQ,
          ssmsp_KE,
          ssmsp_containing_shower_id,
	  ssmsp_containing_shower_ke,
	  ssmsp_containing_shower_flag,
          ssm_kine_reco_Enu,
          ssm_kine_reco_add_energy,
          ssm_kine_energy_particle,
          ssm_kine_energy_info,
          ssm_kine_particle_type,
          ssm_kine_energy_included,
          ssm_kine_pio_mass,
          ssm_kine_pio_flag,
          ssm_kine_pio_vtx_dis,
          ssm_kine_pio_energy_1,
          ssm_kine_pio_theta_1,
          ssm_kine_pio_phi_1,
          ssm_kine_pio_dis_1,
          ssm_kine_pio_energy_2,
          ssm_kine_pio_theta_2,
          ssm_kine_pio_phi_2,
          ssm_kine_pio_dis_2,
          ssm_kine_pio_angle,
          numu_cc_flag,
	  cosmict_flag_1, // fiducial volume vertex
          cosmict_flag_2,  // single muon
          cosmict_flag_3,  // single muon (long)
          cosmict_flag_4,  // kinematics muon
          cosmict_flag_5, // kinematics muon (long)
          cosmict_flag_6, // special ...
          cosmict_flag_7,  // muon+ michel
          cosmict_flag_8,  // muon + michel + special
          cosmict_flag_9,  // this tagger is relevant for nueCC, see "cosmic tagger ones, one case of cosmics ..." (frist one ...)
          cosmict_flag_10,  // front upstream (dirt)
          cosmict_flag
  };

  nsm::NuSelectionBDT::SPID _SPID_init = {
          shw_sp_num_mip_tracks,
          shw_sp_num_muons,
          shw_sp_num_pions,
          shw_sp_num_protons,
          shw_sp_proton_length_1,
          shw_sp_proton_dqdx_1,
          shw_sp_proton_energy_1,
          shw_sp_proton_length_2,
          shw_sp_proton_dqdx_2,
          shw_sp_proton_energy_2,
          shw_sp_n_good_showers,
          shw_sp_n_20mev_showers,
          shw_sp_n_br1_showers,
          shw_sp_n_br2_showers,
          shw_sp_n_br3_showers,
          shw_sp_n_br4_showers,
          shw_sp_n_20br1_showers,
          shw_sp_20mev_showers,
          shw_sp_br1_showers,
          shw_sp_br2_showers,
          shw_sp_br3_showers,
          shw_sp_br4_showers,
          shw_sp_shw_vtx_dis,
          shw_sp_max_shw_dis
  };

  nsm::NuSelectionBDT::SPSHWID1 _SPSHWID1_init = {
          shw_sp_filled,
          shw_sp_flag,
          shw_sp_energy,
          shw_sp_vec_dQ_dx_0,
          shw_sp_vec_dQ_dx_1,
          shw_sp_max_dQ_dx_sample,
          shw_sp_n_below_threshold,
          shw_sp_n_below_zero,
          shw_sp_n_lowest,
          shw_sp_n_highest,
          shw_sp_lowest_dQ_dx,
          shw_sp_highest_dQ_dx,
          shw_sp_medium_dQ_dx,
          shw_sp_stem_length,
          shw_sp_length_main,
          shw_sp_length_total,
          shw_sp_angle_beam,
          shw_sp_iso_angle,
          shw_sp_n_vertex,
          shw_sp_n_good_tracks,
          shw_sp_E_indirect_max_energy,
          shw_sp_flag_all_above,
          shw_sp_min_dQ_dx_5,
          shw_sp_n_other_vertex,
          shw_sp_n_stem_size,
          shw_sp_flag_stem_trajectory,
          shw_sp_min_dis
  };

  nsm::NuSelectionBDT::SPSHWID2 _SPSHWID2_init = {
          shw_sp_vec_median_dedx,
          shw_sp_vec_mean_dedx,
          shw_sp_vec_dQ_dx_2,
          shw_sp_vec_dQ_dx_3,
          shw_sp_vec_dQ_dx_4,
          shw_sp_vec_dQ_dx_5,
          shw_sp_vec_dQ_dx_6,
          shw_sp_vec_dQ_dx_7,
          shw_sp_vec_dQ_dx_8,
          shw_sp_vec_dQ_dx_9,
          shw_sp_vec_dQ_dx_10,
          shw_sp_vec_dQ_dx_11,
          shw_sp_vec_dQ_dx_12,
          shw_sp_vec_dQ_dx_13,
          shw_sp_vec_dQ_dx_14,
          shw_sp_vec_dQ_dx_15,
          shw_sp_vec_dQ_dx_16,
          shw_sp_vec_dQ_dx_17,
          shw_sp_vec_dQ_dx_18,
          shw_sp_vec_dQ_dx_19
  };

  nsm::NuSelectionBDT::SPPi0Tagger1 _SPPi0Tagger1_init = {
          shw_sp_pio_filled,
          shw_sp_pio_flag,
          shw_sp_pio_mip_id,
          shw_sp_pio_flag_pio,
          shw_sp_pio_1_flag,
          shw_sp_pio_1_mass,
          shw_sp_pio_1_pio_type,
          shw_sp_pio_1_energy_1,
          shw_sp_pio_1_energy_2,
          shw_sp_pio_1_dis_1,
          shw_sp_pio_1_dis_2,
          shw_sp_pio_2_v_flag,
          shw_sp_pio_2_v_dis2,
          shw_sp_pio_2_v_angle2,
          shw_sp_pio_2_v_acc_length
  };

  nsm::NuSelectionBDT::SPLowEMichel _SPLowEMichel_init = {
          shw_sp_lem_flag,
          shw_sp_lem_shower_total_length,
          shw_sp_lem_shower_main_length,
          shw_sp_lem_n_3seg,
          shw_sp_lem_e_charge,
          shw_sp_lem_e_dQdx,
          shw_sp_lem_shower_num_segs,
          shw_sp_lem_shower_num_main_segs
  };

  nsm::NuSelectionBDT::SPBadReco1 _SPBadReco1_init = {
          shw_sp_br_filled,
          shw_sp_br1_flag,
          shw_sp_br1_1_flag,
          shw_sp_br1_1_shower_type,
          shw_sp_br1_1_vtx_n_segs,
          shw_sp_br1_1_energy,
          shw_sp_br1_1_n_segs,
          shw_sp_br1_1_flag_sg_topology,
          shw_sp_br1_1_flag_sg_trajectory,
          shw_sp_br1_1_sg_length,
          shw_sp_br1_2_flag,
          shw_sp_br1_2_energy,
          shw_sp_br1_2_n_connected,
          shw_sp_br1_2_max_length,
          shw_sp_br1_2_n_connected_1,
          shw_sp_br1_2_vtx_n_segs,
          shw_sp_br1_2_n_shower_segs,
          shw_sp_br1_2_max_length_ratio,
          shw_sp_br1_2_shower_length,
          shw_sp_br1_3_flag,
          shw_sp_br1_3_energy,
          shw_sp_br1_3_n_connected_p,
          shw_sp_br1_3_max_length_p,
          shw_sp_br1_3_n_shower_segs,
          shw_sp_br1_3_flag_sg_topology,
          shw_sp_br1_3_flag_sg_trajectory,
          shw_sp_br1_3_n_shower_main_segs,
          shw_sp_br1_3_sg_length
  };
  nsm::NuSelectionBDT::SPBadReco2 _SPBadReco2_init = {
          shw_sp_br_filled,
          shw_sp_br2_flag,
          shw_sp_br2_flag_single_shower,
          shw_sp_br2_num_valid_tracks,
          shw_sp_br2_energy,
          shw_sp_br2_angle1,
          shw_sp_br2_angle2,
          shw_sp_br2_angle,
          shw_sp_br2_angle3,
          shw_sp_br2_n_shower_main_segs,
          shw_sp_br2_max_angle,
          shw_sp_br2_sg_length,
          shw_sp_br2_flag_sg_trajectory
  };
  nsm::NuSelectionBDT::SPBadReco3 _SPBadReco3_init = {
          shw_sp_br_filled,
          shw_sp_br3_flag,
          shw_sp_br3_1_flag,
          shw_sp_br3_1_energy,
          shw_sp_br3_1_n_shower_segments,
          shw_sp_br3_1_sg_flag_trajectory,
          shw_sp_br3_1_sg_direct_length,
          shw_sp_br3_1_sg_length,
          shw_sp_br3_1_total_main_length,
          shw_sp_br3_1_total_length,
          shw_sp_br3_1_iso_angle,
          shw_sp_br3_1_sg_flag_topology,
          shw_sp_br3_2_flag,
          shw_sp_br3_2_n_ele,
          shw_sp_br3_2_n_other,
          shw_sp_br3_2_energy,
          shw_sp_br3_2_total_main_length,
          shw_sp_br3_2_total_length,
          shw_sp_br3_2_other_fid,
          shw_sp_br3_3_v_flag,
          shw_sp_br3_3_v_energy,
          shw_sp_br3_3_v_angle,
          shw_sp_br3_3_v_dir_length,
          shw_sp_br3_3_v_length,
          shw_sp_br3_4_flag,
          shw_sp_br3_4_acc_length,
          shw_sp_br3_4_total_length,
          shw_sp_br3_4_energy,
          shw_sp_br3_5_v_flag,
          shw_sp_br3_5_v_dir_length,
          shw_sp_br3_5_v_total_length,
          shw_sp_br3_5_v_flag_avoid_muon_check,
          shw_sp_br3_5_v_n_seg,
          shw_sp_br3_5_v_angle,
          shw_sp_br3_5_v_sg_length,
          shw_sp_br3_5_v_energy,
          shw_sp_br3_5_v_n_main_segs,
          shw_sp_br3_5_v_n_segs,
          shw_sp_br3_5_v_shower_main_length,
          shw_sp_br3_5_v_shower_total_length,
          shw_sp_br3_6_v_flag,
          shw_sp_br3_6_v_angle,
          shw_sp_br3_6_v_angle1,
          shw_sp_br3_6_v_flag_shower_trajectory,
          shw_sp_br3_6_v_direct_length,
          shw_sp_br3_6_v_length,
          shw_sp_br3_6_v_n_other_vtx_segs,
          shw_sp_br3_6_v_energy,
          shw_sp_br3_7_flag,
          shw_sp_br3_7_energy,
          shw_sp_br3_7_min_angle,
          shw_sp_br3_7_sg_length,
          shw_sp_br3_7_shower_main_length,
          shw_sp_br3_8_flag,
          shw_sp_br3_8_max_dQ_dx,
          shw_sp_br3_8_energy,
          shw_sp_br3_8_n_main_segs,
          shw_sp_br3_8_shower_main_length,
          shw_sp_br3_8_shower_length
  };
  nsm::NuSelectionBDT::SPBadReco4 _SPBadReco4_init = {
          shw_sp_br_filled,
          shw_sp_br4_flag,
          shw_sp_br4_1_flag,
          shw_sp_br4_1_shower_main_length,
          shw_sp_br4_1_shower_total_length,
          shw_sp_br4_1_min_dis,
          shw_sp_br4_1_energy,
          shw_sp_br4_1_flag_avoid_muon_check,
          shw_sp_br4_1_n_vtx_segs,
          shw_sp_br4_1_n_main_segs,
          shw_sp_br4_2_flag,
          shw_sp_br4_2_ratio_45,
          shw_sp_br4_2_ratio_35,
          shw_sp_br4_2_ratio_25,
          shw_sp_br4_2_ratio_15,
          shw_sp_br4_2_energy,
          shw_sp_br4_2_ratio1_45,
          shw_sp_br4_2_ratio1_35,
          shw_sp_br4_2_ratio1_25,
          shw_sp_br4_2_ratio1_15,
          shw_sp_br4_2_iso_angle,
          shw_sp_br4_2_iso_angle1,
          shw_sp_br4_2_angle
  };

  nsm::NuSelectionBDT::SPHighEoverlap _SPHighEoverlap_init = {
          shw_sp_hol_flag,
          shw_sp_hol_1_flag,
          shw_sp_hol_1_n_valid_tracks,
          shw_sp_hol_1_min_angle,
          shw_sp_hol_1_energy,
          shw_sp_hol_1_flag_all_shower,
          shw_sp_hol_1_min_length,
          shw_sp_hol_2_flag,
          shw_sp_hol_2_min_angle,
          shw_sp_hol_2_medium_dQ_dx,
          shw_sp_hol_2_ncount,
          shw_sp_hol_2_energy
  };
  nsm::NuSelectionBDT::SPLowEoverlap _SPLowEoverlap_init = {
          shw_sp_lol_flag,
          shw_sp_lol_1_v_flag,
          shw_sp_lol_1_v_energy,
          shw_sp_lol_1_v_vtx_n_segs,
          shw_sp_lol_1_v_nseg,
          shw_sp_lol_1_v_angle,
          shw_sp_lol_2_v_flag,
          shw_sp_lol_2_v_length,
          shw_sp_lol_2_v_angle,
          shw_sp_lol_2_v_type,
          shw_sp_lol_2_v_vtx_n_segs,
          shw_sp_lol_2_v_energy,
          shw_sp_lol_2_v_shower_main_length,
          shw_sp_lol_2_v_flag_dir_weak,
          shw_sp_lol_3_flag,
          shw_sp_lol_3_angle_beam,
          shw_sp_lol_3_n_valid_tracks,
          shw_sp_lol_3_min_angle,
          shw_sp_lol_3_vtx_n_segs,
          shw_sp_lol_3_energy,
          shw_sp_lol_3_shower_main_length,
          shw_sp_lol_3_n_out,
          shw_sp_lol_3_n_sum
  };


  nsm::NuSelectionBDT::CosmicTagger _CosmicTagger_init = {
          cosmic_filled,
          cosmic_flag,
          cosmic_n_solid_tracks,
          cosmic_energy_main_showers,
          cosmic_energy_direct_showers,
          cosmic_energy_indirect_showers,
          cosmic_n_direct_showers,
          cosmic_n_indirect_showers,
          cosmic_n_main_showers
  };
  nsm::NuSelectionBDT::GapID _GapID_init = {
          gap_filled,
          gap_flag,
          gap_flag_prolong_u,
          gap_flag_prolong_v,
          gap_flag_prolong_w,
          gap_flag_parallel,
          gap_n_points,
          gap_n_bad,
          gap_energy,
          gap_num_valid_tracks,
          gap_flag_single_shower
  };
  nsm::NuSelectionBDT::MipCheck _MipCheck_init = {
          mip_quality_filled,
          mip_quality_flag,
          mip_quality_energy,
          mip_quality_overlap,
          mip_quality_n_showers,
          mip_quality_n_tracks,
          mip_quality_flag_inside_pi0,
          mip_quality_n_pi0_showers,
          mip_quality_shortest_length,
          mip_quality_acc_length,
          mip_quality_shortest_angle,
          mip_quality_flag_proton
  };
  nsm::NuSelectionBDT::MipID1 _MipID1_init = {
          mip_filled,
          mip_flag,
          mip_energy,
          mip_n_end_reduction,
          mip_n_first_mip,
          mip_n_first_non_mip,
          mip_n_first_non_mip_1,
          mip_n_first_non_mip_2,
          mip_vec_dQ_dx_0,
          mip_vec_dQ_dx_1,
          mip_max_dQ_dx_sample,
          mip_n_below_threshold,
          mip_n_below_zero,
          mip_n_lowest,
          mip_n_highest,
          mip_lowest_dQ_dx,
          mip_highest_dQ_dx,
          mip_medium_dQ_dx,
          mip_stem_length,
          mip_length_main,
          mip_length_total,
          mip_angle_beam,
          mip_iso_angle,
          mip_n_vertex,
          mip_n_good_tracks,
          mip_E_indirect_max_energy,
          mip_flag_all_above,
          mip_min_dQ_dx_5,
          mip_n_other_vertex,
          mip_n_stem_size,
          mip_flag_stem_trajectory,
          mip_min_dis
  };
  nsm::NuSelectionBDT::MipID2 _MipID2_init = {
          mip_vec_dQ_dx_2,
          mip_vec_dQ_dx_3,
          mip_vec_dQ_dx_4,
          mip_vec_dQ_dx_5,
          mip_vec_dQ_dx_6,
          mip_vec_dQ_dx_7,
          mip_vec_dQ_dx_8,
          mip_vec_dQ_dx_9,
          mip_vec_dQ_dx_10,
          mip_vec_dQ_dx_11,
          mip_vec_dQ_dx_12,
          mip_vec_dQ_dx_13,
          mip_vec_dQ_dx_14,
          mip_vec_dQ_dx_15,
          mip_vec_dQ_dx_16,
          mip_vec_dQ_dx_17,
          mip_vec_dQ_dx_18,
          mip_vec_dQ_dx_19
  };
  nsm::NuSelectionBDT::Pi0Tagger1 _Pi0Tagger1_init = {
          pio_filled,
          pio_flag,
          pio_mip_id,
          pio_flag_pio,
          pio_1_flag,
          pio_1_mass,
          pio_1_pio_type,
          pio_1_energy_1,
          pio_1_energy_2,
          pio_1_dis_1,
          pio_1_dis_2,
          pio_2_v_flag,
          pio_2_v_dis2,
          pio_2_v_angle2,
          pio_2_v_acc_length
  };
  nsm::NuSelectionBDT::Pi0Tagger2 _Pi0Tagger2_init = {
          sig_flag,
          sig_1_v_flag,
          sig_1_v_angle,
          sig_1_v_flag_single_shower,
          sig_1_v_energy,
          sig_1_v_energy_1,
          sig_2_v_flag,
          sig_2_v_energy,
          sig_2_v_shower_angle,
          sig_2_v_flag_single_shower,
          sig_2_v_medium_dQ_dx,
          sig_2_v_start_dQ_dx
  };
  nsm::NuSelectionBDT::MultiGamma1 _MultiGamma1_init = {
          mgo_flag,
          mgo_energy,
          mgo_max_energy,
          mgo_total_energy,
          mgo_n_showers,
          mgo_max_energy_1,
          mgo_max_energy_2,
          mgo_total_other_energy,
          mgo_n_total_showers,
          mgo_total_other_energy_1
  };
  nsm::NuSelectionBDT::MultiGamma2 _MultiGamma2_init = {
          mgt_flag,
          mgt_flag_single_shower,
          mgt_max_energy,
          mgt_energy,
          mgt_total_other_energy,
          mgt_max_energy_1,
          mgt_e_indirect_max_energy,
          mgt_e_direct_max_energy,
          mgt_n_direct_showers,
          mgt_e_direct_total_energy,
          mgt_flag_indirect_max_pio,
          mgt_e_indirect_total_energy
  };
  nsm::NuSelectionBDT::SingleGamma1 _SingleGamma1_init = {
          stw_flag,
          stw_1_flag,
          stw_1_energy,
          stw_1_dis,
          stw_1_dQ_dx,
          stw_1_flag_single_shower,
          stw_1_n_pi0,
          stw_1_num_valid_tracks,
          stw_2_v_flag,
          stw_2_v_medium_dQ_dx,
          stw_2_v_energy,
          stw_2_v_angle,
          stw_2_v_dir_length,
          stw_2_v_max_dQ_dx,
          stw_3_v_flag,
          stw_3_v_angle,
          stw_3_v_dir_length,
          stw_3_v_energy,
          stw_3_v_medium_dQ_dx,
          stw_4_v_flag,
          stw_4_v_angle,
          stw_4_v_dis,
          stw_4_v_energy
  };
  nsm::NuSelectionBDT::SingleGamma2 _SingleGamma2_init = {
          spt_flag,
          spt_flag_single_shower,
          spt_energy,
          spt_shower_main_length,
          spt_shower_total_length,
          spt_angle_beam,
          spt_angle_vertical,
          spt_max_dQ_dx,
          spt_angle_beam_1,
          spt_angle_drift,
          spt_angle_drift_1,
          spt_num_valid_tracks,
          spt_n_vtx_segs,
          spt_max_length
  };
  nsm::NuSelectionBDT::StemLen _StemLen_init = {
          stem_len_flag,
          stem_len_energy,
          stem_len_length,
          stem_len_flag_avoid_muon_check,
          stem_len_num_daughters,
          stem_len_daughter_length
  };
  nsm::NuSelectionBDT::LowEMichel _LowEMichel_init = {
          lem_flag,
          lem_shower_total_length,
          lem_shower_main_length,
          lem_n_3seg,
          lem_e_charge,
          lem_e_dQdx,
          lem_shower_num_segs,
          lem_shower_num_main_segs
  };
  nsm::NuSelectionBDT::BrokenMuon _BrokenMuon_init = {
          brm_flag,
          brm_n_mu_segs,
          brm_Ep,
          brm_energy,
          brm_acc_length,
          brm_shower_total_length,
          brm_connected_length,
          brm_n_size,
          brm_acc_direct_length,
          brm_n_shower_main_segs,
          brm_n_mu_main
  };
  nsm::NuSelectionBDT::MuEnergy _MuEnergy_init = {
          cme_flag,
          cme_mu_energy,
          cme_energy,
          cme_mu_length,
          cme_length,
          cme_angle_beam
  };
  nsm::NuSelectionBDT::ShowerAngle _ShowerAngle_init = {
          anc_flag,
          anc_energy,
          anc_angle,
          anc_max_angle,
          anc_max_length,
          anc_acc_forward_length,
          anc_acc_backward_length,
          anc_acc_forward_length1,
          anc_shower_main_length,
          anc_shower_total_length,
          anc_flag_main_outside
  };
  nsm::NuSelectionBDT::BadStem _BadStem_init = {
          stem_dir_filled,
          stem_dir_flag,
          stem_dir_flag_single_shower,
          stem_dir_angle,
          stem_dir_energy,
          stem_dir_angle1,
          stem_dir_angle2,
          stem_dir_angle3,
          stem_dir_ratio
  };
  nsm::NuSelectionBDT::VtxInShw _VtxInShw_init = {
          vis_flag,
          vis_1_filled,
          vis_1_flag,
          vis_1_n_vtx_segs,
          vis_1_energy,
          vis_1_num_good_tracks,
          vis_1_max_angle,
          vis_1_max_shower_angle,
          vis_1_tmp_length1,
          vis_1_tmp_length2,
          vis_1_particle_type,
          vis_2_filled,
          vis_2_flag,
          vis_2_n_vtx_segs,
          vis_2_min_angle,
          vis_2_min_weak_track,
          vis_2_angle_beam,
          vis_2_min_angle1,
          vis_2_iso_angle1,
          vis_2_min_medium_dQ_dx,
          vis_2_min_length,
          vis_2_sg_length,
          vis_2_max_angle,
          vis_2_max_weak_track
  };
  nsm::NuSelectionBDT::BadReco1 _BadReco1_init = {
          br_filled,
          br1_flag,
          br1_1_flag,
          br1_1_shower_type,
          br1_1_vtx_n_segs,
          br1_1_energy,
          br1_1_n_segs,
          br1_1_flag_sg_topology,
          br1_1_flag_sg_trajectory,
          br1_1_sg_length,
          br1_2_flag,
          br1_2_energy,
          br1_2_n_connected,
          br1_2_max_length,
          br1_2_n_connected_1,
          br1_2_vtx_n_segs,
          br1_2_n_shower_segs,
          br1_2_max_length_ratio,
          br1_2_shower_length,
          br1_3_flag,
          br1_3_energy,
          br1_3_n_connected_p,
          br1_3_max_length_p,
          br1_3_n_shower_segs,
          br1_3_flag_sg_topology,
          br1_3_flag_sg_trajectory,
          br1_3_n_shower_main_segs,
          br1_3_sg_length
  };
  nsm::NuSelectionBDT::BadReco2 _BadReco2_init = {
          br_filled,
          br2_flag,
          br2_flag_single_shower,
          br2_num_valid_tracks,
          br2_energy,
          br2_angle1,
          br2_angle2,
          br2_angle,
          br2_angle3,
          br2_n_shower_main_segs,
          br2_max_angle,
          br2_sg_length,
          br2_flag_sg_trajectory
  };
  nsm::NuSelectionBDT::BadReco3 _BadReco3_init = {
          br_filled,
          br3_flag,
          br3_1_flag,
          br3_1_energy,
          br3_1_n_shower_segments,
          br3_1_sg_flag_trajectory,
          br3_1_sg_direct_length,
          br3_1_sg_length,
          br3_1_total_main_length,
          br3_1_total_length,
          br3_1_iso_angle,
          br3_1_sg_flag_topology,
          br3_2_flag,
          br3_2_n_ele,
          br3_2_n_other,
          br3_2_energy,
          br3_2_total_main_length,
          br3_2_total_length,
          br3_2_other_fid,
          br3_3_v_flag,
          br3_3_v_energy,
          br3_3_v_angle,
          br3_3_v_dir_length,
          br3_3_v_length,
          br3_4_flag,
          br3_4_acc_length,
          br3_4_total_length,
          br3_4_energy,
          br3_5_v_flag,
          br3_5_v_dir_length,
          br3_5_v_total_length,
          br3_5_v_flag_avoid_muon_check,
          br3_5_v_n_seg,
          br3_5_v_angle,
          br3_5_v_sg_length,
          br3_5_v_energy,
          br3_5_v_n_main_segs,
          br3_5_v_n_segs,
          br3_5_v_shower_main_length,
          br3_5_v_shower_total_length,
          br3_6_v_flag,
          br3_6_v_angle,
          br3_6_v_angle1,
          br3_6_v_flag_shower_trajectory,
          br3_6_v_direct_length,
          br3_6_v_length,
          br3_6_v_n_other_vtx_segs,
          br3_6_v_energy,
          br3_7_flag,
          br3_7_energy,
          br3_7_min_angle,
          br3_7_sg_length,
          br3_7_shower_main_length,
          br3_8_flag,
          br3_8_max_dQ_dx,
          br3_8_energy,
          br3_8_n_main_segs,
          br3_8_shower_main_length,
          br3_8_shower_length
  };
  nsm::NuSelectionBDT::BadReco4 _BadReco4_init = {
          br_filled,
          br4_flag,
          br4_1_flag,
          br4_1_shower_main_length,
          br4_1_shower_total_length,
          br4_1_min_dis,
          br4_1_energy,
          br4_1_flag_avoid_muon_check,
          br4_1_n_vtx_segs,
          br4_1_n_main_segs,
          br4_2_flag,
          br4_2_ratio_45,
          br4_2_ratio_35,
          br4_2_ratio_25,
          br4_2_ratio_15,
          br4_2_energy,
          br4_2_ratio1_45,
          br4_2_ratio1_35,
          br4_2_ratio1_25,
          br4_2_ratio1_15,
          br4_2_iso_angle,
          br4_2_iso_angle1,
          br4_2_angle
  };
  nsm::NuSelectionBDT::TrackOverCluster _TrackOverCluster_init = {
          tro_flag,
          tro_1_v_flag,
          tro_1_v_particle_type,
          tro_1_v_flag_dir_weak,
          tro_1_v_min_dis,
          tro_1_v_sg1_length,
          tro_1_v_shower_main_length,
          tro_1_v_max_n_vtx_segs,
          tro_1_v_tmp_length,
          tro_1_v_medium_dQ_dx,
          tro_1_v_dQ_dx_cut,
          tro_1_v_flag_shower_topology,
          tro_2_v_flag,
          tro_2_v_energy,
          tro_2_v_stem_length,
          tro_2_v_iso_angle,
          tro_2_v_max_length,
          tro_2_v_angle,
          tro_3_flag,
          tro_3_stem_length,
          tro_3_n_muon_segs,
          tro_3_energy,
          tro_4_v_flag,
          tro_4_v_dir2_mag,
          tro_4_v_angle,
          tro_4_v_angle1,
          tro_4_v_angle2,
          tro_4_v_length,
          tro_4_v_length1,
          tro_4_v_medium_dQ_dx,
          tro_4_v_end_dQ_dx,
          tro_4_v_energy,
          tro_4_v_shower_main_length,
          tro_4_v_flag_shower_trajectory,
          tro_5_v_flag,
          tro_5_v_max_angle,
          tro_5_v_min_angle,
          tro_5_v_max_length,
          tro_5_v_iso_angle,
          tro_5_v_n_vtx_segs,
          tro_5_v_min_count,
          tro_5_v_max_count,
          tro_5_v_energy
  };
  nsm::NuSelectionBDT::HighEoverlap _HighEoverlap_init = {
          hol_flag,
          hol_1_flag,
          hol_1_n_valid_tracks,
          hol_1_min_angle,
          hol_1_energy,
          hol_1_flag_all_shower,
          hol_1_min_length,
          hol_2_flag,
          hol_2_min_angle,
          hol_2_medium_dQ_dx,
          hol_2_ncount,
          hol_2_energy
  };
  nsm::NuSelectionBDT::LowEoverlap _LowEoverlap_init = {
          lol_flag,
          lol_1_v_flag,
          lol_1_v_energy,
          lol_1_v_vtx_n_segs,
          lol_1_v_nseg,
          lol_1_v_angle,
          lol_2_v_flag,
          lol_2_v_length,
          lol_2_v_angle,
          lol_2_v_type,
          lol_2_v_vtx_n_segs,
          lol_2_v_energy,
          lol_2_v_shower_main_length,
          lol_2_v_flag_dir_weak,
          lol_3_flag,
          lol_3_angle_beam,
          lol_3_n_valid_tracks,
          lol_3_min_angle,
          lol_3_vtx_n_segs,
          lol_3_energy,
          lol_3_shower_main_length,
          lol_3_n_out,
          lol_3_n_sum
  };
  nsm::NuSelectionBDT::MajorCosmicTagger _MajorCosmicTagger_init = {
	  cosmict_flag_1, // fiducial volume vertex
	  cosmict_flag_2,  // single muon
	  cosmict_flag_3,  // single muon (long)
	  cosmict_flag_4,  // kinematics muon
	  cosmict_flag_5, // kinematics muon (long)
	  cosmict_flag_6, // special ...
	  cosmict_flag_7,  // muon+ michel
	  cosmict_flag_8,  // muon + michel + special
	  cosmict_flag_9,  // this tagger is relevant for nueCC, see "cosmic tagger ones, one case of cosmics ..." (frist one ...)
	  cosmict_flag_10,  // front upstream (dirt)
	  cosmict_flag,
	  cosmict_2_filled,
	  cosmict_2_particle_type,
	  cosmict_2_n_muon_tracks,
	  cosmict_2_total_shower_length,
	  cosmict_2_flag_inside,
	  cosmict_2_angle_beam,
	  cosmict_2_flag_dir_weak,
	  cosmict_2_dQ_dx_end,
	  cosmict_2_dQ_dx_front,
	  cosmict_2_theta,
	  cosmict_2_phi,
	  cosmict_2_valid_tracks,
	  cosmict_3_filled,
	  cosmict_3_flag_inside,
	  cosmict_3_angle_beam,
	  cosmict_3_flag_dir_weak,
	  cosmict_3_dQ_dx_end,
	  cosmict_3_dQ_dx_front,
	  cosmict_3_theta,
	  cosmict_3_phi,
	  cosmict_3_valid_tracks,
	  cosmict_4_filled,
	  cosmict_4_flag_inside,
	  cosmict_4_angle_beam,
	  cosmict_4_connected_showers,  // need to be careful about the nueCC ...
	  cosmict_5_filled,
	  cosmict_5_flag_inside,
	  cosmict_5_angle_beam,
	  cosmict_5_connected_showers,
	  cosmict_6_filled,
	  cosmict_6_flag_dir_weak,
	  cosmict_6_flag_inside,
	  cosmict_6_angle,
	  cosmict_7_filled,
	  cosmict_7_flag_sec,
	  cosmict_7_n_muon_tracks,
	  cosmict_7_total_shower_length,
	  cosmict_7_flag_inside,
	  cosmict_7_angle_beam,
	  cosmict_7_flag_dir_weak,
	  cosmict_7_dQ_dx_end,
	  cosmict_7_dQ_dx_front,
	  cosmict_7_theta,
	  cosmict_7_phi,
	  cosmict_8_filled,
	  cosmict_8_flag_out,
	  cosmict_8_muon_length,
	  cosmict_8_acc_length,
	  cosmict_10_flag_inside,
	  cosmict_10_vtx_z,
	  cosmict_10_flag_shower,
	  cosmict_10_flag_dir_weak,
	  cosmict_10_angle_beam,
	  cosmict_10_length
  };
  nsm::NuSelectionBDT::NumuCCTagger _NumuCCTagger_init = {
	  numu_cc_flag,
	  numu_cc_flag_1,
	  numu_cc_1_particle_type,
	  numu_cc_1_length,
	  numu_cc_1_medium_dQ_dx,
	  numu_cc_1_dQ_dx_cut,
	  numu_cc_1_direct_length,
	  numu_cc_1_n_daughter_tracks,
	  numu_cc_1_n_daughter_all,
	  numu_cc_flag_2,
	  numu_cc_2_length,
	  numu_cc_2_total_length,
	  numu_cc_2_n_daughter_tracks,
	  numu_cc_2_n_daughter_all,
	  numu_cc_flag_3,
	  numu_cc_3_particle_type,
	  numu_cc_3_max_length,
	  numu_cc_3_acc_track_length,
	  numu_cc_3_max_length_all,
	  numu_cc_3_max_muon_length,
	  numu_cc_3_n_daughter_tracks,
	  numu_cc_3_n_daughter_all,
  };
  nsm::NuSelectionBDT::BDTscores _BDTscores_init = {
	  cosmict_2_4_score,
	  cosmict_3_5_score,
	  cosmict_6_score,
	  cosmict_7_score,
	  cosmict_8_score,
	  cosmict_10_score,
	  numu_1_score,
	  numu_2_score,
	  numu_3_score,
	  cosmict_score,
	  numu_score,
	  mipid_score,
	  gap_score,
	  hol_lol_score,
	  cme_anc_score,
	  mgo_mgt_score,
	  br1_score,
	  br3_score,
	  br3_3_score,
	  br3_5_score,
	  br3_6_score,
	  stemdir_br2_score,
	  trimuon_score,
	  br4_tro_score,
	  mipquality_score,
	  pio_1_score,
	  pio_2_score,
	  stw_spt_score,
	  vis_1_score,
	  vis_2_score,
	  stw_2_score,
	  stw_3_score,
	  stw_4_score,
	  sig_1_score,
	  sig_2_score,
	  lol_1_score,
	  lol_2_score,
	  tro_1_score,
	  tro_2_score,
	  tro_4_score,
	  tro_5_score,
	  nue_score
  };

  std::cout<<"T_tagger size: "<<tree3->GetEntries()<<std::endl;
  if(tree3->GetEntries()==0){ std::cout<<"Empty T_tagger"<<std::endl;}
  else{

    // set
    if(f_ssmBDT){ 
      nsmbdt.Setstkdar(_stkdar_init); 
    }
    std::cout<<"numu_cc_flag "<<numu_cc_flag<<std::endl;
    bool flag_fill_bdt = true;
    if(numu_cc_flag<0){ flag_fill_bdt = false; }
    if (flag_fill_bdt) {
      nsmbdt.SetSPID(_SPID_init);
      nsmbdt.SetSPSHWID1(_SPSHWID1_init);
      nsmbdt.SetSPSHWID2(_SPSHWID2_init);
      nsmbdt.SetSPPi0Tagger1(_SPPi0Tagger1_init);
      nsmbdt.SetSPLowEMichel(_SPLowEMichel_init);
      nsmbdt.SetSPBadReco1(_SPBadReco1_init);
      nsmbdt.SetSPBadReco2(_SPBadReco2_init);
      nsmbdt.SetSPBadReco3(_SPBadReco3_init);
      nsmbdt.SetSPBadReco4(_SPBadReco4_init);
      nsmbdt.SetSPHighEoverlap(_SPHighEoverlap_init);
      nsmbdt.SetSPLowEoverlap(_SPLowEoverlap_init);
      nsmbdt.SetCosmicTagger(_CosmicTagger_init);
      nsmbdt.SetGapID(_GapID_init);
      nsmbdt.SetMipCheck(_MipCheck_init);
      nsmbdt.SetMipID1(_MipID1_init);
      nsmbdt.SetMipID2(_MipID2_init);
      nsmbdt.SetPi0Tagger1(_Pi0Tagger1_init);
      nsmbdt.SetPi0Tagger2(_Pi0Tagger2_init);
      nsmbdt.SetMultiGamma1(_MultiGamma1_init);
      nsmbdt.SetMultiGamma2(_MultiGamma2_init);
      nsmbdt.SetSingleGamma1(_SingleGamma1_init);
      nsmbdt.SetSingleGamma2(_SingleGamma2_init);
      nsmbdt.SetStemLen(_StemLen_init);
      nsmbdt.SetLowEMichel(_LowEMichel_init);
      nsmbdt.SetBrokenMuon(_BrokenMuon_init);
      nsmbdt.SetMuEnergy(_MuEnergy_init);
      nsmbdt.SetShowerAngle(_ShowerAngle_init);
      nsmbdt.SetBadStem(_BadStem_init);
      nsmbdt.SetVtxInShw(_VtxInShw_init);
      nsmbdt.SetBadReco1(_BadReco1_init);
      nsmbdt.SetBadReco2(_BadReco2_init);
      nsmbdt.SetBadReco3(_BadReco3_init);
      nsmbdt.SetBadReco4(_BadReco4_init);
      nsmbdt.SetTrackOverCluster(_TrackOverCluster_init);
      nsmbdt.SetHighEoverlap(_HighEoverlap_init);
      nsmbdt.SetLowEoverlap(_LowEoverlap_init);
      nsmbdt.SetMajorCosmicTagger(_MajorCosmicTagger_init);
      nsmbdt.SetNumuCCTagger(_NumuCCTagger_init);
      nsmbdt.SetBDTscores(_BDTscores_init);
    }

    outputBDTvars->push_back(nsmbdt);
  }
  }//end if(tree3)
  else {
    mf::LogError("WireCellPF") <<"TTree "<< fInput_tree3 <<" not found in file " << fInput <<"\n";
  }

  e.put(std::move(outputBDTvars));

  }//else
}

if(f_KINEport){
  auto outputKINEvars = std::make_unique< std::vector<nsm::NuSelectionKINE> >();
  if(badinput){
	e.put(std::move(outputKINEvars));
	//return;
  }
  else{
  nsm::NuSelectionKINE nsmkine;

  TTree *tree4 = (TTree*)fin->Get(fInput_tree4.c_str());
  if(tree4){

	  float kine_reco_Enu; // kinetic energy  + additional energy ...
	  float kine_reco_add_energy;  // mass, binding energy ...
	  std::vector<float> *kine_energy_particle = new std::vector<float>;  // energy of each particle
	  std::vector<int> *kine_energy_info = new std::vector<int>; // what kind of energy reconstruction?
	  std::vector<int> *kine_particle_type = new std::vector<int>;
	  std::vector<int> *kine_energy_included = new std::vector<int>; // included in the neutrino energy calculation?
	  float kine_pio_mass; // mass
	  int kine_pio_flag; // 0 not filled, 1, with vertex: CCpio, 2 without vertex: NCpi0
	  float kine_pio_vtx_dis;
	  float kine_pio_energy_1;
	  float kine_pio_theta_1;
	  float kine_pio_phi_1;
	  float kine_pio_dis_1;
	  float kine_pio_energy_2;
	  float kine_pio_theta_2;
	  float kine_pio_phi_2;
	  float kine_pio_dis_2;
	  float kine_pio_angle;

	  tree4->SetBranchAddress("kine_reco_Enu", &kine_reco_Enu); // kinetic energy  + additional energy ...
	  tree4->SetBranchAddress("kine_reco_add_energy", &kine_reco_add_energy);  // mass, binding energy ...
	  tree4->SetBranchAddress("kine_energy_particle", &kine_energy_particle);  // energy of each particle
	  tree4->SetBranchAddress("kine_energy_info", &kine_energy_info); // what kind of energy reconstruction?
	  tree4->SetBranchAddress("kine_particle_type", &kine_particle_type);
	  tree4->SetBranchAddress("kine_energy_included", &kine_energy_included); // included in the neutrino energy calculation?
	  tree4->SetBranchAddress("kine_pio_mass", &kine_pio_mass); // mass
	  tree4->SetBranchAddress("kine_pio_flag", &kine_pio_flag); // 0 not filled, 1, with vertex: CCpio, 2 without vertex: NCpi0
	  tree4->SetBranchAddress("kine_pio_vtx_dis", &kine_pio_vtx_dis);
	  tree4->SetBranchAddress("kine_pio_energy_1", &kine_pio_energy_1);
	  tree4->SetBranchAddress("kine_pio_theta_1", &kine_pio_theta_1);
	  tree4->SetBranchAddress("kine_pio_phi_1", &kine_pio_phi_1);
	  tree4->SetBranchAddress("kine_pio_dis_1", &kine_pio_dis_1);
	  tree4->SetBranchAddress("kine_pio_energy_2", &kine_pio_energy_2);
	  tree4->SetBranchAddress("kine_pio_theta_2", &kine_pio_theta_2);
	  tree4->SetBranchAddress("kine_pio_phi_2", &kine_pio_phi_2);
	  tree4->SetBranchAddress("kine_pio_dis_2", &kine_pio_dis_2);
	  tree4->SetBranchAddress("kine_pio_angle", &kine_pio_angle);

	  //read and port
  	  tree4->GetEntry(0);
	  nsm::NuSelectionKINE::KineInfo _KineInfo_init = {
		  kine_reco_Enu,
		  kine_reco_add_energy,
		  kine_energy_particle,
		  kine_energy_info,
		  kine_particle_type,
		  kine_energy_included,
		  kine_pio_mass,
		  kine_pio_flag,
		  kine_pio_vtx_dis,
		  kine_pio_energy_1,
		  kine_pio_theta_1,
		  kine_pio_phi_1,
		  kine_pio_dis_1,
		  kine_pio_energy_2,
		  kine_pio_theta_2,
		  kine_pio_phi_2,
		  kine_pio_dis_2,
		  kine_pio_angle
	  };
          bool flag_fill_kine = true;
          if(numu_cc_flag<0 && f_BDTport){ flag_fill_kine = false; }
          if (flag_fill_kine) { nsmkine.SetKineInfo(_KineInfo_init); }
	  outputKINEvars->push_back(nsmkine);

  }
  else{
    mf::LogError("WireCellPF") <<"TTree "<< fInput_tree4 <<" not found in file " << fInput <<"\n";
  }

  e.put(std::move(outputKINEvars));

  }
}


  if(fin) fin->Close();
  return;

}

void nsm::WireCellPF::beginJob()
{
  // Implementation of optional member function here.
  mf::LogInfo("WireCellPF") << "Job begins: "
                            << "Number of PF ported: " << NumberOfPF <<"\n";
}

void nsm::WireCellPF::endJob()
{
  // Implementation of optional member function here.
  mf::LogInfo("WireCellPF") << "Job ends: "
                            << "Number of PF ported: " << NumberOfPF <<"\n";
}

DEFINE_ART_MODULE(nsm::WireCellPF)
