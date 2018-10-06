////////////////////////////////////////////////////////////////////////
// Class:       ShowerTruthStudies
// Plugin Type: analyzer (art v2_11_03)
// File:        ShowerTruthStudies_module.cc
//
// Generated at Sat Oct  6 14:03:01 2018 by David Caratelli using cetskelgen
// from cetlib version v3_03_01.
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

#include "TTree.h"
#include "TVector3.h"

#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardata/Utilities/AssociationUtil.h"


class ShowerTruthStudies;


class ShowerTruthStudies : public art::EDAnalyzer {
public:
  explicit ShowerTruthStudies(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ShowerTruthStudies(ShowerTruthStudies const &) = delete;
  ShowerTruthStudies(ShowerTruthStudies &&) = delete;
  ShowerTruthStudies & operator = (ShowerTruthStudies const &) = delete;
  ShowerTruthStudies & operator = (ShowerTruthStudies &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  void Reset();

// Match each MCshower to its best matching RCshower
  std::vector<int> Match(const std::vector<sim::MCShower>& mc_shower_v,
			 const std::vector<recob::Shower>& shr_v);
// Match each RCshower to its best matching MCshower
  std::vector<int> Match(const std::vector<recob::Shower>& shr_v,
			 const std::vector<sim::MCShower>& mc_shower_v);
			 

  void ClearMCRC();

  void SetTTree();

  // shower-by-shower TTree
  TTree* _mcshr_tree;
  TTree* _rcshr_tree;

  // variables common to both ttrees
  int _run, _sub, _evt;

  int _n_reco_showers;

  double _nu_e, _pi0_e;
  double _emc;

  double _mc_vtx_x, _mc_vtx_y, _mc_vtx_z;
  double _rc_vtx_x, _rc_vtx_y, _rc_vtx_z;

  // shower-by-shower variables [reco]
  double _rc_shr_x,  _rc_shr_y,  _rc_shr_z;
  double _rc_shr_px, _rc_shr_py, _rc_shr_pz;
  double _rc_shr_e;
  double _rc_shr_dedx;
  double _rcradlen;
  // shower-by-shower variables [truth]
  double _mc_shr_x,  _mc_shr_y,  _mc_shr_z;
  double _mc_shr_px, _mc_shr_py, _mc_shr_pz;
  double _mc_shr_e, _mc_shr_edep;
  double _mc_shr_dedx;
  double _mcradlen;
  // shower-by-shower variables, comparisons
  double _angle;
  double _strtdiff;
  double _dwallmin;

  std::string fShrProducer;


};


ShowerTruthStudies::ShowerTruthStudies(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  fShrProducer = p.get<std::string>("ShrProducer");
  SetTTree();
}

void ShowerTruthStudies::analyze(art::Event const & e)
{

  _run = e.run();
  _sub = e.subRun();
  _evt = e.event();

  // load input showers
  auto const& shr_h = e.getValidHandle<std::vector<recob::Shower>>(fShrProducer);
  // load mcshowers
  auto const& mcs_h = e.getValidHandle<std::vector<sim::MCShower>>("mcreco");
  // load MCTruth
  auto const& mct_h = e.getValidHandle<std::vector<simb::MCTruth> >("generator");

  // Store MC and RC showers in vectors
  _n_reco_showers = shr_h->size();
  std::vector<recob::Shower> reco_shower_v;
  for (size_t i=0; i < shr_h->size(); i++)
    reco_shower_v.push_back( shr_h->at(i) );
  std::vector<sim::MCShower> mc_shower_v;

  auto mct = mct_h->at(0);

  auto neutrino = mct.GetNeutrino();

  auto nu     = neutrino.Nu();
  //auto lepton = neutrino.Lepton();

  _mc_vtx_x = nu.Trajectory().X( nu.Trajectory().size() - 1 );
  _mc_vtx_y = nu.Trajectory().Y( nu.Trajectory().size() - 1 );
  _mc_vtx_z = nu.Trajectory().Z( nu.Trajectory().size() - 1 );

  // search all MCShowers produced @ the vertex
  for (size_t i=0; i < mcs_h->size(); i++){
    auto const& mcs = mcs_h->at(i);
    // distance from vertex                                                                
    double x = mcs.Start().X();
    double y = mcs.Start().Y();
    double z = mcs.Start().Z();
    double d = sqrt( ( (_mc_vtx_x - x) * (_mc_vtx_x - x) ) +
		     ( (_mc_vtx_y - y) * (_mc_vtx_y - y) ) +
		     ( (_mc_vtx_z - z) * (_mc_vtx_z - z) ) );
    if ( d < 0.01 )
      mc_shower_v.push_back( mcs );
  }// for all MCShowers                                                               
    
  
  // MC <-> RC matching
  auto MCRCmatch_v = Match(mc_shower_v, reco_shower_v);
  
  for (size_t mcidx = 0; mcidx < MCRCmatch_v.size(); mcidx++) {
    
    ClearMCRC();
    
    auto mcshr = mc_shower_v.at(mcidx);
    
    _mc_shr_e  = mcshr.Start().E();
    _mc_shr_edep  = mcshr.DetProfile().E();
    _mc_shr_x  = mcshr.DetProfile().X();
    _mc_shr_y  = mcshr.DetProfile().Y();
    _mc_shr_z  = mcshr.DetProfile().Z();
    
    double mommc = mcshr.Start().Momentum().Vect().Mag();
    _mc_shr_px = mcshr.Start().Px() / mommc;
    _mc_shr_py = mcshr.Start().Py() / mommc;
    _mc_shr_pz = mcshr.Start().Pz() / mommc;
    
    _mcradlen = sqrt( ( (_mc_shr_x - _mc_vtx_x) * (_mc_shr_x - _mc_vtx_x) ) +
		      ( (_mc_shr_y - _mc_vtx_y) * (_mc_shr_y - _mc_vtx_y) ) +
		      ( (_mc_shr_z - _mc_vtx_z) * (_mc_shr_z - _mc_vtx_z) ) );
    
    if (MCRCmatch_v.at(mcidx) == -1) {
      _rc_shr_e = 0;
      _angle    = -5;
      _mcshr_tree->Fill();
      continue;
    }
    
    auto const& rcshr = shr_h->at( MCRCmatch_v.at(mcidx) );
    
    _rc_shr_e = rcshr.Energy()[2];
    _rc_shr_x = rcshr.ShowerStart().X();
    _rc_shr_y = rcshr.ShowerStart().Y();
    _rc_shr_z = rcshr.ShowerStart().Z();
    _rc_shr_dedx = rcshr.dEdx()[2];
    
    double momrc = rcshr.Direction().Mag();      
    _rc_shr_px = rcshr.Direction().X() / momrc;
    _rc_shr_py = rcshr.Direction().Y() / momrc;
    _rc_shr_pz = rcshr.Direction().Z() / momrc;
    
    _rcradlen = sqrt( ( (_rc_shr_x - _rc_vtx_x) * (_rc_shr_x - _rc_vtx_x) ) +
		      ( (_rc_shr_y - _rc_vtx_y) * (_rc_shr_y - _rc_vtx_y) ) +
		      ( (_rc_shr_z - _rc_vtx_z) * (_rc_shr_z - _rc_vtx_z) ) );
    
    _angle  = rcshr.Direction().Angle( mcshr.Start().Momentum().Vect() );
    
    _emc = mcshr.Start().E();
    
    _strtdiff = (rcshr.ShowerStart() - mcshr.DetProfile().Position().Vect()).Mag();
    
    _mcshr_tree->Fill();
    
  }// loop through MC showers
  
  // RC <-> MC matching
  auto RCMCmatch_v = Match(reco_shower_v, mc_shower_v);
  
  for (size_t rcidx = 0; rcidx < RCMCmatch_v.size(); rcidx++) {
    
    ClearMCRC();
    
    auto const& rcshr = reco_shower_v.at(rcidx);
    
    _rc_shr_e = rcshr.Energy()[2];
    _rc_shr_x = rcshr.ShowerStart().X();
    _rc_shr_y = rcshr.ShowerStart().Y();
    _rc_shr_z = rcshr.ShowerStart().Z();
    _rc_shr_dedx = rcshr.dEdx()[2];

    double mom = rcshr.Direction().Mag();    
    _rc_shr_px = rcshr.Direction().X() / mom;
    _rc_shr_py = rcshr.Direction().Y() / mom;
    _rc_shr_pz = rcshr.Direction().Z() / mom;
    
    _rcradlen = sqrt( ( (_rc_shr_x - _rc_vtx_x) * (_rc_shr_x - _rc_vtx_x) ) +
		      ( (_rc_shr_y - _rc_vtx_y) * (_rc_shr_y - _rc_vtx_y) ) +
		      ( (_rc_shr_z - _rc_vtx_z) * (_rc_shr_z - _rc_vtx_z) ) );

    
    if (RCMCmatch_v.at(rcidx) == -1) {
      _mc_shr_e = 0;
      _angle    = -5;
      _rcshr_tree->Fill();
      continue;
    }
    
    auto const& mcshr = mc_shower_v.at( RCMCmatch_v.at(rcidx) );

    _mc_shr_e  = mcshr.Start().E();
    _mc_shr_edep  = mcshr.DetProfile().E();
    _mc_shr_x  = mcshr.DetProfile().X();
    _mc_shr_y  = mcshr.DetProfile().Y();
    _mc_shr_z  = mcshr.DetProfile().Z();

    _mcradlen = sqrt( ( (_mc_shr_x - _mc_vtx_x) * (_mc_shr_x - _mc_vtx_x) ) +
		      ( (_mc_shr_y - _mc_vtx_y) * (_mc_shr_y - _mc_vtx_y) ) +
		      ( (_mc_shr_z - _mc_vtx_z) * (_mc_shr_z - _mc_vtx_z) ) );
    
    double mommc = mcshr.Start().Momentum().Vect().Mag();
    _mc_shr_px = mcshr.Start().Px() / mommc;
    _mc_shr_py = mcshr.Start().Py() / mommc;
    _mc_shr_pz = mcshr.Start().Pz() / mommc;
    
    _angle  = rcshr.Direction().Angle( mcshr.Start().Momentum().Vect() );
    
    _emc = mcshr.Start().E();
    
    _strtdiff = (rcshr.ShowerStart() - mcshr.DetProfile().Position().Vect()).Mag();
    
    _rcshr_tree->Fill();
    
  }// loop through RC showers
  
  return;
}

// Match each MCshower to its best matching RCshower
std::vector<int> ShowerTruthStudies::Match(const std::vector<sim::MCShower>& mc_shower_v,
					   const std::vector<recob::Shower>& shr_v) {


  // now match to true showers 
  std::vector<int> matched_indices;

  // find best matching RC shower for each MC shower               
  for (auto const& mcshr : mc_shower_v){

    double dotmax = -1.;
    int idxmax = -1;
    
    // loop through reco showers
    for (size_t i=0; i < shr_v.size(); i++){
      auto const& rcshr = shr_v.at(i);
      
      double dot = rcshr.Direction().Dot( mcshr.Start().Momentum().Vect() );
      dot /= mcshr.Start().Momentum().Vect().Mag();
      dot /= rcshr.Direction().Mag();
      
      if (dot > dotmax) { dotmax = dot; idxmax = i; }
      
    }// for all reco indices                                                 

    matched_indices.push_back( idxmax );
  }
  
  return matched_indices;
  
}// end of function   

// Match each RCshower to its best matching MCshower
std::vector<int> ShowerTruthStudies::Match(const std::vector<recob::Shower>& shr_v,
					   const std::vector<sim::MCShower>& mc_shower_v) {


  // now match to true showers                                                                                                                                                                                           
  std::vector<int> matched_indices;

  // find best matching MC shower for each RC shower       
  for (auto const& rcshr : shr_v) {
    
    double dotmax = -1.;
    int idxmax = -1;
    
    // loop through reco showers
    for (size_t i=0; i < mc_shower_v.size(); i++) {

      auto const& mcshr = mc_shower_v.at(i);
      
      double dot = rcshr.Direction().Dot( mcshr.Start().Momentum().Vect() );
      dot /= mcshr.Start().Momentum().Vect().Mag();
      dot /= rcshr.Direction().Mag();
      
      if (dot > dotmax) { dotmax = dot; idxmax = i; }
      
    }// for all reco indices                                                 
    
    matched_indices.push_back( idxmax );
  }
  
  return matched_indices;
  
}// end of function   

void ShowerTruthStudies::ClearMCRC() {

  _rc_shr_x=  _rc_shr_y=  _rc_shr_z= 0;
  _rc_shr_px= _rc_shr_py= _rc_shr_pz= 0;
  _rc_shr_e= 0;
  _rc_shr_dedx = 0;

  _mc_shr_x=  _mc_shr_y=  _mc_shr_z= 0;
  _mc_shr_px= _mc_shr_py= _mc_shr_pz= 0;
  _mc_shr_e= _mc_shr_edep = 0;
  _mc_shr_dedx = 0;

  return;
}

void ShowerTruthStudies::SetTTree() {

  art::ServiceHandle<art::TFileService> tfs;

  // MC shower-by-shower TTree
  _mcshr_tree = tfs->make<TTree>("_mcshr_tree","Pi0 Tree TTree");

  _mcshr_tree->Branch("_run",&_run,"run/I");
  _mcshr_tree->Branch("_sub",&_sub,"sub/I");
  _mcshr_tree->Branch("_evt",&_evt,"evt/I");

  _mcshr_tree->Branch("_n_reco_showers",&_n_reco_showers,"n_reco_showers/I");
  // vertex info                                                                                           
  _mcshr_tree->Branch("_mc_vtx_x",&_mc_vtx_x,"mc_vtx_x/D");
  _mcshr_tree->Branch("_mc_vtx_y",&_mc_vtx_y,"mc_vtx_y/D");
  _mcshr_tree->Branch("_mc_vtx_z",&_mc_vtx_z,"mc_vtx_z/D");
  _mcshr_tree->Branch("_rc_vtx_x",&_rc_vtx_x,"rc_vtx_x/D");
  _mcshr_tree->Branch("_rc_vtx_y",&_rc_vtx_y,"rc_vtx_y/D");
  _mcshr_tree->Branch("_rc_vtx_z",&_rc_vtx_z,"rc_vtx_z/D");

  // reco shower info [rec] 
  _mcshr_tree->Branch("_rc_shr_x",&_rc_shr_x,"rc_shr_x/D");
  _mcshr_tree->Branch("_rc_shr_y",&_rc_shr_y,"rc_shr_y/D");
  _mcshr_tree->Branch("_rc_shr_z",&_rc_shr_z,"rc_shr_z/D");
  _mcshr_tree->Branch("_rc_shr_e",&_rc_shr_e,"rc_shr_e/D");
  _mcshr_tree->Branch("_rc_shr_dedx",&_rc_shr_dedx,"rc_shr_dedx/D");
  _mcshr_tree->Branch("_rc_shr_px",&_rc_shr_px,"rc_shr_px/D");
  _mcshr_tree->Branch("_rc_shr_py",&_rc_shr_py,"rc_shr_py/D");
  _mcshr_tree->Branch("_rc_shr_pz",&_rc_shr_pz,"rc_shr_pz/D");
  _mcshr_tree->Branch("_rcradlen",&_rcradlen,"rcradlen/D");
  // reco shower info [mc]
  _mcshr_tree->Branch("_mc_shr_x",&_mc_shr_x,"mc_shr_x/D");
  _mcshr_tree->Branch("_mc_shr_y",&_mc_shr_y,"mc_shr_y/D");
  _mcshr_tree->Branch("_mc_shr_z",&_mc_shr_z,"mc_shr_z/D");
  _mcshr_tree->Branch("_mc_shr_e",&_mc_shr_e,"mc_shr_e/D");
  _mcshr_tree->Branch("_mc_shr_edep",&_mc_shr_edep,"mc_shr_edep/D");
  _mcshr_tree->Branch("_mc_shr_dedx",&_mc_shr_dedx,"mc_shr_dedx/D");
  _mcshr_tree->Branch("_mc_shr_px",&_mc_shr_px,"mc_shr_px/D");
  _mcshr_tree->Branch("_mc_shr_py",&_mc_shr_py,"mc_shr_py/D");
  _mcshr_tree->Branch("_mc_shr_pz",&_mc_shr_pz,"mc_shr_pz/D");
  _mcshr_tree->Branch("_mcradlen",&_mcradlen,"mcradlen/D");

  // MC -> RC shower comparisons                                          
  _mcshr_tree->Branch("_angle",&_angle,"angle/D");
  _mcshr_tree->Branch("_strtdiff",&_strtdiff,"strtdiff/D");
  _mcshr_tree->Branch("_emc",&_emc,"emc/D");

  // MC shower-by-shower TTree
  _rcshr_tree = tfs->make<TTree>("_rcshr_tree","Pi0 Tree TTree");

  _rcshr_tree->Branch("_run",&_run,"run/I");
  _rcshr_tree->Branch("_sub",&_sub,"sub/I");
  _rcshr_tree->Branch("_evt",&_evt,"evt/I");

  _rcshr_tree->Branch("_n_reco_showers",&_n_reco_showers,"n_reco_showers/I");
  // vertex info                                                                                           
  _rcshr_tree->Branch("_mc_vtx_x",&_mc_vtx_x,"mc_vtx_x/D");
  _rcshr_tree->Branch("_mc_vtx_y",&_mc_vtx_y,"mc_vtx_y/D");
  _rcshr_tree->Branch("_mc_vtx_z",&_mc_vtx_z,"mc_vtx_z/D");
  _rcshr_tree->Branch("_rc_vtx_x",&_rc_vtx_x,"rc_vtx_x/D");
  _rcshr_tree->Branch("_rc_vtx_y",&_rc_vtx_y,"rc_vtx_y/D");
  _rcshr_tree->Branch("_rc_vtx_z",&_rc_vtx_z,"rc_vtx_z/D");

  // reco shower info [rec] 
  _rcshr_tree->Branch("_rc_shr_x",&_rc_shr_x,"rc_shr_x/D");
  _rcshr_tree->Branch("_rc_shr_y",&_rc_shr_y,"rc_shr_y/D");
  _rcshr_tree->Branch("_rc_shr_z",&_rc_shr_z,"rc_shr_z/D");
  _rcshr_tree->Branch("_rc_shr_e",&_rc_shr_e,"rc_shr_e/D");
  _rcshr_tree->Branch("_rc_shr_dedx",&_rc_shr_dedx,"rc_shr_dedx/D");
  _rcshr_tree->Branch("_rc_shr_px",&_rc_shr_px,"rc_shr_px/D");
  _rcshr_tree->Branch("_rc_shr_py",&_rc_shr_py,"rc_shr_py/D");
  _rcshr_tree->Branch("_rc_shr_pz",&_rc_shr_pz,"rc_shr_pz/D");
  _rcshr_tree->Branch("_rcradlen",&_rcradlen,"rcradlen/D");
  // reco shower info [mc]
  _rcshr_tree->Branch("_mc_shr_x",&_mc_shr_x,"mc_shr_x/D");
  _rcshr_tree->Branch("_mc_shr_y",&_mc_shr_y,"mc_shr_y/D");
  _rcshr_tree->Branch("_mc_shr_z",&_mc_shr_z,"mc_shr_z/D");
  _rcshr_tree->Branch("_mc_shr_e",&_mc_shr_e,"mc_shr_e/D");
  _rcshr_tree->Branch("_mc_shr_edep",&_mc_shr_edep,"mc_shr_edep/D");
  _rcshr_tree->Branch("_mc_shr_dedx",&_mc_shr_dedx,"mc_shr_dedx/D");
  _rcshr_tree->Branch("_mc_shr_px",&_mc_shr_px,"mc_shr_px/D");
  _rcshr_tree->Branch("_mc_shr_py",&_mc_shr_py,"mc_shr_py/D");
  _rcshr_tree->Branch("_mc_shr_pz",&_mc_shr_pz,"mc_shr_pz/D");
  _rcshr_tree->Branch("_mcradlen",&_mcradlen,"mcradlen/D");

  // MC -> RC shower comparisons                                          
  _rcshr_tree->Branch("_angle",&_angle,"angle/D");
  _rcshr_tree->Branch("_strtdiff",&_strtdiff,"strtdiff/D");
  _rcshr_tree->Branch("_emc",&_emc,"emc/D");

  return;
}


void ShowerTruthStudies::beginJob()
{
  // Implementation of optional member function here.
}

void ShowerTruthStudies::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ShowerTruthStudies)
