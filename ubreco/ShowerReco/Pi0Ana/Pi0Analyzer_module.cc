////////////////////////////////////////////////////////////////////////
// Class:       Pi0Analyzer
// Plugin Type: analyzer (art v2_09_06)
// File:        Pi0Analyzer_module.cc
//
// Generated at Wed Mar 14 09:23:33 2018 by David Caratelli using cetskelgen
// from cetlib version v3_01_03.
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

#include "Selection/SelectionAlg.h"

#include <memory>

#include "TTree.h"
#include "TVector3.h"

#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

class Pi0Analyzer;


class Pi0Analyzer : public art::EDAnalyzer {
public:
  explicit Pi0Analyzer(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Pi0Analyzer(Pi0Analyzer const &) = delete;
  Pi0Analyzer(Pi0Analyzer &&) = delete;
  Pi0Analyzer & operator = (Pi0Analyzer const &) = delete;
  Pi0Analyzer & operator = (Pi0Analyzer &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  void Reset();

// Match each MCshower to its best matching RCshower
  std::vector<int> Match(const std::vector<sim::MCShower>& pi0_shower_v,
			 const std::vector<recob::Shower>& shr_v);
// Match each RCshower to its best matching MCshower
  std::vector<int> Match(const std::vector<recob::Shower>& shr_v,
			 const std::vector<sim::MCShower>& pi0_shower_v);
			 
  
  // unused double _w2cm, _t2cm;

  selection::SelectionAlg _pi0selection;

  void ClearMCRC();

  void SetTTree();

  // shower-by-shower TTree
  TTree* _mcshr_tree;
  TTree* _rcshr_tree;
  // pi0-by-pi0 ttree
  TTree* _pi0_tree;

  // variables common to both ttrees
  int _run, _sub, _evt;

  int _n_reco_showers;

  double _nu_e, _pi0_e;

  double _mc_vtx_x, _mc_vtx_y, _mc_vtx_z;
  double _rc_vtx_x, _rc_vtx_y, _rc_vtx_z;

  double _mc_shr1_x,  _mc_shr1_y,  _mc_shr1_z;
  double _mc_shr1_px, _mc_shr1_py, _mc_shr1_pz;
  double _mc_shr1_e, _mc_shr1_edep;
  double _mc_shr1_cont;
  double _mc_shr2_x,  _mc_shr2_y,  _mc_shr2_z;
  double _mc_shr2_px, _mc_shr2_py, _mc_shr2_pz;
  double _mc_shr2_e, _mc_shr2_edep;
  double _mc_shr2_cont;
  double _mcradlen1, _mcradlen2;

  double _mcangle;
  double _emc;
  double _mcmass, _mcmass_edep;

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

  // pi0 tree information
  double _rc_shr1_x,  _rc_shr1_y,  _rc_shr1_z;
  double _rc_shr1_px, _rc_shr1_py, _rc_shr1_pz;
  double _rc_shr1_e;
  double _rc_shr1_dedx;
  double _rcradlen1;
  double _angle1;
  double _strtdiff1;
  // unused double _dwallmin1;

  double _rc_shr2_x,  _rc_shr2_y,  _rc_shr2_z;
  double _rc_shr2_px, _rc_shr2_py, _rc_shr2_pz;
  double _rc_shr2_e;
  double _rc_shr2_dedx;
  double _rcradlen2;
  double _angle2;
  double _strtdiff2;
  // unused double _dwallmin2;

  double _rcmass;
  double _rcangle;

  std::string fShrProducer;

};


Pi0Analyzer::Pi0Analyzer(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  fShrProducer = p.get<std::string>("ShrProducer");
  SetTTree();
}

void Pi0Analyzer::analyze(art::Event const & e)
{

  _run = e.run();
  _sub = e.subRun();
  _evt = e.event();

  // load input showers
  auto const& shr_h = e.getValidHandle<std::vector<recob::Shower>>(fShrProducer);
  // load input tracks
  //auto const& trk_h = e.getValidHandle<std::vector<recob::Track>>("pandoraCosmic");
  // load input vertices
  auto const& vtx_h = e.getValidHandle<std::vector<recob::Vertex>>("ccvertex");
  // load mcshowers
  auto const& mcs_h = e.getValidHandle<std::vector<sim::MCShower>>("mcreco");

  // store reco vertex info
  if (vtx_h->size() == 1){
    Double_t rcxyz[3] = {};
    auto const& vtx = vtx_h->at(0);
    vtx.XYZ(rcxyz);
    _rc_vtx_x = rcxyz[0];
    _rc_vtx_y = rcxyz[1];
    _rc_vtx_z = rcxyz[2];
  }

  auto const& mct_h = e.getValidHandle<std::vector<simb::MCTruth> >("generator");

  // Store MC and RC showers in vectors
  _n_reco_showers = shr_h->size();
  std::vector<recob::Shower> reco_shower_v;
  for (size_t i=0; i < shr_h->size(); i++)
    reco_shower_v.push_back( shr_h->at(i) );
  std::vector<sim::MCShower> pi0_shower_v;

  auto mct = mct_h->at(0);
  size_t npart = mct.NParticles();

  bool foundPi0 = false;

  for (size_t i=0; i < npart; i++){
    auto const& part = mct.GetParticle(i);
    if ( (part.PdgCode() == 111) and (part.StatusCode() == 1) ){
      _mc_vtx_x = part.Trajectory().X(0);
      _mc_vtx_y = part.Trajectory().Y(0);
      _mc_vtx_z = part.Trajectory().Z(0);
      foundPi0 = true;
      break;
    }
  }

  if (foundPi0 == true) {
    
    size_t idx_1 = 0;
    size_t idx_2 = 0;
    size_t n_found = 0;
    for (size_t i=0; i < mcs_h->size(); i++){
      auto const& mcs = mcs_h->at(i);
      // distance from vertex                                                                
      double x = mcs.Start().X();
      double y = mcs.Start().Y();
      double z = mcs.Start().Z();
      double d = sqrt( ( (_mc_vtx_x - x) * (_mc_vtx_x - x) ) +
		       ( (_mc_vtx_y - y) * (_mc_vtx_y - y) ) +
		       ( (_mc_vtx_z - z) * (_mc_vtx_z - z) ) );
      if ( d < 0.01 ){
	if (n_found == 0){
	  idx_1 = i;
	  n_found += 1;
	}
	else if (n_found == 1){
	  idx_2 = i;
	  n_found += 1;
	}
	else
	  n_found += 1;
      }// if mother is a Pi0   
    }// for all MCShowers                                                               
    
    
    size_t idxLARGE = idx_1;
    size_t idxSMALL = idx_2;
    
    if (mcs_h->at(idx_1).Start().E() < mcs_h->at(idx_2).Start().E() )
      { idxLARGE = idx_2; idxSMALL = idx_1; }
    
    auto const& mcshr1 = mcs_h->at(idxLARGE);
    auto const& mcshr2 = mcs_h->at(idxSMALL);
    
    _pi0_e  = mcshr1.MotherEnd().E();
    
    _mc_shr1_e  = mcshr1.Start().E();
    _mc_shr1_edep  = mcshr1.DetProfile().E();
    _mc_shr1_x  = mcshr1.DetProfile().X();
    _mc_shr1_y  = mcshr1.DetProfile().Y();
    _mc_shr1_z  = mcshr1.DetProfile().Z();

    double mom1 = mcshr1.Start().Momentum().Vect().Mag();    
    _mc_shr1_px = mcshr1.Start().Px() / mom1;
    _mc_shr1_py = mcshr1.Start().Py() / mom1;
    _mc_shr1_pz = mcshr1.Start().Pz() / mom1;
    
    _mc_shr2_e  = mcshr2.Start().E();
    _mc_shr2_edep  = mcshr2.DetProfile().E();
    _mc_shr2_x  = mcshr2.DetProfile().X();
    _mc_shr2_y  = mcshr2.DetProfile().Y();
    _mc_shr2_z  = mcshr2.DetProfile().Z();

    double mom2 = mcshr2.Start().Momentum().Vect().Mag();    
    _mc_shr2_px = mcshr2.Start().Px() / mom2;
    _mc_shr2_py = mcshr2.Start().Py() / mom2;
    _mc_shr2_pz = mcshr2.Start().Pz() / mom2;
    
    _mcangle  = mcshr1.Start().Momentum().Vect().Angle( mcshr2.Start().Momentum().Vect() );
    
    _mcradlen1 = sqrt( ( (_mc_shr1_x - _mc_vtx_x) * (_mc_shr1_x - _mc_vtx_x) ) +
		       ( (_mc_shr1_y - _mc_vtx_y) * (_mc_shr1_y - _mc_vtx_y) ) +
		       ( (_mc_shr1_z - _mc_vtx_z) * (_mc_shr1_z - _mc_vtx_z) ) );
    
    _mcradlen2 = sqrt( ( (_mc_shr2_x - _mc_vtx_x) * (_mc_shr2_x - _mc_vtx_x) ) +
		       ( (_mc_shr2_y - _mc_vtx_y) * (_mc_shr2_y - _mc_vtx_y) ) +
		       ( (_mc_shr2_z - _mc_vtx_z) * (_mc_shr2_z - _mc_vtx_z) ) );
    
    _mcmass      = sqrt( 2 * _mc_shr1_e    * _mc_shr2_e    * (1 - cos(_mcangle) ) );
    _mcmass_edep = sqrt( 2 * _mc_shr1_edep * _mc_shr2_edep * (1 - cos(_mcangle) ) );

    pi0_shower_v = {mcshr1,mcshr2};

    // MC <-> RC matching
    auto MCRCmatch_v = Match(pi0_shower_v, reco_shower_v);
    
    for (size_t mcidx = 0; mcidx < MCRCmatch_v.size(); mcidx++) {

      ClearMCRC();
      
      auto mcshr = pi0_shower_v.at(mcidx);

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
  }// if a pi0 was found
  
  else {
    
    _mcmass  = -1;
    _mcmass_edep = -1;
    _mcangle = -5;
    
  }// if a pi0 was not found in the event


  // RC <-> MC matching
  auto RCMCmatch_v = Match(reco_shower_v, pi0_shower_v);
  
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
    
    auto const& mcshr = pi0_shower_v.at( RCMCmatch_v.at(rcidx) );

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
  
  // apply pi0 selection to event showers
  auto pi0candidate = _pi0selection.ApplySelection(shr_h);
  
  if (pi0candidate.mass < 0) {
    _rcmass  = pi0candidate.mass;
    _rcangle = -5;
  }
  else {

    _rcmass = pi0candidate.mass;
    _rcangle = pi0candidate.angle;
    _rc_shr1_e = pi0candidate.e1;
    _rc_shr2_e = pi0candidate.e2;
    _rc_shr1_dedx = pi0candidate.dedx1;
    _rc_shr2_dedx = pi0candidate.dedx2;

    auto const& rcshr1 = shr_h->at(pi0candidate.idx1);
    auto const& rcshr2 = shr_h->at(pi0candidate.idx2);

    _rc_shr1_x = rcshr1.ShowerStart().X();
    _rc_shr1_y = rcshr1.ShowerStart().Y();
    _rc_shr1_z = rcshr1.ShowerStart().Z();
    double momrc1 = rcshr1.Direction().Mag();      
    _rc_shr1_px = rcshr1.Direction().X() / momrc1;
    _rc_shr1_py = rcshr1.Direction().Y() / momrc1;
    _rc_shr1_pz = rcshr1.Direction().Z() / momrc1;

    _rc_shr2_x = rcshr2.ShowerStart().X();
    _rc_shr2_y = rcshr2.ShowerStart().Y();
    _rc_shr2_z = rcshr2.ShowerStart().Z();
    double momrc2 = rcshr2.Direction().Mag();      
    _rc_shr2_px = rcshr2.Direction().X() / momrc2;
    _rc_shr2_py = rcshr2.Direction().Y() / momrc2;
    _rc_shr2_pz = rcshr2.Direction().Z() / momrc2;

    _rcradlen1 = sqrt( ( (_rc_shr1_x - _rc_vtx_x) * (_rc_shr1_x - _rc_vtx_x) ) +
		       ( (_rc_shr1_y - _rc_vtx_y) * (_rc_shr1_y - _rc_vtx_y) ) +
		       ( (_rc_shr1_z - _rc_vtx_z) * (_rc_shr1_z - _rc_vtx_z) ) );

    _rcradlen2 = sqrt( ( (_rc_shr2_x - _rc_vtx_x) * (_rc_shr2_x - _rc_vtx_x) ) +
		       ( (_rc_shr2_y - _rc_vtx_y) * (_rc_shr2_y - _rc_vtx_y) ) +
		       ( (_rc_shr2_z - _rc_vtx_z) * (_rc_shr2_z - _rc_vtx_z) ) );

  }

  _pi0_tree->Fill();
  
  return;
}

void Pi0Analyzer::beginJob()
{
  // Implementation of optional member function here.
}

void Pi0Analyzer::endJob()
{
  // Implementation of optional member function here.
}

void Pi0Analyzer::Reset() {

  _n_reco_showers = 0;
  _nu_e = 0;
  _pi0_e = 0;

  _mc_vtx_x= _mc_vtx_y= _mc_vtx_z= 0;
  _rc_vtx_x= _rc_vtx_y= _rc_vtx_z= 0;

  _mc_shr1_x=  _mc_shr1_y=  _mc_shr1_z= 0;
  _mc_shr1_px= _mc_shr1_py= _mc_shr1_pz= 0;
  _mc_shr1_e= _mc_shr1_edep = 0;

  _mc_shr2_x=  _mc_shr2_y=  _mc_shr2_z= 0;
  _mc_shr2_px= _mc_shr2_py= _mc_shr2_pz= 0;
  _mc_shr2_e= _mc_shr2_edep = 0;

  _rc_shr1_x=  _rc_shr1_y=  _rc_shr1_z= 0;
  _rc_shr1_px= _rc_shr1_py= _rc_shr1_pz= 0;
  _rc_shr1_e= 0;

  _rc_shr2_x=  _rc_shr2_y=  _rc_shr2_z= 0;
  _rc_shr2_px= _rc_shr2_py= _rc_shr2_pz= 0;
  _rc_shr2_e= 0;

  _rc_shr_x=  _rc_shr_y=  _rc_shr_z= 0;
  _rc_shr_px= _rc_shr_py= _rc_shr_pz= 0;
  _rc_shr_e= 0;

  _mc_shr_x=  _mc_shr_y=  _mc_shr_z= 0;
  _mc_shr_px= _mc_shr_py= _mc_shr_pz= 0;
  _mc_shr_e=  _mc_shr_edep = 0;

  return;
}

void Pi0Analyzer::ClearMCRC() {

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



// Match each MCshower to its best matching RCshower
std::vector<int> Pi0Analyzer::Match(const std::vector<sim::MCShower>& pi0_shower_v,
				    const std::vector<recob::Shower>&   shr_v) {


  // now match to true showers                                                                                                                                                                                           
  std::vector<int> matched_indices;

  // find best matching RC shower for each MC shower               
  for (auto const& mcshr : pi0_shower_v){

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
std::vector<int> Pi0Analyzer::Match(const std::vector<recob::Shower>& shr_v,
				    const std::vector<sim::MCShower>& pi0_shower_v) {


  // now match to true showers                                                                                                                                                                                           
  std::vector<int> matched_indices;

  // find best matching MC shower for each RC shower       
  for (auto const& rcshr : shr_v) {
    
    double dotmax = -1.;
    int idxmax = -1;
    
    // loop through reco showers
    for (size_t i=0; i < pi0_shower_v.size(); i++) {

      auto const& mcshr = pi0_shower_v.at(i);
      
      double dot = rcshr.Direction().Dot( mcshr.Start().Momentum().Vect() );
      dot /= mcshr.Start().Momentum().Vect().Mag();
      dot /= rcshr.Direction().Mag();
      
      if (dot > dotmax) { dotmax = dot; idxmax = i; }
      
    }// for all reco indices                                                 
    
    matched_indices.push_back( idxmax );
  }
  
  return matched_indices;
  
}// end of function   


void Pi0Analyzer::SetTTree() {

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
  // mc shower info                                          
  _mcshr_tree->Branch("_mc_shr1_x",&_mc_shr1_x,"mc_shr1_x/D");
  _mcshr_tree->Branch("_mc_shr1_y",&_mc_shr1_y,"mc_shr1_y/D");
  _mcshr_tree->Branch("_mc_shr1_z",&_mc_shr1_z,"mc_shr1_z/D");
  _mcshr_tree->Branch("_mc_shr1_e",&_mc_shr1_e,"mc_shr1_e/D");
  _mcshr_tree->Branch("_mc_shr1_edep",&_mc_shr1_edep,"mc_shr1_edep/D");
  _mcshr_tree->Branch("_mc_shr1_cont",&_mc_shr1_cont,"mc_shr1_cont/D");
  _mcshr_tree->Branch("_mc_shr1_px",&_mc_shr1_px,"mc_shr1_px/D");
  _mcshr_tree->Branch("_mc_shr1_py",&_mc_shr1_py,"mc_shr1_py/D");
  _mcshr_tree->Branch("_mc_shr1_pz",&_mc_shr1_pz,"mc_shr1_pz/D");
  _mcshr_tree->Branch("_mcradlen1",&_mcradlen1,"mcradlen1/D");
  _mcshr_tree->Branch("_mc_shr2_x",&_mc_shr2_x,"mc_shr2_x/D");
  _mcshr_tree->Branch("_mc_shr2_y",&_mc_shr2_y,"mc_shr2_y/D");
  _mcshr_tree->Branch("_mc_shr2_z",&_mc_shr2_z,"mc_shr2_z/D");
  _mcshr_tree->Branch("_mc_shr2_e",&_mc_shr2_e,"mc_shr2_e/D");
  _mcshr_tree->Branch("_mc_shr2_edep",&_mc_shr2_edep,"mc_shr2_edep/D");
  _mcshr_tree->Branch("_mc_shr2_cont",&_mc_shr2_cont,"mc_shr2_cont/D");
  _mcshr_tree->Branch("_mc_shr2_px",&_mc_shr2_px,"mc_shr2_px/D");
  _mcshr_tree->Branch("_mc_shr2_py",&_mc_shr2_py,"mc_shr2_py/D");
  _mcshr_tree->Branch("_mc_shr2_pz",&_mc_shr2_pz,"mc_shr2_pz/D");
  _mcshr_tree->Branch("_mcradlen2",&_mcradlen2,"mcradlen2/D");

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

  // pi0 related MC information                                  
  _mcshr_tree->Branch("_nu_e",&_nu_e,"nu_e/D");
  _mcshr_tree->Branch("_pi0_e",&_pi0_e,"pi0_e/D");
  _mcshr_tree->Branch("_mcangle",&_mcangle,"mcangle/D");
  _mcshr_tree->Branch("_mcmass"  ,&_mcmass  ,"mcmass/D" );
  _mcshr_tree->Branch("_mcmass_edep"  ,&_mcmass_edep  ,"mcmass_edep/D" );
  _mcshr_tree->Branch("_dwallmin",&_dwallmin,"dwallmin/D");

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
  // mc shower info                                          
  _rcshr_tree->Branch("_mc_shr1_x",&_mc_shr1_x,"mc_shr1_x/D");
  _rcshr_tree->Branch("_mc_shr1_y",&_mc_shr1_y,"mc_shr1_y/D");
  _rcshr_tree->Branch("_mc_shr1_z",&_mc_shr1_z,"mc_shr1_z/D");
  _rcshr_tree->Branch("_mc_shr1_e",&_mc_shr1_e,"mc_shr1_e/D");
  _rcshr_tree->Branch("_mc_shr1_cont",&_mc_shr1_cont,"mc_shr1_cont/D");
  _rcshr_tree->Branch("_mc_shr1_px",&_mc_shr1_px,"mc_shr1_px/D");
  _rcshr_tree->Branch("_mc_shr1_py",&_mc_shr1_py,"mc_shr1_py/D");
  _rcshr_tree->Branch("_mc_shr1_pz",&_mc_shr1_pz,"mc_shr1_pz/D");
  _rcshr_tree->Branch("_mcradlen1",&_mcradlen1,"mcradlen1/D");
  _rcshr_tree->Branch("_mc_shr2_x",&_mc_shr2_x,"mc_shr2_x/D");
  _rcshr_tree->Branch("_mc_shr2_y",&_mc_shr2_y,"mc_shr2_y/D");
  _rcshr_tree->Branch("_mc_shr2_z",&_mc_shr2_z,"mc_shr2_z/D");
  _rcshr_tree->Branch("_mc_shr2_e",&_mc_shr2_e,"mc_shr2_e/D");
  _rcshr_tree->Branch("_mc_shr2_cont",&_mc_shr2_cont,"mc_shr2_cont/D");
  _rcshr_tree->Branch("_mc_shr2_px",&_mc_shr2_px,"mc_shr2_px/D");
  _rcshr_tree->Branch("_mc_shr2_py",&_mc_shr2_py,"mc_shr2_py/D");
  _rcshr_tree->Branch("_mc_shr2_pz",&_mc_shr2_pz,"mc_shr2_pz/D");
  _rcshr_tree->Branch("_mcradlen2",&_mcradlen2,"mcradlen2/D");

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

  // pi0 related MC information                                  
  _rcshr_tree->Branch("_nu_e",&_nu_e,"nu_e/D");
  _rcshr_tree->Branch("_pi0_e",&_pi0_e,"pi0_e/D");
  _rcshr_tree->Branch("_mcangle",&_mcangle,"mcangle/D");
  _rcshr_tree->Branch("_mcmass"  ,&_mcmass  ,"mcmass/D" );
  _rcshr_tree->Branch("_mcmass_edep"  ,&_mcmass_edep  ,"mcmass_edep/D" );
  _rcshr_tree->Branch("_dwallmin",&_dwallmin,"dwallmin/D");

  
  // pi0 ttree
  _pi0_tree = tfs->make<TTree>("_pi0_tree","Pi0 Tree TTree");

  _pi0_tree->Branch("_run",&_run,"run/I");
  _pi0_tree->Branch("_sub",&_sub,"sub/I");
  _pi0_tree->Branch("_evt",&_evt,"evt/I");

  _pi0_tree->Branch("_n_reco_showers",&_n_reco_showers,"n_reco_showers/I");
  // vertex info                                                                                           
  _pi0_tree->Branch("_mc_vtx_x",&_mc_vtx_x,"mc_vtx_x/D");
  _pi0_tree->Branch("_mc_vtx_y",&_mc_vtx_y,"mc_vtx_y/D");
  _pi0_tree->Branch("_mc_vtx_z",&_mc_vtx_z,"mc_vtx_z/D");
  _pi0_tree->Branch("_rc_vtx_x",&_rc_vtx_x,"rc_vtx_x/D");
  _pi0_tree->Branch("_rc_vtx_y",&_rc_vtx_y,"rc_vtx_y/D");
  _pi0_tree->Branch("_rc_vtx_z",&_rc_vtx_z,"rc_vtx_z/D");
  // mc shower info                                          
  _pi0_tree->Branch("_mc_shr1_x",&_mc_shr1_x,"mc_shr1_x/D");
  _pi0_tree->Branch("_mc_shr1_y",&_mc_shr1_y,"mc_shr1_y/D");
  _pi0_tree->Branch("_mc_shr1_z",&_mc_shr1_z,"mc_shr1_z/D");
  _pi0_tree->Branch("_mc_shr1_e",&_mc_shr1_e,"mc_shr1_e/D");
  _pi0_tree->Branch("_mc_shr1_edep",&_mc_shr1_edep,"mc_shr1_edep/D");
  _pi0_tree->Branch("_mc_shr1_cont",&_mc_shr1_cont,"mc_shr1_cont/D");
  _pi0_tree->Branch("_mc_shr1_px",&_mc_shr1_px,"mc_shr1_px/D");
  _pi0_tree->Branch("_mc_shr1_py",&_mc_shr1_py,"mc_shr1_py/D");
  _pi0_tree->Branch("_mc_shr1_pz",&_mc_shr1_pz,"mc_shr1_pz/D");
  _pi0_tree->Branch("_mc_shr2_x",&_mc_shr2_x,"mc_shr2_x/D");
  _pi0_tree->Branch("_mc_shr2_y",&_mc_shr2_y,"mc_shr2_y/D");
  _pi0_tree->Branch("_mc_shr2_z",&_mc_shr2_z,"mc_shr2_z/D");
  _pi0_tree->Branch("_mc_shr2_e",&_mc_shr2_e,"mc_shr2_e/D");
  _pi0_tree->Branch("_mc_shr2_edep",&_mc_shr2_edep,"mc_shr2_edep/D");
  _pi0_tree->Branch("_mc_shr2_cont",&_mc_shr2_cont,"mc_shr2_cont/D");
  _pi0_tree->Branch("_mc_shr2_px",&_mc_shr2_px,"mc_shr2_px/D");
  _pi0_tree->Branch("_mc_shr2_py",&_mc_shr2_py,"mc_shr2_py/D");
  _pi0_tree->Branch("_mc_shr2_pz",&_mc_shr2_pz,"mc_shr2_pz/D");

  // reco shower info      
  _pi0_tree->Branch("_rc_shr1_x",&_rc_shr1_x,"rc_shr1_x/D");
  _pi0_tree->Branch("_rc_shr1_y",&_rc_shr1_y,"rc_shr1_y/D");
  _pi0_tree->Branch("_rc_shr1_z",&_rc_shr1_z,"rc_shr1_z/D");
  _pi0_tree->Branch("_rc_shr1_e",&_rc_shr1_e,"rc_shr1_e/D");
  _pi0_tree->Branch("_rc_shr1_dedx",&_rc_shr1_dedx,"rc_shr1_dedx/D");
  _pi0_tree->Branch("_rc_shr1_px",&_rc_shr1_px,"rc_shr1_px/D");
  _pi0_tree->Branch("_rc_shr1_py",&_rc_shr1_py,"rc_shr1_py/D");
  _pi0_tree->Branch("_rc_shr1_pz",&_rc_shr1_pz,"rc_shr1_pz/D");
  _pi0_tree->Branch("_rc_shr2_x",&_rc_shr2_x,"rc_shr2_x/D");
  _pi0_tree->Branch("_rc_shr2_y",&_rc_shr2_y,"rc_shr2_y/D");
  _pi0_tree->Branch("_rc_shr2_z",&_rc_shr2_z,"rc_shr2_z/D");
  _pi0_tree->Branch("_rc_shr2_e",&_rc_shr2_e,"rc_shr2_e/D");
  _pi0_tree->Branch("_rc_shr2_dedx",&_rc_shr2_dedx,"rc_shr2_dedx/D");
  _pi0_tree->Branch("_rc_shr2_px",&_rc_shr2_px,"rc_shr2_px/D");
  _pi0_tree->Branch("_rc_shr2_py",&_rc_shr2_py,"rc_shr2_py/D");
  _pi0_tree->Branch("_rc_shr2_pz",&_rc_shr2_pz,"rc_shr2_pz/D");

  _pi0_tree->Branch("_rcradlen1",&_rcradlen1,"rcradlen1/D");
  _pi0_tree->Branch("_rcradlen2",&_rcradlen2,"rcradlen2/D");
  _pi0_tree->Branch("_mcradlen1",&_mcradlen1,"mcradlen1/D");
  _pi0_tree->Branch("_mcradlen2",&_mcradlen2,"mcradlen2/D");

  _pi0_tree->Branch("_rcmass",&_rcmass,"rcmass/D");
  _pi0_tree->Branch("_rcangle",&_rcangle,"rcangle/D");

  // MC -> RC shower comparisons                                          
  _pi0_tree->Branch("_angle1",&_angle1,"angle1/D");
  _pi0_tree->Branch("_strtdiff1",&_strtdiff1,"strtdiff1/D");
  _pi0_tree->Branch("_angle2",&_angle2,"angle2/D");
  _pi0_tree->Branch("_strtdiff2",&_strtdiff2,"strtdiff2/D");

  _pi0_tree->Branch("_emc",&_emc,"emc/D");

  // pi0 related MC information                                  
  _pi0_tree->Branch("_nu_e",&_nu_e,"nu_e/D");
  _pi0_tree->Branch("_pi0_e",&_pi0_e,"pi0_e/D");
  _pi0_tree->Branch("_mcangle",&_mcangle,"mcangle/D");
  _pi0_tree->Branch("_mcmass"  ,&_mcmass  ,"mcmass/D" );
  _pi0_tree->Branch("_mcmass_edep"  ,&_mcmass_edep  ,"mcmass_edep/D" );
  _pi0_tree->Branch("_dwallmin",&_dwallmin,"dwallmin/D");

  return;
}

DEFINE_ART_MODULE(Pi0Analyzer)
