////////////////////////////////////////////////////////////////////////
// Class:       Pi0AnalyzerDATA
// Plugin Type: analyzer (art v2_09_06)
// File:        Pi0AnalyzerDATA_module.cc
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
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

class Pi0AnalyzerDATA;


class Pi0AnalyzerDATA : public art::EDAnalyzer {
public:
  explicit Pi0AnalyzerDATA(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Pi0AnalyzerDATA(Pi0AnalyzerDATA const &) = delete;
  Pi0AnalyzerDATA(Pi0AnalyzerDATA &&) = delete;
  Pi0AnalyzerDATA & operator = (Pi0AnalyzerDATA const &) = delete;
  Pi0AnalyzerDATA & operator = (Pi0AnalyzerDATA &&) = delete;

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

  // pi0-by-pi0 ttree
  TTree* _pi0_tree;
  TTree* _rcshr_tree;
  TTree* _trk_tree;

  // variables common to both ttrees
  int _run, _sub, _evt;

  int _n_reco_showers;


  double _rc_vtx_x, _rc_vtx_y, _rc_vtx_z;


  // shower-by-shower variables [reco]
  double _rc_shr_x,  _rc_shr_y,  _rc_shr_z;
  double _rc_shr_px, _rc_shr_py, _rc_shr_pz;
  double _rc_shr_e;
  double _rc_shr_dedx;
  double _rcradlen;
  // shower-by-shower variables, comparisons
  // unused double _angle;
  // unused double _strtdiff;
  // unused double _dwallmin;

  // pi0 tree information
  double _rc_shr1_x,  _rc_shr1_y,  _rc_shr1_z;
  double _rc_shr1_px, _rc_shr1_py, _rc_shr1_pz;
  double _rc_shr1_e;
  double _rc_shr1_dedx;
  double _rcradlen1;
  // unused double _angle1;
  // unused double _strtdiff1;
  // unused double _dwallmin1;

  double _rc_shr2_x,  _rc_shr2_y,  _rc_shr2_z;
  double _rc_shr2_px, _rc_shr2_py, _rc_shr2_pz;
  double _rc_shr2_e;
  double _rc_shr2_dedx;
  double _rcradlen2;
  // unused double _angle2;
  // unused double _strtdiff2;
  // unused double _dwallmin2;

  double _rcmass;
  double _rcangle;

  std::vector<float> _dqdx_v, _dedx_v, _rr_v;
  double _dvtx, _length;
  double _xs, _ys, _zs, _xe, _ye, _ze;
  int _pl;

  std::string fShrProducer;
  std::string fTrkProducer;
  std::string fCaloProducer;

};


Pi0AnalyzerDATA::Pi0AnalyzerDATA(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  fShrProducer  = p.get<std::string>("ShrProducer" );
  fTrkProducer  = p.get<std::string>("TrkProducer" );
  fCaloProducer = p.get<std::string>("CaloProducer");
  SetTTree();
}

void Pi0AnalyzerDATA::analyze(art::Event const & e)
{

  _run = e.run();
  _sub = e.subRun();
  _evt = e.event();

  // load input showers
  auto const& shr_h = e.getValidHandle<std::vector<recob::Shower>>(fShrProducer);
  // load input tracks
  auto const& trk_h = e.getValidHandle<std::vector<recob::Track>>(fTrkProducer);
  // load calorimetry associated to tracks
  art::FindManyP<anab::Calorimetry> trk_calo_assn_v(trk_h, e, fCaloProducer);
  // load input vertices
  auto const& vtx_h = e.getValidHandle<std::vector<recob::Vertex>>("ccvertex");

  // store reco vertex info
  if (vtx_h->size() == 1){
    Double_t rcxyz[3] = {};
    auto const& vtx = vtx_h->at(0);
    vtx.XYZ(rcxyz);
    _rc_vtx_x = rcxyz[0];
    _rc_vtx_y = rcxyz[1];
    _rc_vtx_z = rcxyz[2];
  }

  // track portion of this analysis
  TVector3 nuvtx(_rc_vtx_x,_rc_vtx_y,_rc_vtx_z);
  for (size_t t=0; t < trk_h->size(); t++) {
    
    auto const& trk = trk_h->at(t);
    auto const& beg = trk.Vertex();
    auto const& end = trk.End();

    _xs = beg.X();
    _ys = beg.Y();
    _zs = beg.Z();
    _xe = end.X();
    _ye = end.Y();
    _ze = end.Z();

    double dvtx = 1000.;
    // which one is closest?
    if ( (nuvtx-beg).Mag() < (nuvtx-end).Mag() ) 
      dvtx = (nuvtx-beg).Mag();
    else 
      dvtx = (nuvtx-end).Mag();

    if (dvtx < 5.0) {

      _length = trk.Length();
      _dvtx   = dvtx;
      
      // fill calorimetry info for this track
      // grab the associated calorimetry object
      //const std::vector<const anab::Calorimetry*>& Calo_v = trk_calo_assn_v.at(t);
      auto Calo_v = trk_calo_assn_v.at(t);
      
      for (size_t pl=0; pl < Calo_v.size(); pl++){
	
	auto const& calo = Calo_v.at(pl);
	
	_pl = calo->PlaneID().Plane;
	
	// grab point-by-point information
	auto const& dqdx_v = calo->dQdx();
	auto const& dedx_v = calo->dEdx();
	auto const& rr_v   = calo->ResidualRange();
	_dqdx_v.clear();
	_dedx_v.clear();
	_rr_v.clear();
	for (auto const& dqdx : dqdx_v)
	  _dqdx_v.push_back((float)dqdx);
	for (auto const& dedx : dedx_v)
	  _dedx_v.push_back((float)dedx);
	for (auto const& rr : rr_v)
	  _rr_v.push_back((float)rr);

	_trk_tree->Fill();
	
      }// for all planes
    }// if within 5 cm of vertex
  }// for all tracks

  // Store MC and RC showers in vectors
  _n_reco_showers = shr_h->size();

  for (size_t rcidx = 0; rcidx < shr_h->size(); rcidx++) {

    auto const& rcshr = shr_h->at(rcidx);
    
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

    _rcshr_tree->Fill();

  }

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

void Pi0AnalyzerDATA::beginJob()
{
  // Implementation of optional member function here.
}

void Pi0AnalyzerDATA::endJob()
{
  // Implementation of optional member function here.
}

void Pi0AnalyzerDATA::Reset() {

  _n_reco_showers = 0;

  _rc_vtx_x= _rc_vtx_y= _rc_vtx_z= 0;

  _rc_shr1_x=  _rc_shr1_y=  _rc_shr1_z= 0;
  _rc_shr1_px= _rc_shr1_py= _rc_shr1_pz= 0;
  _rc_shr1_e= 0;

  _rc_shr2_x=  _rc_shr2_y=  _rc_shr2_z= 0;
  _rc_shr2_px= _rc_shr2_py= _rc_shr2_pz= 0;
  _rc_shr2_e= 0;

  _rc_shr_x=  _rc_shr_y=  _rc_shr_z= 0;
  _rc_shr_px= _rc_shr_py= _rc_shr_pz= 0;
  _rc_shr_e= 0;

  return;
}

void Pi0AnalyzerDATA::ClearMCRC() {

  _rc_shr_x=  _rc_shr_y=  _rc_shr_z= 0;
  _rc_shr_px= _rc_shr_py= _rc_shr_pz= 0;
  _rc_shr_e= 0;
  _rc_shr_dedx = 0;

  return;
}



// Match each MCshower to its best matching RCshower
std::vector<int> Pi0AnalyzerDATA::Match(const std::vector<sim::MCShower>& pi0_shower_v,
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
std::vector<int> Pi0AnalyzerDATA::Match(const std::vector<recob::Shower>& shr_v,
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


void Pi0AnalyzerDATA::SetTTree() {

  art::ServiceHandle<art::TFileService> tfs;

  // MC shower-by-shower TTree
  _rcshr_tree = tfs->make<TTree>("_rcshr_tree","Pi0 Tree TTree");

  _rcshr_tree->Branch("_run",&_run,"run/I");
  _rcshr_tree->Branch("_sub",&_sub,"sub/I");
  _rcshr_tree->Branch("_evt",&_evt,"evt/I");

  _rcshr_tree->Branch("_n_reco_showers",&_n_reco_showers,"n_reco_showers/I");
  // vertex info
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


  // track tree
  _trk_tree = tfs->make<TTree>("_trk_tree","track tree");
  _trk_tree->Branch("_run",&_run,"run/I");
  _trk_tree->Branch("_sub",&_sub,"sub/I");
  _trk_tree->Branch("_evt",&_evt,"evt/I");
  _trk_tree->Branch("_dqdx_v","std::vector<float>",&_dqdx_v);
  _trk_tree->Branch("_dedx_v","std::vector<float>",&_dedx_v);
  _trk_tree->Branch("_rr_v"  ,"std::vector<float>",&_rr_v  );
  _trk_tree->Branch("_pl",&_pl,"pl/I");
  _trk_tree->Branch("_dvtx",&_dvtx,"dvtx/D");
  _trk_tree->Branch("_length",&_length,"length/D");
  _trk_tree->Branch("_xs",&_xs,"xs/D");
  _trk_tree->Branch("_ys",&_ys,"ys/D");
  _trk_tree->Branch("_zs",&_zs,"zs/D");
  _trk_tree->Branch("_xe",&_xe,"xe/D");
  _trk_tree->Branch("_ye",&_ye,"ye/D");
  _trk_tree->Branch("_ze",&_ze,"ze/D");
  
  // pi0 ttree
  _pi0_tree = tfs->make<TTree>("_pi0_tree","Pi0 Tree TTree");

  _pi0_tree->Branch("_run",&_run,"run/I");
  _pi0_tree->Branch("_sub",&_sub,"sub/I");
  _pi0_tree->Branch("_evt",&_evt,"evt/I");

  _pi0_tree->Branch("_n_reco_showers",&_n_reco_showers,"n_reco_showers/I");
  // vertex info                                                                                           
  _pi0_tree->Branch("_rc_vtx_x",&_rc_vtx_x,"rc_vtx_x/D");
  _pi0_tree->Branch("_rc_vtx_y",&_rc_vtx_y,"rc_vtx_y/D");
  _pi0_tree->Branch("_rc_vtx_z",&_rc_vtx_z,"rc_vtx_z/D");
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

  _pi0_tree->Branch("_rcmass",&_rcmass,"rcmass/D");
  _pi0_tree->Branch("_rcangle",&_rcangle,"rcangle/D");


  return;
}

DEFINE_ART_MODULE(Pi0AnalyzerDATA)
