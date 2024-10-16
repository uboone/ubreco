////////////////////////////////////////////////////////////////////////
// Class:       ThroughMu
// Plugin Type: analyzer (art v2_11_03)
// File:        ThroughMu_module.cc
//
// Generated at Sun Oct 21 21:38:26 2018 by David Caratelli using cetskelgen
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

#include "art/Persistency/Common/PtrMaker.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "art_root_io/TFileService.h"
#include "TTree.h"
#include "TVector3.h"

class ThroughMu;


class ThroughMu : public art::EDAnalyzer {
public:
  explicit ThroughMu(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ThroughMu(ThroughMu const &) = delete;
  ThroughMu(ThroughMu &&) = delete;
  ThroughMu & operator = (ThroughMu const &) = delete;
  ThroughMu & operator = (ThroughMu &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  bool ACPT();

  // Declare member data here.

  std::string fTrkProducer, fCaloProducer;

  TTree* _trk_tree;
  int   _pl;
  int   _run, _subrun, _event;
  float _trk_len;
  float _trk_start_x, _trk_start_y, _trk_start_z;
  float _trk_end_x, _trk_end_y, _trk_end_z;
  std::vector<float> _dqdx_v, _rr_v;

};


ThroughMu::ThroughMu(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{

  fTrkProducer  = p.get<std::string>("TrkProducer" );
  fCaloProducer = p.get<std::string>("CaloProducer");

  art::ServiceHandle<art::TFileService> tfs;

  _trk_tree = tfs->make<TTree>("_trk_tree","Track Tree");
  _trk_tree->Branch("_run",&_run,"run/I");
  _trk_tree->Branch("_event",&_event,"event/I");
  _trk_tree->Branch("_subrun",&_subrun,"subrun/I");
  _trk_tree->Branch("_trk_len",&_trk_len,"trk_len/F");
  _trk_tree->Branch("_pl",&_pl,"pl/I");
  _trk_tree->Branch("_trk_start_x",&_trk_start_x,"trk_start_x/F");
  _trk_tree->Branch("_trk_start_y",&_trk_start_y,"trk_start_y/F");
  _trk_tree->Branch("_trk_start_z",&_trk_start_z,"trk_start_z/F");
  _trk_tree->Branch("_trk_end_x",  &_trk_end_x,  "trk_end_x/F"  );
  _trk_tree->Branch("_trk_end_y",  &_trk_end_y,  "trk_end_y/F"  );
  _trk_tree->Branch("_trk_end_z",  &_trk_end_z,  "trk_end_z/F"  );
  _trk_tree->Branch("_dqdx_v","std::vector<float>",&_dqdx_v);
  _trk_tree->Branch("_rr_v",  "std::vector<float>",&_rr_v  );

}

void ThroughMu::analyze(art::Event const & e)
{

  _run = e.run();
  _subrun = e.subRun();
  _event = e.event();

  // load input tracks
  auto const& trk_h = e.getValidHandle<std::vector<recob::Track>>(fTrkProducer);
  // grab calorimetry objects associated to tracks
  art::FindMany<anab::Calorimetry> trk_calo_assn_v(trk_h, e, fCaloProducer);

  for (size_t t=0; t < trk_h->size(); t++) {
    
    auto const& trk = trk_h->at(t);
    auto const& beg = trk.Vertex();
    auto const& end = trk.End();
    
    _trk_len = trk.Length();
    _trk_start_x = beg.X();
    _trk_start_y = beg.Y();
    _trk_start_z = beg.Z();
    _trk_end_x   = end.X();
    _trk_end_y   = end.Y();
    _trk_end_z   = end.Z();
    
    // fill calorimetry info for this track
    // grab the associated calorimetry object
    const std::vector<const anab::Calorimetry*>& Calo_v = trk_calo_assn_v.at(t);
    
    for (size_t pl=0; pl < Calo_v.size(); pl++){
      
      auto const& calo = Calo_v.at(pl);
      
      auto const& plane = calo->PlaneID().Plane;
      
      _pl = plane;

      if (_pl != 2) continue;

      // apply cuts on start and end point for ACPT in beam cuts
      if (ACPT() == false) continue;
      
      // grab point-by-point information
      auto const& dqdx_v = calo->dQdx();
      auto const& rr_v   = calo->ResidualRange();
      _dqdx_v.clear();
      _rr_v.clear();
      for (size_t n=0; n < dqdx_v.size(); n++) {
	_dqdx_v.push_back((float)dqdx_v[n]);
	_rr_v.push_back(  (float)rr_v[n]  );
      }
      
      _trk_tree->Fill();
      
    }// for all planes
  }// for all tracks

  return;
}

bool ThroughMu::ACPT() {

  if ( (_trk_start_y < 95) && (_trk_end_y < -95) && (_trk_start_x > 250) && (_trk_start_x < 270) )
    return true;
  if ( (_trk_start_y < 95) && (_trk_end_y < -95) && (_trk_start_x > -10) && (_trk_start_x < 10) )
    return true;
  if ( (_trk_end_y > -95) && (_trk_start_y > 95) && (_trk_end_x > 250) && (_trk_end_x < 270) )
    return true;
  if ( (_trk_end_y > -95) && (_trk_start_y > 95) && (_trk_end_x > -10) && (_trk_end_x < 10) )
    return true;

  return false;
}

void ThroughMu::beginJob()
{
  // Implementation of optional member function here.
}

void ThroughMu::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ThroughMu)
