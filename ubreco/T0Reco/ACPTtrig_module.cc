////////////////////////////////////////////////////////////////////////
// Class:       ACPTtrig
// Plugin Type: producer (art v3_00_00)
// File:        ACPTtrig_module.cc
//
// Generated at Wed Jan  2 13:44:06 2019 by David Caratelli using cetskelgen
// from cetlib version v3_04_00.
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

#include <memory>

// data-products
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include <TTree.h>

class ACPTtrig;


class ACPTtrig : public art::EDProducer {
public:
  explicit ACPTtrig(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ACPTtrig(ACPTtrig const&) = delete;
  ACPTtrig(ACPTtrig&&) = delete;
  ACPTtrig& operator=(ACPTtrig const&) = delete;
  ACPTtrig& operator=(ACPTtrig&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  double fCathodeMin, fCathodeMax;
  double fAnodeMin, fAnodeMax;
  double fYMin, fYMax;
  double fBeamSpillStart, fBeamSpillEnd;
  std::string fTrackProducer, fFlashProducer;
  bool fSaveTree;

  TTree* _tree;
  double _trk_beg_x, _trk_beg_y, _trk_beg_z;
  double _trk_end_x, _trk_end_y, _trk_end_z;
  double _flash_t, _flash_pe, _flash_yw, _flash_yc, _flash_zw, _flash_zc;
  int    _nflash;
  std::vector<double> _flash_pe_v;


};


ACPTtrig::ACPTtrig(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{

  produces< std::vector< anab::T0 > >();
  produces< art::Assns <recob::Track, anab::T0> >();
  produces< art::Assns <recob::Track, recob::OpFlash> >();

  fTrackProducer  = p.get<std::string>("TrackProducer");
  fFlashProducer  = p.get<std::string>("FlashProducer");
  fCathodeMin     = p.get<double>("CathodeMin");
  fCathodeMax     = p.get<double>("CathodeMax");
  fAnodeMin       = p.get<double>("AnodeMin");
  fAnodeMax       = p.get<double>("AnodeMax");
  fYMin           = p.get<double>("YMin");
  fYMax           = p.get<double>("YMax");
  fBeamSpillStart = p.get<double>("BeamSpillStart");
  fBeamSpillEnd   = p.get<double>("BeamSpillEnd");
  fSaveTree       = p.get<bool>("SaveTree",false);

  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("_tree","ACPT trigger tree");
  _tree->Branch("_trk_beg_x",&_trk_beg_x,"trk_beg_x/D");
  _tree->Branch("_trk_beg_y",&_trk_beg_y,"trk_beg_y/D");
  _tree->Branch("_trk_beg_z",&_trk_beg_z,"trk_beg_z/D");
  _tree->Branch("_trk_end_x",&_trk_end_x,"trk_end_x/D");
  _tree->Branch("_trk_end_y",&_trk_end_y,"trk_end_y/D");
  _tree->Branch("_trk_end_z",&_trk_end_z,"trk_end_z/D");
  _tree->Branch("_nflash",&_nflash,"nflash/I");
  _tree->Branch("_flash_t",&_flash_t,"flash_t/D");
  _tree->Branch("_flash_pe",&_flash_pe,"flash_pe/D");
  _tree->Branch("_flash_yw",&_flash_yw,"flash_yw/D");
  _tree->Branch("_flash_zw",&_flash_zw,"flash_zw/D");
  _tree->Branch("_flash_yc",&_flash_yc,"flash_yc/D");
  _tree->Branch("_flash_zc",&_flash_zc,"flash_zc/D");
  
}

void ACPTtrig::produce(art::Event& e)
{

  // produce OpFlash data-product to be filled within module
  std::unique_ptr< std::vector<anab::T0> > T0_v(new std::vector<anab::T0>);
  std::unique_ptr< art::Assns <recob::Track, anab::T0> >       trk_t0_assn_v   ( new art::Assns<recob::Track, anab::T0>       );
  std::unique_ptr< art::Assns <recob::Track, recob::OpFlash> > trk_flash_assn_v( new art::Assns<recob::Track, recob::OpFlash> );

  // load Flash
  art::Handle<std::vector<recob::OpFlash> > flash_h;
  e.getByLabel(fFlashProducer,flash_h);

  // make sure flash look good
  if(!flash_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate Flash!"<<std::endl;
    throw std::exception();
  }

  //  select largest flash within beam-spill window, if available
  size_t flash_ctr = 0;
  double PEmax = 0;
  size_t flash_idx = 0;
  for (size_t f=0; f < flash_h->size(); f++) {
    auto const& flash = flash_h->at(f);
    if ( (flash.Time() > fBeamSpillStart) && (flash.Time() < fBeamSpillEnd) ) {
      flash_ctr += 1;
      if (flash.TotalPE() > PEmax) {
	PEmax = flash.TotalPE();
	flash_idx = f;
      }// if largest flash so far
    }// if within beam spill window
  }// for all flashes

  // load tracks previously created for which T0 reconstruction should occur
  art::Handle<std::vector<recob::Track> > track_h;
  e.getByLabel(fTrackProducer,track_h);

  // make sure tracks look good
  if(!track_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate Track!"<<std::endl;
    throw std::exception();
  }

  std::vector<art::Ptr<recob::Track> > TrkVec;
  art::fill_ptr_vector(TrkVec, track_h);
  
  // loop through tracks
  for (auto& track : TrkVec){

    auto const& beg = track->Vertex();
    auto const& end = track->End();

    bool tagged = false;

    // cathode-side tracks
    // for tracks that enter from the side:
    if ( (beg.X() > end.X()) && (beg.Y() < fYMax) && (end.Y() < fYMin) ) {
      if ( (beg.X() > fCathodeMin) && (beg.X() < fCathodeMax) ) {
	tagged = true;
      }// if within cathode X constraint
    }// for tracks that enter from the side
    // for tracks that exit from the side:
    if ( (beg.X() < end.X()) && (beg.Y() > fYMax) && (end.Y() > fYMin) ) {
      if ( (end.X() > fCathodeMin) && (end.X() < fCathodeMax) ) {
	tagged = true;
      }// if within cathode X constraint
    }// for tracks that exit from the side

    // anode-side tracks
    // for tracks that enter from the side:
    if ( (beg.X() < end.X()) && (beg.Y() < fYMax) && (end.Y() < fYMin) ) {
      if ( (beg.X() > fAnodeMin) && (beg.X() < fAnodeMax) ) {
	tagged = true;
      }// if within anode X constraint
    }// for tracks that enter from the side
    // for tracks that exit from the side:
    if ( (beg.X() > end.X()) && (beg.Y() > fYMax) && (end.Y() > fYMin) ) {
      if ( (end.X() > fAnodeMin) && (end.X() < fAnodeMax) ) {
	tagged = true;
      }// if within anode X constraint
    }// for tracks that exit from the side
    
    // has the track been tagged?
    if (tagged == true) {

      double time = 0.;
      // associate  flash if a flash exits
      if (flash_ctr > 0) {
	time = flash_h->at(flash_idx).Time();
	// get pointer to associated flash
	const art::Ptr<recob::OpFlash> flash_ptr(flash_h, flash_idx);
	trk_flash_assn_v->addSingle(track, flash_ptr);
      }// if there is an associated flash

      // create T0 object with this information!
      anab::T0 t0(time, 0, 0);
      T0_v->emplace_back(t0);
      util::CreateAssn(*this, e, *T0_v, track, *trk_t0_assn_v);

      // save TTree
      if (fSaveTree) {
	_trk_beg_x = beg.X();
	_trk_end_x = end.X();
	_trk_beg_y = beg.y();
	_trk_end_y = end.y();
	_trk_beg_z = beg.Z();
	_trk_end_z = end.Z();
	_nflash    = flash_ctr;
	_flash_t   = time;
	if (flash_ctr > 0) {
	  auto const& flash = flash_h->at(flash_idx);
	  _flash_pe   = flash.TotalPE();
	  _flash_pe_v = flash.PEs();
	  _flash_yc   = flash.YCenter();
	  _flash_yw   = flash.YWidth();
	  _flash_zc   = flash.ZCenter();
	  _flash_zw   = flash.ZWidth();
	}// if a flash was found
	_tree->Fill();
      }// if we are saving to the TTree

    }// if the track was tagged

  }// for all tracks

  e.put(std::move(T0_v));
  e.put(std::move(trk_t0_assn_v));
  e.put(std::move(trk_flash_assn_v));

  return;
}

void ACPTtrig::beginJob()
{
  // Implementation of optional member function here.
}

void ACPTtrig::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ACPTtrig)
