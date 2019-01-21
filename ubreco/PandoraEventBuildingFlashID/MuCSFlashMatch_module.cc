////////////////////////////////////////////////////////////////////////
// Class:       MuCSFlashMatch
// Plugin Type: analyzer (art v3_00_00)
// File:        MuCSFlashMatch_module.cc
//
// Generated at Sun Jan 13 18:33:43 2019 by David Caratelli using cetskelgen
// from cetlib version v3_04_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Utilities/make_tool.h"
#include "FlashMatchingToolBase_tool.h"

#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Hit.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/FindManyInChainP.h"

#include "lardataobj/RecoBase/PFParticleMetadata.h"

#include "TTree.h"

// art TOOLS
#include "art/Utilities/ToolMacros.h"


class MuCSFlashMatch;


class MuCSFlashMatch : public art::EDAnalyzer {
public:
  explicit MuCSFlashMatch(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MuCSFlashMatch(MuCSFlashMatch const&) = delete;
  MuCSFlashMatch(MuCSFlashMatch&&) = delete;
  MuCSFlashMatch& operator=(MuCSFlashMatch const&) = delete;
  MuCSFlashMatch& operator=(MuCSFlashMatch&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  TTree* _tree;
  float _score;
  int _mucs;

  TTree* _evt_tree;
  float _best_score;
  float _mucs_score;
  int _best_mucs;

  std::string fPFPproducer, fTrackproducer, fSpacePointproducer, fCTagproducer;
  
  bool fOnlyTagged;
  
  std::unique_ptr<flashmatch::FlashMatchingToolBase> _flashmatchTool;             ///< The slice id tool

};


MuCSFlashMatch::MuCSFlashMatch(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fPFPproducer         = p.get<std::string>("PFPproducer");
  fTrackproducer       = p.get<std::string>("Trackproducer");
  fSpacePointproducer  = p.get<std::string>("SpacePointproducer");
  fCTagproducer        = p.get<std::string>("CTagproducer");
  fOnlyTagged          = p.get<bool       >("OnlyTagged");

  const fhicl::ParameterSet& flashconfig = p.get<fhicl::ParameterSet>("SliceTool");  
  
  _flashmatchTool = art::make_tool<flashmatch::FlashMatchingToolBase>(flashconfig);

  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("flashmatch","flashmatching tree");
  _tree->Branch("_score",&_score,"score/F");
  _tree->Branch("_mucs",&_mucs,"mucs/I");

  _evt_tree = tfs->make<TTree>("evtflashmatch","evt flashmatching tree");
  _evt_tree->Branch("_best_score",&_best_score,"best_score/F");
  _evt_tree->Branch("_mucs_score",&_mucs_score,"mucs_score/F");
  _evt_tree->Branch("_best_mucs",&_best_mucs,"best_mucs/I");

}

void MuCSFlashMatch::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  // grab pfp in event
  auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle> >(fPFPproducer);

  // grab tracks in event
  auto const& trk_h = e.getValidHandle<std::vector<recob::Track> >(fTrackproducer);
  
  // grab cosmic score association (to find MuCS tagged tracks)
  art::FindManyP<anab::CosmicTag> trk_ctag_assn_v(trk_h, e, fCTagproducer);

  // grab tracks associated to pfparticles
  art::FindManyP<recob::Track> pfp_trk_assn_v(pfp_h, e, fPFPproducer);
  
  // grab spacepoints associated with PFParticles
  art::FindManyP<recob::SpacePoint> pfp_spacepoint_assn_v(pfp_h, e, fPFPproducer);
  
  auto const& spacepoint_h = e.getValidHandle<std::vector<recob::SpacePoint> >(fSpacePointproducer);
  
  art::FindManyP<recob::Hit> spacepoint_hit_assn_v(spacepoint_h, e, fSpacePointproducer);

  _best_score = 10000.;
  _mucs_score = -1; // this way we know if a MuCS track was tagged in the event

  // loop through PFParticles
  for (unsigned int p=0; p < pfp_h->size(); p++){
    
    auto const& pfp = pfp_h->at(p);
    const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, p);
    auto pfp_key = pfp_ptr.key();

    // only primary PFParticles
    if (pfp.IsPrimary() == false) continue;
    
    // grab the associated track
    const std::vector< art::Ptr<recob::Track> >& trk_ptr_v = pfp_trk_assn_v.at(pfp_key);

    // are there assocaited tracks? if no skip
    if (trk_ptr_v.size() != 1) continue;

    // grab the track key
    auto trk_key = trk_ptr_v.at(0).key();

    const art::Ptr<recob::Track> trk_ptr(trk_h, trk_key);

    // associations
    const std::vector< art::Ptr<recob::SpacePoint> >& spacepoint_ptr_v = pfp_spacepoint_assn_v.at(pfp_key);
    const std::vector< art::Ptr<anab::CosmicTag> >& ctag_ptr_v = trk_ctag_assn_v.at(trk_key);
    
    _mucs = ctag_ptr_v.size();

    if (fOnlyTagged && (_mucs != 1) ) continue;

    
    std::vector< art::Ptr<recob::Hit> > hit_ptr_v;
    
    for (size_t sp=0; sp < spacepoint_ptr_v.size(); sp++) {
      auto const& spkey = spacepoint_ptr_v.at(sp).key();
      const std::vector< art::Ptr<recob::Hit> > this_hit_ptr_v = spacepoint_hit_assn_v.at( spkey );
      for (size_t h=0; h < this_hit_ptr_v.size(); h++) {
	hit_ptr_v.push_back( this_hit_ptr_v.at( h ) );
      }// for all hits associated to this spacepoint
    }// fpr all spacepoints
    
    
    _score = _flashmatchTool->ClassifyTrack(e, spacepoint_ptr_v, hit_ptr_v);

    if ( (_score < _best_score) && (_score > 0) ) { _best_score = _score; _best_mucs = _mucs; }
    if (_mucs) { _mucs_score = _score; }
    
    _tree->Fill();
    
  }// for all tracks

  _evt_tree->Fill();
  
  return;
}

void MuCSFlashMatch::beginJob()
{
  // Implementation of optional member function here.
}

void MuCSFlashMatch::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(MuCSFlashMatch)
