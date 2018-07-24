////////////////////////////////////////////////////////////////////////
// Class:       CosmicBackgrounds
// Plugin Type: filter (art v2_09_06)
// File:        CosmicBackgrounds_module.cc
//
// Generated at Sat Feb 24 08:47:57 2018 by David Caratelli using cetskelgen
// from cetlib version v3_01_03.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Services/Optional/TFileService.h"

#include <memory>

#include "TTree.h"
#include "TVector3.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"

class CosmicBackgrounds;

class CosmicBackgrounds : public art::EDFilter {
public:
  explicit CosmicBackgrounds(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CosmicBackgrounds(CosmicBackgrounds const &) = delete;
  CosmicBackgrounds(CosmicBackgrounds &&) = delete;
  CosmicBackgrounds & operator = (CosmicBackgrounds const &) = delete;
  CosmicBackgrounds & operator = (CosmicBackgrounds &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  TTree* _trkangle_tree;
  float  _trkangle;

};


CosmicBackgrounds::CosmicBackgrounds(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  
}

bool CosmicBackgrounds::filter(art::Event & e)
{

  // load input tracks
  auto const& trk_h = e.getValidHandle<std::vector<recob::Track>>("pandoraCosmic");
  // load input vertices
  auto const& vtx_h = e.getValidHandle<std::vector<recob::Vertex>>("ccvertex");

  // check whether vertex in the middle of a cosmic track
  if (vtx_h->size() != 1) return false;

  auto const& vtx = vtx_h->at(0);
  Double_t xyz[3] = {};
  vtx.XYZ(xyz);
  TVector3 vtxpt(xyz[0],xyz[1],xyz[2]);
  _trkangle = 0.;

  // loop trhough tracks. is the angle ~180?
  for (size_t t=0; t < trk_h->size(); t++) {

    auto const& trk = trk_h->at(t);

    if (trk.Length() < 50) continue;

    auto const& start = trk.Vertex();
    auto const& end   = trk.End();

    if ( ( (start-vtxpt).Mag() < 10) && ((end-vtxpt).Mag() < 10) ) continue;

    float a = (start-vtxpt).Angle(end-vtxpt);
    if (a > _trkangle)
      _trkangle = a;
  }// for all tracks

  _trkangle_tree->Fill();

  return true;
}

void CosmicBackgrounds::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  _trkangle_tree = tfs->make<TTree>("_trkangle_tree","Track Angle TTree");
  _trkangle_tree->Branch("_trkangle",&_trkangle,"trkangle/F");

}

void CosmicBackgrounds::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(CosmicBackgrounds)
