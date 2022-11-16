////////////////////////////////////////////////////////////////////////
// Class:       FilteredHitsTruthProducer
// Plugin Type: producer (art v3_06_03)
// File:        FilteredHitsTruthProducer_module.cc
//
// Generated at Tue May 25 10:39:19 2021 by Giuseppe Cerati using cetskelgen
// from cetlib version v3_11_01.
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
#include "canvas/Persistency/Common/FindManyP.h"

#include <memory>

#include "canvas/Persistency/Common/Assns.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

class FilteredHitsTruthProducer;

using HitParticleAssociations = art::Assns<simb::MCParticle, recob::Hit, anab::BackTrackerHitMatchingData>;

class FilteredHitsTruthProducer : public art::EDProducer {
public:
  explicit FilteredHitsTruthProducer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FilteredHitsTruthProducer(FilteredHitsTruthProducer const&) = delete;
  FilteredHitsTruthProducer(FilteredHitsTruthProducer&&) = delete;
  FilteredHitsTruthProducer& operator=(FilteredHitsTruthProducer const&) = delete;
  FilteredHitsTruthProducer& operator=(FilteredHitsTruthProducer&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.
  art::InputTag fHitLabel;
  art::InputTag fFiltHitLabel;
  art::InputTag fHitTruthLabel;
  float fRMSCut;
};


FilteredHitsTruthProducer::FilteredHitsTruthProducer(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fHitLabel(p.get<art::InputTag>("HitLabel")),
  fFiltHitLabel(p.get<art::InputTag>("FiltHitLabel")),
  fHitTruthLabel(p.get<art::InputTag>("HitTruthLabel")),
  fRMSCut(p.get<float>("RMScut",0.001))
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  produces<HitParticleAssociations>();

  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void FilteredHitsTruthProducer::produce(art::Event& e)
{
  // Implementation of required member function here.

  auto outputHitPartAssns = std::make_unique<HitParticleAssociations>();

  art::Handle< std::vector< recob::Hit > > hitListHandle;
  e.getByLabel(fHitLabel, hitListHandle);
  std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> hittruth;
  hittruth = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> >(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(hitListHandle, e, fHitTruthLabel));

  art::Handle< std::vector< recob::Hit > > filtHitListHandle;
  e.getByLabel(fFiltHitLabel, filtHitListHandle);
  std::vector<std::pair<float,size_t> > map[3][3456] = {std::vector<std::pair<float, size_t> >()};
  for (size_t ifilthit=0; ifilthit<filtHitListHandle->size();ifilthit++) {
    art::Ptr<recob::Hit> filthit(filtHitListHandle,ifilthit);
    map[filthit->WireID().Plane][filthit->WireID().Wire].push_back({filthit->PeakTime(),ifilthit});
  }

  for (size_t ihit=0; ihit<hitListHandle->size();ihit++) {
    art::Ptr<recob::Hit> hit(hitListHandle,ihit);

    for (auto ifilt : map[hit->WireID().Plane][hit->WireID().Wire]) {
      if ( std::fabs(hit->PeakTime()-ifilt.first)>fRMSCut*hit->RMS() ) continue;
      //
      const art::Ptr<recob::Hit> ahp(filtHitListHandle,ifilt.second);
      std::vector<art::Ptr<simb::MCParticle> > particle_vec = hittruth->at(hit.key());
      std::vector<anab::BackTrackerHitMatchingData const *> match_vec = hittruth->data(hit.key());
      for (size_t i_p = 0; i_p < particle_vec.size(); ++i_p) {
	outputHitPartAssns->addSingle(particle_vec[i_p],ahp,*match_vec[i_p]);
      }
      //
      break;
    }

  }

  e.put(std::move(outputHitPartAssns));

}

DEFINE_ART_MODULE(FilteredHitsTruthProducer)
