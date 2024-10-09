////////////////////////////////////////////////////////////////////////
// Class:       FilteredHitsProducer
// Plugin Type: producer (art v3_06_03)
// File:        FilteredHitsProducer_module.cc
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
#include "lardataobj/AnalysisBase/MVAOutput.h"

class FilteredHitsProducer;

using HitParticleAssociations = art::Assns<simb::MCParticle, recob::Hit, anab::BackTrackerHitMatchingData>;

class FilteredHitsProducer : public art::EDProducer {
public:
  explicit FilteredHitsProducer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FilteredHitsProducer(FilteredHitsProducer const&) = delete;
  FilteredHitsProducer(FilteredHitsProducer&&) = delete;
  FilteredHitsProducer& operator=(FilteredHitsProducer const&) = delete;
  FilteredHitsProducer& operator=(FilteredHitsProducer&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.
  art::InputTag fHitLabel;
  art::InputTag fMatchHitLabel;
  art::InputTag fHitScoreLabel;
  art::InputTag fHitSemanticLabel;
  art::InputTag fHitTruthLabel;
  float fRMSCut;
  float fScoreCut;
  bool fIsMC;
};


FilteredHitsProducer::FilteredHitsProducer(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fHitLabel(p.get<art::InputTag>("HitLabel")),
  fMatchHitLabel(p.get<art::InputTag>("MatchHitLabel")),
  fHitScoreLabel(p.get<art::InputTag>("HitScoreLabel")),
  fHitSemanticLabel(p.get<art::InputTag>("HitSemanticLabel")),
  fHitTruthLabel(p.get<art::InputTag>("HitTruthLabel")),
  fRMSCut(p.get<float>("RMScut")),
  fScoreCut(p.get<float>("ScoreCut")),
  fIsMC(p.get<bool>("isMC"))
{
  produces<std::vector<recob::Hit>>();
  if (fIsMC) produces<HitParticleAssociations>();

  if (fScoreCut>=0) produces<std::vector<anab::FeatureVector<5> > >("semantic");
}

void FilteredHitsProducer::produce(art::Event& e)
{
  // Implementation of required member function here.

  auto outputHits = std::make_unique<std::vector<recob::Hit> >();
  auto outputHitPartAssns = std::make_unique<HitParticleAssociations>();
  std::unique_ptr<std::vector<anab::FeatureVector<5>>> semtcol(new std::vector<anab::FeatureVector<5> >());
  art::PtrMaker<recob::Hit> hitPtrMaker(e);

  art::Handle< std::vector< recob::Hit > > hitListHandle;
  e.getByLabel(fHitLabel, hitListHandle);
  std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> hittruth;
  if (fIsMC) hittruth = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> >(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(hitListHandle, e, fHitTruthLabel));

  art::Handle< std::vector< recob::Hit > > mtchHitListHandle;
  std::vector<float> map[3][3456] = {{std::vector<float>()}};
  if (fRMSCut>=0) {
    e.getByLabel(fMatchHitLabel, mtchHitListHandle);
    for (size_t imtchhit=0; imtchhit<mtchHitListHandle->size();imtchhit++) {
      art::Ptr<recob::Hit> mtchhit(mtchHitListHandle,imtchhit);
      map[mtchhit->WireID().Plane][mtchhit->WireID().Wire].push_back(mtchhit->PeakTime());
    }
  }

  art::Handle< std::vector< anab::FeatureVector<1> > > ngFiltOutputHandle;
  art::Handle< std::vector< anab::FeatureVector<5> > > ngSemtOutputHandle;
  if (fScoreCut>=0) {
    e.getByLabel(fHitScoreLabel, ngFiltOutputHandle);
    e.getByLabel(fHitSemanticLabel, ngSemtOutputHandle);
  }

  for (size_t ihit=0; ihit<hitListHandle->size();ihit++) {
    if (fScoreCut>=0 && ngFiltOutputHandle->size()==0) break;
    art::Ptr<recob::Hit> hit(hitListHandle,ihit);

    bool keepHit = false;
    for (auto mtcht : map[hit->WireID().Plane][hit->WireID().Wire]) {
      if ( std::fabs(hit->PeakTime()-mtcht)>fRMSCut*hit->RMS() ) continue;
      keepHit = true;
      break;
    }
    if (fScoreCut>=0 && ngFiltOutputHandle->at(ihit).at(0)>=fScoreCut) keepHit = true;
    if (keepHit == false) continue;
    //
    outputHits->emplace_back(*hit);
    if (fScoreCut>=0) semtcol->emplace_back(ngSemtOutputHandle->at(ihit));
    //
    if (fIsMC ) {
      const art::Ptr<recob::Hit> ahp = hitPtrMaker(outputHits->size() - 1);
      std::vector<art::Ptr<simb::MCParticle> > particle_vec = hittruth->at(hit.key());
      std::vector<anab::BackTrackerHitMatchingData const *> match_vec = hittruth->data(hit.key());
      for (size_t i_p = 0; i_p < particle_vec.size(); ++i_p) {
	outputHitPartAssns->addSingle(particle_vec[i_p],ahp,*match_vec[i_p]);
      }
    }

  }

  e.put(std::move(outputHits));
  if (fIsMC) e.put(std::move(outputHitPartAssns));
  if (fScoreCut>=0) e.put(std::move(semtcol), "semantic");
}

DEFINE_ART_MODULE(FilteredHitsProducer)
