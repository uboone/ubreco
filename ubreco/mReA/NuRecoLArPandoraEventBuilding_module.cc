////////////////////////////////////////////////////////////////////////
// Class:       NuRecoLArPandoraEventBuilding
// Plugin Type: producer (art v3_01_02)
// File:        NuRecoLArPandoraEventBuilding_module.cc
//
// Generated at Wed Nov 18 19:59:00 2020 by Matthew Rosenberg using cetskelgen
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

#include "lardata/Utilities/AssociationUtil.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "larpandora/LArPandoraEventBuilding/Slice.h"
#include "larpandora/LArPandoraEventBuilding/SliceIdBaseTool.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraEvent.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

#include <memory>

class NuRecoLArPandoraEventBuilding;


class NuRecoLArPandoraEventBuilding : public art::EDProducer {
public:
  explicit NuRecoLArPandoraEventBuilding(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NuRecoLArPandoraEventBuilding(NuRecoLArPandoraEventBuilding const&) = delete;
  NuRecoLArPandoraEventBuilding(NuRecoLArPandoraEventBuilding&&) = delete;
  NuRecoLArPandoraEventBuilding& operator=(NuRecoLArPandoraEventBuilding const&) = delete;
  NuRecoLArPandoraEventBuilding& operator=(NuRecoLArPandoraEventBuilding&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.
  std::string   m_inputProducerLabel;  ///< Label for the Pandora instance that produced the collections we want to consolidate
  std::string   m_trackProducerLabel;  ///< Label for the track producer using the Pandora instance that produced the collections we want to consolidate
  std::string   m_showerProducerLabel; ///< Label for the shower producer using the Pandora instance that produced the collections we want to consolidate
  std::string   m_hitProducerLabel;    ///< Label for the hit producer that was used as input to the Pandora instance specified
  bool          m_shouldProduceT0s;    ///< If we should produce T0s (relevant when stitching over multiple drift volumes)
  art::InputTag m_pandoraTag;          ///< The input tag for the pandora producer

};


NuRecoLArPandoraEventBuilding::NuRecoLArPandoraEventBuilding(fhicl::ParameterSet const& p)
  : EDProducer{p},
  m_inputProducerLabel(p.get<std::string>("InputProducerLabel")),
  m_trackProducerLabel(p.get<std::string>("TrackProducerLabel")),
  m_showerProducerLabel(p.get<std::string>("ShowerProducerLabel")),
  m_hitProducerLabel(p.get<std::string>("HitProducerLabel")),
  m_shouldProduceT0s(p.get<bool>("ShouldProduceT0s")),
  m_pandoraTag(art::InputTag(m_inputProducerLabel))
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  produces< std::vector<recob::PFParticle> >();
  produces< std::vector<recob::SpacePoint> >();
  produces< std::vector<recob::Cluster> >();
  produces< std::vector<recob::Vertex> >();
  produces< std::vector<recob::Slice> >();
  produces< std::vector<recob::Track> >(); 
  produces< std::vector<recob::Shower> >();
  produces< std::vector<recob::PCAxis> >();
  produces< std::vector<larpandoraobj::PFParticleMetadata> >();

  produces< art::Assns<recob::PFParticle, recob::SpacePoint> >();
  produces< art::Assns<recob::PFParticle, recob::Cluster> >();
  produces< art::Assns<recob::PFParticle, recob::Vertex> >();
  produces< art::Assns<recob::PFParticle, recob::Slice> >();
  produces< art::Assns<recob::PFParticle, recob::Track> >();
  produces< art::Assns<recob::PFParticle, recob::Shower> >();
  produces< art::Assns<recob::PFParticle, recob::PCAxis> >();
  produces< art::Assns<recob::PFParticle, larpandoraobj::PFParticleMetadata> >();
  produces< art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta> >();
  produces< art::Assns<recob::Shower, recob::Hit> >();
  produces< art::Assns<recob::Shower, recob::PCAxis> >();
  produces< art::Assns<recob::SpacePoint, recob::Hit> >();
  produces< art::Assns<recob::Cluster, recob::Hit> >();
  produces< art::Assns<recob::Slice, recob::Hit> >();

  if(m_shouldProduceT0s){
    produces< std::vector<anab::T0> >();
    produces< art::Assns<recob::PFParticle, anab::T0> >();
  }
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}


void NuRecoLArPandoraEventBuilding::produce(art::Event& e)
{
  // Implementation of required member function here.
  art::Handle<std::vector<recob::PFParticle> > pfParticleHandle;
  e.getByLabel(m_pandoraTag, pfParticleHandle);

  lar_pandora::PFParticleVector particles;  
  for(unsigned int i = 0; i < pfParticleHandle->size(); ++i){
    const art::Ptr<recob::PFParticle> part(pfParticleHandle, i);
    particles.push_back(part);
  }

  const lar_pandora::LArPandoraEvent::Labels labels(m_inputProducerLabel, m_trackProducerLabel,
   m_showerProducerLabel, m_hitProducerLabel); 
  const lar_pandora::LArPandoraEvent consolidatedEvent(lar_pandora::LArPandoraEvent(this, &e,
   labels, m_shouldProduceT0s), particles);

  consolidatedEvent.WriteToEvent();
}


DEFINE_ART_MODULE(NuRecoLArPandoraEventBuilding)
