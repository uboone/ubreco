////////////////////////////////////////////////////////////////////////
// Class:       OneTrackFilter
// Plugin Type: filter (art v3_01_02)
// File:        OneTrackFilter_module.cc
//
// Generated at Mon Feb 24 15:58:08 2020 by Matthew Rosenberg using cetskelgen
// Modified by C Thorpe July 2021
// from cetlib version v3_05_01.
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

#include "canvas/Persistency/Common/FindManyP.h"

#include "larcore/Geometry/Geometry.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"

#include "Pandora/PdgTable.h"

#include <memory>

class OneTrackFilter;


class OneTrackFilter : public art::EDFilter {
   public:
      explicit OneTrackFilter(fhicl::ParameterSet const& p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      OneTrackFilter(OneTrackFilter const&) = delete;
      OneTrackFilter(OneTrackFilter&&) = delete;
      OneTrackFilter& operator=(OneTrackFilter const&) = delete;
      OneTrackFilter& operator=(OneTrackFilter&&) = delete;

      // Required functions.
      bool filter(art::Event& e) override;

   private:

      // Declare member data here.
      std::string m_pfp_producer;
      std::string m_track_producer;

};


OneTrackFilter::OneTrackFilter(fhicl::ParameterSet const& p)
   : EDFilter{p},
   m_pfp_producer(p.get<std::string>("PFParticleLabel","pandora")),
   m_track_producer(p.get<std::string>("TrackLabel","pandora"))  
   // More initializers here.
{
   // Call appropriate produces<>() functions here.
   // Call appropriate consumes<>() for any products to be retrieved by this module.
}


bool OneTrackFilter::filter(art::Event& e)
{

   art::Handle<std::vector<recob::PFParticle>> pfParticleHandle;

   if(!e.getByLabel(m_pfp_producer, pfParticleHandle))
      throw cet::exception("OneTrackFilter") << "Repeated PFParticles in input!" << std::endl;

   // Fill vector of PFPs
   std::vector<art::Ptr<recob::PFParticle>> pfParticleVect; 
   art::fill_ptr_vector(pfParticleVect,pfParticleHandle);

   art::FindManyP<recob::Track> trackAssoc(pfParticleVect,e,m_track_producer);

   // Find the neutrino candidate
   size_t neutrinoID = 99999;

   for(const art::Ptr<recob::PFParticle> &pfp : pfParticleVect){

      if((pfp->IsPrimary() && (std::abs(pfp->PdgCode()) == 14 || std::abs(pfp->PdgCode()) == 12 ))){

         neutrinoID = pfp->Self();

      }   

   }

   // Find the longest track amongst its children

   int n_tracks=0;  

   for(const art::Ptr<recob::PFParticle> &pfp : pfParticleVect){

      if(pfp->Parent() != neutrinoID) continue;

      std::vector< art::Ptr<recob::Track> > pfpTracks = trackAssoc.at(pfp.key());
    
      if(pfpTracks.size()) n_tracks++;

   }

   return n_tracks > 0;
}


DEFINE_ART_MODULE(OneTrackFilter)
