////////////////////////////////////////////////////////////////////////
// Class:       MuonRemoval
// Plugin Type: producer (art v3_01_02)
// File:        MuonRemoval_module.cc
//
// Generated at Wed Mar  4 09:54:37 2020 by Matthew Rosenberg using cetskelgen
// from cetlib version v3_05_01.
// Modified by C Thorpe July 2021
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

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "ubreco/mReA/helpers/FindParticleHits.h"

#include "Pandora/PdgTable.h"

#include <utility>
#include <algorithm>
#include <memory>

class MuonRemoval;


class MuonRemoval : public art::EDProducer {
public:
  explicit MuonRemoval(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MuonRemoval(MuonRemoval const&) = delete;
  MuonRemoval(MuonRemoval&&) = delete;
  MuonRemoval& operator=(MuonRemoval const&) = delete;
  MuonRemoval& operator=(MuonRemoval&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.
  std::string m_pfp_producer;
  std::string m_track_producer;
  //std::string m_muon_producer;
  std::string m_slice_producer;
  std::string m_hit_producer;
  std::string m_cluster_producer;
  std::string m_pfpToCls_producer;
  std::string m_clsToHit_producer;
  std::string m_dataReco1_process;
  std::string m_simReco1_process;
  std::string m_mcp_producer;
  std::string m_backTrack_producer;
  std::string m_mcp_truthProcess;
  std::string m_hit_truthProcess;
  std::string m_backTrack_truthProcess;

  bool m_removeMCMuHits;
  bool m_useNuSliceOnly;

  mReA::FindParticleHits hitFinder;

};


MuonRemoval::MuonRemoval(fhicl::ParameterSet const& p)
  : EDProducer{p},
  m_pfp_producer(p.get<std::string>("PFParticleLabel","pandora")),
  m_track_producer(p.get<std::string>("TrackLabel","pandora")),
  //m_muon_producer(p.get<std::string>("MuonFilterLabel","NuCCproducer")),
  m_slice_producer(p.get<std::string>("SliceLabel","pandora")),
  m_hit_producer(p.get<std::string>("HitLabel", "gaushit")),
  m_cluster_producer(p.get<std::string>("ClusterLabel","pandora")),
  m_pfpToCls_producer(p.get<std::string>("PFPToClusterAssnLabel","pandora")),
  m_clsToHit_producer(p.get<std::string>("ClusterToHitAssnLabel","pandora")),
  m_dataReco1_process(p.get<std::string>("DataReco1Process", "DataRecoStage1Test")),
  m_simReco1_process(p.get<std::string>("SimReco1Process", "mReARecoStage1a")),
  m_mcp_producer(p.get<std::string>("MCParticleLabel","largeant")),
  m_backTrack_producer(p.get<std::string>("BackTrackerLabel","gaushitTruthMatch")),
  m_mcp_truthProcess(p.get<std::string>("MCPartTruthProcess", "G4EDep")),
  m_hit_truthProcess(p.get<std::string>("HitTruthProcess", "OverlayStage1a")),
  m_backTrack_truthProcess(p.get<std::string>("BackTrackerTruthProcess", "OverlayRecoStage1b")),
  m_removeMCMuHits(p.get<bool>("RemoveMCMuHits", false)),
  m_useNuSliceOnly(p.get<bool>("UseNeutrinoSliceOnly", true))
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  produces< std::vector<recob::Hit> >();
  produces< art::Assns<recob::Wire,recob::Hit> >();
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  hitFinder.Initialize(m_pfp_producer, m_cluster_producer, m_hit_producer,
   m_pfpToCls_producer, m_clsToHit_producer, m_mcp_producer, m_backTrack_producer,
   m_hit_truthProcess, m_mcp_truthProcess, m_backTrack_truthProcess);
}


void MuonRemoval::produce(art::Event& e)
{
  // Implementation of required member function here.

   // Use a different muon ID - find the longest track (for now)

   // Get reco muon PFParticle
   art::Handle< std::vector<recob::PFParticle> > pfParticleHandle;
   e.getByLabel(m_pfp_producer, pfParticleHandle);

   // Fill vector of PFPs
   std::vector< art::Ptr<recob::PFParticle>>pfParticleVect; 
   art::fill_ptr_vector(pfParticleVect,pfParticleHandle);

   // Get the associated tracks
   art::FindManyP<recob::Track> trackAssoc(pfParticleVect,e,m_track_producer);

   // Find the neutrino candidate
   size_t neutrinoID = 99999;
   for(const art::Ptr<recob::PFParticle> &pfp : pfParticleVect){

      if(pfp->IsPrimary() && (std::abs(pfp->PdgCode()) == 14 || std::abs(pfp->PdgCode()) == 12))
         neutrinoID = pfp->Self();

   }

   // Find the longest track amongst its children

   int recoMuKey=-1;
   double l_max=0.0;

   for(const art::Ptr<recob::PFParticle> &pfp : pfParticleVect){

      if(pfp->Parent() != neutrinoID) continue;

      std::vector< art::Ptr<recob::Track> > pfpTracks = trackAssoc.at(pfp.key());

      for(const art::Ptr<recob::Track> &trk : pfpTracks){

         if(trk->Length() > l_max){ 
            recoMuKey =  pfp->Self(); 
            l_max = trk->Length(); 
         }
      }  
   }

  if(recoMuKey == -1) throw cet::exception("MuonReplacement") <<
    "No tracks in primary hierarchy" << std::endl;

  std::map<size_t, art::Ptr<recob::PFParticle> > pfParticleIdMap;

  for(unsigned int iPF = 0; iPF < pfParticleHandle->size(); ++iPF){
    const art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle, iPF);
    if(!pfParticleIdMap.emplace(pfParticle->Self(), pfParticle).second)
      throw cet::exception("MuonRemoval") << "Repeated PFParticles in input!" << std::endl;
  }

  const art::Ptr<recob::PFParticle>& recoMuon = pfParticleIdMap[recoMuKey];

  //Get hits associated with reco muon
  std::vector< art::Ptr<recob::Hit> > recoMuonHits;
  std::vector< art::Ptr<recob::Hit> > mcMuonHits;
  hitFinder.getPFPHits(pfParticleIdMap, e, recoMuon, recoMuonHits);
  unsigned int nMuonHits = recoMuonHits.size();

  //get hits associated with MC muon for MC input, if requested
  if(m_removeMCMuHits){
    art::Handle< std::vector<simb::MCParticle> > mcParticleHandle;
    art::InputTag mcpTag(m_mcp_producer, "", m_mcp_truthProcess);
    e.getByLabel(mcpTag, mcParticleHandle);
    int mcMuonTrackId = -1;

    std::map<int, art::Ptr<simb::MCParticle> > mcParticleIdMap;
    for (unsigned int i = 0; i < mcParticleHandle->size(); ++i){
      const art::Ptr<simb::MCParticle> mcParticle(mcParticleHandle, i);
      if(!mcParticleIdMap.emplace(mcParticle->TrackId(), mcParticle).second)
        throw cet::exception("MuonRemoval") << "Repeated MCParticles in input!" << std::endl;
      if(std::abs(mcParticle->PdgCode()) == pandora::MU_MINUS && mcParticle->Mother() == 0)
        mcMuonTrackId = mcParticle->TrackId();
    }

    if(mcMuonTrackId < 0) throw cet::exception("MuonRemoval")
      << "Couldn't find primary MC muon. RemoveMCMuHits should only be set "
      << "to true for true CC nu_mu MC input" <<  std::endl;
    const art::Ptr<simb::MCParticle>& mcMuon = mcParticleIdMap.at(mcMuonTrackId);

    hitFinder.getMCPHits(mcParticleIdMap, e, mcMuon, mcMuonHits);
    nMuonHits = mcMuonHits.size();
  }

  //Get neutrino slice
  const art::Ptr<recob::PFParticle>& recoNu = pfParticleIdMap[recoMuon->Parent()];
  art::FindManyP<recob::Slice> pfpToSliceAssoc(pfParticleHandle, e, m_slice_producer);
  const auto nuSlices = pfpToSliceAssoc.at( recoNu.key() );
  if(nuSlices.size() != 1) throw cet::exception("MuonRemoval") <<
   "Could not find a single slice associated with the reconstructed neutrino" << std::endl;
  const art::Ptr<recob::Slice>& nuSlice = nuSlices[0];

  //create output hit collection
  auto outputHits = std::make_unique<std::vector<recob::Hit> >();
  auto outputWireHitAssn = std::make_unique< art::Assns<recob::Wire, recob::Hit> >();

  //add data hits from neutrino slice only (by default), excluding muon hits 
  art::Handle<std::vector<recob::Hit>> HitHandle;
  e.getByLabel(m_hit_producer, HitHandle);
  art::FindManyP<recob::Slice> hitToSliceAssoc(HitHandle, e, m_slice_producer);
  art::FindManyP<recob::Wire> dataHitToWireAssoc(HitHandle, e,m_hit_producer);
  unsigned int removed_hit_count_slice = 0;
  unsigned int removed_hit_count_mu = 0;

  for(unsigned int iH = 0; iH < HitHandle->size(); ++iH){

    const art::Ptr<recob::Hit> hitPtr(HitHandle,iH);
    //we'll need to dereference hitPtr later, so to be safe...
    if(hitPtr.isNull()) throw cet::exception("MuonRemoval") <<
     "null hit Ptr in collection from module: " << m_hit_producer << std::endl;

    if(m_useNuSliceOnly){
      const auto associatedSlices = hitToSliceAssoc.at( hitPtr.key() );
      if(associatedSlices.size() == 0){
        mf::LogWarning("MuonRemoval") << "Encountered hit from original event "
         << "that is not associated with any slices. Will not write hit to output." << std::endl;
        ++removed_hit_count_slice;
        continue;
      }
      if(associatedSlices.size() > 1) throw cet::exception("MuonRemoval")
       << "Hit matched to multiple slices." << std::endl;
      if(associatedSlices[0] != nuSlice){ ++removed_hit_count_slice;  continue; }
    }

    if(m_removeMCMuHits){
      auto muHitItr = std::find(mcMuonHits.begin(), mcMuonHits.end(), hitPtr);
      if(muHitItr != mcMuonHits.end()){
        mcMuonHits.erase(muHitItr); //speed up the next find call a bit
        ++removed_hit_count_mu;
        continue;
      }
    }
    else{
      auto muHitItr = std::find(recoMuonHits.begin(), recoMuonHits.end(), hitPtr);
      if(muHitItr != recoMuonHits.end()){
        recoMuonHits.erase(muHitItr); //speed up the next find call a bit
        ++removed_hit_count_mu;
        continue;
      }
    }

    const auto associatedWires = dataHitToWireAssoc.at( hitPtr.key() );
    if(associatedWires.size() != 1) throw cet::exception("MuonRemoval")
     << "Could not match hit to a single wire." << std::endl;

    const recob::Hit& hitRef = *hitPtr;
    recob::Hit outHit(hitRef);
    outputHits->emplace_back(std::move(outHit));
    outputWireHitAssn->addSingle(associatedWires[0], hitPtr);
  }


  mf::LogDebug("MuonRemoval") << "Wrote " << outputHits->size() <<
   " data hits to mReA hit collection. Of original hit count (" << HitHandle->size() <<
   "), " << removed_hit_count_slice << " hits were removed by neutrino slice cut and " <<
   removed_hit_count_mu << " of " << nMuonHits << " muon hits were removed." << std::endl;

  e.put(std::move(outputHits));
  e.put(std::move(outputWireHitAssn));

}

DEFINE_ART_MODULE(MuonRemoval)
