////////////////////////////////////////////////////////////////////////
// Class:       NuSliceHitsProducerWithRecovery
// Plugin Type: producer (art v3_06_03)
// File:        NuSliceHitsProducerWithRecovery_module.cc
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
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "nusimdata/SimulationBase/MCParticle.h"

class NuSliceHitsProducerWithRecovery;

using HitParticleAssociations =
  art::Assns<simb::MCParticle, recob::Hit, anab::BackTrackerHitMatchingData>;

class NuSliceHitsProducerWithRecovery : public art::EDProducer {
public:
  explicit NuSliceHitsProducerWithRecovery(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NuSliceHitsProducerWithRecovery(NuSliceHitsProducerWithRecovery const&) = delete;
  NuSliceHitsProducerWithRecovery(NuSliceHitsProducerWithRecovery&&) = delete;
  NuSliceHitsProducerWithRecovery& operator=(NuSliceHitsProducerWithRecovery const&) = delete;
  NuSliceHitsProducerWithRecovery& operator=(NuSliceHitsProducerWithRecovery&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:
  // Declare member data here.
  art::InputTag fPfpLabel;
  art::InputTag fSliceLabel;
  art::InputTag fHitLabel;
  art::InputTag fHitTruthLabel;
  art::InputTag fAllOutcomesLabel;
  bool fRecoverHighestNuScoreSlice;
  bool fRecover2ndShower;
  float fVtxDistCut2;
  int fMaxHitCut;
};

NuSliceHitsProducerWithRecovery::NuSliceHitsProducerWithRecovery(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fPfpLabel(p.get<art::InputTag>("PfpLabel", "pandora"))
  , fSliceLabel(p.get<art::InputTag>("SliceLabel", "pandora"))
  , fHitLabel(p.get<art::InputTag>("HitLabel", "gaushit"))
  , fHitTruthLabel(p.get<art::InputTag>("HitTruthLabel", ""))
  , fAllOutcomesLabel(p.get<art::InputTag>("AllOutcomesLabel", art::InputTag("pandoraPatRec","allOutcomes")))
  , fRecoverHighestNuScoreSlice(p.get<bool>("RecoverHighestNuScoreSlice"))
  , fRecover2ndShower(p.get<bool>("Recover2ndShower"))
  , fVtxDistCut2(p.get<float>("VtxDistCut"))
  , fMaxHitCut(p.get<int>("MaxHitCut"))
// More initializers here.
{
  fVtxDistCut2 *= fVtxDistCut2;

  // Call appropriate produces<>() functions here.
  produces<std::vector<recob::Hit>>();
  if (!fHitTruthLabel.empty()) produces<HitParticleAssociations>();

  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void NuSliceHitsProducerWithRecovery::produce(art::Event& e)
{
  // Implementation of required member function here.
  std::cout << "NuSliceHitsProducerWithRecovery --" << e.id() << std::endl;

  size_t nAssocHits = 0;
  auto outputHits = std::make_unique<std::vector<recob::Hit>>();
  auto outputHitPartAssns = std::make_unique<HitParticleAssociations>();
  art::PtrMaker<recob::Hit> hitPtrMaker(e);

  art::ValidHandle<std::vector<recob::PFParticle>> inputPfp =
    e.getValidHandle<std::vector<recob::PFParticle>>(fPfpLabel);
  auto assocPfpSlice = std::unique_ptr<art::FindManyP<recob::Slice>>(
    new art::FindManyP<recob::Slice>(inputPfp, e, fPfpLabel));
  auto assocPfpMetadata = e.getValidHandle<std::vector<larpandoraobj::PFParticleMetadata>>(fPfpLabel);
  auto assocPfpVertex = std::unique_ptr<art::FindManyP<recob::Vertex>>(
    new art::FindManyP<recob::Vertex>(inputPfp, e, fPfpLabel));
  auto assocPfpCluster = std::unique_ptr<art::FindManyP<recob::Cluster>>(
    new art::FindManyP<recob::Cluster>(inputPfp, e, fPfpLabel));

  art::ValidHandle<std::vector<recob::Slice>> inputSlice =
    e.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);
  auto assocSliceHit = std::unique_ptr<art::FindManyP<recob::Hit>>(
    new art::FindManyP<recob::Hit>(inputSlice, e, fSliceLabel));

  auto const& cluster_h = e.getValidHandle<std::vector<recob::Cluster> >(fPfpLabel);
  art::FindManyP<recob::Hit> cluster_hit_assn_v(cluster_h, e, fPfpLabel);

  art::Handle<std::vector<recob::Hit>> hitListHandle;
  e.getByLabel(fHitLabel, hitListHandle);
  std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> hittruth;
  if (!fHitTruthLabel.empty()) {
    hittruth = std::make_unique<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(
      hitListHandle, e, fHitTruthLabel);
  }

  //select the neutrino slice, or the one with largest NuScore if there is no neutino slice
  size_t primaryIdx = inputPfp->size();
  size_t primarySelf = inputPfp->size();
  float maxNuScore = std::numeric_limits<float>::lowest();
  bool foundNuSlice = false;
  recob::Vertex::Point_t nuvtx;
  for (size_t ipfp = 0; ipfp < inputPfp->size(); ipfp++) {
    art::Ptr<recob::PFParticle> pfp(inputPfp, ipfp);
    if (pfp->IsPrimary() == false) continue;
    //
    auto pfParticleMetadata = assocPfpMetadata->at(ipfp);
    auto pfParticlePropertiesMap = pfParticleMetadata.GetPropertiesMap();
    //
    auto PDG = fabs(pfp->PdgCode());
    if (PDG == 12 || PDG == 14) {
      //
      std::cout << "ipfp=" << ipfp << " self=" << pfp->Self() << " primary=" << pfp->IsPrimary() << " pdg=" << pfp->PdgCode() << std::endl;
      for (auto it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); it++) {
	if (it->first == "NuScore") {
	  std::cout << "   meta=" << it->first << " " << it->second << std::endl;
	}
      }
      //
      primaryIdx = ipfp;
      primarySelf = pfp->Self();
      foundNuSlice = true;
      nuvtx = assocPfpVertex->at(pfp.key()).at(0)->position();
      break;
    }
    //
    if (fRecoverHighestNuScoreSlice==false) continue;

    if (!pfParticlePropertiesMap.empty()) {
      auto it = pfParticlePropertiesMap.begin();
      while (it != pfParticlePropertiesMap.end()) {
	//std::cout << "ipfp=" << ipfp << " primary=" << pfp->IsPrimary() << " pdg=" << pfp->PdgCode() << " meta=" << it->first << " " << it->second << std::endl;
	if (it->first == "NuScore") {
	  std::cout << "ipfp=" << ipfp << " self=" << pfp->Self() << " primary=" << pfp->IsPrimary() << " pdg=" << pfp->PdgCode() << " meta=" << it->first << " " << it->second << std::endl;
	  if (pfParticlePropertiesMap.at(it->first)>maxNuScore) {
	    primaryIdx = ipfp;
	    primarySelf = pfp->Self();
	    maxNuScore = pfParticlePropertiesMap.at(it->first);
	    nuvtx = assocPfpVertex->at(pfp.key()).at(0)->position();//fixme: this is not the "neutrino" vertex since this is a cosmic slice
	  }
	}
	it++;
      }
    } // if PFP metadata exists!
  }

  //now save the hits for this slice
  if (primaryIdx < inputPfp->size()) {
    art::Ptr<recob::PFParticle> pfp(inputPfp, primaryIdx);

    auto assocSlice = assocPfpSlice->at(pfp.key());
    auto sliceHits = assocSliceHit->at(assocSlice[0].key());

    for (size_t ihit = 0; ihit < sliceHits.size(); ++ihit) {
      auto hit = sliceHits.at(ihit);
      outputHits->emplace_back(*hit);

      if (!hittruth) continue;
      std::vector<art::Ptr<simb::MCParticle>> particle_vec = hittruth->at(hit.key());
      if (particle_vec.size()>0) nAssocHits++;
      std::vector<anab::BackTrackerHitMatchingData const*> match_vec = hittruth->data(hit.key());
      const art::Ptr<recob::Hit> ahp = hitPtrMaker(outputHits->size() - 1);
      for (size_t i_p = 0; i_p < particle_vec.size(); ++i_p) {
        outputHitPartAssns->addSingle(particle_vec[i_p], ahp, *match_vec[i_p]);
      }
    }

    //now try to recover hits from potential 2nd showers in other slices
    if (fRecover2ndShower) {
      //
      //get the slice with largest NuScore, from the "allOutcomes" output so that we can look at it under the neutrino interpretation
      auto const& slice_ao_h = e.getValidHandle<std::vector<recob::Slice> >(fAllOutcomesLabel);
      auto const& pfp_ao_h = e.getValidHandle<std::vector<recob::PFParticle> >(fAllOutcomesLabel);
      auto const& cluster_ao_h = e.getValidHandle<std::vector<recob::Cluster> >(fAllOutcomesLabel);
      art::FindOneP<larpandoraobj::PFParticleMetadata> assocPfpMetadata_AO(pfp_ao_h, e, fAllOutcomesLabel);
      art::FindManyP<recob::PFParticle> assocSlicePfp_AO(slice_ao_h, e, fAllOutcomesLabel);
      art::FindManyP<recob::Cluster> pfp_cluster_ao_assn_v(pfp_ao_h, e, fAllOutcomesLabel);
      art::FindManyP<recob::Vertex> pfp_vertex_ao_assn_v(pfp_ao_h, e, fAllOutcomesLabel);
      art::FindManyP<recob::Hit> cluster_hit_ao_assn_v(cluster_ao_h, e, fAllOutcomesLabel);
      //
      //std::cout << "assocSlicePfp.size()=" << assocSlicePfp_AO.size() << std::endl;
      std::vector<art::Ptr<recob::Slice> > slice_ao_ptr_v;
      std::vector<art::Ptr<recob::PFParticle> > pfp_ao_ptr_v;
      art::fill_ptr_vector(slice_ao_ptr_v,slice_ao_h);
      art::fill_ptr_vector(pfp_ao_ptr_v,pfp_ao_h);
      std::map<int, art::Ptr<recob::PFParticle>> pfpmap_ao;
      lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfp_ao_ptr_v,pfpmap_ao);
      //
      for (auto& slice : slice_ao_ptr_v) {
	const std::vector<art::Ptr<recob::PFParticle> >& nuSlicePFPs_AO(assocSlicePfp_AO.at(slice.key()));
	for (auto& pfp_AO : nuSlicePFPs_AO) {
	  auto pfp_AO_parent = lar_pandora::LArPandoraHelper::GetParentPFParticle(pfpmap_ao,pfp_AO);
	  if (!lar_pandora::LArPandoraHelper::IsNeutrino(pfp_AO_parent) ) continue;
	  if (pfp_AO->IsPrimary() || pfp_AO_parent->Self()==primarySelf) continue;
	  //std::cout << "self=" << pfp_AO->Self() << " " << " parent=" << pfp_AO->Parent() << " primary=" << primarySelf << " parent2=" << lar_pandora::LArPandoraHelper::GetParentPFParticle(pfpmap_ao,pfp_AO)->Self() << std::endl;
	  // only pfps within fVtxDistCut of the neutrino vertex
	  if (pfp_vertex_ao_assn_v.at(pfp_AO.key()).size()==0) continue;
	  auto pfvtx = pfp_vertex_ao_assn_v.at(pfp_AO.key()).at(0)->position();
	  //std::cout << "pfp vtx dist=" << std::sqrt( (nuvtx-pfvtx).Mag2() ) << std::endl;
	  if ( (nuvtx-pfvtx).Mag2() > fVtxDistCut2 ) continue;
	  //
	  const std::vector< art::Ptr<recob::Cluster> > this_cluster_ptr_v = pfp_cluster_ao_assn_v.at( pfp_AO.key() );
	  //
	  int nhits = 0;
	  for (auto cluster_ptr : this_cluster_ptr_v) nhits += cluster_hit_ao_assn_v.at( cluster_ptr.key() ).size();
	  //
	  //consider only showers, or tracks with a limited number of hits consistent with a pi0 shower
	  if (pfp_AO->PdgCode()==13 && nhits>fMaxHitCut) continue;
	  //
	  std::cout << "out-of-slice recovery: adding hits from this pfp pdg=" << pfp_AO->PdgCode() << " nhits=" << nhits << " vtx dist=" << std::sqrt( (nuvtx-pfvtx).Mag2() ) << std::endl;
	  for (auto cluster_ptr : this_cluster_ptr_v) {
	    const std::vector< art::Ptr<recob::Hit> > this_hit_ptr_v = cluster_hit_ao_assn_v.at( cluster_ptr.key() );
	    std::cout << "plane hits size=" << this_hit_ptr_v.size() << std::endl;
	    for (auto hit : this_hit_ptr_v) {
	      //std::cout << "hit plane=" << hit->WireID().Plane << " wire=" << hit->WireID().Wire << " time=" << hit->PeakTime() << std::endl;
	      outputHits->emplace_back(*hit);
	      //
	      if (!hittruth) continue;
	      std::vector<art::Ptr<simb::MCParticle>> particle_vec = hittruth->at(hit.key());
	      if (particle_vec.size()>0) nAssocHits++;
	      std::vector<anab::BackTrackerHitMatchingData const*> match_vec = hittruth->data(hit.key());
	      const art::Ptr<recob::Hit> ahp = hitPtrMaker(outputHits->size() - 1);
	      for (size_t i_p = 0; i_p < particle_vec.size(); ++i_p) {
		outputHitPartAssns->addSingle(particle_vec[i_p], ahp, *match_vec[i_p]);
	      }
	    }
	  }
	}
      }
    }
  }

  std::cout << "NuSliceHitProducer nhits=" << outputHits->size() << " nAssocHits=" << nAssocHits << " assns=" << outputHitPartAssns->size() << " foundNuSlice=" << foundNuSlice << std::endl;
  e.put(std::move(outputHits));
  if (!fHitTruthLabel.empty()) e.put(std::move(outputHitPartAssns));
}

DEFINE_ART_MODULE(NuSliceHitsProducerWithRecovery)
