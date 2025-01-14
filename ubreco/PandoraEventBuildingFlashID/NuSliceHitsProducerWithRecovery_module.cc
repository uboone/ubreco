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
  produces<std::vector<recob::PFParticle>>();
  produces<std::vector<recob::Vertex>>();
  produces<std::vector<recob::Cluster>>();
  produces<art::Assns<recob::PFParticle,recob::Vertex,void>>();
  produces<art::Assns<recob::PFParticle,recob::Cluster,void>>();
  produces<art::Assns<recob::Cluster,recob::Hit,void>>();

  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void NuSliceHitsProducerWithRecovery::produce(art::Event& e)
{
  // Implementation of required member function here.
  std::cout << "NuSliceHitsProducerWithRecovery --" << e.id() << std::endl;

  size_t nAssocHits = 0;
  auto outputHits = std::make_unique<std::vector<recob::Hit>>();
  auto outputHitPartAssns = std::make_unique<HitParticleAssociations>();
  auto outputPFP        = std::make_unique<std::vector<recob::PFParticle>>();
  auto outputVertex    = std::make_unique<std::vector<recob::Vertex>>();
  auto outPFPVtxAssns = std::make_unique<art::Assns<recob::PFParticle,recob::Vertex,void>>();
  auto outputCluster    = std::make_unique<std::vector<recob::Cluster>>();
  auto outPFPCluAssns = std::make_unique<art::Assns<recob::PFParticle,recob::Cluster,void>>();
  auto outCluHitAssns = std::make_unique<art::Assns<recob::Cluster,recob::Hit,void>>();

  art::PtrMaker<recob::Hit> hitPtrMaker(e);
  art::PtrMaker<recob::Cluster> cluPtrMaker(e);
  art::PtrMaker<recob::Vertex> vtxPtrMaker(e);
  art::PtrMaker<recob::PFParticle> pfpPtrMaker(e);

  auto const& inputSlice = e.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);
  auto const& inputPfp = e.getValidHandle<std::vector<recob::PFParticle>>(fPfpLabel);
  art::FindOneP<larpandoraobj::PFParticleMetadata> assocPfpMetadata(inputPfp, e, fPfpLabel);
  art::FindManyP<recob::PFParticle> assocSlicePfp(inputSlice, e, fPfpLabel);
  art::FindManyP<recob::Vertex> assocPfpVertex(inputPfp, e, fPfpLabel);
  art::FindManyP<recob::Cluster> assocPfpCluster(inputPfp, e, fPfpLabel);
  art::FindManyP<recob::Hit> assocSliceHit(inputSlice, e, fSliceLabel);

  auto const& cluster_h = e.getValidHandle<std::vector<recob::Cluster> >(fPfpLabel);
  art::FindManyP<recob::Hit> cluster_hit_assn_v(cluster_h, e, fPfpLabel);

  auto const& hitListHandle = e.getValidHandle<std::vector<recob::Hit>>(fHitLabel);
  std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> hittruth;
  if (!fHitTruthLabel.empty()) {
    hittruth = std::make_unique<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(hitListHandle, e, fHitTruthLabel);
  }

  std::vector<art::Ptr<recob::Slice> > slice_ptr_v;
  std::vector<art::Ptr<recob::PFParticle> > pfp_ptr_v;
  art::fill_ptr_vector(slice_ptr_v,inputSlice);
  art::fill_ptr_vector(pfp_ptr_v,inputPfp);
  std::map<int, art::Ptr<recob::PFParticle>> pfpmap;
  lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfp_ptr_v,pfpmap);

  //select the neutrino slice, or the one with largest NuScore if there is no neutino slice
  std::vector< art::Ptr<recob::Hit> > slice_hit_ptr_v;
  std::vector<art::Ptr<recob::Hit>> hit_cluster_v;
  size_t primarySelf = inputPfp->size();
  bool foundNuSlice = false;
  recob::Vertex::Point_t nuvtx;
  //
  for (auto& slice : slice_ptr_v) {
    const std::vector<art::Ptr<recob::PFParticle> >& nuSlicePFPs(assocSlicePfp.at(slice.key()));
    for (auto& pfp : nuSlicePFPs) {
      auto pfp_parent = lar_pandora::LArPandoraHelper::GetParentPFParticle(pfpmap,pfp);
      if (!lar_pandora::LArPandoraHelper::IsNeutrino(pfp_parent) ) continue;
      //
      if (pfp->IsPrimary()) {
	//
	primarySelf = pfp->Self();
	foundNuSlice = true;
	nuvtx = assocPfpVertex.at(pfp.key()).at(0)->position();
	slice_hit_ptr_v = assocSliceHit.at(slice.key());
	//
	auto pfParticleMetadata = assocPfpMetadata.at(pfp.key());
	auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
	std::cout << "ipfp self=" << pfp->Self() << " primary=" << pfp->IsPrimary() << " pdg=" << pfp->PdgCode() << std::endl;
	for (auto it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); it++) {
	  if (it->first == "NuScore") {
	    std::cout << "   meta=" << it->first << " " << it->second << std::endl;
	  }
	}
	//
	outputPFP->push_back(*pfp);
	art::Ptr<recob::PFParticle> pfp_ptr = pfpPtrMaker(outputPFP->size()-1);
	outputVertex->push_back(*assocPfpVertex.at(pfp.key()).at(0));
	auto vtx_ptr = vtxPtrMaker(outputVertex->size()-1);
	outPFPVtxAssns->addSingle(pfp_ptr,vtx_ptr);
	//
      } else {
	//std::cout << "ipfp self=" << pfp->Self() << " primary=" << pfp->IsPrimary() << " pdg=" << pfp->PdgCode() << " nvtx=" << assocPfpVertex.at(pfp.key()).size() << " nclu=" << assocPfpCluster.at(pfp.key()).size() << std::endl;
	outputPFP->push_back(*pfp);
	art::Ptr<recob::PFParticle> pfp_ptr = pfpPtrMaker(outputPFP->size()-1);
	if (assocPfpVertex.at(pfp.key()).size() > 0) {
	  outputVertex->push_back(*assocPfpVertex.at(pfp.key()).at(0));
	  auto vtx_ptr = vtxPtrMaker(outputVertex->size()-1);
	  outPFPVtxAssns->addSingle(pfp_ptr,vtx_ptr);
	}
	auto clusters = assocPfpCluster.at(pfp.key());	  
	for (auto clu : clusters) {
	  outputCluster->push_back(*clu);
	  auto clu_ptr = cluPtrMaker(outputCluster->size()-1);
	  outPFPCluAssns->addSingle(pfp_ptr,clu_ptr);
	  auto hits = cluster_hit_assn_v.at(clu.key());
	  for (auto hit : hits) {
	    //std::cout << "hit plane=" << hit->WireID().Plane << " wire=" << hit->WireID().Wire << " time=" << hit->PeakTime() << std::endl;
	    hit_cluster_v.push_back(hit);
	    outputHits->emplace_back(*hit);
	    const art::Ptr<recob::Hit> ahp = hitPtrMaker(outputHits->size() - 1);
	    outCluHitAssns->addSingle(clu_ptr,ahp);
	    //
	    if (!hittruth) continue;
	    std::vector<art::Ptr<simb::MCParticle>> particle_vec = hittruth->at(hit.key());
	    if (particle_vec.size()>0) nAssocHits++;
	    std::vector<anab::BackTrackerHitMatchingData const*> match_vec = hittruth->data(hit.key());
	    for (size_t i_p = 0; i_p < particle_vec.size(); ++i_p) {
	      outputHitPartAssns->addSingle(particle_vec[i_p], ahp, *match_vec[i_p]);
	    }
	  }
	}
      }
    }
  }
  // std::cout << "slice_hit_ptr_v size=" << slice_hit_ptr_v.size() << std::endl;
  // std::cout << "hits at before=" << outputHits->size() << std::endl;

  //get the slice with largest NuScore, from the "allOutcomes" output so that we can look at it under the neutrino interpretation
  auto const& slice_ao_h = e.getValidHandle<std::vector<recob::Slice> >(fAllOutcomesLabel);
  auto const& pfp_ao_h = e.getValidHandle<std::vector<recob::PFParticle> >(fAllOutcomesLabel);
  art::FindOneP<larpandoraobj::PFParticleMetadata> assocPfpMetadata_AO(pfp_ao_h, e, fAllOutcomesLabel);
  art::FindManyP<recob::PFParticle> assocSlicePfp_AO(slice_ao_h, e, fAllOutcomesLabel);
  art::FindManyP<recob::Vertex> pfp_vertex_ao_assn_v(pfp_ao_h, e, fAllOutcomesLabel);
  art::FindManyP<recob::Hit> assocSliceHit_AO(slice_ao_h, e, fAllOutcomesLabel);
  //
  std::vector<art::Ptr<recob::Slice> > slice_ao_ptr_v;
  std::vector<art::Ptr<recob::PFParticle> > pfp_ao_ptr_v;
  art::fill_ptr_vector(slice_ao_ptr_v,slice_ao_h);
  art::fill_ptr_vector(pfp_ao_ptr_v,pfp_ao_h);
  std::map<int, art::Ptr<recob::PFParticle>> pfpmap_ao;
  lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfp_ao_ptr_v,pfpmap_ao);
  //
  auto const& cluster_ao_h = e.getValidHandle<std::vector<recob::Cluster> >(fAllOutcomesLabel);
  art::FindManyP<recob::Cluster> pfp_cluster_ao_assn_v(pfp_ao_h, e, fAllOutcomesLabel);
  art::FindManyP<recob::Hit> cluster_hit_ao_assn_v(cluster_ao_h, e, fAllOutcomesLabel);
  //
  if (foundNuSlice==false && fRecoverHighestNuScoreSlice) {
    //
    float maxNuScore = std::numeric_limits<float>::lowest();
    std::vector<art::Ptr<recob::PFParticle> > bestSlicePFPs_AO;
    for (auto& slice : slice_ao_ptr_v) {
      const std::vector<art::Ptr<recob::PFParticle> >& nuSlicePFPs_AO(assocSlicePfp_AO.at(slice.key()));
      for (auto& pfp_AO : nuSlicePFPs_AO) {
	auto pfp_AO_parent = lar_pandora::LArPandoraHelper::GetParentPFParticle(pfpmap_ao,pfp_AO);
	if (!lar_pandora::LArPandoraHelper::IsNeutrino(pfp_AO_parent) ) continue;
	if (pfp_AO->IsPrimary()==false) continue;
	//
	auto pfParticleMetadata = assocPfpMetadata_AO.at(pfp_AO.key());
	auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
	for (auto it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); it++)  {
	  //std::cout << "ipfp=" << ipfp << " primary=" << pfp->IsPrimary() << " pdg=" << pfp->PdgCode() << " meta=" << it->first << " " << it->second << std::endl;
	  if (it->first != "NuScore") continue;
	  std::cout << "pfp  self=" << pfp_AO->Self() << " primary=" << pfp_AO->IsPrimary() << " pdg=" << pfp_AO->PdgCode() << " meta=" << it->first << " " << it->second << std::endl;
	  if (pfParticlePropertiesMap.at(it->first)>maxNuScore) {
	    primarySelf = pfp_AO->Self();
	    maxNuScore = pfParticlePropertiesMap.at(it->first);
	    nuvtx = pfp_vertex_ao_assn_v.at(pfp_AO.key()).at(0)->position();
	    slice_hit_ptr_v = assocSliceHit_AO.at(slice.key());
	    bestSlicePFPs_AO.clear();
	    for (auto& pfp : nuSlicePFPs_AO) {
	      auto pfp_parent = lar_pandora::LArPandoraHelper::GetParentPFParticle(pfpmap_ao,pfp);
	      if (!lar_pandora::LArPandoraHelper::IsNeutrino(pfp_parent) ) continue;
	      if (pfp->IsPrimary()) continue;
	      bestSlicePFPs_AO.push_back(pfp);
	    }
	  }
	}
      }
    }
    //
    for (auto pfp : bestSlicePFPs_AO) {
      outputPFP->push_back(*pfp);
      art::Ptr<recob::PFParticle> pfp_ptr = pfpPtrMaker(outputPFP->size()-1);
      if (pfp_vertex_ao_assn_v.at(pfp.key()).size() > 0) {
	outputVertex->push_back(*pfp_vertex_ao_assn_v.at(pfp.key()).at(0));
	auto vtx_ptr = vtxPtrMaker(outputVertex->size()-1);
	outPFPVtxAssns->addSingle(pfp_ptr,vtx_ptr);
      }
      auto clusters = pfp_cluster_ao_assn_v.at(pfp.key());	  
      for (auto clu : clusters) {
	outputCluster->push_back(*clu);
	auto clu_ptr = cluPtrMaker(outputCluster->size()-1);
	outPFPCluAssns->addSingle(pfp_ptr,clu_ptr);
	auto hits = cluster_hit_ao_assn_v.at(clu.key());
	for (auto hit : hits) {
	  //std::cout << "hit plane=" << hit->WireID().Plane << " wire=" << hit->WireID().Wire << " time=" << hit->PeakTime() << std::endl;
	  hit_cluster_v.push_back(hit);
	  outputHits->emplace_back(*hit);
	  const art::Ptr<recob::Hit> ahp = hitPtrMaker(outputHits->size() - 1);
	  outCluHitAssns->addSingle(clu_ptr,ahp);
	  //
	  if (!hittruth) continue;
	  std::vector<art::Ptr<simb::MCParticle>> particle_vec = hittruth->at(hit.key());
	  if (particle_vec.size()>0) nAssocHits++;
	  std::vector<anab::BackTrackerHitMatchingData const*> match_vec = hittruth->data(hit.key());
	  for (size_t i_p = 0; i_p < particle_vec.size(); ++i_p) {
	    outputHitPartAssns->addSingle(particle_vec[i_p], ahp, *match_vec[i_p]);
	  }
	}
      }
    }
    //
  }

 // std::cout << "slice_hit_ptr_v size=" << slice_hit_ptr_v.size() << std::endl;
  // std::cout << "hit_cluster_v size=" << hit_cluster_v.size() << std::endl;
  // std::cout << "output hits at before=" << outputHits->size() << std::endl;

  // now save the hits for this slice that are not clustered
  // Sort both vectors
  sort(slice_hit_ptr_v.begin(), slice_hit_ptr_v.end());
  sort(hit_cluster_v.begin(), hit_cluster_v.end());
  // Create a new vector to store the difference
  std::vector<art::Ptr<recob::Hit>> hit_unclustered_v;
  // Use set_difference
  std::set_difference(slice_hit_ptr_v.begin(), slice_hit_ptr_v.end(), hit_cluster_v.begin(),
		      hit_cluster_v.end(), std::back_inserter(hit_unclustered_v));
  // std::cout << "hit_unclustered_v size=" << hit_unclustered_v.size() << std::endl;
  for (size_t ihit = 0; ihit < hit_unclustered_v.size(); ++ihit) {
    auto hit = hit_unclustered_v.at(ihit);
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
  // std::cout << "output hits at after=" << outputHits->size() << std::endl;

  //now try to recover hits from potential 2nd showers in other slices
  if (fRecover2ndShower) {
    //
    for (auto& slice : slice_ao_ptr_v) {
      const std::vector<art::Ptr<recob::PFParticle> >& nuSlicePFPs_AO(assocSlicePfp_AO.at(slice.key()));
      for (auto& pfp_AO : nuSlicePFPs_AO) {
	auto pfp_AO_parent = lar_pandora::LArPandoraHelper::GetParentPFParticle(pfpmap_ao,pfp_AO);
	if (!lar_pandora::LArPandoraHelper::IsNeutrino(pfp_AO_parent) ) continue;
	if (pfp_AO->IsPrimary() || pfp_AO_parent->Self()==primarySelf) continue;
	std::cout << "self=" << pfp_AO->Self() << " " << " parent=" << pfp_AO->Parent() << " primary=" << primarySelf
		  << " parent2=" << lar_pandora::LArPandoraHelper::GetParentPFParticle(pfpmap_ao,pfp_AO)->Self() << std::endl;
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
	std::cout << "out-of-slice recovery: adding hits from this pfp pdg=" << pfp_AO->PdgCode() << " nhits=" << nhits
		  << " vtx dist=" << std::sqrt( (nuvtx-pfvtx).Mag2() ) << std::endl;
	//
	outputPFP->push_back(*pfp_AO);
	art::Ptr<recob::PFParticle> pfp_ptr = pfpPtrMaker(outputPFP->size()-1);
	outputVertex->push_back(*pfp_vertex_ao_assn_v.at(pfp_AO.key()).at(0));
	auto vtx_ptr = vtxPtrMaker(outputVertex->size()-1);
	outPFPVtxAssns->addSingle(pfp_ptr,vtx_ptr);
	//
	for (auto cluster_ptr : this_cluster_ptr_v) {
	  outputCluster->push_back(*cluster_ptr);
	  auto clu_ptr = cluPtrMaker(outputCluster->size()-1);
	  outPFPCluAssns->addSingle(pfp_ptr,clu_ptr);
	  const std::vector< art::Ptr<recob::Hit> > this_hit_ptr_v = cluster_hit_ao_assn_v.at( cluster_ptr.key() );
	  std::cout << "plane hits size=" << this_hit_ptr_v.size() << std::endl;
	  for (auto hit : this_hit_ptr_v) {
	    //std::cout << "hit plane=" << hit->WireID().Plane << " wire=" << hit->WireID().Wire << " time=" << hit->PeakTime() << std::endl;
	    outputHits->emplace_back(*hit);
	    const art::Ptr<recob::Hit> ahp = hitPtrMaker(outputHits->size() - 1);
	    outCluHitAssns->addSingle(clu_ptr,ahp);
	    //
	    if (!hittruth) continue;
	    std::vector<art::Ptr<simb::MCParticle>> particle_vec = hittruth->at(hit.key());
	    if (particle_vec.size()>0) nAssocHits++;
	    std::vector<anab::BackTrackerHitMatchingData const*> match_vec = hittruth->data(hit.key());
	    for (size_t i_p = 0; i_p < particle_vec.size(); ++i_p) {
	      outputHitPartAssns->addSingle(particle_vec[i_p], ahp, *match_vec[i_p]);
	    }
	  }
	}
      }
    }
    //
  }

  std::cout << "NuSliceHitProducer nhits=" << outputHits->size() << " nAssocHits=" << nAssocHits << " assns=" << outputHitPartAssns->size() << " foundNuSlice=" << foundNuSlice << std::endl;
  e.put(std::move(outputHits));
  if (!fHitTruthLabel.empty()) e.put(std::move(outputHitPartAssns));
  e.put(std::move(outputPFP));
  e.put(std::move(outputVertex));
  e.put(std::move(outPFPVtxAssns));
  e.put(std::move(outputCluster));
  e.put(std::move(outPFPCluAssns));
  e.put(std::move(outCluHitAssns));
  //
}

DEFINE_ART_MODULE(NuSliceHitsProducerWithRecovery)
