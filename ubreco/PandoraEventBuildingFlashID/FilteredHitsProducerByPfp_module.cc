////////////////////////////////////////////////////////////////////////
// Class:       FilteredHitsProducerByPfp
// Plugin Type: producer (art v3_06_03)
// File:        FilteredHitsProducerByPfp_module.cc
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
#include "canvas/Persistency/Common/FindOneP.h"

#include <memory>

#include "canvas/Persistency/Common/Assns.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/MVAOutput.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

class FilteredHitsProducerByPfp;

class FilteredHitsProducerByPfp : public art::EDProducer {
public:
  explicit FilteredHitsProducerByPfp(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FilteredHitsProducerByPfp(FilteredHitsProducerByPfp const&) = delete;
  FilteredHitsProducerByPfp(FilteredHitsProducerByPfp&&) = delete;
  FilteredHitsProducerByPfp& operator=(FilteredHitsProducerByPfp const&) = delete;
  FilteredHitsProducerByPfp& operator=(FilteredHitsProducerByPfp&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.
  art::InputTag fHitLabel;
  art::InputTag fPfpLabel;
  art::InputTag fHitScoreLabel;
  art::InputTag fHitSemanticLabel;
  float fScoreCut;
  float fFracCut;
};


FilteredHitsProducerByPfp::FilteredHitsProducerByPfp(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fHitLabel(p.get<art::InputTag>("HitLabel")),
  fPfpLabel(p.get<art::InputTag>("PfpLabel")),
  fHitScoreLabel(p.get<art::InputTag>("HitScoreLabel")),
  fHitSemanticLabel(p.get<art::InputTag>("HitSemanticLabel")),
  fScoreCut(p.get<float>("ScoreCut")),
  fFracCut(p.get<float>("FracCut"))
{
  produces<std::vector<recob::Hit>>();
  produces<std::vector<anab::FeatureVector<5> > >("semantic");
  produces<std::vector<anab::FeatureVector<1> > >("filter");
}

void FilteredHitsProducerByPfp::produce(art::Event& e)
{
  // Implementation of required member function here.

  auto outputHits = std::make_unique<std::vector<recob::Hit> >();
  std::unique_ptr<std::vector<anab::FeatureVector<5>>> semtcol(new std::vector<anab::FeatureVector<5> >());
  std::unique_ptr<std::vector<anab::FeatureVector<1>>> filtcol(new std::vector<anab::FeatureVector<1> >());
  art::PtrMaker<recob::Hit> hitPtrMaker(e);

  art::Handle< std::vector< anab::FeatureVector<1> > > ngFiltOutputHandle;
  art::Handle< std::vector< anab::FeatureVector<5> > > ngSemtOutputHandle;
  e.getByLabel(fHitScoreLabel, ngFiltOutputHandle);
  e.getByLabel(fHitSemanticLabel, ngSemtOutputHandle);

  art::Handle< std::vector< recob::Hit > > hitListHandle;
  e.getByLabel(fHitLabel, hitListHandle);

  auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle> >(fPfpLabel);
  art::FindManyP<recob::Cluster> pfp_cluster_assn_v(pfp_h, e, fPfpLabel);
  art::FindManyP<recob::Vertex> pfp_vertex_assn_v(pfp_h, e, fPfpLabel);
  auto const& cluster_h = e.getValidHandle<std::vector<recob::Cluster> >(fPfpLabel);
  art::FindManyP<recob::Hit> cluster_hit_assn_v(cluster_h, e, fPfpLabel);

  // Step 1: Fill vector of all slice hits
  std::vector<art::Ptr<recob::Hit>> hit_slice_v;
  art::fill_ptr_vector(hit_slice_v,hitListHandle);

  // Step 2: Fill vector of clustered slice hits
  std::vector<std::vector<art::Ptr<recob::Hit>>> hit_pfp_v_v;
  std::vector<art::Ptr<recob::Hit>> hit_cluster_v;
  recob::tracking::Point_t nuvtx;
  std::vector<art::Ptr<recob::Vertex>> vtx_pfp_v;
  for (unsigned int p=0; p < pfp_h->size(); p++){
    const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, p);
    const std::vector< art::Ptr<recob::Vertex> > this_vertex_ptr_v = pfp_vertex_assn_v.at( pfp_ptr.key() );
    if (pfp_ptr->PdgCode()==12 || pfp_ptr->PdgCode()==14) {
      nuvtx = this_vertex_ptr_v[0]->position();
      continue;
    } else {
      if (this_vertex_ptr_v.size()>0) vtx_pfp_v.push_back(this_vertex_ptr_v[0]);
      else vtx_pfp_v.push_back(art::Ptr<recob::Vertex>());
    }
    std::vector<art::Ptr<recob::Hit>> hit_pfp_v;
    const std::vector< art::Ptr<recob::Cluster> > this_cluster_ptr_v = pfp_cluster_assn_v.at( pfp_ptr.key() );
    for (auto cluster_ptr : this_cluster_ptr_v) {
      const std::vector< art::Ptr<recob::Hit> > this_hit_ptr_v = cluster_hit_assn_v.at( cluster_ptr.key() );
      for (auto hit_ptr : this_hit_ptr_v) {
	hit_cluster_v.push_back(hit_ptr);
	hit_pfp_v.push_back(hit_ptr);
      }
    }
    hit_pfp_v_v.push_back(hit_pfp_v);
  }

  // Step 3: Sort both vectors
  sort(hit_slice_v.begin(), hit_slice_v.end());
  sort(hit_cluster_v.begin(), hit_cluster_v.end());

  // Step 4: Create a new vector to store the difference
  std::vector<art::Ptr<recob::Hit>> hit_unclustered_v;

  // Step 5: Use set_difference
  std::set_difference(hit_slice_v.begin(), hit_slice_v.end(), hit_cluster_v.begin(),
		      hit_cluster_v.end(), std::back_inserter(hit_unclustered_v));

  // Step 6: Filter unclustered hits
  for (unsigned int ih=0; ih<hit_unclustered_v.size(); ih++) {
    art::Ptr<recob::Hit> hit_ptr = hit_unclustered_v[ih];
    if (ngFiltOutputHandle->at(hit_ptr.key()).at(0)>=fScoreCut) {
      outputHits->emplace_back(*hit_ptr);
      semtcol->emplace_back(ngSemtOutputHandle->at(hit_ptr.key()));
      filtcol->emplace_back(ngFiltOutputHandle->at(hit_ptr.key()));
    }
  }

  // Step 7: Filter clustered hits based on the overall PFP filtering
  for (size_t j=0; j<hit_pfp_v_v.size(); j++) {
    auto hit_pfp_v = hit_pfp_v_v[j];
    unsigned int pass = 0, fail = 0;
    unsigned int hips = 0;
    if (hit_pfp_v.size()==0) continue;
    for (auto hit_ptr : hit_pfp_v) {
      if (ngFiltOutputHandle->at(hit_ptr.key()).at(0)>=fScoreCut) {
	pass++;
      } else {
	fail++;
      }
      if (ngSemtOutputHandle->at(hit_ptr.key()).at(1)>0.5) hips++;//fixme
    }
    // if (vtx_pfp_v[j].isNull()) std::cout << "pass=" << pass << " fail=" << fail << " hips=" << hips << std::endl;
    // else std::cout << "pass=" << pass << " fail=" << fail << " hips=" << hips
    // 	      << " dr=" << std::sqrt((vtx_pfp_v[j]->position()-nuvtx).Mag2())
    // 	      << " nuvtx=" << nuvtx << " vtx=" << vtx_pfp_v[j]->position()
    // 	      << std::endl;
    // keep if pass is above cut. Include also a protection to keep protons at vertex
    if ( pass < ((pass+fail)*fFracCut) && (vtx_pfp_v[j].isNull() || (hips < ((pass+fail)*fFracCut) || (vtx_pfp_v[j]->position()-nuvtx).Mag2()>1.0)) ) continue;
    for (auto hit_ptr : hit_pfp_v) {
      outputHits->emplace_back(*hit_ptr);
      semtcol->emplace_back(ngSemtOutputHandle->at(hit_ptr.key()));
      filtcol->emplace_back(ngFiltOutputHandle->at(hit_ptr.key()));
    }
  }

  std::cout << "FilteredHitProducerByPfp nhits in=" << hitListHandle->size() << " out=" << outputHits->size() << " diff=" << hitListHandle->size()-outputHits->size() << std::endl;
  e.put(std::move(outputHits));
  e.put(std::move(semtcol), "semantic");
  e.put(std::move(filtcol), "filter");
}

DEFINE_ART_MODULE(FilteredHitsProducerByPfp)
