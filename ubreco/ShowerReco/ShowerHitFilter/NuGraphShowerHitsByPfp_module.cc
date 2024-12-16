////////////////////////////////////////////////////////////////////////
// Class:       PerfectClustering
// Plugin Type: producer (art v2_09_06)
// File:        NuGraphShowerHitsByPfp_module.cc
//
// Generated at Mon Aug 21 2024 by David Caratelli using cetskelgen
// from cetlib version v3_01_03.
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
#include "art_root_io/TFileService.h"

// larsoft data-products
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/AnalysisBase/MVAOutput.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"

#include "art/Persistency/Common/PtrMaker.h"

// ROOT
#include <TTree.h>

#include <vector>
#include <iostream>
#include <memory>
#include <algorithm>
#include <limits>
#include <cstdint>
#include <stdint.h>

class NuGraphShowerHitsByPfp;


class NuGraphShowerHitsByPfp : public art::EDProducer {
public:
  explicit NuGraphShowerHitsByPfp(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NuGraphShowerHitsByPfp(NuGraphShowerHitsByPfp const &) = delete;
  NuGraphShowerHitsByPfp(NuGraphShowerHitsByPfp &&) = delete;
  NuGraphShowerHitsByPfp & operator = (NuGraphShowerHitsByPfp const &) = delete;
  NuGraphShowerHitsByPfp & operator = (NuGraphShowerHitsByPfp &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  std::string fHitProducer;
  std::string fPandoraProducer;

  template <typename T, typename A>
  int arg_max(std::vector<T, A> const& vec)
  {
    return static_cast<int>(std::distance(vec.begin(), max_element(vec.begin(), vec.end())));
  }

  void AddDaughters(const art::Ptr<recob::PFParticle>& pfp_ptr,
                    const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h,
                    std::vector<art::Ptr<recob::PFParticle> > &pfp_v,
		    std::map<unsigned int, unsigned int>& pfpmap);

};


NuGraphShowerHitsByPfp::NuGraphShowerHitsByPfp(fhicl::ParameterSet const & p)
: EDProducer(p)
// Initialize member data here.
{
  
  fHitProducer = p.get<std::string>("HitProducer");
  fPandoraProducer = p.get<std::string>("PandoraProducer");

  produces<std::vector<recob::Hit> >();
  produces<std::vector<recob::PFParticle> >();

}

void NuGraphShowerHitsByPfp::produce(art::Event & e)
{

  std::unique_ptr< std::vector<recob::Hit> > Hit_v(new std::vector<recob::Hit>);
  auto outputPFP = std::make_unique<std::vector<recob::PFParticle>>();

  //load pfps and the associated slice (as done in NuSliceHitsProducer)
  //load hits associated to slice and to pfp clusters (as done in NuSliceHitsProducer, except we need to hits in clusters)
  //take all shower hits that are not clustered in pfps
  //take all hits in pfps that have a majority of hits labeled as shower
  //keep track of the PFPs that are removed

  auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle> >(fPandoraProducer);
  art::FindManyP<recob::Cluster> pfp_cluster_assn_v(pfp_h, e, fPandoraProducer);
  art::FindManyP<recob::Slice>   pfp_slice_assn_v(pfp_h, e, fPandoraProducer);
  auto const& cluster_h = e.getValidHandle<std::vector<recob::Cluster> >(fPandoraProducer);
  art::FindManyP<recob::Hit> cluster_hit_assn_v(cluster_h, e, fPandoraProducer);
  auto const& slice_h = e.getValidHandle<std::vector<recob::Slice> >(fPandoraProducer);
  art::FindManyP<recob::Hit> slice_hit_assn_v(slice_h, e, fPandoraProducer);

  std::map<unsigned int, unsigned int> pfpmap;
  for (unsigned int p=0; p < pfp_h->size(); p++) pfpmap[pfp_h->at(p).Self()] = p;

  std::vector<art::Ptr<recob::PFParticle> > pfp_ptr_v;
  std::vector<std::vector<art::Ptr<recob::Hit>>> hit_pfp_v_v;
  std::vector<art::Ptr<recob::Hit>> hit_cluster_v;
  std::vector<art::Ptr<recob::Hit>> hit_slice_v;
  for (unsigned int p=0; p < pfp_h->size(); p++){

    const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, p);

    // start from primary neutrino PFParticles
    //std::cout << "pfp pdg=" << pfp_ptr->PdgCode() << " primary=" << pfp_ptr->IsPrimary()<< std::endl;
    if (pfp_ptr->IsPrimary() == false) continue;
    if (pfp_ptr->PdgCode()!=12 && pfp_ptr->PdgCode()!=14) {
      //std::cout << "to be removed" << std::endl;
      outputPFP->push_back(*pfp_ptr);
      continue;
    }

    const std::vector< art::Ptr<recob::Slice> > this_slice_ptr_v = pfp_slice_assn_v.at( pfp_ptr.key() );
    for (auto slice_ptr : this_slice_ptr_v) {
      const std::vector< art::Ptr<recob::Hit> > this_hit_ptr_v = slice_hit_assn_v.at( slice_ptr.key() );
      for (auto hit_ptr : this_hit_ptr_v) {
	hit_slice_v.push_back(hit_ptr);
      }
    }

    AddDaughters(pfp_ptr, pfp_h, pfp_ptr_v,pfpmap);

    for (unsigned int q=0; q < pfp_ptr_v.size(); q++) {
      std::vector<art::Ptr<recob::Hit>> hit_pfp_v;
      const std::vector< art::Ptr<recob::Cluster> > this_cluster_ptr_v = pfp_cluster_assn_v.at( pfp_ptr_v[q].key() );
      for (auto cluster_ptr : this_cluster_ptr_v) {
	const std::vector< art::Ptr<recob::Hit> > this_hit_ptr_v = cluster_hit_assn_v.at( cluster_ptr.key() );
	for (auto hit_ptr : this_hit_ptr_v) {
	  hit_cluster_v.push_back(hit_ptr);
	  hit_pfp_v.push_back(hit_ptr);
	}
      }
      hit_pfp_v_v.push_back(hit_pfp_v);
    }

  } // pfp loop

   // Step 3: Sort both vectors
  sort(hit_slice_v.begin(), hit_slice_v.end());
  sort(hit_cluster_v.begin(), hit_cluster_v.end());

  // Step 4: Create a new vector to store the difference
  std::vector<art::Ptr<recob::Hit>> hit_unclustered_v;

    // Step 5: Use set_difference
  std::set_difference(hit_slice_v.begin(), hit_slice_v.end(), hit_cluster_v.begin(),
		      hit_cluster_v.end(), std::back_inserter(hit_unclustered_v));

  // load NuGraph hit semantic scores
  auto const& hitsWithScores = proxy::getCollection<std::vector<recob::Hit>>(e,art::InputTag(fHitProducer), //tag of the hit collection we ran the GNN on
									     proxy::withParallelData<anab::FeatureVector<5>>(art::InputTag(fHitProducer, "semantic")));

  for (unsigned int ih=0; ih<hit_unclustered_v.size(); ih++) {
    art::Ptr<recob::Hit> hit_ptr = hit_unclustered_v[ih];
    auto scores = hitsWithScores[hit_ptr.key()].get<anab::FeatureVector<5>>();    
    std::vector<float> ng2semscores;
    for (size_t i=0;i<scores.size();i++) ng2semscores.push_back(scores[i]);
    unsigned int sem_label = arg_max(ng2semscores);
    if (sem_label == 2){
      Hit_v->emplace_back(*hit_ptr);
    }
  }

  for (size_t j=0; j<hit_pfp_v_v.size(); j++) {
    auto hit_pfp_v = hit_pfp_v_v[j];
    unsigned int lbl[5] = {0,0,0,0,0};
    //std::cout << "pfp pdg=" << pfp_ptr_v[j]->PdgCode() << " primary=" << pfp_ptr_v[j]->IsPrimary() << " nhits="<< hit_pfp_v.size()<< std::endl;
    if (hit_pfp_v.size()==0) continue;
    for (auto hit_ptr : hit_pfp_v) {
      auto scores = hitsWithScores[hit_ptr.key()].get<anab::FeatureVector<5>>();    
      std::vector<float> ng2semscores;
      for (size_t i=0;i<scores.size();i++) ng2semscores.push_back(scores[i]);
      unsigned int sem_label = arg_max(ng2semscores);
      lbl[sem_label]++;
    }
    if (lbl[2]<lbl[0] || lbl[2]<lbl[1] || lbl[2]<lbl[3] || lbl[2]<lbl[4]) continue;
    for (auto hit_ptr : hit_pfp_v) {
      Hit_v->emplace_back(*hit_ptr);
    }
    //std::cout << "to be removed" << std::endl;
    outputPFP->push_back(*pfp_ptr_v[j]);
  }

  e.put(std::move(Hit_v));
  e.put(std::move(outputPFP));
  
}

void NuGraphShowerHitsByPfp::AddDaughters(const art::Ptr<recob::PFParticle>& pfp_ptr,  const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h, std::vector<art::Ptr<recob::PFParticle> > &pfp_v,std::map<unsigned int, unsigned int>& pfpmap) {

  auto daughters = pfp_ptr->Daughters();

  pfp_v.push_back(pfp_ptr);

  for(auto const& daughterid : daughters) {

    if (pfpmap.find(daughterid) == pfpmap.end()) {
      std::cout << "Did not find DAUGHTERID in map! error"<< std::endl;
      continue;
    }

    const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, pfpmap.at(daughterid) );

    AddDaughters(pfp_ptr, pfp_h, pfp_v, pfpmap);

  }// for all daughters

  return;
}

void NuGraphShowerHitsByPfp::beginJob()
{

}

void NuGraphShowerHitsByPfp::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(NuGraphShowerHitsByPfp)
