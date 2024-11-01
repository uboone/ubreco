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

#include <memory>

#include "canvas/Persistency/Common/Assns.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/MVAOutput.h"
// #include "lardata/Utilities/AssociationUtil.h"
// #include "lardata/RecoBaseProxy/ProxyBase.h"

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
  art::InputTag fMatchHitLabel;
  art::InputTag fPfpLabel;
  art::InputTag fHitScoreLabel;
  art::InputTag fHitSemanticLabel;
  float fScoreCut;
  float fFracCut;

  void AddDaughters(const art::Ptr<recob::PFParticle>& pfp_ptr,
                    const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h,
                    std::vector<art::Ptr<recob::PFParticle> > &pfp_v,
		    std::map<unsigned int, unsigned int>& pfpmap);

};


FilteredHitsProducerByPfp::FilteredHitsProducerByPfp(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fHitLabel(p.get<art::InputTag>("HitLabel")),
  fMatchHitLabel(p.get<art::InputTag>("MatchHitLabel")),
  fPfpLabel(p.get<art::InputTag>("PfpLabel")),
  fHitScoreLabel(p.get<art::InputTag>("HitScoreLabel")),
  fHitSemanticLabel(p.get<art::InputTag>("HitSemanticLabel")),
  fScoreCut(p.get<float>("ScoreCut")),
  fFracCut(p.get<float>("FracCut"))
{
  produces<std::vector<recob::Hit>>();
  produces<std::vector<anab::FeatureVector<5> > >("semantic");
}

void FilteredHitsProducerByPfp::produce(art::Event& e)
{
  // Implementation of required member function here.

  auto outputHits = std::make_unique<std::vector<recob::Hit> >();
  std::unique_ptr<std::vector<anab::FeatureVector<5>>> semtcol(new std::vector<anab::FeatureVector<5> >());
  art::PtrMaker<recob::Hit> hitPtrMaker(e);

  art::Handle< std::vector< anab::FeatureVector<1> > > ngFiltOutputHandle;
  art::Handle< std::vector< anab::FeatureVector<5> > > ngSemtOutputHandle;
  e.getByLabel(fHitScoreLabel, ngFiltOutputHandle);
  e.getByLabel(fHitSemanticLabel, ngSemtOutputHandle);

  art::Handle< std::vector< recob::Hit > > hitListHandle;
  e.getByLabel(fHitLabel, hitListHandle);
  art::Handle< std::vector< recob::Hit > > mtchHitListHandle;
  e.getByLabel(fMatchHitLabel, mtchHitListHandle);
  std::map<size_t,size_t> hmap;
  for (size_t imtchhit=0; imtchhit<mtchHitListHandle->size();imtchhit++) {
    art::Ptr<recob::Hit> mtchhit(mtchHitListHandle,imtchhit);
    for (size_t ihit=0; ihit<hitListHandle->size();ihit++) {
      art::Ptr<recob::Hit> hit(hitListHandle,ihit);
      if (hit->WireID().Plane==mtchhit->WireID().Plane &&
	  hit->WireID().Wire==mtchhit->WireID().Wire &&
	  std::fabs(hit->PeakTime()-mtchhit->PeakTime())<0.0001 &&
	  std::fabs(hit->RMS()-mtchhit->RMS())<0.0001) {
	hmap[mtchhit.key()] = hit.key();
	break;
      }
    }
  }

  ///
  auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle> >(fPfpLabel);
  art::FindManyP<recob::Cluster> pfp_cluster_assn_v(pfp_h, e, fPfpLabel);
  art::FindManyP<recob::Slice>   pfp_slice_assn_v(pfp_h, e, fPfpLabel);
  auto const& cluster_h = e.getValidHandle<std::vector<recob::Cluster> >(fPfpLabel);
  art::FindManyP<recob::Hit> cluster_hit_assn_v(cluster_h, e, fPfpLabel);
  auto const& slice_h = e.getValidHandle<std::vector<recob::Slice> >(fPfpLabel);
  art::FindManyP<recob::Hit> slice_hit_assn_v(slice_h, e, fPfpLabel);
  auto assocPfpMetadata = e.getValidHandle<std::vector<larpandoraobj::PFParticleMetadata>>(fPfpLabel);

  std::cout << "N slices=" << slice_h->size() << std::endl;

  std::map<unsigned int, unsigned int> pfpmap;
  for (unsigned int p=0; p < pfp_h->size(); p++) pfpmap[pfp_h->at(p).Self()] = p;

  std::vector<art::Ptr<recob::PFParticle> > pfp_ptr_v;
  std::vector<std::vector<art::Ptr<recob::Hit>>> hit_pfp_v_v;
  std::vector<art::Ptr<recob::Hit>> hit_cluster_v;
  std::vector<art::Ptr<recob::Hit>> hit_slice_v;
  bool foundNuSlice = false;
  for (unsigned int p=0; p < pfp_h->size(); p++){

    const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, p);

    // start from primary neutrino PFParticles
    if (pfp_ptr->IsPrimary() == false) continue;
    if (pfp_ptr->PdgCode()!=12 && pfp_ptr->PdgCode()!=14) continue;
    foundNuSlice = true;

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

  // std::cout << __FILE__ << " " << __LINE__ << std::endl;
  // std::cout << "hit_slice_v.size()=" << hit_slice_v.size() << std::endl;
  // std::cout << "hit_cluster_v.size()=" << hit_cluster_v.size() << std::endl;

  if (!foundNuSlice) {
    //get the slice with largest NuScore, but don't bother looking at PFPs as they are from the cosmic reco path for now
    size_t primaryIdx = pfp_h->size();
    float maxNuScore = std::numeric_limits<float>::lowest();
    for (unsigned int p=0; p < pfp_h->size(); p++){

      const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, p);
      if (pfp_ptr->IsPrimary() == false) continue;
    
      auto pfParticleMetadata = assocPfpMetadata->at(p);
      auto pfParticlePropertiesMap = pfParticleMetadata.GetPropertiesMap();
      if (!pfParticlePropertiesMap.empty()) {
	auto it = pfParticlePropertiesMap.begin();
	while (it != pfParticlePropertiesMap.end()) {
	  //std::cout << "ipfp=" << p << " primary=" << pfp_ptr->IsPrimary() << " pdg=" << pfp_ptr->PdgCode() << " meta=" << it->first << " " << it->second << std::endl;
	  if (it->first == "NuScore") {
	    std::cout << "ipfp=" << p << " primary=" << pfp_ptr->IsPrimary() << " pdg=" << pfp_ptr->PdgCode() << " meta=" << it->first << " " << it->second << std::endl;
	    if (pfParticlePropertiesMap.at(it->first)>maxNuScore) {
	      primaryIdx = p;
	      maxNuScore = pfParticlePropertiesMap.at(it->first);
	    }
	  }
	  it++;
	}
      } // if PFP metadata exists!
    }
    if (primaryIdx < pfp_h->size()) {
      const std::vector< art::Ptr<recob::Slice> > this_slice_ptr_v = pfp_slice_assn_v.at( primaryIdx );
      for (auto slice_ptr : this_slice_ptr_v) {
	const std::vector< art::Ptr<recob::Hit> > this_hit_ptr_v = slice_hit_assn_v.at( slice_ptr.key() );
	for (auto hit_ptr : this_hit_ptr_v) {
	  hit_slice_v.push_back(hit_ptr);
	}
      }
    }
  }

  // std::cout << __FILE__ << " " << __LINE__ << std::endl;
  // std::cout << "hit_slice_v.size()=" << hit_slice_v.size() << std::endl;
  // std::cout << "hit_cluster_v.size()=" << hit_cluster_v.size() << std::endl;

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
    if (ngFiltOutputHandle->at(hmap[hit_ptr.key()]).at(0)>=fScoreCut) {
      outputHits->emplace_back(*hit_ptr);
      semtcol->emplace_back(ngSemtOutputHandle->at(hmap[hit_ptr.key()]));
    }
  }

  // Step 7: Filter clustered hits based on the overall PFP filtering
  for (size_t j=0; j<hit_pfp_v_v.size(); j++) {
    auto hit_pfp_v = hit_pfp_v_v[j];
    unsigned int pass = 0, fail = 0;
    //std::cout << "pfp pdg=" << pfp_ptr_v[j]->PdgCode() << " primary=" << pfp_ptr_v[j]->IsPrimary() << " nhits="<< hit_pfp_v.size()<< std::endl;
    if (hit_pfp_v.size()==0) continue;
    for (auto hit_ptr : hit_pfp_v) {
      if (ngFiltOutputHandle->at(hmap[hit_ptr.key()]).at(0)>=fScoreCut) {
	pass++;
      } else {
	fail++;
      }
    }
    if ( pass < ((pass+fail)*fFracCut) ) continue;
    for (auto hit_ptr : hit_pfp_v) {
      outputHits->emplace_back(*hit_ptr);
      semtcol->emplace_back(ngSemtOutputHandle->at(hmap[hit_ptr.key()]));
    }
  }

  std::cout << "FilteredHitProducerByPfp nhits=" << outputHits->size() << std::endl;
  e.put(std::move(outputHits));
  e.put(std::move(semtcol), "semantic");
}

void FilteredHitsProducerByPfp::AddDaughters(const art::Ptr<recob::PFParticle>& pfp_ptr,  const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h, std::vector<art::Ptr<recob::PFParticle> > &pfp_v,std::map<unsigned int, unsigned int>& pfpmap) {

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

DEFINE_ART_MODULE(FilteredHitsProducerByPfp)
