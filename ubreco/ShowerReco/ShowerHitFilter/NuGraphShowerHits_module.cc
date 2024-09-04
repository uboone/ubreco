////////////////////////////////////////////////////////////////////////
// Class:       PerfectClustering
// Plugin Type: producer (art v2_09_06)
// File:        NuGraphShowerHits_module.cc
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
#include "art/Framework/Services/Optional/TFileService.h"

// larsoft data-products
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/AnalysisBase/MVAOutput.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"

#include "art/Persistency/Common/PtrMaker.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

// ROOT
#include <TTree.h>

#include <vector>
#include <iostream>
#include <memory>
#include <algorithm>
#include <limits>
#include <cstdint>
#include <stdint.h>

class NuGraphShowerHits;


class NuGraphShowerHits : public art::EDProducer {
public:
  explicit NuGraphShowerHits(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NuGraphShowerHits(NuGraphShowerHits const &) = delete;
  NuGraphShowerHits(NuGraphShowerHits &&) = delete;
  NuGraphShowerHits & operator = (NuGraphShowerHits const &) = delete;
  NuGraphShowerHits & operator = (NuGraphShowerHits &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  std::string fHitProducer;

  template <typename T, typename A>
  int arg_max(std::vector<T, A> const& vec)
  {
    return static_cast<int>(std::distance(vec.begin(), max_element(vec.begin(), vec.end())));
  }

};


NuGraphShowerHits::NuGraphShowerHits(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  
  fHitProducer = p.get<std::string>("HitProducer");

  produces<std::vector<recob::Hit> >();

}

void NuGraphShowerHits::produce(art::Event & e)
{

  std::unique_ptr< std::vector<recob::Hit> > Hit_v(new std::vector<recob::Hit>);

  // load existing hits
  auto const& hit_h = e.getValidHandle<std::vector<recob::Hit> >(fHitProducer);
  std::cout << "there are " << hit_h->size() << " hits in the input file. " << std::endl;

  // load NuGraph hit semantic scores
  auto const& hitsWithScores = proxy::getCollection<std::vector<recob::Hit>>(
									     e,
									     art::InputTag(fHitProducer), //tag of the hit collection we ran the GNN on
									     //proxy::withParallelData<anab::FeatureVector<1>>(art::InputTag("NuGraph", "filter")),
									     proxy::withParallelData<anab::FeatureVector<5>>(art::InputTag(fHitProducer, "semantic")));

  std::vector<int> ng2semclucounts(5,0);
  std::vector<int> ng2semslccounts(5,0);
  for (auto& h : hitsWithScores) {
    auto scores = h.get<anab::FeatureVector<5>>();
    std::vector<float> ng2semscores;
    for (size_t i=0;i<scores.size();i++) ng2semscores.push_back(scores[i]);
    unsigned int sem_label = arg_max(ng2semscores);
    ng2semslccounts[sem_label]++;

    //std::cout << "[NuGraph] hit key is " << h.key() << std::endl;
    std::cout << "[NuGraph] hit index is " << h.index() << std::endl;
    std::cout << "[NuGraph] hit label is " << sem_label << std::endl;

    if (sem_label == 2){
      recob::Hit hit = hit_h->at(h.index());
      Hit_v->emplace_back(hit);
    }
  }

  
  //loop through hits, save them in a cluster if they are associated to one of the found showers
  //  for (size_t h=0; h < hit_h->size(); h++) {
    /*
    auto scores = hitsWithScores[h].get<anab::FeatureVector<5>>();
    std::vector<float> ng2semscores;
    for (size_t i=0;i<scores.size();i++) ng2semscores.push_back(scores[i]);
    unsigned int sem_label = arg_max(ng2semscores);
    ng2sempfpcounts[sem_label]++;
    ng2semclucounts[sem_label]++;
    */
  //  }// for all hits
 
  e.put(std::move(Hit_v));
  
}

void NuGraphShowerHits::beginJob()
{

}

void NuGraphShowerHits::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(NuGraphShowerHits)
