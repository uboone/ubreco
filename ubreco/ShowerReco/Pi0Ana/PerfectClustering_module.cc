////////////////////////////////////////////////////////////////////////
// Class:       PerfectClustering
// Plugin Type: producer (art v2_09_06)
// File:        PerfectClustering_module.cc
//
// Generated at Mon Mar 26 15:15:36 2018 by David Caratelli using cetskelgen
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

#include <memory>

// larsoft data-products
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "art/Persistency/Common/PtrMaker.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

// ROOT
#include <TTree.h>

class PerfectClustering;


class PerfectClustering : public art::EDProducer {
public:
  explicit PerfectClustering(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PerfectClustering(PerfectClustering const &) = delete;
  PerfectClustering(PerfectClustering &&) = delete;
  PerfectClustering & operator = (PerfectClustering const &) = delete;
  PerfectClustering & operator = (PerfectClustering &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  std::string fHitProducer;

  std::vector<size_t> AssociatedMCShowers(const size_t& h,
					  const art::FindManyP<simb::MCParticle>& hit_mcp_assn_v,
					  const std::map<size_t, std::vector<unsigned int> >& event_shower_map);

  TTree* _tree;
  double _shr_etot;
  double _shr_edep;
  double _pl2_integral;
  double _pl1_integral;
  double _pl0_integral;

};


PerfectClustering::PerfectClustering(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  
  fHitProducer = p.get<std::string>("HitProducer");

  produces<std::vector<recob::Cluster> >();
  produces<std::vector<recob::PFParticle> >();
  produces<art::Assns <recob::Cluster, recob::Hit> >();
  //produces<art::Assns <recob::Cluster, sim::MCShower> >();
  produces<art::Assns <recob::PFParticle, recob::Cluster> >();
  produces<art::Assns <recob::PFParticle, recob::Hit    > >();

}

void PerfectClustering::produce(art::Event & e)
{

  std::unique_ptr< std::vector<recob::Cluster> >                    Cluster_v                (new std::vector<recob::Cluster>                 );
  std::unique_ptr< art::Assns <recob::Cluster, recob::Hit> >        Cluster_Hit_assn_v       (new art::Assns<recob::Cluster,recob::Hit>       );
  //std::unique_ptr< art::Assns <recob::Cluster, sim::MCShower> >     Cluster_MCShower_assn_v  (new art::Assns<recob::Cluster,sim::MCShower>    );
  std::unique_ptr< std::vector<recob::PFParticle> >                 PFParticle_v             (new std::vector<recob::PFParticle>              );
  std::unique_ptr< art::Assns <recob::PFParticle, recob::Cluster> > PFParticle_Cluster_assn_v(new art::Assns<recob::PFParticle,recob::Cluster>);
  std::unique_ptr< art::Assns <recob::PFParticle, recob::Hit> >     PFParticle_Hit_assn_v    (new art::Assns<recob::PFParticle,recob::Hit>    );

  // load mcshowers & mctruth
  auto const& mcs_h = e.getValidHandle<std::vector<sim::MCShower>>("mcreco");
  auto const& mct_h = e.getValidHandle<std::vector<simb::MCTruth> >("generator");
  // load hits and mcparticles which are associated
  auto const& hit_h = e.getValidHandle<std::vector<recob::Hit> >(fHitProducer);
  // we will need the gaushits no matter what
  auto const& hit_g = e.getValidHandle<std::vector<recob::Hit> >("gaushit");
  // and associated mcparticles
  art::FindManyP<simb::MCParticle> hit_mcp_assn_v(hit_g, e, "gaushitTruthMatch");

  std::cout << "building hit to gaushit map " << std::endl;

  //if the hits are not gaushit, create a map connecting the hit index to the gaushit hit index.
  std::map<size_t,size_t> HitHitMap;
  if (fHitProducer != "gaushit") {
    for (size_t h0=0; h0 < hit_h->size(); h0++){
      auto const& hit0 = hit_h->at(h0);
      bool foundmatch = false;
      for (size_t h1=0; h1 < hit_g->size(); h1++){
	auto const& hit1 = hit_g->at(h1);
	if ( (hit0.PeakTime() == hit1.PeakTime()) && (hit1.WireID().Wire == hit0.WireID().Wire) ){
	  HitHitMap[h0] = h1;
	  foundmatch = true;
	  break;
	}
      }
      if (foundmatch == false)
	std::cout << "\t ERROR no match found!" << std::endl;
    }
  }
  else {
    for (size_t h0=0; h0 < hit_h->size(); h0++)
      HitHitMap[h0] = h0;
  }

  std::cout << "there are " << hit_h->size() << " hits and " << hit_mcp_assn_v.size() << " hit <-> mcp associations" << std::endl;
  std::cout << "there are " << HitHitMap.size() << " map elements" << std::endl;

  auto mct = mct_h->at(0);
  size_t npart = mct.NParticles();

  Double_t xyz[3] = {};

  // save the trackID of the pi0
  unsigned int pi0trkId = 0;
  int npi0 = 0;

  for (size_t i=0; i < npart; i++){
    auto const& part = mct.GetParticle(i);
    if ( (part.PdgCode() == 111) and (part.StatusCode() == 1) ){
      npi0 += 1;
      pi0trkId = part.TrackId();
      xyz[0] = part.Trajectory().X(0);
      xyz[1] = part.Trajectory().Y(0);
      xyz[2] = part.Trajectory().Z(0);
      break;
    }
  }

  std::cout << "pi0 track id : " << pi0trkId << std::endl;

  if (npi0 != 1){
    e.put(std::move(PFParticle_v));
    e.put(std::move(PFParticle_Cluster_assn_v));
    e.put(std::move(PFParticle_Hit_assn_v));
    e.put(std::move(Cluster_v));
    e.put(std::move(Cluster_Hit_assn_v));
    //e.put(std::move(Cluster_MCShower_assn_v));
    return;
  }

  // loop through MCShowers and identify those originating from the pi0
  // map connecting mcshower index to track ID vector for all e+/e- in MCShower
  std::map<size_t, std::vector<unsigned int> > event_shower_map;
  // map connecting e+/e- trackID in mcshower to mcshower index
  //std::map<unsigned int, size_t> event_mcpart_map;

  for (size_t i=0; i < mcs_h->size(); i++) {
    auto const& mcs = mcs_h->at(i);

    double x = mcs.Start().X();
    double y = mcs.Start().Y();
    double z = mcs.Start().Z();
    double d = sqrt( ( (xyz[0] - x) * (xyz[0] - x) ) +
		     ( (xyz[1] - y) * (xyz[1] - y) ) +
		     ( (xyz[2] - z) * (xyz[2] - z) ) );

    if (d < 0.01) {
      std::vector<unsigned int> shrtrackIDs = mcs.DaughterTrackID();
      shrtrackIDs.push_back( mcs.TrackID() );
      event_shower_map[ i ] = shrtrackIDs;
    }// if mcshower matched to pi0
  }// for all mcshowers
  std::cout << "found " << event_shower_map.size() << " mcshowers associated to the pi0" << std::endl;

  // save vector of clusters, one per shower, one per plane
  // for now only hit indices will go in
  // map connects MCShower index to vector of vector of clusters, one per plane.
  std::map< size_t, std::vector< std::vector<size_t> > > shower_cluster_v;
  for (auto const& showerinfo : event_shower_map) {
    shower_cluster_v[ showerinfo.first ] = std::vector<std::vector<size_t> >();
    shower_cluster_v[ showerinfo.first ].resize(3);
    for (size_t pl=0; pl < 3; pl++)
      shower_cluster_v[ showerinfo.first ].at(pl).clear();
  } 
  // shower -> mcshower index vector
  std::vector<size_t> shower_mcshower_v;

  // loop through hits, save them in a cluster if they are associated to one of the found showers
  for (size_t h=0; h < hit_h->size(); h++) {
    
    int plane = hit_h->at(h).WireID().Plane;

    // figure out which MCShower(s) this hit is associated to
    auto ass_mcs_v = AssociatedMCShowers(HitHitMap[h],hit_mcp_assn_v,event_shower_map);

    if (ass_mcs_v.size() == 1)
      shower_cluster_v[ass_mcs_v.at(0)][plane].push_back( h );
    
  }// for all mcparticles
  
  for (auto const& showerinfo : shower_cluster_v){
    for (size_t pl=0; pl < 3; pl++) 
      std::cout << "shower j on plane " << pl << " has " << showerinfo.second.at(pl).size() << " hits" << std::endl;
  }

  // now lets create some actual clusters
  auto const* geom = ::lar::providerFrom<geo::Geometry>();
  // cluster pointer maker for later to create associations
  art::PtrMaker<recob::Cluster> ClusPtrMaker(e, *this);
  // pfp pointer maker
  art::PtrMaker<recob::PFParticle> PFPPtrMaker(e, *this);

  int ctr = 0;
  for (auto const& showerinfo : shower_cluster_v) {
    //for (size_t j=0; j < shower_cluster_v.size(); j++){

    // count how many clusters added?
    // if 2+ then save a pfp as well
    int nclus = 0;
    // and index list of associated clusters
    std::vector<size_t> ass_clusters;

    for (size_t pl=0; pl < 3; pl++) {
      auto const& hit_idx_v = showerinfo.second.at(pl);
      if (hit_idx_v.size() == 0) continue;

      auto planeid = geo::PlaneID(0,0,pl);
      recob::Cluster clus(0., 0., 0., 0., 0., 0., 0., 
			  0., 0., 0., 0., 0., 0., 0., 
			  0., 0., 0., 0., 
			  hit_idx_v.size(), 0., 0., ctr*3+pl,
			  geom->View(planeid),
			  planeid);
      Cluster_v->emplace_back( clus );

      nclus += 1;

      art::Ptr<recob::Cluster> const ClusPtr = ClusPtrMaker(Cluster_v->size()-1);
      ass_clusters.push_back( Cluster_v->size() - 1);

      _shr_etot = mcs_h->at(showerinfo.first).Start().E();
      _shr_edep = mcs_h->at(showerinfo.first).DetProfile().E();
      _pl0_integral = 0.;
      _pl1_integral = 0.;
      _pl2_integral = 0.;

      //const art::Ptr<sim::MCShower> MCSPtr(mcs_h, showerinfo.first );
      //Cluster_MCShower_assn_v->addSingle(ClusPtr, MCSPtr);

      for (size_t h=0; h < showerinfo.second.at(pl).size(); h++) {
	const art::Ptr<recob::Hit> HitPtr(hit_h, showerinfo.second.at(pl).at(h) );
	Cluster_Hit_assn_v->addSingle(ClusPtr, HitPtr );
	if (HitPtr->WireID().Plane == 0)
	  _pl0_integral += HitPtr->Integral();
	if (HitPtr->WireID().Plane == 1)
	  _pl1_integral += HitPtr->Integral();
	if (HitPtr->WireID().Plane == 2)
	  _pl2_integral += HitPtr->Integral();
      }

    }// for all planes

    _tree->Fill();

    // if 2 or more clusters, save pfp as well
    if (nclus >= 2) {

      recob::PFParticle pfp(11,0,0,std::vector<size_t>());
      PFParticle_v->emplace_back(pfp);

      art::Ptr<recob::PFParticle> const PFPPtr = PFPPtrMaker(PFParticle_v->size()-1);

      // PFP -> Cluster associations
      for (auto const& clusidx : ass_clusters) {
	art::Ptr<recob::Cluster> const ClusPtr = ClusPtrMaker(clusidx);
	PFParticle_Cluster_assn_v->addSingle(PFPPtr,ClusPtr);
      }
      // PFP -> Hit associations
      for (size_t pl=0; pl < 3; pl++) {
	auto const& hit_idx_v = showerinfo.second.at(pl);
	for (auto const& hit_idx : hit_idx_v) {
	  const art::Ptr<recob::Hit> HitPtr(hit_h, hit_idx );
	  PFParticle_Hit_assn_v->addSingle(PFPPtr, HitPtr );
	}
      }

    }// if 2 or more clusters
    
    ctr += 1;
  }// for all MCShowers

  e.put(std::move(PFParticle_v));
  e.put(std::move(PFParticle_Cluster_assn_v));
  e.put(std::move(PFParticle_Hit_assn_v));
  e.put(std::move(Cluster_v));
  e.put(std::move(Cluster_Hit_assn_v));
  //e.put(std::move(Cluster_MCShower_assn_v));

}

std::vector<size_t> PerfectClustering::AssociatedMCShowers(const size_t& h,
							   const art::FindManyP<simb::MCParticle>& hit_mcp_assn_v,
							   const std::map<size_t, std::vector<unsigned int> >& event_shower_map) {
  
  std::vector<size_t> associated_mcshowers;

  // grab mcparticle associated to the hit
  auto ass_mcp_v = hit_mcp_assn_v.at(h);
  if (ass_mcp_v.size() == 0) 
    return associated_mcshowers;

  for (auto const& mcp : ass_mcp_v) {

    // does this trackID match any shower?
    bool matched = false;
    int nshrmatch = 0;
    size_t ctr = 0;
    for (auto const& mcsinfo : event_shower_map){
      for (auto const& trkid : mcsinfo.second)
	if ((unsigned int)mcp->TrackId() == trkid) {
	  matched = true;
	  nshrmatch = mcsinfo.first;
	  break;
	}// if we found a matching mcshower for the hit!
      ctr += 1;
    }// for all mcshowers associated to pi0

    if (matched == false) continue;

    // add to vector of associated mcshowers if this mcshower had not been added yet
    bool alreadyadded = false;
    for (auto const& mcshidx : associated_mcshowers)
      if (nshrmatch == (int)mcshidx) { alreadyadded = true; break; }
    
    if (alreadyadded == false) { associated_mcshowers.push_back( nshrmatch ); }
    
  }// for all associated mcparticles to hit    

  return associated_mcshowers;
}

void PerfectClustering::beginJob()
{

  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("_tree","Perfect clustering performance");
  _tree->Branch("_shr_etot",&_shr_etot,"shr_etot/D");
  _tree->Branch("_shr_edep",&_shr_edep,"shr_edep/D");
  _tree->Branch("_pl0_integral",&_pl0_integral,"pl0_integral/D");
  _tree->Branch("_pl1_integral",&_pl1_integral,"pl1_integral/D");
  _tree->Branch("_pl2_integral",&_pl2_integral,"pl2_integral/D");

}

void PerfectClustering::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(PerfectClustering)
