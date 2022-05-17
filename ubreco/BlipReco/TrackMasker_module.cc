////////////////////////////////////////////////////////////////////////
// Class:       TrackMasker
// Plugin Type: producer (art v2_05_01)
// File:        TrackMasker_module.cc
//
// W. Foreman
// March 2022
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes
#include <memory>
#include <cstring>
#include <vector>
#include <map>
#include <utility>
#include <iterator>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <typeinfo>
#include <cmath>

// Larsoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/Utilities/AssociationUtil.h"

// BlipReco-specific utility includes
#include "ubreco/BlipReco/Utils/BlipUtils.h"

class TrackMasker;


//###################################################
// Class Definition
//###################################################
class TrackMasker : public art::EDProducer {
public:
  explicit TrackMasker(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TrackMasker(TrackMasker const &) = delete;
  TrackMasker(TrackMasker &&) = delete;
  TrackMasker & operator = (TrackMasker const &) = delete;
  TrackMasker & operator = (TrackMasker &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;


private:
  
  // Producer labels
  std::string   fHitProducer;
  std::string   fTrkProducer;
  
  // Options
  double fMinTrkLength;
  double fMinTrkLengthVeto;
  double fVetoRadius;

  // Vector to hold hit indices to be removed
  std::vector<size_t> _vetohits;
  
  
  // TODO: remove these hard-coded conversion values
  float pitch        = 0.3;    // cm
  float samplePeriod = 0.5;    // us
  float driftSpeed   = 0.1041; // cm/us

};


//###################################################
// Constructor
//###################################################
TrackMasker::TrackMasker(fhicl::ParameterSet const & p)
{
  produces< std::vector< recob::Hit > >();
  
  fHitProducer            = p.get<std::string>("HitProducer");
  fTrkProducer            = p.get<std::string>("TrkProducer");
  fMinTrkLength           = p.get<double>     ("MinTrkLength");
  fMinTrkLengthVeto       = p.get<double>     ("MinTrkLengthVeto");
  fVetoRadius             = p.get<double>     ("VetoRadius");

  // min means min!
  fMinTrkLength = std::min(fMinTrkLength, fMinTrkLengthVeto);
}



//###################################################
// Main producer function
//###################################################
void TrackMasker::produce(art::Event & e)
{
  std::cout<<"\n"
  <<"=========== TrackMasker =========================\n"
  <<"Event "<<e.id().event()<<" / run "<<e.id().run()<<"\n";
  
  // *****************************
  // Grab data products from file
  // *****************************

  // grab tracks
  auto const& trk_h = e.getValidHandle<std::vector<recob::Track>>(fTrkProducer);
  std::vector<art::Ptr<recob::Track> > trklist;
  art::fill_ptr_vector(trklist,trk_h);
  
  // grab hits
  auto const& hit_h = e.getValidHandle<std::vector<recob::Hit>>(fHitProducer);
  std::vector<art::Ptr<recob::Hit> > hitlist;
  art::fill_ptr_vector(hitlist, hit_h);

  // grab hits associated to track
  art::FindManyP<recob::Track> hit_trk_assn_v(hit_h, e, fTrkProducer);
  
  // produce hits
  std::unique_ptr< std::vector<recob::Hit> > Hit_v(new std::vector<recob::Hit>);
  
  // clear track-like hit vector
  _vetohits.clear();

  std::cout<<"Found "<<hit_h->size()<<" hits from "<<fHitProducer<<"\n";
  std::cout<<"Found "<<trk_h->size()<<" tracks from "<<fTrkProducer<<"\n";
  
  
  
      

  //***************************************************
  // First go through tracks and designate which ones
  // pass our track length cuts by saving a map based
  // on each track's "ID".
  //***************************************************
  
  // Map of track ID to track length
  std::map<size_t,double> trklmap;
  
  for( size_t t=0; t < trklist.size(); t++ ) {
    
    auto const& trk = trklist[t];
    
    // must be some minimum length
    if (trk->Length() < fMinTrkLength) continue;

    // and track length must be < twice start-end distance
    if (trk->Length() > 2 * (trk->Vertex()-trk->End()).R() ) continue;

    trklmap[trk->ID()] = (double)trk->Length();
  }

  std::cout<<" --> "<<trklmap.size()<<" of those tracks pass trk length cut (> "<<fMinTrkLength<<" cm)\n";



  //***************************************************
  // Now categorize hits as 'tracked' or 'untracked'
  //***************************************************
  
  // boolean designating each hit as tracked
  std::vector<bool> hitIsVetoed( hitlist.size(), false);
  
  // map of (un-vetoed) hits to planes
  std::map<size_t,std::vector<size_t>> planehitmap;
  
  // 2D wire-time coordinates for each hit
  std::vector<TVector2> wtpoint( hitlist.size());

  // hits that are in tracks satisfying the veto region reqs
  std::vector<size_t> flaggedhits;

  for (size_t h=0; h < hitlist.size(); h++) {
    
    // find associated track
    auto const& trk_v = hit_trk_assn_v.at(h);
    if( trk_v.size() ) {
      
      size_t trkID = trk_v.at(0)->ID();
      if( trklmap.find(trkID) != trklmap.end() ) {
        hitIsVetoed[h] = true;
        _vetohits.push_back(h);
        if( trklmap[trkID] >= fMinTrkLengthVeto ){
          flaggedhits.push_back(h);
        }
      }
    
    }//endif track association exists
    
    // calculate wire-time coordinates of every hit so we aren't
    // repeating these calculations a bunch of times later on
    float w = hitlist[h]->WireID().Wire * pitch;
    float t = hitlist[h]->PeakTime() * samplePeriod * driftSpeed;
    wtpoint.at(h).Set(w,t);
    
    if( !hitIsVetoed[h] ) planehitmap[hitlist[h]->WireID().Plane].push_back(h);
  
  }//endloop over hits
  
  std::cout<<" --> removed "<<_vetohits.size()<<" tracked hits\n";
  std::cout<<" --> flagged "<<flaggedhits.size()<<" hits in tracks long enough to apply veto radius\n";
  size_t veto_hits_trk = _vetohits.size();
  
  //***************************************************
  // Loop through all tracked hits and for each one, 
  // check all un-tracked hits to determine if any
  // are within the veto radius.
  //**************************************************

  for( auto const& h : flaggedhits ) {
   
    for(auto const& hh : planehitmap[hitlist.at(h)->WireID().Plane] ) {

      // skip same hit (duh)
      //if( hh == h ) continue;
        
      // skip hits that are already vetoed
      if ( hitIsVetoed[hh] ) continue;
      
      float dw = fabs(wtpoint.at(hh).X()-wtpoint.at(h).X());
      float dt = fabs(wtpoint.at(hh).Y()-wtpoint.at(h).Y());
      
      // skip hits on far-away wires
      if( dw > fVetoRadius ) continue;

      // skip hits that are sufficiently separated in time
      if( dt > fVetoRadius ) continue;
      
      // finally, check 2D proximity
      if( sqrt( pow(dw,2)+pow(dt,2) ) > fVetoRadius ) continue;
      
      _vetohits.push_back(hh);
      hitIsVetoed[hh] = true;

    }
  }
  
  std::cout<<" --> vetoed additional "<<_vetohits.size()-veto_hits_trk<<" hits within radius\n";



  // ********************************************
  // finally, save hits that weren't marked veto
  // ********************************************
  for (size_t h=0; h < hit_h->size(); h++) {
    if (std::find(_vetohits.begin(),_vetohits.end(),h) == _vetohits.end() )
      Hit_v->emplace_back(hit_h->at(h));
  }

  
  printf("input hits  : %lu\n",hit_h->size());
  printf("output hits : %lu\n",Hit_v->size());

  e.put(std::move(Hit_v));
        
}

  
//###################################################
void TrackMasker::beginJob()
{
}


//###################################################
void TrackMasker::endJob()
{
}


DEFINE_ART_MODULE(TrackMasker)
