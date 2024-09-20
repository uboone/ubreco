////////////////////////////////////////////////////////////////////////
// Class:       RawDigitSmoother
// Plugin Type: producer (art v2_11_03)
// File:        RawDigitSmoother_module.cc
//
// UNDER CONSTRUCTION -- not yet functional as of Aug 2024.
//
// This module attempts to average together neighboring channels and creates
// a new collection of raw::Digits, which can then be passed through the 
// rest of the signal processing chain.
//
// W. Foreman
// Feb 2024
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft libraries
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "cetlib/search_path.h"

// MicroBooNE-specific includes
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

namespace raw {
  class RawDigitSmoother;
}

namespace raw{

  class RawDigitSmoother : public art::EDProducer 
  {
    public:
    explicit RawDigitSmoother(fhicl::ParameterSet const & p);
    RawDigitSmoother(RawDigitSmoother const &) = delete;
    RawDigitSmoother(RawDigitSmoother &&) = delete;
    RawDigitSmoother & operator = (RawDigitSmoother const &) = delete;
    RawDigitSmoother & operator = (RawDigitSmoother &&) = delete;

    // Required functions.
    void produce(art::Event & e) override;

    private:
    std::string fRawDigitInput;

  };
  
  RawDigitSmoother::RawDigitSmoother(fhicl::ParameterSet const & pset)
  {
    fRawDigitInput = pset.get<std::string>("daq");
    produces< std::vector<  raw::RawDigit > >();
  }

}



//###################################################
//  Main produce function
//###################################################
void raw::RawDigitSmoother::produce(art::Event & evt)
{
  std::cout
  <<"\n******************** RawDigitSmoother ****************************\n"
  <<  "Event "<<evt.id().event()<<" / run "<<evt.id().run()<<"\n";

  //================
  // Services
  //auto const& chanFilt = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
  
  //==========================
  // What we will produce
  std::unique_ptr< std::vector< raw::RawDigit > > RawDigitOutput_v(new std::vector<raw::RawDigit>);
  

  //==================================================
  // Get RawDigits from input file
  //art::Handle< std::vector<raw::RawDigit> > rdHandle;
  //std::vector< art::Ptr< raw::RawDigit> > rdVec;
  //if( evt.getByLabel("daq",rdHandle) ) 
  //  art::fill_ptr_vector(rdVec, rdHandle);
  //std::cout<<"Found "<<rdVec.size()<<" rawdigits\n";
 

  auto rdHandle = evt.getValidHandle<std::vector<raw::RawDigit>>("daq");
  std::cout<<"Found "<<rdHandle->size()<<" raw digits\n";   

  //if( chanFilt.Status(chan) < 4 ) chanisBad
  // chan<2 --> chan is dead?

  /*
  //===========================================
  // Make recob::SpacePoints out of the blip::Blips
  //===========================================
  for(size_t i=0; i<fBlipAlg->blips.size(); i++){
    auto& b = fBlipAlg->blips[i];
    
    Double32_t xyz[3];
    Double32_t xyz_err[6];
    Double32_t chiSquare = 0;
    Double32_t err = 0.; //b.MaxIntersectDiff?
    xyz[0]      = b.Position.X();
    xyz[1]      = b.Position.Y();
    xyz[2]      = b.Position.Z();
    xyz_err[0]  = err;
    xyz_err[1]  = err;
    xyz_err[2]  = err;
    xyz_err[3]  = err;
    xyz_err[4]  = err;
    xyz_err[5]  = err;
    
    recob::SpacePoint newpt(xyz,xyz_err,chiSquare);
    SpacePoint_v->emplace_back(newpt);
    
    // Hit associations 
    for(auto& hc : b.clusters ) {
      for(auto& ihit : hc.HitIDs ) {
        auto& hitptr = hitlist[ihit];
        util::CreateAssn(*this, evt, *SpacePoint_v, hitptr, *assn_hit_sps_v);
      }
    }
  
  }
    
  //===========================================
  // Put them on the event
  //===========================================
  evt.put(std::move(SpacePoint_v));
  evt.put(std::move(assn_hit_sps_v));
  */

  
  evt.put(std::move(RawDigitOutput_v));

}//END EVENT LOOP

DEFINE_ART_MODULE(raw::RawDigitSmoother)
