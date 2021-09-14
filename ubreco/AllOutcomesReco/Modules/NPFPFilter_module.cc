////////////////////////////////////////////////////////////////////////
// Class:       NPFPFilter
// Plugin Type: filter (art v3_01_02)
// File:        NPFPFilter_module.cc
//
// Generated at Tue Sep 14 08:44:52 2021 by Christopher Thorpe using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

#include <memory>

#include "lardataobj/RecoBase/PFParticle.h"

namespace alloutcomes {
   class NPFPFilter;
}


class alloutcomes::NPFPFilter : public art::EDFilter {
   public:
      explicit NPFPFilter(fhicl::ParameterSet const& p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      NPFPFilter(NPFPFilter const&) = delete;
      NPFPFilter(NPFPFilter&&) = delete;
      NPFPFilter& operator=(NPFPFilter const&) = delete;
      NPFPFilter& operator=(NPFPFilter&&) = delete;

      // Required functions.
      bool filter(art::Event& e) override;

   private:

      // fhicl parameters
      std::string f_PFPLabel; 
      int f_MinPFPs;

};


alloutcomes::NPFPFilter::NPFPFilter(fhicl::ParameterSet const& p)
   : EDFilter{p},  // ,
   f_PFPLabel(p.get<std::string>("PFPLabel","pandora")),
   f_MinPFPs(p.get<int>("MinPFPs",3)) 
{
   std::cout << "PFPLabel = " << f_PFPLabel << std::endl;
   std::cout << "Removing events with fewer than = " << f_MinPFPs << " PFPs in the primary hierarchy" << std::endl;

   if(f_MinPFPs < 1)
      throw cet::exception("NPFPFilter") << "You must set MinPFPs > 1 (or just don't run this filter)" << std::endl;

}

bool alloutcomes::NPFPFilter::filter(art::Event& e)
{

   art::Handle<std::vector<recob::PFParticle>> Handle_PFParticle;
   std::vector<art::Ptr<recob::PFParticle>> Vect_PFParticle;

   if(!e.getByLabel(f_PFPLabel,Handle_PFParticle)) 
      throw cet::exception("NPFPFilter") << "PFParticle Data Products Not Found!" << std::endl;

   art::fill_ptr_vector(Vect_PFParticle,Handle_PFParticle);

   if(!Vect_PFParticle.size()) return false;

   // Get the neutrino candidate
   size_t neutrinoID = 99999;
   for(const art::Ptr<recob::PFParticle> &pfp : Vect_PFParticle)
      if(pfp->IsPrimary() && (pfp->PdgCode() == 12 || pfp->PdgCode() == 14))
         neutrinoID = pfp->Self();

   // No reco'd neutrino 
   if(neutrinoID == 99999) return false;

   // Count the number of reco'd particles in the hierarchy
   int NDaughters = 0;
   for(const art::Ptr<recob::PFParticle> &pfp : Vect_PFParticle)
      if(pfp->Parent() == neutrinoID) NDaughters++;

  
  return NDaughters >= f_MinPFPs;
}

DEFINE_ART_MODULE(alloutcomes::NPFPFilter)
