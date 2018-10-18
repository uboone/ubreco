////////////////////////////////////////////////////////////////////////
// Class:       Pi0ShrFilter
// Plugin Type: filter (art v2_11_03)
// File:        Pi0ShrFilter_module.cc
//
// Generated at Wed Oct 17 18:51:47 2018 by David Caratelli using cetskelgen
// from cetlib version v3_03_01.
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

#include "lardataobj/RecoBase/Vertex.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardata/Utilities/AssociationUtil.h"

#include <memory>

class Pi0ShrFilter;


class Pi0ShrFilter : public art::EDFilter {
public:
  explicit Pi0ShrFilter(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Pi0ShrFilter(Pi0ShrFilter const &) = delete;
  Pi0ShrFilter(Pi0ShrFilter &&) = delete;
  Pi0ShrFilter & operator = (Pi0ShrFilter const &) = delete;
  Pi0ShrFilter & operator = (Pi0ShrFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

};


Pi0ShrFilter::Pi0ShrFilter(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
}

bool Pi0ShrFilter::filter(art::Event & e)
{

  auto const& mct_h = e.getValidHandle<std::vector<simb::MCTruth> >("generator");

  auto mct = mct_h->at(0);
  size_t npart = mct.NParticles();

  // save the trackID of the pi0
  int npi0 = 0;

  for (size_t i=0; i < npart; i++){
    auto const& part = mct.GetParticle(i);
    if ( (part.PdgCode() == 111) and (part.StatusCode() == 1) ){
      npi0 += 1;
      break;
    }
  }// for all particles

  if (npi0 == 0)
    return false;

  return true;
  
}

void Pi0ShrFilter::beginJob()
{
  // Implementation of optional member function here.
}

void Pi0ShrFilter::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(Pi0ShrFilter)
