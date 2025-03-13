////////////////////////////////////////////////////////////////////////
// Class:       EmptyFilter
// Plugin Type: filter (Unknown Unknown)
// File:        EmptyFilter_module.cc
//
// Generated at Mon Mar 10 12:25:31 2025 by Erin Yandel using cetskelgen
// from cetlib version 3.18.02.
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

#include <memory>

class EmptyFilter;


class EmptyFilter : public art::EDFilter {
public:
  explicit EmptyFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  EmptyFilter(EmptyFilter const&) = delete;
  EmptyFilter(EmptyFilter&&) = delete;
  EmptyFilter& operator=(EmptyFilter const&) = delete;
  EmptyFilter& operator=(EmptyFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

private:

  // Declare member data here.

};


EmptyFilter::EmptyFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

bool EmptyFilter::filter(art::Event& e)
{
  // Implementation of required member function here.
  return true;
}

DEFINE_ART_MODULE(EmptyFilter)
