////////////////////////////////////////////////////////////////////////
// Class:       DummyRawDigitCreator
// Plugin Type: producer (art v3_01_02)
// File:        DummyRawDigitCreator_module.cc
//
// Generated at Thu Jan  2 11:00:03 2020 by Wesley Ketchum using cetskelgen
// from cetlib version v3_05_01.
//
// Doing this because we need to make a dummy association.
//
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

#include <memory>

#include "lardataobj/RawData/RawDigit.h"

namespace util {
  class DummyRawDigitCreator;
}


class util::DummyRawDigitCreator : public art::EDProducer {
public:
  explicit DummyRawDigitCreator(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DummyRawDigitCreator(DummyRawDigitCreator const&) = delete;
  DummyRawDigitCreator(DummyRawDigitCreator&&) = delete;
  DummyRawDigitCreator& operator=(DummyRawDigitCreator const&) = delete;
  DummyRawDigitCreator& operator=(DummyRawDigitCreator&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.
  size_t fNChannels;
};


util::DummyRawDigitCreator::DummyRawDigitCreator(fhicl::ParameterSet const& p)
  : 
  EDProducer{p},
  fNChannels(p.get<size_t>("NChannels",8256))
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
    produces< std::vector< raw::RawDigit > >();
  
}

void util::DummyRawDigitCreator::produce(art::Event& e)
{
  std::unique_ptr< std::vector<raw::RawDigit> > rawDigitPtrVec(new std::vector<raw::RawDigit>());
  
  raw::RawDigit::ADCvector_t empty_wvfm;

  for(size_t i_ch=0; i_ch<fNChannels; ++i_ch)
    rawDigitPtrVec->emplace_back(i_ch,0,empty_wvfm);

  e.put(std::move(rawDigitPtrVec));
}

DEFINE_ART_MODULE(util::DummyRawDigitCreator)
