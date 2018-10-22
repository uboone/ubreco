////////////////////////////////////////////////////////////////////////
// Class:       StopMu
// Plugin Type: analyzer (art v2_11_03)
// File:        StopMu_module.cc
//
// Generated at Sun Oct 21 21:38:30 2018 by David Caratelli using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

class StopMu;


class StopMu : public art::EDAnalyzer {
public:
  explicit StopMu(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  StopMu(StopMu const &) = delete;
  StopMu(StopMu &&) = delete;
  StopMu & operator = (StopMu const &) = delete;
  StopMu & operator = (StopMu &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

};


StopMu::StopMu(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

void StopMu::analyze(art::Event const & e)
{
  // Implementation of required member function here.
}

void StopMu::beginJob()
{
  // Implementation of optional member function here.
}

void StopMu::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(StopMu)
