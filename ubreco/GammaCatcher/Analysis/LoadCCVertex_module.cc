////////////////////////////////////////////////////////////////////////
// Class:       LoadCCVertex
// Plugin Type: analyzer (art v2_11_03)
// File:        LoadCCVertex_module.cc
//
// Generated at Tue Jan  8 12:25:39 2019 by David Caratelli using cetskelgen
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

#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardata/Utilities/AssociationUtil.h"

class LoadCCVertex;


class LoadCCVertex : public art::EDAnalyzer {
public:
  explicit LoadCCVertex(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  LoadCCVertex(LoadCCVertex const &) = delete;
  LoadCCVertex(LoadCCVertex &&) = delete;
  LoadCCVertex & operator = (LoadCCVertex const &) = delete;
  LoadCCVertex & operator = (LoadCCVertex &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // producers
  std::string fAssnProducer;

};


LoadCCVertex::LoadCCVertex(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  fAssnProducer = p.get<std::string>("AssnProducer");
}

void LoadCCVertex::analyze(art::Event const & e)
{

  Double_t xyz[3] = {};

  art::Handle< art::Assns<recob::Vertex,recob::Track,void> > numuCCassn_h;
  e.getByLabel(fAssnProducer,numuCCassn_h);
  if (numuCCassn_h->size() != 1){
    std::cout << "Number of vertices != 1 -> ERROR ERROR ERROR" << std::endl;
    return;
  }
  
  numuCCassn_h->at(0).first->XYZ(xyz);

  std::cout << "Vertex coordinates are : " << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << std::endl;

  auto ccmuon = numuCCassn_h->at(0).second;

  for (unsigned int i = 0; i < ccmuon->NumberTrajectoryPoints(); i++) {
    // project a point into 2D:
      if (ccmuon->HasValidPoint(i)) {
	auto loc = ccmuon->LocationAtPoint(i);
	std::cout << "muon track point : [ " << loc.X() << ", " << loc.Y() << ", " << loc.Z() << " ]" << std::endl;
	  }// if valid point
  }// for all points
  
  return;
}

void LoadCCVertex::beginJob()
{
  // Implementation of optional member function here.
}

void LoadCCVertex::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(LoadCCVertex)
