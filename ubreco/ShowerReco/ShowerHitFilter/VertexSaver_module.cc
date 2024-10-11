////////////////////////////////////////////////////////////////////////
// Class:       PerfectClustering
// Plugin Type: producer (art v2_09_06)
// File:        VertexSaver_module.cc
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
#include "art_root_io/TFileService.h"

// larsoft data-products
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/AnalysisBase/MVAOutput.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "lardata/Utilities/AssociationUtil.h"
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

class VertexSaver;


class VertexSaver : public art::EDProducer {
public:
  explicit VertexSaver(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  VertexSaver(VertexSaver const &) = delete;
  VertexSaver(VertexSaver &&) = delete;
  VertexSaver & operator = (VertexSaver const &) = delete;
  VertexSaver & operator = (VertexSaver &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  std::string fPFPproducer;
  std::string fVTXproducer;

};


VertexSaver::VertexSaver(fhicl::ParameterSet const & p)
: EDProducer(p)
// Initialize member data here.
{
  
  fPFPproducer = p.get<std::string>("PFPproducer");
  fVTXproducer = p.get<std::string>("VTXproducer");

  produces<std::vector<recob::Vertex> >();

}

void VertexSaver::produce(art::Event & e)
{

  std::unique_ptr< std::vector<recob::Vertex> > Vertex_v(new std::vector<recob::Vertex>);

  // load PFParticles
  auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
  // load associated Vertices
  art::FindManyP<recob::Vertex> pfp_vtx_assn_v(pfp_h, e, fVTXproducer);
  
  // loop through PFParticles
  for (size_t p=0; p < pfp_h->size(); p++) {
    std::cout << "this is PFParticle " << p << std::endl;
    auto const& pfp = pfp_h->at(p);
    if ((pfp.PdgCode() == 12) || (pfp.PdgCode() == 14)){
      // get associated vertex
      auto vtx_v = pfp_vtx_assn_v.at(p);
      for (size_t v=0; v < vtx_v.size(); v++) {
	//auto const& vtx = vtx_v.at(v);
	auto const& vtx = vtx_v.at(v);
	Double_t xyz[3] = {};
	vtx->XYZ(xyz);
	recob::Vertex newvtx(xyz);
	Vertex_v->emplace_back(newvtx);
	// make a copy and replace
      }// for all associated vertices
    }// if PDG == 14 or 12
  }// for all PFParticles
  
  e.put(std::move(Vertex_v));
  
}

void VertexSaver::beginJob()
{

}

void VertexSaver::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(VertexSaver)
