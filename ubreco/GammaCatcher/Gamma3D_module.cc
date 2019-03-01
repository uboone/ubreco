////////////////////////////////////////////////////////////////////////
// Class:       Gamma3D
// Plugin Type: producer (art v2_11_03)
// File:        Gamma3D_module.cc
//
// Generated at Wed Feb 27 14:08:46 2019 by David Caratelli using cetskelgen
// from cetlib version v3_03_01.
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

#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Track.h"

#include "TTree.h"

#include <memory>

class Gamma3D;


class Gamma3D : public art::EDProducer {
public:
  explicit Gamma3D(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Gamma3D(Gamma3D const &) = delete;
  Gamma3D(Gamma3D &&) = delete;
  Gamma3D & operator = (Gamma3D const &) = delete;
  Gamma3D & operator = (Gamma3D &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // configurable inputs to reconstruction
  float f2DcutY;
  float f2DcutUV;
  float deltaY;

  TTree* _tree;
  int   _yhits;
  float _ycharge;
  float _ywire;
  float _ytime;
  int   _umacthes;
  int   _vmatches;
  int   _deltaYmin;

};


Gamma3D::Gamma3D(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  produces< std::vector< recob::SpacePoint > >();
  produces< art::Assns <recob::Cluster, recob::SpacePoint> >();


  f2DcutY  = p.get<float>("2DcutY");
  f2DcutUV = p.get<float>("2DcutUV");
  fdeltaY  = p.get<float>("deltaY");
  
}

void Gamma3D::produce(art::Event & e)
{
  
  std::unique_ptr< std::vector< recob::SpacePoint> > SpacePoint_v(new std::vector<recob::SpacePoint>);
  std::unique_ptr< art::Assn <recob::Cluster, recob::SpacePoint> >(new art::Assn<recob::Cluster,recob::SpacePoint> );

  // code goes here...
  
  e.put(std::move(SpacePoint_v));
  e.put(std::move(sps_clus_assn_v));

}

void Gamma3D::beginJob()
{
  // Implementation of optional member function here.
}

void Gamma3D::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(Gamma3D)
