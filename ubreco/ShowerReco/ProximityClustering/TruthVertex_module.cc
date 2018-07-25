////////////////////////////////////////////////////////////////////////
// Class:       TruthVertex
// Plugin Type: filter (art v2_09_06)
// File:        TruthVertex_module.cc
//
// Generated at Mon Mar  5 14:14:25 2018 by David Caratelli using cetskelgen
// from cetlib version v3_01_03.
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

#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "art/Framework/Services/Optional/TFileService.h"

#include "lardataobj/RecoBase/Vertex.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

class TruthVertex;


class TruthVertex : public art::EDFilter {
public:
  explicit TruthVertex(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TruthVertex(TruthVertex const &) = delete;
  TruthVertex(TruthVertex &&) = delete;
  TruthVertex & operator = (TruthVertex const &) = delete;
  TruthVertex & operator = (TruthVertex &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  // producers
  std::string fMCTruthProducer;

  double _time2cm, _wire2cm;

};


TruthVertex::TruthVertex(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  produces< std::vector< recob::Vertex > >();
  fMCTruthProducer = p.get<std::string>("MCTruthProducer");

  auto const* geom = ::lar::providerFrom<geo::Geometry>();
  auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  _wire2cm = geom->WirePitch(0,1,0);
  _time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );
}

bool TruthVertex::filter(art::Event & e)
{

  // produce vertex
  std::unique_ptr< std::vector<recob::Vertex> > Vtx_v(new std::vector<recob::Vertex>);

  // space charge
  auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

  Double_t xyz[3] = {};
  double tvtx;

  auto const& mct_h = e.getValidHandle<std::vector<simb::MCTruth> >(fMCTruthProducer);
  auto mct = mct_h->at(0);
  size_t npart = mct.NParticles();

  bool foundPi0 = false;

  for (size_t i=0; i < npart; i++){
    auto const& part = mct.GetParticle(i);
    if ( (part.PdgCode() == 111) and (part.StatusCode() == 1) ){
      xyz[0] = part.Trajectory().X(0);
      xyz[1] = part.Trajectory().Y(0);
      xyz[2] = part.Trajectory().Z(0);
      tvtx   = part.Trajectory().T(0);
      foundPi0 = true;
      break;
    }
  }

  if (foundPi0 == false) {
    e.put(std::move(Vtx_v));
    return false;
  }

  std::cout << "truth vertex @ [" << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << " ]" << std::endl;
  
  // get position offsets:
  auto offset = SCE->GetPosOffsets(geo::Point_t(xyz[0],xyz[1],xyz[2]));
  xyz[0] += (-offset.X() + 0.7);
  xyz[1] += offset.Y();
  xyz[2] += offset.Z();
  std::cout << "Offset. X : " << offset.X() << ", Y : " << offset.Y() << ", Z : " << offset.Z() << std::endl;
  
  auto vtxtick   = (tvtx/1000.) * 2.;
  auto vtxtimecm = vtxtick * _time2cm;

  xyz[0] += vtxtimecm;

  std::cout << "truth vertex @ [" << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << " ]" << std::endl;

  if ( (xyz[0] < 1) || (xyz[0] > 249) || (xyz[1] < -115) || (xyz[1] > 115) || (xyz[2] < 1) || (xyz[2] > 1035) ) {
    e.put(std::move(Vtx_v));
    return false;
  }
  
  recob::Vertex vtx(xyz);
  Vtx_v->emplace_back(vtx);
  
  e.put(std::move(Vtx_v));

  return true;
}

void TruthVertex::beginJob()
{
  // Implementation of optional member function here.
}

void TruthVertex::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(TruthVertex)
