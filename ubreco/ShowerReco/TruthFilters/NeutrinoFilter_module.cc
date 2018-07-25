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

#include "TTree.h"

class NeutrinoFilter;


class NeutrinoFilter : public art::EDFilter {
public:
  explicit NeutrinoFilter(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NeutrinoFilter(NeutrinoFilter const &) = delete;
  NeutrinoFilter(NeutrinoFilter &&) = delete;
  NeutrinoFilter & operator = (NeutrinoFilter const &) = delete;
  NeutrinoFilter & operator = (NeutrinoFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  double _time2cm, _wire2cm;

  double fECut;

  TTree* _tree;
  double _nu_energy, _e_energy, _p_energy;

};


NeutrinoFilter::NeutrinoFilter(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  produces< std::vector< recob::Vertex > >();

  auto const* geom = ::lar::providerFrom<geo::Geometry>();
  auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  _wire2cm = geom->WirePitch(0,1,0);
  _time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );

  fECut = p.get<double>("ECut");

  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("_tree","Nue Truth TTree");
  _tree->Branch("_nu_energy",&_nu_energy,"nu_energy/D");
  _tree->Branch("_e_energy" ,&_e_energy ,"e_energy/D" );
  _tree->Branch("_p_energy" ,&_p_energy ,"p_energy/D" );

}

bool NeutrinoFilter::filter(art::Event & e)
{

  auto const& mct_h = e.getValidHandle<std::vector<simb::MCTruth> >("generator");

  // produce vertex
  std::unique_ptr< std::vector<recob::Vertex> > Vtx_v(new std::vector<recob::Vertex>);
  // space charge
  auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  
  auto mct = mct_h->at(0);
  auto neutrino = mct.GetNeutrino();

  auto nu     = neutrino.Nu();
  auto lepton = neutrino.Lepton();
  auto ccnc   = neutrino.CCNC();

  _nu_energy = nu.Trajectory().E(0);

  std::cout << "Nu PdgCode : " << nu.PdgCode()
	    << "\t Lepton PdgCode : " << lepton.PdgCode() 
	    << "\t CC/NC : " << ccnc << std::endl;

  if ( ccnc == 1 ) {
    e.put(std::move(Vtx_v));
    return false;
  }

  _e_energy = 0.;
  _p_energy = 0.;
  for (size_t i=0; i < (size_t)(mct.NParticles()); i++){
    auto const& part = mct.GetParticle(i);
    if ( (fabs(part.PdgCode()) == 11) and (part.StatusCode() == 1) ){
      if (part.Trajectory().E(0) > _e_energy)
	_e_energy = part.Trajectory().E(0);
    }// if electron
    if ( (fabs(part.PdgCode()) == 2212) and (part.StatusCode() == 1) ){
      if (part.Trajectory().E(0) > _e_energy)
	_p_energy = part.Trajectory().E(0);
    }// if proton
  }

  _tree->Fill();

  if (_nu_energy < fECut) {
    e.put(std::move(Vtx_v));
    return false;
  }

  Double_t xyz[3] = {};
  double tvtx;

  // vertex coordinates from neutrino end point
  xyz[0] = lepton.Trajectory().X(0);
  xyz[1] = lepton.Trajectory().Y(0);
  xyz[2] = lepton.Trajectory().Z(0);
  tvtx   = lepton.Trajectory().T(0);

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


  Vtx_v->emplace_back(recob::Vertex(xyz));

  e.put(std::move(Vtx_v));

  return true;
}

void NeutrinoFilter::beginJob()
{
  // Implementation of optional member function here.
}

void NeutrinoFilter::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(NeutrinoFilter)
