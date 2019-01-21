////////////////////////////////////////////////////////////////////////
// Class:       PhotonVis
// Plugin Type: analyzer (art v3_00_00)
// File:        PhotonVis_module.cc
//
// Generated at Sat Jan 19 20:21:57 2019 by David Caratelli using cetskelgen
// from cetlib version v3_04_00.
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

#include "larsim/PhotonPropagation/PhotonVisibilityService.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

class PhotonVis;


class PhotonVis : public art::EDAnalyzer {
public:
  explicit PhotonVis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PhotonVis(PhotonVis const&) = delete;
  PhotonVis(PhotonVis&&) = delete;
  PhotonVis& operator=(PhotonVis const&) = delete;
  PhotonVis& operator=(PhotonVis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  std::vector<double> _pmt_x_pos, _pmt_y_pos, _pmt_z_pos;

};


PhotonVis::PhotonVis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{

  _pmt_x_pos = p.get<std::vector<double>>("X");
  _pmt_y_pos = p.get<std::vector<double>>("Y");
  _pmt_z_pos = p.get<std::vector<double>>("Z");
}

void PhotonVis::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  art::ServiceHandle<phot::PhotonVisibilityService> vis;

  for (size_t pmt=0; pmt < _pmt_x_pos.size(); pmt++) {

    /*
    auto const& x = _pmt_x_pos[pmt] + 10.;
    auto const& y = _pmt_y_pos[pmt];
    auto const& z = _pmt_z_pos[pmt];
    */

    Double_t xyz[3] ={};//  {x,y,z};


    ::art::ServiceHandle<geo::Geometry> geo;
    auto opdet = geo->OpDetFromOpChannel(pmt);
    geo->OpDetGeoFromOpChannel(pmt).GetCenter(xyz); 
    
    xyz[0] += 10;

    double x = xyz[0];
    double y = xyz[1];
    double z = xyz[2];

    std::cout << "PMT ch " << pmt << " with OpDet " << opdet <<
      " @ [ " << x << ", " << y << ", " << z << " ]" << std::endl
	      << "\t with distance to OpChan : " << vis->DistanceToOpDet(xyz,pmt) << std::endl
	      << "\t has visibility " << vis->GetVisibility(xyz, pmt) << std::endl
	      << "\t with distance to OpDet : " << vis->DistanceToOpDet(xyz,opdet) << std::endl
	      << "\t has visibility " << vis->GetVisibility(xyz, opdet) << std::endl << std::endl;

  }// for all PMTs
  return;
}

void PhotonVis::beginJob()
{
  // Implementation of optional member function here.
}

void PhotonVis::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(PhotonVis)
