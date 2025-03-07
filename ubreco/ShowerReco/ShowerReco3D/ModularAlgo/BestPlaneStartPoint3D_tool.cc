#ifndef BESTPLANESTARTPOINT3D_CXX
#define BESTPLANESTARTPOINT3D_CXX

#include <iostream>
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcore/Geometry/WireReadout.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "ubreco/ShowerReco/ShowerReco3D/Base/ShowerRecoModuleBase.h"
/**
   \class ShowerRecoModuleBase
   User defined class ShowerRecoModuleBase ... these comments are used to generate
   doxygen documentation!
 */

namespace showerreco {

  class BestPlaneStartPoint3D : public ShowerRecoModuleBase {
    
  public:
    
    /// Default constructor
    BestPlaneStartPoint3D(const fhicl::ParameterSet& pset);

    /// Default destructor
    ~BestPlaneStartPoint3D() {}
    
    /// Inherited/overloaded function from ShowerRecoModuleBase
    void do_reconstruction(util::GeometryUtilities const&,
                           const ::protoshower::ProtoShower &, Shower_t &);

  private:

    double _wire2cm, _time2cm;
    
  };
  
  BestPlaneStartPoint3D::BestPlaneStartPoint3D(const fhicl::ParameterSet& pset)
  {
    _name = "BestPlaneStartPoint3D"; 
    auto const& channelMap = art::ServiceHandle<geo::WireReadout>()->Get();
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();
    auto const detp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);
    _wire2cm = channelMap.Plane(geo::PlaneID{0,0,0}).WirePitch();
    _time2cm = sampling_rate(clockData) / 1000.0 * detp.DriftVelocity( detp.Efield(), detp.Temperature());
  }

  void BestPlaneStartPoint3D::do_reconstruction(util::GeometryUtilities const&,
                                             const ::protoshower::ProtoShower & proto_shower,
					      Shower_t& resultShower)
{

    //if the module does not have 2D cluster info -> fail the reconstruction
    if (!proto_shower.hasCluster2D()) {
        std::stringstream ss;
        ss << "Fail @ algo " << this->name() << " due to missing 2D cluster";
        throw ShowerRecoException(ss.str());
    }

    if (proto_shower.hasVertex() == false){
      std::stringstream ss;
      ss << "Fail @ algo " << this->name() << " due to != 1 vertex";
      throw ShowerRecoException(ss.str());
      return;
    }
    
    // get the proto-shower 3D vertex
    auto const& vtx3D = proto_shower.vertex();
    auto const vtx3D_pt = geo::vect::toPoint(vtx3D);
    auto const& dir3D = resultShower.fDCosStart;
    float d2D = -1;
    int bestPlane = -1;

    auto & clusters = proto_shower.clusters();

    auto const& channelMap = art::ServiceHandle<geo::WireReadout>()->Get();
    for (size_t i=0; i < clusters.size(); i++) {

      auto const& clus = clusters.at(i);

      // plane
      int pl = clus._plane;

      // project vertex onto this plane
      auto wire = channelMap.Plane(geo::PlaneID(0,0,pl)).WireCoordinate(vtx3D_pt) * _wire2cm;
      auto time = vtx3D.X();
      auto const& start = clus._start;

      util::PxPoint vtx2D(pl,wire,time);
      
      // calculate 2D distance between vertex and cluster start point on plane
      float tmpd2D = sqrt( (start.w - vtx2D.w) * (start.w - vtx2D.w) + (start.t - vtx2D.t) * (start.t - vtx2D.t) );
      std::cout << "plane=" << pl << " d2D=" << tmpd2D << " cldist=" << abs(clus._start.w - clus._end.w) << std::endl;

      if (tmpd2D<d2D || d2D<0.) {//do we need a more sophisticated way to choose the best plane?
	d2D = tmpd2D;
	bestPlane = pl;
      }

    }// for all clusters

    if (d2D<0) {
      std::stringstream ss;
      ss << "Fail @ algo " << this->name() << " due to negative 2D distance";
      throw ShowerRecoException(ss.str());
    }

    std::cout << "BestPlaneStartPoint3D -- best plane=" << bestPlane << std::endl;
    TVector3 wdir(channelMap.Plane(geo::PlaneID(0,0,bestPlane)).MiddleWire().Direction().X(),
		  channelMap.Plane(geo::PlaneID(0,0,bestPlane)).MiddleWire().Direction().Y(),
		  channelMap.Plane(geo::PlaneID(0,0,bestPlane)).MiddleWire().Direction().Z());
    auto dotp = dir3D.Dot(wdir);
    //std::cout << channelMap.Plane(geo::PlaneID(0,0,bestPlane)).MiddleWire().Direction() << " " << dotp << " " << (1 - dir3D[1]*dir3D[1] ) << " " << (1 - dotp*dotp)<< std::endl;

    double f   = (1 - dotp*dotp);
    double d3D = d2D / f;

    // extend by this amount in 3D
    auto start3D = vtx3D + d3D * dir3D;

    resultShower.fXYZStart = start3D;

}

  DEFINE_ART_CLASS_TOOL(BestPlaneStartPoint3D)
} //showerreco

#endif
