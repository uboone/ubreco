#ifndef YPLANESTARTPOINT3D_CXX
#define YPLANESTARTPOINT3D_CXX

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

  class YPlaneStartPoint3D : public ShowerRecoModuleBase {
    
  public:
    
    /// Default constructor
    YPlaneStartPoint3D(const fhicl::ParameterSet& pset);

    /// Default destructor
    ~YPlaneStartPoint3D() {}
    
    /// Inherited/overloaded function from ShowerRecoModuleBase
    void do_reconstruction(util::GeometryUtilities const&,
                           const ::protoshower::ProtoShower &, Shower_t &);

  private:

    double _wire2cm, _time2cm;
    
  };
  
  YPlaneStartPoint3D::YPlaneStartPoint3D(const fhicl::ParameterSet& pset)
  {
    _name = "YPlaneStartPoint3D"; 
    auto const& channelMap = art::ServiceHandle<geo::WireReadout>()->Get();
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();
    auto const detp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);
    _wire2cm = channelMap.Plane(geo::PlaneID{0,0,0}).WirePitch();
    _time2cm = sampling_rate(clockData) / 1000.0 * detp.DriftVelocity( detp.Efield(), detp.Temperature());
  }

  void YPlaneStartPoint3D::do_reconstruction(util::GeometryUtilities const&,
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
    double d2D;

    auto & clusters = proto_shower.clusters();

    // STEP 1:
    // identify the 2D start point on the collection plane.
    // given the shower's reconstructed 3D direction vector
    // we will find where this point lies in 3D.

    ::util::PxHit strtpt0;
    int    pl0 = 0;
    
    auto const& channelMap = art::ServiceHandle<geo::WireReadout>()->Get();
    for (size_t i=0; i < clusters.size(); i++) {

      auto const& clus = clusters.at(i);

      // plane
      auto const& pl = clus._plane;

      if (pl != 2) continue;
      
      // project vertex onto this plane
      //auto const& vtx2D = util::PxPoint(pl,0,0);//geomH->Get2DPointProjection(vtx,pl);
      auto wire = channelMap.Plane(geo::PlaneID(0,0,pl)).WireCoordinate(vtx3D_pt) * _wire2cm;
      auto time = vtx3D.X();
      auto const& start = clus._start;

      util::PxPoint vtx2D(pl,wire,time);
      
      // calculate 2D distance between vertex and cluster start point on plane
      d2D = sqrt( (start.w - vtx2D.w) * (start.w - vtx2D.w) + (start.t - vtx2D.t) * (start.t - vtx2D.t) );


      pl0 = pl;
      strtpt0 = start;
      
    }// for all clusters

    if (pl0 != 2) {
      std::stringstream ss;
      ss << "Fail @ algo " << this->name() << " due to missing Collection-Plane 2D cluster";
      throw ShowerRecoException(ss.str());
    }

    double f   = (1 - dir3D[1]*dir3D[1] );
    double d3D = d2D / f;//(1-fabs(dir3D[1]));

    // extend by this amount in 3D
    auto start3D = vtx3D + d3D * dir3D;

    resultShower.fXYZStart = start3D;

}

  DEFINE_ART_CLASS_TOOL(YPlaneStartPoint3D)
} //showerreco

#endif
