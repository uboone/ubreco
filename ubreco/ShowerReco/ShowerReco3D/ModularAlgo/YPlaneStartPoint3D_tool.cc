#ifndef YPLANESTARTPOINT3D_CXX
#define YPLANESTARTPOINT3D_CXX

#include <iostream>
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
    void do_reconstruction(const ::protoshower::ProtoShower &, Shower_t &);

  private:

    double _wire2cm, _time2cm;
    
  };
  
  YPlaneStartPoint3D::YPlaneStartPoint3D(const fhicl::ParameterSet& pset)
  {
    _name = "YPlaneStartPoint3D"; 
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    _wire2cm = geom->WirePitch(0,1,0);
    _time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );

  }

  void YPlaneStartPoint3D::do_reconstruction( const ::protoshower::ProtoShower & proto_shower,
					      Shower_t& resultShower)
{


    //if the module does not have 2D cluster info -> fail the reconstruction
    if (!proto_shower.hasCluster2D()) {
        std::stringstream ss;
        ss << "Fail @ algo " << this->name() << " due to missing 2D cluster";
        throw ShowerRecoException(ss.str());
    }

    if (proto_shower.hasVertex() == false){
      std::cout << "Number of vertices is not one!" << std::endl;
      return;
    }
    
    // get the proto-shower 3D vertex
    auto const& vtx3D = proto_shower.vertex();
    auto const& dir3D = resultShower.fDCosStart;
    double d2D;

    std::cout << std::endl << " START 3D Start Point CALCULATION" << std::endl;

    std::cout << " Vtx 3D : [" << vtx3D[0] << ", " << vtx3D[1] << ", " << vtx3D[2] << " ]" << std::endl;
    std::cout << " dir 3D : [" << dir3D[0] << ", " << dir3D[1] << ", " << dir3D[2] << " ]" << std::endl;
    
    auto & clusters = proto_shower.clusters();

    // STEP 1:
    // identify the 2D start point on the collection plane.
    // given the shower's reconstructed 3D direction vector
    // we will find where this point lies in 3D.

    ::util::PxHit strtpt0;
    int    pl0 = 0;
    
    for (size_t i=0; i < clusters.size(); i++) {

      auto const& clus = clusters.at(i);

      // plane
      auto const& pl = clus._plane;

      if (pl != 2) continue;
      
      // project vertex onto this plane
      //auto const& vtx2D = util::PxPoint(pl,0,0);//geomH->Get2DPointProjection(vtx,pl);
      auto const* geom = ::lar::providerFrom<geo::Geometry>();
      auto wire = geom->WireCoordinate(vtx3D[1],vtx3D[2],geo::PlaneID(0,0,pl)) * _wire2cm;
      auto time = vtx3D[0];
      auto const& start = clus._start;

      util::PxPoint vtx2D(pl,wire,time);
      
      // calculate 2D distance between vertex and cluster start point on plane
      d2D = sqrt( (start.w - vtx2D.w) * (start.w - vtx2D.w) + (start.t - vtx2D.t) * (start.t - vtx2D.t) );

      std::cout << " Strt 2D : [" << start.w << ", " << start.t << " ]" << std::endl;
      std::cout << " Vtx  2D : [" << vtx2D.w << ", " << vtx2D.t << " ]" << std::endl;

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

    std::cout << " d2D : " << d2D << std::endl;
    std::cout << " d3D : " << d3D << std::endl;

    // extend by this amount in 3D
    auto start3D = vtx3D + d3D * dir3D;

    std::cout << " Strt 3D : [" << start3D[0] << ", " << start3D[1] << ", " << start3D[2] << " ]" << std::endl;

    std::cout << " END 3D Start Point CALCULATION" << std::endl << std::endl;
    
    resultShower.fXYZStart = start3D;

    //std::cout << "DONE " << std::endl << std::endl;
    
}

  DEFINE_ART_CLASS_TOOL(YPlaneStartPoint3D)
} //showerreco

#endif
