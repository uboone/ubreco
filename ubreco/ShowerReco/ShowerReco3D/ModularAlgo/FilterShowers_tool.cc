#ifndef FILTERSHOWERS_CXX
#define FILTERSHOWERS_CXX

#include <iostream>

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcore/Geometry/WireReadout.h"
#include "ubreco/ShowerReco/ShowerReco3D/Base/ShowerRecoModuleBase.h"

#include <sstream>

/**
   \class ShowerRecoModuleBase
   User defined class ShowerRecoModuleBase ... these comments are used to generate
   doxygen documentation!
 */

namespace showerreco {

  class FilterShowers : public ShowerRecoModuleBase {
    
  public:
    
    /// Default constructor
    FilterShowers(const fhicl::ParameterSet& pset);
    
    /// Default destructor
    ~FilterShowers() {}

    void configure(const fhicl::ParameterSet& pset);
    
    
    void do_reconstruction(util::GeometryUtilities const&,
                           const ::protoshower::ProtoShower &, Shower_t &);
    
  private:
    
    double _anglecut;

    double _wire2cm, _time2cm;
    
  };

  FilterShowers::FilterShowers(const fhicl::ParameterSet& pset)
  {
    _name = "FilterShowers";
    auto const& channelMap = art::ServiceHandle<geo::WireReadout>()->Get();
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();
    auto const detp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);
    _wire2cm = channelMap.Plane(geo::PlaneID{0,0,0}).WirePitch();
    _time2cm = sampling_rate(clockData) / 1000.0 * detp.DriftVelocity( detp.Efield(), detp.Temperature() );
    configure(pset);
  }

  void FilterShowers::configure(const fhicl::ParameterSet& pset)
  {
    _anglecut = pset.get<double>("anglecut");
  }

void FilterShowers::do_reconstruction(util::GeometryUtilities const&,
                                      const ::protoshower::ProtoShower & proto_shower,
                                    Shower_t& resultShower) {

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
  
  // take the vtx -> start point direction as the 3D direction
  auto const vtx = geo::vect::toPoint(proto_shower.vertex());
  
  //auto const& geomH = ::util::GeometryUtilities::GetME();
  

  
  auto & clusters = proto_shower.clusters();
  
  // get 3D shower direction projected on collection-plane
  double slope3D = resultShower.fDCosStart[0] / resultShower.fDCosStart[2];
  slope3D       /= sqrt( ( resultShower.fDCosStart[0] * resultShower.fDCosStart[0] ) +
			 ( resultShower.fDCosStart[2] * resultShower.fDCosStart[2] ) );
  

  // angle between projected 3D direction on plane (slope3D)
  // and average hit direction
  double clusterhitangle = 0.;
  
  // is there a collection-plane shower?
  bool collectioncluster = false;
  
  auto const& channelMap = art::ServiceHandle<geo::WireReadout>()->Get();
  for (auto const& clus : clusters) {

    if (clus._plane != 2) continue;


    // project vertex onto this plane
    auto wire = channelMap.Plane(geo::PlaneID(0,0,clus._plane)).WireCoordinate(vtx) * _wire2cm;
    auto time = vtx.X();
    auto const& vtx2D = util::PxPoint(2,wire,time);// geomH->Get2DPointProjection(vtx,2);
    
    // get hits and start-point, calculate average angle of hits w.r.t. projected
    // shower direction

    collectioncluster = true;

    // cluster hit closest to vertex is the point from which to compute hit angle
    double sw = 0;
    double st = 0;
    double vtxdmin = 10000000.;

    for (auto const& hit : clus._hits) {
      double dd = ( (hit.w-vtx2D.w)*(hit.w-vtx2D.w) + (hit.t-vtx2D.t)*(hit.t-vtx2D.t) );
      if (dd < vtxdmin) { vtxdmin = dd; sw = hit.w; st = hit.t; }
    }    
    
    // get cluster start point
    //auto start = clus.start_point;

    if (_verbose)
      std::cout << "start pt @ " << sw << ", " << st << std::endl;

    for (auto const& hit : clus._hits) {

      if ( (hit.t == st) || (hit.w == sw) ) {
	if (_verbose)
	  std::cout << "found overlapping hit!" << std::endl;
	continue;
      }

      double hitslope = (hit.t - st) / (hit.w - sw );
      double hitangle = fabs( atan( ( hitslope - slope3D ) / ( 1 + slope3D * hitslope ) ) ); 
      clusterhitangle += hitangle;
      //std::cout << "\t hit @ [ " << hit.w << ", " << hit.t << "] -> [" << (hit.w-sw) << ", " << (hit.t-st) << "] has angle : " << hitangle * 180. / 3.14 << std::endl;
    }

    clusterhitangle /= ( clus._hits.size() - 1);
    clusterhitangle = fabs(clusterhitangle);

  }// for all clusters

  //if the module does not have 2D cluster info -> fail the reconstruction
  if ( collectioncluster == false ) {
    std::stringstream ss;
    ss << "Fail @ algo " << this->name() << " due to missing 2D collection-plane cluster";
    throw ShowerRecoException(ss.str());
  }

  clusterhitangle *= (180./3.14);

  if (_verbose) {
    std::cout << "Hit Angle = " << clusterhitangle << std::endl;
    std::cout << "angle cut : " << _anglecut << std::endl;
    //std::cout << "Opening Angle = " << resultShower.fOpeningAngle * 180 / 3.14 << std::endl;
  }
  
  if (clusterhitangle > _anglecut) {
  //if ( (clusterhitangle > (resultShower.fOpeningAngle * 180 / 3.14)) and (clusterhitangle > _anglecut) ) {
    std::stringstream ss;
    //ss << "Fail @ algo " << this->name() << " Hit angle w.r.t. 3D dir on Y plane is too large : " << clusterhitangle;
    ss << "aaa";
    //throw std::exception();
    throw ShowerRecoException(ss.str());
  }
}

  DEFINE_ART_CLASS_TOOL(FilterShowers)
} //showerreco

#endif
