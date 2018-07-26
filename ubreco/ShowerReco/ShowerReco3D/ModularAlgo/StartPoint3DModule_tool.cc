#ifndef STARTPOINT3DMODULE_CXX
#define STARTPOINT3DMODULE_CXX

#include <iostream>
#include "ubreco/ShowerReco/ShowerReco3D/Base/ShowerRecoModuleBase.h"
/**
   \class ShowerRecoModuleBase
   User defined class ShowerRecoModuleBase ... these comments are used to generate
   doxygen documentation!
 */

namespace showerreco {

  class StartPoint3DModule : public ShowerRecoModuleBase {
    
  public:
    
    /// Default constructor
    StartPoint3DModule(const fhicl::ParameterSet& pset);

    
    /// Default destructor
    ~StartPoint3DModule() {}
    
    /// Inherited/overloaded function from ShowerRecoModuleBase
    void do_reconstruction(const ::protoshower::ProtoShower &, Shower_t &);
    
  };

  StartPoint3DModule::StartPoint3DModule(const fhicl::ParameterSet& pset)
    {
      _name = "StartPoint3DModule"; 
    }
  
  void StartPoint3DModule::do_reconstruction( const ::protoshower::ProtoShower & proto_shower,
					      Shower_t& resultShower)
{

    //if the module does not have 2D cluster info -> fail the reconstruction
    if (!proto_shower.hasCluster2D()) {
        std::stringstream ss;
        ss << "Fail @ algo " << this->name() << " due to missing 2D cluster";
        throw ShowerRecoException(ss.str());
    }

    auto & clusters = proto_shower.clusters();

    /**
    *  Use 2 best planes to calculate Y, Z of start point.
    *  Use their time average to calculate X.
    *  First calculate the worst plane and don't use that guy.
    */
    float minClusDist      =  999999 ;
    int worstPlane  = -1      ;
    int planeTemp   =  0      ;

    if ( clusters.size() > 2 ) {
        for ( auto const c : clusters ) {
            float distTemp = abs ( c._start.w - c._end.w );

            if ( distTemp < minClusDist ) {
                minClusDist = distTemp ;
                worstPlane = planeTemp;
            }

            planeTemp ++ ;
        }
    }

    /**
    *  Store planes and wire start, calculate average time for best 2 clusters.
    */
    double sX = 0 ;
    double sY ;
    double sZ ;

    std::vector<unsigned int> wireStarts(0) ;
    std::vector<int> planes(0) ;
    
    //auto const& geomH = ::util::GeometryUtilities::GetME();
    auto const* geom  = ::lar::providerFrom<geo::Geometry>();

    for ( auto const c : clusters ) {

        if ((int) c._plane != worstPlane) {
            wireStarts.emplace_back(
                int(c._start.w / geom -> WirePitch(0) ) ) ;
            planes.emplace_back( c._plane ) ;
            sX += c._start.t;
        }
    }

    /**
    *  Caluclate intersection point in Y,Z of the 2 best planes.
    *  Calculate their average start time
    */
    if (planes.size() < 2){
      std::stringstream ss;
      ss << "Fail @ algo " << this->name() << " due to too few (1) 2D start points";
      throw ShowerRecoException(ss.str());
    }

    // check start/end point range
    // wireStarts is a vector of unsigned ints, so the first check is always true
    //if ( (wireStarts.at(0) < 0) or (wireStarts.at(0) > geom->Nwires(planes.at(0))) or
 //	 (wireStarts.at(1) < 0) or (wireStarts.at(1) > geom->Nwires(planes.at(1))) ) {
    if (  (wireStarts.at(0) > geom->Nwires(planes.at(0))) or
	  (wireStarts.at(1) > geom->Nwires(planes.at(1))) ) {
      std::stringstream ss;
      ss << "Fail @ algo " << this->name() << " due to wires out of range";
      throw ShowerRecoException(ss.str());
    } 

    geom->IntersectionPoint(wireStarts.at(0) , wireStarts.at(1) , planes.at(0), planes.at(1), 0, 0, sY, sZ );

    // check if reconstructed start point is outside of TPC volume    
    /*
    if (geomH->ContainedYZ(sY, sZ) == false){
      std::stringstream ss;
      ss << "Fail @ algo " << this->name() << " due to start point reconstructed out of TPC bouns";
      throw ShowerRecoException(ss.str());
    }
    */

    resultShower.fXYZStart = {sX / 2, sY, sZ} ;

}
  DEFINE_ART_CLASS_TOOL(StartPoint3DModule)
} //showerreco

#endif
