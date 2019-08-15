#ifndef STARTPOINT3DFROMVTX_CXX
#define STARTPOINT3DFROMVTX_CXX

#include <iostream>

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "ubreco/ShowerReco/ShowerReco3D/Base/ShowerRecoModuleBase.h"
/**
   \class ShowerRecoModuleBase
   User defined class ShowerRecoModuleBase ... these comments are used to generate
   doxygen documentation!
 */

namespace showerreco {

  class StartPoint3DfromVtx : public ShowerRecoModuleBase {
    
  public:
    
    /// Default constructor
    StartPoint3DfromVtx(const fhicl::ParameterSet& pset);

    /// Default destructor
    ~StartPoint3DfromVtx() {}
    
    /// Inherited/overloaded function from ShowerRecoModuleBase
    void do_reconstruction(const ::protoshower::ProtoShower &, Shower_t &);

  private:

    double _wire2cm, _time2cm;
    
  };
  
  StartPoint3DfromVtx::StartPoint3DfromVtx(const fhicl::ParameterSet& pset)
  {
    _name = "StartPoint3DfromVtx"; 
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    _wire2cm = geom->WirePitch(0,0,0);
    _time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );

  }

  void StartPoint3DfromVtx::do_reconstruction( const ::protoshower::ProtoShower & proto_shower,
					      Shower_t& resultShower)
{

    if (proto_shower.hasVertex() == false){
      std::cout << "Number of vertices is not one!" << std::endl;
      return;
    }
    
    // get the proto-shower 3D vertex
    auto const& vtx3D = proto_shower.vertex();

    auto start3D = vtx3D;

    std::cout << " Strt 3D : [" << start3D[0] << ", " << start3D[1] << ", " << start3D[2] << " ]" << std::endl;

    resultShower.fXYZStart = start3D;

    //std::cout << "DONE " << std::endl << std::endl;
    
}

  DEFINE_ART_CLASS_TOOL(StartPoint3DfromVtx)
} //showerreco

#endif
