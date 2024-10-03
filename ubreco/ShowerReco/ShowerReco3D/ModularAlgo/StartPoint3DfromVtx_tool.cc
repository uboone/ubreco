#ifndef STARTPOINT3DFROMVTX_CXX
#define STARTPOINT3DFROMVTX_CXX

#include <iostream>

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
    void do_reconstruction(const util::GeometryUtilities&,
                           const ::protoshower::ProtoShower &, Shower_t &);

  };
  
  StartPoint3DfromVtx::StartPoint3DfromVtx(const fhicl::ParameterSet& pset)
  {
    _name = "StartPoint3DfromVtx"; 
  }

  void StartPoint3DfromVtx::do_reconstruction(const util::GeometryUtilities&,
                                              const ::protoshower::ProtoShower & proto_shower,
					      Shower_t& resultShower)
{

    if (proto_shower.hasVertex() == false){
      std::cout << "Number of vertices is not one!" << std::endl;
      return;
    }
    
    // get the proto-shower 3D vertex
    auto const& vtx3D = proto_shower.vertex();

    auto start3D = vtx3D;

    resultShower.fXYZStart = start3D;

    //std::cout << "DONE " << std::endl << std::endl;
    
}

  DEFINE_ART_CLASS_TOOL(StartPoint3DfromVtx)
} //showerreco

#endif
