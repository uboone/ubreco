#ifndef EMPTYMODULE_CXX
#define EMPTYMODULE_CXX

#include <iostream>
#include "ubreco/ShowerReco/ShowerReco3D/Base/ShowerRecoModuleBase.h"
/**
   \class ShowerRecoModuleBase
   User defined class ShowerRecoModuleBase ... these comments are used to generate
   doxygen documentation!
 */

namespace showerreco {

  class EmptyModule : public ShowerRecoModuleBase {
    
  public:
    
    /// Default constructor
    EmptyModule(const fhicl::ParameterSet& pset);
    
    /// Default destructor
    ~EmptyModule(){}
    
    
    void do_reconstruction(const ::protoshower::ProtoShower &, Shower_t &);
    
  private:
    
  };

  EmptyModule::EmptyModule(const fhicl::ParameterSet& pset)
  {
    _name = "EmptyModule";
  }
  
  void EmptyModule::do_reconstruction(const ::protoshower::ProtoShower & proto_shower,
				      Shower_t& resultShower) {
    
    // This function takes the shower cluster set and computes the best fit 3D axis
    // and then assigns it to the shower
    
  }
  
  DEFINE_ART_CLASS_TOOL(EmptyModule)
} //showerreco

#endif
