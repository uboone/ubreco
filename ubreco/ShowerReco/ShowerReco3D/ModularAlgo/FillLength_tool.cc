#ifndef FILLLENGTH_CXX
#define FILLLENGTH_CXX

#include <iostream>
#include "ubreco/ShowerReco/ShowerReco3D/Base/ShowerRecoModuleBase.h"
/**
   \class ShowerRecoModuleBase
   User defined class ShowerRecoModuleBase ... these comments are used to generate
   doxygen documentation!
 */

namespace showerreco {

  class FillLength : public ShowerRecoModuleBase {
    
  public:
    
    /// Default constructor
    FillLength(const fhicl::ParameterSet& pset);
    
    /// Default destructor
    ~FillLength() {}
    
    void do_reconstruction(const ::protoshower::ProtoShower &, Shower_t &);
    
  private:
    
  };
  
  FillLength::FillLength(const fhicl::ParameterSet& pset)
  {
    _name        = "FillLength";
  }
  
  void FillLength::do_reconstruction(const ::protoshower::ProtoShower & proto_shower,
				     Shower_t& resultShower) {
    
    //if the module does not have 2D cluster info -> fail the reconstruction
    if (!proto_shower.hasCluster2D()){
      std::stringstream ss;
      ss << "Fail @ algo " << this->name() << " due to missing 2D cluster";
      throw ShowerRecoException(ss.str());
    }
    
    auto & clusters = proto_shower.clusters();
    
    // use collection-plane cluster to compute length
    
    bool filled = false;
    
    for (auto const& clus : clusters) {
      
      if (clus._plane != 2) continue;
      
      auto sW  = clus._start.w;
      auto eW  = clus._end.w;
      
      // length on Y-plane projection
      auto dLY = fabs(eW - sW);
      
      resultShower.fLength = dLY / fabs(resultShower.fDCosStart[2]) ;
      
      // require a minimum length of 20 cm
      if (resultShower.fLength < 20) resultShower.fLength = 20.;
      
      resultShower.fOpeningAngle = clus._opening_angle;
      
      filled = true;
      
    }// for all clusters    
    
    // if collection-plane info not found -> skip
    if (filled == false) {
      std::stringstream ss;
      ss << "Fail @ algo " << this->name() << " due to missing collection-plane cluster";
      throw ShowerRecoException(ss.str());
    }
    
    return;
  }
  DEFINE_ART_CLASS_TOOL(FillLength)
}// showerreco

#endif

