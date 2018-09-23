/**
 * \file ProtoShowerAlgBase.h
 *
 * \ingroup ProtoShower
 *
 * \brief Class def header for a class ProtoShowerAlgBase
 *
 * @author david caratelli
 */

/** \addtogroup ProtoShower

    @{*/
#ifndef PROTOSHOWERALGBASE_H
#define PROTOSHOWERALGBASE_H

#include <iostream>

#include "art/Framework/Principal/Event.h"

#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/FindManyInChainP.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "ProtoShower.h"

// art TOOLS
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/ParameterSet.h"

/**
   \class ProtoShowerAlgBase
   User defined class ProtoShowerAlgBase ... these comments are used to generate
   doxygen documentation!
 */

namespace protoshower {

class ProtoShowerAlgBase {

public:

  /// Default constructor
  virtual ~ProtoShowerAlgBase() noexcept = default;

  void configure(const fhicl::ParameterSet& pset){};
    
  /// Default destructor
  //virtual ~ProtoShowerAlgBase() {}

  /**
     @brief Generate ProtoShower objects starting with from art event
     @input art:Event e -> event information
     @input fPFPproducer -> producer name for PFParticle
     @input fClusterproducer -> producer for Clusters associated to PFParticle
     @input fVtxproducer -> producer for vertex which notes shower origin
     @input proto_shower_v -> vector of protoshowers passed by reference. Filled by function
   */
  virtual void GenerateProtoShowers(::art::Event & e,
				    std::vector<protoshower::ProtoShower> & proto_shower_v) = 0;



  std::string name() { return _name; }
  
  /**
     Set calorimetry conversion for ADC to MeV
     ADC refers to Hit Integral()
   */
  void setCalorimetry(std::vector<double> c) { _ADC_to_MeV = c; }

protected:

  std::string _name;

  double _wire2cm, _time2cm;

  // conversion from ADC to MeV
  std::vector<double> _ADC_to_MeV;

};

}// namespace

#endif
/** @} */ // end of doxygen group

