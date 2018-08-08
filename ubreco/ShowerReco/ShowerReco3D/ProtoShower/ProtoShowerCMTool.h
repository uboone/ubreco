/**
 * \file ProtoShowerCMTool.h
 *
 * \ingroup ProtoShower
 *
 * \brief Class def header for a class ProtoShowerCMTool
 *
 * @author david caratelli
 */

/** \addtogroup ProtoShower

    @{*/
#ifndef PROTOSHOWER_PROTOSHOWERALGBASE_H
#define PROTOSHOWER_PROTOSHOWERALGBASE_H

#include <iostream>

#include "ProtoShowerAlgBase.h"

/**
   \class ProtoShowerCMTool
   User defined class ProtoShowerCMTool ... these comments are used to generate
   doxygen documentation!
 */

namespace protoshower {

  class ProtoShowerCMTool : public ProtoShowerAlgBase {

public:

  /// Default constructor
  ProtoShowerCMTool() { _name = "ProtoShowerCMTool"; }

  /// Default destructor
  virtual ~ProtoShowerCMTool() {}

  /**
     @brief Generate ProtoShower objects starting with from art event
     @input art:Event e -> event information
     @input fPFPproducer -> producer name for PFParticle
     @input fClusterproducer -> producer for Clusters associated to PFParticle
     @input fVtxproducer -> producer for vertex which notes shower origin
     @input proto_shower_v -> vector of protoshowers passed by reference. Filled by function
   */
  void GenerateProtoShowers(::art::Event & e,
			    const std::string& fPFPproducer,
			    const std::string& fClusterproducer,
			    const std::string& fVtxproducer,
			    std::vector<protoshower::ProtoShower> & proto_shower_v);
  

  std::string name() { return _name; }

  //void setVertexProducer(const std::string& producer) { fVertexProducer = producer; }

protected:

  std::string _name;

};

}// namespace

#endif
/** @} */ // end of doxygen group

