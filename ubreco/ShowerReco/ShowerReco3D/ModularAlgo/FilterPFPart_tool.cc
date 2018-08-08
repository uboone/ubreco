#ifndef FILTERPFPART_CXX
#define FILTERPFPART_CXX

#include <iostream>
#include "ubreco/ShowerReco/ShowerReco3D/Base/ShowerRecoModuleBase.h"

#include <sstream>

/**
   \class ShowerRecoModuleBase
   User defined class ShowerRecoModuleBase ... these comments are used to generate
   doxygen documentation!
 */

namespace showerreco {

  class FilterPFPart : public ShowerRecoModuleBase {
    
  public:
    
    /// Default constructor
    FilterPFPart(const fhicl::ParameterSet& pset);
    
    /// Default destructor
    ~FilterPFPart() {}

    void configure(const fhicl::ParameterSet& pset);
    
    void do_reconstruction(const ::protoshower::ProtoShower &, Shower_t &);
    
    void setMinNHitsAbsolute(int n) { _min_nhits_absolute = n; }
    void setMinNHitsLargest (int n) { _min_nhits_largest  = n; }
    
  private:
    
    size_t _min_nhits_absolute, _min_nhits_largest;
    
  };

  FilterPFPart::FilterPFPart(const fhicl::ParameterSet& pset)
  {
    configure(pset);
    _name = "FilterPFPart";
  }
  
  void FilterPFPart::configure(const fhicl::ParameterSet& pset)
  {
    _verbose   = pset.get<bool>("verbose",false);
    _min_nhits_absolute = pset.get<size_t>("min_nhits_absolute");
    _min_nhits_largest  = pset.get<size_t>("min_nhits_largest" );
  }

  void FilterPFPart::do_reconstruction(const ::protoshower::ProtoShower & proto_shower,
				       Shower_t& resultShower) {
    
    //if the module does not have 2D cluster info -> fail the reconstruction
    if (!proto_shower.hasCluster2D()) {
      std::stringstream ss;
      ss << "Fail @ algo " << this->name() << " due to missing 2D cluster";
      throw ShowerRecoException(ss.str());
    }
    
    auto & clusters = proto_shower.clusters();
    
    // keep track of largest cluster
    size_t largest_cluster_nhits = 0;
    
    // loop through clusters, find number of hits in each
    for (size_t n = 0; n < clusters.size(); n++) {
      
      size_t nhits = clusters.at(n)._hits.size();
      
      if (nhits < _min_nhits_absolute) {
	std::stringstream ss;
	ss << "Fail @ algo " << this->name() << " : input cluster has too few hits";
	throw ShowerRecoException(ss.str());
      }
      
      if (nhits > largest_cluster_nhits)
	largest_cluster_nhits = nhits;
      
    }// for all clusters
    
    if (largest_cluster_nhits < _min_nhits_largest) {
      std::stringstream ss;
      ss << "Fail @ algo " << this->name() << " : largest cluster too small";
      throw ShowerRecoException(ss.str());
    }
    
    // This function takes the shower cluster set and computes the best fit 3D axis
    // and then assigns it to the shower
    
  }

  DEFINE_ART_CLASS_TOOL(FilterPFPart)
} //showerreco

#endif
