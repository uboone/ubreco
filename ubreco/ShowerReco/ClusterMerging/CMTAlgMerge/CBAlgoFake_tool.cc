#include <iostream>
#include "ubreco/ShowerReco/ClusterMerging/CMToolBase/CBoolAlgoBase.h"

namespace clusmtool {

  /**
     \class CBAlgoFake
     An abstract fake class for merging algorithm. Having this fake class helps
     to have a better overall design of various merging for iterative approach.
     The algorithms are run through CMergeManager.
  */
  class CBAlgoFake : public CBoolAlgoBase {
    
  public:
    
    /// Default constructor
    explicit CBAlgoFake(const fhicl::ParameterSet& pset);
    
    /// Default destructor
    ~CBAlgoFake(){};
 
    /**
       Core function: given the ClusterParamsAlg input, return whether a cluster should be
       merged or not.
    */
    bool Bool(const ::cluster::Cluster &cluster1,
	      const ::cluster::Cluster &cluster2);

    std::vector<std::vector<size_t> > Merge(const std::vector<::cluster::Cluster>& clus_v);


    /// Function to reset the algorithm instance ... maybe implemented via child class
    void Reset(){}

  protected:

    bool _flip;
    int _ctr;
  };

  //----------------------------------------
  CBAlgoFake::CBAlgoFake(const fhicl::ParameterSet& pset)
  //----------------------------------------
  {
    _flip = false;
    _ctr  = 0;
    // Nothing to be done in the base class
  }

  //--------------------------------------------------------
  bool CBAlgoFake::Bool(const ::cluster::Cluster &cluster1,
                        const ::cluster::Cluster &cluster2)
  //--------------------------------------------------------
  {

    std::cout << "Pair-wise called" << std::endl;
    
    if(cluster1.size() && cluster2.size()) {
      _ctr++;
      if( (_ctr%64) == 0)
        _flip = (!_flip);
      return _flip;
    }
    else return false;
  }

  std::vector< std::vector<size_t> > CBAlgoFake::Merge(const std::vector<::cluster::Cluster>& clus_v) {

    std::cout << "Full view merging called" << std::endl;

    std::vector< std::vector<size_t> > merge_result{ {}, {}, {} };

    for (size_t clus_idx = 0; clus_idx < clus_v.size(); clus_idx++) {

      auto const& clus = clus_v.at(clus_idx);

      merge_result[clus._plane].push_back(clus_idx);
      
    }

    return merge_result;
  }

DEFINE_ART_CLASS_TOOL(CBAlgoFake)  
}

