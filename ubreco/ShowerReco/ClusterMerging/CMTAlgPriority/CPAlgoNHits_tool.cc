#include <map>
#include <algorithm>

#include <iostream>
#include "ubreco/ShowerReco/ClusterMerging/CMToolBase/CPriorityAlgoBase.h"

namespace clusmtool {

  /**
     \class CPAlgoNHits
     Simple algorithm to determine priority based on # of hits.
     If # hits < set cut value by a user, returns -1.
  */
  class CPAlgoNHits : public CPriorityAlgoBase {
    
  public:
    
    /// Default constructor
    explicit CPAlgoNHits(const fhicl::ParameterSet& pset);
    
    /// Default destructor
    ~CPAlgoNHits(){};

    void configure(const fhicl::ParameterSet& pset);

    /**
       Core function: given the CPAN input, return a float which indicates 
       the user-defined priority for analysis.
    */
    /*
    float Priority(const ::cluster::Cluster &cluster) 
    {
      std::cout << "OVERRIDE!" << std::endl;
      return 100;
    }
    */

  protected:

    size_t _min_hits;

  };

  //----------------------------------------------
  CPAlgoNHits::CPAlgoNHits(const fhicl::ParameterSet& pset)
  //----------------------------------------------
  {
    _name = "CPAlgoNHits";
    configure(pset);
  }

  void CPAlgoNHits::configure(const fhicl::ParameterSet& pset)
  {
    std::cout << "MIN HITS :  " << pset.get<size_t>("min_hits") << std::endl;
    _min_hits = pset.get<size_t>("min_hits");
  }

  /*
  //------------------------------------------------------------------------
  float CPAlgoNHits::Priority(const ::cluster::Cluster &cluster)
  //------------------------------------------------------------------------
  {
    auto nhit = cluster.size();

    std::cout << "\t\t size is " << cluster.size() << std::endl;

    return ( nhit < _min_hits ? -1 : (float)nhit );
  }
  */

  DEFINE_ART_CLASS_TOOL(CPAlgoNHits)      
}
