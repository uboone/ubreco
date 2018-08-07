#include <iostream>
#include "ubreco/ShowerReco/ClusterMerging/CMToolBase/CPriorityAlgoBase.h"

namespace clusmtool {

  /**
     \class CPAlgoQSum
     Simple algorithm to determine priority based on charge sum
     If charge sum < set cut value by a user, returns -1.     
  */
  class CPAlgoQSum : public CPriorityAlgoBase {
    
  public:
    
    /// Default constructor
    explicit CPAlgoQSum(const fhicl::ParameterSet& pset);
    
    /// Default destructor
    ~CPAlgoQSum(){};

    void configure(const fhicl::ParameterSet& pset);

    /**
       Core function: given the CPAN input, return a float which indicates 
       the user-defined priority for analysis.
    */
    //float Priority(const ::cluster::Cluster &cluster);

  protected:

    double _qsum_cut;

  };

  //----------------------------------------------------------
  CPAlgoQSum::CPAlgoQSum(const fhicl::ParameterSet& pset)
  //----------------------------------------------------------
  {
    configure(pset);
  }

  void CPAlgoQSum::configure(const fhicl::ParameterSet& pset)
  {
    _qsum_cut = pset.get<double>("qsum_cut");
  }

  /*
  //------------------------------------------------------------------------------
  float CPAlgoQSum::Priority(const ::cluster::Cluster &cluster)
  //------------------------------------------------------------------------------
  {
    //if(cluster._sum_charge < _qsum_cut) return -1;

    return cluster._sum_charge;
  }
  */
  
  DEFINE_ART_CLASS_TOOL(CPAlgoQSum)        
}
