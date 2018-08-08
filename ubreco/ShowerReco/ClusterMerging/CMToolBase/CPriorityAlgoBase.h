/**
 * \file CPriorityAlgoBase.h *
 * \ingroup Clusmtool
 * 
 * \brief Class def header for a class CPriorityAlgoBase
 *
 * @author kazuhiro
 */

/** \addtogroup Clusmtool

    @{*/
#ifndef RECOTOOL_CPRIORITYALGOBASE_H
#define RECOTOOL_CPRIORITYALGOBASE_H

#include "CMAlgoBase.h"

namespace clusmtool {

  /**
     \class CPriorityAlgoBase
     An abstract base class for CMatchManager and CMergeManager to determine 
     cluster "priority" for matching and merging action respectively.
  */
  class CPriorityAlgoBase : public CMAlgoBase {
    
  public:

    CPriorityAlgoBase(){}
    
    /// Default destructor
    virtual ~CPriorityAlgoBase() = default;

    /**
       Core function: given the CPAN input, return whether a cluster should be
       merged or not.
    */
    float Priority(const ::cluster::Cluster &cluster) 
    {
      return ( cluster.size() < 10 ? -1 : (float)cluster.size() );
      //std::cout << "PRIORITY BEING DEFINED "  << std::endl;
      //return cluster.size() * 1000;
      //
      //return -1;
    };

  };

}

#endif
/** @} */ // end of doxygen group 

