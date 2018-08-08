/**
 * \file CFloatAlgoBase.h
 *
 * \ingroup Clusmtool
 * 
 * \brief Class def header for a class CFloatAlgoBase
 *
 * @author kazuhiro
 */

/** \addtogroup Clusmtool

    @{*/
#ifndef RECOTOOL_CFLOATALGOBASE_H
#define RECOTOOL_CFLOATALGOBASE_H

#include "CMAlgoBase.h"

namespace clusmtool {

  /**
     \class CFloatAlgoBase
     An abstract base class for mathcing algorithm. Having this base class helps
     to have a better overall design of various merging for iterative approach.
     The algorithms are run through CMergeManager.
  */
  class CFloatAlgoBase : public CMAlgoBase {
    
  public:
    
    /// Default destructor
    virtual ~CFloatAlgoBase(){};

    /**
       Core function: given the CPAN input, return whether a cluster should be
       merged or not.
    */
    virtual float Float(const std::vector<const ::cluster::Cluster*> &clusters)
    { 
      if(clusters.size()) return 1;
      else return -1;
    }

  };

}

#endif
/** @} */ // end of doxygen group 

