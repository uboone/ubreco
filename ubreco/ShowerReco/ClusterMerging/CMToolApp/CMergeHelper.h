/**
 * \file CMergeHelper.h
 *
 * \ingroup ClusmtoolApp
 * 
 * \brief Class def header for a class CMergeHelper
 *
 * @author kazuhiro
 */

/** \addtogroup ClusmtoolApp

    @{*/
#ifndef CMERGEHELPER_H
#define CMERGEHELPER_H

#include <iostream>
#include "ubreco/ShowerReco/ClusterMerging/CMToolBase/CMergeManager.h"
#include "ubreco/ShowerReco/ClusterMerging/CMToolBase/ClusterMaker.h"

namespace clusmtool {
  /**
     \class CMergeHelper
     User defined class CMergeHelper ... these comments are used to generate
     doxygen documentation!
  */
  class CMergeHelper{
    
  public:
    
    /// Default constructor
    CMergeHelper(){}
    
    /// Default destructor
    //virtual ~CMergeHelper(){}

    //CMergeManager& GetManager(size_t mgr_id=0);
    CMergeManager& GetManager();

    void SetAnaFile(TFile* fout);

    void Process(const std::vector< ::cluster::Cluster >& clusters);

    //size_t size() const { return _mgr_v.size(); }

    const CMergeBookKeeper& GetResult() const { return _bk; }

    const std::vector< ::cluster::Cluster>& GetClusters() const;

  protected:

    ::clusmtool::CMergeManager _mgr;
    //std::vector< ::clusmtool::CMergeManager> _mgr_v;

    CMergeBookKeeper _bk;
    
  };
}

#endif
/** @} */ // end of doxygen group 

