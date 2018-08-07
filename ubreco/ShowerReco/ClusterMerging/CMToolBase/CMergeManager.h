/**
 * \file CMergeManager.h
 *
 * \ingroup Clusmtool
 * 
 * \brief Class def header for a class CMergeManager
 *
 * @author kazuhiro
 */

/** \addtogroup Clusmtool

    @{*/
#ifndef RECOTOOL_CMERGEMANAGER_H
#define RECOTOOL_CMERGEMANAGER_H

#include <iostream>

#include <vector>

#include <typeinfo>

#include "CMManagerBase.h"
#include "CMergeBookKeeper.h"
#include "CBoolAlgoBase.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/make_tool.h"

namespace clusmtool {

  /**
     \class CMergeManager
     A class that instantiates merging algorithm(s) and run.
     The book-keeping of merged cluster sets are done by CMergeBookKeeper.
  */
  class CMergeManager : public CMManagerBase {

  public:
    
    /// Default constructor
    CMergeManager();
    
    /// copy constructor
    //CMergeManager(const CMergeManager&) = delete;
    // CMergeManager& operator=(const CMergeManager&) = delete;
    
    /// Method to reset itself
    virtual void Reset();

    /// A simple method to add an algorithm for merging
    void AddMergeAlgo(std::unique_ptr<clusmtool::CBoolAlgoBase> algo)
    //void AddMergeAlgo(const fhicl::ParameterSet& pset) 
    {
      std::cout << "DD \t\t\t currently there are "  << _merge_algo_v.size() << " algorithms" << std::endl;
      if (_merge_algo_v.size()) {
	std::cout << "DD \t\t\t report last algo's name " << std::endl;
	std::cout << "DD \t\t\t last algo added : " << _merge_algo_v.back()->Name() << std::endl;
      }
      _merge_algo_v.push_back(std::move(algo));
      //_merge_algo_v.push_back( art::make_tool<clusmtool::CBoolAlgoBase>(pset) ); 
      std::cout << "DD \t\t\t and now there are   "  << _merge_algo_v.size() << " algorithms" << std::endl;
    }

    /// A method to obtain output clusters
    const std::vector<::cluster::Cluster>& GetClusters() const { return _out_clusters; }

    /// A method to give up clusters
#ifndef __CINT__
    void PassOutputClusters(CMergeManager& rhs)
    {
      rhs.SetClusters(std::move(_out_clusters));
    }
#endif

    /// A method to obtain book keeper
    const CMergeBookKeeper& GetBookKeeper() const { return _book_keeper; }

    /**
       Report algorithm chain
     */
    void ReportAlgoChain();


  protected:
    
    //
    // FMWK functions override
    //

    /// FMWK function called @ beginning of Process()
    virtual void EventBegin();

    /// FMWK function called @ beginning of iterative loop inside Process()
    virtual void IterationBegin();

    /// FMWK function called @ iterative loop inside Process()
    virtual bool IterationProcess();

    /// FMWK function called @ end of iterative loop inside Process()
    virtual void IterationEnd();
    
    /// FMWK function called @ end of Process()
    virtual void EventEnd();

    void RunMerge(const int& algo_idx,
		  const std::vector<::cluster::Cluster > &in_clusters,
		  CMergeBookKeeper &book_keeper) const;

    void RunMerge(const int& algo_idx,
		  const std::vector<::cluster::Cluster > &in_clusters,
		  const std::vector<bool> &merge_flag,
		  CMergeBookKeeper &book_keeper) const;

    /// Output clusters
    std::vector<cluster::Cluster> _out_clusters;

    /// Book keeper instance
    CMergeBookKeeper _book_keeper;

    size_t _iter_ctr = 0;

    std::vector<CMergeBookKeeper> _book_keeper_v;

    std::vector<std::vector<unsigned short> > _tmp_merged_indexes;

    std::vector<::cluster::Cluster> _tmp_merged_clusters;

    /// Merging algorithm
    std::vector<std::unique_ptr<::clusmtool::CBoolAlgoBase> > _merge_algo_v;
    

  };
}

#endif
/** @} */ // end of doxygen group 

