#ifndef CMERGEHELPER_CXX
#define CMERGEHELPER_CXX

#include "CMergeHelper.h"

namespace clusmtool {

  ::clusmtool::CMergeManager& CMergeHelper::GetManager()//size_t mgr_id)
  {
    return _mgr;
    /*
      if(_mgr_v.size() <= mgr_id){
      while (_mgr_v.size() <= mgr_id)
	_mgr_v.emplace_back();
    }
    //if(_mgr_v.size() <= mgr_id) _mgr_v.resize(mgr_id+1);
    return _mgr_v[mgr_id];
    */
  }

  /*
  ::clusmtool::CMergeManager& CMergeHelper::GetManager(size_t mgr_id)
  {
    if(_mgr_v.size() <= mgr_id) _mgr_v.resize(mgr_id+1);
    return _mgr_v[mgr_id];
  }
  */

  void CMergeHelper::SetAnaFile(TFile* fout)
  {
    _mgr.SetAnaFile(fout);
    //for(auto& mgr : _mgr_v) mgr.SetAnaFile(fout);
  }

  void CMergeHelper::Process(const std::vector< ::cluster::Cluster >& clusters) 
  {
    _bk = ::clusmtool::CMergeBookKeeper(clusters.size());

    //for(size_t i=0; i<_mgr_v.size(); ++i) {
      
    //auto& mgr = _mgr_v[i];

    _mgr.Reset();
    _mgr.SetClusters(clusters);
    
    //if(!i) _mgr.SetClusters(clusters);
    //else _mgr_v[i-1].PassOutputClusters(mgr);
    
    _mgr.Process();
    
    auto const& new_bk = _mgr.GetBookKeeper();
    
    _bk = new_bk;

    //if(!i) _bk = new_bk;
    //else if(new_bk.GetResult().size() < new_bk.size())
    //  _bk.Combine(new_bk);
    
  }
  
  const std::vector< ::cluster::Cluster>& CMergeHelper::GetClusters() const
  {
    //if(!(_mgr_v.size())) throw CMTException("No manager = no output clusters...");
    return _mgr.GetClusters();
  }
}

#endif
