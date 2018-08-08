// base class
#include <iostream>
#include "ubreco/ShowerReco/ClusterMerging/CMToolBase/CFloatAlgoBase.h"

#include "TTree.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"

namespace clusmtool {


  class CFAlgoTimeOverlap : public CFloatAlgoBase {
    
  public:
    
    /// Default constructor
    explicit CFAlgoTimeOverlap(const fhicl::ParameterSet& pset);
    
    /// Default destructor
    ~CFAlgoTimeOverlap(){};

    /**This algorithm calculates the difference between start and end times for merged clusters,
		and compares across planes to form matches. 
    */
    float Float(const std::vector<const cluster::Cluster*> &clusters);

    void Report();
    
    void Reset();

  protected:

    void getMinMaxTime(const cluster::Cluster* cluster, double& min, double& max);

    void configure(const fhicl::ParameterSet& pset);

    /**
       Is a hit matched on the other plane?
     */
    int Matched(const float& time, const std::vector<cluster::pt>& otherhits);
    
    float _timetolerance;

    /*
    TTree* _match_tree;
    int _induction; // which induction plane are we matching to?
    double _compat; // compatibility score
    */

  };
  
  
  //-------------------------------------------------------
  CFAlgoTimeOverlap::CFAlgoTimeOverlap(const fhicl::ParameterSet& pset) 
  //-------------------------------------------------------
  {
    _name = "CFAlgoTimeOverlap";
    configure(pset);
  }

  //--------------------------------------------------------
  void CFAlgoTimeOverlap::configure(const fhicl::ParameterSet& pset)
  //--------------------------------------------------------
  {
    _timetolerance = pset.get<float>("timetolerance");
    _verbose       = pset.get<bool> ("verbose",false);

    std::cout << "TTTTTTT configuring" << std::endl;

    /*
    art::ServiceHandle<art::TFileService> tfs;
    _match_tree = tfs->make<TTree>("match_tree","Match TTree");
    _match_tree->Branch("_induction",&_induction,"induction/I");
    _match_tree->Branch("_compat"   ,&_compat   ,"compat/D"   );
    */

    return;
  }

  //-----------------------------
  void CFAlgoTimeOverlap::Reset()
  //-----------------------------
  {

  }

  //----------------------------------------------------------------------------------------------
  float CFAlgoTimeOverlap::Float(const std::vector<const cluster::Cluster*> &clusters)
  //----------------------------------------------------------------------------------------------
  {


    std::cout << " BLABLABLA" << std::endl;

    // if 3 clusters -> skip
    if (clusters.size() != 2) return -1;

    // identify collection plane
    // if no collection-plane cluster -> ignore match
    size_t collectionPlane = 0;
    size_t inductionPlane  = 0;
    bool   collectionPresent = false;
    if (clusters[0]->_plane == 2) {
      collectionPresent = true;
      collectionPlane   = 0;
      inductionPlane    = 1;
    }
    if (clusters[1]->_plane == 2) {
      collectionPresent = true;
      collectionPlane   = 1;
      inductionPlane    = 0;
    }
    
    if (collectionPresent == false)
      return -1;

    double t_min_abs = 9600; // smallest start point of the 3
    double t_max_abs = 0;    // largest start point of the three

    
    for(auto const& c : clusters){

      double min,max;

      getMinMaxTime(c,min,max);

      if (_verbose) {
	std::cout << "cluster bounds :  [ " << min << ", " << max << "]" 
		  << " cluster size   : " << c->size() 
		  << std::endl;
      }

      if ( min < t_min_abs )
	t_min_abs = min;
      if ( max   > t_max_abs )
	t_max_abs = max;
      
    }
    
    if (clusters.size() < 2) return -1;

    // do the clusters overlap at all?
    bool overlap = true;

   double t_min_common, t_max_common;

    getMinMaxTime(clusters[0],t_min_common,t_max_common);

    for (size_t i=1; i < clusters.size(); i++){

      if (_verbose)
	std::cout << "PLANES " << clusters[0]->_plane
		  << " AND " << clusters[i]->_plane
		  << std::endl;

      double t_min, t_max; // min and max for current cluster
      getMinMaxTime(clusters[i],t_min,t_max);

      if (_verbose) {
	std::cout << "Current t_min,t_max    = [" <<  t_min_common << ", " << t_max_common << " ]"<< std::endl
		  << "Current cluster bounds = [" <<  t_min        << ", " << t_max        << " ]"<< std::endl;
	  }
      
      // now find overlap
      if (t_max < t_min_common){
	overlap = false;
	break;
      }

      if (t_min > t_max_common){
	overlap = false;
	break;
      }

      if (t_min > t_min_common) t_min_common = t_min;
      if (t_max < t_max_common) t_max_common = t_max;
      
      if (_verbose) {
	std::cout << "updated t_min,t_max    = [" <<  t_min_common << ", " << t_max_common << " ]"<< std::endl;
      }
      
    }// for all clusters
    
    if (overlap == false) return -1;

    // if the clusters overlap to some degree, check compatibility on hit-by-hit level
    double compatibility = 0.;
    auto const& hitsCollection = clusters[collectionPlane]->GetHits();
    auto const& hitsInduction  = clusters[inductionPlane ]->GetHits();
    
    for (auto const& hitC : hitsCollection) 
      compatibility += Matched(hitC._t,hitsInduction);

    compatibility /= (double)(hitsCollection.size());

    //_induction = clusters[inductionPlane]->_plane;
    //_compat    = compatibility;
    //std::cout << "Induction @ plane " << _induction << " with compatibility = " << _compat << std::endl;

    //_match_tree->Fill();
    
    return compatibility;
  }
  
  
  //------------------------------
  void CFAlgoTimeOverlap::Report()
  //------------------------------
  {
  }
  

  int CFAlgoTimeOverlap::Matched(const float& time, const std::vector<cluster::pt>& otherhits) {

    for (auto const& hit : otherhits) {

      if ( fabs(hit._t - time) < _timetolerance) 
	return 1;
      
    }// for all hits in other cluster

    return 0;
  }
  
  void CFAlgoTimeOverlap::getMinMaxTime(const cluster::Cluster* cluster, double& min, double& max)
  {
    
    auto const& hits = cluster->GetHits();

    min = 9600;
    max = 0;

    for (auto const& hit : hits){
      
      if (hit._t > max) max = hit._t;
      if (hit._t < min) min = hit._t;

    }// fo rall hits
    
    return;
  }
  
  DEFINE_ART_CLASS_TOOL(CFAlgoTimeOverlap)    
}

