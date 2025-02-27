// base class
#include "ubreco/ShowerReco/ClusterMerging/CMToolBase/CFloatAlgoBase.h"



namespace clusmtool {


  class CFAlgoIoU : public CFloatAlgoBase {
    
  public:
    
    /// Default constructor
    explicit CFAlgoIoU(const fhicl::ParameterSet& pset);
    
    /// Default destructor
    ~CFAlgoIoU(){};

    /**This algorithm calculates the difference between start and end times for merged clusters,
		and compares across planes to form matches. 
    */
    float Float(const std::vector<const cluster::Cluster*> &clusters);

    void Report();
    
    void Reset();

  protected:

    void getMinMaxTime(const cluster::Cluster* cluster, double& min, double& max);

    void configure(const fhicl::ParameterSet& pset);
    
    float _iou_min;
  };
  
  
  //-------------------------------------------------------
  CFAlgoIoU::CFAlgoIoU(const fhicl::ParameterSet& pset) 
  //-------------------------------------------------------
  {
    _name = "CFAlgoIoU";
    configure(pset);
  }

  //--------------------------------------------------------
  void CFAlgoIoU::configure(const fhicl::ParameterSet& pset)
  //--------------------------------------------------------
  {
    _iou_min = pset.get<float>("iou_min");
    _verbose = pset.get<bool> ("verbose",false);
    return;
  }

  //-----------------------------
  void CFAlgoIoU::Reset()
  //-----------------------------
  {

  }

  //----------------------------------------------------------------------------------------------
  float CFAlgoIoU::Float(const std::vector<const cluster::Cluster*> &clusters)
  //----------------------------------------------------------------------------------------------
  {



    // // if 3 clusters -> skip
    // if (clusters.size() != 2) return -1;

    if (clusters.size()==3) {
      std::vector<const cluster::Cluster*> vec01 = {clusters[0],clusters[1]};
      float score01 = this->Float(vec01);
      if (score01<0.) return -1.;
      std::vector<const cluster::Cluster*> vec02 = {clusters[0],clusters[2]};
      float score02 = this->Float(vec02);
      if (score02<0.) return -1.;
      std::vector<const cluster::Cluster*> vec12 = {clusters[1],clusters[2]};
      float score12 = this->Float(vec12);
      if (score12<0.) return -1.;

      return score01+score02+score12;
    }

    // // require collection plane
    // if ( (clusters[0]->_plane != 2) && (clusters[1]->_plane != 2) ) return -1;

    if ( (clusters[0]->size() < 10) || (clusters[1]->size() < 10) ) return -1;

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
    
    if (_verbose) {
      std::cout << "T common interval : [" <<  t_min_common << ", " << t_max_common << " ]"<< std::endl
		<< "T total  interval : [" <<  t_min_abs    << ", " << t_max_abs    << " ]"<< std::endl;
    }	
    
    // calculate overlap
    double iou = (t_max_common - t_min_common) / (t_max_abs - t_min_abs);

    if (_verbose)
      std::cout << "Cluster IoU : " << iou
		<< std::endl << std::endl;

    if (iou < _iou_min) return -1;

    return iou;
  }
  
  
  //------------------------------
  void CFAlgoIoU::Report()
  //------------------------------
  {
  }
  
  
  void CFAlgoIoU::getMinMaxTime(const cluster::Cluster* cluster, double& min, double& max)
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
  
  DEFINE_ART_CLASS_TOOL(CFAlgoIoU)    
}
