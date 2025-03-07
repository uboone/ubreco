#include <iostream>
#include "ubreco/ShowerReco/ClusterMerging/CMToolBase/CBoolAlgoBase.h"

namespace clusmtool {

  class CBAlgoPolar : public CBoolAlgoBase {
    
  public:
    
    /// Default constructor
    explicit CBAlgoPolar(const fhicl::ParameterSet& pset);
    
    /// Default destructor
    ~CBAlgoPolar(){};

    void configure(const fhicl::ParameterSet& pset);
 
    /**
       Core function: given the ClusterParamsAlg input, return whether a cluster should be
       merged or not.
    */
    bool Bool(const ::cluster::Cluster &cluster1,
	      const ::cluster::Cluster &cluster2);


    /// Function to reset the algorithm instance ... maybe implemented via child class
    void Reset(){}

  protected:

    float _buffer;
    
    float ClusterDistance(const ::cluster::Cluster& c1, const ::cluster::Cluster& c2);
    
  };

  //----------------------------------------
  CBAlgoPolar::CBAlgoPolar(const fhicl::ParameterSet& pset) 
  //----------------------------------------
  {
    _name = "CBAlgoPolar";
    configure(pset);
  }

  void CBAlgoPolar::configure(const fhicl::ParameterSet& pset)
  {
    _buffer              = pset.get<float>("buffer");
    _pair_wise           = pset.get<bool >("pair_wise");
    _merge_till_converge = pset.get<bool >("mergetillconverge");
    _verbose             = pset.get<bool> ("verbose",false);

    std::cout << "DD \t BUFFER IS " << _buffer << std::endl;
    std::cout << "DD \t verbose is " << _verbose << std::endl;

  }

  //--------------------------------------------------------
  bool CBAlgoPolar::Bool(const ::cluster::Cluster &cluster1,
			 const ::cluster::Cluster &cluster2)
  //--------------------------------------------------------
  {

    if (cluster1._plane != cluster2._plane)
      std::cout << "Different planes!" << std::endl;

    auto clusDistance = ClusterDistance(cluster1,cluster2);

    if (_verbose) {
      std::cout << "cluster1 start R = " << cluster1._start_pt._r
		<< "\t Qtot = " << cluster1._sum_charge
		<< "\t Nhits = " << cluster1.size() << std::endl
		<< "cluster2 start R = " << cluster2._start_pt._r
		<< "\t Qtot = " << cluster2._sum_charge
		<< "\t Nhits = " << cluster2.size() << std::endl
		<< "Angle span 1 : [" << cluster1._angle_span._amin << ", " << cluster1._angle_span._amax << " ]" << std::endl
		<< "Angle span 2 : [" << cluster2._angle_span._amin << ", " << cluster2._angle_span._amax << " ]" << std::endl
		<< "Angle 1 : " << cluster1._angle << std::endl
		<< "Angle 2 : " << cluster2._angle << std::endl
		<< "Length 1 : " << cluster1.Length() << std::endl
		<< "Length 2 : " << cluster2.Length() << std::endl
		<< "cluster distance = " << clusDistance << std::endl
		<< "plane = " << cluster1._plane << std::endl
		<< std::endl;
    }

    if (clusDistance < 0) {
      std::cout << "Error. Distance is  " << clusDistance << std::endl;
      return false;
    }

    // if cluster 1 upstream of cluster 2 and larger
    if ( (cluster1._start_pt._r < cluster2._start_pt._r) &&
	 (cluster1._sum_charge > cluster2._sum_charge) ) {

      // if the smaller cluster is within the span of the larger -> merge
      if ( cluster1._angle_span.inRange(cluster2._angle, _buffer) == true ) {

	if (clusDistance < cluster1.Length() ) {
	  
	  if (_verbose) { std::cout << "merge 2 with 1 " << std::endl << std::endl; }
	  
	  return true;

	}// if cluster distance not larger than largest cluster's length
	
      }
      
    }// END : if cluster 1 upstream of cluster 2 and larger
    else if ( (cluster2._start_pt._r < cluster1._start_pt._r) &&
	      (cluster2._sum_charge > cluster1._sum_charge) ) {
      
      // if the smaller cluster is within the span of the larger -> merge
      if ( cluster2._angle_span.inRange(cluster1._angle, _buffer) == true ) {

	if (clusDistance < cluster2.Length() ) {
	  
	  if (_verbose) { std::cout << "merge 1 with 2 " << std::endl << std::endl; }
	  
	  return true;

	}// if cluster distance not larger than largest cluster's length
	
      }
      
    }// if cluster 2 upstream of cluster 1 and larger

    return false;
  }

  // calculate maximum gap between clusters
  float CBAlgoPolar::ClusterDistance(const ::cluster::Cluster& c1,
				     const ::cluster::Cluster& c2) {

    if ( (c1._start_pt._r >= c2._start_pt._r) && (c1._end_pt._r <= c2._end_pt._r) )
      return 0.;

    if ( (c2._start_pt._r >= c1._start_pt._r) && (c2._end_pt._r <= c1._end_pt._r) )
      return 0.;

    if ( (c1._start_pt._r >= c2._start_pt._r) && (c1._start_pt._r < c2._end_pt._r) )
      return 0.;

    if ( (c2._start_pt._r >= c1._start_pt._r) && (c2._start_pt._r < c1._end_pt._r) )
      return 0.;

    if (c1._start_pt._r > c2._end_pt._r) return c1._start_pt._r - c2._end_pt._r;
    
    if (c2._start_pt._r > c1._end_pt._r) return c2._start_pt._r - c1._end_pt._r;

    std::cout << "C1 [ " << c1._start_pt._r << ", " << c1._end_pt._r << " ]" << std::endl;
    std::cout << "C2 [ " << c2._start_pt._r << ", " << c2._end_pt._r << " ]" << std::endl;
    
    return -1;
  }

DEFINE_ART_CLASS_TOOL(CBAlgoPolar)
}

