#include <map>
#include <algorithm>

#include <iostream>
#include "ubreco/ShowerReco/ClusterMerging/CMToolBase/CBoolAlgoBase.h"

#include <Eigen/Dense>
#include <cmath>

namespace clusmtool {

  /**
     \class CBAlgoVtxAlign
     An abstract fake class for merging algorithm. Having this fake class helps
     to have a better overall design of various merging for iterative approach.
     The algorithms are run through CMergeManager.
  */

  class CBAlgoVtxAlign : public CBoolAlgoBase {
    
  public:
    
    /// Default constructor
    CBAlgoVtxAlign(const fhicl::ParameterSet& pset);
    
    /// Default destructor
    virtual ~CBAlgoVtxAlign(){};

    void configure(const fhicl::ParameterSet& pset);
 
    /**
       Core function: given the ClusterParamsAlg input, return whether a cluster should be
       merged or not.
    */
    std::vector<std::vector<size_t> > Merge(const std::vector<::cluster::Cluster>& clus_v);


    /// Function to reset the algorithm instance ... maybe implemented via child class
    void Reset(){}

  protected:

    float ClusterDistance(const ::cluster::Cluster& c1, const ::cluster::Cluster& c2);
    float ComputeClusterPCA(const ::cluster::Cluster& c);

    float _max_angle_diff_merge;
    float _min_gammagamma_oangle;
    float _max_merge_dist; // fraction of length of larger shower
    size_t _min_nhits;

    bool _run_shower_merge_algo;
  };

  //----------------------------------------
  CBAlgoVtxAlign::CBAlgoVtxAlign(const fhicl::ParameterSet& pset) 
  //----------------------------------------
  {
    _name = "CBAlgoVtxAlgin";
    configure(pset);

  }

  void CBAlgoVtxAlign::configure(const fhicl::ParameterSet& pset)
  {
    _max_angle_diff_merge  = pset.get<float>("max_angle_diff_merge");
    _min_gammagamma_oangle = pset.get<float>("min_gammagamma_oangle");
    _min_nhits             = pset.get<size_t>("min_nhits");
    _max_merge_dist        = pset.get<float>("max_merge_dist");
    _pair_wise             = pset.get<bool >("pair_wise");
    _merge_till_converge   = pset.get<bool >("mergetillconverge");
    _verbose               = pset.get<bool> ("verbose",false);
    _run_shower_merge_algo = pset.get<bool> ("run_shower_merge_algo",true);

    std::cout << "DD \t MIN NHITS "  << _min_nhits << std::endl;
    std::cout << "DD \t verbose is " << _verbose   << std::endl;

    return;
  }

  std::vector< std::vector<size_t> > CBAlgoVtxAlign::Merge(const std::vector<::cluster::Cluster>& clus_v) {

    std::vector< std::vector<size_t> > merge_result{};

    for (size_t pl=0; pl < 3; pl++) {

      if (_verbose)
	std::cout << std::endl << std::endl << std::endl << "****** PLANE " << pl << " *****" << std::endl << std::endl << std::endl << std::endl;

      // split clusters per plane:
      // and order by size
      std::vector< std::pair<size_t,size_t> > plane_clus_pairs_v; // pair [size, index]
      for (size_t i=0; i < clus_v.size(); i++) {
	if (clus_v[i]._plane == pl) {
	  plane_clus_pairs_v.push_back( std::make_pair( clus_v[i].size(), i ) );
	}
      }

      std::sort(plane_clus_pairs_v.rbegin(), plane_clus_pairs_v.rend());

      std::vector<size_t> plane_clus_idx_v;
      for (size_t i=0; i < plane_clus_pairs_v.size(); i++)
	plane_clus_idx_v.push_back( plane_clus_pairs_v[i].second );


      
      // pair links cluster index to pair of ( NHits, angle )
      // object will hold cluster angle for cluster with largest number of hits per plane
      std::pair< size_t, std::pair<size_t, double> > gamma_00(std::make_pair( 0, std::make_pair( 0, 0.) ) );
      std::pair< size_t, std::pair<size_t, double> > gamma_01(std::make_pair( 0, std::make_pair( 0, 0.) ) );


      // map linking cluster index and angle
      std::map<size_t, double> clus_angle_map;

      /*
      // sort clusters by size
      // map contains nhits -> index
      std::map<int,int> clus_size_map;

      for (auto const& idx : plane_clus_idx_v)
	clus_size_map[ clus_v[idx].size() ] = idx;
      */

      for (auto const& idx : plane_clus_idx_v) {
	
	//std::map<int,int>::reverse_iterator it;
	//for (it = clus_size_map.rbegin(); it != clus_size_map.rend(); it++) {

	//auto idx = it->second;
	auto const& cluster = clus_v.at(idx);

	auto nhits = cluster.size();
	auto angle = cluster._angle;

	clus_angle_map[ idx ] = angle;

	if (_verbose){
	  std::cout << "New cluster @ Plane " << pl
		    << "\t Nhits = " << nhits
		    << "\t angle = " << cluster._angle << std::endl
		    << "\t radius start / end = " << cluster._start_pt._r << " / " << cluster._end_pt._r << std::endl
		    << "Current gamma00 angle : " << gamma_00.second.second << "\t gamma01 angle : " << gamma_01.second.second << std::endl;
	}

	if ( nhits > gamma_00.second.first && nhits > _min_nhits ) {
	  if (_verbose) std::cout << "\t updating gamma00" << std::endl;
	  gamma_00.second.first  = nhits; // update # of hits for largest cluster in plane
	  gamma_00.second.second = angle; // update direction from PCA for largest cluster in plane
	  gamma_00.first         = idx;   // update index
	}
	
	if ( ( nhits > gamma_01.second.first ) &&
	     ( nhits > _min_nhits ) &&
	     ( fabs(angle-gamma_00.second.second ) > _min_gammagamma_oangle) &&
	     ( gamma_00.second.second != 0 ) ) {
	  if (_verbose) std::cout << "\t updating gamma01" << std::endl;
	  gamma_01.second.first  = nhits; // update # of hits for largest cluster in plane
	  gamma_01.second.second = angle; // update direction from PCA for largest cluster in plane
	  gamma_01.first         = idx;   // update index
	}
	
	if (_verbose) std::cout << std::endl << std::endl;

      }// for all clusters on this plane

      if (_verbose)
	std::cout << " gamma00 has " << gamma_00.second.first << " hits"
		  << " and " << gamma_00.second.second << " angle" << std::endl;
      if (_verbose)
	std::cout << " gamma01 has " << gamma_01.second.first << " hits"
		  << " and " << gamma_01.second.second << " angle" << std::endl;

      // "big" cluster 00
      auto gammacluster_00 = clus_v.at( gamma_00.first );
      auto gammaangle_00   = gamma_00.second.second;
      
      // "big" cluster 01
      auto gammacluster_01 = clus_v.at( gamma_01.first );
      auto gammaangle_01   = gamma_01.second.second;
      
      // loop through indices. which to merge?
      // keep track in a vector
      std::vector<size_t> indices_to_merge_00;
      std::vector<size_t> indices_to_merge_01;
      
      for (auto const& clus : clus_angle_map) {
	
	auto clus_idx   = clus.first;
	auto clus_angle = clus.second;
	auto cluster    = clus_v.at(clus_idx);

	// skip clusters with less then threshold of hits
	if (cluster.size() < _min_nhits) continue;
	
	// don't merge cluster with itself!
	if ( (clus.first == gamma_00.first) || ((clus.first == gamma_01.first) && (gammaangle_01 != 0)) ) continue;
	
	// find cluster with respect to which angle is smaller
	double angle_small = fabs(clus_angle - gammaangle_00);
	double angle_large = fabs(clus_angle - gammaangle_01);
	int merge_with = 0; // which gamma to merge with? [0 -> 00, 1 -> 01]
	// flip if we got it wrong
	if (angle_small > angle_large){
	  angle_small = fabs(clus_angle - gammaangle_01);
	  angle_large = fabs(clus_angle - gammaangle_00);
	  merge_with = 1;
	}
	
	if (_verbose)
	  std::cout << "new cluster with " << cluster.size()
		    << " hits and angle " << clus_angle << std::endl
		    << "\t angle small : " << angle_small << "\t angle large : " << angle_large
		    << "\t merge with=" << merge_with
		    << std::endl;
	
	// NOTE
	// below IF statements are not correct
	// when gammaangle_01 is undefined the angle_large requirement may cause
	// merging not to happen when it really should
	
	// Now the IF statements are fixed!
	if (merge_with == 0) {
	  if ( (angle_small < _max_angle_diff_merge) and
	       (gammaangle_00 != 0) and
	       ((angle_large > _min_gammagamma_oangle) || (gammaangle_01 == 0)) ){
	    if ( ClusterDistance(cluster, gammacluster_00) > _max_merge_dist * gammacluster_00.Length() ) continue;
	    indices_to_merge_00.push_back( clus_idx );
	    if (_verbose) std::cout << "MERGE WITH 00" << std::endl;
	  }
	}
	
	if (merge_with == 1) {
	  if ( (angle_small < _max_angle_diff_merge) and
	       (gammaangle_01 != 0) and
	       (angle_large > _min_gammagamma_oangle) ){
	    if ( ClusterDistance(cluster, gammacluster_01) > _max_merge_dist * gammacluster_01.Length() ) continue;
	    indices_to_merge_01.push_back( clus_idx );
	    if (_verbose) std::cout << "MERGE WITH 01" << std::endl;
	  }
	}
	
      }// for all clusters
      
      if (_verbose) {
	std::cout << "identified " << indices_to_merge_00.size() << " clusters to merge with gamma00" << std::endl;
	std::cout << "identified " << indices_to_merge_01.size() << " clusters to merge with gamma01" << std::endl;
      }

      if (indices_to_merge_00.size()) {
	std::vector<size_t> indices_to_merge_v {gamma_00.first};
	for (auto const& idx : indices_to_merge_00) indices_to_merge_v.push_back( idx );
	merge_result.push_back( indices_to_merge_v );
      }
      
      if (indices_to_merge_01.size()) {
	std::vector<size_t> indices_to_merge_v {gamma_01.first};
	for (auto const& idx : indices_to_merge_01) indices_to_merge_v.push_back( idx );
	merge_result.push_back( indices_to_merge_v );
      }

      // merging the start of the shower with its trunk
      if (_verbose) std::cout << "merge_result size: " << merge_result.size() << std::endl;
      if ( (gamma_00.second.second == 0) || (gamma_01.second.second == 0) ) continue;
      size_t gamma00_idx = gamma_00.first;
      size_t gamma01_idx = gamma_01.first;

      size_t gammastart_idx = gamma00_idx;
      size_t gammatrunk_idx = gamma01_idx;
      if (clus_v.at(gamma00_idx)._start_pt._r > clus_v.at(gamma01_idx)._start_pt._r){
        gammastart_idx = gamma01_idx;
        gammatrunk_idx = gamma00_idx;
      }
      bool gamma00isstart = (gammastart_idx == gamma00_idx);
      bool gamma01isstart = (gammastart_idx == gamma01_idx);
      bool gamma00istrunk = (gammatrunk_idx == gamma00_idx);
      bool gamma01istrunk = (gammatrunk_idx == gamma01_idx);
      float starttrunk_dist = std::sqrt(
        pow((clus_v.at(gammastart_idx)._end_pt._w - clus_v.at(gammatrunk_idx)._start_pt._w), 2) +
        pow((clus_v.at(gammastart_idx)._end_pt._t - clus_v.at(gammatrunk_idx)._start_pt._t), 2));

      // compute cluster PCA angle
      float anglePCA_gamma00 = ComputeClusterPCA(clus_v.at(gamma00_idx));
      float anglePCA_gamma01 = ComputeClusterPCA(clus_v.at(gamma01_idx));
      float diff_anglePCA = fabs(anglePCA_gamma00 - anglePCA_gamma01);

      if (_verbose) {
        std::cout<<std::endl;
        std::cout<<"****************************************"<<std::endl;
        std::cout<<"Plane:\t"<<pl<<std::endl;
        std::cout<<"\t\tgamma00\t\t"<<"gamma01"<<std::endl;
        std::cout<<"n hits:\t\t"<<clus_v.at(gamma00_idx).size()<<"\t\t"
                                <<clus_v.at(gamma01_idx).size()<<std::endl;
        std::cout<<"start w:\t" <<clus_v.at(gamma00_idx)._start_pt._w<<"\t\t"
                                <<clus_v.at(gamma01_idx)._start_pt._w<<std::endl;
        std::cout<<"start t:\t" <<clus_v.at(gamma00_idx)._start_pt._t<<"\t\t"
                                <<clus_v.at(gamma01_idx)._start_pt._t<<std::endl;
        std::cout<<"end w:\t\t" <<clus_v.at(gamma00_idx)._end_pt._w<<"\t\t"
                                <<clus_v.at(gamma01_idx)._end_pt._w<<std::endl;
        std::cout<<"end t:\t\t" <<clus_v.at(gamma00_idx)._end_pt._t<<"\t\t"
                                <<clus_v.at(gamma01_idx)._end_pt._t<<std::endl;
        std::cout<<"start r:\t" <<clus_v.at(gamma00_idx)._start_pt._r<<"\t\t"
                                <<clus_v.at(gamma01_idx)._start_pt._r<<std::endl;
        std::cout<<"end r:\t\t" <<clus_v.at(gamma00_idx)._end_pt._r<<"\t\t"
                                <<clus_v.at(gamma01_idx)._end_pt._r<<std::endl;
        std::cout<<"is start:\t"<<gamma00isstart<<"\t\t"<<gamma01isstart<<std::endl;
        std::cout<<"is trunk:\t"<<gamma00istrunk<<"\t\t"<<gamma01istrunk<<std::endl;
        std::cout<<"start to trunk dist:\t" << starttrunk_dist << std::endl;
        std::cout<<"angle:\t\t" <<clus_v.at(gamma00_idx)._angle<<"\t\t"
                                <<clus_v.at(gamma01_idx)._angle<<std::endl;
        std::cout<<"anglePCA:\t"<<anglePCA_gamma00<<"\t\t"
                                <<anglePCA_gamma01<<std::endl;
        std::cout<<"angle span min:\t"<<clus_v.at(gamma00_idx)._angle_span._amin<<"\t\t"
                                      <<clus_v.at(gamma01_idx)._angle_span._amin<<std::endl;
        std::cout<<"angle span max:\t"<<clus_v.at(gamma00_idx)._angle_span._amax<<"\t\t"
                                      <<clus_v.at(gamma01_idx)._angle_span._amax<<std::endl;
        std::cout<<"****************************************"<<std::endl;
        std::cout<<std::endl;
      }

      if (_run_shower_merge_algo) {
        //if ( (clus_v.at(gamma00_idx)._start_pt._r < 5) || (clus_v.at(gamma01_idx)._start_pt._r < 5) )
        if ((starttrunk_dist < 10 && diff_anglePCA < 15) or (starttrunk_dist < 20 && diff_anglePCA < 10)) {

          if(_verbose) std::cout << "need to merge gamma00 and gamma01: starttrunk_dist = " << starttrunk_dist << std::endl;
          if ( (indices_to_merge_00.size() == 0) && (indices_to_merge_01.size() == 0) ) {
            merge_result.push_back( std::vector<size_t> {gamma00_idx, gamma01_idx} );
          }
          else {
            std::vector<size_t> indices_to_merge {gamma_00.first, gamma_01.first};
            indices_to_merge.insert(indices_to_merge.end(), indices_to_merge_00.begin(), indices_to_merge_00.end());
            indices_to_merge.insert(indices_to_merge.end(), indices_to_merge_01.begin(), indices_to_merge_01.end());
            if ( (indices_to_merge_00.size() != 0) ^ (indices_to_merge_01.size() != 0) ) {
              merge_result.pop_back();
            }
            else if ( (indices_to_merge_00.size() != 0) && (indices_to_merge_01.size() != 0) ) {
              merge_result.pop_back();
              merge_result.pop_back();
            }
            merge_result.push_back(indices_to_merge);
          }
          if(_verbose){
            std::cout << "merge_result size: " << merge_result.size() << std::endl;
            std::cout << "identified " << merge_result.back().size() <<
            " clusters (including gamma00 and gamma01) to merge together " << std::endl;
          }

        }
      }


    }// for all 3 planes



    return merge_result;
  }


    // calculate maximum gap between clusters
  float CBAlgoVtxAlign::ClusterDistance(const ::cluster::Cluster& c1,
					const ::cluster::Cluster& c2) {


    if ( (c1._start_pt._r > c2._start_pt._r) && (c1._start_pt._r < c2._end_pt._r) )
      return 0.;

    if ( (c2._start_pt._r > c1._start_pt._r) && (c2._start_pt._r < c1._end_pt._r) )
      return 0.;

    if (c1._start_pt._r > c2._end_pt._r) return c1._start_pt._r - c2._end_pt._r;
    
    if (c2._start_pt._r > c1._end_pt._r) return c2._start_pt._r - c1._end_pt._r;

    return -1;
  }

  // compute cluster PCA angle
  float CBAlgoVtxAlign::ComputeClusterPCA(const ::cluster::Cluster& c) {

    auto hits = c.GetHits();

    float sumWeights = 0;
    float xx = 0;
    float yy = 0;
    float xy = 0;

    float sum_charge = 0;

    float ShowerCentre_c[2] = {0,0};
    for (auto const& ht : hits) {
      ShowerCentre_c[0] += ht._w*ht._q;
      ShowerCentre_c[1] += ht._t*ht._q;
      sum_charge += ht._q;
    }
    for (int i=0; i<2; i++) ShowerCentre_c[i] *= (1. / sum_charge);
    Eigen::Vector2f ShowerCenter(ShowerCentre_c[0], ShowerCentre_c[1]);

    for (auto const& ht : hits) {
      float wht = 1;
      wht *= std::sqrt(ht._q / sum_charge);

      // Normalise the hit point position.
      Eigen::Vector2f ht_pos(ht._w, ht._t);
      auto const ht_pos_norm = ht_pos - ShowerCenter;

      xx += ht_pos_norm.x() * ht_pos_norm.x() * wht;
      yy += ht_pos_norm.y() * ht_pos_norm.y() * wht;
      xy += ht_pos_norm.x() * ht_pos_norm.y() * wht;
      sumWeights += wht;
    }

    // Using Eigen package
    Eigen::Matrix2f matrix;

    // Construct covariance matrix
    matrix << xx, xy, xy, yy;

    // Normalise from the sum of weights
    matrix /= sumWeights;

    // Run the PCA
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2f> eigenMatrix(matrix);

    // Eigen::Vector2f eigenValuesVector = eigenMatrix.eigenvalues();
    Eigen::Matrix2f eigenVectorsMatrix = eigenMatrix.eigenvectors();
    Eigen::Vector2f primaryEigenVector = eigenVectorsMatrix.col(1);

    Eigen::Vector2f StartToEnd(c._end_pt._w-c._start_pt._w, c._end_pt._t-c._start_pt._t);
    float dotProd = primaryEigenVector.dot(StartToEnd);
    float mag_vec_A = primaryEigenVector.norm();
    float mag_vec_B = StartToEnd.norm();
    float cos_dir = dotProd / (mag_vec_A * mag_vec_B);

    float PCA = 0;
    //If the PCA axis is opposite to the shower cluster direction, flip the PCA axis.
    if (cos_dir >= 0) {PCA = std::atan2(primaryEigenVector.y(), primaryEigenVector.x());}
    else {PCA = std::atan2(-primaryEigenVector.y(), -primaryEigenVector.x());}
    float PCAdegrees = PCA * 180.0 / M_PI;

    float anglePCA = PCAdegrees + 180.0;

    return anglePCA;

  }

DEFINE_ART_CLASS_TOOL(CBAlgoVtxAlign)  
}
