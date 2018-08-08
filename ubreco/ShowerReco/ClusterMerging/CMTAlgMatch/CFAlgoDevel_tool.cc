// base class
#include "ubreco/ShowerReco/ClusterMerging/CMToolBase/CFloatAlgoBase.h"

#include "TTree.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"

namespace clusmtool {


  class CFAlgoDevel : public CFloatAlgoBase {
    
  public:
    
    /// Default constructor
    explicit CFAlgoDevel(const fhicl::ParameterSet& pset);
    
    /// Default destructor
    ~CFAlgoDevel(){};

    /**This algorithm calculates the difference between start and end times for merged clusters,
		and compares across planes to form matches. 
    */
    float Float(const std::vector<const cluster::Cluster*> &clusters);

    void Report();
    
    void Reset();

  protected:

    void configure(const fhicl::ParameterSet& pset);

    /**
       Is a hit matched on the other plane?
     */
    int Isochronous(const float& time, const std::vector<cluster::pt>& otherhits);

    float _timetolerance; // in cm, max disagreement in hit peak-time
    float _minoverlap; 

    TTree* _match_tree;
    double _compat;
    int    _induction;
    int    _nhit_induction;
    int    _nhit_collection;

  };
  
  
  //-------------------------------------------------------
  CFAlgoDevel::CFAlgoDevel(const fhicl::ParameterSet& pset) 
  //-------------------------------------------------------
  {
    _name = "CFAlgoDevel";
    configure(pset);

    art::ServiceHandle<art::TFileService> tfs;
    _match_tree = tfs->make<TTree>("match_tree","Match TTree");
    _match_tree->Branch("_induction",&_induction,"induction/I");
    _match_tree->Branch("_compat"   ,&_compat   ,"compat/D"   );
    _match_tree->Branch("_nhit_induction" ,&_nhit_induction  ,"nhit_induction/I" );
    _match_tree->Branch("_nhit_collection",&_nhit_collection,"nhit_collection/I");

  }

  //--------------------------------------------------------
  void CFAlgoDevel::configure(const fhicl::ParameterSet& pset)
  //--------------------------------------------------------
  {

    _timetolerance = pset.get<float>("timetolerance");
    _minoverlap    = pset.get<float>("minoverlap"   );
    _verbose       = pset.get<bool> ("verbose",false);
    return;
  }

  //-----------------------------
  void CFAlgoDevel::Reset()
  //-----------------------------
  {

  }

  //----------------------------------------------------------------------------------------------
  float CFAlgoDevel::Float(const std::vector<const cluster::Cluster*> &clusters)
  //----------------------------------------------------------------------------------------------
  {

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

    // if the clusters overlap to some degree, check compatibility on hit-by-hit level
    double compatibility = 0.;
    auto const& hitsCollection = clusters[collectionPlane]->GetHits();
    auto const& hitsInduction  = clusters[inductionPlane ]->GetHits();

    _nhit_induction  = hitsInduction.size();
    _nhit_collection = hitsCollection.size();
    
    for (auto const& hitC : hitsCollection) 
      compatibility += Isochronous(hitC._t,hitsInduction);
    
    compatibility /= (double)(hitsCollection.size());
    
    _induction = clusters[inductionPlane]->_plane;
    _compat    = compatibility;
    std::cout << "Induction @ plane " << _induction << " with compatibility = " << _compat << std::endl;
    
    _match_tree->Fill();

    if (_compat < _minoverlap) return -1;
    
    return compatibility;

  }

  //------------------------------
  void CFAlgoDevel::Report()
  //------------------------------
  {
  }

  int CFAlgoDevel::Isochronous(const float& time, const std::vector<cluster::pt>& otherhits) {

    for (auto const& hit : otherhits) {

      if ( fabs(hit._t - time) < _timetolerance) 
	return 1;
      
    }// for all hits in other cluster

    return 0;
  }

  
  
  DEFINE_ART_CLASS_TOOL(CFAlgoDevel)    
}
