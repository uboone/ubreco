#include "ProtoShowerAlgBase.h"

namespace protoshower {
  
  class ProtoShowerCMTool : public ProtoShowerAlgBase {

  public:
    
    /// default constructor
    explicit ProtoShowerCMTool(const fhicl::ParameterSet& pset);
    
    /// default destructor
    ~ProtoShowerCMTool(){};
    
    void configure(const fhicl::ParameterSet& pset);

  /**
     @brief function which takes recob::Cluster and vector of recob::Hits to create cluster2d::Cluster2D object
     @input art::Ptr to cluster
     @input vector of art::Ptr to hits associated to the cluster
  */
    cluster2d::Cluster2D MakeCluster2D(const detinfo::DetectorClocksData& clockData,
                                       const art::Ptr<recob::Cluster>& clus, const std::vector< art::Ptr<recob::Hit> >& hit_v);
   
  /**
     @brief Generate ProtoShower objects starting with from art event
     @input art:Event e -> event information
     @input fPFPproducer -> producer name for PFParticle
     @input fClusterproducer -> producer for Clusters associated to PFParticle
     @input fVtxproducer -> producer for vertex which notes shower origin
     @input proto_shower_v -> vector of protoshowers passed by reference. Filled by function
   */ 
    void GenerateProtoShowers(::art::Event & e,
			      std::vector<protoshower::ProtoShower> & proto_shower_v);

  private:

    std::string fPFPproducer, fClusterproducer, fVtxproducer;

  };

  ProtoShowerCMTool::ProtoShowerCMTool(const fhicl::ParameterSet& pset) {
    _name = "ProtoShowerCMTool";
    configure(pset);
  }

  void ProtoShowerCMTool::configure(const fhicl::ParameterSet& pset) {
    fPFPproducer     = pset.get<std::string>("PFPproducer"    );
    fClusterproducer = pset.get<std::string>("ClusterProducer");
    fVtxproducer     = pset.get<std::string>("Vtxproducer"    );

    // get detector specific properties
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();
    auto const detp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);
    _wire2cm = geom->WirePitch(0,0,0);
    _time2cm = sampling_rate(clockData) / 1000.0 * detp.DriftVelocity( detp.Efield(), detp.Temperature() );
  }

  void ProtoShowerCMTool::GenerateProtoShowers(::art::Event & e,
					       std::vector<protoshower::ProtoShower> & proto_shower_v) {

    // grab PFParticles in event
    auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle> >(fPFPproducer);

    // grab clusters associated with PFParticles
    art::FindManyP<recob::Cluster> pfp_clus_assn_v(pfp_h, e, fPFPproducer);

    // ADDITION FROM PETRILLO
    e.getValidHandle<std::vector<recob::Cluster>>(fClusterproducer);

    // grab the hits associated to the PFParticles
    auto pfp_hit_assn_v = lar::FindManyInChainP<recob::Hit, recob::Cluster>::find(pfp_h, e, fPFPproducer);

    // load event vertex associated to tagged neutrino interaction
    auto const& vertex_h = e.getValidHandle<std::vector<recob::Vertex> >(fVtxproducer);

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
    // loop through PFParticles
    for (size_t p=0; p < pfp_h->size(); p++) {

      ::protoshower::ProtoShower proto_shower;
      proto_shower.Reset();

      const recob::PFParticle pfp = pfp_h->at(p);

      // associated clusters
      const std::vector< art::Ptr<recob::Cluster> >& clus_v = pfp_clus_assn_v.at(p);

      // associated hits
      const std::vector< art::Ptr<recob::Hit> >& hit_v = pfp_hit_assn_v.at(p);

      // set number of clusters for protoshower
      proto_shower._clusters.resize(clus_v.size());

      // loop through clusters
      for (size_t c=0; c < clus_v.size(); c++) {

	auto const clus = clus_v.at(c);

	// find the hits associated to this cluster
	std::vector< art::Ptr<recob::Hit> > clusterhits;
	for (auto const& hit : hit_v) {
	  // if hit in same plane as cluster -> add to vector
	  if (hit->WireID().Plane == clus->Plane().Plane )
	    clusterhits.push_back( hit );
	}//for all hits associated to PFParticle

        proto_shower._clusters.at(c) = MakeCluster2D( clockData, clus, clusterhits );
	
	proto_shower.hasCluster2D(true);

      }// for all clusters
      
      // require a single vertex!
      if (vertex_h->size() == 1) {
	
	auto const vtx = vertex_h->at(0);
	Double_t xyz[3] = {};
	vtx.XYZ(xyz);
	proto_shower._vertex = TVector3(xyz[0],xyz[1],xyz[2]);
	proto_shower.hasVertex(true);
	
	std::cout << "\t\t DD VTX @ [ " << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << " ]" << std::endl;
	
      }// if there is only one vertex
      
      proto_shower_v.push_back( proto_shower );
      
    }// for all PFParticles
    
  }// GenerateProtoShower end

  ::cluster2d::Cluster2D ProtoShowerCMTool::MakeCluster2D(const detinfo::DetectorClocksData& clockData,
                                                          const art::Ptr<recob::Cluster>& clus,
							    const std::vector< art::Ptr<recob::Hit> >& hit_v) 
  {
    ::cluster2d::Cluster2D clus2d;
    clus2d.Reset();
    
    for (auto const hit : hit_v) {

      // create PxHit
      ::util::PxHit hit2d( clus->Plane().Plane,
			   hit->WireID().Wire * _wire2cm,
                           (hit->PeakTime() - trigger_offset(clockData)) * _time2cm,
			   hit->Integral() * _ADC_to_MeV[clus->Plane().Plane],
			   hit->Integral() * _ADC_to_MeV[clus->Plane().Plane],
			   hit->PeakAmplitude() );
      
      clus2d._hits.push_back( hit2d );
    }// for all hits

    std::cout << "\t\t done looping through hits" << std::endl;
    
    clus2d._plane = clus->Plane().Plane;

    // Missing : plane offset w.r.t. origin coordinates

    auto const& sw = clus->StartWire() * _wire2cm;
    auto const& ew = clus->EndWire()   * _wire2cm;
    auto const& st = (clus->StartTick() - trigger_offset(clockData)) * _time2cm;
    auto const& et = (clus->EndTick()   - trigger_offset(clockData)) * _time2cm;
    
    clus2d._start = ::util::PxHit(clus->Plane().Plane, sw, st, 0., 0., 0.);
    clus2d._end   = ::util::PxHit(clus->Plane().Plane, ew, et, 0., 0., 0.);
    
    clus2d._angle_2d = clus->StartAngle();

    std::cout << "\t\t returning cluster2D" << std::endl;
    
    return clus2d;
    
  }


DEFINE_ART_CLASS_TOOL(ProtoShowerCMTool)  
}// namespace
