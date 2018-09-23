#include "ProtoShowerAlgBase.h"

#include "larpandora/LArPandoraObjects/PFParticleMetadata.h"

namespace protoshower {

  class ProtoShowerPandora : public ProtoShowerAlgBase {

  public:
    
    /// default constructor
    explicit ProtoShowerPandora(const fhicl::ParameterSet& pset);
    
    /// default destructor
    ~ProtoShowerPandora(){};
    
    void configure(const fhicl::ParameterSet& pset);

  /**
     @brief function which takes recob::Cluster and vector of recob::Hits to create cluster2d::Cluster2D object
     @input art::Ptr to cluster
     @input vector of art::Ptr to hits associated to the cluster
  */
    cluster2d::Cluster2D MakeCluster2D( const art::Ptr<recob::Cluster>& clus, const std::vector< art::Ptr<recob::Hit> >& hit_v);

   
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

    bool fNeutrino;
    double fNeutrinoScoreMin;
    double fTrackScoreMax;

  };

  ProtoShowerPandora::ProtoShowerPandora(const fhicl::ParameterSet& pset) {
    _name = "ProtoShowerPandora";
    configure(pset);
  }

  void ProtoShowerPandora::configure(const fhicl::ParameterSet& pset) {
    fPFPproducer     = pset.get<std::string>("PFPproducer"    );
    fClusterproducer = pset.get<std::string>("ClusterProducer");
    fVtxproducer     = pset.get<std::string>("Vtxproducer"    );
    fNeutrino        = pset.get<bool       >("Neutrino"       );
    fNeutrinoScoreMin  = pset.get<double     >("NeutrinoScoreMin" );
    fTrackScoreMax     = pset.get<double     >("TrackScoreMax"    );

    // get detector specific properties
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    _wire2cm = geom->WirePitch(0,0,0);
    _time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );

  }

  void ProtoShowerPandora::GenerateProtoShowers(::art::Event & e,
						std::vector<protoshower::ProtoShower> & proto_shower_v) {

    // grab PFParticles in event
    auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle> >(fPFPproducer);

    // grab clusters associated with PFParticles
    art::FindManyP<recob::Cluster> pfp_clus_assn_v(pfp_h, e, fPFPproducer);

    // grab associated metadata
    art::FindManyP< larpandoraobj::PFParticleMetadata > pfPartToMetadataAssoc(pfp_h, e, fPFPproducer);
    //art::FindManyP< larpandoraobj::PFParticleMetadata > pfp_meta_assn_v(pfp_h, e, fPFPproducer);

    // grab vertices associated
    art::FindManyP< recob::Vertex > pfp_vtx_assn_v(pfp_h, e, fVtxproducer);

    // ADDITION FROM PETRILLO
    e.getValidHandle<std::vector<recob::Cluster>>(fClusterproducer);

    // grab the hits associated to the PFParticles
    auto pfp_hit_assn_v = lar::FindManyInChainP<recob::Hit, recob::Cluster>::find(pfp_h, e, fPFPproducer);

    std::cout << "Grabbed necessary products " << std::endl;

    /*
    // Get mapping from ID to PFParticle
    std::map<size_t, art::Ptr<recob::PFParticle> > pfp_map;
    for (unsigned int i = 0; i < pfp_h->size(); ++i) {
      const art::Ptr<recob::PFParticle> pfp(pfp_h, i);
      pfp_map.emplace(pfp->Self(), pfp);
      //throw cet::exception("WorkshopAnalyzer") << "Repeated PFParticles!" << std::endl;
    } 
    */
    
    std::cout << "Start event loop " << std::endl;

    // loop through PFParticles
    for (size_t p=0; p < pfp_h->size(); p++) {


      const recob::PFParticle pfp = pfp_h->at(p);

      std::cout << "Now I have a PFP! " << std::endl;

      // get metadata for this PFP
      const std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> > &pfParticleMetadataList(pfPartToMetadataAssoc.at(p));

      // we need to decide if to skip this particle based on user settings
      bool skip = false;
      
      if (!pfParticleMetadataList.empty()) {
	std::cout << "metadata!" << std::endl;
	
	for (unsigned int j=0; j<pfParticleMetadataList.size(); ++j)
	  {
	    const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
	    const pandora::PropertiesMap &pfParticlePropertiesMap(pfParticleMetadata->GetPropertiesMap());
	    if (!pfParticlePropertiesMap.empty())
	      std::cout << " Found PFParticle " << pfp.Self() << " with: " << std::endl;
	    for (pandora::PropertiesMap::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it) {
	      std::cout << "  - " << it->first << " = " << it->second << std::endl;
	      if ( (it->first == "IsClearCosmic") && (it->second == 1) && (fNeutrino == true) ) {
		std::cout << "\t SKIPPING because ClearCosmic " << std::endl;
		skip = true;
	      }// if this is  not a neutrino
	      if ( (it->first == "TrackScore") && (it->second > fTrackScoreMax) ){
		std::cout << "\t SKIPPING because TrackScore is  " << it->second << std::endl;
		skip = true;
	      }// if this is not a shower
	    }
	  }
      }// if PFP metadata exists!
      
      if (skip == true)
	continue;

      /*
      const std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> > &pfp_meta_v(pfp_meta_assn_v.at(p));
      if (pfp_meta_v.empty() == true) continue;
      for (unsigned int j=0; j < pfp_meta_v.size(); j++) {
	const art::Ptr<larpandoraobj::PFParticleMetadata> &pfp_meta(pfp_meta_v.at(j));
	const pandora::PropertiesMap &pfp_map(pfp_meta->GetPropertiesMap());
	if (pfp_map.empty() == true) continue;
	std::cout << " Found PFParticle " << pfp.Self() << " with: " << std::endl;
	for (pandora::PropertiesMap::const_iterator it = pfp_map.begin(); it != pfp_map.end(); ++it)
	  std::cout << "  - " << it->first << " = " << it->second << std::endl;
      }
      */

      // Find the parent particle
      auto parentIdx = pfp.Parent();
      if (parentIdx > pfp_h->size() ) continue;

      auto parent = pfp_h->at( parentIdx );

      // get metadata for parent
      const std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> > &parentMetadataList(pfPartToMetadataAssoc.at(parent.Self()));

      if (!parentMetadataList.empty()) {
	std::cout << "parent metadata!" << std::endl;
	
	for (unsigned int j=0; j<parentMetadataList.size(); ++j)
	  {
	    const art::Ptr<larpandoraobj::PFParticleMetadata> &parentMetadata(parentMetadataList.at(j));
	    const pandora::PropertiesMap &parentPropertiesMap(parentMetadata->GetPropertiesMap());
	    if (!parentPropertiesMap.empty())
	      std::cout << " Found PFParticle " << parent.Self() << " with: " << std::endl;
	    for (pandora::PropertiesMap::const_iterator it = parentPropertiesMap.begin(); it != parentPropertiesMap.end(); ++it) {
	      std::cout << "  - " << it->first << " = " << it->second << std::endl;
	      if ( (it->first == "NuScore") && (it->second <= fNeutrinoScoreMin) && (fNeutrino == true) ) {
		std::cout << "\t SKIPPING because IsNeutrino Score is " << it->second << std::endl;
		skip = true;
	      }// if this is  not a neutrino
	    }
	  }
      }// if PFP metadata exists!
      
      // find vertex defined as the vertex of the primary

      const std::vector< art::Ptr<recob::Vertex> >& vtx_v = pfp_vtx_assn_v.at( parentIdx );//pfp_parent_iter->first );

      std::cout << "grabbed vertex" << std::endl;

      ::protoshower::ProtoShower proto_shower;
      proto_shower.Reset();
      
      
      // require a single vertex!
      if (vtx_v.size() == 1) {
	
	auto const vtx = vtx_v.at(0);
	Double_t xyz[3] = {};
	vtx->XYZ(xyz);
	proto_shower._vertex = TVector3(xyz[0],xyz[1],xyz[2]);
	proto_shower.hasVertex(true);
	
	std::cout << "\t\t DD VTX @ [ " << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << " ]" << std::endl;
	
      }// if there is only one vertex
      else
	continue;
      
      // associated clusters
      std::cout << "grabbed cluster" << std::endl;
      const std::vector< art::Ptr<recob::Cluster> >& clus_v = pfp_clus_assn_v.at(p);
      
      // associated hits
      std::cout << "grabbed hits" << std::endl;
      const std::vector< art::Ptr<recob::Hit> >& hit_v = pfp_hit_assn_v.at(p);
      
      // set number of clusters for protoshower
      std::cout << "PFP associated with " << clus_v.size() << " clusters..." << std::endl;
      proto_shower._clusters.resize(clus_v.size());

      // loop through clusters
      for (size_t c=0; c < clus_v.size(); c++) {

	auto const clus = clus_v.at(c);
	std::cout << " new cluster" << std::endl;

	// find the hits associated to this cluster
	std::vector< art::Ptr<recob::Hit> > clusterhits;
	for (auto const& hit : hit_v) {
	  // if hit in same plane as cluster -> add to vector
	  if (hit->WireID().Plane == clus->Plane().Plane )
	    clusterhits.push_back( hit );
	}//for all hits associated to PFParticle

	std::cout << " make cluster2d for this cluster with " << clusterhits.size() << " hits" << std::endl;
	proto_shower._clusters.at(c) = MakeCluster2D( clus, clusterhits );
	
	proto_shower.hasCluster2D(true);

      }// for all clusters
      
      std::cout << "done with hit loop" << std::endl;
      
      proto_shower_v.push_back( proto_shower );
      
    }// for all PFParticles
    
  }// GenerateProtoShower end


  ::cluster2d::Cluster2D ProtoShowerPandora::MakeCluster2D( const art::Ptr<recob::Cluster>& clus, 
							    const std::vector< art::Ptr<recob::Hit> >& hit_v) 
  {
    
    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    ::cluster2d::Cluster2D clus2d;
    clus2d.Reset();

    for (auto const hit : hit_v) {

      // create PxHit
      ::util::PxHit hit2d( clus->Plane().Plane,
			   hit->WireID().Wire * _wire2cm,
			   (hit->PeakTime() - detp->TriggerOffset()) * _time2cm,
			   hit->Integral() * _ADC_to_MeV[clus->Plane().Plane],
			   hit->Integral() * _ADC_to_MeV[clus->Plane().Plane],
			   hit->PeakAmplitude() );

      clus2d._hits.push_back( hit2d );
    }// for all hits

    clus2d._plane = clus->Plane().Plane;

    // Missing : plane offset w.r.t. origin coordinates

    auto const& sw = clus->StartWire() * _wire2cm;
    auto const& ew = clus->EndWire()   * _wire2cm;
    auto const& st = (clus->StartTick() - detp->TriggerOffset()) * _time2cm;
    auto const& et = (clus->EndTick()   - detp->TriggerOffset()) * _time2cm;
    
    clus2d._start = ::util::PxHit(clus->Plane().Plane, sw, st, 0., 0., 0.);
    clus2d._end   = ::util::PxHit(clus->Plane().Plane, ew, et, 0., 0., 0.);
    
    clus2d._angle_2d = clus->StartAngle();

    return clus2d;
    
  }



DEFINE_ART_CLASS_TOOL(ProtoShowerPandora)
}// namespace

