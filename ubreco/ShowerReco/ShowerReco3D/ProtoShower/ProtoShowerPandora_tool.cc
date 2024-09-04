#include "ProtoShowerAlgBase.h"

#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/MVAOutput.h"

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

    /*
    // grab event hits and associations to NuGraph semantic score
    auto const& hitsWithScores = proxy::getCollection<std::vector<recob::Hit>>(
									       e,
									       art::InputTag("gaushit"), //tag of the hit collection we ran the GNN on
									       //proxy::withParallelData<anab::FeatureVector<1>>(art::InputTag("NuGraph", "filter")),
									       proxy::withParallelData<anab::FeatureVector<5>>(art::InputTag("gaushit", "semantic")));
    */


    // loop through PFParticles
    for (size_t p=0; p < pfp_h->size(); p++) {

      //std::cout << std::endl;
      //std::cout << "NEW PFP" << std::endl;

      const recob::PFParticle pfp = pfp_h->at(p);

      if ( (fabs(pfp.PdgCode()) ==12) || (fabs(pfp.PdgCode()) ==14) ) {
	continue;
      }
	  
      // get metadata for this PFP
      const std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> > &pfParticleMetadataList(pfPartToMetadataAssoc.at(p));
      
      // we need to decide if to skip this particle based on user settings
      bool skip = false;
      
      if (!pfParticleMetadataList.empty()) {
	
	for (unsigned int j=0; j<pfParticleMetadataList.size(); ++j)
	  {
	    const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
	    auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
	    //const larpandoraobj::PropertiesMap &pfParticlePropertiesMap(pfParticleMetadata->GetPropertiesMap());
	    if (!pfParticlePropertiesMap.empty())
	      //std::cout << " Found PFParticle " << pfp.Self() << " with PDG code " << pfp.PdgCode() << std::endl;
	    for (std::map<std::string, float>::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it) {
	      //std::cout << "  - " << it->first << " = " << it->second << std::endl;
	      if ( (it->first == "IsClearCosmic") && (it->second == 1) && (fNeutrino == true) ) {
		//std::cout << "\t SKIPPING because ClearCosmic " << std::endl;
		skip = true;
	      }// if this is  not a neutrino
	      if ( (it->first == "TrackScore") && (it->second > fTrackScoreMax) ){
		//std::cout << "\t SKIPPING because TrackScore is  " << it->second << std::endl;
		skip = true;
	      }// if this is not a shower
	      if ( (it->first == "NuScore") && (it->second <= fNeutrinoScoreMin) && (fNeutrino == true) ) {
		//std::cout << "\t SKIPPING because IsNeutrino Score is " << it->second << std::endl;
		skip = true;
	      }// if neutrino score too low
	    }
	  }
      }// if PFP metadata exists!
      
      if (skip == true){
	//std::cout << "\t SKIPPING after PFP metadata query " << std::endl;
	continue;
      }
      
      /*
      // Find the parent particle
      if (pfp.IsPrimary() == true) {
	//std::cout << "\t SKIPPING because primary PFP" << std::endl;
	continue;
      }
      */

      /*
      size_t parentIdx = pfp.Parent();
      if (parentIdx > pfp_h->size() ) {
	//std::cout << "\t SKIPPING because parent PFP idx is " << parentIdx << std::endl;
	continue;
      }
      */

      // save vertex for PFP
      std::vector< art::Ptr<recob::Vertex> > vtx_v = pfp_vtx_assn_v.at( p );;

      bool foundParent = false;
      // find parent based on pfp Self
      art::Ptr<recob::PFParticle> parent;
      for (size_t pj=0; pj < pfp_h->size(); pj++) {
	auto jpfp = art::Ptr<recob::PFParticle>(pfp_h,pj);
	if (jpfp->Self()==pfp.Parent()) {
	  foundParent = true;
	  parent = jpfp;
	  break;
	}
      }

      // if there is a parent, check if neutrino-like
      if (foundParent == true) {
	
	// get metadata for parent
	const std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> > &parentMetadataList(pfPartToMetadataAssoc.at(parent.key()));
	
	if (!parentMetadataList.empty()) {

	  //vtx_v = pfp_vtx_assn_v.at( parent.key() );
	  
	  for (unsigned int j=0; j<parentMetadataList.size(); ++j)
	    {
	      const art::Ptr<larpandoraobj::PFParticleMetadata> &parentMetadata(parentMetadataList.at(j));
	      auto parentPropertiesMap = parentMetadata->GetPropertiesMap();
	      //const larpandoraobj::PropertiesMap &parentPropertiesMap(parentMetadata->GetPropertiesMap());
	      if (!parentPropertiesMap.empty())
		//std::cout << " Found PFParticle " << parent->Self() << " with: " << std::endl;
	      for (std::map<std::string, float>::const_iterator it = parentPropertiesMap.begin(); it != parentPropertiesMap.end(); ++it) {
		//std::cout << "  - " << it->first << " = " << it->second << std::endl;
		if ( (it->first == "NuScore") && (it->second <= fNeutrinoScoreMin) && (fNeutrino == true) ) {
		  //std::cout << "\t SKIPPING because IsNeutrino Score is " << it->second << std::endl;
		  skip = true;
		}// if this is  not a neutrino
	      }
	    }
	}// if PFP metadata exists!

      }// if there is a parent
      /*
      else {
      }// if there is no parent
      */
      if (vtx_v.size() == 0) {
	std::cout << "\t SKIPPING because no vertex was found" << std::endl;
	continue;
      }

      //std::cout << "\t\t RECO THIS PFP" << std::endl;
      
      ::protoshower::ProtoShower proto_shower;
      proto_shower.Reset();
      proto_shower.SetIndex(p);
      
      // require a single vertex!
      if (vtx_v.size() == 1) {
	
	auto const vtx = vtx_v.at(0);
	Double_t xyz[3] = {};
	vtx->XYZ(xyz);
	proto_shower._vertex = TVector3(xyz[0],xyz[1],xyz[2]);
	proto_shower.hasVertex(true);
	
      }// if there is only one vertex
      else
	continue;
      
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

	proto_shower._clusters.at(c) = MakeCluster2D( clus, clusterhits );
	
	proto_shower.hasCluster2D(true);

      }// for all clusters

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

